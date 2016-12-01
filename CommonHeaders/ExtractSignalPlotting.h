// provided by Gamma Conversion Group, $ALICE_ROOT/PWGGA/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion

/*
//************************************************************************************
//************ Return correct trigger name based on trigger cutnumber ****************
//************************************************************************************
TString ReturnTriggerName(Int_t trigger){
    cout << trigger << endl;
    if (trigger == 1 || trigger == 3 || trigger == 0){  // INT1
        return "INT1";
    } else if (trigger == 10 ){                         // INT7
        return "INT7";
    } else if (trigger == 11){                          // INT7
        return "INT8";
    } else if (trigger == 51){                          // EMC1
        return "EMC1";
    } else if (trigger == 52){                          // EMC7
        return "EMC7";
    } else if (trigger == 81){                          // EGA
        return "EGA";
    } else if (trigger == 85){                          // EG2
        return "EG2";
    } else if (trigger == 83){                          // EG1  
        return "EG1";
    }    
    return "";
}*/



//________________________________ Plotting Invariant mass in a single p_t bin ____________________________________________
void PlotInvMassSinglePtBin(    TH1D* fHistoMappingGGInvMassPtBinPlot, 
                                TH1D* fHistoMappingBackNormInvMassPtBinPlot, 
                                TString namePlot, 
                                TString nameCanvas, 
                                Double_t* fPlottingRangeMeson, 
                                TString decayChannel){
    TCanvas * canvasDataSpectra = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    canvasDataSpectra->SetTopMargin(0.12);
    canvasDataSpectra->SetBottomMargin(0.15);
    canvasDataSpectra->SetRightMargin(0.05);
    canvasDataSpectra->SetLeftMargin(0.15);


    fHistoMappingGGInvMassPtBinPlot->GetXaxis()->SetRangeUser(10, 10000);
    DrawGammaHisto( fHistoMappingGGInvMassPtBinPlot,
                    "2.4 GeV/#it{c} < #it{p}_{T} < 2.6 GeV/#it{c}",
                    Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                    fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);

    DrawGammaHisto( fHistoMappingBackNormInvMassPtBinPlot,
                "2.4 GeV/#it{c} < #it{p}_{T} < 2.6 GeV/#it{c}",
                Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                fPlottingRangeMeson[0],fPlottingRangeMeson[1],1);


    canvasDataSpectra->Print(namePlot.Data());
    delete canvasDataSpectra;
}

//_______________________ Plotting Invariant mass with BG and subtraction in a single p_t bin __________________________________
void PlotExampleInvMassBins(    TH1D* histoInvMassSignalWithBG, 
                                TH1D* histoInvMassSubtracted, 
                                TH1D* histoInvMassBG, 
                                TF1* fitSignal, 
                                Int_t exampleBin, 
                                TString outputDir, 
                                TString suffix, 
                                Double_t* fPlottingRangeMeson, 
                                Float_t* pictDrawingCoordinatesDummy, 
                                Double_t fNumberOfEvents, 
                                TString dateDummy, 
                                TString fMesonType, 
                                TString fSimulation, 
                                TString fPlottingType, 
                                TString fCollisionSystemDummy, 
                                Double_t* fRangeBinsPt, 
                                TString decayChannel,
                                TString detectionChannel = "" ){

    TCanvas* canvasLonelyBin            = new TCanvas("canvasLonelyBin","",200,10,700,500);  // gives the page size
    DrawGammaCanvasSettings( canvasLonelyBin,  0.12, 0.02, 0.07, 0.1);

    Double_t startPt                    = fRangeBinsPt[exampleBin];
    Double_t endPt                      = fRangeBinsPt[exampleBin+1];

    Bool_t kDalitz                      = kFALSE;
    if(decayChannel.CompareTo("e^{+}e^{-}#gamma")==0){
        kDalitz=kTRUE;
    }
    
    canvasLonelyBin->cd();
    DrawGammaHistoBigger( histoInvMassSignalWithBG,
                Form("       %3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);

    DrawGammaHistoBigger( histoInvMassBG,
                Form("       %3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                fPlottingRangeMeson[0],fPlottingRangeMeson[1],1);

    TLegend* legendExamplePlotInvMass   = new TLegend(0.6,0.4,0.9,0.5);
    legendExamplePlotInvMass->SetTextSize(0.04);
    legendExamplePlotInvMass->SetFillColor(0);
    legendExamplePlotInvMass->AddEntry(histoInvMassSignalWithBG,"Data","ep");
    legendExamplePlotInvMass->AddEntry(histoInvMassBG,"combinatorial BG","l");
    legendExamplePlotInvMass->Draw();

    Bool_t lable ;
    Bool_t mcFile;
    if(fMesonType.CompareTo("Pi0") == 0 || fMesonType.CompareTo("Pi0EtaBinning") == 0){
        lable                           = kTRUE;
    } else { 
        lable                           = kFALSE; 
    }
    if(fSimulation.CompareTo("data")==0){ 
        mcFile                          = kFALSE;
    } else {
        mcFile                          = kTRUE;
    }

    if (fPlottingType.CompareTo("thesis") != 0)DrawAliceLogoPi0PerformanceExtract(pictDrawingCoordinatesDummy[0], pictDrawingCoordinatesDummy[1], pictDrawingCoordinatesDummy[2], pictDrawingCoordinatesDummy[3], pictDrawingCoordinatesDummy[4], pictDrawingCoordinatesDummy[5], pictDrawingCoordinatesDummy[6], pictDrawingCoordinatesDummy[7], fNumberOfEvents ,fCollisionSystemDummy ,mcFile, kFALSE, lable ,700 ,500 ,dateDummy,kDalitz, detectionChannel);

    canvasLonelyBin->Update();
    canvasLonelyBin->Print(Form("%s/%s_%s_InvMassDistributionWithBG_PtBin_%i.%s",outputDir.Data(),fMesonType.Data(),fSimulation.Data(),exampleBin,suffix.Data()));

    TCanvas* canvasLonelyBin2           = new TCanvas("canvasLonelyBin2","",200,10,700,500);  // gives the page size
    DrawGammaCanvasSettings( canvasLonelyBin2, 0.12, 0.02, 0.07, 0.1);

    canvasLonelyBin2->cd();
    DrawGammaHistoBigger( histoInvMassSubtracted,
                Form("       %3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                fPlottingRangeMeson[0],fPlottingRangeMeson[1],0    );
    DrawGammaHistoBigger( histoInvMassSubtracted,
                Form("       %3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                fPlottingRangeMeson[0],fPlottingRangeMeson[1],2);
    histoInvMassSubtracted->SetMarkerColor(kRed+2);
    histoInvMassSubtracted->SetLineColor(kRed+2);
    histoInvMassSubtracted->Draw("e1,same");
    if (fitSignal!=0x00){
        fitSignal->SetLineColor(kCyan+3);
        fitSignal->SetLineWidth(2);
        fitSignal->DrawCopy("same");
    }
    TLegend* legendExamplePlotInvMass2;
    if(fMesonType.CompareTo("Pi0") == 0 || fMesonType.CompareTo("Pi0EtaBinning") == 0){
        legendExamplePlotInvMass2       = new TLegend(0.6,0.4,0.95,0.5);
    } else {
        legendExamplePlotInvMass2       = new TLegend(0.6,0.52,0.95,0.62);
    }
    legendExamplePlotInvMass2->SetTextSize(0.033);legendExamplePlotInvMass2->SetFillColor(0);
    legendExamplePlotInvMass2->AddEntry(histoInvMassSubtracted,"Data (BG subtracted)","ep");
    legendExamplePlotInvMass2->AddEntry(fitSignal,"Fit","l");
    legendExamplePlotInvMass2->Draw();

    if (fPlottingType.CompareTo("thesis") != 0)  DrawAliceLogoPi0PerformanceExtract(pictDrawingCoordinatesDummy[0], pictDrawingCoordinatesDummy[1], pictDrawingCoordinatesDummy[2], pictDrawingCoordinatesDummy[3], pictDrawingCoordinatesDummy[4], pictDrawingCoordinatesDummy[5], pictDrawingCoordinatesDummy[6], pictDrawingCoordinatesDummy[7], fNumberOfEvents ,fCollisionSystemDummy, mcFile, kFALSE, lable, 700, 500, dateDummy,kDalitz, detectionChannel);

    canvasLonelyBin2->Update();
    canvasLonelyBin2->Print(Form("%s/%s_%s_InvMassDistributionSubtracted_PtBin_%i.%s",outputDir.Data(),fMesonType.Data(),fSimulation.Data(),exampleBin,suffix.Data()));
}

//_______________________ Plotting Invariant mass with BG and subtraction in a single p_t bin __________________________________
void PlotExampleInvMassBinsV2(  TH1D* histoInvMassSignalWithBG, 
                                TH1D* histoInvMassSubtracted, 
                                TH1D* histoInvMassBG, 
                                TF1* fitSignal, 
                                Int_t exampleBin, 
                                TString outputDir, 
                                TString suffix, 
                                Double_t* fPlottingRangeMeson, 
                                Float_t* pictDrawingCoordinatesDummy, 
                                Double_t fNumberOfEvents, 
                                TString dateDummy, 
                                TString fMesonType, 
                                TString fSimulation, 
                                TString fPlottingType, 
                                TString fCollisionSystemDummy, 
                                Double_t* fRangeBinsPt, 
                                TString decayChannel,
                                TString detectionChannel                = "",
                                Int_t triggerSet                        = 0,
                                Double_t scaleFacSignal                 = 1.0,
                                Int_t detMode                           = 0
                             ){
    
    TString triggerStr              = Form("%s triggered", ((TString)ReturnTriggerName(triggerSet)).Data());
    TString triggerStr2             = ReturnTriggerName(triggerSet);
    TString methodStr               = ReturnTextReconstructionProcess(detMode);
    
    TH1D* histoPi0InvMassSigPlusBG;
    TH1D* histoPi0InvMassSig;
    TH1D* histoPi0InvMassSigRemBG;
    TH1D* histoPi0InvMassSigRemBGSub;
    TH1D* histoPi0InvMassBG;
    TH1D* histoPi0InvMassRemBG;
    TH1D* histoPi0InvMassBGTot;
    TF1* fitPi0InvMassSig;
    TF1* fitPi0InvMassSigRemBG;
    TF1* fitPi0InvMassBG;
    histoPi0InvMassSig               = (TH1D*)histoInvMassSubtracted->Clone("InvMassSig_PtBin07");
    histoPi0InvMassSigRemBG          = (TH1D*)histoInvMassSubtracted->Clone("InvMassSigPlRemBG_PtBin07");
    histoPi0InvMassSigPlusBG         = (TH1D*)histoInvMassSignalWithBG->Clone("InvMassSigPlusBG_PtBin07");
    histoPi0InvMassBG                = (TH1D*)histoInvMassBG->Clone("InvMassBG_PtBin07");
    fitPi0InvMassSig                 = (TF1*)fitSignal->Clone("FitInvMassSig_PtBin07");
    fitPi0InvMassSigRemBG            = (TF1*)fitSignal->Clone("FitInvMassOrig_PtBin07");
    
    histoPi0InvMassSig->Fit(fitPi0InvMassSig,"QRME0");
    for (Int_t l=0; l < 6; l++){
        cout << fitPi0InvMassSig->GetParameter(l) << "\t +- " << fitPi0InvMassSig->GetParError(l) << endl;
    }
    if(fMesonType.CompareTo("Pi0") == 0 || fMesonType.CompareTo("Pi0EtaBinning") == 0){
        fitPi0InvMassBG                                  = new TF1("Linearpp","[0]+[1]*x",0.02,0.25);
    } else {
        fitPi0InvMassBG                                  = new TF1("Linearpp","[0]+[1]*x",0.00,0.3);
    }

    fitPi0InvMassBG->SetParameter(0, fitPi0InvMassSig->GetParameter(4));
    fitPi0InvMassBG->SetParameter(1, fitPi0InvMassSig->GetParameter(5));
    TVirtualFitter * fitter                             = TVirtualFitter::GetFitter();
    Int_t nFreePar                                      = fitPi0InvMassSig->GetNumberFreeParameters();
    double * covMatrix                                  = fitter->GetCovarianceMatrix();

    histoPi0InvMassRemBG                             = (TH1D*)histoPi0InvMassBG->Clone("Pi0_InvMassRemBG_Example");
    for (Int_t j = 1; j < histoPi0InvMassRemBG->GetNbinsX()+1; j++){
        histoPi0InvMassRemBG->SetBinContent(j,0);
        histoPi0InvMassRemBG->SetBinError(j,0);
    }
    if(fMesonType.CompareTo("Pi0") == 0 || fMesonType.CompareTo("Pi0EtaBinning") == 0){
        for (Int_t j = histoPi0InvMassSig->GetXaxis()->FindBin(0.01); j < histoPi0InvMassSig->GetXaxis()->FindBin(0.30)+1; j++){
            Double_t startBinEdge                                   = histoPi0InvMassSig->GetXaxis()->GetBinLowEdge(j);
            Double_t endBinEdge                                     = histoPi0InvMassSig->GetXaxis()->GetBinUpEdge(j);
            Double_t intLinearBack                                  = fitPi0InvMassBG->Integral(startBinEdge, endBinEdge)/(endBinEdge-startBinEdge) ;
            Double_t errorLinearBck                                 = pow(( pow( (endBinEdge-startBinEdge)*fitPi0InvMassSig->GetParError(4),2) +
            pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fitPi0InvMassSig->GetParError(5),2)
            +2*covMatrix[nFreePar*nFreePar-2]*(endBinEdge-startBinEdge)*0.5*
            (endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5)/(endBinEdge-startBinEdge);
            histoPi0InvMassRemBG->SetBinContent(j,intLinearBack);
            histoPi0InvMassRemBG->SetBinError(j,errorLinearBck);
        }
    } else {
        for (Int_t j = histoPi0InvMassSig->GetXaxis()->FindBin(0.30); j < histoPi0InvMassSig->GetXaxis()->FindBin(0.70)+1; j++){
            Double_t startBinEdge                                   = histoPi0InvMassSig->GetXaxis()->GetBinLowEdge(j);
            Double_t endBinEdge                                     = histoPi0InvMassSig->GetXaxis()->GetBinUpEdge(j);
            Double_t intLinearBack                                  = fitPi0InvMassBG->Integral(startBinEdge, endBinEdge)/(endBinEdge-startBinEdge) ;
            Double_t errorLinearBck                                 = pow(( pow( (endBinEdge-startBinEdge)*fitPi0InvMassSig->GetParError(4),2) +
            pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fitPi0InvMassSig->GetParError(5),2)
            +2*covMatrix[nFreePar*nFreePar-2]*(endBinEdge-startBinEdge)*0.5*
            (endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5)/(endBinEdge-startBinEdge);
            histoPi0InvMassRemBG->SetBinContent(j,intLinearBack);
            histoPi0InvMassRemBG->SetBinError(j,errorLinearBck);
        }
    }

    histoPi0InvMassBGTot         = (TH1D*)histoPi0InvMassBG->Clone("Pi0_InvMassTotBG_Example");
    histoPi0InvMassBGTot->Sumw2();
    histoPi0InvMassBGTot->Add(histoPi0InvMassRemBG);
    histoPi0InvMassSigRemBGSub   = (TH1D*)histoPi0InvMassSig->Clone("Pi0_InvMassSigRemBGSub_Example");
    histoPi0InvMassSigRemBGSub->Sumw2();
    histoPi0InvMassSigRemBGSub->Add(histoPi0InvMassRemBG,-1);

    fitPi0InvMassSig->SetParameter(4, 0);
    fitPi0InvMassSig->SetParameter(5, 0);
    histoPi0InvMassSigRemBGSub->Scale(scaleFacSignal);
    histoPi0InvMassSigRemBG->Scale(scaleFacSignal);
    
    Double_t textSizeLabelsPixel                 = 100*3/5;
    TCanvas* canvasInvMassSamplePlot    = new TCanvas("canvasInvMassSamplePlotNew","",0,0,1500,1500);  // gives the page size
    DrawGammaCanvasSettings( canvasInvMassSamplePlot,  0.09, 0.012, 0.035, 0.08);

    Double_t startPt                    = fRangeBinsPt[exampleBin];
    Double_t endPt                      = fRangeBinsPt[exampleBin+1];

    Style_t markerStyleInvMassSGBG      = 0;
    Size_t markerSizeInvMassSGBG        = 0;
    Color_t markerColorInvMassSGBG      = kBlack;
    Style_t markerStyleInvMassMBG       = 24;
    Size_t markerSizeInvMassMBG         = 1.5;
    Color_t markerColorInvMassMBG       = kGray+2;
    Color_t markerColorInvMassMBG1      = kGray+3;
    Color_t markerColorInvMassMBG2      = kGray+1;
    Style_t markerStyleInvMassBG        = 20;
    Size_t markerSizeInvMassBG          = 2;
    Color_t markerColorInvMassBG        = kBlack;
    Style_t markerStyleInvMassSG        = 20;
    Size_t markerSizeInvMassSG          = 3;
    Color_t markerColorInvMassSG        = kRed+2;
    Color_t fitColorInvMassSG           = kAzure+2;

    Double_t textsizeLabelsPP       = 0.04;
    Double_t marginInvMass          = 0.1*1500;
    Double_t textsizeLabelsInvMass  = 0;
    Double_t textsizeFacInvMass     = 0;
    if (canvasInvMassSamplePlot->XtoPixel(canvasInvMassSamplePlot->GetX2()) < canvasInvMassSamplePlot->YtoPixel(canvasInvMassSamplePlot->GetY1())){
        textsizeLabelsInvMass       = (Double_t)textSizeLabelsPixel/canvasInvMassSamplePlot->XtoPixel(canvasInvMassSamplePlot->GetX2()) ;
        textsizeFacInvMass          = (Double_t)1./canvasInvMassSamplePlot->XtoPixel(canvasInvMassSamplePlot->GetX2()) ;
    } else {
        textsizeLabelsInvMass       = (Double_t)textSizeLabelsPixel/canvasInvMassSamplePlot->YtoPixel(canvasInvMassSamplePlot->GetY1());
        textsizeFacInvMass          = (Double_t)1./canvasInvMassSamplePlot->YtoPixel(canvasInvMassSamplePlot->GetY1());
    }

    TH1F * histo1DInvMassDummy;
    if(fMesonType.CompareTo("Pi0") == 0 || fMesonType.CompareTo("Pi0EtaBinning") == 0){
        histo1DInvMassDummy             = new TH1F("histo1DInvMass2","histo1DInvMass2",11000,0.02,0.255);
        SetStyleHistoTH1ForGraphs(histo1DInvMassDummy, Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()),"Counts",0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,
                                0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,0.88, 0.115/(textsizeFacInvMass*marginInvMass));
    } else {
        histo1DInvMassDummy             = new TH1F("histo1DInvMass2","histo1DInvMass2",11000,0.35,0.695);
        SetStyleHistoTH1ForGraphs(histo1DInvMassDummy, Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()),"Counts",0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,
                                0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,0.88, 0.115/(textsizeFacInvMass*marginInvMass));
    }


    // Set range for fits and labels
    TLatex *labelInvMassPtRange;
    if(fMesonType.CompareTo("Pi0") == 0 || fMesonType.CompareTo("Pi0EtaBinning") == 0){
        labelInvMassPtRange = new TLatex(0.95,0.9, Form("#pi^{0}: %3.1f GeV/#it{c} < #it{p}_{T} < %3.1f GeV/#it{c}",startPt,endPt));
        fitPi0InvMassSig->SetRange(0,0.255);
        fitPi0InvMassSigRemBG->SetRange(0,0.255);
    } else {
        labelInvMassPtRange = new TLatex(0.95,0.9, Form("#eta: %3.1f GeV/#it{c} < #it{p}_{T} < %3.1f GeV/#it{c}",startPt,endPt));
        fitPi0InvMassSig->SetRange(0.35,0.695);
        fitPi0InvMassSigRemBG->SetRange(0.35,0.695);
    }
    // Set fit colors
    fitPi0InvMassSig->SetNpx(10000);
    fitPi0InvMassSig->SetLineColor(fitColorInvMassSG);
    fitPi0InvMassSig->SetLineWidth(4);
    fitPi0InvMassSigRemBG->SetNpx(10000);
    fitPi0InvMassSigRemBG->SetLineColor(fitColorInvMassSG);
    fitPi0InvMassSigRemBG->SetLineWidth(4);
        
    TH1D* histoFit  = (TH1D*)fitPi0InvMassSig->GetHistogram();
    histoFit->SetTitle("");
    histoFit->Scale(scaleFacSignal);
    histoFit->SetLineWidth(4);
    TH1D* histoFitWBG  = (TH1D*)fitPi0InvMassSigRemBG->GetHistogram();
    histoFitWBG->SetTitle("");
    histoFitWBG->Scale(scaleFacSignal);
    histoFitWBG->SetLineWidth(4);

    canvasInvMassSamplePlot->cd();
    histo1DInvMassDummy->GetYaxis()->SetRangeUser(histoPi0InvMassSigRemBGSub->GetMinimum(),1.15*histoPi0InvMassSigPlusBG->GetMaximum());
    histo1DInvMassDummy->Draw("AXIS");
    
    DrawGammaSetMarker(histoPi0InvMassSigPlusBG, markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
    histoPi0InvMassSigPlusBG->SetLineWidth(1);
    histoPi0InvMassSigPlusBG->Draw("hist,e,same");
    DrawGammaSetMarker(histoPi0InvMassBGTot, markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG, markerColorInvMassMBG);
    histoPi0InvMassBGTot->Draw("same");

    Int_t nLegendLines      = 5;
    if (scaleFacSignal == 1.0){
        DrawGammaSetMarker(histoPi0InvMassSigRemBGSub, markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
        histoPi0InvMassSigRemBGSub->Draw("same");    
        fitPi0InvMassSig->Draw("same");
    } else {
        DrawGammaSetMarker(histoPi0InvMassSigRemBGSub, markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
        histoPi0InvMassSigRemBGSub->Draw("same");
        
        histoFit->SetLineColor(fitColorInvMassSG);
        histoFit->SetLineWidth(4);
        histoFit->Draw("c,same");
        nLegendLines++;
    }

    TLatex *labelALICE      = new TLatex(0.135,0.9,"ALICE performance");
    SetStyleTLatex( labelALICE, 0.85*textSizeLabelsPixel,4);
    labelALICE->SetTextFont(43);
    labelALICE->Draw();
    
    TLatex *labelInvMassEnergy      = new TLatex(0.135,0.9-0.9*textsizeLabelsPP,fCollisionSystemDummy.Data());
    SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4);
    labelInvMassEnergy->SetTextFont(43);
    labelInvMassEnergy->Draw();

    TLatex *labelTrigger  = new TLatex(0.135,0.9-2*0.9*textsizeLabelsPP,triggerStr.Data());
    SetStyleTLatex( labelTrigger, 0.85*textSizeLabelsPixel,4);
    labelTrigger->SetTextFont(43);
    labelTrigger->Draw();
    
    TLatex *labelInvMassReco  = new TLatex(0.135,0.9-3*0.9*textsizeLabelsPP, methodStr);
    SetStyleTLatex( labelInvMassReco, 0.85*textSizeLabelsPixel,4);
    labelInvMassReco->SetTextFont(43);
    labelInvMassReco->Draw();

    SetStyleTLatex( labelInvMassPtRange, 0.85*textSizeLabelsPixel,4);
    labelInvMassPtRange->SetTextAlign(31);
    labelInvMassPtRange->SetTextFont(43);
    labelInvMassPtRange->Draw();

    TLegend* legendInvMass  = GetAndSetLegend2(0.67, 0.87-nLegendLines*0.75*textsizeLabelsPP, 0.9, 0.87, 0.85*textSizeLabelsPixel);
    legendInvMass->SetMargin(0.25);
    legendInvMass->AddEntry(histoPi0InvMassSigPlusBG,"Raw real events","l");
    legendInvMass->AddEntry(histoPi0InvMassBGTot,"Mixed event +","p");
    legendInvMass->AddEntry((TObject*)0,"rem. BG","");
    legendInvMass->AddEntry(histoPi0InvMassSigRemBGSub,"BG subtracted","p");
    if (scaleFacSignal != 1.0){
        legendInvMass->AddEntry((TObject*)0,Form("scaled by %2.1f",scaleFacSignal),"");
    }    
    legendInvMass->AddEntry(fitPi0InvMassSig, "Fit","l");
    legendInvMass->Draw();

    canvasInvMassSamplePlot->SaveAs(Form("%s/%s_%s_InvMassBin%s_%s.%s",outputDir.Data(),fMesonType.Data(),fSimulation.Data(), methodStr.Data(), triggerStr2.Data(), suffix.Data()));
    
    canvasInvMassSamplePlot->cd();
    histo1DInvMassDummy->Draw();
    histo1DInvMassDummy->GetYaxis()->SetRangeUser(histoPi0InvMassSigRemBGSub->GetMinimum(),1.15*histoPi0InvMassSigPlusBG->GetMaximum());
    histo1DInvMassDummy->Draw("AXIS");
    
    DrawGammaSetMarker(histoPi0InvMassSigPlusBG, markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
    histoPi0InvMassSigPlusBG->SetLineWidth(1);
    histoPi0InvMassSigPlusBG->Draw("hist,e,same");
    DrawGammaSetMarker(histoPi0InvMassBG, markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG1, markerColorInvMassMBG1);
    histoPi0InvMassBG->Draw("same");
    DrawGammaSetMarker(histoPi0InvMassRemBG, markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG2, markerColorInvMassMBG2);
    histoPi0InvMassRemBG->Draw("same");

    
    if (scaleFacSignal == 1.0){
        histoPi0InvMassSigRemBGSub->Draw("same");    
        fitPi0InvMassSig->Draw("same");
    } else {
        histoPi0InvMassSigRemBGSub->Draw("same");    
        histoFit->Draw("c,same");
    }

    labelALICE->Draw();
    labelInvMassEnergy->Draw();
    labelTrigger->Draw();
    labelInvMassReco->Draw();
    labelInvMassPtRange->Draw();

    TLegend* legendInvMass2  = GetAndSetLegend2(0.67, 0.87-nLegendLines*0.75*textsizeLabelsPP, 0.9, 0.87, 0.85*textSizeLabelsPixel);
    legendInvMass2->SetMargin(0.25);
    legendInvMass2->AddEntry(histoPi0InvMassSigPlusBG,"Raw real events","l");
    legendInvMass2->AddEntry(histoPi0InvMassBG,"Mixed event BG","p");
    legendInvMass2->AddEntry(histoPi0InvMassRemBG,"Rem. BG","p");
    legendInvMass2->AddEntry(histoPi0InvMassSigRemBGSub,"BG subtracted","p");
    if (scaleFacSignal != 1.0){
        legendInvMass2->AddEntry((TObject*)0,Form("scaled by %2.1f",scaleFacSignal),"");
    }    
    legendInvMass2->AddEntry(fitPi0InvMassSig, "Fit","l");
    legendInvMass2->Draw();

    canvasInvMassSamplePlot->SaveAs(Form("%s/%s_%s_InvMassBinBGFurtherSplit%s_%s.%s",outputDir.Data(),fMesonType.Data(),fSimulation.Data(), methodStr.Data(), triggerStr2.Data(),  suffix.Data()));

    // Plot Invariant mass sample bin with linear BG in Signal
    canvasInvMassSamplePlot->cd();
    histo1DInvMassDummy->Draw();
    histo1DInvMassDummy->GetYaxis()->SetRangeUser(histoPi0InvMassSigRemBG->GetMinimum(),1.15*histoPi0InvMassSigPlusBG->GetMaximum());
    histo1DInvMassDummy->Draw("AXIS");
    
    DrawGammaSetMarker(histoPi0InvMassSigPlusBG, markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
    histoPi0InvMassSigPlusBG->SetLineWidth(1);
    histoPi0InvMassSigPlusBG->Draw("hist,e,same");
    DrawGammaSetMarker(histoPi0InvMassBG, markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG, markerColorInvMassMBG);
    histoPi0InvMassBG->Draw("same");
    
    if (scaleFacSignal == 1.0){
        DrawGammaSetMarker(histoPi0InvMassSigRemBG, markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
        histoPi0InvMassSigRemBG->Draw("same");    
        fitPi0InvMassSigRemBG->Draw("same");
    } else {
        DrawGammaSetMarker(histoPi0InvMassSigRemBG, markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
        histoPi0InvMassSigRemBG->Draw("same");
        
        histoFitWBG->SetLineColor(fitColorInvMassSG);
        histoFitWBG->SetLineWidth(4);
        histoFitWBG->Draw("c,same");
    }

    labelALICE->Draw();
    labelInvMassEnergy->Draw();
    labelTrigger->Draw();
    labelInvMassReco->Draw();
    labelInvMassPtRange->Draw();

    TLegend* legendInvMass3  = GetAndSetLegend2(0.62, 0.87-nLegendLines*0.75*textsizeLabelsPP, 0.9, 0.87, 0.85*textSizeLabelsPixel);
    legendInvMass3->SetMargin(0.22);
    legendInvMass3->AddEntry(histoPi0InvMassSigPlusBG,"Raw real events","l");
    legendInvMass3->AddEntry(histoPi0InvMassBG,"Mixed event BG","p");
    legendInvMass3->AddEntry(histoPi0InvMassSigRemBG,"Mixed evt. BG sub.","p");
    if (scaleFacSignal != 1.0){
        legendInvMass3->AddEntry((TObject*)0,Form("scaled by %2.1f",scaleFacSignal),"");
    }    
    legendInvMass3->AddEntry(fitPi0InvMassSigRemBG, "Signal fit +","l");
    legendInvMass3->AddEntry((TObject*)0,"linear BG fit","");
    legendInvMass3->Draw();

    canvasInvMassSamplePlot->SaveAs(Form("%s/%s_%s_InvMassBinBGInFit%s_%s.%s",outputDir.Data(),fMesonType.Data(),fSimulation.Data(), methodStr.Data(), triggerStr2.Data(),  suffix.Data()));
    
}


//_______________________ Plotting Invariant mass with BG and subtraction in a single p_t bin __________________________________
void PlotExampleInvMassBinsMC(  TH1D* fHistoTrueSignal, 
                                TH1D* fHistoTrueSignalPhotons, 
                                TH1D* fHistoTrueSignalElectrons, 
                                TH1D* fHistoTrueSignalConvPhotons, 
                                TH1D* fHistoTrueSignalMixed, 
                                Int_t exampleBin, 
                                TString outputDir, 
                                TString suffix, 
                                Double_t* fPlottingRangeMeson, 
                                Float_t* pictDrawingCoordinatesDummy, 
                                Double_t fNumberOfEvents, 
                                TString dateDummy, 
                                TString fMesonType, 
                                TString fSimulation, 
                                TString fPlottingType, 
                                TString fCollisionSystemDummy, 
                                Double_t* fRangeBinsPt, 
                                TString decayChannel,
                                TString detectionChannel = "",
                                Int_t mode = 0 ){

    Double_t startPt                = fRangeBinsPt[exampleBin];
    Double_t endPt                  = fRangeBinsPt[exampleBin+1];

    Bool_t kDalitz                  = kFALSE;
    if(decayChannel.CompareTo("e^{+}e^{-}#gamma")==0){
        kDalitz                     = kTRUE;
    }

    Bool_t lable ;
    Bool_t mcFile;
    if(fMesonType.CompareTo("Pi0") == 0 || fMesonType.CompareTo("Pi0EtaBinning") == 0){ 
        lable                       = kTRUE;
    } else { 
        lable                       = kFALSE; 
    }
    if(fSimulation.CompareTo("data")==0){ 
        mcFile                      = kFALSE;
    } else {
        mcFile                      = kTRUE;
    }

    TCanvas* canvasLonelyBin2       = new TCanvas("canvasLonelyBin2","",200,10,700,500);  // gives the page size
    DrawGammaCanvasSettings( canvasLonelyBin2, 0.12, 0.02, 0.07, 0.1);

    canvasLonelyBin2->cd();
        DrawGammaHistoBigger(   fHistoTrueSignal,
                                Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
        DrawGammaHistoBigger(   fHistoTrueSignal,
                                Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                fPlottingRangeMeson[0],fPlottingRangeMeson[1],2);

        fHistoTrueSignalPhotons->SetMarkerColor(kRed+2);
        fHistoTrueSignalPhotons->SetMarkerStyle(21);
        fHistoTrueSignalPhotons->SetMarkerSize(0.7);
        fHistoTrueSignalPhotons->SetLineColor(kRed+2);
        fHistoTrueSignalPhotons->DrawCopy("same,p,e1");
        
        if (fHistoTrueSignalElectrons != NULL){
            fHistoTrueSignalElectrons->SetMarkerColor(kBlue+2);
            fHistoTrueSignalElectrons->SetMarkerStyle(33);
            fHistoTrueSignalElectrons->SetMarkerSize(0.9);
            fHistoTrueSignalElectrons->SetLineColor(kBlue+2);
            fHistoTrueSignalElectrons->DrawCopy("same,p,e1");
        }
        
        fHistoTrueSignalConvPhotons->SetMarkerColor(kCyan+2);
        fHistoTrueSignalConvPhotons->SetMarkerStyle(27);
        fHistoTrueSignalConvPhotons->SetMarkerSize(0.7);
        fHistoTrueSignalConvPhotons->SetLineColor(kCyan+2);
        fHistoTrueSignalConvPhotons->DrawCopy("same,p,e1");
        
        if (fHistoTrueSignalMixed != NULL){
            fHistoTrueSignalMixed->SetMarkerColor(kViolet+2);
            fHistoTrueSignalMixed->SetMarkerStyle(24);
            fHistoTrueSignalMixed->SetMarkerSize(0.7);
            fHistoTrueSignalMixed->SetLineColor(kViolet+2);
            fHistoTrueSignalMixed->DrawCopy("same,p,e1");
        }

        Double_t nSignals           = 3;
        if (fHistoTrueSignalMixed != NULL) nSignals++;
        if (fHistoTrueSignalElectrons != NULL) nSignals++;
        TLegend* legendMC           = new TLegend(0.7,0.5-(nSignals*0.033),0.95,0.5);
        legendMC->SetTextSize(0.033);
        legendMC->SetTextFont(62);
        legendMC->SetLineColor(0);
        legendMC->SetLineWidth(0);
        legendMC->SetFillStyle(0);
        legendMC->SetFillColor(0);
        legendMC->SetMargin(0.15);
        legendMC->AddEntry(fHistoTrueSignal,"validated meson","ep");
        if (mode == 4 || mode == 5){
            legendMC->AddEntry(fHistoTrueSignalPhotons,"val. #gamma#gamma","ep");
            if (fHistoTrueSignalElectrons != NULL)legendMC->AddEntry(fHistoTrueSignalElectrons,"val. e^{#pm}e^{#pm}","ep");
            legendMC->AddEntry(fHistoTrueSignalConvPhotons,"val. #gamma_{conv}#gamma_{conv}","ep");
            if (fHistoTrueSignalMixed != NULL) legendMC->AddEntry(fHistoTrueSignalMixed,"val. #gamma#gamma_{conv}","ep");
        } else if (mode == 2 || mode == 3) {
            legendMC->AddEntry(fHistoTrueSignalPhotons,"val. #gamma_{conv}#gamma","ep");
            if (fHistoTrueSignalElectrons != NULL)legendMC->AddEntry(fHistoTrueSignalElectrons,"val. #gamma_{conv}e^{#pm}","ep");
            legendMC->AddEntry(fHistoTrueSignalConvPhotons,"val. #gamma_{conv}#gamma_{conv}","ep");
        }
        legendMC->Draw();

    if (fPlottingType.CompareTo("thesis") != 0)  DrawAliceLogoPi0PerformanceExtract(pictDrawingCoordinatesDummy[0], pictDrawingCoordinatesDummy[1], pictDrawingCoordinatesDummy[2], pictDrawingCoordinatesDummy[3], pictDrawingCoordinatesDummy[4], pictDrawingCoordinatesDummy[5], pictDrawingCoordinatesDummy[6], pictDrawingCoordinatesDummy[7], fNumberOfEvents ,fCollisionSystemDummy, mcFile, kFALSE, lable, 700, 500, dateDummy,kDalitz, detectionChannel);

    canvasLonelyBin2->Update();
    canvasLonelyBin2->Print(Form("%s/%s_%s_TrueInvMassDistributionDisentangled_PtBin_%i.%s",outputDir.Data(),fMesonType.Data(),fSimulation.Data(),exampleBin,suffix.Data()));
}




//__________________________________________ Plotting all Invariant Mass bins _______________________________________________
void PlotInvMassInPtBins(   TH1D** fHistoMappingGGInvMassPtBinPlot, 
                            TH1D** fHistoMappingBackNormInvMassPtBinPlot, 
                            TString namePlot, 
                            TString nameCanvas, 
                            TString namePad, 
                            Double_t* fPlottingRangeMeson, 
                            TString dateDummy, 
                            TString fMesonType, 
                            Int_t fRowPlot, 
                            Int_t fColumnPlot, 
                            Int_t fStartBinPtRange, 
                            Int_t fNumberPtBins, 
                            Double_t* fRangeBinsPt, 
                            TString fDecayChannel, 
                            Bool_t fMonteCarloInfo, 
                            TString decayChannel,  
                            TString fDetectionChannel, 
                            TString fEnergy){
    TGaxis::SetMaxDigits(3);

    TCanvas * canvasDataSpectra     = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    canvasDataSpectra->SetTopMargin(0.0);
    canvasDataSpectra->SetBottomMargin(0.0);
    canvasDataSpectra->SetRightMargin(0.0);
    canvasDataSpectra->SetLeftMargin(0.0);

    TPad * padDataSpectra           = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    padDataSpectra->SetFillColor(0);
    padDataSpectra->GetFrame()->SetFillColor(0);
    padDataSpectra->SetBorderMode(0);
    padDataSpectra->SetLogy(0);
    padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
    padDataSpectra->Draw();

    cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
    cout << fMesonType.Data() << endl;
    Int_t place                     = 0;
    for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
        cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
        Double_t startPt            = fRangeBinsPt[iPt];
        Double_t endPt              = fRangeBinsPt[iPt+1];

        place                       = place + 1; //give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataSpectra->cd(place);

            cout << "entered ALICE plotting" << endl;
            TString textAlice       = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo){
                textEvents          = "MC";
            } else {
                textEvents          = "Data";
            }
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
            TLatex *latexDate       = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy          = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
            TLatex *process         = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess      = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
            TLatex *events          = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

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

            TLegend* legendData     = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+2.)*differenceText);
            legendData->SetTextSize(textHeight);
            legendData->SetTextFont(62);
            legendData->SetFillColor(0);
            legendData->SetFillStyle(0);
            legendData->SetLineWidth(0);
            legendData->SetLineColor(0);
            legendData->SetMargin(0.15);
            Size_t markersize       = fHistoMappingGGInvMassPtBinPlot[iPt]->GetMarkerSize();
            fHistoMappingGGInvMassPtBinPlot[iPt]->SetMarkerSize(3*markersize);
            legendData->AddEntry(fHistoMappingGGInvMassPtBinPlot[iPt],"same evt. #it{M}_{#gamma#gamma} (BG+Signal)","ep");
            Size_t linesize         = fHistoMappingBackNormInvMassPtBinPlot[iPt]->GetLineWidth();
            fHistoMappingBackNormInvMassPtBinPlot[iPt]->SetLineWidth(5*linesize);
            legendData->AddEntry(fHistoMappingBackNormInvMassPtBinPlot[iPt],"mixed evt. #it{M}_{#gamma#gamma}","l");
            legendData->Draw();
        } else {
            
            padDataSpectra->cd(place);
            padDataSpectra->cd(place)->SetTopMargin(0.12);
            padDataSpectra->cd(place)->SetBottomMargin(0.15);
            padDataSpectra->cd(place)->SetRightMargin(0.02);
            cout << "place: " << place << endl;
            int remaining           = (place-1)%fColumnPlot;
            if (remaining > 0) padDataSpectra->cd(place)->SetLeftMargin(0.15);
            else padDataSpectra->cd(place)->SetLeftMargin(0.25);
            cout << "here" << endl;
            DrawGammaHisto( fHistoMappingGGInvMassPtBinPlot[iPt],
                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
            cout << "here" << endl;
            DrawGammaHisto( fHistoMappingBackNormInvMassPtBinPlot[iPt],
                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],1);
            cout << "here" << endl;
            Double_t fBGFitRangeLow     = fBGFitRange[0];
            Double_t fBGFitRangeHigh    = fBGFitRange[1];
            if (namePlot.Contains("Left")){
                fBGFitRangeLow          = fBGFitRangeLeft[0];
                fBGFitRangeHigh         = fBGFitRangeLeft[1];
            }  
            TBox *box               = new TBox(fBGFitRangeLow,fHistoMappingGGInvMassPtBinPlot[iPt]->GetMaximum()*0.93,fBGFitRangeHigh,fHistoMappingGGInvMassPtBinPlot[iPt]->GetMaximum()*0.91);
            box->SetFillStyle(1001);
            box->SetFillColor(kAzure+9);
            box->Draw("same");
        }
    }
    canvasDataSpectra->Print(namePlot.Data());
    delete padDataSpectra;
    delete canvasDataSpectra;
}

//________________________________________ Plot Invariant Mass Bin With Secondary Contribution _______________________________
void PlotInvMassSecondaryInPtBins(  TH1D** fHistoMappingGGInvMassPtBinPlot, 
                                    TH1D** fHistoMappingSecondaryTotalInvMassPtBinPlot, 
                                    TH1D** fHistoMappingSecondaryK0sInvMassPtBinPlot, 
                                    TString namePlot, 
                                    TString nameCanvas, 
                                    TString namePad, 
                                    Double_t* fPlottingRangeMeson, 
                                    TString dateDummy, 
                                    TString fMesonType,  
                                    Int_t fRowPlot, 
                                    Int_t fColumnPlot, 
                                    Int_t fStartBinPtRange, 
                                    Int_t fNumberPtBins, 
                                    Double_t* fRangeBinsPt, 
                                    TString fDecayChannel, 
                                    Bool_t fMonteCarloInfo, 
                                    TString decayChannel,  
                                    TString fDetectionChannel, 
                                    TString fEnergy){
    TGaxis::SetMaxDigits(3);

    TCanvas * canvasDataSpectra         = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    canvasDataSpectra->SetTopMargin(0.0);
    canvasDataSpectra->SetBottomMargin(0.0);
    canvasDataSpectra->SetRightMargin(0.0);
    canvasDataSpectra->SetLeftMargin(0.0);

    TPad * padDataSpectra               = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    padDataSpectra->SetFillColor(0);
    padDataSpectra->GetFrame()->SetFillColor(0);
    padDataSpectra->SetBorderMode(0);
    padDataSpectra->SetLogy(0);
    padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
    padDataSpectra->Draw();

    cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
    cout << fMesonType.Data() << endl;
    Int_t place                         = 0;
    for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
        cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
        Double_t startPt                = fRangeBinsPt[iPt];
        Double_t endPt                  = fRangeBinsPt[iPt+1];

        place                           = place + 1; //give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataSpectra->cd(place);

            TString textAlice           = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo){
                textEvents              = "MC";   
            } else {
                textEvents              = "Data";
            }
            Double_t nPixels            = 13;
            Double_t textHeight         = 0.08;
            if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
                textHeight              = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
            } else {
                textHeight              = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
            }
            Double_t startTextX         = 0.10;
            Double_t startTextY         = 0.8;
            Double_t differenceText     = textHeight*1.25;

            TLatex *alice               = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
            TLatex *latexDate           = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy              = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
            TLatex *process             = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess          = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
            TLatex *events              = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

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
            int remaining = (place-1)%fColumnPlot;
            if (remaining > 0) padDataSpectra->cd(place)->SetLeftMargin(0.15);
            else padDataSpectra->cd(place)->SetLeftMargin(0.25);

            DrawGammaHisto( fHistoMappingGGInvMassPtBinPlot[iPt],
                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);

            DrawGammaHisto( fHistoMappingSecondaryTotalInvMassPtBinPlot[iPt],
                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],1);
            
            
            if (fHistoMappingSecondaryK0sInvMassPtBinPlot!=NULL){
                if (fHistoMappingSecondaryK0sInvMassPtBinPlot[iPt]!=NULL){
                    DrawGammaHistoColored( fHistoMappingSecondaryK0sInvMassPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0,2);
                }
            }
        }
    }
    canvasDataSpectra->Print(namePlot.Data());
    delete padDataSpectra;
    delete canvasDataSpectra;

}

//________________________________________ Plot Invariant Mass Bin With GG contamination for Dalitz _______________________________
void PlotInvMassBckGGInPtBins(  TH1D** fHistoMappingGGInvMassPtBinPlot, 
                                TH1D** fHistoMappingSecondaryTotalInvMassPtBinPlot,
                                TString namePlot, 
                                TString nameCanvas, 
                                TString namePad, 
                                Double_t* fPlottingRangeMeson, 
                                TString dateDummy, 
                                TString fMesonType, 
                                Int_t fRowPlot, 
                                Int_t fColumnPlot, 
                                Int_t fStartBinPtRange,
                                Int_t fNumberPtBins, 
                                Double_t* fRangeBinsPt, 
                                TString fDecayChannel, 
                                Bool_t fMonteCarloInfo, 
                                TString decayChannel, 
                                TString fDetectionChannel, 
                                TString fEnergy){
    TGaxis::SetMaxDigits(3);

    TCanvas * canvasDataSpectra     = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    canvasDataSpectra->SetTopMargin(0.00);
    canvasDataSpectra->SetBottomMargin(0.00);
    canvasDataSpectra->SetRightMargin(0.0);
    canvasDataSpectra->SetLeftMargin(0.00);

    TPad * padDataSpectra           = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    padDataSpectra->SetFillColor(0);
    padDataSpectra->GetFrame()->SetFillColor(0);
    padDataSpectra->SetBorderMode(0);
    padDataSpectra->SetLogy(0);
    padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
    padDataSpectra->Draw();

    cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
    cout << fMesonType.Data() << endl;

    Int_t place                     = 0;
    for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
        cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
        Double_t startPt            = fRangeBinsPt[iPt];
        Double_t endPt              = fRangeBinsPt[iPt+1];

        place                       = place + 1; //give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataSpectra->cd(place);

            TString textAlice       = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo){
                textEvents          = "MC";
            } else {
                textEvents          = "Data";
            }
            Double_t nPixels        = 13;
            Double_t textHeight     = 0.08;
            if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
                textHeight          = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
            } else {
                textHeight          = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
            }
            Double_t startTextX     = 0.10;
            Double_t startTextY     = 0.8;
            Double_t differenceText = textHeight*1.25;

            TLatex *alice           = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
            TLatex *latexDate       = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy          = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
            TLatex *process         = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess      = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
            TLatex *events          = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

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
            int remaining           = (place-1)%fColumnPlot;
            if (remaining > 0) padDataSpectra->cd(place)->SetLeftMargin(0.15);
            else padDataSpectra->cd(place)->SetLeftMargin(0.25);

            DrawGammaHisto( fHistoMappingGGInvMassPtBinPlot[iPt],
                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);

            DrawGammaHisto( fHistoMappingSecondaryTotalInvMassPtBinPlot[iPt],
                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],1);
            }
    }
    canvasDataSpectra->Print(namePlot.Data());
    delete padDataSpectra;
    delete canvasDataSpectra;

}


//__________________________________ Plotting Invariant Mass as Ratio _________________________________________________________
void PlotInvMassRatioInPtBins(  TH1D** fHistoMappingGGInvMassPtBinPlot, 
                                TString namePlot, 
                                TString nameCanvas, 
                                TString namePad, 
                                Double_t* fPlottingRangeMeson, 
                                TString dateDummy, 
                                TString fMesonType, 
                                Int_t fRowPlot,
                                Int_t fColumnPlot,
                                Int_t fStartBinPtRange,
                                Int_t fNumberPtBins,
                                Double_t* fRangeBinsPt, 
                                TString fDecayChannel,
                                Bool_t fMonteCarloInfo,
                                TString decayChannel,
                                TString fDetectionChannel, 
                                TString fEnergy){
    TCanvas * canvasDataSpectra     = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    canvasDataSpectra->SetTopMargin(0.00);
    canvasDataSpectra->SetBottomMargin(0.00);
    canvasDataSpectra->SetRightMargin(0.00);
    canvasDataSpectra->SetLeftMargin(0.00);

    TPad * padDataSpectra           = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    padDataSpectra->SetFillColor(0);
    padDataSpectra->GetFrame()->SetFillColor(0);
    padDataSpectra->SetBorderMode(0);
    padDataSpectra->SetLogy(1);
    padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
    padDataSpectra->Draw();

    cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
    cout << fMesonType.Data() << endl;
    Int_t place                     = 0;
    for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
        cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
        Double_t startPt            = fRangeBinsPt[iPt];
        Double_t endPt              = fRangeBinsPt[iPt+1];
        place                       = place + 1;//give the right place in the page

        if (place == fColumnPlot){
            iPt--;
            padDataSpectra->cd(place);

            TString textAlice       = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo){
                textEvents          = "MC";
            } else { 
                textEvents          = "Data";
            }
            Double_t nPixels        = 13;
            Double_t textHeight     = 0.08;
            if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
                textHeight          = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
            } else {
                textHeight          = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
            }
            Double_t startTextX     = 0.10;
            Double_t startTextY     = 0.8;
            Double_t differenceText = textHeight*1.25;

            TLatex *alice           = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
            TLatex *latexDate       = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy          = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
            TLatex *process         = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess      = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
            TLatex *events          = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

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
            int remaining           = (place-1)%fColumnPlot;
            if (remaining > 0) padDataSpectra->cd(place)->SetLeftMargin(0.15);
            else padDataSpectra->cd(place)->SetLeftMargin(0.25);

            DrawGammaHisto( fHistoMappingGGInvMassPtBinPlot[iPt],
                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);

            DrawGammaHisto( fHistoMappingGGInvMassPtBinPlot[iPt],
                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],2);
        }
    }
    canvasDataSpectra->Print(namePlot.Data());
    delete padDataSpectra;
    delete canvasDataSpectra;
}

//____________________________ Plotting Invariant Mass with Subtraction for Single Bin ______________________________________
void PlotWithFitSubtractedInvMassSinglePtBin(   TH1D * fHistoMappingSignalInvMassPtBinPlot, 
                                                TH1D** fHistoMappingTrueMesonInvMassPtBinsPlot,
                                                TF1 * fFitSignalInvMassPtBinPlot, 
                                                TString namePlot,
                                                TString nameCanvas,
                                                Double_t* fPlottingRangeMeson, 
                                                Bool_t fMonteCarloInfo,
                                                TString decayChannel){
    TCanvas *canvasDataFit      = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    canvasDataFit->SetTopMargin(0.12);
    canvasDataFit->SetBottomMargin(0.15);
    canvasDataFit->SetRightMargin(0.05);
    canvasDataFit->SetLeftMargin(0.15);

    fHistoMappingSignalInvMassPtBinPlot->SetAxisRange(fPlottingRangeMeson[0],fPlottingRangeMeson[1]);
    //     cout<<"Maximum::"<<fHistoMappingSignalInvMassPtBinPlot->GetMaximum()<<endl;
    DrawGammaHisto( fHistoMappingSignalInvMassPtBinPlot,
                    "2.4 GeV/#it{c} < #it{p}_{T} < 2.6 GeV/#it{c}",
                    Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                    fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
    if(fMonteCarloInfo){
        fHistoMappingTrueMesonInvMassPtBinsPlot[14]->SetLineColor(kRed);
        fHistoMappingTrueMesonInvMassPtBinsPlot[14]->SetLineWidth(1);
        fHistoMappingTrueMesonInvMassPtBinsPlot[14]->DrawCopy("same");
    }
    if(fMonteCarloInfo) fHistoMappingTrueMesonInvMassPtBinsPlot[14]->DrawCopy("same");
    if (fFitSignalInvMassPtBinPlot!=0x00){
        fFitSignalInvMassPtBinPlot->SetLineColor(kCyan+3);
        fFitSignalInvMassPtBinPlot->SetLineWidth(1);
        fFitSignalInvMassPtBinPlot->DrawCopy("same");
    } 
    canvasDataFit->Print(namePlot.Data());
    delete canvasDataFit;
}

//____________________________ Plotting Invariant Mass with Subtraction for Single Bin ______________________________________
void PlotWithFitSubtractedInvMassSinglePtBin2(  TH1D * fHistoMappingSignalInvMassPtBinPlot, 
                                                TH1D * fHistoMappingSignalInvMassPtBinPlot2,
                                                TH1D * fHistoMappingSignalInvMassPtBinPlot3, 
                                                TF1 * fFitSignalInvMassPtBinPlot, 
                                                TF1 * fFitBGInvMassPtBinPlot, 
                                                TString namePlot, 
                                                TString nameCanvas, 
                                                Double_t* fPlottingRangeMeson, 
                                                Bool_t fMonteCarloInfo, 
                                                TString decayChannel) {
    TCanvas *canvasDataFit      = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    canvasDataFit->SetTopMargin(0.12);
    canvasDataFit->SetBottomMargin(0.15);
    canvasDataFit->SetRightMargin(0.05);
    canvasDataFit->SetLeftMargin(0.15);
    fHistoMappingSignalInvMassPtBinPlot->SetAxisRange(fPlottingRangeMeson[0],fPlottingRangeMeson[1]);
//     cout<<"Maximum::"<<fHistoMappingSignalInvMassPtBinPlot->GetMaximum()<<endl;
    DrawGammaHisto( fHistoMappingSignalInvMassPtBinPlot,
                    "2.4 GeV/#it{c} < #it{p}_{T} < 2.6 GeV/#it{c}",
                    Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
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
        fFitBGInvMassPtBinPlot->SetLineWidth(1);
        fFitBGInvMassPtBinPlot->DrawCopy("same");
    } 
    if (fFitSignalInvMassPtBinPlot!=0x00){
        fFitSignalInvMassPtBinPlot->SetLineColor(kBlack);
        fFitSignalInvMassPtBinPlot->SetLineWidth(1);
        fFitSignalInvMassPtBinPlot->DrawCopy("same");
    } 

    canvasDataFit->Print(namePlot.Data());
    delete canvasDataFit;

        if(fMonteCarloInfo)fMonteCarloInfo=kFALSE; 
}

//____________________________ Plotting Invariant Mass Subtracted for all bins ________________________________________________________
void PlotWithFitSubtractedInvMassInPtBins( TH1D ** fHistoMappingSignalInvMassPtBinPlot, 
                                           TH1D** fHistoMappingTrueMesonInvMassPtBinsPlot, 
                                           TF1 ** fFitSignalInvMassPtBinPlot, 
                                           TString namePlot, 
                                           TString nameCanvas, 
                                           TString namePad, 
                                           Double_t* fPlottingRangeMeson, 
                                           TString dateDummy, 
                                           TString fMesonType, 
                                           Int_t fRowPlot, 
                                           Int_t fColumnPlot, 
                                           Int_t fStartBinPtRange, 
                                           Int_t fNumberPtBins, 
                                           Double_t* fRangeBinsPt, 
                                           TString fDecayChannel,
                                           Bool_t fMonteCarloInfo, 
                                           TString decayChannel, 
                                           TString fDetectionChannel, 
                                           TString fEnergy, 
                                           TString fTextMCvalidated ="",
                                           Bool_t labelData = kTRUE, 
                                           TString fTextFit = "Fit", 
                                           TString fTextMGammaGamma ="mixed evt. subtr. #it{M}_{#gamma#gamma}"){
    TCanvas *canvasDataFit          = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    canvasDataFit->SetTopMargin(0.00);
    canvasDataFit->SetBottomMargin(0.00);
    canvasDataFit->SetRightMargin(0.00);
    canvasDataFit->SetLeftMargin(0.00);

    TPad * padDataFit               = new TPad(namePad.Data(),"",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
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
    Int_t place                     = 0;
    for(Int_t iPt = fStartBinPtRange; iPt < fNumberPtBins; iPt++){
        Double_t startPt            = fRangeBinsPt[iPt];
        Double_t endPt              = fRangeBinsPt[iPt+1];

        place                       = place + 1; //give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataFit->cd(place);

            TString textAlice       = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo){
                textEvents          = "MC";
            } else {
                textEvents          = "Data";
            }
            Double_t nPixels        = 13;
            Double_t textHeight     = 0.08;
            if (padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) < padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1())){
                textHeight          = (Double_t)nPixels/padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) ;
            } else {
                textHeight          = (Double_t)nPixels/padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1());
            }
            Double_t startTextX     = 0.10;
            Double_t startTextY     = 0.9;
            Double_t differenceText = textHeight*1.25;

            TLatex *alice           = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
            TLatex *latexDate       = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy          = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
            TLatex *process         = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess      = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
            TLatex *events          = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

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
                Double_t totalheightLeg     = 2.;
                if (fTextMCvalidated.CompareTo("") != 0 && labelData){
                    totalheightLeg          = 3.;
                }
                TLegend* legendMC           = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+totalheightLeg)*differenceText);
                legendMC->SetTextSize(textHeight);
                legendMC->SetTextFont(62);
                legendMC->SetLineColor(0);
                legendMC->SetLineWidth(0);
                legendMC->SetFillStyle(0);
                legendMC->SetFillColor(0);
                legendMC->SetMargin(0.15);
                if (labelData){
                    Size_t markersize       = fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMarkerSize();
                    fHistoMappingSignalInvMassPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legendMC->AddEntry(fHistoMappingSignalInvMassPtBinPlot[iPt],fTextMGammaGamma.Data(),"ep");
                }
                if (fTextMCvalidated.CompareTo("") != 0){
                    Size_t markersize2      = fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMarkerSize();
                    if (labelData)fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetMarkerSize(3*markersize2);
                    legendMC->AddEntry(fHistoMappingTrueMesonInvMassPtBinsPlot[iPt],fTextMCvalidated.Data(),"ep");
                }
                if (fFitSignalInvMassPtBinPlot[iPt]!=0x00){
                    Size_t linewidth        = fFitSignalInvMassPtBinPlot[iPt]->GetLineWidth();
                    fFitSignalInvMassPtBinPlot[iPt]->SetLineWidth(5*linewidth);
                    legendMC->AddEntry(fFitSignalInvMassPtBinPlot[iPt],fTextFit.Data(),"l");
                }
                legendMC->Draw();
            }else {
                TLegend* legendData         = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+2.)*differenceText);
                legendData->SetTextSize(textHeight);
                legendData->SetTextFont(62);
                legendData->SetFillColor(0);
                legendData->SetFillStyle(0);
                legendData->SetLineWidth(0);
                legendData->SetLineColor(0);
                legendData->SetMargin(0.15);
                Size_t markersize           = fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMarkerSize();
                fHistoMappingSignalInvMassPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                legendData->AddEntry(fHistoMappingSignalInvMassPtBinPlot[iPt],fTextMGammaGamma.Data(),"ep");
                if (fFitSignalInvMassPtBinPlot[iPt]!=0x00){
                    Size_t linewidth        = fFitSignalInvMassPtBinPlot[iPt]->GetLineWidth();
                    fFitSignalInvMassPtBinPlot[iPt]->SetLineWidth(5*linewidth);
                    legendData->AddEntry(fFitSignalInvMassPtBinPlot[iPt],fTextFit.Data(),"l");
                }
                legendData->Draw();
            }

        } else {

            padDataFit->cd(place);
            padDataFit->cd(place)->SetTopMargin(0.12);
            padDataFit->cd(place)->SetBottomMargin(0.15);
            padDataFit->cd(place)->SetRightMargin(0.02);
            int remaining = (place-1)%fColumnPlot;
            if (remaining > 0) padDataFit->cd(place)->SetLeftMargin(0.15);
            else padDataFit->cd(place)->SetLeftMargin(0.25);

            cout << startPt << "-" << endPt <<endl;

            if (labelData) {
                DrawGammaHisto( fHistoMappingSignalInvMassPtBinPlot[iPt],
                                Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
                DrawGammaHisto( fHistoMappingSignalInvMassPtBinPlot[iPt],
                                Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                fPlottingRangeMeson[0],fPlottingRangeMeson[1],2);
                if (fHistoMappingSignalInvMassPtBinPlot[iPt]!=0x00){
                        TString nameOfPlot = fHistoMappingSignalInvMassPtBinPlot[iPt]->GetName();
                        Double_t mass = fMesonMass[iPt];
                        if (nameOfPlot.Contains("Left")) 
                            mass                        = fMesonMassLeft[iPt];
                        if (nameOfPlot.Contains("True")) 
                            mass                        = fMesonTrueMass[iPt];
                        Double_t intRangeLow            = mass + fMesonIntDeltaRange[0];
                        Double_t intRangeWideLow        = mass + fMesonIntDeltaRangeWide[0];
                        Double_t intRangeNarrowLow      = mass + fMesonIntDeltaRangeNarrow[0];
                        Double_t intRangeHigh           = mass + fMesonIntDeltaRange[1];
                        Double_t intRangeWideHigh       = mass + fMesonIntDeltaRangeWide[1];
                        Double_t intRangeNarrowHigh     = mass + fMesonIntDeltaRangeNarrow[1];

                        Double_t normalLow              = intRangeLow-(intRangeLow-fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->GetBinLowEdge(fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->FindBin(intRangeLow)));
                        Double_t normalUp               = intRangeHigh+(fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->GetBinUpEdge(fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->FindBin(intRangeHigh))-intRangeHigh);
                        Double_t wideLow                = intRangeWideLow-(intRangeWideLow-fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->GetBinLowEdge(fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->FindBin(intRangeWideLow)));
                        Double_t wideUp                 = intRangeWideHigh+(fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->GetBinUpEdge(fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->FindBin(intRangeWideHigh))-intRangeWideHigh);
                        Double_t narrowLow              = intRangeNarrowLow-(intRangeNarrowLow-fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->GetBinLowEdge(fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->FindBin(intRangeNarrowLow)));
                        Double_t narrowUp               = intRangeNarrowHigh+(fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->GetBinUpEdge(fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->FindBin(intRangeNarrowHigh))-intRangeNarrowHigh);

                        TLine * lmassPos                = new TLine (mass,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMinimum(),mass,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMaximum()*0.4);
                        lmassPos->SetLineColor(kRed+2);
                        lmassPos->SetLineStyle(1);
                        lmassPos->SetLineWidth(1);
                        lmassPos->Draw("same");
                        TLine * l1a                     = new TLine (normalLow,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMinimum(),normalLow,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMaximum());
                        l1a->SetLineColor(kGray+1);
                        l1a->SetLineStyle(1);
                        l1a->SetLineWidth(1);
                        l1a->Draw("same");
                        TLine * l1b                     = new TLine (normalUp,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMinimum(),normalUp,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMaximum());
                        l1b->SetLineColor(kGray+1);
                        l1b->SetLineStyle(1);
                        l1b->SetLineWidth(1);
                        l1b->Draw("same");
                        TLine * l2a                     = new TLine (wideLow,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMinimum(),wideLow,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMaximum());
                        l2a->SetLineColor(kGray+1);
                        l2a->SetLineStyle(2);
                        l2a->SetLineWidth(1);
                        l2a->Draw("same");
                        TLine * l2b                     = new TLine (wideUp,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMinimum(),wideUp,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMaximum());
                        l2b->SetLineColor(kGray+1);
                        l2b->SetLineStyle(2);
                        l2b->SetLineWidth(1);
                        l2b->Draw("same");
                        TLine * l3a                     = new TLine (narrowLow,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMinimum(),narrowLow,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMaximum());
                        l3a->SetLineColor(kGray+1);
                        l3a->SetLineStyle(3);
                        l3a->SetLineWidth(1);
                        l3a->Draw("same");
                        TLine * l3b                     = new TLine (narrowUp,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMinimum(),narrowUp,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMaximum());
                        l3b->SetLineColor(kGray+1);
                        l3b->SetLineStyle(3);
                        l3b->SetLineWidth(1);
                        l3b->Draw("same");
                } 
                DrawGammaHisto( fHistoMappingSignalInvMassPtBinPlot[iPt],
                                Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                fPlottingRangeMeson[0],fPlottingRangeMeson[1],2);
                if(fMonteCarloInfo){
                    fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetMarkerColor(kRed+2);
                    fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetLineColor(kRed+2);
                    fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetLineWidth(1);
                    fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->DrawCopy("same");
                }
            } else {
                DrawGammaHisto( fHistoMappingTrueMesonInvMassPtBinsPlot[iPt],
                                Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
                DrawGammaHisto( fHistoMappingTrueMesonInvMassPtBinsPlot[iPt],
                                Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                fPlottingRangeMeson[0],fPlottingRangeMeson[1],2);
                if (fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]!=0x00){
                        TString nameOfPlot = fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetName();
                        Double_t mass = fMesonMass[iPt];
                        if (nameOfPlot.Contains("Left")) 
                            mass                        = fMesonMassLeft[iPt];
                        if (nameOfPlot.Contains("True"))
                            mass                        = fMesonTrueMass[iPt];
                        
                        Double_t intRangeLow            = mass + fMesonIntDeltaRange[0];
                        Double_t intRangeWideLow        = mass + fMesonIntDeltaRangeWide[0];
                        Double_t intRangeNarrowLow      = mass + fMesonIntDeltaRangeNarrow[0];
                        Double_t intRangeHigh           = mass + fMesonIntDeltaRange[1];
                        Double_t intRangeWideHigh       = mass + fMesonIntDeltaRangeWide[1];
                        Double_t intRangeNarrowHigh     = mass + fMesonIntDeltaRangeNarrow[1];

                        Double_t normalLow              = intRangeLow-(intRangeLow-fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->GetBinLowEdge(fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->FindBin(intRangeLow)));
                        Double_t normalUp               = intRangeHigh+(fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->GetBinUpEdge(fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->FindBin(intRangeHigh))-intRangeHigh);
                        Double_t wideLow                = intRangeWideLow-(intRangeWideLow-fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->GetBinLowEdge(fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->FindBin(intRangeWideLow)));
                        Double_t wideUp                 = intRangeWideHigh+(fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->GetBinUpEdge(fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->FindBin(intRangeWideHigh))-intRangeWideHigh);
                        Double_t narrowLow              = intRangeNarrowLow-(intRangeNarrowLow-fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->GetBinLowEdge(fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->FindBin(intRangeNarrowLow)));
                        Double_t narrowUp               = intRangeNarrowHigh+(fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->GetBinUpEdge(fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->FindBin(intRangeNarrowHigh))-intRangeNarrowHigh);

                        TLine * lmassPos                = new TLine (mass,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMinimum(),mass,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMaximum()*0.4);
                        lmassPos->SetLineColor(kRed+2);
                        lmassPos->SetLineStyle(1);
                        lmassPos->SetLineWidth(1);
                        lmassPos->Draw("same");
                        
                        TLine * l1a                     = new TLine (normalLow,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMinimum(),normalLow,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMaximum());
                        l1a->SetLineColor(kGray+1);
                        l1a->SetLineStyle(1);
                        l1a->SetLineWidth(1);
                        l1a->Draw("same");
                        TLine * l1b                     = new TLine (normalUp,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMinimum(),normalUp,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMaximum());
                        l1b->SetLineColor(kGray+1);
                        l1b->SetLineStyle(1);
                        l1b->SetLineWidth(1);
                        l1b->Draw("same");
                        TLine * l2a                     = new TLine (wideLow,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMinimum(),wideLow,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMaximum());
                        l2a->SetLineColor(kGray+1);
                        l2a->SetLineStyle(2);
                        l2a->SetLineWidth(1);
                        l2a->Draw("same");
                        TLine * l2b                     = new TLine (wideUp,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMinimum(),wideUp,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMaximum());
                        l2b->SetLineColor(kGray+1);
                        l2b->SetLineStyle(2);
                        l2b->SetLineWidth(1);
                        l2b->Draw("same");
                        TLine * l3a                     = new TLine (narrowLow,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMinimum(),narrowLow,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMaximum());
                        l3a->SetLineColor(kGray+1);
                        l3a->SetLineStyle(3);
                        l3a->SetLineWidth(1);
                        l3a->Draw("same");
                        TLine * l3b                     = new TLine (narrowUp,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMinimum(),narrowUp,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMaximum());
                        l3b->SetLineColor(kGray+1);
                        l3b->SetLineStyle(3);
                        l3b->SetLineWidth(1);
                        l3b->Draw("same");
                } 
                
                
                fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetMarkerColor(kRed+2);
                fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetMarkerSize(0.5);
                fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetLineColor(kRed+2);
//                 fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetLineWidth(1);
                fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->DrawCopy("same,p,e1");
                
            }
            if (fFitSignalInvMassPtBinPlot[iPt]!=0x00){
                fFitSignalInvMassPtBinPlot[iPt]->SetLineColor(kCyan+3);
                fFitSignalInvMassPtBinPlot[iPt]->SetLineWidth(1);
                fFitSignalInvMassPtBinPlot[iPt]->DrawCopy("same");
            } 
        }
    }
    canvasDataFit->Print(namePlot.Data());
    delete padDataFit;
    delete canvasDataFit;
}

//____________________________ Plotting Invariant Mass Subtracted for all bins ________________________________________________________
void PlotWithManyFitSubtractedInvMassInPtBins(  TH1D ** fHistoMappingSignalInvMassPtBinPlot, 
                                                TF1 ** fFitSignalInvMassPtBinPlot, 
                                                TF1 *** fFitSigWithOtherBGInvMassPtBinPlot, 
                                                Int_t nFits,
                                                TString* fTextFitAdd,
                                                TString namePlot, 
                                                TString nameCanvas, 
                                                TString namePad, 
                                                Double_t* fPlottingRangeMeson, 
                                                TString dateDummy, 
                                                TString fMesonType, 
                                                Int_t fRowPlot, 
                                                Int_t fColumnPlot, 
                                                Int_t fStartBinPtRange, 
                                                Int_t fNumberPtBins, 
                                                Double_t* fRangeBinsPt, 
                                                TString fDecayChannel,
                                                Bool_t fMonteCarloInfo, 
                                                TString decayChannel, 
                                                TString fDetectionChannel, 
                                                TString fEnergy, 
                                                TString fTextMCvalidated    = "",
                                                Bool_t labelData            = kTRUE, 
                                                TString fTextFit            = "Fit", 
                                                TString fTextMGammaGamma    = "mixed evt. subtr. #it{M}_{#gamma#gamma}" ){
    
    TCanvas *canvasDataFit          = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    canvasDataFit->SetTopMargin(0.00);
    canvasDataFit->SetBottomMargin(0.00);
    canvasDataFit->SetRightMargin(0.00);
    canvasDataFit->SetLeftMargin(0.00);

    TPad * padDataFit               = new TPad(namePad.Data(),"",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
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
    Int_t place                     = 0;
    for(Int_t iPt = fStartBinPtRange; iPt < fNumberPtBins; iPt++){
        Double_t startPt            = fRangeBinsPt[iPt];
        Double_t endPt              = fRangeBinsPt[iPt+1];

        place                       = place + 1; //give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataFit->cd(place);

            TString textAlice       = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo){
                textEvents          = "MC";
            } else {
                textEvents          = "Data";
            }
            Double_t nPixels        = 13;
            Double_t textHeight     = 0.08;
            if (padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) < padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1())){
                textHeight          = (Double_t)nPixels/padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) ;
            } else {
                textHeight          = (Double_t)nPixels/padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1());
            }
            Double_t startTextX     = 0.10;
            Double_t startTextY     = 0.9;
            Double_t differenceText = textHeight*1.25;

            TLatex *alice           = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
            TLatex *latexDate       = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy          = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
            TLatex *process         = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess      = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
            TLatex *events          = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

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
                Double_t totalheightLeg     = 2.;
                totalheightLeg = totalheightLeg+nFits;
                
                TLegend* legendMC           = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+totalheightLeg)*differenceText);
                legendMC->SetTextSize(textHeight);
                legendMC->SetTextFont(62);
                legendMC->SetLineColor(0);
                legendMC->SetLineWidth(0);
                legendMC->SetFillStyle(0);
                legendMC->SetFillColor(0);
                legendMC->SetMargin(0.15);
                if (labelData){
                    Size_t markersize       = fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMarkerSize();
                    fHistoMappingSignalInvMassPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legendMC->AddEntry(fHistoMappingSignalInvMassPtBinPlot[iPt],fTextMGammaGamma.Data(),"ep");
                }
                if (fFitSignalInvMassPtBinPlot[iPt]!=0x00){
                    Size_t linewidth        = fFitSignalInvMassPtBinPlot[iPt]->GetLineWidth();
                    fFitSignalInvMassPtBinPlot[iPt]->SetLineWidth(3*linewidth);
                    legendMC->AddEntry(fFitSignalInvMassPtBinPlot[iPt],fTextFit.Data(),"l");
                }
                for (Int_t m = 0; (m < nFits && m < 3); m++){
                    if (fFitSigWithOtherBGInvMassPtBinPlot[m][iPt]!=0x00){
                        Size_t linewidth        = fFitSigWithOtherBGInvMassPtBinPlot[m][iPt]->GetLineWidth();
                        fFitSigWithOtherBGInvMassPtBinPlot[m][iPt]->SetLineWidth(3*linewidth);
                        legendMC->AddEntry(fFitSigWithOtherBGInvMassPtBinPlot[m][iPt],fTextFitAdd[m].Data(),"l");
                    }
                }    
                
                legendMC->Draw();
            }else {
                TLegend* legendData         = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+2.+nFits)*differenceText);
                legendData->SetTextSize(textHeight);
                legendData->SetTextFont(62);
                legendData->SetFillColor(0);
                legendData->SetFillStyle(0);
                legendData->SetLineWidth(0);
                legendData->SetLineColor(0);
                legendData->SetMargin(0.15);
                Size_t markersize           = fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMarkerSize();
                fHistoMappingSignalInvMassPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                legendData->AddEntry(fHistoMappingSignalInvMassPtBinPlot[iPt],fTextMGammaGamma.Data(),"ep");
                if (fFitSignalInvMassPtBinPlot[iPt]!=0x00){
                    Size_t linewidth        = fFitSignalInvMassPtBinPlot[iPt]->GetLineWidth();
                    fFitSignalInvMassPtBinPlot[iPt]->SetLineWidth(3*linewidth);
                    legendData->AddEntry(fFitSignalInvMassPtBinPlot[iPt],fTextFit.Data(),"l");
                }
                for (Int_t m = 0; (m < nFits && m < 3); m++){
                    if (fFitSigWithOtherBGInvMassPtBinPlot[m][iPt]!=0x00){
                        Size_t linewidth        = fFitSigWithOtherBGInvMassPtBinPlot[m][iPt]->GetLineWidth();
                        fFitSigWithOtherBGInvMassPtBinPlot[m][iPt]->SetLineWidth(3*linewidth);
                        legendData->AddEntry(fFitSigWithOtherBGInvMassPtBinPlot[m][iPt],fTextFitAdd[m].Data(),"l");
                    }
                }    
                legendData->Draw();
            }

        } else {

            padDataFit->cd(place);
            padDataFit->cd(place)->SetTopMargin(0.12);
            padDataFit->cd(place)->SetBottomMargin(0.15);
            padDataFit->cd(place)->SetRightMargin(0.02);
            int remaining = (place-1)%fColumnPlot;
            if (remaining > 0) padDataFit->cd(place)->SetLeftMargin(0.15);
            else padDataFit->cd(place)->SetLeftMargin(0.25);

            cout << startPt << "-" << endPt <<endl;
            DrawGammaHisto( fHistoMappingSignalInvMassPtBinPlot[iPt],
                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);

            if (fFitSignalInvMassPtBinPlot[iPt]!=0x00){
                fFitSignalInvMassPtBinPlot[iPt]->SetLineColor(kCyan+3);
                fFitSignalInvMassPtBinPlot[iPt]->SetLineWidth(1.5);
                fFitSignalInvMassPtBinPlot[iPt]->DrawCopy("same");
            } 
            Color_t colorFit[3]     = {kRed+1, kAzure+2, 807};
            Style_t styleFit[3]     = {7, 3, 6};
            
            for (Int_t m = 0; (m < nFits && m < 3); m++){
                if (fFitSigWithOtherBGInvMassPtBinPlot[m][iPt]){
                    cout << m << "\t"<< iPt << "\t" << fFitSigWithOtherBGInvMassPtBinPlot[m][iPt] << endl;
                    fFitSigWithOtherBGInvMassPtBinPlot[m][iPt]->SetLineColor(colorFit[m]);
                    fFitSigWithOtherBGInvMassPtBinPlot[m][iPt]->SetLineStyle(styleFit[m]);
                    fFitSigWithOtherBGInvMassPtBinPlot[m][iPt]->SetLineWidth(1.5);
                    fFitSigWithOtherBGInvMassPtBinPlot[m][iPt]->DrawCopy("same");    
                }   
            }    
        }
    }
    canvasDataFit->Print(namePlot.Data());
    delete padDataFit;
    delete canvasDataFit;
}


//____________________________ Plotting Invariant Mass Subtracted for all bins ________________________________________________________
void PlotWith2FitsSubtractedInvMassInPtBins( TH1D ** fHistoMappingSignalInvMassPtBinPlot, 
                                           TH1D** fHistoMappingTrueMesonInvMassPtBinsPlot, 
                                           TF1 ** fFitSignalInvMassPtBinPlot, 
                                           TF1 ** fFitLinearBck,
                                           TString namePlot, 
                                           TString nameCanvas, 
                                           TString namePad, 
                                           Double_t* fPlottingRangeMeson, 
                                           TString dateDummy, 
                                           TString fMesonType, 
                                           Int_t fRowPlot, 
                                           Int_t fColumnPlot, 
                                           Int_t fStartBinPtRange, 
                                           Int_t fNumberPtBins, 
                                           Double_t* fRangeBinsPt, 
                                           TString fDecayChannel,
                                           Bool_t fMonteCarloInfo, 
                                           TString decayChannel, 
                                           TString fDetectionChannel, 
                                           TString fEnergy, 
                                           TString fTextMCvalidated ="",
                                           Bool_t labelData = kTRUE, 
                                           TString fTextFit = "Fit", 
                                           TString fTextMGammaGamma ="mixed evt. subtr. #it{M}_{#gamma#gamma}"){
    TCanvas *canvasDataFit          = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    canvasDataFit->SetTopMargin(0.00);
    canvasDataFit->SetBottomMargin(0.00);
    canvasDataFit->SetRightMargin(0.00);
    canvasDataFit->SetLeftMargin(0.00);

    TPad * padDataFit               = new TPad(namePad.Data(),"",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
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
    Int_t place                     = 0;
    for(Int_t iPt = fStartBinPtRange; iPt < fNumberPtBins; iPt++){
        Double_t startPt            = fRangeBinsPt[iPt];
        Double_t endPt              = fRangeBinsPt[iPt+1];

        place                       = place + 1; //give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataFit->cd(place);

            TString textAlice       = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo){
                textEvents          = "MC";
            } else {
                textEvents          = "Data";
            }
            Double_t nPixels        = 13;
            Double_t textHeight     = 0.08;
            if (padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) < padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1())){
                textHeight          = (Double_t)nPixels/padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) ;
            } else {
                textHeight          = (Double_t)nPixels/padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1());
            }
            Double_t startTextX     = 0.10;
            Double_t startTextY     = 0.9;
            Double_t differenceText = textHeight*1.25;

            TLatex *alice           = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
            TLatex *latexDate       = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy          = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
            TLatex *process         = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess      = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
            TLatex *events          = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

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
                Double_t totalheightLeg     = 3.;
                if (fTextMCvalidated.CompareTo("") != 0 && labelData){
                    totalheightLeg          = 4.;
                }
                TLegend* legendMC           = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+totalheightLeg)*differenceText);
                legendMC->SetTextSize(textHeight);
                legendMC->SetTextFont(62);
                legendMC->SetLineColor(0);
                legendMC->SetLineWidth(0);
                legendMC->SetFillStyle(0);
                legendMC->SetFillColor(0);
                legendMC->SetMargin(0.15);
                if (labelData){
                    Size_t markersize       = fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMarkerSize();
                    fHistoMappingSignalInvMassPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legendMC->AddEntry(fHistoMappingSignalInvMassPtBinPlot[iPt],fTextMGammaGamma.Data(),"ep");
                }
                if (fTextMCvalidated.CompareTo("") != 0){
                    Size_t markersize2      = fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMarkerSize();
                    if (labelData)fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetMarkerSize(3*markersize2);
                    legendMC->AddEntry(fHistoMappingTrueMesonInvMassPtBinsPlot[iPt],fTextMCvalidated.Data(),"ep");
                }
                if (fFitSignalInvMassPtBinPlot[iPt]!=0x00){
                    Size_t linewidth        = fFitSignalInvMassPtBinPlot[iPt]->GetLineWidth();
                    fFitSignalInvMassPtBinPlot[iPt]->SetLineWidth(5*linewidth);
                    legendMC->AddEntry(fFitSignalInvMassPtBinPlot[iPt],fTextFit.Data(),"l");

                    fFitLinearBck[iPt]->SetLineWidth(5*linewidth);
                    legendMC->AddEntry(fFitLinearBck[iPt],"BG fit","l");
                    
                }
                legendMC->Draw();
            }else {
                TLegend* legendData         = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+2.)*differenceText);
                legendData->SetTextSize(textHeight);
                legendData->SetTextFont(62);
                legendData->SetFillColor(0);
                legendData->SetFillStyle(0);
                legendData->SetLineWidth(0);
                legendData->SetLineColor(0);
                legendData->SetMargin(0.15);
                Size_t markersize           = fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMarkerSize();
                fHistoMappingSignalInvMassPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                legendData->AddEntry(fHistoMappingSignalInvMassPtBinPlot[iPt],fTextMGammaGamma.Data(),"ep");
                if (fFitSignalInvMassPtBinPlot[iPt]!=0x00){
                    Size_t linewidth        = fFitSignalInvMassPtBinPlot[iPt]->GetLineWidth();
                    fFitSignalInvMassPtBinPlot[iPt]->SetLineWidth(5*linewidth);
                    legendData->AddEntry(fFitSignalInvMassPtBinPlot[iPt],fTextFit.Data(),"l");
                    
                    fFitLinearBck[iPt]->SetLineWidth(5*linewidth);
                    legendData->AddEntry(fFitLinearBck[iPt],"BG fit","l");
                  
                }
                legendData->Draw();
            }

        } else {

            padDataFit->cd(place);
            padDataFit->cd(place)->SetTopMargin(0.12);
            padDataFit->cd(place)->SetBottomMargin(0.15);
            padDataFit->cd(place)->SetRightMargin(0.02);
            int remaining = (place-1)%fColumnPlot;
            if (remaining > 0) padDataFit->cd(place)->SetLeftMargin(0.15);
            else padDataFit->cd(place)->SetLeftMargin(0.25);

            cout << startPt << "-" << endPt <<endl;

            if (labelData) {
                DrawGammaHisto( fHistoMappingSignalInvMassPtBinPlot[iPt],
                                Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
                DrawGammaHisto( fHistoMappingSignalInvMassPtBinPlot[iPt],
                                Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                fPlottingRangeMeson[0],fPlottingRangeMeson[1],2);
                if (fHistoMappingSignalInvMassPtBinPlot[iPt]!=0x00){
                        TString nameOfPlot = fHistoMappingSignalInvMassPtBinPlot[iPt]->GetName();
                        Double_t mass = fMesonMass[iPt];
                        if (nameOfPlot.Contains("Left")) 
                            mass                        = fMesonMassLeft[iPt];
                        if (nameOfPlot.Contains("True")) 
                            mass                        = fMesonTrueMass[iPt];
                        Double_t intRangeLow            = mass + fMesonIntDeltaRange[0];
                        Double_t intRangeWideLow        = mass + fMesonIntDeltaRangeWide[0];
                        Double_t intRangeNarrowLow      = mass + fMesonIntDeltaRangeNarrow[0];
                        Double_t intRangeHigh           = mass + fMesonIntDeltaRange[1];
                        Double_t intRangeWideHigh       = mass + fMesonIntDeltaRangeWide[1];
                        Double_t intRangeNarrowHigh     = mass + fMesonIntDeltaRangeNarrow[1];

                        Double_t normalLow              = intRangeLow-(intRangeLow-fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->GetBinLowEdge(fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->FindBin(intRangeLow)));
                        Double_t normalUp               = intRangeHigh+(fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->GetBinUpEdge(fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->FindBin(intRangeHigh))-intRangeHigh);
                        Double_t wideLow                = intRangeWideLow-(intRangeWideLow-fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->GetBinLowEdge(fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->FindBin(intRangeWideLow)));
                        Double_t wideUp                 = intRangeWideHigh+(fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->GetBinUpEdge(fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->FindBin(intRangeWideHigh))-intRangeWideHigh);
                        Double_t narrowLow              = intRangeNarrowLow-(intRangeNarrowLow-fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->GetBinLowEdge(fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->FindBin(intRangeNarrowLow)));
                        Double_t narrowUp               = intRangeNarrowHigh+(fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->GetBinUpEdge(fHistoMappingSignalInvMassPtBinPlot[iPt]->GetXaxis()->FindBin(intRangeNarrowHigh))-intRangeNarrowHigh);

                        TLine * lmassPos                = new TLine (mass,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMinimum(),mass,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMaximum()*0.4);
                        lmassPos->SetLineColor(kRed+2);
                        lmassPos->SetLineStyle(1);
                        lmassPos->SetLineWidth(1);
                        lmassPos->Draw("same");
                        TLine * l1a                     = new TLine (normalLow,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMinimum(),normalLow,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMaximum());
                        l1a->SetLineColor(kGray+1);
                        l1a->SetLineStyle(1);
                        l1a->SetLineWidth(1);
                        l1a->Draw("same");
                        TLine * l1b                     = new TLine (normalUp,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMinimum(),normalUp,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMaximum());
                        l1b->SetLineColor(kGray+1);
                        l1b->SetLineStyle(1);
                        l1b->SetLineWidth(1);
                        l1b->Draw("same");
                        TLine * l2a                     = new TLine (wideLow,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMinimum(),wideLow,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMaximum());
                        l2a->SetLineColor(kGray+1);
                        l2a->SetLineStyle(2);
                        l2a->SetLineWidth(1);
                        l2a->Draw("same");
                        TLine * l2b                     = new TLine (wideUp,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMinimum(),wideUp,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMaximum());
                        l2b->SetLineColor(kGray+1);
                        l2b->SetLineStyle(2);
                        l2b->SetLineWidth(1);
                        l2b->Draw("same");
                        TLine * l3a                     = new TLine (narrowLow,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMinimum(),narrowLow,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMaximum());
                        l3a->SetLineColor(kGray+1);
                        l3a->SetLineStyle(3);
                        l3a->SetLineWidth(1);
                        l3a->Draw("same");
                        TLine * l3b                     = new TLine (narrowUp,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMinimum(),narrowUp,fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMaximum());
                        l3b->SetLineColor(kGray+1);
                        l3b->SetLineStyle(3);
                        l3b->SetLineWidth(1);
                        l3b->Draw("same");
                } 
                DrawGammaHisto( fHistoMappingSignalInvMassPtBinPlot[iPt],
                                Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                fPlottingRangeMeson[0],fPlottingRangeMeson[1],2);
                if(fMonteCarloInfo){
                    fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetMarkerColor(kRed+2);
                    fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetLineColor(kRed+2);
                    fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetLineWidth(1);
                    fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->DrawCopy("same");
                }
            } else {
                DrawGammaHisto( fHistoMappingTrueMesonInvMassPtBinsPlot[iPt],
                                Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
                DrawGammaHisto( fHistoMappingTrueMesonInvMassPtBinsPlot[iPt],
                                Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                fPlottingRangeMeson[0],fPlottingRangeMeson[1],2);
                if (fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]!=0x00){
                        TString nameOfPlot = fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetName();
                        Double_t mass = fMesonMass[iPt];
                        if (nameOfPlot.Contains("Left")) 
                            mass                        = fMesonMassLeft[iPt];
                        if (nameOfPlot.Contains("True"))
                            mass                        = fMesonTrueMass[iPt];
                        
                        Double_t intRangeLow            = mass + fMesonIntDeltaRange[0];
                        Double_t intRangeWideLow        = mass + fMesonIntDeltaRangeWide[0];
                        Double_t intRangeNarrowLow      = mass + fMesonIntDeltaRangeNarrow[0];
                        Double_t intRangeHigh           = mass + fMesonIntDeltaRange[1];
                        Double_t intRangeWideHigh       = mass + fMesonIntDeltaRangeWide[1];
                        Double_t intRangeNarrowHigh     = mass + fMesonIntDeltaRangeNarrow[1];

                        Double_t normalLow              = intRangeLow-(intRangeLow-fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->GetBinLowEdge(fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->FindBin(intRangeLow)));
                        Double_t normalUp               = intRangeHigh+(fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->GetBinUpEdge(fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->FindBin(intRangeHigh))-intRangeHigh);
                        Double_t wideLow                = intRangeWideLow-(intRangeWideLow-fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->GetBinLowEdge(fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->FindBin(intRangeWideLow)));
                        Double_t wideUp                 = intRangeWideHigh+(fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->GetBinUpEdge(fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->FindBin(intRangeWideHigh))-intRangeWideHigh);
                        Double_t narrowLow              = intRangeNarrowLow-(intRangeNarrowLow-fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->GetBinLowEdge(fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->FindBin(intRangeNarrowLow)));
                        Double_t narrowUp               = intRangeNarrowHigh+(fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->GetBinUpEdge(fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetXaxis()->FindBin(intRangeNarrowHigh))-intRangeNarrowHigh);

                        TLine * lmassPos                = new TLine (mass,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMinimum(),mass,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMaximum()*0.4);
                        lmassPos->SetLineColor(kRed+2);
                        lmassPos->SetLineStyle(1);
                        lmassPos->SetLineWidth(1);
                        lmassPos->Draw("same");
                        
                        TLine * l1a                     = new TLine (normalLow,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMinimum(),normalLow,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMaximum());
                        l1a->SetLineColor(kGray+1);
                        l1a->SetLineStyle(1);
                        l1a->SetLineWidth(1);
                        l1a->Draw("same");
                        TLine * l1b                     = new TLine (normalUp,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMinimum(),normalUp,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMaximum());
                        l1b->SetLineColor(kGray+1);
                        l1b->SetLineStyle(1);
                        l1b->SetLineWidth(1);
                        l1b->Draw("same");
                        TLine * l2a                     = new TLine (wideLow,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMinimum(),wideLow,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMaximum());
                        l2a->SetLineColor(kGray+1);
                        l2a->SetLineStyle(2);
                        l2a->SetLineWidth(1);
                        l2a->Draw("same");
                        TLine * l2b                     = new TLine (wideUp,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMinimum(),wideUp,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMaximum());
                        l2b->SetLineColor(kGray+1);
                        l2b->SetLineStyle(2);
                        l2b->SetLineWidth(1);
                        l2b->Draw("same");
                        TLine * l3a                     = new TLine (narrowLow,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMinimum(),narrowLow,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMaximum());
                        l3a->SetLineColor(kGray+1);
                        l3a->SetLineStyle(3);
                        l3a->SetLineWidth(1);
                        l3a->Draw("same");
                        TLine * l3b                     = new TLine (narrowUp,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMinimum(),narrowUp,fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->GetMaximum());
                        l3b->SetLineColor(kGray+1);
                        l3b->SetLineStyle(3);
                        l3b->SetLineWidth(1);
                        l3b->Draw("same");
                } 
                
                
                fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetMarkerColor(kRed+2);
                fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetMarkerSize(0.5);
                fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetLineColor(kRed+2);
//                 fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetLineWidth(1);
                fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->DrawCopy("same,p,e1");
                
            }
            if (fFitSignalInvMassPtBinPlot[iPt]!=0x00){
                fFitSignalInvMassPtBinPlot[iPt]->SetLineColor(kCyan+3);
                fFitSignalInvMassPtBinPlot[iPt]->SetLineWidth(1);
                fFitSignalInvMassPtBinPlot[iPt]->DrawCopy("same");
                
                fFitLinearBck[iPt]->SetLineColor(kBlue);   
                fFitLinearBck[iPt]->SetLineWidth(1);
                fFitLinearBck[iPt]->DrawCopy("same");
            } 
        }
    }
    canvasDataFit->Print(namePlot.Data());
    delete padDataFit;
    delete canvasDataFit;
}



//____________________________ Plotting Invariant Mass Subtracted for all bins ________________________________________________________
void PlotWithBGFitSubtractedInvMassInPtBins(    TH1D ** fHistoMappingSignalPlusBGInvMassPtBinPlot, 
                                                TH1D ** fHistoMappingBG, 
                                                TH1D** fHistoMappingSignal, 
                                                TF1 ** fFitBGInvMassPtBinPlot,
                                                TString namePlot,
                                                TString nameCanvas,
                                                TString namePad, 
                                                Double_t* fPlottingRangeMeson, 
                                                TString dateDummy, 
                                                TString fMesonType, 
                                                Int_t fRowPlot, 
                                                Int_t fColumnPlot,
                                                Int_t fStartBinPtRange, 
                                                Int_t fNumberPtBins, 
                                                Double_t* fRangeBinsPt, 
                                                TString fDecayChannel, 
                                                Bool_t fMonteCarloInfo, 
                                                TString decayChannel, 
                                                TString fDetectionChannel, 
                                                TString fEnergy){
    TCanvas *canvasDataFit          = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    canvasDataFit->SetTopMargin(0.00);
    canvasDataFit->SetBottomMargin(0.00);
    canvasDataFit->SetRightMargin(0.00);
    canvasDataFit->SetLeftMargin(0.00);

    TPad * padDataFit               = new TPad(namePad.Data(),"",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
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
    Int_t place                     = 0;
    for(Int_t iPt = fStartBinPtRange; iPt < fNumberPtBins; iPt++){
        Double_t startPt            = fRangeBinsPt[iPt];
        Double_t endPt              = fRangeBinsPt[iPt+1];

        place                       = place + 1; //give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataFit->cd(place);

            TString textAlice           = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo){
                textEvents              = "MC";
            } else {
                textEvents              = "Data";
            }
            Double_t nPixels            = 13;
            Double_t textHeight         = 0.08;
            if (padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) < padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1())){
                textHeight              = (Double_t)nPixels/padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) ;
            } else {
                textHeight              = (Double_t)nPixels/padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1());
            }
            Double_t startTextX         = 0.10;
            Double_t startTextY         = 0.9;
            Double_t differenceText     = textHeight*1.25;

            TLatex *alice               = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
            TLatex *latexDate           = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy              = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
            TLatex *process             = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess          = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
            TLatex *events              = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

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

            TLegend* legendData         = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+4.)*differenceText);
            legendData->SetTextSize(textHeight);
            legendData->SetTextFont(62);
            legendData->SetFillColor(0);
            legendData->SetFillStyle(0);
            legendData->SetLineWidth(0);
            legendData->SetLineColor(0);
            legendData->SetMargin(0.15);
            Size_t markersize           = fHistoMappingSignalPlusBGInvMassPtBinPlot[iPt]->GetMarkerSize();
            fHistoMappingSignalPlusBGInvMassPtBinPlot[iPt]->SetMarkerSize(3*markersize);
            legendData->AddEntry(fHistoMappingSignalPlusBGInvMassPtBinPlot[iPt],"mixed evt. subtr. #it{M}_{#gamma#gamma}","ep");
            fHistoMappingBG[iPt]->SetMarkerSize(3*markersize);
            legendData->AddEntry(fHistoMappingBG[iPt],"fitted add. BG","ep");
            fHistoMappingSignal[iPt]->SetMarkerSize(5*markersize);
            legendData->AddEntry(fHistoMappingSignal[iPt],"signal","ep");
            if (fFitBGInvMassPtBinPlot[iPt]!=0x00){
                Size_t linewidth        = fFitBGInvMassPtBinPlot[iPt]->GetLineWidth();
                fFitBGInvMassPtBinPlot[iPt]->SetLineWidth(5*linewidth);
                legendData->AddEntry(fFitBGInvMassPtBinPlot[iPt],"BG fit","l");
            }
            legendData->Draw();

        } else {

            padDataFit->cd(place);
            padDataFit->cd(place)->SetTopMargin(0.12);
            padDataFit->cd(place)->SetBottomMargin(0.15);
            padDataFit->cd(place)->SetRightMargin(0.02);
            int remaining               = (place-1)%fColumnPlot;
            if (remaining > 0) padDataFit->cd(place)->SetLeftMargin(0.15);
            else padDataFit->cd(place)->SetLeftMargin(0.25);

            cout << startPt << "-" << endPt <<endl;

            DrawGammaHisto( fHistoMappingSignalPlusBGInvMassPtBinPlot[iPt],
                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
            DrawGammaHisto( fHistoMappingSignalPlusBGInvMassPtBinPlot[iPt],
                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],2);

            fHistoMappingBG[iPt]->SetMarkerSize(0.2);
            fHistoMappingBG[iPt]->SetMarkerColor(kCyan+1);
            fHistoMappingBG[iPt]->SetLineColor(kCyan+1);
            fHistoMappingBG[iPt]->SetLineWidth(1);
            fHistoMappingBG[iPt]->DrawCopy("same,pe1");
            
            fHistoMappingSignal[iPt]->SetMarkerSize(0.2);
            fHistoMappingSignal[iPt]->SetMarkerColor(kRed+2);
            fHistoMappingSignal[iPt]->SetLineColor(kRed+2);
            fHistoMappingSignal[iPt]->SetLineWidth(1);
            fHistoMappingSignal[iPt]->DrawCopy("same,pe1");

            if (fFitBGInvMassPtBinPlot[iPt]!=0x00){
                fFitBGInvMassPtBinPlot[iPt]->SetLineColor(kCyan+3);
                fFitBGInvMassPtBinPlot[iPt]->SetLineWidth(1);
                fFitBGInvMassPtBinPlot[iPt]->DrawCopy("same");
                TLine * l1 = new TLine (fFitBGInvMassPtBinPlot[iPt]->GetParameter(1),fHistoMappingSignalPlusBGInvMassPtBinPlot[iPt]->GetMinimum(),fFitBGInvMassPtBinPlot[iPt]->GetParameter(1),fHistoMappingSignalPlusBGInvMassPtBinPlot[iPt]->GetMaximum());
                l1->SetLineColor(4);
                l1->SetLineWidth(1);
            } 

        }
    }
    canvasDataFit->Print(namePlot.Data());
    delete padDataFit;
    delete canvasDataFit;
}



//______________________________ Invariant Mass with Fit to Peak Position for all Pt bins __________________________________
void PlotWithFitPeakPosInvMassInPtBins( TH1D ** fHistoMappingSignalInvMassPtBinPlot, 
                                        TH1D **fHistoMappingTrueMesonInvMassPtBinsPlot,
                                        TF1 ** fFitSignalInvMassPtBinPlot, 
                                        TString namePlot, 
                                        TString nameCanvas, 
                                        TString namePad,
                                        Double_t* fPlottingRangeMeson, 
                                        TString dateDummy, 
                                        TString fMesonType, 
                                        Int_t fRowPlot, 
                                        Int_t fColumnPlot, 
                                        Int_t fStartBinPtRange, 
                                        Int_t fNumberPtBins,
                                        Double_t* fRangeBinsPt, 
                                        TString fDecayChannel, 
                                        Bool_t fMonteCarloInfo, 
                                        TString decayChannel,
                                        TString fDetectionChannel,
                                        TString fEnergy){
    TCanvas *canvasDataFit          = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    canvasDataFit->SetTopMargin(0.02);
    canvasDataFit->SetBottomMargin(0.02);
    canvasDataFit->SetRightMargin(0.02);
    canvasDataFit->SetLeftMargin(0.02);

    TPad * padDataFit               = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    padDataFit->SetFillColor(0);
    padDataFit->GetFrame()->SetFillColor(0);
    padDataFit->SetBorderMode(0);
    padDataFit->Divide(fColumnPlot,fRowPlot,0.0,0.0);
    padDataFit->SetLeftMargin(0.2);
    padDataFit->Draw();

    cout << fMesonType.Data() << endl;
    Int_t place                     = 0;
    for(Int_t iPt = fStartBinPtRange; iPt < fNumberPtBins; iPt++){
        Double_t startPt            = fRangeBinsPt[iPt];
        Double_t endPt              = fRangeBinsPt[iPt+1];

        place                       = place + 1; //give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataFit->cd(place);

            TString textAlice       = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo){
                textEvents          = "MC";
            } else { 
                textEvents          = "Data";
            }
            Double_t nPixels        = 13;
            Double_t textHeight     = 0.08;
            if (padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) < padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1())){
                textHeight          = (Double_t)nPixels/padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) ;
            } else {
                textHeight          = (Double_t)nPixels/padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1());
            }
            Double_t startTextX     = 0.10;
            Double_t startTextY     = 0.8;
            Double_t differenceText = textHeight*1.25;

            TLatex *alice           = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
            TLatex *latexDate       = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy          = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
            TLatex *process         = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess      = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
            TLatex *events          = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

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
            int remaining           = (place-1)%fColumnPlot;
            if (remaining > 0) padDataFit->cd(place)->SetLeftMargin(0.15);
            else padDataFit->cd(place)->SetLeftMargin(0.25);

            cout << startPt << "-" << endPt <<endl;

            DrawGammaHisto( fHistoMappingSignalInvMassPtBinPlot[iPt],
                            Form("%3.2f GeV/#it{c} < E < %3.2f GeV/#it{c}",startPt,endPt),
                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
            DrawGammaHisto( fHistoMappingSignalInvMassPtBinPlot[iPt],
                            Form("%3.2f GeV/#it{c} < E < %3.2f GeV/#it{c}",startPt,endPt),
                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],2);

            if(fMonteCarloInfo){
                fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetLineColor(kRed);
                fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetLineWidth(1);
                fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->DrawCopy("same");
            }
            if (fFitSignalInvMassPtBinPlot[iPt]!=0x00){
                fFitSignalInvMassPtBinPlot[iPt]->SetLineColor(kCyan+3);
                fFitSignalInvMassPtBinPlot[iPt]->SetLineWidth(1);
                fFitSignalInvMassPtBinPlot[iPt]->DrawCopy("same");
                TLine * l1          = new TLine (fFitSignalInvMassPtBinPlot[iPt]->GetParameter(1),fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMinimum(),fFitSignalInvMassPtBinPlot[iPt]->GetParameter(1),fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMaximum());
                l1->SetLineColor(4);
                l1->SetLineWidth(1);
            } 

        }
    }
    canvasDataFit->Print(namePlot.Data());
    delete padDataFit;
    delete canvasDataFit;
}



//____________________________ Plotting Invariant Mass validated splitted in different categories for all bins ________________________________________________________
void PlotTrueInvMassSplittedInPhotonAndElectronInPtBins(    TH1D ** fHistoTrueSignal, 
                                                            TH1D** fHistoTrueSignalPhotons, 
                                                            TH1D** fHistoTrueSignalElectrons, 
                                                            TH1D** fHistoTrueSignalConvPhotons, 
                                                            TH1D** fHistoTrueSignalMixed, 
                                                            TString namePlot, 
                                                            TString nameCanvas, 
                                                            TString namePad, 
                                                            Double_t* fPlottingRangeMeson,
                                                            TString dateDummy,
                                                            TString fMesonType, 
                                                            Int_t fRowPlot, 
                                                            Int_t fColumnPlot,
                                                            Int_t fStartBinPtRange, 
                                                            Int_t fNumberPtBins, 
                                                            Double_t* fRangeBinsPt, 
                                                            TString fDecayChannel, 
                                                            Bool_t fMonteCarloInfo, 
                                                            TString decayChannel, 
                                                            TString fDetectionChannel, 
                                                            TString fEnergy, 
                                                            Int_t mode){
    TCanvas *canvasDataFit          = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    canvasDataFit->SetTopMargin(0.00);
    canvasDataFit->SetBottomMargin(0.00);
    canvasDataFit->SetRightMargin(0.00);
    canvasDataFit->SetLeftMargin(0.00);

    TPad * padDataFit               = new TPad(namePad.Data(),"",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
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
    Int_t place                     = 0;
    for(Int_t iPt = fStartBinPtRange; iPt < fNumberPtBins; iPt++){
        Double_t startPt            = fRangeBinsPt[iPt];
        Double_t endPt              = fRangeBinsPt[iPt+1];

        place                       = place + 1; //give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataFit->cd(place);

            TString textAlice       = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo){
                textEvents          = "MC";
            } else {
                textEvents          = "Data";
            }
            Double_t nPixels        = 10;
            Double_t textHeight     = 0.08;
            if (padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) < padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1())){
                textHeight          = (Double_t)nPixels/padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) ;
            } else {
                textHeight          = (Double_t)nPixels/padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1());
            }
            Double_t startTextX     = 0.10;
            Double_t startTextY     = 0.9;
            Double_t differenceText = textHeight*1.25;

            TLatex *alice           = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
            TLatex *latexDate       = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy          = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
            TLatex *process         = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess      = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
            TLatex *events          = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

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
            Double_t nSignals       = 3;
            if (fHistoTrueSignalMixed != NULL) nSignals++;
            if (fHistoTrueSignalElectrons != NULL) nSignals++;
            TLegend* legendMC       = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+nSignals)*differenceText);
            legendMC->SetTextSize(textHeight);
            legendMC->SetTextFont(62);
            legendMC->SetLineColor(0);
            legendMC->SetLineWidth(0);
            legendMC->SetFillStyle(0);
            legendMC->SetFillColor(0);
            legendMC->SetMargin(0.15);
            Size_t markersize2      = fHistoTrueSignal[iPt-1]->GetMarkerSize();
            fHistoTrueSignal[iPt-1]->SetMarkerSize(3*markersize2);
            fHistoTrueSignalPhotons[iPt-1]->SetMarkerSize(3*markersize2);
            if (fHistoTrueSignalElectrons != NULL)fHistoTrueSignalElectrons[iPt-1]->SetMarkerSize(3*markersize2);
            fHistoTrueSignalConvPhotons[iPt-1]->SetMarkerSize(3*markersize2);
            if (fHistoTrueSignalMixed != NULL)fHistoTrueSignalMixed[iPt-1]->SetMarkerSize(3*markersize2);
            legendMC->AddEntry(fHistoTrueSignal[iPt-1],"validated meson","ep");
            if (mode == 4 || mode == 5){
                legendMC->AddEntry(fHistoTrueSignalPhotons[iPt-1],"val. #gamma#gamma","ep");
                if (fHistoTrueSignalElectrons != NULL)legendMC->AddEntry(fHistoTrueSignalElectrons[iPt-1],"val. e^{#pm}e^{#pm}","ep");
                legendMC->AddEntry(fHistoTrueSignalConvPhotons[iPt-1],"val. #gamma_{conv}#gamma_{conv}","ep");
                if (fHistoTrueSignalMixed != NULL) legendMC->AddEntry(fHistoTrueSignalMixed[iPt-1],"val. #gamma#gamma_{conv}","ep");
            } else if (mode == 2 || mode == 3) {
                legendMC->AddEntry(fHistoTrueSignalPhotons[iPt-1],"val. #gamma_{conv}#gamma","ep");
                if (fHistoTrueSignalElectrons != NULL)legendMC->AddEntry(fHistoTrueSignalElectrons[iPt-1],"val. #gamma_{conv}e^{#pm}","ep");
                legendMC->AddEntry(fHistoTrueSignalConvPhotons[iPt-1],"val. #gamma_{conv}#gamma_{conv}","ep");
            }
            legendMC->Draw();
            
        } else {

            padDataFit->cd(place);
            padDataFit->cd(place)->SetTopMargin(0.12);
            padDataFit->cd(place)->SetBottomMargin(0.15);
            padDataFit->cd(place)->SetRightMargin(0.02);
            int remaining           = (place-1)%fColumnPlot;
            if (remaining > 0) padDataFit->cd(place)->SetLeftMargin(0.15);
            else padDataFit->cd(place)->SetLeftMargin(0.25);

            cout << startPt << "-" << endPt <<endl;

            DrawGammaHisto( fHistoTrueSignal[iPt],
                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
            DrawGammaHisto( fHistoTrueSignal[iPt],
                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],2);

            fHistoTrueSignalPhotons[iPt]->SetMarkerColor(kRed+2);
            fHistoTrueSignalPhotons[iPt]->SetMarkerStyle(21);
            fHistoTrueSignalPhotons[iPt]->SetMarkerSize(0.3);
            fHistoTrueSignalPhotons[iPt]->SetLineColor(kRed+2);
            fHistoTrueSignalPhotons[iPt]->DrawCopy("same,p,e1");
            
            if (fHistoTrueSignalElectrons != NULL){
                fHistoTrueSignalElectrons[iPt]->SetMarkerColor(kBlue+2);
                fHistoTrueSignalElectrons[iPt]->SetMarkerStyle(33);
                fHistoTrueSignalElectrons[iPt]->SetMarkerSize(0.5);
                fHistoTrueSignalElectrons[iPt]->SetLineColor(kBlue+2);
                fHistoTrueSignalElectrons[iPt]->DrawCopy("same,p,e1");
            }
            
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
void PlotTrueInvMassSplittedInMergedInPtBins(   TH1D ** fHistoTrueSignal,
                                                TH1D** fHistoTrueSignalMerged, 
                                                TH1D** fHistoTrueSignalMergedPartConv, 
                                                TString namePlot, 
                                                TString nameCanvas, 
                                                TString namePad, 
                                                Double_t* fPlottingRangeMeson, 
                                                TString dateDummy, 
                                                TString fMesonType,  
                                                Int_t fRowPlot, 
                                                Int_t fColumnPlot, 
                                                Int_t fStartBinPtRange, 
                                                Int_t fNumberPtBins, 
                                                Double_t* fRangeBinsPt,
                                                TString fDecayChannel, 
                                                Bool_t fMonteCarloInfo,
                                                TString decayChannel,
                                                TString fDetectionChannel,
                                                TString fEnergy, 
                                                Int_t mode){
    TCanvas *canvasDataFit          = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    canvasDataFit->SetTopMargin(0.00);
    canvasDataFit->SetBottomMargin(0.00);
    canvasDataFit->SetRightMargin(0.00);
    canvasDataFit->SetLeftMargin(0.00);

    TPad * padDataFit               = new TPad(namePad.Data(),"",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
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
    Int_t place                     = 0;
    for(Int_t iPt = fStartBinPtRange; iPt < fNumberPtBins; iPt++){
        Double_t startPt            = fRangeBinsPt[iPt];
        Double_t endPt              = fRangeBinsPt[iPt+1];

        place                       = place + 1;//give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataFit->cd(place);

            TString textAlice       = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo){
                textEvents          = "MC";
            } else {
                textEvents          = "Data";
            }
            Double_t nPixels        = 10;
            Double_t textHeight     = 0.08;
            if (padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) < padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1())){
                textHeight          = (Double_t)nPixels/padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) ;
            } else {
                textHeight          = (Double_t)nPixels/padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1());
            }
            Double_t startTextX     = 0.10;
            Double_t startTextY     = 0.9;
            Double_t differenceText = textHeight*1.25;

            TLatex *alice           = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
            TLatex *latexDate       = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy          = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
            TLatex *process         = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess      = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
            TLatex *events          = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

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
            Double_t nSignals       = 3;
            TLegend* legendMC       = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+nSignals)*differenceText);
            legendMC->SetTextSize(textHeight);
            legendMC->SetTextFont(62);
            legendMC->SetLineColor(0);
            legendMC->SetLineWidth(0);
            legendMC->SetFillStyle(0);
            legendMC->SetFillColor(0);
            legendMC->SetMargin(0.15);
            Size_t markersize2      = fHistoTrueSignal[iPt-1]->GetMarkerSize();
            fHistoTrueSignal[iPt-1]->SetMarkerSize(3*markersize2);
            fHistoTrueSignalMerged[iPt-1]->SetMarkerSize(3*markersize2);
            fHistoTrueSignalMergedPartConv[iPt-1]->SetMarkerSize(3*markersize2);
            legendMC->AddEntry(fHistoTrueSignal[iPt-1],"validated meson","ep");
            if (mode == 4 || mode == 5){
                legendMC->AddEntry(fHistoTrueSignalMerged[iPt-1],"at least 1 E#it{M}_{merged}","ep");
                legendMC->AddEntry(fHistoTrueSignalMergedPartConv[iPt-1],"at least 1 E#it{M}_{merged}, part conv","ep");
            } else {
                legendMC->AddEntry(fHistoTrueSignalMerged[iPt-1],"val. #gamma_{conv}E#it{M}_{merged}","ep");
                legendMC->AddEntry(fHistoTrueSignalMergedPartConv[iPt-1],"val. #gamma_{conv}E#it{M}_{merged}, part conv","ep");
            }
            legendMC->Draw();
            
        } else {

            padDataFit->cd(place);
            padDataFit->cd(place)->SetTopMargin(0.12);
            padDataFit->cd(place)->SetBottomMargin(0.15);
            padDataFit->cd(place)->SetRightMargin(0.02);
            int remaining           = (place-1)%fColumnPlot;
            if (remaining > 0) padDataFit->cd(place)->SetLeftMargin(0.15);
            else padDataFit->cd(place)->SetLeftMargin(0.25);

            cout << startPt << "-" << endPt <<endl;

            DrawGammaHisto( fHistoTrueSignal[iPt],
                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
            DrawGammaHisto( fHistoTrueSignal[iPt],
                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
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
//________________________________________ Plot Invariant Mass Bin With GG contamination for Dalitz _______________________________
void PlotSignalInvMassW0TrueMesonInPtBins(  TH1D** fHistoMappingSignalInvMassW0TruePi0PtBinsPlot, 
                                            TString namePlot,
                                            TString nameCanvas,
                                            TString namePad,
                                            Double_t* fPlottingRangeMeson,
                                            TString dateDummy, 
                                            TString fMesonType,  
                                            Int_t fRowPlot, 
                                            Int_t fColumnPlot, 
                                            Int_t fStartBinPtRange,
                                            Int_t fNumberPtBins,
                                            Double_t* fRangeBinsPt,
                                            TString fDecayChannel,
                                            Bool_t fMonteCarloInfo,
                                            TString decayChannel,
                                            TString fDetectionChannel,
                                            TString fEnergy){
    TGaxis::SetMaxDigits(3);

    TCanvas * canvasDataSpectra     = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    canvasDataSpectra->SetTopMargin(0.00);
    canvasDataSpectra->SetBottomMargin(0.00);
    canvasDataSpectra->SetRightMargin(0.0);
    canvasDataSpectra->SetLeftMargin(0.00);

    TPad * padDataSpectra           = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    padDataSpectra->SetFillColor(0);
    padDataSpectra->GetFrame()->SetFillColor(0);
    padDataSpectra->SetBorderMode(0);
    padDataSpectra->SetLogy(0);
    padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
    padDataSpectra->Draw();

    cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
    cout << fMesonType.Data() << endl;

    Int_t place                     = 0;
    for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
        cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
        Double_t startPt            = fRangeBinsPt[iPt];
        Double_t endPt              = fRangeBinsPt[iPt+1];

        place                       = place + 1;//give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataSpectra->cd(place);

            TString textAlice       = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo){
                textEvents          = "MC";
            } else {
                textEvents          = "Data";
            }
            Double_t nPixels        = 13;
            Double_t textHeight     = 0.08;
            if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
                textHeight          = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
            } else {
                textHeight          = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
            }
            Double_t startTextX     = 0.10;
            Double_t startTextY     = 0.8;
            Double_t differenceText = textHeight*1.25;
            
            TLatex *alice           = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
            TLatex *latexDate       = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy          = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
            TLatex *process         = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess      = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
            TLatex *events          = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

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
            padDataSpectra->cd(place)->SetRightMargin(0.05);
            padDataSpectra->cd(place)->SetLeftMargin(0.15);

            DrawGammaHisto( fHistoMappingSignalInvMassW0TruePi0PtBinsPlot[iPt],
                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);

            }
    }
    canvasDataSpectra->Print(namePlot.Data());
    delete padDataSpectra;
    delete canvasDataSpectra;

}


void PlotInvMassW0TrueMesonInPtBins(    TH1D** fHistoMappingGGInvMassW0TruePi0PtBinsPlot, 
                                        TH1D** fHistoMappingBackNormInvMassPtBinPlot, 
                                        TString namePlot,
                                        TString nameCanvas,
                                        TString namePad, 
                                        Double_t* fPlottingRangeMeson, 
                                        TString dateDummy, 
                                        TString fMesonType,
                                        Int_t fRowPlot, 
                                        Int_t fColumnPlot, 
                                        Int_t fStartBinPtRange, 
                                        Int_t fNumberPtBins,
                                        Double_t* fRangeBinsPt, 
                                        TString fDecayChannel,
                                        Bool_t fMonteCarloInfo, 
                                        TString decayChannel, 
                                        TString fDetectionChannel,
                                        TString fEnergy){
    TGaxis::SetMaxDigits(3);

    TCanvas * canvasDataSpectra     = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    canvasDataSpectra->SetTopMargin(0.0);
    canvasDataSpectra->SetBottomMargin(0.0);
    canvasDataSpectra->SetRightMargin(0.0);
    canvasDataSpectra->SetLeftMargin(0.0);

    TPad * padDataSpectra           = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    padDataSpectra->SetFillColor(0);
    padDataSpectra->GetFrame()->SetFillColor(0);
    padDataSpectra->SetBorderMode(0);
    padDataSpectra->SetLogy(0);
    padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
    padDataSpectra->Draw();

    cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
    cout << fMesonType.Data() << endl;
    Int_t place                     = 0;
    for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
        cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
        Double_t startPt            = fRangeBinsPt[iPt];
        Double_t endPt              = fRangeBinsPt[iPt+1];

        place = place + 1;//give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataSpectra->cd(place);

            TString textAlice       = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo){
                textEvents          = "MC";
            } else {
                textEvents          = "Data";
            }
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
            TLatex *latexDate       = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy          = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
            TLatex *process         = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess      = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
            TLatex *events          = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

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

            TLegend* legendData     = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+2.)*differenceText);
            legendData->SetTextSize(textHeight);
            legendData->SetTextFont(62);
            legendData->SetFillColor(0);
            legendData->SetFillStyle(0);
            legendData->SetLineWidth(0);
            legendData->SetLineColor(0);
            legendData->SetMargin(0.15);
            Size_t markersize       = fHistoMappingGGInvMassW0TruePi0PtBinsPlot[iPt-1]->GetMarkerSize();
            fHistoMappingGGInvMassW0TruePi0PtBinsPlot[iPt-1]->SetMarkerSize(3*markersize);
            legendData->AddEntry(fHistoMappingGGInvMassW0TruePi0PtBinsPlot[iPt-1],"same evt. #it{M}_{#gamma#gamma} - True #pi^{0} (#gamma#gamma + Dalitz channels)","ep");
            Size_t linesize         = fHistoMappingBackNormInvMassPtBinPlot[iPt-1]->GetLineWidth();
            fHistoMappingBackNormInvMassPtBinPlot[iPt-1]->SetLineWidth(5*linesize);
            legendData->AddEntry(fHistoMappingBackNormInvMassPtBinPlot[iPt-1],"mixed evt. #it{M}_{#gamma#gamma}","l");
            legendData->Draw();

        } else {

            padDataSpectra->cd(place);
            padDataSpectra->cd(place)->SetTopMargin(0.12);
            padDataSpectra->cd(place)->SetBottomMargin(0.15);
            padDataSpectra->cd(place)->SetRightMargin(0.02);
            padDataSpectra->cd(place)->SetLeftMargin(0.15);

            DrawGammaHisto( fHistoMappingGGInvMassW0TruePi0PtBinsPlot[iPt],
                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);

            DrawGammaHisto( fHistoMappingBackNormInvMassPtBinPlot[iPt],
                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],1);
        }
    }
    canvasDataSpectra->Print(namePlot.Data());
    delete padDataSpectra;
    delete canvasDataSpectra;
}


//________________________________________ Plot Invariant Mass Bin With Secondary Contribution _______________________________
void PlotInvMassDoubleCountingInPtBins( TH1D** fHistoMappingGGInvMassPtBinPlot,
                                        TH1D** fHistoMappingDCInvMassPtBinPlot, 
                                        TString namePlot, 
                                        TString nameCanvas, 
                                        TString namePad, 
                                        Double_t* fPlottingRangeMeson,
                                        TString dateDummy,
                                        TString fMesonType,  
                                        Int_t fRowPlot,
                                        Int_t fColumnPlot,
                                        Int_t fStartBinPtRange, 
                                        Int_t fNumberPtBins, 
                                        Double_t* fRangeBinsPt, 
                                        TString fDecayChannel, 
                                        Bool_t fMonteCarloInfo, 
                                        TString decayChannel,
                                        TString fDetectionChannel, 
                                        TString fEnergy){
    TGaxis::SetMaxDigits(3);

    TCanvas* canvasDataSpectra      = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    canvasDataSpectra->SetTopMargin(0.0);
    canvasDataSpectra->SetBottomMargin(0.0);
    canvasDataSpectra->SetRightMargin(0.0);
    canvasDataSpectra->SetLeftMargin(0.0);

    TPad* padDataSpectra            = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    padDataSpectra->SetFillColor(0);
    padDataSpectra->GetFrame()->SetFillColor(0);
    padDataSpectra->SetBorderMode(0);
    padDataSpectra->SetLogy(0);
    padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
    padDataSpectra->Draw();

    cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
    cout << fMesonType.Data() << endl;
    Int_t place                     = 0;
    for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
        cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
        Double_t startPt            = fRangeBinsPt[iPt];
        Double_t endPt              = fRangeBinsPt[iPt+1];

        place                       = place + 1; //give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataSpectra->cd(place);

            TString textAlice       = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo){
                textEvents          = "MC";
            } else {
                textEvents          = "Data";
            }
            Double_t nPixels        = 13;
            Double_t textHeight     = 0.08;
            if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
                textHeight          = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
            } else {
                textHeight          = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
            }
            Double_t startTextX     = 0.10;
            Double_t startTextY     = 0.8;
            Double_t differenceText = textHeight*1.25;

            TLatex *alice           = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
            TLatex *latexDate       = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy          = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
            TLatex *process         = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess      = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
            TLatex *events          = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

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
            int remaining           = (place-1)%fColumnPlot;
            if (remaining > 0) padDataSpectra->cd(place)->SetLeftMargin(0.15);
            else padDataSpectra->cd(place)->SetLeftMargin(0.25);

            DrawGammaHisto( fHistoMappingGGInvMassPtBinPlot[iPt],
                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
            cout << fHistoMappingGGInvMassPtBinPlot[iPt]->GetNbinsX() << "\t" << fHistoMappingDCInvMassPtBinPlot[iPt]->GetNbinsX() << endl;
            DrawGammaHisto( fHistoMappingDCInvMassPtBinPlot[iPt],
                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],1);
            cout << fHistoMappingDCInvMassPtBinPlot[iPt]->GetNbinsX() << endl;
            
            
        }
    }
    canvasDataSpectra->Print(namePlot.Data());
    delete padDataSpectra;
    delete canvasDataSpectra;

}


//__________________________________________ Plotting all Invariant Mass bins _______________________________________________
void PlotInvMassMergedInPtBins( TH1D**      fHistoInvMassPtBinPlot, 
                                TString     namePlot, 
                                TString     nameCanvas, 
                                TString     namePad, 
                                Double_t*   fPlottingRangeMeson, 
                                TString     dateDummy, 
                                TString     fMesonType, 
                                Int_t       fRowPlot, 
                                Int_t       fColumnPlot, 
                                Int_t       fStartBinPtRange, 
                                Int_t       fNumberPtBins, 
                                Double_t*   fRangeBinsPt, 
                                TString     fDecayChannel, 
                                Bool_t      fMonteCarloInfo, 
                                TString     decayChannel,  
                                TString     fDetectionChannel, 
                                TString     fEnergy
                              ){
    TGaxis::SetMaxDigits(3);
    TCanvas * canvasDataSpectra     = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    DrawGammaCanvasSettings(canvasDataSpectra, 0.0, 0., 0., 0.);

    TPad * padDataSpectra           = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    padDataSpectra->SetFillColor(0);
    padDataSpectra->GetFrame()->SetFillColor(0);
    padDataSpectra->SetBorderMode(0);
    padDataSpectra->SetLogy(0);
    padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
    padDataSpectra->Draw();

    cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
    cout << fMesonType.Data() << endl;
    Int_t place                     = 0;
    for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
        cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
        Double_t startPt            = fRangeBinsPt[iPt];
        Double_t endPt              = fRangeBinsPt[iPt+1];

        place                       = place + 1;  //give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataSpectra->cd(place);

            cout << "entered ALICE plotting" << endl;
            TString textAlice       = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo){
                textEvents          = "MC";
            } else {
                textEvents          = "Data";
            }
            Double_t nPixels = 13;
            Double_t textHeight = 0.08;
            if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
                textHeight          = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
            } else {
                textHeight          = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
            }   
            
            Double_t startTextX     = 0.10;
            Double_t startTextY     = 0.9;
            Double_t differenceText = textHeight*1.25;

            TLatex *alice           = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
            TLatex *latexDate       = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy          = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
            TLatex *process         = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess      = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
            TLatex *events          = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

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
//             cout << "place: " << place << "\t "<< startPt << " < #it{p}_{T} < " << endPt  << endl;
            int remaining           = (place-1)%fColumnPlot;
            if (remaining > 0) padDataSpectra->cd(place)->SetLeftMargin(0.15);
            else padDataSpectra->cd(place)->SetLeftMargin(0.25);
//             cout << "here" << endl;
            DrawGammaHisto( fHistoInvMassPtBinPlot[iPt],
                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0,0.6);
            cout << "here" << endl;
        }
    }
    canvasDataSpectra->Print(namePlot.Data());
    delete padDataSpectra;
    delete canvasDataSpectra;
}

//__________________________________________ Plotting all M02 bins _______________________________________________
void PlotM02MergedInPtBins(     TH1D**      fHistoM02PtBinPlot, 
                                TString     namePlot, 
                                TString     nameCanvas, 
                                TString     namePad, 
                                Double_t*   fPlottingM02Range, 
                                TString     dateDummy, 
                                TString     fMesonType, 
                                Int_t       fRowPlot, 
                                Int_t       fColumnPlot, 
                                Int_t       fStartBinPtRange, 
                                Int_t       fNumberPtBins, 
                                Double_t*   fRangeBinsPt, 
                                TString     fDecayChannel, 
                                Bool_t      fMonteCarloInfo, 
                                TString     decayChannel,  
                                TString     fDetectionChannel, 
                                TString     fEnergy
                              ){
    TGaxis::SetMaxDigits(3);
    TCanvas * canvasDataSpectra         = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    DrawGammaCanvasSettings(canvasDataSpectra, 0.0, 0., 0., 0.);

    TPad * padDataSpectra               = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    padDataSpectra->SetFillColor(0);
    padDataSpectra->GetFrame()->SetFillColor(0);
    padDataSpectra->SetBorderMode(0);
    padDataSpectra->SetLogy(1);
    padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
    padDataSpectra->Draw();

    cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
    cout << fMesonType.Data() << endl;
    Int_t place                         = 0;
    for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
        cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
        Double_t startPt                = fRangeBinsPt[iPt];
        Double_t endPt                  = fRangeBinsPt[iPt+1];

        place                           = place + 1;  //give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataSpectra->cd(place);

            cout << "entered ALICE plotting" << endl;
            TString textAlice           = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo){
                textEvents              = "MC";
            } else {
                textEvents              = "Data";
            }
            Double_t nPixels            = 13;
            Double_t textHeight         = 0.08;
            if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
                textHeight              = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
            } else {
                textHeight              = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
            }   
            
            Double_t startTextX         = 0.10;
            Double_t startTextY         = 0.9;
            Double_t differenceText     = textHeight*1.25;

            TLatex *alice               = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
            TLatex *latexDate           = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy              = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
            TLatex *process             = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess          = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
            TLatex *events              = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

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
//             cout << "place: " << place << "\t "<< startPt << " < #it{p}_{T} < " << endPt  << endl;
            int remaining               = (place-1)%fColumnPlot;
            if (remaining > 0) padDataSpectra->cd(place)->SetLeftMargin(0.15);
            else padDataSpectra->cd(place)->SetLeftMargin(0.25);
            if (fHistoM02PtBinPlot[iPt]){
                    DrawGammaHisto( fHistoM02PtBinPlot[iPt],
                                    Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}", startPt, endPt),
                                    "#it{#sigma}_{long}^{2}", Form("dN_{%s}/d#it{#sigma}_{long}^{2}", decayChannel.Data()),
                                    fPlottingM02Range[0], fPlottingM02Range[1], 0, 0.7);
            } else {
                 continue;
            }
        }
    }
    canvasDataSpectra->Print(namePlot.Data());
    delete padDataSpectra;
    delete canvasDataSpectra;
}


//__________________________________________ Plotting all Invariant Mass bins _______________________________________________
void PlotInvMassMergedMCInPtBins(   TH1D**      fHistoMergedRecPtBinPlot, 
                                    TH1D**      fHistoMergedValPtBinPlot, 
                                    TH1D**      fHistoMergedPartConvValPtBinPlot, 
                                    TH1D**      fHistoBGPtBinPlot,
                                    TH1D**      fHistoValGammaPtBinPlot,
                                    TH1D**      fHistoValElectronPtBinPlot,
                                    TString     namePlot, 
                                    TString     nameCanvas, 
                                    TString     namePad, 
                                    Double_t*   fPlottingRangeMeson, 
                                    TString     dateDummy, 
                                    TString     fMesonType, 
                                    Int_t       fRowPlot, 
                                    Int_t       fColumnPlot, 
                                    Int_t       fStartBinPtRange, 
                                    Int_t       fNumberPtBins, 
                                    Double_t*   fRangeBinsPt, 
                                    TString     fDecayChannel, 
                                    Bool_t      fMonteCarloInfo, 
                                    TString     decayChannel,  
                                    TString     fDetectionChannel, 
                                    TString     fEnergy
                              ){
    TGaxis::SetMaxDigits(3);
    TCanvas * canvasDataSpectra         = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    DrawGammaCanvasSettings(canvasDataSpectra, 0.0, 0., 0., 0.);

    TPad * padDataSpectra               = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    padDataSpectra->SetFillColor(0);
    padDataSpectra->GetFrame()->SetFillColor(0);
    padDataSpectra->SetBorderMode(0);
    padDataSpectra->SetLogy(0);
    padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
    padDataSpectra->Draw();

    cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
    cout << fMesonType.Data() << endl;
    Int_t place                         = 0;
    for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
        cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
        Double_t startPt                = fRangeBinsPt[iPt];
        Double_t endPt                  = fRangeBinsPt[iPt+1];

        place                           = place + 1;  //give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataSpectra->cd(place);

            cout << "entered ALICE plotting" << endl;
            TString textAlice           = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo){
                textEvents              = "MC";
            } else {
                textEvents              = "Data";
            }
            Double_t nPixels            = 13;
            Double_t textHeight         = 0.08;
            if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
                textHeight              = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
            } else {
                textHeight              = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
            }   
            
            Double_t startTextX         = 0.10;
            Double_t startTextY         = 0.9;
            Double_t differenceText     = textHeight*1.25;

            TLatex *alice               = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
            TLatex *latexDate           = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy              = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
            TLatex *process             = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess          = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
            TLatex *events              = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

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

            Int_t nLegend               = 0;
            if (fHistoMergedRecPtBinPlot) nLegend++;
            if (fHistoMergedValPtBinPlot) nLegend++;
            if (fHistoMergedPartConvValPtBinPlot) nLegend++;
            if (fHistoValGammaPtBinPlot) nLegend++;
            if (fHistoValElectronPtBinPlot) nLegend++;
            if (fHistoBGPtBinPlot) nLegend++;
            
            TLegend* legend             = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+nLegend)*differenceText);
            legend->SetTextSize(textHeight);
            legend->SetTextFont(62);
            legend->SetFillColor(0);
            legend->SetFillStyle(0);
            legend->SetLineWidth(0);
            legend->SetLineColor(0);
            legend->SetMargin(0.15);
            Size_t markersize           = 0;
            Size_t linesize             = 0;
            if (fHistoMergedRecPtBinPlot){
                if (fHistoMergedRecPtBinPlot[iPt]){
                    markersize          = fHistoMergedRecPtBinPlot[iPt]->GetMarkerSize();
                    fHistoMergedRecPtBinPlot[iPt]->SetMarkerSize(1*markersize);
                    legend->AddEntry(fHistoMergedRecPtBinPlot[iPt],"merged clus. accept","ep");
                }
            }
            if (fHistoMergedValPtBinPlot){
                if (fHistoMergedValPtBinPlot[iPt]){
                    markersize          = fHistoMergedValPtBinPlot[iPt]->GetMarkerSize();
                    fHistoMergedValPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoMergedValPtBinPlot[iPt],"merged clus. val.","ep");
                }
            }
            if (fHistoMergedPartConvValPtBinPlot){
                if (fHistoMergedPartConvValPtBinPlot[iPt]){
                    markersize          = fHistoMergedPartConvValPtBinPlot[iPt]->GetMarkerSize();
                    fHistoMergedPartConvValPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoMergedPartConvValPtBinPlot[iPt],"merged clus. val. p. conv","ep");
                }
            }
            if (fHistoValGammaPtBinPlot){
                if (fHistoValGammaPtBinPlot[iPt]){
                    markersize          = fHistoValGammaPtBinPlot[iPt]->GetMarkerSize();
                    fHistoValGammaPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoValGammaPtBinPlot[iPt],"#gamma","ep");
                }
            }
            if (fHistoValElectronPtBinPlot){
                if (fHistoValElectronPtBinPlot[iPt]){
                    markersize          = fHistoValElectronPtBinPlot[iPt]->GetMarkerSize();
                    fHistoValElectronPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoValElectronPtBinPlot[iPt],"e^{#pm}","ep");
                }
            }
            if (fHistoBGPtBinPlot){
                if (fHistoBGPtBinPlot[iPt]){
                    linesize            = fHistoBGPtBinPlot[iPt]->GetLineWidth();
                    fHistoBGPtBinPlot[iPt]->SetLineWidth(2*linesize);
                    legend->AddEntry(fHistoBGPtBinPlot[iPt],"background","l");
                }
            }
            legend->Draw(); 
            
        } else {            
            padDataSpectra->cd(place);
            padDataSpectra->cd(place)->SetTopMargin(0.12);
            padDataSpectra->cd(place)->SetBottomMargin(0.15);
            padDataSpectra->cd(place)->SetRightMargin(0.02);
//             cout << "place: " << place << "\t "<< startPt << " < #it{p}_{T} < " << endPt  << endl;
            int remaining               = (place-1)%fColumnPlot;
            if (remaining > 0) padDataSpectra->cd(place)->SetLeftMargin(0.15);
            else padDataSpectra->cd(place)->SetLeftMargin(0.25);
//             cout << "here" << endl;
            if (fHistoMergedRecPtBinPlot){
                if (fHistoMergedRecPtBinPlot[iPt]){
                    DrawGammaHisto( fHistoMergedRecPtBinPlot[iPt],
                                    Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                    Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                    fPlottingRangeMeson[0],fPlottingRangeMeson[1],0,0.6);
                } else {
                    continue;
                }
            } else {
                cout << "ABORT: plotting no default histos" << endl;
                return;
            }
            if (fHistoMergedValPtBinPlot){
                if (fHistoMergedValPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoMergedValPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0,kBlue+1, 24);
                }
            }
            if (fHistoMergedPartConvValPtBinPlot){
                if (fHistoMergedPartConvValPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoMergedPartConvValPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0,kGreen+2, 24);
                }
            }
            if (fHistoValGammaPtBinPlot){
                if (fHistoValGammaPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoValGammaPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0,kRed+2, 25);
                }
            }
            if (fHistoValElectronPtBinPlot){
                if (fHistoValElectronPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoValElectronPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0,807, 25);
                }
            }
            if (fHistoBGPtBinPlot){
                if (fHistoBGPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoBGPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0,kGray+2);
                }
            }
        }
    }
    canvasDataSpectra->Print(namePlot.Data());
    delete padDataSpectra;
    delete canvasDataSpectra;
}

//__________________________________________ Plotting all M02 bins _______________________________________________
void PlotM02ExoticsInEBins(     TH1D**      fHistoM02PtBinPlot, 
                                TString     namePlot, 
                                TString     nameCanvas, 
                                TString     namePad, 
                                Double_t*   fPlottingM02Range, 
                                TString     dateDummy, 
                                TString     fMesonType, 
                                Int_t       fRowPlot, 
                                Int_t       fColumnPlot, 
                                Int_t       fStartBinPtRange, 
                                Int_t       fNumberPtBins, 
                                Double_t*   fRangeBinsPt, 
                                TString     fDecayChannel, 
                                Bool_t      fMonteCarloInfo, 
                                TString     decayChannel,  
                                TString     fDetectionChannel, 
                                TString     fEnergy
                              ){
    TGaxis::SetMaxDigits(3);
    TCanvas * canvasDataSpectra         = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    DrawGammaCanvasSettings(canvasDataSpectra, 0.0, 0., 0., 0.);

    TPad * padDataSpectra               = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    padDataSpectra->SetFillColor(0);
    padDataSpectra->GetFrame()->SetFillColor(0);
    padDataSpectra->SetBorderMode(0);
    padDataSpectra->SetLogy(1);
    padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
    padDataSpectra->Draw();

    cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
    cout << fMesonType.Data() << endl;
    Int_t place                         = 0;
    for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
        cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
        Double_t startPt                = fRangeBinsPt[iPt];
        Double_t endPt                  = fRangeBinsPt[iPt+1];

        place                           = place + 1;  //give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataSpectra->cd(place);

            cout << "entered ALICE plotting" << endl;
            TString textAlice           = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo){
                textEvents              = "MC";
            } else {
                textEvents              = "Data";
            }
            Double_t nPixels            = 13;
            Double_t textHeight         = 0.08;
            if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
                textHeight              = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
            } else {
                textHeight              = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
            }   
            
            Double_t startTextX         = 0.10;
            Double_t startTextY         = 0.9;
            Double_t differenceText     = textHeight*1.25;

            TLatex *alice               = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
            TLatex *latexDate           = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy              = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
            TLatex *process             = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess          = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
            TLatex *events              = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

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
//             cout << "place: " << place << "\t "<< startPt << " < #it{p}_{T} < " << endPt  << endl;
            int remaining               = (place-1)%fColumnPlot;
            if (remaining > 0) padDataSpectra->cd(place)->SetLeftMargin(0.15);
            else padDataSpectra->cd(place)->SetLeftMargin(0.25);
            if (fHistoM02PtBinPlot[iPt]){
                    DrawGammaHisto( fHistoM02PtBinPlot[iPt],
                                    Form("%3.2f GeV/#it{c} < #it{E} < %3.2f GeV", startPt, endPt),
                                    "#it{#sigma}_{long}^{2}", Form("dN_{%s}/d#it{#sigma}_{long}^{2}", decayChannel.Data()),
                                    0, 0.5, 0.0001, 0.7);
            } else {
                 continue;
            }
        }
    }
    canvasDataSpectra->Print(namePlot.Data());
    delete padDataSpectra;
    delete canvasDataSpectra;
}


//__________________________________________ Plotting all Invariant Mass bins _______________________________________________
void PlotInvMassMergedTrueInPtBins( TH1D**      fHistoMergedValPtBinPlot, 
                                    TH1D**      fHistoMergedPartConvValPtBinPlot, 
                                    TH1D**      fHistoPi0ValPtBinPlot,
                                    TH1D**      fHistoPi0PartConvValPtBinPlot,
                                    TH1D**      fHistoEtaValPtBinPlot,
                                    TH1D**      fHistoEtaPartConvValPtBinPlot,
                                    TString     namePlot, 
                                    TString     nameCanvas, 
                                    TString     namePad, 
                                    Double_t*   fPlottingRangeMeson, 
                                    TString     dateDummy, 
                                    TString     fMesonType, 
                                    Int_t       fRowPlot, 
                                    Int_t       fColumnPlot, 
                                    Int_t       fStartBinPtRange, 
                                    Int_t       fNumberPtBins, 
                                    Double_t*   fRangeBinsPt, 
                                    TString     fDecayChannel, 
                                    Bool_t      fMonteCarloInfo, 
                                    TString     decayChannel,  
                                    TString     fDetectionChannel, 
                                    TString     fEnergy
                              ){
    TGaxis::SetMaxDigits(3);
    TCanvas * canvasDataSpectra         = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    DrawGammaCanvasSettings(canvasDataSpectra, 0.0, 0., 0., 0.);

    TPad * padDataSpectra               = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    padDataSpectra->SetFillColor(0);
    padDataSpectra->GetFrame()->SetFillColor(0);
    padDataSpectra->SetBorderMode(0);
    padDataSpectra->SetLogy(0);
    padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
    padDataSpectra->Draw();

    cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
    cout << fMesonType.Data() << endl;
    Int_t place                         = 0;
    for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
        cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
        Double_t startPt                = fRangeBinsPt[iPt];
        Double_t endPt                  = fRangeBinsPt[iPt+1];

        place                           = place + 1;  //give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataSpectra->cd(place);

            cout << "entered ALICE plotting" << endl;
            TString textAlice           = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo){
                textEvents              = "MC";
            } else {
                textEvents              = "Data";
            }
            Double_t nPixels            = 13;
            Double_t textHeight         = 0.08;
            if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
                textHeight              = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
            } else {
                textHeight              = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
            }   
            
            Double_t startTextX         = 0.10;
            Double_t startTextY         = 0.9;
            Double_t differenceText     = textHeight*1.25;

            TLatex *alice               = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
            TLatex *latexDate           = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy              = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
            TLatex *process             = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess          = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
            TLatex *events              = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

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
            
            Int_t nLegend               = 0;
            if (fHistoMergedValPtBinPlot) nLegend++;
            if (fHistoMergedPartConvValPtBinPlot) nLegend++;
            if (fHistoPi0ValPtBinPlot) nLegend++;
            if (fHistoPi0PartConvValPtBinPlot) nLegend++;
            if (fHistoEtaValPtBinPlot) nLegend++;
            if (fHistoEtaPartConvValPtBinPlot) nLegend++;
            
            TLegend* legend             = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+nLegend)*differenceText);
            legend->SetTextSize(textHeight);
            legend->SetTextFont(62);
            legend->SetFillColor(0);
            legend->SetFillStyle(0);
            legend->SetLineWidth(0);
            legend->SetLineColor(0);
            legend->SetMargin(0.15);
            Size_t markersize           = 0;
            Size_t linesize             = 0;
            if (fHistoMergedValPtBinPlot){
                if (fHistoMergedValPtBinPlot[iPt]){
                    markersize          = fHistoMergedValPtBinPlot[iPt]->GetMarkerSize();
                    fHistoMergedValPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoMergedValPtBinPlot[iPt],"merged clus.","ep");
                }
            }
            if (fHistoMergedPartConvValPtBinPlot){
                if (fHistoMergedPartConvValPtBinPlot[iPt]){
                    markersize          = fHistoMergedPartConvValPtBinPlot[iPt]->GetMarkerSize();
                    fHistoMergedPartConvValPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoMergedPartConvValPtBinPlot[iPt],"merged clus. p. conv","ep");
                }
            }
            if (fHistoPi0ValPtBinPlot){
                if (fHistoPi0ValPtBinPlot[iPt]){
                    markersize          = fHistoPi0ValPtBinPlot[iPt]->GetMarkerSize();
                    fHistoPi0ValPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoPi0ValPtBinPlot[iPt],"merged #pi^{0}","ep");
                }
            }
            if (fHistoPi0PartConvValPtBinPlot){
                if (fHistoPi0PartConvValPtBinPlot[iPt]){
                    markersize          = fHistoPi0PartConvValPtBinPlot[iPt]->GetMarkerSize();
                    fHistoPi0PartConvValPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoPi0PartConvValPtBinPlot[iPt],"merged #pi^{0} p. conv","ep");
                }
            }
            if (fHistoEtaValPtBinPlot){
                if (fHistoEtaValPtBinPlot[iPt]){
                    markersize          = fHistoEtaValPtBinPlot[iPt]->GetMarkerSize();
                    fHistoEtaValPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoEtaValPtBinPlot[iPt],"merged #eta","ep");
                }
            }
            if (fHistoEtaPartConvValPtBinPlot){
                if (fHistoEtaPartConvValPtBinPlot[iPt]){
                    markersize          = fHistoEtaPartConvValPtBinPlot[iPt]->GetMarkerSize();
                    fHistoEtaPartConvValPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoEtaPartConvValPtBinPlot[iPt],"merged #eta p. conv","ep");
                }
            }
            legend->Draw(); 
        } else {            
            padDataSpectra->cd(place);
            padDataSpectra->cd(place)->SetTopMargin(0.12);
            padDataSpectra->cd(place)->SetBottomMargin(0.15);
            padDataSpectra->cd(place)->SetRightMargin(0.02);
//             cout << "place: " << place << "\t "<< startPt << " < #it{p}_{T} < " << endPt  << endl;
            int remaining               = (place-1)%fColumnPlot;
            if (remaining > 0) padDataSpectra->cd(place)->SetLeftMargin(0.15);
            else padDataSpectra->cd(place)->SetLeftMargin(0.25);
//             cout << "here" << endl;
            if (fHistoMergedValPtBinPlot){
                if (fHistoMergedValPtBinPlot[iPt]){
                    DrawGammaHisto( fHistoMergedValPtBinPlot[iPt],
                                    Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                    Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                    fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
                } else {
                    continue;
                }
            } else {
                cout << "ABORT: plotting no default histos" << endl;
                return;
            }
            if (fHistoMergedPartConvValPtBinPlot){
                if (fHistoMergedPartConvValPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoMergedPartConvValPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0,kGray+2, 20);
                }
            }
            if (fHistoPi0ValPtBinPlot){
                if (fHistoPi0ValPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoPi0ValPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0,kBlue+1, 24);
                }
            }
            if (fHistoPi0PartConvValPtBinPlot){
                if (fHistoPi0PartConvValPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoPi0PartConvValPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0,kAzure-4, 24);
                }
            }
            if (fHistoEtaValPtBinPlot){
                if (fHistoEtaValPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoEtaValPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0,kGreen+2, 25);
                }
            }
            if (fHistoEtaPartConvValPtBinPlot){
                if (fHistoEtaPartConvValPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoEtaPartConvValPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                                            Form("#it{M}_{%s} (GeV/#it{c}^{2})",decayChannel.Data()), Form("dN_{%s}/d#it{M}_{%s}",decayChannel.Data(), decayChannel.Data()),
                                            fPlottingRangeMeson[0],fPlottingRangeMeson[1],0,kSpring+5, 25);
                }
            }
        }
    }
    canvasDataSpectra->Print(namePlot.Data());
    delete padDataSpectra;
    delete canvasDataSpectra;
}

//__________________________________________ Plotting all Invariant Mass bins _______________________________________________
void PlotM02MergedMCInPtBins(       TH1D**      fHistoMergedRecPtBinPlot, 
                                    TH1D**      fHistoMergedValPtBinPlot, 
                                    TH1D**      fHistoMergedPartConvValPtBinPlot, 
                                    TH1D**      fHistoBGPtBinPlot,
                                    TH1D**      fHistoValGammaPtBinPlot,
                                    TH1D**      fHistoValElectronPtBinPlot,
                                    TString     namePlot, 
                                    TString     nameCanvas, 
                                    TString     namePad, 
                                    Double_t*   fPlottingM02Range, 
                                    TString     dateDummy, 
                                    TString     fMesonType, 
                                    Int_t       fRowPlot, 
                                    Int_t       fColumnPlot, 
                                    Int_t       fStartBinPtRange, 
                                    Int_t       fNumberPtBins, 
                                    Double_t*   fRangeBinsPt, 
                                    TString     fDecayChannel, 
                                    Bool_t      fMonteCarloInfo, 
                                    TString     decayChannel,  
                                    TString     fDetectionChannel, 
                                    TString     fEnergy
                              ){
    TGaxis::SetMaxDigits(3);
    TCanvas * canvasDataSpectra         = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    DrawGammaCanvasSettings(canvasDataSpectra, 0.0, 0., 0., 0.);

    TPad * padDataSpectra               = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    padDataSpectra->SetFillColor(0);
    padDataSpectra->GetFrame()->SetFillColor(0);
    padDataSpectra->SetBorderMode(0);
    padDataSpectra->SetLogy(1);
    padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
    padDataSpectra->Draw();

    cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
    cout << fMesonType.Data() << endl;
    Int_t place                         = 0;
    for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
        cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
        Double_t startPt                = fRangeBinsPt[iPt];
        Double_t endPt                  = fRangeBinsPt[iPt+1];

        place                           = place + 1;  //give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataSpectra->cd(place);

            cout << "entered ALICE plotting" << endl;
            TString textAlice           = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo){
                textEvents              = "MC";
            } else {
                textEvents              = "Data";
            }
            Double_t nPixels            = 13;
            Double_t textHeight         = 0.08;
            if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
                textHeight              = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
            } else {
                textHeight              = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
            }   
            
            Double_t startTextX         = 0.10;
            Double_t startTextY         = 0.9;
            Double_t differenceText     = textHeight*1.25;

            TLatex *alice               = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
            TLatex *latexDate           = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy              = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
            TLatex *process             = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess          = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
            TLatex *events              = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

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

            Int_t nLegend               = 0;
            if (fHistoMergedRecPtBinPlot) nLegend++;
            if (fHistoMergedValPtBinPlot) nLegend++;
            if (fHistoMergedPartConvValPtBinPlot) nLegend++;
            if (fHistoValGammaPtBinPlot) nLegend++;
            if (fHistoValElectronPtBinPlot) nLegend++;
            if (fHistoBGPtBinPlot) nLegend++;
            
            TLegend* legend             = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+nLegend)*differenceText);
            legend->SetTextSize(textHeight);
            legend->SetTextFont(62);
            legend->SetFillColor(0);
            legend->SetFillStyle(0);
            legend->SetLineWidth(0);
            legend->SetLineColor(0);
            legend->SetMargin(0.15);
            Size_t markersize           = 0;
            Size_t linesize             = 0;
            if (fHistoMergedRecPtBinPlot){
                if (fHistoMergedRecPtBinPlot[iPt]){
                    markersize          = fHistoMergedRecPtBinPlot[iPt]->GetMarkerSize();
                    fHistoMergedRecPtBinPlot[iPt]->SetMarkerSize(1*markersize);
                    legend->AddEntry(fHistoMergedRecPtBinPlot[iPt],"merged clus. accept","ep");
                }
            }
            if (fHistoMergedValPtBinPlot){
                if (fHistoMergedValPtBinPlot[iPt]){
                    markersize          = fHistoMergedValPtBinPlot[iPt]->GetMarkerSize();
                    fHistoMergedValPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoMergedValPtBinPlot[iPt],"merged clus. val.","ep");
                }
            }
            if (fHistoMergedPartConvValPtBinPlot){
                if (fHistoMergedPartConvValPtBinPlot[iPt]){
                    markersize          = fHistoMergedPartConvValPtBinPlot[iPt]->GetMarkerSize();
                    fHistoMergedPartConvValPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoMergedPartConvValPtBinPlot[iPt],"merged clus. val. p. conv","ep");
                }
            }
            if (fHistoValGammaPtBinPlot){
                if (fHistoValGammaPtBinPlot[iPt]){
                    markersize          = fHistoValGammaPtBinPlot[iPt]->GetMarkerSize();
                    fHistoValGammaPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoValGammaPtBinPlot[iPt],"#gamma","ep");
                }
            }
            if (fHistoValElectronPtBinPlot){
                if (fHistoValElectronPtBinPlot[iPt]){
                    markersize          = fHistoValElectronPtBinPlot[iPt]->GetMarkerSize();
                    fHistoValElectronPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoValElectronPtBinPlot[iPt],"e^{#pm}","ep");
                }
            }
            if (fHistoBGPtBinPlot){
                if (fHistoBGPtBinPlot[iPt]){
                    linesize            = fHistoBGPtBinPlot[iPt]->GetLineWidth();
                    fHistoBGPtBinPlot[iPt]->SetLineWidth(2*linesize);
                    legend->AddEntry(fHistoBGPtBinPlot[iPt],"background","l");
                }
            }
            legend->Draw(); 
            
        } else {            
            padDataSpectra->cd(place);
            padDataSpectra->cd(place)->SetTopMargin(0.12);
            padDataSpectra->cd(place)->SetBottomMargin(0.15);
            padDataSpectra->cd(place)->SetRightMargin(0.02);
//             cout << "place: " << place << "\t "<< startPt << " < #it{p}_{T} < " << endPt  << endl;
            int remaining               = (place-1)%fColumnPlot;
            if (remaining > 0) padDataSpectra->cd(place)->SetLeftMargin(0.15);
            else padDataSpectra->cd(place)->SetLeftMargin(0.25);
//             cout << "here" << endl;
            if (fHistoMergedRecPtBinPlot){
                if (fHistoMergedRecPtBinPlot[iPt]){
                    DrawGammaHisto( fHistoMergedRecPtBinPlot[iPt],
                                    Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}", startPt, endPt),
                                    "#it{#sigma}_{long}^{2}", Form("dN_{%s}/d#it{#sigma}_{long}^{2}", decayChannel.Data()),
                                    fPlottingM02Range[0],fPlottingM02Range[1],0,0.8);
                } else {
                    continue;
                }
            } else {
                cout << "ABORT: plotting no default histos" << endl;
                return;
            }
            if (fHistoMergedValPtBinPlot){
                if (fHistoMergedValPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoMergedValPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}", startPt, endPt),
                                            "#it{#sigma}_{long}^{2}", Form("dN_{%s}/d#it{#sigma}_{long}^{2}", decayChannel.Data()),
                                            fPlottingM02Range[0],fPlottingM02Range[1],0,kBlue+1, 24);
                }
            }
            if (fHistoMergedPartConvValPtBinPlot){
                if (fHistoMergedPartConvValPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoMergedPartConvValPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}", startPt, endPt),
                                            "#it{#sigma}_{long}^{2}", Form("dN_{%s}/d#it{#sigma}_{long}^{2}", decayChannel.Data()),
                                            fPlottingM02Range[0],fPlottingM02Range[1],0,kGreen+2, 24);
                }
            }
            if (fHistoValGammaPtBinPlot){
                if (fHistoValGammaPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoValGammaPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}", startPt, endPt),
                                            "#it{#sigma}_{long}^{2}", Form("dN_{%s}/d#it{#sigma}_{long}^{2}", decayChannel.Data()),
                                            fPlottingM02Range[0],fPlottingM02Range[1],0,kRed+2, 25);
                }
            }
            if (fHistoValElectronPtBinPlot){
                if (fHistoValElectronPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoValElectronPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}", startPt, endPt),
                                            "#it{#sigma}_{long}^{2}", Form("dN_{%s}/d#it{#sigma}_{long}^{2}", decayChannel.Data()),
                                            fPlottingM02Range[0],fPlottingM02Range[1],0,807, 25);
                }
            }
            if (fHistoBGPtBinPlot){
                if (fHistoBGPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoBGPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}", startPt, endPt),
                                            "#it{#sigma}_{long}^{2}", Form("dN_{%s}/d#it{#sigma}_{long}^{2}", decayChannel.Data()),
                                            fPlottingM02Range[0],fPlottingM02Range[1],0,kGray+2);
                }
            }
        }
    }
    canvasDataSpectra->Print(namePlot.Data());
    delete padDataSpectra;
    delete canvasDataSpectra;
}

//__________________________________________ Plotting all Invariant Mass bins _______________________________________________
void PlotM02MergedTrueInPtBins(     TH1D**      fHistoMergedValPtBinPlot, 
                                    TH1D**      fHistoMergedPartConvValPtBinPlot, 
                                    TH1D**      fHistoPi0ValPtBinPlot,
                                    TH1D**      fHistoPi0PartConvValPtBinPlot,
                                    TH1D**      fHistoEtaValPtBinPlot,
                                    TH1D**      fHistoEtaPartConvValPtBinPlot,
                                    TString     namePlot, 
                                    TString     nameCanvas, 
                                    TString     namePad, 
                                    Double_t*   fPlottingM02Range, 
                                    TString     dateDummy, 
                                    TString     fMesonType, 
                                    Int_t       fRowPlot, 
                                    Int_t       fColumnPlot, 
                                    Int_t       fStartBinPtRange, 
                                    Int_t       fNumberPtBins, 
                                    Double_t*   fRangeBinsPt, 
                                    TString     fDecayChannel, 
                                    Bool_t      fMonteCarloInfo, 
                                    TString     decayChannel,  
                                    TString     fDetectionChannel, 
                                    TString     fEnergy
                              ){
    TGaxis::SetMaxDigits(3);
    TCanvas * canvasDataSpectra         = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    DrawGammaCanvasSettings(canvasDataSpectra, 0.0, 0., 0., 0.);

    TPad * padDataSpectra               = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    padDataSpectra->SetFillColor(0);
    padDataSpectra->GetFrame()->SetFillColor(0);
    padDataSpectra->SetBorderMode(0);
    padDataSpectra->SetLogy(1);
    padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
    padDataSpectra->Draw();

    cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
    cout << fMesonType.Data() << endl;
    Int_t place                         = 0;
    for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
        cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
        Double_t startPt                = fRangeBinsPt[iPt];
        Double_t endPt                  = fRangeBinsPt[iPt+1];

        place                           = place + 1;  //give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataSpectra->cd(place);

            cout << "entered ALICE plotting" << endl;
            TString textAlice           = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo){
                textEvents              = "MC";
            } else {
                textEvents              = "Data";
            }
            Double_t nPixels            = 13;
            Double_t textHeight         = 0.08;
            if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
                textHeight              = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
            } else {
                textHeight              = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
            }   
            
            Double_t startTextX         = 0.10;
            Double_t startTextY         = 0.9;
            Double_t differenceText     = textHeight*1.25;

            TLatex *alice               = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
            TLatex *latexDate           = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy              = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
            TLatex *process             = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess          = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
            TLatex *events              = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

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
            
            Int_t nLegend               = 0;
            if (fHistoMergedValPtBinPlot) nLegend++;
            if (fHistoMergedPartConvValPtBinPlot) nLegend++;
            if (fHistoPi0ValPtBinPlot) nLegend++;
            if (fHistoPi0PartConvValPtBinPlot) nLegend++;
            if (fHistoEtaValPtBinPlot) nLegend++;
            if (fHistoEtaPartConvValPtBinPlot) nLegend++;
            
            TLegend* legend = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+nLegend)*differenceText);
            legend->SetTextSize(textHeight);
            legend->SetTextFont(62);
            legend->SetFillColor(0);
            legend->SetFillStyle(0);
            legend->SetLineWidth(0);
            legend->SetLineColor(0);
            legend->SetMargin(0.15);
            Size_t markersize           = 0;
            Size_t linesize             = 0;
            if (fHistoMergedValPtBinPlot){
                if (fHistoMergedValPtBinPlot[iPt]){
                    markersize          = fHistoMergedValPtBinPlot[iPt]->GetMarkerSize();
                    fHistoMergedValPtBinPlot[iPt]->SetMarkerSize(1.*markersize);
                    legend->AddEntry(fHistoMergedValPtBinPlot[iPt],"merged clus.","ep");
                }
            }
            if (fHistoMergedPartConvValPtBinPlot){
                if (fHistoMergedPartConvValPtBinPlot[iPt]){
                    markersize          = fHistoMergedPartConvValPtBinPlot[iPt]->GetMarkerSize();
                    fHistoMergedPartConvValPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoMergedPartConvValPtBinPlot[iPt],"merged clus. p. conv","ep");
                }
            }
            if (fHistoPi0ValPtBinPlot){
                if (fHistoPi0ValPtBinPlot[iPt]){
                    markersize          = fHistoPi0ValPtBinPlot[iPt]->GetMarkerSize();
                    fHistoPi0ValPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoPi0ValPtBinPlot[iPt],"merged #pi^{0}","ep");
                }
            }
            if (fHistoPi0PartConvValPtBinPlot){
                if (fHistoPi0PartConvValPtBinPlot[iPt]){
                    markersize          = fHistoPi0PartConvValPtBinPlot[iPt]->GetMarkerSize();
                    fHistoPi0PartConvValPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoPi0PartConvValPtBinPlot[iPt],"merged #pi^{0} p. conv","ep");
                }
            }
            if (fHistoEtaValPtBinPlot){
                if (fHistoEtaValPtBinPlot[iPt]){
                    markersize          = fHistoEtaValPtBinPlot[iPt]->GetMarkerSize();
                    fHistoEtaValPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoEtaValPtBinPlot[iPt],"merged #eta","ep");
                }
            }
            if (fHistoEtaPartConvValPtBinPlot){
                if (fHistoEtaPartConvValPtBinPlot[iPt]){
                    markersize          = fHistoEtaPartConvValPtBinPlot[iPt]->GetMarkerSize();
                    fHistoEtaPartConvValPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoEtaPartConvValPtBinPlot[iPt],"merged #eta p. conv","ep");
                }
            }
            legend->Draw(); 
        } else {            
            padDataSpectra->cd(place);
            padDataSpectra->cd(place)->SetTopMargin(0.12);
            padDataSpectra->cd(place)->SetBottomMargin(0.15);
            padDataSpectra->cd(place)->SetRightMargin(0.02);
//             cout << "place: " << place << "\t "<< startPt << " < #it{p}_{T} < " << endPt  << endl;
            int remaining               = (place-1)%fColumnPlot;
            if (remaining > 0) padDataSpectra->cd(place)->SetLeftMargin(0.15);
            else padDataSpectra->cd(place)->SetLeftMargin(0.25);
//             cout << "here" << endl;
            if (fHistoMergedValPtBinPlot){
                if (fHistoMergedValPtBinPlot[iPt]){
                    DrawGammaHisto( fHistoMergedValPtBinPlot[iPt],
                                    Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}", startPt, endPt),
                                    "#it{#sigma}_{long}^{2}", Form("dN_{%s}/d#it{#sigma}_{long}^{2}", decayChannel.Data()),
                                    fPlottingM02Range[0],fPlottingM02Range[1],0,0.8);
                } else {
                    continue;
                }
            } else {
                cout << "ABORT: plotting no default histos" << endl;
                return;
            }
            if (fHistoMergedPartConvValPtBinPlot){
                if (fHistoMergedPartConvValPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoMergedPartConvValPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}", startPt, endPt),
                                            "#it{#sigma}_{long}^{2}", Form("dN_{%s}/d#it{#sigma}_{long}^{2}", decayChannel.Data()),
                                            fPlottingM02Range[0],fPlottingM02Range[1],0,kGray+2, 20);
                }
            }
            if (fHistoPi0ValPtBinPlot){
                if (fHistoPi0ValPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoPi0ValPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}", startPt, endPt),
                                            "#it{#sigma}_{long}^{2}", Form("dN_{%s}/d#it{#sigma}_{long}^{2}", decayChannel.Data()),
                                            fPlottingM02Range[0],fPlottingM02Range[1],0,kBlue+1, 24);
                }
            }
            if (fHistoPi0PartConvValPtBinPlot){
                if (fHistoPi0PartConvValPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoPi0PartConvValPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}", startPt, endPt),
                                            "#it{#sigma}_{long}^{2}", Form("dN_{%s}/d#it{#sigma}_{long}^{2}", decayChannel.Data()),
                                            fPlottingM02Range[0],fPlottingM02Range[1],0,kAzure-4, 24);
                }
            }
            if (fHistoEtaValPtBinPlot){
                if (fHistoEtaValPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoEtaValPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}", startPt, endPt),
                                            "#it{#sigma}_{long}^{2}", Form("dN_{%s}/d#it{#sigma}_{long}^{2}", decayChannel.Data()),
                                            fPlottingM02Range[0],fPlottingM02Range[1],0,kGreen+2, 25);
                }
            }
            if (fHistoEtaPartConvValPtBinPlot){
                if (fHistoEtaPartConvValPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoEtaPartConvValPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}", startPt, endPt),
                                            "#it{#sigma}_{long}^{2}", Form("dN_{%s}/d#it{#sigma}_{long}^{2}", decayChannel.Data()),
                                            fPlottingM02Range[0],fPlottingM02Range[1],0,kSpring+5,25);
                }
            }
        }
    }
    canvasDataSpectra->Print(namePlot.Data());
    delete padDataSpectra;
    delete canvasDataSpectra;
}

//__________________________________________ Plotting all Invariant Mass bins _______________________________________________
void PlotM02MergedTruePrimSecInPtBins(      TH1D**      fHistoPi0ValPtBinPlot, 
                                            TH1D**      fHistoPi0PrimValPtBinPlot,
                                            TH1D**      fHistoPi0SecValPtBinPlot,
                                            TH1D**      fHistoPi0SecFK0sValPtBinPlot,
                                            TH1D**      fHistoPi0SecFLambdaValPtBinPlot,
                                            TString     namePlot, 
                                            TString     nameCanvas, 
                                            TString     namePad, 
                                            Double_t*   fPlottingM02Range, 
                                            TString     dateDummy, 
                                            TString     fMesonType, 
                                            Int_t       fRowPlot, 
                                            Int_t       fColumnPlot, 
                                            Int_t       fStartBinPtRange, 
                                            Int_t       fNumberPtBins, 
                                            Double_t*   fRangeBinsPt, 
                                            TString     fDecayChannel, 
                                            Bool_t      fMonteCarloInfo, 
                                            TString     decayChannel,  
                                            TString     fDetectionChannel, 
                                            TString     fEnergy
                                     ){
    TGaxis::SetMaxDigits(3);
    TCanvas * canvasDataSpectra         = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    DrawGammaCanvasSettings(canvasDataSpectra, 0.0, 0., 0., 0.);

    TPad * padDataSpectra               = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    padDataSpectra->SetFillColor(0);
    padDataSpectra->GetFrame()->SetFillColor(0);
    padDataSpectra->SetBorderMode(0);
    padDataSpectra->SetLogy(1);
    padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
    padDataSpectra->Draw();

    cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
    cout << fMesonType.Data() << endl;
    Int_t place                         = 0;
    for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
        cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
        Double_t startPt                = fRangeBinsPt[iPt];
        Double_t endPt                  = fRangeBinsPt[iPt+1];

        place                           = place + 1;  //give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataSpectra->cd(place);

            cout << "entered ALICE plotting" << endl;
            TString textAlice           = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo){
                textEvents              = "MC";
            } else {
                textEvents              = "Data";
            }
            Double_t nPixels            = 13;
            Double_t textHeight         = 0.08;
            if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
                textHeight              = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
            } else {
                textHeight              = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
            }   
            
            Double_t startTextX         = 0.10;
            Double_t startTextY         = 0.9;
            Double_t differenceText     = textHeight*1.25;

            TLatex *alice               = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
            TLatex *latexDate           = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy              = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
            TLatex *process             = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess          = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
            TLatex *events              = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

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
            
            Int_t nLegend               = 0;
            if (fHistoPi0ValPtBinPlot) nLegend++;
            if (fHistoPi0PrimValPtBinPlot) nLegend++;
            if (fHistoPi0SecValPtBinPlot) nLegend++;
            if (fHistoPi0SecFK0sValPtBinPlot) nLegend++;
            if (fHistoPi0SecFLambdaValPtBinPlot) nLegend++;
            
            TLegend* legend             = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+nLegend)*differenceText);
            legend->SetTextSize(textHeight);
            legend->SetTextFont(62);
            legend->SetFillColor(0);
            legend->SetFillStyle(0);
            legend->SetLineWidth(0);
            legend->SetLineColor(0);
            legend->SetMargin(0.15);
            Size_t markersize           = 0;
            Size_t linesize             = 0;
            if (fHistoPi0ValPtBinPlot){
                if (fHistoPi0ValPtBinPlot[iPt]){
                    markersize          = fHistoPi0ValPtBinPlot[iPt]->GetMarkerSize();
                    fHistoPi0ValPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoPi0ValPtBinPlot[iPt],"merged #pi^{0}","ep");
                }
            }
            if (fHistoPi0PrimValPtBinPlot){
                if (fHistoPi0PrimValPtBinPlot[iPt]){
                    markersize          = fHistoPi0PrimValPtBinPlot[iPt]->GetMarkerSize();
                    fHistoPi0PrimValPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoPi0PrimValPtBinPlot[iPt],"merged prim. #pi^{0}","ep");
                }
            }
            if (fHistoPi0SecValPtBinPlot){
                if (fHistoPi0SecValPtBinPlot[iPt]){
                    markersize          = fHistoPi0SecValPtBinPlot[iPt]->GetMarkerSize();
                    fHistoPi0SecValPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoPi0SecValPtBinPlot[iPt],"merged sec #pi^{0}","ep");
                }
            }
            if (fHistoPi0SecFK0sValPtBinPlot){
                if (fHistoPi0SecFK0sValPtBinPlot[iPt]){
                    markersize          = fHistoPi0SecFK0sValPtBinPlot[iPt]->GetMarkerSize();
                    fHistoPi0SecFK0sValPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoPi0SecFK0sValPtBinPlot[iPt],"merged sec #pi^{0} f. K^{0}_{s}","ep");
                }
            }
            if (fHistoPi0SecFLambdaValPtBinPlot){
                if (fHistoPi0SecFLambdaValPtBinPlot[iPt]){
                    markersize          = fHistoPi0SecFLambdaValPtBinPlot[iPt]->GetMarkerSize();
                    fHistoPi0SecFLambdaValPtBinPlot[iPt]->SetMarkerSize(3*markersize);
                    legend->AddEntry(fHistoPi0SecFLambdaValPtBinPlot[iPt],"merged sec #pi^{0} f. #Lambda","ep");
                }
            }
            legend->Draw(); 
        } else {            
            padDataSpectra->cd(place);
            padDataSpectra->cd(place)->SetTopMargin(0.12);
            padDataSpectra->cd(place)->SetBottomMargin(0.15);
            padDataSpectra->cd(place)->SetRightMargin(0.02);
//             cout << "place: " << place << "\t "<< startPt << " < #it{p}_{T} < " << endPt  << endl;
            int remaining               = (place-1)%fColumnPlot;
            if (remaining > 0) padDataSpectra->cd(place)->SetLeftMargin(0.15);
            else padDataSpectra->cd(place)->SetLeftMargin(0.25);
//             cout << "here" << endl;
            if (fHistoPi0ValPtBinPlot){
                if (fHistoPi0ValPtBinPlot[iPt]){
                    DrawGammaHisto( fHistoPi0ValPtBinPlot[iPt],
                                    Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}", startPt, endPt),
                                    "#it{#sigma}_{long}^{2}", Form("dN_{%s}/d#it{#sigma}_{long}^{2}", decayChannel.Data()),
                                    fPlottingM02Range[0],fPlottingM02Range[1],0);
                } else {
                    continue;
                }
            } else {
                cout << "ABORT: plotting no default histos" << endl;
                return;
            }
            if (fHistoPi0PrimValPtBinPlot){
                if (fHistoPi0PrimValPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoPi0PrimValPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}", startPt, endPt),
                                            "#it{#sigma}_{long}^{2}", Form("dN_{%s}/d#it{#sigma}_{long}^{2}", decayChannel.Data()),
                                            fPlottingM02Range[0],fPlottingM02Range[1],0,kRed+1, 24);
                }
            }
            if (fHistoPi0SecValPtBinPlot){
                if (fHistoPi0SecValPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoPi0SecValPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}", startPt, endPt),
                                            "#it{#sigma}_{long}^{2}", Form("dN_{%s}/d#it{#sigma}_{long}^{2}", decayChannel.Data()),
                                            fPlottingM02Range[0],fPlottingM02Range[1],0,kBlue+1, 24);
                }
            }
            if (fHistoPi0SecFK0sValPtBinPlot){
                if (fHistoPi0SecFK0sValPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoPi0SecFK0sValPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}", startPt, endPt),
                                            "#it{#sigma}_{long}^{2}", Form("dN_{%s}/d#it{#sigma}_{long}^{2}", decayChannel.Data()),
                                            fPlottingM02Range[0],fPlottingM02Range[1],0,kAzure-4, 25);
                }
            }
            if (fHistoPi0SecFLambdaValPtBinPlot){
                if (fHistoPi0SecFLambdaValPtBinPlot[iPt]){
                    DrawGammaHistoColored(  fHistoPi0SecFLambdaValPtBinPlot[iPt],
                                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}", startPt, endPt),
                                            "#it{#sigma}_{long}^{2}", Form("dN_{%s}/d#it{#sigma}_{long}^{2}", decayChannel.Data()),
                                            fPlottingM02Range[0],fPlottingM02Range[1],0,kGreen+2, 26);
                }
            }
        }
    }
    canvasDataSpectra->Print(namePlot.Data());
    delete padDataSpectra;
    delete canvasDataSpectra;
}

