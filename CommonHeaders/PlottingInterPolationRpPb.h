

void SetStyleTGraphErrorsForGraphs(TGraphErrors* graph, TString XTitle, TString YTitle, Size_t xLableSize, Size_t xTitleSize, Size_t yLableSize, Size_t yTitleSize, Float_t xTitleOffset = 1, Float_t yTitleOffset = 1, Int_t xNDivisions = 510, Int_t yNDivisions = 510){
	graph->GetXaxis()->SetTitle(XTitle);
	graph->GetYaxis()->SetTitle(YTitle);
	graph->SetTitle("");
	
	graph->GetXaxis()->SetLabelSize(xLableSize);
	graph->GetXaxis()->SetTitleSize(xTitleSize);
	graph->GetXaxis()->SetTitleOffset(xTitleOffset);
	graph->GetXaxis()->SetNdivisions(xNDivisions,kTRUE);
	
	graph->GetXaxis()->SetLabelFont(42);
	graph->GetYaxis()->SetLabelFont(42); 
	graph->GetXaxis()->SetTitleFont(62);
	graph->GetYaxis()->SetTitleFont(62);

	
	graph->GetYaxis()->SetDecimals();
	graph->GetYaxis()->SetLabelSize(yLableSize);
	graph->GetYaxis()->SetTitleSize(yTitleSize);
	graph->GetYaxis()->SetTitleOffset(yTitleOffset);
	graph->GetYaxis()->SetNdivisions(yNDivisions,kTRUE);	
}

void SetStyleTGraphAsymmErrorsForGraphs(TGraphAsymmErrors* graph, TString XTitle, TString YTitle, Size_t xLableSize, Size_t xTitleSize, Size_t yLableSize, Size_t yTitleSize, Float_t xTitleOffset = 1, Float_t yTitleOffset = 1, Int_t xNDivisions = 510, Int_t yNDivisions = 510){
	graph->GetXaxis()->SetTitle(XTitle);
	graph->GetYaxis()->SetTitle(YTitle);
	graph->SetTitle("");
	
	graph->GetXaxis()->SetLabelSize(xLableSize);
	graph->GetXaxis()->SetTitleSize(xTitleSize);
	graph->GetXaxis()->SetTitleOffset(xTitleOffset);
	graph->GetXaxis()->SetNdivisions(xNDivisions,kTRUE);
	
	graph->GetXaxis()->SetLabelFont(42);
	graph->GetYaxis()->SetLabelFont(42); 
	graph->GetXaxis()->SetTitleFont(62);
	graph->GetYaxis()->SetTitleFont(62);

	
	graph->GetYaxis()->SetDecimals();
	graph->GetYaxis()->SetLabelSize(yLableSize);
	graph->GetYaxis()->SetTitleSize(yTitleSize);
	graph->GetYaxis()->SetTitleOffset(yTitleOffset);
	graph->GetYaxis()->SetNdivisions(yNDivisions,kTRUE);	
}




void PlotInterpolationPtBins(TGraphErrors** gPtvSqrts,TGraphErrors** gPtvsEnergies, TF1** fPowerlaw, TGraphErrors* gRpPb,Int_t fColumnPlot, Int_t fRowPlot,TString namePlot){

    TGaxis::SetMaxDigits(3);
    
    TString nameCanvas = "";
    TString namePad    = "";
    //TString namePlot   = "";
    
    TCanvas * canvasPtvsSqrts     = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    canvasPtvsSqrts->SetTopMargin(0.00);
    canvasPtvsSqrts->SetBottomMargin(0.00);
    canvasPtvsSqrts->SetRightMargin(0.0);
    canvasPtvsSqrts->SetLeftMargin(0.00);

    TPad * padPtvsSqrts           = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
       
    padPtvsSqrts->SetFillColor(0);
    padPtvsSqrts->GetFrame()->SetFillColor(0);
    padPtvsSqrts->SetBorderMode(0);
    padPtvsSqrts->Divide(fColumnPlot,fRowPlot,0.0,0.0);
    padPtvsSqrts->SetLeftMargin(0.);
    padPtvsSqrts->SetRightMargin(0.);
    padPtvsSqrts->SetTopMargin(0.);
    padPtvsSqrts->SetBottomMargin(0.);
    padPtvsSqrts->Draw();
    
    
    

    cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
    
    Int_t place = 0;
    
    for(Int_t iPt = 0; iPt <gRpPb->GetN(); iPt++){
      
      cout<<"fallo"<<endl;
      
      place++;
      
      
      padPtvsSqrts->cd(place);
      padPtvsSqrts->cd(place)->SetTopMargin(0.12);
      padPtvsSqrts->cd(place)->SetBottomMargin(0.15);
      padPtvsSqrts->cd(place)->SetRightMargin(0.035);
      //padPtvsSqrts->cd(place)->SetLeftMargin(0.25);
      
      int remaining           = (place-1)%fColumnPlot;
      if (remaining > 0) padPtvsSqrts->cd(place)->SetLeftMargin(0.15);
      else padPtvsSqrts->cd(place)->SetLeftMargin(0.25);

      
      
      DrawGammaSetMarkerTGraphErr(gPtvSqrts[iPt],21,1.5, kRed , kRed);
      DrawGammaSetMarkerTGraphErr(gPtvsEnergies[iPt],20,1.5, kBlack , kBlack);
      
      gPtvsEnergies[iPt]->GetXaxis()->SetTitle("#sqrt{s}");
      gPtvsEnergies[iPt]->GetXaxis()->SetTitleSize(0.065);
      gPtvsEnergies[iPt]->GetXaxis()->SetLabelSize(0.06);
      gPtvsEnergies[iPt]->GetYaxis()->SetTitle("invariant cross section");      
      gPtvsEnergies[iPt]->GetYaxis()->SetTitleSize(0.065);
      gPtvsEnergies[iPt]->GetYaxis()->SetLabelSize(0.06);
      gPtvsEnergies[iPt]->GetXaxis()->SetNdivisions(308,kTRUE);
      gPtvsEnergies[iPt]->GetYaxis()->SetNdivisions(304,kTRUE);
      gPtvsEnergies[iPt]->GetXaxis()->SetLabelOffset(0.015);
      gPtvsEnergies[iPt]->GetYaxis()->SetLabelOffset(0.01);
      
      //gPtvsEnergies[iPt]->Set
      
      
      //gPtvSqrts[iPt]->SetLineColor(kRed);
      //gPtvSqrts[iPt]->Draw("E1psame");
	//graphAsymmErrorsPi0DSS5000->Draw("same,E3");
      gPtvsEnergies[iPt]->Draw("ap");
      //gPtvSqrts[iPt]->Draw("same");
      gPtvSqrts[iPt]->Draw("same, p");
      fPowerlaw[iPt]->SetLineColor(kBlue);
      fPowerlaw[iPt]->SetLineWidth(2);
      fPowerlaw[iPt]->Draw("same");
      
      TString Title = Form("#it{p}_{T} = %3.2f GeV/#it{c} ",gRpPb->GetX()[iPt]);
      
      if(Title.Length() > 0){
		gPtvsEnergies[iPt]->SetTitle("");
		
		Double_t xMin,yMin;
		yMin = 0.78;
		if ( remaining > 0 ) xMin = 0.20;
		else xMin = 0.30;
		
		TLatex *alice = new TLatex(xMin,yMin,Form("%s",Title.Data())); // Bo: this was 
		alice->SetNDC();
		alice->SetTextColor(1);
		alice->SetTextSize(0.062);
		alice->Draw();		
      }
      
       
      
    }

    canvasPtvsSqrts->Print(namePlot.Data());
    delete padPtvsSqrts;
    delete canvasPtvsSqrts;
    
}

//void SetStyleTGraphErrorsForGraphs(TGraphErrors* graph, TString XTitle, TString YTitle, Size_t xLableSize, Size_t xTitleSize, Size_t yLableSize, Size_t yTitleSize, Float_t xTitleOffset = 1, Float_t yTitleOffset = 1, Int_t xNDivisions = 510, Int_t yNDivisions = 510){

void PlotAlphavsPt(TGraphErrors* gAlpha, TString method, TString thesisPlotLabel, TString namePlot){
  
    TCanvas * canvasAlphavsPt     = new TCanvas("AlphavsPt","",1400,900);  // gives the page size
    DrawGammaCanvasSettings( canvasAlphavsPt,  0.1, 0.02, 0.03, 0.1);
	
    SetStyleTGraphErrorsForGraphs(gAlpha,pTLabel.Data(),"#alpha", 0.04,0.04, 0.04,0.04, 1.,1.2, 512, 512); 

    
    
    DrawGammaSetMarkerTGraphErr(gAlpha,21,1.5, kRed , kRed);
    
    
    
    gAlpha->GetYaxis()->CenterTitle();
   
    
    gAlpha->Draw("ap");
    
    cout<<"Aqui debe decir Te "<<thesisPlotLabel.Data()<<endl;
    
    TLatex *labelThesis = new TLatex(0.15,0.90,thesisPlotLabel.Data());
    SetStyleTLatex( labelThesis, 0.04,4);
    labelThesis->Draw();
	
    TLatex *labelEnergy  =  new TLatex(0.15,0.85, energyLabel.Data());
    SetStyleTLatex( labelEnergy, 0.04,4);
    labelEnergy->Draw();
    
    TString methodLabel = "";
    
    if( method.CompareTo("PCM") == 0 ) methodLabel = PCMLabel;
    else if ( method.CompareTo("Dalitz") == 0) methodLabel = DalitzLabel;
     else  methodLabel = CaloLabel;
    TLatex *labelMethod = new TLatex(0.15,0.80,methodLabel.Data());
    SetStyleTLatex( labelMethod, 0.03,4);
    labelMethod->Draw();

  
    canvasAlphavsPt->Print(namePlot.Data());
    
    delete canvasAlphavsPt;
  
}

 
void PlotErrors(TGraphErrors* interpolationYieldStatErr,  TGraphErrors* interpolationYieldSystErr, 
		TGraphErrors* InvYieldPi07TeVStatErr,     TGraphErrors* InvYieldPi07TeVSystErr,
		TGraphErrors* InvYieldPi02760GeVStatErr,  TGraphErrors* InvYieldPi02760GeVSystErr,
		TString method, TString thesisPlotLabel,  TString outputDir, TString suffix){
  
  
    TCanvas * canvasYieldStatErr     = new TCanvas("canvasYieldStatErr","",1400,900);  // gives the page size
    DrawGammaCanvasSettings( canvasYieldStatErr,  0.1, 0.02, 0.03, 0.1);
    
    canvasYieldStatErr->SetLogx();
    
    TH2F * dummyStatErr;
	dummyStatErr = new TH2F("dummyStatErr","dummyStatErr",1000,0.23,30.,1000,0,50);
	SetStyleHistoTH2ForGraphs(dummyStatErr, "#it{p}_{T} (GeV/#it{c})", "Statistical Err (%)", 0.032,0.04, 0.04,0.04, 1,1.2);
       dummyStatErr->DrawCopy(); 

    
    
    DrawGammaSetMarkerTGraphErr(interpolationYieldStatErr,20,2.5,  kRed+1 , 	kRed+1);
    DrawGammaSetMarkerTGraphErr(InvYieldPi02760GeVStatErr,22,2.5,  kGreen+2,    kGreen+2);
    DrawGammaSetMarkerTGraphErr(InvYieldPi07TeVStatErr,   34,2.5,  kBlue+1,     kBlue+1);
    InvYieldPi07TeVStatErr->SetFillColor(0);
    InvYieldPi02760GeVStatErr->SetFillColor(0);
    interpolationYieldStatErr->SetFillColor(0);
    
    
    
    interpolationYieldStatErr->GetYaxis()->CenterTitle();
   
    
    interpolationYieldStatErr->Draw("p,same");
    InvYieldPi02760GeVStatErr->Draw("p,same");
    InvYieldPi07TeVStatErr->Draw("p, same");
    
   
    
    TLatex *labelThesis = new TLatex(0.15,0.90,thesisPlotLabel.Data());
    SetStyleTLatex( labelThesis, 0.04,4);
    labelThesis->Draw();
	
    TLatex *labelEnergy  =  new TLatex(0.15,0.85, energyLabel.Data());
    SetStyleTLatex( labelEnergy, 0.04,4);
    labelEnergy->Draw();
    
    TString methodLabel = "";
    
    if( method.CompareTo("PCM") == 0 ) methodLabel = PCMLabel;
    else if ( method.CompareTo("Dalitz") == 0) methodLabel = DalitzLabel;
    
    TLatex *labelMethod = new TLatex(0.15,0.80,methodLabel.Data());
    SetStyleTLatex( labelMethod, 0.03,4);
    labelMethod->Draw();
    
    
    TLegend* legendInvYields = new TLegend(0.5,0.75,0.9,0.95);
	legendInvYields->SetFillColor(0);
	legendInvYields->SetLineColor(0);
	legendInvYields->SetNColumns(1);
	legendInvYields->SetTextSize(0.04);
	legendInvYields->AddEntry(InvYieldPi07TeVStatErr,"calculated pp #sqrt{#it{s}} = 7 TeV","p");
	legendInvYields->AddEntry(interpolationYieldStatErr,"pp reference #sqrt{#it{s}} = 5.02 TeV","p");
	legendInvYields->AddEntry(InvYieldPi02760GeVStatErr,"calculated pp #sqrt{#it{s}} = 2.76 TeV","p");
	
	legendInvYields->Draw();
    

    TString namePlotStatErr = Form("%s/StatErrors_vs_Pt_%s.%s",outputDir.Data(),method.Data(),suffix.Data());

    canvasYieldStatErr->Print(namePlotStatErr.Data());
    
    
    
    TCanvas * canvasYieldSystErr     = new TCanvas("canvasYieldSystErr","",1400,900);  // gives the page size
    DrawGammaCanvasSettings( canvasYieldSystErr,  0.1, 0.02, 0.03, 0.1);
    
    canvasYieldSystErr->SetLogx();
    
    TH2F * dummySystErr;
    dummySystErr = new TH2F("dummySystErr","dummySystErr",1000,0.23,30.,1000,0,30);
    SetStyleHistoTH2ForGraphs(dummySystErr, "#it{p}_{T} (GeV/#it{c})", "Systmeatical Err (%)", 0.032,0.04, 0.04,0.04, 1,1.2);
    dummySystErr->DrawCopy(); 

    
    
    DrawGammaSetMarkerTGraphErr(interpolationYieldSystErr,20,2.5,   kRed+1,    kRed+1);
    DrawGammaSetMarkerTGraphErr(InvYieldPi02760GeVSystErr,22,2.5,   kGreen+2,  kGreen+2);
    DrawGammaSetMarkerTGraphErr(InvYieldPi07TeVSystErr,   34,2.5,   kBlue+1,   kBlue+1);
    
    InvYieldPi07TeVSystErr->SetFillColor(0);
    InvYieldPi02760GeVSystErr->SetFillColor(0);
    interpolationYieldSystErr->SetFillColor(0);
	
    
    
    interpolationYieldSystErr->GetYaxis()->CenterTitle();
   
    
    interpolationYieldSystErr->Draw("p,same");
    InvYieldPi02760GeVSystErr->Draw("p,same");
    InvYieldPi07TeVSystErr->Draw("p, same");
    
    
       
    labelThesis->Draw();
    labelEnergy->Draw();
     
    labelMethod->Draw();

    TString namePlotSystErr = Form("%s/SystErrors_vs_Pt_%s.%s",outputDir.Data(),method.Data(),suffix.Data());

    legendInvYields->Draw();
  
    canvasYieldSystErr->Print(namePlotSystErr.Data());
    
    
    
    delete canvasYieldStatErr;
  
}
