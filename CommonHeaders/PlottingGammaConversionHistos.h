/***********************************************************************************************
*** provided by Gamma Conversion Group, PWG4, 									******
***	    Friederike Bock, fbock@physi.uni-heidelberg.de ***							******
************************************************************************************************/


/************************************************************************************************
************************************************************************************************
* This header was created to make things easier with making plots for the gamma conversion group.
it offers you several functions for drawing and styling your histogramms.
************************************************************************************************
************************************************************************************************

The functions are 
- StyleSettingsThesis
- StyleSettings
- GammaScalingHistogramm

- DrawAutoGammaHistos
- DrawAutoGammaHisto
- DrawAutoGammaHisto2D

- DrawRatioGammaHisto

- DrawCutGammaHisto
- DrawCutGammaHistos

- DrawGammaLines
*/

#include <Riostream.h>

// extern TRandom *kgRandom;
// extern TBenchmark *kgBenchmark;
// extern TSystem *kgSystem;


// ---------------------------- Function definiton --------------------------------------------------------------------------------------------


/* StyleSettingsThesis will make some standard settings for gStyle 
*/
void StyleSettingsThesis( TString format = ""){
    //gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);   //show day and time
    gStyle->SetOptStat(0);  //show statistic
    gStyle->SetPalette(1,0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetTitleFillColor(0);
    gStyle->SetTextSize(0.5);
    gStyle->SetLabelSize(0.03,"xyz");
    gStyle->SetLabelOffset(0.002,"xyz");
    gStyle->SetTitleFontSize(0.04);
    gStyle->SetTitleOffset(1,"y");
    gStyle->SetTitleOffset(0.7,"x");		
    gStyle->SetCanvasColor(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLineWidth(1);
    
    gStyle->SetPadTopMargin(0.03);
    gStyle->SetPadBottomMargin(0.09);
    gStyle->SetPadRightMargin(0.03);
    gStyle->SetPadLeftMargin(0.13);
    
    
    TGaxis::SetMaxDigits(3);
    gErrorIgnoreLevel=kError;
    
    if (format.CompareTo("eps") == 0 ||format.CompareTo("pdf") == 0  ) gStyle->SetLineScalePS(1);
}


/* StyleSettings will make some standard settings for gStyle 
*/
void StyleSettings(){
    //gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);   //show day and time
    gStyle->SetOptStat(0);  //show statistic
    gStyle->SetPalette(1,0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetTitleFillColor(0);
    gStyle->SetTextSize(0.5);
    gStyle->SetLabelSize(0.03,"xyz");
    gStyle->SetLabelOffset(0.002,"xyz");
    gStyle->SetTitleFontSize(0.04);
    gStyle->SetTitleOffset(1,"y");
    gStyle->SetTitleOffset(0.7,"x");		
    gStyle->SetCanvasColor(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLineWidth(1);
    
    gStyle->SetPadTopMargin(0.1);
    gStyle->SetPadBottomMargin(0.1);
    gStyle->SetPadRightMargin(0.08);
    gStyle->SetPadLeftMargin(0.12);
    
    gErrorIgnoreLevel=kError;
    
    TGaxis::SetMaxDigits(3);
}


void SetPlotStyle() {
// 	const Int_t nRGBs = 7;
    const Int_t nRGBs = 5;
    const Int_t nCont = 255;
    
    Double_t stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[nRGBs] = { 0.31, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[nRGBs]  = { 0.51, 1., 0.12, 0.00, 0.00};

// 	Double_t stops[nRGBs] = {  0.34, 0.61, 0.84, 1.00 };
// 	Double_t red[nRGBs]   = {  0.00, 0.87, 1.00, 0.51 };
// 	Double_t green[nRGBs] = {  0.81, 1.00, 0.20, 0.00 };
// // 	Double_t blue[nRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
// 	Double_t blue[nRGBs]  = {  1., 0.00, 0.00, 0.00};

// 	Double_t blue[nRGBs]  = { 1.00, 0.97, 0.97, 0.00, 0.00, 0.00, 0.00 };
    TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);
    gStyle->SetNumberContours(nCont);
}

void SetPlotStyleNConts(    Int_t nCont = 255) {
    const Int_t nRGBs = 5;
    Double_t stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[nRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[nRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);
    gStyle->SetNumberContours(nCont);
}


void DrawGammaSetMarker(    TH1* histo1, 
                            Style_t markerStyle, 
                            Size_t markerSize, 
                            Color_t markerColor, 
                            Color_t lineColor ) {
    histo1->SetMarkerStyle(markerStyle);
    histo1->SetMarkerSize(markerSize);
    histo1->SetMarkerColor(markerColor);
    histo1->SetLineColor(lineColor);	
    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);
}


// GammaScalingHistogram will scale the histogram by "Factor" 
void GammaScalingHistogramm(TH1 *histo, Double_t Factor){
    histo->Sumw2();
    histo->Scale(Factor);
} 

// GammaScalingHistogram will scale the histogram by "Factor" 
void GammaScalingHistogramm(TH2 *histo, Double_t Factor){
    histo->Sumw2();
    histo->Scale(Factor);
} 

void StylingSliceHistos(TH1 *histo, Float_t markersize){
    histo->SetMarkerStyle(20);
    histo->SetMarkerSize(markersize);
}

void ConvGammaRebinWithBinCorrection(TH1 *histo, Int_t rebinFactor, Int_t bin = 3){
    histo->Sumw2();
    histo->Rebin(rebinFactor);
    Double_t binWidth= histo->GetXaxis()->GetBinWidth(bin);
    for (Int_t i = 1; i < histo->GetNbinsX()+1; i++){
        histo->SetBinContent(i,histo->GetBinContent(i)/binWidth);	
        histo->SetBinError(i,histo->GetBinError(i)/binWidth);	
    }
}

void ConvGammaRebinWithBinCorrection2D(TH2 *histo, Int_t rebinFactor1, Int_t rebinFactor2, Int_t bin = 3){
// 	histo->Sumw2();
    histo->Rebin2D(rebinFactor1,rebinFactor2);
    Double_t binWidthY= histo->GetYaxis()->GetBinWidth(bin);
    Double_t binWidthX= histo->GetXaxis()->GetBinWidth(bin);
    histo->Scale(1/binWidthY*1/binWidthX);
}

void ConvGammaRebinWithBinCorrection2DSumw2(TH2 *histo, Int_t rebinFactor1, Int_t rebinFactor2, Int_t bin = 3){
    histo->Sumw2();
    histo->Rebin2D(rebinFactor1,rebinFactor2);
    Double_t binWidthY= histo->GetYaxis()->GetBinWidth(bin);
    Double_t binWidthX= histo->GetXaxis()->GetBinWidth(bin);
    histo->Scale(1/binWidthY*1/binWidthX);
}

/* DrawAutoGammaHistos is function used for styling the histograms of the gamma conversion group for two histos and standart settings
* histo1 - first histogram (Data)
* histo2 - second histogram (MC)
* Title - histogram title
* XTitle - X-axis title
* YTitle - Y-axis title
* YRangeMax 	= kTRUE will scale by Maximum and Minimum Range in Y
*YMaxFactor - will MaximumY by this factor if YRangeMay = kTRUE 
*YMinimum - this will be used if YRangeMax is set
*YRange  	= kTRUE will Cut y-axis by YMin and YMax 
- will be set to kFAlSE if YRangeMax is set
*YMin - minimum Y
*YMax - maximum Y
*XRange 	= kTRUE will Cut x-axis by XMin and XMax
*XMin - minimum Y
*XMax - maximum Y
*/ 

void DrawAutoGammaHistos(   TH1* histo1, 
                            TH1*histo2, 
                            TString Title, 
                            TString XTitle, 
                            TString YTitle,
                            Bool_t YRangeMax, 
                            Double_t YMaxFactor, 
                            Double_t YMinimum,
                            Bool_t YRange, 
                            Double_t YMin, 
                            Double_t YMax,
                            Bool_t XRange, 
                            Double_t XMin, 
                            Double_t XMax,
                            Float_t xOffset=1., 
                            Float_t yOffset=1.7) {

    if (YRangeMax && !XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        if(maxRangeR < histo2->GetMaximum()){
            maxRangeR = histo2->GetMaximum();
        }
        Double_t minRangeR = histo1->GetMinimum();		
        if(minRangeR > histo2->GetMinimum()){
            minRangeR = histo2->GetMinimum();
        }
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
    }
    if (YRangeMax && XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        if(maxRangeR < histo2->GetMaximum()){
            maxRangeR = histo2->GetMaximum();
        }
        Double_t minRangeR = histo1->GetMinimum();		
        if(minRangeR > histo2->GetMinimum()){
            minRangeR = histo2->GetMinimum();
        }
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (YRange && XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);	
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (!YRangeMax && !YRange && XRange){
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (YRange && !XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);
    }
    
    if(Title.CompareTo("") != 0){
        histo1->SetTitle(Title.Data());
    }else{histo1->SetTitle("");
    histo2->SetTitle("");}
    if(XTitle.CompareTo("") != 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.CompareTo("") != 0){
        histo1->SetYTitle(YTitle.Data());
    }
    DrawGammaSetMarker(histo1, 20, 0.5, kBlack, kBlack);
// 	DrawGammaSetMarker(histo2, 21, 0.8, kRed, kRed);
    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);

    histo1->GetYaxis()->SetLabelSize(0.035);
    histo1->GetYaxis()->SetTitleSize(0.04);	
    histo1->GetYaxis()->SetDecimals();
    histo1->GetYaxis()->SetTitleOffset(yOffset);
    histo1->GetXaxis()->SetTitleSize(0.04);	
    histo1->GetXaxis()->SetLabelSize(0.035);
    histo1->GetXaxis()->SetTitleOffset(xOffset);
    histo1->DrawCopy("e2,p");
    histo2->SetLineColor(2);
    histo2->DrawCopy("e,hist,same");
    TLegend* leg1 = new TLegend( 0.7,0.87,0.97,0.97);
    leg1->SetTextSize(0.04);			
    leg1->SetFillColor(0);
    leg1->AddEntry(histo1,("Data"));
    leg1->AddEntry(histo2,("MC"));
    leg1->Draw();
    
}

void DrawAutoGammaHistosMaterial(   TH1* histo1, 
                                    TH1*histo2, 
                                    TString Title, 
                                    TString XTitle, 
                                    TString YTitle,
                                    Bool_t YRangeMax, 
                                    Float_t YMaxFactor, 
                                    Float_t YMinimum, 
                                    Bool_t YRange, 
                                    Float_t YMin, 
                                    Float_t YMax, 
                                    Bool_t XRange, 
                                    Float_t XMin, 
                                    Float_t XMax) {

    if (YRangeMax && !XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        if(maxRangeR < histo2->GetMaximum()){
            maxRangeR = histo2->GetMaximum();
        }
        Double_t minRangeR = histo1->GetMinimum();		
        if(minRangeR > histo2->GetMinimum()){
            minRangeR = histo2->GetMinimum();
        }
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
    }
    if (YRangeMax && XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        if(maxRangeR < histo2->GetMaximum()){
            maxRangeR = histo2->GetMaximum();
        }
        Double_t minRangeR = histo1->GetMinimum();		
        if(minRangeR > histo2->GetMinimum()){
            minRangeR = histo2->GetMinimum();
        }
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (YRange && XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);	
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (!YRangeMax && !YRange && XRange){
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (YRange && !XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);
    }
    
    if(Title.CompareTo("") != 0){
        histo1->SetTitle(Title.Data());
    } else{histo1->SetTitle("");
    histo2->SetTitle("");}
    if(XTitle.CompareTo("") != 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.CompareTo("") != 0){
        histo1->SetYTitle(YTitle.Data());
    }
    DrawGammaSetMarker(histo1, 20, 0.5, kBlack, kBlack);
    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);

// 	DrawGammaSetMarker(histo2, 21, 0.8, kRed, kRed);
    histo1->GetYaxis()->SetLabelSize(0.035);
    histo1->GetYaxis()->SetTitleSize(0.04);	
    histo1->GetYaxis()->SetDecimals();
    histo1->GetYaxis()->SetTitleOffset(1.5);
    histo1->GetXaxis()->SetTitleSize(0.04);	
    histo1->GetXaxis()->SetLabelSize(0.035);
    histo1->DrawCopy("e2,hist");
    histo2->SetLineColor(kRed);
    histo2->DrawCopy("e,hist,same");
    TLegend* leg1 = new TLegend( 0.7,0.87,0.97,0.97);
    leg1->SetTextSize(0.04);			
    leg1->SetFillColor(0);
    leg1->AddEntry(histo1,("Data"));
    leg1->AddEntry(histo2,("MC"));
    leg1->Draw();
}

void DrawAutoGammaHistosMaterialP(  TH1* histo1, 
                                    TH1*histo2, 
                                    TString Title, 
                                    TString XTitle, 
                                    TString YTitle,
                                    Bool_t YRangeMax, 
                                    Float_t YMaxFactor,
                                    Float_t YMinimum, 
                                    Bool_t YRange, 
                                    Float_t YMin, 
                                    Float_t YMax, 
                                    Bool_t XRange, 
                                    Float_t XMin, 
                                    Float_t XMax) {
    if (YRangeMax && !XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        if(maxRangeR < histo2->GetMaximum()){
            maxRangeR = histo2->GetMaximum();
        }
        Double_t minRangeR = histo1->GetMinimum();		
        if(minRangeR > histo2->GetMinimum()){
            minRangeR = histo2->GetMinimum();
        }
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
    }
    if (YRangeMax && XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        if(maxRangeR < histo2->GetMaximum()){
            maxRangeR = histo2->GetMaximum();
        }
        Double_t minRangeR = histo1->GetMinimum();		
        if(minRangeR > histo2->GetMinimum()){
            minRangeR = histo2->GetMinimum();
        }
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (YRange && XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);	
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (!YRangeMax && !YRange && XRange){
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (YRange && !XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);
    }

    if(Title.CompareTo("") != 0){
        histo1->SetTitle(Title.Data());
    }else{histo1->SetTitle("");
    histo2->SetTitle("");}
    if(XTitle.CompareTo("") != 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.CompareTo("") != 0){
        histo1->SetYTitle(YTitle.Data());
    }
    DrawGammaSetMarker(histo1, 20, 0.5, kBlack, kBlack);
    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);

    // 	DrawGammaSetMarker(histo2, 21, 0.8, kRed, kRed);
    histo1->GetYaxis()->SetLabelSize(0.035);
    histo1->GetYaxis()->SetTitleSize(0.04);	
    histo1->GetYaxis()->SetDecimals();
    histo1->GetYaxis()->SetTitleOffset(1.5);
    histo1->GetXaxis()->SetTitleSize(0.04);	
    histo1->GetXaxis()->SetLabelSize(0.035);
    histo1->DrawCopy("e2,hist");
    histo2->SetLineColor(kRed);
    histo2->DrawCopy("e,histo,same");
    TLegend* leg1 = new TLegend( 0.7,0.87,0.97,0.96);
    leg1->SetTextSize(0.04);			
    leg1->SetFillColor(0);
    // 	leg1->SetLineColor(0);
    leg1->AddEntry(histo1,("Data"));
    leg1->AddEntry(histo2,("MC"));
    leg1->Draw();

}


void DrawAutoGamma3Histos(  TH1* histo1, 
                            TH1* histo2, 
                            TH1* histo3, 
                            TString Title, 
                            TString XTitle, 
                            TString YTitle,
                            Bool_t YRangeMax, 
                            Float_t YMaxFactor, 
                            Float_t YMinimum, 
                            Bool_t YRange, 
                            Float_t YMin, 
                            Float_t YMax, 
                            Bool_t XRange, 
                            Float_t XMin, 
                            Float_t XMax) {

    if (YRangeMax && !XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        if(maxRangeR < histo2->GetMaximum()){
            maxRangeR = histo2->GetMaximum();
        }
        Double_t minRangeR = histo1->GetMinimum();		
        if(minRangeR > histo2->GetMinimum()){
            minRangeR = histo2->GetMinimum();
        }
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
    }
    if (YRangeMax && XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        if(maxRangeR < histo2->GetMaximum()){
            maxRangeR = histo2->GetMaximum();
        }
        Double_t minRangeR = histo1->GetMinimum();		
        if(minRangeR > histo2->GetMinimum()){
            minRangeR = histo2->GetMinimum();
        }
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (YRange && XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);	
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (!YRangeMax && !YRange && XRange){
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (YRange && !XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);
    }

    if(Title.CompareTo("") != 0){
        histo1->SetTitle(Title.Data());
    }else{histo1->SetTitle("");
    histo2->SetTitle("");}
    if(XTitle.CompareTo("") != 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.CompareTo("") != 0){
        histo1->SetYTitle(YTitle.Data());
    }
    DrawGammaSetMarker(histo1, 20, 0.3, kBlack, kBlack);
    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);

    // 	DrawGammaSetMarker(histo2, 21, 0.8, kRed, kRed);
    // 	DrawGammaSetMarker(histo3, 22, 0.8, kBlue, kBlue);
    histo1->GetYaxis()->SetLabelSize(0.035);
    histo1->GetYaxis()->SetTitleSize(0.04);	
    histo1->GetYaxis()->SetDecimals();
    histo1->GetYaxis()->SetTitleOffset(1.5);
    histo1->GetXaxis()->SetTitleSize(0.04);	
    histo1->GetXaxis()->SetLabelSize(0.035);
    histo1->DrawCopy("e,hist");
    histo2->SetLineColor(kRed);
    histo2->DrawCopy("e,hist,same");
    histo3->SetLineColor(kYellow-7);
    histo3->SetFillColor(kYellow-7);
    histo3->Draw("same,hist");
    histo2->DrawCopy("e,hist,same");
    histo1->DrawCopy("e,hist,same");
    histo1->DrawCopy("same,axis");
    TLegend* leg1 = new TLegend( 0.7,0.82,0.97,0.97);
    leg1->SetTextSize(0.04);			
    leg1->SetFillColor(0);
    leg1->AddEntry(histo1,("Data"));
    leg1->AddEntry(histo2,("MC"));
    leg1->AddEntry(histo3,"True conversion","f");
    leg1->Draw();

}


void DrawAutoGammaHistosWOLeg(  TH1* histo1, 
                                TH1*histo2, 
                                TString Title, 
                                TString XTitle, 
                                TString YTitle,
                                Bool_t YRangeMax, 
                                Float_t YMaxFactor, 
                                Float_t YMinimum, 
                                Bool_t YRange, 
                                Float_t YMin, 
                                Float_t YMax, 
                                Bool_t XRange, 
                                Float_t XMin, 
                                Float_t XMax) {
    if (YRangeMax && !XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        if(maxRangeR < histo2->GetMaximum()){
            maxRangeR = histo2->GetMaximum();
        }
        Double_t minRangeR = histo1->GetMinimum();		
        if(minRangeR > histo2->GetMinimum()){
            minRangeR = histo2->GetMinimum();
        }
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
    }
    if (YRangeMax && XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        if(maxRangeR < histo2->GetMaximum()){
            maxRangeR = histo2->GetMaximum();
        }
        Double_t minRangeR = histo1->GetMinimum();		
        if(minRangeR > histo2->GetMinimum()){
            minRangeR = histo2->GetMinimum();
        }
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (YRange && XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);	
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (!YRangeMax && !YRange && XRange){
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (YRange && !XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);
    }

    if(Title.CompareTo("") != 0){
        histo1->SetTitle(Title.Data());
    }else{histo1->SetTitle("");
    histo2->SetTitle("");}
    if(XTitle.CompareTo("") != 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.CompareTo("") != 0){
        histo1->SetYTitle(YTitle.Data());
    }
    DrawGammaSetMarker(histo1, 20, 0.8, kBlack, kBlack);
    // 	DrawGammaSetMarker(histo2, 21, 0.8, kRed, kRed);
    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);

    histo1->GetYaxis()->SetLabelSize(0.035);
    histo1->GetYaxis()->SetTitleSize(0.04);	
    histo1->GetYaxis()->SetDecimals();
    histo1->GetYaxis()->SetTitleOffset(1.4);
    histo1->GetXaxis()->SetTitleSize(0.04);	
    histo1->GetXaxis()->SetLabelSize(0.035);
    histo1->DrawCopy("e2,p");
    histo2->SetLineColor(2);
    histo2->DrawCopy("e,hist,same");

}


/* DrawAutoGammaHisto is function used for styling a histograma of the gamma conversion group with standart settings
* histo1 - first histogram (Data)
* Title - histogram title
* XTitle - X-axis title
* YTitle - Y-axis title
* YRangeMax 	= kTRUE will scale by Maximum and Minimum Range in Y
*YMaxFactor - will MaximumY by this factor if YRangeMay = kTRUE 
*YMinimum - this will be used if YRangeMax is set
*YRange  	= kTRUE will Cut y-axis by YMin and YMax 
- will be set to kFAlSE if YRangeMax is set
*YMin - minimum Y
*YMax - maximum Y
*XRange 	= kTRUE will Cut x-axis by XMin and XMax
*XMin - minimum Y
*XMax - maximum Y
*/ 
void DrawAutoGammaHisto(    TH1* histo1,
                            TString Title, 
                            TString XTitle, 
                            TString YTitle,
                            Bool_t YRangeMax, 
                            Double_t YMaxFactor, 
                            Double_t YMinimum,
                            Bool_t YRange, 
                            Double_t YMin, 
                            Double_t YMax,
                            Bool_t XRange, 
                            Double_t XMin, 
                            Double_t XMax, 
                            Double_t yOffset=1.) {
    if (YRangeMax && !XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        Double_t minRangeR = histo1->GetMinimum();
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);
    }
    if (YRangeMax && XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        Double_t minRangeR = histo1->GetMinimum();
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);
    }
    if (YRange && XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);
    }
    if (!YRangeMax && !YRange && XRange){
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);
    }
    
    if (YRange && !XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);
    }
    
    
    histo1->SetTitle(Title.Data());
    
    if(XTitle.CompareTo("") != 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.CompareTo("") != 0){
        histo1->SetYTitle(YTitle.Data());
    }

    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);
    histo1->GetYaxis()->SetLabelSize(0.035);
    histo1->GetYaxis()->SetTitleSize(0.043);
    histo1->GetYaxis()->SetDecimals();
    histo1->GetYaxis()->SetTitleOffset(yOffset);
    histo1->GetXaxis()->SetTitleSize(0.043);
    histo1->GetXaxis()->SetLabelSize(0.035);
    histo1->DrawCopy("e,hist");
}

/*DrawAutoGammaHisto2D is a function for drawing a 2D-histogram of the gamma conversion group
* histo - histogramm which need to be drawn
* Title - histogram title
* XTitle - X- axis-title
* YTitle - Y-axis-title
* Input - Legend 
* YRange - if kTRUE will scale by YMin and YMay
* YMin  - Y minimum
* YMax - Y maximum
* XRange - if kTRUE will scale by XMin and XMax
* XMin - X minimum
* XMax - X maximum
*/
void DrawAutoGammaHisto2D(  TH2 *histo,  
                            TString Title, 
                            TString XTitle, 
                            TString YTitle, 
                            TString Input,
                            Bool_t YRange, 
                            Float_t YMin, 
                            Float_t YMax, 
                            Bool_t XRange, 
                            Float_t XMin, 
                            Float_t XMax, 
                            Double_t titleOffsetX = 1.2, 
                            Double_t titleOffsetY = 1.4,
                            Size_t labelSizeX = 0.035, 
                            Size_t titleSizeX = 0.043, 
                            Size_t labelSizeY = 0.035, 
                            Size_t titleSizeY = 0.043){
	
	
    if (YRange && XRange){
        histo->GetYaxis()->SetRangeUser(YMin, YMax);	
        histo->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if ( !YRange && XRange){
        histo->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    
    if (YRange && !XRange){
        histo->GetYaxis()->SetRangeUser(YMin, YMax);
    }
    
    if(Title.CompareTo("") != 0){
        histo->SetTitle(Title.Data());
    }
    if(XTitle.CompareTo("") != 0){
        histo->SetXTitle(XTitle.Data());
    }
    if(YTitle.CompareTo("") != 0){
        histo->SetYTitle(YTitle.Data());
    }
    histo->GetYaxis()->SetLabelFont(42);
    histo->GetXaxis()->SetLabelFont(42);
    histo->GetYaxis()->SetTitleFont(62);
    histo->GetXaxis()->SetTitleFont(62);

    histo->GetYaxis()->SetLabelSize(labelSizeY);
    histo->GetYaxis()->SetTitleSize(titleSizeY);	
    histo->GetYaxis()->SetTitleOffset(titleOffsetY);
    histo->GetYaxis()->SetDecimals();
    
    histo->GetXaxis()->SetLabelSize(labelSizeX);
    histo->GetXaxis()->SetTitleSize(titleSizeX);
    histo->GetXaxis()->SetTitleOffset(titleOffsetX);
    histo->DrawCopy("colz");
    if(Input.CompareTo("") != 0){
        TLegend* leg2 = new TLegend(0.6,0.82,0.83,0.9);
        leg2->SetTextSize(0.04);			
        leg2->SetFillColor(0);
        leg2->AddEntry(histo,(Input.Data()));
        leg2->Draw("same");
    }
}


/* DrawRatioGammaHisto is function used for styling the ratio-histograms of the gamma conversion group
* histo1 - histogram 
* Title - histogram title
* XTitle - X-axis title
* YTitle - Y-axis title
* YRangeMax 	= kTRUE will scale by Maximum and Minimum Range in Y
*YMaxFactor - will MaximumY by this factor if YRangeMay = kTRUE 
*YMinimum - this will be used if YRangeMax is set
*YRange  	= kTRUE will Cut y-axis by YMin and YMax 
- will be set to kFAlSE if YRangeMax is set
*YMin - minimum Y
*YMax - maximum Y
*XRange 	= kTRUE will Cut x-axis by XMin and XMax
*XMin - minimum Y
*XMax - maximum Y
*/ 
void DrawRatioGammaHisto(   TH1* histo1, 
                            TString Title, 
                            TString XTitle, 
                            TString YTitle,
                            Bool_t YRangeMax, 
                            Float_t YMaxFactor, 
                            Float_t YMinimum,
                            Bool_t YRange, 
                            Float_t YMin, 
                            Float_t YMax,  
                            Bool_t XRange, 
                            Float_t XMin, 
                            Float_t XMax) {
    if (YRangeMax && !XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        Double_t minRangeR = histo1->GetMinimum();		
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
    }
    if (YRangeMax && XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        Double_t minRangeR = histo1->GetMinimum();		
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (YRange && XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);	
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (!YRangeMax && !YRange && XRange){
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    
    if (YRange && !XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);
    }
    
    if(Title.CompareTo("") != 0){	histo1->SetTitle(Title.Data());
    }else{	histo1->SetTitle("");}
    
    if(XTitle.CompareTo("") != 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.CompareTo("") != 0){
        histo1->SetYTitle(YTitle.Data());
    }
    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);

    histo1->GetYaxis()->SetTitleSize(0.04);	
    histo1->GetYaxis()->SetLabelSize(0.03);
    histo1->GetYaxis()->SetDecimals();
    histo1->GetYaxis()->SetTitleOffset(0.9);
    histo1->GetXaxis()->SetTitleOffset(0.9);
    histo1->GetXaxis()->SetTitleSize(0.04);
    histo1->GetXaxis()->SetLabelSize(0.03);	
    histo1->SetLineColor(kBlue-5);
    histo1->SetMarkerStyle(20);
    histo1->SetMarkerSize(0.5);
    histo1->DrawCopy("hist,e");
    histo1->DrawCopy("same,p");
}

/* DrawCutGammaHistos is function used for styling the Cut-histograms of the gamma conversion group for 4 histos combined
* histo1 - histogram Data
* histo2 - histogram Data Comparision
* histo3 - histogram MC
* histo4 - histogram MC Comparision
* Title - histogram title
* XTitle - X-axis title
* YTitle - Y-axis title
* Legend1 - additional Legend for histo2
* Legend2 - additional Legend for histo4	
* YRangeMax 	= kTRUE will scale by Maximum and Minimum Range in Y
*YMaxFactor - will MaximumY by this factor if YRangeMay = kTRUE 
*YMinimum - this will be used if YRangeMax is set
*YRange  	= kTRUE will Cut y-axis by YMin and YMax 
- will be set to kFAlSE if YRangeMax is set
*YMin - minimum Y
*YMax - maximum Y
*XRange 	= kTRUE will Cut x-axis by XMin and XMax
*XMin - minimum Y
*XMax - maximum Y
*/ 
void DrawCutGammaHistos(    TH1* histo1, 
                            TH1* histo2, 
                            TH1* histo3, 
                            TH1*histo4, 
                            TString Title, 
                            TString XTitle, 
                            TString YTitle, 
                            const char *Legend1, 
                            const char *Legend2,
                            Bool_t YRangeMax, 
                            Float_t YMaxFactor, 
                            Float_t YMinimum,
                            Bool_t YRange, 
                            Float_t YMin, 
                            Float_t YMax,  
                            Bool_t XRange, 
                            Float_t XMin, 
                            Float_t XMax) {
    if (YRangeMax && !XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo2->GetMaximum();
        if(maxRangeR < histo4->GetMaximum()){
            maxRangeR = histo4->GetMaximum();
        }
        Double_t minRangeR = histo2->GetMinimum();		
        if(minRangeR > histo4->GetMinimum()){
            minRangeR = histo4->GetMinimum();
        }
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
    }
    if (YRangeMax && XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo2->GetMaximum();
        if(maxRangeR < histo4->GetMaximum()){
            maxRangeR = histo4->GetMaximum();
        }
        Double_t minRangeR = histo2->GetMinimum();		
        if(minRangeR > histo4->GetMinimum()){
            minRangeR = histo4->GetMinimum();
        }
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (YRange && XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);	
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (!YRangeMax && !YRange && XRange){
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (YRange && !XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);
    }
    
    if(Title.CompareTo("") != 0){
        histo1->SetTitle(Title.Data());
    }
    if(XTitle.CompareTo("") != 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.CompareTo("") != 0){
        histo1->SetYTitle(YTitle.Data());
    }
    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);

    histo1->GetYaxis()->SetTitleSize(0.025);		
    histo1->GetYaxis()->SetLabelSize(0.02);
    histo1->GetYaxis()->SetDecimals();
    histo1->GetYaxis()->SetTitleOffset(1.8);
    histo1->GetXaxis()->SetLabelSize(0.02);
    histo1->GetXaxis()->SetTitleSize(0.025);	
    histo1->Draw("e,hist");
    
    histo2->SetLineColor(15);
    histo2->Draw("e,hist,same");
    
    histo3->SetLineColor(2);
    histo3->Draw("e,hist,same");
    
    histo4->SetLineColor(46);		
    histo4->Draw("e,hist,same");
    
    TLegend* leg1 = new TLegend( 0.6,0.82,0.92,0.9);
    leg1->SetTextSize(0.02);			
    leg1->SetFillColor(0);
    leg1->AddEntry(histo1,("Data"));
    leg1->AddEntry(histo2,(Legend1));
    leg1->AddEntry(histo3,("MC"));
    leg1->AddEntry(histo4,(Legend2));
    
    leg1->Draw();
}

/* DrawCutGammaHisto is function used for styling the Cut-histograms of the gamma conversion group for 2 histos combined
* histo1 - histogram Data
* histo2 - histogram Data Comparision
* Title - histogram title
* XTitle - X-axis title
* YTitle - Y-axis title
* Legend - additional Legend for histo2
* YRangeMax 	= kTRUE will scale by Maximum and Minimum Range in Y
*YMaxFactor - will MaximumY by this factor if YRangeMay = kTRUE 
*YMinimum - this will be used if YRangeMax is set
*YRange  	= kTRUE will Cut y-axis by YMin and YMax 
- will be set to kFAlSE if YRangeMax is set
*YMin - minimum Y
*YMax - maximum Y
*XRange 	= kTRUE will Cut x-axis by XMin and XMax
*XMin - minimum Y
*XMax - maximum Y
*/ 
void DrawCutGammaHisto( TH1* histo1, 
                        TH1* histo2, 
                        TString Title, 
                        TString XTitle, 
                        TString YTitle, 
                        const char *Legend,
                        Bool_t YRangeMax, 
                        Float_t YMaxFactor, 
                        Float_t YMinimum,
                        Bool_t YRange, 
                        Float_t YMin, 
                        Float_t YMax,  
                        Bool_t XRange, 
                        Float_t XMin, 
                        Float_t XMax) {
    if (YRangeMax && !XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo2->GetMaximum();
        Double_t minRangeR = histo2->GetMinimum();		
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
    }
    if (YRangeMax && XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo2->GetMaximum();
        Double_t minRangeR = histo2->GetMinimum();				
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (YRange && XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);	
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (!YRangeMax && !YRange && XRange){
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (YRange && !XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);
    }
    
    if(Title.CompareTo("") != 0){
        histo1->SetTitle(Title.Data());
    }
    if(XTitle.CompareTo("") != 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.CompareTo("") != 0){
        histo1->SetYTitle(YTitle.Data());
    }
    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);

    histo1->GetYaxis()->SetTitleSize(0.025);	
    histo1->GetYaxis()->SetLabelSize(0.02);
    histo1->GetYaxis()->SetDecimals();
    histo1->GetYaxis()->SetTitleOffset(1.8);
    histo1->GetXaxis()->SetLabelSize(0.02);
    histo1->GetXaxis()->SetTitleSize(0.025);	
    histo1->Draw("e,hist");
    
    histo2->SetLineColor(15);
    histo2->Draw("e,hist,same");
    
    TLegend* leg1 = new TLegend( 0.6,0.82,0.92,0.9);
    leg1->SetTextSize(0.04);			
    leg1->SetFillColor(0);
    leg1->AddEntry(histo1,("Data"));
    leg1->AddEntry(histo2,(Legend));
    leg1->Draw();	
}

/* DrawRatioGammaHisto is function used for styling the ratio-histograms of the gamma conversion group
* histo1 - histogram 
* Title - histogram title
* XTitle - X-axis title
* YTitle - Y-axis title
* YRangeMax 	= kTRUE will scale by Maximum and Minimum Range in Y
*YMaxFactor - will MaximumY by this factor if YRangeMay = kTRUE 
*YMinimum - this will be used if YRangeMax is set
*YRange  	= kTRUE will Cut y-axis by YMin and YMax 
- will be set to kFAlSE if YRangeMax is set
*YMin - minimum Y
*YMax - maximum Y
*XRange 	= kTRUE will Cut x-axis by XMin and XMax
*XMin - minimum Y
*XMax - maximum Y
*/ 
void DrawResolutionGammaHisto(  TH1* histo1, 
                                TString Title, 
                                TString XTitle, 
                                TString YTitle,
                                Bool_t YRangeMax, 
                                Float_t YMaxFactor, 
                                Float_t YMinimum,
                                Bool_t YRange, 
                                Float_t YMin, 
                                Float_t YMax,  
                                Bool_t XRange, 
                                Float_t XMin, 
                                Float_t XMax) {
    if (YRangeMax && !XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        Double_t minRangeR = histo1->GetMinimum();		
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
    }
    if (YRangeMax && XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        Double_t minRangeR = histo1->GetMinimum();		
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (YRange && XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);	
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (!YRangeMax && !YRange && XRange){
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    
    if (YRange && !XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);
    }
    
    if(Title.CompareTo("") != 0){
        histo1->SetTitle(Title.Data());
    }else { histo1->SetTitle("");}
    
    if(XTitle.CompareTo("") != 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.CompareTo("") != 0){
        histo1->SetYTitle(YTitle.Data());
    }
    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);

    histo1->GetYaxis()->SetTitleSize(0.055);	
    histo1->GetYaxis()->SetLabelSize(0.045);
    histo1->GetYaxis()->SetDecimals();
    histo1->GetYaxis()->SetTitleOffset(0.9);
    histo1->GetXaxis()->SetTitleOffset(0.85);
    histo1->GetXaxis()->SetTitleSize(0.055);
    histo1->GetXaxis()->SetLabelSize(0.045);	
    histo1->DrawCopy("e1");
}

/*DrawAutoGammaHisto2Dres is a function for drawing a resolution 2D-histogram of the gamma conversion group
* histo - histogramm which need to be drawn
* Title - histogram title
* XTitle - X- axis-title
* YTitle - Y-axis-title
* Input - Legend 
* YRange - if kTRUE will scale by YMin and YMay
* YMin  - Y minimum
* YMax - Y maximum
* XRange - if kTRUE will scale by XMin and XMax
* XMin - X minimum
* XMax - X maximum
*/
void DrawAutoGammaHisto2DRes(   TH2 *histo,  
                                TString Title, 
                                TString XTitle, 
                                TString YTitle, 
                                TString Input,
                                Bool_t YRange,
                                Float_t YMin, 
                                Float_t YMax, 
                                Bool_t XRange, 
                                Float_t XMin, 
                                Float_t XMax) {
	
    if (YRange && XRange){
        histo->GetYaxis()->SetRangeUser(YMin, YMax);	
        histo->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if ( !YRange && XRange){
        histo->GetXaxis()->SetRangeUser(XMin, XMax);	
    }

    if (YRange && !XRange){
        histo->GetYaxis()->SetRangeUser(YMin, YMax);
    }

    if(Title.CompareTo("") != 0){
        histo->SetTitle(Title.Data());
    }
    if(XTitle.CompareTo("") != 0){
        histo->SetXTitle(XTitle.Data());
    }
    if(YTitle.CompareTo("") != 0){
        histo->SetYTitle(YTitle.Data());
    }
    histo->GetYaxis()->SetLabelFont(42);
    histo->GetXaxis()->SetLabelFont(42);
    histo->GetYaxis()->SetTitleFont(62);
    histo->GetXaxis()->SetTitleFont(62);

    histo->GetYaxis()->SetTitleSize(0.045);	
    histo->GetYaxis()->SetLabelSize(0.03);
    histo->GetXaxis()->SetLabelSize(0.03);
    histo->GetYaxis()->SetDecimals();
    histo->GetYaxis()->SetTitleOffset(1.5);
    histo->GetXaxis()->SetTitleSize(0.045);	
    histo->GetYaxis()->SetTitleOffset(1.5);
    histo->DrawCopy("colz");
    if(Input.CompareTo("") != 0){
        TLegend* leg2 = new TLegend(0.6,0.82,0.83,0.9);
        leg2->SetTextSize(0.04);			
        leg2->SetFillColor(0);
        leg2->AddEntry(histo,(Input.Data()));
        leg2->Draw("same");
    }
}


/* DrawAutoGammaMesonHistos is function used for styling the histograms of the gamma conversion group for two histos and standart settings
* histo1 - first histogram
* Title - histogram title
* XTitle - X-axis title
* YTitle - Y-axis title
* YRangeMax 	= kTRUE will scale by Maximum and Minimum Range in Y
*YMaxFactor - will MaximumY by this factor if YRangeMay = kTRUE 
*YMinimum - this will be used if YRangeMax is set
*YRange  	= kTRUE will Cut y-axis by YMin and YMax 
- will be set to kFAlSE if YRangeMax is set
*YMin - minimum Y
*YMax - maximum Y
*XRange 	= kTRUE will Cut x-axis by XMin and XMax
*XMin - minimum Y
*XMax - maximum Y
*/ 

void DrawAutoGammaMesonHistos(  TH1* histo1, 
                                TString Title, 
                                TString XTitle, 
                                TString YTitle, 
                                Bool_t YRangeMax, 
                                Float_t YMaxFactor, 
                                Float_t YMinimum, 
                                Bool_t ScaleByMaxPtBin, 
                                Bool_t YRange, 
                                Float_t YMin, 
                                Float_t YMax, 
                                Bool_t XRange, 
                                Float_t XMin, 
                                Float_t XMax,
                                Style_t textFontTitle = 62, 
                                Size_t textSizeTitle = 0.04, 
                                Style_t textFontLabel = 42, 
                                Size_t textSizeLabel = 0.03,
                                Double_t offsetTitleX = 0.9, 
                                Double_t offsetTitleY = 1.2) {
    if (YRangeMax && !XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        Double_t minRangeR;
        if (ScaleByMaxPtBin) { 
            minRangeR = 0.05*histo1->GetBinContent(histo1->GetNbinsX());
        } else { 
            minRangeR = 0.1*histo1->GetBinContent(histo1->GetMinimumBin());
        }
        cout << maxRangeR << "\t" << minRangeR << endl;
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
    }
    if (YRangeMax && XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        Double_t minRangeR = histo1->GetBinContent(histo1->GetMinimumBin())/10.;
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (YRange && XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);	
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (!YRangeMax && !YRange && XRange){
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (YRange && !XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);
    }
    
    if(Title.CompareTo("") != 0){
        histo1->SetTitle(Title.Data());
    }else{
        histo1->SetTitle("");
    }
    if(XTitle.CompareTo("") != 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.CompareTo("") != 0){
        histo1->SetYTitle(YTitle.Data());
    }
    histo1->GetYaxis()->SetLabelFont(textFontLabel);
    histo1->GetXaxis()->SetLabelFont(textFontLabel);
    histo1->GetYaxis()->SetTitleFont(textFontTitle);
    histo1->GetXaxis()->SetTitleFont(textFontTitle);

    histo1->GetYaxis()->SetLabelSize(textSizeLabel);
    histo1->GetYaxis()->SetTitleSize(textSizeTitle);	
    histo1->GetYaxis()->SetDecimals();
    histo1->GetYaxis()->SetTitleOffset(offsetTitleY);
    
    histo1->GetXaxis()->SetTitleSize(textSizeTitle);	
    histo1->GetXaxis()->SetLabelSize(textSizeLabel);
    histo1->GetXaxis()->SetTitleOffset(offsetTitleX);
    
    histo1->DrawCopy("e1,p");
	
}


void DrawGammaCanvasSettings( TCanvas* c1, 
                              Double_t leftMargin, 
                              Double_t rightMargin, 
                              Double_t topMargin, 
                              Double_t bottomMargin){
    c1->SetTickx();
    c1->SetTicky();
    c1->SetGridx(0);
    c1->SetGridy(0);
    c1->SetLogy(0);	
    c1->SetLeftMargin(leftMargin);
    c1->SetRightMargin(rightMargin);
    c1->SetTopMargin(topMargin);				
    c1->SetBottomMargin(bottomMargin);				
    c1->SetFillColor(0);
}

void DrawGammaPadSettings( TPad* pad1, 
                           Double_t leftMargin, 
                           Double_t rightMargin, 
                           Double_t topMargin, 
                           Double_t bottomMargin){
    pad1->SetFillColor(0);
    pad1->GetFrame()->SetFillColor(0);
    pad1->SetBorderMode(0);
    pad1->SetLeftMargin(leftMargin);
    pad1->SetBottomMargin(bottomMargin);
    pad1->SetRightMargin(rightMargin);
    pad1->SetTopMargin(topMargin);
    pad1->SetTickx();
    pad1->SetTicky();
}

void DrawGammaSetMarkerTGraph(  TGraph* graph, 
                                Style_t markerStyle, 
                                Size_t markerSize, 
                                Color_t markerColor, 
                                Color_t lineColor, 
                                Width_t lineWidth = 1, 
                                Style_t lineStyle = 1, 
                                Bool_t boxes = kFALSE, 
                                Color_t fillColor = 0) {
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(markerSize);
    graph->SetMarkerColor(markerColor);
    graph->SetLineColor(lineColor);
    graph->SetLineWidth(lineWidth);
    graph->SetLineWidth(lineStyle);
    if (boxes){
        graph->SetFillColor(fillColor);
        if (fillColor!=0){
            graph->SetFillStyle(1001);
        } else {
            graph->SetFillStyle(0);
        }
        
    }
}

void DrawGammaSetMarkerTGraphErr(   TGraphErrors* graph, 
                                    Style_t markerStyle, 
                                    Size_t markerSize, 
                                    Color_t markerColor, 
                                    Color_t lineColor, 
                                    Width_t lineWidth=1, 
                                    Bool_t boxes = kFALSE, 
                                    Color_t fillColor = 0) {
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(markerSize);
    graph->SetMarkerColor(markerColor);
    graph->SetLineColor(lineColor);	
    graph->SetLineWidth(lineWidth);
    if (boxes){
        graph->SetFillColor(fillColor);
        if (fillColor!=0){
            graph->SetFillStyle(1001);
        } else {
            graph->SetFillStyle(0);
        }
    }
}

void DrawGammaNLOTGraph( TGraph* graph, 
                         Width_t lineWidth, 
                         Style_t lineStyle, 
                         Color_t lineColor){
    graph->SetLineWidth(lineWidth);
    graph->SetLineColor(lineColor);
    graph->SetLineStyle(lineStyle);
}

void SetStyleGammaNLOTGraphWithBand( TGraph* graph, 
                                     Width_t lineWidth, 
                                     Style_t lineStyle, 
                                     Color_t lineColor, 
                                     Style_t fillStyle, 
                                     Color_t fillColor, 
                                     Size_t markerSize){
    graph->SetMarkerColor(lineColor);
    graph->SetLineColor(lineColor);
    graph->SetFillColor(fillColor);
    graph->SetFillStyle(fillStyle);
    graph->SetLineWidth(lineWidth);
    graph->SetLineStyle(lineStyle);
    graph->SetMarkerSize(markerSize);	
}


void DrawGammaSetMarkerTGraphAsym(  TGraphAsymmErrors* graph, 
                                    Style_t markerStyle, 
                                    Size_t markerSize, 
                                    Color_t markerColor, 
                                    Color_t lineColor,  
                                    Width_t lineWidth=1, 
                                    Bool_t boxes = kFALSE, 
                                    Color_t fillColor = 0) {
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(markerSize);
    graph->SetMarkerColor(markerColor);
    graph->SetLineColor(lineColor);	
    graph->SetLineWidth(lineWidth);
    if (boxes){
        graph->SetFillColor(fillColor);
        if (fillColor!=0){
            graph->SetFillStyle(1001);
        } else {
            graph->SetFillStyle(0);
        }
        
    }
}


void DrawGammaSetMarkerTF1( TF1* fit1, 
                            Style_t lineStyle, 
                            Size_t lineWidth, 
                            Color_t lineColor ) {
    fit1->SetLineColor(lineColor);	
    fit1->SetLineStyle(lineStyle);	
    fit1->SetLineWidth(lineWidth);	
}

void SetStyleTLatex( TLatex* text, 
                     Size_t textSize, 
                     Width_t lineWidth, 
                     Color_t textColor = 1, 
                     Font_t textFont = 42, 
                     Bool_t kNDC = kTRUE){
    if (kNDC) {text->SetNDC();}
    text->SetTextFont(textFont);
    text->SetTextColor(textColor);
    text->SetTextSize(textSize);
    text->SetLineWidth(lineWidth);
}

void SetStyleHisto( TH1* histo, 
                    Width_t lineWidth, 
                    Style_t lineStyle, 
                    Color_t lineColor) { 
    histo->SetLineWidth(lineWidth);
    histo->SetLineStyle(lineStyle);
    histo->SetLineColor(lineColor);
}

void SetStyleFit(   TF1* fit, 
                    Double_t xRangeStart, 
                    Double_t xRangeEnd, 
                    Width_t lineWidth, 
                    Style_t lineStyle, 
                    Color_t lineColor) { 
    fit->SetRange(xRangeStart,xRangeEnd);
    fit->SetLineWidth(lineWidth);
    fit->SetLineStyle(lineStyle);
    fit->SetLineColor(lineColor);
}

void SetStyleHistoTH2ForGraphs( TH2* histo, 
                                TString XTitle, 
                                TString YTitle, 
                                Size_t xLableSize, 
                                Size_t xTitleSize, 
                                Size_t yLableSize, 
                                Size_t yTitleSize, 
                                Float_t xTitleOffset = 1, 
                                Float_t yTitleOffset = 1, 
                                Int_t xNDivisions = 510, 
                                Int_t yNDivisions = 510){
    histo->SetXTitle(XTitle);
    histo->SetYTitle(YTitle);
    histo->SetTitle("");

    histo->GetXaxis()->SetLabelSize(xLableSize);
    histo->GetXaxis()->SetTitleSize(xTitleSize);
    histo->GetXaxis()->SetTitleOffset(xTitleOffset);
    histo->GetXaxis()->SetNdivisions(xNDivisions,kTRUE);

    histo->GetXaxis()->SetLabelFont(42);
    histo->GetYaxis()->SetLabelFont(42); 
    histo->GetXaxis()->SetTitleFont(62);
    histo->GetYaxis()->SetTitleFont(62);


    histo->GetYaxis()->SetDecimals();
    histo->GetYaxis()->SetLabelSize(yLableSize);
    histo->GetYaxis()->SetTitleSize(yTitleSize);
    histo->GetYaxis()->SetTitleOffset(yTitleOffset);
    histo->GetYaxis()->SetNdivisions(yNDivisions,kTRUE);	
}

void SetStyleHistoTH1ForGraphs( TH1* histo, 
                                TString XTitle, 
                                TString YTitle, 
                                Size_t xLableSize, 
                                Size_t xTitleSize, 
                                Size_t yLableSize, 
                                Size_t yTitleSize, 
                                Float_t xTitleOffset = 1, 
                                Float_t yTitleOffset = 1, 
                                Int_t xNDivisions = 510, 
                                Int_t yNDivisions = 510){
    histo->SetXTitle(XTitle);
    histo->SetYTitle(YTitle);
    histo->SetTitle("");

    histo->GetYaxis()->SetLabelFont(42);
    histo->GetXaxis()->SetLabelFont(42);
    histo->GetYaxis()->SetTitleFont(62);
    histo->GetXaxis()->SetTitleFont(62);

    
    histo->GetXaxis()->SetLabelSize(xLableSize);
    histo->GetXaxis()->SetTitleSize(xTitleSize);
    histo->GetXaxis()->SetTitleOffset(xTitleOffset);
    histo->GetXaxis()->SetNdivisions(xNDivisions,kTRUE);
    
    histo->GetYaxis()->SetDecimals();
    histo->GetYaxis()->SetLabelSize(yLableSize);
    histo->GetYaxis()->SetTitleSize(yTitleSize);
    histo->GetYaxis()->SetTitleOffset(yTitleOffset);
    histo->GetYaxis()->SetNdivisions(yNDivisions,kTRUE);  
}

void DrawGammaHistoWithTitleAndFit(     TH1* histo1, 
                                        TH1* histo2, 
                                        TF1* fit1, 
                                        TF1* fit2, 
                                        TString Title, 
                                        TString XTitle, 
                                        TString YTitle,
                                        Float_t xMin,
                                        Float_t xMax, 
                                        Float_t yMin) {

    histo1->GetXaxis()->SetRangeUser(xMin, xMax);
    histo1->GetYaxis()->SetRangeUser(yMin, 2.5*histo1->GetMaximum());
    if(Title.Length() > 0){
        histo1->SetTitle("");
        TLatex *alice = new TLatex(0.1,0.95,Form("%s",Title.Data())); // Bo: this was 
        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextSize(0.062);
        alice->Draw();		
    }
    if(XTitle.Length() > 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.Length() > 0){
        histo1->SetYTitle(YTitle.Data());
    }
    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);

    histo1->GetYaxis()->SetLabelSize(0.02);
    histo1->GetYaxis()->SetTitleSize(0.025);
    histo1->GetYaxis()->SetDecimals();
    histo1->GetXaxis()->SetTitleSize(0.025);
    histo1->GetXaxis()->SetLabelSize(0.02);
    histo1->SetMarkerStyle(20);
    histo1->SetMarkerColor(1);
    histo1->SetLineColor(1);
    histo1->SetLineWidth(1);
    histo1->SetMarkerSize(0.3);
    histo1->SetMarkerStyle(20);
    histo1->SetTitleOffset(1.4,"xy");		
    histo1->SetTitleSize(0.05,"xy");		
    histo1->GetYaxis()->SetLabelSize(0.05);
    histo1->GetXaxis()->SetLabelSize(0.05);
    histo1->GetXaxis()->SetNdivisions(507,kTRUE);
    histo1->DrawCopy("p,e1");
    if(Title.Length() > 0){
        histo1->SetTitle("");
        TLatex *alice = new TLatex(0.1,0.95,Form("%s",Title.Data())); // Bo: this was 
        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextSize(0.062);
        alice->Draw();		
    }
    histo2->SetLineStyle(1);		
    histo2->SetLineColor(2);
    histo2->SetMarkerColor(2);
    histo2->SetMarkerSize(0.3);
    histo2->SetMarkerStyle(20);
    histo2->SetLineWidth(1);
    histo2->DrawCopy("p,e1,same");
    if (fit1 != 0x0){
        fit1->SetLineColor(4);
        fit1->SetLineWidth(1.);
        fit1->Draw("same");
    }
    if (fit2 != 0x0){
        fit2->SetLineColor(kGreen+2);
        fit2->SetLineWidth(1.);
        fit2->Draw("same");
    }
}

void DrawGammaHistoWithTitle2(  TH1* histo1,
                                TString Title, 
                                TString XTitle, 
                                TString YTitle,
                                Float_t xMin,
                                Float_t xMax,
                                Float_t yMin) {
    histo1->GetXaxis()->SetRangeUser(xMin, xMax);
    histo1->GetYaxis()->SetRangeUser(yMin, 1.5*histo1->GetMaximum());
    if(Title.Length() > 0){
        histo1->SetTitle("");
        TLatex *alice = new TLatex(0.1,0.95,Form("%s",Title.Data())); // Bo: this was 
        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextSize(0.062);
        alice->Draw();		
    }
    if(XTitle.Length() > 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.Length() > 0){
        histo1->SetYTitle(YTitle.Data());
    }
    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);

    histo1->GetYaxis()->SetLabelSize(0.02);
    histo1->GetYaxis()->SetTitleSize(0.025);
    histo1->GetYaxis()->SetDecimals();
    histo1->GetXaxis()->SetTitleSize(0.025);
    histo1->GetXaxis()->SetLabelSize(0.02);
    histo1->SetMarkerStyle(20);
    histo1->SetMarkerColor(1);
    histo1->SetLineColor(1);
    histo1->SetLineWidth(1);
    histo1->SetMarkerSize(0.3);
    histo1->SetMarkerStyle(20);
    histo1->SetTitleOffset(1.4,"xy");		
    histo1->SetTitleSize(0.05,"xy");		
    histo1->GetYaxis()->SetLabelSize(0.05);
    histo1->GetXaxis()->SetLabelSize(0.05);
    histo1->GetXaxis()->SetNdivisions(507,kTRUE);
    histo1->DrawCopy("p,e1");
    if(Title.Length() > 0){
        histo1->SetTitle("");
        TLatex *alice = new TLatex(0.1,0.95,Form("%s",Title.Data())); // Bo: this was 
        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextSize(0.062);
        alice->Draw();		
    }
}

void DrawGammaHistoRatioLowerPanel(     TH1* histo1, 
                                        TString yTitle,
                                        Float_t yMin, 
                                        Float_t yMax, 
                                        Int_t nDivisionsY=505,
                                        Double_t yLabelSize = 0.08, 
                                        Double_t yTitleSize= 0.1,
                                        Double_t yTitleOffset = 0.42,
                                        Double_t xLabelSize = 0.08, 
                                        Double_t xTitleSize= 0.11, 
                                        Double_t xTitleOffset = 1.) {
    cout << "here" << endl;
    histo1->SetYTitle(yTitle.Data());	
    cout << "here" << endl;
    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);

    histo1->GetYaxis()->SetRangeUser(yMin,yMax);
    cout << "here" << endl;
    histo1->GetYaxis()->SetNdivisions(nDivisionsY);
    cout << "here" << endl;
    histo1->GetYaxis()->SetLabelSize(yLabelSize);
    cout << "here" << endl;
    histo1->GetYaxis()->SetTitleSize(yTitleSize);	
    cout << "here" << endl;
    histo1->GetYaxis()->SetTitleOffset(yTitleOffset);
    cout << "here" << endl;
    histo1->GetYaxis()->SetDecimals();
    cout << "here" << endl;
    histo1->GetXaxis()->SetLabelSize(xLabelSize);	
    cout << "here" << endl;
    histo1->GetXaxis()->SetTitleSize(xTitleSize);	
    cout << "here" << endl;
    histo1->GetXaxis()->SetTitleOffset(xTitleOffset);	
    cout << "here" << endl;
}


void DrawGammaHistoWithTitle(   TH1* histo1, 
                                TH1* histo2, 
                                TString Title,
                                TString XTitle,
                                TString YTitle,
                                Float_t xMin, 
                                Float_t xMax, 
                                Float_t yMin) {
	
    histo1->GetXaxis()->SetRangeUser(xMin, xMax);
    histo1->GetYaxis()->SetRangeUser(yMin, 2.5*histo1->GetMaximum());
    if(Title.Length() > 0){
        histo1->SetTitle("");
        TLatex *alice = new TLatex(0.1,0.95,Form("%s",Title.Data())); // Bo: this was 
        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextSize(0.062);
        alice->Draw();		
    }
    if(XTitle.Length() > 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.Length() > 0){
        histo1->SetYTitle(YTitle.Data());
    }
    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);

    histo1->GetYaxis()->SetLabelSize(0.02);
    histo1->GetYaxis()->SetTitleSize(0.025);
    histo1->GetYaxis()->SetDecimals();
    histo1->GetXaxis()->SetTitleSize(0.025);
    histo1->GetXaxis()->SetLabelSize(0.02);
    histo1->SetMarkerStyle(20);
    histo1->SetMarkerColor(1);
    histo1->SetLineColor(1);
    histo1->SetLineWidth(1.);
    histo1->SetMarkerSize(0.2);
    histo1->SetTitleOffset(1.4,"xy");		
    histo1->SetTitleSize(0.05,"xy");		
    histo1->GetYaxis()->SetLabelSize(0.05);
    histo1->GetXaxis()->SetLabelSize(0.05);
    histo1->GetXaxis()->SetNdivisions(507,kTRUE);
    histo1->DrawCopy("hist");
    if(Title.Length() > 0){
        histo1->SetTitle("");
        TLatex *alice = new TLatex(0.1,0.95,Form("%s",Title.Data())); // Bo: this was 
        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextSize(0.062);
        alice->Draw();		
    }
    histo2->SetLineStyle(1);		
    histo2->SetLineColor(2);
    histo2->SetMarkerColor(2);
    histo2->SetMarkerSize(0.3);
    histo2->SetMarkerStyle(20);
    histo2->SetLineWidth(1);
    histo2->DrawCopy("p,e1,same");
}

void DrawFitResultsTwoSpecies(  TH1* histo1, 
                                TH1* histo2, 
                                TH1* histo3, 
                                TH1* histo4,
                                TString Title, 
                                TString XTitle, 
                                TString YTitle, 
                                TString legendString1, 
                                TString legendString2,
                                Bool_t YRange, 
                                Float_t YMin, 
                                Float_t YMax,  
                                Bool_t XRange, 
                                Float_t XMin,
                                Float_t XMax) {
    if (YRange && XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);	
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (YRange && !XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);
    }
    if (XRange && !YRange){
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);
    }
    if(Title.CompareTo("") != 0){
        histo1->SetTitle(Title.Data());
    }else { histo1->SetTitle("");}

    if(XTitle.CompareTo("") != 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.CompareTo("") != 0){
        histo1->SetYTitle(YTitle.Data());
    }
    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);

    histo1->GetYaxis()->SetTitleSize(0.04);	
    histo1->GetYaxis()->SetLabelSize(0.03);
    histo1->GetYaxis()->SetDecimals();
    histo1->GetYaxis()->SetTitleOffset(1.3);
    histo1->GetXaxis()->SetTitleOffset(1.1);
    histo1->GetXaxis()->SetTitleSize(0.04);
    histo1->GetXaxis()->SetLabelSize(0.03);	
    histo1->SetLineStyle(1);		
    histo1->SetLineColor(kRed);
    histo1->SetMarkerColor(kRed);
    histo1->SetMarkerSize(0.7);
    histo1->SetMarkerStyle(20);
    histo1->DrawCopy("e1,p");
    histo2->SetLineStyle(1);		
    histo2->SetLineColor(kRed-7);
    histo2->SetMarkerColor(kRed-7);
    histo2->SetMarkerSize(0.7);
    histo2->SetMarkerStyle(24);
    histo2->DrawCopy("e1,p,same");
    histo3->SetLineStyle(1);		
    histo3->SetLineColor(kBlue);
    histo3->SetMarkerColor(kBlue);
    histo3->SetMarkerSize(0.5);
    histo3->SetMarkerStyle(21);
    histo3->DrawCopy("e1,p,same");
    histo4->SetLineStyle(1);		
    histo4->SetLineColor(kBlue-7);
    histo4->SetMarkerColor(kBlue-7);
    histo4->SetMarkerSize(0.5);
    histo4->SetMarkerStyle(25);
    histo4->DrawCopy("e1,p,same");


    TLegend* leg2 = new TLegend(0.7,0.82,0.97,0.97);
    leg2->SetTextSize(0.035);			
    leg2->SetFillColor(0);
    leg2->AddEntry(histo1,Form("Data %s",legendString1.Data()),"pe");
    leg2->AddEntry(histo2,Form("MC %s",legendString1.Data()),"pe");
    leg2->AddEntry(histo3,Form("Data %s",legendString2.Data()),"pe");
    leg2->AddEntry(histo4,Form("MC %s",legendString2.Data()),"pe");
    leg2->Draw("same");

}


/*DrawAutoGammaHisto2D is a function for drawing a 2D-histogram of the gamma conversion group
* histo - histogramm which need to be drawn
* Title - histogram title
* XTitle - X- axis-title
* YTitle - Y-axis-title
* Input - Legend 
* YRange - if kTRUE will scale by YMin and YMay
* YMin  - Y minimum
* YMax - Y maximum
* XRange - if kTRUE will scale by XMin and XMax
* XMin - X minimum
* XMax - X maximum
*/
void DrawHistoCorrelationSurf2D(    TH2 *histo,  
                                    TString Title, 
                                    TString XTitle, 
                                    TString YTitle, 
                                    TString ZTitle, 
                                    Bool_t YRange, 
                                    Float_t YMin, 
                                    Float_t YMax, 
                                    Bool_t XRange, 
                                    Float_t XMin, 
                                    Float_t XMax, 
                                    TString optionDraw = "SURF2Z",
                                    Size_t labelSizeX = 0.035, 
                                    Size_t titleSizeX = 0.043, 
                                    Double_t titleOffsetX = 1.2,
                                    Size_t labelSizeY = 0.035, 
                                    Size_t titleSizeY = 0.043, 
                                    Double_t titleOffsetY = 1.4
                                ) {
    histo->SetTitle(Title.Data());	
    if (YRange && XRange){
        histo->GetYaxis()->SetRangeUser(YMin, YMax);	
        histo->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if ( !YRange && XRange){
        histo->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    
    if (YRange && !XRange){
        histo->GetYaxis()->SetRangeUser(YMin, YMax);
    }
    
    if(XTitle.CompareTo("") != 0){
        histo->SetXTitle(XTitle.Data());
    }
    if(YTitle.CompareTo("") != 0){
        histo->SetYTitle(YTitle.Data());
    }
    if(ZTitle.CompareTo("") != 0){
        histo->SetZTitle(ZTitle.Data());
    }
    histo->GetYaxis()->SetLabelFont(42);
    histo->GetXaxis()->SetLabelFont(42);
    histo->GetYaxis()->SetTitleFont(62);
    histo->GetXaxis()->SetTitleFont(62);

    histo->GetYaxis()->SetLabelSize(labelSizeY);
    histo->GetYaxis()->SetTitleSize(titleSizeY);	
    histo->GetYaxis()->SetTitleOffset(titleOffsetY);
    histo->GetYaxis()->SetDecimals();
    
    histo->GetXaxis()->SetLabelSize(labelSizeX);
    histo->GetXaxis()->SetTitleSize(titleSizeX);	
    histo->GetXaxis()->SetTitleOffset(titleOffsetX);
    
    histo->GetZaxis()->SetTitleOffset(2.);
    histo->DrawCopy(optionDraw.Data());
}

/* DrawAutoGammaMesonHistos is function used for styling the histograms of the gamma conversion group for two histos and standart settings
* histo1 - first histogram
* Title - histogram title
* XTitle - X-axis title
* YTitle - Y-axis title
* YRangeMax 	= kTRUE will scale by Maximum and Minimum Range in Y
*YMaxFactor - will MaximumY by this factor if YRangeMay = kTRUE 
*YMinimum - this will be used if YRangeMax is set
*YRange  	= kTRUE will Cut y-axis by YMin and YMax 
- will be set to kFAlSE if YRangeMax is set
*YMin - minimum Y
*YMax - maximum Y
*XRange 	= kTRUE will Cut x-axis by XMin and XMax
*XMin - minimum Y
*XMax - maximum Y
*/ 

void DrawCorrelationHisto1D(    TH1* histo1, 
                                TString Title, 
                                TString XTitle, 
                                TString YTitle, 
                                Bool_t YRangeMax, 
                                Float_t YMaxFactor, 
                                Float_t YMinimum,
                                Bool_t YRange, 
                                Float_t YMin,
                                Float_t YMax, 
                                Bool_t XRange, 
                                Float_t XMin, 
                                Float_t XMax) {
    if (YRangeMax && !XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        Double_t minRangeR = histo1->GetBinContent(histo1->GetMinimumBin())/10.;
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
    }
    if (YRangeMax && XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        Double_t minRangeR = histo1->GetBinContent(histo1->GetMinimumBin())/10.;
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);	
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (YRange && XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);	
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (!YRangeMax && !YRange && XRange){
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);	
    }
    if (YRange && !XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);
    }

    if(Title.CompareTo("") != 0){
        histo1->SetTitle(Title.Data());
    }else{
        histo1->SetTitle("");
    }
    if(XTitle.CompareTo("") != 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.CompareTo("") != 0){
        histo1->SetYTitle(YTitle.Data());
    }
    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);

    histo1->GetYaxis()->SetLabelSize(0.03);
    histo1->GetYaxis()->SetTitleSize(0.04);	
    histo1->GetYaxis()->SetDecimals();
    histo1->GetYaxis()->SetTitleOffset(1.2);
    histo1->GetXaxis()->SetTitleSize(0.04);	
    histo1->GetXaxis()->SetLabelSize(0.03);

    histo1->DrawCopy("e1,p");

}

void ReturnCorrectValuesForCanvasScaling(   Int_t sizeX, 
                                            Int_t sizeY, 
                                            Int_t nCols, 
                                            Int_t nRows, 
                                            Double_t leftMargin, 
                                            Double_t rightMargin, 
                                            Double_t upperMargin, 
                                            Double_t lowerMargin, 
                                            Double_t* arrayBoundariesX, 
                                            Double_t* arrayBoundariesY,
                                            Double_t* relativeMarginsX, 
                                            Double_t* relativeMarginsY){
    Int_t realsizeX             = sizeX- (Int_t)(sizeX*leftMargin)- (Int_t)(sizeX*rightMargin);
    Int_t realsizeY             = sizeY- (Int_t)(sizeY*upperMargin)- (Int_t)(sizeY*lowerMargin);
    
    Int_t nPixelsLeftColumn     = (Int_t)(sizeX*leftMargin);
    Int_t nPixelsRightColumn    = (Int_t)(sizeX*rightMargin);
    Int_t nPixelsUpperColumn    = (Int_t)(sizeY*upperMargin);
    Int_t nPixelsLowerColumn    = (Int_t)(sizeY*lowerMargin);
    
    Int_t nPixelsSinglePlotX    = (Int_t) (realsizeX/nCols);
    Int_t nPixelsSinglePlotY    = (Int_t) (realsizeY/nRows);
    cout << realsizeX << "\t" << nPixelsSinglePlotX << endl;
    cout << realsizeY << "\t" << nPixelsSinglePlotY << endl;
    
    cout << nPixelsLeftColumn << "\t" << nPixelsRightColumn  << "\t" << nPixelsLowerColumn << "\t" << nPixelsUpperColumn << endl;
    
    Int_t pixel = 0;
    cout << "boundaries X" << endl;
    for (Int_t i = 0; i < nCols+1; i++){
        if (i == 0){
            arrayBoundariesX[i] = 0.;
            pixel = pixel+nPixelsLeftColumn+nPixelsSinglePlotX;
        } else if (i == nCols){
            arrayBoundariesX[i] = 1.;
            pixel = pixel+nPixelsRightColumn;
        } else {
            arrayBoundariesX[i] = (Double_t)pixel/sizeX;
            pixel = pixel+nPixelsSinglePlotX;
        }   
        cout << i << "\t" << arrayBoundariesX[i] << "\t" << pixel<<endl;
    }   
    
    cout << "boundaries Y" << endl;
    pixel = sizeY;
    for (Int_t i = 0; i < nRows+1; i++){
        if (i == 0){
            arrayBoundariesY[i] = 1.;
            pixel = pixel-nPixelsUpperColumn-nPixelsSinglePlotY;
        } else if (i == nRows){
            arrayBoundariesY[i] = 0.;
            pixel = pixel-nPixelsLowerColumn;
        } else {
            arrayBoundariesY[i] = (Double_t)pixel/sizeY;
            pixel = pixel-nPixelsSinglePlotY;
        }   
        cout << i << "\t" << arrayBoundariesY[i] <<"\t" << pixel<<endl;
    }   
    
    relativeMarginsX[0]         = (Double_t)nPixelsLeftColumn/(nPixelsLeftColumn+nPixelsSinglePlotX);
    relativeMarginsX[1]         = 0;
    relativeMarginsX[2]         = (Double_t)nPixelsRightColumn/(nPixelsRightColumn+nPixelsSinglePlotX);;
    
    relativeMarginsY[0]         = (Double_t)nPixelsUpperColumn/(nPixelsUpperColumn+nPixelsSinglePlotY);
    relativeMarginsY[1]         = 0;
    relativeMarginsY[2]         = (Double_t)nPixelsLowerColumn/(nPixelsLowerColumn+nPixelsSinglePlotY);;
    
    return;
}   

void ReturnCorrectValuesTextSize(   TPad * pad, 
                                    Double_t &textsizeLabels, 
                                    Double_t &textsizeFac, 
                                    Int_t textSizeLabelsPixel, 
                                    Double_t margin){
    
    textsizeLabels = 0;
    textsizeFac = 0;
    if (pad->XtoPixel(pad->GetX2()) < pad->YtoPixel(pad->GetY1())){
        textsizeLabels = (Double_t)textSizeLabelsPixel/pad->XtoPixel(pad->GetX2()) ;
        textsizeFac = (Double_t)1./pad->XtoPixel(pad->GetX2()) ;
    } else {
        textsizeLabels = (Double_t)textSizeLabelsPixel/pad->YtoPixel(pad->GetY1());
        textsizeFac = (Double_t)1./pad->YtoPixel(pad->GetY1());
    }
    cout << textsizeLabels << endl;
    cout << textsizeFac << endl;
    
    return;
    
}
