#include <Riostream.h>
#include "TMath.h"
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TPaveLabel.h>
#include <TSystem.h>
#include <TFrame.h>
#include <TStyle.h>
#include <TString.h>
#include "TGaxis.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"

extern TRandom*	gRandom;
extern TBenchmark*	gBenchmark;
extern TSystem*	gSystem;
extern TMinuit*  	gMinuit;
const Double_t kMean=0.136 ; //Approximate peak position to facilitate error estimate


//_______________________ Plotting Invariant mass with BG and subtraction in a single p_t bin __________________________________
void PlotExampleInvMass(  TH1D* histoInvMassSignalWithBG,
                                TH1D* histoInvMassSubtracted,
                                TH1D* histoInvMassBG,
                                TF1* fitSignal,
                                Int_t exampleBin,
                                Double_t* fRangeBinsPt,
                                TString outputDir,
                                TString suffix,
                                TString triggerStr,
                                Double_t* fPlottingRangeMeson,
                                Float_t* pictDrawingCoordinatesDummy,
                                TString dateDummy,
                                TString fMesonType,
                                TString fSimulation,
                                TString fPlottingType,
                                TString fCollisionSystemDummy,
                                TString decayChannel,
                                TString detectionChannel                = "",
                                Double_t scaleFacSignal                 = 1.0,
                                Int_t detMode                           = 0
                             ){

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
    if(fMesonType.CompareTo("Pi0") == 0){
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
    if(fMesonType.CompareTo("Pi0") == 0){
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
    if(fMesonType.CompareTo("Pi0") == 0){
        labelInvMassPtRange = new TLatex(0.95,0.9, Form("#pi^{0}: %3.1f GeV/#it{c} < #it{p}_{T} < %3.1f GeV/#it{c}",startPt,endPt));
        fitPi0InvMassSig->SetRange(0,0.255);
        fitPi0InvMassSigRemBG->SetRange(0,0.255);
    } else {
        labelInvMassPtRange = new TLatex(0.95,0.9, Form("#eta: %3.1f GeV/#it{c} < #it{p}_{T} < %3.1f GeV/#it{c}",startPt,endPt));
        fitPi0InvMassSig->SetRange(0.35,0.75);
        fitPi0InvMassSigRemBG->SetRange(0.35,0.695);
    }
    // Set fit colors
    fitPi0InvMassSig->SetNpx(10000);
    fitPi0InvMassSig->SetLineColor(fitColorInvMassSG);
    fitPi0InvMassSig->SetLineWidth(3);
    fitPi0InvMassSigRemBG->SetNpx(10000);
    fitPi0InvMassSigRemBG->SetLineColor(fitColorInvMassSG);
    fitPi0InvMassSigRemBG->SetLineWidth(3);

    TH1D* histoFit  = (TH1D*)fitPi0InvMassSig->GetHistogram();
    histoFit->SetTitle("");
    histoFit->Scale(scaleFacSignal);
    TH1D* histoFitWBG  = (TH1D*)fitPi0InvMassSigRemBG->GetHistogram();
    histoFitWBG->SetTitle("");
    histoFitWBG->Scale(scaleFacSignal);

    canvasInvMassSamplePlot->cd();
    if(fMesonType.CompareTo("Pi0") == 0){
        histo1DInvMassDummy->GetYaxis()->SetRangeUser(histoPi0InvMassSigRemBGSub->GetMinimum(),1.1*histoPi0InvMassSigPlusBG->GetMaximum());
    } else {
      histo1DInvMassDummy->GetYaxis()->SetRangeUser(1.2*histoPi0InvMassSigRemBGSub->GetMinimum(),1.1*histoPi0InvMassSigPlusBG->GetMaximum());
    }
    histo1DInvMassDummy->GetXaxis()->SetRangeUser(fPlottingRangeMeson[0],fPlottingRangeMeson[1]);
    histo1DInvMassDummy->Draw("AXIS");

    DrawGammaSetMarker(histoPi0InvMassBGTot, /*markerStyleInvMassMBG*/0, /*markerSizeInvMassMBG*/0, markerColorInvMassMBG, markerColorInvMassMBG);
    histoPi0InvMassBGTot->SetLineWidth(2);
    histoPi0InvMassBGTot->SetLineStyle(2);
    histoPi0InvMassBGTot->Draw("hist,e,same");
    DrawGammaSetMarker(histoPi0InvMassSigPlusBG, markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
    histoPi0InvMassSigPlusBG->SetLineWidth(1);
    histoPi0InvMassSigPlusBG->Draw("hist,e,same");

    Int_t nLegendLines      = 5; //4
    histoPi0InvMassBGTot->GetXaxis()->SetRangeUser(fPlottingRangeMeson[0],fPlottingRangeMeson[1]);
    histoPi0InvMassSigRemBGSub->GetXaxis()->SetRangeUser(fPlottingRangeMeson[0],fPlottingRangeMeson[1]);
    if (scaleFacSignal == 1.0){
        fitPi0InvMassSig->Draw("same");
        DrawGammaSetMarker(histoPi0InvMassSigRemBGSub, markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
        histoPi0InvMassSigRemBGSub->Draw("same");
    } else {
        histoFit->SetLineColor(fitColorInvMassSG);
        histoFit->SetLineWidth(4);
        histoFit->Draw("same");
        nLegendLines++;

        DrawGammaSetMarker(histoPi0InvMassSigRemBGSub, markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
        histoPi0InvMassSigRemBGSub->Draw("same");
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

    TLatex *labelInvMassReco  = new TLatex(0.135,0.9-3*0.9*textsizeLabelsPP, detectionChannel.Data());
    SetStyleTLatex( labelInvMassReco, 0.85*textSizeLabelsPixel,4);
    labelInvMassReco->SetTextFont(43);
    labelInvMassReco->Draw();

    SetStyleTLatex( labelInvMassPtRange, 0.85*textSizeLabelsPixel,4);
    labelInvMassPtRange->SetTextAlign(31);
    labelInvMassPtRange->SetTextFont(43);
    labelInvMassPtRange->Draw();

    TLegend* legendInvMass  = GetAndSetLegend2(0.66, 0.87-nLegendLines*0.85*textsizeLabelsPP, 0.9, 0.87, 0.85*textSizeLabelsPixel);
    legendInvMass->SetMargin(0.25);
    legendInvMass->AddEntry(histoPi0InvMassSigPlusBG,"Raw real events","l");
    legendInvMass->AddEntry(histoPi0InvMassBGTot,"Mixed event BG","l");
//     legendInvMass->AddEntry((TObject*)0,"+ rem. BG","");
//     legendInvMass->AddEntry(histoPi0InvMassSigRemBGSub,"BG subtracted","p");
    legendInvMass->AddEntry(histoPi0InvMassSigRemBGSub,"Signal after BG","p");
    legendInvMass->AddEntry((TObject*)0,"subtraction","");
    if (scaleFacSignal != 1.0){
        legendInvMass->AddEntry((TObject*)0,Form("scaled by %2.f",scaleFacSignal),"");
    }
    legendInvMass->AddEntry(fitPi0InvMassSig, "Fit","l");
    legendInvMass->Draw();

    histo1DInvMassDummy->Draw("AXIS,same");
    canvasInvMassSamplePlot->SaveAs(Form("%s/%s_%s_InvMassBin_%s_%s.%s",outputDir.Data(),fMesonType.Data(),fSimulation.Data(),detectionChannel.Data(),fPlottingType.Data(), suffix.Data()));

}



//_______________________ Plotting Invariant mass with BG and subtraction in a single p_t bin __________________________________
void PlotExampleInvMassEMCal(  TH1D* histoInvMassSignalWithBG,
                                TF1* fitSignal,
                                TF1* fitBG,
                                Int_t exampleBin,
                                Double_t* fRangeBinsPt,
                                TString outputDir,
                                TString suffix,
                                TString triggerStr,
                                Double_t* fPlottingRangeMeson,
                                Float_t* pictDrawingCoordinatesDummy,
                                TString dateDummy,
                                TString fMesonType,
                                TString fSimulation,
                                TString fPlottingType,
                                TString fCollisionSystemDummy,
                                TString decayChannel,
                                TString detectionChannel                = "",
                                Double_t scaleFacSignal                 = 1.0,
                                Int_t detMode                           = 0
                             ){

    TH1D* histoPi0InvMassSigPlusBG;
    TH1D* histoPi0InvMassSigRemBG;
    TH1D* histoPi0InvMassSigRemBGSub;
    TH1D* histoPi0InvMassRemBG;
    TF1* fitPi0InvMassSig;
    TF1* fitPi0InvMassSigRemBG;
    TF1* fitPi0InvMassBG;

    histoPi0InvMassSigPlusBG         = (TH1D*)histoInvMassSignalWithBG->Clone("InvMassSigPlusBG_PtBin07");
    histoPi0InvMassSigPlusBG->SetStats(0);
    fitPi0InvMassSig                 = (TF1*)fitSignal->Clone("FitInvMassSig_PtBin07");
    fitPi0InvMassBG                  = (TF1*)fitBG->Clone("FitInvMassBG_PtBin07");

    histoPi0InvMassSigRemBGSub = (TH1D*)histoPi0InvMassSigPlusBG->Clone("histoPi0InvMassSigRemBGSub");
    histoPi0InvMassSigRemBGSub->Add(fitPi0InvMassBG,-1.);

    TH1D* histoPi0InvMassBGTot = (TH1D*) fitPi0InvMassBG->GetHistogram();
    DrawGammaSetMarker(histoPi0InvMassBGTot, /*markerStyleInvMassMBG*/0, /*markerSizeInvMassMBG*/0, kGray+2, kGray+2);
    histoPi0InvMassBGTot->SetLineWidth(2);
    histoPi0InvMassBGTot->SetLineStyle(2);

    if(fMesonType.CompareTo("Pi0") == 0 ){
      fitPi0InvMassSigRemBG  = new TF1("gaus","gaus",0.08,0.22);
      histoPi0InvMassSigRemBGSub->Fit(fitPi0InvMassSigRemBG,"","",0.08,0.22);
    } else {
      fitPi0InvMassSigRemBG  = new TF1("gaus","gaus",0.4,0.7);
      histoPi0InvMassSigRemBGSub->Fit(fitPi0InvMassSigRemBG,"","",0.4,0.7);
    }
//     histoPi0InvMassSigPlusBG->SetMinimum(-20.);


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
    if(fMesonType.CompareTo("Pi0") == 0){
        labelInvMassPtRange = new TLatex(0.95,0.9, Form("#pi^{0}: %3.1f GeV/#it{c} < #it{p}_{T} < %3.1f GeV/#it{c}",startPt,endPt));
        fitSignal->SetRange(0,0.255);
    } else {
        labelInvMassPtRange = new TLatex(0.95,0.9, Form("#eta: %3.1f GeV/#it{c} < #it{p}_{T} < %3.1f GeV/#it{c}",startPt,endPt));
        fitSignal->SetRange(0.35,0.695);
    }
    // Set fit colors
    fitSignal->SetNpx(10000);
    fitSignal->SetLineColor(fitColorInvMassSG);
    fitSignal->SetLineWidth(3);

    TH1D* histoFit  = (TH1D*)fitPi0InvMassSigRemBG->GetHistogram();
    histoFit->SetTitle("");
    histoFit->Scale(scaleFacSignal);

    canvasInvMassSamplePlot->cd();
    histo1DInvMassDummy->GetYaxis()->SetRangeUser(1.45*histoPi0InvMassSigRemBGSub->GetMinimum(),1.1*histoPi0InvMassSigPlusBG->GetMaximum());
    histo1DInvMassDummy->GetXaxis()->SetRangeUser(fPlottingRangeMeson[0],fPlottingRangeMeson[1]);
    histo1DInvMassDummy->Draw("AXIS");

    histoPi0InvMassBGTot->GetXaxis()->SetRangeUser(fPlottingRangeMeson[0],fPlottingRangeMeson[1]);
    histoPi0InvMassBGTot->Draw("hist,same");
    DrawGammaSetMarker(histoPi0InvMassSigPlusBG, markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
    histoPi0InvMassSigPlusBG->SetLineWidth(1);
    histoPi0InvMassSigPlusBG->Draw("hist,e,same");

    Int_t nLegendLines      = 4;//5;
    histoPi0InvMassSigRemBGSub->GetXaxis()->SetRangeUser(fPlottingRangeMeson[0],fPlottingRangeMeson[1]);
    if (scaleFacSignal == 1.0){
        DrawGammaSetMarker(histoPi0InvMassSigRemBGSub, markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
        histoPi0InvMassSigRemBGSub->Draw("p,same");
        fitPi0InvMassSigRemBG->SetLineColor(fitColorInvMassSG);
        fitPi0InvMassSigRemBG->Draw("same");
    } else {
        DrawGammaSetMarker(histoPi0InvMassSigRemBGSub, markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
        histoPi0InvMassSigRemBGSub->Draw("p,same");

        histoFit->SetLineColor(fitColorInvMassSG);
        histoFit->SetLineWidth(4);
        histoFit->Draw("same");
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

    TLatex *labelInvMassReco  = new TLatex(0.135,0.9-3*0.9*textsizeLabelsPP, detectionChannel.Data());
    SetStyleTLatex( labelInvMassReco, 0.85*textSizeLabelsPixel,4);
    labelInvMassReco->SetTextFont(43);
    labelInvMassReco->Draw();

    SetStyleTLatex( labelInvMassPtRange, 0.85*textSizeLabelsPixel,4);
    labelInvMassPtRange->SetTextAlign(31);
    labelInvMassPtRange->SetTextFont(43);
    labelInvMassPtRange->Draw();

    TLegend* legendInvMass  = GetAndSetLegend2(0.67, 0.87-nLegendLines*textsizeLabelsPP, 0.9, 0.87, 0.85*textSizeLabelsPixel);
    legendInvMass->SetMargin(0.25);
    legendInvMass->AddEntry(histoPi0InvMassSigPlusBG,"Raw real events","l");
    legendInvMass->AddEntry(histoPi0InvMassBGTot,"Mixed event BG","l");
//     legendInvMass->AddEntry((TObject*)0,"+ rem. BG","");
//     legendInvMass->AddEntry(histoPi0InvMassSigRemBGSub,"BG subtracted","p");
    legendInvMass->AddEntry(histoPi0InvMassSigRemBGSub,"Signal after BG","p");
    legendInvMass->AddEntry((TObject*)0,"subtraction","");
    if (scaleFacSignal != 1.0){
        legendInvMass->AddEntry((TObject*)0,Form("scaled by %2.1f",scaleFacSignal),"");
    }
    legendInvMass->AddEntry(fitPi0InvMassSigRemBG, "Fit","l");
    legendInvMass->Draw();

    histo1DInvMassDummy->Draw("AXIS,same");
    canvasInvMassSamplePlot->SaveAs(Form("%s/%s_%s_InvMassBin_%s_%s.%s",outputDir.Data(),fMesonType.Data(),fSimulation.Data(),detectionChannel.Data(),fPlottingType.Data(), suffix.Data()));

}

void PlotInvMassEMCalAndPCMTogether(Int_t binPbPb = 4){


	gROOT->Reset();
	gROOT->SetStyle("Plain");
	StyleSettingsThesis();
	SetPlotStyle();

	TString dateForOutput = ReturnDateStringForOutput();
	TString outputDir = Form("pdf/%s/InvMassBinPlots",dateForOutput.Data());
	gSystem->Exec("mkdir -p "+outputDir);

    Double_t fBinsPi0HIPtLHC11h[27] = { 0.0, 0.4, 0.6, 0.8, 1.0,
                                        1.2, 1.4, 1.6, 1.8, 2.0,
                                        2.2, 2.4, 2.6, 3.0, 3.5,
                                        4.0, 5.0, 6.0, 8.0, 10.0,
                                        12.0, 14.0, 16.0, 18.0, 20.0,
                                        25.0, 30.0};
    Double_t fBinsEtaHIPtLHC11hLessBins[13]  = { 0.0, 0.5, 1.0, 1.5, 2.0, 3.0,
                                                  4., 6.0, 8.0, 10., 12.,
                                                  14., 19.};


    Double_t fMesonMassPlotRangePi0[2]     = {0.09,0.215};
    Double_t fMesonMassPlotRangeEta[2]     = {0.4,0.745};
    Float_t  pictDrawingCoordinatesFWHM[9]  = {0.6, 0.8, 0.30, 0.04, 0.15,0.7, 0.1, 0.035,0};

    Double_t scaleFactorPi0PCM = 10;
    Double_t scaleFactorEtaPCM = 40;
    Double_t noScaling = 1;
    Double_t ptMinEMCal = 12;
    Double_t ptMaxEMCal = 14;

    TString fDetectionProcessPCM = "PCM";
    TString fDetectionProcessEMCal = "EMCal";
    TString fCollisionSystem = "Pb#font[122]{-}Pb, #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";
    TString fDecayChannelPi0 = "#pi^{0} #rightarrow #gamma#gamma";
    TString fDecayChannelEta = "#eta #rightarrow #gamma#gamma";

//===============================================================================================

    TString cutselection0010 = "50100013_00200009247602008850404000_0652501500000000";
    TFile *filePi0PCMPbPb0010           = new TFile(Form("%s/PbPb_2.76TeV/Pi0_data_GammaConvV1WithoutCorrection_%s.root",cutselection0010.Data(),cutselection0010.Data()));
	TH1D* histoPi0PCMSignalPlusBG0010   = (TH1D*)filePi0PCMPbPb0010->Get(Form("Mapping_GG_InvMass_in_Pt_Bin%02d",binPbPb));
    TH1D* histoPi0PCMBG0010             = (TH1D*)filePi0PCMPbPb0010->Get(Form("Mapping_BckNorm_InvMass_in_Pt_Bin%02d",binPbPb));
	TH1D* histoPi0PCMSignal0010         = (TH1D*)filePi0PCMPbPb0010->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d",binPbPb));
	TF1* fitPi0PCMSignal0010            = (TF1*)filePi0PCMPbPb0010->Get(Form("Signal_InvMassFit_in_Pt_Bin%02d",binPbPb));

    PlotExampleInvMass(histoPi0PCMSignalPlusBG0010, histoPi0PCMSignal0010, histoPi0PCMBG0010, fitPi0PCMSignal0010,
                       binPbPb,fBinsPi0HIPtLHC11h,outputDir.Data(),"pdf", "0-10%", fMesonMassPlotRangePi0, pictDrawingCoordinatesFWHM, dateForOutput.Data(), "Pi0", "Data","0010",
                       fCollisionSystem, fDecayChannelPi0, fDetectionProcessPCM, scaleFactorPi0PCM, 0);


	TFile * fileEtaPCMPbPb0010          = new TFile(Form("%s/PbPb_2.76TeV/Eta_data_GammaConvV1WithoutCorrection_%s.root",cutselection0010.Data(),cutselection0010.Data()));
	TH1D* histoEtaPCMSignalPlusBG0010   = (TH1D*)fileEtaPCMPbPb0010->Get(Form("Mapping_GG_InvMass_in_Pt_Bin%02d",binPbPb));
    TH1D* histoEtaPCMBG0010             = (TH1D*)fileEtaPCMPbPb0010->Get(Form("Mapping_BckNorm_InvMass_in_Pt_Bin%02d",binPbPb));
	TH1D* histoEtaPCMSignal0010         = (TH1D*)fileEtaPCMPbPb0010->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d",binPbPb));
	TF1* fitEtaPCMSignal0010            = (TF1*)fileEtaPCMPbPb0010->Get(Form("Signal_InvMassFit_in_Pt_Bin%02d",binPbPb));

    PlotExampleInvMass(histoEtaPCMSignalPlusBG0010, histoEtaPCMSignal0010, histoEtaPCMBG0010, fitEtaPCMSignal0010,
                       binPbPb, fBinsEtaHIPtLHC11hLessBins, outputDir.Data(),"pdf", "0-10%",fMesonMassPlotRangeEta, pictDrawingCoordinatesFWHM, dateForOutput.Data(), "Eta", "Data", "0010",
                       fCollisionSystem, fDecayChannelEta, fDetectionProcessPCM, scaleFactorEtaPCM, 0);


    TString cutselection2050 = "52500013_00200009247602008850404000_0652501500000000";
    TFile *filePi0PCMPbPb2050           = new TFile(Form("%s/PbPb_2.76TeV/Pi0_data_GammaConvV1WithoutCorrection_%s.root",cutselection2050.Data(),cutselection2050.Data()));
    TH1D* histoPi0PCMSignalPlusBG2050   = (TH1D*)filePi0PCMPbPb2050->Get(Form("Mapping_GG_InvMass_in_Pt_Bin%02d",binPbPb));
    TH1D* histoPi0PCMBG2050             = (TH1D*)filePi0PCMPbPb2050->Get(Form("Mapping_BckNorm_InvMass_in_Pt_Bin%02d",binPbPb));
    TH1D* histoPi0PCMSignal2050         = (TH1D*)filePi0PCMPbPb2050->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d",binPbPb));
    TF1* fitPi0PCMSignal2050            = (TF1*)filePi0PCMPbPb2050->Get(Form("Signal_InvMassFit_in_Pt_Bin%02d",binPbPb));

    PlotExampleInvMass(histoPi0PCMSignalPlusBG2050, histoPi0PCMSignal2050, histoPi0PCMBG2050, fitPi0PCMSignal2050,
                       binPbPb,fBinsPi0HIPtLHC11h,outputDir.Data(),"pdf", "20-50%",fMesonMassPlotRangePi0, pictDrawingCoordinatesFWHM, dateForOutput.Data(), "Pi0", "Data","2050",
                       fCollisionSystem, fDecayChannelPi0, fDetectionProcessPCM, noScaling, 0);


    TFile * fileEtaPCMPbPb2050          = new TFile(Form("%s/PbPb_2.76TeV/Eta_data_GammaConvV1WithoutCorrection_%s.root",cutselection2050.Data(),cutselection2050.Data()));
    TH1D* histoEtaPCMSignalPlusBG2050   = (TH1D*)fileEtaPCMPbPb2050->Get(Form("Mapping_GG_InvMass_in_Pt_Bin%02d",binPbPb));
    TH1D* histoEtaPCMBG2050             = (TH1D*)fileEtaPCMPbPb2050->Get(Form("Mapping_BckNorm_InvMass_in_Pt_Bin%02d",binPbPb));
    TH1D* histoEtaPCMSignal2050         = (TH1D*)fileEtaPCMPbPb2050->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d",binPbPb));
    TF1* fitEtaPCMSignal2050            = (TF1*)fileEtaPCMPbPb2050->Get(Form("Signal_InvMassFit_in_Pt_Bin%02d",binPbPb));

    PlotExampleInvMass(histoEtaPCMSignalPlusBG2050, histoEtaPCMSignal2050, histoEtaPCMBG2050, fitEtaPCMSignal2050,
                       binPbPb, fBinsEtaHIPtLHC11hLessBins, outputDir.Data(),"pdf","20-50%", fMesonMassPlotRangeEta, pictDrawingCoordinatesFWHM, dateForOutput.Data(), "Eta", "Data", "2050",
                       fCollisionSystem, fDecayChannelEta, fDetectionProcessPCM, noScaling, 0);

//===============================================================================================


    TFile* filePi0EMCalPbPb0010         = TFile::Open("/home/admin1/leardini/alicepietapaper/pbpbrootfiles/InvMassFiles/EMCAL_PbPb_LHC11h_InvMassPi0.root");
    TH1D* histoPi0EMCalSignalPlusBG0010 = (TH1D*)filePi0EMCalPbPb0010->Get("hDiff1D");
    TH1D* histoPi0EMCalSignal0010       = (TH1D*)histoPi0EMCalSignalPlusBG0010->Clone("histoEMCalSignalPi00010");
    TF1*  fitPi0EMCalSignal0010         = (TF1*)filePi0EMCalPbPb0010->Get("fitPolyPi");
    TF1*  fitPi0EMCalBack0010           = (TF1*)filePi0EMCalPbPb0010->Get("fitFullBackpi0");

    PlotExampleInvMassEMCal(histoPi0EMCalSignalPlusBG0010, fitPi0EMCalSignal0010, fitPi0EMCalBack0010,
                       20, fBinsPi0HIPtLHC11h, outputDir.Data(),"pdf", "0-10%",fMesonMassPlotRangePi0, pictDrawingCoordinatesFWHM, dateForOutput.Data(), "Pi0", "Data", "0010",
                       fCollisionSystem, fDecayChannelPi0, fDetectionProcessEMCal, noScaling, 0);


    TFile* fileEtaEMCalPbPb0010         = TFile::Open("/home/admin1/leardini/alicepietapaper/pbpbrootfiles/InvMassFiles/EMCAL_PbPb_LHC11h_InvMassEta0.root");
    TH1D* histoEtaEMCalSignalPlusBG0010 = (TH1D*)fileEtaEMCalPbPb0010->Get("hDiff1DEta__1");
    TH1D* histoEtaEMCalSignal0010       = (TH1D*)histoEtaEMCalSignalPlusBG0010->Clone("histoEMCalSignalEta0010");
    TF1*  fitEtaEMCalSignal0010         = (TF1*)fileEtaEMCalPbPb0010->Get("fitPolyEta");
    TF1*  fitEtaEMCalBack0010           = (TF1*)fileEtaEMCalPbPb0010->Get("fitFullBackEta");

    PlotExampleInvMassEMCal(histoEtaEMCalSignalPlusBG0010, fitEtaEMCalSignal0010, fitEtaEMCalBack0010,
                       10, fBinsEtaHIPtLHC11hLessBins, outputDir.Data(),"pdf", "0-10%",fMesonMassPlotRangeEta, pictDrawingCoordinatesFWHM, dateForOutput.Data(), "Eta", "Data", "0010",
                       fCollisionSystem, fDecayChannelEta, fDetectionProcessEMCal, noScaling, 0);


    TFile* filePi0EMCalPbPb2050         = TFile::Open("/home/admin1/leardini/alicepietapaper/pbpbrootfiles/InvMassFiles/EMCAL_PbPb_LHC11h_InvMassPi2.root");
    TH1D* histoPi0EMCalSignalPlusBG2050 = (TH1D*)filePi0EMCalPbPb2050->Get("hDiff1D");
    TH1D* histoPi0EMCalSignal2050       = (TH1D*)histoPi0EMCalSignalPlusBG2050->Clone("histoEMCalSignalPi02050");
    TF1*  fitPi0EMCalSignal2050         = (TF1*)filePi0EMCalPbPb2050->Get("fitPolyPi");
    TF1*  fitPi0EMCalBack2050           = (TF1*)filePi0EMCalPbPb2050->Get("fitFullBackpi0");

    PlotExampleInvMassEMCal(histoPi0EMCalSignalPlusBG2050, fitPi0EMCalSignal2050, fitPi0EMCalBack2050,
                       20, fBinsPi0HIPtLHC11h, outputDir.Data(),"pdf", "20-50%",fMesonMassPlotRangePi0, pictDrawingCoordinatesFWHM, dateForOutput.Data(), "Pi0", "Data", "2050",
                       fCollisionSystem, fDecayChannelPi0, fDetectionProcessEMCal, noScaling, 0);


    TFile* fileEtaEMCalPbPb2050         = TFile::Open("/home/admin1/leardini/alicepietapaper/pbpbrootfiles/InvMassFiles/EMCAL_PbPb_LHC11h_InvMassEta2.root");
    TH1D* histoEtaEMCalSignalPlusBG2050 = (TH1D*)fileEtaEMCalPbPb2050->Get("hDiff1DEta__1");
    TH1D* histoEtaEMCalSignal2050       = (TH1D*)histoEtaEMCalSignalPlusBG2050->Clone("histoEMCalSignalEta2050");
    TF1*  fitEtaEMCalSignal2050         = (TF1*)fileEtaEMCalPbPb2050->Get("fitPolyEta");
    TF1*  fitEtaEMCalBack2050           = (TF1*)fileEtaEMCalPbPb2050->Get("fitFullBackEta");

    PlotExampleInvMassEMCal(histoEtaEMCalSignalPlusBG2050, fitEtaEMCalSignal2050, fitEtaEMCalBack2050,
                       10, fBinsEtaHIPtLHC11hLessBins, outputDir.Data(),"pdf", "20-50%",fMesonMassPlotRangeEta, pictDrawingCoordinatesFWHM, dateForOutput.Data(), "Eta", "Data", "2050",
                       fCollisionSystem, fDecayChannelEta, fDetectionProcessEMCal, noScaling, 0);


/*


    TFile *fEMC2 = TFile::Open("/home/admin1/leardini/alicepietapaper/pbpbrootfiles/InvMassFiles/EMCAL_PbPb_LHC11h_InvMassEta2.root");
    TFile *fEMCPI0 = TFile::Open("/home/admin1/leardini/alicepietapaper/pbpbrootfiles/InvMassFiles/EMCAL_PbPb_LHC11h_InvMassPi0.root");
    TFile *fEMCPI2 = TFile::Open("/home/admin1/leardini/alicepietapaper/pbpbrootfiles/InvMassFiles/EMCAL_PbPb_LHC11h_InvMassPi2.root");

    TH1F *histoEMCalSignalPlusBG0010;
    TF1 *fitEMCalSignal0010;
    TF1 *FitBack0010;
    TH1F *histoEMCalSignal0010;

    TF1 * fitpeak0010  = new TF1("gaus","gaus",0.08,0.22);
    DrawGammaSetMarkerTF1(fitpeak0010,1,1,kBlue+1);

    histoEMCalSignalPlusBG0010     = (TH1F*)fEMCPI0->Get("hDiff1D"); //hMassPi0
    histoEMCalSignalPlusBG0010->SetStats(0);
    fitEMCalSignal0010        = (TF1*)fEMCPI0->Get("fitPolyPi");    //fFitMassPi0
    FitBack0010       = (TF1*)fEMCPI0->Get("fitFullBackpi0");       //fFitBackPi0
    histoEMCalSignal0010 = (TH1F*)histoEMCalSignalPlusBG0010->Clone("histoEMCalSignal0010");
    histoEMCalSignal0010->Add(FitBack0010,-1.);
    histoEMCalSignal0010->Fit(fitpeak0010,"","",0.08,0.22);
    histoEMCalSignalPlusBG0010->SetMinimum(-20.);
//     fitpeak0010 = (TF1*)histoEMCalSignal0010->GetFunction("gaus");



//===============================================================================================

    TH1F *histoEMCalSignal2050;
    TH1F *histoEMCalSignalPlusBG2050;
    TF1 *fitEMCalSignal2050;
    TF1 *FitBack2050;

    TF1 * fitpeak2050  = new TF1("gaus","gaus",0.08,0.22);
    DrawGammaSetMarkerTF1(fitpeak2050,1,1,kBlue+1);

    histoEMCalSignalPlusBG2050     = (TH1F*)fEMCPI2->Get("hDiff1D");
    histoEMCalSignalPlusBG2050->SetStats(0);
    fitEMCalSignal2050        = (TF1*)fEMCPI2->Get("fitPolyPi");
    FitBack2050       = (TF1*)fEMCPI2->Get("fitFullBackpi0");
    histoEMCalSignal2050 = (TH1F*)histoEMCalSignalPlusBG2050->Clone("histoEMCalSignal2050");
    histoEMCalSignal2050->Add(FitBack2050,-1.);
    histoEMCalSignal2050->Fit(fitpeak2050,"","",0.08,0.22);
    histoEMCalSignalPlusBG2050->SetMinimum(-20.);
//     TF1 *fitpeak2050 = (TF1*)histoEMCalSignal2050->GetFunction("gaus");

//===============================================================================================

//===============================================================================================

    TH1F *histoEMCalSignalPlusBGEta2050;
    TF1 *fitEMCalSignalEta2050;
    TF1 *FitBackEta2050;
    TH1F *histoEMCalSignalEta2050;

    TF1 * fitpeakEta2050  = new TF1("gaus","gaus",0.45,0.7);
    DrawGammaSetMarkerTF1(fitpeakEta2050,1,1,kBlue+1);

    histoEMCalSignalPlusBGEta2050     = (TH1F*)fEMC2->Get("hDiff1DEta__1");
    histoEMCalSignalPlusBGEta2050->SetStats(0);

    fitEMCalSignalEta2050     = (TF1*)fEMC2->Get("fitPolyEta");
    FitBackEta2050        = (TF1*)fEMC2->Get("fitFullBackEta");
    histoEMCalSignalEta2050 = (TH1F*)histoEMCalSignalPlusBGEta2050->Clone("histoEMCalSignalEta2050");
    histoEMCalSignalEta2050->Add(FitBackEta2050,-1.);
    histoEMCalSignalEta2050->Fit(fitpeakEta2050,"","",0.4,0.75);
    histoEMCalSignalPlusBGEta2050->SetMinimum(-20.);
//     TF1 *fitpeakEta2050 = (TF1*)histoEMCalSignalEta2050->GetFunction("gaus");

*/

}



