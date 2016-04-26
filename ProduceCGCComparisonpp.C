
/****************************************************************************************************************************
******      provided by Gamma Conversion Group, GA,                                        *****
******      Ana Marin, a.marin@gsi.de                                      *****
******      Friederike Bock, friederike.bock@cern.ch                                      *****
******      Annika Passfeld, annikapassfeld@uni-muenster.de  
******      Pedro Gonzalez, pedro.gonzalez.zamora@cern.ch 
******      29.04.2014
*****************************************************************************************************************************/

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
#include "TGraph.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h" 
#include "TGaxis.h"
#include "TMarker.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CombineMesonMeasurements.h"


extern TRandom*   gRandom;
extern TBenchmark*   gBenchmark;
extern TSystem*   gSystem;
extern TMinuit*   gMinuit;





TGraph* ScaleGraph (TGraph* graph, Double_t scaleFac){
   TGraph* dummyGraph = (TGraph*)graph->Clone(Form("%s_Scaled",graph->GetName()));
   Double_t * xValue = dummyGraph->GetX();
   Double_t * yValue = dummyGraph->GetY();
   
   Int_t nPoints = dummyGraph->GetN();
   for (Int_t i = 0; i < nPoints; i++){
      yValue[i] = yValue[i]*scaleFac;
   }
   TGraph* returnGraph = new TGraph(nPoints,xValue,yValue);
   return returnGraph;
}

void ProduceCGCComparisonpp(){
  //  gROOT->Reset(); 
  gROOT->SetStyle("Plain");
  
  StyleSettingsThesis();  
  SetPlotStyle();

  Double_t    xSection2760GeV =    62.8*1e-3;
  Double_t    xSection7000GeV =    73.2*1e-3;
  Double_t    recalcBarn =      1e12; //NLO in pbarn!!!!



  Int_t index = 0;
 //*---CGC- T. Lappi and H.M"antysaari  arxiv1309.6963 (kt-factorization)---*//

  TGraph * graphInvYieldCGC7000GeV = new TGraph(66);
  Double_t pt_cgc7000GeV[66];
  Double_t dndydpt_cgc7000GeV[66];
  ifstream file_cgc7000("ExternalInput/Theory/ALICECGCcalcInclusiveYieldPi0Lappi7000GeV.dat");

  if (file_cgc7000.is_open()){
    while(!file_cgc7000.eof()){
      //cout<<index<<endl;
      file_cgc7000 >> pt_cgc7000GeV[index] >> dndydpt_cgc7000GeV[index];
      graphInvYieldCGC7000GeV->SetPoint(index,pt_cgc7000GeV[index],
					dndydpt_cgc7000GeV[index]);
      index++;
    }
    file_cgc7000.close();
    index = 0;
  }
  //  graphInvYieldCGC7000GeV->Print();
  TGraph * graphInvCrossSecCGC7000GeV = ScaleGraph(graphInvYieldCGC7000GeV,(xSection7000GeV*recalcBarn));
  //  graphInvYieldCGC7000GeV->Print();


  TGraph * graphInvYieldCGC7000GeV_mvgamma = new TGraph(54);
  Double_t pt_cgc7000GeV_mvgamma[54];
  Double_t dndydpt_cgc7000GeV_mvgamma[54];
  ifstream file_cgc7000_mvgamma("ExternalInput/Theory/ALICECGCcalcInclusiveYieldPi0Lappi7000GeV_mvgamma.dat");

  
  if (file_cgc7000_mvgamma.is_open()){
    while(!file_cgc7000_mvgamma.eof()){
      //cout<<index<<endl;
      file_cgc7000_mvgamma >> pt_cgc7000GeV_mvgamma[index] >> dndydpt_cgc7000GeV_mvgamma[index];
      graphInvYieldCGC7000GeV_mvgamma->SetPoint(index,pt_cgc7000GeV_mvgamma[index],
					dndydpt_cgc7000GeV_mvgamma[index]);
      index++;
    }
    file_cgc7000_mvgamma.close();
    index = 0;
  }
  
  //  graphInvYieldCGC7000GeV_mvgamma->Print();
  TGraph * graphInvCrossSecCGC7000GeV_mvgamma = ScaleGraph(graphInvYieldCGC7000GeV_mvgamma,(xSection7000GeV*recalcBarn));


  TGraph * graphInvYieldCGC2760GeV = new TGraph(50);
  Double_t pt_cgc2760GeV[50];
  Double_t dndydpt_cgc2760GeV[50];
  ifstream file_cgc2760("ExternalInput/Theory/ALICECGCcalcInclusiveYieldPi0Lappi2760GeV.dat");

  if (file_cgc2760.is_open()){
    while(!file_cgc2760.eof()){ 
      // cout<<index<<endl;
      file_cgc2760 >> pt_cgc2760GeV[index] >> dndydpt_cgc2760GeV[index];
      graphInvYieldCGC2760GeV->SetPoint(index,pt_cgc2760GeV[index],
					dndydpt_cgc2760GeV[index]);
      index++;
    }
    file_cgc2760.close();
    index = 0;
  }


  // graphInvYieldCGC2760GeV->Print();
  TGraph * graphInvCrossSecCGC2760GeV = ScaleGraph(graphInvYieldCGC2760GeV,(xSection2760GeV*recalcBarn));
  TFile fileCGCGraphsPP("ExternalInput/CGCCompilationPP.root","RECREATE");
  graphInvYieldCGC2760GeV->Write("graphInvYieldCGC2760GeV");
  graphInvYieldCGC7000GeV->Write("graphInvYieldCGC7000GeV");
  graphInvCrossSecCGC2760GeV->Write("graphInvCrossSecCGC2760GeV");
  graphInvCrossSecCGC7000GeV->Write("graphInvCrossSecCGC7000GeV");
  graphInvCrossSecCGC7000GeV_mvgamma->Write("graphInvCrossSecCGC7000GeV_mvgamma");

  fileCGCGraphsPP.Close();


  TFile fileppCombined("CombinedResultsPaperX_18_Feb_2014.root");

  TGraphAsymmErrors* graphInvCrossSectionPi0Comb2760GeV = 
    (TGraphAsymmErrors*)fileppCombined.Get("graphInvCrossSectionPi0Comb2760GeV");
  TGraphAsymmErrors* graphInvCrossSectionPi0Comb7TeV = 
    (TGraphAsymmErrors*)fileppCombined.Get("graphInvCrossSectionPi0Comb7TeV");

  TGraphAsymmErrors* graphInvCrossSectionPi0Comb2760GeVSys = 
    (TGraphAsymmErrors*)fileppCombined.Get("graphInvCrossSectionPi0Comb2760GeVSysErr");
  TGraphAsymmErrors* graphInvCrossSectionPi0Comb7TeVSys = 
    (TGraphAsymmErrors*)fileppCombined.Get("graphInvCrossSectionPi0Comb7TeVSysErr");
                                            
  TGraphAsymmErrors* graphInvCrossSectionPi0Comb2760GeVSta = 
    (TGraphAsymmErrors*)fileppCombined.Get("graphInvCrossSectionPi0Comb2760GeVStatErr");
  TGraphAsymmErrors* graphInvCrossSectionPi0Comb7TeVSta = 
    (TGraphAsymmErrors*)fileppCombined.Get("graphInvCrossSectionPi0Comb7TeVStatErr");

  graphInvCrossSectionPi0Comb2760GeVSta->RemovePoint(18);
  graphInvCrossSectionPi0Comb2760GeVSys->RemovePoint(18);
  graphInvCrossSectionPi0Comb2760GeV->RemovePoint(18);

  TCanvas* canvasDummy2 = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
  DrawGammaCanvasSettings( canvasDummy2,  0.15, 0.01, 0.015, 0.08);
  canvasDummy2->SetLogy();
  canvasDummy2->SetLogx();
  TH2F * histo2DDummy2;
  histo2DDummy2 = new TH2F("histo2DDummy2","histo2DDummy2",1000,0.23,30.,1000,5e1,2e12);
  SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 0.8,1.55);
  histo2DDummy2->DrawCopy(); 
  
  DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb7TeVSta, markerStyleCommmonSpectrumPi07TeV,0.5*markerSizeCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kTRUE);
  DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb7TeVSys, markerStyleCommmonSpectrumPi07TeV,0.5*markerSizeCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes,kTRUE,colorCommonSpectrumPi07TeVBox );

  graphInvCrossSectionPi0Comb7TeV->SetLineWidth(widthCommonErrors);	
  //  graphInvCrossSectionPi0Comb7TeV->Draw("p,E2same");
  graphInvCrossSectionPi0Comb7TeVSta->SetMarkerColor(colorCommonSpectrumPi07TeV);
  graphInvCrossSectionPi0Comb7TeVSys->SetMarkerColor(colorCommonSpectrumPi07TeV);
  graphInvCrossSectionPi0Comb7TeVSys->Draw("p,E2same");
  graphInvCrossSectionPi0Comb7TeVSta->Draw("p,same");

  // DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb2760GeV, markerStyleCommmonSpectrumPi07TeV,markerSizeCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kTRUE);

  graphInvCrossSectionPi0Comb2760GeV->SetLineWidth(widthCommonErrors);	
  //  graphInvCrossSectionPi0Comb2760GeV->SetMarkerColor(kMagenta+1);
  //graphInvCrossSectionPi0Comb2760GeV->Draw("p,E2same");

  DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb2760GeVSta, markerStyleCommmonSpectrumPi07TeV,0.5*markerSizeCommonSpectrumPi07TeV, kMagenta+1, kMagenta+1, widthCommonSpectrumBoxes, kTRUE);
  DrawGammaSetMarkerTGraphAsym(graphInvCrossSectionPi0Comb2760GeVSys, markerStyleCommmonSpectrumPi07TeV,0.5*markerSizeCommonSpectrumPi07TeV, kMagenta+1, kMagenta+1, widthCommonSpectrumBoxes,kTRUE,kMagenta-8 );

  graphInvCrossSectionPi0Comb2760GeVSta->Draw("p,same");
  graphInvCrossSectionPi0Comb2760GeVSys->Draw("p,E2same");



  DrawGammaNLOTGraph(graphInvCrossSecCGC2760GeV , widthCommonFit, 2,kMagenta+2 );
  graphInvCrossSecCGC2760GeV->SetLineWidth(3);

  graphInvCrossSecCGC2760GeV->Draw("same,c");


  DrawGammaNLOTGraph(graphInvCrossSecCGC7000GeV , widthCommonFit, 2,colorCommonSpectrumPi07TeV);
  graphInvCrossSecCGC7000GeV->SetLineWidth(3);
  graphInvCrossSecCGC7000GeV->Draw("same,c");

  DrawGammaNLOTGraph(graphInvCrossSecCGC7000GeV_mvgamma , widthCommonFit, 3,colorCommonSpectrumPi07TeV);
  graphInvCrossSecCGC7000GeV_mvgamma->SetLineWidth(3);
  graphInvCrossSecCGC7000GeV_mvgamma->Draw("same,c");


  TLegend * leg= new TLegend(0.18,0.15,0.4,0.35);
  leg->SetTextSize(0.0325);
  leg->AddEntry(graphInvCrossSectionPi0Comb7TeV,"ALICE pp #sqrt{#it{s}} = 7 TeV","p");
  leg->AddEntry(graphInvCrossSectionPi0Comb2760GeV,"ALICE pp #sqrt{#it{s}} = 2.76 TeV","p");
  leg->AddEntry(graphInvCrossSecCGC7000GeV,"CGC MV^{e} #sqrt{#it{s}} = 7 TeV, Phys. Rev. D 88, 114020 (2013)","l");
  leg->AddEntry(graphInvCrossSecCGC7000GeV_mvgamma,"CGC MV^{#gamma} #sqrt{#it{s}} = 7 TeV","l");
  leg->AddEntry(graphInvCrossSecCGC2760GeV,"CGC MV^{e} #sqrt{#it{s}} = 2.76 TeV","l");

  leg->SetFillColor(0);
  leg->SetLineColor(0);

  leg->Draw();

  canvasDummy2->Print("CGCComparisonpp.pdf");
  canvasDummy2->Print("CGCComparisonpp.eps");
  
  TF1* fitInvCrossSectionPi0Comb7TeV = FitObject("l","fitInvCrossSectionPi0Comb7TeV","Pi0");
  fitInvCrossSectionPi0Comb7TeV->SetParameter(0,7.5e10);
  fitInvCrossSectionPi0Comb7TeV->SetParameter(1,7.5);
  fitInvCrossSectionPi0Comb7TeV->SetParameter(2,0.152);
  fitInvCrossSectionPi0Comb7TeV->SetRange(0.,24.);

  graphInvCrossSectionPi0Comb7TeV->Fit(fitInvCrossSectionPi0Comb7TeV,"QNRME+","",0.4,24.);
  //fitInvCrossSectionPi0Comb7TeV->Draw("same");

  Double_t paramPi0pp2760GeV[3]={1.1e11,5.8,0.13};

  TF1* fitInvCrossSectionPi0Comb2760GeV = FitObject("l","fit2760","Pi0",graphInvCrossSectionPi0Comb2760GeV,0.4,10.,paramPi0pp2760GeV,"QNRME+");
  

  //  fit2760->Draw("same");
  


  TGraph * graphRatioCombCGC2760GeV = (TGraph*)graphInvCrossSecCGC2760GeV->Clone("graphRatioCombCGC2760GeV");

  	
  graphRatioCombCGC2760GeV = CalculateGraphRatioToFit (graphRatioCombCGC2760GeV, 
							  fitInvCrossSectionPi0Comb2760GeV); 


  TGraph * graphRatioCombCGC7TeV = (TGraph*)graphInvCrossSecCGC7000GeV->Clone("graphRatioCombCGC7TeV");
	
  graphRatioCombCGC7TeV = CalculateGraphRatioToFit (graphRatioCombCGC7TeV, 
						       fitInvCrossSectionPi0Comb7TeV); 

  TGraph * graphRatioCombCGC7TeV_mvgamma = (TGraph*)graphInvCrossSecCGC7000GeV_mvgamma->Clone("graphRatioCombCGC7TeV_mvgamma");
	
  graphRatioCombCGC7TeV_mvgamma = CalculateGraphRatioToFit (graphRatioCombCGC7TeV_mvgamma, 
						       fitInvCrossSectionPi0Comb7TeV); 



  TGraphAsymmErrors * graphRatioCombCombFit7TeV = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb7TeV->Clone("graphRatioCombCombFit7TeV");
	
  graphRatioCombCombFit7TeV = CalculateGraphErrRatioToFit (graphInvCrossSectionPi0Comb7TeV, 
						       fitInvCrossSectionPi0Comb7TeV); 

  TGraphAsymmErrors * graphRatioCombCombFit7TeVSys = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb7TeVSys->Clone();
  graphRatioCombCombFit7TeVSys = CalculateGraphErrRatioToFit(graphRatioCombCombFit7TeVSys, fitInvCrossSectionPi0Comb7TeV); 


  TGraphAsymmErrors * graphRatioCombCombFit7TeVSta = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb7TeVSta->Clone();
  graphRatioCombCombFit7TeVSta = CalculateGraphErrRatioToFit(graphRatioCombCombFit7TeVSta, fitInvCrossSectionPi0Comb7TeV); 




  DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFit7TeVSta, markerStyleCommmonSpectrumPi07TeV,0.8*markerSizeCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kFALSE);
  graphRatioCombCombFit7TeVSta->SetLineWidth(widthCommonErrors);


  DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFit7TeVSys, markerStyleCommmonSpectrumPi07TeV,0.8*markerSizeCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kTRUE,colorCommonSpectrumPi07TeVBox);
  graphRatioCombCombFit7TeVSys->SetLineWidth(widthCommonErrors);








  TGraphAsymmErrors * graphRatioCombCombFit2760GeVSys = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb2760GeVSys->Clone();
  graphRatioCombCombFit2760GeVSys = CalculateGraphErrRatioToFit(graphRatioCombCombFit2760GeVSys, fitInvCrossSectionPi0Comb2760GeV); 


  TGraphAsymmErrors * graphRatioCombCombFit2760GeVSta = (TGraphAsymmErrors*)graphInvCrossSectionPi0Comb2760GeVSta->Clone();
  graphRatioCombCombFit2760GeVSta = CalculateGraphErrRatioToFit(graphRatioCombCombFit2760GeVSta, fitInvCrossSectionPi0Comb2760GeV); 




  DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFit2760GeVSta, markerStyleCommmonSpectrumPi07TeV,0.8*markerSizeCommonSpectrumPi07TeV, kMagenta+2, kMagenta+2, widthCommonSpectrumBoxes, kFALSE);
  graphRatioCombCombFit2760GeVSta->SetLineWidth(widthCommonErrors);


  DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFit2760GeVSys, markerStyleCommmonSpectrumPi07TeV,0.8*markerSizeCommonSpectrumPi07TeV, kMagenta+2, kMagenta+2, widthCommonSpectrumBoxes, kTRUE,kMagenta-8);
  graphRatioCombCombFit2760GeVSys->SetLineWidth(widthCommonErrors);









  //  graphRatioCombCombFit7TeV->SetMarkerStyle(20);
 
  
  //*************** Size factors ********************
  Double_t scaleMarkerNLO = 1.5;
  Double_t scaleLineWidthNLO = 1.;
  gStyle->SetOptTitle(kFALSE);
  TCanvas* canvasInvXSectionCCG = new TCanvas("canvasInvXSectionCGC","",200,10,1200,1200);  // gives the page size
    DrawGammaCanvasSettings( canvasInvXSectionCGC,  0.15, 0.01, 0.015, 0.08);
  //DrawGammaCanvasSettings( canvasInvXSectionCGC,  0.1, 0.02, 0.05, 0.05);

    //  TPad* padXSectionCGCPi07TeV = new TPad("padXSectionCGCPi07TeV", "", 0., 0.5, 1., 0.9,-1, -1, -2);
  TPad* padXSectionCGCPi07TeV = new TPad("padXSectionCGCPi07TeV", "", 0., 0.5, 1., 1.,-1, -1, -2);
  //  DrawGammaPadSettings( padXSectionCGCPi07TeV, 0.15, 0.01, 0.15, 0.);
  DrawGammaPadSettings( padXSectionCGCPi07TeV, 0.12, 0.02, 0.15, 0.00);
  padXSectionCGCPi07TeV->Draw();
	
  TPad* padXSectionCGCPi02760GeV = new TPad("padXSectionCGCPi02760GeV", "", 0., 0., 1., 0.5,-1, -1, -2);
  //  DrawGammaPadSettings( padXSectionCGCPi02760GeV, 0.15, 0.01, 0., 0.15);
  DrawGammaPadSettings( padXSectionCGCPi02760GeV, 0.12, 0.02, 0.0, 0.15);
  padXSectionCGCPi02760GeV->Draw();
	
	

  //  canvasInvXSectionCCG->SetLogx();
  padXSectionCGCPi07TeV->SetLogx();
  padXSectionCGCPi02760GeV->SetLogx();
  padXSectionCGCPi07TeV->cd();

  TH2F * ratio2DInvXSectionCGCPi0 = new TH2F("ratio2DInvXSectionCGCPi0","ratio2DInvXSectionCGCPi0",1000,0.23,30.,1000,0.3,3.4);
  SetStyleHistoTH2ForGraphs(ratio2DInvXSectionCGCPi0, "#it{p}_{T} (GeV/#it{c})","#frac{CGC}{data fit}", 0.07,0.07, 0.07 , 0.07, 0.9, 0.8, 512, 505);
  ratio2DInvXSectionCGCPi0->GetYaxis()->SetMoreLogLabels(kTRUE);
  ratio2DInvXSectionCGCPi0->GetYaxis()->SetNdivisions(515);
  ratio2DInvXSectionCGCPi0->GetYaxis()->SetNoExponent(kTRUE);
  ratio2DInvXSectionCGCPi0->GetXaxis()->SetTickLength(0.07);
  ratio2DInvXSectionCGCPi0->SetXTitle("#it{p}_{T}(GeV/#it{c})");
  ratio2DInvXSectionCGCPi0->SetYTitle("#frac{CGC}{data fit}");
  ratio2DInvXSectionCGCPi0->DrawCopy(); 
  graphRatioCombCombFit7TeVSta->SetMarkerSize(1.2);
  graphRatioCombCombFit7TeVSta->Draw("p,same");
  graphRatioCombCombFit7TeVSys->Draw("2,same");
  graphRatioCombCombFit7TeVSta->Draw("p,same");


 



  DrawGammaNLOTGraph(graphRatioCombCGC7TeV, widthCommonFit, 2,colorCommonSpectrumPi07TeV);
  graphRatioCombCGC7TeV->SetLineWidth(3);
  graphRatioCombCGC7TeV->Draw("same,c");
  

  DrawGammaNLOTGraph(graphRatioCombCGC7TeV_mvgamma, widthCommonFit, 3,colorCommonSpectrumPi07TeV);
  graphRatioCombCGC7TeV_mvgamma->SetLineWidth(3);
  graphRatioCombCGC7TeV_mvgamma->Draw("same,c");
 


  padXSectionCGCPi02760GeV->cd();

  ratio2DInvXSectionCGCPi0->DrawCopy(); 
  graphRatioCombCombFit2760GeVSta->SetMarkerSize(1.2);
  graphRatioCombCombFit2760GeVSta->Draw("p,same");
  graphRatioCombCombFit2760GeVSys->Draw("2,same");
  graphRatioCombCombFit2760GeVSta->Draw("p,same");

  DrawGammaNLOTGraph(graphRatioCombCGC2760GeV , widthCommonFit, 2, kMagenta+2 );
  graphRatioCombCGC2760GeV->SetLineWidth(3);
  graphRatioCombCGC2760GeV->Draw("same,c");

  TLegend * leg1= new TLegend(0.2,0.55,0.3,0.88);
  leg1->SetTextSize(0.06);
  leg1->AddEntry(graphRatioCombCombFit7TeVSta,"ALICE pp #sqrt{#it{s}} = 7 TeV","p");
  leg1->AddEntry(graphRatioCombCombFit2760GeVSta,"ALICE pp #sqrt{#it{s}} = 2.76 TeV","p");
  leg1->AddEntry(graphRatioCombCGC7TeV,"CGC MV^{e} #sqrt{#it{s}} = 7 TeV, Phys. Rev. D 88, 114020 (2013)","l");
  leg1->AddEntry(graphRatioCombCGC7TeV_mvgamma,"CGC MV^{#gamma} #sqrt{#it{s}} = 7 TeV ","l");
  leg1->AddEntry(graphRatioCombCGC2760GeV,"CGC MV^{e} #sqrt{#it{s}} = 2.76 TeV","l");

  leg1->SetFillColor(0);
  leg1->SetLineColor(0);

  leg1->Draw();

  canvasInvXSectionCCG->Print("ratioCGCpp7TeV2760GeV.eps");
  canvasInvXSectionCCG->Print("ratioCGCpp7TeV2760GeV.pdf");


}
