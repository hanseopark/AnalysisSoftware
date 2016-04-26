
/****************************************************************************************************************************
******          provided by Gamma Conversion Group, PWG-GA                                                                                                    *****
******          Ana Marin, marin@physi.uni-heidelberg.de                                                                                                        *****                                                                                                     *****
******          Friederike Bock, friederike.bock@cern.ch                                                                                                        *****
******          Annika Passfeld, annikapassfeld@uni-muenster,de
******          Pedro Gonzalez, pedro.gonzalez.zamora@cern.ch
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
#include "TArrow.h"
#include "TGraphAsymmErrors.h" 
#include "TGaxis.h"
#include "TMarker.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"
#include "TFitResultPtr.h"

/*#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctions.h"*/

extern TRandom*         gRandom;	
extern TBenchmark*      gBenchmark;
extern TSystem*         gSystem;
extern TMinuit*         gMinuit;
// TGraphErrors* RemovePointsFromGraph(TGraphErrors*, Int_t);
// TGraphAsymmErrors* GetChargeParticlesRpPb();
// TGraphErrors *GetInterpolSpectrum3D(TGraphErrors *g1, TGraphErrors *g2,TGraphErrors *g3, Double_t d1, Double_t d2, Double_t d3, Double_t dSqrts);
// TGraphErrors *GetInterpolSpectrum2D(TGraphErrors *g1, TGraphErrors *g2,Double_t d1, Double_t d2,Double_t dSqrts);


void InterpolationV3(TString fileNameNeutralPionPCMResultsPP="",TString fileNameNeutralPionPCMResultspPb="",TString fileNameNeutralPionDalitzResultspPb="", TString suffix="eps", TString outputDir="eps/OutputRpPbV3")
{



  //  gROOT->Reset(); 
        gROOT->SetStyle("Plain");

        TH1::AddDirectory(kFALSE);
	
	StyleSettings();


        TFile*  fileNeutralPionPCMResultsPP = new TFile( fileNameNeutralPionPCMResultsPP.Data() );
	
	
	Double_t	xSection2760GeVpp    =  62.8*1e-3;
        Double_t	xSection2760GeVErrpp = 	3.9;
	Double_t        xSectionpPb5023GeV   =  70*1e-3;
	Double_t	xSection7TeV         =	73.2*1e-3;
	Double_t        recalcBarn           =  1e12; //NLO in pbarn!!!!
	Double_t        factorToInel         =  1;
	
	

        TGraphAsymmErrors*  graphInvYieldPi0PCMStatPP7TeV      = (TGraphAsymmErrors*)fileNeutralPionPCMResultsPP->Get("graphInvCrossSectionPi0Comb7TeV");
        TGraphAsymmErrors*  graphInvYieldPi0PCMStatPP2760GeV   = (TGraphAsymmErrors*)fileNeutralPionPCMResultsPP->Get("graphInvCrossSectionPi0Comb2760GeV");	

		
    //     TGraphAsymmErrors*  graphInvYieldPi0PCMStatPP7TeV      = (TGraphAsymmErrors*)fileNeutralPionPCMResultsPP->Get("graphInvCrossSectionPi0PCMStat7TeV");
// 	TGraphAsymmErrors*  graphInvYieldPi0PCMStatPP2760GeV   = (TGraphAsymmErrors*)fileNeutralPionPCMResultsPP->Get("graphInvCrossSectionPi0PCM2760GeVStatErr");
	
        //TGraphAsymmErrors*  graphInvYieldPi0PCMStatPP2760GeV   = (TGraphAsymmErrors*)fileNeutralPionPCMResultsPP->Get("graphInvCrossSectionPi0PCMStat2760GeV_YShifted");


	
	
	graphInvYieldPi0PCMStatPP2760GeV   =  ScaleGraph(graphInvYieldPi0PCMStatPP2760GeV,1./(xSection2760GeVpp*recalcBarn)*factorToInel);
	graphInvYieldPi0PCMStatPP7TeV      =  ScaleGraph(graphInvYieldPi0PCMStatPP7TeV,1./(xSection7TeV*recalcBarn)*factorToInel);

	
	
	Int_t nPointsInter=10;        
        Int_t nOffset7000TeV=0;	
	Int_t nOffset2760TeV=0;
	
	
	cout<<" 7Tev" << endl;
	
	graphInvYieldPi0PCMStatPP7TeV->Print();
	
	nPointsInter = graphInvYieldPi0PCMStatPP7TeV->GetN();
	
	
		
	/////////////////Convert TGraphAsymmErrors to TGrahErrors////////////
	
	cout<<"*****************Creating TGraphErrors from TGraphAsymmErrors of 7 TeV**************************"<<endl;
	
	TGraphErrors* graphErrorsInvYieldPi0PCMStat7TeV = new TGraphErrors( nPointsInter );
	
	
	for(Int_t iPoint=0; iPoint < nPointsInter; iPoint++){
	  
	    Double_t x,y,errorX,errorY;
	    graphInvYieldPi0PCMStatPP7TeV->GetPoint(iPoint+nOffset7000TeV,x,y);
	    errorX = graphInvYieldPi0PCMStatPP7TeV->GetErrorX(iPoint+nOffset7000TeV);
	    errorY = graphInvYieldPi0PCMStatPP7TeV->GetErrorY(iPoint+nOffset7000TeV);
	    graphErrorsInvYieldPi0PCMStat7TeV->SetPoint(iPoint,x,y);
	    graphErrorsInvYieldPi0PCMStat7TeV->SetPointError(iPoint,errorX,errorY);
            cout<<iPoint<<" "<< x<< " " << y<<" " << errorX<< " "<< errorY<<endl;   
	    
	 }
	cout<< endl;
	
	
	cout<<"7 TeV TGraphErrors"<<endl;
	
	graphErrorsInvYieldPi0PCMStat7TeV->Print();
	
	cout<<"*************************************************************************************************"<<endl;
	
	cout<<"2.76 TeV Bin shifted"<<endl;    
	graphInvYieldPi0PCMStatPP2760GeV->Print();
	
	//////////////////////Fitting pp@2.76 InvYield spectra to extend the binning to pp@7TeV///////////////
	
	
	cout<<" 2.76TeV Fitted" << endl;
	cout<<"Fitting 2.76 TeV ..."<<endl;
	
	Float_t xMin = 0.30;
	Float_t xMax = graphInvYieldPi0PCMStatPP2760GeV->GetXaxis()->GetBinUpEdge(graphInvYieldPi0PCMStatPP2760GeV->GetXaxis()->GetNbins());
	
	cout<<"fitTsallisPi0PCM2760 xMin: "<<xMin<<"   xMax:  "<<xMax<<endl;
	
	Double_t Parameters[3]={1.0,7.,0.13}; //2.,5.,0.18
	TF1* fitTsallisPi0PCM2760 = FitObject("l","fitTsallisPi0PCM2760","Pi0");
	
	
	fitTsallisPi0PCM2760->SetRange(xMin,xMax);
	fitTsallisPi0PCM2760->SetParameters(Parameters[0],Parameters[1],Parameters[2]); // standard
	
	//SQNRME+ , SEMO
	
	TFitResultPtr resultfitTsallisPi0PCM2760 = graphInvYieldPi0PCMStatPP2760GeV->Fit(fitTsallisPi0PCM2760,"SQNRME+","",xMin,xMax);
	
	
	
	
	
	
	nOffset7000TeV =  0;
	nOffset2760TeV = -2;
	
	//nPointsInter = graphErrorsInvYieldPi0PCMStat7TeV->GetN();  
	
	nPointsInter =  graphInvYieldPi0PCMStatPP7TeV->GetN();
	
	cout<<"*****************Creating TGraphError from TGraphAsymmErrors and Fit of 2.760 TeV***************"<<endl;
	
        TGraphErrors*  graphErrorsInvYieldPi0PCMStat2760GeVFitted = new TGraphErrors( nPointsInter );
	
	for(Int_t iPoint = 0; iPoint < nPointsInter; iPoint++){
	  
	 
	  Double_t x7TeV,y7TeV,errorXlow7TeV,errorXhigh7TeV;
	  Double_t x2760GeV, y2760GeV;
	  Double_t xE2760GeV,yE2760GeV;
	  
	  graphInvYieldPi0PCMStatPP7TeV->GetPoint(iPoint,x7TeV,y7TeV);
	  
	  errorXlow7TeV  = graphInvYieldPi0PCMStatPP7TeV->GetErrorXlow(iPoint+nOffset7000TeV);
	  errorXhigh7TeV = graphInvYieldPi0PCMStatPP7TeV->GetErrorXhigh(iPoint+nOffset7000TeV);
	 
	  
	  
	
	 //cout<<"iPoint: "<<iPoint<<" 7TeV Bin:  xVal: "<<x7TeV <<" yVal: "<<y7TeV<<endl;
	  
	  if( iPoint < 3 || iPoint > 12 ) {
	    
	    //Double_t ptMin = x7TeV - errorXlow7TeV;
	    //Double_t ptMax = x7TeV + errorXhigh7TeV;
	    //Double_t binWidth = ptMax - ptMin;

	      x2760GeV = x7TeV;
	      Double_t err[1];
	      Double_t val[1];
	      val[0] = x2760GeV;
	    
	    //Double_t yInt =      fitTsallisPi0PCM2760->Integral(ptMin,ptMax,resultfitTsallisPi0PCM2760->GetParams()) / binWidth;
	    //Double_t yIntError = fitTsallisPi0PCM2760->IntegralError(ptMin,ptMax, resultfitTsallisPi0PCM2760->GetParams(), resultfitTsallisPi0PCM2760->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
	      Double_t yInt = fitTsallisPi0PCM2760->Eval(val[0]);
	      resultfitTsallisPi0PCM2760->GetConfidenceIntervals(1, 1, 1, val, err, 0.95, false);  
	      Double_t yIntError = err[0];
	    
	    
	
	     cout<<"From Fit "<<"iPoint: "<<iPoint<<"2760 GeV Bin:  xVal: "<<x2760GeV<< " yVal: "<<yInt<<" errorY: "<<yIntError<<endl;
	     
	     graphErrorsInvYieldPi0PCMStat2760GeVFitted->SetPoint(iPoint,x2760GeV,yInt);
	     graphErrorsInvYieldPi0PCMStat2760GeVFitted->SetPointError(iPoint,errorXhigh7TeV,yIntError);
	     
		    
	  }
	  else {
	    
	    
		graphInvYieldPi0PCMStatPP2760GeV->GetPoint(iPoint+nOffset2760TeV,x2760GeV,y2760GeV);
		xE2760GeV = graphInvYieldPi0PCMStatPP2760GeV->GetErrorX(iPoint+nOffset2760TeV);
	        yE2760GeV = graphInvYieldPi0PCMStatPP2760GeV->GetErrorY(iPoint+nOffset2760TeV);
		cout<<"From TGraphAsymmErrors "<<"iPoint: "<<iPoint<<" 2760 GeV Bin:  xVal: "<<x2760GeV<< " yVal: "<<y2760GeV<<" errorY: "<<yE2760GeV<<endl;
		graphErrorsInvYieldPi0PCMStat2760GeVFitted->SetPoint(iPoint,x2760GeV,y2760GeV);
		graphErrorsInvYieldPi0PCMStat2760GeVFitted->SetPointError(iPoint,xE2760GeV,yE2760GeV);
			    
	  }
	  
	}
	cout<<endl;
	 
	
	
	
	////////////////////Interpolate pp@5.023 TeV//////////////////// 
	 
        TGraphErrors* graphErrosInterPolation5023GeV2TPoints = GetInterpolSpectrum2D(graphErrorsInvYieldPi0PCMStat2760GeVFitted,graphErrorsInvYieldPi0PCMStat7TeV,2760,7000,5023);
	       
	
	
	/////////////Drawing InvYield pp@2.76 GeV//////////////
	
	
	TCanvas* cInvYield2760GeV = new TCanvas("cInvYield2760GeV","Invariant Yield pp@2.76 TeV",200,10,700,550);
    	DrawGammaCanvasSettings( cInvYield2760GeV,  0.15, 0.02, 0.03, 0.1);
	cInvYield2760GeV->SetLogx();
	cInvYield2760GeV->SetLogy();
	
	TPad* padComparisonInvYield2760GeV = new TPad("padComparisonInvYield2760GeV", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padComparisonInvYield2760GeV, 0.15, 0.02, 0.03, 0.1);
	padComparisonInvYield2760GeV->Draw();
	
	
	TH2F * histoInvYield2760GeV;
	histoInvYield2760GeV = new TH2F("histoInvYield2760GeV","histoInvYield2760GeV",1000,0.23,30.,1000,2e-9,10e0);
	SetStyleHistoTH2ForGraphs(histoInvYield2760GeV, "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoInvYield2760GeV->DrawCopy(); 
	

	
	DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0PCMStatPP2760GeV, 20, 1, kRed, kRed);
	DrawGammaSetMarkerTGraphErr(graphErrorsInvYieldPi0PCMStat2760GeVFitted, 20,1,kBlue,kBlue);
	
	
	
	
	
	graphInvYieldPi0PCMStatPP2760GeV->Draw("p,same,e1");
	graphErrorsInvYieldPi0PCMStat2760GeVFitted->Draw("p,same,e1");
	fitTsallisPi0PCM2760->Draw("same");
	
	Double_t columnsLegend[4] 	= {0.,0.18,0.47,0.75};
	Double_t rowsLegend[6] 		= {0.88,0.75,0.57,0.4,0.22,0.05}; //with EMCAL {0.88,0.75,0.57,0.4,0.22,0.05};
		//******************* Text sizes *******************
	  Size_t textSizeLeftColumn	= 0.13;
	  Size_t textSizeTopRow		= 0.13; 
	  Size_t textSizeSecondRow 	= 0.11;
		//******************* Offsets ***********************
	  Double_t offsetSystColumn 	= 0.15;
	  Double_t offsetMarkerX		= 0.41;
	  Double_t offsetMarkerY		= 0.05;
	  Double_t offsetBoxSizeY		= 0.05;
	  Double_t offsetFit		= 0.04;
		//****************** Scale factors ******************
	  Double_t scaleWidthLine 		= 0.8;
		
		TPad* padInvYieldPP2760GeVLegend = new TPad("padInvYieldPP2760GeVLegend", "", 0.57, 0.7, 0.85, 0.95,-1, -1, -2);  
		DrawGammaPadSettings( padInvYieldPP2760GeVLegend, 0., 0., 0., 0.);
		padInvYieldPP2760GeVLegend->Draw();
		padInvYieldPP2760GeVLegend->cd();
	
	
	Double_t	 columnsLegend2[4] 	= {0.,0.18,0.47,0.75};
	Double_t	      rowsLegend2[6] 		= {0.88,0.75,0.57,0.4,0.22,0.05}; //with EMCAL {0.88,0.75,0.57,0.4,0.22,0.05};
	
		
	
	        TLatex *textpp2760GeV = new TLatex(columnsLegend2[0],rowsLegend2[2],"pp@2.76 TeV");
			
		
		SetStyleTLatex( textpp2760GeV, textSizeLeftColumn,4);
		
		textpp2760GeV->Draw();
		TLatex *textpp2760GeVFitted = new TLatex(columnsLegend2[0],rowsLegend2[3],"pp@2.76 TeV Fitted");
		SetStyleTLatex( textpp2760GeVFitted, textSizeLeftColumn,4);
		textpp2760GeVFitted->Draw();
		TLatex *textTsallisFit = new TLatex(columnsLegend2[0],rowsLegend2[4],"Tsallis fit");
		SetStyleTLatex( textTsallisFit, textSizeLeftColumn,4);
		textTsallisFit->Draw();
		
		
		TMarker* markerPi0PCMStatPP2760GeV= CreateMarkerFromGraph(graphInvYieldPi0PCMStatPP2760GeV,columnsLegend2[2]+ offsetMarkerX ,rowsLegend2[2]+ offsetMarkerY ,scaleWidthLine);
		markerPi0PCMStatPP2760GeV->DrawMarker(columnsLegend2[2]+ offsetMarkerX ,rowsLegend2[2]+ offsetMarkerY);
		TMarker* markerPi0PCMStatPP2760GeVFitted = CreateMarkerFromGraph(graphErrorsInvYieldPi0PCMStat2760GeVFitted, columnsLegend2[2]+ offsetMarkerX ,rowsLegend2[3]+ offsetMarkerY ,scaleWidthLine);
		markerPi0PCMStatPP2760GeVFitted->DrawMarker(columnsLegend2[2]+ offsetMarkerX ,rowsLegend2[3]+ offsetMarkerY);
		TLine * lineFit27600GeV = CreateLineFromFit(fitTsallisPi0PCM2760, columnsLegend2[2]+ offsetMarkerX- offsetFit, rowsLegend2[4]+offsetMarkerY, columnsLegend2[2]+ offsetMarkerX+ offsetFit, rowsLegend2[4]+offsetMarkerY, scaleWidthLine);
		lineFit27600GeV->Draw("same");
		
		TBox* boxPi0PCMStatPP2760GeV = CreateBoxFromGraph(graphInvYieldPi0PCMStatPP2760GeV,columnsLegend2[2]+offsetSystColumn+offsetMarkerX-offsetFit, rowsLegend2[2]+ offsetMarkerY- offsetBoxSizeY, columnsLegend2[1]+ offsetSystColumn + offsetMarkerX + offsetFit, rowsLegend2[2]+ offsetMarkerY+ offsetBoxSizeY);
		boxPi0PCMStatPP2760GeV->Draw("l");
		
		TBox* boxPi0PCMStatPP2760GeVFitted = CreateBoxFromGraph(graphErrorsInvYieldPi0PCMStat2760GeVFitted,columnsLegend2[2]+offsetSystColumn+offsetMarkerX-offsetFit, rowsLegend2[3]+ offsetMarkerY- offsetBoxSizeY, columnsLegend2[1]+ offsetSystColumn+ offsetMarkerX+ offsetFit, rowsLegend2[3]+ offsetMarkerY+ offsetBoxSizeY);
		boxPi0PCMStatPP2760GeVFitted->Draw("l");
		
		//TBox* boxPi0PCMStatPP2760GeVFitted = CreateBoxFromGraph(graphErrorsInvYieldPi0PCMStat2760GeVFitted,columnsLegend2[2]+offsetSystColumn+offsetMarkerX-offsetFit, rowsLegend2[3]+ offsetMarkerY- offsetBoxSizeY, columnsLegend2[1]+ offsetSystColumn+ offsetMarkerX+ offsetFit, rowsLegend2[3]+ offsetMarkerY+ offsetBoxSizeY);
		//boxPi0PCMStatPP2760GeVFitted->Draw("l");
		
		
		
	
	
		//cout<<endl;
		//graphInvYieldPi0PCMStatPP2760GeV->Print();
		//cout<<endl;
		//graphErrorsInvYieldPi0PCMStat2760GeVFitted->Print();
		
		cInvYield2760GeV->Update();
		cInvYield2760GeV->SaveAs(Form("%s/InvYieldPP2760GeVCombined.%s",outputDir.Data(),suffix.Data()));
	
	
	////////////////////////////Scaling graphsErrors to be painted////////////////////
	
	
	TGraphErrors* graphErrorsInvYieldPi0PCMStat7TeVScale           = new TGraphErrors( *graphErrorsInvYieldPi0PCMStat7TeV );
	TGraphErrors* graphErrorsInvYieldPi0PCMStat2760GeVFittedScale  = new TGraphErrors( *graphErrorsInvYieldPi0PCMStat2760GeVFitted);
	TGraphErrors* graphErrosInterPolation5023GeV2TPointsScale      = new TGraphErrors( *graphErrosInterPolation5023GeV2TPoints );
	
	
	
	
	graphErrorsInvYieldPi0PCMStat7TeVScale          = ScaleGraph(graphErrorsInvYieldPi0PCMStat7TeVScale,4);
	graphErrosInterPolation5023GeV2TPointsScale     = ScaleGraph(graphErrosInterPolation5023GeV2TPointsScale,2);
	graphErrorsInvYieldPi0PCMStat2760GeVFittedScale = ScaleGraph(graphErrorsInvYieldPi0PCMStat2760GeVFittedScale,1);
	
        		
													 
	TCanvas* cInvYieldComparison = new TCanvas("cInvYieldComparison","Comparison Invariant Yields",200,10,700,500);
		
	DrawGammaCanvasSettings( cInvYieldComparison,  0.15, 0.02, 0.03, 0.1);
	cInvYieldComparison->SetLogx();
	cInvYieldComparison->SetLogy();
	
	TPad* padInvYieldComparison = new TPad("padInvYieldComparison", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padInvYieldComparison, 0.15, 0.02, 0.03, 0.1);
	padInvYieldComparison->Draw();
	
	
	TH2F * histoInvYieldComparison;
	histoInvYieldComparison = new TH2F("histoInvYieldComparison","histoInvYieldComparison",1000,0.23,30.,1000,2e-10,5e1);
	SetStyleHistoTH2ForGraphs(histoInvYieldComparison, "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoInvYieldComparison->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphErr(graphErrorsInvYieldPi0PCMStat7TeVScale,     	20, 1, kRed+1,   kRed+1);
	DrawGammaSetMarkerTGraphErr(graphErrosInterPolation5023GeV2TPointsScale, 	20, 1, kGreen+2, kGreen+2);
	DrawGammaSetMarkerTGraphErr(graphErrorsInvYieldPi0PCMStat2760GeVFittedScale,  	20, 1, kBlue+1,  kBlue+1);
	
	
	graphErrorsInvYieldPi0PCMStat7TeVScale->Draw("p,same,e1");
	graphErrosInterPolation5023GeV2TPointsScale->Draw("p,same,e1");
	graphErrorsInvYieldPi0PCMStat2760GeVFittedScale->Draw("p,same,e1");
	
	TLatex *LabelpPb = new TLatex(0.65,0.9,"ALICE work in progress");
	SetStyleTLatex( LabelpPb, 0.03,4);
	LabelpPb->Draw();
	TLatex *labelSpectraPi0LabelPCMpPb = new TLatex(0.65,0.85,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-} (PCM)");
	SetStyleTLatex( labelSpectraPi0LabelPCMpPb, 0.03,4);
	labelSpectraPi0LabelPCMpPb->Draw();

	TLatex *labelSpectraPi0LabelpPb = new TLatex(0.65,0.8,"|y_{#pi^{0},lab}| < 0.8");
	SetStyleTLatex( labelSpectraPi0LabelpPb, 0.03,4);
	labelSpectraPi0LabelpPb->Draw();  

	graphErrorsInvYieldPi0PCMStat7TeVScale->SetFillColor(0);
	graphErrosInterPolation5023GeV2TPointsScale->SetFillColor(0);
	graphErrorsInvYieldPi0PCMStat2760GeVFittedScale->SetFillColor(0);
	TLegend* legendInvYields = new TLegend(0.2,0.15,0.6,0.35);
	legendInvYields->SetFillColor(0);
	legendInvYields->SetLineColor(0);
	legendInvYields->SetNColumns(1);
	legendInvYields->SetTextSize(0.03);
	legendInvYields->AddEntry(graphErrorsInvYieldPi0PCMStat7TeVScale,"pp@7 TeV (x4)","pef");
	legendInvYields->AddEntry(graphErrosInterPolation5023GeV2TPointsScale,"pp@5.023 TeV - Interpolation (x2)","pef");
	legendInvYields->AddEntry(graphErrorsInvYieldPi0PCMStat2760GeVFittedScale,"pp@2.76 TeV (x1)","pef");
	
	legendInvYields->Draw();

		
	
	cout<<"TGraphErrors pp@5.023 TeV"<<endl;
	graphErrosInterPolation5023GeV2TPoints->Print();
	
	cInvYieldComparison->Update();
	cInvYieldComparison->SaveAs(Form("%s/InvYield_Comparison_PP_2760GeV_5023GeV_7TeV.%s",outputDir.Data(),suffix.Data()));
	
	
	 ///////////////////////////////////////////////////////////////////////////////
	 
	 
	///////////Getting the InvYield spectra p-Pb @ 5.023 TeV from PCM results
	
        TFile*  fileNeutralPionPCMResultspPb = new TFile(fileNameNeutralPionPCMResultspPb.Data());
	
	
	TDirectory* fNeutralPionPCMResultspPbContainer = (TDirectory*) fileNeutralPionPCMResultspPb->GetDirectory("Pi0_pPb_5.023TeV_0-100%");
	
	
	if( ! fNeutralPionPCMResultspPbContainer ) {cout<<"TList fNeutralPionPCMResultspPbContainer does not exist: "<<endl; return;}
	 
	 
	 	 
	TH1F* histoInvYieldPi0PCMpPb5023GeV  = (TH1F*)fNeutralPionPCMResultspPbContainer->Get("CorrectedYieldPi0");
   TGraphAsymmErrors* graphPCMYieldPi0SysErrpPb=    (TGraphAsymmErrors*)fNeutralPionPCMResultspPbContainer->Get("Pi0SystError"); 
 
   //   TGraphAsymmErrors* graphPCMYieldPi0SysErrpPbCopy = (TGraphAsymmErrors*) graphPCMYieldPi0SysErrpPb->Clone("Pi0SystErrorClone");	


	TGraphAsymmErrors* graphInvYieldPi0PCMpPb5023GeV = new TGraphAsymmErrors(histoInvYieldPi0PCMpPb5023GeV);


    TGraphAsymmErrors* graphStatCopy = (TGraphAsymmErrors*)graphInvYieldPi0PCMpPb5023GeV->Clone("graphStatCopy");
    TGraphAsymmErrors* graphSysCopy = (TGraphAsymmErrors*)graphPCMYieldPi0SysErrpPb->Clone("graphSysCopy");
    Double_t* xValue = graphStatCopy->GetX();
    Double_t* xErrorLow = graphStatCopy->GetEXlow();
    Double_t* xErrorHigh = graphStatCopy->GetEXhigh();
    Double_t* yValueStat = graphStatCopy->GetY();
    Double_t* yErrorLowStat = graphStatCopy->GetEYlow();
    Double_t* yErrorHighStat = graphStatCopy->GetEYhigh();
    Double_t* yErrorLowSys = graphSysCopy->GetEYlow();
    Double_t* yErrorHighSys = graphSysCopy->GetEYhigh();
    Int_t nPoints = graphStatCopy->GetN();
    cout << "nPoints: " << nPoints<< endl; 
    Double_t yErrorLowComb[31];
    Double_t yErrorHighComb[31];
     for (Int_t i = 0; i < nPoints; i++){
       yErrorLowComb[i] = TMath::Sqrt((TMath::Power(yErrorLowStat[i],2) + TMath::Power(yErrorLowSys[i],2)));
       yErrorHighComb[i] = TMath::Sqrt((TMath::Power(yErrorHighStat[i],2) + TMath::Power(yErrorHighSys[i],2)));
    
     }
     TGraphAsymmErrors* graphSysAndComb = new TGraphAsymmErrors(nPoints,xValue,yValueStat,xErrorLow,xErrorHigh,yErrorLowComb,yErrorHighComb);
 	//------------------------------------


 	graphSysAndComb->RemovePoint(0);
graphInvYieldPi0PCMpPb5023GeV->RemovePoint(0);
//		TGraphAsymmErrors* graphInvYieldPi0PCMpPb5023GeVXShifted = graphInvYieldPi0PCMpPb5023GeV->Clone();
 TGraphAsymmErrors* graphInvYieldPi0PCMpPb5023GeVXShifted =graphSysAndComb->Clone();
	
	//Double_t paramPCMpPb5023GeV[3]={1.0, 7., 0.13}; //1.0e10,7.,0.13
	
	Double_t paramPCMpPb5023GeV[10];
	
	ReturnParameterSetFittingPbPb("8000",paramPCMpPb5023GeV);
	
	
	TF1* fitTsallisPi0PCMpPb5023GeVPtMult = FitObject("rad","fitTsallisPi0PCMpPb5023GeVPtMult","Pi0");
	
	SetParametersLimitsForFit (fitTsallisPi0PCMpPb5023GeVPtMult, 5, paramPCMpPb5023GeV);
	
	 //Double_t Param[5] = { 1.27524, 6.64925, 0.155715, 3.02966, 5.9877};
	 
	  Double_t Param[5] = { 0.97524, 2.64925, 0.155715, 3.02966, 5.9877};
   
          fitTsallisPi0PCMpPb5023GeVPtMult->SetParameters(Param);
	
	// fitTsallisPi0PCMpPb5023GeVPtMult->SetRange(0.4,20);
	 /*for(Int_t i=0; i<10; i++){
	   fitTsallisPi0PCMpPb5023GeVPtMult->SetParameter(i,paramPCMpPb5023GeV[i]) ; // standard parameter optimize if necessary
	   cout<<"Param "<<i<<" "<<paramPCMpPb5023GeV[i]<<endl;
	 }*/
	     
	
        ////Applying BinShifX
	 
	 //graphInvYieldPi0PCMpPb5023GeVXShifted->Fit(fitTsallisPi0PCMpPb5023GeVPtMult,"QNRME+","",0.4,20);
	graphInvYieldPi0PCMpPb5023GeVXShifted  = ApplyXshift(graphInvYieldPi0PCMpPb5023GeVXShifted ,fitTsallisPi0PCMpPb5023GeVPtMult);
	
	
	cout<<"Parameters"<<endl;
	for(Int_t i=0; i<5; i++){
	cout<<"Param "<<i<<" "<<fitTsallisPi0PCMpPb5023GeVPtMult->GetParameter(i)<<endl;
	}
	
	     
	     
	cout<<"PCM pPb X shifted"<<endl;
	graphInvYieldPi0PCMpPb5023GeVXShifted->Print();
	
	
	
	
	TCanvas* cInvYieldPi0PCMpPb5023GeVXShifted = new TCanvas("cInvYieldPi0PCMpPb5023GeVXShifted","Pi0 PCM pPb at 5.23 TeV X shifted",200,10,700,500);
		
	DrawGammaCanvasSettings( cInvYieldPi0PCMpPb5023GeVXShifted,  0.15, 0.02, 0.03, 0.1);
	cInvYieldPi0PCMpPb5023GeVXShifted->SetLogx();
	cInvYieldPi0PCMpPb5023GeVXShifted->SetLogy();
	
	TPad* padInvYieldPi0PCMpPb5023GeVXShifted = new TPad("padInvYieldPi0PCMpPb5023GeVXShifted", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padInvYieldPi0PCMpPb5023GeVXShifted, 0.15, 0.02, 0.03, 0.1);
	padInvYieldPi0PCMpPb5023GeVXShifted->Draw();
	
	
	TH2F * histoInvYieldPi0PCMpPb5023GeVXShifted;
	histoInvYieldPi0PCMpPb5023GeVXShifted = new TH2F("histoInvYieldPi0PCMpPb5023GeVXShifted","histoInvYieldPi0PCMpPb5023GeVXShifted",1000,0.23,30.,1000,2e-9,5e1);
	SetStyleHistoTH2ForGraphs(histoInvYieldPi0PCMpPb5023GeVXShifted, "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoInvYieldPi0PCMpPb5023GeVXShifted->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0PCMpPb5023GeVXShifted,     	20, 1, kRed,   kRed);
	DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0PCMpPb5023GeV,  	        20, 1, kBlue,  kBlue);
	
	
	
	graphInvYieldPi0PCMpPb5023GeV->Draw("p,e1");
	graphInvYieldPi0PCMpPb5023GeVXShifted->Draw("p,same,e1");
	//fitTsallisPi0PCMpPb5023GeVPtMult->Draw("same");
	//fitTsallisPi0PCMpPb5023GeVPtMult->Draw("same");
	
	TLegend* legendInvYieldPi0PCMpPb5023GeVXShifted = new TLegend(0.20,0.25,0.70,0.40);
	legendInvYieldPi0PCMpPb5023GeVXShifted->SetFillColor(0);
	legendInvYieldPi0PCMpPb5023GeVXShifted->SetLineColor(0);
	legendInvYieldPi0PCMpPb5023GeVXShifted->SetTextSize(0.03);
	legendInvYieldPi0PCMpPb5023GeVXShifted->AddEntry(graphInvYieldPi0PCMpPb5023GeV,"#pi^{0} pPb@5.02 TeV (PCM)");
	legendInvYieldPi0PCMpPb5023GeVXShifted->AddEntry(graphInvYieldPi0PCMpPb5023GeVXShifted,"#pi^{0} pPb@5.02 TeV XBin shifted (PCM)");
	//legendInvYieldPi0PCMpPb5023GeVXShifted->AddEntry(fitTsallisPi0PCMpPb5023GeVPtMult,"Tsallis Fit");
	legendInvYieldPi0PCMpPb5023GeVXShifted->Draw("same");
	
		
	
	cInvYieldPi0PCMpPb5023GeVXShifted->Update();
	cInvYieldPi0PCMpPb5023GeVXShifted->SaveAs(Form("%s/InvYield_XShifted_PCM_pPb_5023GeV.%s",outputDir.Data(),suffix.Data()));
	
	
	
	
	
	
	
	
	//Fitting interpolated pp@5.023 TeV spectra
	
	cout<<"Fitting pp@5.023 TeV ..."<<endl;
	
	xMin = 0.2;
	xMax = graphErrosInterPolation5023GeV2TPoints->GetXaxis()->GetBinUpEdge(graphErrosInterPolation5023GeV2TPoints->GetXaxis()->GetNbins());
	
	
	
	cout<<"fitTsallisPi0PCMPP5023GeV xMin: "<<xMin<<"   xMax:  "<<xMax<<endl;
	
	Double_t ParametersPP5023GeV[3]={1.0,7.,0.13}; //1.0e10,7.,0.13
	
	
	TF1* fitTsallisPi0PCMPP5023GeV = FitObject("tmpt","fitTsallisPi0PCMPP5023GeV","Pi0");
	
	fitTsallisPi0PCMPP5023GeV->SetRange(xMin,xMax);
	fitTsallisPi0PCMPP5023GeV->SetParameters(ParametersPP5023GeV[0],ParametersPP5023GeV[1],ParametersPP5023GeV[2]); // standard
	
	TFitResultPtr resultfitTsallisPi0PCMPP5023GeV = graphErrosInterPolation5023GeV2TPoints->Fit(fitTsallisPi0PCMPP5023GeV,"SQNRME+","",xMin,xMax);
	
	
	
	
	nPointsInter = graphInvYieldPi0PCMpPb5023GeV->GetN();
			
			
	
			
	
	
	//Extending interpolated pp@5.023 TeV to the 7 TeV binning
	
	cout<<"*******************Extending interpolated pp@5.023 TeV to the 7 TeV binning ********************"<<endl;
	
	
	TGraphAsymmErrors* graphErrorsPi0PCMPP5023GeVFitted = new TGraphAsymmErrors( nPointsInter);
	
	
	  Double_t* xPoint     = graphInvYieldPi0PCMpPb5023GeVXShifted->GetX();
	  Double_t* errorXlow  = graphInvYieldPi0PCMpPb5023GeVXShifted->GetEXlow();
	  Double_t* errorXhigh = graphInvYieldPi0PCMpPb5023GeVXShifted->GetEXhigh();
	  
	
	cout<<"*******************llego aqui"<<endl;
	
	
	for(Int_t iPoint = 0; iPoint < nPointsInter; iPoint++){
	  
	 // cout<<iPoint<<" "<<xPoint[iPoint] <<" "<<errorXlow[iPoint]<<" "<<errorXhigh[iPoint]<<endl;
	  
	  //Double_t ptMin = xPoint[iPoint] - errorXlow[iPoint];
	  //Double_t ptMax = xPoint[iPoint] + errorXhigh[iPoint];
	  //Double_t binWidth = ptMax - ptMin;
														   
	 // Double_t yInt      = fitTsallisPi0PCMPP5023GeV->Integral(ptMin,ptMax,resultfitTsallisPi0PCMPP5023GeV->GetParams()) / binWidth;
	 // Double_t yIntError = fitTsallisPi0PCMPP5023GeV->IntegralError(ptMin,ptMax, resultfitTsallisPi0PCMPP5023GeV->GetParams(), resultfitTsallisPi0PCMPP5023GeV->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
	  
	  Double_t point[1] = {xPoint[iPoint]};
	  Double_t err[1];
	  Double_t yInt = fitTsallisPi0PCMPP5023GeV->Eval(point[0]);
	  resultfitTsallisPi0PCMPP5023GeV->GetConfidenceIntervals(1, 1, 1, point, err, 0.95, false);  
	  Double_t yIntError = err[0];
	  
	  	  
	  cout<<"Point "<<iPoint<<" xVal "<<xPoint[iPoint] <<" xErrlow "<<errorXlow[iPoint]<< " xErrhigh "<< errorXhigh[iPoint] <<" yVal "<<yInt<<" yErrlow "<<yIntError<<" yErrhigh "<<yIntError<<endl;
	  
	  graphErrorsPi0PCMPP5023GeVFitted->SetPoint(iPoint,xPoint[iPoint],yInt);
	  graphErrorsPi0PCMPP5023GeVFitted->SetPointError(iPoint,errorXlow[iPoint],errorXhigh[iPoint],yIntError,yIntError);
	  
	  
	  
	  
	}
	
	
	//Computing RpPb using the interpolated pp@5.023 TeV and PCM pp@7 TeV spectra
	
	TGraphAsymmErrors* graphRpPbPCM = new TGraphAsymmErrors( nPointsInter);
	
	Double_t *xBinspPb    = graphInvYieldPi0PCMpPb5023GeVXShifted->GetX();
	//Double_t *xBinsPP     = graphErrorsPi0PCMPP5023GeVFitted->GetX();
	Double_t *yBinspPb    = graphInvYieldPi0PCMpPb5023GeVXShifted->GetY();
	Double_t *yBinsPP     = graphErrorsPi0PCMPP5023GeVFitted->GetY();
	
	
	Double_t *xBinsErrpPblow  = graphInvYieldPi0PCMpPb5023GeVXShifted->GetEXlow();
	//Double_t *xBinsErrpPbhigh = graphInvYieldPi0PCMpPb5023GeVXShifted->GetEXhigh();
	Double_t *yBinsErrpPblow  = graphInvYieldPi0PCMpPb5023GeVXShifted->GetEYlow();
	//Double_t *yBinsErrpPbhigh = graphInvYieldPi0PCMpPb5023GeVXShifted->GetEYhigh();
	
	
	
	//Double_t *xBinsErrPPlow  = graphErrorsPi0PCMPP5023GeVFitted->GetEXlow();
	//Double_t *xBinsErrPPhigh = graphErrorsPi0PCMPP5023GeVFitted->GetEXhigh();
	Double_t *yBinsErrPPlow  = graphErrorsPi0PCMPP5023GeVFitted->GetEYlow();
	//Double_t *yBinsErrPPhigh = graphErrorsPi0PCMPP5023GeVFitted->GetEYhigh();
	
	
	
	
	
	Double_t fNcoll = 6.9;
	Double_t fTpPb = fNcoll/xSectionpPb5023GeV;// 0.0983/1e-3;
	
	
	for(Int_t iPoint = 0; iPoint < nPointsInter; iPoint++){
	  
	 graphRpPbPCM->SetPoint( iPoint, xBinspPb[iPoint], yBinspPb[iPoint] /(fNcoll* yBinsPP[iPoint]) );
	 
	 Double_t errYStat = pow( pow( yBinsErrpPblow[iPoint]/yBinspPb[iPoint],2. ) + pow( yBinsErrPPlow[iPoint]/yBinsPP[iPoint],2.), 0.5)*yBinspPb[iPoint] /(fNcoll* yBinsPP[iPoint]);
	 
	 graphRpPbPCM->SetPointError(iPoint, xBinsErrpPblow[iPoint],xBinsErrpPblow[iPoint],errYStat,errYStat);
		 
	 
	}
	  
	
	
	
	
	  //////////////////////Dalitz////////////////////
	  
	  
	//Getting Inv Yield spectro pPb@5.023 from Dalitz results
	
        TFile*  fileNeutralPionDalitzResultspPb = new TFile(fileNameNeutralPionDalitzResultspPb.Data());
	
	TDirectory* fNeutralPionDalitzResultspPbContainer = (TDirectory*) fileNeutralPionDalitzResultspPb->GetDirectory("Pi0_pPb_5.023TeV_0-100%");
	
	 if( ! fNeutralPionDalitzResultspPbContainer ) {cout<<"TList fNeutralPionDalitzResultspPbContainer does not exist: "<<endl; return;}
	 
	 	 
	TH1F* histoInvYieldPi0DalitzpPb5023GeV  = (TH1F*)fNeutralPionDalitzResultspPbContainer->Get("CorrectedYieldPi0");
	
	
	//cout<<"Histo p-Pb@5.023 TeV NBins "<<histoInvYieldPi0DalitzpPb5023GeV->GetNbinsX()<<endl;
	
	
        
	TGraphAsymmErrors* graphInvYieldPi0DalitzpPb5023GeV = new TGraphAsymmErrors( histoInvYieldPi0DalitzpPb5023GeV );
	graphInvYieldPi0DalitzpPb5023GeV->RemovePoint(0);
	graphInvYieldPi0DalitzpPb5023GeV->RemovePoint(0);
	
	TGraphAsymmErrors* graphInvYieldPi0DalitzpPb5023GeVXShifted =(TGraphAsymmErrors*) graphInvYieldPi0DalitzpPb5023GeV->Clone();
	
	
	
	
	 //cout<<"graphInvYieldPi0DalitzpPb5023GeV  Bin Shifted**************************************"<<endl;
	
	// TGraphAsymmErrors* graphInvYieldPi0DalitzpPb5023GeV = ( TGraphAsymmErrors* ) fNeutralPionDalitzResultspPbContainer->Get("Pi0SystErrorBinShifted");
	
	graphInvYieldPi0DalitzpPb5023GeVXShifted  = ApplyXshift(graphInvYieldPi0DalitzpPb5023GeVXShifted ,fitTsallisPi0PCMpPb5023GeVPtMult);
	
	
	graphInvYieldPi0DalitzpPb5023GeVXShifted->Print();
	
	
	
	TCanvas* cInvYieldPi0DalitzpPb5023GeVXShifted = new TCanvas("cInvYieldPi0DalitzpPb5023GeVXShifted","Pi0 Dalitz pPb at 5.23 TeV X shifted",200,10,700,500);
		
	DrawGammaCanvasSettings( cInvYieldPi0DalitzpPb5023GeVXShifted,  0.15, 0.02, 0.03, 0.1);
	cInvYieldPi0DalitzpPb5023GeVXShifted->SetLogx();
	cInvYieldPi0DalitzpPb5023GeVXShifted->SetLogy();
	
	TPad* padInvYieldPi0DalitzpPb5023GeVXShifted = new TPad("padInvYieldPi0DalitzpPb5023GeVXShifted", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padInvYieldPi0DalitzpPb5023GeVXShifted, 0.15, 0.02, 0.03, 0.1);
	padInvYieldPi0DalitzpPb5023GeVXShifted->Draw();
	
	
	TH2F * histoInvYieldPi0DalitzpPb5023GeVXShifted;
	histoInvYieldPi0DalitzpPb5023GeVXShifted = new TH2F("histoInvYieldPi0DalitzpPb5023GeVXShifted","histoInvYieldPi0DalitzpPb5023GeVXShifted",1000,0.23,30.,1000,2e-9,5e1);
	SetStyleHistoTH2ForGraphs(histoInvYieldPi0DalitzpPb5023GeVXShifted, "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoInvYieldPi0DalitzpPb5023GeVXShifted->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0DalitzpPb5023GeVXShifted,     	20, 1, kRed,   kRed);
	DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0DalitzpPb5023GeV,  	        20, 1, kBlue,  kBlue);
	
	
	
	graphInvYieldPi0DalitzpPb5023GeV->Draw("p,e1");
	graphInvYieldPi0DalitzpPb5023GeVXShifted->Draw("p,same,e1");
	//fitTsallisPi0DalitzpPb5023GeVPtMult->Draw("same");
	//fitTsallisPi0DalitzpPb5023GeVPtMult->Draw("same");
	
	TLegend* legendInvYieldPi0DalitzpPb5023GeVXShifted = new TLegend(0.20,0.25,0.70,0.40);
	legendInvYieldPi0DalitzpPb5023GeVXShifted->SetFillColor(0);
	legendInvYieldPi0DalitzpPb5023GeVXShifted->SetLineColor(0);
	legendInvYieldPi0DalitzpPb5023GeVXShifted->SetTextSize(0.03);
	legendInvYieldPi0DalitzpPb5023GeVXShifted->AddEntry(graphInvYieldPi0DalitzpPb5023GeV,"#pi^{0} pPb@5.02 TeV (Dalitz)");
	legendInvYieldPi0DalitzpPb5023GeVXShifted->AddEntry(graphInvYieldPi0DalitzpPb5023GeVXShifted,"#pi^{0} pPb@5.02 TeV XBin shifted (Dalitz)");
	//legendInvYieldPi0DalitzpPb5023GeVXShifted->AddEntry(fitTsallisPi0DalitzpPb5023GeVPtMult,"Tsallis Fit");
	legendInvYieldPi0DalitzpPb5023GeVXShifted->Draw("same");
	
		
	
	cInvYieldPi0DalitzpPb5023GeVXShifted->Update();
	cInvYieldPi0DalitzpPb5023GeVXShifted->SaveAs(Form("%s/InvYield_XShifted_Dalitz_pPb_5023GeV.%s",outputDir.Data(),suffix.Data()));
	
	

	
	
	//Extending the pp@5.023 TeV to the Dalitz pPb@5.023 binning
	
	cout<<"*******************Extending interpolated pp@5.023 TeV to the Dalitz pPb@5.023 binning ********************"<<endl;
	
	
	nPointsInter = graphInvYieldPi0DalitzpPb5023GeVXShifted->GetN();
	
	TGraphAsymmErrors* graphErrorsPi0DalitzBinPP5023GeVFitted = new TGraphAsymmErrors( nPointsInter );
	
	Double_t* xPointDal     = graphInvYieldPi0DalitzpPb5023GeVXShifted->GetX();
	Double_t* errorXlowDal  = graphInvYieldPi0DalitzpPb5023GeVXShifted->GetEXlow();
	Double_t* errorXhighDal = graphInvYieldPi0DalitzpPb5023GeVXShifted->GetEXhigh();
	
	for(Int_t iPoint = 0; iPoint < nPointsInter; iPoint++){
	  
	  //cout<<iPoint<<" "<<xPoint[iPoint] <<" "<<errorXlow[iPoint]<<" "<<errorXhigh[iPoint]<<endl;
	  
	  //Double_t ptMin = xPointDal[iPoint] - errorXlowDal[iPoint];
	  //Double_t ptMax = xPointDal[iPoint] + errorXhighDal[iPoint];
	  //Double_t binWidth = ptMax - ptMin;
														   
	  //Double_t yInt      = fitTsallisPi0PCMPP5023GeV->Integral(ptMin,ptMax,resultfitTsallisPi0PCMPP5023GeV->GetParams()) / binWidth;
	 // Double_t yIntError = fitTsallisPi0PCMPP5023GeV->IntegralError(ptMin,ptMax, resultfitTsallisPi0PCMPP5023GeV->GetParams(), resultfitTsallisPi0PCMPP5023GeV->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
			
	  
	  
	  Double_t point[1] = {xPointDal[iPoint]};
	  Double_t   err[1];
	  Double_t yInt = fitTsallisPi0PCMPP5023GeV->Eval(point[0]);
	  resultfitTsallisPi0PCMPP5023GeV->GetConfidenceIntervals(1, 1, 1, point, err, 0.95, false);  
	  
	  Double_t yIntError = err[0];
	  
	 
	   cout<<"Point "<<iPoint<<" xVal "<<xPointDal[iPoint] <<" xErrlow "<<errorXlowDal[iPoint]<<" xErrhigh "<<errorXhighDal<<" yVal "<<yInt<<" yErrlow "<<yIntError<<" yErrhigh "<<yIntError<<endl;
	
	  
	  graphErrorsPi0DalitzBinPP5023GeVFitted->SetPoint(iPoint,xPointDal[iPoint],yInt);
	  graphErrorsPi0DalitzBinPP5023GeVFitted->SetPointError(iPoint,errorXlowDal[iPoint],errorXhighDal[iPoint],yIntError,yIntError);
	  
	  
  
	}
	
	
	//Computing RpPb using the interpolated pp@5.023 TeV and Dalitz pPb@5.023 TeV spectra.
	
	
	TGraphAsymmErrors* graphDalitzRpPb = new TGraphAsymmErrors( nPointsInter);
	
	Double_t *xBinsDalitzpPb    = graphInvYieldPi0DalitzpPb5023GeVXShifted->GetX();
	Double_t *xBinsDalitzPP     = graphErrorsPi0DalitzBinPP5023GeVFitted->GetX();
	Double_t *yBinsDalitzpPb    = graphInvYieldPi0DalitzpPb5023GeVXShifted->GetY();
	Double_t *yBinsDalitzPP     = graphErrorsPi0DalitzBinPP5023GeVFitted->GetY();
	
	
	Double_t *xBinsErrDalitzpPblow  = graphInvYieldPi0DalitzpPb5023GeVXShifted->GetEXlow();
	//Double_t *xBinsErrDalitzpPbhigh = graphInvYieldPi0DalitzpPb5023GeVXShifted->GetEXhigh();
	
	Double_t *yBinsErrDalitzpPblow = graphInvYieldPi0DalitzpPb5023GeV->GetEYlow();
	//Double_t *xBinsErrDalitzPP  = graphErrorsPi0DalitzBinPP5023GeVFitted->GetEX();
	Double_t *yBinsErrDalitzPPlow  = graphErrorsPi0DalitzBinPP5023GeVFitted->GetEYlow();
	
	cout<<graphInvYieldPi0DalitzpPb5023GeV->GetN()<<" "<<graphErrorsPi0DalitzBinPP5023GeVFitted->GetN()<<endl;
	
	
	
	
	for(Int_t iPoint = 0; iPoint < nPointsInter; iPoint++){
	  
	 graphDalitzRpPb->SetPoint( iPoint, xBinsDalitzpPb[iPoint], yBinsDalitzpPb[iPoint] /(fNcoll* yBinsDalitzPP[iPoint]) );
	 
	 Double_t errYStat = pow( pow( yBinsErrDalitzpPblow[iPoint]/yBinsDalitzpPb[iPoint],2. ) + pow( yBinsErrDalitzPPlow[iPoint]/yBinsDalitzPP[iPoint],2.), 0.5)*yBinsDalitzpPb[iPoint] /(fNcoll* yBinsDalitzPP[iPoint]);
	 
	 graphDalitzRpPb->SetPointError(iPoint, xBinsErrDalitzpPblow[iPoint], xBinsErrDalitzpPblow[iPoint], errYStat,errYStat);
		 
	 
	}
	
	
	TCanvas* canvasPi0InvYieldPP5023GeVPCMBinning = new TCanvas("canvasPi0InvYieldPP5023GeVPCMBinning","PP reference with PCM binning and fit",200,10,700,500);
		
	DrawGammaCanvasSettings( canvasPi0InvYieldPP5023GeVPCMBinning,  0.15, 0.02, 0.03, 0.1);
	canvasPi0InvYieldPP5023GeVPCMBinning->SetLogx();
	canvasPi0InvYieldPP5023GeVPCMBinning->SetLogy();
	
	TPad* padPi0InvYieldPP5023GeVPCMBinning = new TPad("padPi0InvYieldPP5023GeVPCMBinning", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padPi0InvYieldPP5023GeVPCMBinning, 0.15, 0.02, 0.03, 0.1);
	padPi0InvYieldPP5023GeVPCMBinning->Draw();
	
	
	TH2F * histoPi0InvYieldPP5023GeVPCMBinning;
	histoPi0InvYieldPP5023GeVPCMBinning = new TH2F("histoPi0InvYieldPP5023GeVPCMBinning","histoPi0InvYieldPP5023GeVPCMBinning",1000,0.23,30.,1000,2e-9,5e1);
	SetStyleHistoTH2ForGraphs(histoPi0InvYieldPP5023GeVPCMBinning, "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoPi0InvYieldPP5023GeVPCMBinning->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphErrorsPi0PCMPP5023GeVFitted,     20, 1, kRed,   kRed);
	DrawGammaSetMarkerTGraphErr(graphErrosInterPolation5023GeV2TPoints, 	20, 1, kGreen, kGreen);
	//DrawGammaSetMarkerTGraphErr(fitTsallisPi0PCMPP5023GeV, 	        20, 1, kRed, kRed);
	graphErrosInterPolation5023GeV2TPoints->SetFillColor(0);
	

	graphErrorsPi0PCMPP5023GeVFitted->Draw("p,same,e1");
	graphErrosInterPolation5023GeV2TPoints->Draw("p,same,e1");
	//fitTsallisPi0PCMPP5023GeV->Draw("same");
	
	
	
	TLatex *textPi0InvYieldPP5023GeVPCMBinning = new TLatex(0.25,0.43,"pp reference  #sqrt{#it{s}} = 5.023 TeV");
	SetStyleTLatex( textPi0InvYieldPP5023GeVPCMBinning,0.035,4); 
	textPi0InvYieldPP5023GeVPCMBinning->Draw();
	
	
	graphErrorsPi0PCMPP5023GeVFitted->SetFillColor(0);
	
	TLegend* legendPi0InvYieldPP5023GeVPCMBinning = new TLegend(0.20,0.25,0.70,0.40);
	legendPi0InvYieldPP5023GeVPCMBinning->SetFillColor(0);
	legendPi0InvYieldPP5023GeVPCMBinning->SetLineColor(0);
	legendPi0InvYieldPP5023GeVPCMBinning->SetTextSize(0.03);
	legendPi0InvYieldPP5023GeVPCMBinning->AddEntry(graphErrorsPi0PCMPP5023GeVFitted,"#pi^{0} #rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-} (PCM) binning","pef");
	legendPi0InvYieldPP5023GeVPCMBinning->AddEntry(graphErrosInterPolation5023GeV2TPoints,"Interpolated pp reference at 5.023");
	legendPi0InvYieldPP5023GeVPCMBinning->AddEntry(fitTsallisPi0PCMPP5023GeV,"Tsallis Fit");
	
	
	legendPi0InvYieldPP5023GeVPCMBinning->Draw();
	
	canvasPi0InvYieldPP5023GeVPCMBinning->Update();
	canvasPi0InvYieldPP5023GeVPCMBinning->SaveAs(Form("%s/PP_Reference_PCM_Binning_and_Fit.%s",outputDir.Data(),suffix.Data()));
	
	
	///////////////////////////////////
	
	
	TCanvas* canvasPi0InvYieldPP5023GeVDalitzBinning = new TCanvas("canvasPi0InvYieldPP5023GeVDalitzBinning","PP reference with Dalitz binning and fit",200,10,700,500);
		
	DrawGammaCanvasSettings( canvasPi0InvYieldPP5023GeVDalitzBinning,  0.15, 0.02, 0.03, 0.1);
	canvasPi0InvYieldPP5023GeVDalitzBinning->SetLogx();
	canvasPi0InvYieldPP5023GeVDalitzBinning->SetLogy();
	
	TPad* padPi0InvYieldPP5023GeVDalitzBinning = new TPad("padPi0InvYieldPP5023GeVDalitzBinning", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padPi0InvYieldPP5023GeVDalitzBinning, 0.15, 0.02, 0.03, 0.1);
	padPi0InvYieldPP5023GeVDalitzBinning->Draw();
	
	
	TH2F * histoPi0InvYieldPP5023GeVDalitzBinning;
	histoPi0InvYieldPP5023GeVDalitzBinning = new TH2F("histoPi0InvYieldPP5023GeVDalitzBinning","histoPi0InvYieldPP5023GeVDalitzBinning",1000,0.23,30.,1000,2e-9,5e1);
	SetStyleHistoTH2ForGraphs(histoPi0InvYieldPP5023GeVDalitzBinning, "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoPi0InvYieldPP5023GeVDalitzBinning->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphErrorsPi0DalitzBinPP5023GeVFitted,     20, 1, kRed,   kRed);
	//DrawGammaSetMarkerTGraphErr(fitTsallisPi0DalitzPP5023GeV, 	        20, 1, kRed, kRed);
	

	graphErrorsPi0DalitzBinPP5023GeVFitted->Draw("p,same,e1");
	fitTsallisPi0PCMPP5023GeV->Draw("same");
	graphErrosInterPolation5023GeV2TPoints->Draw("p,same,e1");
	
	
	
	TLatex *textPi0InvYieldPP5023GeVDalitzBinning = new TLatex(0.25,0.40,"pp reference  #sqrt{#it{s}} = 5.023 TeV");
	SetStyleTLatex( textPi0InvYieldPP5023GeVDalitzBinning,0.035,4); 
	textPi0InvYieldPP5023GeVDalitzBinning->Draw();
	
	
	//graphErrorsPi0DalitzPP5023GeVFitted->SetFillColor(0);
	graphErrorsPi0DalitzBinPP5023GeVFitted->SetFillColor(0);
	TLegend* legendPi0InvYieldPP5023GeVDalitzBinning = new TLegend(0.20,0.25,0.70,0.40);
	legendPi0InvYieldPP5023GeVDalitzBinning->SetFillColor(0);
	legendPi0InvYieldPP5023GeVDalitzBinning->SetLineColor(0);
	legendPi0InvYieldPP5023GeVDalitzBinning->SetTextSize(0.03);
	legendPi0InvYieldPP5023GeVDalitzBinning->AddEntry(graphErrorsPi0DalitzBinPP5023GeVFitted,"#pi^{0} #rightarrow e^{+}e^{-}#gamma #rightarrow e^{+}e^{-}e^{+}e^{-} (Dalitz) binning","pef");
	legendPi0InvYieldPP5023GeVDalitzBinning->AddEntry(graphErrosInterPolation5023GeV2TPoints,"Interpolated pp reference at 5.023");
	legendPi0InvYieldPP5023GeVDalitzBinning->AddEntry(fitTsallisPi0PCMPP5023GeV,"Tsallis Fit");
	
	
	legendPi0InvYieldPP5023GeVDalitzBinning->Draw();
	
	canvasPi0InvYieldPP5023GeVDalitzBinning->Update();
	
	
	canvasPi0InvYieldPP5023GeVDalitzBinning->SaveAs(Form("%s/PP_Reference_Dalitz_Binning_and_Fit.%s",outputDir.Data(),suffix.Data()));
	
	
	
	
	////////////////////////////////////////////////////////////
	
	
	TCanvas* canvasPi0InvYieldPCMIterPP5023GeV = new TCanvas("canvasPi0InvYieldPCMIterPP5023GeV","PCM Interpolation pp@5.023 TeV",200,10,700,500);
		
	DrawGammaCanvasSettings( canvasPi0InvYieldPCMIterPP5023GeV,  0.15, 0.02, 0.03, 0.1);
	canvasPi0InvYieldPCMIterPP5023GeV->SetLogx();
	canvasPi0InvYieldPCMIterPP5023GeV->SetLogy();
	
	TPad* padInvYieldPi0PCM = new TPad("padInvYieldPi0PCM", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padInvYieldPi0PCM, 0.15, 0.02, 0.03, 0.1);
	padInvYieldPi0PCM->Draw();
	
	
	TH2F * histoInvYieldPi0PCM;
	histoInvYieldPi0PCM = new TH2F("histoInvYieldPi0PCM","histoInvYieldPi0PCM",1000,0.23,30.,1000,2e-9,5e1);
	SetStyleHistoTH2ForGraphs(histoInvYieldPi0PCM, "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoInvYieldPi0PCM->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0PCMpPb5023GeVXShifted,     	20, 1, kBlue,   kBlue);
	DrawGammaSetMarkerTGraphAsym(graphErrorsPi0PCMPP5023GeVFitted, 	        20, 1, kRed, kRed);
	
		
	
	graphInvYieldPi0PCMpPb5023GeVXShifted->Draw("p,same,e1");
	graphErrorsPi0PCMPP5023GeVFitted->Draw("p,same,e1");
	canvasPi0InvYieldPCMIterPP5023GeV->Update();
	
	TLatex *textInvPPPCM = new TLatex(0.25,0.40,"p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV, (PCM)");
	SetStyleTLatex( textInvPPPCM,0.035,4); 
	textInvPPPCM->Draw();
	
	//graphDalitzRpPb->SetFillColor(0);
	
	graphErrorsPi0PCMPP5023GeVFitted->SetFillColor(0);
	graphErrorsPi0DalitzBinPP5023GeVFitted->SetFillColor(0);
	TLegend* legendInvYiedlPPPCM = new TLegend(0.20,0.30,0.70,0.35);
	legendInvYiedlPPPCM->SetFillColor(0);
	legendInvYiedlPPPCM->SetLineColor(0);
	//legendInvYiedlPPPCM->SetNColumns(3);
	legendInvYiedlPPPCM->SetTextSize(0.03);
	legendInvYiedlPPPCM->AddEntry(graphErrorsPi0PCMPP5023GeVFitted,"PP reference, |Y| < 0.8","pef");
	legendInvYiedlPPPCM->AddEntry(graphInvYieldPi0PCMpPb5023GeVXShifted,"|Y| < 0.4","pef");
	
	
	legendInvYiedlPPPCM->Draw();
	
	canvasPi0InvYieldPCMIterPP5023GeV->SaveAs(Form("%s/Reference_PP_PCM.%s",outputDir.Data(),suffix.Data()));
	
	
	
	TCanvas* canvasPi0InvYieldDalitzIterPP5023GeV = new TCanvas("canvasPi0InvYieldDalitzIterPP5023GeV","Dalitz Interpolation pp@5.023 TeV",200,10,700,500);
		
	DrawGammaCanvasSettings( canvasPi0InvYieldDalitzIterPP5023GeV,  0.15, 0.02, 0.03, 0.1);
	canvasPi0InvYieldDalitzIterPP5023GeV->SetLogx();
	canvasPi0InvYieldDalitzIterPP5023GeV->SetLogy();
	
	TPad* padInvYieldPi0Dalitz = new TPad("padInvYieldPi0Dalitz", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padInvYieldPi0Dalitz, 0.15, 0.02, 0.03, 0.1);
	padInvYieldPi0Dalitz->Draw();
	
	
	TH2F * histoInvYieldPi0Dalitz;
	histoInvYieldPi0Dalitz = new TH2F("histoInvYieldPi0Dalitz","histoInvYieldPi0Dalitz",1000,0.23,30.,1000,2e-9,5e1);
	SetStyleHistoTH2ForGraphs(histoInvYieldPi0Dalitz, "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoInvYieldPi0Dalitz->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0DalitzpPb5023GeVXShifted,     	20, 1, kBlue,   kBlue);
	DrawGammaSetMarkerTGraphAsym(graphErrorsPi0DalitzBinPP5023GeVFitted, 	        20, 1, kRed, kRed);
	
	//graphErrorsPi0DalitzPP5023GeVFitted->SetFillColor(2);
	
	
	graphInvYieldPi0DalitzpPb5023GeV->Draw("p,same,e1");
	graphErrorsPi0DalitzBinPP5023GeVFitted->Draw("p,same,e1");
	//graphErrorsInvYieldPi0DalitzStat2760GeVFittedScale->Draw("p,same,e1");
	canvasPi0InvYieldDalitzIterPP5023GeV->Update();
	
	TLatex *textInvPPDalitz = new TLatex(0.25,0.40,"p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV, (Dalitz)");
	SetStyleTLatex( textInvPPDalitz,0.035,4); 
	textInvPPDalitz->Draw();
	
	//graphDalitzRpPb->SetFillColor(0);
	
	//graphErrorsPi0DalitzBinPP5023GeVFitted->RemovePoint(0);
	//graphErrorsPi0DalitzBinPP5023GeVFitted->RemovePoint(0);
	
	//graphInvYieldPi0DalitzpPb5023GeVXShifted->RemovePoint(0);
	//graphInvYieldPi0DalitzpPb5023GeVXShifted->RemovePoint(0);
	
	
	
	graphInvYieldPi0DalitzpPb5023GeVXShifted->SetFillColor(0);
	graphErrorsPi0DalitzBinPP5023GeVFitted->SetFillColor(0);
	TLegend* legendInvYiedlPPDalitz = new TLegend(0.20,0.30,0.70,0.35);
	legendInvYiedlPPDalitz->SetFillColor(0);
	legendInvYiedlPPDalitz->SetLineColor(0);
	//legendInvYiedlPPDalitz->SetNColumns(3);
	legendInvYiedlPPDalitz->SetTextSize(0.03);
	legendInvYiedlPPDalitz->AddEntry(graphErrorsPi0DalitzBinPP5023GeVFitted,"PP reference, |Y| < 0.8","pef");
	legendInvYiedlPPDalitz->AddEntry(graphInvYieldPi0DalitzpPb5023GeVXShifted,"|Y| < 0.4","pef");
	
	
	legendInvYiedlPPDalitz->Draw();
	
	canvasPi0InvYieldDalitzIterPP5023GeV->SaveAs(Form("%s/Reference_PP_Dalitz.%s",outputDir.Data(),suffix.Data()));
	
	
	
	
	
	
	////////////////////////////////////////////////////////////
	
	//	TGraphErrors*	graphRpPbPCM_2=RemovePointsFromGraph(graphRpPbPCM,2);		
	

	graphRpPbPCM->RemovePoint(0);	
	graphRpPbPCM->RemovePoint(0);	
	
	TCanvas* canvasRpPbPCM = new TCanvas("canvasRpPbPCM","R_{pPb} PCM",200,10,700,500);
	
	DrawGammaCanvasSettings( canvasRpPbPCM,  0.1, 0.01, 0.015, 0.13);
	
	TH2F * histo2DRpPbPCM = new TH2F("histo2DRpPbPCM","histo2DRpPbPCM",1000,0.,15.,1000,0.3,1.5);
	SetStyleHistoTH2ForGraphs(histo2DRpPbPCM, "p_{T} (GeV/c)","R_{pPb}", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 
	histo2DRpPbPCM->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphRpPbPCM,   	20, 1, kRed+1,   kRed+1);
 	graphRpPbPCM->Draw("p,same,e");
	graphRpPbPCM->SetFillColor(0);
	DrawGammaLines(0., 15.,1., 1.,2.,kBlack);
	
	TLatex *textRpPbPCM = new TLatex(0.3,0.35,"p-Pb minimum bias #sqrt{#it{s}_{_{NN}}} = 5.02 TeV");
	SetStyleTLatex( textRpPbPCM,0.035,4); 
	textRpPbPCM->Draw();
	
	TLatex *LabelpPb2 = new TLatex(0.3,0.4,"ALICE work in progress");
	SetStyleTLatex( LabelpPb2, 0.035,4);
	LabelpPb2->Draw();
	
	TLegend* legendRpPbPCM = new TLegend(0.3,0.25,0.85,0.3);
	legendRpPbPCM->SetFillColor(0);
	legendRpPbPCM->SetLineColor(0);
	legendRpPbPCM->SetNColumns(3);
	legendRpPbPCM->SetTextSize(0.04);
	legendRpPbPCM->AddEntry(graphRpPbPCM,"#pi^{0} #rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-} (PCM)","pef");
	
	legendRpPbPCM->Draw();
	
	

	canvasRpPbPCM->Update(); 
	
	canvasRpPbPCM->SaveAs(Form("%s/RpPb_PCM.%s",outputDir.Data(),suffix.Data()));
	
	
	TCanvas* canvasRpPbDalitz = new TCanvas("canvasRpPbDalitz","R_{pPb} Dalitz",200,10,1200,700);
	
	DrawGammaCanvasSettings( canvasRpPbDalitz,  0.1, 0.01, 0.015, 0.13);
	
	TH2F * histo2DRpPbDalitz = new TH2F("histo2DRpPbDalitz","histo2DRpPbDalitz",1000,0.,15.,1000,0.3,1.5);
	SetStyleHistoTH2ForGraphs(histo2DRpPbDalitz, "p_{T} (GeV/c)","R_{pPb}", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 
	histo2DRpPbDalitz->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphDalitzRpPb,     	20, 1, kBlue,   kBlue);
 	graphDalitzRpPb->Draw("p,same,e");
	
	DrawGammaLines(0., 15.,1., 1.,2.,kBlack);
	
	TLatex *textRpPbDalitz = new TLatex(0.25,0.40,"p-Pb minimum bias #sqrt{#it{s}_{_{NN}}} = 5.02 TeV");
	SetStyleTLatex( textRpPbDalitz,0.035,4); 
	textRpPbDalitz->Draw();
	
	graphDalitzRpPb->SetFillColor(0);
	TLegend* legendRpPbDalitz = new TLegend(0.23,0.30,0.85,0.35);
	legendRpPbDalitz->SetFillColor(0);
	legendRpPbDalitz->SetLineColor(0);
	legendRpPbDalitz->SetNColumns(3);
	legendRpPbDalitz->SetTextSize(0.04);
	legendRpPbDalitz->AddEntry(graphDalitzRpPb,"#pi^{0} #rightarrow e^{+}e^{-}#gamma #rightarrow e^{+}e^{-}e^{+}e^{-} (Dalitz)","pef");
	
	legendRpPbDalitz->Draw();
	
	canvasRpPbDalitz->Update();
	
	canvasRpPbDalitz->SaveAs(Form("%s/RpPb_Dalitz.%s",outputDir.Data(),suffix.Data()));
	

	
	
	TCanvas* canvasRpPbCombine = new TCanvas("canvasRpPbCombine","Combined R_{pPb}",200,10,1200,700);
	
	DrawGammaCanvasSettings( canvasRpPbCombine,  0.1, 0.01, 0.015, 0.13);
	
	TH2F * histo2DRpPbCombine = new TH2F("histo2DRpPbCombine","histo2DRpPbCombine",1000,0.,15.,1000,0.3,1.5);
	SetStyleHistoTH2ForGraphs(histo2DRpPbCombine, "p_{T} (GeV/c)","R^{#pi^{0}}_{pPb}", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 
	histo2DRpPbCombine->DrawCopy(); 
	
	//DrawGammaSetMarkerTGraphErr(graphDalitzRpPb,     	20, 1, kBlue,   kBlue);
 	graphDalitzRpPb->Draw("p,same,e");
	graphRpPbPCM->Draw("p,same,e");
	
	DrawGammaLines(0., 15.,1., 1.,2.,kBlack);
	
	TLatex *textRpPbCombine = new TLatex(0.25,0.40,"p-Pb minimum bias #sqrt{#it{s}_{_{NN}}} = 5.02 TeV");
	SetStyleTLatex( textRpPbCombine,0.035,4); 
	textRpPbCombine->Draw();
	
	
	//graphDalitzRpPb->SetFillColor(0);
	//graphRpPbPCM->SetFillColor(0);
	
	TLegend* legendRpPbCombine = new TLegend(0.23,0.30,0.85,0.35);
	legendRpPbCombine->SetFillColor(0);
	legendRpPbCombine->SetLineColor(0);
	//egendRpPbCombine->SetNColumns(3);
	legendRpPbCombine->SetTextSize(0.03);
	//legendRpPbCombine->SetEntrySeparation(0.02);
	legendRpPbCombine->AddEntry(graphDalitzRpPb,"#pi^{0} #rightarrow e^{+}e^{-}#gamma #rightarrow e^{+}e^{-}e^{+}e^{-} (Dalitz)","pef");
	legendRpPbCombine->AddEntry(graphRpPbPCM,"#pi^{0} #rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-} (PCM)","pef");
	
	legendRpPbCombine->Draw();
	
	
	canvasRpPbCombine->SaveAs(Form("%s/RpPb_Combined_PCM_Dalitz.%s",outputDir.Data(),suffix.Data()));
	
		
	
	
	TCanvas* canvasRpPbChargedParCombined = new TCanvas("canvasRpPbChargedParCombined","Combined R_{pPb} with charged particles",200,10,1200,700);
	
	TGraphAsymmErrors* graphAsymmChargedParticlesRpPb = GetChargeParticlesRpPb();
	
	DrawGammaCanvasSettings( canvasRpPbChargedParCombined,  0.1, 0.01, 0.015, 0.13);
	
	TH2F * histo2DRpPbChargedParCombined = new TH2F("histo2DRpPbChargedParCombined","histo2DRpPbChargedParCombined",1000,0.,15.,1000,0.3,1.5);
	SetStyleHistoTH2ForGraphs(histo2DRpPbChargedParCombined, "p_{T} (GeV/c)","R_{pPb}", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 
	histo2DRpPbChargedParCombined->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphAsymmChargedParticlesRpPb, 20, 1, kGreen,   kGreen);
	graphAsymmChargedParticlesRpPb->SetFillColor(0);
	
	//	graphDalitzRpPb->Draw("p,same,e");
	graphRpPbPCM->Draw("p,same,e");
	graphAsymmChargedParticlesRpPb->Draw("p,same,e");
	DrawGammaLines(0., 15.,1., 1.,2.,kBlack);
	
	TLatex *textRpPbCombined = new TLatex(0.25,0.40,"p-Pb minimum bias #sqrt{#it{s}_{_{NN}}} = 5.02 TeV");
	SetStyleTLatex( textRpPbCombined,0.035,4); 
	textRpPbCombined->Draw();
	
	

	
	TLegend* legendRpPbCombined = new TLegend(0.23,0.20,0.85,0.35);
	legendRpPbCombined->SetFillColor(0);
	legendRpPbCombined->SetLineColor(0);
	//egendRpPbCombine->SetNColumns(3);
	legendRpPbCombined->SetTextSize(0.03);
	//legendRpPbCombine->SetEntrySeparation(0.02);
	//	legendRpPbCombined->AddEntry(graphDalitzRpPb,"#pi^{0} #rightarrow e^{+}e^{-}#gamma #rightarrow e^{+}e^{-}e^{+}e^{-} (Dalitz)","pef");
	legendRpPbCombined->AddEntry(graphRpPbPCM,"#pi^{0} #rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-} (PCM)","pef");
	legendRpPbCombined->AddEntry(graphAsymmChargedParticlesRpPb,"ALICE,NSD,Charged particles,|#eta_{cms}| < 0.3","pef");
	
	legendRpPbCombined->Draw();
	
	
	/*TLegend* legendRpPbChargedPar = new TLegend(0.65,0.10,0.85,0.31);
	legendRpPbChargedPar->SetFillColor(0);
	legendRpPbChargedPar->SetLineColor(0);
	//egendRpPbCombine->SetNColumns(3);
	legendRpPbChargedPar->SetTextSize(0.03);
	//legendRpPbCombine->SetEntrySeparation(0.02);
	legendRpPbChargedPar->AddEntry(graphAsymmChargedParticlesRpPb,"ALICE,NSD,Charged particles,|#eta_{cms}| < 0.3","pef");
	
	
	legendRpPbChargedPar->Draw("same");*/
	
	
	
	
	
	
	canvasRpPbChargedParCombined->SaveAs(Form("%s/RpPb_Combined_PCM_ChargedParticles.%s",outputDir.Data(),suffix.Data()));
	
	
	
	
	
	
	
	
	
									    //  R_pPb5000_y0_eps09s_fdss_mb.dat
        const char *fileNameEPS09sPi0AKK = "ExternalInputpPb/R_pi0_pPb_y0_eps09s/R_pPb5000_y0_eps09s_akk_mb.dat";
	const char *fileNameEPS09sPi0DSS = "ExternalInputpPb/R_pi0_pPb_y0_eps09s/R_pPb5000_y0_eps09s_fdss_mb.dat";
	const char *fileNameEPS09sPi0KKP = "ExternalInputpPb/R_pi0_pPb_y0_eps09s/R_pPb5000_y0_eps09s_kkp_mb.dat";
	
	const char *fileNameEPS09sPi0CGC = "ExternalInputpPb/ColorGlassCondensate.dat";
	
	//****************************** extracting EPS09s predictions**************************************
	
	ifstream 		inDSS;
	ifstream                inAKK;
	ifstream                inKKP;
	ifstream 		inCGC;	
	
	
	Int_t nlinesEPSsPi0fDSS = 0;
	
	inDSS.open(fileNameEPS09sPi0DSS,ios_base::in);
	
	Double_t xEPSsPi0fDSS[100],yEPSsPi0fDSS[100];
	
	Double_t xUpErrorEPSsPi0DSS[100],xDownErrorEPSsPi0DSS[100];
	Double_t yUpErrorEPSsPi0DSS[100],yDownErrorEPSsPi0DSS[100];
	
	
	while(!inDSS.eof()){
			nlinesEPSsPi0fDSS++;
			inDSS >> xEPSsPi0fDSS[nlinesEPSsPi0fDSS]  >> yEPSsPi0fDSS[nlinesEPSsPi0fDSS] >> yUpErrorEPSsPi0DSS[nlinesEPSsPi0fDSS]>>yDownErrorEPSsPi0DSS[nlinesEPSsPi0fDSS];
			yUpErrorEPSsPi0DSS[nlinesEPSsPi0fDSS] = ( yUpErrorEPSsPi0DSS[nlinesEPSsPi0fDSS] - yEPSsPi0fDSS[nlinesEPSsPi0fDSS] ) /yEPSsPi0fDSS[nlinesEPSsPi0fDSS];
			yDownErrorEPSsPi0DSS[nlinesEPSsPi0fDSS] = -1*( yDownErrorEPSsPi0DSS[nlinesEPSsPi0fDSS] - yEPSsPi0fDSS[nlinesEPSsPi0fDSS] ) /yEPSsPi0fDSS[nlinesEPSsPi0fDSS];
			xUpErrorEPSsPi0DSS[nlinesEPSsPi0fDSS]   = 0;
			xDownErrorEPSsPi0DSS[nlinesEPSsPi0fDSS] = 0;
			
			    
			cout << nlinesEPSsPi0fDSS << "         "  << xEPSsPi0fDSS[nlinesEPSsPi0fDSS] << "         "  <<yEPSsPi0fDSS[nlinesEPSsPi0fDSS]<<"         "<<yUpErrorEPSsPi0DSS[nlinesEPSsPi0fDSS]<<"          "<<yDownErrorEPSsPi0DSS[nlinesEPSsPi0fDSS]<< endl;
	
	}
	inDSS.close();
	
	TGraph* graphPi0DSS5000 = new TGraph(nlinesEPSsPi0fDSS,xEPSsPi0fDSS,yEPSsPi0fDSS);
	TGraphAsymmErrors* graphAsymmErrorsPi0DSS5000 = new TGraphAsymmErrors(nlinesEPSsPi0fDSS,xEPSsPi0fDSS,yEPSsPi0fDSS,xDownErrorEPSsPi0DSS,xUpErrorEPSsPi0DSS,yDownErrorEPSsPi0DSS,yUpErrorEPSsPi0DSS);
	
	
	
	
	Int_t nlinesPi0AKK = 0;
		
	inAKK.open(fileNameEPS09sPi0AKK,ios_base::in);
	
	Double_t xESPsPi0AKK[100], yESPsPi0AKK[100];
	
		
	while(!inAKK.eof()){
			nlinesPi0AKK++;
			inAKK >> xESPsPi0AKK[nlinesPi0AKK]  >> yESPsPi0AKK[nlinesPi0AKK];
			cout << nlinesPi0AKK << "         "  << xESPsPi0AKK[nlinesPi0AKK] << "         "  <<xESPsPi0AKK[nlinesPi0AKK]<<endl;
	
	}
	inAKK.close();
	
	TGraph* graphPi0ESP09sPi0AKK = new TGraph(nlinesPi0AKK,xESPsPi0AKK,yESPsPi0AKK);	
	
	
	Int_t nlinesPi0KKP = 0;
		
	inKKP.open(fileNameEPS09sPi0KKP,ios_base::in);
	
	Double_t xESPsPi0KKP[100], yESPsPi0KKP[100];
	
		
	while(!inKKP.eof()){
			nlinesPi0KKP++;
			inKKP >> xESPsPi0KKP[nlinesPi0KKP]  >> yESPsPi0KKP[nlinesPi0KKP];
			cout << nlinesPi0KKP << "         "  << xESPsPi0KKP[nlinesPi0KKP] << "         "  <<xESPsPi0KKP[nlinesPi0KKP]<<endl;
	
	}
	inKKP.close();
	
	TGraph* graphPi0ESP09sPi0KKP = new TGraph(nlinesPi0KKP,xESPsPi0KKP,yESPsPi0KKP);	
	
	
	
	
	
	
	
	
	TCanvas* canvasRpPbEPS09s = new TCanvas("canvasRpPbEPS09s","Combined predictions for R^{#pi0}_{pPb}",200,10,1200,700);
	
	
	DrawGammaCanvasSettings( canvasRpPbEPS09s,  0.1, 0.01, 0.015, 0.13);
	
	
	/*TPad* padcanvasRpPbEPS09s = new TPad("padcanvasRpPbEPS09s", "", 0., 0.05, 1.0, 1.0,0.9, -1, -2);
	DrawGammaPadSettings( padcanvasRpPbEPS09s, 0.09, 0.02, 0.03, 0.12);
	padcanvasRpPbEPS09s->SetFillStyle(3001);
	padcanvasRpPbEPS09s->Draw();
	
	padcanvasRpPbEPS09s->cd();*/
	
	
	
	TH2F * histo2DRpPbESP09sCompared = new TH2F("histo2DRpPbESP09sCompared","histo2DRpPbESP09sCompared",1500,0.,15.,1000,0.3,1.5);
	SetStyleHistoTH2ForGraphs(histo2DRpPbESP09sCompared, "p_{T} (GeV/c)","R^{#pi^{0}}_{pPb}", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 
	histo2DRpPbESP09sCompared->DrawCopy(); 
	
	graphPi0ESP09sPi0KKP->RemovePoint(0);
	graphPi0ESP09sPi0AKK->RemovePoint(0);
	graphPi0DSS5000->RemovePoint(0);
	graphAsymmErrorsPi0DSS5000->RemovePoint(0);
	
	graphPi0ESP09sPi0KKP->SetLineColor(kRed);
	graphPi0ESP09sPi0KKP->SetLineWidth(3);
	graphPi0ESP09sPi0KKP->SetLineStyle(10);
	
	
	//graphPi0ESP09sPi0KKP->RemovePoint(0);
	graphPi0ESP09sPi0AKK->SetLineColor(kBlue);
	graphPi0ESP09sPi0AKK->SetLineWidth(3);
	graphPi0ESP09sPi0AKK->SetLineStyle(7);
	
	
	
	graphPi0DSS5000->SetLineColor(8);
	graphPi0DSS5000->SetLineWidth(3);
	
	
	
	
	DrawGammaLines(0., 15.,1., 1.,2.,kBlack);
	graphAsymmErrorsPi0DSS5000->SetFillColor(kBlue);
	graphAsymmErrorsPi0DSS5000->SetFillStyle(3018);
	
	
	//	
	
		
	
        TLatex *textRpPbCompared = new TLatex(0.25,0.40,"p-Pb min. bias #sqrt{#it{s}_{_{NN}}} = 5.02 TeV");
	SetStyleTLatex( textRpPbCompared,0.035,4); 
	
	TLatex *textRpPbEPS09s = new TLatex(0.65,0.40,"p-Pb min. bias #sqrt{#it{s}_{_{NN}}} = 5.0 TeV");
	SetStyleTLatex( textRpPbEPS09s,0.035,4); 
	
	TLatex *textRpPbEPS09sPre = new TLatex(0.65,0.36,"Model calculations");
	SetStyleTLatex( textRpPbEPS09sPre,0.035,4); 
	
	//textRpPbCompared->Draw();
	
	
	//graphDalitzRpPb->SetFillColor(0);
	//graphRpPbPCM->SetFillColor(0);
	
	TLegend* legendRpPbCompared = new TLegend(0.23,0.28,0.85,0.35);
	legendRpPbCompared->SetFillColor(0);
	legendRpPbCompared->SetLineColor(0);
	//egendRpPbCombine->SetNColumns(3);
	legendRpPbCompared->SetTextSize(0.03);
	//legendRpPbCombine->SetEntrySeparation(0.02);
	//	legendRpPbCompared->AddEntry(graphDalitzRpPb,"#pi^{0} #rightarrow e^{+}e^{-}#gamma #rightarrow e^{+}e^{-}e^{+}e^{-} (Dalitz)","pef");
	legendRpPbCompared->AddEntry(graphRpPbPCM,"#pi^{0} #rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-} (PCM)","pef");
	
	TLatex *LabelpPb3 = new TLatex(0.25,0.45,"ALICE work in progress");
	SetStyleTLatex( LabelpPb3, 0.035,4);
	//	LabelpPb3->Draw("same");	
	
	
	TLegend* legendRpPbPred = new TLegend(0.65,0.20,0.85,0.31);
	legendRpPbPred->SetFillColor(0);
	legendRpPbPred->SetLineColor(0);
	//egendRpPbCombine->SetNColumns(3);
	legendRpPbPred->SetTextSize(0.03);
	//legendRpPbCombine->SetEntrySeparation(0.02);
	legendRpPbPred->AddEntry(graphPi0ESP09sPi0KKP,"EPS09s KKP NLO");
	legendRpPbPred->AddEntry(graphPi0ESP09sPi0AKK,"EPS09s AKK NLO");
	legendRpPbPred->AddEntry(graphPi0DSS5000,"EPS09s fDSS NLO");
	legendRpPbPred->AddEntry(graphAsymmErrorsPi0DSS5000,"EPS09s fDSS errors","pef");
	
	
	graphPi0ESP09sPi0KKP->SetFillColor(0);
	graphPi0ESP09sPi0AKK->SetFillColor(0);
	graphPi0DSS5000->SetFillColor(0);
	
	
	
	
	histo2DRpPbESP09sCompared->DrawCopy();
	
	textRpPbCompared->Draw("sames");
	textRpPbEPS09s->Draw("sames");
	textRpPbEPS09sPre->Draw("sames");
	
	legendRpPbCompared->Draw("sames");
	legendRpPbPred->Draw("sames");
	
	
	graphPi0DSS5000->Draw("same,p,l");
	graphPi0ESP09sPi0AKK->Draw("same,p,l");
	graphPi0ESP09sPi0KKP->Draw("same,p,l");
	graphAsymmErrorsPi0DSS5000->Draw("same,E3");
	//	graphDalitzRpPb->Draw("p,same,e");
	graphRpPbPCM->Draw("p,same,e");
	LabelpPb3->Draw("same");	
	//graphAsymmChargedParticlesRpPb->Draw("p,same,e");
	
	DrawGammaLines(0., 15.,1., 1.,2.,kBlack);
	
	

	canvasRpPbEPS09s->Update(); 
	
	canvasRpPbEPS09s->SaveAs(Form("%s/RpPb_Comparison_PCM_Models.%s",outputDir.Data(),suffix.Data()));



	Int_t nlinesPi0CGC = 0;	
	inCGC.open(fileNameEPS09sPi0CGC,ios_base::in);
	
	Double_t xESPsPi0CGC[100], yESPsPi0CGC[100];
	
		
	while(!inCGC.eof()){
			nlinesPi0CGC++;
			inCGC >> xESPsPi0CGC[nlinesPi0CGC]  >> yESPsPi0CGC[nlinesPi0CGC];
			cout << nlinesPi0CGC << "         "  << xESPsPi0CGC[nlinesPi0CGC] << "         "  <<yESPsPi0CGC[nlinesPi0CGC]<<endl;
	
	}
	inCGC.close();
	

	TGraph* graphPi0CGC = new TGraph(nlinesPi0CGC,xESPsPi0CGC,yESPsPi0CGC);	
	
	TCanvas* canvasCGC = new TCanvas("canvasCGC","Combined predictions for R^{#pi0}_{pPb}",200,10,1200,700);
	
	
	DrawGammaCanvasSettings( canvasCGC,  0.1, 0.01, 0.015, 0.13);

	
	
	
	TH2F * histo2DRpPbESP09sCompared2 = new TH2F("histo2DRpPbESP09sCompared2","histo2DRpPbESP09sCompared2",1500,0.,15.,1000,0.3,1.5);
	SetStyleHistoTH2ForGraphs(histo2DRpPbESP09sCompared2, "p_{T} (GeV/c)","R^{#pi^{0}}_{pPb}", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 
	histo2DRpPbESP09sCompared2->DrawCopy(); 
	
	
	DrawGammaLines(0., 15.,1., 1.,2.,kBlack);
		
	
        TLatex *textRpPbCompared2 = new TLatex(0.25,0.40,"p-Pb min. bias #sqrt{#it{s}_{_{NN}}} = 5.02 TeV");
	SetStyleTLatex( textRpPbCompared2,0.035,4); 
	
	TLatex *textRpPbEPS09s2 = new TLatex(0.65,0.55,"p-Pb min. bias #sqrt{#it{s}_{_{NN}}} = 5.0 TeV");
	SetStyleTLatex( textRpPbEPS09s2,0.035,4); 
	
	TLatex *textRpPbEPS09sPre2 = new TLatex(0.65,0.5,"pQCD calculations");
	SetStyleTLatex( textRpPbEPS09sPre2,0.035,4); 
	
	TLatex *textCGC = new TLatex(0.65,0.40,"p-Pb min. bias #sqrt{#it{s}_{_{NN}}} = 5.02 TeV");
	SetStyleTLatex( textCGC,0.035,4); 
	
	TLatex *textCGC2 = new TLatex(0.65,0.23,"CGC calculations");
	SetStyleTLatex( textCGC2,0.035,4); 	
	//textRpPbCompared->Draw();
	
	
	//graphDalitzRpPb->SetFillColor(0);
	//graphRpPbPCM->SetFillColor(0);
	
	TLegend* legendRpPbCompared2 = new TLegend(0.23,0.28,0.45,0.35);
	legendRpPbCompared2->SetFillColor(0);
	legendRpPbCompared2->SetLineColor(0);
	//egendRpPbCombine->SetNColumns(3);
	legendRpPbCompared2->SetTextSize(0.03);
	//legendRpPbCombine->SetEntrySeparation(0.02);
	//	legendRpPbCompared->AddEntry(graphDalitzRpPb,"#pi^{0} #rightarrow e^{+}e^{-}#gamma #rightarrow e^{+}e^{-}e^{+}e^{-} (Dalitz)","pef");

	legendRpPbCompared2->AddEntry(graphRpPbPCM,"#pi^{0} #rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-} (PCM)","pef");
	TLatex *text2 = new TLatex(0.285,0.26,"including stat. and sys. errors");
	SetStyleTLatex( text2,0.025,3); 	

	TLatex *LabelpPb4 = new TLatex(0.25,0.45,"ALICE work in progress");
	SetStyleTLatex( LabelpPb4, 0.035,4);
	//	LabelpPb3->Draw("same");	
	
	
	TLegend* legendRpPbPred2 = new TLegend(0.65,0.3,0.85,0.48);
	legendRpPbPred2->SetFillColor(0);
	legendRpPbPred2->SetLineColor(0);
	//egendRpPbCombine->SetNColumns(3);
	legendRpPbPred2->SetTextSize(0.03);
	//legendRpPbCombine->SetEntrySeparation(0.02);
	legendRpPbPred2->AddEntry(graphPi0ESP09sPi0KKP,"EPS09s KKP NLO");
	legendRpPbPred2->AddEntry(graphPi0ESP09sPi0AKK,"EPS09s AKK NLO");
	legendRpPbPred2->AddEntry(graphPi0DSS5000,"EPS09s fDSS NLO");
	legendRpPbPred2->AddEntry(graphAsymmErrorsPi0DSS5000,"EPS09s fDSS errors","pef");
	//	legendRpPbPred2->AddEntry(graphPi0CGC,"Color Glas Condensate","p");	
	
	TLegend* legendRpPbPred3 = new TLegend(0.65,0.15,0.85,0.23);
	legendRpPbPred3->SetFillColor(0);
	legendRpPbPred3->SetLineColor(0);
	//egendRpPbCombine->SetNColumns(3);
	legendRpPbPred3->SetTextSize(0.03);
	//legendRpPbCombine->SetEntrySeparation(0.02);
	legendRpPbPred3->AddEntry(graphPi0CGC,"Color Glas Condensate","p");

	graphPi0ESP09sPi0KKP->SetFillColor(0);
	graphPi0ESP09sPi0AKK->SetFillColor(0);
	graphPi0DSS5000->SetFillColor(0);

	graphPi0CGC->SetFillColor(0);
	graphPi0CGC->SetMarkerColor(kGreen+1);
	graphPi0CGC->SetMarkerStyle(27);
	graphPi0CGC->SetMarkerSize(1);
	graphPi0CGC->SetLineColor(kGreen+3);
	graphPi0CGC->SetLineWidth(3);
	graphPi0CGC->SetLineStyle(1);
	//	graphPi0CGC->SetLineStyle(10);	
	

	histo2DRpPbESP09sCompared2->DrawCopy();
	legendRpPbPred3->Draw("sames");		
	textRpPbCompared2->Draw("sames");
	//	textRpPbEPS09s2->Draw("sames");
	textRpPbEPS09sPre2->Draw("sames");
 	textCGC2->Draw("sames");
	text2->Draw("sames");
	
	legendRpPbCompared2->Draw("sames");
	legendRpPbPred2->Draw("sames");
	//	legendRpPbPred3->Draw("sames");
	
	graphPi0CGC->Draw("same,p");
	graphPi0DSS5000->Draw("same,p,l");
	graphPi0ESP09sPi0AKK->Draw("same,p,l");
	graphPi0ESP09sPi0KKP->Draw("same,p,l");
	graphAsymmErrorsPi0DSS5000->Draw("same,E3");
	//	graphDalitzRpPb->Draw("p,same,e");
	graphRpPbPCM->Draw("p,same,e");
	LabelpPb4->Draw("same");	
	//graphAsymmChargedParticlesRpPb->Draw("p,same,e");
	
	DrawGammaLines(0., 15.,1., 1.,2.,kBlack);
	
	

	canvasCGC->Update(); 
	
	canvasCGC->SaveAs(Form("%s/RpPb_Comparison_PCM_CGC.%s",outputDir.Data(),suffix.Data()));

	
	



}

TGraphErrors *GetInterpolSpectrum2D(TGraphErrors *g1, TGraphErrors *g2,Double_t d1, Double_t d2,Double_t dSqrts= 10000)
{
         if(!g1) return 0x0;
         if(!g2) return 0x0;

         TGraphErrors *gInterpol = new TGraphErrors(g1->GetN());

         for(Int_t i = 0; i < g1->GetN(); i++)
         {

                 TGraphErrors *grint = new TGraphErrors(1);
                 grint->SetPoint(0, dSqrts, 0);
                 TGraphErrors *gToFit = new TGraphErrors(2);
                 gToFit->SetPoint(0, d1, g1->GetY()[i]);
                 gToFit->SetPointError(0, 0, g1->GetEY()[i]);

                 gToFit->SetPoint(1, d2, g2->GetY()[i]);
                 gToFit->SetPointError(1, 0, g2->GetEY()[i]);



                 TF1 *fPowerlaw = new TF1("fPowerlaw","[0]*x^([1])", 0,10000);
                 fPowerlaw->SetParameters(0, 0.1);
                 fPowerlaw->SetParameters(1, 2.0);

		 //                 for(Int_t l = 0; l < 10; l++) gToFit->Fit(fPowerlaw,"0Q");
                 for(Int_t l = 0; l < 10; l++) gToFit->Fit(fPowerlaw,"Q");

		 cout<<" in interpol"<< i<<" " << d1<< " " << d2<< " "<<dSqrts<< " "<<
		   g1->GetY()[i]<<" "<<g2->GetY()[i]<<" "<< fPowerlaw->Eval(dSqrts) << endl;                 



                    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint, 0.68);
                    gInterpol->SetPoint(i, g1->GetX()[i],fPowerlaw->Eval(dSqrts));
		       gInterpol->SetPointError(i, 0, grint->GetEY()[0]);
		    //    gInterpol->SetPointError(i, g1->GetEX()[i], grint->GetEY()[0]);

                 delete grint;
                 delete fPowerlaw;
         }

         return gInterpol;
}

TGraphErrors *GetInterpolSpectrum3D(TGraphErrors *g1, TGraphErrors *g2,TGraphErrors *g3, Double_t d1, Double_t d2, Double_t d3, Double_t dSqrts= 10000)
{
         if(!g1) return 0x0;
         if(!g2) return 0x0;
         if(!g3) return 0x0;

         TGraphErrors *gInterpol = new TGraphErrors(g1->GetN());

         for(Int_t i = 0; i < g1->GetN(); i++)
         {
                 TGraphErrors *grint = new TGraphErrors(1);
                 grint->SetPoint(0, dSqrts, 0);
                 TGraphErrors *gToFit = new TGraphErrors(3);
                 gToFit->SetPoint(0, d1, g1->GetY()[i]);
                 gToFit->SetPointError(0, 0, g1->GetEY()[i]);
                 gToFit->SetPoint(1, d2, g2->GetY()[i]);
                 gToFit->SetPointError(1, 0, g2->GetEY()[i]);
                 gToFit->SetPoint(2, d3, g3->GetY()[i]);
                 gToFit->SetPointError(2, 0, g3->GetEY()[i]);

                 TF1 *fPowerlaw = new TF1("fPowerlaw","[0]*x^([1])", 0,10000);
                 fPowerlaw->SetParameters(0, 0.1);
                 fPowerlaw->SetParameters(1, 2.0);

                 for(Int_t l = 0; l < 10; l++) gToFit->Fit(fPowerlaw,"0Q");

                    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint, 0.68);
                    gInterpol->SetPoint(i, g1->GetX()[i],
                    fPowerlaw->Eval(dSqrts));
                    gInterpol->SetPointError(i, 0, grint->GetEY()[0]);

                 delete grint;
                 delete fPowerlaw;
         }

         return gInterpol;
}


TGraphAsymmErrors* GetChargeParticlesRpPb(){

  // Plot: p8424_d3x1y1
  double p8424_d3x1y1_xval[] = { 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 
    0.975, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 
    1.95, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 
    3.9, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 
    10.5, 11.5, 12.5, 13.5, 15.0, 18.0 };
  double p8424_d3x1y1_xerrminus[] = { 0.025000000000000022, 0.02499999999999991, 0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.02499999999999991, 0.025000000000000022, 0.025000000000000022, 
    0.025000000000000022, 0.050000000000000044, 0.04999999999999982, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.04999999999999982, 0.050000000000000044, 0.050000000000000044, 
    0.050000000000000044, 0.10000000000000009, 0.09999999999999964, 0.10000000000000009, 0.10000000000000009, 0.10000000000000009, 0.10000000000000009, 0.09999999999999964, 0.10000000000000009, 0.10000000000000009, 
    0.10000000000000009, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 
    0.5, 0.5, 0.5, 0.5, 1.0, 2.0 };
  double p8424_d3x1y1_xerrplus[] = { 0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.02499999999999991, 0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.02499999999999991, 
    0.025000000000000022, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.04999999999999982, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.04999999999999982, 
    0.050000000000000044, 0.10000000000000009, 0.10000000000000009, 0.10000000000000009, 0.09999999999999964, 0.10000000000000009, 0.10000000000000009, 0.10000000000000009, 0.10000000000000009, 0.09999999999999964, 
    0.10000000000000009, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 
    0.5, 0.5, 0.5, 0.5, 1.0, 2.0 };

    double p8424_d3x1y1_yval[] = { 0.5824, 0.5987, 0.6139, 0.6337, 0.6514, 0.6716, 0.6878, 0.6989, 0.7118, 
    0.722, 0.7485, 0.7822, 0.8046, 0.8306, 0.8539, 0.8811, 0.8947, 0.9145, 0.9431, 
    0.9584, 0.9844, 1.01, 1.034, 1.051, 1.051, 1.08, 1.083, 1.104, 1.102, 
    1.111, 1.1, 1.091, 1.111, 1.058, 1.035, 1.087, 1.064, 1.064, 1.11, 
    0.9183, 1.144, 1.061, 1.17, 0.9558, 1.121 };
  double p8424_d3x1y1_yerrminus[] = { 0.05600723167591842, 0.05670881765651617, 0.056910631695668255, 0.058712264476853564, 0.060314011639087645, 0.06221575363201831, 0.06371765846294103, 0.06472232999514155, 0.0659245781177248, 
    0.06692697512961422, 0.06931623186527092, 0.07251992829560713, 0.07452684080249208, 0.07693438762997987, 0.07914271918502674, 0.08165512843661445, 0.08296969326205805, 0.08478519918004557, 0.08740583504549339, 
    0.08892429364352579, 0.09128335006998813, 0.09413288479590966, 0.09618731725128839, 0.09725224933131367, 0.09732933781753578, 0.10031948963187563, 0.1004987562112089, 0.10259142264341595, 0.10270345661174213, 
    0.1039471019317037, 0.10270345661174213, 0.10196568050084304, 0.10932977636490435, 0.10413932974625868, 0.10235721762533408, 0.1078563859954523, 0.10412012293500235, 0.1077032961426901, 0.11764777940955791, 
    0.10736125930707034, 0.13917614738165443, 0.14509651959988565, 0.1776091213873882, 0.14606854555310667, 0.17462244987400674 };
  double p8424_d3x1y1_yerrplus[] = { 0.05600723167591842, 0.05670881765651617, 0.056910631695668255, 0.058712264476853564, 0.060314011639087645, 0.06221575363201831, 0.06371765846294103, 0.06472232999514155, 0.0659245781177248, 
    0.06692697512961422, 0.06931623186527092, 0.07251992829560713, 0.07452684080249208, 0.07693438762997987, 0.07914271918502674, 0.08165512843661445, 0.08296969326205805, 0.08478519918004557, 0.08740583504549339, 
    0.08892429364352579, 0.09128335006998813, 0.09413288479590966, 0.09618731725128839, 0.09725224933131367, 0.09732933781753578, 0.10031948963187563, 0.1004987562112089, 0.10259142264341595, 0.10270345661174213, 
    0.1039471019317037, 0.10270345661174213, 0.10196568050084304, 0.10932977636490435, 0.10413932974625868, 0.10235721762533408, 0.1078563859954523, 0.10412012293500235, 0.1077032961426901, 0.11764777940955791, 
    0.10736125930707034, 0.13917614738165443, 0.14509651959988565, 0.1776091213873882, 0.14606854555310667, 0.17462244987400674 };
  double p8424_d3x1y1_ystatminus[] = { 9.0E-4, 0.001, 0.0011, 0.0012, 0.0013, 0.0014, 0.0015, 0.0017, 0.0018, 
    0.0019, 0.0015, 0.0017, 0.002, 0.0023, 0.0026, 0.003, 0.0034, 0.0038, 0.0043, 
    0.0047, 0.0039, 0.005, 0.006, 0.007, 0.008, 0.008, 0.01, 0.011, 0.012, 
    0.014, 0.012, 0.014, 0.017, 0.021, 0.026, 0.032, 0.029, 0.04, 0.055, 
    0.064, 0.089, 0.107, 0.141, 0.1159, 0.138 };
  double p8424_d3x1y1_ystatplus[] = { 9.0E-4, 0.001, 0.0011, 0.0012, 0.0013, 0.0014, 0.0015, 0.0017, 0.0018, 
    0.0019, 0.0015, 0.0017, 0.002, 0.0023, 0.0026, 0.003, 0.0034, 0.0038, 0.0043, 
    0.0047, 0.0039, 0.005, 0.006, 0.007, 0.008, 0.008, 0.01, 0.011, 0.012, 
    0.014, 0.012, 0.014, 0.017, 0.021, 0.026, 0.032, 0.029, 0.04, 0.055, 
    0.064, 0.089, 0.107, 0.141, 0.1159, 0.138 };
  int p8424_d3x1y1_numpoints = 45;
  TGraphAsymmErrors* p8424_d3x1y1 = new TGraphAsymmErrors(p8424_d3x1y1_numpoints, p8424_d3x1y1_xval, p8424_d3x1y1_yval, p8424_d3x1y1_xerrminus, p8424_d3x1y1_xerrplus, p8424_d3x1y1_yerrminus, p8424_d3x1y1_yerrplus);
  p8424_d3x1y1->SetName("/HepData/8424/d3x1y1");
  p8424_d3x1y1->SetTitle("/HepData/8424/d3x1y1");
  
  return p8424_d3x1y1;
  
  
}


	
TGraphErrors* RemovePointsFromGraph(TGraphErrors *graph, Int_t NbPoints){

   TGraphErrors* dummyGraph = (TGraphErrors*)	graph->Clone();
	Double_t * xValue = dummyGraph->GetX();
	Double_t * yValue = dummyGraph->GetY();
	Double_t* xError = dummyGraph->GetEX();
	Double_t* yError = dummyGraph->GetEY();

	Int_t nPoints = dummyGraph->GetN();
	Int_t nPoints_new = nPoints  - NbPoints;
	for (Int_t i = NbPoints; i < nPoints; i++){
		yValue[i-NbPoints] = yValue[i];
		xValue[i-NbPoints] = xValue[i];
		xError[i-NbPoints] = xError[i];
		yError[i-NbPoints] = yError[i];
	}
	TGraphErrors* returnGraph = new TGraphErrors(nPoints_new,xValue,yValue, xError,yError);
	return returnGraph;
}


