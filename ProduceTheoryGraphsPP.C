/*****************************************************************************
******         provided by Gamma Conversion Group, PWG4,                ******
******        Ana Marin, marin@physi.uni-heidelberg.de                  ******
******           Kathrin Koch, kkoch@physi.uni-heidelberg.de            ******
******        Friederike Bock, friederike.bock@cern.ch                  ******
******        Lucia Leardini, lucia.leardini@cern.ch                    ******
*****************************************************************************/

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

extern TRandom*    gRandom;
extern TBenchmark*    gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*      gMinuit;

Double_t    xSection900GeV =     47.78*1e-3;
Double_t    xSection2760GeV =    55.416*1e-3;
Double_t    xSection7000GeV =    62.22*1e-3;
Double_t    xSection8000GeV =    55.74*1e-3;
Double_t    recalcBarn =      1e12; //NLO in pbarn!!!!


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

TGraphAsymmErrors* ScaleGraphAsym (TGraphAsymmErrors* graph, Double_t scaleFac){
	TGraphAsymmErrors* dummyGraph = (TGraphAsymmErrors*)graph->Clone(Form("%s_Scaled",graph->GetName()));

	Double_t * xValue = dummyGraph->GetX(); 
	Double_t * yValue = dummyGraph->GetY();
	Double_t* xErrorLow = dummyGraph->GetEXlow();
	Double_t* xErrorHigh = dummyGraph->GetEXhigh();
	Double_t* yErrorLow = dummyGraph->GetEYlow();
	Double_t* yErrorHigh = dummyGraph->GetEYhigh();
	Int_t nPoints = dummyGraph->GetN();
	for (Int_t i = 0; i < nPoints; i++){
		yValue[i] = yValue[i]*scaleFac;
		yErrorLow[i] = yErrorLow[i]*scaleFac;
		yErrorHigh[i] = yErrorHigh[i]*scaleFac;
	}
	TGraphAsymmErrors* returnGraph =  new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh); 
	return returnGraph;
}

void ProduceTheoryGraphsPP(){    
    
	StyleSettingsThesis();    
	SetPlotStyle();


	Double_t       ptNLOEta900GeV[100];
	Double_t       muHalfEta900GeV[100];
	Double_t       muOneEta900GeV[100];
	Double_t       muTwoEta900GeV[100];
	Int_t       nlinesNLOEta900GeV =       0;
	
	TString fileNameNLOEta900GeV = "ExternalInput/Theory/ALICENLOcalcEtaVogelsang900GeV.dat";
	ifstream  fileNLOEta900GeV;
	fileNLOEta900GeV.open(fileNameNLOEta900GeV,ios_base::in);
	cout << fileNameNLOEta900GeV << endl;
	
	while(!fileNLOEta900GeV.eof()){
		nlinesNLOEta900GeV++;
		fileNLOEta900GeV >> ptNLOEta900GeV[nlinesNLOEta900GeV] >> muHalfEta900GeV[nlinesNLOEta900GeV] >> muOneEta900GeV[nlinesNLOEta900GeV] >> muTwoEta900GeV[nlinesNLOEta900GeV]; 
		cout << nlinesNLOEta900GeV << "         "  << ptNLOEta900GeV[nlinesNLOEta900GeV] << "         "  << muHalfEta900GeV[nlinesNLOEta900GeV] << "         "  << muOneEta900GeV[nlinesNLOEta900GeV] << "         "  << muTwoEta900GeV[nlinesNLOEta900GeV] << endl;;
	}
	fileNLOEta900GeV.close();
	TGraph* graphNLOCalcInvSecEtaMuHalf900GeV = new TGraph(nlinesNLOEta900GeV-1,ptNLOEta900GeV,muHalfEta900GeV); 
	TGraph* graphNLOCalcInvSecEtaMuOne900GeV = new TGraph(nlinesNLOEta900GeV-1,ptNLOEta900GeV,muOneEta900GeV); 
	TGraph* graphNLOCalcInvSecEtaMuTwo900GeV = new TGraph(nlinesNLOEta900GeV-1,ptNLOEta900GeV,muTwoEta900GeV); 
	graphNLOCalcInvSecEtaMuHalf900GeV->RemovePoint(0);
	graphNLOCalcInvSecEtaMuHalf900GeV->RemovePoint(0);
	graphNLOCalcInvSecEtaMuOne900GeV->RemovePoint(0);
	graphNLOCalcInvSecEtaMuTwo900GeV->RemovePoint(0);
	
	TGraph* graphNLOCalcInvYieldEtaMuHalf900GeV = ScaleGraph(graphNLOCalcInvSecEtaMuHalf900GeV, 1/(xSection900GeV*recalcBarn));
	TGraph* graphNLOCalcInvYieldEtaMuOne900GeV =  ScaleGraph(graphNLOCalcInvSecEtaMuOne900GeV, 1/(xSection900GeV*recalcBarn));
	TGraph* graphNLOCalcInvYieldEtaMuTwo900GeV =  ScaleGraph(graphNLOCalcInvSecEtaMuTwo900GeV, 1/(xSection900GeV*recalcBarn));
	
	
	Double_t       ptNLOPi0900GeV[100];
	Double_t       muHalfPi0900GeV[100];
	Double_t       muOnePi0900GeV[100];
	Double_t       muTwoPi0900GeV[100];
	Int_t       nlinesNLOPi0900GeV =       0;
	
	TString fileNameNLOPi0900GeV = "ExternalInput/Theory/ALICENLOcalcPi0Vogelsang900Gev.dat";
	ifstream  fileNLOPi0900GeV;
	fileNLOPi0900GeV.open(fileNameNLOPi0900GeV,ios_base::in);
	cout << fileNameNLOPi0900GeV << endl;
	
	while(!fileNLOPi0900GeV.eof()){
		nlinesNLOPi0900GeV++;
		fileNLOPi0900GeV >> ptNLOPi0900GeV[nlinesNLOPi0900GeV] >> muHalfPi0900GeV[nlinesNLOPi0900GeV] >> muOnePi0900GeV[nlinesNLOPi0900GeV] >> muTwoPi0900GeV[nlinesNLOPi0900GeV]; 
		cout << nlinesNLOPi0900GeV << "         "  << ptNLOPi0900GeV[nlinesNLOPi0900GeV] << "         "  << muHalfPi0900GeV[nlinesNLOPi0900GeV] << "         "  << muOnePi0900GeV[nlinesNLOPi0900GeV] << "         "  << muTwoPi0900GeV[nlinesNLOPi0900GeV] << endl;;
	}
	fileNLOPi0900GeV.close();
	TGraph* graphNLOCalcInvSecPi0MuHalf900GeV = new TGraph(nlinesNLOPi0900GeV,ptNLOPi0900GeV,muHalfPi0900GeV); 
	graphNLOCalcInvSecPi0MuHalf900GeV->RemovePoint(0);
	graphNLOCalcInvSecPi0MuHalf900GeV->RemovePoint(0);
	TGraph* graphNLOCalcInvSecPi0MuOne900GeV = new TGraph(nlinesNLOPi0900GeV,ptNLOPi0900GeV,muOnePi0900GeV); 
	graphNLOCalcInvSecPi0MuOne900GeV->RemovePoint(0);
	TGraph* graphNLOCalcInvSecPi0MuTwo900GeV = new TGraph(nlinesNLOPi0900GeV,ptNLOPi0900GeV,muTwoPi0900GeV); 
	graphNLOCalcInvSecPi0MuTwo900GeV->RemovePoint(0);

	TGraph* graphNLOCalcInvYieldPi0MuHalf900GeV = ScaleGraph(graphNLOCalcInvSecPi0MuHalf900GeV, 1/(xSection900GeV*recalcBarn));
	TGraph* graphNLOCalcInvYieldPi0MuOne900GeV =  ScaleGraph(graphNLOCalcInvSecPi0MuOne900GeV, 1/(xSection900GeV*recalcBarn));
	TGraph* graphNLOCalcInvYieldPi0MuTwo900GeV =  ScaleGraph(graphNLOCalcInvSecPi0MuTwo900GeV, 1/(xSection900GeV*recalcBarn));

	Double_t       ptNLOPi0900GeVBKK[100];
	Double_t       muTwoPi0900GeVBKK[100];
	Int_t       nlinesNLOPi0900GeVBKK =       0;
	
	TString fileNameNLOPi0900GeVBKK = "ExternalInput/Theory/lhc_900_CTEQ5M_BKK_20.dat";
	ifstream  fileNLOPi0900GeVBKK;
	fileNLOPi0900GeVBKK.open(fileNameNLOPi0900GeVBKK,ios_base::in);
	cout << fileNameNLOPi0900GeVBKK << endl;
	
	while(!fileNLOPi0900GeVBKK.eof()){
		nlinesNLOPi0900GeVBKK++;
		fileNLOPi0900GeVBKK >> ptNLOPi0900GeVBKK[nlinesNLOPi0900GeVBKK] >> muTwoPi0900GeVBKK[nlinesNLOPi0900GeVBKK];  
	}
	fileNLOPi0900GeVBKK.close();
	TGraph* graphNLOCalcBKKInvSecPi0MuTwo900GeV = new TGraph(nlinesNLOPi0900GeVBKK,ptNLOPi0900GeVBKK,muTwoPi0900GeVBKK); 
	graphNLOCalcBKKInvSecPi0MuTwo900GeV->RemovePoint(0);
	TGraph* graphNLOCalcBKKInvYieldPi0MuTwo900GeV =  ScaleGraph(graphNLOCalcBKKInvSecPi0MuTwo900GeV, 1/(xSection900GeV*recalcBarn));
	
	
	Double_t       ptNLOEta2760GeV[100];
	Double_t       muHalfEta2760GeV[100];
	Double_t       muOneEta2760GeV[100];
	Double_t       muTwoEta2760GeV[100];
	Int_t       nlinesNLOEta2760GeV =       0;
	
	TString fileNameNLOEta2760GeV = "ExternalInput/Theory/ALICENLOcalcEtaVogelsang2760GeV.dat";
	ifstream  fileNLOEta2760GeV;
	fileNLOEta2760GeV.open(fileNameNLOEta2760GeV,ios_base::in);
	cout << fileNameNLOEta2760GeV << endl;
	
	while(!fileNLOEta2760GeV.eof()){
		nlinesNLOEta2760GeV++;
		fileNLOEta2760GeV >> ptNLOEta2760GeV[nlinesNLOEta2760GeV] >> muHalfEta2760GeV[nlinesNLOEta2760GeV] >> muOneEta2760GeV[nlinesNLOEta2760GeV] >> muTwoEta2760GeV[nlinesNLOEta2760GeV]; 
		cout << nlinesNLOEta2760GeV << "         "  << ptNLOEta2760GeV[nlinesNLOEta2760GeV] << "         "  << muHalfEta2760GeV[nlinesNLOEta2760GeV] << "         "  << muOneEta2760GeV[nlinesNLOEta2760GeV] << "         "  << muTwoEta2760GeV[nlinesNLOEta2760GeV] << endl;;
	}
	fileNLOEta2760GeV.close();
	TGraph* graphNLOCalcInvSecEtaMuHalf2760GeV = new TGraph(nlinesNLOEta2760GeV,ptNLOEta2760GeV,muHalfEta2760GeV); 
	graphNLOCalcInvSecEtaMuHalf2760GeV->RemovePoint(0);
	graphNLOCalcInvSecEtaMuHalf2760GeV->RemovePoint(0);
	TGraph* graphNLOCalcInvSecEtaMuOne2760GeV = new TGraph(nlinesNLOEta2760GeV,ptNLOEta2760GeV,muOneEta2760GeV); 
	graphNLOCalcInvSecEtaMuOne2760GeV->RemovePoint(0);
	TGraph* graphNLOCalcInvSecEtaMuTwo2760GeV = new TGraph(nlinesNLOEta2760GeV,ptNLOEta2760GeV,muTwoEta2760GeV); 
	graphNLOCalcInvSecEtaMuTwo2760GeV->RemovePoint(0);

	TGraph* graphNLOCalcInvYieldEtaMuHalf2760GeV = ScaleGraph(graphNLOCalcInvSecEtaMuHalf2760GeV, 1/(xSection2760GeV*recalcBarn));
	TGraph* graphNLOCalcInvYieldEtaMuOne2760GeV =  ScaleGraph(graphNLOCalcInvSecEtaMuOne2760GeV, 1/(xSection2760GeV*recalcBarn));
	TGraph* graphNLOCalcInvYieldEtaMuTwo2760GeV =  ScaleGraph(graphNLOCalcInvSecEtaMuTwo2760GeV, 1/(xSection2760GeV*recalcBarn));
	
	Double_t       ptNLOPi02760GeV[100];
	Double_t       muHalfPi02760GeV[100];
	Double_t       muOnePi02760GeV[100];
	Double_t       muTwoPi02760GeV[100];
	Int_t       nlinesNLOPi02760GeV =       0;
	
	TString fileNameNLOPi02760GeV = "ExternalInput/Theory/ALICENLOcalcPi0Vogelsang2760Gev.dat";
	ifstream  fileNLOPi02760GeV;
	fileNLOPi02760GeV.open(fileNameNLOPi02760GeV,ios_base::in);
	cout << fileNameNLOPi02760GeV << endl;
	
	while(!fileNLOPi02760GeV.eof()){
		nlinesNLOPi02760GeV++;
		fileNLOPi02760GeV >> ptNLOPi02760GeV[nlinesNLOPi02760GeV] >> muHalfPi02760GeV[nlinesNLOPi02760GeV] >> muOnePi02760GeV[nlinesNLOPi02760GeV] >> muTwoPi02760GeV[nlinesNLOPi02760GeV]; 
		cout << nlinesNLOPi02760GeV << "         "  << ptNLOPi02760GeV[nlinesNLOPi02760GeV] << "         "  << muHalfPi02760GeV[nlinesNLOPi02760GeV] << "         "  << muOnePi02760GeV[nlinesNLOPi02760GeV] << "         "  << muTwoPi02760GeV[nlinesNLOPi02760GeV] << endl;;
	}
	fileNLOPi02760GeV.close();
	TGraph* graphNLOCalcInvSecPi0MuHalf2760GeV = new TGraph(nlinesNLOPi02760GeV,ptNLOPi02760GeV,muHalfPi02760GeV); 
	graphNLOCalcInvSecPi0MuHalf2760GeV->RemovePoint(0);
	graphNLOCalcInvSecPi0MuHalf2760GeV->RemovePoint(0);
	TGraph* graphNLOCalcInvSecPi0MuOne2760GeV = new TGraph(nlinesNLOPi02760GeV,ptNLOPi02760GeV,muOnePi02760GeV); 
	graphNLOCalcInvSecPi0MuOne2760GeV->RemovePoint(0);
	TGraph* graphNLOCalcInvSecPi0MuTwo2760GeV = new TGraph(nlinesNLOPi02760GeV,ptNLOPi02760GeV,muTwoPi02760GeV); 
	graphNLOCalcInvSecPi0MuTwo2760GeV->RemovePoint(0);

	TGraph* graphNLOCalcInvYieldPi0MuHalf2760GeV = ScaleGraph(graphNLOCalcInvSecPi0MuHalf2760GeV, 1/(xSection2760GeV*recalcBarn));
	TGraph* graphNLOCalcInvYieldPi0MuOne2760GeV =  ScaleGraph(graphNLOCalcInvSecPi0MuOne2760GeV, 1/(xSection2760GeV*recalcBarn));
	TGraph* graphNLOCalcInvYieldPi0MuTwo2760GeV =  ScaleGraph(graphNLOCalcInvSecPi0MuTwo2760GeV, 1/(xSection2760GeV*recalcBarn));
	
	Double_t       ptNLOEta7000GeV[100];
	Double_t       muHalfEta7000GeV[100];
	Double_t       muOneEta7000GeV[100];
	Double_t       muTwoEta7000GeV[100];
	Int_t       nlinesNLOEta7000GeV =       0;
	
	TString fileNameNLOEta7000GeV = "ExternalInput/Theory/ALICENLOcalcEtaVogelsang7TeV.dat";
	ifstream  fileNLOEta7000GeV;
	fileNLOEta7000GeV.open(fileNameNLOEta7000GeV,ios_base::in);
	cout << fileNameNLOEta7000GeV << endl;
	
	while(!fileNLOEta7000GeV.eof()){
		nlinesNLOEta7000GeV++;
		fileNLOEta7000GeV >> ptNLOEta7000GeV[nlinesNLOEta7000GeV] >> muHalfEta7000GeV[nlinesNLOEta7000GeV] >> muOneEta7000GeV[nlinesNLOEta7000GeV] >> muTwoEta7000GeV[nlinesNLOEta7000GeV]; 
		cout << nlinesNLOEta7000GeV << "         "  << ptNLOEta7000GeV[nlinesNLOEta7000GeV] << "         "  << muHalfEta7000GeV[nlinesNLOEta7000GeV] << "         "  << muOneEta7000GeV[nlinesNLOEta7000GeV] << "         "  << muTwoEta7000GeV[nlinesNLOEta7000GeV] << endl;;
	}
	fileNLOEta7000GeV.close();
	TGraph* graphNLOCalcInvSecEtaMuHalf7000GeV = new TGraph(nlinesNLOEta7000GeV,ptNLOEta7000GeV,muHalfEta7000GeV); 
	graphNLOCalcInvSecEtaMuHalf7000GeV->RemovePoint(0);
	graphNLOCalcInvSecEtaMuHalf7000GeV->RemovePoint(0);
	TGraph* graphNLOCalcInvSecEtaMuOne7000GeV = new TGraph(nlinesNLOEta7000GeV,ptNLOEta7000GeV,muOneEta7000GeV); 
	graphNLOCalcInvSecEtaMuOne7000GeV->RemovePoint(0);
	TGraph* graphNLOCalcInvSecEtaMuTwo7000GeV = new TGraph(nlinesNLOEta7000GeV,ptNLOEta7000GeV,muTwoEta7000GeV); 
	graphNLOCalcInvSecEtaMuTwo7000GeV->RemovePoint(0);

	TGraph* graphNLOCalcInvYieldEtaMuHalf7000GeV = ScaleGraph(graphNLOCalcInvSecEtaMuHalf7000GeV, 1/(xSection7000GeV*recalcBarn));
	TGraph* graphNLOCalcInvYieldEtaMuOne7000GeV =  ScaleGraph(graphNLOCalcInvSecEtaMuOne7000GeV, 1/(xSection7000GeV*recalcBarn));
	TGraph* graphNLOCalcInvYieldEtaMuTwo7000GeV =  ScaleGraph(graphNLOCalcInvSecEtaMuTwo7000GeV, 1/(xSection7000GeV*recalcBarn));
	
	Double_t       ptNLOPi07000GeV[100];
	Double_t       muHalfPi07000GeV[100];
	Double_t       muOnePi07000GeV[100];
	Double_t       muTwoPi07000GeV[100];
	Int_t       nlinesNLOPi07000GeV =       0;
	
	TString fileNameNLOPi07000GeV = "ExternalInput/Theory/ALICENLOcalcPi0Vogelsang7Tev.dat";
	ifstream  fileNLOPi07000GeV;
	fileNLOPi07000GeV.open(fileNameNLOPi07000GeV,ios_base::in);
	cout << fileNameNLOPi07000GeV << endl;
	
	while(!fileNLOPi07000GeV.eof()){
		nlinesNLOPi07000GeV++;
		fileNLOPi07000GeV >> ptNLOPi07000GeV[nlinesNLOPi07000GeV] >> muHalfPi07000GeV[nlinesNLOPi07000GeV] >> muOnePi07000GeV[nlinesNLOPi07000GeV] >> muTwoPi07000GeV[nlinesNLOPi07000GeV]; 
		cout << nlinesNLOPi07000GeV << "         "  << ptNLOPi07000GeV[nlinesNLOPi07000GeV] << "         "  << muHalfPi07000GeV[nlinesNLOPi07000GeV] << "         "  << muOnePi07000GeV[nlinesNLOPi07000GeV] << "         "  << muTwoPi07000GeV[nlinesNLOPi07000GeV] << endl;;
	}
	fileNLOPi07000GeV.close();
	TGraph* graphNLOCalcInvSecPi0MuHalf7000GeV = new TGraph(nlinesNLOPi07000GeV,ptNLOPi07000GeV,muHalfPi07000GeV); 
	graphNLOCalcInvSecPi0MuHalf7000GeV->RemovePoint(0);
	graphNLOCalcInvSecPi0MuHalf7000GeV->RemovePoint(0);
	TGraph* graphNLOCalcInvSecPi0MuOne7000GeV = new TGraph(nlinesNLOPi07000GeV,ptNLOPi07000GeV,muOnePi07000GeV); 
	graphNLOCalcInvSecPi0MuOne7000GeV->RemovePoint(0);
	TGraph* graphNLOCalcInvSecPi0MuTwo7000GeV = new TGraph(nlinesNLOPi07000GeV,ptNLOPi07000GeV,muTwoPi07000GeV); 
	graphNLOCalcInvSecPi0MuTwo7000GeV->RemovePoint(0);
	
	TGraph* graphNLOCalcInvYieldPi0MuHalf7000GeV = ScaleGraph(graphNLOCalcInvSecPi0MuHalf7000GeV, 1/(xSection7000GeV*recalcBarn));
	TGraph* graphNLOCalcInvYieldPi0MuOne7000GeV =  ScaleGraph(graphNLOCalcInvSecPi0MuOne7000GeV, 1/(xSection7000GeV*recalcBarn));
	TGraph* graphNLOCalcInvYieldPi0MuTwo7000GeV =  ScaleGraph(graphNLOCalcInvSecPi0MuTwo7000GeV, 1/(xSection7000GeV*recalcBarn));

	Double_t       ptNLOPi07000GeVBKK[100];
	Double_t       muTwoPi07000GeVBKK[100];
	Int_t       nlinesNLOPi07000GeVBKK =       0;
	
	TString fileNameNLOPi07000GeVBKK = "ExternalInput/Theory/lhc_7000_CTEQ5M_BKK_20.dat";
	ifstream  fileNLOPi07000GeVBKK;
	fileNLOPi07000GeVBKK.open(fileNameNLOPi07000GeVBKK,ios_base::in);
	cout << fileNameNLOPi07000GeVBKK << endl;
	
	while(!fileNLOPi07000GeVBKK.eof()){
		nlinesNLOPi07000GeVBKK++;
		fileNLOPi07000GeVBKK >> ptNLOPi07000GeVBKK[nlinesNLOPi07000GeVBKK] >> muTwoPi07000GeVBKK[nlinesNLOPi07000GeVBKK];  
	}
	fileNLOPi07000GeVBKK.close();
	TGraph* graphNLOCalcBKKInvSecPi0MuTwo7000GeV = new TGraph(nlinesNLOPi07000GeVBKK,ptNLOPi07000GeVBKK,muTwoPi07000GeVBKK); 
	graphNLOCalcBKKInvSecPi0MuTwo7000GeV->RemovePoint(0);
	TGraph* graphNLOCalcBKKInvYieldPi0MuTwo7000GeV =  ScaleGraph(graphNLOCalcBKKInvSecPi0MuTwo7000GeV, 1/(xSection7000GeV*recalcBarn));

	Double_t       ptNLOPi07000GeVDSS[100];
	Double_t       muTwoPi07000GeVDSS[100];
	Double_t       energyPi07000GeVDSS[100];
	Int_t       nlinesNLOPi07000GeVDSS =       0;
		TString fileNameNLOPi07000GeVDSS = "ExternalInput/Theory/ALICE-DPT-DETA-PI0-7000.RES";
	ifstream  fileNLOPi07000GeVDSS;
	fileNLOPi07000GeVDSS.open(fileNameNLOPi07000GeVDSS,ios_base::in);
	
	while(!fileNLOPi07000GeVDSS.eof()){
			nlinesNLOPi07000GeVDSS++;
			fileNLOPi07000GeVDSS >> ptNLOPi07000GeVDSS[nlinesNLOPi07000GeVDSS]  >> energyPi07000GeVDSS[nlinesNLOPi07000GeVDSS] >> muTwoPi07000GeVDSS[nlinesNLOPi07000GeVDSS]; 
			
	
	}
	fileNLOPi07000GeVDSS.close();
	TGraph* graphNLOCalcDSSInvSecPi0MuTwo7000GeV = new TGraph(nlinesNLOPi07000GeVDSS,ptNLOPi07000GeVDSS,muTwoPi07000GeVDSS);   
	graphNLOCalcDSSInvSecPi0MuTwo7000GeV->RemovePoint(0);
	TGraph* graphNLOCalcDSSInvYieldPi0MuTwo7000GeV =  ScaleGraph(graphNLOCalcDSSInvSecPi0MuTwo7000GeV, 1/(xSection7000GeV*recalcBarn));

	Double_t       ptNLOPi02760GeVDSS[100];
	Double_t       muTwoPi02760GeVDSS[100];
	Double_t       energyPi02760GeVDSS[100];
	Int_t       nlinesNLOPi02760GeVDSS =       0;
		TString fileNameNLOPi02760GeVDSS = "ExternalInput/Theory/ALICE-DPT-DETA-PI0-2760.RES";
	ifstream  fileNLOPi02760GeVDSS;
	fileNLOPi02760GeVDSS.open(fileNameNLOPi02760GeVDSS,ios_base::in);
	
	while(!fileNLOPi02760GeVDSS.eof()){
			nlinesNLOPi02760GeVDSS++;
			fileNLOPi02760GeVDSS >> ptNLOPi02760GeVDSS[nlinesNLOPi02760GeVDSS]  >> energyPi02760GeVDSS[nlinesNLOPi02760GeVDSS] >> muTwoPi02760GeVDSS[nlinesNLOPi02760GeVDSS]; 
			
	
	}
	fileNLOPi02760GeVDSS.close();
	TGraph* graphNLOCalcDSSInvSecPi0MuTwo2760GeV = new TGraph(nlinesNLOPi02760GeVDSS,ptNLOPi02760GeVDSS,muTwoPi02760GeVDSS);   
	graphNLOCalcDSSInvSecPi0MuTwo2760GeV->RemovePoint(0);
	TGraph* graphNLOCalcDSSInvYieldPi0MuTwo2760GeV =  ScaleGraph(graphNLOCalcDSSInvSecPi0MuTwo2760GeV, 1/(xSection2760GeV*recalcBarn));

	Double_t       ptNLOPi0900GeVDSS[100];
	Double_t       muTwoPi0900GeVDSS[100];
	Double_t       energyPi0900GeVDSS[100];
	Int_t       nlinesNLOPi0900GeVDSS =       0;
		TString fileNameNLOPi0900GeVDSS = "ExternalInput/Theory/ALICE-DPT-DETA-PI0-900.RES";
	ifstream  fileNLOPi0900GeVDSS;
	fileNLOPi0900GeVDSS.open(fileNameNLOPi0900GeVDSS,ios_base::in);
	
	while(!fileNLOPi0900GeVDSS.eof()){
			nlinesNLOPi0900GeVDSS++;
			fileNLOPi0900GeVDSS >> ptNLOPi0900GeVDSS[nlinesNLOPi0900GeVDSS]  >> energyPi0900GeVDSS[nlinesNLOPi0900GeVDSS] >> muTwoPi0900GeVDSS[nlinesNLOPi0900GeVDSS]; 
			
	
	}
	fileNLOPi0900GeVDSS.close();
	TGraph* graphNLOCalcDSSInvSecPi0MuTwo900GeV = new TGraph(nlinesNLOPi0900GeVDSS,ptNLOPi0900GeVDSS,muTwoPi0900GeVDSS);   
	graphNLOCalcDSSInvSecPi0MuTwo900GeV->RemovePoint(0);
	TGraph* graphNLOCalcDSSInvYieldPi0MuTwo900GeV =  ScaleGraph(graphNLOCalcDSSInvSecPi0MuTwo900GeV, 1/(xSection900GeV*recalcBarn));
	
	
	Double_t* valueNLOMuHalfEta = graphNLOCalcInvSecEtaMuHalf7000GeV->GetY();
	Double_t* valueNLOMuOneEta =  graphNLOCalcInvSecEtaMuOne7000GeV->GetY();
	Double_t* valueNLOMuTwoEta =  graphNLOCalcInvSecEtaMuTwo7000GeV->GetY();
	Double_t* valueNLOMuHalfPi0 = graphNLOCalcInvSecPi0MuHalf7000GeV->GetY();
	Double_t* valueNLOMuOnePi0 =  graphNLOCalcInvSecPi0MuOne7000GeV->GetY();
	Double_t* valueNLOMuTwoPi0 =  graphNLOCalcInvSecPi0MuTwo7000GeV->GetY();
	Double_t* xValueNLO =      graphNLOCalcInvSecPi0MuOne7000GeV->GetX();
	Int_t    xNBins =          graphNLOCalcInvSecPi0MuOne7000GeV->GetN();
	Double_t    valueNLOEtaToPi0NLOMuHalf[50];
	Double_t    valueNLOEtaToPi0NLOMuOne[50];
	Double_t    valueNLOEtaToPi0NLOMuTwo[50];
	
	
	for ( Int_t n = 0; n < xNBins+1; n++){
		if (n == 0){
			valueNLOEtaToPi0NLOMuHalf[n] = 0.;
		} else { 
			if (valueNLOMuHalfPi0[n] != 0){
				valueNLOEtaToPi0NLOMuHalf[n] = valueNLOMuHalfEta[n]/valueNLOMuHalfPi0[n];
			} else {
				valueNLOEtaToPi0NLOMuHalf[n] = 0.;
			}
		}
		if (valueNLOMuOnePi0[n] != 0){
			valueNLOEtaToPi0NLOMuOne[n] = valueNLOMuOneEta[n]/valueNLOMuOnePi0[n];
		} else {
			valueNLOEtaToPi0NLOMuOne[n] = 0.;
		}
		if (valueNLOMuTwoPi0[n] != 0){
			valueNLOEtaToPi0NLOMuTwo[n] = valueNLOMuTwoEta[n]/valueNLOMuTwoPi0[n];
		} else {
			valueNLOEtaToPi0NLOMuTwo[n] = 0.;
		}
	}  
	TGraph* graphEtaToPi0NLOMuHalf7TeV =  new TGraph(xNBins,xValueNLO,valueNLOEtaToPi0NLOMuHalf);  
	graphEtaToPi0NLOMuHalf7TeV->RemovePoint(0);
	TGraph* graphEtaToPi0NLOMuOne7TeV =   new TGraph(xNBins,xValueNLO,valueNLOEtaToPi0NLOMuOne);   
	TGraph* graphEtaToPi0NLOMuTwo7TeV =   new TGraph(xNBins,xValueNLO,valueNLOEtaToPi0NLOMuTwo);   
	
	Double_t* valueNLOMuHalfEta2760GeV =   graphNLOCalcInvSecEtaMuHalf2760GeV->GetY();
	Double_t* valueNLOMuOneEta2760GeV =    graphNLOCalcInvSecEtaMuOne2760GeV->GetY();
	Double_t* valueNLOMuTwoEta2760GeV =    graphNLOCalcInvSecEtaMuTwo2760GeV->GetY();
	Double_t* valueNLOMuHalfPi02760GeV =   graphNLOCalcInvSecPi0MuHalf2760GeV->GetY();
	Double_t* valueNLOMuOnePi02760GeV =    graphNLOCalcInvSecPi0MuOne2760GeV->GetY();
	Double_t* valueNLOMuTwoPi02760GeV =    graphNLOCalcInvSecPi0MuTwo2760GeV->GetY();
	Double_t* xValueNLO2760GeV =        graphNLOCalcInvSecPi0MuOne2760GeV->GetX();
	Int_t    xNBins2760GeV =         graphNLOCalcInvSecPi0MuOne2760GeV->GetN();
	Double_t    valueNLOEtaToPi0NLOMuHalf2760GeV[100];
	Double_t    valueNLOEtaToPi0NLOMuOne2760GeV[100];
	Double_t    valueNLOEtaToPi0NLOMuTwo2760GeV[100];
	
	for ( Int_t n = 0; n < xNBins2760GeV+1; n++){
		//cout << valueNLOMuHalfPi02760GeV[n]  << "\t" << valueNLOMuOnePi02760GeV[n] << "\t" << valueNLOMuTwoPi02760GeV[n]<<endl;
		if (n == 0){
			valueNLOEtaToPi0NLOMuHalf2760GeV[n] = 0.;
		} else { 
			if (valueNLOMuHalfPi0[n] != 0){
				valueNLOEtaToPi0NLOMuHalf2760GeV[n] = valueNLOMuHalfEta2760GeV[n]/valueNLOMuHalfPi02760GeV[n];
			} else {
				valueNLOEtaToPi0NLOMuHalf2760GeV[n] = 0.;
			}
		}
		if (valueNLOMuOnePi02760GeV[n] != 0){
			valueNLOEtaToPi0NLOMuOne2760GeV[n] = valueNLOMuOneEta2760GeV[n]/valueNLOMuOnePi02760GeV[n];
		} else {
			valueNLOEtaToPi0NLOMuOne2760GeV[n] = 0.;
		}
		if (valueNLOMuTwoPi02760GeV[n] != 0){
			valueNLOEtaToPi0NLOMuTwo2760GeV[n] = valueNLOMuTwoEta2760GeV[n]/valueNLOMuTwoPi02760GeV[n];
		} else {
			valueNLOEtaToPi0NLOMuTwo2760GeV[n] = 0.;
		}
	}  
	TGraph* graphEtaToPi0NLOMuHalf2760GeV =  new TGraph(xNBins2760GeV,xValueNLO2760GeV,valueNLOEtaToPi0NLOMuHalf2760GeV);  
	graphEtaToPi0NLOMuHalf2760GeV->RemovePoint(0);
	TGraph* graphEtaToPi0NLOMuOne2760GeV =   new TGraph(xNBins2760GeV,xValueNLO2760GeV,valueNLOEtaToPi0NLOMuOne2760GeV); 
	TGraph* graphEtaToPi0NLOMuTwo2760GeV =   new TGraph(xNBins2760GeV,xValueNLO2760GeV,valueNLOEtaToPi0NLOMuTwo2760GeV); 

	Double_t* valueNLOMuHalfEta900GeV =   graphNLOCalcInvSecEtaMuHalf900GeV->GetY();
	Double_t* valueNLOMuOneEta900GeV =    graphNLOCalcInvSecEtaMuOne900GeV->GetY();
	Double_t* valueNLOMuTwoEta900GeV =    graphNLOCalcInvSecEtaMuTwo900GeV->GetY();
	Double_t* valueNLOMuHalfPi0900GeV =   graphNLOCalcInvSecPi0MuHalf900GeV->GetY();
	Double_t* valueNLOMuOnePi0900GeV =    graphNLOCalcInvSecPi0MuOne900GeV->GetY();
	Double_t* valueNLOMuTwoPi0900GeV =    graphNLOCalcInvSecPi0MuTwo900GeV->GetY();
	Double_t* xValueNLO900GeV =        graphNLOCalcInvSecPi0MuOne900GeV->GetX();
	Int_t    xNBins900GeV =         18;
	Double_t    valueNLOEtaToPi0NLOMuHalf900GeV[100];
	Double_t    valueNLOEtaToPi0NLOMuOne900GeV[100];
	Double_t    valueNLOEtaToPi0NLOMuTwo900GeV[100];

	cout << "eta" << endl;
	graphNLOCalcInvSecEtaMuOne900GeV->Print();
	cout << "pi0" << endl;
	graphNLOCalcInvSecPi0MuOne900GeV->Print();


	
	for ( Int_t n = 0; n < xNBins900GeV+1; n++){
		//cout << valueNLOMuHalfPi0900GeV[n]  << "\t" << valueNLOMuOnePi0900GeV[n] << "\t" << valueNLOMuTwoPi0900GeV[n]<<endl;
		if (n == 0){
			valueNLOEtaToPi0NLOMuHalf900GeV[n] = 0.;
		} else { 
			if (valueNLOMuHalfPi0900GeV[n] != 0){
				
				valueNLOEtaToPi0NLOMuHalf900GeV[n] = valueNLOMuHalfEta900GeV[n]/valueNLOMuHalfPi0900GeV[n];
			} else {
				valueNLOEtaToPi0NLOMuHalf900GeV[n] = 0.;
			}
		}
		if (valueNLOMuOnePi0900GeV[n] != 0){
			valueNLOEtaToPi0NLOMuOne900GeV[n] = valueNLOMuOneEta900GeV[n]/valueNLOMuOnePi0900GeV[n];
		} else {
			valueNLOEtaToPi0NLOMuOne900GeV[n] = 0.;
		}
		if (valueNLOMuTwoPi0900GeV[n] != 0){
			valueNLOEtaToPi0NLOMuTwo900GeV[n] = valueNLOMuTwoEta900GeV[n]/valueNLOMuTwoPi0900GeV[n];
		} else {
			valueNLOEtaToPi0NLOMuTwo900GeV[n] = 0.;
		}
	}  
	TGraph* graphEtaToPi0NLOMuHalf900GeV =  new TGraph(xNBins900GeV,xValueNLO900GeV,valueNLOEtaToPi0NLOMuHalf900GeV);  
	graphEtaToPi0NLOMuHalf900GeV->RemovePoint(0);
	TGraph* graphEtaToPi0NLOMuOne900GeV =   new TGraph(xNBins900GeV,xValueNLO900GeV,valueNLOEtaToPi0NLOMuOne900GeV); 
	TGraph* graphEtaToPi0NLOMuTwo900GeV =   new TGraph(xNBins900GeV,xValueNLO900GeV,valueNLOEtaToPi0NLOMuTwo900GeV); 
	

	Double_t       ptNLOEta8000GeV[100];
	Double_t       muHalfEta8000GeV[100];
	Double_t       muOneEta8000GeV[100];
	Double_t       muTwoEta8000GeV[100];
	Int_t       nlinesNLOEta8000GeV =       0;
	
	TString fileNameNLOEta8000GeV = "ExternalInput/Theory/ALICENLOcalcEtaVogelsang8Tev.dat";
	ifstream  fileNLOEta8000GeV;
	fileNLOEta8000GeV.open(fileNameNLOEta8000GeV,ios_base::in);
	cout << fileNameNLOEta8000GeV << endl;
	
	while(!fileNLOEta8000GeV.eof() && nlinesNLOEta8000GeV < 100){
		nlinesNLOEta8000GeV++;
		fileNLOEta8000GeV >> ptNLOEta8000GeV[nlinesNLOEta8000GeV] >> muHalfEta8000GeV[nlinesNLOEta8000GeV] >> muOneEta8000GeV[nlinesNLOEta8000GeV] >> muTwoEta8000GeV[nlinesNLOEta8000GeV]; 
		cout << nlinesNLOEta8000GeV << "         "  << ptNLOEta8000GeV[nlinesNLOEta8000GeV] << "         "  << muHalfEta8000GeV[nlinesNLOEta8000GeV] << "         "  << muOneEta8000GeV[nlinesNLOEta8000GeV] << "         "  << muTwoEta8000GeV[nlinesNLOEta8000GeV] << endl;;
	}
	fileNLOEta8000GeV.close();
	TGraph* graphNLOCalcInvSecEtaMuHalf8000GeV = new TGraph(nlinesNLOEta8000GeV,ptNLOEta8000GeV,muHalfEta8000GeV); 
	graphNLOCalcInvSecEtaMuHalf8000GeV->RemovePoint(0);
	graphNLOCalcInvSecEtaMuHalf8000GeV->RemovePoint(0);
	TGraph* graphNLOCalcInvSecEtaMuOne8000GeV = new TGraph(nlinesNLOEta8000GeV,ptNLOEta8000GeV,muOneEta8000GeV); 
	graphNLOCalcInvSecEtaMuOne8000GeV->RemovePoint(0);
	TGraph* graphNLOCalcInvSecEtaMuTwo8000GeV = new TGraph(nlinesNLOEta8000GeV,ptNLOEta8000GeV,muTwoEta8000GeV); 
	graphNLOCalcInvSecEtaMuTwo8000GeV->RemovePoint(0);

	TGraph* graphNLOCalcInvYieldEtaMuHalf8000GeV = ScaleGraph(graphNLOCalcInvSecEtaMuHalf8000GeV, 1/(xSection8000GeV*recalcBarn));
	TGraph* graphNLOCalcInvYieldEtaMuOne8000GeV =  ScaleGraph(graphNLOCalcInvSecEtaMuOne8000GeV, 1/(xSection8000GeV*recalcBarn));
	TGraph* graphNLOCalcInvYieldEtaMuTwo8000GeV =  ScaleGraph(graphNLOCalcInvSecEtaMuTwo8000GeV, 1/(xSection8000GeV*recalcBarn));
	
	Double_t       ptNLOPi08000GeV[100];
	Double_t       muHalfPi08000GeV[100];
	Double_t       muOnePi08000GeV[100];
	Double_t       muTwoPi08000GeV[100];
	Int_t       nlinesNLOPi08000GeV =       0;
	
	TString fileNameNLOPi08000GeV = "ExternalInput/Theory/ALICENLOcalcPi0Vogelsang8Tev.dat";
	ifstream  fileNLOPi08000GeV;
	fileNLOPi08000GeV.open(fileNameNLOPi08000GeV,ios_base::in);
	cout << fileNameNLOPi08000GeV << endl;
	
	while(!fileNLOPi08000GeV.eof() && nlinesNLOPi08000GeV< 150){
		nlinesNLOPi08000GeV++;
		fileNLOPi08000GeV >> ptNLOPi08000GeV[nlinesNLOPi08000GeV] >> muHalfPi08000GeV[nlinesNLOPi08000GeV] >> muOnePi08000GeV[nlinesNLOPi08000GeV] >> muTwoPi08000GeV[nlinesNLOPi08000GeV]; 
		cout << nlinesNLOPi08000GeV << "         "  << ptNLOPi08000GeV[nlinesNLOPi08000GeV] << "         "  << muHalfPi08000GeV[nlinesNLOPi08000GeV] << "         "  << muOnePi08000GeV[nlinesNLOPi08000GeV] << "         "  << muTwoPi08000GeV[nlinesNLOPi08000GeV] << endl;;
	}
	fileNLOPi08000GeV.close();
	TGraph* graphNLOCalcInvSecPi0MuHalf8000GeV = new TGraph(nlinesNLOPi08000GeV,ptNLOPi08000GeV,muHalfPi08000GeV); 
	graphNLOCalcInvSecPi0MuHalf8000GeV->RemovePoint(0);
	graphNLOCalcInvSecPi0MuHalf8000GeV->RemovePoint(0);
	TGraph* graphNLOCalcInvSecPi0MuOne8000GeV = new TGraph(nlinesNLOPi08000GeV,ptNLOPi08000GeV,muOnePi08000GeV); 
	graphNLOCalcInvSecPi0MuOne8000GeV->RemovePoint(0);
	TGraph* graphNLOCalcInvSecPi0MuTwo8000GeV = new TGraph(nlinesNLOPi08000GeV,ptNLOPi08000GeV,muTwoPi08000GeV); 
	graphNLOCalcInvSecPi0MuTwo8000GeV->RemovePoint(0);
	
	TGraph* graphNLOCalcInvYieldPi0MuHalf8000GeV = ScaleGraph(graphNLOCalcInvSecPi0MuHalf8000GeV, 1/(xSection8000GeV*recalcBarn));
	TGraph* graphNLOCalcInvYieldPi0MuOne8000GeV =  ScaleGraph(graphNLOCalcInvSecPi0MuOne8000GeV, 1/(xSection8000GeV*recalcBarn));
	TGraph* graphNLOCalcInvYieldPi0MuTwo8000GeV =  ScaleGraph(graphNLOCalcInvSecPi0MuTwo8000GeV, 1/(xSection8000GeV*recalcBarn));

	Double_t* valueNLOMuHalfEta8000GeV = graphNLOCalcInvSecEtaMuHalf8000GeV->GetY();
	Double_t* valueNLOMuOneEta8000GeV =  graphNLOCalcInvSecEtaMuOne8000GeV->GetY();
	Double_t* valueNLOMuTwoEta8000GeV =  graphNLOCalcInvSecEtaMuTwo8000GeV->GetY();
	Double_t* valueNLOMuHalfPi08000GeV = graphNLOCalcInvSecPi0MuHalf8000GeV->GetY();
	Double_t* valueNLOMuOnePi08000GeV =  graphNLOCalcInvSecPi0MuOne8000GeV->GetY();
	Double_t* valueNLOMuTwoPi08000GeV =  graphNLOCalcInvSecPi0MuTwo8000GeV->GetY();
	Double_t* xValueNLO8000GeV =      graphNLOCalcInvSecPi0MuOne8000GeV->GetX();
	Int_t    xNBins8000GeV =          graphNLOCalcInvSecPi0MuOne8000GeV->GetN();
	Double_t    valueNLOEtaToPi0NLOMuHalf8000GeV[100];
	Double_t    valueNLOEtaToPi0NLOMuOne8000GeV[100];
	Double_t    valueNLOEtaToPi0NLOMuTwo8000GeV[100];
	
	for ( Int_t n = 0; n < xNBins8000GeV+1; n++){
		if (n == 0){
			valueNLOEtaToPi0NLOMuHalf8000GeV[n] = 0.;
		} else { 
			if (valueNLOMuHalfPi08000GeV[n] != 0){
				valueNLOEtaToPi0NLOMuHalf8000GeV[n] = valueNLOMuHalfEta8000GeV[n]/valueNLOMuHalfPi08000GeV[n];
			} else {
				valueNLOEtaToPi0NLOMuHalf8000GeV[n] = 0.;
			}
		}
		if (valueNLOMuOnePi08000GeV[n] != 0){
			valueNLOEtaToPi0NLOMuOne8000GeV[n] = valueNLOMuOneEta8000GeV[n]/valueNLOMuOnePi08000GeV[n];
		} else {
			valueNLOEtaToPi0NLOMuOne8000GeV[n] = 0.;
		}
		if (valueNLOMuTwoPi08000GeV[n] != 0){
			valueNLOEtaToPi0NLOMuTwo8000GeV[n] = valueNLOMuTwoEta8000GeV[n]/valueNLOMuTwoPi08000GeV[n];
		} else {
			valueNLOEtaToPi0NLOMuTwo8000GeV[n] = 0.;
		}
	}  
	TGraph* graphEtaToPi0NLOMuHalf8TeV =  new TGraph(xNBins8000GeV,xValueNLO8000GeV,valueNLOEtaToPi0NLOMuHalf8000GeV);  
	graphEtaToPi0NLOMuHalf8TeV->RemovePoint(0);
	graphEtaToPi0NLOMuHalf8TeV->Print();
	TGraph* graphEtaToPi0NLOMuOne8TeV =   new TGraph(xNBins8000GeV,xValueNLO8000GeV,valueNLOEtaToPi0NLOMuOne8000GeV);   
	graphEtaToPi0NLOMuOne8TeV->Print();
	TGraph* graphEtaToPi0NLOMuTwo8TeV =   new TGraph(xNBins8000GeV,xValueNLO8000GeV,valueNLOEtaToPi0NLOMuTwo8000GeV);   
		graphEtaToPi0NLOMuTwo8TeV->Print();
	
	TFile* file2760GeV_Pythia8_fixedBinning = TFile::Open("ExternalInput/Theory/Pythia/pp_pi0_2760GeV_pythia8_tune4C_comp_fixed_binning_2014-02-04-1330.root");
	TH1F* histoPythia8Spec2760GeV = (TH1F*)file2760GeV_Pythia8_fixedBinning->Get("hPythiaComb");
	TH1F* histoPythia8InvYield2760GeV = (TH1F*)histoPythia8Spec2760GeV->Clone("histoPythia8InvYield2760GeV");
	histoPythia8InvYield2760GeV->Scale(1./(xSection2760GeV*recalcBarn));
	
	TFile* filePythia8_2760GeV = TFile::Open("ExternalInput/Theory/Pythia/pp_pi0_2760GeV_pythia8_tune4C_comp_variable_binning.root");
	TH1F* histoPythia8Spec2760GeVVarBinning = (TH1F*)filePythia8_2760GeV->Get("hPythiaComb");
	TH1F* histoPythia8InvYield2760GeVVarBinning = (TH1F*)histoPythia8Spec2760GeVVarBinning->Clone("histoPythia8InvYield2760GeV");
	histoPythia8InvYield2760GeVVarBinning->Scale(1./(xSection2760GeV*recalcBarn));

	
	//*******************************************************************************************************
	//****************************** Pythia 8 Tune 4C *******************************************************
	//*******************************************************************************************************
	Color_t colorMB         = kBlack;
	Color_t colorptHard2     = kBlue+1;
	Color_t colorptHard5     = kRed+1;
	Color_t colorptHard10     = kGreen+1;
	
	Double_t scalingFactor2760GeVPy8_4C_MB                = 1.096e3/(2361698525/4e8)*1/4e8*1e9;
	Double_t scalingFactor2760GeVPy8_4C_ptHard2            = 2.690e2/(671454494/1e8)*1/1e8*1e9;
	Double_t scalingFactor2760GeVPy8_4C_ptHard5            = 1.522e1/(702687699/1e8)*1/1e8*1e9;
	Double_t scalingFactor2760GeVPy8_4C_ptHard10        = 1.194e0/(659556684/8.6e7)*1/8.6e7*1e9;
	
//     Double_t xSection2760GeVPy8_4CMB                    = 1.09e3;
//     Double_t xSection2760GeVPy8_4CMBErr                    = 3.077e-2;
//     Double_t nEvt2760GeVPy8_4CMB                        = 400e6;
//     Double_t xSection2760GeVPy8_4CptHard5                = 1.522e1;
//     Double_t xSection2760GeVPy8_4CptHard5Err            = 8.484e-4;
//     Double_t nEvt2760GeVPy8_4CptHard5                    = 100e6;
//     Double_t xSection2760GeVPy8_4CptHard10                = 1.194e0;
//     Double_t xSection2760GeVPy8_4CptHard10Err            = 6.968e-5;
//     Double_t nEvt2760GeVPy8_4CptHard10                    = 86e6;
	
	TFile*     file2760GeV_Pythia8_Tune4CMB                 = TFile::Open("ExternalInput/Theory/Pythia/hist_Pythia8_1_Tune4C_MB_400Mio.root");
	TH1F*     histoPi02760GeV_Pythia8_Tune4CMB             = (TH1F*)file2760GeV_Pythia8_Tune4CMB->Get("ptPi0W0K");
	TH1F*     histoPi0FromEta2760GeV_Pythia8_Tune4CMB     = (TH1F*)file2760GeV_Pythia8_Tune4CMB->Get("ptPi0FromEta");
	TH1F*     histoPi0FromLambda2760GeV_Pythia8_Tune4CMB     = (TH1F*)file2760GeV_Pythia8_Tune4CMB->Get("ptPi0FromLambda");
	TH1F*     histoPi0FromK0s2760GeV_Pythia8_Tune4CMB     = (TH1F*)file2760GeV_Pythia8_Tune4CMB->Get("ptPi0");
	histoPi0FromK0s2760GeV_Pythia8_Tune4CMB->Sumw2();
	histoPi0FromK0s2760GeV_Pythia8_Tune4CMB->Add(histoPi02760GeV_Pythia8_Tune4CMB,-1);
	TH1F*     histoEta2760GeV_Pythia8_Tune4CMB             = (TH1F*)file2760GeV_Pythia8_Tune4CMB->Get("ptEta");
	TH1F*     histoPiPl2760GeV_Pythia8_Tune4CMB             = (TH1F*)file2760GeV_Pythia8_Tune4CMB->Get("ptPiPlW0K");
	TH1F*     histoPiMi2760GeV_Pythia8_Tune4CMB             = (TH1F*)file2760GeV_Pythia8_Tune4CMB->Get("ptPiMiW0K");
	TH1F*     histoPi0UB2760GeV_Pythia8_Tune4CMB             = (TH1F*)file2760GeV_Pythia8_Tune4CMB->Get("ptPi0W0KUB");
	TH1F*     histoPi0FromEtaUB2760GeV_Pythia8_Tune4CMB     = (TH1F*)file2760GeV_Pythia8_Tune4CMB->Get("ptPi0FromEtaUB");
	TH1F*     histoPi0FromLambdaUB2760GeV_Pythia8_Tune4CMB     = (TH1F*)file2760GeV_Pythia8_Tune4CMB->Get("ptPi0FromLambdaUB");
	TH1F*     histoPi0FromK0sUB2760GeV_Pythia8_Tune4CMB     = (TH1F*)file2760GeV_Pythia8_Tune4CMB->Get("ptPi0UB");
	histoPi0FromK0sUB2760GeV_Pythia8_Tune4CMB->Sumw2();
	histoPi0FromK0sUB2760GeV_Pythia8_Tune4CMB->Add(histoPi0UB2760GeV_Pythia8_Tune4CMB,-1);
	TH1F*     histoEtaUB2760GeV_Pythia8_Tune4CMB             = (TH1F*)file2760GeV_Pythia8_Tune4CMB->Get("ptEtaUB");
	TH1F*     histoPiPlUB2760GeV_Pythia8_Tune4CMB             = (TH1F*)file2760GeV_Pythia8_Tune4CMB->Get("ptPiPlW0KUB");
	TH1F*     histoPiMiUB2760GeV_Pythia8_Tune4CMB             = (TH1F*)file2760GeV_Pythia8_Tune4CMB->Get("ptPiMiW0KUB");
	
	histoPi02760GeV_Pythia8_Tune4CMB->Scale(scalingFactor2760GeVPy8_4C_MB);
	histoPi0FromLambda2760GeV_Pythia8_Tune4CMB->Scale(scalingFactor2760GeVPy8_4C_MB);
	histoPi0FromEta2760GeV_Pythia8_Tune4CMB->Scale(scalingFactor2760GeVPy8_4C_MB);
	histoPi0FromK0s2760GeV_Pythia8_Tune4CMB->Scale(scalingFactor2760GeVPy8_4C_MB);
	histoEta2760GeV_Pythia8_Tune4CMB->Scale(scalingFactor2760GeVPy8_4C_MB);
	histoPiPl2760GeV_Pythia8_Tune4CMB->Scale(scalingFactor2760GeVPy8_4C_MB);
	histoPiMi2760GeV_Pythia8_Tune4CMB->Scale(scalingFactor2760GeVPy8_4C_MB);
	histoPi0UB2760GeV_Pythia8_Tune4CMB->Scale(scalingFactor2760GeVPy8_4C_MB);
	histoPi0FromLambdaUB2760GeV_Pythia8_Tune4CMB->Scale(scalingFactor2760GeVPy8_4C_MB);
	histoPi0FromEtaUB2760GeV_Pythia8_Tune4CMB->Scale(scalingFactor2760GeVPy8_4C_MB);
	histoPi0FromK0sUB2760GeV_Pythia8_Tune4CMB->Scale(scalingFactor2760GeVPy8_4C_MB);
	histoEtaUB2760GeV_Pythia8_Tune4CMB->Scale(scalingFactor2760GeVPy8_4C_MB);
	histoPiPlUB2760GeV_Pythia8_Tune4CMB->Scale(scalingFactor2760GeVPy8_4C_MB);
	histoPiMiUB2760GeV_Pythia8_Tune4CMB->Scale(scalingFactor2760GeVPy8_4C_MB);

	TFile*     file2760GeV_Pythia8_Tune4CptHard2                 = TFile::Open("ExternalInput/Theory/Pythia/hist_Pythia8_1_Tune4C_pTHard2GeV_100Mio.root");
	TH1F*     histoPi02760GeV_Pythia8_Tune4CptHard2             = (TH1F*)file2760GeV_Pythia8_Tune4CptHard2->Get("ptPi0W0K");
	TH1F*     histoPi0FromEta2760GeV_Pythia8_Tune4CptHard2     = (TH1F*)file2760GeV_Pythia8_Tune4CptHard2->Get("ptPi0FromEta");
	TH1F*     histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard2 = (TH1F*)file2760GeV_Pythia8_Tune4CptHard2->Get("ptPi0FromLambda");
	TH1F*     histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard2     = (TH1F*)file2760GeV_Pythia8_Tune4CptHard2->Get("ptPi0");
	histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard2->Sumw2();
	histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard2->Add(histoPi02760GeV_Pythia8_Tune4CptHard2,-1);
	TH1F*     histoEta2760GeV_Pythia8_Tune4CptHard2             = (TH1F*)file2760GeV_Pythia8_Tune4CptHard2->Get("ptEta");
	TH1F*     histoPiPl2760GeV_Pythia8_Tune4CptHard2             = (TH1F*)file2760GeV_Pythia8_Tune4CptHard2->Get("ptPiPlW0K");
	TH1F*     histoPiMi2760GeV_Pythia8_Tune4CptHard2             = (TH1F*)file2760GeV_Pythia8_Tune4CptHard2->Get("ptPiMiW0K");
	TH1F*     histoPi0UB2760GeV_Pythia8_Tune4CptHard2         = (TH1F*)file2760GeV_Pythia8_Tune4CptHard2->Get("ptPi0W0KUB");
	TH1F*     histoPi0FromEtaUB2760GeV_Pythia8_Tune4CptHard2     = (TH1F*)file2760GeV_Pythia8_Tune4CptHard2->Get("ptPi0FromEtaUB");
	TH1F*     histoPi0FromLambdaUB2760GeV_Pythia8_Tune4CptHard2     = (TH1F*)file2760GeV_Pythia8_Tune4CptHard2->Get("ptPi0FromLambdaUB");
	TH1F*     histoPi0FromK0sUB2760GeV_Pythia8_Tune4CptHard2     = (TH1F*)file2760GeV_Pythia8_Tune4CptHard2->Get("ptPi0UB");
	histoPi0FromK0sUB2760GeV_Pythia8_Tune4CptHard2->Sumw2();
	histoPi0FromK0sUB2760GeV_Pythia8_Tune4CptHard2->Add(histoPi0UB2760GeV_Pythia8_Tune4CptHard2,-1);
	TH1F*     histoEtaUB2760GeV_Pythia8_Tune4CptHard2         = (TH1F*)file2760GeV_Pythia8_Tune4CptHard2->Get("ptEtaUB");
	TH1F*     histoPiPlUB2760GeV_Pythia8_Tune4CptHard2         = (TH1F*)file2760GeV_Pythia8_Tune4CptHard2->Get("ptPiPlW0KUB");
	TH1F*     histoPiMiUB2760GeV_Pythia8_Tune4CptHard2         = (TH1F*)file2760GeV_Pythia8_Tune4CptHard2->Get("ptPiMiW0KUB");
	
	histoPi02760GeV_Pythia8_Tune4CptHard2->Scale(scalingFactor2760GeVPy8_4C_ptHard2);
	histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard2->Scale(scalingFactor2760GeVPy8_4C_ptHard2);
	histoPi0FromEta2760GeV_Pythia8_Tune4CptHard2->Scale(scalingFactor2760GeVPy8_4C_ptHard2);
	histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard2->Scale(scalingFactor2760GeVPy8_4C_ptHard2);
	histoEta2760GeV_Pythia8_Tune4CptHard2->Scale(scalingFactor2760GeVPy8_4C_ptHard2);
	histoPiPl2760GeV_Pythia8_Tune4CptHard2->Scale(scalingFactor2760GeVPy8_4C_ptHard2);
	histoPiMi2760GeV_Pythia8_Tune4CptHard2->Scale(scalingFactor2760GeVPy8_4C_ptHard2);
	histoPi0UB2760GeV_Pythia8_Tune4CptHard2->Scale(scalingFactor2760GeVPy8_4C_ptHard2);
	histoPi0FromLambdaUB2760GeV_Pythia8_Tune4CptHard2->Scale(scalingFactor2760GeVPy8_4C_ptHard2);
	histoPi0FromEtaUB2760GeV_Pythia8_Tune4CptHard2->Scale(scalingFactor2760GeVPy8_4C_ptHard2);
	histoPi0FromK0sUB2760GeV_Pythia8_Tune4CptHard2->Scale(scalingFactor2760GeVPy8_4C_ptHard2);
	histoEtaUB2760GeV_Pythia8_Tune4CptHard2->Scale(scalingFactor2760GeVPy8_4C_ptHard2);
	histoPiPlUB2760GeV_Pythia8_Tune4CptHard2->Scale(scalingFactor2760GeVPy8_4C_ptHard2);
	histoPiMiUB2760GeV_Pythia8_Tune4CptHard2->Scale(scalingFactor2760GeVPy8_4C_ptHard2);
	
	
	TFile*     file2760GeV_Pythia8_Tune4CptHard5                 = TFile::Open("ExternalInput/Theory/Pythia/hist_Pythia8_1_Tune4C_pTHard5GeV_100Mio.root");
	TH1F*     histoPi02760GeV_Pythia8_Tune4CptHard5             = (TH1F*)file2760GeV_Pythia8_Tune4CptHard5->Get("ptPi0W0K");
	TH1F*     histoPi0FromEta2760GeV_Pythia8_Tune4CptHard5     = (TH1F*)file2760GeV_Pythia8_Tune4CptHard5->Get("ptPi0FromEta");
	TH1F*     histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard5 = (TH1F*)file2760GeV_Pythia8_Tune4CptHard5->Get("ptPi0FromLambda");
	TH1F*     histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard5     = (TH1F*)file2760GeV_Pythia8_Tune4CptHard5->Get("ptPi0");
	histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard5->Sumw2();
	histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard5->Add(histoPi02760GeV_Pythia8_Tune4CptHard5,-1);
	TH1F*     histoEta2760GeV_Pythia8_Tune4CptHard5             = (TH1F*)file2760GeV_Pythia8_Tune4CptHard5->Get("ptEta");
	TH1F*     histoPiPl2760GeV_Pythia8_Tune4CptHard5             = (TH1F*)file2760GeV_Pythia8_Tune4CptHard5->Get("ptPiPlW0K");
	TH1F*     histoPiMi2760GeV_Pythia8_Tune4CptHard5             = (TH1F*)file2760GeV_Pythia8_Tune4CptHard5->Get("ptPiMiW0K");
	TH1F*     histoPi0UB2760GeV_Pythia8_Tune4CptHard5         = (TH1F*)file2760GeV_Pythia8_Tune4CptHard5->Get("ptPi0W0KUB");
	TH1F*     histoPi0FromEtaUB2760GeV_Pythia8_Tune4CptHard5     = (TH1F*)file2760GeV_Pythia8_Tune4CptHard5->Get("ptPi0FromEtaUB");
	TH1F*     histoPi0FromLambdaUB2760GeV_Pythia8_Tune4CptHard5     = (TH1F*)file2760GeV_Pythia8_Tune4CptHard5->Get("ptPi0FromLambdaUB");
	TH1F*     histoPi0FromK0sUB2760GeV_Pythia8_Tune4CptHard5     = (TH1F*)file2760GeV_Pythia8_Tune4CptHard5->Get("ptPi0UB");
	histoPi0FromK0sUB2760GeV_Pythia8_Tune4CptHard5->Sumw2();
	histoPi0FromK0sUB2760GeV_Pythia8_Tune4CptHard5->Add(histoPi0UB2760GeV_Pythia8_Tune4CptHard5,-1);
	TH1F*     histoEtaUB2760GeV_Pythia8_Tune4CptHard5         = (TH1F*)file2760GeV_Pythia8_Tune4CptHard5->Get("ptEtaUB");
	TH1F*     histoPiPlUB2760GeV_Pythia8_Tune4CptHard5         = (TH1F*)file2760GeV_Pythia8_Tune4CptHard5->Get("ptPiPlW0KUB");
	TH1F*     histoPiMiUB2760GeV_Pythia8_Tune4CptHard5         = (TH1F*)file2760GeV_Pythia8_Tune4CptHard5->Get("ptPiMiW0KUB");
	
	histoPi02760GeV_Pythia8_Tune4CptHard5->Scale(scalingFactor2760GeVPy8_4C_ptHard5);
	histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard5->Scale(scalingFactor2760GeVPy8_4C_ptHard5);
	histoPi0FromEta2760GeV_Pythia8_Tune4CptHard5->Scale(scalingFactor2760GeVPy8_4C_ptHard5);
	histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard5->Scale(scalingFactor2760GeVPy8_4C_ptHard5);
	histoEta2760GeV_Pythia8_Tune4CptHard5->Scale(scalingFactor2760GeVPy8_4C_ptHard5);
	histoPiPl2760GeV_Pythia8_Tune4CptHard5->Scale(scalingFactor2760GeVPy8_4C_ptHard5);
	histoPiMi2760GeV_Pythia8_Tune4CptHard5->Scale(scalingFactor2760GeVPy8_4C_ptHard5);
	histoPi0UB2760GeV_Pythia8_Tune4CptHard5->Scale(scalingFactor2760GeVPy8_4C_ptHard5);
	histoPi0FromLambdaUB2760GeV_Pythia8_Tune4CptHard5->Scale(scalingFactor2760GeVPy8_4C_ptHard5);
	histoPi0FromEtaUB2760GeV_Pythia8_Tune4CptHard5->Scale(scalingFactor2760GeVPy8_4C_ptHard5);
	histoPi0FromK0sUB2760GeV_Pythia8_Tune4CptHard5->Scale(scalingFactor2760GeVPy8_4C_ptHard5);
	histoEtaUB2760GeV_Pythia8_Tune4CptHard5->Scale(scalingFactor2760GeVPy8_4C_ptHard5);
	histoPiPlUB2760GeV_Pythia8_Tune4CptHard5->Scale(scalingFactor2760GeVPy8_4C_ptHard5);
	histoPiMiUB2760GeV_Pythia8_Tune4CptHard5->Scale(scalingFactor2760GeVPy8_4C_ptHard5);

	TFile*     file2760GeV_Pythia8_Tune4CptHard10                 = TFile::Open("ExternalInput/Theory/Pythia/hist_Pythia8_1_Tune4C_ptHard10GeV_Backup86.0Mio.root");
	TH1F*     histoPi02760GeV_Pythia8_Tune4CptHard10             = (TH1F*)file2760GeV_Pythia8_Tune4CptHard10->Get("ptPi0W0KClone");
	TH1F*     histoPi0FromEta2760GeV_Pythia8_Tune4CptHard10     = (TH1F*)file2760GeV_Pythia8_Tune4CptHard10->Get("ptPi0FromEtaClone");
	TH1F*     histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard10     = (TH1F*)file2760GeV_Pythia8_Tune4CptHard10->Get("ptPi0FromLambdaClone");
	TH1F*     histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard10     = (TH1F*)file2760GeV_Pythia8_Tune4CptHard10->Get("ptPi0Clone");
	histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard10->Sumw2();
	histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard10->Add(histoPi02760GeV_Pythia8_Tune4CptHard10,-1);
	TH1F*     histoEta2760GeV_Pythia8_Tune4CptHard10             = (TH1F*)file2760GeV_Pythia8_Tune4CptHard10->Get("ptEtaClone");
	TH1F*     histoPiPl2760GeV_Pythia8_Tune4CptHard10             = (TH1F*)file2760GeV_Pythia8_Tune4CptHard10->Get("ptPiPlW0KClone");
	TH1F*     histoPiMi2760GeV_Pythia8_Tune4CptHard10             = (TH1F*)file2760GeV_Pythia8_Tune4CptHard10->Get("ptPiMiW0KClone");
	TH1F*     histoPi0UB2760GeV_Pythia8_Tune4CptHard10             = (TH1F*)file2760GeV_Pythia8_Tune4CptHard10->Get("ptPi0W0KUBClone");
	TH1F*     histoPi0FromEtaUB2760GeV_Pythia8_Tune4CptHard10     = (TH1F*)file2760GeV_Pythia8_Tune4CptHard10->Get("ptPi0FromEtaUBClone");
	TH1F*     histoPi0FromLambdaUB2760GeV_Pythia8_Tune4CptHard10     = (TH1F*)file2760GeV_Pythia8_Tune4CptHard10->Get("ptPi0FromLambdaUBClone");
	TH1F*     histoPi0FromK0sUB2760GeV_Pythia8_Tune4CptHard10     = (TH1F*)file2760GeV_Pythia8_Tune4CptHard10->Get("ptPi0UBClone");
	histoPi0FromK0sUB2760GeV_Pythia8_Tune4CptHard10->Sumw2();
	histoPi0FromK0sUB2760GeV_Pythia8_Tune4CptHard10->Add(histoPi0UB2760GeV_Pythia8_Tune4CptHard10,-1);
	TH1F*     histoEtaUB2760GeV_Pythia8_Tune4CptHard10             = (TH1F*)file2760GeV_Pythia8_Tune4CptHard10->Get("ptEtaUBClone");
	TH1F*     histoPiPlUB2760GeV_Pythia8_Tune4CptHard10             = (TH1F*)file2760GeV_Pythia8_Tune4CptHard10->Get("ptPiPlW0KUBClone");
	TH1F*     histoPiMiUB2760GeV_Pythia8_Tune4CptHard10             = (TH1F*)file2760GeV_Pythia8_Tune4CptHard10->Get("ptPiMiW0KUBClone");
	
	histoPi02760GeV_Pythia8_Tune4CptHard10->Scale(scalingFactor2760GeVPy8_4C_ptHard10);
	histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard10->Scale(scalingFactor2760GeVPy8_4C_ptHard10);
	histoPi0FromEta2760GeV_Pythia8_Tune4CptHard10->Scale(scalingFactor2760GeVPy8_4C_ptHard10);
	histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard10->Scale(scalingFactor2760GeVPy8_4C_ptHard10);
	histoEta2760GeV_Pythia8_Tune4CptHard10->Scale(scalingFactor2760GeVPy8_4C_ptHard10);
	histoPiPl2760GeV_Pythia8_Tune4CptHard10->Scale(scalingFactor2760GeVPy8_4C_ptHard10);
	histoPiMi2760GeV_Pythia8_Tune4CptHard10->Scale(scalingFactor2760GeVPy8_4C_ptHard10);
	histoPi0UB2760GeV_Pythia8_Tune4CptHard10->Scale(scalingFactor2760GeVPy8_4C_ptHard10);
	histoPi0FromLambdaUB2760GeV_Pythia8_Tune4CptHard10->Scale(scalingFactor2760GeVPy8_4C_ptHard10);
	histoPi0FromEtaUB2760GeV_Pythia8_Tune4CptHard10->Scale(scalingFactor2760GeVPy8_4C_ptHard10);
	histoPi0FromK0sUB2760GeV_Pythia8_Tune4CptHard10->Scale(scalingFactor2760GeVPy8_4C_ptHard10);
	histoEtaUB2760GeV_Pythia8_Tune4CptHard10->Scale(scalingFactor2760GeVPy8_4C_ptHard10);
	histoPiPlUB2760GeV_Pythia8_Tune4CptHard10->Scale(scalingFactor2760GeVPy8_4C_ptHard10);
	histoPiMiUB2760GeV_Pythia8_Tune4CptHard10->Scale(scalingFactor2760GeVPy8_4C_ptHard10);

	TH1D* ratioPi02760GeV_Pythia8_Tune4CpTHard2DivMB = (TH1D*) histoPi02760GeV_Pythia8_Tune4CptHard2->Clone("ratioPi02760GeV_Pythia8_Tune4CpTHard2DivMB");
	ratioPi02760GeV_Pythia8_Tune4CpTHard2DivMB->Divide(ratioPi02760GeV_Pythia8_Tune4CpTHard2DivMB, histoPi02760GeV_Pythia8_Tune4CMB,1.,1.,"");
	TH1D* ratioPi02760GeV_Pythia8_Tune4CpTHard5DivpTHard2 = (TH1D*) histoPi02760GeV_Pythia8_Tune4CptHard5->Clone("ratioPi02760GeV_Pythia8_Tune4CpTHard5DivpTHard2");
	ratioPi02760GeV_Pythia8_Tune4CpTHard5DivpTHard2->Divide(ratioPi02760GeV_Pythia8_Tune4CpTHard5DivpTHard2, histoPi02760GeV_Pythia8_Tune4CptHard2,1.,1.,"");
	TH1D* ratioPi02760GeV_Pythia8_Tune4CpTHard10DivpTHard5 = (TH1D*) histoPi02760GeV_Pythia8_Tune4CptHard10->Clone("ratioPi02760GeV_Pythia8_Tune4CpTHard10DivpTHard5");
	ratioPi02760GeV_Pythia8_Tune4CpTHard10DivpTHard5->Divide(ratioPi02760GeV_Pythia8_Tune4CpTHard10DivpTHard5, histoPi02760GeV_Pythia8_Tune4CptHard5,1.,1.,"");    
	TH1D* ratioPi0UB2760GeV_Pythia8_Tune4CpTHard2DivMB = (TH1D*) histoPi0UB2760GeV_Pythia8_Tune4CptHard2->Clone("ratioPi0UB2760GeV_Pythia8_Tune4CpTHard2DivMB");
	ratioPi0UB2760GeV_Pythia8_Tune4CpTHard2DivMB->Divide(ratioPi0UB2760GeV_Pythia8_Tune4CpTHard2DivMB, histoPi0UB2760GeV_Pythia8_Tune4CMB,1.,1.,"");
	TH1D* ratioPi0UB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2 = (TH1D*) histoPi0UB2760GeV_Pythia8_Tune4CptHard5->Clone("ratioPi0UB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2");
	ratioPi0UB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2->Divide(ratioPi0UB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2, histoPi0UB2760GeV_Pythia8_Tune4CptHard2,1.,1.,"");
	TH1D* ratioPi0UB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5 = (TH1D*) histoPi0UB2760GeV_Pythia8_Tune4CptHard10->Clone("ratioPi0UB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5");
	ratioPi0UB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5->Divide(ratioPi0UB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5, histoPi0UB2760GeV_Pythia8_Tune4CptHard5,1.,1.,"");    
	DrawGammaSetMarker(ratioPi02760GeV_Pythia8_Tune4CpTHard2DivMB, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(ratioPi02760GeV_Pythia8_Tune4CpTHard5DivpTHard2, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(ratioPi02760GeV_Pythia8_Tune4CpTHard10DivpTHard5, 20, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(ratioPi0UB2760GeV_Pythia8_Tune4CpTHard2DivMB, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(ratioPi0UB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(ratioPi0UB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5, 20, 1, colorptHard10, colorptHard10);

	TH1D* ratioPiPl2760GeV_Pythia8_Tune4CpTHard2DivMB = (TH1D*) histoPiPl2760GeV_Pythia8_Tune4CptHard2->Clone("ratioPiPl2760GeV_Pythia8_Tune4CpTHard2DivMB");
	ratioPiPl2760GeV_Pythia8_Tune4CpTHard2DivMB->Divide(ratioPiPl2760GeV_Pythia8_Tune4CpTHard2DivMB, histoPiPl2760GeV_Pythia8_Tune4CMB,1.,1.,"");
	TH1D* ratioPiPl2760GeV_Pythia8_Tune4CpTHard5DivpTHard2 = (TH1D*) histoPiPl2760GeV_Pythia8_Tune4CptHard5->Clone("ratioPiPl2760GeV_Pythia8_Tune4CpTHard5DivpTHard2");
	ratioPiPl2760GeV_Pythia8_Tune4CpTHard5DivpTHard2->Divide(ratioPiPl2760GeV_Pythia8_Tune4CpTHard5DivpTHard2, histoPiPl2760GeV_Pythia8_Tune4CptHard2,1.,1.,"");
	TH1D* ratioPiPl2760GeV_Pythia8_Tune4CpTHard10DivpTHard5 = (TH1D*) histoPiPl2760GeV_Pythia8_Tune4CptHard10->Clone("ratioPiPl2760GeV_Pythia8_Tune4CpTHard10DivpTHard5");
	ratioPiPl2760GeV_Pythia8_Tune4CpTHard10DivpTHard5->Divide(ratioPiPl2760GeV_Pythia8_Tune4CpTHard10DivpTHard5, histoPiPl2760GeV_Pythia8_Tune4CptHard5,1.,1.,"");    
	TH1D* ratioPiPlUB2760GeV_Pythia8_Tune4CpTHard2DivMB = (TH1D*) histoPiPlUB2760GeV_Pythia8_Tune4CptHard2->Clone("ratioPiPlUB2760GeV_Pythia8_Tune4CpTHard2DivMB");
	ratioPiPlUB2760GeV_Pythia8_Tune4CpTHard2DivMB->Divide(ratioPiPlUB2760GeV_Pythia8_Tune4CpTHard2DivMB, histoPiPlUB2760GeV_Pythia8_Tune4CMB,1.,1.,"");
	TH1D* ratioPiPlUB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2 = (TH1D*) histoPiPlUB2760GeV_Pythia8_Tune4CptHard5->Clone("ratioPiPlUB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2");
	ratioPiPlUB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2->Divide(ratioPiPlUB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2, histoPiPlUB2760GeV_Pythia8_Tune4CptHard2,1.,1.,"");
	TH1D* ratioPiPlUB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5 = (TH1D*) histoPiPlUB2760GeV_Pythia8_Tune4CptHard10->Clone("ratioPiPlUB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5");
	ratioPiPlUB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5->Divide(ratioPiPlUB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5, histoPiPlUB2760GeV_Pythia8_Tune4CptHard5,1.,1.,"");    
	DrawGammaSetMarker(ratioPiPl2760GeV_Pythia8_Tune4CpTHard2DivMB, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(ratioPiPl2760GeV_Pythia8_Tune4CpTHard5DivpTHard2, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(ratioPiPl2760GeV_Pythia8_Tune4CpTHard10DivpTHard5, 20, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(ratioPiPlUB2760GeV_Pythia8_Tune4CpTHard2DivMB, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(ratioPiPlUB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(ratioPiPlUB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5, 20, 1, colorptHard10, colorptHard10);

	TH1D* ratioPiMi2760GeV_Pythia8_Tune4CpTHard2DivMB = (TH1D*) histoPiMi2760GeV_Pythia8_Tune4CptHard2->Clone("ratioPiMi2760GeV_Pythia8_Tune4CpTHard2DivMB");
	ratioPiMi2760GeV_Pythia8_Tune4CpTHard2DivMB->Divide(ratioPiMi2760GeV_Pythia8_Tune4CpTHard2DivMB, histoPiMi2760GeV_Pythia8_Tune4CMB,1.,1.,"");
	TH1D* ratioPiMi2760GeV_Pythia8_Tune4CpTHard5DivpTHard2 = (TH1D*) histoPiMi2760GeV_Pythia8_Tune4CptHard5->Clone("ratioPiMi2760GeV_Pythia8_Tune4CpTHard5DivpTHard2");
	ratioPiMi2760GeV_Pythia8_Tune4CpTHard5DivpTHard2->Divide(ratioPiMi2760GeV_Pythia8_Tune4CpTHard5DivpTHard2, histoPiMi2760GeV_Pythia8_Tune4CptHard2,1.,1.,"");
	TH1D* ratioPiMi2760GeV_Pythia8_Tune4CpTHard10DivpTHard5 = (TH1D*) histoPiMi2760GeV_Pythia8_Tune4CptHard10->Clone("ratioPiMi2760GeV_Pythia8_Tune4CpTHard10DivpTHard5");
	ratioPiMi2760GeV_Pythia8_Tune4CpTHard10DivpTHard5->Divide(ratioPiMi2760GeV_Pythia8_Tune4CpTHard10DivpTHard5, histoPiMi2760GeV_Pythia8_Tune4CptHard5,1.,1.,"");    
	TH1D* ratioPiMiUB2760GeV_Pythia8_Tune4CpTHard2DivMB = (TH1D*) histoPiMiUB2760GeV_Pythia8_Tune4CptHard2->Clone("ratioPiMiUB2760GeV_Pythia8_Tune4CpTHard2DivMB");
	ratioPiMiUB2760GeV_Pythia8_Tune4CpTHard2DivMB->Divide(ratioPiMiUB2760GeV_Pythia8_Tune4CpTHard2DivMB, histoPiMiUB2760GeV_Pythia8_Tune4CMB,1.,1.,"");
	TH1D* ratioPiMiUB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2 = (TH1D*) histoPiMiUB2760GeV_Pythia8_Tune4CptHard5->Clone("ratioPiMiUB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2");
	ratioPiMiUB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2->Divide(ratioPiMiUB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2, histoPiMiUB2760GeV_Pythia8_Tune4CptHard2,1.,1.,"");
	TH1D* ratioPiMiUB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5 = (TH1D*) histoPiMiUB2760GeV_Pythia8_Tune4CptHard10->Clone("ratioPiMiUB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5");
	ratioPiMiUB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5->Divide(ratioPiMiUB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5, histoPiMiUB2760GeV_Pythia8_Tune4CptHard5,1.,1.,"");    
	DrawGammaSetMarker(ratioPiMi2760GeV_Pythia8_Tune4CpTHard2DivMB, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(ratioPiMi2760GeV_Pythia8_Tune4CpTHard5DivpTHard2, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(ratioPiMi2760GeV_Pythia8_Tune4CpTHard10DivpTHard5, 20, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(ratioPiMiUB2760GeV_Pythia8_Tune4CpTHard2DivMB, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(ratioPiMiUB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(ratioPiMiUB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5, 20, 1, colorptHard10, colorptHard10);

	
	TH1D* ratioEta2760GeV_Pythia8_Tune4CpTHard2DivMB = (TH1D*) histoEta2760GeV_Pythia8_Tune4CptHard2->Clone("ratioEta2760GeV_Pythia8_Tune4CpTHard2DivMB");
	ratioEta2760GeV_Pythia8_Tune4CpTHard2DivMB->Divide(ratioEta2760GeV_Pythia8_Tune4CpTHard2DivMB, histoEta2760GeV_Pythia8_Tune4CMB,1.,1.,"");
	TH1D* ratioEta2760GeV_Pythia8_Tune4CpTHard5DivpTHard2 = (TH1D*) histoEta2760GeV_Pythia8_Tune4CptHard5->Clone("ratioEta2760GeV_Pythia8_Tune4CpTHard5DivpTHard2");
	ratioEta2760GeV_Pythia8_Tune4CpTHard5DivpTHard2->Divide(ratioEta2760GeV_Pythia8_Tune4CpTHard5DivpTHard2, histoEta2760GeV_Pythia8_Tune4CptHard2,1.,1.,"");
	TH1D* ratioEta2760GeV_Pythia8_Tune4CpTHard10DivpTHard5 = (TH1D*) histoEta2760GeV_Pythia8_Tune4CptHard10->Clone("ratioEta2760GeV_Pythia8_Tune4CpTHard10DivpTHard5");
	ratioEta2760GeV_Pythia8_Tune4CpTHard10DivpTHard5->Divide(ratioEta2760GeV_Pythia8_Tune4CpTHard10DivpTHard5, histoEta2760GeV_Pythia8_Tune4CptHard5,1.,1.,"");    
	TH1D* ratioEtaUB2760GeV_Pythia8_Tune4CpTHard2DivMB = (TH1D*) histoEtaUB2760GeV_Pythia8_Tune4CptHard2->Clone("ratioEtaUB2760GeV_Pythia8_Tune4CpTHard2DivMB");
	ratioEtaUB2760GeV_Pythia8_Tune4CpTHard2DivMB->Divide(ratioEtaUB2760GeV_Pythia8_Tune4CpTHard2DivMB, histoEtaUB2760GeV_Pythia8_Tune4CMB,1.,1.,"");
	TH1D* ratioEtaUB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2 = (TH1D*) histoEtaUB2760GeV_Pythia8_Tune4CptHard5->Clone("ratioEtaUB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2");
	ratioEtaUB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2->Divide(ratioEtaUB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2, histoEtaUB2760GeV_Pythia8_Tune4CptHard2,1.,1.,"");
	TH1D* ratioEtaUB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5 = (TH1D*) histoEtaUB2760GeV_Pythia8_Tune4CptHard10->Clone("ratioEtaUB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5");
	ratioEtaUB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5->Divide(ratioEtaUB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5, histoEtaUB2760GeV_Pythia8_Tune4CptHard5,1.,1.,"");    
	DrawGammaSetMarker(ratioEta2760GeV_Pythia8_Tune4CpTHard2DivMB, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(ratioEta2760GeV_Pythia8_Tune4CpTHard5DivpTHard2, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(ratioEta2760GeV_Pythia8_Tune4CpTHard10DivpTHard5, 20, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(ratioEtaUB2760GeV_Pythia8_Tune4CpTHard2DivMB, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(ratioEtaUB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(ratioEtaUB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5, 20, 1, colorptHard10, colorptHard10);
	
	DrawGammaSetMarker(histoPi02760GeV_Pythia8_Tune4CMB, 20, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoPi0FromK0s2760GeV_Pythia8_Tune4CMB, 24, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoPi0UB2760GeV_Pythia8_Tune4CMB, 20, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoPi0FromK0sUB2760GeV_Pythia8_Tune4CMB, 24, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoPi0FromEta2760GeV_Pythia8_Tune4CMB, 25, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoPi0FromEtaUB2760GeV_Pythia8_Tune4CMB, 25, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoPi0FromLambda2760GeV_Pythia8_Tune4CMB, 29, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoPi0FromLambdaUB2760GeV_Pythia8_Tune4CMB, 29, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoEta2760GeV_Pythia8_Tune4CMB, 20, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoEtaUB2760GeV_Pythia8_Tune4CMB, 20, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoPiPl2760GeV_Pythia8_Tune4CMB, 20, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoPiPlUB2760GeV_Pythia8_Tune4CMB, 20, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoPiMi2760GeV_Pythia8_Tune4CMB, 20, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoPiMiUB2760GeV_Pythia8_Tune4CMB, 20, 1, colorMB, colorMB);

	DrawGammaSetMarker(histoPi02760GeV_Pythia8_Tune4CptHard2, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoPi0UB2760GeV_Pythia8_Tune4CptHard2, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard2, 24, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoPi0FromK0sUB2760GeV_Pythia8_Tune4CptHard2, 24, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoPi0FromEta2760GeV_Pythia8_Tune4CptHard2, 25, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoPi0FromEtaUB2760GeV_Pythia8_Tune4CptHard2, 25, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard2, 29, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoPi0FromLambdaUB2760GeV_Pythia8_Tune4CptHard2, 29, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoEta2760GeV_Pythia8_Tune4CptHard2, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoEtaUB2760GeV_Pythia8_Tune4CptHard2, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoPiPl2760GeV_Pythia8_Tune4CptHard2, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoPiPlUB2760GeV_Pythia8_Tune4CptHard2, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoPiMi2760GeV_Pythia8_Tune4CptHard2, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoPiMiUB2760GeV_Pythia8_Tune4CptHard2, 20, 1, colorptHard2, colorptHard2);
	
	DrawGammaSetMarker(histoPi02760GeV_Pythia8_Tune4CptHard5, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard5, 24, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoPi0FromK0sUB2760GeV_Pythia8_Tune4CptHard5, 24, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoPi0UB2760GeV_Pythia8_Tune4CptHard5, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoPi0FromEta2760GeV_Pythia8_Tune4CptHard5, 25, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoPi0FromEtaUB2760GeV_Pythia8_Tune4CptHard5, 25, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard5, 29, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoPi0FromLambdaUB2760GeV_Pythia8_Tune4CptHard5, 29, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoEta2760GeV_Pythia8_Tune4CptHard5, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoEtaUB2760GeV_Pythia8_Tune4CptHard5, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoPiPl2760GeV_Pythia8_Tune4CptHard5, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoPiPlUB2760GeV_Pythia8_Tune4CptHard5, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoPiMi2760GeV_Pythia8_Tune4CptHard5, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoPiMiUB2760GeV_Pythia8_Tune4CptHard5, 20, 1, colorptHard5, colorptHard5);

	DrawGammaSetMarker(histoPi0FromK0sUB2760GeV_Pythia8_Tune4CptHard10, 24, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard10, 24, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoPi02760GeV_Pythia8_Tune4CptHard10, 20, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoPi0UB2760GeV_Pythia8_Tune4CptHard10, 20, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoPi0FromEta2760GeV_Pythia8_Tune4CptHard10, 25, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoPi0FromEtaUB2760GeV_Pythia8_Tune4CptHard10, 25, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard10, 29, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoPi0FromLambdaUB2760GeV_Pythia8_Tune4CptHard10, 29, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoEta2760GeV_Pythia8_Tune4CptHard10, 20, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoEtaUB2760GeV_Pythia8_Tune4CptHard10, 20, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoPiPl2760GeV_Pythia8_Tune4CptHard10, 20, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoPiPlUB2760GeV_Pythia8_Tune4CptHard10, 20, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoPiMi2760GeV_Pythia8_Tune4CptHard10, 20, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoPiMiUB2760GeV_Pythia8_Tune4CptHard10, 20, 1, colorptHard10, colorptHard10);

	
	//*******************************************************************************************************
	//****************************** Pythia 8.2 Tune Monach *******************************************************
	//*******************************************************************************************************    
	Double_t scalingFactor2760GeVPy8_Monach_MB                = 6.311e02/(2713849945/4e8)*1/4e8*1e9;
	Double_t scalingFactor2760GeVPy8_Monach_ptHard2            = 1.912e02/(508785747/8e7)*1/8e7*1e9;
	Double_t scalingFactor2760GeVPy8_Monach_ptHard5            = 1.181e01/((199972446+136588301)/(3e7+2.05e7))*1/(3e7+2.05e7)*1e9;
	Double_t scalingFactor2760GeVPy8_Monach_ptHard10        = 9.901e-01/(337055749/5e7)*1/5e7*1e9;
	
//     Double_t xSection2760GeVPy8_MonachMB                    = 1.09e3;
//     Double_t xSection2760GeVPy8_MonachMBErr                    = 3.077e-2;
//     Double_t nEvt2760GeVPy8_MonachMB                        = 400e6;
//     Double_t xSection2760GeVPy8_MonachptHard5                = 1.522e1;
//     Double_t xSection2760GeVPy8_MonachptHard5Err            = 8.484e-4;
//     Double_t nEvt2760GeVPy8_MonachptHard5                    = 100e6;
//     Double_t xSection2760GeVPy8_MonachptHard10                = 1.194e0;
//     Double_t xSection2760GeVPy8_MonachptHard10Err            = 6.968e-5;
//     Double_t nEvt2760GeVPy8_MonachptHard10                    = 86e6;
	
	TFile*     file2760GeV_Pythia8_TuneMonachMB                 = TFile::Open("ExternalInput/Theory/Pythia/hist_Pythia8_2_Monash_MB_Backup400.0Mio.root");
	TH1F*     histoPi02760GeV_Pythia8_TuneMonachMB             = (TH1F*)file2760GeV_Pythia8_TuneMonachMB->Get("ptPi0W0KClone");
	TH1F*     histoPi0FromEta2760GeV_Pythia8_TuneMonachMB     = (TH1F*)file2760GeV_Pythia8_TuneMonachMB->Get("ptPi0FromEtaClone");
	TH1F*     histoPi0FromLambda2760GeV_Pythia8_TuneMonachMB     = (TH1F*)file2760GeV_Pythia8_TuneMonachMB->Get("ptPi0FromLambdaClone");
	TH1F*     histoPi0FromK0s2760GeV_Pythia8_TuneMonachMB     = (TH1F*)file2760GeV_Pythia8_TuneMonachMB->Get("ptPi0Clone");
	histoPi0FromK0s2760GeV_Pythia8_TuneMonachMB->Sumw2();
	histoPi0FromK0s2760GeV_Pythia8_TuneMonachMB->Add(histoPi02760GeV_Pythia8_TuneMonachMB,-1);
	TH1F*     histoEta2760GeV_Pythia8_TuneMonachMB             = (TH1F*)file2760GeV_Pythia8_TuneMonachMB->Get("ptEtaClone");
	TH1F*     histoPiPl2760GeV_Pythia8_TuneMonachMB             = (TH1F*)file2760GeV_Pythia8_TuneMonachMB->Get("ptPiPlW0KClone");
	TH1F*     histoPiMi2760GeV_Pythia8_TuneMonachMB             = (TH1F*)file2760GeV_Pythia8_TuneMonachMB->Get("ptPiMiW0KClone");
	TH1F*     histoPi0UB2760GeV_Pythia8_TuneMonachMB             = (TH1F*)file2760GeV_Pythia8_TuneMonachMB->Get("ptPi0W0KUBClone");
	TH1F*     histoPi0FromEtaUB2760GeV_Pythia8_TuneMonachMB     = (TH1F*)file2760GeV_Pythia8_TuneMonachMB->Get("ptPi0FromEtaUBClone");
	TH1F*     histoPi0FromLambdaUB2760GeV_Pythia8_TuneMonachMB     = (TH1F*)file2760GeV_Pythia8_TuneMonachMB->Get("ptPi0FromLambdaUBClone");
	TH1F*     histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachMB     = (TH1F*)file2760GeV_Pythia8_TuneMonachMB->Get("ptPi0UBClone");
	histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachMB->Sumw2();
	histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachMB->Add(histoPi0UB2760GeV_Pythia8_TuneMonachMB,-1);
	TH1F*     histoEtaUB2760GeV_Pythia8_TuneMonachMB             = (TH1F*)file2760GeV_Pythia8_TuneMonachMB->Get("ptEtaUBClone");
	TH1F*     histoPiPlUB2760GeV_Pythia8_TuneMonachMB         = (TH1F*)file2760GeV_Pythia8_TuneMonachMB->Get("ptPiPlW0KUBClone");
	TH1F*     histoPiMiUB2760GeV_Pythia8_TuneMonachMB         = (TH1F*)file2760GeV_Pythia8_TuneMonachMB->Get("ptPiMiW0KUBClone");
	
	histoPi02760GeV_Pythia8_TuneMonachMB->Scale(scalingFactor2760GeVPy8_Monach_MB);
	histoPi0FromLambda2760GeV_Pythia8_TuneMonachMB->Scale(scalingFactor2760GeVPy8_Monach_MB);
	histoPi0FromEta2760GeV_Pythia8_TuneMonachMB->Scale(scalingFactor2760GeVPy8_Monach_MB);
	histoPi0FromK0s2760GeV_Pythia8_TuneMonachMB->Scale(scalingFactor2760GeVPy8_Monach_MB);
	histoEta2760GeV_Pythia8_TuneMonachMB->Scale(scalingFactor2760GeVPy8_Monach_MB);
	histoPiPl2760GeV_Pythia8_TuneMonachMB->Scale(scalingFactor2760GeVPy8_Monach_MB);
	histoPiMi2760GeV_Pythia8_TuneMonachMB->Scale(scalingFactor2760GeVPy8_Monach_MB);
	histoPi0UB2760GeV_Pythia8_TuneMonachMB->Scale(scalingFactor2760GeVPy8_Monach_MB);
	histoPi0FromLambdaUB2760GeV_Pythia8_TuneMonachMB->Scale(scalingFactor2760GeVPy8_Monach_MB);
	histoPi0FromEtaUB2760GeV_Pythia8_TuneMonachMB->Scale(scalingFactor2760GeVPy8_Monach_MB);
	histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachMB->Scale(scalingFactor2760GeVPy8_Monach_MB);
	histoEtaUB2760GeV_Pythia8_TuneMonachMB->Scale(scalingFactor2760GeVPy8_Monach_MB);
	histoPiPlUB2760GeV_Pythia8_TuneMonachMB->Scale(scalingFactor2760GeVPy8_Monach_MB);
	histoPiMiUB2760GeV_Pythia8_TuneMonachMB->Scale(scalingFactor2760GeVPy8_Monach_MB);

	TFile*     file2760GeV_Pythia8_TuneMonachptHard2                 = TFile::Open("ExternalInput/Theory/Pythia/hist_Pythia8_2_Monash_pTHard2_Backup80.0Mio.root");
	TH1F*     histoPi02760GeV_Pythia8_TuneMonachptHard2             = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard2->Get("ptPi0W0KClone");
	TH1F*     histoPi0FromEta2760GeV_Pythia8_TuneMonachptHard2     = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard2->Get("ptPi0FromEtaClone");
	TH1F*     histoPi0FromLambda2760GeV_Pythia8_TuneMonachptHard2 = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard2->Get("ptPi0FromLambdaClone");
	TH1F*     histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard2     = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard2->Get("ptPi0Clone");
	histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard2->Sumw2();
	histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard2->Add(histoPi02760GeV_Pythia8_TuneMonachptHard2,-1);
	TH1F*     histoEta2760GeV_Pythia8_TuneMonachptHard2             = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard2->Get("ptEtaClone");
	TH1F*     histoPiPl2760GeV_Pythia8_TuneMonachptHard2             = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard2->Get("ptPiPlW0KClone");
	TH1F*     histoPiMi2760GeV_Pythia8_TuneMonachptHard2             = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard2->Get("ptPiMiW0KClone");
	TH1F*     histoPi0UB2760GeV_Pythia8_TuneMonachptHard2         = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard2->Get("ptPi0W0KUBClone");
	TH1F*     histoPi0FromEtaUB2760GeV_Pythia8_TuneMonachptHard2     = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard2->Get("ptPi0FromEtaUBClone");
	TH1F*     histoPi0FromLambdaUB2760GeV_Pythia8_TuneMonachptHard2     = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard2->Get("ptPi0FromLambdaUBClone");
	TH1F*     histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachptHard2     = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard2->Get("ptPi0UBClone");
	histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachptHard2->Sumw2();
	histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachptHard2->Add(histoPi0UB2760GeV_Pythia8_TuneMonachptHard2,-1);
	TH1F*     histoEtaUB2760GeV_Pythia8_TuneMonachptHard2         = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard2->Get("ptEtaUBClone");
	TH1F*     histoPiPlUB2760GeV_Pythia8_TuneMonachptHard2         = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard2->Get("ptPiPlW0KUBClone");
	TH1F*     histoPiMiUB2760GeV_Pythia8_TuneMonachptHard2         = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard2->Get("ptPiMiW0KUBClone");
	
	histoPi02760GeV_Pythia8_TuneMonachptHard2->Scale(scalingFactor2760GeVPy8_Monach_ptHard2);
	histoPi0FromLambda2760GeV_Pythia8_TuneMonachptHard2->Scale(scalingFactor2760GeVPy8_Monach_ptHard2);
	histoPi0FromEta2760GeV_Pythia8_TuneMonachptHard2->Scale(scalingFactor2760GeVPy8_Monach_ptHard2);
	histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard2->Scale(scalingFactor2760GeVPy8_Monach_ptHard2);
	histoEta2760GeV_Pythia8_TuneMonachptHard2->Scale(scalingFactor2760GeVPy8_Monach_ptHard2);
	histoPiPl2760GeV_Pythia8_TuneMonachptHard2->Scale(scalingFactor2760GeVPy8_Monach_ptHard2);
	histoPiMi2760GeV_Pythia8_TuneMonachptHard2->Scale(scalingFactor2760GeVPy8_Monach_ptHard2);
	histoPi0UB2760GeV_Pythia8_TuneMonachptHard2->Scale(scalingFactor2760GeVPy8_Monach_ptHard2);
	histoPi0FromLambdaUB2760GeV_Pythia8_TuneMonachptHard2->Scale(scalingFactor2760GeVPy8_Monach_ptHard2);
	histoPi0FromEtaUB2760GeV_Pythia8_TuneMonachptHard2->Scale(scalingFactor2760GeVPy8_Monach_ptHard2);
	histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachptHard2->Scale(scalingFactor2760GeVPy8_Monach_ptHard2);
	histoEtaUB2760GeV_Pythia8_TuneMonachptHard2->Scale(scalingFactor2760GeVPy8_Monach_ptHard2);
	histoPiPlUB2760GeV_Pythia8_TuneMonachptHard2->Scale(scalingFactor2760GeVPy8_Monach_ptHard2);
	histoPiMiUB2760GeV_Pythia8_TuneMonachptHard2->Scale(scalingFactor2760GeVPy8_Monach_ptHard2);
	
	
	TFile*     file2760GeV_Pythia8_TuneMonachptHard5                 = TFile::Open("ExternalInput/Theory/Pythia/hist_Pythia8_2_Monash_pTHard5_50.5Mio.root");
	TH1F*     histoPi02760GeV_Pythia8_TuneMonachptHard5             = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard5->Get("ptPi0W0KClone");
	TH1F*     histoPi0FromEta2760GeV_Pythia8_TuneMonachptHard5     = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard5->Get("ptPi0FromEtaClone");
	TH1F*     histoPi0FromLambda2760GeV_Pythia8_TuneMonachptHard5 = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard5->Get("ptPi0FromLambdaClone");
	TH1F*     histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard5     = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard5->Get("ptPi0Clone");
	histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard5->Sumw2();
	histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard5->Add(histoPi02760GeV_Pythia8_TuneMonachptHard5,-1);
	TH1F*     histoEta2760GeV_Pythia8_TuneMonachptHard5             = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard5->Get("ptEtaClone");
	TH1F*     histoPiPl2760GeV_Pythia8_TuneMonachptHard5             = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard5->Get("ptPiPlW0KClone");
	TH1F*     histoPiMi2760GeV_Pythia8_TuneMonachptHard5             = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard5->Get("ptPiMiW0KClone");
	TH1F*     histoPi0UB2760GeV_Pythia8_TuneMonachptHard5         = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard5->Get("ptPi0W0KUBClone");
	TH1F*     histoPi0FromEtaUB2760GeV_Pythia8_TuneMonachptHard5     = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard5->Get("ptPi0FromEtaUBClone");
	TH1F*     histoPi0FromLambdaUB2760GeV_Pythia8_TuneMonachptHard5     = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard5->Get("ptPi0FromLambdaUBClone");
	TH1F*     histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachptHard5     = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard5->Get("ptPi0UBClone");
	histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachptHard5->Sumw2();
	histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachptHard5->Add(histoPi0UB2760GeV_Pythia8_TuneMonachptHard5,-1);
	TH1F*     histoEtaUB2760GeV_Pythia8_TuneMonachptHard5         = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard5->Get("ptEtaUBClone");
	TH1F*     histoPiPlUB2760GeV_Pythia8_TuneMonachptHard5         = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard5->Get("ptPiPlW0KUBClone");
	TH1F*     histoPiMiUB2760GeV_Pythia8_TuneMonachptHard5         = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard5->Get("ptPiMiW0KUBClone");
	
	histoPi02760GeV_Pythia8_TuneMonachptHard5->Scale(scalingFactor2760GeVPy8_Monach_ptHard5);
	histoPi0FromLambda2760GeV_Pythia8_TuneMonachptHard5->Scale(scalingFactor2760GeVPy8_Monach_ptHard5);
	histoPi0FromEta2760GeV_Pythia8_TuneMonachptHard5->Scale(scalingFactor2760GeVPy8_Monach_ptHard5);
	histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard5->Scale(scalingFactor2760GeVPy8_Monach_ptHard5);
	histoEta2760GeV_Pythia8_TuneMonachptHard5->Scale(scalingFactor2760GeVPy8_Monach_ptHard5);
	histoPiPl2760GeV_Pythia8_TuneMonachptHard5->Scale(scalingFactor2760GeVPy8_Monach_ptHard5);
	histoPiMi2760GeV_Pythia8_TuneMonachptHard5->Scale(scalingFactor2760GeVPy8_Monach_ptHard5);
	histoPi0UB2760GeV_Pythia8_TuneMonachptHard5->Scale(scalingFactor2760GeVPy8_Monach_ptHard5);
	histoPi0FromLambdaUB2760GeV_Pythia8_TuneMonachptHard5->Scale(scalingFactor2760GeVPy8_Monach_ptHard5);
	histoPi0FromEtaUB2760GeV_Pythia8_TuneMonachptHard5->Scale(scalingFactor2760GeVPy8_Monach_ptHard5);
	histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachptHard5->Scale(scalingFactor2760GeVPy8_Monach_ptHard5);
	histoEtaUB2760GeV_Pythia8_TuneMonachptHard5->Scale(scalingFactor2760GeVPy8_Monach_ptHard5);
	histoPiPlUB2760GeV_Pythia8_TuneMonachptHard5->Scale(scalingFactor2760GeVPy8_Monach_ptHard5);
	histoPiMiUB2760GeV_Pythia8_TuneMonachptHard5->Scale(scalingFactor2760GeVPy8_Monach_ptHard5);

	TFile*     file2760GeV_Pythia8_TuneMonachptHard10                 = TFile::Open("ExternalInput/Theory/Pythia/hist_Pythia8_2_Monash_pTHard10_50Mio.root");
	TH1F*     histoPi02760GeV_Pythia8_TuneMonachptHard10             = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard10->Get("ptPi0W0K");
	TH1F*     histoPi0FromEta2760GeV_Pythia8_TuneMonachptHard10     = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard10->Get("ptPi0FromEta");
	TH1F*     histoPi0FromLambda2760GeV_Pythia8_TuneMonachptHard10     = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard10->Get("ptPi0FromLambda");
	TH1F*     histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard10     = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard10->Get("ptPi0");
	histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard10->Sumw2();
	histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard10->Add(histoPi02760GeV_Pythia8_TuneMonachptHard10,-1);
	TH1F*     histoEta2760GeV_Pythia8_TuneMonachptHard10             = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard10->Get("ptEta");
	TH1F*     histoPiPl2760GeV_Pythia8_TuneMonachptHard10             = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard10->Get("ptPiPlW0K");
	TH1F*     histoPiMi2760GeV_Pythia8_TuneMonachptHard10             = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard10->Get("ptPiMiW0K");
	TH1F*     histoPi0UB2760GeV_Pythia8_TuneMonachptHard10             = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard10->Get("ptPi0W0KUB");
	TH1F*     histoPi0FromEtaUB2760GeV_Pythia8_TuneMonachptHard10     = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard10->Get("ptPi0FromEtaUB");
	TH1F*     histoPi0FromLambdaUB2760GeV_Pythia8_TuneMonachptHard10     = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard10->Get("ptPi0FromLambdaUB");
	TH1F*     histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachptHard10     = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard10->Get("ptPi0UB");
	histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachptHard10->Sumw2();
	histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachptHard10->Add(histoPi0UB2760GeV_Pythia8_TuneMonachptHard10,-1);
	TH1F*     histoEtaUB2760GeV_Pythia8_TuneMonachptHard10             = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard10->Get("ptEtaUB");
	TH1F*     histoPiPlUB2760GeV_Pythia8_TuneMonachptHard10             = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard10->Get("ptPiPlW0KUB");
	TH1F*     histoPiMiUB2760GeV_Pythia8_TuneMonachptHard10             = (TH1F*)file2760GeV_Pythia8_TuneMonachptHard10->Get("ptPiMiW0KUB");
	
	histoPi02760GeV_Pythia8_TuneMonachptHard10->Scale(scalingFactor2760GeVPy8_Monach_ptHard10);
	histoPi0FromLambda2760GeV_Pythia8_TuneMonachptHard10->Scale(scalingFactor2760GeVPy8_Monach_ptHard10);
	histoPi0FromEta2760GeV_Pythia8_TuneMonachptHard10->Scale(scalingFactor2760GeVPy8_Monach_ptHard10);
	histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard10->Scale(scalingFactor2760GeVPy8_Monach_ptHard10);
	histoEta2760GeV_Pythia8_TuneMonachptHard10->Scale(scalingFactor2760GeVPy8_Monach_ptHard10);
	histoPiPl2760GeV_Pythia8_TuneMonachptHard10->Scale(scalingFactor2760GeVPy8_Monach_ptHard10);
	histoPiMi2760GeV_Pythia8_TuneMonachptHard10->Scale(scalingFactor2760GeVPy8_Monach_ptHard10);
	histoPi0UB2760GeV_Pythia8_TuneMonachptHard10->Scale(scalingFactor2760GeVPy8_Monach_ptHard10);
	histoPi0FromLambdaUB2760GeV_Pythia8_TuneMonachptHard10->Scale(scalingFactor2760GeVPy8_Monach_ptHard10);
	histoPi0FromEtaUB2760GeV_Pythia8_TuneMonachptHard10->Scale(scalingFactor2760GeVPy8_Monach_ptHard10);
	histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachptHard10->Scale(scalingFactor2760GeVPy8_Monach_ptHard10);
	histoEtaUB2760GeV_Pythia8_TuneMonachptHard10->Scale(scalingFactor2760GeVPy8_Monach_ptHard10);
	histoPiPlUB2760GeV_Pythia8_TuneMonachptHard10->Scale(scalingFactor2760GeVPy8_Monach_ptHard10);
	histoPiMiUB2760GeV_Pythia8_TuneMonachptHard10->Scale(scalingFactor2760GeVPy8_Monach_ptHard10);

	TH1D* ratioPi02760GeV_Pythia8_TuneMonachpTHard2DivMB = (TH1D*) histoPi02760GeV_Pythia8_TuneMonachptHard2->Clone("ratioPi02760GeV_Pythia8_TuneMonachpTHard2DivMB");
	ratioPi02760GeV_Pythia8_TuneMonachpTHard2DivMB->Divide(ratioPi02760GeV_Pythia8_TuneMonachpTHard2DivMB, histoPi02760GeV_Pythia8_TuneMonachMB,1.,1.,"");
	TH1D* ratioPi02760GeV_Pythia8_TuneMonachpTHard5DivpTHard2 = (TH1D*) histoPi02760GeV_Pythia8_TuneMonachptHard5->Clone("ratioPi02760GeV_Pythia8_TuneMonachpTHard5DivpTHard2");
	ratioPi02760GeV_Pythia8_TuneMonachpTHard5DivpTHard2->Divide(ratioPi02760GeV_Pythia8_TuneMonachpTHard5DivpTHard2, histoPi02760GeV_Pythia8_TuneMonachptHard2,1.,1.,"");
	TH1D* ratioPi02760GeV_Pythia8_TuneMonachpTHard10DivpTHard5 = (TH1D*) histoPi02760GeV_Pythia8_TuneMonachptHard10->Clone("ratioPi02760GeV_Pythia8_TuneMonachpTHard10DivpTHard5");
	ratioPi02760GeV_Pythia8_TuneMonachpTHard10DivpTHard5->Divide(ratioPi02760GeV_Pythia8_TuneMonachpTHard10DivpTHard5, histoPi02760GeV_Pythia8_TuneMonachptHard5,1.,1.,"");    
	TH1D* ratioPi0UB2760GeV_Pythia8_TuneMonachpTHard2DivMB = (TH1D*) histoPi0UB2760GeV_Pythia8_TuneMonachptHard2->Clone("ratioPi0UB2760GeV_Pythia8_TuneMonachpTHard2DivMB");
	ratioPi0UB2760GeV_Pythia8_TuneMonachpTHard2DivMB->Divide(ratioPi0UB2760GeV_Pythia8_TuneMonachpTHard2DivMB, histoPi0UB2760GeV_Pythia8_TuneMonachMB,1.,1.,"");
	TH1D* ratioPi0UB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2 = (TH1D*) histoPi0UB2760GeV_Pythia8_TuneMonachptHard5->Clone("ratioPi0UB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2");
	ratioPi0UB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2->Divide(ratioPi0UB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2, histoPi0UB2760GeV_Pythia8_TuneMonachptHard2,1.,1.,"");
	TH1D* ratioPi0UB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5 = (TH1D*) histoPi0UB2760GeV_Pythia8_TuneMonachptHard10->Clone("ratioPi0UB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5");
	ratioPi0UB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5->Divide(ratioPi0UB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5, histoPi0UB2760GeV_Pythia8_TuneMonachptHard5,1.,1.,"");    
	DrawGammaSetMarker(ratioPi02760GeV_Pythia8_TuneMonachpTHard2DivMB, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(ratioPi02760GeV_Pythia8_TuneMonachpTHard5DivpTHard2, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(ratioPi02760GeV_Pythia8_TuneMonachpTHard10DivpTHard5, 20, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(ratioPi0UB2760GeV_Pythia8_TuneMonachpTHard2DivMB, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(ratioPi0UB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(ratioPi0UB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5, 20, 1, colorptHard10, colorptHard10);

	TH1D* ratioPiPl2760GeV_Pythia8_TuneMonachpTHard2DivMB = (TH1D*) histoPiPl2760GeV_Pythia8_TuneMonachptHard2->Clone("ratioPiPl2760GeV_Pythia8_TuneMonachpTHard2DivMB");
	ratioPiPl2760GeV_Pythia8_TuneMonachpTHard2DivMB->Divide(ratioPiPl2760GeV_Pythia8_TuneMonachpTHard2DivMB, histoPiPl2760GeV_Pythia8_TuneMonachMB,1.,1.,"");
	TH1D* ratioPiPl2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2 = (TH1D*) histoPiPl2760GeV_Pythia8_TuneMonachptHard5->Clone("ratioPiPl2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2");
	ratioPiPl2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2->Divide(ratioPiPl2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2, histoPiPl2760GeV_Pythia8_TuneMonachptHard2,1.,1.,"");
	TH1D* ratioPiPl2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5 = (TH1D*) histoPiPl2760GeV_Pythia8_TuneMonachptHard10->Clone("ratioPiPl2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5");
	ratioPiPl2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5->Divide(ratioPiPl2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5, histoPiPl2760GeV_Pythia8_TuneMonachptHard5,1.,1.,"");    
	TH1D* ratioPiPlUB2760GeV_Pythia8_TuneMonachpTHard2DivMB = (TH1D*) histoPiPlUB2760GeV_Pythia8_TuneMonachptHard2->Clone("ratioPiPlUB2760GeV_Pythia8_TuneMonachpTHard2DivMB");
	ratioPiPlUB2760GeV_Pythia8_TuneMonachpTHard2DivMB->Divide(ratioPiPlUB2760GeV_Pythia8_TuneMonachpTHard2DivMB, histoPiPlUB2760GeV_Pythia8_TuneMonachMB,1.,1.,"");
	TH1D* ratioPiPlUB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2 = (TH1D*) histoPiPlUB2760GeV_Pythia8_TuneMonachptHard5->Clone("ratioPiPlUB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2");
	ratioPiPlUB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2->Divide(ratioPiPlUB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2, histoPiPlUB2760GeV_Pythia8_TuneMonachptHard2,1.,1.,"");
	TH1D* ratioPiPlUB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5 = (TH1D*) histoPiPlUB2760GeV_Pythia8_TuneMonachptHard10->Clone("ratioPiPlUB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5");
	ratioPiPlUB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5->Divide(ratioPiPlUB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5, histoPiPlUB2760GeV_Pythia8_TuneMonachptHard5,1.,1.,"");    
	DrawGammaSetMarker(ratioPiPl2760GeV_Pythia8_TuneMonachpTHard2DivMB, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(ratioPiPl2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(ratioPiPl2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5, 20, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(ratioPiPlUB2760GeV_Pythia8_TuneMonachpTHard2DivMB, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(ratioPiPlUB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(ratioPiPlUB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5, 20, 1, colorptHard10, colorptHard10);

	TH1D* ratioPiMi2760GeV_Pythia8_TuneMonachpTHard2DivMB = (TH1D*) histoPiMi2760GeV_Pythia8_TuneMonachptHard2->Clone("ratioPiMi2760GeV_Pythia8_TuneMonachpTHard2DivMB");
	ratioPiMi2760GeV_Pythia8_TuneMonachpTHard2DivMB->Divide(ratioPiMi2760GeV_Pythia8_TuneMonachpTHard2DivMB, histoPiMi2760GeV_Pythia8_TuneMonachMB,1.,1.,"");
	TH1D* ratioPiMi2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2 = (TH1D*) histoPiMi2760GeV_Pythia8_TuneMonachptHard5->Clone("ratioPiMi2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2");
	ratioPiMi2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2->Divide(ratioPiMi2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2, histoPiMi2760GeV_Pythia8_TuneMonachptHard2,1.,1.,"");
	TH1D* ratioPiMi2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5 = (TH1D*) histoPiMi2760GeV_Pythia8_TuneMonachptHard10->Clone("ratioPiMi2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5");
	ratioPiMi2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5->Divide(ratioPiMi2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5, histoPiMi2760GeV_Pythia8_TuneMonachptHard5,1.,1.,"");    
	TH1D* ratioPiMiUB2760GeV_Pythia8_TuneMonachpTHard2DivMB = (TH1D*) histoPiMiUB2760GeV_Pythia8_TuneMonachptHard2->Clone("ratioPiMiUB2760GeV_Pythia8_TuneMonachpTHard2DivMB");
	ratioPiMiUB2760GeV_Pythia8_TuneMonachpTHard2DivMB->Divide(ratioPiMiUB2760GeV_Pythia8_TuneMonachpTHard2DivMB, histoPiMiUB2760GeV_Pythia8_TuneMonachMB,1.,1.,"");
	TH1D* ratioPiMiUB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2 = (TH1D*) histoPiMiUB2760GeV_Pythia8_TuneMonachptHard5->Clone("ratioPiMiUB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2");
	ratioPiMiUB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2->Divide(ratioPiMiUB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2, histoPiMiUB2760GeV_Pythia8_TuneMonachptHard2,1.,1.,"");
	TH1D* ratioPiMiUB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5 = (TH1D*) histoPiMiUB2760GeV_Pythia8_TuneMonachptHard10->Clone("ratioPiMiUB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5");
	ratioPiMiUB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5->Divide(ratioPiMiUB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5, histoPiMiUB2760GeV_Pythia8_TuneMonachptHard5,1.,1.,"");    
	DrawGammaSetMarker(ratioPiMi2760GeV_Pythia8_TuneMonachpTHard2DivMB, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(ratioPiMi2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(ratioPiMi2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5, 20, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(ratioPiMiUB2760GeV_Pythia8_TuneMonachpTHard2DivMB, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(ratioPiMiUB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(ratioPiMiUB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5, 20, 1, colorptHard10, colorptHard10);

	
	TH1D* ratioEta2760GeV_Pythia8_TuneMonachpTHard2DivMB = (TH1D*) histoEta2760GeV_Pythia8_TuneMonachptHard2->Clone("ratioEta2760GeV_Pythia8_TuneMonachpTHard2DivMB");
	ratioEta2760GeV_Pythia8_TuneMonachpTHard2DivMB->Divide(ratioEta2760GeV_Pythia8_TuneMonachpTHard2DivMB, histoEta2760GeV_Pythia8_TuneMonachMB,1.,1.,"");
	TH1D* ratioEta2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2 = (TH1D*) histoEta2760GeV_Pythia8_TuneMonachptHard5->Clone("ratioEta2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2");
	ratioEta2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2->Divide(ratioEta2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2, histoEta2760GeV_Pythia8_TuneMonachptHard2,1.,1.,"");
	TH1D* ratioEta2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5 = (TH1D*) histoEta2760GeV_Pythia8_TuneMonachptHard10->Clone("ratioEta2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5");
	ratioEta2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5->Divide(ratioEta2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5, histoEta2760GeV_Pythia8_TuneMonachptHard5,1.,1.,"");    
	TH1D* ratioEtaUB2760GeV_Pythia8_TuneMonachpTHard2DivMB = (TH1D*) histoEtaUB2760GeV_Pythia8_TuneMonachptHard2->Clone("ratioEtaUB2760GeV_Pythia8_TuneMonachpTHard2DivMB");
	ratioEtaUB2760GeV_Pythia8_TuneMonachpTHard2DivMB->Divide(ratioEtaUB2760GeV_Pythia8_TuneMonachpTHard2DivMB, histoEtaUB2760GeV_Pythia8_TuneMonachMB,1.,1.,"");
	TH1D* ratioEtaUB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2 = (TH1D*) histoEtaUB2760GeV_Pythia8_TuneMonachptHard5->Clone("ratioEtaUB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2");
	ratioEtaUB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2->Divide(ratioEtaUB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2, histoEtaUB2760GeV_Pythia8_TuneMonachptHard2,1.,1.,"");
	TH1D* ratioEtaUB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5 = (TH1D*) histoEtaUB2760GeV_Pythia8_TuneMonachptHard10->Clone("ratioEtaUB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5");
	ratioEtaUB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5->Divide(ratioEtaUB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5, histoEtaUB2760GeV_Pythia8_TuneMonachptHard5,1.,1.,"");    
	DrawGammaSetMarker(ratioEta2760GeV_Pythia8_TuneMonachpTHard2DivMB, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(ratioEta2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(ratioEta2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5, 20, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(ratioEtaUB2760GeV_Pythia8_TuneMonachpTHard2DivMB, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(ratioEtaUB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(ratioEtaUB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5, 20, 1, colorptHard10, colorptHard10);
	
	DrawGammaSetMarker(histoPi02760GeV_Pythia8_TuneMonachMB, 20, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoPi0FromK0s2760GeV_Pythia8_TuneMonachMB, 24, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoPi0UB2760GeV_Pythia8_TuneMonachMB, 20, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachMB, 24, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoPi0FromEta2760GeV_Pythia8_TuneMonachMB, 25, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoPi0FromEtaUB2760GeV_Pythia8_TuneMonachMB, 25, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoPi0FromLambda2760GeV_Pythia8_TuneMonachMB, 29, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoPi0FromLambdaUB2760GeV_Pythia8_TuneMonachMB, 29, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoEta2760GeV_Pythia8_TuneMonachMB, 20, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoEtaUB2760GeV_Pythia8_TuneMonachMB, 20, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoPiPl2760GeV_Pythia8_TuneMonachMB, 20, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoPiPlUB2760GeV_Pythia8_TuneMonachMB, 20, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoPiMi2760GeV_Pythia8_TuneMonachMB, 20, 1, colorMB, colorMB);
	DrawGammaSetMarker(histoPiMiUB2760GeV_Pythia8_TuneMonachMB, 20, 1, colorMB, colorMB);

	DrawGammaSetMarker(histoPi02760GeV_Pythia8_TuneMonachptHard2, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoPi0UB2760GeV_Pythia8_TuneMonachptHard2, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard2, 24, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachptHard2, 24, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoPi0FromEta2760GeV_Pythia8_TuneMonachptHard2, 25, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoPi0FromEtaUB2760GeV_Pythia8_TuneMonachptHard2, 25, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoPi0FromLambda2760GeV_Pythia8_TuneMonachptHard2, 29, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoPi0FromLambdaUB2760GeV_Pythia8_TuneMonachptHard2, 29, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoEta2760GeV_Pythia8_TuneMonachptHard2, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoEtaUB2760GeV_Pythia8_TuneMonachptHard2, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoPiPl2760GeV_Pythia8_TuneMonachptHard2, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoPiPlUB2760GeV_Pythia8_TuneMonachptHard2, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoPiMi2760GeV_Pythia8_TuneMonachptHard2, 20, 1, colorptHard2, colorptHard2);
	DrawGammaSetMarker(histoPiMiUB2760GeV_Pythia8_TuneMonachptHard2, 20, 1, colorptHard2, colorptHard2);
	
	DrawGammaSetMarker(histoPi02760GeV_Pythia8_TuneMonachptHard5, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard5, 24, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachptHard5, 24, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoPi0UB2760GeV_Pythia8_TuneMonachptHard5, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoPi0FromEta2760GeV_Pythia8_TuneMonachptHard5, 25, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoPi0FromEtaUB2760GeV_Pythia8_TuneMonachptHard5, 25, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoPi0FromLambda2760GeV_Pythia8_TuneMonachptHard5, 29, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoPi0FromLambdaUB2760GeV_Pythia8_TuneMonachptHard5, 29, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoEta2760GeV_Pythia8_TuneMonachptHard5, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoEtaUB2760GeV_Pythia8_TuneMonachptHard5, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoPiPl2760GeV_Pythia8_TuneMonachptHard5, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoPiPlUB2760GeV_Pythia8_TuneMonachptHard5, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoPiMi2760GeV_Pythia8_TuneMonachptHard5, 20, 1, colorptHard5, colorptHard5);
	DrawGammaSetMarker(histoPiMiUB2760GeV_Pythia8_TuneMonachptHard5, 20, 1, colorptHard5, colorptHard5);

	DrawGammaSetMarker(histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachptHard10, 24, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard10, 24, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoPi02760GeV_Pythia8_TuneMonachptHard10, 20, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoPi0UB2760GeV_Pythia8_TuneMonachptHard10, 20, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoPi0FromEta2760GeV_Pythia8_TuneMonachptHard10, 25, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoPi0FromEtaUB2760GeV_Pythia8_TuneMonachptHard10, 25, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoPi0FromLambda2760GeV_Pythia8_TuneMonachptHard10, 29, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoPi0FromLambdaUB2760GeV_Pythia8_TuneMonachptHard10, 29, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoEta2760GeV_Pythia8_TuneMonachptHard10, 20, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoEtaUB2760GeV_Pythia8_TuneMonachptHard10, 20, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoPiPl2760GeV_Pythia8_TuneMonachptHard10, 20, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoPiPlUB2760GeV_Pythia8_TuneMonachptHard10, 20, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoPiMi2760GeV_Pythia8_TuneMonachptHard10, 20, 1, colorptHard10, colorptHard10);
	DrawGammaSetMarker(histoPiMiUB2760GeV_Pythia8_TuneMonachptHard10, 20, 1, colorptHard10, colorptHard10);
	
	//************************************************************************************************
	//************************ Plotting xSection *****************************************************
	//************************************************************************************************
	
	TString collisionSystem2760GeV             = "pp, #sqrt{#it{s}} = 2.76 TeV";   

	TCanvas* canvasXSectionPi0 = new TCanvas("canvasXSectionPi0","",200,10,1350,1350*1.15);  // gives the page size
	DrawGammaCanvasSettings( canvasXSectionPi0, 0.14, 0.02, 0.02, 0.09);
	canvasXSectionPi0->SetLogx();
	canvasXSectionPi0->SetLogy();
	
	TH2F * histo2DXSectionPi0;
	histo2DXSectionPi0 = new TH2F("histo2DXSectionPi0","histo2DXSectionPi0",11000,0.23,70.,1000,2e-2,10e11);
	SetStyleHistoTH2ForGraphs(histo2DXSectionPi0, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",0.035,0.04, 0.035,0.04, 1.,1.45);
	histo2DXSectionPi0->GetXaxis()->SetMoreLogLabels();
	histo2DXSectionPi0->GetXaxis()->SetLabelOffset(-0.01);
	histo2DXSectionPi0->Draw("copy");

		histoPi02760GeV_Pythia8_Tune4CMB->Draw("pE1same");
		histoPi02760GeV_Pythia8_Tune4CptHard2->Draw("pE1same");
		histoPi02760GeV_Pythia8_Tune4CptHard5->Draw("pE1same");
		histoPi02760GeV_Pythia8_Tune4CptHard10->Draw("pE1same");

		TLatex *labelEnergyXSectionPi0 = new TLatex(0.64,0.92,collisionSystem2760GeV.Data());
		SetStyleTLatex( labelEnergyXSectionPi0, 0.035,4);
		labelEnergyXSectionPi0->Draw();
		TLatex *labelDetSysXSectionPi0 = new TLatex(0.64,0.88,"#pi^{0} #rightarrow #gamma#gamma");
		SetStyleTLatex( labelDetSysXSectionPi0, 0.035,4);
		labelDetSysXSectionPi0->Draw();
		TLatex *labelPythiaTuneA = new TLatex(0.64,0.84,"Pythia 8.1, Tune 4C");
		SetStyleTLatex( labelPythiaTuneA, 0.035,4);
		labelPythiaTuneA->Draw();

		TLegend* legendXSectionPi0 = new TLegend(0.62,0.66,0.9,0.82);
		legendXSectionPi0->SetFillColor(0);
		legendXSectionPi0->SetLineColor(0);
		legendXSectionPi0->SetTextFont(42);
		legendXSectionPi0->SetTextSize(0.035);
		legendXSectionPi0->AddEntry(histoPi02760GeV_Pythia8_Tune4CMB,"MB","p");
		legendXSectionPi0->AddEntry(histoPi02760GeV_Pythia8_Tune4CptHard2,"p_{T,hard} > 2 GeV/c","p");
		legendXSectionPi0->AddEntry(histoPi02760GeV_Pythia8_Tune4CptHard5,"p_{T,hard} > 5 GeV/c","p");
		legendXSectionPi0->AddEntry(histoPi02760GeV_Pythia8_Tune4CptHard10,"p_{T,hard} > 10 GeV/c","p");
		legendXSectionPi0->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/Tune4C_Pi0XSectionDiffptHardBins.eps"));

	histo2DXSectionPi0->Draw("copy");
		histoPi02760GeV_Pythia8_TuneMonachMB->Draw("pE1same");
		histoPi02760GeV_Pythia8_TuneMonachptHard2->Draw("pE1same");
		histoPi02760GeV_Pythia8_TuneMonachptHard5->Draw("pE1same");
		histoPi02760GeV_Pythia8_TuneMonachptHard10->Draw("pE1same");

		labelEnergyXSectionPi0->Draw();
		labelDetSysXSectionPi0->Draw();
		TLatex *labelPythiaTuneB = new TLatex(0.64,0.84,"Pythia 8.2, Monach");
		SetStyleTLatex( labelPythiaTuneB, 0.035,4);
		labelPythiaTuneB->Draw();

		legendXSectionPi0->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/TuneMonach_Pi0XSectionDiffptHardBins.eps"));

	histo2DXSectionPi0->Draw("copy");
		histoPi0UB2760GeV_Pythia8_Tune4CMB->Draw("pE1same");
		histoPi0UB2760GeV_Pythia8_Tune4CptHard2->Draw("pE1same");
		histoPi0UB2760GeV_Pythia8_Tune4CptHard5->Draw("pE1same");
		histoPi0UB2760GeV_Pythia8_Tune4CptHard10->Draw("pE1same");

		labelEnergyXSectionPi0->Draw();
		labelDetSysXSectionPi0->Draw();
		labelPythiaTuneA->Draw();
		legendXSectionPi0->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/Tune4C_Pi0UBXSectionDiffptHardBins.eps"));

	histo2DXSectionPi0->Draw("copy");
		histoPi0UB2760GeV_Pythia8_TuneMonachMB->Draw("pE1same");
		histoPi0UB2760GeV_Pythia8_TuneMonachptHard2->Draw("pE1same");
		histoPi0UB2760GeV_Pythia8_TuneMonachptHard5->Draw("pE1same");
		histoPi0UB2760GeV_Pythia8_TuneMonachptHard10->Draw("pE1same");

		labelEnergyXSectionPi0->Draw();
		labelDetSysXSectionPi0->Draw();
		labelPythiaTuneB->Draw();
		legendXSectionPi0->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/TuneMonach_Pi0UBXSectionDiffptHardBins.eps"));
	
	histo2DXSectionPi0->Draw("copy");
		histoPi02760GeV_Pythia8_Tune4CMB->Draw("pE1same");
		histoPi0FromEta2760GeV_Pythia8_Tune4CMB->Draw("pE1same");
//         histoPi0FromK0s2760GeV_Pythia8_Tune4CMB->Draw("pE1same");
		histoPi0FromLambda2760GeV_Pythia8_Tune4CMB->Draw("pE1same");

		labelEnergyXSectionPi0->Draw();
		labelDetSysXSectionPi0->Draw();
		labelPythiaTuneA->Draw();
		
		TLegend* legendXSectionPi02MB = new TLegend(0.62,0.72,0.9,0.82);
		legendXSectionPi02MB->SetFillColor(0);
		legendXSectionPi02MB->SetLineColor(0);
		legendXSectionPi02MB->SetTextFont(42);
		legendXSectionPi02MB->SetTextSize(0.035);
		legendXSectionPi02MB->AddEntry(histoPi02760GeV_Pythia8_Tune4CMB,"X (not K)#rightarrow #pi^{0}","p");
//         legendXSectionPi02MB->AddEntry(histoPi02760GeV_Pythia8_Tune4CptHard2,"K #rightarrow #pi^{0}","p");
		legendXSectionPi02MB->AddEntry(histoPi0FromEta2760GeV_Pythia8_Tune4CMB,"#eta #rightarrow #pi^{0}","p");
		legendXSectionPi02MB->AddEntry(histoPi0FromLambda2760GeV_Pythia8_Tune4CMB,"#Lambda #rightarrow #pi^{0}","p");
		legendXSectionPi02MB->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/Tune4C_Pi0XSectionDiffSourcesMB.eps"));

	histo2DXSectionPi0->Draw("copy");
		TH1D* histoPi0All2760GeV_Pythia8_Tune4CMB = (TH1D*)histoPi02760GeV_Pythia8_Tune4CMB->Clone("histoPi0All2760GeV_Pythia8_Tune4CMB");
		histoPi0All2760GeV_Pythia8_Tune4CMB->Sumw2();
		histoPi0All2760GeV_Pythia8_Tune4CMB->Add(histoPi0FromK0s2760GeV_Pythia8_Tune4CMB);

		DrawGammaSetMarker(histoPi0All2760GeV_Pythia8_Tune4CMB, 20, 1, colorMB, colorMB);

		histoPi0All2760GeV_Pythia8_Tune4CMB->Draw("pE1same");
		histoPi0FromK0s2760GeV_Pythia8_Tune4CMB->Draw("pE1same");
		histoPi0FromEta2760GeV_Pythia8_Tune4CMB->Draw("pE1same");
		histoPi0FromLambda2760GeV_Pythia8_Tune4CMB->Draw("pE1same");

		labelEnergyXSectionPi0->Draw();
		labelDetSysXSectionPi0->Draw();
		labelPythiaTuneA->Draw();
		TLegend* legendXSectionPi03MB = new TLegend(0.62,0.66,0.9,0.82);
		legendXSectionPi03MB->SetFillColor(0);
		legendXSectionPi03MB->SetLineColor(0);
		legendXSectionPi03MB->SetTextFont(42);
		legendXSectionPi03MB->SetTextSize(0.035);
		legendXSectionPi03MB->AddEntry(histoPi0All2760GeV_Pythia8_Tune4CMB,"X #rightarrow #pi^{0}","p");
		legendXSectionPi03MB->AddEntry(histoPi0FromK0s2760GeV_Pythia8_Tune4CMB,"K #rightarrow #pi^{0}","p");
		legendXSectionPi03MB->AddEntry(histoPi0FromEta2760GeV_Pythia8_Tune4CMB,"#eta #rightarrow #pi^{0}","p");
		legendXSectionPi03MB->AddEntry(histoPi0FromLambda2760GeV_Pythia8_Tune4CMB,"#Lambda #rightarrow #pi^{0}","p");
		legendXSectionPi03MB->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/Tune4C_Pi0AllXSectionDiffSourcesMB.eps"));


	histo2DXSectionPi0->Draw("copy");
		histoPi02760GeV_Pythia8_TuneMonachMB->Draw("pE1same");
		histoPi0FromEta2760GeV_Pythia8_TuneMonachMB->Draw("pE1same");
//         histoPi0FromK0s2760GeV_Pythia8_TuneMonachMB->Draw("pE1same");
		histoPi0FromLambda2760GeV_Pythia8_TuneMonachMB->Draw("pE1same");

		labelEnergyXSectionPi0->Draw();
		labelDetSysXSectionPi0->Draw();
		labelPythiaTuneB->Draw();
		
		legendXSectionPi02MB->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/TuneMonach_Pi0XSectionDiffSourcesMB.eps"));

	histo2DXSectionPi0->Draw("copy");
		TH1D* histoPi0All2760GeV_Pythia8_TuneMonachMB = (TH1D*)histoPi02760GeV_Pythia8_TuneMonachMB->Clone("histoPi0All2760GeV_Pythia8_TuneMonachMB");
		histoPi0All2760GeV_Pythia8_TuneMonachMB->Sumw2();
		histoPi0All2760GeV_Pythia8_TuneMonachMB->Add(histoPi0FromK0s2760GeV_Pythia8_TuneMonachMB);

		DrawGammaSetMarker(histoPi0All2760GeV_Pythia8_TuneMonachMB, 20, 1, colorMB, colorMB);

		histoPi0All2760GeV_Pythia8_TuneMonachMB->Draw("pE1same");
		histoPi0FromK0s2760GeV_Pythia8_TuneMonachMB->Draw("pE1same");
		histoPi0FromEta2760GeV_Pythia8_TuneMonachMB->Draw("pE1same");
		histoPi0FromLambda2760GeV_Pythia8_TuneMonachMB->Draw("pE1same");

		labelEnergyXSectionPi0->Draw();
		labelDetSysXSectionPi0->Draw();
		labelPythiaTuneB->Draw();
		legendXSectionPi03MB->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/TuneMonach_Pi0AllXSectionDiffSourcesMB.eps"));
	
	histo2DXSectionPi0->Draw("copy");
		histoPi02760GeV_Pythia8_Tune4CptHard2->Draw("pE1same");
		histoPi0FromEta2760GeV_Pythia8_Tune4CptHard2->Draw("pE1same");
//         histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard2->Draw("pE1same");
		histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard2->Draw("pE1same");

		labelEnergyXSectionPi0->Draw();
		labelDetSysXSectionPi0->Draw();
		labelPythiaTuneA->Draw();
		TLegend* legendXSectionPi02ptHard2 = new TLegend(0.62,0.72,0.9,0.82);
		legendXSectionPi02ptHard2->SetFillColor(0);
		legendXSectionPi02ptHard2->SetLineColor(0);
		legendXSectionPi02ptHard2->SetTextFont(42);
		legendXSectionPi02ptHard2->SetTextSize(0.035);
		legendXSectionPi02ptHard2->AddEntry(histoPi02760GeV_Pythia8_Tune4CptHard2,"X (not K)#rightarrow #pi^{0}","p");
//         legendXSectionPi02ptHard2->AddEntry(histoPi02760GeV_Pythia8_Tune4CptHard2,"K #rightarrow #pi^{0}","p");
		legendXSectionPi02ptHard2->AddEntry(histoPi0FromEta2760GeV_Pythia8_Tune4CptHard2,"#eta #rightarrow #pi^{0}","p");
		legendXSectionPi02ptHard2->AddEntry(histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard2,"#Lambda #rightarrow #pi^{0}","p");
		legendXSectionPi02ptHard2->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/Tune4C_Pi0XSectionDiffSourcesptHard2.eps"));

	histo2DXSectionPi0->Draw("copy");
		TH1D* histoPi0All2760GeV_Pythia8_Tune4CptHard2 = (TH1D*)histoPi02760GeV_Pythia8_Tune4CptHard2->Clone("histoPi0All2760GeV_Pythia8_Tune4CptHard2");
		histoPi0All2760GeV_Pythia8_Tune4CptHard2->Sumw2();
		histoPi0All2760GeV_Pythia8_Tune4CptHard2->Add(histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard2);

		DrawGammaSetMarker(histoPi0All2760GeV_Pythia8_Tune4CptHard2, 20, 1, colorptHard2, colorptHard2);

		histoPi0All2760GeV_Pythia8_Tune4CptHard2->Draw("pE1same");
		histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard2->Draw("pE1same");
		histoPi0FromEta2760GeV_Pythia8_Tune4CptHard2->Draw("pE1same");
		histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard2->Draw("pE1same");

		labelEnergyXSectionPi0->Draw();
		labelDetSysXSectionPi0->Draw();
		labelPythiaTuneA->Draw();
		
		TLegend* legendXSectionPi03ptHard2 = new TLegend(0.62,0.66,0.9,0.82);
		legendXSectionPi03ptHard2->SetFillColor(0);
		legendXSectionPi03ptHard2->SetLineColor(0);
		legendXSectionPi03ptHard2->SetTextFont(42);
		legendXSectionPi03ptHard2->SetTextSize(0.035);
		legendXSectionPi03ptHard2->AddEntry(histoPi0All2760GeV_Pythia8_Tune4CptHard2,"X #rightarrow #pi^{0}","p");
		legendXSectionPi03ptHard2->AddEntry(histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard2,"K #rightarrow #pi^{0}","p");
		legendXSectionPi03ptHard2->AddEntry(histoPi0FromEta2760GeV_Pythia8_Tune4CptHard2,"#eta #rightarrow #pi^{0}","p");
		legendXSectionPi03ptHard2->AddEntry(histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard2,"#Lambda #rightarrow #pi^{0}","p");
		legendXSectionPi03ptHard2->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/Tune4C_Pi0AllXSectionDiffSourcesptHard2.eps"));

	histo2DXSectionPi0->Draw("copy");
		histoPi02760GeV_Pythia8_TuneMonachptHard2->Draw("pE1same");
		histoPi0FromEta2760GeV_Pythia8_TuneMonachptHard2->Draw("pE1same");
//         histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard2->Draw("pE1same");
		histoPi0FromLambda2760GeV_Pythia8_TuneMonachptHard2->Draw("pE1same");

		labelEnergyXSectionPi0->Draw();
		labelDetSysXSectionPi0->Draw();
		labelPythiaTuneB->Draw();
		legendXSectionPi02ptHard2->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/TuneMonach_Pi0XSectionDiffSourcesptHard2.eps"));

	histo2DXSectionPi0->Draw("copy");
		TH1D* histoPi0All2760GeV_Pythia8_TuneMonachptHard2 = (TH1D*)histoPi02760GeV_Pythia8_TuneMonachptHard2->Clone("histoPi0All2760GeV_Pythia8_TuneMonachptHard2");
		histoPi0All2760GeV_Pythia8_TuneMonachptHard2->Sumw2();
		histoPi0All2760GeV_Pythia8_TuneMonachptHard2->Add(histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard2);

		DrawGammaSetMarker(histoPi0All2760GeV_Pythia8_TuneMonachptHard2, 20, 1, colorptHard2, colorptHard2);

		histoPi0All2760GeV_Pythia8_TuneMonachptHard2->Draw("pE1same");
		histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard2->Draw("pE1same");
		histoPi0FromEta2760GeV_Pythia8_TuneMonachptHard2->Draw("pE1same");
		histoPi0FromLambda2760GeV_Pythia8_TuneMonachptHard2->Draw("pE1same");

		labelEnergyXSectionPi0->Draw();
		labelDetSysXSectionPi0->Draw();
		labelPythiaTuneB->Draw();
		legendXSectionPi03ptHard2->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/TuneMonach_Pi0AllXSectionDiffSourcesptHard2.eps"));
	
	histo2DXSectionPi0->Draw("copy");
		histoPi02760GeV_Pythia8_Tune4CptHard5->Draw("pE1same");
		histoPi0FromEta2760GeV_Pythia8_Tune4CptHard5->Draw("pE1same");
//         histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard5->Draw("pE1same");
		histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard5->Draw("pE1same");

		labelEnergyXSectionPi0->Draw();
		labelDetSysXSectionPi0->Draw();
		labelPythiaTuneA->Draw();
		
		TLegend* legendXSectionPi02ptHard5 = new TLegend(0.62,0.72,0.9,0.82);
		legendXSectionPi02ptHard5->SetFillColor(0);
		legendXSectionPi02ptHard5->SetLineColor(0);
		legendXSectionPi02ptHard5->SetTextFont(42);
		legendXSectionPi02ptHard5->SetTextSize(0.035);
		legendXSectionPi02ptHard5->AddEntry(histoPi02760GeV_Pythia8_Tune4CptHard5,"X (not K)#rightarrow #pi^{0}","p");
//         legendXSectionPi02ptHard5->AddEntry(histoPi02760GeV_Pythia8_Tune4CptHard2,"K #rightarrow #pi^{0}","p");
		legendXSectionPi02ptHard5->AddEntry(histoPi0FromEta2760GeV_Pythia8_Tune4CptHard5,"#eta #rightarrow #pi^{0}","p");
		legendXSectionPi02ptHard5->AddEntry(histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard5,"#Lambda #rightarrow #pi^{0}","p");
		legendXSectionPi02ptHard5->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/Tune4C_Pi0XSectionDiffSourcesptHard5.eps"));

	histo2DXSectionPi0->Draw("copy");
		TH1D* histoPi0All2760GeV_Pythia8_Tune4CptHard5 = (TH1D*)histoPi02760GeV_Pythia8_Tune4CptHard5->Clone("histoPi0All2760GeV_Pythia8_Tune4CptHard5");
		histoPi0All2760GeV_Pythia8_Tune4CptHard5->Sumw2();
		histoPi0All2760GeV_Pythia8_Tune4CptHard5->Add(histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard5);

		DrawGammaSetMarker(histoPi0All2760GeV_Pythia8_Tune4CptHard5, 20, 1, colorptHard5, colorptHard5);

		histoPi0All2760GeV_Pythia8_Tune4CptHard5->Draw("pE1same");
		histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard5->Draw("pE1same");
		histoPi0FromEta2760GeV_Pythia8_Tune4CptHard5->Draw("pE1same");
		histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard5->Draw("pE1same");

		labelEnergyXSectionPi0->Draw();
		labelDetSysXSectionPi0->Draw();
		labelPythiaTuneA->Draw();
		TLegend* legendXSectionPi03ptHard5 = new TLegend(0.62,0.66,0.9,0.82);
		legendXSectionPi03ptHard5->SetFillColor(0);
		legendXSectionPi03ptHard5->SetLineColor(0);
		legendXSectionPi03ptHard5->SetTextFont(42);
		legendXSectionPi03ptHard5->SetTextSize(0.035);
		legendXSectionPi03ptHard5->AddEntry(histoPi0All2760GeV_Pythia8_Tune4CptHard5,"X #rightarrow #pi^{0}","p");
		legendXSectionPi03ptHard5->AddEntry(histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard5,"K #rightarrow #pi^{0}","p");
		legendXSectionPi03ptHard5->AddEntry(histoPi0FromEta2760GeV_Pythia8_Tune4CptHard5,"#eta #rightarrow #pi^{0}","p");
		legendXSectionPi03ptHard5->AddEntry(histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard5,"#Lambda #rightarrow #pi^{0}","p");
		legendXSectionPi03ptHard5->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/Tune4C_Pi0AllXSectionDiffSourcesptHard5.eps"));

	histo2DXSectionPi0->Draw("copy");
		histoPi02760GeV_Pythia8_TuneMonachptHard5->Draw("pE1same");
		histoPi0FromEta2760GeV_Pythia8_TuneMonachptHard5->Draw("pE1same");
//         histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard5->Draw("pE1same");
		histoPi0FromLambda2760GeV_Pythia8_TuneMonachptHard5->Draw("pE1same");

		labelEnergyXSectionPi0->Draw();
		labelDetSysXSectionPi0->Draw();
		labelPythiaTuneB->Draw();   
		legendXSectionPi02ptHard5->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/TuneMonach_Pi0XSectionDiffSourcesptHard5.eps"));

	histo2DXSectionPi0->Draw("copy");
		TH1D* histoPi0All2760GeV_Pythia8_TuneMonachptHard5 = (TH1D*)histoPi02760GeV_Pythia8_TuneMonachptHard5->Clone("histoPi0All2760GeV_Pythia8_TuneMonachptHard5");
		histoPi0All2760GeV_Pythia8_TuneMonachptHard5->Sumw2();
		histoPi0All2760GeV_Pythia8_TuneMonachptHard5->Add(histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard5);

		DrawGammaSetMarker(histoPi0All2760GeV_Pythia8_TuneMonachptHard5, 20, 1, colorptHard5, colorptHard5);

		histoPi0All2760GeV_Pythia8_TuneMonachptHard5->Draw("pE1same");
		histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard5->Draw("pE1same");
		histoPi0FromEta2760GeV_Pythia8_TuneMonachptHard5->Draw("pE1same");
		histoPi0FromLambda2760GeV_Pythia8_TuneMonachptHard5->Draw("pE1same");

		labelEnergyXSectionPi0->Draw();
		labelDetSysXSectionPi0->Draw();
		labelPythiaTuneB->Draw();   
		legendXSectionPi03ptHard5->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/TuneMonach_Pi0AllXSectionDiffSourcesptHard5.eps"));
	
	histo2DXSectionPi0->Draw("copy");
		histoPi02760GeV_Pythia8_Tune4CptHard10->Draw("pE1same");
		histoPi0FromEta2760GeV_Pythia8_Tune4CptHard10->Draw("pE1same");
//         histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard10->Draw("pE1same");
		histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard10->Draw("pE1same");

		labelEnergyXSectionPi0->Draw();
		labelDetSysXSectionPi0->Draw();
		labelPythiaTuneA->Draw();
		
		TLegend* legendXSectionPi02ptHard10 = new TLegend(0.62,0.72,0.9,0.82);
		legendXSectionPi02ptHard10->SetFillColor(0);
		legendXSectionPi02ptHard10->SetLineColor(0);
		legendXSectionPi02ptHard10->SetTextFont(42);
		legendXSectionPi02ptHard10->SetTextSize(0.035);
		legendXSectionPi02ptHard10->AddEntry(histoPi02760GeV_Pythia8_Tune4CptHard10,"X (not K)#rightarrow #pi^{0}","p");
//         legendXSectionPi02ptHard10->AddEntry(histoPi02760GeV_Pythia8_Tune4CptHard2,"K #rightarrow #pi^{0}","p");
		legendXSectionPi02ptHard10->AddEntry(histoPi0FromEta2760GeV_Pythia8_Tune4CptHard10,"#eta #rightarrow #pi^{0}","p");
		legendXSectionPi02ptHard10->AddEntry(histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard10,"#Lambda #rightarrow #pi^{0}","p");
		legendXSectionPi02ptHard10->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/Tune4C_Pi0XSectionDiffSourcesptHard10.eps"));

	histo2DXSectionPi0->Draw("copy");
		TH1D* histoPi0All2760GeV_Pythia8_Tune4CptHard10 = (TH1D*)histoPi02760GeV_Pythia8_Tune4CptHard10->Clone("histoPi0All2760GeV_Pythia8_Tune4CptHard10");
		histoPi0All2760GeV_Pythia8_Tune4CptHard10->Sumw2();
		histoPi0All2760GeV_Pythia8_Tune4CptHard10->Add(histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard10);

		DrawGammaSetMarker(histoPi0All2760GeV_Pythia8_Tune4CptHard10, 20, 1, colorptHard10, colorptHard10);

		histoPi0All2760GeV_Pythia8_Tune4CptHard10->Draw("pE1same");
		histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard10->Draw("pE1same");
		histoPi0FromEta2760GeV_Pythia8_Tune4CptHard10->Draw("pE1same");
		histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard10->Draw("pE1same");

		labelEnergyXSectionPi0->Draw();
		labelDetSysXSectionPi0->Draw();
		labelPythiaTuneA->Draw();
		
		TLegend* legendXSectionPi03ptHard10 = new TLegend(0.62,0.66,0.9,0.82);
		legendXSectionPi03ptHard10->SetFillColor(0);
		legendXSectionPi03ptHard10->SetLineColor(0);
		legendXSectionPi03ptHard10->SetTextFont(42);
		legendXSectionPi03ptHard10->SetTextSize(0.035);
		legendXSectionPi03ptHard10->AddEntry(histoPi0All2760GeV_Pythia8_Tune4CptHard10,"X #rightarrow #pi^{0}","p");
		legendXSectionPi03ptHard10->AddEntry(histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard10,"K #rightarrow #pi^{0}","p");
		legendXSectionPi03ptHard10->AddEntry(histoPi0FromEta2760GeV_Pythia8_Tune4CptHard10,"#eta #rightarrow #pi^{0}","p");
		legendXSectionPi03ptHard10->AddEntry(histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard10,"#Lambda #rightarrow #pi^{0}","p");
		legendXSectionPi03ptHard10->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/Tune4C_Pi0AllXSectionDiffSourcesptHard10.eps"));

	histo2DXSectionPi0->Draw("copy");
		histoPi02760GeV_Pythia8_TuneMonachptHard10->Draw("pE1same");
		histoPi0FromEta2760GeV_Pythia8_TuneMonachptHard10->Draw("pE1same");
//         histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard10->Draw("pE1same");
		histoPi0FromLambda2760GeV_Pythia8_TuneMonachptHard10->Draw("pE1same");

		labelEnergyXSectionPi0->Draw();
		labelDetSysXSectionPi0->Draw();
		labelPythiaTuneB->Draw();
		legendXSectionPi02ptHard10->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/TuneMonach_Pi0XSectionDiffSourcesptHard10.eps"));

	histo2DXSectionPi0->Draw("copy");
		TH1D* histoPi0All2760GeV_Pythia8_TuneMonachptHard10 = (TH1D*)histoPi02760GeV_Pythia8_TuneMonachptHard10->Clone("histoPi0All2760GeV_Pythia8_TuneMonachptHard10");
		histoPi0All2760GeV_Pythia8_TuneMonachptHard10->Sumw2();
		histoPi0All2760GeV_Pythia8_TuneMonachptHard10->Add(histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard10);

		DrawGammaSetMarker(histoPi0All2760GeV_Pythia8_TuneMonachptHard10, 20, 1, colorptHard10, colorptHard10);

		histoPi0All2760GeV_Pythia8_TuneMonachptHard10->Draw("pE1same");
		histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard10->Draw("pE1same");
		histoPi0FromEta2760GeV_Pythia8_TuneMonachptHard10->Draw("pE1same");
		histoPi0FromLambda2760GeV_Pythia8_TuneMonachptHard10->Draw("pE1same");

		labelEnergyXSectionPi0->Draw();
		labelDetSysXSectionPi0->Draw();
		labelPythiaTuneA->Draw();
		legendXSectionPi03ptHard10->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/TuneMonach_Pi0AllXSectionDiffSourcesptHard10.eps"));
	
	
	histo2DXSectionPi0->Draw("copy");

		histoEta2760GeV_Pythia8_Tune4CMB->Draw("pE1same");
		histoEta2760GeV_Pythia8_Tune4CptHard2->Draw("pE1same");
		histoEta2760GeV_Pythia8_Tune4CptHard5->Draw("pE1same");
		histoEta2760GeV_Pythia8_Tune4CptHard10->Draw("pE1same");

		TLatex *labelEnergyXSectionEta = new TLatex(0.64,0.92,collisionSystem2760GeV.Data());
		SetStyleTLatex( labelEnergyXSectionEta, 0.035,4);
		labelEnergyXSectionEta->Draw();
		TLatex *labelDetSysXSectionEta = new TLatex(0.64,0.88,"#pi^{0} #rightarrow #gamma#gamma");
		SetStyleTLatex( labelDetSysXSectionEta, 0.035,4);
		labelDetSysXSectionEta->Draw();
		labelPythiaTuneA->Draw();
		
		TLegend* legendXSectionEta = new TLegend(0.62,0.66,0.9,0.82);
		legendXSectionEta->SetFillColor(0);
		legendXSectionEta->SetLineColor(0);
		legendXSectionEta->SetTextFont(42);
		legendXSectionEta->SetTextSize(0.035);
		legendXSectionEta->AddEntry(histoEta2760GeV_Pythia8_Tune4CMB,"MB","p");
		legendXSectionEta->AddEntry(histoEta2760GeV_Pythia8_Tune4CptHard2,"p_{T,hard} > 2 GeV/c","p");
		legendXSectionEta->AddEntry(histoEta2760GeV_Pythia8_Tune4CptHard5,"p_{T,hard} > 5 GeV/c","p");
		legendXSectionEta->AddEntry(histoEta2760GeV_Pythia8_Tune4CptHard10,"p_{T,hard} > 10 GeV/c","p");
		legendXSectionEta->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/Tune4C_EtaXSectionDiffptHardBins.eps"));
	
	histo2DXSectionPi0->Draw("copy");
		histoEtaUB2760GeV_Pythia8_Tune4CMB->Draw("pE1same");
		histoEtaUB2760GeV_Pythia8_Tune4CptHard2->Draw("pE1same");
		histoEtaUB2760GeV_Pythia8_Tune4CptHard5->Draw("pE1same");
		histoEtaUB2760GeV_Pythia8_Tune4CptHard10->Draw("pE1same");

		labelEnergyXSectionEta->Draw();
		labelDetSysXSectionEta->Draw();
		labelPythiaTuneA->Draw();
		legendXSectionEta->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/Tune4C_EtaUBXSectionDiffptHardBins.eps"));

	histo2DXSectionPi0->Draw("copy");

		histoEta2760GeV_Pythia8_TuneMonachMB->Draw("pE1same");
		histoEta2760GeV_Pythia8_TuneMonachptHard2->Draw("pE1same");
		histoEta2760GeV_Pythia8_TuneMonachptHard5->Draw("pE1same");
		histoEta2760GeV_Pythia8_TuneMonachptHard10->Draw("pE1same");

		labelEnergyXSectionEta->Draw();
		labelDetSysXSectionEta->Draw();
		labelPythiaTuneB->Draw();
		legendXSectionEta->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/TuneMonach_EtaXSectionDiffptHardBins.eps"));
	
	histo2DXSectionPi0->Draw("copy");
		histoEtaUB2760GeV_Pythia8_TuneMonachMB->Draw("pE1same");
		histoEtaUB2760GeV_Pythia8_TuneMonachptHard2->Draw("pE1same");
		histoEtaUB2760GeV_Pythia8_TuneMonachptHard5->Draw("pE1same");
		histoEtaUB2760GeV_Pythia8_TuneMonachptHard10->Draw("pE1same");

		labelEnergyXSectionEta->Draw();
		labelDetSysXSectionEta->Draw();
		labelPythiaTuneB->Draw();
		legendXSectionEta->Draw();

	canvasXSectionPi0->SaveAs(Form("ExternalInput/Theory/Pythia/TuneMonach_EtaUBXSectionDiffptHardBins.eps"));
	
	
	TCanvas* canvasRatioXSection = new TCanvas("canvasRatioXSection","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioXSection, 0.1, 0.02, 0.035, 0.09);
	canvasRatioXSection->SetLogx();

	TH2F * histo2DRatioXSection;
	histo2DRatioXSection = new TH2F("histo2DRatioXSection","histo2DRatioXSection",11000,0.23,70.,1000,0,2);
	SetStyleHistoTH2ForGraphs(histo2DRatioXSection, "#it{p}_{T} (GeV/#it{c})","#frac{#sigma_{A}}{#sigma_{B}}",0.035,0.04, 0.035,0.04, 1.,1.);
	histo2DRatioXSection->GetXaxis()->SetMoreLogLabels();
	histo2DRatioXSection->GetXaxis()->SetLabelOffset(-0.01);
//     histo2DRatioXSection->GetYaxis()->SetRangeUser(-10,10);
	histo2DRatioXSection->Draw("copy");

		ratioPi02760GeV_Pythia8_Tune4CpTHard2DivMB->Draw("pE1same");
		ratioPi02760GeV_Pythia8_Tune4CpTHard5DivpTHard2->Draw("pE1same");
		ratioPi02760GeV_Pythia8_Tune4CpTHard10DivpTHard5->Draw("pE1same");

		TF1* const1Pi0Tune4C = new TF1("const1Pi0Tune4C","[0]",4.0,20);
		TF1* const2Pi0Tune4C = new TF1("const2Pi0Tune4C","[0]",7.0,30);
		TF1* const3Pi0Tune4C = new TF1("const3Pi0Tune4C","[0]",12.,50);
		
		ratioPi02760GeV_Pythia8_Tune4CpTHard2DivMB->Fit(const1Pi0Tune4C,"NRME+","",4,15);
		ratioPi02760GeV_Pythia8_Tune4CpTHard5DivpTHard2->Fit(const2Pi0Tune4C,"NRME+","",7,30);
		ratioPi02760GeV_Pythia8_Tune4CpTHard10DivpTHard5->Fit(const3Pi0Tune4C,"NRME+","",12,50);
		
		SetStyleFit(const1Pi0Tune4C, 4.0, 20, 2, 3, colorptHard2+1); 
		SetStyleFit(const2Pi0Tune4C, 7,30, 2, 5, colorptHard5+1);
		SetStyleFit(const3Pi0Tune4C, 12,50, 2, 7, colorptHard10+1);
		
		const1Pi0Tune4C->Draw("same");
		const2Pi0Tune4C->Draw("same");
		const3Pi0Tune4C->Draw("same");
		
		TLatex *labelEnergyRatioXSectionPi0 = new TLatex(0.74,0.90,collisionSystem2760GeV.Data());
		SetStyleTLatex( labelEnergyRatioXSectionPi0, 0.035,4);
		labelEnergyRatioXSectionPi0->Draw();
		TLatex *labelDetSysRatioXSectionPi0 = new TLatex(0.74,0.86,"#pi^{0} #rightarrow #gamma#gamma");
		SetStyleTLatex( labelDetSysRatioXSectionPi0, 0.035,4);
		labelDetSysRatioXSectionPi0->Draw();
		TLatex *labelRatioPythiaTuneA = new TLatex(0.74,0.82,"Pythia 8.1, Tune 4C");
		SetStyleTLatex( labelRatioPythiaTuneA, 0.035,4);
		labelRatioPythiaTuneA->Draw();

		TLegend* legendXSectionRatioPi0 = new TLegend(0.12,0.8,0.6,0.94);
		legendXSectionRatioPi0->SetFillColor(0);
		legendXSectionRatioPi0->SetLineColor(0);
		legendXSectionRatioPi0->SetTextFont(42);
		legendXSectionRatioPi0->SetTextSize(0.035);
		legendXSectionRatioPi0->SetNColumns(2);
		legendXSectionRatioPi0->SetMargin(0.1);
		legendXSectionRatioPi0->AddEntry(ratioPi02760GeV_Pythia8_Tune4CpTHard2DivMB,"(p_{T,hard} > 2 GeV/c)/(MB)","p");
		legendXSectionRatioPi0->AddEntry(const1Pi0Tune4C,Form("%0.2f",const1Pi0Tune4C->GetParameter(0)),"l");
		legendXSectionRatioPi0->AddEntry(ratioPi02760GeV_Pythia8_Tune4CpTHard5DivpTHard2,"(p_{T,hard} > 5 GeV/c)/((p_{T,hard} > 2 GeV/c))","p");
		legendXSectionRatioPi0->AddEntry(const2Pi0Tune4C,Form("%0.2f",const2Pi0Tune4C->GetParameter(0)),"l");
		legendXSectionRatioPi0->AddEntry(ratioPi02760GeV_Pythia8_Tune4CpTHard10DivpTHard5,"(p_{T,hard} > 10 GeV/c)/((p_{T,hard} > 5 GeV/c))","p");
		legendXSectionRatioPi0->AddEntry(const3Pi0Tune4C,Form("%0.2f",const3Pi0Tune4C->GetParameter(0)),"l");
		legendXSectionRatioPi0->Draw();

		
	canvasRatioXSection->SaveAs("ExternalInput/Theory/Pythia/Tune4C_RatioPi0XSectionDiffptHardBins.eps");

	histo2DRatioXSection->Draw("copy");

		ratioPi0UB2760GeV_Pythia8_Tune4CpTHard2DivMB->Draw("pE1same");
		ratioPi0UB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2->Draw("pE1same");
		ratioPi0UB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5->Draw("pE1same");

		const1Pi0Tune4C->Draw("same");
		const2Pi0Tune4C->Draw("same");
		const3Pi0Tune4C->Draw("same");
		
		labelEnergyRatioXSectionPi0->Draw();
		labelDetSysRatioXSectionPi0->Draw();
		labelRatioPythiaTuneA->Draw();
		
		legendXSectionRatioPi0->Draw();   
	canvasRatioXSection->SaveAs("ExternalInput/Theory/Pythia/Tune4C_RatioPi0UBXSectionDiffptHardBins.eps");
	
	histo2DRatioXSection->Draw("copy");

		ratioPi02760GeV_Pythia8_TuneMonachpTHard2DivMB->Draw("pE1same");
		ratioPi02760GeV_Pythia8_TuneMonachpTHard5DivpTHard2->Draw("pE1same");
		ratioPi02760GeV_Pythia8_TuneMonachpTHard10DivpTHard5->Draw("pE1same");

		TF1* const1Pi0TuneMonach = new TF1("const1Pi0TuneMonach","[0]",4.0,20);
		TF1* const2Pi0TuneMonach = new TF1("const2Pi0TuneMonach","[0]",7.0,30);
		TF1* const3Pi0TuneMonach = new TF1("const3Pi0TuneMonach","[0]",12.,50);
		
		ratioPi02760GeV_Pythia8_TuneMonachpTHard2DivMB->Fit(const1Pi0TuneMonach,"NRME+","",4,15);
		ratioPi02760GeV_Pythia8_TuneMonachpTHard5DivpTHard2->Fit(const2Pi0TuneMonach,"NRME+","",7,30);
		ratioPi02760GeV_Pythia8_TuneMonachpTHard10DivpTHard5->Fit(const3Pi0TuneMonach,"NRME+","",12,50);
		
		SetStyleFit(const1Pi0TuneMonach, 4.0, 20, 2, 3, colorptHard2+1); 
		SetStyleFit(const2Pi0TuneMonach, 7,30, 2, 5, colorptHard5+1);
		SetStyleFit(const3Pi0TuneMonach, 12,50, 2, 7, colorptHard10+1);
		
		const1Pi0TuneMonach->Draw("same");
		const2Pi0TuneMonach->Draw("same");
		const3Pi0TuneMonach->Draw("same");
		
		labelEnergyRatioXSectionPi0->Draw();
		labelDetSysRatioXSectionPi0->Draw();
		TLatex *labelRatioPythiaTuneB = new TLatex(0.74,0.82,"Pythia 8.2, Monach");
		SetStyleTLatex( labelRatioPythiaTuneB, 0.035,4);
		labelRatioPythiaTuneB->Draw();

		TLegend* legendXSectionRatioPi0TuneMonach = new TLegend(0.12,0.8,0.6,0.94);
		legendXSectionRatioPi0TuneMonach->SetFillColor(0);
		legendXSectionRatioPi0TuneMonach->SetLineColor(0);
		legendXSectionRatioPi0TuneMonach->SetTextFont(42);
		legendXSectionRatioPi0TuneMonach->SetTextSize(0.035);
		legendXSectionRatioPi0TuneMonach->SetNColumns(2);
		legendXSectionRatioPi0TuneMonach->SetMargin(0.1);
		legendXSectionRatioPi0TuneMonach->AddEntry(ratioPi02760GeV_Pythia8_TuneMonachpTHard2DivMB,"(p_{T,hard} > 2 GeV/c)/(MB)","p");
		legendXSectionRatioPi0TuneMonach->AddEntry(const1Pi0TuneMonach,Form("%0.2f",const1Pi0TuneMonach->GetParameter(0)),"l");
		legendXSectionRatioPi0TuneMonach->AddEntry(ratioPi02760GeV_Pythia8_TuneMonachpTHard5DivpTHard2,"(p_{T,hard} > 5 GeV/c)/((p_{T,hard} > 2 GeV/c))","p");
		legendXSectionRatioPi0TuneMonach->AddEntry(const2Pi0TuneMonach,Form("%0.2f",const2Pi0TuneMonach->GetParameter(0)),"l");
		legendXSectionRatioPi0TuneMonach->AddEntry(ratioPi02760GeV_Pythia8_TuneMonachpTHard10DivpTHard5,"(p_{T,hard} > 10 GeV/c)/((p_{T,hard} > 5 GeV/c))","p");
		legendXSectionRatioPi0TuneMonach->AddEntry(const3Pi0TuneMonach,Form("%0.2f",const3Pi0TuneMonach->GetParameter(0)),"l");
		legendXSectionRatioPi0TuneMonach->Draw();

		
	canvasRatioXSection->SaveAs("ExternalInput/Theory/Pythia/TuneMonach_RatioPi0XSectionDiffptHardBins.eps");

	histo2DRatioXSection->Draw("copy");

		ratioPi0UB2760GeV_Pythia8_TuneMonachpTHard2DivMB->Draw("pE1same");
		ratioPi0UB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2->Draw("pE1same");
		ratioPi0UB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5->Draw("pE1same");

		const1Pi0TuneMonach->Draw("same");
		const2Pi0TuneMonach->Draw("same");
		const3Pi0TuneMonach->Draw("same");
		
		labelEnergyRatioXSectionPi0->Draw();
		labelDetSysRatioXSectionPi0->Draw();
		labelRatioPythiaTuneB->Draw();
		
		legendXSectionRatioPi0TuneMonach->Draw();   
	canvasRatioXSection->SaveAs("ExternalInput/Theory/Pythia/TuneMonach_RatioPi0UBXSectionDiffptHardBins.eps");
	
	
	
	Double_t    xSection2760GeVINEL =    62.8*1e-3;
	Double_t    xSection7000GeVINEL =    73.2*1e-3;

	index = 0;
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
	TGraph * graphInvCrossSecCGC7000GeV = ScaleGraph(graphInvYieldCGC7000GeV,(xSection7000GeVINEL*recalcBarn));
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
	TGraph * graphInvCrossSecCGC7000GeV_mvgamma = ScaleGraph(graphInvYieldCGC7000GeV_mvgamma,(xSection7000GeVINEL*recalcBarn));

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
	TGraph * graphInvCrossSecCGC2760GeV = ScaleGraph(graphInvYieldCGC2760GeV,(xSection2760GeVINEL*recalcBarn));

	Double_t       ptNLODSS14Pi07000GeV[100];
	Double_t       ptErrNLODSS14Pi07000GeV[100];
	Double_t       muHalfDSS14Pi07000GeV[100];
	Double_t       muOneDSS14Pi07000GeV[100];
	Double_t       muOneErrDSS14Pi07000GeV[100];
	Double_t       muTwoDSS14Pi07000GeV[100];
	Int_t       nlinesNLODSS14Pi07000GeV =       0;
	
	TString fileNameNLODSS14Pi07000GeV = "ExternalInput/Theory/pi0dss14-7000gev-dsigmadpt-eta060-pt40.dat";
	ifstream  fileNLODSS14Pi07000GeV;
	fileNLODSS14Pi07000GeV.open(fileNameNLODSS14Pi07000GeV,ios_base::in);
	cout << fileNameNLODSS14Pi07000GeV << endl;
	
	while(!fileNLODSS14Pi07000GeV.eof()){
		nlinesNLODSS14Pi07000GeV++;
		TString garbage;
		fileNLODSS14Pi07000GeV >> ptNLODSS14Pi07000GeV[nlinesNLODSS14Pi07000GeV] >> garbage >> muHalfDSS14Pi07000GeV[nlinesNLODSS14Pi07000GeV]  >> muTwoDSS14Pi07000GeV[nlinesNLODSS14Pi07000GeV]; 
		muOneDSS14Pi07000GeV[nlinesNLODSS14Pi07000GeV] = (muHalfDSS14Pi07000GeV[nlinesNLODSS14Pi07000GeV]+muTwoDSS14Pi07000GeV[nlinesNLODSS14Pi07000GeV])/ (2*2*TMath::Pi()*ptNLODSS14Pi07000GeV[nlinesNLODSS14Pi07000GeV]);
		muOneErrDSS14Pi07000GeV[nlinesNLODSS14Pi07000GeV] = (muHalfDSS14Pi07000GeV[nlinesNLODSS14Pi07000GeV]-muTwoDSS14Pi07000GeV[nlinesNLODSS14Pi07000GeV])/ (2*2*TMath::Pi()*ptNLODSS14Pi07000GeV[nlinesNLODSS14Pi07000GeV]);
		ptErrNLODSS14Pi07000GeV[nlinesNLODSS14Pi07000GeV] = 0.5;
		cout << nlinesNLODSS14Pi07000GeV << "         "  << ptNLODSS14Pi07000GeV[nlinesNLODSS14Pi07000GeV] << "         "  << muHalfDSS14Pi07000GeV[nlinesNLODSS14Pi07000GeV] << "         "  << muOneDSS14Pi07000GeV[nlinesNLODSS14Pi07000GeV] << "         "  << muTwoDSS14Pi07000GeV[nlinesNLODSS14Pi07000GeV] << endl;;
		
		
	}
	fileNLODSS14Pi07000GeV.close();
	TGraphAsymmErrors* graphNLOCalcDSS14InvSecPi07000GeV = new TGraphAsymmErrors(nlinesNLODSS14Pi07000GeV, ptNLODSS14Pi07000GeV, muOneDSS14Pi07000GeV, ptErrNLODSS14Pi07000GeV, ptErrNLODSS14Pi07000GeV,
																				muOneErrDSS14Pi07000GeV, muOneErrDSS14Pi07000GeV);
	TGraphAsymmErrors* graphNLOCalcDSS14InvYieldPi07000GeV = ScaleGraphAsym(graphNLOCalcDSS14InvSecPi07000GeV, 1/(xSection7000GeV*recalcBarn));
	
	Double_t       ptNLODSS14Pi02760GeV[100];
	Double_t       ptErrNLODSS14Pi02760GeV[100];
	Double_t       muHalfDSS14Pi02760GeV[100];
	Double_t       muOneDSS14Pi02760GeV[100];
	Double_t       muOneErrDSS14Pi02760GeV[100];
	Double_t       muTwoDSS14Pi02760GeV[100];
	Int_t       nlinesNLODSS14Pi02760GeV =       0;
	
	TString fileNameNLODSS14Pi02760GeV = "ExternalInput/Theory/pi0dss14-2760gev-dsigmadpt-eta060-pt40.dat";
	ifstream  fileNLODSS14Pi02760GeV;
	fileNLODSS14Pi02760GeV.open(fileNameNLODSS14Pi02760GeV,ios_base::in);
	cout << fileNameNLODSS14Pi02760GeV << endl;
	
	while(!fileNLODSS14Pi02760GeV.eof()){
		nlinesNLODSS14Pi02760GeV++;
		TString garbage;
		fileNLODSS14Pi02760GeV >> ptNLODSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] >> garbage >> muHalfDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV]  >> muTwoDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV]; 
		muOneDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] = (muHalfDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV]+muTwoDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV])/ (2*2*TMath::Pi()*ptNLODSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV]);
		muOneErrDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] = (muHalfDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV]-muTwoDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV])/(2*2*TMath::Pi()*ptNLODSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV]);
		ptErrNLODSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] = 0.5;
		cout << nlinesNLODSS14Pi02760GeV << "         "  << ptNLODSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] << "         "  << muHalfDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] << "         "  << muOneDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] << "         "  << muTwoDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] << endl;;
		
		
	}
	fileNLODSS14Pi02760GeV.close();
	TGraphAsymmErrors* graphNLOCalcDSS14InvSecPi02760GeV = new TGraphAsymmErrors(nlinesNLODSS14Pi02760GeV, ptNLODSS14Pi02760GeV, muOneDSS14Pi02760GeV, ptErrNLODSS14Pi02760GeV, ptErrNLODSS14Pi02760GeV,
																				muOneErrDSS14Pi02760GeV, muOneErrDSS14Pi02760GeV);
	TGraphAsymmErrors* graphNLOCalcDSS14InvYieldPi02760GeV = ScaleGraphAsym(graphNLOCalcDSS14InvSecPi02760GeV, 1/(xSection2760GeV*recalcBarn));


	TFile fileTheoryGraphsPP("ExternalInput/TheoryCompilationPP.root","RECREATE");

		graphNLOCalcInvSecPi0MuHalf8000GeV->Write("graphNLOCalcInvSecPi0MuHalf8000GeV");
		graphNLOCalcInvSecPi0MuOne8000GeV->Write("graphNLOCalcInvSecPi0MuOne8000GeV");
		graphNLOCalcInvSecPi0MuTwo8000GeV->Write("graphNLOCalcInvSecPi0MuTwo8000GeV");

		graphNLOCalcInvSecPi0MuHalf7000GeV->Write("graphNLOCalcInvSecPi0MuHalf7000GeV");
		graphNLOCalcInvSecPi0MuOne7000GeV->Write("graphNLOCalcInvSecPi0MuOne7000GeV");
		graphNLOCalcInvSecPi0MuTwo7000GeV->Write("graphNLOCalcInvSecPi0MuTwo7000GeV");

		graphNLOCalcInvSecPi0MuHalf2760GeV->Write("graphNLOCalcInvSecPi0MuHalf2760GeV");
		graphNLOCalcInvSecPi0MuOne2760GeV->Write("graphNLOCalcInvSecPi0MuOne2760GeV");
		graphNLOCalcInvSecPi0MuTwo2760GeV->Write("graphNLOCalcInvSecPi0MuTwo2760GeV");

		graphNLOCalcInvSecPi0MuHalf900GeV->Write("graphNLOCalcInvSecPi0MuHalf900GeV");
		graphNLOCalcInvSecPi0MuOne900GeV->Write("graphNLOCalcInvSecPi0MuOne900GeV");
		graphNLOCalcInvSecPi0MuTwo900GeV->Write("graphNLOCalcInvSecPi0MuTwo900GeV");

		graphNLOCalcBKKInvSecPi0MuTwo7000GeV->Write("graphNLOCalcBKKInvSecPi0MuTwo7000GeV");
		graphNLOCalcBKKInvSecPi0MuTwo900GeV->Write("graphNLOCalcBKKInvSecPi0MuTwo900GeV");

		graphNLOCalcDSSInvSecPi0MuTwo7000GeV->Write("graphNLOCalcDSSInvSecPi0MuTwo7000GeV");
		graphNLOCalcDSSInvSecPi0MuTwo2760GeV->Write("graphNLOCalcDSSInvSecPi0MuTwo2760GeV");
		graphNLOCalcDSSInvSecPi0MuTwo900GeV->Write("graphNLOCalcDSSInvSecPi0MuTwo900GeV");

		graphNLOCalcInvSecEtaMuHalf8000GeV->Write("graphNLOCalcInvSecEtaMuHalf8000GeV");
		graphNLOCalcInvSecEtaMuOne8000GeV->Write("graphNLOCalcInvSecEtaMuOne8000GeV");
		graphNLOCalcInvSecEtaMuTwo8000GeV->Write("graphNLOCalcInvSecEtaMuTwo8000GeV");
		
		graphNLOCalcInvSecEtaMuHalf7000GeV->Write("graphNLOCalcInvSecEtaMuHalf7000GeV");
		graphNLOCalcInvSecEtaMuOne7000GeV->Write("graphNLOCalcInvSecEtaMuOne7000GeV");
		graphNLOCalcInvSecEtaMuTwo7000GeV->Write("graphNLOCalcInvSecEtaMuTwo7000GeV");

		graphNLOCalcInvSecEtaMuHalf2760GeV->Write("graphNLOCalcInvSecEtaMuHalf2760GeV");
		graphNLOCalcInvSecEtaMuOne2760GeV->Write("graphNLOCalcInvSecEtaMuOne2760GeV");
		graphNLOCalcInvSecEtaMuTwo2760GeV->Write("graphNLOCalcInvSecEtaMuTwo2760GeV");

		graphNLOCalcInvSecEtaMuHalf900GeV->Write("graphNLOCalcInvSecEtaMuHalf900GeV");
		graphNLOCalcInvSecEtaMuOne900GeV->Write("graphNLOCalcInvSecEtaMuOne900GeV");
		graphNLOCalcInvSecEtaMuTwo900GeV->Write("graphNLOCalcInvSecEtaMuTwo900GeV");

		graphNLOCalcInvYieldPi0MuHalf8000GeV->Write("graphNLOCalcInvYieldPi0MuHalf8000GeV");
		graphNLOCalcInvYieldPi0MuOne8000GeV->Write("graphNLOCalcInvYieldPi0MuOne8000GeV");
		graphNLOCalcInvYieldPi0MuTwo8000GeV->Write("graphNLOCalcInvYieldPi0MuTwo8000GeV");

		graphNLOCalcInvYieldPi0MuHalf7000GeV->Write("graphNLOCalcInvYieldPi0MuHalf7000GeV");
		graphNLOCalcInvYieldPi0MuOne7000GeV->Write("graphNLOCalcInvYieldPi0MuOne7000GeV");
		graphNLOCalcInvYieldPi0MuTwo7000GeV->Write("graphNLOCalcInvYieldPi0MuTwo7000GeV");

		graphNLOCalcInvYieldPi0MuHalf2760GeV->Write("graphNLOCalcInvYieldPi0MuHalf2760GeV");
		graphNLOCalcInvYieldPi0MuOne2760GeV->Write("graphNLOCalcInvYieldPi0MuOne2760GeV");
		graphNLOCalcInvYieldPi0MuTwo2760GeV->Write("graphNLOCalcInvYieldPi0MuTwo2760GeV");

		graphNLOCalcInvYieldPi0MuHalf900GeV->Write("graphNLOCalcInvYieldPi0MuHalf900GeV");
		graphNLOCalcInvYieldPi0MuOne900GeV->Write("graphNLOCalcInvYieldPi0MuOne900GeV");
		graphNLOCalcInvYieldPi0MuTwo900GeV->Write("graphNLOCalcInvYieldPi0MuTwo900GeV");

		graphNLOCalcBKKInvYieldPi0MuTwo7000GeV->Write("graphNLOCalcBKKInvYieldPi0MuTwo7000GeV");
		graphNLOCalcBKKInvYieldPi0MuTwo900GeV->Write("graphNLOCalcBKKInvYieldPi0MuTwo900GeV");

		graphNLOCalcDSSInvYieldPi0MuTwo7000GeV->Write("graphNLOCalcDSSInvYieldPi0MuTwo7000GeV");
		graphNLOCalcDSSInvYieldPi0MuTwo2760GeV->Write("graphNLOCalcDSSInvYieldPi0MuTwo2760GeV");
		graphNLOCalcDSSInvYieldPi0MuTwo900GeV->Write("graphNLOCalcDSSInvYieldPi0MuTwo900GeV");

		graphNLOCalcInvYieldEtaMuHalf8000GeV->Write("graphNLOCalcInvYieldEtaMuHalf8000GeV");
		graphNLOCalcInvYieldEtaMuOne8000GeV->Write("graphNLOCalcInvYieldEtaMuOne8000GeV");
		graphNLOCalcInvYieldEtaMuTwo8000GeV->Write("graphNLOCalcInvYieldEtaMuTwo8000GeV");
		
		graphNLOCalcInvYieldEtaMuHalf7000GeV->Write("graphNLOCalcInvYieldEtaMuHalf7000GeV");
		graphNLOCalcInvYieldEtaMuOne7000GeV->Write("graphNLOCalcInvYieldEtaMuOne7000GeV");
		graphNLOCalcInvYieldEtaMuTwo7000GeV->Write("graphNLOCalcInvYieldEtaMuTwo7000GeV");

		graphNLOCalcInvYieldEtaMuHalf2760GeV->Write("graphNLOCalcInvYieldEtaMuHalf2760GeV");
		graphNLOCalcInvYieldEtaMuOne2760GeV->Write("graphNLOCalcInvYieldEtaMuOne2760GeV");
		graphNLOCalcInvYieldEtaMuTwo2760GeV->Write("graphNLOCalcInvYieldEtaMuTwo2760GeV");

		graphNLOCalcInvYieldEtaMuHalf900GeV->Write("graphNLOCalcInvYieldEtaMuHalf900GeV");
		graphNLOCalcInvYieldEtaMuOne900GeV->Write("graphNLOCalcInvYieldEtaMuOne900GeV");
		graphNLOCalcInvYieldEtaMuTwo900GeV->Write("graphNLOCalcInvYieldEtaMuTwo900GeV");

		graphEtaToPi0NLOMuHalf8TeV->Write("graphNLOCalcEtaOverPi0MuHalf8000GeV");
		graphEtaToPi0NLOMuOne8TeV->Write("graphNLOCalcEtaOverPi0MuOne8000GeV");
		graphEtaToPi0NLOMuTwo8TeV->Write("graphNLOCalcEtaOverPi0MuTwo8000GeV");
		
		graphEtaToPi0NLOMuHalf7TeV->Write("graphNLOCalcEtaOverPi0MuHalf7000GeV");
		graphEtaToPi0NLOMuOne7TeV->Write("graphNLOCalcEtaOverPi0MuOne7000GeV");
		graphEtaToPi0NLOMuTwo7TeV->Write("graphNLOCalcEtaOverPi0MuTwo7000GeV");

		graphEtaToPi0NLOMuHalf2760GeV->Write("graphNLOCalcEtaOverPi0MuHalf2760GeV");
		graphEtaToPi0NLOMuOne2760GeV->Write("graphNLOCalcEtaOverPi0MuOne2760GeV");
		graphEtaToPi0NLOMuTwo2760GeV->Write("graphNLOCalcEtaOverPi0MuTwo2760GeV");

		graphEtaToPi0NLOMuHalf900GeV->Write("graphNLOCalcEtaOverPi0MuHalf900GeV");
		graphEtaToPi0NLOMuOne900GeV->Write("graphNLOCalcEtaOverPi0MuOne900GeV");
		graphEtaToPi0NLOMuTwo900GeV->Write("graphNLOCalcEtaOverPi0MuTwo900GeV");

		histoPythia8Spec2760GeV->Write("histoInvSecPythia8Spec2760GeV");
		histoPythia8InvYield2760GeV->Write("histoInvYieldPythia8Spec2760GeV");
		histoPythia8Spec2760GeVVarBinning->Write("histoInvSecPythia8Spec2760GeVVarBinning");
		histoPythia8InvYield2760GeVVarBinning->Write("histoInvYieldPythia8Spec2760GeVVarBinning");
		
		graphInvYieldCGC2760GeV->Write("graphNLOCalcCGCInvYield2760GeV");
		graphInvYieldCGC7000GeV->Write("graphNLOCalcCGCInvYield7000GeV");
		graphInvCrossSecCGC2760GeV->Write("graphNLOCalcCGCInvCrossSec2760GeV");
		graphInvCrossSecCGC7000GeV->Write("graphNLOCalcCGCInvCrossSec7000GeV");
		graphInvCrossSecCGC7000GeV_mvgamma->Write("graphNLOCalcCGCInvCrossSec7000GeV_mvgamma");

		graphNLOCalcDSS14InvSecPi02760GeV->Write("graphNLOCalcDSS14InvCrossSec2760GeV");
		graphNLOCalcDSS14InvSecPi07000GeV->Write("graphNLOCalcDSS14InvCrossSec7000GeV");
		graphNLOCalcDSS14InvYieldPi02760GeV->Write("graphNLOCalcDSS14InvYield2760GeV");
		graphNLOCalcDSS14InvYieldPi07000GeV->Write("graphNLOCalcDSS14InvYield7000GeV");
		
		histoPi02760GeV_Pythia8_Tune4CMB->Write("histoInvYieldPi02760GeV_Pythia8_Tune4CMB");
		histoPi02760GeV_Pythia8_Tune4CptHard2->Write("histoInvYieldPi02760GeV_Pythia8_Tune4CptHard2");
		histoPi02760GeV_Pythia8_Tune4CptHard5->Write("histoInvYieldPi02760GeV_Pythia8_Tune4CptHard5");
		histoPi02760GeV_Pythia8_Tune4CptHard10->Write("histoInvYieldPi02760GeV_Pythia8_Tune4CptHard10");
		histoPi0UB2760GeV_Pythia8_Tune4CMB->Write("histoInvYieldPi02760GeV_Pythia8_Tune4CMB_UB");
		histoPi0UB2760GeV_Pythia8_Tune4CptHard2->Write("histoInvYieldPi02760GeV_Pythia8_Tune4CptHard2_UB");
		histoPi0UB2760GeV_Pythia8_Tune4CptHard5->Write("histoInvYieldPi02760GeV_Pythia8_Tune4CptHard5_UB");
		histoPi0UB2760GeV_Pythia8_Tune4CptHard10->Write("histoInvYieldPi02760GeV_Pythia8_Tune4CptHard10_UB");

		histoPi0FromK0s2760GeV_Pythia8_Tune4CMB->Write("histoInvYieldPi0FromK2760GeV_Pythia8_Tune4CMB");
		histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard2->Write("histoInvYieldPi0FromK2760GeV_Pythia8_Tune4CptHard2");
		histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard5->Write("histoInvYieldPi0FromK2760GeV_Pythia8_Tune4CptHard5");
		histoPi0FromK0s2760GeV_Pythia8_Tune4CptHard10->Write("histoInvYieldPi0FromK2760GeV_Pythia8_Tune4CptHard10");
		histoPi0FromK0sUB2760GeV_Pythia8_Tune4CMB->Write("histoInvYieldPi0FromK2760GeV_Pythia8_Tune4CMB_UB");
		histoPi0FromK0sUB2760GeV_Pythia8_Tune4CptHard2->Write("histoInvYieldPi0FromK2760GeV_Pythia8_Tune4CptHard2_UB");
		histoPi0FromK0sUB2760GeV_Pythia8_Tune4CptHard5->Write("histoInvYieldPi0FromK2760GeV_Pythia8_Tune4CptHard5_UB");
		histoPi0FromK0sUB2760GeV_Pythia8_Tune4CptHard10->Write("histoInvYieldPi0FromK2760GeV_Pythia8_Tune4CptHard10_UB");

		histoPi0FromEta2760GeV_Pythia8_Tune4CMB->Write("histoInvYieldPi0FromEta2760GeV_Pythia8_Tune4CMB");
		histoPi0FromEta2760GeV_Pythia8_Tune4CptHard2->Write("histoInvYieldPi0FromEta2760GeV_Pythia8_Tune4CptHard2");
		histoPi0FromEta2760GeV_Pythia8_Tune4CptHard5->Write("histoInvYieldPi0FromEta2760GeV_Pythia8_Tune4CptHard5");
		histoPi0FromEta2760GeV_Pythia8_Tune4CptHard10->Write("histoInvYieldPi0FromEta2760GeV_Pythia8_Tune4CptHard10");
		histoPi0FromEtaUB2760GeV_Pythia8_Tune4CMB->Write("histoInvYieldPi0FromEta2760GeV_Pythia8_Tune4CMB_UB");
		histoPi0FromEtaUB2760GeV_Pythia8_Tune4CptHard2->Write("histoInvYieldPi0FromEta2760GeV_Pythia8_Tune4CptHard2_UB");
		histoPi0FromEtaUB2760GeV_Pythia8_Tune4CptHard5->Write("histoInvYieldPi0FromEta2760GeV_Pythia8_Tune4CptHard5_UB");
		histoPi0FromEtaUB2760GeV_Pythia8_Tune4CptHard10->Write("histoInvYieldPi0FromEta2760GeV_Pythia8_Tune4CptHard10_UB");

		histoPi0FromLambda2760GeV_Pythia8_Tune4CMB->Write("histoInvYieldPi0FromLambda2760GeV_Pythia8_Tune4CMB");
		histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard2->Write("histoInvYieldPi0FromLambda2760GeV_Pythia8_Tune4CptHard2");
		histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard5->Write("histoInvYieldPi0FromLambda2760GeV_Pythia8_Tune4CptHard5");
		histoPi0FromLambda2760GeV_Pythia8_Tune4CptHard10->Write("histoInvYieldPi0FromLambda2760GeV_Pythia8_Tune4CptHard10");
		histoPi0FromLambdaUB2760GeV_Pythia8_Tune4CMB->Write("histoInvYieldPi0FromLambda2760GeV_Pythia8_Tune4CMB_UB");
		histoPi0FromLambdaUB2760GeV_Pythia8_Tune4CptHard2->Write("histoInvYieldPi0FromLambda2760GeV_Pythia8_Tune4CptHard2_UB");
		histoPi0FromLambdaUB2760GeV_Pythia8_Tune4CptHard5->Write("histoInvYieldPi0FromLambda2760GeV_Pythia8_Tune4CptHard5_UB");
		histoPi0FromLambdaUB2760GeV_Pythia8_Tune4CptHard10->Write("histoInvYieldPi0FromLambda2760GeV_Pythia8_Tune4CptHard10_UB");

		histoEta2760GeV_Pythia8_Tune4CMB->Write("histoInvYieldEta2760GeV_Pythia8_Tune4CMB");
		histoEta2760GeV_Pythia8_Tune4CptHard2->Write("histoInvYieldEta2760GeV_Pythia8_Tune4CptHard2");
		histoEta2760GeV_Pythia8_Tune4CptHard5->Write("histoInvYieldEta2760GeV_Pythia8_Tune4CptHard5");
		histoEta2760GeV_Pythia8_Tune4CptHard10->Write("histoInvYieldEta2760GeV_Pythia8_Tune4CptHard10");
		histoEtaUB2760GeV_Pythia8_Tune4CMB->Write("histoInvYieldEta2760GeV_Pythia8_Tune4CMB_UB");
		histoEtaUB2760GeV_Pythia8_Tune4CptHard2->Write("histoInvYieldEta2760GeV_Pythia8_Tune4CptHard2_UB");
		histoEtaUB2760GeV_Pythia8_Tune4CptHard5->Write("histoInvYieldEta2760GeV_Pythia8_Tune4CptHard5_UB");
		histoEtaUB2760GeV_Pythia8_Tune4CptHard10->Write("histoInvYieldEta2760GeV_Pythia8_Tune4CptHard10_UB");

		histoPiPl2760GeV_Pythia8_Tune4CMB->Write("histoInvYieldPiPl2760GeV_Pythia8_Tune4CMB");
		histoPiPl2760GeV_Pythia8_Tune4CptHard2->Write("histoInvYieldPiPl2760GeV_Pythia8_Tune4CptHard2");
		histoPiPl2760GeV_Pythia8_Tune4CptHard5->Write("histoInvYieldPiPl2760GeV_Pythia8_Tune4CptHard5");
		histoPiPl2760GeV_Pythia8_Tune4CptHard10->Write("histoInvYieldPiPl2760GeV_Pythia8_Tune4CptHard10");
		histoPiPlUB2760GeV_Pythia8_Tune4CMB->Write("histoInvYieldPiPl2760GeV_Pythia8_Tune4CMB_UB");
		histoPiPlUB2760GeV_Pythia8_Tune4CptHard2->Write("histoInvYieldPiPl2760GeV_Pythia8_Tune4CptHard2_UB");
		histoPiPlUB2760GeV_Pythia8_Tune4CptHard5->Write("histoInvYieldPiPl2760GeV_Pythia8_Tune4CptHard5_UB");
		histoPiPlUB2760GeV_Pythia8_Tune4CptHard10->Write("histoInvYieldPiPl2760GeV_Pythia8_Tune4CptHard10_UB");

		histoPiMi2760GeV_Pythia8_Tune4CMB->Write("histoInvYieldPiMi2760GeV_Pythia8_Tune4CMB");
		histoPiMi2760GeV_Pythia8_Tune4CptHard2->Write("histoInvYieldPiMi2760GeV_Pythia8_Tune4CptHard2");
		histoPiMi2760GeV_Pythia8_Tune4CptHard5->Write("histoInvYieldPiMi2760GeV_Pythia8_Tune4CptHard5");
		histoPiMi2760GeV_Pythia8_Tune4CptHard10->Write("histoInvYieldPiMi2760GeV_Pythia8_Tune4CptHard10");
		histoPiMiUB2760GeV_Pythia8_Tune4CMB->Write("histoInvYieldPiMi2760GeV_Pythia8_Tune4CMB_UB");
		histoPiMiUB2760GeV_Pythia8_Tune4CptHard2->Write("histoInvYieldPiMi2760GeV_Pythia8_Tune4CptHard2_UB");
		histoPiMiUB2760GeV_Pythia8_Tune4CptHard5->Write("histoInvYieldPiMi2760GeV_Pythia8_Tune4CptHard5_UB");
		histoPiMiUB2760GeV_Pythia8_Tune4CptHard10->Write("histoInvYieldPiMi2760GeV_Pythia8_Tune4CptHard10_UB");
		
		ratioPi02760GeV_Pythia8_Tune4CpTHard2DivMB->Write();
		ratioPi02760GeV_Pythia8_Tune4CpTHard5DivpTHard2->Write();
		ratioPi02760GeV_Pythia8_Tune4CpTHard10DivpTHard5->Write();
		ratioPi0UB2760GeV_Pythia8_Tune4CpTHard2DivMB->Write();
		ratioPi0UB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2->Write();
		ratioPi0UB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5->Write();

		ratioPiPl2760GeV_Pythia8_Tune4CpTHard2DivMB->Write();
		ratioPiPl2760GeV_Pythia8_Tune4CpTHard5DivpTHard2->Write();
		ratioPiPl2760GeV_Pythia8_Tune4CpTHard10DivpTHard5->Write();
		ratioPiPlUB2760GeV_Pythia8_Tune4CpTHard2DivMB->Write();
		ratioPiPlUB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2->Write();
		ratioPiPlUB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5->Write();

		ratioPiMi2760GeV_Pythia8_Tune4CpTHard2DivMB->Write();
		ratioPiMi2760GeV_Pythia8_Tune4CpTHard5DivpTHard2->Write();
		ratioPiMi2760GeV_Pythia8_Tune4CpTHard10DivpTHard5->Write();
		ratioPiMiUB2760GeV_Pythia8_Tune4CpTHard2DivMB->Write();
		ratioPiMiUB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2->Write();
		ratioPiMiUB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5->Write();

		ratioEta2760GeV_Pythia8_Tune4CpTHard2DivMB->Write();
		ratioEta2760GeV_Pythia8_Tune4CpTHard5DivpTHard2->Write();
		ratioEta2760GeV_Pythia8_Tune4CpTHard10DivpTHard5->Write();
		ratioEtaUB2760GeV_Pythia8_Tune4CpTHard2DivMB->Write();
		ratioEtaUB2760GeV_Pythia8_Tune4CpTHard5DivpTHard2->Write();
		ratioEtaUB2760GeV_Pythia8_Tune4CpTHard10DivpTHard5->Write();

		histoPi02760GeV_Pythia8_TuneMonachMB->Write("histoInvYieldPi02760GeV_Pythia8_TuneMonachMB");
		histoPi02760GeV_Pythia8_TuneMonachptHard2->Write("histoInvYieldPi02760GeV_Pythia8_TuneMonachptHard2");
		histoPi02760GeV_Pythia8_TuneMonachptHard5->Write("histoInvYieldPi02760GeV_Pythia8_TuneMonachptHard5");
		histoPi02760GeV_Pythia8_TuneMonachptHard10->Write("histoInvYieldPi02760GeV_Pythia8_TuneMonachptHard10");
		histoPi0UB2760GeV_Pythia8_TuneMonachMB->Write("histoInvYieldPi02760GeV_Pythia8_TuneMonachMB_UB");
		histoPi0UB2760GeV_Pythia8_TuneMonachptHard2->Write("histoInvYieldPi02760GeV_Pythia8_TuneMonachptHard2_UB");
		histoPi0UB2760GeV_Pythia8_TuneMonachptHard5->Write("histoInvYieldPi02760GeV_Pythia8_TuneMonachptHard5_UB");
		histoPi0UB2760GeV_Pythia8_TuneMonachptHard10->Write("histoInvYieldPi02760GeV_Pythia8_TuneMonachptHard10_UB");

		histoPi0FromK0s2760GeV_Pythia8_TuneMonachMB->Write("histoInvYieldPi0FromK2760GeV_Pythia8_TuneMonachMB");
		histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard2->Write("histoInvYieldPi0FromK2760GeV_Pythia8_TuneMonachptHard2");
		histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard5->Write("histoInvYieldPi0FromK2760GeV_Pythia8_TuneMonachptHard5");
		histoPi0FromK0s2760GeV_Pythia8_TuneMonachptHard10->Write("histoInvYieldPi0FromK2760GeV_Pythia8_TuneMonachptHard10");
		histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachMB->Write("histoInvYieldPi0FromK2760GeV_Pythia8_TuneMonachMB_UB");
		histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachptHard2->Write("histoInvYieldPi0FromK2760GeV_Pythia8_TuneMonachptHard2_UB");
		histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachptHard5->Write("histoInvYieldPi0FromK2760GeV_Pythia8_TuneMonachptHard5_UB");
		histoPi0FromK0sUB2760GeV_Pythia8_TuneMonachptHard10->Write("histoInvYieldPi0FromK2760GeV_Pythia8_TuneMonachptHard10_UB");

		histoPi0FromEta2760GeV_Pythia8_TuneMonachMB->Write("histoInvYieldPi0FromEta2760GeV_Pythia8_TuneMonachMB");
		histoPi0FromEta2760GeV_Pythia8_TuneMonachptHard2->Write("histoInvYieldPi0FromEta2760GeV_Pythia8_TuneMonachptHard2");
		histoPi0FromEta2760GeV_Pythia8_TuneMonachptHard5->Write("histoInvYieldPi0FromEta2760GeV_Pythia8_TuneMonachptHard5");
		histoPi0FromEta2760GeV_Pythia8_TuneMonachptHard10->Write("histoInvYieldPi0FromEta2760GeV_Pythia8_TuneMonachptHard10");
		histoPi0FromEtaUB2760GeV_Pythia8_TuneMonachMB->Write("histoInvYieldPi0FromEta2760GeV_Pythia8_TuneMonachMB_UB");
		histoPi0FromEtaUB2760GeV_Pythia8_TuneMonachptHard2->Write("histoInvYieldPi0FromEta2760GeV_Pythia8_TuneMonachptHard2_UB");
		histoPi0FromEtaUB2760GeV_Pythia8_TuneMonachptHard5->Write("histoInvYieldPi0FromEta2760GeV_Pythia8_TuneMonachptHard5_UB");
		histoPi0FromEtaUB2760GeV_Pythia8_TuneMonachptHard10->Write("histoInvYieldPi0FromEta2760GeV_Pythia8_TuneMonachptHard10_UB");

		histoPi0FromLambda2760GeV_Pythia8_TuneMonachMB->Write("histoInvYieldPi0FromLambda2760GeV_Pythia8_TuneMonachMB");
		histoPi0FromLambda2760GeV_Pythia8_TuneMonachptHard2->Write("histoInvYieldPi0FromLambda2760GeV_Pythia8_TuneMonachptHard2");
		histoPi0FromLambda2760GeV_Pythia8_TuneMonachptHard5->Write("histoInvYieldPi0FromLambda2760GeV_Pythia8_TuneMonachptHard5");
		histoPi0FromLambda2760GeV_Pythia8_TuneMonachptHard10->Write("histoInvYieldPi0FromLambda2760GeV_Pythia8_TuneMonachptHard10");
		histoPi0FromLambdaUB2760GeV_Pythia8_TuneMonachMB->Write("histoInvYieldPi0FromLambda2760GeV_Pythia8_TuneMonachMB_UB");
		histoPi0FromLambdaUB2760GeV_Pythia8_TuneMonachptHard2->Write("histoInvYieldPi0FromLambda2760GeV_Pythia8_TuneMonachptHard2_UB");
		histoPi0FromLambdaUB2760GeV_Pythia8_TuneMonachptHard5->Write("histoInvYieldPi0FromLambda2760GeV_Pythia8_TuneMonachptHard5_UB");
		histoPi0FromLambdaUB2760GeV_Pythia8_TuneMonachptHard10->Write("histoInvYieldPi0FromLambda2760GeV_Pythia8_TuneMonachptHard10_UB");

		histoEta2760GeV_Pythia8_TuneMonachMB->Write("histoInvYieldEta2760GeV_Pythia8_TuneMonachMB");
		histoEta2760GeV_Pythia8_TuneMonachptHard2->Write("histoInvYieldEta2760GeV_Pythia8_TuneMonachptHard2");
		histoEta2760GeV_Pythia8_TuneMonachptHard5->Write("histoInvYieldEta2760GeV_Pythia8_TuneMonachptHard5");
		histoEta2760GeV_Pythia8_TuneMonachptHard10->Write("histoInvYieldEta2760GeV_Pythia8_TuneMonachptHard10");
		histoEtaUB2760GeV_Pythia8_TuneMonachMB->Write("histoInvYieldEta2760GeV_Pythia8_TuneMonachMB_UB");
		histoEtaUB2760GeV_Pythia8_TuneMonachptHard2->Write("histoInvYieldEta2760GeV_Pythia8_TuneMonachptHard2_UB");
		histoEtaUB2760GeV_Pythia8_TuneMonachptHard5->Write("histoInvYieldEta2760GeV_Pythia8_TuneMonachptHard5_UB");
		histoEtaUB2760GeV_Pythia8_TuneMonachptHard10->Write("histoInvYieldEta2760GeV_Pythia8_TuneMonachptHard10_UB");

		histoPiPl2760GeV_Pythia8_TuneMonachMB->Write("histoInvYieldPiPl2760GeV_Pythia8_TuneMonachMB");
		histoPiPl2760GeV_Pythia8_TuneMonachptHard2->Write("histoInvYieldPiPl2760GeV_Pythia8_TuneMonachptHard2");
		histoPiPl2760GeV_Pythia8_TuneMonachptHard5->Write("histoInvYieldPiPl2760GeV_Pythia8_TuneMonachptHard5");
		histoPiPl2760GeV_Pythia8_TuneMonachptHard10->Write("histoInvYieldPiPl2760GeV_Pythia8_TuneMonachptHard10");
		histoPiPlUB2760GeV_Pythia8_TuneMonachMB->Write("histoInvYieldPiPl2760GeV_Pythia8_TuneMonachMB_UB");
		histoPiPlUB2760GeV_Pythia8_TuneMonachptHard2->Write("histoInvYieldPiPl2760GeV_Pythia8_TuneMonachptHard2_UB");
		histoPiPlUB2760GeV_Pythia8_TuneMonachptHard5->Write("histoInvYieldPiPl2760GeV_Pythia8_TuneMonachptHard5_UB");
		histoPiPlUB2760GeV_Pythia8_TuneMonachptHard10->Write("histoInvYieldPiPl2760GeV_Pythia8_TuneMonachptHard10_UB");

		histoPiMi2760GeV_Pythia8_TuneMonachMB->Write("histoInvYieldPiMi2760GeV_Pythia8_TuneMonachMB");
		histoPiMi2760GeV_Pythia8_TuneMonachptHard2->Write("histoInvYieldPiMi2760GeV_Pythia8_TuneMonachptHard2");
		histoPiMi2760GeV_Pythia8_TuneMonachptHard5->Write("histoInvYieldPiMi2760GeV_Pythia8_TuneMonachptHard5");
		histoPiMi2760GeV_Pythia8_TuneMonachptHard10->Write("histoInvYieldPiMi2760GeV_Pythia8_TuneMonachptHard10");
		histoPiMiUB2760GeV_Pythia8_TuneMonachMB->Write("histoInvYieldPiMi2760GeV_Pythia8_TuneMonachMB_UB");
		histoPiMiUB2760GeV_Pythia8_TuneMonachptHard2->Write("histoInvYieldPiMi2760GeV_Pythia8_TuneMonachptHard2_UB");
		histoPiMiUB2760GeV_Pythia8_TuneMonachptHard5->Write("histoInvYieldPiMi2760GeV_Pythia8_TuneMonachptHard5_UB");
		histoPiMiUB2760GeV_Pythia8_TuneMonachptHard10->Write("histoInvYieldPiMi2760GeV_Pythia8_TuneMonachptHard10_UB");
		
		ratioPi02760GeV_Pythia8_TuneMonachpTHard2DivMB->Write();
		ratioPi02760GeV_Pythia8_TuneMonachpTHard5DivpTHard2->Write();
		ratioPi02760GeV_Pythia8_TuneMonachpTHard10DivpTHard5->Write();
		ratioPi0UB2760GeV_Pythia8_TuneMonachpTHard2DivMB->Write();
		ratioPi0UB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2->Write();
		ratioPi0UB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5->Write();

		ratioPiPl2760GeV_Pythia8_TuneMonachpTHard2DivMB->Write();
		ratioPiPl2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2->Write();
		ratioPiPl2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5->Write();
		ratioPiPlUB2760GeV_Pythia8_TuneMonachpTHard2DivMB->Write();
		ratioPiPlUB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2->Write();
		ratioPiPlUB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5->Write();

		ratioPiMi2760GeV_Pythia8_TuneMonachpTHard2DivMB->Write();
		ratioPiMi2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2->Write();
		ratioPiMi2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5->Write();
		ratioPiMiUB2760GeV_Pythia8_TuneMonachpTHard2DivMB->Write();
		ratioPiMiUB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2->Write();
		ratioPiMiUB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5->Write();

		ratioEta2760GeV_Pythia8_TuneMonachpTHard2DivMB->Write();
		ratioEta2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2->Write();
		ratioEta2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5->Write();
		ratioEtaUB2760GeV_Pythia8_TuneMonachpTHard2DivMB->Write();
		ratioEtaUB2760GeV_Pythia8_TuneMonachpTHard5DivpTHard2->Write();
		ratioEtaUB2760GeV_Pythia8_TuneMonachpTHard10DivpTHard5->Write();
		
	fileTheoryGraphsPP.Close();

}







