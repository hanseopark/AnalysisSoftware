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

void ProduceTheoryGraphsPbPb(TString specifier = ""){    
    
	StyleSettingsThesis();    
	SetPlotStyle();

	//********************************************************************************************************************************************        
    //***************************************************** Cracow model for LHC11h yields *******************************************************
	// PRC 90, 014906 (2014)
    // columns in the file represent:  nBin p_{T} [GeV/c]  p_{T} [GeV/c]_[MIN]   p_{T} [GeV/c]_[MAX]    dN/(2 #pi p_{T} dp_{T} dy)    dN/(2 #pi p_{T} dp_{T} dy)_[ERROR]
	
	const Int_t MaxNPtBins = 15; //instead of 25 for the full pt range, the theory guys asked us to stop at 3 Gev/c
	Int_t nbins;
	Double_t BinLowerEdge, BinUpperEdge;

	Double_t ptPi0LowPtNonEq_0010[15];    
	Double_t yieldPi0LowPtNonEq_0010[15];    
	Double_t yieldErrPi0LowPtNonEq_0010[15];    
	Double_t Pi0LowErr_0010[15];
	Double_t Pi0HighErr_0010[15];

	TString fileNamePi0LowPtNonEq_0010 = "ExternalInputPbPb/Theory/CracowModel/lowPt-chem-non-equilibrium_15Oct2015/pion0-10.txt";
	ifstream  filePi0LowPtNonEq_0010;
	filePi0LowPtNonEq_0010.open(fileNamePi0LowPtNonEq_0010,ios_base::in);
	cout << fileNamePi0LowPtNonEq_0010 << endl;

	Int_t nlinesLowPt = 0;
	while(!filePi0LowPtNonEq_0010.eof() && nlinesLowPt < MaxNPtBins){
		filePi0LowPtNonEq_0010 >> nbins >> ptPi0LowPtNonEq_0010[nlinesLowPt] >> BinLowerEdge >> BinUpperEdge >> yieldPi0LowPtNonEq_0010[nlinesLowPt] >> yieldErrPi0LowPtNonEq_0010[nlinesLowPt];
		cout << nlinesLowPt << "\t "  << ptPi0LowPtNonEq_0010[nlinesLowPt] << "\t "  << yieldPi0LowPtNonEq_0010[nlinesLowPt] << " +- "  << yieldErrPi0LowPtNonEq_0010[nlinesLowPt] << endl;

		// error of the yields is not taken (is put to zero) because they have trouble with the centrality range and it is not a "real" error (just the cent. bin width)
		Pi0LowErr_0010[nlinesLowPt] = 0;
		Pi0HighErr_0010[nlinesLowPt] = 0;
		nlinesLowPt++;
	}
	filePi0LowPtNonEq_0010.close();
    
    Double_t ptEtaLowPtNonEq_0010[15];    
    Double_t yieldEtaLowPtNonEq_0010[15];    
    Double_t yieldErrEtaLowPtNonEq_0010[15];    
    Double_t EtaLowErr_0010[15];
    Double_t EtaHighErr_0010[15];

    Double_t ratioEtaPi0LowPtNonEq_0010[15];    

    TString fileNameEtaLowPtNonEq_0010 = "ExternalInputPbPb/Theory/CracowModel/lowPt-chem-non-equilibrium_15Oct2015/eta0-10.txt";
    ifstream  fileEtaLowPtNonEq_0010;
    fileEtaLowPtNonEq_0010.open(fileNameEtaLowPtNonEq_0010,ios_base::in);
    cout << fileNameEtaLowPtNonEq_0010 << endl;
    
    nlinesLowPt = 0;
    while(!fileEtaLowPtNonEq_0010.eof() && nlinesLowPt < MaxNPtBins){
		fileEtaLowPtNonEq_0010 >> nbins >> ptEtaLowPtNonEq_0010[nlinesLowPt] >> BinLowerEdge >> BinUpperEdge >> yieldEtaLowPtNonEq_0010[nlinesLowPt] >> yieldErrEtaLowPtNonEq_0010[nlinesLowPt];
		cout << nlinesLowPt << "\t "  << ptEtaLowPtNonEq_0010[nlinesLowPt] << "\t "  << yieldEtaLowPtNonEq_0010[nlinesLowPt] << " +- "  << yieldErrEtaLowPtNonEq_0010[nlinesLowPt] << endl;
		
		// error of the yields is not taken (is put to zero) because they have trouble with the centrality range and it is not a "real" error (just the cent. bin width)
		EtaLowErr_0010[nlinesLowPt] = 0;
		EtaHighErr_0010[nlinesLowPt] = 0;
        
        ratioEtaPi0LowPtNonEq_0010[nlinesLowPt] = yieldEtaLowPtNonEq_0010[nlinesLowPt]/yieldPi0LowPtNonEq_0010[nlinesLowPt];
		cout << "Eta/Pi0 ratio: " << endl;
        cout << nlinesLowPt << "\t "  << ptEtaLowPtNonEq_0010[nlinesLowPt] << "\t "  << ratioEtaPi0LowPtNonEq_0010[nlinesLowPt] << endl;
		
		nlinesLowPt++;
    }
    fileEtaLowPtNonEq_0010.close();
    
	
	Double_t ptPi0LowPtNonEq_2050[15];    
	Double_t yieldPi0LowPtNonEq_2050[15];    
	Double_t yieldErrPi0LowPtNonEq_2050[15];    
	Double_t Pi0LowErr_2050[15];
	Double_t Pi0HighErr_2050[15];

	TString fileNamePi0LowPtNonEq_2050 = "ExternalInputPbPb/Theory/CracowModel/lowPt-chem-non-equilibrium_15Oct2015/pion20-50.txt";
	ifstream  filePi0LowPtNonEq_2050;
	filePi0LowPtNonEq_2050.open(fileNamePi0LowPtNonEq_2050,ios_base::in);
	cout << fileNamePi0LowPtNonEq_2050 << endl;

	nlinesLowPt = 0;
	while(!filePi0LowPtNonEq_2050.eof() && nlinesLowPt < MaxNPtBins){
		filePi0LowPtNonEq_2050 >> nbins >> ptPi0LowPtNonEq_2050[nlinesLowPt] >> BinLowerEdge >> BinUpperEdge >> yieldPi0LowPtNonEq_2050[nlinesLowPt] >> yieldErrPi0LowPtNonEq_2050[nlinesLowPt];
		cout << nlinesLowPt << "\t "  << ptPi0LowPtNonEq_2050[nlinesLowPt] << "\t "  << yieldPi0LowPtNonEq_2050[nlinesLowPt] << " +- "  << yieldErrPi0LowPtNonEq_2050[nlinesLowPt] << endl;

		// error of the yields is not taken (is put to zero) because they have trouble with the centrality range and it is not a "real" error (just the cent. bin width)
		Pi0LowErr_2050[nlinesLowPt] = 0;
		Pi0HighErr_2050[nlinesLowPt] = 0;
		nlinesLowPt++;
	}
	filePi0LowPtNonEq_2050.close();
    
	
    Double_t ptEtaLowPtNonEq_2050[15];    
    Double_t yieldEtaLowPtNonEq_2050[15];    
    Double_t yieldErrEtaLowPtNonEq_2050[15];    
    Double_t EtaLowErr_2050[15];
    Double_t EtaHighErr_2050[15];

    Double_t ratioEtaPi0LowPtNonEq_2050[15];    

    TString fileNameEtaLowPtNonEq_2050 = "ExternalInputPbPb/Theory/CracowModel/lowPt-chem-non-equilibrium_15Oct2015/eta20-50.txt";
    ifstream  fileEtaLowPtNonEq_2050;
    fileEtaLowPtNonEq_2050.open(fileNameEtaLowPtNonEq_2050,ios_base::in);
    cout << fileNameEtaLowPtNonEq_2050 << endl;
    
    nlinesLowPt = 0;
    while(!fileEtaLowPtNonEq_2050.eof() && nlinesLowPt < MaxNPtBins){
        fileEtaLowPtNonEq_2050 >> nbins >> ptEtaLowPtNonEq_2050[nlinesLowPt] >> BinLowerEdge >> BinUpperEdge >> yieldEtaLowPtNonEq_2050[nlinesLowPt] >> yieldErrEtaLowPtNonEq_2050[nlinesLowPt];
        cout << nlinesLowPt << "\t "  << ptEtaLowPtNonEq_2050[nlinesLowPt] << "\t "  << yieldEtaLowPtNonEq_2050[nlinesLowPt] << " +- "  << yieldErrEtaLowPtNonEq_2050[nlinesLowPt] << endl;
		
		// error of the yields is not taken (is put to zero) because they have trouble with the centrality range and it is not a "real" error (just the cent. bin width)
		EtaLowErr_2050[nlinesLowPt] = 0;
		EtaHighErr_2050[nlinesLowPt] = 0;
        
        ratioEtaPi0LowPtNonEq_2050[nlinesLowPt] = yieldEtaLowPtNonEq_2050[nlinesLowPt]/yieldPi0LowPtNonEq_2050[nlinesLowPt];
		cout << "Eta/Pi0 ratio: " << endl;
        cout << nlinesLowPt << "\t "  << ptEtaLowPtNonEq_2050[nlinesLowPt] << "\t "  << ratioEtaPi0LowPtNonEq_2050[nlinesLowPt] << endl;
        nlinesLowPt++;
    }
    fileEtaLowPtNonEq_2050.close();

	
    Double_t ptChargedLowPtNonEq_0010[15];    
    Double_t yieldChargedPionLowPtNonEq_0010[15];    
    Double_t yieldChargedKaonLowPtNonEq_0010[15];    
    Double_t ratioKaonsToPionsLowPtNonEq_0010[15];    
    Double_t ChargedPionLowErr_0010[15];
    Double_t ChargedPionHighErr_0010[15];

    Double_t yieldNeutralToChargedPionLowPtNonEq_0010[15];    
    
    TString fileNameChargedPionsKaonsLowPtNonEq_0010 = "ExternalInputPbPb/Theory/CracowModel/lowPt-chem-non-equilibrium_15Oct2015/KaonsAndPions0-10.txt";
//     #dN / ( 2 \pi p_T dp_T dy) for pions (\pi^+ + \pi^- ) / 2 and Kaons (K^+ + K^-) / 2, 
//     and the K/pi ratio for 0-10% centrality as the function of p_T, which is (p_T^{MIN} + p_T^{MAX} ) / 2 
//     for the parameters from http://arxiv.org/abs/1503.04040
    ifstream  fileChargedPionsKaonsLowPtNonEq_0010;
    fileChargedPionsKaonsLowPtNonEq_0010.open(fileNameChargedPionsKaonsLowPtNonEq_0010,ios_base::in);
    cout << fileNameChargedPionsKaonsLowPtNonEq_0010 << endl;
    
    nlinesLowPt = 0;
    while(!fileChargedPionsKaonsLowPtNonEq_0010.eof() && nlinesLowPt < MaxNPtBins){
        fileChargedPionsKaonsLowPtNonEq_0010 >> ptChargedLowPtNonEq_0010[nlinesLowPt] >> yieldChargedPionLowPtNonEq_0010[nlinesLowPt] >> yieldChargedKaonLowPtNonEq_0010[nlinesLowPt] >> ratioKaonsToPionsLowPtNonEq_0010[nlinesLowPt]; 
        cout << nlinesLowPt << "\t " << ptChargedLowPtNonEq_0010[nlinesLowPt] << "\t pions: "  << yieldChargedPionLowPtNonEq_0010[nlinesLowPt] << "\t kaons:"  << yieldChargedKaonLowPtNonEq_0010[nlinesLowPt] << "\t kaons/pions:"  << ratioKaonsToPionsLowPtNonEq_0010[nlinesLowPt] << endl;
        
        yieldNeutralToChargedPionLowPtNonEq_0010[nlinesLowPt] = yieldPi0LowPtNonEq_0010[nlinesLowPt]/yieldChargedPionLowPtNonEq_0010[nlinesLowPt];
        ChargedPionLowErr_0010[nlinesLowPt] = 0.;
        ChargedPionHighErr_0010[nlinesLowPt] = 0.;
        nlinesLowPt++;
    }
    fileChargedPionsKaonsLowPtNonEq_0010.close();
    
    TGraphAsymmErrors *TheoryCracowPi0LowPt_0010 = new TGraphAsymmErrors(MaxNPtBins, ptPi0LowPtNonEq_0010, yieldPi0LowPtNonEq_0010, 0, 0, Pi0LowErr_0010, Pi0HighErr_0010);
    TGraphAsymmErrors *TheoryCracowEtaLowPt_0010 = new TGraphAsymmErrors(MaxNPtBins, ptEtaLowPtNonEq_0010, yieldEtaLowPtNonEq_0010, 0, 0, EtaLowErr_0010, EtaHighErr_0010);
    TGraphAsymmErrors *TheoryCracowEtaToPi0LowPt_0010 = new TGraphAsymmErrors(MaxNPtBins, ptEtaLowPtNonEq_0010, ratioEtaPi0LowPtNonEq_0010, 0, 0, 0, 0);

	TGraphAsymmErrors *TheoryCracowPi0LowPt_2050 = new TGraphAsymmErrors(MaxNPtBins, ptPi0LowPtNonEq_2050, yieldPi0LowPtNonEq_2050, 0, 0, Pi0LowErr_2050, Pi0HighErr_2050);
    TGraphAsymmErrors *TheoryCracowEtaLowPt_2050 = new TGraphAsymmErrors(MaxNPtBins, ptEtaLowPtNonEq_2050, yieldEtaLowPtNonEq_2050, 0, 0, EtaLowErr_2050, EtaHighErr_2050);
    TGraphAsymmErrors *TheoryCracowEtaToPi0LowPt_2050 = new TGraphAsymmErrors(MaxNPtBins, ptEtaLowPtNonEq_2050, ratioEtaPi0LowPtNonEq_2050, 0, 0, 0, 0);

    TGraphAsymmErrors *TheoryCracowChargedPionLowPt_0010 = new TGraphAsymmErrors(MaxNPtBins, ptChargedLowPtNonEq_0010, yieldChargedPionLowPtNonEq_0010, 0, 0, ChargedPionLowErr_0010, ChargedPionHighErr_0010);
    TGraphAsymmErrors *TheoryCracowChargedKaonLowPt_0010 = new TGraphAsymmErrors(MaxNPtBins, ptChargedLowPtNonEq_0010, yieldChargedKaonLowPtNonEq_0010, 0, 0, ChargedPionLowErr_0010, ChargedPionHighErr_0010);

	TGraphAsymmErrors *TheoryCracowKaonsToPionsLowPt_0010 = new TGraphAsymmErrors(MaxNPtBins, ptChargedLowPtNonEq_0010, ratioKaonsToPionsLowPtNonEq_0010, 0, 0, ChargedPionLowErr_0010, ChargedPionHighErr_0010);
    TGraphAsymmErrors *TheoryCracowNeutralToChargedPionLowPt_0010 = new TGraphAsymmErrors(MaxNPtBins, ptChargedLowPtNonEq_0010, yieldNeutralToChargedPionLowPtNonEq_0010, 0, 0, 0, 0);

	
	//********************************************************************************************************************************************    
    //*************************************** Raa theory - Djordjevic pred for Pi0 2011  *********************************************************
    // citing M Djordjevic, M. Djordjevic and B. Blagojevic, Phys. Lett. B 737 (2014) 298-302
    // and M Djordjevic and M. Djordjevic, Phys. Lett. B 734 (2014) 286-289
    // prediction for 0-10% for LHC11h data
    
    Double_t ptPi0Djordjevic[27];
    
    Double_t pi0RaaDjordjevic_0010[27];
    Double_t pi0RaaLowEdge_0010[27];
    Double_t pi0RaaHighEdge_0010[27];
    Double_t pi0RaaLowEdgeError_0010[27];
    Double_t pi0RaaHighEdgeError_0010[27];
	
    Double_t ptError[27];
    Int_t nlinesDjordjevic = 0;
    
    TString fileNamePi0Djordjevic_0010 = "ExternalInputPbPb/Theory/DjordjevicPi02011/PiRaa_0-10.txt";
    ifstream  filePi0Djordjevic_0010;
    filePi0Djordjevic_0010.open(fileNamePi0Djordjevic_0010,ios_base::in);
    cout << fileNamePi0Djordjevic_0010 << endl;
    
    while(!filePi0Djordjevic_0010.eof() && nlinesDjordjevic < 27){
        filePi0Djordjevic_0010 >> ptPi0Djordjevic[nlinesDjordjevic] >> pi0RaaLowEdge_0010[nlinesDjordjevic] >> pi0RaaHighEdge_0010[nlinesDjordjevic];
        cout << nlinesDjordjevic << "\t "  << ptPi0Djordjevic[nlinesDjordjevic] << "\t "  << pi0RaaLowEdge_0010[nlinesDjordjevic] << "\t "  << pi0RaaHighEdge_0010[nlinesDjordjevic] << endl;
        ptError[nlinesDjordjevic] = 0.;
        pi0RaaDjordjevic_0010[nlinesDjordjevic] = (pi0RaaHighEdge_0010[nlinesDjordjevic] + pi0RaaLowEdge_0010[nlinesDjordjevic])/2;
        pi0RaaLowEdgeError_0010[nlinesDjordjevic] = pi0RaaDjordjevic_0010[nlinesDjordjevic] - pi0RaaLowEdge_0010[nlinesDjordjevic];
        pi0RaaHighEdgeError_0010[nlinesDjordjevic] = pi0RaaHighEdge_0010[nlinesDjordjevic] - pi0RaaDjordjevic_0010[nlinesDjordjevic];
        nlinesDjordjevic++;

    }
    filePi0Djordjevic_0010.close();
    TGraphAsymmErrors *graphPi0Djordjevic_0010 = new TGraphAsymmErrors(nlinesDjordjevic-1, ptPi0Djordjevic, pi0RaaDjordjevic_0010, ptError, ptError, pi0RaaLowEdgeError_0010, pi0RaaHighEdgeError_0010);    
    
    TString fileNamePi0Djordjevic_2050 = "ExternalInputPbPb/Theory/DjordjevicPi02011/PiRaa_20-50.txt";
    ifstream  filePi0Djordjevic_2050;
    filePi0Djordjevic_2050.open(fileNamePi0Djordjevic_2050,ios_base::in);
    cout << fileNamePi0Djordjevic_2050 << endl;
    
    Double_t pi0RaaDjordjevic_2050[27];
    Double_t pi0RaaLowEdge_2050[27];
    Double_t pi0RaaHighEdge_2050[27];    
    Double_t pi0RaaLowEdgeError_2050[27];
    Double_t pi0RaaHighEdgeError_2050[27];

    nlinesDjordjevic = 0;
    while(!filePi0Djordjevic_2050.eof()  && nlinesDjordjevic < 27){
        filePi0Djordjevic_2050 >> ptPi0Djordjevic[nlinesDjordjevic] >> pi0RaaLowEdge_2050[nlinesDjordjevic] >> pi0RaaHighEdge_2050[nlinesDjordjevic];
        cout << nlinesDjordjevic << "\t "  << ptPi0Djordjevic[nlinesDjordjevic] << "\t "  << pi0RaaLowEdge_2050[nlinesDjordjevic] << "\t "  << pi0RaaHighEdge_2050[nlinesDjordjevic] << endl;
        ptError[nlinesDjordjevic] = 0.;
        pi0RaaDjordjevic_2050[nlinesDjordjevic] = (pi0RaaHighEdge_2050[nlinesDjordjevic] + pi0RaaLowEdge_2050[nlinesDjordjevic])/2;
        pi0RaaLowEdgeError_2050[nlinesDjordjevic] = pi0RaaDjordjevic_2050[nlinesDjordjevic] - pi0RaaLowEdge_2050[nlinesDjordjevic];
        pi0RaaHighEdgeError_2050[nlinesDjordjevic] = pi0RaaHighEdge_2050[nlinesDjordjevic] - pi0RaaDjordjevic_2050[nlinesDjordjevic];
        nlinesDjordjevic++;

    }
    filePi0Djordjevic_2050.close();
    TGraphAsymmErrors *graphPi0Djordjevic_2050 = new TGraphAsymmErrors(nlinesDjordjevic-1, ptPi0Djordjevic, pi0RaaDjordjevic_2050, ptError, ptError, pi0RaaLowEdgeError_2050, pi0RaaHighEdgeError_2050);    


	//*********************************************************************************************************************************
    //***************************************************** Jet quenching  ************************************************************
    // prediction for centrality 0-10% for LHC11h data from arXiv:1506.00838
	
	// Pi0 Raa 
    Double_t ptPi0JetQuenching18_0010[100];
    Double_t ptPi0JetQuenching22_0010[100];
    Double_t ptPi0JetQuenching26_0010[100];
    
    Double_t pi0Raa18_0010[100];
    Double_t pi0Raa22_0010[100];
    Double_t pi0Raa26_0010[100];
    
    Double_t pi0RaaLowErr_0010[100];
    Double_t pi0RaaHighErr_0010[100];
    
    Double_t ptErr[100];

    Int_t nlinesJetQuenching = 0;
    
    TString fileNamePi0JetQuenching18_0010 = "ExternalInputPbPb/Theory/JetQuenching/0.9piraa/0.9piraa1.8.dat";
    ifstream  filePi0JetQuenching18_0010;
    filePi0JetQuenching18_0010.open(fileNamePi0JetQuenching18_0010,ios_base::in);
    cout << fileNamePi0JetQuenching18_0010 << endl;
    
    while(!filePi0JetQuenching18_0010.eof() && nlinesJetQuenching < 100){
        filePi0JetQuenching18_0010 >> ptPi0JetQuenching18_0010[nlinesJetQuenching] >> pi0Raa18_0010[nlinesJetQuenching]; 
        cout << nlinesJetQuenching << "\t "  << ptPi0JetQuenching18_0010[nlinesJetQuenching] << "\t "  << pi0Raa18_0010[nlinesJetQuenching] << endl;
        nlinesJetQuenching++;
    }
    filePi0JetQuenching18_0010.close();
    TGraph* graphPi0JetQuenching18_0010 = new TGraph(nlinesJetQuenching,ptPi0JetQuenching18_0010,pi0Raa18_0010); 
    
	
    nlinesJetQuenching = 0;
    TString fileNamePi0JetQuenching22_0010 = "ExternalInputPbPb/Theory/JetQuenching/0.9piraa/0.9piraa2.2.dat";
    ifstream  filePi0JetQuenching22_0010;
    filePi0JetQuenching22_0010.open(fileNamePi0JetQuenching22_0010,ios_base::in);
    cout << fileNamePi0JetQuenching22_0010 << endl;
    
    while(!filePi0JetQuenching22_0010.eof() && nlinesJetQuenching < 100){
        filePi0JetQuenching22_0010 >> ptPi0JetQuenching22_0010[nlinesJetQuenching] >> pi0Raa22_0010[nlinesJetQuenching]; 
        cout << nlinesJetQuenching << "\t "  << ptPi0JetQuenching22_0010[nlinesJetQuenching] << "\t "  << pi0Raa22_0010[nlinesJetQuenching] << endl;    
        nlinesJetQuenching++;
    }
    filePi0JetQuenching22_0010.close();
    TGraph* graphPi0JetQuenching22_0010 = new TGraph(nlinesJetQuenching,ptPi0JetQuenching22_0010,pi0Raa22_0010); 
    
    
    nlinesJetQuenching = 0;
    TString fileNamePi0JetQuenching26_0010 = "ExternalInputPbPb/Theory/JetQuenching/0.9piraa/0.9piraa2.6.dat";
    ifstream  filePi0JetQuenching26_0010;
    filePi0JetQuenching26_0010.open(fileNamePi0JetQuenching26_0010,ios_base::in);
    cout << fileNamePi0JetQuenching26_0010 << endl;
    
    while(!filePi0JetQuenching26_0010.eof() && nlinesJetQuenching < 100){
        filePi0JetQuenching26_0010 >> ptPi0JetQuenching26_0010[nlinesJetQuenching] >> pi0Raa26_0010[nlinesJetQuenching]; 
        cout << nlinesJetQuenching << "\t "  << ptPi0JetQuenching26_0010[nlinesJetQuenching] << "\t "  << pi0Raa26_0010[nlinesJetQuenching] << endl;
        nlinesJetQuenching++;
    }
    filePi0JetQuenching26_0010.close();    
    TGraph* graphPi0JetQuenching26_0010 = new TGraph(nlinesJetQuenching,ptPi0JetQuenching26_0010,pi0Raa26_0010); 

	
    //Putting the three together:
    for(Int_t k=0; k<nlinesJetQuenching; k++){
        ptErr[k] = 0.;
        pi0RaaHighErr_0010[k] = pi0Raa18_0010[k] - pi0Raa22_0010[k];
        pi0RaaLowErr_0010[k] = pi0Raa22_0010[k] - pi0Raa26_0010[k];
    }
    TGraphAsymmErrors *graphPi0RAAJetQuenching_0010 = new TGraphAsymmErrors(nlinesJetQuenching, ptPi0JetQuenching22_0010, pi0Raa22_0010, ptErr, ptErr, pi0RaaLowErr_0010, pi0RaaHighErr_0010);

    
	// Eta Raa
    Double_t ptEtaJetQuenching18_0010[100];
    Double_t ptEtaJetQuenching22_0010[100];
    Double_t ptEtaJetQuenching26_0010[100];
    
    Double_t etaRaa18_0010[100];
    Double_t etaRaa22_0010[100];
    Double_t etaRaa26_0010[100];

    Double_t etaRaaLowErr_0010[100];
    Double_t etaRaaHighErr_0010[100];
    
    TString fileNameEtaJetQuenching18_0010 = "ExternalInputPbPb/Theory/JetQuenching/0.9etaraa/0.9etaraa1.8.dat";
    ifstream  fileEtaJetQuenching18_0010;
    fileEtaJetQuenching18_0010.open(fileNameEtaJetQuenching18_0010,ios_base::in);
    cout << fileNameEtaJetQuenching18_0010 << endl;
    
    nlinesJetQuenching = 0;
    while(!fileEtaJetQuenching18_0010.eof() && nlinesJetQuenching < 100){
        fileEtaJetQuenching18_0010 >> ptEtaJetQuenching18_0010[nlinesJetQuenching] >> etaRaa18_0010[nlinesJetQuenching]; 
        cout << nlinesJetQuenching << "\t " << ptEtaJetQuenching18_0010[nlinesJetQuenching] << "\t "  << etaRaa18_0010[nlinesJetQuenching] << endl;;
        nlinesJetQuenching++;
    }
    fileEtaJetQuenching18_0010.close();
    TGraph* graphEtaJetQuenching18_0010 = new TGraph(nlinesJetQuenching,ptEtaJetQuenching18_0010,etaRaa18_0010); 

    
    nlinesJetQuenching = 0;
    TString fileNameEtaJetQuenching22_0010 = "ExternalInputPbPb/Theory/JetQuenching/0.9etaraa/0.9etaraa2.2.dat";
    ifstream  fileEtaJetQuenching22_0010;
    fileEtaJetQuenching22_0010.open(fileNameEtaJetQuenching22_0010,ios_base::in);
    cout << fileNameEtaJetQuenching22_0010 << endl;
    
    while(!fileEtaJetQuenching22_0010.eof() && nlinesJetQuenching < 100){
        fileEtaJetQuenching22_0010 >> ptEtaJetQuenching22_0010[nlinesJetQuenching] >> etaRaa22_0010[nlinesJetQuenching]; 
        cout << nlinesJetQuenching << "\t "  << ptEtaJetQuenching22_0010[nlinesJetQuenching] << "\t "  << etaRaa22_0010[nlinesJetQuenching] << endl;;
        nlinesJetQuenching++;
    }
    fileEtaJetQuenching22_0010.close();
    TGraph* graphEtaJetQuenching22_0010 = new TGraph(nlinesJetQuenching,ptEtaJetQuenching22_0010,etaRaa22_0010); 
    
    
    nlinesJetQuenching = 0;
    TString fileNameEtaJetQuenching26_0010 = "ExternalInputPbPb/Theory/JetQuenching/0.9etaraa/0.9etaraa2.6.dat";
    ifstream  fileEtaJetQuenching26_0010;
    fileEtaJetQuenching26_0010.open(fileNameEtaJetQuenching26_0010,ios_base::in);
    cout << fileNameEtaJetQuenching26_0010 << endl;
    
    while(!fileEtaJetQuenching26_0010.eof() && nlinesJetQuenching < 100){
        fileEtaJetQuenching26_0010 >> ptEtaJetQuenching26_0010[nlinesJetQuenching] >> etaRaa26_0010[nlinesJetQuenching]; 
        cout << nlinesJetQuenching << "\t "  << ptEtaJetQuenching26_0010[nlinesJetQuenching] << "\t "  << etaRaa26_0010[nlinesJetQuenching] << endl;;
        nlinesJetQuenching++;
    }
    fileEtaJetQuenching26_0010.close();    
    TGraph* graphEtaJetQuenching26_0010 = new TGraph(nlinesJetQuenching,ptEtaJetQuenching26_0010,etaRaa26_0010); 
    
    
    //Putting the three together:
    for(Int_t k=0; k<nlinesJetQuenching; k++){
        etaRaaLowErr_0010[k] = etaRaa22_0010[k] - etaRaa26_0010[k];
        etaRaaHighErr_0010[k] = etaRaa18_0010[k] - etaRaa22_0010[k];
    }
    TGraphAsymmErrors *graphEtaRAAJetQuenching_0010 = new TGraphAsymmErrors(nlinesJetQuenching, ptEtaJetQuenching22_0010, etaRaa22_0010, ptErr, ptErr, etaRaaLowErr_0010, etaRaaHighErr_0010);

    
	// Eta to Pi0 ratio
    Double_t ptetaTopi0Ratio18_0010[100];
    Double_t ptetaTopi0Ratio22_0010[100];
    Double_t ptetaTopi0Ratio26_0010[100];

    Double_t etaTopi0Ratio18_0010[100];
    Double_t etaTopi0Ratio22_0010[100];
    Double_t etaTopi0Ratio26_0010[100];
    
    Double_t etaTopi0LowErr_0010[100];
    Double_t etaTopi0HighErr_0010[100];

    TString fileNameEtatoPi0RatioJetQuenching18_0010 = "ExternalInputPbPb/Theory/JetQuenching/0.9ratio/0.9ratio1.8.dat";
    ifstream  fileEtatoPi0RatioJetQuenching18_0010;
    fileEtatoPi0RatioJetQuenching18_0010.open(fileNameEtatoPi0RatioJetQuenching18_0010,ios_base::in);
    cout << fileNameEtatoPi0RatioJetQuenching18_0010 << endl;
    
    nlinesJetQuenching = 0;
    while(!fileEtatoPi0RatioJetQuenching18_0010.eof() && nlinesJetQuenching < 100){
        fileEtatoPi0RatioJetQuenching18_0010 >> ptetaTopi0Ratio18_0010[nlinesJetQuenching] >> etaTopi0Ratio18_0010[nlinesJetQuenching]; 
        cout << nlinesJetQuenching << "	 "  << ptetaTopi0Ratio18_0010[nlinesJetQuenching] << "	 "  << etaTopi0Ratio18_0010[nlinesJetQuenching] << endl;;
        nlinesJetQuenching++;
    }
    fileEtatoPi0RatioJetQuenching18_0010.close();
    TGraph* graphEtatoPi0RatioJetQuenching18_0010 = new TGraph(nlinesJetQuenching,ptetaTopi0Ratio18_0010,etaTopi0Ratio18_0010); 

    
    nlinesJetQuenching = 0;
    TString fileNameEtatoPi0RatioJetQuenching22_0010 = "ExternalInputPbPb/Theory/JetQuenching/0.9ratio/0.9ratio2.2.dat";
    ifstream  fileEtatoPi0RatioJetQuenching22_0010;
    fileEtatoPi0RatioJetQuenching22_0010.open(fileNameEtatoPi0RatioJetQuenching22_0010,ios_base::in);
    cout << fileNameEtatoPi0RatioJetQuenching22_0010 << endl;
    
    while(!fileEtatoPi0RatioJetQuenching22_0010.eof() && nlinesJetQuenching < 100){
        fileEtatoPi0RatioJetQuenching22_0010 >> ptetaTopi0Ratio22_0010[nlinesJetQuenching] >> etaTopi0Ratio22_0010[nlinesJetQuenching]; 
        cout << nlinesJetQuenching << "	 "  << ptetaTopi0Ratio22_0010[nlinesJetQuenching] << "	 "  << etaTopi0Ratio22_0010[nlinesJetQuenching] << endl;;
        nlinesJetQuenching++;
    }
    fileEtatoPi0RatioJetQuenching22_0010.close();
    TGraph* graphEtatoPi0RatioJetQuenching22_0010 = new TGraph(nlinesJetQuenching,ptetaTopi0Ratio22_0010,etaTopi0Ratio22_0010); 
    
    
    nlinesJetQuenching = 0;
    TString fileNameEtatoPi0RatioJetQuenching26_0010 = "ExternalInputPbPb/Theory/JetQuenching/0.9ratio/0.9ratio2.6.dat";
    ifstream  fileEtatoPi0RatioJetQuenching26_0010;
    fileEtatoPi0RatioJetQuenching26_0010.open(fileNameEtatoPi0RatioJetQuenching26_0010,ios_base::in);
    cout << fileNameEtatoPi0RatioJetQuenching26_0010 << endl;
    
    while(!fileEtatoPi0RatioJetQuenching26_0010.eof() && nlinesJetQuenching < 100){
        fileEtatoPi0RatioJetQuenching26_0010 >> ptetaTopi0Ratio22_0010[nlinesJetQuenching] >> etaTopi0Ratio26_0010[nlinesJetQuenching]; 
        cout << nlinesJetQuenching << "	 "  << ptetaTopi0Ratio22_0010[nlinesJetQuenching] << "	 "  << etaTopi0Ratio26_0010[nlinesJetQuenching] << endl;;
        nlinesJetQuenching++;
    }
    fileEtatoPi0RatioJetQuenching26_0010.close();    
    TGraph* graphEtatoPi0RatioJetQuenching26_0010 = new TGraph(nlinesJetQuenching,ptetaTopi0Ratio22_0010,etaTopi0Ratio26_0010); 
    
    
    //Putting the three together:
    for(Int_t k=0; k<nlinesJetQuenching; k++){
        etaTopi0LowErr_0010[k] = etaTopi0Ratio22_0010[k] - etaTopi0Ratio18_0010[k];
        etaTopi0HighErr_0010[k] = etaTopi0Ratio26_0010[k] - etaTopi0Ratio22_0010[k];
    }
    TGraphAsymmErrors *graphEtaToPi0JetQuenching_0010 = new TGraphAsymmErrors(nlinesJetQuenching, ptetaTopi0Ratio22_0010, etaTopi0Ratio22_0010, ptErr, ptErr, etaTopi0LowErr_0010, etaTopi0HighErr_0010);

	
    
    //********************************************************************************************************************************************    
    //*************************************** Raa theory Xiao-Fang *******************************************************************************
    Int_t index = 0;
    TGraphErrors* Xiao_Raa_0020 = new TGraphErrors(16);
    Double_t pT_0020[17];
    Double_t high_0020[17];
    Double_t low_0020[17];
    ifstream file_0020_low ("ExternalInputPbPb/Theory/Xiao/0-20-low.dat", ios::app);
    ifstream file_0020_high ("ExternalInputPbPb/Theory/Xiao/0-20-up.dat");
    TGraphErrors* Xiao_Raa_2040 = new TGraphErrors(16);
    Double_t pT_2040[17];
    Double_t high_2040[17];
    Double_t low_2040[17];
    ifstream file_2040_low ("ExternalInputPbPb/Theory/Xiao/20-40-low.dat");
    ifstream file_2040_high ("ExternalInputPbPb/Theory/Xiao/20-40-up.dat");
    TGraphErrors* Xiao_Raa_4060 = new TGraphErrors(16);
    Double_t pT_4060[17];
    Double_t high_4060[17];
    Double_t low_4060[17];
    ifstream file_4060_low ("ExternalInputPbPb/Theory/Xiao/40-60-low.dat");
    ifstream file_4060_high ("ExternalInputPbPb/Theory/Xiao/40-60-up.dat");
    TGraphErrors* Xiao_Raa_6080 = new TGraphErrors(16);
    Double_t pT_6080[17];
    Double_t high_6080[17];
    Double_t low_6080[17];
    ifstream file_6080_low ("ExternalInputPbPb/Theory/Xiao/60-80-low.dat");
    ifstream file_6080_high ("ExternalInputPbPb/Theory/Xiao/60-80-up.dat");

    if (file_0020_low.is_open()){
        while(!file_0020_low.eof()){
            file_0020_low >> pT_0020[index] >> low_0020[index];
            index++;
        }
        file_0020_low.close();
        index = 0;
    }
    if (file_0020_high.is_open()){
        while(!file_0020_high.eof()){
            file_0020_high >> pT_0020[index] >> high_0020[index];
            index++;
        }
        file_0020_high.close();
        index = 0;
    }
    
    if (file_2040_low.is_open()){
        while(!file_2040_low.eof()){
            file_2040_low >> pT_2040[index] >> low_2040[index];
            index++;
        }
        file_2040_low.close();
        index = 0;
    }
    if (file_2040_high.is_open()){
        while(!file_2040_high.eof()){
            file_2040_high >> pT_2040[index] >> high_2040[index];
            index++;
        }
        file_2040_high.close();
        index = 0;
    }
    
    if (file_4060_low.is_open()){
        while(!file_4060_low.eof()){
            file_4060_low >> pT_4060[index] >> low_4060[index];
            index++;
        }
        file_4060_low.close();
        index = 0;
    }
    if (file_4060_high.is_open()){
        while(!file_4060_high.eof()){
            file_4060_high >> pT_4060[index] >> high_4060[index];
            index++;
        }
        file_4060_high.close();
        index = 0;
    }
    
    if (file_6080_low.is_open()){
        while(!file_6080_low.eof()){
            file_6080_low >> pT_6080[index] >> low_6080[index];
            index++;
        }
        file_6080_low.close();
        index = 0;
    }
    if (file_6080_high.is_open()){
        while(!file_6080_high.eof()){
            file_6080_high >> pT_6080[index] >> high_6080[index];
            index++;
        }
        file_6080_high.close();
        index = 0;
    }
    
    for(Int_t i=0; i<16; i++){
        Xiao_Raa_0020->SetPoint(i, pT_0020[i], (high_0020[i] + low_0020[i])/2.);
        Xiao_Raa_0020->SetPointError(i, 0.0000000001, (high_0020[i] - low_0020[i])/2.);
        Xiao_Raa_2040->SetPoint(i, pT_2040[i], (high_2040[i] + low_2040[i])/2.);
        Xiao_Raa_2040->SetPointError(i, 0.0000000001, (high_2040[i] - low_2040[i])/2.);
        Xiao_Raa_4060->SetPoint(i, pT_4060[i], (high_4060[i] + low_4060[i])/2.);
        Xiao_Raa_4060->SetPointError(i, 0.0000000001, (high_4060[i] - low_4060[i])/2.);
        Xiao_Raa_6080->SetPoint(i, pT_6080[i], (high_6080[i] + low_6080[i])/2.);
        Xiao_Raa_6080->SetPointError(i, 0.0000000001, (high_6080[i] - low_6080[i])/2.);
    }

    // *************************************************************************************
    // **************************** Vitev **************************************************
    TGraphErrors* Vitev_Bas_Raa_0020 = new TGraphErrors(35);
    Double_t pTVitev_Bas_0020[40];
    Double_t highVitev_Bas_0020[40];
    Double_t lowVitev_Bas_0020[40];
    ifstream fileVitev_Bas_0020_low ("ExternalInputPbPb/Theory/Vitev/R-PbPb2760pi0.dnbasPI0");
    ifstream fileVitev_Bas_0020_high ("ExternalInputPbPb/Theory/Vitev/R-PbPb2760pi0.upbasPI0");
    TGraphErrors* Vitev_ShlSel_Raa_0020 = new TGraphErrors(35);
    Double_t pTVitev_ShlSel_0020[40];
    Double_t highVitev_ShlSel_0020[40];
    Double_t lowVitev_ShlSel_0020[40];
    ifstream fileVitev_ShlSel_0020_low ("ExternalInputPbPb/Theory/Vitev/R-PbPb2760pi0.dnShISelPI0");
    ifstream fileVitev_ShlSel_0020_high ("ExternalInputPbPb/Theory/Vitev/R-PbPb2760pi0.upShISelPI0");
    
   index = 0;
   if (fileVitev_Bas_0020_low.is_open()){
      while(!fileVitev_Bas_0020_low.eof()){
         fileVitev_Bas_0020_low >> pTVitev_Bas_0020[index] >> lowVitev_Bas_0020[index];
         index++;
      }
      fileVitev_Bas_0020_low.close();
      index = 0;
   }
   if (fileVitev_Bas_0020_high.is_open()){
      while(!fileVitev_Bas_0020_high.eof()){
         fileVitev_Bas_0020_high >> pTVitev_Bas_0020[index] >> highVitev_Bas_0020[index];
         index++;
      }
      fileVitev_Bas_0020_high.close();
      index = 0;
   }
   
   
   if (fileVitev_ShlSel_0020_low.is_open()){
      while(!fileVitev_ShlSel_0020_low.eof()){
         fileVitev_ShlSel_0020_low >> pTVitev_ShlSel_0020[index] >> lowVitev_ShlSel_0020[index];
         index++;
      }
      fileVitev_ShlSel_0020_low.close();
      index = 0;
   }
   if (fileVitev_ShlSel_0020_high.is_open()){
      while(!fileVitev_ShlSel_0020_high.eof()){
         fileVitev_ShlSel_0020_high >> pTVitev_ShlSel_0020[index] >> highVitev_ShlSel_0020[index];
         index++;
      }
      fileVitev_ShlSel_0020_high.close();
      index = 0;
   }

   for(Int_t i=0; i<35; i++){
      Vitev_Bas_Raa_0020->SetPoint(i, pTVitev_Bas_0020[i], (highVitev_Bas_0020[i] + lowVitev_Bas_0020[i])/2.);
      Vitev_Bas_Raa_0020->SetPointError(i, 0.0000000001, (highVitev_Bas_0020[i] - lowVitev_Bas_0020[i])/2.);
//       cout << "Vitev Bas \t " << pTVitev_Bas_0020[i] << "\t" << (highVitev_Bas_0020[i] + lowVitev_Bas_0020[i])/2. << "\t" << (highVitev_Bas_0020[i] - lowVitev_Bas_0020[i])/2. << endl;
   }
   for(Int_t i=0; i<35; i++){
      Vitev_ShlSel_Raa_0020->SetPoint(i, pTVitev_ShlSel_0020[i], (highVitev_ShlSel_0020[i] + lowVitev_ShlSel_0020[i])/2.);
      Vitev_ShlSel_Raa_0020->SetPointError(i, 0.0000000001, (highVitev_ShlSel_0020[i] - lowVitev_ShlSel_0020[i])/2.);
//       cout << "Vitev ShlSel \t " << pTVitev_ShlSel_0020[i] << "\t" << (highVitev_ShlSel_0020[i] + lowVitev_ShlSel_0020[i])/2. << "\t" << (highVitev_ShlSel_0020[i] - lowVitev_ShlSel_0020[i])/2. << endl;
      
   }

   TGraphErrors* Vitev_Bas_Raa_0005 = new TGraphErrors(275);
   Double_t pTVitev_Bas_0005[275];
   Double_t highVitev_Bas_0005[275];
   Double_t lowVitev_Bas_0005[275];
   ifstream fileVitev_Bas_0005_low ("ExternalInputPbPb/Theory/Vitev/R-PbPb2760pi0.05.dn");
   ifstream fileVitev_Bas_0005_high ("ExternalInputPbPb/Theory/Vitev/R-PbPb2760pi0.05.up");
   index = 0;
   if (fileVitev_Bas_0005_low.is_open()){
      while(!fileVitev_Bas_0005_low.eof()){
         fileVitev_Bas_0005_low >> pTVitev_Bas_0005[index] >> lowVitev_Bas_0005[index];
         index++;
      }
      fileVitev_Bas_0005_low.close();
      index = 0;
   }
   if (fileVitev_Bas_0005_high.is_open()){
      while(!fileVitev_Bas_0005_high.eof()){
         fileVitev_Bas_0005_high >> pTVitev_Bas_0005[index] >> highVitev_Bas_0005[index];
         index++;
      }
      fileVitev_Bas_0005_high.close();
      index = 0;
   }
   for(Int_t i=0; i<275; i++){
      Vitev_Bas_Raa_0005->SetPoint(i, pTVitev_Bas_0005[i], (highVitev_Bas_0005[i] + lowVitev_Bas_0005[i])/2.);
      Vitev_Bas_Raa_0005->SetPointError(i, 0.0000000001, (highVitev_Bas_0005[i] - lowVitev_Bas_0005[i])/2.);
//       cout << "Vitev Bas \t " << pTVitev_Bas_0005[i] << "\t" << (highVitev_Bas_0005[i] + lowVitev_Bas_0005[i])/2. << "\t" << (highVitev_Bas_0005[i] - lowVitev_Bas_0005[i])/2. << endl;
   }
   
   TGraphErrors* Vitev_Bas_Raa_0510 = new TGraphErrors(275);
   Double_t pTVitev_Bas_0510[275];
   Double_t highVitev_Bas_0510[275];
   Double_t lowVitev_Bas_0510[275];
   ifstream fileVitev_Bas_0510_low ("ExternalInputPbPb/Theory/Vitev/R-PbPb2760pi0.510.dn");
   ifstream fileVitev_Bas_0510_high ("ExternalInputPbPb/Theory/Vitev/R-PbPb2760pi0.510.up");
   index = 0;
   if (fileVitev_Bas_0510_low.is_open()){
      while(!fileVitev_Bas_0510_low.eof()){
         fileVitev_Bas_0510_low >> pTVitev_Bas_0510[index] >> lowVitev_Bas_0510[index];
         index++;
      }
      fileVitev_Bas_0510_low.close();
      index = 0;
   }
   if (fileVitev_Bas_0510_high.is_open()){
      while(!fileVitev_Bas_0510_high.eof()){
         fileVitev_Bas_0510_high >> pTVitev_Bas_0510[index] >> highVitev_Bas_0510[index];
         index++;
      }
      fileVitev_Bas_0510_high.close();
      index = 0;
   }
   for(Int_t i=0; i<275; i++){
      Vitev_Bas_Raa_0510->SetPoint(i, pTVitev_Bas_0510[i], (highVitev_Bas_0510[i] + lowVitev_Bas_0510[i])/2.);
      Vitev_Bas_Raa_0510->SetPointError(i, 0.0000000001, (highVitev_Bas_0510[i] - lowVitev_Bas_0510[i])/2.);
//       cout << "Vitev Bas \t " << pTVitev_Bas_0510[i] << "\t" << (highVitev_Bas_0510[i] + lowVitev_Bas_0510[i])/2. << "\t" << (highVitev_Bas_0510[i] - lowVitev_Bas_0510[i])/2. << endl;
   }
   
   TGraphErrors* Vitev_Bas_Raa_1020 = new TGraphErrors(275);
   Double_t pTVitev_Bas_1020[275];
   Double_t highVitev_Bas_1020[275];
   Double_t lowVitev_Bas_1020[275];
   ifstream fileVitev_Bas_1020_low ("ExternalInputPbPb/Theory/Vitev/R-PbPb2760pi0.1020.dn");
   ifstream fileVitev_Bas_1020_high ("ExternalInputPbPb/Theory/Vitev/R-PbPb2760pi0.1020.up");
   index = 0;
   if (fileVitev_Bas_1020_low.is_open()){
      while(!fileVitev_Bas_1020_low.eof()){
         fileVitev_Bas_1020_low >> pTVitev_Bas_1020[index] >> lowVitev_Bas_1020[index];
         index++;
      }
      fileVitev_Bas_1020_low.close();
      index = 0;
   }
   if (fileVitev_Bas_1020_high.is_open()){
      while(!fileVitev_Bas_1020_high.eof()){
         fileVitev_Bas_1020_high >> pTVitev_Bas_1020[index] >> highVitev_Bas_1020[index];
         index++;
      }
      fileVitev_Bas_1020_high.close();
      index = 0;
   }
   for(Int_t i=0; i<275; i++){
      Vitev_Bas_Raa_1020->SetPoint(i, pTVitev_Bas_1020[i], (highVitev_Bas_1020[i] + lowVitev_Bas_1020[i])/2.);
      Vitev_Bas_Raa_1020->SetPointError(i, 0.0000000001, (highVitev_Bas_1020[i] - lowVitev_Bas_1020[i])/2.);
//       cout << "Vitev Bas \t " << pTVitev_Bas_1020[i] << "\t" << (highVitev_Bas_1020[i] + lowVitev_Bas_1020[i])/2. << "\t" << (highVitev_Bas_1020[i] - lowVitev_Bas_1020[i])/2. << endl;
   }
   
   TGraphErrors* Vitev_Bas_Raa_2040 = new TGraphErrors(275);
   Double_t pTVitev_Bas_2040[275];
   Double_t highVitev_Bas_2040[275];
   Double_t lowVitev_Bas_2040[275];
   ifstream fileVitev_Bas_2040_low ("ExternalInputPbPb/Theory/Vitev/R-PbPb2760pi0.2040.dn");
   ifstream fileVitev_Bas_2040_high ("ExternalInputPbPb/Theory/Vitev/R-PbPb2760pi0.2040.up");
   index = 0;
   if (fileVitev_Bas_2040_low.is_open()){
      while(!fileVitev_Bas_2040_low.eof()){
         fileVitev_Bas_2040_low >> pTVitev_Bas_2040[index] >> lowVitev_Bas_2040[index];
         index++;
      }
      fileVitev_Bas_2040_low.close();
      index = 0;
   }
   if (fileVitev_Bas_2040_high.is_open()){
      while(!fileVitev_Bas_2040_high.eof()){
         fileVitev_Bas_2040_high >> pTVitev_Bas_2040[index] >> highVitev_Bas_2040[index];
         index++;
      }
      fileVitev_Bas_2040_high.close();
      index = 0;
   }
   for(Int_t i=0; i<275; i++){
      Vitev_Bas_Raa_2040->SetPoint(i, pTVitev_Bas_2040[i], (highVitev_Bas_2040[i] + lowVitev_Bas_2040[i])/2.);
      Vitev_Bas_Raa_2040->SetPointError(i, 0.0000000001, (highVitev_Bas_2040[i] - lowVitev_Bas_2040[i])/2.);
//       cout << "Vitev Bas \t " << pTVitev_Bas_2040[i] << "\t" << (highVitev_Bas_2040[i] + lowVitev_Bas_2040[i])/2. << "\t" << (highVitev_Bas_2040[i] - lowVitev_Bas_2040[i])/2. << endl;
   }
   
   TGraphErrors* Vitev_Bas_Raa_4060 = new TGraphErrors(275);
   Double_t pTVitev_Bas_4060[275];
   Double_t highVitev_Bas_4060[275];
   Double_t lowVitev_Bas_4060[275];
   ifstream fileVitev_Bas_4060_low ("ExternalInputPbPb/Theory/Vitev/R-PbPb2760pi0.4060.dn");
   ifstream fileVitev_Bas_4060_high ("ExternalInputPbPb/Theory/Vitev/R-PbPb2760pi0.4060.up");
   index = 0;
   if (fileVitev_Bas_4060_low.is_open()){
      while(!fileVitev_Bas_4060_low.eof()){
         fileVitev_Bas_4060_low >> pTVitev_Bas_4060[index] >> lowVitev_Bas_4060[index];
         index++;
      }
      fileVitev_Bas_4060_low.close();
      index = 0;
   }
   if (fileVitev_Bas_4060_high.is_open()){
      while(!fileVitev_Bas_4060_high.eof()){
         fileVitev_Bas_4060_high >> pTVitev_Bas_4060[index] >> highVitev_Bas_4060[index];
         index++;
      }
      fileVitev_Bas_4060_high.close();
      index = 0;
   }
   for(Int_t i=0; i<275; i++){
      Vitev_Bas_Raa_4060->SetPoint(i, pTVitev_Bas_4060[i], (highVitev_Bas_4060[i] + lowVitev_Bas_4060[i])/2.);
      Vitev_Bas_Raa_4060->SetPointError(i, 0.0000000001, (highVitev_Bas_4060[i] - lowVitev_Bas_4060[i])/2.);
//       cout << "Vitev Bas \t " << pTVitev_Bas_4060[i] << "\t" << (highVitev_Bas_4060[i] + lowVitev_Bas_4060[i])/2. << "\t" << (highVitev_Bas_4060[i] - lowVitev_Bas_4060[i])/2. << endl;
   }
   
   TGraphErrors* Vitev_Bas_Raa_6080 = new TGraphErrors(275);
   Double_t pTVitev_Bas_6080[275];
   Double_t highVitev_Bas_6080[275];
   Double_t lowVitev_Bas_6080[275];
   ifstream fileVitev_Bas_6080_low("ExternalInputPbPb/Theory/Vitev/R-PbPb2760pi0.6080.dn");
   ifstream fileVitev_Bas_6080_high("ExternalInputPbPb/Theory/Vitev/R-PbPb2760pi0.6080.up");

   index = 0;
   if (fileVitev_Bas_6080_low.is_open()){
      while(!fileVitev_Bas_6080_low.eof()){
         fileVitev_Bas_6080_low >> pTVitev_Bas_6080[index] >> lowVitev_Bas_6080[index];
         index++;
      }
      fileVitev_Bas_6080_low.close();
      index = 0;
   }
   if (fileVitev_Bas_6080_high.is_open()){
      while(!fileVitev_Bas_6080_high.eof()){
         fileVitev_Bas_6080_high >> pTVitev_Bas_6080[index] >> highVitev_Bas_6080[index];
         index++;
      }
      fileVitev_Bas_6080_high.close();
      index = 0;
   }
   for(Int_t i=0; i<275; i++){
      Vitev_Bas_Raa_6080->SetPoint(i, pTVitev_Bas_6080[i], (highVitev_Bas_6080[i] + lowVitev_Bas_6080[i])/2.);
      Vitev_Bas_Raa_6080->SetPointError(i, 0.0000000001, (highVitev_Bas_6080[i] - lowVitev_Bas_6080[i])/2.);
//       cout << "Vitev Bas \t " << pTVitev_Bas_6080[i] << "\t" << (highVitev_Bas_6080[i] + lowVitev_Bas_6080[i])/2. << "\t" << (highVitev_Bas_6080[i] - lowVitev_Bas_6080[i])/2. << endl;
   }
              
          
    //************************************************************************************************************************
    //*************************************** Raa theory WHDG ****************************************************************
    //from S. Wicks, W. Horowitz, M. Djordjevic and M. Gyulassy, ``Elastic, inelastic, and path length fluctuations in jet tomography,'' Nucl. Phys. A 784, 426 (2007) [nucl-th/0512076].
	//and ``The Surprising Transparency of the sQGP at LHC,'' Nucl. Phys. A 872, 265 (2011) [arXiv:1104.4958 [hep-ph]]

    //*****************************************************     Pi0     ************************************************************************
    // 0-5%, 5-10%, 0-10%, 10-20%, 0-20%, 20-40%, 20-50%, 40-60%, 60-80%
    
    TGraphAsymmErrors* gWHDG_Raa_0005 = new TGraphAsymmErrors(4);
    Double_t WHDG_pT_0005[39];
    Double_t WHDG_Raa_0005[39];
    Double_t WHDG_high_0005[39];
    Double_t WHDG_low_0005[39];
    ifstream WHDG_file_0005("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA0005b.dat");
    ifstream WHDG_file_0005_low("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA0005l.dat");
    ifstream WHDG_file_0005_high("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA0005u.dat");
    TGraphAsymmErrors* gWHDG_Raa_0510 = new TGraphAsymmErrors(4);
    Double_t WHDG_pT_0510[39];
    Double_t WHDG_Raa_0510[39];
    Double_t WHDG_high_0510[39];
    Double_t WHDG_low_0510[39];
    ifstream WHDG_file_0510("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA0510b.dat");
    ifstream WHDG_file_0510_low("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA0510l.dat");
    ifstream WHDG_file_0510_high("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA0510u.dat");
    TGraphAsymmErrors* gWHDG_Raa_0010 = new TGraphAsymmErrors(4);
    Double_t WHDG_pT_0010[39];
    Double_t WHDG_Raa_0010[39];
    Double_t WHDG_high_0010[39];
    Double_t WHDG_low_0010[39];
    ifstream WHDG_file_0010("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA0010b.dat");
    ifstream WHDG_file_0010_low("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA0010l.dat");
    ifstream WHDG_file_0010_high("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA0010u.dat");
    TGraphAsymmErrors* gWHDG_Raa_1020 = new TGraphAsymmErrors(4);
    Double_t WHDG_pT_1020[39];
    Double_t WHDG_Raa_1020[39];
    Double_t WHDG_high_1020[39];
    Double_t WHDG_low_1020[39];
    ifstream WHDG_file_1020("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA1020b.dat");
    ifstream WHDG_file_1020_low("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA1020l.dat");
    ifstream WHDG_file_1020_high("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA1020u.dat");
    
    TGraphAsymmErrors* gWHDG_Raa_0020 = new TGraphAsymmErrors(4);
    Double_t WHDG_pT_0020[39];
    Double_t WHDG_Raa_0020[39];
    Double_t WHDG_high_0020[39];
    Double_t WHDG_low_0020[39];
    ifstream WHDG_file_0020("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA0020b.dat");
    ifstream WHDG_file_0020_low("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA0020l.dat");
    ifstream WHDG_file_0020_high("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA0020u.dat");
    TGraphAsymmErrors* gWHDG_Raa_2040 = new TGraphAsymmErrors(4);
    Double_t WHDG_pT_2040[39];
    Double_t WHDG_Raa_2040[39];
    Double_t WHDG_high_2040[39];
    Double_t WHDG_low_2040[39];
    ifstream WHDG_file_2040("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA2040b.dat");
    ifstream WHDG_file_2040_low("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA2040l.dat");
    ifstream WHDG_file_2040_high("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA2040u.dat");
    TGraphAsymmErrors* gWHDG_Raa_2050 = new TGraphAsymmErrors(4);
    Double_t WHDG_pT_2050[6];
    Double_t WHDG_Raa_2050[6];
    Double_t WHDG_high_2050[6];
    Double_t WHDG_low_2050[6];
    ifstream WHDG_file_2050("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA2050b.dat");
    ifstream WHDG_file_2050_low("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA2050l.dat");
    ifstream WHDG_file_2050_high("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA2050u.dat");
    TGraphAsymmErrors* gWHDG_Raa_4060 = new TGraphAsymmErrors(4);
    Double_t WHDG_pT_4060[39];
    Double_t WHDG_Raa_4060[39];
    Double_t WHDG_high_4060[39];
    Double_t WHDG_low_4060[39];
    ifstream WHDG_file_4060("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA4060b.dat");
    ifstream WHDG_file_4060_low("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA4060l.dat");
    ifstream WHDG_file_4060_high("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA4060u.dat");
    TGraphAsymmErrors* gWHDG_Raa_6080 = new TGraphAsymmErrors(4);
    Double_t WHDG_pT_6080[39];
    Double_t WHDG_Raa_6080[39];
    Double_t WHDG_high_6080[39];
    Double_t WHDG_low_6080[39];
    ifstream WHDG_file_6080("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA6080b.dat");
    ifstream WHDG_file_6080_low("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA6080l.dat");
    ifstream WHDG_file_6080_high("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA6080u.dat");    
    
    if (WHDG_file_0005.is_open()){
        while(!WHDG_file_0005.eof()){
            WHDG_file_0005 >> WHDG_pT_0005[index] >> WHDG_Raa_0005[index];
         cout << WHDG_pT_0005[index] << "\t" << WHDG_Raa_0005[index] << endl; 
            index++;
        }
        WHDG_file_0005.close();
        index = 0;
    }
    if (WHDG_file_0005_low.is_open()){
        while(!WHDG_file_0005_low.eof()){
            WHDG_file_0005_low >> WHDG_pT_0005[index] >> WHDG_low_0005[index];
         cout << WHDG_pT_0005[index] << "\t" << WHDG_low_0005[index]<< endl; 
            index++;
        }
        WHDG_file_0005_low.close();
        index = 0;
    }
    if (WHDG_file_0005_high.is_open()){
        while(!WHDG_file_0005_high.eof()){
            WHDG_file_0005_high >> WHDG_pT_0005[index] >> WHDG_high_0005[index];
            index++;
        }
        WHDG_file_0005_high.close();
        index = 0;
    }
    
    if (WHDG_file_0510.is_open()){
        while(!WHDG_file_0510.eof()){
            WHDG_file_0510 >> WHDG_pT_0510[index] >> WHDG_Raa_0510[index];
            index++;
        }
        WHDG_file_0510.close();
        index = 0;
    }
    if (WHDG_file_0510_low.is_open()){
        while(!WHDG_file_0510_low.eof()){
            WHDG_file_0510_low >> WHDG_pT_0510[index] >> WHDG_low_0510[index];
            index++;
        }
        WHDG_file_0510_low.close();
        index = 0;
    }
    if (WHDG_file_0510_high.is_open()){
        while(!WHDG_file_0510_high.eof()){
            WHDG_file_0510_high >> WHDG_pT_0510[index] >> WHDG_high_0510[index];
            index++;
        }
        WHDG_file_0510_high.close();
        index = 0;
    }
    if (WHDG_file_0010.is_open()){
        while(!WHDG_file_0010.eof()){
            WHDG_file_0010 >> WHDG_pT_0010[index] >> WHDG_Raa_0010[index];
            index++;
        }
        WHDG_file_0010.close();
        index = 0;
    }
    if (WHDG_file_0010_low.is_open()){
        while(!WHDG_file_0010_low.eof()){
            WHDG_file_0010_low >> WHDG_pT_0010[index] >> WHDG_low_0010[index];
            index++;
        }
        WHDG_file_0010_low.close();
        index = 0;
    }
    if (WHDG_file_0010_high.is_open()){
        while(!WHDG_file_0010_high.eof()){
            WHDG_file_0010_high >> WHDG_pT_0010[index] >> WHDG_high_0010[index];
            index++;
        }
        WHDG_file_0010_high.close();
        index = 0;
    }
    if (WHDG_file_1020.is_open()){
        while(!WHDG_file_1020.eof()){
            WHDG_file_1020 >> WHDG_pT_1020[index] >> WHDG_Raa_1020[index];
            index++;
        }
        WHDG_file_1020.close();
        index = 0;
    }
    if (WHDG_file_1020_low.is_open()){
        while(!WHDG_file_1020_low.eof()){
            WHDG_file_1020_low >> WHDG_pT_1020[index] >> WHDG_low_1020[index];
            index++;
        }
        WHDG_file_1020_low.close();
        index = 0;
    }
    if (WHDG_file_1020_high.is_open()){
        while(!WHDG_file_1020_high.eof()){
            WHDG_file_1020_high >> WHDG_pT_1020[index] >> WHDG_high_1020[index];
            index++;
        }
        WHDG_file_1020_high.close();
        index = 0;
    }
    
    
    if (WHDG_file_0020.is_open()){
        while(!WHDG_file_0020.eof()){
            WHDG_file_0020 >> WHDG_pT_0020[index] >> WHDG_Raa_0020[index];
            index++;
        }
        WHDG_file_0020.close();
        index = 0;
    }
    if (WHDG_file_0020_low.is_open()){
        while(!WHDG_file_0020_low.eof()){
            WHDG_file_0020_low >> WHDG_pT_0020[index] >> WHDG_low_0020[index];
            index++;
        }
        WHDG_file_0020_low.close();
        index = 0;
    }
    if (WHDG_file_0020_high.is_open()){
        while(!WHDG_file_0020_high.eof()){
            WHDG_file_0020_high >> WHDG_pT_0020[index] >> WHDG_high_0020[index];
            index++;
        }
        WHDG_file_0020_high.close();
        index = 0;
    }
    
    if (WHDG_file_2040.is_open()){
        while(!WHDG_file_2040.eof()){
            WHDG_file_2040 >> WHDG_pT_2040[index] >> WHDG_Raa_2040[index];
            index++;
        }
        WHDG_file_2040.close();
        index = 0;
    }
    if (WHDG_file_2040_low.is_open()){
        while(!WHDG_file_2040_low.eof()){
            WHDG_file_2040_low >> WHDG_pT_2040[index] >> WHDG_low_2040[index];
            index++;
        }
        WHDG_file_2040_low.close();
        index = 0;
    }
    if (WHDG_file_2040_high.is_open()){
        while(!WHDG_file_2040_high.eof()){
            WHDG_file_2040_high >> WHDG_pT_2040[index] >> WHDG_high_2040[index];
            index++;
        }
        WHDG_file_2040_high.close();
        index = 0;
    }
    if (WHDG_file_2050.is_open()){
        while(!WHDG_file_2050.eof()){
            WHDG_file_2050 >> WHDG_pT_2050[index] >> WHDG_Raa_2050[index];
            index++;
        }
        WHDG_file_2050.close();
        index = 0;
    }

    if (WHDG_file_2050_low.is_open()){
        while(!WHDG_file_2050_low.eof()){
            WHDG_file_2050_low >> WHDG_pT_2050[index] >> WHDG_low_2050[index];
            index++;
        }
        WHDG_file_2050_low.close();
        index = 0;
    }
    if (WHDG_file_2050_high.is_open()){
        while(!WHDG_file_2050_high.eof()){
            WHDG_file_2050_high >> WHDG_pT_2050[index] >> WHDG_high_2050[index];
            index++;
        }
        WHDG_file_2050_high.close();
        index = 0;
    }
    
    if (WHDG_file_4060.is_open()){
        while(!WHDG_file_4060.eof()){
            WHDG_file_4060 >> WHDG_pT_4060[index] >> WHDG_Raa_4060[index];
            index++;
        }
        WHDG_file_4060.close();
        index = 0;
    }
    if (WHDG_file_4060_low.is_open()){
        while(!WHDG_file_4060_low.eof()){
            WHDG_file_4060_low >> WHDG_pT_4060[index] >> WHDG_low_4060[index];
            index++;
        }
        WHDG_file_4060_low.close();
        index = 0;
    }
    if (WHDG_file_4060_high.is_open()){
        while(!WHDG_file_4060_high.eof()){
            WHDG_file_4060_high >> WHDG_pT_4060[index] >> WHDG_high_4060[index];
            index++;
        }
        WHDG_file_4060_high.close();
        index = 0;
    }
    
    if (WHDG_file_6080.is_open()){
        while(!WHDG_file_6080.eof()){
            WHDG_file_6080 >> WHDG_pT_6080[index] >> WHDG_Raa_6080[index];
            index++;
        }
        WHDG_file_6080.close();
        index = 0;
    }
    if (WHDG_file_6080_low.is_open()){
        while(!WHDG_file_6080_low.eof()){
            WHDG_file_6080_low >> WHDG_pT_6080[index] >> WHDG_low_6080[index];
            index++;
        }
        WHDG_file_6080_low.close();
        index = 0;
    }
    if (WHDG_file_6080_high.is_open()){
        while(!WHDG_file_6080_high.eof()){
            WHDG_file_6080_high >> WHDG_pT_6080[index] >> WHDG_high_6080[index];
            index++;
        }
        WHDG_file_6080_high.close();
        index = 0;
    }    
    
    for(Int_t i=0; i<4; i++){
        gWHDG_Raa_0005->SetPoint(i, WHDG_pT_0005[i],  WHDG_Raa_0005[i]);
        gWHDG_Raa_0005->SetPointError(i, 0.000001, 0.000001, WHDG_Raa_0005[i] - WHDG_low_0005[i], WHDG_high_0005[i] - WHDG_Raa_0005[i]);
        gWHDG_Raa_0510->SetPoint(i, WHDG_pT_0510[i],  WHDG_Raa_0510[i]);
        gWHDG_Raa_0510->SetPointError(i, 0.000001, 0.000001, WHDG_Raa_0510[i] - WHDG_low_0510[i], WHDG_high_0510[i] - WHDG_Raa_0510[i]);
        gWHDG_Raa_0010->SetPoint(i, WHDG_pT_0010[i],  WHDG_Raa_0010[i]);
        gWHDG_Raa_0010->SetPointError(i, 0.000001, 0.000001, WHDG_Raa_0010[i] - WHDG_low_0010[i], WHDG_high_0010[i] - WHDG_Raa_0010[i]);
        gWHDG_Raa_1020->SetPoint(i, WHDG_pT_1020[i],  WHDG_Raa_1020[i]);
        gWHDG_Raa_1020->SetPointError(i, 0.000001, 0.000001, WHDG_Raa_1020[i] - WHDG_low_1020[i], WHDG_high_1020[i] - WHDG_Raa_1020[i]);
        gWHDG_Raa_0020->SetPoint(i, WHDG_pT_0020[i],  WHDG_Raa_0020[i]);
        gWHDG_Raa_0020->SetPointError(i, 0.000001, 0.000001, WHDG_Raa_0020[i] - WHDG_low_0020[i], WHDG_high_0020[i] - WHDG_Raa_0020[i]);
        gWHDG_Raa_2040->SetPoint(i, WHDG_pT_2040[i],  WHDG_Raa_2040[i]);
        gWHDG_Raa_2040->SetPointError(i, 0.000001, 0.000001, WHDG_Raa_2040[i] - WHDG_low_2040[i], WHDG_high_2040[i] - WHDG_Raa_2040[i]);
        gWHDG_Raa_2050->SetPoint(i, WHDG_pT_2050[i],  WHDG_Raa_2050[i]);
        gWHDG_Raa_2050->SetPointError(i, 0.000001, 0.000001, WHDG_Raa_2050[i] - WHDG_low_2050[i], WHDG_high_2050[i] - WHDG_Raa_2050[i]);
        gWHDG_Raa_4060->SetPoint(i, WHDG_pT_4060[i],  WHDG_Raa_4060[i]);
        gWHDG_Raa_4060->SetPointError(i, 0.000001, 0.000001, WHDG_Raa_4060[i] - WHDG_low_4060[i], WHDG_high_4060[i] - WHDG_Raa_4060[i]);
        gWHDG_Raa_6080->SetPoint(i, WHDG_pT_6080[i],  WHDG_Raa_6080[i]);
        gWHDG_Raa_6080->SetPointError(i, 0.000001, 0.000001, WHDG_Raa_6080[i] - WHDG_low_6080[i], WHDG_high_6080[i] - WHDG_Raa_6080[i]);
    }
    
    
    //*****************************************************     Eta     ************************************************************************
    // 0-10%, (available 10-20%), 20-50%, (available 50-80%)
    TGraphAsymmErrors* gWHDG_Eta_Raa_0010 = new TGraphAsymmErrors(4);
    Double_t WHDG_Eta_pT_0010[7];
    Double_t WHDG_Eta_Raa_0010[7];
    Double_t WHDG_Eta_high_0010[7];
    Double_t WHDG_Eta_low_0010[7];
    ifstream WHDG_Eta_file_0010("ExternalInputPbPb/Theory/WHDG2760etaRAA/WHDGetaRAA0010b.dat");
    ifstream WHDG_Eta_file_0010_low("ExternalInputPbPb/Theory/WHDG2760etaRAA/WHDGetaRAA0010l.dat");
    ifstream WHDG_Eta_file_0010_high("ExternalInputPbPb/Theory/WHDG2760etaRAA/WHDGetaRAA0010u.dat");

    TGraphAsymmErrors* gWHDG_Eta_Raa_2050 = new TGraphAsymmErrors(4);
    Double_t WHDG_Eta_pT_2050[7];
    Double_t WHDG_Eta_Raa_2050[7];
    Double_t WHDG_Eta_high_2050[7];
    Double_t WHDG_Eta_low_2050[7];
    ifstream WHDG_Eta_file_2050("ExternalInputPbPb/Theory/WHDG2760etaRAA/WHDGetaRAA2050b.dat");
    ifstream WHDG_Eta_file_2050_low("ExternalInputPbPb/Theory/WHDG2760etaRAA/WHDGetaRAA2050l.dat");
    ifstream WHDG_Eta_file_2050_high("ExternalInputPbPb/Theory/WHDG2760etaRAA/WHDGetaRAA2050u.dat");

    if (WHDG_Eta_file_0010.is_open()){
        while(!WHDG_Eta_file_0010.eof()){
            WHDG_Eta_file_0010 >> WHDG_Eta_pT_0010[index] >> WHDG_Eta_Raa_0010[index];
            index++;
        }
        WHDG_Eta_file_0010.close();
        index = 0;
    }
    if (WHDG_Eta_file_0010_low.is_open()){
        while(!WHDG_Eta_file_0010_low.eof()){
            WHDG_Eta_file_0010_low >> WHDG_Eta_pT_0010[index] >> WHDG_Eta_low_0010[index];
            index++;
        }
        WHDG_Eta_file_0010_low.close();
        index = 0;
    }
    if (WHDG_Eta_file_0010_high.is_open()){
        while(!WHDG_Eta_file_0010_high.eof()){
            WHDG_Eta_file_0010_high >> WHDG_Eta_pT_0010[index] >> WHDG_Eta_high_0010[index];
            index++;
        }
        WHDG_Eta_file_0010_high.close();
        index = 0;
    }

    if (WHDG_Eta_file_2050.is_open()){
        while(!WHDG_Eta_file_2050.eof()){
            WHDG_Eta_file_2050 >> WHDG_Eta_pT_2050[index] >> WHDG_Eta_Raa_2050[index];
            index++;
        }
        WHDG_Eta_file_2050.close();
        index = 0;
    }
    if (WHDG_Eta_file_2050_low.is_open()){
        while(!WHDG_Eta_file_2050_low.eof()){
            WHDG_Eta_file_2050_low >> WHDG_Eta_pT_2050[index] >> WHDG_Eta_low_2050[index];
            index++;
        }
        WHDG_Eta_file_2050_low.close();
        index = 0;
    }
    if (WHDG_Eta_file_2050_high.is_open()){
        while(!WHDG_Eta_file_2050_high.eof()){
            WHDG_Eta_file_2050_high >> WHDG_Eta_pT_2050[index] >> WHDG_Eta_high_2050[index];
            index++;
        }
        WHDG_Eta_file_2050_high.close();
        index = 0;
    }
    
    for(Int_t i=0; i<4; i++){
        gWHDG_Eta_Raa_0010->SetPoint(i, WHDG_Eta_pT_0010[i],  WHDG_Eta_Raa_0010[i]);
        gWHDG_Eta_Raa_0010->SetPointError(i, 0.000001, 0.000001, WHDG_Eta_Raa_0010[i] - WHDG_Eta_low_0010[i], WHDG_Eta_high_0010[i] - WHDG_Eta_Raa_0010[i]);

        gWHDG_Eta_Raa_2050->SetPoint(i, WHDG_Eta_pT_2050[i],  WHDG_Eta_Raa_2050[i]);
        gWHDG_Eta_Raa_2050->SetPointError(i, 0.000001, 0.000001, WHDG_Eta_Raa_2050[i] - WHDG_Eta_low_2050[i], WHDG_Eta_high_2050[i] - WHDG_Eta_Raa_2050[i]);
    }
    
    
    const Int_t nFilesEpos = 6;
    TString label[nFilesEpos] = {"Epos_pi0_PbPb_2760GeV_00_to_05", "Epos_pi0_PbPb_2760GeV_05_to_10", "Epos_pi0_PbPb_2760GeV_10_to_20", 
                    "Epos_pi0_PbPb_2760GeV_20_to_40", "Epos_pi0_PbPb_2760GeV_40_to_60", "Epos_pi0_PbPb_2760GeV_60_to_80"};
    TString labelOut[nFilesEpos] = {"graphEPOSSpec0005", "graphEPOSSpec0510", "graphEPOSSpec1020","graphEPOSSpec2040", "graphEPOSSpec4060", "graphEPOSSpec6080"};
    TString labelOutWOErr[nFilesEpos] = {"graphEPOSSpecWOErr0005", "graphEPOSSpecWOErr0510", "graphEPOSSpecWOErr1020","graphEPOSSpecWOErr2040", "graphEPOSSpecWOErr4060", "graphEPOSSpecWOErr6080"};
    
    TGraphErrors* gEpos[nFilesEpos];
    TGraph* gEposWOErr[nFilesEpos];
    
    for (int iFile=0; iFile<nFilesEpos; iFile++) {

        const Int_t nPoints = 70;
        Double_t pt[nPoints], dndptpy[nPoints], dndptpy_staterr[nPoints];

        TString filename = "ExternalInputPbPb/Theory/EPOS/" + label[iFile] + ".txt";
        cout << filename << endl;

        ifstream fEposTxt(filename);
        
        Int_t iPoint = 0;
        
        while(!fEposTxt.eof() && iPoint < nPoints){ 
        fEposTxt >> pt[iPoint] >> dndptpy[iPoint] >> dndptpy_staterr[iPoint];
        // cout << pt[iPoint] << "   " << dndptpy[iPoint] << "   " << dndptpy_staterr[iPoint] << endl;
        dndptpy[iPoint]= dndptpy[iPoint]/(pt[iPoint]*2*TMath::Pi());
        dndptpy_staterr[iPoint] = dndptpy_staterr[iPoint]/(pt[iPoint]*2*TMath::Pi());
        if (fEposTxt.fail()) {
            fEposTxt.clear();
            fEposTxt.ignore(0xFFFFFF,'\n');
            continue;
        }

        iPoint++;
        }

        fEposTxt.close();
    
        gEpos[iFile] = new TGraphErrors(nPoints, pt, dndptpy, 0, dndptpy_staterr);
        gEposWOErr[iFile] = new TGraph(nPoints, pt, dndptpy);
        Int_t i = nPoints-1 ;
        cout <<  pt[i] << endl;
        while (pt[i] > 12.){
        cout  << pt[i] << endl;
            gEpos[iFile]->RemovePoint(  gEpos[iFile]->GetN()-1);
        gEposWOErr[iFile]->RemovePoint(  gEposWOErr[iFile]->GetN()-1);
            i--;
        }
        gEpos[iFile]->Print();
    }    

	
	//*********************************************************************************************************************************
	//************************************************************** EPOS 3 ***********************************************************
    // Calculation for 0-10% only for LHC11h
	
    Double_t Pi0pt[70];
	Double_t Etapt[70];
	Double_t Pi0dndptpy[70];
	Double_t Pi0dndptpy_staterr[70];
	Double_t Etadndptpy[70];
	Double_t Etadndptpy_staterr[70];
	
    TString filenamePi0 = "ExternalInputPbPb/Theory/EPOS/pion0-10.txt";
    ifstream fEposPi0Txt;
	fEposPi0Txt.open(filenamePi0,ios_base::in);
    cout << filenamePi0 << endl;
	
    Int_t iPoint = 0;
    while(!fEposPi0Txt.eof() && iPoint < 70){ 
        fEposPi0Txt >> Pi0pt[iPoint] >> Pi0dndptpy[iPoint] >> Pi0dndptpy_staterr[iPoint];
        iPoint++;
    }
    fEposPi0Txt.close();
    TGraphErrors* graphEPOSPi0_0010 = new TGraphErrors(iPoint, Pi0pt, Pi0dndptpy, 0, Pi0dndptpy_staterr);

	
    TString filenameEta = "ExternalInputPbPb/Theory/EPOS/eta0-10.txt";
    ifstream fEposEtaTxt;
	fEposEtaTxt.open(filenameEta,ios_base::in);
    cout << filenameEta << endl;

    iPoint = 0;
    while(!fEposEtaTxt.eof() && iPoint < 70){ 
        fEposEtaTxt >> Etapt[iPoint] >> Etadndptpy[iPoint] >> Etadndptpy_staterr[iPoint];
        iPoint++;
    }
    fEposEtaTxt.close();
    TGraphErrors* graphEPOSEta_0010 = new TGraphErrors(iPoint, Etapt, Etadndptpy, 0, Etadndptpy_staterr);
    
    
    
	//*********************************************************************************************************************************
	//********************************************************** Nemchick *************************************************************
    
    Double_t pTAdditional[43] = {  6.00,  6.25,  6.50,  6.75,  7.00,
                                7.25,  7.50,  7.75,  8.00,  8.25,
                                8.50,  8.75,  9.00,  9.25,  9.50,
                                9.75, 10.00, 10.25, 10.50, 10.75,
                                11.00, 11.25, 11.50, 11.75, 12.00,
                                12.25, 12.50, 12.75, 13.00, 13.25,
                                13.50, 13.75, 14.00, 14.25, 14.50,
                                14.75, 15.00, 15.50, 16.00, 17.00,
                                18.00, 19.00, 20.00};
    
    Double_t ncoll0005 = 1684.4;
    Double_t ncoll0510 = 1316;
    Double_t ncoll0010 = 1500;
    Double_t ncoll1020 = 921.2;
    Double_t ncoll0020 = 1210.6;
    Double_t ncoll2040 = 438.4;
    Double_t ncoll4060 = 127.7;
    Double_t ncoll4050 = 171.25;
    Double_t ncoll5060 = 84.28;
    Double_t ncoll6080 = 26.71;
    Double_t ncoll6070 = 37.855;
    Double_t ncoll7080 = 15.575;
    Double_t ncoll8090 = 6.293;
    Double_t ncoll7590 = 8.219;

    TF1 *spectrum2760GeV = new TF1("Levy_Dummy","[0] / ( 2 * TMath::Pi())*([1]-1.)*([1]-2.) / ([1]*[2]*([1]*[2]+0.134977*([1]-2.)))  * pow(1.+(sqrt(x*x+0.134977*0.134977)-0.134977)/([1]*[2]), -[1])");
    spectrum2760GeV->SetParameter(0,1.6450861188);
    spectrum2760GeV->SetParameter(1,7.2981303064);
    spectrum2760GeV->SetParameter(2,0.1402413223);
    
    TGraph* Kopeliovich_YieldHydro_0005 = new TGraph(97);
    TGraph* Kopeliovich_YieldHydroSmoothed_0005 = new TGraph(140);
    TGraph* Kopeliovich_YieldELoss_0005 = new TGraph(140);
    TGraph* Kopeliovich_YieldTotal_0005 = new TGraph(140);
    Double_t pTKopeliovich_Yield_0005[140];
    Double_t yKopeliovich_YieldHydro_0005[97];
    Double_t yELossKopeliovich_Yield_0005[140];
    Double_t yTotalKopeliovich_Yield_0005[140];
    ifstream file_Kopeliovich_Yield_0005 ("ExternalInputPbPb/Theory/Nemchick/y_AA-pi0-hydro-0-5L.dat");
    index = 0;
    if (file_Kopeliovich_Yield_0005.is_open()){
        while(!file_Kopeliovich_Yield_0005.eof()){
            file_Kopeliovich_Yield_0005 >> pTKopeliovich_Yield_0005[index] >> yKopeliovich_YieldHydro_0005[index];
            index++;
        }
        file_Kopeliovich_Yield_0005.close();
        index = 0;
    }
    for(Int_t i=0; i<97; i++){
        Kopeliovich_YieldHydro_0005->SetPoint(i, pTKopeliovich_Yield_0005[i], yKopeliovich_YieldHydro_0005[i]);
    }
    for (Int_t i= 0; i < 43; i++){
        pTKopeliovich_Yield_0005[i+97] = pTAdditional[i];
    }   
    
    TGraph* Kopeliovich_RAA_0005 = new TGraph(24);
    Double_t pTKopeliovich_RAA_0005[24];
    Double_t yKopeliovich_RAA_0005[24];
    
    ifstream file_Kopeliovich_RAA_0005 ("ExternalInputPbPb/Theory/Nemchick/r_AA-pQCD-0-5L.dat");
    index = 0;
    if (file_Kopeliovich_RAA_0005.is_open()){
        while(!file_Kopeliovich_RAA_0005.eof()){
            file_Kopeliovich_RAA_0005 >> pTKopeliovich_RAA_0005[index] >> yKopeliovich_RAA_0005[index];
            index++;
        }
        file_Kopeliovich_RAA_0005.close();
        index = 0;
    }
    for(Int_t i=0; i<24; i++){
        cout << pTKopeliovich_RAA_0005[i] << endl;
        Kopeliovich_RAA_0005->SetPoint(i, pTKopeliovich_RAA_0005[i], yKopeliovich_RAA_0005[i]);
    }
    Kopeliovich_RAA_0005->Print();
    TGraph* Kopeliovich_RAA_0005_finerBinning = new TGraph(97);
    for(Int_t i=0; i<140; i++){
        Double_t evaluatedRAA = Kopeliovich_RAA_0005->Eval(pTKopeliovich_Yield_0005[i], 0, "S");
    //       Double_t evaluatedHydro = Kopeliovich_YieldHydro_0005->Eval(pTKopeliovich_Yield_0005[i], 0, "S");
        Double_t evaluatedPP = spectrum2760GeV->Eval(pTKopeliovich_Yield_0005[i]);
        Kopeliovich_RAA_0005_finerBinning->SetPoint(i, pTKopeliovich_Yield_0005[i], evaluatedRAA);

            yELossKopeliovich_Yield_0005[i] = ncoll0005*evaluatedPP*evaluatedRAA; 
            if (i < 97 ){
                yTotalKopeliovich_Yield_0005[i] = yKopeliovich_YieldHydro_0005[i] + ncoll0005*evaluatedPP*evaluatedRAA;    
            } else {
                yTotalKopeliovich_Yield_0005[i] = ncoll0005*evaluatedPP*evaluatedRAA; 
            }   
    //       Kopeliovich_YieldHydroSmoothed_0005->SetPoint(i, pTKopeliovich_Yield_0005[i], evaluatedHydro);
        Kopeliovich_YieldELoss_0005->SetPoint(i, pTKopeliovich_Yield_0005[i], yELossKopeliovich_Yield_0005[i]);
        Kopeliovich_YieldTotal_0005->SetPoint(i, pTKopeliovich_Yield_0005[i], yTotalKopeliovich_Yield_0005[i]);
    }                          

    Int_t j = 0;
    while (pTKopeliovich_Yield_0005[j] < 3){
        j++;
        Kopeliovich_YieldELoss_0005->RemovePoint(0);
        Kopeliovich_YieldTotal_0005->RemovePoint(0);
    }
    

    TGraph* Kopeliovich_YieldHydro_2040 = new TGraph(97);
    TGraph* Kopeliovich_YieldHydroSmoothed_2040 = new TGraph(140);
    TGraph* Kopeliovich_YieldELoss_2040 = new TGraph(140);
    TGraph* Kopeliovich_YieldTotal_2040 = new TGraph(140);
    Double_t pTKopeliovich_Yield_2040[140];
    Double_t yKopeliovich_YieldHydro_2040[97];
    Double_t yELossKopeliovich_Yield_2040[140];
    Double_t yTotalKopeliovich_Yield_2040[140];
    ifstream file_Kopeliovich_Yield_2040 ("ExternalInputPbPb/Theory/Nemchick/y_AA-pi0-hydro-20-40L.dat");
    index = 0;
    if (file_Kopeliovich_Yield_2040.is_open()){
        while(!file_Kopeliovich_Yield_2040.eof()){
            file_Kopeliovich_Yield_2040 >> pTKopeliovich_Yield_2040[index] >> yKopeliovich_YieldHydro_2040[index];
            index++;
        }
        file_Kopeliovich_Yield_2040.close();
        index = 0;
    }
    for(Int_t i=0; i<97; i++){
        Kopeliovich_YieldHydro_2040->SetPoint(i, pTKopeliovich_Yield_2040[i], yKopeliovich_YieldHydro_2040[i]);
    }
    for (Int_t i= 0; i < 43; i++){
        pTKopeliovich_Yield_2040[i+97] = pTAdditional[i];
    }   
    
    TGraph* Kopeliovich_RAA_2040 = new TGraph(24);
    Double_t pTKopeliovich_RAA_2040[24];
    Double_t yKopeliovich_RAA_2040[24];
    
    ifstream file_Kopeliovich_RAA_2040 ("ExternalInputPbPb/Theory/Nemchick/r_AA-pQCD-20-40L.dat");
    index = 0;
    if (file_Kopeliovich_RAA_2040.is_open()){
        while(!file_Kopeliovich_RAA_2040.eof()){
            file_Kopeliovich_RAA_2040 >> pTKopeliovich_RAA_2040[index] >> yKopeliovich_RAA_2040[index];
            index++;
        }
        file_Kopeliovich_RAA_2040.close();
        index = 0;
    }
    for(Int_t i=0; i<24; i++){
        cout << pTKopeliovich_RAA_2040[i] << endl;
        Kopeliovich_RAA_2040->SetPoint(i, pTKopeliovich_RAA_2040[i], yKopeliovich_RAA_2040[i]);
    }
    Kopeliovich_RAA_2040->Print();
    TGraph* Kopeliovich_RAA_2040_finerBinning = new TGraph(97);
    for(Int_t i=0; i<140; i++){
        Double_t evaluatedRAA = Kopeliovich_RAA_2040->Eval(pTKopeliovich_Yield_2040[i], 0, "S");
    //       Double_t evaluatedHydro = Kopeliovich_YieldHydro_2040->Eval(pTKopeliovich_Yield_2040[i], 0, "S");
        Double_t evaluatedPP = spectrum2760GeV->Eval(pTKopeliovich_Yield_2040[i]);
        Kopeliovich_RAA_2040_finerBinning->SetPoint(i, pTKopeliovich_Yield_2040[i], evaluatedRAA);

            yELossKopeliovich_Yield_2040[i] = ncoll2040*evaluatedPP*evaluatedRAA; 
            if (i < 97 ){
                yTotalKopeliovich_Yield_2040[i] = yKopeliovich_YieldHydro_2040[i] + ncoll2040*evaluatedPP*evaluatedRAA;    
            } else {
                yTotalKopeliovich_Yield_2040[i] = ncoll2040*evaluatedPP*evaluatedRAA; 
            }   
    //       Kopeliovich_YieldHydroSmoothed_2040->SetPoint(i, pTKopeliovich_Yield_2040[i], evaluatedHydro);
        Kopeliovich_YieldELoss_2040->SetPoint(i, pTKopeliovich_Yield_2040[i], yELossKopeliovich_Yield_2040[i]);
        Kopeliovich_YieldTotal_2040->SetPoint(i, pTKopeliovich_Yield_2040[i], yTotalKopeliovich_Yield_2040[i]);
    }                          

    j = 0;
    while (pTKopeliovich_Yield_2040[j] < 3){
        j++;
        Kopeliovich_YieldELoss_2040->RemovePoint(0);
        Kopeliovich_YieldTotal_2040->RemovePoint(0);
    }
    
    TGraph* Kopeliovich_YieldHydro_6080 = new TGraph(79);
    TGraph* Kopeliovich_YieldHydroSmoothed_6080 = new TGraph(122);
    TGraph* Kopeliovich_YieldELoss_6080 = new TGraph(122);
    TGraph* Kopeliovich_YieldTotal_6080 = new TGraph(122);
    Double_t pTKopeliovich_Yield_6080[122];
    Double_t yKopeliovich_YieldHydro_6080[97];
    Double_t yELossKopeliovich_Yield_6080[122];
    Double_t yTotalKopeliovich_Yield_6080[122];
    ifstream file_Kopeliovich_Yield_6080 ("ExternalInputPbPb/Theory/Nemchick/y_AA-pi0-hydro-60-80L.dat");
    index = 0;
    if (file_Kopeliovich_Yield_6080.is_open()){
        while(!file_Kopeliovich_Yield_6080.eof()){
            file_Kopeliovich_Yield_6080 >> pTKopeliovich_Yield_6080[index] >> yKopeliovich_YieldHydro_6080[index];
            index++;
        }
        file_Kopeliovich_Yield_6080.close();
        index = 0;
    }
    for(Int_t i=0; i<79; i++){
        Kopeliovich_YieldHydro_6080->SetPoint(i, pTKopeliovich_Yield_6080[i], yKopeliovich_YieldHydro_6080[i]);
    }
    for (Int_t i= 0; i < 43; i++){
        pTKopeliovich_Yield_6080[i+79] = pTAdditional[i];
    }   
    Kopeliovich_YieldHydro_6080->Print();
    
    TGraph* Kopeliovich_RAA_6080 = new TGraph(24);
    Double_t pTKopeliovich_RAA_6080[24];
    Double_t yKopeliovich_RAA_6080[24];
    
    ifstream file_Kopeliovich_RAA_6080 ("ExternalInputPbPb/Theory/Nemchick/r_AA-pQCD-60-80L.dat");
    index = 0;
    if (file_Kopeliovich_RAA_6080.is_open()){
        while(!file_Kopeliovich_RAA_6080.eof()){
            file_Kopeliovich_RAA_6080 >> pTKopeliovich_RAA_6080[index] >> yKopeliovich_RAA_6080[index];
            index++;
        }
        file_Kopeliovich_RAA_6080.close();
        index = 0;
    }
    for(Int_t i=0; i<24; i++){
        cout << pTKopeliovich_RAA_6080[i] << endl;
        Kopeliovich_RAA_6080->SetPoint(i, pTKopeliovich_RAA_6080[i], yKopeliovich_RAA_6080[i]);
    }
    Kopeliovich_RAA_6080->Print();
    TGraph* Kopeliovich_RAA_6080_finerBinning = new TGraph(79);
    for(Int_t i=0; i<79+43; i++){
        Double_t evaluatedRAA = Kopeliovich_RAA_6080->Eval(pTKopeliovich_Yield_6080[i], 0, "S");
    //       Double_t evaluatedHydro = Kopeliovich_YieldHydro_6080->Eval(pTKopeliovich_Yield_6080[i], 0, "S");
        Double_t evaluatedPP = spectrum2760GeV->Eval(pTKopeliovich_Yield_6080[i]);
        Kopeliovich_RAA_6080_finerBinning->SetPoint(i, pTKopeliovich_Yield_6080[i], evaluatedRAA);
        
            yELossKopeliovich_Yield_6080[i] = ncoll6080*evaluatedPP*evaluatedRAA; 
            if (i < 79 ){
                yTotalKopeliovich_Yield_6080[i] = yKopeliovich_YieldHydro_6080[i] + ncoll6080*evaluatedPP*evaluatedRAA;    
            } else {
                yTotalKopeliovich_Yield_6080[i] = ncoll6080*evaluatedPP*evaluatedRAA; 
            }   
    //       Kopeliovich_YieldHydroSmoothed_6080->SetPoint(i, pTKopeliovich_Yield_6080[i], evaluatedHydro);
        cout << pTKopeliovich_Yield_6080[i] << "\t" << evaluatedPP << "\t" << evaluatedRAA << "\t" << yELossKopeliovich_Yield_6080[i] << endl;
        Kopeliovich_YieldELoss_6080->SetPoint(i, pTKopeliovich_Yield_6080[i], yELossKopeliovich_Yield_6080[i]);
        Kopeliovich_YieldTotal_6080->SetPoint(i, pTKopeliovich_Yield_6080[i], yTotalKopeliovich_Yield_6080[i]);
    }                          
    Kopeliovich_RAA_6080_finerBinning->Print();
    j = 0;
    while (pTKopeliovich_Yield_6080[j] < 3){
        j++;
        Kopeliovich_YieldELoss_6080->RemovePoint(0);
        Kopeliovich_YieldTotal_6080->RemovePoint(0);
    }
    


    //**************************************************************************************************************
    //************************************* pQCD, pp 2.76TeV Pi0 DSS14 *********************************************
    
    Double_t xSection2760GeVppINEL = 62.8*1e9;
    Double_t       ptNLODSS14Pi02760GeV[39];
    Double_t       ptErrNLODSS14Pi02760GeV[39];
    Double_t       muHalfDSS14Pi02760GeV[39];
    Double_t       muHalfDSS14Pi02760GeVforGraph[39];
    Double_t       muHalfErrDSS14Pi02760GeV[39];
    Double_t       muOneDSS14Pi02760GeV[39];
    Double_t       muOneErrDSS14Pi02760GeV[39];
    Double_t       muTwoDSS14Pi02760GeV[39];
    Double_t       muTwoDSS14Pi02760GeVforGraph[39];
    Double_t       muTwoErrDSS14Pi02760GeV[39];
    Int_t          nlinesNLODSS14Pi02760GeV = 0;
    
    TString fileNameNLODSS14Pi02760GeV = "ExternalInput/Theory/pi0dss14-2760gev-dsigmadpt-eta060-pt40.dat";
    ifstream  fileNLODSS14Pi02760GeV;
    fileNLODSS14Pi02760GeV.open(fileNameNLODSS14Pi02760GeV,ios_base::in);
    cout << fileNameNLODSS14Pi02760GeV << endl;
    
    while(!fileNLODSS14Pi02760GeV.eof() && nlinesNLODSS14Pi02760GeV < 39){
        TString garbage;
        fileNLODSS14Pi02760GeV >> ptNLODSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] >> garbage >> muHalfDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV]  >> muTwoDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV]; 
        
		muOneDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] = (muHalfDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] + muTwoDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV])/(2*2*TMath::Pi()*ptNLODSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV]);
        muOneErrDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] = (muHalfDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV]-muTwoDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV])/(2*2*TMath::Pi()*ptNLODSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV]);
        ptErrNLODSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] = 0.;
        cout << nlinesNLODSS14Pi02760GeV << "\t "  << ptNLODSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] << "\t "  << muHalfDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] << "\t "  << muOneDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] << "         "  << muTwoDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] << endl;;
        
        muHalfDSS14Pi02760GeVforGraph[nlinesNLODSS14Pi02760GeV] = muHalfDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV]/(2*TMath::Pi()*ptNLODSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV]);
        muTwoDSS14Pi02760GeVforGraph[nlinesNLODSS14Pi02760GeV] = muTwoDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV]/(2*TMath::Pi()*ptNLODSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV]);
        
        muHalfErrDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] = 0.;
        muTwoErrDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] = 0.;
        nlinesNLODSS14Pi02760GeV++;
        
    }
    fileNLODSS14Pi02760GeV.close();
	
    TGraphAsymmErrors* graphNLOCalcmuHalfDSS14InvSecPi02760GeV = new TGraphAsymmErrors(nlinesNLODSS14Pi02760GeV, ptNLODSS14Pi02760GeV, muHalfDSS14Pi02760GeVforGraph, ptErrNLODSS14Pi02760GeV,
																					   ptErrNLODSS14Pi02760GeV, muHalfErrDSS14Pi02760GeV, muHalfErrDSS14Pi02760GeV);
    TGraphAsymmErrors* graphNLOCalcmuTwoDSS14InvSecPi02760GeV = new TGraphAsymmErrors(nlinesNLODSS14Pi02760GeV, ptNLODSS14Pi02760GeV, muTwoDSS14Pi02760GeVforGraph, ptErrNLODSS14Pi02760GeV,
																					  ptErrNLODSS14Pi02760GeV,muTwoErrDSS14Pi02760GeV, muTwoErrDSS14Pi02760GeV);
    TGraphAsymmErrors* graphNLOCalcDSS14InvSecPi02760GeV = new TGraphAsymmErrors(nlinesNLODSS14Pi02760GeV, ptNLODSS14Pi02760GeV, muOneDSS14Pi02760GeV, ptErrNLODSS14Pi02760GeV, ptErrNLODSS14Pi02760GeV,
                                                                                muOneErrDSS14Pi02760GeV, muOneErrDSS14Pi02760GeV);
    TGraphAsymmErrors* graphNLOCalcDSS14InvYieldPi02760GeV = ScaleGraphAsym(graphNLOCalcDSS14InvSecPi02760GeV,  1./xSection2760GeVppINEL);
	

    // **************************************************************************************************************
    // ************************************* pQCD, pp 2.76TeV Eta DSS07 *********************************************
    // **************************************************************************************************************
    Double_t       ptNLODSS07Eta2760GeV[24];
    Double_t       ptErrNLODSS07Eta2760GeV[24];
    Double_t       muHalfDSS07Eta2760GeV[24];
    Double_t       muHalfDSS07Eta2760GeVforGraph[24];
    Double_t       muHalfErrDSS07Eta2760GeV[24];
    Double_t       muOneDSS07Eta2760GeV[24];
    Double_t       muOneErrDSS07Eta2760GeV[24];
    Double_t       muTwoDSS07Eta2760GeV[24];
    Double_t       muTwoDSS07Eta2760GeVforGraph[24];
    Double_t       muTwoErrDSS07Eta2760GeV[24];
    Int_t       nlinesNLODSS07Eta2760GeV = 0;

    TString fileNameNLODSS07Eta2760GeV = "ExternalInput/Theory/etadss07-2760gev-dsigmadpt-eta060-pt40.dat";
    ifstream  fileNLODSS07Eta2760GeV;
    fileNLODSS07Eta2760GeV.open(fileNameNLODSS07Eta2760GeV,ios_base::in);
    cout << fileNameNLODSS07Eta2760GeV << endl;
    
    while(!fileNLODSS07Eta2760GeV.eof() && nlinesNLODSS07Eta2760GeV < 24){
        TString garbage;
        fileNLODSS07Eta2760GeV >> ptNLODSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV] >> garbage >> muHalfDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV]  >> muTwoDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV];
		
        muOneDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV] = (muHalfDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV]+muTwoDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV])/ (2*2*TMath::Pi()*ptNLODSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV]);
        muOneErrDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV] = (muHalfDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV]-muTwoDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV])/(2*2*TMath::Pi()*ptNLODSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV]);
        ptErrNLODSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV] = 0.;
        cout << nlinesNLODSS07Eta2760GeV << "\t "  << ptNLODSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV] << "\t "  << muHalfDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV] << "\t "  << muOneDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV] << "\t "  << muTwoDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV] << endl;;
        
        muHalfDSS07Eta2760GeVforGraph[nlinesNLODSS07Eta2760GeV] = muHalfDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV]/(2*TMath::Pi()*ptNLODSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV]);
        muTwoDSS07Eta2760GeVforGraph[nlinesNLODSS07Eta2760GeV] = muTwoDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV]/(2*TMath::Pi()*ptNLODSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV]);
        muHalfErrDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV] = 0.;
        muTwoErrDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV] = 0.;
        nlinesNLODSS07Eta2760GeV++;
        
    }
    fileNLODSS07Eta2760GeV.close();
	
    TGraphAsymmErrors* graphNLOCalcmuHalfDSS07InvSecEta2760GeV = new TGraphAsymmErrors(nlinesNLODSS07Eta2760GeV, ptNLODSS07Eta2760GeV, muHalfDSS07Eta2760GeVforGraph, ptErrNLODSS07Eta2760GeV,
																					   ptErrNLODSS07Eta2760GeV, muHalfErrDSS07Eta2760GeV, muHalfErrDSS07Eta2760GeV);
    TGraphAsymmErrors* graphNLOCalcmuTwoDSS07InvSecEta2760GeV = new TGraphAsymmErrors(nlinesNLODSS07Eta2760GeV, ptNLODSS07Eta2760GeV, muTwoDSS07Eta2760GeVforGraph, ptErrNLODSS07Eta2760GeV,
																					  ptErrNLODSS07Eta2760GeV, muTwoErrDSS07Eta2760GeV, muTwoErrDSS07Eta2760GeV);
    TGraphAsymmErrors* graphNLOCalcDSS07InvSecEta2760GeV = new TGraphAsymmErrors(nlinesNLODSS07Eta2760GeV, ptNLODSS07Eta2760GeV, muOneDSS07Eta2760GeV, ptErrNLODSS07Eta2760GeV,
																				 ptErrNLODSS07Eta2760GeV, muOneErrDSS07Eta2760GeV, muOneErrDSS07Eta2760GeV);
    TGraphAsymmErrors* graphNLOCalcDSS07InvYieldEta2760GeV = ScaleGraphAsym(graphNLOCalcDSS07InvSecEta2760GeV, 1./xSection2760GeVppINEL);


	
	TFile *fileTheoryGraphsPbPb = new TFile(Form("ExternalInputPbPb/Theory/TheoryCompilationPbPb%s.root",specifier.Data()),"RECREATE");
        
        Vitev_Bas_Raa_0020->Write("graphVitevBasRAA0020");
        Vitev_Bas_Raa_0005->Write("graphVitevBasRAA0005");
        Vitev_Bas_Raa_0510->Write("graphVitevBasRAA0510");
        Vitev_Bas_Raa_1020->Write("graphVitevBasRAA1020");
        Vitev_Bas_Raa_2040->Write("graphVitevBasRAA2040");
        Vitev_Bas_Raa_4060->Write("graphVitevBasRAA4060");
        Vitev_Bas_Raa_6080->Write("graphVitevBasRAA6080");
        Vitev_ShlSel_Raa_0020->Write("graphVitevShlSelRAA0020");
        Xiao_Raa_0020->Write("graphXiaoRAA0020");
        Xiao_Raa_2040->Write("graphXiaoRAA2040");
        Xiao_Raa_4060->Write("graphXiaoRAA4060");
        Xiao_Raa_6080->Write("graphXiaoRAA6080");
        gWHDG_Raa_0005->Write("graphWHDGRAA0005");
        gWHDG_Raa_0510->Write("graphWHDGRAA0510");
        gWHDG_Raa_1020->Write("graphWHDGRAA1020");
        gWHDG_Raa_0010->Write("graphWHDGRAA0010");
        gWHDG_Raa_0020->Write("graphWHDGRAA0020");
        gWHDG_Raa_2040->Write("graphWHDGRAA2040");
        gWHDG_Raa_2050->Write("graphWHDGRAA2050");
        gWHDG_Raa_4060->Write("graphWHDGRAA4060");
        gWHDG_Raa_6080->Write("graphWHDGRAA6080");
        gWHDG_Eta_Raa_0010->Write("graphWHDGetaRAA0010");
        gWHDG_Eta_Raa_2050->Write("graphWHDGetaRAA2050");

        for (int iFile=0; iFile<nFilesEpos; iFile++) {
            gEpos[iFile]->Write(labelOut[iFile]);
			gEposWOErr[iFile]->Write(labelOutWOErr[iFile]);
        }
        Kopeliovich_YieldHydro_0005->Write("graphKopeliovichHydroYield0005");
        Kopeliovich_YieldELoss_0005->Write("graphKopeliovichELossYield0005");
        Kopeliovich_YieldTotal_0005->Write("graphKopeliovichTotalYield0005");
        Kopeliovich_RAA_0005->Write("graphKopeliovichRAA0005");
        Kopeliovich_RAA_0005_finerBinning->Write("graphKopeliovichRAA0005_finerBinning");
        Kopeliovich_YieldHydro_2040->Write("graphKopeliovichHydroYield2040");
        Kopeliovich_YieldELoss_2040->Write("graphKopeliovichELossYield2040");
        Kopeliovich_YieldTotal_2040->Write("graphKopeliovichTotalYield2040");
        Kopeliovich_RAA_2040->Write("graphKopeliovichRAA2040");
        Kopeliovich_RAA_2040_finerBinning->Write("graphKopeliovichRAA2040_finerBinning");
        Kopeliovich_YieldHydro_6080->Write("graphKopeliovichHydroYield6080");
        Kopeliovich_YieldELoss_6080->Write("graphKopeliovichELossYield6080");
        Kopeliovich_YieldTotal_6080->Write("graphKopeliovichTotalYield6080");
        Kopeliovich_RAA_6080->Write("graphKopeliovichRAA6080");
        Kopeliovich_RAA_6080_finerBinning->Write("graphKopeliovichRAA6080_finerBinning");

        TheoryCracowPi0LowPt_0010->Write("TheoryCracowPi0LowPt_0010");
        TheoryCracowEtaLowPt_0010->Write("TheoryCracowEtaLowPt_0010");
        TheoryCracowEtaToPi0LowPt_0010->Write("TheoryCracowEtaToPi0LowPt_0010");
        TheoryCracowPi0LowPt_2050->Write("TheoryCracowPi0LowPt_2050");
        TheoryCracowEtaLowPt_2050->Write("TheoryCracowEtaLowPt_2050");
        TheoryCracowEtaToPi0LowPt_2050->Write("TheoryCracowEtaToPi0LowPt_2050");
        TheoryCracowChargedPionLowPt_0010->Write("TheoryCracowChargedPionLowPt_0010");
        TheoryCracowChargedKaonLowPt_0010->Write("TheoryCracowChargedKaonLowPt_0010");
        TheoryCracowKaonsToPionsLowPt_0010->Write("TheoryCracowKaonsToPionsLowPt_0010");
        TheoryCracowNeutralToChargedPionLowPt_0010->Write("TheoryCracowNeutralToChargedPionLowPt_0010");
		
		graphPi0Djordjevic_0010->Write("graphPi0Djordjevic_0010");
		graphPi0Djordjevic_2050->Write("graphPi0Djordjevic_2050");
		
        graphPi0JetQuenching18_0010->Write("graphPi0JetQuenching18_0010");
        graphPi0JetQuenching22_0010->Write("graphPi0JetQuenching22_0010");
        graphPi0JetQuenching26_0010->Write("graphPi0JetQuenching26_0010");
        graphEtaJetQuenching18_0010->Write("graphEtaJetQuenching18_0010");
        graphEtaJetQuenching22_0010->Write("graphEtaJetQuenching22_0010");
        graphEtaJetQuenching26_0010->Write("graphEtaJetQuenching26_0010");
        
        graphEtatoPi0RatioJetQuenching18_0010->Write("graphEtatoPi0RatioJetQuenching18_0010");
        graphEtatoPi0RatioJetQuenching22_0010->Write("graphEtatoPi0RatioJetQuenching22_0010");
        graphEtatoPi0RatioJetQuenching26_0010->Write("graphEtatoPi0RatioJetQuenching26_0010");
        
        graphPi0RAAJetQuenching_0010->Write("graphPi0RAAJetQuenching_0010");
        graphEtaRAAJetQuenching_0010->Write("graphEtaRAAJetQuenching_0010");
        graphEtaToPi0JetQuenching_0010->Write("graphEtaToPi0JetQuenching_0010");

        graphEPOSPi0_0010->Write("epos_pi0_pt_cent0_10");
        graphEPOSEta_0010->Write("epos_eta_pt_cent0_10");
		
		graphNLOCalcmuHalfDSS14InvSecPi02760GeV->Write("graphNLOCalcmuHalfDSS14InvSecPi02760GeV");
		graphNLOCalcmuTwoDSS14InvSecPi02760GeV->Write("graphNLOCalcmuTwoDSS14InvSecPi02760GeV");
		graphNLOCalcDSS14InvSecPi02760GeV->Write("graphNLOCalcDSS14InvCrossSecPi02760GeV");
		graphNLOCalcDSS14InvYieldPi02760GeV->Write("graphNLOCalcDSS14InvYieldPi02760GeV");

		graphNLOCalcmuHalfDSS07InvSecEta2760GeV->Write("graphNLOCalcmuHalfDSS07InvSecEta2760GeV");
		graphNLOCalcmuTwoDSS07InvSecEta2760GeV->Write("graphNLOCalcmuTwoDSS07InvSecEta2760GeV");
		graphNLOCalcDSS07InvSecEta2760GeV->Write("graphNLOCalcDSS07InvCrossSecEta2760GeV");
		graphNLOCalcDSS07InvYieldEta2760GeV->Write("graphNLOCalcDSS07InvYieldEta2760GeV");
        
    fileTheoryGraphsPbPb->Close();


}