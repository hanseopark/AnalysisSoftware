/*****************************************************************************
******         provided by Gamma Conversion Group, PWGGA,               ******
******        Friederike Bock, friederike.bock@cern.ch                  ******
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

Double_t xSection900GeV         = 47.78*1e-3;
Double_t xSection2760GeV        = 55.416*1e-3;
Double_t xSection7000GeV        = 62.22*1e-3;
Double_t xSection8000GeV        = 55.74*1e-3;
Double_t recalcBarn             = 1e12; //NLO in pbarn!!!!
Double_t xSection2760GeVINEL    = 62.8*1e-3;
Double_t xSection7000GeVINEL    = 73.2*1e-3;


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
    
    //**********************************************************************************************************************
    //******************************** 900 GeV Pi0 and Eta calc*************************************************************
    //**********************************************************************************************************************
    
    // Vogelsang
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
    
    // DSS
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

    // Eta/Pi0
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
    

    
    //**********************************************************************************************************************
    //******************************** 2.76TeV GeV Pi0 and Eta calc*********************************************************
    //**********************************************************************************************************************    

    // Vogelsang
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
    
    // DSS
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

    // Eta/pi0
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
            if (valueNLOMuHalfPi02760GeV[n] != 0){
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
    
    //*---CGC- T. Lappi and H.M"antysaari  arxiv1309.6963 (kt-factorization)---*//
    Int_t index = 0;
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
    


    // DSS14
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

    //**********************************************************************************************************************
    //******************************** 7 TeV Pi0 and Eta calc***************************************************************
    //**********************************************************************************************************************

    // Vogelsang
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

    
    // BKK
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

    // DSS
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

    // Eta/Pi0
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

    // DSS14
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
    
    
    //**********************************************************************************************************************
    //******************************** 8 TeV Pi0 and Eta calc***************************************************************
    //**********************************************************************************************************************
    
    // Vogelsang
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

    Double_t* valueNLOMuHalfEta8000GeV  = graphNLOCalcInvSecEtaMuHalf8000GeV->GetY();
    Double_t* valueNLOMuOneEta8000GeV   = graphNLOCalcInvSecEtaMuOne8000GeV->GetY();
    Double_t* valueNLOMuTwoEta8000GeV   = graphNLOCalcInvSecEtaMuTwo8000GeV->GetY();
    Double_t* valueNLOMuHalfPi08000GeV  = graphNLOCalcInvSecPi0MuHalf8000GeV->GetY();
    Double_t* valueNLOMuOnePi08000GeV   = graphNLOCalcInvSecPi0MuOne8000GeV->GetY();
    Double_t* valueNLOMuTwoPi08000GeV   = graphNLOCalcInvSecPi0MuTwo8000GeV->GetY();
    Double_t* xValueNLO8000GeV          = graphNLOCalcInvSecPi0MuOne8000GeV->GetX();
    Int_t xNBins8000GeV                 = graphNLOCalcInvSecPi0MuOne8000GeV->GetN();
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
    TGraph* graphEtaToPi0NLOMuHalf8TeV      = new TGraph(xNBins8000GeV,xValueNLO8000GeV,valueNLOEtaToPi0NLOMuHalf8000GeV);  
    graphEtaToPi0NLOMuHalf8TeV->RemovePoint(0);
    graphEtaToPi0NLOMuHalf8TeV->Print();
    TGraph* graphEtaToPi0NLOMuOne8TeV       = new TGraph(xNBins8000GeV,xValueNLO8000GeV,valueNLOEtaToPi0NLOMuOne8000GeV);   
    graphEtaToPi0NLOMuOne8TeV->Print();
    TGraph* graphEtaToPi0NLOMuTwo8TeV       = new TGraph(xNBins8000GeV,xValueNLO8000GeV,valueNLOEtaToPi0NLOMuTwo8000GeV);   
    graphEtaToPi0NLOMuTwo8TeV->Print();
    

    //**********************************************************************************************************************
    //***************************** Pythia calculations 2.76TeV ************************************************************
    //**********************************************************************************************************************
    
    // Klaus: 8.1 Tune 4C
    TFile* file2760GeV_Pythia8_fixedBinning     = TFile::Open("ExternalInput/Theory/Pythia/pp_pi0_2760GeV_pythia8_tune4C_comp_fixed_binning_2014-02-04-1330.root");
    TH1F* histoPythia8Spec2760GeV               = (TH1F*)file2760GeV_Pythia8_fixedBinning->Get("hPythiaComb");
    TH1F* histoPythia8InvYield2760GeV           = (TH1F*)histoPythia8Spec2760GeV->Clone("histoPythia8InvYield2760GeV");
    histoPythia8InvYield2760GeV->Scale(1./(xSection2760GeV*recalcBarn));
    
    TFile* filePythia8_2760GeV                  = TFile::Open("ExternalInput/Theory/Pythia/pp_pi0_2760GeV_pythia8_tune4C_comp_variable_binning.root");
    TH1F* histoPythia8Spec2760GeVVarBinning     = (TH1F*)filePythia8_2760GeV->Get("hPythiaComb");
    TH1F* histoPythia8InvYield2760GeVVarBinning = (TH1F*)histoPythia8Spec2760GeVVarBinning->Clone("histoPythia8InvYield2760GeV");
    histoPythia8InvYield2760GeVVarBinning->Scale(1./(xSection2760GeV*recalcBarn));

    
    TFile* filePythia8Monash2013_2760GeV        = TFile::Open("ExternalInput/Theory/Pythia/PYTHIA8_Monash2013Tune_2760GeV.root");
    TH1F* histoPi0Pythia8MonashInvSec           = (TH1F*)filePythia8Monash2013_2760GeV->Get("fHistInvXsec_Pi0");
    TH1F* histoEtaPythia8MonashInvSec           = (TH1F*)filePythia8Monash2013_2760GeV->Get("fHistInvXsec_Eta");
    
    //**********************************************************************************************************************
    //********************************* Write graphs and histos to compilation file for pp *********************************
    //**********************************************************************************************************************
    TFile fileTheoryGraphsPP("ExternalInput/Theory/TheoryCompilationPP.root","UPDATE");

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
        histoPi0Pythia8MonashInvSec->Write("histoInvSecPythia8Monash2013Pi02760GeV");
        histoEtaPythia8MonashInvSec->Write("histoInvSecPythia8Monash2013Eta2760GeV");

        
        graphNLOCalcDSS14InvSecPi02760GeV->Write("graphNLOCalcDSS14InvCrossSec2760GeV");
        graphNLOCalcDSS14InvSecPi07000GeV->Write("graphNLOCalcDSS14InvCrossSec7000GeV");
        graphNLOCalcDSS14InvYieldPi02760GeV->Write("graphNLOCalcDSS14InvYield2760GeV");
        graphNLOCalcDSS14InvYieldPi07000GeV->Write("graphNLOCalcDSS14InvYield7000GeV");
        
        
    fileTheoryGraphsPP.Close();

}
