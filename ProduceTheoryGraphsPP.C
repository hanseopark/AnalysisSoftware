/*****************************************************************************
******        provided by Gamma Conversion Group, PWGGA,                ******
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

Double_t xSection900GeV         = 47.78*1e-3;       // V0OR
Double_t xSection2760GeV        = 55.416*1e-3;      // V0OR
// reference for 5 TeV: https://aliceinfo.cern.ch/ArtSubmission/sites/aliceinfo.cern.ch.ArtSubmission/files/draft/mgagliar/2016-Aug-22-paper_draft-vdmNote_5TeV.pdf
Double_t xSection5023GeV        = 51.2*1e-3;        // V0AND - 51.2±1.2 
Double_t xSection7000GeV        = 62.22*1e-3;       // V0OR
Double_t xSection8000GeV        = 55.8*1e-3;        // V0AND
Double_t recalcBarn             = 1e12; //NLO in pbarn!!!!
Double_t xSection2760GeVINEL    = 62.8*1e-3;
Double_t xSection7000GeVINEL    = 73.2*1e-3;

TGraphErrors* ScaleGraph (TGraphErrors* graph, Double_t scaleFac){
    TGraphErrors* dummyGraph    = (TGraphErrors*)graph->Clone(Form("%s_Scaled",graph->GetName()));

    Double_t * xValue           = dummyGraph->GetX(); 
    Double_t * yValue           = dummyGraph->GetY();
    Double_t* xError            = dummyGraph->GetEX();
    Double_t* yError            = dummyGraph->GetEY();
    Int_t nPoints               = dummyGraph->GetN();
    for (Int_t i = 0; i < nPoints; i++){
        yValue[i]               = yValue[i]*scaleFac;
        yError[i]               = yError[i]*scaleFac;
    }
    TGraphErrors* returnGraph =  new TGraphErrors(nPoints,xValue,yValue,xError,yError); 
    return returnGraph;
}


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


TGraphAsymmErrors* CombineMuScales( Int_t nPoints,
                                    Double_t* pt,
                                    Double_t* mu,
                                    Double_t* muHalf,
                                    Double_t* muTwo
                                  ){
    
    
    Double_t yValue[400];
    Double_t yErrorHigh[400];
    Double_t yErrorLow[400];
    Double_t ptErr[400];
    
    for (Int_t i = 0; i< nPoints; i++){
        yValue[i]   = mu[i];
        ptErr[i]    = 0;
        if (muHalf[i] < mu[i] && muHalf[i] != 0){
            yErrorLow[i]    = TMath::Abs(mu[i]-muHalf[i]);
            yErrorHigh[i]   = TMath::Abs(muTwo[i]-mu[i]);
        } else {
            yErrorLow[i]    = TMath::Abs(mu[i]-muTwo[i]);
            yErrorHigh[i]   = TMath::Abs(muHalf[i]-mu[i]);
        }    
    }
    TGraphAsymmErrors* graphReturn = new TGraphAsymmErrors(nPoints, pt, yValue, ptErr, ptErr, yErrorLow, yErrorHigh);
    return graphReturn;
}

//**********************************************************************************************************************
//*********************************** Main function ********************************************************************
//**********************************************************************************************************************
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
//     graphNLOCalcInvSecEtaMuOne900GeV->Print();
    cout << "pi0" << endl;
//     graphNLOCalcInvSecPi0MuOne900GeV->Print();


    
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
    TGraph* graphNLOCalcInvSecEtaMuHalf2760GeV          = new TGraph(nlinesNLOEta2760GeV,ptNLOEta2760GeV,muHalfEta2760GeV); 
    graphNLOCalcInvSecEtaMuHalf2760GeV->RemovePoint(0);
    graphNLOCalcInvSecEtaMuHalf2760GeV->RemovePoint(0);
    TGraph* graphNLOCalcInvSecEtaMuOne2760GeV           = new TGraph(nlinesNLOEta2760GeV,ptNLOEta2760GeV,muOneEta2760GeV); 
    graphNLOCalcInvSecEtaMuOne2760GeV->RemovePoint(0);
    TGraph* graphNLOCalcInvSecEtaMuTwo2760GeV           = new TGraph(nlinesNLOEta2760GeV,ptNLOEta2760GeV,muTwoEta2760GeV); 
    graphNLOCalcInvSecEtaMuTwo2760GeV->RemovePoint(0);
    TGraphAsymmErrors* graphNLOCalcInvSecEta2760GeV     = CombineMuScales(nlinesNLOEta2760GeV, ptNLOEta2760GeV, muOneEta2760GeV, muHalfEta2760GeV, muTwoEta2760GeV);
    graphNLOCalcInvSecEta2760GeV->RemovePoint(0);
    graphNLOCalcInvSecEta2760GeV->RemovePoint(0);
    
    
    TGraph* graphNLOCalcInvYieldEtaMuHalf2760GeV        = ScaleGraph(graphNLOCalcInvSecEtaMuHalf2760GeV, 1/(xSection2760GeV*recalcBarn));
    TGraph* graphNLOCalcInvYieldEtaMuOne2760GeV         = ScaleGraph(graphNLOCalcInvSecEtaMuOne2760GeV, 1/(xSection2760GeV*recalcBarn));
    TGraph* graphNLOCalcInvYieldEtaMuTwo2760GeV         = ScaleGraph(graphNLOCalcInvSecEtaMuTwo2760GeV, 1/(xSection2760GeV*recalcBarn));
    TGraphAsymmErrors* graphNLOCalcInvYieldEta2760GeV   = ScaleGraphAsym(graphNLOCalcInvSecEta2760GeV, 1/(xSection2760GeV*recalcBarn));
    
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

    TGraph* graphNLOCalcInvYieldPi0MuHalf2760GeV        = ScaleGraph(graphNLOCalcInvSecPi0MuHalf2760GeV, 1/(xSection2760GeV*recalcBarn));
    TGraph* graphNLOCalcInvYieldPi0MuOne2760GeV         = ScaleGraph(graphNLOCalcInvSecPi0MuOne2760GeV, 1/(xSection2760GeV*recalcBarn));
    TGraph* graphNLOCalcInvYieldPi0MuTwo2760GeV         = ScaleGraph(graphNLOCalcInvSecPi0MuTwo2760GeV, 1/(xSection2760GeV*recalcBarn));
    
    TGraphAsymmErrors* graphNLOCalcInvSecPi02760GeV     = CombineMuScales(nlinesNLOPi02760GeV, ptNLOPi02760GeV, muOnePi02760GeV, muHalfPi02760GeV, muTwoPi02760GeV);
    graphNLOCalcInvSecPi02760GeV->RemovePoint(0);
    TGraphAsymmErrors* graphNLOCalcInvYieldPi02760GeV   = ScaleGraphAsym(graphNLOCalcInvSecPi02760GeV, 1/(xSection2760GeV*recalcBarn));
    
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
    Double_t    valueNLOEtaToPi0NLOMuHalf2760GeV[100];
    Double_t    valueNLOEtaToPi0NLOMuOne2760GeV[100];
    Double_t    valueNLOEtaToPi0NLOMuTwo2760GeV[100];
    
    for ( Int_t n = 0; n < nlinesNLOEta2760GeV+1; n++){
        if (ptNLOEta2760GeV[n] == ptNLOPi02760GeV[n]){
            //cout << valueNLOMuHalfPi02760GeV[n]  << "\t" << valueNLOMuOnePi02760GeV[n] << "\t" << valueNLOMuTwoPi02760GeV[n]<<endl;
            if (n == 0){
                valueNLOEtaToPi0NLOMuHalf2760GeV[n] = 0.;
            } else { 
                if (muHalfPi02760GeV[n] != 0){
                    valueNLOEtaToPi0NLOMuHalf2760GeV[n] = muHalfEta2760GeV[n]/muHalfPi02760GeV[n];
                } else {
                    valueNLOEtaToPi0NLOMuHalf2760GeV[n] = 0.;
                }
            }
            if (muOnePi02760GeV[n] != 0){
                valueNLOEtaToPi0NLOMuOne2760GeV[n] = muOneEta2760GeV[n]/muOnePi02760GeV[n];
            } else {
                valueNLOEtaToPi0NLOMuOne2760GeV[n] = 0.;
            }
            if (muTwoPi02760GeV[n] != 0){
                valueNLOEtaToPi0NLOMuTwo2760GeV[n] = muTwoEta2760GeV[n]/muTwoPi02760GeV[n];
            } else {
                valueNLOEtaToPi0NLOMuTwo2760GeV[n] = 0.;
            }
        } else {
            cout << "missmatch in pt between eta and pi0: " << ptNLOEta2760GeV[n] << "\t" << ptNLOPi02760GeV[n] << endl;
        }     
    }  
    TGraph* graphEtaToPi0NLOMuHalf2760GeV           = new TGraph(nlinesNLOEta2760GeV,ptNLOEta2760GeV,valueNLOEtaToPi0NLOMuHalf2760GeV);  
    graphEtaToPi0NLOMuHalf2760GeV->RemovePoint(0);
    graphEtaToPi0NLOMuHalf2760GeV->RemovePoint(0);
    TGraph* graphEtaToPi0NLOMuOne2760GeV            = new TGraph(nlinesNLOEta2760GeV,ptNLOEta2760GeV,valueNLOEtaToPi0NLOMuOne2760GeV); 
    graphEtaToPi0NLOMuOne2760GeV->RemovePoint(0);
    TGraph* graphEtaToPi0NLOMuTwo2760GeV            = new TGraph(nlinesNLOEta2760GeV,ptNLOEta2760GeV,valueNLOEtaToPi0NLOMuTwo2760GeV); 
    graphEtaToPi0NLOMuTwo2760GeV->RemovePoint(0);
    
    TGraphAsymmErrors* graphNLOCalcEtaToPi02760GeV  = CombineMuScales(nlinesNLOEta2760GeV, ptNLOEta2760GeV, valueNLOEtaToPi0NLOMuOne2760GeV, valueNLOEtaToPi0NLOMuHalf2760GeV, valueNLOEtaToPi0NLOMuTwo2760GeV);
    graphNLOCalcEtaToPi02760GeV->RemovePoint(0);
    graphNLOCalcEtaToPi02760GeV->RemovePoint(0);
    
    
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
    //******************************** 5.023 TeV Pi0 and Eta calc***************************************************************
    //**********************************************************************************************************************

    //*********************************************************
    // eta Vogelsang - FF: AESSS, PDF: CT10
    Double_t ptNLOEta5023GeV[200];
    Double_t muHalfEta5023GeV[200];
    Double_t muOneEta5023GeV[200];
    Double_t muTwoEta5023GeV[200];
    Int_t nlinesNLOEta5023GeV                       = 0;
    Bool_t fillAllMuScalesEta5023GeV                = kFALSE;
    
    TString fileNameNLOEta5023GeV                   = "ExternalInput/Theory/ALICENLOcalcEtaVogelsang5023GeV_AESSS_CT10.dat"; // currently only mu = 1pt
    ifstream  fileNLOEta5023GeV;
    fileNLOEta5023GeV.open(fileNameNLOEta5023GeV,ios_base::in);
    cout << fileNameNLOEta5023GeV << endl;
    
    // read eta input theory file for 5.023TeV
    while(!fileNLOEta5023GeV.eof() && nlinesNLOEta5023GeV < 200){
        if (fillAllMuScalesEta5023GeV){
            fileNLOEta5023GeV >> ptNLOEta5023GeV[nlinesNLOEta5023GeV] >> muHalfEta5023GeV[nlinesNLOEta5023GeV] >> muOneEta5023GeV[nlinesNLOEta5023GeV] >> muTwoEta5023GeV[nlinesNLOEta5023GeV]; 
        } else {    
            fileNLOEta5023GeV >> ptNLOEta5023GeV[nlinesNLOEta5023GeV] >> muOneEta5023GeV[nlinesNLOEta5023GeV] ; 
            muHalfEta5023GeV[nlinesNLOEta5023GeV]   = muOneEta5023GeV[nlinesNLOEta5023GeV] ; 
            muTwoEta5023GeV[nlinesNLOEta5023GeV]    = muOneEta5023GeV[nlinesNLOEta5023GeV] ; 
        }    
        nlinesNLOEta5023GeV++;
    }
    fileNLOEta5023GeV.close();
    // reduce number of lines by 1
    nlinesNLOEta5023GeV--;
    
    // generate graphs
    TGraph* graphNLOCalcInvSecEtaMuHalf5023GeV      = NULL;
    TGraph* graphNLOCalcInvSecEtaMuOne5023GeV       = NULL;
    TGraph* graphNLOCalcInvSecEtaMuTwo5023GeV       = NULL;
    TGraph* graphNLOCalcInvYieldEtaMuHalf5023GeV    = NULL;
    TGraph* graphNLOCalcInvYieldEtaMuOne5023GeV     = NULL;
    TGraph* graphNLOCalcInvYieldEtaMuTwo5023GeV     = NULL;

    // fill x-section graph for default mu
    graphNLOCalcInvSecEtaMuOne5023GeV               = new TGraph(nlinesNLOEta5023GeV,ptNLOEta5023GeV,muOneEta5023GeV); 
    graphNLOCalcInvSecEtaMuOne5023GeV->Print();
    // fill yield graph for default mu
    graphNLOCalcInvYieldEtaMuOne5023GeV             = ScaleGraph(graphNLOCalcInvSecEtaMuOne5023GeV, 1/(xSection5023GeV*recalcBarn));
    
    if (fillAllMuScalesEta5023GeV){
        // fill x-section graphs for mu variations
        graphNLOCalcInvSecEtaMuHalf5023GeV          = new TGraph(nlinesNLOEta5023GeV,ptNLOEta5023GeV,muHalfEta5023GeV); 
        graphNLOCalcInvSecEtaMuHalf5023GeV->RemovePoint(0);
        graphNLOCalcInvSecEtaMuTwo5023GeV           = new TGraph(nlinesNLOEta5023GeV,ptNLOEta5023GeV,muTwoEta5023GeV); 
        // fill yield graphs for mu variations
        graphNLOCalcInvYieldEtaMuHalf5023GeV        = ScaleGraph(graphNLOCalcInvSecEtaMuHalf5023GeV, 1/(xSection5023GeV*recalcBarn));
        graphNLOCalcInvYieldEtaMuTwo5023GeV         = ScaleGraph(graphNLOCalcInvSecEtaMuTwo5023GeV, 1/(xSection5023GeV*recalcBarn));
    }
    // combine all mu scales into 1 graph for eventual plotting
    TGraphAsymmErrors* graphNLOCalcInvSecEta5023GeV     = CombineMuScales(nlinesNLOEta5023GeV, ptNLOEta5023GeV, muOneEta5023GeV, muHalfEta5023GeV, muTwoEta5023GeV);
    TGraphAsymmErrors* graphNLOCalcInvYieldEta5023GeV   = ScaleGraphAsym(graphNLOCalcInvSecEta5023GeV, 1/(xSection5023GeV*recalcBarn));

    //*********************************************************
    // pi0 Vogelsang - FF: DSS14, PDF: CT10
    Double_t ptNLOPi05023GeV[200];
    Double_t muHalfPi05023GeV[200];
    Double_t muOnePi05023GeV[200];
    Double_t muTwoPi05023GeV[200];
    Int_t nlinesNLOPi05023GeV                       = 0;
    
    TString fileNameNLOPi05023GeV                   = "ExternalInput/Theory/ALICENLOCalcPi0Vogelsang5023GeV_DSS14_CT10.dat";
    ifstream  fileNLOPi05023GeV;
    fileNLOPi05023GeV.open(fileNameNLOPi05023GeV,ios_base::in);
    cout << fileNameNLOPi05023GeV << endl;
    
    while(!fileNLOPi05023GeV.eof() && nlinesNLOPi05023GeV < 200){
        fileNLOPi05023GeV >> ptNLOPi05023GeV[nlinesNLOPi05023GeV] >> muHalfPi05023GeV[nlinesNLOPi05023GeV] >> muOnePi05023GeV[nlinesNLOPi05023GeV] >> muTwoPi05023GeV[nlinesNLOPi05023GeV]; 
        cout << nlinesNLOPi05023GeV << "         "  << ptNLOPi05023GeV[nlinesNLOPi05023GeV] << "         "  << muHalfPi05023GeV[nlinesNLOPi05023GeV] << "         "  << muOnePi05023GeV[nlinesNLOPi05023GeV] << "         "  << muTwoPi05023GeV[nlinesNLOPi05023GeV] << endl;
        nlinesNLOPi05023GeV++;
    }
    fileNLOPi05023GeV.close();
    // reduce number of lines by 1
    nlinesNLOPi05023GeV--;
    
    // fill x-section graphs
    TGraph* graphNLOCalcInvSecPi0MuHalf5023GeV      = new TGraph(nlinesNLOPi05023GeV,ptNLOPi05023GeV,muHalfPi05023GeV); 
    graphNLOCalcInvSecPi0MuHalf5023GeV->RemovePoint(0);
    TGraph* graphNLOCalcInvSecPi0MuOne5023GeV       = new TGraph(nlinesNLOPi05023GeV,ptNLOPi05023GeV,muOnePi05023GeV); 
    TGraph* graphNLOCalcInvSecPi0MuTwo5023GeV       = new TGraph(nlinesNLOPi05023GeV,ptNLOPi05023GeV,muTwoPi05023GeV); 
    // fill yield graphs
    TGraph* graphNLOCalcInvYieldPi0MuHalf5023GeV    = ScaleGraph(graphNLOCalcInvSecPi0MuHalf5023GeV, 1/(xSection5023GeV*recalcBarn));
    TGraph* graphNLOCalcInvYieldPi0MuOne5023GeV     = ScaleGraph(graphNLOCalcInvSecPi0MuOne5023GeV, 1/(xSection5023GeV*recalcBarn));
    TGraph* graphNLOCalcInvYieldPi0MuTwo5023GeV     = ScaleGraph(graphNLOCalcInvSecPi0MuTwo5023GeV, 1/(xSection5023GeV*recalcBarn));
    // combine all mu scales into 1 graph for eventual plotting
    TGraphAsymmErrors* graphNLOCalcInvSecPi05023GeV     = CombineMuScales(nlinesNLOPi05023GeV, ptNLOPi05023GeV, muOnePi05023GeV, muHalfPi05023GeV, muTwoPi05023GeV);
    TGraphAsymmErrors* graphNLOCalcInvYieldPi05023GeV   = ScaleGraphAsym(graphNLOCalcInvSecPi05023GeV, 1/(xSection5023GeV*recalcBarn));

    //*********************************************************
    // Eta/Pi0 Vogelsang FF pi0: DSS14, FF eta: AESSS, PDF: CT10
    Double_t* valueNLOMuHalfEta5023GeV              = NULL;
    Double_t* valueNLOMuOneEta5023GeV               = graphNLOCalcInvSecEtaMuOne5023GeV->GetY();
    Double_t* valueNLOMuTwoEta5023GeV               = NULL;
    if (fillAllMuScalesEta5023GeV){
        valueNLOMuHalfEta5023GeV                    = graphNLOCalcInvSecEtaMuHalf5023GeV->GetY();
        valueNLOMuTwoEta5023GeV                     = graphNLOCalcInvSecEtaMuTwo5023GeV->GetY();
    }
    Double_t* valueNLOMuHalfPi05023GeV              = graphNLOCalcInvSecPi0MuHalf5023GeV->GetY();
    Double_t* valueNLOMuOnePi05023GeV               = graphNLOCalcInvSecPi0MuOne5023GeV->GetY();
    Double_t* valueNLOMuTwoPi05023GeV               = graphNLOCalcInvSecPi0MuTwo5023GeV->GetY();
    Double_t* xValueNLO5023GeV                      = graphNLOCalcInvSecPi0MuOne5023GeV->GetX();
    Int_t xNBins5023GeV                             = graphNLOCalcInvSecPi0MuOne5023GeV->GetN();
    Double_t valueNLOEtaToPi0NLOMuHalf5023GeV[200];
    Double_t valueNLOEtaToPi0NLOMuOne5023GeV[200];
    Double_t valueNLOEtaToPi0NLOMuTwo5023GeV[200];
    
    // calculate eta/pi0 ratios for different mu scales
    for ( Int_t n = 0; n < xNBins5023GeV+1; n++){
        if (valueNLOMuOnePi05023GeV[n] != 0){
            valueNLOEtaToPi0NLOMuOne5023GeV[n] = valueNLOMuOneEta5023GeV[n]/valueNLOMuOnePi05023GeV[n];
        } else {
            valueNLOEtaToPi0NLOMuOne5023GeV[n] = 0.;
        }

        if (fillAllMuScalesEta5023GeV){
            if (n == 0){
                valueNLOEtaToPi0NLOMuHalf5023GeV[n] = 0.;
            } else { 
                if (valueNLOMuHalfPi05023GeV[n] != 0){
                    valueNLOEtaToPi0NLOMuHalf5023GeV[n] = valueNLOMuHalfEta5023GeV[n]/valueNLOMuHalfPi05023GeV[n];
                } else {
                    valueNLOEtaToPi0NLOMuHalf5023GeV[n] = 0.;
                }
            }
            if (valueNLOMuTwoPi05023GeV[n] != 0){
                valueNLOEtaToPi0NLOMuTwo5023GeV[n] = valueNLOMuTwoEta5023GeV[n]/valueNLOMuTwoPi05023GeV[n];
            } else {
                valueNLOEtaToPi0NLOMuTwo5023GeV[n] = 0.;
            }
        } else {
            valueNLOEtaToPi0NLOMuHalf5023GeV[n] = valueNLOEtaToPi0NLOMuOne5023GeV[n];
            valueNLOEtaToPi0NLOMuTwo5023GeV[n] = valueNLOEtaToPi0NLOMuOne5023GeV[n];
        }    
    }  
    // fill graphs
    TGraph* graphEtaToPi0NLOMuHalf5023GeV           = NULL;
    TGraph* graphEtaToPi0NLOMuTwo5023GeV            = NULL;
    TGraph* graphEtaToPi0NLOMuOne5023GeV            = new TGraph(xNBins5023GeV,xValueNLO5023GeV,valueNLOEtaToPi0NLOMuOne5023GeV);   
    if (fillAllMuScalesEta5023GeV){
        graphEtaToPi0NLOMuHalf5023GeV               = new TGraph(xNBins5023GeV,xValueNLO5023GeV,valueNLOEtaToPi0NLOMuHalf5023GeV);  
        graphEtaToPi0NLOMuHalf5023GeV->RemovePoint(0);
        graphEtaToPi0NLOMuTwo5023GeV                = new TGraph(xNBins5023GeV,xValueNLO5023GeV,valueNLOEtaToPi0NLOMuTwo5023GeV);   
    }    
    // combine all mu scales into 1 graph for plotting
    TGraphAsymmErrors* graphNLOCalcEtaToPi05023GeV  = CombineMuScales(xNBins5023GeV, xValueNLO5023GeV, valueNLOEtaToPi0NLOMuOne5023GeV, valueNLOEtaToPi0NLOMuHalf5023GeV, valueNLOEtaToPi0NLOMuTwo5023GeV);
    
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
    // Stratmann
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
    
    //**************************************************************************************************
    //********************** Pythia 6 Perugia2011 MC spectra and ratio *********************************
    //**************************************************************************************************
    // file generated with TaskV1/ExtractMCInputSpectraFromFile.C++ based on PCM only inputs
    TFile* file7TeVPyth6Perugia                 = new TFile("ExternalInput/Theory/MCInputCompilationLHC14j4[b-f]_pp7TeV_0.root");
    TH1D* histoPi07TeVPythia6                   = (TH1D*)file7TeVPyth6Perugia->Get("MC_Pi0_Pt");
    TH1D* histoPi07TeVPythia6Reb                = (TH1D*)file7TeVPyth6Perugia->Get("MC_Pi0_Pt_Rebinned");
    TH1D* histoEta7TeVPythia6                   = (TH1D*)file7TeVPyth6Perugia->Get("MC_Eta_Pt");
    TH1D* histoEta7TeVPythia6Reb                = (TH1D*)file7TeVPyth6Perugia->Get("MC_Eta_Pt_Rebinned");
    TH1D* histoEtaToPi07TeVPythia6              = (TH1D*)file7TeVPyth6Perugia->Get("MCEtaToPi0");
    
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

    TGraphAsymmErrors* graphNLOCalcInvSecEta8000GeV     = CombineMuScales(nlinesNLOEta8000GeV, ptNLOEta8000GeV, muOneEta8000GeV, muHalfEta8000GeV, muTwoEta8000GeV);
    graphNLOCalcInvSecEta8000GeV->RemovePoint(0);
    
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
//     graphEtaToPi0NLOMuHalf8TeV->Print();
    TGraph* graphEtaToPi0NLOMuOne8TeV       = new TGraph(xNBins8000GeV,xValueNLO8000GeV,valueNLOEtaToPi0NLOMuOne8000GeV);   
//     graphEtaToPi0NLOMuOne8TeV->Print();
    TGraph* graphEtaToPi0NLOMuTwo8TeV       = new TGraph(xNBins8000GeV,xValueNLO8000GeV,valueNLOEtaToPi0NLOMuTwo8000GeV);   
//     graphEtaToPi0NLOMuTwo8TeV->Print();

    TGraphAsymmErrors* graphNLOCalcEtaToPi08000GeV  = CombineMuScales(xNBins8000GeV, xValueNLO8000GeV, graphEtaToPi0NLOMuOne8TeV->GetY(), graphEtaToPi0NLOMuHalf8TeV->GetY(), graphEtaToPi0NLOMuTwo8TeV->GetY());

    //**********************************************************************************************************************
    //******************************** 8TeV Pi0 NLO DSS14 calc *************************************************************
    //**********************************************************************************************************************
    
    TString fileNameTheoryPi0DSS14                    = "ExternalInput/Theory/pp8TeV_NLO_Pi0_DSS14.root";
    TFile* fileTheorypp8TeVPi0DSS14                   = new TFile(fileNameTheoryPi0DSS14.Data());
    TGraphAsymmErrors* graphNLOCalcInvSecPi08000GeV   = (TGraphAsymmErrors*)fileTheorypp8TeVPi0DSS14->Get("fGraphInvXsec_Pi0");
    //scale errors 1/2 mu to 2 mu are given by errors of TGraphAsymmErrors

    //fix normalization problem since factor 2 was missing:
    for (int i=0;i<graphNLOCalcInvSecPi08000GeV->GetN();i++) graphNLOCalcInvSecPi08000GeV->GetY()[i] /= 2.;

    //**************************************************************************************************
    //********************** 8TeV Pythia 8 MC spectra and ratio ****************************************
    //**************************************************************************************************
    // file generated with TaskV1/ExtractMCInputSpectraFromFile.C++ based on PCM only inputs
    TFile* file8TeVPyth8                        = new TFile("ExternalInput/Theory/MCInputCompilationLHC15h1_pp8TeV_0.root");
    TH1D* histoPi08TeVPythia8                   = (TH1D*)file8TeVPyth8->Get("MC_Pi0_Pt");
    TH1D* histoPi08TeVPythia8Reb                = (TH1D*)file8TeVPyth8->Get("MC_Pi0_Pt_Rebinned");
    TH1D* histoEta8TeVPythia8                   = (TH1D*)file8TeVPyth8->Get("MC_Eta_Pt");
    TH1D* histoEta8TeVPythia8Reb                = (TH1D*)file8TeVPyth8->Get("MC_Eta_Pt_Rebinned");
    TH1D* histoEtaToPi08TeVPythia8              = (TH1D*)file8TeVPyth8->Get("MCEtaToPi0");
    
    //**************************************************************************************************
    //********************** 8TeV Phojet MC spectra and ratio ******************************************
    //**************************************************************************************************
    // file generated with TaskV1/ExtractMCInputSpectraFromFile.C++ based on PCM only inputs
    TFile* file8TeVPhojet                       = new TFile("ExternalInput/Theory/MCInputCompilationLHC15h2_pp8TeV_0.root");
    TH1D* histoPi08TeVPhojet                    = (TH1D*)file8TeVPhojet->Get("MC_Pi0_Pt");
    TH1D* histoPi08TeVPhojetReb                 = (TH1D*)file8TeVPhojet->Get("MC_Pi0_Pt_Rebinned");
    TH1D* histoEta8TeVPhojet                    = (TH1D*)file8TeVPhojet->Get("MC_Eta_Pt");
    TH1D* histoEta8TeVPhojetReb                 = (TH1D*)file8TeVPhojet->Get("MC_Eta_Pt_Rebinned");
    TH1D* histoEtaToPi08TeVPhojet               = (TH1D*)file8TeVPhojet->Get("MCEtaToPi0");
    

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
    
//     TFile* filePythia8Monash2013_2760GeV        = TFile::Open("ExternalInput/Theory/Pythia/PYTHIA8_Monash2013Tune_2760GeV.root");
//     TH1F* histoPi0Pythia8MonashInvSec2760GeV    = (TH1F*)filePythia8Monash2013_2760GeV->Get("fHistInvXsec_Pi0");
//     TH1F* histoEtaPythia8MonashInvSec2760GeV    = (TH1F*)filePythia8Monash2013_2760GeV->Get("fHistInvXsec_Eta");
//     TH1F* histoEtaToPi0RatioPythia8Monash2760GeV= (TH1F*)histoEtaPythia8MonashInvSec2760GeV->Clone("histoEtaToPi0RatioPythia8Monash2760GeV");
//     histoEtaToPi0RatioPythia8Monash2760GeV->Sumw2();
//     histoEtaToPi0RatioPythia8Monash2760GeV->Divide(histoEtaToPi0RatioPythia8Monash2760GeV,histoPi0Pythia8MonashInvSec2760GeV);

    TFile* filePythia8Monash2013_2760GeVLego        = TFile::Open("ExternalInput/Theory/Pythia/Pythia8_Monash2013_2760GeV_11120Mio.root");
    TH1F* histoPi0Pythia8MonashInvSec2760GeVLego    = (TH1F*)filePythia8Monash2013_2760GeVLego->Get("hPt_Pi0_MB_XSec");
    TH1F* histoEtaPythia8MonashInvSec2760GeVLego    = (TH1F*)filePythia8Monash2013_2760GeVLego->Get("hPt_Eta_MB_XSec");
    TH1F* histoEtaToPi0RatioPythia8Monash2760GeVLego= (TH1F*)histoEtaPythia8MonashInvSec2760GeVLego->Clone("histoEtaToPi0RatioPythia8Monash2760GeVLego");
    histoEtaToPi0RatioPythia8Monash2760GeVLego->Sumw2();
    histoEtaToPi0RatioPythia8Monash2760GeVLego->Divide(histoEtaToPi0RatioPythia8Monash2760GeVLego,histoPi0Pythia8MonashInvSec2760GeVLego);
    
    //**********************************************************************************************************************
    //***************************** Pythia calculations 0.9TeV ************************************************************
    //**********************************************************************************************************************    
//     TFile* filePythia8Monash2013_900GeV        = TFile::Open("ExternalInput/Theory/Pythia/PYTHIA8_Monash2013Tune_pp900GeV.root");
//     TH1F* histoPi0Pythia8MonashInvSec900GeV    = (TH1F*)filePythia8Monash2013_900GeV->Get("fHistInvXsec_Pi0");
//     TH1F* histoEtaPythia8MonashInvSec900GeV    = (TH1F*)filePythia8Monash2013_900GeV->Get("fHistInvXsec_Eta");
//     TH1F* histoEtaToPi0RatioPythia8Monash900GeV= (TH1F*)histoEtaPythia8MonashInvSec900GeV->Clone("histoEtaToPi0RatioPythia8Monash900GeV");
//     histoEtaToPi0RatioPythia8Monash900GeV->Sumw2();
//     histoEtaToPi0RatioPythia8Monash900GeV->Divide(histoEtaToPi0RatioPythia8Monash900GeV,histoPi0Pythia8MonashInvSec900GeV);

    TFile* filePythia8Monash2013_900GeVLego        = TFile::Open("ExternalInput/Theory/Pythia/Pythia8_Monash2013_900GeV_1598Mio.root");
    TH1F* histoPi0Pythia8MonashInvSec900GeVLego    = (TH1F*)filePythia8Monash2013_900GeVLego->Get("hPt_Pi0_MB_XSec");
    TH1F* histoEtaPythia8MonashInvSec900GeVLego    = (TH1F*)filePythia8Monash2013_900GeVLego->Get("hPt_Eta_MB_XSec");
    TH1F* histoEtaToPi0RatioPythia8Monash900GeVLego= (TH1F*)histoEtaPythia8MonashInvSec900GeVLego->Clone("histoEtaToPi0RatioPythia8Monash900GeVLego");
    histoEtaToPi0RatioPythia8Monash900GeVLego->Sumw2();
    histoEtaToPi0RatioPythia8Monash900GeVLego->Divide(histoEtaToPi0RatioPythia8Monash900GeVLego,histoPi0Pythia8MonashInvSec900GeVLego);
    
    //**********************************************************************************************************************
    //***************************** Pythia calculations 5.023TeV ************************************************************
    //**********************************************************************************************************************    
//     TFile* filePythia8Monash2013_5023GeV        = TFile::Open("ExternalInput/Theory/Pythia/PYTHIA8_Monash2013Tune_pp5023GeV.root");
//     TH1F* histoPi0Pythia8MonashInvSec5023GeV    = (TH1F*)filePythia8Monash2013_5023GeV->Get("fHistInvXsec_Pi0");
//     TH1F* histoEtaPythia8MonashInvSec5023GeV    = (TH1F*)filePythia8Monash2013_5023GeV->Get("fHistInvXsec_Eta");
//     TH1F* histoEtaToPi0RatioPythia8Monash5023GeV= (TH1F*)histoEtaPythia8MonashInvSec5023GeV->Clone("histoEtaToPi0RatioPythia8Monash5023GeV");
//     histoEtaToPi0RatioPythia8Monash5023GeV->Sumw2();
//     histoEtaToPi0RatioPythia8Monash5023GeV->Divide(histoEtaToPi0RatioPythia8Monash5023GeV,histoPi0Pythia8MonashInvSec5023GeV);

    TFile* filePythia8Monash2013_5TeVLego        = TFile::Open("ExternalInput/Theory/Pythia/Pythia8_Monash2013_5TeV_1502Mio.root");
    TH1F* histoPi0Pythia8MonashInvSec5TeVLego    = (TH1F*)filePythia8Monash2013_5TeVLego->Get("hPt_Pi0_MB_XSec");
    TH1F* histoEtaPythia8MonashInvSec5TeVLego    = (TH1F*)filePythia8Monash2013_5TeVLego->Get("hPt_Eta_MB_XSec");
    TH1F* histoEtaToPi0RatioPythia8Monash5TeVLego= (TH1F*)histoEtaPythia8MonashInvSec5TeVLego->Clone("histoEtaToPi0RatioPythia8Monash5TeVLego");
    histoEtaToPi0RatioPythia8Monash5TeVLego->Sumw2();
    histoEtaToPi0RatioPythia8Monash5TeVLego->Divide(histoEtaToPi0RatioPythia8Monash5TeVLego,histoPi0Pythia8MonashInvSec5TeVLego);
    
    //**********************************************************************************************************************
    //***************************** Pythia calculations 7TeV ************************************************************
    //**********************************************************************************************************************    
//     TFile* filePythia8Monash2013_7TeV           = TFile::Open("ExternalInput/Theory/Pythia/PYTHIA8_Monash2013Tune_pp7000GeV.root");
//     TH1F* histoPi0Pythia8MonashInvSec7TeV       = (TH1F*)filePythia8Monash2013_7TeV->Get("fHistInvXsec_Pi0");
//     TH1F* histoEtaPythia8MonashInvSec7TeV       = (TH1F*)filePythia8Monash2013_7TeV->Get("fHistInvXsec_Eta");
//     TH1F* histoEtaToPi0RatioPythia8Monash7TeV   = (TH1F*)histoEtaPythia8MonashInvSec7TeV->Clone("histoEtaToPi0RatioPythia8Monash7TeV");
//     histoEtaToPi0RatioPythia8Monash7TeV->Sumw2();
//     histoEtaToPi0RatioPythia8Monash7TeV->Divide(histoEtaToPi0RatioPythia8Monash7TeV,histoPi0Pythia8MonashInvSec7TeV);

    TFile* filePythia8Monash2013_7TeVLego        = TFile::Open("ExternalInput/Theory/Pythia/Pythia8_Monash2013_7TeV_1650Mio.root");
    TH1F* histoPi0Pythia8MonashInvSec7TeVLego    = (TH1F*)filePythia8Monash2013_7TeVLego->Get("hPt_Pi0_MB_XSec");
    TH1F* histoEtaPythia8MonashInvSec7TeVLego    = (TH1F*)filePythia8Monash2013_7TeVLego->Get("hPt_Eta_MB_XSec");
    TH1F* histoEtaToPi0RatioPythia8Monash7TeVLego= (TH1F*)histoEtaPythia8MonashInvSec7TeVLego->Clone("histoEtaToPi0RatioPythia8Monash7TeVLego");
    histoEtaToPi0RatioPythia8Monash7TeVLego->Sumw2();
    histoEtaToPi0RatioPythia8Monash7TeVLego->Divide(histoEtaToPi0RatioPythia8Monash7TeVLego,histoPi0Pythia8MonashInvSec7TeVLego);

    
    //**********************************************************************************************************************
    //***************************** Pythia calculations 8TeV ************************************************************
    //**********************************************************************************************************************    
    TFile* filePythia8Monash2013_8TeV           = TFile::Open("ExternalInput/Theory/Pythia/PYTHIA8_Monash2013Tune_pp8000GeV.root");
    TH1F* histoPi0Pythia8MonashInvSec8TeV       = (TH1F*)filePythia8Monash2013_8TeV->Get("fHistInvXsec_Pi0");
    TH1F* histoEtaPythia8MonashInvSec8TeV       = (TH1F*)filePythia8Monash2013_8TeV->Get("fHistInvXsec_Eta");
    TH1F* histoEtaToPi0RatioPythia8Monash8TeV   = (TH1F*)histoEtaPythia8MonashInvSec8TeV->Clone("histoEtaToPi0RatioPythia8Monash8TeV");
    histoEtaToPi0RatioPythia8Monash8TeV->Sumw2();
    histoEtaToPi0RatioPythia8Monash8TeV->Divide(histoEtaToPi0RatioPythia8Monash8TeV,histoPi0Pythia8MonashInvSec8TeV);

    TFile* filePythia8Monash2013_8TeVLego        = TFile::Open("ExternalInput/Theory/Pythia/Pythia8_Monash2013_8TeV_3640Mio.root");
    TH1F* histoPi0Pythia8MonashInvSec8TeVLego    = (TH1F*)filePythia8Monash2013_8TeVLego->Get("hPt_Pi0_MB_XSec");
    TH1F* histoEtaPythia8MonashInvSec8TeVLego    = (TH1F*)filePythia8Monash2013_8TeVLego->Get("hPt_Eta_MB_XSec");
    TH1F* histoEtaToPi0RatioPythia8Monash8TeVLego= (TH1F*)histoEtaPythia8MonashInvSec8TeVLego->Clone("histoEtaToPi0RatioPythia8Monash8TeVLego");
    histoEtaToPi0RatioPythia8Monash8TeVLego->Sumw2();
    histoEtaToPi0RatioPythia8Monash8TeVLego->Divide(histoEtaToPi0RatioPythia8Monash8TeVLego,histoPi0Pythia8MonashInvSec8TeVLego);

    TFile* filePythia8Tune4C_8TeVLego            = TFile::Open("ExternalInput/Theory/Pythia/Pythia8_Tune4C_8TeV_3508Mio.root");
    TH1F* histoPi0Pythia8Tune4CInvSec8TeVLego    = (TH1F*)filePythia8Tune4C_8TeVLego->Get("hPt_Pi0_MB_XSec");
    TH1F* histoEtaPythia8Tune4CInvSec8TeVLego    = (TH1F*)filePythia8Tune4C_8TeVLego->Get("hPt_Eta_MB_XSec");
    TH1F* histoEtaToPi0RatioPythia8Tune4C8TeVLego= (TH1F*)histoEtaPythia8Tune4CInvSec8TeVLego->Clone("histoEtaToPi0RatioPythia8Tune4C8TeVLego");
    histoEtaToPi0RatioPythia8Tune4C8TeVLego->Sumw2();
    histoEtaToPi0RatioPythia8Tune4C8TeVLego->Divide(histoEtaToPi0RatioPythia8Tune4C8TeVLego,histoPi0Pythia8Tune4CInvSec8TeVLego);
    
    
    //**********************************************************************************************************************
    //******************************** Plot Pythia x-section at different energies *****************************************
    //**********************************************************************************************************************
    TCanvas* canvasXSection  = new TCanvas("canvasXSection","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasXSection, 0.14, 0.02, 0.02, 0.09);
    canvasXSection->SetLogx();
    canvasXSection->SetLogy();
    
    TH2F * histo2DXSection;
    histo2DXSection          = new TH2F("histo2DXSection","histo2DXSection",11000,0.23,60.,1000,0.01,9e11);
    SetStyleHistoTH2ForGraphs(histo2DXSection, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",0.035,0.04, 0.035,0.04, 0.9,1.45);
    histo2DXSection->GetXaxis()->SetMoreLogLabels();
    histo2DXSection->GetXaxis()->SetNoExponent(kTRUE);
//    histo2DXSection->GetXaxis()->SetLabelOffset(-0.01);
    histo2DXSection->Draw("copy");

        DrawGammaSetMarker(histoPi0Pythia8MonashInvSec8TeVLego, 1, 1.5, GetColorDefaultColor("8TeV","","") , GetColorDefaultColor("8TeV","",""));
        histoPi0Pythia8MonashInvSec8TeVLego->SetLineWidth(3);
        histoPi0Pythia8MonashInvSec8TeVLego->Draw("same,hist,l");

        DrawGammaSetMarker(histoPi0Pythia8MonashInvSec7TeVLego, 1, 1.5, GetColorDefaultColor("7TeV","","") , GetColorDefaultColor("7TeV","",""));
        histoPi0Pythia8MonashInvSec7TeVLego->SetLineWidth(3);
        histoPi0Pythia8MonashInvSec7TeVLego->Draw("same,hist,l");

        DrawGammaSetMarker(histoPi0Pythia8MonashInvSec5TeVLego, 1, 1.5, GetColorDefaultColor("5TeV","","") , GetColorDefaultColor("5TeV","",""));
        histoPi0Pythia8MonashInvSec5TeVLego->SetLineWidth(3);
        histoPi0Pythia8MonashInvSec5TeVLego->Draw("same,hist,l");
        
        DrawGammaSetMarker(histoPi0Pythia8MonashInvSec2760GeVLego, 1, 1.5, GetColorDefaultColor("2.76TeV","","") , GetColorDefaultColor("2.76TeV","",""));
        histoPi0Pythia8MonashInvSec2760GeVLego->SetLineWidth(3);
        histoPi0Pythia8MonashInvSec2760GeVLego->Draw("same,hist,l");

        DrawGammaSetMarker(histoPi0Pythia8MonashInvSec900GeVLego, 1, 1.5, GetColorDefaultColor("900GeV","","") , GetColorDefaultColor("900GeV","",""));
        histoPi0Pythia8MonashInvSec900GeVLego->SetLineWidth(3);
        histoPi0Pythia8MonashInvSec900GeVLego->Draw("same,hist,l");

        TLatex *labelDetSysXSectionGen      = new TLatex(0.94,0.92,"Pythia 8.2, Monash 2013");
        SetStyleTLatex( labelDetSysXSectionGen, 0.035,4, 1, 42, kTRUE, 31);
        labelDetSysXSectionGen->Draw();
        
        TLatex *labelDetSysXSectionPi0      = new TLatex(0.94,0.88,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysXSectionPi0, 0.035,4, 1, 42, kTRUE, 31);
        labelDetSysXSectionPi0->Draw();

        TLegend* legendXSectionPi0DiffFits          = GetAndSetLegend2(0.18, 0.13, 0.40 , 0.13+0.03*5, 32); 
        legendXSectionPi0DiffFits->AddEntry(histoPi0Pythia8MonashInvSec8TeVLego,"pp, #sqrt{s} = 8 TeV","l");
        legendXSectionPi0DiffFits->AddEntry(histoPi0Pythia8MonashInvSec7TeVLego,"pp, #sqrt{s} = 7 TeV","l");
        legendXSectionPi0DiffFits->AddEntry(histoPi0Pythia8MonashInvSec5TeVLego,"pp, #sqrt{s} = 5 TeV","l");
        legendXSectionPi0DiffFits->AddEntry(histoPi0Pythia8MonashInvSec2760GeVLego,"pp, #sqrt{s} = 2.76 TeV","l");
        legendXSectionPi0DiffFits->AddEntry(histoPi0Pythia8MonashInvSec900GeVLego,"pp, #sqrt{s} = 0.9 TeV","l");
        legendXSectionPi0DiffFits->Draw();
        
    canvasXSection->SaveAs("ExternalInput/Theory/Pi0_Pythia8_Monash_InvXSection_comp.eps");
    
    histo2DXSection->Draw("copy");
        DrawGammaSetMarker(histoEtaPythia8MonashInvSec8TeVLego, 1, 1.5, GetColorDefaultColor("8TeV","","") , GetColorDefaultColor("8TeV","",""));
        histoEtaPythia8MonashInvSec8TeVLego->SetLineWidth(3);
        histoEtaPythia8MonashInvSec8TeVLego->Draw("same,hist,l");

        DrawGammaSetMarker(histoEtaPythia8MonashInvSec7TeVLego, 1, 1.5, GetColorDefaultColor("7TeV","","") , GetColorDefaultColor("7TeV","",""));
        histoEtaPythia8MonashInvSec7TeVLego->SetLineWidth(3);
        histoEtaPythia8MonashInvSec7TeVLego->Draw("same,hist,l");

        DrawGammaSetMarker(histoEtaPythia8MonashInvSec5TeVLego, 1, 1.5, GetColorDefaultColor("5TeV","","") , GetColorDefaultColor("5TeV","",""));
        histoEtaPythia8MonashInvSec5TeVLego->SetLineWidth(3);
        histoEtaPythia8MonashInvSec5TeVLego->Draw("same,hist,l");
        
        DrawGammaSetMarker(histoEtaPythia8MonashInvSec2760GeVLego, 1, 1.5, GetColorDefaultColor("2.76TeV","","") , GetColorDefaultColor("2.76TeV","",""));
        histoEtaPythia8MonashInvSec2760GeVLego->SetLineWidth(3);
        histoEtaPythia8MonashInvSec2760GeVLego->Draw("same,hist,l");

        DrawGammaSetMarker(histoEtaPythia8MonashInvSec900GeVLego, 1, 1.5, GetColorDefaultColor("900GeV","","") , GetColorDefaultColor("900GeV","",""));
        histoEtaPythia8MonashInvSec900GeVLego->SetLineWidth(3);
        histoEtaPythia8MonashInvSec900GeVLego->Draw("same,hist,l");

        labelDetSysXSectionGen->Draw();
        
        TLatex *labelDetSysXSectionEta      = new TLatex(0.94,0.88,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysXSectionEta, 0.035,4, 1, 42, kTRUE, 31);
        labelDetSysXSectionEta->Draw();

        legendXSectionPi0DiffFits->Draw();
        
    canvasXSection->SaveAs("ExternalInput/Theory/Eta_Pythia8_Monash_InvXSection_comp.eps");
    
    //**********************************************************************************************************************
    //********************************* Write graphs and histos to compilation file for pp *********************************
    //**********************************************************************************************************************
    TFile fileTheoryGraphsPP("ExternalInput/Theory/TheoryCompilationPP.root","UPDATE");

        //***********************************************************************
        // write  calculations for 900 GeV
        //***********************************************************************
        // pi0 Vogelsang PDF: CT10, FF DSS07
        graphNLOCalcInvSecPi0MuHalf900GeV->Write("graphNLOCalcInvSecPi0MuHalf900GeV", TObject::kOverwrite);
        graphNLOCalcInvSecPi0MuOne900GeV->Write("graphNLOCalcInvSecPi0MuOne900GeV", TObject::kOverwrite);
        graphNLOCalcInvSecPi0MuTwo900GeV->Write("graphNLOCalcInvSecPi0MuTwo900GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldPi0MuHalf900GeV->Write("graphNLOCalcInvYieldPi0MuHalf900GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldPi0MuOne900GeV->Write("graphNLOCalcInvYieldPi0MuOne900GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldPi0MuTwo900GeV->Write("graphNLOCalcInvYieldPi0MuTwo900GeV", TObject::kOverwrite);
        // pi0 INCNLO FF:BKK, DSS07
        graphNLOCalcBKKInvSecPi0MuTwo900GeV->Write("graphNLOCalcBKKInvSecPi0MuTwo900GeV", TObject::kOverwrite);
        graphNLOCalcBKKInvYieldPi0MuTwo900GeV->Write("graphNLOCalcBKKInvYieldPi0MuTwo900GeV", TObject::kOverwrite);
        graphNLOCalcDSSInvSecPi0MuTwo900GeV->Write("graphNLOCalcDSSInvSecPi0MuTwo900GeV", TObject::kOverwrite);
        graphNLOCalcDSSInvYieldPi0MuTwo900GeV->Write("graphNLOCalcDSSInvYieldPi0MuTwo900GeV", TObject::kOverwrite);
        // eta Vogelsang PDF: CT10, FF AESSS
        graphNLOCalcInvSecEtaMuHalf900GeV->Write("graphNLOCalcInvSecEtaMuHalf900GeV", TObject::kOverwrite);
        graphNLOCalcInvSecEtaMuOne900GeV->Write("graphNLOCalcInvSecEtaMuOne900GeV", TObject::kOverwrite);
        graphNLOCalcInvSecEtaMuTwo900GeV->Write("graphNLOCalcInvSecEtaMuTwo900GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldEtaMuHalf900GeV->Write("graphNLOCalcInvYieldEtaMuHalf900GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldEtaMuOne900GeV->Write("graphNLOCalcInvYieldEtaMuOne900GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldEtaMuTwo900GeV->Write("graphNLOCalcInvYieldEtaMuTwo900GeV", TObject::kOverwrite);
        // eta/pi0 Vogelsang PDF: CT10, FF eta- AESSS, FF pi0 - DSS07
        graphEtaToPi0NLOMuHalf900GeV->Write("graphNLOCalcEtaOverPi0MuHalf900GeV", TObject::kOverwrite);
        graphEtaToPi0NLOMuOne900GeV->Write("graphNLOCalcEtaOverPi0MuOne900GeV", TObject::kOverwrite);
        graphEtaToPi0NLOMuTwo900GeV->Write("graphNLOCalcEtaOverPi0MuTwo900GeV", TObject::kOverwrite);
//         // pi0, eta, eta/pi0 Pythia 8.1 Monash 2013
//         histoPi0Pythia8MonashInvSec900GeV->Write("histoInvSecPythia8Monash2013Pi0900GeV", TObject::kOverwrite);
//         histoEtaPythia8MonashInvSec900GeV->Write("histoInvSecPythia8Monash2013Eta900GeV", TObject::kOverwrite);
//         histoEtaToPi0RatioPythia8Monash900GeV->Write("histoEtaToPi0RatioPythia8Monash900GeV", TObject::kOverwrite);
        // pi0, eta, eta/pi0 Pythia 8.1 Monash 2013 (central alice)
        histoPi0Pythia8MonashInvSec900GeVLego->Write("histoInvSecPythia8Monash2013LegoPi0900GeV", TObject::kOverwrite);
        histoEtaPythia8MonashInvSec900GeVLego->Write("histoInvSecPythia8Monash2013LegoEta900GeV", TObject::kOverwrite);
        histoEtaToPi0RatioPythia8Monash900GeVLego->Write("histoEtaToPi0RatioPythia8Monash2013Lego900GeV", TObject::kOverwrite);

        
        //***********************************************************************
        // write  calculations for 2.76TeV
        //***********************************************************************
        // pi0 Vogelsang PDF: CT10, FF DSS07
        graphNLOCalcInvSecPi0MuHalf2760GeV->Write("graphNLOCalcInvSecPi0MuHalf2760GeV", TObject::kOverwrite);
        graphNLOCalcInvSecPi0MuOne2760GeV->Write("graphNLOCalcInvSecPi0MuOne2760GeV", TObject::kOverwrite);
        graphNLOCalcInvSecPi0MuTwo2760GeV->Write("graphNLOCalcInvSecPi0MuTwo2760GeV", TObject::kOverwrite);
        graphNLOCalcInvSecPi02760GeV->Write("graphNLOCalcDSS07InvSecPi02760GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldPi0MuHalf2760GeV->Write("graphNLOCalcInvYieldPi0MuHalf2760GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldPi0MuOne2760GeV->Write("graphNLOCalcInvYieldPi0MuOne2760GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldPi0MuTwo2760GeV->Write("graphNLOCalcInvYieldPi0MuTwo2760GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldPi02760GeV->Write("graphNLOCalcDSS07InvYieldPi02760GeV", TObject::kOverwrite);
        // pi0 INCNLO FF:DSS07
        graphNLOCalcDSSInvSecPi0MuTwo2760GeV->Write("graphNLOCalcDSSInvSecPi0MuTwo2760GeV", TObject::kOverwrite);
        graphNLOCalcDSSInvYieldPi0MuTwo2760GeV->Write("graphNLOCalcDSSInvYieldPi0MuTwo2760GeV", TObject::kOverwrite);
        // eta Vogelsang PDF: CT10, FF AESSS
        graphNLOCalcInvSecEtaMuHalf2760GeV->Write("graphNLOCalcInvSecEtaMuHalf2760GeV", TObject::kOverwrite);
        graphNLOCalcInvSecEtaMuOne2760GeV->Write("graphNLOCalcInvSecEtaMuOne2760GeV", TObject::kOverwrite);
        graphNLOCalcInvSecEtaMuTwo2760GeV->Write("graphNLOCalcInvSecEtaMuTwo2760GeV", TObject::kOverwrite);
        graphNLOCalcInvSecEta2760GeV->Write("graphNLOCalcAESSSInvSecEta2760GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldEtaMuHalf2760GeV->Write("graphNLOCalcInvYieldEtaMuHalf2760GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldEtaMuOne2760GeV->Write("graphNLOCalcInvYieldEtaMuOne2760GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldEtaMuTwo2760GeV->Write("graphNLOCalcInvYieldEtaMuTwo2760GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldEta2760GeV->Write("graphNLOCalcAESSSInvYieldEta2760GeV", TObject::kOverwrite);
        // eta/pi0 Vogelsang PDF: CT10, FF eta- AESSS, FF pi0 - DSS07
        graphEtaToPi0NLOMuHalf2760GeV->Write("graphNLOCalcEtaOverPi0MuHalf2760GeV", TObject::kOverwrite);
        graphEtaToPi0NLOMuOne2760GeV->Write("graphNLOCalcEtaOverPi0MuOne2760GeV", TObject::kOverwrite);
        graphEtaToPi0NLOMuTwo2760GeV->Write("graphNLOCalcEtaOverPi0MuTwo2760GeV", TObject::kOverwrite);
        graphNLOCalcEtaToPi02760GeV->Write("graphNLOCalcEtaOverPi02760GeV_AESSS_DSS07", TObject::kOverwrite);
        graphNLOCalcEtaToPi02760GeV->Print();
        // pi0 Stratmann PDF: CT10, FF: DSS14
        graphNLOCalcDSS14InvSecPi02760GeV->Write("graphNLOCalcDSS14InvCrossSec2760GeV", TObject::kOverwrite);
        graphNLOCalcDSS14InvYieldPi02760GeV->Write("graphNLOCalcDSS14InvYield2760GeV", TObject::kOverwrite);
        // pi0 CGC calc
        graphInvYieldCGC2760GeV->Write("graphNLOCalcCGCInvYield2760GeV", TObject::kOverwrite);
        graphInvCrossSecCGC2760GeV->Write("graphNLOCalcCGCInvCrossSec2760GeV", TObject::kOverwrite);
//         // pi0, eta, eta/pi0 Pythia 8.1 Monash 2013 (Satoshi)
//         histoPi0Pythia8MonashInvSec2760GeV->Write("histoInvSecPythia8Monash2013Pi02760GeV", TObject::kOverwrite);
//         histoEtaPythia8MonashInvSec2760GeV->Write("histoInvSecPythia8Monash2013Eta2760GeV", TObject::kOverwrite);
//         histoEtaToPi0RatioPythia8Monash2760GeV->Write("histoEtaToPi0RatioPythia8Monash20132760GeV", TObject::kOverwrite);
        // pi0, eta, eta/pi0 Pythia 8.1 Monash 2013 (central alice)
        histoPi0Pythia8MonashInvSec2760GeVLego->Write("histoInvSecPythia8Monash2013LegoPi02760GeV", TObject::kOverwrite);
        histoEtaPythia8MonashInvSec2760GeVLego->Write("histoInvSecPythia8Monash2013LegoEta2760GeV", TObject::kOverwrite);
        histoEtaToPi0RatioPythia8Monash2760GeVLego->Write("histoEtaToPi0RatioPythia8Monash2013Lego2760GeV", TObject::kOverwrite);
        // pi0, eta, eta/pi0 Pythia 8.0 Tune4C (Klaus Reygers)
        histoPythia8Spec2760GeV->Write("histoInvSecPythia8Spec2760GeV", TObject::kOverwrite);
        histoPythia8InvYield2760GeV->Write("histoInvYieldPythia8Spec2760GeV", TObject::kOverwrite);
        histoPythia8Spec2760GeVVarBinning->Write("histoInvSecPythia8Spec2760GeVVarBinning", TObject::kOverwrite);
        histoPythia8InvYield2760GeVVarBinning->Write("histoInvYieldPythia8Spec2760GeVVarBinning", TObject::kOverwrite);
        

        //***********************************************************************
        // write  calculations for 5.023TeV
        //***********************************************************************
        // pi0 Vogelsang PDF: CT10, FF: DSS14
        graphNLOCalcInvSecPi0MuHalf5023GeV->Write("graphNLOCalcDSS14InvSecPi0MuHalf5023GeV", TObject::kOverwrite);
        graphNLOCalcInvSecPi0MuOne5023GeV->Write("graphNLOCalcDSS14InvSecPi0MuOne5023GeV", TObject::kOverwrite);
        graphNLOCalcInvSecPi0MuTwo5023GeV->Write("graphNLOCalcDSS14InvSecPi0MuTwo5023GeV", TObject::kOverwrite);
        graphNLOCalcInvSecPi05023GeV->Write("graphNLOCalcDSS14InvSecPi05023GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldPi0MuHalf5023GeV->Write("graphNLOCalcDSS14InvYieldPi0MuHalf5023GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldPi0MuOne5023GeV->Write("graphNLOCalcDSS14InvYieldPi0MuOne5023GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldPi0MuTwo5023GeV->Write("graphNLOCalcDSS14InvYieldPi0MuTwo5023GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldPi05023GeV->Write("graphNLOCalcDSS14InvYieldPi05023GeV", TObject::kOverwrite);
        // eta Vogelsang PDF: CT10, FF AESSS
        if (graphNLOCalcInvSecEtaMuHalf5023GeV) graphNLOCalcInvSecEtaMuHalf5023GeV->Write("graphNLOCalcAESSSInvSecEtaMuHalf5023GeV", TObject::kOverwrite);
        if (graphNLOCalcInvSecEtaMuOne5023GeV) graphNLOCalcInvSecEtaMuOne5023GeV->Write("graphNLOCalcAESSSInvSecEtaMuOne5023GeV", TObject::kOverwrite);
        if (graphNLOCalcInvSecEtaMuTwo5023GeV) graphNLOCalcInvSecEtaMuTwo5023GeV->Write("graphNLOCalcAESSSInvSecEtaMuTwo5023GeV", TObject::kOverwrite);
        graphNLOCalcInvSecEta5023GeV->Write("graphNLOCalcAESSSInvSecEta5023GeV", TObject::kOverwrite);
        if (graphNLOCalcInvYieldEtaMuHalf5023GeV) graphNLOCalcInvYieldEtaMuHalf5023GeV->Write("graphNLOCalcAESSSInvYieldEtaMuHalf5023GeV", TObject::kOverwrite);
        if (graphNLOCalcInvYieldEtaMuOne5023GeV) graphNLOCalcInvYieldEtaMuOne5023GeV->Write("graphNLOCalcAESSSInvYieldEtaMuOne5023GeV", TObject::kOverwrite);
        if (graphNLOCalcInvYieldEtaMuTwo5023GeV) graphNLOCalcInvYieldEtaMuTwo5023GeV->Write("graphNLOCalcAESSSInvYieldEtaMuTwo5023GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldEta5023GeV->Write("graphNLOCalcAESSSInvYieldEta5023GeV", TObject::kOverwrite);
        // eta/pi0 Vogelsang PDF: CT10, FF eta- AESSS, FF pi0 - DSS14
        if (graphEtaToPi0NLOMuHalf5023GeV) graphEtaToPi0NLOMuHalf5023GeV->Write("graphNLOCalcEtaOverPi0MuHalf5023GeV", TObject::kOverwrite);
        if (graphEtaToPi0NLOMuOne5023GeV) graphEtaToPi0NLOMuOne5023GeV->Write("graphNLOCalcEtaOverPi0MuOne5023GeV", TObject::kOverwrite);
        if (graphEtaToPi0NLOMuTwo5023GeV) graphEtaToPi0NLOMuTwo5023GeV->Write("graphNLOCalcEtaOverPi0MuTwo5023GeV", TObject::kOverwrite);
        graphNLOCalcEtaToPi05023GeV->Write("graphNLOCalcEtaOverPi05023GeV_AESSS_DSS14", TObject::kOverwrite);
//         // pi0, eta, eta/pi0 Pythia 8.1 Monash 2013 (Satoshi)
//         histoPi0Pythia8MonashInvSec5023GeV->Write("histoInvSecPythia8Monash2013Pi05023GeV", TObject::kOverwrite);
//         histoEtaPythia8MonashInvSec5023GeV->Write("histoInvSecPythia8Monash2013Eta5023GeV", TObject::kOverwrite);
//         histoEtaToPi0RatioPythia8Monash5023GeV->Write("histoEtaToPi0RatioPythia8Monash5023GeV", TObject::kOverwrite);
        // pi0, eta, eta/pi0 Pythia 8.1 Monash 2013 (central alice)
        histoPi0Pythia8MonashInvSec5TeVLego->Write("histoInvSecPythia8Monash2013LegoPi05TeV", TObject::kOverwrite);
        histoEtaPythia8MonashInvSec5TeVLego->Write("histoInvSecPythia8Monash2013LegoEta5TeV", TObject::kOverwrite);
        histoEtaToPi0RatioPythia8Monash5TeVLego->Write("histoEtaToPi0RatioPythia8Monash2013Lego5TeV", TObject::kOverwrite);

        //***********************************************************************
        // write  calculations for 7TeV
        //***********************************************************************
        // pi0 Vogelsang PDF: CT10, FF DSS07
        graphNLOCalcInvSecPi0MuHalf7000GeV->Write("graphNLOCalcInvSecPi0MuHalf7000GeV", TObject::kOverwrite);
        graphNLOCalcInvSecPi0MuOne7000GeV->Write("graphNLOCalcInvSecPi0MuOne7000GeV", TObject::kOverwrite);
        graphNLOCalcInvSecPi0MuTwo7000GeV->Write("graphNLOCalcInvSecPi0MuTwo7000GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldPi0MuHalf7000GeV->Write("graphNLOCalcInvYieldPi0MuHalf7000GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldPi0MuOne7000GeV->Write("graphNLOCalcInvYieldPi0MuOne7000GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldPi0MuTwo7000GeV->Write("graphNLOCalcInvYieldPi0MuTwo7000GeV", TObject::kOverwrite); 
        // pi0 INCNLO FF:BKK, DSS07
        graphNLOCalcBKKInvSecPi0MuTwo7000GeV->Write("graphNLOCalcBKKInvSecPi0MuTwo7000GeV", TObject::kOverwrite);
        graphNLOCalcDSSInvSecPi0MuTwo7000GeV->Write("graphNLOCalcDSSInvSecPi0MuTwo7000GeV", TObject::kOverwrite);
        graphNLOCalcBKKInvYieldPi0MuTwo7000GeV->Write("graphNLOCalcBKKInvYieldPi0MuTwo7000GeV", TObject::kOverwrite);
        graphNLOCalcDSSInvYieldPi0MuTwo7000GeV->Write("graphNLOCalcDSSInvYieldPi0MuTwo7000GeV", TObject::kOverwrite);
        // eta Vogelsang PDF: CT10, FF AESSS
        graphNLOCalcInvSecEtaMuHalf7000GeV->Write("graphNLOCalcInvSecEtaMuHalf7000GeV", TObject::kOverwrite);
        graphNLOCalcInvSecEtaMuOne7000GeV->Write("graphNLOCalcInvSecEtaMuOne7000GeV", TObject::kOverwrite);
        graphNLOCalcInvSecEtaMuTwo7000GeV->Write("graphNLOCalcInvSecEtaMuTwo7000GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldEtaMuHalf7000GeV->Write("graphNLOCalcInvYieldEtaMuHalf7000GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldEtaMuOne7000GeV->Write("graphNLOCalcInvYieldEtaMuOne7000GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldEtaMuTwo7000GeV->Write("graphNLOCalcInvYieldEtaMuTwo7000GeV", TObject::kOverwrite);
        // eta/pi0 Vogelsang PDF: CT10, FF eta- AESSS, FF pi0 - DSS07
        graphEtaToPi0NLOMuHalf7TeV->Write("graphNLOCalcEtaOverPi0MuHalf7000GeV", TObject::kOverwrite);
        graphEtaToPi0NLOMuOne7TeV->Write("graphNLOCalcEtaOverPi0MuOne7000GeV", TObject::kOverwrite);
        graphEtaToPi0NLOMuTwo7TeV->Write("graphNLOCalcEtaOverPi0MuTwo7000GeV", TObject::kOverwrite);
        // pi0 Stratmann FF:DSS14
        graphNLOCalcDSS14InvSecPi07000GeV->Write("graphNLOCalcDSS14InvCrossSec7000GeV", TObject::kOverwrite);
        graphNLOCalcDSS14InvYieldPi07000GeV->Write("graphNLOCalcDSS14InvYield7000GeV", TObject::kOverwrite);
        // pi0 CGC calc
        graphInvYieldCGC7000GeV->Write("graphNLOCalcCGCInvYield7000GeV", TObject::kOverwrite);
        graphInvCrossSecCGC7000GeV->Write("graphNLOCalcCGCInvCrossSec7000GeV", TObject::kOverwrite);
        graphInvCrossSecCGC7000GeV_mvgamma->Write("graphNLOCalcCGCInvCrossSec7000GeV_mvgamma", TObject::kOverwrite);
//         // pi0, eta, eta/pi0 Pythia 8.1 Monash 2013 (Satoshi)
//         histoPi0Pythia8MonashInvSec7TeV->Write("histoInvSecPythia8Monash2013Pi07TeV", TObject::kOverwrite);
//         histoEtaPythia8MonashInvSec7TeV->Write("histoInvSecPythia8Monash2013Eta7TeV", TObject::kOverwrite);
//         histoEtaToPi0RatioPythia8Monash7TeV->Write("histoEtaToPi0RatioPythia8Monash7TeV", TObject::kOverwrite);
        // pi0, eta, eta/pi0 Pythia 8.1 Monash 2013 (central alice)
        histoPi0Pythia8MonashInvSec7TeVLego->Write("histoInvSecPythia8Monash2013LegoPi07TeV", TObject::kOverwrite);
        histoEtaPythia8MonashInvSec7TeVLego->Write("histoInvSecPythia8Monash2013LegoEta7TeV", TObject::kOverwrite);
        histoEtaToPi0RatioPythia8Monash7TeVLego->Write("histoEtaToPi0RatioPythia8Monash2013Lego7TeV", TObject::kOverwrite);
        
        // pi0, eta, eta/pi0 Pythia 6.4 Perugia2011
        histoPi07TeVPythia6->Write("histoPi0Pythia6_Perugia2011_7000TeV", TObject::kOverwrite);
        histoEta7TeVPythia6->Write("histoEtaPythia6_Perugia2011_7000TeV", TObject::kOverwrite);
        histoPi07TeVPythia6Reb->Write("histoPi0Pythia6_Perugia2011_7000TeV_Reb", TObject::kOverwrite);
        histoEta7TeVPythia6Reb->Write("histoEtaPythia6_Perugia2011_7000TeV_Reb", TObject::kOverwrite);
        histoEtaToPi07TeVPythia6->Write("histoEtaToPi0Pythia6_Perugia2011_7000TeV", TObject::kOverwrite);

        //***********************************************************************
        // write  calculations for 8TeV
        //***********************************************************************
        // pi0 Vogelsang PDF: CT10, FF DSS07
        graphNLOCalcInvSecPi0MuHalf8000GeV->Write("graphNLOCalcInvSecPi0MuHalf8000GeV", TObject::kOverwrite);
        graphNLOCalcInvSecPi0MuOne8000GeV->Write("graphNLOCalcInvSecPi0MuOne8000GeV", TObject::kOverwrite);
        graphNLOCalcInvSecPi0MuTwo8000GeV->Write("graphNLOCalcInvSecPi0MuTwo8000GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldPi0MuHalf8000GeV->Write("graphNLOCalcInvYieldPi0MuHalf8000GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldPi0MuOne8000GeV->Write("graphNLOCalcInvYieldPi0MuOne8000GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldPi0MuTwo8000GeV->Write("graphNLOCalcInvYieldPi0MuTwo8000GeV", TObject::kOverwrite);
        // pi0 INCNLO FF:BKK, DSS07
        // eta Vogelsang PDF: CT10, FF AESSS
        graphNLOCalcInvSecEtaMuHalf8000GeV->Write("graphNLOCalcInvSecEtaMuHalf8000GeV", TObject::kOverwrite);
        graphNLOCalcInvSecEtaMuOne8000GeV->Write("graphNLOCalcInvSecEtaMuOne8000GeV", TObject::kOverwrite);
        graphNLOCalcInvSecEtaMuTwo8000GeV->Write("graphNLOCalcInvSecEtaMuTwo8000GeV", TObject::kOverwrite);
        graphNLOCalcInvSecEta8000GeV->Write("graphNLOCalcAESSSInvSecEta8000GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldEtaMuHalf8000GeV->Write("graphNLOCalcInvYieldEtaMuHalf8000GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldEtaMuOne8000GeV->Write("graphNLOCalcInvYieldEtaMuOne8000GeV", TObject::kOverwrite);
        graphNLOCalcInvYieldEtaMuTwo8000GeV->Write("graphNLOCalcInvYieldEtaMuTwo8000GeV", TObject::kOverwrite);
        // eta/pi0 Vogelsang PDF: CT10, FF eta- AESSS, FF pi0 - DSS07
        graphEtaToPi0NLOMuHalf8TeV->Write("graphNLOCalcEtaOverPi0MuHalf8000GeV", TObject::kOverwrite);
        graphEtaToPi0NLOMuOne8TeV->Write("graphNLOCalcEtaOverPi0MuOne8000GeV", TObject::kOverwrite);
        graphEtaToPi0NLOMuTwo8TeV->Write("graphNLOCalcEtaOverPi0MuTwo8000GeV", TObject::kOverwrite);
        graphNLOCalcEtaToPi08000GeV->Write("graphNLOCalcEtaOverPi08000GeV_AESSS_DSS07",TObject::kOverwrite);
        // pi0 Stratmann PDF: CT10, FF: DSS14, 
        graphNLOCalcInvSecPi08000GeV->Write("graphNLOCalcDSS14InvSecPi08000GeV", TObject::kOverwrite);
        // pi0, eta, eta/pi0 Pythia 8.1 Monash 2013 (Satoshi)
        histoPi0Pythia8MonashInvSec8TeV->Write("histoInvSecPythia8Monash2013Pi08TeV", TObject::kOverwrite);
        histoEtaPythia8MonashInvSec8TeV->Write("histoInvSecPythia8Monash2013Eta8TeV", TObject::kOverwrite);
        histoEtaToPi0RatioPythia8Monash8TeV->Write("histoEtaToPi0RatioPythia8Monash8TeV", TObject::kOverwrite);
        // pi0, eta, eta/pi0 Pythia 8.1 Monash 2013 (central)
        histoPi0Pythia8MonashInvSec8TeVLego->Write("histoInvSecPythia8Monash2013LegoPi08TeV", TObject::kOverwrite);
        histoEtaPythia8MonashInvSec8TeVLego->Write("histoInvSecPythia8Monash2013LegoEta8TeV", TObject::kOverwrite);
        histoEtaToPi0RatioPythia8Monash8TeVLego->Write("histoEtaToPi0RatioPythia8Monash2013Lego8TeV", TObject::kOverwrite);
        // pi0, eta, eta/pi0 Pythia 8.1 Tune4C (Satoshi)
        histoPi08TeVPythia8->Write("histoPi0Pythia8_8000TeV", TObject::kOverwrite);
        histoEta8TeVPythia8->Write("histoEtaPythia8_8000TeV", TObject::kOverwrite);
        histoPi08TeVPythia8Reb->Write("histoPi0Pythia8_8000TeV_Reb", TObject::kOverwrite);
        histoEta8TeVPythia8Reb->Write("histoEtaPythia8_8000TeV_Reb", TObject::kOverwrite);
        histoEtaToPi08TeVPythia8->Write("histoEtaToPi0Pythia8_8000TeV", TObject::kOverwrite);
        // pi0, eta, eta/pi0 Pythia 8.1 Tune4C  (central)
        histoPi0Pythia8Tune4CInvSec8TeVLego->Write("histoInvSecPythia8Tune4CLegoPi08TeV", TObject::kOverwrite);
        histoEtaPythia8Tune4CInvSec8TeVLego->Write("histoInvSecPythia8Tune4CLegoEta8TeV", TObject::kOverwrite);
        histoEtaToPi0RatioPythia8Tune4C8TeVLego->Write("histoEtaToPi0RatioPythia8Tune4CLego8TeV", TObject::kOverwrite);
        // pi0, eta, eta/pi0 Phojet
        histoPi08TeVPhojet->Write("histoPi0Phojet_8000TeV", TObject::kOverwrite);
        histoEta8TeVPhojet->Write("histoEtaPhojet_8000TeV", TObject::kOverwrite);
        histoPi08TeVPhojetReb->Write("histoPi0Phojet_8000TeV_Reb", TObject::kOverwrite);
        histoEta8TeVPhojetReb->Write("histoEtaPhojet_8000TeV_Reb", TObject::kOverwrite);
        histoEtaToPi08TeVPhojet->Write("histoEtaToPi0Phojet_8000TeV", TObject::kOverwrite);
        
        
    fileTheoryGraphsPP.Close();

}
