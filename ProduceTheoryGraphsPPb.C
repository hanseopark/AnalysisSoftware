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
#include "CommonHeaders/ConversionFunctions.h"
#include "TSpline.h"

extern TRandom*    gRandom;
extern TBenchmark*    gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*      gMinuit;

Double_t xSection5023GeVINEL        = 67.6*1e-3;    // from https://aliceinfo.cern.ch/Notes/sites/aliceinfo.cern.ch.Notes/files/notes/analysis/stripath/2017-Jun-14-analysis_note-INEL_norm.pdf
Double_t xSection5023GeVINELErr     = 0.6;          // from https://aliceinfo.cern.ch/Notes/sites/aliceinfo.cern.ch.Notes/files/notes/analysis/stripath/2017-Jun-14-analysis_note-INEL_norm.pdf

Double_t xSection5023GeVINELpPb     = 70*1e-3;
Double_t ncollpPb5023GeV            = 6.9;
Double_t recalcBarn                 = 1e12; //NLO in pbarn!!!!


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

//**********************************************************************************************************
// Combine graphs from different muo scales into 1 accordingly
//**********************************************************************************************************
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

//**********************************************************************************************************
// Combine graphs from different muo scales into 1 accordingly
//**********************************************************************************************************
TGraphAsymmErrors* CombineMuScales( TGraph* mu,
                                    TGraph* muHalf,
                                    TGraph* muTwo
){
    cout << "combining different mu scales" << endl;
    Int_t nPoints               = mu->GetN();
    Double_t* pt                = mu->GetX();
    Double_t yValue[400];
    Double_t yErrorHigh[400];
    Double_t yErrorLow[400];
    Double_t ptErr[400];

    for (Int_t i = 0; i< nPoints; i++){
        if (pt[i] == mu->GetX()[i] && pt[i] == muHalf->GetX()[i] && pt[i] == muTwo->GetX()[i] ){
            yValue[i]   = mu->GetY()[i];
            ptErr[i]    = 0;
            if (muHalf->GetY()[i] < mu->GetY()[i] && muHalf->GetY()[i] != 0){
                yErrorLow[i]    = TMath::Abs(mu->GetY()[i]-muHalf->GetY()[i]);
                yErrorHigh[i]   = TMath::Abs(muTwo->GetY()[i]-mu->GetY()[i]);
            } else {
                yErrorLow[i]    = TMath::Abs(mu->GetY()[i]-muTwo->GetY()[i]);
                yErrorHigh[i]   = TMath::Abs(muHalf->GetY()[i]-mu->GetY()[i]);
            }
        } else {
            cout << "*************************************************************" << endl;
            cout << "********** ATTENTION mismatching x values *******************" << endl;
            cout << mu->GetX()[i] << "\t" << muHalf->GetX()[i]  << "\t" << muTwo->GetX()[i] << endl;
            cout << "*************************************************************" << endl;
            yValue[i]           = 0;
            yErrorLow[i]        = 0;
            yErrorHigh[i]       = 0;
        }
    }
    TGraphAsymmErrors* graphReturn = new TGraphAsymmErrors(nPoints, pt, yValue, ptErr, ptErr, yErrorLow, yErrorHigh);
    return graphReturn;
}


// //**********************************************************************************************************
// // Calculates the ratio of two error graphs
// //**********************************************************************************************************
// TGraph* CalculateGraphRatioToGraph(TGraph* graphA, TGraph* graphB){

//     TGraphErrors* graphACopy        = (TGraphErrors*)graphA->Clone("GraphCopy");
//     Double_t* xValue                = graphACopy->GetX();
//     Double_t* yValue                = graphACopy->GetY();
//     Int_t nPoints                   = graphACopy->GetN();

//     for (Int_t i = 0; i < nPoints; i++){
//         yValue[i]                   = yValue[i]/graphB->GetY()[i];
//     }
//     TGraph* returnGraph       = new TGraph(nPoints,xValue,yValue);
//     return returnGraph;
// }


//**********************************************************************************************************
// Calculates the ratio of two error graphs
//**********************************************************************************************************
TGraphErrors* CalculateGraphRatioToGraphErrs(TGraphErrors* graphA, TGraphErrors* graphB){

    TGraphErrors* graphACopy        = (TGraphErrors*)graphA->Clone("GraphCopy");
    Double_t* xValue                = graphACopy->GetX();
    Double_t* yValue                = graphACopy->GetY();
    Double_t* xError                = graphACopy->GetEX();
    Double_t* yError                = graphACopy->GetEY();
    Int_t nPoints                   = graphACopy->GetN();

    for (Int_t i = 0; i < nPoints; i++){
        yValue[i]                   = yValue[i]/graphB->GetY()[i];
        cout<<i<<" A EY "<< graphA->GetEY()[i]<<" A Y "<<graphA->GetY()[i]<<" B EY "<<graphB->GetEY()[i]<<" B Y"<<graphB->GetY()[i]<<endl;
        Double_t yErrorRatio        = yValue[i]*TMath::Sqrt( TMath::Power(graphA->GetEY()[i]/graphA->GetY()[i],2) + TMath::Power(graphB->GetEY()[i]/graphB->GetY()[i],2));
        yError[i] = TMath::Abs(yErrorRatio);
    }
    TGraphErrors* returnGraph       = new TGraphErrors(nPoints,xValue,yValue,xError,yError);
    return returnGraph;
}

//**********************************************************************************************************
// Calculates the ratio of two error graphs
//**********************************************************************************************************
TGraphAsymmErrors* CalculateGraphRatioToGraphAsym(TGraphAsymmErrors* graphA, TGraphAsymmErrors* graphB){

    TGraphErrors* graphACopy        = (TGraphErrors*)graphA->Clone("GraphCopy");
    Double_t* xValue                = graphACopy->GetX();
    Double_t* yValue                = graphACopy->GetY();
    Double_t* xErrorHigh            = graphACopy->GetEXhigh();
    Double_t* yErrorHigh            = graphACopy->GetEYhigh();
    Double_t* xErrorLow             = graphACopy->GetEXlow();
    Double_t* yErrorLow             = graphACopy->GetEYlow();
    Int_t nPoints                   = graphACopy->GetN();

    for (Int_t i = 0; i < nPoints; i++){
        yValue[i]                   = yValue[i]/graphB->GetY()[i];
        Double_t yErrorRatioLow     = yValue[i]*TMath::Sqrt( TMath::Power(graphA->GetEYlow()[i]/graphA->GetY()[i],2) + TMath::Power(graphB->GetEYlow()[i]/graphB->GetY()[i],2));
        Double_t yErrorRatioHigh    = yValue[i]*TMath::Sqrt( TMath::Power(graphA->GetEYhigh()[i]/graphA->GetY()[i],2) + TMath::Power(graphB->GetEYhigh()[i]/graphB->GetY()[i],2));
        yErrorLow[i]                = TMath::Abs(yErrorRatioLow);
        yErrorHigh[i]               = TMath::Abs(yErrorRatioHigh);
    }
    TGraphAsymmErrors* returnGraph  = new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh, yErrorLow, yErrorHigh);
    return returnGraph;
}

//**********************************************************************************************************
// Calculates the ratio of two error graphs
//**********************************************************************************************************
TGraphAsymmErrors* CalculateGraphRatioToGraphAsym2(TGraphAsymmErrors* graphA, TGraph* graphB){

    TGraphErrors* graphACopy        = (TGraphErrors*)graphA->Clone("GraphCopy");
    Double_t* xValue                = graphACopy->GetX();
    Double_t* yValue                = graphACopy->GetY();
    Double_t* xErrorHigh            = graphACopy->GetEXhigh();
    Double_t* yErrorHigh            = graphACopy->GetEYhigh();
    Double_t* xErrorLow             = graphACopy->GetEXlow();
    Double_t* yErrorLow             = graphACopy->GetEYlow();
    Int_t nPoints                   = graphACopy->GetN();

    for (Int_t i = 0; i < nPoints; i++){
        yValue[i]                   = yValue[i]/graphB->GetY()[i];
        Double_t yErrorRatioLow     = yValue[i]*TMath::Sqrt( TMath::Power(graphA->GetEYlow()[i]/graphA->GetY()[i],2));
        Double_t yErrorRatioHigh    = yValue[i]*TMath::Sqrt( TMath::Power(graphA->GetEYhigh()[i]/graphA->GetY()[i],2));
        yErrorLow[i]                = TMath::Abs(yErrorRatioLow);
        yErrorHigh[i]               = TMath::Abs(yErrorRatioHigh);
    }
    TGraphAsymmErrors* returnGraph  = new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh, yErrorLow, yErrorHigh);
    return returnGraph;
}

//**************************************************************************************************
//**************************** Main function *******************************************************
//**************************************************************************************************
void ProduceTheoryGraphsPPb(){

    StyleSettingsThesis();
    SetPlotStyle();

    const char *fileNameEPS09sPi0AKK    = "ExternalInputpPb/Theory/R_pi0_pPb_y0_eps09s/R_pPb5000_y0_eps09s_akk_mb.dat";
    const char *fileNameEPS09sPi0DSS    = "ExternalInputpPb/Theory/R_pi0_pPb_y0_eps09s/R_pPb5000_y0_eps09s_fdss_mb.dat";
    const char *fileNameEPS09sPi0KKP    = "ExternalInputpPb/Theory/R_pi0_pPb_y0_eps09s/R_pPb5000_y0_eps09s_kkp_mb.dat";
    const char *fileNameEPS09sPi0CGC    = "ExternalInputpPb/Theory/ColorGlassCondensate.dat";

    //**************************************************************************************************
    //****************************** extracting EPS09s predictions**************************************
    //**************************************************************************************************
    //     @article{Helenius:2012wd,
    //       author         = "Helenius, Ilkka and Eskola, Kari J. and Honkanen, Heli
    //                         and Salgado, Carlos A.",
    //       title          = "{Impact-Parameter Dependent Nuclear Parton Distribution
    //                         Functions: EPS09s and EKS98s and Their Applications in
    //                         Nuclear Hard Processes}",
    //       journal        = "JHEP",
    //       volume         = "07",
    //       year           = "2012",
    //       pages          = "073",
    //       doi            = "10.1007/JHEP07(2012)073",
    //       eprint         = "1205.5359",
    //       archivePrefix  = "arXiv",
    //       primaryClass   = "hep-ph",
    //       reportNumber   = "CERN-PH-TH-2012-145",
    //       SLACcitation   = "%%CITATION = ARXIV:1205.5359;%%"
    //     }

    // FF DSS07
    ifstream inDSS;
    Int_t nlinesEPSsPi0fDSS = 0;
    inDSS.open(fileNameEPS09sPi0DSS,ios_base::in);
    Double_t xEPSsPi0fDSS[100],yEPSsPi0fDSS[100];
    Double_t xUpErrorEPSsPi0DSS[100],xDownErrorEPSsPi0DSS[100];
    Double_t yUpErrorEPSsPi0DSS[100],yDownErrorEPSsPi0DSS[100];

    while(!inDSS.eof()){

        inDSS >> xEPSsPi0fDSS[nlinesEPSsPi0fDSS]  >> yEPSsPi0fDSS[nlinesEPSsPi0fDSS] >> yUpErrorEPSsPi0DSS[nlinesEPSsPi0fDSS]>>yDownErrorEPSsPi0DSS[nlinesEPSsPi0fDSS];
        yUpErrorEPSsPi0DSS[nlinesEPSsPi0fDSS] = ( yUpErrorEPSsPi0DSS[nlinesEPSsPi0fDSS] - yEPSsPi0fDSS[nlinesEPSsPi0fDSS] ) /yEPSsPi0fDSS[nlinesEPSsPi0fDSS];
        yDownErrorEPSsPi0DSS[nlinesEPSsPi0fDSS] = -1*( yDownErrorEPSsPi0DSS[nlinesEPSsPi0fDSS] - yEPSsPi0fDSS[nlinesEPSsPi0fDSS] ) /yEPSsPi0fDSS[nlinesEPSsPi0fDSS];
        xUpErrorEPSsPi0DSS[nlinesEPSsPi0fDSS]   = 0;
        xDownErrorEPSsPi0DSS[nlinesEPSsPi0fDSS] = 0;

        cout << nlinesEPSsPi0fDSS << "         "  << xEPSsPi0fDSS[nlinesEPSsPi0fDSS] << "         "  <<yEPSsPi0fDSS[nlinesEPSsPi0fDSS]<<"         "<<yUpErrorEPSsPi0DSS[nlinesEPSsPi0fDSS]<<"          "<<yDownErrorEPSsPi0DSS[nlinesEPSsPi0fDSS]<< endl;
        nlinesEPSsPi0fDSS++;
    }
    inDSS.close();

    TGraph* graphPi0RpAEPS09sDSS                        = new TGraph(nlinesEPSsPi0fDSS-1,xEPSsPi0fDSS,yEPSsPi0fDSS);
    TGraphAsymmErrors* graphPi0RpAAsymmErrEPS09sDSS     = new TGraphAsymmErrors(nlinesEPSsPi0fDSS-1,xEPSsPi0fDSS,yEPSsPi0fDSS,xDownErrorEPSsPi0DSS,xUpErrorEPSsPi0DSS,yDownErrorEPSsPi0DSS,yUpErrorEPSsPi0DSS);
//     graphPi0RpAEPS09sDSS->RemovePoint(0);
//     graphPi0RpAAsymmErrEPS09sDSS->RemovePoint(0);


    // FF AKK
    ifstream inAKK;
    Int_t nlinesPi0AKK = 0;
    inAKK.open(fileNameEPS09sPi0AKK,ios_base::in);
    Double_t xESPsPi0AKK[100], yESPsPi0AKK[100];

    while(!inAKK.eof()){
        inAKK >> xESPsPi0AKK[nlinesPi0AKK]  >> yESPsPi0AKK[nlinesPi0AKK];
        cout << nlinesPi0AKK << "         "  << xESPsPi0AKK[nlinesPi0AKK] << "         "  <<xESPsPi0AKK[nlinesPi0AKK]<<endl;
        nlinesPi0AKK++;
    }
    inAKK.close();

    TGraph* graphPi0RpAEPS09sAKK                        = new TGraph(nlinesPi0AKK,xESPsPi0AKK,yESPsPi0AKK);
    graphPi0RpAEPS09sAKK->RemovePoint(0);

    // FF KKP
    ifstream inKKP;
    Int_t nlinesPi0KKP = 0;
    inKKP.open(fileNameEPS09sPi0KKP,ios_base::in);
    Double_t xESPsPi0KKP[100], yESPsPi0KKP[100];

    while(!inKKP.eof()){
        inKKP >> xESPsPi0KKP[nlinesPi0KKP]  >> yESPsPi0KKP[nlinesPi0KKP];
        cout << nlinesPi0KKP << "         "  << xESPsPi0KKP[nlinesPi0KKP] << "         "  <<xESPsPi0KKP[nlinesPi0KKP]<<endl;
        nlinesPi0KKP++;
    }
    inKKP.close();
    TGraph* graphPi0RpAEPS09sKKP = new TGraph(nlinesPi0KKP,xESPsPi0KKP,yESPsPi0KKP);

    //**************************************************************************************************
    //****************************** extracting CGC predictions*****************************************
    //**************************************************************************************************
    //     @article{Lappi:2013zma,
    //         author         = "Lappi, T. and Mäntysaari, H.",
    //         title          = "{Single inclusive particle production at high energy from
    //                             HERA data to proton-nucleus collisions}",
    //         journal        = "Phys. Rev.",
    //         volume         = "D88",
    //         year           = "2013",
    //         pages          = "114020",
    //         doi            = "10.1103/PhysRevD.88.114020",
    //         eprint         = "1309.6963",
    //         archivePrefix  = "arXiv",
    //         primaryClass   = "hep-ph",
    //         SLACcitation   = "%%CITATION = ARXIV:1309.6963;%%"
    //     }

    ifstream inCGC;
    Int_t nlinesPi0CGC = 0;
    inCGC.open(fileNameEPS09sPi0CGC,ios_base::in);
    Double_t xESPsPi0CGC[100], yESPsPi0CGC[100];

    while(!inCGC.eof()){
        inCGC >> xESPsPi0CGC[nlinesPi0CGC]  >> yESPsPi0CGC[nlinesPi0CGC];
        cout << nlinesPi0CGC << "         "  << xESPsPi0CGC[nlinesPi0CGC] << "         "  <<yESPsPi0CGC[nlinesPi0CGC]<<endl;
        nlinesPi0CGC++;
    }
    inCGC.close();
    TGraph* graphPi0RpACGC = new TGraph(nlinesPi0CGC,xESPsPi0CGC,yESPsPi0CGC);

    //**************************************************************************************************
    //********************** extract EPOS3 spectra *****************************************************
    //**************************************************************************************************
    //     @article{Werner:2013tya,
    //       author         = "Werner, K. and Guiot, B. and Karpenko, Iu. and Pierog,
    //                         T.",
    //       title          = "{Analysing radial flow features in p-Pb and p-p
    //                         collisions at several TeV by studying identified particle
    //                         production in EPOS3}",
    //       journal        = "Phys. Rev.",
    //       volume         = "C89",
    //       year           = "2014",
    //       number         = "6",
    //       pages          = "064903",
    //       doi            = "10.1103/PhysRevC.89.064903",
    //       eprint         = "1312.1233",
    //       archivePrefix  = "arXiv",
    //       primaryClass   = "nucl-th",
    //       SLACcitation   = "%%CITATION = ARXIV:1312.1233;%%"
    //     }
    // file from Klaus Werner private communication
    Double_t ptBinningEPOS[59]      = { 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                        1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                                        2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
                                        3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8,
                                        5.0, 5.4, 5.8, 6.2, 6.6, 7.0, 7.5, 8.0, 8.5, 9.0,
                                        10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0};
    Int_t nBinsXEPOS                = 58;

    TFile* fileEPOS3                = new TFile("ExternalInputpPb/Theory/pi0_eta_EPOS3.root");
    TH1D* histoPi0EPOS3             = (TH1D*)fileEPOS3->Get("pi0_pt");
    TH1D* histoEtaEPOS3             = (TH1D*)fileEPOS3->Get("eta_pt");
    histoPi0EPOS3->Sumw2();
    histoPi0EPOS3->Scale(1/(2*TMath::Pi()));
    histoEtaEPOS3->Sumw2();
    histoEtaEPOS3->Scale(1/(2*TMath::Pi()));

    TH1D* histoPi0EPOS3RebPrep      = (TH1D*)histoPi0EPOS3->Clone("bla1");
    TH1D* histoEtaEPOS3RebPrep      = (TH1D*)histoEtaEPOS3->Clone("bla2");
    for (Int_t i = 1; i< histoPi0EPOS3RebPrep->GetNbinsX()+1; i++){
        histoPi0EPOS3RebPrep->SetBinContent(i, histoPi0EPOS3RebPrep->GetBinContent(i)*histoPi0EPOS3RebPrep->GetBinWidth(i));
        histoPi0EPOS3RebPrep->SetBinError(i, histoPi0EPOS3RebPrep->GetBinError(i)*histoPi0EPOS3RebPrep->GetBinWidth(i));
    }
    for (Int_t i = 1; i< histoEtaEPOS3RebPrep->GetNbinsX()+1; i++){
        histoEtaEPOS3RebPrep->SetBinContent(i, histoEtaEPOS3RebPrep->GetBinContent(i)*histoEtaEPOS3RebPrep->GetBinWidth(i));
        histoEtaEPOS3RebPrep->SetBinError(i, histoEtaEPOS3RebPrep->GetBinError(i)*histoEtaEPOS3RebPrep->GetBinWidth(i));
    }
    TH1D* histoPi0EPOS3Reb          = (TH1D*)histoPi0EPOS3RebPrep->Rebin(nBinsXEPOS,"EPOSpi0_Reb",ptBinningEPOS);
    TH1D* histoEtaEPOS3Reb          = (TH1D*)histoEtaEPOS3RebPrep->Rebin(nBinsXEPOS,"EPOSeta_Reb",ptBinningEPOS);

    for (Int_t i = 1; i< histoPi0EPOS3->GetNbinsX()+1; i++){
        histoPi0EPOS3->SetBinContent(i, histoPi0EPOS3->GetBinContent(i)/histoPi0EPOS3->GetBinCenter(i));
        histoPi0EPOS3->SetBinError(i, histoPi0EPOS3->GetBinError(i)/histoPi0EPOS3->GetBinCenter(i));
    }
    for (Int_t i = 1; i< histoEtaEPOS3->GetNbinsX()+1; i++){
        histoEtaEPOS3->SetBinContent(i, histoEtaEPOS3->GetBinContent(i)/histoEtaEPOS3->GetBinCenter(i));
        histoEtaEPOS3->SetBinError(i, histoEtaEPOS3->GetBinError(i)/histoEtaEPOS3->GetBinCenter(i));
    }
    for (Int_t i = 1; i< histoPi0EPOS3Reb->GetNbinsX()+1; i++){
        histoPi0EPOS3Reb->SetBinContent(i, histoPi0EPOS3Reb->GetBinContent(i)/(histoPi0EPOS3Reb->GetBinCenter(i)*histoPi0EPOS3Reb->GetBinWidth(i)));
        histoPi0EPOS3Reb->SetBinError(i, histoPi0EPOS3Reb->GetBinError(i)/(histoPi0EPOS3Reb->GetBinCenter(i)*histoPi0EPOS3Reb->GetBinWidth(i)));
    }
    for (Int_t i = 1; i< histoEtaEPOS3Reb->GetNbinsX()+1; i++){
        histoEtaEPOS3Reb->SetBinContent(i, histoEtaEPOS3Reb->GetBinContent(i)/(histoEtaEPOS3Reb->GetBinCenter(i)*histoEtaEPOS3Reb->GetBinWidth(i)));
        histoEtaEPOS3Reb->SetBinError(i, histoEtaEPOS3Reb->GetBinError(i)/(histoEtaEPOS3Reb->GetBinCenter(i)*histoEtaEPOS3Reb->GetBinWidth(i)));
    }

    TH1D* histoEtaToPi0EPOS3Reb     = (TH1D* )histoEtaEPOS3Reb->Clone("EtaToPi0EPOS3Reb");
    histoEtaToPi0EPOS3Reb->Sumw2();
    histoEtaToPi0EPOS3Reb->Divide(histoEtaToPi0EPOS3Reb,histoPi0EPOS3Reb);
    TH1D* histoEtaToPi0EPOS3        = (TH1D* )histoEtaEPOS3->Clone("EtaToPi0EPOS3");
    histoEtaToPi0EPOS3->Sumw2();
    histoEtaToPi0EPOS3->Divide(histoEtaToPi0EPOS3,histoPi0EPOS3);

    //**************************************************************************************************
    //********************** DPMJet MC spectra and ratio ***********************************************
    //**************************************************************************************************
    TString centrality[5]           = {"", "0-20%", "20-40%", "40-60%", "60-100%"};
    // file generated with TaskV1/ExtractMCInputSpectraFromFile.C++ based on PCM-EMC inputs
    TH1D* histoPi0DPMJet[5]         = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histoPi0DPMJetReb[5]      = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histoPiChDPMJet[5]        = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histoKChDPMJet[5]         = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histoEtaDPMJet[5]         = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histoEtaDPMJetReb[5]      = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histoEtaToPi0DPMJet[5]    = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histoPi0ToPiCHDPMJet[5]   = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histoEtaToKCHDPMJet[5]    = {NULL, NULL, NULL, NULL, NULL};
    TDirectory* directoryDPMJetCent[5]  = {NULL, NULL, NULL, NULL, NULL};

    TFile* fileDPMJet               = new TFile("ExternalInputpPb/Theory/MCInputCompilationLHC13b2_efix_pPb5TeV_2.root");
    histoPi0DPMJet[0]               = (TH1D*)fileDPMJet->Get("MC_Pi0_Pt");
    histoPi0DPMJetReb[0]            = (TH1D*)fileDPMJet->Get("MC_Pi0_Pt_Rebinned");
    histoPiChDPMJet[0]              = (TH1D*)fileDPMJet->Get("MC_PiCh_All_Pt");
    histoKChDPMJet[0]               = (TH1D*)fileDPMJet->Get("MC_KCh_All_Pt");
    histoEtaDPMJet[0]               = (TH1D*)fileDPMJet->Get("MC_Eta_Pt");
    histoEtaDPMJetReb[0]            = (TH1D*)fileDPMJet->Get("MC_Eta_Pt_Rebinned");
    histoEtaToPi0DPMJet[0]          = (TH1D*)fileDPMJet->Get("MCEtaToPi0");
    histoPi0ToPiCHDPMJet[0]         = (TH1D*)fileDPMJet->Get("MCPi0ToPiCh");
    histoEtaToKCHDPMJet[0]          = (TH1D*)fileDPMJet->Get("MCEtaToKCh");

    TFile* fileDPMJetCent           = new TFile("ExternalInputpPb/Theory/MCInputCompilationLHC13b2_efix_pPb5TeV_0.root");
    for (Int_t cent = 1; cent < 5; cent++){
        cout << Form("%spPb_5.023TeV",centrality[cent].Data()) << endl;
        directoryDPMJetCent[cent]       = (TDirectory*)fileDPMJetCent->Get(Form("%spPb_5.023TeV",centrality[cent].Data()));
        histoPi0DPMJet[cent]            = (TH1D*)directoryDPMJetCent[cent]->Get("MC_Pi0_Pt");
        histoPi0DPMJetReb[cent]         = (TH1D*)directoryDPMJetCent[cent]->Get("MC_Pi0_Pt_Rebinned");
        histoPiChDPMJet[cent]           = (TH1D*)directoryDPMJetCent[cent]->Get("MC_PiCh_All_Pt");
        histoKChDPMJet[cent]            = (TH1D*)directoryDPMJetCent[cent]->Get("MC_KCh_All_Pt");
        histoEtaDPMJet[cent]            = (TH1D*)directoryDPMJetCent[cent]->Get("MC_Eta_Pt");
        histoEtaDPMJetReb[cent]         = (TH1D*)directoryDPMJetCent[cent]->Get("MC_Eta_Pt_Rebinned");
        histoEtaToPi0DPMJet[cent]       = (TH1D*)directoryDPMJetCent[cent]->Get("MCEtaToPi0");
        histoPi0ToPiCHDPMJet[cent]      = (TH1D*)directoryDPMJetCent[cent]->Get("MCPi0ToPiCh");
        histoEtaToKCHDPMJet[cent]       = (TH1D*)directoryDPMJetCent[cent]->Get("MCEtaToKCh");
    }

    TFile* fileDPMJet2              = new TFile("ExternalInputpPb/Theory/DPMJet_pPb_5023GeV_990Mio.root");
    TH1D* histoPi0DPMJet2               = (TH1D*)fileDPMJet2->Get("hPtPi0DPMJetpPb");
    TH1D* histoPiPlDPMJet2              = (TH1D*)fileDPMJet2->Get("hPtPiPlDPMJetpPb");
    TH1D* histoPiMiDPMJet2              = (TH1D*)fileDPMJet2->Get("hPtPiMiDPMJetpPb");
    TH1D* histoPiChDPMJet2              = (TH1D*)histoPiPlDPMJet2->Clone("hPtPiChDPMJetpPb");
    histoPiChDPMJet2->Add(histoPiMiDPMJet2);
    histoPiChDPMJet2->Scale(1/2.);
    TH1D* histoEtaDPMJet2               = (TH1D*)fileDPMJet2->Get("hPtEtaDPMJetpPb");
    TH1D* histoEtaToPi0DPMJet2          = (TH1D*)histoEtaDPMJet2->Clone("MCEtaToPi0DPMJet");
    histoEtaToPi0DPMJet2->Divide(histoPi0DPMJet2);
    TH1D* histoPi0ToPiCHDPMJet2         = (TH1D*)histoPi0DPMJet2->Clone("MCPi0ToPiChDPMJet");
    histoPi0ToPiCHDPMJet2->Divide(histoPiChDPMJet2);

    //**************************************************************************************************
    //********************** EPOSLHC MC spectra and ratio ***********************************************
    //**************************************************************************************************
    // file generated with TaskV1/ExtractMCInputSpectraFromFile.C++ based on PCM inputs
    TH1D* histoPi0EPOSLHC[5]         = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histoPi0EPOSLHCReb[5]      = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histoPiChEPOSLHC[5]        = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histoKChEPOSLHC[5]         = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histoEtaEPOSLHC[5]         = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histoEtaEPOSLHCReb[5]      = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histoEtaToPi0EPOSLHC[5]    = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histoPi0ToPiCHEPOSLHC[5]   = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histoEtaToKCHEPOSLHC[5]    = {NULL, NULL, NULL, NULL, NULL};
    TDirectory* directoryEPOSLHCCent[5]  = {NULL, NULL, NULL, NULL, NULL};

    TFile* fileEPOSLHCCent           = new TFile("ExternalInputpPb/Theory/MCInputCompilationLHC17f2a_pPb5TeV_0.root");
    for (Int_t cent = 0; cent < 5; cent++){
        cout << Form("%spPb_5.023TeV",centrality[cent].Data()) << endl;
        directoryEPOSLHCCent[cent]       = (TDirectory*)fileEPOSLHCCent->Get(Form("%spPb_5.023TeV",centrality[cent].Data()));
        histoPi0EPOSLHC[cent]            = (TH1D*)directoryEPOSLHCCent[cent]->Get("MC_Pi0_Pt");
        histoPi0EPOSLHCReb[cent]         = (TH1D*)directoryEPOSLHCCent[cent]->Get("MC_Pi0_Pt_Rebinned");
        histoPiChEPOSLHC[cent]           = (TH1D*)directoryEPOSLHCCent[cent]->Get("MC_PiCh_All_Pt");
        histoKChEPOSLHC[cent]            = (TH1D*)directoryEPOSLHCCent[cent]->Get("MC_KCh_All_Pt");
        histoEtaEPOSLHC[cent]            = (TH1D*)directoryEPOSLHCCent[cent]->Get("MC_Eta_Pt");
        histoEtaEPOSLHCReb[cent]         = (TH1D*)directoryEPOSLHCCent[cent]->Get("MC_Eta_Pt_Rebinned");
        histoEtaToPi0EPOSLHC[cent]       = (TH1D*)directoryEPOSLHCCent[cent]->Get("MCEtaToPi0");
        histoPi0ToPiCHEPOSLHC[cent]      = (TH1D*)directoryEPOSLHCCent[cent]->Get("MCPi0ToPiCh");
        histoEtaToKCHEPOSLHC[cent]       = (TH1D*)directoryEPOSLHCCent[cent]->Get("MCEtaToKCh");
    }

    TFile* fileEPOSLHC2              = new TFile("ExternalInputpPb/Theory/EPOSLHC_pPb_5023GeV_99Mio.root");
    TH1D* histoPi0EPOSLHC2               = (TH1D*)fileEPOSLHC2->Get("hPtPi0EPOSLHCpPb");
    TH1D* histoPiPlEPOSLHC2              = (TH1D*)fileEPOSLHC2->Get("hPtPiPlEPOSLHCpPb");
    TH1D* histoPiMiEPOSLHC2              = (TH1D*)fileEPOSLHC2->Get("hPtPiMiEPOSLHCpPb");
    TH1D* histoPiChEPOSLHC2              = (TH1D*)histoPiPlEPOSLHC2->Clone("hPtPiChEPOSLHCpPb");
    histoPiChEPOSLHC2->Add(histoPiMiEPOSLHC2);
    histoPiChEPOSLHC2->Scale(1/2.);
    TH1D* histoEtaEPOSLHC2               = (TH1D*)fileEPOSLHC2->Get("hPtEtaEPOSLHCpPb");
    TH1D* histoEtaToPi0EPOSLHC2          = (TH1D*)histoEtaEPOSLHC2->Clone("MCEtaToPi0EPOSLHC");
    histoEtaToPi0EPOSLHC2->Divide(histoPi0EPOSLHC2);
    TH1D* histoPi0ToPiCHEPOSLHC2         = (TH1D*)histoPi0EPOSLHC2->Clone("MCPi0ToPiChEPOSLHC");
    histoPi0ToPiCHEPOSLHC2->Divide(histoPiChEPOSLHC2);


    //**************************************************************************************************
    //********************** HIJING MC spectra and ratio ***********************************************
    //**************************************************************************************************
    // file generated with TaskV1/ExtractMCInputSpectraFromFile.C++ based on PCM-EMC inputs
    TFile* fileHIJING               = new TFile("ExternalInputpPb/Theory/MCInputCompilationLHC13e7_pPb5TeV_2.root");
    TH1D* histoPi0HIJING            = (TH1D*)fileHIJING->Get("MC_Pi0_Pt");
    TH1D* histoPi0HIJINGReb         = (TH1D*)fileHIJING->Get("MC_Pi0_Pt_Rebinned");
    TH1D* histoPiChHIJING           = (TH1D*)fileHIJING->Get("MC_PiCh_All_Pt");
    TH1D* histoKChHIJING            = (TH1D*)fileHIJING->Get("MC_KCh_All_Pt");
    TH1D* histoEtaHIJING            = (TH1D*)fileHIJING->Get("MC_Eta_Pt");
    TH1D* histoEtaHIJINGReb         = (TH1D*)fileHIJING->Get("MC_Eta_Pt_Rebinned");
    TH1D* histoEtaToPi0HIJING       = (TH1D*)fileHIJING->Get("MCEtaToPi0");
    TH1D* histoEtaToKCHHIJING       = (TH1D*)fileHIJING->Get("MCEtaToKCh");
    TH1D* histoPi0ToPiCHHIJING      = (TH1D*)fileHIJING->Get("MCPi0ToPiCh");

    TFile* fileHIJING2              = new TFile("ExternalInputpPb/Theory/HIJING_pPb_5023GeV_97Mio.root");
    TH1D* histoPi0HIJING2               = (TH1D*)fileHIJING2->Get("hPtPi0HIJINGpPb");
    TH1D* histoPiPlHIJING2              = (TH1D*)fileHIJING2->Get("hPtPiPlHIJINGpPb");
    TH1D* histoPiMiHIJING2              = (TH1D*)fileHIJING2->Get("hPtPiMiHIJINGpPb");
    TH1D* histoPiChHIJING2              = (TH1D*)histoPiPlHIJING2->Clone("hPtPiChHIJINGpPb");
    histoPiChHIJING2->Add(histoPiMiHIJING2);
    histoPiChHIJING2->Scale(1/2.);
    TH1D* histoEtaHIJING2               = (TH1D*)fileHIJING2->Get("hPtEtaHIJINGpPb");
    TH1D* histoEtaToPi0HIJING2          = (TH1D*)histoEtaHIJING2->Clone("MCEtaToPi0HIJING");
    histoEtaToPi0HIJING2->Divide(histoPi0HIJING2);
    TH1D* histoPi0ToPiCHHIJING2         = (TH1D*)histoPi0HIJING2->Clone("MCPi0ToPiChHIJING");
    histoPi0ToPiCHHIJING2->Divide(histoPiChHIJING2);

    //**************************************************************************************************
    //****************************** extracting McGill predictions**************************************
    //**************************************************************************************************
//     @article{Shen:2016zpp,
//         author         = "Shen, Chun and Paquet, Jean-François and Denicol,
//                             Gabriel S. and Jeon, Sangyong and Gale, Charles",
//         title          = "{Collectivity and electromagnetic radiation in small
//                             systems}",
//         journal        = "Phys. Rev.",
//         volume         = "C95",
//         year           = "2017",
//         pages          = "014906",
//         doi            = "10.1103/PhysRevC.95.014906",
//         eprint         = "1609.02590",
//         archivePrefix  = "arXiv",
//         primaryClass   = "nucl-th",
//         SLACcitation   = "%%CITATION = ARXIV:1609.02590;%%"
//     }

    // *************************************
    // read eta spectra
    ifstream inMCGillEta;
    Int_t nlinesEtaMCGill = 0;
    inMCGillEta.open("ExternalInputpPb/Theory/McGill/spectra_eta_hydro.dat",ios_base::in);
    Double_t xPtEtaMCGill[100], xPtErrEtaMCGill[0], yYieldEtaMCGill[100], yYieldErrEtaMCGill[100];

    while(!inMCGillEta.eof()){
        inMCGillEta >> xPtEtaMCGill[nlinesEtaMCGill]  >> yYieldEtaMCGill[nlinesEtaMCGill] >> yYieldErrEtaMCGill[nlinesEtaMCGill];
        cout << nlinesEtaMCGill << "         "  << xPtEtaMCGill[nlinesEtaMCGill] << "         "  <<yYieldEtaMCGill[nlinesEtaMCGill]<<endl;
        xPtErrEtaMCGill[nlinesEtaMCGill] = 0.1;
        nlinesEtaMCGill++;

    }
    inMCGillEta.close();
    TGraphErrors* graphEtaSpecMcGill = new TGraphErrors(nlinesEtaMCGill-1,xPtEtaMCGill,yYieldEtaMCGill, xPtErrEtaMCGill, yYieldErrEtaMCGill );

    // *************************************
    // read eta v2
    ifstream inMCGillEtaV2;
    nlinesEtaMCGill = 0;
    inMCGillEtaV2.open("ExternalInputpPb/Theory/McGill/v2_eta_hydro.dat",ios_base::in);
    Double_t yV2EtaMCGill[100], yV2ErrEtaMCGill[100];

    while(!inMCGillEtaV2.eof()){
        inMCGillEtaV2 >> xPtEtaMCGill[nlinesEtaMCGill]  >> yV2EtaMCGill[nlinesEtaMCGill] >> yV2ErrEtaMCGill[nlinesEtaMCGill];
        cout << nlinesEtaMCGill << "         "  << xPtEtaMCGill[nlinesEtaMCGill] << "         "  <<yV2EtaMCGill[nlinesEtaMCGill]<<endl;
        xPtErrEtaMCGill[nlinesEtaMCGill] = 0.1;
        nlinesEtaMCGill++;

    }
    inMCGillEtaV2.close();
    TGraphErrors* graphEtaV2McGill = new TGraphErrors(nlinesEtaMCGill-1,xPtEtaMCGill,yV2EtaMCGill, xPtErrEtaMCGill, yV2ErrEtaMCGill );

    // *************************************
    // read pi0 spectra
    ifstream inMCGillPi0;
    Int_t nlinesPi0MCGill = 0;
    inMCGillPi0.open("ExternalInputpPb/Theory/McGill/spectra_pion_p_hydro.dat",ios_base::in);
    Double_t xPtPi0MCGill[100], xPtErrPi0MCGill[0], yYieldPi0MCGill[100], yYieldErrPi0MCGill[100];

    while(!inMCGillPi0.eof()){
        inMCGillPi0 >> xPtPi0MCGill[nlinesPi0MCGill]  >> yYieldPi0MCGill[nlinesPi0MCGill] >> yYieldErrPi0MCGill[nlinesPi0MCGill];
        cout << nlinesPi0MCGill << "         "  << xPtPi0MCGill[nlinesPi0MCGill] << "         "  <<yYieldPi0MCGill[nlinesPi0MCGill]<<endl;
        xPtErrPi0MCGill[nlinesPi0MCGill] = 0.1;
        nlinesPi0MCGill++;

    }
    inMCGillPi0.close();
    TGraphErrors* graphPi0SpecMcGill = new TGraphErrors(nlinesPi0MCGill-1,xPtPi0MCGill,yYieldPi0MCGill, xPtErrPi0MCGill, yYieldErrPi0MCGill );
    // *************************************
    // read pi0 v2
    ifstream inMCGillPi0V2;
    nlinesPi0MCGill = 0;
    inMCGillPi0V2.open("ExternalInputpPb/Theory/McGill/v2_pion_p_hydro.dat",ios_base::in);
    Double_t yV2Pi0MCGill[100], yV2ErrPi0MCGill[100];

    while(!inMCGillPi0V2.eof()){
        inMCGillPi0V2 >> xPtPi0MCGill[nlinesPi0MCGill]  >> yV2Pi0MCGill[nlinesPi0MCGill] >> yV2ErrPi0MCGill[nlinesPi0MCGill];
        cout << nlinesPi0MCGill << "         "  << xPtPi0MCGill[nlinesPi0MCGill] << "         "  <<yV2Pi0MCGill[nlinesPi0MCGill]<<endl;
        xPtErrPi0MCGill[nlinesPi0MCGill] = 0.1;
        nlinesPi0MCGill++;

    }
    inMCGillPi0V2.close();
    TGraphErrors* graphPi0V2McGill = new TGraphErrors(nlinesPi0MCGill-1,xPtPi0MCGill,yV2Pi0MCGill, xPtErrPi0MCGill, yV2ErrPi0MCGill );

    // *************************************
    // calculate eta/pi0
    TGraphErrors* graphEtaToPi0McGill = CalculateGraphRatioToGraphErrs(graphEtaSpecMcGill, graphPi0SpecMcGill);

    //**************************************************************************************************
    //********************* extracting McGill predictions including other observables ******************
    //**************************************************************************************************
    const Int_t nParticles              = 14;
    const Int_t nMeasurements           = 15;
    const Int_t nCent                   = 3;
    TString centNames[nCent]            = {"0-100", "0-5", "5-10"};
    TString centNamesOut[nCent]         = {"", "0-5_", "5-10_"};
    TString particleNames[nParticles]   = {"pion_p", "pion_m", "kaon_p", "kaon_m", "proton",
                                           "anti_proton", "charged_hadron", "Lambda", "anti_Lambda", "Omega",
                                           "anti_Omega", "Xi_m", "anti_Xi_p", "phi"};
    TString particleNamesOut[nParticles]= {"PiPlus", "PiMinus", "KaonPlus", "KaonMinus", "Proton",
                                           "AntiProton", "ChargedHadron", "Lambda", "AntiLambda", "Omega",
                                           "AntiOmega", "Xi", "AntiXi", "Phi"};
    vector<Double_t> ***valuesMcGill    = new vector<Double_t>**[nCent];     // iCent x iParticle x nMeasurements matrix with theory curves
    for (Int_t iCent = 0; iCent < nCent; iCent++){
        valuesMcGill[iCent]             = new vector<Double_t>*[nParticles];
        for(Int_t iParticle=0; iParticle<nParticles; iParticle++){
            valuesMcGill[iCent][iParticle]  = new vector<Double_t>[nMeasurements];
        }
    }
    Int_t nPtPoint[nCent][nParticles];


    //  read from file
    for (Int_t iCent = 0; iCent < nCent; iCent++){
        cout << "reading cent: " << centNames[iCent].Data() << endl;
        for(Int_t iParticle=0; iParticle<nParticles; iParticle++){
            nPtPoint[iCent][iParticle]      = 0;
            ifstream fileMcGillInput;
            TString fileName                = Form("ExternalInputpPb/Theory/McGill/MCGlbpPb5020_C%s_vis010_Ttr155_new/%s_differential_observables_ALICE.dat", centNames[iCent].Data(), particleNames[iParticle].Data());
            fileMcGillInput.open(fileName.Data(),ios_base::in);
            cout << "opening: " << fileName.Data() << endl;
            Int_t iPtCurrent    = 0;
            Int_t nCurrMeas     = 0;
            string line;
            while (getline(fileMcGillInput, line) && iPtCurrent < 100) {
                TString temp        = "";
                TString tempBin     = "";
                istringstream cs(line); // controll stream
                cs >> temp;

                if (!temp.Contains("#")){
                    istringstream ss(line);
                    Int_t iMeasurement  = 0;
                    while(ss && iMeasurement < nMeasurements){
                        ss >> temp;
                        if(!temp.IsNull() && temp.CompareTo("nan") != 0){
                            valuesMcGill[iCent][iParticle][iMeasurement].push_back(temp.Atof());
                            cout << temp.Data() << "\t ";
                            iMeasurement++;
                        } else {
                            valuesMcGill[iCent][iParticle][iMeasurement].push_back(-1.0e-6);
                            cout << -1.0e-6 << "\t ";
                            iMeasurement++;
                        }
                    }
                    nCurrMeas = iMeasurement;
                    cout << endl;
                    iPtCurrent++;
                } else {
                    cout << "first line contains comments" << endl;
                }
            }
            cout << "Number of pT bins: "<< iPtCurrent << "\t number of measurement points: "<< nCurrMeas <<  endl;
            nPtPoint[iCent][iParticle]         = iPtCurrent;
            fileMcGillInput.close();
        }
    }

    TGraphErrors* graphSpectraMcGill[nCent][nParticles];
    TGraphErrors* graphVnMcGill[nCent][nParticles][6];
    for (Int_t iCent = 0; iCent < nCent; iCent++){
        for(Int_t iParticle=0; iParticle<nParticles; iParticle++){
            graphSpectraMcGill[iCent][iParticle]       = NULL;
            for (Int_t iVn = 0; iVn < 6; iVn++){
                graphVnMcGill[iCent][iParticle][iVn]   = NULL;
            }
        }
    }
    for (Int_t iCent = 0; iCent < nCent; iCent++){
        for(Int_t iParticle=0; iParticle<nParticles; iParticle++){
            graphSpectraMcGill[iCent][iParticle]       = new TGraphErrors(nPtPoint[iCent][iParticle]);
            for (Int_t iPt = 0; iPt < nPtPoint[iCent][iParticle]; iPt++){
                graphSpectraMcGill[iCent][iParticle]->SetPoint(iPt, valuesMcGill[iCent][iParticle][0].at(iPt), valuesMcGill[iCent][iParticle][1].at(iPt));
                graphSpectraMcGill[iCent][iParticle]->SetPointError(iPt, 0.01, valuesMcGill[iCent][iParticle][2].at(iPt) );
            }
            graphSpectraMcGill[iCent][iParticle]->SetName(Form("graph%sSpecMcGill5023GeV", particleNamesOut[iParticle].Data()));
            for (Int_t iVn = 0; iVn < 6; iVn++){
                graphVnMcGill[iCent][iParticle][iVn]   = new TGraphErrors(nPtPoint[iCent][iParticle]);
                for (Int_t iPt = 0; iPt < nPtPoint[iCent][iParticle]; iPt++){
                    graphVnMcGill[iCent][iParticle][iVn]->SetPoint(iPt, valuesMcGill[iCent][iParticle][0].at(iPt), valuesMcGill[iCent][iParticle][3+2*iVn].at(iPt));
                    graphVnMcGill[iCent][iParticle][iVn]->SetPointError(iPt, 0.01, valuesMcGill[iCent][iParticle][4+2*iVn].at(iPt) );
                }
                graphVnMcGill[iCent][iParticle][iVn]->SetName(Form("graph%sV%iMcGill5023GeV", particleNamesOut[iParticle].Data(),iVn+1));
            }
        }
    }

    //**********************************************************************************************************************
    //********************************pp 5.023 TeV Pi0 and Eta calc*********************************************************
    //**********************************************************************************************************************

    //*********************************************************
    // eta Vogelsang - FF: AESSS, PDF: CT10
    Double_t ptNLOEta5023GeV[200];
    Double_t muHalfEta5023GeV[200];
    Double_t muOneEta5023GeV[200];
    Double_t muTwoEta5023GeV[200];
    Int_t nlinesNLOEta5023GeV                       = 0;
    Bool_t fillAllMuScalesEta5023GeV                = kTRUE;

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
    TGraph* graphNLOCalcInvYieldEtaMuHalfpPb5023GeV = NULL;
    TGraph* graphNLOCalcInvYieldEtaMuOnepPb5023GeV  = NULL;
    TGraph* graphNLOCalcInvYieldEtaMuTwopPb5023GeV  = NULL;
    TGraph* graphNLOCalcInvYieldEtaMuHalfpp5023GeV = NULL;
    TGraph* graphNLOCalcInvYieldEtaMuOnepp5023GeV  = NULL;
    TGraph* graphNLOCalcInvYieldEtaMuTwopp5023GeV  = NULL;

    // fill x-section graph for default mu
    graphNLOCalcInvSecEtaMuOne5023GeV               = new TGraph(nlinesNLOEta5023GeV,ptNLOEta5023GeV,muOneEta5023GeV);
    graphNLOCalcInvSecEtaMuOne5023GeV->Print();
    // fill yield graph for default mu
    graphNLOCalcInvYieldEtaMuOnepPb5023GeV          = ScaleGraph(graphNLOCalcInvSecEtaMuOne5023GeV, 1/(xSection5023GeVINEL*recalcBarn)*ncollpPb5023GeV);
    graphNLOCalcInvYieldEtaMuOnepp5023GeV           = ScaleGraph(graphNLOCalcInvSecEtaMuOne5023GeV, 1/(xSection5023GeVINEL*recalcBarn));

    if (fillAllMuScalesEta5023GeV){
        // fill x-section graphs for mu variations
        graphNLOCalcInvSecEtaMuHalf5023GeV          = new TGraph(nlinesNLOEta5023GeV,ptNLOEta5023GeV,muHalfEta5023GeV);
        graphNLOCalcInvSecEtaMuHalf5023GeV->RemovePoint(0);
        graphNLOCalcInvSecEtaMuTwo5023GeV           = new TGraph(nlinesNLOEta5023GeV,ptNLOEta5023GeV,muTwoEta5023GeV);
        // fill yield graphs for mu variations
        graphNLOCalcInvYieldEtaMuHalfpPb5023GeV     = ScaleGraph(graphNLOCalcInvSecEtaMuHalf5023GeV, 1/(xSection5023GeVINEL*recalcBarn)*ncollpPb5023GeV);
        graphNLOCalcInvYieldEtaMuTwopPb5023GeV      = ScaleGraph(graphNLOCalcInvSecEtaMuTwo5023GeV, 1/(xSection5023GeVINEL*recalcBarn)*ncollpPb5023GeV);
        graphNLOCalcInvYieldEtaMuHalfpp5023GeV      = ScaleGraph(graphNLOCalcInvSecEtaMuHalf5023GeV, 1/(xSection5023GeVINEL*recalcBarn));
        graphNLOCalcInvYieldEtaMuTwopp5023GeV       = ScaleGraph(graphNLOCalcInvSecEtaMuTwo5023GeV, 1/(xSection5023GeVINEL*recalcBarn));
    }
    // combine all mu scales into 1 graph for eventual plotting
    TGraphAsymmErrors* graphNLOCalcInvSecEta5023GeV     = CombineMuScales(nlinesNLOEta5023GeV, ptNLOEta5023GeV, muOneEta5023GeV, muHalfEta5023GeV, muTwoEta5023GeV);
    TGraphAsymmErrors* graphNLOCalcInvYieldEtapPb5023GeV= ScaleGraphAsym(graphNLOCalcInvSecEta5023GeV, 1/(xSection5023GeVINEL*recalcBarn)*ncollpPb5023GeV);
    TGraphAsymmErrors* graphNLOCalcInvYieldEtapp5023GeV = ScaleGraphAsym(graphNLOCalcInvSecEta5023GeV, 1/(xSection5023GeVINEL*recalcBarn));

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
    TGraph* graphNLOCalcInvYieldPi0MuHalfpPb5023GeV     = ScaleGraph(graphNLOCalcInvSecPi0MuHalf5023GeV, 1/(xSection5023GeVINEL*recalcBarn)*ncollpPb5023GeV);
    TGraph* graphNLOCalcInvYieldPi0MuOnepPb5023GeV      = ScaleGraph(graphNLOCalcInvSecPi0MuOne5023GeV, 1/(xSection5023GeVINEL*recalcBarn)*ncollpPb5023GeV);
    TGraph* graphNLOCalcInvYieldPi0MuTwopPb5023GeV      = ScaleGraph(graphNLOCalcInvSecPi0MuTwo5023GeV, 1/(xSection5023GeVINEL*recalcBarn)*ncollpPb5023GeV);
    TGraph* graphNLOCalcInvYieldPi0MuHalfpp5023GeV      = ScaleGraph(graphNLOCalcInvSecPi0MuHalf5023GeV, 1/(xSection5023GeVINEL*recalcBarn));
    TGraph* graphNLOCalcInvYieldPi0MuOnepp5023GeV       = ScaleGraph(graphNLOCalcInvSecPi0MuOne5023GeV, 1/(xSection5023GeVINEL*recalcBarn));
    TGraph* graphNLOCalcInvYieldPi0MuTwopp5023GeV       = ScaleGraph(graphNLOCalcInvSecPi0MuTwo5023GeV, 1/(xSection5023GeVINEL*recalcBarn));
    // combine all mu scales into 1 graph for eventual plotting
    TGraphAsymmErrors* graphNLOCalcInvSecPi05023GeV     = CombineMuScales(nlinesNLOPi05023GeV, ptNLOPi05023GeV, muOnePi05023GeV, muHalfPi05023GeV, muTwoPi05023GeV);
    TGraphAsymmErrors* graphNLOCalcInvYieldPi0pPb5023GeV= ScaleGraphAsym(graphNLOCalcInvSecPi05023GeV, 1/(xSection5023GeVINEL*recalcBarn)*ncollpPb5023GeV);
    TGraphAsymmErrors* graphNLOCalcInvYieldPi0pp5023GeV = ScaleGraphAsym(graphNLOCalcInvSecPi05023GeV, 1/(xSection5023GeVINEL*recalcBarn));
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


    //*********************************************************
    // eta Vogelsang - FF: AESSS, nPDF: nCTEQ
    Double_t ptNLOEtanPDF5023GeV[200];
    Double_t muHalfEtanPDF5023GeV[200];
    Double_t muOneEtanPDF5023GeV[200];
    Double_t muTwoEtanPDF5023GeV[200];
    Int_t nlinesNLOEtanPDF5023GeV                       = 0;
    Bool_t fillAllMuScalesEtanPDF5023GeV                = kTRUE;

    TString fileNameNLOEtanPDF5023GeV                   = "ExternalInputpPb/Theory/Vogelsang/ALICENLOcalcEtaVogelsang5023GeV_AESSS_nCTEQ.dat"; // currently only mu = 1pt
    ifstream  fileNLOEtanPDF5023GeV;
    fileNLOEtanPDF5023GeV.open(fileNameNLOEtanPDF5023GeV,ios_base::in);
    cout << fileNameNLOEtanPDF5023GeV << endl;

    // read eta input theory file for 5.023TeV
    while(!fileNLOEtanPDF5023GeV.eof() && nlinesNLOEtanPDF5023GeV < 200){
        if (fillAllMuScalesEtanPDF5023GeV){
            fileNLOEtanPDF5023GeV >> ptNLOEtanPDF5023GeV[nlinesNLOEtanPDF5023GeV] >> muHalfEtanPDF5023GeV[nlinesNLOEtanPDF5023GeV] >> muOneEtanPDF5023GeV[nlinesNLOEtanPDF5023GeV] >> muTwoEtanPDF5023GeV[nlinesNLOEtanPDF5023GeV];
        } else {
            fileNLOEtanPDF5023GeV >> ptNLOEtanPDF5023GeV[nlinesNLOEtanPDF5023GeV] >> muOneEtanPDF5023GeV[nlinesNLOEtanPDF5023GeV] ;
            muHalfEtanPDF5023GeV[nlinesNLOEtanPDF5023GeV]   = muOneEtanPDF5023GeV[nlinesNLOEtanPDF5023GeV] ;
            muTwoEtanPDF5023GeV[nlinesNLOEtanPDF5023GeV]    = muOneEtanPDF5023GeV[nlinesNLOEtanPDF5023GeV] ;
        }
        nlinesNLOEtanPDF5023GeV++;
    }
    fileNLOEtanPDF5023GeV.close();
    // reduce number of lines by 1
    nlinesNLOEtanPDF5023GeV--;

    // generate graphs
    TGraph* graphNLOCalcInvSecEtanPDFMuHalf5023GeV      = NULL;
    TGraph* graphNLOCalcInvSecEtanPDFMuOne5023GeV       = NULL;
    TGraph* graphNLOCalcInvSecEtanPDFMuTwo5023GeV       = NULL;
    TGraph* graphNLOCalcInvYieldEtanPDFMuHalfpPb5023GeV = NULL;
    TGraph* graphNLOCalcInvYieldEtanPDFMuOnepPb5023GeV  = NULL;
    TGraph* graphNLOCalcInvYieldEtanPDFMuTwopPb5023GeV  = NULL;

    // fill x-section graph for default mu
    graphNLOCalcInvSecEtanPDFMuOne5023GeV               = new TGraph(nlinesNLOEtanPDF5023GeV,ptNLOEtanPDF5023GeV,muOneEtanPDF5023GeV);
    graphNLOCalcInvSecEtanPDFMuOne5023GeV->Print();
    // fill yield graph for default mu
    graphNLOCalcInvYieldEtanPDFMuOnepPb5023GeV          = ScaleGraph(graphNLOCalcInvSecEtanPDFMuOne5023GeV, 1/(xSection5023GeVINEL*recalcBarn)*ncollpPb5023GeV);

    if (fillAllMuScalesEtanPDF5023GeV){
        // fill x-section graphs for mu variations
        graphNLOCalcInvSecEtanPDFMuHalf5023GeV          = new TGraph(nlinesNLOEtanPDF5023GeV,ptNLOEtanPDF5023GeV,muHalfEtanPDF5023GeV);
        graphNLOCalcInvSecEtanPDFMuHalf5023GeV->RemovePoint(0);
        graphNLOCalcInvSecEtanPDFMuTwo5023GeV           = new TGraph(nlinesNLOEtanPDF5023GeV,ptNLOEtanPDF5023GeV,muTwoEtanPDF5023GeV);
        // fill yield graphs for mu variations
        graphNLOCalcInvYieldEtanPDFMuHalfpPb5023GeV     = ScaleGraph(graphNLOCalcInvSecEtanPDFMuHalf5023GeV, 1/(xSection5023GeVINEL*recalcBarn)*ncollpPb5023GeV);
        graphNLOCalcInvYieldEtanPDFMuTwopPb5023GeV      = ScaleGraph(graphNLOCalcInvSecEtanPDFMuTwo5023GeV, 1/(xSection5023GeVINEL*recalcBarn)*ncollpPb5023GeV);

    }
    // combine all mu scales into 1 graph for eventual plotting
    TGraphAsymmErrors* graphNLOCalcInvSecEtanPDF5023GeV     = CombineMuScales(nlinesNLOEtanPDF5023GeV, ptNLOEtanPDF5023GeV, muOneEtanPDF5023GeV, muHalfEtanPDF5023GeV, muTwoEtanPDF5023GeV);
    TGraphAsymmErrors* graphNLOCalcInvYieldEtanPDFpPb5023GeV= ScaleGraphAsym(graphNLOCalcInvSecEtanPDF5023GeV, 1/(xSection5023GeVINEL*recalcBarn)*ncollpPb5023GeV);

    //*********************************************************
    // pi0 Vogelsang - FF: DSS14, PDF:  nPDF: nCTEQ
    Double_t ptNLOPi0nPDF5023GeV[200];
    Double_t muHalfPi0nPDF5023GeV[200];
    Double_t muOnePi0nPDF5023GeV[200];
    Double_t muTwoPi0nPDF5023GeV[200];
    Int_t nlinesNLOPi0nPDF5023GeV                       = 0;
    Bool_t fillAllMuScalesPi0nPDF5023GeV                = kTRUE;

    TString fileNameNLOPi0nPDF5023GeV                   = "ExternalInputpPb/Theory/Vogelsang/ALICENLOcalcPi0Vogelsang5023GeV_DSS14_nCTEQ.dat";
    ifstream  fileNLOPi0nPDF5023GeV;
    fileNLOPi0nPDF5023GeV.open(fileNameNLOPi0nPDF5023GeV,ios_base::in);
    cout << fileNameNLOPi0nPDF5023GeV << endl;

    // read eta input theory file for 5.023TeV
    while(!fileNLOPi0nPDF5023GeV.eof() && nlinesNLOPi0nPDF5023GeV < 200){
        if (fillAllMuScalesPi0nPDF5023GeV){
            fileNLOPi0nPDF5023GeV >> ptNLOPi0nPDF5023GeV[nlinesNLOPi0nPDF5023GeV] >> muHalfPi0nPDF5023GeV[nlinesNLOPi0nPDF5023GeV] >> muOnePi0nPDF5023GeV[nlinesNLOPi0nPDF5023GeV] >> muTwoPi0nPDF5023GeV[nlinesNLOPi0nPDF5023GeV];
        } else {
            fileNLOPi0nPDF5023GeV >> ptNLOPi0nPDF5023GeV[nlinesNLOPi0nPDF5023GeV] >> muOnePi0nPDF5023GeV[nlinesNLOPi0nPDF5023GeV] ;
            muHalfPi0nPDF5023GeV[nlinesNLOPi0nPDF5023GeV]   = muOnePi0nPDF5023GeV[nlinesNLOPi0nPDF5023GeV] ;
            muTwoPi0nPDF5023GeV[nlinesNLOPi0nPDF5023GeV]    = muOnePi0nPDF5023GeV[nlinesNLOPi0nPDF5023GeV] ;
        }
        nlinesNLOPi0nPDF5023GeV++;
    }
    fileNLOPi0nPDF5023GeV.close();
    // reduce number of lines by 1
    nlinesNLOPi0nPDF5023GeV--;

    // generate graphs
    TGraph* graphNLOCalcInvSecPi0nPDFMuHalf5023GeV      = NULL;
    TGraph* graphNLOCalcInvSecPi0nPDFMuOne5023GeV       = NULL;
    TGraph* graphNLOCalcInvSecPi0nPDFMuTwo5023GeV       = NULL;
    TGraph* graphNLOCalcInvYieldPi0nPDFMuHalfpPb5023GeV = NULL;
    TGraph* graphNLOCalcInvYieldPi0nPDFMuOnepPb5023GeV  = NULL;
    TGraph* graphNLOCalcInvYieldPi0nPDFMuTwopPb5023GeV  = NULL;

    // fill x-section graph for default mu
    graphNLOCalcInvSecPi0nPDFMuOne5023GeV               = new TGraph(nlinesNLOPi0nPDF5023GeV,ptNLOPi0nPDF5023GeV,muOnePi0nPDF5023GeV);
    graphNLOCalcInvSecPi0nPDFMuOne5023GeV->Print();
    // fill yield graph for default mu
    graphNLOCalcInvYieldPi0nPDFMuOnepPb5023GeV          = ScaleGraph(graphNLOCalcInvSecPi0nPDFMuOne5023GeV, 1/(xSection5023GeVINEL*recalcBarn)*ncollpPb5023GeV);

    if (fillAllMuScalesPi0nPDF5023GeV){
        // fill x-section graphs for mu variations
        graphNLOCalcInvSecPi0nPDFMuHalf5023GeV          = new TGraph(nlinesNLOPi0nPDF5023GeV,ptNLOPi0nPDF5023GeV,muHalfPi0nPDF5023GeV);
        graphNLOCalcInvSecPi0nPDFMuHalf5023GeV->RemovePoint(0);
        graphNLOCalcInvSecPi0nPDFMuTwo5023GeV           = new TGraph(nlinesNLOPi0nPDF5023GeV,ptNLOPi0nPDF5023GeV,muTwoPi0nPDF5023GeV);
        // fill yield graphs for mu variations
        graphNLOCalcInvYieldPi0nPDFMuHalfpPb5023GeV     = ScaleGraph(graphNLOCalcInvSecPi0nPDFMuHalf5023GeV, 1/(xSection5023GeVINEL*recalcBarn)*ncollpPb5023GeV);
        graphNLOCalcInvYieldPi0nPDFMuTwopPb5023GeV      = ScaleGraph(graphNLOCalcInvSecPi0nPDFMuTwo5023GeV, 1/(xSection5023GeVINEL*recalcBarn)*ncollpPb5023GeV);

    }
    // combine all mu scales into 1 graph for eventual plotting
    TGraphAsymmErrors* graphNLOCalcInvSecPi0nPDF5023GeV     = CombineMuScales(nlinesNLOPi0nPDF5023GeV, ptNLOPi0nPDF5023GeV, muOnePi0nPDF5023GeV, muHalfPi0nPDF5023GeV, muTwoPi0nPDF5023GeV);
    TGraphAsymmErrors* graphNLOCalcInvYieldPi0nPDFpPb5023GeV= ScaleGraphAsym(graphNLOCalcInvSecPi0nPDF5023GeV, 1/(xSection5023GeVINEL*recalcBarn)*ncollpPb5023GeV);

    //*********************************************************
    // Eta/Pi0 Vogelsang FF pi0: DSS14, FF eta: AESSS,  nPDF: nCTEQ

    Double_t* valueNLOMuHalfEtanPDF5023GeV              = NULL;
    Double_t* valueNLOMuOneEtanPDF5023GeV               = graphNLOCalcInvSecEtanPDFMuOne5023GeV->GetY();
    Double_t* valueNLOMuTwoEtanPDF5023GeV               = NULL;
    if (fillAllMuScalesEtanPDF5023GeV){
        valueNLOMuHalfEtanPDF5023GeV                    = graphNLOCalcInvSecEtanPDFMuHalf5023GeV->GetY();
        valueNLOMuTwoEtanPDF5023GeV                     = graphNLOCalcInvSecEtanPDFMuTwo5023GeV->GetY();
    }
    Double_t* valueNLOMuHalfPi0nPDF5023GeV              = graphNLOCalcInvSecPi0nPDFMuHalf5023GeV->GetY();
    Double_t* valueNLOMuOnePi0nPDF5023GeV               = graphNLOCalcInvSecPi0nPDFMuOne5023GeV->GetY();
    Double_t* valueNLOMuTwoPi0nPDF5023GeV               = graphNLOCalcInvSecPi0nPDFMuTwo5023GeV->GetY();
    Double_t* xValuenPDFNLO5023GeV                      = graphNLOCalcInvSecPi0nPDFMuOne5023GeV->GetX();
    Int_t xNBinsnPdf5023GeV                             = graphNLOCalcInvSecPi0nPDFMuOne5023GeV->GetN();
    Double_t valueNLOEtaToPi0nPDFNLOMuHalf5023GeV[200];
    Double_t valueNLOEtaToPi0nPDFNLOMuOne5023GeV[200];
    Double_t valueNLOEtaToPi0nPDFNLOMuTwo5023GeV[200];

    // calculate eta/pi0 ratios for different mu scales
    for ( Int_t n = 0; n < xNBinsnPdf5023GeV+1; n++){
        if (valueNLOMuOnePi0nPDF5023GeV[n] != 0){
            valueNLOEtaToPi0nPDFNLOMuOne5023GeV[n] = valueNLOMuOneEtanPDF5023GeV[n]/valueNLOMuOnePi0nPDF5023GeV[n];
        } else {
            valueNLOEtaToPi0nPDFNLOMuOne5023GeV[n] = 0.;
        }

        if (fillAllMuScalesEtanPDF5023GeV){
            if (n == 0){
                valueNLOEtaToPi0nPDFNLOMuHalf5023GeV[n] = 0.;
            } else {
                if (valueNLOMuHalfPi0nPDF5023GeV[n] != 0){
                    valueNLOEtaToPi0nPDFNLOMuHalf5023GeV[n] = valueNLOMuHalfEtanPDF5023GeV[n]/valueNLOMuHalfPi0nPDF5023GeV[n];
                } else {
                    valueNLOEtaToPi0nPDFNLOMuHalf5023GeV[n] = 0.;
                }
            }
            if (valueNLOMuTwoPi0nPDF5023GeV[n] != 0){
                valueNLOEtaToPi0nPDFNLOMuTwo5023GeV[n] = valueNLOMuTwoEtanPDF5023GeV[n]/valueNLOMuTwoPi0nPDF5023GeV[n];
            } else {
                valueNLOEtaToPi0nPDFNLOMuTwo5023GeV[n] = 0.;
            }
        } else {
            valueNLOEtaToPi0nPDFNLOMuHalf5023GeV[n] = valueNLOEtaToPi0nPDFNLOMuOne5023GeV[n];
            valueNLOEtaToPi0nPDFNLOMuTwo5023GeV[n] = valueNLOEtaToPi0nPDFNLOMuOne5023GeV[n];
        }
    }
    // fill graphs
    TGraph* graphEtaToPi0nPDFNLOMuHalf5023GeV           = NULL;
    TGraph* graphEtaToPi0nPDFNLOMuTwo5023GeV            = NULL;
    TGraph* graphEtaToPi0nPDFNLOMuOne5023GeV            = new TGraph(xNBinsnPdf5023GeV,xValuenPDFNLO5023GeV,valueNLOEtaToPi0nPDFNLOMuOne5023GeV);
    if (fillAllMuScalesEta5023GeV){
        graphEtaToPi0nPDFNLOMuHalf5023GeV               = new TGraph(xNBinsnPdf5023GeV,xValuenPDFNLO5023GeV,valueNLOEtaToPi0nPDFNLOMuHalf5023GeV);
        graphEtaToPi0nPDFNLOMuHalf5023GeV->RemovePoint(0);
        graphEtaToPi0nPDFNLOMuTwo5023GeV                = new TGraph(xNBinsnPdf5023GeV,xValuenPDFNLO5023GeV,valueNLOEtaToPi0nPDFNLOMuTwo5023GeV);
    }
    // combine all mu scales into 1 graph for plotting
    TGraphAsymmErrors* graphNLOCalcEtaToPi0nPDF5023GeV  = CombineMuScales(xNBinsnPdf5023GeV, xValuenPDFNLO5023GeV, valueNLOEtaToPi0nPDFNLOMuOne5023GeV, valueNLOEtaToPi0nPDFNLOMuHalf5023GeV, valueNLOEtaToPi0nPDFNLOMuTwo5023GeV);


    // Calculate RpA for 5TeV based on NLOs for pPb and pp 5TeV
    graphNLOCalcInvYieldPi0nPDFpPb5023GeV->Print();
    graphNLOCalcInvYieldPi0pPb5023GeV->Print();
    TGraphAsymmErrors* graphRpANLOPi0nCTEQ              = CalculateGraphRatioToGraphAsym(graphNLOCalcInvYieldPi0nPDFpPb5023GeV,graphNLOCalcInvYieldPi0pPb5023GeV);
    TGraphAsymmErrors* graphRpANLOPi0nCTEQ_diffMuOne    = CalculateGraphRatioToGraphAsym2(graphNLOCalcInvYieldPi0nPDFpPb5023GeV,graphNLOCalcInvYieldPi0MuOnepPb5023GeV);
    TGraph* graphRpANLOPi0nCTEQ_muOne                   = CalculateGraphRatioToGraph(graphNLOCalcInvYieldPi0nPDFMuOnepPb5023GeV,graphNLOCalcInvYieldPi0MuOnepPb5023GeV);
    graphRpANLOPi0nCTEQ_muOne->RemovePoint(0);
    TGraph* graphRpANLOPi0nCTEQ_muHalf                  = CalculateGraphRatioToGraph(graphNLOCalcInvYieldPi0nPDFMuHalfpPb5023GeV,graphNLOCalcInvYieldPi0MuHalfpPb5023GeV);
    TGraph* graphRpANLOPi0nCTEQ_muTwo                   = CalculateGraphRatioToGraph(graphNLOCalcInvYieldPi0nPDFMuTwopPb5023GeV,graphNLOCalcInvYieldPi0MuTwopPb5023GeV);
    graphRpANLOPi0nCTEQ_muTwo->RemovePoint(0);

    TGraphAsymmErrors* graphRpANLOPi0nCTEQSepCalc       = CombineMuScales(graphRpANLOPi0nCTEQ_muOne, graphRpANLOPi0nCTEQ_muHalf, graphRpANLOPi0nCTEQ_muTwo);

    TGraphAsymmErrors* graphRpANLOEtanCTEQ              = CalculateGraphRatioToGraphAsym(graphNLOCalcInvYieldEtanPDFpPb5023GeV,graphNLOCalcInvYieldEtapPb5023GeV);
    TGraphAsymmErrors* graphRpANLOEtanCTEQ_diffMuOne    = CalculateGraphRatioToGraphAsym2(graphNLOCalcInvYieldEtanPDFpPb5023GeV,graphNLOCalcInvYieldEtaMuOnepPb5023GeV);
    TGraph* graphRpANLOEtanCTEQ_muOne                   = CalculateGraphRatioToGraph(graphNLOCalcInvYieldEtanPDFMuOnepPb5023GeV,graphNLOCalcInvYieldEtaMuOnepPb5023GeV);
    graphRpANLOEtanCTEQ_muOne->RemovePoint(0);
    TGraph* graphRpANLOEtanCTEQ_muHalf                  = CalculateGraphRatioToGraph(graphNLOCalcInvYieldEtanPDFMuHalfpPb5023GeV,graphNLOCalcInvYieldEtaMuHalfpPb5023GeV);
    TGraph* graphRpANLOEtanCTEQ_muTwo                   = CalculateGraphRatioToGraph(graphNLOCalcInvYieldEtanPDFMuTwopPb5023GeV,graphNLOCalcInvYieldEtaMuTwopPb5023GeV);
    graphRpANLOEtanCTEQ_muTwo->RemovePoint(0);
    TGraphAsymmErrors* graphRpANLOEtanCTEQSepCalc       = CombineMuScales(graphRpANLOEtanCTEQ_muOne, graphRpANLOEtanCTEQ_muHalf, graphRpANLOEtanCTEQ_muTwo);

    //*********************************************************
    // EPPS16+DSS14 pi0 RpPb
    //*********************************************************
    Int_t nPtPointsRpPbEPPS16_DSS14                     = 0;
    ifstream fileRpPbEPPS16_DSS14;
    TString fileNameRpPbEPPS16_DSS14                    = "ExternalInputpPb/Theory/RpPb5020_pi0_ct14_epps16_dss14.dat";
    fileRpPbEPPS16_DSS14.open(fileNameRpPbEPPS16_DSS14.Data(),ios_base::in);
    cout << "opening: " << fileNameRpPbEPPS16_DSS14.Data() << endl;
    Int_t iPtCurrent    = 0;
    string line;
    vector<Double_t> *valuesEPPS16    = new vector<Double_t>[4];     // iCent x iParticle x nMeasurements matrix with theory curves
    while (getline(fileRpPbEPPS16_DSS14, line) && iPtCurrent < 100) {
        TString temp        = "";
        TString tempBin     = "";
        istringstream cs(line); // controll stream
        cs >> temp;

        if (!temp.Contains("#")){
            istringstream ss(line);
            Int_t iMeasurement  = 0;
            while(ss && iMeasurement < 4){
                ss >> temp;
                if(!temp.IsNull() && temp.CompareTo("nan") != 0){
                    valuesEPPS16[iMeasurement].push_back(temp.Atof());
                    cout << temp.Data() << "\t ";
                    iMeasurement++;
                } else {
                    valuesEPPS16[iMeasurement].push_back(-1.0e-6);
                    cout << -1.0e-6 << "\t ";
                    iMeasurement++;
                }
            }
            cout << endl;
            iPtCurrent++;
        } else {
            cout << "line contains comments" << endl;
        }
    }
    fileRpPbEPPS16_DSS14.close();
    cout << "Number of pT bins: "<< iPtCurrent << "\t number of measurement points: "<<  endl;
    nPtPointsRpPbEPPS16_DSS14                   = iPtCurrent;
    TGraphAsymmErrors* graphRpPbEPPS16_DSS14    = new TGraphAsymmErrors(nPtPointsRpPbEPPS16_DSS14);
    TGraph* graphRpPbEPPS16_DSS14_muOne         = new TGraph(nPtPointsRpPbEPPS16_DSS14);
    for (Int_t iPt = 0; iPt < nPtPointsRpPbEPPS16_DSS14; iPt++){
        graphRpPbEPPS16_DSS14->SetPoint(iPt, valuesEPPS16[0].at(iPt), valuesEPPS16[1].at(iPt));
        graphRpPbEPPS16_DSS14_muOne->SetPoint(iPt, valuesEPPS16[0].at(iPt), valuesEPPS16[1].at(iPt));
        Double_t errUp          = valuesEPPS16[2].at(iPt)-valuesEPPS16[1].at(iPt);
        Double_t errDown        = valuesEPPS16[1].at(iPt)-valuesEPPS16[3].at(iPt);
        graphRpPbEPPS16_DSS14->SetPointError(iPt, 0.01, 0.01, errDown, errUp);
    }
    graphRpPbEPPS16_DSS14->Print();


    //////////////////////////////////////////ColorGlassCondensate predictions for pi0 spectrum ///////////////////////////////
    //Predictions sent by Heikki Mäntysaari [heikki.mantysaari@jyu.fi]
    TString fileNameCGCPi0y0pA5020          = "ExternalInputpPb/Theory/CGC/pA_pi0_y0_sqrts_5020";
    ifstream inCGCPi0;

    inCGCPi0.open(fileNameCGCPi0y0pA5020.Data(),ios_base::in);
    cout<<"*********************CGC spectrum******************"<<fileNameCGCPi0y0pA5020.Data()<<endl;
    Int_t nlinesCGCPi0 = 0;
    Double_t xValue[100];
    Double_t yValue[100];

    string currentline;
    while( getline(inCGCPi0, currentline) && nlinesCGCPi0 < 100 ) {
        TString temp1        = "";
        TString temp2        = "";
        //cout<<currentline<<" hola "<<nlinesCGCPi0<<endl;
        istringstream cs(currentline); // controll stream
        cs >> temp1>>temp2;

        if( !temp1.Contains("#") ) {
            xValue[nlinesCGCPi0] = temp1.Atof();
            yValue[nlinesCGCPi0] = temp2.Atof();
            cout << nlinesCGCPi0 <<"\t"<<xValue[nlinesCGCPi0]<<"\t"<<yValue[nlinesCGCPi0]<<endl;
            nlinesCGCPi0++;
        }
    }

    TGraph* grahCGCPi0y0pA5020              = new TGraph(nlinesCGCPi0,xValue,yValue);

    /*********************************************************************************************************************/
    /* Reading theory from  Ilkka with ct14, epps16, dss14 ***************/
    /* at 5.02 TeV in -0.335 < y_CM < 1.265 (|y_LAB|<0.8).
    /* Calculations in NLO pQCD with CT14+EPPS16+DSS14, based on Nucl.Phys.B883
    /* (2014) 615-628 but with updated (n)PDFs and FFs. From I. Helenius*/

    TString fileNameIlkkapPb5020_pi0[4]     =  {"pPb5020_pi0_ct14_epps16_dss14_scale_err.dat","pPb5020_pi0_ct14_epps16_dss14_err.dat",
                                                "pPb5020_pi0_ct14_epps16_err_dss14.dat","pPb5020_pi0_ct14_errSym_epps16_dss14.dat"};
    Int_t iIlkkaFiles                       = 4;
    Int_t currentIlkkaFile = 0;
    Double_t xValueIlkka[4][100];
    Double_t yValueIlkka[4][100];
    Double_t xErrHighIlkka[4][100];
    Double_t xErrLowIlkka[4][100];
    Double_t yErrHighIlkka[4][100];
    Double_t yErrLowIlkka[4][100];
    TGraphAsymmErrors* graphPi0NLOpQCDct14epps16dss14[4];
    TGraph* graphPi0NLOpQCDct14epps16dss14_muOne            = NULL;
    cout<<"****************************************Ilkka*************************************"<<endl;
    while ( currentIlkkaFile < iIlkkaFiles) {
        Int_t nlinesIlkkaPi0                    = 0;
        TString currentfileNameIlkkaName        = Form("ExternalInputpPb/Theory/%s",fileNameIlkkapPb5020_pi0[currentIlkkaFile].Data());
        ifstream         inIlkkaPi0;

        inIlkkaPi0.open(currentfileNameIlkkaName.Data(),ios_base::in);
        cout<<"*********************Ilkka NLO spectrum******************"<<currentfileNameIlkkaName.Data()<<endl;
        string currentline;
        while( getline(inIlkkaPi0, currentline) && nlinesIlkkaPi0 < 100 ) {
            TString temp1        = "";
            TString temp2        = "";
            TString temp3        = "";
            TString temp4        = "";
            istringstream cs(currentline); // controll stream
            cs >> temp1>>temp2>>temp3>>temp4;
            if( !temp1.Contains("#") ) {
                xValueIlkka[currentIlkkaFile][nlinesIlkkaPi0]   = temp1.Atof();
                yValueIlkka[currentIlkkaFile][nlinesIlkkaPi0]   = temp2.Atof();
                xErrHighIlkka[currentIlkkaFile][nlinesIlkkaPi0] = 0.0;
                xErrLowIlkka[currentIlkkaFile][nlinesIlkkaPi0]  = 0.0;
                yErrHighIlkka[currentIlkkaFile][nlinesIlkkaPi0] = temp3.Atof();
                yErrHighIlkka[currentIlkkaFile][nlinesIlkkaPi0] = yErrHighIlkka[currentIlkkaFile][nlinesIlkkaPi0] - yValueIlkka[currentIlkkaFile][nlinesIlkkaPi0];
                yErrLowIlkka[currentIlkkaFile][nlinesIlkkaPi0]  = temp4.Atof();
                yErrLowIlkka[currentIlkkaFile][nlinesIlkkaPi0]  = yValueIlkka[currentIlkkaFile][nlinesIlkkaPi0] - yErrLowIlkka[currentIlkkaFile][nlinesIlkkaPi0];

                cout << nlinesIlkkaPi0 <<"\t"<<xValueIlkka[currentIlkkaFile][nlinesIlkkaPi0]<<"\t"<<yValueIlkka[currentIlkkaFile][nlinesIlkkaPi0]<<
                "\t"<<yErrHighIlkka[currentIlkkaFile][nlinesIlkkaPi0]<<"\t"<<yErrLowIlkka[currentIlkkaFile][nlinesIlkkaPi0]<<endl;
                nlinesIlkkaPi0++;
            }
        }
        graphPi0NLOpQCDct14epps16dss14[currentIlkkaFile]        = new TGraphAsymmErrors(nlinesIlkkaPi0, xValueIlkka[currentIlkkaFile], yValueIlkka[currentIlkkaFile],
                                                                                        xErrLowIlkka[currentIlkkaFile], xErrHighIlkka[currentIlkkaFile],
                                                                                        yErrLowIlkka[currentIlkkaFile], yErrHighIlkka[currentIlkkaFile]);
        if (currentIlkkaFile == 0){
            graphPi0NLOpQCDct14epps16dss14_muOne                = new TGraph(nlinesIlkkaPi0, xValueIlkka[0], yValueIlkka[0]);
        }
        currentIlkkaFile++;
    }
    //Produce Graph with summ errors

    TGraphAsymmErrors* graphPi0NLOpQCDct14epps16dss14_sumerr    = (TGraphAsymmErrors*)graphPi0NLOpQCDct14epps16dss14[0]->Clone();

    Int_t     nPointsIlkkaSumErr        = graphPi0NLOpQCDct14epps16dss14_sumerr->GetN();
    Double_t *xValueIlkkaSumErr         = graphPi0NLOpQCDct14epps16dss14_sumerr->GetX();
    Double_t *yErrLowIlkkaSumErr        = graphPi0NLOpQCDct14epps16dss14_sumerr->GetEYlow();
    Double_t *yErrHighIlkkaSumErr       = graphPi0NLOpQCDct14epps16dss14_sumerr->GetEYhigh();
    Double_t *yValueIlkkaSumErr         = graphPi0NLOpQCDct14epps16dss14_sumerr->GetY();

    for(Int_t iPoint = 0;  iPoint < nPointsIlkkaSumErr; iPoint++){
        Double_t yELow                  = 0;
        Double_t yEHigh                 = 0;
        for(Int_t nHisto = 0; nHisto < 4; nHisto++) {
            yELow                       += TMath::Power( graphPi0NLOpQCDct14epps16dss14[nHisto]->GetErrorYlow(iPoint)/yValueIlkkaSumErr[iPoint]  , 2 );
            yEHigh                      +=  TMath::Power( graphPi0NLOpQCDct14epps16dss14[nHisto]->GetErrorYhigh(iPoint)/yValueIlkkaSumErr[iPoint] , 2 );
        }
        yELow                           = TMath::Sqrt( yELow  );
        yEHigh                          = TMath::Sqrt( yEHigh );
        yErrLowIlkkaSumErr[iPoint]      = yELow*yValueIlkkaSumErr[iPoint];
        yErrHighIlkkaSumErr[iPoint]     = yEHigh*yValueIlkkaSumErr[iPoint];
    }

    Double_t v0ANDxSection = 2.09; //barns
    Double_t v0ANDxSectionNSD = v0ANDxSection / 0.964;
    Double_t factorHelenius = 208;

    for(Int_t nGraph = 0; nGraph < 4; nGraph++) {
        graphPi0NLOpQCDct14epps16dss14[nGraph] = (TGraphAsymmErrors*)ScaleGraphAsym(graphPi0NLOpQCDct14epps16dss14[nGraph],factorHelenius/(v0ANDxSectionNSD*recalcBarn));
    }
    graphPi0NLOpQCDct14epps16dss14_muOne        = (TGraph*)ScaleGraph(graphPi0NLOpQCDct14epps16dss14_muOne,factorHelenius/(v0ANDxSectionNSD*recalcBarn));
    graphPi0NLOpQCDct14epps16dss14_sumerr       = (TGraphAsymmErrors*)ScaleGraphAsym(graphPi0NLOpQCDct14epps16dss14_sumerr,factorHelenius/(v0ANDxSectionNSD*recalcBarn));

    // **********************************************************************************************************************
    // Load cocktail for pPb pure mt-scaling
    // **********************************************************************************************************************
    TString nameCocktailFileMtScaling       = "CocktailInput/GammaCocktail_pPb_MB_pureMT.root";
    TFile* fileCockailPureMtScaling         = new TFile(nameCocktailFileMtScaling.Data());
    TH1D* histoPi0PureMtScaling             = (TH1D*)fileCockailPureMtScaling->Get("Pi0_Pt_OrBin");
    TH1D* histoEtaPureMtScaling             = (TH1D*)fileCockailPureMtScaling->Get("Eta_Pt_OrBin");
    for (Int_t i = 1; i < histoPi0PureMtScaling->GetNbinsX()+1; i++){
        histoPi0PureMtScaling->SetBinContent(i, histoPi0PureMtScaling->GetBinContent(i)/(histoPi0PureMtScaling->GetBinCenter(i) *2* TMath::Pi()) );
        histoPi0PureMtScaling->SetBinError(i, histoPi0PureMtScaling->GetBinError(i)/(histoPi0PureMtScaling->GetBinCenter(i)*2* TMath::Pi()));
    }
    histoPi0PureMtScaling->GetXaxis()->SetRangeUser(0,50);
    for (Int_t i = 1; i < histoEtaPureMtScaling->GetNbinsX()+1; i++){
        histoEtaPureMtScaling->SetBinContent(i, histoEtaPureMtScaling->GetBinContent(i)/(histoEtaPureMtScaling->GetBinCenter(i) *2* TMath::Pi()) );
        histoEtaPureMtScaling->SetBinError(i, histoEtaPureMtScaling->GetBinError(i)/(histoEtaPureMtScaling->GetBinCenter(i)*2* TMath::Pi()));
    }
    histoEtaPureMtScaling->GetXaxis()->SetRangeUser(0,50);
    TH1D* histoEtaPi0PureMtScaling          = (TH1D*)histoEtaPureMtScaling->Clone("histoEtaPi0PureMtScaling");
    histoEtaPi0PureMtScaling->Divide(histoEtaPi0PureMtScaling,histoPi0PureMtScaling);

    //**************************************************************************************************
    //********************** EPOS+JJ MC spectra and ratio ***********************************************
    //**************************************************************************************************
    // file generated with TaskV1/ExtractMCInputSpectraFromFile.C++ based on PCM-EMC inputs
    TFile* fileEPOSJJpPb8TeV               = new TFile("ExternalInputpPb/Theory/MCInputCompilationLHC18b9bc_pPb8TeV_10.root");
    TH1D* histoPi0EPOSJJpPb8TeV            = (TH1D*)fileEPOSJJpPb8TeV->Get("pPb_8TeV/MC_Pi0_Pt");
    TH1D* histoPi0EPOSJJpPb8TeVReb         = (TH1D*)fileEPOSJJpPb8TeV->Get("pPb_8TeV/MC_Pi0_Pt_Rebinned");
    TH1D* histoPiChEPOSJJpPb8TeV           = (TH1D*)fileEPOSJJpPb8TeV->Get("pPb_8TeV/MC_PiCh_All_Pt");
    TH1D* histoKChEPOSJJpPb8TeV            = (TH1D*)fileEPOSJJpPb8TeV->Get("pPb_8TeV/MC_KCh_All_Pt");
    TH1D* histoEtaEPOSJJpPb8TeV            = (TH1D*)fileEPOSJJpPb8TeV->Get("pPb_8TeV/MC_Eta_Pt");
    TH1D* histoEtaEPOSJJpPb8TeVReb         = (TH1D*)fileEPOSJJpPb8TeV->Get("pPb_8TeV/MC_Eta_Pt_Rebinned");
    TH1D* histoEtaToPi0EPOSJJpPb8TeV       = (TH1D*)fileEPOSJJpPb8TeV->Get("pPb_8TeV/MCEtaToPi0");
    TH1D* histoEtaToKCHEPOSJJpPb8TeV       = (TH1D*)fileEPOSJJpPb8TeV->Get("pPb_8TeV/MCEtaToKCh");
    TH1D* histoPi0ToPiCHEPOSJJpPb8TeV      = (TH1D*)fileEPOSJJpPb8TeV->Get("pPb_8TeV/MCPi0ToPiCh");

    TFile* fileDPMJETpPb8TeV               = new TFile("ExternalInputpPb/Theory/MCInputCompilationLHC18f3bc_pPb8TeV_4.root");
    TH1D* histoPi0DPMJETpPb8TeV            = (TH1D*)fileDPMJETpPb8TeV->Get("pPb_8TeV/MC_Pi0_Pt");
    TH1D* histoPi0DPMJETpPb8TeVReb         = (TH1D*)fileDPMJETpPb8TeV->Get("pPb_8TeV/MC_Pi0_Pt_Rebinned");
    TH1D* histoPiChDPMJETpPb8TeV           = (TH1D*)fileDPMJETpPb8TeV->Get("pPb_8TeV/MC_PiCh_All_Pt");
    TH1D* histoKChDPMJETpPb8TeV            = (TH1D*)fileDPMJETpPb8TeV->Get("pPb_8TeV/MC_KCh_All_Pt");
    TH1D* histoEtaDPMJETpPb8TeV            = (TH1D*)fileDPMJETpPb8TeV->Get("pPb_8TeV/MC_Eta_Pt");
    TH1D* histoEtaDPMJETpPb8TeVReb         = (TH1D*)fileDPMJETpPb8TeV->Get("pPb_8TeV/MC_Eta_Pt_Rebinned");
    TH1D* histoEtaToPi0DPMJETpPb8TeV       = (TH1D*)fileDPMJETpPb8TeV->Get("pPb_8TeV/MCEtaToPi0");
    TH1D* histoEtaToKCHDPMJETpPb8TeV       = (TH1D*)fileDPMJETpPb8TeV->Get("pPb_8TeV/MCEtaToKCh");
    TH1D* histoPi0ToPiCHDPMJETpPb8TeV      = (TH1D*)fileDPMJETpPb8TeV->Get("pPb_8TeV/MCPi0ToPiCh");

    //**************************************************************************************************
    //********************** Pythia8 EPPS16 spectra ***********************************************
    //**************************************************************************************************
    TFile* filePythia8EPPS16               = new TFile("ExternalInputpPb/Theory/pythia8_8160GeV_pPb_EPPS16.root");
    TH1D* histoPi0Pythia8EPPS16pPb8TeV            = (TH1D*)filePythia8EPPS16->Get("h_invXsec_pi0primary_etaEMCal");
    histoPi0Pythia8EPPS16pPb8TeV->Scale(1/(2.1*1e12));
    TH1D* histoEtaPythia8EPPS16pPb8TeV            = (TH1D*)filePythia8EPPS16->Get("h_invXsec_eta_etaEMCal");
    histoEtaPythia8EPPS16pPb8TeV->Scale(1/(2.1*1e12));

    //**************************************************************************************************
    //********************** Pythia8 EPPS16 spectra ***********************************************
    //**************************************************************************************************
    TFile* filePythia8nCTEQ15               = new TFile("ExternalInputpPb/Theory/pythia8_8160GeV_pPb_nCTEQ15.root");
    TH1D* histoPi0Pythia8nCTEQ15pPb8TeV            = (TH1D*)filePythia8nCTEQ15->Get("h_invXsec_pi0primary_etaEMCal");
    histoPi0Pythia8nCTEQ15pPb8TeV->Scale(1/(2.1*1e12));
    TH1D* histoEtaPythia8nCTEQ15pPb8TeV            = (TH1D*)filePythia8nCTEQ15->Get("h_invXsec_eta_etaEMCal");
    histoEtaPythia8nCTEQ15pPb8TeV->Scale(1/(2.1*1e12));

    //**************************************************************************************************
    //********************** Hijing spectra ***********************************************
    //**************************************************************************************************
    TFile* fileHIJING8160GeV              = new TFile("ExternalInputpPb/Theory/Hijing_pPb_8160GeV_1034Mio.root");
    TH1D* histoPi0HIJING8160GeV               = (TH1D*)fileHIJING8160GeV->Get("hPtPi0HijingpPb");
    TH1D* histoPiPlHIJING8160GeV              = (TH1D*)fileHIJING8160GeV->Get("hPtPiPlHijingpPb");
    TH1D* histoPiMiHIJING8160GeV              = (TH1D*)fileHIJING8160GeV->Get("hPtPiMiHijingpPb");
    TH1D* histoPiChHIJING8160GeV              = (TH1D*)histoPiPlHIJING8160GeV->Clone("hPtPiChHijingpPb");
    histoPiChHIJING8160GeV->Add(histoPiMiHIJING8160GeV);
    histoPiChHIJING8160GeV->Scale(1/2.);
    TH1D* histoEtaHIJING8160GeV               = (TH1D*)fileHIJING8160GeV->Get("hPtEtaHijingpPb");
    TH1D* histoEtaToPi0HIJING8160GeV          = (TH1D*)histoEtaHIJING8160GeV->Clone("MCEtaToPi0Hijing");
    histoEtaToPi0HIJING8160GeV->Divide(histoPi0HIJING8160GeV);
    TH1D* histoPi0ToPiCHHIJING8160GeV         = (TH1D*)histoPi0HIJING8160GeV->Clone("MCPi0ToPiChHijing");
    histoPi0ToPiCHHIJING8160GeV->Divide(histoPiChHIJING8160GeV);

    //**********************************************************************************************************************
    //******************************** Electron cross section from weak Bosons *********************************************
    //**********************************************************************************************************************
    // using POWHEG processes for single boson production ("POWHEG-BOX/W", "POWHEG-BOX/Z")
    // using NNPDF2.3NLO as PDF, the NLO version of the PDF used in Pythia8 Monash 2013 tune
    TFile *filepPb8TeVPowhegElecFromWeakBosons  = new TFile("ExternalInput/Theory/Powheg/powheg_electronsFromWeakBoson_8TeV_pdf244600_graphs.root");
    TGraphAsymmErrors* graphElecFromWeakBosonPowhegXSecpPb8TeVEMCal = (TGraphAsymmErrors*)filepPb8TeVPowhegElecFromWeakBosons->Get("gae_XSec_Elec_all");
    // using Pythia8 JJ
    TFile *filepPb8TeVPythiaElecFromJJ = new TFile("ExternalInput/Theory/Pythia/pythia8_8TeV_compilation_poppenborg_onlyJJ.root");
    TGraphAsymmErrors* graphElecFromJJXSecpPb8TeVEMCal = (TGraphAsymmErrors*)filepPb8TeVPythiaElecFromJJ->Get("gae_electrons_etaEMCal");

    // fitting ratio (electrons from weak bosons)/(electrons from JJ MC)
    graphElecFromWeakBosonPowhegXSecpPb8TeVEMCal->RemovePoint(0);
    graphElecFromJJXSecpPb8TeVEMCal->RemovePoint(0);
    TGraphAsymmErrors* graphRatioElecFromWeakBoson   = CalculateAsymGraphRatioToGraph(graphElecFromWeakBosonPowhegXSecpPb8TeVEMCal, graphElecFromJJXSecpPb8TeVEMCal);
    TSpline3* splineRatioElecFromWeakBoson = new TSpline3("splineRatioElecFromWeakBoson_pPb8TeV", graphRatioElecFromWeakBoson,"b2e2",0,0);

    //**********************************************************************************************************************
    //********************************* Write graphs and histos to compilation file for pPb ********************************
    //**********************************************************************************************************************
    TFile fileTheoryGraphsPPb("ExternalInputpPb/Theory/TheoryCompilationPPb.root","UPDATE");

        // create subdirectories in file if needed
        TDirectoryFile* directory5TeV = (TDirectoryFile*)fileTheoryGraphsPPb.Get("pPb_5.023TeV");
        if (!directory5TeV){
            fileTheoryGraphsPPb.mkdir("pPb_5.023TeV");
        }
        TDirectoryFile* directory5TeV0005 = (TDirectoryFile*)fileTheoryGraphsPPb.Get("0-5_pPb_5.023TeV");
        if (!directory5TeV0005){
            fileTheoryGraphsPPb.mkdir("0-5_pPb_5.023TeV");
        }
        TDirectoryFile* directory5TeV0510 = (TDirectoryFile*)fileTheoryGraphsPPb.Get("5-10_pPb_5.023TeV");
        if (!directory5TeV0510){
            fileTheoryGraphsPPb.mkdir("5-10_pPb_5.023TeV");
        }

        TDirectoryFile* directory5TeVCents[5]   = {NULL, NULL, NULL, NULL, NULL};
        for (Int_t cent = 1; cent < 5; cent++){
            TString currentName = Form("%spPb_5.023TeV", centrality[cent].Data());
            directory5TeVCents[cent] = (TDirectoryFile*)fileTheoryGraphsPPb.Get(currentName.Data());
            if (!directory5TeVCents[cent]){
                fileTheoryGraphsPPb.mkdir(currentName.Data());
            }
            fileTheoryGraphsPPb.cd(currentName.Data());
            histoPi0DPMJet[cent]->Write("histoPi0SpecDPMJet", TObject::kOverwrite);
            histoEtaDPMJet[cent]->Write("histoEtaSpecDPMJet", TObject::kOverwrite);
            histoPiChDPMJet[cent]->Write("histoPiChSpecDPMJet", TObject::kOverwrite);
            histoKChDPMJet[cent]->Write("histoKChSpecDPMJet", TObject::kOverwrite);
            histoPi0DPMJetReb[cent]->Write("histoPi0SpecDPMJet_Reb", TObject::kOverwrite);
            histoEtaDPMJetReb[cent]->Write("histoEtaSpecDPMJet_Reb", TObject::kOverwrite);
            histoEtaToPi0DPMJet[cent]->Write("histoEtaToPi0DPMJet", TObject::kOverwrite);
            histoEtaToKCHDPMJet[cent]->Write("histoEtaToKChDPMJet", TObject::kOverwrite);
            histoPi0ToPiCHDPMJet[cent]->Write("histoPi0ToPiChDPMJet", TObject::kOverwrite);

            histoPi0EPOSLHC[cent]->Write("histoPi0SpecEPOSLHC", TObject::kOverwrite);
            histoEtaEPOSLHC[cent]->Write("histoEtaSpecEPOSLHC", TObject::kOverwrite);
            histoPiChEPOSLHC[cent]->Write("histoPiChSpecEPOSLHC", TObject::kOverwrite);
            histoKChEPOSLHC[cent]->Write("histoKChSpecEPOSLHC", TObject::kOverwrite);
            histoPi0EPOSLHCReb[cent]->Write("histoPi0SpecEPOSLHC_Reb", TObject::kOverwrite);
            histoEtaEPOSLHCReb[cent]->Write("histoEtaSpecEPOSLHC_Reb", TObject::kOverwrite);
            histoEtaToPi0EPOSLHC[cent]->Write("histoEtaToPi0EPOSLHC", TObject::kOverwrite);
            histoEtaToKCHEPOSLHC[cent]->Write("histoEtaToKChEPOSLHC", TObject::kOverwrite);
            histoPi0ToPiCHEPOSLHC[cent]->Write("histoPi0ToPiChEPOSLHC", TObject::kOverwrite);
        }

        // write MB calcs
        fileTheoryGraphsPPb.cd("pPb_5.023TeV");
            // pi0 EPS09 with old FFs RpPb
            graphPi0RpAEPS09sKKP->Write("graphPi0RpAEPS09sKKP5023GeV", TObject::kOverwrite);
            graphPi0RpAEPS09sAKK->Write("graphPi0RpAEPS09sAKK5023GeV", TObject::kOverwrite);
            graphPi0RpAEPS09sDSS->Write("graphPi0RpAEPS09sDSS5023GeV", TObject::kOverwrite);
            graphPi0RpAAsymmErrEPS09sDSS->Write("graphPi0RpAAsymmErrEPS09sDSS5023GeV", TObject::kOverwrite);
            // pi0, eta, eta/pi0 EP03
            histoPi0EPOS3->Write("histoPi0SpecEPOS35023GeV", TObject::kOverwrite);
            histoEtaEPOS3->Write("histoEtaSpecEPOS35023GeV", TObject::kOverwrite);
            histoPi0EPOS3Reb->Write("histoPi0SpecEPOS35023GeV_Reb", TObject::kOverwrite);
            histoEtaEPOS3Reb->Write("histoEtaSpecEPOS35023GeV_Reb", TObject::kOverwrite);
            histoEtaToPi0EPOS3->Write("histoEtaToPi0EPOS35023GeV", TObject::kOverwrite);
            histoEtaToPi0EPOS3Reb->Write("histoEtaToPi0EPOS35023GeV_Reb", TObject::kOverwrite);
            // pi0, eta, eta/pi0 EPOSLHC
            histoPi0EPOSLHC[0]->Write("histoPi0SpecEPOSLHC", TObject::kOverwrite);
            histoEtaEPOSLHC[0]->Write("histoEtaSpecEPOSLHC", TObject::kOverwrite);
            histoPiChEPOSLHC[0]->Write("histoPiChSpecEPOSLHC", TObject::kOverwrite);
            histoKChEPOSLHC[0]->Write("histoKChSpecEPOSLHC", TObject::kOverwrite);
            histoPi0EPOSLHCReb[0]->Write("histoPi0SpecEPOSLHC_Reb", TObject::kOverwrite);
            histoEtaEPOSLHCReb[0]->Write("histoEtaSpecEPOSLHC_Reb", TObject::kOverwrite);
            histoEtaToPi0EPOSLHC[0]->Write("histoEtaToPi0EPOSLHC", TObject::kOverwrite);
            histoEtaToKCHEPOSLHC[0]->Write("histoEtaToKChEPOSLHC", TObject::kOverwrite);
            histoPi0ToPiCHEPOSLHC[0]->Write("histoPi0ToPiChEPOSLHC", TObject::kOverwrite);
            histoPi0EPOSLHC2->Write("histoPi0SpecEPOSLHC_MCGen", TObject::kOverwrite);
            histoPiChEPOSLHC2->Write("histoPiChSpecEPOSLHC_MCGen", TObject::kOverwrite);
            histoEtaEPOSLHC2->Write("histoEtaSpecEPOSLHC_MCGen", TObject::kOverwrite);
            histoEtaToPi0EPOSLHC2->Write("histoEtaToPi0EPOSLHC_MCGen", TObject::kOverwrite);
            histoPi0ToPiCHEPOSLHC2->Write("histoPi0ToPiChEPOSLHC_MCGen", TObject::kOverwrite);

            // pi0, eta, eta/pi0 DPMJet
            histoPi0DPMJet[0]->Write("histoPi0SpecDPMJet", TObject::kOverwrite);
            histoEtaDPMJet[0]->Write("histoEtaSpecDPMJet", TObject::kOverwrite);
            histoPiChDPMJet[0]->Write("histoPiChSpecDPMJet", TObject::kOverwrite);
            histoKChDPMJet[0]->Write("histoKChSpecDPMJet", TObject::kOverwrite);
            histoPi0DPMJetReb[0]->Write("histoPi0SpecDPMJet_Reb", TObject::kOverwrite);
            histoEtaDPMJetReb[0]->Write("histoEtaSpecDPMJet_Reb", TObject::kOverwrite);
            histoEtaToPi0DPMJet[0]->Write("histoEtaToPi0DPMJet", TObject::kOverwrite);
            histoEtaToKCHDPMJet[0]->Write("histoEtaToKChDPMJet", TObject::kOverwrite);
            histoPi0ToPiCHDPMJet[0]->Write("histoPi0ToPiChDPMJet", TObject::kOverwrite);
            histoPi0DPMJet2->Write("histoPi0SpecDPMJet_MCGen", TObject::kOverwrite);
            histoPiChDPMJet2->Write("histoPiChSpecDPMJet_MCGen", TObject::kOverwrite);
            histoEtaDPMJet2->Write("histoEtaSpecDPMJet_MCGen", TObject::kOverwrite);
            histoEtaToPi0DPMJet2->Write("histoEtaToPi0DPMJet_MCGen", TObject::kOverwrite);
            histoPi0ToPiCHDPMJet2->Write("histoPi0ToPiChDPMJet_MCGen", TObject::kOverwrite);

            // pi0, eta, eta/pi0 HIJING
            histoPi0HIJING->Write("histoPi0SpecHIJING", TObject::kOverwrite);
            histoEtaHIJING->Write("histoEtaSpecHIJING", TObject::kOverwrite);
            histoPiChHIJING->Write("histoPiChSpecHIJING", TObject::kOverwrite);
            histoKChHIJING->Write("histoKChSpecHIJING", TObject::kOverwrite);
            histoPi0HIJINGReb->Write("histoPi0SpecHIJING_Reb", TObject::kOverwrite);
            histoEtaHIJINGReb->Write("histoEtaSpecHIJING_Reb", TObject::kOverwrite);
            histoEtaToPi0HIJING->Write("histoEtaToPi0HIJING", TObject::kOverwrite);
            histoEtaToKCHHIJING->Write("histoEtaToKChHIJING", TObject::kOverwrite);
            histoPi0ToPiCHHIJING->Write("histoPi0ToPiChHIJING", TObject::kOverwrite);
            histoPi0HIJING2->Write("histoPi0SpecHIJING_MCGen", TObject::kOverwrite);
            histoPiChHIJING2->Write("histoPiChSpecHIJING_MCGen", TObject::kOverwrite);
            histoEtaHIJING2->Write("histoEtaSpecHIJING_MCGen", TObject::kOverwrite);
            histoEtaToPi0HIJING2->Write("histoEtaToPi0HIJING_MCGen", TObject::kOverwrite);
            histoPi0ToPiCHHIJING2->Write("histoPi0ToPiChHIJING_MCGen", TObject::kOverwrite);
            // pi0, eta, eta/pi0 McGill Hydro
            graphEtaSpecMcGill->Write("graphEtaSpecMcGill5023GeV", TObject::kOverwrite);
            graphPi0SpecMcGill->Write("graphPi0SpecMcGill5023GeV", TObject::kOverwrite);
            graphEtaV2McGill->Write("graphEtaV2McGill5023GeV", TObject::kOverwrite);
            graphPi0V2McGill->Write("graphPi0V2McGill5023GeV", TObject::kOverwrite);
            graphEtaToPi0McGill->Write("graphEtaToPi0McGill5023GeV", TObject::kOverwrite);

            // pi0 Vogelsang PDF: CT10, FF: DSS14 scaled to pPb
            graphNLOCalcInvYieldPi0MuHalfpPb5023GeV->Write("graphNLOCalcDSS14InvYieldPi0MuHalf5023GeV", TObject::kOverwrite);
            graphNLOCalcInvYieldPi0MuOnepPb5023GeV->Write("graphNLOCalcDSS14InvYieldPi0MuOne5023GeV", TObject::kOverwrite);
            graphNLOCalcInvYieldPi0MuTwopPb5023GeV->Write("graphNLOCalcDSS14InvYieldPi0MuTwo5023GeV", TObject::kOverwrite);
            graphNLOCalcInvYieldPi0pPb5023GeV->Write("graphNLOCalcDSS14InvYieldPi05023GeV", TObject::kOverwrite);
            graphNLOCalcInvYieldPi0MuHalfpp5023GeV->Write("graphNLOCalcDSS14InvYieldPi0MuHalfPP5023GeV", TObject::kOverwrite);
            graphNLOCalcInvYieldPi0MuOnepp5023GeV->Write("graphNLOCalcDSS14InvYieldPi0MuOnePP5023GeV", TObject::kOverwrite);
            graphNLOCalcInvYieldPi0MuTwopp5023GeV->Write("graphNLOCalcDSS14InvYieldPi0MuTwoPP5023GeV", TObject::kOverwrite);
            graphNLOCalcInvYieldPi0pp5023GeV->Write("graphNLOCalcDSS14InvYieldPi0PP5023GeV", TObject::kOverwrite);
            // eta Vogelsang PDF: CT10, FF AESSS
            if (graphNLOCalcInvYieldEtaMuHalfpPb5023GeV) graphNLOCalcInvYieldEtaMuHalfpPb5023GeV->Write("graphNLOCalcAESSSInvYieldEtaMuHalf5023GeV", TObject::kOverwrite);
            if (graphNLOCalcInvYieldEtaMuOnepPb5023GeV) graphNLOCalcInvYieldEtaMuOnepPb5023GeV->Write("graphNLOCalcAESSSInvYieldEtaMuOne5023GeV", TObject::kOverwrite);
            if (graphNLOCalcInvYieldEtaMuTwopPb5023GeV) graphNLOCalcInvYieldEtaMuTwopPb5023GeV->Write("graphNLOCalcAESSSInvYieldEtaMuTwo5023GeV", TObject::kOverwrite);
            graphNLOCalcInvYieldEtapPb5023GeV->Write("graphNLOCalcAESSSInvYieldEta5023GeV", TObject::kOverwrite);
            if (graphNLOCalcInvYieldEtaMuHalfpp5023GeV) graphNLOCalcInvYieldEtaMuHalfpp5023GeV->Write("graphNLOCalcAESSSInvYieldEtaMuHalfPP5023GeV", TObject::kOverwrite);
            if (graphNLOCalcInvYieldEtaMuOnepp5023GeV) graphNLOCalcInvYieldEtaMuOnepp5023GeV->Write("graphNLOCalcAESSSInvYieldEtaMuOnePP5023GeV", TObject::kOverwrite);
            if (graphNLOCalcInvYieldEtaMuTwopp5023GeV) graphNLOCalcInvYieldEtaMuTwopp5023GeV->Write("graphNLOCalcAESSSInvYieldEtaMuTwoPP5023GeV", TObject::kOverwrite);
            graphNLOCalcInvYieldEtapp5023GeV->Write("graphNLOCalcAESSSInvYieldEtaPP5023GeV", TObject::kOverwrite);
            // eta/pi0 Vogelsang PDF: CT10, FF eta- AESSS, FF pi0 - DSS14
            if (graphEtaToPi0NLOMuHalf5023GeV) graphEtaToPi0NLOMuHalf5023GeV->Write("graphNLOCalcEtaOverPi0MuHalf5023GeV", TObject::kOverwrite);
            if (graphEtaToPi0NLOMuOne5023GeV) graphEtaToPi0NLOMuOne5023GeV->Write("graphNLOCalcEtaOverPi0MuOne5023GeV", TObject::kOverwrite);
            if (graphEtaToPi0NLOMuTwo5023GeV) graphEtaToPi0NLOMuTwo5023GeV->Write("graphNLOCalcEtaOverPi0MuTwo5023GeV", TObject::kOverwrite);
            graphNLOCalcEtaToPi05023GeV->Write("graphNLOCalcEtaOverPi05023GeV_AESSS_DSS14", TObject::kOverwrite);

            // pi0 Vogelsang nPDF: nCTEQ, FF: DSS14 scaled to pPb
            if (graphNLOCalcInvYieldPi0nPDFMuHalfpPb5023GeV) graphNLOCalcInvYieldPi0nPDFMuHalfpPb5023GeV->Write("graphNLOCalcDSS14InvYieldPi0MuHalf5023GeV_nCTEQ", TObject::kOverwrite);
            if (graphNLOCalcInvYieldPi0nPDFMuOnepPb5023GeV) graphNLOCalcInvYieldPi0nPDFMuOnepPb5023GeV->Write("graphNLOCalcDSS14InvYieldPi0MuOne5023GeV_nCTEQ", TObject::kOverwrite);
            if (graphNLOCalcInvYieldPi0nPDFMuTwopPb5023GeV) graphNLOCalcInvYieldPi0nPDFMuTwopPb5023GeV->Write("graphNLOCalcDSS14InvYieldPi0MuTwo5023GeV_nCTEQ", TObject::kOverwrite);
            graphNLOCalcInvYieldPi0nPDFpPb5023GeV->Write("graphNLOCalcDSS14InvYieldPi05023GeV_nCTEQ", TObject::kOverwrite);
            // eta Vogelsang nPDF: nCTEQ, FF AESSS
            if (graphNLOCalcInvYieldEtanPDFMuHalfpPb5023GeV) graphNLOCalcInvYieldEtanPDFMuHalfpPb5023GeV->Write("graphNLOCalcAESSSInvYieldEtaMuHalf5023GeV_nCTEQ", TObject::kOverwrite);
            if (graphNLOCalcInvYieldEtanPDFMuOnepPb5023GeV) graphNLOCalcInvYieldEtanPDFMuOnepPb5023GeV->Write("graphNLOCalcAESSSInvYieldEtaMuOne5023GeV_nCTEQ", TObject::kOverwrite);
            if (graphNLOCalcInvYieldEtanPDFMuTwopPb5023GeV) graphNLOCalcInvYieldEtanPDFMuTwopPb5023GeV->Write("graphNLOCalcAESSSInvYieldEtaMuTwo5023GeV_nCTEQ", TObject::kOverwrite);
            graphNLOCalcInvYieldEtanPDFpPb5023GeV->Write("graphNLOCalcAESSSInvYieldEta5023GeV_nCTEQ", TObject::kOverwrite);
            // eta/pi0 Vogelsang nPDF: nCTEQ, FF eta- AESSS, FF pi0 - DSS14
            if (graphEtaToPi0nPDFNLOMuHalf5023GeV) graphEtaToPi0nPDFNLOMuHalf5023GeV->Write("graphNLOCalcEtaOverPi0MuHalf5023GeV_nCTEQ", TObject::kOverwrite);
            if (graphEtaToPi0nPDFNLOMuOne5023GeV) graphEtaToPi0nPDFNLOMuOne5023GeV->Write("graphNLOCalcEtaOverPi0MuOne5023GeV_nCTEQ", TObject::kOverwrite);
            if (graphEtaToPi0nPDFNLOMuTwo5023GeV) graphEtaToPi0nPDFNLOMuTwo5023GeV->Write("graphNLOCalcEtaOverPi0MuTwo5023GeV_nCTEQ", TObject::kOverwrite);
            graphNLOCalcEtaToPi0nPDF5023GeV->Write("graphNLOCalcEtaOverPi05023GeV_AESSS_DSS14_nCTEQ", TObject::kOverwrite);

            graphRpANLOPi0nCTEQ->Write("graphNLOCalcDSS14RpAPi05023GeV_nCTEQ", TObject::kOverwrite);
            graphRpANLOPi0nCTEQ_diffMuOne->Write("graphNLOCalcDSS14RpAPi05023GeV_nCTEQ_onlypPbErrs", TObject::kOverwrite);
            graphRpANLOPi0nCTEQ_muOne->Write("graphNLOCalcDSS14RpAPi05023GeV_muOne_nCTEQ", TObject::kOverwrite);
            graphRpANLOPi0nCTEQ_muHalf->Write("graphNLOCalcDSS14RpAPi05023GeV_muHalf_nCTEQ", TObject::kOverwrite);
            graphRpANLOPi0nCTEQ_muTwo->Write("graphNLOCalcDSS14RpAPi05023GeV_muTwo_nCTEQ", TObject::kOverwrite);
            graphRpANLOPi0nCTEQSepCalc->Write("graphNLOCalcDSS14RpAPi05023GeV_nCTEQ_SepCalc", TObject::kOverwrite);

            graphRpANLOEtanCTEQ->Write("graphNLOCalcAESSSRpAEta5023GeV_nCTEQ", TObject::kOverwrite);
            graphRpANLOEtanCTEQ_diffMuOne->Write("graphNLOCalcAESSSRpAEta5023GeV_nCTEQ_onlypPbErrs", TObject::kOverwrite);
            graphRpANLOEtanCTEQ_muOne->Write("graphNLOCalcAESSSRpAEta5023GeV_muOne_nCTEQ", TObject::kOverwrite);
            graphRpANLOEtanCTEQ_muHalf->Write("graphNLOCalcAESSSRpAEta5023GeV_muHalf_nCTEQ", TObject::kOverwrite);
            graphRpANLOEtanCTEQ_muTwo->Write("graphNLOCalcAESSSRpAEta5023GeV_muTwo_nCTEQ", TObject::kOverwrite);
            graphRpANLOEtanCTEQSepCalc->Write("graphNLOCalcAESSSRpAEta5023GeV_nCTEQ_SepCalc", TObject::kOverwrite);

            graphRpPbEPPS16_DSS14->Write("graphNLOCalcDSS14RpAPi05023GeV_EPPS16", TObject::kOverwrite);
            graphRpPbEPPS16_DSS14_muOne->Write("graphNLOCalcDSS14RpAPi05023GeV_muOne_EPPS16", TObject::kOverwrite);

            // pi0 CGC RpPb
            graphPi0RpACGC->Write("graphPi0RpACGC5023GeV", TObject::kOverwrite);
            //pi0 CGC
            grahCGCPi0y0pA5020->Write("graphPi0SpecCGC5023GeV", TObject::kOverwrite);
            grahCGCPi0y0pA5020->Print();
            //Pi0 NLO pQCD
            if(graphPi0NLOpQCDct14epps16dss14_sumerr) graphPi0NLOpQCDct14epps16dss14_sumerr->Write("graphNLOpQCDPi0_ct14_epps16_dss14_sumerr", TObject::kOverwrite);
            if(graphPi0NLOpQCDct14epps16dss14_muOne) graphPi0NLOpQCDct14epps16dss14_muOne->Write("graphNLOpQCDPi0_ct14_epps16_dss14_muOne", TObject::kOverwrite);
            if( graphPi0NLOpQCDct14epps16dss14[0] ) graphPi0NLOpQCDct14epps16dss14[0]->Write("graphNLOpQCDPi0_ct14_epps16_dss14_scale_err", TObject::kOverwrite);
            if( graphPi0NLOpQCDct14epps16dss14[1] ) graphPi0NLOpQCDct14epps16dss14[1]->Write("graphNLOpQCDPi0_ct14_epps16_dss14_err", TObject::kOverwrite);
            if( graphPi0NLOpQCDct14epps16dss14[2] ) graphPi0NLOpQCDct14epps16dss14[2]->Write("graphNLOpQCDPi0_ct14_epps16_err_dss14", TObject::kOverwrite);
            if( graphPi0NLOpQCDct14epps16dss14[3] ) graphPi0NLOpQCDct14epps16dss14[3]->Write("graphNLOpQCDPi0_ct14_errSym_epps16_dss14", TObject::kOverwrite);

            graphPi0NLOpQCDct14epps16dss14_muOne->Print();

            histoEtaPi0PureMtScaling->GetYaxis()->SetTitle("#eta/#pi^{0}");
            histoEtaPi0PureMtScaling->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            histoEtaPi0PureMtScaling->Write("histoEtaPi0PureMtScaling_ALICECombPi0", TObject::kOverwrite);

        // write McGill calc for different cents and particles
        for (Int_t iCent = 0; iCent < nCent; iCent++){
            // go to correct directory
            fileTheoryGraphsPPb.cd(Form("%spPb_5.023TeV",centNamesOut[iCent].Data()));
            // particle loop
            for(Int_t iParticle=0; iParticle<nParticles; iParticle++){
                // write spectra
                if (graphSpectraMcGill[iCent][iParticle])
                    graphSpectraMcGill[iCent][iParticle]->Write(Form("graph%sSpecMcGill5023GeV", particleNamesOut[iParticle].Data()), TObject::kOverwrite);
                // write vN
                for (Int_t iVn = 0; iVn < 6; iVn++){
                    if (graphVnMcGill[iCent][iParticle][iVn])
                        graphVnMcGill[iCent][iParticle][iVn]->Write(Form("graph%sV%iMcGill5023GeV", particleNamesOut[iParticle].Data(),iVn+1), TObject::kOverwrite);
                }
            }
        }

        fileTheoryGraphsPPb.mkdir("pPb_8.16TeV");
        fileTheoryGraphsPPb.cd("pPb_8.16TeV");
        // electrons from Powheg single weak boson production
        graphElecFromWeakBosonPowhegXSecpPb8TeVEMCal->Write("graphElecFromWeakBosonPowheg_pPb8TeV", TObject::kOverwrite);
        graphElecFromJJXSecpPb8TeVEMCal->Write("graphElecFromJJ_pPb8TeV", TObject::kOverwrite);
        graphRatioElecFromWeakBoson->Write("graphRatioElecFromWeakBoson_pPb8TeV", TObject::kOverwrite);
        splineRatioElecFromWeakBoson->Write("splineRatioElecFromWeakBoson_pPb8TeV", TObject::kOverwrite);

        // pi0, eta, eta/pi0 EPOSJJpPb8TeV
            histoPi0EPOSJJpPb8TeV->Write("histoPi0SpecEPOSJJpPb8TeV", TObject::kOverwrite);
            histoEtaEPOSJJpPb8TeV->Write("histoEtaSpecEPOSJJpPb8TeV", TObject::kOverwrite);
            histoPiChEPOSJJpPb8TeV->Write("histoPiChSpecEPOSJJpPb8TeV", TObject::kOverwrite);
            histoKChEPOSJJpPb8TeV->Write("histoKChSpecEPOSJJpPb8TeV", TObject::kOverwrite);
            histoPi0EPOSJJpPb8TeVReb->Write("histoPi0SpecEPOSJJpPb8TeV_Reb", TObject::kOverwrite);
            histoEtaEPOSJJpPb8TeVReb->Write("histoEtaSpecEPOSJJpPb8TeV_Reb", TObject::kOverwrite);
            histoEtaToPi0EPOSJJpPb8TeV->Write("histoEtaToPi0EPOSJJpPb8TeV", TObject::kOverwrite);
            histoEtaToKCHEPOSJJpPb8TeV->Write("histoEtaToKChEPOSJJpPb8TeV", TObject::kOverwrite);
            histoPi0ToPiCHEPOSJJpPb8TeV->Write("histoPi0ToPiChEPOSJJpPb8TeV", TObject::kOverwrite);
            histoPi0DPMJETpPb8TeV->Write("histoPi0SpecDPMJETpPb8TeV", TObject::kOverwrite);
            histoEtaDPMJETpPb8TeV->Write("histoEtaSpecDPMJETpPb8TeV", TObject::kOverwrite);
            histoPiChDPMJETpPb8TeV->Write("histoPiChSpecDPMJETpPb8TeV", TObject::kOverwrite);
            histoKChDPMJETpPb8TeV->Write("histoKChSpecDPMJETpPb8TeV", TObject::kOverwrite);
            histoPi0DPMJETpPb8TeVReb->Write("histoPi0SpecDPMJETpPb8TeV_Reb", TObject::kOverwrite);
            histoEtaDPMJETpPb8TeVReb->Write("histoEtaSpecDPMJETpPb8TeV_Reb", TObject::kOverwrite);
            histoEtaToPi0DPMJETpPb8TeV->Write("histoEtaToPi0DPMJETpPb8TeV", TObject::kOverwrite);
            histoEtaToKCHDPMJETpPb8TeV->Write("histoEtaToKChDPMJETpPb8TeV", TObject::kOverwrite);
            histoPi0ToPiCHDPMJETpPb8TeV->Write("histoPi0ToPiChDPMJETpPb8TeV", TObject::kOverwrite);

            histoPi0Pythia8EPPS16pPb8TeV->Write("histoPi0Pythia8EPPS16pPb8TeV", TObject::kOverwrite);
            histoEtaPythia8EPPS16pPb8TeV->Write("histoEtaPythia8EPPS16pPb8TeV", TObject::kOverwrite);
            histoPi0Pythia8nCTEQ15pPb8TeV->Write("histoPi0Pythia8nCTEQ15pPb8TeV", TObject::kOverwrite);
            histoEtaPythia8nCTEQ15pPb8TeV->Write("histoEtaPythia8nCTEQ15pPb8TeV", TObject::kOverwrite);


            histoPi0HIJING8160GeV->Write("histoPi0SpecHIJING_MCGenpPb8TeV", TObject::kOverwrite);
            histoPiChHIJING8160GeV->Write("histoPiChSpecHIJING_MCGenpPb8TeV", TObject::kOverwrite);
            histoEtaHIJING8160GeV->Write("histoEtaSpecHIJING_MCGenpPb8TeV", TObject::kOverwrite);
            histoEtaToPi0HIJING8160GeV->Write("histoEtaToPi0HIJING_MCGenpPb8TeV", TObject::kOverwrite);
            histoPi0ToPiCHHIJING8160GeV->Write("histoPi0ToPiChHIJING_MCGenpPb8TeV", TObject::kOverwrite);
    fileTheoryGraphsPPb.Close();

}
