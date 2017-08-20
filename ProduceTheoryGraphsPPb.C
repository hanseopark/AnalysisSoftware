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

Double_t xSection5023GeVINEL        = 67.6*1e-3;    // from https://aliceinfo.cern.ch/Notes/sites/aliceinfo.cern.ch.Notes/files/notes/analysis/stripath/2017-Jun-14-analysis_note-INEL_norm.pdf
Double_t xSection5023GeVINELErr     = 0.6;          // from https://aliceinfo.cern.ch/Notes/sites/aliceinfo.cern.ch.Notes/files/notes/analysis/stripath/2017-Jun-14-analysis_note-INEL_norm.pdf

Double_t xSection5023GeVINELpPb     = 70*1e-3;
Double_t ncollpPb5023GeV            = 6.9;
Double_t recalcBarn                 = 1e12; //NLO in pbarn!!!!

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


//**********************************************************************************************************
// Calculates the ratio of two error graphs
//**********************************************************************************************************
TGraphErrors* CalculateGraphRatioToGraph(TGraphErrors* graphA, TGraphErrors* graphB){

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
    // file generated with TaskV1/ExtractMCInputSpectraFromFile.C++ based on PCM-EMC inputs
    TFile* fileDPMJet               = new TFile("ExternalInputpPb/Theory/MCInputCompilationLHC13b2_efix_pPb5TeV_2.root");
    TH1D* histoPi0DPMJet            = (TH1D*)fileDPMJet->Get("MC_Pi0_Pt");
    TH1D* histoPi0DPMJetReb         = (TH1D*)fileDPMJet->Get("MC_Pi0_Pt_Rebinned");
    TH1D* histoPiChDPMJet           = (TH1D*)fileDPMJet->Get("MC_PiCh_All_Pt");
    TH1D* histoKChDPMJet            = (TH1D*)fileDPMJet->Get("MC_KCh_All_Pt");
    TH1D* histoEtaDPMJet            = (TH1D*)fileDPMJet->Get("MC_Eta_Pt");
    TH1D* histoEtaDPMJetReb         = (TH1D*)fileDPMJet->Get("MC_Eta_Pt_Rebinned");
    TH1D* histoEtaToPi0DPMJet       = (TH1D*)fileDPMJet->Get("MCEtaToPi0");
    TH1D* histoPi0ToPiCHDPMJet      = (TH1D*)fileDPMJet->Get("MCPi0ToPiCh");
    TH1D* histoEtaToKCHDPMJet       = (TH1D*)fileDPMJet->Get("MCEtaToKCh");

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
    TGraphErrors* graphEtaToPi0McGill = CalculateGraphRatioToGraph(graphEtaSpecMcGill, graphPi0SpecMcGill);

    //**************************************************************************************************
    //********************* extracting McGill predictions including other observables ******************
    //**************************************************************************************************
    Int_t nParticles                    = 14;
    Int_t nMeasurements                 = 15;
    Int_t nCent                         = 3;
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
        valuesMcGill[iCent]             = new vector<Double_t>*[nMeasurements];
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
    graphNLOCalcInvYieldEtaMuOnepPb5023GeV          = ScaleGraph(graphNLOCalcInvSecEtaMuOne5023GeV, 1/(xSection5023GeVINELpPb*recalcBarn)*ncollpPb5023GeV);
    graphNLOCalcInvYieldEtaMuOnepp5023GeV           = ScaleGraph(graphNLOCalcInvSecEtaMuOne5023GeV, 1/(xSection5023GeVINEL*recalcBarn));

    if (fillAllMuScalesEta5023GeV){
        // fill x-section graphs for mu variations
        graphNLOCalcInvSecEtaMuHalf5023GeV          = new TGraph(nlinesNLOEta5023GeV,ptNLOEta5023GeV,muHalfEta5023GeV);
        graphNLOCalcInvSecEtaMuHalf5023GeV->RemovePoint(0);
        graphNLOCalcInvSecEtaMuTwo5023GeV           = new TGraph(nlinesNLOEta5023GeV,ptNLOEta5023GeV,muTwoEta5023GeV);
        // fill yield graphs for mu variations
        graphNLOCalcInvYieldEtaMuHalfpPb5023GeV     = ScaleGraph(graphNLOCalcInvSecEtaMuHalf5023GeV, 1/(xSection5023GeVINELpPb*recalcBarn)*ncollpPb5023GeV);
        graphNLOCalcInvYieldEtaMuTwopPb5023GeV      = ScaleGraph(graphNLOCalcInvSecEtaMuTwo5023GeV, 1/(xSection5023GeVINELpPb*recalcBarn)*ncollpPb5023GeV);
        graphNLOCalcInvYieldEtaMuHalfpp5023GeV      = ScaleGraph(graphNLOCalcInvSecEtaMuHalf5023GeV, 1/(xSection5023GeVINEL*recalcBarn));
        graphNLOCalcInvYieldEtaMuTwopp5023GeV       = ScaleGraph(graphNLOCalcInvSecEtaMuTwo5023GeV, 1/(xSection5023GeVINEL*recalcBarn));
    }
    // combine all mu scales into 1 graph for eventual plotting
    TGraphAsymmErrors* graphNLOCalcInvSecEta5023GeV     = CombineMuScales(nlinesNLOEta5023GeV, ptNLOEta5023GeV, muOneEta5023GeV, muHalfEta5023GeV, muTwoEta5023GeV);
    TGraphAsymmErrors* graphNLOCalcInvYieldEtapPb5023GeV= ScaleGraphAsym(graphNLOCalcInvSecEta5023GeV, 1/(xSection5023GeVINELpPb*recalcBarn)*ncollpPb5023GeV);
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
    TGraph* graphNLOCalcInvYieldPi0MuHalfpPb5023GeV     = ScaleGraph(graphNLOCalcInvSecPi0MuHalf5023GeV, 1/(xSection5023GeVINELpPb*recalcBarn)*ncollpPb5023GeV);
    TGraph* graphNLOCalcInvYieldPi0MuOnepPb5023GeV      = ScaleGraph(graphNLOCalcInvSecPi0MuOne5023GeV, 1/(xSection5023GeVINELpPb*recalcBarn)*ncollpPb5023GeV);
    TGraph* graphNLOCalcInvYieldPi0MuTwopPb5023GeV      = ScaleGraph(graphNLOCalcInvSecPi0MuTwo5023GeV, 1/(xSection5023GeVINELpPb*recalcBarn)*ncollpPb5023GeV);
    TGraph* graphNLOCalcInvYieldPi0MuHalfpp5023GeV      = ScaleGraph(graphNLOCalcInvSecPi0MuHalf5023GeV, 1/(xSection5023GeVINEL*recalcBarn));
    TGraph* graphNLOCalcInvYieldPi0MuOnepp5023GeV       = ScaleGraph(graphNLOCalcInvSecPi0MuOne5023GeV, 1/(xSection5023GeVINEL*recalcBarn));
    TGraph* graphNLOCalcInvYieldPi0MuTwopp5023GeV       = ScaleGraph(graphNLOCalcInvSecPi0MuTwo5023GeV, 1/(xSection5023GeVINEL*recalcBarn));
    // combine all mu scales into 1 graph for eventual plotting
    TGraphAsymmErrors* graphNLOCalcInvSecPi05023GeV     = CombineMuScales(nlinesNLOPi05023GeV, ptNLOPi05023GeV, muOnePi05023GeV, muHalfPi05023GeV, muTwoPi05023GeV);
    TGraphAsymmErrors* graphNLOCalcInvYieldPi0pPb5023GeV= ScaleGraphAsym(graphNLOCalcInvSecPi05023GeV, 1/(xSection5023GeVINELpPb*recalcBarn)*ncollpPb5023GeV);
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
    Double_t ptNLOEtanPdf5023GeV[200];
    Double_t muOneEtanPdf5023GeV[200];
    Int_t nlinesNLOEtanPdf5023GeV                       = 0;

    TString fileNameNLOEtanPdf5023GeV                   = "ExternalInputpPb/Theory/Vogelsang/ALICENLOcalcEtaVogelsang5023GeV_AESSS_nCTEQ.dat"; // currently only mu = 1pt
    ifstream  fileNLOEtanPdf5023GeV;
    fileNLOEtanPdf5023GeV.open(fileNameNLOEtanPdf5023GeV,ios_base::in);
    cout << fileNameNLOEtanPdf5023GeV << endl;

    // read eta input theory file for 5.023TeV
    while(!fileNLOEtanPdf5023GeV.eof() && nlinesNLOEtanPdf5023GeV < 200){
        Double_t dummy                              = 0;
        fileNLOEtanPdf5023GeV >> dummy >> dummy >> ptNLOEtanPdf5023GeV[nlinesNLOEtanPdf5023GeV] >> muOneEtanPdf5023GeV[nlinesNLOEtanPdf5023GeV] ;
        nlinesNLOEtanPdf5023GeV++;
    }
    fileNLOEtanPdf5023GeV.close();
    // reduce number of lines by 1
    nlinesNLOEtanPdf5023GeV--;

    // generate graphs
    TGraph* graphNLOCalcInvSecEtaMuOnenPdf5023GeV       = NULL;
    TGraph* graphNLOCalcInvYieldEtaMuOnepPbnPdf5023GeV  = NULL;

    // fill x-section graph for default mu
    graphNLOCalcInvSecEtaMuOnenPdf5023GeV               = new TGraph(nlinesNLOEtanPdf5023GeV,ptNLOEtanPdf5023GeV,muOneEtanPdf5023GeV);
    graphNLOCalcInvSecEtaMuOnenPdf5023GeV->Print();
    // fill yield graph for default mu
    graphNLOCalcInvYieldEtaMuOnepPbnPdf5023GeV          = ScaleGraph(graphNLOCalcInvSecEtaMuOnenPdf5023GeV, 1/(xSection5023GeVINELpPb*recalcBarn)*ncollpPb5023GeV);

    //*********************************************************
    // pi0 Vogelsang - FF: DSS14, PDF:  nPDF: nCTEQ
    Double_t ptNLOPi0nPdf5023GeV[200];
    Double_t muOnePi0nPdf5023GeV[200];
    Int_t nlinesNLOPi0nPdf5023GeV                       = 0;

    TString fileNameNLOPi0nPdf5023GeV                   = "ExternalInputpPb/Theory/Vogelsang/ALICENLOcalcPi0Vogelsang5023GeV_DSS14_nCTEQ.dat"; // currently only mu = 1pt
    ifstream  fileNLOPi0nPdf5023GeV;
    fileNLOPi0nPdf5023GeV.open(fileNameNLOPi0nPdf5023GeV,ios_base::in);
    cout << fileNameNLOPi0nPdf5023GeV << endl;

    while(!fileNLOPi0nPdf5023GeV.eof() && nlinesNLOPi0nPdf5023GeV < 200){
        Double_t dummy                              = 0;
        fileNLOPi0nPdf5023GeV >> dummy >> dummy >> ptNLOPi0nPdf5023GeV[nlinesNLOPi0nPdf5023GeV] >>  muOnePi0nPdf5023GeV[nlinesNLOPi0nPdf5023GeV] ;
        nlinesNLOPi0nPdf5023GeV++;
    }
    fileNLOPi0nPdf5023GeV.close();
    // reduce number of lines by 1
    nlinesNLOPi0nPdf5023GeV--;

    // fill x-section graphs
    TGraph* graphNLOCalcInvSecPi0MuOnenPdf5023GeV       = new TGraph(nlinesNLOPi0nPdf5023GeV,ptNLOPi0nPdf5023GeV,muOnePi0nPdf5023GeV);
    // fill yield graphs
    TGraph* graphNLOCalcInvYieldPi0MuOnepPbnPdf5023GeV  = ScaleGraph(graphNLOCalcInvSecPi0MuOnenPdf5023GeV, 1/(xSection5023GeVINELpPb*recalcBarn)*ncollpPb5023GeV);

    //*********************************************************
    // Eta/Pi0 Vogelsang FF pi0: DSS14, FF eta: AESSS,  nPDF: nCTEQ
    Double_t* valueNLOMuOneEtanPdf5023GeV               = graphNLOCalcInvSecEtaMuOnenPdf5023GeV->GetY();
    Double_t* valueNLOMuOnePi0nPdf5023GeV               = graphNLOCalcInvSecPi0MuOnenPdf5023GeV->GetY();
    Double_t* xValueNLOnPdf5023GeV                      = graphNLOCalcInvSecPi0MuOnenPdf5023GeV->GetX();
    Int_t xNBinsnPdf5023GeV                             = graphNLOCalcInvSecPi0MuOnenPdf5023GeV->GetN();
    Double_t valueNLOEtaToPi0NLOMuOnenPdf5023GeV[200];

    // calculate eta/pi0 ratios for different mu scales
    for ( Int_t n = 0; n < xNBinsnPdf5023GeV+1; n++){
        if (valueNLOMuOnePi0nPdf5023GeV[n] != 0){
            valueNLOEtaToPi0NLOMuOnenPdf5023GeV[n] = valueNLOMuOneEtanPdf5023GeV[n]/valueNLOMuOnePi0nPdf5023GeV[n];
        } else {
            valueNLOEtaToPi0NLOMuOnenPdf5023GeV[n] = 0.;
        }
    }
    // fill graphs
    TGraph* graphEtaToPi0NLOMuOnenPdf5023GeV            = new TGraph(xNBinsnPdf5023GeV,xValueNLOnPdf5023GeV,valueNLOEtaToPi0NLOMuOnenPdf5023GeV);




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

        // write MB calcs
        fileTheoryGraphsPPb.cd("pPb_5.023TeV");
            // pi0 CGC RpPb
            graphPi0RpACGC->Write("graphPi0RpACGC5023GeV", TObject::kOverwrite);
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
            // pi0, eta, eta/pi0 DPMJet
            histoPi0DPMJet->Write("histoPi0SpecDPMJet5023GeV", TObject::kOverwrite);
            histoEtaDPMJet->Write("histoEtaSpecDPMJet5023GeV", TObject::kOverwrite);
            histoPiChDPMJet->Write("histoPiChSpecDPMJet5023GeV", TObject::kOverwrite);
            histoKChDPMJet->Write("histoKChSpecDPMJet5023GeV", TObject::kOverwrite);
            histoPi0DPMJetReb->Write("histoPi0SpecDPMJet5023GeV_Reb", TObject::kOverwrite);
            histoEtaDPMJetReb->Write("histoEtaSpecDPMJet5023GeV_Reb", TObject::kOverwrite);
            histoEtaToPi0DPMJet->Write("histoEtaToPi0DPMJet5023GeV", TObject::kOverwrite);
            histoEtaToKCHDPMJet->Write("histoEtaToKChDPMJet5023GeV", TObject::kOverwrite);
            histoPi0ToPiCHDPMJet->Write("histoPi0ToPiChDPMJet5023GeV", TObject::kOverwrite);
            // pi0, eta, eta/pi0 HIJING
            histoPi0HIJING->Write("histoPi0SpecHIJING5023GeV", TObject::kOverwrite);
            histoEtaHIJING->Write("histoEtaSpecHIJING5023GeV", TObject::kOverwrite);
            histoPiChHIJING->Write("histoPiChSpecHIJING5023GeV", TObject::kOverwrite);
            histoKChHIJING->Write("histoKChSpecHIJING5023GeV", TObject::kOverwrite);
            histoPi0HIJINGReb->Write("histoPi0SpecHIJING5023GeV_Reb", TObject::kOverwrite);
            histoEtaHIJINGReb->Write("histoEtaSpecHIJING5023GeV_Reb", TObject::kOverwrite);
            histoEtaToPi0HIJING->Write("histoEtaToPi0HIJING5023GeV", TObject::kOverwrite);
            histoEtaToKCHHIJING->Write("histoEtaToKChHIJING5023GeV", TObject::kOverwrite);
            histoPi0ToPiCHHIJING->Write("histoPi0ToPiChHIJING5023GeV", TObject::kOverwrite);
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
            graphNLOCalcInvYieldPi0MuOnepPbnPdf5023GeV->Write("graphNLOCalcDSS14InvYieldPi0MuOne5023GeV_nCTEQ", TObject::kOverwrite);
            // eta Vogelsang nPDF: nCTEQ, FF AESSS
            graphNLOCalcInvYieldEtaMuOnepPbnPdf5023GeV->Write("graphNLOCalcAESSSInvYieldEtaMuOne5023GeV_nCTEQ", TObject::kOverwrite);
            // eta/pi0 Vogelsang nPDF: nCTEQ, FF eta- AESSS, FF pi0 - DSS14
            graphEtaToPi0NLOMuOnenPdf5023GeV->Write("graphNLOCalcEtaOverPi0MuOne5023GeV_AESSS_DSS14_nCTEQ", TObject::kOverwrite);

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
    fileTheoryGraphsPPb.Close();

}
