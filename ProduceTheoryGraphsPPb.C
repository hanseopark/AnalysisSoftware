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
    
    TGraph* graphPi0RpAEPS09sDSS                        = new TGraph(nlinesEPSsPi0fDSS,xEPSsPi0fDSS,yEPSsPi0fDSS);
    TGraphAsymmErrors* graphPi0RpAAsymmErrEPS09sDSS     = new TGraphAsymmErrors(nlinesEPSsPi0fDSS,xEPSsPi0fDSS,yEPSsPi0fDSS,xDownErrorEPSsPi0DSS,xUpErrorEPSsPi0DSS,yDownErrorEPSsPi0DSS,yUpErrorEPSsPi0DSS);
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
    TGraph* graphPi0RpACGC = new TGraph(nlinesPi0CGC+1,xESPsPi0CGC,yESPsPi0CGC);	    

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
    TH1D* histoEtaDPMJet            = (TH1D*)fileDPMJet->Get("MC_Eta_Pt");
    TH1D* histoEtaDPMJetReb         = (TH1D*)fileDPMJet->Get("MC_Eta_Pt_Rebinned");
    TH1D* histoEtaToPi0DPMJet       = (TH1D*)fileDPMJet->Get("MCEtaToPi0");
    
    //**************************************************************************************************
    //********************** HIJING MC spectra and ratio ***********************************************
    //**************************************************************************************************
    // file generated with TaskV1/ExtractMCInputSpectraFromFile.C++ based on PCM-EMC inputs
    TFile* fileHIJING               = new TFile("ExternalInputpPb/Theory/MCInputCompilationLHC13e7_pPb5TeV_2.root");
    TH1D* histoPi0HIJING            = (TH1D*)fileHIJING->Get("MC_Pi0_Pt");
    TH1D* histoPi0HIJINGReb         = (TH1D*)fileHIJING->Get("MC_Pi0_Pt_Rebinned");
    TH1D* histoEtaHIJING            = (TH1D*)fileHIJING->Get("MC_Eta_Pt");
    TH1D* histoEtaHIJINGReb         = (TH1D*)fileHIJING->Get("MC_Eta_Pt_Rebinned");
    TH1D* histoEtaToPi0HIJING       = (TH1D*)fileHIJING->Get("MCEtaToPi0");
    
    
    //**********************************************************************************************************************
    //********************************* Write graphs and histos to compilation file for pPb ********************************
    //**********************************************************************************************************************
    TFile fileTheoryGraphsPPb("ExternalInputpPb/Theory/TheoryCompilationPPb.root","UPDATE");

        graphPi0RpACGC->Write("graphPi0RpACGC5023GeV", TObject::kOverwrite);
        graphPi0RpAEPS09sKKP->Write("graphPi0RpAEPS09sKKP5023GeV", TObject::kOverwrite);
        graphPi0RpAEPS09sAKK->Write("graphPi0RpAEPS09sAKK5023GeV", TObject::kOverwrite);
        graphPi0RpAEPS09sDSS->Write("graphPi0RpAEPS09sDSS5023GeV", TObject::kOverwrite);
        graphPi0RpAAsymmErrEPS09sDSS->Write("graphPi0RpAAsymmErrEPS09sDSS5023GeV", TObject::kOverwrite);
        histoPi0EPOS3->Write("histoPi0EPOS35023TeV", TObject::kOverwrite);
        histoEtaEPOS3->Write("histoEtaEPOS35023TeV", TObject::kOverwrite);
        histoPi0EPOS3Reb->Write("histoPi0EPOS35023TeV_Reb", TObject::kOverwrite);
        histoEtaEPOS3Reb->Write("histoEtaEPOS35023TeV_Reb", TObject::kOverwrite);
        histoEtaToPi0EPOS3->Write("histoEtaToPi0EPOS35023TeV", TObject::kOverwrite);
        histoEtaToPi0EPOS3Reb->Write("histoEtaToPi0EPOS35023TeV_Reb", TObject::kOverwrite);
        histoPi0DPMJet->Write("histoPi0DPMJet5023TeV", TObject::kOverwrite);
        histoEtaDPMJet->Write("histoEtaDPMJet5023TeV", TObject::kOverwrite);
        histoPi0DPMJetReb->Write("histoPi0DPMJet5023TeV_Reb", TObject::kOverwrite);
        histoEtaDPMJetReb->Write("histoEtaDPMJet5023TeV_Reb", TObject::kOverwrite);
        histoEtaToPi0DPMJet->Write("histoEtaToPi0DPMJet5023TeV", TObject::kOverwrite);
        histoPi0HIJING->Write("histoPi0HIJING5023TeV", TObject::kOverwrite);
        histoEtaHIJING->Write("histoEtaHIJING5023TeV", TObject::kOverwrite);
        histoPi0HIJINGReb->Write("histoPi0HIJING5023TeV_Reb", TObject::kOverwrite);
        histoEtaHIJINGReb->Write("histoEtaHIJING5023TeV_Reb", TObject::kOverwrite);
        histoEtaToPi0HIJING->Write("histoEtaToPi0HIJING5023TeV", TObject::kOverwrite);
        
    fileTheoryGraphsPPb.Close();

}
