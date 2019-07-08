#include <stdlib.h>
#include <iostream>
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
#include "TH3F.h"
#include "TF1.h"
#include "TExec.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TTree.h"
#include "TMinuit.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TVectorT.h"
#include "TArc.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/CombinationFunctions.h"
typedef TVectorT<double> TVectorD;
typedef TVectorT<float> TVectorF;

using namespace std;

void PalBW()
{
   static Int_t  colors[50];
   static Bool_t initialized = kFALSE;

    Double_t stopsBW[5] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t redBW[5]   = { 1.00, 0.84, 0.61, 0.34, 0.00 };
    Double_t greenBW[5] = { 1.00, 0.84, 0.61, 0.34, 0.00 };
    Double_t blueBW[5]  = { 1.00, 0.84, 0.61, 0.34, 0.00 };

   if(!initialized){
      Int_t FI = TColor::CreateGradientColorTable(5, stopsBW, redBW, greenBW, blueBW, 50);
      for (int i=0; i<50; i++) colors[i] = FI+i;
      initialized = kTRUE;
      return;
   }
   gStyle->SetPalette(50,colors);
}

void PalColor()
{
   static Int_t  colors[50];
   static Bool_t initialized = kFALSE;

    Double_t stops[5] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[5]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[5] = { 0.31, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[5]  = { 0.51, 1., 0.12, 0.00, 0.00};

   if(!initialized){
      Int_t FI = TColor::CreateGradientColorTable(5, stops, red, green, blue, 50);
      for (int i=0; i<50; i++) colors[i] = FI+i;
      initialized = kTRUE;
      return;
   }
   gStyle->SetPalette(50,colors);
}

void MakeClusterQAPlots(
    TString fileNameData    = "data.root",
    TString fileNameMC      = "MC.root",
    TString fileNameMCXT    = "MCXT.root",
    TString energy          = "5TeV2017",
    TString suffix          = "pdf",
    Double_t minEPlot          = 20.,
    Double_t maxEPlot          = 100.,
    Int_t minSM             = 0,
    Int_t maxSM             = 20,
    Int_t rebinFac             = 2,
    TString labelMC1 = "JJ+GJ std",
    TString labelMC2 = "JJ+GJ XTalk"
    // TString labelMC1 = "JJ XTalk";
    // TString labelMC2 = "JJ+GJ XTalk";
    // TString labelMC1 = "MC";
    // TString labelMC2 = "MC + XTalk";
    // TString labelMC1 = "GJ std";
    // TString labelMC2 = "GJ XTalk";
    // TString labelMC1 = "JJ std";
    // TString labelMC2 = "JJ XTalk";
){

    gROOT->Reset();
    gROOT->SetStyle("Plain");

    StyleSettingsThesis();
    SetPlotStyle();


    //********************************************************************************
    //*            File definition/ loading                                          *
    //********************************************************************************

    TFile *inputFile[3]   = {NULL};
    inputFile[0]   = new TFile(fileNameData.Data());
    if (inputFile[0]->IsZombie()) {
        cout <<fileNameData.Data() <<" file does not exist" << endl;
        inputFile[0]->Close();
        delete inputFile[0];
        return;
    }

    inputFile[1]   = new TFile(fileNameMC.Data());
    if (inputFile[1]->IsZombie()) {
        cout <<fileNameData.Data() <<" file does not exist" << endl;
        inputFile[1]->Close();
        delete inputFile[1];
        return;
    }

    inputFile[2]   = new TFile(fileNameMCXT.Data());
    if (inputFile[2]->IsZombie()) {
        cout <<fileNameData.Data() <<" file does not exist" << endl;
        inputFile[2]->Close();
        delete inputFile[2];
        return;
    }

  //___________________________________ Declaration of files _____________________________________________
    TString dateForOutput                           = ReturnDateStringForOutput();
    cout << dateForOutput.Data() << endl;
    TString collisionSystem                    = ReturnFullCollisionsSystem(energy);
    TString labelALICEforPlots                      = "ALICE";
    TString outputDir                               = Form("%s/%s/MakeClusterQAPlots_%2.1f_%2.1f",suffix.Data(),dateForOutput.Data(), minEPlot, maxEPlot);
    cout << outputDir.Data() << endl;

    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec("mkdir -p "+outputDir + "/Ratios");
    gSystem->Exec("mkdir -p "+outputDir + "/Ratios/withFits");
    // gSystem->Exec("mkdir -p "+outputDir + "/Chi2PsiPair");

    //********************************************************************************
    //*      Definition of histograms for reconstructed Conversion Points            *
    //********************************************************************************
    // Color_t  colorSMdata[20] ={ kBlack,  kBlue+2, kViolet+2, kMagenta+2, kPink+2, kRed+2, kOrange+2, kYellow+2, kSpring+2, kGreen+2, kTeal+2, kCyan+2, kAzure+2, kBlue+3, kMagenta+3, kRed+3, kYellow+3, kGreen+3, kCyan+3, kGray+3};
    Color_t  colorSMdata[20] ={ kBlack,  kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack};
    Color_t  colorSMMC[20] ={ kGreen+2,  kGreen+2, kGreen+2, kGreen+2, kGreen+2, kGreen+2, kGreen+2, kGreen+2, kGreen+2, kGreen+2, kGreen+2, kGreen+2, kGreen+2, kGreen+2, kGreen+2, kGreen+2, kGreen+2, kGreen+2, kGreen+2, kGreen+2};
    Color_t  colorSMMCXT[20] ={ kRed+2,  kRed+2, kRed+2, kRed+2, kRed+2, kRed+2, kRed+2, kRed+2, kRed+2, kRed+2, kRed+2, kRed+2, kRed+2, kRed+2, kRed+2, kRed+2, kRed+2, kRed+2, kRed+2, kRed+2};
    // Color_t  colorSMMC[20] ={ kGreen+2, kViolet-8, kMagenta-8, kPink-8, kRed-8, kOrange-8, kYellow-8, kSpring-8, kGreen-8, kTeal-8, kCyan-8, kAzure-8, kBlue-4, kMagenta-4, kRed-4, kYellow-4, kGreen-4, kCyan-4, kGray-4, kGray+2};
    // Color_t  colorSMMCXT[20]   ={ kRed+2, kBlue-8, kViolet-8, kMagenta-8, kPink-8, kRed-8, kOrange-8, kYellow-8, kSpring-8, kGreen-8, kTeal-8, kCyan-8, kAzure-8, kBlue-4, kMagenta-4, kRed-4, kYellow-4, kGreen-4, kCyan-4, kGray-4};
    Color_t  colorTruePart[20]      ={ kBlack,  kBlue+2, kViolet+2, kMagenta+2, kPink+2, kRed+2, kOrange+2, kYellow+2, kSpring+2, kGreen+2, kTeal+2, kCyan+2, kAzure+2, kBlue+3, kMagenta+3, kRed+3, kYellow+3, kGreen+3, kCyan+3, kGray+3};
    Color_t  colorTruePart2[20]     ={ kRed+2, kBlue-8, kViolet-8, kMagenta-8, kPink-8, kRed-8, kOrange-8, kYellow-8, kSpring-8, kGreen-8, kTeal-8, kCyan-8, kAzure-8, kBlue-4, kMagenta-4, kRed-4, kYellow-4, kGreen-4, kCyan-4, kGray-4};
    Style_t  markerSMdata[20]={20, 21, 23, 29, 33, 47, 39, 41, 45, 47, 20, 21, 23, 29, 33, 47, 39, 41, 45, 47};
    Style_t  markerSMMC[20]={24, 32, 30, 27, 46, 37, 40, 44, 46, 24, 25, 32, 30, 27, 46, 37, 40, 44, 46, 24};
    Style_t  markerSMMCXT[20]  ={24, 25, 32, 30, 27, 46, 37, 40, 44, 46, 24, 25, 32, 30, 27, 46, 37, 40, 44, 46};

    const Int_t nBinsE = 14;
    Double_t projBinnin[nBinsE] = {
        0.00, 5.00, 8.00, 10.0, 15.0,
        20.0, 30.0, 40.0, 50.0, 75.0,
        100., 125., 150., 200.};

    TH1F* histoM02allSM[3]              = {NULL};
    TH1F* histoM02allSMPlot[3]              = {NULL};
    TH1F* histoM02RatioallSM[2]         = {NULL};
    TH1F* histoM02RatioallSMMC         = NULL;
    TH1D* histoM02perSM[20][3]          = {{NULL}};
    TH1F* histoM02RatioperSM[20][2]     = {{NULL}};
    TH2F* histoM02perSMvsE[20][3]       = {{NULL}};
    TH1F* histoM20perSM[20][3]          = {{NULL}};
    TH2F* histoM20perSMvsE[20][3]       = {{NULL}};
    TH2F* histoM02vsM20perSM[20][3]     = {{NULL}};
    TH3F* histoM02vsM20perSMvsE[20][3]  = {{NULL}};

    TH2F* histoEFracFirstLabelvsE[2]       = {NULL};
    TH2F* histoEFracLeadPi0vsE[2]          = {NULL};

    TH1F* histoM02TrueMergedPi0allSM[2]                 = {NULL};
    TH1F* histoM02FracTrueMergedPi0allSM[2]                 = {NULL};
    TH1D* histoM02TrueMergedPi0perSM[20][2]                 = {{NULL}};
    // TH1F* histoM02TrueMergedEtaallSM[2]                 = {NULL};
    // TH1F* histoM02FracTrueMergedEtaallSM[2]                 = {NULL};
    // TH1D* histoM02TrueMergedEtaperSM[20][2]                 = {{NULL}};
    TH1F* histoM02TrueMergedPartConvPi0allSM[2]                 = {NULL};
    TH1F* histoM02FracTrueMergedPartConvPi0allSM[2]                 = {NULL};
    TH1D* histoM02TrueMergedPartConvPi0perSM[20][2]         = {{NULL}};
    // TH1F* histoM02TrueMergedPartConvEtaallSM[2]                 = {NULL};
    // TH1F* histoM02FracTrueMergedPartConvEtaallSM[2]                 = {NULL};
    // TH1D* histoM02TrueMergedPartConvEtaperSM[20][2]         = {{NULL}};
    TH1F* histoM02TrueGammaFromPi0allSM[2]                 = {NULL};
    TH1F* histoM02FracTrueGammaFromPi0allSM[2]                 = {NULL};
    TH1D* histoM02TrueGammaFromPi0perSM[20][2]              = {{NULL}};
    // TH1F* histoM02TrueGammaFromEtaallSM[2]                 = {NULL};
    // TH1F* histoM02FracTrueGammaFromEtaallSM[2]                 = {NULL};
    // TH1D* histoM02TrueGammaFromEtaperSM[20][2]              = {{NULL}};
    TH1F* histoM02TrueElectronFromPi0allSM[2]                 = {NULL};
    TH1F* histoM02FracTrueElectronFromPi0allSM[2]                 = {NULL};
    TH1D* histoM02TrueElectronFromPi0perSM[20][2]           = {{NULL}};
    // TH1F* histoM02TrueElectronFromEtaallSM[2]                 = {NULL};
    // TH1F* histoM02FracTrueElectronFromEtaallSM[2]                 = {NULL};
    // TH1D* histoM02TrueElectronFromEtaperSM[20][2]           = {{NULL}};
    TH1F* histoM02TrueBackgroundallSM[2]                 = {NULL};
    TH1F* histoM02FracTrueBackgroundallSM[2]                 = {NULL};
    TH1D* histoM02TrueBackgroundperSM[20][2]                = {{NULL}};

    TH2F* histoM02TrueMergedPi0perSMvsE[20][2]             = {{NULL}};
    // TH2F* histoM02TrueMergedEtaperSMvsE[20][2]             = {{NULL}};
    TH2F* histoM02TrueMergedPartConvPi0perSMvsE[20][2]     = {{NULL}};
    // TH2F* histoM02TrueMergedPartConvEtaperSMvsE[20][2]     = {{NULL}};
    TH2F* histoM02TrueGammaFromPi0perSMvsE[20][2]          = {{NULL}};
    // TH2F* histoM02TrueGammaFromEtaperSMvsE[20][2]          = {{NULL}};
    TH2F* histoM02TrueElectronFromPi0perSMvsE[20][2]       = {{NULL}};
    // TH2F* histoM02TrueElectronFromEtaperSMvsE[20][2]       = {{NULL}};
    TH2F* histoM02TrueBackgroundperSMvsE[20][2]            = {{NULL}};

    Double_t projectMin;
    Double_t projectMax;
    for(Int_t kMC=0;kMC<3;kMC++){

        for(Int_t iSM=minSM; iSM<maxSM; iSM++){

            histoM02perSMvsE[iSM][kMC]             = (TH2F*)inputFile[kMC]->Get(Form("histoM02perSMvsE%d",iSM));
            projectMin            = histoM02perSMvsE[iSM][kMC]->GetYaxis()->FindBin(minEPlot+0.001);
            projectMax            = histoM02perSMvsE[iSM][kMC]->GetYaxis()->FindBin(maxEPlot-0.001);

            histoM02perSM[iSM][kMC]             = (TH1D*) histoM02perSMvsE[iSM][kMC]->ProjectionX(Form("histoProjM02perSM%d%d",iSM,kMC),projectMin,projectMax);
            histoM02perSM[iSM][kMC]->Sumw2();
            histoM02perSM[iSM][kMC]->Rebin(rebinFac);
            if(iSM==minSM){
                histoM02allSM[kMC] = (TH1F*)histoM02perSM[iSM][kMC]->Clone(Form("histoM02allSM_%d",kMC));
                histoM02allSMPlot[kMC] = (TH1F*)histoM02perSM[iSM][kMC]->Clone(Form("histoM02allSMPlot_%d",kMC));
            } else {
                histoM02allSM[kMC]->Add(histoM02perSM[iSM][kMC]);
                histoM02allSMPlot[kMC]->Add(histoM02perSM[iSM][kMC]);
            }
            histoM02perSM[iSM][kMC]->Scale(1/histoM02perSM[iSM][kMC]->Integral());
            if(kMC>0){
                histoM02RatioperSM[iSM][kMC-1] = (TH1F*)histoM02perSM[iSM][0]->Clone(Form("histoM02RatioperSM%d%d",iSM,kMC));
                histoM02RatioperSM[iSM][kMC-1]->Sumw2();
                histoM02RatioperSM[iSM][kMC-1]->Divide(histoM02perSM[iSM][kMC]);
            }
        }

        if (kMC>0){
            histoEFracFirstLabelvsE[kMC-1]             = (TH2F*)inputFile[kMC]->Get("histoEFracFirstLabelvsE");
            histoEFracLeadPi0vsE[kMC-1]                = (TH2F*)inputFile[kMC]->Get("histoEFracLeadPi0vsE");

            for(Int_t iSM=minSM; iSM<maxSM; iSM++){
                // true merged pi0
                histoM02TrueMergedPi0perSMvsE[iSM][kMC-1]  = (TH2F*)inputFile[kMC]->Get(Form("histoM02TrueMergedPi0perSMvsE%d",iSM));
                projectMin            = histoM02TrueMergedPi0perSMvsE[iSM][kMC-1]->GetYaxis()->FindBin(minEPlot+0.001);
                projectMax            = histoM02TrueMergedPi0perSMvsE[iSM][kMC-1]->GetYaxis()->FindBin(maxEPlot-0.001);
                histoM02TrueMergedPi0perSM[iSM][kMC-1] = (TH1D*) histoM02TrueMergedPi0perSMvsE[iSM][kMC-1]->ProjectionX(Form("histoM02TrueMergedPi0perSM%d%d%d",iSM,kMC-1,iSM),projectMin,projectMax);
                histoM02TrueMergedPi0perSM[iSM][kMC-1]->Rebin(rebinFac);
                histoM02TrueMergedPi0perSM[iSM][kMC-1]->Sumw2();
                if(iSM==minSM){
                    histoM02TrueMergedPi0allSM[kMC-1] = (TH1F*)histoM02TrueMergedPi0perSM[iSM][kMC-1]->Clone(Form("histoM02TrueMergedPi0allSM%d",kMC-1));
                } else {
                    histoM02TrueMergedPi0allSM[kMC-1]->Add(histoM02TrueMergedPi0perSM[iSM][kMC-1]);
                }
                histoM02TrueMergedPi0perSM[iSM][kMC-1]->Scale(1/histoM02TrueMergedPi0perSM[iSM][kMC-1]->Integral());

                // true merged eta
                // histoM02TrueMergedEtaperSMvsE[iSM][kMC-1]  = (TH2F*)inputFile[kMC]->Get(Form("histoM02TrueMergedEtaperSMvsE%d",iSM));
                // projectMin            = histoM02TrueMergedEtaperSMvsE[iSM][kMC-1]->GetYaxis()->FindBin(minEPlot+0.001);
                // projectMax            = histoM02TrueMergedEtaperSMvsE[iSM][kMC-1]->GetYaxis()->FindBin(maxEPlot-0.001);
                // histoM02TrueMergedEtaperSM[iSM][kMC-1] = (TH1D*) histoM02TrueMergedEtaperSMvsE[iSM][kMC-1]->ProjectionX(Form("histoM02TrueMergedEtaperSM%d%d%d",iSM,kMC-1,iSM),projectMin,projectMax);
                // histoM02TrueMergedEtaperSM[iSM][kMC-1]->Rebin(rebinFac);
                // histoM02TrueMergedEtaperSM[iSM][kMC-1]->Sumw2();
                // if(iSM==minSM){
                //     histoM02TrueMergedEtaallSM[kMC-1] = (TH1F*)histoM02TrueMergedEtaperSM[iSM][kMC-1]->Clone(Form("histoM02TrueMergedEtaallSM%d",kMC-1));
                // } else {
                //     histoM02TrueMergedEtaallSM[kMC-1]->Add(histoM02TrueMergedEtaperSM[iSM][kMC-1]);
                // }
                // histoM02TrueMergedEtaperSM[iSM][kMC-1]->Scale(1/histoM02TrueMergedEtaperSM[iSM][kMC-1]->Integral());

                // true part. conv pi0
                histoM02TrueMergedPartConvPi0perSMvsE[iSM][kMC-1]  = (TH2F*)inputFile[kMC]->Get(Form("histoM02TrueMergedPartConvPi0perSMvsE%d",iSM));
                projectMin            = histoM02TrueMergedPartConvPi0perSMvsE[iSM][kMC-1]->GetYaxis()->FindBin(minEPlot+0.001);
                projectMax            = histoM02TrueMergedPartConvPi0perSMvsE[iSM][kMC-1]->GetYaxis()->FindBin(maxEPlot-0.001);
                histoM02TrueMergedPartConvPi0perSM[iSM][kMC-1] = (TH1D*) histoM02TrueMergedPartConvPi0perSMvsE[iSM][kMC-1]->ProjectionX(Form("histoM02TrueMergedPartConvPi0perSM%d%d%d",iSM,kMC-1,iSM),projectMin,projectMax);
                histoM02TrueMergedPartConvPi0perSM[iSM][kMC-1]->Rebin(rebinFac);
                histoM02TrueMergedPartConvPi0perSM[iSM][kMC-1]->Sumw2();
                if(iSM==minSM){
                    histoM02TrueMergedPartConvPi0allSM[kMC-1] = (TH1F*)histoM02TrueMergedPartConvPi0perSM[iSM][kMC-1]->Clone(Form("histoM02TrueMergedPartConvPi0allSM%d",kMC-1));
                } else {
                    histoM02TrueMergedPartConvPi0allSM[kMC-1]->Add(histoM02TrueMergedPartConvPi0perSM[iSM][kMC-1]);
                }
                histoM02TrueMergedPartConvPi0perSM[iSM][kMC-1]->Scale(1/histoM02TrueMergedPartConvPi0perSM[iSM][kMC-1]->Integral());

                // true gamma from pi0
                histoM02TrueGammaFromPi0perSMvsE[iSM][kMC-1]  = (TH2F*)inputFile[kMC]->Get(Form("histoM02TrueGammaFromPi0perSMvsE%d",iSM));
                projectMin            = histoM02TrueGammaFromPi0perSMvsE[iSM][kMC-1]->GetYaxis()->FindBin(minEPlot+0.001);
                projectMax            = histoM02TrueGammaFromPi0perSMvsE[iSM][kMC-1]->GetYaxis()->FindBin(maxEPlot-0.001);
                histoM02TrueGammaFromPi0perSM[iSM][kMC-1] = (TH1D*) histoM02TrueGammaFromPi0perSMvsE[iSM][kMC-1]->ProjectionX(Form("histoM02TrueGammaFromPi0perSM%d%d%d",iSM,kMC-1,iSM),projectMin,projectMax);
                histoM02TrueGammaFromPi0perSM[iSM][kMC-1]->Rebin(rebinFac);
                histoM02TrueGammaFromPi0perSM[iSM][kMC-1]->Sumw2();
                if(iSM==minSM){
                    histoM02TrueGammaFromPi0allSM[kMC-1] = (TH1F*)histoM02TrueGammaFromPi0perSM[iSM][kMC-1]->Clone(Form("histoM02TrueGammaFromPi0allSM%d",kMC-1));
                } else {
                    histoM02TrueGammaFromPi0allSM[kMC-1]->Add(histoM02TrueGammaFromPi0perSM[iSM][kMC-1]);
                }
                histoM02TrueGammaFromPi0perSM[iSM][kMC-1]->Scale(1/histoM02TrueGammaFromPi0perSM[iSM][kMC-1]->Integral());

                // true electron from pi0
                histoM02TrueElectronFromPi0perSMvsE[iSM][kMC-1]  = (TH2F*)inputFile[kMC]->Get(Form("histoM02TrueElectronFromPi0perSMvsE%d",iSM));
                projectMin            = histoM02TrueElectronFromPi0perSMvsE[iSM][kMC-1]->GetYaxis()->FindBin(minEPlot+0.001);
                projectMax            = histoM02TrueElectronFromPi0perSMvsE[iSM][kMC-1]->GetYaxis()->FindBin(maxEPlot-0.001);
                histoM02TrueElectronFromPi0perSM[iSM][kMC-1] = (TH1D*) histoM02TrueElectronFromPi0perSMvsE[iSM][kMC-1]->ProjectionX(Form("histoM02TrueElectronFromPi0perSM%d%d%d",iSM,kMC-1,iSM),projectMin,projectMax);
                histoM02TrueElectronFromPi0perSM[iSM][kMC-1]->Rebin(rebinFac);
                histoM02TrueElectronFromPi0perSM[iSM][kMC-1]->Sumw2();
                if(iSM==minSM){
                    histoM02TrueElectronFromPi0allSM[kMC-1] = (TH1F*)histoM02TrueElectronFromPi0perSM[iSM][kMC-1]->Clone(Form("histoM02TrueElectronFromPi0allSM%d",kMC-1));
                } else {
                    histoM02TrueElectronFromPi0allSM[kMC-1]->Add(histoM02TrueElectronFromPi0perSM[iSM][kMC-1]);
                }
                histoM02TrueElectronFromPi0perSM[iSM][kMC-1]->Scale(1/histoM02TrueElectronFromPi0perSM[iSM][kMC-1]->Integral());
            }
            histoM02FracTrueMergedPi0allSM[kMC-1]           = (TH1F*)histoM02TrueMergedPi0allSM[kMC-1]->Clone(Form("histoM02FracTrueMergedPi0allSM%d",kMC-1));
            histoM02FracTrueMergedPi0allSM[kMC-1]->Divide(histoM02allSM[kMC]);
            histoM02FracTrueGammaFromPi0allSM[kMC-1]       = (TH1F*)histoM02TrueGammaFromPi0allSM[kMC-1]->Clone(Form("histoM02FracTrueGammaFromPi0allSM%d",kMC-1));
            histoM02FracTrueGammaFromPi0allSM[kMC-1]->Divide(histoM02allSM[kMC]);
            histoM02FracTrueMergedPartConvPi0allSM[kMC-1]  = (TH1F*)histoM02TrueMergedPartConvPi0allSM[kMC-1]->Clone(Form("histoM02FracTrueMergedPartConvPi0allSM%d",kMC-1));
            histoM02FracTrueMergedPartConvPi0allSM[kMC-1]->Divide(histoM02allSM[kMC]);
            histoM02FracTrueElectronFromPi0allSM[kMC-1]    = (TH1F*)histoM02TrueElectronFromPi0allSM[kMC-1]->Clone(Form("histoM02FracTrueElectronFromPi0allSM%d",kMC-1));
            histoM02FracTrueElectronFromPi0allSM[kMC-1]->Divide(histoM02allSM[kMC]);
        }
    }
    
    



    TLatex* labelpTrange;


    Int_t textSizeLabelsPixel                   = 1200*0.04;

    TCanvas* canvasM02Plots = new TCanvas("canvasM02Plots","",200,10,1350,1200);  // gives the page size
    DrawGammaCanvasSettings( canvasM02Plots, 0.08, 0.02, 0.02, 0.09);
    // canvasM02Plots->SetLogy();
    canvasM02Plots->SetLogz();

    // TH2F * histo2DM02Dummy = new TH2F("histo2DM02Dummy","histo2DM02Dummy",220,0.0, 1.25, 500,0.0008,0.06*rebinFac);
    TH2F * histo2DM02Dummy = new TH2F("histo2DM02Dummy","histo2DM02Dummy",250,0.0, 1.25, 500,-0.0001,0.06*rebinFac);
    // TH2F * histo2DM02Dummy = new TH2F("histo2DM02Dummy","histo2DM02Dummy",220,0.0, 1.5, 500,-0.1,2.0);
    SetStyleHistoTH2ForGraphs(histo2DM02Dummy, "M_{02}","norm. counts",0.035,0.04, 0.035,0.04, 0.98,0.9);
    TLegend* legendM02Dummy;
    TLatex* labelEnergy                  = new TLatex(0.94,0.90,Form("%s, %s",labelALICEforPlots.Data(),collisionSystem.Data()));
    SetStyleTLatex( labelEnergy, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    TLatex* labelSupMod;
    TLatex* labelClusterE = new TLatex(0.12,0.90,Form("%2.1f < #it{E}_{cls} < %2.1f GeV",minEPlot, maxEPlot));
    SetStyleTLatex( labelClusterE, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
    canvasM02Plots->cd();
    for(Int_t iSM=minSM; iSM<maxSM; iSM++){
        histo2DM02Dummy->Draw("copy");
        DrawGammaSetMarker(histoM02perSM[iSM][0], markerSMdata[iSM], 2.0, colorSMdata[iSM] , colorSMdata[iSM]);
        histoM02perSM[iSM][0]->Draw("p,same,e");
        DrawGammaSetMarker(histoM02perSM[iSM][1], markerSMMC[iSM], 2.0, colorSMMC[iSM] , colorSMMC[iSM]);
        histoM02perSM[iSM][1]->Draw("p,same,e");
        DrawGammaSetMarker(histoM02perSM[iSM][2], markerSMMCXT[iSM], 2.0, colorSMMCXT[iSM] , colorSMMCXT[iSM]);
        histoM02perSM[iSM][2]->Draw("p,same,e");

        legendM02Dummy  = GetAndSetLegend2(0.70, 0.84-(0.035*3*1.35), 0.98, 0.84, 0.85*textSizeLabelsPixel);
        legendM02Dummy->AddEntry(histoM02perSM[iSM][0], "data","pe");
        legendM02Dummy->AddEntry(histoM02perSM[iSM][1], labelMC1.Data(),"pe");
        legendM02Dummy->AddEntry(histoM02perSM[iSM][2], labelMC2.Data(),"pe");
        legendM02Dummy->Draw();
        labelEnergy->Draw();
        labelClusterE->Draw();
        labelSupMod = new TLatex(0.94,0.86,Form("Supermodule %d",iSM));
        SetStyleTLatex( labelSupMod, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
        labelSupMod->Draw();
        histo2DM02Dummy->Draw("axis,same");
        canvasM02Plots->SaveAs(Form("%s/M02inSM%d.%s",outputDir.Data(),iSM,suffix.Data()));
    }

    histo2DM02Dummy->Draw("copy");
    DrawGammaSetMarker(histoM02allSM[0], markerSMdata[0], 2.0, colorSMdata[0] , colorSMdata[0]);
    histoM02allSM[0]->Scale(1/histoM02allSM[0]->Integral());
    histoM02allSM[0]->Draw("p,same,e");
    DrawGammaSetMarker(histoM02allSM[1], markerSMMC[0], 2.0, colorSMMC[0] , colorSMMC[0]);
    DrawGammaSetMarker(histoM02allSMPlot[1], markerSMMC[0], 2.0, colorSMMC[0] , colorSMMC[0]);
    Double_t MCscalingForTrue[2] = {histoM02allSM[1]->Integral(),histoM02allSM[2]->Integral()};
    histoM02allSM[1]->Scale(1/MCscalingForTrue[0]);
    histoM02allSMPlot[1]->Scale(1/MCscalingForTrue[0]);
    histoM02allSM[1]->Draw("p,same,e");
    DrawGammaSetMarker(histoM02allSM[2], markerSMMCXT[0], 2.0, colorSMMCXT[0] , colorSMMCXT[0]);
    histoM02allSM[2]->Scale(1/MCscalingForTrue[0]);
    histoM02allSMPlot[2]->Scale(1/MCscalingForTrue[1]);
    histoM02allSM[2]->Draw("p,same,e");
    legendM02Dummy  = GetAndSetLegend2(0.70, 0.84-(0.035*3*1.35), 0.98, 0.84, 0.85*textSizeLabelsPixel);
    legendM02Dummy->AddEntry(histoM02allSM[0], "data","pe");
    legendM02Dummy->AddEntry(histoM02allSM[1], labelMC1.Data(),"pe");
    legendM02Dummy->AddEntry(histoM02allSM[2], labelMC2.Data(),"pe");
    legendM02Dummy->Draw();
    labelEnergy->Draw();
    labelSupMod = new TLatex(0.94,0.86,Form("Supermodules %d-%d",minSM,maxSM-1));
    SetStyleTLatex( labelSupMod, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelSupMod->Draw();
    labelClusterE->Draw();
    histo2DM02Dummy->Draw("axis,same");
    canvasM02Plots->SaveAs(Form("%s/M02allSM.%s",outputDir.Data(),suffix.Data()));

    histo2DM02Dummy->Draw("copy");
    histoM02allSM[0]->Draw("p,same,e");
    histoM02allSM[1]->Draw("p,same,e");
    histoM02allSM[2]->Draw("p,same,e");
    legendM02Dummy->Draw();
    labelEnergy->Draw();
    Double_t integratedCounts[2] = {histoM02allSM[1]->Integral(20/rebinFac,60/rebinFac),histoM02allSM[2]->Integral(20/rebinFac,60/rebinFac)};
    // Double_t integratedCounts[2] = {histoM02allSM[1]->Integral(54/rebinFac,300/rebinFac),histoM02allSM[2]->Integral(54/rebinFac,300/rebinFac)};
    // histoM02allSM[1]->SetBinContent(20/rebinFac,0.1);
    // histoM02allSM[1]->SetBinContent(60/rebinFac,0.1);
    DrawGammaLines( 0.1, 0.1, 0.0, 0.1, 3, kGray+1, 9);
    DrawGammaLines( 0.3, 0.3, 0.0, 0.1, 3, kGray+1, 9);
    TLatex* labelIntegralMC1 = new TLatex(0.94,0.60,Form("%s: #sum_{M_{0}^{2} = 0.1}^{0.3}=%1.3f",labelMC1.Data(), integratedCounts[0]));
    SetStyleTLatex( labelIntegralMC1, 0.85*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelIntegralMC1->Draw();
    TLatex* labelIntegralMC2 = new TLatex(0.94,0.50,Form("%s: #sum_{M_{0}^{2} = 0.1}^{0.3}=%1.3f",labelMC2.Data(), integratedCounts[1]));
    SetStyleTLatex( labelIntegralMC2, 0.85*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelIntegralMC2->Draw();
    // TLatex* labelIntegralMC3 = new TLatex(0.94,0.40,Form("Ratio of MCs: %1.3f", integratedCounts[0]/integratedCounts[1]));
    TLatex* labelIntegralMC3 = new TLatex(0.94,0.40,Form("Ratio of sums: %1.2f %%", (1-(integratedCounts[0]/integratedCounts[1]))*100));
    SetStyleTLatex( labelIntegralMC3, 0.85*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelIntegralMC3->Draw();
    labelSupMod->Draw();
    labelClusterE->Draw();
    histo2DM02Dummy->Draw("axis,same");
    canvasM02Plots->SaveAs(Form("%s/M02allSM_withLowIntegral.%s",outputDir.Data(),suffix.Data()));

    histo2DM02Dummy->Draw("copy");
    histoM02allSM[0]->Draw("p,same,e");
    histoM02allSM[1]->Draw("p,same,e");
    histoM02allSM[2]->Draw("p,same,e");
    legendM02Dummy->Draw();
    labelEnergy->Draw();
    Double_t integratedCounts2[2] = {histoM02allSM[1]->Integral(54/rebinFac,240/rebinFac),histoM02allSM[2]->Integral(54/rebinFac,240/rebinFac)};
    // histoM02allSM[1]->SetBinContent(54/rebinFac,0.1);
    // histoM02allSM[1]->SetBinContent(240/rebinFac,0.1);
    DrawGammaLines( 0.27, 0.27, 0.0, 0.1, 3, kGray+1, 9);
    DrawGammaLines( 1.2, 1.2, 0.0, 0.03, 3, kGray+1, 9);
    labelIntegralMC1 = new TLatex(0.94,0.60,Form("%s: #sum_{M_{0}^{2} = 0.27}^{1.2}=%1.3f",labelMC1.Data(), integratedCounts2[0]));
    SetStyleTLatex( labelIntegralMC1, 0.85*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelIntegralMC1->Draw();
    labelIntegralMC2 = new TLatex(0.94,0.50,Form("%s: #sum_{M_{0}^{2} = 0.27}^{1.2}=%1.3f",labelMC2.Data(), integratedCounts2[1]));
    SetStyleTLatex( labelIntegralMC2, 0.85*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelIntegralMC2->Draw();
    labelIntegralMC3 = new TLatex(0.94,0.40,Form("Ratio of sums: %1.2f %%", (1-(integratedCounts2[0]/integratedCounts2[1]))*100));
    SetStyleTLatex( labelIntegralMC3, 0.85*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelIntegralMC3->Draw();
    labelSupMod->Draw();
    labelClusterE->Draw();
    histo2DM02Dummy->Draw("axis,same");
    canvasM02Plots->SaveAs(Form("%s/M02allSM_withHighIntegral.%s",outputDir.Data(),suffix.Data()));


    DrawGammaSetMarker(histoM02allSM[1], markerSMMC[0], 2.0, kGray+1 , kGray+1);
    DrawGammaSetMarker(histoM02allSM[2], markerSMMCXT[0], 2.0, kGray+3 , kGray+3);
    // DrawGammaSetMarker(histoM02allSM[1], markerSMMC[0], 2.0, colorSMMCXT[0] , colorSMMCXT[0]);
    // DrawGammaSetMarker(histoM02allSM[2], markerSMMCXT[0], 2.0, colorSMMCXT[0] , colorSMMCXT[0]);

    histo2DM02Dummy->Draw("copy");
    histoM02allSM[0]->Draw("p,same,e");
    histoM02allSM[1]->Draw("p,same,e");
    DrawGammaSetMarker(histoM02TrueMergedPi0allSM[0], markerSMMC[4], 2.0, colorSMMC[4] , colorSMMC[4]);
    histoM02TrueMergedPi0allSM[0]->Scale(1/MCscalingForTrue[0]);
    histoM02TrueMergedPi0allSM[0]->Draw("p,same,e");
    DrawGammaSetMarker(histoM02TrueGammaFromPi0allSM[0], markerSMMC[1], 2.0, colorSMMC[1] , colorSMMC[1]);
    histoM02TrueGammaFromPi0allSM[0]->Scale(1/MCscalingForTrue[0]);
    histoM02TrueGammaFromPi0allSM[0]->Draw("p,same,e");
    DrawGammaSetMarker(histoM02TrueMergedPartConvPi0allSM[0], markerSMMC[2], 2.0, colorSMMC[2] , colorSMMC[2]);
    histoM02TrueMergedPartConvPi0allSM[0]->Scale(1/MCscalingForTrue[0]);
    histoM02TrueMergedPartConvPi0allSM[0]->Draw("p,same,e");
    DrawGammaSetMarker(histoM02TrueElectronFromPi0allSM[0], markerSMMC[3], 2.0, colorSMMC[3] , colorSMMC[3]);
    histoM02TrueElectronFromPi0allSM[0]->Scale(1/MCscalingForTrue[0]);
    histoM02TrueElectronFromPi0allSM[0]->Draw("p,same,e");
    // DrawGammaSetMarker(histoM02TrueMergedEtaallSM[0], markerSMMC[4], 2.0, colorSMMC[4] , colorSMMC[4]);
    // histoM02TrueMergedEtaallSM[0]->Scale(1/MCscalingForTrue[0]);
    // histoM02TrueMergedEtaallSM[0]->Draw("p,same,e");
    legendM02Dummy  = GetAndSetLegend2(0.60, 0.84-(0.035*5*1.35), 0.88, 0.84, 0.85*textSizeLabelsPixel);
    legendM02Dummy->AddEntry(histoM02allSM[0], "data","pe");
    legendM02Dummy->AddEntry(histoM02allSM[1], labelMC1.Data(),"pe");
    legendM02Dummy->AddEntry(histoM02TrueMergedPi0allSM[0], "merged #pi^{0}","pe");
    // legendM02Dummy->AddEntry(histoM02TrueMergedEtaallSM[0], "merged #eta","pe");
    legendM02Dummy->AddEntry(histoM02TrueMergedPartConvPi0allSM[0], "merged #pi^{0} part. conv.","pe");
    legendM02Dummy->AddEntry(histoM02TrueGammaFromPi0allSM[0], "#gamma from #pi^{0}","pe");
    legendM02Dummy->AddEntry(histoM02TrueElectronFromPi0allSM[0], "e^{#pm} from #pi^{0}","pe");
    legendM02Dummy->Draw();
    labelEnergy->Draw();
    labelSupMod = new TLatex(0.94,0.86,Form("Supermodules %d-%d",minSM,maxSM-1));
    SetStyleTLatex( labelSupMod, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelSupMod->Draw();
    labelClusterE->Draw();
    histo2DM02Dummy->Draw("axis,same");
    canvasM02Plots->SaveAs(Form("%s/M02allSMTrueDecomp_standardMC.%s",outputDir.Data(),suffix.Data()));

    histo2DM02Dummy->Draw("copy");
    histoM02allSM[0]->Draw("p,same,e");
    histoM02allSM[2]->Draw("p,same,e");
    DrawGammaSetMarker(histoM02TrueMergedPi0allSM[1], markerSMMC[4], 2.0, colorSMMCXT[4] , colorSMMCXT[4]);
    histoM02TrueMergedPi0allSM[1]->Scale(1/MCscalingForTrue[1]);
    histoM02TrueMergedPi0allSM[1]->Draw("p,same,e");
    DrawGammaSetMarker(histoM02TrueGammaFromPi0allSM[1], markerSMMC[1], 2.0, colorSMMCXT[1] , colorSMMCXT[1]);
    histoM02TrueGammaFromPi0allSM[1]->Scale(1/MCscalingForTrue[1]);
    histoM02TrueGammaFromPi0allSM[1]->Draw("p,same,e");
    DrawGammaSetMarker(histoM02TrueMergedPartConvPi0allSM[1], markerSMMC[2], 2.0, colorSMMCXT[2] , colorSMMCXT[2]);
    histoM02TrueMergedPartConvPi0allSM[1]->Scale(1/MCscalingForTrue[1]);
    histoM02TrueMergedPartConvPi0allSM[1]->Draw("p,same,e");
    DrawGammaSetMarker(histoM02TrueElectronFromPi0allSM[1], markerSMMC[3], 2.0, colorSMMCXT[3] , colorSMMCXT[3]);
    histoM02TrueElectronFromPi0allSM[1]->Scale(1/MCscalingForTrue[1]);
    histoM02TrueElectronFromPi0allSM[1]->Draw("p,same,e");
    // DrawGammaSetMarker(histoM02TrueMergedEtaallSM[1], markerSMMC[4], 2.0, colorSMMCXT[4] , colorSMMCXT[4]);
    // histoM02TrueMergedEtaallSM[1]->Scale(1/MCscalingForTrue[1]);
    // histoM02TrueMergedEtaallSM[1]->Draw("p,same,e");
    legendM02Dummy  = GetAndSetLegend2(0.60, 0.84-(0.035*5*1.35), 0.88, 0.84, 0.85*textSizeLabelsPixel);
    legendM02Dummy->AddEntry(histoM02allSM[0], "data","pe");
    legendM02Dummy->AddEntry(histoM02allSM[2], labelMC2.Data(),"pe");
    legendM02Dummy->AddEntry(histoM02TrueMergedPi0allSM[1], "merged #pi^{0}","pe");
    // legendM02Dummy->AddEntry(histoM02TrueMergedEtaallSM[1], "merged #eta","pe");
    legendM02Dummy->AddEntry(histoM02TrueMergedPartConvPi0allSM[1], "merged #pi^{0} part. conv.","pe");
    legendM02Dummy->AddEntry(histoM02TrueGammaFromPi0allSM[1], "#gamma from #pi^{0}","pe");
    legendM02Dummy->AddEntry(histoM02TrueElectronFromPi0allSM[1], "e^{#pm} from #pi^{0}","pe");
    legendM02Dummy->Draw();
    labelEnergy->Draw();
    labelSupMod = new TLatex(0.94,0.86,Form("Supermodules %d-%d",minSM,maxSM-1));
    SetStyleTLatex( labelSupMod, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelSupMod->Draw();
    labelClusterE->Draw();
    histo2DM02Dummy->Draw("axis,same");
    canvasM02Plots->SaveAs(Form("%s/M02allSMTrueDecomp_xtalkMC.%s",outputDir.Data(),suffix.Data()));

    histo2DM02Dummy->Draw("copy");
    histoM02allSM[0]->Draw("p,same,e");
    histoM02allSM[1]->Draw("p,same,e");
    histoM02allSM[2]->Draw("p,same,e");
    histoM02TrueMergedPi0allSM[0]->Draw("p,same,e");
    histoM02TrueGammaFromPi0allSM[0]->Draw("p,same,e");
    histoM02TrueMergedPartConvPi0allSM[0]->Draw("p,same,e");
    histoM02TrueElectronFromPi0allSM[0]->Draw("p,same,e");
    // histoM02TrueMergedEtaallSM[0]->Draw("p,same,e");
    histoM02TrueMergedPi0allSM[1]->Draw("p,same,e");
    histoM02TrueGammaFromPi0allSM[1]->Draw("p,same,e");
    histoM02TrueMergedPartConvPi0allSM[1]->Draw("p,same,e");
    histoM02TrueElectronFromPi0allSM[1]->Draw("p,same,e");
    // histoM02TrueMergedEtaallSM[1]->Draw("p,same,e");
    legendM02Dummy  = GetAndSetLegend2(0.47, 0.84-(0.035*5*1.35), 0.95, 0.84, 0.85*textSizeLabelsPixel,2);
    legendM02Dummy->AddEntry(histoM02allSM[0], "data","pe");
    legendM02Dummy->AddEntry((TObject*)0, "","");
    legendM02Dummy->AddEntry(histoM02allSM[1], labelMC1.Data(),"pe");
    legendM02Dummy->AddEntry(histoM02allSM[2], labelMC2.Data(),"pe");
    legendM02Dummy->AddEntry(histoM02TrueMergedPi0allSM[0], "merged #pi^{0}","pe");
    legendM02Dummy->AddEntry(histoM02TrueMergedPi0allSM[1], "merged #pi^{0}","pe");
    // legendM02Dummy->AddEntry(histoM02TrueMergedEtaallSM[1], "merged #eta","pe");
    legendM02Dummy->AddEntry(histoM02TrueMergedPartConvPi0allSM[0], "merged #pi^{0} p.c.","pe");
    legendM02Dummy->AddEntry(histoM02TrueMergedPartConvPi0allSM[1], "merged #pi^{0} p.c.","pe");
    legendM02Dummy->AddEntry(histoM02TrueGammaFromPi0allSM[0], "#gamma from #pi^{0}","pe");
    legendM02Dummy->AddEntry(histoM02TrueGammaFromPi0allSM[1], "#gamma from #pi^{0}","pe");
    legendM02Dummy->AddEntry(histoM02TrueElectronFromPi0allSM[0], "e^{#pm} from #pi^{0}","pe");
    legendM02Dummy->AddEntry(histoM02TrueElectronFromPi0allSM[1], "e^{#pm} from #pi^{0}","pe");
    legendM02Dummy->Draw();
    labelEnergy->Draw();
    labelSupMod = new TLatex(0.94,0.86,Form("Supermodules %d-%d",minSM,maxSM-1));
    SetStyleTLatex( labelSupMod, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelSupMod->Draw();
    labelClusterE->Draw();
    histo2DM02Dummy->Draw("axis,same");
    canvasM02Plots->SaveAs(Form("%s/M02allSMTrueDecomp_bothMC.%s",outputDir.Data(),suffix.Data()));



    DrawGammaSetMarker(histoM02allSM[1], markerSMMC[0], 2.0, colorSMMC[0] , colorSMMC[0]);
    DrawGammaSetMarker(histoM02allSM[2], markerSMMCXT[0], 2.0, colorSMMCXT[0] , colorSMMCXT[0]);

    TCanvas* canvasM02RatioPlots = new TCanvas("canvasM02RatioPlots","",200,10,1350,1200);  // gives the page size
    DrawGammaCanvasSettings( canvasM02RatioPlots, 0.08, 0.02, 0.02, 0.09);
    // canvasM02RatioPlots->SetLogy();
    canvasM02RatioPlots->SetLogz();

    TH2F * histo2DM02RatioDummy = new TH2F("histo2DM02RatioDummy","histo2DM02RatioDummy",250,0.0, 1.25, 500,-0.1,3.0);
    SetStyleHistoTH2ForGraphs(histo2DM02RatioDummy, "M_{02}","ratio",0.035,0.04, 0.035,0.04, 0.98,0.9);
    // histo2DM02RatioDummy->GetZaxis()->SetRangeUser(6,6e2);
    labelEnergy                  = new TLatex(0.25,0.90,Form("%s, %s",labelALICEforPlots.Data(),collisionSystem.Data()));
    SetStyleTLatex( labelEnergy, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
    labelClusterE = new TLatex(0.12,0.13,Form("%2.1f < #it{E}_{cls} < %2.1f GeV",minEPlot, maxEPlot));
    SetStyleTLatex( labelClusterE, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
    canvasM02RatioPlots->cd();
    for(Int_t iSM=minSM; iSM<maxSM; iSM++){
        histo2DM02RatioDummy->Draw("copy");
        DrawGammaSetMarker(histoM02RatioperSM[iSM][0], markerSMMC[iSM], 2.0, colorSMMC[iSM] , colorSMMC[iSM]);
        histoM02RatioperSM[iSM][0]->Draw("p,same,e");
        DrawGammaSetMarker(histoM02RatioperSM[iSM][1], markerSMMCXT[iSM], 2.0, colorSMMCXT[iSM] , colorSMMCXT[iSM]);
        histoM02RatioperSM[iSM][1]->Draw("p,same,e");
        labelEnergy->Draw();
        labelSupMod = new TLatex(0.25,0.86,Form("Supermodule %d",iSM));
        SetStyleTLatex( labelSupMod, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
        labelSupMod->Draw();
        labelClusterE->Draw();
        legendM02Dummy  = GetAndSetLegend2(0.24, 0.85-(0.035*2*1.35), 0.52, 0.85, 0.85*textSizeLabelsPixel);
        legendM02Dummy->AddEntry(histoM02RatioperSM[iSM][0], Form("data / %s", labelMC1.Data()),"pe");
        legendM02Dummy->AddEntry(histoM02RatioperSM[iSM][1], Form("data / %s", labelMC2.Data()),"pe");
        legendM02Dummy->Draw();
        histo2DM02RatioDummy->Draw("axis,same");
        canvasM02RatioPlots->SaveAs(Form("%s/Ratios/M02RatioinSM%d.%s",outputDir.Data(),iSM,suffix.Data()));
    }

    // histo2DM02RatioDummy->GetYaxis()->SetRangeUser(0.,2.0);
    histo2DM02RatioDummy->Draw("copy");
    histoM02RatioallSM[0] = (TH1F*)histoM02allSM[0]->Clone(Form("histoM02RatioallSM_MC%d",0));
    histoM02RatioallSM[0]->Sumw2();
    histoM02RatioallSM[0]->Divide(histoM02allSM[1]);
    // histoM02RatioallSM[0]->Divide(histoM02RatioallSM[0],histoM02allSM[1],1,1,"B");
    DrawGammaSetMarker(histoM02RatioallSM[0], markerSMMC[0], 2.0, colorSMMC[0] , colorSMMC[0]);
    histoM02RatioallSM[1] = (TH1F*)histoM02allSM[0]->Clone(Form("histoM02RatioallSM_MC%d",0));
    histoM02RatioallSM[1]->Sumw2();
    histoM02RatioallSM[1]->Divide(histoM02allSM[2]);
    // histoM02RatioallSM[1]->Divide(histoM02RatioallSM[1],histoM02allSM[2],1,1,"B");
    DrawGammaSetMarker(histoM02RatioallSM[1], markerSMMCXT[0], 2.0, colorSMMCXT[0] , colorSMMCXT[0]);
    labelEnergy->Draw();
    labelSupMod = new TLatex(0.25,0.86,Form("Supermodules %d-%d",minSM,maxSM-1));
    SetStyleTLatex( labelSupMod, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
    labelSupMod->Draw();
    labelClusterE->Draw();
    legendM02Dummy  = GetAndSetLegend2(0.24, 0.85-(0.035*2*1.35), 0.52, 0.85, 0.85*textSizeLabelsPixel);
    legendM02Dummy->AddEntry(histoM02RatioallSM[0], Form("data / %s", labelMC1.Data()),"pe");
    legendM02Dummy->AddEntry(histoM02RatioallSM[1], Form("data / %s", labelMC2.Data()),"pe");
    legendM02Dummy->Draw();
    histoM02RatioallSM[0]->Draw("p,same,e");
    histoM02RatioallSM[1]->Draw("p,same,e");
    histo2DM02RatioDummy->Draw("axis,same");
    canvasM02RatioPlots->SaveAs(Form("%s/Ratios/M02RatioinallSM.%s",outputDir.Data(),suffix.Data()));

    // histo2DM02RatioDummy->GetYaxis()->SetRangeUser(0.,2.0);
    histo2DM02RatioDummy->Draw("copy");
    labelEnergy->Draw();
    histoM02RatioallSMMC = (TH1F*)histoM02allSM[1]->Clone(Form("histoM02RatioallSMMC%d",0));
    histoM02RatioallSMMC->Sumw2();
    // histoM02RatioallSMMC->Divide(histoM02allSM[2]);
    histoM02RatioallSMMC->Divide(histoM02RatioallSMMC,histoM02allSM[2],1,1,"B");
    DrawGammaSetMarker(histoM02RatioallSMMC, markerSMdata[0], 2.0, kBlue+2 , kBlue+2);
    labelEnergy->Draw();
    labelSupMod->Draw();
    labelClusterE->Draw();
    legendM02Dummy  = GetAndSetLegend2(0.24, 0.85-(0.035*3*1.35), 0.52, 0.85, 0.85*textSizeLabelsPixel);
    legendM02Dummy->AddEntry(histoM02RatioallSM[0], Form("data / %s", labelMC1.Data()),"pe");
    legendM02Dummy->AddEntry(histoM02RatioallSM[1], Form("data / %s", labelMC2.Data()),"pe");
    legendM02Dummy->AddEntry(histoM02RatioallSMMC, Form("%s / %s", labelMC1.Data(),labelMC2.Data()),"pe");
    legendM02Dummy->Draw();
    histoM02RatioallSM[0]->Draw("p,same,e");
    histoM02RatioallSM[1]->Draw("p,same,e");
    histoM02RatioallSMMC->Draw("p,same,e");
    histo2DM02RatioDummy->Draw("axis,same");
    canvasM02RatioPlots->SaveAs(Form("%s/Ratios/M02RatioinallSM_withMCratio.%s",outputDir.Data(),suffix.Data()));

    TF1* pol0[20][2];
    // histo2DM02RatioDummy->GetYaxis()->SetRangeUser(0.,3.0);
    labelClusterE = new TLatex(0.30,0.13,Form("%2.1f < #it{E}_{cls} < %2.1f GeV",minEPlot, maxEPlot));
    SetStyleTLatex( labelClusterE, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
    for(Int_t iSM=minSM; iSM<maxSM; iSM++){
        histo2DM02RatioDummy->Draw("copy");
        histoM02RatioperSM[iSM][0]->Draw("p,same,e");
        histoM02RatioperSM[iSM][1]->Draw("p,same,e");
        DrawGammaLines( 0.1, 0.1, -0.1, 2.0, 3, kGreen+3, 9);
        DrawGammaLines( 0.27, 0.27, -0.1, 2.0, 3, kRed+3, 9);
        DrawGammaLines( 0.0, 1.45, 1.0, 1.0, 2, kGray+1, 2);
        pol0[iSM][0] = new TF1("pol0","[0]",0.27,1.5);//
        histoM02RatioperSM[iSM][0]->Fit(pol0[iSM][0],"NRMEX0+","",0.27,1.5);
        pol0[iSM][0]->SetLineColor(colorSMMC[iSM]+1);
        TBox* boxPol0FitError = CreateBoxConv(colorSMMC[iSM], 0.27,  pol0[iSM][0]->GetParameter(0)-pol0[iSM][0]->GetParError(0), 1.5, pol0[iSM][0]->GetParameter(0)+pol0[iSM][0]->GetParError(0));
        boxPol0FitError->SetLineWidth(8);
        boxPol0FitError->Draw();
        pol0[iSM][0]->Draw("same");
        pol0[iSM][1] = new TF1("pol0","[0]",0.27,1.5);//
        histoM02RatioperSM[iSM][1]->Fit(pol0[iSM][1],"NRMEX0+","",0.27,1.5);
        pol0[iSM][1]->SetLineColor(colorSMMCXT[iSM]+1);
        TBox* boxPol0FitError2 = CreateBoxConv(colorSMMCXT[iSM], 0.27,  pol0[iSM][1]->GetParameter(0)-pol0[iSM][1]->GetParError(0), 1.5, pol0[iSM][1]->GetParameter(0)+pol0[iSM][1]->GetParError(0));
        boxPol0FitError2->SetLineWidth(8);
        boxPol0FitError2->Draw();
        pol0[iSM][1]->Draw("same");
        histoM02RatioperSM[iSM][0]->Draw("p,same,e");
        histoM02RatioperSM[iSM][1]->Draw("p,same,e");
        labelEnergy->Draw();
        labelClusterE->Draw();
        labelSupMod = new TLatex(0.25,0.86,Form("Supermodule %d",iSM));
        SetStyleTLatex( labelSupMod, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
        labelSupMod->Draw();
        legendM02Dummy  = GetAndSetLegend2(0.24, 0.85-(0.035*2*1.35), 0.52, 0.85, 0.85*textSizeLabelsPixel);
        legendM02Dummy->AddEntry(histoM02RatioperSM[iSM][0], labelMC1.Data(),"pe");
        legendM02Dummy->AddEntry(histoM02RatioperSM[iSM][1], labelMC2.Data(),"pe");
        legendM02Dummy->Draw();
        histo2DM02RatioDummy->Draw("axis,same");
        canvasM02RatioPlots->SaveAs(Form("%s/Ratios/withFits/fittedM02RatioinSM%d.%s",outputDir.Data(),iSM,suffix.Data()));
    }

    // histo2DM02RatioDummy->GetYaxis()->SetRangeUser(0.,2.0);
    histo2DM02RatioDummy->Draw("copy");
    histoM02RatioallSM[0]->Draw("p,same,e");
    histoM02RatioallSM[1]->Draw("p,same,e");
    TF1* pol0All = new TF1("pol0","[0]",0.27,1.5);//
    histoM02RatioallSM[0]->Fit(pol0All,"NRMEX0+","",0.27,1.5);
    pol0All->SetLineColor(kGray+2);
    DrawGammaLines( 0.1, 0.1, -0.1, 2.0, 3, kGreen+3, 9);
    DrawGammaLines( 0.27, 0.27, -0.1, 2.0, 3, kRed+3, 9);
    DrawGammaLines( 0.0, 1.45, 1.0, 1.0, 2, kGray+1, 2);
    TBox* boxPol0FitError = CreateBoxConv(kGray+1, 0.27,  pol0All->GetParameter(0)-pol0All->GetParError(0), 1.5, pol0All->GetParameter(0)+pol0All->GetParError(0));
    boxPol0FitError->SetLineWidth(8);
    boxPol0FitError->Draw();
    histoM02RatioallSM[0]->Draw("p,same,e");
    pol0All->Draw("same");
    labelEnergy->Draw();
    labelClusterE->Draw();
    legendM02Dummy  = GetAndSetLegend2(0.24, 0.85-(0.035*2*1.35), 0.52, 0.85, 0.85*textSizeLabelsPixel);
    legendM02Dummy->AddEntry(histoM02RatioallSM[0], Form("data / %s", labelMC1.Data()),"pe");
    legendM02Dummy->AddEntry(histoM02RatioallSM[1], Form("data / %s", labelMC2.Data()),"pe");
    legendM02Dummy->Draw();
    labelSupMod = new TLatex(0.25,0.86,Form("Supermodules %d-%d",minSM,maxSM-1));
    SetStyleTLatex( labelSupMod, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
    labelSupMod->Draw();
    histo2DM02RatioDummy->Draw("axis,same");
    canvasM02RatioPlots->SaveAs(Form("%s/Ratios/withFits/fittedM02RatioinallSM.%s",outputDir.Data(),suffix.Data()));

    // histo2DM02RatioDummy->GetYaxis()->SetRangeUser(0.,2.0);
    histo2DM02RatioDummy->Draw("copy");
    histoM02RatioallSM[0]->Draw("p,same,e");
    histoM02RatioallSM[1]->Draw("p,same,e");
    DrawGammaLines( 0.1, 0.1, -0.1, 2.0, 3, kGreen+3, 9);
    DrawGammaLines( 0.27, 0.27, -0.1, 2.0, 3, kRed+3, 9);
    DrawGammaLines( 0.0, 1.45, 1.0, 1.0, 2, kGray+1, 2);
    histoM02RatioallSMMC->Draw("p,same,e");
    labelEnergy->Draw();
    labelClusterE->Draw();
    legendM02Dummy  = GetAndSetLegend2(0.24, 0.85-(0.035*3*1.35), 0.52, 0.85, 0.85*textSizeLabelsPixel);
    legendM02Dummy->AddEntry(histoM02RatioallSM[0], Form("data / %s", labelMC1.Data()),"pe");
    legendM02Dummy->AddEntry(histoM02RatioallSM[1], Form("data / %s", labelMC2.Data()),"pe");
    legendM02Dummy->AddEntry(histoM02RatioallSMMC, Form("%s / %s", labelMC1.Data(),labelMC2.Data()),"pe");
    legendM02Dummy->Draw();
    labelSupMod->Draw();
    histo2DM02RatioDummy->Draw("axis,same");
    canvasM02RatioPlots->SaveAs(Form("%s/Ratios/withFits/fittedM02RatioinallSMWithMCratio.%s",outputDir.Data(),suffix.Data()));


    histo2DM02RatioDummy->GetYaxis()->SetRangeUser(0.5,2.5);
    histo2DM02RatioDummy->Draw("copy");
    DrawGammaLines( 0.1, 0.1, 0.5, 2.0, 3, kGreen+3, 9);
    DrawGammaLines( 0.27, 0.27, 0.5, 2.0, 3, kRed+3, 9);
    // DrawGammaLines( 0.1, 0.1, -0.1, 2.0, 3, kGreen+3, 9);
    // DrawGammaLines( 0.27, 0.27, -0.1, 2.0, 3, kRed+3, 9);
    DrawGammaLines( 0.0, 1.45, 1.0, 1.0, 2, kGray+1, 2);
    histoM02RatioallSMMC->Draw("p,same,e");
    labelEnergy->Draw();
    labelClusterE->Draw();
    legendM02Dummy  = GetAndSetLegend2(0.24, 0.85-(0.035*1*1.35), 0.52, 0.85, 0.85*textSizeLabelsPixel);
    legendM02Dummy->AddEntry(histoM02RatioallSMMC, Form("%s / %s", labelMC1.Data(),labelMC2.Data()),"pe");
    legendM02Dummy->Draw();
    labelSupMod->Draw();
    histo2DM02RatioDummy->Draw("axis,same");
    canvasM02RatioPlots->SaveAs(Form("%s/Ratios/withFits/M02RatioinallSMOnlyMCratio.%s",outputDir.Data(),suffix.Data()));







    histo2DM02RatioDummy->GetYaxis()->SetRangeUser(-0.1,0.7);
    histo2DM02RatioDummy->Draw("copy");


    DrawGammaSetMarker(histoM02FracTrueMergedPi0allSM[0], markerSMMC[4], 2.0, colorTruePart[4] , colorTruePart[4]);
    DrawGammaSetMarker(histoM02FracTrueGammaFromPi0allSM[0], markerSMMC[1], 2.0, colorTruePart[1] , colorTruePart[1]);
    DrawGammaSetMarker(histoM02FracTrueMergedPartConvPi0allSM[0], markerSMMC[2], 2.0, colorTruePart[2] , colorTruePart[2]);
    DrawGammaSetMarker(histoM02FracTrueElectronFromPi0allSM[0], markerSMMC[3], 2.0, colorTruePart[3] , colorTruePart[3]);


    histoM02FracTrueMergedPi0allSM[0]->Draw("p,same,e");
    histoM02FracTrueGammaFromPi0allSM[0]->Draw("p,same,e");
    histoM02FracTrueMergedPartConvPi0allSM[0]->Draw("p,same,e");
    histoM02FracTrueElectronFromPi0allSM[0]->Draw("p,same,e");
    // histoM02RatioallSM[1]->Draw("p,same,e");
    // TF1* pol0All = new TF1("pol0","[0]",0.27,1.5);//
    // histoM02RatioallSM[0]->Fit(pol0All,"NRMEX0+","",0.27,1.5);
    // pol0All->SetLineColor(kGray+2);
    DrawGammaLines( 0.1, 0.1, -0.1, 2.0, 3, kGreen+3, 9);
    DrawGammaLines( 0.27, 0.27, -0.1, 2.0, 3, kRed+3, 9);
    // TBox* boxPol0FitError = CreateBoxConv(kGray+1, 0.27,  pol0All->GetParameter(0)-pol0All->GetParError(0), 1.5, pol0All->GetParameter(0)+pol0All->GetParError(0));
    // boxPol0FitError->SetLineWidth(8);
    // boxPol0FitError->Draw();
    // histoM02RatioallSM[0]->Draw("p,same,e");
    // pol0All->Draw("same");
    labelEnergy->Draw();
    labelClusterE->Draw();
    // legendM02Dummy  = GetAndSetLegend2(0.24, 0.85-(0.035*2*1.35), 0.52, 0.85, 0.85*textSizeLabelsPixel);
    // legendM02Dummy->AddEntry(histoM02RatioallSM[0], labelMC1.Data(),"pe");
    // legendM02Dummy->AddEntry(histoM02RatioallSM[1], labelMC2.Data(),"pe");
    // legendM02Dummy->Draw();
    labelSupMod = new TLatex(0.25,0.86,Form("Supermodules %d-%d",minSM,maxSM-1));
    SetStyleTLatex( labelSupMod, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
    labelSupMod->Draw();
    histo2DM02RatioDummy->Draw("axis,same");
    canvasM02RatioPlots->SaveAs(Form("%s/Ratios/MCtrueDecomp_stdMC.%s",outputDir.Data(),suffix.Data()));



    histo2DM02RatioDummy->Draw("copy");
    DrawGammaSetMarker(histoM02FracTrueMergedPi0allSM[1], markerSMMC[4], 2.0, colorTruePart2[4] , colorTruePart2[4]);
    DrawGammaSetMarker(histoM02FracTrueGammaFromPi0allSM[1], markerSMMC[1], 2.0, colorTruePart2[1] , colorTruePart2[1]);
    DrawGammaSetMarker(histoM02FracTrueMergedPartConvPi0allSM[1], markerSMMC[2], 2.0, colorTruePart2[2] , colorTruePart2[2]);
    DrawGammaSetMarker(histoM02FracTrueElectronFromPi0allSM[1], markerSMMC[3], 2.0, colorTruePart2[3] , colorTruePart2[3]);


    histoM02FracTrueMergedPi0allSM[1]->Draw("p,same,e");
    histoM02FracTrueGammaFromPi0allSM[1]->Draw("p,same,e");
    histoM02FracTrueMergedPartConvPi0allSM[1]->Draw("p,same,e");
    histoM02FracTrueElectronFromPi0allSM[1]->Draw("p,same,e");
    // histoM02RatioallSM[1]->Draw("p,same,e");
    // TF1* pol0All = new TF1("pol0","[0]",0.27,1.5);//
    // histoM02RatioallSM[0]->Fit(pol0All,"NRMEX0+","",0.27,1.5);
    // pol0All->SetLineColor(kGray+2);
    DrawGammaLines( 0.1, 0.1, -0.1, 2.0, 3, kGreen+3, 9);
    DrawGammaLines( 0.27, 0.27, -0.1, 2.0, 3, kRed+3, 9);
    // TBox* boxPol0FitError = CreateBoxConv(kGray+1, 0.27,  pol0All->GetParameter(0)-pol0All->GetParError(0), 1.5, pol0All->GetParameter(0)+pol0All->GetParError(0));
    // boxPol0FitError->SetLineWidth(8);
    // boxPol0FitError->Draw();
    // histoM02RatioallSM[0]->Draw("p,same,e");
    // pol0All->Draw("same");
    labelEnergy->Draw();
    labelClusterE->Draw();
    // legendM02Dummy  = GetAndSetLegend2(0.24, 0.85-(0.035*2*1.35), 0.52, 0.85, 0.85*textSizeLabelsPixel);
    // legendM02Dummy->AddEntry(histoM02RatioallSM[0], labelMC1.Data(),"pe");
    // legendM02Dummy->AddEntry(histoM02RatioallSM[1], labelMC2.Data(),"pe");
    // legendM02Dummy->Draw();
    labelSupMod = new TLatex(0.25,0.86,Form("Supermodules %d-%d",minSM,maxSM-1));
    SetStyleTLatex( labelSupMod, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
    labelSupMod->Draw();
    histo2DM02RatioDummy->Draw("axis,same");
    canvasM02RatioPlots->SaveAs(Form("%s/Ratios/MCtrueDecomp_XTMC.%s",outputDir.Data(),suffix.Data()));


    histo2DM02RatioDummy->Draw("copy");
    histoM02FracTrueMergedPi0allSM[0]->Draw("p,same,e");
    histoM02FracTrueGammaFromPi0allSM[0]->Draw("p,same,e");
    histoM02FracTrueMergedPartConvPi0allSM[0]->Draw("p,same,e");
    histoM02FracTrueElectronFromPi0allSM[0]->Draw("p,same,e");
    histoM02FracTrueMergedPi0allSM[1]->Draw("p,same,e");
    histoM02FracTrueGammaFromPi0allSM[1]->Draw("p,same,e");
    histoM02FracTrueMergedPartConvPi0allSM[1]->Draw("p,same,e");
    histoM02FracTrueElectronFromPi0allSM[1]->Draw("p,same,e");
    // histoM02RatioallSM[1]->Draw("p,same,e");
    // TF1* pol0All = new TF1("pol0","[0]",0.27,1.5);//
    // histoM02RatioallSM[0]->Fit(pol0All,"NRMEX0+","",0.27,1.5);
    // pol0All->SetLineColor(kGray+2);
    DrawGammaLines( 0.1, 0.1, -0.1, 2.0, 3, kGreen+3, 9);
    DrawGammaLines( 0.27, 0.27, -0.1, 2.0, 3, kRed+3, 9);
    // TBox* boxPol0FitError = CreateBoxConv(kGray+1, 0.27,  pol0All->GetParameter(0)-pol0All->GetParError(0), 1.5, pol0All->GetParameter(0)+pol0All->GetParError(0));
    // boxPol0FitError->SetLineWidth(8);
    // boxPol0FitError->Draw();
    // histoM02RatioallSM[0]->Draw("p,same,e");
    // pol0All->Draw("same");
    labelEnergy->Draw();
    labelClusterE->Draw();
    // legendM02Dummy  = GetAndSetLegend2(0.24, 0.85-(0.035*2*1.35), 0.52, 0.85, 0.85*textSizeLabelsPixel);
    // legendM02Dummy->AddEntry(histoM02RatioallSM[0], labelMC1.Data(),"pe");
    // legendM02Dummy->AddEntry(histoM02RatioallSM[1], labelMC2.Data(),"pe");
    // legendM02Dummy->Draw();
    labelSupMod = new TLatex(0.25,0.86,Form("Supermodules %d-%d",minSM,maxSM-1));
    SetStyleTLatex( labelSupMod, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
    labelSupMod->Draw();
    histo2DM02RatioDummy->Draw("axis,same");
    canvasM02RatioPlots->SaveAs(Form("%s/Ratios/MCtrueDecomp_bothMC.%s",outputDir.Data(),suffix.Data()));


}
