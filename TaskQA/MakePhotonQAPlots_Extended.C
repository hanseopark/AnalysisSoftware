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

void MakePhotonQAPlots_Extended(
    TString fileName   = "GammaConvV1_QA_5460001022092970003190000000.root",
    TString cutnumber   = "0005314140",
    Bool_t isMC = kFALSE,
    TString suffix   = "pdf"
){

    gROOT->Reset();
    gROOT->SetStyle("Plain");

    StyleSettingsThesis();
    SetPlotStyle();


    //********************************************************************************
    //*            File definition/ loading                                          *
    //********************************************************************************

    TFile *f                        = (TFile*)gROOT->GetListOfFiles()->FindObject(fileName.Data());
    if (!f) {
        f                           = new TFile(fileName.Data());
    }
    if (!f) cout << "main List not found" << endl;
    if (f->IsZombie()) {
        cout <<fileName.Data() <<" file does not exist" << endl;
        f->Close();
        delete f;
        return;
    }
    TString nameDirectory           = Form("GammaConv_%s",  cutnumber.Data());
    TDirectory* directoryConv           = (TDirectory*)f->Get(nameDirectory.Data());

  //___________________________________ Declaration of files _____________________________________________
    TString dateForOutput                           = ReturnDateStringForOutput();
    cout << dateForOutput.Data() << endl;
    TString collisionSystem                     = "pp, #sqrt{#it{s}} = 8 TeV";
    TString labelALICEforPlots                      = "ALICE";
    TString outputDir                               = Form("%s/%s/MakePhotonQAPlots_%s",suffix.Data(),dateForOutput.Data(),cutnumber.Data());
    cout << outputDir.Data() << endl;

    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec("mkdir -p "+outputDir + "/AlphaQt");
    gSystem->Exec("mkdir -p "+outputDir + "/Chi2PsiPair");

    //********************************************************************************
    //*      Definition of histograms for reconstructed Conversion Points            *
    //********************************************************************************
    Double_t projBinnin[16] = {
        0.0, 0.1, 0.2, 0.4, 0.7,
        1.0, 1.5, 2.0, 3.0, 4.0,
        6.0, 10., 15., 20., 30.,
        50.};
    TH3F* histoElectrondEdxEtaP                 = NULL;
    TH3F* histoPositrondEdxEtaP                 = NULL;
    TH3F* histoElectronNSigmadEdxEtaP           = NULL;
    TH3F* histoPositronNSigmadEdxEtaP           = NULL;
    TH3F* histoElectronTOFEtaP                  = NULL;
    TH3F* histoPositronTOFEtaP                  = NULL;
    TH3F* histoElectronNSigmaTOFEtaP            = NULL;
    TH3F* histoPositronNSigmaTOFEtaP            = NULL;
    TH3F* histoElectronNSigmaITSEtaP            = NULL;
    TH3F* histoPositronNSigmaITSEtaP            = NULL;
    TH3F* histoElectronITSdEdxEtaP              = NULL;
    TH3F* histoPositronITSdEdxEtaP              = NULL;
    TH3F* histonSigdEdxElnSigdEdxPosGammaP      = NULL;
    TH2F* histoElectronNSigmadEdxPhi            = NULL;
    TH2F* histoPositronNSigmadEdxPhi            = NULL;
    TH2F* histoElectronNSigmadEdxPhiEtaNeg      = NULL;
    TH2F* histoPositronNSigmadEdxPhiEtaNeg      = NULL;
    TH2F* histoElectronNSigmadEdxPhiEtaPos      = NULL;
    TH2F* histoPositronNSigmadEdxPhiEtaPos      = NULL;
    TH2F* histoElectronFClPt                    = NULL;
    TH2F* histoPositronFClPt                    = NULL;
    TH2F* histoElectronClPt                     = NULL;
    TH2F* histoPositronClPt                     = NULL;
    TH2F* histoGammaChi2NDFPt                   = NULL;
    TH2F* histoGammaChi2NDFPtEtaNeg             = NULL;
    TH2F* histoGammaChi2NDFPtEtaPos             = NULL;
    TH2F* histoGammaChi2NDFR                    = NULL;
    TH2F* histoElectronEtaPt                    = NULL;
    TH2F* histoPositronEtaPt                    = NULL;
    TH2F* histoElectronITSClPt                  = NULL;
    TH2F* histoPositronITSClPt                  = NULL;
    TH2F* histoElectronITSClEta                 = NULL;
    TH2F* histoPositronITSClEta                 = NULL;
    TH2F* histoPositronITSClR                   = NULL;
    TH2F* histoElectronITSClR                   = NULL;
    TH2F* histoPositronITSClPhi                 = NULL;
    TH2F* histoElectronITSClPhi                 = NULL;
    TH2F* histoPositronClR                      = NULL;
    TH2F* histoElectronClR                      = NULL;
    TH2F* histoPositronClPhi                    = NULL;
    TH2F* histoElectronClPhi                    = NULL;
    TH1F* histoGammaPhiEtaNeg                   = NULL;
    TH1F* histoGammaPhiEtaPos                   = NULL;
    TH1F* histoPositronNSigmadEdxCut            = NULL;
    TH1F* histoElectronNSigmadEdxCut            = NULL;
    TH1F *histoVertexZ                          = NULL;
    TH1F *histoContrVertexZ                     = NULL;
    TH1F *histoGoodESDTracks                    = NULL;
    TH1F *histoVertexZdummy                     = NULL;
    TH1F *histoContrVertexZdummy                = NULL;
    TH1F *histoGoodESDTracksdummy               = NULL;
    TH2F* histoGammaAsymP                       = NULL;
    TH2F* histoGammaAsymR                       = NULL;
    TH2F* histoGammaAlphaQt                     = NULL;
    TH3F* histoGammaAlphaQtPt                   = NULL;
    TH3F* histoGammaAlphaQtPtCopy               = NULL;
    TH2F* histoGammaAlphaQtPtSliced[15]         = {NULL};
    TH2F* histoGammaAlphaPt                     = NULL;
    TH2F* histoGammaAlphaR                      = NULL;
    TH2F* histoGammaQtPt                        = NULL;
    TH2F* histoGammaQtR                         = NULL;
    TH2F* histoGammaEtaPt                       = NULL;
    TH2F* histoGammaEtaR                        = NULL;
    TH2F* histoGammaPsiPairPt                   = NULL;
    TH2F* histoGammaPsiPairPtEtaNeg             = NULL;
    TH2F* histoGammaPsiPairPtEtaPos             = NULL;
    TH2F* histoGammaPsiPairR                    = NULL;
    TH2F* histoGammaCosPointPt                  = NULL;
    TH2F* histoGammaCosPointR                   = NULL;
    TH2F* histoGammaCosPointChi2                = NULL;
    TH2F* histoGammaInvMassPt                   = NULL;
    TH2F* histoGammaInvMassR                    = NULL;
    TH2F* histoGammaCosPointPsiPair             = NULL;
    TH2F* histoGammaChi2PsiPair                 = NULL;
    TH3F* histoGammaChi2PsiPairPt               = NULL;
    TH3F* histoGammaChi2PsiPairPtCopy           = NULL;
    TH2F* histoGammaChi2PsiPairPtSliced[15]     = {NULL};
    TH2F* histoGammaChi2Pt                      = NULL;
    TH1F* histoGammaPhi                         = NULL;
    TH2F* histoGammaPhiR                        = NULL;
    TH2F* histoGammaZR                          = NULL;
    TH2F* histoGammaXY                          = NULL;

    TH2F* histoTruePrimGammaChi2NDFPt           = NULL;
    TH2F* histoTruePrimGammaAlphaQt             = NULL;
    TH3F* histoTruePrimGammaAlphaQtPt           = NULL;
    TH2F* histoTruePrimGammaAlphaPt             = NULL;
    TH2F* histoTruePrimGammaQtPt                = NULL;
    TH2F* histoTruePrimGammaQtR                 = NULL;
    TH2F* histoTruePrimGammaPsiPairPt           = NULL;
    TH2F* histoTruePrimGammaCosPointPt          = NULL;
    TH2F* histoTruePrimGammaAsymP               = NULL;
    TH2F* histoTrueMCKindQtPt[20]               = {NULL};
    TH2F* histoTrueMCKindAlphaPt[20]            = {NULL};
    TH2F* histoTrueMCKindPsiPairPt[20]          = {NULL};
    TH2F* histoTrueMCKindChi2Pt[20]             = {NULL};
    TH2F* histoTrueMCKindCosPointPt[20]         = {NULL};
    TH2F* histoTrueMCKindChi2PsiPair[20]        = {NULL};
    TH3F* histoTrueMCKindChi2PsiPairPt[20]      = {NULL};
    TH3F* histoTrueMCKindChi2PsiPairPtCopy[20]  = {NULL};
    TH2F* histoTrueMCKindChi2PsiPairPtSliced[20][15]    = {NULL};
    TH2F* histoTrueMCKindAlphaQt[20]            = {NULL};
    TH3F* histoTrueMCKindAlphaQtPt[20]          = {NULL};
    TH3F* histoTrueMCKindAlphaQtPtCopy[20]      = {NULL};
    TH2F* histoTrueMCKindAlphaQtPtSliced[20][15]= {NULL};
    TH2F* histoTrueSecGammaChi2NDFPt            = NULL;
    TH2F* histoTrueSecGammaAlphaQt              = NULL;
    TH2F* histoTrueSecGammaQtPt                 = NULL;
    TH2F* histoTrueSecGammaPsiPairPt            = NULL;
    TH2F* histoTrueSecGammaCosPointPt           = NULL;
    TH2F* histoTrueSecGammaAsymP                = NULL;
    TH2F* histoTrueDalitzGammaChi2NDFPt         = NULL;
    TH2F* histoTrueDalitzGammaAlphaQt           = NULL;
    TH2F* histoTrueDalitzGammaQtPt              = NULL;
    TH2F* histoTrueDalitzGammaPsiPairPt         = NULL;
    TH2F* histoTrueDalitzGammaCosPointPt        = NULL;
    TH2F* histoTrueDalitzGammaAsymP             = NULL;
    TH2F* histoTrueGammaCosPointChi2            = NULL;
    TH2F* histoTrueGammaInvMassPt               = NULL;
    TH2F* histoTrueGammaCosPointPsiPair         = NULL;
    TH2F* histoTrueGammaChi2PsiPair             = NULL;
    TH3F* histoTrueGammaChi2PsiPairPt           = NULL;
    TH2F* histoTrueDalitzGammaCosPointChi2      = NULL;
    TH2F* histoTrueDalitzGammaInvMassPt         = NULL;
    TH2F* histoTrueDalitzGammaCosPointPsiPair   = NULL;
    TH2F* histoTrueDalitzGammaChi2PsiPair       = NULL;
    TH2F* histoTrueElectronNSigmadEdxTPCP       = NULL;
    TH2F* histoTruePositronNSigmadEdxTPCP       = NULL;
    TH2F* histoTrueElectronNSigmadEdxITSP       = NULL;
    TH2F* histoTruePositronNSigmadEdxITSP       = NULL;
    TH2F* histoTrueElectronNSigmaTOFP           = NULL;
    TH2F* histoTruePositronNSigmaTOFP           = NULL;
    TH2F* histoTruePositronITSClR               = NULL;
    TH2F* histoTrueElectronITSClR               = NULL;

    histoElectrondEdxEtaP                   = (TH3F*)directoryConv->Get("histoElectrondEdxEtaP");
    histoPositrondEdxEtaP                   = (TH3F*)directoryConv->Get("histoPositrondEdxEtaP");
    histoElectronNSigmadEdxEtaP             = (TH3F*)directoryConv->Get("histoElectronNSigmadEdxEtaP");
    histoPositronNSigmadEdxEtaP             = (TH3F*)directoryConv->Get("histoPositronNSigmadEdxEtaP");
    histoElectronTOFEtaP                    = (TH3F*)directoryConv->Get("histoElectronTOFEtaP");
    histoPositronTOFEtaP                    = (TH3F*)directoryConv->Get("histoPositronTOFEtaP");
    histoElectronNSigmaTOFEtaP              = (TH3F*)directoryConv->Get("histoElectronNSigmaTOFEtaP");
    histoPositronNSigmaTOFEtaP              = (TH3F*)directoryConv->Get("histoPositronNSigmaTOFEtaP");
    histoElectronNSigmadEdxPhi              = (TH2F*)directoryConv->Get("histoElectronNSigmadEdxPhi");
    histoPositronNSigmadEdxPhi              = (TH2F*)directoryConv->Get("histoPositronNSigmadEdxPhi");
    histoElectronNSigmadEdxPhiEtaNeg        = (TH2F*)directoryConv->Get("histoElectronNSigmadEdxPhiEtaNeg");
    histoPositronNSigmadEdxPhiEtaNeg        = (TH2F*)directoryConv->Get("histoPositronNSigmadEdxPhiEtaNeg");
    histoElectronNSigmadEdxPhiEtaPos        = (TH2F*)directoryConv->Get("histoElectronNSigmadEdxPhiEtaPos");
    histoPositronNSigmadEdxPhiEtaPos        = (TH2F*)directoryConv->Get("histoPositronNSigmadEdxPhiEtaPos");
    histoElectronFClPt                      = (TH2F*)directoryConv->Get("histoElectronFClPt");
    histoPositronFClPt                      = (TH2F*)directoryConv->Get("histoPositronFClPt");
    histoElectronClPt                       = (TH2F*)directoryConv->Get("histoElectronClPt");
    histoPositronClPt                       = (TH2F*)directoryConv->Get("histoPositronClPt");
    histoGammaChi2NDFPt                     = (TH2F*)directoryConv->Get("histoGammaChi2NDFPt");
    histoGammaChi2NDFPtEtaNeg               = (TH2F*)directoryConv->Get("histoGammaChi2NDFPtEtaNeg");
    histoGammaChi2NDFPtEtaPos               = (TH2F*)directoryConv->Get("histoGammaChi2NDFPtEtaPos");
    histoGammaChi2NDFR                      = (TH2F*)directoryConv->Get("histoGammaChi2NDFR");
    histoElectronEtaPt                      = (TH2F*)directoryConv->Get("histoElectronEtaPt");
    histoPositronEtaPt                      = (TH2F*)directoryConv->Get("histoPositronEtaPt");
    histoElectronITSClR                     = (TH2F*)directoryConv->Get("histoElectronITSClR");
    histoPositronITSClR                     = (TH2F*)directoryConv->Get("histoPositronITSClR");
    histoElectronITSClPhi                   = (TH2F*)directoryConv->Get("histoElectronITSClPhi");
    histoPositronITSClPhi                   = (TH2F*)directoryConv->Get("histoPositronITSClPhi");
    histoElectronClR                        = (TH2F*)directoryConv->Get("histoElectronClR");
    histoPositronClR                        = (TH2F*)directoryConv->Get("histoPositronClR");
    histoElectronClPhi                      = (TH2F*)directoryConv->Get("histoElectronClPhi");
    histoPositronClPhi                      = (TH2F*)directoryConv->Get("histoPositronClPhi");
    histoGammaPhiEtaNeg                     = (TH1F*)directoryConv->Get("histoGammaPhiEtaNeg");
    histoGammaPhiEtaPos                     = (TH1F*)directoryConv->Get("histoGammaPhiEtaPos");
    histoPositronNSigmadEdxCut              = (TH1F*)directoryConv->Get("histoPositronNSigmadEdxCut");
    histoElectronNSigmadEdxCut              = (TH1F*)directoryConv->Get("histoElectronNSigmadEdxCut");
    histoGammaEtaPt                         = (TH2F*)directoryConv->Get("histoGammaEtaPt");
    histoGammaEtaR                          = (TH2F*)directoryConv->Get("histoGammaEtaR");
    histoGammaAlphaQt                       = (TH2F*)directoryConv->Get("histoGammaAlphaQt");
    histoGammaAlphaQtPt                     = (TH3F*)directoryConv->Get("histoGammaAlphaQtPt");
    histoGammaAlphaPt                       = (TH2F*)directoryConv->Get("histoGammaAlphaPt");
    histoGammaAlphaR                        = (TH2F*)directoryConv->Get("histoGammaAlphaR");
    histoGammaQtPt                          = (TH2F*)directoryConv->Get("histoGammaQtPt");
    histoGammaQtR                           = (TH2F*)directoryConv->Get("histoGammaQtR");
    histoGammaPhi                           = (TH1F*)directoryConv->Get( "histoGammaPhi");
    histoGammaPhiR                          = (TH2F*)directoryConv->Get( "histoGammaPhiR");
    histoGammaPsiPairPt                     = (TH2F*)directoryConv->Get( "histoGammaPsiPairPt");
    histoGammaPsiPairPtEtaNeg               = (TH2F*)directoryConv->Get( "histoGammaPsiPairPtEtaNeg");
    histoGammaPsiPairPtEtaPos               = (TH2F*)directoryConv->Get( "histoGammaPsiPairPtEtaPos");
    histoGammaPsiPairR                      = (TH2F*)directoryConv->Get( "histoGammaPsiPairR");
    histoGammaCosPointPt                    = (TH2F*)directoryConv->Get( "histoGammaCosPointPt");
    histoGammaCosPointR                     = (TH2F*)directoryConv->Get( "histoGammaCosPointR");
    histoGammaAsymP                         = (TH2F*)directoryConv->Get("histoGammaAsymP");
    histoGammaAsymR                         = (TH2F*)directoryConv->Get("histoGammaAsymR");
    histoGammaChi2PsiPair                   = (TH2F*)directoryConv->Get("histoGammaChi2PsiPair");
    histoGammaChi2PsiPairPt                 = (TH3F*)directoryConv->Get("histoGammaChi2PsiPairPt");
    histoGammaChi2Pt                        = (TH2F*)directoryConv->Get("histoGammaChi2Pt");
    histoGammaCosPointChi2                  = (TH2F*)directoryConv->Get("histoGammaCosPointChi2");
    histoGammaCosPointPsiPair               = (TH2F*)directoryConv->Get("histoGammaCosPointPsiPair");
    histoGammaInvMassPt                     = (TH2F*)directoryConv->Get("histoGammaInvMassPt");
    histoGammaInvMassR                      = (TH2F*)directoryConv->Get("histoGammaInvMassR");
    histoGammaZR                            = (TH2F*)directoryConv->Get("histoGammaZR");
    histoGammaXY                            = (TH2F*)directoryConv->Get("histoGammaXY");
    histoGammaChi2PsiPairPtCopy             = (TH3F*)histoGammaChi2PsiPairPt->Clone("histoGammaChi2PsiPairPtCopy");
    histoGammaAlphaQtPtCopy                 = (TH3F*)histoGammaAlphaQtPt->Clone("histoGammaAlphaQtPtCopy");
    for(Int_t j=0; j<15; j++){
        histoGammaChi2PsiPairPtCopy->GetZaxis()->SetRangeUser(projBinnin[j],projBinnin[j+1]);
        histoGammaChi2PsiPairPtSliced[j]   = (TH2F*)histoGammaChi2PsiPairPtCopy->Project3D(Form("yx%d",j));
        histoGammaAlphaQtPtCopy->GetZaxis()->SetRangeUser(projBinnin[j],projBinnin[j+1]);
        histoGammaAlphaQtPtSliced[j]       = (TH2F*)histoGammaAlphaQtPtCopy->Project3D(Form("yx%d",j));
    }
    //ITS add
    histoElectronITSdEdxEtaP                = (TH3F*)directoryConv->Get("histoElectronITSdEdxEtaP");
    histoPositronITSdEdxEtaP                = (TH3F*)directoryConv->Get("histoPositronITSdEdxEtaP");
    histonSigdEdxElnSigdEdxPosGammaP        = (TH3F*)directoryConv->Get("histonSigdEdxElnSigdEdxPosGammaP"); // new sept 7 2015
    histoElectronNSigmaITSEtaP              = (TH3F*)directoryConv->Get("histoElectronNSigmaITSEtaP");
    histoPositronNSigmaITSEtaP              = (TH3F*)directoryConv->Get("histoPositronNSigmaITSEtaP");
    histoElectronITSClPt                    = (TH2F*)directoryConv->Get("histoElectronITSClPt");
    histoPositronITSClPt                    = (TH2F*)directoryConv->Get("histoPositronITSClPt");
    histoElectronITSClEta                   = (TH2F*)directoryConv->Get("histoElectronITSClEta");
    histoPositronITSClEta                   = (TH2F*)directoryConv->Get("histoPositronITSClEta");

    if (isMC){
        histoTrueElectronNSigmadEdxITSP     = (TH2F*)directoryConv->Get("histoTrueElectronNSigmadEdxITSP");
        histoTrueElectronNSigmadEdxTPCP     = (TH2F*)directoryConv->Get("histoTrueElectronNSigmadEdxTPCP");
        histoTrueElectronNSigmaTOFP         = (TH2F*)directoryConv->Get("histoTrueElectronNSigmaTOFP");
        histoTrueElectronITSClR             = (TH2F*)directoryConv->Get("histoTrueElectronITSClR");
        histoTruePositronNSigmadEdxITSP     = (TH2F*)directoryConv->Get("histoTruePositronNSigmadEdxITSP");
        histoTruePositronNSigmadEdxTPCP     = (TH2F*)directoryConv->Get("histoTruePositronNSigmadEdxTPCP");
        histoTruePositronNSigmaTOFP         = (TH2F*)directoryConv->Get("histoTruePositronNSigmaTOFP");
        histoTruePositronITSClR             = (TH2F*)directoryConv->Get("histoTruePositronITSClR");
        histoTruePrimGammaChi2NDFPt         = (TH2F*)directoryConv->Get("histoTruePrimGammaChi2NDFPt");
        histoTruePrimGammaAlphaQt           = (TH2F*)directoryConv->Get("histoTruePrimGammaAlphaQt");
        histoTruePrimGammaAlphaQtPt         = (TH3F*)directoryConv->Get("histoTruePrimGammaAlphaQtPt");
        histoTruePrimGammaAlphaPt           = (TH2F*)directoryConv->Get("histoTruePrimGammaAlphaPt");
        histoTruePrimGammaQtPt              = (TH2F*)directoryConv->Get("histoTruePrimGammaQtPt");
        histoTruePrimGammaQtR               = (TH2F*)directoryConv->Get("histoTruePrimGammaQtR");
        histoTruePrimGammaPsiPairPt         = (TH2F*)directoryConv->Get("histoTruePrimGammaPsiPairPt");
        histoTruePrimGammaCosPointPt        = (TH2F*)directoryConv->Get("histoTruePrimGammaCosPointPt");
        histoTrueSecGammaChi2NDFPt          = (TH2F*)directoryConv->Get("histoTrueSecGammaChi2NDFPt");
        histoTrueSecGammaAlphaQt            = (TH2F*)directoryConv->Get("histoTrueSecGammaAlphaQt");
        histoTrueSecGammaQtPt               = (TH2F*)directoryConv->Get("histoTrueSecGammaQtPt");
        for(Int_t i=0;i<20;i++){
            if(i!=0 && i!=10 && i!=11 && i!=13) continue;
            histoTrueMCKindQtPt[i]          = (TH2F*)directoryConv->Get(Form("histoTrueMCKindQtPt_kind%d",i));
            histoTrueMCKindAlphaPt[i]       = (TH2F*)directoryConv->Get(Form("histoTrueMCKindAlphaPt_kind%d",i));
            histoTrueMCKindPsiPairPt[i]     = (TH2F*)directoryConv->Get(Form("histoTrueMCKindPsiPairPt_kind%d",i));
            histoTrueMCKindChi2Pt[i]        = (TH2F*)directoryConv->Get(Form("histoTrueMCKindChi2Pt_kind%d",i));
            histoTrueMCKindCosPointPt[i]    = (TH2F*)directoryConv->Get(Form("histoTrueMCKindCosPointPt_kind%d",i));

            histoTrueMCKindChi2PsiPair[i]   = (TH2F*)directoryConv->Get(Form("histoTrueMCKindChi2PsiPair_kind%d",i));
            histoTrueMCKindChi2PsiPairPt[i] = (TH3F*)directoryConv->Get(Form("histoTrueMCKindChi2PsiPairPt_kind%d",i));
            histoTrueMCKindAlphaQt[i]       = (TH2F*)directoryConv->Get(Form("histoTrueMCKindAlphaQt_kind%d",i));
            histoTrueMCKindAlphaQtPt[i]     = (TH3F*)directoryConv->Get(Form("histoTrueMCKindAlphaQtPt_kind%d",i));

            histoTrueMCKindChi2PsiPairPtCopy[i] = (TH3F*)histoTrueMCKindChi2PsiPairPt[i]->Clone(Form("histoTrueMCKindChi2PsiPairPtCopy%d",i));
            histoTrueMCKindAlphaQtPtCopy[i]     = (TH3F*)histoTrueMCKindAlphaQtPt[i]->Clone(Form("histoTrueMCKindAlphaQtPtCopy%d",i));
            for(Int_t j=0; j<15; j++){
                histoTrueMCKindChi2PsiPairPtCopy[i]->GetZaxis()->SetRangeUser(projBinnin[j],projBinnin[j+1]);
                histoTrueMCKindChi2PsiPairPtSliced[i][j]   = (TH2F*)histoTrueMCKindChi2PsiPairPtCopy[i]->Project3D(Form("yx%d",j));
                histoTrueMCKindAlphaQtPtCopy[i]->GetZaxis()->SetRangeUser(projBinnin[j],projBinnin[j+1]);
                histoTrueMCKindAlphaQtPtSliced[i][j]       = (TH2F*)histoTrueMCKindAlphaQtPtCopy[i]->Project3D(Form("yx%d",j));
            }
        }
        histoTrueSecGammaPsiPairPt          = (TH2F*)directoryConv->Get("histoTrueSecGammaPsiPairPt");
        histoTrueSecGammaCosPointPt         = (TH2F*)directoryConv->Get("histoTrueSecGammaCosPointPt");
        histoTrueDalitzGammaChi2NDFPt       = (TH2F*)directoryConv->Get("histoTrueDalitzGammaChi2NDFPt");
        histoTrueDalitzGammaAlphaQt         = (TH2F*)directoryConv->Get("histoTrueDalitzGammaAlphaQt");
        histoTrueDalitzGammaQtPt            = (TH2F*)directoryConv->Get("histoTrueDalitzGammaQtPt");
        histoTrueDalitzGammaPsiPairPt       = (TH2F*)directoryConv->Get("histoTrueDalitzGammaPsiPairPt");
        histoTrueDalitzGammaCosPointPt      = (TH2F*)directoryConv->Get("histoTrueDalitzGammaCosPointPt");
        histoTruePrimGammaAsymP             = (TH2F*)directoryConv->Get("histoTruePrimGammaAsymP");
        histoTrueSecGammaAsymP              = (TH2F*)directoryConv->Get("histoTrueSecGammaAsymP");
        histoTrueDalitzGammaAsymP           = (TH2F*)directoryConv->Get("histoTrueDalitzGammaAsymP");
        histoTrueGammaChi2PsiPair           = (TH2F*)directoryConv->Get("histoTrueGammaChi2PsiPair");
        histoTrueGammaChi2PsiPairPt         = (TH3F*)directoryConv->Get("histoTrueGammaChi2PsiPairPt");
        histoTrueGammaCosPointChi2          = (TH2F*)directoryConv->Get("histoTrueGammaCosPointChi2");
        histoTrueGammaCosPointPsiPair       = (TH2F*)directoryConv->Get("histoTrueGammaCosPointPsiPair");
        histoTrueGammaInvMassPt             = (TH2F*)directoryConv->Get("histoTrueGammaInvMassPt");
        histoTrueDalitzGammaChi2PsiPair     = (TH2F*)directoryConv->Get("histoTrueDalitzGammaChi2PsiPair");
        histoTrueDalitzGammaCosPointChi2    = (TH2F*)directoryConv->Get("histoTrueDalitzGammaCosPointChi2");
        histoTrueDalitzGammaCosPointPsiPair = (TH2F*)directoryConv->Get("histoTrueDalitzGammaCosPointPsiPair");
        histoTrueDalitzGammaInvMassPt       = (TH2F*)directoryConv->Get("histoTrueDalitzGammaInvMassPt");
    }

    TLatex* labelpTrange;

    //    ____ _______
    //   / __ \__   __|
    //  | |  | | | |
    //  | |  | | | |
    //  | |__| | | |
    //   \___\_\ |_|

    Int_t textSizeLabelsPixel                   = 1200*0.04;

    TCanvas* canvasQtPlots = new TCanvas("canvasQtPlots","",200,10,1350,1200);  // gives the page size
    DrawGammaCanvasSettings( canvasQtPlots, 0.08, 0.02, 0.02, 0.09);
    canvasQtPlots->SetLogy();
    canvasQtPlots->SetLogz();

    TH2F * histo2DQtDummy = new TH2F("histo2DQtDummy","histo2DQtDummy",100,0,0.1,1000,0.04,40);
    SetStyleHistoTH2ForGraphs(histo2DQtDummy, "#it{q}_{T}^{#gamma}","#it{p}_{T} (GeV/#it{c})",0.035,0.04, 0.035,0.04, 0.98,0.9);
    histo2DQtDummy->GetZaxis()->SetRangeUser(6,6e2);
    canvasQtPlots->cd();
    histo2DQtDummy->Draw("copy");
    TExec *ex1 = NULL;
    TExec *ex2 = NULL;
    if (isMC){
        // draw true gamma->ee
        histoTrueMCKindQtPt[0]->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTrueMCKindQtPt[0]->Draw("col,same");
        // draw ee combinatorics
        histoTrueMCKindQtPt[11]->Draw("col,same");
        ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        ex2->Draw();
        histoTrueMCKindQtPt[11]->Draw("col,same");
        // draw pipi combinatorics
        histoTrueMCKindQtPt[13]->Draw("col,same");
        ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        ex2->Draw();
        histoTrueMCKindQtPt[13]->Draw("col,same");
    } else {
        histoGammaQtPt->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoGammaQtPt->Draw("col,same");
    }
    TF1 *funcPtDepQtCut_std = new TF1("funcPtDepQtCut_std","x/[0]",0.005,1.);
    DrawGammaSetMarkerTF1( funcPtDepQtCut_std, 1, 3, kMagenta+2);
    funcPtDepQtCut_std->SetParameter(0,0.125);
    funcPtDepQtCut_std->Draw("same");
    TF1 *funcPtDepQtCut_hard = new TF1("funcPtDepQtCut_hard","x/[0]",0.005,1.);
    DrawGammaSetMarkerTF1( funcPtDepQtCut_hard, 2, 3, kMagenta+4);
    funcPtDepQtCut_hard->SetParameter(0,0.11);
    funcPtDepQtCut_hard->Draw("same");
    TF1 *funcPtDepQtCut_soft = new TF1("funcPtDepQtCut_soft","x/[0]",0.005,1.);
    DrawGammaSetMarkerTF1( funcPtDepQtCut_soft, 2, 3, kMagenta-4);
    funcPtDepQtCut_soft->SetParameter(0,0.14);
    funcPtDepQtCut_soft->Draw("same");
    TF1 *funcPtDepQtCut_soft2 = new TF1("funcPtDepQtCut_soft2","x/[0]",0.005,1.);
    DrawGammaSetMarkerTF1( funcPtDepQtCut_soft2, 2, 3, kMagenta-1);
    funcPtDepQtCut_soft2->SetParameter(0,0.16);
    funcPtDepQtCut_soft2->Draw("same");
    DrawGammaLines(0.05, 0.05, 0.04, 8, 3, kGreen-8, 9);
    TF1 *funcOldCutDummy = new TF1("funcOldCutDummy","x/[0]",0.005,1.);
    DrawGammaSetMarkerTF1( funcOldCutDummy, 9, 3, kGreen-8);
    TLegend* legendQtPlotFits  = GetAndSetLegend2(0.67, 0.14, 0.95, 0.14+(0.035*5*1.35), 0.85*textSizeLabelsPixel);
    legendQtPlotFits->AddEntry(funcOldCutDummy, "#it{q}_{T}^{max} = 0.05","l");
    legendQtPlotFits->AddEntry(funcPtDepQtCut_hard, "#it{q}_{T}^{max} = 0.110*#it{p}_{T}","l");
    legendQtPlotFits->AddEntry(funcPtDepQtCut_std,  "#it{q}_{T}^{max} = 0.125*#it{p}_{T}","l");
    legendQtPlotFits->AddEntry(funcPtDepQtCut_soft, "#it{q}_{T}^{max} = 0.140*#it{p}_{T}","l");
    legendQtPlotFits->AddEntry(funcPtDepQtCut_soft2,"#it{q}_{T}^{max} = 0.160*#it{p}_{T}","l");
    legendQtPlotFits->Draw();
    TLatex *labelColor                  = new TLatex(0.12,0.90, isMC ? "#gamma rec. from e^{+}e^{-} (color)" : "");
    SetStyleTLatex( labelColor, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
    labelColor->Draw();
    TLatex *labelBW                     = new TLatex(0.12,0.86, isMC ? "#gamma rec. from e^{#pm}#pi^{#pm}/#pi^{+}#pi^{-} (BAW)" : "");
    SetStyleTLatex( labelBW, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
    labelBW->Draw();
    TLatex *labelEnergy                  = new TLatex(0.95,0.90,Form("%s, %s",labelALICEforPlots.Data(),collisionSystem.Data()));
    SetStyleTLatex( labelEnergy, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelEnergy->Draw();
    TLatex *labelProcess                     = new TLatex(0.95,0.86,isMC ? "#gamma candidates (MC rec.)" : "#gamma candidates (data)");
    SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelProcess->Draw();
    histo2DQtDummy->Draw("axis,same");
    canvasQtPlots->SaveAs(Form("%s/Qt_vs_Pt_Final_withCuts_%s.%s",outputDir.Data(),isMC ? "MC" : "data",suffix.Data()));

    histo2DQtDummy->Draw("copy");
    if (isMC){
        histoTrueMCKindQtPt[0]->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTrueMCKindQtPt[0]->Draw("col,same");
    }
    DrawGammaLines(0.05, 0.05, 0.04, 8, 3, kGreen-8, 9);
    legendQtPlotFits  = GetAndSetLegend2(0.67, 0.14, 0.95, 0.14+(0.035*1*1.35), 0.85*textSizeLabelsPixel);
    legendQtPlotFits->AddEntry(funcOldCutDummy, "#it{q}_{T}^{max} = 0.05","l");
    legendQtPlotFits->Draw();
    labelColor->Draw();
    labelEnergy->Draw();
    labelProcess->Draw();
    histo2DQtDummy->Draw("axis,same");
    if (isMC)
        canvasQtPlots->SaveAs(Form("%s/Qt_vs_Pt_trueGamma_woCuts.%s",outputDir.Data(),suffix.Data()));



    histo2DQtDummy->Draw("copy");
    if (isMC){
        histoTrueMCKindQtPt[0]->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTrueMCKindQtPt[0]->Draw("col,same");
        histoTrueMCKindQtPt[11]->Draw("col,same");
        ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        ex2->Draw();
        histoTrueMCKindQtPt[11]->Draw("col,same");
    }
    DrawGammaLines(0.05, 0.05, 0.04, 8, 3, kGreen-8, 9);
    legendQtPlotFits  = GetAndSetLegend2(0.67, 0.14, 0.95, 0.14+(0.035*1*1.35), 0.85*textSizeLabelsPixel);
    legendQtPlotFits->AddEntry(funcOldCutDummy, "#it{q}_{T}^{max} = 0.05","l");
    legendQtPlotFits->Draw();
    labelBW                     = new TLatex(0.12,0.86,"#gamma rec. from #pi^{+}#pi^{-} (BAW)");
    SetStyleTLatex( labelBW, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
    labelBW->Draw();
    labelColor->Draw();
    labelEnergy->Draw();
    labelProcess->Draw();
    histo2DQtDummy->Draw("axis,same");
    if (isMC)
        canvasQtPlots->SaveAs(Form("%s/Qt_vs_Pt_trueGammaAndComb_woCuts.%s",outputDir.Data(),suffix.Data()));


    histo2DQtDummy->Draw("copy");
    if (isMC){
        histoTrueMCKindQtPt[0]->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTrueMCKindQtPt[0]->Draw("col,same");
        histoTrueMCKindQtPt[11]->Draw("col,same");
        ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        ex2->Draw();
        histoTrueMCKindQtPt[11]->Draw("col,same");
        histoTrueMCKindQtPt[13]->Draw("col,same");
        ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        ex2->Draw();
        histoTrueMCKindQtPt[13]->Draw("col,same");
    }
    DrawGammaLines(0.05, 0.05, 0.04, 8, 3, kGreen-8, 9);
    legendQtPlotFits  = GetAndSetLegend2(0.67, 0.14, 0.95, 0.14+(0.035*1*1.35), 0.85*textSizeLabelsPixel);
    legendQtPlotFits->AddEntry(funcOldCutDummy, "#it{q}_{T}^{max} = 0.05","l");
    legendQtPlotFits->Draw();
    labelBW                     = new TLatex(0.12,0.86,isMC ? "#gamma rec. from e^{#pm}#pi^{#pm}/#pi^{+}#pi^{-} (BAW)" : "");
    SetStyleTLatex( labelBW, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
    labelBW->Draw();
    labelColor->Draw();
    labelEnergy->Draw();
    labelProcess->Draw();
    histo2DQtDummy->Draw("axis,same");
    if (isMC)
        canvasQtPlots->SaveAs(Form("%s/Qt_vs_Pt_trueGammaAndCombPlus_woCuts.%s",outputDir.Data(),suffix.Data()));



    histo2DQtDummy->Draw("copy");
    histoGammaQtPt->Draw("col,same");
    ex1 = new TExec("ex1","PalColor();");
    ex1->Draw();
    histoGammaQtPt->Draw("col,same");
    DrawGammaLines(0.05, 0.05, 0.04, 8, 3, kGreen-8, 9);
    legendQtPlotFits  = GetAndSetLegend2(0.67, 0.14, 0.95, 0.14+(0.035*1*1.35), 0.85*textSizeLabelsPixel);
    legendQtPlotFits->AddEntry(funcOldCutDummy, "#it{q}_{T}^{max} = 0.05","l");
    legendQtPlotFits->Draw();
    labelEnergy->Draw();
    labelProcess                     = new TLatex(0.95,0.86,isMC ? "#gamma candidates (MC rec.)" : "#gamma candidates (data)");
    SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelProcess->Draw();
    histo2DQtDummy->Draw("axis,same");
    canvasQtPlots->SaveAs(Form("%s/Qt_vs_Pt_AllRec_%s.%s",outputDir.Data(),isMC ? "MC" : "data",suffix.Data()));


    //             _______     ____  __ __  __ ______ _______ _______     __
    //      /\    / ____\ \   / /  \/  |  \/  |  ____|__   __|  __ \ \   / /
    //     /  \  | (___  \ \_/ /| \  / | \  / | |__     | |  | |__) \ \_/ /
    //    / /\ \  \___ \  \   / | |\/| | |\/| |  __|    | |  |  _  / \   /
    //   / ____ \ ____) |  | |  | |  | | |  | | |____   | |  | | \ \  | |
    //  /_/    \_\_____/   |_|  |_|  |_|_|  |_|______|  |_|  |_|  \_\ |_|

    textSizeLabelsPixel                   = 1200*0.04;
    TCanvas* canvasAsymPlots = new TCanvas("canvasAsymPlots","",200,10,1350,1200);  // gives the page size
    DrawGammaCanvasSettings( canvasAsymPlots, 0.08, 0.02, 0.02, 0.09);
    canvasAsymPlots->SetLogy();
    canvasAsymPlots->SetLogz();
    TH2F * histo2DAlphaDummy = new TH2F("histo2DAlphaDummy","histo2DAlphaDummy",100,-1.07,1.07,1000,0.04,200);
    SetStyleHistoTH2ForGraphs(histo2DAlphaDummy, "#alpha^{#gamma} = (#it{p}^{+}_{L}-#it{p}^{-}_{L})/(#it{p}^{+}_{L}+#it{p}^{-}_{L})","#it{p}_{T} (GeV/#it{c})",0.035,0.04, 0.035,0.04, 0.98,0.9);
    histo2DAlphaDummy->GetZaxis()->SetRangeUser(1,6e2);
    canvasAsymPlots->cd();

    // RECONSTRUCTED PLOT
    histo2DAlphaDummy->Draw("copy");
    histoGammaAlphaPt->Draw("col,same");
    ex1 = new TExec("ex1","PalColor();");
    ex1->Draw();
    histoGammaAlphaPt->Draw("col,same");
    labelEnergy                  = new TLatex(0.95,0.90,Form("%s, %s",labelALICEforPlots.Data(),collisionSystem.Data()));
    SetStyleTLatex( labelEnergy, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelEnergy->Draw();
    labelProcess                     = new TLatex(0.95,0.86,isMC ? "#gamma candidates (MC rec.)" : "#gamma candidates (data)");
    SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelProcess->Draw();
    histo2DAlphaDummy->Draw("axis,same");
    canvasAsymPlots->SaveAs(Form("%s/Alpha_vs_Pt_AllRec_%s.%s",outputDir.Data(),isMC ? "MC" : "data",suffix.Data()));


    // TRUE GAMMAS PLOT
    histo2DAlphaDummy->Draw("copy");
    // draw true gamma->ee
    if (isMC){
        histoTrueMCKindAlphaPt[0]->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTrueMCKindAlphaPt[0]->Draw("col,same");
    }
    labelColor                  = new TLatex(0.12,0.90, isMC ? "#gamma rec. from e^{+}e^{-} (color)" : "");
    SetStyleTLatex( labelColor, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
    labelColor->Draw();
    labelProcess                     = new TLatex(0.95,0.86,isMC ? "#gamma candidates (MC true)" : "#gamma candidates (data)");
    SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelEnergy->Draw();
    labelProcess->Draw();
    histo2DAlphaDummy->Draw("axis,same");
    if (isMC)
        canvasAsymPlots->SaveAs(Form("%s/Alpha_vs_Pt_trueGamma_woCuts.%s",outputDir.Data(),suffix.Data()));


    // TRUE GAMMAS AND COMBINATORIAL GAMMA CANDIDATES PLOT
    histo2DAlphaDummy->Draw("copy");
    if (isMC){
        // draw true gamma->ee
        histoTrueMCKindAlphaPt[0]->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTrueMCKindAlphaPt[0]->Draw("col,same");
        // draw ee combinatorics
        histoTrueMCKindAlphaPt[11]->Draw("col,same");
        ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        ex2->Draw();
        histoTrueMCKindAlphaPt[11]->Draw("col,same");
        // draw pipi combinatorics
        histoTrueMCKindAlphaPt[13]->Draw("col,same");
        ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        ex2->Draw();
        histoTrueMCKindAlphaPt[13]->Draw("col,same");
    }
    labelColor->Draw();
    labelBW                     = new TLatex(0.12,0.86,isMC ? "#gamma rec. from e^{#pm}#pi^{#pm}/#pi^{+}#pi^{-} (BAW)" : "");
    SetStyleTLatex( labelBW, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
    labelBW->Draw();
    labelEnergy->Draw();
    labelProcess->Draw();
    histo2DAlphaDummy->Draw("axis,same");
    if (isMC)
        canvasAsymPlots->SaveAs(Form("%s/Alpha_vs_Pt_trueGammaAndCombPlus_woCuts.%s",outputDir.Data(),suffix.Data()));


    // TRUE GAMMAS AND COMBINATORIAL GAMMA CANDIDATES PLOT
    histo2DAlphaDummy->Draw("copy");
    if (isMC){
        // draw true gamma->ee
        histoTrueMCKindAlphaPt[0]->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTrueMCKindAlphaPt[0]->Draw("col,same");
        // draw ee combinatorics
        histoTrueMCKindAlphaPt[11]->Draw("col,same");
        ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        ex2->Draw();
        histoTrueMCKindAlphaPt[11]->Draw("col,same");
        // draw pipi combinatorics
        histoTrueMCKindAlphaPt[13]->Draw("col,same");
        ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        ex2->Draw();
        histoTrueMCKindAlphaPt[13]->Draw("col,same");
    } else {
        histoGammaAlphaPt->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoGammaAlphaPt->Draw("col,same");
    }
    // TF1 *funcPtDepAlphaCut_std = new TF1("funcPtDepAlphaCut_std","[0]*TMath::Exp(2*TMath::Abs(x))",-1,1.);
    // DrawGammaSetMarkerTF1( funcPtDepAlphaCut_std, 1, 3, kMagenta+2);
    // funcPtDepAlphaCut_std->SetParameter(0,0.07);
    // funcPtDepAlphaCut_std->Draw("same");
    // TF1 *funcPtDepAlphaCut_hard = new TF1("funcPtDepAlphaCut_hard","[0]*TMath::Exp(2.1*TMath::Abs(x))",-1,1.);
    // DrawGammaSetMarkerTF1( funcPtDepAlphaCut_hard, 2, 3, kMagenta+4);
    // funcPtDepAlphaCut_hard->SetParameter(0,0.07);
    // funcPtDepAlphaCut_hard->Draw("same");
    // TF1 *funcPtDepAlphaCut_soft = new TF1("funcPtDepAlphaCut_soft","[0]*TMath::Exp(1.9*TMath::Abs(x))",-1,1.);
    // DrawGammaSetMarkerTF1( funcPtDepAlphaCut_soft, 2, 3, kMagenta-4);
    // funcPtDepAlphaCut_soft->SetParameter(0,0.07);
    // funcPtDepAlphaCut_soft->Draw("same");
    // TF1 *funcPtDepAlphaCut_soft2 = new TF1("funcPtDepAlphaCut_soft2","[0]*TMath::Exp(1.9*TMath::Abs(x))",-1,1.);
    // DrawGammaSetMarkerTF1( funcPtDepAlphaCut_soft2, 2, 3, kMagenta-1);
    // funcPtDepAlphaCut_soft2->SetParameter(0,0.06);
    // funcPtDepAlphaCut_soft2->Draw("same");
    DrawGammaLines(-0.95, -0.95, 0.04, 50, 3, kGreen+3, 9);
    DrawGammaLines( 0.95,  0.95, 0.04, 50, 3, kGreen+3, 9);
    funcOldCutDummy = new TF1("funcOldCutDummy","x/[0]",0.005,1.);
    DrawGammaSetMarkerTF1( funcOldCutDummy, 9, 3, kGreen+3);
    // DrawGammaLines(-0.99, -0.99, 0.04, 50, 3, kGreen-8, 9);
    // DrawGammaLines( 0.99,  0.99, 0.04, 50, 3, kGreen-8, 9);
    // TF1* funcOldCutDummy1 = new TF1("funcOldCutDummy","x/[0]",0.005,1.);
    // DrawGammaSetMarkerTF1( funcOldCutDummy1, 9, 3, kGreen-8);
    TLegend* legendAlphaPlotFits  = GetAndSetLegend2(0.67, 0.13, 0.95, 0.13+(0.035*1*1.35), 0.85*textSizeLabelsPixel);
    legendAlphaPlotFits->AddEntry(funcOldCutDummy, "|#alpha^{max}| = 0.95","l");
    // legendAlphaPlotFits->AddEntry(funcOldCutDummy1, "|#alpha^{max}| = 0.99","l");
    // legendAlphaPlotFits->AddEntry(funcPtDepAlphaCut_hard, "#it{q}_{T}^{max} = 0.110*#it{p}_{T}","l");
    // legendAlphaPlotFits->AddEntry(funcPtDepAlphaCut_std,  "#it{q}_{T}^{max} = 0.125*#it{p}_{T}","l");
    // legendAlphaPlotFits->AddEntry(funcPtDepAlphaCut_soft, "#it{q}_{T}^{max} = 0.140*#it{p}_{T}","l");
    // legendAlphaPlotFits->AddEntry(funcPtDepAlphaCut_soft2,"#it{q}_{T}^{max} = 0.160*#it{p}_{T}","l");
    legendAlphaPlotFits->Draw();
    TLatex* labelHighpTSignal    = new TLatex(0.90,0.78,"high #it{p}_{T} signal #rightarrow");
    SetStyleTLatex( labelHighpTSignal, 0.95*textSizeLabelsPixel,4,kBlue+2,43,kTRUE,31);
    labelHighpTSignal->Draw();
    labelColor->Draw();
    labelBW->Draw();
    labelEnergy->Draw();
    labelProcess->Draw();
    histo2DAlphaDummy->Draw("axis,same");
    canvasAsymPlots->SaveAs(Form("%s/Alpha_vs_Pt_Final_withCuts_%s.%s",outputDir.Data(),isMC ? "MC" : "data",suffix.Data()));




    //    _____ _    _ _____ ___
    //   / ____| |  | |_   _|__ \
    //  | |    | |__| | | |    ) |
    //  | |    |  __  | | |   / /
    //  | |____| |  | |_| |_ / /_
    //   \_____|_|  |_|_____|____|


    textSizeLabelsPixel                   = 1200*0.04;
    TCanvas* canvasChi2Plots = new TCanvas("canvasChi2Plots","",200,10,1350,1200);  // gives the page size
    DrawGammaCanvasSettings( canvasChi2Plots, 0.08, 0.02, 0.02, 0.09);
    canvasChi2Plots->SetLogy();
    canvasChi2Plots->SetLogz();
    TH2F * histo2DChi2Dummy = new TH2F("histo2DChi2Dummy","histo2DChi2Dummy",200,0.00,50,1000,0.04,200);
    SetStyleHistoTH2ForGraphs(histo2DChi2Dummy, "#chi^{2}/NDF","#it{p}_{T} (GeV/#it{c})",0.035,0.04, 0.035,0.04, 0.98,0.9);
    histo2DChi2Dummy->GetZaxis()->SetRangeUser(1,6e2);
    canvasChi2Plots->cd();

    // RECONSTRUCTED PLOT
    histo2DChi2Dummy->Draw("copy");
    histoGammaChi2Pt->Draw("col,same");
    ex1 = new TExec("ex1","PalColor();");
    ex1->Draw();
    histoGammaChi2Pt->Draw("col,same");
    labelEnergy                  = new TLatex(0.95,0.90,Form("%s, %s",labelALICEforPlots.Data(),collisionSystem.Data()));
    SetStyleTLatex( labelEnergy, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelEnergy->Draw();
    labelProcess                     = new TLatex(0.95,0.86,isMC ? "#gamma candidates (MC rec.)" : "#gamma candidates (data)");
    SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelProcess->Draw();
    histo2DChi2Dummy->Draw("axis,same");
    canvasChi2Plots->SaveAs(Form("%s/Chi2_vs_Pt_AllRec_%s.%s",outputDir.Data(),isMC ? "MC" : "data",suffix.Data()));


    // TRUE GAMMAS PLOT
    histo2DChi2Dummy->Draw("copy");
    if (isMC){
        // draw true gamma->ee
        histoTrueMCKindChi2Pt[0]->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTrueMCKindChi2Pt[0]->Draw("col,same");
    }
    labelColor                  = new TLatex(0.12,0.90, isMC ? "#gamma rec. from e^{+}e^{-} (color)" : "");
    SetStyleTLatex( labelColor, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
    labelColor->Draw();
    labelProcess                     = new TLatex(0.95,0.86,isMC ? "#gamma candidates (MC true)" : "#gamma candidates (data)");
    SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelEnergy->Draw();
    labelProcess->Draw();
    histo2DChi2Dummy->Draw("axis,same");
    if (isMC)
        canvasChi2Plots->SaveAs(Form("%s/Chi2_vs_Pt_trueGamma_woCuts.%s",outputDir.Data(),suffix.Data()));


    // TRUE GAMMAS AND COMBINATORIAL GAMMA CANDIDATES PLOT
    histo2DChi2Dummy->Draw("copy");
    if (isMC){
        // draw true gamma->ee
        histoTrueMCKindChi2Pt[0]->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTrueMCKindChi2Pt[0]->Draw("col,same");
        // // draw ee combinatorics
        // histoTrueMCKindChi2Pt[11]->Draw("col,same");
        // ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        // ex2->Draw();
        // histoTrueMCKindChi2Pt[11]->Draw("col,same");
        // draw pipi combinatorics
        histoTrueMCKindChi2Pt[13]->Draw("col,same");
        ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        ex2->Draw();
        histoTrueMCKindChi2Pt[13]->Draw("col,same");
    }
    labelColor->Draw();
    labelBW                     = new TLatex(0.12,0.86,"#gamma rec. from #pi^{+}#pi^{-} (BAW)");
    SetStyleTLatex( labelBW, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
    labelBW->Draw();
    labelEnergy->Draw();
    labelProcess->Draw();
    histo2DChi2Dummy->Draw("axis,same");
    if (isMC)
        canvasChi2Plots->SaveAs(Form("%s/Chi2_vs_Pt_trueGammaAndCombPlus_woCuts.%s",outputDir.Data(),suffix.Data()));


    // TRUE GAMMAS AND COMBINATORIAL GAMMA CANDIDATES PLOT
    histo2DChi2Dummy->Draw("copy");
    // draw true gamma->ee
    if (isMC){
        histoTrueMCKindChi2Pt[0]->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTrueMCKindChi2Pt[0]->Draw("col,same");
        // // draw ee combinatorics
        // histoTrueMCKindChi2Pt[11]->Draw("col,same");
        // ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        // ex2->Draw();
        // histoTrueMCKindChi2Pt[11]->Draw("col,same");
        // // draw pipi combinatorics
        // histoTrueMCKindChi2Pt[13]->Draw("col,same");
        // ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        // ex2->Draw();
        // histoTrueMCKindChi2Pt[13]->Draw("col,same");
    } else {
        histoGammaChi2Pt->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoGammaChi2Pt->Draw("col,same");
    }
    // TF1 *funcPtDepAlphaCut_std = new TF1("funcPtDepAlphaCut_std","[0]*TMath::Exp(2*TMath::Abs(x))",-1,1.);
    // DrawGammaSetMarkerTF1( funcPtDepAlphaCut_std, 1, 3, kMagenta+2);
    // funcPtDepAlphaCut_std->SetParameter(0,0.07);
    // funcPtDepAlphaCut_std->Draw("same");
    // TF1 *funcPtDepAlphaCut_hard = new TF1("funcPtDepAlphaCut_hard","[0]*TMath::Exp(2.1*TMath::Abs(x))",-1,1.);
    // DrawGammaSetMarkerTF1( funcPtDepAlphaCut_hard, 2, 3, kMagenta+4);
    // funcPtDepAlphaCut_hard->SetParameter(0,0.07);
    // funcPtDepAlphaCut_hard->Draw("same");
    // TF1 *funcPtDepAlphaCut_soft = new TF1("funcPtDepAlphaCut_soft","[0]*TMath::Exp(1.9*TMath::Abs(x))",-1,1.);
    // DrawGammaSetMarkerTF1( funcPtDepAlphaCut_soft, 2, 3, kMagenta-4);
    // funcPtDepAlphaCut_soft->SetParameter(0,0.07);
    // funcPtDepAlphaCut_soft->Draw("same");
    // TF1 *funcPtDepAlphaCut_soft2 = new TF1("funcPtDepAlphaCut_soft2","[0]*TMath::Exp(1.9*TMath::Abs(x))",-1,1.);
    // DrawGammaSetMarkerTF1( funcPtDepAlphaCut_soft2, 2, 3, kMagenta-1);
    // funcPtDepAlphaCut_soft2->SetParameter(0,0.06);
    // funcPtDepAlphaCut_soft2->Draw("same");
    DrawGammaLines( 30.,  30., 0.04, 15, 3, kMagenta+2, 9);
    funcOldCutDummy = new TF1("funcOldCutDummy","x/[0]",0.005,1.);
    DrawGammaSetMarkerTF1( funcOldCutDummy, 9, 3, kMagenta+2);
    DrawGammaLines(40.,  40., 0.04, 15, 3, kMagenta+1, 7);
    TF1* funcOldCutDummy1 = new TF1("funcOldCutDummy","x/[0]",0.005,1.);
    DrawGammaSetMarkerTF1( funcOldCutDummy1, 7, 3, kMagenta+1);
    TLegend* legendChi2PlotFits  = GetAndSetLegend2(0.74, 0.75, 0.95, 0.75+(0.035*2*1.35), 0.85*textSizeLabelsPixel);
    legendChi2PlotFits->AddEntry(funcOldCutDummy, "#chi^{2}_{max} = 30","l");
    legendChi2PlotFits->AddEntry(funcOldCutDummy1, "#chi^{2}_{max} = 40","l");
    legendChi2PlotFits->Draw();
    // TLatex* labelHighpTSignal    = new TLatex(0.90,0.78,"high #it{p}_{T} signal #rightarrow");
    // SetStyleTLatex( labelHighpTSignal, 0.95*textSizeLabelsPixel,4,kBlue+2,43,kTRUE,31);
    // labelHighpTSignal->Draw();
    labelColor->Draw();
    // labelBW->Draw();
    labelEnergy->Draw();
    labelProcess->Draw();
    histo2DChi2Dummy->Draw("axis,same");
    canvasChi2Plots->SaveAs(Form("%s/Chi2_vs_Pt_Final_withCuts_%s.%s",outputDir.Data(),isMC ? "MC" : "data",suffix.Data()));



    //   _____   _____ _____   _____        _____ _____
    //  |  __ \ / ____|_   _| |  __ \ /\   |_   _|  __ \
    //  | |__) | (___   | |   | |__) /  \    | | | |__) |
    //  |  ___/ \___ \  | |   |  ___/ /\ \   | | |  _  /
    //  | |     ____) |_| |_  | |  / ____ \ _| |_| | \ \
    //  |_|    |_____/|_____| |_| /_/    \_\_____|_|  \_\



    textSizeLabelsPixel                   = 1200*0.04;
    TCanvas* canvasPsiPairPlots = new TCanvas("canvasPsiPairPlots","",200,10,1350,1200);  // gives the page size
    DrawGammaCanvasSettings( canvasPsiPairPlots, 0.08, 0.02, 0.02, 0.09);
    canvasPsiPairPlots->SetLogy();
    canvasPsiPairPlots->SetLogz();
    TH2F * histo2DPsiPairDummy = new TH2F("histo2DPsiPairDummy","histo2DPsiPairDummy",200,-0.16,0.16,1000,0.04,200);
    SetStyleHistoTH2ForGraphs(histo2DPsiPairDummy, "#Psi_{pair}","#it{p}_{T} (GeV/#it{c})",0.035,0.04, 0.035,0.04, 0.98,0.9);
    histo2DPsiPairDummy->GetZaxis()->SetRangeUser(1,6e2);
    canvasPsiPairPlots->cd();

    // RECONSTRUCTED PLOT
    histo2DPsiPairDummy->Draw("copy");
    histoGammaPsiPairPt->Draw("col,same");
    ex1 = new TExec("ex1","PalColor();");
    ex1->Draw();
    histoGammaPsiPairPt->Draw("col,same");
    labelEnergy                  = new TLatex(0.95,0.90,Form("%s, %s",labelALICEforPlots.Data(),collisionSystem.Data()));
    SetStyleTLatex( labelEnergy, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelEnergy->Draw();
    labelProcess                     = new TLatex(0.95,0.86,isMC ? "#gamma candidates (MC rec.)" : "#gamma candidates (data)");
    SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelProcess->Draw();
    histo2DPsiPairDummy->Draw("axis,same");
    canvasPsiPairPlots->SaveAs(Form("%s/PsiPair_vs_Pt_AllRec_%s.%s",outputDir.Data(),isMC ? "MC" : "data",suffix.Data()));


    // TRUE GAMMAS PLOT
    histo2DPsiPairDummy->Draw("copy");
    if (isMC){
        // draw true gamma->ee
        histoTrueMCKindPsiPairPt[0]->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTrueMCKindPsiPairPt[0]->Draw("col,same");
    }
    labelColor                  = new TLatex(0.12,0.90, isMC ? "#gamma rec. from e^{+}e^{-} (color)" : "");
    SetStyleTLatex( labelColor, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
    labelColor->Draw();
    labelProcess                     = new TLatex(0.95,0.86,isMC ? "#gamma candidates (MC true)" : "#gamma candidates (data)");
    SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelEnergy->Draw();
    labelProcess->Draw();
    histo2DPsiPairDummy->Draw("axis,same");
    if (isMC)
        canvasPsiPairPlots->SaveAs(Form("%s/PsiPair_vs_Pt_trueGamma_woCuts.%s",outputDir.Data(),suffix.Data()));


    // TRUE GAMMAS AND COMBINATORIAL GAMMA CANDIDATES PLOT
    histo2DPsiPairDummy->Draw("copy");
    if (isMC){
        // draw true gamma->ee
        histoTrueMCKindPsiPairPt[0]->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTrueMCKindPsiPairPt[0]->Draw("col,same");
        // // draw ee combinatorics
        // histoTrueMCKindPsiPairPt[11]->Draw("col,same");
        // ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        // ex2->Draw();
        // histoTrueMCKindPsiPairPt[11]->Draw("col,same");
        // draw pipi combinatorics
        histoTrueMCKindPsiPairPt[13]->Draw("col,same");
        ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        ex2->Draw();
        histoTrueMCKindPsiPairPt[13]->Draw("col,same");
    }
    labelColor->Draw();
    labelBW                     = new TLatex(0.12,0.86,"#gamma rec. from #pi^{+}#pi^{-} (BAW)");
    SetStyleTLatex( labelBW, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
    labelBW->Draw();
    labelEnergy->Draw();
    labelProcess->Draw();
    histo2DPsiPairDummy->Draw("axis,same");
    if (isMC)
        canvasPsiPairPlots->SaveAs(Form("%s/PsiPair_vs_Pt_trueGammaAndCombPlus_woCuts.%s",outputDir.Data(),suffix.Data()));


    // TRUE GAMMAS AND COMBINATORIAL GAMMA CANDIDATES PLOT
    histo2DPsiPairDummy->Draw("copy");
    if (isMC){
        // draw true gamma->ee
        histoTrueMCKindPsiPairPt[0]->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTrueMCKindPsiPairPt[0]->Draw("col,same");
        // // draw ee combinatorics
        // histoTrueMCKindPsiPairPt[11]->Draw("col,same");
        // ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        // ex2->Draw();
        // histoTrueMCKindPsiPairPt[11]->Draw("col,same");
        // // draw pipi combinatorics
        // histoTrueMCKindPsiPairPt[13]->Draw("col,same");
        // ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        // ex2->Draw();
        // histoTrueMCKindPsiPairPt[13]->Draw("col,same");
    } else {
        histoGammaPsiPairPt->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoGammaPsiPairPt->Draw("col,same");
    }

    // TF1 *funcPtDepAlphaCut_std = new TF1("funcPtDepAlphaCut_std","[0]*TMath::Exp(2*TMath::Abs(x))",-1,1.);
    // DrawGammaSetMarkerTF1( funcPtDepAlphaCut_std, 1, 3, kMagenta+2);
    // funcPtDepAlphaCut_std->SetParameter(0,0.07);
    // funcPtDepAlphaCut_std->Draw("same");
    // TF1 *funcPtDepAlphaCut_hard = new TF1("funcPtDepAlphaCut_hard","[0]*TMath::Exp(2.1*TMath::Abs(x))",-1,1.);
    // DrawGammaSetMarkerTF1( funcPtDepAlphaCut_hard, 2, 3, kMagenta+4);
    // funcPtDepAlphaCut_hard->SetParameter(0,0.07);
    // funcPtDepAlphaCut_hard->Draw("same");
    // TF1 *funcPtDepAlphaCut_soft = new TF1("funcPtDepAlphaCut_soft","[0]*TMath::Exp(1.9*TMath::Abs(x))",-1,1.);
    // DrawGammaSetMarkerTF1( funcPtDepAlphaCut_soft, 2, 3, kMagenta-4);
    // funcPtDepAlphaCut_soft->SetParameter(0,0.07);
    // funcPtDepAlphaCut_soft->Draw("same");
    // TF1 *funcPtDepAlphaCut_soft2 = new TF1("funcPtDepAlphaCut_soft2","[0]*TMath::Exp(1.9*TMath::Abs(x))",-1,1.);
    // DrawGammaSetMarkerTF1( funcPtDepAlphaCut_soft2, 2, 3, kMagenta-1);
    // funcPtDepAlphaCut_soft2->SetParameter(0,0.06);
    // funcPtDepAlphaCut_soft2->Draw("same");
    DrawGammaLines( 0.1,  0.1, 0.04, 15, 3, kMagenta+2, 9);
    DrawGammaLines( -0.1,  -0.1, 0.04, 15, 3, kMagenta+2, 9);
    funcOldCutDummy = new TF1("funcOldCutDummy","x/[0]",0.005,1.);
    DrawGammaSetMarkerTF1( funcOldCutDummy, 9, 3, kMagenta+2);
    DrawGammaLines(0.15,  0.15, 0.04, 15, 3, kMagenta+1, 7);
    DrawGammaLines(-0.15,  -0.15, 0.04, 15, 3, kMagenta+1, 7);
    funcOldCutDummy1 = new TF1("funcOldCutDummy","x/[0]",0.005,1.);
    DrawGammaSetMarkerTF1( funcOldCutDummy1, 7, 3, kMagenta+1);
    TLegend* legendPsiPairPlotFits  = GetAndSetLegend2(0.74, 0.75, 0.95, 0.75+(0.035*2*1.35), 0.85*textSizeLabelsPixel);
    legendPsiPairPlotFits->AddEntry(funcOldCutDummy, "|#Psi_{pair}| = 0.1","l");
    legendPsiPairPlotFits->AddEntry(funcOldCutDummy1, "|#Psi_{pair}| = 0.15","l");
    legendPsiPairPlotFits->Draw();
    // TLatex* labelHighpTSignal    = new TLatex(0.90,0.78,"high #it{p}_{T} signal #rightarrow");
    // SetStyleTLatex( labelHighpTSignal, 0.95*textSizeLabelsPixel,4,kBlue+2,43,kTRUE,31);
    // labelHighpTSignal->Draw();
    labelColor->Draw();
    // labelBW->Draw();
    labelEnergy->Draw();
    labelProcess->Draw();
    histo2DPsiPairDummy->Draw("axis,same");
    canvasPsiPairPlots->SaveAs(Form("%s/PsiPair_vs_Pt_Final_withCuts_%s.%s",outputDir.Data(),isMC ? "MC" : "data",suffix.Data()));




    //   _____   _____ _____   _____        _____ _____         _____ _    _ _____ ___
    //  |  __ \ / ____|_   _| |  __ \ /\   |_   _|  __ \       / ____| |  | |_   _|__ \
    //  | |__) | (___   | |   | |__) /  \    | | | |__) | ___ | |    | |__| | | |    ) |
    //  |  ___/ \___ \  | |   |  ___/ /\ \   | | |  _  /  ___ | |    |  __  | | |   / /
    //  | |     ____) |_| |_  | |  / ____ \ _| |_| | \ \      | |____| |  | |_| |_ / /_
    //  |_|    |_____/|_____| |_| /_/    \_\_____|_|  \_\      \_____|_|  |_|_____|____|

    textSizeLabelsPixel                   = 1200*0.04;
    TCanvas* canvasChi2PsiPairPlots = new TCanvas("canvasChi2PsiPairPlots","",200,10,1350,1200);  // gives the page size
    DrawGammaCanvasSettings( canvasChi2PsiPairPlots, 0.1, 0.02, 0.02, 0.09);
    // canvasChi2PsiPairPlots->SetLogy();
    canvasChi2PsiPairPlots->SetLogz();
    TH2F * histo2DChi2PsiPairDummy = new TH2F("histo2DChi2PsiPairDummy","histo2DChi2PsiPairDummy",200,0.0,50,200,-0.215,0.215);
    SetStyleHistoTH2ForGraphs(histo2DChi2PsiPairDummy, "#chi^{2}/NDF", "#Psi_{pair}",0.035,0.04, 0.035,0.04, 0.98,1.2);
    histo2DChi2PsiPairDummy->GetZaxis()->SetRangeUser(8,6e2);
    canvasChi2PsiPairPlots->cd();

    // RECONSTRUCTED PLOT
    histo2DChi2PsiPairDummy->Draw("copy");
    histoGammaChi2PsiPair->Draw("col,same");
    ex1 = new TExec("ex1","PalColor();");
    ex1->Draw();
    histoGammaChi2PsiPair->Draw("col,same");
    labelEnergy                  = new TLatex(0.95,0.90,Form("%s, %s",labelALICEforPlots.Data(),collisionSystem.Data()));
    SetStyleTLatex( labelEnergy, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelEnergy->Draw();
    labelProcess                     = new TLatex(0.95,0.86,isMC ? "#gamma candidates (MC rec.)" : "#gamma candidates (data)");
    SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelProcess->Draw();
    histo2DChi2PsiPairDummy->Draw("axis,same");
    canvasChi2PsiPairPlots->SaveAs(Form("%s/Chi2PsiPair_vs_Pt_AllRec_%s.%s",outputDir.Data(),isMC ? "MC" : "data",suffix.Data()));

    for(Int_t bin=0; bin<15; bin++){
        histo2DChi2PsiPairDummy->GetZaxis()->SetRangeUser(1,6e2);
        histo2DChi2PsiPairDummy->Draw("copy");
        histoGammaChi2PsiPairPtSliced[bin]->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoGammaChi2PsiPairPtSliced[bin]->Draw("col,same");
        labelpTrange    = new TLatex(0.95,0.82,Form("%2.1f < #it{p}_{T} < %2.1f GeV/#it{c}",projBinnin[bin],projBinnin[bin+1]));
        SetStyleTLatex( labelpTrange, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
        labelpTrange->Draw();
        labelEnergy->Draw();
        labelProcess->Draw();
        histo2DChi2PsiPairDummy->Draw("axis,same");
        canvasChi2PsiPairPlots->SaveAs(Form("%s/Chi2PsiPair/%dbin_Chi2PsiPair_vs_Pt_AllRec_%s.%s",outputDir.Data(),bin,isMC ? "MC" : "data",suffix.Data()));
    }
    histo2DChi2PsiPairDummy->GetZaxis()->SetRangeUser(6,6e2);
    // TRUE GAMMAS PLOT
    histo2DChi2PsiPairDummy->Draw("copy");
    // draw true gamma->ee
    if (isMC){
        histoTrueGammaChi2PsiPair->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTrueGammaChi2PsiPair->Draw("col,same");
    }
    labelColor                  = new TLatex(0.12,0.90, isMC ? "#gamma rec. from e^{+}e^{-} (color)" : "");
    SetStyleTLatex( labelColor, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
    labelColor->Draw();
    labelProcess                     = new TLatex(0.95,0.86,isMC ? "#gamma candidates (MC true)" : "#gamma candidates (data)");
    SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelEnergy->Draw();
    labelProcess->Draw();
    histo2DChi2PsiPairDummy->Draw("axis,same");
    if (isMC)
        canvasChi2PsiPairPlots->SaveAs(Form("%s/Chi2PsiPair_vs_Pt_trueGamma_woCuts.%s",outputDir.Data(),suffix.Data()));

    histo2DChi2PsiPairDummy->GetZaxis()->SetRangeUser(1,6e2);
    for(Int_t bin=0; bin<15; bin++){
        if (!isMC) continue;
        histo2DChi2PsiPairDummy->Draw("copy");
        // draw true gamma->ee
        histoTrueMCKindChi2PsiPairPtSliced[0][bin]->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTrueMCKindChi2PsiPairPtSliced[0][bin]->Draw("col,same");
        labelColor->Draw();
        labelEnergy->Draw();
        labelProcess->Draw();
        labelpTrange    = new TLatex(0.95,0.82,Form("%2.1f < #it{p}_{T} < %2.1f GeV/#it{c}",projBinnin[bin],projBinnin[bin+1]));
        SetStyleTLatex( labelpTrange, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
        labelpTrange->Draw();
        histo2DChi2PsiPairDummy->Draw("axis,same");
        canvasChi2PsiPairPlots->SaveAs(Form("%s/Chi2PsiPair/%dbin_Chi2PsiPair_vs_Pt_trueGamma.%s",outputDir.Data(),bin,suffix.Data()));
    }

    labelBW                     = new TLatex(0.12,0.86,"#gamma rec. from e^{#pm}#pi^{#pm} (BAW)");
    SetStyleTLatex( labelBW, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);

    for(Int_t bin=0; bin<15; bin++){
        if (!isMC) continue;
        histo2DChi2PsiPairDummy->Draw("copy");
        // draw true gamma->ee
        histoTrueMCKindChi2PsiPairPtSliced[0][bin]->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTrueMCKindChi2PsiPairPtSliced[0][bin]->Draw("col,same");
        // draw ee combinatorics
        histoTrueMCKindChi2PsiPairPtSliced[11][bin]->Draw("col,same");
        ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        ex2->Draw();
        histoTrueMCKindChi2PsiPairPtSliced[11][bin]->Draw("col,same");
        labelColor->Draw();
        labelBW->Draw();
        labelEnergy->Draw();
        labelProcess->Draw();
        labelpTrange    = new TLatex(0.95,0.82,Form("%2.1f < #it{p}_{T} < %2.1f GeV/#it{c}",projBinnin[bin],projBinnin[bin+1]));
        SetStyleTLatex( labelpTrange, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
        labelpTrange->Draw();
        histo2DChi2PsiPairDummy->Draw("axis,same");
        canvasChi2PsiPairPlots->SaveAs(Form("%s/Chi2PsiPair/%dbin_Chi2PsiPair_vs_Pt_trueGammaAndComb.%s",outputDir.Data(),bin,suffix.Data()));
    }


    histo2DChi2PsiPairDummy->GetZaxis()->SetRangeUser(6,6e2);
    // TRUE GAMMAS AND COMBINATORICS PLOT
    histo2DChi2PsiPairDummy->Draw("copy");
    if (isMC){
        // draw true gamma->ee
        histoTrueGammaChi2PsiPair->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTrueGammaChi2PsiPair->Draw("col,same");
        // draw ee combinatorics
        // histoTrueMCKindChi2PsiPair[11]->Draw("col,same");
        // ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        // ex2->Draw();
        // histoTrueMCKindChi2PsiPair[11]->Draw("col,same");
        // draw pipi combinatorics
        histoTrueMCKindChi2PsiPair[13]->Draw("col,same");
        ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        ex2->Draw();
        histoTrueMCKindChi2PsiPair[13]->Draw("col,same");
    }
    labelColor                  = new TLatex(0.12,0.90, isMC ? "#gamma rec. from e^{+}e^{-} (color)" : "");
    SetStyleTLatex( labelColor, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
    labelColor->Draw();
    labelBW                     = new TLatex(0.12,0.86,"#gamma rec. from e^{#pm}#pi^{#pm} (BAW)");
    SetStyleTLatex( labelBW, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
    labelBW->Draw();
    labelProcess                     = new TLatex(0.95,0.86,isMC ? "#gamma candidates (MC true)" : "#gamma candidates (data)");
    SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelEnergy->Draw();
    labelProcess->Draw();
    histo2DChi2PsiPairDummy->Draw("axis,same");
    if (isMC)
        canvasChi2PsiPairPlots->SaveAs(Form("%s/Chi2PsiPair_vs_Pt_trueGammaAndComb_woCuts.%s",outputDir.Data(),suffix.Data()));


    // TRUE GAMMAS AND COMBINATORIAL GAMMA CANDIDATES PLOT
    histo2DChi2PsiPairDummy->Draw("copy");
    if (isMC){
        // draw true gamma->ee
        histoTrueGammaChi2PsiPair->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTrueGammaChi2PsiPair->Draw("col,same");
        // // draw ee combinatorics
        // histoTrueMCKindChi2PsiPairPt[11]->Draw("col,same");
        // ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        // ex2->Draw();
        // histoTrueMCKindChi2PsiPairPt[11]->Draw("col,same");
        // // draw pipi combinatorics
        // histoTrueMCKindChi2PsiPairPt[13]->Draw("col,same");
        // ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        // ex2->Draw();
        // histoTrueMCKindChi2PsiPairPt[13]->Draw("col,same");
    } else {
        histoGammaChi2PsiPair->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoGammaChi2PsiPair->Draw("col,same");
    }
    TF1 *funcPtDepAlphaCut_std = new TF1("funcPtDepAlphaCut_std","-0.1+(0.1/30)*x",0,30);
    DrawGammaSetMarkerTF1( funcPtDepAlphaCut_std, 1, 3, kMagenta+2);
    funcPtDepAlphaCut_std->Draw("same");
    TF1 *funcPtDepAlphaCut_std2 = new TF1("funcPtDepAlphaCut_std","0.1-(0.1/30)*x",0,30);
    DrawGammaSetMarkerTF1( funcPtDepAlphaCut_std2, 1, 3, kMagenta+2);
    funcPtDepAlphaCut_std2->Draw("same");
    TF1 *funcPtDepAlphaCut_hard = new TF1("funcPtDepAlphaCut_hard","-0.15+0.003*x",0,50);
    DrawGammaSetMarkerTF1( funcPtDepAlphaCut_hard, 2, 3, kMagenta+4);
    funcPtDepAlphaCut_hard->Draw("same");
    TF1 *funcPtDepAlphaCut_hard2 = new TF1("funcPtDepAlphaCut_hard","0.15-0.003*x",0,50);
    DrawGammaSetMarkerTF1( funcPtDepAlphaCut_hard2, 2, 3, kMagenta+4);
    funcPtDepAlphaCut_hard2->Draw("same");
    TF1 *funcPtDepAlphaCut_soft = new TF1("funcPtDepAlphaCut_soft","0.18*TMath::Exp(-0.055*x)",0,50);
    DrawGammaSetMarkerTF1( funcPtDepAlphaCut_soft, 2, 3, kMagenta-4);
    funcPtDepAlphaCut_soft->Draw("same");
    TF1 *funcPtDepAlphaCut_soft2 = new TF1("funcPtDepAlphaCut_soft","-(0.18*TMath::Exp(-0.055*x))",0,50);
    DrawGammaSetMarkerTF1( funcPtDepAlphaCut_soft2, 2, 3, kMagenta-4);
    funcPtDepAlphaCut_soft2->Draw("same");
    TF1 *funcPtDepAlphaCut_verysoft = new TF1("funcPtDepAlphaCut_verysoft","0.20*TMath::Exp(-0.050*x)",0,50);
    DrawGammaSetMarkerTF1( funcPtDepAlphaCut_verysoft, 2, 3, kMagenta-8);
    funcPtDepAlphaCut_verysoft->Draw("same");
    TF1 *funcPtDepAlphaCut_verysoft2 = new TF1("funcPtDepAlphaCut_verysoft2","-(0.20*TMath::Exp(-0.050*x))",0,50);
    DrawGammaSetMarkerTF1( funcPtDepAlphaCut_verysoft2, 2, 3, kMagenta-8);
    funcPtDepAlphaCut_verysoft2->Draw("same");
    TF1 *funcPtDepAlphaCut_veryhard = new TF1("funcPtDepAlphaCut_veryhard","0.15*TMath::Exp(-0.065*x)",0,50);
    DrawGammaSetMarkerTF1( funcPtDepAlphaCut_veryhard, 2, 3, kMagenta+1);
    funcPtDepAlphaCut_veryhard->Draw("same");
    TF1 *funcPtDepAlphaCut_veryhard2 = new TF1("funcPtDepAlphaCut_veryhard2","-(0.15*TMath::Exp(-0.065*x))",0,50);
    DrawGammaSetMarkerTF1( funcPtDepAlphaCut_veryhard2, 2, 3, kMagenta+1);
    funcPtDepAlphaCut_veryhard2->Draw("same");
    TLegend* legendChi2PsiPairPlotFits  = GetAndSetLegend2(0.50, 0.13, 0.95, 0.13+(0.035*5*1.35), 0.85*textSizeLabelsPixel);
    legendChi2PsiPairPlotFits->AddEntry(funcPtDepAlphaCut_std, "|#Psi_{pair}| < 0.10/30#chi^{2}+0.10","l");
    legendChi2PsiPairPlotFits->AddEntry(funcPtDepAlphaCut_hard, "|#Psi_{pair}| < 0.15/50#chi^{2}+0.15","l");
    legendChi2PsiPairPlotFits->AddEntry(funcPtDepAlphaCut_veryhard, "|#Psi_{pair}| < 0.15*exp(-0.065#chi^{2})","l");
    legendChi2PsiPairPlotFits->AddEntry(funcPtDepAlphaCut_soft, "|#Psi_{pair}| < 0.18*exp(-0.055#chi^{2})","l");
    legendChi2PsiPairPlotFits->AddEntry(funcPtDepAlphaCut_verysoft2, "|#Psi_{pair}| < 0.20*exp(-0.050#chi^{2})","l");
    legendChi2PsiPairPlotFits->Draw();
    labelEnergy->Draw();
    labelProcess->Draw();
    histo2DChi2PsiPairDummy->Draw("axis,same");
    canvasChi2PsiPairPlots->SaveAs(Form("%s/Chi2PsiPair_vs_Pt_Final_withCuts_%s.%s",outputDir.Data(),isMC ? "MC" : "data",suffix.Data()));





    //    ____ _______                   _______     ____  __ __  __ ______ _______ _______     __
    //   / __ \__   __|           /\    / ____\ \   / /  \/  |  \/  |  ____|__   __|  __ \ \   / /
    //  | |  | | | |             /  \  | (___  \ \_/ /| \  / | \  / | |__     | |  | |__) \ \_/ /
    //  | |  | | | |            / /\ \  \___ \  \   / | |\/| | |\/| |  __|    | |  |  _  / \   /
    //  | |__| | | |           / ____ \ ____) |  | |  | |  | | |  | | |____   | |  | | \ \  | |
    //   \___\_\ |_|          /_/    \_\_____/   |_|  |_|  |_|_|  |_|______|  |_|  |_|  \_\ |_|
    //


    textSizeLabelsPixel                   = 1200*0.04;
    TCanvas* canvasAlphaQtPlots = new TCanvas("canvasAlphaQtPlots","",200,10,1350,1200);  // gives the page size
    DrawGammaCanvasSettings( canvasAlphaQtPlots, 0.1, 0.02, 0.02, 0.09);
    // canvasAlphaQtPlots->SetLogy();
    canvasAlphaQtPlots->SetLogz();
    TH2F * histo2DAlphaQtDummy = new TH2F("histo2DAlphaQtDummy","histo2DAlphaQtDummy",200,-1.05,1.05,200,0.,0.15);
    SetStyleHistoTH2ForGraphs(histo2DAlphaQtDummy, "#alpha^{#gamma} = (#it{p}^{+}_{L}-#it{p}^{-}_{L})/(#it{p}^{+}_{L}+#it{p}^{-}_{L})", "#it{q}_{T}^{#gamma}",0.035,0.04, 0.035,0.04, 0.98,1.2);
    histo2DAlphaQtDummy->GetZaxis()->SetRangeUser(1,6e2);
    canvasAlphaQtPlots->cd();

    // RECONSTRUCTED PLOT
    histo2DAlphaQtDummy->Draw("copy");
    histoGammaAlphaQt->Draw("col,same");
    ex1 = new TExec("ex1","PalColor();");
    ex1->Draw();
    histoGammaAlphaQt->Draw("col,same");
    labelEnergy                  = new TLatex(0.95,0.90,Form("%s, %s",labelALICEforPlots.Data(),collisionSystem.Data()));
    SetStyleTLatex( labelEnergy, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelEnergy->Draw();
    labelProcess                     = new TLatex(0.95,0.86,isMC ? "#gamma candidates (MC rec.)" : "#gamma candidates (data)");
    SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelProcess->Draw();
    histo2DAlphaQtDummy->Draw("axis,same");
    canvasAlphaQtPlots->SaveAs(Form("%s/AlphaQt_vs_Pt_AllRec_%s.%s",outputDir.Data(),isMC ? "MC" : "data",suffix.Data()));

    for(Int_t bin=0; bin<15; bin++){
        // RECONSTRUCTED PLOT
        histo2DAlphaQtDummy->Draw("copy");
        histoGammaAlphaQtPtSliced[bin]->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoGammaAlphaQtPtSliced[bin]->Draw("col,same");
        labelEnergy->Draw();
        labelProcess->Draw();
        labelpTrange                     = new TLatex(0.95,0.82,Form("%2.1f < #it{p}_{T} < %2.1f GeV/#it{c}",projBinnin[bin],projBinnin[bin+1]));
        SetStyleTLatex( labelpTrange, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
        labelpTrange->Draw();
        histo2DAlphaQtDummy->Draw("axis,same");
        canvasAlphaQtPlots->SaveAs(Form("%s/AlphaQt/%dbin_AlphaQt_vs_Pt_AllRec_%s.%s",outputDir.Data(), bin,isMC ? "MC" : "data",suffix.Data()));
    }

    // TRUE GAMMAS PLOT
    histo2DAlphaQtDummy->Draw("copy");
    if (isMC){
        // draw true gamma->ee
        histoTruePrimGammaAlphaQt->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTruePrimGammaAlphaQt->Draw("col,same");
    }
    labelColor                  = new TLatex(0.12,0.90, isMC ? "#gamma rec. from e^{+}e^{-} (color)" : "");
    SetStyleTLatex( labelColor, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);
    labelColor->Draw();
    labelProcess                     = new TLatex(0.95,0.86,isMC ? "#gamma candidates (MC true)" : "#gamma candidates (data)");
    SetStyleTLatex( labelProcess, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
    labelEnergy->Draw();
    labelProcess->Draw();
    histo2DAlphaQtDummy->Draw("axis,same");
    if (isMC)
        canvasAlphaQtPlots->SaveAs(Form("%s/AlphaQt_vs_Pt_trueGamma_woCuts.%s",outputDir.Data(),suffix.Data()));

    labelBW                     = new TLatex(0.12,0.86,"#gamma rec. from e^{#pm}#pi^{#pm} (BAW)");
    SetStyleTLatex( labelBW, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,11);

    for(Int_t bin=0; bin<15; bin++){
        if (!isMC) continue;
        // RECONSTRUCTED PLOT
        histo2DAlphaQtDummy->Draw("copy");
        histoTrueMCKindAlphaQtPtSliced[0][bin]->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTrueMCKindAlphaQtPtSliced[0][bin]->Draw("col,same");
        labelColor->Draw();
        labelEnergy->Draw();
        labelProcess->Draw();
        labelpTrange    = new TLatex(0.95,0.82,Form("%2.1f < #it{p}_{T} < %2.1f GeV/#it{c}",projBinnin[bin],projBinnin[bin+1]));
        SetStyleTLatex( labelpTrange, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
        labelpTrange->Draw();
        histo2DAlphaQtDummy->Draw("axis,same");
        canvasAlphaQtPlots->SaveAs(Form("%s/AlphaQt/%dbin_AlphaQt_vs_Pt_trueGamma_%s.%s",outputDir.Data(), bin,isMC ? "MC" : "data",suffix.Data()));
    }
    for(Int_t bin=0; bin<15; bin++){
        if (!isMC) continue;
        // RECONSTRUCTED PLOT
        histo2DAlphaQtDummy->Draw("copy");
        histoTrueMCKindAlphaQtPtSliced[0][bin]->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTrueMCKindAlphaQtPtSliced[0][bin]->Draw("col,same");
        histoTrueMCKindAlphaQtPtSliced[11][bin]->Draw("col,same");
        ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        ex2->Draw();
        histoTrueMCKindAlphaQtPtSliced[11][bin]->Draw("col,same");
        labelEnergy->Draw();
        labelProcess->Draw();
        labelColor->Draw();
        labelBW->Draw();
        labelpTrange    = new TLatex(0.95,0.82,Form("%2.1f < #it{p}_{T} < %2.1f GeV/#it{c}",projBinnin[bin],projBinnin[bin+1]));
        SetStyleTLatex( labelpTrange, 0.95*textSizeLabelsPixel,4,1,43,kTRUE,31);
        labelpTrange->Draw();
        histo2DAlphaQtDummy->Draw("axis,same");
        canvasAlphaQtPlots->SaveAs(Form("%s/AlphaQt/%dbin_AlphaQt_vs_Pt_trueGammaAndComb_%s.%s",outputDir.Data(), bin,isMC ? "MC" : "data",suffix.Data()));
    }

    // TRUE GAMMAS AND COMBINATORIAL GAMMA CANDIDATES PLOT
    histo2DAlphaQtDummy->Draw("copy");
    if (isMC){
        // draw true gamma->ee
        histoTruePrimGammaAlphaQt->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTruePrimGammaAlphaQt->Draw("col,same");
        // // draw ee combinatorics
        // histoTrueMCKindAlphaQtPt[11]->Draw("col,same");
        // ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        // ex2->Draw();
        // histoTrueMCKindAlphaQtPt[11]->Draw("col,same");
        // // draw pipi combinatorics
        // histoTrueMCKindAlphaQtPt[13]->Draw("col,same");
        // ex2 = new TExec("ex2","gStyle->SetPalette(9);");
        // ex2->Draw();
        // histoTrueMCKindAlphaQtPt[13]->Draw("col,same");
    } else {
        histoGammaAlphaQt->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoGammaAlphaQt->Draw("col,same");
    }
    // TF1 *funcPtDepAlphaQtCut_std = new TF1("funcPtDepAlphaQtCut_std","-0.1+(0.1/30)*x",0,30);
    // DrawGammaSetMarkerTF1( funcPtDepAlphaQtCut_std, 1, 3, kMagenta+2);
    // funcPtDepAlphaQtCut_std->Draw("same");
    // TF1 *funcPtDepAlphaQtCut_std2 = new TF1("funcPtDepAlphaQtCut_std","0.1-(0.1/30)*x",0,30);
    // DrawGammaSetMarkerTF1( funcPtDepAlphaQtCut_std2, 1, 3, kMagenta+2);
    // funcPtDepAlphaQtCut_std2->Draw("same");
    // TF1 *funcPtDepAlphaQtCut_hard = new TF1("funcPtDepAlphaQtCut_hard","-0.15+0.003*x",0,50);
    // DrawGammaSetMarkerTF1( funcPtDepAlphaQtCut_hard, 2, 3, kMagenta+4);
    // funcPtDepAlphaQtCut_hard->Draw("same");
    // TF1 *funcPtDepAlphaQtCut_hard2 = new TF1("funcPtDepAlphaQtCut_hard","0.15-0.003*x",0,50);
    // DrawGammaSetMarkerTF1( funcPtDepAlphaQtCut_hard2, 2, 3, kMagenta+4);
    // funcPtDepAlphaQtCut_hard2->Draw("same");
    // TF1 *funcPtDepAlphaQtCut_soft = new TF1("funcPtDepAlphaQtCut_soft","0.18*TMath::Exp(-0.055*x)",0,50);
    // DrawGammaSetMarkerTF1( funcPtDepAlphaQtCut_soft, 2, 3, kMagenta-4);
    // funcPtDepAlphaQtCut_soft->Draw("same");
    // TF1 *funcPtDepAlphaQtCut_soft2 = new TF1("funcPtDepAlphaQtCut_soft","-(0.18*TMath::Exp(-0.055*x))",0,50);
    // DrawGammaSetMarkerTF1( funcPtDepAlphaQtCut_soft2, 2, 3, kMagenta-4);
    // funcPtDepAlphaQtCut_soft2->Draw("same");
    // TLegend* legendAlphaQtPlotFits  = GetAndSetLegend2(0.50, 0.15, 0.95, 0.15+(0.035*3*1.35), 0.85*textSizeLabelsPixel);
    // legendAlphaQtPlotFits->AddEntry(funcPtDepAlphaQtCut_std, "|#Psi_{pair}| < -0.10/30#chi^{2}+0.10","l");
    // legendAlphaQtPlotFits->AddEntry(funcPtDepAlphaQtCut_hard, "|#Psi_{pair}| < -0.15/50#chi^{2}+0.15","l");
    // legendAlphaQtPlotFits->AddEntry(funcPtDepAlphaQtCut_soft, "|#Psi_{pair}| < -0.18*exp(-0.055#chi^{2})","l");
    // legendAlphaQtPlotFits->Draw();
    // TLatex* labelHighpTSignal    = new TLatex(0.90,0.78,"high #it{p}_{T} signal #rightarrow");
    // SetStyleTLatex( labelHighpTSignal, 0.95*textSizeLabelsPixel,4,kBlue+2,43,kTRUE,31);
    // labelHighpTSignal->Draw();
    // labelColor->Draw();
    // labelBW->Draw();
    labelEnergy->Draw();
    labelProcess->Draw();
    histo2DAlphaQtDummy->Draw("axis,same");
    canvasAlphaQtPlots->SaveAs(Form("%s/AlphaQt_vs_Pt_Final_withCuts_%s.%s",outputDir.Data(),isMC ? "MC" : "data",suffix.Data()));





}
