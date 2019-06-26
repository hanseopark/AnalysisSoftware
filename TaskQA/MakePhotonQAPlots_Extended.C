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
#include "TASImage.h"
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
    TString fileName        = "GammaConvV1_QA_5460001022092970003190000000.root",
    TString cutnumber       = "0005314140",
    Bool_t isMC             = kFALSE,
    TString suffix          = "pdf",
    TString optionEnergy    = "8TeV",
    TString centCutNumber   = "000",
    TString addName         = ""
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
    TString nameDirectory           = Form("GammaConvV1_QA_%s",  cutnumber.Data());
    TDirectory* directoryConv       = (TDirectory*)f->Get(nameDirectory.Data());
    if (!directoryConv){
        cout << "missing directory: " << nameDirectory.Data() << ", aborting ..."<< endl;
        return;
    }
  //___________________________________ Declaration of files _____________________________________________
    TString dateForOutput                           = ReturnDateStringForOutput();
    cout << dateForOutput.Data() << endl;
    TString collisionSystem                         = ReturnFullCollisionsSystem(optionEnergy);
    TString centralityString                        = GetCentralityString(centCutNumber);
    if (centralityString.CompareTo("pp")==0){
        centralityString    = "";
    } else {
        if ( !centralityString.Contains("0-100%") )
            collisionSystem = Form("%s %s", centralityString.Data(), collisionSystem.Data());
    }

    TString labelALICEforPlots                      = "ALICE";
    TString outputDir                               = Form("%s/%s/MakePhotonQAPlots_%s%s",suffix.Data(),dateForOutput.Data(),cutnumber.Data(),addName.Data());
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
    TH3F* histoGammaAlphaQtPt                   = NULL;
    TH3F* histoGammaAlphaQtPtCopy               = NULL;
    TH3F* histoGammaChi2PsiPairPt               = NULL;
    TH3F* histoGammaChi2PsiPairPtCopy           = NULL;
    // projections
    TH2D* histoGammaAlphaQt                     = NULL;
    TH2D* histoGammaAlphaQtPtSliced[15]         = {NULL};
    TH2D* histoGammaAlphaPt                     = NULL;
    TH2D* histoGammaQtPt                        = NULL;
    TH2D* histoGammaPsiPairPt                   = NULL;
    TH2D* histoGammaChi2PsiPair                 = NULL;
    TH2D* histoGammaChi2PsiPairPtSliced[15]     = {NULL};
    TH2D* histoGammaChi2Pt                      = NULL;

    TH3F* histoTrueMCKindChi2PsiPairPt[20]      = {NULL};
    TH3F* histoTrueMCKindChi2PsiPairPtCopy[20]  = {NULL};
    TH3F* histoTrueMCKindAlphaQtPt[20]          = {NULL};
    TH3F* histoTrueMCKindAlphaQtPtCopy[20]      = {NULL};
    TH3F* histoTrueGammaChi2PsiPairPt           = NULL;

    TH2D* histoTrueMCKindQtPt[20]               = {NULL};
    TH2D* histoTrueMCKindAlphaPt[20]            = {NULL};
    TH2D* histoTrueMCKindPsiPairPt[20]          = {NULL};
    TH2D* histoTrueMCKindChi2Pt[20]             = {NULL};
    TH2D* histoTrueMCKindChi2PsiPair[20]        = {NULL};
    TH2D* histoTrueMCKindChi2PsiPairPtSliced[20][15]    = {NULL};
    TH2D* histoTrueMCKindAlphaQt[20]            = {NULL};
    TH2D* histoTrueMCKindAlphaQtPtSliced[20][15]= {NULL};
    TH2D* histoTrueGammaChi2PsiPair             = NULL;

    histoGammaAlphaQtPt                     = (TH3F*)directoryConv->Get("histoGammaAlphaQtPt");
    histoGammaChi2PsiPairPt                 = (TH3F*)directoryConv->Get("histoGammaChi2PsiPairPt");

    histoGammaAlphaQt                       = (TH2D*)histoGammaAlphaQtPt->Project3D("yx");
    histoGammaAlphaPt                       = (TH2D*)histoGammaAlphaQtPt->Project3D("zx");
    histoGammaQtPt                          = (TH2D*)histoGammaAlphaQtPt->Project3D("zy");
    histoGammaPsiPairPt                     = (TH2D*)histoGammaChi2PsiPairPt->Project3D("zy");
    histoGammaChi2PsiPair                   = (TH2D*)histoGammaChi2PsiPairPt->Project3D("yx");
    histoGammaChi2Pt                        = (TH2D*)histoGammaChi2PsiPairPt->Project3D("zx");
    // copies for projecting
    histoGammaChi2PsiPairPtCopy             = (TH3F*)histoGammaChi2PsiPairPt->Clone("histoGammaChi2PsiPairPtCopy");
    histoGammaAlphaQtPtCopy                 = (TH3F*)histoGammaAlphaQtPt->Clone("histoGammaAlphaQtPtCopy");
    for(Int_t j=0; j<15; j++){
        histoGammaChi2PsiPairPtCopy->GetZaxis()->SetRangeUser(projBinnin[j],projBinnin[j+1]);
        histoGammaChi2PsiPairPtSliced[j]   = (TH2D*)histoGammaChi2PsiPairPtCopy->Project3D(Form("yx%d",j));
        histoGammaAlphaQtPtCopy->GetZaxis()->SetRangeUser(projBinnin[j],projBinnin[j+1]);
        histoGammaAlphaQtPtSliced[j]       = (TH2D*)histoGammaAlphaQtPtCopy->Project3D(Form("yx%d",j));
    }
    if (isMC){
        for(Int_t i=0;i<20;i++){
            if(i!=0 && i!=10 && i!=11 && i!=13) continue;
            histoTrueMCKindChi2PsiPairPt[i] = (TH3F*)directoryConv->Get(Form("histoTrueMCKindChi2PsiPairPt_kind%d",i));
            histoTrueMCKindAlphaQtPt[i]     = (TH3F*)directoryConv->Get(Form("histoTrueMCKindAlphaQtPt_kind%d",i));

            histoTrueMCKindQtPt[i]          = (TH2D*)histoTrueMCKindAlphaQtPt[i]->Project3D(Form("zy%d",i));
            histoTrueMCKindAlphaPt[i]       = (TH2D*)histoTrueMCKindAlphaQtPt[i]->Project3D(Form("zx%d",i));
            histoTrueMCKindAlphaQt[i]       = (TH2D*)histoTrueMCKindAlphaQtPt[i]->Project3D(Form("yx%d",i));

            histoTrueMCKindPsiPairPt[i]     = (TH2D*)histoTrueMCKindChi2PsiPairPt[i]->Project3D(Form("zy%d",i));
            histoTrueMCKindChi2Pt[i]        = (TH2D*)histoTrueMCKindChi2PsiPairPt[i]->Project3D(Form("zx%d",i));
            histoTrueMCKindChi2PsiPair[i]   = (TH2D*)histoTrueMCKindChi2PsiPairPt[i]->Project3D(Form("yx%d",i));

            histoTrueMCKindChi2PsiPairPtCopy[i] = (TH3F*)histoTrueMCKindChi2PsiPairPt[i]->Clone(Form("histoTrueMCKindChi2PsiPairPtCopy%d",i));
            histoTrueMCKindAlphaQtPtCopy[i]     = (TH3F*)histoTrueMCKindAlphaQtPt[i]->Clone(Form("histoTrueMCKindAlphaQtPtCopy%d",i));
            for(Int_t j=0; j<15; j++){
                histoTrueMCKindChi2PsiPairPtCopy[i]->GetZaxis()->SetRangeUser(projBinnin[j],projBinnin[j+1]);
                histoTrueMCKindChi2PsiPairPtSliced[i][j]   = (TH2D*)histoTrueMCKindChi2PsiPairPtCopy[i]->Project3D(Form("yx%d",j));
                histoTrueMCKindAlphaQtPtCopy[i]->GetZaxis()->SetRangeUser(projBinnin[j],projBinnin[j+1]);
                histoTrueMCKindAlphaQtPtSliced[i][j]       = (TH2D*)histoTrueMCKindAlphaQtPtCopy[i]->Project3D(Form("yx%d",j));
            }
        }
        histoTrueGammaChi2PsiPairPt         = (TH3F*)directoryConv->Get("histoTrueGammaChi2PsiPairPt");
        histoTrueGammaChi2PsiPair           = (TH2D*)histoTrueGammaChi2PsiPairPt->Project3D("yx");
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
    // initialize cuts for QT vs pT
    Double_t funcParamQtPt[4][2]    = {{0.11,0.04},{0.125,0.05},{0.14,0.06},{0.16,0.07}};
    if (optionEnergy.CompareTo("XeXe_5.44TeV") == 0 || optionEnergy.CompareTo("13TeVLowB") == 0){
        funcParamQtPt[0][0]    = 0.18;
        funcParamQtPt[0][1]    = 0.032;
        funcParamQtPt[1][0]    = 0.2;
        funcParamQtPt[1][1]    = 0.035;
        funcParamQtPt[2][0]    = 0.25;
        funcParamQtPt[2][1]    = 0.04;
        funcParamQtPt[3][0]    = 0.3;
        funcParamQtPt[3][1]    = 0.045;
    }
    Color_t colorQtFunc[5]          = {kMagenta+2, kMagenta+4, kMagenta-4, kMagenta-8, kMagenta+1};
    Style_t styleQtFunc[5]          = {1,7,8,9,6};

    DrawGammaLines(0.05, 0.05, 0.04, 8, 3, kGreen-8, 9);
    TF1* funcPtDepQtCut[4]      = {NULL};
    for (Int_t k = 0; k< 4; k++){
        funcPtDepQtCut[k] = new TF1(Form("funcPtDepQtCut_%d",k),"x/[0]",0.005,funcParamQtPt[k][1]);
        DrawGammaSetMarkerTF1( funcPtDepQtCut[k], styleQtFunc[k+1], 2, colorQtFunc[k+1]);
        funcPtDepQtCut[k]->SetParameter(0,funcParamQtPt[k][0]);
        funcPtDepQtCut[k]->Draw("same");
        DrawGammaLines(funcParamQtPt[k][1], funcParamQtPt[k][1], funcParamQtPt[k][1]/funcParamQtPt[k][0], 8, 2, colorQtFunc[k+1], styleQtFunc[k+1]);
    }
    TF1 *funcOldCutDummy = new TF1("funcOldCutDummy","x/[0]",0.005,1.);
    DrawGammaSetMarkerTF1( funcOldCutDummy, 9, 3, kGreen-8);
    TLegend* legendQtPlotFits  = GetAndSetLegend2(0.58, 0.13, 0.95, 0.13+(0.035*5*1.15), 0.75*textSizeLabelsPixel,1,"",43,0.1);
    legendQtPlotFits->AddEntry(funcOldCutDummy, "#it{q}_{T}^{max} = 0.05","l");
    for (Int_t k = 0; k< 4; k++){
        legendQtPlotFits->AddEntry(funcPtDepQtCut[k], Form("#it{q}_{T}^{max} = %0.3f#it{p}_{T}, #it{q}_{T}^{max} = %0.3f",funcParamQtPt[k][0],funcParamQtPt[k][1] ),"l");
    }
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
    for (Int_t k = 0; k< 4; k++){
        funcPtDepQtCut[k]->Draw("same");
        DrawGammaLines(funcParamQtPt[k][1], funcParamQtPt[k][1], funcParamQtPt[k][1]/funcParamQtPt[k][0], 8, 2, colorQtFunc[k+1], styleQtFunc[k+1]);
    }
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

    Double_t funcParamChi2PsiPairLin[2][2]  = {{0.1,30}, {0.15, 50}};
    Double_t funcParamChi2PsiPair[3][2]     = {{0.15,-0.065}, {0.18, -0.055}, { 0.20, -0.050}};
    if (optionEnergy.CompareTo("XeXe_5.44TeV") == 0 || optionEnergy.CompareTo("13TeVLowB") == 0){
        funcParamChi2PsiPair[0][0]     = 0.3;
        funcParamChi2PsiPair[0][1]     = -0.085;
        funcParamChi2PsiPair[1][0]     = 0.35;
        funcParamChi2PsiPair[1][1]     = -0.075;
        funcParamChi2PsiPair[2][0]     = 0.4;
        funcParamChi2PsiPair[2][1]     = -0.065;
    }
    Color_t colorTriFunc[2]                 = {kMagenta+2, kMagenta+4};
    Color_t colorModFunc[3]                 = {kMagenta-4, kMagenta-8, kMagenta+1};

    TF1 *funcChi2PsiPairTri[4]      = {NULL};
    for (Int_t k = 0; k< 2; k++){
        funcChi2PsiPairTri[k*2] =  new TF1(Form("funcChi2PsiPairTri_%d",k*2),Form("-%0.2f+(%0.2f/%f*x)",funcParamChi2PsiPairLin[k][0],funcParamChi2PsiPairLin[k][0],funcParamChi2PsiPairLin[k][1]),0,funcParamChi2PsiPairLin[k][1]);
        DrawGammaSetMarkerTF1( funcChi2PsiPairTri[k*2], 2, 3, colorTriFunc[k]);
        funcChi2PsiPairTri[k*2]->Draw("same");
        funcChi2PsiPairTri[k*2+1] =   new TF1(Form("funcChi2PsiPairTri_%d",k*2+1),Form("%0.2f-(%0.2f/%f*x)",funcParamChi2PsiPairLin[k][0],funcParamChi2PsiPairLin[k][0],funcParamChi2PsiPairLin[k][1]),0,funcParamChi2PsiPairLin[k][1]);
        DrawGammaSetMarkerTF1( funcChi2PsiPairTri[k*2+1], 2, 3, colorTriFunc[k]);
        funcChi2PsiPairTri[k*2+1]->Draw("same");
    }
    TF1 *funcChi2PsiPairExp[6]      = {NULL};
    for (Int_t k = 0; k< 3; k++){
        funcChi2PsiPairExp[k*2] =  new TF1(Form("funcChi2PsiPairExp_%d",k*2),Form("%0.2f*TMath::Exp(%0.3f*x)",funcParamChi2PsiPair[k][0],funcParamChi2PsiPair[k][1]),0,50);
        DrawGammaSetMarkerTF1( funcChi2PsiPairExp[k*2], 2, 3, colorModFunc[k]);
        funcChi2PsiPairExp[k*2]->Draw("same");
        funcChi2PsiPairExp[k*2+1] =  new TF1(Form("funcChi2PsiPairExp_%d",k*2+1),Form("-%0.2f*TMath::Exp(%0.3f*x)",funcParamChi2PsiPair[k][0],funcParamChi2PsiPair[k][1]),0,50);
        DrawGammaSetMarkerTF1( funcChi2PsiPairExp[k*2+1], 2, 3, colorModFunc[k]);
        funcChi2PsiPairExp[k*2+1]->Draw("same");

    }
    TLegend* legendChi2PsiPairPlotFits  = GetAndSetLegend2(0.65, 0.13, 0.95, 0.13+(0.035*5*1.15), 0.75*textSizeLabelsPixel,1, "", 43, 0.1);
    for (Int_t k = 0; k< 2; k++){
        legendChi2PsiPairPlotFits->AddEntry(funcChi2PsiPairTri[k*2], Form("|#Psi_{pair}| < %0.2f/%.0f#chi^{2}+%0.2f",funcParamChi2PsiPairLin[k][0], funcParamChi2PsiPairLin[k][1],funcParamChi2PsiPairLin[k][0]),"l");
    }
    for (Int_t k = 0; k< 3; k++){
        legendChi2PsiPairPlotFits->AddEntry(funcChi2PsiPairExp[k*2], Form("|#Psi_{pair}| < %0.2fexp(%0.3f#chi^{2})",funcParamChi2PsiPair[k][0],funcParamChi2PsiPair[k][1]),"l");
    }
    legendChi2PsiPairPlotFits->Draw();
    labelEnergy->Draw();
    labelProcess->Draw();
    histo2DChi2PsiPairDummy->Draw("axis,same");
    canvasChi2PsiPairPlots->SaveAs(Form("%s/Chi2PsiPair_vs_Pt_Final_withCuts_%s.%s",outputDir.Data(),isMC ? "MC" : "data",suffix.Data()));

    // TRUE GAMMAS AND COMBINATORIAL GAMMA CANDIDATES PLOT
    histo2DChi2PsiPairDummy->Draw("copy");
    if (isMC){
        // draw true gamma->ee
        histoTrueGammaChi2PsiPair->Draw("col,same");
        for (Int_t k = 0; k< 2; k++){
            funcChi2PsiPairTri[k*2]->Draw("same");
            funcChi2PsiPairTri[k*2+1]->Draw("same");
        }
        for (Int_t k = 0; k< 3; k++){
            funcChi2PsiPairExp[k*2]->Draw("same");
            funcChi2PsiPairExp[k*2+1]->Draw("same");
        }
        legendChi2PsiPairPlotFits->Draw();
        labelEnergy->Draw();
        labelProcess->Draw();
        histo2DChi2PsiPairDummy->Draw("axis,same");
        canvasChi2PsiPairPlots->SaveAs(Form("%s/Chi2PsiPair_vs_Pt_trueGamma.%s",outputDir.Data(), suffix.Data()));
    }



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
        histoTrueMCKindAlphaQt[0]->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTrueMCKindAlphaQt[0]->Draw("col,same");
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
        histoTrueMCKindAlphaQt[0]->Draw("col,same");
        ex1 = new TExec("ex1","PalColor();");
        ex1->Draw();
        histoTrueMCKindAlphaQt[0]->Draw("col,same");
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
