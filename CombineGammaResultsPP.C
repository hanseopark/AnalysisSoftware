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
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"
#include "CombineGammaResultsPP.h"

extern TRandom*    gRandom;
extern TBenchmark* gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*    gMinuit;

struct SysErrorConversion {
    Double_t value;
    Double_t error;
    //    TString name;
};

void drawLatexAdd(TString latextext, Double_t textcolumn, Double_t textrow, Double_t textSizePixel,Bool_t setFont = kFALSE, Bool_t setFont2 = kFALSE, Bool_t alignRight = kFALSE, Color_t textcolor = kBlack){
    TLatex *latexDummy                  = new TLatex(textcolumn ,textrow,latextext);
    SetStyleTLatex( latexDummy, textSizePixel,4);
    if(setFont)
        latexDummy->SetTextFont(62);
    if(setFont2)
        latexDummy->SetTextFont(43);
    if(alignRight)
        latexDummy->SetTextAlign(31);
    latexDummy->SetTextColor(textcolor);
    latexDummy->Draw();
}
void drawhMarker(TH1D* histoDummy, Double_t column, Double_t row, Double_t markerScale){
    TMarker* markerdummy;
    markerdummy= CreateMarkerFromHisto(histoDummy,column ,row ,markerScale);
    markerdummy->DrawMarker(column ,row);
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
void plotLuminosity (Double_t ymin, Double_t ymax, TString ylabel,TString plotName){
    TH1D* histoData[5];
    TH1D* histoData2[7];
    TGraph *graphAliceLuminosity2010[5];
    TGraph *graphAliceLuminosity2012[7];
    graphAliceLuminosity2010[0] = new TGraph("ThesisQAInputAllEnergies/LumiInput/lumi10c900g","%lg %*lg %*lg %*lg %*lg %*lg %lg %*lg %*lg");
    graphAliceLuminosity2010[1] = new TGraph("ThesisQAInputAllEnergies/LumiInput/lumi10b","%lg %*lg %*lg %*lg %*lg %*lg %lg %*lg %*lg");
    graphAliceLuminosity2010[2] = new TGraph("ThesisQAInputAllEnergies/LumiInput/lumi10c","%lg %*lg %*lg %*lg %*lg %*lg %lg %*lg %*lg");
    graphAliceLuminosity2010[3] = new TGraph("ThesisQAInputAllEnergies/LumiInput/lumi10d","%lg %*lg %*lg %*lg %*lg %*lg %lg %*lg %*lg");
    graphAliceLuminosity2010[4] = new TGraph("ThesisQAInputAllEnergies/LumiInput/lumi10e","%lg %*lg %*lg %*lg %*lg %*lg %lg %*lg %*lg");
    graphAliceLuminosity2010[5] = new TGraph("ThesisQAInputAllEnergies/LumiInput/lumi10f","%lg %*lg %*lg %*lg %*lg %*lg %lg %*lg %*lg");
    graphAliceLuminosity2012[0] = new TGraph("ThesisQAInputAllEnergies/LumiInput/Lumi12a.txt","%lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %lg");
    graphAliceLuminosity2012[1] = new TGraph("ThesisQAInputAllEnergies/LumiInput/Lumi12b.txt","%lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %lg");
    graphAliceLuminosity2012[2] = new TGraph("ThesisQAInputAllEnergies/LumiInput/Lumi12c.txt","%lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %lg");
    graphAliceLuminosity2012[3] = new TGraph("ThesisQAInputAllEnergies/LumiInput/Lumi12d.txt","%lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %lg");
    graphAliceLuminosity2012[4] = new TGraph("ThesisQAInputAllEnergies/LumiInput/Lumi12f.txt","%lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %lg");
    graphAliceLuminosity2012[5] = new TGraph("ThesisQAInputAllEnergies/LumiInput/Lumi12h.txt","%lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %lg");
    graphAliceLuminosity2012[6] = new TGraph("ThesisQAInputAllEnergies/LumiInput/Lumi12i.txt","%lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %lg");

    // for( Int_t i = 0; i<7; i++){
    //     histoData2[i]      = new TH1D(graphAliceLuminosity2012[i]);
    //     if(i<5)
    //     histoData[i]      = new TH1D(graphAliceLuminosity2010[i]);
    // }

    Double_t textSizeLabelsPixel                = 48;
    Double_t arrayBoundariesX1_4[3];
    Double_t arrayBoundariesY1_4[2];
    Double_t relativeMarginsX[3];
    Double_t relativeMarginsY[3];
    TGaxis::SetMaxDigits(6);
    TCanvas* canvasDummy     = new TCanvas("canvasDummy","",0,0,1650,900);  // gives the page size
    DrawGammaCanvasSettings( canvasDummy,  0., 0., 0., 0.0);

    TPad* pad900and7               = new TPad("pad900and7", "",0, 0, 0.5, 1,-1, -1, -2);
    DrawGammaPadSettings( pad900and7, 0.12, 0.,           0.02, 0.07);
    //                                       ^rechts     ^oben
    pad900and7->Draw();

    TPad* pad8                = new TPad("pad8", "", 0.5, 0, 1,1,-1, -1, -2);
    DrawGammaPadSettings( pad8, 0, 0.01, 0.02, 0.07);
    pad8->Draw();


    Double_t textsizeLabelsPP = 0;
    Double_t textsizeFacPP= 0;
    if (pad900and7->XtoPixel(pad900and7->GetX2()) <pad900and7->YtoPixel(pad900and7->GetY1()) ){
        textsizeLabelsPP = (Double_t)textSizeLabelsPixel/pad900and7->XtoPixel(pad900and7->GetX2()) ;
        textsizeFacPP = (Double_t)1./pad900and7->XtoPixel(pad900and7->GetX2()) ;
    } else {
        textsizeLabelsPP = (Double_t)textSizeLabelsPixel/pad900and7->YtoPixel(pad900and7->GetY1());
        textsizeFacPP = (Double_t)1./pad900and7->YtoPixel(pad900and7->GetY1());
    }

    pad900and7->cd();
    pad900and7->SetLogy();
    TH2F * histoPlottingDummy;
    histoPlottingDummy                       = new TH2F("histoPlottingDummy","histoPlottingDummy",1000,950,1449,1000,ymin,ymax   );
    SetStyleHistoTH2ForGraphs(histoPlottingDummy, "",ylabel.Data(), 0.55*textsizeLabelsPP, 0.85*textsizeLabelsPP,
                              0.65*textsizeLabelsPP,0.8*textsizeLabelsPP, 0.4, 1.15, 80511, 508);
    histoPlottingDummy->Draw("copy");
    // Data 900 GeV
    DrawGammaSetMarkerTGraph(graphAliceLuminosity2010[0], 34 ,2, kRed+2, kRed+2);
    graphAliceLuminosity2010[0]->Draw("p,same,e");

    // Data 7 TeV
    for (Int_t i = 1; i < 6; i++){
        DrawGammaSetMarkerTGraph(graphAliceLuminosity2010[i], 20 ,1.8, kBlue+2, kBlue+2);
        graphAliceLuminosity2010[i]->Draw("p,same,e");
    }
    TLegend* legend7900           = GetAndSetLegend2(0.15, 0.81, 0.43, 0.81+2*textsizeLabelsPP,42);
                legend7900->AddEntry(graphAliceLuminosity2010[0],nameMeasGlobal[0].Data(),"p");
                legend7900->AddEntry(graphAliceLuminosity2010[1],nameMeasGlobal[1].Data(),"p");
        legend7900->Draw();
    pad8->cd();
    pad8->SetLogy();
    TH2F * histoPlottingDummy2;
    histoPlottingDummy2                       = new TH2F("histoPlottingDummy2","histoPlottingDummy2",1000,2401,3529,1000,ymin,ymax    );
    SetStyleHistoTH2ForGraphs(histoPlottingDummy2, "fill number","", 0.55*textsizeLabelsPP, 0.75*textsizeLabelsPP,
                              0*textsizeLabelsPP,0.85*textsizeLabelsPP, 0.7, 0.95, 80509, 508);
    histoPlottingDummy2->Draw("copy");
    // Data 8 TeV
    for (Int_t i = 0; i < 7; i++){
        DrawGammaSetMarkerTGraph(graphAliceLuminosity2012[i], 29 ,2, kGreen+2, kGreen+2);
        graphAliceLuminosity2012[i]->Draw("p,same,e");
    }
    TLegend* legend8           = GetAndSetLegend2(0.55, 0.13, 0.93, 0.13+textsizeLabelsPP,42);
                legend8->AddEntry(graphAliceLuminosity2012[2],nameMeasGlobal[2].Data(),"p");
        legend8->Draw();

    canvasDummy->SaveAs(Form("%s/%s.%s",outputDir.Data(),plotName.Data(),suffix.Data()));
}
void plotHistograms (TH1D* histoData[], TH1D* histoMC[], Double_t ymin, Double_t ymax, TString ylabel,TString plotName){
    Double_t textSizeLabelsPixel                = 48;
    Double_t arrayBoundariesX1_4[3];
    Double_t arrayBoundariesY1_4[2];
    Double_t relativeMarginsX[3];
    Double_t relativeMarginsY[3];
    TGaxis::SetMaxDigits(6);
    TCanvas* canvasDummy     = new TCanvas("canvasDummy","",0,0,1650,900);  // gives the page size
    DrawGammaCanvasSettings( canvasDummy,  0., 0., 0., 0.0);

    TPad* pad900and7               = new TPad("pad900and7", "",0, 0, 0.5, 1,-1, -1, -2);
    DrawGammaPadSettings( pad900and7, 0.14, 0.,           0.09, 0.07);
    //                                       ^rechts     ^oben
    pad900and7->Draw();

    TPad* pad8                = new TPad("pad8", "", 0.5, 0, 1,1,-1, -1, -2);
    DrawGammaPadSettings( pad8, 0, 0.01, 0.09, 0.07);
    pad8->Draw();

    TPad* padMassLegend1            = new TPad("padMassLegend1", "", 0.07, 0.88, 0.98, 1.0,-1, -1, -2);
    DrawGammaPadSettings( padMassLegend1, 0., 0., 0., 0.);
    padMassLegend1->SetFillStyle(0);
    padMassLegend1->Draw();

    Double_t textsizeLabelsPP = 0;
    Double_t textsizeFacPP= 0;
    if (pad900and7->XtoPixel(pad900and7->GetX2()) <pad900and7->YtoPixel(pad900and7->GetY1()) ){
        textsizeLabelsPP = (Double_t)textSizeLabelsPixel/pad900and7->XtoPixel(pad900and7->GetX2()) ;
        textsizeFacPP = (Double_t)1./pad900and7->XtoPixel(pad900and7->GetX2()) ;
    } else {
        textsizeLabelsPP = (Double_t)textSizeLabelsPixel/pad900and7->YtoPixel(pad900and7->GetY1());
        textsizeFacPP = (Double_t)1./pad900and7->YtoPixel(pad900and7->GetY1());
    }

    pad900and7->cd();
    TH2F * histoPlottingDummy;
    histoPlottingDummy                       = new TH2F("histoPlottingDummy","histoPlottingDummy",1000,113001,136899,1000,ymin,ymax   );
    SetStyleHistoTH2ForGraphs(histoPlottingDummy, "",ylabel.Data(), 0.55*textsizeLabelsPP, 0.85*textsizeLabelsPP,
                              0.65*textsizeLabelsPP,0.8*textsizeLabelsPP, 0.4, 1.45, 80511, 508);
    histoPlottingDummy->Draw("copy");
    // MC 900GeV
    DrawGammaSetMarker(histoMC[0], 28 ,2, kRed-8, kRed-8);
    histoMC[0]->Draw("p,same,e");
    // Data 900 GeV
    DrawGammaSetMarker(histoData[0], 34 ,2, kRed+2, kRed+2);
    histoData[0]->Draw("p,same,e");

    // MC 7 TeV
    DrawGammaSetMarker(histoMC[1], 24 ,1.8, kBlue-6, kBlue-6);
    histoMC[1]->Draw("p,same,e");
    // Data 7 TeV
    for (Int_t i = 1; i < 6; i++){
        DrawGammaSetMarker(histoData[i], 20 ,1.8, kBlue+2, kBlue+2);
        histoData[i]->Draw("p,same,e");
    }

    pad8->cd();
    TH2F * histoPlottingDummy2;
    histoPlottingDummy2                       = new TH2F("histoPlottingDummy2","histoPlottingDummy2",1000,175001,194999,1000,ymin,ymax    );
    SetStyleHistoTH2ForGraphs(histoPlottingDummy2, "run number","", 0.55*textsizeLabelsPP, 0.75*textsizeLabelsPP,
                              0*textsizeLabelsPP,0.85*textsizeLabelsPP, 0.55, 0.95, 80509, 508);
    histoPlottingDummy2->Draw("copy");
    // MC Pythia 8 TeV
    DrawGammaSetMarker(histoMC[2], 30 ,2, kGreen-8, kGreen-8);
    histoMC[2]->Draw("p,same,e");
    // MC Phojet 8 TeV
    DrawGammaSetMarker(histoMC[3], 30 ,2, kGreen-9, kGreen-9);
    histoMC[3]->Draw("p,same,e");
    // Data 8 TeV
    for (Int_t i = 6; i < nSets; i++){
        DrawGammaSetMarker(histoData[i], 29 ,2, kGreen+2, kGreen+2);
        histoData[i]->Draw("p,same,e");
    }
    // Double_t textsizeLegendPixel        = 30;
    Double_t textsizeLegendPixel        = 0.32;
    Double_t columnsLegendMass2[7]      = {0.,.11,.25,.35,.6,.7,.8};
        // Double_t rowsLegendMass2[5] = {0.8,0.6,0.4,0.2,0.01};
        // Double_t rowsLegendMass2[6] = {0.84,0.66,0.50,0.33,0.16,0.01};
          Double_t  rowsLegendMass2[2]= {0.7,0.4};
        //******************* Offsets ***********************
        Double_t offsetMarkerXMass2         = 0.012;
        Double_t offsetMarkerYMass2         = 0.11;
        //****************** Scale factors ******************
        Double_t scaleMarkerMass2           = 1.2;

        padMassLegend1->cd();
        //****************** first Column **************************************************

        drawLatexAdd("0.9 TeV:",columnsLegendMass2[0],rowsLegendMass2[0],textsizeLegendPixel,kTRUE);
        drawLatexAdd("Data",columnsLegendMass2[1],rowsLegendMass2[0],textsizeLegendPixel,kTRUE);
        drawLatexAdd("Pythia 6",columnsLegendMass2[1],rowsLegendMass2[1],textsizeLegendPixel,kTRUE);

        drawLatexAdd("7 TeV:",columnsLegendMass2[2],rowsLegendMass2[0],textsizeLegendPixel,kTRUE);
        drawLatexAdd("Data",columnsLegendMass2[3],rowsLegendMass2[0],textsizeLegendPixel,kTRUE);
        drawLatexAdd("Pythia 6",columnsLegendMass2[3],rowsLegendMass2[1],textsizeLegendPixel,kTRUE);

        drawLatexAdd("8 TeV:",columnsLegendMass2[4],rowsLegendMass2[0],textsizeLegendPixel,kTRUE);
        drawLatexAdd("Data",columnsLegendMass2[5],rowsLegendMass2[0],textsizeLegendPixel,kTRUE);
        drawLatexAdd("Pythia 8",columnsLegendMass2[6],rowsLegendMass2[0],textsizeLegendPixel,kTRUE);
        drawLatexAdd("Phojet",columnsLegendMass2[6],rowsLegendMass2[1],textsizeLegendPixel,kTRUE);

        drawhMarker(histoData[0],columnsLegendMass2[1]-offsetMarkerXMass2,rowsLegendMass2[0]+ offsetMarkerYMass2,scaleMarkerMass2);
        drawhMarker(histoMC[0],columnsLegendMass2[1]-offsetMarkerXMass2,rowsLegendMass2[1]+ offsetMarkerYMass2,scaleMarkerMass2);
        drawhMarker(histoData[1],columnsLegendMass2[3]-offsetMarkerXMass2,rowsLegendMass2[0]+ offsetMarkerYMass2,scaleMarkerMass2);
        drawhMarker(histoMC[1],columnsLegendMass2[3]-offsetMarkerXMass2,rowsLegendMass2[1]+ offsetMarkerYMass2,scaleMarkerMass2);
        drawhMarker(histoData[6],columnsLegendMass2[5]-offsetMarkerXMass2,rowsLegendMass2[0]+ offsetMarkerYMass2,scaleMarkerMass2);
        drawhMarker(histoMC[2],columnsLegendMass2[6]-offsetMarkerXMass2,rowsLegendMass2[0]+ offsetMarkerYMass2,scaleMarkerMass2);
        drawhMarker(histoMC[3],columnsLegendMass2[6]-offsetMarkerXMass2,rowsLegendMass2[1]+ offsetMarkerYMass2,scaleMarkerMass2);

    canvasDummy->SaveAs(Form("%s/%s.%s",outputDir.Data(),plotName.Data(),suffix.Data()));
}

void plotYield(TH1D* yieldHisto[], Double_t yMin, Double_t yMax, TString meson, TString plotName)
    {


    Double_t textSizeLabelsPixel             = 55;
    Double_t textSizeLabelsRel      = 55./1200;

    TCanvas* canvasDummy       = new TCanvas("canvasDummy", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasDummy,  0.1, 0.01, 0.015, 0.095);
    canvasDummy->SetLogy(1);
    canvasDummy->SetLogx(1);

        TH2F * histoDummy;
            histoDummy                = new TH2F("histoDummy", "histoDummy",1000, 0.23,  19.9, 1000, yMin, yMax );
        SetStyleHistoTH2ForGraphs( histoDummy, "#it{p}_{T} (GeV/#it{c})", "Raw Yield per Event",
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.08);//(#times #epsilon_{pur})
        histoDummy->GetYaxis()->SetLabelOffset(0.001);
        histoDummy->GetXaxis()->SetLabelOffset(-0.01);
        histoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
        histoDummy->DrawCopy();

        if(!plotName.CompareTo("GammaRawYield")){
            yieldHisto[0] ->Scale(1);
            yieldHisto[3] ->Scale(10);
            yieldHisto[1] ->Scale(100);
            yieldHisto[2] ->Scale(1000);
        }else{
            yieldHisto[0] ->Scale(1);
            yieldHisto[1] ->Scale(10);
            yieldHisto[2] ->Scale(100);
        }

        DrawGammaSetMarker(yieldHisto[0], markerstyles[0], markersize[0], colorData[0] , colorData[0]);
        yieldHisto[0]->Draw("p,same,e");
        if(!plotName.CompareTo("GammaRawYield")){
            DrawGammaSetMarker(yieldHisto[3], markerstyles[3], markersize[3], colorData[3] , colorData[3]);
            yieldHisto[3]->Draw("p,same,e");
        }
        DrawGammaSetMarker(yieldHisto[1], markerstyles[1], markersize[1], colorData[1] , colorData[1]);
        yieldHisto[1]->Draw("p,same,e");

        DrawGammaSetMarker(yieldHisto[2], markerstyles[2], markersize[2], colorData[2] , colorData[2]);
        yieldHisto[2]->Draw("p,same,e");

        TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.15, 0.13, 0.43, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel);
        if(!plotName.CompareTo("GammaRawYield")){
            legendEffiAccPi0           = GetAndSetLegend2(0.15, 0.13, 0.43, 0.13+(4.5*textSizeLabelsRel),textSizeLabelsPixel);
            legendEffiAccPi0->AddEntry(yieldHisto[2],Form("%s x10^{3}",nameMeasGlobal[2].Data()),"p");
            legendEffiAccPi0->AddEntry(yieldHisto[1],Form("%s x10^{2}",nameMeasGlobal[1].Data()),"p");
            legendEffiAccPi0->AddEntry(yieldHisto[3],Form("%s x10",nameMeasGlobal[3].Data()),"p");
            legendEffiAccPi0->AddEntry(yieldHisto[0],nameMeasGlobal[0].Data(),"p");
        }else{
            legendEffiAccPi0->AddEntry(yieldHisto[2],Form("%s x10^{2}",nameMeasGlobal[2].Data()),"p");
            legendEffiAccPi0->AddEntry(yieldHisto[1],Form("%s x10",nameMeasGlobal[1].Data()),"p");
            legendEffiAccPi0->AddEntry(yieldHisto[0],nameMeasGlobal[0].Data(),"p");
        }


        legendEffiAccPi0->Draw();

        drawLatexAdd("ALICE",0.9,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd("pp, PCM",0.9,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if(!meson.CompareTo("Pi0"))
            drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.9,0.82,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        else if(!meson.CompareTo("Eta"))
            drawLatexAdd("#eta #rightarrow #gamma#gamma",0.9,0.82,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        else
            drawLatexAdd("#gamma #rightarrow e^{+}e^{-}",0.9,0.82,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

    canvasDummy->Update();
    canvasDummy->Print(Form("%s/%s.%s",outputDir.Data(),plotName.Data(),suffix.Data()));

    }

void plotEverything(TH1D* histoInput[], Double_t yMin, Double_t yMax, Int_t logY, TString yLabel, TString plotName)
    {


    Double_t textSizeLabelsPixel             = 55;
    Double_t textSizeLabelsRel      = 55./1200;

    TCanvas* canvasDummy       = new TCanvas("canvasDummy", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasDummy,  0.1, 0.01, 0.015, 0.095);
    if(logY)
        canvasDummy->SetLogy(1);
    canvasDummy->SetLogx(1);

        TH2F * histoDummy;
            histoDummy                = new TH2F("histoDummy", "histoDummy",1000, 0.23,  19.9, 1000, yMin, yMax );
        SetStyleHistoTH2ForGraphs( histoDummy, "#it{p}_{T} (GeV/#it{c})", yLabel.Data(),
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.06);//(#times #epsilon_{pur})
        histoDummy->GetYaxis()->SetLabelOffset(0.001);
        histoDummy->GetXaxis()->SetLabelOffset(-0.01);
        histoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
        histoDummy->GetYaxis()->SetMoreLogLabels(kTRUE);
        histoDummy->GetYaxis()->SetNoExponent(kTRUE);
        histoDummy->GetYaxis()->SetNdivisions(805);
        if(!plotName.CompareTo("PileupCorrectionFactorPi0"))
            histoDummy->GetYaxis()->SetNdivisions(510);
        histoDummy->DrawCopy();


        DrawGammaSetMarker(histoInput[0], markerstyles[0], markersize[0], colorData[0] , colorData[0]);
        histoInput[0]->Draw("p,same,e");
        if(!plotName.CompareTo("PileupContamination")||!plotName.CompareTo("PhotonPurity")||!plotName.CompareTo("PhotonEfficiency")||!plotName.CompareTo("PileupCorrectionFactorPi0")||!plotName.CompareTo("GammaPileUpCorrectionFactor")){
            DrawGammaSetMarker(histoInput[3], markerstyles[3], markersize[3], colorData[3] , colorData[3]);
            histoInput[3]->Draw("p,same,e");
        }
        DrawGammaSetMarker(histoInput[1], markerstyles[1], markersize[1], colorData[1] , colorData[1]);
        histoInput[1]->Draw("p,same,e");

        DrawGammaSetMarker(histoInput[2], markerstyles[2], markersize[2], colorData[2] , colorData[2]);
        histoInput[2]->Draw("p,same,e");

        TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.55, 0.13, 0.93, 0.13+(3*textSizeLabelsRel),textSizeLabelsPixel);
        if(!plotName.CompareTo("BinShift"))
            legendEffiAccPi0           = GetAndSetLegend2(0.55, 0.81, 0.93, 0.81+(3*textSizeLabelsRel),textSizeLabelsPixel);
        if(!plotName.CompareTo("PileupContamination"))
            legendEffiAccPi0           = GetAndSetLegend2(0.55, 0.81-textSizeLabelsRel, 0.93, 0.81+(3*textSizeLabelsRel),textSizeLabelsPixel);
        if(!plotName.CompareTo("PhotonPurity")||!plotName.CompareTo("PhotonEfficiency")||!plotName.CompareTo("PileupCorrectionFactorPi0")||!plotName.CompareTo("GammaPileUpCorrectionFactor"))
            legendEffiAccPi0           = GetAndSetLegend2(0.55, 0.13, 0.93, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel);
                legendEffiAccPi0->AddEntry(histoInput[2],nameMeasGlobal[2].Data(),"p");
                legendEffiAccPi0->AddEntry(histoInput[1],nameMeasGlobal[1].Data(),"p");
            if(!plotName.CompareTo("PileupContamination")||!plotName.CompareTo("PhotonPurity")||!plotName.CompareTo("PhotonEfficiency")||!plotName.CompareTo("PileupCorrectionFactorPi0")||!plotName.CompareTo("GammaPileUpCorrectionFactor"))
                legendEffiAccPi0->AddEntry(histoInput[3],nameMeasGlobal[3].Data(),"p");
                legendEffiAccPi0->AddEntry(histoInput[0],nameMeasGlobal[0].Data(),"p");


        legendEffiAccPi0->Draw();

        if(!plotName.CompareTo("PhotonPurity")||!plotName.CompareTo("BinShift"))
            DrawGammaLines(0.23, 19.9,1.0, 1.0, 1, kGray+2, 2);

        drawLatexAdd("ALICE",0.15,0.92,textSizeLabelsRel,kFALSE,kFALSE);
        // drawLatexAdd("ALICE",0.15,0.92,textSizeLabelsRel,kFALSE,kFALSE);
        drawLatexAdd("pp, PCM",0.15,0.87,textSizeLabelsRel,kFALSE,kFALSE);
        if(!plotName.CompareTo("PileupContamination")||!plotName.CompareTo("PileupCorrectionFactorPi0"))
            drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,0.82,textSizeLabelsRel);

    canvasDummy->Update();
    canvasDummy->Print(Form("%s/%s.%s",outputDir.Data(),plotName.Data(),suffix.Data()));

    }

void plotPeriodwisePileup(TH1D* histoInput[], Double_t yMin, Double_t yMax, Int_t logY, TString yLabel, TString plotName)
    {
    TString periods[7]                 = {"a","b","c","d","f","h","i"};
    TFile* inputFilePWDCA[7];
    TH1D* DCAPileupEstimate[10];
    for(Int_t i=0;i<7;i++)
    {
        inputFilePWDCA[i]                      = new TFile(Form("ThesisQAInputAllEnergies/PileupStudies/LHC12%s_00010113_00200009227300008250404000_0152103500000000/8TeV/Pi0_data_GammaConvV1Correction_00010113_00200009227300008250404000_0152103500000000.root",periods[i].Data()));
        DCAPileupEstimate[i]              = (TH1D*)inputFilePWDCA[i]->Get("BGEstimateFromPileup");
    }
    // TFile* inputFilePWDCAGamma[7];
    // TH1D* DCAPileupEstimateGamma[10];
    // for(Int_t i=0;i<7;i++)
    // {
    //     inputFilePWDCAGamma[i]                      = new TFile(Form("ThesisQAInputAllEnergies/PileupStudies/LHC12%s_00010113_00200009227300008250404000_0152103500000000/8TeV/Gamma_Pi0_data_GammaConvV1Correction_00010113_00200009227300008250404000_0152103500000000.root",periods[i].Data()));
    //     DCAPileupEstimateGamma[i]              = (TH1D*)inputFilePWDCAGamma[i]->Get("PileUpCorrectionFactor");
    // }
    Double_t textSizeLabelsPixel             = 55;
    Double_t textSizeLabelsRel      = 55./1200;

    TCanvas* canvasDummy       = new TCanvas("canvasDummy", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasDummy,  0.1, 0.01, 0.015, 0.095);
    if(logY)
        canvasDummy->SetLogy(1);
    canvasDummy->SetLogx(1);

        TH2F * histoDummy;
            histoDummy                = new TH2F("histoDummy", "histoDummy",1000, 0.23,  15.9, 1000, yMin, yMax );
        SetStyleHistoTH2ForGraphs( histoDummy, "#it{p}_{T} (GeV/#it{c})", yLabel.Data(),
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.06);//(#times #epsilon_{pur})
        histoDummy->GetYaxis()->SetLabelOffset(0.001);
        histoDummy->GetXaxis()->SetLabelOffset(-0.01);
        histoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
        histoDummy->GetYaxis()->SetMoreLogLabels(kTRUE);
        histoDummy->GetYaxis()->SetNoExponent(kTRUE);
            histoDummy->GetYaxis()->SetNdivisions(510);
        histoDummy->DrawCopy();

        Color_t colorDCA[8]                        ={kRed+2,kOrange+2,kYellow+2,kGreen+2,kCyan+2,kBlue+2,kMagenta+2,kGray+2};

    for(Int_t i=0;i<7;i++)
    {
        DrawGammaSetMarker(DCAPileupEstimate[i], 24, markersize[1], colorDCA[i] , colorDCA[i]);
        DCAPileupEstimate[i]->Draw("p,same,e");
    }

        DrawGammaSetMarker(histoInput[2], markerstyles[2], markersize[2], colorDCA[7] , colorDCA[7]);
        histoInput[2]->Draw("p,same,e");

        TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.55, 0.13, 0.93, 0.13+(3*textSizeLabelsRel),textSizeLabelsPixel);
        legendEffiAccPi0           = GetAndSetLegend2(0.65, 0.13, 0.95, 0.13+(8*textSizeLabelsRel),textSizeLabelsPixel);
        legendEffiAccPi0->AddEntry(histoInput[2],"LHC12[a-i]","p");
        for(Int_t i=0;i<7;i++)
        {
            legendEffiAccPi0->AddEntry(DCAPileupEstimate[i],Form("LHC12%s",periods[i].Data()),"p");
        }

        legendEffiAccPi0->Draw();


        drawLatexAdd("ALICE",0.15,0.92,textSizeLabelsRel,kFALSE,kFALSE);
        drawLatexAdd("pp, PCM",0.15,0.87,textSizeLabelsRel,kFALSE,kFALSE);
            drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma #rightarrow 4e",0.15,0.82,textSizeLabelsRel);

    canvasDummy->Update();
    canvasDummy->Print(Form("%s/%s.%s",outputDir.Data(),plotName.Data(),suffix.Data()));

    // histoDummy->DrawCopy();
    // for(Int_t i=0;i<7;i++)
    // {
    //     DrawGammaSetMarker(DCAPileupEstimateGamma[i], markerstyles[1], markersize[1], colorDCA[i] , colorDCA[i]);
    //     DCAPileupEstimateGamma[i]->Draw("p,same,e");
    // }
    // DrawGammaSetMarker(PileUpCorrectionFactor[2], markerstyles[2], markersize[2], colorDCA[7] , colorDCA[7]);
    // PileUpCorrectionFactor[2]->Draw("p,same,e");
    // canvasDummy->Update();
    // canvasDummy->Print(Form("%s/%s.%s",outputDir.Data(),"GammaPileupPeriodwise",suffix.Data()));

    }

void pi0andEtacombiedPlots(TH1D* histoPi0[],TH1D* histoEta[], TString yLabel, Double_t yMin, Double_t yMax, TString plotname)
    {
    Double_t textSizeLabelsPixel             = 55;
    Double_t textSizeLabelsRel      = 55./1200;

    TCanvas* canvasDummy       = new TCanvas("canvasDummy", "", 200, 10, 1200, 1100);  // gives the page size
    // DrawGammaCanvasSettings( canvasDummy,  0.1, 0.01, 0.015, 0.095);

    TPad* padMassPi0                = new TPad("padMassPi0", "", 0,0, 1, 1,-1, -1, -2);
    DrawGammaPadSettings( padMassPi0, 0.1, 0.01, 0.015, 0.095);
    padMassPi0->Draw();

    TPad* padMassLegend1            = new TPad("padMassLegend1", "", 0.53, 0.12, 0.95, 0.35,-1, -1, -2);
    DrawGammaPadSettings( padMassLegend1, 0., 0., 0., 0.);
    padMassLegend1->SetFillStyle(0);
    padMassLegend1->Draw();

    padMassPi0->cd();
    if(plotname.CompareTo("AcceptanceComb") )
        padMassPi0->SetLogy(1);
    padMassPi0->SetLogx(1);
        TH2F * histoDummy;
            histoDummy                = new TH2F("histoDummy", "histoDummy",1000, 0.23,  19.9, 1000, yMin, yMax );
        SetStyleHistoTH2ForGraphs( histoDummy, "#it{p}_{T} (GeV/#it{c})", yLabel.Data(),
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.08);//(#times #epsilon_{pur})
        histoDummy->GetYaxis()->SetLabelOffset(0.001);
        histoDummy->GetXaxis()->SetLabelOffset(-0.01);
        histoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
        if(!plotname.CompareTo("AcceptanceComb"))
            histoDummy->GetYaxis()->SetMoreLogLabels(kTRUE);
        if(!plotname.CompareTo("AcceptanceComb"))
            histoDummy->GetYaxis()->SetNoExponent(kTRUE);
        histoDummy->DrawCopy();

        // histoPi0[0] ->Scale(1);
        // histoPi0[1] ->Scale(10);
        // histoPi0[2] ->Scale(100);


        DrawGammaSetMarker(histoPi0[1], markerstyles[1], markersize[1], colorData[1] , colorData[1]);
        histoPi0[1]->Draw("p,same,e");
        DrawGammaSetMarker(histoPi0[2], markerstyles[2], markersize[2], colorData[2] , colorData[2]);
        histoPi0[2]->Draw("p,same,e");
        DrawGammaSetMarker(histoPi0[0], markerstyles[0], markersize[0], colorData[0] , colorData[0]);
        histoPi0[0]->Draw("p,same,e");

        Int_t addColor = 9;
        DrawGammaSetMarker(histoEta[1], markerstylesMC[1], markersize[1], colorMC[1]+addColor , colorMC[1]+addColor);
        histoEta[1]->Draw("p,same,e");
        DrawGammaSetMarker(histoEta[2], markerstylesMC[2], markersize[2], colorMC[2]+addColor , colorMC[2]+addColor);
        histoEta[2]->Draw("p,same,e");
        DrawGammaSetMarker(histoEta[0], markerstylesMC[0], markersize[0], colorMC[0]+addColor , colorMC[0]+addColor);
        histoEta[0]->Draw("p,same,e");

        drawLatexAdd("ALICE simulation",0.15,0.92,textSizeLabelsRel,kFALSE,kFALSE);
        drawLatexAdd("pp, PCM",0.15,0.87,textSizeLabelsRel,kFALSE,kFALSE);
            drawLatexAdd("#pi^{0}, #eta #rightarrow #gamma#gamma",0.15,0.82,textSizeLabelsRel,kFALSE,kFALSE);

       //********************************** Defintion of the Legend **************************************************
        Double_t columnsLegendMass2[3]      = {0.,0.62,0.89};
        // Double_t rowsLegendMass2[5] = {0.8,0.6,0.4,0.2,0.01};
        // Double_t rowsLegendMass2[6] = {0.84,0.66,0.50,0.33,0.16,0.01};
          Double_t  rowsLegendMass2[7]= {0.80,0.58,0.35,0.13,0.01,0.16};
        //******************* Offsets ***********************
        Double_t offsetMarkerXMass2         = 0.1;
        Double_t offsetMarkerYMass2         = 0.08;
        //****************** Scale factors ******************
        Double_t scaleMarkerMass2           = 1.2;

        padMassLegend1->cd();
        //****************** first Column **************************************************
        TLatex *textMassPCM[10];
        for (Int_t i = 0; i < 3; i++){
                textMassPCM[i]                  = new TLatex(columnsLegendMass2[0],rowsLegendMass2[i+1],nameMeasGlobal[i].Data());
                SetStyleTLatex( textMassPCM[i], textSizeLabelsPixel,4);
                textMassPCM[i]->SetTextFont(43);
                textMassPCM[i]->Draw();
        }

        //****************** second Column *************************************************
        TLatex *textMassData                = new TLatex(columnsLegendMass2[1]+0.065,rowsLegendMass2[0] ,"#pi^{0}");
        SetStyleTLatex( textMassData, textSizeLabelsPixel,4);
        textMassData->SetTextFont(43);
        textMassData->Draw();
        TLatex *textMassMC                  = new TLatex(columnsLegendMass2[2]+0.035 ,rowsLegendMass2[0],"#eta");
        SetStyleTLatex( textMassMC, textSizeLabelsPixel,4);
        textMassMC->SetTextFont(43);
        textMassMC->Draw();

        TMarker* markerPCMPi0Mass[10];
        TMarker* markerPCMPi0MassMC[10];
        for (Int_t i = 0; i < 3; i++){
                markerPCMPi0Mass[i]             = CreateMarkerFromHisto(histoPi0[i],columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
                markerPCMPi0Mass[i]->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2);
                markerPCMPi0MassMC[i]           = CreateMarkerFromHisto(histoEta[i],columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
                markerPCMPi0MassMC[i]->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2);
        }


    canvasDummy->Update();
    canvasDummy->Print(Form("%s/%s.%s",outputDir.Data(),plotname.Data(),suffix.Data()));

    }

void plotSecondaries(Int_t index, TH1D* histo1[],TH1D* histo2[],TH1D* histo3[],TH1D* histo4[],TH1D* histo5[], Double_t yMin, Double_t yMax, TString yLabel, Int_t logY, TString plotName)
    {


    Double_t textSizeLabelsPixel             = 55;
    Double_t textSizeLabelsRel      = 55./1200;

    TCanvas* canvasDummy       = new TCanvas("canvasDummy", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasDummy,  0.1, 0.01, 0.015, 0.095);
    if(logY)
        canvasDummy->SetLogy(1);
    canvasDummy->SetLogx(1);

        TH2F * histoDummy;
            histoDummy                = new TH2F("histoDummy", "histoDummy",1000, 0.23,  19.9, 1000, yMin, yMax );
        SetStyleHistoTH2ForGraphs( histoDummy, "#it{p}_{T} (GeV/#it{c})", yLabel.Data(),
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.08);//(#times #epsilon_{pur})
        histoDummy->GetYaxis()->SetLabelOffset(0.001);
        histoDummy->GetXaxis()->SetLabelOffset(-0.01);
        histoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
        histoDummy->DrawCopy();
        TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.5, 0.13, 0.88, 0.13+(5*textSizeLabelsRel),textSizeLabelsPixel);
        if(!plotName.CompareTo("SecondaryYield")||!plotName.CompareTo("FracSecYield"))
            legendEffiAccPi0           = GetAndSetLegend2(0.15, 0.13, 0.53, 0.13+(5*textSizeLabelsRel),textSizeLabelsPixel);


        Style_t markerstylesPlot[5]                         ={34,20,33,29,22};
        Size_t markersizePlot[5]                        ={2.3,2.3,3.4,2.8,2.8};
        Color_t colorDataPlot[5]                        ={kBlack,kRed-3,kOrange+2,kCyan-3,kBlue-3};

        if(histo1[index]){
            DrawGammaSetMarker(histo1[index], markerstylesPlot[0], markersizePlot[0], colorDataPlot[0] , colorDataPlot[0]);
            histo1[index]->Draw("p,same,e");
            legendEffiAccPi0->AddEntry(histo1[index],"val. primary","p");
        }
        if(histo2[index]){
            DrawGammaSetMarker(histo2[index], markerstylesPlot[1], markersizePlot[1], colorDataPlot[1] , colorDataPlot[1]);
            histo2[index]->Draw("p,same,e");
            legendEffiAccPi0->AddEntry(histo2[index],"val. #pi^{0} from K^{0}_{S}","p");
        }
        if(histo3[index]){
            DrawGammaSetMarker(histo3[index], markerstylesPlot[2], markersizePlot[2], colorDataPlot[2] , colorDataPlot[2]);
            histo3[index]->Draw("p,same,e");
            legendEffiAccPi0->AddEntry(histo3[index],"val. #pi^{0} from #Lambda","p");
        }
        if(histo4[index]){
            DrawGammaSetMarker(histo4[index], markerstylesPlot[3], markersizePlot[3], colorDataPlot[3] , colorDataPlot[3]);
            histo4[index]->Draw("p,same,e");
            legendEffiAccPi0->AddEntry(histo4[index],"val. #pi^{0} from K^{0}_{L}","p");
        }
        if(histo5[index]){
            DrawGammaSetMarker(histo5[index], markerstylesPlot[4], markersizePlot[4], colorDataPlot[4] , colorDataPlot[4]);
            histo5[index]->Draw("p,same,e");
            legendEffiAccPi0->AddEntry(histo5[index],"val. #pi^{0} from Rest","p");
        }
        legendEffiAccPi0->Draw();
        if(!plotName.CompareTo("SecondaryYield")||!plotName.CompareTo("FracSecYield")){
            drawLatexAdd("ALICE simulation",0.92,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd(Form("pp, %s",nameMeasGlobal[index].Data()),0.92,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        }else{
            drawLatexAdd("ALICE simulation",0.15,0.92,textSizeLabelsRel,kFALSE,kFALSE);
            drawLatexAdd(Form("pp, %s",nameMeasGlobal[index].Data()),0.15,0.87,textSizeLabelsRel,kFALSE,kFALSE);
        }

    canvasDummy->Update();
    canvasDummy->Print(Form("%s/%s.%s",outputDir.Data(),plotName.Data(),suffix.Data()));

    }
void plotDCA(Int_t index, TH1D* histo1[],TH1D* histo2[],TH1D* histo3[],TH1D* histo4[],TH1D* histo5[], TH1D* histo6[], Double_t yMin, Double_t yMax, TString yLabel, Int_t logY, TString plotName)
    {


    Double_t textSizeLabelsPixel             = 55;
    Double_t textSizeLabelsRel      = 55./1200;

    TCanvas* canvasDummy       = new TCanvas("canvasDummy", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasDummy,  0.1, 0.01, 0.015, 0.095);
    if(logY)
        canvasDummy->SetLogy(1);
    canvasDummy->SetLogx(1);

        TH2F * histoDummy;
            histoDummy                = new TH2F("histoDummy", "histoDummy",1000, 0.23,  19.9, 1000, yMin, yMax );
        SetStyleHistoTH2ForGraphs( histoDummy, "#it{p}_{T} (GeV/#it{c})", yLabel.Data(),
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.0);//(#times #epsilon_{pur})
        histoDummy->GetYaxis()->SetLabelOffset(0.001);
        histoDummy->GetXaxis()->SetLabelOffset(-0.01);
        histoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
        histoDummy->DrawCopy();
        TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.63, 0.58, 0.93, 0.58+(6*textSizeLabelsRel),textSizeLabelsPixel);


        Style_t markerstylesPlot[6]                         ={34,20,33,29,20,20};
        Size_t markersizePlot[6]                        ={2.3,2.3,3.4,2.8,2.3,2.3};
        Color_t colorDataPlot[6]                        ={kBlack,kRed-3,kOrange+2,kCyan-3,kBlue-3,kPink};

        if(histo1[index]){
            DrawGammaSetMarker(histo1[index], markerstylesPlot[0], markersizePlot[0], colorDataPlot[0] , colorDataPlot[0]);
            histo1[index]->Draw("p,same,e");
            legendEffiAccPi0->AddEntry(histo1[index],"Category 1","p");
        }
        if(histo2[index]){
            DrawGammaSetMarker(histo2[index], markerstylesPlot[1], markersizePlot[1], colorDataPlot[1] , colorDataPlot[1]);
            histo2[index]->Draw("p,same,e");
            legendEffiAccPi0->AddEntry(histo2[index],"Category 2","p");
        }
        if(histo3[index]){
            DrawGammaSetMarker(histo3[index], markerstylesPlot[2], markersizePlot[2], colorDataPlot[2] , colorDataPlot[2]);
            histo3[index]->Draw("p,same,e");
            legendEffiAccPi0->AddEntry(histo3[index],"Category 3","p");
        }
        if(histo4[index]){
            DrawGammaSetMarker(histo4[index], markerstylesPlot[3], markersizePlot[3], colorDataPlot[3] , colorDataPlot[3]);
            histo4[index]->Draw("p,same,e");
            legendEffiAccPi0->AddEntry(histo4[index],"Category 4","p");
        }
        if(histo5[index]){
            DrawGammaSetMarker(histo5[index], markerstylesPlot[4], markersizePlot[4], colorDataPlot[4] , colorDataPlot[4]);
            histo5[index]->Draw("p,same,e");
            legendEffiAccPi0->AddEntry(histo5[index],"Category 5","p");
        }
        if(histo6[index]){
            DrawGammaSetMarker(histo6[index], markerstylesPlot[5], markersizePlot[5], colorDataPlot[5] , colorDataPlot[5]);
            histo6[index]->Draw("p,same,e");
            legendEffiAccPi0->AddEntry(histo6[index],"Category 6","p");
        }
        legendEffiAccPi0->Draw();
            drawLatexAdd("ALICE",0.92,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd(Form("pp, %s",nameMeasGlobal[index].Data()),0.92,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);


    canvasDummy->Update();
    canvasDummy->Print(Form("%s/%s.%s",outputDir.Data(),plotName.Data(),suffix.Data()));

    }
void plotCrossSection(TGraphAsymmErrors* csGraphs[],TGraphAsymmErrors* csGraphsSys[], Double_t yMin, Double_t yMax, TString meson, TString plotName)
    {


    Double_t textSizeLabelsPixel             = 50;
    Double_t textSizeLabelsRel      = 50./1200;

    TCanvas* canvasDummy       = new TCanvas("canvasDummy", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasDummy,  0.14, 0.01, 0.015, 0.095);
    canvasDummy->SetLogy(1);
    canvasDummy->SetLogx(1);
    Double_t xMax = 19.9;
    Double_t xMin = 0.23;
    if(!meson.CompareTo("Eta")){
        xMin= 0.29;
        xMax= 14;
    }
        TH2F * histoDummy;
            histoDummy                = new TH2F("histoDummy", "histoDummy",1000, xMin,  xMax, 1000, yMin, yMax );
        SetStyleHistoTH2ForGraphs( histoDummy, "#it{p}_{T} (GeV/#it{c})", "#it{E}#frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2}#it{c}^{3})",
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.4);//(#times #epsilon_{pur})
        if(!meson.CompareTo("Gamma"))
            SetStyleHistoTH2ForGraphs( histoDummy, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c})",
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.4);//(#times #epsilon_{pur})

        histoDummy->GetYaxis()->SetLabelOffset(0.001);
        histoDummy->GetXaxis()->SetLabelOffset(-0.01);
        histoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
        histoDummy->DrawCopy();

        TH1F * histoBlack               = (TH1F*)histoPi0InvCrossSection[1]->Clone("histoBlack") ;
        histoBlack->SetLineColor(kBlack);
        histoBlack->SetMarkerStyle(21) ;
        histoBlack->SetMarkerColor(kBlack) ;
        histoBlack->SetMarkerSize(1.) ;

        TGraphAsymmErrors * graphGrey   = (TGraphAsymmErrors*)csGraphsSys[1]->Clone("graphGrey") ;
        graphGrey   ->SetFillColor(kGray);
        graphGrey   ->SetLineColor(kGray);
        graphGrey   ->SetFillStyle(1001);


        Style_t markerstylesFULL[4]                     ={34,20,33,29};
        Size_t markersizeFULL[4]                        ={2.3,2.3,3.4,3.2};
        Color_t colorDataFULL[4]                        ={kRed+2,kBlue+2,kGreen+2,kMagenta+2};
        // systematics
        if(!plotName.CompareTo("InclusivePhotons")){
            for (int j=0;j<csGraphsSys[2]->GetN();j++){
                csGraphsSys[2]->GetY()[j] *= 1000;
                csGraphsSys[2]->GetEYhigh()[j] *= 1000;
                csGraphsSys[2]->GetEYlow()[j] *= 1000;
            }
            DrawGammaSetMarkerTGraphAsym(csGraphsSys[2], markerstylesFULL[2], markersizeFULL[2], colorDataFULL[2], colorDataFULL[2], widthLinesBoxes, kTRUE);
            for (int j=0;j<csGraphsSys[1]->GetN();j++){
                csGraphsSys[1]->GetY()[j] *= 100;
                csGraphsSys[1]->GetEYhigh()[j] *= 100;
                csGraphsSys[1]->GetEYlow()[j] *= 100;
            }
            DrawGammaSetMarkerTGraphAsym(csGraphsSys[1], markerstylesFULL[1], markersizeFULL[1], colorDataFULL[1], colorDataFULL[1], widthLinesBoxes, kTRUE);

            for (int j=0;j<csGraphsSys[3]->GetN();j++){
                csGraphsSys[3]->GetY()[j] *= 10;
                csGraphsSys[3]->GetEYhigh()[j] *= 10;
                csGraphsSys[3]->GetEYlow()[j] *= 10;
            }
            DrawGammaSetMarkerTGraphAsym(csGraphsSys[3], markerstylesFULL[3], markersizeFULL[3], colorDataFULL[3], colorDataFULL[3], widthLinesBoxes, kTRUE);

            DrawGammaSetMarkerTGraphAsym(csGraphsSys[0], markerstylesFULL[0], markersizeFULL[0], colorDataFULL[0], colorDataFULL[0], widthLinesBoxes, kTRUE);


            // statistics
            for (int j=0;j<csGraphs[2]->GetN();j++){
                csGraphs[2]->GetY()[j] *= 1000;
                csGraphs[2]->GetEYhigh()[j] *= 1000;
                csGraphs[2]->GetEYlow()[j] *= 1000;
            }
            DrawGammaSetMarkerTGraph(csGraphs[2], markerstylesFULL[2], 0.7*markersizeFULL[2], colorDataFULL[2] , colorDataFULL[2]);

            for (int j=0;j<csGraphs[1]->GetN();j++){
                csGraphs[1]->GetY()[j] *= 100;
                csGraphs[1]->GetEYhigh()[j] *= 100;
                csGraphs[1]->GetEYlow()[j] *= 100;
            }
            DrawGammaSetMarkerTGraph(csGraphs[1], markerstylesFULL[1], 0.7*markersizeFULL[1], colorDataFULL[1] , colorDataFULL[1]);

            for (int j=0;j<csGraphs[3]->GetN();j++){
                csGraphs[3]->GetY()[j] *= 10;
                csGraphs[3]->GetEYhigh()[j] *= 10;
                csGraphs[3]->GetEYlow()[j] *= 10;
            }
            DrawGammaSetMarkerTGraph(csGraphs[3], markerstylesFULL[3], 0.7*markersizeFULL[3], colorDataFULL[3] , colorDataFULL[3]);


            DrawGammaSetMarkerTGraph(csGraphs[0], markerstylesFULL[0], 0.7*markersizeFULL[0], colorDataFULL[0] , colorDataFULL[0]);
            csGraphsSys[2]     ->Draw("E2same");
            csGraphsSys[1]     ->Draw("E2same");
            csGraphsSys[3]     ->Draw("E2same");
            csGraphsSys[0]     ->Draw("E2same");
            csGraphs[2]->Draw("p,same,z");
            csGraphs[1]->Draw("p,same,z");
            csGraphs[3]->Draw("p,same,z");
            csGraphs[0]->Draw("p,same,z");

            TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.18, 0.13, 0.43, 0.13+(4.5*textSizeLabelsRel),textSizeLabelsPixel);
                    legendEffiAccPi0->AddEntry(csGraphsSys[2],Form("%s x10^{3}",nameMeasGlobal[2].Data()),"pf");
                    legendEffiAccPi0->AddEntry(csGraphsSys[1],Form("%s x10^{2}",nameMeasGlobal[1].Data()),"pf");
                    legendEffiAccPi0->AddEntry(csGraphsSys[3],Form("%s x10",nameMeasGlobal[3].Data()),"pf");
                    legendEffiAccPi0->AddEntry(csGraphsSys[0],nameMeasGlobal[0].Data(),"pf");


            legendEffiAccPi0->Draw();
        }else{
            for (int j=0;j<csGraphsSys[2]->GetN();j++){
                csGraphsSys[2]->GetY()[j] *= 100;
                csGraphsSys[2]->GetEYhigh()[j] *= 100;
                csGraphsSys[2]->GetEYlow()[j] *= 100;
            }
            DrawGammaSetMarkerTGraphAsym(csGraphsSys[2], markerstylesFULL[2], markersizeFULL[2], colorDataFULL[2], colorDataFULL[2], widthLinesBoxes, kTRUE);

            for (int j=0;j<csGraphsSys[1]->GetN();j++){
                csGraphsSys[1]->GetY()[j] *= 10;
                csGraphsSys[1]->GetEYhigh()[j] *= 10;
                csGraphsSys[1]->GetEYlow()[j] *= 10;
            }
            DrawGammaSetMarkerTGraphAsym(csGraphsSys[1], markerstylesFULL[1], markersizeFULL[1], colorDataFULL[1], colorDataFULL[1], widthLinesBoxes, kTRUE);

            DrawGammaSetMarkerTGraphAsym(csGraphsSys[0], markerstylesFULL[0], markersizeFULL[0], colorDataFULL[0], colorDataFULL[0], widthLinesBoxes, kTRUE);


            // statistics
            for (int j=0;j<csGraphs[2]->GetN();j++){
                csGraphs[2]->GetY()[j] *= 100;
                csGraphs[2]->GetEYhigh()[j] *= 100;
                csGraphs[2]->GetEYlow()[j] *= 100;
            }
            DrawGammaSetMarkerTGraph(csGraphs[2], markerstylesFULL[2], 0.7*markersizeFULL[2], colorDataFULL[2] , colorDataFULL[2]);

            for (int j=0;j<csGraphs[1]->GetN();j++){
                csGraphs[1]->GetY()[j] *= 10;
                csGraphs[1]->GetEYhigh()[j] *= 10;
                csGraphs[1]->GetEYlow()[j] *= 10;
            }
            DrawGammaSetMarkerTGraph(csGraphs[1], markerstylesFULL[1], 0.7*markersizeFULL[1], colorDataFULL[1] , colorDataFULL[1]);


            DrawGammaSetMarkerTGraph(csGraphs[0], markerstylesFULL[0], 0.7*markersizeFULL[0], colorDataFULL[0] , colorDataFULL[0]);
            csGraphsSys[2]     ->Draw("E2same");
            csGraphsSys[1]     ->Draw("E2same");
            csGraphsSys[0]     ->Draw("E2same");
            csGraphs[2]->Draw("p,same,z");
            csGraphs[1]->Draw("p,same,z");
            csGraphs[0]->Draw("p,same,z");

            TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.18, 0.13, 0.43, 0.13+(3.5*textSizeLabelsRel),textSizeLabelsPixel);
                    legendEffiAccPi0->AddEntry(csGraphsSys[2],Form("%s x10^{2}",nameMeasGlobal[2].Data()),"pf");
                    legendEffiAccPi0->AddEntry(csGraphsSys[1],Form("%s x10",nameMeasGlobal[1].Data()),"pf");
                    legendEffiAccPi0->AddEntry(csGraphsSys[0],nameMeasGlobal[0].Data(),"pf");


            legendEffiAccPi0->Draw();
        }


        // TLegend* legendPi0Err2 = GetAndSetLegend2(0.72, 0.72, 0.98, 0.72+(2*textSizeLabelsRel),0.85*textSizeLabelsPixel);
        // legendPi0Err2->AddEntry(histoBlack, "stat. Err.","ple");
        // legendPi0Err2->AddEntry(graphGrey,  "syst. Err.","f");
        // legendPi0Err2->Draw();

        drawLatexAdd("ALICE",0.9,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        drawLatexAdd("pp, PCM",0.9,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        if(meson.CompareTo("Gamma")){
            if(!meson.CompareTo("Pi0"))
                drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.9,0.82,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            else
                drawLatexAdd("#eta #rightarrow #gamma#gamma",0.9,0.82,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        }
        histoDummy->Draw("sameaxis");
    canvasDummy->Update();
    canvasDummy->Print(Form("%s/%s.%s",outputDir.Data(),plotName.Data(),suffix.Data()));
    /*
        if(!plotName.CompareTo("InclusivePhotons")){
            TCanvas* canvasDummyRatio       = new TCanvas("canvasDummyRatio", "", 200, 10, 1200, 1100);  // gives the page size
            DrawGammaCanvasSettings( canvasDummyRatio,  0.14, 0.01, 0.015, 0.095);
    //         canvasDummyRatio->SetLogy(1);
            canvasDummyRatio->cd();
            canvasDummyRatio->SetLogx(1);
            TH2F * histoDummyRatio;
                histoDummyRatio                = new TH2F("histoDummyRatio", "histoDummyRatio",1000, xMin,  xMax, 1000, 1.0, 1.8 );
                SetStyleHistoTH2ForGraphs( histoDummyRatio, "#it{p}_{T} (GeV/#it{c})","#gamma_{incl.}^{8 TeV}/#gamma_{incl.}^{7 TeV}",
                                    0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.4);//(#times #epsilon_{pur})

            histoDummyRatio->GetYaxis()->SetLabelOffset(0.001);
            histoDummyRatio->GetXaxis()->SetLabelOffset(-0.01);
            histoDummyRatio->GetXaxis()->SetMoreLogLabels(kTRUE);
            histoDummyRatio->DrawCopy();
    //             TH1D* ratiohisto = (TH1D*) csGraphs[2]->Clone("ratioinclusivehisto");
            // csGraphsSys[2]->RemovePoint(csGraphsSys[2]->GetN());
            // for (int j=0;j<csGraphsSys[2]->GetN();j++){
            //     csGraphsSys[2]->GetY()[j] *= 1/1000;
            //     csGraphsSys[2]->GetEYhigh()[j] *= 1/1000;
            //     csGraphsSys[2]->GetEYlow()[j] *= 1/1000;
            // }
            for (int j=0;j<csGraphsSys[2]->GetN();j++){
                csGraphsSys[2]->GetY()[j] /= 1000;
                csGraphsSys[2]->GetEYhigh()[j] /= 1000;
                csGraphsSys[2]->GetEYlow()[j] /= 1000;
            }
            for (int j=0;j<csGraphs[2]->GetN();j++){
                csGraphs[2]->GetY()[j] /= 1000;
                csGraphs[2]->GetEYhigh()[j] /= 1000;
                csGraphs[2]->GetEYlow()[j] /= 1000;
            }
            for (int j=0;j<csGraphsSys[1]->GetN();j++){
                csGraphsSys[1]->GetY()[j] /= 100;
                csGraphsSys[1]->GetEYhigh()[j] /= 100;
                csGraphsSys[1]->GetEYlow()[j] /= 100;
            }
            for (int j=0;j<csGraphs[1]->GetN();j++){
                csGraphs[1]->GetY()[j] /= 100;
                csGraphs[1]->GetEYhigh()[j] /= 100;
                csGraphs[1]->GetEYlow()[j] /= 100;
            }
            // csGraphsSys[2]->SetPoint(20,14,-1);
            TGraphAsymmErrors* ratiohistosys = CalculateAsymGraphRatioToGraph(csGraphsSys[2],csGraphsSys[1]);
            // ratiohistosys->RemovePoint(ratiohistosys->GetN()-1);
            DrawGammaSetMarkerTGraphAsym(ratiohistosys, markerstylesFULL[1], markersizeFULL[1], colorDataFULL[1], colorDataFULL[1], widthLinesBoxes, kTRUE);
            ratiohistosys->Draw("E2same");

            TGraphAsymmErrors* ratiohisto = CalculateAsymGraphRatioToGraph(csGraphs[2],csGraphs[1]);
            DrawGammaSetMarkerTGraph(ratiohisto, markerstylesFULL[1], 0.7*markersizeFULL[1], colorDataFULL[1] , colorDataFULL[1]);
            ratiohisto->Draw("p,same,z");
            canvasDummyRatio->Update();
            canvasDummyRatio->Print(Form("%s/%s.%s",outputDir.Data(),"InclusivePhotonsRatio",suffix.Data()));

        }
        */
    }

void plotEtaToPi0(TGraphAsymmErrors* csGraphs[],TGraphAsymmErrors* csGraphsSys[], Double_t yMin, Double_t yMax, TString yTitle, TString plotName){

    // ***************************************************************************************************************
    // ******************************* Plotting eta/pi0 ratio for single measurements ********************************
    // ***************************************************************************************************************
    cout << "PLOTTING: Eta/Pi0 ratio" << endl;
    Double_t textSizeLabelsPixel                 = 54;
    TCanvas* canvasEtatoPi0combo       = new TCanvas("canvasEtatoPi0combo","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasEtatoPi0combo, 0.1, 0.01, 0.01, 0.125);
    canvasEtatoPi0combo->SetLogx();

        Double_t textsizeLabelsEtaToPi0 = 0;
        Double_t textsizeFacEtaToPi0    = 0;
        if (canvasEtatoPi0combo->XtoPixel(canvasEtatoPi0combo->GetX2()) <canvasEtatoPi0combo->YtoPixel(canvasEtatoPi0combo->GetY1()) ){
            textsizeLabelsEtaToPi0      = (Double_t)textSizeLabelsPixel/canvasEtatoPi0combo->XtoPixel(canvasEtatoPi0combo->GetX2()) ;
            textsizeFacEtaToPi0         = (Double_t)1./canvasEtatoPi0combo->XtoPixel(canvasEtatoPi0combo->GetX2()) ;
        } else {
            textsizeLabelsEtaToPi0      = (Double_t)textSizeLabelsPixel/canvasEtatoPi0combo->YtoPixel(canvasEtatoPi0combo->GetY1());
            textsizeFacEtaToPi0         = (Double_t)1./canvasEtatoPi0combo->YtoPixel(canvasEtatoPi0combo->GetY1());
        }
        Double_t xMin = 0.33;
        Double_t xMax = 15;
    if(!plotName.CompareTo("InclusiveRatio")||!plotName.CompareTo("CocktailRatio")||!plotName.CompareTo("DoubleRatio")){
        xMin=0.23;
        xMax=19.9;
    }
    TH2F * histo2DEtatoPi0combo;
    histo2DEtatoPi0combo               = new TH2F("histo2DEtatoPi0combo","histo2DEtatoPi0combo",1000,xMin,xMax,1000,yMin,yMax    );
    SetStyleHistoTH2ForGraphs(histo2DEtatoPi0combo, "#it{p}_{T} (GeV/#it{c})",yTitle, 0.85*textsizeLabelsEtaToPi0, textsizeLabelsEtaToPi0,
                              0.85*textsizeLabelsEtaToPi0,1.1*textsizeLabelsEtaToPi0, 0.9, 0.65, 510, 510);
    histo2DEtatoPi0combo->GetXaxis()->SetMoreLogLabels();
    histo2DEtatoPi0combo->GetXaxis()->SetLabelOffset(-0.01);
    // histo2DEtatoPi0combo->GetYaxis()->SetRangeUser(0.0,1.05);
    histo2DEtatoPi0combo->Draw();
        // plotting systematics graphs

    if(!plotName.CompareTo("DoubleRatio")||!plotName.CompareTo("CocktailRatio"))
            DrawGammaLines(0.23, 19.9 , 1., 1.,0.1, kGray);

    Double_t markerScale = 1.3;
        for (Int_t i = 0; i < 4; i++){
            if(csGraphsSys[i]){
                DrawGammaSetMarkerTGraphAsym(csGraphsSys[i], markerstyles[i], markerScale*markersize[i], colorData[i] , colorData[i], widthLinesBoxes, kTRUE);
                if(plotName.CompareTo("CocktailRatio"))
                  csGraphsSys[i]->Draw("E2same");
            }else if(graphEtaToPi0ShiftSys[i]){
                DrawGammaSetMarkerTGraphAsym(graphEtaToPi0ShiftSys[i], markerstyles[i], markerScale*markersize[i], colorData[i] , colorData[i], widthLinesBoxes, kTRUE);
                graphEtaToPi0ShiftSys[i]->Draw("E2same");
            }
        }

        // plotting statistics graphs
        for (Int_t i = 0; i < 4; i++){
            if(csGraphs[i]){
                DrawGammaSetMarkerTGraphAsym(csGraphs[i], markerstyles[i], markerScale*markersize[i], colorData[i] , colorData[i]);
                csGraphs[i]->Draw("p,same,e");
            }
            else if(graphEtaToPi0ShiftStat[i]){
                DrawGammaSetMarkerTGraphAsym(graphEtaToPi0ShiftStat[i], markerstyles[i], markerScale*markersize[i], colorData[i] , colorData[i]);
                graphEtaToPi0ShiftStat[i]->Draw("p,same,e");
            }
        }

        TLegend* legendEtaToPi0;
        if(!plotName.CompareTo("InclusiveRatio"))
            legendEtaToPi0 = GetAndSetLegend2(0.15, 0.15, 0.38, 0.15+(textsizeLabelsEtaToPi0*4*0.9), textSizeLabelsPixel);
        else if(!plotName.CompareTo("CocktailRatio"))
            legendEtaToPi0 = GetAndSetLegend2(0.25, 0.73-textsizeLabelsEtaToPi0, 0.48, 0.73+(textsizeLabelsEtaToPi0*4*0.9), 0.9*textSizeLabelsPixel);
        else
            legendEtaToPi0 = GetAndSetLegend2(0.67, 0.15, 0.9, 0.15+(textsizeLabelsEtaToPi0*4*0.9), textSizeLabelsPixel);
        for (Int_t i = 2; i > -1; i--){
            if(csGraphsSys[i]){
                legendEtaToPi0->AddEntry(csGraphsSys[i],nameMeasGlobal[i],"pf");
                if(i==1 && csGraphsSys[3])
                    legendEtaToPi0->AddEntry(csGraphsSys[3],nameMeasGlobal[3],"pf");
            }
            else if(graphEtaToPi0ShiftSys[i]){
                legendEtaToPi0->AddEntry(graphEtaToPi0ShiftSys[i],nameMeasGlobal[i],"pf");
                if(i==1 && csGraphsSys[3])
                    legendEtaToPi0->AddEntry(csGraphsSys[3],nameMeasGlobal[3],"pf");
            }
        }
        legendEtaToPi0->Draw();



        if(!plotName.CompareTo("InclusiveRatio")){
            drawLatexAdd("ALICE",0.92,0.92,0.85*textsizeLabelsEtaToPi0,kFALSE,kFALSE,kTRUE);
            // drawLatexAdd("#eta/#pi^{0} #rightarrow #gamma#gamma",0.92, 0.92-(1*textsizeLabelsEtaToPi0*0.85),0.85*textsizeLabelsEtaToPi0,kFALSE,kFALSE,kTRUE);
            // drawLatexAdd("#gamma's rec. with PCM",0.13, 0.92-(2*textsizeLabelsEtaToPi0*0.9),textsizeLabelsEtaToPi0,kFALSE);
            drawLatexAdd("pp, PCM",0.92, 0.92-(1*textsizeLabelsEtaToPi0*0.9),0.85*textsizeLabelsEtaToPi0,kFALSE,kFALSE,kTRUE);
        }else if(!plotName.CompareTo("CocktailRatio")){
            drawLatexAdd("ALICE",0.92,0.92,0.85*textsizeLabelsEtaToPi0,kFALSE,kFALSE,kTRUE);
            drawLatexAdd("pp, PCM",0.92, 0.92-(1*textsizeLabelsEtaToPi0*0.85),0.85*textsizeLabelsEtaToPi0,kFALSE,kFALSE,kTRUE);
            // drawLatexAdd("#eta/#pi^{0} #rightarrow #gamma#gamma",0.92, 0.92-(1*textsizeLabelsEtaToPi0*0.85),0.85*textsizeLabelsEtaToPi0,kFALSE,kFALSE,kTRUE);
            // drawLatexAdd("#gamma's rec. with PCM",0.13, 0.92-(2*textsizeLabelsEtaToPi0*0.9),textsizeLabelsEtaToPi0,kFALSE);
            // drawLatexAdd("pp, PCM",0.92, 0.92-(1*textsizeLabelsEtaToPi0*0.9),0.85*textsizeLabelsEtaToPi0,kFALSE,kFALSE,kTRUE);
        } else {
            drawLatexAdd("ALICE",0.13, 0.92,0.85*textsizeLabelsEtaToPi0,kFALSE);
            drawLatexAdd("#eta/#pi^{0} #rightarrow #gamma#gamma",0.13, 0.92-(1*textsizeLabelsEtaToPi0*0.85),0.85*textsizeLabelsEtaToPi0,kFALSE);
            // drawLatexAdd("#gamma's rec. with PCM",0.13, 0.92-(2*textsizeLabelsEtaToPi0*0.9),textsizeLabelsEtaToPi0,kFALSE);
            drawLatexAdd("pp, PCM",0.13, 0.92-(2*textsizeLabelsEtaToPi0*0.9),0.85*textsizeLabelsEtaToPi0,kFALSE);
        }

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/%s.%s",outputDir.Data(),plotName.Data(), suffix.Data()));
}

void plotChargedComparison(TGraphAsymmErrors* neutralStat[],TGraphAsymmErrors* neutralSys[],TGraphAsymmErrors* chargedStat[],TGraphAsymmErrors* chargedSys[], Double_t yMin, Double_t yMax, TString yTitle, TString plotName){

        Double_t xSection8TeVx         = 55.8*1e-3;
        Double_t xSection7TeVx         = 62.22*1e-3;
        Double_t xSection900GeVx         = 47.78*1e-3;
    Double_t recalcBarnx             = 1e12;

    TFile* input8TeVfile = new TFile("ThesisQAInputAllEnergies/CombinedResultsPaperPP8TeV_2017_01_10.root");
    TDirectory* directoryPi08TeV                     = (TDirectory*)input8TeVfile->Get("Pi08TeV");
    TFile* input7TeVfile = new TFile("ThesisQAInputAllEnergies/CombinedResultsPaperPP7TeV_2017_02_25.root");
    TDirectory* directoryPi07TeV                     = (TDirectory*)input7TeVfile->Get("Pi07TeV");
    TFile* input900GeVfile = new TFile("ThesisQAInputAllEnergies/CombinedResultsPaperPP900GeV_2017_02_25.root");
    TDirectory* directoryPi0900GeV                     = (TDirectory*)input900GeVfile->Get("Pi0900GeV");

    neutralStat[2]     = (TGraphAsymmErrors*)directoryPi08TeV->Get("graphInvCrossSectionPi0Comb8TeVAStatErr");
    neutralSys[2]     = (TGraphAsymmErrors*)directoryPi08TeV->Get("graphInvCrossSectionPi0Comb8TeVASysErr");
    neutralStat[1]     = (TGraphAsymmErrors*)directoryPi07TeV->Get("graphCombPi0InvCrossSectionStatPCMEMCPHOS");
    neutralSys[1]     = (TGraphAsymmErrors*)directoryPi07TeV->Get("graphCombPi0InvCrossSectionSysPCMEMCPHOS");
    neutralStat[0]     = (TGraphAsymmErrors*)directoryPi0900GeV->Get("graphCombPi0InvCrossSectionStatPCMPHOS");
    neutralSys[0]     = (TGraphAsymmErrors*)directoryPi0900GeV->Get("graphCombPi0InvCrossSectionSysPCMPHOS");

    for (int j=0;j<neutralStat[1]->GetN();j++){
            neutralStat[1]->GetY()[j] *= 1/(xSection7TeVx*recalcBarnx);
            neutralStat[1]->GetEYhigh()[j] *=  1/(xSection7TeVx*recalcBarnx);
            neutralStat[1]->GetEYlow()[j] *=  1/(xSection7TeVx*recalcBarnx);
        }
    for (int j=0;j<neutralSys[1]->GetN();j++){
            neutralSys[1]->GetY()[j] *=  1/(xSection7TeVx*recalcBarnx);
            neutralSys[1]->GetEYhigh()[j] *=  1/(xSection7TeVx*recalcBarnx);
            neutralSys[1]->GetEYlow()[j] *= 1/(xSection7TeVx*recalcBarnx);
        }
    for (int j=0;j<neutralStat[2]->GetN();j++){
            neutralStat[2]->GetY()[j] *= 1/(xSection8TeVx*recalcBarnx);
            neutralStat[2]->GetEYhigh()[j] *=  1/(xSection8TeVx*recalcBarnx);
            neutralStat[2]->GetEYlow()[j] *=  1/(xSection8TeVx*recalcBarnx);
        }
    for (int j=0;j<neutralSys[2]->GetN();j++){
            neutralSys[2]->GetY()[j] *=  1/(xSection8TeVx*recalcBarnx);
            neutralSys[2]->GetEYhigh()[j] *=  1/(xSection8TeVx*recalcBarnx);
            neutralSys[2]->GetEYlow()[j] *= 1/(xSection8TeVx*recalcBarnx);
        }

    Style_t markerstylesFULL[3]                         ={34,20,33};
    Size_t markersizeFULL[3]                        ={2.3,2.3,3.4};
    Color_t colorDataFULL[3]                        ={kRed+2,kBlue+2,kGreen+2};


    Int_t testhisto = 0;

    //     cout << "charged stat: " << endl;
    // for (int j=0;j<chargedStat[testhisto]->GetN();j++){
    //         cout << j << "\t" << chargedStat[testhisto]->GetX()[j] << "\t" << chargedStat[testhisto]->GetEXlow()[j] << "\t" << chargedStat[testhisto]->GetY()[j] << "\t" << chargedStat[testhisto]->GetEYhigh()[j] << endl;
    //     }
    // cout << "neutral stat: " << endl;
    // for (int j=0;j<neutralStat[testhisto]->GetN();j++){
    //         cout << j << "\t" << neutralStat[testhisto]->GetX()[j] << "\t" << neutralStat[testhisto]->GetEXlow()[j] << "\t" << neutralStat[testhisto]->GetY()[j] << "\t" << neutralStat[testhisto]->GetEYhigh()[j] << endl;
    //     }
    //     cout << "charged sys: " << endl;
    // for (int j=0;j<chargedSys[testhisto]->GetN();j++){
    //         cout << j << "\t" << chargedSys[testhisto]->GetX()[j] << "\t" << chargedSys[testhisto]->GetEXlow()[j] << "\t" << chargedSys[testhisto]->GetY()[j] << "\t" << chargedSys[testhisto]->GetEYhigh()[j] << endl;
    //     }
    // cout << "neutral sys: " << endl;
    // for (int j=0;j<neutralSys[testhisto]->GetN();j++){
    //         cout << j << "\t" << neutralSys[testhisto]->GetX()[j] << "\t" << neutralSys[testhisto]->GetEXlow()[j] << "\t" << neutralSys[testhisto]->GetY()[j] << "\t" << neutralSys[testhisto]->GetEYhigh()[j] << endl;
    //     }

    TGraphErrors* graphChargedPionStatPt[4] = {NULL, NULL, NULL, NULL};
    TGraphErrors* graphChargedPionSysPt[4] = {NULL, NULL, NULL, NULL};
    TGraphErrors* graphPi0YieldStatRebinned[4] = {NULL, NULL, NULL, NULL};
    TGraphErrors* graphPi0YieldSysRebinned[4] = {NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphChargedRatio[4];

    //charged ratio 8 TeV
    double xvalue8[] = { 0.35, 0.45, 0.55, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 3.9, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 17, 19 };
  double xerrlow8[] = { 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1};
  double xerrhigh8[] = { 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1 };
  double yvalue8[] = { 2.7990888134, 1.7489632083, 1.4217942916, 1.2011610685, 1.159936424, 1.1150930719, 1.0593225873, 1.055222081, 1.0550366723, 1.034590305, 1.0315893211, 1.0277596473, 1.0628357862, 1.0796797731, 1.0159511322, 1.06676613, 1.0247474547, 1.0956403179, 1.060344443, 1.0761547589, 1.0594697939, 1.0582718752, 1.0756200019, 1.0406366366, 0.9687301165, 0.9927929636, 1.0340942671, 0.998576203, 1.0222866626, 1.0672124727, 1.0185608416, 1.0638358707, 1.0836731508, 1.1648966151, 1.1131937444, 1.1023324602, 1.364049935  };
  double yerrlow8[] = {0.8930700793, 0.2821731653, 0.1859285338, 0.1386565879, 0.1091139116, 0.0943297868, 0.0726186381, 0.0773348449, 0.0752540141, 0.0735848841, 0.0795975459, 0.0809531955, 0.0926670666, 0.100380364, 0.0960686142, 0.1207798111, 0.1167041964, 0.1263137186, 0.1216740786, 0.1237886418, 0.1214890693, 0.122055556, 0.1254152343, 0.1241613593, 0.1323322081, 0.1358679811, 0.1064463377, 0.102614983, 0.1333092058, 0.1429513545, 0.1385226939, 0.1503214205, 0.1574792584, 0.179412332, 0.1784284604, 0.1421370969, 0.2502848113 };
  double yerrhigh8[] = {0.8930700793, 0.2821731653, 0.1859285338, 0.1386565879, 0.1091139116, 0.0943297868, 0.0726186381, 0.0773348449, 0.0752540141, 0.0735848841, 0.0795975459, 0.0809531955, 0.0926670666, 0.100380364, 0.0960686142, 0.1207798111, 0.1167041964, 0.1263137186, 0.1216740786, 0.1237886418, 0.1214890693, 0.122055556, 0.1254152343, 0.1241613593, 0.1323322081, 0.1358679811, 0.1064463377, 0.102614983, 0.1333092058, 0.1429513545, 0.1385226939, 0.1503214205, 0.1574792584, 0.179412332, 0.1784284604, 0.1421370969, 0.2502848113};
  int numpoints8 = 37;
  graphChargedRatio[2] = new TGraphAsymmErrors(numpoints8, xvalue8, yvalue8, xerrlow8, xerrhigh8, yerrlow8, yerrhigh8);
   //charged ratio 7 TeV

    double xvalue[] = { 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 3.9, 4.5, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 15, 17, 19 };
  double xerrlow[] = { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 1 };
  double xerrhigh[] = { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 1, 2.5 };
  double yvalue[] = { 1.5273321273, 1.3111161613, 1.3008043299, 1.3243456842, 1.2849443505, 1.3579264888, 1.3279891048, 1.3413767951, 1.3517032065, 1.3291600832, 1.2755194674, 1.269901704, 1.2423049934, 1.2349146128, 1.2244323667, 1.2478470823, 1.256122208, 1.2484454809, 1.2270401803, 1.2367246501, 1.2521691547, 1.2190456315, 1.2670823531, 1.2469543137, 1.2295417452, 1.2557074719, 1.2385189238, 1.2907898298, 1.2385395503, 1.2517063152, 1.3146614665, 1.2110698468, 1.2306774421, 1.1755271363, 1.2732085399, 1.2199398271, 1.2799817627, 1.2925853063, 1.320708072, 1.1655989741, 1.2445714492, 1.1686875622  };
  double yerrlow[] = {0.3688455418, 0.1537116558, 0.1320668759, 0.1393712628, 0.1309519484, 0.1238274755, 0.1089425067, 0.1154368962, 0.1152109002, 0.1036342263, 0.0986959511, 0.0976647495, 0.0947355126, 0.0937187883, 0.0930076026, 0.0951211859, 0.0962782843, 0.0741945666, 0.0743441455, 0.0885400853, 0.0963645529, 0.0947365658, 0.1186637445, 0.1170374635, 0.1149022509, 0.1193845652, 0.1168508256, 0.095723645, 0.1170164228, 0.1189977794, 0.1589438842, 0.1495062839, 0.1218411837, 0.119183507, 0.1316273218, 0.1355035073, 0.1475556308, 0.1557265238, 0.1640688181, 0.1301834204, 0.1684801552, 0.1783764276 };
  double yerrhigh[] = { 0.3688455418, 0.1537116558, 0.1320668759, 0.1393712628, 0.1309519484, 0.1238274755, 0.1089425067, 0.1154368962, 0.1152109002, 0.1036342263, 0.0986959511, 0.0976647495, 0.0947355126, 0.0937187883, 0.0930076026, 0.0951211859, 0.0962782843, 0.0741945666, 0.0743441455, 0.0885400853, 0.0963645529, 0.0947365658, 0.1186637445, 0.1170374635, 0.1149022509, 0.1193845652, 0.1168508256, 0.095723645, 0.1170164228, 0.1189977794, 0.1589438842, 0.1495062839, 0.1218411837, 0.119183507, 0.1316273218, 0.1355035073, 0.1475556308, 0.1557265238, 0.1640688181, 0.1301834204, 0.1684801552, 0.1783764276};
  int numpoints = 42;
  graphChargedRatio[1] = new TGraphAsymmErrors(numpoints, xvalue, yvalue, xerrlow, xerrhigh, yerrlow, yerrhigh);

  for (int j=0;j<graphChargedRatio[1]->GetN();j++){
            graphChargedRatio[1]->GetY()[j] *=  0.849;
            graphChargedRatio[1]->GetEYhigh()[j] *=  0.849;
            graphChargedRatio[1]->GetEYlow()[j] *= 0.849;
        }

  double xvalue900[] = { 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.8 };
  double xerrlow900[] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2 };
  double xerrhigh900[] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2 };
  double yvalue900[] = {1.6057453767, 1.3575522945, 1.3500473328, 1.3080317232, 1.3695986108, 1.3618072345, 1.3006519784 };
  double yerrlow900[] = {0.4207057842, 0.2206787688, 0.1706113723, 0.1475450737, 0.1540479183, 0.1553724475, 0.1453547031 };
  double yerrhigh900[] = {0.4207057842, 0.2206787688, 0.1706113723, 0.1475450737, 0.1540479183, 0.1553724475, 0.1453547031};
  int numpoints900 = 7;
  graphChargedRatio[0] = new TGraphAsymmErrors(numpoints900, xvalue900, yvalue900, xerrlow900, xerrhigh900, yerrlow900, yerrhigh900);

  for (int j=0;j<graphChargedRatio[0]->GetN();j++){
            graphChargedRatio[0]->GetY()[j] *=  0.908;
            graphChargedRatio[0]->GetEYhigh()[j] *=  0.908;
            graphChargedRatio[0]->GetEYlow()[j] *= 0.908;
        }
    // ***************************************************************************************************************
    // ******************************* Plotting eta/pi0 ratio for single measurements ********************************
    // ***************************************************************************************************************
    cout << "PLOTTING: charged ratio" << endl;
    Double_t textSizeLabelsPixel                 = 54;
    Double_t textSizeLabelsRel                 = 54/900;
    TCanvas* cavnaschargedcomp       = new TCanvas("cavnaschargedcomp","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( cavnaschargedcomp, 0.09, 0.01, 0.01, 0.125);
    cavnaschargedcomp->SetLogx();
    // cavnaschargedcomp->SetLogy();

        Double_t textsizelabelschargedcomp = 0;
        Double_t textsizefacchargedcomp    = 0;
        if (cavnaschargedcomp->XtoPixel(cavnaschargedcomp->GetX2()) <cavnaschargedcomp->YtoPixel(cavnaschargedcomp->GetY1()) ){
            textsizelabelschargedcomp      = (Double_t)textSizeLabelsPixel/cavnaschargedcomp->XtoPixel(cavnaschargedcomp->GetX2()) ;
            textsizefacchargedcomp         = (Double_t)1./cavnaschargedcomp->XtoPixel(cavnaschargedcomp->GetX2()) ;
        } else {
            textsizelabelschargedcomp      = (Double_t)textSizeLabelsPixel/cavnaschargedcomp->YtoPixel(cavnaschargedcomp->GetY1());
            textsizefacchargedcomp         = (Double_t)1./cavnaschargedcomp->YtoPixel(cavnaschargedcomp->GetY1());
        }

    TH2F * histo2Dchargedcomp;
    histo2Dchargedcomp               = new TH2F("histo2Dchargedcomp","histo2Dchargedcomp",1000,0.23,25.9,1000,0.72,2.25    );
    SetStyleHistoTH2ForGraphs(histo2Dchargedcomp, "#it{p}_{T} (GeV/#it{c})","#pi^{0}/#pi^{#pm}", 0.85*textsizelabelschargedcomp, textsizelabelschargedcomp,
                              0.85*textsizelabelschargedcomp,1.1*textsizelabelschargedcomp, 0.9, 0.65, 510, 510);
    histo2Dchargedcomp->GetXaxis()->SetMoreLogLabels();
    histo2Dchargedcomp->GetXaxis()->SetLabelOffset(-0.01);
    // histo2Dchargedcomp->GetYaxis()->SetRangeUser(0.0,1.05);
    histo2Dchargedcomp->Draw();
        // plotting systematics graphs
    Double_t markerScale = 1.3;

                DrawGammaSetMarkerTGraphAsym(graphChargedRatio[1], markerstyles[1], markersize[1], colorData[1] , colorData[1], widthLinesBoxes, kTRUE);
                graphChargedRatio[1]->Draw("E2same,p");
                // DrawGammaSetMarkerTGraphAsym(graphChargedRatio[2], markerstyles[2], markersize[2], colorData[2] , colorData[2], widthLinesBoxes, kTRUE);
                // graphChargedRatio[2]->Draw("E2same,p");
                DrawGammaSetMarkerTGraphAsym(graphChargedRatio[0], markerstyles[0], markersize[0], colorData[0] , colorData[0], widthLinesBoxes, kTRUE);
                graphChargedRatio[0]->Draw("E2same,p");

        DrawGammaLines(0.23, 25.9 , 1., 1.,0.1, kGray);


        TGraphAsymmErrors * graphBlackBox               = (TGraphAsymmErrors*)graphChargedRatio[1]->Clone("graphBlackBox") ;
        DrawGammaSetMarkerTGraphAsym(graphBlackBox, markerstyles[1], markersize[1], kBlack , kBlack, widthLinesBoxes, kTRUE);

        TLegend* legendXsectionPaper    = GetAndSetLegend2(0.66, 0.72, 0.96, 0.82+0.05, 0.9*textSizeLabelsPixel);
        // TLegend* legendXsectionPaper    = GetAndSetLegend2(0.66, 0.68, 0.96, 0.78+0.05*2, 0.9*textSizeLabelsPixel);
        legendXsectionPaper->SetNColumns(1);
        legendXsectionPaper->SetMargin(0.2);
        // legendXsectionPaper->AddEntry(graphChargedRatio[2],Form("%s",nameMeasGlobal[2].Data()),"pf");
        legendXsectionPaper->AddEntry(graphChargedRatio[1],Form("%s",nameMeasGlobal[1].Data()),"pf");
        legendXsectionPaper->AddEntry(graphChargedRatio[0],Form("%s",nameMeasGlobal[0].Data()),"pf");
        legendXsectionPaper->Draw();

        TLegend* legenuncertaint    = GetAndSetLegend2(0.15, 0.15, 0.44, 0.15+0.065, 0.9*textSizeLabelsPixel);
        legenuncertaint->SetNColumns(1);
        legenuncertaint->SetMargin(0.2);
        // legenuncertaint->AddEntry(graphBlackBox, "Uncertainties: stat. #oplus sys.", "pf");
        legenuncertaint->AddEntry(graphBlackBox, "stat. #oplus sys.", "pf");
        legenuncertaint->Draw();

        drawLatexAdd("ALICE",0.93,0.91,textsizelabelschargedcomp,kFALSE,kFALSE,kTRUE);
            // drawLatexAdd(Form("pp, %s",nameMeasGlobal[index].Data()),0.92,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);

    cavnaschargedcomp->Update();
    cavnaschargedcomp->SaveAs(Form("%s/ChargedComparison.%s",outputDir.Data(), suffix.Data()));
}

void plotGammaSecondaries(Int_t index, TH1D* histo1[],TH1D* histo2[],TH1D* histo3[],TH1D* histo4[], Double_t yMin, Double_t yMax, TString yLabel, Int_t logY, TString plotName)
    {


    Double_t textSizeLabelsPixel             = 55;
    Double_t textSizeLabelsRel      = 55./1200;

    TCanvas* canvasDummy       = new TCanvas("canvasDummy", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasDummy,  0.1, 0.01, 0.015, 0.095);
    if(!plotName.CompareTo("SecGammaConvProb"))
        DrawGammaCanvasSettings( canvasDummy,  0.11, 0.01, 0.015, 0.095);
    if(logY)
        canvasDummy->SetLogy(1);
    canvasDummy->SetLogx(1);

        TH2F * histoDummy;
            histoDummy                = new TH2F("histoDummy", "histoDummy",1000, 0.23,  19.9, 1000, yMin, yMax );
        SetStyleHistoTH2ForGraphs( histoDummy, "#it{p}_{T} (GeV/#it{c})", yLabel.Data(),
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.06);//(#times #epsilon_{pur})
        if(!plotName.CompareTo("SecGammaConvProb"))
            SetStyleHistoTH2ForGraphs( histoDummy, "#it{p}_{T} (GeV/#it{c})", yLabel.Data(),
                                    0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.15);//(#times #epsilon_{pur})
        histoDummy->GetYaxis()->SetLabelOffset(0.001);
        histoDummy->GetXaxis()->SetLabelOffset(-0.01);
        histoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
        histoDummy->DrawCopy();
        TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.55, 0.77, 0.92, 0.77+(4*textSizeLabelsRel),textSizeLabelsPixel);
        if(histo4[index])
            legendEffiAccPi0           = GetAndSetLegend2(0.55, 0.74, 0.92, 0.755+(4.3*textSizeLabelsRel),textSizeLabelsPixel);
        if(!plotName.CompareTo("SecGammaYield"))
            legendEffiAccPi0           = GetAndSetLegend2(0.15, 0.13, 0.53, 0.13+(5*textSizeLabelsRel),textSizeLabelsPixel);


        Style_t markerstylesPlot[5]                         ={20,33,29,22};
        Size_t markersizePlot[5]                        ={2.3,3.4,2.8,2.8};
        Color_t colorDataPlot[5]                        ={kRed-3,kOrange+2,kCyan-3,kBlue-3};

        if(histo1[index]){
            DrawGammaSetMarker(histo1[index], markerstylesPlot[0], markersizePlot[0], colorDataPlot[0] , colorDataPlot[0]);
            histo1[index]->Draw("p,same,e");
            legendEffiAccPi0->AddEntry(histo1[index],"sec. #gamma from K^{0}_{S}","p");
        }
        if(histo4[index]){
            DrawGammaSetMarker(histo4[index], markerstylesPlot[3], markersizePlot[3], colorDataPlot[3] , colorDataPlot[3]);
            histo4[index]->Draw("p,same,e");
            legendEffiAccPi0->AddEntry(histo4[index],"sec. #gamma from Rest","p");
        }
        if(histo2[index]){
            DrawGammaSetMarker(histo2[index], markerstylesPlot[1], markersizePlot[1], colorDataPlot[1] , colorDataPlot[1]);
            histo2[index]->Draw("p,same,e");
            legendEffiAccPi0->AddEntry(histo2[index],"sec. #gamma from K^{0}_{L}","p");
        }
        if(histo3[index]){
            DrawGammaSetMarker(histo3[index], markerstylesPlot[2], markersizePlot[2], colorDataPlot[2] , colorDataPlot[2]);
            histo3[index]->Draw("p,same,e");
            legendEffiAccPi0->AddEntry(histo3[index],"sec. #gamma from #Lambda","p");
        }


        legendEffiAccPi0->Draw();
        if(!plotName.CompareTo("SecGammaYield")||!plotName.CompareTo("FracSecYield")){
            drawLatexAdd("ALICE simulation",0.92,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd(Form("pp, %s",nameMeasGlobal[index].Data()),0.92,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        }else{
            drawLatexAdd("ALICE simulation",0.15,0.92,textSizeLabelsRel,kFALSE,kFALSE);
            drawLatexAdd(Form("pp, %s",nameMeasGlobal[index].Data()),0.15,0.87,textSizeLabelsRel,kFALSE,kFALSE);
        }

    canvasDummy->Update();
    canvasDummy->Print(Form("%s/%s_%d.%s",outputDir.Data(),plotName.Data(),index,suffix.Data()));

    }

void plotGammaPileup(Int_t index, TH1D* histo1[],TH1D* histo2[],TH1D* histo3[],TH1D* histo4[], Double_t yMin, Double_t yMax, TString yLabel, Int_t logY, TString plotName)
    {


    Double_t textSizeLabelsPixel             = 55;
    Double_t textSizeLabelsRel      = 55./1200;

    TCanvas* canvasDummy       = new TCanvas("canvasDummy", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasDummy,  0.1, 0.01, 0.015, 0.095);
    if(logY)
        canvasDummy->SetLogy(1);

        TH2F * histoDummy;
            histoDummy                = new TH2F("histoDummy", "histoDummy",1000, -5.99,  5.99, 1000, yMin, yMax );
        SetStyleHistoTH2ForGraphs( histoDummy, "DCA#it{z} (cm)", yLabel.Data(),
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.06);//(#times #epsilon_{pur})
        histoDummy->GetYaxis()->SetLabelOffset(0.001);
        // histoDummy->GetXaxis()->SetLabelOffset(-0.01);
        histoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
        histoDummy->DrawCopy();
        TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.13, 0.81-2*textSizeLabelsRel, 0.42, 0.81+(3*textSizeLabelsRel),0.9*textSizeLabelsPixel);
            // legendEffiAccPi0           = GetAndSetLegend2(0.55, 0.74, 0.92, 0.755+(4.3*textSizeLabelsRel),textSizeLabelsPixel);
            // legendEffiAccPi0           = GetAndSetLegend2(0.15, 0.13, 0.53, 0.13+(5*textSizeLabelsRel),textSizeLabelsPixel);


        Style_t markerstylesPlot[6]                         ={20,24,29,22,30,26};
        Size_t markersizePlot[5]                        ={1.5,2,2.,2.8};
        Color_t colorDataPlot[6]                        ={kBlack,kGray+2,kCyan-3,kBlue-3,kRed-2,kRed+2};

        Double_t maxValueDCAData    = 0;
    if(HistDCAZUnderMesonCat_1_MesonPt_0304[index]&&fHistDCAZUnderMesonBGEstimateCat_1_MesonPt_0304[index]){
        for(Int_t i=0;i<HistDCAZUnderMesonCat_1_MesonPt_0304[index]->GetNbinsX();i++){
            if(HistDCAZUnderMesonCat_1_MesonPt_0304[index]->GetBinContent(i)>maxValueDCAData)
                maxValueDCAData = HistDCAZUnderMesonCat_1_MesonPt_0304[index]->GetBinContent(i);
        }
    }
    cout << "##############################################################\n\tmaxValueDCAData: " << maxValueDCAData << endl;

        if(HistDCAZUnderMesonCat_1_MesonPt_0304MC[index]&&HistDCAZUnderMesonCat_1_MesonPt_0304MCTRUE[index]){
            DrawGammaSetMarker(HistDCAZUnderMesonCat_1_MesonPt_0304MC[index], markerstylesPlot[3], markersizePlot[1], colorDataPlot[4] , colorDataPlot[4]);
            HistDCAZUnderMesonCat_1_MesonPt_0304MC[index]->Sumw2();
            maxValueDCAData    = 0;
            for(Int_t i=0;i<HistDCAZUnderMesonCat_1_MesonPt_0304MC[index]->GetNbinsX();i++){
                if(HistDCAZUnderMesonCat_1_MesonPt_0304MC[index]->GetBinContent(i)>maxValueDCAData)
                    maxValueDCAData = HistDCAZUnderMesonCat_1_MesonPt_0304MC[index]->GetBinContent(i);
            }
            HistDCAZUnderMesonCat_1_MesonPt_0304MC[index]->Scale(1/maxValueDCAData);
            // HistDCAZUnderMesonCat_1_MesonPt_0304MC[index]->Scale(1/HistDCAZUnderMesonCat_1_MesonPt_0304MC[index]->GetMaximum());
            // HistDCAZUnderMesonCat_1_MesonPt_0304MC[index]->Scale((Double_t)1/fNEventsDCAMC[index]);
            HistDCAZUnderMesonCat_1_MesonPt_0304MC[index]->Draw("p,same,e");
            legendEffiAccPi0->AddEntry(HistDCAZUnderMesonCat_1_MesonPt_0304MC[index],"MC DCAz cat. 1","p");
        }
        if(HistDCAZUnderMesonCat_1_MesonPt_0304[index]&&fHistDCAZUnderMesonBGEstimateCat_1_MesonPt_0304[index]){
            // signal
            DrawGammaSetMarker(HistDCAZUnderMesonCat_1_MesonPt_0304[index], markerstylesPlot[0], markersizePlot[0], colorDataPlot[0] , colorDataPlot[0]);
            cout << "##############################################################\n\t\tNEVENTS: " << fNEventsDCA[index] << endl;
            HistDCAZUnderMesonCat_1_MesonPt_0304[index]->Sumw2();
            maxValueDCAData = 0;
            for(Int_t i=0;i<HistDCAZUnderMesonCat_1_MesonPt_0304[index]->GetNbinsX();i++){
                if(HistDCAZUnderMesonCat_1_MesonPt_0304[index]->GetBinContent(i)>maxValueDCAData)
                    maxValueDCAData = HistDCAZUnderMesonCat_1_MesonPt_0304[index]->GetBinContent(i);
            }
            HistDCAZUnderMesonCat_1_MesonPt_0304[index]->Scale((Double_t)1/maxValueDCAData);
            // HistDCAZUnderMesonCat_1_MesonPt_0304[index]->Scale((Double_t)1/fNEventsDCA[index]);
            // HistDCAZUnderMesonCat_1_MesonPt_0304[index]->Scale(1/HistDCAZUnderMesonCat_1_MesonPt_0304[index]->GetMaximum());
            legendEffiAccPi0->AddEntry(HistDCAZUnderMesonCat_1_MesonPt_0304[index],"data DCAz cat. 1","p");
            HistDCAZUnderMesonCat_1_MesonPt_0304[index]->Draw("p,same,e");

            //background
            DrawGammaSetMarker(fHistDCAZUnderMesonBGEstimateCat_1_MesonPt_0304[index], markerstylesPlot[3], markersizePlot[1], colorDataPlot[3] , colorDataPlot[3]);
            fHistDCAZUnderMesonBGEstimateCat_1_MesonPt_0304[index]->Sumw2();
            fHistDCAZUnderMesonBGEstimateCat_1_MesonPt_0304[index]->Scale((Double_t)1/maxValueDCAData);
            // fHistDCAZUnderMesonBGEstimateCat_1_MesonPt_0304[index]->Scale((Double_t)1/fNEventsDCA[index]);
            legendEffiAccPi0->AddEntry(fHistDCAZUnderMesonBGEstimateCat_1_MesonPt_0304[index],"Estimated Pileup","p");
            fHistDCAZUnderMesonBGEstimateCat_1_MesonPt_0304[index]->Draw("p,same");

            //background subtracted signal
            TH1D* histoSubtracted = (TH1D*) HistDCAZUnderMesonCat_1_MesonPt_0304[index]->Clone("subtractedDCA");
            histoSubtracted->Sumw2();
            histoSubtracted->Add(fHistDCAZUnderMesonBGEstimateCat_1_MesonPt_0304[index],-1);
            DrawGammaSetMarker(histoSubtracted, markerstylesPlot[4], markersizePlot[0], kCyan , kCyan);
            legendEffiAccPi0->AddEntry(histoSubtracted,"Pileup Subtracted","p");
            // histoSubtracted->Scale(1/histoSubtracted->GetMaximum());
            histoSubtracted->Draw("pl,same");
        }

      //
    //     if(histo4[index]){
    //         DrawGammaSetMarker(histo4[index], markerstylesPlot[4], markersizePlot[0], colorDataPlot[5] , colorDataPlot[5]);
    //         histo4[index]->Draw("p,same,e");
    //         legendEffiAccPi0->AddEntry(histo4[index],"MC True DCAz cat. 1","p");
    //     }
    //     if(ESDGammaDCAzAll3MC[index]){
    //         DrawGammaSetMarker(ESDGammaDCAzAll3MC[index], markerstylesPlot[1], markersizePlot[1], colorDataPlot[4] , colorDataPlot[4]);
    //         ESDGammaDCAzAll3MC[index]->Draw("p,same,e");
    //         legendEffiAccPi0->AddEntry(ESDGammaDCAzAll3MC[index],"MC True DCAz cat. 3","p");
    //     }
    //     if(histo1[index]){
    //         DrawGammaSetMarker(histo1[index], markerstylesPlot[0], markersizePlot[0], colorDataPlot[0] , colorDataPlot[0]);
    //         histo1[index]->Draw("p,same,e");
    //         legendEffiAccPi0->AddEntry(histo1[index],"data DCAz cat. 1","p");
    //     }
    //    if(ESDGammaDCAzAll3[index]){
    //         DrawGammaSetMarker(ESDGammaDCAzAll3[index], markerstylesPlot[1], markersizePlot[0], colorDataPlot[1] , colorDataPlot[1]);
    //         ESDGammaDCAzAll3[index]->Draw("p,same,e");
    //         legendEffiAccPi0->AddEntry(ESDGammaDCAzAll3[index],"data DCAz cat. 3","p");
    //     }
    //   if(histo2[index]){
    //         DrawGammaSetMarker(histo2[index], markerstylesPlot[3], markersizePlot[1], colorDataPlot[3] , colorDataPlot[3]);
    //         histo2[index]->Draw("p,same,e");
    //         legendEffiAccPi0->AddEntry(histo2[index],"Estimated Pileup","p");
    //     }
    //     if(histo3[index]){
    //         DrawGammaSetMarker(histo3[index], markerstylesPlot[2], markersizePlot[2], colorDataPlot[2] , colorDataPlot[2]);
    //         histo3[index]->Draw("p,same,e");
    //         legendEffiAccPi0->AddEntry(histo3[index],"Pileup Subtracted","p");
    //     }
      //

        legendEffiAccPi0->Draw();
        // if(!plotName.CompareTo("SecGammaYield")||!plotName.CompareTo("FracSecYield")){
            // drawLatexAdd("0.4 < #it{p}_{T} < 0.6 GeV/#it{c}",0.92,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd("0.5 < #it{p}_{T} < 0.6 GeV/#it{c}",0.92,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd("ALICE",0.92,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd(Form("pp, %s",nameMeasGlobal[index].Data()),0.92,0.82,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
        // }else{
        //     drawLatexAdd("ALICE this thesis",0.15,0.92,textSizeLabelsRel,kFALSE,kFALSE);
        //     drawLatexAdd(Form("pp, %s",nameMeasGlobal[index].Data()),0.15,0.87,textSizeLabelsRel,kFALSE,kFALSE);
        // }

    canvasDummy->Update();
    canvasDummy->Print(Form("%s/%s_%d_0.5-0.6.%s",outputDir.Data(),plotName.Data(),index,suffix.Data()));

    }

void plotGammaDCAMultiple(Int_t index, TH1D* histo1[],TH1D* histo2[],TH1D* histo3[],TH1D* histo4[], Double_t yMin, Double_t yMax, TString yLabel, Int_t logY, TString plotName)
    {
     TFile* inputfileDCAPastFuture[2];
     TFile* inputfileDCAPastFutureMC[2];
    inputfileDCAPastFuture[0]                      = new TFile("/home/nschmidt/AnalysisSoftware/ThesisQAInputAllEnergies/PileupStudies/8TeV_20170329Trainoutput_withDCA/Pi0_Data_GammaConvV1DCATestAnalysed.root");
    inputfileDCAPastFuture[1]                      = new TFile("/home/nschmidt/AnalysisSoftware/ThesisQAInputAllEnergies/PileupStudies/8TeV_newWithDCA/Pi0_Data_GammaConvV1DCATestAnalysed.root");
    inputfileDCAPastFutureMC[0]                      = new TFile("/home/nschmidt/AnalysisSoftware/ThesisQAInputAllEnergies/PileupStudies/8TeV_20170329Trainoutput_withDCA/Pi0_MC_GammaConvV1DCATestAnalysed.root");
    inputfileDCAPastFutureMC[1]                      = new TFile("/home/nschmidt/AnalysisSoftware/ThesisQAInputAllEnergies/PileupStudies/8TeV_newWithDCA/Pi0_MC_GammaConvV1DCATestAnalysed.root");
    TH1D* fHistDCAdata[2];
    TH1D* fHistDCAMC[2];
    TH1D* fHistDCAMCTrue[2];
    TH1D* fHistDCABackground[2];
    TH1D* fHistEventQuality[2];
    TH1D* fHistEventQualityMC[2];
    Int_t fNEventsDCAPF[2];
    Int_t fNEventsDCAPFMC[2];
    for(Int_t i=0;i<2;i++){
        fHistDCAdata[i]              = (TH1D*)inputfileDCAPastFuture[i]->Get("HistDCAZUnderMesonCat_1_MesonPt_0.50-0.60");
        fHistDCABackground[i]   = (TH1D*)inputfileDCAPastFuture[i]->Get("fHistDCAZUnderMesonBGEstimateCat_1_MesonPt_0.50-0.60");
        fHistDCAMCTrue[i]        = (TH1D*)inputfileDCAPastFutureMC[i]->Get("HistDCAZTruePrimaryMesonGammaGammaCat_1_MesonPt_0.50-0.60");
        fHistDCAMC[i]              = (TH1D*)inputfileDCAPastFutureMC[i]->Get("HistDCAZUnderMesonCat_1_MesonPt_0.50-0.60");
        fHistEventQuality[i]                          = (TH1D*)inputfileDCAPastFuture[i]->Get("NEvents");
        if(fHistEventQuality[i])
            fNEventsDCAPF[i]                               = GetNEvents(fHistEventQuality[i]);
        else
            fNEventsDCAPF[i]                          = -1;

        fHistEventQualityMC[i]                          = (TH1D*)inputfileDCAPastFutureMC[i]->Get("NEvents");
        if(fHistEventQualityMC[i])
            fNEventsDCAPFMC[i]                               = GetNEvents(fHistEventQualityMC[i]);
        else
            fNEventsDCAPFMC[i]                          = -1;
    }

    Double_t textSizeLabelsPixel             = 55;
    Double_t textSizeLabelsRel      = 55./1200;

    TCanvas* canvasDummy       = new TCanvas("canvasDummy", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasDummy,  0.1, 0.01, 0.015, 0.095);
    if(logY)
        canvasDummy->SetLogy(1);

        TH2F * histoDummy;
            histoDummy                = new TH2F("histoDummy", "histoDummy",1000, -5.99,  5.99, 1000, yMin, yMax );
        SetStyleHistoTH2ForGraphs( histoDummy, "DCA#it{z} (cm)", yLabel.Data(),
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.06);//(#times #epsilon_{pur})
        histoDummy->GetYaxis()->SetLabelOffset(0.001);
        // histoDummy->GetXaxis()->SetLabelOffset(-0.01);
        histoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
        histoDummy->DrawCopy();
        // TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.13, 0.81-4*textSizeLabelsRel, 0.42, 0.81+(3*textSizeLabelsRel),0.7*textSizeLabelsPixel);
        TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.13, 0.81-1*textSizeLabelsRel, 0.42, 0.81+(3*textSizeLabelsRel),0.7*textSizeLabelsPixel);
        Style_t markerstylesPlot[6]                         ={20,24,29,29,30,26};
        Size_t markersizePlot[5]                        ={1.5,2,2.,2.8};
        Color_t colorDataPlot[6]                        ={kBlack,kGray+2,kCyan-3,kBlue-3,kRed-2,kRed+2};


        Double_t maxValueDCAData[2]    = {0,0};
        for(Int_t i=0;i<2;i++){
            if(fHistDCAdata[i]&&fHistDCABackground[i]){
                for(Int_t j=0;j<fHistDCAdata[i]->GetNbinsX();j++){
                    if(fHistDCAdata[i]->GetBinContent(j)>maxValueDCAData[i])
                        maxValueDCAData[i] = fHistDCAdata[i]->GetBinContent(j);
                }
            }
            cout << "##############################################################\n\tmaxValueDCAData[i]: " << maxValueDCAData[i] << endl;
        }

       for(Int_t i=0;i<1;i++){
            if(fHistDCAMC[i]&&fHistDCAMCTrue[i]){
                DrawGammaSetMarker(fHistDCAMC[i], markerstylesPlot[3], markersizePlot[1], colorDataPlot[i+4] , colorDataPlot[i+4]);
                fHistDCAMC[i]->Sumw2();
                maxValueDCAData[i]    = 0;

           for(Int_t j=0;j<fHistDCAMC[i]->GetNbinsX();j++){
                    if(fHistDCAMC[i]->GetBinContent(j)>maxValueDCAData[i])
                        maxValueDCAData[i] = fHistDCAMC[i]->GetBinContent(j);
                }

           fHistDCAMC[i]->Scale(1/maxValueDCAData[i]);
                fHistDCAMC[i]->Draw("p,same,e");
                legendEffiAccPi0->AddEntry(fHistDCAMC[i],"MC DCAz cat. 1","p");
            }
        }

       for(Int_t i=0;i<2;i++){
            if(fHistDCAdata[i]&&fHistDCABackground[i]){
                // signal
                DrawGammaSetMarker(fHistDCAdata[i], markerstylesPlot[0], markersizePlot[0], colorDataPlot[i] , colorDataPlot[i]);
                cout << "##############################################################\n\t\tNEVENTS: " << fNEventsDCA[i] << endl;
                fHistDCAdata[i]->Sumw2();
                maxValueDCAData[i] = 0;
                for(Int_t j=0;j<fHistDCAdata[i]->GetNbinsX();j++){
                    if(fHistDCAdata[i]->GetBinContent(j)>maxValueDCAData[i])
                        maxValueDCAData[i] = fHistDCAdata[i]->GetBinContent(j);
                }
                fHistDCAdata[i]->Scale((Double_t)1/maxValueDCAData[i]);
                // if(i==0)
                // legendEffiAccPi0->AddEntry(fHistDCAdata[i],"data DCAz cat. 1","p");
                // else
                // legendEffiAccPi0->AddEntry(fHistDCAdata[i],"PF data DCAz cat. 1","p");
                // fHistDCAdata[i]->Draw("p,same,e");


                //background
                DrawGammaSetMarker(fHistDCABackground[i], markerstylesPlot[3], markersizePlot[1], colorDataPlot[3]+i*3 , colorDataPlot[3]+i*3);
                fHistDCABackground[i]->Sumw2();
                fHistDCABackground[i]->Scale((Double_t)1/maxValueDCAData[i]);
                // if(i==0)
                // legendEffiAccPi0->AddEntry(fHistDCABackground[i],"Estimated Pileup","p");
                // else
                // legendEffiAccPi0->AddEntry(fHistDCABackground[i],"PF Estimated Pileup","p");
                // fHistDCABackground[i]->Draw("p,same");


                //background subtracted signal
                TH1D* histoSubtracted = (TH1D*) fHistDCAdata[i]->Clone("subtractedDCA");
                histoSubtracted->Sumw2();
                histoSubtracted->Add(fHistDCABackground[i],-1);
                DrawGammaSetMarker(histoSubtracted, markerstylesPlot[4], markersizePlot[0], kCyan+3*i , kCyan+3*i);
                if(i==0)
                legendEffiAccPi0->AddEntry(histoSubtracted,"Pileup Subtracted","p");
                else
                legendEffiAccPi0->AddEntry(histoSubtracted,"PF Pileup Subtracted","p");
                histoSubtracted->Draw("pl,same");
            }
        }



        legendEffiAccPi0->Draw();
            drawLatexAdd("0.4 < #it{p}_{T} < 0.5 GeV/#it{c}",0.92,0.92,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd("ALICE",0.92,0.87,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);
            drawLatexAdd(Form("pp, %s",nameMeasGlobal[2].Data()),0.92,0.82,textSizeLabelsRel,kFALSE,kFALSE,kTRUE);


    canvasDummy->Update();
    canvasDummy->Print(Form("%s/%s.%s",outputDir.Data(),plotName.Data(),suffix.Data()));

    }

void plotPhotonCombBG(Int_t index, TH1D* histo1[][17], Double_t yMin, Double_t yMax, TString yLabel, Int_t logY, TString plotName)
    {


    Double_t textSizeLabelsPixel             = 55;
    Double_t textSizeLabelsRel      = 55./1200;

    TCanvas* canvasDummy       = new TCanvas("canvasDummy", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasDummy,  0.1, 0.01, 0.015, 0.095);
    if(!plotName.CompareTo("SecGammaConvProb"))
        DrawGammaCanvasSettings( canvasDummy,  0.11, 0.01, 0.015, 0.095);
    if(logY)
        canvasDummy->SetLogy(1);
    canvasDummy->SetLogx(1);

        TH2F * histoDummy;
            histoDummy                = new TH2F("histoDummy", "histoDummy",1000, 0.23,  19.9, 1000, yMin, yMax );
        SetStyleHistoTH2ForGraphs( histoDummy, "#it{p}_{T} (GeV/#it{c})", yLabel.Data(),
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.06);//(#times #epsilon_{pur})
        if(!plotName.CompareTo("SecGammaConvProb"))
            SetStyleHistoTH2ForGraphs( histoDummy, "#it{p}_{T} (GeV/#it{c})", yLabel.Data(),
                                    0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.15);//(#times #epsilon_{pur})
        histoDummy->GetYaxis()->SetLabelOffset(0.001);
        histoDummy->GetXaxis()->SetLabelOffset(-0.01);
        histoDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
        histoDummy->DrawCopy();
        TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.15, 0.82, 0.92, 0.82+(3*textSizeLabelsRel),textSizeLabelsPixel*0.8,3);

            for(Int_t i = 0;i<9;i++){
                    DrawGammaSetMarker(histo1[index][i], markersCombinatorics[i], 1.5, colorsCombinatorics[i], colorsCombinatorics[i]);
                    histo1[index][i]->DrawCopy("e1same");
                legendEffiAccPi0->AddEntry(histo1[index][i],combinatorics[i]);
            }

        legendEffiAccPi0->Draw();

            drawLatexAdd("ALICE simulation",0.15,0.77,textSizeLabelsRel,kFALSE,kFALSE);
            drawLatexAdd(Form("pp, %s",nameMeasGlobal[index].Data()),0.15,0.72,textSizeLabelsRel,kFALSE,kFALSE);

    canvasDummy->Update();
    canvasDummy->Print(Form("%s/%s_%d.%s",outputDir.Data(),plotName.Data(),index,suffix.Data()));

    }

void plotTriggeredRawYields(){
    // SINGLE PLOTS FOR SLIDES
    TFile* inputFilePi0[3];
    TFile* inputFileEta[3];
    inputFilePi0[0] = new TFile("/home/nschmidt/AnalysisSoftware/ThesisQAInputAllEnergies/8tevtriggered/Pi0_data_GammaConvV1WithoutCorrection_00010113_00200009227300008250404000_0152103500000000.root");
    inputFilePi0[1] = new TFile("/home/nschmidt/AnalysisSoftware/ThesisQAInputAllEnergies/8tevtriggered/Pi0_data_GammaConvV1WithoutCorrection_00052113_00200009227300008250404000_0152103500000000.root");
    inputFilePi0[2] = new TFile("/home/nschmidt/AnalysisSoftware/ThesisQAInputAllEnergies/8tevtriggered/Pi0_data_GammaConvV1WithoutCorrection_00081113_00200009227300008250404000_0152103500000000.root");
    inputFileEta[0] = new TFile("/home/nschmidt/AnalysisSoftware/ThesisQAInputAllEnergies/8tevtriggered/Eta_data_GammaConvV1WithoutCorrection_00010113_00200009227300008250404000_0152103500000000.root");
    inputFileEta[1] = new TFile("/home/nschmidt/AnalysisSoftware/ThesisQAInputAllEnergies/8tevtriggered/Eta_data_GammaConvV1WithoutCorrection_00052113_00200009227300008250404000_0152103500000000.root");
    inputFileEta[2] = new TFile("/home/nschmidt/AnalysisSoftware/ThesisQAInputAllEnergies/8tevtriggered/Eta_data_GammaConvV1WithoutCorrection_00081113_00200009227300008250404000_0152103500000000.root");

    TGraphAsymmErrors* csGraphsPi0[3];
    TGraphAsymmErrors* csGraphsEta[3];
    for(Int_t i=0;i<3;i++)
    {
        csGraphsPi0[i] = new TGraphAsymmErrors((TH1D*)inputFilePi0[i]->Get("histoYieldMesonPerEvent"));
        csGraphsEta[i] = new TGraphAsymmErrors((TH1D*)inputFileEta[i]->Get("histoYieldMesonPerEvent"));
    }


    Double_t arrayBoundariesX1_XSec2[2];
    Double_t arrayBoundariesY1_XSec2[3];
    Double_t relativeMarginsXXSec2[3];
    Double_t relativeMarginsYXSec2[3];
    
    Double_t textSizeLabelsPixel = 48;
    ReturnCorrectValuesForCanvasScaling(1250,1300, 1, 1,0.135, 0.005, 0.003,0.09,arrayBoundariesX1_XSec2,arrayBoundariesY1_XSec2,relativeMarginsXXSec2,relativeMarginsYXSec2);

    TCanvas* canvasInvSectionPaperSingle      = new TCanvas("canvasInvSectionPaperSingle","",0,0,1250,1250);  // gives the page size
    DrawGammaCanvasSettings( canvasInvSectionPaperSingle,  0.13, 0.02, 0.03, 0.06);

    TPad* padInvSectionSpecSingle             = new TPad("padInvSectionSpecSingle", "", arrayBoundariesX1_XSec2[0], arrayBoundariesY1_XSec2[3], arrayBoundariesX1_XSec2[1], arrayBoundariesY1_XSec2[0],-1, -1, -2);
    DrawGammaPadSettings( padInvSectionSpecSingle, relativeMarginsXXSec2[0], relativeMarginsXXSec2[2], relativeMarginsYXSec2[0], relativeMarginsYXSec2[2]);
    padInvSectionSpecSingle->Draw();
    Double_t marginXSec                 = relativeMarginsXXSec2[0]*1250;
    Double_t textsizeLabelsXSecUp       = 0;
    Double_t textsizeFacXSecUp          = 0;
    if (padInvSectionSpecSingle->XtoPixel(padInvSectionSpecSingle->GetX2()) < padInvSectionSpecSingle->YtoPixel(padInvSectionSpecSingle->GetY1())){
        textsizeLabelsXSecUp            = (Double_t)textSizeLabelsPixel/padInvSectionSpecSingle->XtoPixel(padInvSectionSpecSingle->GetX2()) ;
        textsizeFacXSecUp               = (Double_t)1./padInvSectionSpecSingle->XtoPixel(padInvSectionSpecSingle->GetX2()) ;
    } else {
        textsizeLabelsXSecUp            = (Double_t)textSizeLabelsPixel/padInvSectionSpecSingle->YtoPixel(padInvSectionSpecSingle->GetY1());
        textsizeFacXSecUp               = (Double_t)1./padInvSectionSpecSingle->YtoPixel(padInvSectionSpecSingle->GetY1());
    }

    padInvSectionSpecSingle->cd();
    padInvSectionSpecSingle->SetLogy(1);
    padInvSectionSpecSingle->SetLogx(1);
    TH2F * histo2DXSectionPi0;
    // histo2DXSectionPi0          = new TH2F("histo2DXSectionPi0","histo2DXSectionPi0",11000,0.23,50.,1000,6,9e11);
    histo2DXSectionPi0          = new TH2F("histo2DXSectionPi0","histo2DXSectionPi0",11000,0.23,29,1000,3.1e-8,1.9e-3);
    SetStyleHistoTH2ForGraphs(histo2DXSectionPi0, "#it{p}_{T} (GeV/#it{c})","Raw Yield per Event",0.035,0.04, 0.035,0.04, 0.9,1.45);
    histo2DXSectionPi0->GetXaxis()->SetMoreLogLabels();
    histo2DXSectionPi0->GetXaxis()->SetNoExponent(kTRUE);

        SetStyleHistoTH2ForGraphs(histo2DXSectionPi0, "#it{p}_{T} (GeV/#it{c})","Raw Yield per Event",
                                0.85*textsizeLabelsXSecUp,textsizeLabelsXSecUp, 0.85*textsizeLabelsXSecUp, textsizeLabelsXSecUp, 1,0.2/(textsizeFacXSecUp*marginXSec));
        histo2DXSectionPi0->GetXaxis()->SetMoreLogLabels();
        histo2DXSectionPi0->GetXaxis()->SetLabelOffset(+0.01);
        histo2DXSectionPi0->Draw();
    histo2DXSectionPi0->Draw();

    Style_t markerstylesFULL[3]                         ={34,20,33};
    Size_t markersizeFULL[3]                        ={2.3,2.3,3.4};
    Color_t colorDataFULL[3]                        ={kGreen+2,kYellow+2,kOrange+2};


    for(Int_t i=0;i<3;i++)
    {
        DrawGammaSetMarkerTGraph(csGraphsPi0[i], markerstylesFULL[i], 0.9*markersizeFULL[i], colorDataFULL[i] , colorDataFULL[i]);
        DrawGammaSetMarkerTGraph(csGraphsEta[i], markerstylesFULL[i], 0.9*markersizeFULL[i], colorDataFULL[i] , colorDataFULL[i]);
    }
        

        csGraphsPi0[2]->Draw("p,same,z");
        csGraphsPi0[1]->Draw("p,same,z");
        csGraphsPi0[0]->Draw("p,same,z");

        histo2DXSectionPi0->Draw("same,axis");

        Double_t rightalignDouble = 0.93;
        drawLatexAdd("ALICE",rightalignDouble,0.91,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE);
        drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",rightalignDouble,0.86,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE);


        TLegend* legendXsectionPaperSingle    = GetAndSetLegend2(0.17, 0.13, 0.5, 0.23+0.05*2, textSizeLabelsPixel);
        legendXsectionPaperSingle->SetNColumns(1);
        legendXsectionPaperSingle->SetMargin(0.2);
        legendXsectionPaperSingle->AddEntry(csGraphsPi0[0],"MinBias","pl");
        legendXsectionPaperSingle->AddEntry(csGraphsPi0[1],"EMC7","pl");
        legendXsectionPaperSingle->AddEntry(csGraphsPi0[2],"EGA","pl");
        legendXsectionPaperSingle->Draw();
        //
        // drawLatexAdd("x10^{2}",0.24,0.94,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE,colorDataFULL[2]);
        // drawLatexAdd("x10^{1}",0.24,0.835,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE,colorDataFULL[1]);

    canvasInvSectionPaperSingle->Print(Form("%s/Pi0_TriggeredRawYields.%s",outputDir.Data(),suffix.Data()));
}

void plotCrossSectionFULL(TGraphAsymmErrors* csGraphs[],TGraphAsymmErrors* csGraphsSys[], Double_t yMin, Double_t yMax, TString meson, TString plotName)
{
      Bool_t usecombinedspectra = kTRUE;
      Bool_t plotTheorycurves   = kFALSE;
    Double_t xSection900GeVx         = 47.78*1e-3;
    Double_t recalcBarnx             = 1e12;

    TFile* input8TeVfile = new TFile("ThesisQAInputAllEnergies/CombinedResultsPaperPP8TeV_2017_07_13.root");
    TDirectory* directoryPi08TeV                     = (TDirectory*)input8TeVfile->Get("Pi08TeV");
    TFile* input7TeVfile = new TFile("ThesisQAInputAllEnergies/CombinedResultsPaperPP7TeV_2017_08_28.root");
    TDirectory* directoryPi07TeV                     = (TDirectory*)input7TeVfile->Get("Pi07TeV");
    TFile* input276TeVfile = new TFile("ThesisQAInputAllEnergies/2.76TeV/CombinedResultsPaperPP2760GeV_2017_03_01_FrediV2Clusterizer.root");
    TDirectory* directoryPi0276TeV                     = (TDirectory*)input276TeVfile->Get("Pi02.76TeV");
    TFile* input900GeVfile = new TFile("ThesisQAInputAllEnergies/CombinedResultsPaperPP900GeV_2017_07_18.root");
    TDirectory* directoryPi0900GeV                     = (TDirectory*)input900GeVfile->Get("Pi0900GeV");
    if(usecombinedspectra){
      csGraphs[3]     = (TGraphAsymmErrors*)directoryPi0276TeV    ->Get("graphInvCrossSectionPi0Comb2760GeVAStatErr");
      csGraphsSys[3]     = (TGraphAsymmErrors*)directoryPi0276TeV ->Get("graphInvCrossSectionPi0Comb2760GeVASysErr");
      csGraphs[2]     = (TGraphAsymmErrors*)directoryPi08TeV      ->Get("graphInvCrossSectionPi0Comb8TeVAStatErr");
      csGraphsSys[2]     = (TGraphAsymmErrors*)directoryPi08TeV   ->Get("graphInvCrossSectionPi0Comb8TeVASysErr");
      csGraphs[1]     = (TGraphAsymmErrors*)directoryPi07TeV      ->Get("graphInvCrossSectionPi0Comb7TeVAStatErr");
      csGraphsSys[1]     = (TGraphAsymmErrors*)directoryPi07TeV   ->Get("graphInvCrossSectionPi0Comb7TeVASysErr");
      csGraphs[0]     = (TGraphAsymmErrors*)directoryPi0900GeV    ->Get("graphInvCrossSectionPi0Comb900GeVAStatErr");
      csGraphsSys[0]     = (TGraphAsymmErrors*)directoryPi0900GeV ->Get("graphInvCrossSectionPi0Comb900GeVASysErr");
    
        } else {
          csGraphs[3]       = (TGraphAsymmErrors*)directoryPi0276TeV->Get("graphInvCrossSectionPi0PCM2760GeVStatErr");
          csGraphsSys[3]    = (TGraphAsymmErrors*)directoryPi0276TeV->Get("graphInvCrossSectionPi0PCM2760GeVSysErr");
          csGraphs[2]       = (TGraphAsymmErrors*)directoryPi08TeV  ->Get("graphInvCrossSectionPi0PCM8TeVStatErr");
          csGraphsSys[2]    = (TGraphAsymmErrors*)directoryPi08TeV  ->Get("graphInvCrossSectionPi0PCM8TeVSysErr");
          csGraphs[1]       = (TGraphAsymmErrors*)directoryPi07TeV  ->Get("graphInvCrossSectionPi0PCM7TeVStatErr");
          csGraphsSys[1]    = (TGraphAsymmErrors*)directoryPi07TeV  ->Get("graphInvCrossSectionPi0PCM7TeVSysErr");
          csGraphs[0]       = (TGraphAsymmErrors*)directoryPi0900GeV->Get("graphInvCrossSectionPi0PCM900GeVStatErr");
          csGraphsSys[0]    = (TGraphAsymmErrors*)directoryPi0900GeV->Get("graphInvCrossSectionPi0PCM900GeVSysErr");
          
        }
    Style_t markerstylesFULL[4]                         ={34,20,33,29};
    Size_t markersizeFULL[4]                        ={2.3,2.3,3.4,3.4};
    Color_t colorDataFULL[4]                        ={kRed+2,kBlue+2,kGreen+2,kMagenta+2};
  
    // systematics

        for (int j=0;j<csGraphsSys[2]->GetN();j++){
            csGraphsSys[2]->GetY()[j] *= 1000;
            csGraphsSys[2]->GetEYhigh()[j] *= 1000;
            csGraphsSys[2]->GetEYlow()[j] *= 1000;
        }
        DrawGammaSetMarkerTGraphAsym(csGraphsSys[2], markerstylesFULL[2], markersizeFULL[2], colorDataFULL[2], colorDataFULL[2], widthLinesBoxes, kTRUE);
        
        for (int j=0;j<csGraphsSys[1]->GetN();j++){
            csGraphsSys[1]->GetY()[j] *= 100;
            csGraphsSys[1]->GetEYhigh()[j] *= 100;
            csGraphsSys[1]->GetEYlow()[j] *= 100;
        }
        DrawGammaSetMarkerTGraphAsym(csGraphsSys[1], markerstylesFULL[1], markersizeFULL[1], colorDataFULL[1], colorDataFULL[1], widthLinesBoxes, kTRUE);
            for (int j=0;j<csGraphsSys[3]->GetN();j++){
                csGraphsSys[3]->GetY()[j] *= 10;
                csGraphsSys[3]->GetEYhigh()[j] *= 10;
                csGraphsSys[3]->GetEYlow()[j] *= 10;
            }
        DrawGammaSetMarkerTGraphAsym(csGraphsSys[3], markerstylesFULL[3], markersizeFULL[3], colorDataFULL[3], colorDataFULL[3], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(csGraphsSys[0], markerstylesFULL[0], markersizeFULL[0], colorDataFULL[0], colorDataFULL[0], widthLinesBoxes, kTRUE);
        //
        // statistics
        
        for (int j=0;j<csGraphs[2]->GetN();j++){
            csGraphs[2]->GetY()[j] *= 1000;
            csGraphs[2]->GetEYhigh()[j] *= 1000;
            csGraphs[2]->GetEYlow()[j] *= 1000;
        }
        DrawGammaSetMarkerTGraph(csGraphs[2], markerstylesFULL[2], 0.7*markersizeFULL[2], colorDataFULL[2] , colorDataFULL[2]);
        
        for (int j=0;j<csGraphs[1]->GetN();j++){
            csGraphs[1]->GetY()[j] *= 100;
            csGraphs[1]->GetEYhigh()[j] *= 100;
            csGraphs[1]->GetEYlow()[j] *= 100;
        }
        DrawGammaSetMarkerTGraph(csGraphs[1], markerstylesFULL[1], 0.7*markersizeFULL[1], colorDataFULL[1] , colorDataFULL[1]);
        for (int j=0;j<csGraphs[3]->GetN();j++){
            csGraphs[3]->GetY()[j] *= 10;
            csGraphs[3]->GetEYhigh()[j] *= 10;
            csGraphs[3]->GetEYlow()[j] *= 10;
        }
        DrawGammaSetMarkerTGraph(csGraphs[3], markerstylesFULL[3], 0.7*markersizeFULL[3], colorDataFULL[3] , colorDataFULL[3]);

        DrawGammaSetMarkerTGraph(csGraphs[0], markerstylesFULL[0], 0.7*markersizeFULL[0], colorDataFULL[0] , colorDataFULL[0]);





    Double_t textSizeLabelsPixel                     = 48;
    TCanvas* canvasRatioPP                  = new TCanvas("canvasRatioPP","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioPP,  0.12, 0.01, 0.01, 0.11);
    canvasRatioPP->cd();
    canvasRatioPP->SetLogx();

        Double_t textsizeLabelsPP                    = 0;
        Double_t textsizeFacPP                       = 0;
        if (canvasRatioPP->XtoPixel(canvasRatioPP->GetX2()) <canvasRatioPP->YtoPixel(canvasRatioPP->GetY1()) ){
            textsizeLabelsPP                = (Double_t)textSizeLabelsPixel/canvasRatioPP->XtoPixel(canvasRatioPP->GetX2()) ;
            textsizeFacPP                   = (Double_t)1./canvasRatioPP->XtoPixel(canvasRatioPP->GetX2()) ;
        } else {
        textsizeLabelsPP                    = (Double_t)textSizeLabelsPixel/canvasRatioPP->YtoPixel(canvasRatioPP->GetY1());
        textsizeFacPP                       = (Double_t)1./canvasRatioPP->YtoPixel(canvasRatioPP->GetY1());
        }
        cout << textsizeLabelsPP << endl;

    Double_t mesonMassExpectPi0                 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    Double_t mesonMassExpectEta                 = TDatabasePDG::Instance()->GetParticle(221)->Mass();
    Color_t  colorComb                          = kMagenta+2;
    Style_t  markerStyleComb                    = 20;
    Size_t   markersizeFULLComb                     = 2;
         Color_t  colorCGC                           = kCyan-8;
    Color_t  colorNLO                           = kAzure-4;
    Style_t  styleMarkerNLOMuHalf               = 24;
    Style_t  styleMarkerNLOMuOne                = 27;
    Style_t  styleMarkerNLOMuTwo                = 30;
    Style_t  styleLineNLOMuHalf                 = 8;
    Style_t  styleLineNLOMuOne                  = 7;
    Style_t  styleLineNLOMuTwo                  = 4;
    Style_t  styleLineNLOMuTwoBKK               = 3;
    Style_t  styleLineNLOMuTwoDSS               = 6;
    Size_t   sizeMarkerNLO                      = 1;
    Width_t  widthLineNLO                       = 2.;



    TString fileNameTheory                      = "ExternalInput/Theory/TheoryCompilationPP.root";
    TFile* fileTheoryCompilation                            = new TFile(fileNameTheory.Data());

        // 8TeV Pythia8 Monash2013:
        histoPythia8InvXSectionPi08TeV                       = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoPi08TeV");
        if(usecombinedspectra)
        histoPythia8InvXSectionPi08TeV->GetXaxis()->SetRangeUser(0.3,25);
        else
        histoPythia8InvXSectionPi08TeV->GetXaxis()->SetRangeUser(0.3,14);
        histoPythia8InvXSectionPi08TeV->Scale(1000);
        // 7TeV Pythia8 Monash2013:
        histoPythia8InvXSectionPi07TeV                       = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013Pi07TeV");
        if(usecombinedspectra)
        histoPythia8InvXSectionPi07TeV->GetXaxis()->SetRangeUser(0.3,25);
        else
        histoPythia8InvXSectionPi07TeV->GetXaxis()->SetRangeUser(0.3,20);
        histoPythia8InvXSectionPi07TeV->Scale(100);
        // 2.76TeV Pythia8 Monash2013:
        histoPythia8InvXSectionPi0276TeV                       = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013Pi02760GeV");
        if(usecombinedspectra)
        histoPythia8InvXSectionPi0276TeV->GetXaxis()->SetRangeUser(0.3,40);
        else
        histoPythia8InvXSectionPi0276TeV->GetXaxis()->SetRangeUser(0.3,8);
        histoPythia8InvXSectionPi0276TeV->Scale(10);
        // 0.9TeV Pythia8 Monash2013:
        histoPythia8InvXSectionPi0900GeV                       = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013Pi0900GeV");
        if(usecombinedspectra)
        histoPythia8InvXSectionPi0900GeV->GetXaxis()->SetRangeUser(0.4,7.);
        else
        histoPythia8InvXSectionPi0900GeV->GetXaxis()->SetRangeUser(0.4,4);

        TGraphErrors* graphPythia8InvXSection8TeV               = new TGraphErrors((TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoPi08TeV"));
        for (int j=0;j<graphPythia8InvXSection8TeV->GetN();j++){
            graphPythia8InvXSection8TeV->GetY()[j] *= 1000;
            graphPythia8InvXSection8TeV->GetEY()[j] *= 1000;
        }
        if(usecombinedspectra){
          while(graphPythia8InvXSection8TeV->GetX()[0] < 0.3) graphPythia8InvXSection8TeV->RemovePoint(0);
          while(graphPythia8InvXSection8TeV->GetX()[graphPythia8InvXSection8TeV->GetN()-1] > 35.) graphPythia8InvXSection8TeV->RemovePoint(graphPythia8InvXSection8TeV->GetN()-1);
        } else {
          while(graphPythia8InvXSection8TeV->GetX()[0] < 0.3) graphPythia8InvXSection8TeV->RemovePoint(0);
          while(graphPythia8InvXSection8TeV->GetX()[graphPythia8InvXSection8TeV->GetN()-1] > 14.) graphPythia8InvXSection8TeV->RemovePoint(graphPythia8InvXSection8TeV->GetN()-1);
        }
        TGraphErrors* graphPythia8InvXSection7TeV               = new TGraphErrors((TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013Pi07TeV"));
        for (int j=0;j<graphPythia8InvXSection7TeV->GetN();j++){
            graphPythia8InvXSection7TeV->GetY()[j] *= 100;
            graphPythia8InvXSection7TeV->GetEY()[j] *= 100;
        }
        if(usecombinedspectra){
          while(graphPythia8InvXSection7TeV->GetX()[0] < 0.3) graphPythia8InvXSection7TeV->RemovePoint(0);
          while(graphPythia8InvXSection7TeV->GetX()[graphPythia8InvXSection7TeV->GetN()-1] > 25.) graphPythia8InvXSection7TeV->RemovePoint(graphPythia8InvXSection7TeV->GetN()-1);
        } else {
          while(graphPythia8InvXSection7TeV->GetX()[0] < 0.3) graphPythia8InvXSection7TeV->RemovePoint(0);
          while(graphPythia8InvXSection7TeV->GetX()[graphPythia8InvXSection7TeV->GetN()-1] > 20.) graphPythia8InvXSection7TeV->RemovePoint(graphPythia8InvXSection7TeV->GetN()-1);
        }
        TGraphErrors* graphPythia8InvXSection276TeV               = new TGraphErrors((TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013Pi02760GeV"));
        for (int j=0;j<graphPythia8InvXSection276TeV->GetN();j++){
            graphPythia8InvXSection276TeV->GetY()[j] *= 10;
            graphPythia8InvXSection276TeV->GetEY()[j] *= 10;
        }
        if(usecombinedspectra){
          while(graphPythia8InvXSection276TeV->GetX()[0] < 0.3) graphPythia8InvXSection276TeV->RemovePoint(0);
          while(graphPythia8InvXSection276TeV->GetX()[graphPythia8InvXSection276TeV->GetN()-1] > 40.) graphPythia8InvXSection276TeV->RemovePoint(graphPythia8InvXSection276TeV->GetN()-1);
        } else {
          while(graphPythia8InvXSection276TeV->GetX()[0] < 0.3) graphPythia8InvXSection276TeV->RemovePoint(0);
          while(graphPythia8InvXSection276TeV->GetX()[graphPythia8InvXSection276TeV->GetN()-1] > 9.) graphPythia8InvXSection276TeV->RemovePoint(graphPythia8InvXSection276TeV->GetN()-1);
        }
        TGraphErrors* graphPythia8InvXSection900GeV               = new TGraphErrors((TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013Pi0900GeV"));
        if(usecombinedspectra){
          while(graphPythia8InvXSection900GeV->GetX()[0] < 0.4) graphPythia8InvXSection900GeV->RemovePoint(0);
          while(graphPythia8InvXSection900GeV->GetX()[graphPythia8InvXSection900GeV->GetN()-1] > 7.) graphPythia8InvXSection900GeV->RemovePoint(graphPythia8InvXSection900GeV->GetN()-1);
        } else {
          while(graphPythia8InvXSection900GeV->GetX()[0] < 0.4) graphPythia8InvXSection900GeV->RemovePoint(0);
          while(graphPythia8InvXSection900GeV->GetX()[graphPythia8InvXSection900GeV->GetN()-1] > 4.) graphPythia8InvXSection900GeV->RemovePoint(graphPythia8InvXSection900GeV->GetN()-1);
        }

        // *******************************************************************************************************
        // NLO calc
        TGraphAsymmErrors* graphPi0DSS148TeV                    = (TGraphAsymmErrors*) fileTheoryCompilation->Get("graphNLOCalcDSS14InvSecPi08000GeV");
        for (int j=0;j<graphPi0DSS148TeV->GetN();j++){
            graphPi0DSS148TeV->GetY()[j] *= 1000;
            graphPi0DSS148TeV->GetEYhigh()[j] *= 1000;
            graphPi0DSS148TeV->GetEYlow()[j] *= 1000;
        }
        if(usecombinedspectra){
          while (graphPi0DSS148TeV->GetX()[graphPi0DSS148TeV->GetN()-1] > 38. ) graphPi0DSS148TeV->RemovePoint(graphPi0DSS148TeV->GetN()-1);
        } else {
          while (graphPi0DSS148TeV->GetX()[graphPi0DSS148TeV->GetN()-1] > 14. ) graphPi0DSS148TeV->RemovePoint(graphPi0DSS148TeV->GetN()-1);
        }
        TGraphAsymmErrors* graphPi0DSS147TeV                    = (TGraphAsymmErrors*) fileTheoryCompilation->Get("graphNLOCalcDSS14InvCrossSec7000GeV");
        for (int j=0;j<graphPi0DSS147TeV->GetN();j++){
            graphPi0DSS147TeV->GetY()[j] *= 100;
            graphPi0DSS147TeV->GetEYhigh()[j] *= 100;
            graphPi0DSS147TeV->GetEYlow()[j] *= 100;
        }
        graphPi0DSS147TeV->RemovePoint(0);
        if(usecombinedspectra){
          while (graphPi0DSS147TeV->GetX()[graphPi0DSS147TeV->GetN()-1] > 30. ) graphPi0DSS147TeV->RemovePoint(graphPi0DSS147TeV->GetN()-1);
        } else {
          while (graphPi0DSS147TeV->GetX()[graphPi0DSS147TeV->GetN()-1] > 20. ) graphPi0DSS147TeV->RemovePoint(graphPi0DSS147TeV->GetN()-1);
        }
        // TGraphAsymmErrors* graphPi0DSS14276TeV                    = (TGraphAsymmErrors*) fileTheoryCompilation->Get("graphNLOCalcDSS07InvSecPi02760GeV");
        TGraphAsymmErrors* graphPi0DSS14276TeV                    = (TGraphAsymmErrors*) fileTheoryCompilation->Get("graphNLOCalcDSS14InvCrossSec2760GeV");
        for (int j=0;j<graphPi0DSS14276TeV->GetN();j++){
            graphPi0DSS14276TeV->GetY()[j] *= 10;
            graphPi0DSS14276TeV->GetEYhigh()[j] *= 10;
            graphPi0DSS14276TeV->GetEYlow()[j] *= 10;
        }
        graphPi0DSS14276TeV->RemovePoint(0);
        if(usecombinedspectra){
          while (graphPi0DSS14276TeV->GetX()[graphPi0DSS14276TeV->GetN()-1] > 30. ) graphPi0DSS14276TeV->RemovePoint(graphPi0DSS14276TeV->GetN()-1);
        } else {
          while (graphPi0DSS14276TeV->GetX()[graphPi0DSS14276TeV->GetN()-1] > 9. ) graphPi0DSS14276TeV->RemovePoint(graphPi0DSS14276TeV->GetN()-1);
        }


        // *******************************************************************************************************

        TGraph* graphNLOCalcPi0MuHalf                       = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuHalf8000GeV");
        TGraph* graphNLOCalcPi0MuOne                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuOne8000GeV");
        TGraph* graphNLOCalcPi0MuTwo                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuTwo8000GeV");



    Double_t arrayBoundariesX1_XSec[2];
    Double_t arrayBoundariesY1_XSec[10];
    Double_t relativeMarginsXXSec[3];
    Double_t relativeMarginsYXSec[3];
    textSizeLabelsPixel = 48;
    ReturnCorrectValuesForCanvasScaling(1250,2500, 1, 9,0.135, 0.005, 0.003,0.04,arrayBoundariesX1_XSec,arrayBoundariesY1_XSec,relativeMarginsXXSec,relativeMarginsYXSec);

    TCanvas* canvasInvSectionPaper      = new TCanvas("canvasInvSectionPaper","",0,0,1250,2500);  // gives the page size
    DrawGammaCanvasSettings( canvasInvSectionPaper,  0.13, 0.02, 0.03, 0.06);

    TPad* padInvSectionSpec             = new TPad("padInvSectionSpec", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[5], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[0],-1, -1, -2);
    DrawGammaPadSettings( padInvSectionSpec, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[0], relativeMarginsYXSec[1]);
    padInvSectionSpec->Draw();
    Double_t marginXSec                 = relativeMarginsXXSec[0]*1250;
    Double_t textsizeLabelsXSecUp       = 0;
    Double_t textsizeFacXSecUp          = 0;
    if (padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) < padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1())){
        textsizeLabelsXSecUp            = (Double_t)textSizeLabelsPixel/padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) ;
        textsizeFacXSecUp               = (Double_t)1./padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) ;
    } else {
        textsizeLabelsXSecUp            = (Double_t)textSizeLabelsPixel/padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1());
        textsizeFacXSecUp               = (Double_t)1./padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1());
    }

    TPad* padInvSectionNLORatio         = new TPad("padInvSectionNLORatio", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[6], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[5],-1, -1, -2);
    DrawGammaPadSettings( padInvSectionNLORatio, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[1]);
    padInvSectionNLORatio->Draw();
    Double_t textsizeLabelsXSecMiddle   = 0;
    Double_t textsizeFacXSecMiddle      = 0;
    if (padInvSectionNLORatio->XtoPixel(padInvSectionNLORatio->GetX2()) < padInvSectionNLORatio->YtoPixel(padInvSectionNLORatio->GetY1())){
        textsizeLabelsXSecMiddle        = (Double_t)textSizeLabelsPixel/padInvSectionNLORatio->XtoPixel(padInvSectionNLORatio->GetX2()) ;
        textsizeFacXSecMiddle           = (Double_t)1./padInvSectionNLORatio->XtoPixel(padInvSectionNLORatio->GetX2()) ;
    } else {
        textsizeLabelsXSecMiddle        = (Double_t)textSizeLabelsPixel/padInvSectionNLORatio->YtoPixel(padInvSectionNLORatio->GetY1());
        textsizeFacXSecMiddle           = (Double_t)1./padInvSectionNLORatio->YtoPixel(padInvSectionNLORatio->GetY1());
    }

    TPad* padInvSectionPythiaRatio      = new TPad("padInvSectionPythiaRatio", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[7], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[6],-1, -1, -2);
    DrawGammaPadSettings( padInvSectionPythiaRatio, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[1]);
    padInvSectionPythiaRatio->Draw();
    Double_t textsizeLabelsXSecDown     = 0;
    Double_t textsizeFacXSecDown        = 0;
    if (padInvSectionPythiaRatio->XtoPixel(padInvSectionPythiaRatio->GetX2()) < padInvSectionPythiaRatio->YtoPixel(padInvSectionPythiaRatio->GetY1())){
        textsizeLabelsXSecDown          = (Double_t)textSizeLabelsPixel/padInvSectionPythiaRatio->XtoPixel(padInvSectionPythiaRatio->GetX2()) ;
        textsizeFacXSecDown             = (Double_t)1./padInvSectionPythiaRatio->XtoPixel(padInvSectionPythiaRatio->GetX2()) ;
    } else {
        textsizeLabelsXSecDown          = (Double_t)textSizeLabelsPixel/padInvSectionPythiaRatio->YtoPixel(padInvSectionPythiaRatio->GetY1());
        textsizeFacXSecDown             = (Double_t)1./padInvSectionPythiaRatio->YtoPixel(padInvSectionPythiaRatio->GetY1());
    }
    TPad* padInvSection2760GeVRatio      = new TPad("padInvSection2760GeVRatio", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[8], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[7],-1, -1, -2);
    DrawGammaPadSettings( padInvSection2760GeVRatio, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[1]);
    padInvSection2760GeVRatio->Draw();
    Double_t textsizeLabelsXSecDown26     = 0;
    Double_t textsizeFacXSecDown26        = 0;
    if (padInvSection2760GeVRatio->XtoPixel(padInvSection2760GeVRatio->GetX2()) < padInvSection2760GeVRatio->YtoPixel(padInvSection2760GeVRatio->GetY1())){
        textsizeLabelsXSecDown26          = (Double_t)textSizeLabelsPixel/padInvSection2760GeVRatio->XtoPixel(padInvSection2760GeVRatio->GetX2()) ;
        textsizeFacXSecDown26             = (Double_t)1./padInvSection2760GeVRatio->XtoPixel(padInvSection2760GeVRatio->GetX2()) ;
    } else {
        textsizeLabelsXSecDown26          = (Double_t)textSizeLabelsPixel/padInvSection2760GeVRatio->YtoPixel(padInvSection2760GeVRatio->GetY1());
        textsizeFacXSecDown26             = (Double_t)1./padInvSection2760GeVRatio->YtoPixel(padInvSection2760GeVRatio->GetY1());
    }
    TPad* padInvSection900GeVRatio      = new TPad("padInvSection900GeVRatio", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[9], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[8],-1, -1, -2);
    DrawGammaPadSettings( padInvSection900GeVRatio, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[2]);
    padInvSection900GeVRatio->Draw();
    Double_t textsizeLabelsXSecDown2     = 0;
    Double_t textsizeFacXSecDown2        = 0;
    if (padInvSection900GeVRatio->XtoPixel(padInvSection900GeVRatio->GetX2()) < padInvSection900GeVRatio->YtoPixel(padInvSection900GeVRatio->GetY1())){
        textsizeLabelsXSecDown2          = (Double_t)textSizeLabelsPixel/padInvSection900GeVRatio->XtoPixel(padInvSection900GeVRatio->GetX2()) ;
        textsizeFacXSecDown2             = (Double_t)1./padInvSection900GeVRatio->XtoPixel(padInvSection900GeVRatio->GetX2()) ;
    } else {
        textsizeLabelsXSecDown2          = (Double_t)textSizeLabelsPixel/padInvSection900GeVRatio->YtoPixel(padInvSection900GeVRatio->GetY1());
        textsizeFacXSecDown2             = (Double_t)1./padInvSection900GeVRatio->YtoPixel(padInvSection900GeVRatio->GetY1());
    }
    Double_t maxX = 22;
    if(usecombinedspectra)
      maxX = 49;

    padInvSectionSpec->cd();
    padInvSectionSpec->SetLogy(1);
    padInvSectionSpec->SetLogx(1);

    TH2F * histo2DXSectionPi0;
    // histo2DXSectionPi0          = new TH2F("histo2DXSectionPi0","histo2DXSectionPi0",11000,0.23,50.,1000,6,9e11);
    if(usecombinedspectra)
    histo2DXSectionPi0          = new TH2F("histo2DXSectionPi0","histo2DXSectionPi0",11000,0.23,maxX,1000,4.1,9e14);
    else
    histo2DXSectionPi0          = new TH2F("histo2DXSectionPi0","histo2DXSectionPi0",11000,0.23,maxX,1000,4.1e2,9e14);
    SetStyleHistoTH2ForGraphs(histo2DXSectionPi0, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",0.035,0.04, 0.035,0.04, 0.9,1.45);
    histo2DXSectionPi0->GetXaxis()->SetMoreLogLabels();
    histo2DXSectionPi0->GetXaxis()->SetNoExponent(kTRUE);

        SetStyleHistoTH2ForGraphs(histo2DXSectionPi0, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",
                                0.85*textsizeLabelsXSecUp,textsizeLabelsXSecUp, 0.85*textsizeLabelsXSecUp, textsizeLabelsXSecUp, 1,0.2/(textsizeFacXSecUp*marginXSec));
        histo2DXSectionPi0->GetXaxis()->SetMoreLogLabels();
        histo2DXSectionPi0->GetXaxis()->SetLabelOffset(+0.01);
        histo2DXSectionPi0->Draw();

        Color_t pythia8color = kRed+2;
        if(plotTheorycurves){
          graphPi0DSS148TeV->SetLineWidth(widthCommonFit);
          graphPi0DSS148TeV->SetLineColor(colorNLO);
          graphPi0DSS148TeV->SetLineStyle(1);
          graphPi0DSS148TeV->SetFillStyle(1001);
          graphPi0DSS148TeV->SetFillColor(colorNLO);
          graphPi0DSS148TeV->Draw("same,e3");
          
          graphPi0DSS147TeV->SetLineWidth(widthCommonFit);
          graphPi0DSS147TeV->SetLineColor(colorNLO);
          graphPi0DSS147TeV->SetLineStyle(1);
          graphPi0DSS147TeV->SetFillStyle(1001);
          graphPi0DSS147TeV->SetFillColor(colorNLO);
          graphPi0DSS147TeV->Draw("same,e3");
          
          graphPi0DSS14276TeV->SetLineWidth(widthCommonFit);
          graphPi0DSS14276TeV->SetLineColor(colorNLO);
          graphPi0DSS14276TeV->SetLineStyle(1);
          graphPi0DSS14276TeV->SetFillStyle(1001);
          graphPi0DSS14276TeV->SetFillColor(colorNLO);
          graphPi0DSS14276TeV->Draw("same,e3");


          DrawGammaSetMarkerTGraphErr(graphPythia8InvXSection8TeV, 0, 0, pythia8color , pythia8color, widthLinesBoxes, kTRUE, pythia8color);
          graphPythia8InvXSection8TeV->Draw("3,same");
          DrawGammaSetMarkerTGraphErr(graphPythia8InvXSection7TeV, 0, 0, pythia8color , pythia8color, widthLinesBoxes, kTRUE, pythia8color);
          graphPythia8InvXSection7TeV->Draw("3,same");
          DrawGammaSetMarkerTGraphErr(graphPythia8InvXSection276TeV, 0, 0, pythia8color , pythia8color, widthLinesBoxes, kTRUE, pythia8color);
          graphPythia8InvXSection276TeV->Draw("3,same");
          DrawGammaSetMarkerTGraphErr(graphPythia8InvXSection900GeV, 0, 0, pythia8color , pythia8color, widthLinesBoxes, kTRUE, pythia8color);
          graphPythia8InvXSection900GeV->Draw("3,same");
          DrawGammaSetMarker(histoPythia8InvXSectionPi08TeV, 24, 1.5, pythia8color , pythia8color);
          histoPythia8InvXSectionPi08TeV->SetLineWidth(widthCommonFit);
          histoPythia8InvXSectionPi08TeV->Draw("same,hist,l");
          DrawGammaSetMarker(histoPythia8InvXSectionPi07TeV, 24, 1.5, pythia8color , pythia8color);
          histoPythia8InvXSectionPi07TeV->SetLineWidth(widthCommonFit);
          histoPythia8InvXSectionPi07TeV->Draw("same,hist,l");
          DrawGammaSetMarker(histoPythia8InvXSectionPi0276TeV, 24, 1.5, pythia8color , pythia8color);
          histoPythia8InvXSectionPi0276TeV->SetLineWidth(widthCommonFit);
          histoPythia8InvXSectionPi0276TeV->Draw("same,hist,l");
          DrawGammaSetMarker(histoPythia8InvXSectionPi0900GeV, 24, 1.5, pythia8color , pythia8color);
          histoPythia8InvXSectionPi0900GeV->SetLineWidth(widthCommonFit);
          histoPythia8InvXSectionPi0900GeV->Draw("same,hist,l");
        }

         Double_t paramTCMPi0New8TeV[5]  = { csGraphs[2]->GetY()[1],0.1,  csGraphs[2]->GetY()[4],0.6,3.0};
         Double_t fitmaxx8TeV = 12;
         if(usecombinedspectra)
          fitmaxx8TeV = 35;
        TF1* fitTCMInvXSectionPi08TeV   = FitObject("tcm","fitTCMInvCrossSectionPi08TeV","Pi0",csGraphs[2],0.3,fitmaxx8TeV ,paramTCMPi0New8TeV,"QNRMEX0+","", kFALSE);

        TF1* fitTCMInvXSectionPi0Plot8TeV = new TF1("twoCompModel_plotting8TeV",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mesonMassExpectPi0,mesonMassExpectPi0,mesonMassExpectPi0));
        fitTCMInvXSectionPi0Plot8TeV->SetRange(0.3,fitmaxx8TeV);
        fitTCMInvXSectionPi0Plot8TeV->SetParameters(fitTCMInvXSectionPi08TeV->GetParameters());
        fitTCMInvXSectionPi0Plot8TeV->SetParErrors(fitTCMInvXSectionPi08TeV->GetParErrors());
        cout << WriteParameterToFile(fitTCMInvXSectionPi08TeV) << endl;
        DrawGammaSetMarkerTF1( fitTCMInvXSectionPi0Plot8TeV, 7, 2, kGray+2);
        fitTCMInvXSectionPi0Plot8TeV->Draw("same");

         Double_t paramTCMPi0New7TeV[5]  = { csGraphs[1]->GetY()[1],0.1,  csGraphs[1]->GetY()[4],0.6,3.0};
         Double_t fitmaxx7TeV = 16;
         if(usecombinedspectra)
          fitmaxx7TeV = 25;
        TF1* fitTCMInvXSectionPi07TeV   = FitObject("tcm","fitTCMInvCrossSectionPi07TeV","Pi0",csGraphs[1],0.3,fitmaxx7TeV ,paramTCMPi0New7TeV,"QNRMEX0+","", kFALSE);

        TF1* fitTCMInvXSectionPi0Plot7TeV = new TF1("twoCompModel_plotting7TeV",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mesonMassExpectPi0,mesonMassExpectPi0,mesonMassExpectPi0));
        fitTCMInvXSectionPi0Plot7TeV->SetRange(0.3,fitmaxx7TeV);
        fitTCMInvXSectionPi0Plot7TeV->SetParameters(fitTCMInvXSectionPi07TeV->GetParameters());
        fitTCMInvXSectionPi0Plot7TeV->SetParErrors(fitTCMInvXSectionPi07TeV->GetParErrors());
        cout << WriteParameterToFile(fitTCMInvXSectionPi07TeV) << endl;
        DrawGammaSetMarkerTF1( fitTCMInvXSectionPi0Plot7TeV, 7, 2, kGray+2);
        fitTCMInvXSectionPi0Plot7TeV->Draw("same");
        
         Double_t paramTCMPi0New276TeV[5]  = { csGraphs[3]->GetY()[1],0.1,csGraphs[3]->GetY()[4],0.6,3.0};
         Double_t fitmaxx276TeV = 8;
         if(usecombinedspectra)
          fitmaxx276TeV = 25;
        TF1* fitTCMInvXSectionPi0276TeV   = FitObject("tcm","fitTCMInvCrossSectionPi0276TeV","Pi0",csGraphs[3],0.3,fitmaxx276TeV ,paramTCMPi0New276TeV,"QNRMEX0+","", kFALSE);

        TF1* fitTCMInvXSectionPi0Plot276TeV = new TF1("twoCompModel_plotting7TeV",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mesonMassExpectPi0,mesonMassExpectPi0,mesonMassExpectPi0));
        fitTCMInvXSectionPi0Plot276TeV->SetRange(0.3,fitmaxx276TeV);
        fitTCMInvXSectionPi0Plot276TeV->SetParameters(fitTCMInvXSectionPi0276TeV->GetParameters());
        fitTCMInvXSectionPi0Plot276TeV->SetParErrors(fitTCMInvXSectionPi0276TeV->GetParErrors());
        cout << WriteParameterToFile(fitTCMInvXSectionPi0276TeV) << endl;
        DrawGammaSetMarkerTF1( fitTCMInvXSectionPi0Plot276TeV, 7, 2, kGray+2);
        fitTCMInvXSectionPi0Plot276TeV->Draw("same");

         Double_t paramTCMPi0New900GeV[5]  = { csGraphs[0]->GetY()[1],0.1,  csGraphs[0]->GetY()[4],0.6,3.0};
         Double_t fitmaxx900GeV = 4;
         if(usecombinedspectra)
          fitmaxx900GeV = 7;
        TF1* fitTCMInvXSectionPi0900GeV   = FitObject("tcm","fitTCMInvCrossSectionPi0900GeV","Pi0",csGraphs[0],0.4,fitmaxx900GeV ,paramTCMPi0New900GeV,"QNRMEX0+","", kFALSE);
        TF1* fitTCMInvXSectionPi0Plot900GeV = new TF1("twoCompModel_plotting900GeV",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mesonMassExpectPi0,mesonMassExpectPi0,mesonMassExpectPi0));
        fitTCMInvXSectionPi0Plot900GeV->SetRange(0.4,fitmaxx900GeV);
        fitTCMInvXSectionPi0Plot900GeV->SetParameters(fitTCMInvXSectionPi0900GeV->GetParameters());
        fitTCMInvXSectionPi0Plot900GeV->SetParErrors(fitTCMInvXSectionPi0900GeV->GetParErrors());
        cout << WriteParameterToFile(fitTCMInvXSectionPi0900GeV) << endl;
        DrawGammaSetMarkerTF1( fitTCMInvXSectionPi0Plot900GeV, 7, 2, kGray+2);
        fitTCMInvXSectionPi0Plot900GeV->Draw("same");


        // Tsallis fit
        Double_t paramGraphPi08[3]                       = {5e11, 6., 0.13};
        TF1* fitInvXSectionPi08TeV                      = FitObject("l","fitInvCrossSectionPi08TeV","Pi0",csGraphs[2],0.3,fitmaxx8TeV,paramGraphPi08,"QNRMEX0+");
        DrawGammaSetMarkerTF1( fitInvXSectionPi08TeV, 3, 2, kGray+1);
        fitInvXSectionPi08TeV->Draw("same");
        cout << WriteParameterToFile(fitInvXSectionPi08TeV) << endl;
        Double_t paramGraphPi07[3]                                   = {5e11, 6., 0.13};
        TF1* fitInvXSectionPi07TeV                      = FitObject("l","fitInvCrossSectionPi07TeV","Pi0",csGraphs[1],0.3,fitmaxx7TeV,paramGraphPi07,"QNRMEX0+");
        DrawGammaSetMarkerTF1( fitInvXSectionPi07TeV, 3, 2, kGray+1);
        fitInvXSectionPi07TeV->Draw("same");
        Double_t paramGraphPi0276[3]                                   = {5e11, 6., 0.13};
        TF1* fitInvXSectionPi0276TeV                      = FitObject("l","fitInvCrossSectionPi0276TeV","Pi0",csGraphs[3],0.3,fitmaxx276TeV,paramGraphPi0276,"QNRMEX0+");
        DrawGammaSetMarkerTF1( fitInvXSectionPi0276TeV, 3, 2, kGray+1);
        fitInvXSectionPi0276TeV->Draw("same");
        cout << WriteParameterToFile(fitInvXSectionPi07TeV) << endl;
        Double_t paramGraphPi0900[3]                                   = {5e11, 6., 0.13};
        TF1* fitInvXSectionPi0900GeV                    = FitObject("l","fitInvCrossSectionPi0900GeV","Pi0",csGraphs[0],0.4,fitmaxx900GeV,paramGraphPi0900,"QNRMEX0+");
        DrawGammaSetMarkerTF1( fitInvXSectionPi0900GeV, 3, 2, kGray+1);
        fitInvXSectionPi0900GeV->Draw("same");
        cout << WriteParameterToFile(fitInvXSectionPi0900GeV) << endl;

        csGraphsSys[2]     ->Draw("E2same");
        csGraphsSys[1]     ->Draw("E2same");
        csGraphsSys[3]     ->Draw("E2same");
        csGraphsSys[0]     ->Draw("E2same");
        csGraphs[2]->Draw("p,same,z");
        csGraphs[1]->Draw("p,same,z");
        csGraphs[3]->Draw("p,same,z");
        csGraphs[0]->Draw("p,same,z");



        Double_t rightalignDouble = 0.93;
        // drawLatexAdd("ALICE",rightalignDouble,0.86,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE);
        drawLatexAdd("ALICE",rightalignDouble,0.91,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE);
        // drawLatexAdd("proton-proton",rightalignDouble,0.86,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE);
        // drawLatexAdd("#gamma's rec. with PCM",rightalignDouble,0.86,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE);
        // drawLatexAdd("#gamma's rec. with PCM",rightalignDouble,0.91,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE);
        // drawLatexAdd("MB triggered",rightalignDouble,0.81,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE);
        // drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",rightalignDouble,0.76,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE);
        drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",rightalignDouble,0.86,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE);

        TBox* boxErrorSigmaRatio = CreateBoxConv(kGray+2, 0.3, 1.-(0.0196 ), 0.35, 1.+(0.0196));
        boxErrorSigmaRatio->SetLineWidth(8);


        TLegend* legendXsectionPaper;
        if(plotTheorycurves)
        legendXsectionPaper    = GetAndSetLegend2(0.17, 0.03, 0.5, 0.12+0.048*6, textSizeLabelsPixel);
        else
        legendXsectionPaper    = GetAndSetLegend2(0.17, 0.03, 0.5, 0.12+0.048*4, textSizeLabelsPixel);
        legendXsectionPaper->SetNColumns(1);
        legendXsectionPaper->SetMargin(0.2);
        legendXsectionPaper->AddEntry(csGraphsSys[2],Form("%s (x 10^{3})",nameMeasGlobal[2].Data()),"pf");
        legendXsectionPaper->AddEntry(csGraphsSys[1],Form("%s (x 10^{2})",nameMeasGlobal[1].Data()),"pf");
        legendXsectionPaper->AddEntry(csGraphsSys[3],Form("%s (x 10)",nameMeasGlobal[3].Data()),"pf");
        legendXsectionPaper->AddEntry(csGraphsSys[0],nameMeasGlobal[0].Data(),"pf");
        legendXsectionPaper->AddEntry(fitTCMInvXSectionPi0Plot8TeV,"TCM fit","l");
        legendXsectionPaper->AddEntry(fitInvXSectionPi08TeV,"Tsallis fit","l");
        if(plotTheorycurves){
          legendXsectionPaper->AddEntry(histoPythia8InvXSectionPi08TeV,"PYTHIA 8.2, Monash 2013","l");
          legendXsectionPaper->AddEntry(graphPi0DSS148TeV,  "NLO, PDF:MSTW08 - FF:DSS14", "f");
        }
        legendXsectionPaper->Draw();

        // drawLatexAdd("x10^{2}",0.24,0.925,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE,colorDataFULL[2]);
        // drawLatexAdd("x 10^{2}",0.45,0.80,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE,colorDataFULL[2]);
        // drawLatexAdd("x 10^{2}",0.25,0.93,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE,colorDataFULL[2]);
        // drawLatexAdd("x10^{1}",0.24,0.81,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE,colorDataFULL[1]);
        //
        // TLatex *labelRatioTheoryPP_Paper   = new TLatex(0.24,0.055,"0.5#it{p}_{T} < #mu < 2#it{p}_{T}");
        // SetStyleTLatex( labelRatioTheoryPP_Paper, 0.8*textsizeLabelsPP,4);
        // labelRatioTheoryPP_Paper->Draw();

    padInvSectionNLORatio->cd();
    padInvSectionNLORatio->SetLogx(1);
        TH2F * ratio8TeVdummy               = new TH2F("ratio8TeVdummy","ratio8TeVdummy",1000,0.23,maxX,1000,0.6,1.95);
        SetStyleHistoTH2ForGraphs(ratio8TeVdummy, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{TCM fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle,
                                  0.85*textsizeLabelsXSecMiddle,0.9*textsizeLabelsXSecMiddle, 1,0.23/(textsizeFacXSecMiddle*marginXSec), 510, 505);
        ratio8TeVdummy->GetYaxis()->SetMoreLogLabels(kTRUE);
        ratio8TeVdummy->GetYaxis()->SetNdivisions(505);
        ratio8TeVdummy->GetYaxis()->SetNoExponent(kTRUE);
        ratio8TeVdummy->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio8TeVdummy->GetXaxis()->SetNoExponent(kTRUE);
        ratio8TeVdummy->GetXaxis()->SetLabelFont(42);
        ratio8TeVdummy->GetYaxis()->SetLabelFont(42);
        ratio8TeVdummy->GetYaxis()->CenterTitle(kTRUE);
        ratio8TeVdummy->GetYaxis()->SetLabelOffset(+0.01);
        ratio8TeVdummy->GetXaxis()->SetTickLength(0.07);
        ratio8TeVdummy->DrawCopy();


        TGraphAsymmErrors* graphRatioPi0DSS148TeV            = (TGraphAsymmErrors*)graphPi0DSS148TeV->Clone();
        graphRatioPi0DSS148TeV                               = CalculateGraphErrRatioToFit (graphRatioPi0DSS148TeV, fitTCMInvXSectionPi0Plot8TeV);

        graphRatioPi0DSS148TeV->SetLineWidth(widthCommonFit);
        graphRatioPi0DSS148TeV->SetLineColor(colorNLO);
        graphRatioPi0DSS148TeV->SetLineStyle(1);
        graphRatioPi0DSS148TeV->SetFillStyle(1001);
        graphRatioPi0DSS148TeV->SetFillColor(colorNLO);
        if(plotTheorycurves)
          graphRatioPi0DSS148TeV->Draw("3,same");

        TH1D* histoRatioPythia8ToFit8TeV                     = (TH1D*) histoPythia8InvXSectionPi08TeV->Clone();
        histoRatioPythia8ToFit8TeV                           = CalculateHistoRatioToFit (histoRatioPythia8ToFit8TeV, fitTCMInvXSectionPi0Plot8TeV);
        DrawGammaSetMarker(histoRatioPythia8ToFit8TeV, 24, 1.5, pythia8color , pythia8color);
        histoRatioPythia8ToFit8TeV->SetLineWidth(widthCommonFit);
        if(plotTheorycurves)
          histoRatioPythia8ToFit8TeV->Draw("same,hist,l");


        TGraphErrors* graphRatioPythia8ToFit8TeV             = (TGraphErrors*) graphPythia8InvXSection8TeV->Clone();
        graphRatioPythia8ToFit8TeV                           = CalculateGraphErrRatioToFit (graphRatioPythia8ToFit8TeV, fitTCMInvXSectionPi0Plot8TeV);
        DrawGammaSetMarkerTGraphErr(graphRatioPythia8ToFit8TeV, 0, 0, pythia8color , pythia8color, widthLinesBoxes, kTRUE, pythia8color);
        if(plotTheorycurves)
          graphRatioPythia8ToFit8TeV->Draw("3,same");


        TGraphAsymmErrors* graphRatioPi0CombCombFitStatA    = (TGraphAsymmErrors*)csGraphs[2]->Clone();
        graphRatioPi0CombCombFitStatA                       = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitStatA, fitTCMInvXSectionPi0Plot8TeV);
        TGraphAsymmErrors* graphRatioPi0CombCombFitStatA_WOXErr = (TGraphAsymmErrors*) graphRatioPi0CombCombFitStatA->Clone("graphRatioPi0CombCombFitStatA_WOXErr");
        ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombCombFitStatA_WOXErr);

        TGraphAsymmErrors* graphRatioPi0CombCombFitSysA     = (TGraphAsymmErrors*)csGraphsSys[2]->Clone();
        graphRatioPi0CombCombFitSysA                        = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitSysA, fitTCMInvXSectionPi0Plot8TeV);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA,markerstylesFULL[2], markersizeFULL[2], colorDataFULL[2], colorDataFULL[2], widthLinesBoxes, kTRUE, 0);
        graphRatioPi0CombCombFitSysA->SetLineWidth(0);
        graphRatioPi0CombCombFitSysA->Draw("2,same");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA_WOXErr, markerstylesFULL[2], markersizeFULL[2], colorDataFULL[2], colorDataFULL[2], widthLinesBoxes, kFALSE);
        graphRatioPi0CombCombFitStatA_WOXErr->SetLineWidth(widthLinesBoxes);
        graphRatioPi0CombCombFitStatA_WOXErr->Draw("p,same");

        // boxErrorSigmaRatio->Draw();
        DrawGammaLines(0.23, maxX , 1., 1.,0.5, kGray+2);

    padInvSectionPythiaRatio->cd();
    padInvSectionPythiaRatio->SetLogx(1);
        TH2F * ratio7TeVdummy               = new TH2F("ratio7TeVdummy","ratio7TeVdummy",1000,0.23,maxX,1000,0.6,1.95);
        SetStyleHistoTH2ForGraphs(ratio7TeVdummy, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{TCM fit}", 0.85*textsizeLabelsXSecDown, 0.8*textsizeLabelsXSecDown,
                                  0.85*textsizeLabelsXSecDown,0.9*textsizeLabelsXSecDown, 1,0.23/(textsizeFacXSecMiddle*marginXSec), 510, 505);
        ratio7TeVdummy->GetYaxis()->SetMoreLogLabels(kTRUE);
        ratio7TeVdummy->GetYaxis()->SetNdivisions(505);
        ratio7TeVdummy->GetYaxis()->SetNoExponent(kTRUE);
        ratio7TeVdummy->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio7TeVdummy->GetXaxis()->SetNoExponent(kTRUE);
        ratio7TeVdummy->GetXaxis()->SetLabelFont(42);
        ratio7TeVdummy->GetYaxis()->SetLabelFont(42);
        ratio7TeVdummy->GetYaxis()->CenterTitle(kTRUE);
        ratio7TeVdummy->GetYaxis()->SetLabelOffset(+0.01);
        ratio7TeVdummy->GetXaxis()->SetTickLength(0.07);
        ratio7TeVdummy->DrawCopy();

        TGraphAsymmErrors* graphRatioPi0DSS147TeV            = (TGraphAsymmErrors*)graphPi0DSS147TeV->Clone();
        graphRatioPi0DSS147TeV                               = CalculateGraphErrRatioToFit (graphRatioPi0DSS147TeV, fitTCMInvXSectionPi0Plot7TeV);

        graphRatioPi0DSS147TeV->SetLineWidth(widthCommonFit);
        graphRatioPi0DSS147TeV->SetLineColor(colorNLO);
        graphRatioPi0DSS147TeV->SetLineStyle(1);
        graphRatioPi0DSS147TeV->SetFillStyle(1001);
        graphRatioPi0DSS147TeV->SetFillColor(colorNLO);
        if(plotTheorycurves)
          graphRatioPi0DSS147TeV->Draw("same,3");

        TH1D* histoRatioPythia8ToFit7TeV                     = (TH1D*) histoPythia8InvXSectionPi07TeV->Clone();
        histoRatioPythia8ToFit7TeV                           = CalculateHistoRatioToFit (histoRatioPythia8ToFit7TeV, fitTCMInvXSectionPi0Plot7TeV);
        DrawGammaSetMarker(histoRatioPythia8ToFit7TeV, 24, 1.5, pythia8color , pythia8color);
        histoRatioPythia8ToFit7TeV->SetLineWidth(widthCommonFit);
        if(plotTheorycurves)
        histoRatioPythia8ToFit7TeV->Draw("same,hist,l");


        TGraphErrors* graphRatioPythia8ToFit7TeV             = (TGraphErrors*) graphPythia8InvXSection7TeV->Clone();
        graphRatioPythia8ToFit7TeV                           = CalculateGraphErrRatioToFit (graphRatioPythia8ToFit7TeV, fitTCMInvXSectionPi0Plot7TeV);
        DrawGammaSetMarkerTGraphErr(graphRatioPythia8ToFit7TeV, 0, 0, pythia8color , pythia8color, widthLinesBoxes, kTRUE, pythia8color);
        if(plotTheorycurves)
          graphRatioPythia8ToFit7TeV->Draw("3,same");


        TGraphAsymmErrors* graphRatioPi0CombCombFitStatA7TeV    = (TGraphAsymmErrors*)csGraphs[1]->Clone();
        graphRatioPi0CombCombFitStatA7TeV                       = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitStatA7TeV, fitTCMInvXSectionPi0Plot7TeV);
        TGraphAsymmErrors* graphRatioPi0CombCombFitStatA7TeV_WOXErr = (TGraphAsymmErrors*) graphRatioPi0CombCombFitStatA7TeV->Clone("graphRatioPi0CombCombFitStatA7TeV_WOXErr");
        ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombCombFitStatA7TeV_WOXErr);

        TGraphAsymmErrors* graphRatioPi0CombCombFitSysA7TeV     = (TGraphAsymmErrors*)csGraphsSys[1]->Clone();
        graphRatioPi0CombCombFitSysA7TeV                        = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitSysA7TeV, fitTCMInvXSectionPi0Plot7TeV);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA7TeV, markerstylesFULL[1], markersizeFULL[1], colorDataFULL[1], colorDataFULL[1], widthLinesBoxes, kTRUE, 0);
        graphRatioPi0CombCombFitSysA7TeV->SetLineWidth(0);
        graphRatioPi0CombCombFitSysA7TeV->Draw("2,same");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA7TeV_WOXErr, markerstylesFULL[1], markersizeFULL[1], colorDataFULL[1], colorDataFULL[1], widthLinesBoxes, kFALSE);
        graphRatioPi0CombCombFitStatA7TeV_WOXErr->SetLineWidth(widthLinesBoxes);
        graphRatioPi0CombCombFitStatA7TeV_WOXErr->Draw("p,same");


        DrawGammaLines(0.23, maxX , 1., 1.,0.5, kGray+2);
        
    padInvSection2760GeVRatio->cd();
    padInvSection2760GeVRatio->SetLogx(1);
        TH2F * ratio276TeVdummy               = new TH2F("ratio276TeVdummy","ratio276TeVdummy",1000,0.23,maxX,1000,0.6,1.95);
        SetStyleHistoTH2ForGraphs(ratio276TeVdummy, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{TCM fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
                                  0.85*textsizeLabelsXSecDown,0.9*textsizeLabelsXSecDown, 1,0.23/(textsizeFacXSecMiddle*marginXSec), 510, 505);
        ratio276TeVdummy->GetYaxis()->SetMoreLogLabels(kTRUE);
        ratio276TeVdummy->GetYaxis()->SetNdivisions(505);
        ratio276TeVdummy->GetYaxis()->SetNoExponent(kTRUE);
        ratio276TeVdummy->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio276TeVdummy->GetXaxis()->SetNoExponent(kTRUE);
        ratio276TeVdummy->GetXaxis()->SetLabelFont(42);
        ratio276TeVdummy->GetYaxis()->SetLabelFont(42);
        ratio276TeVdummy->GetYaxis()->CenterTitle(kTRUE);
        ratio276TeVdummy->GetYaxis()->SetLabelOffset(+0.01);
        ratio276TeVdummy->GetXaxis()->SetTickLength(0.07);
        ratio276TeVdummy->DrawCopy();

        TGraphAsymmErrors* graphRatioPi0DSS07276TeV            = (TGraphAsymmErrors*)graphPi0DSS14276TeV->Clone();
        graphRatioPi0DSS07276TeV                               = CalculateGraphErrRatioToFit (graphRatioPi0DSS07276TeV, fitTCMInvXSectionPi0Plot276TeV);

        graphRatioPi0DSS07276TeV->SetLineWidth(widthCommonFit);
        graphRatioPi0DSS07276TeV->SetLineColor(colorNLO);
        graphRatioPi0DSS07276TeV->SetLineStyle(1);
        graphRatioPi0DSS07276TeV->SetFillStyle(1001);
        graphRatioPi0DSS07276TeV->SetFillColor(colorNLO);
        if(plotTheorycurves)
          graphRatioPi0DSS07276TeV->Draw("same,3");

        TH1D* histoRatioPythia8ToFit276TeV                     = (TH1D*) histoPythia8InvXSectionPi0276TeV->Clone();
        histoRatioPythia8ToFit276TeV                           = CalculateHistoRatioToFit (histoRatioPythia8ToFit276TeV, fitTCMInvXSectionPi0Plot276TeV);
        DrawGammaSetMarker(histoRatioPythia8ToFit276TeV, 24, 1.5, pythia8color , pythia8color);
        histoRatioPythia8ToFit276TeV->SetLineWidth(widthCommonFit);
        if(plotTheorycurves)
          histoRatioPythia8ToFit276TeV->Draw("same,hist,l");


        TGraphErrors* graphRatioPythia8ToFit276TeV             = (TGraphErrors*) graphPythia8InvXSection276TeV->Clone();
        graphRatioPythia8ToFit276TeV                           = CalculateGraphErrRatioToFit (graphRatioPythia8ToFit276TeV, fitTCMInvXSectionPi0Plot276TeV);
        DrawGammaSetMarkerTGraphErr(graphRatioPythia8ToFit276TeV, 0, 0, pythia8color , pythia8color, widthLinesBoxes, kTRUE, pythia8color);
        if(plotTheorycurves)
          graphRatioPythia8ToFit276TeV->Draw("3,same");


        TGraphAsymmErrors* graphRatioPi0CombCombFitStatA276TeV    = (TGraphAsymmErrors*)csGraphs[3]->Clone();
        graphRatioPi0CombCombFitStatA276TeV                       = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitStatA276TeV, fitTCMInvXSectionPi0Plot276TeV);
        TGraphAsymmErrors* graphRatioPi0CombCombFitStatA276TeV_WOXErr = (TGraphAsymmErrors*) graphRatioPi0CombCombFitStatA276TeV->Clone("graphRatioPi0CombCombFitStatA276TeV_WOXErr");
        ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombCombFitStatA276TeV_WOXErr);

        TGraphAsymmErrors* graphRatioPi0CombCombFitSysA276TeV     = (TGraphAsymmErrors*)csGraphsSys[3]->Clone();
        graphRatioPi0CombCombFitSysA276TeV                        = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitSysA276TeV, fitTCMInvXSectionPi0Plot276TeV);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA276TeV, markerstylesFULL[3], markersizeFULL[3], colorDataFULL[3], colorDataFULL[3], widthLinesBoxes, kTRUE, 0);
        graphRatioPi0CombCombFitSysA276TeV->SetLineWidth(0);
        graphRatioPi0CombCombFitSysA276TeV->Draw("2,same");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA276TeV_WOXErr, markerstylesFULL[3], markersizeFULL[3], colorDataFULL[3], colorDataFULL[3], widthLinesBoxes, kFALSE);
        graphRatioPi0CombCombFitStatA276TeV_WOXErr->SetLineWidth(widthLinesBoxes);
        graphRatioPi0CombCombFitStatA276TeV_WOXErr->Draw("p,same");


        DrawGammaLines(0.23, maxX , 1., 1.,0.5, kGray+2);

    padInvSection900GeVRatio->cd();
    padInvSection900GeVRatio->SetLogx(1);
        TH2F * ratio900GeVdummy            = new TH2F("ratio900GeVdummy","ratio900GeVdummy",1000,0.23,maxX,1000,0.6,1.95);
        SetStyleHistoTH2ForGraphs(ratio900GeVdummy, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{Tsallis fit}", 0.85*textsizeLabelsXSecDown2, textsizeLabelsXSecDown2,
                                  0.85*textsizeLabelsXSecDown2,0.9*textsizeLabelsXSecDown2, 0.9,0.23/(textsizeFacXSecDown2*marginXSec), 510, 505);
        ratio900GeVdummy->GetYaxis()->SetMoreLogLabels(kTRUE);
        ratio900GeVdummy->GetYaxis()->SetNdivisions(505);
        ratio900GeVdummy->GetYaxis()->SetNoExponent(kTRUE);
        ratio900GeVdummy->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio900GeVdummy->GetXaxis()->SetNoExponent(kTRUE);
        ratio900GeVdummy->GetXaxis()->SetLabelFont(42);
        ratio900GeVdummy->GetYaxis()->SetLabelFont(42);
        ratio900GeVdummy->GetYaxis()->CenterTitle(kTRUE);
        ratio900GeVdummy->GetYaxis()->SetLabelOffset(+0.01);
        ratio900GeVdummy->GetXaxis()->SetTickLength(0.06);
        ratio900GeVdummy->GetYaxis()->SetTickLength(0.04);
        ratio900GeVdummy->DrawCopy();


        TH1D* histoRatioPythia8ToFit900GeV                     = (TH1D*) histoPythia8InvXSectionPi0900GeV->Clone();
        histoRatioPythia8ToFit900GeV                           = CalculateHistoRatioToFit (histoRatioPythia8ToFit900GeV, fitInvXSectionPi0900GeV);
        DrawGammaSetMarker(histoRatioPythia8ToFit900GeV, 24, 1.5, pythia8color , pythia8color);
        histoRatioPythia8ToFit900GeV->SetLineWidth(widthCommonFit);
        if(plotTheorycurves)
          histoRatioPythia8ToFit900GeV->Draw("same,hist,l");


        TGraphErrors* graphRatioPythia8ToFit900GeV             = (TGraphErrors*) graphPythia8InvXSection900GeV->Clone();
        graphRatioPythia8ToFit900GeV                           = CalculateGraphErrRatioToFit (graphRatioPythia8ToFit900GeV, fitInvXSectionPi0900GeV);
        DrawGammaSetMarkerTGraphErr(graphRatioPythia8ToFit900GeV, 0, 0, pythia8color , pythia8color, widthLinesBoxes, kTRUE, pythia8color);
        if(plotTheorycurves)
          graphRatioPythia8ToFit900GeV->Draw("3,same");


        TGraphAsymmErrors* graphRatioPi0CombCombFitStatA900GeV    = (TGraphAsymmErrors*)csGraphs[0]->Clone();
        graphRatioPi0CombCombFitStatA900GeV                       = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitStatA900GeV, fitInvXSectionPi0900GeV);
        TGraphAsymmErrors* graphRatioPi0CombCombFitStatA900GeV_WOXErr = (TGraphAsymmErrors*) graphRatioPi0CombCombFitStatA900GeV->Clone("graphRatioPi0CombCombFitStatA900GeV_WOXErr");
        ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombCombFitStatA900GeV_WOXErr);

        TGraphAsymmErrors* graphRatioPi0CombCombFitSysA900GeV     = (TGraphAsymmErrors*)csGraphsSys[0]->Clone();
        graphRatioPi0CombCombFitSysA900GeV                        = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitSysA900GeV, fitInvXSectionPi0900GeV);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA900GeV, markerstylesFULL[0], markersizeFULL[0], colorDataFULL[0], colorDataFULL[0], widthLinesBoxes, kTRUE, 0);
        graphRatioPi0CombCombFitSysA900GeV->SetLineWidth(0);
        graphRatioPi0CombCombFitSysA900GeV->Draw("2,same");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA900GeV_WOXErr, markerstylesFULL[0], markersizeFULL[0], colorDataFULL[0], colorDataFULL[0], widthLinesBoxes, kFALSE);
        graphRatioPi0CombCombFitStatA900GeV_WOXErr->SetLineWidth(widthLinesBoxes);
        graphRatioPi0CombCombFitStatA900GeV_WOXErr->Draw("p,same");

        DrawGammaLines(0.23, maxX , 1., 1.,0.5, kGray+2);

    canvasInvSectionPaper->Print(Form("%s/Pi0_InvCrossSectionWithRatios.%s",outputDir.Data(),suffix.Data()));



    // SINGLE PLOTS FOR SLIDES

    Double_t arrayBoundariesX1_XSec2[2];
    Double_t arrayBoundariesY1_XSec2[3];
    Double_t relativeMarginsXXSec2[3];
    Double_t relativeMarginsYXSec2[3];
    textSizeLabelsPixel = 48;
    ReturnCorrectValuesForCanvasScaling(1250,1300, 1, 1,0.135, 0.005, 0.003,0.09,arrayBoundariesX1_XSec2,arrayBoundariesY1_XSec2,relativeMarginsXXSec2,relativeMarginsYXSec2);

    TCanvas* canvasInvSectionPaperSingle      = new TCanvas("canvasInvSectionPaperSingle","",0,0,1250,1250);  // gives the page size
    DrawGammaCanvasSettings( canvasInvSectionPaperSingle,  0.13, 0.02, 0.03, 0.06);

    TPad* padInvSectionSpecSingle             = new TPad("padInvSectionSpecSingle", "", arrayBoundariesX1_XSec2[0], arrayBoundariesY1_XSec2[3], arrayBoundariesX1_XSec2[1], arrayBoundariesY1_XSec2[0],-1, -1, -2);
    DrawGammaPadSettings( padInvSectionSpecSingle, relativeMarginsXXSec2[0], relativeMarginsXXSec2[2], relativeMarginsYXSec2[0], relativeMarginsYXSec2[2]);
    padInvSectionSpecSingle->Draw();
    marginXSec                 = relativeMarginsXXSec2[0]*1250;
    textsizeLabelsXSecUp       = 0;
    textsizeFacXSecUp          = 0;
    if (padInvSectionSpecSingle->XtoPixel(padInvSectionSpecSingle->GetX2()) < padInvSectionSpecSingle->YtoPixel(padInvSectionSpecSingle->GetY1())){
        textsizeLabelsXSecUp            = (Double_t)textSizeLabelsPixel/padInvSectionSpecSingle->XtoPixel(padInvSectionSpecSingle->GetX2()) ;
        textsizeFacXSecUp               = (Double_t)1./padInvSectionSpecSingle->XtoPixel(padInvSectionSpecSingle->GetX2()) ;
    } else {
        textsizeLabelsXSecUp            = (Double_t)textSizeLabelsPixel/padInvSectionSpecSingle->YtoPixel(padInvSectionSpecSingle->GetY1());
        textsizeFacXSecUp               = (Double_t)1./padInvSectionSpecSingle->YtoPixel(padInvSectionSpecSingle->GetY1());
    }

     padInvSectionSpecSingle->cd();
    padInvSectionSpecSingle->SetLogy(1);
    padInvSectionSpecSingle->SetLogx(1);
    histo2DXSectionPi0->Draw();
      if(plotTheorycurves){
        graphPi0DSS148TeV->Draw("same,e3");
        graphPi0DSS147TeV->Draw("same,e3");
        graphPi0DSS14276TeV->Draw("same,e3");
        graphPythia8InvXSection8TeV->Draw("3,same");
        graphPythia8InvXSection7TeV->Draw("3,same");
        graphPythia8InvXSection276TeV->Draw("3,same");
        graphPythia8InvXSection900GeV->Draw("3,same");
        histoPythia8InvXSectionPi08TeV->Draw("same,hist,l");
        histoPythia8InvXSectionPi07TeV->Draw("same,hist,l");
        histoPythia8InvXSectionPi0276TeV->Draw("same,hist,l");
        histoPythia8InvXSectionPi0900GeV->Draw("same,hist,l");
      }

        fitTCMInvXSectionPi0Plot8TeV->Draw("same");

        fitTCMInvXSectionPi0Plot7TeV->Draw("same");
        fitTCMInvXSectionPi0Plot276TeV->Draw("same");

        fitTCMInvXSectionPi0Plot900GeV->Draw("same");


        // Tsallis fit

        fitInvXSectionPi08TeV->Draw("same");

        fitInvXSectionPi07TeV->Draw("same");
        fitInvXSectionPi0276TeV->Draw("same");

        fitInvXSectionPi0900GeV->Draw("same");


        // TF1* fitexp8TeV   = FitObject("e","fitexponential8TeV","Pi0",csGraphs[2],2.,30. ,NULL,"QNRMEX0+");
        // fitexp8TeV->Draw("same");

        csGraphsSys[2]     ->Draw("E2same");
        csGraphsSys[1]     ->Draw("E2same");
        csGraphsSys[3]     ->Draw("E2same");
        csGraphsSys[0]     ->Draw("E2same");
        csGraphs[2]->Draw("p,same,z");
        csGraphs[1]->Draw("p,same,z");
        csGraphs[3]->Draw("p,same,z");
        csGraphs[0]->Draw("p,same,z");



        drawLatexAdd("ALICE",rightalignDouble,0.91,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE);
        drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",rightalignDouble,0.86,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE);


        TLegend* legendXsectionPaperSingle    = GetAndSetLegend2(0.17, 0.13, 0.5, 0.23+0.05*4, textSizeLabelsPixel);
        legendXsectionPaperSingle->SetNColumns(1);
        legendXsectionPaperSingle->SetMargin(0.2);
        legendXsectionPaperSingle->AddEntry(csGraphsSys[2],Form("%s (x 10^{3})",nameMeasGlobal[2].Data()),"pf");
        legendXsectionPaperSingle->AddEntry(csGraphsSys[1],Form("%s (x 10^{2})",nameMeasGlobal[1].Data()),"pf");
        legendXsectionPaperSingle->AddEntry(csGraphsSys[3],Form("%s (x 10)",nameMeasGlobal[3].Data()),"pf");
        legendXsectionPaperSingle->AddEntry(csGraphsSys[0],nameMeasGlobal[0].Data(),"pf");
        legendXsectionPaperSingle->AddEntry(fitTCMInvXSectionPi0Plot8TeV,"TCM fit","l");
        legendXsectionPaperSingle->AddEntry(fitInvXSectionPi08TeV,"Tsallis fit","l");

        // legendXsectionPaperSingle->AddEntry(histoPythia8InvXSectionPi08TeV,"PYTHIA 8.2, Monash 2013","l");
        // legendXsectionPaperSingle->AddEntry(graphPi0DSS148TeV,  "NLO, PDF:MSTW08 - FF:DSS14", "f");
        legendXsectionPaperSingle->Draw();
        //
        // drawLatexAdd("x10^{2}",0.24,0.94,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE,colorDataFULL[2]);
        // drawLatexAdd("x10^{1}",0.24,0.835,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE,colorDataFULL[1]);

    canvasInvSectionPaperSingle->Print(Form("%s/Pi0_InvCrossSectionSingle.%s",outputDir.Data(),suffix.Data()));



     Double_t arrayBoundariesX1_XSec3[2];
    Double_t arrayBoundariesY1_XSec3[5];
    Double_t relativeMarginsXXSec3[3];
    Double_t relativeMarginsYXSec3[3];
    textSizeLabelsPixel = 48;
    ReturnCorrectValuesForCanvasScaling(1250,1300, 1, 4,0.135, 0.005, 0.005,0.082,arrayBoundariesX1_XSec3,arrayBoundariesY1_XSec3,relativeMarginsXXSec3,relativeMarginsYXSec3);

    TCanvas* canvasInvSectionPaperSingleRatio      = new TCanvas("canvasInvSectionPaperSingleRatio","",0,0,1250,1250);  // gives the page size
    DrawGammaCanvasSettings( canvasInvSectionPaperSingleRatio,  0.13, 0.02, 0.03, 0.06);

    TPad* padInvSectionNLORatioSingle         = new TPad("padInvSectionNLORatioSingle", "", arrayBoundariesX1_XSec3[0], arrayBoundariesY1_XSec3[1], arrayBoundariesX1_XSec3[1], arrayBoundariesY1_XSec3[0],-1, -1, -2);
    DrawGammaPadSettings( padInvSectionNLORatioSingle, relativeMarginsXXSec3[0], relativeMarginsXXSec3[2], relativeMarginsYXSec3[0], relativeMarginsYXSec3[1]);
    padInvSectionNLORatioSingle->Draw();
    textsizeLabelsXSecMiddle   = 0;
    textsizeFacXSecMiddle      = 0;
    if (padInvSectionNLORatioSingle->XtoPixel(padInvSectionNLORatioSingle->GetX2()) < padInvSectionNLORatioSingle->YtoPixel(padInvSectionNLORatioSingle->GetY1())){
        textsizeLabelsXSecMiddle        = (Double_t)textSizeLabelsPixel/padInvSectionNLORatioSingle->XtoPixel(padInvSectionNLORatioSingle->GetX2()) ;
        textsizeFacXSecMiddle           = (Double_t)1./padInvSectionNLORatioSingle->XtoPixel(padInvSectionNLORatioSingle->GetX2()) ;
    } else {
        textsizeLabelsXSecMiddle        = (Double_t)textSizeLabelsPixel/padInvSectionNLORatioSingle->YtoPixel(padInvSectionNLORatioSingle->GetY1());
        textsizeFacXSecMiddle           = (Double_t)1./padInvSectionNLORatioSingle->YtoPixel(padInvSectionNLORatioSingle->GetY1());
    }

    TPad* padInvSectionPythiaRatioSignle      = new TPad("padInvSectionPythiaRatioSignle", "", arrayBoundariesX1_XSec3[0], arrayBoundariesY1_XSec3[2], arrayBoundariesX1_XSec3[1], arrayBoundariesY1_XSec3[1],-1, -1, -2);
    DrawGammaPadSettings( padInvSectionPythiaRatioSignle, relativeMarginsXXSec3[0], relativeMarginsXXSec3[2], relativeMarginsYXSec3[1], relativeMarginsYXSec3[1]);
    padInvSectionPythiaRatioSignle->Draw();
    textsizeLabelsXSecDown     = 0;
    textsizeFacXSecDown        = 0;
    if (padInvSectionPythiaRatioSignle->XtoPixel(padInvSectionPythiaRatioSignle->GetX2()) < padInvSectionPythiaRatioSignle->YtoPixel(padInvSectionPythiaRatioSignle->GetY1())){
        textsizeLabelsXSecDown          = (Double_t)textSizeLabelsPixel/padInvSectionPythiaRatioSignle->XtoPixel(padInvSectionPythiaRatioSignle->GetX2()) ;
        textsizeFacXSecDown             = (Double_t)1./padInvSectionPythiaRatioSignle->XtoPixel(padInvSectionPythiaRatioSignle->GetX2()) ;
    } else {
        textsizeLabelsXSecDown          = (Double_t)textSizeLabelsPixel/padInvSectionPythiaRatioSignle->YtoPixel(padInvSectionPythiaRatioSignle->GetY1());
        textsizeFacXSecDown             = (Double_t)1./padInvSectionPythiaRatioSignle->YtoPixel(padInvSectionPythiaRatioSignle->GetY1());
    }

    TPad* padInvSection276TeVRatioSignle      = new TPad("padInvSection276TeVRatioSignle", "", arrayBoundariesX1_XSec3[0], arrayBoundariesY1_XSec3[3], arrayBoundariesX1_XSec3[1], arrayBoundariesY1_XSec3[2],-1, -1, -2);
    DrawGammaPadSettings( padInvSection276TeVRatioSignle, relativeMarginsXXSec3[0], relativeMarginsXXSec3[2], relativeMarginsYXSec3[1], relativeMarginsYXSec3[1]);
    padInvSection276TeVRatioSignle->Draw();
    textsizeLabelsXSecDown26     = 0;
    textsizeFacXSecDown26        = 0;
    if (padInvSection276TeVRatioSignle->XtoPixel(padInvSection276TeVRatioSignle->GetX2()) < padInvSection276TeVRatioSignle->YtoPixel(padInvSection276TeVRatioSignle->GetY1())){
        textsizeLabelsXSecDown26          = (Double_t)textSizeLabelsPixel/padInvSection276TeVRatioSignle->XtoPixel(padInvSection276TeVRatioSignle->GetX2()) ;
        textsizeFacXSecDown26             = (Double_t)1./padInvSection276TeVRatioSignle->XtoPixel(padInvSection276TeVRatioSignle->GetX2()) ;
    } else {
        textsizeLabelsXSecDown26          = (Double_t)textSizeLabelsPixel/padInvSection276TeVRatioSignle->YtoPixel(padInvSection276TeVRatioSignle->GetY1());
        textsizeFacXSecDown26             = (Double_t)1./padInvSection276TeVRatioSignle->YtoPixel(padInvSection276TeVRatioSignle->GetY1());
    }
    TPad* padInvSection900GeVRatioSingle      = new TPad("padInvSection900GeVRatioSingle", "", arrayBoundariesX1_XSec3[0], arrayBoundariesY1_XSec3[4], arrayBoundariesX1_XSec3[1], arrayBoundariesY1_XSec3[3],-1, -1, -2);
    DrawGammaPadSettings( padInvSection900GeVRatioSingle, relativeMarginsXXSec3[0], relativeMarginsXXSec3[2], relativeMarginsYXSec3[1], relativeMarginsYXSec3[2]);
    padInvSection900GeVRatioSingle->Draw();
    textsizeLabelsXSecDown2     = 0;
    textsizeFacXSecDown2        = 0;
    if (padInvSection900GeVRatioSingle->XtoPixel(padInvSection900GeVRatioSingle->GetX2()) < padInvSection900GeVRatioSingle->YtoPixel(padInvSection900GeVRatioSingle->GetY1())){
        textsizeLabelsXSecDown2          = (Double_t)textSizeLabelsPixel/padInvSection900GeVRatioSingle->XtoPixel(padInvSection900GeVRatioSingle->GetX2()) ;
        textsizeFacXSecDown2             = (Double_t)1./padInvSection900GeVRatioSingle->XtoPixel(padInvSection900GeVRatioSingle->GetX2()) ;
    } else {
        textsizeLabelsXSecDown2          = (Double_t)textSizeLabelsPixel/padInvSection900GeVRatioSingle->YtoPixel(padInvSection900GeVRatioSingle->GetY1());
        textsizeFacXSecDown2             = (Double_t)1./padInvSection900GeVRatioSingle->YtoPixel(padInvSection900GeVRatioSingle->GetY1());
    }

    padInvSectionNLORatioSingle->cd();
    padInvSectionNLORatioSingle->SetLogx(1);
    ratio8TeVdummy->DrawCopy();
    ratio8TeVdummy->DrawCopy();
    if(plotTheorycurves){
      graphRatioPi0DSS148TeV->Draw("same,3");
      histoRatioPythia8ToFit8TeV->Draw("same,hist,l");
      graphRatioPythia8ToFit8TeV->Draw("3,same");
    }
    graphRatioPi0CombCombFitSysA->Draw("2,same");
    graphRatioPi0CombCombFitStatA_WOXErr->Draw("p,same");
    
    TLegend* legend8TeVsingleratio    = GetAndSetLegend2(0.96, 0.9, 0.9, 0.93, textSizeLabelsPixel);
    legend8TeVsingleratio->SetNColumns(1);
    legend8TeVsingleratio->SetMargin(0.2);
    legend8TeVsingleratio->SetTextAlign(33);
    legend8TeVsingleratio->AddEntry((TObject*)0,"pp, #sqrt{s} = 8 TeV","");
    legend8TeVsingleratio->Draw();
    
    DrawGammaLines(0.23, maxX , 1., 1.,0.5, kGray+2);

    padInvSectionPythiaRatioSignle->cd();
    padInvSectionPythiaRatioSignle->SetLogx(1);
    ratio7TeVdummy->DrawCopy();
    if(plotTheorycurves){
      graphRatioPi0DSS147TeV->Draw("same,3");
      histoRatioPythia8ToFit7TeV->Draw("same,hist,l");
      graphRatioPythia8ToFit7TeV->Draw("3,same");
    }
    /*
    // pass 2 input BEGIN
    TString fileNameDoubleRatioPass2            = "ExternalInput/PCM/data_GammaConversionResultsFullCorrectionNoBinShifting_PCM_020712.root";
    TFile* fileGammasPass2                            = new TFile(fileNameDoubleRatioPass2.Data());
    TH1D* pass2invcrosssectionstat                       = (TH1D*)fileGammasPass2->Get("Pi07TeV/InvCrossSectionPi0");
    TGraphAsymmErrors* graphpass2invcrosssectionstat                      = new TGraphAsymmErrors(pass2invcrosssectionstat);
    ProduceGraphAsymmWithoutXErrors(graphpass2invcrosssectionstat);
    TGraphAsymmErrors* graphpass2invcrosssectionsys     = (TGraphAsymmErrors*)fileGammasPass2->Get("Pi07TeV/InvCrossSectionPi0Sys");
    Color_t pass2Xseccolor = kGray+2;
    TGraphAsymmErrors* graphpass2invXsecratioStat    = (TGraphAsymmErrors*)graphpass2invcrosssectionstat->Clone();
    for (int j=0;j<graphpass2invXsecratioStat->GetN();j++){
        graphpass2invXsecratioStat->GetY()[j] *= 10;
        graphpass2invXsecratioStat->GetEYhigh()[j] *= 10;
        graphpass2invXsecratioStat->GetEYlow()[j] *= 10;
    }
    graphpass2invXsecratioStat                       = CalculateGraphErrRatioToFit(graphpass2invXsecratioStat, fitTCMInvXSectionPi0Plot7TeV);
    TGraphAsymmErrors* graphpass2invXsecratioStat_noXErr = (TGraphAsymmErrors*) graphpass2invXsecratioStat->Clone("graphpass2invXsecratioStat_noXErr");
    ProduceGraphAsymmWithoutXErrors(graphpass2invXsecratioStat_noXErr);

    TGraphAsymmErrors* graphpass2invXsecratioSys     = (TGraphAsymmErrors*)graphpass2invcrosssectionsys->Clone();
    for (int j=0;j<graphpass2invXsecratioSys->GetN();j++){
        graphpass2invXsecratioSys->GetY()[j] *= 10;
        graphpass2invXsecratioSys->GetEYhigh()[j] *= 10;
        graphpass2invXsecratioSys->GetEYlow()[j] *= 10;
    }
    graphpass2invXsecratioSys                        = CalculateGraphErrRatioToFit(graphpass2invXsecratioSys, fitTCMInvXSectionPi0Plot7TeV);
    DrawGammaSetMarkerTGraphAsym(graphpass2invXsecratioSys, 24, markersizeFULL[1], pass2Xseccolor, pass2Xseccolor, widthLinesBoxes, kTRUE, 0);
    graphpass2invXsecratioSys->SetLineWidth(0);
    // graphpass2invXsecratioSys->Draw("2,same");
    // graphpass2invXsecratioSys->Draw("2,same");
    DrawGammaSetMarkerTGraphAsym(graphpass2invXsecratioStat_noXErr, 24, markersizeFULL[1], pass2Xseccolor, pass2Xseccolor, widthLinesBoxes, kFALSE);
    graphpass2invXsecratioStat_noXErr->SetLineWidth(widthLinesBoxes);
    // graphpass2invXsecratioStat_noXErr->Draw("p,same");
    // pass2 input END
    */
    // histoRatioPythia8ToFit7TeV->Draw("same,hist,l");
    // graphRatioPythia8ToFit7TeV->Draw("3,same");
    graphRatioPi0CombCombFitSysA7TeV->Draw("2,same");
    graphRatioPi0CombCombFitStatA7TeV_WOXErr->Draw("p,same");
    
    TLegend* legend7TeVsingleratio    = GetAndSetLegend2(0.96, 0.9, 0.9, 0.93, textSizeLabelsPixel);
    legend7TeVsingleratio->SetNColumns(1);
    legend7TeVsingleratio->SetMargin(0.2);
    legend7TeVsingleratio->SetTextAlign(33);
    legend7TeVsingleratio->AddEntry((TObject*)0,"pp, #sqrt{s} = 7 TeV","");
    legend7TeVsingleratio->Draw();

    DrawGammaLines(0.23, maxX , 1., 1.,0.5, kGray+2);

    padInvSection276TeVRatioSignle->cd();
    padInvSection276TeVRatioSignle->SetLogx(1);
    ratio276TeVdummy->DrawCopy();
    if(plotTheorycurves){
      histoRatioPythia8ToFit276TeV->Draw("same,hist,l");
      graphRatioPythia8ToFit276TeV->Draw("3,same");
    }
    graphRatioPi0CombCombFitSysA276TeV->Draw("2,same");

    graphRatioPi0CombCombFitStatA276TeV_WOXErr->Draw("p,same");
    
    TLegend* legend276TeVsingleratio    = GetAndSetLegend2(0.96, 0.9, 0.9, 0.93, textSizeLabelsPixel);
    legend276TeVsingleratio->SetNColumns(1);
    legend276TeVsingleratio->SetMargin(0.2);
    legend276TeVsingleratio->SetTextAlign(33);
    legend276TeVsingleratio->AddEntry((TObject*)0,"pp, #sqrt{s} = 2.76 TeV","");
    legend276TeVsingleratio->Draw();
    
    DrawGammaLines(0.23, maxX , 1., 1.,0.5, kGray+2);
    padInvSection900GeVRatioSingle->cd();
    padInvSection900GeVRatioSingle->SetLogx(1);
    ratio900GeVdummy->DrawCopy();
    if(plotTheorycurves){
      histoRatioPythia8ToFit900GeV->Draw("same,hist,l");
      graphRatioPythia8ToFit900GeV->Draw("3,same");
    }

    graphRatioPi0CombCombFitSysA900GeV->Draw("2,same");

    graphRatioPi0CombCombFitStatA900GeV_WOXErr->Draw("p,same");

    TLegend* legend900GeVsingleratio    = GetAndSetLegend2(0.96, 0.9, 0.9, 0.95, textSizeLabelsPixel);
    legend900GeVsingleratio->SetNColumns(1);
    legend900GeVsingleratio->SetMargin(0.2);
    legend900GeVsingleratio->SetTextAlign(33);
    legend900GeVsingleratio->AddEntry((TObject*)0,"pp, #sqrt{s} = 0.9 TeV","");
    legend900GeVsingleratio->Draw();
    
    DrawGammaLines(0.23, maxX , 1., 1.,0.5, kGray+2);




    canvasInvSectionPaperSingleRatio->Print(Form("%s/Pi0_InvCrossSectionRatioSingle.%s",outputDir.Data(),suffix.Data()));
    // canvasInvSectionPaperSingleRatio->Print(Form("%s/Pi0_InvCrossSectionRatioSingle.%s",outputDir.Data(),suffix.Data()));









    //########################################################################################
    //########################################################################################
    //                                   COMBINATION RATIOS
    //########################################################################################
    //########################################################################################
    TString nameMeasGlobalx[11]                  = {"PCM", "PHOS", "EMCal", "PCM-EMCal"};
    TString nameMeasGlobal2[11]                  = {"PCM", "PHOS", "EMCAL", "PCMEMCAL"};
    Color_t colorDet[11];
    Style_t markerStyleDet[11];
    Size_t  markerSizeDet[11];

    for (Int_t i = 0; i < 11; i++){
        colorDet[i]                             = GetDefaultColorDiffDetectors(nameMeasGlobalx[i].Data(), kFALSE, kFALSE, kTRUE);
        markerStyleDet[i]                       = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobalx[i].Data(), kFALSE);
        markerSizeDet[i]                        = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobalx[i].Data(), kFALSE);
    }


    TGraphAsymmErrors * combgraphsstat8TeV[5];
    TGraphAsymmErrors * combgraphssys8TeV[5];
    TGraphAsymmErrors * combgraphsstat7TeV[5];
    TGraphAsymmErrors * combgraphssys7TeV[5];
    TGraphAsymmErrors * combgraphsstat900GeV[5];
    TGraphAsymmErrors * combgraphssys900GeV[5];
     for(Int_t i=0; i<4; i++){
        combgraphsstat8TeV[i]     = (TGraphAsymmErrors*)directoryPi08TeV->Get(Form("graphInvCrossSectionPi0%s8TeVStatErr",nameMeasGlobal2[i].Data()));
        combgraphssys8TeV[i]      = (TGraphAsymmErrors*)directoryPi08TeV->Get(Form("graphInvCrossSectionPi0%s8TeVSysErr",nameMeasGlobal2[i].Data()));
        for (int j=0;j<combgraphsstat8TeV[i]->GetN();j++){
            combgraphsstat8TeV[i]->GetY()[j] *= 100;
            combgraphsstat8TeV[i]->GetEYhigh()[j] *= 100;
            combgraphsstat8TeV[i]->GetEYlow()[j] *= 100;
        }
        for (int j=0;j<combgraphssys8TeV[i]->GetN();j++){
            combgraphssys8TeV[i]->GetY()[j] *= 100;
            combgraphssys8TeV[i]->GetEYhigh()[j] *= 100;
            combgraphssys8TeV[i]->GetEYlow()[j] *= 100;
        }
     }
     for(Int_t i=0; i<4; i++){
        combgraphsstat7TeV[i]     = (TGraphAsymmErrors*)directoryPi07TeV->Get(Form("graphInvCrossSectionPi0%s7TeVStatErr",nameMeasGlobal2[i].Data()));
        combgraphssys7TeV[i]      = (TGraphAsymmErrors*)directoryPi07TeV->Get(Form("graphInvCrossSectionPi0%s7TeVSysErr",nameMeasGlobal2[i].Data()));
        for (int j=0;j<combgraphsstat7TeV[i]->GetN();j++){
            combgraphsstat7TeV[i]->GetY()[j] *= 10;
            combgraphsstat7TeV[i]->GetEYhigh()[j] *= 10;
            combgraphsstat7TeV[i]->GetEYlow()[j] *= 10;
        }
        for (int j=0;j<combgraphssys7TeV[i]->GetN();j++){
            combgraphssys7TeV[i]->GetY()[j] *= 10;
            combgraphssys7TeV[i]->GetEYhigh()[j] *= 10;
            combgraphssys7TeV[i]->GetEYlow()[j] *= 10;
        }
     }
     for(Int_t i=0; i<2; i++){
        combgraphsstat900GeV[i]     = (TGraphAsymmErrors*)directoryPi0900GeV->Get(Form("graphInvCrossSectionPi0%s900GeVStatErr",nameMeasGlobal2[i].Data()));
        combgraphssys900GeV[i]      = (TGraphAsymmErrors*)directoryPi0900GeV->Get(Form("graphInvCrossSectionPi0%s900GeVSysErr",nameMeasGlobal2[i].Data()));
        for (int j=0;j<combgraphsstat900GeV[i]->GetN();j++){
            combgraphsstat900GeV[i]->GetY()[j] *= xSection900GeVx*recalcBarnx;
            combgraphsstat900GeV[i]->GetEYhigh()[j] *= xSection900GeVx*recalcBarnx;
            combgraphsstat900GeV[i]->GetEYlow()[j] *= xSection900GeVx*recalcBarnx;
        }
        for (int j=0;j<combgraphssys900GeV[i]->GetN();j++){
            combgraphssys900GeV[i]->GetY()[j] *= xSection900GeVx*recalcBarnx;
            combgraphssys900GeV[i]->GetEYhigh()[j] *= xSection900GeVx*recalcBarnx;
            combgraphssys900GeV[i]->GetEYlow()[j] *= xSection900GeVx*recalcBarnx;
        }
     }

    TGraphAsymmErrors* graphRatioPi0CombFitStatA8TeV[5];
    TGraphAsymmErrors* graphRatioPi0CombFitStatA8TeV_woxerr[5];
    TGraphAsymmErrors* graphRatioPi0CombFitSysA8TeV[5];
    TGraphAsymmErrors* graphRatioPi0CombFitStatA7TeV[5];
    TGraphAsymmErrors* graphRatioPi0CombFitStatA7TeV_woxerr[5];
    TGraphAsymmErrors* graphRatioPi0CombFitSysA7TeV[5];
    TGraphAsymmErrors* graphRatioPi0CombFitStatA900GeV[5];
    TGraphAsymmErrors* graphRatioPi0CombFitStatA900GeV_woxerr[5];
    TGraphAsymmErrors* graphRatioPi0CombFitSysA900GeV[5];

    padInvSectionNLORatioSingle->cd();
    padInvSectionNLORatioSingle->SetLogx(1);
    ratio8TeVdummy->GetYaxis()->SetTitle("#frac{Data}{Comb Fit}");
    ratio8TeVdummy->DrawCopy();
    for(Int_t i=0; i<4; i++){
        graphRatioPi0CombFitStatA8TeV[i]    = (TGraphAsymmErrors*)combgraphsstat8TeV[i]->Clone();
        graphRatioPi0CombFitStatA8TeV[i]                       = CalculateGraphErrRatioToFit(graphRatioPi0CombFitStatA8TeV[i], fitTCMInvXSectionPi08TeV);
        graphRatioPi0CombFitStatA8TeV_woxerr[i] = (TGraphAsymmErrors*) graphRatioPi0CombFitStatA8TeV[i]->Clone("graphRatioPi0CombFitStatA8TeV_woxerr");
        ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombFitStatA8TeV_woxerr[i]);

        graphRatioPi0CombFitSysA8TeV[i]     = (TGraphAsymmErrors*)combgraphssys8TeV[i]->Clone();
        graphRatioPi0CombFitSysA8TeV[i]                        = CalculateGraphErrRatioToFit(graphRatioPi0CombFitSysA8TeV[i], fitTCMInvXSectionPi08TeV);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombFitSysA8TeV[i], markerStyleDet[i], markerSizeDet[i], colorDet[i], colorDet[i], widthLinesBoxes, kTRUE, 0);
        graphRatioPi0CombFitSysA8TeV[i]->SetLineWidth(0);
        graphRatioPi0CombFitSysA8TeV[i]->Draw("2,same");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombFitStatA8TeV_woxerr[i], markerStyleDet[i], markerSizeDet[i], colorDet[i], colorDet[i], widthLinesBoxes, kFALSE);
        graphRatioPi0CombFitStatA8TeV_woxerr[i]->SetLineWidth(widthLinesBoxes);
        graphRatioPi0CombFitStatA8TeV_woxerr[i]->Draw("p,same");
    }
    DrawGammaLines(0.23, maxX , 1., 1.,0.5, kGray+2);
    drawLatexAdd(nameMeasGlobal[2].Data(),0.27,0.8,0.135,kFALSE,kFALSE);


    padInvSectionPythiaRatioSignle->cd();
    padInvSectionPythiaRatioSignle->SetLogx(1);
    ratio7TeVdummy->GetYaxis()->SetTitle("#frac{Data}{Comb Fit}");
    ratio7TeVdummy->DrawCopy();

    for(Int_t i=0; i<4; i++){
        graphRatioPi0CombFitStatA7TeV[i]    = (TGraphAsymmErrors*)combgraphsstat7TeV[i]->Clone();
        graphRatioPi0CombFitStatA7TeV[i]                       = CalculateGraphErrRatioToFit(graphRatioPi0CombFitStatA7TeV[i], fitTCMInvXSectionPi07TeV);
        graphRatioPi0CombFitStatA7TeV_woxerr[i] = (TGraphAsymmErrors*) graphRatioPi0CombFitStatA7TeV[i]->Clone("graphRatioPi0CombFitStatA7TeV_woxerr");
        ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombFitStatA7TeV_woxerr[i]);

        graphRatioPi0CombFitSysA7TeV[i]     = (TGraphAsymmErrors*)combgraphssys7TeV[i]->Clone();
        graphRatioPi0CombFitSysA7TeV[i]                        = CalculateGraphErrRatioToFit(graphRatioPi0CombFitSysA7TeV[i], fitTCMInvXSectionPi07TeV);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombFitSysA7TeV[i], markerStyleDet[i], markerSizeDet[i], colorDet[i], colorDet[i], widthLinesBoxes, kTRUE, 0);
        graphRatioPi0CombFitSysA7TeV[i]->SetLineWidth(0);
        graphRatioPi0CombFitSysA7TeV[i]->Draw("2,same");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombFitStatA7TeV_woxerr[i], markerStyleDet[i], markerSizeDet[i], colorDet[i], colorDet[i], widthLinesBoxes, kFALSE);
        graphRatioPi0CombFitStatA7TeV_woxerr[i]->SetLineWidth(widthLinesBoxes);
        graphRatioPi0CombFitStatA7TeV_woxerr[i]->Draw("p,same");
    }

    DrawGammaLines(0.23, maxX , 1., 1.,0.5, kGray+2);
    drawLatexAdd(nameMeasGlobal[1].Data(),0.27,0.8,0.135,kFALSE,kFALSE);
    // drawLatexAdd("ALICE",rightalignDouble,0.91,0.1,kFALSE,kFALSE,kTRUE);
        // drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",rightalignDouble,0.86,0.1,kFALSE,kFALSE,kTRUE);

    padInvSection900GeVRatioSingle->cd();
    padInvSection900GeVRatioSingle->SetLogx(1);
    ratio900GeVdummy->GetYaxis()->SetTitle("#frac{Data}{Comb Fit}");
    ratio900GeVdummy->DrawCopy();


    for(Int_t i=0; i<2; i++){
        graphRatioPi0CombFitStatA900GeV[i]    = (TGraphAsymmErrors*)combgraphsstat900GeV[i]->Clone();
        graphRatioPi0CombFitStatA900GeV[i]                       = CalculateGraphErrRatioToFit(graphRatioPi0CombFitStatA900GeV[i], fitTCMInvXSectionPi0900GeV);
        graphRatioPi0CombFitStatA900GeV_woxerr[i] = (TGraphAsymmErrors*) graphRatioPi0CombFitStatA900GeV[i]->Clone("graphRatioPi0CombFitStatA900GeV_woxerr");
        ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombFitStatA900GeV_woxerr[i]);

        graphRatioPi0CombFitSysA900GeV[i]     = (TGraphAsymmErrors*)combgraphssys900GeV[i]->Clone();
        graphRatioPi0CombFitSysA900GeV[i]                        = CalculateGraphErrRatioToFit(graphRatioPi0CombFitSysA900GeV[i], fitTCMInvXSectionPi0900GeV);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombFitSysA900GeV[i], markerStyleDet[i], markerSizeDet[i], colorDet[i], colorDet[i], widthLinesBoxes, kTRUE, 0);
        graphRatioPi0CombFitSysA900GeV[i]->SetLineWidth(0);
        graphRatioPi0CombFitSysA900GeV[i]->Draw("2,same");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombFitStatA900GeV_woxerr[i], markerStyleDet[i], markerSizeDet[i], colorDet[i], colorDet[i], widthLinesBoxes, kFALSE);
        graphRatioPi0CombFitStatA900GeV_woxerr[i]->SetLineWidth(widthLinesBoxes);
        graphRatioPi0CombFitStatA900GeV_woxerr[i]->Draw("p,same");
    }

    DrawGammaLines(0.23, maxX , 1., 1.,0.5, kGray+2);
    drawLatexAdd(nameMeasGlobal[0].Data(),0.27,0.82,0.11,kFALSE,kFALSE);

    Double_t reldiff = 0.13;
    Double_t rowsLegendOnlyPi0Ratio[6]      = {0.91,0.83,0.75,0.67,0.59,0.55};
        // Double_t rowsLegendOnlyPi0RatioAbs[6]   = {0.91,2.2,2.1,2.0,1.9,1.8};
        Double_t rowsLegendOnlyPi0RatioAbs[6]   = {0.91,1.7,1.7-reldiff,1.7-2*reldiff,1.7-3*reldiff,1.9};
        Double_t columnsLegendOnlyPi0Ratio[3]   = {0.65,0.82, 0.88};
        Double_t columnsLegendOnlyPi0RatioAbs[3]= {0.15,20.04, 27.37};
        Double_t lengthBox                      = 1.5;
        Double_t heightBox                      = 0.08/2;
        //****************** first Column **************************************************
        TLatex *textSingleMeasRatioPi0[10];
        for (Int_t i = 0; i < 4; i++){
            textSingleMeasRatioPi0[i]           = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[i+1],nameMeasGlobalx[i].Data());
            SetStyleTLatex( textSingleMeasRatioPi0[i], 0.85*textSizeLabelsPixel,4);
            textSingleMeasRatioPi0[i]->SetTextFont(43);
            textSingleMeasRatioPi0[i]->Draw();
        }

        //****************** second Column *************************************************
        TLatex *textStatOnlyRatioPi0            = new TLatex(columnsLegendOnlyPi0Ratio[1],rowsLegendOnlyPi0Ratio[0] ,"stat");
        SetStyleTLatex( textStatOnlyRatioPi0, 0.85*textSizeLabelsPixel,4);
        textStatOnlyRatioPi0->SetTextFont(43);
        textStatOnlyRatioPi0->Draw();
        TLatex *textSysOnlyRatioPi0             = new TLatex(columnsLegendOnlyPi0Ratio[2] ,rowsLegendOnlyPi0Ratio[0],"syst");
        SetStyleTLatex( textSysOnlyRatioPi0, 0.85*textSizeLabelsPixel,4);
        textSysOnlyRatioPi0->SetTextFont(43);
        textSysOnlyRatioPi0->Draw();

        TMarker* markerPi0OnlyRatio[10];
        for (Int_t i = 0; i < 4; i++){
            markerPi0OnlyRatio[i]               = CreateMarkerFromGraph(graphRatioPi0CombFitSysA7TeV[i],columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[i+1],1);
            markerPi0OnlyRatio[i]->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[i+1]);
        }

        TBox* boxPi0OnlyRatio[10];
        for (Int_t i = 0; i < 4; i++){
            boxPi0OnlyRatio[i]                  = CreateBoxFromGraph(graphRatioPi0CombFitSysA7TeV[i], columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[i+1]- heightBox,
                                                        columnsLegendOnlyPi0RatioAbs[2]+ 3*lengthBox, rowsLegendOnlyPi0RatioAbs[i+1]+ heightBox);
            boxPi0OnlyRatio[i]->Draw("l");
        }

           drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.77,0.28,0.135,kFALSE,kFALSE);


    canvasInvSectionPaperSingleRatio->Print(Form("%s/CombinationRatiosPi0.%s",outputDir.Data(),suffix.Data()));






    }


void plotCrossSectionFULLEta(TGraphAsymmErrors* csGraphs[],TGraphAsymmErrors* csGraphsSys[], Double_t yMin, Double_t yMax, TString meson, TString plotName)
{
    Bool_t plotTheorycurves   = kFALSE;
    Double_t xSection900GeVx         = 47.78*1e-3;
    Double_t recalcBarnx             = 1e12;

    TFile* input8TeVfile = new TFile("ThesisQAInputAllEnergies/CombinedResultsPaperPP8TeV_2017_07_13.root");
    TDirectory* directoryEta8TeV                     = (TDirectory*)input8TeVfile->Get("Eta8TeV");
    TFile* input7TeVfile = new TFile("ThesisQAInputAllEnergies/CombinedResultsPaperPP7TeV_2017_08_28.root");
    TDirectory* directoryEta7TeV                     = (TDirectory*)input7TeVfile->Get("Eta7TeV");
    TFile* input276TeVfile = new TFile("ThesisQAInputAllEnergies/2.76TeV/CombinedResultsPaperPP2760GeV_2017_03_01_FrediV2Clusterizer.root");
    TDirectory* directoryEta276TeV                     = (TDirectory*)input276TeVfile->Get("Eta2.76TeV");
    TFile* input900GeVfile = new TFile("ThesisQAInputAllEnergies/CombinedResultsPaperPP900GeV_2017_07_18.root");
    TDirectory* directoryEta900GeV                     = (TDirectory*)input900GeVfile->Get("Eta900GeV");
    //
    csGraphs[3]     = (TGraphAsymmErrors*)directoryEta276TeV      ->Get("graphInvCrossSectionEtaComb2760GeVAStatErr");
    csGraphsSys[3]     = (TGraphAsymmErrors*)directoryEta276TeV   ->Get("graphInvCrossSectionEtaComb2760GeVASysErr");
    csGraphs[2]     = (TGraphAsymmErrors*)directoryEta8TeV        ->Get("graphInvCrossSectionEtaComb8TeVAStatErr");
    csGraphsSys[2]     = (TGraphAsymmErrors*)directoryEta8TeV     ->Get("graphInvCrossSectionEtaComb8TeVASysErr");
        csGraphs[1]     = (TGraphAsymmErrors*)directoryEta7TeV    ->Get("graphInvCrossSectionEtaComb7TeVAStatErr");
        csGraphsSys[1]     = (TGraphAsymmErrors*)directoryEta7TeV ->Get("graphInvCrossSectionEtaComb7TeVASysErr");
        csGraphs[0]     = (TGraphAsymmErrors*)directoryEta900GeV  ->Get("graphInvCrossSectionEtaPCM900GeVStatErr");
        csGraphsSys[0]    = (TGraphAsymmErrors*)directoryEta900GeV->Get("graphInvCrossSectionEtaPCM900GeVSysErr");
    for(Int_t i = 0; i< 4; i++){
        if(!csGraphs[i])
          cout << "graphstat for " << i << " not found" << endl;
        if(!csGraphsSys[i])
          cout << "graphsys for " << i << " not found" << endl;
    }
    //     for (int j=0;j<csGraphs[0]->GetN();j++){
    //             csGraphs[0]->GetY()[j] *= xSection900GeVx*recalcBarnx;
    //             csGraphs[0]->GetEYhigh()[j] *= xSection900GeVx*recalcBarnx;
    //             csGraphs[0]->GetEYlow()[j] *= xSection900GeVx*recalcBarnx;
    //         }
    //     for (int j=0;j<csGraphsSys[0]->GetN();j++){
    //             csGraphsSys[0]->GetY()[j] *= xSection900GeVx*recalcBarnx;
    //             csGraphsSys[0]->GetEYhigh()[j] *= xSection900GeVx*recalcBarnx;
    //             csGraphsSys[0]->GetEYlow()[j] *= xSection900GeVx*recalcBarnx;
    //         }

    Style_t markerstylesFULL[4]                         ={34,20,33,29};
    Size_t markersizeFULL[4]                        ={2.3,2.3,3.4,3.4};
    Color_t colorDataFULL[4]                        ={kRed+2,kBlue+2,kGreen+2,kMagenta+2};
  
    //
  
            for (int j=0;j<csGraphsSys[2]->GetN();j++){
                csGraphsSys[2]->GetY()[j] *= 1000;
                csGraphsSys[2]->GetEYhigh()[j] *= 1000;
                csGraphsSys[2]->GetEYlow()[j] *= 1000;
            }
            DrawGammaSetMarkerTGraphAsym(csGraphsSys[2], markerstylesFULL[2], markersizeFULL[2], colorDataFULL[2], colorDataFULL[2], widthLinesBoxes, kTRUE);
            
            for (int j=0;j<csGraphsSys[1]->GetN();j++){
                csGraphsSys[1]->GetY()[j] *= 100;
                csGraphsSys[1]->GetEYhigh()[j] *= 100;
                csGraphsSys[1]->GetEYlow()[j] *= 100;
            }
            DrawGammaSetMarkerTGraphAsym(csGraphsSys[1], markerstylesFULL[1], markersizeFULL[1], colorDataFULL[1], colorDataFULL[1], widthLinesBoxes, kTRUE);
                for (int j=0;j<csGraphsSys[3]->GetN();j++){
                    csGraphsSys[3]->GetY()[j] *= 10;
                    csGraphsSys[3]->GetEYhigh()[j] *= 10;
                    csGraphsSys[3]->GetEYlow()[j] *= 10;
                }
            DrawGammaSetMarkerTGraphAsym(csGraphsSys[3], markerstylesFULL[3], markersizeFULL[3], colorDataFULL[3], colorDataFULL[3], widthLinesBoxes, kTRUE);
            DrawGammaSetMarkerTGraphAsym(csGraphsSys[0], markerstylesFULL[0], markersizeFULL[0], colorDataFULL[0], colorDataFULL[0], widthLinesBoxes, kTRUE);
            //
            // statistics
            
            for (int j=0;j<csGraphs[2]->GetN();j++){
                csGraphs[2]->GetY()[j] *= 1000;
                csGraphs[2]->GetEYhigh()[j] *= 1000;
                csGraphs[2]->GetEYlow()[j] *= 1000;
            }
            DrawGammaSetMarkerTGraph(csGraphs[2], markerstylesFULL[2], 0.7*markersizeFULL[2], colorDataFULL[2] , colorDataFULL[2]);
            
            for (int j=0;j<csGraphs[1]->GetN();j++){
                csGraphs[1]->GetY()[j] *= 100;
                csGraphs[1]->GetEYhigh()[j] *= 100;
                csGraphs[1]->GetEYlow()[j] *= 100;
            }
            DrawGammaSetMarkerTGraph(csGraphs[1], markerstylesFULL[1], 0.7*markersizeFULL[1], colorDataFULL[1] , colorDataFULL[1]);
            for (int j=0;j<csGraphs[3]->GetN();j++){
                csGraphs[3]->GetY()[j] *= 10;
                csGraphs[3]->GetEYhigh()[j] *= 10;
                csGraphs[3]->GetEYlow()[j] *= 10;
            }
            DrawGammaSetMarkerTGraph(csGraphs[3], markerstylesFULL[3], 0.7*markersizeFULL[3], colorDataFULL[3] , colorDataFULL[3]);

            DrawGammaSetMarkerTGraph(csGraphs[0], markerstylesFULL[0], 0.7*markersizeFULL[0], colorDataFULL[0] , colorDataFULL[0]);





    Double_t textSizeLabelsPixel                     = 48;
    TCanvas* canvasRatioPP                  = new TCanvas("canvasRatioPP","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioPP,  0.12, 0.01, 0.01, 0.11);
    canvasRatioPP->cd();
    canvasRatioPP->SetLogx();

    Double_t textsizeLabelsPP                    = 0;
    Double_t textsizeFacPP                       = 0;
    if (canvasRatioPP->XtoPixel(canvasRatioPP->GetX2()) <canvasRatioPP->YtoPixel(canvasRatioPP->GetY1()) ){
        textsizeLabelsPP                = (Double_t)textSizeLabelsPixel/canvasRatioPP->XtoPixel(canvasRatioPP->GetX2()) ;
        textsizeFacPP                   = (Double_t)1./canvasRatioPP->XtoPixel(canvasRatioPP->GetX2()) ;
    } else {
        textsizeLabelsPP                    = (Double_t)textSizeLabelsPixel/canvasRatioPP->YtoPixel(canvasRatioPP->GetY1());
        textsizeFacPP                       = (Double_t)1./canvasRatioPP->YtoPixel(canvasRatioPP->GetY1());
    }
    cout << textsizeLabelsPP << endl;

    Double_t mesonMassExpectPi0                 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    Double_t mesonMassExpectEta                 = TDatabasePDG::Instance()->GetParticle(221)->Mass();
    Color_t  colorComb                          = kMagenta+2;
    Style_t  markerStyleComb                    = 20;
    Size_t   markersizeFULLComb                     = 2;
    Color_t  colorCGC                           = kCyan-8;
    Color_t  colorNLO                           = kAzure-4;
    Style_t  styleMarkerNLOMuHalf               = 24;
    Style_t  styleMarkerNLOMuOne                = 27;
    Style_t  styleMarkerNLOMuTwo                = 30;
    Style_t  styleLineNLOMuHalf                 = 8;
    Style_t  styleLineNLOMuOne                  = 7;
    Style_t  styleLineNLOMuTwo                  = 4;
    Style_t  styleLineNLOMuTwoBKK               = 3;
    Style_t  styleLineNLOMuTwoDSS               = 6;
    Size_t   sizeMarkerNLO                      = 1;
    Width_t  widthLineNLO                       = 2.;


    TString fileNameTheory                      = "ExternalInput/Theory/TheoryCompilationPP.root";
    TFile* fileTheoryCompilation                            = new TFile(fileNameTheory.Data());

    // 8TeV Pythia8 Monash2013:
    TH1F* histoPythia8InvXSectionEta8TeV                       = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoEta8TeV");
    histoPythia8InvXSectionEta8TeV->GetXaxis()->SetRangeUser(0.4,35);
    histoPythia8InvXSectionEta8TeV->Scale(1000);
    // 7TeV Pythia8 Monash2013:
    TH1F* histoPythia8InvXSectionEta7TeV                       = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013Eta7TeV");
    histoPythia8InvXSectionEta7TeV->GetXaxis()->SetRangeUser(0.4,35);
    histoPythia8InvXSectionEta7TeV->Scale(100);
    // 7TeV Pythia8 Monash2013:
    TH1F* histoPythia8InvXSectionEta276TeV                       = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013Eta2760GeV");
    histoPythia8InvXSectionEta276TeV->GetXaxis()->SetRangeUser(0.4,35);
    histoPythia8InvXSectionEta276TeV->Scale(10);
    // 0.9TeV Pythia8 Monash2013:
    TH1F* histoPythia8InvXSectionEta900GeV                       = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013Eta900GeV");
    histoPythia8InvXSectionEta900GeV->GetXaxis()->SetRangeUser(0.8,3.);

    TGraphErrors* graphPythia8InvXSection8TeV               = new TGraphErrors((TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoEta8TeV"));
    for (int j=0;j<graphPythia8InvXSection8TeV->GetN();j++){
        graphPythia8InvXSection8TeV->GetY()[j] *= 1000;
        graphPythia8InvXSection8TeV->GetEY()[j] *= 1000;
    }
    while(graphPythia8InvXSection8TeV->GetX()[0] < 0.4) graphPythia8InvXSection8TeV->RemovePoint(0);
    while(graphPythia8InvXSection8TeV->GetX()[graphPythia8InvXSection8TeV->GetN()-1] > 35.) graphPythia8InvXSection8TeV->RemovePoint(graphPythia8InvXSection8TeV->GetN()-1);

    TGraphErrors* graphPythia8InvXSection7TeV               = new TGraphErrors((TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013Eta7TeV"));
    for (int j=0;j<graphPythia8InvXSection7TeV->GetN();j++){
        graphPythia8InvXSection7TeV->GetY()[j] *= 100;
        graphPythia8InvXSection7TeV->GetEY()[j] *= 100;
    }
    while(graphPythia8InvXSection7TeV->GetX()[0] < 0.4) graphPythia8InvXSection7TeV->RemovePoint(0);
    while(graphPythia8InvXSection7TeV->GetX()[graphPythia8InvXSection7TeV->GetN()-1] > 35.) graphPythia8InvXSection7TeV->RemovePoint(graphPythia8InvXSection7TeV->GetN()-1);
    
    TGraphErrors* graphPythia8InvXSection276TeV               = new TGraphErrors((TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013Eta2760GeV"));
    for (int j=0;j<graphPythia8InvXSection276TeV->GetN();j++){
        graphPythia8InvXSection276TeV->GetY()[j] *= 10;
        graphPythia8InvXSection276TeV->GetEY()[j] *= 10;
    }
    while(graphPythia8InvXSection276TeV->GetX()[0] < 0.4) graphPythia8InvXSection276TeV->RemovePoint(0);
    while(graphPythia8InvXSection276TeV->GetX()[graphPythia8InvXSection276TeV->GetN()-1] > 35.) graphPythia8InvXSection276TeV->RemovePoint(graphPythia8InvXSection276TeV->GetN()-1);

    TGraphErrors* graphPythia8InvXSection900GeV               = new TGraphErrors((TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013Eta900GeV"));
    while(graphPythia8InvXSection900GeV->GetX()[0] < 0.8) graphPythia8InvXSection900GeV->RemovePoint(0);
    while(graphPythia8InvXSection900GeV->GetX()[graphPythia8InvXSection900GeV->GetN()-1] > 3.) graphPythia8InvXSection900GeV->RemovePoint(graphPythia8InvXSection900GeV->GetN()-1);


    // *******************************************************************************************************
    // NLO calc
    TGraph* graphNLOCalcEtaMuHalf8TeV                       = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuHalf8000GeV");
    TGraph* graphNLOCalcEtaMuOne8TeV                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuOne8000GeV");
    TGraph* graphNLOCalcEtaMuTwo8TeV                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuTwo8000GeV");
    for (int j=0;j<graphNLOCalcEtaMuHalf8TeV->GetN();j++){
        graphNLOCalcEtaMuHalf8TeV->GetY()[j] *= 1000;
        graphNLOCalcEtaMuOne8TeV->GetY()[j] *= 1000;
        graphNLOCalcEtaMuTwo8TeV->GetY()[j] *= 1000;
    }
    while (graphNLOCalcEtaMuHalf8TeV->GetX()[graphNLOCalcEtaMuHalf8TeV->GetN()-1] > 34. )
        graphNLOCalcEtaMuHalf8TeV->RemovePoint(graphNLOCalcEtaMuHalf8TeV->GetN()-1);
    while (graphNLOCalcEtaMuOne8TeV->GetX()[graphNLOCalcEtaMuOne8TeV->GetN()-1] > 34. )
        graphNLOCalcEtaMuOne8TeV->RemovePoint(graphNLOCalcEtaMuOne8TeV->GetN()-1);
    while (graphNLOCalcEtaMuTwo8TeV->GetX()[graphNLOCalcEtaMuTwo8TeV->GetN()-1] > 34. )
        graphNLOCalcEtaMuTwo8TeV->RemovePoint(graphNLOCalcEtaMuTwo8TeV->GetN()-1);

    TGraph* graphNLOCalcEtaMuHalf7TeV                       = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuHalf7000GeV");
    TGraph* graphNLOCalcEtaMuOne7TeV                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuOne7000GeV");
    TGraph* graphNLOCalcEtaMuTwo7TeV                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuTwo7000GeV");
    for (int j=0;j<graphNLOCalcEtaMuHalf7TeV->GetN();j++){
        graphNLOCalcEtaMuHalf7TeV->GetY()[j] *= 100;
        graphNLOCalcEtaMuOne7TeV->GetY()[j] *= 100;
        graphNLOCalcEtaMuTwo7TeV->GetY()[j] *= 100;
    }
    while (graphNLOCalcEtaMuHalf7TeV->GetX()[graphNLOCalcEtaMuHalf7TeV->GetN()-1] > 20. )
        graphNLOCalcEtaMuHalf7TeV->RemovePoint(graphNLOCalcEtaMuHalf7TeV->GetN()-1);
    while (graphNLOCalcEtaMuOne7TeV->GetX()[graphNLOCalcEtaMuOne7TeV->GetN()-1] > 20. )
        graphNLOCalcEtaMuOne7TeV->RemovePoint(graphNLOCalcEtaMuOne7TeV->GetN()-1);
    while (graphNLOCalcEtaMuTwo7TeV->GetX()[graphNLOCalcEtaMuTwo7TeV->GetN()-1] > 20. )
        graphNLOCalcEtaMuTwo7TeV->RemovePoint(graphNLOCalcEtaMuTwo7TeV->GetN()-1);

    TGraph* graphNLOCalcEtaMuHalf900GeV                       = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuHalf900GeV");
    TGraph* graphNLOCalcEtaMuOne900GeV                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuOne900GeV");
    TGraph* graphNLOCalcEtaMuTwo900GeV                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuTwo900GeV");
    while (graphNLOCalcEtaMuHalf900GeV->GetX()[graphNLOCalcEtaMuHalf900GeV->GetN()-1] > 3.5 )
        graphNLOCalcEtaMuHalf900GeV->RemovePoint(graphNLOCalcEtaMuHalf900GeV->GetN()-1);
    while (graphNLOCalcEtaMuOne900GeV->GetX()[graphNLOCalcEtaMuOne900GeV->GetN()-1] > 3.5 )
        graphNLOCalcEtaMuOne900GeV->RemovePoint(graphNLOCalcEtaMuOne900GeV->GetN()-1);
    while (graphNLOCalcEtaMuTwo900GeV->GetX()[graphNLOCalcEtaMuTwo900GeV->GetN()-1] > 3.5 )
        graphNLOCalcEtaMuTwo900GeV->RemovePoint(graphNLOCalcEtaMuTwo900GeV->GetN()-1);



  
    Double_t arrayBoundariesX1_XSec[2];
    Double_t arrayBoundariesY1_XSec[10];
    Double_t relativeMarginsXXSec[3];
    Double_t relativeMarginsYXSec[3];
    textSizeLabelsPixel = 48;
    ReturnCorrectValuesForCanvasScaling(1250,2500, 1, 9,0.135, 0.005, 0.003,0.04,arrayBoundariesX1_XSec,arrayBoundariesY1_XSec,relativeMarginsXXSec,relativeMarginsYXSec);

    TCanvas* canvasInvSectionPaper      = new TCanvas("canvasInvSectionPaper","",0,0,1250,2500);  // gives the page size
    DrawGammaCanvasSettings( canvasInvSectionPaper,  0.13, 0.02, 0.03, 0.06);

    TPad* padInvSectionSpec             = new TPad("padInvSectionSpec", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[5], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[0],-1, -1, -2);
    DrawGammaPadSettings( padInvSectionSpec, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[0], relativeMarginsYXSec[1]);
    padInvSectionSpec->Draw();
    Double_t marginXSec                 = relativeMarginsXXSec[0]*1250;
    Double_t textsizeLabelsXSecUp       = 0;
    Double_t textsizeFacXSecUp          = 0;
    if (padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) < padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1())){
        textsizeLabelsXSecUp            = (Double_t)textSizeLabelsPixel/padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) ;
        textsizeFacXSecUp               = (Double_t)1./padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) ;
    } else {
        textsizeLabelsXSecUp            = (Double_t)textSizeLabelsPixel/padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1());
        textsizeFacXSecUp               = (Double_t)1./padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1());
    }

    TPad* padInvSectionNLORatio         = new TPad("padInvSectionNLORatio", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[6], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[5],-1, -1, -2);
    DrawGammaPadSettings( padInvSectionNLORatio, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[1]);
    padInvSectionNLORatio->Draw();
    Double_t textsizeLabelsXSecMiddle   = 0;
    Double_t textsizeFacXSecMiddle      = 0;
    if (padInvSectionNLORatio->XtoPixel(padInvSectionNLORatio->GetX2()) < padInvSectionNLORatio->YtoPixel(padInvSectionNLORatio->GetY1())){
        textsizeLabelsXSecMiddle        = (Double_t)textSizeLabelsPixel/padInvSectionNLORatio->XtoPixel(padInvSectionNLORatio->GetX2()) ;
        textsizeFacXSecMiddle           = (Double_t)1./padInvSectionNLORatio->XtoPixel(padInvSectionNLORatio->GetX2()) ;
    } else {
        textsizeLabelsXSecMiddle        = (Double_t)textSizeLabelsPixel/padInvSectionNLORatio->YtoPixel(padInvSectionNLORatio->GetY1());
        textsizeFacXSecMiddle           = (Double_t)1./padInvSectionNLORatio->YtoPixel(padInvSectionNLORatio->GetY1());
    }

    TPad* padInvSectionPythiaRatio      = new TPad("padInvSectionPythiaRatio", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[7], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[6],-1, -1, -2);
    DrawGammaPadSettings( padInvSectionPythiaRatio, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[1]);
    padInvSectionPythiaRatio->Draw();
    Double_t textsizeLabelsXSecDown     = 0;
    Double_t textsizeFacXSecDown        = 0;
    if (padInvSectionPythiaRatio->XtoPixel(padInvSectionPythiaRatio->GetX2()) < padInvSectionPythiaRatio->YtoPixel(padInvSectionPythiaRatio->GetY1())){
        textsizeLabelsXSecDown          = (Double_t)textSizeLabelsPixel/padInvSectionPythiaRatio->XtoPixel(padInvSectionPythiaRatio->GetX2()) ;
        textsizeFacXSecDown             = (Double_t)1./padInvSectionPythiaRatio->XtoPixel(padInvSectionPythiaRatio->GetX2()) ;
    } else {
        textsizeLabelsXSecDown          = (Double_t)textSizeLabelsPixel/padInvSectionPythiaRatio->YtoPixel(padInvSectionPythiaRatio->GetY1());
        textsizeFacXSecDown             = (Double_t)1./padInvSectionPythiaRatio->YtoPixel(padInvSectionPythiaRatio->GetY1());
    }
    TPad* padInvSection2760GeVRatio      = new TPad("padInvSection2760GeVRatio", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[8], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[7],-1, -1, -2);
    DrawGammaPadSettings( padInvSection2760GeVRatio, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[1]);
    padInvSection2760GeVRatio->Draw();
    Double_t textsizeLabelsXSecDown26     = 0;
    Double_t textsizeFacXSecDown26        = 0;
    if (padInvSection2760GeVRatio->XtoPixel(padInvSection2760GeVRatio->GetX2()) < padInvSection2760GeVRatio->YtoPixel(padInvSection2760GeVRatio->GetY1())){
        textsizeLabelsXSecDown26          = (Double_t)textSizeLabelsPixel/padInvSection2760GeVRatio->XtoPixel(padInvSection2760GeVRatio->GetX2()) ;
        textsizeFacXSecDown26             = (Double_t)1./padInvSection2760GeVRatio->XtoPixel(padInvSection2760GeVRatio->GetX2()) ;
    } else {
        textsizeLabelsXSecDown26          = (Double_t)textSizeLabelsPixel/padInvSection2760GeVRatio->YtoPixel(padInvSection2760GeVRatio->GetY1());
        textsizeFacXSecDown26             = (Double_t)1./padInvSection2760GeVRatio->YtoPixel(padInvSection2760GeVRatio->GetY1());
    }
    TPad* padInvSection900GeVRatio      = new TPad("padInvSection900GeVRatio", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[9], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[8],-1, -1, -2);
    DrawGammaPadSettings( padInvSection900GeVRatio, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[2]);
    padInvSection900GeVRatio->Draw();
    Double_t textsizeLabelsXSecDown2     = 0;
    Double_t textsizeFacXSecDown2        = 0;
    if (padInvSection900GeVRatio->XtoPixel(padInvSection900GeVRatio->GetX2()) < padInvSection900GeVRatio->YtoPixel(padInvSection900GeVRatio->GetY1())){
        textsizeLabelsXSecDown2          = (Double_t)textSizeLabelsPixel/padInvSection900GeVRatio->XtoPixel(padInvSection900GeVRatio->GetX2()) ;
        textsizeFacXSecDown2             = (Double_t)1./padInvSection900GeVRatio->XtoPixel(padInvSection900GeVRatio->GetX2()) ;
    } else {
        textsizeLabelsXSecDown2          = (Double_t)textSizeLabelsPixel/padInvSection900GeVRatio->YtoPixel(padInvSection900GeVRatio->GetY1());
        textsizeFacXSecDown2             = (Double_t)1./padInvSection900GeVRatio->YtoPixel(padInvSection900GeVRatio->GetY1());
    }
    
    Double_t maxX = 49;

    padInvSectionSpec->cd();
    padInvSectionSpec->SetLogy(1);
    padInvSectionSpec->SetLogx(1);

    TH2F * histo2DXSectionPi0;
    //     histo2DXSectionPi0          = new TH2F("histo2DXSectionPi0","histo2DXSectionPi0",11000,0.23,50.,1000,6,9e11);
    histo2DXSectionPi0          = new TH2F("histo2DXSectionPi0","histo2DXSectionPi0",11000,0.33,maxX,1000,11,2.3e13);
    SetStyleHistoTH2ForGraphs(histo2DXSectionPi0, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",0.035,0.04, 0.035,0.04, 0.9,1.45);
    histo2DXSectionPi0->GetXaxis()->SetMoreLogLabels();
    histo2DXSectionPi0->GetXaxis()->SetNoExponent(kTRUE);

    SetStyleHistoTH2ForGraphs(histo2DXSectionPi0, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",
                            0.85*textsizeLabelsXSecUp,textsizeLabelsXSecUp, 0.85*textsizeLabelsXSecUp, textsizeLabelsXSecUp, 1,0.2/(textsizeFacXSecUp*marginXSec));
    histo2DXSectionPi0->GetXaxis()->SetMoreLogLabels();
    histo2DXSectionPi0->GetXaxis()->SetLabelOffset(+0.01);
    histo2DXSectionPi0->Draw();
    if(plotTheorycurves){
      DrawGammaNLOTGraph( graphNLOCalcEtaMuHalf8TeV, widthCommonFit, styleLineNLOMuHalf, colorNLO);
      graphNLOCalcEtaMuHalf8TeV->Draw("same,c");
      DrawGammaNLOTGraph( graphNLOCalcEtaMuOne8TeV, widthCommonFit, styleLineNLOMuOne, colorNLO);
      graphNLOCalcEtaMuOne8TeV->Draw("same,c");
      DrawGammaNLOTGraph( graphNLOCalcEtaMuTwo8TeV, widthCommonFit, styleLineNLOMuTwo, colorNLO);
      graphNLOCalcEtaMuTwo8TeV->Draw("same,c");

      DrawGammaNLOTGraph( graphNLOCalcEtaMuHalf7TeV, widthCommonFit, styleLineNLOMuHalf, colorNLO);
      graphNLOCalcEtaMuHalf7TeV->Draw("same,c");
      DrawGammaNLOTGraph( graphNLOCalcEtaMuOne7TeV, widthCommonFit, styleLineNLOMuOne, colorNLO);
      graphNLOCalcEtaMuOne7TeV->Draw("same,c");
      DrawGammaNLOTGraph( graphNLOCalcEtaMuTwo7TeV, widthCommonFit, styleLineNLOMuTwo, colorNLO);
      graphNLOCalcEtaMuTwo7TeV->Draw("same,c");

      DrawGammaNLOTGraph( graphNLOCalcEtaMuHalf900GeV, widthCommonFit, styleLineNLOMuHalf, colorNLO);
      graphNLOCalcEtaMuHalf900GeV->Draw("same,c");
      DrawGammaNLOTGraph( graphNLOCalcEtaMuOne900GeV, widthCommonFit, styleLineNLOMuOne, colorNLO);
      graphNLOCalcEtaMuOne900GeV->Draw("same,c");
      DrawGammaNLOTGraph( graphNLOCalcEtaMuTwo900GeV, widthCommonFit, styleLineNLOMuTwo, colorNLO);
      graphNLOCalcEtaMuTwo900GeV->Draw("same,c");
    }

    Color_t pythia8color = kRed+2;
    if(plotTheorycurves){
      DrawGammaSetMarkerTGraphErr(graphPythia8InvXSection8TeV, 0, 0, pythia8color , pythia8color, widthLinesBoxes, kTRUE, pythia8color);
      graphPythia8InvXSection8TeV->Draw("3,same");
      DrawGammaSetMarkerTGraphErr(graphPythia8InvXSection7TeV, 0, 0, pythia8color , pythia8color, widthLinesBoxes, kTRUE, pythia8color);
      graphPythia8InvXSection7TeV->Draw("3,same");
      DrawGammaSetMarkerTGraphErr(graphPythia8InvXSection900GeV, 0, 0, pythia8color , pythia8color, widthLinesBoxes, kTRUE, pythia8color);
      graphPythia8InvXSection900GeV->Draw("3,same");
      DrawGammaSetMarker(histoPythia8InvXSectionEta8TeV, 24, 1.5, pythia8color , pythia8color);
      histoPythia8InvXSectionEta8TeV->SetLineWidth(widthCommonFit);
      histoPythia8InvXSectionEta8TeV->Draw("same,hist,l");
      DrawGammaSetMarker(histoPythia8InvXSectionEta7TeV, 24, 1.5, pythia8color , pythia8color);
      histoPythia8InvXSectionEta7TeV->SetLineWidth(widthCommonFit);
      histoPythia8InvXSectionEta7TeV->Draw("same,hist,l");
      DrawGammaSetMarker(histoPythia8InvXSectionEta900GeV, 24, 1.5, pythia8color , pythia8color);
      histoPythia8InvXSectionEta900GeV->SetLineWidth(widthCommonFit);
      histoPythia8InvXSectionEta900GeV->Draw("same,hist,l");
    }
    Double_t paramTCMPi0New8TeV[5]  = { csGraphs[2]->GetY()[1],0.1, csGraphs[2]->GetY()[4],0.6,3.0};
    TF1* fitTCMInvXSectionPi08TeV   = FitObject("tcm","fitTCMInvCrossSectionEta8TeV","Pi0",csGraphs[2],0.4,35. ,paramTCMPi0New8TeV,"QNRMEX0+","", kFALSE);

    TF1* fitTCMInvXSectionPi0Plot8TeV = new TF1("twoCompModel_plottingEta8TeV",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mesonMassExpectPi0,mesonMassExpectPi0,mesonMassExpectPi0));
    fitTCMInvXSectionPi0Plot8TeV->SetRange(0.4,35.);
    fitTCMInvXSectionPi0Plot8TeV->SetParameters(fitTCMInvXSectionPi08TeV->GetParameters());
    fitTCMInvXSectionPi0Plot8TeV->SetParErrors(fitTCMInvXSectionPi08TeV->GetParErrors());
    cout << WriteParameterToFile(fitTCMInvXSectionPi08TeV) << endl;
    DrawGammaSetMarkerTF1( fitTCMInvXSectionPi0Plot8TeV, 7, 2, kGray+2);
    // fitTCMInvXSectionPi0Plot8TeV->Draw("same");
    
    Double_t paramTCMPi0New7TeV[5]  = { csGraphs[1]->GetY()[1],0.1, csGraphs[1]->GetY()[4],0.6,3.0};
    TF1* fitTCMInvXSectionPi07TeV   = FitObject("tcm","fitTCMInvCrossSectionEta7TeV","Pi0",csGraphs[1],0.4,25. ,paramTCMPi0New7TeV,"QNRMEX0+","", kFALSE);
    
    TF1* fitTCMInvXSectionPi0Plot7TeV = new TF1("twoCompModel_plottingEta7TeV",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mesonMassExpectPi0,mesonMassExpectPi0,mesonMassExpectPi0));
    fitTCMInvXSectionPi0Plot7TeV->SetRange(0.4,35.);
    fitTCMInvXSectionPi0Plot7TeV->SetParameters(fitTCMInvXSectionPi07TeV->GetParameters());
    fitTCMInvXSectionPi0Plot7TeV->SetParErrors(fitTCMInvXSectionPi07TeV->GetParErrors());
    cout << WriteParameterToFile(fitTCMInvXSectionPi07TeV) << endl;
    DrawGammaSetMarkerTF1( fitTCMInvXSectionPi0Plot7TeV, 7, 2, kGray+2);
    fitTCMInvXSectionPi0Plot7TeV->Draw("same");
    
    Double_t paramTCMPi0New276TeV[5]  = { csGraphs[3]->GetY()[1],0.1, csGraphs[3]->GetY()[4],0.6,3.0};
    TF1* fitTCMInvXSectionPi0276TeV   = FitObject("tcm","fitTCMInvCrossSectionEta7TeV","Pi0",csGraphs[3],0.5,20. ,paramTCMPi0New276TeV,"QNRMEX0+","", kFALSE);
    
    TF1* fitTCMInvXSectionPi0Plot276TeV = new TF1("twoCompModel_plottingEta7TeV",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mesonMassExpectPi0,mesonMassExpectPi0,mesonMassExpectPi0));
    fitTCMInvXSectionPi0Plot276TeV->SetRange(0.5,21.);
    fitTCMInvXSectionPi0Plot276TeV->SetParameters(fitTCMInvXSectionPi0276TeV->GetParameters());
    fitTCMInvXSectionPi0Plot276TeV->SetParErrors(fitTCMInvXSectionPi0276TeV->GetParErrors());
    cout << WriteParameterToFile(fitTCMInvXSectionPi0276TeV) << endl;
    DrawGammaSetMarkerTF1( fitTCMInvXSectionPi0Plot276TeV, 7, 2, kGray+2);
    fitTCMInvXSectionPi0Plot276TeV->Draw("same");

    Double_t paramTCMPi0New900GeV[5]  = { csGraphs[0]->GetY()[1],0.1, csGraphs[0]->GetY()[4],0.6,3.0};
    TF1* fitTCMInvXSectionPi0900GeV   = FitObject("tcm","fitTCMInvCrossSectionEta900GeV","Pi0",csGraphs[0],0.8,3. ,paramTCMPi0New900GeV,"QNRMEX0+","", kFALSE);
    TF1* fitTCMInvXSectionPi0Plot900GeV = new TF1("twoCompModel_plotting900GeV",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mesonMassExpectPi0,mesonMassExpectPi0,mesonMassExpectPi0));
    fitTCMInvXSectionPi0Plot900GeV->SetRange(0.8,3.);
    fitTCMInvXSectionPi0Plot900GeV->SetParameters(fitTCMInvXSectionPi0900GeV->GetParameters());
    fitTCMInvXSectionPi0Plot900GeV->SetParErrors(fitTCMInvXSectionPi0900GeV->GetParErrors());
    //     cout << WriteParameterToFile(fitTCMInvXSectionPi0900GeV) << endl;
    DrawGammaSetMarkerTF1( fitTCMInvXSectionPi0Plot900GeV, 7, 2, kGray+2);
    //         fitTCMInvXSectionPi0Plot900GeV->Draw("same");


    // Tsallis fit
    Double_t paramGraphPi08[3]                       = {5e11, 6., 0.13};
    TF1* fitInvXSectionPi08TeV                      = FitObject("l","fitInvCrossSectionEta8TeV","Pi0",csGraphs[2],0.4,35.,paramGraphPi08,"QNRMEX0+");
    DrawGammaSetMarkerTF1( fitInvXSectionPi08TeV, 3, 2, kGray+1);
    fitInvXSectionPi08TeV->Draw("same");
    cout << WriteParameterToFile(fitInvXSectionPi08TeV) << endl;
    Double_t paramGraphPi07[3]                                   = {5e11, 6., 0.13};
    TF1* fitInvXSectionPi07TeV                      = FitObject("l","fitInvCrossSectionEta7TeV","Pi0",csGraphs[1],0.4,35.,paramGraphPi07,"QNRMEX0+");
    DrawGammaSetMarkerTF1( fitInvXSectionPi07TeV, 3, 2, kGray+1);
    fitInvXSectionPi07TeV->Draw("same");
    cout << WriteParameterToFile(fitInvXSectionPi07TeV) << endl;
    Double_t paramGraphPi0276[3]                                   = {5e11, 6., 0.13};
    TF1* fitInvXSectionPi0276TeV                      = FitObject("l","fitInvCrossSectionEta7TeV","Pi0",csGraphs[3],0.5,21.,paramGraphPi0276,"QNRMEX0+");
    DrawGammaSetMarkerTF1( fitInvXSectionPi0276TeV, 3, 2, kGray+1);
    fitInvXSectionPi0276TeV->Draw("same");
    cout << WriteParameterToFile(fitInvXSectionPi07TeV) << endl;
    Double_t paramGraphPi0900[3]                                   = {5e11, 6., 0.13};
    TF1* fitInvXSectionPi0900GeV                    = FitObject("l","fitInvCrossSectionEta900GeV","Pi0",csGraphs[0],0.8,3,paramGraphPi0900,"QNRMEX0+");
    DrawGammaSetMarkerTF1( fitInvXSectionPi0900GeV, 3, 2, kGray+1);
    fitInvXSectionPi0900GeV->Draw("same");
    cout << WriteParameterToFile(fitInvXSectionPi0900GeV) << endl;

    csGraphsSys[2]     ->Draw("E2same");
    csGraphsSys[1]     ->Draw("E2same");
    csGraphsSys[3]     ->Draw("E2same");
    csGraphsSys[0]     ->Draw("E2same");
    csGraphs[2]->Draw("p,same,z");
    csGraphs[1]->Draw("p,same,z");
    csGraphs[3]->Draw("p,same,z");
    csGraphs[0]->Draw("p,same,z");


    Double_t rightalignDouble = 0.93;
    drawLatexAdd("ALICE",rightalignDouble,0.91,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#eta #rightarrow #gamma#gamma",rightalignDouble,0.86,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE);

    TBox* boxErrorSigmaRatio = CreateBoxConv(kGray+2, 0.3, 1.-(0.0196 ), 0.35, 1.+(0.0196));
    boxErrorSigmaRatio->SetLineWidth(8);


    TLegend* legendXsectionPaper;
    if(plotTheorycurves)
      legendXsectionPaper    = GetAndSetLegend2(0.17, 0.074, 0.5, 0.144+0.05*5, textSizeLabelsPixel);
      else
      legendXsectionPaper    = GetAndSetLegend2(0.17, 0.03, 0.5, 0.104+0.05*4, textSizeLabelsPixel);
    legendXsectionPaper->SetNColumns(1);
    legendXsectionPaper->SetMargin(0.2);
    legendXsectionPaper->AddEntry(csGraphsSys[2],Form("%s (x 10^{3})",nameMeasGlobal[2].Data()),"pf");
    legendXsectionPaper->AddEntry(csGraphsSys[1],Form("%s (x 10^{2})",nameMeasGlobal[1].Data()),"pf");
    legendXsectionPaper->AddEntry(csGraphsSys[3],Form("%s (x 10)",nameMeasGlobal[3].Data()),"pf");
    legendXsectionPaper->AddEntry(csGraphsSys[0],nameMeasGlobal[0].Data(),"pf");
    legendXsectionPaper->AddEntry(fitTCMInvXSectionPi0Plot8TeV,"TCM fit","l");
    legendXsectionPaper->AddEntry(fitInvXSectionPi08TeV,"Tsallis fit","l");
    if(plotTheorycurves){
      legendXsectionPaper->AddEntry(histoPythia8InvXSectionEta8TeV,"PYTHIA 8.2, Monash 2013","l");
      legendXsectionPaper->AddEntry((TObject*)0, "NLO, PDF:CTEQ6M5 - FF:AESSS", "");
      TLegend* legendXsectionPaperNLO    = GetAndSetLegend2(0.23, 0.03, 0.83, 0.08, textSizeLabelsPixel);
      legendXsectionPaperNLO->SetNColumns(3);
      legendXsectionPaperNLO->AddEntry(graphNLOCalcEtaMuHalf8TeV, "#mu = 0.5 #it{p}_{T}", "l");
      legendXsectionPaperNLO->AddEntry(graphNLOCalcEtaMuOne8TeV,  "#mu = #it{p}_{T}", "l");
      legendXsectionPaperNLO->AddEntry(graphNLOCalcEtaMuTwo8TeV,  "#mu = 2 #it{p}_{T}", "l");
      legendXsectionPaperNLO->Draw();
    }
    legendXsectionPaper->Draw();

    //     drawLatexAdd("x10^{2}",0.24,0.905,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE,colorDataFULL[2]);
    //     drawLatexAdd("x10^{1}",0.24,0.78,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE,colorDataFULL[1]);

    Double_t yrangeHighRatio = 2.98;
    Double_t yrangeLowRatio = 0.41;

    padInvSectionNLORatio->cd();
    padInvSectionNLORatio->SetLogx(1);
    TH2F * ratio8TeVdummy               = new TH2F("ratio8TeVdummy","ratio8TeVdummy",1000,0.33,maxX,1000,-10.01,5);
    SetStyleHistoTH2ForGraphs(ratio8TeVdummy, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{Tsallis fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle,
                            0.85*textsizeLabelsXSecMiddle,0.9*textsizeLabelsXSecMiddle, 1,0.2/(textsizeFacXSecMiddle*marginXSec), 510, 505);
    ratio8TeVdummy->GetYaxis()->SetMoreLogLabels(kTRUE);
    ratio8TeVdummy->GetYaxis()->SetNdivisions(505);
    ratio8TeVdummy->GetYaxis()->SetNoExponent(kTRUE);
    ratio8TeVdummy->GetXaxis()->SetMoreLogLabels(kTRUE);
    ratio8TeVdummy->GetXaxis()->SetNoExponent(kTRUE);
    ratio8TeVdummy->GetXaxis()->SetLabelFont(42);
    ratio8TeVdummy->GetYaxis()->SetLabelFont(42);
    ratio8TeVdummy->GetYaxis()->CenterTitle(kTRUE);
    ratio8TeVdummy->GetYaxis()->SetLabelOffset(+0.01);
    ratio8TeVdummy->GetYaxis()->SetRangeUser(yrangeLowRatio,yrangeHighRatio);
    ratio8TeVdummy->GetXaxis()->SetTickLength(0.07);
    ratio8TeVdummy->DrawCopy();

    TGraph* graphRatioEtaCombNLOMuHalf8TeV                  = (TGraph*)graphNLOCalcEtaMuHalf8TeV->Clone();
    TGraph* graphRatioEtaCombNLOMuOne8TeV                   = (TGraph*)graphNLOCalcEtaMuOne8TeV->Clone();
    TGraph* graphRatioEtaCombNLOMuTwo8TeV                   = (TGraph*)graphNLOCalcEtaMuTwo8TeV->Clone();
    graphRatioEtaCombNLOMuHalf8TeV                          = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuHalf8TeV, fitInvXSectionPi08TeV);
    graphRatioEtaCombNLOMuOne8TeV                           = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuOne8TeV, fitInvXSectionPi08TeV);
    graphRatioEtaCombNLOMuTwo8TeV                           = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuTwo8TeV, fitInvXSectionPi08TeV);
    if(plotTheorycurves){
      DrawGammaNLOTGraph( graphRatioEtaCombNLOMuHalf8TeV, widthCommonFit, styleLineNLOMuHalf, colorNLO);
      graphRatioEtaCombNLOMuHalf8TeV->Draw("same,c");
      DrawGammaNLOTGraph( graphRatioEtaCombNLOMuOne8TeV, widthCommonFit, styleLineNLOMuOne, colorNLO);
      graphRatioEtaCombNLOMuOne8TeV->Draw("same,c");
      DrawGammaNLOTGraph( graphRatioEtaCombNLOMuTwo8TeV, widthCommonFit, styleLineNLOMuTwo, colorNLO);
      graphRatioEtaCombNLOMuTwo8TeV->Draw("same,c");
    }
    TH1D* histoRatioPythia8ToFit8TeV                     = (TH1D*) histoPythia8InvXSectionEta8TeV->Clone();
    histoRatioPythia8ToFit8TeV                           = CalculateHistoRatioToFit (histoRatioPythia8ToFit8TeV, fitInvXSectionPi08TeV);
    DrawGammaSetMarker(histoRatioPythia8ToFit8TeV, 24, 1.5, pythia8color , pythia8color);
    histoRatioPythia8ToFit8TeV->SetLineWidth(widthCommonFit);
    if(plotTheorycurves)
    histoRatioPythia8ToFit8TeV->Draw("same,hist,l");


    TGraphErrors* graphRatioPythia8ToFit8TeV             = (TGraphErrors*) graphPythia8InvXSection8TeV->Clone();
    graphRatioPythia8ToFit8TeV                           = CalculateGraphErrRatioToFit (graphRatioPythia8ToFit8TeV, fitInvXSectionPi08TeV);
    DrawGammaSetMarkerTGraphErr(graphRatioPythia8ToFit8TeV, 0, 0, pythia8color , pythia8color, widthLinesBoxes, kTRUE, pythia8color);
    if(plotTheorycurves)
    graphRatioPythia8ToFit8TeV->Draw("3,same");


    TGraphAsymmErrors* graphRatioPi0CombCombFitStatA    = (TGraphAsymmErrors*)csGraphs[2]->Clone();
    graphRatioPi0CombCombFitStatA                       = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitStatA, fitInvXSectionPi08TeV);
    TGraphAsymmErrors* graphRatioPi0CombCombFitStatA_WOXErr = (TGraphAsymmErrors*) graphRatioPi0CombCombFitStatA->Clone("graphRatioPi0CombCombFitStatA_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombCombFitStatA_WOXErr);

    TGraphAsymmErrors* graphRatioPi0CombCombFitSysA     = (TGraphAsymmErrors*)csGraphsSys[2]->Clone();
    graphRatioPi0CombCombFitSysA                        = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitSysA, fitInvXSectionPi08TeV);
    DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA,markerstylesFULL[2], markersizeFULL[2], colorDataFULL[2], colorDataFULL[2], widthLinesBoxes, kTRUE, 0);
    graphRatioPi0CombCombFitSysA->SetLineWidth(0);
    graphRatioPi0CombCombFitSysA->Draw("2,same");
    DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA_WOXErr, markerstylesFULL[2], markersizeFULL[2], colorDataFULL[2], colorDataFULL[2], widthLinesBoxes, kFALSE);
    graphRatioPi0CombCombFitStatA_WOXErr->SetLineWidth(widthLinesBoxes);
    graphRatioPi0CombCombFitStatA_WOXErr->Draw("p,same");
    //
    //         boxErrorSigmaRatio->Draw();
    DrawGammaLines(0.33, maxX , 1., 1.,0.5, kGray+2);

    padInvSectionPythiaRatio->cd();
    padInvSectionPythiaRatio->SetLogx(1);
    TH2F * ratio7TeVdummy               = new TH2F("ratio7TeVdummy","ratio7TeVdummy",1000,0.33,maxX,1000,-1.01,5);
    SetStyleHistoTH2ForGraphs(ratio7TeVdummy, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{Tsallis fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
                            0.85*textsizeLabelsXSecDown,0.9*textsizeLabelsXSecDown, 1,0.2/(textsizeFacXSecMiddle*marginXSec), 510, 505);
    ratio7TeVdummy->GetYaxis()->SetMoreLogLabels(kTRUE);
    ratio7TeVdummy->GetYaxis()->SetNdivisions(505);
    ratio7TeVdummy->GetYaxis()->SetNoExponent(kTRUE);
    ratio7TeVdummy->GetXaxis()->SetMoreLogLabels(kTRUE);
    ratio7TeVdummy->GetXaxis()->SetNoExponent(kTRUE);
    ratio7TeVdummy->GetXaxis()->SetLabelFont(42);
    ratio7TeVdummy->GetYaxis()->SetLabelFont(42);
    ratio7TeVdummy->GetYaxis()->CenterTitle(kTRUE);
    ratio7TeVdummy->GetYaxis()->SetLabelOffset(+0.01);
    ratio7TeVdummy->GetYaxis()->SetRangeUser(yrangeLowRatio,yrangeHighRatio);
    ratio7TeVdummy->GetXaxis()->SetTickLength(0.07);
    ratio7TeVdummy->DrawCopy();

    TGraph* graphRatioEtaCombNLOMuHalf7TeV                  = (TGraph*)graphNLOCalcEtaMuHalf7TeV->Clone();
    TGraph* graphRatioEtaCombNLOMuOne7TeV                   = (TGraph*)graphNLOCalcEtaMuOne7TeV->Clone();
    TGraph* graphRatioEtaCombNLOMuTwo7TeV                   = (TGraph*)graphNLOCalcEtaMuTwo7TeV->Clone();
    graphRatioEtaCombNLOMuHalf7TeV                          = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuHalf7TeV, fitInvXSectionPi07TeV);
    graphRatioEtaCombNLOMuOne7TeV                           = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuOne7TeV, fitInvXSectionPi07TeV);
    graphRatioEtaCombNLOMuTwo7TeV                           = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuTwo7TeV, fitInvXSectionPi07TeV);
    if(plotTheorycurves){
      DrawGammaNLOTGraph( graphRatioEtaCombNLOMuHalf7TeV, widthCommonFit, styleLineNLOMuHalf, colorNLO);
      graphRatioEtaCombNLOMuHalf7TeV->Draw("same,c");
      DrawGammaNLOTGraph( graphRatioEtaCombNLOMuOne7TeV, widthCommonFit, styleLineNLOMuOne, colorNLO);
      graphRatioEtaCombNLOMuOne7TeV->Draw("same,c");
      DrawGammaNLOTGraph( graphRatioEtaCombNLOMuTwo7TeV, widthCommonFit, styleLineNLOMuTwo, colorNLO);
      graphRatioEtaCombNLOMuTwo7TeV->Draw("same,c");
    }

    TH1D* histoRatioPythia8ToFit7TeV                     = (TH1D*) histoPythia8InvXSectionEta7TeV->Clone();
    histoRatioPythia8ToFit7TeV                           = CalculateHistoRatioToFit (histoRatioPythia8ToFit7TeV, fitInvXSectionPi07TeV);
    DrawGammaSetMarker(histoRatioPythia8ToFit7TeV, 24, 1.5, pythia8color , pythia8color);
    histoRatioPythia8ToFit7TeV->SetLineWidth(widthCommonFit);
    if(plotTheorycurves)
    histoRatioPythia8ToFit7TeV->Draw("same,hist,l");


    TGraphErrors* graphRatioPythia8ToFit7TeV             = (TGraphErrors*) graphPythia8InvXSection7TeV->Clone();
    graphRatioPythia8ToFit7TeV                           = CalculateGraphErrRatioToFit (graphRatioPythia8ToFit7TeV, fitInvXSectionPi07TeV);
    DrawGammaSetMarkerTGraphErr(graphRatioPythia8ToFit7TeV, 0, 0, pythia8color , pythia8color, widthLinesBoxes, kTRUE, pythia8color);
    if(plotTheorycurves)
    graphRatioPythia8ToFit7TeV->Draw("3,same");


    TGraphAsymmErrors* graphRatioPi0CombCombFitStatA7TeV    = (TGraphAsymmErrors*)csGraphs[1]->Clone();
    graphRatioPi0CombCombFitStatA7TeV                       = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitStatA7TeV, fitInvXSectionPi07TeV);
    TGraphAsymmErrors* graphRatioPi0CombCombFitStatA7TeV_WOXErr = (TGraphAsymmErrors*) graphRatioPi0CombCombFitStatA7TeV->Clone("graphRatioPi0CombCombFitStatA7TeV_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombCombFitStatA7TeV_WOXErr);

    TGraphAsymmErrors* graphRatioPi0CombCombFitSysA7TeV     = (TGraphAsymmErrors*)csGraphsSys[1]->Clone();
    graphRatioPi0CombCombFitSysA7TeV                        = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitSysA7TeV, fitInvXSectionPi07TeV);
    DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA7TeV, markerstylesFULL[1], markersizeFULL[1], colorDataFULL[1], colorDataFULL[1], widthLinesBoxes, kTRUE, 0);
    graphRatioPi0CombCombFitSysA7TeV->SetLineWidth(0);
    graphRatioPi0CombCombFitSysA7TeV->Draw("2,same");
    DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA7TeV_WOXErr, markerstylesFULL[1], markersizeFULL[1], colorDataFULL[1], colorDataFULL[1], widthLinesBoxes, kFALSE);
    graphRatioPi0CombCombFitStatA7TeV_WOXErr->SetLineWidth(widthLinesBoxes);
    graphRatioPi0CombCombFitStatA7TeV_WOXErr->Draw("p,same");


    DrawGammaLines(0.33, maxX , 1., 1.,0.5, kGray+2);
    
    padInvSection2760GeVRatio->cd();
    padInvSection2760GeVRatio->SetLogx(1);
    TH2F * ratio276TeVdummy               = new TH2F("ratio276TeVdummy","ratio276TeVdummy",1000,0.33,maxX,1000,yrangeLowRatio,yrangeHighRatio);
    SetStyleHistoTH2ForGraphs(ratio276TeVdummy, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{Tsallis fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
                              0.85*textsizeLabelsXSecDown,0.9*textsizeLabelsXSecDown, 1,0.2/(textsizeFacXSecMiddle*marginXSec), 510, 505);
    ratio276TeVdummy->GetYaxis()->SetMoreLogLabels(kTRUE);
    ratio276TeVdummy->GetYaxis()->SetNdivisions(505);
    ratio276TeVdummy->GetYaxis()->SetNoExponent(kTRUE);
    ratio276TeVdummy->GetXaxis()->SetMoreLogLabels(kTRUE);
    ratio276TeVdummy->GetXaxis()->SetNoExponent(kTRUE);
    ratio276TeVdummy->GetXaxis()->SetLabelFont(42);
    ratio276TeVdummy->GetYaxis()->SetLabelFont(42);
    ratio276TeVdummy->GetYaxis()->CenterTitle(kTRUE);
    ratio276TeVdummy->GetYaxis()->SetLabelOffset(+0.01);
    ratio276TeVdummy->GetXaxis()->SetTickLength(0.07);
    ratio276TeVdummy->DrawCopy();
    
    
    TGraphAsymmErrors* graphRatioPi0CombCombFitStatA276TeV    = (TGraphAsymmErrors*)csGraphs[3]->Clone();
    graphRatioPi0CombCombFitStatA276TeV                       = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitStatA276TeV, fitInvXSectionPi0276TeV);
    TGraphAsymmErrors* graphRatioPi0CombCombFitStatA276TeV_WOXErr = (TGraphAsymmErrors*) graphRatioPi0CombCombFitStatA276TeV->Clone("graphRatioPi0CombCombFitStatA276TeV_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombCombFitStatA276TeV_WOXErr);

    TGraphAsymmErrors* graphRatioPi0CombCombFitSysA276TeV     = (TGraphAsymmErrors*)csGraphsSys[3]->Clone();
    graphRatioPi0CombCombFitSysA276TeV                        = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitSysA276TeV, fitInvXSectionPi0276TeV);
    DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA276TeV, markerstylesFULL[3], markersizeFULL[3], colorDataFULL[3], colorDataFULL[3], widthLinesBoxes, kTRUE, 0);
    graphRatioPi0CombCombFitSysA276TeV->SetLineWidth(0);
    graphRatioPi0CombCombFitSysA276TeV->Draw("2,same");
    DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA276TeV_WOXErr, markerstylesFULL[3], markersizeFULL[3], colorDataFULL[3], colorDataFULL[3], widthLinesBoxes, kFALSE);
    graphRatioPi0CombCombFitStatA276TeV_WOXErr->SetLineWidth(widthLinesBoxes);
    graphRatioPi0CombCombFitStatA276TeV_WOXErr->Draw("p,same");
    
    DrawGammaLines(0.33, maxX , 1., 1.,0.5, kGray+2);
    
    padInvSection900GeVRatio->cd();
    padInvSection900GeVRatio->SetLogx(1);
    TH2F * ratio900GeVdummy            = new TH2F("ratio900GeVdummy","ratio900GeVdummy",1000,0.33,maxX,1000,-1.01,5.);
    SetStyleHistoTH2ForGraphs(ratio900GeVdummy, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{Tsallis fit}", 0.85*textsizeLabelsXSecDown2, textsizeLabelsXSecDown2,
                            0.85*textsizeLabelsXSecDown2,0.9*textsizeLabelsXSecDown2, 0.9,0.2/(textsizeFacXSecDown2*marginXSec), 510, 505);
    ratio900GeVdummy->GetYaxis()->SetMoreLogLabels(kTRUE);
    ratio900GeVdummy->GetYaxis()->SetNdivisions(505);
    ratio900GeVdummy->GetYaxis()->SetNoExponent(kTRUE);
    ratio900GeVdummy->GetXaxis()->SetMoreLogLabels(kTRUE);
    ratio900GeVdummy->GetXaxis()->SetNoExponent(kTRUE);
    ratio900GeVdummy->GetXaxis()->SetLabelFont(42);
    ratio900GeVdummy->GetYaxis()->SetLabelFont(42);
    ratio900GeVdummy->GetYaxis()->CenterTitle(kTRUE);
    ratio900GeVdummy->GetYaxis()->SetLabelOffset(+0.01);
    ratio900GeVdummy->GetYaxis()->SetRangeUser(yrangeLowRatio,yrangeHighRatio);
    ratio900GeVdummy->GetXaxis()->SetTickLength(0.06);
    ratio900GeVdummy->GetYaxis()->SetTickLength(0.04);
    ratio900GeVdummy->DrawCopy();

    TGraph* graphRatioEtaCombNLOMuHalf900GeV                  = (TGraph*)graphNLOCalcEtaMuHalf900GeV->Clone();
    TGraph* graphRatioEtaCombNLOMuOne900GeV                   = (TGraph*)graphNLOCalcEtaMuOne900GeV->Clone();
    TGraph* graphRatioEtaCombNLOMuTwo900GeV                   = (TGraph*)graphNLOCalcEtaMuTwo900GeV->Clone();
    graphRatioEtaCombNLOMuHalf900GeV                          = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuHalf900GeV, fitInvXSectionPi0900GeV);
    graphRatioEtaCombNLOMuOne900GeV                           = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuOne900GeV, fitInvXSectionPi0900GeV);
    graphRatioEtaCombNLOMuTwo900GeV                           = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuTwo900GeV, fitInvXSectionPi0900GeV);
    if(plotTheorycurves){
      DrawGammaNLOTGraph( graphRatioEtaCombNLOMuHalf900GeV, widthCommonFit, styleLineNLOMuHalf, colorNLO);
      graphRatioEtaCombNLOMuHalf900GeV->Draw("same,c");
      DrawGammaNLOTGraph( graphRatioEtaCombNLOMuOne900GeV, widthCommonFit, styleLineNLOMuOne, colorNLO);
      graphRatioEtaCombNLOMuOne900GeV->Draw("same,c");
      DrawGammaNLOTGraph( graphRatioEtaCombNLOMuTwo900GeV, widthCommonFit, styleLineNLOMuTwo, colorNLO);
      graphRatioEtaCombNLOMuTwo900GeV->Draw("same,c");
    }


    TH1D* histoRatioPythia8ToFit900GeV                     = (TH1D*) histoPythia8InvXSectionEta900GeV->Clone();
    histoRatioPythia8ToFit900GeV                           = CalculateHistoRatioToFit (histoRatioPythia8ToFit900GeV, fitInvXSectionPi0900GeV);
    DrawGammaSetMarker(histoRatioPythia8ToFit900GeV, 24, 1.5, pythia8color , pythia8color);
    histoRatioPythia8ToFit900GeV->SetLineWidth(widthCommonFit);
    if(plotTheorycurves)
    histoRatioPythia8ToFit900GeV->Draw("same,hist,l");


    TGraphErrors* graphRatioPythia8ToFit900GeV             = (TGraphErrors*) graphPythia8InvXSection900GeV->Clone();
    graphRatioPythia8ToFit900GeV                           = CalculateGraphErrRatioToFit (graphRatioPythia8ToFit900GeV, fitInvXSectionPi0900GeV);
    DrawGammaSetMarkerTGraphErr(graphRatioPythia8ToFit900GeV, 0, 0, pythia8color , pythia8color, widthLinesBoxes, kTRUE, pythia8color);
    if(plotTheorycurves)
    graphRatioPythia8ToFit900GeV->Draw("3,same");


    TGraphAsymmErrors* graphRatioPi0CombCombFitStatA900GeV    = (TGraphAsymmErrors*)csGraphs[0]->Clone();
    graphRatioPi0CombCombFitStatA900GeV                       = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitStatA900GeV, fitInvXSectionPi0900GeV);
    TGraphAsymmErrors* graphRatioPi0CombCombFitStatA900GeV_WOXErr = (TGraphAsymmErrors*) graphRatioPi0CombCombFitStatA900GeV->Clone("graphRatioPi0CombCombFitStatA900GeV_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombCombFitStatA900GeV_WOXErr);

    TGraphAsymmErrors* graphRatioPi0CombCombFitSysA900GeV     = (TGraphAsymmErrors*)csGraphsSys[0]->Clone();
    graphRatioPi0CombCombFitSysA900GeV                        = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitSysA900GeV, fitInvXSectionPi0900GeV);
    DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA900GeV, markerstylesFULL[0], markersizeFULL[0], colorDataFULL[0], colorDataFULL[0], widthLinesBoxes, kTRUE, 0);
    graphRatioPi0CombCombFitSysA900GeV->SetLineWidth(0);
    graphRatioPi0CombCombFitSysA900GeV->Draw("2,same");
    DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA900GeV_WOXErr, markerstylesFULL[0], markersizeFULL[0], colorDataFULL[0], colorDataFULL[0], widthLinesBoxes, kFALSE);
    graphRatioPi0CombCombFitStatA900GeV_WOXErr->SetLineWidth(widthLinesBoxes);
    graphRatioPi0CombCombFitStatA900GeV_WOXErr->Draw("p,same");

    DrawGammaLines(0.33, maxX , 1., 1.,0.5, kGray+2);

    canvasInvSectionPaper->Print(Form("%s/Eta_InvCrossSectionWithRatios.%s",outputDir.Data(),suffix.Data()));




    // SINGLE PLOTS FOR SLIDES

    Double_t arrayBoundariesX1_XSec2[2];
    Double_t arrayBoundariesY1_XSec2[3];
    Double_t relativeMarginsXXSec2[3];
    Double_t relativeMarginsYXSec2[3];
    textSizeLabelsPixel = 48;
    ReturnCorrectValuesForCanvasScaling(1250,1300, 1, 1,0.135, 0.005, 0.003,0.09,arrayBoundariesX1_XSec2,arrayBoundariesY1_XSec2,relativeMarginsXXSec2,relativeMarginsYXSec2);

    TCanvas* canvasInvSectionPaperSingle      = new TCanvas("canvasInvSectionPaperSingle","",0,0,1250,1250);  // gives the page size
    DrawGammaCanvasSettings( canvasInvSectionPaperSingle,  0.13, 0.02, 0.03, 0.06);

    TPad* padInvSectionSpecSingle             = new TPad("padInvSectionSpecSingle", "", arrayBoundariesX1_XSec2[0], arrayBoundariesY1_XSec2[3], arrayBoundariesX1_XSec2[1], arrayBoundariesY1_XSec2[0],-1, -1, -2);
    DrawGammaPadSettings( padInvSectionSpecSingle, relativeMarginsXXSec2[0], relativeMarginsXXSec2[2], relativeMarginsYXSec2[0], relativeMarginsYXSec2[2]);
    padInvSectionSpecSingle->Draw();
    marginXSec                 = relativeMarginsXXSec2[0]*1250;
    textsizeLabelsXSecUp       = 0;
    textsizeFacXSecUp          = 0;
    if (padInvSectionSpecSingle->XtoPixel(padInvSectionSpecSingle->GetX2()) < padInvSectionSpecSingle->YtoPixel(padInvSectionSpecSingle->GetY1())){
        textsizeLabelsXSecUp            = (Double_t)textSizeLabelsPixel/padInvSectionSpecSingle->XtoPixel(padInvSectionSpecSingle->GetX2()) ;
        textsizeFacXSecUp               = (Double_t)1./padInvSectionSpecSingle->XtoPixel(padInvSectionSpecSingle->GetX2()) ;
    } else {
        textsizeLabelsXSecUp            = (Double_t)textSizeLabelsPixel/padInvSectionSpecSingle->YtoPixel(padInvSectionSpecSingle->GetY1());
        textsizeFacXSecUp               = (Double_t)1./padInvSectionSpecSingle->YtoPixel(padInvSectionSpecSingle->GetY1());
    }

    padInvSectionSpecSingle->cd();
    padInvSectionSpecSingle->SetLogy(1);
    padInvSectionSpecSingle->SetLogx(1);

    histo2DXSectionPi0->Draw();

    graphNLOCalcEtaMuHalf8TeV->Draw("same,c");
    graphNLOCalcEtaMuOne8TeV->Draw("same,c");
    graphNLOCalcEtaMuTwo8TeV->Draw("same,c");

    graphNLOCalcEtaMuHalf7TeV->Draw("same,c");
    graphNLOCalcEtaMuOne7TeV->Draw("same,c");
    graphNLOCalcEtaMuTwo7TeV->Draw("same,c");

    graphNLOCalcEtaMuHalf900GeV->Draw("same,c");
    graphNLOCalcEtaMuOne900GeV->Draw("same,c");
    graphNLOCalcEtaMuTwo900GeV->Draw("same,c");

    graphPythia8InvXSection8TeV->Draw("3,same");
    graphPythia8InvXSection7TeV->Draw("3,same");
    graphPythia8InvXSection900GeV->Draw("3,same");

    histoPythia8InvXSectionEta8TeV->Draw("same,hist,l");
    histoPythia8InvXSectionEta7TeV->Draw("same,hist,l");
    histoPythia8InvXSectionEta900GeV->Draw("same,hist,l");
    // fitTCMInvXSectionPi0Plot8TeV->Draw("same");
    fitTCMInvXSectionPi0Plot7TeV->Draw("same");

    // Tsallis fit
    fitInvXSectionPi08TeV->Draw("same");
    fitInvXSectionPi07TeV->Draw("same");
    fitInvXSectionPi0900GeV->Draw("same");

    csGraphsSys[2]     ->Draw("E2same");
    csGraphsSys[1]     ->Draw("E2same");
    csGraphsSys[0]     ->Draw("E2same");
    csGraphs[2]->Draw("p,same,z");
    csGraphs[1]->Draw("p,same,z");
    csGraphs[0]->Draw("p,same,z");

    drawLatexAdd("ALICE",rightalignDouble,0.91,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#eta #rightarrow #gamma#gamma",rightalignDouble,0.86,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE);

    TLegend* legendXsectionPaperSingle    = GetAndSetLegend2(0.17, 0.162, 0.5, 0.232+0.05*5, textSizeLabelsPixel);
    legendXsectionPaperSingle->SetNColumns(1);
    legendXsectionPaperSingle->SetMargin(0.2);
    legendXsectionPaperSingle->AddEntry(csGraphsSys[2],Form("%s (x 10^{2})",nameMeasGlobal[2].Data()),"pf");
    legendXsectionPaperSingle->AddEntry(csGraphsSys[1],Form("%s (x 10)",nameMeasGlobal[1].Data()),"pf");
    legendXsectionPaperSingle->AddEntry(csGraphsSys[0],nameMeasGlobal[0].Data(),"pf");
    legendXsectionPaperSingle->AddEntry(fitTCMInvXSectionPi0Plot8TeV,"TCM fit","l");
    legendXsectionPaperSingle->AddEntry(fitInvXSectionPi08TeV,"Tsallis fit","l");

    legendXsectionPaperSingle->AddEntry(histoPythia8InvXSectionEta8TeV,"PYTHIA 8.2, Monash 2013","l");
    legendXsectionPaperSingle->AddEntry((TObject*)0, "NLO, PDF:CTEQ6M5 - FF:AESSS", "");
    TLegend* legendXsectionPaperNLOSingle    = GetAndSetLegend2(0.23, 0.12, 0.83, 0.17, textSizeLabelsPixel);
    legendXsectionPaperNLOSingle->SetNColumns(3);
    legendXsectionPaperNLOSingle->AddEntry(graphNLOCalcEtaMuHalf8TeV, "#mu = 0.5 #it{p}_{T}", "l");
    legendXsectionPaperNLOSingle->AddEntry(graphNLOCalcEtaMuOne8TeV,  "#mu = #it{p}_{T}", "l");
    legendXsectionPaperNLOSingle->AddEntry(graphNLOCalcEtaMuTwo8TeV,  "#mu = 2 #it{p}_{T}", "l");
    legendXsectionPaperSingle->Draw();
    legendXsectionPaperNLOSingle->Draw();

    canvasInvSectionPaperSingle->Print(Form("%s/Eta_InvCrossSectionSingle.%s",outputDir.Data(),suffix.Data()));



    Double_t arrayBoundariesX1_XSec3[2];
    Double_t arrayBoundariesY1_XSec3[4];
    Double_t relativeMarginsXXSec3[3];
    Double_t relativeMarginsYXSec3[3];
    textSizeLabelsPixel = 48;
    ReturnCorrectValuesForCanvasScaling(1250,1300, 1, 3,0.135, 0.005, 0.003,0.072,arrayBoundariesX1_XSec3,arrayBoundariesY1_XSec3,relativeMarginsXXSec3,relativeMarginsYXSec3);

    TCanvas* canvasInvSectionPaperSingleRatio      = new TCanvas("canvasInvSectionPaperSingleRatio","",0,0,1250,1250);  // gives the page size
    DrawGammaCanvasSettings( canvasInvSectionPaperSingleRatio,  0.13, 0.02, 0.03, 0.06);

    TPad* padInvSectionNLORatioSingle         = new TPad("padInvSectionNLORatioSingle", "", arrayBoundariesX1_XSec3[0], arrayBoundariesY1_XSec3[1], arrayBoundariesX1_XSec3[1], arrayBoundariesY1_XSec3[0],-1, -1, -2);
    DrawGammaPadSettings( padInvSectionNLORatioSingle, relativeMarginsXXSec3[0], relativeMarginsXXSec3[2], relativeMarginsYXSec3[0], relativeMarginsYXSec3[1]);
    padInvSectionNLORatioSingle->Draw();
    textsizeLabelsXSecMiddle   = 0;
    textsizeFacXSecMiddle      = 0;
    if (padInvSectionNLORatioSingle->XtoPixel(padInvSectionNLORatioSingle->GetX2()) < padInvSectionNLORatioSingle->YtoPixel(padInvSectionNLORatioSingle->GetY1())){
        textsizeLabelsXSecMiddle        = (Double_t)textSizeLabelsPixel/padInvSectionNLORatioSingle->XtoPixel(padInvSectionNLORatioSingle->GetX2()) ;
        textsizeFacXSecMiddle           = (Double_t)1./padInvSectionNLORatioSingle->XtoPixel(padInvSectionNLORatioSingle->GetX2()) ;
    } else {
        textsizeLabelsXSecMiddle        = (Double_t)textSizeLabelsPixel/padInvSectionNLORatioSingle->YtoPixel(padInvSectionNLORatioSingle->GetY1());
        textsizeFacXSecMiddle           = (Double_t)1./padInvSectionNLORatioSingle->YtoPixel(padInvSectionNLORatioSingle->GetY1());
    }

    TPad* padInvSectionPythiaRatioSignle      = new TPad("padInvSectionPythiaRatioSignle", "", arrayBoundariesX1_XSec3[0], arrayBoundariesY1_XSec3[2], arrayBoundariesX1_XSec3[1], arrayBoundariesY1_XSec3[1],-1, -1, -2);
    DrawGammaPadSettings( padInvSectionPythiaRatioSignle, relativeMarginsXXSec3[0], relativeMarginsXXSec3[2], relativeMarginsYXSec3[1], relativeMarginsYXSec3[1]);
    padInvSectionPythiaRatioSignle->Draw();
    textsizeLabelsXSecDown     = 0;
    textsizeFacXSecDown        = 0;
    if (padInvSectionPythiaRatioSignle->XtoPixel(padInvSectionPythiaRatioSignle->GetX2()) < padInvSectionPythiaRatioSignle->YtoPixel(padInvSectionPythiaRatioSignle->GetY1())){
        textsizeLabelsXSecDown          = (Double_t)textSizeLabelsPixel/padInvSectionPythiaRatioSignle->XtoPixel(padInvSectionPythiaRatioSignle->GetX2()) ;
        textsizeFacXSecDown             = (Double_t)1./padInvSectionPythiaRatioSignle->XtoPixel(padInvSectionPythiaRatioSignle->GetX2()) ;
    } else {
        textsizeLabelsXSecDown          = (Double_t)textSizeLabelsPixel/padInvSectionPythiaRatioSignle->YtoPixel(padInvSectionPythiaRatioSignle->GetY1());
        textsizeFacXSecDown             = (Double_t)1./padInvSectionPythiaRatioSignle->YtoPixel(padInvSectionPythiaRatioSignle->GetY1());
    }
    TPad* padInvSection900GeVRatioSingle      = new TPad("padInvSection900GeVRatioSingle", "", arrayBoundariesX1_XSec3[0], arrayBoundariesY1_XSec3[3], arrayBoundariesX1_XSec3[1], arrayBoundariesY1_XSec3[2],-1, -1, -2);
    DrawGammaPadSettings( padInvSection900GeVRatioSingle, relativeMarginsXXSec3[0], relativeMarginsXXSec3[2], relativeMarginsYXSec3[1], relativeMarginsYXSec3[2]);
    padInvSection900GeVRatioSingle->Draw();
    textsizeLabelsXSecDown2     = 0;
    textsizeFacXSecDown2        = 0;
    if (padInvSection900GeVRatioSingle->XtoPixel(padInvSection900GeVRatioSingle->GetX2()) < padInvSection900GeVRatioSingle->YtoPixel(padInvSection900GeVRatioSingle->GetY1())){
        textsizeLabelsXSecDown2          = (Double_t)textSizeLabelsPixel/padInvSection900GeVRatioSingle->XtoPixel(padInvSection900GeVRatioSingle->GetX2()) ;
        textsizeFacXSecDown2             = (Double_t)1./padInvSection900GeVRatioSingle->XtoPixel(padInvSection900GeVRatioSingle->GetX2()) ;
    } else {
        textsizeLabelsXSecDown2          = (Double_t)textSizeLabelsPixel/padInvSection900GeVRatioSingle->YtoPixel(padInvSection900GeVRatioSingle->GetY1());
        textsizeFacXSecDown2             = (Double_t)1./padInvSection900GeVRatioSingle->YtoPixel(padInvSection900GeVRatioSingle->GetY1());
    }

    padInvSectionNLORatioSingle->cd();
    padInvSectionNLORatioSingle->SetLogx(1);
    ratio8TeVdummy->DrawCopy();
    graphRatioEtaCombNLOMuHalf8TeV->Draw("same,c");
    graphRatioEtaCombNLOMuOne8TeV->Draw("same,c");
    graphRatioEtaCombNLOMuTwo8TeV->Draw("same,c");
    histoRatioPythia8ToFit8TeV->Draw("same,hist,l");
    graphRatioPythia8ToFit8TeV->Draw("3,same");
    graphRatioPi0CombCombFitSysA->Draw("2,same");
    graphRatioPi0CombCombFitStatA_WOXErr->Draw("p,same");
    DrawGammaLines(0.33, maxX , 1., 1.,0.5, kGray+2);

    padInvSectionPythiaRatioSignle->cd();
    padInvSectionPythiaRatioSignle->SetLogx(1);
    ratio7TeVdummy->DrawCopy();
    graphRatioEtaCombNLOMuHalf7TeV->Draw("same,c");
    graphRatioEtaCombNLOMuOne7TeV->Draw("same,c");
    graphRatioEtaCombNLOMuTwo7TeV->Draw("same,c");
    histoRatioPythia8ToFit7TeV->Draw("same,hist,l");
    graphRatioPythia8ToFit7TeV->Draw("3,same");
    graphRatioPi0CombCombFitSysA7TeV->Draw("2,same");
    graphRatioPi0CombCombFitStatA7TeV_WOXErr->Draw("p,same");
    DrawGammaLines(0.33, maxX , 1., 1.,0.5, kGray+2);

    padInvSection900GeVRatioSingle->cd();
    padInvSection900GeVRatioSingle->SetLogx(1);
    ratio900GeVdummy->DrawCopy();
    graphRatioEtaCombNLOMuHalf900GeV->Draw("same,c");
    graphRatioEtaCombNLOMuOne900GeV->Draw("same,c");
    graphRatioEtaCombNLOMuTwo900GeV->Draw("same,c");
    histoRatioPythia8ToFit900GeV->Draw("same,hist,l");
    graphRatioPythia8ToFit900GeV->Draw("3,same");
    graphRatioPi0CombCombFitSysA900GeV->Draw("2,same");
    graphRatioPi0CombCombFitStatA900GeV_WOXErr->Draw("p,same");
    DrawGammaLines(0.33, maxX , 1., 1.,0.5, kGray+2);

    canvasInvSectionPaperSingleRatio->Print(Form("%s/Eta_InvCrossSectionRatioSingle.%s",outputDir.Data(),suffix.Data()));


}

void plotDoubleRatio(TGraphAsymmErrors* csGraphs[],TGraphAsymmErrors* csGraphsSys[], Double_t yMin, Double_t yMax, TString meson, TString plotName){

    Double_t maxX = 19.9;
    Double_t yrangeLowRatio = 0.71;
    Double_t yrangeHighRatio = 1.59;

    Style_t markerstylesFULL[4]                         ={34,20,33,29};
    Size_t markersizeFULL[4]                        ={2.3,2.3,3.4,3.2};
    Color_t colorDataFULL[4]                        ={kRed+2,kBlue+2,kGreen+2,kMagenta+2};

    Double_t arrayBoundariesX1_XSec3[2];
    Double_t arrayBoundariesY1_XSec3[4];
    Double_t relativeMarginsXXSec3[3];
    Double_t relativeMarginsYXSec3[3];
    Double_t textSizeLabelsPixel = 48;
    //     ReturnCorrectValuesForCanvasScaling(1250,1300, 1, 3,0.135, 0.005, 0.005,0.082,arrayBoundariesX1_XSec3,arrayBoundariesY1_XSec3,relativeMarginsXXSec3,relativeMarginsYXSec3);
    ReturnCorrectValuesForCanvasScaling(1250,1300, 1, 3,0.1, 0.005, 0.005,0.082,arrayBoundariesX1_XSec3,arrayBoundariesY1_XSec3,relativeMarginsXXSec3,relativeMarginsYXSec3);

    TCanvas* canvasInvSectionPaperSingleRatio      = new TCanvas("canvasInvSectionPaperSingleRatio","",0,0,1250,1250);  // gives the page size
    DrawGammaCanvasSettings( canvasInvSectionPaperSingleRatio,  0.13, 0.02, 0.03, 0.06);

    TPad* padInvSectionNLORatioSingle         = new TPad("padInvSectionNLORatioSingle", "", arrayBoundariesX1_XSec3[0], arrayBoundariesY1_XSec3[1], arrayBoundariesX1_XSec3[1], arrayBoundariesY1_XSec3[0],-1, -1, -2);
    DrawGammaPadSettings( padInvSectionNLORatioSingle, relativeMarginsXXSec3[0], relativeMarginsXXSec3[2], relativeMarginsYXSec3[0], relativeMarginsYXSec3[1]);
    padInvSectionNLORatioSingle->Draw();
    Double_t textsizeLabelsXSecMiddle   = 0;
    Double_t textsizeFacXSecMiddle      = 0;
    if (padInvSectionNLORatioSingle->XtoPixel(padInvSectionNLORatioSingle->GetX2()) < padInvSectionNLORatioSingle->YtoPixel(padInvSectionNLORatioSingle->GetY1())){
        textsizeLabelsXSecMiddle        = (Double_t)textSizeLabelsPixel/padInvSectionNLORatioSingle->XtoPixel(padInvSectionNLORatioSingle->GetX2()) ;
        textsizeFacXSecMiddle           = (Double_t)1./padInvSectionNLORatioSingle->XtoPixel(padInvSectionNLORatioSingle->GetX2()) ;
    } else {
        textsizeLabelsXSecMiddle        = (Double_t)textSizeLabelsPixel/padInvSectionNLORatioSingle->YtoPixel(padInvSectionNLORatioSingle->GetY1());
        textsizeFacXSecMiddle           = (Double_t)1./padInvSectionNLORatioSingle->YtoPixel(padInvSectionNLORatioSingle->GetY1());
    }

    TPad* padInvSectionPythiaRatioSignle      = new TPad("padInvSectionPythiaRatioSignle", "", arrayBoundariesX1_XSec3[0], arrayBoundariesY1_XSec3[2], arrayBoundariesX1_XSec3[1], arrayBoundariesY1_XSec3[1],-1, -1, -2);
    DrawGammaPadSettings( padInvSectionPythiaRatioSignle, relativeMarginsXXSec3[0], relativeMarginsXXSec3[2], relativeMarginsYXSec3[1], relativeMarginsYXSec3[1]);
    padInvSectionPythiaRatioSignle->Draw();
    Double_t textsizeLabelsXSecDown     = 0;
    Double_t textsizeFacXSecDown        = 0;
    if (padInvSectionPythiaRatioSignle->XtoPixel(padInvSectionPythiaRatioSignle->GetX2()) < padInvSectionPythiaRatioSignle->YtoPixel(padInvSectionPythiaRatioSignle->GetY1())){
        textsizeLabelsXSecDown          = (Double_t)textSizeLabelsPixel/padInvSectionPythiaRatioSignle->XtoPixel(padInvSectionPythiaRatioSignle->GetX2()) ;
        textsizeFacXSecDown             = (Double_t)1./padInvSectionPythiaRatioSignle->XtoPixel(padInvSectionPythiaRatioSignle->GetX2()) ;
    } else {
        textsizeLabelsXSecDown          = (Double_t)textSizeLabelsPixel/padInvSectionPythiaRatioSignle->YtoPixel(padInvSectionPythiaRatioSignle->GetY1());
        textsizeFacXSecDown             = (Double_t)1./padInvSectionPythiaRatioSignle->YtoPixel(padInvSectionPythiaRatioSignle->GetY1());
    }
    TPad* padInvSection900GeVRatioSingle      = new TPad("padInvSection900GeVRatioSingle", "", arrayBoundariesX1_XSec3[0], arrayBoundariesY1_XSec3[3], arrayBoundariesX1_XSec3[1], arrayBoundariesY1_XSec3[2],-1, -1, -2);
    DrawGammaPadSettings( padInvSection900GeVRatioSingle, relativeMarginsXXSec3[0], relativeMarginsXXSec3[2], relativeMarginsYXSec3[1], relativeMarginsYXSec3[2]);
    padInvSection900GeVRatioSingle->Draw();
    Double_t textsizeLabelsXSecDown2     = 0;
    Double_t textsizeFacXSecDown2        = 0;
    if (padInvSection900GeVRatioSingle->XtoPixel(padInvSection900GeVRatioSingle->GetX2()) < padInvSection900GeVRatioSingle->YtoPixel(padInvSection900GeVRatioSingle->GetY1())){
        textsizeLabelsXSecDown2          = (Double_t)textSizeLabelsPixel/padInvSection900GeVRatioSingle->XtoPixel(padInvSection900GeVRatioSingle->GetX2()) ;
        textsizeFacXSecDown2             = (Double_t)1./padInvSection900GeVRatioSingle->XtoPixel(padInvSection900GeVRatioSingle->GetX2()) ;
    } else {
        textsizeLabelsXSecDown2          = (Double_t)textSizeLabelsPixel/padInvSection900GeVRatioSingle->YtoPixel(padInvSection900GeVRatioSingle->GetY1());
        textsizeFacXSecDown2             = (Double_t)1./padInvSection900GeVRatioSingle->YtoPixel(padInvSection900GeVRatioSingle->GetY1());
    }
    //     Double_t marginXSec                 = relativeMarginsXXSec3[0]*1250;
    Double_t marginXSec                 = relativeMarginsXXSec3[0]*2000;

    TFile* fileTheoryCompilation                   = new TFile("ExternalInput/Theory/TheoryCompilationPP.root");
    TDirectoryFile*   directoryGamma                          = (TDirectoryFile*)fileTheoryCompilation->Get("DirectPhoton");

    TGraphAsymmErrors* graphDirectPhotonNLO[4];
    TGraphAsymmErrors* graphPromptPhotonNLO[4];
    TGraphAsymmErrors* graphThermalPhotonNLO[4];

    graphDirectPhotonNLO[0]                = (TGraphAsymmErrors*)directoryGamma->Get("graphDirectPhotonNLOVogelsangInvYieldINT1_7TeV");
    graphPromptPhotonNLO[0]                = (TGraphAsymmErrors*)directoryGamma->Get("graphPromptPhotonNLOVogelsangInvYieldINT1_7TeV");
    graphThermalPhotonNLO[0]         = (TGraphAsymmErrors*)directoryGamma->Get("graphThermalPhotonNLOVogelsangInvYieldINT1_7TeV");

    graphDirectPhotonNLO[1]            = (TGraphAsymmErrors*)directoryGamma->Get("graphDirectPhotonNLOVogelsangInvYieldINT1_7TeV");
    graphPromptPhotonNLO[1]            = (TGraphAsymmErrors*)directoryGamma->Get("graphPromptDirectPhotonLiuWernerInvYieldINT1_7TeV");
    graphThermalPhotonNLO[1]     = (TGraphAsymmErrors*)directoryGamma->Get("graphThermalAndPromptDirectPhotonLiuWernerInvYieldINT1_7TeV");

    graphDirectPhotonNLO[2]                = (TGraphAsymmErrors*)directoryGamma->Get("graphDirectPhotonNLOVogelsangInvYield_8TeV");
    graphPromptPhotonNLO[2]                = (TGraphAsymmErrors*)directoryGamma->Get("graphPromptPhotonNLOVogelsangInvYield_8TeV");
    graphThermalPhotonNLO[2]         = (TGraphAsymmErrors*)directoryGamma->Get("graphThermalPhotonNLOVogelsangInvYield_8TeV");


    graphDirectPhotonNLO[3]                = (TGraphAsymmErrors*)directoryGamma->Get("graphDirectPhotonNLOVogelsangInvYieldINT1_2760GeV");
    graphPromptPhotonNLO[3]                = (TGraphAsymmErrors*)directoryGamma->Get("graphPromptPhotonNLOVogelsangInvYieldINT1_2760GeV");
    graphThermalPhotonNLO[3]         = (TGraphAsymmErrors*)directoryGamma->Get("graphThermalPhotonNLOVogelsangInvYieldINT1_2760GeV");




    TFile* cocktailFile[4];
    cocktailFile[0] = new TFile("ThesisQAInputAllEnergies/900GeV_gamma_00000113_00200009227302008250404000_0152103500000000/900GeV/GammaCocktail_0.80_00000113_00200009227302008250404000_0152103500000000.root");
    cocktailFile[1] = new TFile("ThesisQAInputAllEnergies/7TeV_gamma_00000113_00200009227302008250404000_0152103500000000/7TeV/GammaCocktail_0.80_00000113_00200009227300008250404000_0152103500000000.root");
    cocktailFile[2] = new TFile("ThesisQAInputAllEnergies/8TeV_gamma_00010113_00200009227300008250404000_0152103500000000/8TeV/GammaCocktail_0.80_00010113_00200009227300008250404000_0152103500000000.root");
    cocktailFile[3] = new TFile("CombinationInputPP/2.76TeV/gamma_00000113_00200009397300008250400000_0163103100900000/2.76TeV/GammaCocktail_0.80_00000113_00200009397300008250400000_0163103100900000.root");
    TH1D* cocktailPi0[4];
    TH1D* cocktailEta[4];
    TH1D* cocktailAllGamma[4];
    TH1D* cocktailAllGammaPi0[4];
    for(Int_t i=0; i<4; i++){
        cocktailPi0[i]                         = (TH1D* )cocktailFile[i]->Get("Pi0_Pt");
        cocktailEta[i]                         = (TH1D* )cocktailFile[i]->Get("Eta_Pt");
        cocktailAllGamma[i]                    = (TH1D* )cocktailFile[i]->Get("Gamma_Pt");
        cocktailAllGammaPi0[i]                 = (TH1D* )cocktailFile[i]->Get("Gamma_From_Pi0_Pt");
    }

    Double_t* xVal                                              = NULL;
    Double_t* xErr                                              = NULL;
    Double_t* yVal                                              = NULL;
    Double_t* yErrUp                                            = NULL;
    Double_t* yErrDown                                          = NULL;

    TH1D* cocktailAllGammaNLO[4];
    TGraphAsymmErrors*   graphDirectPhotonNLOCopy[4];
    TGraphAsymmErrors*   graphNLODoubleRatio[6];
    TGraphAsymmErrors*   graphNLODirGammaSpectra[4];
    TGraphAsymmErrors*   graphDirectPhotonThermalNLOCopy[4];
    TGraphAsymmErrors*   graphNLOThermalDoubleRatio[4];
    TGraphAsymmErrors*   graphNLOThermalDirGammaSpectra[4];
    TGraphAsymmErrors*   graphDirectPhotonPromptNLOCopy[4];
    TGraphAsymmErrors*   graphNLOPromptDoubleRatio[4];
    TGraphAsymmErrors*   graphNLOPromptDirGammaSpectra[4];
    for(Int_t i=0; i<4; i++){
        cocktailAllGammaNLO[i]                     = (TH1D*) cocktailAllGamma[i]->Clone("cocktailAllGammaNLO");
        graphDirectPhotonNLOCopy[i]                = (TGraphAsymmErrors*)graphDirectPhotonNLO[i]->Clone("graphNLOCalcCopy");

        xVal                                    = graphDirectPhotonNLOCopy[i]->GetX();
        xErr                                    = graphDirectPhotonNLOCopy[i]->GetEX();
        yVal                                    = graphDirectPhotonNLOCopy[i]->GetY();
        yErrUp                                  = graphDirectPhotonNLOCopy[i]->GetEYhigh();
        yErrDown                                = graphDirectPhotonNLOCopy[i]->GetEYlow();

        TString cocktailFit                     = "xqcd";
        TString fitOptions                      = "QNRME+";//"IQNRME+";

        TF1*   fitCocktailAllGammaForNLO               = (TF1*)FitObject(cocktailFit,"cocktailFit","Pi0",cocktailAllGammaNLO[i],2.0,16,NULL,fitOptions);

        for (Int_t bin=0; bin<graphDirectPhotonNLOCopy[i]->GetN(); bin++) {
            yVal[bin]                           = (1 + ( yVal[bin] / (fitCocktailAllGammaForNLO->Eval(xVal[bin]))));
        }

        // ------------------------------ NLO Calculations -----------------------------
        graphNLODoubleRatio[i]                     = new TGraphAsymmErrors(graphDirectPhotonNLOCopy[i]->GetN(), xVal, yVal, xErr, xErr, yErrDown, yErrUp);
        graphNLODoubleRatio[i]->SetName("graphNLODoubleRatio");
        graphNLODirGammaSpectra[i]                 = (TGraphAsymmErrors*)graphDirectPhotonNLO[i]->Clone("graphNLODirGammaSpectra");

        graphNLODoubleRatio[i]->SetLineColor(kBlack);
        // graphNLODoubleRatio[i]->SetLineColor(kAzure);
        graphNLODoubleRatio[i]->SetFillColor(kBlack);
        graphNLODoubleRatio[i]->SetLineWidth(3.0);
        graphNLODoubleRatio[i]->SetMarkerSize(0);
    }
    /*for(Int_t i=1; i<2; i++){
        cocktailAllGammaNLO[i]                     = (TH1D*) cocktailAllGamma[i]->Clone("cocktailAllGammaNLO");
        graphDirectPhotonPromptNLOCopy[i]                = (TGraphAsymmErrors*)graphPromptPhotonNLO[i]->Clone("graphNLOPromptCalcCopy");
        graphDirectPhotonThermalNLOCopy[i]                = (TGraphAsymmErrors*)graphThermalPhotonNLO[i]->Clone("graphNLOThermalCalcCopy");

        xVal                                    = graphDirectPhotonPromptNLOCopy[i]->GetX();
        xErr                                    = graphDirectPhotonPromptNLOCopy[i]->GetEX();
        yVal                                    = graphDirectPhotonPromptNLOCopy[i]->GetY();
        yErrUp                                  = graphDirectPhotonPromptNLOCopy[i]->GetEYhigh();
        yErrDown                                = graphDirectPhotonPromptNLOCopy[i]->GetEYlow();

        TString cocktailFit                     = "xqcd";
        TString fitOptions                      = "QNRME+";//"IQNRME+";

        TF1*   fitCocktailAllGammaForNLO               = (TF1*)FitObject(cocktailFit,"cocktailFit","Pi0",cocktailAllGammaNLO[i],2.0,16,NULL,fitOptions);

        for (Int_t bin=0; bin<graphDirectPhotonPromptNLOCopy[i]->GetN(); bin++) {
            yVal[bin]                           = (1 + ( yVal[bin] / (fitCocktailAllGammaForNLO->Eval(xVal[bin]))));
        }

        // ------------------------------ NLO Calculations -----------------------------
        graphNLODoubleRatio[4]                     = new TGraphAsymmErrors(graphDirectPhotonPromptNLOCopy[i]->GetN(), xVal, yVal, xErr, xErr, yErrDown, yErrUp);
        graphNLODoubleRatio[4]->SetName("graphNLODoubleRatio");
        graphNLODirGammaSpectra[4]                 = (TGraphAsymmErrors*)graphDirectPhotonNLO[i]->Clone("graphNLODirGammaSpectra");

        // graphNLODoubleRatio[4]->SetLineColor(kBlack);
        graphNLODoubleRatio[4]->SetLineColor(kAzure);
        graphNLODoubleRatio[4]->SetMarkerColor(kAzure);
        // graphNLODoubleRatio[4]->SetFillColor(kBlack);
        graphNLODoubleRatio[4]->SetFillColor(kAzure);
        graphNLODoubleRatio[4]->SetLineWidth(3.0);
        graphNLODoubleRatio[4]->SetMarkerSize(0);

        xVal                                    = graphDirectPhotonThermalNLOCopy[i]->GetX();
        xErr                                    = graphDirectPhotonThermalNLOCopy[i]->GetEX();
        yVal                                    = graphDirectPhotonThermalNLOCopy[i]->GetY();
        yErrUp                                  = graphDirectPhotonThermalNLOCopy[i]->GetEYhigh();
        yErrDown                                = graphDirectPhotonThermalNLOCopy[i]->GetEYlow();

        cocktailFit                     = "xqcd";
        fitOptions                      = "QNRME+";//"IQNRME+";

        fitCocktailAllGammaForNLO               = (TF1*)FitObject(cocktailFit,"cocktailFit","Pi0",cocktailAllGammaNLO[i],2.0,16,NULL,fitOptions);

        for (Int_t bin=0; bin<graphDirectPhotonThermalNLOCopy[i]->GetN(); bin++) {
            yVal[bin]                           = (1 + ( yVal[bin] / (fitCocktailAllGammaForNLO->Eval(xVal[bin]))));
        }

        // ------------------------------ NLO Calculations -----------------------------
        graphNLODoubleRatio[5]                     = new TGraphAsymmErrors(graphDirectPhotonThermalNLOCopy[i]->GetN(), xVal, yVal, xErr, xErr, yErrDown, yErrUp);
        graphNLODoubleRatio[5]->SetName("graphNLODoubleRatio");
        graphNLODirGammaSpectra[5]                 = (TGraphAsymmErrors*)graphDirectPhotonNLO[i]->Clone("graphNLODirGammaSpectra");

        // graphNLODoubleRatio[5]->SetLineColor(kBlack);
        graphNLODoubleRatio[5]->SetLineColor(kRed+1);
        graphNLODoubleRatio[5]->SetMarkerColor(kRed+1);
        // graphNLODoubleRatio[5]->SetFillColor(kBlack);
        graphNLODoubleRatio[5]->SetFillColor(kRed+1);
        graphNLODoubleRatio[5]->SetLineWidth(3.0);
        graphNLODoubleRatio[5]->SetMarkerSize(0);
    }
    */
    //############### PAD 8 TeV  ###########################
    padInvSectionNLORatioSingle->cd();
    padInvSectionNLORatioSingle->SetLogx(1);
    TH2F * ratio8TeVdummy               = new TH2F("ratio8TeVdummy","ratio8TeVdummy",1000,0.23,maxX,1000,yrangeLowRatio,yrangeHighRatio);
    SetStyleHistoTH2ForGraphs(ratio8TeVdummy, "#it{p}_{T} (GeV/#it{c})","", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle,
                            0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.2/(textsizeFacXSecMiddle*marginXSec), 510, 505);
    ratio8TeVdummy->GetYaxis()->SetMoreLogLabels(kTRUE);
    ratio8TeVdummy->GetYaxis()->SetNdivisions(505);
    ratio8TeVdummy->GetYaxis()->SetNoExponent(kTRUE);
    ratio8TeVdummy->GetXaxis()->SetMoreLogLabels(kTRUE);
    ratio8TeVdummy->GetXaxis()->SetNoExponent(kTRUE);
    ratio8TeVdummy->GetXaxis()->SetLabelFont(42);
    ratio8TeVdummy->GetYaxis()->SetLabelFont(42);
    ratio8TeVdummy->GetYaxis()->SetLabelOffset(+0.01);
    ratio8TeVdummy->GetXaxis()->SetTickLength(0.07);
    ratio8TeVdummy->DrawCopy();

    while(graphNLODoubleRatio[2]->GetX()[graphNLODoubleRatio[2]->GetN()-1] > 18.) graphNLODoubleRatio[2]->RemovePoint(graphNLODoubleRatio[2]->GetN()-1);
    graphNLODoubleRatio[2]->Draw("lp3");

    TGraphAsymmErrors* graphDoubleRatioNoXErrStat8TeV    = (TGraphAsymmErrors*)csGraphs[2]->Clone();
    ProduceGraphAsymmWithoutXErrors(graphDoubleRatioNoXErrStat8TeV);
    TGraphAsymmErrors* graphDoubleRatioSys8TeV     = (TGraphAsymmErrors*)csGraphsSys[2]->Clone();
    DrawGammaSetMarkerTGraphAsym(graphDoubleRatioSys8TeV, markerstylesFULL[2], markersizeFULL[2], colorDataFULL[2], colorDataFULL[2], widthLinesBoxes, kTRUE, 0);
    graphDoubleRatioSys8TeV->SetLineWidth(0);
    graphDoubleRatioSys8TeV->Draw("2,same");
    DrawGammaSetMarkerTGraphAsym(graphDoubleRatioNoXErrStat8TeV, markerstylesFULL[2], markersizeFULL[2], colorDataFULL[2], colorDataFULL[2], widthLinesBoxes, kFALSE);
    graphDoubleRatioNoXErrStat8TeV->SetLineWidth(widthLinesBoxes);
    graphDoubleRatioNoXErrStat8TeV->Draw("p,same");

    TLegend* legendDR8TeV    = GetAndSetLegend2(0.14, 0.67, 0.5, 0.93, textSizeLabelsPixel);
    legendDR8TeV->SetNColumns(1);
    legendDR8TeV->SetMargin(0.2);
    //     legendDR8TeV->AddEntry(graphDoubleRatioSys8TeV,"Direct photon double ratio","pf");
    legendDR8TeV->AddEntry(graphDoubleRatioSys8TeV,Form("pp %s",nameMeasGlobal[2].Data()),"pf");
    //     legendDR8TeV->AddEntry(graphDoubleRatioSys8TeV,"PCM #gamma_{dir} double ratio","pf");
    legendDR8TeV->AddEntry(graphNLODoubleRatio[2],"NLO pQCD PDF: CTEQ6M5 FF: GRV","l");
    legendDR8TeV->Draw();

    drawLatexAdd(Form("pp %s",nameMeasGlobal[2].Data()),0.14,0.84,textsizeLabelsXSecMiddle,kFALSE,kFALSE);
    drawLatexAdd("ALICE",0.95,0.84,textsizeLabelsXSecMiddle,kFALSE,kFALSE,kTRUE);
    //     drawLatexAdd("#gamma's rec. with PCM",0.95,0.74,textsizeLabelsXSecMiddle,kFALSE,kFALSE,kTRUE);
    DrawGammaLines(0.23, maxX , 1., 1.,0.5, kGray+2);

    //############### PAD 7 TeV  ###########################
    padInvSectionPythiaRatioSignle->cd();
    padInvSectionPythiaRatioSignle->SetLogx(1);
    TH2F * ratio7TeVdummy               = new TH2F("ratio7TeVdummy","ratio7TeVdummy",1000,0.23,maxX,1000,yrangeLowRatio,yrangeHighRatio);
    //     SetStyleHistoTH2ForGraphs(ratio7TeVdummy, "#it{p}_{T} (GeV/#it{c})","R_{#gamma}=(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
    SetStyleHistoTH2ForGraphs(ratio7TeVdummy, "#it{p}_{T} (GeV/#it{c})","#it{R}_{#gamma}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
                            0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 1,0.2/(textsizeFacXSecMiddle*marginXSec), 510, 505);
    ratio7TeVdummy->GetYaxis()->SetMoreLogLabels(kTRUE);
    ratio7TeVdummy->GetYaxis()->SetNdivisions(505);
    ratio7TeVdummy->GetYaxis()->CenterTitle(kTRUE);
    ratio7TeVdummy->GetYaxis()->SetNoExponent(kTRUE);
    ratio7TeVdummy->GetXaxis()->SetMoreLogLabels(kTRUE);
    ratio7TeVdummy->GetXaxis()->SetNoExponent(kTRUE);
    ratio7TeVdummy->GetXaxis()->SetLabelFont(42);
    ratio7TeVdummy->GetYaxis()->SetLabelFont(42);
    ratio7TeVdummy->GetYaxis()->SetLabelOffset(+0.01);
    ratio7TeVdummy->GetXaxis()->SetTickLength(0.07);
    ratio7TeVdummy->DrawCopy();
/*
    TString fileNameDoubleRatioPass2            = "ExternalInput/PCM/Gamma_PCMResults_pp.root";
    TFile* fileGammasPass2                            = new TFile(fileNameDoubleRatioPass2.Data());
    DoubleRatioStatError                       = (TH1D*)fileGammasPass2->Get("Gamma_7TeV_pp/DoubleRatioStatError");
    graphDoubleRatioStatError                       = new TGraphAsymmErrors(DoubleRatioStatError);
    ProduceGraphAsymmWithoutXErrors(graphDoubleRatioStatError);
    DoubleRatioSystError     = (TGraphAsymmErrors*)fileGammasPass2->Get("Gamma_7TeV_pp/DoubleRatioSystError");
    Color_t pass2DRcolor = kBlue-8;
    DrawGammaSetMarkerTGraphAsym(DoubleRatioSystError, 24, markersizeFULL[1], pass2DRcolor, pass2DRcolor, widthLinesBoxes, kTRUE, 0);
    DoubleRatioSystError->SetLineWidth(0);
    //     DoubleRatioSystError->Draw("2,same");
    DrawGammaSetMarkerTGraphAsym(graphDoubleRatioStatError, 24, markersizeFULL[1], pass2DRcolor, pass2DRcolor, widthLinesBoxes, kFALSE);
    graphDoubleRatioStatError->SetLineWidth(widthLinesBoxes);
    //     graphDoubleRatioStatError->Draw("p,same");
*/

    while(graphNLODoubleRatio[1]->GetX()[graphNLODoubleRatio[1]->GetN()-1] > 18.) graphNLODoubleRatio[1]->RemovePoint(graphNLODoubleRatio[1]->GetN()-1);
    // while(graphNLODoubleRatio[4]->GetX()[graphNLODoubleRatio[4]->GetN()-1] > 18.) graphNLODoubleRatio[4]->RemovePoint(graphNLODoubleRatio[4]->GetN()-1);
    // while(graphNLODoubleRatio[5]->GetX()[graphNLODoubleRatio[5]->GetN()-1] > 18.) graphNLODoubleRatio[5]->RemovePoint(graphNLODoubleRatio[5]->GetN()-1);
    graphNLODoubleRatio[1]->Draw("lp3");
    // graphNLODoubleRatio[4]->Draw("lp3");
    // graphNLODoubleRatio[5]->Draw("lp3");

    TGraphAsymmErrors* graphDoubleRatioNoXErrStat7TeV    = (TGraphAsymmErrors*)csGraphs[1]->Clone();
    ProduceGraphAsymmWithoutXErrors(graphDoubleRatioNoXErrStat7TeV);
    TGraphAsymmErrors* graphDoubleRatioSys7TeV     = (TGraphAsymmErrors*)csGraphsSys[1]->Clone();
    DrawGammaSetMarkerTGraphAsym(graphDoubleRatioSys7TeV, markerstylesFULL[1], markersizeFULL[1], colorDataFULL[1], colorDataFULL[1], widthLinesBoxes, kTRUE, 0);
    graphDoubleRatioSys7TeV->SetLineWidth(0);
    graphDoubleRatioSys7TeV->Draw("2,same");
    DrawGammaSetMarkerTGraphAsym(graphDoubleRatioNoXErrStat7TeV, markerstylesFULL[1], markersizeFULL[1], colorDataFULL[1], colorDataFULL[1], widthLinesBoxes, kFALSE);
    graphDoubleRatioNoXErrStat7TeV->SetLineWidth(widthLinesBoxes);
    graphDoubleRatioNoXErrStat7TeV->Draw("p,same");

    // TLegend* legendDR7TeV    = GetAndSetLegend2(0.14, 0.52, 0.5, 0.8, textSizeLabelsPixel);
    TLegend* legendDR7TeV    = GetAndSetLegend2(0.14, 0.77, 0.5, 0.92, textSizeLabelsPixel);
    legendDR7TeV->SetNColumns(1);
    legendDR7TeV->SetMargin(0.2);
    //     legendDR7TeV->AddEntry(graphDoubleRatioSys7TeV,"PCM pass4","pf");
    //     legendDR7TeV->AddEntry(DoubleRatioSystError,"PCM pass2","pf");
    legendDR7TeV->AddEntry(graphDoubleRatioSys7TeV,Form("pp %s",nameMeasGlobal[1].Data()),"pf");
    // legendDR7TeV->AddEntry(graphNLODoubleRatio[1],"NLO prediction","l");
    legendDR7TeV->Draw();

    // drawLatexAdd(Form("pp %s",nameMeasGlobal[1].Data()),0.14,0.84,textsizeLabelsXSecDown,kFALSE,kFALSE);
    DrawGammaLines(0.23, maxX , 1., 1.,0.5, kGray+2);

    //############### PAD 900 GeV  ###########################
    padInvSection900GeVRatioSingle->cd();
    padInvSection900GeVRatioSingle->SetLogx(1);
    TH2F * ratio900GeVdummy            = new TH2F("ratio900GeVdummy","ratio900GeVdummy",1000,0.23,maxX,1000,yrangeLowRatio,yrangeHighRatio);
    SetStyleHistoTH2ForGraphs(ratio900GeVdummy, "#it{p}_{T} (GeV/#it{c})","", 0.85*textsizeLabelsXSecDown2, textsizeLabelsXSecDown2,
                            0.85*textsizeLabelsXSecDown2,textsizeLabelsXSecDown2, 0.9,0.2/(textsizeFacXSecDown2*marginXSec), 510, 505);
    ratio900GeVdummy->GetYaxis()->SetMoreLogLabels(kTRUE);
    ratio900GeVdummy->GetYaxis()->SetNdivisions(505);
    ratio900GeVdummy->GetYaxis()->SetNoExponent(kTRUE);
    ratio900GeVdummy->GetXaxis()->SetMoreLogLabels(kTRUE);
    ratio900GeVdummy->GetXaxis()->SetNoExponent(kTRUE);
    ratio900GeVdummy->GetXaxis()->SetLabelFont(42);
    ratio900GeVdummy->GetYaxis()->SetLabelFont(42);
    ratio900GeVdummy->GetYaxis()->SetLabelOffset(+0.01);
    ratio900GeVdummy->GetXaxis()->SetTickLength(0.06);
    ratio900GeVdummy->GetYaxis()->SetTickLength(0.04);
    ratio900GeVdummy->DrawCopy();


    while(graphNLODoubleRatio[0]->GetX()[graphNLODoubleRatio[0]->GetN()-1] > 5.) graphNLODoubleRatio[0]->RemovePoint(graphNLODoubleRatio[0]->GetN()-1);
    //     graphNLODoubleRatio[0]->Draw("lp3");

    TGraphAsymmErrors* graphDoubleRatioNoXErrStat900GeV    = (TGraphAsymmErrors*)csGraphs[0]->Clone();
    ProduceGraphAsymmWithoutXErrors(graphDoubleRatioNoXErrStat900GeV);
    TGraphAsymmErrors* graphDoubleRatioSys900GeV     = (TGraphAsymmErrors*)csGraphsSys[0]->Clone();
    DrawGammaSetMarkerTGraphAsym(graphDoubleRatioSys900GeV, markerstylesFULL[0], markersizeFULL[0], colorDataFULL[0], colorDataFULL[0], widthLinesBoxes, kTRUE, 0);
    graphDoubleRatioSys900GeV->SetLineWidth(0);
    graphDoubleRatioSys900GeV->Draw("2,same");
    DrawGammaSetMarkerTGraphAsym(graphDoubleRatioNoXErrStat900GeV, markerstylesFULL[0], markersizeFULL[0], colorDataFULL[0], colorDataFULL[0], widthLinesBoxes, kFALSE);
    graphDoubleRatioNoXErrStat900GeV->SetLineWidth(widthLinesBoxes);
    graphDoubleRatioNoXErrStat900GeV->Draw("p,same");

    TLegend* legendDR900GeV    = GetAndSetLegend2(0.14, 0.81, 0.5, 0.93, textSizeLabelsPixel);
    //     TLegend* legendDR900GeV    = GetAndSetLegend2(0.18, 0.58, 0.5, 0.83, textSizeLabelsPixel);
    legendDR900GeV->SetNColumns(1);
    legendDR900GeV->SetMargin(0.2);
    //     legendDR8TeV->AddEntry(graphDoubleRatioSys8TeV,"Direct photon double ratio","pf");
    legendDR900GeV->AddEntry(graphDoubleRatioSys900GeV,Form("pp %s",nameMeasGlobal[0].Data()),"pf");
    //     legendDR900GeV->AddEntry(graphDoubleRatioSys900GeV,"PCM #gamma_{dir} double ratio","pf");
    //     legendDR900GeV->AddEntry(graphNLODoubleRatio[0],"NLO prediction","l");
    legendDR900GeV->Draw();

    // drawLatexAdd(Form("pp %s",nameMeasGlobal[0].Data()),0.14,0.87,textsizeLabelsXSecDown2,kFALSE,kFALSE);
    DrawGammaLines(0.23, maxX , 1., 1.,0.5, kGray+2);

    canvasInvSectionPaperSingleRatio->Print(Form("%s/DoubleRatioFinal.%s",outputDir.Data(),suffix.Data()));













    // #############################################################################################################################
    // #############################################################################################################################
    // #############################################################################################################################
    // #############################################################################################################################



        Double_t arrayBoundariesX1_DR4Pad[2];
    Double_t arrayBoundariesY1_DR4Pad[5];
    Double_t relativeMarginsXDR4Pad[3];
    Double_t relativeMarginsYDR4Pad[4];
    Double_t textSizeLabelsPixelDR4Pad = 48;
    //     ReturnCorrectValuesForCanvasScaling(1250,1300, 1, 3,0.135, 0.005, 0.005,0.082,arrayBoundariesX1_DR4Pad,arrayBoundariesY1_DR4Pad,relativeMarginsXDR4Pad,relativeMarginsYDR4Pad);
    ReturnCorrectValuesForCanvasScaling(1250,1750, 1, 4,0.1, 0.005, 0.005,0.062,arrayBoundariesX1_DR4Pad,arrayBoundariesY1_DR4Pad,relativeMarginsXDR4Pad,relativeMarginsYDR4Pad);

    TCanvas* canvasDR4Pad      = new TCanvas("canvasDR4Pad","",0,0,1250,1750);  // gives the page size
    DrawGammaCanvasSettings( canvasDR4Pad,  0.13, 0.02, 0.03, 0.06);

    TPad* padDR8TeV         = new TPad("padDR8TeV", "", arrayBoundariesX1_DR4Pad[0], arrayBoundariesY1_DR4Pad[1], arrayBoundariesX1_DR4Pad[1], arrayBoundariesY1_DR4Pad[0],-1, -1, -2);
    DrawGammaPadSettings( padDR8TeV, relativeMarginsXDR4Pad[0], relativeMarginsXDR4Pad[2], relativeMarginsYDR4Pad[0], relativeMarginsYDR4Pad[1]);
    padDR8TeV->Draw();
    Double_t textsizeFacDR8TeV   = 0;
    Double_t textsizeFacDR8TeVMid      = 0;
    if (padDR8TeV->XtoPixel(padDR8TeV->GetX2()) < padDR8TeV->YtoPixel(padDR8TeV->GetY1())){
        textsizeFacDR8TeV        = (Double_t)textSizeLabelsPixelDR4Pad/padDR8TeV->XtoPixel(padDR8TeV->GetX2()) ;
        textsizeFacDR8TeVMid           = (Double_t)1./padDR8TeV->XtoPixel(padDR8TeV->GetX2()) ;
    } else {
        textsizeFacDR8TeV        = (Double_t)textSizeLabelsPixelDR4Pad/padDR8TeV->YtoPixel(padDR8TeV->GetY1());
        textsizeFacDR8TeVMid           = (Double_t)1./padDR8TeV->YtoPixel(padDR8TeV->GetY1());
    }

    TPad* padDR7TeV      = new TPad("padDR7TeV", "", arrayBoundariesX1_DR4Pad[0], arrayBoundariesY1_DR4Pad[2], arrayBoundariesX1_DR4Pad[1], arrayBoundariesY1_DR4Pad[1],-1, -1, -2);
    DrawGammaPadSettings( padDR7TeV, relativeMarginsXDR4Pad[0], relativeMarginsXDR4Pad[2], relativeMarginsYDR4Pad[1], relativeMarginsYDR4Pad[1]);
    padDR7TeV->Draw();
    Double_t textsizeFacDR7TeV     = 0;
    Double_t textsizeFacDR7TeVMid        = 0;
    if (padDR7TeV->XtoPixel(padDR7TeV->GetX2()) < padDR7TeV->YtoPixel(padDR7TeV->GetY1())){
        textsizeFacDR7TeV          = (Double_t)textSizeLabelsPixelDR4Pad/padDR7TeV->XtoPixel(padDR7TeV->GetX2()) ;
        textsizeFacDR7TeVMid             = (Double_t)1./padDR7TeV->XtoPixel(padDR7TeV->GetX2()) ;
    } else {
        textsizeFacDR7TeV          = (Double_t)textSizeLabelsPixelDR4Pad/padDR7TeV->YtoPixel(padDR7TeV->GetY1());
        textsizeFacDR7TeVMid             = (Double_t)1./padDR7TeV->YtoPixel(padDR7TeV->GetY1());
    }
    TPad* padDR276TeV      = new TPad("padDR276TeV", "", arrayBoundariesX1_DR4Pad[0], arrayBoundariesY1_DR4Pad[3], arrayBoundariesX1_DR4Pad[1], arrayBoundariesY1_DR4Pad[2],-1, -1, -2);
    DrawGammaPadSettings( padDR276TeV, relativeMarginsXDR4Pad[0], relativeMarginsXDR4Pad[2], relativeMarginsYDR4Pad[1], relativeMarginsYDR4Pad[1]);
    padDR276TeV->Draw();
    Double_t textsizeFacDR276TeV     = 0;
    Double_t textsizeFacDR276TeVMid        = 0;
    if (padDR276TeV->XtoPixel(padDR276TeV->GetX2()) < padDR276TeV->YtoPixel(padDR276TeV->GetY1())){
        textsizeFacDR276TeV          = (Double_t)textSizeLabelsPixelDR4Pad/padDR276TeV->XtoPixel(padDR276TeV->GetX2()) ;
        textsizeFacDR276TeVMid             = (Double_t)1./padDR276TeV->XtoPixel(padDR276TeV->GetX2()) ;
    } else {
        textsizeFacDR276TeV          = (Double_t)textSizeLabelsPixelDR4Pad/padDR276TeV->YtoPixel(padDR276TeV->GetY1());
        textsizeFacDR276TeVMid             = (Double_t)1./padDR276TeV->YtoPixel(padDR276TeV->GetY1());
    }

    TPad* padDR900GeV      = new TPad("padDR900GeV", "", arrayBoundariesX1_DR4Pad[0], arrayBoundariesY1_DR4Pad[4], arrayBoundariesX1_DR4Pad[1], arrayBoundariesY1_DR4Pad[3],-1, -1, -2);
    DrawGammaPadSettings( padDR900GeV, relativeMarginsXDR4Pad[0], relativeMarginsXDR4Pad[2], relativeMarginsYDR4Pad[1], relativeMarginsYDR4Pad[2]);
    padDR900GeV->Draw();
    Double_t textsizeFacDR900GeV     = 0;
    Double_t textsizeFacDR900GeVMid        = 0;
    if (padDR900GeV->XtoPixel(padDR900GeV->GetX2()) < padDR900GeV->YtoPixel(padDR900GeV->GetY1())){
        textsizeFacDR900GeV          = (Double_t)textSizeLabelsPixelDR4Pad/padDR900GeV->XtoPixel(padDR900GeV->GetX2()) ;
        textsizeFacDR900GeVMid             = (Double_t)1./padDR900GeV->XtoPixel(padDR900GeV->GetX2()) ;
    } else {
        textsizeFacDR900GeV          = (Double_t)textSizeLabelsPixelDR4Pad/padDR900GeV->YtoPixel(padDR900GeV->GetY1());
        textsizeFacDR900GeVMid             = (Double_t)1./padDR900GeV->YtoPixel(padDR900GeV->GetY1());
    }



    //############### PAD 8 TeV  ###########################
    padDR8TeV->cd();
    padDR8TeV->SetLogx(1);
    ratio8TeVdummy->DrawCopy();

    graphNLODoubleRatio[2]->Draw("lp3");

    graphDoubleRatioSys8TeV->Draw("2,same");
    graphDoubleRatioNoXErrStat8TeV->Draw("p,same");

    legendDR8TeV->Draw();

    // drawLatexAdd(Form("pp %s",nameMeasGlobal[2].Data()),0.14,0.84,textsizeLabelsXSecMiddle,kFALSE,kFALSE);
    // drawLatexAdd("ALICE",0.14,0.84,textsizeLabelsXSecMiddle,kFALSE,kFALSE);
    // drawLatexAdd("PDF:CTEQ6M5 FF:GRV",0.95,0.69,textsizeLabelsXSecMiddle,kFALSE,kFALSE,kTRUE);
    //     drawLatexAdd("#gamma's rec. with PCM",0.95,0.74,textsizeLabelsXSecMiddle,kFALSE,kFALSE,kTRUE);
    DrawGammaLines(0.23, maxX, 1., 1., 1.2, kGray+2, 7);

    //############### PAD 7 TeV  ###########################
    padDR7TeV->cd();
    padDR7TeV->SetLogx(1);
    ratio7TeVdummy->GetYaxis()->SetTitle("");
    ratio7TeVdummy->DrawCopy();


    // load PHOS input
    TFile* PHOSfile = new TFile("/home/nschmidt/AnalysisResults/pp/7TeV/PHOS/Results_gamma_7TeV_Podist_20170823.root");
    TH1D* hDoubleRatioPHOSStat = (TH1D* )PHOSfile->Get("hDoubleRatio_stat");
    TH1D* hDoubleRatioPHOSSys  = (TH1D* )PHOSfile->Get("hDoubleRatio_sys");
    TGraphAsymmErrors* graphDoubleRatioPHOSStat                       = new TGraphAsymmErrors(hDoubleRatioPHOSStat);
    TGraphAsymmErrors* graphDoubleRatioPHOSSys                         = new TGraphAsymmErrors(hDoubleRatioPHOSSys);
    TGraphAsymmErrors* graphDoubleRatioNoXErrStat7TeVPHOS    = (TGraphAsymmErrors*)graphDoubleRatioPHOSStat->Clone();
    ProduceGraphAsymmWithoutXErrors(graphDoubleRatioNoXErrStat7TeVPHOS);
    Style_t markerPHOS = 24;
    Color_t colorPHOS = kGray+2;
    Style_t sizePHOS = markersizeFULL[1];
    TGraphAsymmErrors* graphDoubleRatioSys7TeVPHOS     = (TGraphAsymmErrors*)graphDoubleRatioPHOSSys->Clone();
    DrawGammaSetMarkerTGraphAsym(graphDoubleRatioSys7TeVPHOS, markerPHOS, sizePHOS, colorPHOS, colorPHOS, widthLinesBoxes, kTRUE, 0);
    graphDoubleRatioSys7TeVPHOS->SetLineWidth(0);
    DrawGammaSetMarkerTGraphAsym(graphDoubleRatioNoXErrStat7TeVPHOS, markerPHOS, sizePHOS, colorPHOS, colorPHOS, widthLinesBoxes, kFALSE);
    graphDoubleRatioNoXErrStat7TeVPHOS->SetLineWidth(widthLinesBoxes);


    graphNLODoubleRatio[1]->Draw("lp3");
    //     while(graphNLODoubleRatio[4]->GetX()[0] < 2) graphNLODoubleRatio[4]->RemovePoint(0);
    //     while(graphNLODoubleRatio[5]->GetX()[0] < 2) graphNLODoubleRatio[5]->RemovePoint(0);
    // graphNLODoubleRatio[4]->Draw("lp3");
    // graphNLODoubleRatio[5]->Draw("lp3");
    graphDoubleRatioSys7TeV->Draw("2,same");
    graphDoubleRatioNoXErrStat7TeV->Draw("p,same");

    graphDoubleRatioSys7TeVPHOS->Draw("2,same");
    graphDoubleRatioNoXErrStat7TeVPHOS->Draw("p,same");


    // legendDR7TeV->Draw();
    // TLegend* legendDRX27TeV    = GetAndSetLegend2(0.14, 0.57, 0.5, 0.92, textSizeLabelsPixel);
    TLegend* legendDRX27TeV    = GetAndSetLegend2(0.14, 0.57, 0.5, 0.92, textSizeLabelsPixel);
    legendDRX27TeV->SetNColumns(1);
    legendDRX27TeV->SetMargin(0.2);
        legendDRX27TeV->AddEntry(graphDoubleRatioSys7TeV,Form("pp %s PCM",nameMeasGlobal[1].Data()),"pf");
        legendDRX27TeV->AddEntry(graphDoubleRatioSys7TeVPHOS,Form("pp %s PHOS",nameMeasGlobal[1].Data()),"pf");
        // legendDRX27TeV->AddEntry(DoubleRatioSystError,"PCM pass2","pf");
    // legendDRX27TeV->AddEntry(graphDoubleRatioSys7TeV,Form("pp %s",nameMeasGlobal[1].Data()),"pf");
    // legendDRX27TeV->AddEntry(graphNLODoubleRatio[4],"Prompt (Phys.Rev.Lett 106, 242301)","l");
    // legendDRX27TeV->AddEntry(graphNLODoubleRatio[5],"Prompt + Thermal","l");
    legendDRX27TeV->Draw();

    // drawLatexAdd(Form("pp %s",nameMeasGlobal[1].Data()),0.14,0.84,textsizeLabelsXSecDown,kFALSE,kFALSE);
    DrawGammaLines(0.23, maxX, 1., 1., 1.2, kGray+2, 7);

    //############### PAD 2.76 TeV  ###########################
    padDR276TeV->cd();
    padDR276TeV->SetLogx(1);


     TH2F * DR276TeVdummy               = new TH2F("DR276TeVdummy","DR276TeVdummy",1000,0.23,maxX,1000,yrangeLowRatio,yrangeHighRatio);
    //     SetStyleHistoTH2ForGraphs(DR276TeVdummy, "#it{p}_{T} (GeV/#it{c})","R_{#gamma}=(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
    SetStyleHistoTH2ForGraphs(DR276TeVdummy, "#it{p}_{T} (GeV/#it{c})","#it{R}_{#gamma}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
                            0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 1,0.2/(textsizeFacXSecMiddle*marginXSec), 510, 505);
    DR276TeVdummy->GetYaxis()->SetMoreLogLabels(kTRUE);
    DR276TeVdummy->GetYaxis()->SetNdivisions(505);
    // DR276TeVdummy->GetYaxis()->CenterTitle(kTRUE);
    DR276TeVdummy->GetYaxis()->SetNoExponent(kTRUE);
    DR276TeVdummy->GetXaxis()->SetMoreLogLabels(kTRUE);
    DR276TeVdummy->GetXaxis()->SetNoExponent(kTRUE);
    DR276TeVdummy->GetXaxis()->SetLabelFont(42);
    DR276TeVdummy->GetYaxis()->SetLabelFont(42);
    DR276TeVdummy->GetYaxis()->SetLabelOffset(+0.01);
    DR276TeVdummy->GetXaxis()->SetTickLength(0.07);
    DR276TeVdummy->DrawCopy();



    // graphNLODoubleRatio[1]->Draw("lp3");
    while(graphNLODoubleRatio[3]->GetX()[graphNLODoubleRatio[3]->GetN()-1] > 18.) graphNLODoubleRatio[3]->RemovePoint(graphNLODoubleRatio[3]->GetN()-1);
    graphNLODoubleRatio[3]->Draw("lp3");

    TGraphAsymmErrors* graphDoubleRatioNoXErrStat276TeV    = (TGraphAsymmErrors*)csGraphs[3]->Clone();
    ProduceGraphAsymmWithoutXErrors(graphDoubleRatioNoXErrStat276TeV);
    TGraphAsymmErrors* graphDoubleRatioSys276TeV     = (TGraphAsymmErrors*)csGraphsSys[3]->Clone();
    DrawGammaSetMarkerTGraphAsym(graphDoubleRatioSys276TeV, markerstylesFULL[3], markersizeFULL[3], colorDataFULL[3], colorDataFULL[3], widthLinesBoxes, kTRUE, 0);
    graphDoubleRatioSys276TeV->SetLineWidth(0);
    graphDoubleRatioSys276TeV->Draw("2,same");
    DrawGammaSetMarkerTGraphAsym(graphDoubleRatioNoXErrStat276TeV, markerstylesFULL[3], markersizeFULL[3], colorDataFULL[3], colorDataFULL[3], widthLinesBoxes, kFALSE);
    graphDoubleRatioNoXErrStat276TeV->SetLineWidth(widthLinesBoxes);
    graphDoubleRatioNoXErrStat276TeV->Draw("p,same");

    TLegend* legendDR276TeV    = GetAndSetLegend2(0.14, 0.77, 0.5, 0.92, textSizeLabelsPixel);
    legendDR276TeV->SetNColumns(1);
    legendDR276TeV->SetMargin(0.2);
    //     legendDR276TeV->AddEntry(graphDoubleRatioSys276TeV,"Direct photon double ratio","pf");
    legendDR276TeV->AddEntry(graphDoubleRatioSys276TeV,Form("pp %s",nameMeasGlobal[3].Data()),"pf");
    //     legendDR276TeV->AddEntry(graphDoubleRatioSys276TeV,"PCM #gamma_{dir} double ratio","pf");
    // legendDR276TeV->AddEntry(graphNLODoubleRatio[2],"NLO prediction!!!","l");
    legendDR276TeV->Draw();

    // drawLatexAdd(Form("pp %s",nameMeasGlobal[3].Data()),0.14,0.84,textsizeLabelsXSecMiddle,kFALSE,kFALSE);
    //     drawLatexAdd("#gamma's rec. with PCM",0.95,0.74,textsizeLabelsXSecMiddle,kFALSE,kFALSE,kTRUE);
    DrawGammaLines(0.23, maxX, 1., 1., 1.2, kGray+2, 7);


    //############### PAD 900 GeV  ###########################
    padDR900GeV->cd();
    padDR900GeV->SetLogx(1);
    ratio900GeVdummy->DrawCopy();

    graphDoubleRatioSys900GeV->Draw("2,same");
    graphDoubleRatioNoXErrStat900GeV->Draw("p,same");

    legendDR900GeV->Draw();
    // drawLatexAdd("ALICE",0.95,0.84,0.8*textsizeLabelsXSecMiddle,kFALSE,kFALSE,kTRUE);

    // drawLatexAdd(Form("pp %s",nameMeasGlobal[0].Data()),0.14,0.87,textsizeLabelsXSecDown2,kFALSE,kFALSE);
    DrawGammaLines(0.23, maxX, 1., 1., 1.2, kGray+2, 7);

    canvasDR4Pad->Print(Form("%s/DoubleRatioFinal4Pad.%s",outputDir.Data(),suffix.Data()));


    //
    // #############################################################################################################################
    // #############################################################################################################################
    // #############################################################################################################################
    // #############################################################################################################################

































    // DIRECT PHOTON SPECTRUM
    // DIRECT PHOTON SPECTRUM
    // DIRECT PHOTON SPECTRUM
    // DIRECT PHOTON SPECTRUM

    Double_t arrayBoundariesX1_XSec2[2];
    Double_t arrayBoundariesY1_XSec2[3];
    Double_t relativeMarginsXXSec2[3];
    Double_t relativeMarginsYXSec2[3];
    textSizeLabelsPixel = 48;
    ReturnCorrectValuesForCanvasScaling(1250,1500, 1, 1,0.135, 0.005, 0.003,0.07,arrayBoundariesX1_XSec2,arrayBoundariesY1_XSec2,relativeMarginsXXSec2,relativeMarginsYXSec2);

    TCanvas* canvasInvSectionPaperSingle      = new TCanvas("canvasInvSectionPaperSingle","",0,0,1250,1500);  // gives the page size
    DrawGammaCanvasSettings( canvasInvSectionPaperSingle,  0.13, 0.02, 0.03, 0.04);

    TPad* padInvSectionSpecSingle             = new TPad("padInvSectionSpecSingle", "", arrayBoundariesX1_XSec2[0], arrayBoundariesY1_XSec2[2], arrayBoundariesX1_XSec2[1], arrayBoundariesY1_XSec2[0],-1, -1, -2);
    DrawGammaPadSettings( padInvSectionSpecSingle, relativeMarginsXXSec2[0], relativeMarginsXXSec2[2], relativeMarginsYXSec2[0], relativeMarginsYXSec2[2]);
    padInvSectionSpecSingle->Draw();
    marginXSec                 = relativeMarginsXXSec2[0]*1250;
    Double_t textsizeLabelsXSecUp       = 0;
    Double_t textsizeFacXSecUp          = 0;
    if (padInvSectionSpecSingle->XtoPixel(padInvSectionSpecSingle->GetX2()) < padInvSectionSpecSingle->YtoPixel(padInvSectionSpecSingle->GetY1())){
        textsizeLabelsXSecUp            = (Double_t)textSizeLabelsPixel/padInvSectionSpecSingle->XtoPixel(padInvSectionSpecSingle->GetX2()) ;
        textsizeFacXSecUp               = (Double_t)1./padInvSectionSpecSingle->XtoPixel(padInvSectionSpecSingle->GetX2()) ;
    } else {
        textsizeLabelsXSecUp            = (Double_t)textSizeLabelsPixel/padInvSectionSpecSingle->YtoPixel(padInvSectionSpecSingle->GetY1());
        textsizeFacXSecUp               = (Double_t)1./padInvSectionSpecSingle->YtoPixel(padInvSectionSpecSingle->GetY1());
    }

    TH2F * histo2DXSectionPi0;
    //     histo2DXSectionPi0          = new TH2F("histo2DXSectionPi0","histo2DXSectionPi0",11000,0.23,50.,1000,6,9e11);
    histo2DXSectionPi0          = new TH2F("histo2DXSectionPi0","histo2DXSectionPi0",11000,0.23,maxX,1000,1.1e-9,9e5);
    SetStyleHistoTH2ForGraphs(histo2DXSectionPi0, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",0.035,0.04, 0.035,0.04, 0.9,1.45);
    histo2DXSectionPi0->GetXaxis()->SetMoreLogLabels();
    histo2DXSectionPi0->GetXaxis()->SetNoExponent(kTRUE);

    SetStyleHistoTH2ForGraphs(histo2DXSectionPi0, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",
                              0.85*textsizeLabelsXSecUp,textsizeLabelsXSecUp, 0.85*textsizeLabelsXSecUp, textsizeLabelsXSecUp, 0.8,0.2/(textsizeFacXSecUp*marginXSec));
    histo2DXSectionPi0->GetXaxis()->SetMoreLogLabels();
    histo2DXSectionPi0->GetXaxis()->SetLabelOffset(0.0);

    padInvSectionSpecSingle->cd();
    padInvSectionSpecSingle->SetLogy(1);
    padInvSectionSpecSingle->SetLogx(1);

    histo2DXSectionPi0->Draw();

    for (int j=0;j<graphDirectPhotonSpectrum[2]->GetN();j++){
        graphDirectPhotonSpectrum[2]->GetY()[j] *= 1000000;
        graphDirectPhotonSpectrum[2]->GetEYhigh()[j] *= 1000000;
        graphDirectPhotonSpectrum[2]->GetEYlow()[j] *= 1000000;
    }
    for (int j=0;j<graphDirectPhotonSpectrum[1]->GetN();j++){
        graphDirectPhotonSpectrum[1]->GetY()[j] *= 10000;
        graphDirectPhotonSpectrum[1]->GetEYhigh()[j] *= 10000;
        graphDirectPhotonSpectrum[1]->GetEYlow()[j] *= 10000;
    }
    for (int j=0;j<graphDirectPhotonSpectrum[3]->GetN();j++){
        graphDirectPhotonSpectrum[3]->GetY()[j] *= 100;
        graphDirectPhotonSpectrum[3]->GetEYhigh()[j] *= 100;
        graphDirectPhotonSpectrum[3]->GetEYlow()[j] *= 100;
    }
    DrawGammaSetMarkerTGraph(graphDirectPhotonSpectrum[3], markerstylesFULL[3], 0.7*markersizeFULL[3], colorDataFULL[3] , colorDataFULL[3]);
    DrawGammaSetMarkerTGraph(graphDirectPhotonSpectrum[2], markerstylesFULL[2], 0.7*markersizeFULL[2], colorDataFULL[2] , colorDataFULL[2]);
    DrawGammaSetMarkerTGraph(graphDirectPhotonSpectrum[1], markerstylesFULL[1], 0.7*markersizeFULL[1], colorDataFULL[1] , colorDataFULL[1]);
    DrawGammaSetMarkerTGraph(graphDirectPhotonSpectrum[0], markerstylesFULL[0], 0.7*markersizeFULL[0], colorDataFULL[0] , colorDataFULL[0]);

    for (int i=0;i<4;i++){
        graphNLODirGamma[i]->SetLineColor(kBlack);
        graphNLODirGamma[i]->SetFillColor(kBlack);
        graphNLODirGamma[i]->SetLineWidth(3.0);
        graphNLODirGamma[i]->SetMarkerSize(0);
    }

    for (int j=0;j<graphNLODirGamma[2]->GetN();j++){
        graphNLODirGamma[2]->GetY()[j] *= 1000000;
        graphNLODirGamma[2]->GetEYhigh()[j] *= 1000000;
        graphNLODirGamma[2]->GetEYlow()[j] *= 1000000;
    }
    for (int j=0;j<graphNLODirGamma[1]->GetN();j++){
        graphNLODirGamma[1]->GetY()[j] *= 10000;
        graphNLODirGamma[1]->GetEYhigh()[j] *= 10000;
        graphNLODirGamma[1]->GetEYlow()[j] *= 10000;
    }
    for (int j=0;j<graphNLODirGamma[3]->GetN();j++){
        graphNLODirGamma[3]->GetY()[j] *= 100;
        graphNLODirGamma[3]->GetEYhigh()[j] *= 100;
        graphNLODirGamma[3]->GetEYlow()[j] *= 100;
    }
    while(graphNLODirGamma[2]->GetX()[graphNLODirGamma[2]->GetN()-1] > 12.) graphNLODirGamma[2]->RemovePoint(graphNLODirGamma[2]->GetN()-1);
    graphNLODirGamma[2]->Draw("lp3");
    while(graphNLODirGamma[1]->GetX()[graphNLODirGamma[1]->GetN()-1] > 17.) graphNLODirGamma[1]->RemovePoint(graphNLODirGamma[1]->GetN()-1);
    graphNLODirGamma[1]->Draw("lp3");
    while(graphNLODirGamma[3]->GetX()[graphNLODirGamma[3]->GetN()-1] > 8.) graphNLODirGamma[3]->RemovePoint(graphNLODirGamma[3]->GetN()-1);
    graphNLODirGamma[3]->Draw("lp3");
    // while(graphNLODirGamma[0]->GetX()[graphNLODirGamma[0]->GetN()-1] > 4.) graphNLODirGamma[0]->RemovePoint(graphNLODirGamma[0]->GetN()-1);
    // graphNLODirGamma[0]->Draw("lp3");


    TArrow *ar8[50];
    for(Int_t i=1;i<graphDirectPhotonSpectrum[2]->GetN()+1;i++)
    {
        if (graphDirectPhotonSpectrum[2]->GetY()[i] > 5e-10){
            ar8[i] = new TArrow(graphDirectPhotonSpectrum[2]->GetX()[i],graphDirectPhotonSpectrum[2]->GetY()[i],graphDirectPhotonSpectrum[2]->GetX()[i],graphDirectPhotonSpectrum[2]->GetY()[i]*0.1,0.02,"|->");
            ar8[i]->SetAngle(40);
            ar8[i]->SetLineWidth(2);
            ar8[i]->SetLineColor(colorDataFULL[2]);
            ar8[i]->SetFillColor(colorDataFULL[2]);
            ar8[i]->Draw();
        }
    }

    TArrow *ar7[50];
    for(Int_t i=1;i<graphDirectPhotonSpectrum[1]->GetN()+1;i++)
    {
        if (graphDirectPhotonSpectrum[1]->GetY()[i] > 5e-10){
            ar7[i] = new TArrow(graphDirectPhotonSpectrum[1]->GetX()[i],graphDirectPhotonSpectrum[1]->GetY()[i],graphDirectPhotonSpectrum[1]->GetX()[i],graphDirectPhotonSpectrum[1]->GetY()[i]*0.1,0.02,"|->");
            ar7[i]->SetAngle(40);
            ar7[i]->SetLineWidth(2);
            ar7[i]->SetLineColor(colorDataFULL[1]);
            ar7[i]->SetFillColor(colorDataFULL[1]);
            ar7[i]->Draw();
        }
    }


    TArrow *ar276[50];
    for(Int_t i=1;i<graphDirectPhotonSpectrum[3]->GetN()+1;i++)
    {
        if (graphDirectPhotonSpectrum[3]->GetY()[i] > 5e-10){
            ar276[i] = new TArrow(graphDirectPhotonSpectrum[3]->GetX()[i],graphDirectPhotonSpectrum[3]->GetY()[i],graphDirectPhotonSpectrum[3]->GetX()[i],graphDirectPhotonSpectrum[3]->GetY()[i]*0.1,0.02,"|->");
            ar276[i]->SetAngle(40);
            ar276[i]->SetLineWidth(2);
            ar276[i]->SetLineColor(colorDataFULL[3]);
            ar276[i]->SetFillColor(colorDataFULL[3]);
            ar276[i]->Draw();
        }
    }


    TArrow *ar900[50];
    for(Int_t i=1;i<graphDirectPhotonSpectrum[0]->GetN()+1;i++)
    {
        if (graphDirectPhotonSpectrum[0]->GetY()[i] > 5e-10){
            ar900[i] = new TArrow(graphDirectPhotonSpectrum[0]->GetX()[i],graphDirectPhotonSpectrum[0]->GetY()[i],graphDirectPhotonSpectrum[0]->GetX()[i],graphDirectPhotonSpectrum[0]->GetY()[i]*0.1,0.02,"|->");
            ar900[i]->SetAngle(40);
            ar900[i]->SetLineWidth(2);
            ar900[i]->SetLineColor(colorDataFULL[0]);
            ar900[i]->SetFillColor(colorDataFULL[0]);
            ar900[i]->Draw();
        }
    }


    drawLatexAdd("ALICE",0.92,0.93,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE);
    // drawLatexAdd("#gamma's rec. with PCM",0.92,0.86,textsizeLabelsXSecUp,kFALSE,kFALSE,kTRUE);

    // TLegend* legendXsectionPaperSingle    = GetAndSetLegend2(0.17, 0.162, 0.5, 0.232+0.05*5, textSizeLabelsPixel);
    TLegend* legendXsectionPaperSingle    = GetAndSetLegend2(0.17, 0.102, 0.5, 0.172+0.045*3, textSizeLabelsPixel);
    legendXsectionPaperSingle->SetNColumns(1);
    legendXsectionPaperSingle->SetMargin(0.2);
    legendXsectionPaperSingle->AddEntry(graphDirectPhotonSpectrum[2],Form("%s (x 10^{6})",nameMeasGlobal[2].Data()),"pf");
    legendXsectionPaperSingle->AddEntry(graphDirectPhotonSpectrum[1],Form("%s (x 10^{4})",nameMeasGlobal[1].Data()),"pf");
    legendXsectionPaperSingle->AddEntry(graphDirectPhotonSpectrum[3],Form("%s (x 10^{2})",nameMeasGlobal[3].Data()),"pf");
    legendXsectionPaperSingle->AddEntry(graphDirectPhotonSpectrum[0],nameMeasGlobal[0].Data(),"pf");
    legendXsectionPaperSingle->AddEntry(graphNLODirGamma[1],"NLO prediction","l");

    legendXsectionPaperSingle->Draw();

    canvasInvSectionPaperSingle->Print(Form("%s/DirectPhotonSpectrumFinal.%s",outputDir.Data(),suffix.Data()));








}

void plotThreePads(){

    // **********************************************************************************************************************
    // ******************************************* Mass and width for pi0 ****************************************
    // **********************************************************************************************************************
    Double_t mesonMassExpectPi0                 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    Double_t mesonMassExpectEta                 = TDatabasePDG::Instance()->GetParticle(221)->Mass();


    Double_t arrayBoundariesX1_4[4];
    Double_t arrayBoundariesY1_4[3];
    Double_t relativeMarginsX[3];
    Double_t relativeMarginsY[3];
    Double_t maxX = 19.9;
    Double_t textSizeLabelsPixel             = 50;
    ReturnCorrectValuesForCanvasScaling(2350,1250, 3, 2,0.051, 0.005, 0.005,0.085,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);

    TCanvas* canvasMassWidthPi0     = new TCanvas("canvasMassWidthPi0","",0,0,2350,1250);  // gives the page size
    DrawGammaCanvasSettings( canvasMassWidthPi0,  0.13, 0.02, 0.03, 0.06);

    TPad* padWidthPi08TeV               = new TPad("padWidthPi08TeV", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padWidthPi08TeV, relativeMarginsX[0], relativeMarginsX[1], relativeMarginsY[0], relativeMarginsY[1]);
    padWidthPi08TeV->Draw();

    TPad* padMassPi08TeV                = new TPad("padMassPi08TeV", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
    DrawGammaPadSettings( padMassPi08TeV, relativeMarginsX[0], relativeMarginsX[1], relativeMarginsY[1], relativeMarginsY[2]);
    padMassPi08TeV->Draw();

    TPad* padWidthPi07TeV               = new TPad("padWidthPi07TeV", "", arrayBoundariesX1_4[1], arrayBoundariesY1_4[1], arrayBoundariesX1_4[2], arrayBoundariesY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padWidthPi07TeV, relativeMarginsX[1], relativeMarginsX[1], relativeMarginsY[0], relativeMarginsY[1]);
    padWidthPi07TeV->Draw();

    TPad* padMassPi07TeV                = new TPad("padMassPi07TeV", "", arrayBoundariesX1_4[1], arrayBoundariesY1_4[2], arrayBoundariesX1_4[2], arrayBoundariesY1_4[1],-1, -1, -2);
    DrawGammaPadSettings( padMassPi07TeV, relativeMarginsX[1], relativeMarginsX[1], relativeMarginsY[1], relativeMarginsY[2]);
    padMassPi07TeV->Draw();

    TPad* padWidthPi0900GeV               = new TPad("padWidthPi0900GeV", "", arrayBoundariesX1_4[2], arrayBoundariesY1_4[1], arrayBoundariesX1_4[3], arrayBoundariesY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padWidthPi0900GeV, relativeMarginsX[1], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[1]);
    padWidthPi0900GeV->Draw();

    TPad* padMassPi0900GeV                = new TPad("padMassPi0900GeV", "", arrayBoundariesX1_4[2], arrayBoundariesY1_4[2], arrayBoundariesX1_4[3], arrayBoundariesY1_4[1],-1, -1, -2);
    DrawGammaPadSettings( padMassPi0900GeV, relativeMarginsX[1], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[2]);
    padMassPi0900GeV->Draw();

    TPad* padMassLegend1            = new TPad("padMassLegend1", "", 0.13, 0.32, 0.52, 0.52,-1, -1, -2);
    DrawGammaPadSettings( padMassLegend1, 0., 0., 0., 0.);
    padMassLegend1->SetFillStyle(0);
    padMassLegend1->Draw();

    padWidthPi08TeV->cd();
    padWidthPi08TeV->SetLogx();

    Double_t margin                 = relativeMarginsX[0]*2.7*2350;
    Double_t textsizeLabelsWidth    = 0;
    Double_t textsizeFacWidth       = 0;
    if (padWidthPi08TeV->XtoPixel(padWidthPi08TeV->GetX2()) < padWidthPi08TeV->YtoPixel(padWidthPi08TeV->GetY1())){
        textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthPi08TeV->XtoPixel(padWidthPi08TeV->GetX2()) ;
        textsizeFacWidth            = (Double_t)1./padWidthPi08TeV->XtoPixel(padWidthPi08TeV->GetX2()) ;
    } else {
        textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthPi08TeV->YtoPixel(padWidthPi08TeV->GetY1());
        textsizeFacWidth            = (Double_t)1./padWidthPi08TeV->YtoPixel(padWidthPi08TeV->GetY1());
    }

    TH2F * histo2DAllPi0FWHM    = new TH2F("histo2DAllPi0FWHM","histo2DAllPi0FWHM", 20, 0.23, maxX ,1000., 0.1, 6.9);
    SetStyleHistoTH2ForGraphs(histo2DAllPi0FWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                            0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.77, 512, 505);
    //         histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(-1.,24.5);
    histo2DAllPi0FWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
    histo2DAllPi0FWHM->GetYaxis()->SetNdivisions(505);
    histo2DAllPi0FWHM->GetYaxis()->SetNoExponent(kTRUE);
    histo2DAllPi0FWHM->GetXaxis()->SetTickLength(0.05);
    histo2DAllPi0FWHM->GetYaxis()->SetTickLength(0.026);
    histo2DAllPi0FWHM->DrawCopy();




    DrawGammaSetMarker(histoPi0TrueFWHMMeV[0], markerstylesMC[0], markersizeMC[0], colorMC[0] , colorMC[0]);
    histoPi0TrueFWHMMeV[0]->Draw("p,same,e");
    DrawGammaSetMarker(histoPi0FWHMMeV[0], markerstyles[0], markersize[0], colorData[0] , colorData[0]);
    histoPi0FWHMMeV[0]->Draw("p,same,e");


    drawLatexAdd("ALICE",0.18,0.87,textSizeLabelsPixel,kFALSE,kTRUE);
    drawLatexAdd("pp, PCM",0.18,0.78,textSizeLabelsPixel,kFALSE,kTRUE);
    drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.18,0.69,textSizeLabelsPixel,kFALSE,kTRUE);
    // drawLatexAdd("a)",0.18,0.06,textSizeLabelsPixel,kFALSE,kTRUE);

    padWidthPi07TeV->cd();
    padWidthPi07TeV->SetLogx();
    histo2DAllPi0FWHM->DrawCopy();
    DrawGammaSetMarker(histoPi0TrueFWHMMeV[1], markerstylesMC[1], markersizeMC[1], colorMC[1] , colorMC[1]);
    histoPi0TrueFWHMMeV[1]->Draw("p,same,e");
    DrawGammaSetMarker(histoPi0FWHMMeV[1], markerstyles[1], markersize[1], colorData[1] , colorData[1]);
    histoPi0FWHMMeV[1]->Draw("p,same,e");

    padWidthPi0900GeV->cd();
    padWidthPi0900GeV->SetLogx();
    histo2DAllPi0FWHM->DrawCopy();
    DrawGammaSetMarker(histoPi0TrueFWHMMeV[2], markerstylesMC[2], markersizeMC[2], colorMC[2] , colorMC[2]);
    histoPi0TrueFWHMMeV[2]->Draw("p,same,e");
    DrawGammaSetMarker(histoPi0FWHMMeV[2], markerstyles[2], markersize[2], colorData[2] , colorData[2]);
    histoPi0FWHMMeV[2]->Draw("p,same,e");


    padMassPi08TeV->cd();
    padMassPi08TeV->SetLogx();


    Double_t textsizeLabelsMass         = 0;
    Double_t textsizeFacMass            = 0;
    if (padMassPi08TeV->XtoPixel(padMassPi08TeV->GetX2()) <padMassPi08TeV->YtoPixel(padMassPi08TeV->GetY1()) ){
        textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padMassPi08TeV->XtoPixel(padMassPi08TeV->GetX2()) ;
        textsizeFacMass                 = (Double_t)1./padMassPi08TeV->XtoPixel(padMassPi08TeV->GetX2()) ;
    } else {
        textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padMassPi08TeV->YtoPixel(padMassPi08TeV->GetY1());
        textsizeFacMass                 = (Double_t)1./padMassPi08TeV->YtoPixel(padMassPi08TeV->GetY1());
    }

    TH2F * histo2DAllPi0Mass            = new TH2F("histo2DAllPi0Mass","histo2DAllPi0Mass",20, 0.23, maxX, 1000., 132.1, 139.9);
    SetStyleHistoTH2ForGraphs(histo2DAllPi0Mass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
                            textsizeLabelsMass, 0.9, 0.9, 512, 505);
    histo2DAllPi0Mass->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo2DAllPi0Mass->GetYaxis()->SetNdivisions(505);
    histo2DAllPi0Mass->GetYaxis()->SetNoExponent(kTRUE);
    histo2DAllPi0Mass->GetXaxis()->SetNoExponent(kTRUE);
    histo2DAllPi0Mass->GetXaxis()->SetTickLength(0.05);
    //         histo2DAllPi0Mass->GetXaxis()->SetLabelOffset(-0.015);
    TH2F * histo2DAllPi0MassCopy = (TH2F*) histo2DAllPi0Mass->Clone("histo2DAllPi0MassCopy");
    histo2DAllPi0MassCopy->GetXaxis()->SetTitle("");
    histo2DAllPi0MassCopy->DrawCopy();
    // 900 GeV PAD
    DrawGammaSetMarker(histoPi0TrueMass[0], markerstylesMC[0], markersizeMC[0], colorMC[0] , colorMC[0]);
    histoPi0TrueMass[0]->Draw("p,same,e");
    DrawGammaSetMarker(histoPi0Mass[0], markerstyles[0], markersize[0], colorData[0] , colorData[0]);
    histoPi0Mass[0]->Draw("p,same,e");

    drawLatexAdd(nameMeasGlobal[0].Data(),0.92,0.87,textSizeLabelsPixel,kFALSE,kTRUE,kTRUE);
    // TLegend* legend900GeV    = GetAndSetLegend2(0.64, 0.67, 0.95, 0.85, textSizeLabelsPixel);
    TLegend* legend900GeV    = GetAndSetLegend2(0.20, 0.22, 0.46, 0.4, textSizeLabelsPixel);
    legend900GeV->SetNColumns(1);
    legend900GeV->SetMargin(0.2);
    legend900GeV->AddEntry(histoPi0Mass[0],"data","p");
    legend900GeV->AddEntry(histoPi0TrueMass[0],"MC","p");
    legend900GeV->Draw();

    DrawGammaLines(0.23, maxX , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1, kGray);
    // 7 TeV PAD
    padMassPi07TeV->cd();
    padMassPi07TeV->SetLogx();

    histo2DAllPi0MassCopy->DrawCopy();

    DrawGammaSetMarker(histoPi0TrueMass[1], markerstylesMC[1], markersizeMC[1], colorMC[1] , colorMC[1]);
    histoPi0TrueMass[1]->Draw("p,same,e");
    DrawGammaSetMarker(histoPi0Mass[1], markerstyles[1], markersize[1], colorData[1] , colorData[1]);
    histoPi0Mass[1]->Draw("p,same,e");

    drawLatexAdd(nameMeasGlobal[1].Data(),0.92,0.87,textSizeLabelsPixel,kFALSE,kTRUE,kTRUE);
    // TLegend* legend7TeV    = GetAndSetLegend2(0.64, 0.67, 0.90, 0.85, textSizeLabelsPixel);
    TLegend* legend7TeV    = GetAndSetLegend2(0.07, 0.22, 0.45, 0.4, textSizeLabelsPixel);
    legend7TeV->SetNColumns(1);
    legend7TeV->SetMargin(0.2);
    legend7TeV->AddEntry(histoPi0Mass[1],"data","p");
    legend7TeV->AddEntry(histoPi0TrueMass[1],"MC","p");
    legend7TeV->Draw();
    DrawGammaLines(0.23, maxX , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1, kGray);

    // 8 TeV PAD
    padMassPi0900GeV->cd();
    padMassPi0900GeV->SetLogx();
    histo2DAllPi0Mass->DrawCopy();

    DrawGammaSetMarker(histoPi0TrueMass[2], markerstylesMC[2], markersizeMC[2], colorMC[2] , colorMC[2]);
    histoPi0TrueMass[2]->Draw("p,same,e");
    DrawGammaSetMarker(histoPi0Mass[2], markerstyles[2], markersize[2], colorData[2] , colorData[2]);
    histoPi0Mass[2]->Draw("p,same,e");

    drawLatexAdd(nameMeasGlobal[2].Data(),0.92,0.87,textSizeLabelsPixel,kFALSE,kTRUE,kTRUE);
    // TLegend* legend8TeV    = GetAndSetLegend2(0.64, 0.67, 0.95, 0.85, textSizeLabelsPixel);
    TLegend* legend8TeV    = GetAndSetLegend2(0.05, 0.22, 0.45, 0.4, textSizeLabelsPixel);
    legend8TeV->SetNColumns(1);
    legend8TeV->SetMargin(0.2);
    legend8TeV->AddEntry(histoPi0Mass[2],"data","p");
    legend8TeV->AddEntry(histoPi0TrueMass[2],"MC","p");
    legend8TeV->Draw();

    DrawGammaLines(0.23, maxX , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1, kGray);


    canvasMassWidthPi0->Update();
    canvasMassWidthPi0->Print(Form("%s/Pi0_MassAndWidthSeparated.%s",outputDir.Data(),suffix.Data()));

}

void plotPi0FWHMMassNEW(){

    // **********************************************************************************************************************
    // ******************************************* Mass and width for pi0 ****************************************
    // **********************************************************************************************************************
    Double_t mesonMassExpectPi0                 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    Double_t mesonMassExpectEta                 = TDatabasePDG::Instance()->GetParticle(221)->Mass();


    Double_t arrayBoundariesX1_4[4];
    Double_t arrayBoundariesY1_4[3];
    Double_t relativeMarginsX[3];
    Double_t relativeMarginsY[3];
    Double_t maxX = 19.9;
    Double_t textSizeLabelsPixel             = 50;
    ReturnCorrectValuesForCanvasScaling(2350,1250, 3, 2,0.051, 0.005, 0.005,0.085,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);

    TCanvas* canvasMassWidthPi0     = new TCanvas("canvasMassWidthPi0","",0,0,2350,1250);  // gives the page size
    DrawGammaCanvasSettings( canvasMassWidthPi0,  0.13, 0.02, 0.03, 0.06);

    TPad* padWidthPi08TeV               = new TPad("padWidthPi08TeV", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padWidthPi08TeV, relativeMarginsX[0], relativeMarginsX[1], relativeMarginsY[0], relativeMarginsY[1]);
    padWidthPi08TeV->Draw();

    TPad* padMassPi08TeV                = new TPad("padMassPi08TeV", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
    DrawGammaPadSettings( padMassPi08TeV, relativeMarginsX[0], relativeMarginsX[1], relativeMarginsY[1], relativeMarginsY[2]);
    padMassPi08TeV->Draw();

    TPad* padWidthPi07TeV               = new TPad("padWidthPi07TeV", "", arrayBoundariesX1_4[1], arrayBoundariesY1_4[1], arrayBoundariesX1_4[2], arrayBoundariesY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padWidthPi07TeV, relativeMarginsX[1], relativeMarginsX[1], relativeMarginsY[0], relativeMarginsY[1]);
    padWidthPi07TeV->Draw();

    TPad* padMassPi07TeV                = new TPad("padMassPi07TeV", "", arrayBoundariesX1_4[1], arrayBoundariesY1_4[2], arrayBoundariesX1_4[2], arrayBoundariesY1_4[1],-1, -1, -2);
    DrawGammaPadSettings( padMassPi07TeV, relativeMarginsX[1], relativeMarginsX[1], relativeMarginsY[1], relativeMarginsY[2]);
    padMassPi07TeV->Draw();

    TPad* padWidthPi0900GeV               = new TPad("padWidthPi0900GeV", "", arrayBoundariesX1_4[2], arrayBoundariesY1_4[1], arrayBoundariesX1_4[3], arrayBoundariesY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padWidthPi0900GeV, relativeMarginsX[1], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[1]);
    padWidthPi0900GeV->Draw();

    TPad* padMassPi0900GeV                = new TPad("padMassPi0900GeV", "", arrayBoundariesX1_4[2], arrayBoundariesY1_4[2], arrayBoundariesX1_4[3], arrayBoundariesY1_4[1],-1, -1, -2);
    DrawGammaPadSettings( padMassPi0900GeV, relativeMarginsX[1], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[2]);
    padMassPi0900GeV->Draw();

    TPad* padMassLegend1            = new TPad("padMassLegend1", "", 0.13, 0.32, 0.52, 0.52,-1, -1, -2);
    DrawGammaPadSettings( padMassLegend1, 0., 0., 0., 0.);
    padMassLegend1->SetFillStyle(0);
    padMassLegend1->Draw();

    padWidthPi08TeV->cd();
    padWidthPi08TeV->SetLogx();

    Double_t margin                 = relativeMarginsX[0]*2.7*2350;
    Double_t textsizeLabelsWidth    = 0;
    Double_t textsizeFacWidth       = 0;
    if (padWidthPi08TeV->XtoPixel(padWidthPi08TeV->GetX2()) < padWidthPi08TeV->YtoPixel(padWidthPi08TeV->GetY1())){
        textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthPi08TeV->XtoPixel(padWidthPi08TeV->GetX2()) ;
        textsizeFacWidth            = (Double_t)1./padWidthPi08TeV->XtoPixel(padWidthPi08TeV->GetX2()) ;
    } else {
        textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthPi08TeV->YtoPixel(padWidthPi08TeV->GetY1());
        textsizeFacWidth            = (Double_t)1./padWidthPi08TeV->YtoPixel(padWidthPi08TeV->GetY1());
    }

    TH2F * histo2DAllPi0FWHM    = new TH2F("histo2DAllPi0FWHM","histo2DAllPi0FWHM", 20, 0.23, maxX ,1000., 0.1, 6.9);
    SetStyleHistoTH2ForGraphs(histo2DAllPi0FWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                            0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.77, 512, 505);
    //         histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(-1.,24.5);
    histo2DAllPi0FWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
    histo2DAllPi0FWHM->GetYaxis()->SetNdivisions(505);
    histo2DAllPi0FWHM->GetYaxis()->SetNoExponent(kTRUE);
    histo2DAllPi0FWHM->GetXaxis()->SetTickLength(0.05);
    histo2DAllPi0FWHM->GetYaxis()->SetTickLength(0.026);
    histo2DAllPi0FWHM->DrawCopy();




    DrawGammaSetMarker(histoPi0TrueFWHMMeV[0], markerstylesMC[0], markersizeMC[0], colorMC[0] , colorMC[0]);
    histoPi0TrueFWHMMeV[0]->Draw("p,same,e");
    DrawGammaSetMarker(histoPi0FWHMMeV[0], markerstyles[0], markersize[0], colorData[0] , colorData[0]);
    histoPi0FWHMMeV[0]->Draw("p,same,e");


    drawLatexAdd("ALICE",0.18,0.87,textSizeLabelsPixel,kFALSE,kTRUE);
    drawLatexAdd("pp, PCM",0.18,0.78,textSizeLabelsPixel,kFALSE,kTRUE);
    drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.18,0.69,textSizeLabelsPixel,kFALSE,kTRUE);
    // drawLatexAdd("a)",0.18,0.06,textSizeLabelsPixel,kFALSE,kTRUE);

    padWidthPi07TeV->cd();
    padWidthPi07TeV->SetLogx();
    histo2DAllPi0FWHM->DrawCopy();
    DrawGammaSetMarker(histoPi0TrueFWHMMeV[1], markerstylesMC[1], markersizeMC[1], colorMC[1] , colorMC[1]);
    histoPi0TrueFWHMMeV[1]->Draw("p,same,e");
    DrawGammaSetMarker(histoPi0FWHMMeV[1], markerstyles[1], markersize[1], colorData[1] , colorData[1]);
    histoPi0FWHMMeV[1]->Draw("p,same,e");

    padWidthPi0900GeV->cd();
    padWidthPi0900GeV->SetLogx();
    histo2DAllPi0FWHM->DrawCopy();
    DrawGammaSetMarker(histoPi0TrueFWHMMeV[2], markerstylesMC[2], markersizeMC[2], colorMC[2] , colorMC[2]);
    histoPi0TrueFWHMMeV[2]->Draw("p,same,e");
    DrawGammaSetMarker(histoPi0FWHMMeV[2], markerstyles[2], markersize[2], colorData[2] , colorData[2]);
    histoPi0FWHMMeV[2]->Draw("p,same,e");


    padMassPi08TeV->cd();
    padMassPi08TeV->SetLogx();


    Double_t textsizeLabelsMass         = 0;
    Double_t textsizeFacMass            = 0;
    if (padMassPi08TeV->XtoPixel(padMassPi08TeV->GetX2()) <padMassPi08TeV->YtoPixel(padMassPi08TeV->GetY1()) ){
        textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padMassPi08TeV->XtoPixel(padMassPi08TeV->GetX2()) ;
        textsizeFacMass                 = (Double_t)1./padMassPi08TeV->XtoPixel(padMassPi08TeV->GetX2()) ;
    } else {
        textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padMassPi08TeV->YtoPixel(padMassPi08TeV->GetY1());
        textsizeFacMass                 = (Double_t)1./padMassPi08TeV->YtoPixel(padMassPi08TeV->GetY1());
    }

    TH2F * histo2DAllPi0Mass            = new TH2F("histo2DAllPi0Mass","histo2DAllPi0Mass",20, 0.23, maxX, 1000., 132.1, 139.9);
    SetStyleHistoTH2ForGraphs(histo2DAllPi0Mass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
                            textsizeLabelsMass, 0.9, 0.9, 512, 505);
    histo2DAllPi0Mass->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo2DAllPi0Mass->GetYaxis()->SetNdivisions(505);
    histo2DAllPi0Mass->GetYaxis()->SetNoExponent(kTRUE);
    histo2DAllPi0Mass->GetXaxis()->SetNoExponent(kTRUE);
    histo2DAllPi0Mass->GetXaxis()->SetTickLength(0.05);
    //         histo2DAllPi0Mass->GetXaxis()->SetLabelOffset(-0.015);
    TH2F * histo2DAllPi0MassCopy = (TH2F*) histo2DAllPi0Mass->Clone("histo2DAllPi0MassCopy");
    histo2DAllPi0MassCopy->GetXaxis()->SetTitle("");
    histo2DAllPi0MassCopy->DrawCopy();
    // 900 GeV PAD
    DrawGammaSetMarker(histoPi0TrueMass[0], markerstylesMC[0], markersizeMC[0], colorMC[0] , colorMC[0]);
    histoPi0TrueMass[0]->Draw("p,same,e");
    DrawGammaSetMarker(histoPi0Mass[0], markerstyles[0], markersize[0], colorData[0] , colorData[0]);
    histoPi0Mass[0]->Draw("p,same,e");

    drawLatexAdd(nameMeasGlobal[0].Data(),0.92,0.87,textSizeLabelsPixel,kFALSE,kTRUE,kTRUE);
    // TLegend* legend900GeV    = GetAndSetLegend2(0.64, 0.67, 0.95, 0.85, textSizeLabelsPixel);
    TLegend* legend900GeV    = GetAndSetLegend2(0.20, 0.22, 0.46, 0.4, textSizeLabelsPixel);
    legend900GeV->SetNColumns(1);
    legend900GeV->SetMargin(0.2);
    legend900GeV->AddEntry(histoPi0Mass[0],"data","p");
    legend900GeV->AddEntry(histoPi0TrueMass[0],"MC","p");
    legend900GeV->Draw();

    DrawGammaLines(0.23, maxX , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1, kGray);
    // 7 TeV PAD
    padMassPi07TeV->cd();
    padMassPi07TeV->SetLogx();

    histo2DAllPi0MassCopy->DrawCopy();

    DrawGammaSetMarker(histoPi0TrueMass[1], markerstylesMC[1], markersizeMC[1], colorMC[1] , colorMC[1]);
    histoPi0TrueMass[1]->Draw("p,same,e");
    DrawGammaSetMarker(histoPi0Mass[1], markerstyles[1], markersize[1], colorData[1] , colorData[1]);
    histoPi0Mass[1]->Draw("p,same,e");

    drawLatexAdd(nameMeasGlobal[1].Data(),0.92,0.87,textSizeLabelsPixel,kFALSE,kTRUE,kTRUE);
    // TLegend* legend7TeV    = GetAndSetLegend2(0.64, 0.67, 0.90, 0.85, textSizeLabelsPixel);
    TLegend* legend7TeV    = GetAndSetLegend2(0.07, 0.22, 0.45, 0.4, textSizeLabelsPixel);
    legend7TeV->SetNColumns(1);
    legend7TeV->SetMargin(0.2);
    legend7TeV->AddEntry(histoPi0Mass[1],"data","p");
    legend7TeV->AddEntry(histoPi0TrueMass[1],"MC","p");
    legend7TeV->Draw();
    DrawGammaLines(0.23, maxX , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1, kGray);

    // 8 TeV PAD
    padMassPi0900GeV->cd();
    padMassPi0900GeV->SetLogx();
    histo2DAllPi0Mass->DrawCopy();

    DrawGammaSetMarker(histoPi0TrueMass[2], markerstylesMC[2], markersizeMC[2], colorMC[2] , colorMC[2]);
    histoPi0TrueMass[2]->Draw("p,same,e");
    DrawGammaSetMarker(histoPi0Mass[2], markerstyles[2], markersize[2], colorData[2] , colorData[2]);
    histoPi0Mass[2]->Draw("p,same,e");

    drawLatexAdd(nameMeasGlobal[2].Data(),0.92,0.87,textSizeLabelsPixel,kFALSE,kTRUE,kTRUE);
    // TLegend* legend8TeV    = GetAndSetLegend2(0.64, 0.67, 0.95, 0.85, textSizeLabelsPixel);
    TLegend* legend8TeV    = GetAndSetLegend2(0.05, 0.22, 0.45, 0.4, textSizeLabelsPixel);
    legend8TeV->SetNColumns(1);
    legend8TeV->SetMargin(0.2);
    legend8TeV->AddEntry(histoPi0Mass[2],"data","p");
    legend8TeV->AddEntry(histoPi0TrueMass[2],"MC","p");
    legend8TeV->Draw();

    DrawGammaLines(0.23, maxX , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1, kGray);


    canvasMassWidthPi0->Update();
    canvasMassWidthPi0->Print(Form("%s/Pi0_MassAndWidthSeparated.%s",outputDir.Data(),suffix.Data()));

}

void plotEtaFWHMMassNEW(){

    // **********************************************************************************************************************
    // ******************************************* Mass and width for pi0 ****************************************
    // **********************************************************************************************************************
    Double_t mesonMassExpectPi0                 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    Double_t mesonMassExpectEta                 = TDatabasePDG::Instance()->GetParticle(221)->Mass();


    Double_t arrayBoundariesX1_4[4];
    Double_t arrayBoundariesY1_4[3];
    Double_t relativeMarginsX[3];
    Double_t relativeMarginsY[3];
    Double_t minX = 0.33;
    Double_t maxX = 12.9;
    Double_t textSizeLabelsPixel             = 50;
    ReturnCorrectValuesForCanvasScaling(2350,1250, 3, 2,0.051, 0.005, 0.005,0.085,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);

    TCanvas* canvasMassWidthEta     = new TCanvas("canvasMassWidthEta","",0,0,2350,1250);  // gives the page size
    DrawGammaCanvasSettings( canvasMassWidthEta,  0.13, 0.02, 0.03, 0.06);

    TPad* padWidthEta8TeV               = new TPad("padWidthEta8TeV", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padWidthEta8TeV, relativeMarginsX[0], relativeMarginsX[1], relativeMarginsY[0], relativeMarginsY[1]);
    padWidthEta8TeV->Draw();

    TPad* padMassEta8TeV                = new TPad("padMassEta8TeV", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
    DrawGammaPadSettings( padMassEta8TeV, relativeMarginsX[0], relativeMarginsX[1], relativeMarginsY[1], relativeMarginsY[2]);
    padMassEta8TeV->Draw();

    TPad* padWidthEta7TeV               = new TPad("padWidthEta7TeV", "", arrayBoundariesX1_4[1], arrayBoundariesY1_4[1], arrayBoundariesX1_4[2], arrayBoundariesY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padWidthEta7TeV, relativeMarginsX[1], relativeMarginsX[1], relativeMarginsY[0], relativeMarginsY[1]);
    padWidthEta7TeV->Draw();

    TPad* padMassEta7TeV                = new TPad("padMassEta7TeV", "", arrayBoundariesX1_4[1], arrayBoundariesY1_4[2], arrayBoundariesX1_4[2], arrayBoundariesY1_4[1],-1, -1, -2);
    DrawGammaPadSettings( padMassEta7TeV, relativeMarginsX[1], relativeMarginsX[1], relativeMarginsY[1], relativeMarginsY[2]);
    padMassEta7TeV->Draw();

    TPad* padWidthEta900GeV               = new TPad("padWidthEta900GeV", "", arrayBoundariesX1_4[2], arrayBoundariesY1_4[1], arrayBoundariesX1_4[3], arrayBoundariesY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padWidthEta900GeV, relativeMarginsX[1], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[1]);
    padWidthEta900GeV->Draw();

    TPad* padMassEta900GeV                = new TPad("padMassEta900GeV", "", arrayBoundariesX1_4[2], arrayBoundariesY1_4[2], arrayBoundariesX1_4[3], arrayBoundariesY1_4[1],-1, -1, -2);
    DrawGammaPadSettings( padMassEta900GeV, relativeMarginsX[1], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[2]);
    padMassEta900GeV->Draw();

    TPad* padMassLegend1            = new TPad("padMassLegend1", "", 0.13, 0.32, 0.52, 0.52,-1, -1, -2);
    DrawGammaPadSettings( padMassLegend1, 0., 0., 0., 0.);
    padMassLegend1->SetFillStyle(0);
    padMassLegend1->Draw();

    padWidthEta8TeV->cd();
    padWidthEta8TeV->SetLogx();

    Double_t margin                 = relativeMarginsX[0]*2.7*2350;
    Double_t textsizeLabelsWidth    = 0;
    Double_t textsizeFacWidth       = 0;
    if (padWidthEta8TeV->XtoPixel(padWidthEta8TeV->GetX2()) < padWidthEta8TeV->YtoPixel(padWidthEta8TeV->GetY1())){
        textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthEta8TeV->XtoPixel(padWidthEta8TeV->GetX2()) ;
        textsizeFacWidth            = (Double_t)1./padWidthEta8TeV->XtoPixel(padWidthEta8TeV->GetX2()) ;
    } else {
        textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthEta8TeV->YtoPixel(padWidthEta8TeV->GetY1());
        textsizeFacWidth            = (Double_t)1./padWidthEta8TeV->YtoPixel(padWidthEta8TeV->GetY1());
    }

    TH2F * histo2DAllEtaFWHM    = new TH2F("histo2DAllEtaFWHM","histo2DAllEtaFWHM", 20, minX, maxX ,1000., 0.01, 13.9);
    SetStyleHistoTH2ForGraphs(histo2DAllEtaFWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                            0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.77, 512, 505);
    //         histo2DAllEtaFWHM->GetYaxis()->SetRangeUser(-1.,24.5);
    histo2DAllEtaFWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
    histo2DAllEtaFWHM->GetYaxis()->SetNdivisions(505);
    histo2DAllEtaFWHM->GetYaxis()->SetNoExponent(kTRUE);
    histo2DAllEtaFWHM->GetXaxis()->SetTickLength(0.05);
    histo2DAllEtaFWHM->GetYaxis()->SetTickLength(0.026);
    histo2DAllEtaFWHM->DrawCopy();




    DrawGammaSetMarker(histoEtaTrueFWHMMeV[0], markerstylesMC[0], markersizeMC[0], colorMC[0] , colorMC[0]);
    histoEtaTrueFWHMMeV[0]->Draw("p,same,e");
    DrawGammaSetMarker(histoEtaFWHMMeV[0], markerstyles[0], markersize[0], colorData[0] , colorData[0]);
    histoEtaFWHMMeV[0]->Draw("p,same,e");


    drawLatexAdd("ALICE",0.18,0.87,textSizeLabelsPixel,kFALSE,kTRUE);
    drawLatexAdd("pp, PCM",0.18,0.78,textSizeLabelsPixel,kFALSE,kTRUE);
    drawLatexAdd("#eta #rightarrow #gamma#gamma",0.18,0.69,textSizeLabelsPixel,kFALSE,kTRUE);
    // drawLatexAdd("a)",0.18,0.06,textSizeLabelsPixel,kFALSE,kTRUE);

    padWidthEta7TeV->cd();
    padWidthEta7TeV->SetLogx();
    histo2DAllEtaFWHM->DrawCopy();
    DrawGammaSetMarker(histoEtaTrueFWHMMeV[1], markerstylesMC[1], markersizeMC[1], colorMC[1] , colorMC[1]);
    histoEtaTrueFWHMMeV[1]->Draw("p,same,e");
    DrawGammaSetMarker(histoEtaFWHMMeV[1], markerstyles[1], markersize[1], colorData[1] , colorData[1]);
    histoEtaFWHMMeV[1]->Draw("p,same,e");

    padWidthEta900GeV->cd();
    padWidthEta900GeV->SetLogx();
    histo2DAllEtaFWHM->DrawCopy();
    DrawGammaSetMarker(histoEtaTrueFWHMMeV[2], markerstylesMC[2], markersizeMC[2], colorMC[2] , colorMC[2]);
    histoEtaTrueFWHMMeV[2]->Draw("p,same,e");
    DrawGammaSetMarker(histoEtaFWHMMeV[2], markerstyles[2], markersize[2], colorData[2] , colorData[2]);
    histoEtaFWHMMeV[2]->Draw("p,same,e");


    padMassEta8TeV->cd();
    padMassEta8TeV->SetLogx();


    Double_t textsizeLabelsMass         = 0;
    Double_t textsizeFacMass            = 0;
    if (padMassEta8TeV->XtoPixel(padMassEta8TeV->GetX2()) <padMassEta8TeV->YtoPixel(padMassEta8TeV->GetY1()) ){
        textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padMassEta8TeV->XtoPixel(padMassEta8TeV->GetX2()) ;
        textsizeFacMass                 = (Double_t)1./padMassEta8TeV->XtoPixel(padMassEta8TeV->GetX2()) ;
    } else {
        textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padMassEta8TeV->YtoPixel(padMassEta8TeV->GetY1());
        textsizeFacMass                 = (Double_t)1./padMassEta8TeV->YtoPixel(padMassEta8TeV->GetY1());
    }

    TH2F * histo2DAllEtaMass            = new TH2F("histo2DAllEtaMass","histo2DAllEtaMass",20, minX, maxX, 1000., 537.9, 559.9);
    SetStyleHistoTH2ForGraphs(histo2DAllEtaMass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
                            textsizeLabelsMass, 0.9, 0.9, 512, 505);
    histo2DAllEtaMass->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo2DAllEtaMass->GetYaxis()->SetNdivisions(505);
    histo2DAllEtaMass->GetYaxis()->SetNoExponent(kTRUE);
    histo2DAllEtaMass->GetXaxis()->SetNoExponent(kTRUE);
    histo2DAllEtaMass->GetXaxis()->SetTickLength(0.05);
    //         histo2DAllEtaMass->GetXaxis()->SetLabelOffset(-0.015);
    TH2F * histo2DAllEtaMassCopy = (TH2F*) histo2DAllEtaMass->Clone("histo2DAllEtaMassCopy");
    histo2DAllEtaMassCopy->GetXaxis()->SetTitle("");
    histo2DAllEtaMassCopy->DrawCopy();
    // 900 GeV PAD
    DrawGammaSetMarker(histoEtaTrueMass[0], markerstylesMC[0], markersizeMC[0], colorMC[0] , colorMC[0]);
    histoEtaTrueMass[0]->Draw("p,same,e");
    DrawGammaSetMarker(histoEtaMass[0], markerstyles[0], markersize[0], colorData[0] , colorData[0]);
    histoEtaMass[0]->Draw("p,same,e");

    drawLatexAdd(nameMeasGlobal[0].Data(),0.92,0.87,textSizeLabelsPixel,kFALSE,kTRUE,kTRUE);
    // TLegend* legend900GeV    = GetAndSetLegend2(0.64, 0.67, 0.95, 0.85, textSizeLabelsPixel);
    TLegend* legend900GeV    = GetAndSetLegend2(0.20, 0.22, 0.46, 0.4, textSizeLabelsPixel);
    legend900GeV->SetNColumns(1);
    legend900GeV->SetMargin(0.2);
    legend900GeV->AddEntry(histoEtaMass[0],"data","p");
    legend900GeV->AddEntry(histoEtaTrueMass[0],"MC","p");
    legend900GeV->Draw();

    DrawGammaLines(minX, maxX , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,0.1, kGray);
    // 7 TeV PAD
    padMassEta7TeV->cd();
    padMassEta7TeV->SetLogx();

    histo2DAllEtaMassCopy->DrawCopy();

    DrawGammaSetMarker(histoEtaTrueMass[1], markerstylesMC[1], markersizeMC[1], colorMC[1] , colorMC[1]);
    histoEtaTrueMass[1]->Draw("p,same,e");
    DrawGammaSetMarker(histoEtaMass[1], markerstyles[1], markersize[1], colorData[1] , colorData[1]);
    histoEtaMass[1]->Draw("p,same,e");

    drawLatexAdd(nameMeasGlobal[1].Data(),0.92,0.87,textSizeLabelsPixel,kFALSE,kTRUE,kTRUE);
    // TLegend* legend7TeV    = GetAndSetLegend2(0.64, 0.67, 0.90, 0.85, textSizeLabelsPixel);
    TLegend* legend7TeV    = GetAndSetLegend2(0.07, 0.22, 0.45, 0.4, textSizeLabelsPixel);
    legend7TeV->SetNColumns(1);
    legend7TeV->SetMargin(0.2);
    legend7TeV->AddEntry(histoEtaMass[1],"data","p");
    legend7TeV->AddEntry(histoEtaTrueMass[1],"MC","p");
    legend7TeV->Draw();
    DrawGammaLines(minX, maxX , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,0.1, kGray);

    // 8 TeV PAD
    padMassEta900GeV->cd();
    padMassEta900GeV->SetLogx();
    histo2DAllEtaMass->DrawCopy();

    DrawGammaSetMarker(histoEtaTrueMass[2], markerstylesMC[2], markersizeMC[2], colorMC[2] , colorMC[2]);
    histoEtaTrueMass[2]->Draw("p,same,e");
    DrawGammaSetMarker(histoEtaMass[2], markerstyles[2], markersize[2], colorData[2] , colorData[2]);
    histoEtaMass[2]->Draw("p,same,e");

    drawLatexAdd(nameMeasGlobal[2].Data(),0.92,0.87,textSizeLabelsPixel,kFALSE,kTRUE,kTRUE);
    // TLegend* legend8TeV    = GetAndSetLegend2(0.64, 0.67, 0.95, 0.85, textSizeLabelsPixel);
    TLegend* legend8TeV    = GetAndSetLegend2(0.05, 0.22, 0.45, 0.4, textSizeLabelsPixel);
    legend8TeV->SetNColumns(1);
    legend8TeV->SetMargin(0.2);
    legend8TeV->AddEntry(histoEtaMass[2],"data","p");
    legend8TeV->AddEntry(histoEtaTrueMass[2],"MC","p");
    legend8TeV->Draw();

    DrawGammaLines(minX, maxX , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,0.1, kGray);


    canvasMassWidthEta->Update();
    canvasMassWidthEta->Print(Form("%s/Eta_MassAndWidthSeparated.%s",outputDir.Data(),suffix.Data()));

}

void plotPi0FWHMMass(){

     // **********************************************************************************************************************
    // ******************************************* Mass and width for pi0 ****************************************
    // **********************************************************************************************************************
    Double_t mesonMassExpectPi0                 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    Double_t mesonMassExpectEta                 = TDatabasePDG::Instance()->GetParticle(221)->Mass();


    Double_t arrayBoundariesX1_4[2];
    Double_t arrayBoundariesY1_4[3];
    Double_t relativeMarginsX[3];
    Double_t relativeMarginsY[3];
    Double_t textSizeLabelsPixel             = 50;
    ReturnCorrectValuesForCanvasScaling(1350,1250, 1, 2,0.09, 0.005, 0.005,0.085,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);

    TCanvas* canvasMassWidthPi0     = new TCanvas("canvasMassWidthPi0","",0,0,1350,1250);  // gives the page size
    DrawGammaCanvasSettings( canvasMassWidthPi0,  0.13, 0.02, 0.03, 0.06);

    TPad* padWidthPi0               = new TPad("padWidthPi0", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padWidthPi0, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[1]);
    padWidthPi0->Draw();

    TPad* padMassPi0                = new TPad("padMassPi0", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
    DrawGammaPadSettings( padMassPi0, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[2]);
    padMassPi0->Draw();

    TPad* padMassLegend1            = new TPad("padMassLegend1", "", 0.13, 0.32, 0.52, 0.52,-1, -1, -2);
    DrawGammaPadSettings( padMassLegend1, 0., 0., 0., 0.);
    padMassLegend1->SetFillStyle(0);
    padMassLegend1->Draw();

    padWidthPi0->cd();
    padWidthPi0->SetLogx();

        Double_t margin                 = relativeMarginsX[0]*2.7*1350;
        Double_t textsizeLabelsWidth    = 0;
        Double_t textsizeFacWidth       = 0;
        if (padWidthPi0->XtoPixel(padWidthPi0->GetX2()) < padWidthPi0->YtoPixel(padWidthPi0->GetY1())){
            textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthPi0->XtoPixel(padWidthPi0->GetX2()) ;
            textsizeFacWidth            = (Double_t)1./padWidthPi0->XtoPixel(padWidthPi0->GetX2()) ;
        } else {
            textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthPi0->YtoPixel(padWidthPi0->GetY1());
            textsizeFacWidth            = (Double_t)1./padWidthPi0->YtoPixel(padWidthPi0->GetY1());
        }

        TH2F * histo2DAllPi0FWHM    = new TH2F("histo2DAllPi0FWHM","histo2DAllPi0FWHM", 20, 0.23, 20. ,1000., 0.1, 9.9);
        SetStyleHistoTH2ForGraphs(histo2DAllPi0FWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                                  0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin), 512, 505);
        // histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(-1.,24.5);
        histo2DAllPi0FWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
        histo2DAllPi0FWHM->GetYaxis()->SetNdivisions(505);
        histo2DAllPi0FWHM->GetYaxis()->SetNoExponent(kTRUE);
        histo2DAllPi0FWHM->GetXaxis()->SetTickLength(0.05);
        histo2DAllPi0FWHM->GetYaxis()->SetTickLength(0.026);
        histo2DAllPi0FWHM->DrawCopy();

                DrawGammaSetMarker(histoPi0TrueFWHMMeV[1], markerstylesMC[1], markersizeMC[1], colorMC[1] , colorMC[1]);
                histoPi0TrueFWHMMeV[1]->Draw("p,same,e");
                DrawGammaSetMarker(histoPi0FWHMMeV[1], markerstyles[1], markersize[1], colorData[1] , colorData[1]);
                histoPi0FWHMMeV[1]->Draw("p,same,e");

                DrawGammaSetMarker(histoPi0TrueFWHMMeV[2], markerstylesMC[2], markersizeMC[2], colorMC[2] , colorMC[2]);
                histoPi0TrueFWHMMeV[2]->Draw("p,same,e");
                DrawGammaSetMarker(histoPi0FWHMMeV[2], markerstyles[2], markersize[2], colorData[2] , colorData[2]);
                histoPi0FWHMMeV[2]->Draw("p,same,e");

                DrawGammaSetMarker(histoPi0TrueFWHMMeV[0], markerstylesMC[0], markersizeMC[0], colorMC[0] , colorMC[0]);
                histoPi0TrueFWHMMeV[0]->Draw("p,same,e");
                DrawGammaSetMarker(histoPi0FWHMMeV[0], markerstyles[0], markersize[0], colorData[0] , colorData[0]);
                histoPi0FWHMMeV[0]->Draw("p,same,e");

        TLatex *labelLegendAMass    = new TLatex(0.13,0.06,"a)");
        SetStyleTLatex( labelLegendAMass, textSizeLabelsPixel,4);
        labelLegendAMass->SetTextFont(43);
        labelLegendAMass->Draw();

        TLatex *labelMassPerf       = new TLatex(0.13,0.87,"ALICE");
        SetStyleTLatex( labelMassPerf, textSizeLabelsPixel,4);
        labelMassPerf->SetTextFont(43);
        labelMassPerf->Draw();
        TLatex *labelMassEnergy     = new TLatex(0.13,0.78,"pp, PCM");
        SetStyleTLatex( labelMassEnergy, textSizeLabelsPixel,4);
        labelMassEnergy->SetTextFont(43);
        labelMassEnergy->Draw();
        TLatex *labelMassPi0        = new TLatex(0.13,0.69,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelMassPi0, textSizeLabelsPixel,4);
        labelMassPi0->SetTextFont(43);
        labelMassPi0->Draw();

    padMassPi0->cd();
    padMassPi0->SetLogx();

        Double_t textsizeLabelsMass         = 0;
        Double_t textsizeFacMass            = 0;
        if (padMassPi0->XtoPixel(padMassPi0->GetX2()) <padMassPi0->YtoPixel(padMassPi0->GetY1()) ){
            textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padMassPi0->XtoPixel(padMassPi0->GetX2()) ;
            textsizeFacMass                 = (Double_t)1./padMassPi0->XtoPixel(padMassPi0->GetX2()) ;
        } else {
            textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padMassPi0->YtoPixel(padMassPi0->GetY1());
            textsizeFacMass                 = (Double_t)1./padMassPi0->YtoPixel(padMassPi0->GetY1());
        }

        TH2F * histo2DAllPi0Mass            = new TH2F("histo2DAllPi0Mass","histo2DAllPi0Mass",20, 0.23, 20., 1000., 130.1, 143.9);
        SetStyleHistoTH2ForGraphs(histo2DAllPi0Mass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
                                  textsizeLabelsMass, 0.9, 0.28/(textsizeFacMass*margin), 512, 505);
        histo2DAllPi0Mass->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DAllPi0Mass->GetYaxis()->SetNdivisions(505);
        histo2DAllPi0Mass->GetYaxis()->SetNoExponent(kTRUE);
        histo2DAllPi0Mass->GetXaxis()->SetTickLength(0.05);
        histo2DAllPi0Mass->GetXaxis()->SetLabelOffset(-0.015);
        histo2DAllPi0Mass->DrawCopy();

                DrawGammaSetMarker(histoPi0TrueMass[1], markerstylesMC[1], markersizeMC[1], colorMC[1] , colorMC[1]);
                histoPi0TrueMass[1]->Draw("p,same,e");
                DrawGammaSetMarker(histoPi0Mass[1], markerstyles[1], markersize[1], colorData[1] , colorData[1]);
                histoPi0Mass[1]->Draw("p,same,e");

                DrawGammaSetMarker(histoPi0TrueMass[2], markerstylesMC[2], markersizeMC[2], colorMC[2] , colorMC[2]);
                histoPi0TrueMass[2]->Draw("p,same,e");
                DrawGammaSetMarker(histoPi0Mass[2], markerstyles[2], markersize[2], colorData[2] , colorData[2]);
                histoPi0Mass[2]->Draw("p,same,e");

                DrawGammaSetMarker(histoPi0TrueMass[0], markerstylesMC[0], markersizeMC[0], colorMC[0] , colorMC[0]);
                histoPi0TrueMass[0]->Draw("p,same,e");
                DrawGammaSetMarker(histoPi0Mass[0], markerstyles[0], markersize[0], colorData[0] , colorData[0]);
                histoPi0Mass[0]->Draw("p,same,e");


        DrawGammaLines(0.23, 20. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1, kGray);

        TLatex *labelLegendBMass            = new TLatex(0.13,0.22,"b)");
        SetStyleTLatex( labelLegendBMass, textSizeLabelsPixel,4);
        labelLegendBMass->SetTextFont(43);
        labelLegendBMass->Draw();

        //********************************** Defintion of the Legend **************************************************
        Double_t columnsLegendMass2[3]      = {0.,0.57,0.84};
        // Double_t rowsLegendMass2[5] = {0.8,0.6,0.4,0.2,0.01};
        // Double_t rowsLegendMass2[6] = {0.84,0.66,0.50,0.33,0.16,0.01};
          Double_t  rowsLegendMass2[7]= {0.84,0.66,0.47,0.28,0.01,0.16};
        //******************* Offsets ***********************
        Double_t offsetMarkerXMass2         = 0.1;
        Double_t offsetMarkerYMass2         = 0.1;
        //****************** Scale factors ******************
        Double_t scaleMarkerMass2           = 1.2;

        padMassLegend1->cd();
        //****************** first Column **************************************************
        TLatex *textMassPCM[10];
        for (Int_t i = 0; i < 3; i++){
                textMassPCM[i]                  = new TLatex(columnsLegendMass2[0],rowsLegendMass2[i+1],nameMeasGlobal[i].Data());
                SetStyleTLatex( textMassPCM[i], textSizeLabelsPixel,4);
                textMassPCM[i]->SetTextFont(43);
                textMassPCM[i]->Draw();
        }
        //****************** second Column *************************************************
        TLatex *textMassData                = new TLatex(columnsLegendMass2[1],rowsLegendMass2[0] ,"Data");
        SetStyleTLatex( textMassData, textSizeLabelsPixel,4);
        textMassData->SetTextFont(43);
        textMassData->Draw();
        TLatex *textMassMC                  = new TLatex(columnsLegendMass2[2] ,rowsLegendMass2[0],"MC");
        SetStyleTLatex( textMassMC, textSizeLabelsPixel,4);
        textMassMC->SetTextFont(43);
        textMassMC->Draw();

        TMarker* markerPCMPi0Mass[10];
        TMarker* markerPCMPi0MassMC[10];
        for (Int_t i = 0; i < 3; i++){
                markerPCMPi0Mass[i]             = CreateMarkerFromHisto(histoPi0Mass[i],columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
                markerPCMPi0Mass[i]->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2);
                markerPCMPi0MassMC[i]           = CreateMarkerFromHisto(histoPi0TrueMass[i],columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
                markerPCMPi0MassMC[i]->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2);
        }

    canvasMassWidthPi0->Update();
    canvasMassWidthPi0->Print(Form("%s/Pi0_MassAndWidth.%s",outputDir.Data(),suffix.Data()));

}

void plotEtaFWHMMass(){

    // **********************************************************************************************************************
    // ******************************************* Mass and width for Eta at 8TeV ****************************************
    // **********************************************************************************************************************
    Double_t mesonMassExpectPi0                 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    Double_t mesonMassExpectEta                 = TDatabasePDG::Instance()->GetParticle(221)->Mass();

    Double_t arrayBoundariesX1_4[2];
    Double_t arrayBoundariesY1_4[3];
    Double_t relativeMarginsX[3];
    Double_t relativeMarginsY[3];
    Double_t textSizeLabelsPixel             = 50;
    ReturnCorrectValuesForCanvasScaling(1350,1250, 1, 2,0.09, 0.005, 0.005,0.085,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);

    TCanvas* canvasMassWidthEta     = new TCanvas("canvasMassWidthEta","",0,0,1350,1250);  // gives the page size
    DrawGammaCanvasSettings( canvasMassWidthEta,  0.13, 0.02, 0.03, 0.06);

    TPad* padWidthEta               = new TPad("padWidthEta", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padWidthEta, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[1]);
    padWidthEta->Draw();

    TPad* padMassEta                = new TPad("padMassEta", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
    DrawGammaPadSettings( padMassEta, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[2]);
    padMassEta->Draw();

    TPad* padMassLegend1Eta            = new TPad("padMassLegend1Eta", "", 0.13, 0.34, 0.52, 0.52,-1, -1, -2);
    DrawGammaPadSettings( padMassLegend1Eta, 0., 0., 0., 0.);
    padMassLegend1Eta->SetFillStyle(0);
    padMassLegend1Eta->Draw();

    padWidthEta->cd();
    padWidthEta->SetLogx();

        Double_t margin                 = relativeMarginsX[0]*2.7*1350;
        Double_t textsizeLabelsWidth    = 0;
        Double_t textsizeFacWidth       = 0;
        if (padWidthEta->XtoPixel(padWidthEta->GetX2()) < padWidthEta->YtoPixel(padWidthEta->GetY1())){
            textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthEta->XtoPixel(padWidthEta->GetX2()) ;
            textsizeFacWidth            = (Double_t)1./padWidthEta->XtoPixel(padWidthEta->GetX2()) ;
        } else {
            textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthEta->YtoPixel(padWidthEta->GetY1());
            textsizeFacWidth            = (Double_t)1./padWidthEta->YtoPixel(padWidthEta->GetY1());
        }

        TH2F * histo2DAllEtaFWHM    = new TH2F("histo2DAllEtaFWHM","histo2DAllEtaFWHM", 20, 0.23, 15. ,1000., -4.99, 20.5);
        SetStyleHistoTH2ForGraphs(histo2DAllEtaFWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                                  0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin), 512, 505);
        // histo2DAllEtaFWHM->GetYaxis()->SetRangeUser(-1.,45.5);
        histo2DAllEtaFWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
        histo2DAllEtaFWHM->GetYaxis()->SetNdivisions(505);
        histo2DAllEtaFWHM->GetYaxis()->SetNoExponent(kTRUE);
        histo2DAllEtaFWHM->GetXaxis()->SetTickLength(0.05);
        histo2DAllEtaFWHM->GetYaxis()->SetTickLength(0.026);
        histo2DAllEtaFWHM->DrawCopy();


                DrawGammaSetMarker(histoEtaTrueFWHMMeV[1], markerstylesMC[1], markersizeMC[1], colorMC[1] , colorMC[1]);
                histoEtaTrueFWHMMeV[1]->Draw("p,same,e");
                DrawGammaSetMarker(histoEtaFWHMMeV[1], markerstyles[1], markersize[1], colorData[1] , colorData[1]);
                histoEtaFWHMMeV[1]->Draw("p,same,e");

                DrawGammaSetMarker(histoEtaTrueFWHMMeV[2], markerstylesMC[2], markersizeMC[2], colorMC[2] , colorMC[2]);
                histoEtaTrueFWHMMeV[2]->Draw("p,same,e");
                DrawGammaSetMarker(histoEtaFWHMMeV[2], markerstyles[2], markersize[2], colorData[2] , colorData[2]);
                histoEtaFWHMMeV[2]->Draw("p,same,e");

                DrawGammaSetMarker(histoEtaTrueFWHMMeV[0], markerstylesMC[0], markersizeMC[0], colorMC[0] , colorMC[0]);
                histoEtaTrueFWHMMeV[0]->Draw("p,same,e");
                DrawGammaSetMarker(histoEtaFWHMMeV[0], markerstyles[0], markersize[0], colorData[0] , colorData[0]);
                histoEtaFWHMMeV[0]->Draw("p,same,e");

        drawLatexAdd("b)",0.13,0.06,textSizeLabelsPixel,kFALSE,kTRUE);
        drawLatexAdd("ALICE",0.13,0.87,textSizeLabelsPixel,kFALSE,kTRUE);
        drawLatexAdd("pp, PCM",0.13,0.78,textSizeLabelsPixel,kFALSE,kTRUE);
        drawLatexAdd("#eta #rightarrow #gamma#gamma",0.13,0.69,textSizeLabelsPixel,kFALSE,kTRUE);

    padMassEta->cd();
    padMassEta->SetLogx();

        Double_t textsizeLabelsMassEta         = 0;
        Double_t textsizeFacMassEta            = 0;
        if (padMassEta->XtoPixel(padMassEta->GetX2()) <padMassEta->YtoPixel(padMassEta->GetY1()) ){
            textsizeLabelsMassEta              = (Double_t)textSizeLabelsPixel/padMassEta->XtoPixel(padMassEta->GetX2()) ;
            textsizeFacMassEta                 = (Double_t)1./padMassEta->XtoPixel(padMassEta->GetX2()) ;
        } else {
            textsizeLabelsMassEta              = (Double_t)textSizeLabelsPixel/padMassEta->YtoPixel(padMassEta->GetY1());
            textsizeFacMassEta                 = (Double_t)1./padMassEta->YtoPixel(padMassEta->GetY1());
        }

        TH2F * histo2DAllEtaMass            = new TH2F("histo2DAllEtaMass","histo2DAllEtaMass",20, 0.23, 15., 1000., 535.1, 575.9);
        SetStyleHistoTH2ForGraphs(histo2DAllEtaMass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMassEta, textsizeLabelsMassEta, 0.85*textsizeLabelsMassEta,
                                  textsizeLabelsMassEta, 0.9, 0.28/(textsizeFacMassEta*margin), 512, 505);
        histo2DAllEtaMass->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DAllEtaMass->GetYaxis()->SetNdivisions(505);
        histo2DAllEtaMass->GetYaxis()->SetNoExponent(kTRUE);
        histo2DAllEtaMass->GetXaxis()->SetTickLength(0.05);
        histo2DAllEtaMass->GetXaxis()->SetLabelOffset(-0.015);
        histo2DAllEtaMass->DrawCopy();

                DrawGammaSetMarker(histoEtaTrueMass[1], markerstylesMC[1], markersizeMC[1], colorMC[1] , colorMC[1]);
                histoEtaTrueMass[1]->Draw("p,same,e");
                DrawGammaSetMarker(histoEtaMass[1], markerstyles[1], markersize[1], colorData[1] , colorData[1]);
                histoEtaMass[1]->Draw("p,same,e");

                DrawGammaSetMarker(histoEtaTrueMass[2], markerstylesMC[2], markersizeMC[2], colorMC[2] , colorMC[2]);
                histoEtaTrueMass[2]->Draw("p,same,e");
                DrawGammaSetMarker(histoEtaMass[2], markerstyles[2], markersize[2], colorData[2] , colorData[2]);
                histoEtaMass[2]->Draw("p,same,e");

                DrawGammaSetMarker(histoEtaTrueMass[0], markerstylesMC[0], markersizeMC[0], colorMC[0] , colorMC[0]);
                histoEtaTrueMass[0]->Draw("p,same,e");
                DrawGammaSetMarker(histoEtaMass[0], markerstyles[0], markersize[0], colorData[0] , colorData[0]);
                histoEtaMass[0]->Draw("p,same,e");

        DrawGammaLines(0.23, 15. , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,0.1, kGray);

        drawLatexAdd("b)",0.13,0.22,textSizeLabelsPixel,kFALSE,kTRUE);


        //********************************** Defintion of the Legend **************************************************
        Double_t columnsLegendMass2Eta[3]      = {0.,0.57,0.84};
        Double_t rowsLegendMass2Eta[5]         = {0.8,0.6,0.4,0.2,0.01};
        //******************* Offsets ***********************
        Double_t offsetMarkerXMass2Eta         = 0.1;
        Double_t offsetMarkerYMass2Eta         = 0.1;
        //****************** Scale factors ******************
        Double_t scaleMarkerMass2Eta           = 1.2;

        padMassLegend1Eta->cd();
        //****************** first Column **************************************************
        TLatex *textMassPCMEta[10];
        Int_t counterEta = 1;
        for (Int_t i = 0; i < 3; i++){
            if(histoEtaMass[i] && histoEtaTrueMass[i] && histoEtaFWHMMeV[i] && histoEtaTrueFWHMMeV[i]&&i!=3){
                textMassPCMEta[i]                  = new TLatex(columnsLegendMass2Eta[0],rowsLegendMass2Eta[counterEta],nameMeasGlobal[i].Data());
                SetStyleTLatex( textMassPCMEta[i], textSizeLabelsPixel,4);
                textMassPCMEta[i]->SetTextFont(43);
                textMassPCMEta[i]->Draw();
                counterEta+=1;
            }
        }
        //****************** second Column *************************************************
        TLatex *textMassDataEta                = new TLatex(columnsLegendMass2Eta[1],rowsLegendMass2Eta[0] ,"Data");
        SetStyleTLatex( textMassDataEta, textSizeLabelsPixel,4);
        textMassDataEta->SetTextFont(43);
        textMassDataEta->Draw();
        TLatex *textMassMCEta                  = new TLatex(columnsLegendMass2Eta[2] ,rowsLegendMass2Eta[0],"MC");
        SetStyleTLatex( textMassMCEta, textSizeLabelsPixel,4);
        textMassMCEta->SetTextFont(43);
        textMassMCEta->Draw();

        TMarker* markerPCMEtaMass[10];
        TMarker* markerPCMEtaMassMC[10];
        counterEta = 1;
        for (Int_t i = 0; i < 3; i++){
            if(histoEtaMass[i] && histoEtaTrueMass[i]&&i!=3){
                markerPCMEtaMass[i]             = CreateMarkerFromHisto(histoEtaMass[i],columnsLegendMass2Eta[1]+ offsetMarkerXMass2Eta ,rowsLegendMass2Eta[counterEta]+ offsetMarkerYMass2Eta ,scaleMarkerMass2Eta);
                markerPCMEtaMass[i]->DrawMarker(columnsLegendMass2Eta[1]+ offsetMarkerXMass2Eta ,rowsLegendMass2Eta[counterEta]+ offsetMarkerYMass2Eta);
                markerPCMEtaMassMC[i]           = CreateMarkerFromHisto(histoEtaTrueMass[i],columnsLegendMass2Eta[2]+ offsetMarkerXMass2Eta ,rowsLegendMass2Eta[counterEta]+ offsetMarkerYMass2Eta ,scaleMarkerMass2Eta);
                markerPCMEtaMassMC[i]->DrawMarker(columnsLegendMass2Eta[2]+ offsetMarkerXMass2Eta-0.04 ,rowsLegendMass2Eta[counterEta]+ offsetMarkerYMass2Eta);
                counterEta+=1;
            }
        }

    canvasMassWidthEta->Update();
    canvasMassWidthEta->Print(Form("%s/Eta_MassAndWidth.%s",outputDir.Data(),suffix.Data()));
}

void CombineGammaResultsPP(

){
  gStyle->SetEndErrorSize(0);
    // define input file names
    TString dataSetFileNames[6];
    TString gammaResultsFileNames[6];
    TString dataGammaFileNames[6];
    TString secondaryFileNames[6];
    TString DCAFileNames2[6];
    TString dataCorrectionHistos[6];
    TString DCAFileNames[6];
    TString DCAFileNamesMC2[6];
    TString DCAFileNamesMC[6];
    TString dataSetFileNamesQA[nSets];
    TString MCFileNames[nSetsMC];
    TString chargedFileName[6];
    // 900 GeV
     gammaResultsFileNames[0]                        = "ThesisQAInputAllEnergies/900GeV_gamma_00000113_00200009227302008250404000_0152103500000000/900GeV/Gamma_Pi0_data_GammaConvV1_InclusiveRatio.root";
     secondaryFileNames[0]                        = "ThesisQAInputAllEnergies/900GeV_meson_00000113_00200009227302008250404000_0152103500000000/900GeV/Pi0_data_GammaConvV1Correction_00000113_00200009227302008250404000_0152103500000000.root";
     dataGammaFileNames[0]                        = "ThesisQAInputAllEnergies/900GeV_gamma_00000113_00200009227302008250404000_0152103500000000/900GeV/Gamma_Pi0_data_GammaConvV1Correction_00000113_00200009227302008250404000_0152103500000000.root";
     dataCorrectionHistos[0]                        = "ThesisQAInputAllEnergies/900GeV_gamma_00000113_00200009227302008250404000_0152103500000000/900GeV/Pi0_MC_GammaConvV1CorrectionHistos_00000113_00200009227302008250404000_0152103500000000.root";
     DCAFileNames[0]                        = "ThesisQAInputAllEnergies/900GeV_meson_00000113_00200009227302008250404000_0152103500000000/900GeV/Pi0_Data_GammaConvV1DCATestAnalysed.root";
     DCAFileNamesMC2[0]                        = "ThesisQAInputAllEnergies/900GeV_meson_00000113_00200009227302008250404000_0152103500000000/900GeV/Pi0_MC_GammaConvV1DCATestAnalysed.root";
     DCAFileNames2[0]                        = "ThesisQAInputAllEnergies/900GeV_gamma_00000113_00200009227302008250404000_0152103500000000/900GeV/Pi0_data_GammaConvV1DCAHistogramms_00000113_00200009227302008250404000_0152103500000000.root";
     DCAFileNamesMC[0]                        = "ThesisQAInputAllEnergies/900GeV_gamma_00000113_00200009227302008250404000_0152103500000000/900GeV/Pi0_MC_GammaConvV1DCAHistogramms_00000113_00200009227302008250404000_0152103500000000.root";
     dataSetFileNames[0]                        = "ThesisQAInputAllEnergies/900GeV_data_PCMResultsFullCorrection_PP_20170225.root";
     dataSetFileNamesQA[0]                        = "ThesisQAInputAllEnergies/00000113_00200009227300008250404000_0152103500000000/900GeV/EventQA/LHC10c_900GeV_EventQARunwise.root";
     MCFileNames[0]                             = "ThesisQAInputAllEnergies/00000113_00200009227300008250404000_0152103500000000/900GeV/EventQA/LHC14j4c_900GeV_EventQARunwise.root";
     chargedFileName[0]                             = "ThesisQAInputAllEnergies/ExtInput/ChargedIdentifiedSpectraPP_15_Feb_2016.root";
     //2.76 TeV
     gammaResultsFileNames[3]                        = "ThesisQAInputAllEnergies/2.76TeV/Gamma_Pi0_data_GammaConvV1_InclusiveRatio.root";
     dataSetFileNames[3]                        = "ThesisQAInputAllEnergies/2.76TeV/CombinedResultsPaperPP2760GeV_2017_03_01_FrediV2Clusterizer.root";
     dataGammaFileNames[3]                        = "CombinationInputPP/2.76TeV/gamma_00000113_00200009397300008250400000_0163103100900000/2.76TeV/Gamma_Pi0_data_GammaConvV1Correction_00000113_00200009397300008250400000_0163103100900000.root";
    dataCorrectionHistos[3]                        = "CombinationInputPP/2.76TeV/gamma_00000113_00200009397300008250400000_0163103100900000/2.76TeV/Pi0_MC_GammaConvV1CorrectionHistos_00000113_00200009397300008250400000_0163103100900000.root";
    secondaryFileNames[3]                        = "CombinationInputPP/2.76TeV/gamma_00000113_00200009397300008250400000_0163103100900000/2.76TeV/Pi0_data_GammaConvV1Correction_00000113_00200009397300008250400000_0163103100900000.root";

    // 7 TeV
     gammaResultsFileNames[1]                        = "ThesisQAInputAllEnergies/7TeV_gamma_00000113_00200009227302008250404000_0152103500000000/7TeV/Gamma_Pi0_data_GammaConvV1_InclusiveRatio.root"; //std binning
//      gammaResultsFileNames[1]                        = "ThesisQAInputAllEnergies/7TeV_gamma_00000113_00200009227302008250404000_0152103500000000/7TeV/Gamma_Pi0_data_GammaConvV1_InclusiveRatio_FINE.root";  // 0.1 binning
     secondaryFileNames[1]                        = "ThesisQAInputAllEnergies/7TeV_meson_00000113_00200009227302008250404000_0152103500000000/7TeV/Pi0_data_GammaConvV1Correction_00000113_00200009227300008250404000_0152103500000000.root";
     dataGammaFileNames[1]                        = "ThesisQAInputAllEnergies/7TeV_gamma_00000113_00200009227302008250404000_0152103500000000/7TeV/Gamma_Pi0_data_GammaConvV1Correction_00000113_00200009227300008250404000_0152103500000000.root";
     dataCorrectionHistos[1]                        = "ThesisQAInputAllEnergies/7TeV_gamma_00000113_00200009227302008250404000_0152103500000000/7TeV/Pi0_MC_GammaConvV1CorrectionHistos_00000113_00200009227300008250404000_0152103500000000.root";
     DCAFileNames[1]                        = "ThesisQAInputAllEnergies/7TeV_meson_00000113_00200009227302008250404000_0152103500000000/7TeV/Pi0_Data_GammaConvV1DCATestAnalysed.root";
     DCAFileNamesMC2[1]                        = "ThesisQAInputAllEnergies/7TeV_meson_00000113_00200009227302008250404000_0152103500000000/7TeV/Pi0_MC_GammaConvV1DCATestAnalysed.root";
     DCAFileNames2[1]                        = "ThesisQAInputAllEnergies/7TeV_gamma_00000113_00200009227302008250404000_0152103500000000/7TeV/Pi0_data_GammaConvV1DCAHistogramms_00000113_00200009227300008250404000_0152103500000000.root";
     DCAFileNamesMC[1]                        = "ThesisQAInputAllEnergies/7TeV_gamma_00000113_00200009227302008250404000_0152103500000000/7TeV/Pi0_MC_GammaConvV1DCAHistogramms_00000113_00200009227300008250404000_0152103500000000.root";
     dataSetFileNames[1]                        = "ThesisQAInputAllEnergies/data_PCMResultsFullCorrection_PP_7TeV_20170714.root";
//      dataSetFileNames[1]                        = "ThesisQAInputAllEnergies/7TeV_data_PCMResultsFullCorrection_PP_20170225.root";
//      dataSetFileNames[1]                        = "ThesisQAInputAllEnergies/data_PCMResultsFullCorrection_PP_7TeV_20170709.root";
     dataSetFileNamesQA[1]                        = "ThesisQAInputAllEnergies/00000113_00200009227302008250404000_0152103500000000/7TeV/EventQA/LHC10b_pass4_EventQARunwise.root";
     dataSetFileNamesQA[2]                        = "ThesisQAInputAllEnergies/00000113_00200009227302008250404000_0152103500000000/7TeV/EventQA/LHC10c_pass4_EventQARunwise.root";
     dataSetFileNamesQA[3]                        = "ThesisQAInputAllEnergies/00000113_00200009227302008250404000_0152103500000000/7TeV/EventQA/LHC10d_pass4_EventQARunwise.root";
     dataSetFileNamesQA[4]                        = "ThesisQAInputAllEnergies/00000113_00200009227302008250404000_0152103500000000/7TeV/EventQA/LHC10e_pass4_EventQARunwise.root";
     dataSetFileNamesQA[5]                        = "ThesisQAInputAllEnergies/00000113_00200009227302008250404000_0152103500000000/7TeV/EventQA/LHC10f_pass4_EventQARunwise.root";
     MCFileNames[1]                             = "ThesisQAInputAllEnergies/00000113_00200009227302008250404000_0152103500000000/7TeV/EventQA/LHC14j4_EventQARunwise.root";
     chargedFileName[1]                             = "ThesisQAInputAllEnergies/ExtInput/ChargedIdentifiedSpectraPP_15_Feb_2016.root";
    // 8 TeV
     gammaResultsFileNames[2]                        = "ThesisQAInputAllEnergies/8TeV_gamma_00010113_00200009227300008250404000_0152103500000000/8TeV/Gamma_Pi0_data_GammaConvV1_InclusiveRatio.root"; //std binning
//      gammaResultsFileNames[2]                        = "ThesisQAInputAllEnergies/8TeV_gamma_00010113_00200009227300008250404000_0152103500000000/8TeV/Gamma_Pi0_data_GammaConvV1_InclusiveRatio_FINE.root";  //0.1 binning
     secondaryFileNames[2]                        = "ThesisQAInputAllEnergies/8TeV_meson_00010113_00200009227300008250404000_0152103500000000/8TeV/Pi0_data_GammaConvV1Correction_00010113_00200009227300008250404000_0152103500000000.root";
     dataGammaFileNames[2]                        = "ThesisQAInputAllEnergies/8TeV_gamma_00010113_00200009227300008250404000_0152103500000000/8TeV/Gamma_Pi0_data_GammaConvV1Correction_00010113_00200009227300008250404000_0152103500000000.root";
     dataCorrectionHistos[2]                        = "ThesisQAInputAllEnergies/8TeV_gamma_00010113_00200009227300008250404000_0152103500000000/8TeV/Pi0_MC_GammaConvV1CorrectionHistos_00010113_00200009227300008250404000_0152103500000000.root";
//      DCAFileNames[2]                        = "ThesisQAInputAllEnergies/pp5TeV_hikari_AnalyseDCATests/Pi0_Data_GammaConvV1DCATestAnalysed.root";
//      DCAFileNamesMC2[2]                        = "ThesisQAInputAllEnergies/pp5TeV_hikari_AnalyseDCATests/Pi0_MC_GammaConvV1DCATestAnalysed.root";
     DCAFileNames[2]                        = "ThesisQAInputAllEnergies/8TeV_meson_00010113_00200009227300008250404000_0152103500000000/8TeV/Pi0_Data_GammaConvV1DCATestAnalysed.root";
     DCAFileNamesMC2[2]                        = "ThesisQAInputAllEnergies/8TeV_meson_00010113_00200009227300008250404000_0152103500000000/8TeV/Pi0_MC_GammaConvV1DCATestAnalysed.root";
     DCAFileNames2[2]                        = "ThesisQAInputAllEnergies/8TeV_gamma_00010113_00200009227300008250404000_0152103500000000/8TeV/Pi0_data_GammaConvV1DCAHistogramms_00010113_00200009227300008250404000_0152103500000000.root";
     DCAFileNamesMC[2]                        = "ThesisQAInputAllEnergies/8TeV_gamma_00010113_00200009227300008250404000_0152103500000000/8TeV/Pi0_MC_GammaConvV1DCAHistogramms_00010113_00200009227300008250404000_0152103500000000.root";
//      dataSetFileNames[2]                        = "ThesisQAInputAllEnergies/8TeV_data_PCMResultsFullCorrection_PP_2017_02_01.root";
//      dataSetFileNames[2]                        = "ThesisQAInputAllEnergies/8TeV_data_PCMResultsFullCorrection_PP_2017_04_25.root";
//      dataSetFileNames[2]                        = "ThesisQAInputAllEnergies/8TeV_data_PCMResultsFullCorrection_PP_2017_06_08.root";
//      dataSetFileNames[2]                        = "ThesisQAInputAllEnergies/data_PCMResultsFullCorrection_PP_TESTaf.root";
//      dataSetFileNames[2]                        = "ThesisQAInputAllEnergies/data_PCMResultsFullCorrection_PP_test0406.root";
//      dataSetFileNames[2]                        = "ThesisQAInputAllEnergies/data_PCMResultsFullCorrection_PP_8TeV_20170612.root";
//      dataSetFileNames[2]                        = "ThesisQAInputAllEnergies/data_PCMResultsFullCorrection_PP_8TeV_optimzedPileup.root";
//      dataSetFileNames[2]                        = "ThesisQAInputAllEnergies/data_PCMResultsFullCorrection_PP_8TeV_nsigmacut.root";
//      dataSetFileNames[2]                        = "ThesisQAInputAllEnergies/data_PCMResultsFullCorrection_PP_8TeV_20170619_pythia8.root";
//      dataSetFileNames[2]                        = "ThesisQAInputAllEnergies/data_PCMResultsFullCorrection_PP_8TeV_newTrainoutput.root";
//      dataSetFileNames[2]                        = "ThesisQAInputAllEnergies/data_PCMResultsFullCorrection_PP_8TeV_nocentralCathode.root";
//      dataSetFileNames[2]                        = "ThesisQAInputAllEnergies/data_PCMResultsFullCorrection_PP_8TeV_20170630.root";
//      dataSetFileNames[2]                        = "ThesisQAInputAllEnergies/data_PCMResultsFullCorrection_PP_8TeV_DCAiter-1.root";
//      dataSetFileNames[2]                        = "ThesisQAInputAllEnergies/data_PCMResultsFullCorrection_PP_8TeV_20170708.root";
//      dataSetFileNames[2]                        = "ThesisQAInputAllEnergies/data_PCMResultsFullCorrection_PP_8TeV_newSystematics_20170709.root";
//      dataSetFileNames[2]                        = "ThesisQAInputAllEnergies/data_PCMResultsFullCorrection_PP_8TeV_20170710_newEta.root";
//      dataSetFileNames[2]                        = "ThesisQAInputAllEnergies/data_PCMResultsFullCorrection_PP_8TeV_20170711_newEtaBinning.root";
//      dataSetFileNames[2]                        = "ThesisQAInputAllEnergies/data_PCMResultsFullCorrection_PP_8TeV_20170711_FINAL.root";
     dataSetFileNames[2]                        = "ThesisQAInputAllEnergies/data_PCMResultsFullCorrection_PP_8TeV_20170712_etaTo7GeV.root";
     dataSetFileNamesQA[6]                        = "ThesisQAInputAllEnergies/00010113_00200009227300008250404000_0152103500000000/8TeV/EventQA/LHC12a_EventQARunwise.root";
     dataSetFileNamesQA[7]                        = "ThesisQAInputAllEnergies/00010113_00200009227300008250404000_0152103500000000/8TeV/EventQA/LHC12b_EventQARunwise.root";
     dataSetFileNamesQA[8]                        = "ThesisQAInputAllEnergies/00010113_00200009227300008250404000_0152103500000000/8TeV/EventQA/LHC12c_EventQARunwise.root";
     dataSetFileNamesQA[9]                        = "ThesisQAInputAllEnergies/00010113_00200009227300008250404000_0152103500000000/8TeV/EventQA/LHC12d_EventQARunwise.root";
     dataSetFileNamesQA[10]                       = "ThesisQAInputAllEnergies/00010113_00200009227300008250404000_0152103500000000/8TeV/EventQA/LHC12f_EventQARunwise.root";
     dataSetFileNamesQA[11]                       = "ThesisQAInputAllEnergies/00010113_00200009227300008250404000_0152103500000000/8TeV/EventQA/LHC12h_EventQARunwise.root";
     dataSetFileNamesQA[12]                       = "ThesisQAInputAllEnergies/00010113_00200009227300008250404000_0152103500000000/8TeV/EventQA/LHC12i_EventQARunwise.root";
     MCFileNames[2]                             = "ThesisQAInputAllEnergies/00010113_00200009227300008250404000_0152103500000000/8TeV/EventQA/LHC15h1_EventQARunwise.root";
     MCFileNames[3]                             = "ThesisQAInputAllEnergies/00010113_00200009227300008250404000_0152103500000000/8TeV/EventQA/LHC15h2_EventQARunwise.root";
     chargedFileName[2]                             = "ThesisQAInputAllEnergies/ExtInput/8TeV_extrapolated_spectra.root";


     chargedFileName[3]                             = "ThesisQAInputAllEnergies/ExtInput/Pi_K_P_7TeV_INELSpectra_Paper2016_20150803.root";

    // load input files
    TFile* inputFileData[15];
    TFile* inputFileGammaResults[15];
    TFile* inputFileSecondary[15];
    TFile* inputFileGamma[15];
    TFile* inputFileCorrectionHistos[15];
    TFile* inputFileDCA[15];
    TFile* inputFileDCAHistos[15];
    TFile* inputFileDCAMC[15];
    TFile* inputFileDCAHistosMC[15];
    TFile* inputFileDataQA[15];
    TFile* inputFileMC[10];
    TFile* chargedFiles[10];
    for(Int_t i=0; i<nSets; i++){
        inputFileDataQA[i]                      = new TFile(dataSetFileNamesQA[i].Data());
        if(i<nSetsMC)
            inputFileMC[i]                      = new TFile(MCFileNames[i].Data());
        if(i<4)
            inputFileData[i]                    = new TFile(dataSetFileNames[i].Data());
        if(i<4)
            inputFileGamma[i]                   = new TFile(dataGammaFileNames[i].Data());
        if(i<4)
            inputFileCorrectionHistos[i]        = new TFile(dataCorrectionHistos[i].Data());
        // Errors if files are not found
        if(i<3 && !inputFileData[i])
            cout << i << " " << dataSetFileNames[i].Data() << " not found x!" << endl;
        if(!inputFileDataQA[i])
            cout << dataSetFileNamesQA[i].Data() << " not found y!" << endl;
        if(i<nSetsMC && !inputFileMC[i])
            cout << MCFileNames[i].Data() << " not found z!" << endl;
        if(i<4)
            inputFileSecondary[i]                      = new TFile(secondaryFileNames[i].Data());
        if(i<3)
            inputFileDCA[i]                      = new TFile(DCAFileNames[i].Data());
        if(i<3)
            inputFileDCAHistos[i]                      = new TFile(DCAFileNames2[i].Data());
        if(i<3)
            inputFileDCAMC[i]                      = new TFile(DCAFileNamesMC2[i].Data());
        if(i<3)
            inputFileDCAHistosMC[i]                      = new TFile(DCAFileNamesMC[i].Data());
        if(i<4)
            inputFileGammaResults[i]                      = new TFile(gammaResultsFileNames[i].Data());
        if(i<4)
            chargedFiles[i]                      = new TFile(chargedFileName[i].Data());
    }

    TString binshiftFileName                    = "ThesisQAInputAllEnergies/BinShift_data_PCMResultsFullCorrection_PP_NoBinShifting.root";
    TFile* inputFileBinShift                      = new TFile(binshiftFileName.Data());

    TString date                                = ReturnDateString();


    gROOT->Reset();
    gROOT->SetStyle("Plain");

    StyleSettingsThesis();
    SetPlotStyle();

    TString dateForOutput                       = ReturnDateStringForOutput();
    cout << dateForOutput.Data() << endl;
    //___________________________________ Declaration of files _____________________________________________
    TString collisionSystem900GeV               = "pp, #sqrt{#it{s}} = 900 GeV";
    TString collisionSystem7TeV                 = "pp, #sqrt{#it{s}} = 7 TeV";
    TString collisionSystem8TeV                 = "pp, #sqrt{#it{s}} = 8 TeV";

    outputDir                           = Form("%s/%s/CombineGammaResultsPP",suffix.Data(),dateForOutput.Data());
    gSystem->Exec("mkdir -p "+outputDir);


    Double_t mesonMassExpectPi0                 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    Double_t mesonMassExpectEta                 = TDatabasePDG::Instance()->GetParticle(221)->Mass();

    for (Int_t i = 0; i < nSets; i++){
        hNEventsFracGoodEvents[i]               = (TH1D*)inputFileDataQA[i]->Get(Form("hNEventsFracGoodEvents_%s%s",nameDataset[i].Data(),nameDatasetAdd[i].Data()));
        hPi0Mass[i]                             = (TH1D*)inputFileDataQA[i]->Get(Form("hPi0Mass_%s%s",nameDataset[i].Data(),nameDatasetAdd[i].Data()));
        hEtaMass[i]                             = (TH1D*)inputFileDataQA[i]->Get(Form("hEtaMass_%s%s",nameDataset[i].Data(),nameDatasetAdd[i].Data()));
        hTracksGoodMean[i]                      = (TH1D*)inputFileDataQA[i]->Get(Form("hTracksGood-Mean_%s%s",nameDataset[i].Data(),nameDatasetAdd[i].Data()));
        hConvNCandidatesQA[i]                   = (TH1D*)inputFileDataQA[i]->Get(Form("hConvNCandidatesQA_%s%s",nameDataset[i].Data(),nameDatasetAdd[i].Data()));
        hFracPileup[i]                          = (TH1D*)inputFileDataQA[i]->Get(Form("hFracPileUp_%s%s",nameDataset[i].Data(),nameDatasetAdd[i].Data()));
        hFracWVtxOutside10cm[i]                 = (TH1D*)inputFileDataQA[i]->Get(Form("hFracWVtxOutside10cm_%s%s",nameDataset[i].Data(),nameDatasetAdd[i].Data()));
        hPi0Frac[i]                             = (TH1D*)inputFileDataQA[i]->Get(Form("hPi0Frac_%s%s",nameDataset[i].Data(),nameDatasetAdd[i].Data()));
        hVertexZMean[i]                         = (TH1D*)inputFileDataQA[i]->Get(Form("hVertexZ-Mean_%s%s",nameDataset[i].Data(),nameDatasetAdd[i].Data()));
        if(i<nSetsMC){
            hNEventsFracGoodEventsMC[i]               = (TH1D*)inputFileMC[i]->Get(Form("hNEventsFracGoodEvents_%s",nameDatasetMCLoad[i].Data()));
            hPi0MassMC[i]                             = (TH1D*)inputFileMC[i]->Get(Form("hPi0Mass_%s",nameDatasetMCLoad[i].Data()));
            hEtaMassMC[i]                             = (TH1D*)inputFileMC[i]->Get(Form("hEtaMass_%s",nameDatasetMCLoad[i].Data()));
            hTracksGoodMeanMC[i]                      = (TH1D*)inputFileMC[i]->Get(Form("hTracksGood-Mean_%s",nameDatasetMCLoad[i].Data()));
            hConvNCandidatesQAMC[i]                   = (TH1D*)inputFileMC[i]->Get(Form("hConvNCandidatesQA_%s",nameDatasetMCLoad[i].Data()));
            hFracPileupMC[i]                          = (TH1D*)inputFileMC[i]->Get(Form("hFracPileUp_%s",nameDatasetMCLoad[i].Data()));
            hFracWVtxOutside10cmMC[i]                 = (TH1D*)inputFileMC[i]->Get(Form("hFracWVtxOutside10cm_%s",nameDatasetMCLoad[i].Data()));
            hPi0FracMC[i]                             = (TH1D*)inputFileMC[i]->Get(Form("hPi0Frac_%s",nameDatasetMCLoad[i].Data()));
            hVertexZMeanMC[i]                         = (TH1D*)inputFileMC[i]->Get(Form("hVertexZ-Mean_%s",nameDatasetMCLoad[i].Data()));
        }
    }


    plotHistograms(hNEventsFracGoodEvents, hNEventsFracGoodEventsMC,0.651,1.0,"N^{good Evt}/N^{MinBias}","hNEventsFracGoodEvents");
    plotHistograms(hPi0Mass, hPi0MassMC,0.13,0.14,"m_{#pi^{0}}","hPi0Mass");
    plotHistograms(hFracPileup, hFracPileupMC,-0.002,0.04,"N^{pileup}/N^{Evt}","hFracPileUp");
    plotHistograms(hConvNCandidatesQA, hConvNCandidatesQAMC,0.0,0.349,"N_{#gamma_{conv}, Candidates}/N_{Events}","hConvNCandidatesQA");
    plotHistograms(hTracksGoodMean, hTracksGoodMeanMC,2.5,8.5,"#bar{#lower[0.1]{N}}_{ Good Tracks}","hTracksGood");
    plotHistograms(hFracWVtxOutside10cm, hFracWVtxOutside10cmMC,-0.01,0.17,"N^{zVTX>10cm}/N^{Evt}","hFracWVtxOutside10cm");
    plotHistograms(hVertexZMean, hVertexZMeanMC,-0.03,2.49,"#bar{#lower[0.1]{z}}_{Vertex} (cm)","hVertexZMean");


    // ANALYSIS PLOTS (MASS; WIDTH; SPECTRA; ...)
    directoryPi0[0]                     = (TDirectory*)inputFileData[0]->Get("Pi0900GeV");
    directoryEta[0]                     = (TDirectory*)inputFileData[0]->Get("Eta900GeV");
    directoryPi0[1]                     = (TDirectory*)inputFileData[1]->Get("Pi07TeV");
    directoryEta[1]                     = (TDirectory*)inputFileData[1]->Get("Eta7TeV");
    directoryPi0[2]                     = (TDirectory*)inputFileData[2]->Get("Pi08TeV");
    directoryEta[2]                     = (TDirectory*)inputFileData[2]->Get("Eta8TeV");
    directoryPi0[3]                     = (TDirectory*)inputFileData[2]->Get("Pi02.76TeV");
    directoryEta[3]                     = (TDirectory*)inputFileData[2]->Get("Eta2.76TeV");

    directoryPi02[0]                     = (TDirectory*)inputFileBinShift->Get("Pi0900GeV");
    directoryPi02[1]                     = (TDirectory*)inputFileBinShift->Get("Pi07TeV");
    directoryPi02[2]                     = (TDirectory*)inputFileBinShift->Get("Pi08TeV");

    TString categoryName = "cat1";
    TString backgroundExtractionMethod = "std";
    for (Int_t i = 0; i < 3; i++){


        histoEventQualityMC[i]                               = (TH1F*)inputFileCorrectionHistos[i]->Get("NEvents");
        if(inputFileCorrectionHistos[i]&&histoEventQualityMC[i])
        nEvtMC[i]                                        = GetNEvents(histoEventQualityMC[i]);

        TruePrimaryConvGamma_Pt[i]                           = (TH1D*)inputFileCorrectionHistos[i]->Get("TruePrimaryConvGamma_Pt");

        for(Int_t j = 0;j<17;j++){
            if(inputFileCorrectionHistos[i])
            histoCombinatorialSpecies_Pt[i][j]               = (TH1D*)inputFileCorrectionHistos[i]->Get(Form("ESD_TrueComb%s_Pt",combinatorics[j].Data()));
            else
            histoCombinatorialSpecies_Pt[i][j]=NULL;
            if(histoCombinatorialSpecies_Pt[i][j]){
                histoCombinatorialSpecies_Pt[i][j]->SetMinimum(1e-10);
                histoCombinatorialSpecies_Pt[i][j]->Scale(1./nEvtMC[i]);


                histoSignalToCombBackgroundRatio[i][j] = (TH1D*) histoCombinatorialSpecies_Pt[i][j]->Clone(Form("ESD_TrueCombRatioSignal_%s_Pt",combinatorics[i].Data()));
                histoSignalToCombBackgroundRatio[i][j]->Scale(nEvtMC[i]);
                histoSignalToCombBackgroundRatio[i][j]->Divide(histoSignalToCombBackgroundRatio[i][j],TruePrimaryConvGamma_Pt[i],1,1,"");
                SetHistogramm(histoSignalToCombBackgroundRatio[i][j],"#it{p}_{T} (GeV/#it{c})","identified BG/ real primary photons",10,5e7);
                histoSignalToCombBackgroundRatio[i][j]->SetMinimum(1e-5);
            }
        }


        ESDGammaDCAzAllMC[i]                           = (TH1D*)inputFileDCAHistosMC[i]->Get(Form("ESD_TrueGammaPtDCAz_Full_%s_perEvent", categoryName.Data()));
        ESDGammaDCAzAll3MC[i]                          = (TH1D*)inputFileDCAHistosMC[i]->Get(Form("ESD_TrueGammaPtDCAz_Full_%s_perEvent", "cat3"));
        fDCAEventQualityMC[i]                          = (TH1D*)inputFileDCAMC[i]->Get("NEvents");
        if(fDCAEventQualityMC[i])
        fNEventsDCAMC[i]                               = GetNEvents(fDCAEventQualityMC[i]);
        else
        fNEventsDCAMC[i]                          = -1;

        ESDGammaDCAzAll[i]                             = (TH1D*)inputFileDCAHistos[i]->Get(Form("ESD_GammaPtDCAzBin_Full_%s_perEvent", categoryName.Data()));
        ESDGammaDCAzBack[i]                            = (TH1D*)inputFileDCAHistos[i]->Get(Form("ESD_GammaPtDCAzBackBin_Full_%s_%s_perEvent", categoryName.Data(), backgroundExtractionMethod.Data()));
        ESDGammaDCAzAllSub[i]                          = (TH1D*)inputFileDCAHistos[i]->Get(Form("ESD_SubGammaPtDCAzBin_Full_%s_%s_perEvent", categoryName.Data(), backgroundExtractionMethod.Data()));

        ESDGammaDCAzAll3[i]                             = (TH1D*)inputFileDCAHistos[i]->Get(Form("ESD_GammaPtDCAzBin_Full_%s_perEvent", "cat3"));
        ESDGammaDCAzBack3[i]                            = (TH1D*)inputFileDCAHistos[i]->Get(Form("ESD_GammaPtDCAzBackBin_Full_%s_%s_perEvent", "cat3", backgroundExtractionMethod.Data()));
        ESDGammaDCAzAllSub3[i]                          = (TH1D*)inputFileDCAHistos[i]->Get(Form("ESD_SubGammaPtDCAzBin_Full_%s_%s_perEvent", "cat3", backgroundExtractionMethod.Data()));
        fDCAEventQuality[i]                          = (TH1D*)inputFileDCA[i]->Get("NEvents");
        if(fDCAEventQuality[i])
        fNEventsDCA[i]                               = GetNEvents(fDCAEventQuality[i]);
        else
        fNEventsDCA[i]                          = -1;


        histoGammaSpecCorrPurity[i]                       = (TH1D*)inputFileGammaResults[i]->Get("histoGammaSpecCorrPurity");
        graphGammaSpecCorrPurity[i]                       = new TGraphAsymmErrors(histoGammaSpecCorrPurity[i]);
        if(i!=1&&i!=2&&i!=0)
        ESD_TruePrimaryConvGammaESD_PtMCPt_SystErr[i]     = (TGraphAsymmErrors*)inputFileGammaResults[i]->Get("ESD_TruePrimaryConvGammaESD_PtMCPt_SystErr");
        else
        ESD_TruePrimaryConvGammaESD_PtMCPt_SystErr[i]     = (TGraphAsymmErrors*)inputFileGammaResults[i]->Get("histoGammaSpecCorrPurity_SystErr");


        IncRatioPurity_trueEff[i]                       = (TH1D*)inputFileGammaResults[i]->Get("IncRatioPurity_trueEff");
        graphIncRatioPurity_trueEff[i]                       = new TGraphAsymmErrors(IncRatioPurity_trueEff[i]);
        IncRatioPurity_trueEff_SystErr[i]     = (TGraphAsymmErrors*)inputFileGammaResults[i]->Get("IncRatioPurity_trueEff_SystErr");
        if(i!=1&&i!=2&&i!=0)
        DoubleRatioConversionTrueEffPurity[i]                       = (TH1D*)inputFileGammaResults[i]->Get("DoubleRatioConversionFitPurity");
        // else
        // DoubleRatioConversionTrueEffPurity[i]                       = (TH1D*)inputFileGammaResults[i]->Get("DoubleRatioTrueEffPurity");
        else
        DoubleRatioConversionTrueEffPurity[i]                       = (TH1D*)inputFileGammaResults[i]->Get("DoubleRatioFitPurity");
        graphDoubleRatioConversionTrueEffPurity[i]                       = new TGraphAsymmErrors(DoubleRatioConversionTrueEffPurity[i]);

        if(i!=1&&i!=2&&i!=0)
        DoubleRatioConversionTrueEffPurity_SystErr[i]     = (TGraphAsymmErrors*)inputFileGammaResults[i]->Get("DoubleRatioConversionFitPurity_SystErr");
        // else
        // DoubleRatioConversionTrueEffPurity_SystErr[i]     = (TGraphAsymmErrors*)inputFileGammaResults[i]->Get("DoubleRatioTrueEffPurity_SystErr");
        else
        DoubleRatioConversionTrueEffPurity_SystErr[i]     = (TGraphAsymmErrors*)inputFileGammaResults[i]->Get("DoubleRatioFitPurity_SystErr");

        histoDirectPhotonSpectrum[i]                       = (TH1D*)inputFileGammaResults[i]->Get("histoDirectPhotonSpectrum");
        graphDirectPhotonSpectrum[i]                       = new TGraphAsymmErrors(histoDirectPhotonSpectrum[i]);
        graphNLODirGamma[i]     = (TGraphAsymmErrors*)inputFileGammaResults[i]->Get("graphNLODirGamma");

        //2.76TeV
        if(i==2){
            histoGammaSpecCorrPurity[i+1]                       = (TH1D*)inputFileGammaResults[i+1]->Get("histoGammaSpecCorrPurity");
            GammaRaw_Pt[i+1]                                      = (TH1D*)inputFileGamma[i+1]->Get("GammaRaw_Pt");
            graphGammaSpecCorrPurity[i+1]                       = new TGraphAsymmErrors(histoGammaSpecCorrPurity[i+1]);
            ESD_TruePrimaryConvGammaESD_PtMCPt_SystErr[i+1]     = (TGraphAsymmErrors*)inputFileGammaResults[i+1]->Get("ESD_TruePrimaryConvGammaESD_PtMCPt_SystErr");


            IncRatioPurity_trueEff[i+1]                       = (TH1D*)inputFileGammaResults[i+1]->Get("IncRatioPurity_trueEff");
            graphIncRatioPurity_trueEff[i+1]                       = new TGraphAsymmErrors(IncRatioPurity_trueEff[i+1]);
            IncRatioPurity_trueEff_SystErr[i+1]     = (TGraphAsymmErrors*)inputFileGammaResults[i+1]->Get("IncRatioPurity_trueEff_SystErr");

            DoubleRatioConversionTrueEffPurity[i+1]                       = (TH1D*)inputFileGammaResults[i+1]->Get("DoubleRatioConversionFitPurity");
            graphDoubleRatioConversionTrueEffPurity[i+1]                       = new TGraphAsymmErrors(DoubleRatioConversionTrueEffPurity[i+1]);
            DoubleRatioConversionTrueEffPurity_SystErr[i+1]     = (TGraphAsymmErrors*)inputFileGammaResults[i+1]->Get("DoubleRatioConversionFitPurity_SystErr");

            histoDirectPhotonSpectrum[i+1]                       = (TH1D*)inputFileGammaResults[i+1]->Get("histoDirectPhotonSpectrum");
            graphDirectPhotonSpectrum[i+1]                       = new TGraphAsymmErrors(histoDirectPhotonSpectrum[i+1]);
            graphNLODirGamma[i+1]     = (TGraphAsymmErrors*)inputFileGammaResults[i+1]->Get("graphNLODirGamma");
            histoGammaPurity_Pt[i+1]                            = (TH1D*)inputFileCorrectionHistos[i+1]->Get("GammaPurity_Pt");
            GammaRecoEff_MCPt[i+1]                              = (TH1D*)inputFileGamma[i+1]->Get("GammaRecoEff_MCPt");

        }

        histoGammaPurity_Pt[i]                            = (TH1D*)inputFileCorrectionHistos[i]->Get("GammaPurity_Pt");
        GammaRecoEff_MCPt[i]                              = (TH1D*)inputFileGamma[i]->Get("GammaRecoEff_MCPt");
        GammaConvProb_Pt[i]                               = (TH1D*)inputFileGamma[i]->Get("GammaConvProb_Pt");
        PileUpCorrectionFactor[i]                            = (TH1D*)inputFileGamma[i]->Get("PileUpCorrectionFactor");
        if(i==2)
        PileUpCorrectionFactor[3]                            = (TH1D*)inputFileGamma[3]->Get("PileUpCorrectionFactor");
        // GAMMA INPUT
        SecondaryGammaFromXFromK0sRecoEff_Pt[i]              = (TH1D*)inputFileGamma[i]->Get("SecondaryGammaFromXFromK0sRecoEff_Pt");
        SecondaryGammaFromXFromK0lRecoEff_Pt[i]              = (TH1D*)inputFileGamma[i]->Get("SecondaryGammaFromXFromK0lRecoEff_Pt");
        SecondaryGammaFromXFromLambdaRecoEff_Pt[i]           = (TH1D*)inputFileGamma[i]->Get("SecondaryGammaFromXFromLambdaRecoEff_Pt");


        SecondaryGammaFromXFromK0sConvProb_MCPt[i]           = (TH1D*)inputFileGamma[i]->Get("SecondaryGammaFromXFromK0sMCGammaConvProb_MCPt");
        SecondaryGammaFromXFromK0lConvProb_MCPt[i]           = (TH1D*)inputFileGamma[i]->Get("SecondaryGammaFromXFromK0lMCGammaConvProb_MCPt");
        SecondaryGammaFromXFromLambdaConvProb_MCPt[i]        = (TH1D*)inputFileGamma[i]->Get("SecondaryGammaFromXFromLambdaConvProb_MCPt");


        GammaRaw_Pt[i]                                       = (TH1D*)inputFileGamma[i]->Get("GammaRaw_Pt");
        GammaRaw_Pt[i]->Sumw2();

        if(i!=1&&i!=2){
            histoGammaTrueSecConvGammaFromXFromK0s_Cocktail_Raw_Pt[i]           = (TH1D*)inputFileGamma[i]->Get("histoGammaTrueSecConvGammaFromXFromK0s_Cocktail_Raw_Pt");
            histoGammaTrueSecConvGammaFromXFromK0l_Cocktail_Raw_Pt[i]           = (TH1D*)inputFileGamma[i]->Get("histoGammaTrueSecConvGammaFromXFromK0l_Cocktail_Raw_Pt");
            histoGammaTrueSecConvGammaFromXFromLambda_Cocktail_Raw_Pt[i]        = (TH1D*)inputFileGamma[i]->Get("histoGammaTrueSecConvGammaFromXFromLambda_Cocktail_Raw_Pt");
        } else{
            histoGammaTrueSecConvGammaFromXFromK0s_Cocktail_Raw_Pt[i]           = (TH1D*)inputFileGamma[i]->Get("histoGammaTrueSecGammaFromXFromK0s_Cocktail_Raw_Pt");
            histoGammaTrueSecConvGammaFromXFromK0l_Cocktail_Raw_Pt[i]           = (TH1D*)inputFileGamma[i]->Get("histoGammaTrueSecGammaFromXFromK0l_Cocktail_Raw_Pt");
            histoGammaTrueSecConvGammaFromXFromLambda_Cocktail_Raw_Pt[i]        = (TH1D*)inputFileGamma[i]->Get("histoGammaTrueSecGammaFromXFromLambda_Cocktail_Raw_Pt");
        }
        histoGammaTrueSecConvGammaFromXFromK0s_Cocktail_Raw_Pt[i]->Sumw2();
        histoGammaTrueSecConvGammaFromXFromK0l_Cocktail_Raw_Pt[i]->Sumw2();
        histoGammaTrueSecConvGammaFromXFromLambda_Cocktail_Raw_Pt[i]->Sumw2();
        histoGammaTrueSecCocktailGammaRest_Pt[i]        = (TH1D*)inputFileGamma[i]->Get("histoGammaTrueSecCocktailGammaRest_Pt");
        histoGammaTrueSecCocktailGammaRest_Pt[i]->Sumw2();

        // calculate fractions
        histoFracAllGammaToSecFromXFromK0s_Cocktail_Pt[i]                  = (TH1D*)GammaRaw_Pt[i]->Clone("FracAllGammaToSecFromXFromK0s");
        histoFracAllGammaToSecFromXFromK0s_Cocktail_Pt[i] ->Divide(histoGammaTrueSecConvGammaFromXFromK0s_Cocktail_Raw_Pt[i] ,histoFracAllGammaToSecFromXFromK0s_Cocktail_Pt[i] ,1,1,"B");


        histoFracAllGammaToSecFromXFromK0l_Cocktail_Pt[i]                   = (TH1D*)GammaRaw_Pt[i]->Clone("FracAllGammaToSecFromXFromK0l");
        histoFracAllGammaToSecFromXFromK0l_Cocktail_Pt[i]->Divide(histoGammaTrueSecConvGammaFromXFromK0l_Cocktail_Raw_Pt[i] ,histoFracAllGammaToSecFromXFromK0l_Cocktail_Pt[i] ,1,1,"B");


        histoFracAllGammaToSecFromXFromLambda_Cocktail_Pt[i]                = (TH1D*)GammaRaw_Pt[i]->Clone("FracAllGammaToSecFromXFromLambda");
        histoFracAllGammaToSecFromXFromLambda_Cocktail_Pt[i] ->Divide(histoGammaTrueSecConvGammaFromXFromLambda_Cocktail_Raw_Pt[i] ,histoFracAllGammaToSecFromXFromLambda_Cocktail_Pt[i] ,1,1,"B");


        histoFracAllGammaToSecRest_Cocktail_Pt[i]                           = (TH1D*)GammaRaw_Pt[i]->Clone("FracAllGammaToSecRest");
        histoFracAllGammaToSecRest_Cocktail_Pt[i] ->Divide(histoGammaTrueSecCocktailGammaRest_Pt[i] ,histoFracAllGammaToSecRest_Cocktail_Pt[i] ,1,1,"B");


        // LOAD INPUT
        histoNumberOfEvents[i]              = (TH1D*)inputFileData[i]->Get("histoNumberOfEvents7TeV");
        BGEstimateFromPileup[i]              = (TH1D*)inputFileSecondary[i]->Get("BGEstimateFromPileup");
        if(i==2)
        BGEstimateFromPileup[3]              = (TH1D*)inputFileSecondary[3]->Get("BGEstimateFromPileup");
        PileupContamination[i]              = (TH1D*)inputFileSecondary[i]->Get("PileupContamination");
        if(i==2)
        PileupContamination[3]              = (TH1D*)inputFileSecondary[3]->Get("PileupContamination");
        fHistFracCat_1_vsPt[i]              = (TH1D*)inputFileDCA[i]->Get("fHistFracCat_1_vsPt");
        fHistFracCat_2_vsPt[i]              = (TH1D*)inputFileDCA[i]->Get("fHistFracCat_2_vsPt");
        fHistFracCat_3_vsPt[i]              = (TH1D*)inputFileDCA[i]->Get("fHistFracCat_3_vsPt");
        fHistFracCat_4_vsPt[i]              = (TH1D*)inputFileDCA[i]->Get("fHistFracCat_4_vsPt");
        fHistFracCat_5_vsPt[i]              = (TH1D*)inputFileDCA[i]->Get("fHistFracCat_5_vsPt");
        fHistFracCat_6_vsPt[i]              = (TH1D*)inputFileDCA[i]->Get("fHistFracCat_6_vsPt");
        fHistFracIntHistBGvsPt_Cat_1_Variant_1[i]              = (TH1D*)inputFileDCA[i]->Get("fHistFracIntHistBGvsPt_Cat_1_Variant_1");
        fHistFracIntHistBGvsPt_Cat_2_Variant_1[i]              = (TH1D*)inputFileDCA[i]->Get("fHistFracIntHistBGvsPt_Cat_2_Variant_1");
        fHistFracIntHistBGvsPt_Cat_3_Variant_1[i]              = (TH1D*)inputFileDCA[i]->Get("fHistFracIntHistBGvsPt_Cat_3_Variant_1");
        fHistFracIntHistBGvsPt_Cat_4_Variant_1[i]              = (TH1D*)inputFileDCA[i]->Get("fHistFracIntHistBGvsPt_Cat_4_Variant_1");
        fHistFracIntHistBGvsPt_Cat_5_Variant_1[i]              = (TH1D*)inputFileDCA[i]->Get("fHistFracIntHistBGvsPt_Cat_5_Variant_1");
        fHistFracIntHistBGvsPt_Cat_6_Variant_1[i]              = (TH1D*)inputFileDCA[i]->Get("fHistFracIntHistBGvsPt_Cat_6_Variant_1");


        HistDCAZUnderMesonCat_1_MesonPt_0304[i]              = (TH1D*)inputFileDCA[i]->Get("HistDCAZUnderMesonCat_1_MesonPt_0.50-0.60");
        fHistDCAZUnderMesonBGEstimateCat_1_MesonPt_0304[i]   = (TH1D*)inputFileDCA[i]->Get("fHistDCAZUnderMesonBGEstimateCat_1_MesonPt_0.50-0.60");
        HistDCAZUnderMesonCat_1_MesonPt_0304MCTRUE[i]        = (TH1D*)inputFileDCAMC[i]->Get("HistDCAZTruePrimaryMesonGammaGammaCat_1_MesonPt_0.50-0.60");
        HistDCAZUnderMesonCat_1_MesonPt_0304MC[i]              = (TH1D*)inputFileDCAMC[i]->Get("HistDCAZUnderMesonCat_1_MesonPt_0.50-0.60");

        // load invariant mass peak positions and widths
        histoPi0Mass[i]                     = (TH1D*)directoryPi0[i]->Get(Form("Pi0_Mass_data_%s",strTrigName[i].Data()));
        histoPi0FWHMMeV[i]                  = (TH1D*)directoryPi0[i]->Get(Form("Pi0_Width_data_%s",strTrigName[i].Data()));
        histoPi0TrueMass[i]                 = (TH1D*)directoryPi0[i]->Get(Form("Pi0_Mass_MC_%s",strTrigName[i].Data()));
        histoPi0TrueFWHMMeV[i]              = (TH1D*)directoryPi0[i]->Get(Form("Pi0_Width_MC_%s",strTrigName[i].Data()));
        histoPi0Mass[i]                     ->Scale(1000);
        histoPi0FWHMMeV[i]                  ->Scale(1000);
        histoPi0TrueMass[i]                 ->Scale(1000);
        histoPi0TrueFWHMMeV[i]              ->Scale(1000);
        BinShiftRatio[i]                    = (TH1D*)directoryPi02[i]->Get("BinShiftRatio");
        histoPi0InvCrossSection[i]          = (TH1D*)directoryPi0[i]->Get("InvCrossSectionPi0");
        histoPi0RawYields[i]                = (TH1D*)directoryPi0[i]->Get(Form("RAWYieldPerEventsPi0_%s",strTrigName[i].Data()));


        // load acceptance and efficiency and calculate acc*eff*y*2pi
        histoPi0Acc[i]                      = (TH1D*)directoryPi0[i]->Get(Form("AcceptancePi0_%s",strTrigName[i].Data()));
        histoPi0TrueEffPt[i]                = (TH1D*)directoryPi0[i]->Get(Form("EfficiencyPi0_%s",strTrigName[i].Data()));
        histoPi0AccEff[i]                   = (TH1D*)histoPi0Acc[i]->Clone("AcceptanceTimesEffPi0");
        histoPi0AccEff[i] ->Multiply(histoPi0TrueEffPt[i]);
        //         histoPi0AccTimesEff[i]              = (TH1D*)histoPi0TrueEffPt[i]->Clone(Form("histoPi0AccTimesEff%s",nameMeasGlobal[i].Data()));
        //         histoPi0AccTimesEff[i]->Multiply(histoPi0Acc[i]);
        //         histoPi0AccTimesEff[i]->Scale(2*TMath::Pi()*rapidityMeas[i]);


        // load cross section systematics and datapoints
        graphPi0InvCrossSectionStat[i]       = (TGraphAsymmErrors*)directoryPi0[i]->Get("graphInvCrossSectionPi0");
        graphPi0InvCrossSectionSys[i]       = (TGraphAsymmErrors*)directoryPi0[i]->Get("InvCrossSectionPi0Sys");

        histCorrectedYieldPi0[i]       = (TGraphAsymmErrors*)directoryPi0[i]->Get("CorrectedYieldPi0");
        graphCorrectedYieldPi0[i]       = (TGraphAsymmErrors*)directoryPi0[i]->Get("graphCorrectedYieldPi0");
        Pi0SystError[i]       = (TGraphAsymmErrors*)directoryPi0[i]->Get("Pi0SystError");

        //             graphPi0InvCrossSectionStat[i]      = new TGraphAsymmErrors(histoPi0InvCrossSection[i]);
        //                 graphPi0InvCrossSectionStat[i]  ->RemovePoint(0);


        // load eta invariant mass peak positions and widths
        if(directoryEta[i]){
            histoEtaRawYields[i]                = (TH1D*)directoryEta[i]->Get(Form("RAWYieldPerEventsEta_%s",strTrigName[i].Data()));
            graphEtaInvCrossSectionStat[i]       = (TGraphAsymmErrors*)directoryEta[i]->Get("graphInvCrossSectionEta");
            graphEtaInvCrossSectionSys[i]       = (TGraphAsymmErrors*)directoryEta[i]->Get("InvCrossSectionEtaSys");
            histoEtaMass[i]                 = (TH1D*)directoryEta[i]->Get(Form("Eta_Mass_data_%s",strTrigName[i].Data()));
            histoEtaFWHMMeV[i]              = (TH1D*)directoryEta[i]->Get(Form("Eta_Width_data_%s",strTrigName[i].Data()));
            histoEtaTrueMass[i]             = (TH1D*)directoryEta[i]->Get(Form("Eta_Mass_MC_%s",strTrigName[i].Data()));
            histoEtaTrueFWHMMeV[i]          = (TH1D*)directoryEta[i]->Get(Form("Eta_Width_MC_%s",strTrigName[i].Data()));
            histoEtaMass[i]                 ->Scale(1000);
            histoEtaFWHMMeV[i]              ->Scale(1000);
            histoEtaTrueMass[i]             ->Scale(1000);
            histoEtaTrueFWHMMeV[i]          ->Scale(1000);
            histoEtaAcc[i]                      = (TH1D*)directoryEta[i]->Get(Form("AcceptanceEta_%s",strTrigName[i].Data()));
            histoEtaTrueEffPt[i]                = (TH1D*)directoryEta[i]->Get(Form("EfficiencyEta_%s",strTrigName[i].Data()));
            histoEtaAccEff[i]                   = (TH1D*)histoEtaAcc[i]->Clone("AcceptanceTimesEffEta");
            histoEtaAccEff[i] ->Multiply(histoEtaTrueEffPt[i]);
            graphEtaToPi0Stat[i]             = (TGraphAsymmErrors*)directoryEta[i]->Get(Form("EtaToPi0StatError_%s",strTrigName[i].Data()));
            graphEtaToPi0Sys[i]             = (TGraphAsymmErrors*)directoryEta[i]->Get(Form("EtaToPi0SystError_%s",strTrigName[i].Data()));
            graphEtaToPi0ShiftStat[i]             = (TGraphAsymmErrors*)directoryEta[i]->Get(Form("EtaToPi0YShiftedStatError_%s",strTrigName[i].Data()));
            graphEtaToPi0ShiftSys[i]             = (TGraphAsymmErrors*)directoryEta[i]->Get(Form("EtaToPi0YShiftedSystError_%s",strTrigName[i].Data()));
        }else{
            histoEtaMass[i]                 = NULL;
            histoEtaFWHMMeV[i]              = NULL;
            histoEtaTrueMass[i]             = NULL;
            histoEtaTrueFWHMMeV[i]          = NULL;
        }
        TrueMesonEffiPt[i]                          = (TH1D*)inputFileSecondary[i]->Get("TrueMesonEffiPt");
        TrueSecFromK0SEffiPt[i]                     = (TH1D*)inputFileSecondary[i]->Get("TrueSecFromK0SEffiPt");
        TrueSecFromK0LEffiPt[i]                     = (TH1D*)inputFileSecondary[i]->Get("TrueSecFromK0LEffiPt");
        TrueSecFromLambdaEffiPt[i]                  = (TH1D*)inputFileSecondary[i]->Get("TrueSecFromLambdaEffiPt");
        TrueSecFromRestEffiPt[i]                    = (TH1D*)inputFileSecondary[i]->Get("TrueSecFromRestEffiPt");

        fMCMesonAcceptPt[i]                          = (TH1D*)inputFileSecondary[i]->Get("fMCMesonAccepPt");
        fMCSecPi0FromK0SAccepPt[i]                     = (TH1D*)inputFileSecondary[i]->Get("fMCSecPi0FromK0SAccepPt");
        fMCSecPi0FromK0LAccepPt[i]                     = (TH1D*)inputFileSecondary[i]->Get("fMCSecPi0FromK0LAccepPt");
        fMCSecPi0FromLambdaAccepPt[i]                  = (TH1D*)inputFileSecondary[i]->Get("fMCSecPi0FromLambdaAccepPt");
        fMCSecPi0FromRestAccepPt[i]                    = (TH1D*)inputFileSecondary[i]->Get("fMCSecPi0FromRestAccepPt");


        histoYieldMesonPerEvent[i]                          = (TH1D*)inputFileSecondary[i]->Get("histoYieldMesonPerEvent");
        SecYieldFromK0SMesonFromCocktail[i]                     = (TH1D*)inputFileSecondary[i]->Get("SecYieldFromK0SMesonFromCocktail");
        SecYieldFromK0LMesonFromCocktail[i]                     = (TH1D*)inputFileSecondary[i]->Get("SecYieldFromK0LMesonFromCocktail");
        SecYieldFromLambdaMesonFromCocktail[i]                  = (TH1D*)inputFileSecondary[i]->Get("SecYieldFromLambdaMesonFromCocktail");
        SecYieldFromRestMeson[i]                    = (TH1D*)inputFileSecondary[i]->Get("SecYieldFromRestMeson");


        RatioToRaw[i]                     = (TH1D*)inputFileSecondary[i]->Get("RatioToRaw");
        RatioSecYieldFromK0SMesonFromCocktailToRaw[i]                     = (TH1D*)inputFileSecondary[i]->Get("RatioSecYieldFromK0SMesonFromCocktailToRaw");
        RatioSecYieldFromK0LMesonFromCocktailToRaw[i]                     = (TH1D*)inputFileSecondary[i]->Get("RatioSecYieldFromK0LMesonFromCocktailToRaw");
        RatioSecYieldFromLambdaMesonFromCocktailToRaw[i]                  = (TH1D*)inputFileSecondary[i]->Get("RatioSecYieldFromLambdaMesonFromCocktailToRaw");
        RatioSecYieldFromRestMesonToRaw[i]                    = (TH1D*)inputFileSecondary[i]->Get("RatioSecYieldFromRestMesonToRaw");
    }




    histoChargedPionStat[0]                       = (TH1D*)chargedFiles[0]->Get("histoChargedPionSpecLowPtStat900GeVALICE");
        graphChargedPionStat[0]                       = new TGraphAsymmErrors(histoChargedPionStat[0]);

    histoChargedPionSys[0]                       = (TH1D*)chargedFiles[0]->Get("histoChargedPionSpecLowPtSys900GeVALICE");
        graphChargedPionSys[0]                       = new TGraphAsymmErrors(histoChargedPionSys[0]);


    histoChargedPionStat[1]                       = (TH1D*)chargedFiles[3]->Get("hstat_pp7_pion_sum");
    histoChargedPionStat[1]                      ->Sumw2();
    histoChargedPionStat[1]                      ->Scale(0.5);
        graphChargedPionStat[1]                       = new TGraphAsymmErrors(histoChargedPionStat[1]);
    histoChargedPionSys[1]                       = (TH1D*)chargedFiles[3]->Get("hsys_pp7_pion_sum");
    histoChargedPionSys[1]                      ->Sumw2();
    histoChargedPionSys[1]                      ->Scale(0.5);
        graphChargedPionSys[1]                       = new TGraphAsymmErrors(histoChargedPionSys[1]);
//     histoChargedPionStat[1]                       = (TH1D*)chargedFiles[1]->Get("histoChargedPionSpecLowPtStat7TeVALICE");
//         graphChargedPionStat[1]                       = new TGraphAsymmErrors(histoChargedPionStat[1]);
//     histoChargedPionSys[1]                       = (TH1D*)chargedFiles[1]->Get("histoChargedPionSpecLowPtSys7TeVALICE");
//         graphChargedPionSys[1]                       = new TGraphAsymmErrors(histoChargedPionSys[1]);


    histoChargedPionStat[2]                       = (TH1D*)chargedFiles[1]->Get("histoChargedPionSpecHighPtStat7TeVALICE");
        graphChargedPionStat[2]                       = new TGraphAsymmErrors(histoChargedPionStat[2]);
    graphChargedPionSys[2]                       = (TGraphAsymmErrors*)chargedFiles[1]->Get("graphChargedPionSpecHighPtSys7TeVALICE");

//     histoChargedPionStat[3]                       = (TH1D*)chargedFiles[2]->Get("histPion8_8TeV");
//         graphChargedPionStat[3]                       = new TGraphAsymmErrors(histoChargedPionStat[3]);
        graphChargedPionStat[3]                 = (TGraphAsymmErrors*)chargedFiles[2]->Get("graphStatErr_Pion_8TeV");
        graphChargedPionStat[3]->Print();
    graphChargedPionSys[3]                      = (TGraphAsymmErrors*)chargedFiles[2]->Get("graphSysErr_Pion_8TeV");
          graphChargedPionSys[3]->Print();

//     TH1D* histoChargedPionStat[10];
//     TH1D* histoChargedPionSys[10];
//     TGraphAsymmErrors* graphChargedPionStat[10];
//     TGraphAsymmErrors* graphChargedPionSys[10];

     TFile* cocktailFile[4];
   cocktailFile[0] = new TFile("ThesisQAInputAllEnergies/900GeV_gamma_00000113_00200009227302008250404000_0152103500000000/900GeV/GammaCocktail_0.80_00000113_00200009227302008250404000_0152103500000000.root");
   cocktailFile[1] = new TFile("ThesisQAInputAllEnergies/7TeV_gamma_00000113_00200009227302008250404000_0152103500000000/7TeV/GammaCocktail_0.80_00000113_00200009227300008250404000_0152103500000000.root");
   cocktailFile[2] = new TFile("ThesisQAInputAllEnergies/8TeV_gamma_00010113_00200009227300008250404000_0152103500000000/8TeV/GammaCocktail_0.80_00010113_00200009227300008250404000_0152103500000000.root");
   cocktailFile[3] = new TFile("CombinationInputPP/2.76TeV/gamma_00000113_00200009397300008250400000_0163103100900000/2.76TeV/GammaCocktail_0.80_00000113_00200009397300008250400000_0163103100900000.root");
    TH1D* cocktailPi0[3];
    TGraphAsymmErrors* graphcocktailPi0[3];
    TH1D* cocktailPi0Ratio[4];
    TGraphAsymmErrors* graphcocktailPi0Ratio[4];
    TGraphAsymmErrors* graphcocktailPi0RatioSys[3];
    TH1D* cocktailEta[3];
    TH1D* cocktailAllGamma[3];
    TH1D* cocktailAllGammaPi0[3];
//     for(Int_t i=0; i<3; i++){
//         cocktailPi0[i]                         = (TH1D* )cocktailFile[i]->Get("Pi0_Pt");
//         graphcocktailPi0[i]                    = new TGraphAsymmErrors(cocktailPi0[i]);
//         cocktailEta[i]                         = (TH1D* )cocktailFile[i]->Get("Eta_Pt");
//         cocktailAllGamma[i]                    = (TH1D* )cocktailFile[i]->Get("Gamma_Pt");
//         cocktailAllGammaPi0[i]                 = (TH1D* )cocktailFile[i]->Get("Gamma_From_Pi0_Pt");
//         cocktailPi0Ratio[i] = (TH1D* ) histCorrectedYieldPi0[i]->Clone("cocktailratio");
//         graphcocktailPi0RatioSys[i] = (TGraphAsymmErrors* ) Pi0SystError[i]->Clone("cocktailratiosys");
//         cocktailPi0Ratio[i]->Divide(cocktailPi0[i]);
//         graphcocktailPi0Ratio[i]= new TGraphAsymmErrors(cocktailPi0Ratio[i]);
//         graphcocktailPi0RatioSys[i] = CalculateAsymGraphRatioToGraph(graphcocktailPi0RatioSys[i],graphcocktailPi0[i]);
//     }
    for(Int_t i=0; i<4; i++){
//         cocktailPi0[i]                         = (TH1D* )cocktailFile[i]->Get("Pi0_Pt");
//         graphcocktailPi0[i]                    = new TGraphAsymmErrors(cocktailPi0[i]);
//         cocktailEta[i]                         = (TH1D* )cocktailFile[i]->Get("Eta_Pt");
//         cocktailAllGamma[i]                    = (TH1D* )cocktailFile[i]->Get("Gamma_Pt");
//         cocktailAllGammaPi0[i]                 = (TH1D* )cocktailFile[i]->Get("Gamma_From_Pi0_Pt");
//         cocktailPi0Ratio[i] = (TH1D* ) histCorrectedYieldPi0[i]->Clone("cocktailratio");
//         graphcocktailPi0RatioSys[i] = (TGraphAsymmErrors* ) Pi0SystError[i]->Clone("cocktailratiosys");
        cocktailPi0Ratio[i]         = (TH1D* )cocktailFile[i]->Get("ratioPi0DataCocktail");
        if(cocktailPi0Ratio[i] )
            cout << "found " << i << endl;
        graphcocktailPi0Ratio[i]= new TGraphAsymmErrors(cocktailPi0Ratio[i]);
//         graphcocktailPi0RatioSys[i] = CalculateAsymGraphRatioToGraph(graphcocktailPi0RatioSys[i],graphcocktailPi0[i]);
    }




//     plotPi0FWHMMass();
    plotPi0FWHMMassNEW();
//     plotEtaFWHMMass();
    plotEtaFWHMMassNEW();

    // RAW YIELD PLOTTING
    plotYield(histoPi0RawYields,1e-8, 2, "Pi0","Pi0RawYield");
    plotYield(GammaRaw_Pt,9e-7, 2000, "Gamma","GammaRawYield");
    plotYield(histoEtaRawYields,6e-8, 0.02, "Eta","EtaRawYield");

    plotCrossSection(graphPi0InvCrossSectionStat,graphPi0InvCrossSectionSys,5e3,4e14, "Pi0","Pi0InvCrossSection");
    plotCrossSection(graphEtaInvCrossSectionStat,graphEtaInvCrossSectionSys,5e4,4e12, "Eta","EtaInvCrossSection");

    plotCrossSectionFULL(graphPi0InvCrossSectionStat,graphPi0InvCrossSectionSys,5e4,4e14, "Pi0","Pi0InvCrossSectionFINAL");
    plotCrossSectionFULLEta(graphEtaInvCrossSectionStat,graphEtaInvCrossSectionSys,5e4,4e14, "Eta","EtaInvCrossSectionFINAL");



    plotCrossSection(graphGammaSpecCorrPurity,ESD_TruePrimaryConvGammaESD_PtMCPt_SystErr,5e-7,9000, "Gamma","InclusivePhotons");

    plotEtaToPi0(graphIncRatioPurity_trueEff,IncRatioPurity_trueEff_SystErr,0.0,1.09, "Ratio Inclusive #gamma/#pi^{0}","InclusiveRatio");

    plotEtaToPi0(graphDoubleRatioConversionTrueEffPurity,DoubleRatioConversionTrueEffPurity_SystErr,0.45,1.55, "(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})","DoubleRatio");
    plotEtaToPi0(graphcocktailPi0Ratio,graphcocktailPi0Ratio,0.41,1.75, "#pi^{0}_{data}/#pi^{0}_{cocktail}","CocktailRatio");
    plotDoubleRatio(graphDoubleRatioConversionTrueEffPurity,DoubleRatioConversionTrueEffPurity_SystErr,0.45,1.55, "(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})","DoubleRatio");




//     plotGammaPileup(2, ESDGammaDCAzAll,ESDGammaDCAzBack, ESDGammaDCAzAllSub,ESDGammaDCAzAllMC,4e-9,2e-5, "counts per event", 1, "GammaPileupBin");
    plotGammaPileup(2, ESDGammaDCAzAll,ESDGammaDCAzBack, ESDGammaDCAzAllSub,ESDGammaDCAzAllMC,3e-4,7, "counts per event", 1, "GammaPileupBin");
    plotGammaPileup(1, ESDGammaDCAzAll,ESDGammaDCAzBack, ESDGammaDCAzAllSub,ESDGammaDCAzAllMC,3e-4,7, "counts per event", 1, "GammaPileupBin");
    plotGammaDCAMultiple(1, ESDGammaDCAzAll,ESDGammaDCAzBack, ESDGammaDCAzAllSub,ESDGammaDCAzAllMC,3e-4,7, "normalized counts", 1, "GammaPileupBinPF");
//     plotGammaPileup(1, ESDGammaDCAzAll3,ESDGammaDCAzBack3, ESDGammaDCAzAllSub3,1.1,7e8, "counts", 1, "GammaPileupBinCat3");

    plotGammaSecondaries(0, SecondaryGammaFromXFromK0sRecoEff_Pt,SecondaryGammaFromXFromK0lRecoEff_Pt, SecondaryGammaFromXFromLambdaRecoEff_Pt,SecondaryGammaFromXFromRestRecoEff_Pt,0.0,1.0, "#epsilon_{eff} in |#eta| < 0.9", 0, "SecGammaEff");
    plotGammaSecondaries(0, SecondaryGammaFromXFromK0sConvProb_MCPt,SecondaryGammaFromXFromK0lConvProb_MCPt, SecondaryGammaFromXFromLambdaConvProb_MCPt,SecondaryGammaFromXFromRestConvProb_MCPt,0.0,8.5e-2, "#it{P}_{conv} in |#eta| < 0.9", 0, "SecGammaConvProb");
    plotGammaSecondaries(1, SecondaryGammaFromXFromK0sRecoEff_Pt,SecondaryGammaFromXFromK0lRecoEff_Pt, SecondaryGammaFromXFromLambdaRecoEff_Pt,SecondaryGammaFromXFromRestRecoEff_Pt,0.0,1.0, "#epsilon_{eff} in |#eta| < 0.9", 0, "SecGammaEff");
    plotGammaSecondaries(1, SecondaryGammaFromXFromK0sConvProb_MCPt,SecondaryGammaFromXFromK0lConvProb_MCPt, SecondaryGammaFromXFromLambdaConvProb_MCPt,SecondaryGammaFromXFromRestConvProb_MCPt,0.0,8.5e-2, "#it{P}_{conv} in |#eta| < 0.9", 0, "SecGammaConvProb");
    plotGammaSecondaries(2, SecondaryGammaFromXFromK0sRecoEff_Pt,SecondaryGammaFromXFromK0lRecoEff_Pt, SecondaryGammaFromXFromLambdaRecoEff_Pt,SecondaryGammaFromXFromRestRecoEff_Pt,0.0,1.0, "#epsilon_{eff} in |#eta| < 0.9", 0, "SecGammaEff");
    plotGammaSecondaries(2, SecondaryGammaFromXFromK0sConvProb_MCPt,SecondaryGammaFromXFromK0lConvProb_MCPt, SecondaryGammaFromXFromLambdaConvProb_MCPt,SecondaryGammaFromXFromRestConvProb_MCPt,0.0,8.5e-2, "#it{P}_{conv} in |#eta| < 0.9", 0, "SecGammaConvProb");


    plotGammaSecondaries(0, histoFracAllGammaToSecFromXFromK0s_Cocktail_Pt,histoFracAllGammaToSecFromXFromK0l_Cocktail_Pt, histoFracAllGammaToSecFromXFromLambda_Cocktail_Pt,histoFracAllGammaToSecRest_Cocktail_Pt,8e-7,9.9, "Effective Secondary #gamma Correction", 1, "SecGammaFractions");
    plotGammaSecondaries(1, histoFracAllGammaToSecFromXFromK0s_Cocktail_Pt,histoFracAllGammaToSecFromXFromK0l_Cocktail_Pt, histoFracAllGammaToSecFromXFromLambda_Cocktail_Pt,histoFracAllGammaToSecRest_Cocktail_Pt,8e-7,8., "Effective Secondary #gamma Correction", 1, "SecGammaFractions");
    plotGammaSecondaries(2, histoFracAllGammaToSecFromXFromK0s_Cocktail_Pt,histoFracAllGammaToSecFromXFromK0l_Cocktail_Pt, histoFracAllGammaToSecFromXFromLambda_Cocktail_Pt,histoFracAllGammaToSecRest_Cocktail_Pt,8e-7,8., "Effective Secondary #gamma Correction", 1, "SecGammaFractions");

    plotGammaSecondaries(0, histoGammaTrueSecConvGammaFromXFromK0s_Cocktail_Raw_Pt,histoGammaTrueSecConvGammaFromXFromK0l_Cocktail_Raw_Pt, histoGammaTrueSecConvGammaFromXFromLambda_Cocktail_Raw_Pt,histoGammaTrueSecCocktailGammaRest_Pt,5e-13,9e-2, "Secondary Converted #gamma Raw Yield", 1, "SecGammaYield");
    plotGammaSecondaries(1, histoGammaTrueSecConvGammaFromXFromK0s_Cocktail_Raw_Pt,histoGammaTrueSecConvGammaFromXFromK0l_Cocktail_Raw_Pt, histoGammaTrueSecConvGammaFromXFromLambda_Cocktail_Raw_Pt,histoGammaTrueSecCocktailGammaRest_Pt,5e-13,9e-2, "Secondary Converted #gamma Raw Yield", 1, "SecGammaYield");
    plotGammaSecondaries(2, histoGammaTrueSecConvGammaFromXFromK0s_Cocktail_Raw_Pt,histoGammaTrueSecConvGammaFromXFromK0l_Cocktail_Raw_Pt, histoGammaTrueSecConvGammaFromXFromLambda_Cocktail_Raw_Pt,histoGammaTrueSecCocktailGammaRest_Pt,5e-13,9e-2, "Secondary Converted #gamma Raw Yield", 1, "SecGammaYield");


    plotPhotonCombBG(1, histoSignalToCombBackgroundRatio,1e-5,1., "identified BG/real primary #gamma", 1, "GammaCombBGRatio");

    plotDCA(1, fHistFracCat_1_vsPt,fHistFracCat_2_vsPt, fHistFracCat_3_vsPt, fHistFracCat_4_vsPt, fHistFracCat_5_vsPt,fHistFracCat_6_vsPt,-0.3,100., "N_{meson per cat}/(N_{meson}) (%)", 0, "FracDCAcat");
    plotDCA(1, fHistFracIntHistBGvsPt_Cat_1_Variant_1,fHistFracIntHistBGvsPt_Cat_2_Variant_1, fHistFracIntHistBGvsPt_Cat_3_Variant_1, fHistFracIntHistBGvsPt_Cat_4_Variant_1, fHistFracIntHistBGvsPt_Cat_5_Variant_1,fHistFracIntHistBGvsPt_Cat_6_Variant_1,-0.3,14.9, "BG/Total (%)", 0, "FracDCABGcat7");
    plotDCA(2, fHistFracIntHistBGvsPt_Cat_1_Variant_1,fHistFracIntHistBGvsPt_Cat_2_Variant_1, fHistFracIntHistBGvsPt_Cat_3_Variant_1, fHistFracIntHistBGvsPt_Cat_4_Variant_1, fHistFracIntHistBGvsPt_Cat_5_Variant_1,fHistFracIntHistBGvsPt_Cat_6_Variant_1,-0.3,27.9, "BG/Total (%)", 0, "FracDCABGcat8");



    plotEverything(BGEstimateFromPileup, 0.541, 1.09,0, "Pileup Correction Factor", "PileupCorrectionFactorPi0");
    plotPeriodwisePileup(BGEstimateFromPileup, 0.451, 1.049,0, "Pileup Correction Factor", "PileupCorrectionFactorPi0Periods");
    plotEverything(PileupContamination, -0.2, 39,0, "Pileup Contamination [%]", "PileupContamination");
    plotEverything(BinShiftRatio, 0.95, 1.34,0, "value / binshifted value", "BinShift");

    plotEverything(histoGammaPurity_Pt, 0.89, 1.03,0, "#epsilon_{pur,#gamma} in |#eta| < 0.9", "PhotonPurity");
    plotEverything(GammaRecoEff_MCPt, 0.21, 1.01,0, "#epsilon_{eff,#gamma} in |#eta| < 0.9", "PhotonEfficiency");
    plotEverything(GammaConvProb_Pt, 0.055, .099,0, "P_{conv,#gamma} in |#eta| < 0.9", "PhotonConvProb");
    plotEverything(PileUpCorrectionFactor, 0.78, 1.04,0, "Pileup Correction Factor", "GammaPileUpCorrectionFactor");
    // plotEverything(Pi0spectraRatios, 0.78, 1.04,0, "8 to 7 TeV Ratio", "Pi0Ratio7to8TeV");
    // plotEverything(GammaspectraRatios, 0.78, 1.04,0, "8 to 7 TeV Ratio", "GammaRatio7to8TeV");


//     plotLuminosity(2e-4,100,"Peak Luminosity [Hz/#mub]","Luminosity");
}
