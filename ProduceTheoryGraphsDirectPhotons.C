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
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"

extern TRandom*    gRandom;
extern TBenchmark*    gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*      gMinuit;

// Double_t    xSection900GeV  = 47.78*1e-3;
// Double_t    xSection2760GeV = 55.416*1e-3;
Double_t    xSection7000GeV = 62.22*1e-3;
Double_t    xSection8000GeV = 55.74*1e-3;
// Double_t    recalcBarn      = 1e12; //NLO in pbarn!!!!


TGraph* ScaleGraph (TGraph* graph, Double_t scaleFac){
    TGraph* dummyGraph  = (TGraph*)graph->Clone(Form("%s_Scaled",graph->GetName()));
    Double_t * xValue   = dummyGraph->GetX();
    Double_t * yValue   = dummyGraph->GetY();
    Int_t nPoints       = dummyGraph->GetN();
    
    for (Int_t i = 0; i < nPoints; i++){
        yValue[i]       = yValue[i]*scaleFac;
    }
    TGraph* returnGraph = new TGraph(nPoints,xValue,yValue);
    return returnGraph;
}

TGraphAsymmErrors* ScaleGraphAsym (TGraphAsymmErrors* graph, Double_t scaleFac){
    TGraphAsymmErrors* dummyGraph   = (TGraphAsymmErrors*)graph->Clone(Form("%s_Scaled",graph->GetName()));

    Double_t * xValue               = dummyGraph->GetX(); 
    Double_t * yValue               = dummyGraph->GetY();
    Double_t* xErrorLow             = dummyGraph->GetEXlow();
    Double_t* xErrorHigh            = dummyGraph->GetEXhigh();
    Double_t* yErrorLow             = dummyGraph->GetEYlow();
    Double_t* yErrorHigh            = dummyGraph->GetEYhigh();
    Int_t nPoints = dummyGraph->GetN();
    for (Int_t i = 0; i < nPoints; i++){
        yValue[i] = yValue[i]*scaleFac;
        yErrorLow[i] = yErrorLow[i]*scaleFac;
        yErrorHigh[i] = yErrorHigh[i]*scaleFac;
    }
    TGraphAsymmErrors* returnGraph  = new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh); 
    return returnGraph;
}

//**********************************************************************************************************
// Calculates the ratio of two asymmerror graphs 
//**********************************************************************************************************
TGraphAsymmErrors* CalculateAsymGraphRatioToGraph(TGraphAsymmErrors* graphA, TGraphAsymmErrors* graphB){
  
    TGraphErrors* graphACopy        = (TGraphErrors*)graphA->Clone("GraphCopy");
    Double_t* xValue                = graphACopy->GetX(); 
    Double_t* yValue                = graphACopy->GetY();
    Double_t* xErrorHigh            = graphACopy->GetEXhigh();
    Double_t* xErrorLow             = graphACopy->GetEXlow();
    Double_t* yErrorHigh            = graphACopy->GetEYhigh();
    Double_t* yErrorLow             = graphACopy->GetEYlow();
   
    Int_t nPoints                   = graphACopy->GetN();
    for (Int_t i = 0; i < nPoints; i++){
        if (abs(xValue[i]-graphB->GetX()[i]) < 0.0001){
            if (graphB->GetY()[i] != 0){
                yValue[i]                   = yValue[i]/graphB->GetY()[i];
                Double_t yErrorRatioHigh    = yValue[i]*TMath::Sqrt( TMath::Power(graphA->GetEYhigh()[i]/graphA->GetY()[i],2) + TMath::Power(graphB->GetEYhigh()[i]/graphB->GetY()[i],2));
                yErrorHigh[i]               = yErrorRatioHigh; 
                Double_t yErrorRatioLow     = yValue[i]*TMath::Sqrt( TMath::Power(graphA->GetEYlow()[i]/graphA->GetY()[i],2) + TMath::Power(graphB->GetEYlow()[i]/graphB->GetY()[i],2));
                yErrorLow[i]                = yErrorRatioLow; 
            } else {
                yValue[i]                   = 0;
                yErrorHigh[i]               = 0;
                yErrorLow[i]                = 0;
            }    
        } else {
            cout << "ERROR: graphs don't have same binning " << endl;
            return NULL;
        }    
    }
    TGraphAsymmErrors* returnGraph      = new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh); 
    return returnGraph;
}

//**********************************************************************************************************
//****************************** Main function for theory compilation **************************************
//**********************************************************************************************************
void ProduceTheoryGraphsDirectPhotons(  Bool_t runPP    = kTRUE,
                                        Bool_t runPPb   = kTRUE,
                                        Bool_t runPbPb  = kTRUE,
                                        TString suffix  = "eps"
                                     ){    
    
    StyleSettingsThesis();    
    SetPlotStyle();

    Color_t colorFrag           = kAzure+6;
    Color_t colorPrompt         = kRed+2;
    Color_t colorDir            = kGreen+2;
    
    TString outputDir           = "DirectPhotonPlots/";
    gSystem->Exec("mkdir -p "+outputDir);
    
    TString collisionSystem2760GeV      = "pp, #sqrt{#it{s}} = 2.76 TeV";
    TString collisionSystem5TeV         = "pp, #sqrt{#it{s}} = 5.023 TeV";
    TString collisionSystem7TeV         = "pp, #sqrt{#it{s}} = 7 TeV";
    TString collisionSystem8TeV         = "pp, #sqrt{#it{s}} = 8 TeV";
    
    //******************************************************************************************************************
    //******************************* Put all direct photon calculations for pp in one file ****************************
    //******************************************************************************************************************
    if (runPP){
        // **********************************************************************************************************************
        // *********************************** direct photon calculations for 2.76TeV *******************************************
        // **********************************************************************************************************************

        TString fileNameNLOPhotonHalf2760GeV    = "ExternalInput/Theory/ALICENLOcalcDirectPhoton2760GeVmuhalf.dat";
        TString fileNameNLOPhotonOne2760GeV     = "ExternalInput/Theory/ALICENLOcalcDirectPhoton2760GeVmu.dat";
        TString fileNameNLOPhotonTwo2760GeV     = "ExternalInput/Theory/ALICENLOcalcDirectPhoton2760GeVtwomu.dat";
        Int_t nlinesNLOTwo2760GeV              = 0;
        Int_t nlinesNLOOne2760GeV              = 0;
        Int_t nlinesNLOHalf2760GeV             = 0;
        Double_t ptNLOPhotonTwo2760GeV[100];
        Double_t ptNLOPhotonTwo2760GeVErr[100];
        Double_t ptNLOPhotonOne2760GeV[100];
        Double_t ptNLOPhotonHalf2760GeV[100];
        
        Double_t muHalfDF2760GeV[100];
        Double_t muHalfD2760GeV[100];
        Double_t muHalfF2760GeV[100];
        Double_t muOneDF2760GeV[100];
        Double_t muOneD2760GeV[100];
        Double_t muOneF2760GeV[100];
        Double_t muTwoDF2760GeV[100];
        Double_t muTwoD2760GeV[100];
        Double_t muTwoF2760GeV[100];
        Double_t gammaValue2760GeV[100];
        Double_t gammaErrUp2760GeV[100];
        Double_t gammaErrDown2760GeV[100];
        Double_t fragGammaValue2760GeV[100];
        Double_t fragGammaErrUp2760GeV[100];
        Double_t fragGammaErrDown2760GeV[100];
        Double_t promptGammaValue2760GeV[100];
        Double_t promptGammaErrUp2760GeV[100];
        Double_t promptGammaErrDown2760GeV[100];
        
        ifstream inHalf2760GeV;
        inHalf2760GeV.open(fileNameNLOPhotonHalf2760GeV,ios_base::in);
        cout << fileNameNLOPhotonHalf2760GeV << endl;
            
        while(!inHalf2760GeV.eof()){
            nlinesNLOHalf2760GeV++;
            //               Pt                              DirectPhoton           FragmentationPhoton         SumPhoton
            inHalf2760GeV >> ptNLOPhotonHalf2760GeV[nlinesNLOHalf2760GeV] >> muHalfD2760GeV[nlinesNLOHalf2760GeV] >> muHalfF2760GeV[nlinesNLOHalf2760GeV] >> muHalfDF2760GeV[nlinesNLOHalf2760GeV];
            
        }
        inHalf2760GeV.close();

        ifstream inOne2760GeV;
        inOne2760GeV.open(fileNameNLOPhotonOne2760GeV,ios_base::in);
        cout << fileNameNLOPhotonOne2760GeV << endl;

        while(!inOne2760GeV.eof()){
            nlinesNLOOne2760GeV++;
            inOne2760GeV >> ptNLOPhotonOne2760GeV[nlinesNLOOne2760GeV] >> muOneD2760GeV[nlinesNLOOne2760GeV] >> muOneF2760GeV[nlinesNLOOne2760GeV] >> muOneDF2760GeV[nlinesNLOOne2760GeV];
        }
        inOne2760GeV.close();

        ifstream inTwo2760GeV;
        inTwo2760GeV.open(fileNameNLOPhotonTwo2760GeV,ios_base::in);
        cout << fileNameNLOPhotonTwo2760GeV << endl;

        while(!inTwo2760GeV.eof()){
            nlinesNLOTwo2760GeV++;
            inTwo2760GeV >> ptNLOPhotonTwo2760GeV[nlinesNLOTwo2760GeV] >> muTwoD2760GeV[nlinesNLOTwo2760GeV] >> muTwoF2760GeV[nlinesNLOTwo2760GeV] >> muTwoDF2760GeV[nlinesNLOTwo2760GeV];
            ptNLOPhotonTwo2760GeVErr[nlinesNLOTwo2760GeV] = 0.25;
        }
        inTwo2760GeV.close();
        
        for (Int_t i = 0; i < nlinesNLOTwo2760GeV; i++){
            cout << ptNLOPhotonHalf2760GeV[i] << "\t" << muHalfDF2760GeV[i] << "\t" << muOneDF2760GeV[i] << "\t" << muTwoDF2760GeV[i] << endl;
            gammaValue2760GeV[i]           = muOneDF2760GeV[i];
            gammaErrUp2760GeV[i]           = muHalfDF2760GeV[i]-muOneDF2760GeV[i];
            gammaErrDown2760GeV[i]         = muOneDF2760GeV[i]-muTwoDF2760GeV[i];
            fragGammaValue2760GeV[i]       = muOneF2760GeV[i];
            fragGammaErrUp2760GeV[i]       = muHalfF2760GeV[i]-muOneF2760GeV[i];
            fragGammaErrDown2760GeV[i]     = muOneF2760GeV[i]-muTwoF2760GeV[i];
            promptGammaValue2760GeV[i]     = muOneD2760GeV[i];
            promptGammaErrUp2760GeV[i]     = muHalfD2760GeV[i]-muOneD2760GeV[i];
            promptGammaErrDown2760GeV[i]   = muOneD2760GeV[i]-muTwoD2760GeV[i];
        }    
        
        TGraphAsymmErrors *graphNLOCalcDirGam2760GeV       = new TGraphAsymmErrors(nlinesNLOTwo2760GeV,ptNLOPhotonTwo2760GeV,gammaValue2760GeV,ptNLOPhotonTwo2760GeVErr,ptNLOPhotonTwo2760GeVErr,gammaErrDown2760GeV,gammaErrUp2760GeV);
        graphNLOCalcDirGam2760GeV->RemovePoint(0);
        TGraphAsymmErrors *graphNLOCalcFragGam2760GeV      = new TGraphAsymmErrors(nlinesNLOTwo2760GeV,ptNLOPhotonTwo2760GeV,fragGammaValue2760GeV,ptNLOPhotonTwo2760GeVErr,ptNLOPhotonTwo2760GeVErr,fragGammaErrDown2760GeV,fragGammaErrUp2760GeV);
        graphNLOCalcFragGam2760GeV->RemovePoint(0);
        TGraphAsymmErrors *graphNLOCalcPromGam2760GeV      = new TGraphAsymmErrors(nlinesNLOTwo2760GeV,ptNLOPhotonTwo2760GeV,promptGammaValue2760GeV,ptNLOPhotonTwo2760GeVErr,ptNLOPhotonTwo2760GeVErr,promptGammaErrDown2760GeV,promptGammaErrUp2760GeV);
        graphNLOCalcPromGam2760GeV->RemovePoint(0);
        
        TGraphAsymmErrors* graphRatioNLOFragGammaDivTot2760GeV      = CalculateAsymGraphRatioToGraph(graphNLOCalcFragGam2760GeV, graphNLOCalcDirGam2760GeV);
        TGraphAsymmErrors* graphRatioNLOPromptGammaDivTot2760GeV    = CalculateAsymGraphRatioToGraph(graphNLOCalcPromGam2760GeV, graphNLOCalcDirGam2760GeV);
        
        
        // -----------------------------------------------------------------------------------------------------------------------
        // ----------------------------- plotting fragmentation and prompt to total direct ---------------------------------------
        // -----------------------------------------------------------------------------------------------------------------------
        TCanvas* canvasRatioDirGammaCalc   = new TCanvas("canvasRatioDirGammaCalc","",200,10,900,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRatioDirGammaCalc, 0.11, 0.01, 0.01, 0.09);
        canvasRatioDirGammaCalc->SetLogx();
    
        TH2F * histo2DRatioGammaCalc;
        histo2DRatioGammaCalc           = new TH2F("histo2DRatioGammaCalc","histo2DRatioGammaCalc",11000,0.23,100.,1000,0.,1.25);
        SetStyleHistoTH2ForGraphs(histo2DRatioGammaCalc, "#it{p}_{T} (GeV/#it{c})","#frac{#gamma_{source}}{#gamma_{dir}}",0.035,0.04, 0.035,0.04, 1.,1.,510,505);
        histo2DRatioGammaCalc->GetXaxis()->SetMoreLogLabels();
        histo2DRatioGammaCalc->GetXaxis()->SetLabelOffset(-0.01);
        histo2DRatioGammaCalc->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphRatioNLOFragGammaDivTot2760GeV, 0, 0, colorFrag, colorFrag, 1, kTRUE, colorFrag);
            graphRatioNLOFragGammaDivTot2760GeV->Draw("3,same");

            DrawGammaSetMarkerTGraphAsym(graphRatioNLOPromptGammaDivTot2760GeV, 0, 0, colorPrompt, colorPrompt, 1, kTRUE, colorPrompt);
            graphRatioNLOPromptGammaDivTot2760GeV->SetFillStyle(3245);
            graphRatioNLOPromptGammaDivTot2760GeV->Draw("3,same");
            
            DrawGammaLines(0.23, 100. , 1., 1.,0.1, kGray+2);
            
            TLegend* legendRatioGammaCalc       = GetAndSetLegend2(0.15, 0.14, 0.3, 0.14+(0.035*2*1.25), 32);
            legendRatioGammaCalc->AddEntry(graphRatioNLOFragGammaDivTot2760GeV,"#gamma_{frag}/#gamma_{dir}");
            legendRatioGammaCalc->AddEntry(graphRatioNLOPromptGammaDivTot2760GeV,"#gamma_{prompt}/#gamma_{dir}");
            legendRatioGammaCalc->Draw();

            TLatex *labelRatioGammaCalc2760GeV   = new TLatex(0.15,0.93,collisionSystem2760GeV.Data());
            SetStyleTLatex( labelRatioGammaCalc2760GeV, 0.85*32,4);
            labelRatioGammaCalc2760GeV->SetTextFont(43);
            labelRatioGammaCalc2760GeV->Draw();
            TLatex *labelRatioGamma      = new TLatex(0.15,0.89,"#gamma_{dir}");
            SetStyleTLatex( labelRatioGamma, 0.85*32,4);
            labelRatioGamma->SetTextFont(43);
            labelRatioGamma->Draw();
            
        canvasRatioDirGammaCalc->SaveAs(Form("%s/GammaNLOCalc_Separation_PP2760GeV.%s",outputDir.Data(),suffix.Data()));

        // **********************************************************************************************************************
        // *********************************** direct photon calculations for 5TeV *******************************************
        // **********************************************************************************************************************
        
        TString fileNameNLOPhotonHalf5TeV    = "ExternalInput/Theory/ALICENLOcalcDirectPhoton5023GeVHalfMu.dat";
        TString fileNameNLOPhotonOne5TeV     = "ExternalInput/Theory/ALICENLOcalcDirectPhoton5023GeVMu.dat";
        TString fileNameNLOPhotonTwo5TeV     = "ExternalInput/Theory/ALICENLOcalcDirectPhoton5023GeVTwoMu.dat";
        Int_t nlinesNLOTwo5TeV              = 0;
        Int_t nlinesNLOOne5TeV              = 0;
        Int_t nlinesNLOHalf5TeV             = 0;
        Double_t ptNLOPhotonTwo5TeV[100];
        Double_t ptNLOPhotonTwo5TeVErr[100];
        Double_t ptNLOPhotonOne5TeV[100];
        Double_t ptNLOPhotonHalf5TeV[100];
        
        Double_t muHalfDF5TeV[100];
        Double_t muHalfD5TeV[100];
        Double_t muHalfF5TeV[100];
        Double_t muOneDF5TeV[100];
        Double_t muOneD5TeV[100];
        Double_t muOneF5TeV[100];
        Double_t muTwoDF5TeV[100];
        Double_t muTwoD5TeV[100];
        Double_t muTwoF5TeV[100];
        Double_t gammaValue5TeV[100];
        Double_t gammaErrUp5TeV[100];
        Double_t gammaErrDown5TeV[100];
        Double_t fragGammaValue5TeV[100];
        Double_t fragGammaErrUp5TeV[100];
        Double_t fragGammaErrDown5TeV[100];
        Double_t promptGammaValue5TeV[100];
        Double_t promptGammaErrUp5TeV[100];
        Double_t promptGammaErrDown5TeV[100];
        
        ifstream inHalf5TeV;
        inHalf5TeV.open(fileNameNLOPhotonHalf5TeV,ios_base::in);
        cout << fileNameNLOPhotonHalf5TeV << endl;
            
        while(!inHalf5TeV.eof()){
            nlinesNLOHalf5TeV++;
            //               Pt                              DirectPhoton           FragmentationPhoton         SumPhoton
            inHalf5TeV >> ptNLOPhotonHalf5TeV[nlinesNLOHalf5TeV] >> muHalfD5TeV[nlinesNLOHalf5TeV] >> muHalfF5TeV[nlinesNLOHalf5TeV] >> muHalfDF5TeV[nlinesNLOHalf5TeV];
            
        }
        inHalf5TeV.close();

        ifstream inOne5TeV;
        inOne5TeV.open(fileNameNLOPhotonOne5TeV,ios_base::in);
        cout << fileNameNLOPhotonOne5TeV << endl;

        while(!inOne5TeV.eof()){
            nlinesNLOOne5TeV++;
            inOne5TeV >> ptNLOPhotonOne5TeV[nlinesNLOOne5TeV] >> muOneD5TeV[nlinesNLOOne5TeV] >> muOneF5TeV[nlinesNLOOne5TeV] >> muOneDF5TeV[nlinesNLOOne5TeV];
        }
        inOne5TeV.close();

        ifstream inTwo5TeV;
        inTwo5TeV.open(fileNameNLOPhotonTwo5TeV,ios_base::in);
        cout << fileNameNLOPhotonTwo5TeV << endl;

        while(!inTwo5TeV.eof()){
            nlinesNLOTwo5TeV++;
            inTwo5TeV >> ptNLOPhotonTwo5TeV[nlinesNLOTwo5TeV] >> muTwoD5TeV[nlinesNLOTwo5TeV] >> muTwoF5TeV[nlinesNLOTwo5TeV] >> muTwoDF5TeV[nlinesNLOTwo5TeV];
            ptNLOPhotonTwo5TeVErr[nlinesNLOTwo5TeV] = 0.25;
        }
        inTwo5TeV.close();
        
        for (Int_t i = 0; i < nlinesNLOTwo5TeV; i++){
            cout << ptNLOPhotonHalf5TeV[i] << "\t" << muHalfDF5TeV[i] << "\t"<< ptNLOPhotonOne5TeV[i] <<"\t" << muOneDF5TeV[i] << "\t"<< ptNLOPhotonTwo5TeV[i] << "\t" << muTwoDF5TeV[i] << endl;
            gammaValue5TeV[i]           = muOneDF5TeV[i];
            gammaErrUp5TeV[i]           = muHalfDF5TeV[i]-muOneDF5TeV[i];
            gammaErrDown5TeV[i]         = muOneDF5TeV[i]-muTwoDF5TeV[i];
            fragGammaValue5TeV[i]       = muOneF5TeV[i];
            fragGammaErrUp5TeV[i]       = muHalfF5TeV[i]-muOneF5TeV[i];
            fragGammaErrDown5TeV[i]     = muOneF5TeV[i]-muTwoF5TeV[i];
            promptGammaValue5TeV[i]     = muOneD5TeV[i];
            promptGammaErrUp5TeV[i]     = muHalfD5TeV[i]-muOneD5TeV[i];
            promptGammaErrDown5TeV[i]   = muOneD5TeV[i]-muTwoD5TeV[i];
        }    
        
        TGraphAsymmErrors *graphNLOCalcDirGam5TeV       = new TGraphAsymmErrors(nlinesNLOTwo5TeV,ptNLOPhotonTwo5TeV,gammaValue5TeV,ptNLOPhotonTwo5TeVErr,ptNLOPhotonTwo5TeVErr,gammaErrDown5TeV,gammaErrUp5TeV);
        graphNLOCalcDirGam5TeV->RemovePoint(0);
        TGraphAsymmErrors *graphNLOCalcFragGam5TeV      = new TGraphAsymmErrors(nlinesNLOTwo5TeV,ptNLOPhotonTwo5TeV,fragGammaValue5TeV,ptNLOPhotonTwo5TeVErr,ptNLOPhotonTwo5TeVErr,fragGammaErrDown5TeV,fragGammaErrUp5TeV);
        graphNLOCalcFragGam5TeV->RemovePoint(0);
        TGraphAsymmErrors *graphNLOCalcPromGam5TeV      = new TGraphAsymmErrors(nlinesNLOTwo5TeV,ptNLOPhotonTwo5TeV,promptGammaValue5TeV,ptNLOPhotonTwo5TeVErr,ptNLOPhotonTwo5TeVErr,promptGammaErrDown5TeV,promptGammaErrUp5TeV);
        graphNLOCalcPromGam5TeV->RemovePoint(0);
        
        TGraphAsymmErrors* graphRatioNLOFragGammaDivTot5TeV      = CalculateAsymGraphRatioToGraph(graphNLOCalcFragGam5TeV, graphNLOCalcDirGam5TeV);
        TGraphAsymmErrors* graphRatioNLOPromptGammaDivTot5TeV    = CalculateAsymGraphRatioToGraph(graphNLOCalcPromGam5TeV, graphNLOCalcDirGam5TeV);
        
        
        // -----------------------------------------------------------------------------------------------------------------------
        // ----------------------------- plotting fragmentation and prompt to total direct ---------------------------------------
        // -----------------------------------------------------------------------------------------------------------------------
        
        canvasRatioDirGammaCalc->cd();    
        histo2DRatioGammaCalc->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphRatioNLOFragGammaDivTot5TeV, 0, 0, colorFrag, colorFrag, 1, kTRUE, colorFrag);
            graphRatioNLOFragGammaDivTot5TeV->Draw("3,same");

            DrawGammaSetMarkerTGraphAsym(graphRatioNLOPromptGammaDivTot5TeV, 0, 0, colorPrompt, colorPrompt, 1, kTRUE, colorPrompt);
            graphRatioNLOPromptGammaDivTot5TeV->SetFillStyle(3245);
            graphRatioNLOPromptGammaDivTot5TeV->Draw("3,same");
            
            DrawGammaLines(0.23, 100. , 1., 1.,0.1, kGray+2);
            
            legendRatioGammaCalc->Draw();

            TLatex *labelRatioGammaCalc5TeV   = new TLatex(0.15,0.93,collisionSystem5TeV.Data());
            SetStyleTLatex( labelRatioGammaCalc5TeV, 0.85*32,4);
            labelRatioGammaCalc5TeV->SetTextFont(43);
            labelRatioGammaCalc5TeV->Draw();
            labelRatioGamma->Draw();
            
        canvasRatioDirGammaCalc->SaveAs(Form("%s/GammaNLOCalc_Separation_PP5TeV.%s",outputDir.Data(),suffix.Data()));
        
        
        // **********************************************************************************************************************
        // *********************************** direct photon calculations for 7TeV *******************************************
        // **********************************************************************************************************************
        
        TString fileNameNLOPhotonHalf7TeV    = "ExternalInput/Theory/ALICENLOcalcDirectPhoton7TeVmuhalf.dat";
        TString fileNameNLOPhotonOne7TeV     = "ExternalInput/Theory/ALICENLOcalcDirectPhoton7TeVmu.dat";
        TString fileNameNLOPhotonTwo7TeV     = "ExternalInput/Theory/ALICENLOcalcDirectPhoton7TeVtwomu.dat";
        Int_t nlinesNLOTwo7TeV              = 0;
        Int_t nlinesNLOOne7TeV              = 0;
        Int_t nlinesNLOHalf7TeV             = 0;
        Double_t ptNLOPhotonTwo7TeV[100];
        Double_t ptNLOPhotonTwo7TeVErr[100];
        Double_t ptNLOPhotonOne7TeV[100];
        Double_t ptNLOPhotonHalf7TeV[100];
        
        Double_t muHalfDF7TeV[100];
        Double_t muHalfD7TeV[100];
        Double_t muHalfF7TeV[100];
        Double_t muOneDF7TeV[100];
        Double_t muOneD7TeV[100];
        Double_t muOneF7TeV[100];
        Double_t muTwoDF7TeV[100];
        Double_t muTwoD7TeV[100];
        Double_t muTwoF7TeV[100];
        Double_t gammaValue7TeV[100];
        Double_t gammaErrUp7TeV[100];
        Double_t gammaErrDown7TeV[100];
        Double_t fragGammaValue7TeV[100];
        Double_t fragGammaErrUp7TeV[100];
        Double_t fragGammaErrDown7TeV[100];
        Double_t promptGammaValue7TeV[100];
        Double_t promptGammaErrUp7TeV[100];
        Double_t promptGammaErrDown7TeV[100];
        
        ifstream inHalf7TeV;
        inHalf7TeV.open(fileNameNLOPhotonHalf7TeV,ios_base::in);
        cout << fileNameNLOPhotonHalf7TeV << endl;
            
        while(!inHalf7TeV.eof()){
            nlinesNLOHalf7TeV++;
            //               Pt                              DirectPhoton           FragmentationPhoton         SumPhoton
            inHalf7TeV >> ptNLOPhotonHalf7TeV[nlinesNLOHalf7TeV] >> muHalfD7TeV[nlinesNLOHalf7TeV] >> muHalfF7TeV[nlinesNLOHalf7TeV] >> muHalfDF7TeV[nlinesNLOHalf7TeV];
            
        }
        inHalf7TeV.close();

        ifstream inOne7TeV;
        inOne7TeV.open(fileNameNLOPhotonOne7TeV,ios_base::in);
        cout << fileNameNLOPhotonOne7TeV << endl;

        while(!inOne7TeV.eof()){
            nlinesNLOOne7TeV++;
            inOne7TeV >> ptNLOPhotonOne7TeV[nlinesNLOOne7TeV] >> muOneD7TeV[nlinesNLOOne7TeV] >> muOneF7TeV[nlinesNLOOne7TeV] >> muOneDF7TeV[nlinesNLOOne7TeV];
        }
        inOne7TeV.close();

        ifstream inTwo7TeV;
        inTwo7TeV.open(fileNameNLOPhotonTwo7TeV,ios_base::in);
        cout << fileNameNLOPhotonTwo7TeV << endl;

        while(!inTwo7TeV.eof()){
            nlinesNLOTwo7TeV++;
            inTwo7TeV >> ptNLOPhotonTwo7TeV[nlinesNLOTwo7TeV] >> muTwoD7TeV[nlinesNLOTwo7TeV] >> muTwoF7TeV[nlinesNLOTwo7TeV] >> muTwoDF7TeV[nlinesNLOTwo7TeV];
            ptNLOPhotonTwo7TeVErr[nlinesNLOTwo7TeV] = 0.25;
        }
        inTwo7TeV.close();
        
        for (Int_t i = 0; i < nlinesNLOTwo7TeV; i++){
            cout << ptNLOPhotonHalf7TeV[i] << "\t" << muHalfDF7TeV[i] << "\t"<< ptNLOPhotonOne7TeV[i] <<"\t" << muOneDF7TeV[i] << "\t"<< ptNLOPhotonTwo7TeV[i] << "\t" << muTwoDF7TeV[i] << endl;
            gammaValue7TeV[i]           = muOneDF7TeV[i];
            gammaErrUp7TeV[i]           = muHalfDF7TeV[i]-muOneDF7TeV[i];
            gammaErrDown7TeV[i]         = muOneDF7TeV[i]-muTwoDF7TeV[i];
            fragGammaValue7TeV[i]       = muOneF7TeV[i];
            fragGammaErrUp7TeV[i]       = muHalfF7TeV[i]-muOneF7TeV[i];
            fragGammaErrDown7TeV[i]     = muOneF7TeV[i]-muTwoF7TeV[i];
            promptGammaValue7TeV[i]     = muOneD7TeV[i];
            promptGammaErrUp7TeV[i]     = muHalfD7TeV[i]-muOneD7TeV[i];
            promptGammaErrDown7TeV[i]   = muOneD7TeV[i]-muTwoD7TeV[i];
        }    
        
        TGraphAsymmErrors *graphNLOCalcDirGam7TeV       = new TGraphAsymmErrors(nlinesNLOTwo7TeV,ptNLOPhotonTwo7TeV,gammaValue7TeV,ptNLOPhotonTwo7TeVErr,ptNLOPhotonTwo7TeVErr,gammaErrDown7TeV,gammaErrUp7TeV);
        graphNLOCalcDirGam7TeV->RemovePoint(0);
        TGraphAsymmErrors *graphNLOCalcFragGam7TeV      = new TGraphAsymmErrors(nlinesNLOTwo7TeV,ptNLOPhotonTwo7TeV,fragGammaValue7TeV,ptNLOPhotonTwo7TeVErr,ptNLOPhotonTwo7TeVErr,fragGammaErrDown7TeV,fragGammaErrUp7TeV);
        graphNLOCalcFragGam7TeV->RemovePoint(0);
        TGraphAsymmErrors *graphNLOCalcPromGam7TeV      = new TGraphAsymmErrors(nlinesNLOTwo7TeV,ptNLOPhotonTwo7TeV,promptGammaValue7TeV,ptNLOPhotonTwo7TeVErr,ptNLOPhotonTwo7TeVErr,promptGammaErrDown7TeV,promptGammaErrUp7TeV);
        graphNLOCalcPromGam7TeV->RemovePoint(0);
        
        TGraphAsymmErrors* graphRatioNLOFragGammaDivTot7TeV      = CalculateAsymGraphRatioToGraph(graphNLOCalcFragGam7TeV, graphNLOCalcDirGam7TeV);
        TGraphAsymmErrors* graphRatioNLOPromptGammaDivTot7TeV    = CalculateAsymGraphRatioToGraph(graphNLOCalcPromGam7TeV, graphNLOCalcDirGam7TeV);
        
        
        // -----------------------------------------------------------------------------------------------------------------------
        // ----------------------------- plotting fragmentation and prompt to total direct ---------------------------------------
        // -----------------------------------------------------------------------------------------------------------------------
        
        canvasRatioDirGammaCalc->cd();    
        histo2DRatioGammaCalc->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphRatioNLOFragGammaDivTot7TeV, 0, 0, colorFrag, colorFrag, 1, kTRUE, colorFrag);
            graphRatioNLOFragGammaDivTot7TeV->Draw("3,same");

            DrawGammaSetMarkerTGraphAsym(graphRatioNLOPromptGammaDivTot7TeV, 0, 0, colorPrompt, colorPrompt, 1, kTRUE, colorPrompt);
            graphRatioNLOPromptGammaDivTot7TeV->SetFillStyle(3245);
            graphRatioNLOPromptGammaDivTot7TeV->Draw("3,same");
            
            DrawGammaLines(0.23, 100. , 1., 1.,0.1, kGray+2);
            
            legendRatioGammaCalc->Draw();

            TLatex *labelRatioGammaCalc7TeV   = new TLatex(0.15,0.93,collisionSystem7TeV.Data());
            SetStyleTLatex( labelRatioGammaCalc7TeV, 0.85*32,4);
            labelRatioGammaCalc7TeV->SetTextFont(43);
            labelRatioGammaCalc7TeV->Draw();
            labelRatioGamma->Draw();
            
        canvasRatioDirGammaCalc->SaveAs(Form("%s/GammaNLOCalc_Separation_PP7TeV.%s",outputDir.Data(),suffix.Data()));

        // **********************************************************************************************************************
        // *********************************** direct photon calculations for 8TeV *******************************************
        // **********************************************************************************************************************
        
        TString fileNameNLOPhotonHalf8TeV    = "ExternalInput/Theory/ALICENLOcalcDirectPhoton8TeVmuhalf.dat";
        TString fileNameNLOPhotonOne8TeV     = "ExternalInput/Theory/ALICENLOcalcDirectPhoton8TeVmu.dat";
        TString fileNameNLOPhotonTwo8TeV     = "ExternalInput/Theory/ALICENLOcalcDirectPhoton8TeVtwomu.dat";
        Int_t nlinesNLOTwo8TeV              = 0;
        Int_t nlinesNLOOne8TeV              = 0;
        Int_t nlinesNLOHalf8TeV             = 0;
        Double_t ptNLOPhotonTwo8TeV[100];
        Double_t ptNLOPhotonTwo8TeVErr[100];
        Double_t ptNLOPhotonOne8TeV[100];
        Double_t ptNLOPhotonHalf8TeV[100];
        
        Double_t muHalfDF8TeV[100];
        Double_t muHalfD8TeV[100];
        Double_t muHalfF8TeV[100];
        Double_t muOneDF8TeV[100];
        Double_t muOneD8TeV[100];
        Double_t muOneF8TeV[100];
        Double_t muTwoDF8TeV[100];
        Double_t muTwoD8TeV[100];
        Double_t muTwoF8TeV[100];
        Double_t gammaValue8TeV[100];
        Double_t gammaErrUp8TeV[100];
        Double_t gammaErrDown8TeV[100];
        Double_t fragGammaValue8TeV[100];
        Double_t fragGammaErrUp8TeV[100];
        Double_t fragGammaErrDown8TeV[100];
        Double_t promptGammaValue8TeV[100];
        Double_t promptGammaErrUp8TeV[100];
        Double_t promptGammaErrDown8TeV[100];
        
        ifstream inHalf8TeV;
        inHalf8TeV.open(fileNameNLOPhotonHalf8TeV,ios_base::in);
        cout << fileNameNLOPhotonHalf8TeV << endl;
            
        while(!inHalf8TeV.eof()){
            nlinesNLOHalf8TeV++;
            //               Pt                              DirectPhoton           FragmentationPhoton         SumPhoton
            inHalf8TeV >> ptNLOPhotonHalf8TeV[nlinesNLOHalf8TeV] >> muHalfD8TeV[nlinesNLOHalf8TeV] >> muHalfF8TeV[nlinesNLOHalf8TeV] >> muHalfDF8TeV[nlinesNLOHalf8TeV];
            
        }
        inHalf8TeV.close();

        ifstream inOne8TeV;
        inOne8TeV.open(fileNameNLOPhotonOne8TeV,ios_base::in);
        cout << fileNameNLOPhotonOne8TeV << endl;

        while(!inOne8TeV.eof()){
            nlinesNLOOne8TeV++;
            inOne8TeV >> ptNLOPhotonOne8TeV[nlinesNLOOne8TeV] >> muOneD8TeV[nlinesNLOOne8TeV] >> muOneF8TeV[nlinesNLOOne8TeV] >> muOneDF8TeV[nlinesNLOOne8TeV];
        }
        inOne8TeV.close();

        ifstream inTwo8TeV;
        inTwo8TeV.open(fileNameNLOPhotonTwo8TeV,ios_base::in);
        cout << fileNameNLOPhotonTwo8TeV << endl;

        while(!inTwo8TeV.eof()){
            nlinesNLOTwo8TeV++;
            inTwo8TeV >> ptNLOPhotonTwo8TeV[nlinesNLOTwo8TeV] >> muTwoD8TeV[nlinesNLOTwo8TeV] >> muTwoF8TeV[nlinesNLOTwo8TeV] >> muTwoDF8TeV[nlinesNLOTwo8TeV];
            ptNLOPhotonTwo8TeVErr[nlinesNLOTwo8TeV] = 0.25;
        }
        inTwo8TeV.close();
        
        for (Int_t i = 0; i < nlinesNLOTwo8TeV; i++){
            cout << ptNLOPhotonHalf8TeV[i] << "\t" << muHalfDF8TeV[i] << "\t"<< ptNLOPhotonOne8TeV[i] <<"\t" << muOneDF8TeV[i] << "\t"<< ptNLOPhotonTwo8TeV[i] << "\t" << muTwoDF8TeV[i] << endl;
            gammaValue8TeV[i]           = muOneDF8TeV[i];
            gammaErrUp8TeV[i]           = muHalfDF8TeV[i]-muOneDF8TeV[i];
            gammaErrDown8TeV[i]         = muOneDF8TeV[i]-muTwoDF8TeV[i];
            fragGammaValue8TeV[i]       = muOneF8TeV[i];
            fragGammaErrUp8TeV[i]       = muHalfF8TeV[i]-muOneF8TeV[i];
            fragGammaErrDown8TeV[i]     = muOneF8TeV[i]-muTwoF8TeV[i];
            promptGammaValue8TeV[i]     = muOneD8TeV[i];
            promptGammaErrUp8TeV[i]     = muHalfD8TeV[i]-muOneD8TeV[i];
            promptGammaErrDown8TeV[i]   = muOneD8TeV[i]-muTwoD8TeV[i];
        }    
        
        TGraphAsymmErrors *graphNLOCalcDirGam8TeV       = new TGraphAsymmErrors(nlinesNLOTwo8TeV,ptNLOPhotonTwo8TeV,gammaValue8TeV,ptNLOPhotonTwo8TeVErr,ptNLOPhotonTwo8TeVErr,gammaErrDown8TeV,gammaErrUp8TeV);
        graphNLOCalcDirGam8TeV->RemovePoint(0);
        TGraphAsymmErrors *graphNLOCalcFragGam8TeV      = new TGraphAsymmErrors(nlinesNLOTwo8TeV,ptNLOPhotonTwo8TeV,fragGammaValue8TeV,ptNLOPhotonTwo8TeVErr,ptNLOPhotonTwo8TeVErr,fragGammaErrDown8TeV,fragGammaErrUp8TeV);
        graphNLOCalcFragGam8TeV->RemovePoint(0);
        TGraphAsymmErrors *graphNLOCalcPromGam8TeV      = new TGraphAsymmErrors(nlinesNLOTwo8TeV,ptNLOPhotonTwo8TeV,promptGammaValue8TeV,ptNLOPhotonTwo8TeVErr,ptNLOPhotonTwo8TeVErr,promptGammaErrDown8TeV,promptGammaErrUp8TeV);
        graphNLOCalcPromGam8TeV->RemovePoint(0);
        
        TGraphAsymmErrors* graphRatioNLOFragGammaDivTot8TeV      = CalculateAsymGraphRatioToGraph(graphNLOCalcFragGam8TeV, graphNLOCalcDirGam8TeV);
        TGraphAsymmErrors* graphRatioNLOPromptGammaDivTot8TeV    = CalculateAsymGraphRatioToGraph(graphNLOCalcPromGam8TeV, graphNLOCalcDirGam8TeV);
        
        
        // -----------------------------------------------------------------------------------------------------------------------
        // ----------------------------- plotting fragmentation and prompt to total direct ---------------------------------------
        // -----------------------------------------------------------------------------------------------------------------------
        
        canvasRatioDirGammaCalc->cd();    
        histo2DRatioGammaCalc->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphRatioNLOFragGammaDivTot8TeV, 0, 0, colorFrag, colorFrag, 1, kTRUE, colorFrag);
            graphRatioNLOFragGammaDivTot8TeV->Draw("3,same");

            DrawGammaSetMarkerTGraphAsym(graphRatioNLOPromptGammaDivTot8TeV, 0, 0, colorPrompt, colorPrompt, 1, kTRUE, colorPrompt);
            graphRatioNLOPromptGammaDivTot8TeV->SetFillStyle(3245);
            graphRatioNLOPromptGammaDivTot8TeV->Draw("3,same");
            
            DrawGammaLines(0.23, 100. , 1., 1.,0.1, kGray+2);
            
            legendRatioGammaCalc->Draw();

            TLatex *labelRatioGammaCalc8TeV   = new TLatex(0.15,0.93,collisionSystem8TeV.Data());
            SetStyleTLatex( labelRatioGammaCalc8TeV, 0.85*32,4);
            labelRatioGammaCalc8TeV->SetTextFont(43);
            labelRatioGammaCalc8TeV->Draw();
            labelRatioGamma->Draw();
            
        canvasRatioDirGammaCalc->SaveAs(Form("%s/GammaNLOCalc_Separation_PP8TeV.%s",outputDir.Data(),suffix.Data()));
        
        //******************************************************************************************************************
        //********************** Prompt photon parametrisation Paquett (private communication)******************************
        //******************************************************************************************************************
        // Yield
        TF1* fitTheoryDirectMcGill2760GeV       = new TF1 ("fitTheoryPromptPaquett",
                                                "[0]/64.*1e-9*TMath::Exp(16.20-3.94*TMath::Log(x)-0.269*TMath::Log(x)**2)*1./(0.865779/TMath::Power(4,0.0694875))",
                                                0.1, 50);  
        fitTheoryDirectMcGill2760GeV->SetParameter(0,1);

        //******************************************************************************************************************
        //************************************** Calculate inv. yield ******************************************************
        //******************************************************************************************************************
        TGraphAsymmErrors* graphNLOCalcInvYieldINT1DirGam2760GeV    = (TGraphAsymmErrors*)graphNLOCalcDirGam2760GeV->Clone("graphNLOCalcInvYieldINT1DirGam2760GeV");
        graphNLOCalcInvYieldINT1DirGam2760GeV                       = (TGraphAsymmErrors*)ScaleGraphAsym(graphNLOCalcInvYieldINT1DirGam2760GeV, 1/recalcBarn/ReturnCorrectXSection("2.76TeV", 0));
        TGraphAsymmErrors* graphNLOCalcInvYieldINT1PromGam2760GeV   = (TGraphAsymmErrors*)graphNLOCalcPromGam2760GeV->Clone("graphNLOCalcInvYieldINT1PromGam2760GeV");
        graphNLOCalcInvYieldINT1PromGam2760GeV                      = (TGraphAsymmErrors*)ScaleGraphAsym(graphNLOCalcInvYieldINT1PromGam2760GeV, 1/recalcBarn/ReturnCorrectXSection("2.76TeV", 0));
        TGraphAsymmErrors* graphNLOCalcInvYieldINT1FragGam2760GeV   = (TGraphAsymmErrors*)graphNLOCalcFragGam2760GeV->Clone("graphNLOCalcInvYieldINT1FragGam2760GeV");
        graphNLOCalcInvYieldINT1FragGam2760GeV                      = (TGraphAsymmErrors*)ScaleGraphAsym(graphNLOCalcInvYieldINT1FragGam2760GeV, 1/recalcBarn/ReturnCorrectXSection("2.76TeV", 0));
        
        TGraphAsymmErrors* graphNLOCalcInvYieldINT7DirGam2760GeV    = (TGraphAsymmErrors*)graphNLOCalcDirGam2760GeV->Clone("graphNLOCalcInvYieldINT7DirGam2760GeV");
        graphNLOCalcInvYieldINT7DirGam2760GeV                       = (TGraphAsymmErrors*)ScaleGraphAsym(graphNLOCalcInvYieldINT7DirGam2760GeV, 1/recalcBarn/ReturnCorrectXSection("2.76TeV", 1));
        TGraphAsymmErrors* graphNLOCalcInvYieldINT7PromGam2760GeV   = (TGraphAsymmErrors*)graphNLOCalcPromGam2760GeV->Clone("graphNLOCalcInvYieldINT7PromGam2760GeV");
        graphNLOCalcInvYieldINT7PromGam2760GeV                      = (TGraphAsymmErrors*)ScaleGraphAsym(graphNLOCalcInvYieldINT7PromGam2760GeV, 1/recalcBarn/ReturnCorrectXSection("2.76TeV", 1));
        TGraphAsymmErrors* graphNLOCalcInvYieldINT7FragGam2760GeV   = (TGraphAsymmErrors*)graphNLOCalcFragGam2760GeV->Clone("graphNLOCalcInvYieldINT7FragGam2760GeV");
        graphNLOCalcInvYieldINT7FragGam2760GeV                      = (TGraphAsymmErrors*)ScaleGraphAsym(graphNLOCalcInvYieldINT7FragGam2760GeV, 1/recalcBarn/ReturnCorrectXSection("2.76TeV", 1));

        // xSection 5TeV?
//        TGraphAsymmErrors* graphNLOCalcInvYieldDirGam5TeV       = (TGraphAsymmErrors*)graphNLOCalcDirGam5TeV->Clone("graphNLOCalcInvYieldDirGam5TeV");
//        graphNLOCalcInvYieldDirGam5TeV                          = (TGraphAsymmErrors*)ScaleGraphAsym(graphNLOCalcInvYieldDirGam5TeV, 1/recalcBarn/ReturnCorrectXSection("5TeV", 0));
//        TGraphAsymmErrors* graphNLOCalcInvYieldPromGam5TeV      = (TGraphAsymmErrors*)graphNLOCalcPromGam5TeV->Clone("graphNLOCalcInvYieldPromGam5TeV");
//        graphNLOCalcInvYieldPromGam5TeV                         = (TGraphAsymmErrors*)ScaleGraphAsym(graphNLOCalcInvYieldPromGam5TeV, 1/recalcBarn/ReturnCorrectXSection("5TeV", 0));
//        TGraphAsymmErrors* graphNLOCalcInvYieldFragGam5TeV      = (TGraphAsymmErrors*)graphNLOCalcFragGam5TeV->Clone("graphNLOCalcInvYieldFragGam5TeV");
//        graphNLOCalcInvYieldFragGam5TeV                         = (TGraphAsymmErrors*)ScaleGraphAsym(graphNLOCalcInvYieldFragGam5TeV, 1/recalcBarn/ReturnCorrectXSection("5TeV", 0));
        
        TGraphAsymmErrors* graphNLOCalcInvYieldINT1DirGam7TeV       = (TGraphAsymmErrors*)graphNLOCalcDirGam7TeV->Clone("graphNLOCalcInvYieldINT1DirGam7TeV");
        graphNLOCalcInvYieldINT1DirGam7TeV                          = (TGraphAsymmErrors*)ScaleGraphAsym(graphNLOCalcInvYieldINT1DirGam7TeV, 1/recalcBarn/ReturnCorrectXSection("7TeV", 0));
        TGraphAsymmErrors* graphNLOCalcInvYieldINT1PromGam7TeV      = (TGraphAsymmErrors*)graphNLOCalcPromGam7TeV->Clone("graphNLOCalcInvYieldINT1PromGam7TeV");
        graphNLOCalcInvYieldINT1PromGam7TeV                         = (TGraphAsymmErrors*)ScaleGraphAsym(graphNLOCalcInvYieldINT1PromGam7TeV, 1/recalcBarn/ReturnCorrectXSection("7TeV", 0));
        TGraphAsymmErrors* graphNLOCalcInvYieldINT1FragGam7TeV      = (TGraphAsymmErrors*)graphNLOCalcFragGam7TeV->Clone("graphNLOCalcInvYieldINT1FragGam7TeV");
        graphNLOCalcInvYieldINT1FragGam7TeV                         = (TGraphAsymmErrors*)ScaleGraphAsym(graphNLOCalcInvYieldINT1FragGam7TeV, 1/recalcBarn/ReturnCorrectXSection("7TeV", 0));

        TGraphAsymmErrors* graphNLOCalcInvYieldINT7DirGam7TeV       = (TGraphAsymmErrors*)graphNLOCalcDirGam7TeV->Clone("graphNLOCalcInvYieldINT7DirGam7TeV");
        graphNLOCalcInvYieldINT7DirGam7TeV                          = (TGraphAsymmErrors*)ScaleGraphAsym(graphNLOCalcInvYieldINT7DirGam7TeV, 1/recalcBarn/ReturnCorrectXSection("7TeV", 1));
        TGraphAsymmErrors* graphNLOCalcInvYieldINT7PromGam7TeV      = (TGraphAsymmErrors*)graphNLOCalcPromGam7TeV->Clone("graphNLOCalcInvYieldINT7PromGam7TeV");
        graphNLOCalcInvYieldINT7PromGam7TeV                         = (TGraphAsymmErrors*)ScaleGraphAsym(graphNLOCalcInvYieldINT7PromGam7TeV, 1/recalcBarn/ReturnCorrectXSection("7TeV", 1));
        TGraphAsymmErrors* graphNLOCalcInvYieldINT7FragGam7TeV      = (TGraphAsymmErrors*)graphNLOCalcFragGam7TeV->Clone("graphNLOCalcInvYieldINT7FragGam7TeV");
        graphNLOCalcInvYieldINT7FragGam7TeV                         = (TGraphAsymmErrors*)ScaleGraphAsym(graphNLOCalcInvYieldINT7FragGam7TeV, 1/recalcBarn/ReturnCorrectXSection("7TeV", 1));

        TGraphAsymmErrors* graphNLOCalcInvYieldDirGam8TeV           = (TGraphAsymmErrors*)graphNLOCalcDirGam8TeV->Clone("graphNLOCalcInvYieldDirGam8TeV");
        graphNLOCalcInvYieldDirGam8TeV                              = (TGraphAsymmErrors*)ScaleGraphAsym(graphNLOCalcInvYieldDirGam8TeV, 1/recalcBarn/ReturnCorrectXSection("8TeV", 1));
        TGraphAsymmErrors* graphNLOCalcInvYieldPromGam8TeV          = (TGraphAsymmErrors*)graphNLOCalcPromGam8TeV->Clone("graphNLOCalcInvYieldPromGam8TeV");
        graphNLOCalcInvYieldPromGam8TeV                             = (TGraphAsymmErrors*)ScaleGraphAsym(graphNLOCalcInvYieldPromGam8TeV, 1/recalcBarn/ReturnCorrectXSection("8TeV", 1));
        TGraphAsymmErrors* graphNLOCalcInvYieldFragGam8TeV          = (TGraphAsymmErrors*)graphNLOCalcFragGam8TeV->Clone("graphNLOCalcInvYieldFragGam8TeV");
        graphNLOCalcInvYieldFragGam8TeV                             = (TGraphAsymmErrors*)ScaleGraphAsym(graphNLOCalcInvYieldFragGam8TeV, 1/recalcBarn/ReturnCorrectXSection("8TeV", 1));

        
        //******************************************************************************************************************
        //************************************** Writing output for pp ***************************************************
        //******************************************************************************************************************
        TFile *fileTheoryGraphsPP   = new TFile("ExternalInput/Theory/TheoryCompilationPP.root","UPDATE");

        
            fileTheoryGraphsPP->mkdir("DirectPhoton");
            TDirectoryFile* directoryGamma = (TDirectoryFile*)fileTheoryGraphsPP->Get("DirectPhoton"); 
            fileTheoryGraphsPP->cd("DirectPhoton");

            graphNLOCalcDirGam2760GeV->GetYaxis()->SetTitle("#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )");
            graphNLOCalcDirGam2760GeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            graphNLOCalcDirGam2760GeV->Write("graphDirectPhotonNLOVogelsang_2760GeV",TObject::kOverwrite);
            graphNLOCalcInvYieldINT1DirGam2760GeV->GetYaxis()->SetTitle("#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}dy} (GeV^{-2}#it{c})");
            graphNLOCalcInvYieldINT1DirGam2760GeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            graphNLOCalcInvYieldINT1DirGam2760GeV->Write("graphDirectPhotonNLOVogelsangInvYieldINT1_2760GeV",TObject::kOverwrite);
            graphNLOCalcInvYieldINT7DirGam2760GeV->GetYaxis()->SetTitle("#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}dy} (GeV^{-2}#it{c})");
            graphNLOCalcInvYieldINT7DirGam2760GeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            graphNLOCalcInvYieldINT7DirGam2760GeV->Write("graphDirectPhotonNLOVogelsangInvYieldINT7_2760GeV",TObject::kOverwrite);
            graphNLOCalcPromGam2760GeV->GetYaxis()->SetTitle("#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )");
            graphNLOCalcPromGam2760GeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            graphNLOCalcPromGam2760GeV->Write("graphPromptPhotonNLOVogelsang_2760GeV",TObject::kOverwrite);
            graphNLOCalcInvYieldINT1PromGam2760GeV->GetYaxis()->SetTitle("#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}dy} (GeV^{-2}#it{c})");
            graphNLOCalcInvYieldINT1PromGam2760GeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            graphNLOCalcInvYieldINT1PromGam2760GeV->Write("graphPromptPhotonNLOVogelsangInvYieldINT1_2760GeV",TObject::kOverwrite);
            graphNLOCalcInvYieldINT7PromGam2760GeV->GetYaxis()->SetTitle("#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}dy} (GeV^{-2}#it{c})");
            graphNLOCalcInvYieldINT7PromGam2760GeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            graphNLOCalcInvYieldINT7PromGam2760GeV->Write("graphPromptPhotonNLOVogelsangInvYieldINT7_2760GeV",TObject::kOverwrite);
            graphNLOCalcFragGam2760GeV->GetYaxis()->SetTitle("#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )");
            graphNLOCalcFragGam2760GeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            graphNLOCalcFragGam2760GeV->Write("graphFragmentationPhotonNLOVogelsang_2760GeV",TObject::kOverwrite);
            graphNLOCalcInvYieldINT1FragGam2760GeV->GetYaxis()->SetTitle("#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}dy} (GeV^{-2}#it{c})");
            graphNLOCalcInvYieldINT1FragGam2760GeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            graphNLOCalcInvYieldINT1FragGam2760GeV->Write("graphFragmentationPhotonNLOVogelsangInvYieldINT1_2760GeV",TObject::kOverwrite);
            graphNLOCalcInvYieldINT7FragGam2760GeV->GetYaxis()->SetTitle("#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}dy} (GeV^{-2}#it{c})");
            graphNLOCalcInvYieldINT7FragGam2760GeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            graphNLOCalcInvYieldINT7FragGam2760GeV->Write("graphFragmentationPhotonNLOVogelsangInvYieldINT7_2760GeV",TObject::kOverwrite);
            graphRatioNLOFragGammaDivTot2760GeV->Write("graphFragPhotonDivDirectNLOVogelsang_2760GeV",TObject::kOverwrite);
            graphRatioNLOPromptGammaDivTot2760GeV->Write("graphPromptPhotonDivDirectNLOVogelsang_2760GeV",TObject::kOverwrite);
            fitTheoryDirectMcGill2760GeV->Write("fitYieldDirectPhotonNLOPaquett_2760GeV",TObject::kOverwrite);

            graphNLOCalcDirGam5TeV->GetYaxis()->SetTitle("#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )");
            graphNLOCalcDirGam5TeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");            
            graphNLOCalcDirGam5TeV->Write("graphDirectPhotonNLOVogelsang_5TeV",TObject::kOverwrite);
            graphNLOCalcPromGam5TeV->GetYaxis()->SetTitle("#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )");
            graphNLOCalcPromGam5TeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            graphNLOCalcPromGam5TeV->Write("graphPromptPhotonNLOVogelsang_5TeV",TObject::kOverwrite);
            graphNLOCalcFragGam5TeV->GetYaxis()->SetTitle("#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )");
            graphNLOCalcFragGam5TeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            graphNLOCalcFragGam5TeV->Write("graphFragmentationPhotonNLOVogelsang_5TeV",TObject::kOverwrite);
            graphRatioNLOFragGammaDivTot5TeV->Write("graphFragPhotonDivDirectNLOVogelsang_5TeV",TObject::kOverwrite);
            graphRatioNLOPromptGammaDivTot5TeV->Write("graphPromptPhotonDivDirectNLOVogelsang_5TeV",TObject::kOverwrite);
            
            graphNLOCalcDirGam7TeV->GetYaxis()->SetTitle("#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )");
            graphNLOCalcDirGam7TeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            graphNLOCalcDirGam7TeV->Write("graphDirectPhotonNLOVogelsang_7TeV",TObject::kOverwrite);
            graphNLOCalcInvYieldINT1DirGam7TeV->GetYaxis()->SetTitle("#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}dy} (GeV^{-2}#it{c})");
            graphNLOCalcInvYieldINT1DirGam7TeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            graphNLOCalcInvYieldINT1DirGam7TeV->Write("graphDirectPhotonNLOVogelsangInvYieldINT1_7TeV",TObject::kOverwrite);
            graphNLOCalcInvYieldINT7DirGam7TeV->GetYaxis()->SetTitle("#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}dy} (GeV^{-2}#it{c})");
            graphNLOCalcInvYieldINT7DirGam7TeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            graphNLOCalcInvYieldINT7DirGam7TeV->Write("graphDirectPhotonNLOVogelsangInvYieldINT7_7TeV",TObject::kOverwrite);
            graphNLOCalcPromGam7TeV->GetYaxis()->SetTitle("#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )");
            graphNLOCalcPromGam7TeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");            
            graphNLOCalcPromGam7TeV->Write("graphPromptPhotonNLOVogelsang_7TeV",TObject::kOverwrite);
            graphNLOCalcInvYieldINT1PromGam7TeV->GetYaxis()->SetTitle("#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}dy} (GeV^{-2}#it{c})");
            graphNLOCalcInvYieldINT1PromGam7TeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            graphNLOCalcInvYieldINT1PromGam7TeV->Write("graphPromptPhotonNLOVogelsangInvYieldINT1_7TeV",TObject::kOverwrite);
            graphNLOCalcInvYieldINT7PromGam7TeV->GetYaxis()->SetTitle("#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}dy} (GeV^{-2}#it{c})");
            graphNLOCalcInvYieldINT7PromGam7TeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            graphNLOCalcInvYieldINT7PromGam7TeV->Write("graphPromptPhotonNLOVogelsangInvYieldINT7_7TeV",TObject::kOverwrite);
            graphNLOCalcFragGam7TeV->GetYaxis()->SetTitle("#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )");
            graphNLOCalcFragGam7TeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");                        
            graphNLOCalcFragGam7TeV->Write("graphFragmentationPhotonNLOVogelsang_7TeV",TObject::kOverwrite);
            graphNLOCalcInvYieldINT1FragGam7TeV->GetYaxis()->SetTitle("#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}dy} (GeV^{-2}#it{c})");
            graphNLOCalcInvYieldINT1FragGam7TeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            graphNLOCalcInvYieldINT1FragGam7TeV->Write("graphFragmentationPhotonNLOVogelsangInvYieldINT1_7TeV",TObject::kOverwrite);
            graphNLOCalcInvYieldINT7FragGam7TeV->GetYaxis()->SetTitle("#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}dy} (GeV^{-2}#it{c})");
            graphNLOCalcInvYieldINT7FragGam7TeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            graphNLOCalcInvYieldINT7FragGam7TeV->Write("graphFragmentationPhotonNLOVogelsangInvYieldINT7_7TeV",TObject::kOverwrite);
            graphRatioNLOFragGammaDivTot7TeV->Write("graphFragPhotonDivDirectNLOVogelsang_7TeV",TObject::kOverwrite);
            graphRatioNLOPromptGammaDivTot7TeV->Write("graphPromptPhotonDivDirectNLOVogelsang_7TeV",TObject::kOverwrite);

            graphNLOCalcDirGam8TeV->GetYaxis()->SetTitle("#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )");
            graphNLOCalcDirGam8TeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");            
            graphNLOCalcDirGam8TeV->Write("graphDirectPhotonNLOVogelsang_8TeV",TObject::kOverwrite);
            graphNLOCalcInvYieldDirGam8TeV->GetYaxis()->SetTitle("#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}dy} (GeV^{-2}#it{c})");
            graphNLOCalcInvYieldDirGam8TeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            graphNLOCalcInvYieldDirGam8TeV->Write("graphDirectPhotonNLOVogelsangInvYield_8TeV",TObject::kOverwrite);
            graphNLOCalcPromGam8TeV->GetYaxis()->SetTitle("#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )");
            graphNLOCalcPromGam8TeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");            
            graphNLOCalcPromGam8TeV->Write("graphPromptPhotonNLOVogelsang_8TeV",TObject::kOverwrite);
            graphNLOCalcInvYieldPromGam8TeV->GetYaxis()->SetTitle("#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}dy} (GeV^{-2}#it{c})");
            graphNLOCalcInvYieldPromGam8TeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            graphNLOCalcInvYieldPromGam8TeV->Write("graphPromptPhotonNLOVogelsangInvYield_8TeV",TObject::kOverwrite);
            graphNLOCalcFragGam8TeV->GetYaxis()->SetTitle("#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )");
            graphNLOCalcFragGam8TeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");            
            graphNLOCalcFragGam8TeV->Write("graphFragmentationPhotonNLOVogelsang_8TeV",TObject::kOverwrite);
            graphNLOCalcInvYieldFragGam8TeV->GetYaxis()->SetTitle("#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{#it{p}_{T}d#it{p}_{T}dy} (GeV^{-2}#it{c})");
            graphNLOCalcInvYieldFragGam8TeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            graphNLOCalcInvYieldFragGam8TeV->Write("graphFragmentationPhotonNLOVogelsangInvYield_8TeV",TObject::kOverwrite);
            graphRatioNLOFragGammaDivTot8TeV->Write("graphFragPhotonDivDirectNLOVogelsang_8TeV",TObject::kOverwrite);
            graphRatioNLOPromptGammaDivTot8TeV->Write("graphPromptPhotonDivDirectNLOVogelsang_8TeV",TObject::kOverwrite);

        fileTheoryGraphsPP->Close();
        delete fileTheoryGraphsPP;
    }
    
    //******************************************************************************************************************
    //******************************* Put all direct photon calculations for PbPb in one file **************************
    //******************************************************************************************************************
    if (runPbPb){
        //******************************************************************************************************************
        //***************************** Direct Photon calculations *********************************************************
        //******************************************************************************************************************
        
        //******************************************************************************************************************
        //*********************************Old MC Gill *********************************************************************
        //******************************************************************************************************************
    //     Double_t theoryMcGillPbPb0020_x[15]     = {    2.000000000000000111e-01, 4.000000000000000222e-01, 5.999999999999999778e-01, 8.000000000000000444e-01, 1.000000000000000000e+00, 
    //                                                 1.199999999999999956e+00, 1.399999999999999911e+00, 1.600000000000000089e+00, 1.800000000000000044e+00, 2.000000000000000000e+00, 
    //                                                 2.200000000000000178e+00, 2.399999999999999911e+00, 2.600000000000000089e+00, 2.799999999999999822e+00, 3.000000000000000000e+00};
    //     Double_t theoryMcGillPbPb0020_xerr[15]     = {    0, 0, 0, 0, 0,
    //                                                 0, 0, 0, 0, 0,
    //                                                 0, 0, 0, 0, 0};
    //     Double_t theoryMcGillPbPb0020_y[15]     = {    1.540859385329232794e+02, 2.044997387338719363e+01, 5.878499923344167044e+00, 2.250958012412696885e+00, 9.991863147603015083e-01,
    //                                                 4.845094433818496471e-01, 2.497293734693046829e-01, 1.349175639294052931e-01, 7.576424621951352578e-02, 4.395437070531096890e-02,
    //                                                 2.638121048976970612e-02, 1.632128463217997344e-02, 1.037788835831024437e-02, 6.782615199120981507e-03, 4.547945132616882172e-03 }; 
    //     Double_t theoryMcGillPbPb0020_yerr[15]     = {    2.416453821976370708e+00, 3.301167073906846050e-01, 9.581685975141712719e-02, 3.707312930049495858e-02, 1.673623054763240595e-02, 
    //                                                 8.318531597347338796e-03, 4.417758588599256936e-03, 2.462353352140851659e-03, 1.423226057051350064e-03, 8.454263715662204250e-04, 
    //                                                 5.166073785903992467e-04, 3.233128079544251221e-04, 2.066024310156839574e-04, 1.349119265242332950e-04, 8.990109325466884435e-05}; 
    //     TGraphErrors* graphTheoryMcGill0020     = new TGraphErrors(15,theoryMcGillPbPb0020_x,theoryMcGillPbPb0020_y,theoryMcGillPbPb0020_xerr, theoryMcGillPbPb0020_yerr);  
    //     TGraphErrors* graphTheoryMcGill0020Plot    = ScaleGraph(graphTheoryMcGill0020,100);
    //     graphTheoryMcGill0020Plot->RemovePoint(0);                                           
    //     graphTheoryMcGill0020Plot->RemovePoint(0);
    //     graphTheoryMcGill0020Plot->RemovePoint(0);
    //     Double_t theoryMcGillPbPb2040_x[15]     = {    2.000000000000000111e-01, 4.000000000000000222e-01, 5.999999999999999778e-01, 8.000000000000000444e-01, 1.000000000000000000e+00, 
    //                                                 1.199999999999999956e+00, 1.399999999999999911e+00, 1.600000000000000089e+00, 1.800000000000000044e+00, 2.000000000000000000e+00, 
    //                                                 2.200000000000000178e+00, 2.399999999999999911e+00, 2.600000000000000089e+00, 2.799999999999999822e+00, 3.000000000000000000e+00}; 
    //     Double_t theoryMcGillPbPb2040_xerr[15]     = {    0, 0, 0, 0, 0,
    //                                                 0, 0, 0, 0, 0,
    //                                                 0, 0, 0, 0, 0};
    //     Double_t theoryMcGillPbPb2040_y[15]     = {    5.996162411828655081e+01, 7.682871584213885718e+00, 2.170336260707449672e+00, 8.182166469176795909e-01, 3.583236917664220922e-01,
    //                                                 1.718390778093310256e-01, 8.781549983209124832e-02, 4.716549515653831182e-02, 2.639893678394553828e-02, 1.530106534174231064e-02,
    //                                                 9.200643808420352898e-03, 5.715716866426879920e-03, 3.655916602971762304e-03, 2.407585539368926539e-03, 1.628794054128191544e-03 }; 
    //     Double_t theoryMcGillPbPb2040_yerr[15]     = {    1.293765778931855293e+00, 1.700623864527457396e-01, 4.903841383176194696e-02, 1.908670016811140485e-02, 8.698233566001163999e-03,
    //                                                 4.353428213951922483e-03, 2.319044963732616853e-03, 1.292520172768387397e-03, 7.456386913794923587e-04, 4.415830500157004066e-04, 
    //                                                 2.693008616739346649e-04, 1.683700730118350920e-04, 1.075702856825843344e-04, 7.032156976682658054e-05, 4.696434116417641080e-05 }; 
    //     TGraphErrors* graphTheoryMcGill2040     = new TGraphErrors(15,theoryMcGillPbPb2040_x,theoryMcGillPbPb2040_y,theoryMcGillPbPb2040_xerr, theoryMcGillPbPb2040_yerr);  
    //     TGraphErrors* graphTheoryMcGill2040Plot    = ScaleGraph(graphTheoryMcGill2040,10);
    //     graphTheoryMcGill2040Plot->RemovePoint(0);                                           
    //     graphTheoryMcGill2040Plot->RemovePoint(0);
    //     graphTheoryMcGill2040Plot->RemovePoint(0);
        
        
        //******************************************************************************************************************
        //*********************************MC Gill (Paquett, IPGlasma, 14.08.2015 - physi-mail) ****************************
        //******************************************************************************************************************    
        Double_t ptMCGill0020           [100];
        Double_t yieldMCGill0020        [100];
        Double_t errYieldMCGill0020     [100];
        Double_t errYieldXMCGill0020    [100];
        Int_t nlinesMCGill0020                      = 0;
        
        TString fileNameMCGill0020                  = "ExternalInputPbPb/Theory/McGill/direct_photons_LHC2760_cent0020_Paquet_McGill.dat";
        ifstream  fileMCGill0020;
        fileMCGill0020.open(fileNameMCGill0020,ios_base::in);
        cout << fileNameMCGill0020 << endl;
        
        while(!fileMCGill0020.eof() && nlinesMCGill0020< 100){
            fileMCGill0020 >> ptMCGill0020[nlinesMCGill0020] >> yieldMCGill0020[nlinesMCGill0020] >> errYieldMCGill0020[nlinesMCGill0020] ; 
            cout << nlinesMCGill0020 << "\t"  << ptMCGill0020[nlinesMCGill0020] << "\t"  << yieldMCGill0020[nlinesMCGill0020] << "\t"  << errYieldMCGill0020[nlinesMCGill0020] << endl;;
            errYieldMCGill0020[nlinesMCGill0020]    = 0;
            errYieldXMCGill0020[nlinesMCGill0020]   = 0;
            nlinesMCGill0020++;
        }
        fileMCGill0020.close();
        TGraphErrors* graphDirectPhotonMCGill0020   = new TGraphErrors(nlinesMCGill0020-1,ptMCGill0020,yieldMCGill0020, errYieldXMCGill0020, errYieldMCGill0020); 

        Double_t ptMCGillPrompt0020         [100];
        Double_t yieldMCGillPrompt0020      [100];
        Double_t errYieldMCGillPrompt0020   [100];
        Double_t errYieldXMCGillPrompt0020  [100];
        Int_t nlinesMCGillPrompt0020                = 0;
            
        TString fileNameMCGillPrompt0020            = "ExternalInputPbPb/Theory/McGill/prompt_photons_LHC2760_cent0020_Paquet_McGill.dat";
        ifstream  fileMCGillPrompt0020;
        fileMCGillPrompt0020.open(fileNameMCGillPrompt0020,ios_base::in);
        cout << fileNameMCGillPrompt0020 << endl;
        
        while(!fileMCGillPrompt0020.eof() && nlinesMCGillPrompt0020< 100){
            fileMCGillPrompt0020 >> ptMCGillPrompt0020[nlinesMCGillPrompt0020] >> yieldMCGillPrompt0020[nlinesMCGillPrompt0020] >> errYieldMCGillPrompt0020[nlinesMCGillPrompt0020] ; 
            cout << nlinesMCGillPrompt0020 << "\t"  << ptMCGillPrompt0020[nlinesMCGillPrompt0020] << "\t"  << yieldMCGillPrompt0020[nlinesMCGillPrompt0020] << "\t"  <<
            errYieldMCGillPrompt0020[nlinesMCGillPrompt0020] << endl;;
            errYieldMCGillPrompt0020[nlinesMCGillPrompt0020]        = 0;
            errYieldXMCGillPrompt0020[nlinesMCGillPrompt0020]       = 0;
            nlinesMCGillPrompt0020++;
        }
        fileMCGillPrompt0020.close();
        TGraphErrors* graphPromptPhotonMCGill0020   = new TGraphErrors(nlinesMCGillPrompt0020-1,ptMCGillPrompt0020,yieldMCGillPrompt0020, errYieldXMCGillPrompt0020, errYieldMCGillPrompt0020); 

        
        Double_t ptMCGill2040           [100];
        Double_t yieldMCGill2040        [100];
        Double_t errYieldMCGill2040     [100];
        Double_t errYieldXMCGill2040    [100];
        Int_t nlinesMCGill2040                      = 0;
        
        TString fileNameMCGill2040                  = "ExternalInputPbPb/Theory/McGill/direct_photons_LHC2760_cent2040_Paquet_McGill.dat";
        ifstream  fileMCGill2040;
        fileMCGill2040.open(fileNameMCGill2040,ios_base::in);
        cout << fileNameMCGill2040 << endl;
        
        while(!fileMCGill2040.eof() && nlinesMCGill2040< 100){
            fileMCGill2040 >> ptMCGill2040[nlinesMCGill2040] >> yieldMCGill2040[nlinesMCGill2040] >> errYieldMCGill2040[nlinesMCGill2040] ; 
            cout << nlinesMCGill2040 << "\t"  << ptMCGill2040[nlinesMCGill2040] << "\t"  << yieldMCGill2040[nlinesMCGill2040] << "\t"  << errYieldMCGill2040[nlinesMCGill2040] << endl;;
            errYieldMCGill2040[nlinesMCGill2040]    = 0;
            errYieldXMCGill2040[nlinesMCGill2040]   = 0;
            nlinesMCGill2040++;
        }
        fileMCGill2040.close();
        TGraphErrors* graphDirectPhotonMCGill2040   = new TGraphErrors(nlinesMCGill2040-1,ptMCGill2040,yieldMCGill2040, errYieldXMCGill2040, errYieldMCGill2040); 

        Double_t ptMCGillPrompt2040         [100];
        Double_t yieldMCGillPrompt2040      [100];
        Double_t errYieldMCGillPrompt2040   [100];
        Double_t errYieldXMCGillPrompt2040  [100];
        Int_t nlinesMCGillPrompt2040                = 0;
            
        TString fileNameMCGillPrompt2040            = "ExternalInputPbPb/Theory/McGill/prompt_photons_LHC2760_cent2040_Paquet_McGill.dat";
        ifstream  fileMCGillPrompt2040;
        fileMCGillPrompt2040.open(fileNameMCGillPrompt2040,ios_base::in);
        cout << fileNameMCGillPrompt2040 << endl;
        
        while(!fileMCGillPrompt2040.eof() && nlinesMCGillPrompt2040< 100){
            fileMCGillPrompt2040 >> ptMCGillPrompt2040[nlinesMCGillPrompt2040] >> yieldMCGillPrompt2040[nlinesMCGillPrompt2040] >> errYieldMCGillPrompt2040[nlinesMCGillPrompt2040] ; 
            cout << nlinesMCGillPrompt2040 << "\t"  << ptMCGillPrompt2040[nlinesMCGillPrompt2040] << "\t"  << yieldMCGillPrompt2040[nlinesMCGillPrompt2040] << "\t"  <<
            errYieldMCGillPrompt2040[nlinesMCGillPrompt2040] << endl;;
            errYieldMCGillPrompt2040[nlinesMCGillPrompt2040]    = 0;
            errYieldXMCGillPrompt2040[nlinesMCGillPrompt2040]     = 0;
            nlinesMCGillPrompt2040++;
        }
        fileMCGillPrompt2040.close();
        TGraphErrors* graphPromptPhotonMCGill2040   = new TGraphErrors(nlinesMCGillPrompt2040-1,ptMCGillPrompt2040,yieldMCGillPrompt2040, errYieldXMCGillPrompt2040, errYieldMCGillPrompt2040); 

        
        //******************************************************************************************************************
        //*********************************PHSD arXiv:1504.05699 [nucl-th], O. Linnyk **************************************
        //******************************************************************************************************************
        Double_t ptPHSD             [100];
        Double_t yieldPHSD0020      [100];
        Double_t yieldPHSD0040      [100];
        Double_t yieldPHSD2040      [100];
        Double_t yieldPHSD4080      [100];
        Double_t errPHSD0020        [100];
        Double_t errPHSD0040        [100];
        Double_t errPHSD2040        [100];
        Double_t errPHSD4080        [100];
        Double_t errXPHSD           [100];
        
        Int_t nlinesPHSD                            = 0;
        
        TString fileNamePHSD                        = "ExternalInputPbPb/Theory/PHSD/Direct_Spcetra_Centralities_prepForRead.dat";
        ifstream  filePHSD;
        filePHSD.open(fileNamePHSD,ios_base::in);
        cout << fileNamePHSD << endl;
        
        while(!filePHSD.eof() && nlinesPHSD< 100){
            filePHSD >> ptPHSD[nlinesPHSD] >> yieldPHSD0020[nlinesPHSD] >> yieldPHSD2040[nlinesPHSD] >> yieldPHSD0040[nlinesPHSD]>> yieldPHSD4080[nlinesPHSD]; 
            cout << nlinesPHSD << "\t"  << ptPHSD[nlinesPHSD] << "\t"  << yieldPHSD0020[nlinesPHSD]<< "\t"  << yieldPHSD2040[nlinesPHSD] << "\t"  << yieldPHSD4080[nlinesPHSD]<< endl;;
            errPHSD0020[nlinesPHSD]                 = 0;//yieldPHSD0020[nlinesPHSD]*0.3;
            errPHSD0040[nlinesPHSD]                 = 0;//yieldPHSD0040[nlinesPHSD]*0.3;
            errPHSD2040[nlinesPHSD]                 = 0;//yieldPHSD2040[nlinesPHSD]*0.3;
            errPHSD4080[nlinesPHSD]                 = 0;//yieldPHSD4080[nlinesPHSD]*0.3;
            errXPHSD[nlinesPHSD]                    = 0;
            nlinesPHSD++;
        }
        filePHSD.close();
        TGraphErrors* graphDirectPhotonPHSD0020     = new TGraphErrors(nlinesPHSD-1,ptPHSD,yieldPHSD0020, errXPHSD, errPHSD0020); 
        TGraphErrors* graphDirectPhotonPHSD2040     = new TGraphErrors(nlinesPHSD-1,ptPHSD,yieldPHSD2040, errXPHSD, errPHSD2040); 
        TGraphErrors* graphDirectPhotonPHSD0040     = new TGraphErrors(nlinesPHSD-1,ptPHSD,yieldPHSD0040, errXPHSD, errPHSD0040); 
        TGraphErrors* graphDirectPhotonPHSD4080     = new TGraphErrors(nlinesPHSD-1,ptPHSD,yieldPHSD4080, errXPHSD, errPHSD4080); 

        //******************************************************************************************************************
        //*********************************Fireball He et al. PRC85(2012)044911 ********************************************
        //******************************************************************************************************************
        Double_t ptHe               [100];
        Double_t yieldHe0020        [100];
    //     Double_t yieldHe0040     [100];
        Double_t yieldHe2040        [100];
        Double_t yieldHe4080        [100];
        Double_t errHe0020          [100];
    //     Double_t errHe0040       [100];
        Double_t errHe2040          [100];
        Double_t errHe4080          [100];
        Double_t errXHe             [100];
        
        Int_t nlinesHe                              = 0;
        
        TString fileNameHe                          = "ExternalInputPbPb/Theory/Rapp/YieldFireball.txt";
        ifstream  fileHe;
        fileHe.open(fileNameHe,ios_base::in);
        cout << fileNameHe << endl;
        
        while(!fileHe.eof() && nlinesHe< 100){
            fileHe >> ptHe[nlinesHe] >> yieldHe0020[nlinesHe] >> yieldHe2040[nlinesHe] >> yieldHe4080[nlinesHe]; 
            cout << nlinesHe << "\t"  << ptHe[nlinesHe] << "\t"  << yieldHe0020[nlinesHe]<< "\t"  << yieldHe2040[nlinesHe] << "\t"  << yieldHe4080[nlinesHe]<< endl;;
            errHe0020[nlinesHe]                     = 0;//yieldHe0020[nlinesHe]*0.3;
    //         errHe0040[nlinesHe]                  = 0;//yieldHe0040[nlinesHe]*0.3;
            errHe2040[nlinesHe]                     = 0;//yieldHe2040[nlinesHe]*0.3;
            errHe4080[nlinesHe]                     = 0;//yieldHe4080[nlinesHe]*0.3;
            errXHe[nlinesHe]                        = 0;
            nlinesHe++;
        }
        fileHe.close();
        TGraphErrors* graphDirectPhotonHe0020       = new TGraphErrors(nlinesHe-1,ptHe,yieldHe0020, errXHe, errHe0020); 
        while (graphDirectPhotonHe0020->GetX()[0] < 0.7) graphDirectPhotonHe0020->RemovePoint(0);

        TGraphErrors* graphDirectPhotonHe2040       = new TGraphErrors(nlinesHe-1,ptHe,yieldHe2040, errXHe, errHe2040); 
        while (graphDirectPhotonHe2040->GetX()[0] < 0.7) graphDirectPhotonHe2040->RemovePoint(0);

    //     TGraphErrors* graphDirectPhotonHe0040 = new TGraphErrors(nlinesHe-1,ptHe,yieldHe0040, errXHe, errHe0040); 
        TGraphErrors* graphDirectPhotonHe4080       = new TGraphErrors(nlinesHe-1,ptHe,yieldHe4080, errXHe, errHe4080); 
        while (graphDirectPhotonHe4080->GetX()[0] < 0.7) graphDirectPhotonHe4080->RemovePoint(0);
        
        //******************************************************************************************************************
        //*********************************Chatterjee http://journals.aps.org/prc/pdf/10.1103/PhysRevC.85.064910 ***********
        //******************************************************************************************************************
        Double_t ptChatterjee0020               [100];
        Double_t yieldThermalChatterjee0020     [100];
        Double_t yieldPromptChatterjee0020      [100];
        Double_t yieldChatterjee0020            [100];
        Double_t errYieldChatterjee0020         [100];
        Double_t errYieldXChatterjee0020        [100];
        Int_t nlinesChatterjee0020                  = 0;
        
        TString fileNameChatterjee0020              = "ExternalInputPbPb/Theory/Chatterjee/direct_0_20_forReading.dat";
        ifstream  fileChatterjee0020;
        fileChatterjee0020.open(fileNameChatterjee0020,ios_base::in);
        cout << fileNameChatterjee0020 << endl;
        
        while(!fileChatterjee0020.eof() && nlinesChatterjee0020< 100){
            fileChatterjee0020 >> ptChatterjee0020[nlinesChatterjee0020] >> yieldThermalChatterjee0020[nlinesChatterjee0020] >> yieldPromptChatterjee0020[nlinesChatterjee0020] 
                            >> yieldChatterjee0020[nlinesChatterjee0020] ; 
            cout << nlinesChatterjee0020 << "\t"  << ptChatterjee0020[nlinesChatterjee0020] << "\t"  << yieldThermalChatterjee0020[nlinesChatterjee0020] << "\t"  
                << yieldPromptChatterjee0020[nlinesChatterjee0020] << "\t"  << yieldChatterjee0020[nlinesChatterjee0020] << endl;;
            errYieldChatterjee0020[nlinesChatterjee0020]    = 0;
            errYieldXChatterjee0020[nlinesChatterjee0020]   = 0;
            nlinesChatterjee0020++;
        }
        fileChatterjee0020.close();
        TGraphErrors* graphDirectPhotonChatterjee0020           = new TGraphErrors(nlinesChatterjee0020-1,ptChatterjee0020,yieldChatterjee0020, errYieldXChatterjee0020, errYieldChatterjee0020); 
        
        TGraphErrors* graphDirectPhotonThermalChatterjee0020    = new TGraphErrors(nlinesChatterjee0020-1,ptChatterjee0020,yieldThermalChatterjee0020, errYieldXChatterjee0020, errYieldChatterjee0020); 
        while (graphDirectPhotonThermalChatterjee0020->GetY()[graphDirectPhotonThermalChatterjee0020->GetN()-1] == 0) 
            graphDirectPhotonThermalChatterjee0020->RemovePoint(graphDirectPhotonThermalChatterjee0020->GetN()-1); 
        
        TGraphErrors* graphDirectPhotonPromptChatterjee0020     = new TGraphErrors(nlinesChatterjee0020-1,ptChatterjee0020,yieldPromptChatterjee0020, errYieldXChatterjee0020, errYieldChatterjee0020); 
        while (graphDirectPhotonPromptChatterjee0020->GetY()[0] == 0) 
            graphDirectPhotonPromptChatterjee0020->RemovePoint(0);

        Double_t ptChatterjee2040               [100];
        Double_t yieldThermalChatterjee2040     [100];
        Double_t yieldPromptChatterjee2040      [100];
        Double_t yieldChatterjee2040            [100];
        Double_t errYieldChatterjee2040         [100];
        Double_t errYieldXChatterjee2040        [100];
        Int_t nlinesChatterjee2040                          = 0;
        
        TString fileNameChatterjee2040                      = "ExternalInputPbPb/Theory/Chatterjee/direct_20_40_forReading.dat";
        ifstream  fileChatterjee2040;
        fileChatterjee2040.open(fileNameChatterjee2040,ios_base::in);
        cout << fileNameChatterjee2040 << endl;
        
        while(!fileChatterjee2040.eof() && nlinesChatterjee2040< 100){
            fileChatterjee2040 >> ptChatterjee2040[nlinesChatterjee2040] >> yieldThermalChatterjee2040[nlinesChatterjee2040] >> yieldPromptChatterjee2040[nlinesChatterjee2040] 
                            >> yieldChatterjee2040[nlinesChatterjee2040] ; 
            cout << nlinesChatterjee2040 << "\t"  << ptChatterjee2040[nlinesChatterjee2040] << "\t"  << yieldThermalChatterjee2040[nlinesChatterjee2040] << "\t"  
                << yieldPromptChatterjee2040[nlinesChatterjee2040] << "\t"  << yieldChatterjee2040[nlinesChatterjee2040] << endl;;
            errYieldChatterjee2040[nlinesChatterjee2040]    = 0;
            errYieldXChatterjee2040[nlinesChatterjee2040]   = 0;
            nlinesChatterjee2040++;
        }
        fileChatterjee2040.close();
        TGraphErrors* graphDirectPhotonChatterjee2040           = new TGraphErrors(nlinesChatterjee2040-1,ptChatterjee2040,yieldChatterjee2040, errYieldXChatterjee2040, errYieldChatterjee2040); 
        
        TGraphErrors* graphDirectPhotonThermalChatterjee2040    = new TGraphErrors(nlinesChatterjee2040-1,ptChatterjee2040,yieldThermalChatterjee2040, errYieldXChatterjee2040, errYieldChatterjee2040); 
        while (graphDirectPhotonThermalChatterjee2040->GetY()[graphDirectPhotonThermalChatterjee2040->GetN()-1] == 0) 
            graphDirectPhotonThermalChatterjee2040->RemovePoint(graphDirectPhotonThermalChatterjee2040->GetN()-1); 
        
        TGraphErrors* graphDirectPhotonPromptChatterjee2040     = new TGraphErrors(nlinesChatterjee2040-1,ptChatterjee2040,yieldPromptChatterjee2040, errYieldXChatterjee2040, errYieldChatterjee2040); 
        while (graphDirectPhotonPromptChatterjee2040->GetY()[0] == 0) 
            graphDirectPhotonPromptChatterjee2040->RemovePoint(0);

        Double_t ptChatterjee4060               [100];
        Double_t yieldThermalChatterjee4060     [100];
        Double_t yieldPromptChatterjee4060      [100];
        Double_t yieldChatterjee4060            [100];
        Double_t errYieldChatterjee4060         [100];
        Double_t errYieldXChatterjee4060        [100];
        Int_t nlinesChatterjee4060                          = 0;
        
        TString fileNameChatterjee4060                      = "ExternalInputPbPb/Theory/Chatterjee/direct_40_60_forReading.dat";
        ifstream  fileChatterjee4060;
        fileChatterjee4060.open(fileNameChatterjee4060,ios_base::in);
        cout << fileNameChatterjee4060 << endl;
        
        while(!fileChatterjee4060.eof() && nlinesChatterjee4060< 100){
            fileChatterjee4060 >> ptChatterjee4060[nlinesChatterjee4060] >> yieldThermalChatterjee4060[nlinesChatterjee4060] >> yieldPromptChatterjee4060[nlinesChatterjee4060] 
                            >> yieldChatterjee4060[nlinesChatterjee4060] ; 
            cout << nlinesChatterjee4060 << "\t"  << ptChatterjee4060[nlinesChatterjee4060] << "\t"  << yieldThermalChatterjee4060[nlinesChatterjee4060] << "\t"  
                << yieldPromptChatterjee4060[nlinesChatterjee4060] << "\t"  << yieldChatterjee4060[nlinesChatterjee4060] << endl;;
            errYieldChatterjee4060[nlinesChatterjee4060]    = 0;
            errYieldXChatterjee4060[nlinesChatterjee4060]   = 0;
            nlinesChatterjee4060++;
        }
        fileChatterjee4060.close();
        TGraphErrors* graphDirectPhotonChatterjee4060           = new TGraphErrors(nlinesChatterjee4060-1,ptChatterjee4060,yieldChatterjee4060, errYieldXChatterjee4060, errYieldChatterjee4060); 
        
        TGraphErrors* graphDirectPhotonThermalChatterjee4060    = new TGraphErrors(nlinesChatterjee4060-1,ptChatterjee4060,yieldThermalChatterjee4060, errYieldXChatterjee4060, errYieldChatterjee4060); 
        while (graphDirectPhotonThermalChatterjee4060->GetY()[graphDirectPhotonThermalChatterjee4060->GetN()-1] == 0) 
            graphDirectPhotonThermalChatterjee4060->RemovePoint(graphDirectPhotonThermalChatterjee4060->GetN()-1); 
        
        TGraphErrors* graphDirectPhotonPromptChatterjee4060     = new TGraphErrors(nlinesChatterjee4060-1,ptChatterjee4060,yieldPromptChatterjee4060, errYieldXChatterjee4060, errYieldChatterjee4060); 
        while (graphDirectPhotonPromptChatterjee4060->GetY()[0] == 0) 
            graphDirectPhotonPromptChatterjee4060->RemovePoint(0);

        //******************************************************************************************************************
        //*********************************Chatterjee thermal from Rupa, pQCD from JHEP 1305 (2013) 030 -0-20% *************
        //******************************************************************************************************************
        Double_t ptChatterjee0020_2[6]                = { 1.3, 2, 3, 4, 5,
                                                            6};
        Double_t yieldThermalChatterjee0020_2[6]      = { 0.231531, 0.0343742, 0.00357161, 0.000518726, 9.30806e-05,
                                                            1.91335e-05};
        Double_t yieldPromptChatterjee0020_2[6]       = { 0.0831538, 0.0160724, 0.00253054, 0.00063178, 0.000216407,
                                                            8.72282e-05};
        Double_t yieldChatterjee0020_2[6]             = { 0.314685, 0.0504466, 0.00610215, 0.00115051, 0.000309488,
                                                            0.000106362};
        Double_t errYieldXChatterjee0020_2[6]         = { 0, 0, 0, 0, 0,
                                                            0  };
        Double_t errYieldYChatterjee0020_2[6]         = { 0, 0, 0, 0, 0,
                                                            0  };
        Int_t nlinesChatterjee0020_2                  = 6;

        TGraphErrors* graphDirectPhotonChatterjee0020_2         = new TGraphErrors(nlinesChatterjee0020_2,ptChatterjee0020_2,yieldChatterjee0020_2, errYieldXChatterjee0020_2, errYieldYChatterjee0020_2); 
        TGraphErrors* graphDirectPhotonThermalChatterjee0020_2  = new TGraphErrors(nlinesChatterjee0020_2,ptChatterjee0020_2,yieldThermalChatterjee0020_2, errYieldXChatterjee0020_2, errYieldYChatterjee0020_2);   
        TGraphErrors* graphDirectPhotonPromptChatterjee0020_2   = new TGraphErrors(nlinesChatterjee0020_2,ptChatterjee0020_2,yieldPromptChatterjee0020_2, errYieldXChatterjee0020_2, errYieldYChatterjee0020_2); 

        //******************************************************************************************************************
        //*********************************Chatterjee thermal from Rupa, pQCD from JHEP 1305 (2013) 030 -20-40% ************
        //******************************************************************************************************************
        Double_t ptChatterjee2040_2[6]                = { 1.3, 2, 3, 4, 5,
                                                            6};
        Double_t yieldThermalChatterjee2040_2[6]      = { 0.0920386, 0.013742, 0.00140182, 0.000194282, 3.29015e-05,
                                                            6.38523e-06};
        Double_t yieldPromptChatterjee2040_2[6]       = { 0.0302371, 0.00574401, 0.000894192, 0.000222269, 7.58425e-05,
                                                            3.05298e-05};
        Double_t yieldChatterjee2040_2[6]             = { 0.122276, 0.019486, 0.00229601, 0.000416551, 0.000108744,
                                                            3.6915e-05};
        Double_t errYieldXChatterjee2040_2[6]         = { 0, 0, 0, 0, 0,
                                                            0  };
        Double_t errYieldYChatterjee2040_2[6]         = { 0, 0, 0, 0, 0,
                                                            0  };
        Int_t nlinesChatterjee2040_2                  = 6;                                               

        TGraphErrors* graphDirectPhotonChatterjee2040_2         = new TGraphErrors(nlinesChatterjee2040_2, ptChatterjee2040_2, yieldChatterjee2040_2, 
                                                                                errYieldXChatterjee2040_2, errYieldYChatterjee2040_2); 
        TGraphErrors* graphDirectPhotonThermalChatterjee2040_2  = new TGraphErrors(nlinesChatterjee2040_2, ptChatterjee2040_2, yieldThermalChatterjee2040_2, 
                                                                                errYieldXChatterjee2040_2, errYieldYChatterjee2040_2);   
        TGraphErrors* graphDirectPhotonPromptChatterjee2040_2   = new TGraphErrors(nlinesChatterjee2040_2, ptChatterjee2040_2, yieldPromptChatterjee2040_2, 
                                                                                errYieldXChatterjee2040_2, errYieldYChatterjee2040_2); 
        
        //******************************************************************************************************************
        //*********************************v.Hees, Rapp NPA933(2015)256 ****************************************************
        //******************************************************************************************************************

        Double_t ptHees0020                 [100];
        Double_t yieldRhoSF0020             [100];
        Double_t yieldQGPHees0020           [100];
        Double_t yieldOmegaHees0020         [100];
        Double_t yieldMesonGasHees0020      [100];
        Double_t yieldPrimordialHees0020    [100];
        Double_t yieldHees0020              [100];
        Double_t errYieldHees0020           [100];
        Double_t errYieldXHees0020          [100];
        Int_t nlinesHees0020                        = 0;
        
        TString fileNameHees0020                    = "ExternalInputPbPb/Theory/Rapp/Yield0020.txt";
        ifstream  fileHees0020;
        fileHees0020.open(fileNameHees0020,ios_base::in);
        cout << fileNameHees0020 << endl;
        
        while(!fileHees0020.eof() && nlinesHees0020< 100){
            fileHees0020 >> ptHees0020[nlinesHees0020] >> yieldRhoSF0020[nlinesHees0020] >>yieldQGPHees0020[nlinesHees0020] >>yieldOmegaHees0020[nlinesHees0020] >> yieldMesonGasHees0020[nlinesHees0020]
                        >> yieldPrimordialHees0020[nlinesHees0020] >> yieldHees0020[nlinesHees0020] ; 
            cout << nlinesHees0020 << "\t"  << ptHees0020[nlinesHees0020] << "\t"  << yieldQGPHees0020[nlinesHees0020] << "\t"  
                << yieldPrimordialHees0020[nlinesHees0020] << "\t"  << yieldHees0020[nlinesHees0020] << endl;;
            errYieldHees0020[nlinesHees0020]        = 0;
            errYieldXHees0020[nlinesHees0020]       = 0;
            nlinesHees0020++;
        }
        fileHees0020.close();
        TGraphErrors* graphDirectPhotonHees0020             = new TGraphErrors(nlinesHees0020-1,ptHees0020,yieldHees0020, errYieldXHees0020, errYieldHees0020); 
        while (graphDirectPhotonHees0020->GetX()[0] < 0.7) graphDirectPhotonHees0020->RemovePoint(0);
        TGraphErrors* graphDirectPhotonRhoSFHees0020        = new TGraphErrors(nlinesHees0020-1,ptHees0020,yieldRhoSF0020, errYieldXHees0020, errYieldHees0020);
        TGraphErrors* graphDirectPhotonQGPHees0020          = new TGraphErrors(nlinesHees0020-1,ptHees0020,yieldQGPHees0020, errYieldXHees0020, errYieldHees0020);
        TGraphErrors* graphDirectPhotonOmegaHees0020        = new TGraphErrors(nlinesHees0020-1,ptHees0020,yieldOmegaHees0020, errYieldXHees0020, errYieldHees0020);
        TGraphErrors* graphDirectPhotonMesonGasHees0020     = new TGraphErrors(nlinesHees0020-1,ptHees0020,yieldMesonGasHees0020, errYieldXHees0020, errYieldHees0020);
        TGraphErrors* graphDirectPhotonPrimordialHees0020   = new TGraphErrors(nlinesHees0020-1,ptHees0020,yieldPrimordialHees0020, errYieldXHees0020, errYieldHees0020); 
        
        Double_t ptHees2040                 [100];
        Double_t yieldRhoSF2040             [100];
        Double_t yieldQGPHees2040           [100];
        Double_t yieldOmegaHees2040         [100];
        Double_t yieldMesonGasHees2040      [100];
        Double_t yieldPrimordialHees2040    [100];
        Double_t yieldHees2040              [100];
        Double_t errYieldHees2040           [100];
        Double_t errYieldXHees2040          [100];
        Int_t nlinesHees2040                        = 0;
        
        TString fileNameHees2040                    = "ExternalInputPbPb/Theory/Rapp/Yield2040.txt";
        ifstream  fileHees2040;
        fileHees2040.open(fileNameHees2040,ios_base::in);
        cout << fileNameHees2040 << endl;
        
        while(!fileHees2040.eof() && nlinesHees2040< 100){
            fileHees2040 >> ptHees2040[nlinesHees2040] >> yieldRhoSF2040[nlinesHees2040] >>yieldQGPHees2040[nlinesHees2040] >>yieldOmegaHees2040[nlinesHees2040] >> yieldMesonGasHees2040[nlinesHees2040]
                        >> yieldPrimordialHees2040[nlinesHees2040] >> yieldHees2040[nlinesHees2040] ; 
            cout << nlinesHees2040 << "\t"  << ptHees2040[nlinesHees2040] << "\t"  << yieldQGPHees2040[nlinesHees2040] << "\t"  
                << yieldPrimordialHees2040[nlinesHees2040] << "\t"  << yieldHees2040[nlinesHees2040] << endl;;
            errYieldHees2040[nlinesHees2040]        = 0;
            errYieldXHees2040[nlinesHees2040]       = 0;
            nlinesHees2040++;
        }
        fileHees2040.close();
        TGraphErrors* graphDirectPhotonHees2040             = new TGraphErrors(nlinesHees2040-1,ptHees2040,yieldHees2040, errYieldXHees2040, errYieldHees2040); 
        while (graphDirectPhotonHees2040->GetX()[0] < 0.7) graphDirectPhotonHees2040->RemovePoint(0);
        TGraphErrors* graphDirectPhotonRhoSFHees2040        = new TGraphErrors(nlinesHees2040-1,ptHees2040,yieldRhoSF2040, errYieldXHees2040, errYieldHees2040);
        TGraphErrors* graphDirectPhotonQGPHees2040          = new TGraphErrors(nlinesHees2040-1,ptHees2040,yieldQGPHees2040, errYieldXHees2040, errYieldHees2040);
        TGraphErrors* graphDirectPhotonOmegaHees2040        = new TGraphErrors(nlinesHees2040-1,ptHees2040,yieldOmegaHees2040, errYieldXHees2040, errYieldHees2040);
        TGraphErrors* graphDirectPhotonMesonGasHees2040     = new TGraphErrors(nlinesHees2040-1,ptHees2040,yieldMesonGasHees2040, errYieldXHees2040, errYieldHees2040);
        TGraphErrors* graphDirectPhotonPrimordialHees2040   = new TGraphErrors(nlinesHees2040-1,ptHees2040,yieldPrimordialHees2040, errYieldXHees2040, errYieldHees2040); 

        Double_t ptHees4080                 [100];
        Double_t yieldRhoSF4080             [100];
        Double_t yieldQGPHees4080           [100];
        Double_t yieldOmegaHees4080         [100];
        Double_t yieldMesonGasHees4080      [100];
        Double_t yieldPrimordialHees4080    [100];
        Double_t yieldHees4080              [100];
        Double_t errYieldHees4080           [100];
        Double_t errYieldXHees4080          [100];
        Int_t nlinesHees4080                        = 0;
        
        TString fileNameHees4080                    = "ExternalInputPbPb/Theory/Rapp/Yield4080.txt";
        ifstream  fileHees4080;
        fileHees4080.open(fileNameHees4080,ios_base::in);
        cout << fileNameHees4080 << endl;
        
        while(!fileHees4080.eof() && nlinesHees4080< 100){
            fileHees4080 >> ptHees4080[nlinesHees4080] >> yieldRhoSF4080[nlinesHees4080] >>yieldQGPHees4080[nlinesHees4080] >>yieldOmegaHees4080[nlinesHees4080] >> yieldMesonGasHees4080[nlinesHees4080]
                        >> yieldPrimordialHees4080[nlinesHees4080] >> yieldHees4080[nlinesHees4080] ; 
            cout << nlinesHees4080 << "\t"  << ptHees4080[nlinesHees4080] << "\t"  << yieldQGPHees4080[nlinesHees4080] << "\t"  
                << yieldPrimordialHees4080[nlinesHees4080] << "\t"  << yieldHees4080[nlinesHees4080] << endl;;
            errYieldHees4080[nlinesHees4080]        = 0;
            errYieldXHees4080[nlinesHees4080]       = 0;
            nlinesHees4080++;
        }
        fileHees4080.close();
        TGraphErrors* graphDirectPhotonHees4080             = new TGraphErrors(nlinesHees4080-1,ptHees4080,yieldHees4080, errYieldXHees4080, errYieldHees4080); 
        while (graphDirectPhotonHees4080->GetX()[0] < 0.7) graphDirectPhotonHees4080->RemovePoint(0);
        TGraphErrors* graphDirectPhotonRhoSFHees4080        = new TGraphErrors(nlinesHees4080-1,ptHees4080,yieldRhoSF4080, errYieldXHees4080, errYieldHees4080);
        TGraphErrors* graphDirectPhotonQGPHees4080          = new TGraphErrors(nlinesHees4080-1,ptHees4080,yieldQGPHees4080, errYieldXHees4080, errYieldHees4080);
        TGraphErrors* graphDirectPhotonOmegaHees4080        = new TGraphErrors(nlinesHees4080-1,ptHees4080,yieldOmegaHees4080, errYieldXHees4080, errYieldHees4080);
        TGraphErrors* graphDirectPhotonMesonGasHees4080     = new TGraphErrors(nlinesHees4080-1,ptHees4080,yieldMesonGasHees4080, errYieldXHees4080, errYieldHees4080);
        TGraphErrors* graphDirectPhotonPrimordialHees4080   = new TGraphErrors(nlinesHees4080-1,ptHees4080,yieldPrimordialHees4080, errYieldXHees4080, errYieldHees4080); 
        
        //******************************************************************************************************************
        //*********************************Holopainen **********************************************************************
        //******************************************************************************************************************
        
        Double_t ptHolopainen0020           [100];
        Double_t v2Holopainen0020           [100];
        Double_t yieldHolopainen0020        [100];
        Double_t errYieldHolopainen0020     [100];
        Double_t errYieldXHolopainen0020    [100];
        Int_t nlinesHolopainen0020                          = 0;
        
        TString fileNameHolopainen0020                      = "ExternalInputPbPb/Theory/Holopainen/yieldAndv2_0020_toRead.txt";
        ifstream  fileHolopainen0020;
        fileHolopainen0020.open(fileNameHolopainen0020,ios_base::in);
        cout << fileNameHolopainen0020 << endl;
        
        while(!fileHolopainen0020.eof() && nlinesHolopainen0020< 100){
            fileHolopainen0020 >> ptHolopainen0020[nlinesHolopainen0020] >> yieldHolopainen0020[nlinesHolopainen0020] >> v2Holopainen0020[nlinesHolopainen0020] ; 
            cout << nlinesHolopainen0020 << "\t"  << ptHolopainen0020[nlinesHolopainen0020] << "\t"  << v2Holopainen0020[nlinesHolopainen0020] << "\t"  << yieldHolopainen0020[nlinesHolopainen0020] << endl;    
            errYieldHolopainen0020[nlinesHolopainen0020]    = 0;
            errYieldXHolopainen0020[nlinesHolopainen0020]   = 0;
            nlinesHolopainen0020++;
        }
        fileHolopainen0020.close();
        TGraphErrors* graphDirectPhotonHolopainen0020       = new TGraphErrors(nlinesHolopainen0020-1,ptHolopainen0020,yieldHolopainen0020, errYieldXHolopainen0020, errYieldHolopainen0020); 
        while (graphDirectPhotonHolopainen0020->GetX()[0] < 0.7) graphDirectPhotonHolopainen0020->RemovePoint(0);
        TGraphErrors* graphDirectPhotonV2Holopainen0020     = new TGraphErrors(nlinesHolopainen0020-1,ptHolopainen0020,v2Holopainen0020, errYieldXHolopainen0020, errYieldHolopainen0020); 
        while (graphDirectPhotonV2Holopainen0020->GetX()[0] < 0.7) graphDirectPhotonV2Holopainen0020->RemovePoint(0);
        
        Double_t ptHolopainen2040           [100];
        Double_t v2Holopainen2040           [100];
        Double_t yieldHolopainen2040        [100];
        Double_t errYieldHolopainen2040     [100];
        Double_t errYieldXHolopainen2040    [100];
        Int_t nlinesHolopainen2040                          = 0;
        
        TString fileNameHolopainen2040 = "ExternalInputPbPb/Theory/Holopainen/yieldAndv2_2040_toRead.txt";
        ifstream  fileHolopainen2040;
        fileHolopainen2040.open(fileNameHolopainen2040,ios_base::in);
        cout << fileNameHolopainen2040 << endl;
        
        while(!fileHolopainen2040.eof() && nlinesHolopainen2040< 100){
            fileHolopainen2040 >> ptHolopainen2040[nlinesHolopainen2040] >> yieldHolopainen2040[nlinesHolopainen2040] >> v2Holopainen2040[nlinesHolopainen2040] ; 
            cout << nlinesHolopainen2040 << "\t"  << ptHolopainen2040[nlinesHolopainen2040] << "\t"  << v2Holopainen2040[nlinesHolopainen2040] << "\t"  << yieldHolopainen2040[nlinesHolopainen2040] << endl;    
            errYieldHolopainen2040[nlinesHolopainen2040]    = 0;
            errYieldXHolopainen2040[nlinesHolopainen2040]   = 0;
            nlinesHolopainen2040++;
        }
        fileHolopainen2040.close();
        TGraphErrors* graphDirectPhotonHolopainen2040       = new TGraphErrors(nlinesHolopainen2040-1,ptHolopainen2040,yieldHolopainen2040, errYieldXHolopainen2040, errYieldHolopainen2040); 
        while (graphDirectPhotonHolopainen2040->GetX()[0] < 0.7) graphDirectPhotonHolopainen2040->RemovePoint(0);
        TGraphErrors* graphDirectPhotonV2Holopainen2040     = new TGraphErrors(nlinesHolopainen2040-1,ptHolopainen2040,v2Holopainen2040, errYieldXHolopainen2040, errYieldHolopainen2040); 
        while (graphDirectPhotonV2Holopainen2040->GetX()[0] < 0.7) graphDirectPhotonV2Holopainen2040->RemovePoint(0);
        
        
        //******************************************************************************************************************
        //************************************** Writing output for PbPb ***************************************************
        //******************************************************************************************************************
        TFile *fileTheoryGraphsPbPb = new TFile("ExternalInputPbPb/Theory/TheoryCompilationPbPb.root","RECREATE");

        
            fileTheoryGraphsPbPb->mkdir("DirectPhoton");
            TDirectoryFile* directoryGamma = (TDirectoryFile*)fileTheoryGraphsPbPb->Get("DirectPhoton"); 
            fileTheoryGraphsPbPb->cd("DirectPhoton");

            graphDirectPhotonMCGill0020->Write("graphDirectPhotonYield_McGill_0020", TObject::kOverwrite);
            graphDirectPhotonMCGill2040->Write("graphDirectPhotonYield_McGill_2040", TObject::kOverwrite);
            graphPromptPhotonMCGill0020->Write("graphPromptPhotonYield_McGill_0020", TObject::kOverwrite);
            graphPromptPhotonMCGill2040->Write("graphPromptPhotonYield_McGill_2040", TObject::kOverwrite);
            
            graphDirectPhotonPHSD0020->Write("graphDirectPhotonYield_PHSD_0020", TObject::kOverwrite);
            graphDirectPhotonPHSD0040->Write("graphDirectPhotonYield_PHSD_0040", TObject::kOverwrite);
            graphDirectPhotonPHSD2040->Write("graphDirectPhotonYield_PHSD_2040", TObject::kOverwrite);
            graphDirectPhotonPHSD4080->Write("graphDirectPhotonYield_PHSD_4080", TObject::kOverwrite);

            graphDirectPhotonHe0020->Write("graphDirectPhotonYield_He_0020", TObject::kOverwrite);
    //         graphDirectPhotonHe0040->Write("graphDirectPhotonYield_He_0040", TObject::kOverwrite);
            graphDirectPhotonHe2040->Write("graphDirectPhotonYield_He_2040", TObject::kOverwrite);
            graphDirectPhotonHe4080->Write("graphDirectPhotonYield_He_4080", TObject::kOverwrite);
            
            graphDirectPhotonChatterjee0020->Write("graphDirectPhotonYield_Chatterjee_0020", TObject::kOverwrite);
            graphDirectPhotonThermalChatterjee0020->Write("graphDirectPhotonThermalYield_Chatterjee_0020", TObject::kOverwrite);
            graphDirectPhotonPromptChatterjee0020->Write("graphDirectPhotonPromptYield_Chatterjee_0020", TObject::kOverwrite);
            graphDirectPhotonChatterjee2040->Write("graphDirectPhotonYield_Chatterjee_2040", TObject::kOverwrite);
            graphDirectPhotonThermalChatterjee2040->Write("graphDirectPhotonThermalYield_Chatterjee_2040", TObject::kOverwrite);
            graphDirectPhotonPromptChatterjee2040->Write("graphDirectPhotonPromptYield_Chatterjee_2040", TObject::kOverwrite);
            graphDirectPhotonChatterjee4060->Write("graphDirectPhotonYield_Chatterjee_4060", TObject::kOverwrite);
            graphDirectPhotonThermalChatterjee4060->Write("graphDirectPhotonThermalYield_Chatterjee_4060", TObject::kOverwrite);
            graphDirectPhotonPromptChatterjee4060->Write("graphDirectPhotonPromptYield_Chatterjee_4060", TObject::kOverwrite);

            graphDirectPhotonChatterjee0020_2->Write("graphDirectPhotonYield_Chatterjee_0020_2", TObject::kOverwrite);
            graphDirectPhotonThermalChatterjee0020_2->Write("graphDirectPhotonThermalYield_Chatterjee_0020_2", TObject::kOverwrite);
            graphDirectPhotonPromptChatterjee0020_2->Write("graphDirectPhotonPromptYield_Chatterjee_0020_2", TObject::kOverwrite);
            graphDirectPhotonChatterjee2040_2->Write("graphDirectPhotonYield_Chatterjee_2040_2", TObject::kOverwrite);
            graphDirectPhotonThermalChatterjee2040_2->Write("graphDirectPhotonThermalYield_Chatterjee_2040_2", TObject::kOverwrite);
            graphDirectPhotonPromptChatterjee2040_2->Write("graphDirectPhotonPromptYield_Chatterjee_2040_2", TObject::kOverwrite);
            
            graphDirectPhotonHees0020->Write("graphDirectPhotonYield_VanHees_0020", TObject::kOverwrite);
            graphDirectPhotonQGPHees0020->Write("graphDirectPhotonQGPYield_VanHees_0020", TObject::kOverwrite);
            graphDirectPhotonPrimordialHees0020->Write("graphDirectPhotonPrimordialYield_VanHees_0020", TObject::kOverwrite);
            graphDirectPhotonOmegaHees0020->Write("graphDirectPhotonOmegaYield_VanHees_0020", TObject::kOverwrite);
            graphDirectPhotonMesonGasHees0020->Write("graphDirectPhotonMesonGasYield_VanHees_0020", TObject::kOverwrite);
            graphDirectPhotonRhoSFHees0020->Write("graphDirectPhotonRhoSFYield_VanHees_0020", TObject::kOverwrite);

            graphDirectPhotonHees2040->Write("graphDirectPhotonYield_VanHees_2040", TObject::kOverwrite);
            graphDirectPhotonQGPHees2040->Write("graphDirectPhotonQGPYield_VanHees_2040", TObject::kOverwrite);
            graphDirectPhotonPrimordialHees2040->Write("graphDirectPhotonPrimordialYield_VanHees_2040", TObject::kOverwrite);
            graphDirectPhotonOmegaHees2040->Write("graphDirectPhotonOmegaYield_VanHees_2040", TObject::kOverwrite);
            graphDirectPhotonMesonGasHees2040->Write("graphDirectPhotonMesonGasYield_VanHees_2040", TObject::kOverwrite);
            graphDirectPhotonRhoSFHees2040->Write("graphDirectPhotonRhoSFYield_VanHees_2040", TObject::kOverwrite);

            graphDirectPhotonHees4080->Write("graphDirectPhotonYield_VanHees_4080", TObject::kOverwrite);
            graphDirectPhotonQGPHees4080->Write("graphDirectPhotonQGPYield_VanHees_4080", TObject::kOverwrite);
            graphDirectPhotonPrimordialHees4080->Write("graphDirectPhotonPrimordialYield_VanHees_4080", TObject::kOverwrite);
            graphDirectPhotonOmegaHees4080->Write("graphDirectPhotonOmegaYield_VanHees_4080", TObject::kOverwrite);
            graphDirectPhotonMesonGasHees4080->Write("graphDirectPhotonMesonGasYield_VanHees_4080", TObject::kOverwrite);
            graphDirectPhotonRhoSFHees4080->Write("graphDirectPhotonRhoSFYield_VanHees_4080", TObject::kOverwrite);
            
            graphDirectPhotonHolopainen0020->Write("graphDirectPhotonYield_Holopainen_0020", TObject::kOverwrite);
            graphDirectPhotonV2Holopainen0020->Write("graphDirectPhotonV2_Holopainen_0020", TObject::kOverwrite);
            graphDirectPhotonHolopainen2040->Write("graphDirectPhotonYield_Holopainen_2040", TObject::kOverwrite);
            graphDirectPhotonV2Holopainen2040->Write("graphDirectPhotonV2_Holopainen_2040", TObject::kOverwrite);
            
        fileTheoryGraphsPbPb->Close();
    }
        
    
    
}