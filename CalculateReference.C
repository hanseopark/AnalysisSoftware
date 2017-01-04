/****************************************************************************************************************************
******          provided by Gamma Conversion Group, PWG-GA                                                          *********
******          Daniel Mühlheim, d.muehlheim@cern.ch                                                                *********
******          Friederike Bock, friederike.bock@cern.ch                                                            *********
*****************************************************************************************************************************/

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
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ExtractSignalBinning.h"
// #include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"
#include "TFitResultPtr.h"

TGraphAsymmErrors *GetInterpolSpectrum2D(Int_t nDataPoints, TGraphAsymmErrors** graphs, Double_t* energies, Double_t dSqrts);
TH1D* ConvertYieldHisto(TH1D* input, Bool_t DivideBy2pi, Bool_t DivideByPt, Bool_t MultiplyBy2pi, Bool_t MultiplyByPt);
TF1* DoFitWithTsallis(TGraph* graph, TString name, TString particle, Double_t p0, Double_t p1, Double_t p2);
TF1* DoFitWithTCM(TGraph* graph, TString name, TString particle, Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4 );
void PlotInterpolationPtBins(TGraphErrors** gPtvSqrts,TGraphErrors** gPtvsEnergies, TF1** fPowerlaw, TGraphAsymmErrors* gRpPb,Int_t fColumnPlot, Int_t fRowPlot,TString namePlot);
void PlotAlphavsPt(TGraphErrors* gAlpha, TString method, TString thesisPlotLabel, TString namePlot);
void PlotWithFit(TCanvas* canvas, TH2F* hist, TH2F* hist2, TGraphAsymmErrors* graph, TGraphAsymmErrors* graphRebin, TF1* fit, TString name, TString outputDir, TString suffix, TString energyCur);
void PlotWithFit(TCanvas *canvas, TH2F* hist, TH2F* hist2, TGraphErrors* graph, TGraphErrors* graphRebin, TF1* fit, TString name, TString outputDir, TString suffix, TString energyCur);
void PlotGraphsOfAllEnergies( TCanvas* canvas, TH2F* hist, Int_t nDataSets, TGraphAsymmErrors** graphs, TF1** fits, TString* energies, TString name, TString outputDir, TString suffix );

void SetStyleTGraphErrorsForGraphs(TGraph* graph, TString XTitle, TString YTitle, Size_t xLableSize, Size_t xTitleSize, Size_t yLableSize, Size_t yTitleSize, Float_t xTitleOffset = 1, Float_t yTitleOffset = 1, Int_t xNDivisions = 510, Int_t yNDivisions = 510);
TGraphAsymmErrors* SetGraphErrors (TGraphAsymmErrors* graph, TGraphAsymmErrors* refGraph);
TGraphAsymmErrors* RemoveSysErrors (TGraphAsymmErrors* graph, TGraphAsymmErrors* refGraph);

void RemoveZerosAndPrint(TGraphErrors* graph, TString d){
  cout << d.Data() << endl;
  while (graph->GetY()[0] == 0. ) graph->RemovePoint(0);
  graph->Print();
}
void RemoveZerosAndPrint(TGraphAsymmErrors* graph, TString d){
  cout << d.Data() << endl;
  while (graph->GetY()[0] == 0. ) graph->RemovePoint(0);
  graph->Print();
}

TGraphErrors* graphAlpha            = 0x0;
TGraphErrors** graphPtvsSqrts       = 0x0;
TGraphErrors** gPtvsEnergiesSystem  = 0x0;
TF1** fPowerlawSystem               = 0x0;

void CalculateStatPlusSysErrors(TH1D* histStat, TH1D* histSys);
TGraphAsymmErrors* CalculateStatPlusSysErrors(TGraphAsymmErrors* graphStat, TGraphAsymmErrors* graphSys);

//________________________________________________________________________________________________________________________
void CalculateReference (   TString configFile      = "",
                            Int_t nDataSets         = 2,
                            TString meson           = "Pi0",
                            TString suffix          = "eps",
                            TString finalEnergy     = "8TeV",
                            Double_t energy         = 8000,
                            Int_t mode              = 20,
                            Bool_t doSpecialBinning = kFALSE,
                            TString binningEnergy   = ""
                        ){

    //*************************************************************************************************
    //*************************** general settings
    //*************************************************************************************************
    TH1::AddDirectory(kFALSE);
    gROOT->Reset();
    gROOT->SetStyle("Plain");

    StyleSettingsThesis(suffix);
    SetPlotStyle();

    TString dateForOutput   = ReturnDateStringForOutput();
    TString outputDir       = Form("%s/%s/CalculateReference",suffix.Data(),dateForOutput.Data());
    gSystem->Exec("mkdir -p "+outputDir);
    TString modeName        = ReturnTextReconstructionProcess(mode);
    
    //*************************************************************************************************
    //***************************** read configurarion file *******************************************
    //*************************************************************************************************
    TString nameFile[5]             = {"", "", "", "", ""};
    TString nameHist[4][5]          = {{"", "", "", "", ""},{"", "", "", "", ""},{"", "", "", "", ""},{"", "", "", "", ""}};
    Int_t isHist[2][5]              = {{0, 0, 0, 0, 0},{0, 0, 0, 0, 0}};
    TString energyIndName[6]        = {"", "", "", "", "", ""};
    Double_t energyInd[6]           = {0, 0, 0, 0, 0, 0};
    Double_t scaleFactor[5]         = {1, 1, 1, 1, 1};
    Int_t nrOfTrigToBeComb          = 0;
    // number of triggers which are really used for the respective analysis
    Int_t nEnergyRead               = 0;
    ifstream in(configFile.Data());
    while(!in.eof() && nEnergyRead<nDataSets ){
        cout << "read line:" <<  nEnergyRead+1 << endl;
        in >> nameFile[nEnergyRead] >> isHist[0][nEnergyRead] >> nameHist[0][nEnergyRead] >> isHist[1][nEnergyRead] >> nameHist[1][nEnergyRead]
        >> energyIndName[nEnergyRead]  >> energyInd[nEnergyRead]  >> scaleFactor[nEnergyRead] >> nameHist[2][nEnergyRead] >> nameHist[3][nEnergyRead];;
        cout << nameFile[nEnergyRead] << "\t" << isHist[0][nEnergyRead] << "\t" << nameHist[0][nEnergyRead].Data() << "\t" << isHist[1][nEnergyRead] << "\t" << nameHist[1][nEnergyRead].Data()
        << "\t" << energyIndName[nEnergyRead].Data()  << "\t" << energyInd[nEnergyRead]  << "\t" << scaleFactor[nEnergyRead] << endl;
        nEnergyRead++;
    }

    Double_t finalBinningPt[100];
    Int_t maxNBins                  = 0;
    if (doSpecialBinning){
        maxNBins                    = GetBinning( finalBinningPt, meson.Data(), binningEnergy.Data(), mode );
        if (maxNBins == 0){
            cout << "The requested binning doesn't exist, aborting!" << endl;
            return;
        }
        cout << "Binning" << endl;        
        for (Int_t i = 0; i<maxNBins; i++){
            cout << i << "\t"<< finalBinningPt[i] << "-" << finalBinningPt[i] <<endl;
        }
    }
    
    //*************************************************************************************************
    //*************************** read input
    //*************************************************************************************************
    TFile* inputFile[5]                 = {NULL, NULL, NULL, NULL, NULL};     
    TH1D* histStat[5]                   = {NULL, NULL, NULL, NULL, NULL};     
    TH1D* histSyst[5]                   = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histComb[5]                   = {NULL, NULL, NULL, NULL, NULL}; 
    TH1D* histRelStat[5]                = {NULL, NULL, NULL, NULL, NULL};     
    TH1D* histRelSyst[5]                = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histRelComb[5]                = {NULL, NULL, NULL, NULL, NULL}; 
    TGraphAsymmErrors* graphStat[5]     = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphSyst[5]     = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphComb[5]     = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombReb[6]  = {NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphSystReb[6]  = {NULL, NULL, NULL, NULL, NULL, NULL};
    TF1* fitComb[6]                     = {NULL, NULL, NULL, NULL, NULL, NULL};
    Double_t minY                       = 1e100;
    Double_t maxY                       = 0;
    for (Int_t i = 0; i < nDataSets; i++){
        inputFile[i]                    = new TFile(nameFile[i].Data());
        if (isHist[0][i] && isHist[1][i]){
            histStat[i]                 = (TH1D*)inputFile[i]->Get(nameHist[0][i].Data());
            histStat[i]->Scale(scaleFactor[i]);
            histSyst[i]                 = (TH1D*)inputFile[i]->Get(nameHist[1][i].Data());
            histSyst[i]->Scale(scaleFactor[i]);
            histComb[i]                 = (TH1D*)histStat[i]->Clone(Form("histComb_%d",i));
            CalculateStatPlusSysErrors(histComb[i],histSyst[i]);
            graphComb[i]                = new TGraphAsymmErrors(histComb[i]);
            graphStat[i]                = new TGraphAsymmErrors(histStat[i]);
            graphSyst[i]                = new TGraphAsymmErrors(histSyst[i]);
            histRelStat[i]              = (TH1D*)histStat[i]->Clone(Form("Rel%s",nameHist[0][i].Data()));
            histRelSyst[i]              = (TH1D*)histSyst[i]->Clone(Form("Rel%s",nameHist[1][i].Data()));
            histRelComb[i]              = (TH1D*)histComb[i]->Clone(Form("Rel_histComb_%d",i));
            for (Int_t j = 1; j< histRelStat[i]->GetNbinsX()+1;j++){
                histRelStat[i]->SetBinContent(j,histRelStat[i]->GetBinError(j)/histRelStat[i]->GetBinContent(j));
                histRelStat[i]->SetBinError(j, 0);
                histRelSyst[i]->SetBinContent(j, histRelSyst[i]->GetBinError(j)/histRelSyst[i]->GetBinContent(j));
                histRelSyst[i]->SetBinError(j, 0);
                histRelComb[i]->SetBinContent(j, histRelComb[i]->GetBinError(j)/histRelComb[i]->GetBinContent(j));
                histRelComb[i]->SetBinError(j, 0);
            }    
        } else if (!isHist[0][i] && !isHist[1][i]){
            graphStat[i]                = (TGraphAsymmErrors*)inputFile[i]->Get(nameHist[0][i].Data());
            ScaleGraph(graphStat[i],scaleFactor[i]);
            graphSyst[i]                = (TGraphAsymmErrors*)inputFile[i]->Get(nameHist[1][i].Data());
            ScaleGraph(graphSyst[i],scaleFactor[i]);
            if (nameHist[2][i].CompareTo("bla")!= 0){
                graphComb[i]            = (TGraphAsymmErrors*)inputFile[i]->Get(nameHist[2][i].Data());
                ScaleGraph(graphComb[i],scaleFactor[i]);
            } else {
                graphComb[i]            = CalculateStatPlusSysErrors(graphStat[i],graphSyst[i]);
            }
            // read binning from graphs
            Double_t binningGraphs[100];
            binningGraphs[0]            = 0;
            for (Int_t j = 0; j < graphStat[i]->GetN(); j++){
                binningGraphs[j+1]      = graphStat[i]->GetX()[j]-graphStat[i]->GetEXlow()[j];
                cout << binningGraphs[j+1] << ",";
                if (j == graphStat[i]->GetN()-1){
                    binningGraphs[j+2]      = graphStat[i]->GetX()[j]+graphStat[i]->GetEXhigh()[j];
                    cout << binningGraphs[j+2];
                }    
            }    
            cout << endl;

            histRelStat[i]              = new TH1D(Form("Rel%s",nameHist[0][i].Data()),Form("Rel%s",nameHist[0][i].Data()),graphStat[i]->GetN()+1, binningGraphs );
            histRelComb[i]              = new TH1D(Form("Rel_histComb_%d",i), Form("Rel_histComb_%d",i), graphStat[i]->GetN()+1, binningGraphs );
            histRelSyst[i]              = new TH1D(Form("Rel%s",nameHist[0][i].Data()),Form("Rel%s",nameHist[0][i].Data()),graphStat[i]->GetN()+1, binningGraphs );
            
            for (Int_t j = 0; j< graphStat[i]->GetN()+1; j++){
                cout << histRelStat[i]->GetBinCenter(j) << "\t" << histRelStat[i]->GetBinWidth(j) << endl;
                if (graphStat[i]->GetY()[j] != 0){
                    histRelStat[i]->SetBinContent(j+2, (graphStat[i]->GetEYhigh()[j]+graphStat[i]->GetEYlow()[j])/(2*graphStat[i]->GetY()[j]));
                    histRelStat[i]->SetBinError(j+2, 0.);
                    histRelSyst[i]->SetBinContent(j+2, (graphSyst[i]->GetEYhigh()[j]+graphSyst[i]->GetEYlow()[j])/(2*graphSyst[i]->GetY()[j]));
                    histRelSyst[i]->SetBinError(j+2, 0.);
                    histRelComb[i]->SetBinContent(j+2, (graphComb[i]->GetEYhigh()[j]+graphComb[i]->GetEYlow()[j])/(2*graphComb[i]->GetY()[j]));
                    histRelComb[i]->SetBinError(j+2, 0.);
                } else {
                    histRelStat[i]->SetBinContent(j+2, 0);
                    histRelStat[i]->SetBinError(j+2, 0.);
                    histRelSyst[i]->SetBinContent(j+2, 0);
                    histRelSyst[i]->SetBinError(j+2, 0.);
                    histRelComb[i]->SetBinContent(j+2, 0);
                    histRelComb[i]->SetBinError(j+2, 0.);    
                }    
            }
        }
        RemoveZerosAndPrint(graphComb[i],Form("graphComb_%d",i));
        RemoveZerosAndPrint(graphStat[i],Form("graphStat_%d",i));
        RemoveZerosAndPrint(graphSyst[i],Form("graphSyst_%d",i));
        
        if (maxY < graphComb[i]->GetY()[0])
            maxY        = graphComb[i]->GetY()[0];
        if (minY > graphComb[i]->GetY()[graphComb[i]->GetN()-1])
            minY        = graphComb[i]->GetY()[graphComb[i]->GetN()-1];
        //*************************************************************************************************
        //*************************** Fits
        //*************************************************************************************************
        if (nameHist[3][i].CompareTo("bla") != 0){
            fitComb[i]                  = (TF1*)inputFile[i]->Get(nameHist[3][i].Data());
        } else {
            fitComb[i]                  = DoFitWithTsallis(graphComb[i],Form("fitComb_%d",i),meson.Data(), graphComb[i]->GetY()[0],7.,0.2);
//             fitComb[i]                  = DoFitWithTCM(graphComb[i],Form("fitComb_%d",i),meson.Data(), graphComb[i]->GetY()[0],7.,0.2,graphComb[i]->GetY()[0]/10,0.3);
        }    
        
        if (doSpecialBinning){
            graphCombReb[i]             = new TGraphAsymmErrors(maxNBins);
            graphSystReb[i]             = new TGraphAsymmErrors(maxNBins);
            for(Int_t j = 0; j<maxNBins; j++){
                Double_t ptCurrent      = (finalBinningPt[j+1]+finalBinningPt[j])/2.;
                
                graphCombReb[i]->SetPoint(j,ptCurrent,fitComb[i]->Eval(ptCurrent));
                Double_t relCombErr     =  histRelComb[i]->Interpolate(ptCurrent);
                Double_t relSystErr     =  histRelSyst[i]->Interpolate(ptCurrent);
                if (relCombErr != 0){
                    graphCombReb[i]->SetPointError(j,(finalBinningPt[j+1]-finalBinningPt[j])/2.,(finalBinningPt[j+1]-finalBinningPt[j])/2.,relCombErr*fitComb[i]->Eval(ptCurrent),relCombErr*fitComb[i]->Eval(ptCurrent));
                } else {
                    graphCombReb[i]->SetPointError(j,(finalBinningPt[j+1]-finalBinningPt[j])/2.,(finalBinningPt[j+1]-finalBinningPt[j])/2.,fitComb[i]->Eval(ptCurrent)*0.2,fitComb[i]->Eval(ptCurrent)*0.2);
                }    
                graphSystReb[i]->SetPoint(j,graphSystReb[i]->GetX()[j],fitComb[i]->Eval(graphSystReb[i]->GetX()[j]));
                if (relSystErr != 0){
                    graphSystReb[i]->SetPointError(j,(finalBinningPt[j+1]-finalBinningPt[j])/2.,(finalBinningPt[j+1]-finalBinningPt[j])/2.,relSystErr*fitComb[i]->Eval(ptCurrent),relSystErr*fitComb[i]->Eval(ptCurrent));
                } else {
                    graphSystReb[i]->SetPointError(j,(finalBinningPt[j+1]-finalBinningPt[j])/2.,(finalBinningPt[j+1]-finalBinningPt[j])/2.,fitComb[i]->Eval(ptCurrent)*0.2,fitComb[i]->Eval(ptCurrent)*0.2);
                }    
            }
        } else {    
            if (i==0 ){
                graphCombReb[i]             = (TGraphAsymmErrors*)graphComb[i]->Clone(Form("graphCombReb_%d",i));
                graphSystReb[i]             = (TGraphAsymmErrors*)graphSyst[i]->Clone(Form("graphSystReb_%d",i));
            } else {
                graphCombReb[i]             = (TGraphAsymmErrors*)graphComb[0]->Clone(Form("graphCombReb_%d",i));
                graphSystReb[i]             = (TGraphAsymmErrors*)graphSyst[0]->Clone(Form("graphSystReb_%d",i));
                for(Int_t j = 0; j<graphCombReb[i]->GetN(); j++){
                    graphCombReb[i]->SetPoint(j,graphCombReb[i]->GetX()[j],fitComb[i]->Eval(graphCombReb[i]->GetX()[j]));
                    Double_t relCombErr     =  histRelComb[i]->Interpolate(graphCombReb[i]->GetX()[j]);
                    Double_t relSystErr     =  histRelSyst[i]->Interpolate(graphCombReb[i]->GetX()[j]);
                    if (relCombErr != 0){
                        graphCombReb[i]->SetPointError(j,graphCombReb[i]->GetEXlow()[j],graphCombReb[i]->GetEXhigh()[j],relCombErr*fitComb[i]->Eval(graphCombReb[i]->GetX()[j]),relCombErr*fitComb[i]->Eval(graphCombReb[i]->GetX()[j]));
                    } else {
                        graphCombReb[i]->SetPointError(j,graphCombReb[i]->GetEXlow()[j],graphCombReb[i]->GetEXhigh()[j],fitComb[i]->Eval(graphCombReb[i]->GetX()[j])*0.2,fitComb[i]->Eval(graphCombReb[i]->GetX()[j])*0.2);
                    }    
                    graphSystReb[i]->SetPoint(j,graphSystReb[i]->GetX()[j],fitComb[i]->Eval(graphSystReb[i]->GetX()[j]));
                    if (relSystErr != 0){
                        graphSystReb[i]->SetPointError(j,graphSystReb[i]->GetEXlow()[j],graphSystReb[i]->GetEXhigh()[j],relSystErr*fitComb[i]->Eval(graphCombReb[i]->GetX()[j]),relSystErr*fitComb[i]->Eval(graphCombReb[i]->GetX()[j]));
                    } else {
                        graphSystReb[i]->SetPointError(j,graphSystReb[i]->GetEXlow()[j],graphSystReb[i]->GetEXhigh()[j],fitComb[i]->Eval(graphSystReb[i]->GetX()[j])*0.2,fitComb[i]->Eval(graphSystReb[i]->GetX()[j])*0.2);
                    }    
                }
            }
        }    
    }

// 
//     TFile* input8TeV              = new TFile("ExternalInput/IdentifiedCharged/spectra8TeV_MarekFirstShot20161129.root");
//     TList *list                   = (TList*)input8TeV->Get("output");
//     TH1D* hNegPionResult8TeV      = (TH1D*)list->FindObject("SpectraFinalPionMinus");
//     TH1D* hPosPionResult8TeV      = (TH1D*)list->FindObject("SpectraFinalPionPlus");
//     hPosPionResult8TeV->Scale(0.5);
//     hPosPionResult8TeV->Add(hNegPionResult8TeV,0.5);
//     ConvertYieldHisto(hPosPionResult8TeV,kTRUE,kTRUE,kFALSE,kFALSE);
//     TGraphErrors* graphChargedPion8TeVStat    = new TGraphErrors(hPosPionResult8TeV);
//     RemoveZerosAndPrint(graphChargedPion8TeVStat,"graphChargedPion8TeVStat");
// 

    //*************************************************************************************************
    //*************************** extrapolate spectra *************************************************
    //*************************************************************************************************
    TGraphAsymmErrors*  graphFinalEnergy    = GetInterpolSpectrum2D(    nDataSets,
                                                                        graphCombReb,
                                                                        energyInd,
                                                                        energy
                                                                    );

    if(graphAlpha && graphPtvsSqrts && gPtvsEnergiesSystem && fPowerlawSystem ){
        Int_t columns   = 2;
        Int_t rows      = 2;
        Int_t counter   = 0;
        while (columns*rows < graphCombReb[0]->GetN()){
            if (counter%2 != 0) rows++;
            else columns++;
            counter++;
        }    
        PlotInterpolationPtBins(graphPtvsSqrts,gPtvsEnergiesSystem,fPowerlawSystem,graphFinalEnergy,columns, rows,Form("%s/%s_%s_Pt_vs_Sqrts.%s",outputDir.Data(),meson.Data(),modeName.Data(), suffix.Data()));
        PlotAlphavsPt(graphAlpha, "pp", meson.Data(), Form("%s/%s_%s_Alpha_vs_Pt.%s", outputDir.Data(),meson.Data(),modeName.Data(), suffix.Data()));
    }else{
        cout << "ERROR: NULL pointer - returning..." << endl;
        cout << graphAlpha << endl;
        cout << graphPtvsSqrts << endl;
        cout << gPtvsEnergiesSystem << endl;
        cout << fPowerlawSystem << endl;
        return;
    }
    graphFinalEnergy->Print();
    
    TF1* fitFinal                       = DoFitWithTsallis(graphFinalEnergy,Form("fitComb_%s",finalEnergy.Data()),meson.Data(), graphFinalEnergy->GetY()[0],7,0.2); 
//     TF1* fitFinal                       = DoFitWithTCM(graphFinalEnergy,Form("fitComb_%s",finalEnergy.Data()),meson.Data(), graphFinalEnergy->GetY()[0],7.,0.2,graphFinalEnergy->GetY()[0]/10,0.3);

    graphCombReb[nDataSets]             = graphFinalEnergy;
    fitComb[nDataSets]                  = fitFinal;
    energyIndName[nDataSets]            = finalEnergy;
    energyInd[nDataSets]                = energy;
    if (maxY < graphFinalEnergy->GetY()[0])
        maxY        = graphFinalEnergy->GetY()[0];
    if (minY > graphFinalEnergy->GetY()[graphFinalEnergy->GetN()-1])
        minY        = graphFinalEnergy->GetY()[graphFinalEnergy->GetN()-1];
    
    //*************************************************************************************************
    //*************************** plotting
    //*************************************************************************************************
    TString axisLabel   = "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (c/GeV)^{2}";
    if (maxY > 1000)
        axisLabel       = "#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )";
    
    TCanvas* canvasDummy2 = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
    DrawGammaCanvasSettings( canvasDummy2,  0.15, 0.01, 0.015, 0.08);
    canvasDummy2->SetLogy();
    canvasDummy2->SetLogx();
    TH2F * histo2DDummy = new TH2F("histo2DDummy","histo2DDummy",1000,0.1,31.,1000,minY/10,maxY*10);
    SetStyleHistoTH2ForGraphs(histo2DDummy, "#it{p}_{T} (GeV/#it{c})",axisLabel, 0.032,0.04, 0.032,0.04, 0.85,1.65);
    histo2DDummy->GetXaxis()->SetMoreLogLabels();
    histo2DDummy->GetXaxis()->SetLabelOffset(-0.008);
    TH2F * histo2DDummy2 = new TH2F("histo2DDummy","histo2DDummy",1000,0.1,31.,1000, 0, 2);
    SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","data/fit", 0.032,0.04, 0.032,0.04, 0.85,1.65);
    histo2DDummy2->GetXaxis()->SetMoreLogLabels();
    histo2DDummy2->GetXaxis()->SetLabelOffset(-0.008);
    
    for (Int_t i = 0; i < nDataSets; i++){
        PlotWithFit(canvasDummy2, histo2DDummy, histo2DDummy2, graphComb[i], graphCombReb[i], fitComb[i], Form("input_%s%s%s_withFit",meson.Data(), modeName.Data(), energyIndName[i].Data()), outputDir, suffix, energyIndName[i]);    
    }    
    PlotWithFit(canvasDummy2, histo2DDummy, histo2DDummy2, graphFinalEnergy, 0x0, fitFinal, Form("compare_%s%s%s_withFit", meson.Data(), modeName.Data(), finalEnergy.Data()), outputDir, suffix, energyIndName[nDataSets]);
    PlotGraphsOfAllEnergies(canvasDummy2, histo2DDummy, nDataSets+1, graphCombReb, fitComb, energyIndName, Form("output%s%s_withFit",meson.Data(), modeName.Data()),  outputDir, suffix);
    

    
    // **************************************************************************************************
    // ************************ prepare histos and graphs for writing ***********************************
    // **************************************************************************************************
    TH1D* histoFinalStat                    = NULL;
    if (maxNBins > 0){
        histoFinalStat                    = new TH1D(Form("histStatErr%s_%s_%s",modeName.Data(),meson.Data(),finalEnergy.Data()),Form("histStatErr_%s_%s",meson.Data(),finalEnergy.Data()),maxNBins,finalBinningPt);
        for(Int_t i=1; i<maxNBins+1; i++){
            histoFinalStat->SetBinContent(i,fitFinal->Eval(histoFinalStat->GetBinCenter(i)));
            histoFinalStat->SetBinError(i,0.1*fitFinal->Eval(histoFinalStat->GetBinCenter(i)));
        }
    }
    
    TGraphAsymmErrors* graphStatFinalEnergy = RemoveSysErrors(graphFinalEnergy,graphSystReb[0]);    
    graphFinalEnergy = SetGraphErrors(graphFinalEnergy,graphSystReb[0]);

    //*************************************************************************************************
    //*************************** write output file ***************************************************
    //*************************************************************************************************
    TFile *fOutput = new TFile(Form("%s/Interpolation.root",outputDir.Data()),"UPDATE");

        graphStatFinalEnergy->Write(Form("graphStatErr%s_%s_%s", modeName.Data(),meson.Data(), finalEnergy.Data()),TObject::kOverwrite);
        graphFinalEnergy->Write(Form("graphSystErr%s_%s_%s",modeName.Data(),meson.Data(),finalEnergy.Data()),TObject::kOverwrite);
        if (histoFinalStat) histoFinalStat->Write(Form("histStatErr%s_%s_%s",modeName.Data(),meson.Data(),finalEnergy.Data()),TObject::kOverwrite);

    fOutput->Write();
    fOutput->Close();

    return;
}


//________________________________________________________________________________________________________________________
TGraphAsymmErrors *GetInterpolSpectrum2D(Int_t nDataPoints, TGraphAsymmErrors** graphs, Double_t* energies, Double_t dSqrts)
{
    if (!graphs){   
        cout << "failed to load array" << endl;
    }    
    for (Int_t j = 0; j< nDataPoints; j++){
        if(!graphs[j]) {
            cout << "couldn't load graph " << j << " properly." << endl;
            return 0x0;
        }    
    }
    // initialize graphs properly
    TGraphAsymmErrors  *gInterpol   = new TGraphAsymmErrors(graphs[0]->GetN());
    TGraphErrors  *gAlpha           = new TGraphErrors(graphs[0]->GetN());
    TGraphErrors** gPtvsSqrts       = new TGraphErrors*[graphs[0]->GetN()];
    TGraphErrors** gPtvsEnergies    = new TGraphErrors*[graphs[0]->GetN()];
    TF1**          fPowerlawFits    = new TF1*[graphs[0]->GetN()];

    for(Int_t i = 0; i < graphs[0]->GetN(); i++){
        TGraphErrors *grint = new TGraphErrors(1);
        grint->SetPoint(0, dSqrts, 0);

        // check if graph are in same binning?
        Bool_t isSameBinning = kTRUE;
        for (Int_t j = 0; j<nDataPoints-1; j++){
            cout << graphs[j]->GetX()[i] << "-" << graphs[j+1]->GetX()[i] << endl;
            if (TMath::Abs(graphs[j]->GetX()[i] - graphs[j+1]->GetX()[i]) > 0.0001 )
                isSameBinning = kFALSE;
        }
        if (!isSameBinning){
            cout << "graphs are not in same binning" << endl;
            return 0x0;
        }
        
        // set data points for interpolation
        TGraphErrors *gToFit = new TGraphErrors(nDataPoints);
        for (Int_t j = 0; j < nDataPoints; j++){
            gToFit->SetPoint(j, energies[j], graphs[j]->GetY()[i]);
            gToFit->SetPointError(j, 0, graphs[j]->GetEYhigh()[i]);
        }    
//         gToFit->Print();
        
        TF1 *fPowerlaw = new TF1("fPowerlaw","[0]*x^([1])", 0,20000);
        if(i==0){
            fPowerlaw->SetParameters(0, 0.005);
            fPowerlaw->SetParameters(1, 0.13);
        }else{
            fPowerlaw->SetParameters(0, 0.1);
            fPowerlaw->SetParameters(1, 2.0);
        }

        for(Int_t l = 0; l < 10; l++) gToFit->Fit(fPowerlaw,"Q");
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint, 0.68);
        Double_t alpha  = fPowerlaw->GetParameter(1);
        cout << "pT: ";
        for (Int_t j = 0; j< nDataPoints; j++){
             cout << graphs[j]->GetX()[i] << " % " ;
        }
        cout << ": " << fPowerlaw->GetParameter(0) << " - " << fPowerlaw->GetParameter(1) << endl;
        
        gInterpol->SetPoint(i, graphs[0]->GetX()[i],fPowerlaw->Eval(dSqrts));
        gInterpol->SetPointError(i, graphs[0]->GetEXlow()[i], graphs[0]->GetEXhigh()[i], grint->GetEY()[0],grint->GetEY()[0]);

        gAlpha->SetPoint(i, graphs[0]->GetX()[i],alpha);
        //gAlpha->SetPointError(i, 0,alphaE);

        gPtvsSqrts[i]= new TGraphErrors(1);
        gPtvsSqrts[i]->SetPoint(0,dSqrts,gInterpol->GetY()[i]);
        //gPtvsSqrts[i]->Sort();

        fPowerlawFits[i] = fPowerlaw;
        gPtvsEnergies[i] = gToFit;

        gToFit->SetPoint(nDataPoints, dSqrts, gInterpol->GetY()[i]);
        gToFit->SetPointError(nDataPoints, 0, gInterpol->GetEYhigh()[i]);

        delete grint;
    }

    graphAlpha          = gAlpha;
    graphPtvsSqrts      = new TGraphErrors*[graphs[0]->GetN()];
    fPowerlawSystem	    = new TF1*[graphs[0]->GetN()];
    gPtvsEnergiesSystem = new TGraphErrors*[graphs[0]->GetN()];

    for ( Int_t i = 0; i < graphs[0]->GetN(); i++ ){
        graphPtvsSqrts[i]       = gPtvsSqrts[i];
        fPowerlawSystem[i]      = fPowerlawFits[i];
        gPtvsEnergiesSystem[i]  = gPtvsEnergies[i];
    }

    return gInterpol;
}


//________________________________________________________________________________________________________________________
TH1D* ConvertYieldHisto(TH1D* input, Bool_t DivideBy2pi, Bool_t DivideByPt, Bool_t MultiplyBy2pi, Bool_t MultiplyByPt){
    if (!input) {
        cout << "Error: Histogram is NULL" << endl;
        return NULL;
    }

    Int_t nBins                 = input->GetNbinsX();
    Double_t newValue           = 0;
    Double_t newErrorValue      = 0;
    Double_t correctionValue    = 1;

    //correct by 2pi if specified
    if (DivideBy2pi) input->Scale(1/(2*TMath::Pi()));
    if (MultiplyBy2pi) input->Scale(2*TMath::Pi());

    for(Int_t i=0;i<nBins;i++){

        //correct by 1/Pt if specified
        if(DivideByPt)    correctionValue  = 1/(input->GetBinCenter(i+1));
        if(MultiplyByPt)  correctionValue  = input->GetBinCenter(i+1);

        //set the value and error of the bin
        input->SetBinContent(i+1,input->GetBinContent(i+1)*correctionValue);
        input->SetBinError(i+1,input->GetBinError(i+1)*correctionValue);
    }

    return input;
}

//________________________________________________________________________________________________________________________
void PlotInterpolationPtBins(TGraphErrors** gPtvSqrts,TGraphErrors** gPtvsEnergies, TF1** fPowerlaw, TGraphAsymmErrors* gRpPb,Int_t fColumnPlot, Int_t fRowPlot,TString namePlot){

    TGaxis::SetMaxDigits(3);
    TString nameCanvas = "";
    TString namePad    = "";

    TCanvas * canvasPtvsSqrts     = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    canvasPtvsSqrts->SetTopMargin(0.00);
    canvasPtvsSqrts->SetBottomMargin(0.00);
    canvasPtvsSqrts->SetRightMargin(0.0);
    canvasPtvsSqrts->SetLeftMargin(0.00);

    TPad * padPtvsSqrts           = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    padPtvsSqrts->SetFillColor(0);
    padPtvsSqrts->GetFrame()->SetFillColor(0);
    padPtvsSqrts->SetBorderMode(0);
    padPtvsSqrts->Divide(fColumnPlot,fRowPlot,0.0,0.0);
    padPtvsSqrts->SetLeftMargin(0.);
    padPtvsSqrts->SetRightMargin(0.);
    padPtvsSqrts->SetTopMargin(0.);
    padPtvsSqrts->SetBottomMargin(0.);
    padPtvsSqrts->Draw();

    cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
    Int_t place = 0;
    for(Int_t iPt = 0; iPt <gRpPb->GetN(); iPt++){
        place++;
        padPtvsSqrts->cd(place);
        padPtvsSqrts->cd(place)->SetTopMargin(0.12);
        padPtvsSqrts->cd(place)->SetBottomMargin(0.15);
        padPtvsSqrts->cd(place)->SetRightMargin(0.035);
        //padPtvsSqrts->cd(place)->SetLeftMargin(0.25);

        int remaining           = (place-1)%fColumnPlot;
        if (remaining > 0) padPtvsSqrts->cd(place)->SetLeftMargin(0.15);
        else padPtvsSqrts->cd(place)->SetLeftMargin(0.25);
        DrawGammaSetMarkerTGraphErr(gPtvSqrts[iPt],21,1.5, kRed , kRed);
        DrawGammaSetMarkerTGraphErr(gPtvsEnergies[iPt],20,1.5, kBlack , kBlack);

        gPtvsEnergies[iPt]->GetXaxis()->SetTitle("#sqrt{s}");
        gPtvsEnergies[iPt]->GetXaxis()->SetTitleSize(0.065);
        gPtvsEnergies[iPt]->GetXaxis()->SetLabelSize(0.06);
        gPtvsEnergies[iPt]->GetYaxis()->SetTitle("invariant cross section");
        gPtvsEnergies[iPt]->GetYaxis()->SetTitleSize(0.065);
        gPtvsEnergies[iPt]->GetYaxis()->SetLabelSize(0.06);
        gPtvsEnergies[iPt]->GetXaxis()->SetNdivisions(308,kTRUE);
        gPtvsEnergies[iPt]->GetYaxis()->SetNdivisions(304,kTRUE);
        gPtvsEnergies[iPt]->GetXaxis()->SetLabelOffset(0.015);
        gPtvsEnergies[iPt]->GetYaxis()->SetLabelOffset(0.01);

        gPtvsEnergies[iPt]->Draw("ap");
        gPtvSqrts[iPt]->Draw("same, p");
        fPowerlaw[iPt]->SetLineColor(kBlue);
        fPowerlaw[iPt]->SetLineWidth(2);
        fPowerlaw[iPt]->Draw("same");

        TString Title = Form("#it{p}_{T} = %3.2f GeV/#it{c} ",gRpPb->GetX()[iPt]);

        if(Title.Length() > 0){
        gPtvsEnergies[iPt]->SetTitle("");

        Double_t xMin,yMin;
        yMin = 0.78;
        if ( remaining > 0 ) xMin = 0.20;
        else xMin = 0.30;

        TLatex *alice = new TLatex(xMin,yMin,Form("%s",Title.Data())); // Bo: this was
        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextSize(0.062);
        alice->Draw();
        }
    }

    canvasPtvsSqrts->Print(namePlot.Data());
    delete padPtvsSqrts;
    delete canvasPtvsSqrts;
}

//________________________________________________________________________________________________________________________
void PlotAlphavsPt(TGraphErrors* gAlpha, TString method, TString thesisPlotLabel, TString namePlot){

    TCanvas * canvasAlphavsPt     = new TCanvas("AlphavsPt","",1400,900);  // gives the page size
    DrawGammaCanvasSettings( canvasAlphavsPt,  0.08, 0.02, 0.03, 0.09);
    SetStyleTGraphErrorsForGraphs(gAlpha,"#it{p}_{T} GeV/#it{c}","#alpha", 0.04,0.04, 0.04,0.04, 1.,1., 512, 512);
    DrawGammaSetMarkerTGraphErr(gAlpha,21,1.5, kRed+2 , kRed+2);
    gAlpha->GetYaxis()->CenterTitle();
    gAlpha->Draw("ap");

    TLatex *labelThesis = new TLatex(0.15,0.90,thesisPlotLabel.Data());
    SetStyleTLatex( labelThesis, 0.04,4);
    labelThesis->Draw();

    canvasAlphavsPt->Print(namePlot.Data());
    delete canvasAlphavsPt;
}


//________________________________________________________________________________________________________________________
TF1* DoFitWithTsallis(TGraph* graph, TString name, TString particle, Double_t p0, Double_t p1, Double_t p2){

    cout << "-----------------------------------" << endl;
    cout << "fit: '" << name.Data() << "' for '" << particle.Data() << "'" << endl;

    Double_t paramFit[3] = {p0, p1, p2};
    TF1* fit = FitObject("l",name.Data(),particle.Data(),graph,graph->GetX()[0],graph->GetX()[graph->GetN()-1],paramFit,"QNRMEX0+");

    cout << "chi2/ndf: " << fit->GetChisquare()/fit->GetNDF() << endl;
    cout << endl;

    paramFit[0]=fit->GetParameter(0);
    paramFit[1]=fit->GetParameter(1);
    paramFit[2]=fit->GetParameter(2);
    fit = FitObject("l",name.Data(),particle.Data(),graph,graph->GetX()[0],graph->GetX()[graph->GetN()-1],paramFit,"QNRMEX0+");

    cout << WriteParameterToFile(fit) << endl;
    cout << "-----------------------------------" << endl;
    cout << endl;
    return fit;
}

//________________________________________________________________________________________________________________________
TF1* DoFitWithTCM(TGraph* graph, TString name, TString particle, Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4 ){

    cout << "-----------------------------------" << endl;
    cout << "fit: '" << name.Data() << "' for '" << particle.Data() << "'" << endl;

    Double_t paramFit[5] = {p0, p1, p2, p3, p4};
    TF1* fit = FitObject("tcm",name.Data(),particle.Data(),graph,graph->GetX()[0],graph->GetX()[graph->GetN()-1],paramFit,"QNRMEX0+");

    cout << "chi2/ndf: " << fit->GetChisquare()/fit->GetNDF() << endl;
    cout << endl;

    paramFit[0]=fit->GetParameter(0);
    paramFit[1]=fit->GetParameter(1);
    paramFit[2]=fit->GetParameter(2);
    paramFit[3]=fit->GetParameter(3);
    paramFit[4]=fit->GetParameter(4);
    fit = FitObject("tcm",name.Data(),particle.Data(),graph,graph->GetX()[0],graph->GetX()[graph->GetN()-1],paramFit,"QNRMEX0+");

    cout << WriteParameterToFile(fit) << endl;
    cout << "-----------------------------------" << endl;
    cout << endl;
    return fit;
}

//________________________________________________________________________________________________________________________
void CalculateStatPlusSysErrors(TH1D* histStat, TH1D* histSys){
    for(Int_t i=1; i<=histStat->GetNbinsX(); i++){
        histStat->SetBinError(i,TMath::Sqrt(TMath::Power(histStat->GetBinError(i),2)+TMath::Power(histSys->GetBinError(i),2)));
    }
    return;
}

//________________________________________________________________________________________________________________________
TGraphAsymmErrors* CalculateStatPlusSysErrors(TGraphAsymmErrors* graphStat, TGraphAsymmErrors* graphSys){
    Double_t* xValue            = graphStat->GetX();
    Double_t* yValue            = graphStat->GetY();
    Double_t* xErrorHigh        = graphStat->GetEXhigh();
    Double_t* yErrorHigh        = graphStat->GetEYhigh();
    Double_t* xErrorLow         = graphStat->GetEXlow();
    Double_t* yErrorLow         = graphStat->GetEYlow();
    Int_t nPoints               = graphStat->GetN();

    TGraphAsymmErrors* graphR   = new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow, yErrorHigh);

    Double_t* yErrorRHigh       = graphR->GetEYhigh();
    Double_t* yErrorRLow        = graphR->GetEYlow();
    Double_t* yErrorSysHigh     = graphSys->GetEYhigh();
    Double_t* yErrorSysLow      = graphSys->GetEYlow();

    for(Int_t i=1; i<nPoints; i++){
        yErrorRHigh[i]           = TMath::Sqrt(TMath::Power(yErrorHigh[i],2)+TMath::Power(yErrorSysHigh[i],2)); 
        yErrorRLow[i]           = TMath::Sqrt(TMath::Power(yErrorLow[i],2)+TMath::Power(yErrorSysLow[i],2)); 
    }
    return graphR;
}

//________________________________________________________________________________________________________________________
void PlotWithFit(TCanvas* canvas, TH2F* hist, TH2F* hist2, TGraphAsymmErrors* graph, TGraphAsymmErrors* graphRebin, TF1* fit, TString name, TString outputDir, TString suffix, TString energyCur){

    canvas->Clear();
    hist->DrawCopy();

    DrawGammaSetMarkerTGraphAsym(graph, GetDefaultMarkerStyle(energyCur.Data(),"", ""), GetDefaultMarkerSize(energyCur.Data(),"", ""), GetColorDefaultColor( energyCur.Data(),"", ""), GetColorDefaultColor( energyCur.Data(),"", ""), 1.4, kTRUE);
    graph->Draw("pEsame");
    if(graphRebin){
        DrawGammaSetMarkerTGraphAsym(graphRebin, GetDefaultMarkerStyle(energyCur.Data(),"MC", ""), GetDefaultMarkerSize(energyCur.Data(),"MC", ""), kGray+1, kGray+1, 1.4, kTRUE);
        graphRebin->Draw("pEsame");
    }

    fit->SetLineColor(GetColorDefaultColor( energyCur.Data(),"", ""));
    fit->Draw("same");

    canvas->Update();
    canvas->Print(Form("%s/%s.%s",outputDir.Data(),name.Data(),suffix.Data()));

    canvas->SetLogy(kFALSE);
    if(graphRebin){
        canvas->Clear();
        hist2->GetYaxis()->SetRangeUser(0.5,1.5);
        hist2->DrawCopy();

        TGraphAsymmErrors* graphRatio        = CalculateGraphErrRatioToFit (graph, fit);
        TGraphAsymmErrors* graphRebinRatio   = CalculateGraphErrRatioToFit (graphRebin, fit);

        DrawGammaSetMarkerTGraphAsym(graphRatio, GetDefaultMarkerStyle(energyCur.Data(),"", ""), GetDefaultMarkerSize(energyCur.Data(),"", ""), GetColorDefaultColor( energyCur.Data(),"", ""), GetColorDefaultColor( energyCur.Data(),"", ""), 1.4, kTRUE);
        graphRatio->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphRebinRatio, GetDefaultMarkerStyle(energyCur.Data(),"MC", ""), GetDefaultMarkerSize(energyCur.Data(),"MC", ""), kGray+1, kGray+1, 1.4, kTRUE);
        graphRebinRatio->Draw("pEsame");
        DrawGammaLines(hist->GetXaxis()->GetBinCenter(1), hist->GetXaxis()->GetBinCenter(hist->GetNbinsX()) , 1.0, 1.0,0.1, kGray, 1);

        canvas->Update();
        canvas->Print(Form("%s/%s_ratio.%s",outputDir.Data(),name.Data(),suffix.Data()));
    } else {
        canvas->Clear();
        hist2->GetYaxis()->SetRangeUser(0,2.);
        hist2->DrawCopy();

        TGraphAsymmErrors* graphRatio       = CalculateGraphErrRatioToFit (graph, fit);

        DrawGammaSetMarkerTGraphAsym(graphRatio, GetDefaultMarkerStyle(energyCur.Data(),"", ""), GetDefaultMarkerSize(energyCur.Data(),"", ""), GetColorDefaultColor( energyCur.Data(),"", ""), GetColorDefaultColor( energyCur.Data(),"", ""), 1.4, kTRUE);
        graphRatio->Draw("pEsame");
        DrawGammaLines(hist->GetXaxis()->GetBinCenter(1), hist->GetXaxis()->GetBinCenter(hist->GetNbinsX()) , 1.0, 1.0,0.1, kGray, 1);

        canvas->Update();
        canvas->Print(Form("%s/%s_ratio.%s",outputDir.Data(),name.Data(),suffix.Data()));
    }
    canvas->SetLogy();

    return;
}

//________________________________________________________________________________________________________________________
void PlotWithFit(TCanvas* canvas, TH2F* hist, TH2F* hist2, TGraphErrors* graph, TGraphErrors* graphRebin, TF1* fit, TString name, TString outputDir, TString suffix, TString energyCur){

    canvas->Clear();
    hist->DrawCopy();

    DrawGammaSetMarkerTGraphErr(graph, GetDefaultMarkerStyle(energyCur.Data(),"", ""), GetDefaultMarkerSize(energyCur.Data(),"", ""), GetColorDefaultColor( energyCur.Data(),"", ""), GetColorDefaultColor( energyCur.Data(),"", ""), 1.4, kTRUE);
    graph->Draw("pEsame");
    if(graphRebin){
        DrawGammaSetMarkerTGraphErr(graphRebin, GetDefaultMarkerStyle(energyCur.Data(),"MC", ""), GetDefaultMarkerSize(energyCur.Data(),"MC", ""), kGray+1, kGray+1, 1.4, kTRUE);
        graphRebin->Draw("pEsame");
    }

    fit->SetLineColor(GetColorDefaultColor( energyCur.Data(),"", ""));
    fit->Draw("same");

    canvas->Update();
    canvas->Print(Form("%s/%s.%s",outputDir.Data(),name.Data(),suffix.Data()));

    canvas->SetLogy(kFALSE);
    if(graphRebin){
        canvas->Clear();
        hist2->GetYaxis()->SetRangeUser(0.5,1.5);
        hist2->DrawCopy();

        TGraphErrors* graphRatio        = CalculateGraphErrRatioToFit (graph, fit);
        TGraphErrors* graphRebinRatio   = CalculateGraphErrRatioToFit (graphRebin, fit);

        DrawGammaSetMarkerTGraphErr(graphRatio, GetDefaultMarkerStyle(energyCur.Data(),"", ""), GetDefaultMarkerSize(energyCur.Data(),"", ""), GetColorDefaultColor( energyCur.Data(),"", ""), GetColorDefaultColor( energyCur.Data(),"", ""), 1.4, kTRUE);
        graphRatio->Draw("pEsame");
        DrawGammaSetMarkerTGraphErr(graphRebinRatio, GetDefaultMarkerStyle(energyCur.Data(),"MC", ""), GetDefaultMarkerSize(energyCur.Data(),"MC", ""), kGray+1, kGray+1, 1.4, kTRUE);
        graphRebinRatio->Draw("pEsame");
        DrawGammaLines(hist->GetXaxis()->GetBinCenter(1), hist->GetXaxis()->GetBinCenter(hist->GetNbinsX()) , 1.0, 1.0,0.1, kGray, 1);

        canvas->Update();
        canvas->Print(Form("%s/%s_ratio.%s",outputDir.Data(),name.Data(),suffix.Data()));
    } else {
        canvas->Clear();
        hist2->GetYaxis()->SetRangeUser(0,2.);
        hist2->DrawCopy();

        TGraphErrors* graphRatio        = CalculateGraphErrRatioToFit (graph, fit);

        DrawGammaSetMarkerTGraphErr(graphRatio, GetDefaultMarkerStyle(energyCur.Data(),"", ""), GetDefaultMarkerSize(energyCur.Data(),"", ""), GetColorDefaultColor( energyCur.Data(),"", ""), GetColorDefaultColor( energyCur.Data(),"", ""), 1.4, kTRUE);
        graphRatio->Draw("pEsame");
        DrawGammaLines(hist->GetXaxis()->GetBinCenter(1), hist->GetXaxis()->GetBinCenter(hist->GetNbinsX()) , 1.0, 1.0,0.1, kGray, 1);

        canvas->Update();
        canvas->Print(Form("%s/%s_ratio.%s",outputDir.Data(),name.Data(),suffix.Data()));
    }
    canvas->SetLogy();

    return;
}

//________________________________________________________________________________________________________________________
void PlotGraphsOfAllEnergies( TCanvas* canvas, 
                              TH2F* hist, 
                              Int_t nDataSets,
                              TGraphAsymmErrors** graphs,
                              TF1** fits, 
                              TString* energies,
                              TString name, 
                              TString outputDir, 
                              TString suffix
                            ){

    canvas->Clear();
    hist->DrawCopy();

    TLegend* legendEnergies     = GetAndSetLegend2(0.18, 0.12, 0.5, 0.12+(0.035*nDataSets), 32, 2, "", 43, 0.2);
    for (Int_t i = 0; i< nDataSets; i++){
        cout << energies[i].Data() << endl;
        DrawGammaSetMarkerTGraphAsym(graphs[i], GetDefaultMarkerStyle(energies[i].Data(),"", ""), GetDefaultMarkerSize(energies[i].Data(),"", ""), GetColorDefaultColor( energies[i].Data(),"", ""), GetColorDefaultColor( energies[i].Data(),"", ""), 1.4, kTRUE);
        graphs[i]->Draw("pEsame");
        fits[i]->SetLineColor(GetColorDefaultColor( energies[i].Data(),"", ""));
        fits[i]->Draw("same");
        legendEnergies->AddEntry(graphs[i],energies[i],"pf");
        legendEnergies->AddEntry(fits[i],"fit","l");
    }     
    legendEnergies->Draw();
    
    canvas->Update();
    canvas->Print(Form("%s/%s.%s",outputDir.Data(),name.Data(),suffix.Data()));

    return;
}

//________________________________________________________________________________________________________________________
void SetStyleTGraphErrorsForGraphs(TGraph* graph, TString XTitle, TString YTitle, Size_t xLableSize, Size_t xTitleSize, Size_t yLableSize, Size_t yTitleSize, Float_t xTitleOffset, Float_t yTitleOffset, Int_t xNDivisions, Int_t yNDivisions){
    graph->GetXaxis()->SetTitle(XTitle);
    graph->GetYaxis()->SetTitle(YTitle);
    graph->SetTitle("");

    graph->GetXaxis()->SetLabelSize(xLableSize);
    graph->GetXaxis()->SetTitleSize(xTitleSize);
    graph->GetXaxis()->SetTitleOffset(xTitleOffset);
    graph->GetXaxis()->SetNdivisions(xNDivisions,kTRUE);

    graph->GetXaxis()->SetLabelFont(42);
    graph->GetYaxis()->SetLabelFont(42);
    graph->GetXaxis()->SetTitleFont(62);
    graph->GetYaxis()->SetTitleFont(62);


    graph->GetYaxis()->SetDecimals();
    graph->GetYaxis()->SetLabelSize(yLableSize);
    graph->GetYaxis()->SetTitleSize(yTitleSize);
    graph->GetYaxis()->SetTitleOffset(yTitleOffset);
    graph->GetYaxis()->SetNdivisions(yNDivisions,kTRUE);
}

//________________________________________________________________________________________________________________________
TGraphAsymmErrors* SetGraphErrors (TGraphAsymmErrors* graph, TGraphAsymmErrors* refGraph){

    Double_t* xValue            = graph->GetX();
    Double_t* yValue            = graph->GetY();
    Double_t* xErrorHigh        = graph->GetEXhigh();
    Double_t* yErrorHigh        = graph->GetEYhigh();
    Double_t* xErrorLow         = graph->GetEXlow();
    Double_t* yErrorLow         = graph->GetEYlow();
    Int_t nPoints               = graph->GetN();

    Double_t* yValueRef         = refGraph->GetY();
    Double_t* yErrorRefHigh     = refGraph->GetEYhigh();
    Double_t* yErrorRefLow      = refGraph->GetEYlow();

    for (Int_t i = 0; i < nPoints; i++){
        yErrorHigh[i]           = (yErrorRefHigh[i]/yValueRef[i])*yValue[i];
        yErrorLow[i]            = (yErrorRefLow[i]/yValueRef[i])*yValue[i];
    }

    graph                       = new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow, yErrorHigh);
    return graph;
}

//________________________________________________________________________________________________________________________
TGraphAsymmErrors* RemoveSysErrors (TGraphAsymmErrors* graph, TGraphAsymmErrors* refGraph){

    Double_t* xValue            = graph->GetX();
    Double_t* yValue            = graph->GetY();
    Double_t* xErrorHigh        = graph->GetEXhigh();
    Double_t* yErrorHigh        = graph->GetEYhigh();
    Double_t* xErrorLow         = graph->GetEXlow();
    Double_t* yErrorLow         = graph->GetEYlow();
    Int_t nPoints               = graph->GetN();

    TGraphAsymmErrors* graphR   = new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow, yErrorHigh);

    Double_t* yErrorRHigh       = graphR->GetEYhigh();
    Double_t* yErrorRLow        = graphR->GetEYlow();
    Int_t nPointsR              = graphR->GetN();
    Double_t* yErrorRefHigh     = refGraph->GetEYhigh();
    Double_t* yErrorRefLow      = refGraph->GetEYlow();
    
    for (Int_t i = 0; i < nPointsR; i++){
      if(((yErrorRLow[i]*yErrorRLow[i])-(yErrorRefLow[i]*yErrorRefLow[i]))<0.)
        yErrorRLow[i]               = yErrorRefLow[i]*yErrorRefLow[i];
      else
        yErrorRLow[i]               = TMath::Sqrt((yErrorRLow[i]*yErrorRLow[i])-(yErrorRefLow[i]*yErrorRefLow[i]));
      if(((yErrorRHigh[i]*yErrorRHigh[i])-(yErrorRefHigh[i]*yErrorRefHigh[i]))<0.)
        yErrorRHigh[i]               = yErrorRefHigh[i]*yErrorRefHigh[i];
      else
        yErrorRHigh[i]               = TMath::Sqrt((yErrorRHigh[i]*yErrorRHigh[i])-(yErrorRefHigh[i]*yErrorRefHigh[i]));
    }

    return graphR;
}

