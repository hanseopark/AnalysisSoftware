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
#include "CombineMesonMeasurementsPP.h"

extern TRandom*    gRandom;
extern TBenchmark* gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*    gMinuit;

void CombineMesonMeasurementsPP()
{
    //__________________________________________ Defining input files
    TString fileName[6];
    /*  900GeV  */  fileName[0]                 = "CombinationInputPP/900GeV/CombinedResultsPaperPP900GeV_2017_07_18.root";
    /*  2.76TeV */  fileName[1]                 = "CombinationInputPP/2.76TeV/CombinedResultsPaperPP2760GeV_2017_07_10_FrediV2Clusterizer.root";
    /*  5TeV    */  fileName[2]                 = "CombinationInputPP/5TeV/";
    /*  7TeV    */  fileName[3]                 = "CombinationInputPP/7TeV/CombinedResultsPaperPP7TeV_2017_07_18.root";
    /*  8TeV    */  fileName[4]                 = "CombinationInputPP/8TeV/CombinedResultsPaperPP8TeV_2017_07_13.root";
    /*  13TeV   */  fileName[5]                 = "CombinationInputPP/13TeV/";
    TString         fileNameTheory              = "ExternalInput/Theory/TheoryCompilationPP.root";
    TString         suffix                      = "eps";
    Int_t includeEnergy[6]                      = {1,1,0,1,1,0};
    Int_t numActiveMeas                         = 4;
    Double_t maxX                               = 49.9;
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    StyleSettingsThesis();
    SetPlotStyle();
    Int_t exampleActiveMeas                     = 0;
    //__________________________________________ Producing date strings
    TString date                                = ReturnDateString();
    TString dateForOutput                       = ReturnDateStringForOutput();
    cout << dateForOutput.Data() << endl;
    Double_t mesonMassExpectPi0                 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    Double_t mesonMassExpectEta                 = TDatabasePDG::Instance()->GetParticle(221)->Mass();

    //__________________________________________ Copying input files to output directory
    TString outputDir                           = Form("%s/%s/CombineMesonMeasurementsPP",suffix.Data(),dateForOutput.Data());
    cout << outputDir.Data() << endl;

    gSystem->Exec("mkdir -p "+outputDir);
    for(Int_t i=0; i<6; i++){
        if(includeEnergy[i]){
            gSystem->Exec(Form("cp %s %s/Input%s.root", fileName[i].Data(),  outputDir.Data(),nameEnergyGlobal[i].Data()));
        }
    }
    gSystem->Exec(Form("cp %s %s/Theory.root",fileNameTheory.Data(),  outputDir.Data()));

    //__________________________________________ Loading input files and directories
    for(Int_t i=0; i<6; i++){
        if(includeEnergy[i]){
            inputFile[i]                        = new TFile(fileName[i].Data());
            directoryPi0[i]                     = (TDirectory*)inputFile[i]->Get(Form("Pi0%s",nameEnergyGlobal[i].Data()));
            directoryEta[i]                     = (TDirectory*)inputFile[i]->Get(Form("Eta%s",nameEnergyGlobal[i].Data()));
            if(inputFile[i])
                cout << "Found " << nameEnergyGlobal[i] << " input and directories " << Form("Pi0%s",nameEnergyGlobal[i].Data()) << " and " << Form("Eta%s",nameEnergyGlobal[i].Data()) << endl;
            exampleActiveMeas                   = i;
        }
    }
    TFile* fileTheoryCompilation                            = new TFile(fileNameTheory.Data());

    //__________________________________________ Definition of colors, styles and markers sizes
    LoadColorsMarkersAndSizes();

    //__________________________________________ Loading graphs from input files || [6] = energy, [11] = reco. method
    TGraphAsymmErrors* graphPi0InvariantCrossSectionStat[6][11];
    TGraphAsymmErrors* graphPi0InvariantCrossSectionSys[6][11];
    TGraphAsymmErrors* graphEtaInvariantCrossSectionStat[6][11];
    TGraphAsymmErrors* graphEtaInvariantCrossSectionSys[6][11];
    TH1D* histoEtaToPi0RatioStat[6][11];
    TGraphAsymmErrors* graphEtaToPi0RatioStat[6][11];
    TGraphAsymmErrors* graphEtaToPi0RatioSys[6][11];
    Double_t minXSpectra[6][11];
    Double_t maxXSpectra[6][11];
    Double_t minXSpectraEta[6][11];
    Double_t maxXSpectraEta[6][11];
    Bool_t recoMethodFoundPi0[6][11];
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            for (Int_t j = 0; j < 11; j++){


              if(j==3){
                graphPi0InvariantCrossSectionStat[i][j] = NULL;
                graphPi0InvariantCrossSectionSys[i][j] = NULL;
                graphEtaInvariantCrossSectionStat[i][j]=NULL;
                graphEtaInvariantCrossSectionSys[i][j]=NULL;
                continue; //skip PCM-PHOS for the moment
              }

                //______________________________ Loading pi0 inv. cross sections
                graphPi0InvariantCrossSectionStat[i][j]   = (TGraphAsymmErrors*)directoryPi0[i]->Get(Form("graphInvCrossSectionPi0%s%s%sStatErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data(),graphNameModifier[j].Data()));
                graphPi0InvariantCrossSectionSys[i][j]   = (TGraphAsymmErrors*)directoryPi0[i]->Get(Form("graphInvCrossSectionPi0%s%s%sSysErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data(),graphNameModifier[j].Data()));
                if(graphPi0InvariantCrossSectionStat[i][j]&&graphPi0InvariantCrossSectionSys[i][j]){
                    cout << "found " << Form("graphInvCrossSectionPi0%s%sStatErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data()) << endl;
                    //__________________________ Reading global xRange for each spectrum
                    minXSpectra[i][j]           = graphPi0InvariantCrossSectionStat[i][j]->GetX()[0]-graphPi0InvariantCrossSectionStat[i][j]->GetEXlow()[0];
                    maxXSpectra[i][j]           = graphPi0InvariantCrossSectionStat[i][j]->GetX()[graphPi0InvariantCrossSectionStat[i][j]->GetN()-1]+graphPi0InvariantCrossSectionStat[i][j]->GetEXhigh()[graphPi0InvariantCrossSectionStat[i][j]->GetN()-1];
                    cout << "minx: " << minXSpectra[i][j] << "\tmaxX: " << maxXSpectra[i][j] << endl;
                }else{
                  graphPi0InvariantCrossSectionStat[i][j]=NULL;
                  graphPi0InvariantCrossSectionSys[i][j]=NULL;
                }

                if(i==0 && j==10){
                  nameMeasGlobal[j]                  = "PCM";
                  nameMeasGlobalshort[j]             = "PCM";
                  graphNameModifier[j]               = "";
                }

                //______________________________ Loading eta inv. cross sections
                graphEtaInvariantCrossSectionStat[i][j]   = (TGraphAsymmErrors*)directoryEta[i]->Get(Form("graphInvCrossSectionEta%s%s%sStatErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data(),graphNameModifier[j].Data()));
                graphEtaInvariantCrossSectionSys[i][j]   = (TGraphAsymmErrors*)directoryEta[i]->Get(Form("graphInvCrossSectionEta%s%s%sSysErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data(),graphNameModifier[j].Data()));
                if(graphEtaInvariantCrossSectionStat[i][j]&&graphEtaInvariantCrossSectionSys[i][j]){
                    cout << "found " << Form("graphInvCrossSectionEta%s%sStatErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data()) << endl;
                    minXSpectraEta[i][j]           = graphEtaInvariantCrossSectionStat[i][j]->GetX()[0]-graphEtaInvariantCrossSectionStat[i][j]->GetEXlow()[0];
                    maxXSpectraEta[i][j]           = graphEtaInvariantCrossSectionStat[i][j]->GetX()[graphEtaInvariantCrossSectionStat[i][j]->GetN()-1]+graphEtaInvariantCrossSectionStat[i][j]->GetEXhigh()[graphEtaInvariantCrossSectionStat[i][j]->GetN()-1];
                }else{
                  graphEtaInvariantCrossSectionStat[i][j]=NULL;
                  graphEtaInvariantCrossSectionSys[i][j]=NULL;
                }
                //______________________________ Loading eta/pi0 ratios
                graphEtaToPi0RatioStat[i][j]   = (TGraphAsymmErrors*)directoryEta[i]->Get(Form("graphRatioEtaToPi0%s%sStatErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data()));
                histoEtaToPi0RatioStat[i][j]   = (TH1D*)directoryEta[i]->Get(Form("histoRatioEtaToPi0%s%sStatErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data()));
                if(histoEtaToPi0RatioStat[i][j])
                    graphEtaToPi0RatioStat[i][j]= new TGraphAsymmErrors(histoEtaToPi0RatioStat[i][j]);
                graphEtaToPi0RatioSys[i][j]   = (TGraphAsymmErrors*)directoryEta[i]->Get(Form("graphRatioEtaToPi0%s%sSysErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data()));
                if(graphEtaToPi0RatioStat[i][j])
                    cout << "found " << Form("graphRatioEtaToPi0%s%sStatErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data()) << endl;
                if(i==0 && j==10){
                  nameMeasGlobal[j]                  = "Comb";
                  nameMeasGlobalshort[j]             = "Comb";
                  graphNameModifier[j]               = "A";
                }
            }
        }
    }
    TString pythiaHistoNames[6]                 = { "histoInvSecPythia8Monash2013Pi0900GeV", "histoInvSecPythia8Monash2013Pi02760GeV", "", "histoInvSecPythia8Monash2013Pi07TeV", "histoInvSecPythia8Monash2013LegoPi08TeV", ""};
    TH1F* histoPythiaInvCrossSectionPi0[6];
    TGraphErrors* graphPythiaInvCrossSectionPi0[6];

    TString pythiaHistoNamesEta[6]                 = { "histoInvSecPythia8Monash2013Eta900GeV", "histoInvSecPythia8Monash2013Eta2760GeV", "", "histoInvSecPythia8Monash2013Eta7TeV", "histoInvSecPythia8Monash2013LegoEta8TeV", ""};
    TH1F* histoPythiaInvCrossSectionEta[6];
    TGraphErrors* graphPythiaInvCrossSectionEta[6];
    Double_t scalingFactorPythia = 1;

    for(Int_t i=0;i<6;i++){
        if(includeEnergy[i]){
            histoPythiaInvCrossSectionPi0[i]    = (TH1F*) fileTheoryCompilation->Get(pythiaHistoNames[i].Data());
            if(histoPythiaInvCrossSectionPi0[i]){
                histoPythiaInvCrossSectionPi0[i]->GetXaxis()->SetRangeUser(0.3,25);
                histoPythiaInvCrossSectionPi0[i]->Scale(scalingFactorPythia);
                DrawGammaSetMarker(histoPythiaInvCrossSectionPi0[i], 24, 1.5, pythia8color , pythia8color);
                histoPythiaInvCrossSectionPi0[i]->SetLineWidth(widthCommonFit);
                graphPythiaInvCrossSectionPi0[i]= new TGraphErrors((TH1F*) fileTheoryCompilation->Get(pythiaHistoNames[i].Data()));
                for (int j=0;j<graphPythiaInvCrossSectionPi0[i]->GetN();j++){
                    graphPythiaInvCrossSectionPi0[i]->GetY()[j] *= scalingFactorPythia;
                    graphPythiaInvCrossSectionPi0[i]->GetEY()[j] *= scalingFactorPythia;
                }
                while(graphPythiaInvCrossSectionPi0[i]->GetX()[0] < 0.3) graphPythiaInvCrossSectionPi0[i]->RemovePoint(0);
                while(graphPythiaInvCrossSectionPi0[i]->GetX()[graphPythiaInvCrossSectionPi0[i]->GetN()-1] > 35.) graphPythiaInvCrossSectionPi0[i]->RemovePoint(graphPythiaInvCrossSectionPi0[i]->GetN()-1);
                DrawGammaSetMarkerTGraphErr(graphPythiaInvCrossSectionPi0[i], 0, 0, pythia8color , pythia8color, widthLinesBoxes, kTRUE, pythia8color);
            }else{
                histoPythiaInvCrossSectionPi0[i]=NULL;
                graphPythiaInvCrossSectionPi0[i]=NULL;
            }
            histoPythiaInvCrossSectionEta[i]    = (TH1F*) fileTheoryCompilation->Get(pythiaHistoNamesEta[i].Data());
            if(histoPythiaInvCrossSectionEta[i]){
                histoPythiaInvCrossSectionEta[i]->GetXaxis()->SetRangeUser(0.3,25);
                histoPythiaInvCrossSectionEta[i]->Scale(scalingFactorPythia);
                DrawGammaSetMarker(histoPythiaInvCrossSectionEta[i], 24, 1.5, pythia8color , pythia8color);
                histoPythiaInvCrossSectionEta[i]->SetLineWidth(widthCommonFit);
                graphPythiaInvCrossSectionEta[i]= new TGraphErrors((TH1F*) fileTheoryCompilation->Get(pythiaHistoNamesEta[i].Data()));
                for (int j=0;j<graphPythiaInvCrossSectionEta[i]->GetN();j++){
                    graphPythiaInvCrossSectionEta[i]->GetY()[j] *= scalingFactorPythia;
                    graphPythiaInvCrossSectionEta[i]->GetEY()[j] *= scalingFactorPythia;
                }
                while(graphPythiaInvCrossSectionEta[i]->GetX()[0] < 0.3) graphPythiaInvCrossSectionEta[i]->RemovePoint(0);
                while(graphPythiaInvCrossSectionEta[i]->GetX()[graphPythiaInvCrossSectionEta[i]->GetN()-1] > 35.) graphPythiaInvCrossSectionEta[i]->RemovePoint(graphPythiaInvCrossSectionEta[i]->GetN()-1);
                DrawGammaSetMarkerTGraphErr(graphPythiaInvCrossSectionEta[i], 0, 0, pythia8color , pythia8color, widthLinesBoxes, kTRUE, pythia8color);
            }else{
                histoPythiaInvCrossSectionEta[i]=NULL;
                graphPythiaInvCrossSectionEta[i]=NULL;
            }
            scalingFactorPythia                 *=10;
        }
    }

    //__________________________________________ Scaling cross sections for plotting
    Double_t scalingFactorXsec = 1;
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            for (Int_t j = 0; j < 11; j++){
                if(graphPi0InvariantCrossSectionSys[i][j]){
                    for (int k=0;k<graphPi0InvariantCrossSectionSys[i][j]->GetN();k++){
                        graphPi0InvariantCrossSectionSys[i][j]->GetY()[k]       *= scalingFactorXsec;
                        graphPi0InvariantCrossSectionSys[i][j]->GetEYhigh()[k]  *= scalingFactorXsec;
                        graphPi0InvariantCrossSectionSys[i][j]->GetEYlow()[k]   *= scalingFactorXsec;
                    }
                    DrawGammaSetMarkerTGraphAsym(graphPi0InvariantCrossSectionSys[i][j], markerStyleMeas[j], markerSizeMeas[j], colorMeas[j], colorMeas[j], widthLinesBoxes, kTRUE);
                    if(j==10)
                        DrawGammaSetMarkerTGraphAsym(graphPi0InvariantCrossSectionSys[i][j], markerStyleEnergy[i], markerSizeEnergy[i], colorEnergy[i], colorEnergy[i], widthLinesBoxes, kTRUE);
                }
                if(graphPi0InvariantCrossSectionStat[i][j]){
                    for (int k=0;k<graphPi0InvariantCrossSectionStat[i][j]->GetN();k++){
                        graphPi0InvariantCrossSectionStat[i][j]->GetY()[k]      *= scalingFactorXsec;
                        graphPi0InvariantCrossSectionStat[i][j]->GetEYhigh()[k] *= scalingFactorXsec;
                        graphPi0InvariantCrossSectionStat[i][j]->GetEYlow()[k]  *= scalingFactorXsec;
                    }
                    DrawGammaSetMarkerTGraph(graphPi0InvariantCrossSectionStat[i][j], markerStyleMeas[j], markerSizeMeas[j], colorMeas[j] , colorMeas[j]);
                    if(j==10)
                        DrawGammaSetMarkerTGraph(graphPi0InvariantCrossSectionStat[i][j], markerStyleEnergy[i], markerSizeEnergy[i], colorEnergy[i] , colorEnergy[i]);
                }
                if(graphEtaInvariantCrossSectionSys[i][j]){
                    for (int k=0;k<graphEtaInvariantCrossSectionSys[i][j]->GetN();k++){
                        graphEtaInvariantCrossSectionSys[i][j]->GetY()[k]       *= scalingFactorXsec;
                        graphEtaInvariantCrossSectionSys[i][j]->GetEYhigh()[k]  *= scalingFactorXsec;
                        graphEtaInvariantCrossSectionSys[i][j]->GetEYlow()[k]   *= scalingFactorXsec;
                    }
                    DrawGammaSetMarkerTGraphAsym(graphEtaInvariantCrossSectionSys[i][j], markerStyleMeas[j], markerSizeMeas[j], colorMeas[j], colorMeas[j], widthLinesBoxes, kTRUE);
                    if(j==10)
                        DrawGammaSetMarkerTGraphAsym(graphEtaInvariantCrossSectionSys[i][j], markerStyleEnergy[i], markerSizeEnergy[i], colorEnergy[i], colorEnergy[i], widthLinesBoxes, kTRUE);
                }
                if(graphEtaInvariantCrossSectionStat[i][j]){
                    for (int k=0;k<graphEtaInvariantCrossSectionStat[i][j]->GetN();k++){
                        graphEtaInvariantCrossSectionStat[i][j]->GetY()[k]      *= scalingFactorXsec;
                        graphEtaInvariantCrossSectionStat[i][j]->GetEYhigh()[k] *= scalingFactorXsec;
                        graphEtaInvariantCrossSectionStat[i][j]->GetEYlow()[k]  *= scalingFactorXsec;
                    }
                    DrawGammaSetMarkerTGraph(graphEtaInvariantCrossSectionStat[i][j], markerStyleMeas[j], markerSizeMeas[j], colorMeas[j] , colorMeas[j]);
                    if(j==10)
                        DrawGammaSetMarkerTGraph(graphEtaInvariantCrossSectionStat[i][j], markerStyleEnergy[i], markerSizeEnergy[i], colorEnergy[i] , colorEnergy[i]);
                }
            }
            scalingFactorXsec                   *= 10;
        }
    }


    //__________________________________________ Define canvas, pads and sizes for Spectrum (top) and Ratio (bottom) plot
    Int_t canvasDividerXSec                     = 5+numActiveMeas;
    Double_t textSizeLabelsPixel                = 48;
    Double_t xRangeMinXSec                      = 0.23;
    Double_t xRangeMaxXSec                      = 49;
    Double_t arrayBoundariesX1_XSec[2];
    Double_t arrayBoundariesY1_XSec[canvasDividerXSec];
    Double_t relativeMarginsXXSec[3];
    Double_t relativeMarginsYXSec[3];
    Double_t textsizeLabelsXSec[numActiveMeas+1];
    Double_t textsizeFacXSec[numActiveMeas+1];
    ReturnCorrectValuesForCanvasScaling(1250,4*330+330*numActiveMeas, 1, canvasDividerXSec-1,0.135, 0.005, 0.001,0.04*6/(3+numActiveMeas),arrayBoundariesX1_XSec,arrayBoundariesY1_XSec,relativeMarginsXXSec,relativeMarginsYXSec);
    Double_t marginXSec                         = relativeMarginsXXSec[0]*1250;
    TPad* padXSec[numActiveMeas+1];
    TH2F * histoXSecDummy[numActiveMeas+1];

    TCanvas* canvasXSec                         = new TCanvas("canvasXSec","",0,0,1250,330*canvasDividerXSec);  // gives the page size
    DrawGammaCanvasSettings( canvasXSec,  0.13, 0.02, 0.03, 0.06);

    //__________________________________________ Set spectrum pad properties
    padXSec[0]                         = new TPad("padXSec_0", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[4], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[0],-1, -1, -2);
    DrawGammaPadSettings( padXSec[0], relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[0], relativeMarginsYXSec[1]);
    padXSec[0]->Draw();

    //__________________________________________ Set ratio pad properties
    for(Int_t i=0;i<numActiveMeas;i++){
        padXSec[i+1]                         = new TPad(Form("padXSecFull_%d",i), "",
                                                      arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[4+i+1],
                                                      arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[4+i],-1, -1, -2);
        if(i==numActiveMeas-1)
            DrawGammaPadSettings( padXSec[i+1], relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[2]);
        else
            DrawGammaPadSettings( padXSec[i+1], relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[1]);
        padXSec[i+1]->Draw();
    }
    //__________________________________________ Set pad text sizes
    for(Int_t i=0;i<numActiveMeas+1;i++){
        if (padXSec[i]->XtoPixel(padXSec[i]->GetX2()) < padXSec[i]->YtoPixel(padXSec[i]->GetY1())){
            textsizeLabelsXSec[i]                  = (Double_t)textSizeLabelsPixel/padXSec[i]->XtoPixel(padXSec[i]->GetX2()) ;
            textsizeFacXSec[i]                   = (Double_t)1./padXSec[i]->XtoPixel(padXSec[i]->GetX2()) ;
        } else {
            textsizeLabelsXSec[i]                  = (Double_t)textSizeLabelsPixel/padXSec[i]->YtoPixel(padXSec[i]->GetY1());
            textsizeFacXSec[i]                   = (Double_t)1./padXSec[i]->YtoPixel(padXSec[i]->GetY1());
        }
    }
    //__________________________________________ Draw in spectrum pad
    padXSec[0]->cd();
    padXSec[0]->SetLogy(1);
    padXSec[0]->SetLogx(1);


    histoXSecDummy[0]                           = new TH2F("histoXSecDummy_0","histoXSecDummy_0",11000,xRangeMinXSec,xRangeMaxXSec,1000,1.1e1,9e14);
    SetStyleHistoTH2ForGraphs(histoXSecDummy[0], "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",0.035,0.04, 0.035,0.04, 0.9,1.45);
    histoXSecDummy[0]->GetXaxis()->SetMoreLogLabels();
    histoXSecDummy[0]->GetXaxis()->SetNoExponent(kTRUE);
    SetStyleHistoTH2ForGraphs(histoXSecDummy[0], "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",
                            0.85*textsizeLabelsXSec[0],textsizeLabelsXSec[0], 0.85*textsizeLabelsXSec[0], textsizeLabelsXSec[0], 1,0.2/(textsizeFacXSec[0]*marginXSec));
    histoXSecDummy[0]->GetXaxis()->SetMoreLogLabels();
    histoXSecDummy[0]->GetXaxis()->SetLabelOffset(+0.01);
    histoXSecDummy[0]->Draw();

    //__________________________________________ Loop over all ratio pads and draw
    for(Int_t i=0;i<numActiveMeas;i++){
        padXSec[i+1]->cd();
        padXSec[i+1]->SetLogx(1);
        histoXSecDummy[i+1]                     = new TH2F(Form("histoXSecDummy_%d",i+1),Form("histoXSecDummy_%d",i+1),1000,xRangeMinXSec,xRangeMaxXSec,1000,0.6,1.95);
        SetStyleHistoTH2ForGraphs(histoXSecDummy[i+1], "#it{p}_{T} (GeV/#it{c})","#frac{Data}{TCM fit}", 0.85*textsizeLabelsXSec[i+1], textsizeLabelsXSec[i+1],
                                  0.85*textsizeLabelsXSec[i+1],textsizeLabelsXSec[i+1], 1,0.2/(textsizeFacXSec[i+1]*marginXSec), 510, 505);
        histoXSecDummy[i+1]->GetYaxis()->SetMoreLogLabels(kTRUE);
        histoXSecDummy[i+1]->GetYaxis()->SetNdivisions(505);
        histoXSecDummy[i+1]->GetYaxis()->SetNoExponent(kTRUE);
        histoXSecDummy[i+1]->GetXaxis()->SetMoreLogLabels(kTRUE);
        histoXSecDummy[i+1]->GetXaxis()->SetNoExponent(kTRUE);
        histoXSecDummy[i+1]->GetXaxis()->SetLabelFont(42);
        histoXSecDummy[i+1]->GetYaxis()->SetLabelFont(42);
        histoXSecDummy[i+1]->GetYaxis()->CenterTitle(kTRUE);
        histoXSecDummy[i+1]->GetYaxis()->SetLabelOffset(+0.01);
        histoXSecDummy[i+1]->GetXaxis()->SetTickLength(0.07);
        histoXSecDummy[i+1]->DrawCopy();
    }

    //__________________________________________ Draw combined spectra in first pad
    TF1* fitTsallisInvCrossSectionPi0Comb[6];
    TF1* fitTCMInvCrossSectionPi0Comb[6];
    TF1* fitTCMInvCrossSectionPi0CombPlot[6];
    padXSec[0]->cd();
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
                graphPi0InvariantCrossSectionSys[i][10]     ->Draw("E2same");
                graphPi0InvariantCrossSectionStat[i][10]    ->Draw("p,same,z");
                if(i== 1){
                  Double_t paramTCMPi0New[5]  = { 703678218.2483119965,0.5727060723,
                                                  68561118949.5456314087,0.4481417931,3.0869940568};
                  fitTCMInvCrossSectionPi0Comb[i]   = FitObject("tcm","fitTCMInvCrossSectionPi0Comb","Pi0",graphPi0InvariantCrossSectionStat[i][10],minXSpectra[i][10],maxXSpectra[i][10] ,paramTCMPi0New,"QNRMEX0+","", kFALSE);
                }else{
                  Double_t paramTCMComb[5]  = { graphPi0InvariantCrossSectionStat[i][10]->GetY()[1],0.1,graphPi0InvariantCrossSectionStat[i][10]->GetY()[4],0.6,3.0};
                  fitTCMInvCrossSectionPi0Comb[i]   = FitObject("tcm","fitTCMInvCrossSectionPi0Comb","Pi0",graphPi0InvariantCrossSectionStat[i][10],minXSpectra[i][10],maxXSpectra[i][10] ,paramTCMComb,"QNRMEX0+","", kFALSE);
                }
                fitTCMInvCrossSectionPi0CombPlot[i] = new TF1("twoCompModel_plotting7TeV",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mesonMassExpectPi0,mesonMassExpectPi0,mesonMassExpectPi0));
                fitTCMInvCrossSectionPi0CombPlot[i]->SetRange(minXSpectra[i][10],maxXSpectra[i][10]);
                fitTCMInvCrossSectionPi0CombPlot[i]->SetParameters(fitTCMInvCrossSectionPi0Comb[i]->GetParameters());
                fitTCMInvCrossSectionPi0CombPlot[i]->SetParErrors(fitTCMInvCrossSectionPi0Comb[i]->GetParErrors());
                cout << WriteParameterToFile(fitTCMInvCrossSectionPi0Comb[i]) << endl;
                DrawGammaSetMarkerTF1( fitTCMInvCrossSectionPi0CombPlot[i], 7, 2, kGray+2);
                fitTCMInvCrossSectionPi0CombPlot[i]->Draw("same");
                Double_t paramTsallisComb[3]    = {5e11, 6., 0.13};
                fitTsallisInvCrossSectionPi0Comb[i] = FitObject("l","fitInvCrossSectionPi08TeV","Pi0",graphPi0InvariantCrossSectionStat[i][10],minXSpectra[i][10],maxXSpectra[i][10],paramTsallisComb,"QNRMEX0+");
                DrawGammaSetMarkerTF1( fitTsallisInvCrossSectionPi0Comb[i], 3, 2, kGray+1);
                fitTsallisInvCrossSectionPi0Comb[i]->Draw("same");
            }
        }
    }
    Double_t rightalignDouble = 0.93;
    drawLatexAdd("ALICE",rightalignDouble,0.91,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",rightalignDouble,0.87,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

    TLegend* legendInvariantCrossSectionPi0    = GetAndSetLegend2(0.17, 0.05, 0.5, 0.05+textsizeLabelsXSec[0]*(numActiveMeas+1)+textsizeLabelsXSec[0], textSizeLabelsPixel);
    legendInvariantCrossSectionPi0->SetNColumns(1);
    legendInvariantCrossSectionPi0->SetMargin(0.2);
    TString legendScalingString[6] = {"", " (x 10)"," (x 10^{2})"," (x 10^{3})"," (x 10^{4})"," (x 10^{5})"};
    Int_t legendRunningIndex = numActiveMeas-1;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            legendInvariantCrossSectionPi0->AddEntry(graphPi0InvariantCrossSectionSys[i][10],Form("%s%s",energyLatex[i].Data(),legendScalingString[legendRunningIndex].Data()),"pf");
            legendRunningIndex-=1;
        }
    }
    legendInvariantCrossSectionPi0->AddEntry(fitTCMInvCrossSectionPi0CombPlot[exampleActiveMeas],"TCM fit","l");
    legendInvariantCrossSectionPi0->AddEntry(fitTsallisInvCrossSectionPi0Comb[exampleActiveMeas],"Tsallis fit","l");
    legendInvariantCrossSectionPi0->Draw();
    //__________________________________________ Draw ratios of individual spectra to combined fit in lower pads
    TGraphAsymmErrors* graphRatioIndivMeasToCombFitPi0Sys[6][11];
    TGraphAsymmErrors* graphRatioIndivMeasToCombFitPi0Stat[6][11];
    TGraphAsymmErrors* graphRatioIndivMeasToCombFitPi0Stat_woXErr[6][11];
    Int_t padCounter = 0;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            padXSec[padCounter+1]->cd();
            for (Int_t j = 0; j < 10; j++){
                if(graphPi0InvariantCrossSectionSys[i][j]&&graphPi0InvariantCrossSectionStat[i][j]){
                    //__________________________ Sys points
                    graphRatioIndivMeasToCombFitPi0Sys[i][j]     = (TGraphAsymmErrors*)graphPi0InvariantCrossSectionSys[i][j]->Clone();
                    graphRatioIndivMeasToCombFitPi0Sys[i][j]     = CalculateGraphErrRatioToFit(graphRatioIndivMeasToCombFitPi0Sys[i][j], fitTCMInvCrossSectionPi0CombPlot[i]);
                    DrawGammaSetMarkerTGraphAsym(graphRatioIndivMeasToCombFitPi0Sys[i][j], markerStyleMeas[j], markerSizeMeas[j], colorMeas[j], colorMeas[j], widthLinesBoxes, kTRUE, 0);
                    graphRatioIndivMeasToCombFitPi0Sys[i][j]->SetLineWidth(0);
                    graphRatioIndivMeasToCombFitPi0Sys[i][j]->Draw("2,same");
                    //__________________________ Stat points
                    graphRatioIndivMeasToCombFitPi0Stat[i][j]    = (TGraphAsymmErrors*)graphPi0InvariantCrossSectionStat[i][j]->Clone();
                    graphRatioIndivMeasToCombFitPi0Stat[i][j]    = CalculateGraphErrRatioToFit(graphRatioIndivMeasToCombFitPi0Stat[i][j], fitTCMInvCrossSectionPi0CombPlot[i]);
                    graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][j] = (TGraphAsymmErrors*) graphRatioIndivMeasToCombFitPi0Stat[i][j]->Clone("graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][j]");
                    ProduceGraphAsymmWithoutXErrors(graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][j]);
                    DrawGammaSetMarkerTGraphAsym(graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][j], markerStyleMeas[j], markerSizeMeas[j], colorMeas[j], colorMeas[j], widthLinesBoxes, kFALSE);
                    graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][j]->SetLineWidth(widthLinesBoxes);
                    graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][j]->Draw("p,same");
                }
            }
            DrawGammaLines(xRangeMinXSec, xRangeMaxXSec , 1., 1., 1.2, kGray+2, 7);
            padCounter+=1;
        }
    }
    canvasXSec->Print(Form("%s/Pi0_XSecFullPlot.%s",outputDir.Data(),suffix.Data()));

    TCanvas* canvasXSectionPi0  = new TCanvas("canvasXSectionPi0","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasXSectionPi0, 0.14, 0.02, 0.02, 0.09);
    canvasXSectionPi0->SetLogx();
    canvasXSectionPi0->SetLogy();
    histoXSecDummy[0]->Draw();
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
                graphPi0InvariantCrossSectionSys[i][10]     ->Draw("E2same");
                graphPi0InvariantCrossSectionStat[i][10]    ->Draw("p,same,z");
                fitTCMInvCrossSectionPi0CombPlot[i]->Draw("same");
                fitTsallisInvCrossSectionPi0Comb[i]->Draw("same");
            }
        }
    }
    drawLatexAdd("ALICE",rightalignDouble,0.91,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",rightalignDouble,0.87,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

    TLegend* legendInvariantCrossSectionPi02    = GetAndSetLegend2(0.17, 0.15, 0.5, 0.15+textsizeLabelsXSec[0]*(numActiveMeas+1)+textsizeLabelsXSec[0], textSizeLabelsPixel);
    legendInvariantCrossSectionPi02->SetNColumns(1);
    legendInvariantCrossSectionPi02->SetMargin(0.2);
    legendRunningIndex = numActiveMeas-1;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            legendInvariantCrossSectionPi02->AddEntry(graphPi0InvariantCrossSectionSys[i][10],Form("%s%s",energyLatex[i].Data(),legendScalingString[legendRunningIndex].Data()),"pf");
            legendRunningIndex-=1;
        }
    }
    legendInvariantCrossSectionPi02->AddEntry(fitTCMInvCrossSectionPi0CombPlot[exampleActiveMeas],"TCM fit","l");
    legendInvariantCrossSectionPi02->AddEntry(fitTsallisInvCrossSectionPi0Comb[exampleActiveMeas],"Tsallis fit","l");
    legendInvariantCrossSectionPi02->Draw();

    canvasXSectionPi0->Print(Form("%s/Pi0_XSec.%s",outputDir.Data(),suffix.Data()));



    //--------------------------------------------------------------------------------------------
    //__________________________________________ Define canvas, pads and sizes for ratio only plot
    //--------------------------------------------------------------------------------------------
    Int_t canvasDividerRatios                   = numActiveMeas+1;
    Double_t xRangeMinRatios                    = 0.23;
    Double_t xRangeMaxRatios                    = 49;
    Double_t arrayBoundariesX1_Ratios[2];
    Double_t arrayBoundariesY1_Ratios[canvasDividerRatios];
    Double_t relativeMarginsXRatios[3];
    Double_t relativeMarginsYRatios[3];
    Double_t textsizeLabelsRatios[canvasDividerRatios];
    Double_t textsizeFacRatios[canvasDividerRatios];
    ReturnCorrectValuesForCanvasScaling(1250,330*canvasDividerRatios, 1, canvasDividerRatios-1,0.135, 0.005, 0.005,0.072,arrayBoundariesX1_Ratios,arrayBoundariesY1_Ratios,relativeMarginsXRatios,relativeMarginsYRatios);
    Double_t marginRatios                         = relativeMarginsXRatios[0]*1250;
    TPad* padRatios[canvasDividerRatios];
    TH2F * histoRatiosDummy[canvasDividerRatios];

    TCanvas* canvasRatios                         = new TCanvas("canvasRatios","",0,0,1250,330*canvasDividerRatios);  // gives the page size
    DrawGammaCanvasSettings( canvasRatios,  0.13, 0.02, 0.03, 0.06);
    //__________________________________________ Set ratio pad properties
    for(Int_t i=0;i<numActiveMeas;i++){
        padRatios[i]                         = new TPad(Form("padRatiosFull_%d",i), "",
                                                      arrayBoundariesX1_Ratios[0], arrayBoundariesY1_Ratios[i+1],
                                                      arrayBoundariesX1_Ratios[1], arrayBoundariesY1_Ratios[i],-1, -1, -2);
        if(i==numActiveMeas-1)
            DrawGammaPadSettings( padRatios[i], relativeMarginsXRatios[0], relativeMarginsXRatios[2], relativeMarginsYRatios[1], relativeMarginsYRatios[2]);
        else
            DrawGammaPadSettings( padRatios[i], relativeMarginsXRatios[0], relativeMarginsXRatios[2], relativeMarginsYRatios[1], relativeMarginsYRatios[1]);
        padRatios[i]->Draw();
    }
   //__________________________________________ Set pad text sizes
    for(Int_t i=0;i<numActiveMeas;i++){
        if (padRatios[i]->XtoPixel(padRatios[i]->GetX2()) < padRatios[i]->YtoPixel(padRatios[i]->GetY1())){
            textsizeLabelsRatios[i]                  = (Double_t)textSizeLabelsPixel/padRatios[i]->XtoPixel(padRatios[i]->GetX2()) ;
            textsizeFacRatios[i]                   = (Double_t)1./padRatios[i]->XtoPixel(padRatios[i]->GetX2()) ;
        } else {
            textsizeLabelsRatios[i]                  = (Double_t)textSizeLabelsPixel/padRatios[i]->YtoPixel(padRatios[i]->GetY1());
            textsizeFacRatios[i]                   = (Double_t)1./padRatios[i]->YtoPixel(padRatios[i]->GetY1());
        }
    }

    //__________________________________________ Loop over all ratio pads and draw
    for(Int_t i=0;i<numActiveMeas;i++){
        padRatios[i]->cd();
        padRatios[i]->SetLogx(1);
        histoRatiosDummy[i]                     = new TH2F(Form("histoRatiosDummy_%d",i),Form("histoRatiosDummy_%d",i),1000,xRangeMinRatios,xRangeMaxRatios,1000,0.6,1.95);
        SetStyleHistoTH2ForGraphs(histoRatiosDummy[i], "#it{p}_{T} (GeV/#it{c})","#frac{Data}{TCM fit}", 0.85*textsizeLabelsRatios[i], textsizeLabelsRatios[i],
                                  0.85*textsizeLabelsRatios[i],textsizeLabelsRatios[i], 1,0.2/(textsizeFacRatios[i]*marginRatios), 510, 505);
        histoRatiosDummy[i]->GetYaxis()->SetMoreLogLabels(kTRUE);
        histoRatiosDummy[i]->GetYaxis()->SetNdivisions(505);
        histoRatiosDummy[i]->GetYaxis()->SetNoExponent(kTRUE);
        histoRatiosDummy[i]->GetXaxis()->SetMoreLogLabels(kTRUE);
        histoRatiosDummy[i]->GetXaxis()->SetNoExponent(kTRUE);
        histoRatiosDummy[i]->GetXaxis()->SetLabelFont(42);
        histoRatiosDummy[i]->GetYaxis()->SetLabelFont(42);
        histoRatiosDummy[i]->GetYaxis()->CenterTitle(kTRUE);
        histoRatiosDummy[i]->GetYaxis()->SetLabelOffset(+0.01);
        histoRatiosDummy[i]->GetXaxis()->SetTickLength(0.07);
        histoRatiosDummy[i]->DrawCopy();
    }
    padCounter = 0;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            padRatios[padCounter]->cd();
            for (Int_t j = 0; j < 10; j++){
                if(graphPi0InvariantCrossSectionSys[i][j]&&graphPi0InvariantCrossSectionStat[i][j]){
                    //__________________________ Sys points
                    graphRatioIndivMeasToCombFitPi0Sys[i][j]->Draw("2,same");
                    //__________________________ Stat points
                    graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][j]->Draw("p,same");
                }
            }
            DrawGammaLines(xRangeMinRatios, xRangeMaxRatios, 1., 1., 1.2, kGray+2, 7);
            padCounter+=1;
        }
    }
    canvasRatios->Print(Form("%s/Pi0_CombinationRatiosPP.%s",outputDir.Data(),suffix.Data()));

    //--------------------------------------------------------------------------------------------
    //__________________________________________ Plot spectra, ratios and theory
    //--------------------------------------------------------------------------------------------
    padXSec[0]->cd();
    histoXSecDummy[0]->Draw();
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
                fitTCMInvCrossSectionPi0CombPlot[i]->Draw("same");
                fitTsallisInvCrossSectionPi0Comb[i]->Draw("same");
                if(graphPythiaInvCrossSectionPi0[i])
                    graphPythiaInvCrossSectionPi0[i]->Draw("3,same");
                if(histoPythiaInvCrossSectionPi0[i])
                    histoPythiaInvCrossSectionPi0[i]->Draw("same,hist,l");
                graphPi0InvariantCrossSectionSys[i][10]     ->Draw("E2same");
                graphPi0InvariantCrossSectionStat[i][10]    ->Draw("p,same,z");
            }
        }
    }
    drawLatexAdd("ALICE",rightalignDouble,0.91,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",rightalignDouble,0.87,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

    TLegend* legendInvariantXSecTheoryPi0    = GetAndSetLegend2(0.17, 0.05, 0.5, 0.05+textsizeLabelsXSec[0]*numActiveMeas+2*textsizeLabelsXSec[0], textSizeLabelsPixel);
    legendInvariantXSecTheoryPi0->SetNColumns(1);
    legendInvariantXSecTheoryPi0->SetMargin(0.2);
    legendRunningIndex = numActiveMeas-1;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            legendInvariantXSecTheoryPi0->AddEntry(graphPi0InvariantCrossSectionSys[i][10],Form("%s%s",energyLatex[i].Data(),legendScalingString[legendRunningIndex].Data()),"pf");
            legendRunningIndex-=1;
        }
    }
    legendInvariantXSecTheoryPi0->AddEntry(fitTCMInvCrossSectionPi0CombPlot[exampleActiveMeas],"TCM fit","l");
    legendInvariantXSecTheoryPi0->AddEntry(fitTsallisInvCrossSectionPi0Comb[exampleActiveMeas],"Tsallis fit","l");
    legendInvariantXSecTheoryPi0->Draw();
    //__________________________________________ Loop over all ratio pads and draw
    for(Int_t i=0;i<numActiveMeas;i++){
        padXSec[i+1]->cd();
        histoXSecDummy[i+1]->DrawCopy();
    }
    padCounter = 0;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            padXSec[padCounter+1]->cd();
            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
                //__________________________ Sys points
                graphRatioIndivMeasToCombFitPi0Sys[i][10]     = (TGraphAsymmErrors*)graphPi0InvariantCrossSectionSys[i][10]->Clone();
                graphRatioIndivMeasToCombFitPi0Sys[i][10]     = CalculateGraphErrRatioToFit(graphRatioIndivMeasToCombFitPi0Sys[i][10], fitTCMInvCrossSectionPi0CombPlot[i]);
                DrawGammaSetMarkerTGraphAsym(graphRatioIndivMeasToCombFitPi0Sys[i][10], markerStyleEnergy[i], markerSizeEnergy[i], colorEnergy[i], colorEnergy[i], widthLinesBoxes, kTRUE, 0);
                graphRatioIndivMeasToCombFitPi0Sys[i][10]->SetLineWidth(0);
                graphRatioIndivMeasToCombFitPi0Sys[i][10]->Draw("2,same");
                //__________________________ Stat points
                graphRatioIndivMeasToCombFitPi0Stat[i][10]    = (TGraphAsymmErrors*)graphPi0InvariantCrossSectionStat[i][10]->Clone();
                graphRatioIndivMeasToCombFitPi0Stat[i][10]    = CalculateGraphErrRatioToFit(graphRatioIndivMeasToCombFitPi0Stat[i][10], fitTCMInvCrossSectionPi0CombPlot[i]);
                graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][10] = (TGraphAsymmErrors*) graphRatioIndivMeasToCombFitPi0Stat[i][10]->Clone("graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][10]");
                ProduceGraphAsymmWithoutXErrors(graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][10]);
                DrawGammaSetMarkerTGraphAsym(graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][10], markerStyleEnergy[i], markerSizeEnergy[i], colorEnergy[i], colorEnergy[i], widthLinesBoxes, kFALSE);
                graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][10]->SetLineWidth(widthLinesBoxes);
                graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][10]->Draw("p,same");
            }
            DrawGammaLines(xRangeMinXSec, xRangeMaxXSec , 1., 1., 1.2, kGray+2, 7);
            padCounter+=1;
        }
    }
    canvasXSec->Print(Form("%s/Pi0_XSecTheoryComparison.%s",outputDir.Data(),suffix.Data()));
    //--------------------------------------------------------------------------------------------
    padXSec[0]->cd();
    histoXSecDummy[0]->Draw();
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
                fitTCMInvCrossSectionPi0CombPlot[i]->Draw("same");
                fitTsallisInvCrossSectionPi0Comb[i]->Draw("same");
                graphPi0InvariantCrossSectionSys[i][10]     ->Draw("E2same");
                graphPi0InvariantCrossSectionStat[i][10]    ->Draw("p,same,z");
            }
        }
    }
    drawLatexAdd("ALICE",rightalignDouble,0.91,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",rightalignDouble,0.87,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

    legendInvariantXSecTheoryPi0->Draw();
    //__________________________________________ Loop over all ratio pads and draw
    for(Int_t i=0;i<numActiveMeas;i++){
        padXSec[i+1]->cd();
        histoXSecDummy[i+1]->DrawCopy();
    }
    padCounter = 0;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            padXSec[padCounter+1]->cd();
            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
                //__________________________ Sys points
                graphRatioIndivMeasToCombFitPi0Sys[i][10]     = (TGraphAsymmErrors*)graphPi0InvariantCrossSectionSys[i][10]->Clone();
                graphRatioIndivMeasToCombFitPi0Sys[i][10]     = CalculateGraphErrRatioToFit(graphRatioIndivMeasToCombFitPi0Sys[i][10], fitTCMInvCrossSectionPi0CombPlot[i]);
                DrawGammaSetMarkerTGraphAsym(graphRatioIndivMeasToCombFitPi0Sys[i][10], markerStyleEnergy[i], markerSizeEnergy[i], colorEnergy[i], colorEnergy[i], widthLinesBoxes, kTRUE, 0);
                graphRatioIndivMeasToCombFitPi0Sys[i][10]->SetLineWidth(0);
                graphRatioIndivMeasToCombFitPi0Sys[i][10]->Draw("2,same");
                //__________________________ Stat points
                graphRatioIndivMeasToCombFitPi0Stat[i][10]    = (TGraphAsymmErrors*)graphPi0InvariantCrossSectionStat[i][10]->Clone();
                graphRatioIndivMeasToCombFitPi0Stat[i][10]    = CalculateGraphErrRatioToFit(graphRatioIndivMeasToCombFitPi0Stat[i][10], fitTCMInvCrossSectionPi0CombPlot[i]);
                graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][10] = (TGraphAsymmErrors*) graphRatioIndivMeasToCombFitPi0Stat[i][10]->Clone("graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][10]");
                ProduceGraphAsymmWithoutXErrors(graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][10]);
                DrawGammaSetMarkerTGraphAsym(graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][10], markerStyleEnergy[i], markerSizeEnergy[i], colorEnergy[i], colorEnergy[i], widthLinesBoxes, kFALSE);
                graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][10]->SetLineWidth(widthLinesBoxes);
                graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][10]->Draw("p,same");
            }
            DrawGammaLines(xRangeMinXSec, xRangeMaxXSec , 1., 1., 1.2, kGray+2, 7);
            padCounter+=1;
        }
    }
    canvasXSec->Print(Form("%s/Pi0_XSecComparison.%s",outputDir.Data(),suffix.Data()));

    canvasRatios->cd();
    //__________________________________________ Set ratio pad properties
    for(Int_t i=0;i<numActiveMeas;i++){
        padRatios[i]->Draw();
    }

    //__________________________________________ Loop over all ratio pads and draw
    for(Int_t i=0;i<numActiveMeas;i++){
        padRatios[i]->cd();
        padRatios[i]->SetLogx(1);
        histoRatiosDummy[i]->DrawCopy();
    }
    padCounter = 0;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            padRatios[padCounter]->cd();

            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
                graphRatioIndivMeasToCombFitPi0Sys[i][10]->Draw("2,same");
                graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][10]->Draw("p,same");
            }
            DrawGammaLines(xRangeMinRatios, xRangeMaxRatios, 1., 1., 1.2, kGray+2, 7);
            padCounter+=1;
        }
    }
    canvasRatios->Print(Form("%s/Pi0_CombinationFullRatiosPP.%s",outputDir.Data(),suffix.Data()));







    //---------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    canvasXSec->cd();
    padXSec[0]->Draw();
    //__________________________________________ Set ratio pad properties
    for(Int_t i=0;i<numActiveMeas;i++){
        padXSec[i+1]->Draw();
    }

    //__________________________________________ Draw combined spectra in first pad
    TF1* fitTsallisInvCrossSectionEtaComb[6];
    TF1* fitTCMInvCrossSectionEtaComb[6];
    TF1* fitTCMInvCrossSectionEtaCombPlot[6];
    padXSec[0]->cd();
    histoXSecDummy[0]->Draw();
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){

            if(graphEtaInvariantCrossSectionSys[i][10]&&graphEtaInvariantCrossSectionStat[i][10]){
                graphEtaInvariantCrossSectionSys[i][10]     ->Draw("E2same");
                graphEtaInvariantCrossSectionStat[i][10]    ->Draw("p,same,z");
                if(i==3){
                  Double_t paramTCMComb[5]  = { graphEtaInvariantCrossSectionStat[i][10]->GetY()[1],0.1,graphEtaInvariantCrossSectionStat[i][10]->GetY()[2],0.6,3.0};
                  fitTCMInvCrossSectionEtaComb[i]   = FitObject("tcm","fitTCMInvCrossSectionEtaComb","Eta",graphEtaInvariantCrossSectionStat[i][10],minXSpectraEta[i][10],maxXSpectraEta[i][10] ,paramTCMComb,"QNRMEX0+","", kFALSE);
                }else{
                  Double_t paramTCMComb[5]  = { graphEtaInvariantCrossSectionStat[i][10]->GetY()[1],0.1,graphEtaInvariantCrossSectionStat[i][10]->GetY()[4],0.6,3.0};
                  fitTCMInvCrossSectionEtaComb[i]   = FitObject("tcm","fitTCMInvCrossSectionEtaComb","Eta",graphEtaInvariantCrossSectionStat[i][10],minXSpectraEta[i][10],maxXSpectraEta[i][10] ,paramTCMComb,"QNRMEX0+","", kFALSE);
                }
                fitTCMInvCrossSectionEtaCombPlot[i] = new TF1("twoCompModel_plotting7TeV",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mesonMassExpectEta,mesonMassExpectEta,mesonMassExpectEta));
                fitTCMInvCrossSectionEtaCombPlot[i]->SetRange(minXSpectraEta[i][10],maxXSpectraEta[i][10]);
                fitTCMInvCrossSectionEtaCombPlot[i]->SetParameters(fitTCMInvCrossSectionEtaComb[i]->GetParameters());
                fitTCMInvCrossSectionEtaCombPlot[i]->SetParErrors(fitTCMInvCrossSectionEtaComb[i]->GetParErrors());
                cout << WriteParameterToFile(fitTCMInvCrossSectionEtaComb[i]) << endl;
                DrawGammaSetMarkerTF1( fitTCMInvCrossSectionEtaCombPlot[i], 7, 2, kGray+2);
                if(i!=0) fitTCMInvCrossSectionEtaCombPlot[i]->Draw("same");
                Double_t paramTsallisComb[3]    = {graphEtaInvariantCrossSectionStat[i][10]->GetY()[1], 6.6, 0.22};
                fitTsallisInvCrossSectionEtaComb[i] = FitObject("l","fitInvCrossSectionEta8TeV","Eta",graphEtaInvariantCrossSectionStat[i][10],minXSpectraEta[i][10],maxXSpectraEta[i][10],paramTsallisComb,"QNRMEX0+");
                DrawGammaSetMarkerTF1( fitTsallisInvCrossSectionEtaComb[i], 3, 2, kGray+1);
                if(i!=0) fitTsallisInvCrossSectionEtaComb[i]->Draw("same");
            }
        }
    }
    drawLatexAdd("ALICE",rightalignDouble,0.91,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#eta #rightarrow #gamma#gamma",rightalignDouble,0.87,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

    TLegend* legendInvariantCrossSectionEta    = GetAndSetLegend2(0.17, 0.05, 0.5, 0.05+textsizeLabelsXSec[0]*(numActiveMeas+1)+textsizeLabelsXSec[0], textSizeLabelsPixel);
    legendInvariantCrossSectionEta->SetNColumns(1);
    legendInvariantCrossSectionEta->SetMargin(0.2);
    legendRunningIndex = numActiveMeas-1;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            legendInvariantCrossSectionEta->AddEntry(graphEtaInvariantCrossSectionSys[i][10],Form("%s%s",energyLatex[i].Data(),legendScalingString[legendRunningIndex].Data()),"pf");
            legendRunningIndex-=1;
        }
    }
    legendInvariantCrossSectionEta->AddEntry(fitTCMInvCrossSectionEtaCombPlot[exampleActiveMeas],"TCM fit","l");
    legendInvariantCrossSectionEta->AddEntry(fitTsallisInvCrossSectionEtaComb[exampleActiveMeas],"Tsallis fit","l");
    legendInvariantCrossSectionEta->Draw();
    //__________________________________________ Draw ratios of individual spectra to combined fit in lower pads
    TGraphAsymmErrors* graphRatioIndivMeasToCombFitEtaSys[6][11];
    TGraphAsymmErrors* graphRatioIndivMeasToCombFitEtaStat[6][11];
    TGraphAsymmErrors* graphRatioIndivMeasToCombFitEtaStat_woXErr[6][11];
    padCounter = 0;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            padXSec[padCounter+1]->cd();
            histoXSecDummy[padCounter+1]->DrawCopy();
            for (Int_t j = 0; j < 10; j++){
                if(graphEtaInvariantCrossSectionSys[i][j]&&graphEtaInvariantCrossSectionStat[i][j]){
                    //__________________________ Sys points
//                    graphRatioIndivMeasToCombFitEtaSys[i][j]     = (TGraphAsymmErrors*)graphEtaInvariantCrossSectionSys[i][j]->Clone();
//                    graphRatioIndivMeasToCombFitEtaSys[i][j]     = CalculateGraphErrRatioToFit(graphRatioIndivMeasToCombFitEtaSys[i][j], fitTCMInvCrossSectionEtaCombPlot[i]);
//                    DrawGammaSetMarkerTGraphAsym(graphRatioIndivMeasToCombFitEtaSys[i][j], markerStyleMeas[j], markerSizeMeas[j], colorMeas[j], colorMeas[j], widthLinesBoxes, kTRUE, 0);
//                    graphRatioIndivMeasToCombFitEtaSys[i][j]->SetLineWidth(0);
//                    graphRatioIndivMeasToCombFitEtaSys[i][j]->Draw("2,same");
//                    //__________________________ Stat points
//                    graphRatioIndivMeasToCombFitEtaStat[i][j]    = (TGraphAsymmErrors*)graphEtaInvariantCrossSectionStat[i][j]->Clone();
//                    graphRatioIndivMeasToCombFitEtaStat[i][j]    = CalculateGraphErrRatioToFit(graphRatioIndivMeasToCombFitEtaStat[i][j], fitTCMInvCrossSectionEtaCombPlot[i]);
//                    graphRatioIndivMeasToCombFitEtaStat_woXErr[i][j] = (TGraphAsymmErrors*) graphRatioIndivMeasToCombFitEtaStat[i][j]->Clone("graphRatioIndivMeasToCombFitEtaStat_woXErr[i][j]");
//                    ProduceGraphAsymmWithoutXErrors(graphRatioIndivMeasToCombFitEtaStat_woXErr[i][j]);
//                    DrawGammaSetMarkerTGraphAsym(graphRatioIndivMeasToCombFitEtaStat_woXErr[i][j], markerStyleMeas[j], markerSizeMeas[j], colorMeas[j], colorMeas[j], widthLinesBoxes, kFALSE);
//                    graphRatioIndivMeasToCombFitEtaStat_woXErr[i][j]->SetLineWidth(widthLinesBoxes);
//                    graphRatioIndivMeasToCombFitEtaStat_woXErr[i][j]->Draw("p,same");
                }
            }
            DrawGammaLines(xRangeMinXSec, xRangeMaxXSec , 1., 1., 1.2, kGray+2, 7);
            padCounter+=1;
        }
    }
    canvasXSec->Print(Form("%s/Eta_XSecFullPlot.%s",outputDir.Data(),suffix.Data()));


    TCanvas* canvasXSectionEta  = new TCanvas("canvasXSectionEta","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasXSectionEta, 0.14, 0.02, 0.02, 0.09);
    canvasXSectionEta->SetLogx();
    canvasXSectionEta->SetLogy();
    histoXSecDummy[0]->GetYaxis()->SetRangeUser(1E3,9E12);
    histoXSecDummy[0]->Draw();
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            if(graphEtaInvariantCrossSectionSys[i][10]&&graphEtaInvariantCrossSectionStat[i][10]){
                graphEtaInvariantCrossSectionSys[i][10]     ->Draw("E2same");
                graphEtaInvariantCrossSectionStat[i][10]    ->Draw("p,same,z");
                if(i!=0) fitTCMInvCrossSectionEtaCombPlot[i]->Draw("same");
                if(i!=0) fitTsallisInvCrossSectionEtaComb[i]->Draw("same");
            }
        }
    }
    drawLatexAdd("ALICE",rightalignDouble,0.91,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#eta #rightarrow #gamma#gamma",rightalignDouble,0.87,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

    TLegend* legendInvariantCrossSectionEta2    = GetAndSetLegend2(0.17, 0.15, 0.5, 0.15+textsizeLabelsXSec[0]*(numActiveMeas+1)+textsizeLabelsXSec[0], textSizeLabelsPixel);
    legendInvariantCrossSectionEta2->SetNColumns(1);
    legendInvariantCrossSectionEta2->SetMargin(0.2);
    legendRunningIndex = numActiveMeas-1;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            legendInvariantCrossSectionEta2->AddEntry(graphEtaInvariantCrossSectionSys[i][10],Form("%s%s",energyLatex[i].Data(),legendScalingString[legendRunningIndex].Data()),"pf");
            legendRunningIndex-=1;
        }
    }
    legendInvariantCrossSectionEta2->AddEntry(fitTCMInvCrossSectionEtaCombPlot[exampleActiveMeas],"TCM fit","l");
    legendInvariantCrossSectionEta2->AddEntry(fitTsallisInvCrossSectionEtaComb[exampleActiveMeas],"Tsallis fit","l");
    legendInvariantCrossSectionEta2->Draw();

    canvasXSectionEta->Print(Form("%s/Eta_XSec.%s",outputDir.Data(),suffix.Data()));
}
