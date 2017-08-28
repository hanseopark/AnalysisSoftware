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
#include "TGraph.h"
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

    //---------------------------------------------------------------------------------------------------------------
    //---------------------------- General setting & global variables -----------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    TString fileName[6];
        // /*  900GeV  */  fileName[0]                 = "CombinationInputPP/900GeV/CombinedResultsPaperPP900GeV_2017_07_18.root";
    /*  900GeV  */  fileName[0]                 = "CombinationInputPP/900GeV/CombinedResultsPaper7TeVand900GeV_IncludingPP2760YShiftedPrelim_Pub2012.root";
    /*  2.76TeV */  fileName[1]                 = "CombinationInputPP/2.76TeV/CombinedResultsPaperPP2760GeV_2017_07_10_FrediV2Clusterizer.root";
    /*  5TeV    */  fileName[2]                 = "CombinationInputPP/5TeV/";
        // /*  7TeV    */  fileName[3]                 = "CombinationInputPP/7TeV/CombinedResultsPaperPP7TeV_2017_07_18.root";
    /*  7TeV    */  fileName[3]                 = "CombinationInputPP/7TeV/CombinedResultsPaper7TeVand900GeV_IncludingPP2760YShiftedPrelim_Pub2012.root";
    /*  8TeV    */  fileName[4]                 = "CombinationInputPP/8TeV/CombinedResultsPaperPP8TeV_2017_07_13.root";
    /*  13TeV   */  fileName[5]                 = "CombinationInputPP/13TeV/";
    TString         fileNameTheory              = "ExternalInput/Theory/TheoryCompilationPP.root";
    TString         suffix                      = "eps";

    Int_t includeEnergy[6]                      = {1,1,0,1,1,0};
    Int_t numActiveMeas                         = 4;

    gROOT->Reset();
    gROOT->SetStyle("Plain");

    StyleSettingsThesis();
    SetPlotStyle();

    gStyle->SetEndErrorSize(0);

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

    fstream fLog;
    fLog.open(Form("%s/CombineMesonPP.log",outputDir.Data()), ios::out);
    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fLog << dateForOutput.Data() << endl;

    //__________________________________________ Loading input files and directories
    for(Int_t i=0; i<6; i++){
        if(includeEnergy[i]){
            inputFile[i]                        = new TFile(fileName[i].Data());
              if(inputFile[i])     cout << "Found " << nameEnergyGlobal[i] << " input"  << endl;
            directoryPi0[i]                     = (TDirectory*)inputFile[i]->Get(Form("Pi0%s",nameEnergyGlobal[i].Data()));
              if(!directoryPi0[i]) directoryPi0[i] = NULL;
              else                 cout << "--> Found directory " << Form("Pi0%s",nameEnergyGlobal[i].Data()) <<  endl;
            directoryEta[i]                     = (TDirectory*)inputFile[i]->Get(Form("Eta%s",nameEnergyGlobal[i].Data()));
              if(!directoryEta[i]) directoryEta[i] = NULL;
              else                 cout << "--> Found directory " << Form("Eta%s",nameEnergyGlobal[i].Data()) << endl;
            exampleActiveMeas                   = i;
        }
    }
    TFile* fileTheoryCompilation                            = new TFile(fileNameTheory.Data());

    //__________________________________________ Definition of colors, styles and markers sizes
    LoadColorsMarkersAndSizes();

    //---------------------------------------------------------------------------------------------------------------
    //---------------------------- Loading graphs from files --------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    TGraphAsymmErrors* graphPi0InvariantCrossSectionStat[6][11];
    TGraphAsymmErrors* graphPi0InvariantCrossSectionSys[6][11];
    TGraphAsymmErrors* graphPi0InvariantCrossSectionTot[6][11];
    TGraphAsymmErrors* graphEtaInvariantCrossSectionStat[6][11];
    TGraphAsymmErrors* graphEtaInvariantCrossSectionSys[6][11];
    TGraphAsymmErrors* graphEtaInvariantCrossSectionTot[6][11];

    TH1D* histoEtaToPi0RatioStat[6][11];
    TGraphAsymmErrors* graphEtaToPi0RatioStat[6][11];
    TGraphAsymmErrors* graphEtaToPi0RatioSys[6][11];
    TGraphAsymmErrors* graphEtaToPi0RatioTot[6][11];

    Double_t minXSpectra[6][11];
    Double_t maxXSpectra[6][11];
    Double_t minXSpectraEta[6][11];
    Double_t maxXSpectraEta[6][11];

    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            for (Int_t j = 0; j < 11; j++){
              // skip PCM-PHOS for the moment
              if(j==3){
                graphPi0InvariantCrossSectionStat[i][j] = NULL;
                graphPi0InvariantCrossSectionSys[i][j]  = NULL;
                graphPi0InvariantCrossSectionTot[i][j]  = NULL;
                graphEtaInvariantCrossSectionStat[i][j] = NULL;
                graphEtaInvariantCrossSectionSys[i][j]  = NULL;
                graphEtaInvariantCrossSectionTot[i][j]  = NULL;
                continue;
              }

                //______________________________ Loading pi0 inv. cross sections
                if(directoryPi0[i]) {
                  graphPi0InvariantCrossSectionStat[i][j]   = (TGraphAsymmErrors*)directoryPi0[i]->Get(Form("graphInvCrossSectionPi0%s%s%sStatErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data(),graphNameModifier[j].Data()));
                  graphPi0InvariantCrossSectionSys[i][j]    = (TGraphAsymmErrors*)directoryPi0[i]->Get(Form("graphInvCrossSectionPi0%s%s%sSysErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data(),graphNameModifier[j].Data()));
                }else{
                  graphPi0InvariantCrossSectionStat[i][j]   = (TGraphAsymmErrors*)inputFile[i]->Get(Form("graphInvCrossSectionPi0%sStat%s%s",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data(),graphNameModifier[j].Data()));
                  graphPi0InvariantCrossSectionSys[i][j]    = (TGraphAsymmErrors*)inputFile[i]->Get(Form("graphInvCrossSectionPi0%sSys%s%s",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data(),graphNameModifier[j].Data()));
                  if(!graphPi0InvariantCrossSectionStat[i][j]&&!graphPi0InvariantCrossSectionSys[i][j]){
                    graphPi0InvariantCrossSectionStat[i][j]   = (TGraphAsymmErrors*)inputFile[i]->Get(Form("graphInvCrossSectionPi0%s%sStatErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data()));
                    graphPi0InvariantCrossSectionSys[i][j]    = (TGraphAsymmErrors*)inputFile[i]->Get(Form("graphInvCrossSectionPi0%s%sSysErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data()));
                  }
                }
                if(graphPi0InvariantCrossSectionStat[i][j]&&graphPi0InvariantCrossSectionSys[i][j]){
                    graphPi0InvariantCrossSectionTot[i][j] = (TGraphAsymmErrors*) graphPi0InvariantCrossSectionStat[i][j]->Clone(Form("graphInvCrossSectionPi0%s%s%sTotErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data(),graphNameModifier[j].Data()));
                    for(Int_t iE=0; iE<graphPi0InvariantCrossSectionTot[i][j]->GetN();iE++){
                      graphPi0InvariantCrossSectionTot[i][j]->GetEYhigh()[iE] = TMath::Sqrt( TMath::Power(graphPi0InvariantCrossSectionStat[i][j]->GetErrorYhigh(iE),2) + TMath::Power(graphPi0InvariantCrossSectionSys[i][j]->GetErrorYhigh(iE),2));
                      graphPi0InvariantCrossSectionTot[i][j]->GetEYlow()[iE] = TMath::Sqrt( TMath::Power(graphPi0InvariantCrossSectionStat[i][j]->GetErrorYlow(iE),2) + TMath::Power(graphPi0InvariantCrossSectionSys[i][j]->GetErrorYlow(iE),2));
                    }
                    cout << "found " << Form("graphInvCrossSectionPi0%s%sStatErr+SysErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data()) << endl;
                    //__________________________ Reading global xRange for each spectrum
                    minXSpectra[i][j]           = graphPi0InvariantCrossSectionStat[i][j]->GetX()[0]-graphPi0InvariantCrossSectionStat[i][j]->GetEXlow()[0];
                    maxXSpectra[i][j]           = graphPi0InvariantCrossSectionStat[i][j]->GetX()[graphPi0InvariantCrossSectionStat[i][j]->GetN()-1]+graphPi0InvariantCrossSectionStat[i][j]->GetEXhigh()[graphPi0InvariantCrossSectionStat[i][j]->GetN()-1];
                    if(j==10) cout << "--> minx: " << minXSpectra[i][j] << "\tmaxX: " << maxXSpectra[i][j] << endl;
                }else{
                  graphPi0InvariantCrossSectionStat[i][j]=NULL;
                  graphPi0InvariantCrossSectionSys[i][j]=NULL;
                  graphPi0InvariantCrossSectionTot[i][j]=NULL;
                }

//                if(i==0 && j==10){
//                  nameMeasGlobal[j]                  = "PCM";
//                  nameMeasGlobalshort[j]             = "PCM";
//                  graphNameModifier[j]               = "";
//                }

                //______________________________ Loading eta inv. cross sections
                if(directoryEta[i]) {
                  graphEtaInvariantCrossSectionStat[i][j]   = (TGraphAsymmErrors*)directoryEta[i]->Get(Form("graphInvCrossSectionEta%s%s%sStatErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data(),graphNameModifier[j].Data()));
                  graphEtaInvariantCrossSectionSys[i][j]    = (TGraphAsymmErrors*)directoryEta[i]->Get(Form("graphInvCrossSectionEta%s%s%sSysErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data(),graphNameModifier[j].Data()));
                }else{
                  graphEtaInvariantCrossSectionStat[i][j]   = (TGraphAsymmErrors*)inputFile[i]->Get(Form("graphInvCrossSectionEta%sStat%s%s",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data(),graphNameModifier[j].Data()));
                  graphEtaInvariantCrossSectionSys[i][j]    = (TGraphAsymmErrors*)inputFile[i]->Get(Form("graphInvCrossSectionEta%sSys%s%s",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data(),graphNameModifier[j].Data()));
                  if(!graphEtaInvariantCrossSectionStat[i][j]&&!graphEtaInvariantCrossSectionSys[i][j]){
                    graphEtaInvariantCrossSectionStat[i][j]   = (TGraphAsymmErrors*)inputFile[i]->Get(Form("graphInvCrossSectionEta%s%sStatErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data()));
                    graphEtaInvariantCrossSectionSys[i][j]    = (TGraphAsymmErrors*)inputFile[i]->Get(Form("graphInvCrossSectionEta%s%sSysErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data()));
                    if(!graphEtaInvariantCrossSectionStat[i][j]&&!graphEtaInvariantCrossSectionSys[i][j]){
                      graphEtaInvariantCrossSectionStat[i][j]   = (TGraphAsymmErrors*)inputFile[i]->Get(Form("graphInvCrossSectionEta%s%sStatErr_PrelimQM2011",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data()));
                      graphEtaInvariantCrossSectionSys[i][j]    = (TGraphAsymmErrors*)inputFile[i]->Get(Form("graphInvCrossSectionEta%s%sSysErr_PrelimQM2011",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data()));
                      if(i==0 && j==0){
                        graphEtaInvariantCrossSectionStat[i][j]   = (TGraphAsymmErrors*)inputFile[i]->Get(Form("graphInvCrossSectionEtaComb%sStatErr_PrelimQM2011",nameEnergyGlobal2[i].Data()));
                        graphEtaInvariantCrossSectionSys[i][j]    = (TGraphAsymmErrors*)inputFile[i]->Get(Form("graphInvCrossSectionEtaComb%sSysErr_PrelimQM2011",nameEnergyGlobal2[i].Data()));
                        graphEtaInvariantCrossSectionStat[i][j]->SetName(Form("graphInvCrossSectionEta%s%sStatErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data()));
                        graphEtaInvariantCrossSectionSys[i][j]->SetName(Form("graphInvCrossSectionEta%s%sSysErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data()));
                      }
                    }
                  }
                }
                if(graphEtaInvariantCrossSectionStat[i][j]&&graphEtaInvariantCrossSectionSys[i][j]){
                    graphEtaInvariantCrossSectionTot[i][j] = (TGraphAsymmErrors*) graphEtaInvariantCrossSectionStat[i][j]->Clone(Form("graphInvCrossSectionEta%s%s%sTotErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data(),graphNameModifier[j].Data()));
                    for(Int_t iE=0; iE<graphEtaInvariantCrossSectionTot[i][j]->GetN();iE++){
                      graphEtaInvariantCrossSectionTot[i][j]->GetEYhigh()[iE] = TMath::Sqrt( TMath::Power(graphEtaInvariantCrossSectionStat[i][j]->GetErrorYhigh(iE),2) + TMath::Power(graphEtaInvariantCrossSectionSys[i][j]->GetErrorYhigh(iE),2));
                      graphEtaInvariantCrossSectionTot[i][j]->GetEYlow()[iE] = TMath::Sqrt( TMath::Power(graphEtaInvariantCrossSectionStat[i][j]->GetErrorYlow(iE),2) + TMath::Power(graphEtaInvariantCrossSectionSys[i][j]->GetErrorYlow(iE),2));
                    }
                    cout << "found " << Form("graphInvCrossSectionEta%s%sStatErr+SysErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data()) << endl;
                    minXSpectraEta[i][j]           = graphEtaInvariantCrossSectionStat[i][j]->GetX()[0]-graphEtaInvariantCrossSectionStat[i][j]->GetEXlow()[0];
                    maxXSpectraEta[i][j]           = graphEtaInvariantCrossSectionStat[i][j]->GetX()[graphEtaInvariantCrossSectionStat[i][j]->GetN()-1]+graphEtaInvariantCrossSectionStat[i][j]->GetEXhigh()[graphEtaInvariantCrossSectionStat[i][j]->GetN()-1];
                    if(j==10) cout << "--> minx: " << minXSpectraEta[i][j] << "\tmaxX: " << maxXSpectraEta[i][j] << endl;
                }else{
                  graphEtaInvariantCrossSectionStat[i][j]=NULL;
                  graphEtaInvariantCrossSectionSys[i][j]=NULL;
                  graphEtaInvariantCrossSectionTot[i][j]=NULL;
                }
                //______________________________ Loading eta/pi0 ratios
                if(directoryEta[i]) {
                  graphEtaToPi0RatioStat[i][j]   = (TGraphAsymmErrors*)directoryEta[i]->Get(Form("graphRatioEtaToPi0%s%sStatErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data()));
                  histoEtaToPi0RatioStat[i][j]   = (TH1D*)directoryEta[i]->Get(Form("histoRatioEtaToPi0%s%sStatErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data()));
                  if(histoEtaToPi0RatioStat[i][j]) graphEtaToPi0RatioStat[i][j]= new TGraphAsymmErrors(histoEtaToPi0RatioStat[i][j]);
                  graphEtaToPi0RatioSys[i][j]   = (TGraphAsymmErrors*)directoryEta[i]->Get(Form("graphRatioEtaToPi0%s%sSysErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data()));
                }else{
                  graphEtaToPi0RatioStat[i][j]   = (TGraphAsymmErrors*)inputFile[i]->Get(Form("graphEtaToPi0%s%sStat",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data()));
                  histoEtaToPi0RatioStat[i][j]   = (TH1D*)inputFile[i]->Get(Form("histoEtaToPi0%s%sStat",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data()));
                  if(histoEtaToPi0RatioStat[i][j]) graphEtaToPi0RatioStat[i][j]= new TGraphAsymmErrors(histoEtaToPi0RatioStat[i][j]);
                  graphEtaToPi0RatioSys[i][j]   = (TGraphAsymmErrors*)inputFile[i]->Get(Form("graphEtaToPi0%s%sSys",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data()));
                }
                if(graphEtaToPi0RatioStat[i][j] && graphEtaToPi0RatioSys[i][j])
                    cout << "found " << Form("graphRatioEtaToPi0%s%sStatErr+SysErr",nameMeasGlobalshort[j].Data(),nameEnergyGlobal2[i].Data()) << endl;

//                if(i==0 && j==10){
//                  nameMeasGlobal[j]                  = "Comb";
//                  nameMeasGlobalshort[j]             = "Comb";
//                  graphNameModifier[j]               = "A";
//                }
            }
        }
    }

    //---------------------------------------------------------------------------------------------------------------
    //---------------------------- Loading pythia predictions from files --------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

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
                histoPythiaInvCrossSectionPi0[i]->Scale(scalingFactorPythia);
                DrawGammaSetMarker(histoPythiaInvCrossSectionPi0[i], 24, 1.5, pythia8color , pythia8color);
                histoPythiaInvCrossSectionPi0[i]->SetLineWidth(widthCommonFit);
                histoPythiaInvCrossSectionPi0[i]->GetXaxis()->SetRangeUser(minXSpectra[i][10]+0.05,maxXSpectra[i][10]);
                graphPythiaInvCrossSectionPi0[i]= new TGraphErrors((TH1F*) fileTheoryCompilation->Get(pythiaHistoNames[i].Data()));
                for (int j=0;j<graphPythiaInvCrossSectionPi0[i]->GetN();j++){
                    graphPythiaInvCrossSectionPi0[i]->GetY()[j] *= scalingFactorPythia;
                    graphPythiaInvCrossSectionPi0[i]->GetEY()[j] *= scalingFactorPythia;
                }
                while(graphPythiaInvCrossSectionPi0[i]->GetX()[0] < minXSpectra[i][10]+0.05) graphPythiaInvCrossSectionPi0[i]->RemovePoint(0);
                while(graphPythiaInvCrossSectionPi0[i]->GetX()[graphPythiaInvCrossSectionPi0[i]->GetN()-1] > maxXSpectra[i][10]) graphPythiaInvCrossSectionPi0[i]->RemovePoint(graphPythiaInvCrossSectionPi0[i]->GetN()-1);
                DrawGammaSetMarkerTGraphErr(graphPythiaInvCrossSectionPi0[i], 0, 0, pythia8color , pythia8color, widthLinesBoxes, kTRUE, pythia8color);
                ProduceGraphErrWithoutXErrors(graphPythiaInvCrossSectionPi0[i]);
            }else{
                histoPythiaInvCrossSectionPi0[i]=NULL;
                graphPythiaInvCrossSectionPi0[i]=NULL;
            }
            histoPythiaInvCrossSectionEta[i]    = (TH1F*) fileTheoryCompilation->Get(pythiaHistoNamesEta[i].Data());
            if(histoPythiaInvCrossSectionEta[i]){
                histoPythiaInvCrossSectionEta[i]->GetXaxis()->SetRangeUser(minXSpectraEta[i][10]+0.05,maxXSpectraEta[i][10]);
                histoPythiaInvCrossSectionEta[i]->Scale(scalingFactorPythia);
                DrawGammaSetMarker(histoPythiaInvCrossSectionEta[i], 24, 1.5, pythia8color , pythia8color);
                histoPythiaInvCrossSectionEta[i]->SetLineWidth(widthCommonFit);
                histoPythiaInvCrossSectionEta[i]->GetXaxis()->SetRangeUser(minXSpectraEta[i][10]+0.05,maxXSpectraEta[i][10]);
                graphPythiaInvCrossSectionEta[i]= new TGraphErrors((TH1F*) fileTheoryCompilation->Get(pythiaHistoNamesEta[i].Data()));
                for (int j=0;j<graphPythiaInvCrossSectionEta[i]->GetN();j++){
                    graphPythiaInvCrossSectionEta[i]->GetY()[j] *= scalingFactorPythia;
                    graphPythiaInvCrossSectionEta[i]->GetEY()[j] *= scalingFactorPythia;
                }
                while(graphPythiaInvCrossSectionEta[i]->GetX()[0] < minXSpectraEta[i][10]+0.05) graphPythiaInvCrossSectionEta[i]->RemovePoint(0);
                while(graphPythiaInvCrossSectionEta[i]->GetX()[graphPythiaInvCrossSectionEta[i]->GetN()-1] > maxXSpectraEta[i][10]) graphPythiaInvCrossSectionEta[i]->RemovePoint(graphPythiaInvCrossSectionEta[i]->GetN()-1);
                DrawGammaSetMarkerTGraphErr(graphPythiaInvCrossSectionEta[i], 0, 0, pythia8color , pythia8color, widthLinesBoxes, kTRUE, pythia8color);
                ProduceGraphErrWithoutXErrors(graphPythiaInvCrossSectionEta[i]);
            }else{
                histoPythiaInvCrossSectionEta[i]=NULL;
                graphPythiaInvCrossSectionEta[i]=NULL;
            }
            scalingFactorPythia                 *=10;
        }
    }

    //---------------------------------------------------------------------------------------------------------------
    //---------------------------- Loading NLO predictions from files -----------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    TString NLO07HistoNamesHalf[6]            = { "graphNLOCalcInvSecPi0MuHalf900GeV", "graphNLOCalcInvSecPi0MuHalf2760GeV", "", "graphNLOCalcInvSecPi0MuHalf7000GeV", "graphNLOCalcInvSecPi0MuHalf8000GeV", ""};
    TString NLO07HistoNamesOne[6]             = { "graphNLOCalcInvSecPi0MuOne900GeV", "graphNLOCalcInvSecPi0MuOne2760GeV", "", "graphNLOCalcInvSecPi0MuOne7000GeV", "graphNLOCalcInvSecPi0MuOne8000GeV", ""};
    TString NLO07HistoNamesTwo[6]             = { "graphNLOCalcInvSecPi0MuTwo900GeV", "graphNLOCalcInvSecPi0MuTwo2760GeV", "", "graphNLOCalcInvSecPi0MuTwo7000GeV", "graphNLOCalcInvSecPi0MuTwo8000GeV", ""};
    TGraph* graphNLODSS07InvCrossSectionPi0Half[6];
    TGraph* graphNLODSS07InvCrossSectionPi0One[6];
    TGraph* graphNLODSS07InvCrossSectionPi0Two[6];
    TString NLO14HistoNames[6]                 = { "", "graphNLOCalcDSS14InvCrossSec2760GeV", "", "graphNLOCalcDSS14InvCrossSec7000GeV", "graphNLOCalcDSS14InvSecPi08000GeV", ""};
    TGraphAsymmErrors* graphNLODSS14InvCrossSectionPi0[6];
    TString NLOAESSSHistoNamesEtaHalf[6]      = { "graphNLOCalcInvSecEtaMuHalf900GeV", "graphNLOCalcInvSecEtaMuHalf2760GeV", "", "graphNLOCalcInvSecEtaMuHalf7000GeV", "graphNLOCalcInvSecEtaMuHalf8000GeV", ""};
    TString NLOAESSSHistoNamesEtaOne[6]       = { "graphNLOCalcInvSecEtaMuOne900GeV", "graphNLOCalcInvSecEtaMuOne2760GeV", "", "graphNLOCalcInvSecEtaMuOne7000GeV", "graphNLOCalcInvSecEtaMuOne8000GeV", ""};
    TString NLOAESSSHistoNamesEtaTwo[6]       = { "graphNLOCalcInvSecEtaMuTwo900GeV", "graphNLOCalcInvSecEtaMuTwo2760GeV", "", "graphNLOCalcInvSecEtaMuTwo7000GeV", "graphNLOCalcInvSecEtaMuTwo8000GeV", ""};
    TGraph* graphNLOAESSSInvCrossSectionEtaHalf[6];
    TGraph* graphNLOAESSSInvCrossSectionEtaOne[6];
    TGraph* graphNLOAESSSInvCrossSectionEtaTwo[6];

    Double_t scalingFactorNLO = 1;

    for(Int_t i=0;i<6;i++){
        if(includeEnergy[i]){
            if(fileTheoryCompilation->Get(NLO07HistoNamesHalf[i].Data()))
              graphNLODSS07InvCrossSectionPi0Half[i] = (TGraph*) ((TGraph*) fileTheoryCompilation->Get(NLO07HistoNamesHalf[i].Data()))->Clone(Form("graphNLODSS07InvCrossSectionPi0Half_%i",i));
            else
              graphNLODSS07InvCrossSectionPi0Half[i] = NULL;

            if(graphNLODSS07InvCrossSectionPi0Half[i]){
                for (int j=0;j<graphNLODSS07InvCrossSectionPi0Half[i]->GetN();j++){
                    graphNLODSS07InvCrossSectionPi0Half[i]->GetY()[j] *= scalingFactorNLO;
                }
                while(i!=0 && graphNLODSS07InvCrossSectionPi0Half[i]->GetX()[0] < minXSpectra[i][10]) graphNLODSS07InvCrossSectionPi0Half[i]->RemovePoint(0);
                while(i!=0 && graphNLODSS07InvCrossSectionPi0Half[i]->GetX()[graphNLODSS07InvCrossSectionPi0Half[i]->GetN()-1] > maxXSpectra[i][10]) graphNLODSS07InvCrossSectionPi0Half[i]->RemovePoint(graphNLODSS07InvCrossSectionPi0Half[i]->GetN()-1);
                DrawGammaNLOTGraph( graphNLODSS07InvCrossSectionPi0Half[i], widthCommonFit, 8, colorNLO);
            }else{
                graphNLODSS07InvCrossSectionPi0Half[i]=NULL;
            }

            if(fileTheoryCompilation->Get(NLO07HistoNamesOne[i].Data()))
              graphNLODSS07InvCrossSectionPi0One[i] = (TGraph*) ((TGraph*) fileTheoryCompilation->Get(NLO07HistoNamesOne[i].Data()))->Clone(Form("graphNLODSS07InvCrossSectionPi0One_%i",i));
            else
              graphNLODSS07InvCrossSectionPi0One[i] = NULL;

            if(graphNLODSS07InvCrossSectionPi0One[i]){
                for (int j=0;j<graphNLODSS07InvCrossSectionPi0One[i]->GetN();j++){
                    graphNLODSS07InvCrossSectionPi0One[i]->GetY()[j] *= scalingFactorNLO;
                }
                while(i!=0 && graphNLODSS07InvCrossSectionPi0One[i]->GetX()[0] < minXSpectra[i][10]) graphNLODSS07InvCrossSectionPi0One[i]->RemovePoint(0);
                while(i!=0 && graphNLODSS07InvCrossSectionPi0One[i]->GetX()[graphNLODSS07InvCrossSectionPi0One[i]->GetN()-1] > maxXSpectra[i][10]) graphNLODSS07InvCrossSectionPi0One[i]->RemovePoint(graphNLODSS07InvCrossSectionPi0One[i]->GetN()-1);
                DrawGammaNLOTGraph( graphNLODSS07InvCrossSectionPi0One[i], widthCommonFit, 7, colorNLO);
            }else{
                graphNLODSS07InvCrossSectionPi0One[i]=NULL;
            }

            if(fileTheoryCompilation->Get(NLO07HistoNamesTwo[i].Data()))
              graphNLODSS07InvCrossSectionPi0Two[i] = (TGraph*) ((TGraph*) fileTheoryCompilation->Get(NLO07HistoNamesTwo[i].Data()))->Clone(Form("graphNLODSS07InvCrossSectionPi0Two_%i",i));
            else
              graphNLODSS07InvCrossSectionPi0Two[i] = NULL;

            if(graphNLODSS07InvCrossSectionPi0Two[i]){
                for (int j=0;j<graphNLODSS07InvCrossSectionPi0Two[i]->GetN();j++){
                    graphNLODSS07InvCrossSectionPi0Two[i]->GetY()[j] *= scalingFactorNLO;
                }
                while(i!=0 && graphNLODSS07InvCrossSectionPi0Two[i]->GetX()[0] < minXSpectra[i][10]) graphNLODSS07InvCrossSectionPi0Two[i]->RemovePoint(0);
                while(i!=0 && graphNLODSS07InvCrossSectionPi0Two[i]->GetX()[graphNLODSS07InvCrossSectionPi0Two[i]->GetN()-1] > maxXSpectra[i][10]) graphNLODSS07InvCrossSectionPi0Two[i]->RemovePoint(graphNLODSS07InvCrossSectionPi0Two[i]->GetN()-1);
                DrawGammaNLOTGraph( graphNLODSS07InvCrossSectionPi0Two[i], widthCommonFit, 4, colorNLO);
            }else{
                graphNLODSS07InvCrossSectionPi0Two[i]=NULL;
            }

            if(fileTheoryCompilation->Get(NLO14HistoNames[i].Data()))
              graphNLODSS14InvCrossSectionPi0[i] = (TGraphAsymmErrors*) ((TGraphAsymmErrors*) fileTheoryCompilation->Get(NLO14HistoNames[i].Data()))->Clone(Form("graphNLODSS14InvCrossSectionPi0_%i",i));
            else
              graphNLODSS14InvCrossSectionPi0[i] = NULL;

            if(graphNLODSS14InvCrossSectionPi0[i]){
                for (int j=0;j<graphNLODSS14InvCrossSectionPi0[i]->GetN();j++){
                    graphNLODSS14InvCrossSectionPi0[i]->GetY()[j] *= scalingFactorNLO;
                    graphNLODSS14InvCrossSectionPi0[i]->GetEYhigh()[j] *= scalingFactorNLO;
                    graphNLODSS14InvCrossSectionPi0[i]->GetEYlow()[j] *= scalingFactorNLO;
                }
                while(graphNLODSS14InvCrossSectionPi0[i]->GetX()[0] < minXSpectra[i][10]) graphNLODSS14InvCrossSectionPi0[i]->RemovePoint(0);
                while(graphNLODSS14InvCrossSectionPi0[i]->GetX()[graphNLODSS14InvCrossSectionPi0[i]->GetN()-1] > maxXSpectra[i][10]) graphNLODSS14InvCrossSectionPi0[i]->RemovePoint(graphNLODSS14InvCrossSectionPi0[i]->GetN()-1);
                DrawGammaSetMarkerTGraphAsym(graphNLODSS14InvCrossSectionPi0[i], 0, 0, colorNLO , colorNLO, widthLinesBoxes, kTRUE, colorNLO);
                ProduceGraphAsymmWithoutXErrors(graphNLODSS14InvCrossSectionPi0[i]);
                graphNLODSS14InvCrossSectionPi0[i]->SetLineWidth(widthCommonFit);
                graphNLODSS14InvCrossSectionPi0[i]->SetLineColor(colorNLO);
                graphNLODSS14InvCrossSectionPi0[i]->SetLineStyle(1);
                graphNLODSS14InvCrossSectionPi0[i]->SetFillStyle(1001);
                graphNLODSS14InvCrossSectionPi0[i]->SetFillColor(colorNLO);
            }else{
                graphNLODSS14InvCrossSectionPi0[i]=NULL;
            }
            if(fileTheoryCompilation->Get(NLOAESSSHistoNamesEtaHalf[i].Data()))
              graphNLOAESSSInvCrossSectionEtaHalf[i] = (TGraph*) ((TGraph*) fileTheoryCompilation->Get(NLOAESSSHistoNamesEtaHalf[i].Data()))->Clone(Form("graphNLOAESSSInvCrossSectionEtaHalf_%i",i));
            else
              graphNLOAESSSInvCrossSectionEtaHalf[i] = NULL;

            if(graphNLOAESSSInvCrossSectionEtaHalf[i]){
                for (int j=0;j<graphNLOAESSSInvCrossSectionEtaHalf[i]->GetN();j++){
                    graphNLOAESSSInvCrossSectionEtaHalf[i]->GetY()[j] *= scalingFactorNLO;
                }
                while(graphNLOAESSSInvCrossSectionEtaHalf[i]->GetX()[0] < minXSpectraEta[i][10]) graphNLOAESSSInvCrossSectionEtaHalf[i]->RemovePoint(0);
                if(i==0) while(graphNLOAESSSInvCrossSectionEtaHalf[i]->GetX()[graphNLOAESSSInvCrossSectionEtaHalf[i]->GetN()-1] > maxXSpectraEta[i][10]+2.) graphNLOAESSSInvCrossSectionEtaHalf[i]->RemovePoint(graphNLOAESSSInvCrossSectionEtaHalf[i]->GetN()-1);
                else while(graphNLOAESSSInvCrossSectionEtaHalf[i]->GetX()[graphNLOAESSSInvCrossSectionEtaHalf[i]->GetN()-1] > maxXSpectraEta[i][10]) graphNLOAESSSInvCrossSectionEtaHalf[i]->RemovePoint(graphNLOAESSSInvCrossSectionEtaHalf[i]->GetN()-1);
                DrawGammaNLOTGraph( graphNLOAESSSInvCrossSectionEtaHalf[i], widthCommonFit, 8, colorNLO);
            }else{
                graphNLOAESSSInvCrossSectionEtaHalf[i]=NULL;
            }

            if(fileTheoryCompilation->Get(NLOAESSSHistoNamesEtaOne[i].Data()))
              graphNLOAESSSInvCrossSectionEtaOne[i] = (TGraph*) ((TGraph*) fileTheoryCompilation->Get(NLOAESSSHistoNamesEtaOne[i].Data()))->Clone(Form("graphNLOAESSSInvCrossSectionEtaOne_%i",i));
            else
              graphNLOAESSSInvCrossSectionEtaOne[i] = NULL;

            if(graphNLOAESSSInvCrossSectionEtaOne[i]){
                for (int j=0;j<graphNLOAESSSInvCrossSectionEtaOne[i]->GetN();j++){
                    graphNLOAESSSInvCrossSectionEtaOne[i]->GetY()[j] *= scalingFactorNLO;
                }
                while(graphNLOAESSSInvCrossSectionEtaOne[i]->GetX()[0] < minXSpectraEta[i][10]) graphNLOAESSSInvCrossSectionEtaOne[i]->RemovePoint(0);
                if(i==0) while(graphNLOAESSSInvCrossSectionEtaOne[i]->GetX()[graphNLOAESSSInvCrossSectionEtaOne[i]->GetN()-1] > maxXSpectraEta[i][10]+2.) graphNLOAESSSInvCrossSectionEtaOne[i]->RemovePoint(graphNLOAESSSInvCrossSectionEtaOne[i]->GetN()-1);
                else while(graphNLOAESSSInvCrossSectionEtaOne[i]->GetX()[graphNLOAESSSInvCrossSectionEtaOne[i]->GetN()-1] > maxXSpectraEta[i][10]) graphNLOAESSSInvCrossSectionEtaOne[i]->RemovePoint(graphNLOAESSSInvCrossSectionEtaOne[i]->GetN()-1);
                DrawGammaNLOTGraph( graphNLOAESSSInvCrossSectionEtaOne[i], widthCommonFit, 7, colorNLO);
            }else{
                graphNLOAESSSInvCrossSectionEtaOne[i]=NULL;
            }

            if(fileTheoryCompilation->Get(NLOAESSSHistoNamesEtaTwo[i].Data()))
              graphNLOAESSSInvCrossSectionEtaTwo[i] = (TGraph*) ((TGraph*) fileTheoryCompilation->Get(NLOAESSSHistoNamesEtaTwo[i].Data()))->Clone(Form("graphNLOAESSSInvCrossSectionEtaTwo_%i",i));
            else
              graphNLOAESSSInvCrossSectionEtaTwo[i] = NULL;

            if(graphNLOAESSSInvCrossSectionEtaTwo[i]){
                for (int j=0;j<graphNLOAESSSInvCrossSectionEtaTwo[i]->GetN();j++){
                    graphNLOAESSSInvCrossSectionEtaTwo[i]->GetY()[j] *= scalingFactorNLO;
                }
                while(graphNLOAESSSInvCrossSectionEtaTwo[i]->GetX()[0] < minXSpectraEta[i][10]) graphNLOAESSSInvCrossSectionEtaTwo[i]->RemovePoint(0);
                if(i==0) while(graphNLOAESSSInvCrossSectionEtaTwo[i]->GetX()[graphNLOAESSSInvCrossSectionEtaTwo[i]->GetN()-1] > maxXSpectraEta[i][10]+2.) graphNLOAESSSInvCrossSectionEtaTwo[i]->RemovePoint(graphNLOAESSSInvCrossSectionEtaTwo[i]->GetN()-1);
                else while(graphNLOAESSSInvCrossSectionEtaTwo[i]->GetX()[graphNLOAESSSInvCrossSectionEtaTwo[i]->GetN()-1] > maxXSpectraEta[i][10]) graphNLOAESSSInvCrossSectionEtaTwo[i]->RemovePoint(graphNLOAESSSInvCrossSectionEtaTwo[i]->GetN()-1);
                DrawGammaNLOTGraph( graphNLOAESSSInvCrossSectionEtaTwo[i], widthCommonFit, 4, colorNLO);
            }else{
                graphNLOAESSSInvCrossSectionEtaTwo[i]=NULL;
            }

            scalingFactorNLO                 *=10;
        }
    }
    //---------------------------------------------------------------------------------------------------------------
    //---------------------------- Scaling cross sections for plotting ----------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

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
                    ProduceGraphAsymmWithoutXErrors(graphPi0InvariantCrossSectionStat[i][j]);
                    DrawGammaSetMarkerTGraph(graphPi0InvariantCrossSectionStat[i][j], markerStyleMeas[j], markerSizeMeas[j], colorMeas[j] , colorMeas[j]);
                    if(j==10)
                        DrawGammaSetMarkerTGraph(graphPi0InvariantCrossSectionStat[i][j], markerStyleEnergy[i], markerSizeEnergy[i], colorEnergy[i] , colorEnergy[i]);
                }
                if(graphPi0InvariantCrossSectionTot[i][j]){
                    for (int k=0;k<graphPi0InvariantCrossSectionTot[i][j]->GetN();k++){
                        graphPi0InvariantCrossSectionTot[i][j]->GetY()[k]      *= scalingFactorXsec;
                        graphPi0InvariantCrossSectionTot[i][j]->GetEYhigh()[k] *= scalingFactorXsec;
                        graphPi0InvariantCrossSectionTot[i][j]->GetEYlow()[k]  *= scalingFactorXsec;
                    }
                    ProduceGraphAsymmWithoutXErrors(graphPi0InvariantCrossSectionTot[i][j]);
                    DrawGammaSetMarkerTGraph(graphPi0InvariantCrossSectionTot[i][j], markerStyleMeas[j], markerSizeMeas[j], colorMeas[j] , colorMeas[j]);
                    if(j==10)
                        DrawGammaSetMarkerTGraph(graphPi0InvariantCrossSectionTot[i][j], markerStyleEnergy[i], markerSizeEnergy[i], colorEnergy[i] , colorEnergy[i]);
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
                    ProduceGraphAsymmWithoutXErrors(graphEtaInvariantCrossSectionStat[i][j]);
                    DrawGammaSetMarkerTGraph(graphEtaInvariantCrossSectionStat[i][j], markerStyleMeas[j], markerSizeMeas[j], colorMeas[j] , colorMeas[j]);
                    if(j==10)
                        DrawGammaSetMarkerTGraph(graphEtaInvariantCrossSectionStat[i][j], markerStyleEnergy[i], markerSizeEnergy[i], colorEnergy[i] , colorEnergy[i]);
                }
                if(graphEtaInvariantCrossSectionTot[i][j]){
                    for (int k=0;k<graphEtaInvariantCrossSectionTot[i][j]->GetN();k++){
                        graphEtaInvariantCrossSectionTot[i][j]->GetY()[k]      *= scalingFactorXsec;
                        graphEtaInvariantCrossSectionTot[i][j]->GetEYhigh()[k] *= scalingFactorXsec;
                        graphEtaInvariantCrossSectionTot[i][j]->GetEYlow()[k]  *= scalingFactorXsec;
                    }
                    ProduceGraphAsymmWithoutXErrors(graphEtaInvariantCrossSectionTot[i][j]);
                    DrawGammaSetMarkerTGraph(graphEtaInvariantCrossSectionTot[i][j], markerStyleMeas[j], markerSizeMeas[j], colorMeas[j] , colorMeas[j]);
                    if(j==10)
                        DrawGammaSetMarkerTGraph(graphEtaInvariantCrossSectionTot[i][j], markerStyleEnergy[i], markerSizeEnergy[i], colorEnergy[i] , colorEnergy[i]);
                }
            }
            scalingFactorXsec                   *= 10;
        }
    }

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------TCM fits -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    TF1* fitTCMInvCrossSectionPi0Comb[6];
    TF1* fitTCMInvCrossSectionPi0CombPlot[6];
    TF1* fitTCMInvCrossSectionEtaComb[6];
    TF1* fitTCMInvCrossSectionEtaCombPlot[6];
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            if(graphPi0InvariantCrossSectionTot[i][10]){
                fLog << "TCM - " << nameEnergyGlobal2[i].Data() << endl;
                fLog << endl;
                if(i == 1){
                  Double_t paramTCMPi0New[5]  = { 703678218.2483119965,0.5727060723,
                                                  68561118949.5456314087,0.4481417931,3.0869940568};
                  fitTCMInvCrossSectionPi0Comb[i]   = FitObject("tcm",Form("fitTCMInvCrossSectionPi0Comb_%s",nameEnergyGlobal2[i].Data()),"Pi0",graphPi0InvariantCrossSectionTot[i][10],minXSpectra[i][10],maxXSpectra[i][10] ,paramTCMPi0New,"QNRMEX0+","", kFALSE);
                }else{
                  Double_t paramTCMComb[5]  = { graphPi0InvariantCrossSectionTot[i][10]->GetY()[1],0.1,graphPi0InvariantCrossSectionTot[i][10]->GetY()[3],0.6,3.0};
                  fitTCMInvCrossSectionPi0Comb[i]   = FitObject("tcm",Form("fitTCMInvCrossSectionPi0Comb_%s",nameEnergyGlobal2[i].Data()),"Pi0",graphPi0InvariantCrossSectionTot[i][10],minXSpectra[i][10],maxXSpectra[i][10] ,paramTCMComb,"QNRMEX0+","", kFALSE);
                }
                fitTCMInvCrossSectionPi0CombPlot[i] = new TF1(Form("twoCompModel_plotting_%s",nameEnergyGlobal2[i].Data()),Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mesonMassExpectPi0,mesonMassExpectPi0,mesonMassExpectPi0));
                fitTCMInvCrossSectionPi0CombPlot[i]->SetRange(minXSpectra[i][10],maxXSpectra[i][10]);
                fitTCMInvCrossSectionPi0CombPlot[i]->SetParameters(fitTCMInvCrossSectionPi0Comb[i]->GetParameters());
                fitTCMInvCrossSectionPi0CombPlot[i]->SetParErrors(fitTCMInvCrossSectionPi0Comb[i]->GetParErrors());
                DrawGammaSetMarkerTF1( fitTCMInvCrossSectionPi0CombPlot[i], 7, 2, kGray+2);

                fLog << WriteParameterToFile(fitTCMInvCrossSectionPi0Comb[i]) << endl;
            }


            if(graphEtaInvariantCrossSectionTot[i][10]){
                fLog << "TCM - " << nameEnergyGlobal2[i].Data() << endl;
                fLog << endl;
                if(i==3){
                  Double_t paramTCMComb[5]  = { graphEtaInvariantCrossSectionTot[i][10]->GetY()[1],0.1,graphEtaInvariantCrossSectionTot[i][10]->GetY()[2],0.6,3.0};
                  fitTCMInvCrossSectionEtaComb[i]   = FitObject("tcm",Form("fitTCMInvCrossSectionEtaComb_%s",nameEnergyGlobal2[i].Data()),"Eta",graphEtaInvariantCrossSectionTot[i][10],minXSpectraEta[i][10],maxXSpectraEta[i][10] ,paramTCMComb,"QNRMEX0+","", kFALSE);
                }else{
                  Double_t paramTCMComb[5]  = { graphEtaInvariantCrossSectionTot[i][10]->GetY()[1],0.1,graphEtaInvariantCrossSectionTot[i][10]->GetY()[4],0.6,3.0};
                  fitTCMInvCrossSectionEtaComb[i]   = FitObject("tcm",Form("fitTCMInvCrossSectionEtaComb_%s",nameEnergyGlobal2[i].Data()),"Eta",graphEtaInvariantCrossSectionTot[i][10],minXSpectraEta[i][10],maxXSpectraEta[i][10] ,paramTCMComb,"QNRMEX0+","", kFALSE);
                }
                fitTCMInvCrossSectionEtaCombPlot[i] = new TF1(Form("twoCompModel_Eta_plotting_%s",nameEnergyGlobal2[i].Data()),Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mesonMassExpectEta,mesonMassExpectEta,mesonMassExpectEta));
                fitTCMInvCrossSectionEtaCombPlot[i]->SetRange(minXSpectraEta[i][10],maxXSpectraEta[i][10]);
                fitTCMInvCrossSectionEtaCombPlot[i]->SetParameters(fitTCMInvCrossSectionEtaComb[i]->GetParameters());
                fitTCMInvCrossSectionEtaCombPlot[i]->SetParErrors(fitTCMInvCrossSectionEtaComb[i]->GetParErrors());
                DrawGammaSetMarkerTF1( fitTCMInvCrossSectionEtaCombPlot[i], 7, 2, kGray+2);

                fLog << WriteParameterToFile(fitTCMInvCrossSectionEtaComb[i]) << endl;
            }
        }
    }

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Tsallis fits -------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    TF1* fitTsallisInvCrossSectionPi0Comb[6];
    TF1* fitTsallisInvCrossSectionEtaComb[6];
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            if(graphPi0InvariantCrossSectionTot[i][10]){
                fLog << "Tsallis - " << nameEnergyGlobal2[i].Data() << endl;
                fLog << endl;
                Double_t paramTsallisComb[3]    = {5e11, 6., 0.13};
                fitTsallisInvCrossSectionPi0Comb[i] = FitObject("l",Form("fitInvCrossSectionPi08TeV_%s",nameEnergyGlobal2[i].Data()),"Pi0",graphPi0InvariantCrossSectionTot[i][10],minXSpectra[i][10],maxXSpectra[i][10],paramTsallisComb,"QNRMEX0+");
                DrawGammaSetMarkerTF1( fitTsallisInvCrossSectionPi0Comb[i], 3, 2, kGray+1);

                fLog << WriteParameterToFile(fitTsallisInvCrossSectionPi0Comb[i]) << endl;
            }

            if(graphEtaInvariantCrossSectionTot[i][10]){
              fLog << "Tsallis - Eta -" << nameEnergyGlobal2[i].Data() << endl;
              fLog << endl;
              if(i==0){
                fitTsallisInvCrossSectionEtaComb[i] = new TF1(Form("fitInvCrossSectionEta8TeV_%s",nameEnergyGlobal2[i].Data()),Form("[0] / ( 2 * TMath::Pi())*([1]-1.)*([1]-2.) / ([1]*[2]*([1]*[2]+%.10f*([1]-2.)))  * pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1])",mesonMassExpectEta,mesonMassExpectEta,mesonMassExpectEta,mesonMassExpectEta));
                fitTsallisInvCrossSectionEtaComb[i]->SetRange(minXSpectraEta[i][10],maxXSpectraEta[i][10]);
                fitTsallisInvCrossSectionEtaComb[i]->SetParameter(0,graphEtaInvariantCrossSectionTot[i][10]->GetY()[1]);
                fitTsallisInvCrossSectionEtaComb[i]->FixParameter(1,7.9);
                fitTsallisInvCrossSectionEtaComb[i]->SetParameter(2,0.18);
                graphEtaInvariantCrossSectionTot[i][10]->Fit(fitTsallisInvCrossSectionEtaComb[i],"QNRMEX0+","",minXSpectraEta[i][10],maxXSpectraEta[i][10]);
              }else{
                Double_t paramTsallisComb[3]    = {graphEtaInvariantCrossSectionTot[i][10]->GetY()[1], 6.6, 0.22};
                fitTsallisInvCrossSectionEtaComb[i] = FitObject("l",Form("fitInvCrossSectionEta8TeV_%s",nameEnergyGlobal2[i].Data()),"Eta",graphEtaInvariantCrossSectionTot[i][10],minXSpectraEta[i][10],maxXSpectraEta[i][10],paramTsallisComb,"QNRMEX0+");
              }
              DrawGammaSetMarkerTF1( fitTsallisInvCrossSectionEtaComb[i], 3, 2, kGray+1);

              fLog << WriteParameterToFile(fitTsallisInvCrossSectionEtaComb[i]) << endl;
            }
        }
    }

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Data to Fit Ratios -------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    TGraphAsymmErrors* graphRatioIndivMeasToCombFitPi0Sys[6][11];
    TGraphAsymmErrors* graphRatioIndivMeasToCombFitPi0Stat[6][11];
    TGraphAsymmErrors* graphRatioIndivMeasToCombFitPi0Stat_woXErr[6][11];
    TGraphAsymmErrors* graphRatioIndivMeasToCombFitEtaSys[6][11];
    TGraphAsymmErrors* graphRatioIndivMeasToCombFitEtaStat[6][11];
    TGraphAsymmErrors* graphRatioIndivMeasToCombFitEtaStat_woXErr[6][11];
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            for (Int_t j = 0; j < 11; j++){
                if(graphPi0InvariantCrossSectionSys[i][j]&&graphPi0InvariantCrossSectionStat[i][j]){
                    //__________________________ Sys points
                    graphRatioIndivMeasToCombFitPi0Sys[i][j]     = (TGraphAsymmErrors*)graphPi0InvariantCrossSectionSys[i][j]->Clone(Form("graphRatioIndivMeasToCombFitPi0Sys%i_%i",i,j));
                    graphRatioIndivMeasToCombFitPi0Sys[i][j]     = CalculateGraphErrRatioToFit(graphRatioIndivMeasToCombFitPi0Sys[i][j], fitTCMInvCrossSectionPi0CombPlot[i]);
                    if(j==10) DrawGammaSetMarkerTGraphAsym(graphRatioIndivMeasToCombFitPi0Sys[i][j], markerStyleEnergy[i], markerSizeEnergy[i], colorEnergy[i], colorEnergy[i], widthLinesBoxes, kTRUE, 0);
                    else      DrawGammaSetMarkerTGraphAsym(graphRatioIndivMeasToCombFitPi0Sys[i][j], markerStyleMeas[j], markerSizeMeas[j], colorMeas[j], colorMeas[j], widthLinesBoxes, kTRUE, 0);
                    graphRatioIndivMeasToCombFitPi0Sys[i][j]->SetLineWidth(0);
                    //__________________________ Stat points
                    graphRatioIndivMeasToCombFitPi0Stat[i][j]    = (TGraphAsymmErrors*)graphPi0InvariantCrossSectionStat[i][j]->Clone(Form("graphRatioIndivMeasToCombFitPi0Stat%i_%i",i,j));
                    graphRatioIndivMeasToCombFitPi0Stat[i][j]    = CalculateGraphErrRatioToFit(graphRatioIndivMeasToCombFitPi0Stat[i][j], fitTCMInvCrossSectionPi0CombPlot[i]);
                    graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][j] = (TGraphAsymmErrors*) graphRatioIndivMeasToCombFitPi0Stat[i][j]->Clone(Form("graphRatioIndivMeasToCombFitPi0Stat_woXErr%i_%i",i,j));
                    ProduceGraphAsymmWithoutXErrors(graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][j]);
                    if(j==10) DrawGammaSetMarkerTGraphAsym(graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][j], markerStyleEnergy[i], markerSizeEnergy[i], colorEnergy[i], colorEnergy[i], widthLinesBoxes, kTRUE, 0);
                    else      DrawGammaSetMarkerTGraphAsym(graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][j], markerStyleMeas[j], markerSizeMeas[j], colorMeas[j], colorMeas[j], widthLinesBoxes, kFALSE);
                    graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][j]->SetLineWidth(widthLinesBoxes);
                }

                if(graphEtaInvariantCrossSectionSys[i][j]&&graphEtaInvariantCrossSectionStat[i][j]){
                    //__________________________ Sys points
                    graphRatioIndivMeasToCombFitEtaSys[i][j]     = (TGraphAsymmErrors*)graphEtaInvariantCrossSectionSys[i][j]->Clone(Form("graphRatioIndivMeasToCombFitEtaSys%i_%i",i,j));
                    if(i==0&&(j==0||j==10)) graphRatioIndivMeasToCombFitEtaSys[i][j]     = CalculateGraphErrRatioToFit(graphRatioIndivMeasToCombFitEtaSys[i][j], fitTsallisInvCrossSectionEtaComb[i]);
                    else graphRatioIndivMeasToCombFitEtaSys[i][j]     = CalculateGraphErrRatioToFit(graphRatioIndivMeasToCombFitEtaSys[i][j], fitTCMInvCrossSectionEtaCombPlot[i]);
                    if(j==10) DrawGammaSetMarkerTGraphAsym(graphRatioIndivMeasToCombFitEtaSys[i][j], markerStyleEnergy[i], markerSizeEnergy[i], colorEnergy[i], colorEnergy[i], widthLinesBoxes, kTRUE, 0);
                    else      DrawGammaSetMarkerTGraphAsym(graphRatioIndivMeasToCombFitEtaSys[i][j], markerStyleMeas[j], markerSizeMeas[j], colorMeas[j], colorMeas[j], widthLinesBoxes, kTRUE, 0);
                    graphRatioIndivMeasToCombFitEtaSys[i][j]->SetLineWidth(0);
                    //__________________________ Stat points
                    graphRatioIndivMeasToCombFitEtaStat[i][j]    = (TGraphAsymmErrors*)graphEtaInvariantCrossSectionStat[i][j]->Clone(Form("graphRatioIndivMeasToCombFitEtaStat%i_%i",i,j));
                    if(i==0&&(j==0||j==10)) graphRatioIndivMeasToCombFitEtaStat[i][j]     = CalculateGraphErrRatioToFit(graphRatioIndivMeasToCombFitEtaStat[i][j], fitTsallisInvCrossSectionEtaComb[i]);
                    else graphRatioIndivMeasToCombFitEtaStat[i][j]    = CalculateGraphErrRatioToFit(graphRatioIndivMeasToCombFitEtaStat[i][j], fitTCMInvCrossSectionEtaCombPlot[i]);
                    graphRatioIndivMeasToCombFitEtaStat_woXErr[i][j] = (TGraphAsymmErrors*) graphRatioIndivMeasToCombFitEtaStat[i][j]->Clone(Form("graphRatioIndivMeasToCombFitEtaStat_woXErr%i_%i",i,j));
                    ProduceGraphAsymmWithoutXErrors(graphRatioIndivMeasToCombFitEtaStat_woXErr[i][j]);
                    if(j==10) DrawGammaSetMarkerTGraphAsym(graphRatioIndivMeasToCombFitEtaStat_woXErr[i][j], markerStyleEnergy[i], markerSizeEnergy[i], colorEnergy[i], colorEnergy[i], widthLinesBoxes, kTRUE, 0);
                    else      DrawGammaSetMarkerTGraphAsym(graphRatioIndivMeasToCombFitEtaStat_woXErr[i][j], markerStyleMeas[j], markerSizeMeas[j], colorMeas[j], colorMeas[j], widthLinesBoxes, kFALSE);
                    graphRatioIndivMeasToCombFitEtaStat_woXErr[i][j]->SetLineWidth(widthLinesBoxes);
                }
            }
        }
    }

    //---------------------------------------------------------------------------------------------------------------
    //------------------------------- Pythia to Fit Ratios ----------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    TGraphErrors* graphRatioPythiaToCombFitPi0[6];
    TH1D*         histoRatioPythiaToCombFitPi0[6];
    TGraphErrors* graphRatioPythiaToCombFitEta[6];
    TH1D*         histoRatioPythiaToCombFitEta[6];
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
          if(graphPythiaInvCrossSectionPi0[i]){
            graphRatioPythiaToCombFitPi0[i]     = (TGraphErrors*)graphPythiaInvCrossSectionPi0[i]->Clone(Form("graphRatioPythiaToCombFitPi0%i",i));
            graphRatioPythiaToCombFitPi0[i]     = CalculateGraphErrRatioToFit(graphRatioPythiaToCombFitPi0[i], fitTCMInvCrossSectionPi0CombPlot[i]);
            DrawGammaSetMarkerTGraphErr(graphRatioPythiaToCombFitPi0[i], 0, 0, pythia8color , pythia8color, widthLinesBoxes, kTRUE, pythia8color);
          }
          if(histoPythiaInvCrossSectionPi0[i]){
            histoRatioPythiaToCombFitPi0[i]     = (TH1D*)histoPythiaInvCrossSectionPi0[i]->Clone(Form("histoRatioPythiaToCombFitPi0%i",i));
            histoRatioPythiaToCombFitPi0[i]     = CalculateHistoRatioToFit(histoRatioPythiaToCombFitPi0[i], fitTCMInvCrossSectionPi0CombPlot[i]);
            DrawGammaSetMarker(histoRatioPythiaToCombFitPi0[i], 24, 1.5, pythia8color , pythia8color);
            histoRatioPythiaToCombFitPi0[i]->SetLineWidth(widthCommonFit);
          }

          if(graphPythiaInvCrossSectionEta[i]){
            graphRatioPythiaToCombFitEta[i]     = (TGraphErrors*)graphPythiaInvCrossSectionEta[i]->Clone(Form("graphRatioPythiaToCombFitEta%i",i));
            if(i==0) graphRatioPythiaToCombFitEta[i]     = CalculateGraphErrRatioToFit(graphRatioPythiaToCombFitEta[i], fitTsallisInvCrossSectionEtaComb[i]);
            else graphRatioPythiaToCombFitEta[i]     = CalculateGraphErrRatioToFit(graphRatioPythiaToCombFitEta[i], fitTCMInvCrossSectionEtaCombPlot[i]);
            DrawGammaSetMarkerTGraphErr(graphRatioPythiaToCombFitEta[i], 0, 0, pythia8color , pythia8color, widthLinesBoxes, kTRUE, pythia8color);

            if(i==4) while(graphRatioPythiaToCombFitEta[i]->GetX()[0] < 1.) graphRatioPythiaToCombFitEta[i]->RemovePoint(0);

          }
          if(histoPythiaInvCrossSectionEta[i]){
            histoRatioPythiaToCombFitEta[i]     = (TH1D*)histoPythiaInvCrossSectionEta[i]->Clone(Form("histoRatioPythiaToCombFitEta%i",i));
            if(i==0) histoRatioPythiaToCombFitEta[i]     = CalculateHistoRatioToFit(histoRatioPythiaToCombFitEta[i], fitTsallisInvCrossSectionEtaComb[i]);
            else histoRatioPythiaToCombFitEta[i]     = CalculateHistoRatioToFit(histoRatioPythiaToCombFitEta[i], fitTCMInvCrossSectionEtaCombPlot[i]);
            DrawGammaSetMarker(histoRatioPythiaToCombFitEta[i], 24, 1.5, pythia8color , pythia8color);
            histoRatioPythiaToCombFitEta[i]->SetLineWidth(widthCommonFit);

            if(i==4) histoRatioPythiaToCombFitEta[i]->GetXaxis()->SetRangeUser(0.8,maxXSpectraEta[i][10]);
          }
        }
    }

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Data to NLO Ratios -------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    TGraph* graphRatioNLODSS07ToCombFitPi0Half[6];
    TGraph* graphRatioNLODSS07ToCombFitPi0One[6];
    TGraph* graphRatioNLODSS07ToCombFitPi0Two[6];
    TGraphAsymmErrors* graphRatioNLODSS14ToCombFitPi0[6];
    TGraph* graphRatioNLOAESSSToCombFitEtaHalf[6];
    TGraph* graphRatioNLOAESSSToCombFitEtaOne[6];
    TGraph* graphRatioNLOAESSSToCombFitEtaTwo[6];
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
          if(graphNLODSS07InvCrossSectionPi0Half[i]){
            graphRatioNLODSS07ToCombFitPi0Half[i]     = (TGraph*)graphNLODSS07InvCrossSectionPi0Half[i]->Clone(Form("graphRatioNLODSS07ToCombFitPi0Half%i",i));
            graphRatioNLODSS07ToCombFitPi0Half[i]     = CalculateGraphRatioToFit(graphRatioNLODSS07ToCombFitPi0Half[i], fitTCMInvCrossSectionPi0CombPlot[i]);
            DrawGammaNLOTGraph( graphRatioNLODSS07ToCombFitPi0Half[i], widthCommonFit, 8, colorNLO);
          }else graphRatioNLODSS07ToCombFitPi0Half[i] = NULL;
          if(graphNLODSS07InvCrossSectionPi0One[i]){
            graphRatioNLODSS07ToCombFitPi0One[i]     = (TGraph*)graphNLODSS07InvCrossSectionPi0One[i]->Clone(Form("graphRatioNLODSS07ToCombFitPi0One%i",i));
            graphRatioNLODSS07ToCombFitPi0One[i]     = CalculateGraphRatioToFit(graphRatioNLODSS07ToCombFitPi0One[i], fitTCMInvCrossSectionPi0CombPlot[i]);
            DrawGammaNLOTGraph( graphRatioNLODSS07ToCombFitPi0One[i], widthCommonFit, 7, colorNLO);
          }else graphRatioNLODSS07ToCombFitPi0One[i] = NULL;
          if(graphNLODSS07InvCrossSectionPi0Two[i]){
            graphRatioNLODSS07ToCombFitPi0Two[i]     = (TGraph*)graphNLODSS07InvCrossSectionPi0Two[i]->Clone(Form("graphRatioNLODSS07ToCombFitPi0Two%i",i));
            graphRatioNLODSS07ToCombFitPi0Two[i]     = CalculateGraphRatioToFit(graphRatioNLODSS07ToCombFitPi0Two[i], fitTCMInvCrossSectionPi0CombPlot[i]);
            DrawGammaNLOTGraph( graphRatioNLODSS07ToCombFitPi0Two[i], widthCommonFit, 4, colorNLO);
          }else graphRatioNLODSS07ToCombFitPi0Two[i] = NULL;

          if(graphNLODSS14InvCrossSectionPi0[i]){
            graphRatioNLODSS14ToCombFitPi0[i]     = (TGraphAsymmErrors*)graphNLODSS14InvCrossSectionPi0[i]->Clone(Form("graphRatioNLODSS14ToCombFitPi0%i",i));
            graphRatioNLODSS14ToCombFitPi0[i]     = CalculateGraphErrRatioToFit(graphRatioNLODSS14ToCombFitPi0[i], fitTCMInvCrossSectionPi0CombPlot[i]);
            DrawGammaSetMarkerTGraphAsym(graphRatioNLODSS14ToCombFitPi0[i], 0, 0, colorNLO , colorNLO, widthLinesBoxes, kTRUE, colorNLO);
            ProduceGraphAsymmWithoutXErrors(graphRatioNLODSS14ToCombFitPi0[i]);
            graphRatioNLODSS14ToCombFitPi0[i]->SetLineWidth(widthCommonFit);
            graphRatioNLODSS14ToCombFitPi0[i]->SetLineColor(colorNLO);
            graphRatioNLODSS14ToCombFitPi0[i]->SetLineStyle(1);
            graphRatioNLODSS14ToCombFitPi0[i]->SetFillStyle(1001);
            graphRatioNLODSS14ToCombFitPi0[i]->SetFillColor(colorNLO);
          }else graphRatioNLODSS14ToCombFitPi0[i] = NULL;

          if(graphNLOAESSSInvCrossSectionEtaHalf[i]){
            graphRatioNLOAESSSToCombFitEtaHalf[i]     = (TGraph*)graphNLOAESSSInvCrossSectionEtaHalf[i]->Clone(Form("graphRatioNLOAESSSToCombFitEtaHalf%i",i));
            if(i==0) graphRatioNLOAESSSToCombFitEtaHalf[i]     = CalculateGraphRatioToFit(graphRatioNLOAESSSToCombFitEtaHalf[i], fitTsallisInvCrossSectionEtaComb[i]);
            else graphRatioNLOAESSSToCombFitEtaHalf[i]     = CalculateGraphRatioToFit(graphRatioNLOAESSSToCombFitEtaHalf[i], fitTCMInvCrossSectionEtaCombPlot[i]);
            DrawGammaNLOTGraph( graphRatioNLOAESSSToCombFitEtaHalf[i], widthCommonFit, 8, colorNLO);
          }else graphRatioNLOAESSSToCombFitEtaHalf[i] = NULL;
          if(graphNLOAESSSInvCrossSectionEtaOne[i]){
            graphRatioNLOAESSSToCombFitEtaOne[i]     = (TGraph*)graphNLOAESSSInvCrossSectionEtaOne[i]->Clone(Form("graphRatioNLOAESSSToCombFitEtaOne%i",i));
            if(i==0) graphRatioNLOAESSSToCombFitEtaOne[i]     = CalculateGraphRatioToFit(graphRatioNLOAESSSToCombFitEtaOne[i], fitTsallisInvCrossSectionEtaComb[i]);
            else graphRatioNLOAESSSToCombFitEtaOne[i]     = CalculateGraphRatioToFit(graphRatioNLOAESSSToCombFitEtaOne[i], fitTCMInvCrossSectionEtaCombPlot[i]);
            DrawGammaNLOTGraph( graphRatioNLOAESSSToCombFitEtaOne[i], widthCommonFit, 7, colorNLO);
          }else graphRatioNLOAESSSToCombFitEtaOne[i] = NULL;
          if(graphNLOAESSSInvCrossSectionEtaTwo[i]){
            graphRatioNLOAESSSToCombFitEtaTwo[i]     = (TGraph*)graphNLOAESSSInvCrossSectionEtaTwo[i]->Clone(Form("graphRatioNLOAESSSToCombFitEtaTwo%i",i));
            if(i==0) graphRatioNLOAESSSToCombFitEtaTwo[i]     = CalculateGraphRatioToFit(graphRatioNLOAESSSToCombFitEtaTwo[i], fitTsallisInvCrossSectionEtaComb[i]);
            else graphRatioNLOAESSSToCombFitEtaTwo[i]     = CalculateGraphRatioToFit(graphRatioNLOAESSSToCombFitEtaTwo[i], fitTCMInvCrossSectionEtaCombPlot[i]);
            DrawGammaNLOTGraph( graphRatioNLOAESSSToCombFitEtaTwo[i], widthCommonFit, 4, colorNLO);
          }else graphRatioNLOAESSSToCombFitEtaTwo[i] = NULL;
        }
    }

    //---------------------------------------------------------------------------------------------------------------
    //------------ Define canvas, pads and sizes for Spectrum (top) and Ratio (bottom) plot -------------------------
    //---------------------------------------------------------------------------------------------------------------

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
    TH2F * histoXSecDummyEta[numActiveMeas+1];

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


    histoXSecDummy[0]                           = new TH2F("histoXSecDummy_0","histoXSecDummy_0",11000,xRangeMinXSec,xRangeMaxXSec,1000,0.5e1,9e14);

    SetStyleHistoTH2ForGraphs(histoXSecDummy[0], "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3})",
                            0.85*textsizeLabelsXSec[0],textsizeLabelsXSec[0], 0.85*textsizeLabelsXSec[0], textsizeLabelsXSec[0], 1,0.2/(textsizeFacXSec[0]*marginXSec));
    histoXSecDummy[0]->GetXaxis()->SetMoreLogLabels();
    histoXSecDummy[0]->GetXaxis()->SetNoExponent(kTRUE);
    histoXSecDummy[0]->GetXaxis()->SetLabelOffset(+0.01);
    histoXSecDummy[0]->Draw();

    //__________________________________________ Loop over all ratio pads and draw
    for(Int_t i=0;i<numActiveMeas;i++){
        padXSec[i+1]->cd();
        padXSec[i+1]->SetLogx(1);
        histoXSecDummy[i+1]                     = new TH2F(Form("histoXSecDummy_%d",i+1),Form("histoXSecDummy_%d",i+1),1000,xRangeMinXSec,xRangeMaxXSec,1000,0.1,4);
        SetStyleHistoTH2ForGraphs(histoXSecDummy[i+1], "#it{p}_{T} (GeV/#it{c})","#frac{Data}{TCM fit}", 0.85*textsizeLabelsXSec[i+1], textsizeLabelsXSec[i+1],
                                  0.85*textsizeLabelsXSec[i+1],textsizeLabelsXSec[i+1], 1,0.2/(textsizeFacXSec[i+1]*marginXSec), 510, 505);
        histoXSecDummy[i+1]->GetYaxis()->SetRangeUser(0.4,1.99);
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

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    //__________________________________________ Draw combined spectra in first pad
    padXSec[0]->cd();
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
              fitTCMInvCrossSectionPi0CombPlot[i]->Draw("same");
              fitTsallisInvCrossSectionPi0Comb[i]->Draw("same");
//              graphPi0InvariantCrossSectionSys[i][10]->Draw("E2same");
//              graphPi0InvariantCrossSectionStat[i][10]->Draw("p,same,z");
            }
            for(Int_t iM = 0; iM<10; iM++){
              if(graphPi0InvariantCrossSectionSys[i][iM]&&graphPi0InvariantCrossSectionStat[i][iM]){
                graphPi0InvariantCrossSectionSys[i][iM]->Draw("E2same");
                graphPi0InvariantCrossSectionStat[i][iM]->Draw("p,same,z");
              }
            }
        }
    }
    Double_t rightalignDouble = 0.93;
    drawLatexAdd("ALICE",rightalignDouble,0.92,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",rightalignDouble,0.88,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

    TLegend* legendInvariantCrossSectionPi0    = GetAndSetLegend2(0.17, 0.03, 0.5, 0.03+textsizeLabelsXSec[0]*(numActiveMeas+1)+textsizeLabelsXSec[0], textSizeLabelsPixel);
    legendInvariantCrossSectionPi0->SetNColumns(1);
    legendInvariantCrossSectionPi0->SetMargin(0.2);
    TString legendScalingString[6] = {"", " (x10)"," (x10^{2})"," (x10^{3})"," (x10^{4})"," (x10^{5})"};
    Int_t legendRunningIndex = numActiveMeas-1;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            legendInvariantCrossSectionPi0->AddEntry(graphPi0InvariantCrossSectionSys[i][10],Form("pp, %s%s",energyLatex[i].Data(),legendScalingString[legendRunningIndex].Data()),"pf");
            legendRunningIndex-=1;
        }
    }
    legendInvariantCrossSectionPi0->AddEntry(fitTCMInvCrossSectionPi0CombPlot[exampleActiveMeas],"TCM fit","l");
    legendInvariantCrossSectionPi0->AddEntry(fitTsallisInvCrossSectionPi0Comb[exampleActiveMeas],"Tsallis fit","l");
    legendInvariantCrossSectionPi0->Draw();

    //__________________________________________ Draw ratios of individual spectra to combined fit in lower pads
    Int_t padCounter = 0;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            padXSec[padCounter+1]->cd();
            DrawGammaLines(xRangeMinXSec, xRangeMaxXSec , 1., 1., 1.2, kGray+2);
            for (Int_t j = 0; j < 10; j++){
                if(graphPi0InvariantCrossSectionSys[i][j]&&graphPi0InvariantCrossSectionStat[i][j]){
                    graphRatioIndivMeasToCombFitPi0Sys[i][j]->Draw("2,same");
                    graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][j]->Draw("p,same");
                }
            }
            padCounter+=1;
        }
    }
    canvasXSec->Print(Form("%s/Pi0_XSec_Ratios.%s",outputDir.Data(),suffix.Data()));

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    padXSec[0]->cd();
    histoXSecDummy[0]->Draw();
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
              if(graphNLODSS14InvCrossSectionPi0[i]) graphNLODSS14InvCrossSectionPi0[i]->Draw("same,e3");
              if(graphPythiaInvCrossSectionPi0[i]) graphPythiaInvCrossSectionPi0[i]->Draw("3,same");
              if(histoPythiaInvCrossSectionPi0[i]) histoPythiaInvCrossSectionPi0[i]->Draw("same,hist,l");
              fitTCMInvCrossSectionPi0CombPlot[i]->Draw("same");
              fitTsallisInvCrossSectionPi0Comb[i]->Draw("same");
              graphPi0InvariantCrossSectionSys[i][10]->Draw("E2same");
              graphPi0InvariantCrossSectionStat[i][10]->Draw("p,same,z");
            }
        }
    }
    drawLatexAdd("ALICE",rightalignDouble,0.92,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",rightalignDouble,0.88,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

    TLegend* legendInvariantXSecTheoryPi0    = GetAndSetLegend2(0.17, 0.03, 0.5, 0.03+textsizeLabelsXSec[0]*numActiveMeas+4*textsizeLabelsXSec[0], textSizeLabelsPixel);
    legendInvariantXSecTheoryPi0->SetNColumns(1);
    legendInvariantXSecTheoryPi0->SetMargin(0.2);
    legendRunningIndex = numActiveMeas-1;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            legendInvariantXSecTheoryPi0->AddEntry(graphPi0InvariantCrossSectionSys[i][10],Form("pp, %s%s",energyLatex[i].Data(),legendScalingString[legendRunningIndex].Data()),"pf");
            legendRunningIndex-=1;
        }
    }
    legendInvariantXSecTheoryPi0->AddEntry(fitTCMInvCrossSectionPi0CombPlot[exampleActiveMeas],"TCM fit","l");
    legendInvariantXSecTheoryPi0->AddEntry(fitTsallisInvCrossSectionPi0Comb[exampleActiveMeas],"Tsallis fit","l");
    legendInvariantXSecTheoryPi0->AddEntry(histoRatioPythiaToCombFitPi0[exampleActiveMeas],"PYTHIA 8.2, Monash 2013","l");
    legendInvariantXSecTheoryPi0->AddEntry(graphNLODSS14InvCrossSectionPi0[exampleActiveMeas],"NLO, PDF: MSTW - FF: DSS14 ","f");
    legendInvariantXSecTheoryPi0->Draw();
    //__________________________________________ Loop over all ratio pads and draw
    for(Int_t i=0;i<numActiveMeas;i++){
        padXSec[i+1]->cd();
        histoXSecDummy[i+1]->GetYaxis()->SetTitle("#frac{Data, Theory}{TCM fit}");
        histoXSecDummy[i+1]->DrawCopy();
        histoXSecDummy[i+1]->GetYaxis()->SetTitle("#frac{Data}{TCM fit}");
    }
    padCounter = 0;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            padXSec[padCounter+1]->cd();
            DrawGammaLines(xRangeMinXSec, xRangeMaxXSec , 1., 1., 1.2, kGray+2);
            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
              if(graphRatioNLODSS14ToCombFitPi0[i]) graphRatioNLODSS14ToCombFitPi0[i]->Draw("same,e3");
              graphRatioPythiaToCombFitPi0[i]->Draw("3,same");
              histoRatioPythiaToCombFitPi0[i]->Draw("same,hist,l");
              graphRatioIndivMeasToCombFitPi0Sys[i][10]->Draw("2,same");
              graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][10]->Draw("p,same");
            }
            padCounter+=1;
        }
    }
    canvasXSec->Print(Form("%s/Pi0_XSec_FullRatios_Pythia_TheoryDSS14.%s",outputDir.Data(),suffix.Data()));

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    padXSec[0]->cd();
    histoXSecDummy[0]->Draw();
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
                if(graphPythiaInvCrossSectionPi0[i]) graphPythiaInvCrossSectionPi0[i]->Draw("3,same");
                if(histoPythiaInvCrossSectionPi0[i]) histoPythiaInvCrossSectionPi0[i]->Draw("same,hist,l");
                fitTCMInvCrossSectionPi0CombPlot[i]->Draw("same");
                fitTsallisInvCrossSectionPi0Comb[i]->Draw("same");
                graphPi0InvariantCrossSectionSys[i][10]->Draw("E2same");
                graphPi0InvariantCrossSectionStat[i][10]->Draw("p,same,z");
            }
        }
    }
    drawLatexAdd("ALICE",rightalignDouble,0.92,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",rightalignDouble,0.88,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

    TLegend* legendInvariantXSecPi0    = GetAndSetLegend2(0.17, 0.03, 0.5, 0.03+textsizeLabelsXSec[0]*numActiveMeas+4*textsizeLabelsXSec[0], textSizeLabelsPixel);
    legendInvariantXSecPi0->SetNColumns(1);
    legendInvariantXSecPi0->SetMargin(0.2);
    legendRunningIndex = numActiveMeas-1;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            legendInvariantXSecPi0->AddEntry(graphPi0InvariantCrossSectionSys[i][10],Form("pp, %s%s",energyLatex[i].Data(),legendScalingString[legendRunningIndex].Data()),"pf");
            legendRunningIndex-=1;
        }
    }
    legendInvariantXSecPi0->AddEntry(fitTCMInvCrossSectionPi0CombPlot[exampleActiveMeas],"TCM fit","l");
    legendInvariantXSecPi0->AddEntry(fitTsallisInvCrossSectionPi0Comb[exampleActiveMeas],"Tsallis fit","l");
    legendInvariantXSecPi0->AddEntry(histoRatioPythiaToCombFitPi0[exampleActiveMeas],"PYTHIA 8.2, Monash 2013","l");
    legendInvariantXSecPi0->Draw();
    //__________________________________________ Loop over all ratio pads and draw
    for(Int_t i=0;i<numActiveMeas;i++){
        padXSec[i+1]->cd();
        histoXSecDummy[i+1]->GetYaxis()->SetTitle("#frac{Data, Theory}{TCM fit}");
        histoXSecDummy[i+1]->DrawCopy();
        histoXSecDummy[i+1]->GetYaxis()->SetTitle("#frac{Data}{TCM fit}");
    }
    padCounter = 0;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            padXSec[padCounter+1]->cd();
            DrawGammaLines(xRangeMinXSec, xRangeMaxXSec , 1., 1., 1.2, kGray+2);
            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
                graphRatioPythiaToCombFitPi0[i]->Draw("3,same");
                histoRatioPythiaToCombFitPi0[i]->Draw("same,hist,l");
                graphRatioIndivMeasToCombFitPi0Sys[i][10]->Draw("2,same");
                graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][10]->Draw("p,same");
            }
            padCounter+=1;
        }
    }
    canvasXSec->Print(Form("%s/Pi0_XSec_FullRatios_Pythia.%s",outputDir.Data(),suffix.Data()));

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    padXSec[0]->cd();
    histoXSecDummy[0]->Draw();
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
              if(graphNLODSS07InvCrossSectionPi0Half[i]) graphNLODSS07InvCrossSectionPi0Half[i]->Draw("same,c");
              if(graphNLODSS07InvCrossSectionPi0One[i]) graphNLODSS07InvCrossSectionPi0One[i]->Draw("same,c");
              if(graphNLODSS07InvCrossSectionPi0Two[i]) graphNLODSS07InvCrossSectionPi0Two[i]->Draw("same,c");
              if(graphPythiaInvCrossSectionPi0[i]) graphPythiaInvCrossSectionPi0[i]->Draw("3,same");
              if(histoPythiaInvCrossSectionPi0[i]) histoPythiaInvCrossSectionPi0[i]->Draw("same,hist,l");
              fitTCMInvCrossSectionPi0CombPlot[i]->Draw("same");
              fitTsallisInvCrossSectionPi0Comb[i]->Draw("same");
              graphPi0InvariantCrossSectionSys[i][10]->Draw("E2same");
              graphPi0InvariantCrossSectionStat[i][10]->Draw("p,same,z");
            }
        }
    }
    drawLatexAdd("ALICE",rightalignDouble,0.92,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",rightalignDouble,0.88,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

    TLegend* legendInvariantCrossSectionPi044    = GetAndSetLegend2(0.696, 0.652, 0.94, 0.652+textsizeLabelsXSec[0]*(numActiveMeas-2)+textsizeLabelsXSec[0], 0.9*textSizeLabelsPixel);
    legendInvariantCrossSectionPi044->SetNColumns(1);
    legendInvariantCrossSectionPi044->SetMargin(0.3);
    legendInvariantCrossSectionPi044->AddEntry((TObject*)0,"","");
    legendInvariantCrossSectionPi044->AddEntry(graphNLODSS07InvCrossSectionPi0Half[exampleActiveMeas],"#mu = 0.5 #it{p}_{T}","l");
    legendInvariantCrossSectionPi044->AddEntry(graphNLODSS07InvCrossSectionPi0One[exampleActiveMeas],"#mu = #it{p}_{T}","l");
    legendInvariantCrossSectionPi044->AddEntry(graphNLODSS07InvCrossSectionPi0Two[exampleActiveMeas],"#mu = 2 #it{p}_{T}","l");
    legendInvariantCrossSectionPi044->Draw();
    drawLatexAdd("NLO,",        0.7,0.81,0.9*textsizeLabelsXSec[0],kFALSE,kFALSE,kFALSE);
    drawLatexAdd("PDF:CTEQ6M5", 0.7,0.78,0.9*textsizeLabelsXSec[0],kFALSE,kFALSE,kFALSE);
    drawLatexAdd("FF:DSS07",    0.7,0.75,0.9*textsizeLabelsXSec[0],kFALSE,kFALSE,kFALSE);

    legendInvariantXSecPi0->Draw();
    //__________________________________________ Loop over all ratio pads and draw
    for(Int_t i=0;i<numActiveMeas;i++){
        padXSec[i+1]->cd();
        histoXSecDummy[i+1]->GetYaxis()->SetTitle("#frac{Data, Theory}{TCM fit}");
        histoXSecDummy[i+1]->GetYaxis()->SetRangeUser(0.4,3.95);
        histoXSecDummy[i+1]->DrawCopy();
        histoXSecDummy[i+1]->GetYaxis()->SetRangeUser(0.4,1.99);
        histoXSecDummy[i+1]->GetYaxis()->SetTitle("#frac{Data}{TCM fit}");
    }
    padCounter = 0;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            padXSec[padCounter+1]->cd();
            DrawGammaLines(xRangeMinXSec, xRangeMaxXSec , 1., 1., 1.2, kGray+2);
            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
              if(graphRatioNLODSS07ToCombFitPi0Half[i]) graphRatioNLODSS07ToCombFitPi0Half[i]->Draw("same,c");
              if(graphRatioNLODSS07ToCombFitPi0One[i]) graphRatioNLODSS07ToCombFitPi0One[i]->Draw("same,c");
              if(graphRatioNLODSS07ToCombFitPi0Two[i]) graphRatioNLODSS07ToCombFitPi0Two[i]->Draw("same,c");
              graphRatioPythiaToCombFitPi0[i]->Draw("3,same");
              histoRatioPythiaToCombFitPi0[i]->Draw("same,hist,l");
              graphRatioIndivMeasToCombFitPi0Sys[i][10]->Draw("2,same");
              graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][10]->Draw("p,same");
            }
            padCounter+=1;
        }
    }
    canvasXSec->Print(Form("%s/Pi0_XSec_FullRatios_Pythia_TheoryDSS07.%s",outputDir.Data(),suffix.Data()));

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    padXSec[0]->cd();
    histoXSecDummy[0]->Draw();
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
                fitTCMInvCrossSectionPi0CombPlot[i]->Draw("same");
                fitTsallisInvCrossSectionPi0Comb[i]->Draw("same");
                graphPi0InvariantCrossSectionSys[i][10]->Draw("E2same");
                graphPi0InvariantCrossSectionStat[i][10]->Draw("p,same,z");
            }
        }
    }
    drawLatexAdd("ALICE",rightalignDouble,0.92,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",rightalignDouble,0.88,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

    legendInvariantCrossSectionPi0->Draw();
    //__________________________________________ Loop over all ratio pads and draw
    for(Int_t i=0;i<numActiveMeas;i++){
        padXSec[i+1]->cd();
        histoXSecDummy[i+1]->DrawCopy();
    }
    padCounter = 0;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            padXSec[padCounter+1]->cd();
            DrawGammaLines(xRangeMinXSec, xRangeMaxXSec , 1., 1., 1.2, kGray+2);
            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
                graphRatioIndivMeasToCombFitPi0Sys[i][10]->Draw("2,same");
                graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][10]->Draw("p,same");
            }
            padCounter+=1;
        }
    }
    canvasXSec->Print(Form("%s/Pi0_XSec_FullRatios.%s",outputDir.Data(),suffix.Data()));

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    TCanvas* canvasXSectionPi0  = new TCanvas("canvasXSectionPi0","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasXSectionPi0, 0.14, 0.01, 0.01, 0.07);
    canvasXSectionPi0->SetLogx();
    canvasXSectionPi0->SetLogy();
    histoXSecDummy[0]->GetXaxis()->SetLabelOffset(0.0);
    histoXSecDummy[0]->GetXaxis()->SetTitleOffset(0.8);
    histoXSecDummy[0]->Draw();
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
              fitTCMInvCrossSectionPi0CombPlot[i]->Draw("same");
              fitTsallisInvCrossSectionPi0Comb[i]->Draw("same");
              graphPi0InvariantCrossSectionSys[i][10]->Draw("E2same");
              graphPi0InvariantCrossSectionStat[i][10]->Draw("p,same,z");
            }
        }
    }
    drawLatexAdd("ALICE",rightalignDouble,0.92,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",rightalignDouble,0.88,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

    TLegend* legendInvariantCrossSectionPi02    = GetAndSetLegend2(0.17, 0.10, 0.5, 0.10+textsizeLabelsXSec[0]*(numActiveMeas+1)+textsizeLabelsXSec[0], textSizeLabelsPixel);
    legendInvariantCrossSectionPi02->SetNColumns(1);
    legendInvariantCrossSectionPi02->SetMargin(0.2);
    legendRunningIndex = numActiveMeas-1;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            legendInvariantCrossSectionPi02->AddEntry(graphPi0InvariantCrossSectionSys[i][10],Form("pp, %s%s",energyLatex[i].Data(),legendScalingString[legendRunningIndex].Data()),"pf");
            legendRunningIndex-=1;
        }
    }
    legendInvariantCrossSectionPi02->AddEntry(fitTCMInvCrossSectionPi0CombPlot[exampleActiveMeas],"TCM fit","l");
    legendInvariantCrossSectionPi02->AddEntry(fitTsallisInvCrossSectionPi0Comb[exampleActiveMeas],"Tsallis fit","l");
    legendInvariantCrossSectionPi02->Draw();

    canvasXSectionPi0->Print(Form("%s/Pi0_XSec.%s",outputDir.Data(),suffix.Data()));

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    histoXSecDummy[0]->Draw();
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
              if(graphPythiaInvCrossSectionPi0[i]) graphPythiaInvCrossSectionPi0[i]->Draw("3,same");
              if(histoPythiaInvCrossSectionPi0[i]) histoPythiaInvCrossSectionPi0[i]->Draw("same,hist,l");
              fitTCMInvCrossSectionPi0CombPlot[i]->Draw("same");
              fitTsallisInvCrossSectionPi0Comb[i]->Draw("same");
              graphPi0InvariantCrossSectionSys[i][10]->Draw("E2same");
              graphPi0InvariantCrossSectionStat[i][10]->Draw("p,same,z");
            }
        }
    }
    drawLatexAdd("ALICE",rightalignDouble,0.92,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",rightalignDouble,0.88,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

    TLegend* legendInvariantCrossSectionPi03    = GetAndSetLegend2(0.17, 0.10, 0.5, 0.10+textsizeLabelsXSec[0]*(numActiveMeas+2)+textsizeLabelsXSec[0], textSizeLabelsPixel);
    legendInvariantCrossSectionPi03->SetNColumns(1);
    legendInvariantCrossSectionPi03->SetMargin(0.2);
    legendRunningIndex = numActiveMeas-1;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            legendInvariantCrossSectionPi03->AddEntry(graphPi0InvariantCrossSectionSys[i][10],Form("pp, %s%s",energyLatex[i].Data(),legendScalingString[legendRunningIndex].Data()),"pf");
            legendRunningIndex-=1;
        }
    }
    legendInvariantCrossSectionPi03->AddEntry(fitTCMInvCrossSectionPi0CombPlot[exampleActiveMeas],"TCM fit","l");
    legendInvariantCrossSectionPi03->AddEntry(fitTsallisInvCrossSectionPi0Comb[exampleActiveMeas],"Tsallis fit","l");
    legendInvariantCrossSectionPi03->AddEntry(histoRatioPythiaToCombFitPi0[exampleActiveMeas],"PYTHIA 8.2, Monash 2013","l");
    legendInvariantCrossSectionPi03->Draw();

    canvasXSectionPi0->Print(Form("%s/Pi0_XSec_Pythia.%s",outputDir.Data(),suffix.Data()));

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    histoXSecDummy[0]->Draw();
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
              if(graphNLODSS07InvCrossSectionPi0Half[i]) graphNLODSS07InvCrossSectionPi0Half[i]->Draw("same,c");
              if(graphNLODSS07InvCrossSectionPi0One[i]) graphNLODSS07InvCrossSectionPi0One[i]->Draw("same,c");
              if(graphNLODSS07InvCrossSectionPi0Two[i]) graphNLODSS07InvCrossSectionPi0Two[i]->Draw("same,c");
              if(graphPythiaInvCrossSectionPi0[i]) graphPythiaInvCrossSectionPi0[i]->Draw("3,same");
              if(histoPythiaInvCrossSectionPi0[i]) histoPythiaInvCrossSectionPi0[i]->Draw("same,hist,l");
              fitTCMInvCrossSectionPi0CombPlot[i]->Draw("same");
              fitTsallisInvCrossSectionPi0Comb[i]->Draw("same");
              graphPi0InvariantCrossSectionSys[i][10]->Draw("E2same");
              graphPi0InvariantCrossSectionStat[i][10]->Draw("p,same,z");
            }
        }
    }
    drawLatexAdd("ALICE",rightalignDouble,0.92,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",rightalignDouble,0.88,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

    legendInvariantCrossSectionPi044->Draw();
    drawLatexAdd("NLO,",        0.7,0.81,0.9*textsizeLabelsXSec[0],kFALSE,kFALSE,kFALSE);
    drawLatexAdd("PDF:CTEQ6M5", 0.7,0.78,0.9*textsizeLabelsXSec[0],kFALSE,kFALSE,kFALSE);
    drawLatexAdd("FF:DSS07",    0.7,0.75,0.9*textsizeLabelsXSec[0],kFALSE,kFALSE,kFALSE);
    legendInvariantCrossSectionPi03->Draw();

    canvasXSectionPi0->Print(Form("%s/Pi0_XSec_Pythia_TheoryDSS07.%s",outputDir.Data(),suffix.Data()));

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    histoXSecDummy[0]->Draw();
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
              if(graphNLODSS14InvCrossSectionPi0[i]) graphNLODSS14InvCrossSectionPi0[i]->Draw("same,e3");
              if(graphPythiaInvCrossSectionPi0[i]) graphPythiaInvCrossSectionPi0[i]->Draw("3,same");
              if(histoPythiaInvCrossSectionPi0[i]) histoPythiaInvCrossSectionPi0[i]->Draw("same,hist,l");
              fitTCMInvCrossSectionPi0CombPlot[i]->Draw("same");
              fitTsallisInvCrossSectionPi0Comb[i]->Draw("same");
              graphPi0InvariantCrossSectionSys[i][10]->Draw("E2same");
              graphPi0InvariantCrossSectionStat[i][10]->Draw("p,same,z");
            }
        }
    }
    drawLatexAdd("ALICE",rightalignDouble,0.92,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",rightalignDouble,0.88,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

    TLegend* legendInvariantCrossSectionPi05    = GetAndSetLegend2(0.17, 0.10, 0.5, 0.10+textsizeLabelsXSec[0]*(numActiveMeas+3)+textsizeLabelsXSec[0], textSizeLabelsPixel);
    legendInvariantCrossSectionPi05->SetNColumns(1);
    legendInvariantCrossSectionPi05->SetMargin(0.2);
    legendRunningIndex = numActiveMeas-1;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            legendInvariantCrossSectionPi05->AddEntry(graphPi0InvariantCrossSectionSys[i][10],Form("pp, %s%s",energyLatex[i].Data(),legendScalingString[legendRunningIndex].Data()),"pf");
            legendRunningIndex-=1;
        }
    }
    legendInvariantCrossSectionPi05->AddEntry(fitTCMInvCrossSectionPi0CombPlot[exampleActiveMeas],"TCM fit","l");
    legendInvariantCrossSectionPi05->AddEntry(fitTsallisInvCrossSectionPi0Comb[exampleActiveMeas],"Tsallis fit","l");
    legendInvariantCrossSectionPi05->AddEntry(histoRatioPythiaToCombFitPi0[exampleActiveMeas],"PYTHIA 8.2, Monash 2013","l");
    legendInvariantCrossSectionPi05->AddEntry(graphNLODSS14InvCrossSectionPi0[exampleActiveMeas],"NLO, PDF: MSTW - FF: DSS14 ","f");
    legendInvariantCrossSectionPi05->Draw();

    canvasXSectionPi0->Print(Form("%s/Pi0_XSec_Pythia_TheoryDSS14.%s",outputDir.Data(),suffix.Data()));
    histoXSecDummy[0]->GetXaxis()->SetLabelOffset(+0.01);
    histoXSecDummy[0]->GetXaxis()->SetTitleOffset(1.0);

    //--------------------------------------------------------------------------------------------
    //--------------------Define canvas, pads and sizes for ratio only plot-----------------------
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
        histoRatiosDummy[i]                     = new TH2F(Form("histoRatiosDummy_%d",i),Form("histoRatiosDummy_%d",i),1000,xRangeMinRatios,xRangeMaxRatios,1000,0.1,4.0);
        SetStyleHistoTH2ForGraphs(histoRatiosDummy[i], "#it{p}_{T} (GeV/#it{c})","#frac{Data}{TCM fit}", 0.85*textsizeLabelsRatios[i], textsizeLabelsRatios[i],
                                  0.85*textsizeLabelsRatios[i],textsizeLabelsRatios[i], 1,0.2/(textsizeFacRatios[i]*marginRatios), 510, 505);
        histoRatiosDummy[i]->GetYaxis()->SetRangeUser(0.4,1.99);
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

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    TLegend* legendRatiosEnergy[6];
    TLegend* legendRatiosPythia = NULL;
    TH1D* tempPythia = (TH1D*) histoRatioPythiaToCombFitPi0[exampleActiveMeas]->Clone("tempPythiaEx");
    tempPythia->SetLineWidth(4);

    padCounter = 0;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            padRatios[padCounter]->cd();
            DrawGammaLines(xRangeMinRatios, xRangeMaxRatios, 1., 1., 1.2, kGray+2);
            if(padCounter==0){
              drawLatexAdd("ALICE",rightalignDouble,0.84,3*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
              drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",rightalignDouble,0.73,3*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
            }
            for (Int_t j = 0; j < 10; j++){
                if(graphPi0InvariantCrossSectionSys[i][j]&&graphPi0InvariantCrossSectionStat[i][j]){
                    graphRatioIndivMeasToCombFitPi0Sys[i][j]->Draw("2,same");
                    graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][j]->Draw("p,same");
                }
            }
            padCounter+=1;
        }
    }
    for(Int_t i=0;i<numActiveMeas;i++) padRatios[i]->RedrawAxis();
    canvasRatios->Print(Form("%s/Pi0_Ratios.%s",outputDir.Data(),suffix.Data()));

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

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
            DrawGammaLines(xRangeMinRatios, xRangeMaxRatios, 1., 1., 1.2, kGray+2);
            if(padCounter==0){
              drawLatexAdd("ALICE",rightalignDouble,0.84,3*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
              drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",rightalignDouble,0.73,3*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
            }
            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
                graphRatioIndivMeasToCombFitPi0Sys[i][10]->Draw("2,same");
                graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][10]->Draw("p,same");
                if(i==0) legendRatiosEnergy[i] = GetAndSetLegend2(0.15, 0.81, 0.5, 0.76+textsizeLabelsXSec[0]*(numActiveMeas+1)+textsizeLabelsXSec[0], textSizeLabelsPixel);
                else legendRatiosEnergy[i] = GetAndSetLegend2(0.15, 0.75, 0.5, 0.75+textsizeLabelsXSec[0]*(numActiveMeas+1)+textsizeLabelsXSec[0], textSizeLabelsPixel);
                legendRatiosEnergy[i]->SetNColumns(1);
                legendRatiosEnergy[i]->SetMargin(0.2);
                legendRatiosEnergy[i]->AddEntry(graphPi0InvariantCrossSectionSys[i][10],Form("pp, %s",energyLatex[i].Data()),"pf");
                legendRatiosEnergy[i]->Draw();
            }
            padCounter+=1;
        }
    }
    for(Int_t i=0;i<numActiveMeas;i++) padRatios[i]->RedrawAxis();
    canvasRatios->Print(Form("%s/Pi0_FullRatios.%s",outputDir.Data(),suffix.Data()));

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    canvasRatios->cd();
    //__________________________________________ Set ratio pad properties
    for(Int_t i=0;i<numActiveMeas;i++){
        padRatios[i]->Draw();
    }
    //__________________________________________ Loop over all ratio pads and draw
    for(Int_t i=0;i<numActiveMeas;i++){
        padRatios[i]->cd();
        padRatios[i]->SetLogx(1);
        histoRatiosDummy[i]->GetYaxis()->SetTitle("#frac{Data, Theory}{TCM fit}");
        histoRatiosDummy[i]->DrawCopy();
    }
    padCounter = 0;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            padRatios[padCounter]->cd();
            DrawGammaLines(xRangeMinRatios, xRangeMaxRatios, 1., 1., 1.2, kGray+2);
            if(padCounter==0){
              drawLatexAdd("ALICE",rightalignDouble,0.84,3*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
              drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",rightalignDouble,0.73,3*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
            }else if(padCounter==3){
              legendRatiosPythia = GetAndSetLegend2(0.7, 0.74, 0.95, 0.74+textsizeLabelsXSec[0]*(numActiveMeas)+textsizeLabelsXSec[0], textSizeLabelsPixel);
              legendRatiosPythia->SetNColumns(1);
              legendRatiosPythia->SetMargin(0.2);
              legendRatiosPythia->AddEntry(tempPythia,"PYTHIA 8.2,","l");
              legendRatiosPythia->AddEntry((TObject*)0,"Monash 2013","");
              legendRatiosPythia->Draw();
            }
            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
              if(i==1){
                graphRatioPythiaToCombFitPi0[i]->RemovePoint(0);
                histoRatioPythiaToCombFitPi0[i]->GetXaxis()->SetRangeUser(0.5,maxXSpectra[i][10]);
              }
              graphRatioPythiaToCombFitPi0[i]->Draw("3,same");
              histoRatioPythiaToCombFitPi0[i]->Draw("same,hist,l");
              graphRatioIndivMeasToCombFitPi0Sys[i][10]->Draw("2,same");
              graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][10]->Draw("p,same");
              if(i==0) legendRatiosEnergy[i] = GetAndSetLegend2(0.15, 0.81, 0.5, 0.76+textsizeLabelsXSec[0]*(numActiveMeas+1)+textsizeLabelsXSec[0], textSizeLabelsPixel);
              else legendRatiosEnergy[i] = GetAndSetLegend2(0.15, 0.75, 0.5, 0.75+textsizeLabelsXSec[0]*(numActiveMeas+1)+textsizeLabelsXSec[0], textSizeLabelsPixel);
              legendRatiosEnergy[i]->SetNColumns(1);
              legendRatiosEnergy[i]->SetMargin(0.2);
              legendRatiosEnergy[i]->AddEntry(graphPi0InvariantCrossSectionSys[i][10],Form("pp, %s",energyLatex[i].Data()),"pf");
              legendRatiosEnergy[i]->Draw();
            }
            padCounter+=1;
        }
    }
    for(Int_t i=0;i<numActiveMeas;i++) padRatios[i]->RedrawAxis();
    canvasRatios->Print(Form("%s/Pi0_FullRatios_Pythia.%s",outputDir.Data(),suffix.Data()));

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    canvasRatios->cd();
    //__________________________________________ Set ratio pad properties
    for(Int_t i=0;i<numActiveMeas;i++){
        padRatios[i]->Draw();
    }
    //__________________________________________ Loop over all ratio pads and draw
    for(Int_t i=0;i<numActiveMeas;i++){
        padRatios[i]->cd();
        padRatios[i]->SetLogx(1);
        histoRatiosDummy[i]->GetYaxis()->SetRangeUser(0.4,2.49);
        histoRatiosDummy[i]->DrawCopy();
        histoRatiosDummy[i]->GetYaxis()->SetRangeUser(0.4,1.99);
    }
    padCounter = 0;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            padRatios[padCounter]->cd();
            DrawGammaLines(xRangeMinRatios, xRangeMaxRatios, 1., 1., 1.2, kGray+2);
            if(padCounter==0){
              drawLatexAdd("ALICE",rightalignDouble,0.84,3*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
              drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",rightalignDouble,0.73,3*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
            }else if(padCounter==3){
              if(legendRatiosPythia) legendRatiosPythia->Draw();
              legendInvariantXSecTheoryPi0 = GetAndSetLegend2(0.52, 0.69, 0.69, 0.69+textsizeLabelsXSec[0]*(numActiveMeas+1)+textsizeLabelsXSec[0], textSizeLabelsPixel);
              legendInvariantXSecTheoryPi0->SetNColumns(1);
              legendInvariantXSecTheoryPi0->SetMargin(0.2);
              legendInvariantXSecTheoryPi0->AddEntry(graphNLODSS14InvCrossSectionPi0[exampleActiveMeas],"NLO,","f");
              legendInvariantXSecTheoryPi0->AddEntry((TObject*)0,"#scale[0.8]{PDF: MSTW}","");
              legendInvariantXSecTheoryPi0->AddEntry((TObject*)0,"#scale[0.8]{FF: DSS14}","");
              legendInvariantXSecTheoryPi0->Draw();
            }
            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
              if(graphRatioNLODSS14ToCombFitPi0[i]) graphRatioNLODSS14ToCombFitPi0[i]->Draw("same,e3");
              graphRatioPythiaToCombFitPi0[i]->Draw("3,same");
              histoRatioPythiaToCombFitPi0[i]->Draw("same,hist,l");
              graphRatioIndivMeasToCombFitPi0Sys[i][10]->Draw("2,same");
              graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][10]->Draw("p,same");
              if(i==0) legendRatiosEnergy[i] = GetAndSetLegend2(0.15, 0.81, 0.5, 0.76+textsizeLabelsXSec[0]*(numActiveMeas+1)+textsizeLabelsXSec[0], textSizeLabelsPixel);
              else legendRatiosEnergy[i] = GetAndSetLegend2(0.15, 0.75, 0.5, 0.75+textsizeLabelsXSec[0]*(numActiveMeas+1)+textsizeLabelsXSec[0], textSizeLabelsPixel);
              legendRatiosEnergy[i]->SetNColumns(1);
              legendRatiosEnergy[i]->SetMargin(0.2);
              legendRatiosEnergy[i]->AddEntry(graphPi0InvariantCrossSectionSys[i][10],Form("pp, %s",energyLatex[i].Data()),"pf");
              legendRatiosEnergy[i]->Draw();
            }
            padCounter+=1;
        }
    }
    for(Int_t i=0;i<numActiveMeas;i++) padRatios[i]->RedrawAxis();
    canvasRatios->Print(Form("%s/Pi0_FullRatios_Pythia_DSS14.%s",outputDir.Data(),suffix.Data()));

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    canvasRatios->cd();
    //__________________________________________ Set ratio pad properties
    for(Int_t i=0;i<numActiveMeas;i++){
        padRatios[i]->Draw();
    }
    //__________________________________________ Loop over all ratio pads and draw
    for(Int_t i=0;i<numActiveMeas;i++){
        padRatios[i]->cd();
        padRatios[i]->SetLogx(1);
        histoRatiosDummy[i]->GetYaxis()->SetRangeUser(0.4,3.95);
        histoRatiosDummy[i]->DrawCopy();
        histoRatiosDummy[i]->GetYaxis()->SetRangeUser(0.4,1.99);
    }
    padCounter = 0;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            padRatios[padCounter]->cd();
            DrawGammaLines(xRangeMinRatios, xRangeMaxRatios, 1., 1., 1.2, kGray+2);
            if(padCounter==0){
              drawLatexAdd("ALICE",rightalignDouble,0.84,3*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
              drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",rightalignDouble,0.73,3*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
            }else if(padCounter==3){
              legendRatiosPythia = GetAndSetLegend2(0.15, 0.6, 0.5, 0.6+textsizeLabelsXSec[0]*(numActiveMeas)+textsizeLabelsXSec[0], textSizeLabelsPixel);
              legendRatiosPythia->SetNColumns(1);
              legendRatiosPythia->SetMargin(0.2);
              legendRatiosPythia->AddEntry(tempPythia,"PYTHIA 8.2,","l");
              legendRatiosPythia->AddEntry((TObject*)0,"Monash 2013","");
              legendRatiosPythia->Draw();

              legendInvariantXSecTheoryPi0 = GetAndSetLegend2(0.75, 0.595, 0.95, 0.595+textsizeLabelsXSec[0]*(numActiveMeas+3)+textsizeLabelsXSec[0], textSizeLabelsPixel);
              legendInvariantXSecTheoryPi0->SetNColumns(1);
              legendInvariantXSecTheoryPi0->SetMargin(0.2);
              legendInvariantXSecTheoryPi0->AddEntry((TObject*)0,"","");
              legendInvariantXSecTheoryPi0->AddEntry(graphNLODSS07InvCrossSectionPi0Half[exampleActiveMeas],"#mu = 0.5 #it{p}_{T}","l");
              legendInvariantXSecTheoryPi0->AddEntry(graphNLODSS07InvCrossSectionPi0One[exampleActiveMeas],"#mu = #it{p}_{T}","l");
              legendInvariantXSecTheoryPi0->AddEntry(graphNLODSS07InvCrossSectionPi0Two[exampleActiveMeas],"#mu = 2 #it{p}_{T}","l");
              legendInvariantXSecTheoryPi0->Draw();
              drawLatexAdd("NLO,",0.83,0.86,2.5*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
              drawLatexAdd("PDF:CTEQ6M5",0.974,0.51,2.2*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
              drawLatexAdd("FF:DSS07",0.9,0.44,2.2*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

            }
            if(graphPi0InvariantCrossSectionSys[i][10]&&graphPi0InvariantCrossSectionStat[i][10]){
              if(graphRatioNLODSS07ToCombFitPi0Half[i]) graphRatioNLODSS07ToCombFitPi0Half[i]->Draw("same,c");
              if(graphRatioNLODSS07ToCombFitPi0One[i]) graphRatioNLODSS07ToCombFitPi0One[i]->Draw("same,c");
              if(graphRatioNLODSS07ToCombFitPi0Two[i]) graphRatioNLODSS07ToCombFitPi0Two[i]->Draw("same,c");
              graphRatioPythiaToCombFitPi0[i]->Draw("3,same");
              histoRatioPythiaToCombFitPi0[i]->Draw("same,hist,l");
              graphRatioIndivMeasToCombFitPi0Sys[i][10]->Draw("2,same");
              graphRatioIndivMeasToCombFitPi0Stat_woXErr[i][10]->Draw("p,same");
              if(i==0) legendRatiosEnergy[i] = GetAndSetLegend2(0.15, 0.81, 0.5, 0.76+textsizeLabelsXSec[0]*(numActiveMeas+1)+textsizeLabelsXSec[0], textSizeLabelsPixel);
              else legendRatiosEnergy[i] = GetAndSetLegend2(0.15, 0.75, 0.5, 0.75+textsizeLabelsXSec[0]*(numActiveMeas+1)+textsizeLabelsXSec[0], textSizeLabelsPixel);
              legendRatiosEnergy[i]->SetNColumns(1);
              legendRatiosEnergy[i]->SetMargin(0.2);
              legendRatiosEnergy[i]->AddEntry(graphPi0InvariantCrossSectionSys[i][10],Form("pp, %s",energyLatex[i].Data()),"pf");
              legendRatiosEnergy[i]->Draw();
            }
            padCounter+=1;
        }
    }
    for(Int_t i=0;i<numActiveMeas;i++) padRatios[i]->RedrawAxis();
    canvasRatios->Print(Form("%s/Pi0_FullRatios_Pythia_DSS07.%s",outputDir.Data(),suffix.Data()));

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------EtaMeson -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------EtaMeson -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------EtaMeson -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------EtaMeson -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------EtaMeson -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------EtaMeson -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------


    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    canvasXSec->cd();
    padXSec[0]->Draw();
    //__________________________________________ Set ratio pad properties
    for(Int_t i=0;i<numActiveMeas;i++){
        padXSec[i+1]->Draw();
    }

    //__________________________________________ Draw combined spectra in first pad
    padXSec[0]->cd();

    histoXSecDummyEta[0]                           = new TH2F("histoXSecDummyEta_0","histoXSecDummyEta_0",11000,xRangeMinXSec,xRangeMaxXSec,1000,1e2,2e13);
    SetStyleHistoTH2ForGraphs(histoXSecDummyEta[0], "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3})",
                            0.85*textsizeLabelsXSec[0],textsizeLabelsXSec[0], 0.85*textsizeLabelsXSec[0], textsizeLabelsXSec[0], 1,0.2/(textsizeFacXSec[0]*marginXSec));
    histoXSecDummyEta[0]->GetXaxis()->SetMoreLogLabels();
    histoXSecDummyEta[0]->GetXaxis()->SetNoExponent(kTRUE);
    histoXSecDummyEta[0]->GetXaxis()->SetLabelOffset(+0.01);
    histoXSecDummyEta[0]->Draw();

//    padCounter = 0;
//    for (Int_t i = 5; i > -1; i--){
//        if(includeEnergy[i]){
//            padXSec[padCounter+1]->cd();
//            if(i==0) histoXSecDummy[padCounter+1]->GetYaxis()->SetTitle("#frac{Data}{Tsallis fit}");
//            histoXSecDummy[padCounter+1]->DrawCopy();
//            if(i==0) histoXSecDummy[padCounter+1]->GetYaxis()->SetTitle("#frac{Data}{TCM fit}");
    for(Int_t i=0;i<numActiveMeas;i++){
        padXSec[i+1]->cd();
        histoXSecDummy[i+1]->GetYaxis()->SetRangeUser(0.3,1.99);
        histoXSecDummy[i+1]->GetYaxis()->SetTitle("#frac{Data}{TCM fit}");
        if(i==3) histoXSecDummy[i+1]->GetYaxis()->SetTitle("#frac{Data}{Tsallis fit}");
        histoXSecDummy[i+1]->DrawCopy();
    }

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    //__________________________________________ Draw combined spectra in first pad
    padXSec[0]->cd();
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            if(graphEtaInvariantCrossSectionSys[i][10]&&graphEtaInvariantCrossSectionStat[i][10]){
              if(i!=0) fitTCMInvCrossSectionEtaCombPlot[i]->Draw("same");
              fitTsallisInvCrossSectionEtaComb[i]->Draw("same");
              //graphEtaInvariantCrossSectionSys[i][10]->Draw("E2same");
              //graphEtaInvariantCrossSectionStat[i][10]->Draw("p,same,z");
            }
            for(Int_t iM = 0; iM<10; iM++){
              if(graphEtaInvariantCrossSectionSys[i][iM]&&graphEtaInvariantCrossSectionStat[i][iM]){
                graphEtaInvariantCrossSectionSys[i][iM]->Draw("E2same");
                graphEtaInvariantCrossSectionStat[i][iM]->Draw("p,same,z");
              }
            }
        }
    }
    drawLatexAdd("ALICE",rightalignDouble,0.92,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#eta #rightarrow #gamma#gamma",rightalignDouble,0.88,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

    TLegend* legendInvariantCrossSectionEta    = GetAndSetLegend2(0.17, 0.03, 0.5, 0.03+textsizeLabelsXSec[0]*(numActiveMeas+1)+textsizeLabelsXSec[0], textSizeLabelsPixel);
    legendInvariantCrossSectionEta->SetNColumns(1);
    legendInvariantCrossSectionEta->SetMargin(0.2);
    legendRunningIndex = numActiveMeas-1;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            legendInvariantCrossSectionEta->AddEntry(graphEtaInvariantCrossSectionSys[i][10],Form("pp, %s%s",energyLatex[i].Data(),legendScalingString[legendRunningIndex].Data()),"pf");
            legendRunningIndex-=1;
        }
    }
    legendInvariantCrossSectionEta->AddEntry(fitTCMInvCrossSectionEtaCombPlot[exampleActiveMeas],"TCM fit","l");
    legendInvariantCrossSectionEta->AddEntry(fitTsallisInvCrossSectionEtaComb[exampleActiveMeas],"Tsallis fit","l");
    legendInvariantCrossSectionEta->Draw();

    //__________________________________________ Draw ratios of individual spectra to combined fit in lower pads
    padCounter = 0;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            padXSec[padCounter+1]->cd();
            DrawGammaLines(xRangeMinXSec, xRangeMaxXSec , 1., 1., 1.2, kGray+2);
            for (Int_t j = 0; j < 10; j++){
                if(graphEtaInvariantCrossSectionSys[i][j]&&graphEtaInvariantCrossSectionStat[i][j]){
                    graphRatioIndivMeasToCombFitEtaSys[i][j]->Draw("2,same");
                    graphRatioIndivMeasToCombFitEtaStat_woXErr[i][j]->Draw("p,same");
                }
            }
            padCounter+=1;
        }
    }
    canvasXSec->Print(Form("%s/Eta_XSec_Ratios.%s",outputDir.Data(),suffix.Data()));

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    padXSec[0]->cd();
    histoXSecDummyEta[0]->Draw();
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            if(graphEtaInvariantCrossSectionSys[i][10]&&graphEtaInvariantCrossSectionStat[i][10]){
                if(graphPythiaInvCrossSectionEta[i]) graphPythiaInvCrossSectionEta[i]->Draw("3,same");
                if(histoPythiaInvCrossSectionEta[i]) histoPythiaInvCrossSectionEta[i]->Draw("same,hist,l");
                if(i!=0) fitTCMInvCrossSectionEtaCombPlot[i]->Draw("same");
                fitTsallisInvCrossSectionEtaComb[i]->Draw("same");
                graphEtaInvariantCrossSectionSys[i][10]->Draw("E2same");
                graphEtaInvariantCrossSectionStat[i][10]->Draw("p,same,z");
            }
        }
    }
    drawLatexAdd("ALICE",rightalignDouble,0.92,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#eta #rightarrow #gamma#gamma",rightalignDouble,0.88,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

    TLegend* legendInvariantXSecEta    = GetAndSetLegend2(0.17, 0.03, 0.5, 0.03+textsizeLabelsXSec[0]*numActiveMeas+4*textsizeLabelsXSec[0], textSizeLabelsPixel);
    legendInvariantXSecEta->SetNColumns(1);
    legendInvariantXSecEta->SetMargin(0.2);
    legendRunningIndex = numActiveMeas-1;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            legendInvariantXSecEta->AddEntry(graphEtaInvariantCrossSectionSys[i][10],Form("pp, %s%s",energyLatex[i].Data(),legendScalingString[legendRunningIndex].Data()),"pf");
            legendRunningIndex-=1;
        }
    }
    legendInvariantXSecEta->AddEntry(fitTCMInvCrossSectionEtaCombPlot[exampleActiveMeas],"TCM fit","l");
    legendInvariantXSecEta->AddEntry(fitTsallisInvCrossSectionEtaComb[exampleActiveMeas],"Tsallis fit","l");
    legendInvariantXSecEta->AddEntry(histoRatioPythiaToCombFitEta[exampleActiveMeas],"PYTHIA 8.2, Monash 2013","l");
    legendInvariantXSecEta->Draw();
    //__________________________________________ Loop over all ratio pads and draw
    for(Int_t i=0;i<numActiveMeas;i++){
        padXSec[i+1]->cd();
        histoXSecDummy[i+1]->GetYaxis()->SetTitle("#frac{Data, Theory}{TCM fit}");
        if(i==3) histoXSecDummy[i+1]->GetYaxis()->SetTitle("#frac{Data, Theory}{Tsallis fit}");
        histoXSecDummy[i+1]->DrawCopy();
        histoXSecDummy[i+1]->GetYaxis()->SetTitle("#frac{Data}{TCM fit}");
        if(i==3) histoXSecDummy[i+1]->GetYaxis()->SetTitle("#frac{Data}{Tsallis fit}");
    }
    padCounter = 0;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            padXSec[padCounter+1]->cd();
            DrawGammaLines(xRangeMinXSec, xRangeMaxXSec , 1., 1., 1.2, kGray+2);
            if(graphEtaInvariantCrossSectionSys[i][10]&&graphEtaInvariantCrossSectionStat[i][10]){
                graphRatioPythiaToCombFitEta[i]->Draw("3,same");
                histoRatioPythiaToCombFitEta[i]->Draw("same,hist,l");
                graphRatioIndivMeasToCombFitEtaSys[i][10]->Draw("2,same");
                graphRatioIndivMeasToCombFitEtaStat_woXErr[i][10]->Draw("p,same");
            }
            padCounter+=1;
        }
    }
    canvasXSec->Print(Form("%s/Eta_XSec_FullRatios_Pythia.%s",outputDir.Data(),suffix.Data()));

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    padXSec[0]->cd();
    histoXSecDummyEta[0]->Draw();
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            if(graphEtaInvariantCrossSectionSys[i][10]&&graphEtaInvariantCrossSectionStat[i][10]){
              if(graphNLOAESSSInvCrossSectionEtaHalf[i]) graphNLOAESSSInvCrossSectionEtaHalf[i]->Draw("same,c");
              if(graphNLOAESSSInvCrossSectionEtaOne[i]) graphNLOAESSSInvCrossSectionEtaOne[i]->Draw("same,c");
              if(graphNLOAESSSInvCrossSectionEtaTwo[i]) graphNLOAESSSInvCrossSectionEtaTwo[i]->Draw("same,c");
              if(graphPythiaInvCrossSectionEta[i]) graphPythiaInvCrossSectionEta[i]->Draw("3,same");
              if(histoPythiaInvCrossSectionEta[i]) histoPythiaInvCrossSectionEta[i]->Draw("same,hist,l");
              if(i!=0) fitTCMInvCrossSectionEtaCombPlot[i]->Draw("same");
              fitTsallisInvCrossSectionEtaComb[i]->Draw("same");
              graphEtaInvariantCrossSectionSys[i][10]->Draw("E2same");
              graphEtaInvariantCrossSectionStat[i][10]->Draw("p,same,z");
            }
        }
    }
    drawLatexAdd("ALICE",rightalignDouble,0.92,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#eta #rightarrow #gamma#gamma",rightalignDouble,0.88,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

    TLegend* legendInvariantCrossSectionEta44    = GetAndSetLegend2(0.696, 0.652, 0.94, 0.652+textsizeLabelsXSec[0]*(numActiveMeas-2)+textsizeLabelsXSec[0], 0.9*textSizeLabelsPixel);
    legendInvariantCrossSectionEta44->SetNColumns(1);
    legendInvariantCrossSectionEta44->SetMargin(0.3);
    legendInvariantCrossSectionEta44->AddEntry((TObject*)0,"","");
    legendInvariantCrossSectionEta44->AddEntry(graphNLOAESSSInvCrossSectionEtaHalf[exampleActiveMeas],"#mu = 0.5 #it{p}_{T}","l");
    legendInvariantCrossSectionEta44->AddEntry(graphNLOAESSSInvCrossSectionEtaOne[exampleActiveMeas],"#mu = #it{p}_{T}","l");
    legendInvariantCrossSectionEta44->AddEntry(graphNLOAESSSInvCrossSectionEtaTwo[exampleActiveMeas],"#mu = 2 #it{p}_{T}","l");
    legendInvariantCrossSectionEta44->Draw();
    legendInvariantXSecEta->Draw();
    drawLatexAdd("NLO,",        0.7,0.81,0.9*textsizeLabelsXSec[0],kFALSE,kFALSE,kFALSE);
    drawLatexAdd("PDF:CTEQ6M5", 0.7,0.78,0.9*textsizeLabelsXSec[0],kFALSE,kFALSE,kFALSE);
    drawLatexAdd("FF:AESSS",    0.7,0.75,0.9*textsizeLabelsXSec[0],kFALSE,kFALSE,kFALSE);
    //__________________________________________ Loop over all ratio pads and draw
    for(Int_t i=0;i<numActiveMeas;i++){
        padXSec[i+1]->cd();
        histoXSecDummy[i+1]->GetYaxis()->SetTitle("#frac{Data, Theory}{TCM fit}");
        if(i==3) histoXSecDummy[i+1]->GetYaxis()->SetTitle("#frac{Data, Theory}{Tsallis fit}");
        histoXSecDummy[i+1]->GetYaxis()->SetRangeUser(0.3,3.95);
        histoXSecDummy[i+1]->DrawCopy();
        histoXSecDummy[i+1]->GetYaxis()->SetRangeUser(0.3,1.99);
        histoXSecDummy[i+1]->GetYaxis()->SetTitle("#frac{Data}{TCM fit}");
        if(i==3) histoXSecDummy[i+1]->GetYaxis()->SetTitle("#frac{Data}{Tsallis fit}");
    }
    padCounter = 0;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            padXSec[padCounter+1]->cd();
            DrawGammaLines(xRangeMinXSec, xRangeMaxXSec , 1., 1., 1.2, kGray+2);
            if(graphEtaInvariantCrossSectionSys[i][10]&&graphEtaInvariantCrossSectionStat[i][10]){
              if(graphRatioNLOAESSSToCombFitEtaHalf[i]) graphRatioNLOAESSSToCombFitEtaHalf[i]->Draw("same,c");
              if(graphRatioNLOAESSSToCombFitEtaOne[i]) graphRatioNLOAESSSToCombFitEtaOne[i]->Draw("same,c");
              if(graphRatioNLOAESSSToCombFitEtaTwo[i]) graphRatioNLOAESSSToCombFitEtaTwo[i]->Draw("same,c");
              graphRatioPythiaToCombFitEta[i]->Draw("3,same");
              if(i==4)histoRatioPythiaToCombFitEta[i]->GetXaxis()->SetRangeUser(0.6,35.);
              histoRatioPythiaToCombFitEta[i]->DrawCopy("same,hist,l");
              if(i==4)histoRatioPythiaToCombFitEta[i]->GetXaxis()->SetRangeUser(0.8,35.);
              graphRatioIndivMeasToCombFitEtaSys[i][10]->Draw("2,same");
              graphRatioIndivMeasToCombFitEtaStat_woXErr[i][10]->Draw("p,same");
            }
            padCounter+=1;
        }
    }
    canvasXSec->Print(Form("%s/Eta_XSec_FullRatios_Pythia_TheoryAESSS.%s",outputDir.Data(),suffix.Data()));

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    padXSec[0]->cd();
    histoXSecDummyEta[0]->Draw();
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            if(graphEtaInvariantCrossSectionSys[i][10]&&graphEtaInvariantCrossSectionStat[i][10]){
                if(i!=0) fitTCMInvCrossSectionEtaCombPlot[i]->Draw("same");
                fitTsallisInvCrossSectionEtaComb[i]->Draw("same");
                graphEtaInvariantCrossSectionSys[i][10]->Draw("E2same");
                graphEtaInvariantCrossSectionStat[i][10]->Draw("p,same,z");
            }
        }
    }
    drawLatexAdd("ALICE",rightalignDouble,0.92,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#eta #rightarrow #gamma#gamma",rightalignDouble,0.88,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

    legendInvariantCrossSectionEta->Draw();
    //__________________________________________ Loop over all ratio pads and draw
    for(Int_t i=0;i<numActiveMeas;i++){
        padXSec[i+1]->cd();
        histoXSecDummy[i+1]->DrawCopy();
    }
    padCounter = 0;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            padXSec[padCounter+1]->cd();
            DrawGammaLines(xRangeMinXSec, xRangeMaxXSec , 1., 1., 1.2, kGray+2);
            if(graphEtaInvariantCrossSectionSys[i][10]&&graphEtaInvariantCrossSectionStat[i][10]){
                graphRatioIndivMeasToCombFitEtaSys[i][10]->Draw("2,same");
                graphRatioIndivMeasToCombFitEtaStat_woXErr[i][10]->Draw("p,same");
            }
            padCounter+=1;
        }
    }
    canvasXSec->Print(Form("%s/Eta_XSec_FullRatios.%s",outputDir.Data(),suffix.Data()));

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    TCanvas* canvasXSectionEta  = new TCanvas("canvasXSectionEta","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasXSectionEta, 0.14, 0.01, 0.01, 0.07);
    canvasXSectionEta->SetLogx();
    canvasXSectionEta->SetLogy();
    histoXSecDummyEta[0]->GetXaxis()->SetLabelOffset(0.0);
    histoXSecDummyEta[0]->GetXaxis()->SetTitleOffset(0.8);
    histoXSecDummyEta[0]->Draw();
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            if(graphEtaInvariantCrossSectionSys[i][10]&&graphEtaInvariantCrossSectionStat[i][10]){
              if(i!=0) fitTCMInvCrossSectionEtaCombPlot[i]->Draw("same");
              fitTsallisInvCrossSectionEtaComb[i]->Draw("same");
              graphEtaInvariantCrossSectionSys[i][10]->Draw("E2same");
              graphEtaInvariantCrossSectionStat[i][10]->Draw("p,same,z");
            }
        }
    }
    drawLatexAdd("ALICE",rightalignDouble,0.92,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#eta #rightarrow #gamma#gamma",rightalignDouble,0.88,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

    TLegend* legendInvariantCrossSectionEta2    = GetAndSetLegend2(0.17, 0.10, 0.5, 0.10+textsizeLabelsXSec[0]*(numActiveMeas+1)+textsizeLabelsXSec[0], textSizeLabelsPixel);
    legendInvariantCrossSectionEta2->SetNColumns(1);
    legendInvariantCrossSectionEta2->SetMargin(0.2);
    legendRunningIndex = numActiveMeas-1;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            legendInvariantCrossSectionEta2->AddEntry(graphEtaInvariantCrossSectionSys[i][10],Form("pp, %s%s",energyLatex[i].Data(),legendScalingString[legendRunningIndex].Data()),"pf");
            legendRunningIndex-=1;
        }
    }
    legendInvariantCrossSectionEta2->AddEntry(fitTCMInvCrossSectionEtaCombPlot[exampleActiveMeas],"TCM fit","l");
    legendInvariantCrossSectionEta2->AddEntry(fitTsallisInvCrossSectionEtaComb[exampleActiveMeas],"Tsallis fit","l");
    legendInvariantCrossSectionEta2->Draw();

    canvasXSectionEta->Print(Form("%s/Eta_XSec.%s",outputDir.Data(),suffix.Data()));

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    histoXSecDummyEta[0]->Draw();
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            if(graphEtaInvariantCrossSectionSys[i][10]&&graphEtaInvariantCrossSectionStat[i][10]){
              if(graphPythiaInvCrossSectionEta[i]) graphPythiaInvCrossSectionEta[i]->Draw("3,same");
              if(histoPythiaInvCrossSectionEta[i]) histoPythiaInvCrossSectionEta[i]->Draw("same,hist,l");
              if(i!=0) fitTCMInvCrossSectionEtaCombPlot[i]->Draw("same");
              fitTsallisInvCrossSectionEtaComb[i]->Draw("same");
              graphEtaInvariantCrossSectionSys[i][10]->Draw("E2same");
              graphEtaInvariantCrossSectionStat[i][10]->Draw("p,same,z");
            }
        }
    }
    drawLatexAdd("ALICE",rightalignDouble,0.92,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#eta #rightarrow #gamma#gamma",rightalignDouble,0.88,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

    TLegend* legendInvariantCrossSectionEta3    = GetAndSetLegend2(0.17, 0.10, 0.5, 0.10+textsizeLabelsXSec[0]*(numActiveMeas+2)+textsizeLabelsXSec[0], textSizeLabelsPixel);
    legendInvariantCrossSectionEta3->SetNColumns(1);
    legendInvariantCrossSectionEta3->SetMargin(0.2);
    legendRunningIndex = numActiveMeas-1;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            legendInvariantCrossSectionEta3->AddEntry(graphEtaInvariantCrossSectionSys[i][10],Form("pp, %s%s",energyLatex[i].Data(),legendScalingString[legendRunningIndex].Data()),"pf");
            legendRunningIndex-=1;
        }
    }
    legendInvariantCrossSectionEta3->AddEntry(fitTCMInvCrossSectionEtaCombPlot[exampleActiveMeas],"TCM fit","l");
    legendInvariantCrossSectionEta3->AddEntry(fitTsallisInvCrossSectionEtaComb[exampleActiveMeas],"Tsallis fit","l");
    legendInvariantCrossSectionEta3->AddEntry(histoRatioPythiaToCombFitEta[exampleActiveMeas],"PYTHIA 8.2, Monash 2013","l");
    legendInvariantCrossSectionEta3->Draw();

    canvasXSectionEta->Print(Form("%s/Eta_XSec_Pythia.%s",outputDir.Data(),suffix.Data()));

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    histoXSecDummyEta[0]->Draw();
    for (Int_t i = 0; i < 6; i++){
        if(includeEnergy[i]){
            if(graphEtaInvariantCrossSectionSys[i][10]&&graphEtaInvariantCrossSectionStat[i][10]){
              if(graphNLOAESSSInvCrossSectionEtaHalf[i]) graphNLOAESSSInvCrossSectionEtaHalf[i]->Draw("same,c");
              if(graphNLOAESSSInvCrossSectionEtaOne[i]) graphNLOAESSSInvCrossSectionEtaOne[i]->Draw("same,c");
              if(graphNLOAESSSInvCrossSectionEtaTwo[i]) graphNLOAESSSInvCrossSectionEtaTwo[i]->Draw("same,c");
              if(graphPythiaInvCrossSectionEta[i]) graphPythiaInvCrossSectionEta[i]->Draw("3,same");
              if(histoPythiaInvCrossSectionEta[i]) histoPythiaInvCrossSectionEta[i]->Draw("same,hist,l");
              if(i!=0) fitTCMInvCrossSectionEtaCombPlot[i]->Draw("same");
              fitTsallisInvCrossSectionEtaComb[i]->Draw("same");
              graphEtaInvariantCrossSectionSys[i][10]->Draw("E2same");
              graphEtaInvariantCrossSectionStat[i][10]->Draw("p,same,z");
            }
        }
    }
    drawLatexAdd("ALICE",rightalignDouble,0.92,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
    drawLatexAdd("#eta #rightarrow #gamma#gamma",rightalignDouble,0.88,textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

    legendInvariantCrossSectionEta3->Draw();
    TLegend* legendInvariantCrossSectionEta4    = GetAndSetLegend2(0.696, 0.652, 0.94, 0.652+textsizeLabelsXSec[0]*(numActiveMeas-2)+textsizeLabelsXSec[0], 0.9*textSizeLabelsPixel);
    legendInvariantCrossSectionEta4->SetNColumns(1);
    legendInvariantCrossSectionEta4->SetMargin(0.2);
    legendInvariantCrossSectionEta4->AddEntry((TObject*)0,"","");
    legendInvariantCrossSectionEta4->AddEntry(graphNLOAESSSInvCrossSectionEtaHalf[exampleActiveMeas],"#mu = 0.5 #it{p}_{T}","l");
    legendInvariantCrossSectionEta4->AddEntry(graphNLOAESSSInvCrossSectionEtaOne[exampleActiveMeas],"#mu = #it{p}_{T}","l");
    legendInvariantCrossSectionEta4->AddEntry(graphNLOAESSSInvCrossSectionEtaTwo[exampleActiveMeas],"#mu = 2 #it{p}_{T}","l");
    legendInvariantCrossSectionEta4->Draw();

    drawLatexAdd("NLO,",        0.7,0.81,0.9*textsizeLabelsXSec[0],kFALSE,kFALSE,kFALSE);
    drawLatexAdd("PDF:CTEQ6M5", 0.7,0.78,0.9*textsizeLabelsXSec[0],kFALSE,kFALSE,kFALSE);
    drawLatexAdd("FF:AESSS",    0.7,0.75,0.9*textsizeLabelsXSec[0],kFALSE,kFALSE,kFALSE);

    canvasXSectionEta->Print(Form("%s/Eta_XSec_Pythia_TheoryAESSS.%s",outputDir.Data(),suffix.Data()));


    canvasRatios->cd();
    //__________________________________________ Set ratio pad properties
    for(Int_t i=0;i<numActiveMeas;i++){
        if(i==numActiveMeas-1)
            DrawGammaPadSettings( padRatios[i], relativeMarginsXRatios[0], relativeMarginsXRatios[2], relativeMarginsYRatios[1], relativeMarginsYRatios[2]);
        else
            DrawGammaPadSettings( padRatios[i], relativeMarginsXRatios[0], relativeMarginsXRatios[2], relativeMarginsYRatios[1], relativeMarginsYRatios[1]);
        padRatios[i]->Draw();
    }

    //__________________________________________ Loop over all ratio pads and draw
    for(Int_t i=0;i<numActiveMeas;i++){
        padRatios[i]->cd();
        padRatios[i]->SetLogx(1);
        histoRatiosDummy[i]                     = new TH2F(Form("histoRatiosDummy_%d",i),Form("histoRatiosDummy_%d",i),1000,xRangeMinRatios,xRangeMaxRatios,1000,0.1,4.0);
        SetStyleHistoTH2ForGraphs(histoRatiosDummy[i], "#it{p}_{T} (GeV/#it{c})","#frac{Data}{TCM fit}", 0.85*textsizeLabelsRatios[i], textsizeLabelsRatios[i],
                                  0.85*textsizeLabelsRatios[i],textsizeLabelsRatios[i], 1,0.2/(textsizeFacRatios[i]*marginRatios), 510, 505);
        if(i==3) histoRatiosDummy[i]->GetYaxis()->SetTitle("#frac{Data}{Tsallis fit}");
        histoRatiosDummy[i]->GetYaxis()->SetRangeUser(0.3,1.99);
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

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    padCounter = 0;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            padRatios[padCounter]->cd();
            DrawGammaLines(xRangeMinRatios, xRangeMaxRatios, 1., 1., 1.2, kGray+2);
            if(padCounter==0){
              drawLatexAdd("ALICE",0.27,0.84,3*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
              drawLatexAdd("#eta #rightarrow #gamma#gamma",0.28,0.73,3*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
            }
            for (Int_t j = 0; j < 10; j++){
                if(graphEtaInvariantCrossSectionSys[i][j]&&graphEtaInvariantCrossSectionStat[i][j]){
                    graphRatioIndivMeasToCombFitEtaSys[i][j]->Draw("2,same");
                    graphRatioIndivMeasToCombFitEtaStat_woXErr[i][j]->Draw("p,same");
                }
            }
            padCounter+=1;
        }
    }
    for(Int_t i=0;i<numActiveMeas;i++) padRatios[i]->RedrawAxis();
    canvasRatios->Print(Form("%s/Eta_Ratios.%s",outputDir.Data(),suffix.Data()));

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

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
            DrawGammaLines(xRangeMinRatios, xRangeMaxRatios, 1., 1., 1.2, kGray+2);
            if(padCounter==0){
              drawLatexAdd("ALICE",rightalignDouble,0.84,3*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
              drawLatexAdd("#eta #rightarrow #gamma#gamma",rightalignDouble,0.73,3*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
            }
            if(graphEtaInvariantCrossSectionSys[i][10]&&graphEtaInvariantCrossSectionStat[i][10]){
                if(i==0) legendRatiosEnergy[i] = GetAndSetLegend2(0.15, 0.81, 0.5, 0.76+textsizeLabelsXSec[0]*(numActiveMeas+1)+textsizeLabelsXSec[0], textSizeLabelsPixel);
                else legendRatiosEnergy[i] = GetAndSetLegend2(0.15, 0.75, 0.5, 0.75+textsizeLabelsXSec[0]*(numActiveMeas+1)+textsizeLabelsXSec[0], textSizeLabelsPixel);
                legendRatiosEnergy[i]->SetNColumns(1);
                legendRatiosEnergy[i]->SetMargin(0.2);
                legendRatiosEnergy[i]->AddEntry(graphEtaInvariantCrossSectionSys[i][10],Form("pp, %s",energyLatex[i].Data()),"pf");
                legendRatiosEnergy[i]->Draw();
                graphRatioIndivMeasToCombFitEtaSys[i][10]->Draw("2,same");
                graphRatioIndivMeasToCombFitEtaStat_woXErr[i][10]->Draw("p,same");
            }
            padCounter+=1;
        }
    }
    for(Int_t i=0;i<numActiveMeas;i++) padRatios[i]->RedrawAxis();
    canvasRatios->Print(Form("%s/Eta_FullRatios.%s",outputDir.Data(),suffix.Data()));

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    canvasRatios->cd();
    //__________________________________________ Set ratio pad properties
    for(Int_t i=0;i<numActiveMeas;i++){
        padRatios[i]->Draw();
    }
    //__________________________________________ Loop over all ratio pads and draw
    for(Int_t i=0;i<numActiveMeas;i++){
        padRatios[i]->cd();
        padRatios[i]->SetLogx(1);
        histoRatiosDummy[i]->GetYaxis()->SetTitle("#frac{Data, Theory}{TCM fit}");
        histoRatiosDummy[i]->DrawCopy();
    }
    padCounter = 0;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            padRatios[padCounter]->cd();
            DrawGammaLines(xRangeMinRatios, xRangeMaxRatios, 1., 1., 1.2, kGray+2);
            if(padCounter==0){
              drawLatexAdd("ALICE",0.27,0.84,3*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
              drawLatexAdd("#eta #rightarrow #gamma#gamma",0.28,0.73,3*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
            }else if(padCounter==3){
              legendRatiosPythia = GetAndSetLegend2(0.15, 0.735, 0.4, 0.735+textsizeLabelsXSec[0]*(numActiveMeas)+textsizeLabelsXSec[0], textSizeLabelsPixel);
              legendRatiosPythia->SetNColumns(1);
              legendRatiosPythia->SetMargin(0.2);
              legendRatiosPythia->AddEntry(tempPythia,"PYTHIA 8.2,","l");
              legendRatiosPythia->AddEntry((TObject*)0,"Monash 2013","");
              legendRatiosPythia->Draw();
            }
            if(graphEtaInvariantCrossSectionSys[i][10]&&graphEtaInvariantCrossSectionStat[i][10]){
              graphRatioPythiaToCombFitEta[i]->Draw("3,same");
              histoRatioPythiaToCombFitEta[i]->Draw("same,hist,l");
              if(i==0) legendRatiosEnergy[i] = GetAndSetLegend2(0.6, 0.81, 0.95, 0.76+textsizeLabelsXSec[0]*(numActiveMeas+1)+textsizeLabelsXSec[0], textSizeLabelsPixel);
              else legendRatiosEnergy[i] = GetAndSetLegend2(0.6, 0.75, 0.95, 0.75+textsizeLabelsXSec[0]*(numActiveMeas+1)+textsizeLabelsXSec[0], textSizeLabelsPixel);
              legendRatiosEnergy[i]->SetNColumns(1);
              legendRatiosEnergy[i]->SetMargin(0.2);
              legendRatiosEnergy[i]->AddEntry(graphEtaInvariantCrossSectionSys[i][10],Form("pp, %s",energyLatex[i].Data()),"pf");
              legendRatiosEnergy[i]->Draw();
              graphRatioIndivMeasToCombFitEtaSys[i][10]->Draw("2,same");
              graphRatioIndivMeasToCombFitEtaStat_woXErr[i][10]->Draw("p,same");
            }
            padCounter+=1;
        }
    }
    for(Int_t i=0;i<numActiveMeas;i++) padRatios[i]->RedrawAxis();
    canvasRatios->Print(Form("%s/Eta_FullRatios_Pythia.%s",outputDir.Data(),suffix.Data()));

    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------Plotting -----------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------

    TLegend* legendInvariantXSecTheoryEta = NULL;
    canvasRatios->cd();
    //__________________________________________ Set ratio pad properties
    for(Int_t i=0;i<numActiveMeas;i++){
        padRatios[i]->Draw();
    }
    //__________________________________________ Loop over all ratio pads and draw
    for(Int_t i=0;i<numActiveMeas;i++){
        padRatios[i]->cd();
        padRatios[i]->SetLogx(1);
        histoRatiosDummy[i]->GetYaxis()->SetRangeUser(0.4,3.95);
        histoRatiosDummy[i]->DrawCopy();
        histoRatiosDummy[i]->GetYaxis()->SetRangeUser(0.4,1.99);
    }
    padCounter = 0;
    for (Int_t i = 5; i > -1; i--){
        if(includeEnergy[i]){
            padRatios[padCounter]->cd();
            DrawGammaLines(xRangeMinRatios, xRangeMaxRatios, 1., 1., 1.2, kGray+2);
            if(padCounter==0){
              drawLatexAdd("ALICE",rightalignDouble,0.84,3*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
              drawLatexAdd("#eta #rightarrow #gamma#gamma",rightalignDouble,0.73,3*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
            }else if(padCounter==3){
              legendRatiosPythia = GetAndSetLegend2(0.15, 0.605, 0.5, 0.605+textsizeLabelsXSec[0]*(numActiveMeas)+textsizeLabelsXSec[0], textSizeLabelsPixel);
              legendRatiosPythia->SetNColumns(1);
              legendRatiosPythia->SetMargin(0.2);
              legendRatiosPythia->AddEntry(tempPythia,"PYTHIA 8.2,","l");
              legendRatiosPythia->AddEntry((TObject*)0,"Monash 2013","");
              legendRatiosPythia->Draw();

              legendInvariantXSecTheoryEta = GetAndSetLegend2(0.75, 0.595, 0.95, 0.595+textsizeLabelsXSec[0]*(numActiveMeas+3)+textsizeLabelsXSec[0], textSizeLabelsPixel);
              legendInvariantXSecTheoryEta->SetNColumns(1);
              legendInvariantXSecTheoryEta->SetMargin(0.2);
              legendInvariantXSecTheoryEta->AddEntry((TObject*)0,"","");
              legendInvariantXSecTheoryEta->AddEntry(graphNLOAESSSInvCrossSectionEtaHalf[exampleActiveMeas],"#mu = 0.5 #it{p}_{T}","l");
              legendInvariantXSecTheoryEta->AddEntry(graphNLOAESSSInvCrossSectionEtaOne[exampleActiveMeas],"#mu = #it{p}_{T}","l");
              legendInvariantXSecTheoryEta->AddEntry(graphNLOAESSSInvCrossSectionEtaTwo[exampleActiveMeas],"#mu = 2 #it{p}_{T}","l");
              legendInvariantXSecTheoryEta->Draw();
              drawLatexAdd("NLO,",0.83,0.86,2.5*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
              drawLatexAdd("PDF:CTEQ6M5",0.968,0.51,2.2*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);
              drawLatexAdd("FF:AESSS",0.9,0.44,2.2*textsizeLabelsXSec[0],kFALSE,kFALSE,kTRUE);

            }
            if(graphEtaInvariantCrossSectionSys[i][10]&&graphEtaInvariantCrossSectionStat[i][10]){
              if(graphRatioNLOAESSSToCombFitEtaHalf[i]) graphRatioNLOAESSSToCombFitEtaHalf[i]->Draw("same,c");
              if(graphRatioNLOAESSSToCombFitEtaOne[i]) graphRatioNLOAESSSToCombFitEtaOne[i]->Draw("same,c");
              if(graphRatioNLOAESSSToCombFitEtaTwo[i]) graphRatioNLOAESSSToCombFitEtaTwo[i]->Draw("same,c");
              graphRatioPythiaToCombFitEta[i]->Draw("3,same");
              if(i==4)histoRatioPythiaToCombFitEta[i]->GetXaxis()->SetRangeUser(0.6,35.);
              histoRatioPythiaToCombFitEta[i]->DrawCopy("same,hist,l");
              if(i==4)histoRatioPythiaToCombFitEta[i]->GetXaxis()->SetRangeUser(0.8,35.);
              graphRatioIndivMeasToCombFitEtaSys[i][10]->Draw("2,same");
              graphRatioIndivMeasToCombFitEtaStat_woXErr[i][10]->Draw("p,same");
              if(i==0) legendRatiosEnergy[i] = GetAndSetLegend2(0.15, 0.81, 0.5, 0.76+textsizeLabelsXSec[0]*(numActiveMeas+1)+textsizeLabelsXSec[0], textSizeLabelsPixel);
              else legendRatiosEnergy[i] = GetAndSetLegend2(0.15, 0.75, 0.5, 0.75+textsizeLabelsXSec[0]*(numActiveMeas+1)+textsizeLabelsXSec[0], textSizeLabelsPixel);
              legendRatiosEnergy[i]->SetNColumns(1);
              legendRatiosEnergy[i]->SetMargin(0.2);
              legendRatiosEnergy[i]->AddEntry(graphEtaInvariantCrossSectionSys[i][10],Form("pp, %s",energyLatex[i].Data()),"pf");
              legendRatiosEnergy[i]->Draw();
            }
            padCounter+=1;
        }
    }
    for(Int_t i=0;i<numActiveMeas;i++) padRatios[i]->RedrawAxis();
    canvasRatios->Print(Form("%s/Eta_FullRatios_Pythia_AESSS.%s",outputDir.Data(),suffix.Data()));
}
