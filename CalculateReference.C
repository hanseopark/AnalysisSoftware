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
#include "TProfile.h"
#include "TRandom3.h"
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
#include "TSpline.h"
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

TGraphAsymmErrors *GetInterpolSpectrum2D(Int_t nDataPoints, TGraphAsymmErrors** graphs, Double_t* energies, Double_t dSqrts, Int_t nIterations, Int_t modus );
TH1D* ConvertYieldHisto(TH1D* input, Bool_t DivideBy2pi, Bool_t DivideByPt, Bool_t MultiplyBy2pi, Bool_t MultiplyByPt);
TF1* DoFitWithTsallis(TGraph* graph, TString name, TString particle, Double_t p0, Double_t p1, Double_t p2);
TF1* DoFitWithTCM(TGraph* graph, TString name, TString particle, Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4 );
void PlotInterpolationPtBins(TGraphErrors** gPtvSqrts,TGraphErrors** gPtvsEnergies, TF1** fPowerlaw, TF1** fPowerlawMC, TGraphAsymmErrors* gRpPb,Int_t fColumnPlot, Int_t fRowPlot,TString namePlot);
void PlotInterpolationSinglePtBin(TGraphErrors* gPtvSqrts,TGraphErrors* gPtvsEnergies, TF1* fPowerlaw, TF1* fPowerlawMC, TGraphAsymmErrors* gRpPb, Int_t ptBin, TString namePlot);
void PlotAlphavsPt(TGraphErrors* gAlpha, TGraphErrors* gAlphaSyst, TGraphErrors* gAlphaMC, TSpline* splineMC, TString method, TString thesisPlotLabel, TString namePlot);
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

TGraphErrors* graphAlpha                = 0x0;
TGraphErrors** graphPtvsSqrts           = 0x0;
TGraphErrors** gPtvsEnergiesSystem      = 0x0;
TF1** fPowerlawSystem                   = 0x0;

TGraphErrors* graphAlphaMC              = 0x0;
TGraphErrors** graphPtvsSqrtsMC         = 0x0;
TGraphErrors** gPtvsEnergiesSystemMC    = 0x0;
TF1** fPowerlawSystemMC                 = 0x0;
TSpline5* splineAlphaMC                 = 0x0;
TGraphErrors* graphDiffMCIntAndReal     = 0x0;

TGraphErrors* graphAlphaStat            = 0x0;
TGraphErrors** graphPtvsSqrtsStat       = 0x0;
TGraphErrors** gPtvsEnergiesSystemStat  = 0x0;
TF1** fPowerlawSystemStat               = 0x0;

TGraphErrors* graphAlphaSyst            = 0x0;
TGraphAsymmErrors* graphInterpolSyst    = 0x0;
TGraphErrors** graphPtvsSqrtsSyst       = 0x0;
TGraphErrors** gPtvsEnergiesSystemSyst  = 0x0;
TF1** fPowerlawSystemSyst               = 0x0;
TF1** fPowerlawSystemSystMC             = 0x0;


void CalculateStatPlusSysErrors(TH1D* histStat, TH1D* histSys);
TGraphAsymmErrors* CalculateStatPlusSysErrors(TGraphAsymmErrors* graphStat, TGraphAsymmErrors* graphSys);

//________________________________________________________________________________________________________________________
void CalculateReference (   TString configFile                  = "",
                            Int_t nDataSets                     = 2,
                            TString meson                       = "Pi0",
                            TString suffix                      = "eps",
                            TString finalEnergy                 = "8TeV",
                            Double_t energy                     = 8000,
                            Int_t mode                          = 20,
                            Bool_t doSpecialBinning             = kFALSE,
                            TString binningEnergy               = "",
                            TString configfileExclusionErrors   = "",
                            TString nameFileMCRef               = "",
                            TString nameHistMCRef               = "",
                            Int_t nTrials                       = 1000
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
    TString modeName        = ReturnTextReconstructionProcess(mode);
    TString mesonString     = ReturnMesonString (meson);
    TString collisionSystem = ReturnFullCollisionsSystem(finalEnergy);
    TString collSystOut      = ReturnCollisionEnergyOutputString(finalEnergy);
    TString detProcess      = ReturnFullTextReconstructionProcess(mode);
    TString outputDir       = Form("%s/%s/CalculateReference",suffix.Data(),dateForOutput.Data());
    TString outputDirPlots  = Form("%s/%s/CalculateReference/%s",suffix.Data(), dateForOutput.Data(), modeName.Data());
    gSystem->Exec("mkdir -p "+outputDirPlots);
    Int_t exampleBin        = 7; 
    Double_t dummyScaleFac  = 1;
    
    //*************************************************************************************************
    //***************************** read configurarion file *******************************************
    //*************************************************************************************************
    TString nameFile[5]             = {"", "", "", "", ""};
    TString nameHist[4][5]          = {{"", "", "", "", ""},{"", "", "", "", ""},{"", "", "", "", ""},{"", "", "", "", ""}};
    TString nameFileSys[5]          = {"", "", "", "", ""};
    TString nameFileMC[5]           = {"", "", "", "", ""};
    TString nameHistMC[5]           = {"", "", "", "", ""};
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
        >> energyIndName[nEnergyRead]  >> energyInd[nEnergyRead]  >> scaleFactor[nEnergyRead] >> nameHist[2][nEnergyRead] >> nameHist[3][nEnergyRead] >> nameFileSys[nEnergyRead] >> nameFileMC[nEnergyRead] >> nameHistMC[nEnergyRead];
        cout << nameFile[nEnergyRead] << "\t" << isHist[0][nEnergyRead] << "\t" << nameHist[0][nEnergyRead].Data() << "\t" << isHist[1][nEnergyRead] << "\t" << nameHist[1][nEnergyRead].Data()
        << "\t" << energyIndName[nEnergyRead].Data()  << "\t" << energyInd[nEnergyRead]  << "\t" << scaleFactor[nEnergyRead] << endl;
        nEnergyRead++;
    }

    Double_t finalBinningPt[100];
    Int_t maxNBins                  = 0;
    Int_t startBin                  = 0;
    if (doSpecialBinning){
        maxNBins                    = GetBinning( finalBinningPt, meson.Data(), binningEnergy.Data(), mode );
        startBin                    = GetStartBin( meson.Data(), binningEnergy.Data(), mode );
//         exampleBin                  = ReturnSingleInvariantMassBinPlotting (meson.Data(), binningEnergy.Data(), mode, 0, dummyScaleFac);
        if (maxNBins == 0){
            cout << "The requested binning doesn't exist, aborting!" << endl;
            return;
        }
        cout << "Binning" << endl;        
        for (Int_t i = 0; i<maxNBins; i++){
            cout << i << "\t"<< finalBinningPt[i] << "-" << finalBinningPt[i+1] <<endl;
        }
    }
    
    vector<TString>* ptSysRemNames  = new vector<TString>[5];
    Bool_t enableSysExcl[5]         = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
    Int_t nTakeOutSys               = 0;
    TString currentString           = "";
    Int_t nCurrColumn               = 0;
    
    if (configfileExclusionErrors.CompareTo("") != 0){
        ifstream in2(configfileExclusionErrors.Data());
        while(!in2.eof() && nTakeOutSys<nDataSets ){
//             cout << "read line:" <<  nTakeOutSys+1 << endl;
            in2 >> currentString;
//             cout << currentString << endl;
            if (currentString.CompareTo("STOP") == 0){
                cout << "found " << ptSysRemNames[nTakeOutSys].size() << " sys errors to be excluded." << endl;
                nTakeOutSys++;
                nCurrColumn     = 0;
            } else {
                if (nCurrColumn == 0){
                    cout << "reading energy: " << currentString.Data() << endl;
                    if (currentString.CompareTo(energyIndName[nTakeOutSys].Data()) != 0){ 
                        cout << "wrong order in file, expecting: " <<  energyIndName[nTakeOutSys].Data() << endl;
                        return;
                    }
                    nCurrColumn++;    
                } else {
                    ptSysRemNames[nTakeOutSys].push_back(currentString);
                    enableSysExcl[nTakeOutSys] = kTRUE;
                    nCurrColumn++;
                }    
            }    
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
    TH1D* histRelSystUncorr[5]          = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histRelSystCorr[5]            = {NULL, NULL, NULL, NULL, NULL};
    TH1D* histRelCombUncorr[5]          = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphStat[5]     = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphSyst[5]     = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphComb[5]     = {NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphStatReb[6]  = {NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombReb[6]  = {NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphSystReb[6]  = {NULL, NULL, NULL, NULL, NULL, NULL};
    TF1* fitComb[6]                     = {NULL, NULL, NULL, NULL, NULL, NULL};
    Bool_t haveDetailedSys[5]           = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
    Int_t nPtBinsSysAvailSingle[5]      = {-1, -1, -1, -1, -1};
    Int_t numberDetSysAvailSingle[5]    = {-1, -1, -1, -1, -1};
    vector<TString>** ptSysDetail       = new vector<TString>*[5];
    vector<Double_t>** ptSysSplit       = new vector<Double_t>*[5];
    TFile* inputFileMC[5]               = {NULL, NULL, NULL, NULL, NULL};
    Bool_t haveMC[5]                    = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};     
    TH1D* histMC[6]                     = {NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphMC[6]       = {NULL, NULL, NULL, NULL, NULL, NULL};
    TF1* fitCombMC[6]                   = {NULL, NULL, NULL, NULL, NULL, NULL};
    
    vector<Bool_t>* isTaken             = new vector<Bool_t>[5];
    for(Int_t iR=0; iR<5; iR++) {
        ptSysDetail[iR] = new vector<TString>[50];
        ptSysSplit[iR]  = new vector<Double_t>[50];
    }    
    Double_t minY                       = 1e100;
    Double_t maxY                       = 0;

    for (Int_t i = 0; i < nDataSets; i++){
        // read external detailed systematics file
        ifstream fileSysErrDetailed;
        fileSysErrDetailed.open(nameFileSys[i],ios_base::in);
        // check if the file exists
        if(fileSysErrDetailed.is_open()) {
            haveDetailedSys[i] = kTRUE;
            gSystem->Exec(Form("cp %s %s/SystematicErrorAveragedSingle_%s%s_%s.txt", nameFileSys[i].Data(), outputDir.Data(),meson.Data(),modeName.Data(), energyIndName[i].Data()));
        } else{
            haveDetailedSys[i] = kFALSE; 
            cout << "No single errors were found" << endl;
        }
            
        // read detailed file    
        if (haveDetailedSys[i]){
            cout << nameFileSys[i] << endl;
            Int_t iPtBin = 0;
            string line;
            Int_t iPtBinColumn = 0;
            // format: ptSysDetail[$NEnergyRead][$PtBin].at($DetailedSys)
            //         $PtBin == 0: names of variations
            //         $DetailedSys == 0 && $PtBin != 0: pt-value
            while (getline(fileSysErrDetailed, line) && iPtBin < 100) {
                istringstream ss(line);
                TString temp        ="";
                iPtBinColumn        = 0;
                while(ss && iPtBinColumn < 100){
                    ss >> temp;
                    if( !(iPtBin==0 && temp.CompareTo("bin")==0) && !temp.IsNull()){
                        ptSysDetail[i][iPtBin].push_back(temp);
                        iPtBinColumn++;
                    }
                }
                if(iPtBin == 0){
                    if (temp.CompareTo("TotalError") != 0){
                        cout << "adding name of for TotalError" << endl;
                        ptSysDetail[i][iPtBin++].push_back("TotalError");
                        iPtBinColumn++;
                    } else {
                        cout << "Nothing to be added" << endl;
                    }
                }else iPtBin++;
            }
            nPtBinsSysAvailSingle[i]        = iPtBin-1;
            fileSysErrDetailed.close();
            
            numberDetSysAvailSingle[i]      = (Int_t)ptSysDetail[i][0].size()-1;
            cout <<  numberDetSysAvailSingle[i] << " of individual errors:"<< endl;
            for (Int_t k = 1; k < numberDetSysAvailSingle[i]+1; k++ ){
                cout << ((TString)ptSysDetail[i][0].at(k)).Data() << "\t";
            }
            cout << endl;
            for (Int_t counter = 0; counter < nPtBinsSysAvailSingle[i]; counter++){   // pt loop
                for (Int_t k = 1; k < numberDetSysAvailSingle[i]+1; k++ ){            // sources loop 
                    cout << ((TString)ptSysDetail[i][counter+1].at(k)) << "\t" ;
                }
                cout <<  endl; 
            }
            for (Int_t k = 1; k < numberDetSysAvailSingle[i]; k++){
                Bool_t isTakenCurrent       = kTRUE;
                cout << (TString)ptSysDetail[i][0].at(k) << endl;
                if (enableSysExcl[i]){
                    for (Int_t b = 0; b < (Int_t)ptSysRemNames[i].size(); b++){
                        if (((TString)ptSysDetail[i][0].at(k)).CompareTo(ptSysRemNames[i].at(b)) == 0){
                            cout << "disabled sys: " << ptSysRemNames[i].at(b) << endl;
                            isTakenCurrent = kFALSE;
                        }
                    }
                }    
                isTaken[i].push_back(isTakenCurrent);    
            }    
            cout << "matrix for sys err disabling" << endl;
            for (Int_t k = 0; k < numberDetSysAvailSingle[i]-1; k++){
                cout << isTaken[i].at(k) << "\t" ;
            }    
            cout << endl;
            
            for (Int_t counter = 1; counter < nPtBinsSysAvailSingle[i]+1; counter++){   // pt loop
                Double_t correlatedError    = 0;
                Double_t uncorrelatedError  = 0;
                for (Int_t k = 1; k < numberDetSysAvailSingle[i]; k++){
//                     cout << ((TString)ptSysDetail[i][0].at(k)) << endl;
                    Double_t currentError   = ((TString)ptSysDetail[i][counter].at(k)).Atof();
                    if (!isTaken[i].at(k-1)){
//                         cout << "correlated" << endl;
                        correlatedError     = correlatedError+ currentError*currentError;
                    } else { 
                        uncorrelatedError   = uncorrelatedError + currentError*currentError;
                    }    
                }
                
                correlatedError             = TMath::Sqrt(correlatedError);
                uncorrelatedError           = TMath::Sqrt(uncorrelatedError);
                Double_t ptCurr             = ((TString)ptSysDetail[i][counter].at(0)).Atof();
                Double_t totalError         = ((TString)ptSysDetail[i][counter].at(numberDetSysAvailSingle[i])).Atof();

                ptSysSplit[i][counter-1].push_back(ptCurr);
                ptSysSplit[i][counter-1].push_back(correlatedError);
                ptSysSplit[i][counter-1].push_back(uncorrelatedError);
                ptSysSplit[i][counter-1].push_back(totalError);
                
                cout << ptCurr<< "\t"<< correlatedError << "\t" << uncorrelatedError << "\t" << totalError << endl;
                cout << ptSysSplit[i][counter-1].at(0) << "\t" << ptSysSplit[i][counter-1].at(1) << "\t"<< ptSysSplit[i][counter-1].at(2)<< "\t" << ptSysSplit[i][counter-1].at(3)<< endl;
                
            }    
        }
        
        // MC file reading if exisiting
        if (nameFileMC[i].CompareTo("bla") != 0){
            inputFileMC[i]              = new TFile(nameFileMC[i]);
            if (!inputFileMC[i]->IsZombie() && nameHistMC[i].CompareTo("bla") != 0){
                histMC[i]               = (TH1D*)inputFileMC[i]->Get(nameHistMC[i].Data());
                if (histMC[i]){
                    haveMC[i]           = kTRUE;
                    graphMC[i]          = new TGraphAsymmErrors(histMC[i]);
                    while (graphMC[i]->GetX()[graphMC[i]->GetN()-1] > 30) graphMC[i]->RemovePoint(graphMC[i]->GetN()-1);
                    while (graphMC[i]->GetX()[0] < 0.1) graphMC[i]->RemovePoint(0);
                    fitCombMC[i]        = DoFitWithTsallis(graphMC[i],Form("fitCombMC_%d",i),meson.Data(), graphMC[i]->GetY()[0],7.,0.2);
                }
            }
        }
        
        // Data file reading
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
            histRelSystUncorr[i]        = (TH1D*)histSyst[i]->Clone(Form("Rel%s_Uncorrelated",nameHist[1][i].Data()));
            histRelSystCorr[i]          = (TH1D*)histSyst[i]->Clone(Form("Rel%s_Correlated",nameHist[1][i].Data()));
            histRelCombUncorr[i]        = (TH1D*)histComb[i]->Clone(Form("Rel_histComb_%d_Uncorrelated",i));
            for (Int_t j = 1; j< histRelStat[i]->GetNbinsX()+1;j++){
                histRelStat[i]->SetBinContent(j,histRelStat[i]->GetBinError(j)/histRelStat[i]->GetBinContent(j));
                histRelStat[i]->SetBinError(j, 0);
                histRelSyst[i]->SetBinContent(j, histRelSyst[i]->GetBinError(j)/histRelSyst[i]->GetBinContent(j));
                histRelSyst[i]->SetBinError(j, 0);
                histRelComb[i]->SetBinContent(j, histRelComb[i]->GetBinError(j)/histRelComb[i]->GetBinContent(j));
                histRelComb[i]->SetBinError(j, 0);
                histRelSystUncorr[i]->SetBinContent(j, histRelSystUncorr[i]->GetBinError(j)/histRelSystUncorr[i]->GetBinContent(j));
                histRelSystUncorr[i]->SetBinError(j, 0);
                histRelSystCorr[i]->SetBinContent(j, 0);
                histRelSystCorr[i]->SetBinError(j, 0);
                histRelCombUncorr[i]->SetBinContent(j, histRelCombUncorr[i]->GetBinError(j)/histRelCombUncorr[i]->GetBinContent(j));
                histRelCombUncorr[i]->SetBinError(j, 0);                
            }    
            if (haveDetailedSys[i] && enableSysExcl[i]){
                for (Int_t iPt = 0; iPt < nPtBinsSysAvailSingle[i]; iPt++){
                    Double_t ptCurr     = ptSysSplit[i][iPt].at(0);
                    Int_t binPt         = histRelSyst[i]->FindBin(ptCurr);
                    Double_t stat       = histRelStat[i]->GetBinContent(binPt);
                    Double_t newSys     = histRelSystUncorr[i]->GetBinContent(binPt)*ptSysSplit[i][iPt].at(2)/ptSysSplit[i][iPt].at(3);
                    Double_t newComb    = TMath::Sqrt(newSys*newSys + stat*stat);
                    histRelSystUncorr[i]->SetBinContent(binPt, newSys);
                    histRelCombUncorr[i]->SetBinContent(binPt, newComb);
                    cout << histRelSyst[i]->GetBinContent(binPt) << "\t" <<histRelComb[i]->GetBinContent(binPt) << "\t" <<  histRelSystUncorr[i]->GetBinContent(binPt) << "\t" << histRelCombUncorr[i]->GetBinContent(binPt) << endl;
                    histRelSystCorr[i]->SetBinContent(binPt, ptSysSplit[i][iPt].at(1)/100.);
                    histRelSystCorr[i]->SetBinError(binPt, 0);
                }    
            }    
        } else if (!isHist[0][i] && !isHist[1][i]){
            graphStat[i]                = (TGraphAsymmErrors*)inputFile[i]->Get(nameHist[0][i].Data());
            ScaleGraph(graphStat[i],scaleFactor[i]);
            cout << "statistical graph" << endl;
            graphStat[i]->Print();
            
            graphSyst[i]                = (TGraphAsymmErrors*)inputFile[i]->Get(nameHist[1][i].Data());
            ScaleGraph(graphSyst[i],scaleFactor[i]);
            cout << "systematics graph" << endl;
            graphSyst[i]->Print();
            
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
            histRelSystCorr[i]          = new TH1D(Form("Rel%s_Correlated",nameHist[0][i].Data()),Form("Rel%s_Uncorrelated",nameHist[0][i].Data()),graphStat[i]->GetN()+1, binningGraphs );
            histRelSystUncorr[i]        = new TH1D(Form("Rel%s_Uncorrelated",nameHist[0][i].Data()),Form("Rel%s_Uncorrelated",nameHist[0][i].Data()),graphStat[i]->GetN()+1, binningGraphs );
            histRelCombUncorr[i]        = new TH1D(Form("Rel_histComb_%d_Uncorrelated",i), Form("Rel_histComb_%d_Uncorrelated",i), graphStat[i]->GetN()+1, binningGraphs );
            
            for (Int_t j = 0; j< graphStat[i]->GetN(); j++){
                Int_t binStat = histRelStat[i]-> FindBin(graphStat[i]->GetX()[j]);
                if (graphStat[i]->GetY()[j] != 0){
                    histRelStat[i]->SetBinContent(binStat, (graphStat[i]->GetEYhigh()[j]+graphStat[i]->GetEYlow()[j])/(2*graphStat[i]->GetY()[j]));
                    histRelStat[i]->SetBinError(binStat, 0.);
                } else {
                    histRelStat[i]->SetBinContent(binStat, 0);
                    histRelStat[i]->SetBinError(binStat, 0.);
                }    
            }
            for (Int_t j = 0; j< graphSyst[i]->GetN(); j++){
                Int_t binSyst = histRelSyst[i]-> FindBin(graphSyst[i]->GetX()[j]);
                if (graphSyst[i]->GetY()[j] != 0){
                    histRelSyst[i]->SetBinContent(binSyst, (graphSyst[i]->GetEYhigh()[j]+graphSyst[i]->GetEYlow()[j])/(2*graphSyst[i]->GetY()[j]));
                    histRelSyst[i]->SetBinError(binSyst, 0.);
                    histRelSystUncorr[i]->SetBinContent(binSyst, (graphSyst[i]->GetEYhigh()[j]+graphSyst[i]->GetEYlow()[j])/(2*graphSyst[i]->GetY()[j]));
                    histRelSystUncorr[i]->SetBinError(binSyst, 0.);
                    histRelSystCorr[i]->SetBinContent(binSyst, 0.);
                    histRelSystCorr[i]->SetBinError(binSyst, 0.);
                } else {
                    histRelSyst[i]->SetBinContent(binSyst, 0);
                    histRelSyst[i]->SetBinError(binSyst, 0.);
                    histRelSystUncorr[i]->SetBinContent(binSyst, 0);
                    histRelSystUncorr[i]->SetBinError(binSyst, 0.);
                    histRelSystCorr[i]->SetBinContent(binSyst, 0.);
                    histRelSystCorr[i]->SetBinError(binSyst, 0.);
                }    
            }
            for (Int_t j = 0; j< graphComb[i]->GetN(); j++){
                Int_t binComb = histRelComb[i]-> FindBin(graphComb[i]->GetX()[j]);
                if (graphComb[i]->GetY()[j] != 0){
                    histRelComb[i]->SetBinContent(binComb, (graphComb[i]->GetEYhigh()[j]+graphComb[i]->GetEYlow()[j])/(2*graphComb[i]->GetY()[j]));
                    histRelComb[i]->SetBinError(binComb, 0.);
                    histRelCombUncorr[i]->SetBinContent(binComb, (graphComb[i]->GetEYhigh()[j]+graphComb[i]->GetEYlow()[j])/(2*graphComb[i]->GetY()[j]));
                    histRelCombUncorr[i]->SetBinError(binComb, 0.);
                } else {
                    histRelComb[i]->SetBinContent(binComb, 0);
                    histRelComb[i]->SetBinError(binComb, 0.);
                    histRelCombUncorr[i]->SetBinContent(binComb, 0);
                    histRelCombUncorr[i]->SetBinError(binComb, 0.);
                }    
            }
            
            
            if (haveDetailedSys[i] && enableSysExcl[i]){
                for (Int_t iPt = 0; iPt < nPtBinsSysAvailSingle[i]; iPt++){
                    Double_t ptCurr     = ptSysSplit[i][iPt].at(0);
                    Int_t binPt         = histRelSyst[i]->FindBin(ptCurr);
                    cout << histRelStat[i]->FindBin(ptCurr) << "\t"<< histRelSyst[i]->FindBin(ptCurr) << "\t" << histRelComb[i]->FindBin(ptCurr) << endl;
                    Double_t stat       = histRelStat[i]->GetBinContent(binPt);
                    Double_t newSys     = histRelSystUncorr[i]->GetBinContent(binPt)*ptSysSplit[i][iPt].at(2)/ptSysSplit[i][iPt].at(3);
                    Double_t newComb    = TMath::Sqrt(newSys*newSys + stat*stat);
                    histRelSystUncorr[i]->SetBinContent(binPt, newSys);
                    histRelCombUncorr[i]->SetBinContent(binPt, newComb);
                    cout << histRelStat[i]->GetBinContent(binPt) << "\t" << histRelSyst[i]->GetBinContent(binPt) << "\t" <<histRelComb[i]->GetBinContent(binPt) << "\t" <<  histRelSystUncorr[i]->GetBinContent(binPt) << "\t" << histRelCombUncorr[i]->GetBinContent(binPt) << endl;
                    histRelSystCorr[i]->SetBinContent(binPt, ptSysSplit[i][iPt].at(1)/100.);
                    histRelSystCorr[i]->SetBinError(binPt, 0);                    
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
            graphCombReb[i]             = new TGraphAsymmErrors(maxNBins-startBin);
            graphSystReb[i]             = new TGraphAsymmErrors(maxNBins-startBin);
            graphStatReb[i]             = new TGraphAsymmErrors(maxNBins-startBin);
            for(Int_t j = startBin; j<maxNBins; j++){
                Double_t ptCurrent      = (finalBinningPt[j+1]+finalBinningPt[j])/2.;
                Double_t relStatErr     =  histRelStat[i]->Interpolate(ptCurrent);
                Double_t relCombErr     =  histRelCombUncorr[i]->Interpolate(ptCurrent);
                Double_t relSystErr     =  histRelSystUncorr[i]->Interpolate(ptCurrent);
                
                graphCombReb[i]->SetPoint(j-startBin,ptCurrent,fitComb[i]->Eval(ptCurrent));
                if (relCombErr != 0){
                    graphCombReb[i]->SetPointError(j-startBin,(finalBinningPt[j+1]-finalBinningPt[j])/2.,(finalBinningPt[j+1]-finalBinningPt[j])/2.,
                                                   relCombErr*fitComb[i]->Eval(ptCurrent),relCombErr*fitComb[i]->Eval(ptCurrent));
                } else {
                    graphCombReb[i]->SetPointError(j-startBin,(finalBinningPt[j+1]-finalBinningPt[j])/2.,(finalBinningPt[j+1]-finalBinningPt[j])/2.,
                                                   fitComb[i]->Eval(ptCurrent)*0.2,fitComb[i]->Eval(ptCurrent)*0.2);
                }    
                graphSystReb[i]->SetPoint(j-startBin,ptCurrent,fitComb[i]->Eval(ptCurrent));
                if (relSystErr != 0){
                    graphSystReb[i]->SetPointError(j-startBin,(finalBinningPt[j+1]-finalBinningPt[j])/2.,(finalBinningPt[j+1]-finalBinningPt[j])/2.,
                                                   relSystErr*fitComb[i]->Eval(ptCurrent),relSystErr*fitComb[i]->Eval(ptCurrent));
                } else {
                    graphSystReb[i]->SetPointError(j-startBin,(finalBinningPt[j+1]-finalBinningPt[j])/2.,(finalBinningPt[j+1]-finalBinningPt[j])/2.,
                                                   fitComb[i]->Eval(ptCurrent)*0.2,fitComb[i]->Eval(ptCurrent)*0.2);
                }
                graphStatReb[i]->SetPoint(j-startBin,ptCurrent,fitComb[i]->Eval(ptCurrent));
                if (relStatErr != 0){
                    graphStatReb[i]->SetPointError(j-startBin,(finalBinningPt[j+1]-finalBinningPt[j])/2.,(finalBinningPt[j+1]-finalBinningPt[j])/2.,
                                                   relStatErr*fitComb[i]->Eval(ptCurrent),relStatErr*fitComb[i]->Eval(ptCurrent));
                } else {
                    graphStatReb[i]->SetPointError(j-startBin,(finalBinningPt[j+1]-finalBinningPt[j])/2.,(finalBinningPt[j+1]-finalBinningPt[j])/2.,
                                                   fitComb[i]->Eval(ptCurrent)*0.2,fitComb[i]->Eval(ptCurrent)*0.2);
                }
            }
            graphCombReb[i]->Print();
        } else {    
            if (i==0 ){
                graphCombReb[i]             = (TGraphAsymmErrors*)graphComb[i]->Clone(Form("graphCombReb_%d",i));
                graphSystReb[i]             = (TGraphAsymmErrors*)graphSyst[i]->Clone(Form("graphSystReb_%d",i));
                graphStatReb[i]             = (TGraphAsymmErrors*)graphStat[i]->Clone(Form("graphStatReb_%d",i));
                for(Int_t j = 0; j<graphCombReb[i]->GetN(); j++){
                    // take out correlated uncertainty for interpolation
                    Double_t relCombErr     =  histRelCombUncorr[i]->GetBinContent(histRelCombUncorr[i]->FindBin(graphCombReb[i]->GetX()[j]));
                    Double_t relSystErr     =  histRelSystUncorr[i]->GetBinContent(histRelSystUncorr[i]->FindBin(graphCombReb[i]->GetX()[j]));
                    if (relCombErr != 0){
                        graphCombReb[i]->SetPointError(j,graphCombReb[i]->GetEXlow()[j],graphCombReb[i]->GetEXhigh()[j],
                                                       relCombErr*graphCombReb[i]->GetY()[j],relCombErr*graphCombReb[i]->GetY()[j]);
                    } else {
                        graphCombReb[i]->SetPointError(j,graphCombReb[i]->GetEXlow()[j],graphCombReb[i]->GetEXhigh()[j],
                                                       graphCombReb[i]->GetY()[j]*0.2,graphCombReb[i]->GetY()[j]*0.2);
                    }    
                    if (relSystErr != 0){
                        graphSystReb[i]->SetPointError(j,graphSystReb[i]->GetEXlow()[j],graphSystReb[i]->GetEXhigh()[j],
                                                       relSystErr*graphSystReb[i]->GetY()[j],relSystErr*graphSystReb[i]->GetY()[j]);
                    } else {
                        graphSystReb[i]->SetPointError(j,graphSystReb[i]->GetEXlow()[j],graphSystReb[i]->GetEXhigh()[j],
                                                       graphSystReb[i]->GetY()[j]*0.2,graphSystReb[i]->GetY()[j]*0.2);
                    }    
                }
            } else {
                graphCombReb[i]             = (TGraphAsymmErrors*)graphComb[0]->Clone(Form("graphCombReb_%d",i));
                graphSystReb[i]             = (TGraphAsymmErrors*)graphSyst[0]->Clone(Form("graphSystReb_%d",i));
                graphStatReb[i]             = (TGraphAsymmErrors*)graphStat[0]->Clone(Form("graphStatReb_%d",i));
                for(Int_t j = 0; j<graphCombReb[i]->GetN(); j++){
                    Double_t relStatErr     =  histRelStat[i]->Interpolate(graphCombReb[i]->GetX()[j]);
                    // take out correlated uncertainty for interpolation
                    Double_t relCombErr     =  histRelCombUncorr[i]->Interpolate(graphCombReb[i]->GetX()[j]);
                    Double_t relSystErr     =  histRelSystUncorr[i]->Interpolate(graphCombReb[i]->GetX()[j]);
                    graphCombReb[i]->SetPoint(j,graphCombReb[i]->GetX()[j],fitComb[i]->Eval(graphCombReb[i]->GetX()[j]));
                    if (relCombErr != 0){
                        graphCombReb[i]->SetPointError(j,graphCombReb[i]->GetEXlow()[j],graphCombReb[i]->GetEXhigh()[j],
                                                       relCombErr*fitComb[i]->Eval(graphCombReb[i]->GetX()[j]),relCombErr*fitComb[i]->Eval(graphCombReb[i]->GetX()[j]));
                    } else {
                        graphCombReb[i]->SetPointError(j,graphCombReb[i]->GetEXlow()[j],graphCombReb[i]->GetEXhigh()[j],
                                                       fitComb[i]->Eval(graphCombReb[i]->GetX()[j])*0.2,fitComb[i]->Eval(graphCombReb[i]->GetX()[j])*0.2);
                    }    
                    graphSystReb[i]->SetPoint(j,graphSystReb[i]->GetX()[j],fitComb[i]->Eval(graphSystReb[i]->GetX()[j]));
                    if (relSystErr != 0){
                        graphSystReb[i]->SetPointError(j,graphSystReb[i]->GetEXlow()[j],graphSystReb[i]->GetEXhigh()[j],
                                                       relSystErr*fitComb[i]->Eval(graphCombReb[i]->GetX()[j]),relSystErr*fitComb[i]->Eval(graphCombReb[i]->GetX()[j]));
                    } else {
                        graphSystReb[i]->SetPointError(j,graphSystReb[i]->GetEXlow()[j],graphSystReb[i]->GetEXhigh()[j],
                                                       fitComb[i]->Eval(graphSystReb[i]->GetX()[j])*0.2,fitComb[i]->Eval(graphSystReb[i]->GetX()[j])*0.2);
                    }    
                    graphStatReb[i]->SetPoint(j,graphStatReb[i]->GetX()[j],fitComb[i]->Eval(graphStatReb[i]->GetX()[j]));
                    if (relStatErr != 0){
                        graphStatReb[i]->SetPointError(j,graphStatReb[i]->GetEXlow()[j],graphStatReb[i]->GetEXhigh()[j],
                                                       relStatErr*fitComb[i]->Eval(graphStatReb[i]->GetX()[j]),relStatErr*fitComb[i]->Eval(graphStatReb[i]->GetX()[j]));
                    } else {
                        graphStatReb[i]->SetPointError(j,graphStatReb[i]->GetEXlow()[j],graphStatReb[i]->GetEXhigh()[j],
                                                       fitComb[i]->Eval(graphStatReb[i]->GetX()[j])*0.2,fitComb[i]->Eval(graphStatReb[i]->GetX()[j])*0.2);
                    }    

                }
            }
        }
    }

    Int_t nDataSetsMC = nDataSets;
    if (haveMC[0]){
        if (nameFileMCRef.CompareTo("") != 0){
            inputFileMC[nDataSets]              = new TFile(nameFileMCRef);
            if (!inputFileMC[nDataSets]->IsZombie() && nameHistMCRef.CompareTo("") != 0){
                histMC[nDataSets]               = (TH1D*)inputFileMC[nDataSets]->Get(nameHistMCRef.Data());
                if (histMC[nDataSets]){
                    nDataSetsMC++;
                    haveMC[nDataSets]           = kTRUE;
                    graphMC[nDataSets]          = new TGraphAsymmErrors(histMC[nDataSets]);
                    while (graphMC[nDataSets]->GetX()[graphMC[nDataSets]->GetN()-1] > 30) graphMC[nDataSets]->RemovePoint(graphMC[nDataSets]->GetN()-1);
                    while (graphMC[nDataSets]->GetX()[0] < 0.1) graphMC[nDataSets]->RemovePoint(0);
                    fitCombMC[nDataSets]        = DoFitWithTsallis(graphMC[nDataSets],Form("fitCombMC_%d",nDataSets),meson.Data(), graphMC[nDataSets]->GetY()[0],7.,0.2);
                    energyIndName[nDataSets]    = finalEnergy;
                    energyInd[nDataSets]        = energy;
                }
            }
        }
    }    
    
    //*************************************************************************************************
    //*************************** extrapolate spectra stat errors only ********************************
    //*************************************************************************************************
    if (haveMC[0]){
        TGraphAsymmErrors*  graphFinalEnergyMC    = GetInterpolSpectrum2D(      nDataSetsMC,
                                                                                graphMC,
                                                                                energyInd,
                                                                                energy,
                                                                                100,
                                                                                3
                                                                            );

        if( graphPtvsSqrtsMC && gPtvsEnergiesSystemMC && fPowerlawSystemMC ){
            Int_t columns   = 2;
            Int_t rows      = 2;
            Int_t counter   = 0;
            while (columns*rows < graphMC[0]->GetN()){
                if (counter%2 != 0) rows++;
                else columns++;
                counter++;
            }    
            PlotInterpolationPtBins(graphPtvsSqrtsMC,gPtvsEnergiesSystemMC,fPowerlawSystemMC,0x0, graphFinalEnergyMC,columns, rows, Form("%s/%s_%s_MC_Pt_vs_Sqrts.%s",outputDirPlots.Data(),meson.Data(),modeName.Data(), suffix.Data()));
            PlotInterpolationSinglePtBin(graphPtvsSqrtsMC[exampleBin+10],gPtvsEnergiesSystemMC[exampleBin+10], fPowerlawSystemMC[exampleBin+10], 0x0, graphFinalEnergyMC, exampleBin+10, Form("%s/%s_%s_MC_Pt_vs_Sqrts_SinglePtBin.%s",outputDirPlots.Data(),meson.Data(),modeName.Data(), suffix.Data()));
            
            splineAlphaMC                       = new TSpline5("alphaMCSpline",graphAlphaMC); 
        } else {
            cout << "ERROR: NULL pointer - returning..." << endl;
            cout << graphAlphaMC << endl;
            cout << graphPtvsSqrtsMC << endl;
            cout << gPtvsEnergiesSystemMC << endl;
            cout << fPowerlawSystemMC << endl;
            return;
        }
        
    }
    //*************************************************************************************************
    //********************* extrapolate spectra combined errors without common*************************
    //*************************************************************************************************
    TGraphAsymmErrors*  graphFinalEnergyComb1    = GetInterpolSpectrum2D(   nDataSets,
                                                                            graphCombReb,
                                                                            energyInd,
                                                                            energy,
                                                                            nTrials,
                                                                            0
                                                                        );

    if( graphPtvsSqrts && gPtvsEnergiesSystem && fPowerlawSystem ){
        Int_t columns   = 2;
        Int_t rows      = 2;
        Int_t counter   = 0;
        while (columns*rows < graphCombReb[0]->GetN()){
            if (counter%2 != 0) rows++;
            else columns++;
            counter++;
        }    
        PlotInterpolationPtBins(graphPtvsSqrts, gPtvsEnergiesSystem, fPowerlawSystem, fPowerlawSystemSystMC, graphFinalEnergyComb1, columns, rows, Form("%s/%s_%s_CombUnCorr_Pt_vs_Sqrts.%s",outputDirPlots.Data(),meson.Data(),modeName.Data(), suffix.Data()));
        PlotInterpolationSinglePtBin(graphPtvsSqrts[exampleBin], gPtvsEnergiesSystem[exampleBin], fPowerlawSystem[exampleBin], fPowerlawSystemSystMC[exampleBin], graphFinalEnergyComb1, exampleBin, Form("%s/%s_%s_CombUnCorr_Pt_vs_Sqrts_SinglePtBin.%s",outputDirPlots.Data(),meson.Data(),modeName.Data(), suffix.Data()));
        
    } else {
        cout << "ERROR: NULL pointer - returning..." << endl;
        cout << graphAlpha << endl;
        cout << graphPtvsSqrts << endl;
        cout << gPtvsEnergiesSystem << endl;
        cout << fPowerlawSystem << endl;
        return;
    }
    graphFinalEnergyComb1->Print();

    //*************************************************************************************************
    //*************************** extrapolate spectra stat errors only ********************************
    //*************************************************************************************************
    TGraphAsymmErrors*  graphFinalEnergyStat    = GetInterpolSpectrum2D(    nDataSets,
                                                                            graphStatReb,
                                                                            energyInd,
                                                                            energy,
                                                                            nTrials,
                                                                            1
                                                                        );

    if( graphPtvsSqrtsStat && gPtvsEnergiesSystemStat && fPowerlawSystemStat ){
        Int_t columns   = 2;
        Int_t rows      = 2;
        Int_t counter   = 0;
        while (columns*rows < graphCombReb[0]->GetN()){
            if (counter%2 != 0) rows++;
            else columns++;
            counter++;
        }    
        PlotInterpolationPtBins(graphPtvsSqrtsStat, gPtvsEnergiesSystemStat, fPowerlawSystemStat, 0x0, graphFinalEnergyStat, columns, rows, Form("%s/%s_%s_Stat_Pt_vs_Sqrts.%s",outputDirPlots.Data(), meson.Data(), modeName.Data(), suffix.Data()));
        PlotInterpolationSinglePtBin(graphPtvsSqrtsStat[exampleBin], gPtvsEnergiesSystemStat[exampleBin], fPowerlawSystemStat[exampleBin], 0x0, graphFinalEnergyStat, exampleBin, Form("%s/%s_%s_Stat_Pt_vs_Sqrts_SinglePtBin.%s", outputDirPlots.Data(), meson.Data(), modeName.Data(), suffix.Data()));
        
    }else{
        cout << "ERROR: NULL pointer - returning..." << endl;
        cout << graphAlphaStat << endl;
        cout << graphPtvsSqrtsStat << endl;
        cout << gPtvsEnergiesSystemStat << endl;
        cout << fPowerlawSystemStat << endl;
        return;
    }
    for (Int_t iPt = 0; iPt < graphFinalEnergyStat->GetN(); iPt++ ){
        Double_t relSysError    = graphFinalEnergyStat->GetEYhigh()[iPt]/graphFinalEnergyStat->GetY()[iPt];
        Double_t newError       = graphFinalEnergyComb1->GetY()[iPt] * relSysError;
        graphFinalEnergyStat->SetPoint(iPt, graphFinalEnergyStat->GetX()[iPt], graphFinalEnergyComb1->GetY()[iPt]);
        graphFinalEnergyStat->SetPointError(iPt, graphFinalEnergyStat->GetEXlow()[iPt], graphFinalEnergyStat->GetEXhigh()[iPt], newError, newError);
    }    
    graphFinalEnergyStat->Print();

    //*************************************************************************************************
    //*************************** extrapolate spectra uncorr syst errors only *************************
    //*************************************************************************************************
    TGraphAsymmErrors*  graphFinalEnergySyst    = GetInterpolSpectrum2D(    nDataSets,
                                                                            graphSystReb,
                                                                            energyInd,
                                                                            energy,
                                                                            nTrials,
                                                                            2
                                                                        );
    TGraphAsymmErrors* graphFinalEnergyIntSyst  = 0x0;
    
    if( graphPtvsSqrtsSyst && gPtvsEnergiesSystemSyst && fPowerlawSystemSyst ){
        Int_t columns   = 2;
        Int_t rows      = 2;
        Int_t counter   = 0;
        while (columns*rows < graphCombReb[0]->GetN()){
            if (counter%2 != 0) rows++;
            else columns++;
            counter++;
        }    
        PlotInterpolationPtBins(graphPtvsSqrtsSyst, gPtvsEnergiesSystemSyst, fPowerlawSystemSyst, 0x0, graphFinalEnergySyst, columns, rows, Form("%s/%s_%s_Syst_Pt_vs_Sqrts.%s", outputDirPlots.Data(), meson.Data(), modeName.Data(), suffix.Data()));
        PlotInterpolationSinglePtBin(graphPtvsSqrtsSyst[exampleBin], gPtvsEnergiesSystemSyst[exampleBin], fPowerlawSystemSyst[exampleBin], 0x0, graphFinalEnergyStat, exampleBin, Form("%s/%s_%s_Syst_Pt_vs_Sqrts_SinglePtBin.%s", outputDirPlots.Data(), meson.Data(), modeName.Data(), suffix.Data()));
        
    }else{
        cout << "ERROR: NULL pointer - returning..." << endl;
        cout << graphAlphaSyst << endl;
        cout << graphPtvsSqrtsSyst << endl;
        cout << gPtvsEnergiesSystemSyst << endl;
        cout << fPowerlawSystemSyst << endl;
        return;
    }
    for (Int_t iPt = 0; iPt < graphFinalEnergySyst->GetN(); iPt++ ){
        Double_t relSysError    = graphFinalEnergySyst->GetEYhigh()[iPt]/graphFinalEnergySyst->GetY()[iPt];
        Double_t newError       = graphFinalEnergyComb1->GetY()[iPt] * relSysError;
        graphFinalEnergySyst->SetPoint(iPt, graphFinalEnergySyst->GetX()[iPt], graphFinalEnergyComb1->GetY()[iPt]);
        graphFinalEnergySyst->SetPointError(iPt, graphFinalEnergySyst->GetEXlow()[iPt], graphFinalEnergySyst->GetEXhigh()[iPt], newError, newError);
    }    
    graphFinalEnergySyst->Print();
    
    cout << "calculating relative interpolation error" << endl;
    if(graphInterpolSyst){
        while (graphDiffMCIntAndReal->GetX()[graphDiffMCIntAndReal->GetN()-1] > graphFinalEnergySyst->GetX()[graphFinalEnergySyst->GetN()-1])
            graphDiffMCIntAndReal->RemovePoint(graphDiffMCIntAndReal->GetN()-1);
        
        graphDiffMCIntAndReal->Print();
        TSpline5* splineRelDiffMC   = new TSpline5("relDiffSpline",graphDiffMCIntAndReal);
        graphFinalEnergyIntSyst     = (TGraphAsymmErrors*)graphFinalEnergySyst->Clone("graphFinalEnergyInterpolationSyst");    
        for (Int_t iPt = 0; iPt < graphFinalEnergyIntSyst->GetN(); iPt++ ){
            Double_t relSysError1   = TMath::Abs(graphInterpolSyst->GetY()[iPt]-graphFinalEnergyComb1->GetY()[iPt])/graphFinalEnergyComb1->GetY()[iPt];
            Double_t relSysError2   = splineRelDiffMC->Eval(graphFinalEnergyIntSyst->GetX()[iPt]);
            cout << graphFinalEnergyIntSyst->GetX()[iPt] << "\t"<< relSysError1 << "\t"<< relSysError2 << endl;
            Double_t relSysError    = TMath::Sqrt(relSysError1*relSysError1+relSysError2*relSysError2);    
            Double_t newError       = graphFinalEnergyComb1->GetY()[iPt] * relSysError;
            graphFinalEnergyIntSyst->SetPoint(iPt, graphFinalEnergyIntSyst->GetX()[iPt], graphFinalEnergyComb1->GetY()[iPt]);
            graphFinalEnergyIntSyst->SetPointError(iPt, graphFinalEnergyIntSyst->GetEXlow()[iPt], graphFinalEnergyIntSyst->GetEXhigh()[iPt], newError, newError);
        }    
        graphFinalEnergyIntSyst->Print();
    }    
        
    TGraphAsymmErrors* graphFinalEnergyCombWOCorr   = CalculateStatPlusSysErrors(graphFinalEnergyStat,graphFinalEnergySyst);
    
    // plot alpha with proper errors
    if (graphAlphaStat && graphAlphaSyst){
        PlotAlphavsPt(graphAlphaStat, graphAlphaSyst, graphAlphaMC, splineAlphaMC, "pp", Form("%s, %s",mesonString.Data(), detProcess.Data()), Form("%s/%s_%s_Alpha_vs_Pt.%s", outputDirPlots.Data(),meson.Data(),modeName.Data(), suffix.Data()));
    }    

    TF1* fitFinal                       = DoFitWithTsallis(graphFinalEnergyCombWOCorr,Form("fitComb_%s",finalEnergy.Data()),meson.Data(), graphFinalEnergyCombWOCorr->GetY()[0],7,0.2); 
//     TF1* fitFinal                       = DoFitWithTCM(graphFinalEnergyCombWOCorr,Form("fitComb_%s",finalEnergy.Data()),meson.Data(), graphFinalEnergyCombWOCorr->GetY()[0],7.,0.2,graphFinalEnergyCombWOCorr->GetY()[0]/10,0.3);

    graphCombReb[nDataSets]             = graphFinalEnergyCombWOCorr;
    fitComb[nDataSets]                  = fitFinal;
    energyIndName[nDataSets]            = finalEnergy;
    energyInd[nDataSets]                = energy;
    if (maxY < graphFinalEnergyCombWOCorr->GetY()[0])
        maxY        = graphFinalEnergyCombWOCorr->GetY()[0];
    if (minY > graphFinalEnergyCombWOCorr->GetY()[graphFinalEnergyCombWOCorr->GetN()-1])
        minY        = graphFinalEnergyCombWOCorr->GetY()[graphFinalEnergyCombWOCorr->GetN()-1];
    
    //*************************************************************************************************
    //*************************** plotting
    //*************************************************************************************************
    TString axisLabel       = "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (c/GeV)^{2}";
    TString labelWriting    = "InvYield";
    if (maxY > 1000){
        axisLabel       = "#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )";
        labelWriting    = "InvXSection";
    }    
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
        PlotWithFit(canvasDummy2, histo2DDummy, histo2DDummy2, graphComb[i], graphCombReb[i], fitComb[i], Form("input_%s%s%s_withFit",meson.Data(), modeName.Data(), energyIndName[i].Data()), outputDirPlots, suffix, energyIndName[i]);    
    }    
    PlotWithFit(canvasDummy2, histo2DDummy, histo2DDummy2, graphFinalEnergyCombWOCorr, 0x0, fitFinal, Form("compare_%s%s%s_withFit", meson.Data(), modeName.Data(), collSystOut.Data()), outputDirPlots, suffix, energyIndName[nDataSets]);
    PlotGraphsOfAllEnergies(canvasDummy2, histo2DDummy, nDataSets+1, graphCombReb, fitComb, energyIndName, Form("output%s%s_withFit",meson.Data(), modeName.Data()),  outputDirPlots, suffix);
    
    // **************************************************************************************************
    // ************************ prepare histos and graphs for writing ***********************************
    // **************************************************************************************************
    TH1D* histoFinalStat                    = NULL;
    if (maxNBins > 0){
        histoFinalStat                      = new TH1D(Form("histStatErr%s_%s_%s",modeName.Data(),meson.Data(),finalEnergy.Data()),Form("histStatErr_%s_%s",meson.Data(),finalEnergy.Data()),maxNBins,finalBinningPt);
        for(Int_t i=0; i<graphFinalEnergyStat->GetN(); i++){
            Int_t ptBin                     = histoFinalStat->FindBin(graphFinalEnergyStat->GetX()[i]);
            histoFinalStat->SetBinContent(ptBin,graphFinalEnergyStat->GetY()[i]);
            histoFinalStat->SetBinError(ptBin,graphFinalEnergyStat->GetEYhigh()[i]);
        }
    }
    graphFinalEnergySyst->Print();
    TGraphAsymmErrors* graphFinalEnergySystFull = (TGraphAsymmErrors*)graphFinalEnergySyst->Clone(Form("graphSystErr%s_%s_%s",modeName.Data(),meson.Data(),finalEnergy.Data()));

    for (Int_t iPt= 0; iPt< graphFinalEnergySystFull->GetN(); iPt++){
        Double_t value                  = graphFinalEnergySyst->GetY()[iPt];
        Double_t oldError               = graphFinalEnergySyst->GetEYhigh()[iPt];
        Double_t currentRelSys          = oldError/value;
        Double_t relCorrErr             = histRelSystCorr[0]->Interpolate(graphFinalEnergySyst->GetX()[iPt]);
        Double_t newRelSys              = TMath::Sqrt(currentRelSys*currentRelSys+relCorrErr*relCorrErr);
        Double_t newError               = newRelSys*graphFinalEnergySyst->GetY()[iPt];
        graphFinalEnergySystFull->SetPointError(iPt, graphFinalEnergySyst->GetEXlow()[iPt], graphFinalEnergySyst->GetEXhigh()[iPt], newError, newError);
    }
    TGraphAsymmErrors* graphFinalEnergyComb     = CalculateStatPlusSysErrors(graphFinalEnergyStat,graphFinalEnergySystFull);
    
    // **************************************************************************************************
    // ************************* Plot relative errors ***************************************************
    // **************************************************************************************************
    Int_t nErrors                           = 3;
    TGraphAsymmErrors* graphStatRel         = CalculateRelErrUpAsymmGraph( graphFinalEnergyStat, Form("relativeStatError_%s_%s_%s", modeName.Data(), meson.Data(), finalEnergy.Data()));
    TGraphAsymmErrors* graphSystUncorrRel   = CalculateRelErrUpAsymmGraph( graphFinalEnergySyst, Form("relativeUncorrSystError_%s_%s_%s", modeName.Data(), meson.Data(), finalEnergy.Data()));
    TGraphAsymmErrors* graphSystFullRel     = CalculateRelErrUpAsymmGraph( graphFinalEnergySystFull, Form("relativeSystFullError_%s_%s_%s", modeName.Data(), meson.Data(), finalEnergy.Data()));
    TGraphAsymmErrors* graphSystInterRel    = 0x0;
    if (graphFinalEnergyIntSyst){    
        graphSystInterRel                   = CalculateRelErrUpAsymmGraph( graphFinalEnergyIntSyst, Form("relativeSystInterError_%s_%s_%s", modeName.Data(), meson.Data(), finalEnergy.Data()));
        nErrors++;
    }
    
    TCanvas* canvasRelErr                   = new TCanvas("canvasRelErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelErr, 0.08, 0.015, 0.02, 0.09);
//     canvasRelErr->SetLogx();
   
    TH2F * histo2DRelErr;
    histo2DRelErr                   = new TH2F("histo2DRelErr","histo2DRelErr", 11000, 0., 
                                               graphFinalEnergyStat->GetX()[graphFinalEnergyStat->GetN()-1]+ graphFinalEnergyStat->GetEXhigh()[graphFinalEnergyStat->GetN()-1], 1000, 0, 40.5);
    SetStyleHistoTH2ForGraphs(histo2DRelErr, "#it{p}_{T} (GeV/#it{c})","Error (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelErr->Draw("copy");

        TLegend* legendRelErr       = GetAndSetLegend2(0.1, 0.94-(0.035*nErrors), 0.3, 0.94, 32);

        DrawGammaSetMarkerTGraph(graphStatRel, 20, 1.4, kRed+2, kRed+2);
        graphStatRel->Draw("p,same,z");
        legendRelErr->AddEntry(graphStatRel,"stat Err","p");
        DrawGammaSetMarkerTGraph(graphSystUncorrRel, 21, 1.4, kBlue+2, kBlue+2);
        graphSystUncorrRel->Draw("p,same,z");
        legendRelErr->AddEntry(graphSystUncorrRel,"uncorr syst Err","p");
        DrawGammaSetMarkerTGraph(graphSystFullRel, 24, 1.4, kGreen+2, kGreen+2);
        graphSystFullRel->Draw("p,same,z");
        legendRelErr->AddEntry(graphSystFullRel,"full syst Err","p");
        
        if (graphSystInterRel){
            DrawGammaSetMarkerTGraph(graphSystInterRel, 25, 1.4, 807, 807);
            graphSystInterRel->Draw("p,same,z,c");                
            legendRelErr->AddEntry(graphSystInterRel,"interpolation syst Err","p");
        }                
        legendRelErr->Draw();

        TLatex *labelRelErrEnergy       = new TLatex(0.95,0.92,collisionSystem.Data());
        SetStyleTLatex( labelRelErrEnergy, 0.035, 4, 1, 42, kTRUE, 31);
        labelRelErrEnergy->Draw();
        TLatex *labelRelErrProcess      = new TLatex(0.95,0.88,Form("%s, %s",mesonString.Data(), detProcess.Data()));
        SetStyleTLatex( labelRelErrProcess, 0.035, 4, 1, 42, kTRUE, 31);
        labelRelErrProcess->Draw();
        
    canvasRelErr->SaveAs(Form("%s/%s%s_RelErr_%s.%s",outputDirPlots.Data(),meson.Data(), modeName.Data(), collSystOut.Data(), suffix.Data()));

    histo2DRelErr->Draw("copy");

        TLegend* legendRelStatErr       = GetAndSetLegend2(0.1, 0.94-(0.035*(nDataSets+1)), 0.3, 0.94, 32);
        for (Int_t j = 0; j< nDataSets; j++){
            histRelStat[j]->Scale(100);
            DrawGammaSetMarker(histRelStat[j], GetDefaultMarkerStyle(energyIndName[j].Data(),"", ""), GetDefaultMarkerSize(energyIndName[j].Data(),"", ""), GetColorDefaultColor( energyIndName[j].Data(),"", ""), GetColorDefaultColor( energyIndName[j].Data(),"", ""));
            histRelStat[j]->Draw("p,same,z");
            legendRelStatErr->AddEntry(histRelStat[j],Form("stat Err, %s", energyIndName[j].Data()),"p");
        }    
            
        DrawGammaSetMarkerTGraph(graphStatRel, GetDefaultMarkerStyle(finalEnergy.Data(),"", ""), GetDefaultMarkerSize(finalEnergy.Data(),"", ""), GetColorDefaultColor( finalEnergy.Data(),"", ""), GetColorDefaultColor( finalEnergy.Data(),"", ""));
        graphStatRel->Draw("p,same,z");
        legendRelStatErr->AddEntry(graphStatRel,Form("stat Err, %s", finalEnergy.Data()),"p");
        legendRelStatErr->Draw();

        labelRelErrEnergy->Draw();
        labelRelErrProcess->Draw();
        
    canvasRelErr->SaveAs(Form("%s/%s%s_RelStatErrDiffEnergies.%s",outputDirPlots.Data(),meson.Data(), modeName.Data(), suffix.Data()));
    histo2DRelErr->Draw("copy");

        TLegend* legendRelSystErr       = GetAndSetLegend2(0.1, 0.94-(0.035*(nDataSets+1)), 0.3, 0.94, 32);
        for (Int_t j = 0; j< nDataSets; j++){
            histRelSystUncorr[j]->Scale(100);
            DrawGammaSetMarker(histRelSystUncorr[j], GetDefaultMarkerStyle(energyIndName[j].Data(),"", ""), GetDefaultMarkerSize(energyIndName[j].Data(),"", ""), GetColorDefaultColor( energyIndName[j].Data(),"", ""), GetColorDefaultColor( energyIndName[j].Data(),"", ""));
            histRelSystUncorr[j]->Draw("p,same,z");
            legendRelSystErr->AddEntry(histRelSystUncorr[j],Form("uncorr. syst Err, %s", energyIndName[j].Data()),"p");
        }    
            
        DrawGammaSetMarkerTGraph(graphSystUncorrRel, GetDefaultMarkerStyle(finalEnergy.Data(),"", ""), GetDefaultMarkerSize(finalEnergy.Data(),"", ""), GetColorDefaultColor( finalEnergy.Data(),"", ""), GetColorDefaultColor( finalEnergy.Data(),"", ""));
        graphSystUncorrRel->Draw("p,same,z");
        legendRelSystErr->AddEntry(graphSystUncorrRel,Form("uncorr. syst Err, %s", finalEnergy.Data()),"p");
        legendRelSystErr->Draw();

        labelRelErrEnergy->Draw();
        labelRelErrProcess->Draw();
        
    canvasRelErr->SaveAs(Form("%s/%s%s_RelUncorrSystErrDiffEnergies.%s",outputDirPlots.Data(),meson.Data(), modeName.Data(), suffix.Data()));
    
    //*************************************************************************************************
    //*************************** write systematics dat-file ***************************************************
    //*************************************************************************************************
    
    const char *SysErrDatnameMean = Form("%s/SystematicErrorAveragedSinglePP_%s%s_%s.dat",outputDir.Data(), meson.Data(),modeName.Data(), ((TString)ReturnCollisionEnergyOutputString(finalEnergy.Data())).Data());
    fstream SysErrDatAver;
    cout << SysErrDatnameMean << endl;
    SysErrDatAver.open(SysErrDatnameMean, ios::out);
    
    SysErrDatAver << "Pt" << "\t" << "UncorrPP" << "\t";
    if (graphSystInterRel) SysErrDatAver << "Interpolation"  << "\t";
    SysErrDatAver << "TotalErrorUncorrPP" << endl;
    
    for (Int_t l=0;l< graphSystUncorrRel->GetN();l++){
        Double_t totErrWOCorr   = graphSystUncorrRel->GetY()[l];
        SysErrDatAver << graphSystUncorrRel->GetX()[l] << "\t" << graphSystUncorrRel->GetY()[l] << "\t";
        if (graphSystInterRel){
            totErrWOCorr        = TMath::Sqrt(totErrWOCorr*totErrWOCorr+graphSystInterRel->GetY()[l]*graphSystInterRel->GetY()[l]);
            SysErrDatAver << graphSystInterRel->GetY()[l] << "\t";
        }
        SysErrDatAver << totErrWOCorr << endl;
    }
    SysErrDatAver.close();

    
    //*************************************************************************************************
    //*************************** write output file ***************************************************
    //*************************************************************************************************
    TFile *fOutput = new TFile(Form("%s/Interpolation.root",outputDir.Data()),"UPDATE");

        graphFinalEnergyStat->Write(Form("graph%sStatErr%s_%s_%s", labelWriting.Data(), modeName.Data(), meson.Data(), finalEnergy.Data()), TObject::kOverwrite);
        graphFinalEnergySyst->Write(Form("graph%sUnCorrSystErr%s_%s_%s", labelWriting.Data(), modeName.Data(), meson.Data(), finalEnergy.Data()), TObject::kOverwrite);
        graphFinalEnergySystFull->Write(Form("graph%sSystErr%s_%s_%s", labelWriting.Data(), modeName.Data(), meson.Data(), finalEnergy.Data()), TObject::kOverwrite);
        if (graphFinalEnergyIntSyst) graphFinalEnergyIntSyst->Write(Form("graph%sInterpolSystErr%s_%s_%s", labelWriting.Data(), modeName.Data(), meson.Data(), finalEnergy.Data()), TObject::kOverwrite);
        graphFinalEnergyCombWOCorr->Write(Form("graph%sUnCorrCombErr%s_%s_%s", labelWriting.Data(), modeName.Data(), meson.Data(), finalEnergy.Data()), TObject::kOverwrite);
        graphFinalEnergyComb->Write(Form("graph%sCombErr%s_%s_%s", labelWriting.Data(), modeName.Data(), meson.Data(), finalEnergy.Data()), TObject::kOverwrite);
        if (histoFinalStat) histoFinalStat->Write(Form("hist%sStatErr%s_%s_%s", labelWriting.Data(), modeName.Data(), meson.Data(), finalEnergy.Data()), TObject::kOverwrite);

        if(graphAlphaStat) graphAlphaStat->Write(Form("graphAlphaStat%s_%s_%s", modeName.Data(), meson.Data(), finalEnergy.Data()), TObject::kOverwrite);
        if(graphAlphaSyst) graphAlphaSyst->Write(Form("graphAlphaSyst%s_%s_%s", modeName.Data(), meson.Data(), finalEnergy.Data()), TObject::kOverwrite);
        if(graphAlphaMC) graphAlphaMC->Write(Form("graphAlphaMC%s_%s_%s", modeName.Data(), meson.Data(), finalEnergy.Data()), TObject::kOverwrite);
        
        graphStatRel->Write(Form("graphRelStatErr%s_%s_%s", modeName.Data(), meson.Data(), finalEnergy.Data()),TObject::kOverwrite);
        graphSystUncorrRel->Write(Form("graphRelUnCorrSystErr%s_%s_%s", modeName.Data(), meson.Data(), finalEnergy.Data()), TObject::kOverwrite);
        graphSystFullRel->Write(Form("graphRelSystErr%s_%s_%s", modeName.Data(), meson.Data(), finalEnergy.Data()), TObject::kOverwrite);
        if (graphSystInterRel) graphSystInterRel->Write(Form("graphRelInterpolSystErr%s_%s_%s", modeName.Data(), meson.Data(), finalEnergy.Data()), TObject::kOverwrite);
        
    fOutput->Write();
    fOutput->Close();

    return;
}


//________________________________________________________________________________________________________________________
TGraphAsymmErrors *GetInterpolSpectrum2D(Int_t nDataPoints, TGraphAsymmErrors** graphs, Double_t* energies, Double_t dSqrts, Int_t nIterations, Int_t modus )
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
    TGraphAsymmErrors *gInterpol    = new TGraphAsymmErrors(graphs[0]->GetN());
    TGraphAsymmErrors *gInterpolSys = new TGraphAsymmErrors(graphs[0]->GetN());
    TGraphErrors  *gAlpha           = new TGraphErrors(graphs[0]->GetN());
    TGraphErrors** gPtvsSqrts       = new TGraphErrors*[graphs[0]->GetN()];
    TGraphErrors** gPtvsEnergies    = new TGraphErrors*[graphs[0]->GetN()];
    TGraphErrors* gDiff             = new TGraphErrors(graphs[0]->GetN());
    TF1** fPowerlawFits             = new TF1*[graphs[0]->GetN()];
    TF1** fPowerlawFitsMC           = new TF1*[graphs[0]->GetN()];

    for(Int_t i = 0; i < graphs[0]->GetN(); i++){        
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
        
        // run pseudo experiments
        TProfile* pseudoExp[nDataPoints];   
        TRandom3* random[nDataPoints];
        for (Int_t j = 0; j < nDataPoints; j++){
            random[j]               = new TRandom3();
            pseudoExp[j]            = new TProfile(Form("resultsInt_%d",j), Form("resultsInt_%d",j), 1, 0, 1, "S");
        }
        TProfile* pseudoExpFinal    = new TProfile("resultsFinal", "resultsFinal", 1, 0, 1, "S");
        TProfile* pseudoExpAlpha    = new TProfile("resultsAlpha", "resultsAlpha", 1, 0, 1, "S");
        
        for (Int_t k = 0; k< nIterations; k++){
            TGraphErrors *gToFitCurrent = new TGraphErrors(nDataPoints);
            for (Int_t j = 0; j < nDataPoints; j++){
                Double_t currentValue       = random[j]->Gaus(graphs[j]->GetY()[i],graphs[j]->GetEYhigh()[i]);
                pseudoExp[j]->Fill(0.5,currentValue);
                gToFitCurrent->SetPoint(j, energies[j], currentValue);
                gToFitCurrent->SetPointError(j, 0, currentValue*0.001);
            }    
            
            TF1 *fPowerlawCurrent = new TF1("fPowerlawCurrent","[0]*x^([1])", 0,20000);
            if(i==0){
                fPowerlawCurrent->SetParameter(0, 0.005);
                fPowerlawCurrent->SetParameter(1, 0.13);
            }else{
                fPowerlawCurrent->SetParameter(0, 0.1);
                fPowerlawCurrent->SetParameter(1, 2.0);
            }
            for(Int_t l = 0; l < 10; l++) gToFitCurrent->Fit(fPowerlawCurrent,"Q");
            pseudoExpFinal->Fill(0.5,fPowerlawCurrent->Eval(dSqrts));
            pseudoExpAlpha->Fill(0.5,fPowerlawCurrent->GetParameter(1));
            delete fPowerlawCurrent;
        }
        
        // set data points for interpolation
        TGraphErrors *gToFit = new TGraphErrors(nDataPoints);
        for (Int_t j = 0; j < nDataPoints; j++){
            gToFit->SetPoint(j, energies[j], graphs[j]->GetY()[i]);
            gToFit->SetPointError(j, 0, graphs[j]->GetEYhigh()[i]);
        }

        if (splineAlphaMC && modus == 0){
            Double_t alphaMC    = splineAlphaMC->Eval(graphs[0]->GetX()[i]);
            cout << "Alpha MC: " << alphaMC << endl;
            TF1* fPowerlawMC    = new TF1("fPowerlawMC","[0]*x^([1])", 0,20000);
            fPowerlawMC->FixParameter(1,alphaMC);
            for(Int_t l = 0; l < 10; l++) gToFit->Fit(fPowerlawMC,"Q");
            
            gInterpolSys->SetPoint(i, graphs[0]->GetX()[i],fPowerlawMC->Eval(dSqrts));
            gInterpolSys->SetPointError(i, graphs[0]->GetEXlow()[i], graphs[0]->GetEXhigh()[i], 0., 0.);
            
            fPowerlawFitsMC[i]  = fPowerlawMC;
        }

        if (modus != 3){
            gToFit->SetPoint(nDataPoints, dSqrts, pseudoExpFinal->GetBinContent(1));
            gToFit->SetPointError(nDataPoints, 0, pseudoExpFinal->GetBinError(1));
        } else {
            Double_t relDiff    = TMath::Abs(graphs[nDataPoints-1]->GetY()[i] -pseudoExpFinal->GetBinContent(1))/graphs[nDataPoints-1]->GetY()[i];
            if (!(isfinite(relDiff))) relDiff = 0;
            cout << "Diff MC: "  <<relDiff << endl;
            gDiff->SetPoint(i, graphs[nDataPoints-1]->GetX()[i], relDiff);
            gDiff->SetPointError(i, 0, 0);
        }
            
        gInterpol->SetPoint(i, graphs[0]->GetX()[i],pseudoExpFinal->GetBinContent(1));
        gInterpol->SetPointError(i, graphs[0]->GetEXlow()[i], graphs[0]->GetEXhigh()[i], pseudoExpFinal->GetBinError(1), pseudoExpFinal->GetBinError(1));
        
        for (Int_t j = 0; j < nDataPoints; j++){
            cout << energies[j] << "\t" << pseudoExp[j]->GetBinContent(1) << "\t" << pseudoExp[j]->GetBinError(1) << "\t"<< graphs[j]->GetY()[i] << "\t"<< graphs[j]->GetEYhigh()[i]<< endl;
            delete random[j];
            delete pseudoExp[j];
        }
                
        TF1 *fPowerlaw = new TF1("fPowerlaw","[0]*x^([1])", 0,20000);
        fPowerlaw->SetParameter(0, 0.1);
        fPowerlaw->SetParameter(1, 2.0);
        fPowerlaw->SetParLimits(1, 0, 100);

        for(Int_t l = 0; l < 10; l++) gToFit->Fit(fPowerlaw,"Q");
        Double_t alpha      = fPowerlaw->GetParameter(1);
        Double_t alphaE     = TMath::Sqrt(fPowerlaw->GetParError(1)*fPowerlaw->GetParError(1)+pseudoExpAlpha->GetBinError(1)*pseudoExpAlpha->GetBinError(1));
        fPowerlaw->FixParameter(1,pseudoExpAlpha->GetBinContent(1));
        gToFit->Fit(fPowerlaw,"Q");
        
        cout << "pT: ";
        for (Int_t j = 0; j< nDataPoints; j++){
             cout << graphs[j]->GetX()[i] << " % " ;
        }
        cout << ": " << fPowerlaw->GetParameter(0) << " - " << fPowerlaw->GetParameter(1) << endl;
        
        gAlpha->SetPoint(i, graphs[0]->GetX()[i],alpha);
        gAlpha->SetPointError(i, 0, alphaE);

        gPtvsSqrts[i]= new TGraphErrors(1);
        gPtvsSqrts[i]->SetPoint(0,dSqrts,pseudoExpFinal->GetBinContent(1));
        gPtvsSqrts[i]->SetPointError(0,0, pseudoExpFinal->GetBinError(1));

        fPowerlawFits[i] = fPowerlaw;
        gPtvsEnergies[i] = gToFit;

        cout << dSqrts << "\t" << pseudoExpFinal->GetBinContent(1) << "\t" << pseudoExpFinal->GetBinError(1) << endl;
        
        delete pseudoExpFinal;
        
    }

    if (modus == 1){
        graphAlphaStat          = gAlpha;
        graphPtvsSqrtsStat      = new TGraphErrors*[graphs[0]->GetN()];
        fPowerlawSystemStat     = new TF1*[graphs[0]->GetN()];
        gPtvsEnergiesSystemStat = new TGraphErrors*[graphs[0]->GetN()];
        for ( Int_t i = 0; i < graphs[0]->GetN(); i++ ){
            graphPtvsSqrtsStat[i]       = gPtvsSqrts[i];
            fPowerlawSystemStat[i]      = fPowerlawFits[i];
            gPtvsEnergiesSystemStat[i]  = gPtvsEnergies[i];
        }
    } else if (modus == 2){ 
        graphAlphaSyst          = gAlpha;
        graphPtvsSqrtsSyst      = new TGraphErrors*[graphs[0]->GetN()];
        fPowerlawSystemSyst     = new TF1*[graphs[0]->GetN()];
        gPtvsEnergiesSystemSyst = new TGraphErrors*[graphs[0]->GetN()];
        for ( Int_t i = 0; i < graphs[0]->GetN(); i++ ){
            graphPtvsSqrtsSyst[i]       = gPtvsSqrts[i];
            fPowerlawSystemSyst[i]      = fPowerlawFits[i];
            gPtvsEnergiesSystemSyst[i]  = gPtvsEnergies[i];
        }
    } else if (modus == 3){ 
        graphDiffMCIntAndReal   = gDiff;
        graphAlphaMC            = gAlpha;
        graphPtvsSqrtsMC        = new TGraphErrors*[graphs[0]->GetN()];
        fPowerlawSystemMC       = new TF1*[graphs[0]->GetN()];
        gPtvsEnergiesSystemMC   = new TGraphErrors*[graphs[0]->GetN()];
        for ( Int_t i = 0; i < graphs[0]->GetN(); i++ ){
            graphPtvsSqrtsMC[i]         = gPtvsSqrts[i];
            fPowerlawSystemMC[i]        = fPowerlawFits[i];
            gPtvsEnergiesSystemMC[i]    = gPtvsEnergies[i];
        }
    } else { 
        graphAlpha              = gAlpha;
        graphInterpolSyst       = gInterpolSys;
        graphPtvsSqrts          = new TGraphErrors*[graphs[0]->GetN()];
        fPowerlawSystem         = new TF1*[graphs[0]->GetN()];
        gPtvsEnergiesSystem     = new TGraphErrors*[graphs[0]->GetN()];
        if (splineAlphaMC)
            fPowerlawSystemSystMC   = new TF1*[graphs[0]->GetN()];
        for ( Int_t i = 0; i < graphs[0]->GetN(); i++ ){
            graphPtvsSqrts[i]           = gPtvsSqrts[i];
            fPowerlawSystem[i]          = fPowerlawFits[i];
            gPtvsEnergiesSystem[i]      = gPtvsEnergies[i];
            if (splineAlphaMC)
                fPowerlawSystemSystMC[i]    = fPowerlawFitsMC[i];
        }
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
void PlotInterpolationPtBins(TGraphErrors** gPtvSqrts,TGraphErrors** gPtvsEnergies, TF1** fPowerlaw, TF1** fPowerlawMC, TGraphAsymmErrors* gRpPb,Int_t fColumnPlot, Int_t fRowPlot,TString namePlot){

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

        gPtvsEnergies[iPt]->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
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

        if (fPowerlawMC){
            if (fPowerlawMC[iPt]){
                fPowerlawMC[iPt]->SetLineColor(kRed+2);
                fPowerlawMC[iPt]->SetLineWidth(2);
                fPowerlawMC[iPt]->SetLineStyle(7);
                fPowerlawMC[iPt]->Draw("same");
            }    
        }
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
void PlotInterpolationSinglePtBin( TGraphErrors* gPtvSqrts,
                                   TGraphErrors* gPtvsEnergies, 
                                   TF1* fPowerlaw, 
                                   TF1* fPowerlawMC, 
                                   TGraphAsymmErrors* gRpPb, 
                                   Int_t ptBin,
                                   TString namePlot
                                  ){

    TGaxis::SetMaxDigits(3);
    TString nameCanvas = "";
    
    TCanvas * canvasPtvsSqrts     = new TCanvas(nameCanvas.Data(),"",1200,900);  // gives the page size
    canvasPtvsSqrts->SetTopMargin(0.1);
    canvasPtvsSqrts->SetBottomMargin(0.1);
    canvasPtvsSqrts->SetRightMargin(0.055);
    canvasPtvsSqrts->SetLeftMargin(0.1);

        DrawGammaSetMarkerTGraphErr(gPtvSqrts,21,1.5, kRed , kRed);
        DrawGammaSetMarkerTGraphErr(gPtvsEnergies,20,1.5, kBlack , kBlack);
        gPtvsEnergies->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
        gPtvsEnergies->SetTitle("");
        gPtvsEnergies->GetXaxis()->SetTitleSize(0.04);
        gPtvsEnergies->GetXaxis()->SetLabelSize(0.035);
        gPtvsEnergies->GetYaxis()->SetTitle("invariant cross section");
        gPtvsEnergies->GetYaxis()->SetTitleSize(0.04);
        gPtvsEnergies->GetYaxis()->SetLabelSize(0.035);
        gPtvsEnergies->GetXaxis()->SetNdivisions(308,kTRUE);
        gPtvsEnergies->GetYaxis()->SetNdivisions(304,kTRUE);
        gPtvsEnergies->GetXaxis()->SetLabelOffset(0.0);
        gPtvsEnergies->GetYaxis()->SetLabelOffset(0.0);
        gPtvsEnergies->GetXaxis()->SetTitleOffset(1.1);
        gPtvsEnergies->GetYaxis()->SetTitleOffset(1.1);

        gPtvsEnergies->Draw("ap");
        gPtvSqrts->Draw("same, p");
        fPowerlaw->SetLineColor(kBlue);
        fPowerlaw->SetLineWidth(2);
        fPowerlaw->Draw("same");

        if (fPowerlawMC){
            fPowerlawMC->SetLineColor(kRed+2);
            fPowerlawMC->SetLineWidth(2);
            fPowerlawMC->SetLineStyle(7);
            fPowerlawMC->Draw("same");
        }
        TString Title = Form("#it{p}_{T} = %3.2f GeV/#it{c} ",gRpPb->GetX()[ptBin]);


        Double_t yMin = 0.94;
        Double_t xMin = 0.45;

        TLatex *alice = new TLatex(xMin,yMin,Form("%s",Title.Data())); // Bo: this was
        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextSize(0.04);
        alice->Draw();

    canvasPtvsSqrts->Print(namePlot.Data());
    delete canvasPtvsSqrts;
}


//________________________________________________________________________________________________________________________
void PlotAlphavsPt(TGraphErrors* gAlpha, TGraphErrors* gAlphaSyst, TGraphErrors* gAlphaMC, TSpline* splineMC, TString method, TString thesisPlotLabel, TString namePlot){

    TCanvas * canvasAlphavsPt     = new TCanvas("AlphavsPt","",1400,900);  // gives the page size
    DrawGammaCanvasSettings( canvasAlphavsPt,  0.08, 0.02, 0.03, 0.09);

    Int_t nlabels = 2;
    SetStyleTGraphErrorsForGraphs(gAlpha,"#it{p}_{T} (GeV/#it{c})","#alpha", 0.04,0.04, 0.04,0.04, 1.,1., 512, 512);
    DrawGammaSetMarkerTGraphErr(gAlpha,21,1.5, kBlue+2 , kBlue+2);
    gAlpha->Draw("ap");
    DrawGammaSetMarkerTGraphErr(gAlphaSyst, 25, 1.5, kBlue+2 , kBlue+2, 1.4, kTRUE, kBlue-9);        
    gAlphaSyst->Draw("same,3");
    
    if (gAlphaMC){
        nlabels++;
        DrawGammaSetMarkerTGraphErr(gAlphaMC,21,1.5, kRed+2 , kRed+2,1.4, kTRUE, kRed-2);
//         gAlphaMC->SetFillStyle(3144);
        gAlphaMC->Draw("same,3");
    }
    if (splineMC){
        splineMC->SetLineColor(kRed+2);
        splineMC->Draw("same");
    }    
    gAlpha->Draw("same,p");
    
    TLatex *labelThesis = new TLatex(0.13,0.90,thesisPlotLabel.Data());
    SetStyleTLatex( labelThesis, 0.04,4);
    labelThesis->Draw();

    TLegend* legendAlpha     = GetAndSetLegend2(0.12, 0.88-(0.035*nlabels), 0.35, 0.88, 32, 1, "", 43, 0.2);
    legendAlpha->AddEntry(gAlpha,"stat only","p");
    legendAlpha->AddEntry(gAlphaSyst,"uncorr syst only","f");
    if (gAlphaMC)legendAlpha->AddEntry(gAlphaMC,"MC based","f");
    legendAlpha->Draw();
    
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

    for(Int_t i=0; i<nPoints; i++){
        yErrorRHigh[i]          = TMath::Sqrt(TMath::Power(yErrorHigh[i],2)+TMath::Power(yErrorSysHigh[i],2)); 
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

