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
#include "TObjString.h"
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
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"

extern TRandom*    gRandom;
extern TBenchmark*    gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*      gMinuit;

Double_t    xSection900GeV =     47.78*1e-3;
Double_t    xSection2760GeV =    55.416*1e-3;
Double_t    xSection7000GeV =    62.22*1e-3;
Double_t    xSection8000GeV =    55.74*1e-3;
Double_t    recalcBarn =      1e12; //NLO in pbarn!!!!

void readFilePaquet(  TString fileName,
                      TFile &outputfile,
                      TString folderName,
                      TString baseOutputName,
                      TString xAxis,
                      TString yAxis,
                      Bool_t verbose
                   ){

    cout << "INFO: You have chosen the following input file:  " << fileName.Data() << endl;
    ifstream inputDataFile;
    inputDataFile.open(fileName,ios_base::in);
    if (!inputDataFile) {
        cout << "ERROR: file " << fileName.Data() << " not found!" << endl;
        return;
    }

    Int_t nPoints = -1;
    Double_t arrayPT[100];
    Double_t arrayYield[11][100];
    TString centArray[11]    = {"0005", "0510", "1020", "2030", "3040", "4050", "5060", "0010", "0020", "2040", "4060"};
    // read settings from file
    for( TString tempLine; tempLine.ReadLine(inputDataFile, kTRUE); ) {
        // check if line should be considered
        if (tempLine.BeginsWith("%") || tempLine.BeginsWith("#")){
            continue;
        }
        if(verbose) cout << tempLine.Data() << endl;
        nPoints++;

        // Separate the string according to tabulators
        TObjArray *tempArr  = tempLine.Tokenize(" ");
        if(tempArr->GetEntries()<1){
            cout << "nothing to be done" << endl;
            delete tempArr;
            continue;
        }


        // Put them to the correct variables
        arrayPT[nPoints]       = (Double_t)atof(((TString)((TObjString*)tempArr->At(0))->GetString()).Data());
        Int_t inter            = TMath::Nint(arrayPT[nPoints]*10);
        Double_t inter2        = (Double_t)inter/10;
        arrayPT[nPoints]       = inter2;
        arrayYield[0][nPoints] = (Double_t)atof(((TString)((TObjString*)tempArr->At(1))->GetString()).Data());
        arrayYield[1][nPoints] = (Double_t)atof(((TString)((TObjString*)tempArr->At(2))->GetString()).Data());
        arrayYield[2][nPoints] = (Double_t)atof(((TString)((TObjString*)tempArr->At(3))->GetString()).Data());
        arrayYield[3][nPoints] = (Double_t)atof(((TString)((TObjString*)tempArr->At(4))->GetString()).Data());
        arrayYield[4][nPoints] = (Double_t)atof(((TString)((TObjString*)tempArr->At(5))->GetString()).Data());
        arrayYield[5][nPoints] = (Double_t)atof(((TString)((TObjString*)tempArr->At(6))->GetString()).Data());
        arrayYield[6][nPoints] = (Double_t)atof(((TString)((TObjString*)tempArr->At(7))->GetString()).Data());
        // 0-10%
        arrayYield[7][nPoints] = (arrayYield[0][nPoints]+arrayYield[1][nPoints])/2;
        // 0-20%
        arrayYield[8][nPoints] = (arrayYield[7][nPoints]+arrayYield[2][nPoints])/2;
        //20-40%
        arrayYield[9][nPoints] = (arrayYield[3][nPoints]+arrayYield[4][nPoints])/2;
        //40-60%
        arrayYield[10][nPoints] = (arrayYield[5][nPoints]+arrayYield[6][nPoints])/2;

        if(verbose) { cout << arrayPT[nPoints]
                           << " " << arrayYield[0][nPoints]
                           << " " << arrayYield[1][nPoints]
                           << " " << arrayYield[2][nPoints]
                           << " " << arrayYield[3][nPoints]
                           << " " << arrayYield[4][nPoints]
                           << " " << arrayYield[5][nPoints]
                           << " " << arrayYield[6][nPoints]
                           << endl; }
        delete tempArr;
    }

    if(nPoints == -1) {
      cout << "ERROR: did not find any valid data in file" << endl;
      return;
    }

    outputfile.cd(folderName.Data());
    for (Int_t i = 0; i< 11; i++){
        TGraph graphYield(nPoints,arrayPT,arrayYield[i]);
        graphYield.GetYaxis()->SetTitle(yAxis.Data());
        graphYield.GetXaxis()->SetTitle(xAxis.Data());
        graphYield.Write(Form("%s_%s",baseOutputName.Data(),centArray[i].Data()), TObject::kOverwrite);
    }
    return;

}

void readFileBegun( TString fileName,
                    TFile &outputfile,
                    TString folderName1,
                    TString folderName2,
                    TString baseOutputName,
                    TString baseOutputName2,
                    TString xAxis,
                    TString yAxis,
                    TString yAxis2,
                    Bool_t verbose
                  ){

    cout << "INFO: You have chosen the following input file:  " << fileName.Data() << endl;
    ifstream inputDataFile;
    inputDataFile.open(fileName,ios_base::in);
    if (!inputDataFile) {
        cout << "ERROR: file " << fileName.Data() << " not found!" << endl;
        return;
    }

    Int_t nPoints = -1;
    Double_t arrayPT[100];
    Double_t arrayPi0Yield[5][100];
    Double_t arrayEtaYield[5][100];
    Double_t arrayEtaToPi0[5][100];
    TString centArray[5]    = {"0010", "1020", "2040", "4060", "6080"};
    // read settings from file
    for( TString tempLine; tempLine.ReadLine(inputDataFile, kTRUE); ) {
        if (tempLine.BeginsWith("%") || tempLine.BeginsWith("#")){
            continue;
        }
        if(verbose) cout << tempLine.Data() << endl;
        nPoints++;

        // Separate the string according to tabulators
        TObjArray *tempArr  = tempLine.Tokenize(" ");
        if(tempArr->GetEntries()<1){
            cout << "nothing to be done" << endl;
            delete tempArr;
            continue;
        }

        // Put them to the correct variables
        arrayPT[nPoints]          = (Double_t)atof(((TString)((TObjString*)tempArr->At(0))->GetString()).Data());
        arrayEtaYield[0][nPoints] = (Double_t)atof(((TString)((TObjString*)tempArr->At(1))->GetString()).Data());
        arrayEtaYield[1][nPoints] = (Double_t)atof(((TString)((TObjString*)tempArr->At(2))->GetString()).Data());
        arrayEtaYield[2][nPoints] = (Double_t)atof(((TString)((TObjString*)tempArr->At(3))->GetString()).Data());
        arrayEtaYield[3][nPoints] = (Double_t)atof(((TString)((TObjString*)tempArr->At(4))->GetString()).Data());
        arrayEtaYield[4][nPoints] = (Double_t)atof(((TString)((TObjString*)tempArr->At(5))->GetString()).Data());
        arrayPi0Yield[0][nPoints] = (Double_t)atof(((TString)((TObjString*)tempArr->At(6))->GetString()).Data());
        arrayPi0Yield[1][nPoints] = (Double_t)atof(((TString)((TObjString*)tempArr->At(7))->GetString()).Data());
        arrayPi0Yield[2][nPoints] = (Double_t)atof(((TString)((TObjString*)tempArr->At(8))->GetString()).Data());
        arrayPi0Yield[3][nPoints] = (Double_t)atof(((TString)((TObjString*)tempArr->At(9))->GetString()).Data());
        arrayPi0Yield[4][nPoints] = (Double_t)atof(((TString)((TObjString*)tempArr->At(10))->GetString()).Data());

        for (Int_t i = 0; i < 5; i++){
            arrayEtaToPi0[i][nPoints] = arrayEtaYield[i][nPoints]/arrayPi0Yield[i][nPoints];
        }

        if(verbose) { cout << arrayPT[nPoints]
                           << " " << arrayEtaYield[0][nPoints]
                           << " " << arrayEtaYield[1][nPoints]
                           << " " << arrayEtaYield[2][nPoints]
                           << " " << arrayEtaYield[3][nPoints]
                           << " " << arrayEtaYield[4][nPoints]
                           << " " << arrayPi0Yield[0][nPoints]
                           << " " << arrayPi0Yield[1][nPoints]
                           << " " << arrayPi0Yield[2][nPoints]
                           << " " << arrayPi0Yield[3][nPoints]
                           << " " << arrayPi0Yield[4][nPoints]
                           << endl; }
        delete tempArr;
    } // end of loop over lines in input file

    if(nPoints == -1) {
      cout << "ERROR: did not find any valid data in file" << endl;
      return;
    }

    // write graphs to file
    outputfile.cd(folderName1.Data());
    for (Int_t i = 0; i < 5; i++){
        TGraph graphPi0Yield(nPoints,arrayPT,arrayPi0Yield[i]);
        graphPi0Yield.GetYaxis()->SetTitle(yAxis.Data());
        graphPi0Yield.GetXaxis()->SetTitle(xAxis.Data());
        graphPi0Yield.Write(Form("%s_%s",baseOutputName.Data(),centArray[i].Data()), TObject::kOverwrite);

    }
    outputfile.cd(folderName2.Data());
    for (Int_t i = 0; i < 5; i++){
        TGraph graphEtaYield(nPoints,arrayPT,arrayEtaYield[i]);
        graphEtaYield.GetYaxis()->SetTitle(yAxis.Data());
        graphEtaYield.GetXaxis()->SetTitle(xAxis.Data());
        graphEtaYield.Write(Form("%s_%s",baseOutputName.Data(),centArray[i].Data()), TObject::kOverwrite);
        TGraph graphEtaToPi0(nPoints,arrayPT,arrayEtaToPi0[i]);
        graphEtaToPi0.GetYaxis()->SetTitle(yAxis2.Data());
        graphEtaToPi0.GetXaxis()->SetTitle(xAxis.Data());
        graphEtaToPi0.Write(Form("%s_%s",baseOutputName2.Data(),centArray[i].Data()), TObject::kOverwrite);

    }
    return;

}

void readFileBegunPerCentBoth ( TString fileName,
                                TFile &outputfile,
                                TString folderName1,
                                TString folderName2,
                                TString outputName1,
                                TString outputName2,
                                TString xAxis,
                                TString yAxis,
                                TString yAxis2,
                                Bool_t verbose
                            ){

    Double_t ptBegun[100];
    Double_t yieldPi0Begun[100];
    Double_t yieldEtaBegun[100];
    Double_t yieldEtaToPi0Begun[100];

    ifstream  fileBegun;
    fileBegun.open(fileName,ios_base::in);
    cout << fileName << endl;

    Int_t nlin = 0;
    while(!fileBegun.eof() && nlin < 100){
        // pt eta yield pi0 yield eta/pi0 ratio
        fileBegun >> ptBegun[nlin] >> yieldEtaBegun[nlin] >> yieldPi0Begun[nlin] >> yieldEtaToPi0Begun[nlin];
        if (verbose) cout << nlin << "\t "  << ptBegun[nlin] << "\t "  << yieldEtaBegun[nlin] << "\t "  << yieldPi0Begun[nlin] << "\t "  << yieldEtaToPi0Begun[nlin] << endl;
        nlin++;
    }
    fileBegun.close();

    outputfile.cd(folderName1.Data());
        TGraph TheoryBegunPi0 = TGraph(nlin-1, ptBegun, yieldPi0Begun);
        TheoryBegunPi0.GetYaxis()->SetTitle(yAxis.Data());
        TheoryBegunPi0.GetXaxis()->SetTitle(xAxis.Data());
        TheoryBegunPi0.Write(outputName1.Data(), TObject::kOverwrite);
        outputfile.cd(folderName2.Data());
        TGraph TheoryBegunEta = TGraph(nlin-1, ptBegun, yieldEtaBegun);
        TheoryBegunEta.GetYaxis()->SetTitle(yAxis.Data());
        TheoryBegunEta.GetXaxis()->SetTitle(xAxis.Data());
        TheoryBegunEta.Write(outputName1.Data(), TObject::kOverwrite);
        TGraph TheoryBegunEtaToPi0 = TGraph(nlin-1, ptBegun, yieldEtaToPi0Begun);
        TheoryBegunEtaToPi0.GetYaxis()->SetTitle(yAxis2.Data());
        TheoryBegunEtaToPi0.GetXaxis()->SetTitle(xAxis.Data());
        TheoryBegunEtaToPi0.Write(outputName2.Data(), TObject::kOverwrite);

    return;
}


void readFileBegunPerCentSep (  TString fileName1,
                                TString fileName2,
                                TFile &outputfile,
                                TString folderName1,
                                TString folderName2,
                                TString outputName1,
                                TString outputName2,
                                TString xAxis,
                                TString yAxis,
                                TString yAxis2,
                                Bool_t verbose
){

    Int_t MaxNPtBins = 100; //the theory guys asked us to stop at 3 Gev/c
    Int_t nbins;
    Double_t BinLowerEdge, BinUpperEdge;

    Double_t ptPi0[100];
    Double_t yieldPi0[100];
    Double_t yieldErrPi0[100];
    Double_t Pi0LowErr[100];
    Double_t Pi0HighErr[100];

    ifstream  filePi0;
    filePi0.open(fileName1,ios_base::in);
    cout << fileName1.Data() << endl;

    Int_t nlinesLowPt = 0;
    while(!filePi0.eof() && nlinesLowPt < MaxNPtBins){
        filePi0 >> nbins >> ptPi0[nlinesLowPt] >> BinLowerEdge >> BinUpperEdge >> yieldPi0[nlinesLowPt] >> yieldErrPi0[nlinesLowPt];
        if (verbose)cout << nlinesLowPt << "\t "  << ptPi0[nlinesLowPt] << "\t "  << yieldPi0[nlinesLowPt] << " +- "  << yieldErrPi0[nlinesLowPt] << endl;

        // error of the yields is not taken (is put to zero) because they have trouble with the centrality range and it is not a "real" error (just the cent. bin width)
        Pi0LowErr[nlinesLowPt] = 0;
        Pi0HighErr[nlinesLowPt] = 0;
        nlinesLowPt++;
    }
    filePi0.close();

    Double_t ptEta[100];
    Double_t yieldEta[100];
    Double_t yieldErrEta[100];
    Double_t EtaLowErr[100];
    Double_t EtaHighErr[100];

    Double_t ratioEtaPi0[100];

    ifstream  fileEta;
    fileEta.open(fileName2,ios_base::in);
    cout << fileName2.Data() << endl;

    nlinesLowPt = 0;
    while(!fileEta.eof() && nlinesLowPt < MaxNPtBins){
        fileEta >> nbins >> ptEta[nlinesLowPt] >> BinLowerEdge >> BinUpperEdge >> yieldEta[nlinesLowPt] >> yieldErrEta[nlinesLowPt];
        if (verbose)cout << nlinesLowPt << "\t "  << ptEta[nlinesLowPt] << "\t "  << yieldEta[nlinesLowPt] << " +- "  << yieldErrEta[nlinesLowPt] << endl;

        // error of the yields is not taken (is put to zero) because they have trouble with the centrality range and it is not a "real" error (just the cent. bin width)
        EtaLowErr[nlinesLowPt] = 0;
        EtaHighErr[nlinesLowPt] = 0;

        ratioEtaPi0[nlinesLowPt] = yieldEta[nlinesLowPt]/yieldPi0[nlinesLowPt];
        if (verbose)cout << "Eta/Pi0 ratio: " << endl;
        if (verbose)cout << nlinesLowPt << "\t "  << ptEta[nlinesLowPt] << "\t "  << ratioEtaPi0[nlinesLowPt] << endl;

        nlinesLowPt++;
    }
    fileEta.close();

    outputfile.cd(folderName1.Data());
        TGraphAsymmErrors TheoryPi0(nlinesLowPt-1, ptPi0, yieldPi0, 0, 0, Pi0LowErr, Pi0HighErr);
        TheoryPi0.GetYaxis()->SetTitle(yAxis.Data());
        TheoryPi0.GetXaxis()->SetTitle(xAxis.Data());
        TheoryPi0.Write(outputName1.Data(), TObject::kOverwrite);
    outputfile.cd(folderName2.Data());
        TGraphAsymmErrors TheoryEta(nlinesLowPt-1, ptEta, yieldEta, 0, 0, EtaLowErr, EtaHighErr);
        TheoryEta.GetYaxis()->SetTitle(yAxis.Data());
        TheoryEta.GetXaxis()->SetTitle(xAxis.Data());
        TheoryEta.Write(outputName1.Data(), TObject::kOverwrite);
        TGraphAsymmErrors TheoryEtaToPi0(nlinesLowPt-1, ptEta, ratioEtaPi0, 0, 0, 0, 0);
        TheoryEtaToPi0.GetYaxis()->SetTitle(yAxis2.Data());
        TheoryEtaToPi0.GetXaxis()->SetTitle(xAxis.Data());
        TheoryEtaToPi0.Write(outputName2.Data(), TObject::kOverwrite);
    return;
}

void readFileDjordjevic(    TString fileName,
                            TFile &outputfile,
                            TString folderName,
                            TString outputName,
                            TString xAxis,
                            TString yAxis,
                            Bool_t verbose){

    cout << "INFO: You have chosen the following input file:  " << fileName.Data() << endl;
    ifstream inputDataFile;
    inputDataFile.open(fileName,ios_base::in);
    if (!inputDataFile) {
        cout << "ERROR: file " << fileName.Data() << " not found!" << endl;
        return;
    }

    Int_t nPoints = -1;
    Double_t arrayPT[100], arrayPTErr[100];
    Double_t arrayPi0RAAMin[100], arrayPi0RAAMax[100], arrayPi0RAAAvg[100], arrayPi0RAAErr[100];

    // read settings from file
    for( TString tempLine; tempLine.ReadLine(inputDataFile, kTRUE); ) {
        if (tempLine.BeginsWith("%") || tempLine.BeginsWith("#")){
            continue;
        }
        if(verbose) cout << tempLine.Data() << endl;
        nPoints++;

        // Separate the string according to tabulators
        TObjArray *tempArr  = tempLine.Tokenize("\t");
        if(tempArr->GetEntries()<1){
            cout << "nothing to be done" << endl;
            delete tempArr;
            continue;
        }

        // Put them to the correct variables
        arrayPT[nPoints]        = (Double_t)atof(((TString)((TObjString*)tempArr->At(0))->GetString()).Data());
        arrayPTErr[nPoints]     = 0;
        arrayPi0RAAMin[nPoints] = (Double_t)atof(((TString)((TObjString*)tempArr->At(1))->GetString()).Data());
        arrayPi0RAAMax[nPoints] = (Double_t)atof(((TString)((TObjString*)tempArr->At(2))->GetString()).Data());
        arrayPi0RAAAvg[nPoints] = (arrayPi0RAAMin[nPoints] + arrayPi0RAAMax[nPoints] ) / 2.0;
        arrayPi0RAAErr[nPoints] = (arrayPi0RAAMax[nPoints] - arrayPi0RAAMin[nPoints] ) / 2.0;

       if(verbose) { cout << arrayPT[nPoints]
                          << "\t" << arrayPi0RAAMin[nPoints]
                          << "\t" << arrayPi0RAAMax[nPoints]
                          << "\t" << arrayPi0RAAAvg[nPoints]
                          << "\t" << arrayPi0RAAErr[nPoints]
                          << endl; }
        delete tempArr;

    } // end of loop over lines in input file

    if(nPoints == -1) {
      cout << "ERROR: did not find any valid data in file" << endl;
      return;
    }

    // write graphs to file
    outputfile.cd(folderName.Data());
    TGraphErrors graphPi0RAA(nPoints,arrayPT,arrayPi0RAAAvg,arrayPTErr,arrayPi0RAAErr);
    graphPi0RAA.GetYaxis()->SetTitle(yAxis.Data());
    graphPi0RAA.GetXaxis()->SetTitle(xAxis.Data());
    graphPi0RAA.Write(outputName.Data(), TObject::kOverwrite);

}

void readFileMinAndMaxSeparately(   TString fileNameMin,
                                    TString fileNameMax,
                                    TFile &outputfile,
                                    TString folderName,
                                    TString outputName,
                                    TString xAxis,
                                    TString yAxis,
                                    Bool_t verbose){

    cout << "INFO: You have chosen the following input file:  " << fileNameMin.Data() << endl;
    ifstream inputDataFileMin;
    ifstream inputDataFileMax;
    inputDataFileMin.open(fileNameMin,ios_base::in);
    inputDataFileMax.open(fileNameMax,ios_base::in);
    if (!inputDataFileMin) {
        cout << "ERROR: file " << fileNameMin.Data() << " not found!" << endl;
        return;
    }
    if (!inputDataFileMax) {
        cout << "ERROR: file " << fileNameMax.Data() << " not found!" << endl;
        return;
    }

    Int_t nPoints = -1;
    Int_t lMax = 4000;
    Double_t arrayPT[4000], arrayPTErr[4000];
    Double_t arrayPi0RAAMin[4000], arrayPi0RAAMax[4000], arrayPi0RAAAvg[4000], arrayPi0RAAErr[4000];

    // read settings from file
    for( TString tempLine; tempLine.ReadLine(inputDataFileMin, kTRUE); ) {
        if(nPoints < lMax){
            if (tempLine.BeginsWith("%") || tempLine.BeginsWith("#")){
                continue;
            }
            if(verbose) cout << tempLine.Data() << endl;
            nPoints++;

            // Separate the string according to tabulators
            TObjArray *tempArr  = tempLine.Tokenize(" ");
            if(tempArr->GetEntries()<1){
                cout << "nothing to be done" << endl;
                delete tempArr;
                continue;
            }
            cout << ((TString)((TObjString*)tempArr->At(0))->GetString()).Data() << "\t" << ((TString)((TObjString*)tempArr->At(1))->GetString()).Data() << endl;

            // Put them to the correct variables
            arrayPT[nPoints]        = (Double_t)atof(((TString)((TObjString*)tempArr->At(0))->GetString()).Data());
            arrayPTErr[nPoints]     = 0;
            arrayPi0RAAMin[nPoints] = (Double_t)atof(((TString)((TObjString*)tempArr->At(1))->GetString()).Data());

            if(verbose) { cout << arrayPT[nPoints]
                               << "\t" << arrayPi0RAAMin[nPoints]
                               << endl; }
            delete tempArr;
        }
    } // end of loop over lines in input file

    if(nPoints == -1) {
      cout << "ERROR: did not find any valid data in file" << endl;
      return;
    }

    cout << "INFO: You have chosen the following input file:  " << fileNameMax.Data() << endl;
    nPoints=-1;
    // read settings from file
    for( TString tempLine; tempLine.ReadLine(inputDataFileMax, kTRUE); ) {
        if(nPoints < lMax){
            if (tempLine.BeginsWith("%") || tempLine.BeginsWith("#")){
                continue;
            }
            if(verbose) cout << tempLine.Data() << endl;
            nPoints++;

            // Separate the string according to tabulators
            TObjArray *tempArr  = tempLine.Tokenize(" ");
            if(tempArr->GetEntries()<1){
                cout << "nothing to be done" << endl;
                delete tempArr;
                continue;
            }

            // Put them to the correct variables
            arrayPT[nPoints]        = (Double_t)atof(((TString)((TObjString*)tempArr->At(0))->GetString()).Data());
            arrayPi0RAAMax[nPoints] = (Double_t)atof(((TString)((TObjString*)tempArr->At(1))->GetString()).Data());

            if(verbose) { cout << arrayPT[nPoints]
                               << "\t" << arrayPi0RAAMax[nPoints]
                               << endl; }
            delete tempArr;
        }
    } // end of loop over lines in input file

    if(nPoints == -1) {
      cout << "ERROR: did not find any valid data in file" << endl;
      return;
    }

    for(Int_t n=0; n<nPoints+1; n++){
        arrayPi0RAAAvg[n] = (arrayPi0RAAMin[n] + arrayPi0RAAMax[n] ) / 2.0;
        arrayPi0RAAErr[n] = TMath::Abs((arrayPi0RAAMax[n] - arrayPi0RAAMin[n] )) / 2.0;

        if(verbose) { cout << arrayPT[n]
            << "\t" << arrayPi0RAAAvg[n]
            << "\t" << arrayPi0RAAErr[n]
                           << endl; }
    }
    TGraphErrors graphPi0RAA(nPoints+1,arrayPT,arrayPi0RAAAvg,arrayPTErr,arrayPi0RAAErr);

    // write graphs to file
    outputfile.cd(folderName.Data());
        graphPi0RAA.GetYaxis()->SetTitle(yAxis.Data());
        graphPi0RAA.GetXaxis()->SetTitle(xAxis.Data());
        graphPi0RAA.Write( outputName.Data(), TObject::kOverwrite);
    if (verbose) graphPi0RAA.Print();
    return;
}


void readFileWHDG(  TString fileNameMain,
                    TString fileNameMax,
                    TString fileNameMin,
                    TFile &outputfile,
                    TString folderName,
                    TString outputName,
                    TString xAxis,
                    TString yAxis,
                    Bool_t verbose
                 ){

    Double_t WHDG_pT[100];
    Double_t WHDG_ptErr[100];
    Double_t WHDG_Raa[100];
    Double_t WHDG_high[100];
    Double_t WHDG_highErr[100];
    Double_t WHDG_low[100];
    Double_t WHDG_lowErr[100];
    for (Int_t i = 0; i< 100; i++){
        WHDG_ptErr[i]           = {0.00001};
    }

    ifstream WHDG_fileMain(fileNameMain.Data());
    ifstream WHDG_fileLow(fileNameMin.Data());
    ifstream WHDG_fileHigh(fileNameMax.Data());

    Int_t index = 0;
    cout << "reading " << fileNameMain.Data() << endl;
    if (WHDG_fileMain.is_open()){
        while(!WHDG_fileMain.eof()){
            WHDG_fileMain >> WHDG_pT[index] >> WHDG_Raa[index];
            if(verbose) cout << WHDG_pT[index] << "\t" << WHDG_Raa[index] << endl;
            index++;
        }
        WHDG_fileMain.close();
        index = 0;
    }
    cout << "reading " << fileNameMin.Data() << endl;
    if (WHDG_fileLow.is_open()){
        while(!WHDG_fileLow.eof()){
            WHDG_fileLow >> WHDG_pT[index] >> WHDG_low[index];
            WHDG_lowErr[index] = WHDG_Raa[index] - WHDG_low[index];
            if(verbose) cout << WHDG_pT[index] << "\t" << WHDG_low[index]<< "\t" << WHDG_lowErr[index]<< endl;

            index++;
        }
        WHDG_fileLow.close();
        index = 0;
    }
    cout << "reading " << fileNameMax.Data() << endl;
    if (WHDG_fileHigh.is_open()){
        while(!WHDG_fileHigh.eof()){
            WHDG_fileHigh >> WHDG_pT[index] >> WHDG_high[index];
            WHDG_highErr[index] = WHDG_high[index] - WHDG_Raa[index] ;
            if(verbose) cout << WHDG_pT[index] << "\t" << WHDG_high[index]<< "\t"<< WHDG_highErr[index] << endl;
            index++;
        }
        WHDG_fileHigh.close();
    }

    TGraphAsymmErrors gWHDG_Raa     = TGraphAsymmErrors(index-1, WHDG_pT, WHDG_Raa, WHDG_ptErr, WHDG_ptErr, WHDG_lowErr, WHDG_highErr );
    if (verbose) gWHDG_Raa.Print();
    outputfile.cd(folderName.Data());
        gWHDG_Raa.GetXaxis()->SetTitle(xAxis.Data());
        gWHDG_Raa.GetYaxis()->SetTitle(yAxis.Data());
        gWHDG_Raa.Write(outputName.Data(),TObject::kOverwrite);
    return;
}

void readJetQuenchingFiles (TString fileName1,
                            TString fileName2,
                            TString fileName3,
                            TFile &outputfile,
                            TString folderName,
                            TString outputName1,
                            TString outputName2,
                            TString outputName3,
                            TString outputName,
                            TString xAxis,
                            TString yAxis,
                            Bool_t verbose
                           ){

    // Pi0 Raa
    Double_t ptPi0JetQuenching18[100];
    Double_t ptPi0JetQuenching22[100];
    Double_t ptPi0JetQuenching26[100];

    Double_t pi0Raa18[100];
    Double_t pi0Raa22[100];
    Double_t pi0Raa26[100];

    Double_t pi0RaaLowErr[100];
    Double_t pi0RaaHighErr[100];

    Double_t ptErr[100];

    Int_t nlinesJetQuenching = 0;

    ifstream  file1;
    file1.open(fileName1,ios_base::in);
    cout << fileName1 << endl;

    while(!file1.eof() && nlinesJetQuenching < 100){
        file1 >> ptPi0JetQuenching18[nlinesJetQuenching] >> pi0Raa18[nlinesJetQuenching];
        cout << nlinesJetQuenching << "\t "  << ptPi0JetQuenching18[nlinesJetQuenching] << "\t "  << pi0Raa18[nlinesJetQuenching] << endl;
        nlinesJetQuenching++;
    }
    file1.close();


    nlinesJetQuenching = 0;
    ifstream  file2;
    file2.open(fileName2,ios_base::in);
    cout << fileName2 << endl;

    while(!file2.eof() && nlinesJetQuenching < 100){
        file2 >> ptPi0JetQuenching22[nlinesJetQuenching] >> pi0Raa22[nlinesJetQuenching];
        cout << nlinesJetQuenching << "\t "  << ptPi0JetQuenching22[nlinesJetQuenching] << "\t "  << pi0Raa22[nlinesJetQuenching] << endl;
        nlinesJetQuenching++;
    }
    file2.close();


    nlinesJetQuenching = 0;
    ifstream  file3;
    file3.open(fileName3,ios_base::in);
    cout << fileName3 << endl;

    while(!file3.eof() && nlinesJetQuenching < 100){
        file3 >> ptPi0JetQuenching26[nlinesJetQuenching] >> pi0Raa26[nlinesJetQuenching];
        cout << nlinesJetQuenching << "\t "  << ptPi0JetQuenching26[nlinesJetQuenching] << "\t "  << pi0Raa26[nlinesJetQuenching] << endl;
        nlinesJetQuenching++;
    }
    file3.close();


    //Putting the three together:
    for(Int_t k=0; k<nlinesJetQuenching; k++){
        ptErr[k] = 0.;
        if (pi0Raa18[k] > pi0Raa22[k]){
            pi0RaaHighErr[k] = pi0Raa18[k] - pi0Raa22[k];
            pi0RaaLowErr[k] = pi0Raa22[k] - pi0Raa26[k];
        } else {
            pi0RaaHighErr[k] = pi0Raa26[k] - pi0Raa22[k];
            pi0RaaLowErr[k] = pi0Raa22[k] - pi0Raa18[k];
        }
    }
    TGraph graphPi0JetQuenching18 = TGraph(nlinesJetQuenching,ptPi0JetQuenching18,pi0Raa18);
    TGraph graphPi0JetQuenching22 = TGraph(nlinesJetQuenching,ptPi0JetQuenching22,pi0Raa22);
    TGraph graphPi0JetQuenching26 = TGraph(nlinesJetQuenching,ptPi0JetQuenching26,pi0Raa26);
    TGraphAsymmErrors graphPi0RAAJetQuenching = TGraphAsymmErrors(nlinesJetQuenching, ptPi0JetQuenching22, pi0Raa22, ptErr, ptErr, pi0RaaLowErr, pi0RaaHighErr);

    outputfile.cd(folderName.Data());
    graphPi0JetQuenching18.GetXaxis()->SetTitle(xAxis.Data());
    graphPi0JetQuenching18.GetYaxis()->SetTitle(yAxis.Data());
    graphPi0JetQuenching22.GetXaxis()->SetTitle(xAxis.Data());
    graphPi0JetQuenching22.GetYaxis()->SetTitle(yAxis.Data());
    graphPi0JetQuenching26.GetXaxis()->SetTitle(xAxis.Data());
    graphPi0JetQuenching26.GetYaxis()->SetTitle(yAxis.Data());
    graphPi0RAAJetQuenching.GetXaxis()->SetTitle(xAxis.Data());
    graphPi0RAAJetQuenching.GetYaxis()->SetTitle(yAxis.Data());
    graphPi0JetQuenching18.Write(outputName1.Data(),TObject::kOverwrite);
    graphPi0JetQuenching22.Write(outputName2.Data(),TObject::kOverwrite);
    graphPi0JetQuenching26.Write(outputName3.Data(),TObject::kOverwrite);
    graphPi0RAAJetQuenching.Write(outputName.Data(),TObject::kOverwrite);
    return;

}

void readEPOS2014Calc(  TString folderBase,
                        TFile &outputfile,
                        TString folderName,
                        TString outputName1,
                        TString outputName,
                        TString xAxis,
                        TString yAxis,
                        Bool_t verbose
                     ){


    Int_t index = 0;
    const Int_t nFilesEpos              = 6;
    TString label[nFilesEpos]           = {"Epos_pi0_PbPb_2760GeV_00_to_05", "Epos_pi0_PbPb_2760GeV_05_to_10", "Epos_pi0_PbPb_2760GeV_10_to_20",
                                            "Epos_pi0_PbPb_2760GeV_20_to_40", "Epos_pi0_PbPb_2760GeV_40_to_60", "Epos_pi0_PbPb_2760GeV_60_to_80"};
    TString centArray[nFilesEpos]       = {"0005", "0510", "1020", "2040", "4060", "6080"};

    TGraphErrors* gEpos[nFilesEpos];
    TGraph* gEposWOErr[nFilesEpos];

    for (int iFile=0; iFile<nFilesEpos; iFile++) {

        const Int_t nPoints = 70;
        Double_t pt[nPoints], dndptpy[nPoints], dndptpy_staterr[nPoints];

        TString filename = folderBase + label[iFile] + ".txt";
        cout << filename << endl;

        ifstream fEposTxt(filename);

        Int_t iPoint = 0;

        while(!fEposTxt.eof() && iPoint < nPoints){
            fEposTxt >> pt[iPoint] >> dndptpy[iPoint] >> dndptpy_staterr[iPoint];
            if(verbose) cout << pt[iPoint] << "   " << dndptpy[iPoint] << "   " << dndptpy_staterr[iPoint] << endl;
            dndptpy[iPoint]= dndptpy[iPoint]/(pt[iPoint]*2*TMath::Pi());
            dndptpy_staterr[iPoint] = dndptpy_staterr[iPoint]/(pt[iPoint]*2*TMath::Pi());
            if (fEposTxt.fail()) {
                fEposTxt.clear();
                fEposTxt.ignore(0xFFFFFF,'\n');
                continue;
            }

            iPoint++;
        }

        fEposTxt.close();

        gEpos[iFile] = new TGraphErrors(nPoints, pt, dndptpy, 0, dndptpy_staterr);
        gEposWOErr[iFile] = new TGraph(nPoints, pt, dndptpy);
        Int_t i = nPoints-1 ;
        cout <<  pt[i] << endl;
        while (pt[i] > 12.){
            cout  << pt[i] << endl;
            gEpos[iFile]->RemovePoint(  gEpos[iFile]->GetN()-1);
            gEposWOErr[iFile]->RemovePoint(  gEposWOErr[iFile]->GetN()-1);
            i--;
        }
        if (verbose) gEpos[iFile]->Print();

        outputfile.cd(folderName.Data());
            gEpos[iFile]->Write(Form("%s_%s", outputName.Data(), centArray[iFile].Data()), TObject::kOverwrite);
            gEposWOErr[iFile]->Write(Form("%s_%s", outputName1.Data(), centArray[iFile].Data()), TObject::kOverwrite);
    }
    return;
}

void readFileNemchick(  TString fileNameHydro,
                        TString fileNameQCD,
                        TFile &outputfile,
                        TString folderName,
                        TString outputNameBase,
                        TString outputNameBase2,
                        TString outputNameAdd,
                        Double_t ncoll,
                        Int_t switchPoint,
                        Bool_t verbose
                    ){
    cout << "entered" << endl;

    Double_t pTAdditional[43] = {   6.00,  6.25,  6.50,  6.75,  7.00,
                                    7.25,  7.50,  7.75,  8.00,  8.25,
                                    8.50,  8.75,  9.00,  9.25,  9.50,
                                    9.75, 10.00, 10.25, 10.50, 10.75,
                                    11.00, 11.25, 11.50, 11.75, 12.00,
                                    12.25, 12.50, 12.75, 13.00, 13.25,
                                    13.50, 13.75, 14.00, 14.25, 14.50,
                                    14.75, 15.00, 15.50, 16.00, 17.00,
                                    18.00, 19.00, 20.00};
    TF1 *spectrum2760GeV = new TF1("Levy_Dummy","[0] / ( 2 * TMath::Pi())*([1]-1.)*([1]-2.) / ([1]*[2]*([1]*[2]+0.134977*([1]-2.)))  * pow(1.+(sqrt(x*x+0.134977*0.134977)-0.134977)/([1]*[2]), -[1])");
    spectrum2760GeV->SetParameter(0,1.6450861188);
    spectrum2760GeV->SetParameter(1,7.2981303064);
    spectrum2760GeV->SetParameter(2,0.1402413223);

    TGraph* YieldHydro = new TGraph(switchPoint);
    TGraph* YieldHydroSmoothed = new TGraph(140);
    TGraph* YieldELoss = new TGraph(140);
    TGraph* YieldTotal = new TGraph(140);
    Double_t pTYield[140];
    Double_t yYieldHydro[switchPoint];
    Double_t yELossYield[140];
    Double_t yTotalYield[140];
    ifstream file_Yield (fileNameHydro.Data());
    Int_t index = 0;
    if (file_Yield.is_open()){
        while(!file_Yield.eof()){
            file_Yield >> pTYield[index] >> yYieldHydro[index];
            index++;
        }
        file_Yield.close();
        index = 0;
    }
    for(Int_t i=0; i<switchPoint; i++){
        YieldHydro->SetPoint(i, pTYield[i], yYieldHydro[i]);
    }
    for (Int_t i= 0; i < 43; i++){
        pTYield[i+switchPoint] = pTAdditional[i];
    }

    TGraph* RAA = new TGraph(24);
    Double_t pTRAA[24];
    Double_t yRAA[24];

    ifstream file_RAA (fileNameQCD.Data());
    index = 0;
    if (file_RAA.is_open()){
        while(!file_RAA.eof()){
            file_RAA >> pTRAA[index] >> yRAA[index];
            index++;
        }
        file_RAA.close();
        index = 0;
    }
    for(Int_t i=0; i<24; i++){
        cout << pTRAA[i] << endl;
        RAA->SetPoint(i, pTRAA[i], yRAA[i]);
    }
    RAA->Print();

    TGraph* RAA_finerBinning = new TGraph(switchPoint);
    for(Int_t i=0; i<140; i++){
        Double_t evaluatedRAA = RAA->Eval(pTYield[i], 0, "S");
        //       Double_t evaluatedHydro = YieldHydro->Eval(pTYield[i], 0, "S");
        Double_t evaluatedPP = spectrum2760GeV->Eval(pTYield[i]);
        RAA_finerBinning->SetPoint(i, pTYield[i], evaluatedRAA);

        yELossYield[i] = ncoll*evaluatedPP*evaluatedRAA;
        if (i < switchPoint ){
            yTotalYield[i] = yYieldHydro[i] + ncoll*evaluatedPP*evaluatedRAA;
        } else {
            yTotalYield[i] = ncoll*evaluatedPP*evaluatedRAA;
        }
        //       YieldHydroSmoothed->SetPoint(i, pTYield[i], evaluatedHydro);
        YieldELoss->SetPoint(i, pTYield[i], yELossYield[i]);
        YieldTotal->SetPoint(i, pTYield[i], yTotalYield[i]);
    }

    Int_t j = 0;
    while (pTYield[j] < 3){
        j++;
        YieldELoss->RemovePoint(0);
        YieldTotal->RemovePoint(0);
    }

    outputfile.cd(folderName.Data());
        YieldHydro->Write(Form("%sHydro_%s", outputNameBase.Data(), outputNameAdd.Data()), TObject::kOverwrite);
        YieldELoss->Write(Form("%sEloss_%s", outputNameBase.Data(), outputNameAdd.Data()), TObject::kOverwrite);
        YieldTotal->Write(Form("%sTotal_%s", outputNameBase.Data(), outputNameAdd.Data()), TObject::kOverwrite);

        RAA->Write(Form("%s_%s", outputNameBase2.Data(), outputNameAdd.Data()), TObject::kOverwrite);
        RAA_finerBinning->Write(Form("%s_%s_finnerBinning", outputNameBase2.Data(), outputNameAdd.Data()), TObject::kOverwrite);

    return;
}



void ProduceTheoryGraphsPbPb(TString specifier = "", TString energy = ""){

    if(energy.CompareTo("PbPb_2.76TeV")==0){

        // Begun
        // Cracow
        // Djordjevic
        // Jet quenching
        // Xiao-Fang
        // Vitev
        // Horowitz
        // EPOS
        // Nemchick
        // pQCD pp

        StyleSettingsThesis();
        SetPlotStyle();

        TString outputFileName = Form("ExternalInputPbPb/Theory/TheoryCompilationPbPb.root");
        TFile fileTheoryGraphsPbPb(outputFileName,"UPDATE");
        TDirectory *directory2760GeVPi0    = (TDirectory*)fileTheoryGraphsPbPb.Get("Pi0_PbPb_2.76TeV");
        if (!directory2760GeVPi0){
            fileTheoryGraphsPbPb.mkdir("Pi0_PbPb_2.76TeV");
        }
        TDirectory *directory2760GeVEta    = (TDirectory*)fileTheoryGraphsPbPb.Get("Eta_PbPb_2.76TeV");
        if (!directory2760GeVEta){
            fileTheoryGraphsPbPb.mkdir("Eta_PbPb_2.76TeV");
        }
        TDirectory *directory2760GeVPiCh    = (TDirectory*)fileTheoryGraphsPbPb.Get("PiCh_PbPb_2.76TeV");
        if (!directory2760GeVPiCh){
            fileTheoryGraphsPbPb.mkdir("PiCh_PbPb_2.76TeV");
        }
        TDirectory *directory2760GeVKCh    = (TDirectory*)fileTheoryGraphsPbPb.Get("KCh_PbPb_2.76TeV");
        if (!directory2760GeVKCh){
            fileTheoryGraphsPbPb.mkdir("KCh_PbPb_2.76TeV");
        }
        //***************************************************************************************
        // V. Begun - Equilibrium model
        //> http://inspirehep.net/record/1267669
        //> http://inspirehep.net/record/1298405
        //> http://inspirehep.net/record/779957
        //> http://inspirehep.net/record/1352138

        readFileBegunPerCentBoth ( "ExternalInputPbPb/Theory/CracowModel/Begun_EQModel_0-10.txt", fileTheoryGraphsPbPb, "Pi0_PbPb_2.76TeV", "Eta_PbPb_2.76TeV", "Spectra_SHM_EQ_0010", "EtaToPi0_SHM_EQ_0010",
                                   "#it{p}_{T} (GeV/#it{c})", "d#it{N}/(2#pi d#it{y} #it{p}_{T} d#it{p}_{T})", "#eta/#pi^{0}", kFALSE);
        readFileBegunPerCentBoth ( "ExternalInputPbPb/Theory/CracowModel/Begun_EQModel_20-50.txt", fileTheoryGraphsPbPb, "Pi0_PbPb_2.76TeV", "Eta_PbPb_2.76TeV", "Spectra_SHM_EQ_2050", "EtaToPi0_SHM_EQ_2050",
                                   "#it{p}_{T} (GeV/#it{c})", "d#it{N}/(2#pi d#it{y} #it{p}_{T} d#it{p}_{T})", "#eta/#pi^{0}", kFALSE);

        //********************************************************************************************************************************************
        //***************************************************** Cracow model for LHC11h yields *******************************************************
        // PRC 90, 014906 (2014)
        // columns in the file represent:  nBin p_{T} [GeV/c]  p_{T} [GeV/c]_[MIN]   p_{T} [GeV/c]_[MAX]    dN/(2 #pi p_{T} dp_{T} dy)    dN/(2 #pi p_{T} dp_{T} dy)_[ERROR]

        readFileBegunPerCentSep (  "ExternalInputPbPb/Theory/CracowModel/lowPt-chem-non-equilibrium_15Oct2015/pion0-10.txt",
                                   "ExternalInputPbPb/Theory/CracowModel/lowPt-chem-non-equilibrium_15Oct2015/eta0-10.txt",
                                   fileTheoryGraphsPbPb, "Pi0_PbPb_2.76TeV", "Eta_PbPb_2.76TeV", "Spectra_SHM_NEQ_0010", "EtaToPi0_SHM_NEQ_0010",
                                   "#it{p}_{T} (GeV/#it{c})", "d#it{N}/(2#pi d#it{y} #it{p}_{T} d#it{p}_{T})", "#eta/#pi^{0}", kFALSE);
        readFileBegunPerCentSep (  "ExternalInputPbPb/Theory/CracowModel/lowPt-chem-non-equilibrium_15Oct2015/pion20-50.txt",
                                   "ExternalInputPbPb/Theory/CracowModel/lowPt-chem-non-equilibrium_15Oct2015/eta20-50.txt",
                                   fileTheoryGraphsPbPb, "Pi0_PbPb_2.76TeV", "Eta_PbPb_2.76TeV", "Spectra_SHM_NEQ_0010", "EtaToPi0_SHM_NEQ_0010",
                                   "#it{p}_{T} (GeV/#it{c})", "d#it{N}/(2#pi d#it{y} #it{p}_{T} d#it{p}_{T})", "#eta/#pi^{0}", kFALSE);
        readFileBegunPerCentBoth ( "ExternalInputPbPb/Theory/CracowModel/lowPt-chem-non-equilibrium_15Oct2015/KaonsAndPions0-10.txt", fileTheoryGraphsPbPb, "PiCh_PbPb_2.76TeV", "KCh_PbPb_2.76TeV",
                                   "Spectra_SHM_NEQ_0010", "KToPi_SHM_NEQ_0010",
                                   "#it{p}_{T} (GeV/#it{c})", "d#it{N}/(2#pi d#it{y} #it{p}_{T} d#it{p}_{T})", "#eta/#pi^{0}", kFALSE);

//         TGraphAsymmErrors *TheoryCracowNeutralToChargedPionLowPt_0010 = new TGraphAsymmErrors(MaxNPtBins, ptChargedLowPtNonEq_0010, yieldNeutralToChargedPionLowPtNonEq_0010, 0, 0, 0, 0);


        //********************************************************************************************************************************************
        //*************************************** Raa theory - Djordjevic pred for Pi0 2011  *********************************************************
        // citing M Djordjevic, M. Djordjevic and B. Blagojevic, Phys. Lett. B 737 (2014) 298-302
        // and M Djordjevic and M. Djordjevic, Phys. Lett. B 734 (2014) 286-289
        // prediction for 0-10% for LHC11h data

        readFileDjordjevic("ExternalInputPbPb/Theory/DjordjevicPi02011/PiRaa_0-10.txt", fileTheoryGraphsPbPb,
                           "Pi0_PbPb_2.76TeV", "RAA_Djordjevic_0010", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kTRUE);
        readFileDjordjevic("ExternalInputPbPb/Theory/DjordjevicPi02011/PiRaa_20-50.txt", fileTheoryGraphsPbPb,
                           "Pi0_PbPb_2.76TeV", "RAA_Djordjevic_2050", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kTRUE);

        //*********************************************************************************************************************************
        //***************************************************** Jet quenching  ************************************************************
        // prediction for centrality 0-10% for LHC11h data from arXiv:1506.00838

        readJetQuenchingFiles ("ExternalInputPbPb/Theory/JetQuenching/0.9piraa/0.9piraa1.8.dat",
                               "ExternalInputPbPb/Theory/JetQuenching/0.9piraa/0.9piraa2.2.dat",
                               "ExternalInputPbPb/Theory/JetQuenching/0.9piraa/0.9piraa2.6.dat",
                               fileTheoryGraphsPbPb,
                               "Pi0_PbPb_2.76TeV", "RAA_JetQuenching18_0010", "RAA_JetQuenching22_0010", "RAA_JetQuenching26_0010", "RAA_JetQuenching_0010",
                               "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kTRUE) ;

        readJetQuenchingFiles ("ExternalInputPbPb/Theory/JetQuenching/0.9etaraa/0.9etaraa1.8.dat",
                               "ExternalInputPbPb/Theory/JetQuenching/0.9etaraa/0.9etaraa2.2.dat",
                               "ExternalInputPbPb/Theory/JetQuenching/0.9etaraa/0.9etaraa2.6.dat",
                                fileTheoryGraphsPbPb,
                                "Eta_PbPb_2.76TeV", "RAA_JetQuenching18_0010", "RAA_JetQuenching22_0010", "RAA_JetQuenching26_0010", "RAA_JetQuenching_0010",
                                "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kTRUE) ;
        readJetQuenchingFiles ("ExternalInputPbPb/Theory/JetQuenching/0.9ratio/0.9ratio1.8.dat",
                               "ExternalInputPbPb/Theory/JetQuenching/0.9ratio/0.9ratio2.2.dat",
                               "ExternalInputPbPb/Theory/JetQuenching/0.9ratio/0.9ratio2.6.dat",
                                fileTheoryGraphsPbPb,
                                "Eta_PbPb_2.76TeV", "EtaToPi0_JetQuenching18_0010", "EtaToPi0_JetQuenching22_0010", "EtaToPi0_JetQuenching26_0010", "EtaToPi0_JetQuenching_0010",
                                "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", kTRUE) ;

        //********************************************************************************************************************************************
        //*************************************** Raa theory Xiao-Fang *******************************************************************************

        readFileMinAndMaxSeparately("ExternalInputPbPb/Theory/Xiao/0-20-low.dat", "ExternalInputPbPb/Theory/Xiao/0-20-up.dat", fileTheoryGraphsPbPb,
                                    "Pi0_PbPb_2.76TeV", "RAA_Xiao_0020", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}",   kFALSE);
        readFileMinAndMaxSeparately("ExternalInputPbPb/Theory/Xiao/20-40-low.dat", "ExternalInputPbPb/Theory/Xiao/20-40-up.dat", fileTheoryGraphsPbPb,
                                    "Pi0_PbPb_2.76TeV", "RAA_Xiao_2040", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}" ,  kFALSE);
        readFileMinAndMaxSeparately("ExternalInputPbPb/Theory/Xiao/40-60-low.dat", "ExternalInputPbPb/Theory/Xiao/40-60-up.dat", fileTheoryGraphsPbPb,
                                    "Pi0_PbPb_2.76TeV", "RAA_Xiao_4060", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}" ,  kFALSE);
        readFileMinAndMaxSeparately("ExternalInputPbPb/Theory/Xiao/60-80-low.dat", "ExternalInputPbPb/Theory/Xiao/60-80-up.dat", fileTheoryGraphsPbPb,
                                    "Pi0_PbPb_2.76TeV", "RAA_Xiao_4060", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}" ,  kFALSE);

        // *************************************************************************************
        // new predictions from Vitev: new effective theory (SCET_G) + new technique of evaluating
        // the production of hadrons in the presence of a QGP based on QCD evolution techniques
        // Phys.Rev.Lett. 114 (2015) no.9, 092002
        // Phys.Rev. D93 (2016) no.7, 074030

        readFileMinAndMaxSeparately("ExternalInputPbPb/Theory/Vitev_2760GeV/R-N.aa_010_cron1.5_eloss0.2760GeVpi.g1.9",
                                    "ExternalInputPbPb/Theory/Vitev_2760GeV/R-N.aa_010_cron1.5_eloss0.2760GeVpi.g2.1", fileTheoryGraphsPbPb,
                                    "Pi0_PbPb_2.76TeV", "RAA_Vitev_0010", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}" ,  kFALSE);
//         readFileMinAndMaxSeparately("ExternalInputPbPb/Theory/Vitev_2760GeV/R-N.aa_2050_cron1.5_eloss0.2760GeVpi.g1.9",
//                                     "ExternalInputPbPb/Theory/Vitev_2760GeV/R-N.aa_2050_cron1.5_eloss0.2760GeVpi.g2.1", fileTheoryGraphsPbPb,
//                                     "Pi0_PbPb_2.76TeV", "RAA_Vitev_2050", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}" ,  kFALSE);

        // **************************** Vitev **************************************************
        readFileMinAndMaxSeparately("ExternalInputPbPb/Theory/Vitev_2760GeV/R-PbPb2760pi0.dnbasPI0", "ExternalInputPbPb/Theory/Vitev_2760GeV/R-PbPb2760pi0.upbasPI0",
                                    fileTheoryGraphsPbPb, "Pi0_PbPb_2.76TeV", "RAA_VitevBas_0020", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}" ,  kTRUE);
        readFileMinAndMaxSeparately("ExternalInputPbPb/Theory/Vitev_2760GeV/R-PbPb2760pi0.dnShISelPI0", "ExternalInputPbPb/Theory/Vitev_2760GeV/R-PbPb2760pi0.upShISelPI0",
                                    fileTheoryGraphsPbPb, "Pi0_PbPb_2.76TeV", "RAA_VitevShlSel_0020", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}" ,  kFALSE);
        readFileMinAndMaxSeparately("ExternalInputPbPb/Theory/Vitev_2760GeV/R-PbPb2760pi0.05.dn", "ExternalInputPbPb/Theory/Vitev_2760GeV/R-PbPb2760pi0.05.up",
                                    fileTheoryGraphsPbPb, "Pi0_PbPb_2.76TeV", "RAA_VitevBas_0005", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}" ,  kFALSE);
        readFileMinAndMaxSeparately("ExternalInputPbPb/Theory/Vitev_2760GeV/R-PbPb2760pi0.510.dn", "ExternalInputPbPb/Theory/Vitev_2760GeV/R-PbPb2760pi0.510.up",
                                    fileTheoryGraphsPbPb, "Pi0_PbPb_2.76TeV", "RAA_VitevBas_0510", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}" ,  kFALSE);
        readFileMinAndMaxSeparately("ExternalInputPbPb/Theory/Vitev_2760GeV/R-PbPb2760pi0.1020.dn", "ExternalInputPbPb/Theory/Vitev_2760GeV/R-PbPb2760pi0.1020.up",
                                    fileTheoryGraphsPbPb, "Pi0_PbPb_2.76TeV", "RAA_VitevBas_1020", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}" ,  kFALSE);
        readFileMinAndMaxSeparately("ExternalInputPbPb/Theory/Vitev_2760GeV/R-PbPb2760pi0.2040.dn", "ExternalInputPbPb/Theory/Vitev_2760GeV/R-PbPb2760pi0.2040.up",
                                    fileTheoryGraphsPbPb, "Pi0_PbPb_2.76TeV", "RAA_VitevBas_2040", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}" ,  kFALSE);
        readFileMinAndMaxSeparately("ExternalInputPbPb/Theory/Vitev_2760GeV/R-PbPb2760pi0.4060.dn", "ExternalInputPbPb/Theory/Vitev_2760GeV/R-PbPb2760pi0.4060.up",
                                    fileTheoryGraphsPbPb, "Pi0_PbPb_2.76TeV", "RAA_VitevBas_4060", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}" ,  kFALSE);
        readFileMinAndMaxSeparately("ExternalInputPbPb/Theory/Vitev_2760GeV/R-PbPb2760pi0.6080.dn", "ExternalInputPbPb/Theory/Vitev_2760GeV/R-PbPb2760pi0.6080.up",
                                    fileTheoryGraphsPbPb, "Pi0_PbPb_2.76TeV", "RAA_VitevBas_6080", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}" ,  kFALSE);

        //************************************************************************************************************************
        //*************************************** Raa theory WHDG ****************************************************************
        //from S. Wicks, W. Horowitz, M. Djordjevic and M. Gyulassy, ``Elastic, inelastic, and path length fluctuations in jet tomography,'' Nucl. Phys. A 784, 426 (2007) [nucl-th/0512076].
        //and ``The Surprising Transparency of the sQGP at LHC,'' Nucl. Phys. A 872, 265 (2011) [arXiv:1104.4958 [hep-ph]]

        //*****************************************************     Pi0     ************************************************************************
        // 0-5%, 5-10%, 0-10%, 10-20%, 0-20%, 20-40%, 20-50%, 40-60%, 60-80%
        readFileWHDG ("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA0005b.dat",
                      "ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA0005u.dat",
                      "ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA0005l.dat",
                      fileTheoryGraphsPbPb, "Pi0_PbPb_2.76TeV", "RAA_WHDG_0005", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kFALSE);
        readFileWHDG ("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA0510b.dat",
                      "ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA0510u.dat",
                      "ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA0510l.dat",
                      fileTheoryGraphsPbPb, "Pi0_PbPb_2.76TeV", "RAA_WHDG_0510", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kFALSE);
        readFileWHDG ("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA0010b.dat",
                      "ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA0010u.dat",
                      "ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA0010l.dat",
                      fileTheoryGraphsPbPb, "Pi0_PbPb_2.76TeV", "RAA_WHDG_0010", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kFALSE);
        readFileWHDG ("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA1020b.dat",
                      "ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA1020u.dat",
                      "ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760Pi0RAA1020l.dat",
                      fileTheoryGraphsPbPb, "Pi0_PbPb_2.76TeV", "RAA_WHDG_1020", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kFALSE);
        readFileWHDG ("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA0020b.dat",
                      "ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA0020u.dat",
                      "ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA0020l.dat",
                      fileTheoryGraphsPbPb, "Pi0_PbPb_2.76TeV", "RAA_WHDG_0020", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kFALSE);
        readFileWHDG ("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA2040b.dat",
                      "ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA2040u.dat",
                      "ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA2040l.dat",
                      fileTheoryGraphsPbPb, "Pi0_PbPb_2.76TeV", "RAA_WHDG_2040", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kFALSE);
        readFileWHDG ("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA2050b.dat",
                      "ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA2050u.dat",
                      "ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA2050l.dat",
                      fileTheoryGraphsPbPb, "Pi0_PbPb_2.76TeV", "RAA_WHDG_2050", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kFALSE);
        readFileWHDG ("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA4060b.dat",
                      "ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA4060u.dat",
                      "ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA4060l.dat",
                      fileTheoryGraphsPbPb, "Pi0_PbPb_2.76TeV", "RAA_WHDG_4060", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kFALSE);
        readFileWHDG ("ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA6080b.dat",
                      "ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA6080u.dat",
                      "ExternalInputPbPb/Theory/WHDG2760pi0RAA/WHDG2760pi0RAA6080l.dat",
                      fileTheoryGraphsPbPb, "Pi0_PbPb_2.76TeV", "RAA_WHDG_6080", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kFALSE);

        //*****************************************************     Eta     ************************************************************************
        // 0-10%, (available 10-20%), 20-50%, (available 50-80%)
        readFileWHDG ("ExternalInputPbPb/Theory/WHDG2760etaRAA/WHDGetaRAA0010b.dat",
                      "ExternalInputPbPb/Theory/WHDG2760etaRAA/WHDGetaRAA0010u.dat",
                      "ExternalInputPbPb/Theory/WHDG2760etaRAA/WHDGetaRAA0010l.dat",
                      fileTheoryGraphsPbPb, "Eta_PbPb_2.76TeV", "RAA_WHDG_0010", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kFALSE);
        readFileWHDG ("ExternalInputPbPb/Theory/WHDG2760etaRAA/WHDGetaRAA1020b.dat",
                      "ExternalInputPbPb/Theory/WHDG2760etaRAA/WHDGetaRAA1020l.dat",
                      "ExternalInputPbPb/Theory/WHDG2760etaRAA/WHDGetaRAA1020u.dat",
                      fileTheoryGraphsPbPb, "Eta_PbPb_2.76TeV", "RAA_WHDG_1020", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kFALSE);
        readFileWHDG ("ExternalInputPbPb/Theory/WHDG2760etaRAA/WHDGetaRAA2050b.dat",
                      "ExternalInputPbPb/Theory/WHDG2760etaRAA/WHDGetaRAA2050u.dat",
                      "ExternalInputPbPb/Theory/WHDG2760etaRAA/WHDGetaRAA2050l.dat",
                      fileTheoryGraphsPbPb, "Eta_PbPb_2.76TeV", "RAA_WHDG_2050", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kFALSE);
        readFileWHDG ("ExternalInputPbPb/Theory/WHDG2760etaRAA/WHDGetaRAA5080b.dat",
                      "ExternalInputPbPb/Theory/WHDG2760etaRAA/WHDGetaRAA5080l.dat",
                      "ExternalInputPbPb/Theory/WHDG2760etaRAA/WHDGetaRAA5080u.dat",
                      fileTheoryGraphsPbPb, "Eta_PbPb_2.76TeV", "RAA_WHDG_5080", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kFALSE);

        //*********************************************************************************************************************************
        //********************************************************** Nemchick *************************************************************
        Double_t ncoll0005 = 1684.4;;
        Double_t ncoll2040 = 438.4;
        Double_t ncoll6080 = 26.71;

        readFileNemchick(  "ExternalInputPbPb/Theory/Nemchick/y_AA-pi0-hydro-0-5L.dat",
                           "ExternalInputPbPb/Theory/Nemchick/r_AA-pQCD-0-5L.dat",
                           fileTheoryGraphsPbPb,
                           "Pi0_PbPb_2.76TeV", "Spectra_Kopeliovich", "RAA_Kopeliovich", "0005", ncoll0005, 97, kFALSE);
        readFileNemchick(  "ExternalInputPbPb/Theory/Nemchick/y_AA-pi0-hydro-20-40L.dat",
                           "ExternalInputPbPb/Theory/Nemchick/r_AA-pQCD-20-40L.dat",
                           fileTheoryGraphsPbPb,
                           "Pi0_PbPb_2.76TeV", "Spectra_Kopeliovich", "RAA_Kopeliovich", "2040", ncoll2040, 97, kFALSE);
        readFileNemchick(  "ExternalInputPbPb/Theory/Nemchick/y_AA-pi0-hydro-60-80L.dat",
                           "ExternalInputPbPb/Theory/Nemchick/r_AA-pQCD-60-80L.dat",
                           fileTheoryGraphsPbPb,
                           "Pi0_PbPb_2.76TeV", "Spectra_Kopeliovich", "RAA_Kopeliovich", "6080", ncoll6080, 79, kFALSE);

        //*********************************************************************************************************************************
        //*************************************************************** EPOS ************************************************************
        readEPOS2014Calc(  "ExternalInputPbPb/Theory/EPOS/", fileTheoryGraphsPbPb, "Pi0_PbPb_2.76TeV", "Spectra_EPOS", "SpectrWOErr_EPOS",
                           "#it{p}_{T} (GeV/#it{c})", "d#it{N}/(2#pi d#it{y} #it{p}_{T} d#it{p}_{T})", kFALSE);


        //EPOS(2?) (file from Anders Knospe)
        //necessary to do -> dndptpy[iPoint] = dndptpy[iPoint]/(pt[iPoint]*2*TMath::Pi())
        TFile* filePredictionEpos = new TFile("ExternalInputPbPb/Theory/EPOS/pi0_eta_Anders19Feb2016_NoUrQMD.root");
        TH1D *histoEPOSPi0_0010 = (TH1D*)filePredictionEpos->Get("rb_pi0_pt_cent0_10");
        for (Int_t i = 1; i < histoEPOSPi0_0010->GetNbinsX()+1 ; i++){
            Double_t newBinContent = histoEPOSPi0_0010->GetBinContent(i)/(2*TMath::Pi()*histoEPOSPi0_0010->GetBinCenter(i));
            Double_t newBinError = histoEPOSPi0_0010->GetBinError(i)/(2*TMath::Pi()*histoEPOSPi0_0010->GetBinCenter(i));
            histoEPOSPi0_0010->SetBinContent(i,newBinContent);
            histoEPOSPi0_0010->SetBinError(i,newBinError);
        }
        TH1D *histoEPOSEta_0010 = (TH1D*)filePredictionEpos->Get("rb_eta_pt_cent0_10");
        for (Int_t i = 1; i < histoEPOSEta_0010->GetNbinsX()+1 ; i++){
            Double_t newBinContent = histoEPOSEta_0010->GetBinContent(i)/(2*TMath::Pi()*histoEPOSEta_0010->GetBinCenter(i));
            Double_t newBinError = histoEPOSEta_0010->GetBinError(i)/(2*TMath::Pi()*histoEPOSEta_0010->GetBinCenter(i));
            histoEPOSEta_0010->SetBinContent(i,newBinContent);
            histoEPOSEta_0010->SetBinError(i,newBinError);
        }
        TH1D *histoEPOSPi0_2050 = (TH1D*)filePredictionEpos->Get("rb_pi0_pt_cent20_50");
        for (Int_t i = 1; i < histoEPOSPi0_2050->GetNbinsX()+1 ; i++){
            Double_t newBinContent = histoEPOSPi0_2050->GetBinContent(i)/(2*TMath::Pi()*histoEPOSPi0_2050->GetBinCenter(i));
            Double_t newBinError = histoEPOSPi0_2050->GetBinError(i)/(2*TMath::Pi()*histoEPOSPi0_2050->GetBinCenter(i));
            histoEPOSPi0_2050->SetBinContent(i,newBinContent);
            histoEPOSPi0_2050->SetBinError(i,newBinError);
        }
        TH1D *histoEPOSEta_2050 = (TH1D*)filePredictionEpos->Get("rb_eta_pt_cent20_50");
        for (Int_t i = 1; i < histoEPOSEta_2050->GetNbinsX()+1 ; i++){
            Double_t newBinContent = histoEPOSEta_2050->GetBinContent(i)/(2*TMath::Pi()*histoEPOSEta_2050->GetBinCenter(i));
            Double_t newBinError = histoEPOSEta_2050->GetBinError(i)/(2*TMath::Pi()*histoEPOSEta_2050->GetBinCenter(i));
            histoEPOSEta_2050->SetBinContent(i,newBinContent);
            histoEPOSEta_2050->SetBinError(i,newBinError);
        }

        TH1D* histoEPOSEtaToPi0Ratio_0010 = (TH1D*)histoEPOSEta_0010->Clone("histoEPOSEtaToPi0Ratio_0010");
        histoEPOSEtaToPi0Ratio_0010->Divide(histoEPOSEta_0010,histoEPOSPi0_0010,1.,1.,"");
        TH1D* histoEPOSEtaToPi0Ratio_2050 = (TH1D*)histoEPOSEta_2050->Clone("histoEPOSEtaToPi0Ratio_2050");
        histoEPOSEtaToPi0Ratio_2050->Divide(histoEPOSEta_2050,histoEPOSPi0_2050,1.,1.,"");

        TGraphErrors *graphEPOSPi0_0010 = new TGraphErrors(histoEPOSPi0_0010);
        TGraphErrors *graphEPOSEta_0010 = new TGraphErrors(histoEPOSEta_0010);
        TGraphErrors *graphEPOSPi0_2050 = new TGraphErrors(histoEPOSPi0_2050);
        TGraphErrors *graphEPOSEta_2050 = new TGraphErrors(histoEPOSEta_2050);
        TGraphErrors *graphEPOSEtaToPi0Ratio_0010 = new TGraphErrors(histoEPOSEtaToPi0Ratio_0010);
        TGraphErrors *graphEPOSEtaToPi0Ratio_2050 = new TGraphErrors(histoEPOSEtaToPi0Ratio_2050);

        //*********************************************************************************************************************************
        //************************************************************** EPOS 3 ***********************************************************
        // Calculation for 0-10% (for LHC11h -> is worse than old EPOS...)

        Double_t Pi0pt[70];
        Double_t Etapt[70];
        Double_t Pi0dndptpy[70];
        Double_t Pi0dndptpy_staterr[70];
        Double_t Etadndptpy[70];
        Double_t Etadndptpy_staterr[70];

        TString filenamePi0 = "ExternalInputPbPb/Theory/EPOS/pion0-10.txt";
        ifstream fEposPi0Txt;
        fEposPi0Txt.open(filenamePi0,ios_base::in);
        cout << filenamePi0 << endl;

        Int_t iPoint = 0;
        while(!fEposPi0Txt.eof() && iPoint < 70){
            fEposPi0Txt >> Pi0pt[iPoint] >> Pi0dndptpy[iPoint] >> Pi0dndptpy_staterr[iPoint];
            Pi0dndptpy[iPoint]= Pi0dndptpy[iPoint]/(Pi0pt[iPoint]*2*TMath::Pi());
            Pi0dndptpy_staterr[iPoint] = Pi0dndptpy_staterr[iPoint]/(Pi0pt[iPoint]*2*TMath::Pi());

            iPoint++;
        }
        fEposPi0Txt.close();
        TGraphErrors* graphEPOS3Pi0_0010 = new TGraphErrors(iPoint, Pi0pt, Pi0dndptpy, 0, Pi0dndptpy_staterr);


        TString filenameEta = "ExternalInputPbPb/Theory/EPOS/eta0-10.txt";
        ifstream fEposEtaTxt;
        fEposEtaTxt.open(filenameEta,ios_base::in);
        cout << filenameEta << endl;

        iPoint = 0;
        while(!fEposEtaTxt.eof() && iPoint < 70){
            fEposEtaTxt >> Etapt[iPoint] >> Etadndptpy[iPoint] >> Etadndptpy_staterr[iPoint];
            Etadndptpy[iPoint]= Etadndptpy[iPoint]/(Etapt[iPoint]*2*TMath::Pi());
            Etadndptpy_staterr[iPoint] = Etadndptpy_staterr[iPoint]/(Etapt[iPoint]*2*TMath::Pi());

            iPoint++;
        }
        fEposEtaTxt.close();
        TGraphErrors* graphEPOS3Eta_0010 = new TGraphErrors(iPoint, Etapt, Etadndptpy, 0, Etadndptpy_staterr);


        //**************************************************************************************************************
        //************************************* pQCD, pp 2.76TeV Pi0 DSS14 *********************************************

        Double_t xSection2760GeVppINEL = 62.8*1e9;
        Double_t       ptNLODSS14Pi02760GeV[39];
        Double_t       ptErrNLODSS14Pi02760GeV[39];
        Double_t       muHalfDSS14Pi02760GeV[39];
        Double_t       muHalfDSS14Pi02760GeVforGraph[39];
        Double_t       muHalfErrDSS14Pi02760GeV[39];
        Double_t       muOneDSS14Pi02760GeV[39];
        Double_t       muOneErrDSS14Pi02760GeV[39];
        Double_t       muTwoDSS14Pi02760GeV[39];
        Double_t       muTwoDSS14Pi02760GeVforGraph[39];
        Double_t       muTwoErrDSS14Pi02760GeV[39];
        Int_t          nlinesNLODSS14Pi02760GeV = 0;

        TString fileNameNLODSS14Pi02760GeV = "ExternalInput/Theory/pi0dss14-2760gev-dsigmadpt-eta060-pt40.dat";
        ifstream  fileNLODSS14Pi02760GeV;
        fileNLODSS14Pi02760GeV.open(fileNameNLODSS14Pi02760GeV,ios_base::in);
        cout << fileNameNLODSS14Pi02760GeV << endl;

        while(!fileNLODSS14Pi02760GeV.eof() && nlinesNLODSS14Pi02760GeV < 39){
            TString garbage;
            fileNLODSS14Pi02760GeV >> ptNLODSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] >> garbage >> muHalfDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV]  >> muTwoDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV];

            muOneDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] = (muHalfDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] + muTwoDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV])/(2*2*TMath::Pi()*ptNLODSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV]);
            muOneErrDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] = (muHalfDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV]-muTwoDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV])/(2*2*TMath::Pi()*ptNLODSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV]);
            ptErrNLODSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] = 0.;
            cout << nlinesNLODSS14Pi02760GeV << "\t "  << ptNLODSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] << "\t "  << muHalfDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] << "\t "  << muOneDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] << "         "  << muTwoDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] << endl;;

            muHalfDSS14Pi02760GeVforGraph[nlinesNLODSS14Pi02760GeV] = muHalfDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV]/(2*TMath::Pi()*ptNLODSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV]);
            muTwoDSS14Pi02760GeVforGraph[nlinesNLODSS14Pi02760GeV] = muTwoDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV]/(2*TMath::Pi()*ptNLODSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV]);

            muHalfErrDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] = 0.;
            muTwoErrDSS14Pi02760GeV[nlinesNLODSS14Pi02760GeV] = 0.;
            nlinesNLODSS14Pi02760GeV++;

        }
        fileNLODSS14Pi02760GeV.close();

        TGraphAsymmErrors* graphNLOCalcmuHalfDSS14InvSecPi02760GeV = new TGraphAsymmErrors(nlinesNLODSS14Pi02760GeV, ptNLODSS14Pi02760GeV, muHalfDSS14Pi02760GeVforGraph, ptErrNLODSS14Pi02760GeV,
                                                                                        ptErrNLODSS14Pi02760GeV, muHalfErrDSS14Pi02760GeV, muHalfErrDSS14Pi02760GeV);
        TGraphAsymmErrors* graphNLOCalcmuTwoDSS14InvSecPi02760GeV = new TGraphAsymmErrors(nlinesNLODSS14Pi02760GeV, ptNLODSS14Pi02760GeV, muTwoDSS14Pi02760GeVforGraph, ptErrNLODSS14Pi02760GeV,
                                                                                        ptErrNLODSS14Pi02760GeV,muTwoErrDSS14Pi02760GeV, muTwoErrDSS14Pi02760GeV);
        TGraphAsymmErrors* graphNLOCalcDSS14InvSecPi02760GeV = new TGraphAsymmErrors(nlinesNLODSS14Pi02760GeV, ptNLODSS14Pi02760GeV, muOneDSS14Pi02760GeV, ptErrNLODSS14Pi02760GeV, ptErrNLODSS14Pi02760GeV,
                                                                                    muOneErrDSS14Pi02760GeV, muOneErrDSS14Pi02760GeV);
        TGraphAsymmErrors* graphNLOCalcDSS14InvYieldPi02760GeV = ScaleGraph(graphNLOCalcDSS14InvSecPi02760GeV,  1./xSection2760GeVppINEL);


        // **************************************************************************************************************
        // ************************************* pQCD, pp 2.76TeV Eta DSS07 *********************************************
        // **************************************************************************************************************
        Double_t       ptNLODSS07Eta2760GeV[24];
        Double_t       ptErrNLODSS07Eta2760GeV[24];
        Double_t       muHalfDSS07Eta2760GeV[24];
        Double_t       muHalfDSS07Eta2760GeVforGraph[24];
        Double_t       muHalfErrDSS07Eta2760GeV[24];
        Double_t       muOneDSS07Eta2760GeV[24];
        Double_t       muOneErrDSS07Eta2760GeV[24];
        Double_t       muTwoDSS07Eta2760GeV[24];
        Double_t       muTwoDSS07Eta2760GeVforGraph[24];
        Double_t       muTwoErrDSS07Eta2760GeV[24];
        Int_t          nlinesNLODSS07Eta2760GeV = 0;

        TString fileNameNLODSS07Eta2760GeV = "ExternalInput/Theory/etadss07-2760gev-dsigmadpt-eta060-pt40.dat";
        ifstream  fileNLODSS07Eta2760GeV;
        fileNLODSS07Eta2760GeV.open(fileNameNLODSS07Eta2760GeV,ios_base::in);
        cout << fileNameNLODSS07Eta2760GeV << endl;

        while(!fileNLODSS07Eta2760GeV.eof() && nlinesNLODSS07Eta2760GeV < 24){
            TString garbage;
            fileNLODSS07Eta2760GeV >> ptNLODSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV] >> garbage >> muHalfDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV]  >> muTwoDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV];

            muOneDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV] = (muHalfDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV]+muTwoDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV])/ (2*2*TMath::Pi()*ptNLODSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV]);
            muOneErrDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV] = (muHalfDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV]-muTwoDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV])/(2*2*TMath::Pi()*ptNLODSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV]);
            ptErrNLODSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV] = 0.;
            cout << nlinesNLODSS07Eta2760GeV << "\t "  << ptNLODSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV] << "\t "  << muHalfDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV] << "\t "  << muOneDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV] << "\t "  << muTwoDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV] << endl;;

            muHalfDSS07Eta2760GeVforGraph[nlinesNLODSS07Eta2760GeV] = muHalfDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV]/(2*TMath::Pi()*ptNLODSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV]);
            muTwoDSS07Eta2760GeVforGraph[nlinesNLODSS07Eta2760GeV] = muTwoDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV]/(2*TMath::Pi()*ptNLODSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV]);
            muHalfErrDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV] = 0.;
            muTwoErrDSS07Eta2760GeV[nlinesNLODSS07Eta2760GeV] = 0.;
            nlinesNLODSS07Eta2760GeV++;

        }
        fileNLODSS07Eta2760GeV.close();

        TGraphAsymmErrors* graphNLOCalcmuHalfDSS07InvSecEta2760GeV = new TGraphAsymmErrors(nlinesNLODSS07Eta2760GeV, ptNLODSS07Eta2760GeV, muHalfDSS07Eta2760GeVforGraph, ptErrNLODSS07Eta2760GeV,
                                                                                        ptErrNLODSS07Eta2760GeV, muHalfErrDSS07Eta2760GeV, muHalfErrDSS07Eta2760GeV);
        TGraphAsymmErrors* graphNLOCalcmuTwoDSS07InvSecEta2760GeV = new TGraphAsymmErrors(nlinesNLODSS07Eta2760GeV, ptNLODSS07Eta2760GeV, muTwoDSS07Eta2760GeVforGraph, ptErrNLODSS07Eta2760GeV,
                                                                                        ptErrNLODSS07Eta2760GeV, muTwoErrDSS07Eta2760GeV, muTwoErrDSS07Eta2760GeV);
        TGraphAsymmErrors* graphNLOCalcDSS07InvSecEta2760GeV = new TGraphAsymmErrors(nlinesNLODSS07Eta2760GeV, ptNLODSS07Eta2760GeV, muOneDSS07Eta2760GeV, ptErrNLODSS07Eta2760GeV,
                                                                                    ptErrNLODSS07Eta2760GeV, muOneErrDSS07Eta2760GeV, muOneErrDSS07Eta2760GeV);
        TGraphAsymmErrors* graphNLOCalcDSS07InvYieldEta2760GeV = ScaleGraph(graphNLOCalcDSS07InvSecEta2760GeV, 1./xSection2760GeVppINEL);


        fileTheoryGraphsPbPb.cd("Pi0_PbPb_2.76TeV");
            //EPOS2
            graphEPOSPi0_0010->Write("graphEPOS_Pi0_0010");
            graphEPOSPi0_2050->Write("graphEPOS_Pi0_2050");
            //EPOS3
            graphEPOS3Pi0_0010->Write("graphEpos3_pi0_pt_cent0_10");

            graphNLOCalcmuHalfDSS14InvSecPi02760GeV->Write("graphNLOCalcmuHalfDSS14InvSecPi02760GeV");
            graphNLOCalcmuTwoDSS14InvSecPi02760GeV->Write("graphNLOCalcmuTwoDSS14InvSecPi02760GeV");
            graphNLOCalcDSS14InvSecPi02760GeV->Write("graphNLOCalcDSS14InvCrossSecPi02760GeV");
            graphNLOCalcDSS14InvYieldPi02760GeV->Write("graphNLOCalcDSS14InvYieldPi02760GeV");

//             TheoryCracowNeutralToChargedPionLowPt_0010->Write("TheoryCracowNeutralToChargedPionLowPt_0010");

        fileTheoryGraphsPbPb.cd("Eta_PbPb_2.76TeV");

            //EPOS2
            graphEPOSEta_0010->Write("graphEPOS_Eta_0010");
            graphEPOSEta_2050->Write("graphEPOS_Eta_2050");
            //EPOS3
            graphEPOS3Eta_0010->Write("graphEpos3_eta_pt_cent0_10");

            graphNLOCalcmuHalfDSS07InvSecEta2760GeV->Write("graphNLOCalcmuHalfDSS07InvSecEta2760GeV");
            graphNLOCalcmuTwoDSS07InvSecEta2760GeV->Write("graphNLOCalcmuTwoDSS07InvSecEta2760GeV");
            graphNLOCalcDSS07InvSecEta2760GeV->Write("graphNLOCalcDSS07InvCrossSecEta2760GeV");
            graphNLOCalcDSS07InvYieldEta2760GeV->Write("graphNLOCalcDSS07InvYieldEta2760GeV");

            graphEPOSEtaToPi0Ratio_0010->Write("graphEPOS_EtaToPi0_0010");
            graphEPOSEtaToPi0Ratio_2050->Write("graphEPOS_EtaToPi0_2050");



        fileTheoryGraphsPbPb.Close();
    } else if(energy.CompareTo("PbPb_5.02TeV")==0){

        TString outputFileName = Form("ExternalInputPbPb/Theory/TheoryCompilationPbPb.root");
        TFile fileTheoryGraphsPbPb(outputFileName,"UPDATE");
        cout << "INFO: You have chosen the following output file: " << outputFileName.Data() << endl;

        TDirectory *directory5TeVPi0    = (TDirectory*)fileTheoryGraphsPbPb.Get("Pi0_PbPb_5.02TeV");
        if (!directory5TeVPi0){
            fileTheoryGraphsPbPb.mkdir("Pi0_PbPb_5.02TeV");
        }
        TDirectory *directory5TeVEta    = (TDirectory*)fileTheoryGraphsPbPb.Get("Eta_PbPb_5.02TeV");
        if (!directory5TeVEta){
            fileTheoryGraphsPbPb.mkdir("Eta_PbPb_5.02TeV");
        }
        TDirectory *directory5TeVPiCh    = (TDirectory*)fileTheoryGraphsPbPb.Get("PiCh_PbPb_5.02TeV");
        if (!directory5TeVPiCh){
            fileTheoryGraphsPbPb.mkdir("PiCh_PbPb_5.02TeV");
        }
        TDirectory *directory5TeVKCh    = (TDirectory*)fileTheoryGraphsPbPb.Get("KCh_PbPb_5.02TeV");
        if (!directory5TeVKCh){
            fileTheoryGraphsPbPb.mkdir("KCh_PbPb_5.02TeV");
        }

        // Paquet Pi0 and Eta spectra at midrapidity
        TString filePaquetPi0 = "ExternalInputPbPb/Theory/McGill_5TeV/spectra_pi0_5020GeV_PbPb_Paquet.dat";
        TString filePaquetEta = "ExternalInputPbPb/Theory/McGill_5TeV/spectra_etas_5020GeV_PbPb_Paquet.dat";
        readFilePaquet(filePaquetPi0, fileTheoryGraphsPbPb, "Pi0_PbPb_5.02TeV", "Spectra_Paquet", "#it{p}_{T} (GeV/#it{c})", "d#it{N}/(2#pi d#it{y} #it{p}_{T} d#it{p}_{T})" ,kFALSE);
        readFilePaquet(filePaquetEta, fileTheoryGraphsPbPb, "Eta_PbPb_5.02TeV", "Spectra_Paquet", "#it{p}_{T} (GeV/#it{c})", "d#it{N}/(2#pi d#it{y} #it{p}_{T} d#it{p}_{T})" ,kFALSE);

        TString centArrayPaquet[11]    = {"0005", "0510", "1020", "2030", "3040", "4050", "5060", "0010", "0020", "2040", "4060"};
        for (Int_t i = 0; i < 11; i++){
            TGraph* graphPi0Curr    = (TGraph*)fileTheoryGraphsPbPb.Get(Form("Pi0_PbPb_5.02TeV/Spectra_Paquet_%s", centArrayPaquet[i].Data()));
            TGraph* graphEtaCurr    = (TGraph*)fileTheoryGraphsPbPb.Get(Form("Eta_PbPb_5.02TeV/Spectra_Paquet_%s", centArrayPaquet[i].Data()));
            if (graphPi0Curr && graphEtaCurr){
                TGraph* etaToPi0RatioCurr   = CalculateGraphRatioToGraph(graphEtaCurr, graphPi0Curr);
                fileTheoryGraphsPbPb.cd("Eta_PbPb_5.02TeV");
                etaToPi0RatioCurr->GetYaxis()->SetTitle("#eta/#pi^{0}");
                etaToPi0RatioCurr->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
                etaToPi0RatioCurr->Write(Form("EtaToPi0_Paquet_%s", centArrayPaquet[i].Data()));
            }
        }
        // Begun Pi0 and Eta spectra
        TString fileBegunEQ  = "ExternalInputPbPb/Theory/CracowModel_5TeV/Begun_EQ.txt";
        TString fileBegunNEQ = "ExternalInputPbPb/Theory/CracowModel_5TeV/Begun_NEQ.txt";
        readFileBegun(fileBegunEQ, fileTheoryGraphsPbPb,  "Pi0_PbPb_5.02TeV", "Eta_PbPb_5.02TeV", "Spectra_SHM_EQ", "EtaToPi0_SHM_EQ", "#it{p}_{T} (GeV/#it{c})", "d#it{N}/(2#pi d#it{y} #it{p}_{T} d#it{p}_{T})", "#eta/#pi^{0}", kFALSE);
        readFileBegun(fileBegunNEQ, fileTheoryGraphsPbPb,  "Pi0_PbPb_5.02TeV", "Eta_PbPb_5.02TeV", "Spectra_SHM_NEQ", "EtaToPi0_SHM_NEQ", "#it{p}_{T} (GeV/#it{c})", "d#it{N}/(2#pi d#it{y} #it{p}_{T} d#it{p}_{T})", "#eta/#pi^{0}", kFALSE);

        // Djordjevic Pi0 RAA  (one file per cent class and model)
        TString fileDjordjevicBjorken0010 = "ExternalInputPbPb/Theory/Djordjevic_5TeV/PionRaa_5TeV_Bjorken_010.txt";
        TString fileDjordjevicBjorken1020 = "ExternalInputPbPb/Theory/Djordjevic_5TeV/PionRaa_5TeV_Bjorken_1020.txt";
        TString fileDjordjevicBjorken2040 = "ExternalInputPbPb/Theory/Djordjevic_5TeV/PionRaa_5TeV_Bjorken_2040.txt";
        TString fileDjordjevicBjorken4060 = "ExternalInputPbPb/Theory/Djordjevic_5TeV/PionRaa_5TeV_Bjorken_4060.txt";
        TString fileDjordjevicBjorken6080 = "ExternalInputPbPb/Theory/Djordjevic_5TeV/PionRaa_5TeV_Bjorken_6080.txt";
        TString fileDjordjevicConstTemp0010 = "ExternalInputPbPb/Theory/Djordjevic_5TeV/PionRaa_5TeV_ConstTemp_010.txt";
        TString fileDjordjevicConstTemp1020 = "ExternalInputPbPb/Theory/Djordjevic_5TeV/PionRaa_5TeV_ConstTemp_1020.txt";
        TString fileDjordjevicConstTemp2040 = "ExternalInputPbPb/Theory/Djordjevic_5TeV/PionRaa_5TeV_ConstTemp_2040.txt";
        TString fileDjordjevicConstTemp4060 = "ExternalInputPbPb/Theory/Djordjevic_5TeV/PionRaa_5TeV_ConstTemp_4060.txt";
        TString fileDjordjevicConstTemp6080 = "ExternalInputPbPb/Theory/Djordjevic_5TeV/PionRaa_5TeV_ConstTemp_6080.txt";
        readFileDjordjevic(fileDjordjevicBjorken0010, fileTheoryGraphsPbPb, "Pi0_PbPb_5.02TeV", "RAA_DjordjevicBjorken_0010", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kFALSE);
        readFileDjordjevic(fileDjordjevicBjorken1020, fileTheoryGraphsPbPb, "Pi0_PbPb_5.02TeV", "RAA_DjordjevicBjorken_1020", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kFALSE);
        readFileDjordjevic(fileDjordjevicBjorken2040, fileTheoryGraphsPbPb, "Pi0_PbPb_5.02TeV", "RAA_DjordjevicBjorken_2040", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kFALSE);
        readFileDjordjevic(fileDjordjevicBjorken4060, fileTheoryGraphsPbPb, "Pi0_PbPb_5.02TeV", "RAA_DjordjevicBjorken_4060", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kFALSE);
        readFileDjordjevic(fileDjordjevicBjorken6080, fileTheoryGraphsPbPb, "Pi0_PbPb_5.02TeV", "RAA_DjordjevicBjorken_6080", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kFALSE);
        readFileDjordjevic(fileDjordjevicConstTemp0010, fileTheoryGraphsPbPb, "Pi0_PbPb_5.02TeV", "RAA_DjordjevicConstTemp_0010", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kFALSE);
        readFileDjordjevic(fileDjordjevicConstTemp1020, fileTheoryGraphsPbPb, "Pi0_PbPb_5.02TeV", "RAA_DjordjevicConstTemp_1020", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kFALSE);
        readFileDjordjevic(fileDjordjevicConstTemp2040, fileTheoryGraphsPbPb, "Pi0_PbPb_5.02TeV", "RAA_DjordjevicConstTemp_2040", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kFALSE);
        readFileDjordjevic(fileDjordjevicConstTemp4060, fileTheoryGraphsPbPb, "Pi0_PbPb_5.02TeV", "RAA_DjordjevicConstTemp_4060", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kFALSE);
        readFileDjordjevic(fileDjordjevicConstTemp6080, fileTheoryGraphsPbPb, "Pi0_PbPb_5.02TeV", "RAA_DjordjevicConstTemp_6080", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}", kFALSE);

        // Vitev Pi0 RAA  (two files, min and max, per cent class)
        TString fileVitevMin0010 = "ExternalInputPbPb/Theory/Vitev_5TeV/R-N.aa_010_cron1.5_eloss0.5100GeVpi.g2.1";
        TString fileVitevMin2040 = "ExternalInputPbPb/Theory/Vitev_5TeV/R-N.aa_2040_cron1.5_eloss0.5100GeVpi.g2.1";
        TString fileVitevMin4060 = "ExternalInputPbPb/Theory/Vitev_5TeV/R-N.aa_4060_cron1.5_eloss0.5100GeVpi.g2.1";
        TString fileVitevMax0010 = "ExternalInputPbPb/Theory/Vitev_5TeV/R-N.aa_010_cron1.5_eloss0.5100GeVpi.g1.9";
        TString fileVitevMax2040 = "ExternalInputPbPb/Theory/Vitev_5TeV/R-N.aa_2040_cron1.5_eloss0.5100GeVpi.g1.9";
        TString fileVitevMax4060 = "ExternalInputPbPb/Theory/Vitev_5TeV/R-N.aa_4060_cron1.5_eloss0.5100GeVpi.g1.9";
        readFileMinAndMaxSeparately(fileVitevMin0010, fileVitevMax0010, fileTheoryGraphsPbPb, "Pi0_PbPb_5.02TeV", "RAA_Vitev_0010", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}" , kFALSE);
        readFileMinAndMaxSeparately(fileVitevMin2040, fileVitevMax2040, fileTheoryGraphsPbPb, "Pi0_PbPb_5.02TeV", "RAA_Vitev_2040", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}" , kFALSE);
        readFileMinAndMaxSeparately(fileVitevMin4060, fileVitevMax4060, fileTheoryGraphsPbPb, "Pi0_PbPb_5.02TeV", "RAA_Vitev_4060", "#it{p}_{T} (GeV/#it{c})", "#it{R}_{AA}" , kFALSE);


        TString centArrayReadHIJING[4]  = {"0-10%", "10-30%", "30-50%", "50-90%"};
        TString centArrayOutHIJING[4]   = {"0010", "1030", "3050", "5090"};
        TFile* fileHIJING               = new TFile("ExternalInputPbPb/Theory/MCInputCompilationLHC18e1_PbPb5TeV_0.root");
        for (Int_t j = 0; j<4; j++){
            TString hijingDir               = Form("%sPbPb_5.02TeV",centArrayReadHIJING[j].Data());
            TH1D* histoPi0HIJING            = (TH1D*)fileHIJING->Get(Form("%s/MC_Pi0_Pt", hijingDir.Data()));
            TH1D* histoPi0ToPiCHHIJING      = (TH1D*)fileHIJING->Get(Form("%s/MCPi0ToPiCh", hijingDir.Data()));
            directory5TeVPi0->cd();
                histoPi0HIJING->Write(Form("Spectra_HIJING_MCGen_%s",centArrayOutHIJING[j].Data()), TObject::kOverwrite);
                histoPi0ToPiCHHIJING->Write(Form("Pi0ToPiCh_HIJING_MCGen_%s",centArrayOutHIJING[j].Data()), TObject::kOverwrite);
//             TH1D* histoPi0HIJINGReb         = (TH1D*)fileHIJING->Get("MC_Pi0_Pt_Rebinned");
            TH1D* histoPiChHIJING           = (TH1D*)fileHIJING->Get(Form("%s/MC_PiCh_All_Pt", hijingDir.Data()));
            directory5TeVPiCh->cd();
                histoPiChHIJING->Write(Form("Spectra_HIJING_MCGen_%s",centArrayOutHIJING[j].Data()), TObject::kOverwrite);

            TH1D* histoKChHIJING            = (TH1D*)fileHIJING->Get(Form("%s/MC_KCh_All_Pt", hijingDir.Data()));
            directory5TeVKCh->cd();
                histoKChHIJING->Write(Form("Spectra_HIJING_MCGen_%s",centArrayOutHIJING[j].Data()), TObject::kOverwrite);

            TH1D* histoEtaHIJING            = (TH1D*)fileHIJING->Get(Form("%s/MC_Eta_Pt", hijingDir.Data()));
            TH1D* histoEtaToPi0HIJING       = (TH1D*)fileHIJING->Get(Form("%s/MCEtaToPi0", hijingDir.Data()));
            TH1D* histoEtaToKChHIJING       = (TH1D*)fileHIJING->Get(Form("%s/MCEtaToKCh", hijingDir.Data()));
            directory5TeVEta->cd();
            histoEtaHIJING->Write(Form("Spectra_HIJING_MCGen_%s",centArrayOutHIJING[j].Data()), TObject::kOverwrite);
                histoEtaToPi0HIJING->Write(Form("EtaToPi0_HIJING_MCGen_%s",centArrayOutHIJING[j].Data()), TObject::kOverwrite);
                histoEtaToKChHIJING->Write(Form("EtaToKCh_HIJING_MCGen_%s",centArrayOutHIJING[j].Data()), TObject::kOverwrite);

            //             TH1D* histoEtaHIJINGReb         = (TH1D*)fileHIJING->Get("MC_Eta_Pt_Rebinned");

        }

    }
}
