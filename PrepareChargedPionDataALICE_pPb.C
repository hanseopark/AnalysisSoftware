/****************************************************************************************************************************
******      provided by Gamma Conversion Group, PWGGA,                                                                  *****
******      Friederike Bock, friederike.bock@cern.ch                                                                    *****
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
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"

extern TRandom* gRandom;
extern TBenchmark* gBenchmark;
extern TSystem* gSystem;
extern TMinuit* gMinuit;

void RemoveZerosAndPrint(TGraphAsymmErrors* graph){
  while (graph->GetY()[0] == 0. ) graph->RemovePoint(0);
  while (graph->GetY()[graph->GetN()-1] == 0. ) graph->RemovePoint(graph->GetN()-1);
  graph->Print();
}

TGraphAsymmErrors* ParseHEPData(TString hepDataFile,
                                Int_t   totalNumberOfColumns,
                                Int_t   columnX,
                                Int_t   columnXErrLow,
                                Int_t   columnXErrHigh,
                                Int_t   columnY,
                                Int_t   columnYErrLow,
                                Int_t   columnYErrHigh,
                                Bool_t  isXErrVal,
                                Bool_t  isYErrVal,
                                Bool_t  debugMode = kFALSE) {

    // create streamer
    ifstream file;
    if (debugMode) cout << "HEP data file: " << hepDataFile.Data() << endl;
    file.open(hepDataFile,ios_base::in);
    if (!file) {
        cout << "ERROR: HEP data file " << hepDataFile.Data() << " not found!" << endl;
        return NULL;
    }

    // check for correct column numbers
    if (columnX<0) {
        cout << "ERROR: columnX set to " << columnX << endl;
        return NULL;
    }
    if (columnY<0) {
        cout << "ERROR: columnY set to " << columnY << endl;
        return NULL;
    }
    if (columnYErrLow<0 || columnYErrHigh<0) {
        cout << "ERROR: columnYErrLow set to " << columnYErrLow << " and columnYErrHigh set to " << columnYErrHigh << endl;
        return NULL;
    }

    // initialize vectors for temporary storage of values
    std::vector<Double_t> xVal;
    std::vector<Double_t> xErrLow;
    std::vector<Double_t> xErrHigh;
    std::vector<Double_t> yVal;
    std::vector<Double_t> yErrLow;
    std::vector<Double_t> yErrHigh;

    // read from file
    TString                 tempString;
    std::vector<TString>    tempStringColumn(totalNumberOfColumns);
    std::string line;
    for( std::string line; getline(file, line); ) {
        file >> tempString;
        if (!tempString.BeginsWith("%") && !tempString.BeginsWith("%") && tempString.CompareTo("")) {
            tempStringColumn[0]     = tempString;
            if (debugMode) cout << tempStringColumn[0].Data() << "\t";
            for (Int_t i=1; i<totalNumberOfColumns; i++) {
                file >> tempStringColumn[i];
                if (debugMode) cout << tempStringColumn[i].Data() << "\t";
            }
            if (debugMode) cout << endl;

            // x value and error
            xVal.push_back(tempStringColumn[columnX].Atof());
            if (columnXErrLow>=0)   xErrLow.push_back(tempStringColumn[columnXErrLow].Atof());
            else                    xErrLow.push_back(-1);
            if (columnXErrHigh>=0)  xErrHigh.push_back(tempStringColumn[columnXErrHigh].Atof());
            else                    xErrHigh.push_back(-1);

            // y value and error
            yVal.push_back(tempStringColumn[columnY].Atof());
            yErrLow.push_back(tempStringColumn[columnYErrLow].Atof());
            yErrHigh.push_back(tempStringColumn[columnYErrHigh].Atof());
        } else
            continue;
    }

    // check for equal number of rows for each column
    Bool_t  isEqualNumberOfRows     = kTRUE;
    Int_t   nRowsTemp[6];
    nRowsTemp[0]                    = xVal.size();
    nRowsTemp[1]                    = xErrLow.size();
    nRowsTemp[2]                    = xErrHigh.size();
    nRowsTemp[3]                    = yVal.size();
    nRowsTemp[4]                    = yErrLow.size();
    nRowsTemp[5]                    = yErrHigh.size();
    for (Int_t i=0; i<5; i++) {
        if (nRowsTemp[i]!=nRowsTemp[i+1]) {
            isEqualNumberOfRows     = kFALSE;
            break;
        }
    }
    if (!isEqualNumberOfRows) {
        cout << "number of rows in " << hepDataFile.Data() << " are not equal for different columns!" << endl;
        return NULL;
    }
    Int_t nRows                     = xVal.size();

    // calculate x errors if necessary (i.e. column numbers set to -1)
    std::vector<Double_t> tempXErr(xVal.size());
    if (columnXErrLow<0 || columnXErrHigh<0) {
        for (Int_t i=0; i<nRows; i++) {

            // calculate x error
            if (i==0)               tempXErr[i] = (xVal[1]-xVal[0])/2;
            else if (i==nRows-1)    tempXErr[i] = xVal[i]-(xVal[i-1] + tempXErr[i-1]);
            else                    tempXErr[i] = (xVal[i]-xVal[i-1])/2;

            // set error
            xErrLow[i]              = tempXErr[i];
            xErrHigh[i]             = tempXErr[i];
        }
    }

    // calculate errors if bin boundaries were given
    if (!isXErrVal && columnXErrLow>=0 && columnXErrHigh>=0) {
        for (Int_t i=0; i<nRows; i++) {
            xErrLow[i]              = TMath::Abs(xVal[i]-xErrLow[i]);
            xErrHigh[i]             = TMath::Abs(xErrHigh[i]-xVal[i]);
        }
    }
    if (!isYErrVal) {
        for (Int_t i=0; i<nRows; i++) {
            yErrLow[i]              = TMath::Abs(yVal[i]-yErrLow[i]);
            yErrHigh[i]             = TMath::Abs(yErrHigh[i]-yVal[i]);
        }
    }

    // set errors to absolute values, direction is taken care of by TGraphAsymmErrors
    for (Int_t i=0; i<nRows; i++) {
        xErrLow[i]                  = TMath::Abs(xErrLow[i]);
        xErrHigh[i]                 = TMath::Abs(xErrHigh[i]);

        yErrLow[i]                  = TMath::Abs(yErrLow[i]);
        yErrHigh[i]                 = TMath::Abs(yErrHigh[i]);
    }

    // cout values (debug mode)
    if (debugMode) {
        cout << "nRows = " << nRows << endl;
        for (Int_t i=0; i<nRows; i++) {
            cout << "x = " << xVal[i] << "\t+ " << xErrHigh[i] << "\t- " << xErrLow[i] << "\t y = " << yVal[i] << "\t+ " << yErrHigh[i] << "\t- " << yErrLow[i] << endl;
        }
    }

    // create TGraphAsymmErrors
    TGraphAsymmErrors* graph        = new TGraphAsymmErrors(nRows);
    for (Int_t i=0; i<nRows; i++) {
        graph->SetPoint(        i, xVal[i], yVal[i]);
        graph->SetPointError(   i, xErrLow[i], xErrHigh[i], TMath::Abs(yErrLow[i]), TMath::Abs(yErrHigh[i]));
    }
    return graph;
}


// -----------------------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------- Main function -----------------------------------------------------------------------------
// ----------------------- This macro is used to compile the pPb external input file for data, i.e. charged hadron, pions, kaons ... --------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------------------------------
void PrepareChargedPionDataALICE_pPb(){

    // *********************************************************************************************************************
    // ************************************** global variable definition ***************************************************
    // *********************************************************************************************************************
    TString dateForOutput                               = ReturnDateStringForOutput();

        // ************************** Low Pt  charged pions pPb************************************************
    TFile *fPionLowPtpPbMinBias = TFile::Open("ExternalInputpPb/20130201.CombinedSpectra_pA_pilot_ycms_0005_minbias.root");

    TH1D* histoChargedPionPlusSpecLowPtStatpPb = (TH1D*)fPionLowPtpPbMinBias->Get("stat_mb_pion_plus");
    TH1D* histoChargedPionPlusSpecLowPtSyspPb = (TH1D*)fPionLowPtpPbMinBias->Get("sys_mb_pion_plus");
    TH1D* histoChargedPionMinusSpecLowPtStatpPb = (TH1D*)fPionLowPtpPbMinBias->Get("stat_mb_pion_minus");
    TH1D* histoChargedPionMinusSpecLowPtSyspPb = (TH1D*)fPionLowPtpPbMinBias->Get("sys_mb_pion_minus");

    TH1D* histoChargedPionSpecLowPtStatpPb = (TH1D*)histoChargedPionMinusSpecLowPtStatpPb->Clone("histoChargedPionSpecLowPtStatpPb");
    histoChargedPionSpecLowPtStatpPb->Add(histoChargedPionPlusSpecLowPtStatpPb);
    histoChargedPionSpecLowPtStatpPb->Scale(0.5);

    TH1D* histoChargedPionSpecLowPtSyspPb = (TH1D*)histoChargedPionMinusSpecLowPtSyspPb->Clone("histoChargedPionSpecLowPtSyspPb");

    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyspPb->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStatpPb->SetBinContent(i, histoChargedPionSpecLowPtStatpPb->GetBinContent(i)/histoChargedPionSpecLowPtStatpPb->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStatpPb->SetBinError(i, i, histoChargedPionSpecLowPtStatpPb->GetBinError(i)/histoChargedPionSpecLowPtStatpPb->GetBinCenter(i)/(2*TMath::Pi()));
        Double_t error = 0;
        if (histoChargedPionMinusSpecLowPtSyspPb->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyspPb->GetBinContent(i) != 0 && histoChargedPionSpecLowPtStatpPb->GetBinContent(i) != 0){
            error= (histoChargedPionMinusSpecLowPtSyspPb->GetBinError(i)/histoChargedPionMinusSpecLowPtSyspPb->GetBinContent(i) +  histoChargedPionPlusSpecLowPtSyspPb->GetBinError(i) /histoChargedPionPlusSpecLowPtSyspPb->GetBinContent(i))/2.*histoChargedPionSpecLowPtStatpPb->GetBinContent(i);
        } else {
            error = 0;
        }
        histoChargedPionSpecLowPtSyspPb->SetBinContent(i, histoChargedPionSpecLowPtStatpPb->GetBinContent(i));
        histoChargedPionSpecLowPtSyspPb->SetBinError(i,error);
    }
    cout << "bis hier" << endl;
    // ************************** Low Pt  charged pions pPb depending on centrality************************************************
    TFile *fPionLowPtpPbCent = TFile::Open("ExternalInputpPb/20130201.CombinedSpectra_pA_pilot_ycms_0005.root");

    TH1D* histoChargedPionPlusSpecLowPtStatpPb0020 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent0_pion_plus");
    TH1D* histoChargedPionPlusSpecLowPtSyspPb0020 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent0_pion_plus");
    TH1D* histoChargedPionMinusSpecLowPtStatpPb0020 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent0_pion_minus");
    TH1D* histoChargedPionMinusSpecLowPtSyspPb0020 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent0_pion_minus");

    TH1D* histoChargedPionSpecLowPtStatpPb0020 = (TH1D*)histoChargedPionMinusSpecLowPtStatpPb0020->Clone("histoChargedPionSpecLowPtStatpPb0020");
    histoChargedPionSpecLowPtStatpPb0020->Add(histoChargedPionPlusSpecLowPtStatpPb0020);
    histoChargedPionSpecLowPtStatpPb0020->Scale(0.5);

    TH1D* histoChargedPionSpecLowPtSyspPb0020 = (TH1D*)histoChargedPionMinusSpecLowPtSyspPb0020->Clone("histoChargedPionSpecLowPtSyspPb0020");

    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyspPb0020->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStatpPb0020->SetBinContent(i, histoChargedPionSpecLowPtStatpPb0020->GetBinContent(i)/histoChargedPionSpecLowPtStatpPb0020->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStatpPb0020->SetBinError(i, i, histoChargedPionSpecLowPtStatpPb0020->GetBinError(i)/histoChargedPionSpecLowPtStatpPb0020->GetBinCenter(i)/(2*TMath::Pi()));
        Double_t error = 0;
        if (histoChargedPionMinusSpecLowPtSyspPb0020->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyspPb0020->GetBinContent(i) != 0 && histoChargedPionSpecLowPtStatpPb0020->GetBinContent(i) != 0){
            error= (histoChargedPionMinusSpecLowPtSyspPb0020->GetBinError(i)/histoChargedPionMinusSpecLowPtSyspPb0020->GetBinContent(i) +  histoChargedPionPlusSpecLowPtSyspPb0020->GetBinError(i) /histoChargedPionPlusSpecLowPtSyspPb0020->GetBinContent(i))/2.*histoChargedPionSpecLowPtStatpPb0020->GetBinContent(i);
        } else {
            error = 0;
        }
        histoChargedPionSpecLowPtSyspPb0020->SetBinContent(i, histoChargedPionSpecLowPtStatpPb0020->GetBinContent(i));
        histoChargedPionSpecLowPtSyspPb0020->SetBinError(i,error);
    }
    cout << "bis hier" << endl;
    TH1D* histoChargedPionPlusSpecLowPtStatpPb2040 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent1_pion_plus");
    TH1D* histoChargedPionPlusSpecLowPtSyspPb2040 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent1_pion_plus");
    TH1D* histoChargedPionMinusSpecLowPtStatpPb2040 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent1_pion_minus");
    TH1D* histoChargedPionMinusSpecLowPtSyspPb2040 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent1_pion_minus");

    TH1D* histoChargedPionSpecLowPtStatpPb2040 = (TH1D*)histoChargedPionMinusSpecLowPtStatpPb2040->Clone("histoChargedPionSpecLowPtStatpPb2040");
    histoChargedPionSpecLowPtStatpPb2040->Add(histoChargedPionPlusSpecLowPtStatpPb2040);
    histoChargedPionSpecLowPtStatpPb2040->Scale(0.5);

    TH1D* histoChargedPionSpecLowPtSyspPb2040 = (TH1D*)histoChargedPionMinusSpecLowPtSyspPb2040->Clone("histoChargedPionSpecLowPtSyspPb2040");

    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyspPb2040->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStatpPb2040->SetBinContent(i, histoChargedPionSpecLowPtStatpPb2040->GetBinContent(i)/histoChargedPionSpecLowPtStatpPb2040->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStatpPb2040->SetBinError(i, i, histoChargedPionSpecLowPtStatpPb2040->GetBinError(i)/histoChargedPionSpecLowPtStatpPb2040->GetBinCenter(i)/(2*TMath::Pi()));
        Double_t error = 0;
        if (histoChargedPionMinusSpecLowPtSyspPb2040->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyspPb2040->GetBinContent(i) != 0 && histoChargedPionSpecLowPtStatpPb2040->GetBinContent(i) != 0){
            error= (histoChargedPionMinusSpecLowPtSyspPb2040->GetBinError(i)/histoChargedPionMinusSpecLowPtSyspPb2040->GetBinContent(i) +  histoChargedPionPlusSpecLowPtSyspPb2040->GetBinError(i) /histoChargedPionPlusSpecLowPtSyspPb2040->GetBinContent(i))/2.*histoChargedPionSpecLowPtStatpPb2040->GetBinContent(i);
        } else {
            error = 0;
        }
        histoChargedPionSpecLowPtSyspPb2040->SetBinContent(i, histoChargedPionSpecLowPtStatpPb2040->GetBinContent(i));
        histoChargedPionSpecLowPtSyspPb2040->SetBinError(i,error);
    }
    TH1D* histoChargedPionPlusSpecLowPtStatpPb4060 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent2_pion_plus");
    TH1D* histoChargedPionPlusSpecLowPtSyspPb4060 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent2_pion_plus");
    TH1D* histoChargedPionMinusSpecLowPtStatpPb4060 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent2_pion_minus");
    TH1D* histoChargedPionMinusSpecLowPtSyspPb4060 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent2_pion_minus");

    TH1D* histoChargedPionSpecLowPtStatpPb4060 = (TH1D*)histoChargedPionMinusSpecLowPtStatpPb4060->Clone("histoChargedPionSpecLowPtStatpPb4060");
    histoChargedPionSpecLowPtStatpPb4060->Add(histoChargedPionPlusSpecLowPtStatpPb4060);
    histoChargedPionSpecLowPtStatpPb4060->Scale(0.5);

    TH1D* histoChargedPionSpecLowPtSyspPb4060 = (TH1D*)histoChargedPionMinusSpecLowPtSyspPb4060->Clone("histoChargedPionSpecLowPtSyspPb4060");

    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyspPb4060->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStatpPb4060->SetBinContent(i, histoChargedPionSpecLowPtStatpPb4060->GetBinContent(i)/histoChargedPionSpecLowPtStatpPb4060->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStatpPb4060->SetBinError(i, i, histoChargedPionSpecLowPtStatpPb4060->GetBinError(i)/histoChargedPionSpecLowPtStatpPb4060->GetBinCenter(i)/(2*TMath::Pi()));
        Double_t error = 0;
        if (histoChargedPionMinusSpecLowPtSyspPb4060->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyspPb4060->GetBinContent(i) != 0 && histoChargedPionSpecLowPtStatpPb4060->GetBinContent(i) != 0){
            error= (histoChargedPionMinusSpecLowPtSyspPb4060->GetBinError(i)/histoChargedPionMinusSpecLowPtSyspPb4060->GetBinContent(i) +  histoChargedPionPlusSpecLowPtSyspPb4060->GetBinError(i) /histoChargedPionPlusSpecLowPtSyspPb4060->GetBinContent(i))/2.*histoChargedPionSpecLowPtStatpPb4060->GetBinContent(i);
        } else {
            error = 0;
        }
        histoChargedPionSpecLowPtSyspPb4060->SetBinContent(i, histoChargedPionSpecLowPtStatpPb4060->GetBinContent(i));
        histoChargedPionSpecLowPtSyspPb4060->SetBinError(i,error);
    }
    cout << "bis hier" << endl;
    TH1D* histoChargedPionPlusSpecLowPtStatpPb6080 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent3_pion_plus");
    TH1D* histoChargedPionPlusSpecLowPtSyspPb6080 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent3_pion_plus");
    TH1D* histoChargedPionMinusSpecLowPtStatpPb6080 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent3_pion_minus");
    TH1D* histoChargedPionMinusSpecLowPtSyspPb6080 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent3_pion_minus");

    TH1D* histoChargedPionSpecLowPtStatpPb6080 = (TH1D*)histoChargedPionMinusSpecLowPtStatpPb6080->Clone("histoChargedPionSpecLowPtStatpPb6080");
    histoChargedPionSpecLowPtStatpPb6080->Add(histoChargedPionPlusSpecLowPtStatpPb6080);
    histoChargedPionSpecLowPtStatpPb6080->Scale(0.5);

    TH1D* histoChargedPionSpecLowPtSyspPb6080 = (TH1D*)histoChargedPionMinusSpecLowPtSyspPb6080->Clone("histoChargedPionSpecLowPtSyspPb6080");

    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyspPb6080->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStatpPb6080->SetBinContent(i, histoChargedPionSpecLowPtStatpPb6080->GetBinContent(i)/histoChargedPionSpecLowPtStatpPb6080->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStatpPb6080->SetBinError(i, i, histoChargedPionSpecLowPtStatpPb6080->GetBinError(i)/histoChargedPionSpecLowPtStatpPb6080->GetBinCenter(i)/(2*TMath::Pi()));
        Double_t error = 0;
        if (histoChargedPionMinusSpecLowPtSyspPb6080->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyspPb6080->GetBinContent(i) != 0 && histoChargedPionSpecLowPtStatpPb6080->GetBinContent(i) != 0){
            error= (histoChargedPionMinusSpecLowPtSyspPb6080->GetBinError(i)/histoChargedPionMinusSpecLowPtSyspPb6080->GetBinContent(i) +  histoChargedPionPlusSpecLowPtSyspPb6080->GetBinError(i) /histoChargedPionPlusSpecLowPtSyspPb6080->GetBinContent(i))/2.*histoChargedPionSpecLowPtStatpPb6080->GetBinContent(i);
        } else {
            error = 0;
        }
        histoChargedPionSpecLowPtSyspPb6080->SetBinContent(i, histoChargedPionSpecLowPtStatpPb6080->GetBinContent(i));
        histoChargedPionSpecLowPtSyspPb6080->SetBinError(i,error);
    }
    TH1D* histoChargedPionPlusSpecLowPtStatpPb80100 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent4_pion_plus");
    TH1D* histoChargedPionPlusSpecLowPtSyspPb80100 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent4_pion_plus");
    TH1D* histoChargedPionMinusSpecLowPtStatpPb80100 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent4_pion_minus");
    TH1D* histoChargedPionMinusSpecLowPtSyspPb80100 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent4_pion_minus");

    TH1D* histoChargedPionSpecLowPtStatpPb80100 = (TH1D*)histoChargedPionMinusSpecLowPtStatpPb80100->Clone("histoChargedPionSpecLowPtStatpPb80100");
    histoChargedPionSpecLowPtStatpPb80100->Add(histoChargedPionPlusSpecLowPtStatpPb80100);
    histoChargedPionSpecLowPtStatpPb80100->Scale(0.5);

    TH1D* histoChargedPionSpecLowPtSyspPb80100 = (TH1D*)histoChargedPionMinusSpecLowPtSyspPb80100->Clone("histoChargedPionSpecLowPtSyspPb80100");

    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyspPb80100->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStatpPb80100->SetBinContent(i, histoChargedPionSpecLowPtStatpPb80100->GetBinContent(i)/histoChargedPionSpecLowPtStatpPb80100->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStatpPb80100->SetBinError(i, i, histoChargedPionSpecLowPtStatpPb80100->GetBinError(i)/histoChargedPionSpecLowPtStatpPb80100->GetBinCenter(i)/(2*TMath::Pi()));
        Double_t error = 0;
        if (histoChargedPionMinusSpecLowPtSyspPb80100->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyspPb80100->GetBinContent(i) != 0 && histoChargedPionSpecLowPtStatpPb80100->GetBinContent(i) != 0){
            error= (histoChargedPionMinusSpecLowPtSyspPb80100->GetBinError(i)/histoChargedPionMinusSpecLowPtSyspPb80100->GetBinContent(i) +  histoChargedPionPlusSpecLowPtSyspPb80100->GetBinError(i) /histoChargedPionPlusSpecLowPtSyspPb80100->GetBinContent(i))/2.*histoChargedPionSpecLowPtStatpPb80100->GetBinContent(i);
        } else {
            error = 0;
        }
        histoChargedPionSpecLowPtSyspPb80100->SetBinContent(i, histoChargedPionSpecLowPtStatpPb80100->GetBinContent(i));
        histoChargedPionSpecLowPtSyspPb80100->SetBinError(i,error);
    }

    // **************************************************************************************
    // ******************************** Reading Charged Pions *******************************
    // publication:
    // twiki: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSPECTRALowToHighPtpA
    // **************************************************************************************
    TFile* fileIdentfiedPub                         = new TFile("ExternalInputpPb/InputRpPb/pPb502.fullpT.INEL.20151204_mb_wo_V0Acorr.root");
    TH1D* histoChargedPionspPbStatErr               = (TH1D*)fileIdentfiedPub->Get("hstat_pPb502_mb_pion_sum");
    histoChargedPionspPbStatErr->Sumw2();
    histoChargedPionspPbStatErr->Scale(0.5);
    TH1D* histoChargedPionspPbSystErr               = (TH1D*)fileIdentfiedPub->Get("hsys_pPb502_mb_pion_sum");
    histoChargedPionspPbSystErr->Sumw2();
    histoChargedPionspPbSystErr->Scale(0.5);
    TGraphAsymmErrors* graphChargedPionspPbSyst     = new TGraphAsymmErrors(histoChargedPionspPbStatErr);
    TGraphAsymmErrors* graphChargedPionspPbStat     = new TGraphAsymmErrors(histoChargedPionspPbSystErr);
    TH1D* histoChargedKaonspPbStatErr               = (TH1D*)fileIdentfiedPub->Get("hstat_pPb502_mb_kaon_sum");
    histoChargedKaonspPbStatErr->Sumw2();
    histoChargedKaonspPbStatErr->Scale(0.5);
    TH1D* histoChargedKaonspPbSystErr               = (TH1D*)fileIdentfiedPub->Get("hsys_pPb502_mb_kaon_sum");
    histoChargedKaonspPbSystErr->Sumw2();
    histoChargedKaonspPbSystErr->Scale(0.5);
    TGraphAsymmErrors* graphChargedKaonspPbSyst     = new TGraphAsymmErrors(histoChargedKaonspPbStatErr);
    TGraphAsymmErrors* graphChargedKaonspPbStat     = new TGraphAsymmErrors(histoChargedKaonspPbSystErr);
    TH1D* histoChargedProtonspPbStatErr             = (TH1D*)fileIdentfiedPub->Get("hstat_pPb502_mb_proton_sum");
    histoChargedProtonspPbStatErr->Sumw2();
    histoChargedProtonspPbStatErr->Scale(0.5);
    TH1D* histoChargedProtonspPbSystErr             = (TH1D*)fileIdentfiedPub->Get("hsys_pPb502_mb_pion_sum");
    histoChargedProtonspPbSystErr->Sumw2();
    histoChargedProtonspPbSystErr->Scale(0.5);
    TGraphAsymmErrors* graphChargedProtonspPbSyst   = new TGraphAsymmErrors(histoChargedProtonspPbSystErr);
    TGraphAsymmErrors* graphChargedProtonspPbStat   = new TGraphAsymmErrors(histoChargedProtonspPbStatErr);

    TH1D* histoChargedPionsStatCent[7]              = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoChargedPionsSystCent[7]              = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoChargedKaonsStatCent[7]              = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoChargedKaonsSystCent[7]              = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoChargedProtonsStatCent[7]            = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoChargedProtonsSystCent[7]            = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoChargedKaonToPionStatCent[9]         = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoChargedKaonToPionSystCent[9]         = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoChargedProtonToPionStatCent[9]       = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoChargedProtonToPionSystCent[9]       = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphChargedPionsStatCent[7]              = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphChargedPionsSystCent[7]              = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphChargedKaonsStatCent[7]              = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphChargedKaonsSystCent[7]              = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphChargedProtonsStatCent[7]            = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphChargedProtonsSystCent[7]            = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphChargedKaonToPionStatCent[9]         = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphChargedKaonToPionStatCentW0XErr[9]   = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphChargedKaonToPionSystCent[9]         = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphChargedKaonToPionTotCent[9]          = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphChargedProtonToPionStatCent[9]       = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphChargedProtonToPionStatCentW0XErr[9] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphChargedProtonToPionSystCent[9]       = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphChargedProtonToPionTotCent[9]        = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TString centStringCharged[7]                    = {"0005","0510","1020","2040","4060","6080","80100"};
    TFile* fileIdentfiedPubCents                    = new TFile("ExternalInputpPb/IdentifiedCharged/pPb502.fullpT.INEL.20151204.root");
    TFile* fileIdentfiedPubCentsRatios              = new TFile("ExternalInputpPb/IdentifiedCharged/pPb502.fullpT.RATIOS.20151204.root");
    for (Int_t i = 0; i<7; i++){
        histoChargedPionsStatCent[i]                = (TH1D*)fileIdentfiedPubCents->Get(Form("hstat_pPb502_%s_pion_sum",centStringCharged[i].Data()));
        histoChargedPionsStatCent[i]->Sumw2();
        histoChargedPionsStatCent[i]->Scale(0.5);
        histoChargedPionsSystCent[i]                = (TH1D*)fileIdentfiedPubCents->Get(Form("hsys_pPb502_%s_pion_sum",centStringCharged[i].Data()));
        histoChargedPionsSystCent[i]->Sumw2();
        histoChargedPionsSystCent[i]->Scale(0.5);

        histoChargedKaonsStatCent[i]                = (TH1D*)fileIdentfiedPubCents->Get(Form("hstat_pPb502_%s_kaon_sum",centStringCharged[i].Data()));
        histoChargedKaonsStatCent[i]->Sumw2();
        histoChargedKaonsStatCent[i]->Scale(0.5);
        histoChargedKaonsSystCent[i]                = (TH1D*)fileIdentfiedPubCents->Get(Form("hsys_pPb502_%s_kaon_sum",centStringCharged[i].Data()));
        histoChargedKaonsSystCent[i]->Sumw2();
        histoChargedKaonsSystCent[i]->Scale(0.5);

        histoChargedProtonsStatCent[i]              = (TH1D*)fileIdentfiedPubCents->Get(Form("hstat_pPb502_%s_proton_sum",centStringCharged[i].Data()));
        histoChargedProtonsStatCent[i]->Sumw2();
        histoChargedProtonsStatCent[i]->Scale(0.5);
        histoChargedProtonsSystCent[i]              = (TH1D*)fileIdentfiedPubCents->Get(Form("hsys_pPb502_%s_proton_sum",centStringCharged[i].Data()));
        histoChargedProtonsSystCent[i]->Sumw2();
        histoChargedProtonsSystCent[i]->Scale(0.5);


        histoChargedKaonToPionStatCent[i]           = (TH1D*)fileIdentfiedPubCentsRatios->Get(Form("hstat_pPb502_%s_kaon_to_pion_sum",centStringCharged[i].Data()));
        histoChargedKaonToPionSystCent[i]           = (TH1D*)fileIdentfiedPubCentsRatios->Get(Form("hsys_pPb502_%s_kaon_to_pion_sum",centStringCharged[i].Data()));
        histoChargedProtonToPionStatCent[i]         = (TH1D*)fileIdentfiedPubCentsRatios->Get(Form("hstat_pPb502_%s_proton_to_pion_sum",centStringCharged[i].Data()));
        histoChargedProtonToPionSystCent[i]         = (TH1D*)fileIdentfiedPubCentsRatios->Get(Form("hsys_pPb502_%s_proton_to_pion_sum",centStringCharged[i].Data()));
    }

    if (histoChargedKaonToPionStatCent[0] && histoChargedKaonToPionStatCent[1] && histoChargedKaonToPionStatCent[2]){
        histoChargedKaonToPionStatCent[7]           = (TH1D*)histoChargedKaonToPionStatCent[0]->Clone("hstat_pPb502_0020_kaon_to_pion_sum");
        histoChargedKaonToPionStatCent[7]->Sumw2();
        histoChargedKaonToPionStatCent[7]->Scale(0.25);
        histoChargedKaonToPionStatCent[7]->Add(histoChargedKaonToPionStatCent[1], 0.25);
        histoChargedKaonToPionStatCent[7]->Add(histoChargedKaonToPionStatCent[2], 0.5);
    }
    if (histoChargedKaonToPionSystCent[0] && histoChargedKaonToPionSystCent[1] && histoChargedKaonToPionSystCent[2]){
        histoChargedKaonToPionSystCent[7]           = (TH1D*)histoChargedKaonToPionSystCent[0]->Clone("hsys_pPb502_0020_kaon_to_pion_sum");
        for (Int_t iPt = 1; iPt < histoChargedKaonToPionSystCent[7]->GetNbinsX()+1; iPt++){
            Double_t val        = histoChargedKaonToPionSystCent[0]->GetBinContent(iPt)*0.25+ histoChargedKaonToPionSystCent[1]->GetBinContent(iPt)*0.25+histoChargedKaonToPionSystCent[2]->GetBinContent(iPt)*0.5 ;
            Double_t error      = ( histoChargedKaonToPionSystCent[0]->GetBinError(iPt)/histoChargedKaonToPionSystCent[0]->GetBinContent(iPt)*0.25+
                                    histoChargedKaonToPionSystCent[1]->GetBinError(iPt)/histoChargedKaonToPionSystCent[1]->GetBinContent(iPt)*0.25+
                                    histoChargedKaonToPionSystCent[2]->GetBinError(iPt)/histoChargedKaonToPionSystCent[2]->GetBinContent(iPt)*0.5 ) * val;
            histoChargedKaonToPionSystCent[7]->SetBinContent(iPt, val);
            histoChargedKaonToPionSystCent[7]->SetBinError(iPt, error);
        }
    }
    if (histoChargedProtonToPionStatCent[0] && histoChargedProtonToPionStatCent[1] && histoChargedProtonToPionStatCent[2]){
        histoChargedProtonToPionStatCent[7]           = (TH1D*)histoChargedProtonToPionStatCent[0]->Clone("hstat_pPb502_0020_proton_to_pion_sum");
        histoChargedProtonToPionStatCent[7]->Sumw2();
        histoChargedProtonToPionStatCent[7]->Scale(0.25);
        histoChargedProtonToPionStatCent[7]->Add(histoChargedProtonToPionStatCent[1], 0.25);
        histoChargedProtonToPionStatCent[7]->Add(histoChargedProtonToPionStatCent[2], 0.5);
    }
    if (histoChargedProtonToPionSystCent[0] && histoChargedProtonToPionSystCent[1] && histoChargedProtonToPionSystCent[2]){
        histoChargedProtonToPionSystCent[7]           = (TH1D*)histoChargedProtonToPionSystCent[0]->Clone("hsys_pPb502_0020_proton_to_pion_sum");
        for (Int_t iPt = 1; iPt < histoChargedProtonToPionSystCent[7]->GetNbinsX()+1; iPt++){
            Double_t val        = histoChargedProtonToPionSystCent[0]->GetBinContent(iPt)*0.25+ histoChargedProtonToPionSystCent[1]->GetBinContent(iPt)*0.25+histoChargedProtonToPionSystCent[2]->GetBinContent(iPt)*0.5 ;
            Double_t error      = ( histoChargedProtonToPionSystCent[0]->GetBinError(iPt)/histoChargedProtonToPionSystCent[0]->GetBinContent(iPt)*0.25+
                                    histoChargedProtonToPionSystCent[1]->GetBinError(iPt)/histoChargedProtonToPionSystCent[1]->GetBinContent(iPt)*0.25+
                                    histoChargedProtonToPionSystCent[2]->GetBinError(iPt)/histoChargedProtonToPionSystCent[2]->GetBinContent(iPt)*0.5 ) * val;
            histoChargedProtonToPionSystCent[7]->SetBinContent(iPt, val);
            histoChargedProtonToPionSystCent[7]->SetBinError(iPt, error);
        }
    }
    if (histoChargedKaonToPionStatCent[5] && histoChargedKaonToPionStatCent[6]){
        histoChargedKaonToPionStatCent[8]           = (TH1D*)histoChargedKaonToPionStatCent[5]->Clone("hstat_pPb502_60100_kaon_to_pion_sum");
        histoChargedKaonToPionStatCent[8]->Sumw2();
        histoChargedKaonToPionStatCent[8]->Scale(0.5);
        histoChargedKaonToPionStatCent[8]->Add(histoChargedKaonToPionStatCent[6], 0.5);
    }
    if (histoChargedKaonToPionSystCent[5] && histoChargedKaonToPionSystCent[6]){
        histoChargedKaonToPionSystCent[8]           = (TH1D*)histoChargedKaonToPionSystCent[5]->Clone("hsys_pPb502_60100_kaon_to_pion_sum");
        for (Int_t iPt = 1; iPt < histoChargedKaonToPionSystCent[8]->GetNbinsX()+1; iPt++){
            Double_t val        = histoChargedKaonToPionSystCent[5]->GetBinContent(iPt)*0.5+histoChargedKaonToPionSystCent[6]->GetBinContent(iPt)*0.5 ;
            Double_t error      = ( histoChargedKaonToPionSystCent[5]->GetBinError(iPt)/histoChargedKaonToPionSystCent[5]->GetBinContent(iPt)*0.5+
                                    histoChargedKaonToPionSystCent[6]->GetBinError(iPt)/histoChargedKaonToPionSystCent[6]->GetBinContent(iPt)*0.5 ) * val;
            histoChargedKaonToPionSystCent[8]->SetBinContent(iPt, val);
            histoChargedKaonToPionSystCent[8]->SetBinError(iPt, error);
        }
    }
    if (histoChargedProtonToPionStatCent[5] && histoChargedProtonToPionStatCent[6]){
        histoChargedProtonToPionStatCent[8]           = (TH1D*)histoChargedProtonToPionStatCent[5]->Clone("hstat_pPb502_60100_kaon_to_pion_sum");
        histoChargedProtonToPionStatCent[8]->Sumw2();
        histoChargedProtonToPionStatCent[8]->Scale(0.5);
        histoChargedProtonToPionStatCent[8]->Add(histoChargedProtonToPionStatCent[6], 0.5);
    }
    if (histoChargedProtonToPionSystCent[5] && histoChargedProtonToPionSystCent[6]){
        histoChargedProtonToPionSystCent[8]           = (TH1D*)histoChargedProtonToPionSystCent[5]->Clone("hsys_pPb502_60100_kaon_to_pion_sum");
        for (Int_t iPt = 1; iPt < histoChargedProtonToPionSystCent[8]->GetNbinsX()+1; iPt++){
            Double_t val        = histoChargedProtonToPionSystCent[5]->GetBinContent(iPt)*0.5+histoChargedProtonToPionSystCent[6]->GetBinContent(iPt)*0.5 ;
            Double_t error      = ( histoChargedProtonToPionSystCent[5]->GetBinError(iPt)/histoChargedProtonToPionSystCent[5]->GetBinContent(iPt)*0.5+
            histoChargedProtonToPionSystCent[6]->GetBinError(iPt)/histoChargedProtonToPionSystCent[6]->GetBinContent(iPt)*0.5 ) * val;
            histoChargedProtonToPionSystCent[8]->SetBinContent(iPt, val);
            histoChargedProtonToPionSystCent[8]->SetBinError(iPt, error);
        }
    }
    for (Int_t i = 0; i< 9; i++){
        if (i < 7){
            graphChargedPionsStatCent[i]            = new TGraphAsymmErrors(histoChargedPionsStatCent[i]);
            graphChargedPionsSystCent[i]            = new TGraphAsymmErrors(histoChargedPionsSystCent[i]);
            graphChargedKaonsStatCent[i]            = new TGraphAsymmErrors(histoChargedKaonsStatCent[i]);
            graphChargedKaonsSystCent[i]            = new TGraphAsymmErrors(histoChargedKaonsSystCent[i]);
            graphChargedProtonsStatCent[i]          = new TGraphAsymmErrors(histoChargedProtonsStatCent[i]);
            graphChargedProtonsSystCent[i]          = new TGraphAsymmErrors(histoChargedProtonsSystCent[i]);
            RemoveZerosAndPrint(graphChargedPionsStatCent[i]);
            RemoveZerosAndPrint(graphChargedPionsSystCent[i]);
            RemoveZerosAndPrint(graphChargedKaonsStatCent[i]);
            RemoveZerosAndPrint(graphChargedKaonsSystCent[i]);
            RemoveZerosAndPrint(graphChargedProtonsStatCent[i]);
            RemoveZerosAndPrint(graphChargedProtonsSystCent[i]);
        }
        graphChargedKaonToPionStatCent[i]           = new TGraphAsymmErrors(histoChargedKaonToPionStatCent[i]);
        graphChargedKaonToPionStatCentW0XErr[i]     = new TGraphAsymmErrors(histoChargedKaonToPionStatCent[i]);
        graphChargedKaonToPionSystCent[i]           = new TGraphAsymmErrors(histoChargedKaonToPionSystCent[i]);
        graphChargedKaonToPionTotCent[i]            = AddErrorsOfGraphsQuadratically(graphChargedKaonToPionStatCent[i],graphChargedKaonToPionSystCent[i]);
        graphChargedProtonToPionStatCent[i]         = new TGraphAsymmErrors(histoChargedProtonToPionStatCent[i]);
        graphChargedProtonToPionStatCentW0XErr[i]   = new TGraphAsymmErrors(histoChargedProtonToPionStatCent[i]);
        graphChargedProtonToPionSystCent[i]         = new TGraphAsymmErrors(histoChargedProtonToPionSystCent[i]);
        graphChargedProtonToPionTotCent[i]          = AddErrorsOfGraphsQuadratically(graphChargedProtonToPionStatCent[i],graphChargedProtonToPionSystCent[i]);
        RemoveZerosAndPrint(graphChargedKaonToPionStatCent[i]);
        RemoveZerosAndPrint(graphChargedKaonToPionSystCent[i]);
        RemoveZerosAndPrint(graphChargedKaonToPionTotCent[i]);
        RemoveZerosAndPrint(graphChargedKaonToPionStatCentW0XErr[i]);
        RemoveZerosAndPrint(graphChargedProtonToPionSystCent[i]);
        RemoveZerosAndPrint(graphChargedProtonToPionStatCent[i]);
        RemoveZerosAndPrint(graphChargedProtonToPionTotCent[i]);
        RemoveZerosAndPrint(graphChargedProtonToPionStatCentW0XErr[i]);
        ProduceGraphAsymmWithoutXErrors(graphChargedKaonToPionStatCentW0XErr[i]);
        ProduceGraphAsymmWithoutXErrors(graphChargedProtonToPionStatCentW0XErr[i]);
        ProduceGraphAsymmWithoutXErrors(graphChargedKaonToPionTotCent[i]);
        ProduceGraphAsymmWithoutXErrors(graphChargedProtonToPionTotCent[i]);


    }

    // **************************************************************************************
    // ******************************** Reading Charged Pions *******************************
    // publication:
    // twiki: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSPECTRALowToHighPtpA
    // **************************************************************************************
    TFile* fileIdentfiedPubRpPb                     = new TFile("ExternalInputpPb/IdentifiedCharged/RpPb_502_PiKp__20151204.root");
    TH1D* histoChargedPionsRpPbStatErr              = (TH1D*)fileIdentfiedPubRpPb->Get("hstat_RpPb_pion");
    TH1D* histoChargedPionsRpPbSystErr              = (TH1D*)fileIdentfiedPubRpPb->Get("hsys_RpPb_pion");
    TH1D* histoChargedKaonsRpPbStatErr              = (TH1D*)fileIdentfiedPubRpPb->Get("hstat_RpPb_kaon");
    TH1D* histoChargedKaonsRpPbSystErr              = (TH1D*)fileIdentfiedPubRpPb->Get("hsys_RpPb_kaon");
    TH1D* histoProtonRpPbStatErr                    = (TH1D*)fileIdentfiedPubRpPb->Get("hstat_RpPb_proton");
    TH1D* histoProtonRpPbSystErr                    = (TH1D*)fileIdentfiedPubRpPb->Get("hsys_RpPb_proton");
    TH1D* histoChargedHadronRpPbStatErr             = (TH1D*)fileIdentfiedPubRpPb->Get("hstat_RpPb_charged"); // /HepData/8550/d4x1y1
    TH1D* histoChargedHadronRpPbSystErr             = (TH1D*)fileIdentfiedPubRpPb->Get("hsys_RpPb_charged"); // /HepData/8550/d4x1y1
    TGraphAsymmErrors* graphChargedPionsRpPbStatErr = new TGraphAsymmErrors(histoChargedPionsRpPbStatErr);
    TGraphAsymmErrors* graphChargedPionsRpPbSystErr = new TGraphAsymmErrors(histoChargedPionsRpPbSystErr);
    TGraphAsymmErrors* graphChargedKaonsRpPbStatErr = new TGraphAsymmErrors(histoChargedKaonsRpPbStatErr);
    TGraphAsymmErrors* graphChargedKaonsRpPbSystErr = new TGraphAsymmErrors(histoChargedKaonsRpPbSystErr);
    TGraphAsymmErrors* graphProtonRpPbStatErr       = new TGraphAsymmErrors(histoProtonRpPbStatErr);
    TGraphAsymmErrors* graphProtonRpPbSystErr       = new TGraphAsymmErrors(histoProtonRpPbSystErr);
    TGraphAsymmErrors* graphChargedHadronRpPbSystErr= new TGraphAsymmErrors(histoChargedHadronRpPbSystErr);
    TGraphAsymmErrors* graphChargedHadronRpPbStatErr= new TGraphAsymmErrors(histoChargedHadronRpPbStatErr);

    RemoveZerosAndPrint(graphChargedPionsRpPbStatErr);
    RemoveZerosAndPrint(graphChargedPionsRpPbSystErr);
    RemoveZerosAndPrint(graphChargedKaonsRpPbStatErr);
    RemoveZerosAndPrint(graphChargedKaonsRpPbSystErr);
    RemoveZerosAndPrint(graphProtonRpPbStatErr);
    RemoveZerosAndPrint(graphProtonRpPbSystErr);
    RemoveZerosAndPrint(graphChargedHadronRpPbStatErr);
    RemoveZerosAndPrint(graphChargedHadronRpPbSystErr);

    TGraphAsymmErrors* graphChargedHadronRpPbStatErrW0XErr     = (TGraphAsymmErrors*)graphChargedHadronRpPbStatErr->Clone("graphChWOXErr");
    ProduceGraphAsymmWithoutXErrors(graphChargedHadronRpPbStatErrW0XErr);
    TGraphAsymmErrors* graphChargedHadronRpPbTotErr          = AddErrorsOfGraphsQuadratically(graphChargedHadronRpPbStatErr,graphChargedHadronRpPbSystErr);
    ProduceGraphAsymmWithoutXErrors(graphChargedHadronRpPbTotErr);


    TGraphAsymmErrors* graphD0QpPbStat[3][4];
    TGraphAsymmErrors* graphD0QpPbStatW0XErr[3][4];
    TGraphAsymmErrors* graphD0QpPbSys[3][4];
    TGraphAsymmErrors* graphD0QpPbTot[3][4];
    TString centNameWide[4]     = {"0020", "2040", "4060", "60100"};
    TString centNameWideOut[4]  = {"0-20%", "20-40%", "40-60%", "60-100%"};
    TString centNameEst[3]      = {"V0A", "CL1", "ZNA"};
    Double_t nCollWide[3][4]    = {{0, 0, 0, 0},{0, 0, 0, 0}, {0, 0, 0, 0}} ;

    for (Int_t est = 0; est < 3; est++){
        for (Int_t cent = 0; cent < 4; cent++){
            TString DmesonHEPDataFile        = Form("ExternalInputpPb/IdentifiedCharged/Dmeson_QpPb_%s_%s.csv", centNameEst[est].Data(), centNameWide[cent].Data()) ;
            graphD0QpPbStat[est][cent]            = ParseHEPData(DmesonHEPDataFile, 8, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
            graphD0QpPbSys[est][cent]             = ParseHEPData(DmesonHEPDataFile, 8, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE, kFALSE);
            nCollWide[est][cent]                  = GetNCollFromName(centNameWide[cent]+"_"+centNameEst[est], "pPb_5.023TeV");
            graphD0QpPbStatW0XErr[est][cent]      = (TGraphAsymmErrors*)graphD0QpPbStat[est][cent]->Clone("graphDWOXErr");
            ProduceGraphAsymmWithoutXErrors(graphD0QpPbStatW0XErr[est][cent]);
            graphD0QpPbTot[est][cent]             = AddErrorsOfGraphsQuadratically(graphD0QpPbStat[est][cent],graphD0QpPbSys[est][cent]);
            ProduceGraphAsymmWithoutXErrors(graphD0QpPbTot[est][cent]);
        }
    }

    TGraphAsymmErrors* graphChHadQpPbStat[3][9];
    TGraphAsymmErrors* graphChHadQpPbStatW0XErr[3][9];
    TGraphAsymmErrors* graphChHadQpPbSys[3][9];
    TGraphAsymmErrors* graphChHadQpPbTot[3][9];
    TString centNameFine[9]     = {"0005", "0510", "1020", "2040", "4060", "6080", "80100", "0020", "60100"};
    TString centNameFineOut[9]  = {"0-5%", "5-10%", "10-20%", "20-40%", "40-60%", "60-80%", "80-100%", "0-20%", "60-100%"};
    Double_t nCollFine[3][7]    = {{0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}};

    for (Int_t est = 0; est < 3; est++){
        for (Int_t cent = 0; cent < 7; cent++){
            TString chargedHadHEPFile           = Form("ExternalInputpPb/IdentifiedCharged/ChargedQpPb_%s_%s.csv", centNameEst[est].Data(), centNameFine[cent].Data()) ;
            graphChHadQpPbStat[est][cent]           = ParseHEPData(chargedHadHEPFile, 12, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
            graphChHadQpPbSys[est][cent]            = ParseHEPData(chargedHadHEPFile, 12, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE, kFALSE);
            nCollFine[est][cent]                    = GetNCollFromName(centNameFine[cent]+"_"+centNameEst[est], "pPb_5.023TeV");
            graphChHadQpPbStatW0XErr[est][cent]     = (TGraphAsymmErrors*)graphChHadQpPbStat[est][cent]->Clone("graphChWOXErr");
            ProduceGraphAsymmWithoutXErrors(graphChHadQpPbStatW0XErr[est][cent]);
            graphChHadQpPbTot[est][cent]            = AddErrorsOfGraphsQuadratically(graphChHadQpPbStat[est][cent],graphChHadQpPbSys[est][cent]);
            ProduceGraphAsymmWithoutXErrors(graphChHadQpPbTot[est][cent]);
        }
        if ( graphChHadQpPbStat[est][0] && graphChHadQpPbStat[est][1] && graphChHadQpPbStat[est][2]){
            graphChHadQpPbStat[est][7]          = (TGraphAsymmErrors*)graphChHadQpPbStat[est][0]->Clone("cent0020Stat");
            for (Int_t i = 0; i < graphChHadQpPbStat[est][7]->GetN(); i++){
                Double_t yVal                   = graphChHadQpPbStat[est][0]->GetY()[i]*0.25*nCollFine[est][0]/nCollWide[est][0]+
                                                  graphChHadQpPbStat[est][1]->GetY()[i]*0.25*nCollFine[est][1]/nCollWide[est][0]+
                                                  graphChHadQpPbStat[est][1]->GetY()[i]*0.5*nCollFine[est][2]/nCollWide[est][0];
                Double_t yErrLow                = graphChHadQpPbStat[est][0]->GetEYlow()[i]*0.25*nCollFine[est][0]/nCollWide[est][0]+
                                                  graphChHadQpPbStat[est][1]->GetEYlow()[i]*0.25*nCollFine[est][1]/nCollWide[est][0]+
                                                  graphChHadQpPbStat[est][1]->GetEYlow()[i]*0.5*nCollFine[est][2]/nCollWide[est][0];
                Double_t yErrHigh               = graphChHadQpPbStat[est][0]->GetEYhigh()[i]*0.25*nCollFine[est][0]/nCollWide[est][0]+
                                                  graphChHadQpPbStat[est][1]->GetEYhigh()[i]*0.25*nCollFine[est][1]/nCollWide[est][0]+
                                                  graphChHadQpPbStat[est][1]->GetEYhigh()[i]*0.5*nCollFine[est][2]/nCollWide[est][0];
                graphChHadQpPbStat[est][7]->SetPoint(i, graphChHadQpPbStat[est][7]->GetX()[i], yVal);
                graphChHadQpPbStat[est][7]->SetPointError(i, graphChHadQpPbStat[est][7]->GetEXlow()[i], graphChHadQpPbStat[est][7]->GetEXlow()[i], yErrLow, yErrHigh);
            }
            graphChHadQpPbStatW0XErr[est][7]     = (TGraphAsymmErrors*)graphChHadQpPbStat[est][7]->Clone("graphChWOXErr");
            ProduceGraphAsymmWithoutXErrors(graphChHadQpPbStatW0XErr[est][7]);
        }
        if ( graphChHadQpPbSys[est][0] && graphChHadQpPbSys[est][1] && graphChHadQpPbSys[est][2]){
            graphChHadQpPbSys[est][7]          = (TGraphAsymmErrors*)graphChHadQpPbSys[est][0]->Clone("cent0020Sys");
            for (Int_t i = 0; i < graphChHadQpPbSys[est][7]->GetN(); i++){
                Double_t yVal                   = graphChHadQpPbSys[est][0]->GetY()[i]*0.25*nCollFine[est][0]/nCollWide[est][0]+
                graphChHadQpPbSys[est][1]->GetY()[i]*0.25*nCollFine[est][1]/nCollWide[est][0]+
                graphChHadQpPbSys[est][1]->GetY()[i]*0.5*nCollFine[est][2]/nCollWide[est][0];
                Double_t yErrLow                = graphChHadQpPbSys[est][0]->GetEYlow()[i]*0.25*nCollFine[est][0]/nCollWide[est][0]+
                graphChHadQpPbSys[est][1]->GetEYlow()[i]*0.25*nCollFine[est][1]/nCollWide[est][0]+
                graphChHadQpPbSys[est][1]->GetEYlow()[i]*0.5*nCollFine[est][2]/nCollWide[est][0];
                Double_t yErrHigh               = graphChHadQpPbSys[est][0]->GetEYhigh()[i]*0.25*nCollFine[est][0]/nCollWide[est][0]+
                graphChHadQpPbSys[est][1]->GetEYhigh()[i]*0.25*nCollFine[est][1]/nCollWide[est][0]+
                graphChHadQpPbSys[est][1]->GetEYhigh()[i]*0.5*nCollFine[est][2]/nCollWide[est][0];
                graphChHadQpPbSys[est][7]->SetPoint(i, graphChHadQpPbSys[est][7]->GetX()[i], yVal);
                graphChHadQpPbSys[est][7]->SetPointError(i, graphChHadQpPbSys[est][7]->GetEXlow()[i], graphChHadQpPbSys[est][7]->GetEXlow()[i], yErrLow, yErrHigh);
            }
            graphChHadQpPbTot[est][7]            = AddErrorsOfGraphsQuadratically(graphChHadQpPbStat[est][7],graphChHadQpPbSys[est][7]);
            ProduceGraphAsymmWithoutXErrors(graphChHadQpPbTot[est][7]);
        }
        if ( graphChHadQpPbStat[est][5] && graphChHadQpPbStat[est][6] ){
            graphChHadQpPbStat[est][8]          = (TGraphAsymmErrors*)graphChHadQpPbStat[est][0]->Clone("cent60100Stat");
            for (Int_t i = 0; i < graphChHadQpPbStat[est][7]->GetN(); i++){
                Double_t yVal                   = graphChHadQpPbStat[est][5]->GetY()[i]*0.5*nCollFine[est][5]/nCollWide[est][3]+
                                                  graphChHadQpPbStat[est][6]->GetY()[i]*0.5*nCollFine[est][6]/nCollWide[est][3];
                Double_t yErrLow                = graphChHadQpPbStat[est][5]->GetEYlow()[i]*0.5*nCollFine[est][5]/nCollWide[est][3]+
                                                  graphChHadQpPbStat[est][6]->GetEYlow()[i]*0.5*nCollFine[est][6] /nCollWide[est][3];
                Double_t yErrHigh               = graphChHadQpPbStat[est][5]->GetEYhigh()[i]*0.5*nCollFine[est][5]/nCollWide[est][3]+
                                                  graphChHadQpPbStat[est][6]->GetEYhigh()[i]*0.5*nCollFine[est][6]/nCollWide[est][3];
                graphChHadQpPbStat[est][8]->SetPoint(i, graphChHadQpPbStat[est][8]->GetX()[i], yVal);
                graphChHadQpPbStat[est][8]->SetPointError(i, graphChHadQpPbStat[est][8]->GetEXlow()[i], graphChHadQpPbStat[est][7]->GetEXlow()[i], yErrLow, yErrHigh);
            }
            graphChHadQpPbStatW0XErr[est][8]     = (TGraphAsymmErrors*)graphChHadQpPbStat[est][8]->Clone("graphChWOXErr");
            ProduceGraphAsymmWithoutXErrors(graphChHadQpPbStatW0XErr[est][8]);
        }
        if ( graphChHadQpPbSys[est][5] && graphChHadQpPbSys[est][6] ){
            graphChHadQpPbSys[est][8]          = (TGraphAsymmErrors*)graphChHadQpPbSys[est][0]->Clone("cent60100Sys");
            for (Int_t i = 0; i < graphChHadQpPbSys[est][7]->GetN(); i++){
                Double_t yVal                   = graphChHadQpPbSys[est][5]->GetY()[i]*0.5*nCollFine[est][5]/nCollWide[est][3]+
                graphChHadQpPbSys[est][6]->GetY()[i]*0.5*nCollFine[est][6]/nCollWide[est][3];
                Double_t yErrLow                = graphChHadQpPbSys[est][5]->GetEYlow()[i]*0.5*nCollFine[est][5]/nCollWide[est][3]+
                graphChHadQpPbSys[est][6]->GetEYlow()[i]*0.5*nCollFine[est][6] /nCollWide[est][3];
                Double_t yErrHigh               = graphChHadQpPbSys[est][5]->GetEYhigh()[i]*0.5*nCollFine[est][5]/nCollWide[est][3]+
                graphChHadQpPbSys[est][6]->GetEYhigh()[i]*0.5*nCollFine[est][6]/nCollWide[est][3];
                graphChHadQpPbSys[est][8]->SetPoint(i, graphChHadQpPbSys[est][8]->GetX()[i], yVal);
                graphChHadQpPbSys[est][8]->SetPointError(i, graphChHadQpPbSys[est][8]->GetEXlow()[i], graphChHadQpPbSys[est][7]->GetEXlow()[i], yErrLow, yErrHigh);
                graphChHadQpPbTot[est][8]            = AddErrorsOfGraphsQuadratically(graphChHadQpPbStat[est][8],graphChHadQpPbSys[est][8]);
                ProduceGraphAsymmWithoutXErrors(graphChHadQpPbTot[est][8]);
            }
        }
    }




    // *********************************************************************************************************************
    // ********************************** Write Output files ***************************************************************
    // *********************************************************************************************************************
    TFile fileChargedPionspPb(Form("ExternalInputpPb/IdentifiedCharged/ChargedIdentifiedSpectrapPb_%s.root",dateForOutput.Data()) ,"RECREATE");
        histoChargedPionSpecLowPtSyspPb->Write("histoChargedPionSpecLowPtSyspPb", TObject::kOverwrite);
        histoChargedPionSpecLowPtStatpPb->Write("histoChargedPionSpecLowPtStatpPb", TObject::kOverwrite);
        histoChargedPionSpecLowPtSyspPb0020->Write("histoChargedPionSpecLowPtSyspPb0020", TObject::kOverwrite);
        histoChargedPionSpecLowPtStatpPb0020->Write("histoChargedPionSpecLowPtStatpPb0020", TObject::kOverwrite);
        histoChargedPionSpecLowPtSyspPb2040->Write("histoChargedPionSpecLowPtSyspPb2040", TObject::kOverwrite);
        histoChargedPionSpecLowPtStatpPb2040->Write("histoChargedPionSpecLowPtStatpPb2040", TObject::kOverwrite);
        histoChargedPionSpecLowPtSyspPb4060->Write("histoChargedPionSpecLowPtSyspPb4060", TObject::kOverwrite);
        histoChargedPionSpecLowPtStatpPb4060->Write("histoChargedPionSpecLowPtStatpPb4060", TObject::kOverwrite);
        histoChargedPionSpecLowPtSyspPb6080->Write("histoChargedPionSpecLowPtSyspPb6080", TObject::kOverwrite);
        histoChargedPionSpecLowPtStatpPb6080->Write("histoChargedPionSpecLowPtStatpPb6080", TObject::kOverwrite);
        histoChargedPionSpecLowPtSyspPb80100->Write("histoChargedPionSpecLowPtSyspPb80100", TObject::kOverwrite);
        histoChargedPionSpecLowPtStatpPb80100->Write("histoChargedPionSpecLowPtStatpPb80100", TObject::kOverwrite);

        TDirectoryFile* centFolder          = (TDirectoryFile*)fileChargedPionspPb.Get("pPb5TeV_0-100%");

        if (!centFolder){
            fileChargedPionspPb.mkdir("pPb5TeV_0-100%");
            centFolder              = (TDirectoryFile*)fileChargedPionspPb.Get("pPb5TeV_0-100%");
        }
        centFolder->cd();


        if (histoChargedPionspPbStatErr) histoChargedPionspPbStatErr->Write("histoChargedPionPubStat", TObject::kOverwrite);
        if (histoChargedPionspPbSystErr) histoChargedPionspPbSystErr->Write("histoChargedPionPubSys", TObject::kOverwrite);
        if (histoChargedKaonspPbStatErr) histoChargedKaonspPbStatErr->Write("histoChargedKaonPubStat", TObject::kOverwrite);
        if (histoChargedKaonspPbSystErr) histoChargedKaonspPbSystErr->Write("histoChargedKaonPubSys", TObject::kOverwrite);
        if (histoChargedProtonspPbStatErr) histoChargedProtonspPbStatErr->Write("histoProtonPubStat", TObject::kOverwrite);
        if (histoChargedProtonspPbSystErr) histoChargedProtonspPbSystErr->Write("histoProtonPubSys", TObject::kOverwrite);
        if (graphChargedPionspPbSyst) graphChargedPionspPbSyst->Write("graphChargedPionPubSys", TObject::kOverwrite);
        if (graphChargedPionspPbStat) graphChargedPionspPbStat->Write("graphChargedPionPubStat", TObject::kOverwrite);
        if (graphChargedKaonspPbStat) graphChargedKaonspPbStat->Write("graphChargedKaonPubStat", TObject::kOverwrite);
        if (graphChargedKaonspPbSyst) graphChargedKaonspPbSyst->Write("graphChargedKaonPubSys", TObject::kOverwrite);
        if (graphChargedProtonspPbStat) graphChargedProtonspPbStat->Write("graphProtonPubStat", TObject::kOverwrite);
        if (graphChargedProtonspPbSyst) graphChargedProtonspPbSyst->Write("graphProtonPubSys", TObject::kOverwrite);

        if (histoChargedPionsRpPbStatErr) histoChargedPionsRpPbStatErr->Write("histoChargedPionPubStat_RpPb");
        if (histoChargedPionsRpPbSystErr) histoChargedPionsRpPbSystErr->Write("histoChargedPionPubSys_RpPb");
        if (histoChargedKaonsRpPbStatErr) histoChargedKaonsRpPbStatErr->Write("histoChargedKaonPubStat_RpPb");
        if (histoChargedKaonsRpPbSystErr) histoChargedKaonsRpPbSystErr->Write("histoChargedKaonPubSys_RpPb");
        if (histoProtonRpPbStatErr) histoProtonRpPbStatErr->Write("histoProtonPubStat_RpPb");
        if (histoProtonRpPbSystErr) histoProtonRpPbSystErr->Write("histoProtonPubSys_RpPb");
        if (histoChargedHadronRpPbStatErr) histoChargedHadronRpPbStatErr->Write("histoChargedHadronPubStat_RpPb");
        if (histoChargedHadronRpPbSystErr) histoChargedHadronRpPbSystErr->Write("histoChargedHadronPubSys_RpPb");
        if (graphChargedPionsRpPbStatErr) graphChargedPionsRpPbStatErr->Write("graphChargedPionPubStat_RpPb");
        if (graphChargedPionsRpPbSystErr) graphChargedPionsRpPbSystErr->Write("graphChargedPionPubSys_RpPb");
        if (graphChargedKaonsRpPbStatErr) graphChargedKaonsRpPbStatErr->Write("graphChargedKaonPubStat_RpPb");
        if (graphChargedKaonsRpPbSystErr) graphChargedKaonsRpPbSystErr->Write("graphChargedKaonPubSys_RpPb");
        if (graphProtonRpPbStatErr) graphProtonRpPbStatErr->Write("graphProtonPubStat_RpPb");
        if (graphProtonRpPbSystErr) graphProtonRpPbSystErr->Write("graphProtonPubSys_RpPb");
        if (graphChargedHadronRpPbStatErr) graphChargedHadronRpPbStatErr->Write("ChargedHadQpPbStat");
        if (graphChargedHadronRpPbSystErr) graphChargedHadronRpPbSystErr->Write("ChargedHadQpPbSys");
        if (graphChargedHadronRpPbStatErrW0XErr) graphChargedHadronRpPbStatErrW0XErr->Write(Form("ChargedHadQpPbStatW0XErr"),TObject::kOverwrite );
        if (graphChargedHadronRpPbTotErr) graphChargedHadronRpPbTotErr->Write(Form("ChargedHadQpPbTot"),TObject::kOverwrite );

        for (Int_t cent = 0; cent < 9; cent++){
            for (Int_t est = 0; est < 3; est++){
                if (cent < 9){
                    cout << cent << "\t" << Form("pPb5TeV_%s_%s", centNameFineOut[cent].Data(), centNameEst[est].Data()) << endl;
                    TDirectoryFile* centFolder          = (TDirectoryFile*)fileChargedPionspPb.Get(Form("pPb5TeV_%s_%s", centNameFineOut[cent].Data(), centNameEst[est].Data()));

                    if (!centFolder){
                        fileChargedPionspPb.mkdir(Form("pPb5TeV_%s_%s",  centNameFineOut[cent].Data(), centNameEst[est].Data()));
                        centFolder              = (TDirectoryFile*)fileChargedPionspPb.Get(Form("pPb5TeV_%s_%s",  centNameFineOut[cent].Data(), centNameEst[est].Data()));
                    }
                    centFolder->cd();
                    if (est == 0 && cent < 7){
                        if (histoChargedPionsStatCent[cent]) histoChargedPionsStatCent[cent]->Write(Form("histoChargedPionPubStat"), TObject::kOverwrite);
                        if (histoChargedPionsSystCent[cent]) histoChargedPionsStatCent[cent]->Write(Form("histoChargedPionPubSys"), TObject::kOverwrite);
                        if (histoChargedKaonsStatCent[cent]) histoChargedKaonsStatCent[cent]->Write(Form("histoChargedKaonPubStat"), TObject::kOverwrite);
                        if (histoChargedKaonsSystCent[cent]) histoChargedKaonsStatCent[cent]->Write(Form("histoChargedKaonPubSys"), TObject::kOverwrite);
                        if (histoChargedProtonsStatCent[cent]) histoChargedProtonsStatCent[cent]->Write(Form("histoProtonPubStat"), TObject::kOverwrite);
                        if (histoChargedProtonsSystCent[cent]) histoChargedProtonsStatCent[cent]->Write(Form("histoProtonPubSys"), TObject::kOverwrite);
                        if (graphChargedPionsStatCent[cent]) graphChargedPionsStatCent[cent]->Write(Form("graphChargedPionPubStat"), TObject::kOverwrite);
                        if (graphChargedPionsSystCent[cent]) graphChargedPionsStatCent[cent]->Write(Form("graphChargedPionPubSys"), TObject::kOverwrite);
                        if (graphChargedKaonsStatCent[cent]) graphChargedKaonsStatCent[cent]->Write(Form("graphChargedKaonPubStat"), TObject::kOverwrite);
                        if (graphChargedKaonsSystCent[cent]) graphChargedKaonsStatCent[cent]->Write(Form("graphChargedKaonPubSys"), TObject::kOverwrite);
                        if (graphChargedProtonsStatCent[cent]) graphChargedProtonsStatCent[cent]->Write(Form("graphProtonPubStat"), TObject::kOverwrite);
                        if (graphChargedProtonsSystCent[cent]) graphChargedProtonsStatCent[cent]->Write(Form("graphProtonPubSys"), TObject::kOverwrite);
                    }
                    if (est == 0){
                        if (graphChargedKaonToPionStatCent[cent]) graphChargedKaonToPionStatCent[cent]->Write(Form("graphChargedKaonToPionPubStat"), TObject::kOverwrite);
                        if (graphChargedKaonToPionStatCentW0XErr[cent]) graphChargedKaonToPionStatCentW0XErr[cent]->Write(Form("graphChargedKaonToPionPubStatW0XErr"), TObject::kOverwrite);
                        if (graphChargedKaonToPionSystCent[cent]) graphChargedKaonToPionSystCent[cent]->Write(Form("graphChargedKaonToPionPubSys"), TObject::kOverwrite);
                        if (graphChargedKaonToPionTotCent[cent]) graphChargedKaonToPionTotCent[cent]->Write(Form("graphChargedKaonToPionPubTot"), TObject::kOverwrite);
                        if (graphChargedProtonToPionStatCent[cent]) graphChargedProtonToPionStatCent[cent]->Write(Form("graphProtonToPionPubStat"), TObject::kOverwrite);
                        if (graphChargedProtonToPionStatCentW0XErr[cent]) graphChargedProtonToPionStatCentW0XErr[cent]->Write(Form("graphProtonToPionPubStatW0XErr"), TObject::kOverwrite);
                        if (graphChargedProtonToPionSystCent[cent]) graphChargedProtonToPionSystCent[cent]->Write(Form("graphProtonToPionPubSys"), TObject::kOverwrite);
                        if (graphChargedProtonToPionTotCent[cent]) graphChargedProtonToPionTotCent[cent]->Write(Form("graphProtonToPionPubTot"), TObject::kOverwrite);
                    }

                    if (graphChHadQpPbStat[est][cent]) graphChHadQpPbStat[est][cent]->Write(Form("ChargedHadQpPbStat"),TObject::kOverwrite );
                    if (graphChHadQpPbStatW0XErr[est][cent]) graphChHadQpPbStatW0XErr[est][cent]->Write(Form("ChargedHadQpPbStatW0XErr"),TObject::kOverwrite );
                    if (graphChHadQpPbSys[est][cent]) graphChHadQpPbSys[est][cent]->Write(Form("ChargedHadQpPbSys"),TObject::kOverwrite );
                    if (graphChHadQpPbTot[est][cent]) graphChHadQpPbTot[est][cent]->Write(Form("ChargedHadQpPbTot"),TObject::kOverwrite );
                }
            }
        }
        for (Int_t est = 0; est < 3; est++){
            for (Int_t cent = 0; cent < 4; cent++){
                TDirectoryFile* centFolder          = (TDirectoryFile*)fileChargedPionspPb.Get(Form("pPb5TeV_%s_%s", centNameWideOut[cent].Data(), centNameEst[est].Data()));

                if (!centFolder){
                    fileChargedPionspPb.mkdir(Form("pPb5TeV_%s_%s",  centNameWideOut[cent].Data(), centNameEst[est].Data()));
                    centFolder              = (TDirectoryFile*)fileChargedPionspPb.Get(Form("pPb5TeV_%s_%s",  centNameWideOut[cent].Data(), centNameEst[est].Data()));

                }
                centFolder->cd();
                if (graphD0QpPbStat[est][cent]) graphD0QpPbStat[est][cent]->Write(Form("DMesonQpPbStat"),TObject::kOverwrite );
                if (graphD0QpPbStatW0XErr[est][cent]) graphD0QpPbStatW0XErr[est][cent]->Write(Form("DMesonQpPbStatW0XErr"),TObject::kOverwrite );
                if (graphD0QpPbSys[est][cent]) graphD0QpPbSys[est][cent]->Write(Form("DMesonQpPbSys"),TObject::kOverwrite );
                if (graphD0QpPbTot[est][cent]) graphD0QpPbTot[est][cent]->Write(Form("DMesonQpPbTot"),TObject::kOverwrite );

            }
        }
    fileChargedPionspPb.Close();

}
