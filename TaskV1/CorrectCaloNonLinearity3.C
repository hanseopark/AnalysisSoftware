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
#include "TProfile2D.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TMath.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/FittingGammaConversion.h"

TF1* FitRecursiveGaussian (TH1* histo, Double_t precision, Double_t correctRange, Double_t fitRangeMin, Double_t fitRangeMax);
TF1* FitExpPlusGaussian(TH1D* histo, Double_t fitRangeMin, Double_t fitRangeMax, Int_t mode);
TF1* FitBckg(TH1* fHisto, Double_t minFit, Double_t maxFit);
TF1* FitDataMC(TH1* fHisto, Double_t minFit, Double_t maxFit, TString selection, Double_t constPar = -1);
Float_t    FunctionNL_kSDM(Float_t e, Float_t p0, Float_t p1, Float_t p2);
template<class ForwardIt>

void SetLogBinningXTH(ForwardIt* histoRebin){
    TAxis *axisafter    = histoRebin->GetXaxis();
    Int_t bins          = axisafter->GetNbins();
    Double_t from       = axisafter->GetXmin();
    Double_t to         = axisafter->GetXmax();
    Double_t *newbins   = new Double_t[bins+1];
    newbins[0]          = from;
    Double_t factor     = TMath::Power(to/from, 1./bins);
    for(Int_t i=1; i<=bins; ++i) 
        newbins[i]      = factor * newbins[i-1];
    axisafter->Set(bins, newbins);
    delete [] newbins;
    return;
}

//****************************************************************************
//************** Function to Correct CaloNonLinearity3 ***********************
//****************************************************************************
void CorrectCaloNonLinearity3(TString select = "LHC11a-Pythia-ConvCalo")
{
    gROOT->Reset();

    StyleSettingsThesis();
    SetPlotStyle();

    TString outputDir           = "CorrectCaloNonLinearity3";
    gSystem->Exec("mkdir -p "+outputDir);
    TString outputDirSample     = Form("CorrectCaloNonLinearity3/%s",select.Data());
    gSystem->Exec("mkdir -p "+outputDirSample);

    // General options
    TString suffix              = "eps";
    TString optionEnergy        = "";
    TString fPlot[6]            = {"", "", "", "", "", ""};
    TString strDataFile[6]      = {"", "", "", "", "", ""};
    TString strMCFile[6]        = {"", "", "", "", "", ""};
    TString dataCut[6]          = {"", "", "", "", "", ""};
    TString mcCut[6]            = {"", "", "", "", "", ""};
    TString dataMainDir[6]      = {"", "", "", "", "", ""};
    TString mcMainDir[6]        = {"", "", "", "", "", ""};
    Int_t mode                  = 2;

    const Int_t nColor          = 13;
    const Int_t nStyle          = 7;
    Color_t color[nColor]       = {kBlack,633,807,/*800,*/418,/*kGreen+4,*/435,601,879,806,852,kCyan+3,426};
    Int_t markerStyle[nStyle]   = {24,25,27,28,29,30,31};
    Double_t massPi0            = 0.1349766;

    // pT range for mass fitting
    Int_t startPtBin            = 0;
    Int_t endPtBin              = 17;
    Int_t firstTriggerBin[5]    = {-1, -1, -1, -1, -1};
    Int_t exampleBin1           = -1;
    Int_t exampleBin2           = -1;
    Int_t fixedOffSet           = -1;
    //*******************************************************************************
    // Choosing data set
    if (select.Contains("LHC11a") || select.Contains("LHC13g")){
        optionEnergy        = "2.76TeV";        
//         fixedOffSet         = 0.977024661;
    } else  if (select.Contains("LHC10") ) {
        optionEnergy        = "7TeV";        
    } else  if (select.Contains("LHC12") ) {
        optionEnergy        = "8TeV";
    } else  if (select.Contains("LHC13bc") ) {
        optionEnergy        = "pPb_5.023TeV";        
    }    
        
    if (select.Contains("ConvCalo") ){
        mode                = 2;
        for (Int_t i = 0; i < 6; i++){
            dataMainDir[i]  = "GammaConvCalo";
            mcMainDir[i]    = "GammaConvCalo";
        }    
    } else if (select.Contains("Calo")){
        mode                = 4;
        for (Int_t i = 0; i < 6; i++){
            dataMainDir[i]  = "GammaCalo";
            mcMainDir[i]    = "GammaCalo";
        }    
    }

    if (select.Contains("ConvCalo") && (select.Contains("LHC11a") || select.Contains("LHC13g")) ){
        startPtBin          = 0;
        endPtBin            = 21;        
        exampleBin1         = 4;
        exampleBin2         = 19;
    } else if (select.Contains("Calo") && (select.Contains("LHC11a") || select.Contains("LHC13g")) ){
        startPtBin          = 0;
        endPtBin            = 17;       
        exampleBin1         = 4;
        exampleBin2         = 16;
    } else if (select.Contains("ConvCalo") && select.Contains("LHC10") ){
        startPtBin          = 0;
        endPtBin            = 20;
        exampleBin1         = 4;
        exampleBin2         = 19;
    } else if (select.Contains("Calo") && select.Contains("LHC10") ){
        startPtBin          = 0;
        endPtBin            = 17;       
        exampleBin1         = 4;
        exampleBin2         = 16;
    } else if (select.Contains("ConvCalo") && select.Contains("LHC12") ){
        startPtBin          = 3;
        endPtBin            = 20;
        exampleBin1         = 4;
        exampleBin2         = 15;
    } else if (select.Contains("Calo") && select.Contains("LHC12") ){
        startPtBin          = 3;
        endPtBin            = 16;
        if(select.Contains("Pythia")) endPtBin=15;
        exampleBin1         = 8;
        exampleBin2         = 14;
    } else if (select.Contains("ConvCalo") && select.Contains("LHC13bc") ){
        startPtBin          = 0;
        endPtBin            = 20;
        exampleBin1         = 4;
        exampleBin2         = 19;
    } else if (select.Contains("Calo") && select.Contains("LHC13bc") ){
        startPtBin          = 0;
        endPtBin            = 17;       
        exampleBin1         = 4;
        exampleBin2         = 16;        
    }    
    

    // Standard configuration for LHC11a with default timing cuts in tender
    if(select.CompareTo("LHC11a-Pythia-Calo")==0){
        strDataFile[0]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20151213-CalibDefaultTiming/GammaCalo_LHC11a-pass4_40.root";
        strMCFile[0]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20151213-CalibDefaultTiming/GammaCalo_MC_LHC12f1a_LHC12i3_40.root";
        dataCut[0]          = "00003113_1111100053032220000_0163103100000050";
        mcCut[0]            = "00003113_1111100053032220000_0163103100000050";

        firstTriggerBin[0]  = 12;
        strDataFile[1]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20151213-CalibDefaultTiming/GammaCalo_LHC11a-pass4_40.root";
        strMCFile[1]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20151213-CalibDefaultTiming/GammaCalo_MC_LHC15g1aFinerPtHardBins_40.root";
        dataCut[1]          = "00051113_1111100053032220000_0163103100000050";
        mcCut[1]            = "00051113_1111100053032220000_0163103100000050";

        fPlot[0]            = "#frac{LHC12f1a & LHC12i3}{LHC11a - INT1}";
        fPlot[1]            = "#frac{LHC15g1a}{LHC11a - EMC1}";
    } else if(select.CompareTo("LHC11a-Pythia-ConvCalo")==0){
        strDataFile[0]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20151213-CalibDefaultTiming/GammaConvCalo_LHC11a-pass4_28.root";
        strMCFile[0]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20151213-CalibDefaultTiming/GammaConvCalo_MC_LHC12f1a_28.root";
        dataCut[0]          = "00003113_00200009327000008250400000_1111100053032230000_0163103100000010";
        mcCut[0]            = "00003113_00200009327000008250400000_1111100053032230000_0163103100000010";

        firstTriggerBin[0]  = 13;
        strDataFile[1]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20151213-CalibDefaultTiming/GammaConvCalo_LHC11a-pass4_28.root";
        strMCFile[1]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20151213-CalibDefaultTiming/GammaConvCalo_MC_LHC15g1aFinerPtHardBins_34.root";
        dataCut[1]          = "00051113_00200009327000008250400000_1111100053032230000_0163103100000010";
        mcCut[1]            = "00051113_00200009327000008250400000_1111100053032230000_0163103100000010";

        fPlot[0]            = "#frac{LHC12f1a}{LHC11a - INT1}";
        fPlot[1]            = "#frac{LHC15g1a}{LHC11a - EMC1}";
    } else if(select.CompareTo("LHC11a-Phojet-Calo")==0){
        strDataFile[0]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20151213-CalibDefaultTiming/GammaCalo_LHC11a-pass4_40.root";
        strMCFile[0]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20151213-CalibDefaultTiming/GammaCalo_MC_LHC12f1b_40.root";
        dataCut[0]          = "00003113_1111100053032220000_0163103100000050";
        mcCut[0]            = "00003113_1111100053032220000_0163103100000050";

        firstTriggerBin[0]  = 12;
        strDataFile[1]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20151213-CalibDefaultTiming/GammaCalo_LHC11a-pass4_40.root";
        strMCFile[1]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20151213-CalibDefaultTiming/GammaCalo_MC_LHC15g1aFinerPtHardBins_40.root";
        dataCut[1]          = "00051113_1111100053032220000_0163103100000050";
        mcCut[1]            = "00051113_1111100053032220000_0163103100000050";

        fPlot[0]            = "#frac{LHC12f1b}{LHC11a - INT1}";
        fPlot[1]            = "#frac{LHC15g1a}{LHC11a - EMC1}";
    } else if(select.CompareTo("LHC11a-Phojet-ConvCalo")==0){
        strDataFile[0]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20151213-CalibDefaultTiming/GammaConvCalo_LHC11a-pass4_28.root";
        strMCFile[0]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20151213-CalibDefaultTiming/GammaConvCalo_MC_LHC12f1b_28.root";
        dataCut[0]          = "00003113_00200009327000008250400000_1111100053032230000_0163103100000010";
        mcCut[0]            = "00003113_00200009327000008250400000_1111100053032230000_0163103100000010";
        dataMainDir[0]      = "GammaConvCalo";
        mcMainDir[0]        = "GammaConvCalo";

        firstTriggerBin[0]  = 13;
        strDataFile[1]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20151213-CalibDefaultTiming/GammaConvCalo_LHC11a-pass4_28.root";
        strMCFile[1]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20151213-CalibDefaultTiming/GammaConvCalo_MC_LHC15g1aFinerPtHardBins_34.root";
        dataCut[1]          = "00051113_00200009327000008250400000_1111100053032230000_0163103100000010";
        mcCut[1]            = "00051113_00200009327000008250400000_1111100053032230000_0163103100000010";

        fPlot[0]            = "#frac{LHC12f1b}{LHC11a - INT1}";
        fPlot[1]            = "#frac{LHC15g1a}{LHC11a - EMC1}";
    } else if(select.CompareTo("LHC11a-JetJet-ConvCalo")==0){
        strDataFile[0]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20151213-CalibDefaultTiming/GammaConvCalo_LHC11a-pass4_28.root";
        strMCFile[0]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20151213-CalibDefaultTiming/GammaConvCalo_MC_LHC15g1aFinerPtHardBins_34.root";
        dataCut[0]          = "00003113_00200009327000008250400000_1111100053032230000_0163103100000010";
        mcCut[0]            = "00003113_00200009327000008250400000_1111100053032230000_0163103100000010";
        dataMainDir[0]      = "GammaConvCalo";
        mcMainDir[0]        = "GammaConvCalo";

        firstTriggerBin[0]  = 13;
        strDataFile[1]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20151213-CalibDefaultTiming/GammaConvCalo_LHC11a-pass4_28.root";
        strMCFile[1]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20151213-CalibDefaultTiming/GammaConvCalo_MC_LHC15g1aFinerPtHardBins_34.root";
        dataCut[1]          = "00051113_00200009327000008250400000_1111100053032230000_0163103100000010";
        mcCut[1]            = "00003113_00200009327000008250400000_1111100053032230000_0163103100000010";

        fPlot[0]            = "#frac{LHC15g1a}{LHC11a - INT1}";
        fPlot[1]            = "#frac{LHC15g1a}{LHC11a - EMC1}";
        
    // Alternate configuration for LHC11a with open timing cuts in tender
    } else if(select.CompareTo("LHC11a-Phojet-Calo-OpenTime")==0){
        startPtBin          = 4;
        strDataFile[0]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_LHC11a-pass4_40.root";
        strMCFile[0]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_MC_LHC12f1b_40.root";
        dataCut[0]          = "00003113_1111100053032220000_0163103100000050";
        mcCut[0]            = "00003113_1111100053032220000_0163103100000050";

        firstTriggerBin[0]  = 13;
        strDataFile[1]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_LHC11a-pass4_40.root";
        strMCFile[1]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_MC_LHC15g1aFinerPtHardBins_40.root";
        dataCut[1]          = "00051113_1111100053032220000_0163103100000050";
        mcCut[1]            = "00003113_1111100053032220000_0163103100000050";

        fPlot[0]            = "#frac{LHC12f1b}{LHC11a - INT1}";
        fPlot[1]            = "#frac{LHC15g1a}{LHC11a - EMC1}";
    } else if(select.CompareTo("LHC11a-Pythia-Calo-OpenTime")==0){
        startPtBin          = 4;
        strDataFile[0]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_LHC11a-pass4_40.root";
        strMCFile[0]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_MC_LHC12f1a_LHC12i3_40.root";
        dataCut[0]          = "00003113_1111100053032220000_0163103100000050";
        mcCut[0]            = "00003113_1111100053032220000_0163103100000050";

        firstTriggerBin[0]  = 13;
        strDataFile[1]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_LHC11a-pass4_40.root";
        strMCFile[1]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_MC_LHC15g1aFinerPtHardBins_40.root";
        dataCut[1]          = "00051113_1111100053032220000_0163103100000050";
        mcCut[1]            = "00003113_1111100053032220000_0163103100000050";

        fPlot[0]            = "#frac{LHC12f1a & LHC12i3}{LHC11a - INT1}";
        fPlot[1]            = "#frac{LHC15g1a}{LHC11a - EMC1}";
    } else if(select.CompareTo("LHC13g-Pythia-Calo-OpenTime")==0){
        startPtBin          = 4;
        endPtBin            = 17;       
        strDataFile[0]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_LHC13g-pass1_42.root";
        strMCFile[0]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_MC_LHC15g2_42.root";
        dataCut[0]          = "00000113_1111100063032220000_0163103100000050";
        mcCut[0]            = "00000113_1111100063032220000_0163103100000050";

        firstTriggerBin[0]  = 11;
        strDataFile[1]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_LHC13g-pass1_42.root";
//         strMCFile[1]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_MC_LHC15g2_42.root";
        strMCFile[1]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_42.root";
        dataCut[1]          = "00052113_1111100063032220000_0163103100000050";
        mcCut[1]            = "00000113_1111100063032220000_0163103100000050";

        firstTriggerBin[1]  = 14;
        strDataFile[2]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_LHC13g-pass1_42.root";
//         strMCFile[2]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_MC_LHC15g2_42.root";
        strMCFile[2]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_42.root";
        dataCut[2]          = "00085113_1111100063032220000_0163103100000050";
        mcCut[2]            = "00000113_1111100063032220000_0163103100000050";

//         firstTriggerBin[2]  = 17;
//         strDataFile[3]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_LHC13g-pass1_42.root";
//         strMCFile[3]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_42.root";
//         dataCut[3]          = "00083113_1111100063032220000_0163103100000050";
//         mcCut[3]            = "00000113_1111100063032220000_0163103100000050";
        
        fPlot[0]            = "#frac{LHC15g2}{LHC13g - INT7}";
        fPlot[1]            = "#frac{LHC15a3a/+}{LHC13g - EMC7}";
//         fPlot[1]            = "#frac{LHC15g2}{LHC13g - EMC7}";
        fPlot[2]            = "#frac{LHC15a3a/+}{LHC13g - EG2}";
//         fPlot[2]            = "#frac{LHC15g2}{LHC13g - EG2}";
//         fPlot[3]            = "#frac{LHC15a3a/+}{LHC13g - EG1}";
    } else if(select.CompareTo("LHC13g-JetJet-Calo-OpenTime")==0){
        startPtBin          = 4;
        endPtBin            = 17;       
        strDataFile[0]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_LHC13g-pass1_42.root";
        strMCFile[0]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_42.root";
        dataCut[0]          = "00000113_1111100063032220000_0163103100000050";
        mcCut[0]            = "00000113_1111100063032220000_0163103100000050";

        firstTriggerBin[0]  = 11;
        strDataFile[1]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_LHC13g-pass1_42.root";
        strMCFile[1]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_42.root";
        dataCut[1]          = "00052113_1111100063032220000_0163103100000050";
        mcCut[1]            = "00000113_1111100063032220000_0163103100000050";

        firstTriggerBin[1]  = 14;
        strDataFile[2]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_LHC13g-pass1_42.root";
        strMCFile[2]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_42.root";
        dataCut[2]          = "00085113_1111100063032220000_0163103100000050";
        mcCut[2]            = "00000113_1111100063032220000_0163103100000050";

        fPlot[0]            = "#frac{LHC15a3a/+}{LHC13g - INT7}";
        fPlot[1]            = "#frac{LHC15a3a/+}{LHC13g - EMC7}";
        fPlot[2]            = "#frac{LHC15a3a/+}{LHC13g - EG2}";
        
    } else if(select.CompareTo("LHC11a-JetJet-Calo-OpenTime")==0){
        startPtBin          = 4;
        strDataFile[0]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_LHC11a-pass4_40.root";
        strMCFile[0]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_MC_LHC15g1aFinerPtHardBins_40.root";
        dataCut[0]          = "00003113_1111100053032220000_0163103100000050";
        mcCut[0]            = "00003113_1111100053032220000_0163103100000050";

        firstTriggerBin[0]  = 13;
        strDataFile[1]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_LHC11a-pass4_40.root";
        strMCFile[1]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaCalo_MC_LHC15g1aFinerPtHardBins_40.root";
        dataCut[1]          = "00051113_1111100053032220000_0163103100000050";
        mcCut[1]            = "00003113_1111100053032220000_0163103100000050";

        fPlot[0]            = "#frac{LHC15g1a}{LHC11a - INT1}";
        fPlot[1]            = "#frac{LHC15g1a}{LHC11a - EMC1}";

        
    } else if(select.CompareTo("LHC11a-Phojet-ConvCalo-OpenTime")==0){
        startPtBin          = 3;
        strDataFile[0]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_LHC11a-pass4_34.root";
        strMCFile[0]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_MC_LHC12f1b_34.root";
        dataCut[0]          = "00003113_00200009327000008250400000_1111100053032230000_0163103100000010";
        mcCut[0]            = "00003113_00200009327000008250400000_1111100053032230000_0163103100000010";

        firstTriggerBin[0]  = 13;
        strDataFile[1]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_LHC11a-pass4_34.root";
        strMCFile[1]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_MC_LHC15g1aFinerPtHardBins_34.root";
        dataCut[1]          = "00051113_00200009327000008250400000_1111100053032230000_0163103100000010";
        mcCut[1]            = "00003113_00200009327000008250400000_1111100053032230000_0163103100000010";

        fPlot[0]            = "#frac{LHC12f1b}{LHC11a - INT1}";
        fPlot[1]            = "#frac{LHC15g1a}{LHC11a - EMC1}";
    } else if(select.CompareTo("LHC11a-Pythia-ConvCalo-OpenTime")==0){
        startPtBin          = 3;
        strDataFile[0]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_LHC11a-pass4_34.root";
        strMCFile[0]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_MC_LHC12f1a_34.root";
        dataCut[0]          = "00003113_00200009327000008250400000_1111100053032230000_0163103100000010";
        mcCut[0]            = "00003113_00200009327000008250400000_1111100053032230000_0163103100000010";

        firstTriggerBin[0]  = 12;
        strDataFile[1]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_LHC11a-pass4_34.root";
        strMCFile[1]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_MC_LHC15g1aFinerPtHardBins_34.root";
        dataCut[1]          = "00051113_00200009327000008250400000_1111100053032230000_0163103100000010";
        mcCut[1]            = "00003113_00200009327000008250400000_1111100053032230000_0163103100000010";

        fPlot[0]            = "#frac{LHC12f1a & LHC12i3}{LHC11a - INT1}";
        fPlot[1]            = "#frac{LHC15g1a}{LHC11a - EMC1}";
    } else if(select.CompareTo("LHC11a-JetJet-ConvCalo-OpenTime")==0){
        startPtBin          = 3;
        strDataFile[0]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_LHC11a-pass4_34.root";
        strMCFile[0]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_MC_LHC15g1aFinerPtHardBins_34.root";
        dataCut[0]          = "00003113_00200009327000008250400000_1111100053032230000_0163103100000010";
        mcCut[0]            = "00003113_00200009327000008250400000_1111100053032230000_0163103100000010";

        firstTriggerBin[0]  = 12;
        strDataFile[1]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_LHC11a-pass4_34.root";
        strMCFile[1]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_MC_LHC15g1aFinerPtHardBins_34.root";
        dataCut[1]          = "00051113_00200009327000008250400000_1111100053032230000_0163103100000010";
        mcCut[1]            = "00003113_00200009327000008250400000_1111100053032230000_0163103100000010";

        fPlot[0]            = "#frac{LHC15g1a}{LHC11a - INT1}";
        fPlot[1]            = "#frac{LHC15g1a}{LHC11a - EMC1}";
    } else if(select.CompareTo("LHC13g-Pythia-ConvCalo-OpenTime")==0){
        startPtBin          = 3;
        strDataFile[0]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_LHC13g-pass3_33.root";
        strMCFile[0]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_MC_LHC15g2_33.root";
        dataCut[0]          = "00000113_00200009327000008250400000_1111100063032230000_0163103100000010";
        mcCut[0]            = "00000113_00200009327000008250400000_1111100063032230000_0163103100000010";

        firstTriggerBin[0]  = 11;
        strDataFile[1]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_LHC13g-pass3_33.root";
        strMCFile[1]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_MC_LHC15g2_33.root";
//         strMCFile[1]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_33.root";
        dataCut[1]          = "00052113_00200009327000008250400000_1111100063032230000_0163103100000010";
        mcCut[1]            = "00000113_00200009327000008250400000_1111100063032230000_0163103100000010";

        firstTriggerBin[1]  = 15;
        strDataFile[2]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_LHC13g-pass3_33.root";
        strMCFile[2]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_MC_LHC15g2_33.root";
//         strMCFile[2]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_33.root";
        dataCut[2]          = "00085113_00200009327000008250400000_1111100063032230000_0163103100000010";
        mcCut[2]            = "00000113_00200009327000008250400000_1111100063032230000_0163103100000010";

        firstTriggerBin[2]  = 18;
        strDataFile[3]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_LHC13g-pass3_33.root";
        strMCFile[3]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_33.root";
        dataCut[3]          = "00083113_00200009327000008250400000_1111100063032230000_0163103100000010";
        mcCut[3]            = "00000113_00200009327000008250400000_1111100063032230000_0163103100000010";
        
        fPlot[0]            = "#frac{LHC15g2}{LHC13g - INT7}";
//         fPlot[1]            = "#frac{LHC15a3a/+}{LHC13g - EMC7}";
        fPlot[1]            = "#frac{LHC15g2}{LHC13g - EMC7}";
//         fPlot[2]            = "#frac{LHC15a3a/+}{LHC13g - EG2}";
        fPlot[2]            = "#frac{LHC15g2}{LHC13g - EG2}";
        fPlot[3]            = "#frac{LHC15a3a/+}{LHC13g - EG1}";
        
//         endPtBin            = 17;        
    } else if(select.CompareTo("LHC13g-JetJet-ConvCalo-OpenTime")==0){
        startPtBin          = 3;
        strDataFile[0]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_LHC13g-pass3_33.root";
        strMCFile[0]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_33.root";
        dataCut[0]          = "00000113_00200009327000008250400000_1111100063032230000_0163103100000010";
        mcCut[0]            = "00000113_00200009327000008250400000_1111100063032230000_0163103100000010";

        firstTriggerBin[0]  = 11;
        strDataFile[1]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_LHC13g-pass3_33.root";
//         strMCFile[1]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_MC_LHC15g2_33.root";
        strMCFile[1]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_33.root";
        dataCut[1]          = "00052113_00200009327000008250400000_1111100063032230000_0163103100000010";
        mcCut[1]            = "00000113_00200009327000008250400000_1111100063032230000_0163103100000010";

        firstTriggerBin[1]  = 15;
        strDataFile[2]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_LHC13g-pass3_33.root";
//         strMCFile[2]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_MC_LHC15g2_33.root";
        strMCFile[2]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_33.root";
        dataCut[2]          = "00085113_00200009327000008250400000_1111100063032230000_0163103100000010";
        mcCut[2]            = "00000113_00200009327000008250400000_1111100063032230000_0163103100000010";

        firstTriggerBin[2]  = 18;
        strDataFile[3]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_LHC13g-pass3_33.root";
        strMCFile[3]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTiming/GammaConvCalo_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_33.root";
        dataCut[3]          = "00083113_00200009327000008250400000_1111100063032230000_0163103100000010";
        mcCut[3]            = "00000113_00200009327000008250400000_1111100063032230000_0163103100000010";
        
        fPlot[0]            = "#frac{LHC15a3a/+}{LHC13g - INT7}";
        fPlot[1]            = "#frac{LHC15a3a/+}{LHC13g - EMC7}";
//         fPlot[1]            = "#frac{LHC15g2}{LHC13g - EMC7}";
        fPlot[2]            = "#frac{LHC15a3a/+}{LHC13g - EG2}";
//         fPlot[2]            = "#frac{LHC15g2}{LHC13g - EG2}";
        fPlot[3]            = "#frac{LHC15a3a/+}{LHC13g - EG1}";
        
//  
    // Alternate configuration for LHC11a with open timing cuts in tender, but delta T cut
    } else if(select.CompareTo("LHC11a-Pythia-ConvCalo-OpenTimeDeltaT")==0){
        startPtBin          = 3;
        strDataFile[0]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTimingDeltaT/GammaConvCalo_LHC11a-pass4_34.root";
        strMCFile[0]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTimingDeltaT/GammaConvCalo_MC_LHC12f1a_34.root";
        dataCut[0]          = "00003113_00200009327000008250400000_1111100053032230000_0163103100000010";
        mcCut[0]            = "00003113_00200009327000008250400000_1111100053032230000_0163103100000010";

        firstTriggerBin[0]  = 13;
        strDataFile[1]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTimingDeltaT/GammaConvCalo_LHC11a-pass4_34.root";
        strMCFile[1]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTimingDeltaT/GammaConvCalo_MC_LHC15g1aFinerPtHardBins_34.root";
        dataCut[1]          = "00051113_00200009327000008250400000_1111100053032230000_0163103100000010";
        mcCut[1]            = "00003113_00200009327000008250400000_1111100053032230000_0163103100000010";

        fPlot[0]            = "#frac{LHC12f1a & LHC12i3}{LHC11a - INT1}";
        fPlot[1]            = "#frac{LHC15g1a}{LHC11a - EMC1}";
    } else if(select.CompareTo("LHC11a-Phojet-Calo-OpenTimeDeltaT")==0){
        startPtBin          = 3;
        strDataFile[0]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTimingDeltaT/GammaCalo_LHC11a-pass4_40.root";
        strMCFile[0]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTimingDeltaT/GammaCalo_MC_LHC12f1b_40.root";
        dataCut[0]          = "00003113_1111100053032220000_0163103100000050";
        mcCut[0]            = "00003113_1111100053032220000_0163103100000050";

        firstTriggerBin[0]  = 12;
        strDataFile[1]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTimingDeltaT/GammaCalo_LHC11a-pass4_40.root";
        strMCFile[1]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTimingDeltaT/GammaCalo_MC_LHC15g1aFinerPtHardBins_40.root";
        dataCut[1]          = "00051113_1111100053032220000_0163103100000050";
        mcCut[1]            = "00003113_1111100053032220000_0163103100000050";

        fPlot[0]            = "#frac{LHC12f1b}{LHC11a - INT1}";
        fPlot[1]            = "#frac{LHC15g1a}{LHC11a - EMC1}";
    } else if(select.CompareTo("LHC11a-Pythia-Calo-OpenTimeDeltaT")==0){
        startPtBin          = 3;
        strDataFile[0]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTimingDeltaT/GammaCalo_LHC11a-pass4_40.root";
        strMCFile[0]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTimingDeltaT/GammaCalo_MC_LHC12f1a_LHC12i3_40.root";
        dataCut[0]          = "00003113_1111100053032220000_0163103100000050";
        mcCut[0]            = "00003113_1111100053032220000_0163103100000050";

        firstTriggerBin[0]  = 12;
        strDataFile[1]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTimingDeltaT/GammaCalo_LHC11a-pass4_40.root";
        strMCFile[1]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTimingDeltaT/GammaCalo_MC_LHC15g1aFinerPtHardBins_40.root";
        dataCut[1]          = "00051113_1111100053032220000_0163103100000050";
        mcCut[1]            = "00003113_1111100053032220000_0163103100000050";

        fPlot[0]            = "#frac{LHC12f1a & LHC12i3}{LHC11a - INT1}";
        fPlot[1]            = "#frac{LHC15g1a}{LHC11a - EMC1}";
    }else if(select.CompareTo("LHC11a-Phojet-ConvCalo-OpenTimeDeltaT")==0){
        startPtBin          = 3;
        strDataFile[0]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTimingDeltaT/GammaConvCalo_LHC11a-pass4_34.root";
        strMCFile[0]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTimingDeltaT/GammaConvCalo_MC_LHC12f1b_34.root";
        dataCut[0]          = "00003113_00200009327000008250400000_1111100053032230000_0163103100000010";
        mcCut[0]            = "00003113_00200009327000008250400000_1111100053032230000_0163103100000010";

        firstTriggerBin[0]  = 13;
        strDataFile[1]      = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTimingDeltaT/GammaConvCalo_LHC11a-pass4_34.root";
        strMCFile[1]        = "/mnt/additionalStorage/OutputLegoTrains/pp/Legotrain-vAN-20160108-CalibOpenTimingDeltaT/GammaConvCalo_MC_LHC15g1aFinerPtHardBins_34.root";
        dataCut[1]          = "00051113_00200009327000008250400000_1111100053032230000_0163103100000010";
        mcCut[1]            = "00003113_00200009327000008250400000_1111100053032230000_0163103100000010";

        fPlot[0]            = "#frac{LHC12f1b}{LHC11a - INT1}";
        fPlot[1]            = "#frac{LHC15g1a}{LHC11a - EMC1}";
    // Default configuration for LHC10 (7TeV) 
    } else if(select.CompareTo("LHC10-Calo")==0){
        strDataFile[0]      = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150922-LHC10_p4-QA/GammaCalo_LHC10bcdef-pass4_201.root";
        strMCFile[0]        = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150922-LHC10_p4-QA/GammaCalo_MC_LHC14j4b_LHC14j4c_LHC14j4d_LHC14j4e_LHC14j4f_201.root";
        dataCut[0]          = "00000113_1111100010032230000_0163103100000050";
        mcCut[0]            = "00000113_1111100010032230000_0163103100000050";

        fPlot[0]            = "#frac{LHC14j4b-f}{LHC10b-f}";
    } else if(select.CompareTo("LHC10-ConvCalo")==0){
        strDataFile[0]      = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150922-LHC10_p4-QA/GammaConvCalo_LHC10bcdef-pass4_201.root";
        strMCFile[0]        = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20150922-LHC10_p4-QA/GammaConvCalo_MC_LHC14j4b_LHC14j4c_LHC14j4d_LHC14j4e_LHC14j4f_201.root";
        dataCut[0]          = "00000113_00200009327000008250400000_1111100013032230000_0163103100000010";
        mcCut[0]            = "00000113_00200009327000008250400000_1111100013032230000_0163103100000010";

        fPlot[0]            = "#frac{LHC14j4b-f}{LHC10b-f}";
    } else if(select.CompareTo("LHC13bc-Calo")==0){
        strDataFile[0]      = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20151005-LHC13bc_NonLin/GammaCalo_LHC13b-pass3_LHC13c-pass3_13.root";
        strMCFile[0]        = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20151005-LHC13bc_NonLin/GammaCalo_MC_LHC13b2_efix_p1_p2_p3_p4_13.root";
        dataCut[0]          = "80000013_1111100050032230000_0163103100000050";
        mcCut[0]            = "80000013_1111100050032230000_0163103100000050";

        fPlot[0]            = "#frac{LHC13b2_efix}{LHC13bc}";
    } else if(select.CompareTo("LHC13bc-Calo2")==0){
        strDataFile[0]      = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20151005-LHC13bc_NonLin/GammaCalo_LHC13b-pass3_LHC13c-pass3_13.root";
        strMCFile[0]        = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20151005-LHC13bc_NonLin/GammaCalo_MC_LHC13b2_efix_p1_p2_p3_p4_13.root";
        dataCut[0]          = "80000013_1111100050032230000_0163103100000050";
        mcCut[0]            = "80000013_1111100050032230000_0163103100000050";

        fPlot[0]            = "#frac{LHC13b2_efix}{LHC13bc}";
        startPtBin          = 2;
        endPtBin            = 17;
    } else if(select.CompareTo("LHC13bc-ConvCalo")==0){
        strDataFile[0]      = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20151005-LHC13bc_NonLin/GammaConvCalo_LHC13b-pass4_LHC13c-pass4_3.root";
        strMCFile[0]        = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20151005-LHC13bc_NonLin/GammaConvCalo_MC_LHC13b2_efix_p1_p2_p3_p4_3.root";
        dataCut[0]          = "80000013_00200009327002008250400000_1111100053032230000_0163103100000010";
        mcCut[0]            = "80000013_00200009327002008250400000_1111100053032230000_0163103100000010";
        fPlot[0]            = "#frac{LHC13b2_efix}{LHC13bc}";
    } else if(select.CompareTo("LHC12-Pythia-Calo")==0){
      strDataFile[0]      = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20160120-8TeV_NL/kINT7-LHC12_GammaCalo_111.root";
      strMCFile[0]        = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20160120-8TeV_NL/kINT7-LHC15h1_GammaCalo_111.root";
      dataCut[0]          = "00000113_1111100063032230000_0163103100000050";
      mcCut[0]            = "00000113_1111100063032230000_0163103100000050";

      firstTriggerBin[0]  = 12;
      strDataFile[1]      = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20160120-8TeV_NL/kEMC7-LHC12_GammaCalo_111.root";
      strMCFile[1]        = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20160120-8TeV_NL/kINT7-LHC15h1_GammaCalo_111.root";
      dataCut[1]          = "00052113_1111100063032230000_0163103100000050";
      mcCut[1]            = "00052113_1111100063032230000_0163103100000050";

      fPlot[0]            = "#frac{LHC15h1}{LHC12 - INT7}";
      fPlot[1]            = "#frac{LHC15h1}{LHC12 - EMC7}";
    } else if(select.CompareTo("LHC12-Pythia-ConvCalo")==0){
      strDataFile[0]      = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20160120-8TeV_NL/kINT7-LHC12_GammaConvCalo_31.root";
      strMCFile[0]        = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20160120-8TeV_NL/kINT7-LHC15h1_GammaConvCalo_31.root";
      dataCut[0]          = "00000113_00200009327000008250400000_1111100063032230000_0163103100000010";
      mcCut[0]            = "00000113_00200009327000008250400000_1111100063032230000_0163103100000010";

      firstTriggerBin[0]  = 12;
      strDataFile[1]      = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20160120-8TeV_NL/kEMC7-LHC12_GammaConvCalo_31.root";
      strMCFile[1]        = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20160120-8TeV_NL/kINT7-LHC15h1_GammaConvCalo_31.root";
      dataCut[1]          = "00052113_00200009327000008250400000_1111100063032230000_0163103100000010";
      mcCut[1]            = "00052113_00200009327000008250400000_1111100063032230000_0163103100000010";

      fPlot[0]            = "#frac{LHC15h1}{LHC12 - INT7}";
      fPlot[1]            = "#frac{LHC15h1}{LHC12 - EMC7}";
    } else if(select.CompareTo("LHC12-Phojet-Calo")==0){
      strDataFile[0]      = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20160120-8TeV_NL/kINT7-LHC12_GammaCalo_111.root";
      strMCFile[0]        = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20160120-8TeV_NL/kINT7-LHC15h2_GammaCalo_111.root";
      dataCut[0]          = "00000113_1111100063032230000_0163103100000050";
      mcCut[0]            = "00000113_1111100063032230000_0163103100000050";

      firstTriggerBin[0]  = 12;
      strDataFile[1]      = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20160120-8TeV_NL/kEMC7-LHC12_GammaCalo_111.root";
      strMCFile[1]        = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20160120-8TeV_NL/kINT7-LHC15h2_GammaCalo_111.root";
      dataCut[1]          = "00052113_1111100063032230000_0163103100000050";
      mcCut[1]            = "00052113_1111100063032230000_0163103100000050";

      fPlot[0]            = "#frac{LHC15h2}{LHC12 - INT7}";
      fPlot[1]            = "#frac{LHC15h1}{LHC12 - EMC7}";
    } else if(select.CompareTo("LHC12-Phojet-ConvCalo")==0){
      strDataFile[0]      = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20160120-8TeV_NL/kINT7-LHC12_GammaConvCalo_31.root";
      strMCFile[0]        = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20160120-8TeV_NL/kINT7-LHC15h2_GammaConvCalo_31.root";
      dataCut[0]          = "00000113_00200009327000008250400000_1111100063032230000_0163103100000010";
      mcCut[0]            = "00000113_00200009327000008250400000_1111100063032230000_0163103100000010";
      dataMainDir[0]      = "GammaConvCalo";
      mcMainDir[0]        = "GammaConvCalo";

      firstTriggerBin[0]  = 12;
      strDataFile[1]      = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20160120-8TeV_NL/kEMC7-LHC12_GammaConvCalo_31.root";
      strMCFile[1]        = "/home/daniel/Desktop/work/Grid/Legotrain-vAN-20160120-8TeV_NL/kINT7-LHC15h2_GammaConvCalo_31.root";
      dataCut[1]          = "00052113_00200009327000008250400000_1111100063032230000_0163103100000010";
      mcCut[1]            = "00052113_00200009327000008250400000_1111100063032230000_0163103100000010";

      fPlot[0]            = "#frac{LHC15h2}{LHC12 - INT7}";
      fPlot[1]            = "#frac{LHC15h2}{LHC12 - EMC7}";
    } else{
        cout << "No valid selection '" << select.Data() << "'' given, returning..." << endl;
        return;
    }
//*******************************************************************************

    // binning for fits Data vs MC mean mass position
    const Int_t fNBins          = 24;
    Double_t fBins[fNBins+1]    = { 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4, 1.6,
                                    1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0,
                                    5.0, 6.0, 8.0, 10.0, 12.0, 16.0,
                                    20.0, 25.0, 30.0};

    Int_t numberOfTriggers      = 0;
    fPlot[0] = Form("#it{E}_{Cluster} < %0.1f GeV : %s",fBins[firstTriggerBin[0]],fPlot[0].Data());
    while (firstTriggerBin[numberOfTriggers]>-1){
        if (firstTriggerBin[numberOfTriggers+1]>-1)
            fPlot[numberOfTriggers+1] = Form("%0.1f GeV #leq #it{E}_{Cluster} < %0.1f GeV : %s",fBins[firstTriggerBin[numberOfTriggers]], fBins[firstTriggerBin[numberOfTriggers+1]],
                                             fPlot[numberOfTriggers+1].Data());
        else 
            fPlot[numberOfTriggers+1] = Form("#it{E}_{Cluster} #geq %0.1f GeV : %s",fBins[firstTriggerBin[numberOfTriggers]],fPlot[numberOfTriggers+1].Data());
        numberOfTriggers++;
    }

    TString recGamma            = "";
    if(select.Contains("-Calo")) 
        recGamma                = "#gamma's rec. with EMCal";
    else if(select.Contains("-ConvCalo")) 
        recGamma                = "#gamma's rec. with PCM, EMCal";

    TString fTextMeasurement    = Form("#pi^{0} #rightarrow #gamma#gamma");
    TString fCollisionSystem    = ReturnFullCollisionsSystem(optionEnergy);
    if (fCollisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }
//*******************************************************************************
// Input
    TFile* dataFile             = new TFile(strDataFile[0].Data(),"READ"); 
    if(dataFile->IsZombie()) {cout << "Info: ROOT file '" << strDataFile[0].Data() << "' could not be openend, return!" << endl; return;}
    TList* dataTopDir           = (TList*) dataFile->Get(dataMainDir[0].Data()); 
    if(dataTopDir == NULL) {cout << "ERROR: dataTopDir not Found"<<endl; return;}
    TList* dataTopContainer     = (TList*) dataTopDir->FindObject(Form("Cut Number %s",dataCut[0].Data())); 
    if(dataTopContainer == NULL) {cout << "ERROR: " << Form("Cut Number '%s'",dataCut[0].Data()) << " not found in Data-File" << endl; return;}
    TList* dataESDContainer     = (TList*) dataTopContainer->FindObject(Form("%s ESD histograms",dataCut[0].Data())); 
    if(dataESDContainer == NULL) {cout << "ERROR: " << Form("'%s' ESD histograms",dataCut[0].Data()) << " not found in Data-File" << endl; return;}

    TFile* mcFile               = new TFile(strMCFile[0].Data(),"READ"); 
    if(mcFile->IsZombie()) {cout << "Info: ROOT file '" << strMCFile[0].Data() << "' could not be openend, return!" << endl; return;}
    TList* mcTopDir             = (TList*) mcFile->Get(mcMainDir[0].Data()); 
    if(mcTopDir == NULL) {cout << "ERROR: mcTopDir not Found"<<endl; return;}
    TList* mcTopContainer       = (TList*) mcTopDir->FindObject(Form("Cut Number %s",mcCut[0].Data())); 
    if(mcTopContainer == NULL) {cout << "ERROR: " << Form("Cut Number '%s'",mcCut[0].Data()) << " not found in MC-File" << endl; return;}
    TList* mcESDContainer       = (TList*) mcTopContainer->FindObject(Form("%s ESD histograms",mcCut[0].Data())); 
    if(mcESDContainer == NULL) {cout << "ERROR: " << Form("'%s' ESD histograms",mcCut[0].Data()) << " not found in MC-File" << endl; return;}
    //*******************************************************************************
    // Output
    TString nameOutput          = "";
    nameOutput                  = Form("%s/CorrectCaloNonLinearity3_%s.root",outputDir.Data(),select.Data());
    TFile* fOutput              = new TFile(nameOutput,"RECREATE");
    //*******************************************************************************
    // Fitting Data+MC
    TH2F* dataInvMassPtAlpha    = NULL;
    TH2F* dataBGInvMassPtAlpha  = NULL;
    TH2F* mcInvMassPtAlpha      = NULL;
    TH2F* mcBGInvMassPtAlpha    = NULL;
    if(mode==2||mode==3){
        cout << "entered Conv Calo mode" << endl;
        dataInvMassPtAlpha      = (TH2F*) dataESDContainer->FindObject("ESD_Mother_InvMass_E_Calib");
        dataBGInvMassPtAlpha    = (TH2F*) dataESDContainer->FindObject("ESD_Background_InvMass_E_Calib");
        mcInvMassPtAlpha        = (TH2F*) mcESDContainer->FindObject("ESD_Mother_InvMass_E_Calib");
        mcBGInvMassPtAlpha      = (TH2F*) mcESDContainer->FindObject("ESD_Background_InvMass_E_Calib");
    }else{
        cout << "entered Calo mode" << endl;
        dataInvMassPtAlpha      = (TH2F*) dataESDContainer->FindObject("ESD_Mother_InvMass_vs_Pt_Alpha");
        dataBGInvMassPtAlpha    = (TH2F*) dataESDContainer->FindObject("ESD_Background_InvMass_vs_Pt_Alpha");
        mcInvMassPtAlpha        = (TH2F*) mcESDContainer->FindObject("ESD_Mother_InvMass_vs_Pt_Alpha");
        mcBGInvMassPtAlpha      = (TH2F*) mcESDContainer->FindObject("ESD_Background_InvMass_vs_Pt_Alpha");
    }
    if(!dataInvMassPtAlpha){cout << "did not find ESD_Mother_InvMass_E_Calib in data" << endl; return;}
    if(!dataBGInvMassPtAlpha){cout << "did not find ESD_Background_InvMass_E_Calib in data" << endl; return;}
    if(!mcInvMassPtAlpha){cout << "did not find ESD_Mother_InvMass_E_Calib in mc" << endl; return;}
    if(!mcBGInvMassPtAlpha){cout << "did not find ESD_Background_InvMass_E_Calib in mc" << endl; return;}
    
    dataInvMassPtAlpha->Write("Data - ESD_Mother_InvMass");
    dataBGInvMassPtAlpha->Write("Data - ESD_Background_InvMass");
    mcInvMassPtAlpha->Write("MC - ESD_Mother_InvMass");
    mcBGInvMassPtAlpha->Write("MC - ESD_Background_InvMass");
    
    TH1D* histMCResults         = NULL;
    TH1D* histDataResults       = NULL;
    TH1D* histDataMCResults     = NULL;

    histMCResults               = new TH1D("Mean mass MC","; #it{E}_{Cluster} (GeV); mean mass MC",fNBins,fBins);
    histDataResults             = new TH1D("Mean mass Data","; #it{E}_{Cluster} (GeV); mean mass Data",fNBins,fBins);
    histDataMCResults           = new TH1D("Mean mass ratio MC/Data","; #it{E}_{Cluster} (GeV); mean mass ratio (MC/Data)",fNBins,fBins);

    histMCResults->SetDirectory(0);
    histDataResults->SetDirectory(0);
    histMCResults->GetXaxis()->SetRangeUser(fBins[startPtBin],fBins[endPtBin]);
    histDataResults->GetXaxis()->SetRangeUser(fBins[startPtBin],fBins[endPtBin]);
    histMCResults->GetYaxis()->SetRangeUser(0.12,0.14);
    histDataResults->GetYaxis()->SetRangeUser(0.12,0.14);

    TF1* fFitReco;
    TF1* fFitMassPos;

    TCanvas *canvas             = new TCanvas("canvas","",200,0,1350,900);  // gives the page size
    DrawGammaCanvasSettings(canvas, 0.1, 0.02, 0.06, 0.1);

    TString dataMC[2]           = {"Data","MC"};
    Int_t triggerSel            = 0;
    for(Int_t iClusterPt=startPtBin; iClusterPt<endPtBin; iClusterPt++){
        if(iClusterPt==firstTriggerBin[triggerSel]){
            // switching to trigger
            cout << endl;
            cout << "-----------------------------------------------------" << endl;
            cout << "\t Closing open files, switching to Trigger!" << endl;
            cout << "bin: " << firstTriggerBin[triggerSel] << endl;
            cout << "-----------------------------------------------------" << endl;

            if(!strDataFile[triggerSel+1].IsNull()){
                dataESDContainer->Clear(); dataTopContainer->Clear(); dataTopDir->Clear(); dataFile->Delete();
                dataFile            = new TFile(strDataFile[triggerSel+1].Data(),"READ"); 
                if(dataFile->IsZombie()) {cout << "Info: ROOT file '" << strDataFile[triggerSel+1].Data() << "' could not be openend, return!" << endl; return;}
                dataTopDir          = (TList*) dataFile->Get(dataMainDir[triggerSel+1].Data()); 
                if(dataTopDir == NULL) {cout << "ERROR: dataTopDir not Found"<<endl; return;}
                dataTopContainer    = (TList*) dataTopDir->FindObject(Form("Cut Number %s",dataCut[triggerSel+1].Data())); 
                if(dataTopContainer == NULL) {cout << "ERROR: " << Form("Cut Number '%s'",dataCut[triggerSel+1].Data()) << " not found in File" << endl; return;}
                dataESDContainer    = (TList*) dataTopContainer->FindObject(Form("%s ESD histograms",dataCut[triggerSel+1].Data())); 
                if(dataESDContainer == NULL) {cout << "ERROR: " << Form("'%s' ESD histograms",dataCut[triggerSel+1].Data()) << " not found in File" << endl; return;}
            }
            if(!strMCFile[triggerSel+1].IsNull()){
                mcESDContainer->Clear(); mcTopContainer->Clear(); mcTopDir->Clear(); mcFile->Delete();
                mcFile              = new TFile(strMCFile[triggerSel+1].Data(),"READ"); 
                if(mcFile->IsZombie()) {cout << "Info: ROOT file '" << strMCFile[triggerSel+1].Data() << "' could not be openend, return!" << endl; return;}
                mcTopDir            = (TList*) mcFile->Get(mcMainDir[triggerSel+1].Data()); 
                if(mcTopDir == NULL) {cout << "ERROR: mcTopDir not Found"<<endl; return;}
                mcTopContainer      = (TList*) mcTopDir->FindObject(Form("Cut Number %s",mcCut[triggerSel+1].Data())); 
                if(mcTopContainer == NULL) {cout << "ERROR: " << Form("Cut Number '%s'",mcCut[triggerSel+1].Data()) << " not found in File" << endl; return;}
                mcESDContainer      = (TList*) mcTopContainer->FindObject(Form("%s ESD histograms",mcCut[triggerSel+1].Data())); 
                if(mcESDContainer == NULL) {cout << "ERROR: " << Form("'%s' ESD histograms",mcCut[triggerSel+1].Data()) << " not found in File" << endl; return;}
            }
            if(mode==2||mode==3){
                dataInvMassPtAlpha      = (TH2F*) dataESDContainer->FindObject("ESD_Mother_InvMass_E_Calib");
                dataBGInvMassPtAlpha    = (TH2F*) dataESDContainer->FindObject("ESD_Background_InvMass_E_Calib");
                mcInvMassPtAlpha        = (TH2F*) mcESDContainer->FindObject("ESD_Mother_InvMass_E_Calib");
                mcBGInvMassPtAlpha      = (TH2F*) mcESDContainer->FindObject("ESD_Background_InvMass_E_Calib");
            }else{
                dataInvMassPtAlpha      = (TH2F*) dataESDContainer->FindObject("ESD_Mother_InvMass_vs_Pt_Alpha");
                dataBGInvMassPtAlpha    = (TH2F*) dataESDContainer->FindObject("ESD_Background_InvMass_vs_Pt_Alpha");
                mcInvMassPtAlpha        = (TH2F*) mcESDContainer->FindObject("ESD_Mother_InvMass_vs_Pt_Alpha");
                mcBGInvMassPtAlpha      = (TH2F*) mcESDContainer->FindObject("ESD_Background_InvMass_vs_Pt_Alpha");
            }
            if(!dataInvMassPtAlpha){cout << "did not find ESD_Mother_InvMass_E_Calib in triggered data" << endl; return;}
            if(!dataBGInvMassPtAlpha){cout << "did not find ESD_Background_InvMass_E_Calib in triggered data" << endl; return;}
            if(!mcInvMassPtAlpha){cout << "did not find ESD_Mother_InvMass_E_Calib in trigger mc" << endl; return;}
            if(!mcBGInvMassPtAlpha){cout << "did not find ESD_Background_InvMass_E_Calib in trigger mc" << endl; return;}

            triggerSel++;
            fOutput->cd();
        }

        cout << endl;
        cout << "-----------------------------------------------------" << endl;
        cout << "\t MC/Data Fitting mass positions" << endl;
        cout << "loop: " << iClusterPt << ", " << fBins[iClusterPt] << " - " << fBins[iClusterPt+1] << " GeV" << endl;
        cout << "-----------------------------------------------------" << endl;

        TH2* Hist2D;
        TH2* HistBG2D;
        for(Int_t iDataMC = 0; iDataMC < 2; iDataMC++)        {
            if(iDataMC==0) {
                Hist2D              = dataInvMassPtAlpha;
                HistBG2D            = dataBGInvMassPtAlpha;
            } else if(iDataMC==1) {
                Hist2D              = mcInvMassPtAlpha;
                HistBG2D            = mcBGInvMassPtAlpha;
            } else {
                cout << "ERROR: data/mc loop, returning..." << endl; return;
            }

            Double_t projectMin;
            Double_t projectMax;
            if(mode==2||mode==3){
              projectMin            = Hist2D->GetYaxis()->FindBin(fBins[iClusterPt]+0.001);
              projectMax            = Hist2D->GetYaxis()->FindBin(fBins[iClusterPt+1]-0.001);
            }else{
              projectMin            = Hist2D->GetYaxis()->FindBin((fBins[iClusterPt]*2)+0.001);
              projectMax            = Hist2D->GetYaxis()->FindBin((fBins[iClusterPt+1]*2)-0.001);
            }
            TH1D* sliceHist         = (TH1D*) Hist2D->ProjectionX(Form("slice%sAlpha_%f-%f",dataMC[iDataMC].Data(),fBins[iClusterPt],fBins[iClusterPt+1]),projectMin,projectMax);
            sliceHist->SetDirectory(0);
            sliceHist->SetTitle(Form("%s - %.01f < #it{E}_{Cluster} < %.01f (GeV)",dataMC[iDataMC].Data(),fBins[iClusterPt],fBins[iClusterPt+1]));
            sliceHist->GetYaxis()->SetTitle("#frac{d#it{M}_{inv}}{dN}");
            sliceHist->Sumw2();
            TH1D* sliceBGHist       = (TH1D*) HistBG2D->ProjectionX(Form("sliceBG%sAlpha_%f-%f",dataMC[iDataMC].Data(),fBins[iClusterPt],fBins[iClusterPt+1]),projectMin,projectMax);
            sliceBGHist->SetDirectory(0);
            sliceBGHist->SetTitle(Form("%s - %.01f < #it{E}_{Cluster} < %.01f (GeV)",dataMC[iDataMC].Data(),fBins[iClusterPt],fBins[iClusterPt+1]));
            sliceBGHist->GetYaxis()->SetTitle("#frac{d#it{M}_{inv}}{dN}");
            sliceBGHist->Sumw2();
            
    //*******************************************************************************
    // Rebin
             if( select.Contains("LHC11a-Pythia-ConvCalo") || select.Contains("LHC11a-Phojet-ConvCalo") || select.Contains("LHC11a-JetJet-ConvCalo") || select.Contains("LHC13g-Pythia-ConvCalo") ||
                 select.Contains("LHC13g-JetJet-ConvCalo") ) {
                if(fBins[iClusterPt]>=20){
                    sliceHist->Rebin(10);
                    sliceBGHist->Rebin(10);
                }else if(fBins[iClusterPt]>=12){
                    sliceHist->Rebin(8);
                    sliceBGHist->Rebin(8);
                } else if(fBins[iClusterPt]>=10){
                    sliceHist->Rebin(4);
                    sliceBGHist->Rebin(4);
//                else if(firstTriggerBin[0]>0 && fBins[iClusterPt]>=3.2 && fBins[iClusterPt]<fBins[firstTriggerBin[0]]) sliceHist->Rebin(4);
                } else {
                    sliceHist->Rebin(2);
                    sliceBGHist->Rebin(2);
                }    
            } else if(select.Contains("LHC11a-Pythia-Calo") || select.Contains("LHC11a-Phojet-Calo") || select.Contains("LHC11a-JetJet-Calo") || select.Contains("LHC13g-Pythia-Calo")  ||
                      select.Contains("LHC13g-JetJet-Calo")) {
                if(fBins[iClusterPt]>=4.5){ 
                    sliceHist->Rebin(8);
                    sliceBGHist->Rebin(8);
                } else if(fBins[iClusterPt]>=2.5){ 
                    sliceHist->Rebin(5);
                    sliceBGHist->Rebin(5);
                } else if(fBins[iClusterPt]>=1.5) {
                    sliceHist->Rebin(4);
                    sliceBGHist->Rebin(4);
                } else { 
                    sliceHist->Rebin(2);
                    sliceBGHist->Rebin(2);
                }    
            } else if(select.Contains("LHC10-Calo") ) {
               if(fBins[iClusterPt]>=5){
                   sliceHist->Rebin(4);
                   sliceBGHist->Rebin(4);
               } else if(fBins[iClusterPt]>=3.2){
                   sliceHist->Rebin(2);
                   sliceBGHist->Rebin(2);
               }    
            } else if(select.Contains("LHC10-ConvCalo") ) {
               if(fBins[iClusterPt]>=8){
                   sliceHist->Rebin(4);
                   sliceBGHist->Rebin(4);
               } else if(fBins[iClusterPt]>=5) {
                   sliceHist->Rebin(2);
                   sliceBGHist->Rebin(2);
               }    
            } else if(select.Contains("LHC12-Pythia-Calo") || select.Contains("LHC12-Phojet-Calo") ) {
               if(fBins[iClusterPt]>=3.6){
                 sliceHist->Rebin(4);
                 sliceBGHist->Rebin(4);
               } else if(fBins[iClusterPt]<=1.0){
                 sliceHist->Rebin(4);
                 sliceBGHist->Rebin(4);
               } else {
                 sliceHist->Rebin(2);
                 sliceBGHist->Rebin(2);
               }
            } else if(select.Contains("LHC12-Pythia-ConvCalo") || select.Contains("LHC12-Phojet-ConvCalo") ) {
               if(fBins[iClusterPt]>=10){
                 sliceHist->Rebin(10);
                 sliceBGHist->Rebin(10);
               } else if(fBins[iClusterPt]>=8){
                 sliceHist->Rebin(4);
                 sliceBGHist->Rebin(4);
               } else {
                 sliceHist->Rebin(2);
                 sliceBGHist->Rebin(2);
               }
            } else if(select.Contains("LHC13bc-Calo")) {
               if(fBins[iClusterPt]>=3.2){
                   sliceHist->Rebin(2);
                   sliceBGHist->Rebin(2);
               }   
            } else if(select.Contains("LHC13bc-ConvCalo")) {
               if(fBins[iClusterPt]>=8){
                   sliceHist->Rebin(2);
                   sliceBGHist->Rebin(2);
               }    
            }
    //             else if(fBins[iClusterPt]>=8) sliceHist->Rebin(5);
    //             else sliceHist->Rebin(4);
            //*******************************************************************************
            // Background subtraction ranges
//             if( ( ( select.Contains("LHC11a-Pythia-ConvCalo") || select.Contains("LHC11a-Phojet-ConvCalo") || select.Contains("LHC11a-JetJet-ConvCalo")) ) //&& fBins[iClusterPt]<2.8 
//                 || ( (select.Contains("LHC11a-Pythia-ConvCalo") || select.Contains("LHC11a-Phojet-ConvCalo") || select.Contains("LHC11a-JetJet-ConvCalo")) ) // && iDataMC==1 && fBins[iClusterPt]>=5
//                 || ( (select.Contains("LHC11a-Pythia-Calo") || select.Contains("LHC11a-Phojet-Calo"))  ) //&& fBins[iClusterPt]<3.6
//                 || ( select.Contains("LHC10-ConvCalo")  && fBins[iClusterPt]<5 )
//                 || ( select.Contains("LHC10-Calo")  && fBins[iClusterPt]<3.6 )
//                 || ( select.Contains("LHC13bc-ConvCalo")  && fBins[iClusterPt]<5 )
//                 || ( select.Contains("LHC13bc-Calo")  && fBins[iClusterPt]<2 )
//                 ){
                Double_t range              = 0.3;
                Double_t rangeMin           = 0.0;
                if(select.Contains("-ConvCalo")) 
                    rangeMin                = 0.0;
                if (select.Contains("JetJet"))
                    range                   = 0.2;
                
                Double_t integralSigAndBG   = sliceHist->Integral(sliceHist->FindBin(0.17), sliceHist->FindBin(0.3));
                Double_t integralBG         = sliceBGHist->Integral(sliceBGHist->FindBin(0.17), sliceBGHist->FindBin(0.3));
                cout << integralSigAndBG << "\t" << integralBG << "\t" << integralSigAndBG/ integralBG << endl;
                
                if (!isnan(integralSigAndBG/ integralBG) && !isinf(integralSigAndBG/ integralBG)){
                    sliceBGHist->Scale( integralSigAndBG/ integralBG );
                }    
                TH1D* sliceHistCopy         = (TH1D*)sliceHist->Clone("SliceCopy");
                if (!isnan(integralSigAndBG/ integralBG) && !isinf(integralSigAndBG/ integralBG)){
                    sliceHist->Add( sliceBGHist, -1);
                }                
//                 TF1* fBckFit                = FitBckg(sliceHist,rangeMin,range);
//                 fBckFit->SetLineColor(kRed+1);
//                 if (fBins[iClusterPt]<4){
//                     sliceHist->GetListOfFunctions()->Add(fBckFit);
//                 }
                sliceHistCopy->GetXaxis()->SetRangeUser(0.0,0.3);
                sliceHistCopy->Write(Form("slice%sAlpha_%f-%f-withWithBG",dataMC[iDataMC].Data(),fBins[iClusterPt],fBins[iClusterPt+1]));
                sliceBGHist->SetLineColor(kGreen+2);
                sliceBGHist->GetXaxis()->SetRangeUser(0.0,0.3);
                sliceBGHist->Write(Form("slice%sAlpha_%f-%f-BG",dataMC[iDataMC].Data(),fBins[iClusterPt],fBins[iClusterPt+1]));                
                sliceHist->SetLineColor(kRed+2);
                sliceHist->GetXaxis()->SetRangeUser(0.0,0.3);
                sliceHist->Write(Form("slice%sAlpha_%f-%f-withRemainingBckg",dataMC[iDataMC].Data(),fBins[iClusterPt],fBins[iClusterPt+1]));
                
                sliceHistCopy->GetYaxis()->SetRangeUser(sliceHist->GetMinimum(),sliceHistCopy->GetMaximum());
                sliceHistCopy->GetXaxis()->SetRangeUser(0.0,0.3);
                sliceHistCopy->DrawCopy();
                sliceBGHist->DrawCopy("same");
                sliceHist->DrawCopy("same");
                DrawGammaLines(0., 0.3, 0, 0, 1, kGray+1, 7);
                canvas->SetLogx(0); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
                if(select.Contains("-ConvCalo")) {
                    if(iClusterPt==exampleBin1 || iClusterPt==exampleBin2){
                        canvas->SaveAs(Form("%s/ConvCalo_slice%sAlpha_%.01f-%.01f-withBckgAndFit.%s",outputDirSample.Data(),dataMC[iDataMC].Data(),fBins[iClusterPt],fBins[iClusterPt+1],suffix.Data()));
                    }    
                    canvas->Write(Form("canvas_ConvCalo_slice%sAlpha_%.01f-%.01f-withBckgAndFit",dataMC[iDataMC].Data(),fBins[iClusterPt],fBins[iClusterPt+1]));
                } else {
                    if(iClusterPt==exampleBin1 || iClusterPt==exampleBin2){
                        canvas->SaveAs(Form("%s/Calo_slice%sAlpha_%.01f-%.01f-withBckgAndFit.%s",outputDirSample.Data(),dataMC[iDataMC].Data(),fBins[iClusterPt],fBins[iClusterPt+1],suffix.Data()));
                    }
                    canvas->Write(Form("canvas_Calo_slice%sAlpha_%.01f-%.01f-withBckgAndFit",dataMC[iDataMC].Data(),fBins[iClusterPt],fBins[iClusterPt+1]));
                }
                canvas->Clear();
                sliceHist->GetXaxis()->SetRangeUser(0.0,0.3);

//                 TF1* fBckg = new TF1(Form("fBckg%s",dataMC[iDataMC].Data()),"[0]+[1]*x",0.05,0.2); //+[2]*x*x
//                 fBckg->SetParameter(0,fBckFit->GetParameter(0));
//                 fBckg->SetParameter(1,fBckFit->GetParameter(1));
//                 fBckg->SetParameter(2,fBckFit->GetParameter(2));

//                 if (fBins[iClusterPt]<4){
//                     sliceHist->Add(fBckg,-1);
//                 }     
//             } else 
//                 sliceHist->GetXaxis()->SetRangeUser(0.05,0.2);

            Double_t sigmaRangeAdjust = 1.5;
            Double_t precision = 0.1;
            Double_t minMax[2]={0.04,0.3};
            
    //*******************************************************************************
    // Adjusting sigma fitting
            if(select.Contains("LHC13bc-")){
              sigmaRangeAdjust = 1;
              if(select.Contains("-Calo")){
//                if(fBins[iClusterPt]>=2.8) sigmaRangeAdjust = 1.5;
              }else if(select.Contains("-ConvCalo")){
                if(fBins[iClusterPt]>=5) sigmaRangeAdjust = 1.5;
              }
            }
    //*******************************************************************************
    // Fit
            if( mode == 4|| mode == 5){
                Double_t min = 0.02*fBins[iClusterPt] - 0.001;
                if (min > 0.04)
                    minMax[0]   = min;
                cout << minMax[0] << endl;
            }
//             fFitReco = FitRecursiveGaussian (sliceHist, precision, sigmaRangeAdjust, minMax[0], minMax[1]);
            fFitReco = FitExpPlusGaussian (sliceHist, minMax[0], minMax[1], mode);

            if(iDataMC==0) {
                histDataResults->SetBinContent(iClusterPt+1,fFitReco->GetParameter(1));
                histDataResults->SetBinError(iClusterPt+1,fFitReco->GetParError(1));
            }
            else if(iDataMC==1) {
                histMCResults->SetBinContent(iClusterPt+1,fFitReco->GetParameter(1));
                histMCResults->SetBinError(iClusterPt+1,fFitReco->GetParError(1));
            }

            sliceHist->GetListOfFunctions()->Add(fFitReco);
            sliceHist->GetXaxis()->SetRangeUser(0.0,0.3);
//             sliceHist->Write();
            
            sliceHist->DrawCopy();
            DrawGammaLines(0., 0.3, 0, 0, 1, kGray+1, 7);
            canvas->SetLogx(0); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
            if(select.Contains("-ConvCalo")) {
                if(iClusterPt==exampleBin1 || iClusterPt==exampleBin2){
                    canvas->SaveAs(Form("%s/ConvCalo_slice%sAlpha_%.01f-%.01f.%s",outputDirSample.Data(),dataMC[iDataMC].Data(),fBins[iClusterPt],fBins[iClusterPt+1],suffix.Data()));
                }    
                canvas->Write(Form("canvas_ConvCalo_slice%sAlpha_%.01f-%.01f",dataMC[iDataMC].Data(),fBins[iClusterPt],fBins[iClusterPt+1]));
            } else { 
                if(iClusterPt==exampleBin1 || iClusterPt==exampleBin2){
                    canvas->SaveAs(Form("%s/Calo_slice%sAlpha_%.01f-%.01f.%s",outputDirSample.Data(),dataMC[iDataMC].Data(),fBins[iClusterPt],fBins[iClusterPt+1],suffix.Data()));
                }    
                canvas->Write(Form("canvas_Calo_slice%sAlpha_%.01f-%.01f",dataMC[iDataMC].Data(),fBins[iClusterPt],fBins[iClusterPt+1]));
            }    
            canvas->Clear();
        }
    }
    delete canvas;

    cout << endl;
    cout << "-----------------------------------------------------" << endl;
    cout << "-----------------------------------------------------" << endl;

    histDataMCResults->Divide(histMCResults,histDataResults,1,1);
    DrawGammaSetMarker(histDataMCResults, 24, 2, kBlack, kBlack);

    Double_t minPlotY = 0.95;
    if(select.Contains("LHC10-Calo")) minPlotY = 0.9;


    //*********************************************************************************************************************************
    //************************************ Write mean mass for MC and data into output file *******************************************
    //*********************************************************************************************************************************    
    SetStyleHistoTH1ForGraphs(histMCResults, "#it{E}_{Cluster} (GeV)","#LT M_{#pi^{0} (MC)} #GT",0.035,0.043, 0.035,0.043, 1.,1.);
    DrawGammaSetMarker(histMCResults, markerStyle[1], 1, color[1], color[1]);
    histMCResults->Write("Mean mass MC");
    SetStyleHistoTH1ForGraphs(histDataResults, "#it{E}_{Cluster} (GeV)","#LT M_{#pi^{0} (data)} #GT",0.035,0.043, 0.035,0.043, 1.,1.);
    DrawGammaSetMarker(histDataResults, markerStyle[0], 1, color[0], color[0]);
    histDataResults->Write("Mean mass Data");
    
    //*********************************************************************************************************************************
    //*********************************** Plotting Mean mass for data and MC vs PDG value *********************************************
    //*********************************************************************************************************************************
    TCanvas *canvasMassPDG = new TCanvas("canvasMassPDG","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings(canvasMassPDG, 0.1, 0.02, 0.06, 0.1);
    canvasMassPDG->SetLogx(1); 
    canvasMassPDG->SetLogy(0); 
  
    TH2F * histoDummyMeanMassVsPDG;
    histoDummyMeanMassVsPDG = new TH2F("histoDummyMeanMassVsPDG","histoDummyMeanMassVsPDG",11000,0.5,fBins[endPtBin]*1.5,1000,0.89,1.1);
    SetStyleHistoTH2ForGraphs(histoDummyMeanMassVsPDG, "#it{E}_{Cluster} (GeV)","#LT M_{#pi^{0} (MC/data)} #GT / M_{#pi^{0} (PDG)}",0.035,0.043, 0.035,0.043, 1.,1.);
    histoDummyMeanMassVsPDG->GetXaxis()->SetMoreLogLabels();
    histoDummyMeanMassVsPDG->GetXaxis()->SetLabelOffset(-0.01);
    histoDummyMeanMassVsPDG->DrawCopy("");
    
    Double_t rangeExponent[2]   = {-0.5, -0.08};
    Double_t rangeMult[2]       = {-0.2, -0.001};
    
    TLegend *legend = GetAndSetLegend2(0.15, 0.95, 0.95, 0.99, 0.043, 2, "", 42);
   
    TH1D* histDataResultsVsPDG =  (TH1D*)histDataResults->Clone("Mean mass data / mass PDG Pi0");
    histDataResultsVsPDG->Scale(1/massPi0);
    SetStyleHistoTH1ForGraphs(histDataResultsVsPDG, "#it{E}_{Cluster} (GeV)","#LT M_{#pi^{0} (data)} #GT / M_{#pi^{0} (PDG)}",0.035,0.043, 0.035,0.043, 1.,1.);
    DrawGammaSetMarker(histDataResultsVsPDG, markerStyle[0], 1, color[0], color[0]);
    TF1* fitMassDataVsPDG = new TF1("fitMassDataVsPDG", "[0] + [1]*pow(x,[2])" ,fBins[startPtBin],fBins[endPtBin]);
    fitMassDataVsPDG->SetParLimits(1, rangeMult[0], rangeMult[1]);
    fitMassDataVsPDG->SetParLimits(2, rangeExponent[0], rangeExponent[1]);
    
    histDataResultsVsPDG->Fit(fitMassDataVsPDG,"QRME0");
    fitMassDataVsPDG->SetLineColor(kBlack);
    fitMassDataVsPDG->Draw("same");
    cout << WriteParameterToFile(fitMassDataVsPDG) << endl;
    
    histDataResultsVsPDG->Write();
    histDataResultsVsPDG->DrawCopy("same");
    legend->AddEntry(histDataResultsVsPDG,"Data");

    TH1D* histMCResultsVsPDG =  (TH1D*)histMCResults->Clone("Mean mass MC / mass PDG Pi0");
    histMCResultsVsPDG->Scale(1/massPi0);
    SetStyleHistoTH1ForGraphs(histMCResultsVsPDG, "#it{E}_{Cluster} (GeV)","#LT M_{#pi^{0} (MC)} #GT / M_{#pi^{0} (PDG)}",0.035,0.043, 0.035,0.043, 1.,1.);
    DrawGammaSetMarker(histMCResultsVsPDG, markerStyle[1], 1, color[1], color[1]);
    TF1* fitMassMCVsPDG = new TF1("fitMassMCVsPDG", "[0] + [1]*pow(x,[2])" ,fBins[startPtBin],fBins[endPtBin]);
    fitMassMCVsPDG->SetParLimits(1, rangeMult[0], rangeMult[1]);
    fitMassMCVsPDG->SetParLimits(2, rangeExponent[0], rangeExponent[1]);

    histMCResultsVsPDG->Fit(fitMassMCVsPDG,"QRME0");
    fitMassMCVsPDG->SetLineColor(kRed+2);
    fitMassMCVsPDG->Draw("same");
    cout << WriteParameterToFile(fitMassMCVsPDG) << endl;
    
    histMCResultsVsPDG->Write();
    histMCResultsVsPDG->DrawCopy("same");
    legend->AddEntry(histMCResultsVsPDG,"MC");
    
    PutProcessLabelAndEnergyOnPlot(0.7, 0.89, 0.03, fCollisionSystem.Data(), fTextMeasurement.Data(), recGamma.Data());
    PutProcessLabelAndEnergyOnPlot(0.2, 0.89, 0.03, fPlot[0].Data(),"", "");
    for (Int_t i = 0; i < numberOfTriggers; i++){
       PutProcessLabelAndEnergyOnPlot(0.2, 0.89-2*0.03*(i+1), 0.03, fPlot[i+1].Data(),"", "");
    }
        
    legend->Draw("same");
    canvasMassPDG->Update();
    canvasMassPDG->SaveAs(Form("%s/MeanMass_Pi0_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
    canvasMassPDG->Clear();
    delete canvasMassPDG;

    //*********************************************************************************************************************************
    //****************************** Fitting ratio of mean mass position in MC/data ***************************************************
    //*********************************************************************************************************************************
    
    Double_t startFit               = 3;
    if(select.Contains("LHC10-ConvCalo")) 
        startFit                    = 3;
    if(select.Contains("-Calo")) 
        startFit                    = 3;
    if(select.Contains("LHC13bc-Calo")) 
        startFit                    = 2.4;
    if(select.Contains("LHC13bc-ConvCalo")) 
        startFit                    = 3;
    if(select.Contains("LHC12-Phojet-ConvCalo"))
        startFit                    = 6;
    if(select.Contains("LHC12-Pythia-ConvCalo"))
        startFit                    = 8;

    TF1* fFitConst = new TF1("DataMCConst", "[0]" ,startFit,fBins[endPtBin]);
    histDataMCResults->Fit(fFitConst,"QRME0");
    Double_t highPtConst            = fixedOffSet;
    if (highPtConst == -1){
        highPtConst = fFitConst->GetParameter(0);
    }
    
    fFitMassPos = FitDataMC(histDataMCResults, fBins[startPtBin], fBins[endPtBin], select,  highPtConst);
    TF1* fFitComposit = new TF1("fFitComposit", "([0] + [1]*pow(x,[2]))/([3] + [4]*pow(x,[5]))" ,fBins[startPtBin],fBins[endPtBin]);
    fFitComposit->SetParameter(0, fitMassMCVsPDG->GetParameter(0) );
    fFitComposit->SetParameter(1, fitMassMCVsPDG->GetParameter(1) );
    fFitComposit->SetParameter(2, fitMassMCVsPDG->GetParameter(2) );
    fFitComposit->SetParameter(3, fitMassDataVsPDG->GetParameter(0) );
    fFitComposit->SetParameter(4, fitMassDataVsPDG->GetParameter(1) );
    fFitComposit->SetParameter(5, fitMassDataVsPDG->GetParameter(2) );

    
    TF1* fFitCompositFitted = new TF1("fFitCompositFitted", "([0] + [1]*pow(x,[2]))/([3] + [4]*pow(x,[5]))" ,fBins[startPtBin],fBins[endPtBin-2]);
    fFitCompositFitted->SetParameter(0, fitMassMCVsPDG->GetParameter(0) );
    fFitCompositFitted->SetParameter(1, fitMassMCVsPDG->GetParameter(1) );
    fFitCompositFitted->FixParameter(2, fitMassMCVsPDG->GetParameter(2) );
    fFitCompositFitted->SetParameter(3, fitMassDataVsPDG->GetParameter(0) );
    fFitCompositFitted->SetParameter(4, fitMassDataVsPDG->GetParameter(1) );
    fFitCompositFitted->FixParameter(5, fitMassDataVsPDG->GetParameter(2) );
    histDataMCResults->Fit(fFitCompositFitted,"QRME0");
    
    TF1* fFitCompositInverted = new TF1("fFitCompositInverted", "([0] + [1]*pow(x,[2]))/([3] + [4]*pow(x,[5]))" ,fBins[startPtBin],fBins[endPtBin]);
    fFitCompositInverted->SetParameter(0, fitMassDataVsPDG->GetParameter(0) );
    fFitCompositInverted->SetParameter(1, fitMassDataVsPDG->GetParameter(1) );
    fFitCompositInverted->SetParameter(2, fitMassDataVsPDG->GetParameter(2) );
    fFitCompositInverted->SetParameter(3, fitMassMCVsPDG->GetParameter(0) );
    fFitCompositInverted->SetParameter(4, fitMassMCVsPDG->GetParameter(1) );
    fFitCompositInverted->SetParameter(5, fitMassMCVsPDG->GetParameter(2) );
    
    TF1* fFitCompositInvertedFitted = new TF1("fFitCompositInvertedFitted", "([0] + [1]*pow(x,[2]))/([3] + [4]*pow(x,[5]))" ,fBins[startPtBin],fBins[endPtBin]);
    fFitCompositInvertedFitted->SetParameter(0, fFitCompositFitted->GetParameter(3) );
    fFitCompositInvertedFitted->SetParameter(1, fFitCompositFitted->GetParameter(4) );
    fFitCompositInvertedFitted->SetParameter(2, fFitCompositFitted->GetParameter(5) );
    fFitCompositInvertedFitted->SetParameter(3, fFitCompositFitted->GetParameter(0) );
    fFitCompositInvertedFitted->SetParameter(4, fFitCompositFitted->GetParameter(1) );
    fFitCompositInvertedFitted->SetParameter(5, fFitCompositFitted->GetParameter(2) );
    
    histDataMCResults->GetYaxis()->SetRangeUser(minPlotY,1.05);
    histDataMCResults->GetXaxis()->SetRangeUser(fBins[startPtBin],fBins[endPtBin]);

    histDataMCResults->Write("MeanMassRatioMCData-noFit");
    fFitMassPos->Write("MeanMassRatioMCData-Fit");
    histDataMCResults->GetListOfFunctions()->Add(fFitMassPos);

    fstream fLog;
    fLog.open(Form("%s/CorrectCaloNonLinearity3_%s.log",outputDir.Data(),select.Data()), ios::out);
    fLog << "FitDataMC results:" << endl;
//     if (select.Contains("JetJet")){
//         fLog << "Par " << 0 << ": " << fFitMassPos->GetParameter(0) << " +- " << fFitMassPos->GetParError(0) << endl;
//     } else {    
        for(Int_t i=0;i<=2;i++) fLog << "Par " << i << ": " << fFitMassPos->GetParameter(i) << " +- " << fFitMassPos->GetParError(i) << endl;
//     }    
    
    fLog << WriteParameterToFile(fFitComposit) << endl;
    fLog << WriteParameterToFile(fFitCompositFitted) << endl;
    fLog << WriteParameterToFile(fFitCompositInverted) << endl;
    fLog.close();
    
    

    //*******************************************************************************
    // plotting mass ratios
    //*******************************************************************************
    TCanvas *canvasMassRatioMCData = new TCanvas("canvasMassPDG","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings(canvasMassRatioMCData, 0.1, 0.02, 0.02, 0.1);
    canvasMassRatioMCData->SetLogx(1); 
    canvasMassRatioMCData->SetLogy(0); 
    
    SetStyleHistoTH1ForGraphs(histDataMCResults, "#it{E}_{Cluster} (GeV)","#LT M_{#pi^{0} (MC)} #GT / #LT M_{#pi^{0} (data)} #GT",0.035,0.043, 0.035,0.043, 1.,0.9);
    DrawGammaSetMarker(histDataMCResults, markerStyle[0], 1, color[0], color[0]);
    histDataMCResults->Draw();
    fFitComposit->SetLineColor(kGreen+2);
    fFitComposit->Draw("same");
//     fFitCompositFitted->SetLineColor(kBlue+2);
//     fFitCompositFitted->Draw("same");

    fFitMassPos->Draw("same");
    PutProcessLabelAndEnergyOnPlot(0.7, 0.89, 0.03, fCollisionSystem.Data(), fTextMeasurement.Data(), recGamma.Data());
    PutProcessLabelAndEnergyOnPlot(0.2, 0.89, 0.03, fPlot[0].Data(),"", "");
    for (Int_t i = 0; i < numberOfTriggers; i++){
       PutProcessLabelAndEnergyOnPlot(0.2, 0.89-2*0.03*(i+1), 0.03, fPlot[i+1].Data(),"", "");
    }

    canvasMassRatioMCData->Update();
    canvasMassRatioMCData->SaveAs(Form("%s/MeanMassRatio_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
    canvasMassRatioMCData->Clear();
    //*******************************************************************************
    // plotting total correction
    //*******************************************************************************
    canvasMassRatioMCData->cd();
    TH1D* totalCorrection = new TH1D("Total Correction","; #it{E}_{Cluster} (GeV); correction factor",1000,0.5,50);

    
    SetStyleHistoTH1ForGraphs(totalCorrection, "#it{E}_{Cluster} (GeV)","correction factor",0.035,0.043, 0.035,0.043, 1.,0.9);
    DrawGammaSetMarker(totalCorrection, 8, 1, kBlack, kBlack);
    totalCorrection->GetYaxis()->SetRangeUser(minPlotY,1.1);
    SetLogBinningXTH(totalCorrection);
 
    for(Int_t iBin = 1; iBin < totalCorrection->GetNbinsX()+1; iBin++){
        Float_t e = totalCorrection->GetXaxis()->GetBinCenter(iBin);
        Float_t factor = 1;
        Float_t p0 = fFitMassPos->GetParameter(0);
//        if(select.Contains("-ConvCalo")) p0*= 0.995*0.9970;
//        else if(select.Contains("-Calo"))p0*= 0.995*0.9981;
        factor /= FunctionNL_kSDM(e,p0,fFitMassPos->GetParameter(1),fFitMassPos->GetParameter(2));
        totalCorrection->SetBinContent(iBin,factor);
    }
    totalCorrection->DrawCopy("p");
    fFitCompositInverted->SetRange(0.5,50);
    fFitCompositInverted->SetLineColor(kGreen+2);
    fFitCompositInverted->Draw("same");
//     fFitCompositInvertedFitted->SetRange(0.5,50);
//     fFitCompositInvertedFitted->SetLineColor(kBlue+2);
//     fFitCompositInvertedFitted->Draw("same");
    
    
    TLegend *legend2 = GetAndSetLegend2(0.2, 0.2, 0.4, 0.29, 0.03, 1, "", 42);
    legend2->AddEntry(totalCorrection,"Correction factor for MC, direct ratio fit","l");
    legend2->AddEntry(fFitCompositInverted,"Correction factor for MC from mass fits","l");
//     legend2->AddEntry(fFitCompositInvertedFitted,"Correction factor for MC from mass fits, refit ratio","l");
    legend2->Draw("same");

    PutProcessLabelAndEnergyOnPlot(0.7, 0.89, 0.03, fCollisionSystem.Data(), fTextMeasurement.Data(), recGamma.Data());
    PutProcessLabelAndEnergyOnPlot(0.2, 0.89, 0.03, fPlot[0].Data(),"", "");
    for (Int_t i = 0; i < numberOfTriggers; i++){
       PutProcessLabelAndEnergyOnPlot(0.2, 0.89-2*0.03*(i+1), 0.03, fPlot[i+1].Data(),"", "");
    }   

    canvasMassRatioMCData->Update();
    canvasMassRatioMCData->SaveAs(Form("%s/TotalCorrection_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
    canvasMassRatioMCData->Clear();
    delete canvasMassRatioMCData;

    cout << "-----------------------------------------------------" << endl;
    cout << "-----------------------------------------------------" << endl;

    fOutput->Write();
    fOutput->Close();
    return;
}

//*******************************************************************************
//*******************************************************************************
//*******************************************************************************

Float_t FunctionNL_kSDM(Float_t e, Float_t p0, Float_t p1, Float_t p2){
    return ( p0 + exp( p1 + ( p2 * e ) ) );
}

Double_t fitExcludeSignal(Double_t *x, Double_t *par)
{
    if (x[0] > 0.08 && x[0] < 0.17) {
      TF1::RejectPoint();
      return 0;
   }
   return par[0] + par[1]*x[0]; //+ par[2]*x[0]*x[0]
}


TF1* FitBckg(TH1* fHisto, Double_t minFit, Double_t maxFit){
    TF1* fFitBckg = new TF1("fFitBckg",fitExcludeSignal,minFit,maxFit,3);
    fFitBckg->SetLineColor(kBlue);
    fFitBckg->SetLineWidth(2);
    fFitBckg->SetLineStyle(1);
    fHisto->Fit(fFitBckg,"QRME0");
    return fFitBckg;
}

TF1* FitDataMC(TH1* fHisto, Double_t minFit, Double_t maxFit, TString selection, Double_t constPar){

//     if (!selection.Contains("JetJet")){
        cout << "running standard fit from " <<  minFit << "\t"<<  maxFit << endl;
        TF1* fFitReco = new TF1("DataMC", "[0]+exp([1]+([2]*x))" ,minFit,maxFit);

        fFitReco->SetParameter(0,1.);
        fFitReco->SetParameter(1,-1.);
        fFitReco->SetParameter(2,-0.5);
        if(constPar!=-1) fFitReco->FixParameter(0,constPar);
//         if(constPar!=-1) fFitReco->FixParameter(1,-2.95303);
//         if(constPar!=-1) fFitReco->FixParameter(2,-1.17524);
        fHisto->Fit(fFitReco,"QRME0");

        fFitReco->SetLineColor(kRed);
        fFitReco->SetLineWidth(2);
        fFitReco->SetLineStyle(1);

        if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("PROBLEMS") == 0){
            cout << "Parameters for DataMC: " << endl;
            for(Int_t i=0;i<=2;i++) cout << "Par " << i << ": " << fFitReco->GetParameter(i) << " +- " << fFitReco->GetParError(i) << endl;
        } else {
            cout << "DataMC fitting failed with status " << gMinuit->fCstatu.Data() <<endl << endl;
        }

        return fFitReco;
//     } else {
//         cout << "running Jet-Jet MC fit" << endl;
//         TF1* fFitReco = new TF1("DataMC", "[0]" ,minFit,maxFit);
// 
//         fFitReco->SetParameter(0,1.);
//         fHisto->Fit(fFitReco,"QRME0");
// 
//         fFitReco->SetLineColor(kRed);
//         fFitReco->SetLineWidth(2);
//         fFitReco->SetLineStyle(1);
// 
//         if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("PROBLEMS") == 0){
//             cout << "Parameters for DataMC: " << endl;
//             for(Int_t i=0;i<=0;i++) cout << "Par " << i << ": " << fFitReco->GetParameter(i) << " +- " << fFitReco->GetParError(i) << endl;
//         } else {
//             cout << "DataMC fitting failed with status " << gMinuit->fCstatu.Data() <<endl << endl;
//         }
// 
//         return fFitReco;
//         
//     }
}

TF1* FitExpPlusGaussian(TH1D* histo, Double_t fitRangeMin, Double_t fitRangeMax, Int_t mode ){

    Double_t mesonAmplitude =histo->GetMaximum();
    Double_t mesonAmplitudeMin;
    Double_t mesonAmplitudeMax;
    
    mesonAmplitudeMin = mesonAmplitude*98./100.;
    mesonAmplitudeMax = mesonAmplitude*400./100.;

    TF1* fFitReco    = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)", 
                               fitRangeMin, fitRangeMax);
    Double_t fMesonMassExpect = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    
    fFitReco->SetParameter(0, mesonAmplitude);
    fFitReco->SetParameter(1, fMesonMassExpect);
    fFitReco->SetParameter(2, 0.01);
    fFitReco->SetParameter(3, 0.012);

    fFitReco->SetParLimits(0, mesonAmplitudeMin, mesonAmplitudeMax);
    if (mode == 4 || mode == 5){
        fFitReco->SetParLimits(1, fMesonMassExpect*0.8, fMesonMassExpect*1.2);
    } else {
        fFitReco->SetParLimits(1, fMesonMassExpect*0.9, fMesonMassExpect*1.1);
    }    
    fFitReco->SetParLimits(2, 0.001, 0.05);
    fFitReco->SetParLimits(3, 0.001, 0.09);
    
    histo->Fit(fFitReco,"QRME0");
    
    fFitReco->SetLineColor(kRed+1);
    fFitReco->SetLineWidth(1);
    fFitReco->SetLineStyle(1);

    if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("PROBLEMS") == 0){
        cout << "Parameter for exponential+Gaussian "<< endl;
        cout << gMinuit->fCstatu.Data() << endl;
        cout << "Gausexp: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<<endl;
    } else {
        cout << "Exp+Gaussian fitting failed in with status " << gMinuit->fCstatu.Data() <<endl << endl;
    }
    return fFitReco;
}


TF1* FitRecursiveGaussian (TH1* histo, Double_t precision, Double_t correctRange, Double_t fitRangeMin, Double_t fitRangeMax ){
    TF1 *f0             = new TF1("f0", "gaus", fitRangeMin,fitRangeMax);
    histo->Fit(f0,"0RMEQ");
    Double_t rp         = f0->GetParameter(2);
    Double_t mp         = f0->GetParameter(1);
    Double_t ymin       = mp -(rp * correctRange);
    Double_t ymax       = mp + (rp * correctRange);
    Double_t deviation  = 100;
    Int_t counter       = 0;
    TF1* f1             = new TF1 ("f1", "gaus", ymin, ymax);
    while(deviation > precision && counter < 100){
        f1->SetRange(ymin,ymax);
        histo->Fit(f1,"0RMEQ");
        Double_t rp2    = f1->GetParameter(2);
        if (rp2>rp){ 
            deviation   = rp2-rp;
        } else { 
            deviation   = rp-rp2 ;
        }
        rp              = rp2 ;
        mp              = f1->GetParameter(1);
        ymin            = mp -(rp * correctRange);
        ymax            = mp +(rp * correctRange);
        counter++;
    }
    delete f0;

    if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("PROBLEMS") == 0){
        cout << "Parameters for FitRecursiveGaussian: " << endl;
        for(Int_t i=0;i<=2;i++) cout << "Par " << i << ": " << f1->GetParameter(i) << " +- " << f1->GetParError(i) << endl;
    } else {
        cout << "FitRecursiveGaussian fitting failed with status " << gMinuit->fCstatu.Data() <<endl << endl;
    }

    return f1;
}
