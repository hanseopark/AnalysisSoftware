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

// -----------------------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------- Main function -----------------------------------------------------------------------------
// ----------------------- This macro is used to compile the pp external input file for data, i.e. charged hadron, pions, kaons ... --------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------------------------------
void PrepareChargedPionDataALICE_pp(){ 
    
    // *********************************************************************************************************************
    // ************************************** global variable definition ***************************************************
    // *********************************************************************************************************************
    TString dateForOutput                               = ReturnDateStringForOutput();

    
    // -----------------------------------------------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------- Read 0.9 TeV spectra -----------------------------------------------------------------------------------------
    // -----------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // *********************************************************************************************************************
    // ********************************** CMS charged pion results 0.9 TeV *************************************************
    // *********************************************************************************************************************
    Double_t chPionCMS900GeV_xval[22]               = { 0.125, 0.175, 0.225, 0.275, 0.325,
                                                        0.375, 0.425, 0.475, 0.525, 0.575,
                                                        0.625, 0.675, 0.725, 0.775, 0.825, 
                                                        0.875, 0.925, 0.975, 1.025, 1.075,
                                                        1.125, 1.175 };
    Double_t ptbins900GeV[23]                       = { 0.1, 0.15, 0.2, 0.25, 0.3,
                                                        0.35, 0.4, 0.45, 0.5, 0.55,
                                                        0.6, 0.65, 0.7, 0.75, 0.8, 
                                                        0.85, 0.9, 0.95, 1.0, 1.05,
                                                        1.1, 1.15, 1.2};
    Double_t chPionPlCMS900GeV_xerrminus[22]        = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
                                                        0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
                                                        0.025, 0.025, 0.025 };
    Double_t chPionPlCMS900GeV_yval[22]             = { 3.094, 3.861, 3.843, 3.494, 3.023,
                                                        2.584, 2.187, 1.836, 1.541, 1.293,
                                                        1.085, 0.9125, 0.77, 0.6521, 0.5514,
                                                        0.4708, 0.4042, 0.347, 0.3008, 0.2591,
                                                        0.2201, 0.1933 };
    Double_t chPionPlCMS900GeV_yErrMin[22]          = { 0.269, 0.124, 0.092, 0.076, 0.063,
                                                        0.053, 0.044, 0.037, 0.031, 0.026,
                                                        0.022, 0.0188, 0.0163, 0.0145, 0.0126,
                                                        0.0111, 0.0101, 0.0100, 0.0099, 0.0095,
                                                        0.0094, 0.0089};
    Double_t chPionPlCMS900GeV_yStatMin[22]         = { 0.002, 0.002, 0.002, 0.002, 0.002,
                                                        0.001, 0.001, 0.001, 0.001, 0.001,
                                                        0.001, 0.001, 9.0E-4, 9.0E-4, 9.0E-4,
                                                        9.0E-4, 9.0E-4, 9.0E-4, 0.001, 0.001,
                                                        0.0012, 0.0017 };
    Double_t chPionMinCMS900GeV_yval[22]            = { 2.952, 3.704, 3.768, 3.43, 2.989, 
                                                        2.546, 2.157, 1.808, 1.519, 1.271,
                                                        1.068, 0.8984, 0.7608, 0.6457, 0.5527,
                                                        0.4681, 0.3991, 0.3386, 0.2913, 0.2646,
                                                        0.2335, 0.2064};
    Double_t chPionMinCMS900GeV_yErrMin[22]         = { 0.271, 0.124, 0.092, 0.073, 0.061,
                                                        0.052, 0.044, 0.037, 0.031, 0.026,
                                                        0.022, 0.0186, 0.0161, 0.0142, 0.0125,
                                                        0.0108, 0.0097, 0.0093, 0.0086, 0.0087,
                                                        0.0091, 0.0081};
    Double_t chPionMinCMS900GeV_yStatMin[22]        = { 0.002, 0.002, 0.002, 0.002, 0.002, 
                                                        0.001, 0.001, 0.001, 0.001, 0.001, 
                                                        0.001, 0.001, 9.0E-4, 9.0E-4, 9.0E-4,
                                                        9.0E-4, 9.0E-4, 9.0E-4, 9.0E-4, 0.001,
                                                        0.0011, 0.0016 };

    TH1D* histoChargedPionPlusSpecLowPtStat900GeVCMS    = new TH1D("histoChargedPionPlusSpecLowPtStat900GeVCMS","histoChargedPionPlusSpecLowPtStat900GeVCMS", 22, ptbins900GeV);
    TH1D* histoChargedPionPlusSpecLowPtSys900GeVCMS     = new TH1D("histoChargedPionPlusSpecLowPtSys900GeVCMS","histoChargedPionPlusSpecLowPtSys900GeVCMS", 22, ptbins900GeV);
    TH1D* histoChargedPionMinusSpecLowPtStat900GeVCMS   = new TH1D("histoChargedPionMinusSpecLowPtStat900GeVCMS","histoChargedPionMinusSpecLowPtStat900GeVCMS", 22, ptbins900GeV);
    TH1D* histoChargedPionMinusSpecLowPtSys900GeVCMS    = new TH1D("histoChargedPionMinusSpecLowPtSys900GeVCMS","histoChargedPionMinusSpecLowPtSys900GeVCMS", 22, ptbins900GeV);
    
    for(Int_t i=1; i<23; i++){
        histoChargedPionPlusSpecLowPtStat900GeVCMS->SetBinContent(i, chPionPlCMS900GeV_yval[i-1]);
        histoChargedPionPlusSpecLowPtStat900GeVCMS->SetBinError(i, chPionPlCMS900GeV_yStatMin[i-1]);
        histoChargedPionPlusSpecLowPtSys900GeVCMS->SetBinContent(i, chPionPlCMS900GeV_yval[i-1]);
        histoChargedPionPlusSpecLowPtSys900GeVCMS->SetBinError(i, chPionPlCMS900GeV_yErrMin[i-1]);
        
        histoChargedPionMinusSpecLowPtStat900GeVCMS->SetBinContent(i, chPionMinCMS900GeV_yval[i-1]);
        histoChargedPionMinusSpecLowPtStat900GeVCMS->SetBinError(i, chPionMinCMS900GeV_yStatMin[i-1]);
        histoChargedPionMinusSpecLowPtSys900GeVCMS->SetBinContent(i, chPionMinCMS900GeV_yval[i-1]);
        histoChargedPionMinusSpecLowPtSys900GeVCMS->SetBinError(i, chPionMinCMS900GeV_yErrMin[i-1]);
    }
    TH1D*   histoChargedPionSpecLowPtSys900GeVCMS       = (TH1D*)histoChargedPionMinusSpecLowPtSys900GeVCMS->Clone("histoChargedPionSpecLowPtSys900GeVCMS");
    histoChargedPionSpecLowPtSys900GeVCMS->Add(histoChargedPionPlusSpecLowPtSys900GeVCMS);
    histoChargedPionSpecLowPtSys900GeVCMS->Scale(0.5);
    for (Int_t i = 0; i < histoChargedPionSpecLowPtSys900GeVCMS->GetNbinsX(); i++){
        Double_t fractionalSystematicError              =  0;
        if (histoChargedPionPlusSpecLowPtSys900GeVCMS->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSys900GeVCMS->GetBinContent(i) !=0){
            fractionalSystematicError                   = ( histoChargedPionPlusSpecLowPtSys900GeVCMS->GetBinError(i)/histoChargedPionPlusSpecLowPtSys900GeVCMS->GetBinContent(i)*100 +
                                                            histoChargedPionMinusSpecLowPtSys900GeVCMS->GetBinError(i)/histoChargedPionMinusSpecLowPtSys900GeVCMS->GetBinContent(i)*100 )/2;
        }
        histoChargedPionSpecLowPtSys900GeVCMS->SetBinError(i, histoChargedPionSpecLowPtSys900GeVCMS->GetBinContent(i)*fractionalSystematicError/100.);
        
    }   
    TH1D*   histoChargedPionSpecLowPtStat900GeVCMS      = (TH1D*)histoChargedPionMinusSpecLowPtStat900GeVCMS->Clone("histoChargedPionSpecLowPtStat900GeVCMS");
    histoChargedPionSpecLowPtStat900GeVCMS->Add(histoChargedPionPlusSpecLowPtStat900GeVCMS);
    histoChargedPionSpecLowPtStat900GeVCMS->Scale(0.5);

    for (Int_t i = 1; i < histoChargedPionSpecLowPtSys900GeVCMS->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtSys900GeVCMS->SetBinContent(i, histoChargedPionSpecLowPtSys900GeVCMS->GetBinContent(i)/histoChargedPionSpecLowPtSys900GeVCMS->GetBinCenter(i)/(2*TMath::Pi())*0.78);
        histoChargedPionSpecLowPtSys900GeVCMS->SetBinError(i, histoChargedPionSpecLowPtSys900GeVCMS->GetBinError(i)/histoChargedPionSpecLowPtSys900GeVCMS->GetBinCenter(i)/(2*TMath::Pi())*0.78);
        histoChargedPionSpecLowPtStat900GeVCMS->SetBinContent(i, histoChargedPionSpecLowPtStat900GeVCMS->GetBinContent(i)/histoChargedPionSpecLowPtStat900GeVCMS->GetBinCenter(i)/(2*TMath::Pi())*0.78);
        histoChargedPionSpecLowPtStat900GeVCMS->SetBinError(i, i, histoChargedPionSpecLowPtStat900GeVCMS->GetBinError(i)/histoChargedPionSpecLowPtStat900GeVCMS->GetBinCenter(i)/(2*TMath::Pi())*0.78);
    }
    
    
    
    // *********************************************************************************************************************
    // *********************************** ALICE charged pion results 0.9 TeV **********************************************
    // *********************************** EPJ C71,1655, arXiv: 1101.4110 **************************************************
    // *********************************************************************************************************************
    Double_t chPionALICE900GeV_xval[33]         = { 0.11, 0.13, 0.15, 0.17, 0.19, 0.225, 0.275, 0.325, 0.375, 
                                                    0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 
                                                    0.925, 0.975, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 
                                                    1.9, 2.1, 2.3, 2.5 };
    Double_t ptbinsALICE900GeV[34]              = { 0.1, 0.12, 0.14, 0.16, 0.18,
                                                    0.20, 0.25, 0.30, 0.35, 0.40,
                                                    0.45, 0.50, 0.55, 0.60, 0.65,
                                                    0.70, 0.75, 0.80, 0.85, 0.9,
                                                    0.95, 1., 1.1, 1.2, 1.3, 
                                                    1.4, 1.5, 1.6, 1.7, 1.8, 
                                                    2., 2.2, 2.4, 2.6};

    Double_t p8038_d1x1y1_xerrminus[33]         = { 0.01, 0.01, 0.01, 0.01, 0.01, 0.025, 0.025, 0.025, 0.025, 
                                                    0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
                                                    0.025, 0.025, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 
                                                    0.10, 0.10, 0.10, 0.10 };
    Double_t chPionPlALICE900GeV_yval[33]       = { 2.658354, 2.888623, 3.064501, 2.9925, 3.115569,
                                                    2.977894, 2.673851, 2.362263, 2.005345, 1.663289,
                                                    1.412473, 1.235371, 1.00173, 0.82576, 0.697456, 
                                                    0.616081, 0.504638, 0.40602, 0.350737, 0.299044, 
                                                    0.258118, 0.206977, 0.143769, 0.109963, 0.083156,
                                                    0.065611, 0.047197, 0.037742, 0.030838, 0.021651, 
                                                    0.012858, 0.008694, 0.006367 };
    Double_t chPionPlALICE900GeV_yErrMin[33]    = { 0.368345, 0.09928, 0.093292, 0.069068, 0.070011, 
                                                    0.059972, 0.055397, 0.048109, 0.040615, 0.033297,
                                                    0.029966, 0.020852, 0.017897, 0.014441, 0.023033,
                                                    0.020305, 0.0167, 0.013316, 0.011571, 0.00985, 
                                                    0.008538, 0.007001, 0.004977, 0.003897, 0.00304, 
                                                    0.002433, 0.001784, 0.001482, 0.001231,8.88E-4,
                                                    5.57E-4, 3.85E-4, 2.93E-4};
    Double_t chPionPlALICE900GeV_yStatMin[33]   = { 0.034355, 0.031949, 0.03253, 0.029677, 0.030266, 0.020327, 0.01895, 0.017561, 0.015952, 
                                                    0.014239, 0.013051, 0.01227, 0.010865, 0.009732, 0.010077, 0.009482, 0.008571, 0.007629, 0.007099, 
                                                    0.006588, 0.006106, 0.00385, 0.003184, 0.002797, 0.002437, 0.002176, 0.001843, 0.001662, 0.001514, 
                                                    8.97E-4, 6.8E-4, 5.69E-4, 4.9E-4 };
   
    Double_t chPionMinALICE900GeV_yval[33]      = { 2.664472, 2.790213, 2.976358, 3.066082, 3.026797, 2.950636, 2.668285, 2.343535, 1.998761, 
                                                    1.716659, 1.436463, 1.221818, 1.024817, 0.81446, 0.69469, 0.595425, 0.513817, 0.410521, 0.346321, 
                                                    0.308068, 0.256686, 0.198254, 0.143069, 0.110292, 0.083025, 0.060515, 0.046562, 0.036817, 0.028274, 
                                                    0.020431, 0.012429, 0.008739, 0.005992 };
    Double_t chPionMinALICE900GeV_yErrMin[33]   = { 0.361721, 0.094854, 0.087912, 0.070138, 0.072816, 
                                                    0.059373, 0.056005, 0.047364, 0.041263, 0.033914,
                                                    0.029768, 0.021121, 0.018372, 0.014214, 0.023628,
                                                    0.020205, 0.017292, 0.013744, 0.011479, 0.010148, 
                                                    0.008487, 0.006695, 0.004908, 0.003895, 0.002985, 
                                                    0.002261, 0.001767, 0.00144, 0.001115,8.36E-4,
                                                    5.29E-4, 3.85E-4, 2.81E-4};
    Double_t chPionMinALICE900GeV_yStatMin[33]  = { 0.034388, 0.031276, 0.031371, 0.029964, 0.029549, 0.020477, 0.018907, 0.01754, 0.015626, 
                                                    0.014267, 0.01317, 0.012069, 0.011089, 0.009505, 0.010239, 0.009452, 0.008735, 0.007768, 0.007097, 
                                                    0.006716, 0.006117, 0.00378, 0.003204, 0.002812, 0.002439, 0.002088, 0.001829, 0.001643, 0.001451, 
                                                    8.66E-4, 6.74E-4, 5.64E-4, 4.75E-4 };
  
    TH1D* histoChargedPionPlusSpecLowPtStat900GeVALICE  = new TH1D("histoChargedPionPlusSpecLowPtStat900GeVALICE","histoChargedPionPlusSpecLowPtStat900GeVALICE", 33, ptbinsALICE900GeV);
    TH1D* histoChargedPionPlusSpecLowPtSys900GeVALICE   = new TH1D("histoChargedPionPlusSpecLowPtSys900GeVALICE","histoChargedPionPlusSpecLowPtSys900GeVALICE", 33, ptbinsALICE900GeV);
    TH1D* histoChargedPionMinusSpecLowPtStat900GeVALICE = new TH1D("histoChargedPionMinusSpecLowPtStat900GeVALICE","histoChargedPionMinusSpecLowPtStat900GeVALICE", 33, ptbinsALICE900GeV);
    TH1D* histoChargedPionMinusSpecLowPtSys900GeVALICE  = new TH1D("histoChargedPionMinusSpecLowPtSys900GeVALICE","histoChargedPionMinusSpecLowPtSys900GeVALICE", 33, ptbinsALICE900GeV);
    
    for(Int_t i=1; i<34; i++){
        histoChargedPionPlusSpecLowPtStat900GeVALICE->SetBinContent(i, chPionPlALICE900GeV_yval[i-1]);
        histoChargedPionPlusSpecLowPtStat900GeVALICE->SetBinError(i, chPionPlALICE900GeV_yStatMin[i-1]);
        histoChargedPionPlusSpecLowPtSys900GeVALICE->SetBinContent(i, chPionPlALICE900GeV_yval[i-1]);
        histoChargedPionPlusSpecLowPtSys900GeVALICE->SetBinError(i, chPionPlALICE900GeV_yErrMin[i-1]);
        
        histoChargedPionMinusSpecLowPtStat900GeVALICE->SetBinContent(i, chPionMinALICE900GeV_yval[i-1]);
        histoChargedPionMinusSpecLowPtStat900GeVALICE->SetBinError(i, chPionMinALICE900GeV_yStatMin[i-1]);
        histoChargedPionMinusSpecLowPtSys900GeVALICE->SetBinContent(i, chPionMinALICE900GeV_yval[i-1]);
        histoChargedPionMinusSpecLowPtSys900GeVALICE->SetBinError(i, chPionMinALICE900GeV_yErrMin[i-1]);
    }
    TH1D* histoChargedPionSpecLowPtSys900GeVALICE       = (TH1D*)histoChargedPionMinusSpecLowPtSys900GeVALICE->Clone("histoChargedPionSpecLowPtSys900GeVALICE");
    histoChargedPionSpecLowPtSys900GeVALICE->Add(histoChargedPionPlusSpecLowPtSys900GeVALICE);
    histoChargedPionSpecLowPtSys900GeVALICE->Scale(0.5);
    TH1D* histoChargedPionSpecLowPtStat900GeVALICE      = (TH1D*)histoChargedPionMinusSpecLowPtStat900GeVALICE->Clone("histoChargedPionSpecLowPtStat900GeVALICE");
    histoChargedPionSpecLowPtStat900GeVALICE->Add(histoChargedPionPlusSpecLowPtStat900GeVALICE);
    histoChargedPionSpecLowPtStat900GeVALICE->Scale(0.5);

    for (Int_t i = 1; i < histoChargedPionSpecLowPtSys900GeVALICE->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtSys900GeVALICE->SetBinContent(i, histoChargedPionSpecLowPtSys900GeVALICE->GetBinContent(i)/histoChargedPionSpecLowPtSys900GeVALICE->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtSys900GeVALICE->SetBinError(i, histoChargedPionSpecLowPtSys900GeVALICE->GetBinError(i)/histoChargedPionSpecLowPtSys900GeVALICE->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStat900GeVALICE->SetBinContent(i, histoChargedPionSpecLowPtStat900GeVALICE->GetBinContent(i)/histoChargedPionSpecLowPtStat900GeVALICE->GetBinCenter(i)/(2*TMath::Pi()));
        Double_t fractionalSystematicError              = 0;
        if (histoChargedPionPlusSpecLowPtSys900GeVALICE->GetBinContent(i) !=  0 && histoChargedPionMinusSpecLowPtSys900GeVALICE->GetBinContent(i)){
            fractionalSystematicError                   = ( histoChargedPionPlusSpecLowPtSys900GeVALICE->GetBinError(i)/histoChargedPionPlusSpecLowPtSys900GeVALICE->GetBinContent(i)*100 +
                                                            histoChargedPionMinusSpecLowPtSys900GeVALICE->GetBinError(i)/histoChargedPionMinusSpecLowPtSys900GeVALICE->GetBinContent(i)*100)/2;
        }
        histoChargedPionSpecLowPtStat900GeVALICE->SetBinError(i, histoChargedPionSpecLowPtStat900GeVALICE->GetBinContent(i)*fractionalSystematicError/100.);
    }

    // *********************************************************************************************************************
    // *********************************** ALICE charged kaon results 0.9 TeV **********************************************
    // *********************************** EPJ C71,1655, arXiv: 1101.4110 **************************************************
    // *********************************************************************************************************************
    
    Double_t chKaonPlALICE_PP900GeV_xval[27]        = { 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 
                                                        0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.05, 1.15, 1.25, 
                                                        1.35, 1.45, 1.55, 1.65, 1.75, 1.9, 2.1, 2.3 };
    Double_t chKaonPlALICE_PP900GeV_xerrminus[27]   = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
                                                        0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.050, 0.050, 0.050, 0.050, 
                                                        0.050, 0.050, 0.050, 0.050, 0.10, 0.10, 0.10 };
    Double_t chKaonPlALICE_PP900GeV_xerrplus[27]    = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
                                                        0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.050, 0.050, 0.050, 0.050, 
                                                        0.050, 0.050, 0.050, 0.050, 0.10, 0.10, 0.10 };
    Double_t chKaonPlALICE_PP900GeV_yval[27]        = { 0.202248, 0.230913, 0.230026, 0.240136, 0.220438, 0.202544, 0.191518, 0.176627, 0.158315, 0.137851, 
                                                        0.130541, 0.11659, 0.098213, 0.091789, 0.081461, 0.06399, 0.05873, 0.045613, 0.030931, 0.029217, 
                                                        0.023534, 0.020722, 0.016497, 0.013606, 0.007975, 0.007493, 0.005139 };
    Double_t chKaonPlALICE_PP900GeV_yerrminus[27]   = { 0.03061919822594968, 0.018455291951090885, 0.013897246058122452, 0.010877166772648104, 0.013138047724072248, 
                                                        0.012436832152923832, 0.010419838482433401, 0.012180917042653233, 0.010781991281762382, 0.00928566771966346, 
                                                        0.008996097153766182, 0.008143350231937713, 0.00720326377415127, 0.006806095209442782, 0.0062881502049489885, 
                                                        0.0052131625717984275, 0.004326051432888888, 0.003618112491341307, 0.002717267929373178, 0.0026539383941606483, 
                                                        0.0023332903805570364, 0.0021739829806141538, 0.0018572291188757514, 0.0016157812351924377, 9.630160954002794E-4, 
                                                        9.526352922288781E-4, 7.318825042313828E-4 };
    Double_t chKaonPlALICE_PP900GeV_yerrplus[27]    = { 0.03061919822594968, 0.018455291951090885, 0.013897246058122452, 0.010877166772648104, 0.013138047724072248, 
                                                        0.012436832152923832, 0.010419838482433401, 0.012180917042653233, 0.010781991281762382, 0.00928566771966346, 
                                                        0.008996097153766182, 0.008143350231937713, 0.00720326377415127, 0.006806095209442782, 0.0062881502049489885, 
                                                        0.0052131625717984275, 0.004326051432888888, 0.003618112491341307, 0.002717267929373178, 0.0026539383941606483, 
                                                        0.0023332903805570364, 0.0021739829806141538, 0.0018572291188757514, 0.0016157812351924377, 9.630160954002794E-4, 
                                                        9.526352922288781E-4, 7.318825042313828E-4 };
    Double_t chKaonPlALICE_PP900GeV_ystatminus[27]  = { 0.00772, 0.008315, 0.007338, 0.010034, 0.009283, 0.005965, 0.007197, 0.007676, 0.006894, 0.006116, 
                                                        0.005808, 0.005363, 0.00482, 0.004534, 0.004192, 0.003642, 0.002436, 0.002113, 0.001696, 0.001642, 
                                                        0.00147, 0.001381, 0.001234, 0.001093, 5.94E-4, 5.65E-4, 4.64E-4 };
    Double_t chKaonPlALICE_PP900GeV_ystatplus[27]   = { 0.00772, 0.008315, 0.007338, 0.010034, 0.009283, 0.005965, 0.007197, 0.007676, 0.006894, 0.006116, 
                                                        0.005808, 0.005363, 0.00482, 0.004534, 0.004192, 0.003642, 0.002436, 0.002113, 0.001696, 0.001642, 
                                                        0.00147, 0.001381, 0.001234, 0.001093, 5.94E-4, 5.65E-4, 4.64E-4 };
    Int_t chKaonPlALICE_PP900GeV_numpoints          = 27;
    
    Double_t chKaonPlALICE_PP900GeV_ysysminus[27];
    Double_t chKaonPlALICE_PP900GeV_ysysplus[27];
    ExtractSystematicFromTotal(chKaonPlALICE_PP900GeV_numpoints, chKaonPlALICE_PP900GeV_yerrminus, chKaonPlALICE_PP900GeV_ystatminus, chKaonPlALICE_PP900GeV_ysysminus );
    ExtractSystematicFromTotal(chKaonPlALICE_PP900GeV_numpoints, chKaonPlALICE_PP900GeV_yerrplus, chKaonPlALICE_PP900GeV_ystatplus, chKaonPlALICE_PP900GeV_ysysplus );
    
    TGraphAsymmErrors* graphChKaonPlALICEStatPP900GeV   = new TGraphAsymmErrors(chKaonPlALICE_PP900GeV_numpoints, chKaonPlALICE_PP900GeV_xval, chKaonPlALICE_PP900GeV_yval, chKaonPlALICE_PP900GeV_xerrminus,
                                                                                chKaonPlALICE_PP900GeV_xerrplus, chKaonPlALICE_PP900GeV_ystatminus, chKaonPlALICE_PP900GeV_ystatplus);
    TGraphAsymmErrors* graphChKaonPlALICESysPP900GeV    = new TGraphAsymmErrors(chKaonPlALICE_PP900GeV_numpoints, chKaonPlALICE_PP900GeV_xval, chKaonPlALICE_PP900GeV_yval, chKaonPlALICE_PP900GeV_xerrminus,
                                                                                chKaonPlALICE_PP900GeV_xerrplus, chKaonPlALICE_PP900GeV_ysysminus, chKaonPlALICE_PP900GeV_ysysplus);

    Double_t chKaonMinALICE_PP900GeV_xval[27]       = { 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 
                                                        0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.05, 1.15, 1.25, 1.35, 
                                                        1.45, 1.55, 1.65, 1.75, 1.9, 2.1, 2.3 };
    Double_t chKaonMinALICE_PP900GeV_xerrminus[27]  = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
                                                        0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.050, 0.050, 0.050, 0.050,
                                                        0.050, 0.050, 0.050, 0.050, 0.10, 0.10, 0.10 };
    Double_t chKaonMinALICE_PP900GeV_xerrplus[27]   = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
                                                        0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.050, 0.050, 0.050, 0.050,
                                                        0.050, 0.050, 0.050, 0.050, 0.10, 0.10, 0.10 };
    Double_t chKaonMinALICE_PP900GeV_yval[27]       = { 0.209256, 0.22594, 0.23524, 0.237882, 0.222334, 0.215391, 0.198406, 0.16838, 0.161305, 0.133114, 
                                                        0.129337, 0.113505, 0.103056, 0.084654, 0.077756, 0.06853, 0.055854, 0.047125, 0.037168, 0.026843, 
                                                        0.020708, 0.019498, 0.014419, 0.010567, 0.009063, 0.00633, 0.004643 };
    Double_t chKaonMinALICE_PP900GeV_yerrminus[27]  = { 0.03142716062580264, 0.018171193466583312, 0.014159715533865785, 0.01091617648263347, 0.012911892386478442, 
                                                        0.013362359671854368, 0.010234025845189175, 0.012562617919844572, 0.011587632760836011, 0.009759357458357595, 
                                                        0.009312138798364209, 0.00838503911738043, 0.007716341425831285, 0.006793307883498289, 0.006281135327311457,
                                                        0.005720049387898675, 0.004268888028515154, 0.0037681970489877514, 0.003210732782403419, 0.00253426320653558, 
                                                        0.0021020413887457116, 0.0020598681996671533, 0.0016787471518963182, 0.0013416948237211024, 0.001065, 
                                                        8.342835249482036E-4, 6.900275356824538E-4 };
    Double_t chKaonMinALICE_PP900GeV_yerrplus[27]   = { 0.03142716062580264, 0.018171193466583312, 0.014159715533865785, 0.01091617648263347, 0.012911892386478442, 
                                                        0.013362359671854368, 0.010234025845189175, 0.012562617919844572, 0.011587632760836011, 0.009759357458357595, 
                                                        0.009312138798364209, 0.00838503911738043, 0.007716341425831285, 0.006793307883498289, 0.006281135327311457, 
                                                        0.005720049387898675, 0.004268888028515154, 0.0037681970489877514, 0.003210732782403419, 0.00253426320653558, 
                                                        0.0021020413887457116, 0.0020598681996671533, 0.0016787471518963182, 0.0013416948237211024, 0.001065, 
                                                        8.342835249482036E-4, 6.900275356824538E-4 };
    Double_t chKaonMinALICE_PP900GeV_ystatminus[27] = { 0.008045, 0.008144, 0.007238, 0.010085, 0.009026, 0.0063, 0.006122, 0.007845, 0.007228, 0.006253, 
                                                        0.005948, 0.005441, 0.005055, 0.004526, 0.004206, 0.003881, 0.002434, 0.002178, 0.001903, 0.001593, 
                                                        0.001383, 0.001344, 0.001136, 9.69E-4, 6.39E-4, 5.23E-4, 4.57E-4 };
    Double_t chKaonMinALICE_PP900GeV_ystatplus[27]  = { 0.008045, 0.008144, 0.007238, 0.010085, 0.009026, 0.0063, 0.006122, 0.007845, 0.007228, 0.006253, 
                                                        0.005948, 0.005441, 0.005055, 0.004526, 0.004206, 0.003881, 0.002434, 0.002178, 0.001903, 0.001593, 
                                                        0.001383, 0.001344, 0.001136, 9.69E-4, 6.39E-4, 5.23E-4, 4.57E-4 };    
    Int_t chKaonMinALICE_PP900GeV_numpoints          = 27;
    
    Double_t chKaonMinALICE_PP900GeV_ysysminus[27];
    Double_t chKaonMinALICE_PP900GeV_ysysplus[27];
    ExtractSystematicFromTotal(chKaonMinALICE_PP900GeV_numpoints, chKaonMinALICE_PP900GeV_yerrminus, chKaonMinALICE_PP900GeV_ystatminus, chKaonMinALICE_PP900GeV_ysysminus );
    ExtractSystematicFromTotal(chKaonMinALICE_PP900GeV_numpoints, chKaonMinALICE_PP900GeV_yerrplus, chKaonMinALICE_PP900GeV_ystatplus, chKaonMinALICE_PP900GeV_ysysplus );
    
    TGraphAsymmErrors* graphChKaonMinALICEStatPP900GeV  = new TGraphAsymmErrors(chKaonMinALICE_PP900GeV_numpoints, chKaonMinALICE_PP900GeV_xval, chKaonMinALICE_PP900GeV_yval, chKaonMinALICE_PP900GeV_xerrminus,
                                                                                chKaonMinALICE_PP900GeV_xerrplus, chKaonMinALICE_PP900GeV_ystatminus, chKaonMinALICE_PP900GeV_ystatplus);
    TGraphAsymmErrors* graphChKaonMinALICESysPP900GeV   = new TGraphAsymmErrors(chKaonMinALICE_PP900GeV_numpoints, chKaonMinALICE_PP900GeV_xval, chKaonMinALICE_PP900GeV_yval, chKaonMinALICE_PP900GeV_xerrminus,
                                                                                chKaonMinALICE_PP900GeV_xerrplus, chKaonMinALICE_PP900GeV_ysysminus, chKaonMinALICE_PP900GeV_ysysplus);

    // *********************************************************************************************************************
    // *********************************** ALICE charged hadron results 0.9 TeV ********************************************
    // *********************************** arXiv: Eur.Phys.J. C73 (2013) no.12, 2662 arXiv:1307.1093 ***********************
    // *********************************************************************************************************************
    
    Double_t chHadALICE_PP900GeV_xval[54]       = { 0.175, 0.225, 0.275, 0.325, 0.375, 
                                                    0.425, 0.475, 0.525, 0.575, 0.625, 
                                                    0.675, 0.725, 0.775, 0.825, 0.875, 
                                                    0.925, 0.975, 1.05, 1.15, 1.25, 
                                                    1.35, 1.45, 1.55, 1.65, 1.75,
                                                    1.85, 1.95, 2.1, 2.3, 2.5, 
                                                    2.7, 2.9, 3.1, 3.3, 3.5, 
                                                    3.7, 3.9, 4.25, 4.75, 5.25, 
                                                    5.75, 6.25, 6.75, 7.5, 8.5,
                                                    9.5, 10.5, 11.5, 12.5, 13.5, 
                                                    14.5, 15.5, 17.0, 19.0 };
    Double_t chHadALICE_PP900GeV_xerrminus[54]  = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
                                                    0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.050, 0.050, 
                                                    0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.100, 0.100, 
                                                    0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.25, 0.25, 
                                                    0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
                                                    0.5, 0.5, 0.5, 1.0, 1.0 };
    Double_t chHadALICE_PP900GeV_xerrplus[54]   = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
                                                    0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.050, 0.050, 
                                                    0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.100, 0.100, 
                                                    0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.25, 0.25, 
                                                    0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
                                                    0.5, 0.5, 0.5, 1.0, 1.0 };
    Double_t chHadALICE_PP900GeV_yval[54]       = { 256.2, 198.9, 154.9, 119.4, 91.68, 70.9, 54.94, 42.71, 33.56, 
                                                    26.51, 21.06, 16.9, 13.54, 11.0, 9.002, 7.431, 6.098, 4.615, 3.201, 
                                                    2.27, 1.625, 1.162, 0.8488, 0.6375, 0.4721, 0.3581, 0.2728, 0.1874, 0.114, 
                                                    0.07148, 0.04502, 0.02896, 0.01986, 0.01325, 0.008893, 0.006293, 0.004496, 0.0026, 0.001199, 
                                                    6.158E-4, 3.394E-4, 1.857E-4, 1.025E-4, 5.206E-5, 2.229E-5, 9.721E-6, 5.3E-6, 3.236E-6, 1.825E-6, 
                                                    5.102E-7, 6.276E-7, 4.19E-7, 1.066E-7, 8.192E-8 };
    Double_t chHadALICE_PP900GeV_yerrminus[54]  = { 20.90215299915298, 13.201515064567399, 10.300485425454474, 7.900632886041472, 6.100663898298282, 
                                                    4.7205190392582885, 3.6504931173746926, 2.8404401067440235, 2.2303587155433093, 
                                                    1.7604544867732308, 1.4003213916812096, 1.1204017136723774, 0.9002221947941519, 0.7302739212103907,
                                                    0.5992703897240377, 0.4942590413942875, 0.40624130759931343, 0.30710421683851885, 0.21308449028495716, 
                                                    0.15108275877809485, 0.10807404868885037, 0.07705841939723394, 0.05656933798445939, 0.04246787020795839, 
                                                    0.03146362979695763, 0.02386063704095094, 0.01815406290613757, 0.012525573839150046, 0.007623647421018368, 
                                                    0.004772221704824704, 0.003009269014229203, 0.0019502051174171398, 0.0013366001645967278, 8.944271909999159E-4,
                                                    6.040744987168388E-4, 4.3123659399452635E-4, 3.1030630029053554E-4, 1.7756407294269863E-4, 8.381527307120105E-5, 
                                                    4.430846420267802E-5, 2.5495097567963924E-5, 1.4977316181479243E-5, 9.068627239003708E-6, 4.476293109259044E-6, 
                                                    2.291920591992663E-6, 1.2701409370617105E-6, 8.405242411733288E-7, 6.143590155601202E-7, 4.273839023641391E-7, 
                                                    2.1184466006958967E-7, 2.2685486549774508E-7, 1.7396700836652906E-7, 6.203104061677507E-8, 4.765967687678967E-8 };
    Double_t chHadALICE_PP900GeV_yerrplus[54]   = { 20.90215299915298, 13.201515064567399, 10.300485425454474, 7.900632886041472, 6.100663898298282, 
                                                    4.7205190392582885, 3.6504931173746926, 2.8404401067440235, 2.2303587155433093, 
                                                    1.7604544867732308, 1.4003213916812096, 1.1204017136723774, 0.9002221947941519, 0.7302739212103907, 
                                                    0.5992703897240377, 0.4942590413942875, 0.40624130759931343, 0.30710421683851885, 0.21308449028495716, 
                                                    0.15108275877809485, 0.10807404868885037, 0.07705841939723394, 0.05656933798445939, 0.04246787020795839, 
                                                    0.03146362979695763, 0.02386063704095094, 0.01815406290613757, 0.012525573839150046, 0.007623647421018368, 
                                                    0.004772221704824704, 0.003009269014229203, 0.0019502051174171398, 0.0013366001645967278, 8.944271909999159E-4,
                                                    6.040744987168388E-4, 4.3123659399452635E-4, 3.1030630029053554E-4, 1.7756407294269863E-4, 8.381527307120105E-5, 
                                                    4.430846420267802E-5, 2.5495097567963924E-5, 1.4977316181479243E-5, 9.068627239003708E-6, 4.476293109259044E-6,
                                                    2.291920591992663E-6, 1.2701409370617105E-6, 8.405242411733288E-7, 6.143590155601202E-7, 4.273839023641391E-7, 
                                                    2.1184466006958967E-7, 2.2685486549774508E-7, 1.7396700836652906E-7, 6.203104061677507E-8, 4.765967687678967E-8 };
    Double_t chHadALICE_PP900GeV_ystatminus[54] = { 0.3, 0.2, 0.1, 0.1, 0.09, 0.07, 0.06, 0.05, 0.04, 
                                                    0.04, 0.03, 0.03, 0.02, 0.02, 0.018, 0.016, 0.014, 0.008, 0.006, 
                                                    0.005, 0.004, 0.003, 0.0028, 0.0024, 0.002, 0.0017, 0.0014, 8.0E-4, 6.0E-4, 
                                                    4.6E-4, 3.4E-4, 2.8E-4, 2.1E-4, 1.6E-4, 1.25E-4, 1.02E-4, 8.3E-5, 4.0E-5, 2.5E-5, 
                                                    1.68E-5, 1.18E-5, 8.4E-6, 6.0E-6, 2.84E-6, 1.75E-6, 1.093E-6, 7.6E-7, 5.74E-7, 4.09E-7, 
                                                    2.09E-7, 2.228E-7, 1.716E-7, 6.16E-8, 4.732E-8 };
    Double_t chHadALICE_PP900GeV_ystatplus[54]  = { 0.3, 0.2, 0.1, 0.1, 0.09, 0.07, 0.06, 0.05, 0.04, 
                                                    0.04, 0.03, 0.03, 0.02, 0.02, 0.018, 0.016, 0.014, 0.008, 0.006, 
                                                    0.005, 0.004, 0.003, 0.0028, 0.0024, 0.002, 0.0017, 0.0014, 8.0E-4, 6.0E-4, 
                                                    4.6E-4, 3.4E-4, 2.8E-4, 2.1E-4, 1.6E-4, 1.25E-4, 1.02E-4, 8.3E-5, 4.0E-5, 2.5E-5, 
                                                    1.68E-5, 1.18E-5, 8.4E-6, 6.0E-6, 2.84E-6, 1.75E-6, 1.093E-6, 7.6E-7, 5.74E-7, 4.09E-7, 
                                                    2.09E-7, 2.228E-7, 1.716E-7, 6.16E-8, 4.732E-8 };
    Int_t chHadALICE_PP900GeV_numpoints         = 54;
    
    Double_t chHadALICE_PP900GeV_ysysminus[54];
    Double_t chHadALICE_PP900GeV_ysysplus[54];
    
    ExtractSystematicFromTotal(chHadALICE_PP900GeV_numpoints, chHadALICE_PP900GeV_yerrminus, chHadALICE_PP900GeV_ystatminus, chHadALICE_PP900GeV_ysysminus );
    ExtractSystematicFromTotal(chHadALICE_PP900GeV_numpoints, chHadALICE_PP900GeV_yerrplus, chHadALICE_PP900GeV_ystatplus, chHadALICE_PP900GeV_ysysplus );
    
    TGraphAsymmErrors* graphChHadALICEStatPP900GeV      = new TGraphAsymmErrors(chHadALICE_PP900GeV_numpoints, chHadALICE_PP900GeV_xval, chHadALICE_PP900GeV_yval, chHadALICE_PP900GeV_xerrminus,
                                                                                chHadALICE_PP900GeV_xerrplus, chHadALICE_PP900GeV_ystatminus, chHadALICE_PP900GeV_ystatplus);
    TGraphAsymmErrors* graphChHadALICESysPP900GeV       = new TGraphAsymmErrors(chHadALICE_PP900GeV_numpoints, chHadALICE_PP900GeV_xval, chHadALICE_PP900GeV_yval, chHadALICE_PP900GeV_xerrminus,
                                                                                chHadALICE_PP900GeV_xerrplus, chHadALICE_PP900GeV_ysysminus, chHadALICE_PP900GeV_ysysplus);
    
    // *********************************************************************************************************************
    // *********************************** CMS charged hadron results 0.9 TeV **********************************************
    // *********************************** JHEP 08,086 arXiv:1104.3547 *****************************************************
    // *********************************************************************************************************************
    
    Double_t chHadCMS_PP900GeV_xval[20]         = { 0.5, 0.7, 0.9, 1.1, 1.4, 1.8, 2.2, 2.8, 3.6, 4.4,
                                                    5.2, 6.0, 6.80, 8.2, 10.2, 12.2, 15.2, 19.2, 23.2, 31.2 };
    Double_t chHadCMS_PP900GeV_xerrminus[20]    = { 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.4, 0.4, 0.4, 
                                                    0.4, 0.4, 0.4, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 6.0 };
    Double_t chHadCMS_PP900GeV_xerrplus[20]     = { 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.4, 0.4, 0.4,
                                                    0.4, 0.4, 0.4, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 6.0 };
    Double_t chHadCMS_PP900GeV_yval[20]         = { 1.12039, 0.43415, 0.184944, 0.085811, 0.030526, 0.0091459, 0.0031714, 7.9086E-4, 1.62026E-4, 4.2808E-5, 
                                                    1.36076E-5, 5.098E-6, 2.09996E-6, 5.8109E-7, 1.19391E-7, 3.3481E-8, 7.9256E-9, 1.5183E-9, 4.632E-10, 4.601E-11 };
    Double_t chHadCMS_PP900GeV_yerrminus[20]    = { 0.048140502697832314, 0.018660324220120077, 0.007952249555943276, 0.0036912060359725247, 0.0013130746361117481, 
                                                    3.9385870562931574E-4, 1.366474295404052E-4, 3.4108066494599194E-5, 7.005212130977905E-6, 1.8608871540208989E-6, 
                                                    5.981088947674996E-7, 2.294484691602888E-7, 9.909224843548562E-8, 2.7742505294223158E-8, 6.880535008267889E-9, 
                                                    2.606266678603707E-9, 7.46213588726445E-10, 2.6744616280664786E-10, 1.3501959117105931E-10, 1.899983420980299E-11 };
    Double_t chHadCMS_PP900GeV_yerrplus[20]     = { 0.048140502697832314, 0.018660324220120077, 0.007952249555943276, 0.0036912060359725247, 0.0013130746361117481,
                                                    3.9385870562931574E-4, 1.366474295404052E-4, 3.4108066494599194E-5, 7.005212130977905E-6, 1.8608871540208989E-6, 
                                                    5.981088947674996E-7, 2.294484691602888E-7, 9.909224843548562E-8, 2.7742505294223158E-8, 6.880535008267889E-9, 
                                                    2.606266678603707E-9, 7.46213588726445E-10, 2.6744616280664786E-10, 1.3501959117105931E-10, 1.899983420980299E-11 };
    Double_t chHadCMS_PP900GeV_ystatminus[20]   = { 2.2E-4, 1.1E-4, 6.3E-5, 3.9E-5, 1.4E-5, 6.8E-6, 3.6E-6, 1.11E-6, 4.46E-7, 2.1E-7, 
                                                    1.095E-7, 6.38E-8, 3.959E-8, 1.171E-8, 4.539E-9, 2.165E-9, 6.624E-10, 2.592E-10, 1.335E-10, 1.889E-11 };
    Double_t chHadCMS_PP900GeV_ystatplus[20]    = { 2.2E-4, 1.1E-4, 6.3E-5, 3.9E-5, 1.4E-5, 6.8E-6, 3.6E-6, 1.11E-6, 4.46E-7, 2.1E-7, 
                                                    1.095E-7, 6.38E-8, 3.959E-8, 1.171E-8, 4.539E-9, 2.165E-9, 6.624E-10, 2.592E-10, 1.335E-10, 1.889E-11 };
    Int_t chHadCMS_PP900GeV_numpoints           = 20;
    Double_t chHadCMS_PP900GeV_ysysminus[20];
    Double_t chHadCMS_PP900GeV_ysysplus[20];
    
    ExtractSystematicFromTotal(chHadCMS_PP900GeV_numpoints, chHadCMS_PP900GeV_yerrminus, chHadCMS_PP900GeV_ystatminus, chHadCMS_PP900GeV_ysysminus );
    ExtractSystematicFromTotal(chHadCMS_PP900GeV_numpoints, chHadCMS_PP900GeV_yerrplus, chHadCMS_PP900GeV_ystatplus, chHadCMS_PP900GeV_ysysplus );
    
    TGraphAsymmErrors* graphChHadCMSStatPP900GeV    = new TGraphAsymmErrors(chHadCMS_PP900GeV_numpoints, chHadCMS_PP900GeV_xval, chHadCMS_PP900GeV_yval, chHadCMS_PP900GeV_xerrminus,
                                                                            chHadCMS_PP900GeV_xerrplus, chHadCMS_PP900GeV_ystatminus, chHadCMS_PP900GeV_ystatplus);
    TGraphAsymmErrors* graphChHadCMSSysPP900GeV     = new TGraphAsymmErrors(chHadCMS_PP900GeV_numpoints, chHadCMS_PP900GeV_xval, chHadCMS_PP900GeV_yval, chHadCMS_PP900GeV_xerrminus,
                                                                            chHadCMS_PP900GeV_xerrplus, chHadCMS_PP900GeV_ysysminus, chHadCMS_PP900GeV_ysysplus);
  
    // -----------------------------------------------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------- Read 7TeV spectra------------------------------------------------------------------------------------
    // -----------------------------------------------------------------------------------------------------------------------------------------------------------------

    // *********************************************************************************************************************    
    // ********************************** CMS charged pion results 7 TeV ***************************************************
    // *********************************************************************************************************************

    Double_t chPionCMS7TeV_xval[22]         = { 0.125, 0.175, 0.225, 0.275, 0.325,
                                                0.375, 0.425, 0.475, 0.525, 0.575,
                                                0.625, 0.675, 0.725, 0.775, 0.825, 
                                                0.875, 0.925, 0.975, 1.025, 1.075,
                                                1.125, 1.175 };
    Double_t ptbins7TeV[23]                 = { 0.1, 0.15, 0.2, 0.25, 0.3,
                                                0.35, 0.4, 0.45, 0.5, 0.55,
                                                0.6, 0.65, 0.7, 0.75, 0.8, 
                                                0.85, 0.9, 0.95, 1.0, 1.05,
                                                1.1, 1.15, 1.2};
    Double_t chPionPlCMS7TeV_xerrminus[22]  = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
                                                0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
                                                0.025, 0.025, 0.025 };
    Double_t chPionPlCMS7TeV_yval[22]       = {4.705, 5.801, 5.673, 5.107, 4.44,
                                                3.818, 3.273, 2.798, 2.392, 2.053, 
                                                1.753, 1.506, 1.304, 1.136, 0.9809,
                                                0.8585, 0.748, 0.6617, 0.5889, 0.5137,
                                                0.4518, 0.3356 };
    Double_t chPionPlCMS7TeV_yErrMin[22]    = { 0.413, 0.179, 0.134, 0.111, 0.094, 
                                                0.078, 0.067, 0.057, 0.048, 0.042, 
                                                0.036, 0.031, 0.028, 0.026, 0.0229,
                                                0.0206, 0.0191, 0.0200, 0.196, 0.0193, 
                                                0.0206, 0.0212};
    Double_t chPionPlCMS7TeV_yStatMin[22]   = { 0.003, 0.003, 0.003, 0.003, 0.002,
                                                0.002, 0.002, 0.002, 0.002, 0.002,
                                                0.002, 0.002, 0.002, 0.001, 0.0015,
                                                0.0015, 0.0016, 0.0016, 0.0017, 0.0019,
                                                0.0022, 0.003 };
    Double_t chPionMinCMS7TeV_yval[22]      = { 4.44, 5.575, 5.588, 5.041, 4.397,
                                                3.794, 3.246, 2.769, 2.375, 2.032,
                                                1.739, 1.499, 1.29, 1.117, 0.977,
                                                0.8601, 0.7628, 0.6559, 0.5884, 0.507,
                                                0.4364, 0.3664 };
    Double_t chPionMinCMS7TeV_yErrMin[22]   = { 0.406, 0.179, 0.134, 0.107, 0.090,
                                                0.077, 0.066, 0.056, 0.048, 0.041,
                                                0.035, 0.031, 0.028, 0.025, 0.0224,
                                                0.0200, 0.0188, 0.0185, 0.0180, 0.0176, 
                                                0.0183, 0.0186};
    Double_t chPionMinCMS7TeV_yStatMin[22]  = { 0.003, 0.003, 0.003, 0.002, 0.002, 
                                                0.002, 0.002, 0.002, 0.002, 0.002,
                                                0.002, 0.002, 0.002, 0.001, 0.0015,
                                                0.0015, 0.0016, 0.0016, 0.0017, 0.0018,
                                                0.002, 0.0028 };
  
    TH1D* histoChargedPionPlusSpecLowPtStat7TeVCMS  = new TH1D("histoChargedPionPlusSpecLowPtStat7TeVCMS","histoChargedPionPlusSpecLowPtStat7TeVCMS", 22, ptbins7TeV);
    TH1D* histoChargedPionPlusSpecLowPtSys7TeVCMS   = new TH1D("histoChargedPionPlusSpecLowPtSys7TeVCMS","histoChargedPionPlusSpecLowPtSys7TeVCMS", 22, ptbins7TeV);
    TH1D* histoChargedPionMinusSpecLowPtStat7TeVCMS = new TH1D("histoChargedPionMinusSpecLowPtStat7TeVCMS","histoChargedPionMinusSpecLowPtStat7TeVCMS", 22, ptbins7TeV);
    TH1D* histoChargedPionMinusSpecLowPtSys7TeVCMS  = new TH1D("histoChargedPionMinusSpecLowPtSys7TeVCMS","histoChargedPionMinusSpecLowPtSys7TeVCMS", 22, ptbins7TeV);
    
    for(Int_t i=1; i<23; i++){
        histoChargedPionPlusSpecLowPtStat7TeVCMS->SetBinContent(i, chPionPlCMS7TeV_yval[i-1]);
        histoChargedPionPlusSpecLowPtStat7TeVCMS->SetBinError(i, chPionPlCMS7TeV_yStatMin[i-1]);
        histoChargedPionPlusSpecLowPtSys7TeVCMS->SetBinContent(i, chPionPlCMS7TeV_yval[i-1]);
        histoChargedPionPlusSpecLowPtSys7TeVCMS->SetBinError(i, chPionPlCMS7TeV_yErrMin[i-1]);
        
        histoChargedPionMinusSpecLowPtStat7TeVCMS->SetBinContent(i, chPionMinCMS7TeV_yval[i-1]);
        histoChargedPionMinusSpecLowPtStat7TeVCMS->SetBinError(i, chPionMinCMS7TeV_yStatMin[i-1]);
        histoChargedPionMinusSpecLowPtSys7TeVCMS->SetBinContent(i, chPionMinCMS7TeV_yval[i-1]);
        histoChargedPionMinusSpecLowPtSys7TeVCMS->SetBinError(i, chPionMinCMS7TeV_yErrMin[i-1]);
    }
    TH1D*   histoChargedPionSpecLowPtSys7TeVCMS     = (TH1D*)histoChargedPionMinusSpecLowPtSys7TeVCMS->Clone("histoChargedPionSpecLowPtSys7TeVCMS");
    histoChargedPionSpecLowPtSys7TeVCMS->Add(histoChargedPionPlusSpecLowPtSys7TeVCMS);
    histoChargedPionSpecLowPtSys7TeVCMS->Scale(0.5);
    TH1D*   histoChargedPionSpecLowPtStat7TeVCMS    = (TH1D*)histoChargedPionMinusSpecLowPtStat7TeVCMS->Clone("histoChargedPionSpecLowPtStat7TeVCMS");
    histoChargedPionSpecLowPtStat7TeVCMS->Add(histoChargedPionPlusSpecLowPtStat7TeVCMS);
    histoChargedPionSpecLowPtStat7TeVCMS->Scale(0.5);

    for (Int_t i = 1; i < histoChargedPionSpecLowPtSys7TeVCMS->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStat7TeVCMS->SetBinContent(i, histoChargedPionSpecLowPtStat7TeVCMS->GetBinContent(i)/histoChargedPionSpecLowPtStat7TeVCMS->GetBinCenter(i)/(2*TMath::Pi())*0.78);
        histoChargedPionSpecLowPtStat7TeVCMS->SetBinError(i, i, histoChargedPionSpecLowPtStat7TeVCMS->GetBinError(i)/histoChargedPionSpecLowPtStat7TeVCMS->GetBinCenter(i)/(2*TMath::Pi())*0.78);
        histoChargedPionSpecLowPtSys7TeVCMS->SetBinContent(i, histoChargedPionSpecLowPtSys7TeVCMS->GetBinContent(i)/histoChargedPionSpecLowPtSys7TeVCMS->GetBinCenter(i)/(2*TMath::Pi())*0.78);
        histoChargedPionSpecLowPtSys7TeVCMS->SetBinError(i, histoChargedPionSpecLowPtSys7TeVCMS->GetBinError(i)/histoChargedPionSpecLowPtSys7TeVCMS->GetBinCenter(i)/(2*TMath::Pi())*0.78);
      
        Double_t fractionalSystematicError          = 0;
        if (histoChargedPionMinusSpecLowPtSys7TeVCMS->GetBinContent(i) != 0 && histoChargedPionPlusSpecLowPtSys7TeVCMS->GetBinContent(i) != 0){
            fractionalSystematicError               = ( histoChargedPionMinusSpecLowPtSys7TeVCMS->GetBinError(i)/histoChargedPionMinusSpecLowPtSys7TeVCMS->GetBinContent(i)*100 +
                                                        histoChargedPionPlusSpecLowPtSys7TeVCMS->GetBinError(i)/histoChargedPionPlusSpecLowPtSys7TeVCMS->GetBinContent(i)*100)/2;
        }
        histoChargedPionSpecLowPtSys7TeVCMS->SetBinError(i, histoChargedPionSpecLowPtSys7TeVCMS->GetBinContent(i)*fractionalSystematicError/100.);
    }
    
    // *********************************************************************************************************************
    // **************************** full pt spectra 7 TeV identified hadrons ***********************************************
    // **************************** arXiv:1504.00024, 1601.03658 ***********************************************************
    // *********************************************************************************************************************
    TFile* fileChargedPionSpectraPublishedPP7TeV        = new TFile("ExternalInput/IdentifiedCharged/pp7TeV.mb.fullpT.INEL.20150803.root");
    TH1D*   histoChargedPionSpecPubStatPP7TeV           = (TH1D*)fileChargedPionSpectraPublishedPP7TeV->Get("hstat_pp7_pion_sum");
    TH1D*   histoChargedPionSpecPubSystPP7TeV           = (TH1D*)fileChargedPionSpectraPublishedPP7TeV->Get("hsys_pp7_pion_sum");
    histoChargedPionSpecPubStatPP7TeV->Scale(0.5);
    histoChargedPionSpecPubSystPP7TeV->Scale(0.5);
    TH1D*   histoChargedKaonSpecPubStatPP7TeV           = (TH1D*)fileChargedPionSpectraPublishedPP7TeV->Get("hstat_pp7_kaon_sum");
    TH1D*   histoChargedKaonSpecPubSystPP7TeV           = (TH1D*)fileChargedPionSpectraPublishedPP7TeV->Get("hsys_pp7_kaon_sum");
    histoChargedKaonSpecPubStatPP7TeV->Scale(0.5);
    histoChargedKaonSpecPubSystPP7TeV->Scale(0.5);
    TH1D*   histoProtonSpecPubStatPP7TeV                = (TH1D*)fileChargedPionSpectraPublishedPP7TeV->Get("hstat_pp7_proton_sum");
    TH1D*   histoProtonSpecPubSystPP7TeV                = (TH1D*)fileChargedPionSpectraPublishedPP7TeV->Get("hsys_pp7_proton_sum");
    histoProtonSpecPubStatPP7TeV->Scale(0.5);
    histoProtonSpecPubSystPP7TeV->Scale(0.5);


    // *********************************************************************************************************************
    // ************************** High Pt  charged pions 7 TeV ALICE *******************************************************
    // ************************** rather old spectra ***********************************************************************
    // *********************************************************************************************************************
    TFile *fPionChargedNew                                          = TFile::Open("ExternalInput/IdentifiedCharged/charged_pion_pectrum_correctedfraction-20120314.root");
    TH1D *histoChargedPionSpecHighPtStat7TeVALICE                   = (TH1D*)fPionChargedNew->Get("hPionSpectrum_PP");
    histoChargedPionSpecHighPtStat7TeVALICE->Scale(0.50);
    TGraphAsymmErrors *graphChargedPionSpecHighPtSys7TeVALICE       = (TGraphAsymmErrors*)fPionChargedNew->Get("hPionGraphAsymUncertainties_PP");
    graphChargedPionSpecHighPtSys7TeVALICE = ScaleGraph(graphChargedPionSpecHighPtSys7TeVALICE, 0.5);
    graphChargedPionSpecHighPtSys7TeVALICE->Print();
    TGraphAsymmErrors *graphChargedPionSpecHighPtSys7TeVALICECopy   = (TGraphAsymmErrors*)fPionChargedNew->Get("hPionGraphAsymUncertainties_PP");
    Double_t* yValue    = graphChargedPionSpecHighPtSys7TeVALICECopy->GetY();
    Int_t counter       = 0;
    while ((abs(yValue[counter] - 0.) )< 1e-10){
        cout << "removed point" << counter << endl;
        counter++;
        graphChargedPionSpecHighPtSys7TeVALICE->RemovePoint(0);   
    }
    
    // *********************************************************************************************************************
    // ************************** High Pt charged pions 7 TeV ALICE 26 Mar 2012 ********************************************
    // ************************** rather old spectra ***********************************************************************
    // *********************************************************************************************************************
    TFile *fPionLowPt = TFile::Open("ExternalInput/IdentifiedCharged/SpectraCombined_13_Marzo_2013.root");
  
    TH1D* histoChargedPionPlusSpecLowPtStat7TeVALICE    = (TH1D*)fPionLowPt->Get("hComb_ITSsa0_TPC0_TOF0_HMPID0");
    TH1D* histoChargedPionPlusSpecLowPtSys7TeVALICE     = (TH1D*)fPionLowPt->Get("Syst_hComb_ITSsa0_TPC0_TOF0_HMPID0");
    TH1D* histoChargedPionMinusSpecLowPtStat7TeVALICE   = (TH1D*)fPionLowPt->Get("hComb_ITSsa3_TPC3_TOF3_HMPID3");
    TH1D* histoChargedPionMinusSpecLowPtSys7TeVALICE    = (TH1D*)fPionLowPt->Get("Syst_hComb_ITSsa3_TPC3_TOF3_HMPID3");
    
    TH1D*   histoChargedPionSpecLowPtStat7TeVALICE      = (TH1D*)histoChargedPionMinusSpecLowPtStat7TeVALICE->Clone("histoChargedPionSpecLowPtStat7TeVALICE");
    histoChargedPionSpecLowPtStat7TeVALICE->Add(histoChargedPionPlusSpecLowPtStat7TeVALICE);
    histoChargedPionSpecLowPtStat7TeVALICE->Scale(0.5);
    TH1D*   histoChargedPionSpecLowPtSys7TeVALICE       = (TH1D*)histoChargedPionMinusSpecLowPtSys7TeVALICE->Clone("histoChargedPionSpecLowPtSys7TeVALICE");  
    
    for (Int_t i = 1; i < histoChargedPionSpecLowPtSys7TeVALICE->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStat7TeVALICE->SetBinContent(i, histoChargedPionSpecLowPtStat7TeVALICE->GetBinContent(i)/histoChargedPionSpecLowPtStat7TeVALICE->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStat7TeVALICE->SetBinError(i, i, histoChargedPionSpecLowPtStat7TeVALICE->GetBinError(i)/histoChargedPionSpecLowPtStat7TeVALICE->GetBinCenter(i)/(2*TMath::Pi()));
        Double_t error                      = ( histoChargedPionMinusSpecLowPtSys7TeVALICE->GetBinContent(i) *histoChargedPionMinusSpecLowPtStat7TeVALICE->GetBinContent(i) +
                                                histoChargedPionPlusSpecLowPtSys7TeVALICE->GetBinContent(i) *histoChargedPionPlusSpecLowPtStat7TeVALICE->GetBinContent(i))/2.;
        histoChargedPionSpecLowPtSys7TeVALICE->SetBinContent(i, histoChargedPionSpecLowPtStat7TeVALICE->GetBinContent(i));
        histoChargedPionSpecLowPtSys7TeVALICE->SetBinError(i,error/histoChargedPionSpecLowPtSys7TeVALICE->GetBinCenter(i)/(2*TMath::Pi()));
        
    }

    
    // *********************************************************************************************************************
    // *********************************** ALICE charged hadron results 7 TeV **********************************************
    // *********************************** Eur.Phys.J. C73 (2013) no.12, 2662 arXiv:1307.1093 ******************************
    // *********************************************************************************************************************
    
    Double_t chHadALICE_PP7TeV_xval[]           = { 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 
                                                    0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.05, 1.15, 
                                                    1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.1, 2.3, 
                                                    2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 3.9, 4.25, 4.75, 
                                                    5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 
                                                    13.5, 14.5, 15.5, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 29.0, 
                                                    31.0, 33.0, 35.0, 38.0, 42.5, 47.5 };
    Double_t chHadALICE_PP7TeV_xerrminus[]      = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
                                                    0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.050, 0.050, 
                                                    0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.100, 0.100, 
                                                    0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.25, 0.25, 
                                                    0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
                                                    0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
                                                    1.0, 1.0, 1.0, 2.0, 2.5, 2.5 };
    Double_t chHadALICE_PP7TeV_xerrplus[]       = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
                                                    0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.050, 0.050, 
                                                    0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.100, 0.100, 
                                                    0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.25, 0.25, 
                                                    0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
                                                    0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
                                                    1.0, 1.0, 1.0, 2.0, 2.5, 2.5 };
    Double_t chHadALICE_PP7TeV_yval[]           = { 490.7, 389.4, 301.1, 232.3, 179.8, 140.3, 110.3, 87.47, 70.0, 
                                                    56.52, 46.04, 37.71, 31.06, 25.77, 21.55, 18.13, 15.33, 11.94, 8.722, 
                                                    6.463, 4.84, 3.65, 2.775, 2.131, 1.662, 1.315, 1.045, 0.7477, 0.4863, 
                                                    0.3226, 0.218, 0.1506, 0.1058, 0.07536, 0.05444, 0.03993, 0.02952, 0.01822, 0.009507, 
                                                    0.005201, 0.003045, 0.001835, 0.001159, 6.256E-4, 2.931E-4, 1.483E-4, 8.094E-5, 4.742E-5, 2.804E-5, 
                                                    1.786E-5, 1.141E-5, 7.72E-6, 4.512E-6, 2.206E-6, 1.211E-6, 6.881E-7, 4.184E-7, 2.687E-7, 1.721E-7, 
                                                    1.071E-7, 8.669E-8, 4.078E-8, 2.977E-8, 1.818E-8, 1.057E-8 };
    Double_t chHadALICE_PP7TeV_yerrminus[]      = { 37.20013440835933, 25.60019531175495, 19.800252523642218, 15.4, 12.4, 10.2, 8.1, 6.390031298827886, 5.1100097847264445, 
                                                    4.130012106519786, 3.360014880919428, 2.750018181758077, 2.2700220263248547, 1.880026595556563, 1.5700318468107584, 
                                                    1.3200378782444087, 1.120044641967453, 0.87, 0.6370031397096878, 
                                                    0.47200423726911606, 0.35400564967243103, 0.26700187265260894, 0.20300246303924493, 0.15600320509528, 0.12100413216084813,
                                                    0.09600520819205591, 0.07600657866263946, 0.05460146518180625, 0.03550126758300328, 
                                                    0.023600847442411893, 0.01590125781188394, 0.011001818031580053, 0.007700649323271383, 0.005500736314349198, 
                                                    0.003980615530291767, 0.0029206163733020464, 0.002160578626201787, 0.001330338302838793, 6.942081820318744E-4, 
                                                    3.801894264705425E-4, 2.2214409737825582E-4, 1.3413426109685773E-4, 8.514693182963201E-5, 4.575784085815239E-5,
                                                    2.145250568115529E-5, 1.0846197490364998E-5, 5.952495275092623E-6, 3.4988569562072696E-6, 2.025956564193813E-6, 
                                                    1.3026895255585651E-6, 8.516454661418682E-7, 5.831517812713942E-7, 3.4065818645674734E-7, 1.7404022523543226E-7,
                                                    9.974467404327913E-8, 6.03675409471017E-8, 3.9062257999250375E-8, 2.6923595599399424E-8, 1.8958902921846507E-8, 
                                                    1.3261221663180206E-8, 1.1226116870939835E-8, 6.842514157822401E-9, 4.163291966701351E-9, 2.7148664792214E-9,
                                                    1.8469975636150688E-9 };
    Double_t chHadALICE_PP7TeV_yerrplus[]       = { 37.20013440835933, 25.60019531175495, 19.800252523642218, 15.4, 12.4, 10.2, 8.1, 6.390031298827886, 5.1100097847264445, 
                                                    4.130012106519786, 3.360014880919428, 2.750018181758077, 2.2700220263248547, 1.880026595556563, 
                                                    1.5700318468107584, 1.3200378782444087, 1.120044641967453, 0.87, 0.6370031397096878, 
                                                    0.47200423726911606, 0.35400564967243103, 0.26700187265260894, 0.20300246303924493, 0.15600320509528,
                                                    0.12100413216084813, 0.09600520819205591, 0.07600657866263946, 0.05460146518180625, 0.03550126758300328, 
                                                    0.023600847442411893, 0.01590125781188394, 0.011001818031580053, 0.007700649323271383, 0.005500736314349198,
                                                    0.003980615530291767, 0.0029206163733020464, 0.002160578626201787, 0.001330338302838793, 6.942081820318744E-4, 
                                                    3.801894264705425E-4, 2.2214409737825582E-4, 1.3413426109685773E-4, 8.514693182963201E-5, 4.575784085815239E-5,
                                                    2.145250568115529E-5, 1.0846197490364998E-5, 5.952495275092623E-6, 3.4988569562072696E-6, 2.025956564193813E-6, 
                                                    1.3026895255585651E-6, 8.516454661418682E-7, 5.831517812713942E-7, 3.4065818645674734E-7, 1.7404022523543226E-7,
                                                    9.974467404327913E-8, 6.03675409471017E-8, 3.9062257999250375E-8, 2.6923595599399424E-8, 1.8958902921846507E-8, 
                                                    1.3261221663180206E-8, 1.1226116870939835E-8, 6.842514157822401E-9, 4.163291966701351E-9, 2.7148664792214E-9, 
                                                    1.8469975636150688E-9 };
    Double_t chHadALICE_PP7TeV_ystatminus[]     = { 0.1, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.02, 0.01, 
                                                    0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.0, 0.002, 
                                                    0.002, 0.002, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 4.0E-4, 3.0E-4, 
                                                    2.0E-4, 2.0E-4, 2.0E-4, 1.0E-4, 9.0E-5, 7.0E-5, 6.0E-5, 5.0E-5, 3.0E-5, 1.7E-5, 
                                                    1.2E-5, 8.0E-6, 6.0E-6, 5.0E-6, 2.3E-6, 1.5E-6, 1.0E-6, 7.1E-7, 5.2E-7, 3.8E-7, 
                                                    2.9E-7, 2.3E-7, 1.79E-7, 9.2E-8, 6.1E-8, 4.3E-8, 3.1E-8, 2.31E-8, 1.78E-8, 1.38E-8, 
                                                    1.05E-8, 9.11E-9, 6.1E-9, 3.49E-9, 2.32E-9, 1.65E-9 };
    Double_t chHadALICE_PP7TeV_ystatplus[]      = { 0.1, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.02, 0.01, 
                                                    0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.0, 0.002, 
                                                    0.002, 0.002, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 4.0E-4, 3.0E-4, 
                                                    2.0E-4, 2.0E-4, 2.0E-4, 1.0E-4, 9.0E-5, 7.0E-5, 6.0E-5, 5.0E-5, 3.0E-5, 1.7E-5, 
                                                    1.2E-5, 8.0E-6, 6.0E-6, 5.0E-6, 2.3E-6, 1.5E-6, 1.0E-6, 7.1E-7, 5.2E-7, 3.8E-7, 
                                                    2.9E-7, 2.3E-7, 1.79E-7, 9.2E-8, 6.1E-8, 4.3E-8, 3.1E-8, 2.31E-8, 1.78E-8, 1.38E-8, 
                                                    1.05E-8, 9.11E-9, 6.1E-9, 3.49E-9, 2.32E-9, 1.65E-9 };
    Int_t chHadALICE_PP7TeV_numpoints = 65;
    
    Double_t chHadALICE_PP7TeV_ysysminus[65];
    Double_t chHadALICE_PP7TeV_ysysplus[65];
    ExtractSystematicFromTotal(chHadALICE_PP7TeV_numpoints, chHadALICE_PP7TeV_yerrminus, chHadALICE_PP7TeV_ystatminus, chHadALICE_PP7TeV_ysysminus );
    ExtractSystematicFromTotal(chHadALICE_PP7TeV_numpoints, chHadALICE_PP7TeV_yerrplus, chHadALICE_PP7TeV_ystatplus, chHadALICE_PP7TeV_ysysplus );
    
    TGraphAsymmErrors* graphChHadALICEStatPP7TeV    = new TGraphAsymmErrors(chHadALICE_PP7TeV_numpoints, chHadALICE_PP7TeV_xval, chHadALICE_PP7TeV_yval, chHadALICE_PP7TeV_xerrminus,
                                                                            chHadALICE_PP7TeV_xerrplus, chHadALICE_PP7TeV_ystatminus, chHadALICE_PP7TeV_ystatplus);
    TGraphAsymmErrors* graphChHadALICESysPP7TeV     = new TGraphAsymmErrors(chHadALICE_PP7TeV_numpoints, chHadALICE_PP7TeV_xval, chHadALICE_PP7TeV_yval, chHadALICE_PP7TeV_xerrminus, 
                                                                            chHadALICE_PP7TeV_xerrplus, chHadALICE_PP7TeV_ysysminus, chHadALICE_PP7TeV_ysysplus);
 
    // *********************************************************************************************************************
    // *********************************** CMS charged hadron results 7 TeV ************************************************
    // *********************************** JHEP 1108 (2011) 086, arXiv:1104.3547 *******************************************
    // *********************************************************************************************************************      
    Double_t chHadCMS_PP7TeV_xval[27]           = { 0.5, 0.7, 0.9, 1.1, 1.4, 1.8, 2.2, 2.8, 3.6, 4.4, 
                                                    5.2, 6.0, 6.8, 8.2, 10.2, 12.2, 15.2, 19.2, 23.2, 31.2, 
                                                    43.2, 55.2, 71.2, 91.2, 111.2, 141.2, 181.2 };
    Double_t chHadCMS_PP7TeV_xerrminus[27]      = { 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.4, 0.4, 0.4, 
                                                    0.4, 0.4, 0.4, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 6.0, 
                                                    6.0, 6.0, 10.0, 10.0, 10.0, 20.0, 20.0 };
    Double_t chHadCMS_PP7TeV_xerrplus[27]       = { 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.4, 0.4, 0.4, 
                                                    0.4, 0.4, 0.4, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 6.0, 
                                                    6.0, 6.0, 10.0, 10.0, 10.0, 20.0, 20.0 };
    Double_t chHadCMS_PP7TeV_yval[27]           = { 1.68022, 0.71916, 0.33621, 0.171979, 0.070213, 0.024796, 0.0099366, 0.0029975, 7.7056E-4, 2.4402E-4, 
                                                    9.0609E-5, 3.8438E-5, 1.78328E-5, 5.7719E-6, 1.52171E-6, 5.1369E-7, 1.3567E-7, 3.3087E-8, 1.01372E-8, 1.68816E-9, 
                                                    2.2497E-10, 4.4431E-11, 8.4116E-12, 1.5299E-12, 3.9626E-13, 8.1866E-14, 1.2708E-14 };
    Double_t chHadCMS_PP7TeV_yerrminus[27]      = { 0.0796102267299874, 0.0340101470152659, 0.015880113349721402, 0.008113084370324273, 0.003308034008289516, 
                                                    0.0011670274204147906, 4.670207275913993E-4, 1.4070799550842874E-4, 3.612621347442879E-5, 1.1435357449594656E-5, 
                                                    4.244435887135062E-6, 1.8018027084006728E-6, 8.368496220946749E-7, 2.7064848050561816E-7, 7.187439043219776E-8, 
                                                    2.4662076149424242E-8, 6.577688727813136E-9, 1.7259855155823295E-9, 5.810505313653882E-10, 1.0056383544793824E-10, 
                                                    1.629800601300662E-11, 3.4921698125950288E-12, 7.393300886613503E-13, 1.3229618286254523E-13, 3.115363542188937E-14, 
                                                    6.9566728397992094E-15, 1.496573753611896E-15 };
    Double_t chHadCMS_PP7TeV_yerrplus[27]       = { 0.0796102267299874, 0.0340101470152659, 0.015880113349721402, 0.008113084370324273, 0.003308034008289516, 
                                                    0.0011670274204147906, 4.670207275913993E-4, 1.4070799550842874E-4, 3.612621347442879E-5, 1.1435357449594656E-5, 
                                                    4.244435887135062E-6, 1.8018027084006728E-6, 8.368496220946749E-7, 2.7064848050561816E-7, 7.187439043219776E-8, 
                                                    2.4662076149424242E-8, 6.577688727813136E-9, 1.7259855155823295E-9, 5.810505313653882E-10, 1.0056383544793824E-10, 
                                                    1.629800601300662E-11, 3.4921698125950288E-12, 7.393300886613503E-13, 1.3229618286254523E-13, 3.115363542188937E-14, 
                                                    6.9566728397992094E-15, 1.496573753611896E-15 };
    Double_t chHadCMS_PP7TeV_ystatminus[27]     = { 1.9E-4, 1.0E-4, 6.0E-5, 3.7E-5, 1.5E-5, 8.0E-6, 4.4E-6, 1.5E-6, 6.7E-7, 3.5E-7, 
                                                    1.94E-7, 1.17E-7, 7.48E-8, 2.38E-8, 1.092E-8, 5.76E-9, 1.745E-9, 7.49E-10, 3.236E-10, 5.201E-11, 
                                                    8.75E-12, 1.845E-12, 4.517E-13, 7.32E-14, 1.039E-14, 2.959E-15, 1.097E-15 };
    Double_t chHadCMS_PP7TeV_ystatplus[27]      = { 1.9E-4, 1.0E-4, 6.0E-5, 3.7E-5, 1.5E-5, 8.0E-6, 4.4E-6, 1.5E-6, 6.7E-7, 3.5E-7, 
                                                    1.94E-7, 1.17E-7, 7.48E-8, 2.38E-8, 1.092E-8, 5.76E-9, 1.745E-9, 7.49E-10, 3.236E-10, 5.201E-11, 
                                                    8.75E-12, 1.845E-12, 4.517E-13, 7.32E-14, 1.039E-14, 2.959E-15, 1.097E-15 };
    Int_t chHadCMS_PP7TeV_numpoints             = 27;
    Double_t chHadCMS_PP7TeV_ysysminus[27];
    Double_t chHadCMS_PP7TeV_ysysplus[27];
    ExtractSystematicFromTotal(chHadCMS_PP7TeV_numpoints, chHadCMS_PP7TeV_yerrminus, chHadCMS_PP7TeV_ystatminus, chHadCMS_PP7TeV_ysysminus );
    ExtractSystematicFromTotal(chHadCMS_PP7TeV_numpoints, chHadCMS_PP7TeV_yerrplus, chHadCMS_PP7TeV_ystatplus, chHadCMS_PP7TeV_ysysplus );
    
    TGraphAsymmErrors* graphChHadCMSStatPP7TeV      = new TGraphAsymmErrors(chHadCMS_PP7TeV_numpoints, chHadCMS_PP7TeV_xval, chHadCMS_PP7TeV_yval, chHadCMS_PP7TeV_xerrminus,
                                                                            chHadCMS_PP7TeV_xerrplus, chHadCMS_PP7TeV_ystatminus, chHadCMS_PP7TeV_ystatplus);
    TGraphAsymmErrors* graphChHadCMSSysPP7TeV       = new TGraphAsymmErrors(chHadCMS_PP7TeV_numpoints, chHadCMS_PP7TeV_xval, chHadCMS_PP7TeV_yval, chHadCMS_PP7TeV_xerrminus, 
                                                                            chHadCMS_PP7TeV_xerrplus, chHadCMS_PP7TeV_ysysminus, chHadCMS_PP7TeV_ysysplus);
    
    // -----------------------------------------------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------- Read 2.76 TeV spectra ----------------------------------------------------------------------------------------
    // -----------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // *********************************************************************************************************************
    // *********************************** CMS charged pion results 2.76 TeV ***********************************************
    //     @article{CMS:2012gva,
    //         author         = "CMS Collaboration",
    //         title          = "{Spectra of identified charged hadrons at sqrt(s) = 0.9, 2.76 and 7 TeV via tracker energy loss}",
    //         collaboration  = "CMS",
    //         year           = "2012",
    //         reportNumber   = "CMS-PAS-FSQ-12-014",
    //         SLACcitation   = "%%CITATION = CMS-PAS-FSQ-12-014;%%"
    //     }
    // *********************************************************************************************************************
    Double_t chPionCMS_xval[22]             = { 0.125, 0.175, 0.225, 0.275, 0.325,
                                                0.375, 0.425, 0.475, 0.525, 0.575,
                                                0.625, 0.675, 0.725, 0.775, 0.825,
                                                0.875, 0.925, 0.975, 1.025, 1.075, 
                                                1.125, 1.175 };
    Double_t chPionPlCMS_xerrminus[22]      = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
                                                0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
                                                0.025, 0.025, 0.025 };
    Double_t ptbins[23]                     = { 0.1, 0.15, 0.2, 0.25, 0.3,
                                                0.35, 0.4, 0.45, 0.5, 0.55,
                                                0.6, 0.65, 0.7, 0.75, 0.8, 
                                                0.85, 0.9, 0.95, 1.0, 1.05,
                                                1.1, 1.15, 1.2};
    Double_t chPionPlCMS_yval[22]           = { 3.828, 4.749, 4.668, 4.201, 3.66, 3.153, 2.69, 2.275, 1.928, 
                                                1.638, 1.389, 1.185, 1.011, 0.8673, 0.7424, 0.6375, 0.553, 0.4952, 0.4314, 
                                                0.3836, 0.3383, 0.2812 };
    Double_t chPionPlCMS_yErrMin[22]        = { 0.339, 0.147, 0.108, 0.090, 0.077, 0.064, 0.055, 0.046, 0.039, 
                                                0.033, 0.028, 0.025, 0.021, 0.0192, 0.0170, 0.0151, 0.0140, 0.0141, 0.0135, 
                                                0.0132, 0.0133, 0.0124 };
    Double_t chPionPlCMS_yStatMin[22]       = { 0.003, 0.003, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 
                                                0.001, 0.001, 0.001, 0.001, 0.0012, 0.0012, 0.0012, 0.0013, 0.0012, 0.0013, 
                                                0.0014, 0.0016, 0.0022 };
    
    Double_t chPionMinCMS_yval[22]          = { 3.596, 4.526, 4.593, 4.159, 3.612, 3.109, 2.648, 2.243, 1.903, 
                                                1.616, 1.375, 1.168, 0.9971, 0.8497, 0.7288, 0.6328, 0.5609, 0.4865, 0.4147, 
                                                0.37, 0.321, 0.2809 };
    Double_t chPionMinCMS_yErrMin[22]       = { 0.328, 0.145, 0.109, 0.089, 0.074, 0.063, 0.054, 0.045, 0.038, 
                                                0.033, 0.028, 0.024, 0.0211, 0.0187, 0.0166, 0.0147, 0.0137, 0.0134, 0.0123, 
                                                0.0122, 0.0120, 0.0112};
    Double_t chPionMinCMS_yStatMin[22]      = { 0.003, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 
                                                0.001, 0.001, 0.001, 0.0012, 0.0012, 0.0012, 0.0012, 0.0012, 0.0012, 0.0013, 
                                                0.0014, 0.0015, 0.0021 };
  
    TH1D* histoChargedPionPlusSpecLowPtStat2760GeVCMS   = new TH1D("histoChargedPionPlusSpecLowPtStat2760GeVCMS","histoChargedPionPlusSpecLowPtStat2760GeVCMS", 22, ptbins);
    TH1D* histoChargedPionPlusSpecLowPtSys2760GeVCMS    = new TH1D("histoChargedPionPlusSpecLowPtSys2760GeVCMS","histoChargedPionPlusSpecLowPtSys2760GeVCMS", 22, ptbins);
    TH1D* histoChargedPionMinusSpecLowPtStat2760GeVCMS  = new TH1D("histoChargedPionMinusSpecLowPtStat2760GeVCMS","histoChargedPionMinusSpecLowPtStat2760GeVCMS", 22, ptbins);
    TH1D* histoChargedPionMinusSpecLowPtSys2760GeVCMS   = new TH1D("histoChargedPionMinusSpecLowPtSys2760GeVCMS","histoChargedPionMinusSpecLowPtSys2760GeVCMS", 22, ptbins);
    
    for(Int_t i=1; i<23; i++){
        histoChargedPionPlusSpecLowPtStat2760GeVCMS->SetBinContent(i, chPionPlCMS_yval[i-1]);
        histoChargedPionPlusSpecLowPtStat2760GeVCMS->SetBinError(i, chPionPlCMS_yStatMin[i-1]);
        histoChargedPionPlusSpecLowPtSys2760GeVCMS->SetBinContent(i, chPionPlCMS_yval[i-1]);
        histoChargedPionPlusSpecLowPtSys2760GeVCMS->SetBinError(i, chPionPlCMS_yErrMin[i-1]);
        
        histoChargedPionMinusSpecLowPtStat2760GeVCMS->SetBinContent(i, chPionMinCMS_yval[i-1]);
        histoChargedPionMinusSpecLowPtStat2760GeVCMS->SetBinError(i, chPionMinCMS_yStatMin[i-1]);
        histoChargedPionMinusSpecLowPtSys2760GeVCMS->SetBinContent(i, chPionMinCMS_yval[i-1]);
        histoChargedPionMinusSpecLowPtSys2760GeVCMS->SetBinError(i, chPionMinCMS_yErrMin[i-1]);
    }
    TH1D*   histoChargedPionSpecLowPtSys2760GeVCMS      = (TH1D*)histoChargedPionMinusSpecLowPtSys2760GeVCMS->Clone("histoChargedPionSpecLowPtSys2760GeVCMS");
    histoChargedPionSpecLowPtSys2760GeVCMS->Add(histoChargedPionPlusSpecLowPtSys2760GeVCMS);
    histoChargedPionSpecLowPtSys2760GeVCMS->Scale(0.5);
    TH1D*   histoChargedPionSpecLowPtStat2760GeVCMS     = (TH1D*)histoChargedPionMinusSpecLowPtStat2760GeVCMS->Clone("histoChargedPionSpecLowPtStat2760GeVCMS");
    histoChargedPionSpecLowPtStat2760GeVCMS->Add(histoChargedPionPlusSpecLowPtStat2760GeVCMS);
    histoChargedPionSpecLowPtStat2760GeVCMS->Scale(0.5);

    for (Int_t i = 1; i < histoChargedPionSpecLowPtSys2760GeVCMS->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStat2760GeVCMS->SetBinError(i, histoChargedPionSpecLowPtSys2760GeVCMS->GetBinError(i)/histoChargedPionSpecLowPtStat2760GeVCMS->GetBinCenter(i)/(2*TMath::Pi())*0.78);
        histoChargedPionSpecLowPtStat2760GeVCMS->SetBinContent(i, histoChargedPionSpecLowPtStat2760GeVCMS->GetBinContent(i)/histoChargedPionSpecLowPtStat2760GeVCMS->GetBinCenter(i)/(2*TMath::Pi())*0.78);
        histoChargedPionSpecLowPtSys2760GeVCMS->SetBinContent(i, histoChargedPionSpecLowPtSys2760GeVCMS->GetBinContent(i)/histoChargedPionSpecLowPtSys2760GeVCMS->GetBinCenter(i)/(2*TMath::Pi())*0.78);
        Double_t fractionalSystematicError              = 0;
        if (histoChargedPionPlusSpecLowPtSys2760GeVCMS->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtStat2760GeVCMS->GetBinContent(i) != 0){
            fractionalSystematicError                   = ( histoChargedPionPlusSpecLowPtSys2760GeVCMS->GetBinError(i)/histoChargedPionPlusSpecLowPtSys2760GeVCMS->GetBinContent(i)*100 +
                                                            histoChargedPionMinusSpecLowPtStat2760GeVCMS->GetBinError(i)/histoChargedPionMinusSpecLowPtStat2760GeVCMS->GetBinContent(i)*100)/2;
        }
        histoChargedPionSpecLowPtSys2760GeVCMS->SetBinError(i, histoChargedPionSpecLowPtSys2760GeVCMS->GetBinContent(i)*fractionalSystematicError/100.);
    }
    
    
    
    // *********************************************************************************************************************
    // **************************** full pt spectra 2.76 TeV ***************************************************************
    // **************************** arXiv: 1401.1250, 1506.07287, 1601.03658 ***********************************************
    // *********************************************************************************************************************
    TFile* fileChargedPionSpectraPublishedPP2760GeV        = new TFile("ExternalInput/IdentifiedCharged/pp276.fullpT.INEL.20140504.root");
    TH1D*   histoChargedPionSpecPubStatPP2760GeV           = (TH1D*)fileChargedPionSpectraPublishedPP2760GeV->Get("hstat_pp276_pion_sum");
    TH1D*   histoChargedPionSpecPubSystPP2760GeV           = (TH1D*)fileChargedPionSpectraPublishedPP2760GeV->Get("hsys_pp276_pion_sum");
    histoChargedPionSpecPubStatPP2760GeV->Scale(0.5);
    histoChargedPionSpecPubSystPP2760GeV->Scale(0.5);
    TH1D*   histoChargedKaonSpecPubStatPP2760GeV           = (TH1D*)fileChargedPionSpectraPublishedPP2760GeV->Get("hstat_pp276_kaon_sum");
    TH1D*   histoChargedKaonSpecPubSystPP2760GeV           = (TH1D*)fileChargedPionSpectraPublishedPP2760GeV->Get("hsys_pp276_kaon_sum");
    histoChargedKaonSpecPubStatPP2760GeV->Scale(0.5);
    histoChargedKaonSpecPubSystPP2760GeV->Scale(0.5);
    TH1D*   histoProtonSpecPubStatPP2760GeV                = (TH1D*)fileChargedPionSpectraPublishedPP2760GeV->Get("hstat_pp276_proton_sum");
    TH1D*   histoProtonSpecPubSystPP2760GeV                = (TH1D*)fileChargedPionSpectraPublishedPP2760GeV->Get("hsys_pp276_proton_sum");
    histoProtonSpecPubStatPP2760GeV->Scale(0.5);
    histoProtonSpecPubSystPP2760GeV->Scale(0.5);
    
    // *********************************************************************************************************************
    // *************************** high pt spectra pion 2.76 TeV ***********************************************************
    // *************************** rather old spectra **********************************************************************
    // *********************************************************************************************************************
    TFile* fileChargedPionSpectraHighPtFinal            = new TFile("ExternalInputPbPb/IdentifiedCharged/SpectraHighPtPionFinal_20131101.root");
    TH1D*   histoChargedPionSpecHighPtStatPP2760GeV     = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrum_pp2760");
    TH1D*   histoChargedPionSpecHighPtSystPP2760GeV     = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrumSyst_pp2760");
    histoChargedPionSpecHighPtStatPP2760GeV->Scale(0.5);
    histoChargedPionSpecHighPtSystPP2760GeV->Scale(0.5);
    
    // *********************************************************************************************************************
    // ************************** high pt spectra kaon 2.76 TeV ************************************************************
    // ************************** rather old spectra ***********************************************************************
    // *********************************************************************************************************************
    TFile* fileChargedKaonSpectraHighPtFinal            = new TFile("ExternalInputPbPb/IdentifiedCharged/SpectraHighPtKaonFinal_20131108.root");
    TH1D* histoChargedKaonSpecHighPtStatPP              = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrum_pp2760");
    TH1D* histoChargedKaonSpecHighPtSystPP              = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrumSyst_pp2760");
    histoChargedKaonSpecHighPtStatPP->Scale(0.5);
    histoChargedKaonSpecHighPtSystPP->Scale(0.5);

    // *********************************************************************************************************************
    // ************************** K0s spectra 2.76 TeV *********************************************************************
    // ************************** published results + Simone Schuchmann's points *******************************************
    // ************************** send 20.4.2016 ***************************************************************************
    // *********************************************************************************************************************
    TFile* fileNeutralStrangeSpectra                    = new TFile("ExternalInput/NeutralStrange/pp_spectra_18_02_16_extrapol_etaErr_2760GeV.root");
    TH1D* histoNeutralKaonSpecStatPP2760GeV             = (TH1D*)fileNeutralStrangeSpectra->Get("K0s_corr_staterr_centpp_extrapolated");
    TH1D* histoNeutralKaonSpecSysPP2760GeV              = (TH1D*)fileNeutralStrangeSpectra->Get("K0s_corr_systerr_centpp_extrapolated");
    for (Int_t i = 1; i < histoNeutralKaonSpecStatPP2760GeV->GetNbinsX()+1; i++){
        histoNeutralKaonSpecStatPP2760GeV->SetBinContent(i, histoNeutralKaonSpecStatPP2760GeV->GetBinContent(i)/histoNeutralKaonSpecStatPP2760GeV->GetBinCenter(i)/(2*TMath::Pi()));
        histoNeutralKaonSpecStatPP2760GeV->SetBinError(i, histoNeutralKaonSpecStatPP2760GeV->GetBinError(i)/histoNeutralKaonSpecStatPP2760GeV->GetBinCenter(i)/(2*TMath::Pi()));
        histoNeutralKaonSpecSysPP2760GeV->SetBinContent(i, histoNeutralKaonSpecSysPP2760GeV->GetBinContent(i)/histoNeutralKaonSpecSysPP2760GeV->GetBinCenter(i)/(2*TMath::Pi()));
        histoNeutralKaonSpecSysPP2760GeV->SetBinError(i, histoNeutralKaonSpecSysPP2760GeV->GetBinError(i)/histoNeutralKaonSpecSysPP2760GeV->GetBinCenter(i)/(2*TMath::Pi()));
    }
    TH1D* histoLambda1115SpecStatPP2760GeV              = (TH1D*)fileNeutralStrangeSpectra->Get("L_corr_staterr_centpp_extrapolated");
    TH1D* histoLambda1115SpecSystPP2760GeV              = (TH1D*)fileNeutralStrangeSpectra->Get("L_corr_systerr_centpp_extrapolated");
    for (Int_t i = 1; i < histoLambda1115SpecStatPP2760GeV->GetNbinsX()+1; i++){
        histoLambda1115SpecStatPP2760GeV->SetBinContent(i, histoLambda1115SpecStatPP2760GeV->GetBinContent(i)/histoLambda1115SpecStatPP2760GeV->GetBinCenter(i)/(2*TMath::Pi()));
        histoLambda1115SpecStatPP2760GeV->SetBinError(i, histoLambda1115SpecStatPP2760GeV->GetBinError(i)/histoLambda1115SpecStatPP2760GeV->GetBinCenter(i)/(2*TMath::Pi()));
        histoLambda1115SpecSystPP2760GeV->SetBinContent(i, histoLambda1115SpecSystPP2760GeV->GetBinContent(i)/histoLambda1115SpecSystPP2760GeV->GetBinCenter(i)/(2*TMath::Pi()));
        histoLambda1115SpecSystPP2760GeV->SetBinError(i, histoLambda1115SpecSystPP2760GeV->GetBinError(i)/histoLambda1115SpecSystPP2760GeV->GetBinCenter(i)/(2*TMath::Pi()));
    }
    TH1D* histoAntiLambda1115SpecStatPP2760GeV           = (TH1D*)fileNeutralStrangeSpectra->Get("AL_corr_staterr_centpp_extrapolated");
    TH1D* histoAntiLambda1115SpecSystPP2760GeV           = (TH1D*)fileNeutralStrangeSpectra->Get("AL_corr_systerr_centpp_extrapolated");
    for (Int_t i = 1; i < histoAntiLambda1115SpecStatPP2760GeV->GetNbinsX()+1; i++){
        histoAntiLambda1115SpecStatPP2760GeV->SetBinContent(i, histoAntiLambda1115SpecStatPP2760GeV->GetBinContent(i)/histoAntiLambda1115SpecStatPP2760GeV->GetBinCenter(i)/(2*TMath::Pi()));
        histoAntiLambda1115SpecStatPP2760GeV->SetBinError(i, histoAntiLambda1115SpecStatPP2760GeV->GetBinError(i)/histoAntiLambda1115SpecStatPP2760GeV->GetBinCenter(i)/(2*TMath::Pi()));
        histoAntiLambda1115SpecSystPP2760GeV->SetBinContent(i, histoAntiLambda1115SpecSystPP2760GeV->GetBinContent(i)/histoAntiLambda1115SpecSystPP2760GeV->GetBinCenter(i)/(2*TMath::Pi()));
        histoAntiLambda1115SpecSystPP2760GeV->SetBinError(i, histoAntiLambda1115SpecSystPP2760GeV->GetBinError(i)/histoAntiLambda1115SpecSystPP2760GeV->GetBinCenter(i)/(2*TMath::Pi()));
    }
    
    // *********************************************************************************************************************
    // ************************** Low Pt pion data 2.76 TeV ****************************************************************
    // ************************** rather old spectra ***********************************************************************
    // *********************************************************************************************************************
    TFile* fileChPionSpectraLowPtPP2760GeVInternal      = new TFile("ExternalInput/IdentifiedCharged/276SpectraCombined_INELstat4072013.root");
    TH1D* histoChargedPionPlusSpecLowPtStatPP2760GeV    = (TH1D*)fileChPionSpectraLowPtPP2760GeVInternal->Get("pion_plus3");
    TH1D* histoChargedPionMinusSpecLowPtStatPP2760GeV   = (TH1D*)fileChPionSpectraLowPtPP2760GeVInternal->Get("pion_minus3");
    TH1D* histoChargedPionSpecLowPtStatPP2760GeV        = (TH1D*)histoChargedPionMinusSpecLowPtStatPP2760GeV->Clone("histoChargedPionSpecLowPtStatPP2760GeV");
    histoChargedPionSpecLowPtStatPP2760GeV->Add(histoChargedPionPlusSpecLowPtStatPP2760GeV);
    histoChargedPionSpecLowPtStatPP2760GeV->Scale(0.5);
    
    TFile* fileChPionSpectraLowPtPP2760GeVInternalSys   = new TFile("ExternalInput/IdentifiedCharged/276SpectraCombined_INEL4072013.root");
    TH1D* histoChargedPionPlusSpecLowPtSysPP2760GeV     = (TH1D*)fileChPionSpectraLowPtPP2760GeVInternalSys->Get("pion_plus3");
    TH1D* histoChargedPionMinusSpecLowPtSysPP2760GeV    = (TH1D*)fileChPionSpectraLowPtPP2760GeVInternalSys->Get("pion_minus3");
    TH1D* histoChargedPionSpecLowPtSysPP2760GeV         = (TH1D*)histoChargedPionPlusSpecLowPtSysPP2760GeV->Clone("histoChargedPionSpecLowPtSysPP2760GeV");
    histoChargedPionSpecLowPtSysPP2760GeV->Add(histoChargedPionMinusSpecLowPtSysPP2760GeV);
    histoChargedPionSpecLowPtSysPP2760GeV->Scale(0.5);
    
    for (Int_t i = 1; i < histoChargedPionSpecLowPtStatPP2760GeV->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStatPP2760GeV->SetBinContent(i, histoChargedPionSpecLowPtStatPP2760GeV->GetBinContent(i)/histoChargedPionSpecLowPtStatPP2760GeV->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStatPP2760GeV->SetBinError(i, histoChargedPionSpecLowPtStatPP2760GeV->GetBinError(i)/histoChargedPionSpecLowPtStatPP2760GeV->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtSysPP2760GeV->SetBinContent(i, histoChargedPionSpecLowPtSysPP2760GeV->GetBinContent(i)/histoChargedPionSpecLowPtSysPP2760GeV->GetBinCenter(i)/(2*TMath::Pi()));
        Double_t fractionalSystematicError                  = 0;
        if ( histoChargedPionMinusSpecLowPtSysPP2760GeV->GetBinContent(i) != 0 && histoChargedPionPlusSpecLowPtSysPP2760GeV->GetBinContent(i) != 0){
            fractionalSystematicError                       = ( histoChargedPionMinusSpecLowPtSysPP2760GeV->GetBinError(i)/histoChargedPionMinusSpecLowPtSysPP2760GeV->GetBinContent(i)*100 +
                                                                histoChargedPionPlusSpecLowPtSysPP2760GeV->GetBinError(i)/histoChargedPionPlusSpecLowPtSysPP2760GeV->GetBinContent(i)*100)/2;
        }
        histoChargedPionSpecLowPtSysPP2760GeV->SetBinError(i, histoChargedPionSpecLowPtSysPP2760GeV->GetBinContent(i)*fractionalSystematicError/100);
    }
    
    cout << "bis hier" << endl;

    // *********************************************************************************************************************
    // *********************************** ALICE charged hadron results 2.76 TeV *******************************************
    // *********************************** arXiv: Eur.Phys.J. C73 (2013) no.12, 2662 arXiv:1307.1093 ***********************
    // *********************************************************************************************************************
    Double_t chHadALICE_PP2760GeV_xval[65]      = { 0.175, 0.225, 0.275, 0.325, 0.375, 
                                                    0.425, 0.475, 0.525, 0.575, 0.625, 
                                                    0.675, 0.725, 0.775, 0.825, 0.875, 
                                                    0.925, 0.975, 1.05, 1.15, 1.25, 
                                                    1.35, 1.45, 1.55, 1.65, 1.75, 
                                                    1.85, 1.95, 2.1, 2.3, 2.5, 
                                                    2.7, 2.9, 3.1, 3.3, 3.5, 
                                                    3.7, 3.9, 4.25, 4.75, 5.25, 
                                                    5.75, 6.25, 6.75, 7.5, 8.5, 
                                                    9.5, 10.5, 11.5, 12.5, 13.5, 
                                                    14.5, 15.5, 17.0, 19.0, 21.0, 
                                                    23.0, 25.0, 27.0, 29.0, 31.0 };
    Double_t chHadALICE_PP2760GeV_xerrminus[65] = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
                                                    0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.05, 0.05, 
                                                    0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.10, 0.1, 
                                                    0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.25, 0.25, 
                                                    0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
                                                    0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
                                                    1.0 };
    Double_t chHadALICE_PP2760GeV_xerrplus[65]  = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
                                                    0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.05, 0.05, 
                                                    0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.10, 0.1, 
                                                    0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.25, 0.25, 
                                                    0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
                                                    0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
                                                    1.0 };
    Double_t chHadALICE_PP2760GeV_yval[65]      = { 357.6, 286.0, 223.9, 172.8, 133.2, 103.3, 80.54, 63.46, 50.38, 
                                                    40.29, 32.43, 26.33, 21.57, 17.65, 14.62, 12.15, 10.13, 7.799, 5.571, 
                                                    4.066, 2.984, 2.224, 1.662, 1.266, 0.9676, 0.7463, 0.5803, 0.4131, 0.2614, 
                                                    0.1667, 0.1104, 0.0749, 0.05148, 0.03607, 0.02516, 0.01807, 0.01326, 0.007935, 0.003978, 
                                                    0.002122, 0.001209, 7.046E-4, 4.346E-4, 2.269E-4, 1.008E-4, 4.915E-5, 2.735E-5, 1.488E-5, 8.995E-6, 
                                                    5.704E-6, 3.021E-6, 2.003E-6, 1.252E-6, 5.819E-7, 3.291E-7, 1.359E-7, 1.163E-7, 5.505E-8, 3.128E-8, 
                                                    2.73E-8 };
    Double_t chHadALICE_PP2760GeV_ysysminus[65] = { 27.2, 18.1, 15.1, 12.40, 10.10, 8.3, 6.39, 4.92, 3.81, 
                                                    2.89, 2.33, 1.89, 1.55, 1.27, 1.05, 0.87, 0.73, 0.56, 0.400, 
                                                    0.292, 0.214, 0.160, 0.119, 0.091, 0.0695, 0.0536, 0.0417, 0.0297, 0.0188, 
                                                    0.0120, 0.0079, 0.00538, 0.00370, 0.00259, 0.00181, 0.00130, 9.5E-4, 5.7E-4, 2.86E-4, 
                                                    1.52E-4, 8.7E-5, 5.06E-5, 3.12E-5, 1.63E-5, 7.2E-6, 3.53E-6, 2.07E-6, 1.13E-6, 6.84E-7, 
                                                    4.36E-7, 2.31E-7, 1.54E-7, 0.96E-7, 4.51E-8, 2.55E-8, 1.06E-8, 0.91E-8, 0.431E-8, 0.246E-8, 
                                                    2.16E-9 };
    Double_t chHadALICE_PP2760GeV_ysysplus[65]  = { 27.2, 18.1, 15.1, 12.40, 10.10, 8.3, 6.39, 4.92, 3.81, 
                                                    2.89, 2.33, 1.89, 1.55, 1.27, 1.05, 0.87, 0.73, 0.56, 0.400, 
                                                    0.292, 0.214, 0.160, 0.119, 0.091, 0.0695, 0.0536, 0.0417, 0.0297, 0.0188, 
                                                    0.0120, 0.0079, 0.00538, 0.00370, 0.00259, 0.00181, 0.00130, 9.5E-4, 5.7E-4, 2.86E-4, 
                                                    1.52E-4, 8.7E-5, 5.06E-5, 3.12E-5, 1.63E-5, 7.2E-6, 3.53E-6, 2.07E-6, 1.13E-6, 6.84E-7, 
                                                    4.36E-7, 2.31E-7, 1.54E-7, 0.96E-7, 4.51E-8, 2.55E-8, 1.06E-8, 0.91E-8, 0.431E-8, 0.246E-8, 
                                                    2.16E-9 };
    Double_t chHadALICE_PP2760GeV_yStatMin[65]  = { 0.5, 0.2, 0.2, 0.1, 0.1, 0.1, 0.07, 0.06, 0.05, 
                                                    0.04, 0.03, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.008, 0.006, 
                                                    0.005, 0.004, 0.004, 0.003, 0.003, 0.002, 0.0019, 0.0015, 0.001, 6.0E-4, 
                                                    5.0E-4, 4.0E-4, 4.2E-4, 2.0E-4, 1.5E-4, 1.1E-4, 8.0E-5, 7.0E-5, 4.3E-5, 2.4E-5, 
                                                    1.5E-5, 1.0E-5, 6.6E-6, 4.8E-6, 2.3E-6, 1.4E-6, 9.0E-7, 6.3E-7, 4.4E-7, 3.25E-7, 
                                                    2.48E-7, 1.74E-7, 1.36E-7, 7.3E-8, 4.69E-8, 3.35E-8, 2.05E-8, 1.82E-8, 1.203E-8, 8.68E-9, 
                                                    7.89E-9 };
    Double_t chHadALICE_PP2760GeV_ystatplus[65] = { 0.5, 0.2, 0.2, 0.1, 0.1, 0.1, 0.07, 0.06, 0.05, 
                                                    0.04, 0.03, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.008, 0.006, 
                                                    0.005, 0.004, 0.004, 0.003, 0.003, 0.002, 0.0019, 0.0015, 0.001, 6.0E-4, 
                                                    5.0E-4, 4.0E-4, 4.2E-4, 2.0E-4, 1.5E-4, 1.1E-4, 8.0E-5, 7.0E-5, 4.3E-5, 2.4E-5, 
                                                    1.5E-5, 1.0E-5, 6.6E-6, 4.8E-6, 2.3E-6, 1.4E-6, 9.0E-7, 6.3E-7, 4.4E-7, 3.25E-7, 
                                                    2.48E-7, 1.74E-7, 1.36E-7, 7.3E-8, 4.69E-8, 3.35E-8, 2.05E-8, 1.82E-8, 1.203E-8, 8.68E-9, 
                                                    7.89E-9 };
    int chHadALICE_PP2760GeV_numpoints          = 60;
    TGraphAsymmErrors* graphChHadALICESysPP2760GeV      = new TGraphAsymmErrors(chHadALICE_PP2760GeV_numpoints, chHadALICE_PP2760GeV_xval, chHadALICE_PP2760GeV_yval, chHadALICE_PP2760GeV_xerrminus,
                                                                                chHadALICE_PP2760GeV_xerrplus, chHadALICE_PP2760GeV_ysysminus, chHadALICE_PP2760GeV_ysysplus);
    TGraphAsymmErrors* graphChHadALICEStatPP2760GeV     = new TGraphAsymmErrors(chHadALICE_PP2760GeV_numpoints, chHadALICE_PP2760GeV_xval, chHadALICE_PP2760GeV_yval, chHadALICE_PP2760GeV_xerrminus,
                                                                                chHadALICE_PP2760GeV_xerrplus, chHadALICE_PP2760GeV_yStatMin, chHadALICE_PP2760GeV_ystatplus);

    
    // *********************************************************************************************************************
    // *********************************** ATLAS charged hadron results 2.76 TeV *******************************************
    // *********************************** arXiv: 1504.04337, JHEP 1509 (2015) 050 (2015-09-09) ****************************
    // *********************************************************************************************************************
    
    Double_t chHadATLAS_PP2760GeV_xval[37]          = { 0.5365, 0.615, 0.705, 0.808, 0.926, 1.0595, 1.21, 1.385, 1.59, 
                                                        1.825, 2.095, 2.405, 2.755, 3.155, 3.62, 4.15, 4.755, 5.455, 6.255, 
                                                        7.17, 8.20, 9.39, 10.75, 12.35, 14.15, 16.2, 18.6, 21.35, 24.45, 
                                                        28.05, 33.85, 42.6, 53.65, 67.55, 85.05, 106.9, 134.5 };
    Double_t chHadATLAS_PP2760GeV_xerrminus[37]     = { 0.0365, 0.042, 0.048, 0.055, 0.063, 0.0705, 0.080, 0.095, 0.11, 
                                                        0.125, 0.145, 0.165, 0.185, 0.215, 0.25, 0.28, 0.325, 0.375, 0.425, 
                                                        0.490, 0.56, 0.610, 0.75, 0.85, 0.95, 1.1, 1.3, 1.45, 1.65, 
                                                        1.95, 3.85, 4.90, 6.15, 7.75, 9.75, 12.1, 15.5 };
    Double_t chHadATLAS_PP2760GeV_xerrplus[37]      = { 0.0365, 0.042, 0.048, 0.055, 0.063, 0.0705, 0.080, 0.095, 0.11, 
                                                        0.125, 0.145, 0.165, 0.185, 0.215, 0.25, 0.28, 0.325, 0.375, 0.425, 
                                                        0.490, 0.56, 0.610, 0.75, 0.85, 0.95, 1.10, 1.3, 1.45, 1.65, 
                                                        1.95, 3.85, 4.90, 6.15, 7.75, 9.75, 12.1, 15.5 };
    Double_t chHadATLAS_PP2760GeV_yval[37]          = { 62.5, 43.7, 29.6, 19.5, 12.4, 7.76, 4.71, 2.74, 1.52, 
                                                        0.812, 0.416, 0.206, 0.0987, 0.0457, 0.0204, 0.00895, 0.00385, 0.00162, 6.74E-4, 
                                                        2.8E-4, 1.16E-4, 4.96E-5, 2.11E-5, 8.58E-6, 3.59E-6, 1.51E-6, 6.22E-7, 2.55E-7, 1.06E-7, 
                                                        4.54E-8, 1.31E-8, 2.86E-9, 5.84E-10, 1.17E-10, 2.29E-11, 4.06E-12, 6.86E-13 };
    Double_t chHadATLAS_PP2760GeV_yerrminus[37]     = { 2.9375265956242846, 2.0976200745914406, 1.4208149232549607, 0.955511461155752, 0.5952087316026202, 
                                                        0.372486798021299, 0.226085023944179, 0.13152369894798427, 0.07448273594974879, 
                                                        0.0389781145259747, 0.01996961236823589, 0.009683232634301419, 0.004738527791965032, 0.0021942853929240837, 
                                                        9.593554135981097E-4, 4.121287930065066E-4, 1.8511414964826432E-4, 7.476067439503204E-5, 3.1871292836030354E-5, 
                                                        1.3305156894978728E-5, 5.448605326136221E-6, 2.2756618377957655E-6, 9.910825205299507E-7, 4.058076968220292E-7, 
                                                        1.6630653204249072E-7, 6.995065832427884E-8, 2.913457561043236E-8, 1.1938793908934018E-8, 5.198324730141433E-9, 
                                                        2.0675700326712032E-9, 5.814393347547101E-10, 1.297129029819316E-10, 2.8890238583992348E-11, 6.306115488951974E-12,
                                                        1.34974084179149E-12, 2.8576170912142867E-13, 6.240715100050635E-14 };
    Double_t chHadATLAS_PP2760GeV_yerrplus[37]      = { 2.9375265956242846, 2.0976200745914406, 1.4208149232549607, 0.955511461155752, 0.5952087316026202,
                                                        0.372486798021299, 0.226085023944179, 0.13152369894798427, 0.07448273594974879, 
                                                        0.0389781145259747, 0.01996961236823589, 0.009683232634301419, 0.004738527791965032, 0.0021942853929240837, 
                                                        9.593554135981097E-4, 4.121287930065066E-4, 1.8511414964826432E-4, 7.476067439503204E-5, 3.1871292836030354E-5, 
                                                        1.3305156894978728E-5, 5.448605326136221E-6, 2.2756618377957655E-6, 9.910825205299507E-7, 4.058076968220292E-7, 
                                                        1.6630653204249072E-7, 6.995065832427884E-8, 2.913457561043236E-8, 1.1938793908934018E-8, 5.198324730141433E-9, 
                                                        2.0675700326712032E-9, 5.814393347547101E-10, 1.297129029819316E-10, 2.8890238583992348E-11, 6.306115488951974E-12, 
                                                        1.34974084179149E-12, 2.8576170912142867E-13, 6.240715100050635E-14 };
    Double_t chHadATLAS_PP2760GeV_ystatminus[37]    = { 0.0125, 0.009177000000000001, 0.006512, 0.00468, 0.003224, 
                                                        0.0022504, 0.0015072000000000002, 9.864E-4, 6.384E-4, 
                                                        4.0600000000000006E-4, 2.5375999999999996E-4, 1.545E-4, 9.376500000000001E-5, 5.484E-5, 
                                                        3.2640000000000006E-5, 1.8795E-5, 1.0780000000000002E-5, 5.9940000000000005E-6, 3.5048E-6, 
                                                        1.96E-6, 1.1019999999999998E-6, 6.448000000000001E-7, 2.0045000000000002E-7, 9.437999999999999E-8, 
                                                        3.949E-8, 1.661E-8, 8.086E-9, 4.08E-9, 1.802E-9, 
                                                        6.81E-10, 1.441E-10, 3.146E-11, 3.971200000000001E-12, 1.1466E-12, 
                                                        3.435E-13, 1.0962E-13, 3.43E-14 };
    Double_t chHadATLAS_PP2760GeV_ystatplus[37]     = { 0.0125, 0.009177000000000001, 0.006512, 0.00468, 0.003224, 
                                                        0.0022504, 0.0015072000000000002, 9.864E-4, 6.384E-4, 
                                                        4.0600000000000006E-4, 2.5375999999999996E-4, 1.545E-4, 9.376500000000001E-5, 5.484E-5,
                                                        3.2640000000000006E-5, 1.8795E-5, 1.0780000000000002E-5, 5.9940000000000005E-6, 3.5048E-6, 
                                                        1.96E-6, 1.1019999999999998E-6, 6.448000000000001E-7, 2.0045000000000002E-7, 9.437999999999999E-8,
                                                        3.949E-8, 1.661E-8, 8.086E-9, 4.08E-9, 1.802E-9, 
                                                        6.81E-10, 1.441E-10, 3.146E-11, 3.971200000000001E-12, 1.1466E-12,
                                                        3.435E-13, 1.0962E-13, 3.43E-14 };
    Int_t chHadATLAS_PP2760GeV_numpoints            = 37;
    Double_t chHadATLAS_PP2760GeV_ysysminus[37];
    Double_t chHadATLAS_PP2760GeV_ysysplus[37];
    ExtractSystematicFromTotal(chHadATLAS_PP2760GeV_numpoints, chHadATLAS_PP2760GeV_yerrminus, chHadATLAS_PP2760GeV_ystatminus, chHadATLAS_PP2760GeV_ysysminus );
    ExtractSystematicFromTotal(chHadATLAS_PP2760GeV_numpoints, chHadATLAS_PP2760GeV_yerrplus, chHadATLAS_PP2760GeV_ystatplus, chHadATLAS_PP2760GeV_ysysplus );
    
    TGraphAsymmErrors* graphChHadATLASStatPP2760GeV = new TGraphAsymmErrors(chHadATLAS_PP2760GeV_numpoints, chHadATLAS_PP2760GeV_xval, chHadATLAS_PP2760GeV_yval, chHadATLAS_PP2760GeV_xerrminus,
                                                                            chHadATLAS_PP2760GeV_xerrplus, chHadATLAS_PP2760GeV_ystatminus, chHadATLAS_PP2760GeV_ystatplus);    
    TGraphAsymmErrors* graphChHadATLASSysPP2760GeV  = new TGraphAsymmErrors(chHadATLAS_PP2760GeV_numpoints, chHadATLAS_PP2760GeV_xval, chHadATLAS_PP2760GeV_yval, chHadATLAS_PP2760GeV_xerrminus,
                                                                            chHadATLAS_PP2760GeV_xerrplus, chHadATLAS_PP2760GeV_ysysminus, chHadATLAS_PP2760GeV_ysysplus);
    
    
    // *********************************************************************************************************************
    // *********************************** CMS charged hadron results 2.76 TeV *********************************************
    // *********************************** arXiv:1202.2554, Eur.Phys.J. C72 (2012) 1945 ************************************
    // *********************************************************************************************************************
    Double_t chHadCMS_PP2760GeV_xval[22]            = { 0.525, 0.675, 0.825, 0.975, 1.125, 1.35, 1.65, 1.95, 2.25, 
                                                        3.0, 4.2, 5.4, 6.6, 9.0, 12.6, 18.0, 25.2, 33.6, 43.2, 
                                                        57.6, 76.8, 99.3 };
    Double_t chHadCMS_PP2760GeV_xerrminus[22]       = { 0.075, 0.075, 0.075, 0.075, 0.075, 0.15, 0.15, 0.15, 0.15, 
                                                        0.60, 0.60, 0.60, 0.60, 1.80, 1.80, 3.6, 3.6, 4.80, 4.80, 
                                                        9.60, 9.60, 12.90 };
    Double_t chHadCMS_PP2760GeV_xerrplus[22]        = { 0.075, 0.075, 0.075, 0.075, 0.075, 0.15, 0.15, 0.15, 0.15, 
                                                        0.60, 0.60, 0.60, 0.60, 1.80, 1.80, 3.60, 3.60, 4.80, 4.80, 
                                                        9.60, 9.60, 12.90 };
    Double_t chHadCMS_PP2760GeV_yval[22]            = { 1.2569, 0.6522, 0.3548, 0.20427, 0.12288, 0.06077, 0.02604, 0.01205, 0.005967, 
                                                        0.0012795, 1.7414E-4, 3.694E-5, 9.993E-6, 1.3496E-6, 1.4955E-7, 1.4511E-8, 1.7416E-9, 2.796E-10, 4.739E-11, 
                                                        6.798E-12, 1.0808E-12, 1.868E-13 };
    Double_t chHadCMS_PP2760GeV_yerrminus[22]       = { 0.058900339557595084, 0.030600163398256552, 0.016600301202086665, 0.009580255737713895, 0.005760217009800933, 
                                                        0.0028500701745746543, 0.0012200409829181968, 5.650566343296926E-4, 2.800446392988089E-4, 
                                                        6.000833275470999E-5, 8.175879157619687E-6, 1.734877517290486E-6, 4.729503145151719E-7, 6.415294225520759E-8, 
                                                        7.55079466016657E-9, 7.942474425517529E-10, 1.2467630087550723E-10, 2.3940133667128928E-11, 3.5638742963241563E-12, 
                                                        5.975717530138117E-13, 1.1667064755113002E-13, 2.920719089539423E-14 };
    Double_t chHadCMS_PP2760GeV_yerrplus[22]        = { 0.058900339557595084, 0.030600163398256552, 0.016600301202086665, 0.009580255737713895, 0.005760217009800933, 
                                                        0.0028500701745746543, 0.0012200409829181968, 5.650566343296926E-4, 2.800446392988089E-4, 
                                                        6.000833275470999E-5, 8.175879157619687E-6, 1.734877517290486E-6, 4.729503145151719E-7, 6.415294225520759E-8, 
                                                        7.55079466016657E-9, 7.942474425517529E-10, 1.2467630087550723E-10, 2.3940133667128928E-11, 3.5638742963241563E-12, 
                                                        5.975717530138117E-13, 1.1667064755113002E-13, 2.920719089539423E-14 };
    Double_t chHadCMS_PP2760GeV_ystatminus[22]      = { 2.0E-4, 1.0E-4, 1.0E-4, 7.0E-5, 5.0E-5, 2.0E-5, 1.0E-5, 8.0E-6, 5.0E-6, 
                                                        1.0E-6, 3.1E-7, 1.3E-7, 6.1E-8, 9.8E-9, 2.73E-9, 4.02E-10, 9.33E-11, 1.92E-11, 1.44E-12, 
                                                        2.24E-13, 7.4E-14, 2.45E-14 };
    Double_t chHadCMS_PP2760GeV_ystatplus[22]       = { 2.0E-4, 1.0E-4, 1.0E-4, 7.0E-5, 5.0E-5, 2.0E-5, 1.0E-5, 8.0E-6, 5.0E-6, 
                                                        1.0E-6, 3.1E-7, 1.3E-7, 6.1E-8, 9.8E-9, 2.73E-9, 4.02E-10, 9.33E-11, 1.92E-11, 1.44E-12, 
                                                        2.24E-13, 7.4E-14, 2.45E-14 };
    int chHadCMS_PP2760GeV_numpoints                = 22;
    Double_t chHadCMS_PP2760GeV_ysysminus[22];
    Double_t chHadCMS_PP2760GeV_ysysplus[22];
    ExtractSystematicFromTotal(chHadCMS_PP2760GeV_numpoints, chHadCMS_PP2760GeV_yerrminus, chHadCMS_PP2760GeV_ystatminus, chHadCMS_PP2760GeV_ysysminus );
    ExtractSystematicFromTotal(chHadCMS_PP2760GeV_numpoints, chHadCMS_PP2760GeV_yerrplus, chHadCMS_PP2760GeV_ystatplus, chHadCMS_PP2760GeV_ysysplus );
    
    TGraphAsymmErrors* graphChHadYieldCMSStatPP2760GeV  = new TGraphAsymmErrors(chHadCMS_PP2760GeV_numpoints, chHadCMS_PP2760GeV_xval, chHadCMS_PP2760GeV_yval, chHadCMS_PP2760GeV_xerrminus,
                                                                                chHadCMS_PP2760GeV_xerrplus, chHadCMS_PP2760GeV_ystatminus, chHadCMS_PP2760GeV_ystatplus);    
    TGraphAsymmErrors* graphChHadYieldCMSSysPP2760GeV   = new TGraphAsymmErrors(chHadCMS_PP2760GeV_numpoints, chHadCMS_PP2760GeV_xval, chHadCMS_PP2760GeV_yval, chHadCMS_PP2760GeV_xerrminus,
                                                                                chHadCMS_PP2760GeV_xerrplus, chHadCMS_PP2760GeV_ysysminus, chHadCMS_PP2760GeV_ysysplus);
    TGraphAsymmErrors* graphChHadCMSStatPP2760GeV       = ScaleGraph(graphChHadYieldCMSStatPP2760GeV, 64.);
    TGraphAsymmErrors* graphChHadCMSSysPP2760GeV        = ScaleGraph(graphChHadYieldCMSSysPP2760GeV, 64.);
    
    
    // *********************************************************************************************************************
    // ************************************ charged particle density dN/dy *************************************************
    // ************************************ http://arxiv.org/abs/1509.07541 ************************************************
    // ************************************ pp 900 GeV, 2.76 TeV, 7 TeV, 8 TeV *********************************************
    // *********************************************************************************************************************
    Double_t chHad_dNdeta_xValHist[21]              = { -2.0,   -1.8,   -1.6,   -1.4,   -1.2,
                                                        -1.0,   -0.8,   -0.6,   -0.4,   -0.2,
                                                        0.0,    0.2,    0.4,    0.6,    0.8,
                                                        1.0,    1.2,    1.4,    1.6,    1.8,
                                                        2.0 };    
    Double_t chHad_dNdeta_xVal[20]                  = { -1.9,   -1.7,   -1.5,   -1.3,   -1.1, 
                                                        -0.9,   -0.7,   -0.5,   -0.3,   -0.1,
                                                        0.1,    0.3,    0.5,    0.7,    0.9,
                                                        1.1,    1.3,    1.5,    1.7,    1.9};
    Double_t chHad_dNdeta_xErr[20]                  = { 0.1,    0.1,    0.1,    0.1,    0.1, 
                                                        0.1,    0.1,    0.1,    0.1,    0.1,
                                                        0.1,    0.1,    0.1,    0.1,    0.1, 
                                                        0.1,    0.1,    0.1,    0.1,    0.1};
    // 900 GeV values
    Double_t chHad_dNdeta_pp900GeVINEL_yVal[20]     = { 3.08, 3.16, 3.15, 3.16, 3.13,
                                                        3.08, 3.02, 2.97, 2.93, 2.90,
                                                        2.91, 2.95, 2.99, 3.02, 3.10,
                                                        3.14, 3.15, 3.17, 3.16, 3.07 };                        
    Double_t chHad_dNdeta_pp900GeVINEL_yErrUp[20]   = { 0.14, 0.14, 0.13, 0.13, 0.12,
                                                        0.12, 0.12, 0.11, 0.11, 0.11,
                                                        0.11, 0.11, 0.11, 0.12, 0.12,
                                                        0.12, 0.13, 0.13, 0.13, 0.13};
    Double_t chHad_dNdeta_pp900GeVINEL_yErrDown[20] = { 0.06, 0.06, 0.06, 0.06, 0.05,
                                                        0.05, 0.05, 0.05, 0.05, 0.05,
                                                        0.05, 0.05, 0.05, 0.05, 0.05,
                                                        0.05, 0.06, 0.06, 0.06, 0.06 };

    TGraphAsymmErrors* graphChHaddNdyINELTotPP900GeV    = new TGraphAsymmErrors(20, chHad_dNdeta_xVal, chHad_dNdeta_pp900GeVINEL_yVal, chHad_dNdeta_xErr,
                                                                                chHad_dNdeta_xErr, chHad_dNdeta_pp900GeVINEL_yErrDown, chHad_dNdeta_pp900GeVINEL_yErrUp);    

    TH1D* histoChHaddNdyINELTotPP900GeV                 = new TH1D("histoChargedHadrondNdEtaALICEPP900GeV", "histoChargedHadrondNdEtaALICEPP900GeV", 20, chHad_dNdeta_xValHist);
    for (Int_t i = 1; i< 20+1; i++){
        histoChHaddNdyINELTotPP900GeV->SetBinContent(i,chHad_dNdeta_pp900GeVINEL_yVal[i-1]);
        histoChHaddNdyINELTotPP900GeV->SetBinError(i,(chHad_dNdeta_pp900GeVINEL_yErrUp[i-1]+chHad_dNdeta_pp900GeVINEL_yErrDown[i-1])/2.);
    }
    
    // 2760 GeV values
    Double_t chHad_dNdeta_pp2760GeVINEL_yVal[20]    = { 4.00, 4.08, 4.08, 4.06, 4.01,
                                                        3.94, 3.86, 3.80, 3.72, 3.71,
                                                        3.71, 3.74, 3.81, 3.87, 3.95,
                                                        4.02, 4.06, 4.08, 4.09, 4.02 };                        
    Double_t chHad_dNdeta_pp2760GeVINEL_yErrUp[20]  = { 0.29, 0.29, 0.29, 0.29, 0.28,
                                                        0.27, 0.27, 0.26, 0.26, 0.26, 
                                                        0.26, 0.26, 0.26, 0.27, 0.27,
                                                        0.28, 0.28, 0.29, 0.29, 0.29 };
    Double_t chHad_dNdeta_pp2760GeVINEL_yErrDown[20]= { 0.17, 0.18, 0.17, 0.17, 0.17,
                                                        0.17, 0.16, 0.16, 0.16, 0.16,
                                                        0.16, 0.16, 0.16, 0.16, 0.17,
                                                        0.17, 0.17, 0.17, 0.18, 0.17 };
    TGraphAsymmErrors* graphChHaddNdyINELTotPP2760GeV   = new TGraphAsymmErrors(20, chHad_dNdeta_xVal, chHad_dNdeta_pp2760GeVINEL_yVal, chHad_dNdeta_xErr,
                                                                                chHad_dNdeta_xErr, chHad_dNdeta_pp2760GeVINEL_yErrDown, chHad_dNdeta_pp2760GeVINEL_yErrUp);    
    
    TH1D* histoChHaddNdyINELTotPP2760GeV                = new TH1D("histoChargedHadrondNdEtaALICEPP2760GeV", "histoChargedHadrondNdEtaALICEPP2760GeV", 20, chHad_dNdeta_xValHist);
    for (Int_t i = 1; i< 20+1; i++){
        histoChHaddNdyINELTotPP2760GeV->SetBinContent(i,chHad_dNdeta_pp2760GeVINEL_yVal[i-1]);
        histoChHaddNdyINELTotPP2760GeV->SetBinError(i,(chHad_dNdeta_pp2760GeVINEL_yErrUp[i-1]+chHad_dNdeta_pp2760GeVINEL_yErrDown[i-1])/2.);
    }

    // 7 TeV values
    Double_t chHad_dNdeta_pp7TeVINEL_yVal[20]       = { 4.94, 5.00, 5.01, 4.98, 4.92,
                                                        4.84, 4.74, 4.66, 4.59, 4.55,
                                                        4.55, 4.59, 4.66, 4.74, 4.84,
                                                        4.91, 4.97, 5.00, 4.99, 4.94 };
    Double_t chHad_dNdeta_pp7TeVINEL_yErrUp[20]     = { 0.38, 0.38, 0.38, 0.37, 0.37,
                                                        0.36, 0.35, 0.35, 0.34, 0.34,
                                                        0.34, 0.34, 0.35, 0.35, 0.36,
                                                        0.37, 0.37, 0.38, 0.38, 0.38 };
    Double_t chHad_dNdeta_pp7TeVINEL_yErrDown[20]   = { 0.19, 0.19, 0.19, 0.19, 0.18,
                                                        0.18, 0.18, 0.17, 0.17, 0.17, 
                                                        0.17, 0.17, 0.17, 0.18, 0.18,
                                                        0.18, 0.19, 0.19, 0.19, 0.20 };
    TGraphAsymmErrors* graphChHaddNdyINELTotPP7TeV      = new TGraphAsymmErrors(20, chHad_dNdeta_xVal, chHad_dNdeta_pp7TeVINEL_yVal, chHad_dNdeta_xErr,
                                                                                chHad_dNdeta_xErr, chHad_dNdeta_pp7TeVINEL_yErrDown, chHad_dNdeta_pp7TeVINEL_yErrUp);    
    TH1D* histoChHaddNdyINELTotPP7TeV                   = new TH1D("histoChargedHadrondNdEtaALICEPP7TeV", "histoChargedHadrondNdEtaALICEPP7TeV", 20, chHad_dNdeta_xValHist);
    for (Int_t i = 1; i< 20+1; i++){
        histoChHaddNdyINELTotPP7TeV->SetBinContent(i,chHad_dNdeta_pp7TeVINEL_yVal[i-1]);
        histoChHaddNdyINELTotPP7TeV->SetBinError(i,(chHad_dNdeta_pp7TeVINEL_yErrUp[i-1]+chHad_dNdeta_pp7TeVINEL_yErrDown[i-1])/2.);
    }
    
    // 8 TeV values
    Double_t chHad_dNdeta_pp8TeVINEL_yVal[20]       = { 6.16, 6.23, 6.25, 6.21, 6.14,
                                                        6.04, 5.92, 5.82, 5.73, 5.68,
                                                        5.68, 5.73, 5.82, 5.91, 6.03,
                                                        6.13, 6.20, 6.23, 6.22, 6.16 };
    Double_t chHad_dNdeta_pp8TeVINEL_yErrUp[20]     = { 0.22, 0.20, 0.19, 0.18, 0.17,
                                                        0.16, 0.16, 0.15, 0.15, 0.15, 
                                                        0.15, 0.15, 0.15, 0.16, 0.16,
                                                        0.17, 0.17, 0.19, 0.20, 0.22 };
    Double_t chHad_dNdeta_pp8TeVINEL_yErrDown[20]   = { 0.16, 0.15, 0.15, 0.14, 0.14,
                                                        0.14, 0.13, 0.13, 0.13, 0.13,
                                                        0.13, 0.13, 0.13, 0.14, 0.14,
                                                        0.14, 0.14, 0.15, 0.15, 0.17 };
    TGraphAsymmErrors* graphChHaddNdyINELTotPP8TeV      = new TGraphAsymmErrors(20, chHad_dNdeta_xVal, chHad_dNdeta_pp8TeVINEL_yVal, chHad_dNdeta_xErr,
                                                                                chHad_dNdeta_xErr, chHad_dNdeta_pp8TeVINEL_yErrDown, chHad_dNdeta_pp8TeVINEL_yErrUp);    
    TH1D* histoChHaddNdyINELTotPP8TeV                   = new TH1D("histoChargedHadrondNdEtaALICEPP8TeV", "histoChargedHadrondNdEtaALICEPP8TeV", 20, chHad_dNdeta_xValHist);
    for (Int_t i = 1; i< 20+1; i++){
        histoChHaddNdyINELTotPP8TeV->SetBinContent(i,chHad_dNdeta_pp8TeVINEL_yVal[i-1]);
        histoChHaddNdyINELTotPP8TeV->SetBinError(i,(chHad_dNdeta_pp8TeVINEL_yErrUp[i-1]+chHad_dNdeta_pp8TeVINEL_yErrDown[i-1])/2.);
    }

    // *********************************************************************************************************************
    // ********************************** Write Output files ***************************************************************
    // *********************************************************************************************************************
    TFile fileChargedHadronspp(Form("ExternalInput/UnidentifiedCharged/ChargedHadronSpectraPP_%s.root",dateForOutput.Data()) ,"RECREATE");
        graphChHadALICEStatPP2760GeV->Write("graphChargedHadronsALICEStatPP2760GeV");
        graphChHadALICESysPP2760GeV->Write("graphChargedHadronsALICESysPP2760GeV");
        graphChHadATLASStatPP2760GeV->Write("graphChargedHadronsATLASStatPP2760GeV");
        graphChHadATLASSysPP2760GeV->Write("graphChargedHadronsATLASSysPP2760GeV");
        graphChHadCMSStatPP2760GeV->Write("graphChargedHadronsCMSStatPP2760GeV");
        graphChHadCMSSysPP2760GeV->Write("graphChargedHadronsCMSSysPP2760GeV");
        graphChHadALICEStatPP900GeV->Write("graphChargedHadronsALICEStatPP900GeV");
        graphChHadALICESysPP900GeV->Write("graphChargedHadronsALICESysPP900GeV");
        graphChHadCMSStatPP900GeV->Write("graphChargedHadronsCMSStatPP900GeV");
        graphChHadCMSSysPP900GeV->Write("graphChargedHadronsCMSSysPP900GeV");
        graphChHadALICEStatPP7TeV->Write("graphChargedHadronsALICEStatPP7TeV");
        graphChHadALICESysPP7TeV->Write("graphChargedHadronsALICESysPP7TeV");
        graphChHadCMSStatPP7TeV->Write("graphChargedHadronsCMSStatPP7TeV");
        graphChHadCMSSysPP7TeV->Write("graphChargedHadronsCMSSysPP7TeV");
        
        graphChHaddNdyINELTotPP900GeV->Write("graphChargedHadrondNdEtaALICEPP900GeV");
        histoChHaddNdyINELTotPP900GeV->Write("histoChargedHadrondNdEtaALICEPP900GeV");
        graphChHaddNdyINELTotPP2760GeV->Write("graphChargedHadrondNdEtaALICEPP2760GeV");
        histoChHaddNdyINELTotPP2760GeV->Write("histoChargedHadrondNdEtaALICEPP2760GeV");
        graphChHaddNdyINELTotPP7TeV->Write("graphChargedHadrondNdEtaALICEPP7TeV");
        histoChHaddNdyINELTotPP7TeV->Write("histoChargedHadrondNdEtaALICEPP7TeV");
        graphChHaddNdyINELTotPP8TeV->Write("graphChargedHadrondNdEtaALICEPP8TeV");
        histoChHaddNdyINELTotPP8TeV->Write("histoChargedHadrondNdEtaALICEPP8TeV");
        
        
    fileChargedHadronspp.Close();

    TFile fileChargedPionspp(Form("ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_%s.root",dateForOutput.Data()) ,"RECREATE");
        histoChargedPionSpecPubStatPP2760GeV->Write("histoChargedPionSpecPubStat2760GeV");
        histoChargedPionSpecPubSystPP2760GeV->Write("histoChargedPionSpecPubSyst2760GeV");
        histoChargedKaonSpecPubStatPP2760GeV->Write("histoChargedKaonSpecPubStat2760GeV");
        histoChargedKaonSpecPubSystPP2760GeV->Write("histoChargedKaonSpecPubSyst2760GeV");
        histoProtonSpecPubStatPP2760GeV->Write("histoProtonSpecPubStat2760GeV");
        histoProtonSpecPubSystPP2760GeV->Write("histoProtonSpecPubSyst2760GeV");
        histoNeutralKaonSpecStatPP2760GeV->Write("histoNeutralKaonSpecStat2760GeV");
        histoNeutralKaonSpecSysPP2760GeV->Write("histoNeutralKaonSpecSyst2760GeV");
        histoLambda1115SpecStatPP2760GeV->Write("histoLambda1115SpecStat2760GeV");
        histoLambda1115SpecSystPP2760GeV->Write("histoLambda1115SpecSyst2760GeV");
        histoAntiLambda1115SpecStatPP2760GeV->Write("histoAntiLambda1115SpecStat2760GeV");
        histoAntiLambda1115SpecSystPP2760GeV->Write("histoAntiLambda1115SpecSyst2760GeV");

        histoChargedPionSpecHighPtStatPP2760GeV->Write("histoChargedPionSpecHighPtStat2760GeV");
        histoChargedPionSpecHighPtSystPP2760GeV->Write("histoChargedPionSpecHighPtSyst2760GeV");
        
        histoChargedPionSpecLowPtStat2760GeVCMS->Write("histoChargedPionSpecLowPtStat2760GeVCMS");
        histoChargedPionSpecLowPtSys2760GeVCMS->Write("histoChargedPionSpecLowPtSys2760GeVCMS");
        histoChargedPionSpecLowPtStatPP2760GeV->Write("histoChargedPionSpecLowPtStatPP2760GeV");
        histoChargedPionSpecLowPtSysPP2760GeV->Write("histoChargedPionSpecLowPtSysPP2760GeV");
        
        histoChargedPionSpecLowPtStat900GeVCMS->Write("histoChargedPionSpecLowPtStat900GeVCMS");
        histoChargedPionSpecLowPtSys900GeVCMS->Write("histoChargedPionSpecLowPtSys900GeVCMS");
        histoChargedPionSpecLowPtStat900GeVALICE->Write("histoChargedPionSpecLowPtStat900GeVALICE");
        histoChargedPionSpecLowPtSys900GeVALICE->Write("histoChargedPionSpecLowPtSys900GeVALICE");
        graphChKaonMinALICEStatPP900GeV->Write("graphNegKaonSpecLowPtStat900GeV");
        graphChKaonMinALICESysPP900GeV->Write("graphNegKaonSpecLowPtSys900GeV");
        graphChKaonPlALICEStatPP900GeV->Write("graphPosKaonSpecLowPtStat900GeV");
        graphChKaonPlALICESysPP900GeV->Write("graphPosKaonSpecLowPtSys900GeV");
        
        histoChargedPionSpecLowPtStat7TeVCMS->Write("histoChargedPionSpecLowPtStat7TeVCMS");
        histoChargedPionSpecLowPtSys7TeVCMS->Write("histoChargedPionSpecLowPtSys7TeVCMS");
        histoChargedPionSpecHighPtStat7TeVALICE->Write("histoChargedPionSpecHighPtStat7TeVALICE");
        graphChargedPionSpecHighPtSys7TeVALICE->Write("graphChargedPionSpecHighPtSys7TeVALICE");

        histoChargedPionSpecLowPtSys7TeVALICE->Write("histoChargedPionSpecLowPtSys7TeVALICE");
        histoChargedPionSpecLowPtStat7TeVALICE->Write("histoChargedPionSpecLowPtStat7TeVALICE");
        
        histoChargedPionSpecPubStatPP7TeV->Write("histoChargedPionSpecPubStat7TeV");
        histoChargedPionSpecPubSystPP7TeV->Write("histoChargedPionSpecPubSyst7TeV");
        histoChargedKaonSpecPubStatPP7TeV->Write("histoChargedKaonSpecPubStat7TeV");
        histoChargedKaonSpecPubSystPP7TeV->Write("histoChargedKaonSpecPubSyst7TeV");
        histoProtonSpecPubStatPP7TeV->Write("histoProtonSpecPubStat7TeV");
        histoProtonSpecPubSystPP7TeV->Write("histoProtonSpecPubSyst7TeV");

        graphChHaddNdyINELTotPP900GeV->Write("graphChargedHadrondNdEtaALICEPP900GeV");
        histoChHaddNdyINELTotPP900GeV->Write("histoChargedHadrondNdEtaALICEPP900GeV");
        graphChHaddNdyINELTotPP2760GeV->Write("graphChargedHadrondNdEtaALICEPP2760GeV");
        histoChHaddNdyINELTotPP2760GeV->Write("histoChargedHadrondNdEtaALICEPP2760GeV");
        graphChHaddNdyINELTotPP7TeV->Write("graphChargedHadrondNdEtaALICEPP7TeV");
        histoChHaddNdyINELTotPP7TeV->Write("histoChargedHadrondNdEtaALICEPP7TeV");
        graphChHaddNdyINELTotPP8TeV->Write("graphChargedHadrondNdEtaALICEPP8TeV");
        histoChHaddNdyINELTotPP8TeV->Write("histoChargedHadrondNdEtaALICEPP8TeV");
        
    fileChargedPionspp.Close();
    
}

