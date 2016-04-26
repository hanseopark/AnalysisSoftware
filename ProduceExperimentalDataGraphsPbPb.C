/****************************************************************************************************************************
****** 		provided by Gamma Conversion Group, PWG4, 													*****
******		Ana Marin, marin@physi.uni-heidelberg.de													*****
******	   	Kathrin Koch, kkoch@physi.uni-heidelberg.de 													*****
******		Friederike Bock, friederike.bock@cern.ch													*****
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

extern TRandom*	gRandom;
extern TBenchmark*	gBenchmark;
extern TSystem*	gSystem;
extern TMinuit*  	gMinuit;

void ProduceExperimentalDataGraphsPbPb(){	
    
    //********************************************************************************************************************************************	
    //*************************************** Raa theory Xiao-Fang *******************************************************************************
    
    //PHENIX 200 GeV
    //pi0 points from http://arxiv.org/pdf/0801.4020v1.pdf
    //eta points from http://arxiv.org/pdf/1005.4916v1.pdf
    Double_t PHENIX200GeVpt_0010[23] =  { 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75, 11, 13, 15, 17, 19};
    Double_t PHENIX200GeVRAA_0010[23] =  {3.640e-01, 4.123e-01, 3.958e-01, 3.159e-01, 2.556e-01, 2.272e-01, 2.103e-01, 2.051e-01, 1.855e-01, 1.854e-01, 1.911e-01, 2.070e-01, 1.964e-01, 1.875e-01, 1.856e-01, 1.901e-01, 2.004e-01, 2.200e-01, 2.099e-01, 1.663e-01, 2.892e-01, 3.739e-01, 3.449e-01};
    Double_t PHENIX200GeVRAAErr_0010[23] = {3.572e-02, 4.197e-02, 4.287e-02, 3.558e-02, 2.952e-02, 2.641e-02, 2.547e-02, 2.483e-02, 2.240e-02, 2.232e-02, 2.299e-02, 2.482e-02, 2.355e-02, 2.251e-02, 2.237e-02, 2.302e-02, 2.452e-02, 2.719e-02, 2.812e-02, 2.704e-02, 5.437e-02, 8.636e-02, 1.534e-01};
    Double_t PHENIX200GeVRAAErrPt_0010[23] = {5.983e-03, 5.426e-03, 5.159e-03, 4.425e-03, 4.123e-03, 4.208e-03, 3.988e-03, 4.284e-03, 4.507e-03, 5.381e-03, 6.178e-03, 5.652e-03, 6.401e-03, 7.132e-03, 8.370e-03, 9.373e-03, 1.170e-02, 1.542e-02, 1.134e-02, 2.028e-02, 4.825e-02, 9.209e-02, 1.869e-01};
    Double_t PHENIX200GeVNormError_0010 = TMath::Sqrt(6.86*6.86+9.7*9.7);

    Double_t PHENIX200GeVRAA_1020[23] = {4.221e-01, 4.551e-01, 4.435e-01, 3.925e-01, 3.345e-01, 3.095e-01, 2.822e-01, 2.694e-01, 2.637e-01, 2.447e-01, 2.528e-01, 2.685e-01, 2.573e-01, 2.589e-01, 2.686e-01, 2.406e-01, 2.842e-01, 2.535e-01, 2.647e-01, 2.670e-01, 4.705e-01, 3.253e-01, 6.186e-01};
    Double_t PHENIX200GeVRAAErr_1020[23] = {4.290e-02, 4.667e-02, 4.811e-02, 4.430e-02, 3.874e-02, 3.612e-02, 3.432e-02, 3.276e-02, 3.203e-02, 2.965e-02, 3.056e-02, 3.243e-02, 3.105e-02, 3.128e-02, 3.261e-02, 2.935e-02, 3.511e-02, 3.194e-02, 3.589e-02, 4.315e-02, 8.998e-02, 8.087e-02, 2.753e-01};
    Double_t PHENIX200GeVRAAErrPt_1020[23] = {6.757e-03, 5.178e-03, 4.882e-03, 4.566e-03, 4.516e-03, 4.817e-03, 4.341e-03, 4.975e-03, 5.649e-03, 6.429e-03, 7.736e-03, 6.831e-03, 7.752e-03, 9.609e-03, 1.142e-02, 1.266e-02, 1.635e-02, 1.818e-02, 1.487e-02, 2.772e-02, 7.497e-02, 9.467e-02, 3.054e-01};
    Double_t PHENIX200GeVNormError_1020 = TMath::Sqrt(6.97*6.97+9.7*9.7);

    Double_t PHENIX200GeVRAA_0020[23] = {3.895e-01,4.347e-01, 4.210e-01, 3.530e-01, 2.915e-01, 2.628e-01, 2.432e-01, 2.325e-01, 2.159e-01, 2.098e-01, 2.151e-01, 2.320e-01, 2.211e-01, 2.129e-01, 2.146e-01, 2.080e-01, 2.288e-01, 2.340e-01, 2.301e-01, 2.014e-01, 3.423e-01, 3.503e-01, 4.195e-01};
    Double_t PHENIX200GeVRAAErr_0020[23] = {3.888e-02, 4.457e-02, 4.546e-02, 3.989e-02, 3.365e-02, 3.060e-02, 2.954e-02, 2.823e-02, 2.613e-02, 2.527e-02, 2.590e-02, 2.793e-02, 2.662e-02, 2.562e-02, 2.583e-02, 2.524e-02, 2.807e-02, 2.921e-02, 3.100e-02, 3.276e-02, 6.465e-02, 8.398e-02, 1.867e-01};
    Double_t PHENIX200GeVRAAErrPt_0020[23] = {4.479e-03, 3.746e-03, 3.546e-03, 3.178e-03, 3.045e-03, 3.169e-03, 2.937e-03, 3.246e-03, 3.523e-03, 4.126e-03, 4.827e-03, 4.355e-03, 4.936e-03, 5.727e-03, 6.751e-03, 7.533e-03, 9.515e-03, 1.176e-02, 9.017e-03, 1.637e-02, 4.057e-02, 6.601e-02, 1.594e-01};
    Double_t PHENIX200GeVNormError_0020 = TMath::Sqrt(6.85*6.85+9.7*9.7);
    
    
    TGraphAsymmErrors* graphPHENIX200GeVRAA_0010 = new TGraphAsymmErrors(23);
    TGraphAsymmErrors* graphPHENIX200GeVRAA_0020 = new TGraphAsymmErrors(23);
    TGraphAsymmErrors* graphPHENIX200GeVRAA_1020 = new TGraphAsymmErrors(23);

    for(Int_t i=0; i<23; i++){
        Double_t quadErr1 = TMath::Sqrt(PHENIX200GeVRAAErr_0010[i]*PHENIX200GeVRAAErr_0010[i] + PHENIX200GeVRAAErrPt_0010[i] * PHENIX200GeVRAAErrPt_0010[i]);
        Double_t quadErr2 = TMath::Sqrt(PHENIX200GeVRAAErr_1020[i]*PHENIX200GeVRAAErr_1020[i] + PHENIX200GeVRAAErrPt_1020[i] * PHENIX200GeVRAAErrPt_1020[i]);
        Double_t quadErr3 = TMath::Sqrt(PHENIX200GeVRAAErr_0020[i]*PHENIX200GeVRAAErr_0020[i] + PHENIX200GeVRAAErrPt_0020[i] * PHENIX200GeVRAAErrPt_0020[i]);
        
        graphPHENIX200GeVRAA_0010->SetPoint(i, PHENIX200GeVpt_0010[i],  PHENIX200GeVRAA_0010[i]);
        graphPHENIX200GeVRAA_0010->SetPointError(i, 0.05, 0.05, quadErr1,quadErr1);
        graphPHENIX200GeVRAA_1020->SetPoint(i, PHENIX200GeVpt_0010[i],  PHENIX200GeVRAA_1020[i]);
        graphPHENIX200GeVRAA_1020->SetPointError(i, 0.05, 0.05, quadErr2,quadErr2);
        graphPHENIX200GeVRAA_0020->SetPoint(i, PHENIX200GeVpt_0010[i],  PHENIX200GeVRAA_0020[i]);
        graphPHENIX200GeVRAA_0020->SetPointError(i, 0.05, 0.05, quadErr3,quadErr3);
    }

    Double_t PHENIX200GeVpt_2030[22] =  { 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75, 11, 13, 15, 17};
    Double_t PHENIX200GeVRAA_2030[22] =  {4.524e-01, 5.142e-01, 5.025e-01, 4.535e-01, 3.990e-01, 3.802e-01, 3.604e-01, 3.516e-01, 3.282e-01, 3.223e-01, 3.351e-01, 3.262e-01, 3.003e-01, 3.102e-01, 3.324e-01, 3.566e-01, 3.445e-01, 3.502e-01, 3.198e-01, 3.356e-01, 4.554e-01, 5.281e-01};
    Double_t PHENIX200GeVRAAErr_2030[22] =  {4.457e-02, 5.178e-02, 5.395e-02, 5.082e-02, 4.596e-02, 4.413e-02, 4.363e-02, 4.256e-02, 3.966e-02, 3.887e-02, 4.031e-02, 3.921e-02, 3.610e-02, 3.738e-02, 4.029e-02, 4.350e-02, 4.235e-02, 4.378e-02, 4.327e-02, 5.544e-02, 9.300e-02, 1.241e-01};
    Double_t PHENIX200GeVRAAErrPt_2030[22] = {6.227e-03, 5.463e-03, 5.220e-03, 5.047e-03, 5.031e-03, 5.654e-03, 5.252e-03, 5.860e-03, 6.770e-03, 8.000e-03, 9.911e-03, 8.409e-03, 9.800e-03, 1.209e-02, 1.535e-02, 1.810e-02, 2.266e-02, 2.903e-02, 2.056e-02, 3.941e-02, 7.952e-02, 1.661e-01};
    Double_t PHENIX200GeVNormError_2030 = TMath::Sqrt(8.09*8.09+9.7*9.7);
    Double_t PHENIX200GeVRAA_3040[22] =  {4.960e-01, 5.663e-01, 5.436e-01, 5.098e-01, 4.554e-01, 4.506e-01, 4.142e-01, 4.240e-01, 4.013e-01, 3.960e-01, 3.895e-01, 3.968e-01, 4.254e-01, 4.089e-01, 4.132e-01, 3.964e-01, 3.817e-01, 5.656e-01, 4.607e-01, 4.302e-01, 4.003e-01, 4.426e-01};
    Double_t PHENIX200GeVRAAErr_3040[22] = {4.826e-02, 5.672e-02, 5.819e-02, 5.701e-02, 5.237e-02, 5.224e-02, 5.011e-02, 5.130e-02, 4.848e-02, 4.774e-02, 4.685e-02, 4.769e-02, 5.109e-02, 4.924e-02, 4.997e-02, 4.864e-02, 4.746e-02, 7.081e-02, 6.236e-02, 7.068e-02, 8.396e-02, 1.121e-01};
    Double_t PHENIX200GeVRAAErrPt_3040[22] = {6.366e-03, 5.802e-03, 5.487e-03, 5.459e-03, 5.580e-03, 6.384e-03, 5.832e-03, 7.060e-03, 8.289e-03, 1.044e-02, 1.208e-02, 1.129e-02, 1.492e-02, 1.844e-02, 2.123e-02, 2.361e-02, 2.830e-02, 4.268e-02, 2.891e-02, 5.306e-02, 9.278e-02, 1.761e-01};
    Double_t PHENIX200GeVNormError_3040 = TMath::Sqrt(8.41*8.41+9.7*9.7);

    TGraphAsymmErrors* graphPHENIX200GeVRAA_2030 = new TGraphAsymmErrors(22);
    TGraphAsymmErrors* graphPHENIX200GeVRAA_3040 = new TGraphAsymmErrors(22);
    TGraphAsymmErrors* graphPHENIX200GeVRAA_2040 = new TGraphAsymmErrors(22);
    
    for(Int_t i=0; i<22; i++){
        Double_t quadErr1 = TMath::Sqrt(PHENIX200GeVRAAErr_2030[i]*PHENIX200GeVRAAErr_2030[i] + PHENIX200GeVRAAErrPt_2030[i] * PHENIX200GeVRAAErrPt_2030[i]);
        Double_t quadErr2 = TMath::Sqrt(PHENIX200GeVRAAErr_3040[i]*PHENIX200GeVRAAErr_3040[i] + PHENIX200GeVRAAErrPt_3040[i] * PHENIX200GeVRAAErrPt_3040[i]);
        graphPHENIX200GeVRAA_2030->SetPoint(i, PHENIX200GeVpt_2030[i],  PHENIX200GeVRAA_2030[i]);
        graphPHENIX200GeVRAA_2030->SetPointError(i, 0.05, 0.05, quadErr1,quadErr1);
        graphPHENIX200GeVRAA_3040->SetPoint(i, PHENIX200GeVpt_2030[i],  PHENIX200GeVRAA_3040[i]);
        graphPHENIX200GeVRAA_3040->SetPointError(i, 0.05, 0.05, quadErr2, quadErr2);
        Double_t val = (PHENIX200GeVRAA_2030[i]+PHENIX200GeVRAA_3040[i])/2.;
        Double_t errorStat = TMath::Sqrt(PHENIX200GeVRAAErrPt_3040[i]*PHENIX200GeVRAAErrPt_3040[i] + PHENIX200GeVRAAErrPt_2030[i]*PHENIX200GeVRAAErrPt_2030[i]);
        Double_t errorSys = (PHENIX200GeVRAAErr_2030[i] + PHENIX200GeVRAAErr_3040[i])/2.;
        Double_t error = TMath::Sqrt(errorStat*errorStat + errorSys*errorSys);
        graphPHENIX200GeVRAA_2040->SetPoint(i, PHENIX200GeVpt_2030[i],  val);
        graphPHENIX200GeVRAA_2040->SetPointError(i, 0.05, 0.05, error,error);
    }
    
    Double_t PHENIX200GeVpt_4050[21] =  { 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75, 11, 13, 15};
    Double_t PHENIX200GeVRAA_4050[21] =  {5.491e-01, 6.287e-01, 6.235e-01, 6.007e-01, 5.724e-01, 5.417e-01, 5.463e-01, 5.217e-01, 5.247e-01, 5.462e-01, 5.258e-01, 5.114e-01, 4.716e-01, 4.980e-01, 5.567e-01, 5.173e-01, 5.330e-01, 5.254e-01, 4.746e-01, 5.588e-01, 6.307e-01};
    Double_t PHENIX200GeVRAAErr_4050[21] = {4.981e-02, 6.237e-02, 6.655e-02, 6.711e-02, 6.576e-02, 6.280e-02, 6.609e-02, 6.313e-02, 6.341e-02, 6.591e-02, 6.330e-02, 6.151e-02, 5.678e-02, 6.036e-02, 6.767e-02, 6.333e-02, 6.705e-02, 6.598e-02, 6.548e-02, 9.548e-02, 1.309e-01};
    Double_t PHENIX200GeVRAAErrPt_4050[21] = {6.093e-03, 6.184e-03, 6.128e-03, 6.325e-03, 6.698e-03, 7.730e-03, 7.688e-03, 9.205e-03, 1.146e-02, 1.506e-02, 1.764e-02, 1.662e-02, 2.060e-02, 2.412e-02, 3.118e-02, 3.621e-02, 4.256e-02, 5.433e-02, 3.761e-02, 7.278e-02, 1.629e-01};
    Double_t PHENIX200GeVNormError_4050 = TMath::Sqrt(9.79*9.79+9.7*9.7);
    Double_t PHENIX200GeVRAA_5060[21] =  {6.159e-01, 6.968e-01, 6.792e-01, 6.665e-01, 6.524e-01, 6.293e-01, 6.568e-01, 6.693e-01, 6.746e-01, 6.435e-01, 6.289e-01, 6.099e-01, 6.658e-01, 6.153e-01, 6.759e-01, 6.030e-01, 5.504e-01, 8.479e-01, 6.345e-01, 6.580e-01, 1.251};
    Double_t PHENIX200GeVRAAErr_5060[21] = {5.569e-02, 6.895e-02, 7.234e-02, 7.429e-02, 7.477e-02, 7.277e-02, 7.924e-02, 8.077e-02, 8.129e-02, 7.736e-02, 7.543e-02, 7.314e-02, 8.038e-02, 7.542e-02, 8.197e-02, 7.437e-02, 6.911e-02, 1.068e-01, 8.744e-02, 1.129e-01, 2.325e-01};
    Double_t PHENIX200GeVRAAErrPt_5060[21] = {6.694e-03, 6.859e-03, 6.502e-03, 6.850e-03, 7.724e-03, 9.208e-03, 9.970e-03, 1.298e-02, 1.682e-02, 1.964e-02, 2.529e-02, 2.379e-02, 2.921e-02, 3.814e-02, 4.656e-02, 5.159e-02, 5.704e-02, 9.234e-02, 5.745e-02, 1.121e-01, 5.687e-01};
    Double_t PHENIX200GeVNormError_5060 =  TMath::Sqrt(15.9*15.9+9.7*9.7);

    TGraphAsymmErrors* graphPHENIX200GeVRAA_4050 = new TGraphAsymmErrors(21);
    TGraphAsymmErrors* graphPHENIX200GeVRAA_5060 = new TGraphAsymmErrors(21);
    TGraphAsymmErrors* graphPHENIX200GeVRAA_4060 = new TGraphAsymmErrors(21);
    
    for(Int_t i=0; i<21; i++){
        Double_t quadErr1 = TMath::Sqrt(PHENIX200GeVRAAErr_4050[i]*PHENIX200GeVRAAErr_4050[i] + PHENIX200GeVRAAErrPt_4050[i] * PHENIX200GeVRAAErrPt_4050[i]);
        Double_t quadErr2 = TMath::Sqrt(PHENIX200GeVRAAErr_5060[i]*PHENIX200GeVRAAErr_5060[i] + PHENIX200GeVRAAErrPt_5060[i] * PHENIX200GeVRAAErrPt_5060[i]);
        graphPHENIX200GeVRAA_4050->SetPoint(i, PHENIX200GeVpt_4050[i],  PHENIX200GeVRAA_4050[i]);
        graphPHENIX200GeVRAA_4050->SetPointError(i, 0.05, 0.05, quadErr1, quadErr1);
        graphPHENIX200GeVRAA_5060->SetPoint(i, PHENIX200GeVpt_4050[i],  PHENIX200GeVRAA_5060[i]);
        graphPHENIX200GeVRAA_5060->SetPointError(i, 0.05, 0.05, quadErr2, quadErr2);
        Double_t val = (PHENIX200GeVRAA_4050[i]+PHENIX200GeVRAA_5060[i])/2.;
        Double_t errorStat = TMath::Sqrt(PHENIX200GeVRAAErrPt_5060[i]*PHENIX200GeVRAAErrPt_5060[i] + PHENIX200GeVRAAErrPt_4050[i]*PHENIX200GeVRAAErrPt_4050[i]);
        Double_t errorSys = (PHENIX200GeVRAAErr_4050[i] + PHENIX200GeVRAAErr_5060[i])/2.;
        Double_t error = TMath::Sqrt(errorStat*errorStat + errorSys*errorSys);
        graphPHENIX200GeVRAA_4060->SetPoint(i, PHENIX200GeVpt_4050[i],  val);
        graphPHENIX200GeVRAA_4060->SetPointError(i, 0.05, 0.05, error,error);
    }
    
    Double_t PHENIX200GeVpt_6070[20] =  {1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75, 11, 13};
    Double_t PHENIX200GeVRAA_6070[20] =  {6.802e-01, 7.472e-01, 7.341e-01, 7.351e-01, 7.309e-01, 7.178e-01, 7.337e-01, 7.286e-01, 7.362e-01, 7.611e-01, 7.287e-01, 7.747e-01, 7.629e-01, 6.946e-01, 7.066e-01, 7.793e-01, 6.991e-01, 8.433e-01, 6.871e-01, 8.026e-01};
    Double_t PHENIX200GeVRAAErr_6070[20] = {6.492e-02, 7.429e-02, 7.818e-02, 8.183e-02, 8.361e-02, 8.282e-02, 8.832e-02, 8.771e-02, 8.851e-02, 9.135e-02, 8.715e-02, 9.315e-02, 1.024e-01, 8.469e-02, 8.688e-02, 9.648e-02, 8.815e-02, 1.077e-01, 9.681e-02, 1.415e-01};
    Double_t PHENIX200GeVRAAErrPt_6070[20] = {7.526e-03, 7.176e-03, 6.863e-03, 7.567e-03, 8.958e-03, 1.125e-02, 1.332e-02, 1.773e-02, 2.207e-02, 3.030e-02, 3.868e-02, 3.576e-02, 4.908e-02, 5.205e-02, 6.345e-02, 8.107e-02, 9.110e-02, 1.273e-01, 8.190e-02, 1.725e-01};
    Double_t PHENIX200GeVNormError_6070 =  TMath::Sqrt(26.5*26.5+9.7*9.7);

    TGraphAsymmErrors* graphPHENIX200GeVRAA_6070 = new TGraphAsymmErrors(20);
    for(Int_t i=0; i<20; i++){
        Double_t quadErr1 = TMath::Sqrt(PHENIX200GeVRAAErr_6070[i]*PHENIX200GeVRAAErr_6070[i] + PHENIX200GeVRAAErrPt_6070[i] * PHENIX200GeVRAAErrPt_6070[i]);
        graphPHENIX200GeVRAA_6070->SetPoint(i, PHENIX200GeVpt_6070[i],  PHENIX200GeVRAA_6070[i]);
        graphPHENIX200GeVRAA_6070->SetPointError(i, 0.05, 0.05, quadErr1,quadErr1);
    }
    
    Double_t PHENIX200GeVpt_7080[19] =  { 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75, 11};
    Double_t PHENIX200GeVRAA_7080[19] =  {7.109e-01, 7.693e-01, 7.679e-01, 7.949e-01, 7.739e-01, 7.685e-01, 8.134e-01, 7.718e-01, 8.385e-01, 8.252e-01, 8.340e-01, 8.073e-01, 8.410e-01, 7.400e-01, 6.670e-01, 7.737e-01, 6.927e-01, 6.069e-01, 6.219e-01};
    Double_t PHENIX200GeVRAAErr_7080[19] = {6.546e-02, 7.607e-02, 8.163e-02, 8.840e-02, 8.851e-02, 8.867e-02, 9.795e-02, 9.298e-02, 1.008e-01, 9.903e-02, 1.000e-01, 9.649e-02, 1.004e-01, 9.264e-02, 8.267e-02, 9.727e-02, 8.793e-02, 7.846e-02, 8.971e-02};
    Double_t PHENIX200GeVRAAErrPt_7080[19] = {7.397e-03, 7.337e-03, 7.270e-03, 8.709e-03, 1.078e-02, 1.434e-02, 1.883e-02, 2.421e-02, 3.440e-02, 4.602e-02, 5.977e-02, 7.890e-02, 9.907e-02, 8.159e-02, 9.088e-02, 1.145e-01, 1.351e-01, 1.533e-01, 1.108e-01};
    Double_t PHENIX200GeVNormError_7080 =  TMath::Sqrt(33.3*33.3+9.7*9.7);

    TGraphAsymmErrors* graphPHENIX200GeVRAA_7080 = new TGraphAsymmErrors(19);
    TGraphAsymmErrors* graphPHENIX200GeVRAA_6080 = new TGraphAsymmErrors(19);
    for(Int_t i=0; i<19; i++){
        Double_t quadErr1 = TMath::Sqrt(PHENIX200GeVRAAErr_7080[i]*PHENIX200GeVRAAErr_7080[i] + PHENIX200GeVRAAErrPt_7080[i] * PHENIX200GeVRAAErrPt_7080[i]);
        graphPHENIX200GeVRAA_7080->SetPoint(i, PHENIX200GeVpt_7080[i],  PHENIX200GeVRAA_7080[i]);
        graphPHENIX200GeVRAA_7080->SetPointError(i, 0.05, 0.05, quadErr1, quadErr1);
        Double_t val = (PHENIX200GeVRAA_6070[i]+PHENIX200GeVRAA_7080[i])/2.;
        Double_t errorStat = TMath::Sqrt(PHENIX200GeVRAAErrPt_7080[i]*PHENIX200GeVRAAErrPt_7080[i] + PHENIX200GeVRAAErrPt_6070[i]*PHENIX200GeVRAAErrPt_6070[i]);
        Double_t errorSys = (PHENIX200GeVRAAErr_6070[i] + PHENIX200GeVRAAErr_7080[i])/2.;
        Double_t error = TMath::Sqrt(errorStat*errorStat + errorSys*errorSys);
        graphPHENIX200GeVRAA_6080->SetPoint(i, PHENIX200GeVpt_6070[i],  val);
        graphPHENIX200GeVRAA_6080->SetPointError(i, 0.05, 0.05, error,error);

    }
    
    
    //======================= eta RAA 0-10% ======================
    Double_t PHENIX200GeVEtapt_0010[9] =  { 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13, 15, 17};
    Double_t PHENIX200GeVEtaRAA_0010[9] =  {1.96e-01, 1.64e-01, 2.14e-01, 2.50e-01, 2.37e-01, 2.36e-01, 3.64e-01, 1.51e-01, 1.70e-01};
    Double_t PHENIX200GeVEtaRAAStatErr_0010[9] = {1.71e-02, 1.88e-02, 2.38e-02, 3.06e-02, 4.13e-02, 4.14e-02, 6.67e-02, 9.41e-02, 1.44e-01};
    Double_t PHENIX200GeVEtaRAASysErrB_0010[9] = {2.64e-02, 2.21e-02, 2.89e-02, 3.37e-02, 3.20e-02, 3.18e-02, 4.92e-02, 2.03e-02, 2.30e-02};
    Double_t PHENIX200GeVEtaRAASysErrC_0010[9] = {2.39e-02, 2.00e-02, 2.61e-02, 3.04e-02, 2.89e-02, 2.87e-02, 4.44e-02, 1.84e-02, 2.07e-02};
    
    //======================= eta RAA 20-40% ======================
    Double_t PHENIX200GeVEtapt_2040[10] =  { 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13, 15, 17, 19};
    Double_t PHENIX200GeVEtaRAA_2040[10] =  {3.92e-01, 3.85e-01, 4.11e-01, 3.62e-01, 4.77e-01, 3.33e-01, 7.00e-01, 1.57e-01, 6.33e-01, 2.11};
    Double_t PHENIX200GeVEtaRAAStatErr_2040[10] = {2.33e-02, 2.56e-02, 3.21e-02, 4.05e-02, 6.25e-02, 6.41e-02, 1.12e-01, 1.37e-01, 2.62e-01, 1.80};
    Double_t PHENIX200GeVEtaRAASysErrB_2040[10] = {5.29e-02, 5.20e-02, 5.55e-02, 4.88e-02, 6.44e-02, 4.50e-02, 9.45e-02, 2.11e-02, 8.55e-02, 2.85e-01};
    Double_t PHENIX200GeVEtaRAASysErrC_2040[10] = {4.78e-02, 4.70e-02, 5.02e-02, 4.41e-02, 5.82e-02, 4.07e-02, 8.54e-02, 1.91e-02, 7.73e-02, 2.57e-01};

    //======================= eta RAA 20-60% ======================
    Double_t PHENIX200GeVEtapt_2060[10] =  { 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13, 15, 17, 19};
    Double_t PHENIX200GeVEtaRAA_2060[10] =  {4.52e-01, 4.47e-01, 4.41e-01, 4.27e-01, 5.54e-01, 3.79e-01, 6.95e-01, 3.60e-01, 7.05e-01, 1.57};
    Double_t PHENIX200GeVEtaRAAStatErr_2060[10] = {2.38e-02, 2.56e-02, 3.08e-02, 3.92e-02, 6.12e-02, 6.38e-02, 1.04e-01, 1.35e-01, 2.47e-01, 8.16e-01};
    Double_t PHENIX200GeVEtaRAASysErrB_2060[10] = {6.10e-02, 6.03e-02, 5.96e-02, 5.76e-02, 7.47e-02, 5.12e-02, 9.39e-02, 4.86e-02, 9.52e-02, 2.11e-01};
    Double_t PHENIX200GeVEtaRAASysErrC_2060[10] = {5.51e-02, 5.45e-02, 5.38e-02, 5.21e-02, 6.75e-02, 4.63e-02, 8.48e-02, 4.39e-02, 8.60e-02, 1.91e-01};
    

    TGraphAsymmErrors* graphPHENIX200GeVEtaRAA_0010 = new TGraphAsymmErrors(9);
    TGraphAsymmErrors* graphPHENIX200GeVEtaRAA_2040 = new TGraphAsymmErrors(10);
    TGraphAsymmErrors* graphPHENIX200GeVEtaRAA_2060 = new TGraphAsymmErrors(10);

    for(Int_t i=0; i<9; i++){
        Double_t quadErrEta_0010 = TMath::Sqrt(PHENIX200GeVEtaRAAStatErr_0010[i]*PHENIX200GeVEtaRAAStatErr_0010[i] + PHENIX200GeVEtaRAASysErrB_0010[i]*PHENIX200GeVEtaRAASysErrB_0010[i] + PHENIX200GeVEtaRAASysErrC_0010[i] * PHENIX200GeVEtaRAASysErrC_0010[i]);
        
        graphPHENIX200GeVEtaRAA_0010->SetPoint(i, PHENIX200GeVEtapt_0010[i],  PHENIX200GeVEtaRAA_0010[i]);
        graphPHENIX200GeVEtaRAA_0010->SetPointError(i, 0.05, 0.05, quadErrEta_0010,quadErrEta_0010);
    }	

    for(Int_t i=0; i<10; i++){
        Double_t quadErrEta_2040 = TMath::Sqrt(PHENIX200GeVEtaRAAStatErr_2040[i]*PHENIX200GeVEtaRAAStatErr_2040[i] + PHENIX200GeVEtaRAASysErrB_2040[i]*PHENIX200GeVEtaRAASysErrB_2040[i] + PHENIX200GeVEtaRAASysErrC_2040[i] * PHENIX200GeVEtaRAASysErrC_2040[i]);
        
        graphPHENIX200GeVEtaRAA_2040->SetPoint(i, PHENIX200GeVEtapt_2040[i],  PHENIX200GeVEtaRAA_2040[i]);
        graphPHENIX200GeVEtaRAA_2040->SetPointError(i, 0.05, 0.05, quadErrEta_2040,quadErrEta_2040);
        
        Double_t quadErrEta_2060 = TMath::Sqrt(PHENIX200GeVEtaRAAStatErr_2060[i]*PHENIX200GeVEtaRAAStatErr_2060[i] + PHENIX200GeVEtaRAASysErrB_2060[i]*PHENIX200GeVEtaRAASysErrB_2060[i] + PHENIX200GeVEtaRAASysErrC_2060[i] * PHENIX200GeVEtaRAASysErrC_2060[i]);
        
        graphPHENIX200GeVEtaRAA_2060->SetPoint(i, PHENIX200GeVEtapt_2060[i],  PHENIX200GeVEtaRAA_2060[i]);
        graphPHENIX200GeVEtaRAA_2060->SetPointError(i, 0.05, 0.05, quadErrEta_2060,quadErrEta_2060);

    }	

    //============================================================	
    
    
    

    Double_t PHENIX39GeVpt_0010[18] = {0.90, 1.10, 1.30, 1.50, 1.70,
                                                1.90, 2.10, 2.30, 2.50, 2.70, 
                                                2.90, 3.20, 3.80, 4.20, 4.80,
                                                5.50, 6.50, 7.50};
    Double_t PHENIX39GeVRAA_0010[18] =  {3.69e-01, 4.47e-01, 4.19e-01, 4.85e-01, 4.70e-01, 4.57e-01, 4.76e-01, 5.09e-01, 5.39e-01, 5.37e-01, 5.63e-01, 5.39e-01, 5.25e-01, 4.96e-01, 4.91e-01, 4.68e-01, 4.28e-01, 3.77e-01};
    Double_t PHENIX39GeVRAAErrPt_0010[18] = {7.54e-05, 1.13e-04, 1.33e-04, 1.98e-04, 2.51e-04, 3.24e-04, 4.49e-04, 6.42e-04, 9.01e-04, 1.21e-03, 1.66e-03, 2.70e-03, 5.01e-03, 8.88e-03, 1.54e-02, 3.18e-02, 8.92e-02, 2.10e-01};
    Double_t PHENIX39GeVRAAErr_0010[18] = {8.85e-02, 8.33e-02, 9.69e-02, 9.40e-02, 9.18e-02, 1.08e-01, 1.13e-01, 1.32e-01, 1.07e-01, 1.05e-01, 9.98e-02, 1.06e-01, 1.01e-01, 1.05e-01, 9.88e-02, 9.35e-02, 1.38e-01};

    TGraphAsymmErrors* graphPHENIX39GeVRAA_0010 = new TGraphAsymmErrors(18);
    for(Int_t i=0; i<18; i++){
        Double_t quadErr = TMath::Sqrt(PHENIX39GeVRAAErr_0010[i]*PHENIX39GeVRAAErr_0010[i] + PHENIX39GeVRAAErrPt_0010[i] * PHENIX39GeVRAAErrPt_0010[i]);
        graphPHENIX39GeVRAA_0010->SetPoint(i, PHENIX39GeVpt_0010[i],  PHENIX39GeVRAA_0010[i]);
        graphPHENIX39GeVRAA_0010->SetPointError(i, 0.0, 0.0, quadErr,quadErr);
    }

    Double_t PHENIX39GeVpt_0020[18] = {0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5};
    Double_t PHENIX39GeVRAA_0020[18] =  {3.40e-01, 4.32e-01, 4.11e-01, 4.81e-01, 4.66e-01, 4.60e-01, 4.88e-01, 5.26e-01, 5.55e-01, 5.66e-01, 5.93e-01, 5.79e-01, 5.80e-01, 5.59e-01, 5.54e-01, 5.26e-01, 4.83e-01, 3.93e-01};
    Double_t PHENIX39GeVRAAErrPt_0020[18] = {4.08e-02, 1.73e-02, 7.80e-03, 3.56e-03, 1.64e-03, 7.98e-04, 5.63e-04, 5.38e-04, 7.11e-04, 9.59e-04, 1.31e-03, 2.16e-03, 4.08e-03, 7.32e-03, 1.28e-02, 2.75e-02, 7.40e-02, 1.67e-01};
    Double_t PHENIX39GeVRAAErr_0020[18] = {6.70e-02, 8.51e-02, 8.09e-02, 9.51e-02, 9.20e-02, 9.09e-02, 1.10e-01, 1.14e-01, 1.34e-01, 1.10e-01, 1.08e-01, 1.04e-01, 1.14e-01, 1.08e-01, 1.09e-01, 1.04e-01, 9.73e-02, 1.08e-01};

    TGraphAsymmErrors* graphPHENIX39GeVRAA_0020 = new TGraphAsymmErrors(18);
    for(Int_t i=0; i<18; i++){
        Double_t quadErr = TMath::Sqrt(PHENIX39GeVRAAErr_0020[i]*PHENIX39GeVRAAErr_0020[i] + PHENIX39GeVRAAErrPt_0020[i] * PHENIX39GeVRAAErrPt_0020[i]);
        graphPHENIX39GeVRAA_0020->SetPoint(i, PHENIX39GeVpt_0020[i],  PHENIX39GeVRAA_0020[i]);
        graphPHENIX39GeVRAA_0020->SetPointError(i, 0.0, 0.0, quadErr,quadErr);
    }

    Double_t PHENIX39GeVpt_2040[19] = {0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 8.5};
    Double_t PHENIX39GeVRAA_2040[19] =  {4.99e-01, 6.06e-01, 5.38e-01, 6.16e-01, 5.92e-01, 5.88e-01, 6.25e-01, 6.76e-01, 7.20e-01, 7.62e-01, 7.77e-01, 7.94e-01, 7.97e-01, 7.52e-01, 7.70e-01, 8.05e-01, 7.61e-01, 9.44e-01, 1.08};
    Double_t PHENIX39GeVRAAErrPt_2040[19] = {1.10e-04, 1.58e-04, 1.78e-04, 2.61e-04, 3.27e-04, 4.26e-04, 5.96e-04, 8.58e-04, 1.21e-03, 1.68e-03, 2.28e-03, 3.85e-03, 7.30e-03, 1.30e-02, 2.32e-02, 5.58e-02, 1.43e-01, 3.96e-01, 7.03e-01};
    Double_t PHENIX39GeVRAAErr_2040[19] = {9.86e-02, 1.20e-01, 1.07e-01, 1.23e-01, 1.18e-01, 1.18e-01, 1.42e-01, 1.50e-01, 1.76e-01, 1.51e-01, 1.45e-01, 1.47e-01, 1.61e-01, 1.53e-01, 1.64e-01, 1.70e-01, 1.66e-01, 3.47e-01, 4.24e-01};

    TGraphAsymmErrors* graphPHENIX39GeVRAA_2040 = new TGraphAsymmErrors(19);
    for(Int_t i=0; i<19; i++){
        Double_t quadErr = TMath::Sqrt(PHENIX39GeVRAAErr_2040[i]*PHENIX39GeVRAAErr_2040[i] + PHENIX39GeVRAAErrPt_2040[i] * PHENIX39GeVRAAErrPt_2040[i]);
        graphPHENIX39GeVRAA_2040->SetPoint(i, PHENIX39GeVpt_2040[i],  PHENIX39GeVRAA_2040[i]);
        graphPHENIX39GeVRAA_2040->SetPointError(i, 0.0, 0.0, quadErr,quadErr);
    }

    Double_t PHENIX39GeVpt_4060[18] = {0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5};
    Double_t PHENIX39GeVRAA_4060[18] =  {5.55e-01, 7.02e-01, 6.26e-01, 7.13e-01, 6.95e-01, 6.93e-01, 7.37e-01, 8.21e-01, 8.64e-01, 9.16e-01, 9.82e-01, 9.81e-01, 1.06, 1.10, 1.20, 1.09, 9.42e-01, 8.26e-01};
    Double_t PHENIX39GeVRAAErrPt_4060[18] = {1.91e-04, 2.85e-04, 3.23e-04, 4.76e-04, 6.04e-04, 7.90e-04, 1.11e-03, 1.62e-03, 2.28e-03, 3.18e-03, 4.42e-03, 7.39e-03, 1.46e-02, 2.73e-02, 5.04e-02, 1.13e-01, 2.78e-01, 6.49e-01};
    Double_t PHENIX39GeVRAAErr_4060[18] = {1.10e-01, 1.39e-01, 1.24e-01, 1.42e-01, 1.39e-01, 1.39e-01, 1.68e-01, 1.82e-01, 2.11e-01, 1.82e-01, 1.84e-01, 1.82e-01, 2.14e-01, 2.23e-01, 2.56e-01, 2.31e-01, 2.06e-01, 3.03e-01};

    TGraphAsymmErrors* graphPHENIX39GeVRAA_4060 = new TGraphAsymmErrors(18);
    for(Int_t i=0; i<18; i++){
        Double_t quadErr = TMath::Sqrt(PHENIX39GeVRAAErr_4060[i]*PHENIX39GeVRAAErr_4060[i] + PHENIX39GeVRAAErrPt_4060[i] * PHENIX39GeVRAAErrPt_4060[i]);
        graphPHENIX39GeVRAA_4060->SetPoint(i, PHENIX39GeVpt_4060[i],  PHENIX39GeVRAA_4060[i]);
        graphPHENIX39GeVRAA_4060->SetPointError(i, 0.0, 0.0, quadErr,quadErr);
    }

    Double_t PHENIX39GeVpt_6086[18] = {0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5};
    Double_t PHENIX39GeVRAA_6086[18] =  {7.62e-01, 9.47e-01, 8.29e-01, 9.31e-01, 8.88e-01, 8.74e-01, 9.53e-01, 1.03, 1.10, 1.16, 1.21, 1.27, 1.22, 1.45, 1.06, 1.34, 8.41e-01, 1.84};
    Double_t PHENIX39GeVRAAErrPt_6086[18] = {4.68e-04, 6.93e-04, 7.81e-04, 1.14e-03, 1.44e-03, 1.87e-03, 2.65e-03, 3.82e-03, 5.41e-03, 7.53e-03, 1.03e-02, 1.77e-02, 3.30e-02, 6.61e-02, 9.98e-02, 2.65e-01, 5.55e-01, 2.05};
    Double_t PHENIX39GeVRAAErr_6086[18] = {1.51e-01, 1.88e-01, 1.65e-01, 1.86e-01, 1.78e-01, 1.75e-01, 2.17e-01, 2.28e-01, 2.68e-01, 2.30e-01, 2.27e-01, 2.36e-01, 2.46e-01, 2.94e-01, 2.26e-01, 2.83e-01, 1.84e-01, 6.77e-01};

    TGraphAsymmErrors* graphPHENIX39GeVRAA_6086 = new TGraphAsymmErrors(18);
    for(Int_t i=0; i<18; i++){
        Double_t quadErr = TMath::Sqrt(PHENIX39GeVRAAErr_6086[i]*PHENIX39GeVRAAErr_6086[i] + PHENIX39GeVRAAErrPt_6086[i] * PHENIX39GeVRAAErrPt_6086[i]);
        graphPHENIX39GeVRAA_6086->SetPoint(i, PHENIX39GeVpt_6086[i],  PHENIX39GeVRAA_6086[i]);
        graphPHENIX39GeVRAA_6086->SetPointError(i, 0.0, 0.0, quadErr,quadErr);
    }

    Double_t PHENIX62GeVpt_0010[23] = {0.90, 1.10, 1.30, 1.50, 1.70, 1.90, 2.10, 2.30, 2.50, 2.70, 2.90, 3.20, 3.80, 4.20, 4.80, 5.20, 5.80, 6.20, 6.8, 7.2, 7.8, 8.5, 9.5};
    Double_t PHENIX62GeVRAA_0010[23] =  {3.57e-01, 3.73e-01, 4.31e-01, 4.75e-01, 4.88e-01, 5.11e-01, 5.03e-01, 4.94e-01, 4.72e-01, 4.56e-01, 4.35e-01, 4.04e-01, 3.82e-01, 3.59e-01, 3.16e-01, 3.18e-01, 2.92e-01, 2.31e-01, 2.62e-01, 2.33e-01, 2.34e-01, 1.89e-01, 2.09e-01};
    Double_t PHENIX62GeVRAAErrPt_0010[23] = {4.59e-05, 5.92e-05, 8.33e-05, 1.15e-04, 1.52e-04, 2.06e-04, 2.64e-04, 3.43e-04, 4.29e-04, 5.40e-04, 6.70e-04, 6.56e-03, 1.10e-02, 1.74e-02, 2.53e-02, 3.78e-02, 4.74e-02, 5.83e-02, 2.77e-02, 3.25e-02, 4.56e-02, 4.56e-02, 1.58e-01};
    Double_t PHENIX62GeVRAAErr_0010[23] = {3.96e-02, 3.72e-02, 4.26e-02, 4.78e-02, 5.00e-02, 5.35e-02, 5.82e-02, 6.26e-02, 6.30e-02, 6.38e-02, 6.39e-02, 6.20e-02, 6.11e-02, 6.08e-02, 5.77e-02, 6.48e-02, 6.81e-02, 6.93e-02, 4.99e-02, 4.68e-02, 4.93e-02, 4.25e-02, 5.11e-02};
    TGraphAsymmErrors* graphPHENIX62GeVRAA_0010 = new TGraphAsymmErrors(23);
    for(Int_t i=0; i<23; i++){
        Double_t quadErr = TMath::Sqrt(PHENIX62GeVRAAErr_0010[i]*PHENIX62GeVRAAErr_0010[i] + PHENIX62GeVRAAErrPt_0010[i] * PHENIX62GeVRAAErrPt_0010[i]);
        graphPHENIX62GeVRAA_0010->SetPoint(i, PHENIX62GeVpt_0010[i],  PHENIX62GeVRAA_0010[i]);
        graphPHENIX62GeVRAA_0010->SetPointError(i, 0.0, 0.0, quadErr,quadErr);
    }

    Double_t PHENIX62GeVpt_0020[23] = {0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.5, 9.5};
    Double_t PHENIX62GeVRAA_0020[23] =  {3.51e-01, 3.75e-01, 4.28e-01, 4.69e-01, 4.87e-01, 5.12e-01, 5.08e-01, 5.05e-01, 4.81e-01, 4.76e-01, 4.50e-01, 4.31e-01, 4.08e-01, 3.95e-01, 3.72e-01, 3.51e-01, 3.30e-01, 2.89e-01, 3.19e-01, 2.42e-01, 2.47e-01, 2.04e-01, 2.03e-01};
    Double_t PHENIX62GeVRAAErrPt_0020[23] = {3.52e-05, 4.56e-05, 6.35e-05, 8.72e-05, 1.16e-04, 1.58e-04, 2.03e-04, 2.65e-04, 3.31e-04, 4.23e-04, 5.22e-04, 6.96e-03, 1.17e-02, 1.90e-02, 2.96e-02, 4.14e-02, 5.31e-02, 7.22e-02, 2.99e-02, 2.97e-02, 4.04e-02, 4.26e-02, 3.74e-02};
    Double_t PHENIX62GeVRAAErr_0020[23] = {3.89e-02, 3.75e-02, 4.24e-02, 4.72e-02, 4.99e-02, 5.37e-02, 5.89e-02, 6.40e-02, 6.42e-02, 6.67e-02, 6.61e-02, 6.61e-02, 6.52e-02, 6.68e-02, 6.79e-02, 7.15e-02, 7.71e-02, 8.68e-02, 6.09e-02, 4.80e-02, 5.15e-02, 4.60e-02, 4.72e-02};
    TGraphAsymmErrors* graphPHENIX62GeVRAA_0020 = new TGraphAsymmErrors(23);
    for(Int_t i=0; i<23; i++){
        Double_t quadErr = TMath::Sqrt(PHENIX62GeVRAAErr_0020[i]*PHENIX62GeVRAAErr_0020[i] + PHENIX62GeVRAAErrPt_0020[i] * PHENIX62GeVRAAErrPt_0020[i]);
        graphPHENIX62GeVRAA_0020->SetPoint(i, PHENIX62GeVpt_0020[i],  PHENIX62GeVRAA_0020[i]);
        graphPHENIX62GeVRAA_0020->SetPointError(i, 0.0, 0.0, quadErr,quadErr);
    }

    Double_t PHENIX62GeVpt_2040[23] = {0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.5, 9.5};
    Double_t PHENIX62GeVRAA_2040[23] =  {4.99e-01, 5.17e-01, 5.57e-01, 5.94e-01, 6.18e-01, 6.56e-01, 6.47e-01, 6.58e-01, 6.42e-01, 6.39e-01, 6.22e-01, 6.06e-01, 5.91e-01, 5.83e-01, 5.71e-01, 5.64e-01, 5.36e-01, 4.43e-01, 4.19e-01, 3.48e-01, 3.47e-01, 2.74e-01, 2.32e-01};
    Double_t PHENIX62GeVRAAErrPt_2040[23] = {6.97e-05, 8.51e-05, 1.13e-04, 1.52e-04, 2.01e-04, 2.74e-04, 3.50e-04, 4.64e-04, 5.86e-04, 7.52e-04, 9.42e-04, 9.83e-03, 1.70e-02, 2.82e-02, 4.56e-02, 6.68e-02, 8.66e-02, 1.11e-01, 4.32e-02, 4.78e-02, 6.48e-02, 7.11e-02, 1.46e-01};
    Double_t PHENIX62GeVRAAErr_2040[23] = {5.53e-02, 5.17e-02, 5.51e-02, 5.98e-02, 6.33e-02, 6.87e-02, 7.49e-02, 8.34e-02, 8.57e-02, 8.96e-02, 9.13e-02, 9.30e-02, 9.46e-02, 9.87e-02, 1.04e-01, 1.15e-01, 1.25e-01, 1.33e-01, 8.01e-02, 7.00e-02, 7.31e-02, 6.17e-02, 5.70e-02};
    TGraphAsymmErrors* graphPHENIX62GeVRAA_2040 = new TGraphAsymmErrors(23);
    for(Int_t i=0; i<23; i++){
        Double_t quadErr = TMath::Sqrt(PHENIX62GeVRAAErr_2040[i]*PHENIX62GeVRAAErr_2040[i] + PHENIX62GeVRAAErrPt_2040[i] * PHENIX62GeVRAAErrPt_2040[i]);
        graphPHENIX62GeVRAA_2040->SetPoint(i, PHENIX62GeVpt_2040[i],  PHENIX62GeVRAA_2040[i]);
        graphPHENIX62GeVRAA_2040->SetPointError(i, 0.0, 0.0, quadErr,quadErr);
    }

    Double_t PHENIX62GeVpt_4060[23] = {0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.5, 9.5};
    Double_t PHENIX62GeVRAA_4060[23] =  {5.63e-01, 6.13e-01, 6.56e-01, 7.00e-01, 7.34e-01, 7.73e-01, 7.74e-01, 7.90e-01, 7.84e-01, 7.70e-01, 7.72e-01, 7.59e-01, 7.65e-01, 7.70e-01, 7.24e-01, 8.13e-01, 8.00e-01, 6.85e-01, 6.12e-01, 4.93e-01, 4.93e-01, 5.23e-01, 5.24e-01};
    Double_t PHENIX62GeVRAAErrPt_4060[23] = {1.25e-04, 1.60e-04, 2.14e-04, 2.88e-04, 3.84e-04, 5.21e-04, 6.70e-04, 8.87e-04, 1.13e-03, 1.44e-03, 1.82e-03, 1.25e-02, 2.23e-02, 3.77e-02, 5.87e-02, 9.76e-02, 1.31e-01, 1.75e-01, 7.54e-02, 8.50e-02, 1.22e-01, 1.47e-01, 1.72e-01};
    Double_t PHENIX62GeVRAAErr_4060[23] = {6.24e-02, 6.13e-02, 6.49e-02, 7.05e-02, 7.53e-02, 8.10e-02, 8.96e-02, 1.00e-01, 1.05e-01, 1.08e-01, 1.13e-01, 1.16e-01, 1.22e-01, 1.30e-01, 1.32e-01, 1.66e-01, 1.87e-01, 2.05e-01, 1.17e-01, 9.91e-02, 1.04e-01, 1.18e-01, 1.29e-01};
    TGraphAsymmErrors* graphPHENIX62GeVRAA_4060 = new TGraphAsymmErrors(23);
    for(Int_t i=0; i<23; i++){
        Double_t quadErr = TMath::Sqrt(PHENIX62GeVRAAErr_4060[i]*PHENIX62GeVRAAErr_4060[i] + PHENIX62GeVRAAErrPt_4060[i] * PHENIX62GeVRAAErrPt_4060[i]);
        graphPHENIX62GeVRAA_4060->SetPoint(i, PHENIX62GeVpt_4060[i],  PHENIX62GeVRAA_4060[i]);
        graphPHENIX62GeVRAA_4060->SetPointError(i, 0.0, 0.0, quadErr,quadErr);
    }

    Double_t PHENIX62GeVpt_6086[22] = {0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.5};
    Double_t PHENIX62GeVRAA_6086[22] =  {7.18e-01, 7.56e-01, 7.85e-01, 8.28e-01, 8.52e-01, 9.14e-01, 9.11e-01, 9.26e-01, 9.25e-01, 9.14e-01, 9.13e-01, 9.40e-01, 9.54e-01, 9.35e-01, 1.02, 9.44e-01, 1.05, 1.12, 8.27e-01, 5.40e-01, 9.01e-01, 6.44e-01};
    Double_t PHENIX62GeVRAAErrPt_6086[22] = {3.00e-04, 3.75e-04, 4.91e-04, 6.58e-04, 8.69e-04, 1.19e-03, 1.53e-03, 2.02e-03, 2.58e-03, 3.30e-03, 4.19e-03, 1.65e-02, 2.95e-02, 4.88e-02, 8.72e-02, 1.21e-01, 1.86e-01, 3.01e-01, 1.56e-01, 2.01e-01, 2.85e-01, 3.55e-01};
    Double_t PHENIX62GeVRAAErr_6086[22] = {7.96e-02, 7.55e-02, 7.76e-02, 8.34e-02, 8.73e-02, 9.58e-02, 1.06e-01, 1.17e-01, 1.24e-01, 1.28e-01, 1.34e-01, 1.44e-01, 1.53e-01, 1.58e-01, 1.86e-01, 1.92e-01, 2.45e-01, 3.36e-01, 1.58e-01, 1.09e-01, 1.90e-01, 1.45e-01};
    TGraphAsymmErrors* graphPHENIX62GeVRAA_6086 = new TGraphAsymmErrors(22);
    for(Int_t i=0; i<22; i++){
        Double_t quadErr = TMath::Sqrt(PHENIX62GeVRAAErr_6086[i]*PHENIX62GeVRAAErr_6086[i] + PHENIX62GeVRAAErrPt_6086[i] * PHENIX62GeVRAAErrPt_6086[i]);
        graphPHENIX62GeVRAA_6086->SetPoint(i, PHENIX62GeVpt_6086[i],  PHENIX62GeVRAA_6086[i]);
        graphPHENIX62GeVRAA_6086->SetPointError(i, 0.0, 0.0, quadErr,quadErr);
    }

    TFile* fileWA98Pi0 = 	new TFile("ExternalInputPbPb/OtherExperiments/wa98_pi0_raa_pC_and_pPb_ref.root");
        TGraphErrors* graphWA09_0013 = (TGraphErrors*)fileWA98Pi0->Get("g_raa_pC_ref_0_13");
        
    // data from:
    // 200 GeV: http://arxiv.org/abs/1208.2254, link to data: http://www.phenix.bnl.gov/phenix/WWW/info/data/ppg133_data.html
    // 39, 62.4 GeV: http://arxiv.org/abs/1204.1526, link to data: http://www.phenix.bnl.gov/phenix/WWW/info/data/ppg138_data.html

    TGraphErrors* g200 = new TGraphErrors(0);
    TGraphErrors* g62 = new TGraphErrors(0);
    TGraphErrors* g39 = new TGraphErrors(0);

    g200->SetPoint(g200->GetN(), 325.8, 0.191494434545626); g200->SetPointError(g200->GetN()-1, 0, 0.0165748162217721);
    // g200->SetPoint( 0-5%, 0.17840589873838); g200->SetPointError(0, 0.0156417727408418);
    g200->SetPoint(g200->GetN(), 236.1, 0.263384119812022); g200->SetPointError(g200->GetN()-1, 0, 0.0227886686616168);
    g200->SetPoint(g200->GetN(), 167.6, 0.336436899849347); g200->SetPointError(g200->GetN()-1, 0, 0.0291548689271431);
    g200->SetPoint(g200->GetN(), 115.5, 0.426075881697806); g200->SetPointError(g200->GetN()-1, 0, 0.0370418089090795);
    g200->SetPoint(g200->GetN(), 76.2, 0.537046084601494); g200->SetPointError(g200->GetN()-1, 0, 0.0469661110954586);
    g200->SetPoint(g200->GetN(), 47.1, 0.664455346234181); g200->SetPointError(g200->GetN()-1, 0, 0.0587351485351021);
    g200->SetPoint(g200->GetN(), 26.7, 0.784509965872809); g200->SetPointError(g200->GetN()-1, 0, 0.0710844439819864);
    g200->SetPoint(g200->GetN(), 13.7, 0.892646228451866); g200->SetPointError(g200->GetN()-1, 0, 0.0854157417473091);
    g200->SetPoint(g200->GetN(), 5.6, 0.925898260451632); g200->SetPointError(g200->GetN()-1, 0, 0.0982863802862068);

    g62->SetPoint(g62->GetN(), 319.6, 0.246133064086216); g62->SetPointError(g62->GetN()-1, 0, 0.0455733674205351);
    g62->SetPoint(g62->GetN(), 229.7, 0.302722085145605); g62->SetPointError(g62->GetN()-1, 0, 0.0594455619587974);
    g62->SetPoint(g62->GetN(), 138.7, 0.378744052011265); g62->SetPointError(g62->GetN()-1, 0, 0.0696310088710886);
    g62->SetPoint(g62->GetN(), 59.7, 0.542710158865597); g62->SetPointError(g62->GetN()-1, 0, 0.110228758494566);
    g62->SetPoint(g62->GetN(), 12.4, 0.632545718550685); g62->SetPointError(g62->GetN()-1, 0, 0.199780167231788);

    g39->SetPoint(g39->GetN(), 319.6, 0.411954122992487); g39->SetPointError(g39->GetN()-1, 0, 0.168439541823257);
    g39->SetPoint(g39->GetN(), 229.7, 0.526984615384615); g39->SetPointError(g39->GetN()-1, 0, 0.22920344506741);
    g39->SetPoint(g39->GetN(), 138.7, 0.795080681242186); g39->SetPointError(g39->GetN()-1, 0, 0.308308947471857);
    g39->SetPoint(g39->GetN(), 59.7, 0.905331401541957); g39->SetPointError(g39->GetN()-1, 0, 0.493815618562743);

    
    
    TFile fileTheoryGraphs("ExternalInputPbPb/OtherExperiments/DataCompilationFromOtherEnergiesPbPbWithEta.root","RECREATE");
        graphWA09_0013->Write("graphWA98RAA_0013");
        graphPHENIX200GeVRAA_0010->Write("graphPHENIX200GeVRAA_0010");
        graphPHENIX200GeVRAA_0020->Write("graphPHENIX200GeVRAA_0020");
        graphPHENIX200GeVRAA_1020->Write("graphPHENIX200GeVRAA_1020");
        graphPHENIX200GeVRAA_2030->Write("graphPHENIX200GeVRAA_2030");
        graphPHENIX200GeVRAA_3040->Write("graphPHENIX200GeVRAA_3040");
        graphPHENIX200GeVRAA_2040->Write("graphPHENIX200GeVRAA_2040");
        graphPHENIX200GeVRAA_4050->Write("graphPHENIX200GeVRAA_4050");
        graphPHENIX200GeVRAA_5060->Write("graphPHENIX200GeVRAA_5060");
        graphPHENIX200GeVRAA_4060->Write("graphPHENIX200GeVRAA_4060");
        graphPHENIX200GeVRAA_6070->Write("graphPHENIX200GeVRAA_6070");
        graphPHENIX200GeVRAA_7080->Write("graphPHENIX200GeVRAA_7080");
        graphPHENIX200GeVRAA_6080->Write("graphPHENIX200GeVRAA_6080");
        graphPHENIX39GeVRAA_0010->Write("graphPHENIX39GeVRAA_0010");
        graphPHENIX39GeVRAA_0020->Write("graphPHENIX39GeVRAA_0020");
        graphPHENIX39GeVRAA_2040->Write("graphPHENIX39GeVRAA_2040");
        graphPHENIX39GeVRAA_4060->Write("graphPHENIX39GeVRAA_4060");
        graphPHENIX39GeVRAA_6086->Write("graphPHENIX39GeVRAA_6086");
        graphPHENIX62GeVRAA_0010->Write("graphPHENIX62GeVRAA_0010");
        graphPHENIX62GeVRAA_0020->Write("graphPHENIX62GeVRAA_0020");
        graphPHENIX62GeVRAA_2040->Write("graphPHENIX62GeVRAA_2040");
        graphPHENIX62GeVRAA_4060->Write("graphPHENIX62GeVRAA_4060");
        graphPHENIX62GeVRAA_6086->Write("graphPHENIX62GeVRAA_6086");
        graphPHENIX200GeVEtaRAA_0010->Write("graphPHENIX200GeVEtaRAA_0010");
        graphPHENIX200GeVEtaRAA_2040->Write("graphPHENIX200GeVEtaRAA_2040");
        graphPHENIX200GeVEtaRAA_2060->Write("graphPHENIX200GeVEtaRAA_2060");
    g200->Write("graphPHENIX200GeVRAAvsNPartAt7GeVc");
    g62->Write("graphPHENIX62GeVRAAvsNPartAt7GeVc");
    g39->Write("graphPHENIX39GeVRAAvsNPartAt7GeVc");
    fileTheoryGraphs.Close();
}
