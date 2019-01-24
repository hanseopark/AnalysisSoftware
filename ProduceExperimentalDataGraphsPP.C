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
#include "CommonHeaders/ConversionFunctions.h"

extern TRandom*    gRandom;
extern TBenchmark*    gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*      gMinuit;

void RemoveUpperLimits(TGraphErrors* tg){
    int N = tg->GetN();
    double *x = tg->GetX();
    double *y = tg->GetY();
    double *ey = tg->GetEY();
    for(int i = 0; i < N; i++){
        if(ey[i]/y[i] > 0.8){
            y[i]=0;
            ey[i]=0;
        }
    }

}


void ProduceExperimentalDataGraphsPP(){

    // LMEE pp
    //==========================================================================================
    //                                         ALICE 7TeV
    //==========================================================================================
    // arXiv:1805.04391
    TString aliceLMEE7TeVHEPDataFile        = "ExternalInput/OtherParticles/LMEE_7TeV_Rgamma_HEPData-ins1672792-v1-Table_11.csv" ;
    TGraphAsymmErrors* RGammaStat_7TeV_LMEE_ALICE   = ParseHEPData(aliceLMEE7TeVHEPDataFile, 8, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
    TGraphAsymmErrors* RGammaSys_7TeV_LMEE_ALICE    = ParseHEPData(aliceLMEE7TeVHEPDataFile, 8, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE, kFALSE);
    TGraphAsymmErrors* RGammaTot_7TeV_LMEE_ALICE    = AddErrorsOfGraphsQuadratically(RGammaStat_7TeV_LMEE_ALICE, RGammaSys_7TeV_LMEE_ALICE);

    // Photons pp
    //==========================================================================================
    //                                         ALICE prelim isolated
    //==========================================================================================
//     Measurement:
//     Isolation Criteria: Charged+neutral Iso <2 GeV/c
//     measured in |eta|<0.27
    Float_t xbins[10]={10., 12., 14., 16., 18., 20., 25., 30., 40., 60.};
    Double_t pt[9];
    Double_t pterr[9];
    for (Int_t i = 0; i < 9; i++){
        pt[i]       = (xbins[i]+xbins[i+1])/2;
        pterr[i]    = (xbins[i]-xbins[i+1])/2;
    }
    Double_t sigma[9]={31.6,19.18,10.7,5.8,4.68,1.98,0.66,0.378,0.059};
    Double_t estatsigma[9]={2.84,1.81,1.07,0.8,0.67,0.26,0.21,0.09,0.024};
    Double_t esystlsigma[9]={4.74,2.58,1.4,0.8,0.8,0.39,0.155,0.08,0.011};
//
    TGraphAsymmErrors *tStat_7TeV_IsoPh_ALICE       = new TGraphAsymmErrors(9, pt, sigma, pterr, pterr, estatsigma, estatsigma);
    TGraphAsymmErrors *tSys_7TeV_IsoPh_ALICE        = new TGraphAsymmErrors(9, pt, sigma, pterr, pterr, esystlsigma, esystlsigma);
    TGraphAsymmErrors* tTot_7TeV_IsoPh_ALICE        = AddErrorsOfGraphsQuadratically(tStat_7TeV_IsoPh_ALICE, tSys_7TeV_IsoPh_ALICE);
    TGraphAsymmErrors *tTotxT_7TeV_IsoPh_ALICE      = xTScalePhoton(tTot_7TeV_IsoPh_ALICE, 7000);

    //==========================================================================================
    //                                         CMS
    //==========================================================================================
    // Phys.Rev. D84 (2011) 052011
    // Phys.Rev.Lett. 106 (2011) 082001
    const int N_CMS = 11;
    double pT_CMS[N_CMS] = {
        22, 24.5, 28.0, 32.5, 37.5, 42.5, 47.5, 55., 72.5, 102.5, 210.
    };
    double y_CMS[N_CMS] = {
        2.17, 1.39, 0.774, 0.402, 0.209, 0.1244, 0.0740, 0.0403, 0.01236, 0.00243, 0.000188
    };
    double ey_CMS[N_CMS] = {
        0.03, 0.02, 0.010, 0.006, 0.004, 0.0028, 0.0021, 0.0010, 0.00035, 0.00012, 0.000013
    };
    double sUPy_CMS[N_CMS] = {
        13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 14.0, 14.0, 13.0
    };
    double sLOWy_CMS[N_CMS] = {
        16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 9.0, 9.0
    };
    double E_CMS = 7000.;

    double sUy_CMS[N_CMS] = {0};
    double sLy_CMS[N_CMS] = {0};
    for(int i = 0; i < N_CMS; i++){
        sUy_CMS[i] = y_CMS[i]*sUPy_CMS[i];
        sLy_CMS[i] = y_CMS[i]*sLOWy_CMS[i];
    }
    TGraphErrors *tE_DirPh_CMS          = new TGraphErrors(N_CMS, pT_CMS, y_CMS, 0, ey_CMS);
    TGraphAsymmErrors *tS_DirPh_CMS     = new TGraphAsymmErrors(N_CMS, pT_CMS, y_CMS, 0, 0, sLy_CMS, sUy_CMS);

    TGraphErrors *tExT_DirPh_CMS        = xTScalePhoton(tE_DirPh_CMS, E_CMS);

    //Phys.Rev.Lett. 106 (2011) 082001, 2011
    TString cms7TeVHEPDataFile        = "ExternalInput/OtherExperiments/CMS_7TeV_HEPData-ins922830-v1-Table_1.csv" ;
    TGraphAsymmErrors* tStat_7TeV_DirPh_CMS   = ParseHEPData(cms7TeVHEPDataFile, 8, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
    TGraphAsymmErrors* tSys_7TeV_DirPh_CMS    = ParseHEPData(cms7TeVHEPDataFile, 8, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE, kFALSE);
    TGraphAsymmErrors* tTot_7TeV_DirPh_CMS    = AddErrorsOfGraphsQuadratically(tStat_7TeV_DirPh_CMS, tSys_7TeV_DirPh_CMS);
    TGraphAsymmErrors *tTotxT_7TeV_DirPh_CMS  = xTScalePhoton(tTot_7TeV_DirPh_CMS, E_CMS);

    //Phys.Lett. B710 (2012) 256-277, 2012
    TString cms2760GeVHEPDataFile        = "ExternalInput/OtherExperiments/CMS_2760GeV_HEPData-ins1084729-v1-Table_3.csv" ;
    TGraphAsymmErrors* tStat_2760GeV_DirPh_CMS   = ParseHEPData(cms2760GeVHEPDataFile, 8, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
    TGraphAsymmErrors* tSys_2760GeV_DirPh_CMS    = ParseHEPData(cms2760GeVHEPDataFile, 8, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE, kFALSE);
    tStat_2760GeV_DirPh_CMS                      = ScaleGraph(tStat_2760GeV_DirPh_CMS, 1e-3);
    tSys_2760GeV_DirPh_CMS                       = ScaleGraph(tSys_2760GeV_DirPh_CMS, 1e-3);
    TGraphAsymmErrors* tTot_2760GeV_DirPh_CMS    = AddErrorsOfGraphsQuadratically(tStat_2760GeV_DirPh_CMS, tSys_2760GeV_DirPh_CMS);

    tTot_2760GeV_DirPh_CMS->Print();
    TGraphAsymmErrors *tTotxT_2760GeV_DirPh_CMS  = xTScalePhoton(tTot_2760GeV_DirPh_CMS, 2760.);
    tTotxT_2760GeV_DirPh_CMS->Print();
    //==========================================================================================
    //                                         ATLAS
    //==========================================================================================
    // Phys.Lett. B706 (2011) 150-167, 2011
    const int N_ATLAS = 8;
    double pT_ATLAS[N_ATLAS] = { 50.0, 62.5, 77.5, 92.5, 112.5, 137.5, 175.0, 300.0 };
    double y_ATLAS[N_ATLAS] = { 83.3, 32.7, 12.3, 5.3, 2.2, 0.8, 0.26, 0.028 };
    double sy_ATLAS[N_ATLAS] = { 7.210409142344143, 2.6851443164195103, 0.938083151964686, 0.4123105625617661, 0.20615528128088306, 0.06324555320336758, 0.022360679774997897, 0.0031304951684997056 };
    double ey_ATLAS[N_ATLAS] = { 0.5, 0.3, 0.2, 0.1, 0.05, 0.03, 0.01, 0.002 };
    double E_ATLAS = 7000.;

    TGraphErrors *tE_DirPh_ATLAS    = new TGraphErrors(N_ATLAS, pT_ATLAS, y_ATLAS, 0, ey_ATLAS);
    TGraphErrors *tS_DirPh_ATLAS    = new TGraphErrors(N_ATLAS, pT_ATLAS, y_ATLAS, 0, sy_ATLAS);
    tE_DirPh_ATLAS                  = ScaleGraph(tE_DirPh_ATLAS, 0.001);
    tS_DirPh_ATLAS                  = ScaleGraph(tS_DirPh_ATLAS, 0.001);
//     styleme(tE_DirPh_ATLAS, 22);
    TGraphErrors *tExT_DirPh_ATLAS  = xTScalePhoton(tE_DirPh_ATLAS, E_ATLAS);

    // Phys.Rev. D83 (2011) 052005, 2011
    TString atlas7TeVHEPDataFile        = "ExternalInput/OtherExperiments/ATLAS_7TeV_HEPData-ins882463-v1-Table_1.csv" ;
    TGraphAsymmErrors* tStat_7TeV_DirPh_ATLAS   = ParseHEPData(atlas7TeVHEPDataFile, 10, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
    TGraphAsymmErrors* tSys_7TeV_DirPh_ATLAS    = ParseHEPData(atlas7TeVHEPDataFile, 10, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE, kFALSE);
    TGraphAsymmErrors* tTot_7TeV_DirPh_ATLAS    = AddErrorsOfGraphsQuadratically(tStat_7TeV_DirPh_ATLAS, tSys_7TeV_DirPh_ATLAS);
    TGraphAsymmErrors *tTotxT_7TeV_DirPh_ATLAS  = xTScalePhoton(tTot_7TeV_DirPh_ATLAS, E_ATLAS);

    //JHEP 1608 (2016) 005, 2016
    TString atlas8TeVHEPDataFile        = "ExternalInput/OtherExperiments/ATLAS_8TeV_HEPData-ins1457605-v1-Table_1.csv" ;
    TGraphAsymmErrors* tStat_8TeV_DirPh_ATLAS   = ParseHEPData(atlas8TeVHEPDataFile, 10, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
    TGraphAsymmErrors* tSys_8TeV_DirPh_ATLAS    = ParseHEPData(atlas8TeVHEPDataFile, 10, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE, kFALSE);
    tStat_8TeV_DirPh_ATLAS                      = ScaleGraph(tStat_8TeV_DirPh_ATLAS, 0.001);
    tSys_8TeV_DirPh_ATLAS                       = ScaleGraph(tSys_8TeV_DirPh_ATLAS, 0.001);
    TGraphAsymmErrors* tTot_8TeV_DirPh_ATLAS    = AddErrorsOfGraphsQuadratically(tStat_8TeV_DirPh_ATLAS, tSys_8TeV_DirPh_ATLAS);
    TGraphAsymmErrors *tTotxT_8TeV_DirPh_ATLAS  = xTScalePhoton(tTot_8TeV_DirPh_ATLAS, 8000.);

    // Phys.Lett. B770 (2017) 473-493
//     TFile* fileatlas13TeV                       = new TFile("ExternalInput/OtherExperiments/ATLAS_13TeV_compilation_Hendrik.root");
//     TGraphAsymmErrors* tTot_13TeV_DirPh_ATLAS   = (TGraphAsymmErrors*)fileatlas13TeV->Get("atlas2017_pt/atlas2017_pt_data_rap1");
//     tTot_13TeV_DirPh_ATLAS                      = ScaleGraph(tTot_13TeV_DirPh_ATLAS, 1e-3);
//     TGraphAsymmErrors *tTotxT_13TeV_DirPh_ATLAS = xTScalePhoton(tTot_13TeV_DirPh_ATLAS, 13000.);

    TFile* fileatlas13TeVStat                   = new TFile("ExternalInput/OtherExperiments/xsections_incphoton_ATLAS_13TeV.root");
    TH1D* histStat_13TeV_DirPh_ATLAS            = (TH1D*)fileatlas13TeVStat->Get("hxsect_eta_0");
    TGraphAsymmErrors* tStat_13TeV_DirPh_ATLAS  = new TGraphAsymmErrors(histStat_13TeV_DirPh_ATLAS);
    RemoveZerosAtBeginningAndEndFromGraph(tStat_13TeV_DirPh_ATLAS);
    TFile* fileatlas13TeVSys                    = new TFile("ExternalInput/OtherExperiments/total_syst_incphoton_ATLAS_13TeV.root");
    TH1D* histRelSysDown_13TeV_DirPh_ATLAS      = (TH1D*)fileatlas13TeVSys->Get("histo_totalup_eta_0");
    TGraphAsymmErrors* tRelSysDown13TeVATLAS    = new TGraphAsymmErrors(histRelSysDown_13TeV_DirPh_ATLAS);
    RemoveZerosAtBeginningAndEndFromGraph(tRelSysDown13TeVATLAS);
    TH1D* histRelSysUp_13TeV_DirPh_ATLAS        = (TH1D*)fileatlas13TeVSys->Get("histo_totaldown_eta_0");
    TGraphAsymmErrors* tRelSysUp13TeVATLAS    = new TGraphAsymmErrors(histRelSysUp_13TeV_DirPh_ATLAS);
    RemoveZerosAtBeginningAndEndFromGraph(tRelSysUp13TeVATLAS);
    TGraphAsymmErrors* tSys_13TeV_DirPh_ATLAS   = (TGraphAsymmErrors*)tStat_13TeV_DirPh_ATLAS->Clone("tSys_13TeV_DirPh_ATLAS");
    for (Int_t i = 0; i< tSys_13TeV_DirPh_ATLAS->GetN(); i++){
        tSys_13TeV_DirPh_ATLAS->SetPointError(i, tSys_13TeV_DirPh_ATLAS->GetEXlow()[i], tSys_13TeV_DirPh_ATLAS->GetEXhigh()[i], tRelSysDown13TeVATLAS->GetY()[i]*tSys_13TeV_DirPh_ATLAS->GetY()[i], tRelSysUp13TeVATLAS->GetY()[i]*tSys_13TeV_DirPh_ATLAS->GetY()[i] );
    }

    TGraphAsymmErrors* tTot_13TeV_DirPh_ATLAS   = AddErrorsOfGraphsQuadratically(tStat_13TeV_DirPh_ATLAS, tSys_13TeV_DirPh_ATLAS);
    tTot_13TeV_DirPh_ATLAS                      = ScaleGraph(tTot_13TeV_DirPh_ATLAS, 1e-3);
    TGraphAsymmErrors *tTotxT_13TeV_DirPh_ATLAS = xTScalePhoton(tTot_13TeV_DirPh_ATLAS, 13000.);

    //==========================================================================================
    //                                         D0
    //==========================================================================================
    // Phys.Rev.Lett. 77 (1996) 5011-5015, 1996
    const int N_D0 = 23;
    double pT_D0[N_D0] = { 10.5, 13.5, 16.5, 19.5, 22.5, 25.5, 28.5, 31.5, 37.4,
    40.5, 43.5, 46.5, 49.5, 52.5, 55.5, 58.5, 61.5, 65.7, 72.0,
    78.0, 85.1, 94.4, 108.4 };
    double y_D0[N_D0] = { 14200.0, 4010.0, 1290.0, 661.0, 315.0, 157.0, 105.0, 63.9, 28.1,
    19.7, 14.2, 11.1, 7.76, 5.83, 4.72, 3.3, 2.98, 2.17, 1.39,
    0.88, 0.67, 0.28, 0.11 };
    double sy_D0[N_D0] = { 5262.052546297879, 1013.61383179197, 261.00766272276377, 113.81124724736128, 40.85339643163099, 21.95449840010015, 15.0, 9.859513172565876, 2.4738633753705965,
    1.772004514666935, 1.2649110640673518, 1.077032961426901, 0.73824115301167, 0.5730619512757762, 0.47169905660283024, 0.35608987629529715, 0.33421549934136807, 0.23021728866442678, 0.1562049935181331,
    0.1131370849898476, 0.08485281374238571, 0.044721359549995794, 0.022360679774997897 };
    double ey_D0[N_D0] = { 291.0, 153.0, 90.0, 67.0, 15.0, 11.0, 9.0, 7.5, 0.6,
    0.5, 0.4, 0.4, 0.31, 0.28, 0.25, 0.22, 0.21, 0.13, 0.1,
    0.08, 0.06, 0.04, 0.02 };
    double E_D0 = 1800.;

    TGraphErrors *tE_DirPh_D0 = new TGraphErrors(N_D0, pT_D0, y_D0, 0, ey_D0);
    TGraphErrors *tS_DirPh_D0 = new TGraphErrors(N_D0, pT_D0, y_D0, 0, sy_D0);
    tE_DirPh_D0               = ScaleGraph(tE_DirPh_D0, 0.001*2*TMath::Pi());
    tS_DirPh_D0               = ScaleGraph(tS_DirPh_D0, 0.001*2*TMath::Pi());
    //     styleme(tE_DirPh_D0, 25);
    TGraphErrors *tExT_DirPh_D0  = xTScalePhoton(tE_DirPh_D0, E_D0);

    //==========================================================================================
    //                                         CDF
    //==========================================================================================
    // Phys.Rev.Lett. 73 (1994) 2662-2666, 1994
    const int N_CDF = 16;
    double pT_CDF[N_CDF] = { 12.3, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 29.0, 31.0,
    33.9, 37.9, 41.9, 48.9, 62.4, 80.8, 114.7 };
    double y_CDF[N_CDF] = { 4460.0, 1300.0, 805.0, 458.0, 308.0, 226.0, 163.0, 106.0, 76.7,
    53.7, 30.9, 20.5, 7.61, 3.09, 0.911, 0.163 };
    double sy_CDF[N_CDF] = { 825.8456272209716, 160.5615146914104, 91.44397191723465, 48.38388161361178, 33.24154027718932, 25.079872407968907, 17.88854381999832, 12.529964086141668, 8.823831367382313,
    5.768882040742383, 3.6878177829171546, 2.6172504656604803, 1.0748023074035522, 0.4455333881989093, 0.18319934497699494, 0.04477722635447622 };
    double ey_CDF[N_CDF] = { 415.0, 38.0, 21.0, 15.0, 12.0, 10.0, 8.0, 6.0, 5.5,
    3.2, 2.4, 1.9, 0.76, 0.32, 0.159, 0.041 };
    double E_CDF = 1800;
    TGraphErrors *tE_DirPh_CDF  = new TGraphErrors(N_CDF, pT_CDF, y_CDF, 0, ey_CDF);
    TGraphErrors *tS_DirPh_CDF  = new TGraphErrors(N_CDF, pT_CDF, y_CDF, 0, sy_CDF);
    tE_DirPh_CDF                = ScaleGraph(tE_DirPh_CDF, 0.001*2*TMath::Pi());
    tS_DirPh_CDF                = ScaleGraph(tS_DirPh_CDF, 0.001*2*TMath::Pi());
    //     styleme(tE_DirPh_CDF, 24);
    TGraphErrors *tExT_DirPh_CDF  = xTScalePhoton(tE_DirPh_CDF, E_CDF);


    //==========================================================================================
    //                                         UA1_0
    //==========================================================================================
    // Phys.Lett. B209 (1988) 385-396
    const int N_UA1_0 = 16;
    double pT_UA1_0[N_UA1_0] = { 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 29.0, 31.5, 34.5,
    37.5, 40.5, 46.0, 55.0, 65.0, 75.0, 90.0 };
    double y_UA1_0[N_UA1_0] = { 6.42, 3.3, 1.54, 0.74, 0.5, 0.381, 0.246, 0.123, 0.056,
    0.051, 0.03, 0.0111, 0.0039, 0.0037, 0.0013, 2.0E-4 };
    double sy_UA1_0[N_UA1_0] = { 1.2567020331009258, 0.6158733636065129, 0.2973213749463701, 0.1140175425099138, 0.07810249675906655, 0.0604400529450463, 0.043046486500061765, 0.023259406699226014, 0.014560219778561038,
    0.013341664064126334, 0.01019803902718557, 0.003522782990761708, 0.0018027756377319946, 0.0016, 9.0E-4, 2.0E-4 };
    double ey_UA1_0[N_UA1_0] = { 0.57, 0.33, 0.2, 0.07, 0.05, 0.047, 0.037, 0.021, 0.014,
    0.013, 0.01, 0.0035, 0.0018, 0.0016, 9.0E-4, 2.0E-4 };
    double E_UA1_0 = 630;
    TGraphErrors *tE_DirPh_UA1_0 = new TGraphErrors(N_UA1_0, pT_UA1_0, y_UA1_0, 0, ey_UA1_0);
    TGraphErrors *tS_DirPh_UA1_0 = new TGraphErrors(N_UA1_0, pT_UA1_0, y_UA1_0, 0, sy_UA1_0);
    RemoveUpperLimits(tE_DirPh_UA1_0);
    RemoveUpperLimits(tS_DirPh_UA1_0);
    //     styleme(tE_DirPh_UA1_0, 27);
    TGraphErrors *tExT_DirPh_UA1_0  = xTScalePhoton(tE_DirPh_UA1_0, E_UA1_0);
    //==========================================================================================
    //                                         UA1_1
    //==========================================================================================
    // CERN-EP-89-138
    const int N_UA1_1 = 6;
    double pT_UA1_1[N_UA1_1] = { 17.0, 19.0, 21.0, 25.0, 34.5, 46.0 };
    double y_UA1_1[N_UA1_1] = { 3.91, 1.74, 1.12, 0.38, 0.049, 0.0084 };
    double sy_UA1_1[N_UA1_1] = { 0.5748912940721924, 0.29410882339705485, 0.21470910553583888, 0.0670820393249937, 0.013152946437965905, 0.006003332407921454 };
    double ey_UA1_1[N_UA1_1] = { 0.37, 0.24, 0.19, 0.06, 0.013, 0.006 };
    double E_UA1_1 = 546;
    TGraphErrors *tE_DirPh_UA1_1 = new TGraphErrors(N_UA1_1, pT_UA1_1, y_UA1_1, 0, ey_UA1_1);
    TGraphErrors *tS_DirPh_UA1_1 = new TGraphErrors(N_UA1_1, pT_UA1_1, y_UA1_1, 0, sy_UA1_1);
    RemoveUpperLimits(tE_DirPh_UA1_1);
    RemoveUpperLimits(tS_DirPh_UA1_1);
    //     styleme(tE_DirPh_UA1_1, 28);
    TGraphErrors *tExT_DirPh_UA1_1   = xTScalePhoton(tE_DirPh_UA1_1, E_UA1_1);
    //==========================================================================================
    //                                         UA2
    //==========================================================================================
    // Phys.Lett. B263 (1991) 544-550
    const int N_UA2 = 13;
    double pT_UA2[N_UA2] = { 15.9, 17.9, 19.9, 21.9, 23.9, 25.9, 28.7, 33.5, 38.6,
    46.3, 54.1, 64.5, 82.3 };
    double y_UA2[N_UA2] = { 7.46, 3.97, 1.79, 0.992, 0.615, 0.366, 0.151, 0.0657, 0.0179,
    0.00694, 0.00231, 4.84E-4, 1.51E-4 };
    double sy_UA2[N_UA2] = { 1.4684004903295285, 0.7164091010030512, 0.337249166047894, 0.17425466995176916, 0.0937469466169432, 0.05783121994217311, 0.024233035303073365, 0.010589357865328755, 0.00403624825797423,
    0.0018672439583514522, 9.989479465918132E-4, 2.780319406111463E-4, 1.0094681768139103E-4 };
    double ey_UA2[N_UA2] = { 0.41, 0.251, 0.156, 0.0713, 0.05, 0.0362, 0.016, 0.00728, 0.00367,
    0.00171, 9.36E-4, 2.72E-4, 9.99E-5 };
    double E_UA2 = 630;
    TGraphErrors *tE_DirPh_UA2 = new TGraphErrors(N_UA2, pT_UA2, y_UA2, 0, ey_UA2);
    TGraphErrors *tS_DirPh_UA2 = new TGraphErrors(N_UA2, pT_UA2, y_UA2, 0, sy_UA2);
    RemoveUpperLimits(tE_DirPh_UA2);
    RemoveUpperLimits(tS_DirPh_UA2);
    //     styleme(tE_DirPh_UA2, 26);
    TGraphErrors *tExT_DirPh_UA2 = xTScalePhoton(tE_DirPh_UA2, E_UA2);

    //==========================================================================================
    //                                         PHENIX
    //==========================================================================================
    // Phys. Rev. D86, 072008 (2012)
    const int N_PHENIX = 24;
    double pT_PHENIX[N_PHENIX] = {1.18449, 1.69126, 2.19773, 2.70364, 3.34880, 4.37494, 5.25,5.75,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75,11.00,13.00,15.00,17.00,19.00,21.00,23.00,25.00};
    double y_PHENIX[N_PHENIX] = {0.001489490e+09, 0.000541715e+09, 0.000111656e+09,  4.84535e+04, 9.89977e+03, 3.03996e+03, 1.14e+03,6.13e+02,3.48e+02,2.31e+02,1.36e+02,9.29e+01,6.70e+01,4.83e+01,3.21e+01,2.04e+01,9.81e+00,2.97e+00,1.06e+00,3.38e-01,1.73e-01,8.82e-02,4.22e-02,2.87e-02};
    double ey_PHENIX[N_PHENIX] = {0.000758668e+09, 0.000163956e+09, 3.04992e+04, 1.18865e+04, 3.26577e+03, 1.01686e+03, 3.04e+01,1.92e+01,1.27e+01,8.50e+00,6.12e+00,4.41e+00,3.22e+00,2.45e+00,1.89e+00,1.46e+00,4.23e-01,1.89e-01,9.85e-02,5.51e-02,3.37e-02,2.06e-02,1.52e-02,1.41e-02};
    double sy_PHENIX[N_PHENIX] = {0.002163150e+09, 0.000304584e+09, 5.79046e+04, 1.65357e+04, 3.86545e+03, 6.28280e+02, 4.78e+02, 2.21e+02, 1.01e+02, 6.24e+01, 3.13e+01, 1.95e+01, 1.34e+01, 9.18e+00, 6.10e+00, 3.68e+00, 1.67e+00, 4.75e-01, 1.69e-01, 5.42e-02, 2.77e-02, 1.50e-02, 7.18e-03,4.30e-03};
    double E_PHENIX = 200;

    TGraphErrors *tE_DirPh_PHENIX = new TGraphErrors(N_PHENIX, pT_PHENIX, y_PHENIX, 0, ey_PHENIX);
    TGraphErrors *tS_DirPh_PHENIX = new TGraphErrors(N_PHENIX, pT_PHENIX, y_PHENIX, 0, sy_PHENIX);
//     styleme(tE_DirPh_PHENIX, 21);
    RemoveUpperLimits(tE_DirPh_PHENIX);
    RemoveUpperLimits(tS_DirPh_PHENIX);
    TGraphErrors *tExT_DirPh_PHENIX = xTScalePhoton(tE_DirPh_PHENIX, E_PHENIX);
    //==========================================================================================
    //                                         R110
    //==========================================================================================
    const int N_R110 = 7;
    double pT_R110[N_R110] = { 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 9.0 };
    double y_R110[N_R110] = { 2.25E-34, 1.41E-34, 6.56E-35, 3.94E-35, 1.96E-35, 7.21E-36, 1.37E-36 };
    double sy_R110[N_R110] = { 3.3E-35, 4.2E-35, 2.42E-35, 6.3E-36, 6.1E-36, 2.72E-36, 8.6E-37 };
    double ey_R110[N_R110] = { 3.3E-35, 4.2E-35, 2.42E-35, 6.3E-36, 6.1E-36, 2.72E-36, 8.6E-37 };
    double E_R110 = 63;
    //==========================================================================================
    //                                         R806
    //==========================================================================================

    //==========================================================================================
    //                                         R807
    //==========================================================================================
    const int N_R807 = 11;
    double pT_R807[N_R807] = { 4.75, 5.25, 5.73, 6.23, 6.74, 7.23, 7.72, 8.22, 8.74,
    9.44, 10.36 };
    double y_R807[N_R807] = { 314.0, 122.0, 59.6, 32.0, 18.7, 9.07, 6.2, 3.81, 2.54,
    1.26, 0.605 };
    double sy_R807[N_R807] = { 79.0, 27.0, 12.2, 6.3, 3.5, 1.82, 1.26, 0.86, 0.63,
    0.32, 0.192 };
    double ey_R807[N_R807] = { 79.0, 27.0, 12.2, 6.3, 3.5, 1.82, 1.26, 0.86, 0.63,
    0.32, 0.192 };
    double E_R807 = 63;

    TGraphErrors *tE_DirPh_R807 = new TGraphErrors(N_R807, pT_R807, y_R807, 0, ey_R807);
    TGraphErrors *tS_DirPh_R807 = new TGraphErrors(N_R807, pT_R807, y_R807, 0, sy_R807);
//     styleme(tE_DirPh_R807, 25);
    RemoveUpperLimits(tE_DirPh_R807);
    RemoveUpperLimits(tS_DirPh_R807);
    TGraphErrors *tExT_DirPh_R807 = xTScalePhoton(tE_DirPh_R807, E_R807);
    //==========================================================================================
    //                                         R108
    //==========================================================================================
    const int N_R108 = 8;
    double pT_R108[N_R108] = { 5.4, 6.4, 7.4, 8.4, 9.4, 10.4, 11.4, 12.4 };
    double y_R108[N_R108] = { 1.3E-34, 1.5E-35, 4.3E-36, 1.9E-36, 5.3E-37, 4.5E-37, 2.4E-37, 1.1E-37 };
    double ey_R108[N_R108] = { 8.0E-35, 1.0E-35, 1.1E-36, 5.0E-37, 1.6E-37, 1.0E-37, 6.0E-38, 4.0E-38 };
    double sy_R108[N_R108] = { 8.0E-35, 1.0E-35, 1.1E-36, 5.0E-37, 1.6E-37, 1.0E-37, 6.0E-38, 4.0E-38 };
    double E_R108 = 62.4;
    TGraphErrors *tE_DirPh_R108 = new TGraphErrors(N_R108, pT_R108, y_R108, 0, ey_R108);
    TGraphErrors *tS_DirPh_R108 = new TGraphErrors(N_R108, pT_R108, y_R108, 0, sy_R108);
    tE_DirPh_R108               = ScaleGraph(tE_DirPh_R108, 1e36);
    tS_DirPh_R108               = ScaleGraph(tS_DirPh_R108, 1e36);
    RemoveUpperLimits(tE_DirPh_R108);
    RemoveUpperLimits(tS_DirPh_R108);
    //     styleme(tE_DirPh_R108, 24);
    TGraphErrors *tExT_DirPh_R108 = xTScalePhoton(tE_DirPh_R108, E_R108);
    //==========================================================================================
    //                                         E706
    //==========================================================================================
    // Phys.Rev. D72 (2005) 032003
    // Phys.Rev. D70 (2004) 092009
    // Phys.Rev. D48 (1993) 5-28
    const int N_E706 = 9;
    double pT_E706[N_E706] = { 3.75, 4.25, 4.75, 5.25, 5.75, 6.5, 7.5, 9.0, 11.0 };
    double y_E706[N_E706] = { 1880.0, 600.0, 225.5, 82.2, 35.0, 9.8, 1.73, 0.339, 0.017 };
    double sy_E706[N_E706] = { 484.1487374764082, 109.2520022699813, 35.27548723972498, 11.84736257569591, 5.011985634456668, 1.4098226838861687, 0.3440930106817051, 0.08013114251026252, 0.01414213562373095 };
    double ey_E706[N_E706] = { 300.0, 44.0, 9.4, 4.4, 2.4, 0.74, 0.28, 0.07, 0.014 };
    double E_E706 = 38.8;
    TGraphErrors *tE_DirPh_E706 = new TGraphErrors(N_E706, pT_E706, y_E706, 0, ey_E706);
    TGraphErrors *tS_DirPh_E706 = new TGraphErrors(N_E706, pT_E706, y_E706, 0, sy_E706);
    RemoveUpperLimits(tE_DirPh_E706);
    RemoveUpperLimits(tS_DirPh_E706);
    //     styleme(tE_DirPh_E706, 32);
    TGraphErrors *tExT_DirPh_E706 = xTScalePhoton(tE_DirPh_E706, E_E706);
    const int N_E706_0 = 8;
    double pT_E706_0[N_E706_0] = { 3.75, 4.25, 4.75, 5.25, 5.75, 6.5, 7.5, 9.0 };
    double y_E706_0[N_E706_0] = { 1160.0, 375.0, 107.6, 33.6, 11.4, 3.43, 0.42, 0.01 };
    double sy_E706_0[N_E706_0] = { 290.68883707497264, 68.87670143089025, 16.62077013859466, 4.919349550499538, 1.7804493814764855, 0.5597320787662612, 0.12083045973594572, 0.015033296378372907 };
    double ey_E706_0[N_E706_0] = { 190.0, 30.0, 4.5, 2.2, 1.1, 0.37, 0.11, 0.015 };
    double E_E706_0 = 31.6;
    TGraphErrors *tE_DirPh_E706_0 = new TGraphErrors(N_E706_0, pT_E706_0, y_E706_0, 0, ey_E706_0);
    TGraphErrors *tS_DirPh_E706_0 = new TGraphErrors(N_E706_0, pT_E706_0, y_E706_0, 0, sy_E706_0);
    RemoveUpperLimits(tE_DirPh_E706_0);
    RemoveUpperLimits(tS_DirPh_E706_0);
    //     styleme(tE_DirPh_E706_0, 31);
    TGraphErrors *tExT_DirPh_E706_0 = xTScalePhoton(tE_DirPh_E706_0, E_E706_0);
    //==========================================================================================
    //                                         E704
    //==========================================================================================
    // Phys.Lett. B345 (1995) 569-575
    const int N_E704 = 5;
    double pT_E704[N_E704] = { 2.59, 2.79, 2.99, 3.24, 3.59 };
    double y_E704[N_E704] = { 4.9, 1.58, 1.15, 0.539, 0.312 };
    double sy_E704[N_E704] = { 0.7824321056807421, 0.45486261662176636, 0.29966648127543394, 0.17053152201279387, 0.09492101980067429 };
    double ey_E704[N_E704] = { 0.61, 0.38, 0.27, 0.16, 0.091 };
    double E_E704 = 19.4;
    TGraphErrors *tE_DirPh_E704 = new TGraphErrors(N_E704, pT_E704, y_E704, 0, ey_E704);
    TGraphErrors *tS_DirPh_E704 = new TGraphErrors(N_E704, pT_E704, y_E704, 0, sy_E704);
    tE_DirPh_E704               = ScaleGraph(tE_DirPh_E704, 1000);
    tS_DirPh_E704               = ScaleGraph(tS_DirPh_E704, 1000);
    RemoveUpperLimits(tE_DirPh_E704);
    RemoveUpperLimits(tS_DirPh_E704);
    TGraphErrors *tExT_DirPh_E704 = xTScalePhoton(tE_DirPh_E704, E_E704);
    tExT_DirPh_E704->Print();
    //==========================================================================================
    //                                         NA24
    //==========================================================================================
    // Phys.Rev. D36 (1987) 8
    const int N_NA24 = 5;
    double pT_NA24[N_NA24] = { 3.25, 3.75, 4.25, 5.0, 6.0 };
    double y_NA24[N_NA24] = { 3.75E-34, 1.21E-34, 2.5E-35, 5.48E-36, 9.5E-37 };
    double sy_NA24[N_NA24] = { 1.9377564346429095E-34, 4.860041152089147E-35, 8.499999999999999E-36, 1.6278820596099707E-36, 3.9458839313897716E-37 };
    double ey_NA24[N_NA24] = { 9.3E-35, 3.9E-35, 4.0E-36, 1.2E-36, 3.9E-37 };
    double E_NA24 = 23.8;
    TGraphErrors *tE_DirPh_NA24 = new TGraphErrors(N_NA24, pT_NA24, y_NA24, 0, ey_NA24);
    TGraphErrors *tS_DirPh_NA24 = new TGraphErrors(N_NA24, pT_NA24, y_NA24, 0, sy_NA24);
    tE_DirPh_NA24               = ScaleGraph(tE_DirPh_NA24, 1e36);
    tS_DirPh_NA24               = ScaleGraph(tS_DirPh_NA24, 1e36);
    RemoveUpperLimits(tE_DirPh_NA24);
    RemoveUpperLimits(tS_DirPh_NA24);
    //     styleme(tE_DirPh_NA24, 21, kGray);
    TGraphErrors *tExT_DirPh_NA24 = xTScalePhoton(tE_DirPh_NA24, E_NA24);

    //==========================================================================================
    //                                         WA70
    //==========================================================================================
    // Z.Phys. C38 (1988) 371
    const int N_WA70 = 5;
    double pT_WA70[N_WA70] = { 4.11, 4.36, 4.7, 5.2, 5.7 };
    double y_WA70[N_WA70] = { 32.6, 24.0, 10.5, 3.92, 0.68 };
    double sy_WA70[N_WA70] = { 8.5,  5.7 , 2.1, 0.94, 0.22 };
    double ey_WA70[N_WA70] = { 3.9, 3.3, 1.4, 0.87, 0.32 };
    double E_WA70 = 23.0;
    TGraphErrors *tE_DirPh_WA70 = new TGraphErrors(N_WA70, pT_WA70, y_WA70, 0, ey_WA70);
    TGraphErrors *tS_DirPh_WA70 = new TGraphErrors(N_WA70, pT_WA70, y_WA70, 0, sy_WA70);
    tE_DirPh_WA70               = ScaleGraph(tE_DirPh_WA70, 1e0);
    tS_DirPh_WA70               = ScaleGraph(tS_DirPh_WA70, 1e0);
    RemoveUpperLimits(tE_DirPh_WA70);
    RemoveUpperLimits(tS_DirPh_WA70);

    //     styleme(tE_DirPh_WA70, 22, kGray);
    TGraphErrors *tExT_DirPh_WA70 = xTScalePhoton(tE_DirPh_WA70, E_WA70);

    //==========================================================================================
    //                                         UA6_0
    //==========================================================================================
    // Phys.Lett. B317 (1993) 243-249
    const int N_UA6_0 = 10;
    double pT_UA6_0[N_UA6_0] = { 4.19, 4.39, 4.59, 4.79, 4.99, 5.19, 5.46, 5.89, 6.32,
    7.07 };
    double y_UA6_0[N_UA6_0] = { 56.3, 40.3, 24.1, 16.7, 7.7, 5.2, 2.37, 0.76, 0.44,
    0.0 };
    double ey_UA6_0[N_UA6_0] = { 4.9, 3.7, 2.7, 2.2, 1.5, 1.2, 0.53, 0.29, 0.15,
    0.0 };
    double sy_UA6_0[N_UA6_0] = { 4.9, 3.7, 2.7, 2.2, 1.5, 1.2, 0.53, 0.29, 0.15,
    0.03 };
    double E_UA6_0 = 24.3;
    TGraphErrors *tE_DirPh_UA6_0 = new TGraphErrors(N_UA6_0, pT_UA6_0, y_UA6_0, 0, ey_UA6_0);
    TGraphErrors *tS_DirPh_UA6_0 = new TGraphErrors(N_UA6_0, pT_UA6_0, y_UA6_0, 0, sy_UA6_0);
    RemoveUpperLimits(tS_DirPh_UA6_0);
    RemoveUpperLimits(tE_DirPh_UA6_0);
    //     styleme(tE_DirPh_UA6_0, 33);
    TGraphErrors *tExT_DirPh_UA6_0 = xTScalePhoton(tE_DirPh_UA6_0, E_UA6_0);

    //==========================================================================================
    //                                         UA6_1
    //==========================================================================================
    // Phys.Lett. B436 (1998) 222-230
    const int N_UA6_1 = 6;
    double pT_UA6_1[N_UA6_1] = { 4.199999999999999, 4.4, 4.6, 4.800000000000001, 5.1, 5.699999999999999 };
    double y_UA6_1[N_UA6_1] = { 119.1, 66.6, 44.6, 26.8, 7.7, 5.4 };
    double sy_UA6_1[N_UA6_1] = { 28.207268566807386, 18.11104635298579, 12.70944530654269, 9.033825324855467, 3.6999999999999997, 1.6124515496597098 };
    double ey_UA6_1[N_UA6_1] = { 21.8, 15.1, 10.8, 8.1, 3.5, 1.4 };
    double E_UA6_1 = 24.3;

    TGraphErrors *tE_DirPh_UA6_1 = new TGraphErrors(N_UA6_1, pT_UA6_1, y_UA6_1, 0, ey_UA6_1);
    TGraphErrors *tS_DirPh_UA6_1 = new TGraphErrors(N_UA6_1, pT_UA6_1, y_UA6_1, 0, sy_UA6_1);
    RemoveUpperLimits(tS_DirPh_UA6_1);
    RemoveUpperLimits(tE_DirPh_UA6_1);
    //     styleme(tE_DirPh_UA6_1, 34);
    TGraphErrors *tExT_DirPh_UA6_1 = xTScalePhoton(tE_DirPh_UA6_1, E_UA6_1);



    TFile fileTheoryGraphs("ExternalInput/OtherExperiments/DataCompilationFromOtherEnergiesPP.root","UPDATE");
        tE_DirPh_ATLAS->Write("ATLAS_tot_Gamma_7TeV", TObject::kOverwrite);
        tTot_7TeV_DirPh_ATLAS->Write("ATLAS_1_tot_Gamma_7TeV", TObject::kOverwrite);
        tTot_8TeV_DirPh_ATLAS->Write("ATLAS_tot_Gamma_8TeV", TObject::kOverwrite);
        tTot_13TeV_DirPh_ATLAS->Write("ATLAS_tot_Gamma_13TeV", TObject::kOverwrite);
        tTot_2760GeV_DirPh_CMS->Write("CMS_tot_Gamma_2.76TeV", TObject::kOverwrite);
        tTot_7TeV_IsoPh_ALICE->Write("ALICE_tot_Gamma_7TeV", TObject::kOverwrite);
        tE_DirPh_CMS->Write("CMS_tot_Gamma_7TeV", TObject::kOverwrite);
        tTot_7TeV_DirPh_CMS->Write("CMS_1_tot_Gamma_7TeV", TObject::kOverwrite);
        tE_DirPh_CDF->Write("CDF_tot_Gamma_1.8TeV", TObject::kOverwrite);
        tE_DirPh_D0->Write("D0_tot_Gamma_1.8TeV", TObject::kOverwrite);
        tE_DirPh_UA1_0->Write("UA1_tot_Gamma_630GeV", TObject::kOverwrite);
        tE_DirPh_UA1_1->Write("UA1_tot_Gamma_546GeV", TObject::kOverwrite);
        tE_DirPh_UA2->Write("UA2_tot_Gamma_630GeV", TObject::kOverwrite);
        tE_DirPh_PHENIX->Write("PHENIX_tot_Gamma_200GeV", TObject::kOverwrite);
        tE_DirPh_R807->Write("R807_tot_Gamma_63GeV", TObject::kOverwrite);
        tE_DirPh_R108->Write("R807_tot_Gamma_62.4GeV", TObject::kOverwrite);
        tE_DirPh_E706->Write("E706_tot_Gamma_38.8GeV", TObject::kOverwrite);
        tE_DirPh_E706_0->Write("E706_tot_Gamma_31.6GeV", TObject::kOverwrite);
        tE_DirPh_E704->Write("E704_tot_Gamma_19.4GeV", TObject::kOverwrite);
        tE_DirPh_NA24->Write("NA24_tot_Gamma_23.8GeV", TObject::kOverwrite);
        tE_DirPh_WA70->Write("WA70_tot_Gamma_23GeV", TObject::kOverwrite);
        tE_DirPh_UA6_0->Write("UA6_0_tot_Gamma_24.3GeV", TObject::kOverwrite);
        tE_DirPh_UA6_1->Write("UA6_1_tot_Gamma_24.3GeV", TObject::kOverwrite);

        tExT_DirPh_ATLAS->Write("ATLAS_tot_Gamma_7TeV_xT", TObject::kOverwrite);
        tTotxT_7TeV_DirPh_ATLAS->Write("ATLAS_1_tot_Gamma_7TeV_xT", TObject::kOverwrite);
        tTotxT_8TeV_DirPh_ATLAS->Write("ATLAS_tot_Gamma_8TeV_xT", TObject::kOverwrite);
        tTotxT_13TeV_DirPh_ATLAS->Write("ATLAS_tot_Gamma_13TeV_xT", TObject::kOverwrite);
        tExT_DirPh_CMS->Write("CMS_tot_Gamma_7TeV_xT", TObject::kOverwrite);
        tTotxT_7TeV_IsoPh_ALICE->Write("ALICE_tot_Gamma_7TeV_xT", TObject::kOverwrite);
        tTotxT_7TeV_DirPh_CMS->Write("CMS_1_tot_Gamma_7TeV_xT", TObject::kOverwrite);
        tTotxT_2760GeV_DirPh_CMS->Write("CMS_tot_Gamma_2.76TeV_xT", TObject::kOverwrite);
        tExT_DirPh_CDF->Write("CDF_tot_Gamma_1.8TeV_xT", TObject::kOverwrite);
        tExT_DirPh_D0->Write("D0_tot_Gamma_1.8TeV_xT", TObject::kOverwrite);
        tExT_DirPh_UA1_0->Write("UA1_tot_Gamma_630GeV_xT", TObject::kOverwrite);
        tExT_DirPh_UA1_1->Write("UA1_tot_Gamma_546GeV_xT", TObject::kOverwrite);
        tExT_DirPh_UA2->Write("UA2_tot_Gamma_630GeV_xT", TObject::kOverwrite);
        tExT_DirPh_PHENIX->Write("PHENIX_tot_Gamma_200GeV_xT", TObject::kOverwrite);
        tExT_DirPh_R807->Write("R807_tot_Gamma_63GeV_xT", TObject::kOverwrite);
        tExT_DirPh_R108->Write("R807_tot_Gamma_62.4GeV_xT", TObject::kOverwrite);
        tExT_DirPh_E706->Write("E706_tot_Gamma_38.8GeV_xT", TObject::kOverwrite);
        tExT_DirPh_E706_0->Write("E706_tot_Gamma_31.6GeV_xT", TObject::kOverwrite);
        tExT_DirPh_E704->Write("E704_tot_Gamma_19.4GeV_xT", TObject::kOverwrite);
        tExT_DirPh_NA24->Write("NA24_tot_Gamma_23.8GeV_xT", TObject::kOverwrite);
        tExT_DirPh_WA70->Write("WA70_tot_Gamma_23GeV_xT", TObject::kOverwrite);
        tExT_DirPh_UA6_0->Write("UA6_0_tot_Gamma_24.3GeV_xT", TObject::kOverwrite);
        tExT_DirPh_UA6_1->Write("UA6_1_tot_Gamma_24.3GeV_xT", TObject::kOverwrite);

        RGammaStat_7TeV_LMEE_ALICE->Write("RGamma_stat_LMEE_7TeV", TObject::kOverwrite);
        RGammaSys_7TeV_LMEE_ALICE->Write("RGamma_sys_LMEE_7TeV", TObject::kOverwrite);
        RGammaTot_7TeV_LMEE_ALICE->Write("RGamma_tot_LMEE_7TeV", TObject::kOverwrite);

    fileTheoryGraphs.Close();
}
