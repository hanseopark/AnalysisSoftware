
/****************************************************************************************************************************
******      provided by Gamma Conversion Group, PWG4,                                        *****
******      Ana Marin, marin@physi.uni-heidelberg.de                                      *****
******         Kathrin Koch, kkoch@physi.uni-heidelberg.de                                      *****
******      Friederike Bock, friederike.bock@cern.ch                                      *****
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
#include "CommonHeaders/CombinationFunctions.h"

extern TRandom*   gRandom;
extern TBenchmark*   gBenchmark;
extern TSystem*   gSystem;
extern TMinuit*   gMinuit;

void CompareChargedAndNeutralPionDataALICEpPb(TString outputDir = "pdf/CombineMesonMeasurementspPbX", TString suffix = "pdf"){

   gROOT->Reset();   
   gROOT->SetStyle("Plain");
   
   StyleSettingsThesis();  
   SetPlotStyle();
   
   Double_t xSection2760GeVpp =     55.416*1e-3;
   Double_t xSection2760GeVErrpp =  3.9;
   Double_t xSection2760GeVppINEL = 62.8*1e9;
   Double_t xSection900GeVppINEL = 52.5*1e9;
   Double_t xSection7TeVppINEL = 73.2*1e9;   
   Double_t recalcBarn =         1e12; //NLO in pbarn!!!!

   gSystem->Exec("mkdir -p "+outputDir);
   Color_t  colorComb0005           = kRed+1;
   Color_t  colorComb0510           = 807;
   Color_t  colorComb1020           = 800;
   Color_t  colorComb2040           = kGreen+2;
   Color_t  colorComb4060           = kCyan+2;
   Color_t  colorComb6080           = kBlue+1;
   Color_t  colorComb60100           = kBlue+3;

   Style_t  markerStyleCommmonSpectrum0005   = 20 ;
   Style_t  markerStyleCommmonSpectrum0510   = 21 ;
   Style_t  markerStyleCommmonSpectrum1020   = 29 ;
   Style_t  markerStyleCommmonSpectrum0020   = 20 ;
   Style_t  markerStyleCommmonSpectrum2040   = 33 ;
   Style_t  markerStyleCommmonSpectrum4060   = 20 ;
   Style_t  markerStyleCommmonSpectrum6080   = 29 ;
   Style_t  markerStyleCommmonSpectrum60100   = 21 ;

   Size_t   markerSizeCommonSpectrum0005  = 2.;
   Size_t   markerSizeCommonSpectrum0510  = 2.;
   Size_t   markerSizeCommonSpectrum1020  = 2.5;
   Size_t   markerSizeCommonSpectrum2040  = 2.5;
   Size_t   markerSizeCommonSpectrum4060  = 2.;
   Size_t   markerSizeCommonSpectrum6080  = 2.;
   Size_t   markerSizeCommonSpectrum60100  = 2.;
   
   Width_t  widthLinesBoxes;

   TString collisionSystem0020 = "0-20% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";     
   TString collisionSystem2040 = "20-40% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";      
   TString collisionSystem4060 = "40-60% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";    
   TString collisionSystem6080 = "60-80% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";     
   TString collisionSystem60100 = "60-100% p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";     
   TString collisionSystempPb = "p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";     
   
   Size_t markerSizeComparison = 1.5;
   TString nameHistoPCM = "CorrectedYieldPi0";
   TString nameGraphPCM = "Pi0SystError";
   
    TFile* fileNeutralPionPCMDatapPbForwardRap = new TFile("data_PCMResults_pPb_fwd_Legotrain-v5-05-20-AN_allcent.root");
   //   TFile* fileNeutralPionPCMDatapPbForwardRap = new TFile("data_GammaConversionResultsFullCorrection_pPb_ForRap_CentBinning_Lego.root");
   //   TFile* fileNeutralPionPCMDatapPbForwardRap = new TFile("data_PCMResults_pPb_20130618_fulleta.root");
   TDirectory* directoryPCMForwardPi0pPb =             (TDirectory*)fileNeutralPionPCMDatapPbForwardRap->Get("Pi0_pPb_5.023TeV_0-100%"); 
   TH1D* histoPCMForwardYieldPi0pPb =            (TH1D*)directoryPCMForwardPi0pPb->Get(nameHistoPCM.Data());
   TGraphAsymmErrors* graphPCMForwardYieldPi0SysErrpPb=    (TGraphAsymmErrors*)directoryPCMForwardPi0pPb->Get(nameGraphPCM.Data()); 
   TDirectory* directoryPCMForwardPi0pPb0020 =            (TDirectory*)fileNeutralPionPCMDatapPbForwardRap->Get("Pi0_pPb_5.023TeV_0-20%"); 
   TH1D* histoPCMForwardYieldPi0pPb0020 =           (TH1D*)directoryPCMForwardPi0pPb0020->Get(nameHistoPCM.Data());
   TGraphAsymmErrors* graphPCMForwardYieldPi0SysErrpPb0020=   (TGraphAsymmErrors*)directoryPCMForwardPi0pPb0020->Get(nameGraphPCM.Data());   
   TDirectory* directoryPCMForwardPi0pPb2040 =            (TDirectory*)fileNeutralPionPCMDatapPbForwardRap->Get("Pi0_pPb_5.023TeV_20-40%"); 
   TH1D* histoPCMForwardYieldPi0pPb2040 =           (TH1D*)directoryPCMForwardPi0pPb2040->Get(nameHistoPCM.Data());
   TGraphAsymmErrors* graphPCMForwardYieldPi0SysErrpPb2040=   (TGraphAsymmErrors*)directoryPCMForwardPi0pPb2040->Get(nameGraphPCM.Data());   
   TDirectory* directoryPCMForwardPi0pPb4060 =            (TDirectory*)fileNeutralPionPCMDatapPbForwardRap->Get("Pi0_pPb_5.023TeV_40-60%"); 
   TH1D* histoPCMForwardYieldPi0pPb4060 =           (TH1D*)directoryPCMForwardPi0pPb4060->Get(nameHistoPCM.Data());
   TGraphAsymmErrors* graphPCMForwardYieldPi0SysErrpPb4060=   (TGraphAsymmErrors*)directoryPCMForwardPi0pPb4060->Get(nameGraphPCM.Data());   
   TDirectory* directoryPCMForwardPi0pPb6080 =            (TDirectory*)fileNeutralPionPCMDatapPbForwardRap->Get("Pi0_pPb_5.023TeV_60-80%"); 
   TH1D* histoPCMForwardYieldPi0pPb6080 =           (TH1D*)directoryPCMForwardPi0pPb6080->Get(nameHistoPCM.Data());
   TGraphAsymmErrors* graphPCMForwardYieldPi0SysErrpPb6080=   (TGraphAsymmErrors*)directoryPCMForwardPi0pPb6080->Get(nameGraphPCM.Data());   
   TDirectory* directoryPCMForwardPi0pPb60100 =            (TDirectory*)fileNeutralPionPCMDatapPbForwardRap->Get("Pi0_pPb_5.023TeV_60-100%"); 
   TH1D* histoPCMForwardYieldPi0pPb60100 =           (TH1D*)directoryPCMForwardPi0pPb60100->Get(nameHistoPCM.Data());
   TGraphAsymmErrors* graphPCMForwardYieldPi0SysErrpPb60100=   (TGraphAsymmErrors*)directoryPCMForwardPi0pPb60100->Get(nameGraphPCM.Data());   
   
//   TFile* fileNeutralPionPCMDatapPbMidRap = new TFile("data_GammaConversionResultsFullCorrection_pPb_MidRap_CentBinning_Lego.root");
//   TFile* fileNeutralPionPCMDatapPbMidRap = new TFile("data_PCMResults_pPb_20130618_midrap.root");
   TFile* fileNeutralPionPCMDatapPbMidRap = new TFile("data_PCMResults_pPb_mid_Legotrain-v5-05-20-AN_allcent.root");
   TDirectory* directoryPCMMidPi0pPb =             (TDirectory*)fileNeutralPionPCMDatapPbMidRap->Get("Pi0_pPb_5.023TeV_0-100%"); 
   TH1D* histoPCMMidYieldPi0pPb =            (TH1D*)directoryPCMMidPi0pPb->Get(nameHistoPCM.Data());
   TGraphAsymmErrors* graphPCMMidYieldPi0SysErrpPb=    (TGraphAsymmErrors*)directoryPCMMidPi0pPb->Get(nameGraphPCM.Data()); 
   TDirectory* directoryPCMMidPi0pPb0020 =            (TDirectory*)fileNeutralPionPCMDatapPbMidRap->Get("Pi0_pPb_5.023TeV_0-20%"); 
   TH1D* histoPCMMidYieldPi0pPb0020 =           (TH1D*)directoryPCMMidPi0pPb0020->Get(nameHistoPCM.Data());
   TGraphAsymmErrors* graphPCMMidYieldPi0SysErrpPb0020=   (TGraphAsymmErrors*)directoryPCMMidPi0pPb0020->Get(nameGraphPCM.Data());   
   TDirectory* directoryPCMMidPi0pPb2040 =            (TDirectory*)fileNeutralPionPCMDatapPbMidRap->Get("Pi0_pPb_5.023TeV_20-40%"); 
   TH1D* histoPCMMidYieldPi0pPb2040 =           (TH1D*)directoryPCMMidPi0pPb2040->Get(nameHistoPCM.Data());
   TGraphAsymmErrors* graphPCMMidYieldPi0SysErrpPb2040=   (TGraphAsymmErrors*)directoryPCMMidPi0pPb2040->Get(nameGraphPCM.Data());   
   TDirectory* directoryPCMMidPi0pPb4060 =            (TDirectory*)fileNeutralPionPCMDatapPbMidRap->Get("Pi0_pPb_5.023TeV_40-60%"); 
   TH1D* histoPCMMidYieldPi0pPb4060 =           (TH1D*)directoryPCMMidPi0pPb4060->Get(nameHistoPCM.Data());
   TGraphAsymmErrors* graphPCMMidYieldPi0SysErrpPb4060=   (TGraphAsymmErrors*)directoryPCMMidPi0pPb4060->Get(nameGraphPCM.Data());   
   TDirectory* directoryPCMMidPi0pPb6080 =            (TDirectory*)fileNeutralPionPCMDatapPbMidRap->Get("Pi0_pPb_5.023TeV_60-80%"); 
   TH1D* histoPCMMidYieldPi0pPb6080 =           (TH1D*)directoryPCMMidPi0pPb6080->Get(nameHistoPCM.Data());
   TGraphAsymmErrors* graphPCMMidYieldPi0SysErrpPb6080=   (TGraphAsymmErrors*)directoryPCMMidPi0pPb6080->Get(nameGraphPCM.Data());   
   TDirectory* directoryPCMMidPi0pPb60100 =            (TDirectory*)fileNeutralPionPCMDatapPbMidRap->Get("Pi0_pPb_5.023TeV_60-100%"); 
   TH1D* histoPCMMidYieldPi0pPb60100 =           (TH1D*)directoryPCMMidPi0pPb60100->Get(nameHistoPCM.Data());
   TGraphAsymmErrors* graphPCMMidYieldPi0SysErrpPb60100=   (TGraphAsymmErrors*)directoryPCMMidPi0pPb60100->Get(nameGraphPCM.Data());   
   
//    TFile*   filePHOSpPb =       new TFile("ExternalInputpPb/PHOS_pPb_08042013.root");
//    TH1D* histoPHOSYieldPi0pPb = (TH1D*)filePHOSpPb->Get("yeild1_GS_Disp_cen0");
//    TFile*   filePHOSpPb =       new TFile("ExternalInputpPb/data_PHOSResultsFullCorrection_pPb_20130826.root");
    TFile*   filePHOSpPb =       new TFile("ExternalInputpPb/data_PHOSResultsFullCorrection_pPb-21112013.root");
   TDirectory *directoryPHOSPi0pPb = 	(TDirectory*)filePHOSpPb->Get("Pi0_pPb_5.023TeV_0-100%"); 
   TH1D* histoPHOSYieldPi0pPb = 	(TH1D*)directoryPHOSPi0pPb->Get(nameHistoPCM.Data());      
 
   
//    TFile* fileChargedPionpPb = new TFile("ExternalInputpPb/ChargedPionSpectrapPb_10_Apr_2013.root");
   
//    TH1D* histoChargedPionSpecLowPtSyspPb= (TH1D*)fileChargedPionpPb->Get("histoChargedPionSpecLowPtSyspPb");
//    TH1D* histoChargedPionSpecLowPtStatpPb= (TH1D*)fileChargedPionpPb->Get("histoChargedPionSpecLowPtStatpPb");
//    TH1D* histoChargedPionSpecLowPtSyspPb0020= (TH1D*)fileChargedPionpPb->Get("histoChargedPionSpecLowPtSyspPb0020");
//    TH1D* histoChargedPionSpecLowPtStatpPb0020= (TH1D*)fileChargedPionpPb->Get("histoChargedPionSpecLowPtStatpPb0020");
//    TH1D* histoChargedPionSpecLowPtSyspPb2040= (TH1D*)fileChargedPionpPb->Get("histoChargedPionSpecLowPtSyspPb2040");
//    TH1D* histoChargedPionSpecLowPtStatpPb2040= (TH1D*)fileChargedPionpPb->Get("histoChargedPionSpecLowPtStatpPb2040");
//    TH1D* histoChargedPionSpecLowPtSyspPb4060= (TH1D*)fileChargedPionpPb->Get("histoChargedPionSpecLowPtSyspPb4060");
//    TH1D* histoChargedPionSpecLowPtStatpPb4060= (TH1D*)fileChargedPionpPb->Get("histoChargedPionSpecLowPtStatpPb4060");
//    TH1D* histoChargedPionSpecLowPtSyspPb6080= (TH1D*)fileChargedPionpPb->Get("histoChargedPionSpecLowPtSyspPb6080");
//    TH1D* histoChargedPionSpecLowPtStatpPb6080= (TH1D*)fileChargedPionpPb->Get("histoChargedPionSpecLowPtStatpPb6080");
//    TH1D* histoChargedPionSpecLowPtSyspPb80100= (TH1D*)fileChargedPionpPb->Get("histoChargedPionSpecLowPtSyspPb80100");
//    TH1D* histoChargedPionSpecLowPtStatpPb80100= (TH1D*)fileChargedPionpPb->Get("histoChargedPionSpecLowPtStatpPb80100");
   
   TFile* fileChargedPionpPb = new TFile("ExternalInputpPb/20130723_CombinedSpectra_pA_ITSsa_TPCTOF_TOF.root");
   // positive pions  
   TH1D* histoChargedPionSpecLowPtSyspPbcent0= (TH1D*)fileChargedPionpPb->Get("sys_cent0_pion_plus"); //0-5% 
   TH1D* histoChargedPionSpecLowPtStatpPbcent0= (TH1D*)fileChargedPionpPb->Get("stat_cent0_pion_plus");

   TH1D* histoChargedPionSpecLowPtSyspPbcent1= (TH1D*)fileChargedPionpPb->Get("sys_cent1_pion_plus"); //5-10%
   TH1D* histoChargedPionSpecLowPtStatpPbcent1= (TH1D*)fileChargedPionpPb->Get("stat_cent1_pion_plus");

   TH1D* histoChargedPionSpecLowPtSyspPbcent2= (TH1D*)fileChargedPionpPb->Get("sys_cent2_pion_plus");//10-20%
   TH1D* histoChargedPionSpecLowPtStatpPbcent2= (TH1D*)fileChargedPionpPb->Get("stat_cent2_pion_plus");

   TH1D* histoChargedPionSpecLowPtSyspPbcent3= (TH1D*)fileChargedPionpPb->Get("sys_cent3_pion_plus");//20-40%
   TH1D* histoChargedPionSpecLowPtStatpPbcent3= (TH1D*)fileChargedPionpPb->Get("stat_cent3_pion_plus");

   TH1D* histoChargedPionSpecLowPtSyspPbcent4= (TH1D*)fileChargedPionpPb->Get("sys_cent4_pion_plus");//40-60%
   TH1D* histoChargedPionSpecLowPtStatpPbcent4= (TH1D*)fileChargedPionpPb->Get("stat_cent4_pion_plus");

   TH1D* histoChargedPionSpecLowPtSyspPbcent5= (TH1D*)fileChargedPionpPb->Get("sys_cent5_pion_plus");//60-80%
   TH1D* histoChargedPionSpecLowPtStatpPbcent5= (TH1D*)fileChargedPionpPb->Get("stat_cent5_pion_plus");

   TH1D* histoChargedPionSpecLowPtSyspPbcent6= (TH1D*)fileChargedPionpPb->Get("sys_cent6_pion_plus");//80-100%
   TH1D* histoChargedPionSpecLowPtStatpPbcent6= (TH1D*)fileChargedPionpPb->Get("stat_cent6_pion_plus");
   // negative pions
   TH1D* histoChargedPionSpecLowPtSyspPbnegcent0= (TH1D*)fileChargedPionpPb->Get("sys_cent0_pion_minus"); //0-5% 
   TH1D* histoChargedPionSpecLowPtStatpPbnegcent0= (TH1D*)fileChargedPionpPb->Get("stat_cent0_pion_minus");

   TH1D* histoChargedPionSpecLowPtSyspPbnegcent1= (TH1D*)fileChargedPionpPb->Get("sys_cent1_pion_minus"); //5-10%
   TH1D* histoChargedPionSpecLowPtStatpPbnegcent1= (TH1D*)fileChargedPionpPb->Get("stat_cent1_pion_minus");

   TH1D* histoChargedPionSpecLowPtSyspPbnegcent2= (TH1D*)fileChargedPionpPb->Get("sys_cent2_pion_minus");//10-20%
   TH1D* histoChargedPionSpecLowPtStatpPbnegcent2= (TH1D*)fileChargedPionpPb->Get("stat_cent2_pion_minus");

   TH1D* histoChargedPionSpecLowPtSyspPbnegcent3= (TH1D*)fileChargedPionpPb->Get("sys_cent3_pion_minus");//20-40%
   TH1D* histoChargedPionSpecLowPtStatpPbnegcent3= (TH1D*)fileChargedPionpPb->Get("stat_cent3_pion_minus");

   TH1D* histoChargedPionSpecLowPtSyspPbnegcent4= (TH1D*)fileChargedPionpPb->Get("sys_cent4_pion_minus");//40-60%
   TH1D* histoChargedPionSpecLowPtStatpPbnegcent4= (TH1D*)fileChargedPionpPb->Get("stat_cent4_pion_minus");

   TH1D* histoChargedPionSpecLowPtSyspPbnegcent5= (TH1D*)fileChargedPionpPb->Get("sys_cent5_pion_minus");//60-80%
   TH1D* histoChargedPionSpecLowPtStatpPbnegcent5= (TH1D*)fileChargedPionpPb->Get("stat_cent5_pion_minus");

   TH1D* histoChargedPionSpecLowPtSyspPbnegcent6= (TH1D*)fileChargedPionpPb->Get("sys_cent6_pion_minus");//80-100%
   TH1D* histoChargedPionSpecLowPtStatpPbnegcent6= (TH1D*)fileChargedPionpPb->Get("stat_cent6_pion_minus");

   histoChargedPionSpecLowPtSyspPbcent0->Add(histoChargedPionSpecLowPtSyspPbnegcent0);
   histoChargedPionSpecLowPtStatpPbcent0->Add(histoChargedPionSpecLowPtStatpPbnegcent0);

   histoChargedPionSpecLowPtSyspPbcent1->Add(histoChargedPionSpecLowPtSyspPbnegcent1);
   histoChargedPionSpecLowPtStatpPbcent1->Add(histoChargedPionSpecLowPtStatpPbnegcent1);

   histoChargedPionSpecLowPtSyspPbcent2->Add(histoChargedPionSpecLowPtSyspPbnegcent2);
   histoChargedPionSpecLowPtStatpPbcent2->Add(histoChargedPionSpecLowPtStatpPbnegcent2);

   histoChargedPionSpecLowPtSyspPbcent3->Add(histoChargedPionSpecLowPtSyspPbnegcent3);
   histoChargedPionSpecLowPtStatpPbcent3->Add(histoChargedPionSpecLowPtStatpPbnegcent3);

   histoChargedPionSpecLowPtSyspPbcent4->Add(histoChargedPionSpecLowPtSyspPbnegcent4);
   histoChargedPionSpecLowPtStatpPbcent4->Add(histoChargedPionSpecLowPtStatpPbnegcent4);

   histoChargedPionSpecLowPtSyspPbcent5->Add(histoChargedPionSpecLowPtSyspPbnegcent5);
   histoChargedPionSpecLowPtStatpPbcent5->Add(histoChargedPionSpecLowPtStatpPbnegcent5);

   histoChargedPionSpecLowPtSyspPbcent6->Add(histoChargedPionSpecLowPtSyspPbnegcent6);
   histoChargedPionSpecLowPtStatpPbcent6->Add(histoChargedPionSpecLowPtStatpPbnegcent6);

   histoChargedPionSpecLowPtSyspPbcent0->Sumw2();
   histoChargedPionSpecLowPtStatpPbcent0->Sumw2();
   histoChargedPionSpecLowPtSyspPbcent1->Sumw2();
   histoChargedPionSpecLowPtStatpPbcent1->Sumw2();
   histoChargedPionSpecLowPtSyspPbcent2->Sumw2();
   histoChargedPionSpecLowPtStatpPbcent2->Sumw2();
   histoChargedPionSpecLowPtSyspPbcent3->Sumw2();
   histoChargedPionSpecLowPtStatpPbcent3->Sumw2();
   histoChargedPionSpecLowPtSyspPbcent4->Sumw2();
   histoChargedPionSpecLowPtStatpPbcent4->Sumw2();
   histoChargedPionSpecLowPtSyspPbcent5->Sumw2();
   histoChargedPionSpecLowPtStatpPbcent5->Sumw2();
   histoChargedPionSpecLowPtSyspPbcent6->Sumw2();
   histoChargedPionSpecLowPtStatpPbcent6->Sumw2();
   // divide pos and neg pions by 2
   histoChargedPionSpecLowPtSyspPbcent0->Scale(0.5);
   histoChargedPionSpecLowPtStatpPbcent0->Scale(0.5);
   histoChargedPionSpecLowPtSyspPbcent1->Scale(0.5);
   histoChargedPionSpecLowPtStatpPbcent1->Scale(0.5);
   histoChargedPionSpecLowPtSyspPbcent2->Scale(0.5);
   histoChargedPionSpecLowPtStatpPbcent2->Scale(0.5);
   histoChargedPionSpecLowPtSyspPbcent3->Scale(0.5);
   histoChargedPionSpecLowPtStatpPbcent3->Scale(0.5);
   histoChargedPionSpecLowPtSyspPbcent4->Scale(0.5);
   histoChargedPionSpecLowPtStatpPbcent4->Scale(0.5);
   histoChargedPionSpecLowPtSyspPbcent5->Scale(0.5);
   histoChargedPionSpecLowPtStatpPbcent5->Scale(0.5);
   histoChargedPionSpecLowPtSyspPbcent6->Scale(0.5);
   histoChargedPionSpecLowPtStatpPbcent6->Scale(0.5);

   for (int i=0; i<histoChargedPionSpecLowPtSyspPbcent0->GetNbinsX();i++){ // divide by pT

     histoChargedPionSpecLowPtSyspPbcent0->SetBinContent(i+1,(histoChargedPionSpecLowPtSyspPbcent0->GetBinContent(i+1)/histoChargedPionSpecLowPtSyspPbcent0->GetBinCenter(i+1)));
     histoChargedPionSpecLowPtSyspPbcent1->SetBinContent(i+1,(histoChargedPionSpecLowPtSyspPbcent1->GetBinContent(i+1)/histoChargedPionSpecLowPtSyspPbcent1->GetBinCenter(i+1)));
     histoChargedPionSpecLowPtSyspPbcent2->SetBinContent(i+1,(histoChargedPionSpecLowPtSyspPbcent2->GetBinContent(i+1)/histoChargedPionSpecLowPtSyspPbcent2->GetBinCenter(i+1)));
     histoChargedPionSpecLowPtSyspPbcent3->SetBinContent(i+1,(histoChargedPionSpecLowPtSyspPbcent3->GetBinContent(i+1)/histoChargedPionSpecLowPtSyspPbcent3->GetBinCenter(i+1)));
     histoChargedPionSpecLowPtSyspPbcent4->SetBinContent(i+1,(histoChargedPionSpecLowPtSyspPbcent4->GetBinContent(i+1)/histoChargedPionSpecLowPtSyspPbcent4->GetBinCenter(i+1)));
     histoChargedPionSpecLowPtSyspPbcent5->SetBinContent(i+1,(histoChargedPionSpecLowPtSyspPbcent5->GetBinContent(i+1)/histoChargedPionSpecLowPtSyspPbcent5->GetBinCenter(i+1)));
     histoChargedPionSpecLowPtSyspPbcent6->SetBinContent(i+1,(histoChargedPionSpecLowPtSyspPbcent6->GetBinContent(i+1)/histoChargedPionSpecLowPtSyspPbcent6->GetBinCenter(i+1)));

     histoChargedPionSpecLowPtStatpPbcent0->SetBinContent(i+1,(histoChargedPionSpecLowPtStatpPbcent0->GetBinContent(i+1)/histoChargedPionSpecLowPtStatpPbcent0->GetBinCenter(i+1)));
     histoChargedPionSpecLowPtStatpPbcent1->SetBinContent(i+1,(histoChargedPionSpecLowPtStatpPbcent1->GetBinContent(i+1)/histoChargedPionSpecLowPtStatpPbcent1->GetBinCenter(i+1)));
     histoChargedPionSpecLowPtStatpPbcent2->SetBinContent(i+1,(histoChargedPionSpecLowPtStatpPbcent2->GetBinContent(i+1)/histoChargedPionSpecLowPtStatpPbcent2->GetBinCenter(i+1)));
     histoChargedPionSpecLowPtStatpPbcent3->SetBinContent(i+1,(histoChargedPionSpecLowPtStatpPbcent3->GetBinContent(i+1)/histoChargedPionSpecLowPtStatpPbcent3->GetBinCenter(i+1)));
     histoChargedPionSpecLowPtStatpPbcent4->SetBinContent(i+1,(histoChargedPionSpecLowPtStatpPbcent4->GetBinContent(i+1)/histoChargedPionSpecLowPtStatpPbcent4->GetBinCenter(i+1)));
     histoChargedPionSpecLowPtStatpPbcent5->SetBinContent(i+1,(histoChargedPionSpecLowPtStatpPbcent5->GetBinContent(i+1)/histoChargedPionSpecLowPtStatpPbcent5->GetBinCenter(i+1)));
     histoChargedPionSpecLowPtStatpPbcent6->SetBinContent(i+1,(histoChargedPionSpecLowPtStatpPbcent6->GetBinContent(i+1)/histoChargedPionSpecLowPtStatpPbcent6->GetBinCenter(i+1)));

   }
   //0-100%
   TH1D *histoChargedPionSpecLowPtSyspPb=(TH1D*)histoChargedPionSpecLowPtSyspPbcent0->Clone() ;
   TH1D *histoChargedPionSpecLowPtStatpPb =(TH1D*)histoChargedPionSpecLowPtStatpPbcent0->Clone();

   histoChargedPionSpecLowPtSyspPb->Add(histoChargedPionSpecLowPtSyspPbcent0,histoChargedPionSpecLowPtSyspPbcent1,0.05,0.05); //weight with the centrality bin width
   histoChargedPionSpecLowPtSyspPb->Add(histoChargedPionSpecLowPtSyspPbcent2,0.1);
   histoChargedPionSpecLowPtSyspPb->Add(histoChargedPionSpecLowPtSyspPbcent3,0.2);
   histoChargedPionSpecLowPtSyspPb->Add(histoChargedPionSpecLowPtSyspPbcent4,0.2);
   histoChargedPionSpecLowPtSyspPb->Add(histoChargedPionSpecLowPtSyspPbcent5,0.2);
   histoChargedPionSpecLowPtSyspPb->Add(histoChargedPionSpecLowPtSyspPbcent6,0.2);
   
   histoChargedPionSpecLowPtStatpPb->Add(histoChargedPionSpecLowPtStatpPbcent0,histoChargedPionSpecLowPtStatpPbcent1,0.05,0.05);
   histoChargedPionSpecLowPtStatpPb->Add(histoChargedPionSpecLowPtStatpPbcent2,0.1);
   histoChargedPionSpecLowPtStatpPb->Add(histoChargedPionSpecLowPtStatpPbcent3,0.2);
   histoChargedPionSpecLowPtStatpPb->Add(histoChargedPionSpecLowPtStatpPbcent4,0.2);
   histoChargedPionSpecLowPtStatpPb->Add(histoChargedPionSpecLowPtStatpPbcent5,0.2);
   histoChargedPionSpecLowPtStatpPb->Add(histoChargedPionSpecLowPtStatpPbcent6,0.2);

   histoChargedPionSpecLowPtSyspPb->Scale(1/(2*TMath::Pi())); // divide by 2*pi
   histoChargedPionSpecLowPtStatpPb->Scale(1/(2*TMath::Pi()));
   //0-20%
    TH1D *histoChargedPionSpecLowPtSyspPb0020=(TH1D*)histoChargedPionSpecLowPtSyspPbcent0->Clone() ;
   TH1D *histoChargedPionSpecLowPtStatpPb0020 =(TH1D*)histoChargedPionSpecLowPtStatpPbcent0->Clone();

   histoChargedPionSpecLowPtSyspPb0020->Add(histoChargedPionSpecLowPtSyspPbcent0,histoChargedPionSpecLowPtSyspPbcent1,0.25,0.25); //weight with the centrality bin width
   histoChargedPionSpecLowPtSyspPb0020->Add(histoChargedPionSpecLowPtSyspPbcent2,0.5);
   
   histoChargedPionSpecLowPtStatpPb0020->Add(histoChargedPionSpecLowPtStatpPbcent0,histoChargedPionSpecLowPtStatpPbcent1,0.25,0.25);
   histoChargedPionSpecLowPtStatpPb0020->Add(histoChargedPionSpecLowPtStatpPbcent2,0.5);

   histoChargedPionSpecLowPtSyspPb0020->Scale(1/(2*TMath::Pi())); // divide by 2*pi
   histoChargedPionSpecLowPtStatpPb0020->Scale(1/(2*TMath::Pi())); 
   //20-40%
    TH1D *histoChargedPionSpecLowPtSyspPb2040=(TH1D*)histoChargedPionSpecLowPtSyspPbcent3->Clone() ;
   TH1D *histoChargedPionSpecLowPtStatpPb2040 =(TH1D*)histoChargedPionSpecLowPtStatpPbcent3->Clone();

   histoChargedPionSpecLowPtSyspPb2040->Scale(1/(2*TMath::Pi())); // divide by 2*pi
   histoChargedPionSpecLowPtStatpPb2040->Scale(1/(2*TMath::Pi())); 
   //40-60%
    TH1D *histoChargedPionSpecLowPtSyspPb4060=(TH1D*)histoChargedPionSpecLowPtSyspPbcent4->Clone() ;
   TH1D *histoChargedPionSpecLowPtStatpPb4060 =(TH1D*)histoChargedPionSpecLowPtStatpPbcent4->Clone();

   histoChargedPionSpecLowPtSyspPb4060->Scale(1/(2*TMath::Pi())); // divide by 2*pi
   histoChargedPionSpecLowPtStatpPb4060->Scale(1/(2*TMath::Pi())); 
   //60-80%
    TH1D *histoChargedPionSpecLowPtSyspPb6080=(TH1D*)histoChargedPionSpecLowPtSyspPbcent5->Clone() ;
   TH1D *histoChargedPionSpecLowPtStatpPb6080 =(TH1D*)histoChargedPionSpecLowPtStatpPbcent5->Clone();

   histoChargedPionSpecLowPtSyspPb6080->Scale(1/(2*TMath::Pi())); // divide by 2*pi
   histoChargedPionSpecLowPtStatpPb6080->Scale(1/(2*TMath::Pi())); 

   //60-100%
    TH1D *histoChargedPionSpecLowPtSyspPb60100=(TH1D*)histoChargedPionSpecLowPtSyspPbcent5->Clone() ;
   TH1D *histoChargedPionSpecLowPtStatpPb60100 =(TH1D*)histoChargedPionSpecLowPtStatpPbcent5->Clone();
 
   histoChargedPionSpecLowPtSyspPb60100->Add(histoChargedPionSpecLowPtSyspPbcent5,histoChargedPionSpecLowPtSyspPbcent6,1,1); //weight with the centrality bin width
   histoChargedPionSpecLowPtStatpPb60100->Add(histoChargedPionSpecLowPtStatpPbcent5,histoChargedPionSpecLowPtStatpPbcent6,0.5,0.5);

   histoChargedPionSpecLowPtSyspPb60100->Scale(1/(2*TMath::Pi())); // divide by 2*pi
   histoChargedPionSpecLowPtStatpPb60100->Scale(1/(2*TMath::Pi())); 

   cout << "*************************************************************************"<< endl;  
   cout << "******************************  pPb *************************************"<< endl;
   cout << "*************************************************************************"<< endl;
   
   TGraphAsymmErrors* graphPCMForwardYieldPi0SysErrpPbCopy = (TGraphAsymmErrors*) graphPCMForwardYieldPi0SysErrpPb->Clone("graphPCMForwardYieldPi0SysErrpPbCopy");
   TGraphAsymmErrors* graphPCMForwardYieldPi0SysErrpPb0020Copy = (TGraphAsymmErrors*) graphPCMForwardYieldPi0SysErrpPb0020->Clone("graphPCMForwardYieldPi0SysErrpPb0020Copy");
   TGraphAsymmErrors* graphPCMForwardYieldPi0SysErrpPb2040Copy = (TGraphAsymmErrors*) graphPCMForwardYieldPi0SysErrpPb2040->Clone("graphPCMForwardYieldPi0SysErrpPb2040Copy");
   TGraphAsymmErrors* graphPCMForwardYieldPi0SysErrpPb4060Copy = (TGraphAsymmErrors*) graphPCMForwardYieldPi0SysErrpPb4060->Clone("graphPCMForwardYieldPi0SysErrpPb4060Copy");
    TGraphAsymmErrors* graphPCMForwardYieldPi0SysErrpPb6080Copy = (TGraphAsymmErrors*) graphPCMForwardYieldPi0SysErrpPb6080->Clone("graphPCMForwardYieldPi0SysErrpPb6080Copy");
   TGraphAsymmErrors* graphPCMForwardYieldPi0SysErrpPb60100Copy = (TGraphAsymmErrors*) graphPCMForwardYieldPi0SysErrpPb60100->Clone("graphPCMForwardYieldPi0SysErrpPb80100Copy");
   TGraphAsymmErrors* graphPCMMidYieldPi0SysErrpPbCopy = (TGraphAsymmErrors*) graphPCMMidYieldPi0SysErrpPb->Clone("graphPCMMidYieldPi0SysErrpPbCopy");
   TGraphAsymmErrors* graphPCMMidYieldPi0SysErrpPb0020Copy = (TGraphAsymmErrors*) graphPCMMidYieldPi0SysErrpPb0020->Clone("graphPCMMidYieldPi0SysErrpPb0020Copy");
   TGraphAsymmErrors* graphPCMMidYieldPi0SysErrpPb2040Copy = (TGraphAsymmErrors*) graphPCMMidYieldPi0SysErrpPb2040->Clone("graphPCMMidYieldPi0SysErrpPb2040Copy");
   TGraphAsymmErrors* graphPCMMidYieldPi0SysErrpPb4060Copy = (TGraphAsymmErrors*) graphPCMMidYieldPi0SysErrpPb4060->Clone("graphPCMMidYieldPi0SysErrpPb4060Copy");
    TGraphAsymmErrors* graphPCMMidYieldPi0SysErrpPb6080Copy = (TGraphAsymmErrors*) graphPCMMidYieldPi0SysErrpPb6080->Clone("graphPCMMidYieldPi0SysErrpPb6080Copy");
    TGraphAsymmErrors* graphPCMMidYieldPi0SysErrpPb60100Copy = (TGraphAsymmErrors*) graphPCMMidYieldPi0SysErrpPb60100->Clone("graphPCMMidYieldPi0SysErrpPb80100Copy");
   
   cout << "*************************************************************************"<< endl;  
   cout << "******************************  pPb MinBias *****************************"<< endl;
   cout << "*************************************************************************"<< endl;

   cout << "PCM Spectrum forward - low Pt" << endl;
   TGraphErrors* graphRatioLowPtChargedPionsPCMForwardpPb = CalculateRatioBetweenSpectraWithDifferentBinning(histoPCMForwardYieldPi0pPb, graphPCMForwardYieldPi0SysErrpPbCopy, histoChargedPionSpecLowPtStatpPb, histoChargedPionSpecLowPtSyspPb,  kTRUE,  kTRUE)  ;
   graphRatioLowPtChargedPionsPCMForwardpPb->Print();
   cout << "PCM Spectrum mid - low Pt" << endl;
   TGraphErrors* graphRatioLowPtChargedPionsPCMMidpPb = CalculateRatioBetweenSpectraWithDifferentBinning(histoPCMMidYieldPi0pPb, graphPCMMidYieldPi0SysErrpPbCopy, histoChargedPionSpecLowPtStatpPb, histoChargedPionSpecLowPtSyspPb,  kTRUE,  kTRUE)  ;
   graphRatioLowPtChargedPionsPCMMidpPb->Print();
   
   cout << "PHOS Spectrum forward - low Pt" << endl;
   TGraphErrors* graphRatioLowPtChargedPionsPHOSpPb = CalculateRatioBetweenSpectraWithDifferentBinning(histoPHOSYieldPi0pPb, histoPHOSYieldPi0pPb, histoChargedPionSpecLowPtStatpPb, histoChargedPionSpecLowPtSyspPb,  kTRUE,  kTRUE)  ;
   graphRatioLowPtChargedPionsPHOSpPb->Print();
   
   cout << "*************************************************************************"<< endl;  
   cout << "******************************  pPb 00-20% ******************************"<< endl;
   cout << "*************************************************************************"<< endl;

   TGraphErrors* graphRatioLowPtChargedPionsPCMForwardpPb0020 = CalculateRatioBetweenSpectraWithDifferentBinning(histoPCMForwardYieldPi0pPb0020, graphPCMForwardYieldPi0SysErrpPb0020Copy, histoChargedPionSpecLowPtStatpPb0020, histoChargedPionSpecLowPtSyspPb0020,  kTRUE,  kTRUE)  ;
   graphRatioLowPtChargedPionsPCMForwardpPb0020->Print();
   TGraphErrors* graphRatioLowPtChargedPionsPCMMidpPb0020 = CalculateRatioBetweenSpectraWithDifferentBinning(histoPCMMidYieldPi0pPb0020, graphPCMMidYieldPi0SysErrpPb0020Copy, histoChargedPionSpecLowPtStatpPb0020, histoChargedPionSpecLowPtSyspPb0020,  kTRUE,  kTRUE)  ;
   graphRatioLowPtChargedPionsPCMMidpPb0020->Print();
   
   cout << "*************************************************************************"<< endl;  
   cout << "******************************  pPb 20-40% ******************************"<< endl;
   cout << "*************************************************************************"<< endl;
   
   TGraphErrors* graphRatioLowPtChargedPionsPCMForwardpPb2040 = CalculateRatioBetweenSpectraWithDifferentBinning(histoPCMForwardYieldPi0pPb2040, graphPCMForwardYieldPi0SysErrpPb2040Copy, histoChargedPionSpecLowPtStatpPb2040, histoChargedPionSpecLowPtSyspPb2040,  kTRUE,  kTRUE)  ;
   graphRatioLowPtChargedPionsPCMForwardpPb2040->Print();
   TGraphErrors* graphRatioLowPtChargedPionsPCMMidpPb2040 = CalculateRatioBetweenSpectraWithDifferentBinning(histoPCMMidYieldPi0pPb2040, graphPCMMidYieldPi0SysErrpPb2040Copy, histoChargedPionSpecLowPtStatpPb2040, histoChargedPionSpecLowPtSyspPb2040,  kTRUE,  kTRUE)  ;
   graphRatioLowPtChargedPionsPCMMidpPb2040->Print();
   
   cout << "*************************************************************************"<< endl;  
   cout << "******************************  pPb 40-60% ******************************"<< endl;
   cout << "*************************************************************************"<< endl;

   TGraphErrors* graphRatioLowPtChargedPionsPCMForwardpPb4060 = CalculateRatioBetweenSpectraWithDifferentBinning(histoPCMForwardYieldPi0pPb4060, graphPCMForwardYieldPi0SysErrpPb4060Copy, histoChargedPionSpecLowPtStatpPb4060, histoChargedPionSpecLowPtSyspPb4060,  kTRUE,  kTRUE)  ;
   graphRatioLowPtChargedPionsPCMForwardpPb4060->Print();
   TGraphErrors* graphRatioLowPtChargedPionsPCMMidpPb4060 = CalculateRatioBetweenSpectraWithDifferentBinning(histoPCMMidYieldPi0pPb4060, graphPCMMidYieldPi0SysErrpPb4060Copy, histoChargedPionSpecLowPtStatpPb4060, histoChargedPionSpecLowPtSyspPb4060,  kTRUE,  kTRUE)  ;
   graphRatioLowPtChargedPionsPCMMidpPb4060->Print();

   cout << "*************************************************************************"<< endl;  
   cout << "******************************  pPb 60-80% ******************************"<< endl;
   cout << "*************************************************************************"<< endl;

   TGraphErrors* graphRatioLowPtChargedPionsPCMForwardpPb6080 = CalculateRatioBetweenSpectraWithDifferentBinning(histoPCMForwardYieldPi0pPb6080, graphPCMForwardYieldPi0SysErrpPb6080Copy, histoChargedPionSpecLowPtStatpPb6080, histoChargedPionSpecLowPtSyspPb6080,  kTRUE,  kTRUE)  ;
   graphRatioLowPtChargedPionsPCMForwardpPb6080->Print();
   TGraphErrors* graphRatioLowPtChargedPionsPCMMidpPb6080 = CalculateRatioBetweenSpectraWithDifferentBinning(histoPCMMidYieldPi0pPb6080, graphPCMMidYieldPi0SysErrpPb6080Copy, histoChargedPionSpecLowPtStatpPb6080, histoChargedPionSpecLowPtSyspPb6080,  kTRUE,  kTRUE)  ;
   graphRatioLowPtChargedPionsPCMMidpPb6080->Print();
// 
   cout << "*************************************************************************"<< endl;  
   cout << "******************************  pPb 60-100% ******************************"<< endl;
   cout << "*************************************************************************"<< endl;

   TGraphErrors* graphRatioLowPtChargedPionsPCMForwardpPb60100 = CalculateRatioBetweenSpectraWithDifferentBinning(histoPCMForwardYieldPi0pPb60100, graphPCMForwardYieldPi0SysErrpPb60100Copy, histoChargedPionSpecLowPtStatpPb60100, histoChargedPionSpecLowPtSyspPb60100,  kTRUE,  kTRUE)  ;
   graphRatioLowPtChargedPionsPCMForwardpPb60100->Print();
   TGraphErrors* graphRatioLowPtChargedPionsPCMMidpPb60100 = CalculateRatioBetweenSpectraWithDifferentBinning(histoPCMMidYieldPi0pPb60100, graphPCMMidYieldPi0SysErrpPb60100Copy, histoChargedPionSpecLowPtStatpPb60100, histoChargedPionSpecLowPtSyspPb60100,  kTRUE,  kTRUE)  ;
   graphRatioLowPtChargedPionsPCMMidpPb60100->Print();

   
   //************************************************************************************************************
   //******************************  plotting just minBias individual measurements ******************************
   //************************************************************************************************************

   TCanvas* canvasCompYieldpPbInd = new TCanvas("canvasCompYieldpPbInd","",200,10,700,500);  // gives the page size
   DrawGammaCanvasSettings( canvasCompYieldpPbInd,  0.12, 0.02, 0.02, 0.12);
   
   canvasCompYieldpPbInd->SetLogx();
   TH2F * histo2DCompCombinedRatio2;
   histo2DCompCombinedRatio2 = new TH2F("histo2DCompCombinedRatio2","histo2DCompCombinedRatio2",1000,0.3,20.,1000,0.2,4.   );
   SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatio2, "p_{T} (GeV/c)","#pi^{0}/#pi^{#pm}", 0.05,0.064, 0.05,0.06, 0.8,0.9, 512, 505);
   histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,15.);
   histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.1,2.1);
   histo2DCompCombinedRatio2->DrawCopy();

      DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMForwardpPb,20,markerSizeComparison, kBlue+2, kBlue+2);
      graphRatioLowPtChargedPionsPCMForwardpPb->Draw("E1psame");
      DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMMidpPb,25,markerSizeComparison, kBlue+2, kBlue+2);
      graphRatioLowPtChargedPionsPCMMidpPb->Draw("E1psame");

      DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOSpPb,20,markerSizeComparison,  kMagenta+1, kMagenta+1);
      graphRatioLowPtChargedPionsPHOSpPb->Draw("E1psame");

      TLatex *labelRatioPi0pPb = new TLatex(0.16,0.9,collisionSystempPb.Data());
      SetStyleTLatex( labelRatioPi0pPb, 0.06,4);
      labelRatioPi0pPb->Draw();

      TLegend* legendPi0CompIndChargedPionspPb = new TLegend(0.13,0.15,0.98,0.21);
      legendPi0CompIndChargedPionspPb->SetFillColor(0);
      legendPi0CompIndChargedPionspPb->SetLineColor(0);
      legendPi0CompIndChargedPionspPb->SetNColumns(2);
      legendPi0CompIndChargedPionspPb->SetNColumns(2);
      legendPi0CompIndChargedPionspPb->SetTextSize(0.038);
      legendPi0CompIndChargedPionspPb->SetMargin(0.14);
      legendPi0CompIndChargedPionspPb->AddEntry(graphRatioLowPtChargedPionsPCMForwardpPb,"#pi^{0}/#pi^{#pm} low pt (PCM, 0.165 < y < 0.765)","p");
      legendPi0CompIndChargedPionspPb->AddEntry(graphRatioLowPtChargedPionsPHOSpPb,"#pi^{0}/#pi^{#pm} low pt (PHOS)","p");
      legendPi0CompIndChargedPionspPb->AddEntry(graphRatioLowPtChargedPionsPCMMidpPb,"#pi^{0}/#pi^{#pm} low pt (PCM, |y| < 0.4)","p");
      legendPi0CompIndChargedPionspPb->Draw();

      legendPi0CompIndChargedPionspPb->Draw();
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
   
   DrawGammaLines(0., 20.,1., 1.,0.1,kGray,2);
   
   
   canvasCompYieldpPbInd->Update();
   canvasCompYieldpPbInd->Print(Form("%s/ComparisonChargedToNeutralInd_pPb.%s",outputDir.Data(),suffix.Data()));

   
   //************************************************************************************************************************************
   //******************************  plotting just centrality dependend comparison individual measurements ******************************
   //************************************************************************************************************************************
   TCanvas * canvas6PartCompChargedIndPions = new TCanvas("canvas6PartCompChargedIndPions","",10,10,1834,1000);  // gives the page size      
   canvas6PartCompChargedIndPions->cd();
   DrawGammaCanvasSettings( canvas6PartCompChargedIndPions, 0.13, 0.0, 0.02, 0.09);
   
   TPad* pad6PartCompChargedIndPions1 = new TPad("pad6PartCompChargedIndPions1", "", 0., 0.52, 0.35, 1.,-1, -1, -2);
   DrawGammaPadSettings( pad6PartCompChargedIndPions1, 0.12, 0.0, 0.02, 0.);
   pad6PartCompChargedIndPions1->Draw();
   TPad* pad6PartCompChargedIndPions2 = new TPad("pad6PartCompChargedIndPions2", "", 0., 0., 0.35, 0.52,-1, -1, -2);
   DrawGammaPadSettings( pad6PartCompChargedIndPions2, 0.12, 0.0, 0., 0.12);
   pad6PartCompChargedIndPions2->Draw();
   
   TPad* pad6PartCompChargedIndPions3 = new TPad("pad6PartCompChargedIndPions3", "", 0.35, 0.52, 0.68, 1.,-1, -1, -2);
   DrawGammaPadSettings( pad6PartCompChargedIndPions3, 0.0, 0.0, 0.02, 0.);
   pad6PartCompChargedIndPions3->Draw();
   TPad* pad6PartCompChargedIndPions4 = new TPad("pad6PartCompChargedIndPions4", "", 0.35, 0., 0.68, 0.52,-1, -1, -2);
   DrawGammaPadSettings( pad6PartCompChargedIndPions4, 0.0, 0.0, 0., 0.12);
   pad6PartCompChargedIndPions4->Draw();

   TPad* pad6PartCompChargedIndPions5 = new TPad("pad6PartCompChargedIndPions5", "", 0.68, 0.52, 1., 1.,-1, -1, -2);
   DrawGammaPadSettings( pad6PartCompChargedIndPions5, 0.0, 0.02, 0.02, 0.);
   pad6PartCompChargedIndPions5->Draw();
   TPad* pad6PartCompChargedIndPions6 = new TPad("pad6PartCompChargedIndPions6", "", 0.68, 0., 1., 0.52,-1, -1, -2);
   DrawGammaPadSettings( pad6PartCompChargedIndPions6, 0.0, 0.02, 0., 0.12);
   pad6PartCompChargedIndPions6->Draw();

   TH2F* histo2DCompCombinedRatio;
   histo2DCompCombinedRatio = new TH2F("histo2DCompCombinedRatio","histo2DCompCombinedRatio",1000,0.3,20.,1000,0.2,4.   );
   histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.1,2.1);
   histo2DCompCombinedRatio->GetXaxis()->SetRangeUser(-0.05,20.);
   SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatio, "p_{T} (GeV/c)","#pi^{0}/#pi^{#pm}", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 

   
   pad6PartCompChargedIndPions1->cd();
   pad6PartCompChargedIndPions1->SetLogx();
   histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,15.);
   histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.1,2.1);
   histo2DCompCombinedRatio2->DrawCopy();
      
      DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMForwardpPb0020,20,markerSizeComparison, kBlue+2, kBlue+2);
      graphRatioLowPtChargedPionsPCMForwardpPb0020->Draw("E1psame");
      DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMMidpPb0020,25,markerSizeComparison, kBlue+2, kBlue+2);
      graphRatioLowPtChargedPionsPCMMidpPb0020->Draw("E1psame");

//       DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOS0005,20,markerSizeComparison, kMagenta+1, kMagenta+1);
//       graphRatioLowPtChargedPionsPHOS0005->Draw("E1psame");
       TLatex *labelPi0CompChargedPions0020 = new TLatex(0.15,0.9,collisionSystem0020.Data());
      SetStyleTLatex( labelPi0CompChargedPions0020, 0.05,4);
      labelPi0CompChargedPions0020->Draw();
      
      TLegend* legendPi0CompChargedIndPionspPb0020 = new TLegend(0.13,0.76,0.98,0.88);
      legendPi0CompChargedIndPionspPb0020->SetFillColor(0);
      legendPi0CompChargedIndPionspPb0020->SetLineColor(0);
      legendPi0CompChargedIndPionspPb0020->SetNColumns(2);
      legendPi0CompChargedIndPionspPb0020->SetTextSize(0.038);
      legendPi0CompChargedIndPionspPb0020->SetMargin(0.14);
      legendPi0CompChargedIndPionspPb0020->AddEntry(graphRatioLowPtChargedPionsPCMForwardpPb,"#pi^{0}/#pi^{#pm} low pt (PCM, 0.165 < y < 0.765)","p");
      legendPi0CompChargedIndPionspPb0020->AddEntry(graphRatioLowPtChargedPionsPHOSpPb,"#pi^{0}/#pi^{#pm} low pt (PHOS)","p");
      legendPi0CompChargedIndPionspPb0020->AddEntry(graphRatioLowPtChargedPionsPCMMidpPb,"#pi^{0}/#pi^{#pm} low pt (PCM, |y| < 0.4)","p");
      legendPi0CompChargedIndPionspPb0020->Draw();
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
   
   histo2DCompCombinedRatio2->Draw("axis,same");
   pad6PartCompChargedIndPions1->Update();
   
   pad6PartCompChargedIndPions2->cd();
   pad6PartCompChargedIndPions2->SetLogx();
   histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.1,2.1);
   histo2DCompCombinedRatio2->DrawCopy();
   
//       DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMForwardpPb6080,20,markerSizeComparison, kBlue+2, kBlue+2);
//       graphRatioLowPtChargedPionsPCMForwardpPb6080->Draw("E1psame");
//       DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMMidpPb6080,25,markerSizeComparison, kBlue+2, kBlue+2);
//       graphRatioLowPtChargedPionsPCMMidpPb6080->Draw("E1psame");
// 
// //       DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOS2040,20,markerSizeComparison, kMagenta+1, kMagenta+1);
// //       graphRatioLowPtChargedPionsPHOS2040->Draw("E1psame");
//    
   DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMForwardpPb6080,20,markerSizeComparison, kBlue+2, kBlue+2);
   graphRatioLowPtChargedPionsPCMForwardpPb6080->Draw("E1psame");
   DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMMidpPb6080,25,markerSizeComparison, kBlue+2, kBlue+2);
   graphRatioLowPtChargedPionsPCMMidpPb6080->Draw("E1psame");

      TLatex *labelPi0CompChargedPionspPb6080 = new TLatex(0.15,0.93,collisionSystem6080.Data());
      SetStyleTLatex( labelPi0CompChargedPionspPb6080, 0.05,4);
      labelPi0CompChargedPionspPb6080->Draw();

      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray); 

   histo2DCompCombinedRatio2->Draw("axis,same");
   pad6PartCompChargedIndPions2->Update();


   pad6PartCompChargedIndPions3->cd();
   pad6PartCompChargedIndPions3->SetLogx();
   histo2DCompCombinedRatio->GetXaxis()->SetRangeUser(-0.25,15.);
   //  histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.1,2.1);
   histo2DCompCombinedRatio->DrawCopy();
         
      DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMForwardpPb2040,20,markerSizeComparison, kBlue+2, kBlue+2);
      graphRatioLowPtChargedPionsPCMForwardpPb2040->Draw("E1psame");
      DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMMidpPb2040,25,markerSizeComparison, kBlue+2, kBlue+2);
      graphRatioLowPtChargedPionsPCMMidpPb2040->Draw("E1psame");

 /*     DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOS0510,20,markerSizeComparison, kMagenta+1, kMagenta+1);
      graphRatioLowPtChargedPionsPHOS0510->Draw("E1psame");
 */     
      TLatex *labelPi0CompChargedPionsPbPb2040 = new TLatex(0.03,0.9,collisionSystem2040.Data());
      SetStyleTLatex( labelPi0CompChargedPionsPbPb2040, 0.05,4);
      labelPi0CompChargedPionsPbPb2040->Draw(); 
      
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
   
   histo2DCompCombinedRatio->Draw("axis,same");
   pad6PartCompChargedIndPions3->Update();
   pad6PartCompChargedIndPions4->cd();
   pad6PartCompChargedIndPions4->SetLogx();
   histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.1,2.1);
   histo2DCompCombinedRatio->DrawCopy();

      DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMForwardpPb60100,20,markerSizeComparison, kBlue+2, kBlue+2);
      graphRatioLowPtChargedPionsPCMForwardpPb60100->Draw("E1psame");
      DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMMidpPb60100,25,markerSizeComparison, kBlue+2, kBlue+2);
      graphRatioLowPtChargedPionsPCMMidpPb60100->Draw("E1psame");

//       DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOS4060,20,markerSizeComparison,  kMagenta+1, kMagenta+1);
//       graphRatioLowPtChargedPionsPHOS4060->Draw("E1psame");
  
      TLatex *labelPi0CompChargedPionsPbPb60100 = new TLatex(0.03,0.93,collisionSystem60100.Data());
      SetStyleTLatex( labelPi0CompChargedPionsPbPb60100, 0.047,4);
      
      labelPi0CompChargedPionsPbPb60100->Draw();    
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);

   histo2DCompCombinedRatio->Draw("axis,same");
   
   pad6PartCompChargedIndPions4->Update();
   pad6PartCompChargedIndPions5->cd();
   pad6PartCompChargedIndPions5->SetLogx();
   histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.1,2.1);
   histo2DCompCombinedRatio->DrawCopy();
   
      DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMForwardpPb4060,20,markerSizeComparison, kBlue+2, kBlue+2);
      graphRatioLowPtChargedPionsPCMForwardpPb4060->Draw("E1psame");
      DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMMidpPb4060,25,markerSizeComparison, kBlue+2, kBlue+2);
      graphRatioLowPtChargedPionsPCMMidpPb4060->Draw("E1psame");

 /*     DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOS1020,20,markerSizeComparison,  kMagenta+1, kMagenta+1);
      graphRatioLowPtChargedPionsPHOS1020->Draw("E1psame");
 */  
      TLatex *labelPi0CompChargedPionsPbPb4060 = new TLatex(0.03,0.9,collisionSystem4060.Data());
      SetStyleTLatex( labelPi0CompChargedPionsPbPb4060, 0.05,4);
      labelPi0CompChargedPionsPbPb4060->Draw();    

      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
      
   histo2DCompCombinedRatio->Draw("axis,same");
   pad6PartCompChargedIndPions5->Update();
   pad6PartCompChargedIndPions6->cd();
   pad6PartCompChargedIndPions6->SetLogx();
   histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.1,2.1);
   histo2DCompCombinedRatio->DrawCopy();
   
      DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMForwardpPb,20,markerSizeComparison, kBlue+2, kBlue+2);
      graphRatioLowPtChargedPionsPCMForwardpPb->Draw("E1psame");
      DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCMMidpPb,25,markerSizeComparison, kBlue+2, kBlue+2);
      graphRatioLowPtChargedPionsPCMMidpPb->Draw("E1psame");

//       DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOSpPb,20,markerSizeComparison,  kMagenta+1, kMagenta+1);
      graphRatioLowPtChargedPionsPHOSpPb->Draw("E1psame");

      TLatex *labelPi0CompChargedPionsPbPb = new TLatex(0.04,0.93,collisionSystempPb.Data());
      SetStyleTLatex( labelPi0CompChargedPionsPbPb, 0.047,4);
      labelPi0CompChargedPionsPbPb->Draw(); 

      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);

   histo2DCompCombinedRatio->Draw("axis,same");
   
   pad6PartCompChargedIndPions6->Update();

   canvas6PartCompChargedIndPions->Update(); 
   canvas6PartCompChargedIndPions->SaveAs(Form("%s/ComparisonChargedToNeutralInd_6Parted_Paper.%s",outputDir.Data(),suffix.Data()));
   delete pad6PartCompChargedIndPions1;   
   delete pad6PartCompChargedIndPions2;   
   delete pad6PartCompChargedIndPions3;   
   delete pad6PartCompChargedIndPions4;   
   delete canvas6PartCompChargedIndPions;

}
