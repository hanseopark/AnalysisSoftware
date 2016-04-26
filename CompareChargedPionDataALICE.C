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

void CompareChargedPionDataALICE( TString suffix = "eps"){

   gROOT->Reset();   
   gROOT->SetStyle("Plain");
   
   StyleSettingsThesis();  
   SetPlotStyle();
   
   
   TString dateForOutput	 			= ReturnDateStringForOutput();
   TString outputDir 					= Form("%s/%s/ComparisonNeutralAndChargedPions",suffix.Data(),dateForOutput.Data());
   TString fileNameChargedPionPbPbOld 	= "ExternalInputPbPb/IdentifiedCharged/ChargedPionSpectraPbPb_4_Aug_2013.root";
   TString fileNameChargedPionPPOld 	= "ExternalInput/IdentifiedCharged/ChargedPionSpectraPP_4_Aug_2013.root";
   TString fileNameChargedPionPbPbNew 	= "ExternalInputPbPb/IdentifiedCharged/ChargedPionSpectraPbPb_5_Nov_2013.root";
   TString fileNameChargedPionPPNew 	= "ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_5_Nov_2013.root";
   TString fileNameCaloPhos2760GeVOld 	= "ExternalInput/PHOS/2.76TeV/LHC11a_PHOS_pi0_pp2760_noBWCorr_FDcorr_20130913.root";
   TString fileNameCaloPhos2760GeVNew 	= "ExternalInput/PHOS/2.76TeV/LHC11a_PHOS_pi0_pp2760_noBWCorr_FDcorr_20131107.root";
   TString fileNameCaloPHOSPbPbOld 		= "ExternalInputPbPb/PHOS/LHC10h_PHOS_pi0_PbPb_28102013.root";
   TString fileNameCaloPHOSPbPbNew 		= "ExternalInputPbPb/PHOS/LHC10h_PHOS_pi0_PbPb_08112013.root";      
   gSystem->Exec("mkdir -p "+outputDir);
   
   Color_t  colorComb0005           = kRed+1;
   Color_t  colorComb0510           = 807;
   Color_t  colorComb1020           = 800;
   Color_t  colorComb2040           = kGreen+2;
   Color_t  colorComb4060           = kCyan+2;
   Color_t  colorComb6080           = kBlue+1;

   Style_t  markerStyleCommmonSpectrum0005   = 20 ;
   Style_t  markerStyleCommmonSpectrum0510   = 21 ;
   Style_t  markerStyleCommmonSpectrum1020   = 29 ;
   Style_t  markerStyleCommmonSpectrum2040   = 33 ;
   Style_t  markerStyleCommmonSpectrum4060   = 20 ;
   Style_t  markerStyleCommmonSpectrum6080   = 21 ;

   Size_t   markerSizeCommonSpectrum0005  = 2.;
   Size_t   markerSizeCommonSpectrum0510  = 2.;
   Size_t   markerSizeCommonSpectrum1020  = 2.5;
   Size_t   markerSizeCommonSpectrum2040  = 2.5;
   Size_t   markerSizeCommonSpectrum4060  = 2.;
   Size_t   markerSizeCommonSpectrum6080  = 2.;
   
   Width_t  widthLinesBoxes;

   TString collisionSystemPP = "pp #sqrt{#it{s}} = 2.76 TeV";    
   TString collisionSystemCent0 = "0-5% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";      
   TString collisionSystemCent1 = "5-10% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";     
   TString collisionSystemCent2 = "10-20% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";    
   TString collisionSystemCent = "0-20% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";      
   TString collisionSystemSemiCent = "20-40% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";    
   TString collisionSystemSemiPer = "40-60% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";     
   TString collisionSystemPer = "60-80% Pb-Pb #sqrt{s_{_{NN}}} = 2.76 TeV";      

   Size_t markerSizeComparison = 1.5;
   
   
   TFile* fileChargedPionInputPbPbOld = new TFile(fileNameChargedPionPbPbOld.Data());
   TH1D* histoChargedPionHighPtOldStat0005 = (TH1D*)fileChargedPionInputPbPbOld->Get("histoChargedPionSpecHighPtStat0005");
   TH1D* histoChargedPionHighPtOldSyst0005 = (TH1D*)fileChargedPionInputPbPbOld->Get("histoChargedPionSpecHighPtSyst0005");
   TH1D* histoChargedPionHighPtOldStat0510 = (TH1D*)fileChargedPionInputPbPbOld->Get("histoChargedPionSpecHighPtStat0510");
   TH1D* histoChargedPionHighPtOldSyst0510 = (TH1D*)fileChargedPionInputPbPbOld->Get("histoChargedPionSpecHighPtSyst0510");
   TH1D* histoChargedPionHighPtOldStat1020 = (TH1D*)fileChargedPionInputPbPbOld->Get("histoChargedPionSpecHighPtStat1020");
   TH1D* histoChargedPionHighPtOldSyst1020 = (TH1D*)fileChargedPionInputPbPbOld->Get("histoChargedPionSpecHighPtSyst1020");
   TH1D* histoChargedPionHighPtOldSyst2040 = (TH1D*)fileChargedPionInputPbPbOld->Get("histoChargedPionSpecHighPtSyst2040");
   TH1D* histoChargedPionHighPtOldStat2040 = (TH1D*)fileChargedPionInputPbPbOld->Get("histoChargedPionSpecHighPtStat2040");
   TH1D* histoChargedPionHighPtOldSyst4060 = (TH1D*)fileChargedPionInputPbPbOld->Get("histoChargedPionSpecHighPtSyst4060");
   TH1D* histoChargedPionHighPtOldStat4060 = (TH1D*)fileChargedPionInputPbPbOld->Get("histoChargedPionSpecHighPtStat4060");
   TH1D* histoChargedPionHighPtOldSyst6080 = (TH1D*)fileChargedPionInputPbPbOld->Get("histoChargedPionSpecHighPtSyst6080");
   TH1D* histoChargedPionHighPtOldStat6080 = (TH1D*)fileChargedPionInputPbPbOld->Get("histoChargedPionSpecHighPtStat6080");

   TFile* fileChargedPionInputppOld = new TFile(fileNameChargedPionPPOld.Data());
   TH1D* histoChargedPionHighPtOldStatPP = (TH1D*)fileChargedPionInputppOld->Get("histoChargedPionSpecHighPtStatPP");
   TH1D* histoChargedPionHighPtOldSystPP = (TH1D*)fileChargedPionInputppOld->Get("histoChargedPionSpecHighPtSystPP");

   TFile* fileChargedPionInputPbPbNew = new TFile(fileNameChargedPionPbPbNew.Data());
   TH1D* histoChargedPionHighPtNewStat0005 = (TH1D*)fileChargedPionInputPbPbNew->Get("histoChargedPionSpecHighPtStat0005");
   TH1D* histoChargedPionHighPtNewSyst0005 = (TH1D*)fileChargedPionInputPbPbNew->Get("histoChargedPionSpecHighPtSyst0005");
   TH1D* histoChargedPionHighPtNewStat0510 = (TH1D*)fileChargedPionInputPbPbNew->Get("histoChargedPionSpecHighPtStat0510");
   TH1D* histoChargedPionHighPtNewSyst0510 = (TH1D*)fileChargedPionInputPbPbNew->Get("histoChargedPionSpecHighPtSyst0510");
   TH1D* histoChargedPionHighPtNewStat1020 = (TH1D*)fileChargedPionInputPbPbNew->Get("histoChargedPionSpecHighPtStat1020");
   TH1D* histoChargedPionHighPtNewSyst1020 = (TH1D*)fileChargedPionInputPbPbNew->Get("histoChargedPionSpecHighPtSyst1020");
   TH1D* histoChargedPionHighPtNewSyst2040 = (TH1D*)fileChargedPionInputPbPbNew->Get("histoChargedPionSpecHighPtSyst2040");
   TH1D* histoChargedPionHighPtNewStat2040 = (TH1D*)fileChargedPionInputPbPbNew->Get("histoChargedPionSpecHighPtStat2040");
   TH1D* histoChargedPionHighPtNewSyst4060 = (TH1D*)fileChargedPionInputPbPbNew->Get("histoChargedPionSpecHighPtSyst4060");
   TH1D* histoChargedPionHighPtNewStat4060 = (TH1D*)fileChargedPionInputPbPbNew->Get("histoChargedPionSpecHighPtStat4060");
   TH1D* histoChargedPionHighPtNewSyst6080 = (TH1D*)fileChargedPionInputPbPbNew->Get("histoChargedPionSpecHighPtSyst6080");
   TH1D* histoChargedPionHighPtNewStat6080 = (TH1D*)fileChargedPionInputPbPbNew->Get("histoChargedPionSpecHighPtStat6080");

   TFile* fileChargedPionInputppNew = new TFile(fileNameChargedPionPPNew.Data());
   TH1D* histoChargedPionHighPtNewStatPP = (TH1D*)fileChargedPionInputppNew->Get("histoChargedPionSpecHighPtStat2760GeV");
   TH1D* histoChargedPionHighPtNewSystPP = (TH1D*)fileChargedPionInputppNew->Get("histoChargedPionSpecHighPtSyst2760GeV");


   TH1D* ratioChargedOldDivNewStat0005 = (TH1D*)histoChargedPionHighPtOldStat0005->Clone("ratioChargedOldDivNewStat0005");
   ratioChargedOldDivNewStat0005->Divide(ratioChargedOldDivNewStat0005,histoChargedPionHighPtNewStat0005,1,1,"B");
   TH1D* ratioChargedOldDivNewStat0510 = (TH1D*)histoChargedPionHighPtOldStat0510->Clone("ratioChargedOldDivNewStat0510");
   ratioChargedOldDivNewStat0510->Divide(ratioChargedOldDivNewStat0510,histoChargedPionHighPtNewStat0510,1,1,"B");
   TH1D* ratioChargedOldDivNewStat1020 = (TH1D*)histoChargedPionHighPtOldStat1020->Clone("ratioChargedOldDivNewStat1020");
   ratioChargedOldDivNewStat1020->Divide(ratioChargedOldDivNewStat1020,histoChargedPionHighPtNewStat1020,1,1,"B");
   TH1D* ratioChargedOldDivNewStat2040 = (TH1D*)histoChargedPionHighPtOldStat2040->Clone("ratioChargedOldDivNewStat2040");
   ratioChargedOldDivNewStat2040->Divide(ratioChargedOldDivNewStat2040,histoChargedPionHighPtNewStat2040,1,1,"B");
   TH1D* ratioChargedOldDivNewStat4060 = (TH1D*)histoChargedPionHighPtOldStat4060->Clone("ratioChargedOldDivNewStat4060");
   ratioChargedOldDivNewStat4060->Divide(ratioChargedOldDivNewStat4060,histoChargedPionHighPtNewStat4060,1,1,"B");
   TH1D* ratioChargedOldDivNewStat6080 = (TH1D*)histoChargedPionHighPtOldStat6080->Clone("ratioChargedOldDivNewStat6080");
   ratioChargedOldDivNewStat6080->Divide(ratioChargedOldDivNewStat6080,histoChargedPionHighPtNewStat6080,1,1,"B");
   TH1D* ratioChargedOldDivNewStatPP2760GeV = (TH1D*)histoChargedPionHighPtOldStatPP->Clone("ratioChargedOldDivNewStatPP2760GeV");
   ratioChargedOldDivNewStatPP2760GeV->Divide(ratioChargedOldDivNewStatPP2760GeV,histoChargedPionHighPtNewStatPP,1,1,"B");
   
   
   TFile* filePhos2760GeVOld =       new TFile(fileNameCaloPhos2760GeVOld);
   TDirectory* directoryPHOSPi02760GeVOld =   (TDirectory*)filePhos2760GeVOld->Get("pp2760"); 
   TH1D* histoPi0Phos2760GeVOld =      (TH1D*)directoryPHOSPi02760GeVOld->Get("hPi02760GeVStat");
   TH1D* histoPi0PhosSys2760GeVOld =   (TH1D*)directoryPHOSPi02760GeVOld->Get("hPi02760GeVSys");
   

   TFile* filePhos2760GeVNew =       new TFile(fileNameCaloPhos2760GeVNew);
   TDirectory* directoryPHOSPi02760GeVNew =   (TDirectory*)filePhos2760GeVNew->Get("pp2760"); 
   TH1D* histoPi0Phos2760GeVNew =      (TH1D*)directoryPHOSPi02760GeVNew->Get("hPi02760GeVStat");
   TH1D* histoPi0PhosSys2760GeVNew =   (TH1D*)directoryPHOSPi02760GeVNew->Get("hPi02760GeVSys");
   
   TH1D* ratioPHOSPP2760GeV = (TH1D*)histoPi0Phos2760GeVOld->Clone("ratioPHOSPP2760GeV");
   ratioPHOSPP2760GeV->Divide(ratioPHOSPP2760GeV,histoPi0Phos2760GeVNew,1,1,"B");
   
   TH1D* relativeSystematicErrorsPHOSPP2760GeVOld = (TH1D*)histoPi0PhosSys2760GeVOld->Clone("relativeSystematicErrorsPHOSPP2760GeVOld");
   TH1D* relativeSystematicErrorsPHOSPP2760GeVNew = (TH1D*)histoPi0PhosSys2760GeVNew->Clone("relativeSystematicErrorsPHOSPP2760GeVNew");
   
   for (Int_t i = 1; i < relativeSystematicErrorsPHOSPP2760GeVOld->GetNbinsX(); i++){
      relativeSystematicErrorsPHOSPP2760GeVOld->SetBinContent(i,relativeSystematicErrorsPHOSPP2760GeVOld->GetBinError(i)/relativeSystematicErrorsPHOSPP2760GeVOld->GetBinContent(i)*100);
      relativeSystematicErrorsPHOSPP2760GeVOld->SetBinError(i,0);
      relativeSystematicErrorsPHOSPP2760GeVNew->SetBinContent(i,relativeSystematicErrorsPHOSPP2760GeVNew->GetBinError(i)/relativeSystematicErrorsPHOSPP2760GeVNew->GetBinContent(i)*100);
      relativeSystematicErrorsPHOSPP2760GeVNew->SetBinError(i,0);
   }  
   
   TFile* filePHOSPbPbOld =       new TFile(fileNameCaloPHOSPbPbOld);
   TDirectory* directoryPHOSPi0PbPb0005Old = (TDirectory*)filePHOSPbPbOld->Get("pi0_PbPb_2760_Centrality_0-5%");
   TDirectory* directoryPHOSPi0PbPb0510Old = (TDirectory*)filePHOSPbPbOld->Get("pi0_PbPb_2760_Centrality_5-10%");
   TDirectory* directoryPHOSPi0PbPb1020Old = (TDirectory*)filePHOSPbPbOld->Get("pi0_PbPb_2760_Centrality_10-20%");
   TDirectory* directoryPHOSPi0PbPb2040Old = (TDirectory*)filePHOSPbPbOld->Get("pi0_PbPb_2760_Centrality_20-40%");
   TDirectory* directoryPHOSPi0PbPb4060Old = (TDirectory*)filePHOSPbPbOld->Get("pi0_PbPb_2760_Centrality_40-60%");
   TDirectory* directoryPHOSPi0PbPb6080Old = (TDirectory*)filePHOSPbPbOld->Get("pi0_PbPb_2760_Centrality_60-80%");
   
   TH1D* histoPi0PHOSPbPb0005Old =     (TH1D*)directoryPHOSPi0PbPb0005Old->Get("hPi0_PbPb_cen0_NoBW_Stat");
   TH1D* histoPi0PHOSSysPbPb0005Old =  (TH1D*)directoryPHOSPi0PbPb0005Old->Get("hPi0_PbPb_cen0_NoBW_Syst");
   TH1D* histoPi0PHOSPbPb0510Old =     (TH1D*)directoryPHOSPi0PbPb0510Old->Get("hPi0_PbPb_cen1_NoBW_Stat");
   TH1D* histoPi0PHOSSysPbPb0510Old =  (TH1D*)directoryPHOSPi0PbPb0510Old->Get("hPi0_PbPb_cen1_NoBW_Syst");
   TH1D* histoPi0PHOSPbPb1020Old =     (TH1D*)directoryPHOSPi0PbPb1020Old->Get("hPi0_PbPb_cen2_NoBW_Stat");
   TH1D* histoPi0PHOSSysPbPb1020Old =  (TH1D*)directoryPHOSPi0PbPb1020Old->Get("hPi0_PbPb_cen2_NoBW_Syst");
   TH1D* histoPi0PHOSPbPb2040Old =     (TH1D*)directoryPHOSPi0PbPb2040Old->Get("hPi0_PbPb_cen3_NoBW_Stat");
   TH1D* histoPi0PHOSSysPbPb2040Old =  (TH1D*)directoryPHOSPi0PbPb2040Old->Get("hPi0_PbPb_cen3_NoBW_Syst");
   TH1D* histoPi0PHOSPbPb4060Old =     (TH1D*)directoryPHOSPi0PbPb4060Old->Get("hPi0_PbPb_cen4_NoBW_Stat");
   TH1D* histoPi0PHOSSysPbPb4060Old =  (TH1D*)directoryPHOSPi0PbPb4060Old->Get("hPi0_PbPb_cen4_NoBW_Syst");
   TH1D* histoPi0PHOSPbPb6080Old =     (TH1D*)directoryPHOSPi0PbPb6080Old->Get("hPi0_PbPb_cen5_NoBW_Stat");
   TH1D* histoPi0PHOSSysPbPb6080Old =  (TH1D*)directoryPHOSPi0PbPb6080Old->Get("hPi0_PbPb_cen5_NoBW_Syst");


   TFile* filePHOSPbPbNew =       new TFile(fileNameCaloPHOSPbPbNew);
   TDirectory* directoryPHOSPi0PbPb0005New = (TDirectory*)filePHOSPbPbNew->Get("pi0_PbPb_2760_Centrality_0-5%");
   TDirectory* directoryPHOSPi0PbPb0510New = (TDirectory*)filePHOSPbPbNew->Get("pi0_PbPb_2760_Centrality_5-10%");
   TDirectory* directoryPHOSPi0PbPb1020New = (TDirectory*)filePHOSPbPbNew->Get("pi0_PbPb_2760_Centrality_10-20%");
   TDirectory* directoryPHOSPi0PbPb2040New = (TDirectory*)filePHOSPbPbNew->Get("pi0_PbPb_2760_Centrality_20-40%");
   TDirectory* directoryPHOSPi0PbPb4060New = (TDirectory*)filePHOSPbPbNew->Get("pi0_PbPb_2760_Centrality_40-60%");
   TDirectory* directoryPHOSPi0PbPb6080New = (TDirectory*)filePHOSPbPbNew->Get("pi0_PbPb_2760_Centrality_60-80%");
   
   TH1D* histoPi0PHOSPbPb0005New =     (TH1D*)directoryPHOSPi0PbPb0005New->Get("hPi0_PbPb_cen0_NoBW_Stat");
   TH1D* histoPi0PHOSSysPbPb0005New =  (TH1D*)directoryPHOSPi0PbPb0005New->Get("hPi0_PbPb_cen0_NoBW_Syst");
   TH1D* histoPi0PHOSPbPb0510New =     (TH1D*)directoryPHOSPi0PbPb0510New->Get("hPi0_PbPb_cen1_NoBW_Stat");
   TH1D* histoPi0PHOSSysPbPb0510New =  (TH1D*)directoryPHOSPi0PbPb0510New->Get("hPi0_PbPb_cen1_NoBW_Syst");
   TH1D* histoPi0PHOSPbPb1020New =     (TH1D*)directoryPHOSPi0PbPb1020New->Get("hPi0_PbPb_cen2_NoBW_Stat");
   TH1D* histoPi0PHOSSysPbPb1020New =  (TH1D*)directoryPHOSPi0PbPb1020New->Get("hPi0_PbPb_cen2_NoBW_Syst");
   TH1D* histoPi0PHOSPbPb2040New =     (TH1D*)directoryPHOSPi0PbPb2040New->Get("hPi0_PbPb_cen3_NoBW_Stat");
   TH1D* histoPi0PHOSSysPbPb2040New =  (TH1D*)directoryPHOSPi0PbPb2040New->Get("hPi0_PbPb_cen3_NoBW_Syst");
   TH1D* histoPi0PHOSPbPb4060New =     (TH1D*)directoryPHOSPi0PbPb4060New->Get("hPi0_PbPb_cen4_NoBW_Stat");
   TH1D* histoPi0PHOSSysPbPb4060New =  (TH1D*)directoryPHOSPi0PbPb4060New->Get("hPi0_PbPb_cen4_NoBW_Syst");
   TH1D* histoPi0PHOSPbPb6080New =     (TH1D*)directoryPHOSPi0PbPb6080New->Get("hPi0_PbPb_cen5_NoBW_Stat");
   TH1D* histoPi0PHOSSysPbPb6080New =  (TH1D*)directoryPHOSPi0PbPb6080New->Get("hPi0_PbPb_cen5_NoBW_Syst");
   
   TH1D* ratioPi0PHOSOldDivNew0005 = (TH1D*) histoPi0PHOSPbPb0005Old->Clone("ratioPi0PHOSOldDivNew0005");
   ratioPi0PHOSOldDivNew0005->Divide(ratioPi0PHOSOldDivNew0005,histoPi0PHOSPbPb0005New,1,1,"B");
   TH1D* ratioPi0PHOSOldDivNew0510 = (TH1D*) histoPi0PHOSPbPb0510Old->Clone("ratioPi0PHOSOldDivNew0510");
   ratioPi0PHOSOldDivNew0510->Divide(ratioPi0PHOSOldDivNew0510,histoPi0PHOSPbPb0510New,1,1,"B");
   TH1D* ratioPi0PHOSOldDivNew1020 = (TH1D*) histoPi0PHOSPbPb1020Old->Clone("ratioPi0PHOSOldDivNew1020");
   ratioPi0PHOSOldDivNew1020->Divide(ratioPi0PHOSOldDivNew1020,histoPi0PHOSPbPb1020New,1,1,"B");
   TH1D* ratioPi0PHOSOldDivNew2040 = (TH1D*) histoPi0PHOSPbPb2040Old->Clone("ratioPi0PHOSOldDivNew2040");
   ratioPi0PHOSOldDivNew2040->Divide(ratioPi0PHOSOldDivNew2040,histoPi0PHOSPbPb2040New,1,1,"B");
   TH1D* ratioPi0PHOSOldDivNew4060 = (TH1D*) histoPi0PHOSPbPb4060Old->Clone("ratioPi0PHOSOldDivNew4060");
   ratioPi0PHOSOldDivNew4060->Divide(ratioPi0PHOSOldDivNew4060,histoPi0PHOSPbPb4060New,1,1,"B");
   TH1D* ratioPi0PHOSOldDivNew6080 = (TH1D*) histoPi0PHOSPbPb6080Old->Clone("ratioPi0PHOSOldDivNew6080");
   ratioPi0PHOSOldDivNew6080->Divide(ratioPi0PHOSOldDivNew6080,histoPi0PHOSPbPb6080New,1,1,"B");
   
   TH1D* relativeSystematicErrorsPHOSPbPb0005New= (TH1D*) histoPi0PHOSSysPbPb0005New->Clone("relativeSystematicErrorsPHOSPbPb0005New");
   TH1D* relativeSystematicErrorsPHOSPbPb0510New= (TH1D*) histoPi0PHOSSysPbPb0510New->Clone("relativeSystematicErrorsPHOSPbPb0510New");
   TH1D* relativeSystematicErrorsPHOSPbPb1020New= (TH1D*) histoPi0PHOSSysPbPb1020New->Clone("relativeSystematicErrorsPHOSPbPb1020New");
   TH1D* relativeSystematicErrorsPHOSPbPb2040New= (TH1D*) histoPi0PHOSSysPbPb2040New->Clone("relativeSystematicErrorsPHOSPbPb2040New");
   TH1D* relativeSystematicErrorsPHOSPbPb4060New= (TH1D*) histoPi0PHOSSysPbPb4060New->Clone("relativeSystematicErrorsPHOSPbPb4060New");
   TH1D* relativeSystematicErrorsPHOSPbPb6080New= (TH1D*) histoPi0PHOSSysPbPb6080New->Clone("relativeSystematicErrorsPHOSPbPb6080New");
   TH1D* relativeSystematicErrorsPHOSPbPb0005Old= (TH1D*) histoPi0PHOSSysPbPb0005Old->Clone("relativeSystematicErrorsPHOSPbPb0005Old");
   TH1D* relativeSystematicErrorsPHOSPbPb0510Old= (TH1D*) histoPi0PHOSSysPbPb0510Old->Clone("relativeSystematicErrorsPHOSPbPb0510Old");
   TH1D* relativeSystematicErrorsPHOSPbPb1020Old= (TH1D*) histoPi0PHOSSysPbPb1020Old->Clone("relativeSystematicErrorsPHOSPbPb1020Old");
   TH1D* relativeSystematicErrorsPHOSPbPb2040Old= (TH1D*) histoPi0PHOSSysPbPb2040Old->Clone("relativeSystematicErrorsPHOSPbPb2040Old");
   TH1D* relativeSystematicErrorsPHOSPbPb4060Old= (TH1D*) histoPi0PHOSSysPbPb4060Old->Clone("relativeSystematicErrorsPHOSPbPb4060Old");
   TH1D* relativeSystematicErrorsPHOSPbPb6080Old= (TH1D*) histoPi0PHOSSysPbPb6080Old->Clone("relativeSystematicErrorsPHOSPbPb6080Old");
   
   for (Int_t i = 3; i < relativeSystematicErrorsPHOSPbPb0005New->GetNbinsX(); i++){
      
      relativeSystematicErrorsPHOSPbPb0005New->SetBinContent(i,relativeSystematicErrorsPHOSPbPb0005New->GetBinError(i)/relativeSystematicErrorsPHOSPbPb0005New->GetBinContent(i)*100);
      relativeSystematicErrorsPHOSPbPb0005New->SetBinError(i,0);
      cout << relativeSystematicErrorsPHOSPbPb0005New->GetBinContent(i) << endl;
      relativeSystematicErrorsPHOSPbPb0510New->SetBinContent(i,relativeSystematicErrorsPHOSPbPb0510New->GetBinError(i)/relativeSystematicErrorsPHOSPbPb0510New->GetBinContent(i)*100);
      relativeSystematicErrorsPHOSPbPb0510New->SetBinError(i,0);
      relativeSystematicErrorsPHOSPbPb1020New->SetBinContent(i,relativeSystematicErrorsPHOSPbPb1020New->GetBinError(i)/relativeSystematicErrorsPHOSPbPb1020New->GetBinContent(i)*100);
      relativeSystematicErrorsPHOSPbPb1020New->SetBinError(i,0);
      relativeSystematicErrorsPHOSPbPb2040New->SetBinContent(i,relativeSystematicErrorsPHOSPbPb2040New->GetBinError(i)/relativeSystematicErrorsPHOSPbPb2040New->GetBinContent(i)*100);
      relativeSystematicErrorsPHOSPbPb2040New->SetBinError(i,0);
      relativeSystematicErrorsPHOSPbPb4060New->SetBinContent(i,relativeSystematicErrorsPHOSPbPb4060New->GetBinError(i)/relativeSystematicErrorsPHOSPbPb4060New->GetBinContent(i)*100);
      relativeSystematicErrorsPHOSPbPb4060New->SetBinError(i,0);
      relativeSystematicErrorsPHOSPbPb6080New->SetBinContent(i,relativeSystematicErrorsPHOSPbPb6080New->GetBinError(i)/relativeSystematicErrorsPHOSPbPb6080New->GetBinContent(i)*100);
      relativeSystematicErrorsPHOSPbPb6080New->SetBinError(i,0);
      relativeSystematicErrorsPHOSPbPb0005Old->SetBinContent(i,relativeSystematicErrorsPHOSPbPb0005Old->GetBinError(i)/relativeSystematicErrorsPHOSPbPb0005Old->GetBinContent(i)*100);
      relativeSystematicErrorsPHOSPbPb0005Old->SetBinError(i,0);
      relativeSystematicErrorsPHOSPbPb0510Old->SetBinContent(i,relativeSystematicErrorsPHOSPbPb0510Old->GetBinError(i)/relativeSystematicErrorsPHOSPbPb0510Old->GetBinContent(i)*100);
      relativeSystematicErrorsPHOSPbPb0510Old->SetBinError(i,0);
      relativeSystematicErrorsPHOSPbPb1020Old->SetBinContent(i,relativeSystematicErrorsPHOSPbPb1020Old->GetBinError(i)/relativeSystematicErrorsPHOSPbPb1020Old->GetBinContent(i)*100);
      relativeSystematicErrorsPHOSPbPb1020Old->SetBinError(i,0);
      relativeSystematicErrorsPHOSPbPb2040Old->SetBinContent(i,relativeSystematicErrorsPHOSPbPb2040Old->GetBinError(i)/relativeSystematicErrorsPHOSPbPb2040Old->GetBinContent(i)*100);
      relativeSystematicErrorsPHOSPbPb2040Old->SetBinError(i,0);
      relativeSystematicErrorsPHOSPbPb4060Old->SetBinContent(i,relativeSystematicErrorsPHOSPbPb4060Old->GetBinError(i)/relativeSystematicErrorsPHOSPbPb4060Old->GetBinContent(i)*100);
      relativeSystematicErrorsPHOSPbPb4060Old->SetBinError(i,0);
      relativeSystematicErrorsPHOSPbPb6080Old->SetBinContent(i,relativeSystematicErrorsPHOSPbPb6080Old->GetBinError(i)/relativeSystematicErrorsPHOSPbPb6080Old->GetBinContent(i)*100);
      relativeSystematicErrorsPHOSPbPb6080Old->SetBinError(i,0);
   }  
   
   
   // ***************************************************************************************************************
   // ************************************ Charged Pion  *********************************************************
   // ***************************************************************************************************************
   TCanvas * canvas6PartCompChargedPions = new TCanvas("canvas6PartCompChargedPions","",10,10,1834,1000);  // gives the page size      
   canvas6PartCompChargedPions->cd();
   DrawGammaCanvasSettings( canvas6PartCompChargedPions, 0.13, 0.0, 0.02, 0.09);
   
   TPad* pad6PartCompChargedPions1 = new TPad("pad6PartCompChargedPions1", "", 0., 0.52, 0.35, 1.,-1, -1, -2);
   DrawGammaPadSettings( pad6PartCompChargedPions1, 0.12, 0.0, 0.02, 0.);
   pad6PartCompChargedPions1->Draw();
   TPad* pad6PartCompChargedPions2 = new TPad("pad6PartCompChargedPions2", "", 0., 0., 0.35, 0.52,-1, -1, -2);
   DrawGammaPadSettings( pad6PartCompChargedPions2, 0.12, 0.0, 0., 0.12);
   pad6PartCompChargedPions2->Draw();
   
   TPad* pad6PartCompChargedPions3 = new TPad("pad6PartCompChargedPions3", "", 0.35, 0.52, 0.68, 1.,-1, -1, -2);
   DrawGammaPadSettings( pad6PartCompChargedPions3, 0.0, 0.0, 0.02, 0.);
   pad6PartCompChargedPions3->Draw();
   TPad* pad6PartCompChargedPions4 = new TPad("pad6PartCompChargedPions4", "", 0.35, 0., 0.68, 0.52,-1, -1, -2);
   DrawGammaPadSettings( pad6PartCompChargedPions4, 0.0, 0.0, 0., 0.12);
   pad6PartCompChargedPions4->Draw();

   TPad* pad6PartCompChargedPions5 = new TPad("pad6PartCompChargedPions5", "", 0.68, 0.52, 1., 1.,-1, -1, -2);
   DrawGammaPadSettings( pad6PartCompChargedPions5, 0.0, 0.02, 0.02, 0.);
   pad6PartCompChargedPions5->Draw();
   TPad* pad6PartCompChargedPions6 = new TPad("pad6PartCompChargedPions6", "", 0.68, 0., 1., 0.52,-1, -1, -2);
   DrawGammaPadSettings( pad6PartCompChargedPions6, 0.0, 0.02, 0., 0.12);
   pad6PartCompChargedPions6->Draw();

   TH2F * histo2DCompCombinedRatio2;
   histo2DCompCombinedRatio2 = new TH2F("histo2DCompCombinedRatio2","histo2DCompCombinedRatio2",1000,0.3,20.,1000,0.2,4.   );
   histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.82,1.15);
   histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,20.);
   histo2DCompCombinedRatio2->GetXaxis()->SetLabelOffset(-0.015);
   SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatio2, "#it{p}_{T} (GeV/#it{c})","ratio", 0.05,0.064, 0.05,0.06, 0.8,0.9, 512, 505);
   
   TH2F* histo2DCompCombinedRatio;
   histo2DCompCombinedRatio = new TH2F("histo2DCompCombinedRatio","histo2DCompCombinedRatio",1000,0.3,20.,1000,0.2,4.   );
   histo2DCompCombinedRatio->GetXaxis()->SetLabelOffset(-0.015);
   histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.82,1.15);
   histo2DCompCombinedRatio->GetXaxis()->SetRangeUser(-0.05,20.);
   SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatio, "#it{p}_{T} (GeV/#it{c})","ratio", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 

   pad6PartCompChargedPions1->cd();
   pad6PartCompChargedPions1->SetLogx();
   histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,20.);
   histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.82,1.15);
   histo2DCompCombinedRatio2->DrawCopy();
      
      DrawGammaSetMarker(ratioChargedOldDivNewStat0005,21,markerSizeComparison, kBlack , kBlack);
      ratioChargedOldDivNewStat0005->Draw("E1psame");
 
      TLatex *labelPi0CompChargedPionsPbPb0005 = new TLatex(0.15,0.9,collisionSystemCent0.Data());
      SetStyleTLatex( labelPi0CompChargedPionsPbPb0005, 0.05,4);
      labelPi0CompChargedPionsPbPb0005->Draw();
            
      TLegend* legendPi0CompChargedPionsPbPb0005 = new TLegend(0.18,0.82,0.9,0.88);
      legendPi0CompChargedPionsPbPb0005->SetFillColor(0);
      legendPi0CompChargedPionsPbPb0005->SetLineColor(0);
      legendPi0CompChargedPionsPbPb0005->SetNColumns(2);
      legendPi0CompChargedPionsPbPb0005->SetTextSize(0.045);
      legendPi0CompChargedPionsPbPb0005->AddEntry(ratioChargedOldDivNewStat0005,"#pi^{#pm} high #it{p}_{T} QM12/final, only stat","p");
      legendPi0CompChargedPionsPbPb0005->Draw();
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
   
   histo2DCompCombinedRatio2->Draw("axis,same");
   pad6PartCompChargedPions1->Update();
   pad6PartCompChargedPions2->cd();
   pad6PartCompChargedPions2->SetLogx();
   histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.82,1.15);
   histo2DCompCombinedRatio2->DrawCopy();
   
      DrawGammaSetMarker(ratioChargedOldDivNewStat2040,21,markerSizeComparison, kBlack , kBlack);
      ratioChargedOldDivNewStat2040->Draw("E1psame");
 
      TLatex *labelPi0CompChargedPionsPbPb2040 = new TLatex(0.15,0.93,collisionSystemSemiCent.Data());
      SetStyleTLatex( labelPi0CompChargedPionsPbPb2040, 0.05,4);
      labelPi0CompChargedPionsPbPb2040->Draw();
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray); 

   histo2DCompCombinedRatio2->Draw("axis,same");
   pad6PartCompChargedPions2->Update();
   pad6PartCompChargedPions3->cd();
   pad6PartCompChargedPions3->SetLogx();
   histo2DCompCombinedRatio->GetXaxis()->SetRangeUser(-0.25,20.);
   histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.82,1.15);
   histo2DCompCombinedRatio->DrawCopy();
         
     DrawGammaSetMarker(ratioChargedOldDivNewStat0510,21,markerSizeComparison, kBlack , kBlack);
      ratioChargedOldDivNewStat0510->Draw("E1psame");
  
      TLatex *labelPi0CompChargedPionsPbPb0510 = new TLatex(0.03,0.9,collisionSystemCent1.Data());
      SetStyleTLatex( labelPi0CompChargedPionsPbPb0510, 0.05,4);
      labelPi0CompChargedPionsPbPb0510->Draw(); 
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
   
   histo2DCompCombinedRatio->Draw("axis,same");
   pad6PartCompChargedPions3->Update();
   pad6PartCompChargedPions4->cd();
   pad6PartCompChargedPions4->SetLogx();
   histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.82,1.15);
   histo2DCompCombinedRatio->DrawCopy();

      DrawGammaSetMarker(ratioChargedOldDivNewStat4060,21,markerSizeComparison, kBlack , kBlack);
      ratioChargedOldDivNewStat4060->Draw("E1psame");
 
      TLatex *labelPi0CompChargedPionsPbPb4060 = new TLatex(0.03,0.93,collisionSystemSemiPer.Data());
      SetStyleTLatex( labelPi0CompChargedPionsPbPb4060, 0.047,4);
      labelPi0CompChargedPionsPbPb4060->Draw();    
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);

   histo2DCompCombinedRatio->Draw("axis,same");
   
   pad6PartCompChargedPions4->Update();
   pad6PartCompChargedPions5->cd();
   pad6PartCompChargedPions5->SetLogx();
   histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.82,1.15);
   histo2DCompCombinedRatio->DrawCopy();
   
      DrawGammaSetMarker(ratioChargedOldDivNewStat1020,21,markerSizeComparison, kBlack , kBlack);
      ratioChargedOldDivNewStat1020->Draw("E1psame");
 
      TLatex *labelPi0CompChargedPionsPbPb1020 = new TLatex(0.03,0.9,collisionSystemCent2.Data());
      SetStyleTLatex( labelPi0CompChargedPionsPbPb1020, 0.05,4);
      labelPi0CompChargedPionsPbPb1020->Draw();    
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
      
   histo2DCompCombinedRatio->Draw("axis,same");
   pad6PartCompChargedPions5->Update();
   pad6PartCompChargedPions6->cd();
   pad6PartCompChargedPions6->SetLogx();
   histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.82,1.15);
   histo2DCompCombinedRatio->DrawCopy();
   
      DrawGammaSetMarker(ratioChargedOldDivNewStat6080,21,markerSizeComparison, kBlack , kBlack);
      ratioChargedOldDivNewStat6080->Draw("E1psame");
 
      TLatex *labelPi0CompChargedPionsPbPb6080 = new TLatex(0.04,0.93,collisionSystemPer.Data());
      SetStyleTLatex( labelPi0CompChargedPionsPbPb6080, 0.047,4);
      labelPi0CompChargedPionsPbPb6080->Draw(); 
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);

   histo2DCompCombinedRatio->Draw("axis,same");
   
   pad6PartCompChargedPions6->Update();

   canvas6PartCompChargedPions->Update(); 
   canvas6PartCompChargedPions->SaveAs(Form("%s/ComparisonChargedQM12DivFinal_PbPb_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));


   pad6PartCompChargedPions1->cd();
   pad6PartCompChargedPions1->SetLogx();
   histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,20.);
   histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.82,1.15);
   histo2DCompCombinedRatio2->DrawCopy();
      
      DrawGammaSetMarker(ratioPi0PHOSOldDivNew0005,21,markerSizeComparison, kBlack , kBlack);
      ratioPi0PHOSOldDivNew0005->Draw("E1psame");
 
      labelPi0CompChargedPionsPbPb0005->Draw();
            
      TLegend* legendPi0PHOS = new TLegend(0.18,0.82,0.9,0.88);
      legendPi0PHOS->SetFillColor(0);
      legendPi0PHOS->SetLineColor(0);
      legendPi0PHOS->SetNColumns(2);
      legendPi0PHOS->SetTextSize(0.045);
      legendPi0PHOS->AddEntry(ratioPi0PHOSOldDivNew0005,"#pi^{#pm} high #it{p}_{T} QM12/final, only stat","p");
      legendPi0PHOS->Draw();
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
   
   histo2DCompCombinedRatio2->Draw("axis,same");
   pad6PartCompChargedPions1->Update();
   pad6PartCompChargedPions2->cd();
   pad6PartCompChargedPions2->SetLogx();
   histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.82,1.15);
   histo2DCompCombinedRatio2->DrawCopy();
   
      DrawGammaSetMarker(ratioPi0PHOSOldDivNew2040,21,markerSizeComparison, kBlack , kBlack);
      ratioPi0PHOSOldDivNew2040->Draw("E1psame");
 
      labelPi0CompChargedPionsPbPb2040->Draw();
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray); 

   histo2DCompCombinedRatio2->Draw("axis,same");
   pad6PartCompChargedPions2->Update();
   pad6PartCompChargedPions3->cd();
   pad6PartCompChargedPions3->SetLogx();
   histo2DCompCombinedRatio->GetXaxis()->SetRangeUser(-0.25,20.);
   histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.82,1.15);
   histo2DCompCombinedRatio->DrawCopy();
         
     DrawGammaSetMarker(ratioPi0PHOSOldDivNew0510,21,markerSizeComparison, kBlack , kBlack);
      ratioPi0PHOSOldDivNew0510->Draw("E1psame");
  
      labelPi0CompChargedPionsPbPb0510->Draw(); 
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
   
   histo2DCompCombinedRatio->Draw("axis,same");
   pad6PartCompChargedPions3->Update();
   pad6PartCompChargedPions4->cd();
   pad6PartCompChargedPions4->SetLogx();
   histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.82,1.15);
   histo2DCompCombinedRatio->DrawCopy();

      DrawGammaSetMarker(ratioPi0PHOSOldDivNew4060,21,markerSizeComparison, kBlack , kBlack);
      ratioPi0PHOSOldDivNew4060->Draw("E1psame");
 
      labelPi0CompChargedPionsPbPb4060->Draw();    
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);

   histo2DCompCombinedRatio->Draw("axis,same");
   
   pad6PartCompChargedPions4->Update();
   pad6PartCompChargedPions5->cd();
   pad6PartCompChargedPions5->SetLogx();
   histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.82,1.15);
   histo2DCompCombinedRatio->DrawCopy();
   
      DrawGammaSetMarker(ratioPi0PHOSOldDivNew1020,21,markerSizeComparison, kBlack , kBlack);
      ratioPi0PHOSOldDivNew1020->Draw("E1psame");
 
      labelPi0CompChargedPionsPbPb1020->Draw();    
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
      
   histo2DCompCombinedRatio->Draw("axis,same");
   pad6PartCompChargedPions5->Update();
   pad6PartCompChargedPions6->cd();
   pad6PartCompChargedPions6->SetLogx();
   histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.82,1.15);
   histo2DCompCombinedRatio->DrawCopy();
   
      DrawGammaSetMarker(ratioPi0PHOSOldDivNew6080,21,markerSizeComparison, kBlack , kBlack);
      ratioPi0PHOSOldDivNew6080->Draw("E1psame");
 
      labelPi0CompChargedPionsPbPb6080->Draw(); 
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);

   histo2DCompCombinedRatio->Draw("axis,same");
   
   pad6PartCompChargedPions6->Update();

   canvas6PartCompChargedPions->Update(); 
   canvas6PartCompChargedPions->SaveAs(Form("%s/ComparisonPHOS_13092013_07112013_PbPb_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));
   
   pad6PartCompChargedPions1->cd();
   pad6PartCompChargedPions1->SetLogx();
   TH2F * histo2DrelativeErrors;
   histo2DrelativeErrors = new TH2F("histo2DrelativeErrors","histo2DrelativeErrors",1000,0.3,20.,1000,0.,50.   );
   histo2DrelativeErrors->GetYaxis()->SetRangeUser(0.,40.);
   histo2DrelativeErrors->GetXaxis()->SetRangeUser(0.,20.);
   histo2DrelativeErrors->GetXaxis()->SetLabelOffset(-0.015);
   SetStyleHistoTH2ForGraphs(histo2DrelativeErrors, "#it{p}_{T} (GeV/#it{c})","rel error (%)", 0.05,0.064, 0.05,0.06, 0.8,0.9, 512, 505);
      
   histo2DrelativeErrors->DrawCopy();
      
      DrawGammaSetMarker(relativeSystematicErrorsPHOSPbPb0005New,21,0.5, kRed+2 , kRed+2);
      relativeSystematicErrorsPHOSPbPb0005New->Draw("E1psame");
      DrawGammaSetMarker(relativeSystematicErrorsPHOSPbPb0005Old,21,0.5, kBlue+2 , kBlue+2);
      relativeSystematicErrorsPHOSPbPb0005Old->Draw("E1psame");
      
      labelPi0CompChargedPionsPbPb0005->Draw();
            
      TLegend* legendRelErrorsPHOSPbPb = new TLegend(0.18,0.82,0.9,0.88);
      legendRelErrorsPHOSPbPb->SetFillColor(0);
      legendRelErrorsPHOSPbPb->SetLineColor(0);
//       legendRelErrorsPHOSPbPb->SetNColumns(2);
      legendRelErrorsPHOSPbPb->SetTextSize(0.045);
      legendRelErrorsPHOSPbPb->AddEntry(relativeSystematicErrorsPHOSPbPb0005Old,"#pi^{0} rel. syst. PHOS 28.10.2013","p");
      legendRelErrorsPHOSPbPb->AddEntry(relativeSystematicErrorsPHOSPbPb0005New,"#pi^{0} rel. syst. PHOS 08.11.2013","p");
      legendRelErrorsPHOSPbPb->Draw();
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
   
//    histo2DrelativeErrors->Draw("axis,same");
   pad6PartCompChargedPions1->Update();
   pad6PartCompChargedPions2->cd();
   pad6PartCompChargedPions2->SetLogx();
   
   histo2DrelativeErrors->DrawCopy();
   
      DrawGammaSetMarker(relativeSystematicErrorsPHOSPbPb2040New,21,0.5, kRed+2 , kRed+2);
      relativeSystematicErrorsPHOSPbPb2040New->Draw("E1psame");
      DrawGammaSetMarker(relativeSystematicErrorsPHOSPbPb2040Old,21,0.5, kBlue+2 , kBlue+2);
      relativeSystematicErrorsPHOSPbPb2040Old->Draw("E1psame");
      
      labelPi0CompChargedPionsPbPb2040->Draw();
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray); 

//    histo2DrelativeErrors->Draw("axis,same");
   pad6PartCompChargedPions2->Update();
   pad6PartCompChargedPions3->cd();
   pad6PartCompChargedPions3->SetLogx();
   histo2DrelativeErrors->DrawCopy();
         
     DrawGammaSetMarker(relativeSystematicErrorsPHOSPbPb0510New,21,0.5, kRed+2 , kRed+2);
      relativeSystematicErrorsPHOSPbPb0510New->Draw("E1psame");
      DrawGammaSetMarker(relativeSystematicErrorsPHOSPbPb0510Old,21,0.5, kBlue+2 , kBlue+2);
      relativeSystematicErrorsPHOSPbPb0510Old->Draw("E1psame");
      
      labelPi0CompChargedPionsPbPb0510->Draw(); 
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
   
//    histo2DrelativeErrors->Draw("axis,same");
   pad6PartCompChargedPions3->Update();
   pad6PartCompChargedPions4->cd();
   pad6PartCompChargedPions4->SetLogx();
   histo2DrelativeErrors->DrawCopy();

      DrawGammaSetMarker(relativeSystematicErrorsPHOSPbPb4060New,21,0.5, kRed+2 , kRed+2);
      relativeSystematicErrorsPHOSPbPb4060New->Draw("E1psame");
      DrawGammaSetMarker(relativeSystematicErrorsPHOSPbPb4060Old,21,0.5, kBlue+2 , kBlue+2);
      relativeSystematicErrorsPHOSPbPb4060Old->Draw("E1psame");
     
      labelPi0CompChargedPionsPbPb4060->Draw();    
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);

//    histo2DrelativeErrors->Draw("axis,same");
   
   pad6PartCompChargedPions4->Update();
   pad6PartCompChargedPions5->cd();
   pad6PartCompChargedPions5->SetLogx();
   histo2DrelativeErrors->DrawCopy();
   
      DrawGammaSetMarker(relativeSystematicErrorsPHOSPbPb1020New,21,0.5, kRed+2 , kRed+2);
      relativeSystematicErrorsPHOSPbPb1020New->Draw("E1psame");
      DrawGammaSetMarker(relativeSystematicErrorsPHOSPbPb1020Old,21,0.5, kBlue+2 , kBlue+2);
      relativeSystematicErrorsPHOSPbPb1020Old->Draw("E1psame");
     
 
      labelPi0CompChargedPionsPbPb1020->Draw();    
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
      
//    histo2DrelativeErrors->Draw("axis,same");
   pad6PartCompChargedPions5->Update();
   pad6PartCompChargedPions6->cd();
   pad6PartCompChargedPions6->SetLogx();
   histo2DrelativeErrors->DrawCopy();
   
      DrawGammaSetMarker(relativeSystematicErrorsPHOSPbPb6080New,21,0.5, kRed+2 , kRed+2);
      relativeSystematicErrorsPHOSPbPb6080New->Draw("E1psame");
      DrawGammaSetMarker(relativeSystematicErrorsPHOSPbPb6080Old,21,0.5, kBlue+2 , kBlue+2);
      relativeSystematicErrorsPHOSPbPb6080Old->Draw("E1psame");
     
      labelPi0CompChargedPionsPbPb6080->Draw(); 
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);

//    histo2DrelativeErrors->Draw("axis,same");
   
   pad6PartCompChargedPions6->Update();

   canvas6PartCompChargedPions->Update(); 
   canvas6PartCompChargedPions->SaveAs(Form("%s/UpdatedRelativeErrorsPHOS_PbPb_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));
   
   delete pad6PartCompChargedPions1;   
   delete pad6PartCompChargedPions2;   
   delete pad6PartCompChargedPions3;   
   delete pad6PartCompChargedPions4;   
   delete canvas6PartCompChargedPions;

   
   
   TCanvas* canvasCompYieldPPComb = new TCanvas("canvasCompYieldPPComb","",200,10,700,500);  // gives the page size
   DrawGammaCanvasSettings( canvasCompYieldPPComb,  0.12, 0.02, 0.02, 0.12);
   
   canvasCompYieldPPComb->SetLogx();
   histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,15.);
   histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.82,1.15);
   histo2DCompCombinedRatio2->DrawCopy();
   
      DrawGammaSetMarker(ratioChargedOldDivNewStatPP2760GeV,21,markerSizeComparison, kBlack , kBlack);
      ratioChargedOldDivNewStatPP2760GeV->Draw("E1psame");
      
      TLatex *labelRatioPi02760GeV = new TLatex(0.16,0.9,"pp #sqrt{#it{s}} = 2.76 TeV");
      SetStyleTLatex( labelRatioPi02760GeV, 0.06,4);
      labelRatioPi02760GeV->Draw();

      TLegend* legendPi0CompChargedPionsPP = new TLegend(0.15,0.75,0.9,0.85);
      legendPi0CompChargedPionsPP->SetFillColor(0);
      legendPi0CompChargedPionsPP->SetLineColor(0);
      legendPi0CompChargedPionsPP->SetNColumns(2);
      legendPi0CompChargedPionsPP->SetTextSize(0.045);
      legendPi0CompChargedPionsPP->SetMargin(0.12);
      legendPi0CompChargedPionsPP->AddEntry(ratioChargedOldDivNewStatPP2760GeV,"#pi^{#pm} high #it{p}_{T} QM12/final, only stat","p");
      legendPi0CompChargedPionsPP->Draw();
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
   
   DrawGammaLines(0., 20.,1., 1.,0.1,kGray,2);
   
   
   canvasCompYieldPPComb->Update();
   canvasCompYieldPPComb->Print(Form("%s/ComparisonChargedQM12DivFinal_PP2760GeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));

   canvasCompYieldPPComb->cd();
   histo2DCompCombinedRatio2->GetXaxis()->SetRangeUser(0.,15.);
   histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.89,1.11);
   histo2DCompCombinedRatio2->DrawCopy();
   
      DrawGammaSetMarker(ratioPHOSPP2760GeV,21,markerSizeComparison, kBlack , kBlack);
      ratioPHOSPP2760GeV->Draw("E1psame");
      
      labelRatioPi02760GeV->Draw();

      TLegend* legendPi0CompPHOSPP = new TLegend(0.15,0.75,0.9,0.85);
      legendPi0CompPHOSPP->SetFillColor(0);
      legendPi0CompPHOSPP->SetLineColor(0);
      legendPi0CompPHOSPP->SetNColumns(2);
      legendPi0CompPHOSPP->SetTextSize(0.045);
      legendPi0CompPHOSPP->SetMargin(0.12);
      legendPi0CompPHOSPP->AddEntry(ratioPHOSPP2760GeV,"#pi^{0} PHOS 13.09.2013/07.11.2013, only stat","p");
      legendPi0CompPHOSPP->Draw();
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
   
   DrawGammaLines(0., 20.,1., 1.,0.1,kGray,2);
   
   
   canvasCompYieldPPComb->Update();
   canvasCompYieldPPComb->Print(Form("%s/ComparisonNeutralPionPHOS_13092013_07112013_PP2760GeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));
 
   
   canvasCompYieldPPComb->cd();
      
      
//       histo2DrelativeErrors->GetYaxis()->SetRangeUser(0.82,1.15);
      histo2DrelativeErrors->DrawCopy();
      
      
      DrawGammaSetMarker(relativeSystematicErrorsPHOSPP2760GeVNew,21,0.5, kRed+2 , kRed+2);
      relativeSystematicErrorsPHOSPP2760GeVNew->Draw("E1psame");
      DrawGammaSetMarker(relativeSystematicErrorsPHOSPP2760GeVOld,21,0.5, kBlue+2 , kBlue+2);
      relativeSystematicErrorsPHOSPP2760GeVOld->Draw("E1psame");
      
      labelRatioPi02760GeV->Draw();

      TLegend* legendRelativeError = new TLegend(0.15,0.75,0.9,0.85);
      legendRelativeError->SetFillColor(0);
      legendRelativeError->SetLineColor(0);
      legendRelativeError->SetTextSize(0.045);
      legendRelativeError->SetMargin(0.12);
      legendRelativeError->AddEntry(relativeSystematicErrorsPHOSPP2760GeVOld,"#pi^{0} rel. syst. PHOS 13.09.2013","p");
      legendRelativeError->AddEntry(relativeSystematicErrorsPHOSPP2760GeVNew,"#pi^{0} rel. syst. PHOS 07.11.2013","p");
      legendRelativeError->Draw();
      DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
   
   DrawGammaLines(0., 20.,1., 1.,0.1,kGray,2);
   
   
   canvasCompYieldPPComb->Update();
   canvasCompYieldPPComb->Print(Form("%s/UpdatedRelativePHOS_PP2760GeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));
 
}