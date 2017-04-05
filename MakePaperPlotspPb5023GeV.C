// #include <Riostream.h>
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
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"

extern TRandom*    gRandom;
extern TBenchmark*    gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*      gMinuit;



void MakePaperPlotspPb5023GeV(Bool_t EPOS=kFALSE, Bool_t mT=kFALSE, Bool_t TAPS=kFALSE, Bool_t CGCPi0=kTRUE, Bool_t HIJINGPi0 = kTRUE, Bool_t DPMJetPi0 = kTRUE){    

  //input files
							        //PaperPlots_Tsallis_2017_02_09.root
  TString fileNamepPbSpectra             = "InputMakePaperPlots/PaperPlots_Tsallis_2017_03_17.root";
  TString fileNameRpPb                   = "InputMakePaperPlots/PaperPlotsRpPb_2017_03_15.root";
  TString fileNamePeaks                  = "InputMakePaperPlots/PaperPlotsPeaks_2016_11_28.root";
  TString fileNameEPOS3                  = "InputMakePaperPlots/pi0_eta_EPOS3.root";
  TString fileNameTAPS                   = "InputMakePaperPlots/TAPS_eta2pi0.root";
  TString fileNameCGC			 = "InputMakePaperPlots/TheoryGraphsCGCPi0.root";
  TString fileNameTheory                 = "ExternalInputpPb/Theory/TheoryCompilationPPb.root";
   
  
  TFile* filepPbSpectra                  = new TFile(fileNamepPbSpectra);

  if( !filepPbSpectra  ){
    cout<<"The file filepPbSpectra  does not exist "<<endl; 
    return;  
  }
  TFile* fileRpPb                          	= new TFile(fileNameRpPb);
  if( !fileRpPb ){
    cout<<"The file fileRpPb does not exist "<<endl; 
    return;  
  }   
  TFile* filePeaks                          	= new TFile(fileNamePeaks);
  if( !filePeaks ){
    cout<<"The file filePeaks does not exist "<<endl; 
    return;  
  } 
    TFile* fileEPOS3                          	= new TFile(fileNameEPOS3);
  if( !fileEPOS3 ){
    cout<<"The file fileEPOS3 does not exist "<<endl; 
    return;  
  }
  TFile* fileTAPS                          	= new TFile(fileNameTAPS);   
  if( !fileTAPS ){
    cout<<"The file fileTAPS does not exist "<<endl; 
    return;  
  }

  TFile* fileCGC			     	= new TFile(fileNameCGC);
  if( !fileCGC ){
      cout<<"The file fileCGC does not exist "<<endl;
  }
  
  TFile* fileTheoryCompilation                            = new TFile(fileNameTheory.Data());
  
  if( ! fileTheoryCompilation ) {
    
      cout<<"The file "<<fileNameTheory.Data()<<" does not exist"<<endl;
      
  }
  
  

  TF1* FitCombPi0                      = (TF1*)filepPbSpectra->Get("FitCombPi0");
  TF1* FitCombEta                      = (TF1*)filepPbSpectra->Get("FitCombEta");
  TGraphAsymmErrors* CombPi0Syst  =(TGraphAsymmErrors*)filepPbSpectra->Get("CombPi0Syst");
  TGraphAsymmErrors* CombPi0Stat =(TGraphAsymmErrors*)filepPbSpectra->Get("CombPi0Stat");
  TGraphAsymmErrors* CombEtaSyst  =(TGraphAsymmErrors*)filepPbSpectra->Get("CombEtaSyst");
  TGraphAsymmErrors* CombEtaStat =(TGraphAsymmErrors*)filepPbSpectra->Get("CombEtaStat");
  TGraphAsymmErrors*RatioTsallisCombPi0Syst  =(TGraphAsymmErrors*)filepPbSpectra->Get("RatioTsallisCombPi0Syst");
  TGraphAsymmErrors* RatioTsallisCombPi0Stat =(TGraphAsymmErrors*)filepPbSpectra->Get("RatioTsallisCombPi0Stat");
  TGraphAsymmErrors* RatioTsallisCombEtaSyst =(TGraphAsymmErrors*)filepPbSpectra->Get("RatioTsallisCombEtaSyst");
  TGraphAsymmErrors* RatioTsallisCombEtaStat =(TGraphAsymmErrors*)filepPbSpectra->Get("RatioTsallisCombEtaStat");
  TGraphAsymmErrors* RatioTsallisPCMSyst =(TGraphAsymmErrors*)filepPbSpectra->Get("RatioTsallisPCMSyst");
  TH1D* RatioTsallisPCMStat =(TH1D*)filepPbSpectra->Get("RatioTsallisPCMStat");
  TGraphAsymmErrors* RatioTsallisDalitzSyst =(TGraphAsymmErrors*)filepPbSpectra->Get("RatioTsallisDalitzSyst");
  TH1D*  RatioTsallisDalitzStat=(TH1D*)filepPbSpectra->Get("RatioTsallisDalitzStat");
  TGraphAsymmErrors* RatioTsallisEMCalSyst =(TGraphAsymmErrors*)filepPbSpectra->Get("RatioTsallisEMCalSyst");
  TH1D* RatioTsallisEMCalStat =(TH1D*)filepPbSpectra->Get("RatioTsallisEMCalStat");
  TGraphAsymmErrors*RatioTsallisPHOSSyst  =(TGraphAsymmErrors*)filepPbSpectra->Get("RatioTsallisPHOSSyst");
  TH1D* RatioTsallisPHOSStat =(TH1D*)filepPbSpectra->Get("RatioTsallisPHOSStat");
  TGraphAsymmErrors* RatioTsallisPCMEMCalSyst =(TGraphAsymmErrors*)filepPbSpectra->Get("RatioTsallisPCMEMCalSyst");
  TH1D* RatioTsallisPCMEMCalStat =(TH1D*)filepPbSpectra->Get("RatioTsallisPCMEMCalStat");
 
  TGraphAsymmErrors* RatioEtaTsallisPCMSyst =(TGraphAsymmErrors*)filepPbSpectra->Get("RatioEtaTsallisPCMSyst");
  TH1D* RatioEtaTsallisPCMStat =(TH1D*)filepPbSpectra->Get("RatioEtaTsallisPCMStat");
  TGraphAsymmErrors* RatioEtaTsallisEMCalSyst =(TGraphAsymmErrors*)filepPbSpectra->Get("RatioEtaTsallisEMCalSyst");
  TH1D* RatioEtaTsallisEMCalStat = (TH1D*)filepPbSpectra->Get("RatioEtaTsallisEMCalStat");
  TGraphAsymmErrors* RatioEtaTsallisPCMEMCalSyst =(TGraphAsymmErrors*)filepPbSpectra->Get("RatioEtaTsallisPCMEMCalSyst");
  TH1D* RatioEtaTsallisPCMEMCalStat =(TH1D*)filepPbSpectra->Get("RatioEtaTsallisPCMEMCalStat");
  
  
  
  //Theory calculations
  
   TH1F* histoDPMJetPi0                                = (TH1F*) fileTheoryCompilation->Get("histoPi0DPMJet5023TeV_Reb");
   TH1F* histoDPMJetEta                                = (TH1F*) fileTheoryCompilation->Get("histoEtaDPMJet5023TeV_Reb");
   TH1F* histoDPMJetEtaToPi0                           = (TH1F*) fileTheoryCompilation->Get("histoEtaToPi0DPMJet5023TeV");
   TH1F* histoHIJINGPi0                                = (TH1F*) fileTheoryCompilation->Get("histoPi0HIJING5023TeV_Reb");
   TH1F* histoHIJINGEta                                = (TH1F*) fileTheoryCompilation->Get("histoEtaHIJING5023TeV_Reb");
   TH1F* histoHIJINGEtaToPi0                           = (TH1F*) fileTheoryCompilation->Get("histoEtaToPi0HIJING5023TeV"); 
  ////
  
    TH1D* histoRatioPi0DPMJetToFit                      = (TH1D*) histoDPMJetPi0->Clone("histoRatioPi0DPMJetToFit"); 
    histoRatioPi0DPMJetToFit                            = CalculateHistoRatioToFit (histoRatioPi0DPMJetToFit, FitCombPi0); 
    histoRatioPi0DPMJetToFit->GetXaxis()->SetRangeUser(0.3, 20);
    TH1D* histoRatioPi0HIJINGToFit                      = (TH1D*) histoHIJINGPi0->Clone("histoRatioPi0HIJINGToFit"); 
    histoRatioPi0HIJINGToFit                            = CalculateHistoRatioToFit (histoRatioPi0HIJINGToFit, FitCombPi0); 
    histoRatioPi0HIJINGToFit->GetXaxis()->SetRangeUser(0.3, 20);
    
    TH1D* histoRatioEtaDPMJetToFit                      = (TH1D*) histoDPMJetEta->Clone("histoRatioEtaDPMJetToFit"); 
    histoRatioEtaDPMJetToFit                            = CalculateHistoRatioToFit (histoRatioEtaDPMJetToFit, FitCombEta); 
    histoRatioEtaDPMJetToFit->GetXaxis()->SetRangeUser(0.6,20);
    TH1D* histoRatioEtaHIJINGToFit                      = (TH1D*) histoHIJINGEta->Clone("histoRatioEtaHIJINGToFit"); 
    histoRatioEtaHIJINGToFit                            = CalculateHistoRatioToFit (histoRatioEtaHIJINGToFit, FitCombEta); 
    histoRatioEtaHIJINGToFit->GetXaxis()->SetRangeUser(0.6,20);
 
  
  
  TF1*  mTScaling=(TF1*)filepPbSpectra->Get("mTScaling");
  TF1*  EtaSpectrum_mT_pPb=(TF1*)filepPbSpectra->Get("EtaSpectrum_mT_pPb");
  TGraphAsymmErrors* EtaPi07TeVStat =(TGraphAsymmErrors*)filepPbSpectra->Get("EtaPi07TeVStat");
  TGraphAsymmErrors* EtaPi07TeVSyst =(TGraphAsymmErrors*)filepPbSpectra->Get("EtaPi07TeVSyst");
  TGraphAsymmErrors* EtaPi0pPbStat =(TGraphAsymmErrors*)filepPbSpectra->Get("EtaPi0pPbStat");
  TGraphAsymmErrors* EtaPi0pPbSyst =(TGraphAsymmErrors*)filepPbSpectra->Get("EtaPi0pPbSyst");
  TGraphAsymmErrors* EtaPi07TeVStat_mT =(TGraphAsymmErrors*)filepPbSpectra->Get("EtaPi0Ratio_vsmT_Stat_pp7TeV");
  TGraphAsymmErrors* EtaPi07TeVSyst_mT =(TGraphAsymmErrors*)filepPbSpectra->Get("EtaPi0Ratio_vsmT_Sys_pp7TeV");
  TGraphAsymmErrors* EtaPi0pPbStat_mT =(TGraphAsymmErrors*)filepPbSpectra->Get("EtaPi0Ratio_vsmT_Stat");
  TGraphAsymmErrors* EtaPi0pPbSyst_mT =(TGraphAsymmErrors*)filepPbSpectra->Get("EtaPi0Ratio_vsmT_Sys");

  //RpPb
  TGraphAsymmErrors*CombinedPi0RpPbSystErr=(TGraphAsymmErrors*)fileRpPb->Get("CombinedPi0RpPbSystErr");
  TGraphAsymmErrors*CombinedPi0RpPbStatErr=(TGraphAsymmErrors*)fileRpPb->Get("CombinedPi0RpPbStatErr");
  TGraph*EPS09s_KKP_NLO=(TGraph*)fileRpPb->Get("EPS09s_KKP_NLO");
   TGraph*	EPS09s_AKK_NLO=(TGraph*)fileRpPb->Get("EPS09s_AKK_NLO");
   TGraph*	EPS09s_fDSS_NLO=(TGraph*)fileRpPb->Get("EPS09s_fDSS_NLO");
   TGraphAsymmErrors*	EPS09s_fDSS_errors=(TGraphAsymmErrors*)fileRpPb->Get("EPS09s_fDSS_errors");
   TGraph*	CGC=(TGraph*)fileRpPb->Get("CGC");
   TGraphAsymmErrors*ppRefStat=(TGraphAsymmErrors*)fileRpPb->Get("ppRefStat");
   TGraphAsymmErrors*ppRefSyst=(TGraphAsymmErrors*)fileRpPb->Get("ppRefSyst");

   //Peaks


   TH1D*     PHOSPi0Width =(TH1D*)filePeaks->Get("PHOSPi0Width"); 
   TH1D*     PHOSPi0WidthMC =(TH1D*)filePeaks->Get("PHOSPi0WidthMC");        
   TH1D*     DalitzPi0Width =(TH1D*)filePeaks->Get("DalitzPi0Width");
   TH1D*     DalitzPi0WidthMC =(TH1D*)filePeaks->Get("DalitzPi0WidthMC");
   TH1D*     PCMPi0Width =(TH1D*)filePeaks->Get("PCMPi0Width");
   TH1D*     PCMPi0WidthMC =(TH1D*)filePeaks->Get("PCMPi0WidthMC");
   TH1D*     EMCalPi0Width =(TH1D*)filePeaks->Get("EMCalPi0Width");
   TH1D*     EMCalPi0WidthMC =(TH1D*)filePeaks->Get("EMCalPi0WidthMC");
   TH1D*     PCMEtaWidth =(TH1D*)filePeaks->Get("PCMEtaWidth");
   TH1D*     PCMEtaWidthMC =(TH1D*)filePeaks->Get("PCMEtaWidthMC");
   TH1D*     EMCalEtaWidth =(TH1D*)filePeaks->Get("EMCalEtaWidth");
   TH1D*     EMCalEtaWidthMC =(TH1D*)filePeaks->Get("EMCalEtaWidthMC");
   TH1D*     PHOSPi0Mass =(TH1D*)filePeaks->Get("PHOSPi0Mass");   
   TH1D*     PHOSPi0MassMC =(TH1D*)filePeaks->Get("PHOSPi0MassMC"); 
   TH1D*     DalitzPi0Mass =(TH1D*)filePeaks->Get("DalitzPi0Mass");
   TH1D*     DalitzPi0MassMC =(TH1D*)filePeaks->Get("DalitzPi0MassMC");
   TH1D*     PCMPi0Mass =(TH1D*)filePeaks->Get("PCMPi0Mass");
   TH1D*     PCMPi0MassMC =(TH1D*)filePeaks->Get("PCMPi0MassMC");
   TH1D*     EMCalPi0Mass =(TH1D*)filePeaks->Get("EMCalPi0Mass");
   TH1D*     EMCalPi0MassMC =(TH1D*)filePeaks->Get("EMCalPi0MassMC");
   TH1D*     PCMEtaMass =(TH1D*)filePeaks->Get("PCMEtaMass");
   TH1D*     PCMEtaMassMC =(TH1D*)filePeaks->Get("PCMEtaMassMC");
   TH1D*     EMCalEtaMass =(TH1D*)filePeaks->Get("EMCalEtaMass");
   TH1D*     EMCalEtaMassMC =(TH1D*)filePeaks->Get("EMCalEtaMassMC");
   //EPOS3
   TH1D*     Pi0EPOS =(TH1D*)fileEPOS3->Get("pi0_pt");Pi0EPOS->Sumw2(); 
   TH1D*     EtaEPOS =(TH1D*)fileEPOS3->Get("eta_pt");   EtaEPOS->Sumw2(); 
   Pi0EPOS=CorrectHistoToBinCenter(Pi0EPOS); 
   Pi0EPOS->Scale(1/(2*TMath::Pi()));
   EtaEPOS=CorrectHistoToBinCenter(EtaEPOS);
   EtaEPOS->Scale(1/(2*TMath::Pi()));

   //   Double_t xbins[47] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.};
   //   Double_t xbins[36] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0,8.0,9.0,10.,12.,14.,16.,18.,20.};
   Double_t xbins[36] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0,8.0,9.0,10.,11.,13.,15.,17.,20.};

   //   TH1D* Pi0EPOS2=(TH1D*)Pi0EPOS->Rebin(46,"pi0",xbins);
   TH1D* Pi0EPOS2=(TH1D*)Pi0EPOS->Rebin(35,"pi0",xbins);//35


   //     TH1D* EtaEPOS2=(TH1D*)  EtaEPOS->Rebin(46,"eta",xbins);
    TH1D* EtaEPOS2=(TH1D*)  EtaEPOS->Rebin(35,"eta",xbins);


    // for (int i=1; i<=46;i++){
    //   x=EtaEPOS2->GetBinCenter(i);
    //   if (x<1.){
    // 	Pi0EPOS2->SetBinContent(i,(Pi0EPOS2->GetBinContent(i)/10));
    // 	EtaEPOS2->SetBinContent(i,(EtaEPOS2->GetBinContent(i)/10));

    //   }
    //   else if(x>1. && x<2.){
    // 	Pi0EPOS2->SetBinContent(i,(Pi0EPOS2->GetBinContent(i)/1));
    // 	EtaEPOS2->SetBinContent(i,(EtaEPOS2->GetBinContent(i)/1));
    //   }
    //   else if(x>2. && x<10.){
    // 	Pi0EPOS2->SetBinContent(i,(Pi0EPOS2->GetBinContent(i)/5));
    // 	EtaEPOS2->SetBinContent(i,(EtaEPOS2->GetBinContent(i)/5));
    //   }
    //   else{
    // 	Pi0EPOS2->SetBinContent(i,(Pi0EPOS2->GetBinContent(i)/10));
    // 	EtaEPOS2->SetBinContent(i,(EtaEPOS2->GetBinContent(i)/10));
    //   }
    // } 
    Double_t x=0.;
    for (int i=1; i<=35;i++){
      x=EtaEPOS2->GetBinCenter(i);
      if (x<1.){
	Pi0EPOS2->SetBinContent(i,(Pi0EPOS2->GetBinContent(i)/10.));
	Pi0EPOS2->SetBinError(i,(Pi0EPOS2->GetBinError(i)/10.));
	EtaEPOS2->SetBinContent(i,(EtaEPOS2->GetBinContent(i)/10.));
	EtaEPOS2->SetBinError(i,(EtaEPOS2->GetBinError(i)/10.));
      }
      else if(x>1. && x<2.){
	Pi0EPOS2->SetBinContent(i,(Pi0EPOS2->GetBinContent(i)/1.));
	Pi0EPOS2->SetBinError(i,(Pi0EPOS2->GetBinError(i)/1.));
	EtaEPOS2->SetBinContent(i,(EtaEPOS2->GetBinContent(i)/1.));
	EtaEPOS2->SetBinError(i,(EtaEPOS2->GetBinError(i)/1.));
      }
      else if(x>2. && x<4.){
	Pi0EPOS2->SetBinContent(i,(Pi0EPOS2->GetBinContent(i)/5.));
	Pi0EPOS2->SetBinError(i,(Pi0EPOS2->GetBinError(i)/5.));
	EtaEPOS2->SetBinContent(i,(EtaEPOS2->GetBinContent(i)/5.));
	EtaEPOS2->SetBinError(i,(EtaEPOS2->GetBinError(i)/5.));
      }     else if(x>4. && x<11.){//10
	Pi0EPOS2->SetBinContent(i,(Pi0EPOS2->GetBinContent(i)/10.));
	Pi0EPOS2->SetBinError(i,(Pi0EPOS2->GetBinError(i)/10.));
	EtaEPOS2->SetBinContent(i,(EtaEPOS2->GetBinContent(i)/10.));
	EtaEPOS2->SetBinError(i,(EtaEPOS2->GetBinError(i)/10));
      }     else if(x>17.){//10
	Pi0EPOS2->SetBinContent(i,(Pi0EPOS2->GetBinContent(i)/30.));
	Pi0EPOS2->SetBinError(i,(Pi0EPOS2->GetBinError(i)/30.));
	EtaEPOS2->SetBinContent(i,(EtaEPOS2->GetBinContent(i)/30.));
	EtaEPOS2->SetBinError(i,(EtaEPOS2->GetBinError(i)/30));
      }
      else{
	Pi0EPOS2->SetBinContent(i,(Pi0EPOS2->GetBinContent(i)/20.));
	Pi0EPOS2->SetBinError(i,(Pi0EPOS2->GetBinError(i)/20.));
	EtaEPOS2->SetBinContent(i,(EtaEPOS2->GetBinContent(i)/20.));
	EtaEPOS2->SetBinError(i,(EtaEPOS2->GetBinError(i)/20.));
      }
    } 
    TGraphAsymmErrors* graphPi0EPOS                       = new TGraphAsymmErrors(Pi0EPOS2);//2
  TGraphAsymmErrors* graphEtaEPOS                       = new TGraphAsymmErrors(EtaEPOS2);
    graphPi0EPOS->RemovePoint(34);
  //  graphPi0EPOS->RemovePoint(33);
  // //   graphPi0EPOS->RemovePoint(32);
  // // graphEtaEPOS->RemovePoint(32);
    graphEtaEPOS->RemovePoint(34);
  //  graphEtaEPOS->RemovePoint(33);
 TH1D* EtaPi0EPOS=(TH1D*)EtaEPOS2->Clone("EtaPi0EPOS");
   EtaPi0EPOS->Divide(Pi0EPOS2);

  TGraphAsymmErrors* graphEtaPi0EPOS    = new TGraphAsymmErrors(EtaPi0EPOS);
    graphEtaPi0EPOS->RemovePoint(34);
  //  graphEtaPi0EPOS->RemovePoint(33);
  // // graphEtaPi0EPOS->RemovePoint(32);


 //read TAPS data

 TGraphErrors*	eta2pi0_pAu=(TGraphErrors*)fileTAPS->Get("eta2pi0_pAu");
 TGraphErrors*	eta2pi0_pBe=(TGraphErrors*)fileTAPS->Get("eta2pi0_pBe");
 
 //read CGC data
  TDirectoryFile* directoryCGC = (TDirectoryFile*)fileCGC->Get("CGC"); 

  TGraphAsymmErrors* graphAsymmErrCGCTheoryPi0y0pA5020 = (TGraphAsymmErrors*)directoryCGC->Get("graphCGCInvSecPi0pPb5020GeV");
 
 //read McGill data
 

  TGraphErrors* graphErrMcGillTheoryPion_p_hydro 	= (TGraphErrors*)fileTheoryCompilation->Get("graphPi0SpecMcGill5023TeV");
  graphErrMcGillTheoryPion_p_hydro->Print();
  TGraphErrors* graphErrMcGillTheoryEta_p_hydro  	= (TGraphErrors*)fileTheoryCompilation->Get("graphEtaSpecMcGill5023TeV");
  graphErrMcGillTheoryEta_p_hydro->Print();
  
  TGraphErrors* graphMcGillEtaToPi0                      = (TGraphErrors*) fileTheoryCompilation->Get("graphEtaToPi0McGill5023TeV");
   
  
  //read Ilkka data
 TDirectoryFile* directoryIlkka = (TDirectoryFile*)fileCGC->Get("Ilkka"); 

 TGraphAsymmErrors* graphAsymmErrIlkkapPb5020_pi0_ct14_epps16_dss14_scale_err = (TGraphAsymmErrors*)directoryIlkka->Get("graphAsymmErr_pi0_ct14_epps16_dss14_scale_sumerr");
 TGraphAsymmErrors* graphAsymmErr_pi0_ct14_epps16_dss14_sumerr                = (TGraphAsymmErrors*)directoryIlkka->Get("graphAsymmErr_pi0_ct14_epps16_dss14_sumerr");
 TGraphAsymmErrors* graphAsymmErr_pi0_ct14_epps16_dss14 		      = (TGraphAsymmErrors*)directoryIlkka->Get("graphAsymmErr_pi0_ct14_epps16_dss14");
 
 // read RpPb data
 
 TDirectoryFile* directoryRpPb = (TDirectoryFile*)fileCGC->Get("RpPb"); 

 TGraphAsymmErrors* graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14 = (TGraphAsymmErrors*)directoryRpPb->Get("graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14");
  
 
 ////////////////////////////////////////////////////////////////////////////////////////////////////
  gROOT->Reset();
  gROOT->SetStyle("Plain");
    
  StyleSettingsThesis();
  SetPlotStyle();

	Double_t 	mesonMassExpectPi0 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
	cout << mesonMassExpectPi0 << endl;
	Double_t 	mesonMassExpectEta = TDatabasePDG::Instance()->GetParticle(221)->Mass();
	cout << mesonMassExpectEta << endl;

  // general settings
  TString dateForOutput                               = ReturnDateStringForOutput();
  TString outputDir=Form("PaperPlots/%s",dateForOutput.Data());
  gSystem->Exec("mkdir -p "+outputDir);
  TString suffix="pdf";
  Double_t CMarginL=0.11;
  Double_t CMarginT=0.02;
  Double_t CMarginR=0.02;
  Double_t CMarginB=0.12;

  Int_t CDimX=500;
  Int_t CDimY=400;

  Double_t TitleOffsetX=1.;
  Double_t TitleOffsetY=1.;
  Double_t TitleSizeX=0.05;
  Double_t TitleSizeY=0.05;

  Double_t LabelOffsetX=1;
  Double_t LabelOffsetY=1.;
  Double_t LabelSizeX=0.05;
  Double_t LabelSizeY=0.05;

  Int_t NDivX=512;
  Int_t NDivY=508;

  Int_t Font=42;
  Double_t TextSize=0.045;
  Double_t TextSizeRpA=0.04;

  Color_t colorPCM   		= kBlack;
  Color_t colorPCMMC   		= kGray+1;
  Color_t colorDalitz   	= kViolet; //kBlue+1;kViolet-4 
  Color_t colorDalitzMC   	= kViolet-4; //kBlue-6;kViolet-4
  Color_t colorPHOS   		= kRed+1;
  Color_t colorPHOSMC  		= kRed-7;
  Color_t colorEMCal   		= kGreen+2;
  Color_t colorEMCalMC   	= kGreen-6 ;
  Color_t colorPCMEMCal   	= kBlue+1;
  Color_t colorPCMEMCalMC 	= kBlue-6;
  Color_t colorCombYieldPi0   	= kBlue+2;
  Color_t colorCombYieldEta 	= kGreen+3;
  Color_t colorEtaPi0RatiopPb   = kBlack;
  Color_t colorEtaPi0Ratiopp   	= kGray+2;
  Color_t colormTScaling  	= kRed+2;
  Color_t RpPb   		= kBlue+2;
  Color_t ppRef   		= kGray+2;
  Color_t  colorDPMJet          = kViolet+2;
  Color_t  colorDPMJetPi0       = kViolet+2;
  Color_t  colorDPMJetEta       = kViolet+2; //-2
  
  Color_t  colorHIJING          = kGreen-2;
  Color_t  colorHIJINGPi0       = kGreen-2;
  Color_t  colorHIJINGEta       = kGreen-2; //+2
  
  
  Color_t  colorMcGill	        = kPink + 2;
  Color_t  colorEPOS            = kAzure+7;
  Color_t  colorEPOSPi0         = kAzure+7;
  Color_t  colorEPOSEta		= kAzure+7;
  
  Color_t  colorMcGillPi0	= kPink + 2;
  Color_t  colorMcGillEta	= kPink + 2; //-6
  
  Color_t  colorIlkka           = 41;
  Color_t  colorIlkkaPi0	= kYellow - 7;//41;
  Color_t  colorEtaPi07TeV	= kBlue - 3;
  
  Style_t  styleLineDPMJet                    = 7;
  Style_t  styleLineHIJING                    = 8;
  Style_t  styleLineCGC			      = 1;
  Style_t  styleLineEPOS3		      = 1;
  Style_t  styleLineMcGill		      = 1;
  Style_t  styleLineIlkka		      = 1;
  
  
  Style_t markerStylePi0	= 20;
  Style_t markerStyleEta	= 21;
  Style_t markerStylePCM	= 20;
  Style_t markerStylePHOS	= 21;
  Style_t markerStyleEMCal 	= 33;
  Style_t markerStyleDalitz	= 29;
  Style_t markerStylePCMEMCal   = 34;
  Style_t markerStylePCMMC	= 24;
  Style_t markerStylePHOSMC	= 25;
  Style_t markerStyleEMCalMC 	= 27;
  Style_t markerStyleDalitzMC	= 30;
  Style_t markerStylePCMEMCalMC = 28;
  Size_t  markerSizePCM=1.;	
  Size_t  markerSizeDalitz=1.2;	
  Size_t  markerSizePHOS=1.;	
  Size_t  markerSizeEMCal=1.2;
  Size_t  markerSizePCMEMCal=1.2;

  Double_t LabelOffsetLog=-0.015;



  //Inv Yields
  TCanvas* c1 = new TCanvas("c1","",200,10,1000,1200);  // gives the page size
  DrawGammaCanvasSettings( c1, 0.3, 0.02, 0.02, 0.16);
  TPad* padHistosEtaandPi0 = new TPad("padHistosEtaandPi0", "", 0., 0.25, 1., 1.,-1, -1, -2);
  DrawGammaPadSettings( padHistosEtaandPi0, 0.18, 0.02, 0.02, 0.);
  padHistosEtaandPi0->Draw();

  TPad* padRatiosEtaandPi0 = new TPad("padRatiosEtaandPi0", "", 0., 0., 1., 0.25,-1, -1, -2);
  DrawGammaPadSettings( padRatiosEtaandPi0, 0.18, 0.02, 0., 0.3);
  padRatiosEtaandPi0->Draw();

  padHistosEtaandPi0->cd(); 
  padHistosEtaandPi0->SetLogx();
  padHistosEtaandPi0->SetLogy();
  TH2F * hist1;
  hist1 = new TH2F("hist1","hist1",1000,0.3,30.,1000,1.2e-9,50 );
  SetStyleHistoTH2ForGraphs(hist1, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2} ",0.04,0.05, 0.045,0.045, 0.7,1.8, 512, 510);
  hist1->DrawCopy(); 
  FitCombPi0->SetLineWidth(2);  
  if (!EPOS)  FitCombPi0->Draw("same"); 
  FitCombEta->SetLineWidth(2);  
   if (!EPOS) FitCombEta->Draw("same");
  // DrawGammaSetMarker(Pi0EPOS,8,1,kRed+2,kRed+2 );
  // DrawGammaSetMarker(EtaEPOS,8,1,kCyan,kCyan ); 



  TGraphAsymmErrors*ppRefStat_scale=(TGraphAsymmErrors*)ScaleGraph(ppRefStat,(0.0983e-9));
  //  TGraphAsymmErrors*ppRefSyst_scale=(TGraphAsymmErrors*)ScaleGraph(ppRefSyst,( 0.0983e-9));
  Double_t* EXh = ppRefStat_scale->GetEXhigh();
  Double_t* EXl = ppRefStat_scale->GetEXlow();
  Double_t* EYh = ppRefStat_scale->GetEYhigh();
  Double_t* EYl = ppRefStat_scale->GetEYlow();
  for(Int_t i = 0; i<ppRefStat_scale->GetN(); i++){
    EXh[i]=0.;
    EXl[i]=0.;
    EYh[i]=0.;
    EYl[i]=0.;
  }
  ppRefStat_scale->Print();
       
  ppRefStat_scale->SetMarkerStyle(1);
  ppRefStat_scale->SetLineStyle(7);
  ppRefStat_scale->SetLineWidth(3);
  ppRefStat_scale->SetLineColor(kGray+1);
 
  if (EPOS){
    graphPi0EPOS->SetLineWidth(2);
    graphEtaEPOS->SetLineWidth(2);
    graphPi0EPOS->SetLineColor(kAzure+7);
    graphPi0EPOS->SetFillColor(kAzure+7);
    graphEtaEPOS->SetLineColor(kGreen-3);
    graphEtaEPOS->SetFillColor(kGreen-3);
    graphPi0EPOS->Draw("C3same"); //C
    graphEtaEPOS->Draw("C3same");
  }
  
  if( CGCPi0 ) {
    
    graphAsymmErrCGCTheoryPi0y0pA5020->SetLineWidth(2);
    graphAsymmErrCGCTheoryPi0y0pA5020->SetLineColor(kOrange+1);
    graphAsymmErrCGCTheoryPi0y0pA5020->SetFillColor(kOrange+1);
    graphAsymmErrCGCTheoryPi0y0pA5020->Draw("C3same");
    
   
    
    
    
  }

  DrawGammaSetMarkerTGraphAsym(CombPi0Syst,markerStylePi0,1,colorCombYieldPi0 ,colorCombYieldPi0, 1, kTRUE);
  CombPi0Syst->Draw("E2,same"); 
  DrawGammaSetMarkerTGraphAsym(CombPi0Stat,markerStylePi0,1,colorCombYieldPi0, colorCombYieldPi0 ); 
  CombPi0Stat->Draw("pz,same"); 
  DrawGammaSetMarkerTGraphAsym(CombEtaSyst,markerStyleEta,1,colorCombYieldEta ,colorCombYieldEta, 1, kTRUE);
  CombEtaSyst->Draw("E2,same");
  DrawGammaSetMarkerTGraphAsym(CombEtaStat,markerStyleEta,1,colorCombYieldEta, colorCombYieldEta ); 
  CombEtaStat->Draw("pz,same"); 

  if (mT){
    EtaSpectrum_mT_pPb->SetLineColor(kGreen-3);
    EtaSpectrum_mT_pPb->SetLineStyle(2);
    EtaSpectrum_mT_pPb->Draw("same"); 
      }
  TLegend* leg1;
  if(EPOS ) leg1= new TLegend(0.25,0.08,0.6,0.38); 
  else if(mT ) leg1= new TLegend(0.25,0.12,0.6,0.38);
  else  leg1 = new TLegend(0.25,0.15,0.6,0.38);
  leg1->SetFillColor(0);
  leg1->SetLineColor(0);
  leg1->SetTextFont(Font);
  leg1->SetTextSize(TextSize);
  leg1->AddEntry(CombPi0Syst,"NSD #pi^{0}","pef");
  leg1->AddEntry(CombEtaSyst,"NSD #eta","pef");
  if (!EPOS)  leg1->AddEntry(FitCombPi0,"Tsallis Fit","l");
  if (EPOS) leg1->AddEntry(graphPi0EPOS,"#pi^{0} EPOS3","l");
  if (EPOS) leg1->AddEntry(graphEtaEPOS,"#eta EPOS3","l");
  if (mT) leg1->AddEntry(EtaSpectrum_mT_pPb,"#eta from #it{m}_{T} scaled #pi^{0}","l");
  if (CGCPi0) leg1->AddEntry(graphAsymmErrCGCTheoryPi0y0pA5020,"#pi^{0} CGC","l");
//  leg1->AddEntry(ppRefStat_scale,"pp Reference x #it{T}_{pA}","l");
  //leg1->AddEntry(ppRefStat_scale,"pp Reference x 1/#sigma_{NN}","l");
  leg1->Draw("same");

 	
  TLatex * lt1 = new TLatex(2.,2.1,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV") ;
  lt1->SetTextColor(kBlack) ;
  lt1->SetTextSize(TextSize) ;
  lt1->SetTextFont(Font) ;
  lt1->Draw() ;

  padRatiosEtaandPi0->cd();
    
  padRatiosEtaandPi0->SetLogx();

  //EPOS
  TGraphAsymmErrors* RatioEPOSFitPi0=(TGraphAsymmErrors*)graphPi0EPOS->Clone("RatioEPOSFitPi0");
  TGraphAsymmErrors* RatioEPOSFitEta=(TGraphAsymmErrors*)graphEtaEPOS->Clone("RatioEPOSFitEta");
  RatioEPOSFitPi0 =(TGraphAsymmErrors*)CalculateGraphErrRatioToFit(RatioEPOSFitPi0,FitCombPi0);
  RatioEPOSFitEta =(TGraphAsymmErrors*)CalculateGraphErrRatioToFit(RatioEPOSFitEta,FitCombEta);
  
  //CGC
  
  TGraphAsymmErrors* RatioCGCFitPi0=(TGraphAsymmErrors*)graphAsymmErrCGCTheoryPi0y0pA5020->Clone("RatioCGCFitPi0");
  RatioCGCFitPi0=(TGraphAsymmErrors*)CalculateGraphErrRatioToFit(RatioCGCFitPi0,FitCombPi0);
  
  //McGill
  
  TGraphErrors* RatioMcGillFitPi0=(TGraphErrors*)graphErrMcGillTheoryPion_p_hydro->Clone("RatioMcGillFitPi0");
  RatioMcGillFitPi0=(TGraphErrors*)CalculateGraphErrRatioToFit(RatioMcGillFitPi0,FitCombPi0);
  
  TGraphErrors* RatioMcGillFitEta=(TGraphErrors*)graphErrMcGillTheoryEta_p_hydro->Clone("RatioMcGillFitEta");
  RatioMcGillFitEta=(TGraphErrors*)CalculateGraphErrRatioToFit(RatioMcGillFitEta,FitCombEta);
  
  
  TGraphAsymmErrors* RatioIlkkaFitPi0 = (TGraphAsymmErrors*)graphAsymmErr_pi0_ct14_epps16_dss14_sumerr->Clone("RatioIlkkaFitPi0");
  RatioIlkkaFitPi0=(TGraphAsymmErrors*)CalculateGraphErrRatioToFit(RatioIlkkaFitPi0,FitCombPi0);
  
  TGraphAsymmErrors* RatioIlkkaFitPi0NoErr = (TGraphAsymmErrors*)graphAsymmErr_pi0_ct14_epps16_dss14->Clone("RatioIlkkaFitPi0NoErr");
  RatioIlkkaFitPi0NoErr=(TGraphAsymmErrors*)CalculateGraphErrRatioToFit(RatioIlkkaFitPi0NoErr,FitCombPi0);
  
  TGraphAsymmErrors* RatioIlkkaFitPi0scaleerr=(TGraphAsymmErrors*)graphAsymmErrIlkkapPb5020_pi0_ct14_epps16_dss14_scale_err->Clone("RatioIlkkaFitPi0scaleerr");
  RatioIlkkaFitPi0scaleerr=(TGraphAsymmErrors*)CalculateGraphErrRatioToFit(RatioIlkkaFitPi0scaleerr,FitCombPi0);
 
  
 
  TH2F * hist1a;
  hist1a = new TH2F("hist1a","histo2DRatioAllppreferencesEtaandPi0",1000,.3,30.,1000,0.31,2.19);//1.89
  SetStyleHistoTH2ForGraphs(hist1a, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.13,0.13, 0.13,0.13, 1.,0.6, 502, 505); 
  hist1a->GetYaxis()->SetLabelOffset(0.005);
  hist1a->GetXaxis()->SetLabelOffset(LabelOffsetLog+0.05);
  hist1a->GetXaxis()->SetTickLength(0.07);
   hist1a->DrawCopy();
   TLine *line1=new TLine(0.,1.,30.,1.);
   line1->SetLineColor(kGray+1);
   line1->Draw("same");

   if (EPOS){
     RatioEPOSFitPi0->SetLineWidth(2);
     //   RatioEPOSFitPi0->SetLineStyle(2); 
     RatioEPOSFitPi0->SetLineColor(kAzure+7);
     RatioEPOSFitPi0->SetFillColor(kAzure+7);
        RatioEPOSFitEta->SetLineWidth(2);
	//	  RatioEPOSFitEta->SetLineStyle(2);
     RatioEPOSFitEta->SetLineColor(kGreen-3);
     RatioEPOSFitEta->SetFillColor(kGreen-3);
     
     RatioEPOSFitEta->SetLineColor(kAzure+7);
     RatioEPOSFitEta->SetFillColor(kAzure+7);
     
     RatioEPOSFitPi0->Draw("C3same");
     RatioEPOSFitEta->Draw("C3same");
     cout<<"bla"<<endl;
     RatioEPOSFitEta->Print();
  }
  
  if( CGCPi0){
      RatioCGCFitPi0->SetLineWidth(2);
      RatioCGCFitPi0->SetLineColor(kOrange+1);
      RatioCGCFitPi0->SetFillColor(kOrange+1);
      RatioCGCFitPi0->Draw("C3same");
  }

   
   //  if (!EPOS){ 
   DrawGammaSetMarkerTGraphAsym(RatioTsallisCombPi0Syst,markerStylePi0,1,colorCombYieldPi0 ,colorCombYieldPi0 , 1, kTRUE);  
   RatioTsallisCombPi0Syst->Draw("E2,same");
   DrawGammaSetMarkerTGraphAsym(RatioTsallisCombPi0Stat,markerStylePi0,1,colorCombYieldPi0, colorCombYieldPi0 );  
   RatioTsallisCombPi0Stat->Draw("Ez,p,same");  

   DrawGammaSetMarkerTGraphAsym(RatioTsallisCombEtaSyst,markerStyleEta,1,colorCombYieldEta ,colorCombYieldEta, 1, kTRUE);  
   RatioTsallisCombEtaSyst->Draw("E2,same");
   DrawGammaSetMarkerTGraphAsym(RatioTsallisCombEtaStat,markerStyleEta,1,colorCombYieldEta ,colorCombYieldEta );  
   RatioTsallisCombEtaStat->Draw("Ez,p,same");
   // }

 // if (mT){ TEMP
 //    Ratio_mT_Fit->SetLineColor(kGreen-3);
 //    Ratio_mT_Fit->SetLineStyle(2);
 //    Ratio_mT_Fit->Draw("same");
 //  }
  c1->Update();
  if(EPOS) c1->Print(Form("%s/MesonYields_EPOS.%s",outputDir.Data(),suffix.Data()));
  else if(mT) c1->Print(Form("%s/MesonYields_mT.%s",outputDir.Data(),suffix.Data()));
  else c1->Print(Form("%s/MesonYields.%s",outputDir.Data(),suffix.Data()));
  //Inv Yields Ratio Pi0
 

 TCanvas* c2 = new TCanvas("c2","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c2, CMarginL, CMarginR,CMarginT ,CMarginB);
  c2->SetLogx();
  TH2F * hist2;
  hist2 = new TH2F("hist2","hist2",1000,0.27,25.,1000,0.41,2.6  );
  SetStyleHistoTH2ForGraphs(hist2, "#it{p}_{T} (GeV/#it{c})","Data/Fit",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY,TitleOffsetX,TitleOffsetY, 512, 508); 
  hist2->GetXaxis()->SetLabelOffset(LabelOffsetLog);
  hist2->DrawCopy(); 

  
  DrawGammaLines(0.0,25.,1.,1.,2.0,kRed+2,2);
  
  

   DrawGammaSetMarkerTGraphAsym(RatioTsallisPHOSSyst,markerStylePHOS,markerSizePHOS, colorPHOS, colorPHOS, 1, kTRUE);  
   RatioTsallisPHOSSyst->Draw("E2same");
   DrawGammaSetMarker(RatioTsallisPHOSStat, markerStylePHOS,markerSizePHOS, colorPHOS, colorPHOS );
   RatioTsallisPHOSStat->Draw("Ez,same") ;
 
   DrawGammaSetMarkerTGraphAsym(RatioTsallisEMCalSyst,markerStyleEMCal,markerSizeEMCal, colorEMCal, colorEMCal, 1, kTRUE);  
   RatioTsallisEMCalSyst->Draw("E2same");
   DrawGammaSetMarker(RatioTsallisEMCalStat,markerStyleEMCal,markerSizeEMCal, colorEMCal, colorEMCal );
   RatioTsallisEMCalStat->Draw("Ez,same") ; 

   DrawGammaSetMarkerTGraphAsym(RatioTsallisPCMSyst,markerStylePCM,markerSizePCM, colorPCM, colorPCM, 1, kTRUE);  
   RatioTsallisPCMSyst->Draw("E2same");
   DrawGammaSetMarker(RatioTsallisPCMStat,markerStylePCM,markerSizePCM, colorPCM, colorPCM );
   RatioTsallisPCMStat->Draw("Ez,same") ;
 
   DrawGammaSetMarkerTGraphAsym(RatioTsallisDalitzSyst,markerStyleDalitz,markerSizeDalitz, colorDalitz, colorDalitz, 1, kTRUE);  
   RatioTsallisDalitzSyst->Draw("E2same");
   DrawGammaSetMarker(RatioTsallisDalitzStat,markerStyleDalitz,markerSizeDalitz, colorDalitz, colorDalitz );
   RatioTsallisDalitzStat->Draw("Ez,same") ; 
  
  DrawGammaSetMarkerTGraphAsym(RatioTsallisPCMEMCalSyst,markerStylePCMEMCal,markerSizePCMEMCal,colorPCMEMCal,colorPCMEMCal,1,kTRUE);
  RatioTsallisPCMEMCalSyst->Draw("E2same");
  DrawGammaSetMarker(RatioTsallisPCMEMCalStat,markerStylePCMEMCal,markerSizePCMEMCal, colorPCMEMCal, colorPCMEMCal );
  RatioTsallisPCMEMCalStat->Draw("Ez,same") ; 
  
  
  
 	
  TLatex * lt2 = new TLatex(2.8,2.4,"ALICE") ;
  lt2->SetTextColor(kBlack) ;
  lt2->SetTextSize(TextSize) ;
  lt2->SetTextFont(Font) ;
  lt2->DrawLatex(2.8,2.25,"p-Pb, NSD, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  lt2->DrawLatex(2.8,2.08,"#pi^{0}");
 
  
  lt2->Draw() ;
  
  //TLatex * lt2 = new TLatex(2.8,2.4,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV") ;
  //lt2->SetTextColor(kBlack) ;
  //lt2->SetTextSize(TextSize) ;
  //lt2->SetTextFont(Font) ;
  //lt2->DrawLatex(2.8,2.25,"ALICE, NSD #pi^{0}");
  
  //lt2->Draw() ;

  TLegend* leg2 = new TLegend(0.18,0.65,0.5,0.95);
  leg2->SetFillColor(0);
  leg2->SetLineColor(0);
  leg2->SetTextFont(Font);
  leg2->SetTextSize(TextSize);
  leg2->AddEntry(RatioTsallisPHOSStat,Form("PHOS"),"pf");
  leg2->AddEntry(RatioTsallisEMCalStat,Form("EMC"),"pf");
  leg2->AddEntry(RatioTsallisPCMStat,Form("PCM"),"pf");
  leg2->AddEntry(RatioTsallisPCMEMCalStat,Form("PCM-EMC"),"pf");
  leg2->AddEntry(RatioTsallisDalitzStat,Form("PCM-#gamma*#gamma"),"pf");
  leg2->Draw("same");





  c2->Update();
  c2->Print(Form("%s/RatioCombFitIndividualPi0.%s",outputDir.Data(),suffix.Data()));



  //Inv Yields Ratio Eta
 


 TCanvas* c3 = new TCanvas("c3","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c3, CMarginL, CMarginR,CMarginT ,CMarginB);
  c3->SetLogx();
  TH2F * hist3;
  hist3 = new TH2F("hist3","hist3",1000,0.30,25.,1000,0.41,2.6  );
  SetStyleHistoTH2ForGraphs(hist3, "#it{p}_{T} (GeV/#it{c})","Data/Fit",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY, TitleOffsetX,TitleOffsetY, 512, 508); 
  hist3->GetXaxis()->SetLabelOffset(LabelOffsetLog);
  hist3->DrawCopy(); 

   DrawGammaLines(0.0,25.,1.,1.,2.0,kRed+2,2);
  
 

 
   DrawGammaSetMarkerTGraphAsym(RatioEtaTsallisEMCalSyst,markerStyleEMCal,markerSizeEMCal, colorEMCal, colorEMCal, 1, kTRUE);  
   DrawGammaSetMarker(RatioEtaTsallisEMCalStat,markerStyleEMCal,markerSizeEMCal, colorEMCal, colorEMCal );
   RatioEtaTsallisEMCalSyst->Draw("E2same");
   RatioEtaTsallisEMCalStat->Draw("Ez,same") ;

   DrawGammaSetMarkerTGraphAsym(RatioEtaTsallisPCMSyst,markerStylePCM,markerSizePCM, colorPCM, colorPCM, 1, kTRUE);  
   DrawGammaSetMarker(RatioEtaTsallisPCMStat,markerStylePCM,markerSizePCM, colorPCM, colorPCM ); 
   RatioEtaTsallisPCMSyst->Draw("E2same");
   RatioEtaTsallisPCMStat->Draw("Ez,same") ;	
   
   DrawGammaSetMarkerTGraphAsym(RatioEtaTsallisPCMEMCalSyst,markerStylePCMEMCal,markerSizePCMEMCal, colorPCMEMCal, colorPCMEMCal, 1, kTRUE);  
   DrawGammaSetMarker(RatioEtaTsallisPCMEMCalStat,markerStylePCMEMCal,markerSizePCMEMCal, colorPCMEMCal, colorPCMEMCal); 
   RatioEtaTsallisPCMEMCalSyst->Draw("E2same");
   RatioEtaTsallisPCMEMCalStat->Draw("Ez,same") ;


  TLegend* leg3 = new TLegend(0.18,0.8,0.5,0.95);
  leg3->SetFillColor(0);
  leg3->SetLineColor(0);
  leg3->SetTextFont(Font);
  leg3->SetTextSize(TextSize);
  leg3->AddEntry(RatioEtaTsallisEMCalStat,Form("EMC"),"pf");
  leg3->AddEntry(RatioEtaTsallisPCMStat,Form("PCM"),"pf");
  leg3->AddEntry(RatioEtaTsallisPCMEMCalStat,Form("PCM-EMC"),"pf");
  leg3->Draw("same");
  
  TLatex * lt3 = new TLatex(2.9,2.4,"ALICE");
  lt3->SetTextColor(kBlack) ;
  lt3->SetTextSize(TextSize) ;
  lt3->SetTextFont(Font) ;
  lt3->DrawLatex(2.9,2.25,"p-Pb, NSD, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  lt3->DrawLatex(2.9,2.08,"#eta");
  
  lt3->Draw() ;

  //  TLatex * lt3 = new TLatex(2.8,2.1,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV") ;
//   TLatex * lt3 = new TLatex(2.8,2.1,"") ;
//   lt3->SetTextColor(kBlack) ;
//   lt3->SetTextSize(TextSize) ;
//   lt3->SetTextFont(Font) ;
//    lt3->DrawLatex(2.8,1.95,"NSD #eta");
//   
//   lt3->Draw() ;



  c3->Update();
  c3->Print(Form("%s/RatioCombFitIndividualEta.%s",outputDir.Data(),suffix.Data()));





  //Eta/Pi0 Ratio

  

  TCanvas* c4 = new TCanvas("c4","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c4, CMarginL, CMarginR,CMarginT ,CMarginB);

  TH2F * hist4;
  if (TAPS)  hist4 = new TH2F("hist4","hist4",1000,-0.2,16.,1000,-0.02,.99   );
  else  hist4 = new TH2F("hist4","hist4",1000,0.3,16.,1000,0.0,.99   );

  SetStyleHistoTH2ForGraphs(hist4, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0} ",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY, TitleOffsetX,TitleOffsetY, 512, 508);
  hist4->DrawCopy(); 
  if (EPOS){
    graphEtaPi0EPOS->SetLineColor(kGreen+1);
    graphEtaPi0EPOS->SetFillColor(kGreen+1);
    //  graphEtaPi0EPOS->SetLineStyle(2);
    graphEtaPi0EPOS->SetLineWidth(2);
    graphEtaPi0EPOS->Draw("C3same");
  }      
  
 // if (!EPOS || (EPOS && TAPS)) mTScaling->Draw("same"); TEMP
 // graphEtaPi0Ratio_mTScaled_pPb->SetMarkerColor(kRed+2);
 // graphEtaPi0Ratio_mTScaled_pPb->SetLineColor(1);
 //  graphEtaPi0Ratio_mTScaled_pPb->SetFillColor(kRed+2);

   if (TAPS){
     DrawGammaSetMarkerTGraphErr(eta2pi0_pAu, 27, 1.4,kGreen, kGreen);
     DrawGammaSetMarkerTGraphErr(eta2pi0_pBe, 27, 1.4,kRed, kRed);
     eta2pi0_pAu->Draw("same,zp");
     eta2pi0_pBe->Draw("same,pz");
   }
   if (!EPOS || (EPOS && TAPS)){ 
  DrawGammaSetMarkerTGraphAsym(EtaPi07TeVStat, 25, 1,kBlue, kBlue);
  EtaPi07TeVStat->Draw("same,pz");
  DrawGammaSetMarkerTGraphAsym(EtaPi07TeVSyst, 25, 1, kBlue, kBlue, 1., kTRUE);
  EtaPi07TeVSyst->Draw("same,E2");
   }  
  DrawGammaSetMarkerTGraphAsym(EtaPi0pPbStat,20,1,1,1);  
  EtaPi0pPbStat->Draw("pz,same");
  DrawGammaSetMarkerTGraphAsym(EtaPi0pPbSyst,20,1, 1,1, 1, kTRUE);  
  EtaPi0pPbSyst->Draw("E2,same");
  TLegend* leg4;
  if (TAPS)   leg4 = new TLegend(0.12,0.6,0.4,0.95);
  else  if (EPOS)   leg4 = new TLegend(0.25,0.2,0.55,0.35);
  else leg4 = new TLegend(0.25,0.15,0.55,0.38);
  leg4->SetFillColor(0);
  leg4->SetLineColor(0);
  leg4->SetTextFont(Font);
  leg4->SetTextSize(0.045);
  leg4->AddEntry(EtaPi0pPbSyst,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  if (!EPOS || (EPOS && TAPS))leg4->AddEntry(EtaPi07TeVSyst,"pp, #sqrt{#it{s}} = 7 TeV","pef");// (PLB717 (2012) 162)","pef"); 
  if (EPOS) leg4->AddEntry(graphEtaPi0EPOS,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV - EPOS3","l"); 
  if (TAPS) {leg4->AddEntry( eta2pi0_pAu,"p-Au,  #sqrt{#it{s}_{NN}} = 29.1 GeV","p");  
    leg4->AddEntry( eta2pi0_pBe,"p-Be,  #sqrt{#it{s}_{NN}} = 29.1 GeV","p");  
       // leg4a->AddEntry(  eta2pi0_pBe,"(Eur. Phys. J. C 4, 249–257 (1998))","l");   
  } 
  if (!EPOS || (EPOS && TAPS))leg4->AddEntry(mTScaling,"#eta from #it{m}_{T} scaled #pi^{0}","pl");  
  
  leg4->Draw("same");


  c4->Update();
  if (EPOS) c4->Print(Form("%s/EtaPi0Ratio_EPOS.%s",outputDir.Data(),suffix.Data()));
  else if (TAPS) c4->Print(Form("%s/EtaPi0Ratio_TAPS.%s",outputDir.Data(),suffix.Data()));
  else c4->Print(Form("%s/EtaPi0Ratio.%s",outputDir.Data(),suffix.Data()));
//Eta/Pi0 Ratio Logx

  

  TCanvas* c4a = new TCanvas("c4a","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c4a, CMarginL, CMarginR,CMarginT ,CMarginB);
  c4a->SetLogx(); 
  // c4a->SetLogy();
  TH2F * hist4a;
  if (TAPS)   hist4a = new TH2F("hist4a","hist4a",1000,0.09,20.,1000,0.01,.99   );
//else   hist4a = new TH2F("hist4a","hist4a",1000,0.6,20.,1000,0.0,.99   );
  else hist4a = new TH2F("hist4a","hist4a",1000,0.3,20.,1000,0.0,.99   );
  //  hist4a = new TH2F("hist4a","hist4a",1000,0.03,7.,1000,0.007,1.   );
  SetStyleHistoTH2ForGraphs(hist4a, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0} ",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY, TitleOffsetX,TitleOffsetY, 512, 508);
  hist4a->GetYaxis()->SetLabelOffset(0.01);
  hist4a->GetXaxis()->SetLabelOffset(LabelOffsetLog); 
  hist4a->DrawCopy(); 
  
  graphMcGillEtaToPi0->SetLineColor(colorMcGill);
   graphMcGillEtaToPi0->SetFillColor(colorMcGill);
   graphMcGillEtaToPi0->SetLineWidth(3);
   graphMcGillEtaToPi0->Draw("C3same");

   if (EPOS){
    graphEtaPi0EPOS->SetLineColor(colorEPOS);
    graphEtaPi0EPOS->SetFillColor(colorEPOS);
    graphEtaPi0EPOS->SetLineWidth(3);
    graphEtaPi0EPOS->Draw("C3same");
  }  
  
   if (!EPOS || (EPOS && TAPS)) mTScaling->Draw("same"); //NOTE TEMP
 
   if (TAPS){
     DrawGammaSetMarkerTGraphErr(eta2pi0_pAu, 27, 1.4,kGreen, kGreen);
     DrawGammaSetMarkerTGraphErr(eta2pi0_pBe, 27, 1.4,kRed, kRed);
     eta2pi0_pAu->Draw("same,zp");
     eta2pi0_pBe->Draw("same,pz");
   }
   if( HIJINGPi0 ){
   SetStyleHisto(histoDPMJetEtaToPi0, 3, styleLineDPMJet, colorDPMJet );  
   histoDPMJetEtaToPi0->Draw("same,hist,l");
   }

   if( DPMJetPi0 ){
   SetStyleHisto(histoHIJINGEtaToPi0,3, 8, colorHIJING);  
   histoHIJINGEtaToPi0->Draw("same,hist,l");
   }
   
   
   

  if (!EPOS || (EPOS && TAPS)){ 
  DrawGammaSetMarkerTGraphAsym(EtaPi07TeVStat, 34, 1, colorEtaPi07TeV, colorEtaPi07TeV);
  EtaPi07TeVStat->Draw("same,zp");
  DrawGammaSetMarkerTGraphAsym(EtaPi07TeVSyst, 34, 1, colorEtaPi07TeV, colorEtaPi07TeV, 1., kTRUE);
  EtaPi07TeVSyst->Draw("same,E2");
   }  
  DrawGammaSetMarkerTGraphAsym(EtaPi0pPbStat,20,1,1,1);  
  EtaPi0pPbStat->Draw("pz,same");
  DrawGammaSetMarkerTGraphAsym(EtaPi0pPbSyst,20,1, 1,1, 1, kTRUE);  
  EtaPi0pPbSyst->Draw("E2,same");
  
  
  
  TLatex * lt4a = new TLatex(.12,0.92,"ALICE");
  lt4a->SetTextColor(kBlack);
  lt4a->SetTextSize(TextSize);
  lt4a->SetTextFont(Font);
  lt4a->Draw() ;
  
  TLegend* leg4a;
  if (EPOS && TAPS)leg4a = new TLegend(0.15,0.56,0.45,0.90);
  else if (EPOS)leg4a = new TLegend(0.15,0.75,0.45,0.95);
  else if (TAPS)leg4a = new TLegend(0.15,0.56,0.45,0.95);
  else leg4a = new TLegend(0.15,0.60,0.45,0.9);
  leg4a->SetFillColor(0);
  leg4a->SetLineColor(0);
  leg4a->SetTextFont(Font);
  leg4a->SetTextSize(0.041);
  leg4a->AddEntry(EtaPi0pPbSyst,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  if (!EPOS || (EPOS && TAPS))leg4a->AddEntry(EtaPi07TeVSyst,"pp, #sqrt{#it{s}} = 7 TeV","pef");// (PLB717 (2012) 162)","pef");  
  if (!EPOS || (EPOS && TAPS)){ 
    leg4a->AddEntry(mTScaling,"#eta from #it{m}_{T} scaled #pi^{0}","pl"); 
  }
  if (EPOS) leg4a->AddEntry(graphEtaPi0EPOS,"EPOS3","l");
   leg4a->AddEntry(graphMcGillEtaToPi0,"VISHNU","l");
  if( HIJINGPi0){
    leg4a->AddEntry(histoHIJINGEtaToPi0,"HIJING","l"); 
  }
  
  
  
  if( DPMJetPi0 ){
     leg4a->AddEntry(histoDPMJetEtaToPi0,"DPMJet","l"); 
  }
   leg4a->Draw("same");
  
  if (TAPS) {
    
     TLatex * lt4b = new TLatex(0.12,0.43,"CERES-TAPS #sqrt{#it{s}_{NN}} = 29.1 GeV");
     lt4b->SetTextColor(kBlack);
     lt4b->SetTextSize(TextSize);
     lt4b->SetTextFont(Font);
     lt4b->Draw();
     
     TLegend* leg4b;
     if (EPOS && TAPS){leg4b = new TLegend(0.15,0.38,0.35,0.47);
     leg4b->SetFillColor(0);
     leg4b->SetLineColor(0);
     leg4b->SetTextFont(Font);
     leg4b->SetTextSize(0.041);
     leg4b->AddEntry( eta2pi0_pAu,"p-Au","p");  
     leg4b->AddEntry( eta2pi0_pBe,"p-Be","p");  
     leg4b->Draw("same");}
  }  
  
  



  c4a->Update();
  if (EPOS && TAPS)  	c4a->Print(Form("%s/EtaPi0Ratio_LogX_EPOS_TAPS.%s",outputDir.Data(),suffix.Data()));
  else if (EPOS)  	c4a->Print(Form("%s/EtaPi0Ratio_LogX_EPOS.%s",outputDir.Data(),suffix.Data()));
  else if (TAPS)  	c4a->Print(Form("%s/EtaPi0Ratio_LogX_TAPS.%s",outputDir.Data(),suffix.Data()));
  else  		c4a->Print(Form("%s/EtaPi0Ratio_LogX.%s",outputDir.Data(),suffix.Data()));

  

  TCanvas* c4b = new TCanvas("c4b","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c4b, CMarginL, CMarginR,CMarginT ,CMarginB);

  TH2F * hist4b;
  hist4b = new TH2F("hist4b","hist4b",1000,0.,20.,1000,0.01,.99   );
  SetStyleHistoTH2ForGraphs(hist4b, "#it{m}_{T} (GeV/#it{c})","#eta/#pi^{0} ",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY, TitleOffsetX,TitleOffsetY, 512, 508);
  hist4b->DrawCopy(); 
    
/* TEMP
  DrawGammaSetMarkerTGraphAsym(EtaPi07TeVStat_mT, 25, 1,kBlue, kBlue);
  EtaPi07TeVStat_mT->Draw("same,zp");
  DrawGammaSetMarkerTGraphAsym(EtaPi07TeVSyst_mT, 25, 1, kBlue, kBlue, 1., kTRUE);
  EtaPi07TeVSyst_mT->Draw("same,E2");
   
  DrawGammaSetMarkerTGraphAsym(EtaPi0pPbStat_mT,20,1,1,1);  
  EtaPi0pPbStat_mT->Draw("pz,same");
  DrawGammaSetMarkerTGraphAsym(EtaPi0pPbSyst_mT,20,1, 1,1, 1, kTRUE);  
  EtaPi0pPbSyst_mT->Draw("E2,same");*/
     

  TLegend* leg4b = new TLegend(0.25,0.15,0.55,0.33);
  leg4b->SetFillColor(0);
  leg4b->SetLineColor(0);
  leg4b->SetTextFont(Font);
  leg4b->SetTextSize(0.045);
  leg4b->AddEntry(EtaPi0pPbSyst,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  leg4b->AddEntry(EtaPi07TeVSyst,"pp, #sqrt{#it{s}} = 7 TeV (PLB717 (2012) 162)","pef");
  // leg4b->AddEntry(mTScaling,"#eta from #it{m}_{T} scaled #pi^{0}","pl");  
  leg4b->Draw("same");
 TLine* line4b=new TLine(0.,0.47,20.,0.47);
  line4b->SetLineColor(kRed+2);
  line4b->SetLineStyle(2);
  line4b->Draw("same");

  c4b->Update();
  c4b->Print(Form("%s/EtaPi0Ratio_mT.%s",outputDir.Data(),suffix.Data()));





  //Eta/Pi0 Ratio vs mT logx

  

  TCanvas* c4c = new TCanvas("c4c","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c4c, CMarginL, CMarginR,CMarginT ,CMarginB);
  c4c->SetLogx();
  TH2F * hist4c;
  hist4c = new TH2F("hist4c","hist4c",1000,0.6,20.,1000,0.0,.99   );
  SetStyleHistoTH2ForGraphs(hist4c, "#it{m}_{T} (GeV/#it{c})","#eta/#pi^{0} ",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY, TitleOffsetX,TitleOffsetY, 512, 508);
  hist4c->GetXaxis()->SetLabelOffset(LabelOffsetLog); 
  hist4c->DrawCopy(); 
    
/*   TEMP
  DrawGammaSetMarkerTGraphAsym(EtaPi07TeVStat_mT, 25, 1,kBlue, kBlue);
  EtaPi07TeVStat_mT->Draw("same,zp");
  DrawGammaSetMarkerTGraphAsym(EtaPi07TeVSyst_mT, 25, 1, kBlue, kBlue, 1., kTRUE);
  EtaPi07TeVSyst_mT->Draw("same,E2");
   
  DrawGammaSetMarkerTGraphAsym(EtaPi0pPbStat_mT,20,1,1,1);  
  EtaPi0pPbStat_mT->Draw("pz,same");
  DrawGammaSetMarkerTGraphAsym(EtaPi0pPbSyst_mT,20,1, 1,1, 1, kTRUE);  
  EtaPi0pPbSyst_mT->Draw("E2,same");*/
     

  TLegend* leg4c = new TLegend(0.15,0.70,0.45,0.90);
  leg4c->SetFillColor(0);
  leg4c->SetLineColor(0);
  leg4c->SetTextFont(Font);
  leg4c->SetTextSize(0.045);
  //leg4c->AddEntry(EtaPi0pPbSyst,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV","pef"); TEMP
  //leg4c->AddEntry(EtaPi07TeVSyst,"pp, #sqrt{#it{s}} = 7 TeV (PLB717 (2012) 162)","pef"); TEMP
  // leg4c->AddEntry(mTScaling,"#eta from #it{m}_{T} scaled #pi^{0}","pl");  
  leg4c->Draw("same");
 TLine* line4c=new TLine(0.6,0.47,20.,0.47);
  line4c->SetLineColor(kRed+2);
  line4c->SetLineStyle(2);
  line4c->Draw("same");

  c4c->Update();
  c4c->Print(Form("%s/EtaPi0Ratio_mT_LogX.%s",outputDir.Data(),suffix.Data()));

  

  TCanvas* c4d = new TCanvas("c4d","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c4d, CMarginL+0.01, CMarginR,CMarginT ,CMarginB);
  c4d->SetLogx();
  TH2F * hist4d;
  hist4d = new TH2F("hist4d","hist4d",1000,0.6,20.,1000,0.21,1.39   );
  SetStyleHistoTH2ForGraphs(hist4d, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}_{measured} / #eta/#pi^{0}_{#it{m}_{T} scaling}",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY, TitleOffsetX,TitleOffsetY+0.05, 512, 508);
  hist4d->GetXaxis()->SetLabelOffset(LabelOffsetLog); 
 
  hist4d->DrawCopy(); 
  
     //TEMP
   TGraphAsymmErrors* EtaPi0Stat_div_mT=(TGraphAsymmErrors*)EtaPi0pPbStat->Clone("EtaPi0Stat_div_mT");
   TGraphAsymmErrors* EtaPi0Syst_div_mT=(TGraphAsymmErrors*)EtaPi0pPbSyst->Clone("EtaPi0Syst_div_mT");
   EtaPi0Stat_div_mT=CalculateGraphErrRatioToFit (EtaPi0pPbStat ,mTScaling  );  
   EtaPi0Syst_div_mT=CalculateGraphErrRatioToFit (EtaPi0pPbSyst ,mTScaling  );  
 
   DrawGammaSetMarkerTGraphAsym(EtaPi0Stat_div_mT, 20, 1,kBlack, kBlack);
   EtaPi0Stat_div_mT->Draw("same,zp");
   DrawGammaSetMarkerTGraphAsym(EtaPi0Syst_div_mT, 20, 1, kBlack, kBlack, 1., kTRUE);
   EtaPi0Syst_div_mT->Draw("same,E2"); //TEMP
   
   


  TLegend* leg4d = new TLegend(0.45,0.25,0.95,0.35);
  leg4d->SetFillColor(0);
  leg4d->SetLineColor(0);
  leg4d->SetTextFont(Font);
  leg4d->SetTextSize(0.045);
 // leg4d->AddEntry(EtaPi0Syst_div_mT,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV","pef"); TEMP
  // leg4d->AddEntry(mTScaling,"#eta from #it{m}_{T} scaled #pi^{0}","pl");  
   leg4d->Draw("same");
   
  DrawGammaLines(0.6,20.,1.,1.,2.0,kRed+2,2);
  
  
  //TLatex *Prelim4d=new TLatex(0.15,0.9,"ALICE Preliminary");
  //Prelim4d->SetNDC();
  //Prelim4d->SetTextColor(1);
  //Prelim4d->SetTextSize(0.05);
  //Prelim4d->SetTextFont(42);
  //  Prelim4d->Draw("same");  
  TLatex * Prelim4d = new TLatex(0.7,1.3,"ALICE");
  Prelim4d->SetTextColor(kBlack) ;
  Prelim4d->SetTextSize(TextSize) ;
  Prelim4d->SetTextFont(Font) ;
  Prelim4d->DrawLatex(0.7,1.22,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  
  Prelim4d->Draw() ;

  
  c4d->Update();
  c4d->Print(Form("%s/EtaPi0Ratio2mTScaling.%s",outputDir.Data(),suffix.Data()));




  //RpPb

  TCanvas* c5 = new TCanvas("c5","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c5, CMarginL, CMarginR,CMarginT ,CMarginB);

  TH2F * hist5;
  hist5 = new TH2F("hist5","hist5",1000,0.,22.,1000,0.3,1.5   );
  SetStyleHistoTH2ForGraphs(hist5, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#pi^{0}}_{p-Pb}",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY,TitleOffsetX,TitleOffsetY, 512, 508);
  hist5->DrawCopy(); 
 
  TLine *line5=new TLine(0.,1.,20.,1.);
  line5->SetLineColor(kGray+1);
  line5->SetLineColor(kRed+2);
  line5->SetLineStyle(2);
  line5->Draw();
     
  EPS09s_fDSS_errors->Draw("same,E3");
  CGC->SetMarkerSize(0.65);
  CGC->Draw("same,p");
  EPS09s_fDSS_NLO->Draw("same,l");
  EPS09s_AKK_NLO->Draw("same,l");
  EPS09s_AKK_NLO->SetLineColor(kGreen+3);
  EPS09s_KKP_NLO->Draw("same,l");
  CombinedPi0RpPbSystErr->SetMarkerSize(0.65);
  CombinedPi0RpPbSystErr->Draw("E2,same");
  CombinedPi0RpPbStatErr->SetMarkerSize(0.65);
  CombinedPi0RpPbStatErr->Draw("pz,same"); 
  cout <<"=========RpPb binning================="<< endl;  
  CombinedPi0RpPbStatErr->Print();
  cout <<"=========================="<< endl;
  TLegend* leg5a = new TLegend(0.14,0.16,0.5,0.39);
  leg5a->SetFillColor(0);
  leg5a->SetLineColor(0);
  leg5a->SetTextFont(Font);
  leg5a->SetTextSize(TextSizeRpA);
  leg5a->AddEntry(EPS09s_KKP_NLO,"EPS09s KKP NLO","l");
  leg5a->AddEntry(EPS09s_AKK_NLO,"EPS09s AKK NLO","l");
  leg5a->AddEntry(EPS09s_fDSS_NLO,"EPS09s fDSS NLO","l");
  leg5a->AddEntry(EPS09s_fDSS_errors,"EPS09s fDSS errors","ef");
  leg5a->AddEntry((TObject*)0,"JHEP 1207 (2012) 073","");
  leg5a->Draw("same");

  TLegend* leg5b = new TLegend(0.51,0.31,0.64,0.39);
  leg5b->SetFillColor(0);
  leg5b->SetLineColor(0);
  leg5b->SetTextFont(Font);
  leg5b->SetTextSize(TextSizeRpA);

  leg5b->AddEntry(CGC,"CGC","p");
  leg5b->AddEntry((TObject*)0,"Phys.Rev. D88 (2013) 114020","");
  //leg5b->AddEntry((TObject*)0,"114020","");
  leg5b->Draw("same");
  
  TLatex * lt4 = new TLatex(1.5,1.4,"ALICE");
  lt4->SetTextColor(kBlack) ;
  lt4->SetTextSize(TextSizeRpA) ;
  lt4->SetTextFont(Font) ;
  //lt4->DrawLatex(2.8,2.25,"ALICE, NSD #eta");
  lt4->DrawLatex(1.5,1.32,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  
  lt4->Draw() ;

  TLegend* leg5c = new TLegend(0.15,0.76,0.50,0.83);
  leg5c->SetFillColor(0);
  leg5c->SetLineColor(0);
  leg5c->SetTextFont(Font);
  leg5c->SetTextSize(TextSizeRpA);

  leg5c->AddEntry(CombinedPi0RpPbSystErr,"NSD #pi^{0}","pef");
  leg5c->Draw("same");

    Double_t NormalizationError =TMath::Sqrt(pow(0.031,2)+pow(0.036,2)+pow(0.036,2));//pPb normalization,TpPb and pp normalization errors
  TBox* BoxNorm =new TBox(0.45, 1.-NormalizationError, 0.7, 1.+NormalizationError);
  BoxNorm->SetFillColor(kBlack);
  BoxNorm->Draw("same");

  c5->Update();
  c5->Print(Form("%s/CombRpA_Models.%s",outputDir.Data(),suffix.Data()));

  //RpPb Version2
  TCanvas* c5a = new TCanvas("c5a","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c5a, CMarginL, CMarginR,CMarginT ,CMarginB);

  TH2F * hist5a;
  hist5a = new TH2F("hist5a","hist5a",1000,0.,21.,1000,0.55,1.65   );
  SetStyleHistoTH2ForGraphs(hist5a, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#pi^{0}}_{p-Pb}",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY,TitleOffsetX,TitleOffsetY, 512, 508);
  hist5a->DrawCopy(); 
 
  TLine *line5a=new TLine(0.,1.,20.,1.);
  line5a->SetLineColor(kGray+1);
  line5a->Draw();
     
  EPS09s_fDSS_errors->Draw("same,E3");
  CGC->SetMarkerSize(0.65);
  CGC->Draw("same,p");
  EPS09s_fDSS_NLO->Draw("same,l");
  EPS09s_AKK_NLO->Draw("same,l");
  EPS09s_AKK_NLO->SetLineColor(kGreen+3);
  EPS09s_KKP_NLO->Draw("same,l");
  CombinedPi0RpPbSystErr->SetMarkerSize(0.65);
  CombinedPi0RpPbSystErr->Draw("E2,same");
  CombinedPi0RpPbStatErr->SetMarkerSize(0.65);
  CombinedPi0RpPbStatErr->Draw("pz,same");   

  TLegend* leg5aa = new TLegend(0.17,0.7,0.61,0.93);
  leg5aa->SetFillColor(0);
  leg5aa->SetLineColor(0);
  leg5aa->SetTextFont(Font);
  leg5aa->SetTextSize(TextSizeRpA);
  leg5aa->AddEntry(EPS09s_KKP_NLO,"EPS09s KKP NLO","l");
  leg5aa->AddEntry(EPS09s_AKK_NLO,"EPS09s AKK NLO","l");
  leg5aa->AddEntry(EPS09s_fDSS_NLO,"EPS09s fDSS NLO","l");
  leg5aa->AddEntry(EPS09s_fDSS_errors,"EPS09s fDSS errors","ef");
  leg5aa->AddEntry((TObject*)0,"JHEP 1207 (2012) 073","");
  leg5aa->Draw("same");

  TLegend* leg5ab = new TLegend(0.6,0.78,0.85,0.93);
  leg5ab->SetFillColor(0);
  leg5ab->SetLineColor(0);
  leg5ab->SetTextFont(Font);
  leg5ab->SetTextSize(TextSizeRpA);

  leg5ab->AddEntry(CGC,"CGC","p");
  leg5ab->AddEntry((TObject*)0,"Phys.Rev. D88 (2013)","");
  leg5ab->AddEntry((TObject*)0,"114020","");
  leg5ab->Draw("same");

  BoxNorm->Draw("same");
  CombinedPi0RpPbSystErr->Draw("E2,same");
  c5a->Update();
  c5a->Print(Form("%s/CombRpA_Models_V2.%s",outputDir.Data(),suffix.Data()));
  
  
  //RpPb version 3
  
  TCanvas* c5aa = new TCanvas("c5aa","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c5aa, CMarginL, CMarginR,CMarginT ,CMarginB);

  TH2F * hist5aa;
  hist5aa = new TH2F("hist5aa","hist5aa",1000,0.,22.,1000,0.50,1.5   );
  SetStyleHistoTH2ForGraphs(hist5aa, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#pi^{0}}_{p-Pb}",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY,TitleOffsetX,TitleOffsetY, 512, 508);
  hist5aa->DrawCopy(); 
  
  TLegend* leg5abcde = new TLegend(0.16,0.77,0.41,0.82);
  leg5abcde->SetFillColor(0);
  leg5abcde->SetLineColor(0);
  leg5abcde->SetTextFont(Font);
  leg5abcde->SetTextSize(TextSizeRpA);
  leg5abcde->AddEntry(CombinedPi0RpPbSystErr,"#pi^{0}","pef");
  //leg5abcde->AddEntry(graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14,"NLO pQCD","ef");
  //leg5abcde->AddEntry(CGC,"CGC","p");
  
  leg5abcde->Draw("same");
 
 
     
  graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14->SetFillColor(kYellow-7);
  graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14->SetFillStyle(1001);
  graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14->Draw("same,E3");
  
  CGC->SetMarkerSize(0.65);
  CGC->Draw("same,p");
  
  
  DrawGammaSetMarkerTGraphAsym(CombinedPi0RpPbStatErr, 20, 1, kBlack, kBlack);
  DrawGammaSetMarkerTGraphAsym(CombinedPi0RpPbSystErr, 20, 1, kBlack, kBlack, 1., kTRUE);
   
  
  //CombinedPi0RpPbSystErr->SetMarkerSize(0.65);
  CombinedPi0RpPbSystErr->Draw("E2,same");
  //CombinedPi0RpPbStatErr->SetMarkerSize(0.65);
  CombinedPi0RpPbStatErr->Draw("pz,same"); 
  
 
  cout <<"=========RpPb binning================="<< endl;  
  CombinedPi0RpPbStatErr->Print();
  cout <<"=========================="<< endl;
  TLegend* leg5abc = new TLegend(0.31,0.25,0.54,0.34);
  leg5abc->SetFillColor(0);
  leg5abc->SetLineColor(0);
  leg5abc->SetTextFont(Font);
  leg5abc->SetTextSize(TextSizeRpA);
  leg5abc->AddEntry(graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14,"NLO pQCD","ef");
  leg5abc->AddEntry(CGC,"CGC","p");
  //leg5abc->AddEntry((TObject*)0,"JHEP 1207 (2012) 073","");
  leg5abc->Draw("same");

  TLegend* leg5abcd = new TLegend(0.51,0.38,0.64,0.41);
  leg5abcd->SetFillColor(0);
  leg5abcd->SetLineColor(0);
  leg5abcd->SetTextFont(Font);
  leg5abcd->SetTextSize(TextSizeRpA);

  //leg5abc->AddEntry(CGC,"CGC","p");
  //leg5abc->AddEntry((TObject*)0,"Phys.Rev. D88 (2013) 114020","");
  //leg5abc->Draw("same");
  
  
  
  
  
  
  
  
  
  
  
  
  TLatex * lt4ab = new TLatex(1.5,1.4,"ALICE");
  lt4ab->SetTextColor(kBlack) ;
  lt4ab->SetTextSize(TextSize) ;
  lt4ab->SetTextFont(Font) ;
  lt4ab->DrawLatex(1.5,1.35,"p-Pb, NSD, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  lt4ab->Draw();

  

  BoxNorm->Draw("same");
  DrawGammaLines(0.,22.,1.,1.,2.0,kRed+2,2);
  
//   TLine *line5aa=new TLine(0.,1.,22.,1.);
//   line5aa->SetLineColor(kGray+1);
//   line5aa->SetLineColor(kRed+2);
//   line5aa->SetLineStyle(2);
//   line5aa->Draw();

  c5aa->Update();
  c5aa->Print(Form("%s/CombRpA_Models_V3.%s",outputDir.Data(),suffix.Data()));
  
    
  
  
  
  


  //Width Pi0
  TCanvas* c10 = new TCanvas("c10","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c10, CMarginL, CMarginR,CMarginT ,CMarginB);
  c10->SetLogx();
  TH2F * hist10;
  hist10 = new TH2F("hist10","hist10",1000,0.3,30.,1000,-2.,25.5);
  SetStyleHistoTH2ForGraphs(hist10, "#it{p}_{T} (GeV/#it{c})","peak width (MeV/#it{c}^{2})",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY,TitleOffsetX,TitleOffsetY, 512, 508);
  hist10->GetXaxis()->SetLabelOffset(LabelOffsetLog);
  hist10->DrawCopy(); 
	 DrawGammaSetMarker(PHOSPi0WidthMC, markerStylePHOSMC,markerSizePHOS , colorPHOSMC , colorPHOSMC);
	 PHOSPi0WidthMC->DrawCopy("same,z,p,x0");      
 	 DrawGammaSetMarker(PHOSPi0Width, markerStylePHOS,markerSizePHOS , colorPHOS, colorPHOS);
	 PHOSPi0Width->DrawCopy("same,z,p,x0"); 

	 DrawGammaSetMarker(PCMPi0WidthMC, markerStylePCMMC, markerSizePCM, colorPCMMC , colorPCMMC);
	 PCMPi0WidthMC->DrawCopy("same,z,p,x0");
	 DrawGammaSetMarker(PCMPi0Width, markerStylePCM,markerSizePCM , colorPCM, colorPCM);
	 PCMPi0Width->DrawCopy("same,z,p,x0"); 

 	 DrawGammaSetMarker(DalitzPi0WidthMC, markerStyleDalitzMC,markerSizeDalitz , colorDalitzMC , colorDalitzMC);
	 DalitzPi0WidthMC->DrawCopy("same,z,p,x0");    
	 DrawGammaSetMarker(DalitzPi0Width, markerStyleDalitz,markerSizeDalitz , colorDalitz, colorDalitz);
	 DalitzPi0Width->DrawCopy("same,z,p,x0"); 

 	 DrawGammaSetMarker(EMCalPi0WidthMC, markerStyleEMCalMC,markerSizeEMCal , colorEMCalMC , colorEMCalMC);
	 EMCalPi0WidthMC->DrawCopy("same,z,p,x0");    
	 DrawGammaSetMarker(EMCalPi0Width, markerStyleEMCal,markerSizeEMCal , colorEMCal, colorEMCal);
	 EMCalPi0Width->DrawCopy("same,z,p,x0"); 



	//********************************** Defintion of the Legend ************************************************** 
	 Double_t textSizeMassWidth =0.049;  
	Double_t columnsLegendFWHM2[4]    = {0.15,0.426,0.515,0.39};
	Double_t columnsLegendFWHM2Abs[4]    = {4,1.9,2.88,12};
	//    Double_t rowsLegendFWHM2[3]       = {0.66,0.33,0.0};
	Double_t rowsLegendFWHM2[7]       = {0.91,0.86,0.81,0.76, 0.71,0.68, 0.63};
	Double_t rowsLegendFWHM2Abs[7]       = {0.2,22.,20.45,18.9, 17.35, 18.4, 17.0};
	//******************* Text sizes *******************
	Size_t textSizeLeftColumnFWHM2 = textSizeMassWidth*0.85;
	Size_t textSizeTopRowFWHM2  = textSizeMassWidth*0.85; 
	//****************** Scale factors ******************
	Double_t scaleMarkerFWHM2      = 1.1;
	
	//    padFWHMLegend1->cd();
	//****************** first Column **************************************************
	TLatex *textFWHM2PCM = new TLatex(columnsLegendFWHM2[0],rowsLegendFWHM2[3],"PCM (FWHM/2.35)");
	SetStyleTLatex( textFWHM2PCM, textSizeLeftColumnFWHM2,4);
	textFWHM2PCM->SetTextFont(42);
	textFWHM2PCM->Draw();

	TLatex *textFWHM2Dalitz = new TLatex(columnsLegendFWHM2[0],rowsLegendFWHM2[4],"PCM-#gamma*#gamma (FWHM/2.35)");
	SetStyleTLatex( textFWHM2Dalitz, textSizeLeftColumnFWHM2,4);
	textFWHM2Dalitz->SetTextFont(42);
	textFWHM2Dalitz->Draw();
	
 	TLatex *textFWHMEMCALEMCAL = new TLatex(columnsLegendFWHM2[0],rowsLegendFWHM2[2],"EMCal (FWHM/2.35)");
 	SetStyleTLatex( textFWHMEMCALEMCAL, textSizeLeftColumnFWHM2,4);
 	textFWHMEMCALEMCAL->SetTextFont(42);
 	textFWHMEMCALEMCAL->Draw();
	TLatex *textFWHM2PHOS = new TLatex(columnsLegendFWHM2[0],rowsLegendFWHM2[1],"PHOS (#sigma)");
	SetStyleTLatex( textFWHM2PHOS, textSizeLeftColumnFWHM2,4);
	textFWHM2PHOS->SetTextFont(42);
	textFWHM2PHOS->Draw();


	
	//****************** second Column *************************************************
	TLatex *textFWHM2Data2 = new TLatex(columnsLegendFWHM2[1],rowsLegendFWHM2[0] ,"data");
	SetStyleTLatex( textFWHM2Data2, textSizeTopRowFWHM2 ,4);
	textFWHM2Data2->SetTextFont(42);
		textFWHM2Data2->Draw();
	TLatex *textFWHM2MC2 = new TLatex(columnsLegendFWHM2[2]
	,rowsLegendFWHM2[0],"MC");
	SetStyleTLatex( textFWHM2MC2, textSizeTopRowFWHM2,4);
	textFWHM2MC2->SetTextFont(42);
	textFWHM2MC2->Draw();
	
	TMarker* markerPCMPi0FWHM2 = CreateMarkerFromHisto(PCMPi0Width,columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2Abs[3] ,scaleMarkerFWHM2);
	markerPCMPi0FWHM2->Draw();
 	TMarker* markerEMCALEMCALPi0FWHM = CreateMarkerFromHisto(EMCalPi0Width,columnsLegendFWHM2Abs[1],rowsLegendFWHM2Abs[2],scaleMarkerFWHM2);
 	markerEMCALEMCALPi0FWHM->Draw();
	TMarker* markerDalitzPi0FWHM2 = CreateMarkerFromHisto(DalitzPi0Width,columnsLegendFWHM2Abs[1],rowsLegendFWHM2Abs[4],scaleMarkerFWHM2);
	markerDalitzPi0FWHM2->Draw();
	TMarker* markerPHOSPi0FWHM2 = CreateMarkerFromHisto(PHOSPi0Width,columnsLegendFWHM2Abs[1],rowsLegendFWHM2Abs[1],scaleMarkerFWHM2);
	markerPHOSPi0FWHM2->Draw();

	
	TMarker* markerPCMPi0FWHM2MC = CreateMarkerFromHisto(PCMPi0WidthMC,columnsLegendFWHM2Abs[2],rowsLegendFWHM2Abs[3],scaleMarkerFWHM2);
	markerPCMPi0FWHM2MC->Draw();
 	TMarker* markerEMCALEMCALPi0FWHMMC = CreateMarkerFromHisto(EMCalPi0WidthMC,columnsLegendFWHM2Abs[2] ,rowsLegendFWHM2Abs[2] ,scaleMarkerFWHM2);
 	markerEMCALEMCALPi0FWHMMC->Draw();
	TMarker* markerDalitzPi0FWHM2MC = CreateMarkerFromHisto(DalitzPi0WidthMC,columnsLegendFWHM2Abs[2] ,rowsLegendFWHM2Abs[4] ,scaleMarkerFWHM2);
	markerDalitzPi0FWHM2MC->Draw();
	TMarker* markerPHOSPi0FWHM2MC = CreateMarkerFromHisto(PHOSPi0WidthMC,columnsLegendFWHM2Abs[2] ,rowsLegendFWHM2Abs[1] ,scaleMarkerFWHM2);
	markerPHOSPi0FWHM2MC->Draw();






  c10->Update();
  c10->Print(Form("%s/Pi0Width.%s",outputDir.Data(),suffix.Data()));

  //WidthEta
  TCanvas* c11 = new TCanvas("c11","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c11, CMarginL, CMarginR,CMarginT ,CMarginB);
  c11->SetLogx();
  TH2F * hist11;
  hist11 = new TH2F("hist11","hist11",1000,0.3,30.,1000,-5.,75.);
  SetStyleHistoTH2ForGraphs(hist11, "#it{p}_{T} (GeV/#it{c})","peak width (MeV/#it{c}^{2})",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY,TitleOffsetX,TitleOffsetY, 512, 508);
  hist11->GetXaxis()->SetLabelOffset(LabelOffsetLog);
  hist11->DrawCopy(); 

	 DrawGammaSetMarker(PCMEtaWidthMC, markerStylePCMMC, markerSizePCM, colorPCMMC , colorPCMMC);
	 PCMEtaWidthMC->SetBinError(PCMEtaWidthMC->FindBin(3.5),1.);//Check later with new output (large error was not observed in previous versions)
	 PCMEtaWidthMC->DrawCopy("same,z,p,x0");      
         DrawGammaSetMarker(PCMEtaWidth, markerStylePCM,markerSizePCM , colorPCM, colorPCM);
	 PCMEtaWidth->DrawCopy("same,z,p,x0"); 

	 DrawGammaSetMarker(EMCalEtaWidthMC, markerStyleEMCalMC,markerSizeEMCal , colorEMCalMC , colorEMCalMC);
	 EMCalEtaWidthMC->DrawCopy("same,z,p,x0");      
	 DrawGammaSetMarker(EMCalEtaWidth, markerStyleEMCal,markerSizeEMCal , colorEMCal, colorEMCal);
	 EMCalEtaWidth->DrawCopy("same,z,p,x0"); 


	 Double_t rowsLegendFWHM2AbsEta[7]       = {0.2,65.,60.3,21.2,19.7, 18.2, 18.2};

	//****************** first Column **************************************************
	TLatex *textFWHM2PCMEta = new TLatex(columnsLegendFWHM2[0],rowsLegendFWHM2[2],"PCM (FWHM/2.35)");
	SetStyleTLatex( textFWHM2PCMEta, textSizeLeftColumnFWHM2,4);
	textFWHM2PCMEta->SetTextFont(42);
	textFWHM2PCMEta->Draw();

	
 	TLatex *textFWHMEMCALEMCALEta = new TLatex(columnsLegendFWHM2[0],rowsLegendFWHM2[1],"EMCal (FWHM/2.35)");
 	SetStyleTLatex( textFWHMEMCALEMCALEta, textSizeLeftColumnFWHM2,4);
 	textFWHMEMCALEMCALEta->SetTextFont(42);
 	textFWHMEMCALEMCALEta->Draw();
	
	//****************** second Column *************************************************
	TLatex *textFWHM2Data2Eta = new TLatex(columnsLegendFWHM2[1],rowsLegendFWHM2[0] ,"data");
	SetStyleTLatex( textFWHM2Data2Eta, textSizeTopRowFWHM2 ,4);
	textFWHM2Data2Eta->SetTextFont(42);
		textFWHM2Data2Eta->Draw();
	TLatex *textFWHM2MC2Eta = new TLatex(columnsLegendFWHM2[2],rowsLegendFWHM2[0],"MC");
	SetStyleTLatex( textFWHM2MC2Eta, textSizeTopRowFWHM2,4);
	textFWHM2MC2Eta->SetTextFont(42);
	textFWHM2MC2Eta->Draw();
	
	TMarker* markerPCMEtaFWHM2Eta = CreateMarkerFromHisto(PCMEtaWidth,columnsLegendFWHM2Abs[1] ,rowsLegendFWHM2AbsEta[2] ,scaleMarkerFWHM2);
	markerPCMEtaFWHM2Eta->Draw();
 	TMarker* markerEMCALEMCALEtaFWHMEta= CreateMarkerFromHisto(EMCalEtaWidth,columnsLegendFWHM2Abs[1],rowsLegendFWHM2AbsEta[1],scaleMarkerFWHM2);
 	markerEMCALEMCALEtaFWHMEta->Draw();
	
	TMarker* markerPCMEtaFWHM2MCEta = CreateMarkerFromHisto(PCMEtaWidthMC,columnsLegendFWHM2Abs[2],rowsLegendFWHM2AbsEta[2],scaleMarkerFWHM2);
	markerPCMEtaFWHM2MCEta->Draw();
 	TMarker* markerEMCALEMCALEtaFWHMMCEta = CreateMarkerFromHisto(EMCalEtaWidthMC,columnsLegendFWHM2Abs[2] ,rowsLegendFWHM2AbsEta[1] ,scaleMarkerFWHM2);
 	markerEMCALEMCALEtaFWHMMCEta->Draw();

  c11->Update();
  c11->Print(Form("%s/EtaWidth.%s",outputDir.Data(),suffix.Data()));

  //Rec Mass Pi0
  TCanvas* c12 = new TCanvas("c12","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c12, CMarginL, CMarginR,CMarginT ,CMarginB);
  c12->SetLogx();
  TH2F * hist12;
  hist12 = new TH2F("hist12","hist12",1000,0.3,30.,1000,126.,149);
  SetStyleHistoTH2ForGraphs(hist12, "#it{p}_{T} (GeV/#it{c})","peak position (MeV/#it{c}^{2})",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY,TitleOffsetX,TitleOffsetY, 512, 508);
  hist12->GetXaxis()->SetLabelOffset(LabelOffsetLog);
  hist12->DrawCopy(); 

  TLine* linePi0=new TLine(0.3,mesonMassExpectPi0*1000,30,mesonMassExpectPi0*1000);
  linePi0->SetLineColor(kGray);
  linePi0->Draw("same");

 	 DrawGammaSetMarker(PHOSPi0MassMC, markerStylePHOSMC,markerSizePHOS , colorPHOSMC , colorPHOSMC);
	 PHOSPi0MassMC->DrawCopy("same,z,p,x0");     
 	 DrawGammaSetMarker(PHOSPi0Mass, markerStylePHOS,markerSizePHOS , colorPHOS, colorPHOS);
	 PHOSPi0Mass->DrawCopy("same,z,p,x0"); 

	 DrawGammaSetMarker(PCMPi0MassMC, markerStylePCMMC, markerSizePCM, colorPCMMC , colorPCMMC);
	 PCMPi0MassMC->DrawCopy("same,z,p,x0");
	 DrawGammaSetMarker(PCMPi0Mass, markerStylePCM,markerSizePCM , colorPCM, colorPCM);
	 PCMPi0Mass->DrawCopy("same,z,p,x0"); 

	 DrawGammaSetMarker(DalitzPi0MassMC, markerStyleDalitzMC,markerSizeDalitz , colorDalitzMC , colorDalitzMC);
	 DalitzPi0MassMC->DrawCopy("same,z,p,x0");     
	 DrawGammaSetMarker(DalitzPi0Mass, markerStyleDalitz,markerSizeDalitz , colorDalitz, colorDalitz);
	 DalitzPi0Mass->DrawCopy("same,z,p,x0"); 

 	 DrawGammaSetMarker(EMCalPi0MassMC, markerStyleEMCalMC,markerSizeEMCal , colorEMCalMC , colorEMCalMC);
	 EMCalPi0MassMC->DrawCopy("same,z,p,x0");    
	 DrawGammaSetMarker(EMCalPi0Mass, markerStyleEMCal,markerSizeEMCal , colorEMCal, colorEMCal);
	 EMCalPi0Mass->DrawCopy("same,z,p,x0"); 


	//********************************** Defintion of the Legend ************************************************** 
	 Double_t textSizeMassMass =0.049;  
	Double_t columnsLegendMass[4]    = {0.15,0.3,0.389,0.3};
	Double_t columnsLegendMassAbs[4]    = {4,.96,1.47,12};
	//    Double_t rowsLegendMass[3]       = {0.66,0.33,0.0};
	Double_t rowsLegendMass[7]       = {0.91,0.86,0.81,0.76, 0.71,0.68, 0.63};
	Double_t rowsLegendMassAbs[7]       = {0.2,146.2,144.9,143.6,142.3 , 18.2, 17.0};
	//******************* Text sizes *******************
	Size_t textSizeLeftColumnMass = textSizeMassMass*0.85;
	Size_t textSizeTopRowMass  = textSizeMassMass*0.85; 
	//****************** Scale factors ******************
	Double_t scaleMarkerMass      = 1.1;
	
	//    padMassLegend1->cd();
	//****************** first Column **************************************************
	TLatex *textMassPCM = new TLatex(columnsLegendMass[0],rowsLegendMass[3],"PCM");
	SetStyleTLatex( textMassPCM, textSizeLeftColumnMass,4);
	textMassPCM->SetTextFont(42);
	textMassPCM->Draw();

	TLatex *textMassDalitz = new TLatex(columnsLegendMass[0],rowsLegendMass[4],"PCM-#gamma*#gamma");
	SetStyleTLatex( textMassDalitz, textSizeLeftColumnMass,4);
	textMassDalitz->SetTextFont(42);
	textMassDalitz->Draw();
	
 	TLatex *textMassEMCALEMCAL = new TLatex(columnsLegendMass[0],rowsLegendMass[2],"EMCal");
 	SetStyleTLatex( textMassEMCALEMCAL, textSizeLeftColumnMass,4);
 	textMassEMCALEMCAL->SetTextFont(42);
 	textMassEMCALEMCAL->Draw();
	TLatex *textMassPHOS = new TLatex(columnsLegendMass[0],rowsLegendMass[1],"PHOS");
	SetStyleTLatex( textMassPHOS, textSizeLeftColumnMass,4);
	textMassPHOS->SetTextFont(42);
	textMassPHOS->Draw();


	
	//****************** second Column *************************************************
	TLatex *textMassData2 = new TLatex(columnsLegendMass[1],rowsLegendMass[0] ,"data");
	SetStyleTLatex( textMassData2, textSizeTopRowMass ,4);
	textMassData2->SetTextFont(42);
		textMassData2->Draw();
	TLatex *textMassMC2 = new TLatex(columnsLegendMass[2]
	,rowsLegendMass[0],"MC");
	SetStyleTLatex( textMassMC2, textSizeTopRowMass,4);
	textMassMC2->SetTextFont(42);
	textMassMC2->Draw();
	
	TMarker* markerPCMPi0Mass = CreateMarkerFromHisto(PCMPi0Mass,columnsLegendMassAbs[1] ,rowsLegendMassAbs[3] ,scaleMarkerMass);
	markerPCMPi0Mass->Draw();
 	TMarker* markerEMCALEMCALPi0Mass = CreateMarkerFromHisto(EMCalPi0Mass,columnsLegendMassAbs[1],rowsLegendMassAbs[2],scaleMarkerMass);
 	markerEMCALEMCALPi0Mass->Draw();
	TMarker* markerDalitzPi0Mass = CreateMarkerFromHisto(DalitzPi0Mass,columnsLegendMassAbs[1],rowsLegendMassAbs[4],scaleMarkerMass);
	markerDalitzPi0Mass->Draw();
	TMarker* markerPHOSPi0Mass = CreateMarkerFromHisto(PHOSPi0Mass,columnsLegendMassAbs[1],rowsLegendMassAbs[1],scaleMarkerMass);
	markerPHOSPi0Mass->Draw();

	
	TMarker* markerPCMPi0MassMC = CreateMarkerFromHisto(PCMPi0MassMC,columnsLegendMassAbs[2],rowsLegendMassAbs[3],scaleMarkerMass);
	markerPCMPi0MassMC->Draw();
 	TMarker* markerEMCALEMCALPi0MassMC = CreateMarkerFromHisto(EMCalPi0MassMC,columnsLegendMassAbs[2] ,rowsLegendMassAbs[2] ,scaleMarkerMass);
 	markerEMCALEMCALPi0MassMC->Draw();
	TMarker* markerDalitzPi0MassMC = CreateMarkerFromHisto(DalitzPi0MassMC,columnsLegendMassAbs[2] ,rowsLegendMassAbs[4] ,scaleMarkerMass);
	markerDalitzPi0MassMC->Draw();
	TMarker* markerPHOSPi0MassMC = CreateMarkerFromHisto(PHOSPi0MassMC,columnsLegendMassAbs[2] ,rowsLegendMassAbs[1] ,scaleMarkerMass);
	markerPHOSPi0MassMC->Draw();


 	
  TLatex * lt12 = new TLatex(2.8,146.5,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV") ;
  lt12->SetTextColor(kBlack) ;
  lt12->SetTextSize(TextSize) ;
  lt12->SetTextFont(Font) ;
  // lt12->DrawLatex(2.8,1.95,"NSD #pi^{0}");
  
  lt12->Draw() ;



  c12->Update();
  c12->Print(Form("%s/Pi0Mass.%s",outputDir.Data(),suffix.Data()));

  //MassEta
  TCanvas* c13 = new TCanvas("c13","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c13, CMarginL, CMarginR,CMarginT ,CMarginB);
  c13->SetLogx();
  TH2F * hist13;
  hist13 = new TH2F("hist13","hist13",1000,0.3,30.,1000,495.,595.);
  SetStyleHistoTH2ForGraphs(hist13, "#it{p}_{T} (GeV/#it{c})","peak position (MeV/#it{c}^{2})",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY,TitleOffsetX,TitleOffsetY, 512, 508);
  hist13->GetXaxis()->SetLabelOffset(LabelOffsetLog);
  hist13->DrawCopy(); 

  TLine* lineEta=new TLine(0.3,mesonMassExpectEta*1000,30,mesonMassExpectEta*1000);
  lineEta->SetLineColor(kGray);
  lineEta->Draw("same");

 	 DrawGammaSetMarker(PCMEtaMassMC, markerStylePCMMC, markerSizePCM, colorPCMMC , colorPCMMC);
	 PCMEtaMassMC->DrawCopy("same,z,p,x0");     
         DrawGammaSetMarker(PCMEtaMass, markerStylePCM,markerSizePCM , colorPCM, colorPCM);
	 PCMEtaMass->DrawCopy("same,z,p,x0"); 

	 DrawGammaSetMarker(EMCalEtaMassMC, markerStyleEMCalMC,markerSizeEMCal , colorEMCalMC , colorEMCalMC);
	 EMCalEtaMassMC->DrawCopy("same,z,p,x0");        
	 DrawGammaSetMarker(EMCalEtaMass, markerStyleEMCal,markerSizeEMCal , colorEMCal, colorEMCal);
	 EMCalEtaMass->DrawCopy("same,z,p,x0"); 


	 Double_t rowsLegendMassAbsEta[7]       = {0.2,582.8,576.9,21.2,19.7, 18.2, 18.2};

	//****************** first Column **************************************************
	TLatex *textMassPCMEta = new TLatex(columnsLegendMass[0],rowsLegendMass[2],"PCM");
	SetStyleTLatex( textMassPCMEta, textSizeLeftColumnMass,4);
	textMassPCMEta->SetTextFont(42);
	textMassPCMEta->Draw();

	
 	TLatex *textMassEMCALEMCALEta = new TLatex(columnsLegendMass[0],rowsLegendMass[1],"EMCal");
 	SetStyleTLatex( textMassEMCALEMCALEta, textSizeLeftColumnMass,4);
 	textMassEMCALEMCALEta->SetTextFont(42);
 	textMassEMCALEMCALEta->Draw();
	
	//****************** second Column *************************************************
	TLatex *textMassData2Eta = new TLatex(columnsLegendMass[1],rowsLegendMass[0] ,"data");
	SetStyleTLatex( textMassData2Eta, textSizeTopRowMass ,4);
	textMassData2Eta->SetTextFont(42);
		textMassData2Eta->Draw();
	TLatex *textMassMC2Eta = new TLatex(columnsLegendMass[2],rowsLegendMass[0],"MC");
	SetStyleTLatex( textMassMC2Eta, textSizeTopRowMass,4);
	textMassMC2Eta->SetTextFont(42);
	textMassMC2Eta->Draw();
	
	TMarker* markerPCMEtaMassEta = CreateMarkerFromHisto(PCMEtaMass,columnsLegendMassAbs[1] ,rowsLegendMassAbsEta[2] ,scaleMarkerMass);
	markerPCMEtaMassEta->Draw();
 	TMarker* markerEMCALEMCALEtaMassEta= CreateMarkerFromHisto(EMCalEtaMass,columnsLegendMassAbs[1],rowsLegendMassAbsEta[1],scaleMarkerMass);
 	markerEMCALEMCALEtaMassEta->Draw();
	
	TMarker* markerPCMEtaMassMCEta = CreateMarkerFromHisto(PCMEtaMassMC,columnsLegendMassAbs[2],rowsLegendMassAbsEta[2],scaleMarkerMass);
	markerPCMEtaMassMCEta->Draw();
 	TMarker* markerEMCALEMCALEtaMassMCEta = CreateMarkerFromHisto(EMCalEtaMassMC,columnsLegendMassAbs[2] ,rowsLegendMassAbsEta[1] ,scaleMarkerMass);
 	markerEMCALEMCALEtaMassMCEta->Draw();



 	
  TLatex * lt13 = new TLatex(2.8,584.5,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV") ;
  lt13->SetTextColor(kBlack);
  lt13->SetTextSize(TextSize);
  lt13->SetTextFont(Font);
   lt13->DrawLatex(2.8,1.95,"NSD #pi^{0}");
  
  lt13->Draw("same") ;


  c13->Update();
  c13->Print(Form("%s/EtaMass.%s",outputDir.Data(),suffix.Data()));
  
  //Definition of the plot
 	
  
  Double_t rowsTheoryLegendPrelAndFinal[8] = {0.9,0.77,0.66,0.50,0.39,0.28,0.17,0.05};
  
  //*************** Size factors ********************
	
  Double_t scaleMarkerTheory    = 1.0;
  Double_t scaleLineWidthTheory = 1.;

  //*************** Label sizes *********************
	
  Double_t textSizeLeftLabelsPrelAndFinal 	= 0.11;
  Double_t textSizeTopLablesPrelAndFinal 	= 0.115;
  Double_t textSizeTopLowerLablesPrelAndFinal  = 0.11;
  //*************** Offsets *************************
	
  Double_t offsetTheoryLegendPrelAndFinalMarkerX = 0.09;
  Double_t offsetTheoryLegendPrelAndFinalMarkerY = 0.03;
  Double_t offsetTheoryLegendPrelAndFinalBox     = 0.05;
  Double_t offsetTheoryLegendPrelAndFinalLine    = 0.06;
		
  Double_t rowsTheoryLegendPrelAndFinalOnlyRatioPi0[8] = {0.9,0.77,0.66,0.50,0.39,0.28,0.17,0.05}; // {0.88,0.75,0.62,0.49,0.36,0.25,0.12}; //{0.9,0.79,0.69,0.54,0.41,0.26,0.12,0.0};
  Double_t columnsTheoryLegendPrelAndFinalOnlyRatio[4] = {0.,0.27,0.61,0.89};
	
  Double_t textSizeLeftLabelsPrelAndFinalRatio = 0.105;
  Double_t textSizeTopLablesPrelAndFinalRatio = 0.11;
  Double_t textSizeTopLowerLablesPrelAndFinalRatio = 0.105;
	
	
	//*************** Column def ***********************
	
  Double_t columnsTheoryLegendPrelAndFinal[4] = {0.,0.27,0.61,0.89};//{0.,0.23,0.46,0.74};
	
    	
  TLatex *textTsallisFitALLEnergies = new TLatex(columnsTheoryLegendPrelAndFinal[0],rowsTheoryLegendPrelAndFinal[3],"Tsallis fit");
  SetStyleTLatex( textTsallisFitALLEnergies, textSizeLeftLabelsPrelAndFinal,4);
  
  TLatex *textEPOS3 = new TLatex(columnsTheoryLegendPrelAndFinal[0],rowsTheoryLegendPrelAndFinal[4],"EPOS3");
  SetStyleTLatex( textEPOS3, textSizeLeftLabelsPrelAndFinal,4);
  TLatex *textCGC = new TLatex(columnsTheoryLegendPrelAndFinal[0],rowsTheoryLegendPrelAndFinal[5],"CGC");
  SetStyleTLatex( textCGC, textSizeLeftLabelsPrelAndFinal,4);
 
  //*************** second Column **********************************************************
  TLatex *textPi0pPb5023GeV = new TLatex(columnsTheoryLegendPrelAndFinal[1],rowsTheoryLegendPrelAndFinal[0],"#pi^{0}, #sqrt{#it{s_{NN}}} = 5.02 TeV");
  SetStyleTLatex( textPi0pPb5023GeV, textSizeTopLablesPrelAndFinal,4);
	
  TLatex *textPi0pPb5023GeVsys = new TLatex(columnsTheoryLegendPrelAndFinal[1],rowsTheoryLegendPrelAndFinal[1],"syst. + stat.");
  SetStyleTLatex( textPi0pPb5023GeVsys, textSizeTopLowerLablesPrelAndFinal,4);
	
  TBox* boxCombinedPi0pPb5023GeV = CreateBoxFromGraph(CombPi0Syst,columnsTheoryLegendPrelAndFinal[1]+offsetTheoryLegendPrelAndFinalMarkerX-offsetTheoryLegendPrelAndFinalLine, rowsTheoryLegendPrelAndFinal[2]+ offsetTheoryLegendPrelAndFinalMarkerY- offsetTheoryLegendPrelAndFinalBox, columnsTheoryLegendPrelAndFinal[1]+offsetTheoryLegendPrelAndFinalMarkerX+offsetTheoryLegendPrelAndFinalLine, rowsTheoryLegendPrelAndFinal[2]+ offsetTheoryLegendPrelAndFinalMarkerY+offsetTheoryLegendPrelAndFinalBox);	
  TMarker* markerCombinedPi0pPb5023GeV = CreateMarkerFromGraph(CombPi0Syst,columnsTheoryLegendPrelAndFinal[1]+offsetTheoryLegendPrelAndFinalMarkerX,rowsTheoryLegendPrelAndFinal[2]+offsetTheoryLegendPrelAndFinalMarkerY ,scaleMarkerTheory);
  
  TLine * lineTsallisFitPi0pPb5023GeV 	= CreateLineFromFit(FitCombPi0, 				columnsTheoryLegendPrelAndFinal[1]+ offsetTheoryLegendPrelAndFinalMarkerX- offsetTheoryLegendPrelAndFinalLine, rowsTheoryLegendPrelAndFinal[3]+offsetTheoryLegendPrelAndFinalMarkerY, columnsTheoryLegendPrelAndFinal[1]+ offsetTheoryLegendPrelAndFinalMarkerX+ offsetTheoryLegendPrelAndFinalLine,rowsTheoryLegendPrelAndFinal[3]+offsetTheoryLegendPrelAndFinalMarkerY, scaleLineWidthTheory);
  TLine * lineEPOS3Pi0pPb5023GeV 	= CreateLineFromGraph(graphPi0EPOS, 				columnsTheoryLegendPrelAndFinal[1]+ offsetTheoryLegendPrelAndFinalMarkerX- offsetTheoryLegendPrelAndFinalLine, rowsTheoryLegendPrelAndFinal[4]+offsetTheoryLegendPrelAndFinalMarkerY, columnsTheoryLegendPrelAndFinal[1]+ offsetTheoryLegendPrelAndFinalMarkerX+ offsetTheoryLegendPrelAndFinalLine,rowsTheoryLegendPrelAndFinal[4]+offsetTheoryLegendPrelAndFinalMarkerY, scaleLineWidthTheory); 
  TLine * lineCGCPi0pPb5023GeV 	        = CreateLineFromGraph(graphAsymmErrCGCTheoryPi0y0pA5020,	columnsTheoryLegendPrelAndFinal[1]+ offsetTheoryLegendPrelAndFinalMarkerX- offsetTheoryLegendPrelAndFinalLine, rowsTheoryLegendPrelAndFinal[5]+offsetTheoryLegendPrelAndFinalMarkerY, columnsTheoryLegendPrelAndFinal[1]+ offsetTheoryLegendPrelAndFinalMarkerX+ offsetTheoryLegendPrelAndFinalLine,rowsTheoryLegendPrelAndFinal[5]+offsetTheoryLegendPrelAndFinalMarkerY, scaleLineWidthTheory); 
//   TLine * lineNLOPi07TeVMuTwoALLEnergies 	= CreateLineFromGraph(graphNLOMuTwoPi07TeV, 			columnsTheoryLegendPrelAndFinal[1]+ offsetTheoryLegendPrelAndFinalMarkerX- offsetTheoryLegendPrelAndFinalLine, rowsTheoryLegendPrelAndFinal[6]+offsetTheoryLegendPrelAndFinalMarkerY, columnsTheoryLegendPrelAndFinal[1]+ offsetTheoryLegendPrelAndFinalMarkerX+ offsetTheoryLegendPrelAndFinalLine,rowsTheoryLegendPrelAndFinal[6]+offsetTheoryLegendPrelAndFinalMarkerY, scaleLineWidthTheory); 
//   TLine * lineNLODSS14MuOnePi07TeVEnergies      = CreateLineFromGraph(graphNLODSS14MuOne7TeV,		        columnsTheoryLegendPrelAndFinal[1]+ offsetTheoryLegendPrelAndFinalMarkerX- offsetTheoryLegendPrelAndFinalLine, rowsTheoryLegendPrelAndFinal[5]+offsetTheoryLegendPrelAndFinalMarkerY, columnsTheoryLegendPrelAndFinal[1]+ offsetTheoryLegendPrelAndFinalMarkerX+ offsetTheoryLegendPrelAndFinalLine,rowsTheoryLegendPrelAndFinal[5]+offsetTheoryLegendPrelAndFinalMarkerY, scaleLineWidthTheory); 
// 	

  //*************** third Column **********************************************************
  TLatex *textEtapPb5023GeV = new TLatex(columnsTheoryLegendPrelAndFinal[2],rowsTheoryLegendPrelAndFinal[0],"#eta, #sqrt{#it{s_{NN}}} = 5.02 TeV");
  SetStyleTLatex( textEtapPb5023GeV, textSizeTopLablesPrelAndFinal,4);
	
  TLatex *textEtapPb5023GeVsys = new TLatex(columnsTheoryLegendPrelAndFinal[2],rowsTheoryLegendPrelAndFinal[1],"syst. + stat.");
  SetStyleTLatex( textEtapPb5023GeVsys, textSizeTopLowerLablesPrelAndFinal,4);
	
  TBox* boxCombinedEtapPb5023GeV       = CreateBoxFromGraph(CombEtaSyst,columnsTheoryLegendPrelAndFinal[2]+offsetTheoryLegendPrelAndFinalMarkerX-offsetTheoryLegendPrelAndFinalLine, rowsTheoryLegendPrelAndFinal[2]+ offsetTheoryLegendPrelAndFinalMarkerY- offsetTheoryLegendPrelAndFinalBox, columnsTheoryLegendPrelAndFinal[2]+offsetTheoryLegendPrelAndFinalMarkerX+offsetTheoryLegendPrelAndFinalLine, rowsTheoryLegendPrelAndFinal[2]+ offsetTheoryLegendPrelAndFinalMarkerY+offsetTheoryLegendPrelAndFinalBox);
  TMarker* markerCombinedEtapPb5023GeV = CreateMarkerFromGraph(CombEtaSyst,columnsTheoryLegendPrelAndFinal[2]+offsetTheoryLegendPrelAndFinalMarkerX,rowsTheoryLegendPrelAndFinal[2]+offsetTheoryLegendPrelAndFinalMarkerY ,scaleMarkerTheory);
	
  TLine * lineTsallisFitEtapPb5023GeV 	=	CreateLineFromFit(FitCombEta, 		columnsTheoryLegendPrelAndFinal[2]+offsetTheoryLegendPrelAndFinalMarkerX-offsetTheoryLegendPrelAndFinalLine, rowsTheoryLegendPrelAndFinal[3]+offsetTheoryLegendPrelAndFinalMarkerY, columnsTheoryLegendPrelAndFinal[2]+ offsetTheoryLegendPrelAndFinalMarkerX+ offsetTheoryLegendPrelAndFinalLine,rowsTheoryLegendPrelAndFinal[3]+offsetTheoryLegendPrelAndFinalMarkerY, scaleLineWidthTheory);
  TLine * lineEPOS3EtapPb5023GeV 	= 	CreateLineFromGraph(graphEtaEPOS,  			columnsTheoryLegendPrelAndFinal[2]+offsetTheoryLegendPrelAndFinalMarkerX-offsetTheoryLegendPrelAndFinalLine, rowsTheoryLegendPrelAndFinal[4]+offsetTheoryLegendPrelAndFinalMarkerY, columnsTheoryLegendPrelAndFinal[2]+ offsetTheoryLegendPrelAndFinalMarkerX+ offsetTheoryLegendPrelAndFinalLine,rowsTheoryLegendPrelAndFinal[4]+offsetTheoryLegendPrelAndFinalMarkerY, scaleLineWidthTheory); 
  
  /////////////////////////////////////////////////////////////////////////
  //Scaling factors
  
  Double_t scaleFacPi0Yield = 1e1;
  Double_t scaleFacEtaYield = 1e-1;
  
  TGraphAsymmErrors* CombPi0ScaledSyst = (TGraphAsymmErrors*)ScaleGraph(CombPi0Syst,scaleFacPi0Yield);
  TGraphAsymmErrors* CombPi0ScaledStat = (TGraphAsymmErrors*)ScaleGraph(CombPi0Stat,scaleFacPi0Yield);
  TH1F* FitCombPi0Scaled = (TH1F*)FitCombPi0->GetHistogram();
  FitCombPi0Scaled->Scale(scaleFacPi0Yield);
  
  TGraphAsymmErrors* CombEtaScaledSyst = (TGraphAsymmErrors*)ScaleGraph(CombEtaSyst,scaleFacEtaYield);
  TGraphAsymmErrors* CombEtaScaledStat = (TGraphAsymmErrors*)ScaleGraph(CombEtaStat,scaleFacEtaYield);
  TH1F* FitCombEtaScaled = (TH1F*)FitCombEta->GetHistogram();
  FitCombEtaScaled->Scale(scaleFacEtaYield);
  
  
  TGraphAsymmErrors* graphPi0EPOSScaled 	= (TGraphAsymmErrors*)ScaleGraph(graphPi0EPOS,scaleFacPi0Yield);
  TGraphAsymmErrors* graphEtaEPOSScaled 	= (TGraphAsymmErrors*)ScaleGraph(graphEtaEPOS,scaleFacEtaYield);
  
  TGraphAsymmErrors* graphAsymmErrCGCTheoryPi0y0pA5020Scaled = (TGraphAsymmErrors*)ScaleGraph(graphAsymmErrCGCTheoryPi0y0pA5020,1e1);
  
  TGraphAsymmErrors* graphErrMcGillTheoryPion_p_hydroScaled = (TGraphAsymmErrors*)ScaleGraph(graphErrMcGillTheoryPion_p_hydro,1e1);
  TGraphAsymmErrors* graphErrMcGillTheoryEta_p_hydroScaled  = (TGraphAsymmErrors*)ScaleGraph(graphErrMcGillTheoryEta_p_hydro,scaleFacEtaYield);
   
  
  TH1F* histoHIJINGPi0Scaled = (TH1F*) histoHIJINGPi0->Clone("histoHIJINGPi0Scaled");
  histoHIJINGPi0Scaled->Scale(scaleFacPi0Yield);
  histoHIJINGPi0Scaled->GetXaxis()->SetRangeUser(0.3, 20);
   
   
  TH1F* histoDPMJetPi0Scaled = (TH1F*) histoDPMJetPi0->Clone("histoDPMJetPi0Scaled");
  histoDPMJetPi0Scaled->Scale(scaleFacPi0Yield);
  histoDPMJetPi0Scaled->GetXaxis()->SetRangeUser(0.3, 20);
  
  
  TH1F* histoHIJINGEtaScaled = (TH1F*) histoHIJINGEta->Clone("histoHIJINGEtaScaled");
  histoHIJINGEtaScaled->Scale(scaleFacEtaYield);
  histoHIJINGEtaScaled->GetXaxis()->SetRangeUser(0.3, 20);
   
   
  TH1F* histoDPMJetEtaScaled = (TH1F*) histoDPMJetEta->Clone("histoDPMJetEtaScaled");
  histoDPMJetEtaScaled->Scale(scaleFacEtaYield);
  histoDPMJetEtaScaled->GetXaxis()->SetRangeUser(0.3, 20);
  
  
 
 
  DrawGammaSetMarkerTGraphAsym(CombPi0ScaledSyst,markerStylePi0,1,colorCombYieldPi0 ,colorCombYieldPi0, 1, kTRUE);
  //CombPi0ScaledSyst->Draw("E2,same"); 
  DrawGammaSetMarkerTGraphAsym(CombPi0ScaledStat,markerStylePi0,1,colorCombYieldPi0, colorCombYieldPi0 ); 
  //CombPi0ScaledStat->Draw("pz,same"); 	
  
  DrawGammaSetMarkerTGraphAsym(CombEtaScaledSyst,markerStyleEta,1,colorCombYieldEta ,colorCombYieldEta, 1, kTRUE);
  //CombPi0ScaledSyst->Draw("E2,same"); 
  DrawGammaSetMarkerTGraphAsym(CombEtaScaledStat,markerStyleEta,1,colorCombYieldEta, colorCombYieldEta ); 
  //CombPi0ScaledStat->Draw("pz,same"); 	
  
  
  
  TCanvas* canvasInvYieldTSallisTheoryOnlyPi0AndEtaSpectra = new TCanvas("canvasInvYieldTSallisTheoryOnlyPi0AndEtaSpectra","",200,10,1000,1200);  // gives the page size
  DrawGammaCanvasSettings( canvasInvYieldTSallisTheoryOnlyPi0AndEtaSpectra,  0.3, 0.02, 0.02, 0.16);
	
  TPad* padComparisonInvYieldTSallisTheoryOnlyPi0AndEtaSpectra = new TPad("padComparisonInvYieldTSallisTheoryOnlyPi0AndEtaSpectra", "", 0., 0., 1., 1.,-1, -1, -2);
  DrawGammaPadSettings( padComparisonInvYieldTSallisTheoryOnlyPi0AndEtaSpectra, 0.15, 0.02, 0.02, 0.06);
  padComparisonInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->Draw();
	
  	
  padComparisonInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->cd();
  padComparisonInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->SetLogy();		
  padComparisonInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->SetLogx();		
	
  TH2F * histo2DInvYieldTSallisTheoryOnlyPi0AndEtaSpectra;
  histo2DInvYieldTSallisTheoryOnlyPi0AndEtaSpectra = new TH2F("histo2DInvYieldTSallisTheoryOnlyPi0AndEtaSpectra","histo2DInvYieldTSallisTheoryOnlyPi0AndEtaSpectra",1000,0.2,30.,1000,4.2e-10,20 );
  SetStyleHistoTH2ForGraphs(histo2DInvYieldTSallisTheoryOnlyPi0AndEtaSpectra, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2} ",0.035,0.035,0.035,0.035, 0.8,1.9, 512, 510);
  histo2DInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->GetXaxis()->SetLabelOffset(-0.009);
  histo2DInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->DrawCopy(); 
	
 
  
  FitCombPi0->SetLineWidth(2);  
  FitCombEtaScaled->SetLineWidth(2);  
  FitCombPi0->SetLineColor(kBlack);  
  FitCombEtaScaled->SetLineColor(kBlack);
  FitCombPi0->SetLineStyle(3);  
  FitCombEtaScaled->SetLineStyle(3);
  
   graphAsymmErrIlkkapPb5020_pi0_ct14_epps16_dss14_scale_err->SetLineWidth(3);
   graphAsymmErrIlkkapPb5020_pi0_ct14_epps16_dss14_scale_err->SetLineStyle(styleLineIlkka);
   graphAsymmErrIlkkapPb5020_pi0_ct14_epps16_dss14_scale_err->SetLineColor(colorIlkkaPi0);
   graphAsymmErrIlkkapPb5020_pi0_ct14_epps16_dss14_scale_err->SetFillColor(colorIlkkaPi0);
   
   //raphAsymmErrRpPb5020_pi0_ct14_epps16_dss14->SetFillColor(kYellow-7);
  //graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14->SetFillStyle(1001);
   graphAsymmErrIlkkapPb5020_pi0_ct14_epps16_dss14_scale_err->Draw("c3same");
   
   
  FitCombPi0->Draw("lsame"); 
  FitCombEtaScaled->Draw("lsame");
 
 
  if (EPOS){
      while(graphPi0EPOS->GetX()[0] < 0.3)
          graphPi0EPOS->RemovePoint(0); 
    graphPi0EPOS->SetLineWidth(3);
    graphPi0EPOS->SetLineStyle(styleLineEPOS3);
    graphPi0EPOS->SetLineColor(colorEPOSPi0);
    graphPi0EPOS->SetFillColor(colorEPOSPi0);
   
     while(graphEtaEPOSScaled->GetX()[0] < 0.7)
          graphEtaEPOSScaled->RemovePoint(0); 
    graphEtaEPOSScaled->SetLineWidth(3);
    graphEtaEPOSScaled->SetLineStyle(styleLineEPOS3);
    graphEtaEPOSScaled->SetLineColor(colorEPOSEta);
    graphEtaEPOSScaled->SetFillColor(colorEPOSEta);
    
    
    graphPi0EPOS->Draw("C3same"); //C
    graphEtaEPOSScaled->Draw("C3same");
    
    while(graphErrMcGillTheoryPion_p_hydro->GetX()[0] < 0.3)
          graphErrMcGillTheoryPion_p_hydro->RemovePoint(0);   
    graphErrMcGillTheoryPion_p_hydro->SetLineWidth(3);
    graphErrMcGillTheoryPion_p_hydro->SetLineStyle(styleLineMcGill);
    graphErrMcGillTheoryPion_p_hydro->SetLineColor(colorMcGillPi0);
    graphErrMcGillTheoryPion_p_hydro->SetFillColor(colorMcGillPi0);
    graphErrMcGillTheoryPion_p_hydro->Draw("C3same");
    
   while(graphErrMcGillTheoryEta_p_hydroScaled->GetX()[0] < 0.7)
          graphErrMcGillTheoryEta_p_hydroScaled->RemovePoint(0);   
    graphErrMcGillTheoryEta_p_hydroScaled->SetLineWidth(3);
    graphErrMcGillTheoryEta_p_hydroScaled->SetLineStyle(styleLineMcGill);
    graphErrMcGillTheoryEta_p_hydroScaled->SetLineColor(colorMcGillEta);
    graphErrMcGillTheoryEta_p_hydroScaled->SetFillColor(colorMcGillEta);
    graphErrMcGillTheoryEta_p_hydroScaled->Draw("c3same");
        
    
  }
  
  
  if( CGCPi0 ) {
    while(graphAsymmErrCGCTheoryPi0y0pA5020->GetX()[0] < 0.3)
          graphAsymmErrCGCTheoryPi0y0pA5020->RemovePoint(0);  
    graphAsymmErrCGCTheoryPi0y0pA5020->SetLineWidth(3);
    graphAsymmErrCGCTheoryPi0y0pA5020->SetLineStyle(styleLineCGC);
    graphAsymmErrCGCTheoryPi0y0pA5020->SetLineColor(kOrange+1);
    graphAsymmErrCGCTheoryPi0y0pA5020->SetFillColor(kOrange+1);
    graphAsymmErrCGCTheoryPi0y0pA5020->Draw("C3same");    
  }
  
  if( HIJINGPi0 ){
    SetStyleHisto(histoHIJINGPi0, 3, styleLineHIJING, colorHIJINGPi0);          
    histoHIJINGPi0->GetXaxis()->SetRangeUser(0.3,20);
    histoHIJINGPi0->Draw("same,hist,l");
    
    SetStyleHisto(histoHIJINGEtaScaled, 3, styleLineHIJING, colorHIJINGEta);    
     histoHIJINGEtaScaled->GetXaxis()->SetRangeUser(0.7,20);
    histoHIJINGEtaScaled->Draw("same,hist,l");
    
    
    
  }
  
  if( DPMJetPi0 ){
    
   SetStyleHisto(histoDPMJetPi0, 3, styleLineDPMJet, colorDPMJetPi0);      
   histoDPMJetPi0->GetXaxis()->SetRangeUser(0.3,20);
   histoDPMJetPi0->Draw("same,hist,l"); 
   
   SetStyleHisto(histoDPMJetEtaScaled, 3, styleLineDPMJet, colorDPMJetEta);     
   histoDPMJetEtaScaled->GetXaxis()->SetRangeUser(0.7,20);
   histoDPMJetEtaScaled->Draw("same,hist,l"); 
   
  }
  CombPi0Syst->Draw("E2,same"); 
  CombPi0Stat->Draw("pz,same"); 
  CombEtaScaledSyst->Draw("E2,same");
  CombEtaScaledStat->Draw("pz,same"); 

  if (mT){
      EtaSpectrum_mT_pPb->Draw("same"); 
  }
 
  
  TLatex * latexInvYieldTSallisTheoryOnlyPi0AndEtaSpectra = new TLatex(2.5,4,"ALICE") ;
  latexInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->SetTextColor(kBlack) ;
  latexInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->SetTextSize(TextSize-0.013) ;
  latexInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->SetTextFont(Font) ;
  latexInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->DrawLatex(2.5,2,"p-Pb, NSD, #sqrt{#it{s}_{NN}} = 5.02 TeV");
 
  latexInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->Draw() ;
  
  

  TLegend* legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra;
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra   = new TLegend(0.20,0.35,0.4,0.58); 
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->SetFillColor(0);
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->SetLineColor(0);
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->SetTextFont(Font);
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->SetTextSize(TextSize-0.015);
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->AddEntry(CombPi0Syst,"#pi^{0}","pef");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->AddEntry(graphPi0EPOS,"EPOS3","l");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->AddEntry(graphErrMcGillTheoryPion_p_hydro,"VISHNU","l");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->AddEntry(histoDPMJetPi0,"DPMJet","l");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->AddEntry(histoHIJINGPi0,"HIJING","l");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->AddEntry(graphAsymmErrCGCTheoryPi0y0pA5020,"CGC MV^{#gamma}","l");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->AddEntry(graphAsymmErrIlkkapPb5020_pi0_ct14_epps16_dss14_scale_err,"NLO pQCD","l");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->AddEntry(FitCombPi0,"Tsallis fit","l");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->Draw();
  
  TLegend* legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2;
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2   = new TLegend(0.20,0.12,0.4,0.30); 
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2->SetFillColor(0);
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2->SetLineColor(0);
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2->SetTextFont(Font);
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2->SetTextSize(TextSize-0.015); 
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2->AddEntry(CombEtaScaledSyst,"#eta #times 10^{-1}","pef");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2->AddEntry(graphEtaEPOSScaled,"EPOS3","l");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2->AddEntry(graphErrMcGillTheoryEta_p_hydroScaled,"VISHNU","l");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2->AddEntry(histoDPMJetEtaScaled,"DPMJet","l");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2->AddEntry(histoHIJINGEtaScaled,"HIJING","l");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2->AddEntry(FitCombPi0,"Tsallis fit","l");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2->Draw();
  
  
   padComparisonInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->Update();

   canvasInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->Update();
	
   if(EPOS) canvasInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->Print(Form("%s/MesonYields_EPOSv2.%s",outputDir.Data(),suffix.Data()));
   else if(mT) canvasInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->Print(Form("%s/MesonYields_mTv2.%s",outputDir.Data(),suffix.Data()));
   else canvasInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->Print(Form("%s/MesonYieldsv2.%s",outputDir.Data(),suffix.Data()));
  
   
   
   
   
   
   
   ///////////////////////////////////////////////////////////////
   
  TCanvas* canvasInvYieldTSallisTheoryOnlyPi0Spectrum = new TCanvas("canvasInvYieldTSallisTheoryOnlyPi0Spectrum","",200,10,1000,1200);  // gives the page size
  DrawGammaCanvasSettings( canvasInvYieldTSallisTheoryOnlyPi0Spectrum,  0.3, 0.02, 0.02, 0.16);
	
  TPad* padComparisonInvYieldTSallisTheoryOnlyPi0Spectrum = new TPad("padComparisonInvYieldTSallisTheoryOnlyPi0Spectrum", "", 0., 0., 1., 1.,-1, -1, -2);
  DrawGammaPadSettings( padComparisonInvYieldTSallisTheoryOnlyPi0Spectrum, 0.15, 0.02, 0.02, 0.09);
  padComparisonInvYieldTSallisTheoryOnlyPi0Spectrum->Draw();
	
  	
  padComparisonInvYieldTSallisTheoryOnlyPi0Spectrum->cd();
  padComparisonInvYieldTSallisTheoryOnlyPi0Spectrum->SetLogy();		
  padComparisonInvYieldTSallisTheoryOnlyPi0Spectrum->SetLogx();		
	
  TH2F * histo2DInvYieldTSallisTheoryOnlyPi0Spectrum;
  histo2DInvYieldTSallisTheoryOnlyPi0Spectrum = new TH2F("histo2DInvYieldTSallisTheoryOnlyPi0Spectrum","histo2DInvYieldTSallisTheoryOnlyPi0Spectrum",1000,0.2,30.,1000,1.2e-8,30 );
  SetStyleHistoTH2ForGraphs(histo2DInvYieldTSallisTheoryOnlyPi0Spectrum, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2} ",0.035,0.035,0.035,0.035, 1.0,1.8, 512, 510);
  histo2DInvYieldTSallisTheoryOnlyPi0Spectrum->DrawCopy(); 
	
 
  
  FitCombPi0->SetLineWidth(2);  
  FitCombEta->SetLineWidth(2);  
  FitCombPi0->SetLineColor(kBlack);  
  FitCombEta->SetLineColor(kBlack);
  FitCombPi0->SetLineStyle(4);  
  FitCombEta->SetLineStyle(4);
  
   
   
  FitCombPi0->Draw("same"); 
  //FitCombEta->Draw("same");
 
 
  if (EPOS){
    graphPi0EPOS->SetLineWidth(3);
    //graphEtaEPOS->SetLineWidth(3);
    //graphEtaEPOS->SetLineColor(kGreen+3);
    //graphEtaEPOS->SetFillColor(kGreen+3);
    
    graphPi0EPOS->Draw("C3same"); //C
    
    //graphEtaEPOS->Draw("C3same");
  }
  
  if( CGCPi0 ) {
    graphAsymmErrCGCTheoryPi0y0pA5020->SetLineWidth(3);
    graphAsymmErrCGCTheoryPi0y0pA5020->SetLineColor(kOrange+1);
    graphAsymmErrCGCTheoryPi0y0pA5020->SetFillColor(kOrange+1);
    graphAsymmErrCGCTheoryPi0y0pA5020->Draw("C3same");
    
  }
  
  if( HIJINGPi0 ){
   // DrawGammaSetMarker(histoHIJINGPi0, 24, 1.5, colorHIJING , colorHIJING);
    SetStyleHisto(histoHIJINGPi0, 3, styleLineHIJING, colorHIJING);          
  
    histoHIJINGPi0->SetLineWidth(3);
    histoHIJINGPi0->Draw("same,hist,l");
    
  }
  
  if( DPMJetPi0 ){
    
   SetStyleHisto(histoDPMJetPi0, 3, styleLineDPMJet, colorDPMJet);          
   histoDPMJetPi0->Draw("same,hist,l"); 
   
  }
  CombPi0Syst->Draw("E2,same"); 
  CombPi0Stat->Draw("pz,same"); 
  //CombEtaSyst->Draw("E2,same");
  //CombEtaStat->Draw("pz,same"); 

 
 
  
  TLatex * latexInvYieldTSallisTheoryOnlyPi0Spectrum = new TLatex(3,9,"ALICE") ;
  latexInvYieldTSallisTheoryOnlyPi0Spectrum->SetTextColor(kBlack) ;
  latexInvYieldTSallisTheoryOnlyPi0Spectrum->SetTextSize(TextSize-0.015) ;
  latexInvYieldTSallisTheoryOnlyPi0Spectrum->SetTextFont(Font) ;
  latexInvYieldTSallisTheoryOnlyPi0Spectrum->DrawLatex(3,4,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
 
  latexInvYieldTSallisTheoryOnlyPi0Spectrum->Draw() ;
  
  

  TLegend* legendInvYieldTSallisTheoryOnlyPi0Spectrum;
  legendInvYieldTSallisTheoryOnlyPi0Spectrum   = new TLegend(0.20,0.12,0.6,0.33); 
  legendInvYieldTSallisTheoryOnlyPi0Spectrum->SetFillColor(0);
  legendInvYieldTSallisTheoryOnlyPi0Spectrum->SetLineColor(0);
  legendInvYieldTSallisTheoryOnlyPi0Spectrum->SetTextFont(Font);
  legendInvYieldTSallisTheoryOnlyPi0Spectrum->SetTextSize(TextSize-0.015);
  legendInvYieldTSallisTheoryOnlyPi0Spectrum->AddEntry(CombPi0Syst,"NSD #pi^{0}","pef");
  legendInvYieldTSallisTheoryOnlyPi0Spectrum->AddEntry(graphPi0EPOS,"#pi^{0} EPOS3","l");
  legendInvYieldTSallisTheoryOnlyPi0Spectrum->AddEntry(graphAsymmErrCGCTheoryPi0y0pA5020,"#pi^{0} CGC MV^{#gamma}","l");
   legendInvYieldTSallisTheoryOnlyPi0Spectrum->AddEntry(histoDPMJetPi0,"DPMJet","l");
  legendInvYieldTSallisTheoryOnlyPi0Spectrum->AddEntry(histoHIJINGPi0,"HIJING","l");
  legendInvYieldTSallisTheoryOnlyPi0Spectrum->AddEntry(FitCombPi0,"Tsallis fit","l");
  legendInvYieldTSallisTheoryOnlyPi0Spectrum->Draw();
  
  
  
   padComparisonInvYieldTSallisTheoryOnlyPi0Spectrum->Update();

   canvasInvYieldTSallisTheoryOnlyPi0Spectrum->Update();
	
   if(EPOS) canvasInvYieldTSallisTheoryOnlyPi0Spectrum->Print(Form("%s/Pi0MesonYield_EPOSv2.%s",outputDir.Data(),suffix.Data()));
   else if(mT) canvasInvYieldTSallisTheoryOnlyPi0Spectrum->Print(Form("%s/Pi0MesonYields_mTv2.%s",outputDir.Data(),suffix.Data()));
   else canvasInvYieldTSallisTheoryOnlyPi0Spectrum->Print(Form("%s/Pi0MesonYieldsv2.%s",outputDir.Data(),suffix.Data()));
  
   
   
  ///////////////////////////////////////////////////////////////
   
  TCanvas* canvasInvYieldTSallisTheoryOnlyEtaSpectrum = new TCanvas("canvasInvYieldTSallisTheoryOnlyEtaSpectrum","",200,10,1000,1200);  // gives the page size
  DrawGammaCanvasSettings( canvasInvYieldTSallisTheoryOnlyEtaSpectrum,  0.3, 0.02, 0.02, 0.16);
	
  TPad* padComparisonInvYieldTSallisTheoryOnlyEtaSpectrum = new TPad("padComparisonInvYieldTSallisTheoryOnlyEtaSpectrum", "", 0., 0., 1., 1.,-1, -1, -2);
  DrawGammaPadSettings( padComparisonInvYieldTSallisTheoryOnlyEtaSpectrum, 0.15, 0.02, 0.02, 0.09);
  padComparisonInvYieldTSallisTheoryOnlyEtaSpectrum->Draw();
	
  	
  padComparisonInvYieldTSallisTheoryOnlyEtaSpectrum->cd();
  padComparisonInvYieldTSallisTheoryOnlyEtaSpectrum->SetLogy();		
  padComparisonInvYieldTSallisTheoryOnlyEtaSpectrum->SetLogx();		
	
  TH2F * histo2DInvYieldTSallisTheoryOnlyEtaSpectrum;
  histo2DInvYieldTSallisTheoryOnlyEtaSpectrum = new TH2F("histo2DInvYieldTSallisTheoryOnlyEtaSpectrum","histo2DInvYieldTSallisTheoryOnlyEtaSpectrum",1000,0.2,30.,1000,1.2e-8,30 );
  SetStyleHistoTH2ForGraphs(histo2DInvYieldTSallisTheoryOnlyEtaSpectrum, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2} ",0.035,0.035,0.035,0.035, 1.0,1.8, 512, 510);
  histo2DInvYieldTSallisTheoryOnlyEtaSpectrum->DrawCopy(); 
	
 
  
  FitCombEta->SetLineWidth(2);  
  FitCombEta->SetLineColor(kBlack);
  FitCombEta->SetLineStyle(4);
  FitCombEta->Draw("same"); 
 
 
  if (EPOS){
    
    // graphEtaEPOS->SetLineColor(kGreen-3);
    //graphEtaEPOS->SetFillColor(kGreen-3);
    graphEtaEPOS->SetLineColor(kAzure+7);
    graphEtaEPOS->SetFillColor(kAzure+7);
    
    graphEtaEPOS->SetLineWidth(3);
    //graphEtaEPOS->SetLineColor(kGreen+3);
   // graphEtaEPOS->SetFillColor(kGreen+3);
    graphEtaEPOS->Draw("C3same"); //C
    
  }
  
  
  
  if( HIJINGPi0 ){
    SetStyleHisto(histoHIJINGEta, 3, styleLineHIJING, colorHIJING);          
    histoHIJINGEta->SetLineWidth(3);
    histoHIJINGEta->Draw("same,hist,l");
    
  }
  
  if( DPMJetPi0 ){
    
   SetStyleHisto(histoDPMJetEta, 3, styleLineDPMJet, colorDPMJet);          
   histoDPMJetEta->Draw("same,hist,l"); 
   
  }
  
  CombEtaSyst->Draw("E2,same"); 
  CombEtaStat->Draw("pz,same"); 
 
 
 
  
  TLatex * latexInvYieldTSallisTheoryOnlyEtaSpectrum = new TLatex(3,9,"ALICE") ;
  latexInvYieldTSallisTheoryOnlyEtaSpectrum->SetTextColor(kBlack) ;
  latexInvYieldTSallisTheoryOnlyEtaSpectrum->SetTextSize(TextSize-0.015) ;
  latexInvYieldTSallisTheoryOnlyEtaSpectrum->SetTextFont(Font) ;
  latexInvYieldTSallisTheoryOnlyEtaSpectrum->DrawLatex(3,4,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latexInvYieldTSallisTheoryOnlyEtaSpectrum->Draw() ;
  
  

  TLegend* legendInvYieldTSallisTheoryOnlyEtaSpectrum;
  legendInvYieldTSallisTheoryOnlyEtaSpectrum   = new TLegend(0.20,0.12,0.6,0.33); 
  legendInvYieldTSallisTheoryOnlyEtaSpectrum->SetFillColor(0);
  legendInvYieldTSallisTheoryOnlyEtaSpectrum->SetLineColor(0);
  legendInvYieldTSallisTheoryOnlyEtaSpectrum->SetTextFont(Font);
  legendInvYieldTSallisTheoryOnlyEtaSpectrum->SetTextSize(TextSize-0.015);
  legendInvYieldTSallisTheoryOnlyEtaSpectrum->AddEntry(CombEtaSyst,"NSD #eta","pef");
  legendInvYieldTSallisTheoryOnlyEtaSpectrum->AddEntry(graphEtaEPOS,"#eta EPOS3","l");
  legendInvYieldTSallisTheoryOnlyEtaSpectrum->AddEntry(histoDPMJetEta,"DPMJet","l");
  legendInvYieldTSallisTheoryOnlyEtaSpectrum->AddEntry(histoHIJINGEta,"HIJING","l");
  legendInvYieldTSallisTheoryOnlyEtaSpectrum->AddEntry(FitCombEta,"Tsallis fit","l");
  legendInvYieldTSallisTheoryOnlyEtaSpectrum->Draw();
  
  
  
   padComparisonInvYieldTSallisTheoryOnlyEtaSpectrum->Update();

   canvasInvYieldTSallisTheoryOnlyEtaSpectrum->Update();
	
   if(EPOS) canvasInvYieldTSallisTheoryOnlyEtaSpectrum->Print(Form("%s/EtaMesonYield_EPOSv2.%s",outputDir.Data(),suffix.Data()));
   else if(mT) canvasInvYieldTSallisTheoryOnlyEtaSpectrum->Print(Form("%s/EtaMesonYields_mTv2.%s",outputDir.Data(),suffix.Data()));
   else canvasInvYieldTSallisTheoryOnlyEtaSpectrum->Print(Form("%s/EtaMesonYieldsv2.%s",outputDir.Data(),suffix.Data()));
  
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
	
   //*************** first Column **********************************************************
   TLatex *textFitTsallisSpectrumOnlyRatioPi0 = new TLatex(columnsTheoryLegendPrelAndFinalOnlyRatio[0],rowsTheoryLegendPrelAndFinalOnlyRatioPi0[3],"Tsallis fit");
   SetStyleTLatex( textFitTsallisSpectrumOnlyRatioPi0, textSizeLeftLabelsPrelAndFinalRatio,4);
	
	
   TLatex *textEPOS3OnlyRatio = new TLatex(columnsTheoryLegendPrelAndFinalOnlyRatio[0],rowsTheoryLegendPrelAndFinalOnlyRatioPi0[4],"EPOS3");
   SetStyleTLatex( textEPOS3OnlyRatio, textSizeLeftLabelsPrelAndFinalRatio,4);
 	
   TLatex *textCGCOnlyRatioPi0 = new TLatex(columnsTheoryLegendPrelAndFinalOnlyRatio[0],rowsTheoryLegendPrelAndFinalOnlyRatioPi0[5],"CGC");
   SetStyleTLatex( textCGCOnlyRatioPi0, textSizeLeftLabelsPrelAndFinalRatio,4);

   //*************** second Column **********************************************************
	
   TLatex *textPi0pPb5023GeVOnlyRatio = new TLatex(columnsTheoryLegendPrelAndFinalOnlyRatio[1],rowsTheoryLegendPrelAndFinalOnlyRatioPi0[0],"#pi^{0}, #sqrt{#it{s_{NN}}} = 5.02 TeV");
   SetStyleTLatex( textPi0pPb5023GeVOnlyRatio, textSizeTopLablesPrelAndFinalRatio,4);
   TLatex *textPi0pPb5023GeVsysOnlyRatio = new TLatex(columnsTheoryLegendPrelAndFinalOnlyRatio[1],rowsTheoryLegendPrelAndFinalOnlyRatioPi0[1],"syst. + stat.");
   SetStyleTLatex( textPi0pPb5023GeVsysOnlyRatio, textSizeTopLowerLablesPrelAndFinalRatio,4);
	
	
   TBox* boxCombinedPi0pPb5023GeVOnlyRatio = CreateBoxFromGraph(RatioTsallisCombPi0Syst,columnsTheoryLegendPrelAndFinalOnlyRatio[1]+offsetTheoryLegendPrelAndFinalMarkerX-offsetTheoryLegendPrelAndFinalLine, rowsTheoryLegendPrelAndFinalOnlyRatioPi0[2]+ offsetTheoryLegendPrelAndFinalMarkerY- offsetTheoryLegendPrelAndFinalBox, columnsTheoryLegendPrelAndFinalOnlyRatio[1]+offsetTheoryLegendPrelAndFinalMarkerX+offsetTheoryLegendPrelAndFinalLine, rowsTheoryLegendPrelAndFinalOnlyRatioPi0[2]+ offsetTheoryLegendPrelAndFinalMarkerY+offsetTheoryLegendPrelAndFinalBox);
   TMarker* markerCombinedPi0pPb5023GeVOnlyRatio = CreateMarkerFromGraph(RatioTsallisCombPi0Syst,columnsTheoryLegendPrelAndFinalOnlyRatio[1]+offsetTheoryLegendPrelAndFinalMarkerX,rowsTheoryLegendPrelAndFinalOnlyRatioPi0[2]+offsetTheoryLegendPrelAndFinalMarkerY ,scaleMarkerTheory);
	
	
		
	TLine * lineTsallisFitPi0pPb5023GeVOnlyRatio = CreateLineFromGraph(RatioTsallisCombPi0Syst, 	columnsTheoryLegendPrelAndFinalOnlyRatio[1]+ offsetTheoryLegendPrelAndFinalMarkerX- offsetTheoryLegendPrelAndFinalLine, rowsTheoryLegendPrelAndFinalOnlyRatioPi0[3]+offsetTheoryLegendPrelAndFinalMarkerY, columnsTheoryLegendPrelAndFinalOnlyRatio[1]+ offsetTheoryLegendPrelAndFinalMarkerX+ offsetTheoryLegendPrelAndFinalLine, rowsTheoryLegendPrelAndFinalOnlyRatioPi0[3]+offsetTheoryLegendPrelAndFinalMarkerY, scaleLineWidthTheory*1.5); 
	TLine * lineEPOS3Pi0pPb5023GeVOnlyRatio      = CreateLineFromGraph(RatioEPOSFitPi0,	        columnsTheoryLegendPrelAndFinalOnlyRatio[1]+ offsetTheoryLegendPrelAndFinalMarkerX- offsetTheoryLegendPrelAndFinalLine, rowsTheoryLegendPrelAndFinalOnlyRatioPi0[4]+offsetTheoryLegendPrelAndFinalMarkerY, columnsTheoryLegendPrelAndFinalOnlyRatio[1]+ offsetTheoryLegendPrelAndFinalMarkerX+ offsetTheoryLegendPrelAndFinalLine, rowsTheoryLegendPrelAndFinalOnlyRatioPi0[4]+offsetTheoryLegendPrelAndFinalMarkerY, scaleLineWidthTheory*1.5); 
	TLine * lineCGCPi0pPb5023GeVOnlyRatio        = CreateLineFromGraph(RatioCGCFitPi0, 	        columnsTheoryLegendPrelAndFinalOnlyRatio[1]+ offsetTheoryLegendPrelAndFinalMarkerX- offsetTheoryLegendPrelAndFinalLine, rowsTheoryLegendPrelAndFinalOnlyRatioPi0[5]+offsetTheoryLegendPrelAndFinalMarkerY, columnsTheoryLegendPrelAndFinalOnlyRatio[1]+ offsetTheoryLegendPrelAndFinalMarkerX+ offsetTheoryLegendPrelAndFinalLine, rowsTheoryLegendPrelAndFinalOnlyRatioPi0[5]+offsetTheoryLegendPrelAndFinalMarkerY, scaleLineWidthTheory*1.5); 
	
	//*************** third Column **********************************************************
	TLatex *textEtapPb5023GeVOnlyRatio = new TLatex(columnsTheoryLegendPrelAndFinalOnlyRatio[2],rowsTheoryLegendPrelAndFinalOnlyRatioPi0[0],"#eta, #sqrt{#it{s_{NN}}} = 5.02 TeV");
	SetStyleTLatex( textEtapPb5023GeVOnlyRatio, textSizeTopLablesPrelAndFinalRatio,4);
	TLatex *textEtapPb5023GeVsysOnlyRatio = new TLatex(columnsTheoryLegendPrelAndFinalOnlyRatio[2],rowsTheoryLegendPrelAndFinalOnlyRatioPi0[1],"syst. + stat.");
	SetStyleTLatex( textEtapPb5023GeVsysOnlyRatio, textSizeTopLowerLablesPrelAndFinalRatio,4);
	
	TBox* boxCombinedEtapPb5023GeVOnlyRatio = CreateBoxFromGraph(RatioTsallisCombEtaSyst,columnsTheoryLegendPrelAndFinalOnlyRatio[2]+offsetTheoryLegendPrelAndFinalMarkerX-offsetTheoryLegendPrelAndFinalLine, rowsTheoryLegendPrelAndFinalOnlyRatioPi0[2]+ offsetTheoryLegendPrelAndFinalMarkerY- offsetTheoryLegendPrelAndFinalBox, columnsTheoryLegendPrelAndFinalOnlyRatio[2]+offsetTheoryLegendPrelAndFinalMarkerX+offsetTheoryLegendPrelAndFinalLine, rowsTheoryLegendPrelAndFinalOnlyRatioPi0[2]+ offsetTheoryLegendPrelAndFinalMarkerY+offsetTheoryLegendPrelAndFinalBox);
	TMarker* markerCombinedEtapPb5023GeVOnlyRatio = CreateMarkerFromGraph(RatioTsallisCombEtaSyst,columnsTheoryLegendPrelAndFinalOnlyRatio[2]+offsetTheoryLegendPrelAndFinalMarkerX,rowsTheoryLegendPrelAndFinalOnlyRatioPi0[2]+offsetTheoryLegendPrelAndFinalMarkerY ,scaleMarkerTheory);
	
	
	
	
	
	TLine * lineTsallisFitEtapPb5023GeVOnlyRatio = CreateLineFromGraph(RatioTsallisCombEtaSyst,  	columnsTheoryLegendPrelAndFinalOnlyRatio[2]+ offsetTheoryLegendPrelAndFinalMarkerX- offsetTheoryLegendPrelAndFinalLine, rowsTheoryLegendPrelAndFinalOnlyRatioPi0[3]+offsetTheoryLegendPrelAndFinalMarkerY, columnsTheoryLegendPrelAndFinalOnlyRatio[2]+ offsetTheoryLegendPrelAndFinalMarkerX+ offsetTheoryLegendPrelAndFinalLine,rowsTheoryLegendPrelAndFinalOnlyRatioPi0[3]+offsetTheoryLegendPrelAndFinalMarkerY, scaleLineWidthTheory*1.5); 
	TLine * lineEPOS3EtapPb5023GeVOnlyRatio  = CreateLineFromGraph(RatioEPOSFitEta,  	columnsTheoryLegendPrelAndFinalOnlyRatio[2]+ offsetTheoryLegendPrelAndFinalMarkerX- offsetTheoryLegendPrelAndFinalLine, rowsTheoryLegendPrelAndFinalOnlyRatioPi0[4]+offsetTheoryLegendPrelAndFinalMarkerY, columnsTheoryLegendPrelAndFinalOnlyRatio[2]+ offsetTheoryLegendPrelAndFinalMarkerX+ offsetTheoryLegendPrelAndFinalLine,rowsTheoryLegendPrelAndFinalOnlyRatioPi0[4]+offsetTheoryLegendPrelAndFinalMarkerY, scaleLineWidthTheory*1.5); 
	
	
	
	
	
	TCanvas* canvasInvYieldTSallisTheoryOnlyRatioPi0AndEtaSpectra 	= new TCanvas("canvasInvYieldTSallisTheoryOnlyRatioPi0AndEtaSpectra","",200,10,1000,800);  // gives the page size
	DrawGammaCanvasSettings( canvasInvYieldTSallisTheoryOnlyRatioPi0AndEtaSpectra,  0.15, 0.02, 0.00, 0.00);  

	TPad* padComparisonInvYieldTSallisTheoryOnlyRatioPi0 		= new TPad("padComparisonInvYieldTSallisTheoryOnlyRatioPi0", "", 0.,0.55, 1., 1,   -1, -1, -2);
	DrawGammaPadSettings( padComparisonInvYieldTSallisTheoryOnlyRatioPi0,  0.15, 0.02, 0.02, 0.0);
	padComparisonInvYieldTSallisTheoryOnlyRatioPi0->Draw();
	
 	TPad* padComparisonInvYieldTSallisTheoryOnlyRatioEta		= new TPad("padComparisonInvYieldTSallisTheoryOnlyRatioEta", "", 0.,0.0,  1., 0.55,-1, -1, -2);
 	DrawGammaPadSettings( padComparisonInvYieldTSallisTheoryOnlyRatioEta,  0.15, 0.02, 0.00,  0.20);
 	padComparisonInvYieldTSallisTheoryOnlyRatioEta->Draw();
 	
 	
 	padComparisonInvYieldTSallisTheoryOnlyRatioPi0->cd();
 	padComparisonInvYieldTSallisTheoryOnlyRatioPi0->SetLogx();
	
	
 	
	TH2F * histo2DInvYieldTSallisTheoryOnlyRatioPi0;
	histo2DInvYieldTSallisTheoryOnlyRatioPi0 = new TH2F("histo2DInvYieldTSallisTheoryOnlyRatioPi0","histo2DRatioAllppreferencesEtaandPi0",1000,.2,40.,1000,0.001,4.2);//1.89
	SetStyleHistoTH2ForGraphs(histo2DInvYieldTSallisTheoryOnlyRatioPi0, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.098,0.098, 0.098,0.098, 1.,0.4, 502, 505); 
	histo2DInvYieldTSallisTheoryOnlyRatioPi0->GetYaxis()->SetLabelOffset(0.005);
	histo2DInvYieldTSallisTheoryOnlyRatioPi0->GetXaxis()->SetLabelOffset(LabelOffsetLog);
	histo2DInvYieldTSallisTheoryOnlyRatioPi0->GetXaxis()->SetTickLength(0.07);
	histo2DInvYieldTSallisTheoryOnlyRatioPi0->GetYaxis()-> CenterTitle();
	histo2DInvYieldTSallisTheoryOnlyRatioPi0->DrawCopy();
	
	TLegend* legendInvYieldTSallisTheoryOnlyRatioPi0;
	legendInvYieldTSallisTheoryOnlyRatioPi0   = new TLegend(0.18,0.81,0.47,0.96); 
	legendInvYieldTSallisTheoryOnlyRatioPi0-> SetNColumns(2);
	legendInvYieldTSallisTheoryOnlyRatioPi0->SetFillColor(0);
	legendInvYieldTSallisTheoryOnlyRatioPi0->SetFillStyle(0);
	legendInvYieldTSallisTheoryOnlyRatioPi0->SetLineColor(0);
	legendInvYieldTSallisTheoryOnlyRatioPi0->SetTextFont(Font);
	legendInvYieldTSallisTheoryOnlyRatioPi0->SetTextSize(TextSize+0.045);
	legendInvYieldTSallisTheoryOnlyRatioPi0->AddEntry(RatioTsallisCombPi0Syst,"#pi^{0}","pef");
        //legendInvYieldTSallisTheoryOnlyRatioPi0->Draw();

	
	DrawGammaLines(0., 40.,1., 1.,2.0,kRed+2,2);
	
	RatioIlkkaFitPi0scaleerr->SetLineWidth(3);
	RatioIlkkaFitPi0scaleerr->SetLineColor(colorIlkkaPi0);
	RatioIlkkaFitPi0scaleerr->SetFillColor(colorIlkkaPi0);
	RatioIlkkaFitPi0scaleerr->Draw("C3same");
	
	
	RatioIlkkaFitPi0NoErr->SetLineWidth(3);
	RatioIlkkaFitPi0NoErr->SetLineColor(41);
	RatioIlkkaFitPi0NoErr->SetFillColor(41);
	//RatioIlkkaFitPi0NoErr->Draw("lE2same");
	
	if (EPOS){
	  while(RatioEPOSFitPi0->GetX()[0] < 0.7)
          RatioEPOSFitPi0->RemovePoint(0);   
	  RatioEPOSFitPi0->Draw("C3same");
	}
  
	if( CGCPi0){
	  while(RatioCGCFitPi0->GetX()[0] < 0.3)
	    RatioCGCFitPi0->RemovePoint(0); 
	    RatioCGCFitPi0->Draw("C3same");
	}
	
	
	RatioMcGillFitPi0->SetLineWidth(2);
	RatioMcGillFitPi0->SetLineColor(colorMcGill);
	RatioMcGillFitPi0->SetFillColor(colorMcGill);
	 while(RatioMcGillFitPi0->GetX()[0] < 0.3)
          RatioMcGillFitPi0->RemovePoint(0);   
	  RatioMcGillFitPi0->Draw("C3same");
	
	
	
	if( HIJINGPi0 ){
           SetStyleHisto(histoRatioPi0HIJINGToFit,3, styleLineHIJING, colorHIJING );  
	   histoRatioPi0HIJINGToFit->GetXaxis()->SetRangeUser(0.3,20.0);
           histoRatioPi0HIJINGToFit->Draw("same,hist,l");  
	}

	if( DPMJetPi0 ){
	   SetStyleHisto(histoRatioPi0DPMJetToFit,3, styleLineDPMJet, colorDPMJet );  
	   histoRatioPi0HIJINGToFit->GetXaxis()->SetRangeUser(0.3,20.0);
	   histoRatioPi0DPMJetToFit->Draw("same,hist,l");  
	}
  
   
	RatioTsallisCombPi0Syst->Draw("E2,same");
	RatioTsallisCombPi0Stat->Draw("Ez,p,same");  
	
	TLatex * latexInvYieldTSallisTheoryOnlyRatioPi0 = new TLatex(2.8,3.7,"ALICE") ;
	latexInvYieldTSallisTheoryOnlyRatioPi0->SetTextColor(kBlack) ;
	latexInvYieldTSallisTheoryOnlyRatioPi0->SetTextSize(TextSize+0.045) ;
	latexInvYieldTSallisTheoryOnlyRatioPi0->SetTextFont(Font) ;
	
	latexInvYieldTSallisTheoryOnlyRatioPi0->DrawLatex(2.8,3.3,"p-Pb, NSD, #sqrt{#it{s}_{NN}} = 5.02 TeV");
	latexInvYieldTSallisTheoryOnlyRatioPi0->Draw() ;
	
	
  
	TLegend* legendInvYieldTSallisTheoryOnlyRatioPi01;
	legendInvYieldTSallisTheoryOnlyRatioPi01   = new TLegend(0.18,0.63,0.38,0.93); 
	//legendInvYieldTSallisTheoryOnlyRatioPi01-> SetNColumns(2);
	legendInvYieldTSallisTheoryOnlyRatioPi01->SetFillColor(0);
	legendInvYieldTSallisTheoryOnlyRatioPi01->SetFillStyle(0);
	legendInvYieldTSallisTheoryOnlyRatioPi01->SetLineColor(0);
	legendInvYieldTSallisTheoryOnlyRatioPi01->SetTextFont(Font);
	legendInvYieldTSallisTheoryOnlyRatioPi01->SetTextSize(TextSize+0.035);
	//legendInvYieldTSallisTheoryOnlyRatioPi01->AddEntry(RatioEPOSFitPi0,"EPOS3","l");
	legendInvYieldTSallisTheoryOnlyRatioPi01->AddEntry(RatioTsallisCombPi0Syst,"#pi^{0}","pef");
       
 	legendInvYieldTSallisTheoryOnlyRatioPi01->AddEntry(RatioCGCFitPi0,"CGC MV^{#gamma}","l");
	//legendInvYieldTSallisTheoryOnlyRatioPi01->AddEntry(RatioIlkkaFitPi0NoErr,"NLO pQCD","l");
	legendInvYieldTSallisTheoryOnlyRatioPi01->AddEntry(RatioIlkkaFitPi0scaleerr,"NLO pQCD","l");
	//legendInvYieldTSallisTheoryOnlyRatioPi01->AddEntry(RatioMcGillFitPi0,"C. Shen et al.","l");
	
	//legendInvYieldTSallisTheoryOnlyRatioPi01->AddEntry(histoRatioPi0HIJINGToFit,"HIJING","l");
	//legendInvYieldTSallisTheoryOnlyRatioPi01->AddEntry(histoRatioPi0DPMJetToFit,"DPMJet","l");
	//legendInvYieldTSallisTheoryOnlyRatioPi01->AddEntry(RatioMcGillFitPi0,"C. Shen et al.","l");
	//legendInvYieldTSallisTheoryOnlyRatioPi01->AddEntry(RatioIlkkaFitPi0scaleerr,"pQCD(epss16,dss14,ct14)","l");

	legendInvYieldTSallisTheoryOnlyRatioPi01->Draw();
	
	
	TLegend* legendInvYieldTSallisTheoryOnlyRatioPi02;
	legendInvYieldTSallisTheoryOnlyRatioPi02   = new TLegend(0.59,0.48,0.73,0.68); 
	legendInvYieldTSallisTheoryOnlyRatioPi02->SetFillColor(0);
	legendInvYieldTSallisTheoryOnlyRatioPi02->SetFillStyle(0);
	legendInvYieldTSallisTheoryOnlyRatioPi02->SetLineColor(0);
	legendInvYieldTSallisTheoryOnlyRatioPi02->SetTextFont(Font);
	legendInvYieldTSallisTheoryOnlyRatioPi02->SetTextSize(TextSize+0.035);
	legendInvYieldTSallisTheoryOnlyRatioPi02->AddEntry(RatioIlkkaFitPi0NoErr,"NLO pQCD","l");
	//legendInvYieldTSallisTheoryOnlyRatioPi02->AddEntry(RatioIlkkaFitPi0scaleerr,"epss16,dss14,ct14,scale","ef");
	
	//legendInvYieldTSallisTheoryOnlyRatioPi02->Draw();

	padComparisonInvYieldTSallisTheoryOnlyRatioEta->cd();
	padComparisonInvYieldTSallisTheoryOnlyRatioEta->SetLogx();
      
      
	TH2F * histo2DInvYieldTSallisTheoryOnlyRatioEta;
	histo2DInvYieldTSallisTheoryOnlyRatioEta = new TH2F("histo2DInvYieldTSallisTheoryOnlyRatioEta","histo2DRatioAllppreferencesEtaandEta",1000,.2,40.,1000,0.001,4.2);//1.89
	SetStyleHistoTH2ForGraphs(histo2DInvYieldTSallisTheoryOnlyRatioEta, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}",0.08,0.08,0.08,0.08,1.,0.5, 502, 505); 
	histo2DInvYieldTSallisTheoryOnlyRatioEta->GetYaxis()->SetLabelOffset(0.005);
	histo2DInvYieldTSallisTheoryOnlyRatioEta->GetXaxis()->SetLabelOffset(LabelOffsetLog);
	histo2DInvYieldTSallisTheoryOnlyRatioEta->GetXaxis()->SetTickLength(0.07);
	histo2DInvYieldTSallisTheoryOnlyRatioEta->GetYaxis()-> CenterTitle();
	histo2DInvYieldTSallisTheoryOnlyRatioEta->DrawCopy();

	
	DrawGammaLines(0., 40.,1., 1.,2.0,kRed+2,2);
	
	if (EPOS){
	 
	  while(RatioEPOSFitEta->GetX()[0] < 0.7)
          RatioEPOSFitEta->RemovePoint(0);   
	  
	 RatioEPOSFitEta->Draw("C3same"); 
	
      	 
	 
	}
  
	if( HIJINGPi0 ){
           SetStyleHisto(histoRatioEtaHIJINGToFit,3, styleLineHIJING, colorHIJING );  
           histoRatioEtaHIJINGToFit->Draw("same,hist,l");  
	   histoRatioEtaHIJINGToFit->GetXaxis()->SetRangeUser(0.7,20);
	}

	if( DPMJetPi0 ){
	   SetStyleHisto(histoRatioEtaDPMJetToFit,3, styleLineDPMJet, colorDPMJet );  
	   
	   histoRatioEtaDPMJetToFit->Draw("same,hist,l");  
	   histoRatioEtaDPMJetToFit->GetXaxis()->SetRangeUser(0.7,20);
	}
	
	
	RatioMcGillFitEta->SetLineWidth(2);
	RatioMcGillFitEta->SetLineColor(colorMcGill);
	RatioMcGillFitEta->SetFillColor(colorMcGill);
	 while(RatioMcGillFitEta->GetX()[0] < 0.7)
          RatioMcGillFitEta->RemovePoint(0);   
	  
	RatioMcGillFitEta->Draw("C3same");
	
  
  
	RatioTsallisCombEtaSyst->Draw("E2,same");
        RatioTsallisCombEtaStat->Draw("Ez,p,same");  
	
	
	TLegend* legendInvYieldTSallisTheoryOnlyRatioEta;
	legendInvYieldTSallisTheoryOnlyRatioEta   = new TLegend(0.18,0.55,0.38,0.95); 
	legendInvYieldTSallisTheoryOnlyRatioEta->SetFillColor(0);
	legendInvYieldTSallisTheoryOnlyRatioEta->SetLineColor(0);
	legendInvYieldTSallisTheoryOnlyRatioEta->SetTextFont(Font);
	legendInvYieldTSallisTheoryOnlyRatioEta->SetTextSize((TextSize+0.035)*0.82);
	//legendInvYieldTSallisTheoryOnlyRatioEta-> SetNColumns(2);
	legendInvYieldTSallisTheoryOnlyRatioEta->AddEntry(RatioTsallisCombEtaSyst,"#eta","pef");
 	legendInvYieldTSallisTheoryOnlyRatioEta->AddEntry(RatioEPOSFitPi0,"EPOS3","l");
 	legendInvYieldTSallisTheoryOnlyRatioEta->AddEntry(RatioMcGillFitPi0,"VISHNU","l");
	legendInvYieldTSallisTheoryOnlyRatioEta->AddEntry(histoRatioPi0HIJINGToFit,"HIJING","l");
	legendInvYieldTSallisTheoryOnlyRatioEta->AddEntry(histoRatioPi0DPMJetToFit,"DPMJet","l");
	
// 	legendInvYieldTSallisTheoryOnlyRatioEta->AddEntry(RatioCGCFitPi0,"CGC MV^{#gamma}","l");
// 	legendInvYieldTSallisTheoryOnlyRatioEta->AddEntry(RatioIlkkaFitPi0scaleerr,"pQCD(epss16,dss14,ct14)","l");

	legendInvYieldTSallisTheoryOnlyRatioEta->Draw();
	
	
	//legendInvYieldTSallisTheoryOnlyRatioPi02->Draw();
  

	canvasInvYieldTSallisTheoryOnlyRatioPi0AndEtaSpectra->Update();
	
	if(EPOS) canvasInvYieldTSallisTheoryOnlyRatioPi0AndEtaSpectra->Print(Form("%s/MesonYields_EPOS_OnlyRatiov2.%s",outputDir.Data(),suffix.Data()));
        else if(mT) canvasInvYieldTSallisTheoryOnlyRatioPi0AndEtaSpectra->Print(Form("%s/MesonYields_mT_OnlyRatiov2.%s",outputDir.Data(),suffix.Data()));
        else canvasInvYieldTSallisTheoryOnlyRatioPi0AndEtaSpectra->Print(Form("%s/MesonYields_OnlyRatiov2.%s",outputDir.Data(),suffix.Data()));
  
	
	
	
	TCanvas* canvasInvYieldTSallisTheoryOnlyRatioPi0Spectrum 	= new TCanvas("canvasInvYieldTSallisTheoryOnlyRatioPi0Spectrum","",200,10,1000,800);  // gives the page size
	DrawGammaCanvasSettings( canvasInvYieldTSallisTheoryOnlyRatioPi0Spectrum,  0.15, 0.02, 0.00, 0.00);  

	TPad* padComparisonInvYieldTSallisTheoryOnlyRatioPi0v2 		= new TPad("padComparisonInvYieldTSallisTheoryOnlyRatioPi0v2", "", 0.,0.55, 1., 1,   -1, -1, -2);
	DrawGammaPadSettings( padComparisonInvYieldTSallisTheoryOnlyRatioPi0v2,  0.15, 0.02, 0.02, 0.0);
	padComparisonInvYieldTSallisTheoryOnlyRatioPi0v2->Draw();
	
 	TPad* padComparisonInvYieldTSallisGenOnlyRatioPi0		= new TPad("padComparisonInvYieldTSallisGenOnlyRatioPi0", "", 0.,0.0,  1., 0.55,-1, -1, -2);
 	DrawGammaPadSettings( padComparisonInvYieldTSallisGenOnlyRatioPi0,  0.15, 0.02, 0.00,  0.20);
 	padComparisonInvYieldTSallisGenOnlyRatioPi0->Draw();
 	
 	
 	padComparisonInvYieldTSallisTheoryOnlyRatioPi0v2->cd();
 	padComparisonInvYieldTSallisTheoryOnlyRatioPi0v2->SetLogx();
	
	
 	
	TH2F * histo2DInvYieldTSallisTheoryOnlyRatioPi0v2;
	histo2DInvYieldTSallisTheoryOnlyRatioPi0v2 = new TH2F("histo2DInvYieldTSallisTheoryOnlyRatioPi0v2","histo2DRatioAllppreferencesEtaandPi0",1000,.2,30.,1000,0.21,2.69);//1.89
	SetStyleHistoTH2ForGraphs(histo2DInvYieldTSallisTheoryOnlyRatioPi0v2, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.098,0.098, 0.098,0.098, 1.,0.6, 502, 505); 
	histo2DInvYieldTSallisTheoryOnlyRatioPi0v2->GetYaxis()->SetLabelOffset(0.005);
	histo2DInvYieldTSallisTheoryOnlyRatioPi0v2->GetXaxis()->SetLabelOffset(LabelOffsetLog);
	histo2DInvYieldTSallisTheoryOnlyRatioPi0v2->GetXaxis()->SetTickLength(0.07);
	histo2DInvYieldTSallisTheoryOnlyRatioPi0v2->GetYaxis()-> CenterTitle();
	histo2DInvYieldTSallisTheoryOnlyRatioPi0v2->DrawCopy();

	
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	if (EPOS){
	  
      	 RatioEPOSFitPi0->Draw("C3same");
	 
	}
  
	if( CGCPi0){
	  
           RatioCGCFitPi0->Draw("C3same");
	}
	
	   
	RatioTsallisCombPi0Syst->Draw("E2,same");
	RatioTsallisCombPi0Stat->Draw("Ez,p,same");  
	
	TLatex * latexInvYieldTSallisTheoryOnlyRatioPi0v2 = new TLatex(3.5,2.4,"ALICE") ;
	latexInvYieldTSallisTheoryOnlyRatioPi0v2->SetTextColor(kBlack) ;
	latexInvYieldTSallisTheoryOnlyRatioPi0v2->SetTextSize(0.098) ;
	latexInvYieldTSallisTheoryOnlyRatioPi0v2->SetTextFont(Font) ;
	
	latexInvYieldTSallisTheoryOnlyRatioPi0v2->DrawLatex(3.5,2.1,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
	latexInvYieldTSallisTheoryOnlyRatioPi0v2->Draw() ;
	
	TLegend* legendInvYieldTSallisTheoryOnlyRatioPi0v2;
	legendInvYieldTSallisTheoryOnlyRatioPi0v2   = new TLegend(0.18,0.52,0.38,0.95); 
	legendInvYieldTSallisTheoryOnlyRatioPi0v2->SetFillColor(0);
	legendInvYieldTSallisTheoryOnlyRatioPi0v2->SetLineColor(0);
	legendInvYieldTSallisTheoryOnlyRatioPi0v2->SetTextFont(Font);
	legendInvYieldTSallisTheoryOnlyRatioPi0v2->SetTextSize(0.098);
	legendInvYieldTSallisTheoryOnlyRatioPi0v2->AddEntry(RatioTsallisCombPi0Syst,"NSD #pi^{0}","pef");
	legendInvYieldTSallisTheoryOnlyRatioPi0v2->AddEntry(RatioEPOSFitPi0,"EPOS3","l");
	legendInvYieldTSallisTheoryOnlyRatioPi0v2->AddEntry(RatioCGCFitPi0,"CGC MV^{#gamma}","l");
	
	legendInvYieldTSallisTheoryOnlyRatioPi0v2->Draw();
  

	padComparisonInvYieldTSallisGenOnlyRatioPi0->cd();
	padComparisonInvYieldTSallisGenOnlyRatioPi0->SetLogx();
      
      
	TH2F * histo2DInvYieldTSallisGenOnlyRatioPi0;
	histo2DInvYieldTSallisGenOnlyRatioPi0 = new TH2F("histo2DInvYieldTSallisGenOnlyRatioPi0","histo2DRatioAllppreferencesEtaandEta",1000,.2,30.,1000,0.21,2.69);//1.89
	SetStyleHistoTH2ForGraphs(histo2DInvYieldTSallisGenOnlyRatioPi0, "#it{p}_{T} (GeV/#it{c})","#frac{Generator, Data}{fit}",0.08,0.08,0.08,0.08,1.,0.7, 502, 505); 
	histo2DInvYieldTSallisGenOnlyRatioPi0->GetYaxis()->SetLabelOffset(0.005);
	histo2DInvYieldTSallisGenOnlyRatioPi0->GetXaxis()->SetLabelOffset(LabelOffsetLog);
	histo2DInvYieldTSallisGenOnlyRatioPi0->GetXaxis()->SetTickLength(0.07);
	histo2DInvYieldTSallisGenOnlyRatioPi0->GetYaxis()-> CenterTitle();
	histo2DInvYieldTSallisGenOnlyRatioPi0->DrawCopy();

	
	DrawGammaLines(0., 30.,1., 1.,0.1);

	
	if( HIJINGPi0 ){
           SetStyleHisto(histoRatioPi0HIJINGToFit,3, styleLineHIJING, colorHIJING );  
           histoRatioPi0HIJINGToFit->Draw("same,hist,l");  
	}

	if( DPMJetPi0 ){
	   SetStyleHisto(histoRatioPi0DPMJetToFit,3, styleLineDPMJet, colorDPMJet );  
	   histoRatioPi0DPMJetToFit->Draw("same,hist,l");  
	}
  
	
	
	
	RatioTsallisCombPi0Syst->Draw("E2,same");
	RatioTsallisCombPi0Stat->Draw("Ez,p,same");  
	
	TLegend* legendInvYieldTSallisGenOnlyRatioPi0;
	legendInvYieldTSallisGenOnlyRatioPi0   = new TLegend(0.18,0.78,0.38,0.95); 
	legendInvYieldTSallisGenOnlyRatioPi0->SetFillColor(0);
	legendInvYieldTSallisGenOnlyRatioPi0->SetLineColor(0);
	legendInvYieldTSallisGenOnlyRatioPi0->SetTextFont(Font);
	legendInvYieldTSallisGenOnlyRatioPi0->SetTextSize(0.08);
	//legendInvYieldTSallisGenOnlyRatioPi0->AddEntry(RatioTsallisCombEtaSyst,"NSD #pi^{0}","pef");
	legendInvYieldTSallisGenOnlyRatioPi0->AddEntry(histoRatioPi0HIJINGToFit,"HIJING","l");
	legendInvYieldTSallisGenOnlyRatioPi0->AddEntry(histoRatioPi0DPMJetToFit,"DPMJet","l");
	legendInvYieldTSallisGenOnlyRatioPi0->Draw();
  

	canvasInvYieldTSallisTheoryOnlyRatioPi0Spectrum->Update();
	
	if(EPOS) canvasInvYieldTSallisTheoryOnlyRatioPi0Spectrum->Print(Form("%s/Pi0MesonYield_EPOS_OnlyRatiov2.%s",outputDir.Data(),suffix.Data()));
        else if(mT) canvasInvYieldTSallisTheoryOnlyRatioPi0Spectrum->Print(Form("%s/Pi0MesonYield_mT_OnlyRatiov2.%s",outputDir.Data(),suffix.Data()));
        else canvasInvYieldTSallisTheoryOnlyRatioPi0Spectrum->Print(Form("%s/Pi0MesonYield_OnlyRatiov2.%s",outputDir.Data(),suffix.Data()));
  
	
	
	TCanvas* canvasInvYieldTSallisTheoryOnlyRatioEtaSpectrum 	= new TCanvas("canvasInvYieldTSallisTheoryOnlyRatioEtaSpectrum","",200,10,1000,800);  // gives the page size
	DrawGammaCanvasSettings( canvasInvYieldTSallisTheoryOnlyRatioEtaSpectrum,  0.15, 0.02, 0.00, 0.00);  

	TPad* padComparisonInvYieldTSallisTheoryOnlyRatioEtav2 		= new TPad("padComparisonInvYieldTSallisTheoryOnlyRatioEtav2", "", 0.,0.55, 1., 1,   -1, -1, -2);
	DrawGammaPadSettings( padComparisonInvYieldTSallisTheoryOnlyRatioEtav2,  0.15, 0.02, 0.02, 0.0);
	padComparisonInvYieldTSallisTheoryOnlyRatioEtav2->Draw();
	
 	TPad* padComparisonInvYieldTSallisGenOnlyRatioEta		= new TPad("padComparisonInvYieldTSallisGenOnlyRatioEta", "", 0.,0.0,  1., 0.55,-1, -1, -2);
 	DrawGammaPadSettings( padComparisonInvYieldTSallisGenOnlyRatioEta,  0.15, 0.02, 0.00,  0.20);
 	padComparisonInvYieldTSallisGenOnlyRatioEta->Draw();
 	
 	
 	padComparisonInvYieldTSallisTheoryOnlyRatioEtav2->cd();
 	padComparisonInvYieldTSallisTheoryOnlyRatioEtav2->SetLogx();
	
	
 	
	TH2F * histo2DInvYieldTSallisTheoryOnlyRatioEtav2;
	histo2DInvYieldTSallisTheoryOnlyRatioEtav2 = new TH2F("histo2DInvYieldTSallisTheoryOnlyRatioEtav2","histo2DRatioAllppreferencesEtaandPi0",1000,.5,30.,1000,0.21,2.69);//1.89
	SetStyleHistoTH2ForGraphs(histo2DInvYieldTSallisTheoryOnlyRatioEtav2, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.098,0.098, 0.098,0.098, 1.,0.6, 502, 505); 
	histo2DInvYieldTSallisTheoryOnlyRatioEtav2->GetYaxis()->SetLabelOffset(0.005);
	histo2DInvYieldTSallisTheoryOnlyRatioEtav2->GetXaxis()->SetLabelOffset(LabelOffsetLog);
	histo2DInvYieldTSallisTheoryOnlyRatioEtav2->GetXaxis()->SetTickLength(0.07);
	histo2DInvYieldTSallisTheoryOnlyRatioEtav2->GetYaxis()-> CenterTitle();
	histo2DInvYieldTSallisTheoryOnlyRatioEtav2->DrawCopy();

	
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	if (EPOS){
	  
      	 RatioEPOSFitEta->Draw("C3same");
	 
	}
  
	   
	RatioTsallisCombEtaSyst->Draw("E2,same");
	RatioTsallisCombEtaStat->Draw("Ez,p,same");  
	
	TLatex * latexInvYieldTSallisTheoryOnlyRatioEtav2 = new TLatex(2.5,2.4,"ALICE") ;
	latexInvYieldTSallisTheoryOnlyRatioEtav2->SetTextColor(kBlack) ;
	latexInvYieldTSallisTheoryOnlyRatioEtav2->SetTextSize(0.098) ;
	latexInvYieldTSallisTheoryOnlyRatioEtav2->SetTextFont(Font) ;
	
	latexInvYieldTSallisTheoryOnlyRatioEtav2->DrawLatex(2.5,2.1,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
	latexInvYieldTSallisTheoryOnlyRatioEtav2->Draw() ;
	
	TLegend* legendInvYieldTSallisTheoryOnlyRatioEtav2;
	legendInvYieldTSallisTheoryOnlyRatioEtav2   = new TLegend(0.18,0.75,0.38,0.95);  
	legendInvYieldTSallisTheoryOnlyRatioEtav2->SetFillColor(0);
	legendInvYieldTSallisTheoryOnlyRatioEtav2->SetLineColor(0);
	legendInvYieldTSallisTheoryOnlyRatioEtav2->SetTextFont(Font);
	legendInvYieldTSallisTheoryOnlyRatioEtav2->SetTextSize(0.098);
	legendInvYieldTSallisTheoryOnlyRatioEtav2->AddEntry(RatioTsallisCombEtaSyst,"NSD #eta","pef");
	legendInvYieldTSallisTheoryOnlyRatioEtav2->AddEntry(RatioEPOSFitEta,"EPOS3","l");
	
	legendInvYieldTSallisTheoryOnlyRatioEtav2->Draw();
  

	padComparisonInvYieldTSallisGenOnlyRatioEta->cd();
	padComparisonInvYieldTSallisGenOnlyRatioEta->SetLogx();
      
      
	TH2F * histo2DInvYieldTSallisGenOnlyRatioEta;
	histo2DInvYieldTSallisGenOnlyRatioEta = new TH2F("histo2DInvYieldTSallisGenOnlyRatioEta","histo2DRatioAllppreferencesEtaandEta",1000,.5,30.,1000,0.21,2.69);//1.89
	SetStyleHistoTH2ForGraphs(histo2DInvYieldTSallisGenOnlyRatioEta, "#it{p}_{T} (GeV/#it{c})","#frac{Generator, Data}{fit}",0.08,0.08,0.08,0.08,1.,0.7, 502, 505); 
	histo2DInvYieldTSallisGenOnlyRatioEta->GetYaxis()->SetLabelOffset(0.005);
	histo2DInvYieldTSallisGenOnlyRatioEta->GetXaxis()->SetLabelOffset(LabelOffsetLog);
	histo2DInvYieldTSallisGenOnlyRatioEta->GetXaxis()->SetTickLength(0.07);
	histo2DInvYieldTSallisGenOnlyRatioEta->GetYaxis()-> CenterTitle();
	histo2DInvYieldTSallisGenOnlyRatioEta->DrawCopy();

	
	DrawGammaLines(0., 30.,1., 1.,0.1);

	
	if( HIJINGPi0 ){
           SetStyleHisto(histoRatioEtaHIJINGToFit,3, styleLineHIJING, colorHIJING );  
           histoRatioEtaHIJINGToFit->Draw("same,hist,l");  
	}

	if( DPMJetPi0 ){
	   SetStyleHisto(histoRatioEtaDPMJetToFit,3, styleLineDPMJet, colorDPMJet );  
	   histoRatioEtaDPMJetToFit->Draw("same,hist,l");  
	}
  
	
	
	
	RatioTsallisCombEtaSyst->Draw("E2,same");
	RatioTsallisCombEtaStat->Draw("Ez,p,same");  
	
	TLegend* legendInvYieldTSallisGenOnlyRatioEta;
	legendInvYieldTSallisGenOnlyRatioEta   = new TLegend(0.18,0.78,0.38,0.95); 
	legendInvYieldTSallisGenOnlyRatioEta->SetFillColor(0);
	legendInvYieldTSallisGenOnlyRatioEta->SetLineColor(0);
	legendInvYieldTSallisGenOnlyRatioEta->SetTextFont(Font);
	legendInvYieldTSallisGenOnlyRatioEta->SetTextSize(0.08);
	legendInvYieldTSallisGenOnlyRatioEta->AddEntry(histoRatioEtaHIJINGToFit,"HIJING","l");
	legendInvYieldTSallisGenOnlyRatioEta->AddEntry(histoRatioEtaDPMJetToFit,"DPMJet","l");
	legendInvYieldTSallisGenOnlyRatioEta->Draw();
  

	canvasInvYieldTSallisTheoryOnlyRatioEtaSpectrum->Update();
	
	if(EPOS) canvasInvYieldTSallisTheoryOnlyRatioEtaSpectrum->Print(Form("%s/EtaMesonYield_EPOS_OnlyRatiov2.%s",outputDir.Data(),suffix.Data()));
        else if(mT) canvasInvYieldTSallisTheoryOnlyRatioEtaSpectrum->Print(Form("%s/EtaMesonYield_mT_OnlyRatiov2.%s",outputDir.Data(),suffix.Data()));
        else canvasInvYieldTSallisTheoryOnlyRatioEtaSpectrum->Print(Form("%s/EtaMesonYield_OnlyRatiov2.%s",outputDir.Data(),suffix.Data()));
  
	
	
	cout << "works!"<<endl;
}
