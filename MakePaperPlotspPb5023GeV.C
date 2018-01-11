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


TGraphErrors* Eta2Pi0dAu();

void SetEx(TGraphAsymmErrors* gae, Double_t Ex)
{
   Int_t np = gae->GetN();
   for (Int_t i=0; i<np; i++) {
      gae->SetPointEXhigh(i,Ex);
      gae->SetPointEXlow(i,Ex);
   }
}



void MakePaperPlotspPb5023GeV(Bool_t EPOS=kFALSE, Bool_t mT=kFALSE, Bool_t TAPS=kFALSE, Bool_t CGCPi0=kTRUE, Bool_t HIJINGPi0 = kTRUE, Bool_t DPMJetPi0 = kTRUE){    

  //input files
							        //PaperPlots_Tsallis_2017_02_09.root
                                                                //18
  TString fileNamepPbSpectra             = "InputMakePaperPlots/PaperPlots_Tsallis_2017_12_18.root";
  TString fileNamepPbOutput              = "InputMakePaperPlots/PaperPlots_Tsallis_2017_12_18_ALL.root";
  //TString fileNameRpPb                   = "InputMakePaperPlots/PaperPlotsRpPb_2017_11_03.root";
  TString fileNameRpPb                   = "InputMakePaperPlots/PaperPlotsRpPb_2017_12_17.root";
  TString fileNamePeaks                  = "InputMakePaperPlots/PaperPlotsPeaks_2016_11_28.root";
  TString fileNameEPOS3                  = "InputMakePaperPlots/pi0_eta_EPOS3.root";
  TString fileNameTAPS                   = "InputMakePaperPlots/TAPS_eta2pi0.root";
  TString fileNameCGC			 = "ExternalInputpPb/Theory/CGC/TheoryGraphsCGCPi0.root";
  TString fileNameTheory                 = "ExternalInputpPb/Theory/TheoryCompilationPPb.root";
  
  
  
  
  
   
  gSystem->Exec(Form("cp %s %s",fileNamepPbSpectra.Data(),fileNamepPbOutput.Data()));
  
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
  
  TDirectoryFile* directoryTheoryCompilation = (TDirectoryFile*)fileTheoryCompilation->Get("pPb_5.023TeV"); 
  
  if( ! directoryTheoryCompilation ) {
    
    cout<<"The directory pPb_5.023TeV does not exist in the file: "<<fileNameTheory.Data()<<endl;
    
  }
  
  
   
  Double_t pi0PtMin =  0.3;
  Double_t pi0PtMax = 20.0;
  
  Double_t etaPtMin =  0.8;
  Double_t etaPtMax = 20.0;
  
  

  TF1* FitCombPi0                      = (TF1*)filepPbSpectra->Get("Pi0pPb/FitCombPi0");
  TF1* FitCombEta                      = (TF1*)filepPbSpectra->Get("EtapPb/FitCombEta");
  TGraphAsymmErrors* CombPi0Syst  =(TGraphAsymmErrors*)filepPbSpectra->Get("Pi0pPb/graphInvYieldPi0CombpPb5023GeVSystErr_xShifted_NSD");
  TGraphAsymmErrors* CombPi0Stat =(TGraphAsymmErrors*)filepPbSpectra->Get("Pi0pPb/graphInvYieldPi0CombpPb5023GeVStatErr_xShifted_NSD");
  TGraphAsymmErrors* CombEtaSyst  =(TGraphAsymmErrors*)filepPbSpectra->Get("EtapPb/graphInvYieldEtaCombpPb5023GeVSystErr_xShifted_NSD");
  TGraphAsymmErrors* CombEtaStat =(TGraphAsymmErrors*)filepPbSpectra->Get("EtapPb/graphInvYieldEtaCombpPb5023GeVStatErr_xShifted_NSD");
  TGraphAsymmErrors* RatioTsallisCombPi0Syst =(TGraphAsymmErrors*)filepPbSpectra->Get("Pi0pPb/RatioTsallisCombPi0Syst");
  TGraphAsymmErrors* RatioTsallisCombPi0Stat =(TGraphAsymmErrors*)filepPbSpectra->Get("Pi0pPb/RatioTsallisCombPi0Stat");
  TGraphAsymmErrors* RatioTsallisCombEtaSyst =(TGraphAsymmErrors*)filepPbSpectra->Get("EtapPb/RatioTsallisCombEtaSyst");
  TGraphAsymmErrors* RatioTsallisCombEtaStat =(TGraphAsymmErrors*)filepPbSpectra->Get("EtapPb/RatioTsallisCombEtaStat");
  TGraphAsymmErrors* RatioTsallisPCMSyst =(TGraphAsymmErrors*)filepPbSpectra->Get("Pi0pPb/RatioTsallisPCMSyst");
  TGraphAsymmErrors* RatioTsallisPCMStat =(TGraphAsymmErrors*)filepPbSpectra->Get("Pi0pPb/RatioTsallisPCMStat");
  TGraphAsymmErrors* RatioTsallisDalitzSyst =(TGraphAsymmErrors*)filepPbSpectra->Get("Pi0pPb/RatioTsallisDalitzSyst");
  TGraphAsymmErrors*  RatioTsallisDalitzStat=(TGraphAsymmErrors*)filepPbSpectra->Get("Pi0pPb/RatioTsallisDalitzStat");
  TGraphAsymmErrors* RatioTsallisEMCalSyst =(TGraphAsymmErrors*)filepPbSpectra->Get("Pi0pPb/RatioTsallisEMCalSyst");
  TGraphAsymmErrors* RatioTsallisEMCalStat =(TGraphAsymmErrors*)filepPbSpectra->Get("Pi0pPb/RatioTsallisEMCalStat");
  TGraphAsymmErrors*RatioTsallisPHOSSyst  =(TGraphAsymmErrors*)filepPbSpectra->Get("Pi0pPb/RatioTsallisPHOSSyst");
  TGraphAsymmErrors* RatioTsallisPHOSStat =(TGraphAsymmErrors*)filepPbSpectra->Get("Pi0pPb/RatioTsallisPHOSStat");
  TGraphAsymmErrors* RatioTsallisPCMEMCalSyst =(TGraphAsymmErrors*)filepPbSpectra->Get("Pi0pPb/RatioTsallisPCMEMCalSyst");
  TGraphAsymmErrors* RatioTsallisPCMEMCalStat =(TGraphAsymmErrors*)filepPbSpectra->Get("Pi0pPb/RatioTsallisPCMEMCalStat");
 
  TGraphAsymmErrors* RatioEtaTsallisPCMSyst =(TGraphAsymmErrors*)filepPbSpectra->Get("EtapPb/RatioEtaTsallisPCMSyst");
  TGraphAsymmErrors* RatioEtaTsallisPCMStat =(TGraphAsymmErrors*)filepPbSpectra->Get("EtapPb/RatioEtaTsallisPCMStat");
  TGraphAsymmErrors* RatioEtaTsallisEMCalSyst =(TGraphAsymmErrors*)filepPbSpectra->Get("EtapPb/RatioEtaTsallisEMCalSyst");
  TGraphAsymmErrors* RatioEtaTsallisEMCalStat = (TGraphAsymmErrors*)filepPbSpectra->Get("EtapPb/RatioEtaTsallisEMCalStat");
  TGraphAsymmErrors* RatioEtaTsallisPCMEMCalSyst =(TGraphAsymmErrors*)filepPbSpectra->Get("EtapPb/RatioEtaTsallisPCMEMCalSyst");
  TGraphAsymmErrors* RatioEtaTsallisPCMEMCalStat =(TGraphAsymmErrors*)filepPbSpectra->Get("EtapPb/RatioEtaTsallisPCMEMCalStat");
  
  
  
  //Theory calculations
  
   TH1F* histoDPMJetPi0                                = (TH1F*) directoryTheoryCompilation->Get("histoPi0SpecDPMJet5023GeV_Reb");
   TH1F* histoDPMJetEta                                = (TH1F*) directoryTheoryCompilation->Get("histoEtaSpecDPMJet5023GeV_Reb");
   TH1F* histoDPMJetEtaToPi0                           = (TH1F*) directoryTheoryCompilation->Get("histoEtaToPi0DPMJet5023GeV");
   TH1F* histoHIJINGPi0                                = (TH1F*) directoryTheoryCompilation->Get("histoPi0SpecHIJING5023GeV_Reb");
   TH1F* histoHIJINGEta                                = (TH1F*) directoryTheoryCompilation->Get("histoEtaSpecHIJING5023GeV_Reb");
   TH1F* histoHIJINGEtaToPi0                           = (TH1F*) directoryTheoryCompilation->Get("histoEtaToPi0HIJING5023GeV"); 

   // -AM   Pi0 Vogelsang spectrum
   TGraphAsymmErrors* graphNLOCalcDSS14InvYieldPi05023GeV_nCTEQ =(TGraphAsymmErrors*)directoryTheoryCompilation->Get("graphNLOCalcDSS14InvYieldPi05023GeV_nCTEQ");
   
   while ( graphNLOCalcDSS14InvYieldPi05023GeV_nCTEQ->GetX()[graphNLOCalcDSS14InvYieldPi05023GeV_nCTEQ->GetN()-1] > pi0PtMax )
       graphNLOCalcDSS14InvYieldPi05023GeV_nCTEQ->RemovePoint(graphNLOCalcDSS14InvYieldPi05023GeV_nCTEQ->GetN()-1);
   
   
   
   
   TGraph * graphNLOCalcDSS14InvYieldPi0MuTwo5023GeV_nCTEQ =(TGraph*)directoryTheoryCompilation->Get("graphNLOCalcDSS14InvYieldPi0MuTwo5023GeV_nCTEQ");
   TGraph * graphNLOCalcDSS14InvYieldPi0MuOne5023GeV_nCTEQ =(TGraph*)directoryTheoryCompilation->Get("graphNLOCalcDSS14InvYieldPi0MuOne5023GeV_nCTEQ");
   TGraph * graphNLOCalcDSS14InvYieldPi0MuHalf5023GeV_nCTEQ =(TGraph*)directoryTheoryCompilation->Get("graphNLOCalcDSS14InvYieldPi0MuHalf5023GeV_nCTEQ");

   // -AM   Eta Vogelsang spectrum
   TGraphAsymmErrors* graphNLOCalcAESSSInvYieldEta5023GeV_nCTEQ =(TGraphAsymmErrors*)directoryTheoryCompilation->Get("graphNLOCalcAESSSInvYieldEta5023GeV_nCTEQ");
   
    while ( graphNLOCalcAESSSInvYieldEta5023GeV_nCTEQ->GetX()[graphNLOCalcAESSSInvYieldEta5023GeV_nCTEQ->GetN()-1] >  etaPtMax )
       graphNLOCalcAESSSInvYieldEta5023GeV_nCTEQ->RemovePoint(graphNLOCalcAESSSInvYieldEta5023GeV_nCTEQ->GetN()-1);
   
   
   TGraph * graphNLOCalcAESSSInvYieldEtaMuTwo5023GeV_nCTEQ =(TGraph*)directoryTheoryCompilation->Get("graphNLOCalcAESSSInvYieldEtaMuTwo5023GeV_nCTEQ");
   TGraph * graphNLOCalcAESSSInvYieldEtaMuOne5023GeV_nCTEQ =(TGraph*)directoryTheoryCompilation->Get("graphNLOCalcAESSS14InvYieldEtaMuOne5023GeV_nCTEQ");
   TGraph * graphNLOCalcAESSSInvYieldEtaMuHalf5023GeV_nCTEQ =(TGraph*)directoryTheoryCompilation->Get("graphNLOCalcAESSS14InvYieldEtaMuHalf5023GeV_nCTEQ");

   cout<< " First addition AM reading NLO lines for pi0 and eta "<< endl;

   // -AM   Pi0 Vogelsang R_pA pi0
   TGraphAsymmErrors*	graphNLOCalcDSS14RpAPi05023GeV_nCTEQ=(TGraphAsymmErrors*)directoryTheoryCompilation->Get("graphNLOCalcDSS14RpAPi05023GeV_nCTEQ");
   TGraphAsymmErrors*	graphNLOCalcDSS14RpAPi05023GeV_nCTEQ_onlypPbErrs=(TGraphAsymmErrors*)directoryTheoryCompilation->Get("graphNLOCalcDSS14RpAPi05023GeV_nCTEQ_onlypPbErrs");
   TGraphAsymmErrors*	graphNLOCalcDSS14RpAPi05023GeV_nCTEQ_SepCalc=(TGraphAsymmErrors*)directoryTheoryCompilation->Get("graphNLOCalcDSS14RpAPi05023GeV_nCTEQ_SepCalc");
   TGraph*	graphNLOCalcDSS14RpAPi05023GeV_muOne_nCTEQ  = (TGraph*)directoryTheoryCompilation->Get("graphNLOCalcDSS14RpAPi05023GeV_muOne_nCTEQ");
   TGraph*	graphNLOCalcDSS14RpAPi05023GeV_muHalf_nCTEQ = (TGraph*)directoryTheoryCompilation->Get("graphNLOCalcDSS14RpAPi05023GeV_muHalf_nCTE");
   TGraph*	graphNLOCalcDSS14RpAPi05023GeV_muTwo_nCTEQ  = (TGraph*)directoryTheoryCompilation->Get("graphNLOCalcDSS14RpAPi05023GeV_muTwo_nCTEQ");
   

   // -AM   Pi0 Vogelsang R_pA Eta
   TGraphAsymmErrors*	graphNLOCalcAESSSRpAEta5023GeV_nCTEQ=(TGraphAsymmErrors*)directoryTheoryCompilation->Get("graphNLOCalcAESSSRpAEta5023GeV_nCTEQ");
   TGraphAsymmErrors*	graphNLOCalcAESSSRpAEta5023GeV_nCTEQ_onlypPbErrs=(TGraphAsymmErrors*)directoryTheoryCompilation->Get("graphNLOCalcAESSSRpAEta5023GeV_nCTEQ_onlypPbErrs");
   TGraphAsymmErrors*	graphNLOCalcAESSSRpAEta5023GeV_nCTEQ_SepCalc=(TGraphAsymmErrors*)directoryTheoryCompilation->Get("graphNLOCalcAESSSRpAEta5023GeV_nCTEQ_SepCalc");
   TGraph*	graphNLOCalcAESSSRpAEta5023GeV_muOne_nCTEQ  = (TGraph*) directoryTheoryCompilation->Get("graphNLOCalcAESSSRpAEta5023GeV_muOne_nCTEQ");  
   TGraph*	graphNLOCalcAESSSRpAEta5023GeV_muHalf_nCTEQ = (TGraph*) directoryTheoryCompilation->Get("graphNLOCalcAESSSRpAEta5023GeV_muHalf_nCTEQ");
   TGraph*	graphNLOCalcAESSSRpAEta5023GeV_muTwo_nCTEQ  = (TGraph*) directoryTheoryCompilation->Get("graphNLOCalcAESSSRpAEta5023GeV_muTwo_nCTEQ");
 
   cout<< " First addition AM reading NLO lines for R_pPb pi0 and eta "<< endl;



  ////
  
    TH1D* histoRatioPi0DPMJetToFit                      = (TH1D*) histoDPMJetPi0->Clone("histoRatioPi0DPMJetToFit"); 
    histoRatioPi0DPMJetToFit                            = CalculateHistoRatioToFit (histoRatioPi0DPMJetToFit, FitCombPi0); 
    histoRatioPi0DPMJetToFit->GetXaxis()->SetRangeUser(pi0PtMin,pi0PtMax);
    TH1D* histoRatioPi0HIJINGToFit                      = (TH1D*) histoHIJINGPi0->Clone("histoRatioPi0HIJINGToFit"); 
    histoRatioPi0HIJINGToFit                            = CalculateHistoRatioToFit (histoRatioPi0HIJINGToFit, FitCombPi0); 
    histoRatioPi0HIJINGToFit->GetXaxis()->SetRangeUser(pi0PtMin,pi0PtMax);
    
    TH1D* histoRatioEtaDPMJetToFit                      = (TH1D*) histoDPMJetEta->Clone("histoRatioEtaDPMJetToFit"); 
    histoRatioEtaDPMJetToFit                            = CalculateHistoRatioToFit (histoRatioEtaDPMJetToFit, FitCombEta); 
    histoRatioEtaDPMJetToFit->GetXaxis()->SetRangeUser(0.6,20);
    TH1D* histoRatioEtaHIJINGToFit                      = (TH1D*) histoHIJINGEta->Clone("histoRatioEtaHIJINGToFit"); 
    histoRatioEtaHIJINGToFit                            = CalculateHistoRatioToFit (histoRatioEtaHIJINGToFit, FitCombEta); 
    histoRatioEtaHIJINGToFit->GetXaxis()->SetRangeUser(0.6,20);
 


  
  
  TF1*  mTScaling=(TF1*)filepPbSpectra->Get("EtapPb/mTScaling");
  TF1*  EtaSpectrum_mT_pPb=(TF1*)filepPbSpectra->Get("EtapPb/EtaSpectrum_mT_pPb");
  TGraphAsymmErrors* EtaPi07TeVStat =(TGraphAsymmErrors*)filepPbSpectra->Get("EtapPb/EtaPi07TeVStat");
  TGraphAsymmErrors* EtaPi07TeVSyst =(TGraphAsymmErrors*)filepPbSpectra->Get("EtapPb/EtaPi07TeVSyst");
  TGraphAsymmErrors* EtaPi0pPbStat =(TGraphAsymmErrors*)filepPbSpectra->Get("EtapPb/EtaPi0pPbStat");
  TGraphAsymmErrors* EtaPi0pPbSyst =(TGraphAsymmErrors*)filepPbSpectra->Get("EtapPb/EtaPi0pPbSyst");
  TGraphAsymmErrors* EtaPi07TeVStat_mT =(TGraphAsymmErrors*)filepPbSpectra->Get("EtaPi0Ratio_vsmT_Stat_pp7TeV");
  TGraphAsymmErrors* EtaPi07TeVSyst_mT =(TGraphAsymmErrors*)filepPbSpectra->Get("EtaPi0Ratio_vsmT_Sys_pp7TeV");
  TGraphAsymmErrors* EtaPi0pPbStat_mT =(TGraphAsymmErrors*)filepPbSpectra->Get("EtapPb/EtaPi0Ratio_vsmT_Stat");
  TGraphAsymmErrors* EtaPi0pPbSyst_mT =(TGraphAsymmErrors*)filepPbSpectra->Get("EtapPb/EtaPi0Ratio_vsmT_Sys");


 TGraphAsymmErrors* graphEtaPi0pPbTot = AddErrorsOfGraphsQuadratically(EtaPi0pPbStat,EtaPi0pPbSyst);

  //RpPb
  TGraphAsymmErrors*CombinedPi0RpPbSystErr=(TGraphAsymmErrors*)fileRpPb->Get("CombinedPi0RpPbSystErr");
  TGraphAsymmErrors*CombinedPi0RpPbStatErr=(TGraphAsymmErrors*)fileRpPb->Get("CombinedPi0RpPbStatErr");
  
  TGraphAsymmErrors*CombinedEtaRpPbSystErr=(TGraphAsymmErrors*)fileRpPb->Get("CombinedEtaRpPbSystErr");
  TGraphAsymmErrors*CombinedEtaRpPbStatErr=(TGraphAsymmErrors*)fileRpPb->Get("CombinedEtaRpPbStatErr");
  

  
  //PPRef
  TGraphAsymmErrors* graphCombPi0InvCrossSectionSystErrPPRef=(TGraphAsymmErrors*)fileRpPb->Get("Pi0PPReferenceSystErr");
  TGraphAsymmErrors* graphCombPi0InvCrossSectionStatErrPPRef=(TGraphAsymmErrors*)fileRpPb->Get("Pi0PPReferenceStatErr");
  TF1* fitTsallisPi0ppRef5023GeVPt = (TF1*)fileRpPb->Get("Pi0PPReferenceTSallisFit");
  TF1* fitTCMPi0ppRef5023GeVPt     = (TF1*)fileRpPb->Get("Pi0PPReferenceTCMFit");
  
  
  
  
  
  
  
  
  TGraphAsymmErrors* graphRatioEtappRefTsallisFitStatErr = (TGraphAsymmErrors*)fileRpPb->Get("graphRatioEtappRefTsallisFitStatErr");   
  TGraphAsymmErrors* graphRatioEtappRefTsallisFitSystErr = (TGraphAsymmErrors*)fileRpPb->Get("graphRatioEtappRefTsallisFitSystErr");   
 
  
  TGraphAsymmErrors* graphRatioEtappRefTCMFitStatErr = (TGraphAsymmErrors*)fileRpPb->Get("graphRatioEtappRefTCMFitStatErr");   
  TGraphAsymmErrors* graphRatioEtappRefTCMFitSystErr = (TGraphAsymmErrors*)fileRpPb->Get("graphRatioEtappRefTCMFitSystErr");   
 
  TGraphAsymmErrors* graphRatioPi0ppRefTCMFitStatErr = (TGraphAsymmErrors*)fileRpPb->Get("graphRatioPi0ppRefTCMFitStatErr");     
  TGraphAsymmErrors* graphRatioPi0ppRefTCMFitSystErr = (TGraphAsymmErrors*)fileRpPb->Get("graphRatioPi0ppRefTCMFitSystErr"); 
  
  TGraphAsymmErrors* graphRatioPi0ppRefTsallisFitStatErr = (TGraphAsymmErrors*)fileRpPb->Get("graphRatioPi0ppRefTsallisFitStatErr");     
  TGraphAsymmErrors* graphRatioPi0ppRefTsallisFitSystErr = (TGraphAsymmErrors*)fileRpPb->Get("graphRatioPi0ppRefTsallisFitSystErr"); 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  TGraphAsymmErrors* graphCombEtaInvCrossSectionSystErrPPRef=(TGraphAsymmErrors*)fileRpPb->Get("EtaPPReferenceSystErr");
  TGraphAsymmErrors* graphCombEtaInvCrossSectionStatErrPPRef=(TGraphAsymmErrors*)fileRpPb->Get("EtaPPReferenceStatErr");
  TF1* fitTsallisEtappRef5023GeVPt = (TF1*)fileRpPb->Get("EtaPPReferenceTSallisFit");
  TF1* fitTCMEtappRef5023GeVPt     = (TF1*)fileRpPb->Get("EtaPPReferenceTCMFit");
  
  
  
  TGraph*EPS09s_KKP_NLO=(TGraph*)fileRpPb->Get("EPS09s_KKP_NLO");
   TGraph*	EPS09s_AKK_NLO=(TGraph*)fileRpPb->Get("EPS09s_AKK_NLO");
   TGraph*	EPS09s_fDSS_NLO=(TGraph*)fileRpPb->Get("EPS09s_fDSS_NLO");
   TGraphAsymmErrors*	EPS09s_fDSS_errors=(TGraphAsymmErrors*)fileRpPb->Get("EPS09s_fDSS_errors");
   TGraph*	CGC=(TGraph*)fileRpPb->Get("CGC");
   TGraph*      CGCline = (TGraph*)CGC->Clone("CGCline");
   
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
 
 

  TGraphErrors* graphErrMcGillTheoryPion_p_hydro 	= (TGraphErrors*)directoryTheoryCompilation->Get("graphPi0SpecMcGill5023GeV");
  graphErrMcGillTheoryPion_p_hydro->Print();
  TGraphErrors* graphErrMcGillTheoryEta_p_hydro  	= (TGraphErrors*)directoryTheoryCompilation->Get("graphEtaSpecMcGill5023GeV");
  graphErrMcGillTheoryEta_p_hydro->Print();
  
  TGraphErrors* graphMcGillEtaToPi0                      = (TGraphErrors*) directoryTheoryCompilation->Get("graphEtaToPi0McGill5023GeV");
   
  
  //read Ilkka data
 TDirectoryFile* directoryIlkka = (TDirectoryFile*)fileCGC->Get("Ilkka"); 

 TGraphAsymmErrors* graphAsymmErrIlkkapPb5020_pi0_ct14_epps16_dss14_scale_err = (TGraphAsymmErrors*)directoryIlkka->Get("graphAsymmErr_pi0_ct14_epps16_dss14_scale_sumerr");
 TGraphAsymmErrors* graphAsymmErr_pi0_ct14_epps16_dss14_sumerr                = (TGraphAsymmErrors*)directoryIlkka->Get("graphAsymmErr_pi0_ct14_epps16_dss14_sumerr");
 TGraphAsymmErrors* graphAsymmErr_pi0_ct14_epps16_dss14 		      = (TGraphAsymmErrors*)directoryIlkka->Get("graphAsymmErr_pi0_ct14_epps16_dss14");
 
 // read RpPb data
 
 TDirectoryFile* directoryRpPb = (TDirectoryFile*)fileCGC->Get("RpPb"); 

 TGraphAsymmErrors* graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14 = (TGraphAsymmErrors*)directoryRpPb->Get("graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14");
 
 /////////////////////read eta2pi0_d-Au of PHENIX
 
 TGraphErrors* graphEta2Pi0dAuPhenix = Eta2Pi0dAu();
 
 
  
 
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
  
  Double_t CMarginEtaToPi0L=0.12;
  Double_t CMarginEtaToPi0T=0.02;
  Double_t CMarginEtaToPi0R=0.02;
  Double_t CMarginEtaToPi0B=0.12;
  
  
  
  

  Int_t CDimX=500;
  Int_t CDimY=400;

  Double_t TitleOffsetX=1.;
  Double_t TitleOffsetY=1.;
  Double_t TitleSizeX=0.05;
  Double_t TitleSizeY=0.05;
  Double_t TitleSizeYEtaToPi0 = 0.07;
  Double_t TitleOffsetYEtaToPi0 = 0.78;

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
  Color_t colorPPRef            = kGray+2;
  Color_t  colorDPMJet          = kViolet+2;
  Color_t  colorDPMJetPi0       = kViolet+2;
  Color_t  colorDPMJetEta       = kViolet+2; //-2
  
  Color_t  colorHIJING          = kGreen-2;
  Color_t  colorHIJINGPi0       = kGreen-2;
  Color_t  colorHIJINGEta       = kGreen-2; //+2
  
  
  Color_t  colorMcGill	        = kPink + 2;
  Color_t  colorMcGillline      = kPink + 3;
  
  Color_t  colorEPOS            = kAzure+7;
  Color_t  colorEPOSline        = kAzure+8;
  Color_t  colorEPOSPi0         = kAzure+7;
  Color_t  colorEPOSEta		= kAzure+7;
  
  Color_t  colorMcGillPi0	= kPink + 2;
  Color_t  colorMcGillEta	= kPink + 2; //-6
  
  Color_t  colorIlkka           = kGray;// kYellow - 10;//41;
  Color_t  colorIlkkaline       = kYellow + 2 ;
  Color_t  colorIlkkaPi0	= kGray;//kYellow - ;//41;
  Color_t  colorEtaPi07TeV	= kBlue - 3;
  
  Color_t colornCTEQPi0Line     = kBlue-5;
  Color_t colornCTEQEtaLine     = kBlue-5;
  Color_t colorDSSnPDFEPPSBand               = kRed-2;
  Width_t lineWidthDSSnPDFEPPSLine = 2;
  
  
   
   
  


  Style_t  styleLineDPMJet                    = 7;
  Style_t  styleLineHIJING                    = 8;
  Style_t  styleLineCGC			      = 1;
  Style_t  styleLineEPOS3		      = 1;
  Style_t  styleLineMcGill		      = 1;
  Style_t  styleLineIlkka		      = 3;
  Style_t  styleLineIlkkaRatio                = 3;
  Style_t  styleLineTsallis                   = 1;
  Style_t  styleLineTCM                       = 2;
  Style_t  styleLineDSSnPDFEPPS               = 2;
  
  Int_t     lineWidthMcGill  = 2;
  Int_t     lineWidthIlkka   = 2;
  Int_t     lineWidthEPOS3   = 2;
  Int_t     lineWidthTsallis = 2;
  Int_t     lineWidthTCM     = 2;

  
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
  Width_t  widthLinesBoxes                    = 1.4;
  



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
  FitCombPi0->SetLineWidth(lineWidthTsallis);  
  if (!EPOS)  FitCombPi0->Draw("same"); 
  FitCombEta->SetLineWidth(lineWidthTsallis);  
  if (!EPOS) FitCombEta->Draw("same");

  TGraphAsymmErrors*ppRefStat_scale=(TGraphAsymmErrors*)ScaleGraph(ppRefStat,(0.0983e-9));
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
  TGraphAsymmErrors* CombPi0Stat_noXerrors = (TGraphAsymmErrors*)CombPi0Stat->Clone();
  SetEx(CombPi0Stat_noXerrors,0.);
  DrawGammaSetMarkerTGraphAsym(CombPi0Stat_noXerrors,markerStylePi0,1,colorCombYieldPi0, colorCombYieldPi0 ); 
 
  
  DrawGammaSetMarkerTGraphAsym(CombEtaSyst,markerStyleEta,1,colorCombYieldEta ,colorCombYieldEta, 1, kTRUE);
  CombEtaSyst->Draw("E2,same");
  DrawGammaSetMarkerTGraphAsym(CombEtaStat,markerStyleEta,1,colorCombYieldEta, colorCombYieldEta ); 
  CombEtaStat->Draw("pz,same"); 
  TGraphAsymmErrors* CombEtaStat_noXerrors = (TGraphAsymmErrors*)CombEtaStat->Clone();
  SetEx(CombEtaStat_noXerrors,0.);
  DrawGammaSetMarkerTGraphAsym(CombEtaStat_noXerrors,markerStyleEta,1,colorCombYieldEta, colorCombYieldEta ); 
 
  
  

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

  // -AM
  // nCTEQ- Vogelsang

  cout<< " second addition AM : Ratio to Fit  for pi0  "<< endl;
  TGraphAsymmErrors* RatioDSS14nCTEQFitPi05023GeV = (TGraphAsymmErrors*)graphNLOCalcDSS14InvYieldPi05023GeV_nCTEQ->Clone("RatioDSS14nCTEQFitPi05023GeV");
  RatioDSS14nCTEQFitPi05023GeV=(TGraphAsymmErrors*)CalculateGraphErrRatioToFit(RatioDSS14nCTEQFitPi05023GeV,FitCombPi0);



  cout<< " second addition AM : Ratio to Fit  for eta "<< endl;  
  TGraphAsymmErrors* RatioAESSSnCTEQFitEta5023GeV = (TGraphAsymmErrors*)graphNLOCalcAESSSInvYieldEta5023GeV_nCTEQ->Clone("RatioAESSSnCTEQFitEta5023GeV");
  RatioAESSSnCTEQFitEta5023GeV=(TGraphAsymmErrors*)CalculateGraphErrRatioToFit(RatioAESSSnCTEQFitEta5023GeV,FitCombEta);



  
 
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
       
       
     RatioEPOSFitPi0->SetLineWidth(lineWidthEPOS3);
     RatioEPOSFitPi0->SetLineStyle(styleLineEPOS3);
     RatioEPOSFitPi0->SetLineColor(colorEPOSline);
     RatioEPOSFitPi0->SetFillColor(colorEPOS);
     
     
     RatioEPOSFitEta->SetLineWidth(lineWidthEPOS3);
     RatioEPOSFitEta->SetLineStyle(styleLineEPOS3);
     RatioEPOSFitEta->SetLineColor(colorEPOSline);
     RatioEPOSFitEta->SetFillColor(colorEPOS);
     
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
   DrawGammaSetMarkerTGraphAsym(RatioTsallisCombPi0Syst,markerStylePi0,1.2,colorCombYieldPi0 ,colorCombYieldPi0 , 1, kTRUE);  
   RatioTsallisCombPi0Syst->Draw("E2same");
   DrawGammaSetMarkerTGraphAsym(RatioTsallisCombPi0Stat,markerStylePi0,1.2,colorCombYieldPi0, colorCombYieldPi0 );  
   SetEx(RatioTsallisCombPi0Stat,0.); //Set the x-errors bars to 0
   RatioTsallisCombPi0Stat->Draw("Ez,p,same");
   
  
  
   DrawGammaSetMarkerTGraphAsym(RatioTsallisCombEtaSyst,markerStyleEta,1.2,colorCombYieldEta ,colorCombYieldEta, 1, kTRUE);  
   RatioTsallisCombEtaSyst->Draw("E2same");
   DrawGammaSetMarkerTGraphAsym(RatioTsallisCombEtaStat,markerStyleEta,1.2,colorCombYieldEta ,colorCombYieldEta );  
   SetEx(RatioTsallisCombEtaStat,0.); //Set the x-errors bars to 0
   RatioTsallisCombEtaStat->Draw("Ez,p,same");
   
  c1->Update();
  if(EPOS) c1->Print(Form("%s/MesonYields_EPOS.%s",outputDir.Data(),suffix.Data()));
  else if(mT) c1->Print(Form("%s/MesonYields_mT.%s",outputDir.Data(),suffix.Data()));
  else c1->Print(Form("%s/MesonYields.%s",outputDir.Data(),suffix.Data()));
  

  
  
  
  
  
  
  


  //Inv Yields Ratio Pi0
 

   TCanvas* c2 = new TCanvas("c2","",200,10,CDimX,CDimY);  // gives the page size
   DrawGammaCanvasSettings( c2, CMarginL, CMarginR,CMarginT ,CMarginB);
   c2->SetLogx();
   TH2F * hist2;
   hist2 = new TH2F("hist2","hist2",1000,0.25,25.,1000,0.34,2.8  );
   SetStyleHistoTH2ForGraphs(hist2, "#it{p}_{T} (GeV/#it{c})","Data/Fit",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY,TitleOffsetX,TitleOffsetY, 512, 508); 
   hist2->GetXaxis()->SetLabelOffset(LabelOffsetLog);
   hist2->DrawCopy(); 

  
   DrawGammaLines(0.0,25.,1.,1.,2.0,kGray+2,2);
  
  

   DrawGammaSetMarkerTGraphAsym(RatioTsallisPHOSSyst,markerStylePHOS,markerSizePHOS, colorPHOS, colorPHOS, 1, kTRUE);  
   RatioTsallisPHOSSyst->Draw("E2same");
   SetEx(RatioTsallisPHOSStat,0);
   DrawGammaSetMarkerTGraphAsym(RatioTsallisPHOSStat, markerStylePHOS,markerSizePHOS, colorPHOS, colorPHOS );
   RatioTsallisPHOSStat->Draw("same,pz") ;
 
   DrawGammaSetMarkerTGraphAsym(RatioTsallisEMCalSyst,markerStyleEMCal,markerSizeEMCal, colorEMCal, colorEMCal, 1, kTRUE);  
   RatioTsallisEMCalSyst->Draw("E2same");
   SetEx(RatioTsallisEMCalStat,0);
   DrawGammaSetMarkerTGraphAsym(RatioTsallisEMCalStat,markerStyleEMCal,markerSizeEMCal, colorEMCal, colorEMCal );
   RatioTsallisEMCalStat->Draw("same,pz") ; 

   DrawGammaSetMarkerTGraphAsym(RatioTsallisPCMSyst,markerStylePCM,markerSizePCM, colorPCM, colorPCM, 1, kTRUE);  
   RatioTsallisPCMSyst->Draw("E2same");
   SetEx(RatioTsallisPCMStat,0);
   DrawGammaSetMarkerTGraphAsym(RatioTsallisPCMStat,markerStylePCM,markerSizePCM, colorPCM, colorPCM );
   RatioTsallisPCMStat->Draw("same,pz") ;
 
   DrawGammaSetMarkerTGraphAsym(RatioTsallisDalitzSyst,markerStyleDalitz,markerSizeDalitz, colorDalitz, colorDalitz, 1, kTRUE);  
   RatioTsallisDalitzSyst->Draw("E2same");
   SetEx(RatioTsallisDalitzStat,0);
   DrawGammaSetMarkerTGraphAsym(RatioTsallisDalitzStat,markerStyleDalitz,markerSizeDalitz, colorDalitz, colorDalitz );
   RatioTsallisDalitzStat->Draw("same,pz") ; 

  
  DrawGammaSetMarkerTGraphAsym(RatioTsallisPCMEMCalSyst,markerStylePCMEMCal,markerSizePCMEMCal,colorPCMEMCal,colorPCMEMCal,1,kTRUE);
  RatioTsallisPCMEMCalSyst->Draw("E2same");
  SetEx(RatioTsallisPCMEMCalStat,0);
  DrawGammaSetMarkerTGraphAsym(RatioTsallisPCMEMCalStat,markerStylePCMEMCal,markerSizePCMEMCal, colorPCMEMCal, colorPCMEMCal );
  RatioTsallisPCMEMCalStat->Draw("same,pz") ; 
  
  
  
 	
  TLatex * lt2 = new TLatex(1.9,2.61,"ALICE") ;
  lt2->SetTextColor(kBlack) ;
  lt2->SetTextSize(TextSize) ;
  lt2->SetTextFont(Font) ;
  lt2->DrawLatex(1.9,2.45,"p-Pb, NSD, #sqrt{#it{s}_{NN}} = 5.02 TeV, #pi^{0}");
  
 
  
  lt2->Draw() ;
  
  TLegend* leg2 = new TLegend(0.13,0.65,0.45,0.95);
  leg2->SetFillColor(0);
  leg2->SetLineColor(0);
  leg2->SetTextFont(Font);
  leg2->SetTextSize(TextSize);
  leg2->AddEntry(RatioTsallisPHOSSyst,Form("PHOS"),"pef");
  leg2->AddEntry(RatioTsallisEMCalSyst,Form("EMC"),"pef");
  leg2->AddEntry(RatioTsallisPCMSyst,Form("PCM"),"pef");
  leg2->AddEntry(RatioTsallisPCMEMCalSyst,Form("PCM-EMC"),"pef");
  leg2->AddEntry(RatioTsallisDalitzSyst,Form("PCM-#gamma*#gamma"),"pef");
  leg2->Draw("same");





  c2->Update();
  c2->Print(Form("%s/RatioCombFitIndividualPi0.%s",outputDir.Data(),suffix.Data()));



  //Inv Yields Ratio Eta
 


  TCanvas* c3 = new TCanvas("c3","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c3, CMarginL, CMarginR,CMarginT ,CMarginB);
  c3->SetLogx();
  TH2F * hist3;
  hist3 = new TH2F("hist3","hist3",1000,0.60,25.,1000,0.34,2.8);
  SetStyleHistoTH2ForGraphs(hist3, "#it{p}_{T} (GeV/#it{c})","Data/Fit",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY, TitleOffsetX,TitleOffsetY, 512, 508); 
  hist3->GetXaxis()->SetLabelOffset(LabelOffsetLog);
  hist3->DrawCopy(); 

  DrawGammaLines(0.0,25.,1.,1.,2.0,kGray+2,2);
  
  DrawGammaSetMarkerTGraphAsym(RatioEtaTsallisEMCalSyst,markerStyleEMCal,markerSizeEMCal, colorEMCal, colorEMCal, 1, kTRUE);
  SetEx(RatioEtaTsallisEMCalStat,0);  
  DrawGammaSetMarkerTGraphAsym(RatioEtaTsallisEMCalStat,markerStyleEMCal,markerSizeEMCal, colorEMCal, colorEMCal );
  RatioEtaTsallisEMCalSyst->Draw("E2same");
  RatioEtaTsallisEMCalStat->Draw("pz,same") ;

  DrawGammaSetMarkerTGraphAsym(RatioEtaTsallisPCMSyst,markerStylePCM,markerSizePCM, colorPCM, colorPCM, 1, kTRUE);
  SetEx(RatioEtaTsallisPCMStat,0);    
  DrawGammaSetMarkerTGraphAsym(RatioEtaTsallisPCMStat,markerStylePCM,markerSizePCM, colorPCM, colorPCM ); 
  RatioEtaTsallisPCMSyst->Draw("E2same");
  RatioEtaTsallisPCMStat->Draw("pz,same") ;	
   
  DrawGammaSetMarkerTGraphAsym(RatioEtaTsallisPCMEMCalSyst,markerStylePCMEMCal,markerSizePCMEMCal, colorPCMEMCal, colorPCMEMCal, 1, kTRUE);  
  SetEx(RatioEtaTsallisPCMEMCalStat,0); 
  DrawGammaSetMarkerTGraphAsym(RatioEtaTsallisPCMEMCalStat,markerStylePCMEMCal,markerSizePCMEMCal, colorPCMEMCal, colorPCMEMCal); 
  RatioEtaTsallisPCMEMCalSyst->Draw("E2same");
  RatioEtaTsallisPCMEMCalStat->Draw("pz,same") ;
  
  RatioEtaTsallisPCMEMCalStat->Print();
  return;


  TLegend* leg3 = new TLegend(0.13,0.8,0.45,0.95);
  leg3->SetFillColor(0);
  leg3->SetLineColor(0);
  leg3->SetTextFont(Font);
  leg3->SetTextSize(TextSize);
  leg3->AddEntry(RatioEtaTsallisEMCalSyst,Form("EMC"),"pef");
  leg3->AddEntry(RatioEtaTsallisPCMSyst,Form("PCM"),"pef");
  leg3->AddEntry(RatioEtaTsallisPCMEMCalSyst,Form("PCM-EMC"),"pef");
  leg3->Draw("same");
  
  TLatex * lt3 = new TLatex(3.1,2.61,"ALICE");
  lt3->SetTextColor(kBlack) ;
  lt3->SetTextSize(TextSize) ;
  lt3->SetTextFont(Font) ;
  lt3->DrawLatex(3.1,2.45,"p-Pb, NSD, #sqrt{#it{s}_{NN}} = 5.02 TeV, #eta");
  lt3->Draw() ;


  c3->Update();
  c3->Print(Form("%s/RatioCombFitIndividualEta.%s",outputDir.Data(),suffix.Data()));





  //Eta/Pi0 Ratio

  

  TCanvas* c4 = new TCanvas("c4","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c4, CMarginEtaToPi0L, CMarginEtaToPi0R,CMarginEtaToPi0T ,CMarginEtaToPi0B);

  TH2F * hist4;
  if (TAPS)  hist4 = new TH2F("hist4","hist4",1000,-0.2,16.,1000,-0.02,.99   );
  else  hist4 = new TH2F("hist4","hist4",1000,0.3,16.,1000,0.0,.99   );

  SetStyleHistoTH2ForGraphs(hist4, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0} ",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeYEtaToPi0, TitleOffsetX,TitleOffsetYEtaToPi0, 512, 508);
  hist4->DrawCopy(); 
  if (EPOS){
      
    graphEtaPi0EPOS->SetLineColor(kGreen+1);
    graphEtaPi0EPOS->SetFillColor(kGreen+1);
    graphEtaPi0EPOS->SetLineWidth(2);
    graphEtaPi0EPOS->Draw("C3same");

  }      
  
 // if (!EPOS || (EPOS && TAPS)) mTScaling->Draw("same"); TEMP
 // graphEtaPi0Ratio_mTScaled_pPb->SetMarkerColor(kRed+2);
 // graphEtaPi0Ratio_mTScaled_pPb->SetLineColor(1);
 //  graphEtaPi0Ratio_mTScaled_pPb->SetFillColor(kRed+2);

   if (TAPS){
       
     //DrawGammaSetMarkerTGraphErr(eta2pi0_pAu, 27, 1.4,kGreen+1, kGreen+1);
     //DrawGammaSetMarkerTGraphErr(eta2pi0_pAu, 27, 1.4,colorCombYieldEta+1, colorCombYieldEta+1);
     DrawGammaSetMarkerTGraphErr(eta2pi0_pAu, 27, 1.4,kBlack,kBlack);
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


  

  TCanvas* c4a = new TCanvas("c4a","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c4a, CMarginEtaToPi0L, CMarginEtaToPi0R,CMarginEtaToPi0T ,CMarginEtaToPi0B);
  c4a->SetLogx(); 
  // c4a->SetLogy();
  TH2F * hist4a;
  if (TAPS)   hist4a = new TH2F("hist4a","hist4a",1000,0.101,25.,1000,-0.05,1.2   );
  else hist4a = new TH2F("hist4a","hist4a",1000,0.3,25.,1000,0.0,.99   );
  SetStyleHistoTH2ForGraphs(hist4a, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0} ",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeYEtaToPi0, TitleOffsetX,TitleOffsetYEtaToPi0, 512, 508);
  hist4a->GetYaxis()->SetLabelOffset(0.01);
  hist4a->GetXaxis()->SetLabelOffset(LabelOffsetLog); 
  hist4a->DrawCopy(); 
  
   graphMcGillEtaToPi0->SetLineStyle(styleLineMcGill);
   graphMcGillEtaToPi0->SetLineColor(colorMcGillline);
   graphMcGillEtaToPi0->SetFillColor(colorMcGill);
   graphMcGillEtaToPi0->SetLineWidth(lineWidthMcGill);
   graphMcGillEtaToPi0->Draw("C3same");
   graphMcGillEtaToPi0->Draw("lXZsame");

   if (EPOS){
    graphEtaPi0EPOS->SetLineColor(colorEPOSline);
    graphEtaPi0EPOS->SetFillColor(colorEPOS);
    graphEtaPi0EPOS->SetLineWidth(3);
    graphEtaPi0EPOS->Draw("C3same");
    graphEtaPi0EPOS->Draw("lXZsame");
  }  
  
   if (!EPOS || (EPOS && TAPS)) mTScaling->Draw("same"); 
   
   
   if( graphEta2Pi0dAuPhenix ) {
       
        DrawGammaSetMarkerTGraphErr(graphEta2Pi0dAuPhenix, 24, 1.2,kRed, kRed);
        graphEta2Pi0dAuPhenix->Draw("same,zp");
       
   }
 
   if (TAPS){
     DrawGammaSetMarkerTGraphErr(eta2pi0_pAu, 27, 1.5,colorCombYieldEta,colorCombYieldEta);
     DrawGammaSetMarkerTGraphErr(eta2pi0_pBe, 27, 1.5,kRed, kRed);
     eta2pi0_pBe->Draw("same,zp");
     eta2pi0_pAu->Draw("same,zp");
     //eta2pi0_pBe->Draw("same,pz");
   }
   if( HIJINGPi0 ){
   SetStyleHisto(histoDPMJetEtaToPi0, 3, styleLineDPMJet, colorDPMJet );  
   histoDPMJetEtaToPi0->Draw("same,hist,l");
   }

   if( DPMJetPi0 ){
   SetStyleHisto(histoHIJINGEtaToPi0,3, 8, colorHIJING);  
   histoHIJINGEtaToPi0->Draw("same,hist,l");
   }
   

   
  TGraphAsymmErrors* EtaPi07TeVStat_noXerrors = (TGraphAsymmErrors*)EtaPi07TeVStat->Clone();
  SetEx(EtaPi07TeVStat_noXerrors,0);
   

  if (!EPOS || (EPOS && TAPS)){ 
    
  
 
    
  DrawGammaSetMarkerTGraphAsym(EtaPi07TeVStat_noXerrors, 34, 1, colorEtaPi07TeV, colorEtaPi07TeV);
  EtaPi07TeVStat_noXerrors->Draw("same,zp");
  DrawGammaSetMarkerTGraphAsym(EtaPi07TeVSyst, 34, 1, colorEtaPi07TeV, colorEtaPi07TeV, 1., kTRUE);
  EtaPi07TeVSyst->Draw("same,E2");
  }  
  TGraphAsymmErrors* EtaPi0pPbStat_noXerrors = (TGraphAsymmErrors*)EtaPi0pPbStat->Clone();
  SetEx(EtaPi0pPbStat_noXerrors,0);
  
  DrawGammaSetMarkerTGraphAsym(EtaPi0pPbStat_noXerrors,20,1,1,1);  
  EtaPi0pPbStat_noXerrors->Draw("pz,same");
  DrawGammaSetMarkerTGraphAsym(EtaPi0pPbSyst,20,1, 1,1, 1, kTRUE);  
  EtaPi0pPbSyst->Draw("E2,same");
  

  // -AM   Fit with a constant the high pT part of the eta/pi0 ratio


  TF1 * fitEtaPi0ConstStat = new TF1("fitEtaPi0ConstStat","[0]",4,20);
  TF1 * fitEtaPi0ConstTot = new TF1("fitEtaPi0ConstTot","[0]",4,20); 

  EtaPi0pPbStat->Fit("fitEtaPi0ConstStat","QRME0","",4,20);
  graphEtaPi0pPbTot->Fit("fitEtaPi0ConstTot","QRME0","",4,20);
  cout<< "Stat"<< endl;
  EtaPi0pPbStat->Print();
  cout<< "Sys"<< endl;
  EtaPi0pPbSyst->Print();
  cout<< "Tot"<< endl;
  graphEtaPi0pPbTot->Print();
  cout<< endl;
  cout<< endl;
  cout<< "////////////////////////////////"<<endl;
  cout<< "Fitting eta/pi0 ratio pt>4 GeV"<< endl;
  cout<< "stat errors:::"<<  fitEtaPi0ConstStat->GetParameter(0) << " +/-"<<   fitEtaPi0ConstStat->GetParError(0)<< endl; 
  cout<< "Tot errors:::" <<  fitEtaPi0ConstTot->GetParameter(0) << " +/-"<<   fitEtaPi0ConstTot->GetParError(0)<< endl; 
  cout<< "Sys errors:::"<< TMath::Power((fitEtaPi0ConstTot->GetParError(0)*fitEtaPi0ConstTot->GetParError(0)-fitEtaPi0ConstStat->GetParError(0)*fitEtaPi0ConstStat->GetParError(0)),0.5);
  cout<< endl;
  cout<< endl;
  cout<< endl;  

  TLatex * lt4a = new TLatex(.12,1.12,"ALICE");
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
  //if (!EPOS || (EPOS && TAPS))leg4a->AddEntry(EtaPi07TeVSyst,"pp, #sqrt{#it{s}} = 7 TeV","pef");// (PLB717 (2012) 162)","pef");  
  if (!EPOS || (EPOS && TAPS)){ 
    leg4a->AddEntry(mTScaling,"#eta from #it{m}_{T} scaled #pi^{0}","pl"); 
  }
  if (EPOS) leg4a->AddEntry(graphEtaPi0EPOS,"EPOS3","fl");
   leg4a->AddEntry(graphMcGillEtaToPi0,"VISHNU","fl");
  if( HIJINGPi0){
    leg4a->AddEntry(histoHIJINGEtaToPi0,"HIJING","l"); 
  }
  
  
  if( DPMJetPi0 ){
     leg4a->AddEntry(histoDPMJetEtaToPi0,"DPMJet","l"); 
  }
   if (!EPOS || (EPOS && TAPS))leg4a->AddEntry(EtaPi07TeVSyst,"pp, #sqrt{#it{s}} = 7 TeV","pef");// (PLB717 (2012) 162)","pef");
  
  
  
  
  
   leg4a->Draw("same");
  
  if (TAPS) {
    
     TLatex * lt4b = new TLatex(0.12,0.53,"CERES-TAPS #sqrt{#it{s}_{NN}} = 29.1 GeV");
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
  
  if( graphEta2Pi0dAuPhenix ) {
      
      
     TLatex * lt4c = new TLatex(2.7,1.12,"PHENIX");
     lt4c->SetTextColor(kBlack);
     lt4c->SetTextSize(TextSize);
     lt4c->SetTextFont(Font);
     lt4c->Draw();
     
     TLegend* leg4c;
     leg4c = new TLegend(0.6,0.85,0.9,0.9);
     leg4c->SetFillColor(0);
     leg4c->SetLineColor(0);
     leg4c->SetTextFont(Font);
     leg4c->SetTextSize(0.041);
     leg4c->AddEntry( graphEta2Pi0dAuPhenix,"d-Au, #sqrt{#it{s}_{NN}} = 200 GeV","p");  
     leg4c->Draw("same");
      
      
  }



  c4a->Update();
  if (EPOS && TAPS)  	c4a->Print(Form("%s/EtaPi0Ratio_LogX_EPOS_TAPS.%s",outputDir.Data(),suffix.Data()));
  else if (EPOS)  	c4a->Print(Form("%s/EtaPi0Ratio_LogX_EPOS.%s",outputDir.Data(),suffix.Data()));
  else if (TAPS)  	c4a->Print(Form("%s/EtaPi0Ratio_LogX_TAPS.%s",outputDir.Data(),suffix.Data()));
  else  		c4a->Print(Form("%s/EtaPi0Ratio_LogX.%s",outputDir.Data(),suffix.Data()));

  

  TCanvas* c4b = new TCanvas("c4b","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c4b, CMarginEtaToPi0L, CMarginEtaToPi0R,CMarginEtaToPi0T ,CMarginEtaToPi0B);

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
  DrawGammaCanvasSettings( c4c, CMarginEtaToPi0L, CMarginEtaToPi0R,CMarginEtaToPi0T ,CMarginEtaToPi0B);
  c4c->SetLogx();
  TH2F * hist4c;
  hist4c = new TH2F("hist4c","hist4c",1000,0.6,20.,1000,0.0,.99   );
  SetStyleHistoTH2ForGraphs(hist4c, "#it{m}_{T} (GeV/#it{c})","#eta/#pi^{0} ",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeYEtaToPi0, TitleOffsetX,TitleOffsetYEtaToPi0, 512, 508);
  hist4c->GetXaxis()->SetLabelOffset(LabelOffsetLog); 
  hist4c->DrawCopy(); 
    

     

  TLegend* leg4c = new TLegend(0.15,0.70,0.45,0.90);
  leg4c->SetFillColor(0);
  leg4c->SetLineColor(0);
  leg4c->SetTextFont(Font);
  leg4c->SetTextSize(0.045);
  leg4c->Draw("same");
  TLine* line4c=new TLine(0.6,0.47,20.,0.47);
  line4c->SetLineColor(kRed+2);
  line4c->SetLineStyle(2);
  line4c->Draw("same");

  c4c->Update();
  c4c->Print(Form("%s/EtaPi0Ratio_mT_LogX.%s",outputDir.Data(),suffix.Data()));

  

  TCanvas* c4d = new TCanvas("c4d","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c4d, CMarginEtaToPi0L+0.02, CMarginEtaToPi0R,CMarginEtaToPi0T ,CMarginEtaToPi0B);
  c4d->SetLogx();
  TH2F * hist4d;
  hist4d = new TH2F("hist4d","hist4d",1000,0.6,25.,1000,0.21,1.39   );
  SetStyleHistoTH2ForGraphs(hist4d, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}_{measured}/#eta/#pi^{0}_{#it{m}_{T} scaling}",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeYEtaToPi0, TitleOffsetX,TitleOffsetYEtaToPi0+0.05, 512, 508);
  hist4d->GetXaxis()->SetLabelOffset(LabelOffsetLog); 
 
  hist4d->DrawCopy(); 
  
  TGraphAsymmErrors* EtaPi0Stat_div_mT=(TGraphAsymmErrors*)EtaPi0pPbStat->Clone("EtaPi0Stat_div_mT");
  TGraphAsymmErrors* EtaPi0Syst_div_mT=(TGraphAsymmErrors*)EtaPi0pPbSyst->Clone("EtaPi0Syst_div_mT");
  EtaPi0Stat_div_mT=CalculateGraphErrRatioToFit (EtaPi0pPbStat ,mTScaling  );  
  EtaPi0Syst_div_mT=CalculateGraphErrRatioToFit (EtaPi0pPbSyst ,mTScaling  );  
   
  TGraphAsymmErrors* EtaPi0Stat_div_mT_noXErrors = (TGraphAsymmErrors*)EtaPi0Stat_div_mT->Clone();
  SetEx(EtaPi0Stat_div_mT_noXErrors,0.);
  DrawGammaSetMarkerTGraphAsym(EtaPi0Stat_div_mT_noXErrors, 20, 1,kBlack, kBlack);
  EtaPi0Stat_div_mT_noXErrors->Draw("same,zp");
   
  DrawGammaSetMarkerTGraphAsym(EtaPi0Syst_div_mT, 20, 1, kBlack, kBlack, 1., kTRUE);
  EtaPi0Syst_div_mT->Draw("same,E2"); //TEMP
   
  TLegend* leg4d = new TLegend(0.45,0.25,0.95,0.35);
  leg4d->SetFillColor(0);
  leg4d->SetLineColor(0);
  leg4d->SetTextFont(Font);
  leg4d->SetTextSize(0.045);
  leg4d->Draw("same");
   
  DrawGammaLines(0.6,25.,1.,1.,2.0,kGray+2,2);
  
  
  TLatex * Prelim4d = new TLatex(0.7,1.3,"ALICE");
  Prelim4d->SetTextColor(kBlack) ;
  Prelim4d->SetTextSize(TextSize) ;
  Prelim4d->SetTextFont(Font) ;
  Prelim4d->DrawLatex(0.7,1.22,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  Prelim4d->Draw() ;

  
  c4d->Update();
  c4d->Print(Form("%s/EtaPi0Ratio2mTScaling.%s",outputDir.Data(),suffix.Data()));


  
  TCanvas* c4e = new TCanvas("c4e","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c4e, CMarginEtaToPi0L, CMarginEtaToPi0R,CMarginEtaToPi0T ,CMarginEtaToPi0B);
  c4e->SetLogx(); 
  
  TH2F* hist4e = new TH2F("hist4e","hist4e",1000,0.101,25.,1000,-0.05,1.2); //0.101,25.,1000,-0.01,1.2
                                                      //          1000,0.09,25.,1000,0.01,1.2
  
  SetStyleHistoTH2ForGraphs(hist4e, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0} ",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeYEtaToPi0, TitleOffsetX,TitleOffsetYEtaToPi0, 512, 508);
  hist4e->GetYaxis()->SetLabelOffset(0.01);
  hist4e->GetXaxis()->SetLabelOffset(LabelOffsetLog); 
  hist4e->DrawCopy(); 

  EtaPi0pPbStat_noXerrors->Draw("pz,same");
  EtaPi0pPbSyst->Draw("E2,same");
  mTScaling->Draw("same"); 
  graphEta2Pi0dAuPhenix->Draw("same,zp");
  eta2pi0_pAu->Draw("same,zp");
  eta2pi0_pBe->Draw("same,pz");
  EtaPi07TeVStat_noXerrors->Draw("same,zp");
  EtaPi07TeVSyst->Draw("same,E2");

  
  
  TLatex * lt4e = new TLatex(.12,1.12,"ALICE");
  lt4e->SetTextColor(kBlack);
  lt4e->SetTextSize(TextSize);
  lt4e->SetTextFont(Font);
  //lt4e->DrawLatex(0.64,0.89,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  lt4e->Draw() ;
  
  
    
  TLegend* leg4e = new TLegend(0.15,0.73,0.45,0.90);
  leg4e->SetFillColor(0);
  leg4e->SetLineColor(0);
  leg4e->SetTextFont(Font);
  leg4e->SetTextSize(0.041);
  leg4e->AddEntry(EtaPi0pPbSyst,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  leg4e->AddEntry(mTScaling,"#eta from #it{m}_{T} scaled #pi^{0}","pl"); 
  leg4e->AddEntry(EtaPi07TeVSyst,"pp, #sqrt{#it{s}} = 7 TeV","pef");// (PLB717 (2012) 162)","pef");
  leg4e->Draw("same");
  
  
  TLatex* lt4ea = new TLatex(0.12,0.63,"CERES-TAPS #sqrt{#it{s}_{NN}} = 29.1 GeV");
  lt4ea->SetTextColor(kBlack); 
  lt4ea->SetTextSize(TextSize);
  lt4ea->SetTextFont(Font);
  lt4ea->Draw();
     
  
  TLegend* leg4ea = new TLegend(0.15,0.45,0.35,0.54);
  leg4ea->SetFillColor(0);
  leg4ea->SetLineColor(0);
  leg4ea->SetTextFont(Font);
  leg4ea->SetTextSize(0.041);
  leg4ea->AddEntry( eta2pi0_pAu,"p-Au","p");  
  leg4ea->AddEntry( eta2pi0_pBe,"p-Be","p");  
  leg4ea->Draw("same");
  
  
  TLatex* lt4eb = new TLatex(2.7,1.12,"PHENIX");
  lt4eb->SetTextColor(kBlack);
  lt4eb->SetTextSize(TextSize);
  lt4eb->SetTextFont(Font);
  lt4eb->Draw();
     
  TLegend* leg4ec = new TLegend(0.6,0.85,0.9,0.9);
  leg4ec->SetFillColor(0);
  leg4ec->SetLineColor(0);
  leg4ec->SetTextFont(Font);
  leg4ec->SetTextSize(0.041);
  leg4ec->AddEntry( graphEta2Pi0dAuPhenix,"d-Au, #sqrt{#it{s}_{NN}} = 200 GeV","p"); 
  leg4ec->Draw("same");
      
  
  c4e->Update();
  c4e->Print(Form("%s/EtaPi0Ratio_WO_Models.%s",outputDir.Data(),suffix.Data()));
  
  
  
  TCanvas* c4f = new TCanvas("c4f","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c4f, CMarginEtaToPi0L, CMarginEtaToPi0R,CMarginEtaToPi0T ,CMarginEtaToPi0B);
  c4f->SetLogx(); 
  
  TH2F* hist4f = new TH2F("hist4f","hist4f",1000,0.101,25.,1000,-0.05,1.2   );
  SetStyleHistoTH2ForGraphs(hist4f, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0} ",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeYEtaToPi0, TitleOffsetX,TitleOffsetYEtaToPi0, 512, 508);
  hist4f->GetYaxis()->SetLabelOffset(0.01);
  hist4f->GetXaxis()->SetLabelOffset(LabelOffsetLog); 
  hist4f->DrawCopy(); 
  
   
   graphMcGillEtaToPi0->Draw("C3same");
   graphMcGillEtaToPi0->Draw("lXZsame");

  
   graphEtaPi0EPOS->Draw("C3same");
   graphEtaPi0EPOS->Draw("lXZsame");
  
   mTScaling->Draw("same"); 
   
   
   
   if( HIJINGPi0 ){
   SetStyleHisto(histoDPMJetEtaToPi0, 3, styleLineDPMJet, colorDPMJet );  
   histoDPMJetEtaToPi0->Draw("same,hist,l");
   }

   if( DPMJetPi0 ){
   SetStyleHisto(histoHIJINGEtaToPi0,3, 8, colorHIJING);  
   histoHIJINGEtaToPi0->Draw("same,hist,l");
   }
   

  EtaPi0pPbStat_noXerrors->Draw("pz,same");
  EtaPi0pPbSyst->Draw("E2,same");
  

  
  TLatex * lt4f = new TLatex(.12,1.12,"ALICE");
  lt4f->SetTextColor(kBlack);
  lt4f->SetTextSize(TextSize);
  lt4f->SetTextFont(Font);
  lt4f->Draw() ;
  
  TLegend* leg4fa = new TLegend(0.15,0.56,0.45,0.90);
  leg4fa->SetFillColor(0);
  leg4fa->SetLineColor(0);
  leg4fa->SetTextFont(Font);
  leg4fa->SetTextSize(0.041);
  leg4fa->AddEntry(EtaPi0pPbSyst,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  leg4fa->AddEntry(mTScaling,"#eta from #it{m}_{T} scaled #pi^{0}","pl"); 
  leg4fa->AddEntry(graphEtaPi0EPOS,"EPOS3","fl");
  leg4fa->AddEntry(graphMcGillEtaToPi0,"VISHNU","fl");
  leg4fa->AddEntry(histoHIJINGEtaToPi0,"HIJING","l"); 
  if( DPMJetPi0 ){
     leg4fa->AddEntry(histoDPMJetEtaToPi0,"DPMJet","l"); 
  }
  leg4fa->Draw("same");
  
  


  c4f->Update();
  c4f->Print(Form("%s/EtaPi0Ratio_LogX_EPOS_TAPSv2.%s",outputDir.Data(),suffix.Data()));
  
  
  


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
  TBox* BoxNorm =new TBox(-0.5, 1.-NormalizationError, -0.15, 1.+NormalizationError);
  BoxNorm->SetFillColor(kBlack);
  BoxNorm->Draw("same");
  
  TBox* BoxNormEta = (TBox*) BoxNorm->Clone();
  TBox* BoxNormPi0 = (TBox*) BoxNorm->Clone();
  
  BoxNormEta->SetFillColor(colorCombYieldEta);
  BoxNormPi0->SetFillColor(colorCombYieldPi0);
  
  
  

  c5->Update();
  c5->Print(Form("%s/CombRpA_Models.%s",outputDir.Data(),suffix.Data()));

  //RpPb Version2
  TCanvas* c5a = new TCanvas("c5a","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c5a, CMarginL, CMarginR,CMarginT ,CMarginB);

  TH2F * hist5a;
  hist5a = new TH2F("hist5a","hist5a",1000,0.,21.,1000,0.55,1.65   );
  SetStyleHistoTH2ForGraphs(hist5a, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#pi^{0}}_{p-Pb}",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY,TitleOffsetX,TitleOffsetY, 512, 508);
  
  hist5a->GetYaxis()->SetLabelOffset(0.01);
  hist5a->GetXaxis()->SetLabelOffset(LabelOffsetLog); 
  
  
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
  //c5aa->SetLogx();
  

  TH2F * hist5aa;
  hist5aa = new TH2F("hist5aa","hist5aa",1000,-1.5,22.,1000,0.37,1.5   );
  SetStyleHistoTH2ForGraphs(hist5aa, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#pi^{0}}_{p-Pb}",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY,TitleOffsetX,TitleOffsetY, 512, 508);
  hist5aa->DrawCopy(); 
  
  TLegend* leg5abcde = new TLegend(0.16,0.79,0.40,0.86);
  leg5abcde->SetFillColor(0);
  leg5abcde->SetLineColor(0);
  leg5abcde->SetTextFont(Font);
  leg5abcde->SetTextSize(TextSizeRpA);
  leg5abcde->AddEntry(CombinedPi0RpPbSystErr,"#pi^{0}","pef");
  leg5abcde->Draw("same");
 
 
  graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14->SetLineWidth(lineWidthIlkka+1);   
  graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14->SetFillColor(colorIlkka);
  graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14->SetLineColor(colorIlkkaline);
  graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14->SetLineStyle(styleLineIlkka);
  graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14->SetFillStyle(1001);
  graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14->Draw("same,E3");
  graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14->Draw("lXYsame");
  
  //CGC->SetMarkerSize(0.65);
  //CGC->Draw("same,p");
  
  //-AM
  cout<< " Drawing  graphNLOCalcDSS14RpAPi05023GeV_nCTEQ_SepCalc"<< endl;
  
  
         
   DrawGammaSetMarkerTGraphAsym(graphNLOCalcDSS14RpAPi05023GeV_nCTEQ_SepCalc, 0, 0, colorDSSnPDFEPPSBand, colorDSSnPDFEPPSBand, widthLinesBoxes, kTRUE, colorDSSnPDFEPPSBand, kTRUE);
   graphNLOCalcDSS14RpAPi05023GeV_nCTEQ_SepCalc->SetLineStyle(styleLineDSSnPDFEPPS);
   graphNLOCalcDSS14RpAPi05023GeV_nCTEQ_SepCalc->Draw("lZXsame");
   graphNLOCalcDSS14RpAPi05023GeV_nCTEQ_SepCalc->Draw("3,same");
  
 
  

  //-AM For the R_pA of eta meson
  cout<< " Drawing graphNLOCalcAESSSRpAEta5023GeV_nCTEQ_SepCalc"<< endl;
  
  TGraphAsymmErrors* CombinedPi0RpPbStatErr_noXerrors = (TGraphAsymmErrors*)CombinedPi0RpPbStatErr->Clone();
  SetEx(CombinedPi0RpPbStatErr_noXerrors,0);
  
  //DrawGammaSetMarkerTGraphAsym(CombinedPi0RpPbStatErr_noXerrors, 20, 1, kBlack, kBlack);
  //DrawGammaSetMarkerTGraphAsym(CombinedPi0RpPbSystErr, 20, 1, kBlack, kBlack, 1., kTRUE);
  
  DrawGammaSetMarkerTGraphAsym(CombinedPi0RpPbSystErr,markerStylePi0,1.2,colorCombYieldPi0 ,colorCombYieldPi0 , 1, kTRUE);  
  DrawGammaSetMarkerTGraphAsym(CombinedPi0RpPbStatErr_noXerrors,markerStylePi0,1.2,colorCombYieldPi0, colorCombYieldPi0 );  
   
  
  
  CombinedPi0RpPbSystErr->Draw("E2,same");
  CombinedPi0RpPbStatErr_noXerrors->Draw("pz,same"); 
  
  graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14->Draw("same,E3");
  graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14->Draw("lXYsame");
  
  //CGC->SetMarkerSize(0.65);
  CGC->Draw("lXYsame");
  
 
  cout <<"=========RpPb binning================="<< endl;  
  CombinedPi0RpPbStatErr->Print();
  cout <<"=========================="<< endl;
  TLegend* leg5abc = new TLegend(0.31,0.20,0.54,0.34);
  leg5abc->SetFillColor(0);
  leg5abc->SetLineColor(0);
  leg5abc->SetTextFont(Font);
  leg5abc->SetTextSize(TextSizeRpA);
  leg5abc->AddEntry(graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14,"NLO: EPPS16, DSS14","lf");
  leg5abc->AddEntry(graphNLOCalcDSS14RpAPi05023GeV_nCTEQ_SepCalc,"NLO: nCTEQ, DSS14","lf");
  leg5abc->AddEntry(CGC,"CGC","p");
  leg5abc->Draw("same");

  TLegend* leg5abcd = new TLegend(0.51,0.38,0.64,0.41);
  leg5abcd->SetFillColor(0);
  leg5abcd->SetLineColor(0);
  leg5abcd->SetTextFont(Font);
  leg5abcd->SetTextSize(TextSizeRpA);

  
  
  
  TLatex * lt4ab = new TLatex(0.2,1.4,"ALICE");
  lt4ab->SetTextColor(kBlack) ;
  lt4ab->SetTextSize(TextSize) ;
  lt4ab->SetTextFont(Font) ;
  lt4ab->DrawLatex(0.2,1.35,"p-Pb, NSD, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  lt4ab->Draw();

  

  BoxNormPi0->Draw("same");
  DrawGammaLines(-1.5,22.,1.,1.,2.0,kGray+2,2);
 
  c5aa->Update();
  c5aa->Print(Form("%s/Comb_Pi0_RpA_Models_V3.%s",outputDir.Data(),suffix.Data()));
  
  
  
  //////////////////////Combined Eta RpPb/////////////////////////////////////
  
  TCanvas* c5ac = new TCanvas("c5ac","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c5ac, CMarginL, CMarginR,CMarginT ,CMarginB);
  

  TH2F * hist5ac;
  hist5ac = new TH2F("hist5ac","hist5ac",1000,-1.5,22.,1000,0.37,1.5   );
  SetStyleHistoTH2ForGraphs(hist5ac, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#eta}_{p-Pb}",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY,TitleOffsetX,TitleOffsetY, 512, 508);
  hist5ac->DrawCopy(); 
  hist5ac->GetYaxis()->SetLabelOffset(0.01);
  hist5ac->GetXaxis()->SetLabelOffset(LabelOffsetLog);
  
  
 
  
  TLegend* leg5abcdf = new TLegend(0.16,0.79,0.40,0.86);
  leg5abcdf->SetFillColor(0);
  leg5abcdf->SetLineColor(0);
  leg5abcdf->SetTextFont(Font);
  leg5abcdf->SetTextSize(TextSizeRpA);
  leg5abcdf->AddEntry(CombinedEtaRpPbSystErr,"#eta","pef");
  leg5abcdf->Draw("same");
 
 
  
         
   DrawGammaSetMarkerTGraphAsym(graphNLOCalcAESSSRpAEta5023GeV_nCTEQ_SepCalc, 0, 0, colorDSSnPDFEPPSBand, colorDSSnPDFEPPSBand, widthLinesBoxes, kTRUE, colorDSSnPDFEPPSBand, kTRUE);
   graphNLOCalcAESSSRpAEta5023GeV_nCTEQ_SepCalc->SetLineStyle(styleLineDSSnPDFEPPS);
   graphNLOCalcAESSSRpAEta5023GeV_nCTEQ_SepCalc->Draw("lZXsame");
   graphNLOCalcAESSSRpAEta5023GeV_nCTEQ_SepCalc->Draw("3,same");
  
 
  

  //-AM For the R_pA of eta meson
  cout<< " Drawing graphNLOCalcAESSSRpAEta5023GeV_nCTEQ_SepCalc"<< endl;
  
  TGraphAsymmErrors* CombinedEtaRpPbStatErr_noXerrors = (TGraphAsymmErrors*)CombinedEtaRpPbStatErr->Clone();
  SetEx(CombinedEtaRpPbStatErr_noXerrors,0);
  
  //DrawGammaSetMarkerTGraphAsym(CombinedEtaRpPbStatErr_noXerrors, 20, 1, kBlack, kBlack);
  //DrawGammaSetMarkerTGraphAsym(CombinedEtaRpPbSystErr, 20, 1, kBlack, kBlack, 1., kTRUE);
  
  
  DrawGammaSetMarkerTGraphAsym(CombinedEtaRpPbSystErr,markerStyleEta,1,colorCombYieldEta ,colorCombYieldEta, 1, kTRUE);
  DrawGammaSetMarkerTGraphAsym(CombinedEtaRpPbStatErr_noXerrors,markerStyleEta,1,colorCombYieldEta, colorCombYieldEta ); 

  
   
  
  CombinedEtaRpPbSystErr->Draw("E2,same");
  CombinedEtaRpPbStatErr_noXerrors->Draw("pz,same"); 
  
 
  TLegend* leg5abe = new TLegend(0.23,0.20,0.43,0.25);
  leg5abe->SetFillColor(0);
  leg5abe->SetLineColor(0);
  leg5abe->SetTextFont(Font);
  leg5abe->SetTextSize(TextSizeRpA);
  leg5abe->AddEntry(graphNLOCalcAESSSRpAEta5023GeV_nCTEQ_SepCalc,"NLO: nCTEQ, AESSS","lf");
  leg5abe->Draw("same");

  TLegend* leg5abcf = new TLegend(0.51,0.38,0.64,0.41);
  leg5abcf->SetFillColor(0);
  leg5abcf->SetLineColor(0);
  leg5abcf->SetTextFont(Font);
  leg5abcf->SetTextSize(TextSizeRpA);

  
  
  
  TLatex * lt4ad = new TLatex(0.2,1.4,"ALICE");
  lt4ad->SetTextColor(kBlack) ;
  lt4ad->SetTextSize(TextSize) ;
  lt4ad->SetTextFont(Font) ;
  lt4ad->DrawLatex(0.2,1.35,"p-Pb, NSD, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  lt4ad->Draw();

  

  BoxNormEta->Draw("same");
  DrawGammaLines(-1.5,22.,1.,1.,2.0,kGray+2,2);
 
  c5ac->Update();
  c5ac->Print(Form("%s/Comb_Eta_RpA_Models_V3.%s",outputDir.Data(),suffix.Data()));
  
  ///////////////////////////////////////////
  
  
  

  TCanvas* c5ad = new TCanvas("c5ad","",200,10,CDimX,CDimY);  // gives the pagesize
  DrawGammaCanvasSettings( c5ad, CMarginL, CMarginR,CMarginT ,CMarginB);
  c5ad->SetLogx();
  

  TH2F * hist5ad;
  hist5ad = new TH2F("hist5ad","hist5ad",1000,0.5,25.,1000,0.37,1.5);
  SetStyleHistoTH2ForGraphs(hist5ad, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#eta}_{p-Pb}",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY,TitleOffsetX,TitleOffsetY, 512, 508);
  hist5ad->GetYaxis()->SetLabelOffset(0.005);
  hist5ad->GetXaxis()->SetLabelOffset(LabelOffsetLog);
  hist5ad->DrawCopy(); 
  
  
  
  
  TLegend* leg5abcdeg = new TLegend(0.16,0.79,0.40,0.86);
  leg5abcdeg->SetFillColor(0);
  leg5abcdeg->SetLineColor(0);
  leg5abcdeg->SetTextFont(Font);
  leg5abcdeg->SetTextSize(TextSizeRpA);
  leg5abcdeg->AddEntry(CombinedEtaRpPbSystErr,"#eta","pef");
  leg5abcdeg->Draw("same");
  
  graphNLOCalcAESSSRpAEta5023GeV_nCTEQ_SepCalc->Draw("lZXsame");
  graphNLOCalcAESSSRpAEta5023GeV_nCTEQ_SepCalc->Draw("3,same");
 
  
  CombinedEtaRpPbSystErr->Draw("E2,same");
  CombinedEtaRpPbStatErr_noXerrors->Draw("pz,same"); 
  
 
  TLegend* leg5abg = new TLegend(0.41,0.20,0.64,0.27);
  leg5abg->SetFillColor(0);
  leg5abg->SetLineColor(0);
  leg5abg->SetTextFont(Font);
  leg5abg->SetTextSize(TextSizeRpA);
  leg5abg->AddEntry(graphNLOCalcAESSSRpAEta5023GeV_nCTEQ_SepCalc,"NLO: nCTEQ, AESSS","lf");
  leg5abg->Draw("same");

  
  TLegend* leg5abch = new TLegend(0.48,0.38,0.61,0.41);
  leg5abch->SetFillColor(0);
  leg5abch->SetLineColor(0);
  leg5abch->SetTextFont(Font);
  leg5abch->SetTextSize(TextSizeRpA);

  
  
  
  TLatex * lt4ae = new TLatex(0.65,1.4,"ALICE");
  lt4ae->SetTextColor(kBlack) ;
  lt4ae->SetTextSize(TextSize) ;
  lt4ae->SetTextFont(Font) ;
  lt4ae->DrawLatex(0.65,1.35,"p-Pb, NSD, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  lt4ae->Draw();

  TBox* BoxNormb =new TBox(0.55, 1.-NormalizationError, 0.60, 1.+NormalizationError);
  BoxNormb->SetFillColor(colorCombYieldEta);
  BoxNormb->Draw("same");
  

  
  DrawGammaLines(0.5,25.,1.,1.,2.0,kGray+2,2);
 
  c5ad->Update();
  c5ad->Print(Form("%s/Comb_Eta_RpA_Models_LogX_V3.%s",outputDir.Data(),suffix.Data()));
  
  
  

  TCanvas* c5ab = new TCanvas("c5ab","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c5ab, CMarginL, CMarginR,CMarginT ,CMarginB);
  c5ab->SetLogx();
  

  TH2F * hist5ab;
  hist5ab = new TH2F("hist5ab","hist5ab",1000,0.2,25.,1000,0.37,1.5);
  SetStyleHistoTH2ForGraphs(hist5ab, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#pi^{0}}_{p-Pb}",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY,TitleOffsetX,TitleOffsetY, 512, 508);

  
  hist5ab->GetYaxis()->SetLabelOffset(0.005);
  hist5ab->GetXaxis()->SetLabelOffset(LabelOffsetLog);
  hist5ab->DrawCopy(); 
  
 
  graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14->Draw("same,E3");
  graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14->Draw("lXYsame");
  
 
  
  graphNLOCalcDSS14RpAPi05023GeV_nCTEQ_SepCalc->Draw("lZXsame");
  graphNLOCalcDSS14RpAPi05023GeV_nCTEQ_SepCalc->Draw("3,same");
  
 
  
  DrawGammaSetMarkerTGraphAsym(CombinedPi0RpPbSystErr,markerStylePi0,1.2,colorCombYieldPi0 ,colorCombYieldPi0 , 1, kTRUE);  
  DrawGammaSetMarkerTGraphAsym(CombinedPi0RpPbStatErr_noXerrors,markerStylePi0,1.2,colorCombYieldPi0, colorCombYieldPi0 );  
 
  
  //DrawGammaSetMarkerTGraphAsym(CombinedPi0RpPbStatErr_noXerrors, 20, 1, kBlack, kBlack);
  //DrawGammaSetMarkerTGraphAsym(CombinedPi0RpPbSystErr, 20, 1, kBlack, kBlack, 1., kTRUE);
   
  
  CombinedPi0RpPbSystErr->Draw("E2,same");
  CombinedPi0RpPbStatErr_noXerrors->Draw("pz,same"); 
  
  CGCline->RemovePoint(0);
  CGCline->SetFillColor(0);
  CGCline->SetMarkerColor(kMagenta-2); 
  CGCline->SetMarkerStyle(21);
  CGCline->SetMarkerSize(1.5);
  CGCline->SetLineColor(kOrange+1);
  CGCline->SetLineWidth(3);
  CGCline->SetLineStyle(1);
  CGCline->Draw("lZXsame");
  CGCline->Draw("3,same");
  
    cout<<"///////////////tetwerewrwerwere////////////////////////////////"<<endl;
    CGCline->Print();
    cout<<"///////////////////////////////////////////////"<<endl;

  
  TLegend* leg5abcef = new TLegend(0.16,0.79,0.40,0.86);
  leg5abcef->SetFillColor(0);
  leg5abcef->SetLineColor(0);
  leg5abcef->SetTextFont(Font);
  leg5abcef->SetTextSize(TextSizeRpA);
  leg5abcef->AddEntry(CombinedPi0RpPbSystErr,"#pi^{0}","pef");
  leg5abcef->Draw("same");
  
  
 
  TLegend* leg5abd = new TLegend(0.59,0.20,0.84,0.34);
  leg5abd->SetFillColor(0);
  leg5abd->SetLineColor(0);
  leg5abd->SetTextFont(Font);
  leg5abd->SetTextSize(TextSizeRpA);
  leg5abd->AddEntry(graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14,"NLO: EPPS16, DSS14","lf");
  leg5abd->AddEntry(graphNLOCalcDSS14RpAPi05023GeV_nCTEQ_SepCalc,"NLO: nCTEQ, DSS14","lf");
  leg5abd->AddEntry(CGCline,"CGC","l");
  leg5abd->Draw("same");

  
  TLegend* leg5abce = new TLegend(0.48,0.38,0.61,0.41);
  leg5abce->SetFillColor(0);
  leg5abce->SetLineColor(0);
  leg5abce->SetTextFont(Font);
  leg5abce->SetTextSize(TextSizeRpA);

  
  
  
  TLatex * lt4ac = new TLatex(0.27,1.4,"ALICE");
  lt4ac->SetTextColor(kBlack) ;
  lt4ac->SetTextSize(TextSize) ;
  lt4ac->SetTextFont(Font) ;
  lt4ac->DrawLatex(0.27,1.35,"p-Pb, NSD, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  lt4ac->Draw();

  TBox* BoxNorma =new TBox(0.22, 1.-NormalizationError, 0.245, 1.+NormalizationError);
  BoxNorma->SetFillColor(colorCombYieldPi0);
  BoxNorma->Draw("same");
  

  
  DrawGammaLines(0.2,25.,1.,1.,2.0,kGray+2,2);
 
  c5ab->Update();
  c5ab->Print(Form("%s/Comb_Pi0_RpA_Models_LogX_V3.%s",outputDir.Data(),suffix.Data()));
  
  //////////////////////////////////////////////////////////////////////////////
  

  TCanvas* c5ae = new TCanvas("c5ae","",200,10,CDimX,CDimY);  // gives the page size
  DrawGammaCanvasSettings( c5ae, CMarginL, CMarginR,CMarginT ,CMarginB);
  c5ae->SetLogx();
  

  TH2F* hist5ae = new TH2F("hist5ae","hist5ae",1000,0.2,25.,1000,0.37,1.5);
  SetStyleHistoTH2ForGraphs(hist5ae, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#pi^{0}}_{p-Pb}",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY,TitleOffsetX,TitleOffsetY, 512, 508);
  hist5ae->GetYaxis()->SetLabelOffset(0.005);
  hist5ae->GetXaxis()->SetLabelOffset(LabelOffsetLog);
  hist5ae->DrawCopy(); 
  
  
  CombinedPi0RpPbSystErr->Draw("E2,same");
  CombinedPi0RpPbStatErr_noXerrors->Draw("pz,same"); 
  
  
  TLegend* leg5ae = new TLegend(0.16,0.79,0.40,0.86);
  leg5ae->SetFillColor(0);
  leg5ae->SetLineColor(0);
  leg5ae->SetTextFont(Font);
  leg5ae->SetTextSize(TextSizeRpA);
  leg5ae->AddEntry(CombinedPi0RpPbSystErr,"#pi^{0}","pef");
  leg5ae->Draw("same");
  
 
  
  TLatex * lt5ae = new TLatex(0.27,1.4,"ALICE");
  lt5ae->SetTextColor(kBlack) ;
  lt5ae->SetTextSize(TextSize) ;
  lt5ae->SetTextFont(Font) ;
  lt5ae->DrawLatex(0.27,1.35,"p-Pb, NSD, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  lt5ae->Draw();

  BoxNorma->Draw("same");
  

  
  DrawGammaLines(0.2,25.,1.,1.,2.0,kGray+2,2);
 
  c5ae->Update();
  c5ae->Print(Form("%s/Comb_Pi0_RpA_WO_Models_LogX_V3.%s",outputDir.Data(),suffix.Data()));
  
  /////////////////////////////////////////////////////////////////////////////////////////
  
  
  TCanvas* c5af = new TCanvas("c5af","",200,10,CDimX,CDimY);  // gives the pagesize
  DrawGammaCanvasSettings( c5af, CMarginL, CMarginR,CMarginT ,CMarginB);
  c5af->SetLogx();
  

  TH2F * hist5af = new TH2F("hist5af","hist5af",1000,0.5,25.,1000,0.37,1.5);
  SetStyleHistoTH2ForGraphs(hist5af, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#eta}_{p-Pb}",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY,TitleOffsetX,TitleOffsetY, 512, 508);
  hist5af->GetYaxis()->SetLabelOffset(0.005);
  hist5af->GetXaxis()->SetLabelOffset(LabelOffsetLog);
  hist5af->DrawCopy(); 
  
  
  
  CombinedEtaRpPbSystErr->Draw("E2,same");
  CombinedEtaRpPbStatErr_noXerrors->Draw("pz,same"); 
  
  TLegend* leg5af = new TLegend(0.16,0.79,0.40,0.86);
  leg5af->SetFillColor(0);
  leg5af->SetLineColor(0);
  leg5af->SetTextFont(Font);
  leg5af->SetTextSize(TextSizeRpA);
  leg5af->AddEntry(CombinedEtaRpPbSystErr,"#eta","pef");
  leg5af->Draw("same");
  
 
  
  
  
  TLatex * lt4af = new TLatex(0.65,1.4,"ALICE");
  lt4af->SetTextColor(kBlack) ;
  lt4af->SetTextSize(TextSize) ;
  lt4af->SetTextFont(Font) ;
  lt4af->DrawLatex(0.65,1.35,"p-Pb, NSD, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  lt4af->Draw();

  BoxNormb->Draw("same");
  

  
  DrawGammaLines(0.5,25.,1.,1.,2.0,kGray+2,2);
 
  c5af->Update();
  c5af->Print(Form("%s/Comb_Eta_RpA_WO_Models_LogX_V3.%s",outputDir.Data(),suffix.Data()));


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
  TGraphAsymmErrors* CombEtaScaledStat_noXerrors = (TGraphAsymmErrors*)ScaleGraph(CombEtaStat_noXerrors,scaleFacEtaYield);
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
  DrawGammaSetMarkerTGraphAsym(CombEtaScaledStat_noXerrors,markerStyleEta,1,colorCombYieldEta, colorCombYieldEta ); 
  //CombPi0ScaledStat->Draw("pz,same"); 	
  
  
  
  TCanvas* canvasInvYieldTSallisTheoryOnlyPi0AndEtaSpectra = new TCanvas("canvasInvYieldTSallisTheoryOnlyPi0AndEtaSpectra","",200,10,1000,1200);  // gives the page size
  DrawGammaCanvasSettings( canvasInvYieldTSallisTheoryOnlyPi0AndEtaSpectra,  0.3, 0.02, 0.02, 0.16);
	
  TPad* padComparisonInvYieldTSallisTheoryOnlyPi0AndEtaSpectra = new TPad("padComparisonInvYieldTSallisTheoryOnlyPi0AndEtaSpectra", "", 0., 0., 1., 1.,-1, -1, -2);
  DrawGammaPadSettings( padComparisonInvYieldTSallisTheoryOnlyPi0AndEtaSpectra, 0.15, 0.02, 0.02, 0.07);
  padComparisonInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->Draw();
	
  	
  padComparisonInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->cd();
  padComparisonInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->SetLogy();		
  padComparisonInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->SetLogx();		
	
  TH2F * histo2DInvYieldTSallisTheoryOnlyPi0AndEtaSpectra;
  histo2DInvYieldTSallisTheoryOnlyPi0AndEtaSpectra = new TH2F("histo2DInvYieldTSallisTheoryOnlyPi0AndEtaSpectra","histo2DInvYieldTSallisTheoryOnlyPi0AndEtaSpectra",1000,0.2,30.,1000,3e-11,20 );
  SetStyleHistoTH2ForGraphs(histo2DInvYieldTSallisTheoryOnlyPi0AndEtaSpectra, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2} ",0.035,0.035,0.035,0.035, 0.8,1.9, 512, 510);
  histo2DInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->GetXaxis()->SetLabelOffset(-0.009);
  histo2DInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->DrawCopy(); 
	
 
  
  FitCombPi0->SetLineWidth(lineWidthTsallis);  
  FitCombEtaScaled->SetLineWidth(lineWidthTsallis);  
  FitCombPi0->SetLineColor(kBlack);  
  FitCombEtaScaled->SetLineColor(kBlack);
  FitCombPi0->SetLineStyle(styleLineTsallis);  
  FitCombEtaScaled->SetLineStyle(styleLineTsallis);
  
   graphAsymmErrIlkkapPb5020_pi0_ct14_epps16_dss14_scale_err->SetLineWidth(lineWidthIlkka);
   graphAsymmErrIlkkapPb5020_pi0_ct14_epps16_dss14_scale_err->SetLineStyle(styleLineIlkka);
   graphAsymmErrIlkkapPb5020_pi0_ct14_epps16_dss14_scale_err->SetLineColor(colorIlkkaline);
   graphAsymmErrIlkkapPb5020_pi0_ct14_epps16_dss14_scale_err->SetFillColor(colorIlkka);
   
   TGraphAsymmErrors* graphAsymmErrIlkkapPb5020_pi0_ct14_epps16_dss14_scale_Noerr = (TGraphAsymmErrors*)graphAsymmErrIlkkapPb5020_pi0_ct14_epps16_dss14_scale_err->Clone();
   
   //raphAsymmErrRpPb5020_pi0_ct14_epps16_dss14->SetFillColor(kYellow-7);
  //graphAsymmErrRpPb5020_pi0_ct14_epps16_dss14->SetFillStyle(1001);
   graphAsymmErrIlkkapPb5020_pi0_ct14_epps16_dss14_scale_err->Draw("c3same");
   graphAsymmErrIlkkapPb5020_pi0_ct14_epps16_dss14_scale_Noerr->Draw("lZXsame");
   
   
  FitCombPi0->Draw("lsame"); 
  FitCombEtaScaled->Draw("lsame");
 
 TGraphAsymmErrors* graphPi0EPOSNoErr;
 TGraphAsymmErrors* graphEtaEPOSScaledNoErr;
 TGraphAsymmErrors* graphErrMcGillTheoryPion_p_hydroNoErr;
  
  if (EPOS){
      while(graphPi0EPOS->GetX()[0] < pi0PtMin)
          graphPi0EPOS->RemovePoint(0); 
    graphPi0EPOS->SetLineWidth(2);
    graphPi0EPOS->SetLineStyle(styleLineEPOS3);
    graphPi0EPOS->SetLineColor(colorEPOSline);
    graphPi0EPOS->SetFillColor(colorEPOS);
    
   graphPi0EPOSNoErr  = (TGraphAsymmErrors*)graphPi0EPOS->Clone();
    
    
    
    
   
     while(graphEtaEPOSScaled->GetX()[0] < 0.7)
          graphEtaEPOSScaled->RemovePoint(0); 
    graphEtaEPOSScaled->SetLineWidth(2);
    graphEtaEPOSScaled->SetLineStyle(styleLineEPOS3);
    graphEtaEPOSScaled->SetLineColor(colorEPOSline);
    graphEtaEPOSScaled->SetFillColor(colorEPOS);
    
    graphEtaEPOSScaledNoErr = (TGraphAsymmErrors*)graphEtaEPOSScaled->Clone();
    
    
    
    graphPi0EPOS->Draw("C3same"); //C
    graphPi0EPOSNoErr->Draw("lZXsame");
    
    graphEtaEPOSScaled->Draw("C3same");
    graphEtaEPOSScaledNoErr->Draw("lZXsame");
    
    
    while(graphErrMcGillTheoryPion_p_hydro->GetX()[0] < pi0PtMin)
          graphErrMcGillTheoryPion_p_hydro->RemovePoint(0);  
    
    graphErrMcGillTheoryPion_p_hydro->SetLineWidth(3);
    graphErrMcGillTheoryPion_p_hydro->SetLineStyle(styleLineMcGill);
    graphErrMcGillTheoryPion_p_hydro->SetLineColor(colorMcGillline);
    graphErrMcGillTheoryPion_p_hydro->SetFillColor(colorMcGill);
    
    graphErrMcGillTheoryPion_p_hydroNoErr = (TGraphAsymmErrors*)graphErrMcGillTheoryPion_p_hydro->Clone();
    
    graphErrMcGillTheoryPion_p_hydro->Draw("C3same");
    graphErrMcGillTheoryPion_p_hydroNoErr->Draw("lXYsame");
    
   while(graphErrMcGillTheoryEta_p_hydroScaled->GetX()[0] < 0.7)
          graphErrMcGillTheoryEta_p_hydroScaled->RemovePoint(0);   
    graphErrMcGillTheoryEta_p_hydroScaled->SetLineWidth(3);
    graphErrMcGillTheoryEta_p_hydroScaled->SetLineStyle(styleLineMcGill);
    graphErrMcGillTheoryEta_p_hydroScaled->SetLineColor(colorMcGillline);
    graphErrMcGillTheoryEta_p_hydroScaled->SetFillColor(colorMcGill);
   
   
    TGraphAsymmErrors* graphErrMcGillTheoryEta_p_hydroScaledNoErr = (TGraphAsymmErrors*)graphErrMcGillTheoryEta_p_hydroScaled->Clone();
   
    graphErrMcGillTheoryEta_p_hydroScaled->Draw("c3same");
    graphErrMcGillTheoryEta_p_hydroScaledNoErr->Draw("lZXsame");
   
   
   
   
   
        
    
  }
  
  
  if( CGCPi0 ) {
    while(graphAsymmErrCGCTheoryPi0y0pA5020->GetX()[0] < pi0PtMin)
          graphAsymmErrCGCTheoryPi0y0pA5020->RemovePoint(0);  
    graphAsymmErrCGCTheoryPi0y0pA5020->SetLineWidth(3);
    graphAsymmErrCGCTheoryPi0y0pA5020->SetLineStyle(styleLineCGC);
    graphAsymmErrCGCTheoryPi0y0pA5020->SetLineColor(kOrange+1);
    graphAsymmErrCGCTheoryPi0y0pA5020->SetFillColor(kOrange+1);
    graphAsymmErrCGCTheoryPi0y0pA5020->Draw("C3same");    
  }
  
  if( HIJINGPi0 ){
    SetStyleHisto(histoHIJINGPi0, 3, styleLineHIJING, colorHIJINGPi0);          
    histoHIJINGPi0->GetXaxis()->SetRangeUser(pi0PtMin,pi0PtMax);
    histoHIJINGPi0->Draw("same,hist,l");
    
    SetStyleHisto(histoHIJINGEtaScaled, 3, styleLineHIJING, colorHIJINGEta);    
     histoHIJINGEtaScaled->GetXaxis()->SetRangeUser(0.7,20);
    histoHIJINGEtaScaled->Draw("same,hist,l");
    
    
    
  }
  
  if( DPMJetPi0 ){
    
   SetStyleHisto(histoDPMJetPi0, 3, styleLineDPMJet, colorDPMJetPi0);      
   histoDPMJetPi0->GetXaxis()->SetRangeUser(pi0PtMin,pi0PtMax);
   histoDPMJetPi0->Draw("same,hist,l"); 
   
   SetStyleHisto(histoDPMJetEtaScaled, 3, styleLineDPMJet, colorDPMJetEta);     
   histoDPMJetEtaScaled->GetXaxis()->SetRangeUser(0.7,20);
   histoDPMJetEtaScaled->Draw("same,hist,l"); 
   
  }
  
  
  //-AM pi0 yield Vogelsang
  // graphNLOCalcDSS14InvYieldPi05023GeV_nCTEQ->Draw("same");
  
  
  
  DrawGammaSetMarkerTGraphAsym(graphNLOCalcDSS14InvYieldPi05023GeV_nCTEQ, 0, 0, colorDSSnPDFEPPSBand, colorDSSnPDFEPPSBand, widthLinesBoxes, kTRUE, colorDSSnPDFEPPSBand, kTRUE);
  
  graphNLOCalcDSS14InvYieldPi05023GeV_nCTEQ->SetLineStyle(styleLineDSSnPDFEPPS);
  
  graphNLOCalcDSS14InvYieldPi05023GeV_nCTEQ->Draw("3,same");
  graphNLOCalcDSS14InvYieldPi05023GeV_nCTEQ->Draw("lZXsame");
  
  


  //-AM eta yield Vogelsang
  TGraphAsymmErrors*  graphNLOCalcAESSSInvYieldEta5023GeV_nCTEQScaled 	= (TGraphAsymmErrors*)ScaleGraph(graphNLOCalcAESSSInvYieldEta5023GeV_nCTEQ,scaleFacEtaYield);
 
  
   DrawGammaSetMarkerTGraphAsym(graphNLOCalcAESSSInvYieldEta5023GeV_nCTEQScaled, 0, 0, colorDSSnPDFEPPSBand, colorDSSnPDFEPPSBand, widthLinesBoxes, kTRUE, colorDSSnPDFEPPSBand, kTRUE);
   
   graphNLOCalcAESSSInvYieldEta5023GeV_nCTEQScaled->SetLineStyle(styleLineDSSnPDFEPPS);
   graphNLOCalcAESSSInvYieldEta5023GeV_nCTEQScaled->Draw("3,same");
   graphNLOCalcAESSSInvYieldEta5023GeV_nCTEQScaled->Draw("lZXsame");
  
  


  CombPi0Syst->Draw("E2,same"); 
  CombPi0Stat_noXerrors->Draw("pz,same"); 
  CombEtaScaledSyst->Draw("E2,same");
  CombEtaScaledStat_noXerrors->Draw("pz,same"); 

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
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra   = new TLegend(0.20,0.38,0.4,0.61); 
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->SetFillColor(0);
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->SetLineColor(0);
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->SetTextFont(Font);
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->SetTextSize(TextSize-0.015);
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->AddEntry(CombPi0Syst,"#pi^{0}","pef");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->AddEntry(graphPi0EPOS,"EPOS3","fl");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->AddEntry(graphErrMcGillTheoryPion_p_hydro,"VISHNU","fl");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->AddEntry(histoDPMJetPi0,"DPMJet","l");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->AddEntry(histoHIJINGPi0,"HIJING","l");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->AddEntry(graphAsymmErrCGCTheoryPi0y0pA5020,"CGC MV^{#gamma}","l");
 
  
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->AddEntry(graphAsymmErrIlkkapPb5020_pi0_ct14_epps16_dss14_scale_err,"NLO: EPPS16, DSS14","fl");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->AddEntry(graphNLOCalcDSS14InvYieldPi05023GeV_nCTEQ,"NLO: nCTEQ, DSS14","fl");
  
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->AddEntry(FitCombPi0,"Tsallis fit","l");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->Draw();
  
  TLegend* legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2;
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2   = new TLegend(0.20,0.12,0.4,0.33); 
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2->SetFillColor(0);
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2->SetLineColor(0);
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2->SetTextFont(Font);
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2->SetTextSize(TextSize-0.015); 
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2->AddEntry(CombEtaScaledSyst,"#eta #times 10^{-1}","pef");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2->AddEntry(graphEtaEPOSScaled,"EPOS3","fl");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2->AddEntry(graphErrMcGillTheoryEta_p_hydroScaled,"VISHNU","fl");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2->AddEntry(histoDPMJetEtaScaled,"DPMJet","l");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2->AddEntry(histoHIJINGEtaScaled,"HIJING","l");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2->AddEntry(graphNLOCalcAESSSInvYieldEta5023GeV_nCTEQScaled,"NLO: nCTEQ, AESSS","fl"); 
  
  
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2->AddEntry(FitCombPi0,"Tsallis fit","l");
  legendInvYieldTSallisTheoryOnlyPi0AndEtaSpectrav2->Draw();
  
  
   padComparisonInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->Update();

   canvasInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->Update();
	
   if(EPOS) canvasInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->Print(Form("%s/MesonYields_EPOSv2.%s",outputDir.Data(),suffix.Data()));
   else if(mT) canvasInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->Print(Form("%s/MesonYields_mTv2.%s",outputDir.Data(),suffix.Data()));
   else canvasInvYieldTSallisTheoryOnlyPi0AndEtaSpectra->Print(Form("%s/MesonYieldsv2.%s",outputDir.Data(),suffix.Data()));
  
   
   
  TCanvas* canvasInvYieldTSallisOnlyPi0AndEtaSpectra = new TCanvas("canvasInvYieldTSallisOnlyPi0AndEtaSpectra","",200,10,1000,1200);  // gives the page size
  DrawGammaCanvasSettings( canvasInvYieldTSallisOnlyPi0AndEtaSpectra,  0.3, 0.02, 0.02, 0.16);
	
  TPad* padComparisonInvYieldTSallisOnlyPi0AndEtaSpectra = new TPad("padComparisonInvYieldTSallisOnlyPi0AndEtaSpectra", "", 0., 0., 1., 1.,-1, -1, -2);
  DrawGammaPadSettings( padComparisonInvYieldTSallisOnlyPi0AndEtaSpectra, 0.15, 0.02, 0.02, 0.06);
  padComparisonInvYieldTSallisOnlyPi0AndEtaSpectra->Draw();
	
  	
  padComparisonInvYieldTSallisOnlyPi0AndEtaSpectra->cd();
  padComparisonInvYieldTSallisOnlyPi0AndEtaSpectra->SetLogy();		
  padComparisonInvYieldTSallisOnlyPi0AndEtaSpectra->SetLogx();		
	
  TH2F * histo2DInvYieldTSallisOnlyPi0AndEtaSpectra;
  histo2DInvYieldTSallisOnlyPi0AndEtaSpectra = new TH2F("histo2DInvYieldTSallisOnlyPi0AndEtaSpectra","histo2DInvYieldTSallisOnlyPi0AndEtaSpectra",1000,0.2,30.,1000,3e-10,20 );
  SetStyleHistoTH2ForGraphs(histo2DInvYieldTSallisOnlyPi0AndEtaSpectra, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2} ",0.035,0.035,0.035,0.035, 0.8,1.9, 512, 510);
  histo2DInvYieldTSallisOnlyPi0AndEtaSpectra->GetXaxis()->SetLabelOffset(-0.009);
  histo2DInvYieldTSallisOnlyPi0AndEtaSpectra->DrawCopy(); 
	
 
  
  FitCombPi0->SetLineWidth(lineWidthTsallis);  
  FitCombEtaScaled->SetLineWidth(lineWidthTsallis);  
  FitCombPi0->SetLineColor(kBlack);  
  FitCombEtaScaled->SetLineColor(kBlack);
  FitCombPi0->SetLineStyle(styleLineTsallis);  
  FitCombEtaScaled->SetLineStyle(styleLineTsallis);
  
   
  FitCombPi0->Draw("lsame"); 
  FitCombEtaScaled->Draw("lsame");
 
  
  


  CombPi0Syst->Draw("E2,same"); 
  CombPi0Stat_noXerrors->Draw("pz,same"); 
  CombEtaScaledSyst->Draw("E2,same");
  CombEtaScaledStat_noXerrors->Draw("pz,same"); 

  
  TLatex * latexInvYieldTSallisOnlyPi0AndEtaSpectra = new TLatex(2.5,4,"ALICE") ;
  latexInvYieldTSallisOnlyPi0AndEtaSpectra->SetTextColor(kBlack) ;
  latexInvYieldTSallisOnlyPi0AndEtaSpectra->SetTextSize(TextSize-0.013) ;
  latexInvYieldTSallisOnlyPi0AndEtaSpectra->SetTextFont(Font) ;
  latexInvYieldTSallisOnlyPi0AndEtaSpectra->DrawLatex(2.5,2,"p-Pb, NSD, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latexInvYieldTSallisOnlyPi0AndEtaSpectra->Draw() ;
  
  

  TLegend* legendInvYieldTSallisOnlyPi0AndEtaSpectra;
  legendInvYieldTSallisOnlyPi0AndEtaSpectra   = new TLegend(0.20,0.35,0.4,0.45); 
  legendInvYieldTSallisOnlyPi0AndEtaSpectra->SetFillColor(0);
  legendInvYieldTSallisOnlyPi0AndEtaSpectra->SetLineColor(0);
  legendInvYieldTSallisOnlyPi0AndEtaSpectra->SetTextFont(Font);
  legendInvYieldTSallisOnlyPi0AndEtaSpectra->SetTextSize(TextSize-0.015);
  legendInvYieldTSallisOnlyPi0AndEtaSpectra->AddEntry(CombPi0Syst,"#pi^{0}","pef");
  legendInvYieldTSallisOnlyPi0AndEtaSpectra->AddEntry(CombEtaScaledSyst,"#eta #times 10^{-1}","pef");
  legendInvYieldTSallisOnlyPi0AndEtaSpectra->AddEntry(FitCombPi0,"Tsallis fit","l");
  legendInvYieldTSallisOnlyPi0AndEtaSpectra->Draw();
  
  
  
  
  padComparisonInvYieldTSallisOnlyPi0AndEtaSpectra->Update();

  canvasInvYieldTSallisOnlyPi0AndEtaSpectra->Update();	
  canvasInvYieldTSallisOnlyPi0AndEtaSpectra->Print(Form("%s/MesonYields_WOModels.%s",outputDir.Data(),suffix.Data()));
   
   
  
  
  



  
  
  
  
  
  
  
  
  
  
  
  
  TCanvas* canvasInvYieldTSallisOnlyPi0AndEtaSpectraPPRef = new TCanvas("canvasInvYieldTSallisOnlyPi0AndEtaSpectraPPRef","",200,10,1000,1200);  // gives the page size
  DrawGammaCanvasSettings( canvasInvYieldTSallisOnlyPi0AndEtaSpectraPPRef,  0.3, 0.02, 0.02, 0.16);
	
  TPad* padComparisonInvYieldTSallisOnlyPi0AndEtaSpectraPPRef = new TPad("padComparisonInvYieldTSallisOnlyPi0AndEtaSpectraPPRef", "", 0., 0., 1., 1.,-1, -1, -2);
  DrawGammaPadSettings( padComparisonInvYieldTSallisOnlyPi0AndEtaSpectraPPRef, 0.15, 0.02, 0.02, 0.06);
  padComparisonInvYieldTSallisOnlyPi0AndEtaSpectraPPRef->Draw();
	
  	
  padComparisonInvYieldTSallisOnlyPi0AndEtaSpectraPPRef->cd();
  padComparisonInvYieldTSallisOnlyPi0AndEtaSpectraPPRef->SetLogy();		
  padComparisonInvYieldTSallisOnlyPi0AndEtaSpectraPPRef->SetLogx();		
	
  TH2F * histo2DInvYieldTSallisOnlyPi0AndEtaSpectraPPRef;
  histo2DInvYieldTSallisOnlyPi0AndEtaSpectraPPRef = new TH2F("histo2DInvYieldTSallisOnlyPi0AndEtaSpectraPPRef","histo2DInvYieldTSallisOnlyPi0AndEtaSpectraPPRef",1000,0.2,30.,1000,1e01,2e11);
  SetStyleHistoTH2ForGraphs(histo2DInvYieldTSallisOnlyPi0AndEtaSpectraPPRef, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2} ",0.035,0.035,0.035,0.035, 0.8,1.9, 512, 510);
  histo2DInvYieldTSallisOnlyPi0AndEtaSpectraPPRef->GetXaxis()->SetLabelOffset(-0.009);
  histo2DInvYieldTSallisOnlyPi0AndEtaSpectraPPRef->DrawCopy(); 
	
   TH1F* fitTsallisEtappRef5023GeVPtScaled = (TH1F*)fitTsallisEtappRef5023GeVPt->GetHistogram();
   fitTsallisEtappRef5023GeVPtScaled->Scale(scaleFacEtaYield);
  
   
 
  
  fitTsallisPi0ppRef5023GeVPt->SetLineWidth(lineWidthTsallis);  
  fitTsallisEtappRef5023GeVPtScaled->SetLineWidth(lineWidthTsallis);  
  fitTsallisPi0ppRef5023GeVPt->SetLineColor(kBlack);  
  fitTsallisEtappRef5023GeVPtScaled->SetLineColor(kBlack);
  fitTsallisPi0ppRef5023GeVPt->SetLineStyle(styleLineTsallis);  
  fitTsallisEtappRef5023GeVPtScaled->SetLineStyle(styleLineTsallis);
  
   
  fitTsallisPi0ppRef5023GeVPt->Draw("lsame"); 
  fitTsallisEtappRef5023GeVPtScaled->Draw("lsame");
 
 
  
  
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionSystErrPPRef,markerStylePi0,1,colorCombYieldPi0 ,colorCombYieldPi0, 1, kTRUE);
  graphCombPi0InvCrossSectionSystErrPPRef->Draw("E2,same"); 
  TGraphAsymmErrors* graphCombPi0InvCrossSectionStatErrPPRef_noXerrors = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionStatErrPPRef->Clone();
  SetEx(graphCombPi0InvCrossSectionStatErrPPRef_noXerrors,0.);
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionStatErrPPRef_noXerrors,markerStylePi0,1,colorCombYieldPi0, colorCombYieldPi0 ); 
  graphCombPi0InvCrossSectionStatErrPPRef_noXerrors->Draw("pz,same");
  
  
  TGraphAsymmErrors* graphCombEtaInvCrossSectionStatErrPPRef_noXerrors = (TGraphAsymmErrors*)graphCombEtaInvCrossSectionStatErrPPRef->Clone();
  SetEx(graphCombEtaInvCrossSectionStatErrPPRef_noXerrors,0.);
  
  TGraphAsymmErrors* graphCombEtaInvCrossSectionSystErrPPRefScaled = (TGraphAsymmErrors*)ScaleGraph(graphCombEtaInvCrossSectionSystErrPPRef,scaleFacEtaYield);
  TGraphAsymmErrors* graphCombEtaInvCrossSectionStatErrPPRefScaled_noXerrors = (TGraphAsymmErrors*)ScaleGraph(graphCombEtaInvCrossSectionStatErrPPRef_noXerrors,scaleFacEtaYield);
  
  
  
  
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionSystErrPPRefScaled,markerStyleEta,1,colorCombYieldEta ,colorCombYieldEta, 1, kTRUE);
  graphCombEtaInvCrossSectionSystErrPPRefScaled->Draw("E2,same");
  
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionStatErrPPRefScaled_noXerrors,markerStyleEta,1,colorCombYieldEta, colorCombYieldEta ); 
  graphCombEtaInvCrossSectionStatErrPPRefScaled_noXerrors->Draw("pz,same");



  
  TLatex * latexInvYieldTSallisOnlyPi0AndEtaSpectraPPRef = new TLatex(2.5,0.8e11,"ALICE") ;
  latexInvYieldTSallisOnlyPi0AndEtaSpectraPPRef->SetTextColor(kBlack) ;
  latexInvYieldTSallisOnlyPi0AndEtaSpectraPPRef->SetTextSize(TextSize-0.013) ;
  latexInvYieldTSallisOnlyPi0AndEtaSpectraPPRef->SetTextFont(Font) ;
  latexInvYieldTSallisOnlyPi0AndEtaSpectraPPRef->DrawLatex(2.5,0.4e11,"p-p reference, #sqrt{#it{s}} = 5.02 TeV");
  latexInvYieldTSallisOnlyPi0AndEtaSpectraPPRef->Draw() ;
  
  

  TLegend* legendInvYieldTSallisOnlyPi0AndEtaSpectraPPRef;
  legendInvYieldTSallisOnlyPi0AndEtaSpectraPPRef   = new TLegend(0.20,0.35,0.4,0.45); 
  legendInvYieldTSallisOnlyPi0AndEtaSpectraPPRef->SetFillColor(0);
  legendInvYieldTSallisOnlyPi0AndEtaSpectraPPRef->SetLineColor(0);
  legendInvYieldTSallisOnlyPi0AndEtaSpectraPPRef->SetTextFont(Font);
  legendInvYieldTSallisOnlyPi0AndEtaSpectraPPRef->SetTextSize(TextSize-0.015);
  legendInvYieldTSallisOnlyPi0AndEtaSpectraPPRef->AddEntry(CombPi0Syst,"#pi^{0}","pef");
  legendInvYieldTSallisOnlyPi0AndEtaSpectraPPRef->AddEntry(CombEtaScaledSyst,"#eta #times 10^{-1}","pef");
  legendInvYieldTSallisOnlyPi0AndEtaSpectraPPRef->AddEntry(FitCombPi0,"Tsallis fit","l");
  legendInvYieldTSallisOnlyPi0AndEtaSpectraPPRef->Draw();
  
  
  
  
  padComparisonInvYieldTSallisOnlyPi0AndEtaSpectraPPRef->Update();

  canvasInvYieldTSallisOnlyPi0AndEtaSpectraPPRef->Update();	
  canvasInvYieldTSallisOnlyPi0AndEtaSpectraPPRef->Print(Form("%s/MesonYields_PPRef_Eta_Pi0_Tsallis.%s",outputDir.Data(),suffix.Data()));
   
   
   
  
  TCanvas* canvasInvYieldTCMOnlyPi0AndEtaSpectraPPRef = new TCanvas("canvasInvYieldTCMOnlyPi0AndEtaSpectraPPRef","",200,10,1000,1200);  // gives the page size
  DrawGammaCanvasSettings( canvasInvYieldTCMOnlyPi0AndEtaSpectraPPRef,  0.3, 0.02, 0.02, 0.16);
	
  TPad* padComparisonInvYieldTCMOnlyPi0AndEtaSpectraPPRef = new TPad("padComparisonInvYieldTCMOnlyPi0AndEtaSpectraPPRef", "", 0., 0., 1., 1.,-1, -1, -2);
  DrawGammaPadSettings( padComparisonInvYieldTCMOnlyPi0AndEtaSpectraPPRef, 0.15, 0.02, 0.02, 0.06);
  padComparisonInvYieldTCMOnlyPi0AndEtaSpectraPPRef->Draw();
	
  	
  padComparisonInvYieldTCMOnlyPi0AndEtaSpectraPPRef->cd();
  padComparisonInvYieldTCMOnlyPi0AndEtaSpectraPPRef->SetLogy();		
  padComparisonInvYieldTCMOnlyPi0AndEtaSpectraPPRef->SetLogx();		
	
  TH2F * histo2DInvYieldTCMOnlyPi0AndEtaSpectraPPRef;
  histo2DInvYieldTCMOnlyPi0AndEtaSpectraPPRef = new TH2F("histo2DInvYieldTCMOnlyPi0AndEtaSpectraPPRef","histo2DInvYieldTCMOnlyPi0AndEtaSpectraPPRef",1000,0.2,30.,1000,1e01,2e11);
  SetStyleHistoTH2ForGraphs(histo2DInvYieldTCMOnlyPi0AndEtaSpectraPPRef, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2} ",0.035,0.035,0.035,0.035, 0.8,1.9, 512, 510);
  histo2DInvYieldTCMOnlyPi0AndEtaSpectraPPRef->GetXaxis()->SetLabelOffset(-0.009);
  histo2DInvYieldTCMOnlyPi0AndEtaSpectraPPRef->DrawCopy(); 
	
   TH1F* fitTCMEtappRef5023GeVPtScaled = (TH1F*)fitTCMEtappRef5023GeVPt->GetHistogram();
   fitTCMEtappRef5023GeVPtScaled->Scale(scaleFacEtaYield);
  
   
 
  
  fitTCMPi0ppRef5023GeVPt->SetLineWidth(lineWidthTsallis);  
  fitTCMEtappRef5023GeVPtScaled->SetLineWidth(lineWidthTsallis);  
  fitTCMPi0ppRef5023GeVPt->SetLineColor(kBlack);  
  fitTCMEtappRef5023GeVPtScaled->SetLineColor(kBlack);
  fitTCMPi0ppRef5023GeVPt->SetLineStyle(styleLineTsallis);  
  fitTCMEtappRef5023GeVPtScaled->SetLineStyle(styleLineTsallis);
  
   
  fitTCMPi0ppRef5023GeVPt->Draw("lsame"); 
  fitTCMEtappRef5023GeVPtScaled->Draw("lsame");
 
 
  
  
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionSystErrPPRef,markerStylePi0,1,colorCombYieldPi0 ,colorCombYieldPi0, 1, kTRUE);
  graphCombPi0InvCrossSectionSystErrPPRef->Draw("E2,same"); 
  //TGraphAsymmErrors* graphCombPi0InvCrossSectionStatErrPPRef_noXerrors = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionStatErrPPRef->Clone();
  //SetEx(graphCombPi0InvCrossSectionStatErrPPRef_noXerrors,0.);
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionStatErrPPRef_noXerrors,markerStylePi0,1,colorCombYieldPi0, colorCombYieldPi0 ); 
  graphCombPi0InvCrossSectionStatErrPPRef_noXerrors->Draw("pz,same");
  
  
  //TGraphAsymmErrors* graphCombEtaInvCrossSectionStatErrPPRef_noXerrors = (TGraphAsymmErrors*)graphCombEtaInvCrossSectionStatErrPPRef->Clone();
  //SetEx(graphCombEtaInvCrossSectionStatErrPPRef_noXerrors,0.);
  
  //TGraphAsymmErrors* graphCombEtaInvCrossSectionSystErrPPRefScaled = (TGraphAsymmErrors*)ScaleGraph(graphCombEtaInvCrossSectionSystErrPPRef,scaleFacEtaYield);
  //TGraphAsymmErrors* graphCombEtaInvCrossSectionStatErrPPRefScaled_noXerrors = (TGraphAsymmErrors*)ScaleGraph(graphCombEtaInvCrossSectionStatErrPPRef_noXerrors,scaleFacEtaYield);
  
  
  
  
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionSystErrPPRefScaled,markerStyleEta,1,colorCombYieldEta ,colorCombYieldEta, 1, kTRUE);
  graphCombEtaInvCrossSectionSystErrPPRefScaled->Draw("E2,same");
  
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionStatErrPPRefScaled_noXerrors,markerStyleEta,1,colorCombYieldEta, colorCombYieldEta ); 
  graphCombEtaInvCrossSectionStatErrPPRefScaled_noXerrors->Draw("pz,same");



  
  TLatex * latexInvYieldTCMOnlyPi0AndEtaSpectraPPRef = new TLatex(2.5,0.8e11,"ALICE") ;
  latexInvYieldTCMOnlyPi0AndEtaSpectraPPRef->SetTextColor(kBlack) ;
  latexInvYieldTCMOnlyPi0AndEtaSpectraPPRef->SetTextSize(TextSize-0.013) ;
  latexInvYieldTCMOnlyPi0AndEtaSpectraPPRef->SetTextFont(Font) ;
  latexInvYieldTCMOnlyPi0AndEtaSpectraPPRef->DrawLatex(2.5,0.4e11,"p-p reference, #sqrt{#it{s}} = 5.02 TeV");
  latexInvYieldTCMOnlyPi0AndEtaSpectraPPRef->Draw() ;
  
  

  TLegend* legendInvYieldTCMOnlyPi0AndEtaSpectraPPRef;
  legendInvYieldTCMOnlyPi0AndEtaSpectraPPRef   = new TLegend(0.20,0.35,0.4,0.45); 
  legendInvYieldTCMOnlyPi0AndEtaSpectraPPRef->SetFillColor(0);
  legendInvYieldTCMOnlyPi0AndEtaSpectraPPRef->SetLineColor(0);
  legendInvYieldTCMOnlyPi0AndEtaSpectraPPRef->SetTextFont(Font);
  legendInvYieldTCMOnlyPi0AndEtaSpectraPPRef->SetTextSize(TextSize-0.015);
  legendInvYieldTCMOnlyPi0AndEtaSpectraPPRef->AddEntry(CombPi0Syst,"#pi^{0}","pef");
  legendInvYieldTCMOnlyPi0AndEtaSpectraPPRef->AddEntry(CombEtaScaledSyst,"#eta #times 10^{-1}","pef");
  legendInvYieldTCMOnlyPi0AndEtaSpectraPPRef->AddEntry(FitCombPi0,"TCM fit","l");
  legendInvYieldTCMOnlyPi0AndEtaSpectraPPRef->Draw();
  
  
  
  
  padComparisonInvYieldTCMOnlyPi0AndEtaSpectraPPRef->Update();

  canvasInvYieldTCMOnlyPi0AndEtaSpectraPPRef->Update();	
  canvasInvYieldTCMOnlyPi0AndEtaSpectraPPRef->Print(Form("%s/MesonYields_PPRef_Eta_Pi0_TCM.%s",outputDir.Data(),suffix.Data()));
   
   
   
   
  
  
  TCanvas* canvasInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef = new TCanvas("canvasInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef","",200,10,1000,1200);  // gives the page size
  DrawGammaCanvasSettings( canvasInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef,  0.3, 0.02, 0.02, 0.16);
	
  TPad* padComparisonInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef = new TPad("padComparisonInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef", "", 0., 0., 1., 1.,-1, -1, -2);
  DrawGammaPadSettings( padComparisonInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef, 0.15, 0.02, 0.02, 0.06);
  padComparisonInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->Draw();
	
  	
  padComparisonInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->cd();
  padComparisonInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->SetLogy();		
  padComparisonInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->SetLogx();		
	
  TH2F * histo2DInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef;
  histo2DInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef = new TH2F("histo2DInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef","histo2DInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef",1000,0.2,30.,1000,3e-10,20 );
  SetStyleHistoTH2ForGraphs(histo2DInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2} ",0.035,0.035,0.035,0.035, 0.8,1.9, 512, 510);
  histo2DInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->GetXaxis()->SetLabelOffset(-0.009);
  histo2DInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->DrawCopy(); 
	
 
  
  FitCombPi0->SetLineWidth(lineWidthTsallis);  
  FitCombEtaScaled->SetLineWidth(lineWidthTsallis);  
  FitCombPi0->SetLineColor(kBlack);  
  FitCombEtaScaled->SetLineColor(kBlack);
  FitCombPi0->SetLineStyle(styleLineTsallis);  
  FitCombEtaScaled->SetLineStyle(styleLineTsallis);
  
  Double_t fTpPb              = 0.0983e3*(1/recalcBarn);
  //Double_t ScalingPPRef = 
  
  TH1D* fitTsallisPi0ppRef5023GeVPtTpPbScaled = (TH1D*)fitTsallisPi0ppRef5023GeVPt->GetHistogram();
  fitTsallisPi0ppRef5023GeVPtTpPbScaled->Scale(fTpPb);
  TH1D* fitTsallisEtappRef5023GeVPtTpPbScaled = (TH1D*)fitTsallisEtappRef5023GeVPtScaled->Clone();
  fitTsallisEtappRef5023GeVPtTpPbScaled->Scale(fTpPb);
  
  
  TGraphAsymmErrors* graphCombPi0InvCrossSectionSystErrPPRefTpPbScaled = (TGraphAsymmErrors*) ScaleGraph(graphCombPi0InvCrossSectionSystErrPPRef,fTpPb);
  TGraphAsymmErrors* graphCombPi0InvCrossSectionStatErrPPRef_noXerrorsTpPbScaled = (TGraphAsymmErrors*)ScaleGraph( graphCombPi0InvCrossSectionStatErrPPRef_noXerrors,fTpPb);
  
  TGraphAsymmErrors* graphCombEtaInvCrossSectionSystErrPPReTpPbScaled = (TGraphAsymmErrors*) ScaleGraph(graphCombEtaInvCrossSectionSystErrPPRefScaled,fTpPb);
  TGraphAsymmErrors* graphCombEtaInvCrossSectionStatErrPPRef_noXerrorsTpPbScaled = (TGraphAsymmErrors*)ScaleGraph( graphCombEtaInvCrossSectionStatErrPPRefScaled_noXerrors,fTpPb);
  
  
  
  fitTsallisPi0ppRef5023GeVPtTpPbScaled->SetLineWidth(lineWidthTsallis);  
  fitTsallisPi0ppRef5023GeVPtTpPbScaled->SetLineColor(kBlack);  
  fitTsallisPi0ppRef5023GeVPtTpPbScaled->SetLineStyle(styleLineTsallis);  
  
  fitTsallisEtappRef5023GeVPtTpPbScaled->SetLineWidth(lineWidthTsallis);  
  fitTsallisEtappRef5023GeVPtTpPbScaled->SetLineColor(kBlack);  
  fitTsallisEtappRef5023GeVPtTpPbScaled->SetLineStyle(styleLineTsallis);  
  
  
  

  CombPi0Syst->Draw("E2,same"); 
  CombPi0Stat_noXerrors->Draw("pz,same");
  
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionSystErrPPRefTpPbScaled,markerStylePi0+4,1,colorCombYieldPi0 ,colorCombYieldPi0, 1, kTRUE);
  graphCombPi0InvCrossSectionSystErrPPRefTpPbScaled->Draw("E2,same"); 
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionStatErrPPRef_noXerrorsTpPbScaled,markerStylePi0+4,1,colorCombYieldPi0, colorCombYieldPi0 ); 
  graphCombPi0InvCrossSectionStatErrPPRef_noXerrorsTpPbScaled->Draw("pz,same"); 
 
  
  
  FitCombPi0->Draw("lsame"); 
  fitTsallisPi0ppRef5023GeVPtTpPbScaled->Draw("lsame");
 
 
  
  CombEtaScaledSyst->Draw("E2,same");
  CombEtaScaledStat_noXerrors->Draw("pz,same"); 
  
  
    
  
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionSystErrPPReTpPbScaled,markerStyleEta+4,1,colorCombYieldEta ,colorCombYieldEta, 1, kTRUE);
  graphCombEtaInvCrossSectionSystErrPPReTpPbScaled->Draw("E2,same");
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionStatErrPPRef_noXerrorsTpPbScaled,markerStyleEta+4,1,colorCombYieldEta, colorCombYieldEta ); 
  graphCombEtaInvCrossSectionStatErrPPRef_noXerrorsTpPbScaled->Draw("pz,same"); 
 
  FitCombEtaScaled->Draw("lsame");
  fitTsallisEtappRef5023GeVPtTpPbScaled->Draw("lsame");
  
 
  
  TLatex * latexInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef = new TLatex(2.5,4,"ALICE") ;
  latexInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->SetTextColor(kBlack) ;
  latexInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->SetTextSize(TextSize-0.013) ;
  latexInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->SetTextFont(Font) ;
  latexInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->DrawLatex(2.5,2,"p-Pb, NSD, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latexInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->Draw() ;
  
  

  TLegend* legendInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef;
  legendInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef   = new TLegend(0.20,0.35,0.4,0.5); 
  legendInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->SetFillColor(0);
  legendInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->SetLineColor(0);
  legendInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->SetTextFont(Font);
  legendInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->SetTextSize(TextSize-0.015);
  legendInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->AddEntry(CombPi0Syst,"#pi^{0}","pef");
  legendInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->AddEntry(CombEtaScaledSyst,"#eta #times 10^{-1}","pef");
  legendInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->AddEntry(graphCombPi0InvCrossSectionSystErrPPRefTpPbScaled,"#pi^{0} pp reference #times <TpPb>","pef");
  legendInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->AddEntry(graphCombEtaInvCrossSectionSystErrPPReTpPbScaled,"#eta  pp reference #times <TpPb> #times 10^{-1}","pef");
  
  legendInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->AddEntry(FitCombPi0,"Tsallis fit","l");
  legendInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->Draw();
  
  
  
  
  padComparisonInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->Update();

  canvasInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->Update();	
  canvasInvYieldTSallisOnlyPi0AndEtaSpectraAndPPRef->Print(Form("%s/MesonYields_And_PPReferences_Tsallis.%s",outputDir.Data(),suffix.Data()));  
   
  
  
  
  
 
  
  
  TCanvas* canvasInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef = new TCanvas("canvasInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef","",200,10,1000,1200);  // gives the page size
  DrawGammaCanvasSettings( canvasInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef,  0.3, 0.02, 0.02, 0.16);
	
  TPad* padComparisonInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef = new TPad("padComparisonInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef", "", 0., 0., 1., 1.,-1, -1, -2);
  DrawGammaPadSettings( padComparisonInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef, 0.15, 0.02, 0.02, 0.07);
  padComparisonInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->Draw();
	
  	
  padComparisonInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->cd();
  padComparisonInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->SetLogy();		
  padComparisonInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->SetLogx();		
	
  TH2F * histo2DInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef;
  histo2DInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef = new TH2F("histo2DInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef","histo2DInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef",1000,0.2,30.,1000,3e-10,20 );
  SetStyleHistoTH2ForGraphs(histo2DInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2} ",0.035,0.035,0.035,0.035, 0.8,1.9, 512, 510);
  histo2DInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->GetXaxis()->SetLabelOffset(-0.009);
  histo2DInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->DrawCopy(); 
	
 
  
  FitCombPi0->SetLineWidth(lineWidthTsallis);  
  FitCombEtaScaled->SetLineWidth(lineWidthTsallis);  
  FitCombPi0->SetLineColor(kBlack);  
  FitCombEtaScaled->SetLineColor(kBlack);
  FitCombPi0->SetLineStyle(styleLineTsallis);  
  FitCombEtaScaled->SetLineStyle(styleLineTsallis);
  
   
  
  TH1D* fitTCMPi0ppRef5023GeVPtTpPbScaled = (TH1D*)fitTCMPi0ppRef5023GeVPt->GetHistogram();
  fitTCMPi0ppRef5023GeVPtTpPbScaled->Scale(fTpPb);
  TH1D* fitTCMEtappRef5023GeVPtTpPbScaled = (TH1D*)fitTCMEtappRef5023GeVPtScaled->Clone();
  fitTCMEtappRef5023GeVPtTpPbScaled->Scale(fTpPb);
  
  
 
  fitTCMPi0ppRef5023GeVPtTpPbScaled->SetLineWidth(lineWidthTCM);  
  fitTCMPi0ppRef5023GeVPtTpPbScaled->SetLineColor(kBlack);  
  fitTCMPi0ppRef5023GeVPtTpPbScaled->SetLineStyle(styleLineTCM);  
  
  fitTCMEtappRef5023GeVPtTpPbScaled->SetLineWidth(lineWidthTCM);  
  fitTCMEtappRef5023GeVPtTpPbScaled->SetLineColor(kBlack);  
  fitTCMEtappRef5023GeVPtTpPbScaled->SetLineStyle(styleLineTCM+1);  
  
  
  

  CombPi0Syst->Draw("E2,same"); 
  CombPi0Stat_noXerrors->Draw("pz,same");
  
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionSystErrPPRefTpPbScaled,markerStylePi0+4,1,colorPPRef ,colorPPRef, 1, kTRUE);
  //graphCombPi0InvCrossSectionSystErrPPRefTpPbScaled->Draw("E2,same"); 
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionStatErrPPRef_noXerrorsTpPbScaled,markerStylePi0+4,1,colorPPRef, colorPPRef ); 
  //graphCombPi0InvCrossSectionStatErrPPRef_noXerrorsTpPbScaled->Draw("pz,same"); 
 
  
  
  FitCombPi0->Draw("lsame"); 
  fitTCMPi0ppRef5023GeVPtTpPbScaled->Draw("lsame");
 
 
  
  CombEtaScaledSyst->Draw("E2,same");
  CombEtaScaledStat_noXerrors->Draw("pz,same"); 
  
  
    
  
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionSystErrPPReTpPbScaled,markerStyleEta+4,1,colorPPRef ,colorPPRef, 1, kTRUE);
  //graphCombEtaInvCrossSectionSystErrPPReTpPbScaled->Draw("E2,same");
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionStatErrPPRef_noXerrorsTpPbScaled,markerStyleEta+4,1,colorPPRef, colorPPRef ); 
  //graphCombEtaInvCrossSectionStatErrPPRef_noXerrorsTpPbScaled->Draw("pz,same"); 
 
  FitCombEtaScaled->Draw("lsame");
  fitTCMEtappRef5023GeVPtTpPbScaled->Draw("lsame");
  
 
  
  TLatex * latexInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef = new TLatex(2.5,4,"ALICE") ;
  latexInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->SetTextColor(kBlack) ;
  latexInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->SetTextSize(TextSize-0.013) ;
  latexInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->SetTextFont(Font) ;
  latexInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->DrawLatex(2.5,2,"p-Pb, NSD, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latexInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->Draw() ;
  
  

  TLegend* legendInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef;
  legendInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef   = new TLegend(0.20,0.30,0.4,0.48); 
  legendInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->SetFillColor(0);
  legendInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->SetLineColor(0);
  legendInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->SetTextFont(Font);
  legendInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->SetTextSize(TextSize-0.015);
  legendInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->AddEntry(CombPi0Syst,"#pi^{0}","pef");
  legendInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->AddEntry(CombEtaScaledSyst,"#eta #times 10^{-1}","pef");
  //legendInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->AddEntry(graphCombPi0InvCrossSectionSystErrPPRefTpPbScaled,"#pi^{0} pp reference #times <TpPb>","pef");
  //legendInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->AddEntry(graphCombEtaInvCrossSectionSystErrPPReTpPbScaled,"#eta  pp reference #times <TpPb> #times 10^{-1}","pef");
  
  legendInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->AddEntry(FitCombPi0,"Tsallis fit","l");
  legendInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->AddEntry(fitTCMEtappRef5023GeVPtTpPbScaled,"TCM fit #pi^{0} pp reference #times <TpPb>","l");
  legendInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->AddEntry(fitTCMPi0ppRef5023GeVPtTpPbScaled,"TCM fit #eta pp reference #times <TpPb> #times 10^{-1}","l");
  
  legendInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->Draw();
  
  
  
  
  padComparisonInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->Update();

  canvasInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->Update();	
  canvasInvYieldTCMOnlyPi0AndEtaSpectraAndPPRef->Print(Form("%s/MesonYields_And_PPReferences_TCM_NoPoints.%s",outputDir.Data(),suffix.Data())); 
  
   
   
   
   
   ////////////////////////////////////////////////////////////////////
   

   
   
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
    graphPi0EPOS->SetLineWidth(lineWidthEPOS3);
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
        //padComparisonInvYieldTSallisTheoryOnlyRatioPi0->SetLogy();
	
	
 	
	TH2F * histo2DInvYieldTSallisTheoryOnlyRatioPi0;
	histo2DInvYieldTSallisTheoryOnlyRatioPi0 = new TH2F("histo2DInvYieldTSallisTheoryOnlyRatioPi0","histo2DRatioAllppreferencesEtaandPi0",1000,.2,25.,1000,0.001,4.2);//1.89
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
        
	
	DrawGammaLines(0., 25.,1., 1.,2.0,kGray+2,2);
	
	RatioIlkkaFitPi0scaleerr->SetLineWidth(3);
        RatioIlkkaFitPi0scaleerr->SetLineStyle(styleLineIlkkaRatio);
        RatioIlkkaFitPi0scaleerr->SetFillColor(colorIlkkaPi0);
	RatioIlkkaFitPi0scaleerr->SetLineColor(colorIlkkaline);
	RatioIlkkaFitPi0scaleerr->Draw("C3same");
	
	
	RatioIlkkaFitPi0NoErr->SetLineWidth(3);
        RatioIlkkaFitPi0NoErr->SetLineStyle(styleLineIlkkaRatio);
	RatioIlkkaFitPi0NoErr->SetLineColor(colorIlkkaline);
	RatioIlkkaFitPi0NoErr->SetFillColor(colorIlkkaPi0);
	RatioIlkkaFitPi0NoErr->Draw("lEsame");
	



	if (EPOS){
	  while(RatioEPOSFitPi0->GetX()[0] < pi0PtMin)
          RatioEPOSFitPi0->RemovePoint(0); 
          
          RatioEPOSFitPi0->Draw("C3same");
          RatioEPOSFitPi0->Draw("lZXsame");
        }
  
	if( CGCPi0){
	  while(RatioCGCFitPi0->GetX()[0] < pi0PtMin)
	    RatioCGCFitPi0->RemovePoint(0); 
	    RatioCGCFitPi0->Draw("C3same");
        }
	
	
	RatioMcGillFitPi0->SetLineWidth(lineWidthMcGill);
        RatioMcGillFitPi0->SetLineStyle(styleLineMcGill);
	RatioMcGillFitPi0->SetLineColor(colorMcGillline);
	RatioMcGillFitPi0->SetFillColor(colorMcGill);
	 while(RatioMcGillFitPi0->GetX()[0] < pi0PtMin)
          RatioMcGillFitPi0->RemovePoint(0);
         
        RatioMcGillFitPi0->Draw("C3same");
        RatioMcGillFitPi0->Draw("lZXsame");
        
	cout<< " third addition AM : Plotting to Fit  for pi0 and eta "<< endl;
        
        //-AM
        
        TGraphAsymmErrors* RatioDSS14nCTEQFitPi05023GeVNoErr = (TGraphAsymmErrors*)RatioDSS14nCTEQFitPi05023GeV->Clone();
        
        RatioDSS14nCTEQFitPi05023GeVNoErr->SetLineWidth(lineWidthDSSnPDFEPPSLine);
        RatioDSS14nCTEQFitPi05023GeVNoErr->SetLineStyle(styleLineDSSnPDFEPPS);
	RatioDSS14nCTEQFitPi05023GeVNoErr->SetLineColor(colorDSSnPDFEPPSBand);
	RatioDSS14nCTEQFitPi05023GeVNoErr->SetFillColor(0);
	RatioDSS14nCTEQFitPi05023GeVNoErr->Draw("lZXsame");
        
        
       
        
        DrawGammaSetMarkerTGraphAsym(RatioDSS14nCTEQFitPi05023GeV, 0, 0, colorDSSnPDFEPPSBand, colorDSSnPDFEPPSBand, widthLinesBoxes, kTRUE, colorDSSnPDFEPPSBand, kTRUE);
        RatioDSS14nCTEQFitPi05023GeV->SetLineWidth(lineWidthDSSnPDFEPPSLine);
        RatioDSS14nCTEQFitPi05023GeV->SetLineStyle(styleLineDSSnPDFEPPS);
        RatioDSS14nCTEQFitPi05023GeV->Draw("3,same");

        
        
        

	
	if( HIJINGPi0 ){
           SetStyleHisto(histoRatioPi0HIJINGToFit,3, styleLineHIJING, colorHIJING );  
	   histoRatioPi0HIJINGToFit->GetXaxis()->SetRangeUser(pi0PtMin,pi0PtMax);
           histoRatioPi0HIJINGToFit->Draw("same,hist,l");  
	}

	if( DPMJetPi0 ){
	   SetStyleHisto(histoRatioPi0DPMJetToFit,3, styleLineDPMJet, colorDPMJet );  
	   histoRatioPi0HIJINGToFit->GetXaxis()->SetRangeUser(pi0PtMin,pi0PtMax);
	   histoRatioPi0DPMJetToFit->Draw("same,hist,l");  
	}
  
   
	RatioTsallisCombPi0Syst->Draw("E2,same");
	RatioTsallisCombPi0Stat->Draw("Ez,p,same"); 
	
	TLatex * latexInvYieldTSallisTheoryOnlyRatioPi0 = new TLatex(2.3,3.7,"ALICE") ;
	latexInvYieldTSallisTheoryOnlyRatioPi0->SetTextColor(kBlack) ;
	latexInvYieldTSallisTheoryOnlyRatioPi0->SetTextSize(TextSize+0.045) ;
	latexInvYieldTSallisTheoryOnlyRatioPi0->SetTextFont(Font) ;
	
	latexInvYieldTSallisTheoryOnlyRatioPi0->DrawLatex(2.3,3.3,"p-Pb, NSD, #sqrt{#it{s}_{NN}} = 5.02 TeV");
	latexInvYieldTSallisTheoryOnlyRatioPi0->Draw() ;
	
	
  
	TLegend* legendInvYieldTSallisTheoryOnlyRatioPi01;
	legendInvYieldTSallisTheoryOnlyRatioPi01   = new TLegend(0.18,0.43,0.38,0.93); 
	legendInvYieldTSallisTheoryOnlyRatioPi01->SetFillColor(0);
	legendInvYieldTSallisTheoryOnlyRatioPi01->SetFillStyle(0);
	legendInvYieldTSallisTheoryOnlyRatioPi01->SetLineColor(0);
	legendInvYieldTSallisTheoryOnlyRatioPi01->SetTextFont(Font);
	legendInvYieldTSallisTheoryOnlyRatioPi01->SetTextSize(TextSize+0.035);
	legendInvYieldTSallisTheoryOnlyRatioPi01->AddEntry(RatioTsallisCombPi0Syst,"#pi^{0}","pef");
       
 	legendInvYieldTSallisTheoryOnlyRatioPi01->AddEntry(RatioCGCFitPi0,"CGC MV^{#gamma}","l");
        legendInvYieldTSallisTheoryOnlyRatioPi01->AddEntry(RatioEPOSFitPi0,"EPOS3","fl");
 	legendInvYieldTSallisTheoryOnlyRatioPi01->AddEntry(RatioMcGillFitPi0,"VISHNU","fl");
        legendInvYieldTSallisTheoryOnlyRatioPi01->AddEntry(RatioIlkkaFitPi0NoErr,"NLO:","fl");
        legendInvYieldTSallisTheoryOnlyRatioPi01->AddEntry((TObject*)0,"EPPS16, DSS14","");
	legendInvYieldTSallisTheoryOnlyRatioPi01->Draw();
        

	padComparisonInvYieldTSallisTheoryOnlyRatioEta->cd();
	padComparisonInvYieldTSallisTheoryOnlyRatioEta->SetLogx();
        //padComparisonInvYieldTSallisTheoryOnlyRatioEta->SetLogy();
      
      
	TH2F * histo2DInvYieldTSallisTheoryOnlyRatioEta;
	histo2DInvYieldTSallisTheoryOnlyRatioEta = new TH2F("histo2DInvYieldTSallisTheoryOnlyRatioEta","histo2DRatioAllppreferencesEtaandEta",1000,.2,25.,1000,0.001,4.4);//1.89
	SetStyleHistoTH2ForGraphs(histo2DInvYieldTSallisTheoryOnlyRatioEta, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}",0.08,0.08,0.08,0.08,1.,0.5, 502, 505); 
	histo2DInvYieldTSallisTheoryOnlyRatioEta->GetYaxis()->SetLabelOffset(0.005);
	histo2DInvYieldTSallisTheoryOnlyRatioEta->GetXaxis()->SetLabelOffset(LabelOffsetLog);
	histo2DInvYieldTSallisTheoryOnlyRatioEta->GetXaxis()->SetTickLength(0.07);
	histo2DInvYieldTSallisTheoryOnlyRatioEta->GetYaxis()-> CenterTitle();
	histo2DInvYieldTSallisTheoryOnlyRatioEta->DrawCopy();

	
	DrawGammaLines(0., 25.,1., 1.,2.0,kGray+2,2);
	
	if (EPOS){
	 
	  while(RatioEPOSFitEta->GetX()[0] < 0.7)
          RatioEPOSFitEta->RemovePoint(0);   
	  
          TGraphAsymmErrors* RatioEPOSFitEtaNoErr = (TGraphAsymmErrors*)RatioEPOSFitEta->Clone();
           
	  RatioEPOSFitEta->Draw("C3same"); 
          RatioEPOSFitEtaNoErr->Draw("lZXsame");
	
      	 
	 
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
	
	
	RatioMcGillFitEta->SetLineWidth(lineWidthMcGill);
        RatioMcGillFitEta->SetLineStyle(styleLineMcGill);
	RatioMcGillFitEta->SetLineColor(colorMcGillline);
	RatioMcGillFitEta->SetFillColor(colorMcGill);
        while(RatioMcGillFitEta->GetX()[0] < 0.7)
        RatioMcGillFitEta->RemovePoint(0);   
	  
	RatioMcGillFitEta->Draw("C3same");
        
        TGraphAsymmErrors* RatioMcGillFitEtaNoErr = (TGraphAsymmErrors*) RatioMcGillFitEta->Clone();
        RatioMcGillFitEta->Draw("C3same");
        RatioMcGillFitEtaNoErr->Draw("lZXsame");
	
	//-AM	
        
        TGraphAsymmErrors* RatioAESSSnCTEQFitEta5023GeVNoErr = (TGraphAsymmErrors*)RatioAESSSnCTEQFitEta5023GeV->Clone();
        
        RatioAESSSnCTEQFitEta5023GeVNoErr->SetLineWidth(3);
        RatioAESSSnCTEQFitEta5023GeVNoErr->SetLineStyle(styleLineIlkkaRatio);
	RatioAESSSnCTEQFitEta5023GeVNoErr->SetLineColor(colorDSSnPDFEPPSBand);
	RatioAESSSnCTEQFitEta5023GeVNoErr->SetFillColor(colorDSSnPDFEPPSBand);
	RatioAESSSnCTEQFitEta5023GeVNoErr->Draw("lZXsame");
        
        DrawGammaSetMarkerTGraphAsym(RatioAESSSnCTEQFitEta5023GeV, 0, 0, colorDSSnPDFEPPSBand, colorDSSnPDFEPPSBand, widthLinesBoxes, kTRUE, colorDSSnPDFEPPSBand, kTRUE);
        RatioAESSSnCTEQFitEta5023GeV->Draw("3,same");

  
  
	RatioTsallisCombEtaSyst->Draw("E2same");
        RatioTsallisCombEtaStat->Draw("Ez,p,same");  
	
	
	TLegend* legendInvYieldTSallisTheoryOnlyRatioEta;
	legendInvYieldTSallisTheoryOnlyRatioEta   = new TLegend(0.18,0.60,0.38,0.95); 
	legendInvYieldTSallisTheoryOnlyRatioEta->SetFillColor(0);
	legendInvYieldTSallisTheoryOnlyRatioEta->SetLineColor(0);
	legendInvYieldTSallisTheoryOnlyRatioEta->SetTextFont(Font);
	legendInvYieldTSallisTheoryOnlyRatioEta->SetTextSize((TextSize+0.035)*0.82);
	legendInvYieldTSallisTheoryOnlyRatioEta->AddEntry(RatioTsallisCombEtaSyst,"#eta","pef");
        legendInvYieldTSallisTheoryOnlyRatioEta->AddEntry(RatioDSS14nCTEQFitPi05023GeV,"NLO:","fl");
        legendInvYieldTSallisTheoryOnlyRatioEta->AddEntry((TObject*)0,"nCTEQ, AESSS","");
 	legendInvYieldTSallisTheoryOnlyRatioEta->AddEntry(histoRatioPi0HIJINGToFit,"HIJING","l");
	legendInvYieldTSallisTheoryOnlyRatioEta->AddEntry(histoRatioPi0DPMJetToFit,"DPMJet","l");
	legendInvYieldTSallisTheoryOnlyRatioEta->Draw();
	

	canvasInvYieldTSallisTheoryOnlyRatioPi0AndEtaSpectra->Update();
	
	if(EPOS) canvasInvYieldTSallisTheoryOnlyRatioPi0AndEtaSpectra->Print(Form("%s/MesonYields_EPOS_OnlyRatiov2.%s",outputDir.Data(),suffix.Data()));
        else if(mT) canvasInvYieldTSallisTheoryOnlyRatioPi0AndEtaSpectra->Print(Form("%s/MesonYields_mT_OnlyRatiov2.%s",outputDir.Data(),suffix.Data()));
        else canvasInvYieldTSallisTheoryOnlyRatioPi0AndEtaSpectra->Print(Form("%s/MesonYields_OnlyRatiov2.%s",outputDir.Data(),suffix.Data()));
        
        
        
        
	TCanvas*    canvasInvYieldTSallisOnlyRatioPi0AndEtaSpectra  	= new TCanvas("canvasInvYieldTSallisOnlyRatioPi0AndEtaSpectra","",200,10,1000,800);  // gives the page size
	DrawGammaCanvasSettings( canvasInvYieldTSallisOnlyRatioPi0AndEtaSpectra,  0.15, 0.02, 0.00, 0.00);  

	TPad* padComparisonInvYieldTSallisOnlyRatioPi0 		= new TPad("padComparisonInvYieldTSallisOnlyRatioPi0", "", 0.,0.55, 1., 1,   -1, -1, -2);
	DrawGammaPadSettings( padComparisonInvYieldTSallisOnlyRatioPi0,  0.15, 0.02, 0.02, 0.0);
	padComparisonInvYieldTSallisOnlyRatioPi0->Draw();
	
 	TPad* padComparisonInvYieldTSallisOnlyRatioEta		= new TPad("padComparisonInvYieldTSallisOnlyRatioEta", "", 0.,0.0,  1., 0.55,-1, -1, -2);
 	DrawGammaPadSettings( padComparisonInvYieldTSallisOnlyRatioEta,  0.15, 0.02, 0.00,  0.20);
 	padComparisonInvYieldTSallisOnlyRatioEta->Draw();
 	
 	
 	padComparisonInvYieldTSallisOnlyRatioPi0->cd();
 	padComparisonInvYieldTSallisOnlyRatioPi0->SetLogx();
	
	
 	
	TH2F * histo2DInvYieldTSallisOnlyRatioPi0;
	histo2DInvYieldTSallisOnlyRatioPi0 = new TH2F("histo2DInvYieldTSallisOnlyRatioPi0","histo2DRatioAllppreferencesEtaandPi0",1000,.2,25.,1000,0.41,2.2);//1.89
	SetStyleHistoTH2ForGraphs(histo2DInvYieldTSallisOnlyRatioPi0, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{fit}", 0.098,0.098, 0.098,0.098, 1.,0.56, 502, 505); 
	histo2DInvYieldTSallisOnlyRatioPi0->GetYaxis()->SetLabelOffset(0.005);
	histo2DInvYieldTSallisOnlyRatioPi0->GetXaxis()->SetLabelOffset(LabelOffsetLog);
	histo2DInvYieldTSallisOnlyRatioPi0->GetXaxis()->SetTickLength(0.07);
	histo2DInvYieldTSallisOnlyRatioPi0->GetYaxis()-> CenterTitle();
	histo2DInvYieldTSallisOnlyRatioPi0->DrawCopy();
	
	
        
	
	DrawGammaLines(0., 25.,1., 1.,2.0,kGray+2,2);
	
	
	
	  
   
	RatioTsallisCombPi0Syst->Draw("E2,same");
	RatioTsallisCombPi0Stat->Draw("Ez,p,same"); 
	
	TLatex * latexInvYieldTSallisOnlyRatioPi0 = new TLatex(2.3,2.0,"ALICE") ;
	latexInvYieldTSallisOnlyRatioPi0->SetTextColor(kBlack) ;
	latexInvYieldTSallisOnlyRatioPi0->SetTextSize(TextSize+0.045) ;
	latexInvYieldTSallisOnlyRatioPi0->SetTextFont(Font) ;
	latexInvYieldTSallisOnlyRatioPi0->DrawLatex(2.3,1.8,"p-Pb, NSD, #sqrt{#it{s}_{NN}} = 5.02 TeV");
	latexInvYieldTSallisOnlyRatioPi0->Draw() ;
	
	
  
	TLegend* legendInvYieldTSallisOnlyRatioPi01;
	legendInvYieldTSallisOnlyRatioPi01   = new TLegend(0.18,0.80,0.38,0.93); 
	legendInvYieldTSallisOnlyRatioPi01->SetFillColor(0);
	legendInvYieldTSallisOnlyRatioPi01->SetFillStyle(0);
	legendInvYieldTSallisOnlyRatioPi01->SetLineColor(0);
	legendInvYieldTSallisOnlyRatioPi01->SetTextFont(Font);
	legendInvYieldTSallisOnlyRatioPi01->SetTextSize(TextSize+0.035);
	legendInvYieldTSallisOnlyRatioPi01->AddEntry(RatioTsallisCombPi0Syst,"#pi^{0}","pef");
        legendInvYieldTSallisOnlyRatioPi01->Draw();
        

	padComparisonInvYieldTSallisOnlyRatioEta->cd();
	padComparisonInvYieldTSallisOnlyRatioEta->SetLogx();
      
      
	TH2F * histo2DInvYieldTSallisOnlyRatioEta;
	histo2DInvYieldTSallisOnlyRatioEta = new TH2F("histo2DInvYieldTSallisOnlyRatioEta","histo2DRatioAllppreferencesEtaandEta",1000,.2,25.,1000,0.41,2.4);//1.89
	SetStyleHistoTH2ForGraphs(histo2DInvYieldTSallisOnlyRatioEta, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{fit}",0.08,0.08,0.08,0.08,1.,0.7, 502, 505); 
	histo2DInvYieldTSallisOnlyRatioEta->GetYaxis()->SetLabelOffset(0.005);
	histo2DInvYieldTSallisOnlyRatioEta->GetXaxis()->SetLabelOffset(LabelOffsetLog);
	histo2DInvYieldTSallisOnlyRatioEta->GetXaxis()->SetTickLength(0.07);
	histo2DInvYieldTSallisOnlyRatioEta->GetYaxis()-> CenterTitle();
	histo2DInvYieldTSallisOnlyRatioEta->DrawCopy();

        
        
	RatioTsallisCombEtaSyst->Draw("E2,same");
	RatioTsallisCombEtaStat->Draw("Ez,p,same"); 
	
	DrawGammaLines(0., 25.,1., 1.,2.0,kGray+2,2);
	
	
	
	
	TLegend* legendInvYieldTSallisOnlyRatioEta;
	legendInvYieldTSallisOnlyRatioEta   = new TLegend(0.18,0.83,0.38,0.95); 
	legendInvYieldTSallisOnlyRatioEta->SetFillColor(0);
	legendInvYieldTSallisOnlyRatioEta->SetLineColor(0);
	legendInvYieldTSallisOnlyRatioEta->SetTextFont(Font);
	legendInvYieldTSallisOnlyRatioEta->SetTextSize((TextSize+0.035)*0.82);
	legendInvYieldTSallisOnlyRatioEta->AddEntry(RatioTsallisCombEtaSyst,"#eta","pef");
	legendInvYieldTSallisOnlyRatioEta->Draw();
	

	canvasInvYieldTSallisOnlyRatioPi0AndEtaSpectra->Update();
	
        canvasInvYieldTSallisOnlyRatioPi0AndEtaSpectra->Print(Form("%s/MesonYields_WOModels_OnlyRatio.%s",outputDir.Data(),suffix.Data()));
       
        
        
        
        
        TCanvas*    canvasInvYieldTSallisOnlyRatioPi0AndEtaSpectraPPRef  	= new TCanvas("canvasInvYieldTSallisOnlyRatioPi0AndEtaSpectraPPRef","",200,10,1000,800);  // gives the page size
	DrawGammaCanvasSettings( canvasInvYieldTSallisOnlyRatioPi0AndEtaSpectraPPRef,  0.15, 0.02, 0.00, 0.00);  

	TPad* padComparisonInvYieldTSallisOnlyRatioPi0PPRef 		= new TPad("padComparisonInvYieldTSallisOnlyRatioPi0PPRef", "", 0.,0.55, 1., 1,   -1, -1, -2);
	DrawGammaPadSettings( padComparisonInvYieldTSallisOnlyRatioPi0PPRef,  0.15, 0.02, 0.02, 0.0);
	padComparisonInvYieldTSallisOnlyRatioPi0PPRef->Draw();
	
 	TPad* padComparisonInvYieldTSallisOnlyRatioEtaPPRef		= new TPad("padComparisonInvYieldTSallisOnlyRatioEtaPPRef", "", 0.,0.0,  1., 0.55,-1, -1, -2);
 	DrawGammaPadSettings( padComparisonInvYieldTSallisOnlyRatioEtaPPRef,  0.15, 0.02, 0.00,  0.20);
 	padComparisonInvYieldTSallisOnlyRatioEtaPPRef->Draw();
 	
 	
 	padComparisonInvYieldTSallisOnlyRatioPi0PPRef->cd();
 	padComparisonInvYieldTSallisOnlyRatioPi0PPRef->SetLogx();
	
	
 	
        
	TH2F * histo2DInvYieldTSallisOnlyRatioPi0PPRef;
	histo2DInvYieldTSallisOnlyRatioPi0PPRef = new TH2F("histo2DInvYieldTSallisOnlyRatioPi0PPRef","histo2DRatioAllppreferencesEtaandPi0",1000,.2,25.,1000,0.41,2.2);//1.89
	SetStyleHistoTH2ForGraphs(histo2DInvYieldTSallisOnlyRatioPi0PPRef, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{fit}", 0.098,0.098, 0.098,0.098, 1.,0.56, 502, 505); 
	histo2DInvYieldTSallisOnlyRatioPi0PPRef->GetYaxis()->SetLabelOffset(0.005);
	histo2DInvYieldTSallisOnlyRatioPi0PPRef->GetXaxis()->SetLabelOffset(LabelOffsetLog);
	histo2DInvYieldTSallisOnlyRatioPi0PPRef->GetXaxis()->SetTickLength(0.07);
	histo2DInvYieldTSallisOnlyRatioPi0PPRef->GetYaxis()-> CenterTitle();
	histo2DInvYieldTSallisOnlyRatioPi0PPRef->DrawCopy();
	
	
        
	
	DrawGammaLines(0., 25.,1., 1.,2.0,kGray+2,2);
	
	
	
	  
   
        
        
        
         DrawGammaSetMarkerTGraphAsym(graphRatioPi0ppRefTsallisFitSystErr,markerStylePi0,1.2,colorCombYieldPi0 ,colorCombYieldPi0 , 1, kTRUE);  
         graphRatioPi0ppRefTsallisFitSystErr->Draw("E2same");
         DrawGammaSetMarkerTGraphAsym(graphRatioPi0ppRefTsallisFitStatErr,markerStylePi0,1.2,colorCombYieldPi0, colorCombYieldPi0 );  
         SetEx(graphRatioPi0ppRefTsallisFitStatErr,0.); //Set the x-errors bars to 0
         graphRatioPi0ppRefTsallisFitStatErr->Draw("Ez,p,same");
   
  
        
        
	
	TLatex * latexInvYieldTSallisOnlyRatioPi0PPRef = new TLatex(2.3,2.0,"ALICE") ;
	latexInvYieldTSallisOnlyRatioPi0PPRef->SetTextColor(kBlack) ;
	latexInvYieldTSallisOnlyRatioPi0PPRef->SetTextSize(TextSize+0.045) ;
	latexInvYieldTSallisOnlyRatioPi0PPRef->SetTextFont(Font) ;
	latexInvYieldTSallisOnlyRatioPi0PPRef->DrawLatex(2.3,1.8,"p-p, reference, #sqrt{#it{s}} = 5.02 TeV");
	latexInvYieldTSallisOnlyRatioPi0PPRef->Draw() ;
	
	
  
	TLegend* legendInvYieldTSallisOnlyRatioPi01PPRef;
	legendInvYieldTSallisOnlyRatioPi01PPRef   = new TLegend(0.18,0.80,0.38,0.93); 
	legendInvYieldTSallisOnlyRatioPi01PPRef->SetFillColor(0);
	legendInvYieldTSallisOnlyRatioPi01PPRef->SetFillStyle(0);
	legendInvYieldTSallisOnlyRatioPi01PPRef->SetLineColor(0);
	legendInvYieldTSallisOnlyRatioPi01PPRef->SetTextFont(Font);
	legendInvYieldTSallisOnlyRatioPi01PPRef->SetTextSize(TextSize+0.035);
	legendInvYieldTSallisOnlyRatioPi01PPRef->AddEntry(RatioTsallisCombPi0Syst,"#pi^{0}","pef");
        legendInvYieldTSallisOnlyRatioPi01PPRef->Draw();
        

	padComparisonInvYieldTSallisOnlyRatioEtaPPRef->cd();
	padComparisonInvYieldTSallisOnlyRatioEtaPPRef->SetLogx();
      
      
	TH2F * histo2DInvYieldTSallisOnlyRatioEtaPPRef;
	histo2DInvYieldTSallisOnlyRatioEtaPPRef = new TH2F("histo2DInvYieldTSallisOnlyRatioEtaPPRef","histo2DRatioAllppreferencesEtaandEta",1000,.2,25.,1000,0.41,2.1);//1.89
	SetStyleHistoTH2ForGraphs(histo2DInvYieldTSallisOnlyRatioEtaPPRef, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{fit}",0.08,0.08,0.08,0.08,1.,0.7, 502, 505); 
	histo2DInvYieldTSallisOnlyRatioEtaPPRef->GetYaxis()->SetLabelOffset(0.005);
	histo2DInvYieldTSallisOnlyRatioEtaPPRef->GetXaxis()->SetLabelOffset(LabelOffsetLog);
	histo2DInvYieldTSallisOnlyRatioEtaPPRef->GetXaxis()->SetTickLength(0.07);
	histo2DInvYieldTSallisOnlyRatioEtaPPRef->GetYaxis()-> CenterTitle();
	histo2DInvYieldTSallisOnlyRatioEtaPPRef->DrawCopy();

        
        
	
        
        
       DrawGammaSetMarkerTGraphAsym(graphRatioEtappRefTsallisFitSystErr,markerStyleEta,1.2,colorCombYieldEta ,colorCombYieldEta, 1, kTRUE);  
       graphRatioEtappRefTsallisFitSystErr->Draw("E2same");
       DrawGammaSetMarkerTGraphAsym(graphRatioEtappRefTsallisFitStatErr,markerStyleEta,1.2,colorCombYieldEta ,colorCombYieldEta );  
       SetEx(graphRatioEtappRefTsallisFitStatErr,0.); //
       graphRatioEtappRefTsallisFitStatErr->Draw("Ez,p,same");
        
        
	
	DrawGammaLines(0., 25.,1., 1.,2.0,kGray+2,2);
	
	
	
	
	TLegend* legendInvYieldTSallisOnlyRatioEtaPPRef;
	legendInvYieldTSallisOnlyRatioEtaPPRef   = new TLegend(0.18,0.83,0.38,0.95); 
	legendInvYieldTSallisOnlyRatioEtaPPRef->SetFillColor(0);
	legendInvYieldTSallisOnlyRatioEtaPPRef->SetLineColor(0);
	legendInvYieldTSallisOnlyRatioEtaPPRef->SetTextFont(Font);
	legendInvYieldTSallisOnlyRatioEtaPPRef->SetTextSize((TextSize+0.035)*0.82);
	legendInvYieldTSallisOnlyRatioEtaPPRef->AddEntry(graphRatioEtappRefTsallisFitSystErr,"#eta","pef");
	legendInvYieldTSallisOnlyRatioEtaPPRef->Draw();
	

	canvasInvYieldTSallisOnlyRatioPi0AndEtaSpectraPPRef->Update();
	
        canvasInvYieldTSallisOnlyRatioPi0AndEtaSpectraPPRef->Print(Form("%s/MesonYields_PPRef_Eta_Pi0_Tsallis_OnlyRatios.%s",outputDir.Data(),suffix.Data()));
       
        
        
        
        
        
        
        TCanvas*    canvasInvYieldTCMOnlyRatioPi0AndEtaSpectraPPRef  	= new TCanvas("canvasInvYieldTCMOnlyRatioPi0AndEtaSpectraPPRef","",200,10,1000,800);  // gives the page size
	DrawGammaCanvasSettings( canvasInvYieldTCMOnlyRatioPi0AndEtaSpectraPPRef,  0.15, 0.02, 0.00, 0.00);  

	TPad* padComparisonInvYieldTCMOnlyRatioPi0PPRef 		= new TPad("padComparisonInvYieldTCMOnlyRatioPi0PPRef", "", 0.,0.55, 1., 1,   -1, -1, -2);
	DrawGammaPadSettings( padComparisonInvYieldTCMOnlyRatioPi0PPRef,  0.15, 0.02, 0.02, 0.0);
	padComparisonInvYieldTCMOnlyRatioPi0PPRef->Draw();
	
 	TPad* padComparisonInvYieldTCMOnlyRatioEtaPPRef		= new TPad("padComparisonInvYieldTCMOnlyRatioEtaPPRef", "", 0.,0.0,  1., 0.55,-1, -1, -2);
 	DrawGammaPadSettings( padComparisonInvYieldTCMOnlyRatioEtaPPRef,  0.15, 0.02, 0.00,  0.20);
 	padComparisonInvYieldTCMOnlyRatioEtaPPRef->Draw();
 	
 	
 	padComparisonInvYieldTCMOnlyRatioPi0PPRef->cd();
 	padComparisonInvYieldTCMOnlyRatioPi0PPRef->SetLogx();
	
	
 	
        
	TH2F * histo2DInvYieldTCMOnlyRatioPi0PPRef;
	histo2DInvYieldTCMOnlyRatioPi0PPRef = new TH2F("histo2DInvYieldTCMOnlyRatioPi0PPRef","histo2DRatioAllppreferencesEtaandPi0",1000,.2,25.,1000,0.41,2.2);//1.89
	SetStyleHistoTH2ForGraphs(histo2DInvYieldTCMOnlyRatioPi0PPRef, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{fit}", 0.098,0.098, 0.098,0.098, 1.,0.56, 502, 505); 
	histo2DInvYieldTCMOnlyRatioPi0PPRef->GetYaxis()->SetLabelOffset(0.005);
	histo2DInvYieldTCMOnlyRatioPi0PPRef->GetXaxis()->SetLabelOffset(LabelOffsetLog);
	histo2DInvYieldTCMOnlyRatioPi0PPRef->GetXaxis()->SetTickLength(0.07);
	histo2DInvYieldTCMOnlyRatioPi0PPRef->GetYaxis()-> CenterTitle();
	histo2DInvYieldTCMOnlyRatioPi0PPRef->DrawCopy();
	
	
        
	
	DrawGammaLines(0., 25.,1., 1.,2.0,kGray+2,2);
	
	
	
	  
   
        
        
        
         DrawGammaSetMarkerTGraphAsym(graphRatioPi0ppRefTCMFitSystErr,markerStylePi0,1.2,colorCombYieldPi0 ,colorCombYieldPi0 , 1, kTRUE);  
         graphRatioPi0ppRefTCMFitSystErr->Draw("E2same");
         DrawGammaSetMarkerTGraphAsym(graphRatioPi0ppRefTCMFitStatErr,markerStylePi0,1.2,colorCombYieldPi0, colorCombYieldPi0 );  
         SetEx(graphRatioPi0ppRefTCMFitStatErr,0.); //Set the x-errors bars to 0
         graphRatioPi0ppRefTCMFitStatErr->Draw("Ez,p,same");
   
  
        
        
	
	TLatex * latexInvYieldTCMOnlyRatioPi0PPRef = new TLatex(2.3,2.0,"ALICE") ;
	latexInvYieldTCMOnlyRatioPi0PPRef->SetTextColor(kBlack) ;
	latexInvYieldTCMOnlyRatioPi0PPRef->SetTextSize(TextSize+0.045) ;
	latexInvYieldTCMOnlyRatioPi0PPRef->SetTextFont(Font) ;
	latexInvYieldTCMOnlyRatioPi0PPRef->DrawLatex(2.3,1.8,"p-p, reference, #sqrt{#it{s}} = 5.02 TeV");
	latexInvYieldTCMOnlyRatioPi0PPRef->Draw() ;
	
	
  
	TLegend* legendInvYieldTCMOnlyRatioPi01PPRef;
	legendInvYieldTCMOnlyRatioPi01PPRef   = new TLegend(0.18,0.80,0.38,0.93); 
	legendInvYieldTCMOnlyRatioPi01PPRef->SetFillColor(0);
	legendInvYieldTCMOnlyRatioPi01PPRef->SetFillStyle(0);
	legendInvYieldTCMOnlyRatioPi01PPRef->SetLineColor(0);
	legendInvYieldTCMOnlyRatioPi01PPRef->SetTextFont(Font);
	legendInvYieldTCMOnlyRatioPi01PPRef->SetTextSize(TextSize+0.035);
	legendInvYieldTCMOnlyRatioPi01PPRef->AddEntry(graphRatioPi0ppRefTCMFitSystErr,"#pi^{0}","pef");
        legendInvYieldTCMOnlyRatioPi01PPRef->Draw();
        

	padComparisonInvYieldTCMOnlyRatioEtaPPRef->cd();
	padComparisonInvYieldTCMOnlyRatioEtaPPRef->SetLogx();
      
      
	TH2F * histo2DInvYieldTCMOnlyRatioEtaPPRef;
	histo2DInvYieldTCMOnlyRatioEtaPPRef = new TH2F("histo2DInvYieldTCMOnlyRatioEtaPPRef","histo2DRatioAllppreferencesEtaandEta",1000,.2,25.,1000,0.41,2.1);//1.89
	SetStyleHistoTH2ForGraphs(histo2DInvYieldTCMOnlyRatioEtaPPRef, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{fit}",0.08,0.08,0.08,0.08,1.,0.7, 502, 505); 
	histo2DInvYieldTCMOnlyRatioEtaPPRef->GetYaxis()->SetLabelOffset(0.005);
	histo2DInvYieldTCMOnlyRatioEtaPPRef->GetXaxis()->SetLabelOffset(LabelOffsetLog);
	histo2DInvYieldTCMOnlyRatioEtaPPRef->GetXaxis()->SetTickLength(0.07);
	histo2DInvYieldTCMOnlyRatioEtaPPRef->GetYaxis()-> CenterTitle();
	histo2DInvYieldTCMOnlyRatioEtaPPRef->DrawCopy();

        
        
	
        
        
       DrawGammaSetMarkerTGraphAsym(graphRatioEtappRefTCMFitSystErr,markerStyleEta,1.2,colorCombYieldEta ,colorCombYieldEta, 1, kTRUE);  
       graphRatioEtappRefTCMFitSystErr->Draw("E2same");
       DrawGammaSetMarkerTGraphAsym(graphRatioEtappRefTCMFitStatErr,markerStyleEta,1.2,colorCombYieldEta ,colorCombYieldEta );  
       SetEx(graphRatioEtappRefTCMFitStatErr,0.); //
       graphRatioEtappRefTCMFitStatErr->Draw("Ez,p,same");
        
        
	
	DrawGammaLines(0., 25.,1., 1.,2.0,kGray+2,2);
	
	
	
	
	TLegend* legendInvYieldTCMOnlyRatioEtaPPRef;
	legendInvYieldTCMOnlyRatioEtaPPRef   = new TLegend(0.18,0.83,0.38,0.95); 
	legendInvYieldTCMOnlyRatioEtaPPRef->SetFillColor(0);
	legendInvYieldTCMOnlyRatioEtaPPRef->SetLineColor(0);
	legendInvYieldTCMOnlyRatioEtaPPRef->SetTextFont(Font);
	legendInvYieldTCMOnlyRatioEtaPPRef->SetTextSize((TextSize+0.035)*0.82);
	legendInvYieldTCMOnlyRatioEtaPPRef->AddEntry(graphRatioEtappRefTCMFitSystErr,"#eta","pef");
	legendInvYieldTCMOnlyRatioEtaPPRef->Draw();
	

	canvasInvYieldTCMOnlyRatioPi0AndEtaSpectraPPRef->Update();
	
        canvasInvYieldTCMOnlyRatioPi0AndEtaSpectraPPRef->Print(Form("%s/MesonYields_PPRef_Eta_Pi0_TCM_OnlyRatios.%s",outputDir.Data(),suffix.Data()));
       
        
        
        
        
        
        
       

    /////////////////////////////////////////////////////////////////////////////////
        
    ///Invariant pi0 yield
  
    TCanvas* canvasInvPi0YieldAndRatios = new TCanvas("canvasInvPi0YieldAndRatios","",200,10,1000,1200);  // gives the page size
    DrawGammaCanvasSettings( canvasInvPi0YieldAndRatios, 0.3, 0.02, 0.02, 0.16);
    TPad* padInvPi0YieldAndRatios = new TPad("padInvPi0YieldAndRatios", "", 0., 0.25, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padInvPi0YieldAndRatios, 0.18, 0.02, 0.02, 0.);
    padInvPi0YieldAndRatios->Draw();

    TPad* padRatiosPi0AndModels = new TPad("padRatiosPi0AndModels", "", 0., 0., 1., 0.25,-1, -1, -2);
    DrawGammaPadSettings( padRatiosPi0AndModels, 0.18, 0.02, 0., 0.3);
    padRatiosPi0AndModels->Draw();

    padInvPi0YieldAndRatios->cd(); 
    padInvPi0YieldAndRatios->SetLogx();
    padInvPi0YieldAndRatios->SetLogy();
    TH2F * histoInvPi0YieldAndModels;
    histoInvPi0YieldAndModels = new TH2F("histoInvPi0YieldAndModels","histoInvPi0YieldAndModels",1000,0.2,30.,1000,1.2e-9,20 );
    SetStyleHistoTH2ForGraphs(histoInvPi0YieldAndModels, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2} ",0.035,0.035, 0.035,0.035, 0.8,1.9, 512, 510);
    histoInvPi0YieldAndModels->GetXaxis()->SetLabelOffset(-0.009);
    histoInvPi0YieldAndModels->DrawCopy(); 

  

  
  graphAsymmErrIlkkapPb5020_pi0_ct14_epps16_dss14_scale_err->Draw("c3same");
  graphAsymmErrIlkkapPb5020_pi0_ct14_epps16_dss14_scale_Noerr->Draw("lZXsame");
  
  FitCombPi0->SetLineWidth(lineWidthTsallis); 
  
  
  
  
  
  graphPi0EPOS->Draw("C3same"); //C
  graphPi0EPOSNoErr->Draw("lZXsame");
  
  
  
  
  graphErrMcGillTheoryPion_p_hydro->Draw("C3same");
  graphErrMcGillTheoryPion_p_hydroNoErr->Draw("lXYsame");
  
 


  //if (EPOS){
  //  graphPi0EPOS->Draw("C3same"); //C
 // }
  
  if( CGCPi0 ) {
    
      graphAsymmErrCGCTheoryPi0y0pA5020->Draw("C3same");
    
  }
  
  
  if( HIJINGPi0 ){
    histoHIJINGPi0->Draw("same,hist,l");
  }
  
  if( DPMJetPi0 ){
    
   histoDPMJetPi0->Draw("same,hist,l"); 
  }
  
  graphNLOCalcDSS14InvYieldPi05023GeV_nCTEQ->Draw("3,same");
  graphNLOCalcDSS14InvYieldPi05023GeV_nCTEQ->Draw("lZXsame");
  
  
  
  CombPi0Syst->Draw("E2,same"); 
  CombPi0Stat_noXerrors->Draw("pz,same"); 
  
  
  
  
  
  
   TLegend* legInvPi0YieldAndRatios;
  legInvPi0YieldAndRatios   = new TLegend(0.25,0.08,0.6,0.41); 
  legInvPi0YieldAndRatios->SetFillColor(0);
  legInvPi0YieldAndRatios->SetLineColor(0);
  legInvPi0YieldAndRatios->SetTextFont(Font);
  legInvPi0YieldAndRatios->SetTextSize(TextSize-0.015);
  legInvPi0YieldAndRatios->AddEntry(CombPi0Syst,"#pi^{0}","pef");
  legInvPi0YieldAndRatios->AddEntry(graphPi0EPOS,"EPOS3","fl");
  legInvPi0YieldAndRatios->AddEntry(graphErrMcGillTheoryPion_p_hydro,"VISHNU","fl");
  legInvPi0YieldAndRatios->AddEntry(histoDPMJetPi0,"DPMJet","l");
  legInvPi0YieldAndRatios->AddEntry(histoHIJINGPi0,"HIJING","l");
  legInvPi0YieldAndRatios->AddEntry(graphAsymmErrCGCTheoryPi0y0pA5020,"CGC MV^{#gamma}","l");
  legInvPi0YieldAndRatios->AddEntry(graphAsymmErrIlkkapPb5020_pi0_ct14_epps16_dss14_scale_err,"NLO: EPPS16, DSS14","fl");
  legInvPi0YieldAndRatios->AddEntry(graphNLOCalcDSS14InvYieldPi05023GeV_nCTEQ,"NLO: nCTEQ, DSS14","fl");
  legInvPi0YieldAndRatios->AddEntry(FitCombPi0,"Tsallis fit","l");
  legInvPi0YieldAndRatios->Draw();
  
  
  
 	
  TLatex * ltInvPi0YieldAndRatios = new TLatex(2.,2.1,"ALICE") ;
  ltInvPi0YieldAndRatios->SetTextColor(kBlack) ;
  ltInvPi0YieldAndRatios->SetTextSize(TextSize-0.015) ;
  ltInvPi0YieldAndRatios->SetTextFont(Font) ;
  ltInvPi0YieldAndRatios->DrawLatex(2.,1,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  ltInvPi0YieldAndRatios->Draw() ;

  ltInvPi0YieldAndRatios->Draw() ;

  padRatiosPi0AndModels->cd();
    
  padRatiosPi0AndModels->SetLogx();

  
 
  
  
  TH2F * canvasRatioPi0YieldAndModels;
  canvasRatioPi0YieldAndModels = new TH2F("canvasRatioPi0YieldAndModels","histo2DRatioAllppreferencesEtaandPi0",1000,.2,30.,1000,0.001,4.2);//1.89
  
  //1000,.2,25.,1000,0.001,4.2
  SetStyleHistoTH2ForGraphs(canvasRatioPi0YieldAndModels, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.1,0.1, 0.1,0.1,1.,0.5, 502, 505); 
  canvasRatioPi0YieldAndModels->GetYaxis()->SetLabelOffset(0.005);
  canvasRatioPi0YieldAndModels->GetXaxis()->SetLabelOffset(LabelOffsetLog+0.05);
  canvasRatioPi0YieldAndModels->GetXaxis()->SetTickLength(0.07);
  canvasRatioPi0YieldAndModels->DrawCopy();
  
  
 	
  

  TLine *lineRatioPi0YieldAndModels=new TLine(0.,1.,30.,1.);
  lineRatioPi0YieldAndModels->SetLineColor(kGray+1);
  lineRatioPi0YieldAndModels->Draw("same");

  RatioIlkkaFitPi0scaleerr->Draw("C3same");
  RatioIlkkaFitPi0NoErr->Draw("lEsame");
	
  if (EPOS){
    RatioEPOSFitPi0->Draw("C3same");
    RatioEPOSFitPi0->Draw("lZXsame");
  }
  
  if( CGCPi0){
	    RatioCGCFitPi0->Draw("C3same");
  }
	
  RatioMcGillFitPi0->Draw("C3same");
  RatioMcGillFitPi0->Draw("lZXsame");
  RatioDSS14nCTEQFitPi05023GeVNoErr->Draw("lZXsame");
  RatioDSS14nCTEQFitPi05023GeV->Draw("3,same");

  if( HIJINGPi0 ){
    histoRatioPi0HIJINGToFit->Draw("same,hist,l");  
  }

  if( DPMJetPi0 ){
    	   histoRatioPi0DPMJetToFit->Draw("same,hist,l");  
  }
  
  RatioTsallisCombPi0Syst->Draw("E2,same");
  RatioTsallisCombPi0Stat->Draw("Ez,p,same"); 
	
   
  
  canvasInvPi0YieldAndRatios->Update();
  canvasInvPi0YieldAndRatios->Print(Form("%s/MesonPi0YieldAndModelswithRatios.%s",outputDir.Data(),suffix.Data()));
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  TCanvas* canvasInvPi0YieldAndRatiosb = new TCanvas("canvasInvPi0YieldAndRatiosb","",200,10,1000,1200);  // gives the page size
  DrawGammaCanvasSettings( canvasInvPi0YieldAndRatiosb, 0.3, 0.02, 0.02, 0.16);
  TPad* padInvPi0YieldAndRatiosb = new TPad("padInvPi0YieldAndRatiosb", "", 0., 0.25, 1., 1.,-1, -1, -2);
  DrawGammaPadSettings( padInvPi0YieldAndRatiosb, 0.18, 0.02, 0.02, 0.);
  padInvPi0YieldAndRatiosb->Draw();

  TPad* padRatiosPi0AndModelsb = new TPad("padRatiosPi0AndModelsb", "", 0., 0., 1., 0.25,-1, -1, -2);
  DrawGammaPadSettings( padRatiosPi0AndModelsb, 0.18, 0.02, 0., 0.3);
  padRatiosPi0AndModelsb->Draw();

  padInvPi0YieldAndRatiosb->cd(); 
  padInvPi0YieldAndRatiosb->SetLogx();
  padInvPi0YieldAndRatiosb->SetLogy();
  TH2F * histoInvPi0YieldAndModelsb;
  histoInvPi0YieldAndModelsb = new TH2F("histoInvPi0YieldAndModelsb","histoInvPi0YieldAndModelsb",1000,0.2,30.,1000,1.2e-8,20 );
  SetStyleHistoTH2ForGraphs(histoInvPi0YieldAndModelsb, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2} ",0.040,0.040, 0.040,0.040, 0.8,1.9, 512, 510);
  histoInvPi0YieldAndModelsb->GetXaxis()->SetLabelOffset(-0.009);
  histoInvPi0YieldAndModelsb->DrawCopy(); 

  

  
  
  FitCombPi0->SetLineWidth(lineWidthTsallis); 
  
 

  
  CombPi0Syst->Draw("E2,same"); 
  CombPi0Stat_noXerrors->Draw("pz,same"); 
  
  FitCombPi0->Draw("same");
  
  
  TLegend* legInvPi0YieldAndRatiosb;
  legInvPi0YieldAndRatiosb   = new TLegend(0.25,0.31,0.45,0.41); 
  legInvPi0YieldAndRatiosb->SetFillColor(0);
  legInvPi0YieldAndRatiosb->SetLineColor(0);
  legInvPi0YieldAndRatiosb->SetTextFont(Font);
  legInvPi0YieldAndRatiosb->SetTextSize(TextSize-0.005);
  legInvPi0YieldAndRatiosb->AddEntry(CombPi0Syst,"#pi^{0}","pef");
  legInvPi0YieldAndRatiosb->AddEntry(FitCombPi0,"Tsallis fit","l");
  legInvPi0YieldAndRatiosb->Draw();
  
  TLatex * ltInvPi0YieldAndRatiosb = new TLatex(3.,6.0,"ALICE") ;
  ltInvPi0YieldAndRatiosb->SetTextColor(kBlack) ;
  ltInvPi0YieldAndRatiosb->SetTextSize(TextSize-0.005) ;
  ltInvPi0YieldAndRatiosb->SetTextFont(Font) ;
  ltInvPi0YieldAndRatiosb->DrawLatex(3.,2,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  ltInvPi0YieldAndRatiosb->Draw() ;

  padRatiosPi0AndModelsb->cd();
    
  padRatiosPi0AndModelsb->SetLogx();

  
 
  
  
  TH2F* canvasRatioPi0YieldAndModelsb = new TH2F("canvasRatioPi0YieldAndModelsb","histo2DRatioAllppreferencesEtaandPi0",1000,.2,30.,1000,0.001,2.2);//1.89
  
  SetStyleHistoTH2ForGraphs(canvasRatioPi0YieldAndModelsb, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.125,0.12,0.125,0.125,1.0,0.5, 502, 505); 
  canvasRatioPi0YieldAndModelsb->GetYaxis()->SetLabelOffset(0.005);
  canvasRatioPi0YieldAndModelsb->GetXaxis()->SetLabelOffset(LabelOffsetLog+0.00);
  canvasRatioPi0YieldAndModelsb->GetXaxis()->SetTickLength(0.07);
  canvasRatioPi0YieldAndModelsb->DrawCopy();
  
  
 	
 
  RatioTsallisCombPi0Syst->Draw("E2,same");
  RatioTsallisCombPi0Stat->Draw("Ez,p,same"); 
	
   
  DrawGammaLines(0.2, 30.,1.,1.,2.0,kGray+2,2);
  
  canvasInvPi0YieldAndRatiosb->Update();
  canvasInvPi0YieldAndRatiosb->Print(Form("%s/MesonPi0YieldAndWOModelswithRatios.%s",outputDir.Data(),suffix.Data()));

    
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

  
  TCanvas* canvasInvEtaYieldAndRatios = new TCanvas("canvasInvEtaYieldAndRatios","",200,10,1000,1200);  // gives the page size
  DrawGammaCanvasSettings( canvasInvEtaYieldAndRatios, 0.3, 0.02, 0.02, 0.16);
  TPad* padInvEtaYieldAndRatios = new TPad("padInvEtaYieldAndRatios", "", 0., 0.25, 1., 1.,-1, -1, -2);
  DrawGammaPadSettings( padInvEtaYieldAndRatios, 0.18, 0.02, 0.02, 0.);
  padInvEtaYieldAndRatios->Draw();

  TPad* padRatiosEtaAndModels = new TPad("padRatiosEtaAndModels", "", 0., 0., 1., 0.25,-1, -1, -2);
  DrawGammaPadSettings( padRatiosEtaAndModels, 0.18, 0.02, 0., 0.3);
  padRatiosEtaAndModels->Draw();

  padInvEtaYieldAndRatios->cd(); 
  padInvEtaYieldAndRatios->SetLogx();
  padInvEtaYieldAndRatios->SetLogy();
  TH2F * histoInvEtaYieldAndModels;
  histoInvEtaYieldAndModels = new TH2F("histoInvEtaYieldAndModels","histoInvEtaYieldAndModels",1000,0.2,30.,1000,1.2e-9,20 );
  SetStyleHistoTH2ForGraphs(histoInvEtaYieldAndModels, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2} ",0.035,0.035, 0.035,0.035, 0.8,1.9, 512, 510);
  histoInvEtaYieldAndModels->GetXaxis()->SetLabelOffset(-0.009);
  histoInvEtaYieldAndModels->DrawCopy(); 

  

  
  
  FitCombEta->SetLineWidth(lineWidthTsallis); 
  
  
    while(graphEtaEPOS->GetX()[0] < 0.7)
          graphEtaEPOS->RemovePoint(0); 
    graphEtaEPOS->SetLineWidth(2);
    graphEtaEPOS->SetLineStyle(styleLineEPOS3);
    graphEtaEPOS->SetLineColor(colorEPOSline);
    graphEtaEPOS->SetFillColor(colorEPOS);
    
   TGraphAsymmErrors* graphEtaEPOSNoErr = (TGraphAsymmErrors*)graphEtaEPOS->Clone();
    
    
    
    
  
  
  
  graphEtaEPOS->Draw("C3same"); //C
  graphEtaEPOSNoErr->Draw("lZXsame");
  
  
  //graphErrMcGillTheoryPion_p_hydro->Draw("C3same");
 // graphErrMcGillTheoryPion_p_hydroNoErr->Draw("lXYsame");
  
 


  //if (EPOS){
  //  graphEtaEPOS->Draw("C3same"); //C
 // }
  
  //if( CGCEta ) {
    
  //    graphAsymmErrCGCTheoryEtay0pA5020->Draw("C3same");
    
  //}
  
  
  if( HIJINGPi0 ){
    histoHIJINGEta->Draw("same,hist,l");
    histoHIJINGEta->GetXaxis()->SetRangeUser(0.7,20);
  }
  
  if( DPMJetPi0 ){
    
   histoDPMJetEta->Draw("same,hist,l"); 
   histoDPMJetEta->GetXaxis()->SetRangeUser(0.7,20);
   
  }
  
  
  
  
   DrawGammaSetMarkerTGraphAsym(graphNLOCalcAESSSInvYieldEta5023GeV_nCTEQ, 0, 0, colorDSSnPDFEPPSBand, colorDSSnPDFEPPSBand, widthLinesBoxes, kTRUE, colorDSSnPDFEPPSBand, kTRUE);
   
   graphNLOCalcAESSSInvYieldEta5023GeV_nCTEQ->SetLineStyle(styleLineDSSnPDFEPPS);
   graphNLOCalcAESSSInvYieldEta5023GeV_nCTEQ->Draw("3,same");
   graphNLOCalcAESSSInvYieldEta5023GeV_nCTEQ->Draw("lZXsame");
  
  
    while(graphErrMcGillTheoryEta_p_hydro->GetX()[0] < 0.7)
          graphErrMcGillTheoryEta_p_hydro->RemovePoint(0);   
    graphErrMcGillTheoryEta_p_hydro->SetLineWidth(3);
    graphErrMcGillTheoryEta_p_hydro->SetLineStyle(styleLineMcGill);
    graphErrMcGillTheoryEta_p_hydro->SetLineColor(colorMcGillline);
    graphErrMcGillTheoryEta_p_hydro->SetFillColor(colorMcGill);
   
   
    TGraphAsymmErrors* graphErrMcGillTheoryEta_p_hydroNoErr = (TGraphAsymmErrors*)graphErrMcGillTheoryEta_p_hydro->Clone();
   
    graphErrMcGillTheoryEta_p_hydro->Draw("c3same");
    graphErrMcGillTheoryEta_p_hydroNoErr->Draw("lZXsame");
   

  
  
  CombEtaSyst->Draw("E2,same"); 
  CombEtaStat_noXerrors->Draw("pz,same"); 
  
  
  
  TLegend* legInvEtaYieldAndRatios;
  legInvEtaYieldAndRatios   = new TLegend(0.25,0.08,0.6,0.41); 
  legInvEtaYieldAndRatios->SetFillColor(0);
  legInvEtaYieldAndRatios->SetLineColor(0);
  legInvEtaYieldAndRatios->SetTextFont(Font);
  legInvEtaYieldAndRatios->SetTextSize(TextSize-0.015);
  legInvEtaYieldAndRatios->AddEntry(CombEtaSyst,"#pi^{0}","pef");
  legInvEtaYieldAndRatios->AddEntry(graphEtaEPOS,"EPOS3","fl");
  legInvEtaYieldAndRatios->AddEntry(graphErrMcGillTheoryPion_p_hydro,"VISHNU","fl");
  legInvEtaYieldAndRatios->AddEntry(histoDPMJetEta,"DPMJet","l");
  legInvEtaYieldAndRatios->AddEntry(histoHIJINGEta,"HIJING","l");
  //legInvEtaYieldAndRatios->AddEntry(graphAsymmErrCGCTheoryEtay0pA5020,"CGC MV^{#gamma}","l");
  //legInvEtaYieldAndRatios->AddEntry(graphAsymmErrIlkkapPb5020_pi0_ct14_epps16_dss14_scale_err,"NLO: EPPS16, DSS14","fl");
  legInvEtaYieldAndRatios->AddEntry(graphNLOCalcAESSSInvYieldEta5023GeV_nCTEQ,"NLO: nCTEQ, DSS14","fl");
  legInvEtaYieldAndRatios->AddEntry(FitCombEta,"Tsallis fit","l");
  legInvEtaYieldAndRatios->Draw();
  
  
  
 	
  TLatex * ltInvEtaYieldAndRatios = new TLatex(2.,2.1,"ALICE") ;
  ltInvEtaYieldAndRatios->SetTextColor(kBlack) ;
  ltInvEtaYieldAndRatios->SetTextSize(TextSize-0.015) ;
  ltInvEtaYieldAndRatios->SetTextFont(Font) ;
  ltInvEtaYieldAndRatios->DrawLatex(2.,1,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  ltInvEtaYieldAndRatios->Draw() ;

  ltInvEtaYieldAndRatios->Draw() ;

  padRatiosEtaAndModels->cd();
    
  padRatiosEtaAndModels->SetLogx();

  
 
  
  
  TH2F * canvasRatioEtaYieldAndModels;
  canvasRatioEtaYieldAndModels = new TH2F("canvasRatioEtaYieldAndModels","histo2DRatioAllppreferencesEtaandEta",1000,.2,30.,1000,0.001,4.2);//1.89
  
  //1000,.2,25.,1000,0.001,4.2
  SetStyleHistoTH2ForGraphs(canvasRatioEtaYieldAndModels, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.1,0.1, 0.1,0.1,1.,0.5, 502, 505); 
  canvasRatioEtaYieldAndModels->GetYaxis()->SetLabelOffset(0.005);
  canvasRatioEtaYieldAndModels->GetXaxis()->SetLabelOffset(LabelOffsetLog+0.05);
  canvasRatioEtaYieldAndModels->GetXaxis()->SetTickLength(0.07);
  canvasRatioEtaYieldAndModels->DrawCopy();
  
  
 	
  

  TLine *lineRatioEtaYieldAndModels=new TLine(0.,1.,30.,1.);
  lineRatioEtaYieldAndModels->SetLineColor(kGray+1);
  lineRatioEtaYieldAndModels->Draw("same");

  //RatioIlkkaFitEtascaleerr->Draw("C3same");
  //RatioIlkkaFitEtaNoErr->Draw("lEsame");
	
  if (EPOS){
    RatioEPOSFitEta->Draw("C3same");
    RatioEPOSFitEta->Draw("lZXsame");
  }
  
  //if( CGCEta){
	  //  RatioCGCFitEta->Draw("C3same");
  //}
	
  //RatioMcGillFitEta->Draw("C3same");
  //RatioMcGillFitEta->Draw("lZXsame");
  //RatioDSS14nCTEQFitEta5023GeVNoErr->Draw("lZXsame");
  //RatioDSS14nCTEQFitEta5023GeV->Draw("3,same");

  if( HIJINGPi0 ){
    histoRatioEtaHIJINGToFit->Draw("same,hist,l");  
  }

  if( DPMJetPi0 ){
    	   histoRatioEtaDPMJetToFit->Draw("same,hist,l");  
  }
  
  RatioMcGillFitEta->Draw("C3same");
  RatioMcGillFitEtaNoErr->Draw("lZXsame");
	
	//-AM	
        
        
  RatioAESSSnCTEQFitEta5023GeVNoErr->Draw("lZXsame");
        
  RatioAESSSnCTEQFitEta5023GeV->Draw("3,same");
  
  
  
  RatioTsallisCombEtaSyst->Draw("E2,same");
  RatioTsallisCombEtaStat->Draw("Ez,p,same"); 
	
   
  
  canvasInvEtaYieldAndRatios->Update();
  canvasInvEtaYieldAndRatios->Print(Form("%s/MesonEtaYieldAndModelswithRatios.%s",outputDir.Data(),suffix.Data()));
   
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  TCanvas* canvasInvEtaYieldAndRatiosb = new TCanvas("canvasInvEtaYieldAndRatiosb","",200,10,1000,1200);  // gives the page size
  DrawGammaCanvasSettings( canvasInvEtaYieldAndRatiosb, 0.3, 0.02, 0.02, 0.16);
  TPad* padInvEtaYieldAndRatiosb = new TPad("padInvEtaYieldAndRatiosb", "", 0., 0.25, 1., 1.,-1, -1, -2);
  DrawGammaPadSettings( padInvEtaYieldAndRatiosb, 0.18, 0.02, 0.02, 0.);
  padInvEtaYieldAndRatiosb->Draw();

  TPad* padRatiosEtaAndModelsb = new TPad("padRatiosEtaAndModelsb", "", 0., 0., 1., 0.25,-1, -1, -2);
  DrawGammaPadSettings( padRatiosEtaAndModelsb, 0.18, 0.02, 0., 0.3);
  padRatiosEtaAndModelsb->Draw();

  padInvEtaYieldAndRatiosb->cd(); 
  padInvEtaYieldAndRatiosb->SetLogx();
  padInvEtaYieldAndRatiosb->SetLogy();
  
  TH2F* histoInvEtaYieldAndModelsb = new TH2F("histoInvEtaYieldAndModelsb","histoInvEtaYieldAndModelsb",1000,0.5,30.,1000,1.2e-8,1 );
  SetStyleHistoTH2ForGraphs(histoInvEtaYieldAndModelsb, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2} ",0.040,0.040, 0.040,0.040, 0.8,1.9, 512, 510);
  histoInvEtaYieldAndModelsb->GetXaxis()->SetLabelOffset(-0.009);
  histoInvEtaYieldAndModelsb->DrawCopy(); 

  
  
  FitCombEta->SetLineWidth(lineWidthTsallis); 

  
  CombEtaSyst->Draw("E2,same"); 
  CombEtaStat_noXerrors->Draw("pz,same"); 
  
  FitCombEta->Draw("same");
  
  
  
  TLegend* legInvEtaYieldAndRatiosb;
  legInvEtaYieldAndRatiosb   = new TLegend(0.25,0.31,0.45,0.41); 
  legInvEtaYieldAndRatiosb->SetFillColor(0);
  legInvEtaYieldAndRatiosb->SetLineColor(0);
  legInvEtaYieldAndRatiosb->SetTextFont(Font);
  legInvEtaYieldAndRatiosb->SetTextSize(TextSize-0.005);
  legInvEtaYieldAndRatiosb->AddEntry(CombEtaSyst,"#pi^{0}","pef");
  legInvEtaYieldAndRatiosb->AddEntry(FitCombEta,"Tsallis fit","l");
  legInvEtaYieldAndRatiosb->Draw();
  
  
  
 	
  TLatex * ltInvEtaYieldAndRatiosb = new TLatex(4.5,0.3,"ALICE") ;
  ltInvEtaYieldAndRatiosb->SetTextColor(kBlack) ;
  ltInvEtaYieldAndRatiosb->SetTextSize(TextSize-0.005) ;
  ltInvEtaYieldAndRatiosb->SetTextFont(Font) ;
  ltInvEtaYieldAndRatiosb->DrawLatex(4.5,0.08,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  ltInvEtaYieldAndRatiosb->Draw() ;


  padRatiosEtaAndModelsb->cd();
    
  padRatiosEtaAndModelsb->SetLogx();

  
 
  
  
  TH2F * canvasRatioEtaYieldAndModelsb;
  canvasRatioEtaYieldAndModelsb = new TH2F("canvasRatioEtaYieldAndModelsb","histo2DRatioAllppreferencesEtaandEta",1000,.5,30.,1000,0.001,2.2);//1.89
  
  //1000,.2,25.,1000,0.001,4.2
  SetStyleHistoTH2ForGraphs(canvasRatioEtaYieldAndModelsb, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.125,0.125, 0.125,0.125,1.,0.5, 502, 505); 
  canvasRatioEtaYieldAndModelsb->GetYaxis()->SetLabelOffset(0.005);
  canvasRatioEtaYieldAndModelsb->GetXaxis()->SetLabelOffset(LabelOffsetLog+0.00);
  canvasRatioEtaYieldAndModelsb->GetXaxis()->SetTickLength(0.07);
  canvasRatioEtaYieldAndModelsb->DrawCopy();
 
  
  
 	
  

  TLine *lineRatioEtaYieldAndModelsb=new TLine(0.,1.,30.,1.);
  lineRatioEtaYieldAndModelsb->SetLineColor(kGray+1);
  lineRatioEtaYieldAndModelsb->Draw("same");

  
  

  
  RatioTsallisCombEtaSyst->Draw("E2,same");
  RatioTsallisCombEtaStat->Draw("Ez,p,same"); 
	
   
  
  canvasInvEtaYieldAndRatiosb->Update();
  canvasInvEtaYieldAndRatiosb->Print(Form("%s/MesonEtaYieldAndWOModelswithRatios.%s",outputDir.Data(),suffix.Data()));
  
  ///////////////////////////////////////////////////////////////////////////////////////////
        
   
  
  
  
  
  
  
   
        
        
        
	
	
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

	
        //DrawGammaLines(0.6,25.,1.,1.,2.0,kGray+2,2);
        DrawGammaLines(0., 25.,1.,1.,2.0,kGray+2,2);

	
	if( HIJINGPi0 ){
           SetStyleHisto(histoRatioPi0HIJINGToFit,3, styleLineHIJING, colorHIJING );  
           histoRatioPi0HIJINGToFit->Draw("same,hist,l");  
	}

	if( DPMJetPi0 ){
	   SetStyleHisto(histoRatioPi0DPMJetToFit,3, styleLineDPMJet, colorDPMJet );  
	   histoRatioPi0DPMJetToFit->Draw("same,hist,l");  
	}
  
	
	RatioTsallisCombPi0Syst->Draw("E2same");
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
        
        
        TFile fPaperPlotsFinal(fileNamepPbOutput.Data(),"UPDATE");
        
        //EPOS
        RatioEPOSFitPi0->Write("RatioEPOSFitPi0");
        RatioEPOSFitEta->Write("RatioEPOSFitEta");
        
  //CGC
  
        RatioCGCFitPi0->Write("RatioCGCFitPi0");
  
  //McGill
  
        RatioMcGillFitPi0->Write("RatioMcGillFitPi0");
        RatioMcGillFitEta->Write("RatioMcGillFitEta");
  
        RatioIlkkaFitPi0->Write("RatioIlkkaFitPi0");
        RatioIlkkaFitPi0NoErr->Write("RatioIlkkaFitPi0NoErr");
        RatioIlkkaFitPi0scaleerr->Write("RatioIlkkaFitPi0scaleerr");
        
        histoRatioPi0DPMJetToFit->Write("histoRatioPi0DPMJetToFit");
        histoRatioPi0HIJINGToFit->Write("histoRatioPi0HIJINGToFit");
   
        histoRatioEtaDPMJetToFit->Write("histoRatioEtaDPMJetToFit");
        histoRatioEtaHIJINGToFit->Write("histoRatioEtaHIJINGToFit");
        
        EtaPi0Stat_div_mT->Write("EtaPi0Stat_div_mT");
        EtaPi0Syst_div_mT->Write("EtaPi0Syst_div_mT");
   
   
        CombinedPi0RpPbSystErr->Write("CombinedPi0RpPbSystErr");
        CombinedPi0RpPbStatErr->Write("CombinedPi0RpPbStatErr");
 
    
        fPaperPlotsFinal.Close();

        
        
}


TGraphErrors* Eta2Pi0dAu(){
    
    
//     Centrality & $\pt$ (GeV/$c$) & $\eta/\pi^0$ & tot. err. & stat. err. + error A & error B & error C \\\hline
//        & 2.25 & 0.420 & 0.038 & 0.028 & 0.025 & 0 \\
//        & 2.75 & 0.472 & 0.044 & 0.033 & 0.028 & 0 \\
//        & 3.25 & 0.383 & 0.045 & 0.039 & 0.023 & 0 \\
//        & 3.75 & 0.472 & 0.034 & 0.018 & 0.028 & 0 \\
//        & 4.25 & 0.478 & 0.033 & 0.017 & 0.029 & 0 \\
//        & 4.75 & 0.483 & 0.035 & 0.020 & 0.029 & 0 \\
// 0-88\% & 5.25 & 0.465 & 0.037 & 0.025 & 0.028 & 0 \\
//  (MB)  & 5.75 & 0.510 & 0.043 & 0.030 & 0.031 & 0 \\
//        & 6.5 & 0.552 & 0.048 & 0.034 & 0.033 & 0 \\
//        & 7.5 & 0.478 & 0.070 & 0.064 & 0.029 & 0 \\
//        & 8.5 & 0.499 & 0.092 & 0.087 & 0.030 & 0 \\
//        & 9.5 & 0.677 & 0.118 & 0.111 & 0.041 & 0 \\
//        & 11 & 0.609 & 0.124 & 0.119 & 0.037 & 0 \\ \hline
//     
     Int_t nPoints = 13;
    
     Double_t xval[] = { 2.25, 2.75, 3.25, 3.75, 4.25, 
                         4.75, 5.25, 5.75, 6.5,  7.5,  
                         8.5,  9.5,  11};
                         
     Double_t yval[] = { 0.420, 0.472, 0.383, 0.472, 0.478,
                         0.483, 0.465, 0.510, 0.552, 0.478,
                         0.499, 0.677, 0.609 };
                         
//      Double_t xErr[] = { 0.25, 0.25, 0.25, 0.25, 0.25,
//                          0.25, 0.25, 0.25, 0.5,  0.5,
//                          0.5,  0.5, 1.0 };
//                          
      Double_t xErr[] = { 0.0, 0.0, 0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0 };
    
     Double_t yErr[] = { 0.038, 0.044, 0.045, 0.034, 0.033,
                         0.035, 0.037, 0.043, 0.048, 0.070,
                         0.092, 0.118, 0.124 };
     
    
 
    TGraphErrors* graphEta2Pi0Ratio = new TGraphErrors(nPoints,xval,yval,xErr,yErr);
    
    
    
    return graphEta2Pi0Ratio;
    
    
}
