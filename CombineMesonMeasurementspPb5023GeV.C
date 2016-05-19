/****************************************************************************************************************************
 ******         provided by Gamma Conversion Group, PWG4,                                                     *****
 ******        Ana Marin, marin@physi.uni-heidelberg.de                                                    *****
 ******           Kathrin Koch, kkoch@physi.uni-heidelberg.de                                                     *****
 ******        Friederike Bock, friederike.bock@cern.ch                                                    *****
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


void DrawPad(Int_t iPad,Double_t a,Double_t b, Double_t c, Double_t d){
  
  Double_t bottom=1. ;
  TPad * sp;
  if(iPad==0)  sp = new TPad(Form("spectrum%d",iPad),"",0.,0,1.,bottom*0.22*(iPad+1)+0.1);
  else  sp = new TPad(Form("spectrum%d",iPad),"",0.,bottom*0.22*iPad+0.1,1.,bottom*0.22*(iPad+1)+0.1);
  //  TPad * sp = new TPad(Form("spectrum%d",iPad),"",0.,bottom*0.25*iPad,1.,bottom*0.25*(iPad+1));
  sp->Draw();
  sp->cd();
  sp->Range(0,0,1,1);
  sp->SetFillColor(0);
  sp->SetFillStyle(0);
  sp->SetBorderSize(1);
  sp->SetTopMargin(a);
  sp->SetBottomMargin(b);
  sp->SetLeftMargin(c);
  sp->SetRightMargin(d);
  sp->SetLogx();
  sp->SetGridy();
  sp->cd() ;
}  

void Draw2Pads(Int_t iPad,Double_t a,Double_t b, Double_t c, Double_t d){
  
  Double_t bottom=1. ;
  TPad * sp;
  if(iPad==0)  sp = new TPad(Form("spectrum%d",iPad),"",0.,0.,1.,0.53);
  else  sp = new TPad(Form("spectrum%d",iPad),"",0.,0.53,1.,0.95);
  //  TPad * sp = new TPad(Form("spectrum%d",iPad),"",0.,bottom*0.25*iPad,1.,bottom*0.25*(iPad+1));
  sp->Draw();
  sp->cd();
  sp->Range(0,0,1,1);
  sp->SetFillColor(0);
  sp->SetFillStyle(0);
  sp->SetBorderSize(1);
  sp->SetTopMargin(a);
  sp->SetBottomMargin(b);
  sp->SetLeftMargin(c);
  sp->SetRightMargin(d);
  sp->SetLogx();
  sp->SetGridy();
  sp->cd() ;
}
struct SysErrorConversion {
  Double_t value;
  Double_t error;
  //    TString name;
};


const Int_t  Ntotal = 31;
const Int_t  NtotalLow = 26;
const Int_t  nPtLimits = Ntotal+1;
const Int_t  nPtLimitsLow = NtotalLow+1;
const Int_t  NtotalEta = 16;
const Int_t  nPtLimitsEta = NtotalEta+1;
const Int_t  NtotalEtaPi0Ratio = 15;
const Int_t  nPtLimitsEtaPi0Ratio = NtotalEtaPi0Ratio+1;

void CombineMesonMeasurementspPb5023GeV(TString FittingType = "Tsallis",Bool_t LowPtCut=kFALSE){    

  TString date                                        = ReturnDateString();
    
  gROOT->Reset();
  gROOT->SetStyle("Plain");
    
  StyleSettingsThesis();
  SetPlotStyle();
    
  TString suffix                                      = "eps";
    
  TString dateForOutput                               = ReturnDateStringForOutput();
  TString outputDir                                   = Form("CombinepPbSpectra/%s/%s/%s",dateForOutput.Data(),FittingType.Data(),suffix.Data());
  if (LowPtCut) outputDir                             = Form("CombinepPbSpectra/%s_LowPtCut/%s/%s",dateForOutput.Data(),FittingType.Data(),suffix.Data());
  gSystem->Exec("mkdir -p "+outputDir);
  fstream fFits;   
  TString nameFits = Form("%s/FittingParameters_%s_%s.dat",outputDir.Data(),dateForOutput.Data(),FittingType.Data());
  fFits.open(nameFits.Data(), ios::out);	
  //___________________________________ Declaration of files _____________________________________________
  TString fileNameNeutralPionDalitz                   = "ExternalInputpPb/PCM/data_PCMResults_Dalitz_pPb_20150806.root";
   TString fileNameNeutralPionPCM                      = "ExternalInputpPb/PCM/data_PCMResults_pPb_20151111_standard_CatErrors.root";
      TString fileNameNeutralPionEMCal                    = "ExternalInputpPb/EMCAL/data_EMCalEMCalResults_160503_pPb.root"; 
      // TString fileNameNeutralPionEMCal                    = "ExternalInputpPb/EMCAL/data_EMCalEMCalResults_160215_pPb_5023_Mike.roo";
  TString fileNameNeutralPionPHOS                     = "ExternalInputpPb/PHOS/data_PHOSResults_pPb_20160208.root";
 

  TString nameHistoPHOS                               = "hCor_stat";
  TString nameHistoPHOSSysErrors                      = "hCor_syst";
 
  TString nameHistoEMCal                              = "CorrectedYieldPi0";
  TString nameHistoEMCalSysErrors                     = "Pi0SystError";
    
    
  TString collisionSystempPb                          = "p-Pb #sqrt{s_{NN}} = 5.02 TeV"; 
    
    
  Double_t pTLimits[nPtLimits]                        = { 0.3, 0.4, 0.5, 0.6, 0.7, 
							  0.8, 1.0, 1.2, 1.4, 1.6, 
							  1.8, 2.0, 2.2, 2.4, 2.6,
							  2.8, 3.0, 3.2, 3.4, 3.6,
							  3.8, 4.0, 4.5, 5.0, 5.5,
							  6.0, 7.0, 8.0, 10.0,12.0,
							  16.0, 20.0};
    
  Double_t pTLimitsLow[nPtLimitsLow]                        = {  1.0, 1.2, 1.4, 1.6, 
							  1.8, 2.0, 2.2, 2.4, 2.6,
							  2.8, 3.0, 3.2, 3.4, 3.6,
							  3.8, 4.0, 4.5, 5.0, 5.5,
							  6.0, 7.0, 8.0, 10.0,12.0,
							  16.0, 20.0};

  Double_t pTLimitsEta[nPtLimits]                        = { 0.7,0.9,1.1,1.4,1.8,
							     2.2,2.6,3.0,3.5,4.0,
							     5.0,6.0,8.0,10.0,12.0,
							     16.0,20.0};
  Double_t pTLimitsEtaPi0Ratio[nPtLimits]                        = { 0.7,0.9,1.1,1.4,1.8,
							     2.2,2.6,3.0,3.5,4.0,
							     5.0,6.0,8.0,10.0,12.0,
							     16.0};
   
   Int_t offSets[11]                                   =  { -1, 5, -1, 0, 0,  2, 0, 0, 0,  0, 0};
  Int_t offSetsSys[11]                                =  {  0, 6,  8, 0, 0, 3, 0, 0, 0, 0, 0};  
   Int_t offSetsLow[11]                                   =  { -7, -1, -7, 0, 0,  -4, 0, 0, 0,  0, 0};
  Int_t offSetsSysLow[11]                                =  {  -6, 0,  2, 0, 0, -3, 0, 0, 0, 0, 0};  
  Int_t offSetsEta[11]                                   =  { -3, 0, -3, 0, 0,  0, 0, 0, 0,  0, 0};
  Int_t offSetsSysEta[11]                                =  {  0, 0,  5, 0, 0, 0, 0, 0, 0, 0, 0};
  Int_t offSetsEtaPi0Ratio[11]                                   =  { 0, 0, 5, 0, 0,  0, 0, 0, 0,  0, 0};
  Int_t offSetsSysEtaPi0Ratio[11]                                =  {  0, 0,  5, 0, 0, 0, 0, 0, 0, 0, 0};
    
    
  TH1D* statErrorCollection[11];
  for (Int_t i = 0; i< 11; i++){
    statErrorCollection[i]                          = NULL;
  }    
    
  TGraphAsymmErrors* sysErrorCollection[11];
  for (Int_t i = 0; i< 11; i++){
    sysErrorCollection[i]                           = NULL;
  }    
    
  TH1D* statErrorCollectionEta[11];
  for (Int_t i = 0; i< 11; i++){
    statErrorCollectionEta[i]                          = NULL;
  }    
    
  TH1D* statErrorCollectionEtaPi0Ratio[11];
  for (Int_t i = 0; i< 11; i++){
    statErrorCollectionEtaPi0Ratio[i]                          = NULL;
  } 
    
  TGraphAsymmErrors* sysErrorCollectionEta[11];
  for (Int_t i = 0; i< 11; i++){
    sysErrorCollectionEta[i]                           = NULL;
  }      
    
  TGraphAsymmErrors* sysErrorCollectionEtaPi0Ratio[11];
  for (Int_t i = 0; i< 11; i++){
    sysErrorCollectionEtaPi0Ratio[i]                           = NULL;
  }    
  // **************************************************************************************
  // ****************************** Reading Dalitz ****************************************
  // **************************************************************************************
  TFile* fileNeutralPionsDalitz                       = new TFile(fileNameNeutralPionDalitz.Data());
  TDirectory* directoryDalitzPi0pPb                   = (TDirectory*) fileNeutralPionsDalitz->GetDirectory("Pi0_pPb_5.023TeV_0-100%");
  if( ! directoryDalitzPi0pPb ){
    cout<<"Dalitz: The directory Pi0_pPb_5.023TeV_0-100% does not exist "<<endl; 
    return;
  }

  TH1D* histoDalitzYieldPi0pPb                        = (TH1D*)directoryDalitzPi0pPb->Get("CorrectedYieldPi0");
  TGraphAsymmErrors* graphDalitzYieldPi0pPbSystErr    = (TGraphAsymmErrors*)directoryDalitzPi0pPb->Get("Pi0SystError");
    
  TGraphAsymmErrors* graphRatioCombDalitzSys          = (TGraphAsymmErrors*)graphDalitzYieldPi0pPbSystErr->Clone();
  TGraphAsymmErrors* temp02                           = new TGraphAsymmErrors(histoDalitzYieldPi0pPb);

  cout<<"Dalitz systematic"<<endl;
  graphRatioCombDalitzSys->Print();
  cout<<"Dalitz statistic"<<endl;
  temp02->Print();
    
  // **************************************************************************************
  // ****************************** Reading PCM *******************************************
  // **************************************************************************************    
  TFile* fileNeutralPionPCM                           = new TFile(fileNameNeutralPionPCM.Data());
  TDirectory* directoryPCMPi0pPb                      = (TDirectory*) fileNeutralPionPCM->GetDirectory("Pi0_pPb_5.023TeV_0-100%");    
  if( ! directoryPCMPi0pPb ){
    cout<<"PCM: The directory Pi0_pPb_5.023TeV_0-100% does not exist "<<endl; 
    return; 
  }
  TDirectory*	directoryPCMEtapPb 			= (TDirectory*)fileNeutralPionPCM->Get("Eta_pPb_5.023TeV_0-100%"); 
  if( ! directoryPCMEtapPb ){
    cout<<"PCM: The directory Eta_pPb_5.023TeV_0-100% does not exist "<<endl; 
    return; 
  }
  TDirectory*	directoryPCMPi0EtaBinningpPb 			= (TDirectory*)fileNeutralPionPCM->Get("Pi0EtaBinning_pPb_5.023TeV_0-100%"); 
  if( ! directoryPCMPi0EtaBinningpPb ){
    cout<<"PCM: The directory Pi0EtaBinning_pPb_5.023TeV_0-100% does not exist "<<endl; 
    return; 
  }
  TH1D* histoPCMYieldPi0pPb                           = (TH1D*)directoryPCMPi0pPb->Get("CorrectedYieldPi0");
  TGraphAsymmErrors* graphPCMYieldPi0pPbSystErr       = (TGraphAsymmErrors*)directoryPCMPi0pPb->Get("Pi0SystError");
  TGraphAsymmErrors* temp01                           = new TGraphAsymmErrors(histoPCMYieldPi0pPb);

  TH1D*  histoPCMYieldEtapPb 				= (TH1D*)directoryPCMEtapPb->Get("CorrectedYieldEta");  
  TGraphAsymmErrors*  graphPCMYieldEtapPb 				= new TGraphAsymmErrors(histoPCMYieldEtapPb);
  TGraphAsymmErrors*	graphPCMYieldEtapPbSystErr	= (TGraphAsymmErrors*)directoryPCMEtapPb->Get("EtaSystError");	
  TGraphAsymmErrors*	graphPCMYieldEtapPbSystErrA	= (TGraphAsymmErrors*)directoryPCMEtapPb->Get("EtaSystErrorA");	
  TH1D*	histoPCMEtaPi0RatiopPb 			= (TH1D*)directoryPCMEtapPb->Get("EtatoPi0Ratio");
  TGraphAsymmErrors*	graphPCMEtaPi0RatiopPbSystErr 	= (TGraphAsymmErrors*)directoryPCMEtapPb->Get("EtatoPi0RatioSys");

  TH1D*  histoPCMYieldPi0EtaBinningpPb 				= (TH1D*)directoryPCMPi0EtaBinningpPb->Get("CorrectedYieldPi0EtaBinning");
   TGraphAsymmErrors*  graphPCMYieldPi0EtaBinningpPb 				= new TGraphAsymmErrors(histoPCMYieldPi0EtaBinningpPb);
  // TGraphAsymmErrors*	graphPCMYieldPi0EtaBinningpPbSystErr	= (TGraphAsymmErrors*)directoryPCMEtapPb->Get("Pi0EtaBinningSystError");
  
    
  cout<<"PCM systematic"<<endl;
  graphPCMYieldPi0pPbSystErr->Print();
  cout<<"PCM statistic"<<endl;
  temp01->Print();
  cout<<"PCM Eta systematic"<<endl;
  graphPCMYieldEtapPbSystErr->Print();   
  cout<<"PCM Pi0Eta systematicA"<<endl;
  graphPCMYieldPi0EtaBinningpPb->RemovePoint(0); 
  graphPCMYieldPi0EtaBinningpPb->RemovePoint(0); 
  graphPCMYieldPi0EtaBinningpPb->RemovePoint(0); 
  graphPCMYieldPi0EtaBinningpPb->Print(); 
  cout<<"PCM Eta "<<endl;
  graphPCMYieldEtapPb->RemovePoint(0);
  graphPCMYieldEtapPb->RemovePoint(0);
  graphPCMYieldEtapPb->RemovePoint(0);
  graphPCMYieldEtapPb->Print(); 
  cout<<"PCM Eta systematicA"<<endl;
  graphPCMYieldEtapPbSystErrA->Print(); 
  TGraphAsymmErrors* graphRatioCombPCMSys             = (TGraphAsymmErrors*)graphPCMYieldPi0pPbSystErr->Clone();;
  TGraphAsymmErrors* graphRatioCombEtaPCMSys             = (TGraphAsymmErrors*)graphPCMYieldEtapPbSystErr->Clone();;
    
    
  // **************************************************************************************
  // ******************************* Reading PHOS *****************************************
  // **************************************************************************************    
  TFile* fileNeutralPionPHOS                          = new TFile(fileNameNeutralPionPHOS);
  //TDirectory*    directoryPHOSPi0pPb =                 (TDirectory*)fileNeutralPionPHOS->Get("Pi0_pPb_5.023TeV_0-100%"); 
  if( ! fileNeutralPionPHOS ){
    cout<<"PHOS: The file fileNeutralPionPHOS does not exist "<<endl; 
    return;  
  }
   
  TH1D* histoPHOSYieldPi0pPbStat                      = (TH1D*)fileNeutralPionPHOS->Get(nameHistoPHOS.Data());
  TH1D* histoPHOSYieldPi0pPbSyst                      = (TH1D*)fileNeutralPionPHOS->Get(nameHistoPHOSSysErrors.Data());
  TGraphAsymmErrors* graphPHOSYieldPi0pPbSystErr      = new TGraphAsymmErrors(histoPHOSYieldPi0pPbSyst);
  graphPHOSYieldPi0pPbSystErr->RemovePoint(0);
  TGraphAsymmErrors* temp03                           = new TGraphAsymmErrors(histoPHOSYieldPi0pPbStat);
    
  cout<<"Systematic PHOS "<<endl;
  graphPHOSYieldPi0pPbSystErr->Print();
  cout<<"Statistic PHOS  "<<endl;
  temp03->Print();
  TH1D * histoRatioCombPHOS                           = (TH1D*)histoPHOSYieldPi0pPbSyst->Clone();
    
  // **************************************************************************************
  // ******************************** Reading EMCal ***************************************
  // **************************************************************************************
  TFile* fileNeutralPionEMCal                         = new TFile(fileNameNeutralPionEMCal);
  TDirectory* directoryEMCalPi0pPb                    = (TDirectory*)fileNeutralPionEMCal->Get("Pi0_pPb_5.023TeV_0-100%");
  if( ! directoryEMCalPi0pPb ){
    cout<<"EMCal: The directory Pi0_pPb_5.023TeV_0-100% does not exist "<<endl; 
    return;  
  }
  TDirectory*	directoryEMCalEtapPb 	      	= (TDirectory*)fileNeutralPionEMCal->Get("Eta_pPb_5.023TeV_0-100%"); 
  if( ! directoryEMCalEtapPb  ){
    cout<<"EMCal: The directory Eta_pPb_5.023TeV_0-100% does not exist "<<endl; 
    return;  
  }  TDirectory*	directoryEMCalPi0EtaBinningpPb 	      	= (TDirectory*)fileNeutralPionEMCal->Get("Pi0EtaBinning_pPb_5.023TeV_0-100%"); 
  if( ! directoryEMCalPi0EtaBinningpPb  ){
    cout<<"EMCal: The directory Pi0EtaBinning_pPb_5.023TeV_0-100% does not exist "<<endl; 
    return;  
  }
  //Pi0
  TH1D* histoEMCalYieldPi0pPbStat                     = (TH1D*)directoryEMCalPi0pPb->Get(nameHistoEMCal.Data());
  TGraphAsymmErrors* graphEMCalYieldPi0pPbSystErr     = (TGraphAsymmErrors*)directoryEMCalPi0pPb->Get(nameHistoEMCalSysErrors.Data());
  TGraphAsymmErrors* graphRatioCombEMCalSys           = (TGraphAsymmErrors*)graphEMCalYieldPi0pPbSystErr->Clone();
  //Eta
  TH1D*	histoEMCalYieldEtapPbStat 			= (TH1D*)directoryEMCalEtapPb->Get("CorrectedYieldEta");
  TGraphAsymmErrors*	graphEMCalYieldEtapPbSystErr		= (TGraphAsymmErrors*)directoryEMCalEtapPb->Get("EtaSystError");
  TGraphAsymmErrors* graphEMCalYieldEtapPb=new TGraphAsymmErrors(histoEMCalYieldEtapPbStat); 
  TGraphAsymmErrors* graphRatioCombEtaEMCalSys           = (TGraphAsymmErrors*)graphEMCalYieldEtapPbSystErr->Clone();
  //Pi0EtaBinning
  TH1D* histoEMCalYieldPi0EtaBinningpPbStat  = (TH1D*)directoryEMCalPi0EtaBinningpPb->Get("CorrectedYieldPi0EtaBinning");
  TGraphAsymmErrors*   graphEMCalYieldPi0EtaBinningpPb=new TGraphAsymmErrors(histoEMCalYieldPi0EtaBinningpPbStat);
  //EtaPi0Ratio	
  TH1D*	histoEMCalEtaPi0RatiopPb 			= (TH1D*)directoryEMCalEtapPb->Get("EtatoPi0Ratio");
  TGraphAsymmErrors*	graphEMCalEtaPi0RatiopPbSystErr 	        = (TGraphAsymmErrors*)directoryEMCalEtapPb->Get("EtatoPi0RatioSys"); 
  cout<<"EMCal Eta systematic"<<endl;
  graphEMCalYieldEtapPbSystErr->Print();
  cout<<"EMCal Eta systematictatistic"<<endl;
  graphEMCalYieldEtapPbSystErr->Print();       
  // **************************************************************************************
  // ********************************* Combine Pi0 spectra ********************************
  // **************************************************************************************
    
  TString fileNameOutputWeightingPi0                  = Form("%s/WeightingOld.dat",outputDir.Data());

  statErrorCollection[0]          = (TH1D*)histoPCMYieldPi0pPb->Clone("statErrPCMPi0");
  statErrorCollection[1]          = (TH1D*)histoPHOSYieldPi0pPbStat->Clone("statErrPHOSPi0");
  statErrorCollection[2]          = (TH1D*)histoEMCalYieldPi0pPbStat->Clone("statErrEMCalPi0");
  statErrorCollection[5]          = (TH1D*)histoDalitzYieldPi0pPb->Clone("statErrDalitzPi0");
    
  sysErrorCollection[0]           = (TGraphAsymmErrors*)graphPCMYieldPi0pPbSystErr->Clone("sysErrPCMPi0");
  sysErrorCollection[1]           = (TGraphAsymmErrors*)graphPHOSYieldPi0pPbSystErr->Clone("sysErrPHOSPi0");
  sysErrorCollection[2]           = (TGraphAsymmErrors*)graphEMCalYieldPi0pPbSystErr->Clone("sysErrEMCalPi0");
  sysErrorCollection[5]           = (TGraphAsymmErrors*)graphDalitzYieldPi0pPbSystErr->Clone("sysErrDalitzPi0");
    
    
    
    
  TGraphAsymmErrors* graphCombPi0InvCrossSectionStatpPb5023GeV= NULL;
  TGraphAsymmErrors* graphCombPi0InvCrossSectionSyspPb5023GeV = NULL;
    
  TGraphAsymmErrors* graphCombPi0InvCrossSectionTotpPb5023GeV=NULL;
  if(!LowPtCut)graphCombPi0InvCrossSectionTotpPb5023GeV= CombinePtPointsSpectraFullCorrMat(    statErrorCollection,    sysErrorCollection,     
												      pTLimits, Ntotal,
												      offSets, offSetsSys,
												      graphCombPi0InvCrossSectionStatpPb5023GeV, graphCombPi0InvCrossSectionSyspPb5023GeV,
												      fileNameOutputWeightingPi0,1
												      );
  else graphCombPi0InvCrossSectionTotpPb5023GeV= CombinePtPointsSpectraFullCorrMat(    statErrorCollection,    sysErrorCollection,     
												      pTLimitsLow, NtotalLow,
												      offSetsLow, offSetsSysLow,
												      graphCombPi0InvCrossSectionStatpPb5023GeV, graphCombPi0InvCrossSectionSyspPb5023GeV,
												      fileNameOutputWeightingPi0,1
												      );
  graphCombPi0InvCrossSectionStatpPb5023GeV->Print();
    
  TGraphAsymmErrors* graphInvYieldPi0CombpPb5023GeVStaClone   = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionStatpPb5023GeV->Clone();
  TGraphAsymmErrors* graphInvYieldPi0CombpPb5023GeVSysClone   = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionSyspPb5023GeV->Clone();
  TGraphAsymmErrors* graphInvYieldPi0CombpPb5023GeVTotClone   = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionTotpPb5023GeV->Clone();
    
  TGraphAsymmErrors* graphRatioCombCombFit                    = (TGraphAsymmErrors*) graphInvYieldPi0CombpPb5023GeVTotClone ->Clone(); 
  TGraphAsymmErrors* graphRatioCombCombFitSta                 = (TGraphAsymmErrors*) graphInvYieldPi0CombpPb5023GeVStaClone ->Clone();
  TGraphAsymmErrors* graphRatioCombCombFitSys                 = (TGraphAsymmErrors*) graphInvYieldPi0CombpPb5023GeVSysClone ->Clone(); 

    
    
  // **************************************************************************************
  // ************************* Plotting yield of different systems ************************
  // **************************************************************************************
  TCanvas* canvasCompYieldpPbInd = new TCanvas("canvasCompYieldpPbInd","",200,10,500,400);  // gives the page size
  DrawGammaCanvasSettings( canvasCompYieldpPbInd, 0.15, 0.02, 0.02, 0.12);
    
  canvasCompYieldpPbInd->SetLogx();
  canvasCompYieldpPbInd->SetLogy();
  TH2F * histo2DCompCombined;
  histo2DCompCombined = new TH2F("histo2DCompCombined","histo2DCompCombined",1000,0.3,30.,1000,1.2e-9,30.   );
  SetStyleHistoTH2ForGraphs(histo2DCompCombined, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.03,0.04, 0.03,0.04, 1.0,1.4, 512, 508);
  // SetStyleHistoTH2ForGraphs(histo2DCompCombined, "#it{p}_{T} (GeV/#it{c})","arbitrary units", 0.03,0.04, 0.03,0.04, 1.0,1.1, 512, 508);
  //   histo2DCompCombined->GetXaxis()->SetRangeUser(0.,30.);
  //  histo2DCompCombined->GetYaxis()->SetLabelColor(0);
  //  histo2DCompCombined->GetYaxis()->SetRangeUser(0.1,2.1);
  histo2DCompCombined->DrawCopy();
 
  TF1* fitCombPi0pPb5023GeVPt;
  Double_t *ParameterspPb  = new Double_t[5];
  Double_t minPt = 0.3;
  Double_t maxPt = 20.0;
  if (FittingType.CompareTo("Tsallis")==0){
    fitCombPi0pPb5023GeVPt = FitObject("l","fitInvCrossSectionPi0","Pi0");
    ParameterspPb[0] = 8.66870e+00;//7.4e+1;
    ParameterspPb[1] = 7.03073e+00;
    ParameterspPb[2] = 1.63511e-01;
    fitCombPi0pPb5023GeVPt->SetRange(minPt,maxPt);
    fitCombPi0pPb5023GeVPt->SetParameters(ParameterspPb[0],ParameterspPb[1],ParameterspPb[2]);
    //  fitCombPi0pPb5023GeVPt->SetParLimits(2,0.17,0.18);
  }else if (FittingType.CompareTo("Bylinkin2")==0)  {
    fitCombPi0pPb5023GeVPt = FitObject("2tcm","fitInvCrossSectionPi0","Pi0");
    ParameterspPb[0] = 8.66870e+00;//7.4e+1;
    ParameterspPb[1] = 7.03073e+00;
    ParameterspPb[2] = 1.63511e-01;
    fitCombPi0pPb5023GeVPt->SetRange(minPt,maxPt);
    fitCombPi0pPb5023GeVPt->SetParameters(ParameterspPb[0],ParameterspPb[1],ParameterspPb[2]);

  }else{
    fitCombPi0pPb5023GeVPt = FitObject("tcm","fitInvCrossSectionPi0","Pi0"); 
    ParameterspPb[0] = 0.01;//7.4e+1;
    ParameterspPb[1] = 0.3;
    ParameterspPb[2] = 0.06;
    ParameterspPb[3] = 0.07;
    ParameterspPb[4] = 400.8;
    fitCombPi0pPb5023GeVPt->SetRange(minPt,maxPt);
    fitCombPi0pPb5023GeVPt->SetParameters(ParameterspPb[0],ParameterspPb[1],ParameterspPb[2],ParameterspPb[3],ParameterspPb[4]); // standard

  }    

  graphInvYieldPi0CombpPb5023GeVTotClone->Fit(fitCombPi0pPb5023GeVPt,"QVNRMEX0+","",minPt,maxPt);
  cout<<"!!!!!!!!Fit!!!!!!!!!!!!"<< endl;
  graphInvYieldPi0CombpPb5023GeVTotClone->Fit(fitCombPi0pPb5023GeVPt,"VNRMEX0+","",minPt,maxPt);
  cout<<"!!!!!!!!Fit!!!!!!!!!!!!"<< endl;    
  histo2DCompCombined->DrawCopy();

  fFits<< FittingType.Data()<<endl;
  fFits<<"Fitting Combined Pi0 spectrum"<< endl; 
  if (FittingType.CompareTo("Tsallis")==0){
    for (int i=0;i < 3;i++){    
      fFits<< fitCombPi0pPb5023GeVPt->GetParName(i) <<":  "<< fitCombPi0pPb5023GeVPt->GetParameter(i)<<" +/- "<< fitCombPi0pPb5023GeVPt->GetParError(i)  << endl; 
    }   
  }else {
    for (int i=0;i < 5;i++){    
      fFits<< fitCombPi0pPb5023GeVPt->GetParName(i) <<":  "<< fitCombPi0pPb5023GeVPt->GetParameter(i)<<" +/- "<< fitCombPi0pPb5023GeVPt->GetParError(i) << endl; 
    }  

  }
  fFits<< "Chi2/NDF:  "<<fitCombPi0pPb5023GeVPt->GetChisquare()/fitCombPi0pPb5023GeVPt->GetNDF() << endl;
  //DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0CombpPb5023GeVSysClone,20,1, kBlue, kBlue, 1, kTRUE, kBlue);  
  DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0CombpPb5023GeVSysClone,20,0.7, 4, 4, 1, kTRUE);  
      
  graphInvYieldPi0CombpPb5023GeVSysClone->Draw("E2,same");

  DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0CombpPb5023GeVStaClone,20,0.7, kBlue, kBlue);  
  graphInvYieldPi0CombpPb5023GeVStaClone->Draw("p,same");

     
  TLegend* legendRpPbCombine = new TLegend(0.23,0.30,0.65,0.40);
  legendRpPbCombine->SetFillColor(0);
  legendRpPbCombine->SetLineColor(0);
  legendRpPbCombine->SetTextSize(0.03);
  legendRpPbCombine->AddEntry(graphInvYieldPi0CombpPb5023GeVSysClone,"#pi^{0}, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  if(FittingType.CompareTo("Tsallis")==0) legendRpPbCombine->AddEntry(fitCombPi0pPb5023GeVPt,"Tsallis Fit","l");
  else legendRpPbCombine->AddEntry(fitCombPi0pPb5023GeVPt,"Bylinkin-Rostovtsev Fit","l");
    
  legendRpPbCombine->Draw();
	
    
    
    
  fitCombPi0pPb5023GeVPt->Draw("same");
    
  canvasCompYieldpPbInd->Print(Form("%s/Comb_pPb.%s",outputDir.Data(),suffix.Data()));

  // **************************************************************************************
  // ************************* Creating ratios of different systems ************************
  // **************************************************************************************

  graphRatioCombPCMSys        = CalculateGraphErrRatioToFit(graphRatioCombPCMSys, fitCombPi0pPb5023GeVPt);
  graphRatioCombDalitzSys     = CalculateGraphErrRatioToFit(graphRatioCombDalitzSys, fitCombPi0pPb5023GeVPt);
  graphRatioCombEMCalSys      = CalculateGraphErrRatioToFit(graphRatioCombEMCalSys, fitCombPi0pPb5023GeVPt);
    
  graphRatioCombCombFit       = CalculateGraphErrRatioToFit(graphInvYieldPi0CombpPb5023GeVTotClone,fitCombPi0pPb5023GeVPt); 
  graphRatioCombCombFitSta    = CalculateGraphErrRatioToFit(graphInvYieldPi0CombpPb5023GeVStaClone,fitCombPi0pPb5023GeVPt); 
  graphRatioCombCombFitSys    = CalculateGraphErrRatioToFit(graphInvYieldPi0CombpPb5023GeVSysClone,fitCombPi0pPb5023GeVPt); 
  //graphRatioCombEMCalSys      = CalculateGraphErrRatioToFit(graphRatioCombEMCalSys,fitCombPi0pPb5023GeVPt);

  histoRatioCombPHOS          = CalculateHistoRatioToFit(histoRatioCombPHOS,fitCombPi0pPb5023GeVPt);
    

  TH1D* histoRatioCombPCM     = CalculateHistoRatioToFit(histoPCMYieldPi0pPb,fitCombPi0pPb5023GeVPt);
  TH1D* histoRatioCombDalitz  = CalculateHistoRatioToFit(histoDalitzYieldPi0pPb,fitCombPi0pPb5023GeVPt);
  TH1D* histoRatioCombEMCal   = CalculateHistoRatioToFit(histoEMCalYieldPi0pPbStat,fitCombPi0pPb5023GeVPt);
  TH1D* histoRatioCombPHOSsys = CalculateHistoRatioToFit(histoPHOSYieldPi0pPbSyst,fitCombPi0pPb5023GeVPt);

  // **************************************************************************************
  // ************************* Plotting ratio of different systems ************************
  // **************************************************************************************

  TCanvas* canvasRatioCompYieldpPbInd = new TCanvas("canvasRatioCompYieldpPbInd","",200,10,900,700);  // gives the page size
  DrawGammaCanvasSettings( canvasRatioCompYieldpPbInd,  0.12, 0.02, 0.02, 0.12);
  canvasRatioCompYieldpPbInd->SetLogx();
  //     canvasRatioCompYieldpPbInd->SetGridx();
  //     canvasRatioCompYieldpPbInd->SetGridy();


  DrawPad(0,0.0,0.3,0.12,0.05);
  TH2F * ratio2DInvXSectionPi0;
  ratio2DInvXSectionPi0 = new TH2F("ratio2DInvXSectionPi0","ratio2DInvXSectionPi0",1000,0.3,25.,1000,0.51,1.49);
    
  //SetStyleHistoTH1ForGraphs(ratio2DInvXSectionPi0, "#it{p}_{T} (GeV/#it{c})", "Data/Fit",12,14,12,13,4.85, 1.5, 505,505);
    
  //#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}
    
  SetStyleHistoTH2ForGraphs(ratio2DInvXSectionPi0, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.09,0.11,0.07,0.1,1.,0.45,505,505);

  ratio2DInvXSectionPi0->GetYaxis()->CenterTitle();
	
	    
  /* 
     ratio2DInvXSectionPi0->GetYaxis()->SetLabelSize(0.055) ;
     ratio2DInvXSectionPi0->GetYaxis()->SetTitleSize(0.1) ;
     ratio2DInvXSectionPi0->GetYaxis()->SetTitleOffset(0.3) ;
     ratio2DInvXSectionPi0->SetYTitle("Data/Fit");
     ratio2DInvXSectionPi0->GetXaxis()->SetLabelSize(0.055) ;
     ratio2DInvXSectionPi0->GetXaxis()->SetTitleSize(0.1) ;
     ratio2DInvXSectionPi0->SetXTitle("#it{p}_{T}(GeV/#it{c})");
  */

  ratio2DInvXSectionPi0->DrawCopy();  
  graphRatioCombCombFitSys->SetPointEYlow(30,0.137146);//attribute EMCal Sys error to 
  graphRatioCombCombFitSys->SetPointEYhigh(30,0.137146);
  DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys,20,1, kBlue, kBlue, 1, kTRUE, kBlue-9);  
  graphRatioCombCombFitSys->Draw("E2,same");

  DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFit,20,1, kBlue, kBlue);  
  graphRatioCombCombFit->Draw("pe,same");

  DrawGammaSetMarkerTGraphAsym(graphRatioCombPCMSys,24,1, 1, 1, 1, kTRUE);  
  graphRatioCombPCMSys->Draw("E2,same");

  DrawGammaSetMarker(histoRatioCombPCM, 24,1, 1 , 1);
  histoRatioCombPCM->Draw("same") ;

  TLatex * lt2 = new TLatex(2.,1.25,"PCM") ;
  lt2->SetTextColor(kBlack) ;
  lt2->SetTextSize(0.11) ;
  lt2->Draw() ;

  //=============Dalitz=====================  
  canvasRatioCompYieldpPbInd->cd() ;
  DrawPad(1,0.0,0.,0.12,0.05);
  SetStyleHistoTH2ForGraphs(ratio2DInvXSectionPi0, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.1,0.15,0.1,0.15,0.8,0.3,505,505);
  ratio2DInvXSectionPi0->DrawCopy(); 
  graphRatioCombCombFitSys->Draw("E2,same") ;

  graphRatioCombCombFit->Draw("pe,same");

  DrawGammaSetMarkerTGraphAsym(graphRatioCombDalitzSys,24,1, kCyan+2, kCyan+2, 1, kTRUE);  
  graphRatioCombDalitzSys->Draw("E2same");
  DrawGammaSetMarker(histoRatioCombDalitz, 24,1 ,kCyan+2 ,kCyan+2);
  histoRatioCombDalitz->Draw("same") ; 
  TLatex * lt = new TLatex(2.,0.7,"PCM"); 
  lt->SetTextSize(0.16) ;
  lt->SetTextColor(kCyan+2) ;
  lt->DrawText(2.,1.25,"Dalitz") ;
 
  // ===========EMCal======================  
  canvasRatioCompYieldpPbInd->cd() ;
  DrawPad(2,0.0,0.,0.12,0.05);
  ratio2DInvXSectionPi0->DrawCopy(); 
  
  graphRatioCombCombFitSys->Draw("E2") ;
  
  graphRatioCombCombFit->Draw("pe,same");

            
  DrawGammaSetMarkerTGraphAsym(graphRatioCombEMCalSys,24,1,kGreen+2 , kGreen+2, 1, kTRUE);  
  graphRatioCombEMCalSys->Draw("E2same");
  DrawGammaSetMarker(histoRatioCombEMCal, 24,1 ,kGreen+2 ,kGreen+2);
  histoRatioCombEMCal->Draw("same") ;
  cout<<"EMCal"<<endl;
  graphRatioCombEMCalSys->Print();
  lt->SetTextColor(kGreen+2) ;
  lt->DrawText(2.,1.25,"EMCal") ;

  //PHOS  
  canvasRatioCompYieldpPbInd->cd() ;
  DrawPad(3,0.,0.,0.12,0.05);
  
  ratio2DInvXSectionPi0->DrawCopy(); 
  graphRatioCombCombFitSys->Draw("E2") ;
  graphRatioCombCombFit->Draw("pe,same");
  graphRatioCombCombFit->Print();
  graphRatioCombCombFitSys->Print();
  histoRatioCombPHOSsys->SetLineColor(kRed+1);
  histoRatioCombPHOSsys->SetFillStyle(0) ;
  histoRatioCombPHOSsys->Draw("E2same") ;
  
  histoRatioCombPHOS->SetMarkerColor(kRed+1);
  histoRatioCombPHOS->SetLineColor(kRed+1);
  histoRatioCombPHOS->Draw("same");


  lt->SetTextColor(kRed+1) ;
  lt->DrawText(2.,1.25,"PHOS") ;

  canvasRatioCompYieldpPbInd->Update();
  canvasRatioCompYieldpPbInd->Print(Form("%s/RatioCompYieldpPb.%s",outputDir.Data(),suffix.Data()));
  // **************************************************************************************
  // ************************* Plotting ratio combined to fit ************************
  // **************************************************************************************

  TCanvas* canvasRatioCompYieldpPb = new TCanvas("canvasRatioCompYieldpPb","",200,10,900,700);  // gives the page size
  DrawGammaCanvasSettings( canvasRatioCompYieldpPb,  0.12, 0.02, 0.02, 0.12);
  canvasRatioCompYieldpPb->SetLogx();
  //     canvasRatioCompYieldpPb->SetGridx();
  //     canvasRatioCompYieldpPb->SetGridy();


  TH2F * ratio2DInvXSectionPi0a;
  ratio2DInvXSectionPi0a = new TH2F("ratio2DInvXSectionPi0a","ratio2DInvXSectionPi0a",1000,0.3,25.,1000,0.21,2.69);

  //SetStyleHistoTH1ForGraphs(ratio2DInvXSectionPi0a, "#it{p}_{T} (GeV/#it{c})", "Data/Fit",12,14,12,13,4.85, 1.5, 505,505);
    
  //#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}
    
  if (FittingType.CompareTo("Tsallis")==0)  SetStyleHistoTH2ForGraphs(ratio2DInvXSectionPi0a, "#it{p}_{T} (GeV/#it{c})","Data/Tsallis Fit", 0.05,0.05,0.05,0.05,1.05,1.,505,505);
  else   SetStyleHistoTH2ForGraphs(ratio2DInvXSectionPi0a, "#it{p}_{T} (GeV/#it{c})","Data/Bylinkin-Rostovtsev Fit", 0.05,0.05,0.05,0.05,1.05,1.,505,505);
  //  ratio2DInvXSectionPi0a->GetYaxis()->CenterTitle();
	
	    
  /* 
     ratio2DInvXSectionPi0a->GetYaxis()->SetLabelSize(0.055) ;
     ratio2DInvXSectionPi0a->GetYaxis()->SetTitleSize(0.1) ;
     ratio2DInvXSectionPi0a->GetYaxis()->SetTitleOffset(0.3) ;
     ratio2DInvXSectionPi0a->SetYTitle("Data/Fit");
     ratio2DInvXSectionPi0a->GetXaxis()->SetLabelSize(0.055) ;
     ratio2DInvXSectionPi0a->GetXaxis()->SetTitleSize(0.1) ;
     ratio2DInvXSectionPi0a->SetXTitle("#it{p}_{T}(GeV/#it{c})");
  */

  ratio2DInvXSectionPi0a->DrawCopy();  

  TLine *lineA=new TLine(0.,1.,25.,1.);
  lineA->SetLineColor(kGray+1);
  lineA->Draw();
  DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys,20,1., kBlue, kBlue, 1.5, kTRUE, 0);  
  graphRatioCombCombFitSys->Draw("E2,same");

  DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFit,20,1., kBlue, kBlue,1.5);  
  graphRatioCombCombFit->Draw("p,same");

  TLatex * lt3 = new TLatex(2.8,2.2,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV") ;
  lt3->SetTextColor(kBlack) ;
  lt3->SetTextSize(0.05) ;
  lt3->DrawLatex(2.8,2.0,"#pi^{0}, ALICE");
  //		lt3->DrawLatex(3.,1.8,"#pi^{0} #rightarrow #gamma#gamma");
  //	lt3->DrawLatex(3.,1.8,"#pi^{0}");
  //	lt3->DrawLatex(2.8,1.8,"#pi^{0} #rightarrow #gamma#gamma, #pi^{0} #rightarrow e^{+}e^{-}#gamma");
  lt3->Draw() ;

 

  canvasRatioCompYieldpPb->Update();
  canvasRatioCompYieldpPb->Print(Form("%s/RatioCombYieldFitpPb.%s",outputDir.Data(),suffix.Data()));


  // **************************************************************************************
  // ************************* Plotting ratio individual spectra to fit ************************
  // **************************************************************************************

  TCanvas* canvasRatioIndYieldpPb = new TCanvas("canvasRatioIndYieldpPb","",200,10,900,700);  // gives the page size
  DrawGammaCanvasSettings( canvasRatioIndYieldpPb,  0.12, 0.02, 0.02, 0.12);
  canvasRatioIndYieldpPb->SetLogx();
  //     canvasRatioIndYieldpPb->SetGridx();
  //     canvasRatioIndYieldpPb->SetGridy();


  TH2F * ratio2DInvXSectionPi0b;
  ratio2DInvXSectionPi0b = new TH2F("ratio2DInvXSectionPi0b","ratio2DInvXSectionPi0b",1000,0.3,25.,1000,0.21,2.69);

  //SetStyleHistoTH1ForGraphs(ratio2DInvXSectionPi0b, "#it{p}_{T} (GeV/#it{c})", "Data/Fit",12,14,12,13,4.85, 1.5, 505,505);
    
  //#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}
    
  if (FittingType.CompareTo("Tsallis")==0)  SetStyleHistoTH2ForGraphs(ratio2DInvXSectionPi0b, "#it{p}_{T} (GeV/#it{c})","Data/Tsallis Fit", 0.05,0.05,0.05,0.05,1.05,1.,505,505);
  else   SetStyleHistoTH2ForGraphs(ratio2DInvXSectionPi0b, "#it{p}_{T} (GeV/#it{c})","Data/Bylinkin-Rostovtsev Fit", 0.05,0.05,0.05,0.05,1.05,1.,505,505);
  //  ratio2DInvXSectionPi0b->GetYaxis()->CenterTitle();
	
	    
  /* 
     ratio2DInvXSectionPi0b->GetYaxis()->SetLabelSize(0.055) ;
     ratio2DInvXSectionPi0b->GetYaxis()->SetTitleSize(0.1) ;
     ratio2DInvXSectionPi0b->GetYaxis()->SetTitleOffset(0.3) ;
     ratio2DInvXSectionPi0b->SetYTitle("Data/Fit");
     ratio2DInvXSectionPi0b->GetXaxis()->SetLabelSize(0.055) ;
     ratio2DInvXSectionPi0b->GetXaxis()->SetTitleSize(0.1) ;
     ratio2DInvXSectionPi0b->SetXTitle("#it{p}_{T}(GeV/#it{c})");
  */

  ratio2DInvXSectionPi0b->DrawCopy();  

  TLine *lineB=new TLine(0.,1.,25.,1.);
  lineB->SetLineColor(kGray+1);
  lineB->Draw();


  DrawGammaSetMarkerTGraphAsym(graphRatioCombPCMSys,20,1, 1, 1, 1, kTRUE);  
  graphRatioCombPCMSys->Draw("E2same");
  DrawGammaSetMarker(histoRatioCombPCM, 20,1, 1 , 1);
  histoRatioCombPCM->Draw("same") ;

  DrawGammaSetMarkerTGraphAsym(graphRatioCombDalitzSys,29,1, kCyan+2, kCyan+2, 1, kTRUE);  
  graphRatioCombDalitzSys->Draw("E2same");
  DrawGammaSetMarker(histoRatioCombDalitz, 29,1.5 ,kCyan+2 ,kCyan+2);
  histoRatioCombDalitz->Draw("same") ; 
 
  DrawGammaSetMarkerTGraphAsym(graphRatioCombEMCalSys,33,1,kGreen+2 , kGreen+2, 1, kTRUE);  
  graphRatioCombEMCalSys->Draw("E2same");
  DrawGammaSetMarker(histoRatioCombEMCal, 33,1.5 ,kGreen+2 ,kGreen+2);
  histoRatioCombEMCal->Draw("same") ;
	
  histoRatioCombPHOSsys->Draw("E2same") ;
	
  histoRatioCombPHOS->SetMarkerStyle(21);
  histoRatioCombPHOS->SetMarkerSize(1);
  histoRatioCombPHOS->Draw("same");
	
  TLatex * lt4 = new TLatex(2.8,2.2,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV") ;
  lt4->SetTextColor(kBlack) ;
  lt4->SetTextSize(0.05) ;
  lt4->DrawLatex(2.8,2.,"#pi^{0}, ALICE");
  //	lt4->DrawLatex(3.,1.8,"#pi^{0} #rightarrow #gamma#gamma");
  //	lt4->DrawLatex(2.8,1.8,"#pi^{0} #rightarrow #gamma#gamma, #pi^{0} #rightarrow e^{+}e^{-}#gamma");
  lt4->Draw() ;

  TLegend* legendSpectraDiffDetMinBiasStrip = new TLegend(0.18,0.7,0.5,0.95);
  legendSpectraDiffDetMinBiasStrip->SetFillColor(0);
  legendSpectraDiffDetMinBiasStrip->SetLineColor(0);
  legendSpectraDiffDetMinBiasStrip->SetTextFont(42);
  legendSpectraDiffDetMinBiasStrip->AddEntry(histoRatioCombPCM,Form("PCM"),"pf");
  legendSpectraDiffDetMinBiasStrip->AddEntry(histoRatioCombDalitz,Form("Dalitz"),"pf");
  legendSpectraDiffDetMinBiasStrip->AddEntry(histoRatioCombEMCal,Form("EMCal"),"pf");
  legendSpectraDiffDetMinBiasStrip->AddEntry(histoRatioCombPHOS,Form("PHOS"),"pf");
  legendSpectraDiffDetMinBiasStrip->Draw();

  canvasRatioIndYieldpPb->Update();
  canvasRatioIndYieldpPb->Print(Form("%s/RatioIndYieldFitpPb.%s",outputDir.Data(),suffix.Data()));

  //------------------ Apply x Shift to Combined Spectrum------------------
  TF1* fitCombPi0pPb5023GeVPtXShift=(TF1*)fitCombPi0pPb5023GeVPt->Clone("FitXShift");
  cout << "X shift!!!!"<<endl;
  TGraphAsymmErrors* graphCombPi0InvCrossSectionStatpPb5023GeVXShifted  = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionStatpPb5023GeV->Clone(); 
  cout << "X shift!!!!"<<endl;	
  TGraphAsymmErrors* graphCombPi0InvCrossSectionSyspPb5023GeVXShifted  = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionSyspPb5023GeV->Clone();
  TGraphAsymmErrors* graphCombPi0InvCrossSectionTotpPb5023GeVXShifted = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionTotpPb5023GeV->Clone();
  graphCombPi0InvCrossSectionTotpPb5023GeVXShifted = ApplyXshift(graphCombPi0InvCrossSectionTotpPb5023GeVXShifted ,fitCombPi0pPb5023GeVPtXShift,"Pi0");
  graphCombPi0InvCrossSectionSyspPb5023GeVXShifted = ApplyXshiftIndividualSpectra(graphCombPi0InvCrossSectionTotpPb5023GeVXShifted,graphCombPi0InvCrossSectionSyspPb5023GeVXShifted , fitCombPi0pPb5023GeVPtXShift, 0, graphCombPi0InvCrossSectionTotpPb5023GeVXShifted->GetN()); 
  graphCombPi0InvCrossSectionStatpPb5023GeVXShifted  = ApplyXshiftIndividualSpectra(graphCombPi0InvCrossSectionTotpPb5023GeVXShifted,graphCombPi0InvCrossSectionStatpPb5023GeVXShifted , fitCombPi0pPb5023GeVPtXShift, 0, graphCombPi0InvCrossSectionTotpPb5023GeVXShifted->GetN());

  TCanvas* canvasXShift = new TCanvas("canvasXShift","",200,10,500,400);  // gives the page size
  DrawGammaCanvasSettings( canvasXShift, 0.15, 0.02, 0.02, 0.12);
    
  canvasXShift->SetLogx();
  canvasXShift->SetLogy();
  TH2F * histo2DCompCombinedXShift;
  histo2DCompCombinedXShift = new TH2F("histo2DCompCombinedXShift","histo2DCompCombinedXShift",1000,0.3,30.,1000,1.2e-9,30.   );
  SetStyleHistoTH2ForGraphs(histo2DCompCombinedXShift, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.03,0.04, 0.03,0.04, 1.0,1.4, 512, 508);



  histo2DCompCombinedXShift->DrawCopy();  
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionTotpPb5023GeV,21,0.7, 4, 4, 1, kTRUE);  
  graphCombPi0InvCrossSectionTotpPb5023GeV->Draw("pe,same");

  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionTotpPb5023GeVXShifted,20,0.7, 2, 2, 1, kTRUE);  
  graphCombPi0InvCrossSectionTotpPb5023GeVXShifted->Draw("pe,same");
  // DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionSyspPb5023GeVXShifted,20,0.7, 4, 4, 1, kTRUE); 
     
      
  // graphCombPi0InvCrossSectionSyspPb5023GeVXShifted->Draw("E2,same");

  // DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionStatpPb5023GeVXShifted,20,0.7, kBlue, kBlue);  
  // graphCombPi0InvCrossSectionStatpPb5023GeVXShifted->Draw("p,same");
  fitCombPi0pPb5023GeVPt->Draw("same");
     
  TLegend* legendRpPbCombineXShift = new TLegend(0.23,0.25,0.65,0.40);
  legendRpPbCombineXShift->SetFillColor(0);
  legendRpPbCombineXShift->SetLineColor(0);
  legendRpPbCombineXShift->SetTextSize(0.03);
  legendRpPbCombineXShift->AddEntry(graphCombPi0InvCrossSectionTotpPb5023GeV,"#pi^{0} , p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  legendRpPbCombineXShift->AddEntry(graphCombPi0InvCrossSectionTotpPb5023GeVXShifted,"#pi^{0} x bin shifted, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  if (FittingType.CompareTo("Tsallis")==0)  legendRpPbCombineXShift->AddEntry(fitCombPi0pPb5023GeVPt,"Tsallis","l");
  else  legendRpPbCombineXShift->AddEntry(fitCombPi0pPb5023GeVPt,"Bylinkin-Rostovtsev Fit","l");
    
  legendRpPbCombineXShift->Draw();
	


  canvasXShift->Update();
  canvasXShift->Print(Form("%s/CombinedSpectrum_XShifted.%s",outputDir.Data(),suffix.Data()));



  Double_t* XPi0Unshifted=graphCombPi0InvCrossSectionTotpPb5023GeV->GetX();
  Double_t* XPi0Shifted=graphCombPi0InvCrossSectionTotpPb5023GeVXShifted->GetX();

  for (int i=0;i<graphCombPi0InvCrossSectionTotpPb5023GeV->GetN();i++){

    fFits <<  XPi0Unshifted[i]-graphCombPi0InvCrossSectionTotpPb5023GeV->GetErrorXlow(i)<< " - " <<  XPi0Unshifted[i]+ graphCombPi0InvCrossSectionTotpPb5023GeV->GetErrorXhigh(i)<< " GeV/$c$ & " << XPi0Unshifted[i] <<" GeV/$c$ & " << XPi0Shifted[i]<<" GeV/$c$ & " << (XPi0Shifted[i]-XPi0Unshifted[i])*1000 <<" MeV/$c$" << endl;
  }








  // **************************************************************************************
  // ********************************* Combine Eta spectra ********************************
  // **************************************************************************************
    
  TString fileNameOutputWeightingEta                 = Form("%s/WeightingEta.dat",outputDir.Data());

  statErrorCollectionEta[0]          = (TH1D*)histoPCMYieldEtapPb->Clone("statErrPCMEta");
  statErrorCollectionEta[2]          = (TH1D*)histoEMCalYieldEtapPbStat->Clone("statErrEMCalEta");
    
  sysErrorCollectionEta[0]           = (TGraphAsymmErrors*)graphPCMYieldEtapPbSystErr->Clone("sysErrPCMEta");
  sysErrorCollectionEta[2]           = (TGraphAsymmErrors*)graphEMCalYieldEtapPbSystErr->Clone("sysErrEMCalEta");
    
    
  TGraphAsymmErrors* graphCombEtaInvCrossSectionStatpPb5023GeV= NULL;
  TGraphAsymmErrors* graphCombEtaInvCrossSectionSyspPb5023GeV = NULL;
    
  TGraphAsymmErrors* graphCombEtaInvCrossSectionTotpPb5023GeV = CombinePtPointsSpectraFullCorrMat(    statErrorCollectionEta,    sysErrorCollectionEta,     
												      pTLimitsEta, NtotalEta,
												      offSetsEta, offSetsSysEta,
												      graphCombEtaInvCrossSectionStatpPb5023GeV, graphCombEtaInvCrossSectionSyspPb5023GeV,
												      fileNameOutputWeightingEta,1
												      );
  graphCombEtaInvCrossSectionStatpPb5023GeV->Print();
    
  TGraphAsymmErrors* graphInvYieldEtaCombpPb5023GeVStaClone   = (TGraphAsymmErrors*) graphCombEtaInvCrossSectionStatpPb5023GeV->Clone();
  TGraphAsymmErrors* graphInvYieldEtaCombpPb5023GeVSysClone   = (TGraphAsymmErrors*) graphCombEtaInvCrossSectionSyspPb5023GeV->Clone();
  TGraphAsymmErrors* graphInvYieldEtaCombpPb5023GeVTotClone   = (TGraphAsymmErrors*) graphCombEtaInvCrossSectionTotpPb5023GeV->Clone();
    
  TGraphAsymmErrors* graphRatioCombEtaFit                    = (TGraphAsymmErrors*) graphInvYieldEtaCombpPb5023GeVTotClone ->Clone(); 
  TGraphAsymmErrors* graphRatioCombEtaFitSta                 = (TGraphAsymmErrors*) graphInvYieldEtaCombpPb5023GeVStaClone ->Clone();
  TGraphAsymmErrors* graphRatioCombEtaFitSys                 = (TGraphAsymmErrors*) graphInvYieldEtaCombpPb5023GeVSysClone ->Clone(); 

    

  //test Tsallis//
	
  TF1* fitTsallisEtapPb5023GeVPt = FitObject("tmpt","tmptEtapPb5023GeVPt","Eta");
	

	
	

   
  //      fitTsallisEtapPb5023GeVPt->SetParameters();
  fitTsallisEtapPb5023GeVPt->SetRange(0.7,20.);
  //   //test Tsallis//


  // // 1  p0           1.18085e+00   4.01565e-01   1.04073e-04   3.99481e-04
  // // 2  p1           8.10783e+00   9.35534e-01   1.86868e-04  -7.95453e-04
  // // 3  p2           1.87850e-01   3.37139e-02   4.78260e-06   2.60935e-03
  // // 4  p3          -4.62234e-09   3.08040e-08   1.67148e-11  -8.25946e+02



  // **************************************************************************************
  // ************************* Plotting Eta yield of different systems ************************
  // **************************************************************************************

  TCanvas* canvasCompEtaYieldpPbInd = new TCanvas("canvasCompEtaYieldpPbInd","",200,10,500,400);  // gives the page size
  DrawGammaCanvasSettings( canvasCompEtaYieldpPbInd, 0.15, 0.02, 0.02, 0.12);
    
  canvasCompEtaYieldpPbInd->SetLogx();
  canvasCompEtaYieldpPbInd->SetLogy();
  TH2F * histo2DCompEtaCombined;
  histo2DCompEtaCombined = new TH2F("histo2DCompEtaCombined","histo2DCompEtaCombined",1000,0.3,30.,1000,1.2e-9,3.   );
  SetStyleHistoTH2ForGraphs(histo2DCompEtaCombined, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.03,0.04, 0.03,0.04, 1.0,1.4, 512, 508);
  // SetStyleHistoTH2ForGraphs(histo2DCompEtaCombined, "#it{p}_{T} (GeV/#it{c})","arbitrary units", 0.03,0.04, 0.03,0.04, 1.0,1.1, 512, 508);
  //   histo2DCompEtaCombined->GetXaxis()->SetRangeUser(0.,30.);
  //  histo2DCompEtaCombined->GetYaxis()->SetLabelColor(0);
  //  histo2DCompEtaCombined->GetYaxis()->SetRangeUser(0.1,2.1);
  histo2DCompEtaCombined->DrawCopy();
    
  //     TF1 * CurrentFit = new TF1();
  TF1* fitCombEtapPb5023GeVPt;
  Double_t minPtEta = 0.7;
  Double_t maxPtEta = 20.0;
  Double_t *ParameterspPbEta  = new Double_t[5];// = { 7.4e+10, 0.3, 1e+09,0.3,8};
  if (FittingType.CompareTo("Tsallis")==0){
    fitCombEtapPb5023GeVPt= FitObject("l","fitInvCrossSectionEta","Eta");
    fitCombEtapPb5023GeVPt->SetRange(minPtEta,maxPtEta);
  } else if (FittingType.CompareTo("Bylinkin2")==0) {
    fitCombEtapPb5023GeVPt= FitObject("2tcm","fitInvCrossSectionEta","Eta"); 
    ParameterspPbEta[0] = 0.01;//7.4e+1;
    ParameterspPbEta[1] = 0.3;
    ParameterspPbEta[2] = .06;
    ParameterspPbEta[3] = 0.07;
    ParameterspPbEta[4] = 400.8;
    fitCombEtapPb5023GeVPt->SetRange(minPtEta,maxPtEta);
    fitCombEtapPb5023GeVPt->SetParameters(ParameterspPbEta[0],ParameterspPbEta[1],ParameterspPbEta[2],ParameterspPbEta[3],ParameterspPbEta[4]); // standard
  }else{
    fitCombEtapPb5023GeVPt= FitObject("tcm","fitInvCrossSectionEta","Eta"); 
    ParameterspPbEta[0] = 0.01;//7.4e+1;
    ParameterspPbEta[1] = 0.3;
    ParameterspPbEta[2] = .06;
    ParameterspPbEta[3] = 0.07;
    ParameterspPbEta[4] = 400.8;
    fitCombEtapPb5023GeVPt->SetRange(minPtEta,maxPtEta);
    fitCombEtapPb5023GeVPt->SetParameters(ParameterspPbEta[0],ParameterspPbEta[1],ParameterspPbEta[2],ParameterspPbEta[3],ParameterspPbEta[4]); // standard
  }    
 
  graphInvYieldEtaCombpPb5023GeVTotClone->Fit(fitCombEtapPb5023GeVPt,"QVNRMEX0+","",minPtEta,maxPtEta);
  cout<<"!!!!!!!!Fit!!!!!!!!!!!!"<< endl;
  graphInvYieldEtaCombpPb5023GeVTotClone->Fit(fitCombEtapPb5023GeVPt,"VNRMEX0+","",minPtEta,maxPtEta);
  cout<<"!!!!!!!!Fit!!!!!!!!!!!!"<< endl; 

  
  fFits<<"Fitting Combined Eta spectrum"<< endl;  
  if (FittingType.CompareTo("Tsallis")==0){
    for (int i=0;i < 3;i++){    
      fFits<< fitCombEtapPb5023GeVPt->GetParName(i) <<":  "<< fitCombEtapPb5023GeVPt->GetParameter(i)<<" +/- "<< fitCombEtapPb5023GeVPt->GetParError(i)  << endl; 
    }   
  }else {
    for (int i=0;i < 5;i++){    
      fFits<< fitCombEtapPb5023GeVPt->GetParName(i) <<":  "<< fitCombEtapPb5023GeVPt->GetParameter(i)<<" +/- "<< fitCombEtapPb5023GeVPt->GetParError(i) << endl; 
    }  

  } 
  fFits<< "Chi2/NDF:  "<<fitCombEtapPb5023GeVPt->GetChisquare()/fitCombEtapPb5023GeVPt->GetNDF() << endl;
   
  // graphInvYieldEtaCombpPb5023GeVTotClone->Fit(fitTsallisEtapPb5023GeVPt,"QVNRMEX0+","",minPtEta,maxPtEta);
  //  cout<<"!!!!!!!!Fit!!!!!!!!!!!!"<< endl;
  //  graphInvYieldEtaCombpPb5023GeVTotClone->Fit(fitTsallisEtapPb5023GeVPt,"VNRMEX0+","",minPtEta,maxPtEta);
  //  cout<<"!!!!!!!!Fit!!!!!!!!!!!!"<< endl;    
  //  histo2DCompEtaCombined->DrawCopy();
    
     
    
  //DrawGammaSetMarkerTGraphAsym(graphInvYieldEtaCombpPb5023GeVSysClone,20,1, kBlue, kBlue, 1, kTRUE, kBlue);  
  DrawGammaSetMarkerTGraphAsym(graphInvYieldEtaCombpPb5023GeVSysClone,20,0.7, 4, 4, 1, kTRUE);  
      
  graphInvYieldEtaCombpPb5023GeVSysClone->Draw("E2,same");

  DrawGammaSetMarkerTGraphAsym(graphInvYieldEtaCombpPb5023GeVStaClone,20,0.7, kBlue, kBlue);  
  graphInvYieldEtaCombpPb5023GeVStaClone->Draw("p,same");

     
  TLegend* legendRpPbCombineEta = new TLegend(0.23,0.30,0.65,0.40);
  legendRpPbCombineEta->SetFillColor(0);
  legendRpPbCombineEta->SetLineColor(0);
  legendRpPbCombineEta->SetTextSize(0.03);
  legendRpPbCombineEta->AddEntry(graphInvYieldEtaCombpPb5023GeVSysClone,"#eta, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  if (FittingType.CompareTo("Tsallis")==0)  legendRpPbCombineEta->AddEntry(fitCombEtapPb5023GeVPt,"Tsallis Fit","l");
  else legendRpPbCombineEta->AddEntry(fitCombEtapPb5023GeVPt,"Bylinkin-Rostovtsev Fit","l");
   
    
  legendRpPbCombineEta->Draw();
	
    
    
    
  fitCombEtapPb5023GeVPt->Draw("same");
    
  canvasCompEtaYieldpPbInd->Print(Form("%s/CombEta_pPb.%s",outputDir.Data(),suffix.Data()));

  // **************************************************************************************
  // ************************* Creating ratios of different systems ************************
  // **************************************************************************************

  graphRatioCombEtaPCMSys        = CalculateGraphErrRatioToFit(graphRatioCombEtaPCMSys, fitCombEtapPb5023GeVPt);
  graphRatioCombEtaEMCalSys      = CalculateGraphErrRatioToFit(graphRatioCombEtaEMCalSys, fitCombEtapPb5023GeVPt);
    
  graphRatioCombEtaFit       = CalculateGraphErrRatioToFit(graphInvYieldEtaCombpPb5023GeVTotClone,fitCombEtapPb5023GeVPt); 
  graphRatioCombEtaFitSta    = CalculateGraphErrRatioToFit(graphInvYieldEtaCombpPb5023GeVStaClone,fitCombEtapPb5023GeVPt); 
  graphRatioCombEtaFitSys    = CalculateGraphErrRatioToFit(graphInvYieldEtaCombpPb5023GeVSysClone,fitCombEtapPb5023GeVPt); 
  //graphRatioCombEtaEMCalSys      = CalculateGraphErrRatioToFit(graphRatioCombEtaEMCalSys,fitCombEtapPb5023GeVPt);

  TH1D* histoRatioCombEtaPCM     = CalculateHistoRatioToFit(histoPCMYieldEtapPb,fitCombEtapPb5023GeVPt);
  TH1D* histoRatioCombEtaEMCal   = CalculateHistoRatioToFit(histoEMCalYieldEtapPbStat,fitCombEtapPb5023GeVPt);
 
  // graphRatioCombEtaPCMSys        = CalculateGraphErrRatioToFit(graphRatioCombEtaPCMSys,fitTsallisEtapPb5023GeVPt );
  // graphRatioCombEtaEMCalSys      = CalculateGraphErrRatioToFit(graphRatioCombEtaEMCalSys,fitTsallisEtapPb5023GeVPt );
    
  // graphRatioCombEtaFit       = CalculateGraphErrRatioToFit(graphInvYieldEtaCombpPb5023GeVTotClone,fitTsallisEtapPb5023GeVPt); 
  // graphRatioCombEtaFitSta    = CalculateGraphErrRatioToFit(graphInvYieldEtaCombpPb5023GeVStaClone,fitTsallisEtapPb5023GeVPt); 
  // graphRatioCombEtaFitSys    = CalculateGraphErrRatioToFit(graphInvYieldEtaCombpPb5023GeVSysClone,fitTsallisEtapPb5023GeVPt); 
  // //graphRatioCombEtaEMCalSys      = CalculateGraphErrRatioToFit(graphRatioCombEtaEMCalSys,fitTsallisEtapPb5023GeVPt);

  //  TH1D* histoRatioCombEtaPCM     = CalculateHistoRatioToFit(histoPCMYieldEtapPb,fitTsallisEtapPb5023GeVPt);
  //  TH1D* histoRatioCombEtaEMCal   = CalculateHistoRatioToFit(histoEMCalYieldEtapPbStat,fitTsallisEtapPb5023GeVPt);
  // **************************************************************************************
  // ************************* Plotting ratio of different systems ************************
  // **************************************************************************************

  TCanvas* canvasRatioCompEtaYieldpPbInd = new TCanvas("canvasRatioCompEtaYieldpPbInd","",200,10,900,700);  // gives the page size
  DrawGammaCanvasSettings( canvasRatioCompEtaYieldpPbInd,  0.12, 0.02, 0.02, 0.12);
  canvasRatioCompEtaYieldpPbInd->SetLogx();
  //     canvasRatioCompEtaYieldpPbInd->SetGridx();
  //     canvasRatioCompEtaYieldpPbInd->SetGridy();


  Draw2Pads(0,0.,0.2,0.1,0.03);
  TH2F * ratio2DInvXSectionEta;
  ratio2DInvXSectionEta = new TH2F("ratio2DInvXSectionEta","ratio2DInvXSectionEta",1000,0.6,25.,1000,0.31,1.69);
    
  //SetStyleHistoTH1ForGraphs(ratio2DInvXSectionEta, "#it{p}_{T} (GeV/#it{c})", "Data/Fit",12,14,12,13,4.85, 1.5, 505,505);
    
  //#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}
    
  SetStyleHistoTH2ForGraphs(ratio2DInvXSectionEta, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.07,0.09,0.055,0.08,0.9,0.55,505,505);
  ratio2DInvXSectionEta->GetYaxis()->CenterTitle();
	
	    
  /* 
     ratio2DInvXSectionEta->GetYaxis()->SetLabelSize(0.055) ;
     ratio2DInvXSectionEta->GetYaxis()->SetTitleSize(0.1) ;
     ratio2DInvXSectionEta->GetYaxis()->SetTitleOffset(0.3) ;
     ratio2DInvXSectionEta->SetYTitle("Data/Fit");
     ratio2DInvXSectionEta->GetXaxis()->SetLabelSize(0.055) ;
     ratio2DInvXSectionEta->GetXaxis()->SetTitleSize(0.1) ;
     ratio2DInvXSectionEta->SetXTitle("#it{p}_{T}(GeV/#it{c})");
  */

  ratio2DInvXSectionEta->DrawCopy();  
  graphRatioCombEtaFitSys->SetPointEYlow(30,0.137146);//attribute EMCal Sys error to 
  graphRatioCombEtaFitSys->SetPointEYhigh(30,0.137146);
  DrawGammaSetMarkerTGraphAsym(graphRatioCombEtaFitSys,20,1, kBlue, kBlue, 1, kTRUE, kBlue-9);  
  graphRatioCombEtaFitSys->Draw("E2,same");

  DrawGammaSetMarkerTGraphAsym(graphRatioCombEtaFit,20,1, kBlue, kBlue);  
  graphRatioCombEtaFit->Draw("pe,same");

  DrawGammaSetMarkerTGraphAsym(graphRatioCombEtaPCMSys,24,1, 1, 1, 1, kTRUE);  
  graphRatioCombEtaPCMSys->Draw("E2,same");

  DrawGammaSetMarker(histoRatioCombEtaPCM, 24,1, 1 , 1);
  histoRatioCombEtaPCM->Draw("same") ;

  TLatex * lt2a = new TLatex(1.,1.3,"PCM") ;
  lt2a->SetTextColor(kBlack) ;
  lt2a->SetTextSize(0.08) ;
  lt2a->Draw() ;


 
  // ===========EMCal======================     
  canvasRatioCompEtaYieldpPbInd->cd() ;
  //      Draw2Pads(1,0.05,0.,0.1,0.03);
  Draw2Pads(1,0.0,0.,0.1,0.03);

  SetStyleHistoTH2ForGraphs(ratio2DInvXSectionEta, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.09,0.11,0.07,0.1,1.,0.45,505,505);
  ratio2DInvXSectionEta->DrawCopy(); 
  
  graphRatioCombEtaFitSys->Draw("E2") ;
  
  graphRatioCombEtaFit->Draw("pe,same");

            
  DrawGammaSetMarkerTGraphAsym(graphRatioCombEtaEMCalSys,24,1,kGreen+2 , kGreen+2, 1, kTRUE);  
  graphRatioCombEtaEMCalSys->Draw("E2same");
  DrawGammaSetMarker(histoRatioCombEtaEMCal, 24,1 ,kGreen+2 ,kGreen+2);
  histoRatioCombEtaEMCal->Draw("same") ;
  cout<<"EMCal"<<endl;
  graphRatioCombEtaEMCalSys->Print();
  lt->SetTextColor(kGreen+2) ;
  lt->SetTextSize(0.1) ;
  lt->DrawText(1.,1.3,"EMCal") ;
  TLatex* labelEta=new TLatex(10.5,1.55,"#eta #rightarrow #gamma#gamma") ;
  labelEta->SetTextSize(0.1);
  labelEta->Draw("same") ;
  canvasRatioCompEtaYieldpPbInd->Update();
  canvasRatioCompEtaYieldpPbInd->Print(Form("%s/RatioCompEtaYieldpPb.%s",outputDir.Data(),suffix.Data()));
  // **************************************************************************************
  // ************************* Plotting Eta ratio combined to fit ************************
  // **************************************************************************************

  TCanvas* canvasRatioCompEtaYieldpPb = new TCanvas("canvasRatioCompEtaYieldpPb","",200,10,900,700);  // gives the page size
  DrawGammaCanvasSettings( canvasRatioCompEtaYieldpPb,  0.12, 0.02, 0.02, 0.12);
  canvasRatioCompEtaYieldpPb->SetLogx();
  //     canvasRatioCompEtaYieldpPb->SetGridx();
  //     canvasRatioCompEtaYieldpPb->SetGridy();


  TH2F * ratio2DInvXSectionEtaa;
  ratio2DInvXSectionEtaa = new TH2F("ratio2DInvXSectionEtaa","ratio2DInvXSectionEtaa",1000,0.5,25.,1000,0.21,2.69);

  //SetStyleHistoTH1ForGraphs(ratio2DInvXSectionEtaa, "#it{p}_{T} (GeV/#it{c})", "Data/Fit",12,14,12,13,4.85, 1.5, 505,505);
    
  //#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}
    
  if (FittingType.CompareTo("Tsallis")==0)      SetStyleHistoTH2ForGraphs(ratio2DInvXSectionEtaa, "#it{p}_{T} (GeV/#it{c})","Data/Tsallis Fit", 0.05,0.05,0.05,0.05,1.05,1.,505,505);
  else       SetStyleHistoTH2ForGraphs(ratio2DInvXSectionEtaa, "#it{p}_{T} (GeV/#it{c})","Data/Bylinkin-Rostovtsev Fit", 0.05,0.05,0.05,0.05,1.05,1.,505,505);
  //  ratio2DInvXSectionEtaa->GetYaxis()->CenterTitle();
	
	    
  /* 
     ratio2DInvXSectionEtaa->GetYaxis()->SetLabelSize(0.055) ;
     ratio2DInvXSectionEtaa->GetYaxis()->SetTitleSize(0.1) ;
     ratio2DInvXSectionEtaa->GetYaxis()->SetTitleOffset(0.3) ;
     ratio2DInvXSectionEtaa->SetYTitle("Data/Fit");
     ratio2DInvXSectionEtaa->GetXaxis()->SetLabelSize(0.055) ;
     ratio2DInvXSectionEtaa->GetXaxis()->SetTitleSize(0.1) ;
     ratio2DInvXSectionEtaa->SetXTitle("#it{p}_{T}(GeV/#it{c})");
  */

  ratio2DInvXSectionEtaa->DrawCopy();  

  TLine *lineAEta=new TLine(0.,1.,25.,1.);
  lineAEta->SetLineColor(kGray+1);
  lineAEta->Draw();
  DrawGammaSetMarkerTGraphAsym(graphRatioCombEtaFitSys,20,1., kBlue, kBlue, 1.5, kTRUE, 0);  
  graphRatioCombEtaFitSys->Draw("E2,same");

  DrawGammaSetMarkerTGraphAsym(graphRatioCombEtaFit,20,1., kBlue, kBlue,1.5);  
  graphRatioCombEtaFit->Draw("p,same");

  TLatex * lt3a = new TLatex(2.8,2.2,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV") ;
  lt3a->SetTextColor(kBlack) ;
  lt3a->SetTextSize(0.05) ;
  lt3a->DrawLatex(2.8,2.0,"#eta, ALICE");
  //		lt3->DrawLatex(3.,1.8,"#pi^{0} #rightarrow #gamma#gamma");
  //	lt3->DrawLatex(3.,1.8,"#pi^{0}");
  //	lt3->DrawLatex(2.8,1.8,"#pi^{0} #rightarrow #gamma#gamma, #pi^{0} #rightarrow e^{+}e^{-}#gamma");
  lt3->Draw() ;

 

  canvasRatioCompEtaYieldpPb->Update();
  canvasRatioCompEtaYieldpPb->Print(Form("%s/RatioCombEtaYieldFitpPb.%s",outputDir.Data(),suffix.Data()));


  // **************************************************************************************
  // ************************* Plotting ratio individual Eta spectra to fit ************************
  // **************************************************************************************

  TCanvas* canvasRatioIndEtaYieldpPb = new TCanvas("canvasRatioIndEtaYieldpPb","",200,10,900,700);  // gives the page size
  DrawGammaCanvasSettings( canvasRatioIndEtaYieldpPb,  0.12, 0.02, 0.02, 0.12);
  canvasRatioIndEtaYieldpPb->SetLogx();
  //     canvasRatioIndEtaYieldpPb->SetGridx();
  //     canvasRatioIndEtaYieldpPb->SetGridy();


  TH2F * ratio2DInvXSectionEtab;
  ratio2DInvXSectionEtab = new TH2F("ratio2DInvXSectionEtab","ratio2DInvXSectionEtab",1000,0.5,25.,1000,0.41,2.69);

  //SetStyleHistoTH1ForGraphs(ratio2DInvXSectionEtab, "#it{p}_{T} (GeV/#it{c})", "Data/Fit",12,14,12,13,4.85, 1.5, 505,505);
    
  //#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}
    
  if (FittingType.CompareTo("Tsallis")==0)      SetStyleHistoTH2ForGraphs(ratio2DInvXSectionEtab, "#it{p}_{T} (GeV/#it{c})","Data/Tsallis Fit", 0.05,0.05,0.05,0.05,1.05,1.,505,505);
  else    SetStyleHistoTH2ForGraphs(ratio2DInvXSectionEtab, "#it{p}_{T} (GeV/#it{c})","Data/Bylinkin-Rostovtsev Fit", 0.05,0.05,0.05,0.05,1.05,1.,505,505);
  //  ratio2DInvXSectionEtab->GetYaxis()->CenterTitle();
	
	    
  /* 
     ratio2DInvXSectionEtab->GetYaxis()->SetLabelSize(0.055) ;
     ratio2DInvXSectionEtab->GetYaxis()->SetTitleSize(0.1) ;
     ratio2DInvXSectionEtab->GetYaxis()->SetTitleOffset(0.3) ;
     ratio2DInvXSectionEtab->SetYTitle("Data/Fit");
     ratio2DInvXSectionEtab->GetXaxis()->SetLabelSize(0.055) ;
     ratio2DInvXSectionEtab->GetXaxis()->SetTitleSize(0.1) ;
     ratio2DInvXSectionEtab->SetXTitle("#it{p}_{T}(GeV/#it{c})");
  */

  ratio2DInvXSectionEtab->DrawCopy();  

  TLine *lineBEta=new TLine(0.,1.,25.,1.);
  lineBEta->SetLineColor(kGray+1);
  lineBEta->Draw();


  DrawGammaSetMarkerTGraphAsym(graphRatioCombEtaPCMSys,20,1, 1, 1, 1, kTRUE);  
  graphRatioCombEtaPCMSys->Draw("E2same");
  DrawGammaSetMarker(histoRatioCombEtaPCM, 20,1, 1 , 1);
  histoRatioCombEtaPCM->Draw("same") ;
 
  DrawGammaSetMarkerTGraphAsym(graphRatioCombEtaEMCalSys,33,1,kGreen+2 , kGreen+2, 1, kTRUE);  
  graphRatioCombEtaEMCalSys->Draw("E2same");
  DrawGammaSetMarker(histoRatioCombEtaEMCal, 33,1.5 ,kGreen+2 ,kGreen+2);
  histoRatioCombEtaEMCal->Draw("same") ;
	
  TLatex * lt4a = new TLatex(2.8,2.2,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV") ;
  lt4a->SetTextColor(kBlack) ;
  lt4a->SetTextSize(0.05) ;
  lt4a->DrawLatex(2.8,2.,"#eta, ALICE");
  //	lt4->DrawLatex(3.,1.8,"#pi^{0} #rightarrow #gamma#gamma");
  //	lt4->DrawLatex(2.8,1.8,"#pi^{0} #rightarrow #gamma#gamma, #pi^{0} #rightarrow e^{+}e^{-}#gamma");
  lt4->Draw() ;

  TLegend* legendSpectraDiffDetMinBiasStripEta = new TLegend(0.18,0.83,0.5,0.95);
  legendSpectraDiffDetMinBiasStripEta->SetFillColor(0);
  legendSpectraDiffDetMinBiasStripEta->SetLineColor(0);
  legendSpectraDiffDetMinBiasStripEta->SetTextFont(42);
  legendSpectraDiffDetMinBiasStripEta->AddEntry(histoRatioCombEtaPCM,Form("PCM"),"pf");
  legendSpectraDiffDetMinBiasStripEta->AddEntry(histoRatioCombEtaEMCal,Form("EMCal"),"pf");
  legendSpectraDiffDetMinBiasStripEta->Draw();

  canvasRatioIndEtaYieldpPb->Update();
  canvasRatioIndEtaYieldpPb->Print(Form("%s/RatioIndEtaYieldFitpPb.%s",outputDir.Data(),suffix.Data()));

  //------------------ Apply x Shift to Combined Eta Spectrum------------------


   
  //    TF1*  fitTsallisEtapPb5023GeVPtXShift=(TF1*)fitTsallisEtapPb5023GeVPt->Clone("FitXShiftEtaTsallis");

  cout << "Par0: !!!!!!!!!"<<      fitTsallisEtapPb5023GeVPt->GetParameter(0)<< endl;
  TF1* fitCombEtapPb5023GeVPtXShift;
  if (FittingType.CompareTo("Tsallis")==0)
    {
      fitCombEtapPb5023GeVPtXShift=FitObject("tmpt","fitInvCrossSectionEtaXShift","Eta");
      fitCombEtapPb5023GeVPtXShift->SetParameters(fitCombEtapPb5023GeVPt->GetParameter(0),fitCombEtapPb5023GeVPt->GetParameter(1),fitCombEtapPb5023GeVPt->GetParameter(2));
      fitCombEtapPb5023GeVPtXShift->SetRange(0.7,20.);
    } else if(FittingType.CompareTo("Bylinkin2")==0)  {
    fitCombEtapPb5023GeVPtXShift=FitObject("2tcmpt","fitInvCrossSectionEtaXShift","Eta");
    fitCombEtapPb5023GeVPtXShift->SetParameters(fitCombEtapPb5023GeVPt->GetParameter(0),fitCombEtapPb5023GeVPt->GetParameter(1),fitCombEtapPb5023GeVPt->GetParameter(2),fitCombEtapPb5023GeVPt->GetParameter(3),fitCombEtapPb5023GeVPt->GetParameter(4));
    fitCombEtapPb5023GeVPtXShift->SetRange(0.7,20.);
    fitCombEtapPb5023GeVPtXShift->SetParLimits(0,0.,5.);
    fitCombEtapPb5023GeVPtXShift->SetParLimits(1,0.,.15);
  } else {fitCombEtapPb5023GeVPtXShift=FitObject("tcmpt","fitInvCrossSectionEtaXShift","Eta");
    fitCombEtapPb5023GeVPtXShift->SetParameters(fitCombEtapPb5023GeVPt->GetParameter(0),fitCombEtapPb5023GeVPt->GetParameter(1),fitCombEtapPb5023GeVPt->GetParameter(2),fitCombEtapPb5023GeVPt->GetParameter(3),fitCombEtapPb5023GeVPt->GetParameter(4));
    fitCombEtapPb5023GeVPtXShift->SetRange(0.7,20.);
    fitCombEtapPb5023GeVPtXShift->SetParLimits(0,0.,5.);
    fitCombEtapPb5023GeVPtXShift->SetParLimits(1,0.,.15);
  }
  cout << "X shift!!!!"<<endl;
  TGraphAsymmErrors* graphCombEtaInvCrossSectionStatpPb5023GeVXShifted    = (TGraphAsymmErrors*) graphCombEtaInvCrossSectionStatpPb5023GeV->Clone(); 
  cout << "X shift!!!!"<<endl;	
  TGraphAsymmErrors* graphCombEtaInvCrossSectionSyspPb5023GeVXShifted     = (TGraphAsymmErrors*) graphCombEtaInvCrossSectionSyspPb5023GeV->Clone();
  TGraphAsymmErrors* graphCombEtaInvCrossSectionTotpPb5023GeVXShifted     = (TGraphAsymmErrors*) graphCombEtaInvCrossSectionTotpPb5023GeV->Clone();
  TGraphAsymmErrors* graphCombEtaInvCrossSectionTotpPb5023GeVXShiftedtest = (TGraphAsymmErrors*) graphCombEtaInvCrossSectionTotpPb5023GeV->Clone(); 

  RemoveScalingWithPtGraph(graphCombEtaInvCrossSectionTotpPb5023GeVXShiftedtest);
  TCanvas *c1=new TCanvas("c","c",900,700);
  c1->SetLogy();
  c1->SetLogx();
  graphCombEtaInvCrossSectionTotpPb5023GeVXShiftedtest->Draw();
  graphCombEtaInvCrossSectionTotpPb5023GeVXShiftedtest->GetXaxis()->SetRangeUser(0.3,20.);
  // graphCombEtaInvCrossSectionTotpPb5023GeVXShiftedtest->Fit(fitTsallisEtapPb5023GeVPtXShift,"R");
  //   fitTsallisEtapPb5023GeVPt->Draw("same");
  graphCombEtaInvCrossSectionTotpPb5023GeVXShiftedtest->Fit(fitCombEtapPb5023GeVPtXShift,"R");
  fitCombEtapPb5023GeVPtXShift->Draw("same");

  c1->Update();     c1->Print(Form("%s/XShiftTest.%s",outputDir.Data(),suffix.Data()));


  graphCombEtaInvCrossSectionTotpPb5023GeVXShifted = ApplyXshift(graphCombEtaInvCrossSectionTotpPb5023GeVXShifted ,fitCombEtapPb5023GeVPtXShift,"Eta");
  graphCombEtaInvCrossSectionSyspPb5023GeVXShifted = ApplyXshiftIndividualSpectra(graphCombEtaInvCrossSectionTotpPb5023GeVXShifted,graphCombEtaInvCrossSectionSyspPb5023GeVXShifted , fitCombEtapPb5023GeVPtXShift, 0, graphCombEtaInvCrossSectionTotpPb5023GeVXShifted->GetN(),"Eta"); 
  graphCombEtaInvCrossSectionStatpPb5023GeVXShifted  = ApplyXshiftIndividualSpectra(graphCombEtaInvCrossSectionTotpPb5023GeVXShifted,graphCombEtaInvCrossSectionStatpPb5023GeVXShifted , fitCombEtapPb5023GeVPtXShift, 0, graphCombEtaInvCrossSectionTotpPb5023GeVXShifted->GetN(),"Eta");

  // graphCombEtaInvCrossSectionTotpPb5023GeVXShifted = ApplyXshift(graphCombEtaInvCrossSectionTotpPb5023GeVXShifted ,fitTsallisEtapPb5023GeVPtXShift,"Eta");
  // graphCombEtaInvCrossSectionSyspPb5023GeVXShifted = ApplyXshiftIndividualSpectra(graphCombEtaInvCrossSectionTotpPb5023GeVXShifted,graphCombEtaInvCrossSectionSyspPb5023GeVXShifted , fitTsallisEtapPb5023GeVPtXShift, 0, graphCombEtaInvCrossSectionTotpPb5023GeVXShifted->GetN(),"Eta"); 
  // graphCombEtaInvCrossSectionStatpPb5023GeVXShifted  = ApplyXshiftIndividualSpectra(graphCombEtaInvCrossSectionTotpPb5023GeVXShifted,graphCombEtaInvCrossSectionStatpPb5023GeVXShifted , fitTsallisEtapPb5023GeVPtXShift, 0, graphCombEtaInvCrossSectionTotpPb5023GeVXShifted->GetN(),"Eta");


  TCanvas* canvasXShiftEta = new TCanvas("canvasXShiftEta","",200,10,500,400);  // gives the page size
  DrawGammaCanvasSettings( canvasXShiftEta, 0.15, 0.02, 0.02, 0.12);
    
  canvasXShiftEta->SetLogx();
  canvasXShiftEta->SetLogy();
  TH2F * histo2DCompCombinedXShiftEta;
  histo2DCompCombinedXShiftEta = new TH2F("histo2DCompCombinedXShiftEta","histo2DCompCombinedXShiftEta",1000,0.3,30.,1000,1.2e-9,30.   );
  SetStyleHistoTH2ForGraphs(histo2DCompCombinedXShiftEta, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.03,0.04, 0.03,0.04, 1.0,1.4, 512, 508);



  histo2DCompCombinedXShiftEta->DrawCopy();  
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionTotpPb5023GeV,21,0.7, 4, 4, 1, kTRUE);  
  graphCombEtaInvCrossSectionTotpPb5023GeV->Draw("pe,same");

  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionTotpPb5023GeVXShifted,20,0.7, 2, 2, 1, kTRUE);  
  graphCombEtaInvCrossSectionTotpPb5023GeVXShifted->Draw("pe,same");
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionSyspPb5023GeVXShifted,20,0.7, 4, 4, 1, kTRUE); 
     
      
  // graphCombEtaInvCrossSectionSyspPb5023GeVXShifted->Draw("E2,same");

  // DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionStatpPb5023GeVXShifted,20,0.7, kBlue+1, kBlue+1);  
  // graphCombEtaInvCrossSectionStatpPb5023GeVXShifted->Draw("p,same");
  // //      fitCombEtapPb5023GeVPtXShift->Draw("same");
  fitCombEtapPb5023GeVPt->Draw("same");
     
  TLegend* legendRpPbCombineXShiftEta = new TLegend(0.23,0.25,0.65,0.40);
  legendRpPbCombineXShiftEta->SetFillColor(0);
  legendRpPbCombineXShiftEta->SetLineColor(0);
  legendRpPbCombineXShiftEta->SetTextSize(0.03);
  legendRpPbCombineXShiftEta->AddEntry(graphCombEtaInvCrossSectionTotpPb5023GeV,"#eta, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  legendRpPbCombineXShiftEta->AddEntry(graphCombEtaInvCrossSectionTotpPb5023GeVXShifted,"#eta x bin shifted, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  if (FittingType.CompareTo("Tsallis")==0) legendRpPbCombineXShiftEta->AddEntry(fitCombEtapPb5023GeVPtXShift,"Tsallis Fit","l");
  else legendRpPbCombineXShiftEta->AddEntry(fitCombEtapPb5023GeVPtXShift,"Bylinkin-Rostovtsev Fit","l");
    
  legendRpPbCombineXShiftEta->Draw();
	


  canvasXShiftEta->Update();
  canvasXShiftEta->Print(Form("%s/CombinedEtaSpectrum_XShifted.%s",outputDir.Data(),suffix.Data()));

  Double_t* XEtaUnshifted=graphCombEtaInvCrossSectionTotpPb5023GeV->GetX();
  Double_t* XEtaShifted=graphCombEtaInvCrossSectionTotpPb5023GeVXShifted->GetX();

  for (int i=0;i<graphCombEtaInvCrossSectionTotpPb5023GeV->GetN();i++){

    fFits <<  XEtaUnshifted[i]-graphCombEtaInvCrossSectionTotpPb5023GeV->GetErrorXlow(i)<< " - " <<  XEtaUnshifted[i]+ graphCombEtaInvCrossSectionTotpPb5023GeV->GetErrorXhigh(i)<< " GeV/$c$ & " << XEtaUnshifted[i] <<" GeV/$c$ & " << XEtaShifted[i]<<" GeV/$c$ & " << (XEtaShifted[i]-XEtaUnshifted[i])*1000 << " MeV/$c$"<< endl;
  }

  //*******************************yShift for eta/pi0*******************************************************

  // //PCM:	

  TGraphAsymmErrors* graphPCMYieldEtapPbYShift=(TGraphAsymmErrors*)graphPCMYieldEtapPb->Clone("graphPCMYieldEtapPbYShift");
   TGraphAsymmErrors* graphPCMYieldPi0EtaBinningpPbYShift=(TGraphAsymmErrors*)graphPCMYieldPi0EtaBinningpPb->Clone("graphPCMYieldPi0EtaBinningpPbYShift");


  graphPCMYieldEtapPbYShift  = ApplyYshiftIndividualSpectra(graphPCMYieldEtapPbYShift,fitCombEtapPb5023GeVPt); 
  graphPCMYieldPi0EtaBinningpPbYShift   = ApplyYshiftIndividualSpectra(graphPCMYieldPi0EtaBinningpPbYShift,fitCombPi0pPb5023GeVPt );

  //EMCal:
 cout <<"!!!!!!!!!!!!!!!!!!!EMCal!!!!!!!!!!!!!!!!!!!!"<< endl;

   graphEMCalYieldEtapPb->RemovePoint(0);
   graphEMCalYieldEtapPb->RemovePoint(0);
   graphEMCalYieldEtapPb->RemovePoint(0);
   graphEMCalYieldEtapPb->RemovePoint(0);
   graphEMCalYieldEtapPb->RemovePoint(0);
   graphEMCalYieldEtapPb->RemovePoint(0);
   graphEMCalYieldEtapPb->RemovePoint(0);
   graphEMCalYieldEtapPb->RemovePoint(0);
   graphEMCalYieldEtapPb->RemovePoint(10);
   graphEMCalYieldPi0EtaBinningpPb->RemovePoint(0);
   graphEMCalYieldPi0EtaBinningpPb->RemovePoint(0);
   graphEMCalYieldPi0EtaBinningpPb->RemovePoint(0);
   graphEMCalYieldPi0EtaBinningpPb->RemovePoint(0);
   graphEMCalYieldPi0EtaBinningpPb->RemovePoint(0);
   graphEMCalYieldPi0EtaBinningpPb->RemovePoint(0);
   graphEMCalYieldPi0EtaBinningpPb->RemovePoint(0);
   graphEMCalYieldPi0EtaBinningpPb->RemovePoint(0);
   graphEMCalYieldPi0EtaBinningpPb->RemovePoint(10);
  graphEMCalYieldEtapPb->Print();
graphEMCalYieldPi0EtaBinningpPb->Print();
  TGraphAsymmErrors* graphEMCalYieldEtapPbYShift=(TGraphAsymmErrors*)graphEMCalYieldEtapPb->Clone("graphEMCalYieldEtapPbYShift");
  TGraphAsymmErrors* graphEMCalYieldPi0EtaBinningpPbYShift=(TGraphAsymmErrors*)graphEMCalYieldPi0EtaBinningpPb->Clone("graphEMCalYieldPi0EtaBinningpPbYShift");

  graphEMCalYieldEtapPbYShift  = ApplyYshiftIndividualSpectra(graphEMCalYieldEtapPbYShift,fitCombEtapPb5023GeVPt); 
  graphEMCalYieldPi0EtaBinningpPbYShift   = ApplyYshiftIndividualSpectra(graphEMCalYieldPi0EtaBinningpPbYShift,fitCombPi0pPb5023GeVPt );

  TCanvas* canvasYShiftEtaPCM = new TCanvas("canvasYShiftEtaPCM","",200,10,500,400);  // gives the page size
  DrawGammaCanvasSettings( canvasYShiftEtaPCM, 0.15, 0.02, 0.02, 0.12);
    
  canvasYShiftEtaPCM->SetLogx();
  canvasYShiftEtaPCM->SetLogy();
  TH2F * histo2DCompCombinedYShiftEtaPCM;
  histo2DCompCombinedYShiftEtaPCM = new TH2F("histo2DCompCombinedYShiftEtaPCM","histo2DCompCombinedYShiftEtaPCM",1000,0.3,30.,1000,1.2e-9,30.   );
  SetStyleHistoTH2ForGraphs(histo2DCompCombinedYShiftEtaPCM, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.03,0.04, 0.03,0.04, 1.0,1.4, 512, 508);

  histo2DCompCombinedYShiftEtaPCM->DrawCopy(); 

  fitCombEtapPb5023GeVPt->Draw("same");
     fitCombPi0pPb5023GeVPt->Draw("same"); 
  DrawGammaSetMarkerTGraphAsym(graphPCMYieldEtapPb,21,0.7, 4, 4, 1, kTRUE);  
  graphPCMYieldEtapPb->Draw("pe,same");
 
  DrawGammaSetMarkerTGraphAsym(graphPCMYieldEtapPbYShift,20,0.7, 2, 2, 1, kTRUE);  
  graphPCMYieldEtapPbYShift->Draw("pe,same");

  
  DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0EtaBinningpPb,25,0.7, 4, 4, 1, kTRUE);  
  graphPCMYieldPi0EtaBinningpPb->Draw("pe,same");
 
  DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0EtaBinningpPbYShift,24,0.7, 2, 2, 1, kTRUE);  
  graphPCMYieldPi0EtaBinningpPbYShift->Draw("pe,same");
     
  TLegend* legendRpPbCombineYShiftEtaPCM = new TLegend(0.23,0.2,0.65,0.40);
  legendRpPbCombineYShiftEtaPCM->SetFillColor(0);
  legendRpPbCombineYShiftEtaPCM->SetLineColor(0);
  legendRpPbCombineYShiftEtaPCM->SetTextSize(0.03);
  legendRpPbCombineYShiftEtaPCM->AddEntry(graphPCMYieldPi0EtaBinningpPb,"PCM #pi^{0}","pef");
  legendRpPbCombineYShiftEtaPCM->AddEntry(graphPCMYieldPi0EtaBinningpPbYShift,"PCM #pi^{0} #it{y}-shifted","pef");
  legendRpPbCombineYShiftEtaPCM->AddEntry(graphPCMYieldEtapPb,"PCM #eta","pef");
  legendRpPbCombineYShiftEtaPCM->AddEntry(graphPCMYieldEtapPbYShift,"PCM #eta #it{y}-shifted","pef");
  if (FittingType.CompareTo("Tsallis")==0)  legendRpPbCombineYShiftEtaPCM->AddEntry(fitCombEtapPb5023GeVPt,"Tsallis Fit","l");
  else  legendRpPbCombineYShiftEtaPCM->AddEntry(fitCombEtapPb5023GeVPt,"Bylinkin-Rostovtsev Fit","l");
    
  legendRpPbCombineYShiftEtaPCM->Draw();
	


  canvasYShiftEtaPCM->Update();
  canvasYShiftEtaPCM->Print(Form("%s/EtaPCMSpectra_YShifted.%s",outputDir.Data(),suffix.Data()));

  TCanvas* canvasYShiftEtaEMCal = new TCanvas("canvasYShiftEtaEMCal","",200,10,500,400);  // gives the page size
  DrawGammaCanvasSettings( canvasYShiftEtaEMCal, 0.15, 0.02, 0.02, 0.12);
    
  canvasYShiftEtaEMCal->SetLogx();
  canvasYShiftEtaEMCal->SetLogy();
  TH2F * histo2DCompCombinedYShiftEtaEMCal;
  histo2DCompCombinedYShiftEtaEMCal = new TH2F("histo2DCompCombinedYShiftEtaEMCal","histo2DCompCombinedYShiftEtaEMCal",1000,0.3,30.,1000,1.2e-9,30.   );
  SetStyleHistoTH2ForGraphs(histo2DCompCombinedYShiftEtaEMCal, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.03,0.04, 0.03,0.04, 1.0,1.4, 512, 508);
  histo2DCompCombinedYShiftEtaEMCal->DrawCopy(); 
      

  fitCombEtapPb5023GeVPt->Draw("same");
    fitCombPi0pPb5023GeVPt->Draw("same");

  DrawGammaSetMarkerTGraphAsym(graphEMCalYieldEtapPb,21,0.7, 4, 4, 1, kTRUE);  
  graphEMCalYieldEtapPb->Draw("pe,same");
 
  DrawGammaSetMarkerTGraphAsym(graphEMCalYieldEtapPbYShift,20,0.7, 2, 2, 1, kTRUE);  
  graphEMCalYieldEtapPbYShift->Draw("pe,same");
  
  DrawGammaSetMarkerTGraphAsym(graphEMCalYieldPi0EtaBinningpPb,25,0.7, 4, 4, 1, kTRUE);  
  graphEMCalYieldPi0EtaBinningpPb->Draw("pe,same");
 
  DrawGammaSetMarkerTGraphAsym(graphEMCalYieldPi0EtaBinningpPbYShift,24,0.7, 2, 2, 1, kTRUE);  
  graphEMCalYieldPi0EtaBinningpPbYShift->Draw("pe,same");
     
  TLegend* legendRpPbCombineYShiftEtaEMCal = new TLegend(0.23,0.2,0.65,0.40);
  legendRpPbCombineYShiftEtaEMCal->SetFillColor(0);
  legendRpPbCombineYShiftEtaEMCal->SetLineColor(0);
  legendRpPbCombineYShiftEtaEMCal->SetTextSize(0.03);
  legendRpPbCombineYShiftEtaEMCal->AddEntry(graphEMCalYieldPi0EtaBinningpPb,"EMCal #pi^{0}","pef");
  legendRpPbCombineYShiftEtaEMCal->AddEntry(graphEMCalYieldPi0EtaBinningpPbYShift,"EMCal #pi^{0} #it{y}-shifted","pef");
  legendRpPbCombineYShiftEtaEMCal->AddEntry(graphEMCalYieldEtapPb,"EMCal #eta","pef");
  legendRpPbCombineYShiftEtaEMCal->AddEntry(graphEMCalYieldEtapPbYShift,"EMCal #eta #it{y}-shifted","pef");
  if (FittingType.CompareTo("Tsallis")==0)  legendRpPbCombineYShiftEtaEMCal->AddEntry(fitCombEtapPb5023GeVPt,"Tsallis Fit","l");
  else  legendRpPbCombineYShiftEtaEMCal->AddEntry(fitCombEtapPb5023GeVPt,"Bylinkin-Rostovtsev Fit","l");
    
  legendRpPbCombineYShiftEtaEMCal->Draw();
	


  canvasYShiftEtaEMCal->Update();
  canvasYShiftEtaEMCal->Print(Form("%s/EtaEMCalSpectra_YShifted.%s",outputDir.Data(),suffix.Data()));


  //******************************* Eta/Pi0 Ratio with Y-Shift ******************************************************

  TGraphAsymmErrors* graphEMCalEtaPi0RatioStatErrYShift=  CalculateGraphErrRatioToOtherTGraphErr(graphEMCalYieldEtapPbYShift,graphEMCalYieldPi0EtaBinningpPbYShift,kTRUE); 
  TGraphAsymmErrors* graphPCMEtaPi0RatioStatErrYShift=  CalculateGraphErrRatioToOtherTGraphErr(graphPCMYieldEtapPbYShift,graphPCMYieldPi0EtaBinningpPbYShift,kTRUE); 
  TGraphAsymmErrors* graphEMCalEtaPi0RatioSystErrYShift =(TGraphAsymmErrors*)graphEMCalEtaPi0RatiopPbSystErr->Clone("graphEMCalEtaPi0RatiopPbSystErrYShift");
  TGraphAsymmErrors* graphPCMEtaPi0RatioSystErrYShift =(TGraphAsymmErrors*)graphPCMEtaPi0RatiopPbSystErr->Clone("graphPCMEtaPi0RatiopPbSystErrYShift");
  graphEMCalEtaPi0RatiopPbSystErr->RemovePoint(10);
  graphEMCalEtaPi0RatiopPbSystErr->RemovePoint(10);
  graphEMCalEtaPi0RatiopPbSystErr->RemovePoint(10);
  graphEMCalEtaPi0RatiopPbSystErr->Print();
  graphPCMEtaPi0RatiopPbSystErr->Print();

  Double_t* PCMY= graphPCMEtaPi0RatioStatErrYShift->GetY();
  Double_t* PCMX= graphPCMEtaPi0RatioStatErrYShift->GetX();
  Double_t* EMCalY= graphEMCalEtaPi0RatioStatErrYShift->GetY();
  Double_t* EMCalX= graphEMCalEtaPi0RatioStatErrYShift->GetX();
 for (int i;i<graphPCMEtaPi0RatiopPbSystErr->GetN();i++){
   graphPCMEtaPi0RatioSystErrYShift->SetPoint(i,PCMX[i],PCMY[i]);
 }
 for (int i;i<graphEMCalEtaPi0RatiopPbSystErr->GetN();i++){
   graphEMCalEtaPi0RatioSystErrYShift->SetPoint(i,EMCalX[i],EMCalY[i]);
 }  

 TH1D *histoPCMEtaPi0RatioStatErrYShift=GraphAsymErrorsToHist_withErrors(graphPCMEtaPi0RatioStatErrYShift,"histoPCMEtaPi0RatiopPbStatErrYShift");
 TH1D *histoEMCalEtaPi0RatioStatErrYShift=GraphAsymErrorsToHist_withErrors(graphEMCalEtaPi0RatioStatErrYShift,"histoEMCalEtaPi0RatiopPbStatErrYShift");


  TString fileNameOutputWeightingEtaPi0Ratio                 = Form("%s/WeightingEtaPi0Ratio.dat",outputDir.Data());

  statErrorCollectionEtaPi0Ratio[0]          = (TH1D*)histoPCMEtaPi0RatioStatErrYShift->Clone("statErrPCMEta");
  statErrorCollectionEtaPi0Ratio[2]          = (TH1D*)histoEMCalEtaPi0RatioStatErrYShift->Clone("statErrEMCalEta");
    
  sysErrorCollectionEtaPi0Ratio[0]           = (TGraphAsymmErrors*)graphPCMEtaPi0RatioSystErrYShift->Clone("sysErrPCMEta");
  sysErrorCollectionEtaPi0Ratio[2]           = (TGraphAsymmErrors*)graphEMCalEtaPi0RatioSystErrYShift->Clone("sysErrEMCalEta");
    
    
  TGraphAsymmErrors* graphCombEtaPi0RatioStatpPb5023GeV= NULL;
  TGraphAsymmErrors* graphCombEtaPi0RatioSyspPb5023GeV = NULL;
    
  TGraphAsymmErrors* graphCombEtaPi0RatioTotpPb5023GeV = CombinePtPointsSpectraFullCorrMat(    statErrorCollectionEtaPi0Ratio,    sysErrorCollectionEtaPi0Ratio,     
												      pTLimitsEtaPi0Ratio, NtotalEtaPi0Ratio,
												      offSetsEtaPi0Ratio, offSetsSysEtaPi0Ratio,
												      graphCombEtaPi0RatioStatpPb5023GeV,graphCombEtaPi0RatioSyspPb5023GeV ,
												      fileNameOutputWeightingEtaPi0Ratio,1
												      );






                                   
  TCanvas* canvasYShiftEtaPi0Ratio = new TCanvas("canvasYShiftEtaPi0Ratio","",200,10,500,400);  // gives the page size
  DrawGammaCanvasSettings( canvasYShiftEtaPi0Ratio, 0.1, 0.02, 0.02, 0.11);

  TH2F * histo2DCompCombinedYShiftEtaPi0Ratio;
  histo2DCompCombinedYShiftEtaPi0Ratio = new TH2F("histo2DCompCombinedYShiftEtaPi0Ratio","histo2DCompCombinedYShiftEtaPi0Ratio",1000,0.3,16.,1000,0.,1.   );
  SetStyleHistoTH2ForGraphs(histo2DCompCombinedYShiftEtaPi0Ratio, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0} ", 0.03,0.04, 0.03,0.05, 1.0,.8, 512, 508);
  histo2DCompCombinedYShiftEtaPi0Ratio->DrawCopy(); 
      


  DrawGammaSetMarkerTGraphAsym(graphPCMEtaPi0RatioStatErrYShift,20,0.7, 1, 1);  
  graphPCMEtaPi0RatioStatErrYShift->Draw("pe,same");
  DrawGammaSetMarkerTGraphAsym(graphPCMEtaPi0RatioSystErrYShift,20,0.7, 1, 1, 1, kTRUE);  
  graphPCMEtaPi0RatioSystErrYShift->Draw("E2,same");
 
  DrawGammaSetMarkerTGraphAsym(graphEMCalEtaPi0RatioStatErrYShift,20,0.7, kGreen+2, kGreen+2);  
  graphEMCalEtaPi0RatioStatErrYShift->Draw("pe,same");
  DrawGammaSetMarkerTGraphAsym(graphEMCalEtaPi0RatioSystErrYShift,20,0.7, kGreen+2, kGreen+2, 1, kTRUE);  
  graphEMCalEtaPi0RatioSystErrYShift->Draw("E2,same");
  
  DrawGammaSetMarkerTGraphAsym(graphCombEtaPi0RatioStatpPb5023GeV,24,0.7,4, 4);  
  graphCombEtaPi0RatioStatpPb5023GeV->Draw("pe,same");
  DrawGammaSetMarkerTGraphAsym(graphCombEtaPi0RatioSyspPb5023GeV,24,0.7, 4, 4, 1, kTRUE);  
  graphCombEtaPi0RatioSyspPb5023GeV->Draw("E2,same");
  
 
     
  TLegend* legendRpPbCombineYShiftEtaPi0Ratio = new TLegend(0.15,0.75,0.5,0.90);
  legendRpPbCombineYShiftEtaPi0Ratio->SetFillColor(0);
  legendRpPbCombineYShiftEtaPi0Ratio->SetLineColor(0);
  legendRpPbCombineYShiftEtaPi0Ratio->SetTextSize(0.03);
  legendRpPbCombineYShiftEtaPi0Ratio->AddEntry(graphCombEtaPi0RatioSyspPb5023GeV,"#eta/#pi^{0} combined","pef");
  legendRpPbCombineYShiftEtaPi0Ratio->AddEntry(graphPCMEtaPi0RatioSystErrYShift,"#eta/#pi^{0} PCM","pef");
  legendRpPbCombineYShiftEtaPi0Ratio->AddEntry(graphEMCalEtaPi0RatioSystErrYShift,"#eta/#pi^{0} EMCal","pef");
 

  
  legendRpPbCombineYShiftEtaPi0Ratio->Draw();
	
  TLatex * lt4b = new TLatex(9.,0.9,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV") ;
  lt4b->SetTextColor(kBlack) ;
  lt4b->SetTextSize(0.04) ;
  lt4b->Draw() ;

  canvasYShiftEtaPi0Ratio->Update();
  canvasYShiftEtaPi0Ratio->Print(Form("%s/EtaPi0Ratio_YShifted.%s",outputDir.Data(),suffix.Data()));



  fFits.close();

  //--------------------Write results to file---------------------

  TFile fResults(Form("%s/ResultspPb_%s_%s.root",outputDir.Data(),FittingType.Data(), dateForOutput.Data()),"RECREATE");
  fitCombPi0pPb5023GeVPt->Write("FitCombinedPi0pPbSpectrum");
  graphCombPi0InvCrossSectionStatpPb5023GeV->Write("CombinedPi0pPbSpectrumStatErr");
  graphCombPi0InvCrossSectionSyspPb5023GeV->Write("CombinedPi0pPbSpectrumSysErr");
  graphCombPi0InvCrossSectionTotpPb5023GeV->Write("CombinedPi0pPbSpectrumTotErr");
  graphCombPi0InvCrossSectionStatpPb5023GeVXShifted->Write("CombinedPi0pPbSpectrumStatErrXShifted");
  graphCombPi0InvCrossSectionSyspPb5023GeVXShifted->Write("CombinedPi0pPbSpectrumSysErrXShifted");
  graphCombPi0InvCrossSectionTotpPb5023GeVXShifted->Write("CombinedPi0pPbSpectrumTotErrXShifted");
  fitCombEtapPb5023GeVPt->Write("FitCombinedEtapPbSpectrum");
  graphCombEtaInvCrossSectionStatpPb5023GeV->Write("CombinedEtapPbSpectrumStatErr");
  graphCombEtaInvCrossSectionSyspPb5023GeV->Write("CombinedEtapPbSpectrumSysErr");
  graphCombEtaInvCrossSectionTotpPb5023GeV->Write("CombinedEtapPbSpectrumTotErr");
  graphCombEtaInvCrossSectionStatpPb5023GeVXShifted->Write("CombinedEtapPbSpectrumStatErrXShifted");
  graphCombEtaInvCrossSectionSyspPb5023GeVXShifted->Write("CombinedEtapPbSpectrumSysErrXShifted");
  graphCombEtaInvCrossSectionTotpPb5023GeVXShifted->Write("CombinedEtapPbSpectrumTotErrXShifted");
  graphCombEtaPi0RatioStatpPb5023GeV->Write("CombinedEtaPi0RatioStatErrYShifted");
  graphCombEtaPi0RatioSyspPb5023GeV->Write("CombinedEtaPi0RatioSysErrYShifted");
  graphCombEtaPi0RatioTotpPb5023GeV->Write("CombinedEtaPi0RatioTotErrYShifted");
}
