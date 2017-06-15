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
TGraphAsymmErrors* ApplyNSDSysError(TGraphAsymmErrors* graphMesonSystErr, Double_t ScalingErr){
  TGraphAsymmErrors*   graphMesonSystErrClone=(TGraphAsymmErrors*) graphMesonSystErr->Clone();
  Int_t nPoints=graphMesonSystErrClone->GetN();
  Double_t* valueY    = graphMesonSystErrClone->GetY();
  Double_t* valueX    = graphMesonSystErrClone->GetX();
  Double_t* ErrYlowSys     = graphMesonSystErrClone->GetEYlow();
  Double_t* ErrYhighSys   = graphMesonSystErrClone->GetEYhigh();
  Double_t* ErrXlowSys     = graphMesonSystErrClone->GetEXlow();
  Double_t* ErrXhighSys   = graphMesonSystErrClone->GetEXhigh();

  for(Int_t i = 0; i < graphMesonSystErrClone->GetN(); i++){


    ErrYlowSys[i]  = TMath::Sqrt(   (ErrYlowSys[i] *ErrYlowSys[i]  ) + ( (valueY[i]*ScalingErr) *(valueY[i]*ScalingErr) ) );
    ErrYhighSys[i]  = TMath::Sqrt(   (ErrYhighSys[i] *ErrYhighSys[i]  ) + ( (valueY[i]*ScalingErr) *(valueY[i]*ScalingErr) ) );
        
  }
  TGraphAsymmErrors* graphMesonSystErrClone1= new TGraphAsymmErrors(nPoints,valueX,valueY,ErrXlowSys,ErrXhighSys,ErrYlowSys,ErrYhighSys);
      
  return graphMesonSystErrClone1;
}

//const Int_t  Ntotal = 31;
const Int_t  Ntotal = 28;
const Int_t  NtotalLow = 25;
const Int_t  nPtLimits = Ntotal+1;
const Int_t  nPtLimitsLow = NtotalLow+1;
const Int_t  NtotalEta = 16;
const Int_t  nPtLimitsEta = NtotalEta+1;
const Int_t  NtotalEtaPi0Ratio = 15;
const Int_t  nPtLimitsEtaPi0Ratio = NtotalEtaPi0Ratio+1;

void CombineMesonMeasurementspPb5023GeV(TString FittingType = "Tsallis",Bool_t IsNSD=kTRUE, Bool_t LowPtCut=kFALSE){    

  TString date                                        = ReturnDateString();
    
  gROOT->Reset();
  gROOT->SetStyle("Plain");
    
  StyleSettingsThesis();
  SetPlotStyle();
    
  TString suffix                                      = "eps";
    
  TString dateForOutput                               = ReturnDateStringForOutput();
  TString outputDir                                   = Form("CombinepPbSpectra/%s/%s/%s",dateForOutput.Data(),FittingType.Data(),suffix.Data());
  if (LowPtCut && IsNSD) outputDir                    = Form("CombinepPbSpectra/%s_NSD_LowPtCut/%s/%s",dateForOutput.Data(),FittingType.Data(),suffix.Data());
  else  if (LowPtCut) outputDir                       = Form("CombinepPbSpectra/%s_LowPtCut/%s/%s",dateForOutput.Data(),FittingType.Data(),suffix.Data());
  else if (IsNSD) outputDir                           = Form("CombinepPbSpectra/%s_NSD/%s/%s",dateForOutput.Data(),FittingType.Data(),suffix.Data());
  gSystem->Exec("mkdir -p "+outputDir);
  fstream fFits;   
  TString nameFits = Form("%s/FittingParameters_%s_%s.dat",outputDir.Data(),dateForOutput.Data(),FittingType.Data());
  fFits.open(nameFits.Data(), ios::out);

  Double_t ScalingErr = 0.031;
  Double_t Scaling = 0.964; 
  Double_t ScalingFit =1./Scaling; 
  //___________________________________ Declaration of files _____________________________________________
  TString fileNameNeutralPionDalitz                   = "ExternalInputpPb/PCM/data_PCMResults_Dalitz_pPb_20170606.root";
  //TString fileNameNeutralPionDalitz                   = "ExternalInputpPb/PCM/data_PCMResults_Dalitz_pPb_20160929.root";
  //TString fileNameNeutralPionPCM                      = "ExternalInputpPb/PCM/data_PCMResults_pPb_20151111_standard_CatErrors.root";
  //TString fileNameNeutralPionPCM                      = "ExternalInputpPb/PCM/data_PCMResults_pPb_20170308.root";
  TString fileNameNeutralPionPCM                      = "ExternalInputpPb/PCM/data_PCMResults_pPb_20170606_oldformat_WithSystematics.root";
  //TString fileNameNeutralPionPCM                      = "ExternalInputpPb/PCM/data_PCMResultsFullCorrection_pPb_20170602.root";
  //TString fileNameNeutralPionPCMOnlyErr               = "ExternalInputpPb/PCM/data_PCMResults_pPb_20170606_oldformat_woEffiCorrection.root";
  //TString fileNameNeutralPionPCM                      = "ExternalInputpPb/PCM/data_PCMResultsFullCorrection_pPb_20170606_woEffiCorrection.root";
  //TString fileNameNeutralPionEMCal                    = "ExternalInputpPb/EMCAL/data_EMCalEMCalResults_160602_newhighptbinningPi0_pPb.root";
  //TString fileNameNeutralPionEMCal 		      = "ExternalInputpPb/EMCAL/data_EMCalEMCalResults_161118_pPb.root";
  //TString fileNameNeutralPionEMCal		      = "ExternalInputpPb/EMCAL/data_EMCalEMCalResults_170602_pPb.root";
  TString fileNameNeutralPionEMCal                    = "ExternalInputpPb/EMCAL/data_EMCalEMCalResults_170607_pPb.root";
 
  //TString fileNamePCMEMCal			      = "ExternalInputpPb/PCM-EMCAL/data_PCM-EMCALResultsFullCorrection_pPb_2016_12_22.root";
  TString fileNamePCMEMCal                            = "ExternalInputpPb/PCM-EMCAL/data_PCM-EMCALResultsFullCorrection_pPb_2017_06_13.root";
 // TString fileNameNeutralPionEMCalEta                 = "ExternalInputpPb/EMCAL/data_EMCalEMCalResults_160503_pPb.root"; 
 // TString fileNameNeutralPionEMCalEta                 = "ExternalInputpPb/EMCAL/data_EMCalEMCalResults_161118_pPb.root";
  //TString fileNameNeutralPionEMCalEta                 = "ExternalInputpPb/EMCAL/data_EMCalEMCalResults_161124_pPb.root";
  					
  TString fileNameNeutralPionEMCalEta                 = "ExternalInputpPb/EMCAL/data_EMCalEMCalResults_170607_pPb.root";
  TString fileNameNeutralPionPHOS                     = "ExternalInputpPb/PHOS/20160601_Pi0InvariantSpectrum_pPb_PHOS.root";//data_PHOSResults_pPb_20160208.root";
  TString fileNameChargedPions                        = "ExternalInputpPb/InputRpPb/pPb502.fullpT.INEL.20151204_mb_wo_V0Acorr.root";//Charged pion pPb spectrum 
  TString fileNamemTScalingpPb                        = "CombinepPbSpectra/mT-Scaling/2017_06_14/mTScaling.root"; 
  TString fileNamemTScalingpp7TeV                     = "CombinepPbSpectra/mT-Scaling/2017_06_14/mTScaling_pp7TeV.root";
  TString fileNameChargedPionsRatioPythia             = "ExternalInputpPb/Pythia6_pionRatio.root";
  
  

  TString nameHistoPHOS                               = "hCor_stat";
  TString nameHistoPHOSSysErrors                      = "hCor_syst";
 
  TString nameHistoEMCal                              = "CorrectedYieldPi0";
  TString nameHistoEMCalSysErrors                     = "Pi0SystError";
  
  //Correlation factors file for PCM and Dalitz
  //TString corrFactorsFile = "ExternalInputpPb/CorrelationFactors_PCM_Dalitz_2016_10_06_pPb_5.023TeV.root";
  TString corrFactorsFile = "ExternalInputpPb/CorrFactorsFiles/pPb5TeV_2017_01_17.root";
    
    
  TString collisionSystempPb                          = "p-Pb #sqrt{s_{NN}} = 5.02 TeV"; 
    
    

//   Double_t pTLimits[nPtLimits]                        = { 0.3, 0.4, 0.5, 0.6, 0.7, 
// 							  0.8, 1.0, 1.2, 1.4, 1.6, 
// 							  1.8, 2.0, 2.2, 2.4, 2.6,
// 							  2.8, 3.0, 3.2, 3.4, 3.6,
// 							  3.8, 4.0, 4.5, 5.0, 5.5,
// 							  6.0, 7.0, 8.0, 10.0,12.0,
// 							  16.0, 20.0};
                                                          
                                                          
  
  Double_t pTLimits[nPtLimits]                        = { 0.6, 0.7, 
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
   
  Int_t offSets[11]                                   =  { -4, 2, -4, 0, -4,  0, 0, 0, 0,  0, 0};
  Int_t offSetsSys[11]                                =  {  0, 3,  5, 0,  2,   0, 0, 0, 0, 0, 0};

  //Int_t offSets[11]                                   =  { -1, 5, -1, 0, -1,  3, 0, 0, 0,  0, 0};
  //Int_t offSetsSys[11]                                =  {  0, 6,  8, 0, 5, 3, 0, 0, 0, 0, 0};   
  
  
  Int_t offSetsLow[11]                                =  { -7, -1, -7, 0, 0,  -4, 0, 0, 0,  0, 0};
  Int_t offSetsSysLow[11]                             =  {  -6, 0,  2, 0, 0, -3, 0, 0, 0, 0, 0};  
  Int_t offSetsEta[11]                                =  { -3, 0, -3, 0, -3,  0, 0, 0, 0,  0, 0};
  Int_t offSetsSysEta[11]                             =  {  0, 0,  5, 0, 2, 0, 0, 0, 0, 0, 0};
  Int_t offSetsEtaPi0Ratio[11]                        =  { 0, 0, 5, 0, -3,  0, 0, 0, 0,  0, 0};
  Int_t offSetsSysEtaPi0Ratio[11]                     =  {  0, 0,  5, 0, 2, 0, 0, 0, 0, 0, 0};
    
    
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
  TFile* fileNeutralPionsDalitz                    = new TFile(fileNameNeutralPionDalitz.Data());
  TDirectory* directoryDalitzPi0pPb                   = (TDirectory*) fileNeutralPionsDalitz->GetDirectory("Pi0_pPb_5.023TeV_0-100%");
  if( ! directoryDalitzPi0pPb ){
    cout<<"Dalitz: The directory Pi0_pPb_5.023TeV_0-100% does not exist "<<endl; 
    return;
  }

  TH1D* histoDalitzYieldPi0pPbStat                    = (TH1D*)directoryDalitzPi0pPb->Get("CorrectedYieldPi0");
  TGraphAsymmErrors* graphDalitzYieldPi0pPbSystErr    = (TGraphAsymmErrors*)directoryDalitzPi0pPb->Get("Pi0SystError");
  
  TGraphAsymmErrors* graphRatioCombDalitzSys          = (TGraphAsymmErrors*)graphDalitzYieldPi0pPbSystErr->Clone();
  TH1D* histoRatioCombDalitzStat          	      = (TH1D*)histoDalitzYieldPi0pPbStat->Clone();
 
  TGraphAsymmErrors* graphRatioCombDalitzStat          = new TGraphAsymmErrors(histoDalitzYieldPi0pPbStat);
    
  graphRatioCombDalitzStat->RemovePoint(0);
  //graphRatioCombDalitzStat->RemovePoint(graphRatioCombDalitzStat->GetN()-1);
  //graphDalitzYieldPi0pPbSystErr->RemovePoint(graphDalitzYieldPi0pPbSystErr->GetN()-1);
  
  TH1D* histoDalitzYieldPi0pPbStat2=GraphAsymErrorsToHist_withErrors(graphRatioCombDalitzStat, "histoDalitzYieldPi0pPbStat2"); //Just for the combined
  

  cout<<"Dalitz systematic"<<endl;
  graphDalitzYieldPi0pPbSystErr->Print();
  cout<<"Dalitz statistic2"<<endl;
  graphRatioCombDalitzStat->Print();
  
  // **************************************************************************************
  // ****************************** Reading PCM *******************************************
  // **************************************************************************************    
  TFile* fileNeutralPionPCM                           = new TFile(fileNameNeutralPionPCM.Data());
  
  TDirectory* directoryPCMPi0pPb                      = (TDirectory*) fileNeutralPionPCM->GetDirectory("Pi0_pPb_5.023TeV_0-100%");    
  if( ! directoryPCMPi0pPb ){
    cout<<"PCM: The directory Pi0pPb 5.023TeV does not exist "<<endl; 
    return; 
  }
  
  TDirectory*	directoryPCMEtapPb 			= (TDirectory*)fileNeutralPionPCM->Get("Eta_pPb_5.023TeV_0-100%"); 
  if( ! directoryPCMEtapPb ){
    cout<<"PCM: The directory Eta_pPb_5.023TeV_0-100% does not exist "<<endl; 
    return; 
  }
  
  
  //TFile* fileNeutralPionPCMOnlyErr                 = new TFile(fileNameNeutralPionPCMOnlyErr.Data());
  
  //TDirectory*	directoryPCMPi0EtaBinningpPb 			= (TDirectory*)fileNeutralPionPCMOnlyErr->Get("Pi0EtaBinning_pPb_5.023TeV_0-100%");
  
  
  TDirectory*   directoryPCMPi0EtaBinningpPb                    = (TDirectory*)fileNeutralPionPCM->Get("Pi0EtaBinning_pPb_5.023TeV_0-100%");
  
  if( ! directoryPCMPi0EtaBinningpPb ){
    cout<<"PCM: The directory Pi0EtaBinning_pPb_5.023TeV_0-100% does not exist "<<endl; 
    return; 
  }
  
  
  
  
  
  TH1D* histoPCMYieldPi0pPbStat                       = (TH1D*)directoryPCMPi0pPb->Get("CorrectedYieldPi0");
  
   ///NOTE Temporal, start for 4 point
  
  histoPCMYieldPi0pPbStat->SetBinContent(2,.0);
  histoPCMYieldPi0pPbStat->SetBinContent(3,.0);
  histoPCMYieldPi0pPbStat->SetBinContent(4,.0);  
  histoPCMYieldPi0pPbStat->SetBinError(2,.0);
  histoPCMYieldPi0pPbStat->SetBinError(3,.0);
  histoPCMYieldPi0pPbStat->SetBinError(4,.0);
  
 
  for(Int_t iPoint = 1; iPoint < histoPCMYieldPi0pPbStat->GetNbinsX(); iPoint++) {
     cout<<iPoint<<" "<<histoPCMYieldPi0pPbStat->GetBinCenter(iPoint)<<" "<<histoPCMYieldPi0pPbStat->GetBinContent(iPoint)<<endl;
      
  }
  
 ///NOTE Temporal, start for 4 point
  
  TGraphAsymmErrors* graphPCMYieldPi0pPbSystErr       = (TGraphAsymmErrors*)directoryPCMPi0pPb->Get("Pi0SystError");
  
  graphPCMYieldPi0pPbSystErr->RemovePoint(0);
  graphPCMYieldPi0pPbSystErr->RemovePoint(0);
  graphPCMYieldPi0pPbSystErr->RemovePoint(0);
  
  
  //cout<<"For temporal, start for 4 point"<<endl;
  //graphPCMYieldPi0pPbSystErr->Print();
  //cout<<"//////////////////////"<<endl;
 // for(Int_t iPoint = 0; iPoint < graphPCMYieldPi0pPbSystErr->GetN(); iPoint++)
 // graphPCMYieldPi0pPbSystErr->SetPointError (Int_t i, Double_t exl, Double_t exh, Double_t eyl, Double_t eyh)
  
  
  
  
  
  
  
  TGraphAsymmErrors* temp01                           = new TGraphAsymmErrors(histoPCMYieldPi0pPbStat);

  TH1D*  histoPCMYieldEtapPbStat 			= (TH1D*)directoryPCMEtapPb->Get("CorrectedYieldEta");  
  TGraphAsymmErrors*  graphPCMYieldEtapPbStat 		= new TGraphAsymmErrors(histoPCMYieldEtapPbStat);
  TGraphAsymmErrors*	graphPCMYieldEtapPbSystErr	= (TGraphAsymmErrors*)directoryPCMEtapPb->Get("EtaSystError");	
  //TGraphAsymmErrors*	graphPCMYieldEtapPbSystErrA	= (TGraphAsymmErrors*)directoryPCMEtapPb->Get("EtaSystErrorA");	
  TH1D*	histoPCMEtaPi0RatiopPb 				= (TH1D*)directoryPCMEtapPb->Get("EtaToPi0StatError");
  
//  TGraphAsymmErrors*	graphPCMEtaPi0RatiopPbSystErr 	= (TGraphAsymmErrors*)directoryPCMEtapPb->Get("EtaToPi0SystError");
  
   TGraphAsymmErrors*	graphPCMEtaPi0RatiopPbSystErr 	= (TGraphAsymmErrors*)directoryPCMEtapPb->Get("EtatoPi0RatioSys");


  TH1D*  histoPCMYieldPi0EtaBinningpPb 			= (TH1D*)directoryPCMPi0EtaBinningpPb->Get("CorrectedYieldPi0EtaBinning");
  TGraphAsymmErrors*  graphPCMYieldPi0EtaBinningpPb 				= new TGraphAsymmErrors(histoPCMYieldPi0EtaBinningpPb);
   
    
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
  graphPCMYieldEtapPbStat->RemovePoint(0);
  graphPCMYieldEtapPbStat->RemovePoint(0);
  graphPCMYieldEtapPbStat->RemovePoint(0);
  graphPCMYieldEtapPbStat->Print(); 
  //cout<<"PCM Eta systematicA"<<endl;
  //graphPCMYieldEtapPbSystErrA->Print(); 

  TGraphAsymmErrors* graphRatioCombPCMSys             = (TGraphAsymmErrors*)graphPCMYieldPi0pPbSystErr->Clone();
  TH1D* histoRatioCombPCMStat                         = (TH1D*)histoPCMYieldPi0pPbStat->Clone();
  TGraphAsymmErrors* graphRatioCombEtaPCMSys          = (TGraphAsymmErrors*)graphPCMYieldEtapPbSystErr->Clone();
  TH1D* histoRatioCombEtaPCMStat                      = (TH1D*)histoPCMYieldEtapPbStat->Clone(); 
  
  
    
    
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

  TGraphAsymmErrors* graphRatioCombPHOSSys            = (TGraphAsymmErrors*)graphPHOSYieldPi0pPbSystErr->Clone();
  TH1D* histoRatioCombPHOSStat                        = (TH1D*)histoPHOSYieldPi0pPbStat->Clone();  
  // **************************************************************************************
  // ******************************** Reading EMCal ***************************************
  // **************************************************************************************
  TFile* fileNeutralPionEMCal                         = new TFile(fileNameNeutralPionEMCal);
  TFile* fileNeutralPionEMCalEta                      = new TFile(fileNameNeutralPionEMCalEta);
  TDirectory* directoryEMCalPi0pPb                    = (TDirectory*)fileNeutralPionEMCal->Get("Pi0_pPb_5.023TeV_0-100%");
  if( ! directoryEMCalPi0pPb ){
    cout<<"EMCal: The directory Pi0_pPb_5.023TeV_0-100% does not exist "<<endl; 
    return;  
  }
  TDirectory*	directoryEMCalEtapPb 	      	= (TDirectory*)fileNeutralPionEMCalEta->Get("Eta_pPb_5.023TeV_0-100%"); 
  if( ! directoryEMCalEtapPb  ){
    cout<<"EMCal: The directory Eta_pPb_5.023TeV_0-100% does not exist "<<endl; 
    return;  
  }  TDirectory*	directoryEMCalPi0EtaBinningpPb 	      	= (TDirectory*)fileNeutralPionEMCal->Get("Pi0EtaBinning_pPb_5.023TeV_0-100%"); 
  if( ! directoryEMCalPi0EtaBinningpPb  ){
    cout<<"EMCal: The directory Pi0EtaBinning_pPb_5.023TeV_0-100% does not exist "<<endl; 
    return;  
  }
  //Pi0
  TH1D* histoEMCalYieldPi0pPbStat1                     = (TH1D*)directoryEMCalPi0pPb->Get(nameHistoEMCal.Data());
  TGraphAsymmErrors* graphEMCalYieldPi0pPbStatErr     =new TGraphAsymmErrors(histoEMCalYieldPi0pPbStat1);

  cout <<"EMCal StatErr"<<endl;
  
 // for(Int_t iPoint = 0; iPoint < 9; iPoint++)
      //graphEMCalYieldPi0pPbStatErr->RemovePoint(0);
  //graphEMCalYieldPi0pPbStatErr->Print();
  graphEMCalYieldPi0pPbStatErr->RemovePoint(31);
  graphEMCalYieldPi0pPbStatErr->Print();
  TH1D* histoEMCalYieldPi0pPbStat=GraphAsymErrorsToHist_withErrors(graphEMCalYieldPi0pPbStatErr, "histoEMCalYieldPi0pPbStat");
  TGraphAsymmErrors* graphEMCalYieldPi0pPbSystErr     = (TGraphAsymmErrors*)directoryEMCalPi0pPb->Get(nameHistoEMCalSysErrors.Data());
  cout<<"EMCal SystErr"<<endl;
  //graphEMCalYieldPi0pPbSystErr->Print();
  graphEMCalYieldPi0pPbSystErr->RemovePoint(22);
  graphEMCalYieldPi0pPbSystErr->RemovePoint(22);
  graphEMCalYieldPi0pPbSystErr->RemovePoint(22);
  graphEMCalYieldPi0pPbSystErr->RemovePoint(22);
  graphEMCalYieldPi0pPbSystErr->Print();
  cout<<"/////////////////////////////"<<endl;
  //Eta
  TH1D*	histoEMCalYieldEtapPbStat 			= (TH1D*)directoryEMCalEtapPb->Get("CorrectedYieldEta");
  TGraphAsymmErrors*	graphEMCalYieldEtapPbSystErr		= (TGraphAsymmErrors*)directoryEMCalEtapPb->Get("EtaSystError");
  
  TGraphAsymmErrors* graphEMCalYieldEtapPb=new TGraphAsymmErrors(histoEMCalYieldEtapPbStat); 

  //Pi0EtaBinning
  TH1D* histoEMCalYieldPi0EtaBinningpPbStat  = (TH1D*)directoryEMCalPi0EtaBinningpPb->Get("CorrectedYieldPi0EtaBinning");
  TGraphAsymmErrors*   graphEMCalYieldPi0EtaBinningpPb=new TGraphAsymmErrors(histoEMCalYieldPi0EtaBinningpPbStat);
  //EtaPi0Ratio	
  TH1D*	histoEMCalEtaPi0RatiopPb 			= (TH1D*)directoryEMCalEtapPb->Get("EtatoPi0Ratio");
  TGraphAsymmErrors*	graphEMCalEtaPi0RatiopPbSystErr 	        = (TGraphAsymmErrors*)directoryEMCalEtapPb->Get("EtatoPi0RatioSys"); 
  cout<<"EMCal Eta systematic"<<endl;
  graphEMCalYieldEtapPbSystErr->Print();
  cout<<"EMCal Eta statistical"<<endl;
  graphEMCalYieldEtapPb->Print();        

  TGraphAsymmErrors* graphRatioCombEMCalSys               = (TGraphAsymmErrors*)graphEMCalYieldPi0pPbSystErr->Clone();
  TH1D* histoRatioCombEMCalStat                        	  = (TH1D*)histoEMCalYieldPi0pPbStat; 
  TGraphAsymmErrors* graphRatioCombEtaEMCalSys            = (TGraphAsymmErrors*)graphEMCalYieldEtapPbSystErr->Clone();
  TH1D* histoRatioCombEtaEMCalStat                        = (TH1D*)histoEMCalYieldPi0pPbStat;  
  
  // **************************************************************************************
  // ******************************* Reading PCM-EMCAL ************************************
  // **************************************************************************************    
  
  TFile* fileNeutralMesonsPCMEMCal                       = new TFile(fileNamePCMEMCal.Data());
  TDirectory* directoryPCMEMCalPi0pPb                    = (TDirectory*)fileNeutralMesonsPCMEMCal->Get("Pi0pPb_5.023TeV");
  if( ! directoryPCMEMCalPi0pPb ){
    cout<<"PCM-EMCal: The directory Pi0pPb_5.023TeV does not exist "<<endl; 
    return;  
  }
  TDirectory* directoryPCMEMCalEtapPb                    = (TDirectory*)fileNeutralMesonsPCMEMCal->Get("EtapPb_5.023TeV");
  if( ! directoryPCMEMCalEtapPb ){
    cout<<"PCM-EMCal: The directory EtapPb_5.023TeV does not exist "<<endl; 
    return;  
  }
  
    //Pi0
  TH1D* histoPCMEMCalYieldPi0pPbStat                  = (TH1D*)directoryPCMEMCalPi0pPb->Get("CorrectedYieldPi0");
  TGraphAsymmErrors* graphPCMEMCalYieldPi0pPbStatErr  = new TGraphAsymmErrors(histoPCMEMCalYieldPi0pPbStat);
   
  TGraphAsymmErrors* graphPCMEMCalYieldPi0pPbSystErr  = (TGraphAsymmErrors*)directoryPCMEMCalPi0pPb->Get("Pi0SystError");
  
  cout<<"PCM-EMCal statistic  erros Pi0"<<endl;
  graphPCMEMCalYieldPi0pPbStatErr->Print();
  cout<<"PCM-EMCal systematic erros Pi0"<<endl;
  graphPCMEMCalYieldPi0pPbSystErr->Print();
  
  TH1D* histoPCMEMCalYieldEtapPbStat                  = (TH1D*)directoryPCMEMCalEtapPb->Get("CorrectedYieldEta");
  TGraphAsymmErrors* graphPCMEMCalYieldEtapPbStatErr  = new TGraphAsymmErrors(histoPCMEMCalYieldEtapPbStat);
  TGraphAsymmErrors* graphPCMEMCalYieldEtapPbSystErr  = (TGraphAsymmErrors*)directoryPCMEMCalEtapPb->Get("EtaSystError");
  
  cout<<"PCM-EMCal statistic  erros Eta"<<endl;
  graphPCMEMCalYieldEtapPbStatErr->Print();
  cout<<"PCM-EMCal systematic erros Eta"<<endl;
  graphPCMEMCalYieldEtapPbSystErr->Print();
  
  //Pi0 Eta binning
 
  TH1D* histoPCMEMCalEtaPi0RatioStatErrYShift                     = (TH1D*)directoryPCMEMCalEtapPb->Get("EtaToPi0YShiftedStatError");
  TGraphAsymmErrors* graphPCMEMCalEtaPi0RatioSystErrYShift        = (TGraphAsymmErrors*)directoryPCMEMCalEtapPb->Get("EtaToPi0YShiftedSystError");
  
  if (!histoPCMEMCalEtaPi0RatioStatErrYShift){
            histoPCMEMCalEtaPi0RatioStatErrYShift                 = (TH1D*)directoryPCMEMCalEtapPb->Get("EtaToPi0StatError");
  } 
  if (!graphPCMEMCalEtaPi0RatioSystErrYShift){
            graphPCMEMCalEtaPi0RatioSystErrYShift                 = (TGraphAsymmErrors*)directoryPCMEMCalEtapPb->Get("EtaToPi0SystError");
  } 
  
  
  TGraphAsymmErrors* graphRatioCombPCMEMCalSys        = (TGraphAsymmErrors*)graphPCMEMCalYieldPi0pPbSystErr->Clone();
  TH1D* histoRatioCombPCMEMCalStat                    = (TH1D*)histoPCMEMCalYieldPi0pPbStat->Clone();
  TGraphAsymmErrors* graphRatioCombEtaPCMEMCalSys     = (TGraphAsymmErrors*)graphPCMEMCalYieldEtapPbSystErr->Clone();
  TH1D* histoRatioCombEtaPCMEMCalStat                 = (TH1D*)histoPCMEMCalYieldEtapPbStat->Clone(); 
  
  // **************************************************************************************
  // ******************************** Reading mT Scaling **********************************
  // **************************************************************************************  
  TFile* filemTScaling                         	       = new TFile(fileNamemTScalingpPb);
  TGraphAsymmErrors* graphEtaPi0Ratio_mTScaled_pPb     =(TGraphAsymmErrors*)filemTScaling->Get("graph_EtaPi0Ratio_mTscaled");
  TGraphAsymmErrors* graphEtaPi0Ratio_vsmT_Sys         =(TGraphAsymmErrors*)filemTScaling->Get("EtaFitPi0Ratio_vs_mT_SysErr");
  TGraphAsymmErrors* graphEtaPi0Ratio_vsmT_Stat        =(TGraphAsymmErrors*)filemTScaling->Get("EtaFitPi0Ratio_vs_mT_StatErr");
  TF1* EtaPi0Ratio_mTScaled_pPb                        =(TF1*)filemTScaling->Get("func_EtaPi0Ratio_mTscaled");
  TF1* EtaSpectrum_mTScaled_pPb                        =(TF1*)filemTScaling->Get("func_EtaSpectrum_mTscaled");
  TFile* filemTScalingpp7TeV                           = new TFile(fileNamemTScalingpp7TeV);
  TGraphAsymmErrors* graphEtaPi0Ratio_mTScaled_pp7TeV  =(TGraphAsymmErrors*)filemTScalingpp7TeV->Get("graph_EtaPi0Ratio_mTscaled");
  TGraphAsymmErrors* graphEtaPi0Ratio_vsmT_Sys_pp7TeV  =(TGraphAsymmErrors*)filemTScalingpp7TeV->Get("EtaFitPi0Ratio_vs_mT_SysErr");
  TGraphAsymmErrors* graphEtaPi0Ratio_vsmT_Stat_pp7TeV =(TGraphAsymmErrors*)filemTScalingpp7TeV->Get("EtaFitPi0Ratio_vs_mT_StatErr");
  TF1* EtaPi0Ratio_mTScaled_pp7TeV                     =(TF1*)filemTScalingpp7TeV->Get("func_EtaPi0Ratio_mTscaled");
  TF1* EtaSpectrum_mTScaled_pp7TeV                     =(TF1*)filemTScalingpp7TeV->Get("func_EtaSpectrum_mTscaled");
  // **************************************************************************************
  // ******************************** Reading Charged Pions *******************************
  // **************************************************************************************  
  TFile* fileChargedPions                         = new TFile(fileNameChargedPions);
  TH1D* histoChargedPionspPbStatErr               = (TH1D*)fileChargedPions->Get("hstat_pPb502_mb_pion_sum");
  histoChargedPionspPbStatErr->Sumw2();
  histoChargedPionspPbStatErr->Scale(0.5);
  TH1D* histoChargedPionspPbSystErr               = (TH1D*)fileChargedPions->Get("hsys_pPb502_mb_pion_sum"); 
  histoChargedPionspPbSystErr->Sumw2();
  histoChargedPionspPbSystErr->Scale(0.5);
  TGraphAsymmErrors* graphChargedPionspPbSyst     =new TGraphAsymmErrors(histoChargedPionspPbSystErr);
  TGraphAsymmErrors* graphChargedPionspPbStat     =new TGraphAsymmErrors(histoChargedPionspPbStatErr);

  // **************************************************************************************
  // ********************************* Combine Pi0 spectra ********************************
  // **************************************************************************************
    
  TString fileNameOutputWeightingPi0                  = Form("%s/WeightingPi0.dat",outputDir.Data());

  statErrorCollection[0]          = (TH1D*)histoPCMYieldPi0pPbStat->Clone("statErrPCMPi0");
  statErrorCollection[1]          = (TH1D*)histoPHOSYieldPi0pPbStat->Clone("statErrPHOSPi0");
  statErrorCollection[2]          = (TH1D*)histoEMCalYieldPi0pPbStat->Clone("statErrEMCalPi0");
  statErrorCollection[4]          = (TH1D*)histoPCMEMCalYieldPi0pPbStat->Clone("statErrPCMEMCalPi0");
  statErrorCollection[5]          = (TH1D*)histoDalitzYieldPi0pPbStat2->Clone("statErrDalitzPi0");
    
  sysErrorCollection[0]           = (TGraphAsymmErrors*)graphPCMYieldPi0pPbSystErr->Clone("sysErrPCMPi0");
  sysErrorCollection[1]           = (TGraphAsymmErrors*)graphPHOSYieldPi0pPbSystErr->Clone("sysErrPHOSPi0");
  sysErrorCollection[2]           = (TGraphAsymmErrors*)graphEMCalYieldPi0pPbSystErr->Clone("sysErrEMCalPi0");
  sysErrorCollection[4]	          = (TGraphAsymmErrors*)graphPCMEMCalYieldPi0pPbSystErr->Clone("sysErrPCMEMCalPi0");
  sysErrorCollection[5]           = (TGraphAsymmErrors*)graphDalitzYieldPi0pPbSystErr->Clone("sysErrDalitzPi0");
    
    
    
    
  TGraphAsymmErrors* graphCombPi0InvCrossSectionStatpPb5023GeV= NULL;
  TGraphAsymmErrors* graphCombPi0InvCrossSectionSyspPb5023GeV = NULL;
   
  cout<<"\n\n */////////////////Combination of Pi0 meson begins here**//////////////"<<endl;
  
  TGraphAsymmErrors* graphCombPi0InvCrossSectionTotpPb5023GeV=NULL;
  if(!LowPtCut)graphCombPi0InvCrossSectionTotpPb5023GeV= CombinePtPointsSpectraFullCorrMat(    statErrorCollection,    sysErrorCollection,     
											       pTLimits, Ntotal,
											       offSets, offSetsSys,
											       graphCombPi0InvCrossSectionStatpPb5023GeV, graphCombPi0InvCrossSectionSyspPb5023GeV,
											       fileNameOutputWeightingPi0,"pPb_5.023TeV","Pi0",kFALSE,NULL,corrFactorsFile.Data());

  else graphCombPi0InvCrossSectionTotpPb5023GeV= CombinePtPointsSpectraFullCorrMat(    statErrorCollection,    sysErrorCollection,     
										       pTLimitsLow, NtotalLow,
										       offSetsLow, offSetsSysLow,
										       graphCombPi0InvCrossSectionStatpPb5023GeV, graphCombPi0InvCrossSectionSyspPb5023GeV,
										       fileNameOutputWeightingPi0,"pPb_5.023TeV","Pi0",kFALSE,NULL,corrFactorsFile.Data());
 // cout<< "CombinedSpectrum:" << endl;
  //graphCombPi0InvCrossSectionTotpPb5023GeV->Print();
  //graphCombPi0InvCrossSectionStatpPb5023GeV->Print();
  //graphCombPi0InvCrossSectionSyspPb5023GeV->Print();
  
  cout<<"\n\n */////////////////Combination of Pi0 meson ends here**//////////////"<<endl;

  if (IsNSD){

   graphCombPi0InvCrossSectionStatpPb5023GeV=ScaleGraph(graphCombPi0InvCrossSectionStatpPb5023GeV,Scaling);
   graphCombPi0InvCrossSectionSyspPb5023GeV=ScaleGraph(graphCombPi0InvCrossSectionSyspPb5023GeV,Scaling);

   graphCombPi0InvCrossSectionSyspPb5023GeV=ApplyNSDSysError(graphCombPi0InvCrossSectionSyspPb5023GeV,ScalingErr);

   graphCombPi0InvCrossSectionTotpPb5023GeV  =  CalculateCombinedSysAndStatError( graphCombPi0InvCrossSectionStatpPb5023GeV ,graphCombPi0InvCrossSectionSyspPb5023GeV );
    
   cout<< "CombinedSpectrum NSD scaled:" << endl;
   graphCombPi0InvCrossSectionTotpPb5023GeV->Print();
   graphCombPi0InvCrossSectionStatpPb5023GeV->Print();
   graphCombPi0InvCrossSectionSyspPb5023GeV->Print();
  }


  TGraphAsymmErrors* graphInvYieldPi0CombpPb5023GeVStaClone   = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionStatpPb5023GeV->Clone();
  TGraphAsymmErrors* graphInvYieldPi0CombpPb5023GeVSysClone   = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionSyspPb5023GeV->Clone();
  TGraphAsymmErrors* graphInvYieldPi0CombpPb5023GeVTotClone   = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionTotpPb5023GeV->Clone();
    
  TGraphAsymmErrors* graphRatioCombCombFitTot                 = (TGraphAsymmErrors*) graphInvYieldPi0CombpPb5023GeVTotClone ->Clone(); 
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
  Double_t minPt = 0.6; //0.3
  Double_t maxPt = 20.0;
  if (FittingType.CompareTo("Tsallis")==0){
    fitCombPi0pPb5023GeVPt = FitObject("l","fitInvCrossSectionPi0","Pi0");
    ParameterspPb[0] = 8.66870e+00;//7.4e+1;
    ParameterspPb[1] = 7.03073e+00;
    ParameterspPb[2] = 1.63511e-01;
    fitCombPi0pPb5023GeVPt->SetRange(minPt,maxPt);
    fitCombPi0pPb5023GeVPt->SetParameters(ParameterspPb[0],ParameterspPb[1],ParameterspPb[2]);
    //  fitCombPi0pPb5023GeVPt->SetParLimits(2,0.17,0.18);
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
    for (Int_t i=0;i < 3;i++){    
      fFits<< fitCombPi0pPb5023GeVPt->GetParName(i) <<":  "<< fitCombPi0pPb5023GeVPt->GetParameter(i)<<" +/- "<< fitCombPi0pPb5023GeVPt->GetParError(i)  << endl; 
    }   
  }else {
    for (Int_t i=0;i < 5;i++){    
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
  if (IsNSD)  legendRpPbCombine->AddEntry(graphInvYieldPi0CombpPb5023GeVSysClone,"NSD #pi^{0}, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  else legendRpPbCombine->AddEntry(graphInvYieldPi0CombpPb5023GeVSysClone,"#pi^{0}, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  if(FittingType.CompareTo("Tsallis")==0) legendRpPbCombine->AddEntry(fitCombPi0pPb5023GeVPt,"Tsallis Fit","l");
  else legendRpPbCombine->AddEntry(fitCombPi0pPb5023GeVPt,"Bylinkin-Rostovtsev Fit","l");
    
  legendRpPbCombine->Draw();
	
    
    
    
  fitCombPi0pPb5023GeVPt->Draw("same");
    
  canvasCompYieldpPbInd->Print(Form("%s/Comb_pPb.%s",outputDir.Data(),suffix.Data()));

  // **************************************************************************************
  // ************************* Creating ratios of different systems ************************
  // **************************************************************************************
  TF1* fitCombPi0pPb5023GeVPt_woNSD=(TF1*)fitCombPi0pPb5023GeVPt->Clone("Fit_woNSD");
  if (IsNSD){
    if (FittingType.CompareTo("Tsallis")==0)fitCombPi0pPb5023GeVPt_woNSD->SetParameter(0,fitCombPi0pPb5023GeVPt_woNSD->GetParameter(0)*ScalingFit);
    else{
      fitCombPi0pPb5023GeVPt_woNSD->SetParameter(0,fitCombPi0pPb5023GeVPt_woNSD->GetParameter(0)*ScalingFit);
      fitCombPi0pPb5023GeVPt_woNSD->SetParameter(2,fitCombPi0pPb5023GeVPt_woNSD->GetParameter(2)*ScalingFit);
    }
  }
  graphRatioCombPCMSys        = CalculateGraphErrRatioToFit(graphRatioCombPCMSys,      fitCombPi0pPb5023GeVPt_woNSD);
  graphRatioCombDalitzSys     = CalculateGraphErrRatioToFit(graphRatioCombDalitzSys,   fitCombPi0pPb5023GeVPt_woNSD);
  graphRatioCombEMCalSys      = CalculateGraphErrRatioToFit(graphRatioCombEMCalSys,    fitCombPi0pPb5023GeVPt_woNSD);
  graphRatioCombPHOSSys       = CalculateGraphErrRatioToFit(graphRatioCombPHOSSys,     fitCombPi0pPb5023GeVPt_woNSD);
  graphRatioCombPCMEMCalSys   = CalculateGraphErrRatioToFit(graphRatioCombPCMEMCalSys, fitCombPi0pPb5023GeVPt_woNSD);
   
  graphRatioCombCombFitTot    = CalculateGraphErrRatioToFit(graphInvYieldPi0CombpPb5023GeVTotClone,fitCombPi0pPb5023GeVPt); 
  graphRatioCombCombFitSta    = CalculateGraphErrRatioToFit(graphInvYieldPi0CombpPb5023GeVStaClone,fitCombPi0pPb5023GeVPt); 
  graphRatioCombCombFitSys    = CalculateGraphErrRatioToFit(graphInvYieldPi0CombpPb5023GeVSysClone,fitCombPi0pPb5023GeVPt); 

  histoRatioCombPCMStat       = CalculateHistoRatioToFit(histoPCMYieldPi0pPbStat,	fitCombPi0pPb5023GeVPt_woNSD);
  histoRatioCombDalitzStat    = CalculateHistoRatioToFit(histoDalitzYieldPi0pPbStat,	fitCombPi0pPb5023GeVPt_woNSD);
  histoRatioCombEMCalStat     = CalculateHistoRatioToFit(histoEMCalYieldPi0pPbStat,	fitCombPi0pPb5023GeVPt_woNSD);
  histoRatioCombPHOSStat      = CalculateHistoRatioToFit(histoRatioCombPHOSStat,	fitCombPi0pPb5023GeVPt_woNSD);
  histoRatioCombPCMEMCalStat  = CalculateHistoRatioToFit(histoRatioCombPCMEMCalStat, 	fitCombPi0pPb5023GeVPt_woNSD);

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
  // graphRatioCombCombFitSys->SetPointEYlow(30,0.137146);//attribute EMCal Sys error to 
  // graphRatioCombCombFitSys->SetPointEYhigh(30,0.137146);
  DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys,20,1, kBlue, kBlue, 1, kTRUE, kBlue-9);  
  graphRatioCombCombFitSys->Draw("E2,same");

  DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSta,20,1, kBlue, kBlue);  
  graphRatioCombCombFitSta->Draw("pe,same");

  DrawGammaSetMarkerTGraphAsym(graphRatioCombPCMSys,24,1, 1, 1, 1, kTRUE);  
  graphRatioCombPCMSys->Draw("E2,same");

  DrawGammaSetMarker(histoRatioCombPCMStat, 24,1, 1 , 1);
  histoRatioCombPCMStat->Draw("E1,same") ;

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

  graphRatioCombCombFitSta->Draw("pe,same");

  DrawGammaSetMarkerTGraphAsym(graphRatioCombDalitzSys,24,1, kCyan+2, kCyan+2, 1, kTRUE);  
  graphRatioCombDalitzSys->Draw("E2same");
  DrawGammaSetMarker(histoRatioCombDalitzStat, 24,1 ,kCyan+2 ,kCyan+2);
  histoRatioCombDalitzStat->Draw("E1,same") ; 
  TLatex * lt = new TLatex(2.,0.7,"PCM"); 
  lt->SetTextSize(0.16) ;
  lt->SetTextColor(kCyan+2) ;
  lt->DrawText(2.,1.25,"Dalitz") ;
 
  // ===========EMCal======================  
  canvasRatioCompYieldpPbInd->cd() ;
  DrawPad(2,0.0,0.,0.12,0.05);
  ratio2DInvXSectionPi0->DrawCopy(); 
  
  graphRatioCombCombFitSys->Draw("E2") ;
  
  graphRatioCombCombFitSta->Draw("pe,same");

            
  DrawGammaSetMarkerTGraphAsym(graphRatioCombEMCalSys,24,1,kGreen+2 , kGreen+2, 1, kTRUE);  
  graphRatioCombEMCalSys->Draw("E2same");
  DrawGammaSetMarker(histoRatioCombEMCalStat, 24,1 ,kGreen+2 ,kGreen+2);
  histoRatioCombEMCalStat->Draw("E1,same") ;
  cout<<"EMCal"<<endl;
  graphRatioCombEMCalSys->Print();
  lt->SetTextColor(kGreen+2) ;
  lt->DrawText(2.,1.25,"EMCal") ;

  //PHOS  
  canvasRatioCompYieldpPbInd->cd() ;
  DrawPad(3,0.,0.,0.12,0.05);
  
  ratio2DInvXSectionPi0->DrawCopy(); 
  graphRatioCombCombFitSys->Draw("E2") ;
  graphRatioCombCombFitSta->Draw("pe,same");

  DrawGammaSetMarkerTGraphAsym(graphRatioCombPHOSSys,24,1, kRed+1, kRed+1, 1, kTRUE);  
  graphRatioCombPHOSSys->Draw("E2same");
  DrawGammaSetMarker(histoRatioCombPHOSStat, 24,1 ,kRed+1 ,kRed+1);
  histoRatioCombPHOSStat->Draw("E1,same") ; 

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
  ratio2DInvXSectionPi0a = new TH2F("ratio2DInvXSectionPi0a","ratio2DInvXSectionPi0a",1000,0.2,25.,1000,0.41,2.29);

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

  DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSta,20,1., kBlue, kBlue,1.5);  
  graphRatioCombCombFitSta->Draw("p,same");

  TLatex * lt3 = new TLatex(2.8,2.,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV") ;
  lt3->SetTextColor(kBlack) ;
  lt3->SetTextSize(0.05) ;
  if (IsNSD)  lt3->DrawLatex(2.8,1.85,"NSD #pi^{0}, ALICE");
  else lt3->DrawLatex(2.8,1.85,"#pi^{0}, ALICE");
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
  ratio2DInvXSectionPi0b = new TH2F("ratio2DInvXSectionPi0b","ratio2DInvXSectionPi0b",1000,0.2,25.,1000,0.41,2.29);

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
  DrawGammaSetMarker(histoRatioCombPCMStat, 20,1, 1 , 1);
  histoRatioCombPCMStat->Draw("pE1,same") ;

  DrawGammaSetMarkerTGraphAsym(graphRatioCombDalitzSys,29,1, kCyan+2, kCyan+2, 1, kTRUE);  
  graphRatioCombDalitzSys->Draw("E2same");
  DrawGammaSetMarker(histoRatioCombDalitzStat, 29,1.5 ,kCyan+2 ,kCyan+2);
  histoRatioCombDalitzStat->Draw("pE1,same") ; 
 
  DrawGammaSetMarkerTGraphAsym(graphRatioCombEMCalSys,33,1,kGreen+2 , kGreen+2, 1, kTRUE);  
  graphRatioCombEMCalSys->Draw("E2same");
  DrawGammaSetMarker(histoRatioCombEMCalStat, 33,1.5 ,kGreen+2 ,kGreen+2);
  histoRatioCombEMCalStat->Draw("pE1,same") ;
	
  graphRatioCombPHOSSys->Draw("E2same") ;
	
  histoRatioCombPHOSStat->SetMarkerStyle(21);
  histoRatioCombPHOSStat->SetMarkerSize(1);
  histoRatioCombPHOSStat->Draw("pE1,same");
	
  TLatex * lt4 = new TLatex(2.8,2.1,"ALICE Preliminary") ;
  lt4->SetTextColor(kBlack) ;
  lt4->SetTextSize(0.05) ;
  if (IsNSD)   lt4->DrawLatex(2.8,1.95,"p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  if (IsNSD)  lt4->DrawLatex(2.8,1.8,"NSD #pi^{0}");
  else lt4->DrawLatex(2.8,1.95,"#pi^{0}, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  //	lt4->DrawLatex(3.,1.8,"#pi^{0} #rightarrow #gamma#gamma");
  //	lt4->DrawLatex(2.8,1.8,"#pi^{0} #rightarrow #gamma#gamma, #pi^{0} #rightarrow e^{+}e^{-}#gamma");
  lt4->Draw() ;

  TLegend* legendSpectraDiffDetMinBiasStrip = new TLegend(0.18,0.7,0.5,0.95);
  legendSpectraDiffDetMinBiasStrip->SetFillColor(0);
  legendSpectraDiffDetMinBiasStrip->SetLineColor(0);
  legendSpectraDiffDetMinBiasStrip->SetTextFont(42);
  legendSpectraDiffDetMinBiasStrip->AddEntry(histoRatioCombPCMStat,Form("PCM"),"pf");
  legendSpectraDiffDetMinBiasStrip->AddEntry(histoRatioCombDalitzStat,Form("Dalitz"),"pf");
  legendSpectraDiffDetMinBiasStrip->AddEntry(histoRatioCombEMCalStat,Form("EMCal"),"pf");
  legendSpectraDiffDetMinBiasStrip->AddEntry(histoRatioCombPHOSStat,Form("PHOS"),"pf");
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
  if (IsNSD) { 
    legendRpPbCombineXShift->AddEntry(graphCombPi0InvCrossSectionTotpPb5023GeV,"NSD #pi^{0} , p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
    legendRpPbCombineXShift->AddEntry(graphCombPi0InvCrossSectionTotpPb5023GeVXShifted,"NSD #pi^{0} x bin shifted, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  } else {
    legendRpPbCombineXShift->AddEntry(graphCombPi0InvCrossSectionTotpPb5023GeV,"#pi^{0} , p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
    legendRpPbCombineXShift->AddEntry(graphCombPi0InvCrossSectionTotpPb5023GeVXShifted,"#pi^{0} x bin shifted, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  }
  if (FittingType.CompareTo("Tsallis")==0)  legendRpPbCombineXShift->AddEntry(fitCombPi0pPb5023GeVPt,"Tsallis","l");
  else  legendRpPbCombineXShift->AddEntry(fitCombPi0pPb5023GeVPt,"Bylinkin-Rostovtsev Fit","l");
    
  legendRpPbCombineXShift->Draw();
	


  canvasXShift->Update();
  canvasXShift->Print(Form("%s/CombinedSpectrum_XShifted.%s",outputDir.Data(),suffix.Data()));



  Double_t* XPi0Unshifted=graphCombPi0InvCrossSectionTotpPb5023GeV->GetX();
  Double_t* XPi0Shifted=graphCombPi0InvCrossSectionTotpPb5023GeVXShifted->GetX();

  for (Int_t i=0;i<graphCombPi0InvCrossSectionTotpPb5023GeV->GetN();i++){

    fFits <<  XPi0Unshifted[i]-graphCombPi0InvCrossSectionTotpPb5023GeV->GetErrorXlow(i)<< " - " <<  XPi0Unshifted[i]+ graphCombPi0InvCrossSectionTotpPb5023GeV->GetErrorXhigh(i)<< " GeV/$c$ & " << XPi0Unshifted[i] <<" GeV/$c$ & " << XPi0Shifted[i]<<" GeV/$c$ & " << (XPi0Shifted[i]-XPi0Unshifted[i])*1000 <<" MeV/$c$" << endl;
  }



  // **************************************************************************************
  // ************************* Plotting x shifted spectrum + ratio ************************
  // **************************************************************************************

  TGraphAsymmErrors*  graphRatioCombXshiftFitSta    = CalculateGraphErrRatioToFit(graphCombPi0InvCrossSectionStatpPb5023GeVXShifted,fitCombPi0pPb5023GeVPt); 
  TGraphAsymmErrors*   graphRatioCombXShiftFitSys    = CalculateGraphErrRatioToFit(graphCombPi0InvCrossSectionSyspPb5023GeVXShifted,fitCombPi0pPb5023GeVPt); 



  TCanvas* CanvasRatioPi0 = new TCanvas("CanvasRatioPi0","",200,10,1200,1200);  // gives the page size
  DrawGammaCanvasSettings(CanvasRatioPi0 , 0.3, 0.02, 0.02, 0.16);
  TPad* padHistos = new TPad("padHistos", "", 0., 0.25, 1., 1.,-1, -1, -2);
  DrawGammaPadSettings( padHistos, 0.15, 0.02, 0.02, 0.);
  padHistos->Draw();

  TPad* padRatios = new TPad("padRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
  DrawGammaPadSettings( padRatios, 0.15, 0.02, 0., 0.2);
  padRatios->Draw();

  padHistos->cd(); 
  padHistos->SetLogx();
  padHistos->SetLogy();
  TH2F * histoRatioPi0;
  histoRatioPi0 = new TH2F("histoRatioPi0","histoRatioPi0",1000,0.2,30.,1000,1.2e-9,30.);
  SetStyleHistoTH2ForGraphs(histoRatioPi0,"#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.04,0.05, 0.045,0.045, 0.7,1.4, 512, 505);
  histoRatioPi0->DrawCopy();
  fitCombPi0pPb5023GeVPt->Draw("same");
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionSyspPb5023GeVXShifted,20,0.7, 4, 4, 1, kTRUE);  
      
  graphCombPi0InvCrossSectionSyspPb5023GeVXShifted->Draw("E2,same");

  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionStatpPb5023GeVXShifted,20,0.7, kBlue, kBlue);  
  graphCombPi0InvCrossSectionStatpPb5023GeVXShifted->Draw("p,same");

  TLatex * ltPre = new TLatex(.35,1e-6,"ALICE Preliminary"); 
  ltPre->SetTextSize(0.05) ;
  //  lt->DrawText(1.,1.45,"2016/06/14") ;
  ltPre->Draw("same");     
  TLegend* legendRpPbCombineRatio = new TLegend(0.23,0.10,0.7,0.250);
  legendRpPbCombineRatio->SetFillColor(0);
  legendRpPbCombineRatio->SetLineColor(0);
  legendRpPbCombineRatio->SetTextSize(0.04);
  if (IsNSD)  legendRpPbCombineRatio->AddEntry(graphInvYieldPi0CombpPb5023GeVSysClone,"NSD #pi^{0}, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  else  legendRpPbCombineRatio->AddEntry(graphInvYieldPi0CombpPb5023GeVSysClone,"#pi^{0}, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  if(FittingType.CompareTo("Tsallis")==0) legendRpPbCombineRatio->AddEntry(fitCombPi0pPb5023GeVPt,"Tsallis Fit","l");
  else legendRpPbCombineRatio->AddEntry(fitCombPi0pPb5023GeVPt,"Bylinkin-Rostovtsev Fit","l");
    
  legendRpPbCombineRatio->Draw();
  
  padRatios->cd();
    
  padRatios->SetLogx();
  //  canvasRatioAllppreferences->SetLogy();
  TH2F * histo2DRatioAllppreferences2;
  histo2DRatioAllppreferences2 = new TH2F("histo2DRatioAllppreferences2","histo2DRatioAllppreferences",1000,.2,30.,1000,0.41,1.89);
  SetStyleHistoTH2ForGraphs(histo2DRatioAllppreferences2, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.13,0.13, 0.13,0.13, 0.5,.5, 502, 503);
  histo2DRatioAllppreferences2->DrawCopy();

  DrawGammaSetMarkerTGraphAsym(graphRatioCombXShiftFitSys,20,0.7, 4, 4, 1, kTRUE);  
      
  graphRatioCombXShiftFitSys->Draw("E2,same");

  DrawGammaSetMarkerTGraphAsym(graphRatioCombXshiftFitSta,20,0.7, kBlue, kBlue);  
  graphRatioCombXshiftFitSta->Draw("p,same");  

 
  TLine *lineA2=new TLine(0.,1.,25.,1.);
  lineA2->SetLineColor(kGray+1);
  lineA2->Draw();   

    
  CanvasRatioPi0->Print(Form("%s/Pi0xShift_ratio.%s",outputDir.Data(),suffix.Data()));

  //------------------ Apply y Shift to Combined Spectrum------------------


  TF1* fitCombPi0pPb5023GeVPtYShift=(TF1*)fitCombPi0pPb5023GeVPt->Clone("FitYShift");
  cout << "Y shift!!!!"<<endl;
  TGraphAsymmErrors* graphCombPi0InvCrossSectionStatpPb5023GeVYShifted  = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionStatpPb5023GeV->Clone(); 
  cout << "Y shift!!!!"<<endl;	
  TGraphAsymmErrors* graphCombPi0InvCrossSectionSyspPb5023GeVYShifted  = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionSyspPb5023GeV->Clone();
  TGraphAsymmErrors* graphCombPi0InvCrossSectionTotpPb5023GeVYShifted = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionTotpPb5023GeV->Clone();
  graphCombPi0InvCrossSectionTotpPb5023GeVYShifted = ApplyYshiftIndividualSpectra(graphCombPi0InvCrossSectionTotpPb5023GeVYShifted ,fitCombPi0pPb5023GeVPtYShift);
  graphCombPi0InvCrossSectionSyspPb5023GeVYShifted = ApplyYshiftIndividualSpectra(graphCombPi0InvCrossSectionSyspPb5023GeVYShifted,fitCombPi0pPb5023GeVPtYShift); 
  graphCombPi0InvCrossSectionStatpPb5023GeVYShifted  = ApplyYshiftIndividualSpectra(graphCombPi0InvCrossSectionStatpPb5023GeVYShifted,fitCombPi0pPb5023GeVPtYShift);

  TCanvas* canvasYShift = new TCanvas("canvasYShift","",200,10,500,400);  // gives the page size
  DrawGammaCanvasSettings( canvasYShift, 0.15, 0.02, 0.02, 0.12);
    
  canvasYShift->SetLogx();
  canvasYShift->SetLogy();
  TH2F * histo2DCompCombinedYShift;
  histo2DCompCombinedYShift = new TH2F("histo2DCompCombinedYShift","histo2DCompCombinedYShift",1000,0.3,30.,1000,1.2e-9,30.   );
  SetStyleHistoTH2ForGraphs(histo2DCompCombinedYShift, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.03,0.04, 0.03,0.04, 1.0,1.4, 512, 508);



  histo2DCompCombinedYShift->DrawCopy();  
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionTotpPb5023GeV,21,0.7, 4, 4, 1, kTRUE);  
  graphCombPi0InvCrossSectionTotpPb5023GeV->Draw("pe,same");

  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionTotpPb5023GeVYShifted,20,0.7, 2, 2, 1, kTRUE);  
  graphCombPi0InvCrossSectionTotpPb5023GeVYShifted->Draw("pe,same");
  // DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionSyspPb5023GeVYShifted,20,0.7, 4, 4, 1, kTRUE); 
     
      
  // graphCombPi0InvCrossSectionSyspPb5023GeVYShifted->Draw("E2,same");

  // DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionStatpPb5023GeVYShifted,20,0.7, kBlue, kBlue);  
  // graphCombPi0InvCrossSectionStatpPb5023GeVYShifted->Draw("p,same");
  fitCombPi0pPb5023GeVPt->Draw("same");
     
  TLegend* legendRpPbCombineYShift = new TLegend(0.23,0.25,0.65,0.40);
  legendRpPbCombineYShift->SetFillColor(0);
  legendRpPbCombineYShift->SetLineColor(0);
  legendRpPbCombineYShift->SetTextSize(0.03);
  if (IsNSD)  legendRpPbCombineYShift->AddEntry(graphCombPi0InvCrossSectionTotpPb5023GeV,"NSD #pi^{0} , p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  else  legendRpPbCombineYShift->AddEntry(graphCombPi0InvCrossSectionTotpPb5023GeVYShifted,"#pi^{0} y bin shifted, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  if (FittingType.CompareTo("Tsallis")==0)  legendRpPbCombineYShift->AddEntry(fitCombPi0pPb5023GeVPt,"Tsallis","l");
  else  legendRpPbCombineYShift->AddEntry(fitCombPi0pPb5023GeVPt,"Bylinkin-Rostovtsev Fit","l");
    
  legendRpPbCombineYShift->Draw();
	


  canvasYShift->Update();
  canvasYShift->Print(Form("%s/CombinedSpectrum_YShifted.%s",outputDir.Data(),suffix.Data()));




  // **************************************************************************************
  // ********************************* Combine Eta spectra ********************************
  // **************************************************************************************
    
  TString fileNameOutputWeightingEta    = Form("%s/WeightingEta.dat",outputDir.Data());

  statErrorCollectionEta[0]          	= (TH1D*)histoPCMYieldEtapPbStat->Clone("statErrPCMEta");
  statErrorCollectionEta[2]          	= (TH1D*)histoEMCalYieldEtapPbStat->Clone("statErrEMCalEta");
  statErrorCollectionEta[4]	     	= (TH1D*)histoPCMEMCalYieldEtapPbStat->Clone("statErrPCMEMCalEta");
    
  sysErrorCollectionEta[0]           	= (TGraphAsymmErrors*)graphPCMYieldEtapPbSystErr->Clone("sysErrPCMEta");
  sysErrorCollectionEta[2]           	= (TGraphAsymmErrors*)graphEMCalYieldEtapPbSystErr->Clone("sysErrEMCalEta");
  sysErrorCollectionEta[4]	     	= (TGraphAsymmErrors*)graphPCMEMCalYieldEtapPbSystErr->Clone("sysErrPCMEMCalEta");
    
    
  TGraphAsymmErrors* graphCombEtaInvCrossSectionStatpPb5023GeV= NULL;
  TGraphAsymmErrors* graphCombEtaInvCrossSectionSyspPb5023GeV = NULL;
  
  cout<<"\n\n */////////////////Combination of Eta meson begins here**//////////////"<<endl;
    
  TGraphAsymmErrors* graphCombEtaInvCrossSectionTotpPb5023GeV = CombinePtPointsSpectraFullCorrMat(    statErrorCollectionEta,    sysErrorCollectionEta,     
												      pTLimitsEta, NtotalEta,
												      offSetsEta, offSetsSysEta,
												      graphCombEtaInvCrossSectionStatpPb5023GeV, graphCombEtaInvCrossSectionSyspPb5023GeV,
												      fileNameOutputWeightingEta,"pPb_5.023TeV", "Eta", kFALSE, NULL, corrFactorsFile.Data()
												 );
  
  cout<<"\n\n */////////////////Combination of Eta meson ends here**//////////////"<<endl;
  //cout<<"Comb Eta Sys"<<endl;
  //graphCombEtaInvCrossSectionSyspPb5023GeV->Print();


  if (IsNSD){

    graphCombEtaInvCrossSectionStatpPb5023GeV=ScaleGraph(graphCombEtaInvCrossSectionStatpPb5023GeV,Scaling);
    graphCombEtaInvCrossSectionSyspPb5023GeV=ScaleGraph(graphCombEtaInvCrossSectionSyspPb5023GeV,Scaling);

     graphCombEtaInvCrossSectionSyspPb5023GeV=ApplyNSDSysError(graphCombEtaInvCrossSectionSyspPb5023GeV,ScalingErr);

    graphCombEtaInvCrossSectionTotpPb5023GeV  =  CalculateCombinedSysAndStatError( graphCombEtaInvCrossSectionStatpPb5023GeV ,graphCombEtaInvCrossSectionSyspPb5023GeV );
    
    cout<< "CombinedSpectrum NSD scaled:" << endl;
    graphCombEtaInvCrossSectionTotpPb5023GeV->Print();
    graphCombEtaInvCrossSectionStatpPb5023GeV->Print();
    graphCombEtaInvCrossSectionSyspPb5023GeV->Print();
  }


 
  TGraphAsymmErrors* graphInvYieldEtaCombpPb5023GeVStaClone   = (TGraphAsymmErrors*) graphCombEtaInvCrossSectionStatpPb5023GeV->Clone();
  TGraphAsymmErrors* graphInvYieldEtaCombpPb5023GeVSysClone   = (TGraphAsymmErrors*) graphCombEtaInvCrossSectionSyspPb5023GeV->Clone();
  TGraphAsymmErrors* graphInvYieldEtaCombpPb5023GeVTotClone   = (TGraphAsymmErrors*) graphCombEtaInvCrossSectionTotpPb5023GeV->Clone();
    
  TGraphAsymmErrors* graphRatioCombEtaFitTot                    = (TGraphAsymmErrors*) graphInvYieldEtaCombpPb5023GeVTotClone ->Clone(); 
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
    for (Int_t i=0;i < 3;i++){    
      fFits<< fitCombEtapPb5023GeVPt->GetParName(i) <<":  "<< fitCombEtapPb5023GeVPt->GetParameter(i)<<" +/- "<< fitCombEtapPb5023GeVPt->GetParError(i)  << endl; 
    }   
  }else {
    for (Int_t i=0;i < 5;i++){    
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
  if (IsNSD)  legendRpPbCombineEta->AddEntry(graphInvYieldEtaCombpPb5023GeVSysClone,"NSD #eta, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  else legendRpPbCombineEta->AddEntry(graphInvYieldEtaCombpPb5023GeVSysClone,"#eta, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  if (FittingType.CompareTo("Tsallis")==0)  legendRpPbCombineEta->AddEntry(fitCombEtapPb5023GeVPt,"Tsallis Fit","l");
  else legendRpPbCombineEta->AddEntry(fitCombEtapPb5023GeVPt,"Bylinkin-Rostovtsev Fit","l");
   
    
  legendRpPbCombineEta->Draw();
	
    
    
    
  fitCombEtapPb5023GeVPt->Draw("same");
    
  canvasCompEtaYieldpPbInd->Print(Form("%s/CombEta_pPb.%s",outputDir.Data(),suffix.Data()));

  // **************************************************************************************
  // ************************* Creating ratios of different systems ************************
  // **************************************************************************************
  TF1* fitCombEtapPb5023GeVPt_woNSD=(TF1*)fitCombEtapPb5023GeVPt->Clone("FitEta_woNSD");
  if (IsNSD){
    if (FittingType.CompareTo("Tsallis")==0)fitCombEtapPb5023GeVPt_woNSD->SetParameter(0,fitCombEtapPb5023GeVPt_woNSD->GetParameter(0)*ScalingFit);
    else{
      fitCombEtapPb5023GeVPt_woNSD->SetParameter(0,fitCombEtapPb5023GeVPt_woNSD->GetParameter(0)*ScalingFit);
      fitCombEtapPb5023GeVPt_woNSD->SetParameter(2,fitCombEtapPb5023GeVPt_woNSD->GetParameter(2)*ScalingFit);
    }
  }
  graphRatioCombEtaPCMSys        = CalculateGraphErrRatioToFit(graphRatioCombEtaPCMSys,      fitCombEtapPb5023GeVPt_woNSD);
  graphRatioCombEtaEMCalSys      = CalculateGraphErrRatioToFit(graphRatioCombEtaEMCalSys,    fitCombEtapPb5023GeVPt_woNSD);
  graphRatioCombEtaPCMEMCalSys   = CalculateGraphErrRatioToFit(graphRatioCombEtaPCMEMCalSys, fitCombEtapPb5023GeVPt_woNSD);
    
  graphRatioCombEtaFitTot       = CalculateGraphErrRatioToFit(graphInvYieldEtaCombpPb5023GeVTotClone,fitCombEtapPb5023GeVPt); 
  graphRatioCombEtaFitSta    = CalculateGraphErrRatioToFit(graphInvYieldEtaCombpPb5023GeVStaClone,fitCombEtapPb5023GeVPt); 
  graphRatioCombEtaFitSys    = CalculateGraphErrRatioToFit(graphInvYieldEtaCombpPb5023GeVSysClone,fitCombEtapPb5023GeVPt); 
 

  histoRatioCombEtaPCMStat      = CalculateHistoRatioToFit(histoPCMYieldEtapPbStat,		fitCombEtapPb5023GeVPt_woNSD);
  histoRatioCombEtaEMCalStat    = CalculateHistoRatioToFit(histoEMCalYieldEtapPbStat,		fitCombEtapPb5023GeVPt_woNSD);
  histoRatioCombEtaPCMEMCalStat = CalculateHistoRatioToFit(histoPCMEMCalYieldEtapPbStat, 	fitCombEtapPb5023GeVPt_woNSD);
 
 
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

  DrawGammaSetMarkerTGraphAsym(graphRatioCombEtaFitSta,20,1, kBlue, kBlue);  
  graphRatioCombEtaFitSta->Draw("pe,same");

  DrawGammaSetMarkerTGraphAsym(graphRatioCombEtaPCMSys,24,1, 1, 1, 1, kTRUE);  
  graphRatioCombEtaPCMSys->Draw("E2,same");

  DrawGammaSetMarker(histoRatioCombEtaPCMStat, 24,1, 1 , 1);
  histoRatioCombEtaPCMStat->Draw("p,E1,same") ;

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
  
  graphRatioCombEtaFitSta->Draw("pe,same");

            
  DrawGammaSetMarkerTGraphAsym(graphRatioCombEtaEMCalSys,24,1,kGreen+2 , kGreen+2, 1, kTRUE);  
  graphRatioCombEtaEMCalSys->Draw("E2same");
  DrawGammaSetMarker(histoRatioCombEtaEMCalStat, 24,1 ,kGreen+2 ,kGreen+2);
  histoRatioCombEtaEMCalStat->Draw("p,E1,same") ;
  cout<<"EMCal"<<endl;
  graphRatioCombEtaEMCalSys->Print();
  lt->SetTextColor(kGreen+2) ;
  lt->SetTextSize(0.1) ;
  lt->DrawText(1.,1.3,"EMCal") ;
  TLatex* labelEta;
  if (IsNSD) labelEta=new TLatex(5,1.55,"NSD #eta #rightarrow #gamma#gamma") ;
  else labelEta=new TLatex(10.5,1.55,"#eta #rightarrow #gamma#gamma") ;
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
  ratio2DInvXSectionEtaa = new TH2F("ratio2DInvXSectionEtaa","ratio2DInvXSectionEtaa",1000,0.2,25.,1000,0.41,2.29);

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

  DrawGammaSetMarkerTGraphAsym(graphRatioCombEtaFitSta,20,1., kBlue, kBlue,1.5);  
  graphRatioCombEtaFitSta->Draw("p,same");

  TLatex * lt3a = new TLatex(2.8,2.,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV") ;
  lt3a->SetTextColor(kBlack) ;
  lt3a->SetTextSize(0.05) ;
  if (IsNSD)  lt3a->DrawLatex(2.8,1.85,"NSD #eta, ALICE");
  else  lt3a->DrawLatex(2.8,1.85,"#eta, ALICE");
  //		lt3->DrawLatex(3.,1.8,"#pi^{0} #rightarrow #gamma#gamma");
  //	lt3->DrawLatex(3.,1.8,"#pi^{0}");
  //	lt3->DrawLatex(2.8,1.8,"#pi^{0} #rightarrow #gamma#gamma, #pi^{0} #rightarrow e^{+}e^{-}#gamma");
  lt3a->Draw("same") ;

 

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
  ratio2DInvXSectionEtab = new TH2F("ratio2DInvXSectionEtab","ratio2DInvXSectionEtab",1000,0.2,25.,1000,0.41,2.29);

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
  DrawGammaSetMarker(histoRatioCombEtaPCMStat, 20,1, 1 , 1);
  histoRatioCombEtaPCMStat->Draw("E1,same") ;
 
  DrawGammaSetMarkerTGraphAsym(graphRatioCombEtaEMCalSys,33,1,kGreen+2 , kGreen+2, 1, kTRUE);  
  graphRatioCombEtaEMCalSys->Draw("E2same");
  DrawGammaSetMarker(histoRatioCombEtaEMCalStat, 33,1.5 ,kGreen+2 ,kGreen+2);
  histoRatioCombEtaEMCalStat->Draw("E1,same") ;
	
  TLatex * lt4a = new TLatex(2.8,2.1,"ALICE Preliminary") ;
  lt4a->SetTextColor(kBlack) ;
  lt4a->SetTextSize(0.05) ; 
  if (IsNSD)lt4a->DrawLatex(2.8,1.95,"p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  if (IsNSD)lt4a->DrawLatex(2.8,1.8,"NSD #eta");
  else  lt4a->DrawLatex(2.8,1.95,"#eta, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  //	lt4->DrawLatex(3.,1.8,"#pi^{0} #rightarrow #gamma#gamma");
  //	lt4->DrawLatex(2.8,1.8,"#pi^{0} #rightarrow #gamma#gamma, #pi^{0} #rightarrow e^{+}e^{-}#gamma");
  lt4a->Draw("same");

  TLegend* legendSpectraDiffDetMinBiasStripEta = new TLegend(0.18,0.75,0.4,0.9);
  legendSpectraDiffDetMinBiasStripEta->SetFillColor(0);
  legendSpectraDiffDetMinBiasStripEta->SetLineColor(0);
  legendSpectraDiffDetMinBiasStripEta->SetTextFont(42);
  legendSpectraDiffDetMinBiasStripEta->AddEntry(histoRatioCombEtaPCMStat,Form("PCM"),"pf");
  legendSpectraDiffDetMinBiasStripEta->AddEntry(histoRatioCombEtaEMCalStat,Form("EMCal"),"pf");
  legendSpectraDiffDetMinBiasStripEta->Draw("same");

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

  c1->Update();   
  c1->Print(Form("%s/XShiftTest.%s",outputDir.Data(),suffix.Data()));


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
  if (IsNSD){
    legendRpPbCombineXShiftEta->AddEntry(graphCombEtaInvCrossSectionTotpPb5023GeV,"NSD #eta, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
    legendRpPbCombineXShiftEta->AddEntry(graphCombEtaInvCrossSectionTotpPb5023GeVXShifted,"NSD #eta x bin shifted, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  }else{
    legendRpPbCombineXShiftEta->AddEntry(graphCombEtaInvCrossSectionTotpPb5023GeV,"#eta, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
    legendRpPbCombineXShiftEta->AddEntry(graphCombEtaInvCrossSectionTotpPb5023GeVXShifted,"#eta x bin shifted, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  }
  if (FittingType.CompareTo("Tsallis")==0) legendRpPbCombineXShiftEta->AddEntry(fitCombEtapPb5023GeVPtXShift,"Tsallis Fit","l");
  else legendRpPbCombineXShiftEta->AddEntry(fitCombEtapPb5023GeVPtXShift,"Bylinkin-Rostovtsev Fit","l");
    
  legendRpPbCombineXShiftEta->Draw();
	


  canvasXShiftEta->Update();
  canvasXShiftEta->Print(Form("%s/CombinedEtaSpectrum_XShifted.%s",outputDir.Data(),suffix.Data()));

  Double_t* XEtaUnshifted=graphCombEtaInvCrossSectionTotpPb5023GeV->GetX();
  Double_t* XEtaShifted=graphCombEtaInvCrossSectionTotpPb5023GeVXShifted->GetX();

  for (Int_t i=0;i<graphCombEtaInvCrossSectionTotpPb5023GeV->GetN();i++){

    fFits <<  XEtaUnshifted[i]-graphCombEtaInvCrossSectionTotpPb5023GeV->GetErrorXlow(i)<< " - " <<  XEtaUnshifted[i]+ graphCombEtaInvCrossSectionTotpPb5023GeV->GetErrorXhigh(i)<< " GeV/$c$ & " << XEtaUnshifted[i] <<" GeV/$c$ & " << XEtaShifted[i]<<" GeV/$c$ & " << (XEtaShifted[i]-XEtaUnshifted[i])*1000 << " MeV/$c$"<< endl;
  }







  // **************************************************************************************
  // ************************* Plotting x shifted Eta spectrum + ratio ************************
  // **************************************************************************************

  TGraphAsymmErrors*  graphEtaRatioCombXshiftFitSta    = CalculateGraphErrRatioToFit(graphCombEtaInvCrossSectionStatpPb5023GeVXShifted,fitCombEtapPb5023GeVPt); 
  TGraphAsymmErrors*   graphEtaRatioCombXShiftFitSys    = CalculateGraphErrRatioToFit(graphCombEtaInvCrossSectionSyspPb5023GeVXShifted,fitCombEtapPb5023GeVPt); 



  TCanvas* CanvasRatioEta = new TCanvas("CanvasRatioEta","",200,10,1200,1200);  // gives the page size
  DrawGammaCanvasSettings(CanvasRatioEta , 0.3, 0.02, 0.02, 0.16);
  TPad* padHistosEta = new TPad("padHistosEta", "", 0., 0.25, 1., 1.,-1, -1, -2);
  DrawGammaPadSettings( padHistosEta, 0.15, 0.02, 0.02, 0.);
  padHistosEta->Draw();

  TPad* padRatiosEta = new TPad("padRatiosEta", "", 0., 0., 1., 0.25,-1, -1, -2);
  DrawGammaPadSettings( padRatiosEta, 0.15, 0.02, 0., 0.2);
  padRatiosEta->Draw();

  padHistosEta->cd(); 
  padHistosEta->SetLogx();
  padHistosEta->SetLogy();
  TH2F * histoRatioEta;
  histoRatioEta = new TH2F("histoRatioEta","histoRatioEta",1000,0.2,30.,1000,1.2e-9,.9);
  SetStyleHistoTH2ForGraphs(histoRatioEta,"#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.04,0.05, 0.045,0.045, 0.7,1.4, 512, 505);
  histoRatioEta->DrawCopy();
  fitCombEtapPb5023GeVPt->Draw("same");
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionSyspPb5023GeVXShifted,20,0.7, 4, 4, 1, kTRUE);  
      
  graphCombEtaInvCrossSectionSyspPb5023GeVXShifted->Draw("E2,same");

  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionStatpPb5023GeVXShifted,20,0.7, kBlue, kBlue);  
  graphCombEtaInvCrossSectionStatpPb5023GeVXShifted->Draw("p,same");

  TLatex * ltPre2 = new TLatex(.35,3.7e-7,"ALICE Preliminary"); 
  ltPre2->SetTextSize(0.05) ;
  //  lt->DrawText(1.,1.45,"2016/06/14") ;
  ltPre2->Draw("same"); 
  
  TLegend* legendRpPbCombineRatioEta = new TLegend(0.23,0.10,0.7,0.250);
  legendRpPbCombineRatioEta->SetFillColor(0);
  legendRpPbCombineRatioEta->SetLineColor(0);
  legendRpPbCombineRatioEta->SetTextSize(0.04);
  if (IsNSD)  legendRpPbCombineRatioEta->AddEntry(graphInvYieldEtaCombpPb5023GeVSysClone,"NSD #eta, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  else  legendRpPbCombineRatioEta->AddEntry(graphInvYieldEtaCombpPb5023GeVSysClone,"#eta, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  if(FittingType.CompareTo("Tsallis")==0) legendRpPbCombineRatioEta->AddEntry(fitCombEtapPb5023GeVPt,"Tsallis Fit","l");
  else legendRpPbCombineRatioEta->AddEntry(fitCombEtapPb5023GeVPt,"Bylinkin-Rostovtsev Fit","l");
    
  legendRpPbCombineRatioEta->Draw();
  
  padRatiosEta->cd();
    
  padRatiosEta->SetLogx();
  //  canvasRatioAllppreferences->SetLogy();
  TH2F * histo2DRatioAllppreferences2Eta;
  histo2DRatioAllppreferences2Eta = new TH2F("histo2DRatioAllppreferences2Eta","histo2DRatioAllppreferencesEta",1000,.2,30.,1000,0.41,1.89);
  SetStyleHistoTH2ForGraphs(histo2DRatioAllppreferences2Eta, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.13,0.13, 0.13,0.13, 0.5,.5, 502, 503);
  histo2DRatioAllppreferences2Eta->DrawCopy();

  DrawGammaSetMarkerTGraphAsym(graphEtaRatioCombXShiftFitSys,20,0.7, 4, 4, 1, kTRUE);  
      
  graphEtaRatioCombXShiftFitSys->Draw("E2,same");

  DrawGammaSetMarkerTGraphAsym(graphEtaRatioCombXshiftFitSta,20,0.7, kBlue, kBlue);  
  graphEtaRatioCombXshiftFitSta->Draw("p,same");  

 
  TLine *lineA2Eta=new TLine(0.,1.,25.,1.);
  lineA2Eta->SetLineColor(kGray+1);
  lineA2Eta->Draw();   

    
  CanvasRatioEta->Print(Form("%s/EtaxShift_ratio.%s",outputDir.Data(),suffix.Data()));




  //*******************************yShift for eta/pi0*******************************************************

  // //PCM:	

  TGraphAsymmErrors* graphPCMYieldEtapPbYShift=(TGraphAsymmErrors*)graphPCMYieldEtapPbStat->Clone("graphPCMYieldEtapPbYShift");
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
  
  
  //PCM-EMCal:
  cout <<"!!!!!!!!!!!!!!!!!!!PCM-EMCal!!!!!!!!!!!!!!!!!!!!"<< endl;

  
  
  
  
  

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
  DrawGammaSetMarkerTGraphAsym(graphPCMYieldEtapPbStat,21,0.7, 4, 4, 1, kTRUE);  
  graphPCMYieldEtapPbStat->Draw("pe,same");
 
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
  legendRpPbCombineYShiftEtaPCM->AddEntry(graphPCMYieldEtapPbStat,"PCM #eta","pef");
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
    
  legendRpPbCombineYShiftEtaEMCal->Draw("same");
	


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


  for (Int_t j=0;j<graphPCMEtaPi0RatioSystErrYShift->GetN();j++){
    graphPCMEtaPi0RatioSystErrYShift->SetPoint(j,PCMX[j],PCMY[j]);
  }
  for (Int_t k=0;k<graphEMCalEtaPi0RatiopPbSystErr->GetN();k++){
    graphEMCalEtaPi0RatioSystErrYShift->SetPoint(k,EMCalX[k],EMCalY[k]);
  }
 

  TH1D *histoPCMEtaPi0RatioStatErrYShift=GraphAsymErrorsToHist_withErrors(graphPCMEtaPi0RatioStatErrYShift,"histoPCMEtaPi0RatiopPbStatErrYShift");
  TH1D *histoEMCalEtaPi0RatioStatErrYShift=GraphAsymErrorsToHist_withErrors(graphEMCalEtaPi0RatioStatErrYShift,"histoEMCalEtaPi0RatiopPbStatErrYShift");


  TString fileNameOutputWeightingEtaPi0Ratio                 = Form("%s/WeightingEtaPi0Ratio.dat",outputDir.Data());
  
  

  statErrorCollectionEtaPi0Ratio[0]          = (TH1D*)histoPCMEtaPi0RatioStatErrYShift->Clone("statErrPCMEta");
  statErrorCollectionEtaPi0Ratio[2]          = (TH1D*)histoEMCalEtaPi0RatioStatErrYShift->Clone("statErrEMCalEta");
  statErrorCollectionEtaPi0Ratio[4]          = (TH1D*)histoPCMEMCalEtaPi0RatioStatErrYShift->Clone("statErrPCMEMCalEta");
    
  sysErrorCollectionEtaPi0Ratio[0]           = (TGraphAsymmErrors*)graphPCMEtaPi0RatioSystErrYShift->Clone("sysErrPCMEta");
  sysErrorCollectionEtaPi0Ratio[2]           = (TGraphAsymmErrors*)graphEMCalEtaPi0RatioSystErrYShift->Clone("sysErrEMCalEta");
  sysErrorCollectionEtaPi0Ratio[4]           = (TGraphAsymmErrors*)graphPCMEMCalEtaPi0RatioSystErrYShift->Clone("sysErrPCMEMCalEta");
    
    
  TGraphAsymmErrors* graphCombEtaPi0RatioStatpPb5023GeV= NULL;
  TGraphAsymmErrors* graphCombEtaPi0RatioSyspPb5023GeV = NULL;
  
  cout<<"\n\n */////////////////Combination of EtatoPi0 ratio  begins here**//////////////"<<endl;
    
  TGraphAsymmErrors* graphCombEtaPi0RatioTotpPb5023GeV = CombinePtPointsSpectraFullCorrMat(    statErrorCollectionEtaPi0Ratio,    sysErrorCollectionEtaPi0Ratio,     
											       pTLimitsEtaPi0Ratio, NtotalEtaPi0Ratio,
											       offSetsEtaPi0Ratio, offSetsSysEtaPi0Ratio,
											       graphCombEtaPi0RatioStatpPb5023GeV,graphCombEtaPi0RatioSyspPb5023GeV ,
											       fileNameOutputWeightingEtaPi0Ratio,"2.76TeV", "EtaToPi0", kFALSE,
											       NULL, corrFactorsFile.Data()
											       );





  cout<<"\n\n */////////////////Combination of EtatoPi0 ratio  ends here**//////////////"<<endl;
  
                                   
  TCanvas* canvasYShiftEtaPi0Ratio = new TCanvas("canvasYShiftEtaPi0Ratio","",200,10,500,400);  // gives the page size
  DrawGammaCanvasSettings( canvasYShiftEtaPi0Ratio, 0.11, 0.02, 0.02, 0.12);

  TH2F * histo2DCompCombinedYShiftEtaPi0Ratio;
  histo2DCompCombinedYShiftEtaPi0Ratio = new TH2F("histo2DCompCombinedYShiftEtaPi0Ratio","histo2DCompCombinedYShiftEtaPi0Ratio",1000,0.3,16.,1000,0.,1.   );
  SetStyleHistoTH2ForGraphs(histo2DCompCombinedYShiftEtaPi0Ratio, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0} ", 0.05,0.05, 0.05,0.06, 1.0,.8, 512, 508);
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
  legendRpPbCombineYShiftEtaPi0Ratio->SetTextSize(0.04);
  legendRpPbCombineYShiftEtaPi0Ratio->AddEntry(graphCombEtaPi0RatioSyspPb5023GeV,"#eta/#pi^{0} combined","pef");
  legendRpPbCombineYShiftEtaPi0Ratio->AddEntry(graphPCMEtaPi0RatioSystErrYShift,"#eta/#pi^{0} PCM","pef");
  legendRpPbCombineYShiftEtaPi0Ratio->AddEntry(graphEMCalEtaPi0RatioSystErrYShift,"#eta/#pi^{0} EMCal","pef");
 

  
  legendRpPbCombineYShiftEtaPi0Ratio->Draw("same");
	
  TLatex * lt4b = new TLatex(9.,0.9,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV") ;
  lt4b->SetTextColor(kBlack) ;
  lt4b->SetTextSize(0.05) ;
  lt4b->Draw("same") ;

  canvasYShiftEtaPi0Ratio->Update();
  canvasYShiftEtaPi0Ratio->Print(Form("%s/EtaPi0Ratio_YShifted.%s",outputDir.Data(),suffix.Data()));

                                   
  TCanvas* canvasYShiftEtaPi0RatioPrelim = new TCanvas("canvasYShiftEtaPi0RatioPrelim","",200,10,500,400);  // gives the page size
  DrawGammaCanvasSettings( canvasYShiftEtaPi0RatioPrelim, 0.11, 0.02, 0.02, 0.12);

  TH2F * histo2DCompCombinedYShiftEtaPi0RatioPrelim;
  histo2DCompCombinedYShiftEtaPi0RatioPrelim = new TH2F("histo2DCompCombinedYShiftEtaPi0RatioPrelim","histo2DCompCombinedYShiftEtaPi0RatioPrelim",1000,0.3,16.,1000,0.,1.   );
  SetStyleHistoTH2ForGraphs(histo2DCompCombinedYShiftEtaPi0RatioPrelim, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0} ", 0.05,0.05, 0.05,0.06, 1.,.8, 512, 508);
  histo2DCompCombinedYShiftEtaPi0RatioPrelim->DrawCopy(); 
      
  
  DrawGammaSetMarkerTGraphAsym(graphCombEtaPi0RatioStatpPb5023GeV,20,0.7,4, 4);  
  graphCombEtaPi0RatioStatpPb5023GeV->Draw("pe,same");
  DrawGammaSetMarkerTGraphAsym(graphCombEtaPi0RatioSyspPb5023GeV,20,0.7, 4, 4, 1, kTRUE);  
  graphCombEtaPi0RatioSyspPb5023GeV->Draw("E2,same");
  
 
     
  TLegend* legendRpPbCombineYShiftEtaPi0RatioPrelim = new TLegend(0.15,0.75,0.5,0.90);
  legendRpPbCombineYShiftEtaPi0RatioPrelim->SetFillColor(0);
  legendRpPbCombineYShiftEtaPi0RatioPrelim->SetLineColor(0);
  legendRpPbCombineYShiftEtaPi0RatioPrelim->SetTextSize(0.03);
  legendRpPbCombineYShiftEtaPi0RatioPrelim->AddEntry(graphCombEtaPi0RatioSyspPb5023GeV,"#eta/#pi^{0} ","pef");

  
  // legendRpPbCombineYShiftEtaPi0RatioPrelim->Draw("same");
	
  TLatex * lt4Pre = new TLatex(9.,0.9,"ALICE Preliminary") ;
  lt4Pre->SetTextColor(kBlack) ;
  lt4Pre->SetTextSize(0.05) ; 
 lt4Pre->DrawLatex(9.,0.8,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  lt4Pre->Draw("same") ;

  canvasYShiftEtaPi0RatioPrelim->Update();
  canvasYShiftEtaPi0RatioPrelim->Print(Form("%s/EtaPi0RatioPrelim_YShifted.%s",outputDir.Data(),suffix.Data()));	
  //*********************** Comparison of Measurements in ALICE ***************************************************
  Double_t fBinsPi07TeVPt[22] =					{0.0, 0.1, 0.2, 0.3, 0.4,
								 0.5, 0.6, 0.8, 1.0, 1.5,  
								 1.9, 2.2, 2.7, 3.2,
								 4.0,  5.0, 6.0, 8.0,
								 10.0, 12.0, 18.0, 22.0};

  TFile* cocktailFileMtScaledEta = new TFile("CocktailInput/cocktail_PbPb_0020c_dmtsallis_MTEta.root");
  TDirectory* cocktailDirMtScaledEta = (TDirectoryFile*) cocktailFileMtScaledEta->Get("cocktail_PbPb_0020c_dmtsallis_MTEta");

  TH1D* cocktailPi0_MtScaled = (TH1D* )cocktailDirMtScaledEta->Get("ptPi0");
  cocktailPi0_MtScaled->Sumw2();
  TH1D* cocktailPi0_MtScaledRebinned = (TH1D* )cocktailPi0_MtScaled->Rebin(21,"ptPionMTScaledRebinned",fBinsPi07TeVPt);
  TH1D* cocktailEta_MtScaled = (TH1D* )cocktailDirMtScaledEta->Get("ptEta");
  cocktailEta_MtScaled->Sumw2();
  TH1D* cocktailEta_MtScaledRebinned = (TH1D* )cocktailEta_MtScaled->Rebin(21,"ptEtaMTScaledRebinned",fBinsPi07TeVPt);
  // TH1D* cocktailOmega_MtScaled = (TH1D* )cocktailDirMtScaledEta->Get("ptOmega");
  // cocktailOmega_MtScaled->Sumw2();
  // TH1D* cocktailOmega_MtScaledRebinned = (TH1D* )cocktailOmega_MtScaled->Rebin(21,"ptOmegaMTScaledRebinned",fBinsPi07TeVPt);

  TH1D* cocktailEtaToPi0Ratio_MtScaled = (TH1D* )cocktailEta_MtScaled->Clone("EtaToPi0Ratio_MtScaled");
  cocktailEtaToPi0Ratio_MtScaled->Sumw2();
  cocktailEtaToPi0Ratio_MtScaled->Divide(cocktailEtaToPi0Ratio_MtScaled,cocktailPi0_MtScaled);
  TH1D* cocktailEtaToPi0Ratio_MtScaledRebinned = (TH1D* )cocktailEta_MtScaledRebinned->Clone("EtaToPi0Ratio_MtScaledRebinned");
  cocktailEtaToPi0Ratio_MtScaledRebinned->Sumw2();
  cocktailEtaToPi0Ratio_MtScaledRebinned->Divide(cocktailEtaToPi0Ratio_MtScaledRebinned,cocktailPi0_MtScaledRebinned);

  TFile* fileCombinedpp 				= new TFile("FinalResults/CombinedResultsPP_ShiftedX_PaperRAA_16_May_2014.root ");
  // 	TFile* fileCombinedpp 				= new TFile("CombinedResultsPaperX_18_Feb_2014.root");
  TGraphAsymmErrors* graphCombEtaToPi0Ratiopp7TeV =         (TGraphAsymmErrors*)fileCombinedpp->Get("graphEtaToPi0Comb7TeVStat");
  TGraphAsymmErrors* graphCombEtaToPi0RatioSysErrpp7TeV=    (TGraphAsymmErrors*)fileCombinedpp->Get("graphEtaToPi0Comb7TeVSys"); 
                                 
  TCanvas* canvasYShiftEtaPi0RatioALICE = new TCanvas("canvasYShiftEtaPi0RatioALICE","",200,10,500,400);  // gives the page size
  DrawGammaCanvasSettings( canvasYShiftEtaPi0RatioALICE, 0.11, 0.02, 0.02, 0.12);

  TH2F * histo2DCompCombinedYShiftEtaPi0RatioALICE;
  histo2DCompCombinedYShiftEtaPi0RatioALICE = new TH2F("histo2DCompCombinedYShiftEtaPi0RatioALICE","histo2DCompCombinedYShiftEtaPi0RatioALICE",1000,0.3,16.,1000,0.,1.   );
  SetStyleHistoTH2ForGraphs(histo2DCompCombinedYShiftEtaPi0RatioALICE, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0} ", 0.05,0.05, 0.05,0.06, 1.0,.8, 512, 508);
  histo2DCompCombinedYShiftEtaPi0RatioALICE->DrawCopy(); 
      
  cocktailEtaToPi0Ratio_MtScaledRebinned->GetXaxis()->SetRangeUser(0.3,16);
  DrawGammaSetMarker(cocktailEtaToPi0Ratio_MtScaledRebinned, 2, 0, kRed+2, kRed+2);
  cocktailEtaToPi0Ratio_MtScaledRebinned->SetLineStyle(2);
  cocktailEtaToPi0Ratio_MtScaledRebinned->SetLineWidth(2);
  cocktailEtaToPi0Ratio_MtScaledRebinned->Draw("same,hist,c");

  DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Ratiopp7TeV, 21, 1, kGray+3, kGray+3);
  graphCombEtaToPi0Ratiopp7TeV->Draw("same,e1p");
  DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RatioSysErrpp7TeV, 21, 1, kGray+3, kGray+3, 1., kTRUE);
  graphCombEtaToPi0RatioSysErrpp7TeV->Draw("same,E2");
   
  DrawGammaSetMarkerTGraphAsym(graphCombEtaPi0RatioStatpPb5023GeV,24,1,4, 4);  
  graphCombEtaPi0RatioStatpPb5023GeV->Draw("pe,same");
  DrawGammaSetMarkerTGraphAsym(graphCombEtaPi0RatioSyspPb5023GeV,24,1, 4, 4, 1, kTRUE);  
  graphCombEtaPi0RatioSyspPb5023GeV->Draw("E2,same");
     
  //  TLegend* legendRpPbCombineYShiftEtaPi0RatioALICE = new TLegend(0.15,0.75,0.5,0.90);
  TLegend* legendRpPbCombineYShiftEtaPi0RatioALICE = new TLegend(0.15,0.78,0.5,0.93);
  legendRpPbCombineYShiftEtaPi0RatioALICE->SetFillColor(0);
  legendRpPbCombineYShiftEtaPi0RatioALICE->SetLineColor(0);
  legendRpPbCombineYShiftEtaPi0RatioALICE->SetTextSize(0.035);
  if(IsNSD)  legendRpPbCombineYShiftEtaPi0RatioALICE->AddEntry(graphCombEtaPi0RatioSyspPb5023GeV,"NSD #eta/#pi^{0} combined PCM+EMCal ","pef");
  else  legendRpPbCombineYShiftEtaPi0RatioALICE->AddEntry(graphCombEtaPi0RatioSyspPb5023GeV,"#eta/#pi^{0} combined PCM+EMCal ","pef");
  legendRpPbCombineYShiftEtaPi0RatioALICE->AddEntry(graphCombEtaToPi0RatioSysErrpp7TeV,"#eta/#pi^{0} PCM - pp 7 TeV","pef");
  legendRpPbCombineYShiftEtaPi0RatioALICE->AddEntry(cocktailEtaToPi0Ratio_MtScaledRebinned,"#eta/#pi^{0}, #eta from #it{m}_{T} scaled #pi^{0}","pl");

  
  legendRpPbCombineYShiftEtaPi0RatioALICE->Draw("same");
	
  //  TLatex * lt4c = new TLatex(9.,0.9,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV") ;
  TLatex * lt4c = new TLatex(7.,0.15,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV") ;
  lt4c->SetTextColor(kBlack) ;
  lt4c->SetTextSize(0.05) ;
  lt4c->Draw("same") ;

  canvasYShiftEtaPi0RatioALICE->Update();
  canvasYShiftEtaPi0RatioALICE->Print(Form("%s/EtaPi0Ratio_ALICE.%s",outputDir.Data(),suffix.Data()));


   //*********************** Comparison of Measurements in ALICE prelim  ***************************************************

  TCanvas* canvasYShiftEtaPi0RatioALICEPrelim = new TCanvas("canvasYShiftEtaPi0RatioALICEPrelim","",200,10,500,400);  // gives the page size
  DrawGammaCanvasSettings( canvasYShiftEtaPi0RatioALICEPrelim, 0.11, 0.02, 0.02, 0.12);

  TH2F * histo2DCompCombinedYShiftEtaPi0RatioALICEPrelim;
  histo2DCompCombinedYShiftEtaPi0RatioALICEPrelim = new TH2F("histo2DCompCombinedYShiftEtaPi0RatioALICEPrelim","histo2DCompCombinedYShiftEtaPi0RatioALICEPrelim",1000,0.3,16.,1000,0.,1.   );
  SetStyleHistoTH2ForGraphs(histo2DCompCombinedYShiftEtaPi0RatioALICEPrelim, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0} ", 0.05,0.05, 0.05,0.06, 1.0,.8, 512, 508);
  histo2DCompCombinedYShiftEtaPi0RatioALICEPrelim->DrawCopy(); 
      
  EtaPi0Ratio_mTScaled_pPb->Draw("same"); //NOTE all
  graphEtaPi0Ratio_mTScaled_pPb->SetMarkerColor(kRed+2);
  graphEtaPi0Ratio_mTScaled_pPb->SetLineColor(1);
  graphEtaPi0Ratio_mTScaled_pPb->SetFillColor(kRed+2);
  
  
  DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Ratiopp7TeV, 21, 1, kGray+3, kGray+3);
  graphCombEtaToPi0Ratiopp7TeV->Draw("same,e1p");
  DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RatioSysErrpp7TeV, 21, 1, kGray+3, kGray+3, 1., kTRUE);
  graphCombEtaToPi0RatioSysErrpp7TeV->Draw("same,E2");
   
  DrawGammaSetMarkerTGraphAsym(graphCombEtaPi0RatioStatpPb5023GeV,24,1,4, 4);  
  graphCombEtaPi0RatioStatpPb5023GeV->Draw("pe,same");
  DrawGammaSetMarkerTGraphAsym(graphCombEtaPi0RatioSyspPb5023GeV,24,1, 4, 4, 1, kTRUE);  
  graphCombEtaPi0RatioSyspPb5023GeV->Draw("E2,same");
  
  
  TLegend* legendRpPbCombineYShiftEtaPi0RatioALICEPrelim = new TLegend(0.22,0.15,0.6,0.38);
  legendRpPbCombineYShiftEtaPi0RatioALICEPrelim->SetFillColor(0);
  legendRpPbCombineYShiftEtaPi0RatioALICEPrelim->SetLineColor(0);
  legendRpPbCombineYShiftEtaPi0RatioALICEPrelim->SetTextSize(0.04);
  legendRpPbCombineYShiftEtaPi0RatioALICEPrelim->AddEntry(graphCombEtaPi0RatioSyspPb5023GeV,"#eta/#pi^{0}, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  legendRpPbCombineYShiftEtaPi0RatioALICEPrelim->AddEntry(graphCombEtaToPi0RatioSysErrpp7TeV,"#eta/#pi^{0}, pp #sqrt{#it{s}} = 7 TeV (PLB717 (2012) 162-172)","pef");
  legendRpPbCombineYShiftEtaPi0RatioALICEPrelim->AddEntry(EtaPi0Ratio_mTScaled_pPb,"#eta/#pi^{0}, eta from #it{m}_{T} scaled #pi^{0}","pl"); //NOTE
   
  legendRpPbCombineYShiftEtaPi0RatioALICEPrelim->Draw("same");

  TLatex * lt4d1 = new TLatex(2.,0.88,"ALICE Preliminary") ;
  lt4d1->SetTextColor(kBlack) ;
  lt4d1->SetTextSize(0.05) ;
  lt4d1->Draw("same") ;

  canvasYShiftEtaPi0RatioALICEPrelim->Update();
  canvasYShiftEtaPi0RatioALICEPrelim->Print(Form("%s/EtaPi0Ratio_ALICEPrelim.%s",outputDir.Data(),suffix.Data()));
 


  // **************************************************************************************
  // ************************* Comparison Charged Pions            ************************
  // **************************************************************************************


  TGraphAsymmErrors* graphCombPi0InvCrossSectionSyspPb5023GeVYShiftedCopy = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionSyspPb5023GeVYShifted->Clone("graphCombPi0InvCrossSectionSyspPb5023GeVYShiftedCopy");
  TGraphAsymmErrors* graphCombPi0InvCrossSectionStatpPb5023GeVYShiftedCopy = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionStatpPb5023GeVYShifted->Clone("graphCombPi0InvCrossSectionStatpPb5023GeVYShiftedCopy");
  TGraphAsymmErrors* graphChargedPionspPbSystCopy = (TGraphAsymmErrors*) graphChargedPionspPbSyst->Clone("graphChargedPionspPbSystCopy");
  TGraphAsymmErrors* graphChargedPionspPbStatCopy = (TGraphAsymmErrors*) graphChargedPionspPbStat->Clone("graphChargedPionspPbStatCopy");
  TH1D*  histoCombPi0InvCrossSectionStatpPb5023GeVYShifted=(TH1D *)GraphAsymErrorsToHist_withErrors(graphCombPi0InvCrossSectionStatpPb5023GeVYShifted, "histoCombPi0InvCrossSectionStatpPb5023GeVYShifted");

//   TGraphErrors* bla1 = NULL;
//   TGraphErrors* bla2 = NULL;
//   TGraphErrors* bla3 = NULL;
//   TGraphErrors* bla4 = NULL;
//   cout << "PCM Spectrum  - Dalitz" << endl;
//   TGraphErrors* graphRatioPi0ChargedPions = CalculateRatioBetweenSpectraWithDifferentBinning(histoCombPi0InvCrossSectionStatpPb5023GeVYShifted,graphCombPi0InvCrossSectionSyspPb5023GeVYShiftedCopy,histoChargedPionspPbStatErr,graphChargedPionspPbSystCopy ,  kTRUE,  kTRUE,&bla1,&bla2,&bla3,&bla4)  ;


  //NOTE splitting errors
  TGraphErrors* graphCombPi0InvCrossSectionStatpPb5023GeVYShiftedRebin = NULL;
  TGraphErrors* graphCombPi0InvCrossSectionSyspPb5023GeVYShiftedRebinOneMat = NULL;
  TGraphErrors* graphChargedPionspPbStatErrRebin = NULL;
  TGraphErrors* graphChargedPionspPbSystErrRebin = NULL;
  
  //TGraphErrors* graphCombPi0InvCrossSectionSyspPb5023GeVYShiftedRebin = NULL;
  cout << "PCM Spectrum  - Dalitz" << endl;
  TGraphErrors* graphRatioPi0ChargedPions = CalculateRatioBetweenSpectraWithDifferentBinning(histoCombPi0InvCrossSectionStatpPb5023GeVYShifted,graphCombPi0InvCrossSectionSyspPb5023GeVYShiftedCopy,histoChargedPionspPbStatErr,graphChargedPionspPbSystCopy ,  kTRUE,  kTRUE,&graphCombPi0InvCrossSectionStatpPb5023GeVYShiftedRebin,&graphCombPi0InvCrossSectionSyspPb5023GeVYShiftedRebinOneMat,&graphChargedPionspPbStatErrRebin,&graphChargedPionspPbSystErrRebin)  ;
  
  
   TGraphErrors* graphCombPi0InvCrossSectionSyspPb5023GeVYShiftedRebin         = (TGraphErrors*)graphCombPi0InvCrossSectionSyspPb5023GeVYShiftedRebinOneMat->Clone("Dummy");
   Double_t * yValue           = graphCombPi0InvCrossSectionSyspPb5023GeVYShiftedRebin->GetY();
   Double_t* yError            = graphCombPi0InvCrossSectionSyspPb5023GeVYShiftedRebin->GetEY();
   Int_t nPoints               = graphCombPi0InvCrossSectionSyspPb5023GeVYShiftedRebin->GetN();
  
  
  for(Int_t i=0; i< nPoints; i++){
    
   
    yError[i] = TMath::Sqrt( TMath::Power((yError[i]/yValue[i]),2) + 0.045*0.045)*yValue[i];
    
	  
    
  }
  
  
  

  TGraphErrors* graphRatioPi0ChargedPionsStatErr = CalculateGraphRatioToGraph( graphCombPi0InvCrossSectionStatpPb5023GeVYShiftedRebin,graphChargedPionspPbStatErrRebin);
  TGraphErrors* graphRatioPi0ChargedPionsSystErr = CalculateGraphRatioToGraph( graphCombPi0InvCrossSectionSyspPb5023GeVYShiftedRebin,graphChargedPionspPbSystErrRebin);

  
  
  
  
  TGraphErrors*  tempTest  = CalculateCombinedSysAndStatError(graphRatioPi0ChargedPionsStatErr,graphRatioPi0ChargedPionsSystErr);

  
  
  
   
  TFile* fileChargedPionsRatioPythia = 0x0;
  TH1F* hPionRatio = 0x0;
  
  fileChargedPionsRatioPythia = new TFile(fileNameChargedPionsRatioPythia.Data(),"READ");
  
  if( fileChargedPionsRatioPythia ){
      hPionRatio = (TH1F*)fileChargedPionsRatioPythia->Get("hPionRatio");
  }
  
  
  
  
  
  
  


  TCanvas* CanvasRatioPi0ChargedPions = new TCanvas("CanvasRatioPi0ChargedPions","",200,10,1200,1200);  // gives the page size
  DrawGammaCanvasSettings(CanvasRatioPi0ChargedPions , 0.3, 0.02, 0.02, 0.16);
  TPad* padHistos1 = new TPad("padHistos1", "", 0., 0.25, 1., 1.,-1, -1, -2);
  DrawGammaPadSettings( padHistos1, 0.15, 0.02, 0.02, 0.);
  padHistos1->Draw();

  TPad* padRatios1 = new TPad("padRatios1", "", 0., 0., 1., 0.25,-1, -1, -2);
  DrawGammaPadSettings( padRatios1, 0.15, 0.02, 0., 0.2);
  padRatios1->Draw();

  padHistos1->cd(); 
  padHistos1->SetLogx();
  padHistos1->SetLogy();
  TH2F * histoPi0ChargedPions;
  histoPi0ChargedPions = new TH2F("histoPi0ChargedPions","histoPi0ChargedPions",1000,0.2,30.,1000,1.2e-9,30.);
  SetStyleHistoTH2ForGraphs(histoPi0ChargedPions,"#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.04,0.05, 0.045,0.045, 0.7,1.4, 512, 505);
  histoPi0ChargedPions->DrawCopy();


  
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionSyspPb5023GeVYShifted,20,0.7, 4, 4, 1, kTRUE);  
  graphCombPi0InvCrossSectionSyspPb5023GeVYShifted->Draw("E2,same");
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionStatpPb5023GeVYShifted,20,0.7, kBlue, kBlue);  
  graphCombPi0InvCrossSectionStatpPb5023GeVYShifted->Draw("p,same");

  DrawGammaSetMarkerTGraphAsym(graphChargedPionspPbSyst,20,0.7,kRed+1 ,kRed+1 , 1, kTRUE);  
  graphChargedPionspPbSyst->Draw("E2,same");
  DrawGammaSetMarkerTGraphAsym(graphChargedPionspPbStat,20,0.7, kRed+1, kRed+1);  
  graphChargedPionspPbStat->Draw("p,same");
     
  TLegend* legendPi0ChargedPions = new TLegend(0.23,0.10,0.7,0.30);
  legendPi0ChargedPions->SetFillColor(0);
  legendPi0ChargedPions->SetLineColor(0); 
  legendPi0ChargedPions->SetTextSize(0.03);
  if (IsNSD)  legendPi0ChargedPions->AddEntry(graphCombPi0InvCrossSectionSyspPb5023GeVYShifted,"NSD #pi^{0}, |#it{y}_{lab}| < 0.8","pef");
  else  legendPi0ChargedPions->AddEntry(graphCombPi0InvCrossSectionSyspPb5023GeVYShifted,"#pi^{0}, |#it{y}_{lab}| < 0.8","pef");
  legendPi0ChargedPions->AddEntry(graphChargedPionspPbSyst,"NSD #pi^{+/-}, 0 < #it{y}_{cms} < 0.5","pef");
  TLatex * lt4d = new TLatex(3,1,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV") ;
  lt4d->SetTextColor(kBlack) ;
  lt4d->SetTextSize(0.04) ;
  lt4d->Draw("same");
  legendPi0ChargedPions->Draw();

  padRatios1->cd();
    
  padRatios1->SetLogx();
  
  DrawGammaSetMarkerTGraphErr(graphRatioPi0ChargedPionsSystErr,20,1.0, kBlue, kBlue, 1, kTRUE);
  DrawGammaSetMarkerTGraphErr(graphRatioPi0ChargedPionsStatErr,20,1.0, kBlue, kBlue);

  
  if( hPionRatio ) {
  DrawGammaSetMarker(hPionRatio,20,1.0, kRed, kRed);
  }


  //  canvasRatioAllppreferences->SetLogy();
  TH2F * histoRatioPi0ChargedPions;
  histoRatioPi0ChargedPions = new TH2F("histoRatioPi0ChargedPions","histo2DRatioAllppreferences",1000,.2,30.,1000,0.51,1.69);
  SetStyleHistoTH2ForGraphs(histoRatioPi0ChargedPions, "#it{p}_{T} (GeV/#it{c})","pi^{0}/pi^{+/-}", 0.13,0.13, 0.13,0.13, 0.5,.5, 502, 505);
 
 
  histoRatioPi0ChargedPions->DrawCopy();
     
  TLine *lineA3=new TLine(0.,1.,25.,1.);
  lineA3->SetLineColor(kGray+1);
  lineA3->Draw("same"); 
     
  //DrawGammaSetMarkerTGraph(graphRatioPi0ChargedPions,20,0.7, kBlue, kBlue);  
  //graphRatioPi0ChargedPions->Draw("p,same");  
  
  //tempTest->Draw("p,same");
  graphRatioPi0ChargedPionsStatErr->Draw("p,same,e");
  //graphRatioPi0ChargedPionsStatErr->SetFillColor(0);
  graphRatioPi0ChargedPionsSystErr->Draw("2same");
  
  if(hPionRatio)hPionRatio->Draw("p,same,e");
  
  
  
  TLegend* legendPi0ChargedPionsRatioPythia = new TLegend(0.2,0.22,0.5,0.4);
  legendPi0ChargedPionsRatioPythia->SetFillColor(0);
  legendPi0ChargedPionsRatioPythia->SetLineColor(0); 
  legendPi0ChargedPionsRatioPythia->SetTextSize(0.08);
  legendPi0ChargedPionsRatioPythia->SetNColumns(2);
  legendPi0ChargedPionsRatioPythia->AddEntry(graphRatioPi0ChargedPionsSystErr,"Data","efp");
  if (hPionRatio) legendPi0ChargedPionsRatioPythia->AddEntry(hPionRatio,"Pythia","efp");
  legendPi0ChargedPionsRatioPythia->Draw("same");
  
  
          
  CanvasRatioPi0ChargedPions->Print(Form("%s/Comparison_ChargedPions.%s",outputDir.Data(),suffix.Data()));


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
  
  fResults.Close();
  
  
  TFile fPaperPlots(Form("%s/PaperPlots_%s_%s.root",outputDir.Data(),FittingType.Data(), dateForOutput.Data()),"RECREATE");
  fitCombEtapPb5023GeVPt->Write("FitCombEta");
  fitCombPi0pPb5023GeVPt->Write("FitCombPi0");
  graphCombPi0InvCrossSectionSyspPb5023GeVXShifted->Write("CombPi0Syst");
  graphCombPi0InvCrossSectionStatpPb5023GeVXShifted->Write("CombPi0Stat");
  graphCombEtaInvCrossSectionSyspPb5023GeVXShifted->Write("CombEtaSyst");
  graphCombEtaInvCrossSectionStatpPb5023GeVXShifted->Write("CombEtaStat");

  graphRatioCombXShiftFitSys->Write("RatioTsallisCombPi0Syst");
  graphRatioCombXshiftFitSta->Write("RatioTsallisCombPi0Stat");  
  graphEtaRatioCombXShiftFitSys->Write("RatioTsallisCombEtaSyst");
  graphEtaRatioCombXshiftFitSta->Write("RatioTsallisCombEtaStat");

  graphRatioCombPCMSys->Write("RatioTsallisPCMSyst");
  histoRatioCombPCMStat->Write("RatioTsallisPCMStat") ;
  graphRatioCombDalitzSys->Write("RatioTsallisDalitzSyst");
  histoRatioCombDalitzStat->Write("RatioTsallisDalitzStat") ; 
  graphRatioCombEMCalSys->Write("RatioTsallisEMCalSyst");
  histoRatioCombEMCalStat->Write("RatioTsallisEMCalStat") ;
  graphRatioCombPHOSSys->Write("RatioTsallisPHOSSyst") ;
  histoRatioCombPHOSStat->Write("RatioTsallisPHOSStat");
  graphRatioCombPCMEMCalSys->Write( "RatioTsallisPCMEMCalSyst");
  histoRatioCombPCMEMCalStat->Write("RatioTsallisPCMEMCalStat");

  graphRatioCombEtaPCMSys->Write("RatioEtaTsallisPCMSyst");
  histoRatioCombEtaPCMStat->Write("RatioEtaTsallisPCMStat") ;
  graphRatioCombEtaEMCalSys->Write("RatioEtaTsallisEMCalSyst");
  histoRatioCombEtaEMCalStat->Write("RatioEtaTsallisEMCalStat") ;
  graphRatioCombEtaPCMEMCalSys->Write("RatioEtaTsallisPCMEMCalSyst");
  histoRatioCombEtaPCMEMCalStat->Write("RatioEtaTsallisPCMEMCalStat");
  
  EtaPi0Ratio_mTScaled_pPb->Write("mTScaling"); 
  graphCombEtaToPi0RatioSysErrpp7TeV->Write("EtaPi07TeVSyst");
  graphCombEtaToPi0Ratiopp7TeV->Write("EtaPi07TeVStat");
  graphCombEtaPi0RatioStatpPb5023GeV->Write("EtaPi0pPbStat");
  graphCombEtaPi0RatioSyspPb5023GeV->Write("EtaPi0pPbSyst");
  
  
  //graphEtaPi0Ratio_vsmT_Sys_pp7TeV->Write("EtaPi0Ratio_vsmT_Sys_pp7TeV");
  //graphEtaPi0Ratio_vsmT_Stat_pp7TeV->Write("EtaPi0Ratio_vsmT_Stat_pp7TeV");
  //graphEtaPi0Ratio_vsmT_Sys->Write("EtaPi0Ratio_vsmT_Sys"); NOTE
  //graphEtaPi0Ratio_vsmT_Stat->Write("EtaPi0Ratio_vsmT_Stat"); NOTE
  EtaSpectrum_mTScaled_pPb->Write("EtaSpectrum_mT_pPb"); 
  
  
  

  fPaperPlots.Close();

  
  
}
