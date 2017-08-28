/*****************************************************************************************************************************
 ******         provided by Gamma Conversion Group, PWG4,                                                                *****
 ******        Ana Marin, marin@physi.uni-heidelberg.de                                                                  *****
 ******        Annika Passfeld, annikapassfeld@uni-muenster.de                                                           *****
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

struct SysErrorConversion {
  Double_t value;
  Double_t error;
  //    TString name;
};

TGraphAsymmErrors* ApplyNSDSysError(TGraphAsymmErrors* graphMesonSystErr, Double_t ScalingErr, Double_t ScalingErrTpPb){
  TGraphAsymmErrors*   graphMesonSystErrClone=(TGraphAsymmErrors*) graphMesonSystErr->Clone();
  Int_t nPoints=graphMesonSystErrClone->GetN();
  Double_t* valueY    = graphMesonSystErrClone->GetY();
  Double_t* valueX    = graphMesonSystErrClone->GetX();
  Double_t* ErrYlowSys     = graphMesonSystErrClone->GetEYlow();
  Double_t* ErrYhighSys   = graphMesonSystErrClone->GetEYhigh();
  Double_t* ErrXlowSys     = graphMesonSystErrClone->GetEXlow();
  Double_t* ErrXhighSys   = graphMesonSystErrClone->GetEXhigh();

  for(Int_t i = 0; i < graphMesonSystErrClone->GetN(); i++){


    ErrYlowSys[i]  = TMath::Sqrt(   (ErrYlowSys[i] *ErrYlowSys[i]  ) + ( (valueY[i]*ScalingErr) *(valueY[i]*ScalingErr) ) + ( (valueY[i]*ScalingErrTpPb) *(valueY[i]*ScalingErrTpPb) ) );
    ErrYhighSys[i]  = TMath::Sqrt(   (ErrYhighSys[i] *ErrYhighSys[i]  ) + ( (valueY[i]*ScalingErr) *(valueY[i]*ScalingErr) ) + ( (valueY[i]*ScalingErrTpPb) *(valueY[i]*ScalingErrTpPb) ) );
        
  }
  TGraphAsymmErrors* graphMesonSystErrClone1= new TGraphAsymmErrors(nPoints,valueX,valueY,ErrXlowSys,ErrXhighSys,ErrYlowSys,ErrYhighSys);
      
  return graphMesonSystErrClone1;
}

TGraphAsymmErrors* ApplyTpPbSysError(TGraphAsymmErrors* graphMesonSystErr, Double_t ScalingErrTpPb){
  TGraphAsymmErrors*   graphMesonSystErrClone=(TGraphAsymmErrors*) graphMesonSystErr->Clone();
  Int_t nPoints=graphMesonSystErrClone->GetN();
  Double_t* valueY    = graphMesonSystErrClone->GetY();
  Double_t* valueX    = graphMesonSystErrClone->GetX();
  Double_t* ErrYlowSys     = graphMesonSystErrClone->GetEYlow();
  Double_t* ErrYhighSys   = graphMesonSystErrClone->GetEYhigh();
  Double_t* ErrXlowSys     = graphMesonSystErrClone->GetEXlow();
  Double_t* ErrXhighSys   = graphMesonSystErrClone->GetEXhigh();

  for(Int_t i = 0; i < graphMesonSystErrClone->GetN(); i++){


    ErrYlowSys[i]  = TMath::Sqrt(   (ErrYlowSys[i] *ErrYlowSys[i]  ) + ( (valueY[i]*ScalingErrTpPb) *(valueY[i]*ScalingErrTpPb) ) );
    ErrYhighSys[i]  = TMath::Sqrt(   (ErrYhighSys[i] *ErrYhighSys[i]  ) + ( (valueY[i]*ScalingErrTpPb) *(valueY[i]*ScalingErrTpPb) ) );
        
  }
  TGraphAsymmErrors* graphMesonSystErrClone1= new TGraphAsymmErrors(nPoints,valueX,valueY,ErrXlowSys,ErrXhighSys,ErrYlowSys,ErrYhighSys);
      
  return graphMesonSystErrClone1;
}
const Int_t  Ntotal = 26;//31
const Int_t  nPtLimits = Ntotal+1;

void CombineRpPb5023GeV(Bool_t IsNSD=kTRUE){    

  TString date                                        = ReturnDateString();
    
  gROOT->Reset();
  gROOT->SetStyle("Plain");
    
  StyleSettingsThesis();
  SetPlotStyle();
    
  TString suffix                                      = "eps";
    
  TString dateForOutput                               = ReturnDateStringForOutput();
  TString outputDir                                   = Form("CombineRpPb/%s/%s",suffix.Data(),dateForOutput.Data());
  if (IsNSD) outputDir                                = Form("CombineRpPb/%s/%s_NSD",suffix.Data(),dateForOutput.Data());

  gSystem->Exec("mkdir -p "+outputDir);
    
  Double_t ScalingErr = 0.031; //pure NSD Err
  Double_t Scaling = 0.964;// NSD Err
  Double_t ScalingErrTpPb =   0.0356; //Normalization Err on TpPb

  Double_t NormalizationError =TMath::Sqrt(pow(0.031,2)+pow(0.036,2)+pow(0.036,2));//pPb normalization,TpPb and pp normalization errors
  //cout << dateForOutput.Data() << endl;
  //___________________________________ Declaration of files _____________________________________________

    
 
  TString fileNameRpPb                        = "ExternalInputpPb/InputRpPb/Pi0RpPb_PCM_2017_08_25.root";
  TString fileNameRpPbCombined                = "ExternalInputpPb/InputRpPb/Pi0RpPb_Comb_2017_08_25.root";
  TString fileNameRpPbPCM                     = "ExternalInputpPb/InputRpPb/Pi0RpPb_PCM_2017_08_25.root";
  TString fileNameRpPbDalitz                  = "ExternalInputpPb/InputRpPb/Pi0RpPb_Dalitz_2017_08_25.root";
  TString fileNameRpPbPHOS                    = "ExternalInputpPb/InputRpPb/Pi0RpPb_PHOS_2017_08_25.root";
  TString fileNameRpPbEMCal                   = "ExternalInputpPb/InputRpPb/Pi0RpPb_EMCal_2017_08_25.root";
  TString fileNameRpPbPCMEMCal                = "ExternalInputpPb/InputRpPb/Pi0RpPb_PCM-EMCal_2017_08_25.root";
  //File to load the correlation factors
  //TString corrFactorsFileName 		      = "ExternalInputpPb/InputRpPb/CorrelationFactors/RpPb_5.023TeV_2016_09_27.root";
  //TString corrFactorsFileName                 = "eps/2017_02_14/ComputeCorrelationFactors_pPb5TeV/pPb5TeV.root";
  //TString corrFactorsFileName                 = "eps/2017_04_05/ComputeCorrelationFactors_pPb5TeV/pPb5TeV.root";
  
   TString corrFactorsFileName                 = "eps/2017_06_07/ComputeCorrelationFactors_pPb5TeV/pPb5TeV.root";
  
  
  TFile* fileNeutralPionRpPb                           	= new TFile(fileNameRpPb.Data());
  TFile* fileNeutralPionRpPbComb                        = new TFile(fileNameRpPbCombined.Data());
  TFile* fileNeutralPionRpPbPCM                         = new TFile(fileNameRpPbPCM.Data());
  TFile* fileNeutralPionRpPbDalitz                      = new TFile(fileNameRpPbDalitz.Data());
  TFile* fileNeutralPionRpPbPHOS                        = new TFile(fileNameRpPbPHOS.Data());
  TFile* fileNeutralPionRpPbEMCal                       = new TFile(fileNameRpPbEMCal.Data());
  TFile* fileNeutralPionRpPbPCMEMCal		        = new TFile(fileNameRpPbPCMEMCal.Data());
    
  TString nameHistoPCM                               = "Pi0_RpPb_PCM_StatErr";
  TString nameHistoPCMSysErrors                      = "Pi0_RpPb_PCM_SystErr";
  TString nameHistoDalitz                            = "Pi0_RpPb_Dalitz_StatErr";
  TString nameHistoDalitzSysErrors                   = "Pi0_RpPb_Dalitz_SystErr";
  TString nameHistoPHOS                              = "Pi0_RpPb_PHOS_StatErr";
  TString nameHistoPHOSSysErrors                     = "Pi0_RpPb_PHOS_SystErr";
  TString nameHistoEMCal                             = "Pi0_RpPb_EMCal_StatErr";
  TString nameHistoEMCalSysErrors                    = "Pi0_RpPb_EMCal_SystErr";
  TString nameHistoPCMEMCal			     = "Pi0_RpPb_PCM-EMCal_StatErr";
  TString nameHistoPCMEMCalSysErrors		     = "Pi0_RpPb_PCM-EMCal_SystErr";
  TString nameHistoCombSysErrors                     = "Pi0_RpPb_Comb_SystErr";
  TString nameHistoComb                              = "Pi0_RpPb_Comb_StatErr";
  TString nameHistoPCMppReferenceStat                = "Pi0_pp_reference_PCMBinning_StatErr";
  TString nameHistoPCMppReferenceSyst                = "Pi0_pp_reference_PCMBinning_SystErr";
  TString nameHistoDalitzppReferenceStat             = "Pi0_pp_reference_DalitzBinning_StatErr";
  TString nameHistoDalitzppReferenceSyst             = "Pi0_pp_reference_DalitzBinning_SystErr";
  TString nameHistoPHOSppReferenceStat               = "Pi0_pp_reference_PHOSBinning_StatErr";
  TString nameHistoPHOSppReferenceSyst               = "Pi0_pp_reference_PHOSBinning_SystErr";
  TString nameHistoEMCalppReferenceStat              = "Pi0_pp_reference_EMCalBinning_StatErr";
  TString nameHistoEMCalppReferenceSyst              = "Pi0_pp_reference_EMCalBinning_SystErr";
  TString nameHistoPCMEMCalppReferenceStat           = "Pi0_pp_reference_PCM-EMCalBinning_StatErr";
  TString nameHistoPCMEMCalppReferenceSyst	     = "Pi0_pp_reference_PCM-EMCalBinning_SystErr";
  TString nameHistoCombppReferenceStat               = "Pi0_pp_reference_CombBinning_StatErr";
  TString nameHistoCombppReferenceSyst               = "Pi0_pp_reference_CombBinning_SystErr";
  TString nameHistoPCMAlpha                          = "Pi0_RpPb_PCM_Alpha";
  TString nameHistoDalitzAlpha                       = "Pi0_RpPb_Dalitz_Alpha";
  TString nameHistoPHOSAlpha                         = "Pi0_RpPb_PHOS_Alpha";
  TString nameHistoEMCalAlpha                        = "Pi0_RpPb_EMCal_Alpha";
  TString nameHistoPCMEMCalAlpha		     = "Pi0_RpPb_PCM-EMCal_Alpha";
  TString nameHistoCombAlpha                         = "Pi0_RpPb_Comb_Alpha";
  
 
   
    
  TString collisionSystempPb                          = "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV"; 
    
    
  // Double_t pTLimits[nPtLimits]                        = { 0.3, 0.4, 0.5, 0.6, 0.7, 
  //                                                         0.8, 1.0, 1.2, 1.4, 1.6, 
  Double_t pTLimits[nPtLimits]                        = { 1.0, 1.2, 1.4, 1.6, 
							  1.8, 2.0, 2.2, 2.4, 2.6,
							  2.8, 3.0, 3.2, 3.4, 3.6,
							  3.8, 4.0, 4.5, 5.0, 5.5,
							  6.0, 7.0, 8.0, 10.0,12.0,
							  16.0, 20.0};
                                
    
  //    Int_t offSets[11]                             =  { 0, 6, 8, 0, 0,  3, 0, 0, 0,  0, 0};
  // Int_t offSetsSys[11]                             =  {  0, 6,  8, 0, 0, 3, 0, 0, 0, 0, 0};
  Int_t offSets[11]                                   =  { -6, 0, 2, 0, -1,  -3, 0, 0, 0,  0, 0};
  Int_t offSetsSys[11]                                =  { -6, 0, 2, 0, -1, -3, 0, 0, 0, 0, 0};
    
    
  TH1D* statErrorCollection[11];
  for (Int_t i = 0; i< 11; i++){
    statErrorCollection[i]                          = NULL;
  }    
    
  TGraphAsymmErrors* sysErrorCollection[11];
  for (Int_t i = 0; i< 11; i++){
    sysErrorCollection[i]                           = NULL;
  }    
  
  
  TGraphAsymmErrors* sysErrorPPReference[11];   //NOTE this histo is to compute the correlation factor as a function of pT
  
    for (Int_t i = 0; i< 11; i++){
      sysErrorPPReference[i]                           = NULL;
    }    
     
 
  // **************************************************************************************
  // ****************************** Reading Dalitz ****************************************
  // **************************************************************************************

  TGraphAsymmErrors* graphDalitzYieldPi0pPb           = (TGraphAsymmErrors*)fileNeutralPionRpPbDalitz->Get(nameHistoDalitz.Data());
  TH1D* histoDalitzYieldPi0pPb                        = (TH1D*)GraphAsymErrorsToHist_withErrors(graphDalitzYieldPi0pPb,nameHistoDalitz.Data());
  TGraphAsymmErrors* graphDalitzYieldPi0pPbSystErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbDalitz->Get(nameHistoDalitzSysErrors.Data());
  TGraphAsymmErrors* graphDalitzppreferenceStatErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbDalitz->Get(nameHistoDalitzppReferenceStat.Data());
  TGraphAsymmErrors* graphDalitzppreferenceSystErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbDalitz->Get(nameHistoDalitzppReferenceSyst.Data());
  TGraphAsymmErrors* graphDalitzAlpha                 = (TGraphAsymmErrors*)fileNeutralPionRpPbDalitz->Get( nameHistoDalitzAlpha.Data());

  // **************************************************************************************
  // ****************************** Reading PCM *******************************************
  // **************************************************************************************    
  TGraphAsymmErrors* graphPCMYieldPi0pPb           = (TGraphAsymmErrors*)fileNeutralPionRpPbPCM->Get(nameHistoPCM.Data());
  TH1D* histoPCMYieldPi0pPb           = (TH1D*)GraphAsymErrorsToHist_withErrors(graphPCMYieldPi0pPb,nameHistoPCM.Data());
  TGraphAsymmErrors* graphPCMYieldPi0pPbSystErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbPCM->Get(nameHistoPCMSysErrors.Data());
  TGraphAsymmErrors* graphPCMppreferenceStatErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbPCM->Get(nameHistoPCMppReferenceStat.Data());
  TGraphAsymmErrors* graphPCMppreferenceSystErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbPCM->Get(nameHistoPCMppReferenceSyst.Data());
  TGraphAsymmErrors* graphPCMAlpha                 = (TGraphAsymmErrors*)fileNeutralPionRpPbPCM->Get( nameHistoPCMAlpha.Data());    
  // **************************************************************************************
  // ******************************* Reading PHOS *****************************************
  // **************************************************************************************    
  
  TGraphAsymmErrors* graphPHOSYieldPi0pPb                      = (TGraphAsymmErrors*)fileNeutralPionRpPbPHOS->Get(nameHistoPHOS.Data()); 
  TH1D* histoPHOSYieldPi0pPb           = GraphAsymErrorsToHist_withErrors(graphPHOSYieldPi0pPb,nameHistoPHOS.Data());
  TGraphAsymmErrors* graphPHOSYieldPi0pPbSystErr                   = (TGraphAsymmErrors*)fileNeutralPionRpPbPHOS->Get(nameHistoPHOSSysErrors.Data());
  TGraphAsymmErrors* graphPHOSppreferenceStatErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbPHOS->Get(nameHistoPHOSppReferenceStat.Data());
  TGraphAsymmErrors* graphPHOSppreferenceSystErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbPHOS->Get(nameHistoPHOSppReferenceSyst.Data()); 
  TGraphAsymmErrors* graphPHOSAlpha                 = (TGraphAsymmErrors*)fileNeutralPionRpPbPHOS->Get( nameHistoPHOSAlpha.Data());
  // **************************************************************************************
  // ******************************** Reading EMCal ***************************************
  // **************************************************************************************
    
  TGraphAsymmErrors* graphEMCalYieldPi0pPb       = (TGraphAsymmErrors*)fileNeutralPionRpPbEMCal->Get(nameHistoEMCal.Data());
  TH1D* histoEMCalYieldPi0pPb           = GraphAsymErrorsToHist_withErrors(graphEMCalYieldPi0pPb,nameHistoEMCal.Data());
  TGraphAsymmErrors* graphEMCalYieldPi0pPbSystErr     = (TGraphAsymmErrors*)fileNeutralPionRpPbEMCal->Get(nameHistoEMCalSysErrors.Data());  
  TGraphAsymmErrors* graphEMCalppreferenceStatErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbEMCal->Get(nameHistoEMCalppReferenceStat.Data());
  TGraphAsymmErrors* graphEMCalppreferenceSystErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbEMCal->Get(nameHistoEMCalppReferenceSyst.Data()); 
  TGraphAsymmErrors* graphEMCalAlpha                 = (TGraphAsymmErrors*)fileNeutralPionRpPbEMCal->Get( nameHistoEMCalAlpha.Data());
  cout <<"stat"<< endl;
  graphEMCalYieldPi0pPb->Print();
  cout <<"sys"<< endl;
  graphEMCalYieldPi0pPbSystErr->Print();
  // ******************************** Reading PCM-EMCal ***************************************
  // **************************************************************************************
  TGraphAsymmErrors* graphPCMEMCalYieldPi0pPb       	= (TGraphAsymmErrors*)fileNeutralPionRpPbPCMEMCal->Get(nameHistoPCMEMCal.Data());
  TH1D* histoPCMEMCalYieldPi0pPb          		= GraphAsymErrorsToHist_withErrors(graphPCMEMCalYieldPi0pPb,nameHistoPCMEMCal.Data());
  TGraphAsymmErrors* graphPCMEMCalYieldPi0pPbSystErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbPCMEMCal->Get(nameHistoPCMEMCalSysErrors.Data());  
  TGraphAsymmErrors* graphPCMEMCalppreferenceStatErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbPCMEMCal->Get(nameHistoPCMEMCalppReferenceStat.Data());
  TGraphAsymmErrors* graphPCMEMCalppreferenceSystErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbPCMEMCal->Get(nameHistoPCMEMCalppReferenceSyst.Data()); 
  TGraphAsymmErrors* graphPCMEMCalAlpha                 = (TGraphAsymmErrors*)fileNeutralPionRpPbPCMEMCal->Get(nameHistoPCMEMCalAlpha.Data());
  cout <<"stat"<< endl;
  graphPCMEMCalYieldPi0pPb->Print();
  cout <<"sys"<< endl;
  graphPCMEMCalYieldPi0pPbSystErr->Print();
  
  // **************************************************************************************
  // ******************************** Reading RpPb from Combined pPb spectrum *************
  // **************************************************************************************
    
  TGraphAsymmErrors* graphCombYieldPi0pPb       = (TGraphAsymmErrors*)fileNeutralPionRpPbComb->Get(nameHistoComb.Data());
  TH1D* histoCombYieldPi0pPb           = GraphAsymErrorsToHist_withErrors(graphCombYieldPi0pPb,nameHistoComb.Data());
  TGraphAsymmErrors* graphCombYieldPi0pPbSystErr     = (TGraphAsymmErrors*)fileNeutralPionRpPbComb->Get(nameHistoCombSysErrors.Data());  
  TGraphAsymmErrors* graphCombppreferenceStatErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbComb->Get(nameHistoCombppReferenceStat.Data());
  TGraphAsymmErrors* graphCombppreferenceSystErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbComb->Get(nameHistoCombppReferenceSyst.Data()); 
  TGraphAsymmErrors* graphCombAlpha                 = (TGraphAsymmErrors*)fileNeutralPionRpPbComb->Get( nameHistoCombAlpha.Data());

  // **************************************************************************************
  // ******************************** Reading Model calculations **************************
  // **************************************************************************************
  TGraphAsymmErrors*  graphAsymmErrorsPi0DSS5000    = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get("EPS09s_fDSS_errors"); 
  TGraph*  graphPi0CGC                     = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get("ColorGlasCondensate"); 
  TGraph*  graphPi0DSS5000                 = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get("EPS09s_fDSS_NLO"); 
  TGraph*  graphPi0ESP09sPi0AKK            = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get("EPS09s_AKK_NLO"); 
  TGraph*  graphPi0ESP09sPi0KKP            = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get("EPS09s_KKP_NLO"); 


  graphAsymmErrorsPi0DSS5000->Draw("same,E3");
  graphPi0CGC->Draw("same,p");
  graphPi0DSS5000->Draw("same,p,l");
  graphPi0ESP09sPi0AKK->Draw("same,p,l");
  graphPi0ESP09sPi0KKP->Draw("same,p,l");
  // **************************************************************************************
  // ******************************** Reading Charged Particles and Pions **************************
  // **************************************************************************************
  TGraphAsymmErrors*  graphRpPbChargedParticlesStatErr  = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get("RpPb_ChargedParticles_StatErr"); 
  TGraphAsymmErrors*  graphRpPbChargedParticlesSystErr  = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get("RpPb_ChargedParticles_SystErr"); 
  TGraphAsymmErrors*  graphRpPbChargedPionsStatErr    	= (TGraphAsymmErrors*)fileNeutralPionRpPb->Get("RpPb_ChargedPions_StatErr"); 
  TGraphAsymmErrors*  graphRpPbChargedPionsSystErr      = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get("RpPb_ChargedPions_SystErr"); 
  // **************************************************************************************
  // ********************************* Combine spectra ************************************
  // **************************************************************************************
    
  TString fileNameOutputWeightingOld                = Form("%s/WeightingRpPb.dat",outputDir.Data());
  
  

  statErrorCollection[0]          = (TH1D*)histoPCMYieldPi0pPb->Clone("statErrPCMPi0");
  statErrorCollection[1]          = (TH1D*)histoPHOSYieldPi0pPb->Clone("statErrPHOSPi0");
  statErrorCollection[2]          = (TH1D*)histoEMCalYieldPi0pPb->Clone("statErrEMCalPi0");
  statErrorCollection[4]	  = (TH1D*)histoPCMEMCalYieldPi0pPb->Clone("statErrPCM-EMCalPi0");
  statErrorCollection[5]          = (TH1D*)histoDalitzYieldPi0pPb->Clone("statErrDalitzPi0");
    
  sysErrorCollection[0]           = (TGraphAsymmErrors*)graphPCMYieldPi0pPbSystErr->Clone("sysErrPCMPi0");
  sysErrorCollection[1]           = (TGraphAsymmErrors*)graphPHOSYieldPi0pPbSystErr->Clone("sysErrPHOSPi0");
  sysErrorCollection[2]           = (TGraphAsymmErrors*)graphEMCalYieldPi0pPbSystErr->Clone("sysErrEMCalPi0");
  sysErrorCollection[4]		  = (TGraphAsymmErrors*)graphPCMEMCalYieldPi0pPbSystErr->Clone("sysErrPCM-EMCalPi0");
  sysErrorCollection[5]           = (TGraphAsymmErrors*)graphDalitzYieldPi0pPbSystErr->Clone("sysErrDalitzPi0");
  cout << "Sys Error PCM:" <<endl;
  graphPCMYieldPi0pPbSystErr->Print();    
  cout << "Sys Error PHOS:" <<endl;
  graphPHOSYieldPi0pPbSystErr->Print();  
  cout << "Sys Error Dalitz:" <<endl;
  graphDalitzYieldPi0pPbSystErr->Print();

  sysErrorPPReference[0]  = (TGraphAsymmErrors*)graphPCMppreferenceSystErr->Clone("sysErrPPRefPCMPi0");
  sysErrorPPReference[1]  = (TGraphAsymmErrors*)graphPHOSppreferenceSystErr->Clone("sysErrPPRefPHOSPi0");
  sysErrorPPReference[2]  = (TGraphAsymmErrors*)graphEMCalppreferenceSystErr->Clone("sysErrPPRefEMCalPi0");
  sysErrorPPReference[4]  = (TGraphAsymmErrors*)graphPCMEMCalppreferenceSystErr->Clone("sysErrPPRefPCM-EMCalPi0");
  sysErrorPPReference[5]  = (TGraphAsymmErrors*)graphDalitzppreferenceSystErr->Clone("sysErrPPRefDalitzPi0");
    
  TGraphAsymmErrors* graphCombPi0InvCrossSectionStatpPb5023GeV= NULL;
  TGraphAsymmErrors* graphCombPi0InvCrossSectionSyspPb5023GeV = NULL;
  
  cout<<"\n\n */////////////////Combination of Pi0 RpPb  begins here**//////////////"<<endl;
 
    
  TGraphAsymmErrors* graphCombPi0InvCrossSectionTotpPb5023GeV = CombinePtPointsSpectraFullCorrMat(   statErrorCollection,    sysErrorCollection,     
												      pTLimits, Ntotal,
												      offSets, offSetsSys,
												      graphCombPi0InvCrossSectionStatpPb5023GeV, graphCombPi0InvCrossSectionSyspPb5023GeV,
												      fileNameOutputWeightingOld,"pPb_5.023GeV_RpPb","Pi0",kFALSE,sysErrorPPReference,corrFactorsFileName.Data()
												      ); 
  cout<<"\n\n */////////////////Combination of Pi0 RpPb  ends here**//////////////"<<endl;
 
  graphCombPi0InvCrossSectionStatpPb5023GeV->Print();

  if (IsNSD){

    graphCombPi0InvCrossSectionStatpPb5023GeV=ScaleGraph(graphCombPi0InvCrossSectionStatpPb5023GeV,Scaling);
    graphCombPi0InvCrossSectionSyspPb5023GeV=ScaleGraph(graphCombPi0InvCrossSectionSyspPb5023GeV,Scaling);

    graphCombPi0InvCrossSectionTotpPb5023GeV  =  CalculateCombinedSysAndStatError( graphCombPi0InvCrossSectionStatpPb5023GeV ,graphCombPi0InvCrossSectionSyspPb5023GeV );
    
    cout<< "Combined RpPb NSD scaled:" << endl;
    graphCombPi0InvCrossSectionTotpPb5023GeV->Print();
    graphCombPi0InvCrossSectionStatpPb5023GeV->Print();
    graphCombPi0InvCrossSectionSyspPb5023GeV->Print();
  }

    
  TGraphAsymmErrors* graphInvYieldPi0CombpPb5023GeVStaClone   = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionStatpPb5023GeV->Clone();
  TGraphAsymmErrors* graphInvYieldPi0CombpPb5023GeVSysClone   = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionSyspPb5023GeV->Clone();
  TGraphAsymmErrors* graphInvYieldPi0CombpPb5023GeVTotClone   = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionTotpPb5023GeV->Clone();
    
  // **************************************************************************************
  // ************************* Plotting Combined R_pPb ************************
  // **************************************************************************************
  TCanvas* canvasCombRpPb = new TCanvas("canvasCombRpPb","",200,10,1200,700);  // gives the page size
  DrawGammaCanvasSettings( canvasCombRpPb, 0.09, 0.02, 0.02, 0.13);
    
  //  canvasCombRpPb->SetLogx();
  // canvasCombRpPb->SetLogy();
  TH2F * histo2DCombined;
  histo2DCombined = new TH2F("histo2DCombined","histo2DCombined",1000,0.,20.,1000,0.3,1.7);
  SetStyleHistoTH2ForGraphs(histo2DCombined, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#pi^{0}}_{p-Pb}", 0.05,0.06, 0.05,0.06, 0.9,0.6, 512, 505);
  histo2DCombined->DrawCopy();
    

     
    

  DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0CombpPb5023GeVSysClone,20,1.5, 4, 4, 1, kTRUE);  
      
  graphInvYieldPi0CombpPb5023GeVSysClone->Draw("E2,same");

  DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0CombpPb5023GeVStaClone,20,1.5, kBlue, kBlue);  
  graphInvYieldPi0CombpPb5023GeVStaClone->Draw("p,same");


  TLatex * lt = new TLatex(1.,1.55,"ALICE Preliminary"); 
  lt->SetTextSize(0.05) ;
  //  lt->DrawText(1.,1.45,"2016/06/14") ;
  lt->Draw("same");
     
  TLegend* legendRpPbCombine = new TLegend(0.1,0.75,0.55,0.85);
  legendRpPbCombine->SetFillColor(0);
  legendRpPbCombine->SetLineColor(0);
  legendRpPbCombine->SetTextSize(0.04);
  legendRpPbCombine->AddEntry(graphInvYieldPi0CombpPb5023GeVSysClone,"NSD, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");

    
  legendRpPbCombine->Draw();
  
  TLine* line =new TLine(0.,1.,20.,1.);
  line->Draw("same");
  
  TBox* BoxNorm =new TBox(0.25, 1.-NormalizationError, 0.5, 1.+NormalizationError);
  BoxNorm->SetFillColor(4);
  BoxNorm->Draw("same");
    
  canvasCombRpPb->Print(Form("%s/Comb_RpPb.%s",outputDir.Data(),suffix.Data()));
    
    
  // **************************************************************************************
  // ************************* Plotting All GA R_pPb ************************
  // **************************************************************************************
  //Apply NSD scaling to individual RpPb

  if (IsNSD){
    graphDalitzYieldPi0pPb=ScaleGraph(graphDalitzYieldPi0pPb,Scaling);
    graphDalitzYieldPi0pPbSystErr=ScaleGraph(graphDalitzYieldPi0pPbSystErr,Scaling);
    //   graphDalitzYieldPi0pPbSystErr=ApplyNSDSysError(graphDalitzYieldPi0pPbSystErr,ScalingErr,ScalingErrTpPb);

    graphPCMYieldPi0pPb=ScaleGraph(graphPCMYieldPi0pPb,Scaling);
    graphPCMYieldPi0pPbSystErr=ScaleGraph(graphPCMYieldPi0pPbSystErr,Scaling);
    //   graphPCMYieldPi0pPbSystErr=ApplyNSDSysError(graphPCMYieldPi0pPbSystErr,ScalingErr,ScalingErrTpPb);

    graphPHOSYieldPi0pPb=ScaleGraph(graphPHOSYieldPi0pPb,Scaling);
    graphPHOSYieldPi0pPbSystErr=ScaleGraph(graphPHOSYieldPi0pPbSystErr,Scaling);
    //  graphPHOSYieldPi0pPbSystErr=ApplyNSDSysError(graphPHOSYieldPi0pPbSystErr,ScalingErr,ScalingErrTpPb);

    graphEMCalYieldPi0pPb=ScaleGraph(graphEMCalYieldPi0pPb,Scaling);
    graphEMCalYieldPi0pPbSystErr=ScaleGraph(graphEMCalYieldPi0pPbSystErr,Scaling);
    //  graphEMCalYieldPi0pPbSystErr=ApplyNSDSysError(graphEMCalYieldPi0pPbSystErr,ScalingErr,ScalingErrTpPb);
  }//  else{
  //   graphDalitzYieldPi0pPbSystErr=ApplyTpPbSysError(graphDalitzYieldPi0pPbSystErr,ScalingErrTpPb);
  //   graphPCMYieldPi0pPbSystErr=ApplyTpPbSysError(graphPCMYieldPi0pPbSystErr,ScalingErrTpPb);
  //   graphPHOSYieldPi0pPbSystErr=ApplyTpPbSysError(graphPHOSYieldPi0pPbSystErr,ScalingErrTpPb);
  //   graphEMCalYieldPi0pPbSystErr=ApplyTpPbSysError(graphEMCalYieldPi0pPbSystErr,ScalingErrTpPb);
  // }  //Normalization uncertainties come as error box!!


  TCanvas* canvasAllGARpPb = new TCanvas("canvasAllGARpPb","",200,10,1200,700);  // gives the page size
  DrawGammaCanvasSettings( canvasAllGARpPb, 0.08, 0.02, 0.02, 0.13);
    
  //  canvasAllGARpPb->SetLogx();
  // canvasAllGARpPb->SetLogy();
  TH2F * histo2DAllGA;
  histo2DAllGA = new TH2F("histo2DAllGA","histo2DAllGA",1000,0.,20.,1000,0.3,2.);
  SetStyleHistoTH2ForGraphs(histo2DAllGA, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#pi^{0}}_{p-Pb}", 0.05,0.06, 0.05,0.05, 0.9,0.7, 512, 505);
  histo2DAllGA->DrawCopy();
    
  DrawGammaSetMarkerTGraphAsym(graphDalitzYieldPi0pPbSystErr,20,1.5, kCyan+2, kCyan+2, 1, kTRUE);
  DrawGammaSetMarkerTGraphAsym(graphDalitzYieldPi0pPb,20,1.5, kCyan+2, kCyan+2); 
  graphDalitzYieldPi0pPbSystErr->Draw("E2,same");  
  graphDalitzYieldPi0pPb->Draw("p,same");  

  DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0pPbSystErr,20,1.5, 1, 1, 1, kTRUE);
  DrawGammaSetMarkerTGraphAsym(graphPCMYieldPi0pPb,20,1.5, 1, 1); 
  graphPCMYieldPi0pPbSystErr->Draw("E2,same");  
  graphPCMYieldPi0pPb->Draw("p,same");  

  DrawGammaSetMarkerTGraphAsym(graphPHOSYieldPi0pPbSystErr,20,1.5, kRed+1,kRed+1, 1, kTRUE);
  DrawGammaSetMarkerTGraphAsym(graphPHOSYieldPi0pPb,20,1.5, kRed+1, kRed+1); 
  graphPHOSYieldPi0pPbSystErr->Draw("E2,same");  
  graphPHOSYieldPi0pPb->Draw("p,same");
 
  DrawGammaSetMarkerTGraphAsym(graphEMCalYieldPi0pPbSystErr,20,1.5, kGreen+2, kGreen+2, 1, kTRUE);
  graphEMCalYieldPi0pPbSystErr->Draw("E2,same");
  DrawGammaSetMarkerTGraphAsym(graphEMCalYieldPi0pPb,20,1.5, kGreen+2, kGreen+2); 
  graphEMCalYieldPi0pPb->Draw("p,same");
  
  DrawGammaSetMarkerTGraphAsym(graphPCMEMCalYieldPi0pPbSystErr,20,1.5, kOrange+1, kOrange+1, 1, kTRUE);
  graphPCMEMCalYieldPi0pPbSystErr->Draw("E2,same");
  DrawGammaSetMarkerTGraphAsym(graphPCMEMCalYieldPi0pPb,20,1.5, kOrange+1, kOrange+1); 
  graphPCMEMCalYieldPi0pPb->Draw("p,same");
    

  graphInvYieldPi0CombpPb5023GeVSysClone->Draw("E2,same");
  graphInvYieldPi0CombpPb5023GeVStaClone->Draw("p,same");


  TLegend* legendRpPbAllGA = new TLegend(0.17,0.65,0.4,0.95);
  legendRpPbAllGA->SetFillColor(0);
  legendRpPbAllGA->SetLineColor(0);
  legendRpPbAllGA->SetTextSize(0.03);
  legendRpPbAllGA->AddEntry(graphInvYieldPi0CombpPb5023GeVSysClone,"Combined","pef");
  legendRpPbAllGA->AddEntry(graphPCMYieldPi0pPbSystErr,"PCM","pef");
  legendRpPbAllGA->AddEntry(graphDalitzYieldPi0pPbSystErr,"Dalitz","pef");
  legendRpPbAllGA->AddEntry(graphPHOSYieldPi0pPbSystErr,"PHOS","pef");
  legendRpPbAllGA->AddEntry(graphEMCalYieldPi0pPbSystErr,"EMCal","pef");
  legendRpPbAllGA->AddEntry(graphPCMEMCalYieldPi0pPbSystErr,"PCM-EMCal","pef");

    
  legendRpPbAllGA->Draw();
  

  line->Draw("same");
 
    
  canvasAllGARpPb->Print(Form("%s/AllGA_RpPb.%s",outputDir.Data(),suffix.Data()));   
   
    
  // **************************************************************************************
  // ************************* PlottingCombined R_pPb and PCM************************
  // **************************************************************************************
  TCanvas* canvasAllGA1RpPb = new TCanvas("canvasAllGA1RpPb","",200,10,1200,700);  // gives the page size
  DrawGammaCanvasSettings( canvasAllGA1RpPb, 0.08, 0.02, 0.02, 0.13);
    
  //  canvasAllGA1RpPb->SetLogx();
  // canvasAllGA1RpPb->SetLogy();
  TH2F * histo2DAllGA1;
  histo2DAllGA1 = new TH2F("histo2DAllGA1","histo2DAllGA1",1000,0.,20.,1000,0.3,1.7);
  SetStyleHistoTH2ForGraphs(histo2DAllGA1, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#pi^{0}}_{p-Pb}", 0.05,0.06, 0.05,0.05, 0.9,0.7, 512, 505);
  histo2DAllGA1->DrawCopy();
  
  graphInvYieldPi0CombpPb5023GeVSysClone->Draw("E2,same");
  graphInvYieldPi0CombpPb5023GeVStaClone->Draw("p,same");
  
  graphPCMYieldPi0pPbSystErr-> SetMarkerStyle(24);
  graphPCMYieldPi0pPb->SetMarkerStyle(24);
  graphPCMYieldPi0pPbSystErr->Draw("E2,same");  
  graphPCMYieldPi0pPb->Draw("p,same");

  TLegend* legendRpPbAllGA1 = new TLegend(0.17,0.2,0.4,0.35);
  legendRpPbAllGA1->SetFillColor(0);
  legendRpPbAllGA1->SetLineColor(0);
  legendRpPbAllGA1->SetTextSize(0.03);
  legendRpPbAllGA1->AddEntry(graphInvYieldPi0CombpPb5023GeVSysClone,"Combined","pef");
  legendRpPbAllGA1->AddEntry(graphPCMYieldPi0pPbSystErr,"PCM","pef");
    
  legendRpPbAllGA1->Draw();
  

  line->Draw("same");
 
    
  canvasAllGA1RpPb->Print(Form("%s/CombinedAndPCM_RpPb.%s",outputDir.Data(),suffix.Data()));   
  // **************************************************************************************
  // ************************* PlottingCombined R_pPb and Dalitz************************
  // **************************************************************************************
  TCanvas* canvasAllGA2RpPb = new TCanvas("canvasAllGA2RpPb","",200,10,1200,700);  // gives the page size
  DrawGammaCanvasSettings( canvasAllGA2RpPb, 0.08, 0.02, 0.02, 0.13);
    
  //  canvasAllGA2RpPb->SetLogx();
  // canvasAllGA2RpPb->SetLogy();
  TH2F * histo2DAllGA2;
  histo2DAllGA2 = new TH2F("histo2DAllGA2","histo2DAllGA2",1000,0.,20.,1000,0.3,1.7);
  SetStyleHistoTH2ForGraphs(histo2DAllGA2, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#pi^{0}}_{p-Pb}", 0.05,0.06, 0.05,0.05, 0.9,0.7, 512, 505);
  histo2DAllGA2->DrawCopy();

  graphInvYieldPi0CombpPb5023GeVSysClone->Draw("E2,same");
  graphInvYieldPi0CombpPb5023GeVStaClone->Draw("p,same");

  graphDalitzYieldPi0pPbSystErr->SetMarkerStyle(24);
  graphDalitzYieldPi0pPb->SetMarkerStyle(24);
  graphDalitzYieldPi0pPbSystErr->Draw("E2,same");  
  graphDalitzYieldPi0pPb->Draw("p,same");  

  TLegend* legendRpPbAllGA2 = new TLegend(0.17,0.2,0.4,0.35);
  legendRpPbAllGA2->SetFillColor(0);
  legendRpPbAllGA2->SetLineColor(0);
  legendRpPbAllGA2->SetTextSize(0.03);
  legendRpPbAllGA2->AddEntry(graphInvYieldPi0CombpPb5023GeVSysClone,"Combined","pef");
  legendRpPbAllGA2->AddEntry(graphDalitzYieldPi0pPbSystErr,"Dalitz","pef");
        
  legendRpPbAllGA2->Draw();
  

  line->Draw("same");
 
    
  canvasAllGA2RpPb->Print(Form("%s/CombinedAndDalitz_RpPb.%s",outputDir.Data(),suffix.Data())); 
  // **************************************************************************************
  // ************************* PlottingCombined R_pPb and PHOS************************
  // **************************************************************************************
  TCanvas* canvasAllGA3RpPb = new TCanvas("canvasAllGA3RpPb","",200,10,1200,700);  // gives the page size
  DrawGammaCanvasSettings( canvasAllGA3RpPb, 0.08, 0.02, 0.02, 0.13);
    
  //  canvasAllGA3RpPb->SetLogx();
  // canvasAllGA3RpPb->SetLogy();
  TH2F * histo2DAllGA3;
  histo2DAllGA3 = new TH2F("histo2DAllGA3","histo2DAllGA3",1000,0.,20.,1000,0.3,1.7);
  SetStyleHistoTH2ForGraphs(histo2DAllGA3, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#pi^{0}}_{p-Pb}", 0.05,0.06, 0.05,0.05, 0.9,0.7, 512, 505);
  histo2DAllGA3->DrawCopy();
    


  graphInvYieldPi0CombpPb5023GeVSysClone->Draw("E2,same");
  graphInvYieldPi0CombpPb5023GeVStaClone->Draw("p,same");

 
  graphPHOSYieldPi0pPbSystErr->SetMarkerStyle(24);
  graphPHOSYieldPi0pPb->SetMarkerStyle(24);
  graphPHOSYieldPi0pPbSystErr->Draw("E2,same");  
  graphPHOSYieldPi0pPb->Draw("p,same"); 




  TLegend* legendRpPbAllGA3 = new TLegend(0.17,0.2,0.4,0.35);
  legendRpPbAllGA3->SetFillColor(0);
  legendRpPbAllGA3->SetLineColor(0);
  legendRpPbAllGA3->SetTextSize(0.03);
  legendRpPbAllGA3->AddEntry(graphInvYieldPi0CombpPb5023GeVSysClone,"Combined","pef");
  legendRpPbAllGA3->AddEntry(graphPHOSYieldPi0pPbSystErr,"PHOS","pef");
        
  legendRpPbAllGA3->Draw();
  

  line->Draw("same");
 
    
  canvasAllGA3RpPb->Print(Form("%s/CombinedAndPHOS_RpPb.%s",outputDir.Data(),suffix.Data())); 

  // **************************************************************************************
  // ************************* PlottingCombined R_pPb and EMCal************************
  // **************************************************************************************
  TCanvas* canvasAllGA4RpPb = new TCanvas("canvasAllGA4RpPb","",200,10,1200,700);  // gives the page size
  DrawGammaCanvasSettings( canvasAllGA4RpPb, 0.08, 0.02, 0.02, 0.13);
    
  //  canvasAllGA4RpPb->SetLogx();
  // canvasAllGA4RpPb->SetLogy();
  TH2F * histo2DAllGA4;
  histo2DAllGA4 = new TH2F("histo2DAllGA4","histo2DAllGA4",1000,0.,20.,1000,0.3,1.7);
  SetStyleHistoTH2ForGraphs(histo2DAllGA4, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#pi^{0}}_{p-Pb}", 0.05,0.06, 0.05,0.05, 0.9,0.7, 512, 505);
  histo2DAllGA4->DrawCopy();

  graphInvYieldPi0CombpPb5023GeVSysClone->Draw("E2,same");
  graphInvYieldPi0CombpPb5023GeVStaClone->Draw("p,same");   

  DrawGammaSetMarkerTGraphAsym(graphEMCalYieldPi0pPbSystErr,24,1.5, kGreen+2, kGreen+2, 1, kTRUE);
  graphEMCalYieldPi0pPb->SetMarkerStyle(24);
  graphEMCalYieldPi0pPbSystErr->Draw("E2,same");  
  graphEMCalYieldPi0pPb->Draw("p,same");  



  TLegend* legendRpPbAllGA4 = new TLegend(0.17,0.2,0.4,0.35);
  legendRpPbAllGA4->SetFillColor(0);
  legendRpPbAllGA4->SetLineColor(0);
  legendRpPbAllGA4->SetTextSize(0.03);
  legendRpPbAllGA4->AddEntry(graphInvYieldPi0CombpPb5023GeVSysClone,"Combined","pef");
  legendRpPbAllGA4->AddEntry(graphEMCalYieldPi0pPbSystErr,"EMCal","pef");
        
  legendRpPbAllGA4->Draw();
  

  line->Draw("same");
 
    
  canvasAllGA4RpPb->Print(Form("%s/CombinedAndEMCal_RpPb.%s",outputDir.Data(),suffix.Data())); 
    // **************************************************************************************
  // ************************* PlottingCombined R_pPb and PCM-EMCal************************
  // **************************************************************************************
  TCanvas* canvasAllGA5RpPb = new TCanvas("canvasAllGA5RpPb","",200,10,1200,700);  // gives the page size
  DrawGammaCanvasSettings( canvasAllGA5RpPb, 0.08, 0.02, 0.02, 0.13);
  TH2F * histo2DAllGA5;
  histo2DAllGA5 = new TH2F("histo2DAllGA5","histo2DAllGA5",1000,0.,20.,1000,0.3,1.7);
  SetStyleHistoTH2ForGraphs(histo2DAllGA5, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#pi^{0}}_{p-Pb}", 0.05,0.06, 0.05,0.05, 0.9,0.7, 512, 505);
  histo2DAllGA5->DrawCopy();

  graphInvYieldPi0CombpPb5023GeVSysClone->Draw("E2,same");
  graphInvYieldPi0CombpPb5023GeVStaClone->Draw("p,same");   

  DrawGammaSetMarkerTGraphAsym(graphPCMEMCalYieldPi0pPbSystErr,24,1.5, kGreen+2, kGreen+2, 1, kTRUE);
  graphPCMEMCalYieldPi0pPb->SetMarkerStyle(24);
  graphPCMEMCalYieldPi0pPbSystErr->Draw("E2,same");  
  graphPCMEMCalYieldPi0pPb->Draw("p,same");  



  TLegend* legendRpPbAllGA5 = new TLegend(0.17,0.2,0.4,0.35);
  legendRpPbAllGA5->SetFillColor(0);
  legendRpPbAllGA5->SetLineColor(0);
  legendRpPbAllGA5->SetTextSize(0.03);
  legendRpPbAllGA5->AddEntry(graphInvYieldPi0CombpPb5023GeVSysClone,"Combined","pef");
  legendRpPbAllGA5->AddEntry(graphPCMEMCalYieldPi0pPbSystErr,"PCM-EMCal","pef");
        
  legendRpPbAllGA5->Draw();
  

  line->Draw("same");
 
    
  canvasAllGA5RpPb->Print(Form("%s/CombinedAndPCM-EMCal_RpPb.%s",outputDir.Data(),suffix.Data())); 
  
  // **************************************************************************************
  // ************************* Plotting Combined R_pPb with models ************************
  // **************************************************************************************
  TCanvas* canvasCombinedWithModelsRpPb = new TCanvas("canvasCombinedWithModelsRpPb","",200,10,1200,700);  // gives the page size
  DrawGammaCanvasSettings( canvasCombinedWithModelsRpPb, 0.09, 0.02, 0.02, 0.13);
    
  //  canvasCombinedWithModelsRpPb->SetLogx();
  // canvasCombinedWithModelsRpPb->SetLogy();
  TH2F * histo2DCombinedWithModels;
  histo2DCombinedWithModels = new TH2F("histo2DCombinedWithModels","histo2DCombinedWithModels",1000,0.,20.,1000,0.3,1.5);
  SetStyleHistoTH2ForGraphs(histo2DCombinedWithModels, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#pi^{0}}_{p-Pb}", 0.05,0.06, 0.05,0.06, 0.9,0.6, 512, 505);
  histo2DCombinedWithModels->DrawCopy();
    
  graphAsymmErrorsPi0DSS5000->Draw("same,E3");
  graphPi0CGC->Draw("same,p");
  graphPi0DSS5000->Draw("same,p,l");
  graphPi0ESP09sPi0AKK->Draw("same,p,l");
  graphPi0ESP09sPi0KKP->Draw("same,p,l");
  graphInvYieldPi0CombpPb5023GeVSysClone->Draw("E2,same");
  graphInvYieldPi0CombpPb5023GeVStaClone->Draw("p,same");
  Float_t offsetYMaxLegendPCM = -0.15;
  Float_t xMinTheoryPCM = 0.30;

  TLegend* legendRpPbESP09sPCM = new TLegend(0.2,0.15,0.45,0.43);
  legendRpPbESP09sPCM->SetFillColor(0);
  legendRpPbESP09sPCM->SetLineColor(0);
  legendRpPbESP09sPCM->SetTextSize(0.04);
  legendRpPbESP09sPCM->SetTextFont(42);
  legendRpPbESP09sPCM->AddEntry(graphPi0ESP09sPi0KKP,"EPS09s KKP NLO");
  legendRpPbESP09sPCM->AddEntry(graphPi0ESP09sPi0AKK,"EPS09s AKK NLO");
  legendRpPbESP09sPCM->AddEntry(graphPi0DSS5000,"EPS09s fDSS NLO");
  legendRpPbESP09sPCM->AddEntry(graphAsymmErrorsPi0DSS5000,"EPS09s fDSS errors","pef");
	
	
  TLegend* legendRpPbCGCPCM = new TLegend(0.5,0.15,0.85,0.25);
  legendRpPbCGCPCM->SetFillColor(0);
  legendRpPbCGCPCM->SetLineColor(0);
  legendRpPbCGCPCM->SetTextSize(0.04);
  legendRpPbCGCPCM->SetTextFont(42);
  legendRpPbCGCPCM->AddEntry(graphPi0CGC,"CGC (Lappi and M\344ntysaari)","p");
  legendRpPbESP09sPCM->Draw("same");
  legendRpPbCGCPCM->Draw("same");

   TLegend* legendRpPbCombine2a= new TLegend(0.1,0.75,0.55,0.85);
  legendRpPbCombine2a->SetFillColor(0);
  legendRpPbCombine2a->SetLineColor(0);
  legendRpPbCombine2a->SetTextSize(0.04);
  legendRpPbCombine2a->AddEntry(graphInvYieldPi0CombpPb5023GeVSysClone,"NSD, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");

    
  legendRpPbCombine->Draw(); 

  line->Draw("same");
  TLatex * lt2 = new TLatex(1.,1.37,"ALICE Preliminary"); 
  lt2->SetTextSize(0.05) ;
  //  lt->DrawText(1.,1.45,"2016/06/14") ;
  lt2->Draw("same");
      BoxNorm->Draw("same");
  canvasCombinedWithModelsRpPb->Print(Form("%s/CombinedWithModels_RpPb.%s",outputDir.Data(),suffix.Data()));
 


  // **************************************************************************************
  // ************************* Plotting Combined R_pPb with Charged Particles ************************
  // **************************************************************************************
  TCanvas* canvasCombinedWithChargedParticlesRpPb = new TCanvas("canvasCombinedWithChargedParticlesRpPb","",200,10,1200,700);  // gives the page size
  DrawGammaCanvasSettings( canvasCombinedWithChargedParticlesRpPb, 0.08, 0.02, 0.02, 0.13);
    
  //  canvasCombinedWithChargedParticlesRpPb->SetLogx();
  // canvasCombinedWithChargedParticlesRpPb->SetLogy();
  TH2F * histo2DCombinedWithChargedParticles;
  histo2DCombinedWithChargedParticles = new TH2F("histo2DCombinedWithChargedParticles","histo2DCombinedWithChargedParticles",1000,0.,20.,1000,0.3,1.5);
  SetStyleHistoTH2ForGraphs(histo2DCombinedWithChargedParticles, "#it{p}_{T} (GeV/#it{c})","#it{R}_{p-Pb}", 0.05,0.06, 0.05,0.05, 0.9,0.7, 512, 505);
  histo2DCombinedWithChargedParticles->DrawCopy();
    
  graphRpPbChargedParticlesSystErr->Draw("E2,same");  
  graphRpPbChargedParticlesStatErr->Draw("p,same");  
  graphInvYieldPi0CombpPb5023GeVSysClone->Draw("E2,same");
  graphInvYieldPi0CombpPb5023GeVStaClone->Draw("p,same");
  
  TLegend* legendRpPbChargedParticles = new TLegend(0.1,0.85,0.55,0.95);
  legendRpPbChargedParticles->SetFillColor(0);
  legendRpPbChargedParticles->SetLineColor(0);
  legendRpPbChargedParticles->SetTextSize(0.03);
  legendRpPbChargedParticles->AddEntry(graphInvYieldPi0CombpPb5023GeVSysClone,"NSD #pi^{0}, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  legendRpPbChargedParticles->AddEntry(graphRpPbChargedParticlesSystErr,"NSD charged particles, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
   
  legendRpPbChargedParticles->Draw();
  BoxNorm->Draw("same");
  line->Draw("same");
 
    
  canvasCombinedWithChargedParticlesRpPb->Print(Form("%s/CombinedWithChargedParticles_RpPb.%s",outputDir.Data(),suffix.Data()));

  // **************************************************************************************
  // ************************* Plotting Combined R_pPb with Charged Pions ************************
  // **************************************************************************************
  TCanvas* canvasCombinedWithChargedPionsRpPb = new TCanvas("canvasCombinedWithChargedPionsRpPb","",200,10,1200,700);  // gives the page size
  DrawGammaCanvasSettings( canvasCombinedWithChargedPionsRpPb, 0.08, 0.02, 0.02, 0.13);
    
  //  canvasCombinedWithChargedPionsRpPb->SetLogx();
  // canvasCombinedWithChargedPionsRpPb->SetLogy();
  TH2F * histo2DCombinedWithChargedPions;
  histo2DCombinedWithChargedPions = new TH2F("histo2DCombinedWithChargedPions","histo2DCombinedWithChargedPions",1000,0.,20.,1000,0.3,1.5);
  SetStyleHistoTH2ForGraphs(histo2DCombinedWithChargedPions, "#it{p}_{T} (GeV/#it{c})","#it{R}_{p-Pb}", 0.05,0.06, 0.05,0.05, 0.9,0.7, 512, 505);
  histo2DCombinedWithChargedPions->DrawCopy();
    
  graphRpPbChargedPionsSystErr->Draw("E2,same");  
  graphRpPbChargedPionsStatErr->Draw("p,same");  
  graphInvYieldPi0CombpPb5023GeVSysClone->Draw("E2,same");
  graphInvYieldPi0CombpPb5023GeVStaClone->Draw("p,same");
  
  TLegend* legendRpPbChargedPions = new TLegend(0.1,0.85,0.55,0.95);
  legendRpPbChargedPions->SetFillColor(0);
  legendRpPbChargedPions->SetLineColor(0);
  legendRpPbChargedPions->SetTextSize(0.03);
  legendRpPbChargedPions->AddEntry(graphInvYieldPi0CombpPb5023GeVSysClone,"NSD #pi^{0}, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
  legendRpPbChargedPions->AddEntry(graphRpPbChargedPionsSystErr,"NSD #pi^{+/-}, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");
   
  legendRpPbChargedPions->Draw();
  BoxNorm->Draw("same");
  line->Draw("same");
 
    
  canvasCombinedWithChargedPionsRpPb->Print(Form("%s/CombinedWithChargedPions_RpPb.%s",outputDir.Data(),suffix.Data()));

  // **************************************************************************************
  // ************************* Plotting all ppreferences ************************
  // **************************************************************************************
  TCanvas* canvasAllppreferences = new TCanvas("canvasAllppreferences","",200,10,1200,700);  // gives the page size
  DrawGammaCanvasSettings(canvasAllppreferences , 0.12, 0.02, 0.02, 0.11);
    
  canvasAllppreferences->SetLogx();
  canvasAllppreferences->SetLogy();
  TH2F * histo2DAllppreferences;
  histo2DAllppreferences = new TH2F("histo2DAllppreferences","histo2DAllppreferences",1000,.3,30.,1000,2e2,10e11);
  SetStyleHistoTH2ForGraphs(histo2DAllppreferences, "#it{p}_{T} (GeV/#it{c})","#it{E}#frac{d^{2}#sigma}{d#it{p}^{3}}(pb GeV^{-2} #it{c}^{3})", 0.05,0.06, 0.05,0.05, 0.7,1., 512, 505);
  histo2DAllppreferences->DrawCopy();
    
  DrawGammaSetMarkerTGraphAsym(graphPCMppreferenceSystErr,20,1.5, 1, 1, 1, kTRUE);
  DrawGammaSetMarkerTGraphAsym(graphPCMppreferenceStatErr,20,1.5, 1, 1);
  DrawGammaSetMarkerTGraphAsym(graphDalitzppreferenceSystErr,20,1.5,  kCyan+2,  kCyan+2, 1, kTRUE);
  DrawGammaSetMarkerTGraphAsym(graphDalitzppreferenceStatErr,20,1.5,  kCyan+2,  kCyan+2);
  DrawGammaSetMarkerTGraphAsym(graphPHOSppreferenceSystErr,20,1.5, kRed+1, kRed+1, 1, kTRUE);
  DrawGammaSetMarkerTGraphAsym(graphPHOSppreferenceStatErr,20,1.5, kRed+1, kRed+1);
  DrawGammaSetMarkerTGraphAsym(graphEMCalppreferenceSystErr,20,1.5, kGreen+2, kGreen+2, 1, kTRUE);
  DrawGammaSetMarkerTGraphAsym(graphEMCalppreferenceStatErr,20,1.5, kGreen+2, kGreen+2);

  graphPCMppreferenceSystErr->Draw("E2,same");  
  graphPCMppreferenceStatErr->Draw("p,same");  
  graphDalitzppreferenceSystErr->Draw("E2,same");  
  graphDalitzppreferenceStatErr->Draw("p,same");  
  graphPHOSppreferenceSystErr->Draw("E2,same");  
  graphPHOSppreferenceStatErr->Draw("p,same");  
  graphEMCalppreferenceSystErr->Draw("E2,same");  
  graphEMCalppreferenceStatErr->Draw("p,same");  

  
  TLegend* legendAllppreferences = new TLegend(0.15,0.15,0.45,0.5);
  legendAllppreferences->SetFillColor(0);
  legendAllppreferences->SetLineColor(0);
  legendAllppreferences->SetTextSize(0.03);
  legendAllppreferences->AddEntry(graphPCMppreferenceSystErr,"PCM pp reference #sqrt{#it{s}} = 5.02 TeV","pef");
  legendAllppreferences->AddEntry(graphDalitzppreferenceSystErr,"Dalitz pp reference #sqrt{#it{s}} = 5.02 TeV","pef");
  legendAllppreferences->AddEntry(graphPHOSppreferenceSystErr,"PHOS pp reference #sqrt{#it{s}} = 5.02 TeV","pef");
  legendAllppreferences->AddEntry(graphEMCalppreferenceSystErr,"EMCal pp reference #sqrt{#it{s}} = 5.02 TeV","pef");
   
  legendAllppreferences->Draw();

  line->Draw("same");
 
    
  canvasAllppreferences->Print(Form("%s/Allppreferences.%s",outputDir.Data(),suffix.Data()));

  // **************************************************************************************
  // ************************* Plotting Ratio all ppreferences ************************
  // **************************************************************************************


 

  TGraphErrors* bla1 = NULL;
  TGraphErrors* bla2 = NULL;
  TGraphErrors* bla3 = NULL;
  TGraphErrors* bla4 = NULL;
  cout << "PCM Spectrum  - Dalitz" << endl;
  TGraphErrors* graphRatioPCMDalitz = CalculateRatioBetweenSpectraWithDifferentBinning(graphDalitzppreferenceStatErr,graphDalitzppreferenceSystErr ,graphPCMppreferenceStatErr,graphPCMppreferenceSystErr,  kTRUE,  kTRUE,&bla1,&bla2,&bla3,&bla4)  ;
  //  graphRatioPCMDalitz->SetPointError(0,0.05,0.34515);

  graphRatioPCMDalitz->Print();
   
  cout << "PCM Spectrum  - PHOS " << endl;
  TGraphErrors* graphRatioPCMPHOS = CalculateRatioBetweenSpectraWithDifferentBinning(graphPHOSppreferenceStatErr,graphPHOSppreferenceSystErr,graphPCMppreferenceStatErr,graphPCMppreferenceSystErr,  kTRUE,  kTRUE,&bla1,&bla2,&bla3,&bla4)  ;
  graphRatioPCMPHOS->Print();
                                                                                
  cout << "PCM Spectrum  - EMCal " << endl;
  TGraphErrors* graphRatioPCMEMCal = CalculateRatioBetweenSpectraWithDifferentBinning(graphEMCalppreferenceStatErr,graphEMCalppreferenceSystErr,graphPCMppreferenceStatErr,graphPCMppreferenceSystErr,  kTRUE,  kTRUE,&bla1,&bla2,&bla3,&bla4)  ;








  TCanvas* canvasRatioAllppreferences = new TCanvas("canvasRatioAllppreferences","",200,10,1200,700);  // gives the page size
  DrawGammaCanvasSettings(canvasRatioAllppreferences , 0.1, 0.02, 0.02, 0.11);
    
  canvasRatioAllppreferences->SetLogx();
  //  canvasRatioAllppreferences->SetLogy();
  TH2F * histo2DRatioAllppreferences;
  histo2DRatioAllppreferences = new TH2F("histo2DRatioAllppreferences","histo2DRatioAllppreferences",1000,.3,30.,1000,0.5,1.5);
  SetStyleHistoTH2ForGraphs(histo2DRatioAllppreferences, "#it{p}_{T} (GeV/#it{c})","ratio pp references", 0.05,0.06, 0.05,0.05, 0.7,.8, 512, 505);
  histo2DRatioAllppreferences->DrawCopy();
    

  DrawGammaSetMarkerTGraphErr(graphRatioPCMDalitz,20,1.5,  kCyan+2,  kCyan+2);
  DrawGammaSetMarkerTGraphErr(graphRatioPCMPHOS,20,1.5, kRed+1, kRed+1);
  DrawGammaSetMarkerTGraphErr(graphRatioPCMEMCal,20,1.5, kGreen+2, kGreen+2);


  graphRatioPCMDalitz->Draw("p,same");  
  graphRatioPCMPHOS->Draw("p,same");  
  graphRatioPCMEMCal->Draw("p,same");  

  
  TLegend* legendRatioAllppreferences = new TLegend(0.15,0.15,0.45,0.35);
  legendRatioAllppreferences->SetFillColor(0);
  legendRatioAllppreferences->SetLineColor(0);
  legendRatioAllppreferences->SetTextSize(0.03);
  legendRatioAllppreferences->AddEntry(graphRatioPCMDalitz,"pp reference #sqrt{#it{s}} = 5.02 TeV, Dalitz/PCM ","pe");
  legendRatioAllppreferences->AddEntry(graphRatioPCMPHOS,"pp reference #sqrt{#it{s}} = 5.02 TeV, PHOS/PCM","pe");
  legendRatioAllppreferences->AddEntry(graphRatioPCMEMCal,"pp reference #sqrt{#it{s}} = 5.02 TeV, EMCal/PCM","pe");
   
  legendRatioAllppreferences->Draw();

  line->Draw("same");
 
    
  canvasRatioAllppreferences->Print(Form("%s/RatioAllppreferences.%s",outputDir.Data(),suffix.Data()));



  // **************************************************************************************
  // ************************* Plotting all Alpha ************************
  // **************************************************************************************
//   TCanvas* canvasAllAlpha = new TCanvas("canvasAllAlpha","",200,10,1200,700);  // gives the page size
//   DrawGammaCanvasSettings(canvasAllAlpha , 0.07, 0.02, 0.02, 0.1);
//     
//   //     canvasAllAlpha->SetLogx();
//   // canvasAllAlpha->SetLogy();
//   TH2F * histo2DAllAlpha;
//   histo2DAllAlpha = new TH2F("histo2DAllAlpha","histo2DAllAlpha",1000,0.,20.,1000,0.2,1.7);
//   SetStyleHistoTH2ForGraphs(histo2DAllAlpha, "#it{p}_{T} (GeV/#it{c})","#alpha", 0.04,0.04, 0.04,0.04, 1.,.6, 512, 505);
//   histo2DAllAlpha->DrawCopy();
//     
//   DrawGammaSetMarkerTGraphAsym(graphPCMAlpha,20,1.5, 1, 1);
//   DrawGammaSetMarkerTGraphAsym(graphDalitzAlpha,20,1.5,  kCyan+2,  kCyan+2);
//   DrawGammaSetMarkerTGraphAsym(graphPHOSAlpha,20,1.5, kRed+1, kRed+1);
//   DrawGammaSetMarkerTGraphAsym(graphEMCalAlpha,20,1.5, kGreen+2, kGreen+2);
// 
// 
//   graphPCMAlpha->Draw("p,same");
//   graphDalitzAlpha->Draw("p,same");  
//   graphPHOSAlpha->Draw("p,same");  
//   graphEMCalAlpha->Draw("p,same");  
// 
//   
//   TLegend* legendAllAlpha = new TLegend(0.65,0.15,0.95,0.35);
//   legendAllAlpha->SetFillColor(0);
//   legendAllAlpha->SetLineColor(0);
//   legendAllAlpha->SetTextSize(0.03);
//   legendAllAlpha->AddEntry(graphPCMAlpha,"PCM","pe");
//   legendAllAlpha->AddEntry(graphDalitzAlpha,"Dalitz","pe");
//   legendAllAlpha->AddEntry(graphPHOSAlpha,"PHOS","pe");
//   legendAllAlpha->AddEntry(graphEMCalAlpha,"EMCal","pe");
//    
//   legendAllAlpha->Draw();
//  
//     
//   canvasAllAlpha->Print(Form("%s/AllAlpha.%s",outputDir.Data(),suffix.Data()));

  // **************************************************************************************
  // ************************* Plotting  R_pPb from Combined pPb spectrum *****************
  // **************************************************************************************

  if (IsNSD){

    graphCombYieldPi0pPbSystErr=ScaleGraph(graphCombYieldPi0pPbSystErr,Scaling);
    graphCombYieldPi0pPb=ScaleGraph(graphCombYieldPi0pPb,Scaling);
    
    //  graphCombYieldPi0pPbSystErr=ApplyNSDSysError(graphCombYieldPi0pPbSystErr,ScalingErr,0.);  
    
  }



  TCanvas* canvasCombRpPb2 = new TCanvas("canvasCombRpPb2","",200,10,1200,700);  // gives the page size
  DrawGammaCanvasSettings( canvasCombRpPb2, 0.08, 0.02, 0.02, 0.13);
    
  //  canvasCombRpPb2->SetLogx();
  // canvasCombRpPb2->SetLogy();
  TH2F * histo2DCombined2;
  histo2DCombined2 = new TH2F("histo2DCombined2","histo2DCombined2",1000,0.,20.,1000,0.3,1.7);
  SetStyleHistoTH2ForGraphs(histo2DCombined2, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#pi^{0}}_{p-Pb}", 0.05,0.06, 0.05,0.05, 0.9,0.7, 512, 505);
  histo2DCombined2->DrawCopy();
    

     
    

  DrawGammaSetMarkerTGraphAsym(graphCombYieldPi0pPbSystErr,20,1.5, 4, 4, 1, kTRUE);  
      
  graphCombYieldPi0pPbSystErr->Draw("E2,same");

  DrawGammaSetMarkerTGraphAsym(graphCombYieldPi0pPb,20,1.5, kBlue, kBlue);  
  graphCombYieldPi0pPb->Draw("p,same");


  TLatex * lt3 = new TLatex(1.,1.55,"ALICE Preliminary"); 
  lt3->SetTextSize(0.05) ;
  //  lt3->DrawText(1.,1.45,"2016/06/14") ;
  lt3->Draw("same");
     
  TLegend* legendRpPbCombine2 = new TLegend(0.1,0.75,0.55,0.85);
  legendRpPbCombine2->SetFillColor(0);
  legendRpPbCombine2->SetLineColor(0);
  legendRpPbCombine2->SetTextSize(0.04);
  legendRpPbCombine2->AddEntry(graphCombYieldPi0pPbSystErr,"p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");

    
  legendRpPbCombine->Draw();
  
  TLine* line2 =new TLine(0.,1.,20.,1.);
  line2->Draw("same");

  BoxNorm->Draw("same");

  canvasCombRpPb2->Print(Form("%s/RpPb_CombinedpPbSpectrum.%s",outputDir.Data(),suffix.Data()));
    
    

  //************************ WriteResultsToFile **************************************************
  TFile fResults(Form("%s/ResultsRpPbpPb_%s.root",outputDir.Data(), dateForOutput.Data()),"RECREATE");
     
  graphInvYieldPi0CombpPb5023GeVSysClone->Write("CombinedPi0RpPbSystErr");
  graphInvYieldPi0CombpPb5023GeVStaClone->Write("CombinedPi0RpPbStatErr");
  
  
  fResults.Close();

  TFile fPaperPlots(Form("%s/PaperPlotsRpPb_%s.root",outputDir.Data(), dateForOutput.Data()),"RECREATE");
     
  graphInvYieldPi0CombpPb5023GeVSysClone->Write("CombinedPi0RpPbSystErr");
  graphInvYieldPi0CombpPb5023GeVStaClone->Write("CombinedPi0RpPbStatErr");
  graphPi0ESP09sPi0KKP->Write("EPS09s_KKP_NLO");
  graphPi0ESP09sPi0AKK->Write("EPS09s_AKK_NLO");
  graphPi0DSS5000->Write("EPS09s_fDSS_NLO");
  graphAsymmErrorsPi0DSS5000->Write("EPS09s_fDSS_errors");  
  graphPi0CGC->Write("CGC");
  graphCombppreferenceStatErr->Write("ppRefStat");
  graphCombppreferenceSystErr->Write("ppRefSyst");


  fPaperPlots.Close();



}
    
