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
const Int_t  Ntotal = 31;//31
const Int_t  nPtLimits = Ntotal+1;

const Int_t  NtotalEta = 16;
const Int_t  nPtLimitsEta = NtotalEta+1;


 
Style_t markerStylePi0	= 20;
Style_t markerStyleEta	= 21;
Style_t markerStylePCM	= 20;
Style_t markerStylePHOS	= 21;
Style_t markerStyleEMCal 	= 33;
Style_t markerStyleDalitz	= 29;
Style_t markerStylePCMEMCal   = 34;



Color_t colorPCM   		= kBlack;
Color_t colorDalitz     	= kViolet; //kBlue+1;kViolet-4 
Color_t colorPHOS   		= kRed+1;
Color_t colorEMCal   		= kGreen+2;
Color_t colorPCMEMCal   	= kBlue+1;
Color_t colorCombYieldPi0   	= kBlue+2;
Color_t colorCombYieldEta 	= kGreen+3;


 Size_t  markerSizePCM=1.;	
  Size_t  markerSizeDalitz=1.2;	
  Size_t  markerSizePHOS=1.;	
  Size_t  markerSizeEMCal=1.2;
  Size_t  markerSizePCMEMCal=1.2;



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

  Double_t LabelOffsetLog=-0.015;





void CombineRpPb5023GeV(Bool_t IsNSD=kFALSE){    

  TString date                                        = ReturnDateString();
    
  gROOT->Reset();
  gROOT->SetStyle("Plain");
    
  StyleSettingsThesis();
  SetPlotStyle();
    
  TString suffix                                      = "pdf";
  TString FittingType                                 = "Tsallis";
    
  TString dateForOutput                               = ReturnDateStringForOutput();
  TString outputDir                                   = Form("CombineRpPb/%s/%s",suffix.Data(),dateForOutput.Data());
  if (IsNSD) outputDir                                = Form("CombineRpPb/%s/%s_NSD",suffix.Data(),dateForOutput.Data());

  gSystem->Exec("mkdir -p "+outputDir);
    
  Double_t ScalingErr = 0.031; //pure NSD Err
  Double_t Scaling = 0.964;// NSD Err
  Double_t ScalingErrTpPb =   0.0356; //Normalization Err on TpPb

  Double_t NormalizationError =TMath::Sqrt(pow(0.031,2)+pow(0.036,2)+pow(0.036,2));//pPb normalization,TpPb and pp
  
  //normalization errors
  //cout << dateForOutput.Data() << endl;
  //___________________________________ Declaration of files _____________________________________________

  cout<<"Normalization error: "<<NormalizationError<<endl;
  
  TString fileNameResultspPb                  = "CombinepPbSpectra/2017_12_18_NSD/Tsallis/eps/PaperPlots_Tsallis_2017_12_18.root";
  TString fileOutputResultspPb                = Form("%s/PaperPlots_Tsallis_%s_ALL.root",outputDir.Data(),dateForOutput.Data());
  
  gSystem->Exec(Form("cp %s %s",fileNameResultspPb.Data(),fileOutputResultspPb.Data()));
  
  fstream fFits;
  TString nameFits = Form("%s/FittingParametersPPReference_%s.dat",outputDir.Data(),dateForOutput.Data());
  fFits.open(nameFits.Data(), ios::out);
    
 
  /*TString fileNameRpPb                        = "ExternalInputpPb/InputRpPb/Pi0RpPb_PCM_2017_10_12.root";
  TString fileNamePi0RpPbCombined             = "ExternalInputpPb/InputRpPb/Pi0RpPb_Comb_2017_10_24.root";
  TString fileNameEtaRpPbCombined             = "ExternalInputpPb/InputRpPb/EtaRpPb_Comb_2017_10_24.root";
  
  
  /////////////////////////////////////////////////Pi0///////////////////////////////////////////
  TString fileNameRpPbPCM                     = "ExternalInputpPb/InputRpPb/Pi0RpPb_PCM_2017_10_12.root";
  TString fileNameRpPbDalitz                  = "ExternalInputpPb/InputRpPb/Pi0RpPb_Dalitz_2017_10_12.root";
  TString fileNameRpPbPHOS                    = "ExternalInputpPb/InputRpPb/Pi0RpPb_PHOS_2017_10_24.root";
  TString fileNameRpPbEMCal                   = "ExternalInputpPb/InputRpPb/Pi0RpPb_EMCal_2017_10_12.root";
  TString fileNameRpPbPCMEMCal                = "ExternalInputpPb/InputRpPb/Pi0RpPb_PCM-EMCal_2017_10_12.root";
  ///////////////////////////////////////////////Eta/////////////////////////////////////////////
  
  TString fileNameRpPbPCMEta                  = "ExternalInputpPb/InputRpPb/EtaRpPb_PCM_2017_10_12.root";
  TString fileNameRpPbEMCalEta                = "ExternalInputpPb/InputRpPb/EtaRpPb_EMCal_2017_10_12.root";
  TString fileNameRpPbPCMEMCalEta             = "ExternalInputpPb/InputRpPb/EtaRpPb_PCM-EMCal_2017_10_12.root";
  
  */
  
  
  TString fileNameRpPb                        = "ExternalInputpPb/InputRpPb/Pi0RpPb_PCM_2017_12_17.root";
  TString fileNamePi0RpPbCombined             = "ExternalInputpPb/InputRpPb/Pi0RpPb_Comb_2017_12_17.root";
  TString fileNameEtaRpPbCombined             = "ExternalInputpPb/InputRpPb/EtaRpPb_Comb_2017_12_17.root";
  
  
  /////////////////////////////////////////////////Pi0///////////////////////////////////////////
  TString fileNameRpPbPCM                     = "ExternalInputpPb/InputRpPb/Pi0RpPb_PCM_2017_12_17.root";
  TString fileNameRpPbDalitz                  = "ExternalInputpPb/InputRpPb/Pi0RpPb_Dalitz_2017_12_17.root";
  TString fileNameRpPbPHOS                    = "ExternalInputpPb/InputRpPb/Pi0RpPb_PHOS_2017_12_17.root";
  TString fileNameRpPbEMCal                   = "ExternalInputpPb/InputRpPb/Pi0RpPb_EMCal_2017_12_17.root";
  TString fileNameRpPbPCMEMCal                = "ExternalInputpPb/InputRpPb/Pi0RpPb_PCM-EMCal_2017_12_17.root";
  ///////////////////////////////////////////////Eta/////////////////////////////////////////////
  
  TString fileNameRpPbPCMEta                  = "ExternalInputpPb/InputRpPb/EtaRpPb_PCM_2017_12_17.root";
  TString fileNameRpPbEMCalEta                = "ExternalInputpPb/InputRpPb/EtaRpPb_EMCal_2017_12_17.root";
  TString fileNameRpPbPCMEMCalEta             = "ExternalInputpPb/InputRpPb/EtaRpPb_PCM-EMCal_2017_12_17.root";
  
  
  
  
  //TString corrFactorsFileName                 = "eps/2017_06_07/ComputeCorrelationFactors_pPb5TeV/pPb5TeV.root";
   TString corrFactorsFileName                 = "/opt/AnalysisSoftware/pPb5TeV.root";
  
  
  TFile* fileNeutralPionRpPb                           	= new TFile(fileNameRpPb.Data());
  TFile* fileNeutralPionRpPbComb                        = new TFile(fileNamePi0RpPbCombined.Data());
  TFile* fileNeutralPionRpPbPCM                         = new TFile(fileNameRpPbPCM.Data());
  TFile* fileNeutralPionRpPbDalitz                      = new TFile(fileNameRpPbDalitz.Data());
  TFile* fileNeutralPionRpPbPHOS                        = new TFile(fileNameRpPbPHOS.Data());
  TFile* fileNeutralPionRpPbEMCal                       = new TFile(fileNameRpPbEMCal.Data());
  TFile* fileNeutralPionRpPbPCMEMCal		        = new TFile(fileNameRpPbPCMEMCal.Data());
  TFile* fileNeutralCombMesonResults                    = new TFile(fileNameResultspPb.Data());
  TDirectory* directoryPionComb                         = (TDirectory*) fileNeutralCombMesonResults->Get("Pi0pPb");
  TDirectory* directoryEtaComb                          = (TDirectory*) fileNeutralCombMesonResults->Get("EtapPb");
  
  
  
  
  
  
  
  
    
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
  TString nameHistoPi0RpPbCombSystErr                = "Pi0_RpPb_Comb_SystErr";
  TString nameHistoPi0RpPbCombStatErr                = "Pi0_RpPb_Comb_StatErr";
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
  TString nameHistoCombPi0ppReferenceStat            = "Pi0_pp_reference_CombBinning_StatErr";
  TString nameHistoCombPi0ppReferenceSyst            = "Pi0_pp_reference_CombBinning_SystErr";
  TString nameHistoPCMAlpha                          = "Pi0_RpPb_PCM_Alpha";
  TString nameHistoDalitzAlpha                       = "Pi0_RpPb_Dalitz_Alpha";
  TString nameHistoPHOSAlpha                         = "Pi0_RpPb_PHOS_Alpha";
  TString nameHistoEMCalAlpha                        = "Pi0_RpPb_EMCal_Alpha";
  TString nameHistoPCMEMCalAlpha		     = "Pi0_RpPb_PCM-EMCal_Alpha";
  TString nameHistoCombAlpha                         = "Pi0_RpPb_Comb_Alpha";
  

  
  
  
  
  TString nameGraphInvYieldPi0PCMpPb5023GeVSystErryShifted = "graphInvYieldPi0PCMpPb5023TeVSystErr_yShifted";
  TString nameGraphInvYieldPi0PCMpPb5023GeVStatErryShifted = "graphInvYieldPi0PCMpPb5023TeVStatErr_yShifted";
  TString nameGraphInvYieldPi0PHOSpPb5023GeVSystErryShifted = "graphInvYieldPi0PHOSpPb5023TeVSystErr_yShifted";
  TString nameGraphInvYieldPi0PHOSpPb5023GeVStatErryShifted = "graphInvYieldPi0PHOSpPb5023TeVStatErr_yShifted";
  TString nameGraphInvYieldPi0EMCalpPb5023GeVSystErryShifted = "graphInvYieldPi0EMCalpPb5023TeVSystErr_yShifted";
  TString nameGraphInvYieldPi0EMCalpPb5023GeVStatErryShifted = "graphInvYieldPi0EMCalpPb5023TeVStatErr_yShifted";
  TString nameGraphInvYieldPi0PCMEMCalpPb5023GeVSystErryShifted = "graphInvYieldPi0PCM-EMCalpPb5023TeVSystErr_yShifted";
  TString nameGraphInvYieldPi0PCMEMCalpPb5023GeVStatErryShifted = "graphInvYieldPi0PCM-EMCalpPb5023TeVStatErr_yShifted";
  TString nameGraphInvYieldPi0DalitzpPb5023GeVSystErryShifted = "graphInvYieldPi0DalitzpPb5023TeVSystErr_yShifted";
  TString nameGraphInvYieldPi0DalitzpPb5023GeVStatErryShifted = "graphInvYieldPi0DalitzpPb5023TeVStatErr_yShifted";
  
  
 
 
  
  
  TFile* fileEtaRpPbPCM                         = new TFile(fileNameRpPbPCMEta.Data());
  TFile* fileEtaRpPbEMCal                       = new TFile(fileNameRpPbEMCalEta.Data());
  TFile* fileEtaRpPbPCMEMCal		        = new TFile(fileNameRpPbPCMEMCalEta.Data());
  TFile* fileEtaRpPbComb                        = new TFile(fileNameEtaRpPbCombined.Data());
 
   
  TString nameHistoPCMEtaRpPbStat                    = "Eta_RpPb_PCM_StatErr";
  TString nameHistoPCMEtaRpPbSyst                    = "Eta_RpPb_PCM_SystErr";
  TString nameHistoEMCalEtaRpPbStat                  = "Eta_RpPb_EMCal_StatErr";
  TString nameHistoEMCalEtaRpPbSyst                  = "Eta_RpPb_EMCal_SystErr";
  TString nameHistoPCMEMCalEtaRpPbStat		     = "Eta_RpPb_PCM-EMCal_StatErr";
  TString nameHistoPCMEMCalEtaRpPbSyst		     = "Eta_RpPb_PCM-EMCal_SystErr";
  TString nameHistoEtaRpPbCombSystErr                = "Eta_RpPb_Comb_SystErr";
  TString nameHistoEtaRpPbCombStatErr                = "Eta_RpPb_Comb_StatErr";
  
  
  
  
  
  
  TString nameHistoPCMEtappReferenceStat             = "Eta_pp_reference_PCMBinning_StatErr";
  TString nameHistoPCMEtappReferenceSyst             = "Eta_pp_reference_PCMBinning_SystErr";
  TString nameHistoEMCalEtappReferenceStat           = "Eta_pp_reference_EMCalBinning_StatErr";
  TString nameHistoEMCalEtappReferenceSyst           = "Eta_pp_reference_EMCalBinning_SystErr";
  TString nameHistoPCMEMCalEtappReferenceStat        = "Eta_pp_reference_PCM-EMCalBinning_StatErr";
  TString nameHistoPCMEMCalEtappReferenceSyst	     = "Eta_pp_reference_PCM-EMCalBinning_SystErr";
  TString nameHistoCombEtappReferenceStat            = "Eta_pp_reference_CombBinning_StatErr";
  TString nameHistoCombEtappReferenceSyst            = "Eta_pp_reference_CombBinning_SystErr";
  TString nameHistoCombEtaAlpha                      = "Eta_RpPb_Comb_Alpha";
 
  
  
  TString nameGraphInvYieldEtaPCMpPb5023GeVSystErryShifted = "graphInvYieldEtaPCMpPb5023TeVSystErr_yShifted";
  TString nameGraphInvYieldEtaPCMpPb5023GeVStatErryShifted = "graphInvYieldEtaPCMpPb5023TeVStatErr_yShifted";
  TString nameGraphInvYieldEtaEMCalpPb5023GeVSystErryShifted = "graphInvYieldEtaEMCalpPb5023TeVSystErr_yShifted";
  TString nameGraphInvYieldEtaEMCalpPb5023GeVStatErryShifted = "graphInvYieldEtaEMCalpPb5023TeVStatErr_yShifted";
  TString nameGraphInvYieldEtaPCMEMCalpPb5023GeVSystErryShifted = "graphInvYieldEtaPCM-EMCalpPb5023TeVSystErr_yShifted";
  TString nameGraphInvYieldEtaPCMEMCalpPb5023GeVStatErryShifted = "graphInvYieldEtaPCM-EMCalpPb5023TeVStatErr_yShifted";
  
  
  
  
  
  
    
  TString collisionSystempPb                          = "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV"; 
  
                                                          
  Double_t pTLimits[nPtLimits]                        = { 0.3, 0.4, 0.5, 0.6, 0.7,
                                                          0.8, 1.0, 1.2, 1.4, 1.6, 
							  1.8, 2.0, 2.2, 2.4, 2.6,
							  2.8, 3.0, 3.2, 3.4, 3.6,
							  3.8, 4.0, 4.5, 5.0, 5.5,
							  6.0, 7.0, 8.0, 10.0,12.0,
							  16.0, 20.0};
                                                          
 
   
  Double_t pTLimitsEta[nPtLimitsEta]                        = { 0.7,0.9,1.1,1.4,1.8,
							     2.2,2.6,3.0,3.5,4.0,
							     5.0,6.0,8.0,10.0,12.0,
							     16.0,20.0};
                                                          
                                                          
                                                        
                                
    
  Int_t offSets[11]                                   =  { 0, 6,  8, 0,  5, 3, 0, 0, 0,  0, 0};
  Int_t offSetsSys[11]                                =  { 0, 6,  8, 0,  5, 3, 0, 0, 0, 0, 0};
  
  Int_t offSetsEta[11]                                =  { 0, 0, 4, 0,  2, 0, 0, 0, 0,  0, 0};
  Int_t offSetsSysEta[11]                             =  { 0, 0, 4, 0,  2, 0, 0, 0, 0, 0, 0};
  
  
  
    
    
    
    
  TH1D* statErrorCollection[11];
  TH1D* statErrorCollectionPi0PPRef[11];
  TH1D* statErrorCollectionEtaPPRef[11];
  TH1D* statErrorCollectionEtaRpPb[11];
  for (Int_t i = 0; i< 11; i++){
    statErrorCollection[i]                          = NULL;
    statErrorCollectionPi0PPRef[i]                  = NULL;
    statErrorCollectionEtaRpPb[i]                   = NULL;
    statErrorCollectionEtaPPRef[i]                  = NULL;
  }    
    
  TGraphAsymmErrors* sysErrorCollection[11];
  TGraphAsymmErrors* sysErrorCollectionPi0PPRef[11];
  TGraphAsymmErrors* sysErrorCollectionEtaPPRef[11];
  TGraphAsymmErrors* sysErrorCollectionEtaRpPb[11];
  
  for (Int_t i = 0; i< 11; i++){
    sysErrorCollection[i]                           = NULL;
    sysErrorCollectionPi0PPRef[i]                   = NULL;
    sysErrorCollectionEtaRpPb[i]                    = NULL;
    sysErrorCollectionEtaPPRef[i]                   = NULL;
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
  TH1D* histoDalitzppreferenceStatErr                 = (TH1D*)GraphAsymErrorsToHist_withErrors(graphDalitzppreferenceStatErr,"histoDalitzppreferenceStatErr");
  TGraphAsymmErrors* graphDalitzAlpha                 = (TGraphAsymmErrors*)fileNeutralPionRpPbDalitz->Get( nameHistoDalitzAlpha.Data());
  
  TGraphAsymmErrors* graphInvYieldPi0DalitzpPb5023GeVStatErryShifted = (TGraphAsymmErrors*) fileNeutralPionRpPbDalitz->Get(nameGraphInvYieldPi0DalitzpPb5023GeVStatErryShifted.Data());  
   TGraphAsymmErrors* graphInvYieldPi0DalitzpPb5023GeVSystErryShifted = (TGraphAsymmErrors*) fileNeutralPionRpPbDalitz->Get(nameGraphInvYieldPi0DalitzpPb5023GeVSystErryShifted.Data());  
  

  // **************************************************************************************
  // ****************************** Reading PCM *******************************************
  // **************************************************************************************    
  TGraphAsymmErrors* graphPCMYieldPi0pPb           = (TGraphAsymmErrors*)fileNeutralPionRpPbPCM->Get(nameHistoPCM.Data());
  TH1D* histoPCMYieldPi0pPb                        = (TH1D*)GraphAsymErrorsToHist_withErrors(graphPCMYieldPi0pPb,nameHistoPCM.Data());
  TGraphAsymmErrors* graphPCMYieldPi0pPbSystErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbPCM->Get(nameHistoPCMSysErrors.Data());
  TGraphAsymmErrors* graphPCMppreferenceStatErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbPCM->Get(nameHistoPCMppReferenceStat.Data());
  TH1D*              histoPCMppreferenceStatErr    = (TH1D*)GraphAsymErrorsToHist_withErrors(graphPCMppreferenceStatErr,"histoPCMppreferenceStatErr");
  TGraphAsymmErrors* graphPCMppreferenceSystErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbPCM->Get(nameHistoPCMppReferenceSyst.Data());
  TGraphAsymmErrors* graphPCMAlpha                 = (TGraphAsymmErrors*)fileNeutralPionRpPbPCM->Get( nameHistoPCMAlpha.Data()); 

   TGraphAsymmErrors* graphInvYieldPi0PCMpPb5023GeVStatErryShifted = (TGraphAsymmErrors*) fileNeutralPionRpPbPCM->Get(nameGraphInvYieldPi0PCMpPb5023GeVStatErryShifted.Data());  
   TGraphAsymmErrors* graphInvYieldPi0PCMpPb5023GeVSystErryShifted = (TGraphAsymmErrors*) fileNeutralPionRpPbPCM->Get(nameGraphInvYieldPi0PCMpPb5023GeVSystErryShifted.Data());  
  
  
  
  
  TGraphAsymmErrors* graphPCMEtaRpPbStat           = (TGraphAsymmErrors*)fileEtaRpPbPCM->Get(nameHistoPCMEtaRpPbStat.Data());
  TH1D* histoPCMEtaRpPbStat                        = (TH1D*)GraphAsymErrorsToHist_withErrors(graphPCMEtaRpPbStat,nameHistoPCMEtaRpPbStat.Data());
  TGraphAsymmErrors* graphPCMEtaRpPbSyst    = (TGraphAsymmErrors*)fileEtaRpPbPCM->Get(nameHistoPCMEtaRpPbSyst.Data());
  
  TGraphAsymmErrors* graphPCMEtappreferenceStatErr    = (TGraphAsymmErrors*)fileEtaRpPbPCM->Get(nameHistoPCMEtappReferenceStat.Data());
  TH1D*              histoPCMEtappreferenceStatErr    = (TH1D*) GraphAsymErrorsToHist_withErrors(graphPCMEtappreferenceStatErr,"histoPCMEtappreferenceStatErr");
  
  TGraphAsymmErrors* graphPCMEtappreferenceSystErr    = (TGraphAsymmErrors*)fileEtaRpPbPCM->Get(nameHistoPCMEtappReferenceSyst.Data());
 
   TGraphAsymmErrors* graphInvYieldEtaPCMpPb5023GeVStatErryShifted = (TGraphAsymmErrors*) fileEtaRpPbPCM->Get(nameGraphInvYieldEtaPCMpPb5023GeVStatErryShifted.Data());  
   TGraphAsymmErrors* graphInvYieldEtaPCMpPb5023GeVSystErryShifted = (TGraphAsymmErrors*) fileEtaRpPbPCM->Get(nameGraphInvYieldEtaPCMpPb5023GeVSystErryShifted.Data()); 
  
  

  // **************************************************************************************
  // ******************************* Reading PHOS *****************************************
  // **************************************************************************************    
  
  TGraphAsymmErrors* graphPHOSYieldPi0pPb                      = (TGraphAsymmErrors*)fileNeutralPionRpPbPHOS->Get(nameHistoPHOS.Data()); 
  TH1D* histoPHOSYieldPi0pPb           = GraphAsymErrorsToHist_withErrors(graphPHOSYieldPi0pPb,nameHistoPHOS.Data());
  TGraphAsymmErrors* graphPHOSYieldPi0pPbSystErr                   = (TGraphAsymmErrors*)fileNeutralPionRpPbPHOS->Get(nameHistoPHOSSysErrors.Data());
  TGraphAsymmErrors* graphPHOSppreferenceStatErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbPHOS->Get(nameHistoPHOSppReferenceStat.Data());
  TGraphAsymmErrors* graphPHOSppreferenceSystErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbPHOS->Get(nameHistoPHOSppReferenceSyst.Data()); 
  TH1D* histoPHOSppreferenceStatErr                 = (TH1D*)GraphAsymErrorsToHist_withErrors(graphPHOSppreferenceStatErr,"histoPHOSppreferenceStatErr");
  
  TGraphAsymmErrors* graphPHOSAlpha                 = (TGraphAsymmErrors*)fileNeutralPionRpPbPHOS->Get( nameHistoPHOSAlpha.Data());
  
    TGraphAsymmErrors* graphInvYieldPi0PHOSpPb5023GeVStatErryShifted = (TGraphAsymmErrors*) fileNeutralPionRpPbPHOS->Get(nameGraphInvYieldPi0PHOSpPb5023GeVStatErryShifted.Data());  
   TGraphAsymmErrors* graphInvYieldPi0PHOSpPb5023GeVSystErryShifted = (TGraphAsymmErrors*) fileNeutralPionRpPbPHOS->Get(nameGraphInvYieldPi0PHOSpPb5023GeVSystErryShifted.Data());  
  // **************************************************************************************
  // ******************************** Reading EMCal ***************************************
  // **************************************************************************************
    
  TGraphAsymmErrors* graphEMCalYieldPi0pPb       = (TGraphAsymmErrors*)fileNeutralPionRpPbEMCal->Get(nameHistoEMCal.Data());
  TH1D* histoEMCalYieldPi0pPb           = GraphAsymErrorsToHist_withErrors(graphEMCalYieldPi0pPb,nameHistoEMCal.Data());
  TGraphAsymmErrors* graphEMCalYieldPi0pPbSystErr     = (TGraphAsymmErrors*)fileNeutralPionRpPbEMCal->Get(nameHistoEMCalSysErrors.Data());  
  TGraphAsymmErrors* graphEMCalppreferenceStatErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbEMCal->Get(nameHistoEMCalppReferenceStat.Data());
  TGraphAsymmErrors* graphEMCalppreferenceSystErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbEMCal->Get(nameHistoEMCalppReferenceSyst.Data()); 
  TH1D* histoEMCalppreferenceStatErr                 = (TH1D*)GraphAsymErrorsToHist_withErrors(graphEMCalppreferenceStatErr,"histoEMCalppreferenceStatErr");
  
  TGraphAsymmErrors* graphEMCalAlpha                 = (TGraphAsymmErrors*)fileNeutralPionRpPbEMCal->Get( nameHistoEMCalAlpha.Data());
  
    TGraphAsymmErrors* graphInvYieldPi0EMCalpPb5023GeVStatErryShifted = (TGraphAsymmErrors*) fileNeutralPionRpPbEMCal->Get(nameGraphInvYieldPi0EMCalpPb5023GeVStatErryShifted.Data());  
   TGraphAsymmErrors* graphInvYieldPi0EMCalpPb5023GeVSystErryShifted = (TGraphAsymmErrors*) fileNeutralPionRpPbEMCal->Get(nameGraphInvYieldPi0EMCalpPb5023GeVSystErryShifted.Data());  
  
  //cout <<"stat"<< endl;
  //graphEMCalYieldPi0pPb->Print();
  //cout <<"sys"<< endl;
  //graphEMCalYieldPi0pPbSystErr->Print();
  
  
  TGraphAsymmErrors* graphEMCalEtaRpPbStat                 = (TGraphAsymmErrors*)fileEtaRpPbEMCal->Get(nameHistoEMCalEtaRpPbStat.Data());
  TH1D* histoEMCalEtaRpPbStat                              = GraphAsymErrorsToHist_withErrors(graphEMCalEtaRpPbStat,nameHistoEMCalEtaRpPbStat.Data());
  TGraphAsymmErrors* graphEMCalEtaRpPbSyst                 = (TGraphAsymmErrors*)fileEtaRpPbEMCal->Get(nameHistoEMCalEtaRpPbSyst.Data());
  
  
  TGraphAsymmErrors* graphEMCalEtappreferenceStatErr          = (TGraphAsymmErrors*)fileEtaRpPbEMCal->Get(nameHistoEMCalEtappReferenceStat);
  TH1D*              histoEMCalEtappreferenceStatErr          = GraphAsymErrorsToHist_withErrors(graphEMCalEtappreferenceStatErr,"histoEMCalEtappreferenceStatErr");
  
  
  TGraphAsymmErrors* graphEMCalEtappreferenceSystErr          = (TGraphAsymmErrors*)fileEtaRpPbEMCal->Get(nameHistoEMCalEtappReferenceSyst.Data());
  
   TGraphAsymmErrors* graphInvYieldEtaEMCalpPb5023GeVStatErryShifted = (TGraphAsymmErrors*) fileEtaRpPbEMCal->Get(nameGraphInvYieldEtaEMCalpPb5023GeVStatErryShifted.Data());  
   TGraphAsymmErrors* graphInvYieldEtaEMCalpPb5023GeVSystErryShifted = (TGraphAsymmErrors*) fileEtaRpPbEMCal->Get(nameGraphInvYieldEtaEMCalpPb5023GeVSystErryShifted.Data()); 
  
  
  
  
  
  // ******************************** Reading PCM-EMCal ***************************************
  // **************************************************************************************
  TGraphAsymmErrors* graphPCMEMCalYieldPi0pPb       	= (TGraphAsymmErrors*)fileNeutralPionRpPbPCMEMCal->Get(nameHistoPCMEMCal.Data());
  TH1D* histoPCMEMCalYieldPi0pPb          		= GraphAsymErrorsToHist_withErrors(graphPCMEMCalYieldPi0pPb,nameHistoPCMEMCal.Data());
  TGraphAsymmErrors* graphPCMEMCalYieldPi0pPbSystErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbPCMEMCal->Get(nameHistoPCMEMCalSysErrors.Data());  
  TGraphAsymmErrors* graphPCMEMCalppreferenceStatErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbPCMEMCal->Get(nameHistoPCMEMCalppReferenceStat.Data());
  TGraphAsymmErrors* graphPCMEMCalppreferenceSystErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbPCMEMCal->Get(nameHistoPCMEMCalppReferenceSyst.Data()); 
  TH1D*              histoPCMEMCalppreferenceStatErr    = (TH1D*) GraphAsymErrorsToHist_withErrors(graphPCMEMCalppreferenceStatErr,"histoPCMEMCalppreferenceStatErr");
  
  TGraphAsymmErrors* graphPCMEMCalAlpha                 = (TGraphAsymmErrors*)fileNeutralPionRpPbPCMEMCal->Get(nameHistoPCMEMCalAlpha.Data());
  cout <<"stat"<< endl;
  graphPCMEMCalYieldPi0pPb->Print();
  cout <<"sys"<< endl;
  graphPCMEMCalYieldPi0pPbSystErr->Print();
  
  TGraphAsymmErrors* graphInvYieldPi0PCMEMCalpPb5023GeVStatErryShifted = (TGraphAsymmErrors*) fileNeutralPionRpPbPCMEMCal->Get(nameGraphInvYieldPi0PCMEMCalpPb5023GeVStatErryShifted.Data());  
   TGraphAsymmErrors* graphInvYieldPi0PCMEMCalpPb5023GeVSystErryShifted = (TGraphAsymmErrors*) fileNeutralPionRpPbPCMEMCal->Get(nameGraphInvYieldPi0PCMEMCalpPb5023GeVSystErryShifted.Data());  
  
  
  
  TGraphAsymmErrors* graphPCMEMCalEtaRpPbStat       	= (TGraphAsymmErrors*)fileEtaRpPbPCMEMCal->Get(nameHistoPCMEMCalEtaRpPbStat.Data());
  TH1D* histoPCMEMCalEtaRpPbStat          		= GraphAsymErrorsToHist_withErrors(graphPCMEMCalEtaRpPbStat,nameHistoPCMEMCalEtaRpPbStat.Data());
  
  
  TGraphAsymmErrors* graphPCMEMCalEtaRpPbSyst    = (TGraphAsymmErrors*)fileEtaRpPbPCMEMCal->Get(nameHistoPCMEMCalEtaRpPbSyst.Data());  
  

  
  TGraphAsymmErrors* graphPCMEMCalEtappreferenceStatErr    = (TGraphAsymmErrors*)fileEtaRpPbPCMEMCal->Get(nameHistoPCMEMCalEtappReferenceStat.Data());
  TH1D* histoPCMEMCalEtappreferenceStatErr  = (TH1D*)GraphAsymErrorsToHist_withErrors(graphPCMEMCalEtappreferenceStatErr, "histoPCMEMCalEtappreferenceStatErr");
  TGraphAsymmErrors* graphPCMEMCalEtappreferenceSystErr    = (TGraphAsymmErrors*)fileEtaRpPbPCMEMCal->Get(nameHistoPCMEMCalEtappReferenceSyst.Data()); 
  
  cout<<"//////////////////////////////////////"<<endl;
  graphPCMEMCalEtappreferenceStatErr->Print();
  cout<<"//////////////////////////////////////"<<endl;
  graphPCMEMCalEtappreferenceSystErr->Print();
  
   TGraphAsymmErrors* graphInvYieldEtaPCMEMCalpPb5023GeVStatErryShifted = (TGraphAsymmErrors*) fileEtaRpPbPCMEMCal->Get(nameGraphInvYieldEtaPCMEMCalpPb5023GeVStatErryShifted.Data());  
   TGraphAsymmErrors* graphInvYieldEtaPCMEMCalpPb5023GeVSystErryShifted = (TGraphAsymmErrors*) fileEtaRpPbPCMEMCal->Get(nameGraphInvYieldEtaPCMEMCalpPb5023GeVSystErryShifted.Data());  
   
 
  // **************************************************************************************
  // ******************************** Reading RpPb from Combined pPb spectrum *************
  // **************************************************************************************
    
  TGraphAsymmErrors* graphCombPi0RpPbStatErr        = (TGraphAsymmErrors*)fileNeutralPionRpPbComb->Get(nameHistoPi0RpPbCombStatErr.Data());
  TGraphAsymmErrors* graphCombPi0RpPbSystErr        = (TGraphAsymmErrors*)fileNeutralPionRpPbComb->Get(nameHistoPi0RpPbCombSystErr.Data());  
  TGraphAsymmErrors* graphCombppreferenceStatErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbComb->Get(nameHistoCombPi0ppReferenceStat.Data());
  TGraphAsymmErrors* graphCombppreferenceSystErr    = (TGraphAsymmErrors*)fileNeutralPionRpPbComb->Get(nameHistoCombPi0ppReferenceSyst.Data()); 
  TGraphAsymmErrors* graphCombAlpha                 = (TGraphAsymmErrors*)fileNeutralPionRpPbComb->Get( nameHistoCombAlpha.Data());
  
  
  
     
  TGraphAsymmErrors* graphCombEtaRpPbStatErr          = (TGraphAsymmErrors*)fileEtaRpPbComb->Get(nameHistoEtaRpPbCombStatErr.Data());
  TGraphAsymmErrors* graphCombEtaRpPbSystErr          = (TGraphAsymmErrors*)fileEtaRpPbComb->Get(nameHistoEtaRpPbCombSystErr.Data());  
  TGraphAsymmErrors* graphCombEtappreferenceStatErr   = (TGraphAsymmErrors*)fileEtaRpPbComb->Get(nameHistoCombEtappReferenceStat.Data());
  TGraphAsymmErrors* graphCombEtappreferenceSystErr   = (TGraphAsymmErrors*)fileEtaRpPbComb->Get(nameHistoCombEtappReferenceSyst.Data()); 
  TGraphAsymmErrors* graphCombEtaAlpha                = (TGraphAsymmErrors*)fileEtaRpPbComb->Get( nameHistoCombEtaAlpha.Data());
   
   
  //**********************************************************************************************
  //**********************************Reading combined pi0 and eta spectra************************
  //**********************************************************************************************
  
  
  
  
  TGraphAsymmErrors* graphAsymmCombPi0pPbSystErr = (TGraphAsymmErrors*) directoryPionComb->Get("graphInvYieldPi0CombpPb5023GeVSystErr_NSD");
  TGraphAsymmErrors* graphAsymmCombPi0pPbStatErr = (TGraphAsymmErrors*) directoryPionComb->Get("graphInvYieldPi0CombpPb5023GeVStatErr_NSD");
  
  
  TGraphAsymmErrors* graphAsymmCombEtapPbSystErr = (TGraphAsymmErrors*) directoryEtaComb->Get("graphInvYieldEtaCombpPb5023GeVSystErr_NSD");
  TGraphAsymmErrors* graphAsymmCombEtapPbStatErr = (TGraphAsymmErrors*) directoryEtaComb->Get("graphInvYieldEtaCombpPb5023GeVStatErr_NSD");
  
  
  
  
  
  
  
  
  
  
  

  // **************************************************************************************
  // ******************************** Reading Model calculations **************************
  // **************************************************************************************
  TGraphAsymmErrors*  graphAsymmErrorsPi0DSS5000    = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get("EPS09s_fDSS_errors"); 
  TGraph*  graphPi0CGC                              = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get("ColorGlasCondensate"); 
  TGraph*  graphPi0DSS5000                          = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get("EPS09s_fDSS_NLO"); 
  TGraph*  graphPi0ESP09sPi0AKK                     = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get("EPS09s_AKK_NLO"); 
  TGraph*  graphPi0ESP09sPi0KKP                     = (TGraphAsymmErrors*)fileNeutralPionRpPb->Get("EPS09s_KKP_NLO"); 


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
												      fileNameOutputWeightingOld,"pPb_5.023TeV","Pi0RpPb",kTRUE,sysErrorPPReference,corrFactorsFileName.Data()
												      ); 
cout<<"\n\n */////////////////Combination of Pi0 RpPb  ends here**//////////////"<<endl;
 


  



  if (IsNSD){

    graphCombPi0InvCrossSectionStatpPb5023GeV=ScaleGraph(graphCombPi0InvCrossSectionStatpPb5023GeV,Scaling);
    graphCombPi0InvCrossSectionSyspPb5023GeV=ScaleGraph(graphCombPi0InvCrossSectionSyspPb5023GeV,Scaling);

    graphCombPi0InvCrossSectionTotpPb5023GeV  =  CalculateCombinedSysAndStatError( graphCombPi0InvCrossSectionStatpPb5023GeV ,graphCombPi0InvCrossSectionSyspPb5023GeV );
    
    cout<< "Combined RpPb NSD scaled:" << endl;
    graphCombPi0InvCrossSectionTotpPb5023GeV->Print();
    graphCombPi0InvCrossSectionStatpPb5023GeV->Print();
    graphCombPi0InvCrossSectionSyspPb5023GeV->Print();
  }
  
   
 Double_t fTpPb              = 0.0983e3*(1/recalcBarn);
 Double_t fTpPbErr           = 0.0035e3*(1/recalcBarn);

TGraphAsymmErrors* graphCombPi0InvCrossSectionStatErrPPRef = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionStatpPb5023GeV->Clone();
TGraphAsymmErrors* graphCombPi0InvCrossSectionSystErrPPRef = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionSyspPb5023GeV->Clone();

  Int_t nPointsPPPi0Ref = graphCombPi0InvCrossSectionStatErrPPRef->GetN();
  
        

  for(Int_t iPoint=0; iPoint < nPointsPPPi0Ref; iPoint++)
  {
      
      
    Double_t value =  graphAsymmCombPi0pPbStatErr->GetY()[iPoint]/(fTpPb*graphCombPi0InvCrossSectionStatErrPPRef->GetY()[iPoint]  );
         
        Double_t  yVal = graphCombPi0InvCrossSectionStatErrPPRef->GetY()[iPoint];
        Double_t  yValPP = graphAsymmCombPi0pPbStatErr->GetY()[iPoint];
        
        Double_t errYStatlow      = graphCombPi0InvCrossSectionStatErrPPRef->GetEYlow()[iPoint];
        Double_t errYStathigh     = graphCombPi0InvCrossSectionStatErrPPRef->GetEYhigh()[iPoint];
        
        Double_t errYStatlowPP    = graphAsymmCombPi0pPbStatErr->GetEYlow()[iPoint];
        Double_t errYStathighPP   = graphAsymmCombPi0pPbStatErr->GetEYhigh()[iPoint];
        
        Double_t errYSystlow      = graphCombPi0InvCrossSectionSystErrPPRef->GetEYlow()[iPoint];
        Double_t errYSysthigh     = graphCombPi0InvCrossSectionSystErrPPRef->GetEYhigh()[iPoint];
        
        Double_t errYSystlowPP    = graphAsymmCombPi0pPbSystErr->GetEYlow()[iPoint];
        Double_t errYSysthighPP   = graphAsymmCombPi0pPbSystErr->GetEYhigh()[iPoint];
        
         
        errYSystlow     = pow( pow( errYSystlow/yVal, 2. ) + pow(  errYSystlowPP/yValPP  ,2. ) ,   0.5) * value;
        errYSysthigh    = pow( pow( errYSysthigh/yVal, 2. ) + pow( errYSysthighPP/yValPP  ,2. ) ,   0.5) * value;  
        
        errYStatlow     = pow( pow( errYStatlow/yVal, 2. ) + pow(  errYStatlowPP/yValPP  ,2. ) ,   0.5) * value;
        errYStathigh    = pow( pow( errYStathigh/yVal, 2. ) + pow( errYStathighPP/yValPP  ,2. ) ,   0.5) * value;  
        
        graphCombPi0InvCrossSectionStatErrPPRef->SetPoint(iPoint,graphCombPi0InvCrossSectionStatErrPPRef->GetX()[iPoint],value);
        graphCombPi0InvCrossSectionStatErrPPRef->SetPointEYhigh(iPoint,errYStathigh);
        graphCombPi0InvCrossSectionStatErrPPRef->SetPointEYlow(iPoint,errYStatlow);
        
        graphCombPi0InvCrossSectionSystErrPPRef->SetPoint(iPoint, graphCombPi0InvCrossSectionSystErrPPRef->GetX()[iPoint], value);
        graphCombPi0InvCrossSectionSystErrPPRef->SetPointEYhigh(iPoint,errYSysthigh);
        graphCombPi0InvCrossSectionSystErrPPRef->SetPointEYlow(iPoint,errYSystlow);
        
        
  }
  

 TGraphAsymmErrors* graphCombPi0InvCrossSectionTotErrPPRef  =  CalculateCombinedSysAndStatError( graphCombPi0InvCrossSectionStatErrPPRef ,graphCombPi0InvCrossSectionSystErrPPRef );
  
  
  
  
  

    
  TGraphAsymmErrors* graphInvYieldPi0CombpPb5023GeVStaClone   = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionStatpPb5023GeV->Clone();
  TGraphAsymmErrors* graphInvYieldPi0CombpPb5023GeVSysClone   = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionSyspPb5023GeV->Clone();
  TGraphAsymmErrors* graphInvYieldPi0CombpPb5023GeVTotClone   = (TGraphAsymmErrors*) graphCombPi0InvCrossSectionTotpPb5023GeV->Clone();
  
  
  
   cout<<"\n\n */////////////////Combination of Eta RpPb  begins here**//////////////"<<endl;
 
  
  
  TString fileNameOutputWeightingEtaRpPb                = Form("%s/Eta_WeightingRpPb.dat",outputDir.Data());
  
  

  statErrorCollectionEtaRpPb[0]          = (TH1D*)histoPCMEtaRpPbStat->Clone("statErrPCMEta");
  statErrorCollectionEtaRpPb[2]          = (TH1D*)histoEMCalEtaRpPbStat->Clone("statErrEMCalPi0");
  statErrorCollectionEtaRpPb[4]	         = (TH1D*)histoPCMEMCalEtaRpPbStat->Clone("statErrPCM-EMCalPi0");
    
  sysErrorCollectionEtaRpPb[0]           = (TGraphAsymmErrors*)graphPCMEtaRpPbSyst->Clone("sysErrPCMEta");
  sysErrorCollectionEtaRpPb[2]           = (TGraphAsymmErrors*)graphEMCalEtaRpPbSyst->Clone("sysErrEMCalEta");
  sysErrorCollectionEtaRpPb[4]		 = (TGraphAsymmErrors*)graphPCMEMCalEtaRpPbSyst->Clone("sysErrPCM-EMCalEta");
  
  TGraphAsymmErrors* graphCombEtaRpPb5023GeVStat = NULL;
  TGraphAsymmErrors* graphCombEtaRpPb5023GeVSyst = NULL;
  
  cout<<"\n\n */////////////////Combination of Eta RpPb  begins here**//////////////"<<endl;
 
    
  TGraphAsymmErrors* graphCombEtaRpPb5023GeVTot = CombinePtPointsSpectraFullCorrMat(   statErrorCollectionEtaRpPb,    sysErrorCollectionEtaRpPb,     
												      pTLimitsEta, NtotalEta,
												      offSetsEta, offSetsSysEta,
												      graphCombEtaRpPb5023GeVStat, graphCombEtaRpPb5023GeVSyst,
                                                                                                      fileNameOutputWeightingEtaRpPb,"pPb_5.023TeV","EtaRpPb",kTRUE,NULL,corrFactorsFileName.Data()
												      ); 
  cout<<"\n\n */////////////////Combination of Pi0 RpPb  ends here**//////////////"<<endl;
 
  graphCombEtaRpPb5023GeVStat->Print();

  if (IsNSD){

    graphCombEtaRpPb5023GeVStat=ScaleGraph(graphCombEtaRpPb5023GeVStat,Scaling);
    graphCombEtaRpPb5023GeVSyst=ScaleGraph(graphCombEtaRpPb5023GeVSyst,Scaling);

    graphCombEtaRpPb5023GeVTot  =  CalculateCombinedSysAndStatError( graphCombEtaRpPb5023GeVStat ,graphCombEtaRpPb5023GeVSyst );
    graphCombEtaRpPb5023GeVTot->Print();
    graphCombEtaRpPb5023GeVStat->Print();
    graphCombEtaRpPb5023GeVSyst->Print();
    
  }

    
  TGraphAsymmErrors* graphCombEtaRpPb5023GeVStatClone   = (TGraphAsymmErrors*) graphCombEtaRpPb5023GeVStat->Clone();
  TGraphAsymmErrors* graphCombEtaRpPb5023GeVSystClone   = (TGraphAsymmErrors*) graphCombEtaRpPb5023GeVSyst->Clone();
  TGraphAsymmErrors* graphCombEtaRpPb5023GeVTotClone    = (TGraphAsymmErrors*) graphCombEtaRpPb5023GeVTot->Clone();
  
  
  TGraphAsymmErrors* graphCombEtaInvCrossSectionStatErrPPRef   = (TGraphAsymmErrors*) graphCombEtaRpPb5023GeVStat->Clone();
  TGraphAsymmErrors* graphCombEtaInvCrossSectionSystErrPPRef   = (TGraphAsymmErrors*) graphCombEtaRpPb5023GeVSyst->Clone();
 
  
  
  Int_t nPointsPPEtaRef = graphCombEtaInvCrossSectionStatErrPPRef->GetN();
  
        

  for(Int_t iPoint=0; iPoint < nPointsPPEtaRef; iPoint++)
  {
      
       // Double_t value =  graphCombEtaInvCrossSectionStatErrPPRef->GetY()[iPoint]/(fTpPb*graphAsymmCombEtapPbStatErr->GetY()[iPoint]  );
        
        Double_t value =  graphAsymmCombEtapPbStatErr->GetY()[iPoint]/(fTpPb*graphCombEtaInvCrossSectionStatErrPPRef->GetY()[iPoint]   );
        
        
        Double_t  yVal = graphCombEtaInvCrossSectionStatErrPPRef->GetY()[iPoint];
        Double_t  yValPP = graphAsymmCombEtapPbStatErr->GetY()[iPoint];
        
        Double_t errYStatlow      = graphCombEtaInvCrossSectionStatErrPPRef->GetEYlow()[iPoint];
        Double_t errYStathigh     = graphCombEtaInvCrossSectionStatErrPPRef->GetEYhigh()[iPoint];
        
        Double_t errYStatlowPP    = graphAsymmCombEtapPbStatErr->GetEYlow()[iPoint];
        Double_t errYStathighPP   = graphAsymmCombEtapPbStatErr->GetEYhigh()[iPoint];
        
        Double_t errYSystlow      = graphCombEtaInvCrossSectionSystErrPPRef->GetEYlow()[iPoint];
        Double_t errYSysthigh     = graphCombEtaInvCrossSectionSystErrPPRef->GetEYhigh()[iPoint];
        
        Double_t errYSystlowPP    = graphAsymmCombEtapPbSystErr->GetEYlow()[iPoint];
        Double_t errYSysthighPP   = graphAsymmCombEtapPbSystErr->GetEYhigh()[iPoint];
        
        
        
        
        
        
        errYSystlow     = pow( pow( errYSystlow/yVal, 2. ) + pow(  errYSystlowPP/yValPP  ,2. ) ,   0.5) * value;
        errYSysthigh    = pow( pow( errYSysthigh/yVal, 2. ) + pow( errYSysthighPP/yValPP  ,2. ) ,   0.5) * value;  
        
        errYStatlow     = pow( pow( errYStatlow/yVal, 2. ) + pow(  errYStatlowPP/yValPP  ,2. ) ,   0.5) * value;
        errYStathigh    = pow( pow( errYStathigh/yVal, 2. ) + pow( errYStathighPP/yValPP  ,2. ) ,   0.5) * value;  
        
        graphCombEtaInvCrossSectionStatErrPPRef->SetPoint(iPoint,graphCombEtaInvCrossSectionStatErrPPRef->GetX()[iPoint], value);
        graphCombEtaInvCrossSectionStatErrPPRef->SetPointEYhigh(iPoint,errYStathigh);
        graphCombEtaInvCrossSectionStatErrPPRef->SetPointEYlow(iPoint,errYStatlow);
        
        graphCombEtaInvCrossSectionSystErrPPRef->SetPoint(iPoint,graphCombEtaInvCrossSectionStatErrPPRef->GetX()[iPoint],value);
        graphCombEtaInvCrossSectionSystErrPPRef->SetPointEYhigh(iPoint,errYSysthigh);
        graphCombEtaInvCrossSectionSystErrPPRef->SetPointEYlow(iPoint,errYSystlow);
        
        
  }
  
  
  
  
 TGraphAsymmErrors* graphCombEtaInvCrossSectionTotErrPPRef  =  CalculateCombinedSysAndStatError( graphCombEtaInvCrossSectionStatErrPPRef ,graphCombEtaInvCrossSectionSystErrPPRef );
  
  
  
  
  
  
  
  
  
  
  
  
  //Combination pp reference
  
  
  TString fileNameOutputWeightingPi0PPRef                = Form("%s/WeightingPi0PPRef.dat",outputDir.Data());
  
  
                                        
  statErrorCollectionPi0PPRef[0]          = (TH1D*)histoPCMppreferenceStatErr->Clone("statErrPCMPi0ppRef");
  statErrorCollectionPi0PPRef[1]          = (TH1D*)histoPHOSppreferenceStatErr->Clone("statErrPHOSPi0ppRef");
  statErrorCollectionPi0PPRef[2]          = (TH1D*)histoEMCalppreferenceStatErr->Clone("statErrEMCalPi0ppRef");
  statErrorCollectionPi0PPRef[4]	  = (TH1D*)histoPCMEMCalppreferenceStatErr->Clone("statErrPCM-EMCalPi0ppRef");
  statErrorCollectionPi0PPRef[5]          = (TH1D*)histoDalitzppreferenceStatErr->Clone("statErrDalitzPi0ppRef");
    
  sysErrorCollectionPi0PPRef[0]           = (TGraphAsymmErrors*)graphPCMppreferenceSystErr->Clone("sysErrPCMPi0ppRef");
  sysErrorCollectionPi0PPRef[1]           = (TGraphAsymmErrors*)graphPHOSppreferenceSystErr->Clone("sysErrPHOSPi0ppRef");
  sysErrorCollectionPi0PPRef[2]           = (TGraphAsymmErrors*)graphEMCalppreferenceSystErr->Clone("sysErrEMCalPi0ppRef");
  sysErrorCollectionPi0PPRef[4]		  = (TGraphAsymmErrors*)graphPCMEMCalppreferenceSystErr->Clone("sysErrPCM-EMCalPi0ppRef");
  sysErrorCollectionPi0PPRef[5]           = (TGraphAsymmErrors*)graphDalitzppreferenceSystErr->Clone("sysErrDalitzPi0ppRef");
  
    
  TGraphAsymmErrors* graphCombPi0InvCrossSectionStatpp5023GeV= NULL;
  TGraphAsymmErrors* graphCombPi0InvCrossSectionSyspp5023GeV = NULL;
  
  cout<<"\n\n */////////////////Combination of Pi0 pp Ref  begins here**//////////////"<<endl;
 
    
  TGraphAsymmErrors* graphCombPi0InvCrossSectionTotpp5023GeV = CombinePtPointsSpectraFullCorrMat(    statErrorCollectionPi0PPRef,    sysErrorCollectionPi0PPRef,     
												      pTLimits, Ntotal,
												      offSets, offSetsSys,
												      graphCombPi0InvCrossSectionStatpp5023GeV, graphCombPi0InvCrossSectionSyspp5023GeV,
												      fileNameOutputWeightingPi0PPRef,"pPb_5.023TeV","Pi0RpPb",kTRUE,NULL,corrFactorsFileName.Data()                                                                             ); 
  
  
  
  
  
  TGraphAsymmErrors* graphCombPi0InvCrossSectionTotErrPPReffitTsallis = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionTotErrPPRef->Clone();
  
  TGraphAsymmErrors* graphCombEtaInvCrossSectionTotErrPPReffitTsallis = (TGraphAsymmErrors*)graphCombEtaInvCrossSectionTotErrPPRef->Clone();
  
  
  TGraphAsymmErrors* graphCombPi0InvCrossSectionTotErrPPReffitTCM = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionTotErrPPRef->Clone();
  
  TGraphAsymmErrors* graphCombEtaInvCrossSectionTotErrPPReffitTCM = (TGraphAsymmErrors*)graphCombEtaInvCrossSectionTotErrPPRef->Clone();


  
  
  
  TF1* fitCombPi0ppRef5023GeVPt;
  TF1* fitTsallisPi0ppRef5023GeVPt;
  TF1* fitTsallisEtappRef5023GeVPt;
  TF1* fitTCMCombPi0ppRef5023GeVPt;
  TF1* fitTCMPi0ppRef5023GeVPt;
  TF1* fitTCMEtappRef5023GeVPt;
  
  
  
  
  Double_t *ParametersTsallisPi0ppRef  = new Double_t[5];
  Double_t *ParametersTCMPi0ppRef      = new Double_t[5];
  
  Double_t minPt = 0.3; 
  Double_t maxPt = 20.0;
  
  fitCombPi0ppRef5023GeVPt        =  FitObject("l","fitInvCrossSectionPi0","Pi0");
  fitTsallisPi0ppRef5023GeVPt     =  FitObject("l","fitInvCrossSectionPi0","Pi0");
  fitTsallisEtappRef5023GeVPt            =  FitObject("l","fitInvCrossSectionEta","Eta");
  ParametersTsallisPi0ppRef[0] = 8.66870e+08;//7.4e+1;
  ParametersTsallisPi0ppRef[1] = 7.03073e+00;
  ParametersTsallisPi0ppRef[2] = 1.63511e-01;
  fitCombPi0ppRef5023GeVPt->SetRange(minPt,maxPt);
  fitCombPi0ppRef5023GeVPt->SetParameters(ParametersTsallisPi0ppRef[0],ParametersTsallisPi0ppRef[1],ParametersTsallisPi0ppRef[2]);
  fitTsallisPi0ppRef5023GeVPt->SetRange(minPt,maxPt);
  fitTsallisPi0ppRef5023GeVPt->SetParameters(ParametersTsallisPi0ppRef[0],ParametersTsallisPi0ppRef[1],ParametersTsallisPi0ppRef[2]);
  fitTsallisEtappRef5023GeVPt->SetRange(0.7,maxPt);
  fitTsallisEtappRef5023GeVPt->SetParameters(ParametersTsallisPi0ppRef[0],ParametersTsallisPi0ppRef[1],ParametersTsallisPi0ppRef[2]);
    
  fitTCMCombPi0ppRef5023GeVPt = FitObject("tcm","fitInvCrossSectionPi0","Pi0");
  fitTCMPi0ppRef5023GeVPt     = FitObject("tcm","fitInvCrossSectionPi0","Pi0");
  fitTCMEtappRef5023GeVPt     = FitObject("tcm","fitInvCrossSectionEta","Eta");
  
    //Eta
    //Parameters8TeV[0] = 1.48e+09;
    //Parameters8TeV[1] = 2.25e-01;
    //Parameters8TeV[2] = 2.98e+09;
    //Parameters8TeV[3] = 8.050e-01;
    //Parameters8TeV[4] = 3.041e+00;
  
  ParametersTCMPi0ppRef[0] = 6.69e+11;
  ParametersTCMPi0ppRef[1] = 1.439e-01;
  ParametersTCMPi0ppRef[2] = 3.44e+10;
  ParametersTCMPi0ppRef[3] = 6.040e-01;
  ParametersTCMPi0ppRef[4] = 3.028e+00;
  
  fitTCMCombPi0ppRef5023GeVPt->SetRange(minPt,maxPt);
  fitTCMCombPi0ppRef5023GeVPt->SetParameters(ParametersTCMPi0ppRef[0],ParametersTCMPi0ppRef[1],ParametersTCMPi0ppRef[2],ParametersTCMPi0ppRef[3],ParametersTCMPi0ppRef[4]); // standard
  fitTCMPi0ppRef5023GeVPt->SetRange(minPt,maxPt);
  fitTCMPi0ppRef5023GeVPt->SetParameters(ParametersTCMPi0ppRef[0],ParametersTCMPi0ppRef[1],ParametersTCMPi0ppRef[2],ParametersTCMPi0ppRef[3],ParametersTCMPi0ppRef[4]);

  fitTCMEtappRef5023GeVPt->SetRange(0.7,maxPt);
  fitTCMEtappRef5023GeVPt->SetParameters(ParametersTCMPi0ppRef[0],ParametersTCMPi0ppRef[1],ParametersTCMPi0ppRef[2],ParametersTCMPi0ppRef[3],ParametersTCMPi0ppRef[4]);

  
  graphCombPi0InvCrossSectionTotpp5023GeV->Fit(fitCombPi0ppRef5023GeVPt,"QVNRMEX0+","",minPt,maxPt);
  cout<<"!!!!!!!!Fit!!!!!!!!!!!!"<< endl;
  graphCombPi0InvCrossSectionTotpp5023GeV->Fit(fitCombPi0ppRef5023GeVPt,"VNRMEX0+","",minPt,maxPt);
  cout<<"!!!!!!!!Fit!!!!!!!!!!!!"<< endl;
  
  TGraphAsymmErrors* graphRatioCombPi0ppRefFitTotErr  = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionTotpp5023GeV->Clone();
  TGraphAsymmErrors* graphRatioCombPi0ppRefFitStatErr = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionStatpp5023GeV->Clone();
  TGraphAsymmErrors* graphRatioCombPi0ppRefFitSystErr = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionSyspp5023GeV->Clone();
  
  graphRatioCombPi0ppRefFitTotErr     = CalculateGraphErrRatioToFit(graphRatioCombPi0ppRefFitTotErr,fitCombPi0ppRef5023GeVPt);
  graphRatioCombPi0ppRefFitStatErr    = CalculateGraphErrRatioToFit(graphRatioCombPi0ppRefFitStatErr,fitCombPi0ppRef5023GeVPt);
  graphRatioCombPi0ppRefFitSystErr    = CalculateGraphErrRatioToFit(graphRatioCombPi0ppRefFitSystErr,fitCombPi0ppRef5023GeVPt);
  
  

  
  TGraphAsymmErrors*   graphRatioPi0PPRefCombPCMSys      = (TGraphAsymmErrors*)graphPCMppreferenceSystErr->Clone("graphRatioPi0PPRefCombPCMSys");
  TGraphAsymmErrors*   graphRatioPi0PPRefCombDalitzSys   = (TGraphAsymmErrors*)graphDalitzppreferenceSystErr->Clone("graphRatioPi0PPRefCombDalitzSys");
  TGraphAsymmErrors*   graphRatioPi0PPRefCombEMCalSys    = (TGraphAsymmErrors*)graphEMCalppreferenceSystErr->Clone("graphRatioPi0PPRefCombEMCalSys");
  TGraphAsymmErrors*   graphRatioPi0PPRefCombPHOSSys     = (TGraphAsymmErrors*)graphPHOSppreferenceSystErr->Clone("graphRatioPi0PPRefCombPHOSSys");
  TGraphAsymmErrors*   graphRatioPi0PPRefCombPCMEMCalSys = (TGraphAsymmErrors*)graphPCMEMCalppreferenceSystErr->Clone("graphRatioPi0PPRefCombPCMEMCalSys");
  
  
  cout<<"PHOS sys"<<endl;
  graphPHOSppreferenceSystErr->Print();
  cout<<"PHOS stat"<<endl;
  graphPHOSppreferenceStatErr->Print();
  
  
  
  TH1D*   histoRatioPi0PPRefCombPCMStat      = (TH1D*)histoPCMppreferenceStatErr->Clone("histoRatioPi0PPRefCombPCMStat");
  TH1D*   histoRatioPi0PPRefCombDalitzStat   = (TH1D*)histoDalitzppreferenceStatErr->Clone("histoRatioPi0PPRefCombDalitzStat");
  TH1D*   histoRatioPi0PPRefCombEMCalStat    = (TH1D*)histoEMCalppreferenceStatErr->Clone("histoRatioPi0PPRefCombEMCalStat");
  TH1D*   histoRatioPi0PPRefCombPHOSStat     = (TH1D*)histoPHOSppreferenceStatErr->Clone("histoRatioPi0PPRefCombPHOSStat");
  TH1D*   histoRatioPi0PPRefCombPCMEMCalStat = (TH1D*)histoPCMEMCalppreferenceStatErr->Clone("histoRatioPi0PPRefCombPCMEMCalStat");
  

  graphRatioPi0PPRefCombPCMSys        = CalculateGraphErrRatioToFit(graphRatioPi0PPRefCombPCMSys,      fitCombPi0ppRef5023GeVPt);
  graphRatioPi0PPRefCombDalitzSys     = CalculateGraphErrRatioToFit(graphRatioPi0PPRefCombDalitzSys,   fitCombPi0ppRef5023GeVPt);
  graphRatioPi0PPRefCombEMCalSys      = CalculateGraphErrRatioToFit(graphRatioPi0PPRefCombEMCalSys,    fitCombPi0ppRef5023GeVPt);
  graphRatioPi0PPRefCombPHOSSys       = CalculateGraphErrRatioToFit(graphRatioPi0PPRefCombPHOSSys,     fitCombPi0ppRef5023GeVPt);
  graphRatioPi0PPRefCombPCMEMCalSys   = CalculateGraphErrRatioToFit(graphRatioPi0PPRefCombPCMEMCalSys, fitCombPi0ppRef5023GeVPt);
  
  
  histoRatioPi0PPRefCombPCMStat       = CalculateHistoRatioToFit(histoRatioPi0PPRefCombPCMStat,	        fitCombPi0ppRef5023GeVPt);
  histoRatioPi0PPRefCombDalitzStat    = CalculateHistoRatioToFit(histoRatioPi0PPRefCombDalitzStat,	fitCombPi0ppRef5023GeVPt);
  histoRatioPi0PPRefCombEMCalStat     = CalculateHistoRatioToFit(histoRatioPi0PPRefCombEMCalStat,	        fitCombPi0ppRef5023GeVPt);
  histoRatioPi0PPRefCombPHOSStat      = CalculateHistoRatioToFit(histoRatioPi0PPRefCombPHOSStat,	        fitCombPi0ppRef5023GeVPt);
  histoRatioPi0PPRefCombPCMEMCalStat  = CalculateHistoRatioToFit(histoRatioPi0PPRefCombPCMEMCalStat, 	fitCombPi0ppRef5023GeVPt);

  
  



  graphCombPi0InvCrossSectionTotErrPPReffitTsallis->Fit(fitTsallisPi0ppRef5023GeVPt,"QVNRMEX0+","",minPt,maxPt);
  cout<<"!!!!!!!!Fit!!!!!!!!!!!!"<< endl;
  graphCombPi0InvCrossSectionTotErrPPReffitTsallis->Fit(fitTsallisPi0ppRef5023GeVPt,"VNRMEX0+","",minPt,maxPt);
  cout<<"!!!!!!!!Fit!!!!!!!!!!!!"<< endl;
  
  
  fFits<<"Fitting Pi0 PP reference spectrum with Tsallis"<< endl;
 
    for (Int_t i=0;i < 3;i++){
      fFits<< fitTsallisPi0ppRef5023GeVPt->GetParName(i) <<":  "<< fitTsallisPi0ppRef5023GeVPt->GetParameter(i)<<" +/- "<< fitTsallisPi0ppRef5023GeVPt->GetParError(i) << endl;
    }
  fFits<< "Chi2/NDF:  "<<fitTsallisPi0ppRef5023GeVPt->GetChisquare()/fitTsallisPi0ppRef5023GeVPt->GetNDF() << endl;
 
  
  
  
  
  TGraphAsymmErrors* graphRatioPi0ppRefFitTotErr  = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionTotErrPPRef->Clone();
  TGraphAsymmErrors* graphRatioPi0ppRefFitStatErr = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionStatErrPPRef->Clone();
  TGraphAsymmErrors* graphRatioPi0ppRefFitSystErr = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionSystErrPPRef->Clone();
  
  graphRatioPi0ppRefFitTotErr     = CalculateGraphErrRatioToFit(graphRatioPi0ppRefFitTotErr,fitTsallisPi0ppRef5023GeVPt);
  graphRatioPi0ppRefFitStatErr    = CalculateGraphErrRatioToFit(graphRatioPi0ppRefFitStatErr,fitTsallisPi0ppRef5023GeVPt);
  graphRatioPi0ppRefFitSystErr    = CalculateGraphErrRatioToFit(graphRatioPi0ppRefFitSystErr,fitTsallisPi0ppRef5023GeVPt);
  
  
  
  graphCombPi0InvCrossSectionTotErrPPReffitTCM->Fit(fitTCMPi0ppRef5023GeVPt,"QVNRMEX0+","",minPt,maxPt);
  cout<<"!!!!!!!!Fit!!!!!!!!!!!!"<< endl;
  graphCombPi0InvCrossSectionTotErrPPReffitTCM->Fit(fitTCMPi0ppRef5023GeVPt,"VNRMEX0+","",minPt,maxPt);
  cout<<"!!!!!!!!Fit!!!!!!!!!!!!"<< endl;
  
  
   fFits<<"Fitting Pi0 PP reference spectrum with TCM"<< endl;
 
    for (Int_t i=0;i < 5;i++){
      fFits<< fitTCMPi0ppRef5023GeVPt->GetParName(i) <<":  "<< fitTCMPi0ppRef5023GeVPt->GetParameter(i)<<" +/- "<< fitTCMPi0ppRef5023GeVPt->GetParError(i) << endl;
    }
  fFits<< "Chi2/NDF:  "<<fitTCMPi0ppRef5023GeVPt->GetChisquare()/fitTCMPi0ppRef5023GeVPt->GetNDF() << endl;
  
  
  TGraphAsymmErrors* graphRatioPi0ppRefFitTCMTotErr  = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionTotErrPPRef->Clone();
  TGraphAsymmErrors* graphRatioPi0ppRefFitTCMStatErr = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionStatErrPPRef->Clone();
  TGraphAsymmErrors* graphRatioPi0ppRefFitTCMSystErr = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionSystErrPPRef->Clone();
  
  graphRatioPi0ppRefFitTCMTotErr     = CalculateGraphErrRatioToFit(graphRatioPi0ppRefFitTCMTotErr,fitTCMPi0ppRef5023GeVPt);
  graphRatioPi0ppRefFitTCMStatErr    = CalculateGraphErrRatioToFit(graphRatioPi0ppRefFitTCMStatErr,fitTCMPi0ppRef5023GeVPt);
  graphRatioPi0ppRefFitTCMSystErr    = CalculateGraphErrRatioToFit(graphRatioPi0ppRefFitTCMSystErr,fitTCMPi0ppRef5023GeVPt);
  
  

  
  graphCombEtaInvCrossSectionTotErrPPReffitTsallis->Fit(fitTsallisEtappRef5023GeVPt,"QVNRMEX0+","",0.7,maxPt);
  cout<<"!!!!!!!!Fit!!!!!!!!!!!!"<< endl;
  graphCombEtaInvCrossSectionTotErrPPReffitTsallis->Fit(fitTsallisEtappRef5023GeVPt,"VNRMEX0+","",0.7,maxPt);
  cout<<"!!!!!!!!Fit!!!!!!!!!!!!"<< endl;
  
   fFits<<"Fitting Eta PP reference spectrum with Tsallis"<< endl;
 
    for (Int_t i=0;i < 3;i++){
      fFits<< fitTsallisEtappRef5023GeVPt->GetParName(i) <<":  "<< fitTsallisEtappRef5023GeVPt->GetParameter(i)<<" +/- "<< fitTsallisEtappRef5023GeVPt->GetParError(i) << endl;
    }
  fFits<< "Chi2/NDF:  "<<fitTsallisEtappRef5023GeVPt->GetChisquare()/fitTsallisEtappRef5023GeVPt->GetNDF() << endl;
 
  
  
  
  
  
  TGraphAsymmErrors* graphRatioEtappRefFitTotErr  = (TGraphAsymmErrors*)graphCombEtaInvCrossSectionTotErrPPRef->Clone();
  TGraphAsymmErrors* graphRatioEtappRefFitStatErr = (TGraphAsymmErrors*)graphCombEtaInvCrossSectionStatErrPPRef->Clone();
  TGraphAsymmErrors* graphRatioEtappRefFitSystErr = (TGraphAsymmErrors*)graphCombEtaInvCrossSectionSystErrPPRef->Clone();
  
  graphRatioEtappRefFitTotErr     = CalculateGraphErrRatioToFit(graphRatioEtappRefFitTotErr,fitTsallisEtappRef5023GeVPt);
  graphRatioEtappRefFitStatErr    = CalculateGraphErrRatioToFit(graphRatioEtappRefFitStatErr,fitTsallisEtappRef5023GeVPt);
  graphRatioEtappRefFitSystErr    = CalculateGraphErrRatioToFit(graphRatioEtappRefFitSystErr,fitTsallisEtappRef5023GeVPt);
  
  
  
  
  
  graphCombEtaInvCrossSectionTotErrPPReffitTCM->Fit(fitTCMEtappRef5023GeVPt,"QVNRMEX0+","",0.7,maxPt);
  cout<<"!!!!!!!!Fit!!!!!!!!!!!!"<< endl;
  graphCombEtaInvCrossSectionTotErrPPReffitTCM->Fit(fitTCMEtappRef5023GeVPt,"VNRMEX0+","",0.7,maxPt);
  cout<<"!!!!!!!!Fit!!!!!!!!!!!!"<< endl; fFits<< FittingType.Data()<<endl;
  
  fFits<<"Fitting Eta PP reference spectrum with TCM"<< endl;
 
    for (Int_t i=0;i < 5;i++){
      fFits<< fitTCMEtappRef5023GeVPt->GetParName(i) <<":  "<< fitTCMEtappRef5023GeVPt->GetParameter(i)<<" +/- "<< fitTCMEtappRef5023GeVPt->GetParError(i) << endl;
    }
  fFits<< "Chi2/NDF:  "<<fitTCMEtappRef5023GeVPt->GetChisquare()/fitTCMEtappRef5023GeVPt->GetNDF() << endl;

  
  
  
  
  
  TGraphAsymmErrors* graphRatioEtappRefFitTCMTotErr  = (TGraphAsymmErrors*)graphCombEtaInvCrossSectionTotErrPPRef->Clone();
  TGraphAsymmErrors* graphRatioEtappRefFitTCMStatErr = (TGraphAsymmErrors*)graphCombEtaInvCrossSectionStatErrPPRef->Clone();
  TGraphAsymmErrors* graphRatioEtappRefFitTCMSystErr = (TGraphAsymmErrors*)graphCombEtaInvCrossSectionSystErrPPRef->Clone();
  
  graphRatioEtappRefFitTCMTotErr     = CalculateGraphErrRatioToFit(graphRatioEtappRefFitTCMTotErr,fitTCMEtappRef5023GeVPt);
  graphRatioEtappRefFitTCMStatErr    = CalculateGraphErrRatioToFit(graphRatioEtappRefFitTCMStatErr,fitTCMEtappRef5023GeVPt);
  graphRatioEtappRefFitTCMSystErr    = CalculateGraphErrRatioToFit(graphRatioEtappRefFitTCMSystErr,fitTCMEtappRef5023GeVPt);
  
  
  
  
  
  TCanvas* canvas0a = new TCanvas("canvas0a","",200,10,1000,1200);  // gives the page size
  DrawGammaCanvasSettings( canvas0a, 0.3, 0.02, 0.02, 0.16);
  TPad* pad0a = new TPad("pad0a", "", 0., 0.25, 1., 1.,-1, -1, -2);
  DrawGammaPadSettings( pad0a, 0.18, 0.02, 0.02, 0.);
  pad0a->Draw();

  TPad* padRatios0a = new TPad("padRatios0a", "", 0., 0., 1., 0.25,-1, -1, -2);
  DrawGammaPadSettings( padRatios0a, 0.18, 0.02, 0., 0.3);
  padRatios0a->Draw();

  pad0a->cd(); 
  pad0a->SetLogx();
  pad0a->SetLogy();
  
  TH2F* histo0a = new TH2F("histo0a","histo0a",1000,0.2,30.,1000,1e02,2e11 );
  SetStyleHistoTH2ForGraphs(histo0a, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2} ",0.040,0.040, 0.040,0.040, 0.8,1.9, 512, 510);
  histo0a->GetXaxis()->SetLabelOffset(-0.009);
  histo0a->DrawCopy(); 

  
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionSystErrPPRef,20,1,kBlue+2 ,kBlue+2, 1, kTRUE);
  graphCombPi0InvCrossSectionSystErrPPRef->Draw("E2,same"); 
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionStatErrPPRef,20,1,kBlue+2 ,kBlue+2  ); 
  graphCombPi0InvCrossSectionStatErrPPRef->Draw("pz,same"); 
  
  
  
  fitTsallisPi0ppRef5023GeVPt->SetLineWidth(2); 

  
  fitTsallisPi0ppRef5023GeVPt->Draw("same");
  
  
  
  TLegend* leg0a;
  leg0a   = new TLegend(0.25,0.31,0.45,0.41); 
  leg0a->SetFillColor(0);
  leg0a->SetLineColor(0);
  leg0a->SetTextFont(42);
  leg0a->SetTextSize(0.045-0.005);
  leg0a->AddEntry(graphCombPi0InvCrossSectionSystErrPPRef,"#pi^{0} pp reference","pef");
  leg0a->AddEntry(fitTsallisPi0ppRef5023GeVPt,"Tsallis fit","l");
  leg0a->Draw();
  
  
  
 	
  TLatex * lt0a = new TLatex(4.5,0.3,"ALICE") ;
  lt0a->SetTextColor(kBlack) ;
  lt0a->SetTextSize(0.045-0.005) ;
  lt0a->SetTextFont(42) ;
  lt0a->DrawLatex(4.5,0.08,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  lt0a->Draw() ;


  padRatios0a->cd();
  padRatios0a->SetLogx();

  
 
  
  
  TH2F * histoRatio0a;
  histoRatio0a = new TH2F("histoRatio0a","histo2DRatioPi0PPRefToFit",1000,.2,30.,1000,0.5,1.5);//1.89
  
  SetStyleHistoTH2ForGraphs(histoRatio0a, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.125,0.125, 0.125,0.125,1.,0.5, 502, 505); 
  histoRatio0a->GetYaxis()->SetLabelOffset(0.005);
  histoRatio0a->GetXaxis()->SetLabelOffset(0.015+0.00);
  histoRatio0a->GetXaxis()->SetTickLength(0.07);
  histoRatio0a->DrawCopy();
 
  
  
 	
  

  TLine *lineRatio0a=new TLine(0.,1.,30.,1.);
  lineRatio0a->SetLineColor(kGray+1);
  lineRatio0a->Draw("same");

  
  
  
  
  DrawGammaSetMarkerTGraphAsym(graphRatioPi0ppRefFitSystErr,20,1,kBlue+2 ,kBlue+2, 1, kTRUE);
  graphRatioPi0ppRefFitSystErr->Draw("E2,same"); 
  DrawGammaSetMarkerTGraphAsym(graphRatioPi0ppRefFitStatErr,20,1,kBlue+2 ,kBlue+2  ); 
  graphRatioPi0ppRefFitStatErr->Draw("pz,same"); 
  
  
	
   
  
  canvas0a->Update();
  canvas0a->Print(Form("%s/Pi0_pp_Ref_Ratio_To_Fit_Tsallis.%s",outputDir.Data(),suffix.Data()));
  
  
  
  
  
  TCanvas* canvas0aTCM = new TCanvas("canvas0aTCM","",200,10,1000,1200);  // gives the page size
  DrawGammaCanvasSettings( canvas0aTCM, 0.3, 0.02, 0.02, 0.16);
  TPad* pad0aTCM = new TPad("pad0aTCM", "", 0., 0.25, 1., 1.,-1, -1, -2);
  DrawGammaPadSettings( pad0aTCM, 0.18, 0.02, 0.02, 0.);
  pad0aTCM->Draw();

  TPad* padRatios0aTCM = new TPad("padRatios0aTCM", "", 0., 0., 1., 0.25,-1, -1, -2);
  DrawGammaPadSettings( padRatios0aTCM, 0.18, 0.02, 0., 0.3);
  padRatios0aTCM->Draw();

  pad0aTCM->cd(); 
  pad0aTCM->SetLogx();
  pad0aTCM->SetLogy();
  
  TH2F* histo0aTCM = new TH2F("histo0aTCM","histo0aTCM",1000,0.2,30.,1000,1e02,2e11 );
  SetStyleHistoTH2ForGraphs(histo0aTCM, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2} ",0.040,0.040, 0.040,0.040, 0.8,1.9, 512, 510);
  histo0aTCM->GetXaxis()->SetLabelOffset(-0.009);
  histo0aTCM->DrawCopy(); 

  
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionSystErrPPRef,20,1,kBlue+2 ,kBlue+2, 1, kTRUE);
  graphCombPi0InvCrossSectionSystErrPPRef->Draw("E2,same"); 
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionStatErrPPRef,20,1,kBlue+2 ,kBlue+2  ); 
  graphCombPi0InvCrossSectionStatErrPPRef->Draw("pz,same"); 
  
  
  
  fitTCMPi0ppRef5023GeVPt->SetLineWidth(2); 

  
  fitTCMPi0ppRef5023GeVPt->Draw("same");
  
  
  
  TLegend* leg0aTCM;
  leg0aTCM   = new TLegend(0.25,0.31,0.45,0.41); 
  leg0aTCM->SetFillColor(0);
  leg0aTCM->SetLineColor(0);
  leg0aTCM->SetTextFont(42);
  leg0aTCM->SetTextSize(0.045-0.005);
  leg0aTCM->AddEntry(graphCombPi0InvCrossSectionSystErrPPRef,"#pi^{0} pp reference","pef");
  leg0aTCM->AddEntry(fitTCMPi0ppRef5023GeVPt,"TCM fit","l");
  leg0aTCM->Draw();
  
  
  
 	
  TLatex * lt0aTCM = new TLatex(4.5,0.3,"ALICE") ;
  lt0aTCM->SetTextColor(kBlack) ;
  lt0aTCM->SetTextSize(0.045-0.005) ;
  lt0aTCM->SetTextFont(42) ;
  lt0aTCM->DrawLatex(4.5,0.08,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  lt0aTCM->Draw() ;


  padRatios0aTCM->cd();
  padRatios0aTCM->SetLogx();

  
 
  
  
  TH2F * histoRatio0aTCM;
  histoRatio0aTCM = new TH2F("histoRatio0aTCM","histo2DRatioPi0PPRefToFit",1000,.2,30.,1000,0.5,1.5);//1.89
  SetStyleHistoTH2ForGraphs(histoRatio0aTCM, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.125,0.125, 0.125,0.125,1.,0.5, 502, 505); 
  histoRatio0aTCM->GetYaxis()->SetLabelOffset(0.005);
  histoRatio0aTCM->GetXaxis()->SetLabelOffset(0.015+0.00);
  histoRatio0aTCM->GetXaxis()->SetTickLength(0.07);
  histoRatio0aTCM->DrawCopy();
 
  
  
 	
  

  TLine *lineRatio0aTCM=new TLine(0.,1.,30.,1.);
  lineRatio0aTCM->SetLineColor(kGray+1);
  lineRatio0aTCM->Draw("same");

  
  
  
  
  DrawGammaSetMarkerTGraphAsym(graphRatioPi0ppRefFitTCMSystErr,20,1,kBlue+2 ,kBlue+2, 1, kTRUE);
  graphRatioPi0ppRefFitTCMSystErr->Draw("E2,same"); 
  DrawGammaSetMarkerTGraphAsym(graphRatioPi0ppRefFitTCMStatErr,20,1,kBlue+2 ,kBlue+2  ); 
  graphRatioPi0ppRefFitTCMStatErr->Draw("pz,same"); 
  
  
	
   
  
  canvas0aTCM->Update();
  canvas0aTCM->Print(Form("%s/Pi0_pp_Ref_Ratio_To_Fit_TCM.%s",outputDir.Data(),suffix.Data()));
  
  
  
  
  
  Int_t CDimX=500;
  Int_t CDimY=400;
  
  Double_t CMarginL=0.11;
  Double_t CMarginT=0.02;
  Double_t CMarginR=0.02;
  Double_t CMarginB=0.12;
  
  
  
  
   TCanvas* c0b = new TCanvas("c2","",200,10,CDimX,CDimY);  // gives the page size
   DrawGammaCanvasSettings( c0b, CMarginL, CMarginR,CMarginT ,CMarginB);
   c0b->SetLogx();
   TH2F * hist0b;
   hist0b = new TH2F("hist0b","hist0b",1000,0.25,25.,1000,0.5,2.0);
   SetStyleHistoTH2ForGraphs(hist0b, "#it{p}_{T} (GeV/#it{c})","Data/ Fit",LabelSizeX,TitleSizeX,LabelSizeY,TitleSizeY,TitleOffsetX,TitleOffsetY, 512, 508); 
   hist0b->GetXaxis()->SetLabelOffset(LabelOffsetLog);
   hist0b->DrawCopy(); 

  
   DrawGammaLines(0.0,25.,1.,1.,2.0,kGray+2,2);
  
  

   DrawGammaSetMarkerTGraphAsym(graphRatioPi0PPRefCombPHOSSys,markerStylePHOS,markerSizePHOS, colorPHOS, colorPHOS, 1, kTRUE);  
   graphRatioPi0PPRefCombPHOSSys->Draw("E2same");
   DrawGammaSetMarker(histoRatioPi0PPRefCombPHOSStat, markerStylePHOS,markerSizePHOS, colorPHOS, colorPHOS );
   histoRatioPi0PPRefCombPHOSStat->Draw("EX0,same") ;
 
   DrawGammaSetMarkerTGraphAsym(graphRatioPi0PPRefCombEMCalSys,markerStyleEMCal,markerSizeEMCal, colorEMCal, colorEMCal, 1, kTRUE);  
   graphRatioPi0PPRefCombEMCalSys->Draw("E2same");
   DrawGammaSetMarker(histoRatioPi0PPRefCombEMCalStat,markerStyleEMCal,markerSizeEMCal, colorEMCal, colorEMCal );
   histoRatioPi0PPRefCombEMCalStat->Draw("EX0,same") ; 

   DrawGammaSetMarkerTGraphAsym(graphRatioPi0PPRefCombPCMSys,markerStylePCM,markerSizePCM, colorPCM, colorPCM, 1, kTRUE);  
   graphRatioPi0PPRefCombPCMSys->Draw("E2same");
   DrawGammaSetMarker(histoRatioPi0PPRefCombPCMStat,markerStylePCM,markerSizePCM, colorPCM, colorPCM );
   histoRatioPi0PPRefCombPCMStat->Draw("EX0,same") ;
 
   DrawGammaSetMarkerTGraphAsym(graphRatioPi0PPRefCombDalitzSys,markerStyleDalitz,markerSizeDalitz, colorDalitz, colorDalitz, 1, kTRUE);  
   graphRatioPi0PPRefCombDalitzSys->Draw("E2same");
   DrawGammaSetMarker(histoRatioPi0PPRefCombDalitzStat,markerStyleDalitz,markerSizeDalitz, colorDalitz, colorDalitz );
   histoRatioPi0PPRefCombDalitzStat->Draw("EX0,same") ; 
  
   DrawGammaSetMarkerTGraphAsym(graphRatioPi0PPRefCombPCMEMCalSys,markerStylePCMEMCal,markerSizePCMEMCal,colorPCMEMCal,colorPCMEMCal,1,kTRUE);
  graphRatioPi0PPRefCombPCMEMCalSys->Draw("E2same");
   DrawGammaSetMarker(histoRatioPi0PPRefCombPCMEMCalStat,markerStylePCMEMCal,markerSizePCMEMCal, colorPCMEMCal, colorPCMEMCal );
  histoRatioPi0PPRefCombPCMEMCalStat->Draw("EX0,same") ; 
  
  
  
 	
  TLatex * lt0b = new TLatex(1.0,2.61,"ALICE") ;
  lt0b->SetTextColor(kBlack) ;
  lt0b->SetTextSize(TextSize) ;
  lt0b->SetTextFont(Font) ;
  lt0b->DrawLatex(1.0,2.45,"p-Pb, NSD, #sqrt{#it{s}_{NN}} = 5.02 TeV, #pi^{0}");
  lt0b->Draw() ;
  
  TLegend* leg0b = new TLegend(0.13,0.65,0.45,0.95);
  leg0b->SetFillColor(0);
  leg0b->SetLineColor(0);
  leg0b->SetTextFont(Font);
  leg0b->SetTextSize(TextSize);
  leg0b->AddEntry(histoRatioPi0PPRefCombPHOSStat,Form("PHOS"),"pf");
  leg0b->AddEntry(histoRatioPi0PPRefCombEMCalStat,Form("EMC"),"pf");
  leg0b->AddEntry(histoRatioPi0PPRefCombPCMStat,Form("PCM"),"pf");
  leg0b->AddEntry(graphRatioPi0PPRefCombPCMEMCalSys,Form("PCM-EMC"),"pf");
  leg0b->AddEntry(histoRatioPi0PPRefCombDalitzStat,Form("PCM-#gamma*#gamma"),"pf");
  leg0b->Draw("same");





  c0b->Update();
  c0b->Print(Form("%s/RatioCombFitIndividualPi0PPRef.%s",outputDir.Data(),suffix.Data()));


  cout<<"\n\n */////////////////Combination of Pi0 pp Ref  ends here**//////////////"<<endl;
 


  cout<<"\n\n */////////////////Combination of Eta pp Ref  begins here**//////////////"<<endl;


  TString fileNameOutputWeightingEtaPPRef                = Form("%s/WeightingEtaPPRef.dat",outputDir.Data());
  
 
  statErrorCollectionEtaPPRef[0]          = (TH1D*)histoPCMEtappreferenceStatErr->Clone("statErrPCMEtappRef");
  statErrorCollectionEtaPPRef[2]          = (TH1D*)histoEMCalEtappreferenceStatErr->Clone("statErrEMCalEtappRef");
  statErrorCollectionEtaPPRef[4]	  = (TH1D*)histoPCMEMCalEtappreferenceStatErr->Clone("statErrPCM-EMCalEtappRef");
    
  sysErrorCollectionEtaPPRef[0]           = (TGraphAsymmErrors*)graphPCMEtappreferenceSystErr->Clone("sysErrPCMEtappRef");
  sysErrorCollectionEtaPPRef[2]           = (TGraphAsymmErrors*)graphEMCalEtappreferenceSystErr->Clone("sysErrEMCalEtappRef");
  sysErrorCollectionEtaPPRef[4]           = (TGraphAsymmErrors*)graphPCMEMCalEtappreferenceSystErr->Clone("sysErrPCM-EMCalEtappRef");
  
    
  TGraphAsymmErrors* graphCombEtaInvCrossSectionStatpp5023GeV = NULL;
  TGraphAsymmErrors* graphCombEtaInvCrossSectionSyspp5023GeV  = NULL;

 
  TGraphAsymmErrors* graphCombEtaInvCrossSectionTotpp5023GeV = CombinePtPointsSpectraFullCorrMat(    statErrorCollectionEtaPPRef,                                                                                                                                                                                                                                                                                                        sysErrorCollectionEtaRpPb,   
                                                                                                        pTLimitsEta, NtotalEta,
                                                                                                        offSetsEta, offSetsSysEta,
                                                                                                        graphCombEtaInvCrossSectionStatpp5023GeV, graphCombEtaInvCrossSectionSyspp5023GeV,
                                                                                                        fileNameOutputWeightingEtaPPRef,"pPb_5.023TeV","EtaRpPb",kTRUE,NULL,corrFactorsFileName.Data()                                                                             );



cout<<"\n\n */////////////////Combination of Eta pp Ref  ends here**//////////////"<<endl;


  /////////////////////////////////Comparison with the combined pp reference////////////////////////


  
  TGraphErrors* bla1a = NULL;
  TGraphErrors* bla2a = NULL;
  TGraphErrors* bla3a = NULL;
  TGraphErrors* bla4a = NULL;
  TGraphErrors* graphRatiotoCombinedPPRefBLUEMethod = CalculateRatioBetweenSpectraWithDifferentBinning(graphCombPi0InvCrossSectionStatpp5023GeV,graphCombPi0InvCrossSectionSyspp5023GeV ,graphCombPi0InvCrossSectionStatErrPPRef,graphCombPi0InvCrossSectionSystErrPPRef,  kTRUE,  kTRUE,&bla1a,&bla2a,&bla3a,&bla4a);
  
  
  TGraphErrors* bla1b = NULL;
  TGraphErrors* bla2b = NULL;
  TGraphErrors* bla3b = NULL;
  TGraphErrors* bla4b = NULL;
  TGraphErrors* graphRatiotoCombinedPPRef = CalculateRatioBetweenSpectraWithDifferentBinning(graphCombppreferenceStatErr,graphCombppreferenceSystErr ,graphCombPi0InvCrossSectionStatErrPPRef,graphCombPi0InvCrossSectionSystErrPPRef,  kTRUE,  kTRUE,&bla1b,&bla2b,&bla3b,&bla4b);
  
  
  

  
 
  
  
  
  TGraphErrors* bla1c = NULL;
  TGraphErrors* bla2c = NULL;
  TGraphErrors* bla3c = NULL;
  TGraphErrors* bla4c = NULL;
  
  TGraphErrors* graphRatioEtaPPRef = CalculateRatioBetweenSpectraWithDifferentBinning(graphCombEtappreferenceStatErr,graphCombEtappreferenceSystErr ,graphCombEtaInvCrossSectionStatpp5023GeV,graphCombEtaInvCrossSectionSyspp5023GeV,  kTRUE,  kTRUE,&bla1c,&bla2c,&bla3c,&bla4c)  ;
  graphRatioEtaPPRef->Print();
  
  

  TCanvas* c1a = new TCanvas("c1a","",200,10,1200,1200);  // gives the page size
  DrawGammaCanvasSettings(c1a , 0.3, 0.02, 0.02, 0.16);
  TPad* pad1a = new TPad("pad1a", "", 0., 0.25, 1., 1.,-1, -1, -2);
  DrawGammaPadSettings( pad1a, 0.15, 0.02, 0.02, 0.);
  pad1a->Draw();

  TPad* padRatios1a = new TPad("padRatios1a", "", 0., 0., 1., 0.25,-1, -1, -2);
  DrawGammaPadSettings( padRatios1a, 0.15, 0.02, 0., 0.2);
  padRatios1a->Draw();

  pad1a->cd();
  pad1a->SetLogx();
  pad1a->SetLogy();
  TH2F * histo1a;
  histo1a = new TH2F("histo1a","histo1a",1000,0.2,30.,1000,2e2,10e11);
  SetStyleHistoTH2ForGraphs(histo1a,"#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.04,0.05, 0.045,0.045, 0.7,1.4, 512, 505);
  histo1a->DrawCopy();


  
  
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionSystErrPPRef,20,0.7, kBlack, kBlack, 1, kTRUE);
  graphCombPi0InvCrossSectionSystErrPPRef->Draw("E2,same");
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionStatErrPPRef,20,0.7, kBlack, kBlack);
  graphCombPi0InvCrossSectionStatErrPPRef->Draw("p,same");
  

  DrawGammaSetMarkerTGraphAsym(graphCombppreferenceSystErr,20,0.7, 4, 4, 1, kTRUE);
  graphCombppreferenceSystErr->Draw("E2,same");
  DrawGammaSetMarkerTGraphAsym(graphCombppreferenceStatErr,20,0.7, kBlue, kBlue);
  graphCombppreferenceStatErr->Draw("p,same");

  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionSyspp5023GeV,20,0.7,kRed+1 ,kRed+1 , 1, kTRUE);
  graphCombPi0InvCrossSectionSyspp5023GeV->Draw("E2,same");
  DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionStatpp5023GeV,20,0.7, kRed+1, kRed+1);
  graphCombPi0InvCrossSectionStatpp5023GeV->Draw("p,same");

  TLegend* legend1a = new TLegend(0.23,0.10,0.7,0.30);
  legend1a->SetFillColor(0);
  legend1a->SetLineColor(0);
  legend1a->SetTextSize(0.03);
   legend1a->AddEntry(graphCombPi0InvCrossSectionSystErrPPRef,"Annika","pef");
  legend1a->AddEntry(graphCombppreferenceSystErr,"pp reference from combined","pef");
  legend1a->AddEntry(graphCombPi0InvCrossSectionSyspp5023GeV,"pp reference with BLUE Method","pef");
  TLatex * lt1a = new TLatex(3,1,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV") ;
  lt1a->SetTextColor(kBlack) ;
  lt1a->SetTextSize(0.04) ;
  lt1a->Draw("same");
  legend1a->Draw();

  padRatios1a->cd();

  padRatios1a->SetLogx();

  DrawGammaSetMarkerTGraphErr(graphRatiotoCombinedPPRef,20,1.0, kBlue, kBlue, 1, kTRUE);
  DrawGammaSetMarkerTGraphErr(graphRatiotoCombinedPPRefBLUEMethod,20,1.0, kRed+1, kRed+1, 1, kTRUE);
  
  TH2F * histoRatio1a;
  histoRatio1a = new TH2F("histoRatio1a","histo2Dppreference",1000,.2,30.,1000,0.81,1.19);
  SetStyleHistoTH2ForGraphs(histoRatio1a, "#it{p}_{T} (GeV/#it{c})","{Combined,BLUE}/Annika", 0.13,0.13, 0.13,0.13, 0.5,.5, 502, 505);


  histoRatio1a->DrawCopy();

  TLine *lineA3=new TLine(0.,1.,25.,1.);
  lineA3->SetLineColor(kGray+1);
  lineA3->Draw("same");

  
  graphRatiotoCombinedPPRef->Draw("p,same,e");
  graphRatiotoCombinedPPRefBLUEMethod->Draw("p,same,e");
 
  
  c1a->Print(Form("%s/Pi0_pp_Ref.%s",outputDir.Data(),suffix.Data()));
  
  
  
  
  TCanvas* c1b = new TCanvas("c1b","",200,10,1200,1200);  // gives the page size
  DrawGammaCanvasSettings(c1b , 0.3, 0.02, 0.02, 0.16);
  TPad* pad1b = new TPad("pad1b", "", 0., 0.25, 1., 1.,-1, -1, -2);
  DrawGammaPadSettings( pad1b, 0.15, 0.02, 0.02, 0.);
  pad1b->Draw();

  TPad* padRatios1b = new TPad("padRatios1b", "", 0., 0., 1., 0.25,-1, -1, -2);
  DrawGammaPadSettings( padRatios1b, 0.15, 0.02, 0., 0.2);
  padRatios1b->Draw();

  pad1b->cd();
  pad1b->SetLogx();
  pad1b->SetLogy();
  TH2F * histo1b;
  histo1b = new TH2F("histo1b","histo1b",1000,0.2,30.,1000,2e2,10e11);
  SetStyleHistoTH2ForGraphs(histo1b,"#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.04,0.05, 0.045,0.045, 0.7,1.4, 512, 505);
  histo1b->DrawCopy();



  DrawGammaSetMarkerTGraphAsym(graphCombEtappreferenceSystErr,20,0.7, 4, 4, 1, kTRUE);
  graphCombEtappreferenceSystErr->Draw("E2,same");
  DrawGammaSetMarkerTGraphAsym(graphCombEtappreferenceStatErr,20,0.7, kBlue, kBlue);
  graphCombEtappreferenceStatErr->Draw("p,same");

  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionSyspp5023GeV,20,0.7,kRed+1 ,kRed+1 , 1, kTRUE);
  graphCombEtaInvCrossSectionSyspp5023GeV->Draw("E2,same");
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionStatpp5023GeV,20,0.7, kRed+1, kRed+1);
  graphCombEtaInvCrossSectionStatpp5023GeV->Draw("p,same");

  TLegend* legend1b = new TLegend(0.23,0.10,0.7,0.30);
  legend1b->SetFillColor(0);
  legend1b->SetLineColor(0);
  legend1b->SetTextSize(0.03);
  legend1b->AddEntry(graphCombEtappreferenceSystErr,"Pi0 pp reference","pef");
  legend1b->AddEntry(graphCombEtaInvCrossSectionSyspp5023GeV,"Pi0 pp reference with BLUE Method","pef");
  
  TLatex * lt1b = new TLatex(3,1,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV") ;
  lt1b->SetTextColor(kBlack) ;
  lt1b->SetTextSize(0.04) ;
  lt1b->Draw("same");
  legend1b->Draw();

  padRatios1b->cd();

  padRatios1b->SetLogx();

  DrawGammaSetMarkerTGraphErr(graphRatioEtaPPRef,20,1.0, kBlue, kBlue, 1, kTRUE);
  
  TH2F * histoRatio1b;
  histoRatio1b = new TH2F("histoRatio1b","histo2Dppreference",1000,.2,30.,1000,0.81,1.19);
  SetStyleHistoTH2ForGraphs(histoRatio1b, "#it{p}_{T} (GeV/#it{c})","pp ref/pp ref BLUE", 0.13,0.13, 0.13,0.13, 0.5,.5, 502, 505);


  histoRatio1b->DrawCopy();

  TLine *lineA1b=new TLine(0.,1.,25.,1.);
  lineA1b->SetLineColor(kGray+1);
  lineA1b->Draw("same");

  
  graphRatioEtaPPRef->Draw("p,same,e");
 
  
  c1b->Print(Form("%s/Eta_pp_Ref.%s",outputDir.Data(),suffix.Data()));
  
  
  
  
  
  
  
  
  
  TCanvas* c1c = new TCanvas("c1c","",200,10,1000,1200);  // gives the page size
  DrawGammaCanvasSettings( c1c, 0.3, 0.02, 0.02, 0.16);
  TPad* pad1c = new TPad("pad1c", "", 0., 0.25, 1., 1.,-1, -1, -2);
  DrawGammaPadSettings( pad1c, 0.18, 0.02, 0.02, 0.);
  pad1c->Draw();

  TPad* padRatios1c = new TPad("padRatios1c", "", 0., 0., 1., 0.25,-1, -1, -2);
  DrawGammaPadSettings( padRatios1c, 0.18, 0.02, 0., 0.3);
  padRatios1c->Draw();

  pad1c->cd(); 
  pad1c->SetLogx();
  pad1c->SetLogy();
  
  TH2F* histo1c = new TH2F("histo1c","histo1c",1000,0.2,30.,1000,1e02,2e11 );
  SetStyleHistoTH2ForGraphs(histo1c, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2} ",0.040,0.040, 0.040,0.040, 0.8,1.9, 512, 510);
  histo1c->GetXaxis()->SetLabelOffset(-0.009);
  histo1c->DrawCopy(); 

  
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionSystErrPPRef,20,1,kBlue+2 ,kBlue+2, 1, kTRUE);
  graphCombEtaInvCrossSectionSystErrPPRef->Draw("E2,same"); 
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionStatErrPPRef,20,1,kBlue+2 ,kBlue+2  ); 
  graphCombEtaInvCrossSectionStatErrPPRef->Draw("pz,same"); 
  
  
  
  fitTsallisEtappRef5023GeVPt->SetLineWidth(2); 

  
  fitTsallisEtappRef5023GeVPt->Draw("same");
  
  
  
  TLegend* leg1c;
  leg1c   = new TLegend(0.25,0.31,0.45,0.41); 
  leg1c->SetFillColor(0);
  leg1c->SetLineColor(0);
  leg1c->SetTextFont(42);
  leg1c->SetTextSize(0.045-0.005);
  leg1c->AddEntry(graphCombEtaInvCrossSectionSystErrPPRef,"#eta pp reference","pef");
  leg1c->AddEntry(fitTsallisEtappRef5023GeVPt,"Tsallis fit","l");
  leg1c->Draw();
  
  
  
 	
  TLatex * lt1c = new TLatex(4.5,0.3,"ALICE") ;
  lt1c->SetTextColor(kBlack) ;
  lt1c->SetTextSize(0.045-0.005) ;
  lt1c->SetTextFont(42) ;
  lt1c->DrawLatex(4.5,0.08,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  lt1c->Draw() ;


  padRatios1c->cd();
  padRatios1c->SetLogx();

  
 
  
  
  TH2F * histoRatio1c;
  histoRatio1c = new TH2F("histoRatio1c","histo2DRatioPi0PPRefToFit",1000,.2,30.,1000,0.5,1.5);//1.89
  
  SetStyleHistoTH2ForGraphs(histoRatio1c, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.125,0.125, 0.125,0.125,1.,0.5, 502, 505); 
  histoRatio1c->GetYaxis()->SetLabelOffset(0.005);
  histoRatio1c->GetXaxis()->SetLabelOffset(0.015+0.00);
  histoRatio1c->GetXaxis()->SetTickLength(0.07);
  histoRatio1c->DrawCopy();
 
  
  
 	
  

  TLine *lineRatio1c=new TLine(0.,1.,30.,1.);
  lineRatio1c->SetLineColor(kGray+1);
  lineRatio1c->Draw("same");

  
  
  
  
  DrawGammaSetMarkerTGraphAsym(graphRatioEtappRefFitSystErr,20,1,kBlue+2 ,kBlue+2, 1, kTRUE);
  graphRatioEtappRefFitSystErr->Draw("E2,same"); 
  DrawGammaSetMarkerTGraphAsym(graphRatioEtappRefFitStatErr,20,1,kBlue+2 ,kBlue+2  ); 
  graphRatioEtappRefFitStatErr->Draw("pz,same"); 
  
  
	
  
  c1c->Update();
  c1c->Print(Form("%s/Eta_pp_Ref_Ratio_To_Fit_Tsallis.%s",outputDir.Data(),suffix.Data()));
  
  
  
  
  TCanvas* c1cTCM = new TCanvas("c1cTCM","",200,10,1000,1200);  // gives the page size
  DrawGammaCanvasSettings( c1cTCM, 0.3, 0.02, 0.02, 0.16);
  TPad* pad1cTCM = new TPad("pad1cTCM", "", 0., 0.25, 1., 1.,-1, -1, -2);
  DrawGammaPadSettings( pad1cTCM, 0.18, 0.02, 0.02, 0.);
  pad1cTCM->Draw();

  TPad* padRatios1cTCM = new TPad("padRatios1cTCM", "", 0., 0., 1., 0.25,-1, -1, -2);
  DrawGammaPadSettings( padRatios1cTCM, 0.18, 0.02, 0., 0.3);
  padRatios1cTCM->Draw();

  pad1cTCM->cd(); 
  pad1cTCM->SetLogx();
  pad1cTCM->SetLogy();
  
  TH2F* histo1cTCM = new TH2F("histo1cTCM","histo1cTCM",1000,0.2,30.,1000,1e02,2e11 );
  SetStyleHistoTH2ForGraphs(histo1cTCM, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2} ",0.040,0.040, 0.040,0.040, 0.8,1.9, 512, 510);
  histo1cTCM->GetXaxis()->SetLabelOffset(-0.009);
  histo1cTCM->DrawCopy(); 

  
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionSystErrPPRef,20,1,kBlue+2 ,kBlue+2, 1, kTRUE);
  graphCombEtaInvCrossSectionSystErrPPRef->Draw("E2,same"); 
  DrawGammaSetMarkerTGraphAsym(graphCombEtaInvCrossSectionStatErrPPRef,20,1,kBlue+2 ,kBlue+2  ); 
  graphCombEtaInvCrossSectionStatErrPPRef->Draw("pz,same"); 
  
  
  
  fitTCMEtappRef5023GeVPt->SetLineWidth(2);   
  fitTCMEtappRef5023GeVPt->Draw("same");
  
  
  
  TLegend* leg1cTCM;
  leg1cTCM   = new TLegend(0.25,0.31,0.45,0.41); 
  leg1cTCM->SetFillColor(0);
  leg1cTCM->SetLineColor(0);
  leg1cTCM->SetTextFont(42);
  leg1cTCM->SetTextSize(0.045-0.005);
  leg1cTCM->AddEntry(graphCombEtaInvCrossSectionSystErrPPRef,"#eta pp reference","pef");
  leg1cTCM->AddEntry(fitTCMEtappRef5023GeVPt,"TCM fit","l");
  leg1cTCM->Draw();
  
  
  
 	
  TLatex * lt1cTCM = new TLatex(4.5,0.3,"ALICE") ;
  lt1cTCM->SetTextColor(kBlack) ;
  lt1cTCM->SetTextSize(0.045-0.005) ;
  lt1cTCM->SetTextFont(42) ;
  lt1cTCM->DrawLatex(4.5,0.08,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
  lt1cTCM->Draw() ;


  padRatios1cTCM->cd();
  padRatios1cTCM->SetLogx();

  
 
  
  
  TH2F * histoRatio1cTCM;
  histoRatio1cTCM = new TH2F("histoRatio1cTCM","histo2DRatioPi0PPRefToFit",1000,.2,30.,1000,0.5,1.5);//1.89
  
  SetStyleHistoTH2ForGraphs(histoRatio1cTCM, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.125,0.125, 0.125,0.125,1.,0.5, 502, 505); 
  histoRatio1cTCM->GetYaxis()->SetLabelOffset(0.005);
  histoRatio1cTCM->GetXaxis()->SetLabelOffset(0.015+0.00);
  histoRatio1cTCM->GetXaxis()->SetTickLength(0.07);
  histoRatio1cTCM->DrawCopy();
 
  
  
 	
  

  TLine *lineRatio1cTCM=new TLine(0.,1.,30.,1.);
  lineRatio1cTCM->SetLineColor(kGray+1);
  lineRatio1cTCM->Draw("same");

  
  
  
  
  DrawGammaSetMarkerTGraphAsym(graphRatioEtappRefFitTCMSystErr,20,1,kBlue+2 ,kBlue+2, 1, kTRUE);
  graphRatioEtappRefFitTCMSystErr->Draw("E2,same"); 
  DrawGammaSetMarkerTGraphAsym(graphRatioEtappRefFitTCMStatErr,20,1,kBlue+2 ,kBlue+2  ); 
  graphRatioEtappRefFitTCMStatErr->Draw("pz,same"); 
  
  
	
   
  
  c1cTCM->Update();
  c1cTCM->Print(Form("%s/Eta_pp_Ref_Ratio_To_Fit_TCM.%s",outputDir.Data(),suffix.Data()));
  
  
  
  
  
  
  
  
  
  
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
  // ************************* Plotting Combined R_pPb ************************
  // **************************************************************************************
  TCanvas* canvasCombEtaRpPb = new TCanvas("canvasCombEtaRpPb","",200,10,1200,700);  // gives the page size
  DrawGammaCanvasSettings( canvasCombEtaRpPb, 0.09, 0.02, 0.02, 0.13);
 
  TH2F * histo2DCombinedEtaRpPb;
  histo2DCombinedEtaRpPb = new TH2F("histo2DCombinedEtaRpPb","histo2DCombinedEtaRpPb",1000,0.,20.,1000,0.3,1.7);
  SetStyleHistoTH2ForGraphs(histo2DCombinedEtaRpPb, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#eta}_{p-Pb}", 0.05,0.06, 0.05,0.06, 0.9,0.6, 512, 505);
  histo2DCombinedEtaRpPb->DrawCopy();
    

     
    

  DrawGammaSetMarkerTGraphAsym(graphCombEtaRpPb5023GeVSystClone,20,1.5, 4, 4, 1, kTRUE);  
      
  graphCombEtaRpPb5023GeVSystClone->Draw("E2,same");

  DrawGammaSetMarkerTGraphAsym(graphCombEtaRpPb5023GeVStatClone,20,1.5, kBlue, kBlue);  
  graphCombEtaRpPb5023GeVStatClone->Draw("p,same");


 // TLatex * lt = new TLatex(1.,1.55,"ALICE Preliminary"); 
 // lt->SetTextSize(0.05) ;
  //  lt->DrawText(1.,1.45,"2016/06/14") ;
  lt->Draw("same");
     
  TLegend* legendEtaRpPbCombine = new TLegend(0.1,0.75,0.55,0.85);
  legendEtaRpPbCombine->SetFillColor(0);
  legendEtaRpPbCombine->SetLineColor(0);
  legendEtaRpPbCombine->SetTextSize(0.04);
  legendEtaRpPbCombine->AddEntry(graphInvYieldPi0CombpPb5023GeVSysClone,"NSD, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");

    
  legendEtaRpPbCombine->Draw();
  
  TLine* lineEtaRpPb =new TLine(0.,1.,20.,1.);
  lineEtaRpPb->Draw("same");
  
  TBox* BoxNormEtaRpPb =new TBox(0.25, 1.-NormalizationError, 0.5, 1.+NormalizationError);
  BoxNormEtaRpPb->SetFillColor(4);
  BoxNormEtaRpPb->Draw("same");
    
  canvasCombEtaRpPb->Print(Form("%s/Comb_Eta_RpPb.%s",outputDir.Data(),suffix.Data()));
  
  
  
  
  
  
  
  
  
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

    graphCombPi0RpPbSystErr=ScaleGraph(graphCombPi0RpPbSystErr,Scaling);
    graphCombPi0RpPbStatErr=ScaleGraph(graphCombPi0RpPbStatErr,Scaling);
    
    //  graphCombPi0RpPbSystErr=ApplyNSDSysError(graphCombPi0RpPbSystErr,ScalingErr,0.);  
    
  }



  TCanvas* canvasCombRpPb2 = new TCanvas("canvasCombRpPb2","",200,10,1200,700);  // gives the page size
  DrawGammaCanvasSettings( canvasCombRpPb2, 0.08, 0.02, 0.02, 0.13);
    
  //  canvasCombRpPb2->SetLogx();
  // canvasCombRpPb2->SetLogy();
  TH2F * histo2DCombined2;
  histo2DCombined2 = new TH2F("histo2DCombined2","histo2DCombined2",1000,0.,20.,1000,0.3,1.7);
  SetStyleHistoTH2ForGraphs(histo2DCombined2, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#pi^{0}}_{p-Pb}", 0.05,0.06, 0.05,0.05, 0.9,0.7, 512, 505);
  histo2DCombined2->DrawCopy();
    

     
    

  DrawGammaSetMarkerTGraphAsym(graphCombPi0RpPbSystErr,20,1.5, 4, 4, 1, kTRUE);  
      
  graphCombPi0RpPbSystErr->Draw("E2,same");

  DrawGammaSetMarkerTGraphAsym(graphCombPi0RpPbStatErr,20,1.5, kBlue, kBlue);  
  graphCombPi0RpPbStatErr->Draw("p,same");


  TLatex * lt3 = new TLatex(1.,1.55,"ALICE Preliminary"); 
  lt3->SetTextSize(0.05) ;
  //  lt3->DrawText(1.,1.45,"2016/06/14") ;
  lt3->Draw("same");
     
  TLegend* legendRpPbCombine2 = new TLegend(0.1,0.75,0.55,0.85);
  legendRpPbCombine2->SetFillColor(0);
  legendRpPbCombine2->SetLineColor(0);
  legendRpPbCombine2->SetTextSize(0.04);
  legendRpPbCombine2->AddEntry(graphCombPi0RpPbSystErr,"p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","pef");

    
  legendRpPbCombine->Draw();
  
  TLine* line2 =new TLine(0.,1.,20.,1.);
  line2->Draw("same");

  BoxNorm->Draw("same");

  canvasCombRpPb2->Print(Form("%s/RpPb_CombinedpPbSpectrum.%s",outputDir.Data(),suffix.Data()));
    
    
  
  fFits.close();

  //************************ WriteResultsToFile **************************************************
  TFile fResults(Form("%s/ResultsRpPbpPb_%s.root",outputDir.Data(), dateForOutput.Data()),"RECREATE");
     
  graphInvYieldPi0CombpPb5023GeVSysClone->Write("CombinedPi0RpPbSystErr");
  graphInvYieldPi0CombpPb5023GeVStaClone->Write("CombinedPi0RpPbStatErr");
  graphCombEtaRpPb5023GeVSystClone->Write("CombinedEtaRpPbSystErr");
  graphCombEtaRpPb5023GeVStatClone->Write("CombinedEtaRpPbSystErr");
  
   
  
  
  fResults.Close();

  TFile fPaperPlots1(Form("%s/PaperPlotsRpPb_%s.root",outputDir.Data(), dateForOutput.Data()),"RECREATE");
  
  
     
  graphInvYieldPi0CombpPb5023GeVSysClone->Write("CombinedPi0RpPbSystErr");
  graphInvYieldPi0CombpPb5023GeVStaClone->Write("CombinedPi0RpPbStatErr");
  graphCombEtaRpPb5023GeVSystClone->Write("CombinedEtaRpPbSystErr");
  graphCombEtaRpPb5023GeVStatClone->Write("CombinedEtaRpPbStatErr");
  
  graphCombEtaInvCrossSectionSystErrPPRef->Write("EtaPPReferenceSystErr");
  graphCombEtaInvCrossSectionStatErrPPRef->Write("EtaPPReferenceStatErr");
  graphCombPi0InvCrossSectionSystErrPPRef->Write("Pi0PPReferenceSystErr");
  graphCombPi0InvCrossSectionStatErrPPRef->Write("Pi0PPReferenceStatErr");
  fitTsallisEtappRef5023GeVPt->Write("EtaPPReferenceTSallisFit");
  fitTCMEtappRef5023GeVPt->Write("EtaPPReferenceTCMFit");
  fitTsallisPi0ppRef5023GeVPt->Write("Pi0PPReferenceTSallisFit");
  fitTCMPi0ppRef5023GeVPt->Write("Pi0PPReferenceTCMFit");
  
 
  graphRatioEtappRefFitTotErr->Write("graphRatioEtappRefTsallisFitTotErr");
  graphRatioEtappRefFitStatErr->Write("graphRatioEtappRefTsallisFitStatErr");   
  graphRatioEtappRefFitSystErr->Write("graphRatioEtappRefTsallisFitSystErr"); 


  graphRatioEtappRefFitTCMTotErr->Write("graphRatioEtappRefTCMFitTotErr");
  graphRatioEtappRefFitTCMStatErr->Write("graphRatioEtappRefTCMFitStatErr");   
  graphRatioEtappRefFitTCMSystErr->Write("graphRatioEtappRefTCMFitSystErr");   
 
  graphRatioPi0ppRefFitTCMTotErr->Write("graphRatioPi0ppRefTCMFitTotErr");     
  graphRatioPi0ppRefFitTCMStatErr->Write("graphRatioPi0ppRefTCMFitStatErr");     
  graphRatioPi0ppRefFitTCMSystErr->Write("graphRatioPi0ppRefTCMFitSystErr"); 
  
  graphRatioPi0ppRefFitTotErr->Write("graphRatioPi0ppRefTsallisFitTotErr");     
  graphRatioPi0ppRefFitStatErr->Write("graphRatioPi0ppRefTsallisFitStatErr");     
  graphRatioPi0ppRefFitSystErr->Write("graphRatioPi0ppRefTsallisFitSystErr"); 
  
  
  
  
  
  
  
  
  graphPi0ESP09sPi0KKP->Write("EPS09s_KKP_NLO");
  graphPi0ESP09sPi0AKK->Write("EPS09s_AKK_NLO");
  graphPi0DSS5000->Write("EPS09s_fDSS_NLO");
  graphAsymmErrorsPi0DSS5000->Write("EPS09s_fDSS_errors");  
  graphPi0CGC->Write("CGC");
  graphCombppreferenceStatErr->Write("ppRefStat");
  graphCombppreferenceSystErr->Write("ppRefSyst");
  
  fPaperPlots1.Close();
  
  TFile fPaperPlots2(fileOutputResultspPb.Data(),"UPDATE");
  
  fPaperPlots2.mkdir("Pi0RpPb");
  fPaperPlots2.mkdir("EtaRpPb");
  fPaperPlots2.cd("Pi0RpPb");
  
  graphInvYieldPi0CombpPb5023GeVSysClone->Write("CombinedPi0RpPbSystErr");
  graphInvYieldPi0CombpPb5023GeVStaClone->Write("CombinedPi0RpPbStatErr");
  
  
  if( statErrorCollection[0] && sysErrorCollection[0] ){
      statErrorCollection[0]->Write("Pi0_RpPb_PCM_StatErr");
      sysErrorCollection[0]->Write("Pi0_RpPb_PCM_SystErr");
  }
  if( statErrorCollection[1] && sysErrorCollection[1] ) {
      statErrorCollection[1]->Write("Pi0_RpPb_PHOS_StatErr");
      sysErrorCollection[1]->Write("Pi0_RpPb_PHOS_SystErr");
  }
  if( statErrorCollection[2] && sysErrorCollection[2] ) {
      statErrorCollection[2]->Write("Pi0_RpPb_EMCAL_StatErr");
      sysErrorCollection[2]->Write("Pi0_RpPb_EMCAL_SystErr");
  }
  if( statErrorCollection[4] && sysErrorCollection[4] ) {
      statErrorCollection[4]->Write("Pi0_RpPb_PCMEMCAL_StatErr");
      sysErrorCollection[4]->Write("Pi0_RpPb_PCMEMCAL_SystErr");
  }
  if( statErrorCollection[5] && sysErrorCollection[5] ) {
      statErrorCollection[5]->Write("Pi0_RpPb_Dalitz_StatErr");
      sysErrorCollection[5]->Write("Pi0_RpPb_Dalitz_SystErr");
  }
  
  
  graphPCMppreferenceStatErr->Write("Pi0_pp_reference_PCMBinning_StatErr");
  graphPCMppreferenceSystErr->Write("Pi0_pp_reference_PCMBinning_SystErr");
  graphPHOSppreferenceStatErr->Write("Pi0_pp_reference_PHOSBinning_StatErr");
  graphPHOSppreferenceSystErr->Write("Pi0_pp_reference_PHOSBinning_SystErr");
  graphEMCalppreferenceStatErr->Write("Pi0_pp_reference_EMCALBinning_StatErr");
  graphEMCalppreferenceSystErr->Write("Pi0_pp_reference_EMCALBinning_SystErr");
  graphPCMEMCalppreferenceStatErr->Write("Pi0_pp_reference_PCMEMCALBinning_StatErr");
  graphPCMEMCalppreferenceSystErr->Write("Pi0_pp_reference_PCMEMCALBinning_SystErr");
  graphDalitzppreferenceStatErr->Write("Pi0_pp_reference_DalitzBinning_StatErr");
  graphDalitzppreferenceSystErr->Write("Pi0_pp_reference_DalitzBinning_SystErr");
  
  

  
  graphInvYieldPi0PCMpPb5023GeVStatErryShifted->Write("graphInvYieldPi0PCMpPb5023GeVStatErr_yShifted_NSD");
  graphInvYieldPi0PCMpPb5023GeVSystErryShifted->Write("graphInvYieldPi0PCMpPb5023GeVSystErr_yShifted_NSD");
  
  graphInvYieldPi0PHOSpPb5023GeVStatErryShifted->Write("graphInvYieldPi0PHOSpPb5023GeVStatErr_yShifted_NSD");
  graphInvYieldPi0PHOSpPb5023GeVSystErryShifted->Write("graphInvYieldPi0PHOSpPb5023GeVSystErr_yShifted_NSD");
  
  graphInvYieldPi0EMCalpPb5023GeVStatErryShifted->Write("graphInvYieldPi0EMCALpPb5023GeVStatErr_yShifted_NSD");
  graphInvYieldPi0EMCalpPb5023GeVSystErryShifted->Write("graphInvYieldPi0EMCALpPb5023GeVSystErr_yShifted_NSD");
  
  graphInvYieldPi0PCMEMCalpPb5023GeVStatErryShifted->Write("graphInvYieldPi0PCMEMCALpPb5023GeVStatErr_yShifted_NSD");
  graphInvYieldPi0PCMEMCalpPb5023GeVSystErryShifted->Write("graphInvYieldPi0PCMEMCALpPb5023GeVSystErr_yShifted_NSD");
  
  
  graphInvYieldPi0DalitzpPb5023GeVStatErryShifted->Write("graphInvYieldPi0DalitzpPb5023GeVStatErr_yShifted_NSD");
  graphInvYieldPi0DalitzpPb5023GeVSystErryShifted->Write("graphInvYieldPi0DalitzpPb5023GeVSystErr_yShifted_NSD");
  
  
  
  fPaperPlots2.cd("EtaRpPb");
  
  graphCombEtaRpPb5023GeVSystClone->Write("CombinedEtaRpPbSystErr");
  graphCombEtaRpPb5023GeVStatClone->Write("CombinedEtaRpPbStatErr");
  
  if( statErrorCollectionEtaRpPb[0] && statErrorCollectionEtaRpPb[0] ) {
      statErrorCollectionEtaRpPb[0]->Write("Eta_RpPb_PCM_StatErr");
      sysErrorCollectionEtaRpPb[0]->Write("Eta_RpPb_EMCal_SystErr");
  }
  if( statErrorCollectionEtaRpPb[2] && sysErrorCollectionEtaRpPb[2] ){
      statErrorCollectionEtaRpPb[2]->Write("Eta_RpPb_EMCAL_StatErr");
      sysErrorCollectionEtaRpPb[2]->Write("Eta_RpPb_EMCAL_SystErr");
  }
  if( statErrorCollectionEtaRpPb[4] && sysErrorCollectionEtaRpPb[4] ){
      statErrorCollectionEtaRpPb[4]->Write("Eta_RpPb_PCMEMCAL_StatErr");
      sysErrorCollectionEtaRpPb[4]->Write("Eta_RpPb_PCMEMCAL_SystErr");
      
    
  }
  
   
  
  
  graphPCMEtappreferenceStatErr->Write("Eta_pp_reference_PCMBinning_StatErr");
  graphPCMEtappreferenceSystErr->Write("Eta_pp_reference_PCMBinning_SystErr");
  graphEMCalEtappreferenceStatErr->Write("Eta_pp_reference_EMCALBinning_StatErr");
  graphEMCalEtappreferenceSystErr->Write("Eta_pp_reference_EMCALBinning_SystErr");
  graphPCMEMCalEtappreferenceStatErr->Write("Eta_pp_reference_PCMEMCALBinning_StatErr");
  graphPCMEMCalEtappreferenceSystErr->Write("Eta_pp_reference_PCMEMCALBinning_SystErr");
  

  
  graphInvYieldEtaPCMpPb5023GeVStatErryShifted->Write("graphInvYieldEtaPCMpPb5023GeVStatErr_yShifted_NSD");
  graphInvYieldEtaPCMpPb5023GeVSystErryShifted->Write("graphInvYieldEtaPCMpPb5023GeVSystErr_yShifted_NSD");

  graphInvYieldEtaEMCalpPb5023GeVStatErryShifted->Write("graphInvYieldEtaEMCALpPb5023GeVStatErr_yShifted_NSD");
  graphInvYieldEtaEMCalpPb5023GeVSystErryShifted->Write("graphInvYieldEtaEMCALpPb5023GeVSystErr_yShifted_NSD");
  
  graphInvYieldEtaPCMEMCalpPb5023GeVStatErryShifted->Write("graphInvYieldEtaPCMEMCALpPb5023GeVStatErr_yShifted_NSD");
  graphInvYieldEtaPCMEMCalpPb5023GeVSystErryShifted->Write("graphInvYieldEtaPCMEMCALpPb5023GeVSystErr_yShifted_NSD");
  
  
  fPaperPlots2.Close();
  
 


}
    
