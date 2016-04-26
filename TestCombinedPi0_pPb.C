/****************************************************************************************************************************
******      provided by Gamma Conversion Group, GA,                                        *****
******      Ana Marin, marin@physi.uni-heidelberg.de                                      *****
******      Friederike Bock, friederike.bock@cern.ch                                      *****
******      Annika Passfeld    
******      Pedro Gonzalez 
******      25.04.2014
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
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"

extern TRandom*   gRandom;
extern TBenchmark*   gBenchmark;
extern TSystem*   gSystem;
extern TMinuit*   gMinuit;


void TestCombinedPi0_pPb(TString outputDir = "Combined_pPb/", TString suffix = "pdf"){


// Style_t 	markerStyleConv		 	= 20 ;
// Style_t 	markerStyleMCEtaToPi0		= 24 ;
// Style_t 	markerStylePHOS		 	= 21 ;
// Style_t 	markerStyleEMCAL		 	= 33 ;
// Style_t 	markerStyleConvMC			= 24 ;
// Style_t 	markerStylePHOSMC		 	= 25 ;
// Style_t 	markerStyleEMCALMC		 	= 27 ;
// Color_t 	colorConv 				= kBlack;
// Color_t 	colorConvMC 				= kGray+1;
// Color_t	colorEMCAL				= kGreen+2;
// Color_t	colorEMCALMC				= kGreen-6;
// Color_t 	colorPHOS					= kBlue+1;
// Color_t 	colorPHOSMC				= kBlue-7;
// Size_t 	markerSizeInvYield			= 1.5;
// Size_t 	markerSizeMass				= 3.2;
// Size_t 	markerSizeMassMC			= 1.5;
// Size_t       markerSizeSpectrum                      = 1.;


  gROOT->Reset();   
  gROOT->SetStyle("Plain");
  
  StyleSettingsThesis();  
  SetPlotStyle();
  gSystem->Exec("mkdir -p "+outputDir);

  gStyle->SetOptTitle(kFALSE);

  TString nameHistoPCM = "CorrectedYieldPi0";
  TString nameGraphPCM = "Pi0SystError";

  
  Width_t         widthLinesBoxes;
	widthLinesBoxes				= 2.3;
// 		widthCommonFit					= 2.6;
// 		widthStatErrBars				= 2.6;
// 		widthCommonErrors				= 2.;
// 		widthCommonSpectrumBoxes			= 2.3;

	TFile* fileNeutralPionPCMDatapPb = new TFile("ExternalInputpPb/PCM/data_PCMResults_pPb_20141023.root");
	TFile* fileNeutralPionDalitzpPb = new TFile("ExternalInputpPb/PCM/data_PCMResults_Dalitz_pPb_2014-09-12-RCut5cm.root");
	TFile*   filePHOSpPb =       new TFile("ExternalInputpPb/PHOS/data_PHOSResultsFullCorrection_pPb-23102014.root");
  TFile*   filePHOSpPbSys =       new TFile("ExternalInputpPb/PHOS/sys_err_PHOS_20140130.root");
  TFile *  fileEMCalpPb = new TFile("ExternalInputpPb/EMCAL/EMCALResults_11Sept.root"); 
  

  cout<< "PCM::"<<endl;
  //****************************************************************************************************
  //************************** Read data for PCM *******************************************************
  //****************************************************************************************************
  TDirectory* directoryPCMPi0pPb =  (TDirectory*)fileNeutralPionPCMDatapPb->Get("Pi0_pPb_5.023TeV_0-100%"); 
  TH1D* histoPCMYieldPi0pPb =       (TH1D*)directoryPCMPi0pPb->Get(nameHistoPCM.Data());
  TGraphAsymmErrors* graphPCMYieldPi0SysErrpPb= (TGraphAsymmErrors*)directoryPCMPi0pPb->Get(nameGraphPCM.Data());

  TGraphAsymmErrors* graphRatioCombPCMSys=(TGraphAsymmErrors*)graphPCMYieldPi0SysErrpPb->Clone();;
  TGraphAsymmErrors* graphRatioCombPCMSta;


  Int_t nP_PCM = graphPCMYieldPi0SysErrpPb->GetN();
  Double_t * eyieldPi0PCMSyst = graphPCMYieldPi0SysErrpPb->GetEYlow();
  Double_t * ptPi0PCM =   graphPCMYieldPi0SysErrpPb->GetX();
  Double_t yieldPi0PCM[30];
  Double_t eyieldPi0PCMStat[30];
  Double_t eyieldPi0PCMTot[30];

  //  Double_t eyieldPi0PCMSyst[30];
  cout<< nP_PCM<<endl;
  for(Int_t ii=1;ii<histoPCMYieldPi0pPb->GetNbinsX();ii++){   // first bin is eempty for the histo
    yieldPi0PCM[ii-1]=histoPCMYieldPi0pPb->GetBinContent(ii+1);
    eyieldPi0PCMStat[ii-1]=histoPCMYieldPi0pPb->GetBinError(ii+1);
    eyieldPi0PCMTot[ii-1]=pow((pow(eyieldPi0PCMStat[ii-1],2)+pow(eyieldPi0PCMSyst[ii-1],2)),0.5);
    // cout<< ii<<" "<<histoPCMYieldPi0pPb->GetBinCenter(ii+1)<<" "<< histoPCMYieldPi0pPb->GetBinContent(ii+1)<<" "<< histoPCMYieldPi0pPb->GetBinError(ii+1) <<" "<<ptPi0PCM[ii-1]<<"  "<< eyieldPi0PCMTot[ii-1]<< " " <<eyieldPi0PCMStat[ii-1]<<" "<< eyieldPi0PCMSyst[ii-1]<<endl;

  }
  
  graphPCMYieldPi0SysErrpPb->Print();
  

  cout<< "Dalitz::"<<endl;
  //****************************************************************************************************
  //************************** Read data for Dalitz *******************************************************
  //****************************************************************************************************  
  
  TDirectory* directoryDalitzPi0pPb = (TDirectory*)fileNeutralPionDalitzpPb->Get("Pi0_pPb_5.023TeV_0-100%"); 
  TH1D* histoDalitzYieldPi0pPb =   (TH1D*)directoryDalitzPi0pPb->Get(nameHistoPCM.Data());
  TGraphAsymmErrors* graphDalitzYieldPi0SysErrpPb=(TGraphAsymmErrors*)directoryDalitzPi0pPb->Get(nameGraphPCM.Data()); 

  TGraphAsymmErrors* graphRatioCombDalitzSys= (TGraphAsymmErrors*)graphDalitzYieldPi0SysErrpPb->Clone();
  TGraphAsymmErrors* graphRatioCombDalitzSta;




  graphDalitzYieldPi0SysErrpPb->Print();


  //  Double_t ptPi0Dalitz[21];
  //Double_t yieldPi0Dalitz[21];
  //Double_t eyieldPi0DalitzStat[21];
  // 



  Int_t nP_Dalitz = graphDalitzYieldPi0SysErrpPb->GetN();
  Double_t * eyieldPi0DalitzSyst = graphDalitzYieldPi0SysErrpPb->GetEYlow();
  Double_t * ptPi0Dalitz =   graphDalitzYieldPi0SysErrpPb->GetX();
  Double_t yieldPi0Dalitz[21];
  Double_t eyieldPi0DalitzStat[21];
  Double_t eyieldPi0DalitzTot[21];

  //  Double_t eyieldPi0PCMSyst[30];
  cout<< "npdalitz::"<<nP_Dalitz<<endl;
  for(Int_t ii=1;ii<histoDalitzYieldPi0pPb->GetNbinsX();ii++){   // first bin is eempty for the histo
    cout<< ii<<" "<<histoDalitzYieldPi0pPb->GetBinCenter(ii+1)<< " " <<histoDalitzYieldPi0pPb->GetBinContent(ii+1)<<endl;
    yieldPi0Dalitz[ii-1]=histoDalitzYieldPi0pPb->GetBinContent(ii+1);
    eyieldPi0DalitzStat[ii-1]=histoDalitzYieldPi0pPb->GetBinError(ii+1);
    eyieldPi0DalitzTot[ii-1]=pow((pow(eyieldPi0DalitzStat[ii-1],2)+pow(eyieldPi0DalitzSyst[ii-1],2)),0.5);
    cout<<"after"<< ii<<" "<<histoDalitzYieldPi0pPb->GetBinCenter(ii+1)<<" "<< histoDalitzYieldPi0pPb->GetBinContent(ii+1)
    <<" "<< histoDalitzYieldPi0pPb->GetBinError(ii+1) <<" "<<ptPi0Dalitz[ii-1]<<"  "<< eyieldPi0DalitzTot[ii-1]<< " " <<eyieldPi0DalitzStat[ii-1]<<" "<< eyieldPi0DalitzSyst[ii-1]<<endl;
  }


  cout<< "PHOS::"<<endl;
  //****************************************************************************************************
  //************************** Read data for PHOS *******************************************************
  //****************************************************************************************************
  
  TDirectory *directoryPHOSPi0pPb = 	(TDirectory*)filePHOSpPb->Get("Pi0_pPb_5.023TeV_0-100%"); 
  TH1D* histoPHOSYieldPi0pPb = 	(TH1D*)directoryPHOSPi0pPb->Get(nameHistoPCM.Data());      
  
  TH1D* histoPHOSYieldPi0pPbSys = 	(TH1D*)filePHOSpPbSys->Get("yeild1_GS_All_cen0");      


  Double_t ptPi0PHOS[22];
  Double_t yieldPi0PHOS[22];
  Double_t eyieldPi0PHOSStat[22];
  Double_t eyieldPi0PHOSSyst[22];
  Double_t eyieldPi0PHOSTot[22];
  Int_t nP_PHOS=22;
  cout<<"PHOS"<<endl;
  for (int i=1; i<8;i++ ){
    histoPHOSYieldPi0pPbSys->SetBinError(i,0.);
    histoPHOSYieldPi0pPbSys->SetBinContent(i,0.);
    histoPHOSYieldPi0pPb->SetBinError(i,0.);
    histoPHOSYieldPi0pPb->SetBinContent(i,0.);
  }

  TH1D * histoRatioCombPHOS= (TH1D*)histoPHOSYieldPi0pPbSys->Clone();
  for(Int_t ii=7;ii<histoPHOSYieldPi0pPb->GetNbinsX();ii++){
    //    cout<<histoPHOSYieldPi0pPb->GetBinCenter(ii+1)<<" "<< histoPHOSYieldPi0pPb->GetBinContent(ii+1)<<" "<< histoPHOSYieldPi0pPb->GetBinError(ii+1)<<" "<< histoPHOSYieldPi0pPbSys->GetBinError(ii+1)<<endl;
    ptPi0PHOS[ii-7]=histoPHOSYieldPi0pPb->GetBinCenter(ii+1);
    yieldPi0PHOS[ii-7]=histoPHOSYieldPi0pPb->GetBinContent(ii+1);
    eyieldPi0PHOSStat[ii-7]=histoPHOSYieldPi0pPb->GetBinError(ii+1);
    eyieldPi0PHOSSyst[ii-7]=histoPHOSYieldPi0pPbSys->GetBinError(ii+1);
    eyieldPi0PHOSTot[ii-7]=pow((pow(eyieldPi0PHOSStat[ii-7],2)+pow(eyieldPi0PHOSSyst[ii-7],2)),0.5);
    cout<<  ptPi0PHOS[ii-7]<<" "<<  yieldPi0PHOS[ii-7]<<" "<< eyieldPi0PHOSStat[ii-7]<< " "<< eyieldPi0PHOSSyst[ii-7]<<" " << eyieldPi0PHOSTot[ii-7]<<" " <<eyieldPi0PHOSStat[ii-7]/ yieldPi0PHOS[ii-7]<<" "<<  eyieldPi0PHOSSyst[ii-7]/yieldPi0PHOS[ii-7]<<" " << eyieldPi0PHOSTot[ii-7]/yieldPi0PHOS[ii-7]<< endl;
  }
  cout<<endl;
  cout<< "EMCAL::"<<endl;
  //****************************************************************************************************
  //************************** Read data for EMCal *******************************************************
  //****************************************************************************************************
  TString nameHistoEMCal = "CorrectedYieldPi0";
  TString nameHistoEMCalSysErrors = "Pi0SystError";
  
  TDirectory*	directoryEMCalPi0pPb = 	(TDirectory*)fileEMCalpPb->Get("Pi05.02TeV_pPb"); 
  TH1D*	histoEMCalYieldPi0pPb    = (TH1D*)directoryEMCalPi0pPb->Get(nameHistoEMCal.Data());
  TH1D* histoEMCalYieldPi0pPbSys = (TH1D*)directoryEMCalPi0pPb->Get(nameHistoEMCalSysErrors.Data());
  TH1D * histoRatioCombEMCAL = (TH1D*)histoEMCalYieldPi0pPbSys->Clone();
  Double_t ptPi0EMCAL[24];
  Double_t yieldPi0EMCAL[24];
  Double_t eyieldPi0EMCALStat[24];
  Double_t eyieldPi0EMCALSyst[24];
  Double_t eyieldPi0EMCALTot[24];

  Int_t nP_EMCAL=24;

  for (Int_t ii=7;ii< histoEMCalYieldPi0pPbSys->GetNbinsX();ii++){
    //    cout<< histoEMCalYieldPi0pPbSys->GetBinCenter(ii+1)<<" "<<  histoEMCalYieldPi0pPbSys->GetBinContent(ii+1)<<" "<<
    // histoEMCalYieldPi0pPbSys->GetBinError(ii+1)<< " "<< histoEMCalYieldPi0pPb->GetBinError(ii+1) <<endl;
    ptPi0EMCAL[ii-7]=histoEMCalYieldPi0pPb->GetBinCenter(ii+1);
    yieldPi0EMCAL[ii-7]=histoEMCalYieldPi0pPb->GetBinContent(ii+1);
    eyieldPi0EMCALStat[ii-7]=histoEMCalYieldPi0pPb->GetBinError(ii+1);
    eyieldPi0EMCALSyst[ii-7]=histoEMCalYieldPi0pPbSys->GetBinError(ii+1);
    eyieldPi0EMCALTot[ii-7]=pow((pow(eyieldPi0EMCALStat[ii-7],2)+pow(eyieldPi0EMCALSyst[ii-7],2)),0.5);
  }

  Double_t xValComb[34]={0.35, 0.45, 0.55, 0.65, 0.75, 0.9,
			 1.1, 1.3, 1.5, 1.7, 1.9,
			 2.1, 2.3, 2.5, 2.7, 2.9,
			 3.1, 3.3, 3.5, 3.7, 3.9,
			 4.25,4.75,5.25,5.75,6.5, 7.5,9,11,13,15,17,19,22.5};

  Double_t exValComb[34]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0,
      0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0,};

  Double_t yieldPi0Comb[34];
  Double_t totErrPi0Comb[34];
  Double_t staErrPi0Comb[34];
  Double_t sysErrPi0Comb[34];
  Double_t totWeightPi0Comb[34];


  Int_t Ntotal=34;

  Int_t NusedPCM=30;
  Int_t NusedDalitz=15;
  Int_t NusedPHOS=22; 
  Int_t NusedEMCAL=22;

  Int_t NoffsetPCM=0;
  Int_t NoffsetDalitz=2;
  Int_t NoffsetPHOS=12;
  Int_t NoffsetEMCAL=6;

  Double_t weight_PCM[30];
  Double_t weight_Dalitz[15];
  Double_t weight_PHOS[22];
  Double_t weight_EMCAL[22];




  // For the moment assume that PCM-gg and PCM-eeg are independent (statistically they are, not fully systematically)
  // So, PCM-gg, and PCM-eeg EMCAL and PHOS  are taken as fully independent


  for (Int_t ii=0;ii<Ntotal;ii++){
    //    cout<<"ii::"<< ii<<endl;
    if(ii<3){                           //PCM only
      //  cout<<"case 1"<<endl;
      yieldPi0Comb[ii]= yieldPi0PCM[ii];
      totErrPi0Comb[ii]=eyieldPi0PCMTot[ii];
      staErrPi0Comb[ii]=eyieldPi0PCMStat[ii];
      sysErrPi0Comb[ii]=eyieldPi0PCMSyst[ii];


    }else if( ii>=3 && ii<6) {          //PCM+Dalitz       

      //      cout<<"Test 1st interval PCM+Dalitz " <<  xValComb[ii]<<" "<<ptPi0PCM[ii]<<" "<< ptPi0Dalitz[ii-NoffsetDalitz]<<endl;
      weight_PCM[ii]=1./pow(eyieldPi0PCMTot[ii],2);
      weight_Dalitz[ii-NoffsetDalitz]=1./pow(eyieldPi0DalitzTot[ii-NoffsetDalitz],2);

      totWeightPi0Comb[ii]=weight_PCM[ii]+weight_Dalitz[ii-NoffsetDalitz];


      yieldPi0Comb[ii]=(yieldPi0PCM[ii]*weight_PCM[ii]+
			yieldPi0Dalitz[ii-NoffsetDalitz]*weight_Dalitz[ii-NoffsetDalitz])/
	(totWeightPi0Comb[ii]);
      
      totErrPi0Comb[ii]= pow(totWeightPi0Comb[ii],-0.5);
    
      staErrPi0Comb[ii]=pow((1./2.*(weight_PCM[ii]/totWeightPi0Comb[ii])*pow(eyieldPi0PCMStat[ii],2)+
			     1./2.*(weight_Dalitz[ii-NoffsetDalitz]/totWeightPi0Comb[ii])*
			     pow(eyieldPi0DalitzStat[ii],2)),0.5);

      sysErrPi0Comb[ii]=pow((1./2.*(weight_PCM[ii]/totWeightPi0Comb[ii])*pow(eyieldPi0PCMSyst[ii],2)+
			     1./2.*(weight_Dalitz[ii-NoffsetDalitz]/totWeightPi0Comb[ii])*
			     pow(eyieldPi0DalitzSyst[ii],2)),0.5);


    }else if( ii>=6 && ii<12){          //PCM+dalitz+EMCAL
      //      cout<<"Test 2nd interval PCM+Dalitz+EMCAL"<< xValComb[ii]<<" "<<ptPi0PCM[ii]<<" "<< ptPi0Dalitz[ii-NoffsetDalitz]<<" "<< ptPi0EMCAL[ii-NoffsetEMCAL]<<endl;
      weight_PCM[ii]=1./pow(eyieldPi0PCMTot[ii],2);
      weight_Dalitz[ii-NoffsetDalitz]=1./pow(eyieldPi0DalitzTot[ii-NoffsetDalitz],2);
      weight_EMCAL[ii-NoffsetEMCAL]=1./pow(eyieldPi0EMCALTot[ii-NoffsetEMCAL],2);
      
      totWeightPi0Comb[ii]=weight_PCM[ii]+weight_Dalitz[ii-NoffsetDalitz]+weight_EMCAL[ii-NoffsetEMCAL];

      
      yieldPi0Comb[ii]=(yieldPi0PCM[ii]*weight_PCM[ii]+
			yieldPi0Dalitz[ii-NoffsetDalitz]*weight_Dalitz[ii-NoffsetDalitz]+
			yieldPi0EMCAL[ii-NoffsetEMCAL]*weight_EMCAL[ii-NoffsetEMCAL])/totWeightPi0Comb[ii];
      
      
      totErrPi0Comb[ii]= pow(totWeightPi0Comb[ii],-0.5);
      
      staErrPi0Comb[ii]=pow((1./3.*(weight_PCM[ii]/totWeightPi0Comb[ii])*pow(eyieldPi0PCMStat[ii],2)+
			     1./3.*(weight_Dalitz[ii-NoffsetDalitz]/totWeightPi0Comb[ii])*
			     pow(eyieldPi0DalitzStat[ii-NoffsetDalitz],2)+
			     1./3.*(weight_EMCAL[ii-NoffsetEMCAL]/totWeightPi0Comb[ii])*
			     pow(eyieldPi0EMCALStat[ii-NoffsetEMCAL],2)),0.5);

      sysErrPi0Comb[ii]=pow((1./3.*(weight_PCM[ii]/totWeightPi0Comb[ii])*pow(eyieldPi0PCMSyst[ii],2)+
			     1./3.*(weight_Dalitz[ii-NoffsetDalitz]/totWeightPi0Comb[ii])*
			     pow(eyieldPi0DalitzSyst[ii-NoffsetDalitz],2)+
			     1./3.*(weight_EMCAL[ii-NoffsetEMCAL]/totWeightPi0Comb[ii])*
			     pow(eyieldPi0EMCALSyst[ii-NoffsetEMCAL],2)),0.5);

      
    }else if(ii>=12 && ii<17){          // PCM+dalitz+EMCAL+PHOS
      // cout<<"Test 3rd interval PCM+Dalitz+EMCAL+PHOS " <<xValComb[ii]<<" "<<  ptPi0PCM[ii]<< "  "<< ptPi0Dalitz[ii-NoffsetDalitz]<<" "<< ptPi0EMCAL[ii-NoffsetEMCAL]<<" "<< ptPi0PHOS[ii-NoffsetPHOS]<< endl;

      
      weight_PCM[ii]=1./pow(eyieldPi0PCMTot[ii],2);
      weight_Dalitz[ii-NoffsetDalitz]=1./pow(eyieldPi0DalitzTot[ii-NoffsetDalitz],2);
      weight_EMCAL[ii-NoffsetEMCAL]=1./pow(eyieldPi0EMCALTot[ii-NoffsetEMCAL],2);
      weight_PHOS[ii-NoffsetPHOS]=1./pow(eyieldPi0PHOSTot[ii-NoffsetPHOS],2);
      totWeightPi0Comb[ii]=weight_PCM[ii]+weight_Dalitz[ii-NoffsetDalitz]+weight_EMCAL[ii-NoffsetEMCAL]+weight_PHOS[ii-NoffsetPHOS];


      yieldPi0Comb[ii]=(yieldPi0PCM[ii]*weight_PCM[ii]+
			yieldPi0Dalitz[ii-NoffsetDalitz]*weight_Dalitz[ii-NoffsetDalitz]+
			yieldPi0EMCAL[ii-NoffsetEMCAL]*weight_EMCAL[ii-NoffsetEMCAL]+
	                yieldPi0PHOS[ii-NoffsetPHOS]*weight_PHOS[ii-NoffsetPHOS])/totWeightPi0Comb[ii];
      
      
      totErrPi0Comb[ii]= pow(totWeightPi0Comb[ii],-0.5);
      
        
      staErrPi0Comb[ii]=pow((1./4.*(weight_PCM[ii]/totWeightPi0Comb[ii])*pow(eyieldPi0PCMStat[ii],2)+
			     1./4.*(weight_Dalitz[ii-NoffsetDalitz]/totWeightPi0Comb[ii])*
			     pow(eyieldPi0DalitzStat[ii-NoffsetDalitz],2)+
			     1./4.*(weight_EMCAL[ii-NoffsetEMCAL]/totWeightPi0Comb[ii])*
			     pow(eyieldPi0EMCALStat[ii-NoffsetEMCAL],2)+
			     1./4.*(weight_EMCAL[ii-NoffsetPHOS]/totWeightPi0Comb[ii])*
			     pow(eyieldPi0PHOSStat[ii-NoffsetPHOS],2)),0.5);
      

      
      sysErrPi0Comb[ii]=pow((1./4.*(weight_PCM[ii]/totWeightPi0Comb[ii])*pow(eyieldPi0PCMSyst[ii],2)+
			     1./4.*(weight_Dalitz[ii-NoffsetDalitz]/totWeightPi0Comb[ii])*
			     pow(eyieldPi0DalitzSyst[ii-NoffsetDalitz],2)+
			     1./4.*(weight_EMCAL[ii-NoffsetEMCAL]/totWeightPi0Comb[ii])*
			     pow(eyieldPi0EMCALSyst[ii-NoffsetEMCAL],2)+
			     1./4.*(weight_PHOS[ii-NoffsetPHOS]/totWeightPi0Comb[ii])*
			     pow(eyieldPi0PHOSSyst[ii-NoffsetPHOS],2)),0.5);
      
  

    
      
    }else if(ii>=17 && ii<28){          // PCM+EMCAL+PHOS
      //      cout<<"Test 4th interval PCM+EMCAL+PHOS " << xValComb[ii]<<" "<< ptPi0PCM[ii]<< "  "<<" "<< ptPi0EMCAL[ii-NoffsetEMCAL]<<" "<< ptPi0PHOS[ii-NoffsetPHOS]<< endl;
      
      weight_PCM[ii]=1./pow(eyieldPi0PCMTot[ii],2);
      weight_EMCAL[ii-NoffsetEMCAL]=1./pow(eyieldPi0EMCALTot[ii-NoffsetEMCAL],2);
      weight_PHOS[ii-NoffsetPHOS]=1./pow(eyieldPi0PHOSTot[ii-NoffsetPHOS],2);
      totWeightPi0Comb[ii]=weight_PCM[ii]+weight_EMCAL[ii-NoffsetEMCAL]+weight_PHOS[ii-NoffsetPHOS];

      yieldPi0Comb[ii]=(yieldPi0PCM[ii]*weight_PCM[ii]+
			yieldPi0EMCAL[ii-NoffsetEMCAL]*weight_EMCAL[ii-NoffsetEMCAL]+
			yieldPi0PHOS[ii-NoffsetPHOS]*weight_PHOS[ii-NoffsetPHOS])/
	totWeightPi0Comb[ii];
      totErrPi0Comb[ii]= pow(totWeightPi0Comb[ii] ,-0.5);
      
     
      staErrPi0Comb[ii]=pow((1./3.*(weight_PCM[ii]/totWeightPi0Comb[ii])*pow(eyieldPi0PCMStat[ii],2)+
			     1./3.*(weight_EMCAL[ii-NoffsetEMCAL]/totWeightPi0Comb[ii])*
			     pow(eyieldPi0EMCALStat[ii-NoffsetEMCAL],2)+
			     1./3.*(weight_PHOS[ii-NoffsetPHOS]/totWeightPi0Comb[ii])*
			     pow(eyieldPi0PHOSStat[ii-NoffsetPHOS],2)),0.5);
      

      
      sysErrPi0Comb[ii]=pow((1./3.*(weight_PCM[ii]/totWeightPi0Comb[ii])*pow(eyieldPi0PCMSyst[ii],2)+
			     1./3.*(weight_EMCAL[ii-NoffsetEMCAL]/totWeightPi0Comb[ii])*
			     pow(eyieldPi0EMCALSyst[ii-NoffsetEMCAL],2)+
			     1./3.*(weight_PHOS[ii-NoffsetPHOS]/totWeightPi0Comb[ii])*
			     pow(eyieldPi0PHOSSyst[ii-NoffsetPHOS],2)),0.5);
      
  


    
    }else if (ii>=28){                  // PHOS
      
      yieldPi0Comb[ii]= yieldPi0PHOS[ii-NoffsetPHOS];
      totErrPi0Comb[ii]=eyieldPi0PHOSTot[ii-NoffsetPHOS];
      staErrPi0Comb[ii]=eyieldPi0PHOSStat[ii-NoffsetPHOS];
      sysErrPi0Comb[ii]=eyieldPi0PHOSSyst[ii-NoffsetPHOS];


    }
  }

  TGraphAsymmErrors * graphInvYieldPi0CombpPb5023GeVTot = new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0Comb,exValComb,exValComb,totErrPi0Comb,totErrPi0Comb);
  TGraphAsymmErrors * graphInvYieldPi0CombpPb5023GeVSta = new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0Comb,exValComb,exValComb,staErrPi0Comb,staErrPi0Comb);
  TGraphAsymmErrors * graphInvYieldPi0CombpPb5023GeVSys = new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0Comb,exValComb,exValComb,sysErrPi0Comb,sysErrPi0Comb);


//   TGraphAsymmErrors* graphInvYieldPi0CombpPb5023GeVXShiftedSta  = (TGraphAsymmErrors*) graphInvYieldPi0CombpPb5023GeVSta->Clone();	
//   TGraphAsymmErrors* graphInvYieldPi0CombpPb5023GeVXShiftedSys  = (TGraphAsymmErrors*) graphInvYieldPi0CombpPb5023GeVSys->Clone();
//   TGraphAsymmErrors* graphInvYieldPi0CombpPb5023GeVXShiftedTot = (TGraphAsymmErrors*) graphInvYieldPi0CombpPb5023GeVTot->Clone();
	

  TGraphAsymmErrors* graphInvYieldPi0CombpPb5023GeVStaClone  = (TGraphAsymmErrors*) graphInvYieldPi0CombpPb5023GeVSta->Clone();	
  TGraphAsymmErrors* graphInvYieldPi0CombpPb5023GeVSysClone  = (TGraphAsymmErrors*) graphInvYieldPi0CombpPb5023GeVSys->Clone();
  TGraphAsymmErrors* graphInvYieldPi0CombpPb5023GeVTotClone = (TGraphAsymmErrors*) graphInvYieldPi0CombpPb5023GeVTot->Clone();
	
  TGraphAsymmErrors* graphRatioCombCombFit = (TGraphAsymmErrors*) graphInvYieldPi0CombpPb5023GeVTotClone ->Clone(); 
  TGraphAsymmErrors* graphRatioCombCombFitSta = (TGraphAsymmErrors*) graphInvYieldPi0CombpPb5023GeVStaClone ->Clone();
  TGraphAsymmErrors* graphRatioCombCombFitSys = (TGraphAsymmErrors*) graphInvYieldPi0CombpPb5023GeVSysClone ->Clone(); 



 cout<<"combined total errors"<<endl;
  graphInvYieldPi0CombpPb5023GeVTot->Print();
  cout<<endl;

  cout<<"combined stat errors"<<endl;
  graphInvYieldPi0CombpPb5023GeVSta->Print();

  cout<<"combined sys errors"<<endl;
  graphInvYieldPi0CombpPb5023GeVSys->Print();

  //  new TCanvas;
  
  TCanvas* canvasCompYieldpPbInd = new TCanvas("canvasCompYieldpPbInd","",200,10,700,500);  // gives the page size
  DrawGammaCanvasSettings( canvasCompYieldpPbInd,  0.12, 0.02, 0.02, 0.12);
  
  canvasCompYieldpPbInd->SetLogx();
  canvasCompYieldpPbInd->SetLogy();
  TH2F * histo2DCompCombined;
  histo2DCompCombined = new TH2F("histo2DCompCombined","histo2DCompCombined",1000,0.3,30.,1000,1.2e-9,30.   );
  SetStyleHistoTH2ForGraphs(histo2DCompCombined, "p_{T} (GeV/c)","Invariant yield", 0.05,0.064, 0.05,0.06, 0.8,0.9, 512, 505);
  histo2DCompCombined->GetXaxis()->SetRangeUser(0.,30.);
  //  histo2DCompCombined->GetYaxis()->SetRangeUser(0.1,2.1);
  histo2DCompCombined->DrawCopy();
  graphInvYieldPi0CombpPb5023GeVTot->SetLineColor(2);
  graphInvYieldPi0CombpPb5023GeVTot->SetMarkerColor(2);
  graphInvYieldPi0CombpPb5023GeVTot->Draw("E1psame");
  graphInvYieldPi0CombpPb5023GeVSta->Draw("E1psame");
  graphInvYieldPi0CombpPb5023GeVSys->Draw("E1psame");

  
  histoEMCalYieldPi0pPbSys->SetLineColor(3);
  histoEMCalYieldPi0pPbSys->SetMarkerColor(3);
  histoEMCalYieldPi0pPbSys->Draw("p,same");


  histoPHOSYieldPi0pPbSys->SetLineColor(5);
  histoPHOSYieldPi0pPbSys->SetMarkerColor(5);
  histoPHOSYieldPi0pPbSys->Draw("p,same");

  
  histoPCMYieldPi0pPb->SetLineColor(4);
  histoPCMYieldPi0pPb->SetMarkerColor(4);
  histoPCMYieldPi0pPb->Draw("p,same");

  
  histoDalitzYieldPi0pPb->SetLineColor(6);
  histoDalitzYieldPi0pPb->SetMarkerColor(6);
  histoDalitzYieldPi0pPb->Draw("p,same");
  canvasCompYieldpPbInd->Print("Comb_pPb.pdf");





  //Xshift and Fitting the combined spectrum





  Double_t paramPi0pPb5023GeV[10];
  
  ReturnParameterSetFittingPbPb("8000",paramPi0pPb5023GeV);
	
	
  TF1* fitTsallisCombPi0pPb5023GeVPt = 
    FitObject("rad","fitTsallisCombPi0pPb5023GeVPt","Pi0",
	      graphInvYieldPi0CombpPb5023GeVTotClone,0.3,24.,paramPi0pPb5023GeV,"QNRME+");
	

  fitTsallisCombPi0pPb5023GeVPt->Draw("same");





  TCanvas* canvasRatioCompYieldpPbInd = new TCanvas("canvasRatioCompYieldpPbInd","",200,10,900,700);  // gives the page size
 DrawGammaCanvasSettings( canvasRatioCompYieldpPbInd,  0.12, 0.02, 0.02, 0.12);
 canvasRatioCompYieldpPbInd->SetLogx();
 canvasRatioCompYieldpPbInd->SetGridx();
 canvasRatioCompYieldpPbInd->SetGridy();


  TH2F * ratio2DInvXSectionPi0;
  ratio2DInvXSectionPi0 = new TH2F("ratio2DInvXSectionPi0","ratio2DInvXSectionPi0",1000,0.1,30.,1000,0.3,2.);
  

  ratio2DInvXSectionPi0->SetXTitle("p_{T}(GeV/c)");
  ratio2DInvXSectionPi0->SetYTitle("Ratio to Combined");

  ratio2DInvXSectionPi0->DrawCopy(); 
  graphRatioCombPCMSys= CalculateGraphErrRatioToFit(graphRatioCombPCMSys,fitTsallisCombPi0pPb5023GeVPt);
  graphRatioCombDalitzSys= CalculateGraphErrRatioToFit(graphRatioCombDalitzSys,fitTsallisCombPi0pPb5023GeVPt);

  graphRatioCombCombFit = CalculateGraphErrRatioToFit(graphInvYieldPi0CombpPb5023GeVTotClone,fitTsallisCombPi0pPb5023GeVPt); 
  graphRatioCombCombFitSta = CalculateGraphErrRatioToFit(graphInvYieldPi0CombpPb5023GeVStaClone,fitTsallisCombPi0pPb5023GeVPt); 
  graphRatioCombCombFitSys = CalculateGraphErrRatioToFit(graphInvYieldPi0CombpPb5023GeVSysClone,fitTsallisCombPi0pPb5023GeVPt); 

  histoRatioCombPHOS = CalculateHistoRatioToFit(histoRatioCombPHOS,fitTsallisCombPi0pPb5023GeVPt);
  histoRatioCombEMCAL = CalculateHistoRatioToFit(histoRatioCombEMCAL,fitTsallisCombPi0pPb5023GeVPt);

  
  //DrawGammaSetMarkerTGraphAsym(graphRatioCombComFit, 20,2, 2, 2, 1, kTRUE);
 
  
 graphRatioCombCombFit->SetMarkerStyle(20);
  graphRatioCombCombFit->SetMarkerColor(kBlack);
  graphRatioCombCombFit->SetMarkerSize(.75);
  
  graphRatioCombCombFit->Draw("pe,same");

  

  graphRatioCombPCMSys->SetMarkerColor(2);
  graphRatioCombPCMSys->SetMarkerStyle(21);
  graphRatioCombPCMSys->SetLineColor(2);
  graphRatioCombPCMSys->SetMarkerSize(0.5);
  graphRatioCombPCMSys->Draw("pE,same");
 cout<<"ratio Combined pcm"<<endl;
 graphRatioCombPCMSys->Print();

   graphRatioCombDalitzSys->SetMarkerColor(4);

   graphRatioCombDalitzSys->SetLineColor(4);
  graphRatioCombDalitzSys->SetMarkerStyle(22);
   graphRatioCombDalitzSys->SetMarkerSize(0.5);
   // graphRatioCombDalitzSys->SetAxisRange(0.3,4.);
  graphRatioCombDalitzSys->Draw("pE,same");


  histoRatioCombPHOS->SetMarkerColor(3);
  histoRatioCombEMCAL->SetMarkerColor(6);
  histoRatioCombPHOS->SetLineColor(3);
  histoRatioCombEMCAL->SetLineColor(6);
  histoRatioCombEMCAL->SetMarkerStyle(23);

  histoRatioCombPHOS->Draw("same");
  histoRatioCombEMCAL->Draw("same");

  cout<< "  ratio combined dalitz"<<endl;
  graphRatioCombDalitzSys->Print();





  /*
 TH1D* histoRatioChargedPionSpecSys= (TH1D*)histoChargedPionSpecSys->Clone();

 histoRatioChargedPionSpecSys = CalculateHistoRatioToFit(histoChargedPionSpecSys,fitTsallisCombPi0pPb5023GeVPt);
 histoRatioChargedPionSpecSys->Draw("same");
 histoRatioChargedPionSpecSys->SetLineColor(7);
 histoRatioChargedPionSpecSys->SetMarkerColor(7);
  */


 TLegend * leg= new TLegend(0.15,0.15,0.4,0.3);

 leg->AddEntry(graphRatioCombCombFit,"Combined spectrum to fit of Combined, total err","l");
 leg->AddEntry(graphRatioCombPCMSys,"PCM to fit of Combined, sys err.","l");
 leg->AddEntry(graphRatioCombDalitzSys,"Dalitz to fit of Combined, sys err.","l");
 leg->AddEntry(histoRatioCombPHOS,"PHOS to fit of Combined, sys err.","l");
 leg->AddEntry(histoRatioCombEMCAL,"EMCAL to fit of Combined, sys err.","l");
 // leg->AddEntry(histoRatioChargedPionSpecSys,"Charged pions to fit of Combined, sys err.","l");
 leg->Draw();


  canvasRatioCompYieldpPbInd->Print("RatioCompYieldpPb.pdf");
}
