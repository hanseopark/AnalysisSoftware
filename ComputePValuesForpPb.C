/****************************************************************************************************************************
****** 		provided by Gamma Conversion Group, PWG4, 													*****
******		Ana Marin, marin@physi.uni-heidelberg.de													*****
******	   	Kathrin Koch, kkoch@physi.uni-heidelberg.de 													*****
******		Friederike Bock, friederike.bock@cern.ch													*****
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
#include "CommonHeaders/Interpolation.h"
#include "TPaveText.h"

extern TRandom*	gRandom;
extern TBenchmark*	gBenchmark;
extern TSystem*	gSystem;
extern TMinuit*  	gMinuit;




TGraphAsymmErrors* GetSpectrumWithSystErrors(TH1D* spectrum, TString fileNameSysErr,TString outputName,Int_t offsetSyst);
TGraphAsymmErrors* CalculateRatio(TGraphAsymmErrors* ObjA_StatErr, TGraphAsymmErrors* ObjA_SystErr, TGraphAsymmErrors* ObjB_StatErr, TGraphAsymmErrors* ObjB_SystErr,Int_t startBin, Int_t endBin);


	
	
TGraphAsymmErrors* g_methodA_stat[6];
TGraphAsymmErrors* g_methodA_sysA[6];
TGraphAsymmErrors* g_methodA_sysB[6];
TGraphAsymmErrors* g_methodA_sysC[6];
TGraphAsymmErrors* g_methodB_stat[6];
TGraphAsymmErrors* g_methodB_sysA[6];
TGraphAsymmErrors* g_methodB_sysB[6];
TGraphAsymmErrors* g_methodB_sysC[6];

enum { kDalitz = 0, kPCM, kPHOS, kEMCAL };
TString methodName[4]={"Dalitz","PCM","PHOS","EMCAL"};
Int_t   ratioAB[6][2];
Int_t   ratioABRange[6][2];
//Double_t pValue[6][2];
const Int_t  Ntotal = 34;
const Int_t  nPtLimits = Ntotal+1;
Double_t pValue[6][nPtLimits+1][2];
TH1F* htest_statistic[6][nPtLimits+1];
Double_t xValComb[Ntotal] = {0};
Double_t xErr[Ntotal]={0};


TRandom rndm;
	



Double_t p_value_to_n_sigma(const Double_t& pvalue) {

  TF1 f("f","1-TMath::Erf(x/sqrt(2.)) - [0]", -100, 100);
  f.SetParameter(0, pvalue);

  // create wrapper function
  ROOT::Math::WrappedTF1 wfrf(f);

  // create root finder
  ROOT::Math::BrentRootFinder brf;

  // set parameters of the method
  brf.SetFunction(wfrf, 0., 50.);
  brf.Solve();
  Double_t nsigma = brf.Root();

  return nsigma;
}

  
Double_t test_statistic(const bool real_data = false, Int_t iRatio=0,Int_t startBin=3,Int_t endBin=14) {
  
  // null hypothesis for ratio methodB/methodA
  Double_t R0 = 1.;

  Double_t test_statistic = 0;

  Double_t eps_b_methodA = 0;
  Double_t eps_c_methodA = 0;
  Double_t eps_b_methodB = 0;
  Double_t eps_c_methodB = 0;

  if (!real_data) {
    eps_b_methodA = rndm.Gaus(0, 1);
    eps_c_methodA = rndm.Gaus(0, 1);
    eps_b_methodB = rndm.Gaus(0, 1);
    eps_c_methodB = rndm.Gaus(0, 1);
  }
  
  
  
    for (Int_t i=startBin;i<=endBin; i++) {

      Double_t Y_methodA = g_methodA_stat[iRatio]->GetY()[i];
      Double_t Y_methodA_ErrStat = g_methodA_stat[iRatio]->GetEYhigh()[i];
      Double_t Y_methodA_ErrSysA = g_methodA_sysA[iRatio]->GetEYhigh()[i];
      Double_t Y_methodA_ErrStatSysA = sqrt(Y_methodA_ErrStat*Y_methodA_ErrStat + Y_methodA_ErrSysA*Y_methodA_ErrSysA);
      Double_t Y_methodA_ErrStatSysA_rel = Y_methodA_ErrStatSysA / Y_methodA;
      Double_t Y_methodA_ErrSysB_rel = g_methodA_sysB[iRatio]->GetEYhigh()[i] / Y_methodA;
      Double_t Y_methodA_ErrSysC_rel = g_methodA_sysC[iRatio]->GetEYhigh()[i] / Y_methodA;
      

      Double_t Y_methodB = g_methodB_stat[iRatio]->GetY()[i];
      Double_t Y_methodB_ErrStat = g_methodB_stat[iRatio]->GetEYhigh()[i];
      Double_t Y_methodB_ErrSysA = g_methodB_sysA[iRatio]->GetEYhigh()[i];
      Double_t Y_methodB_ErrStatSysA = sqrt(Y_methodB_ErrStat*Y_methodB_ErrStat + Y_methodB_ErrSysA*Y_methodB_ErrSysA);
      Double_t Y_methodB_ErrStatSysA_rel = Y_methodB_ErrStatSysA / Y_methodB;
      Double_t Y_methodB_ErrSysB_rel = g_methodB_sysB[iRatio]->GetEYhigh()[i] / Y_methodB;
      Double_t Y_methodB_ErrSysC_rel = g_methodB_sysC[iRatio]->GetEYhigh()[i] / Y_methodB;

      Double_t sigma_rel = sqrt(Y_methodA_ErrStatSysA_rel*Y_methodA_ErrStatSysA_rel + Y_methodB_ErrStatSysA_rel*Y_methodB_ErrStatSysA_rel);

      Double_t R = 0; // methodB/methodA for real data or pseudo data
      if (!real_data) {

	// generate pseudo data for methodB/methodA ratio
	Double_t f_smear_methodA = (1. + Y_methodA_ErrSysB_rel * eps_b_methodA) * (1. + Y_methodA_ErrSysC_rel * eps_c_methodA);
	Double_t f_smear_methodB = (1. + Y_methodB_ErrSysB_rel * eps_b_methodB) * (1. + Y_methodB_ErrSysC_rel * eps_c_methodB);
	Double_t R_mod = R0 * f_smear_methodB/f_smear_methodA;

	// cout << "f_smear_methodA = " << f_smear_methodA << endl;
	// cout << "f_smear_methodB = " << f_smear_methodB << endl;
	// cout << "Y_methodB_ErrSysB_rel = " << Y_methodB_ErrSysB_rel << endl;
	// cout << "Y_methodB_ErrSysC_rel = " << Y_methodB_ErrSysC_rel << endl;
	// cout << "sigma_rel = " << sigma_rel << endl;
	// cout << "R_mod = " << R_mod << endl;
	
	R = rndm.Gaus(R_mod, R_mod * sigma_rel);

	// cout << "R = " << R << endl;
      }
      else {
	R = Y_methodB / Y_methodA;
      }
      
      test_statistic += pow((R - R0)/(R0 * sigma_rel), 2);
		
    }
 
  return test_statistic;

}



void ComputePValue(Bool_t binBybin = kFALSE){
  
  const Int_t nevt = 100000;
  
  // chi2 histogram
  //TH1F htest_statistic("htest_statistic","htest_statistic", 3000, 0., 3000.);
  //htest_statistic.SetXTitle("test statistic t");
  //htest_statistic.SetYTitle("counts");
  
  Int_t startBin = nPtLimits;
  
  if(binBybin) startBin = 0;
  
  
  
  for(Int_t ii = 0; ii < 6; ii++){
    
     TString detA = methodName[ratioAB[ii][0]];
     TString detB = methodName[ratioAB[ii][1]];
     
    /*
    if( binBybin == kFALSE ){
      
      htest_statistic[ii][nPtLimits] = new TH1F(Form("%s_%s_statistic_%d",detA.Data(),detB.Data(),ii),Form("%s_%s_statistic_%d",detA.Data(),detB.Data(),ii), 3000, 0., 3000.);
      htest_statistic[ii][nPtLimits]->SetXTitle("test statistic t");
      htest_statistic[ii][nPtLimits]->SetYTitle("counts");
    
    } else {*/
      
        for(Int_t iBin = startBin; iBin < nPtLimits+1; iBin++){
	  
	  htest_statistic[ii][iBin] = new TH1F(Form("%s_%s_statistic_%d",detA.Data(),detB.Data(),iBin),Form("%s_%s_statistic_%d",detA.Data(),detB.Data(),iBin), 3000, 0., 3000.);
	  htest_statistic[ii][iBin]->SetXTitle("test statistic t");
	  htest_statistic[ii][iBin]->SetYTitle("counts");
	  
	}      
    //}
    
    cout<<"Llego aqui"<<endl;
    
    
    bool real_data_flag;
    
    //cout<<"Method A"<<endl;
    //g_methodA_stat[ii]->Print();
    //g_methodA_sysA[ii]->Print();
    //g_methodA_sysB[ii]->Print();
    //g_methodA_sysC[ii]->Print();
    
    
    
    //cout<<"Computing the p-Value for the ratio: "<<methodName[ratioAB[ii][0]].Data()<<"/"<<methodName[ratioAB[ii][1]].Data()<<endl;
    
        
	  for (Int_t ievt=0; ievt<nevt; ievt++) {

		real_data_flag = false;
		
		htest_statistic[ii][nPtLimits]->Fill(test_statistic(real_data_flag,ii,ratioABRange[ii][0],ratioABRange[ii][1]));
		
		if( binBybin == kTRUE ){
		  
		   for(Int_t iBin = ratioABRange[ii][0]; iBin <= ratioABRange[ii][1]; iBin++){
		     
		     htest_statistic[ii][iBin]->Fill(test_statistic(real_data_flag,ii,iBin,iBin));
		     
		   }
		  
		}
		
	  }
	  
	  cout<<"Paso aqui"<<endl;
	  
	  Double_t test_statistic_data[nPtLimits+1];
	 
	  
	   real_data_flag = true;
	   test_statistic_data[nPtLimits] = test_statistic(real_data_flag,ii,ratioABRange[ii][0],ratioABRange[ii][1]);
	  
	    if( binBybin == kTRUE ){
	      
	       for(Int_t iBin = ratioABRange[ii][0]; iBin <= ratioABRange[ii][1]; iBin++){
		     
		     test_statistic_data[iBin] = test_statistic(real_data_flag,ii,iBin,iBin);
		     
		}
	      
	    }
	  
	   
	  
	  //cout << "test_statistic_data = " << test_statistic_data << endl;
	  
	  Double_t int_total = htest_statistic[ii][nPtLimits]->GetEntries();
	  Double_t ib_test_statistic_data = htest_statistic[ii][nPtLimits]->FindBin(test_statistic_data[nPtLimits]);
	  Double_t pval = 1 - htest_statistic[ii][nPtLimits]->Integral(1, ib_test_statistic_data) / int_total;

	  //TString s = "p-value (comparison Dalitz/PCM inlusive photon spectrum) = ";
	  //if (takeDR) s = "p-value (comparison PCM/PHOS double ratio) = ";
  
	  pValue[ii][nPtLimits][0] = pval;
	  pValue[ii][nPtLimits][1] = p_value_to_n_sigma(pval);
	  
	  if( binBybin ){
	    
	      for(Int_t iBin = ratioABRange[ii][0]; iBin <= ratioABRange[ii][1]; iBin++){
	    
	      int_total = htest_statistic[ii][iBin]->GetEntries();
	      ib_test_statistic_data = htest_statistic[ii][iBin]->FindBin(test_statistic_data[iBin]);
	      pval = 1 - htest_statistic[ii][iBin]->Integral(1, ib_test_statistic_data) / int_total;

	  //TString s = "p-value (comparison Dalitz/PCM inlusive photon spectrum) = ";
	  //if (takeDR) s = "p-value (comparison PCM/PHOS double ratio) = ";
  
	      pValue[ii][iBin][0] = pval;
	      pValue[ii][iBin][1] = p_value_to_n_sigma(pval);
	      
	      //cout <<"("<<methodName[ratioAB[ii][0]].Data()<<"/"<<methodName[ratioAB[ii][1]].Data()<<"), pT range: ("<< xValComb[iBin]<<" - "<< xValComb[iBin]<<" GeV)  pVal "<< pValue[ii][iBin][0]  << ", nsigma = " <<  pValue[ii][iBin][1] << endl;
	 
	      
	      }
	    
	  }
	  
	  
	  
	  
	  //cout <<"("<<methodName[ratioAB[ii][0]].Data()<<"/"<<methodName[ratioAB[ii][1]].Data()<<"), pT range: ("<< xValComb[ratioABRange[ii][0]]<<" - "<< xValComb[ratioABRange[ii][1]]<<" GeV)  pVal "<<  pValue[ii][nPtLimits][0] << ", nsigma = " << pValue[ii][nPtLimits][1] << endl;
	  
	  
    
  }
    
    
    
 }

  
  



void ComputePValuesForpPb(TString suffix="pdf",Bool_t binBybin = kFALSE){

	gROOT->Reset();	
	gROOT->SetStyle("Plain");
	
	
	
	StyleSettingsThesis();	
	SetPlotStyle();
	
	
	TString dateForOutput 				= ReturnDateStringForOutput();
	TString outputDir 				= Form("PValuesForpPb/%s/%s",suffix.Data(),dateForOutput.Data());
	gSystem->Exec("mkdir -p "+outputDir);
	
	
	
	//TString fileNameNeutralPionDalitz       	= "ExternalInputpPb/PCM/data_PCMResults_Dalitz_pPb_2015-06-28.root";
	TString fileNameNeutralPionDalitz		= "ExternalInputpPb/PCM/data_PCMResults_Dalitz_pPb_20150806.root";
	
	//TString fileNameNeutralPionPCM          	= "ExternalInputpPb/PCM/data_PCMResults_pPb_20150623_standard_dc4.root";
	TString fileNameNeutralPionPCM                  = "ExternalInputpPb/PCM/data_PCMResults_pPb_20150804_standard_CatErrors.root";
	
	
	TString	fileNameNeutralPionPHOS 		= "ExternalInputpPb/PHOS/data_PHOSResults_pPb_20150623.root";
	TString nameHistoPHOS 				= "hCor_stat";
	TString nameHistoPHOSSysErrors 			= "hCor_syst";
	
	
	
	TString fileNameNeutralPionEMCAL 		= "ExternalInputpPb/EMCAL/EMCALResults_11Sept.root"; 
	TString nameHistoEMCal 				= "CorrectedYieldPi0";
	TString nameHistoEMCalSysErrors 		= "Pi0SystError";
	
	
	TString collisionSystempPb = "p-Pb, MB #sqrt{s_{NN}} = 5.02 TeV";     
	
	cout<<"kDalitz: "<<kDalitz<<endl;
	
	
	ratioAB[0][0]=kDalitz; ratioAB[0][1]=kPCM;//Dalitz/PCM
	ratioAB[1][0]=kDalitz; ratioAB[1][1]=kPHOS; //Dalitz/PHOS
	ratioAB[2][0]=kDalitz; ratioAB[2][1]=kEMCAL; //Dalitz/EMCAL
	ratioAB[3][0]=kPCM;    ratioAB[3][1]=kPHOS; //PCM/PHOS
	ratioAB[4][0]=kPCM;    ratioAB[4][1]=kEMCAL; //PCM/EMCAL
	ratioAB[5][0]=kPHOS;   ratioAB[5][1]=kEMCAL; //PHOS/EMCAL
	
	
	ratioABRange[0][0]=3; ratioABRange[0][1]=16;//Dalitz/PCM
	ratioABRange[1][0]=6; ratioABRange[1][1]=16; //Dalitz/PHOS
	ratioABRange[2][0]=6; ratioABRange[2][1]=16; //Dalitz/EMCAL
	ratioABRange[3][0]=6; ratioABRange[3][1]=29; //PCM/PHOS
	ratioABRange[4][0]=6; ratioABRange[4][1]=27; //PCM/EMCAL
	ratioABRange[5][0]=6; ratioABRange[5][1]=27; //PHOS/EMCAL
	
	
	
	Double_t materialErr = 4.5;
	Double_t dalitzBRErr = 2.98;
	
	
	
	
	Double_t pTLimits[nPtLimits] =			       {0.3, 0.4, 0.5, 0.6, 0.7, 
								0.8, 1.0, 1.2, 1.4, 1.6, 
								1.8, 2.0, 2.2, 2.4, 2.6,
								2.8, 3.0, 3.2, 3.4, 3.6,
								3.8, 4.0, 4.5, 5.0, 5.5,
								6.0, 7.0, 8.0, 10.0,12.0,
								14.0, 16.0, 18.0, 20.0, 25.0};
	
	
	
								
								
	for(Int_t ii=0; ii < nPtLimits-1; ii++){
	  
	  //cout<<(pTLimits[ii]+((pTLimits[ii+1] - pTLimits[ii])/2))<<" ";
	  xValComb[ii] = pTLimits[ii] + (pTLimits[ii+1]-pTLimits[ii])/2;
	  xErr[ii]     = (pTLimits[ii+1] -pTLimits[ii] ) / 2;
	  //cout<<xValComb[ii]<<" +/- "<<xErr[ii]<<" ";
	  
	  //if(((ii+1)%5)==0) cout<<"\n";
	  
	}
	
	
	
	const Int_t nP_Dalitz = 14;   
	  
	Int_t nDalitzOffset=3 ;
	  
	Double_t ptPi0Dalitz[Ntotal]={0};
	Double_t xPi0DalitzErr[Ntotal]={0};
	Double_t yieldPi0Dalitz[Ntotal]={0};
	Double_t yPi0DalitzStatErr[Ntotal]={0};
	Double_t yPi0DalitzSystErr[Ntotal]={0};
	Double_t yPi0DalitzSystErrMatCanceled[Ntotal]={0};
	Double_t yPi0DalitzTotErr[Ntotal]={0};
	Double_t yPi0DalitzSystErrA[Ntotal]={0};
	Double_t yPi0DalitzSystErrB[Ntotal]={0};
	Double_t yPi0DalitzSystErrC[Ntotal]={0};
	Double_t yPi0DalitzSystErrCMatCanceled[Ntotal]={0};
	
	
	const Int_t nP_PCM = 30;
	Int_t nPCMOffset = 0;
	  
	Double_t ptPi0PCM[Ntotal]={0};
	Double_t xPi0PCMErr[Ntotal]={0};
	Double_t yieldPi0PCM[Ntotal]={0};
	Double_t yPi0PCMStatErr[Ntotal]={0};
	Double_t yPi0PCMSystErr[Ntotal]={0};
	Double_t yPi0PCMSystErrMatCanceled[Ntotal]={0};
	Double_t yPi0PCMTotErr[Ntotal]={0};
	Double_t yPi0PCMSystErrA[Ntotal]={0};
	Double_t yPi0PCMSystErrB[Ntotal]={0};
	Double_t yPi0PCMSystErrC[Ntotal]={0};
	Double_t yPi0PCMSystErrCMatCanceled[Ntotal]={0};
	
	const Int_t nP_PHOS = 27;
	Int_t nPHOSOffset = 6;
	
	Double_t ptPi0PHOS[Ntotal]={0};
	Double_t xPi0PHOSErr[Ntotal]={0};
	Double_t yieldPi0PHOS[Ntotal]={0};
	Double_t yPi0PHOSStatErr[Ntotal]={0};
	Double_t yPi0PHOSSystErr[Ntotal]={0};
	Double_t yPi0PHOSTotErr[Ntotal]={0};
	Double_t yPi0PHOSSystErrA[Ntotal]={0};
	Double_t yPi0PHOSSystErrB[Ntotal]={0};
	Double_t yPi0PHOSSystErrC[Ntotal]={0};
	Double_t yPi0PHOSSystErrCW0Mat[Ntotal]={0};
	
	const Int_t nP_EMCAL = 22;
	Int_t nEMCALOffset = 6;
  
	
	Double_t ptPi0EMCAL[Ntotal]={0};
	Double_t xPi0EMCALErr[Ntotal]={0};
	Double_t yieldPi0EMCAL[Ntotal]={0};
	Double_t yPi0EMCALStatErr[Ntotal]={0};
	Double_t yPi0EMCALSystErr[Ntotal]={0};
	Double_t yPi0EMCALTotErr[Ntotal]={0};
	Double_t yPi0EMCALSystErrA[Ntotal]={0};
	Double_t yPi0EMCALSystErrB[Ntotal]={0};
	Double_t yPi0EMCALSystErrC[Ntotal]={0};
	Double_t yPi0EMCALSystErrCW0Mat[Ntotal]={0};
	
	
	
	
	TFile*  fileNeutralPionsDalitz = new TFile(fileNameNeutralPionDalitz.Data());
	TDirectory* directoryDalitzPi0pPb = (TDirectory*) fileNeutralPionsDalitz->GetDirectory("Pi0_pPb_5.023TeV_0-100%");
	
	if( ! directoryDalitzPi0pPb ){
	  cout<<"Dalitz: The directory Pi0_pPb_5.023TeV_0-100% does not exist "<<endl; 
	  return;
	}
	
	TFile*  fileNeutralPionPCM = new TFile(fileNameNeutralPionPCM.Data());
	TDirectory* directoryPCMPi0pPb = (TDirectory*) fileNeutralPionPCM->GetDirectory("Pi0_pPb_5.023TeV_0-100%");
	
	if( ! directoryPCMPi0pPb ){
	  
	  cout<<"PCM: The directory Pi0_pPb_5.023TeV_0-100% does not exist "<<endl; 
	  return;
	  
	}
	
	TFile*		fileNeutralPionEMCAL = 				new TFile(fileNameNeutralPionEMCAL);
	TDirectory*	directoryEMCalPi0pPb = 				(TDirectory*)fileNeutralPionEMCAL->Get("Pi05.02TeV_pPb");
	
	if( ! directoryEMCalPi0pPb ){
	  cout<<"EMCAL: The directory Pi0_pPb_5.023TeV_0-100% does not exist "<<endl; 
	  return;	  
	}
	
	
	TFile*		fileNeutralPionPHOS = 				new TFile(fileNameNeutralPionPHOS);
	//TDirectory*	directoryPHOSPi0pPb = 				(TDirectory*)fileNeutralPionPHOS->Get("Pi0_pPb_5.023TeV_0-100%");
	
	if( ! fileNeutralPionPHOS ){
	  cout<<"PHOS: The file fileNeutralPionPHOS does not exist "<<endl; 
	  return;	  
	}
	
	

	
		
	
        TH1D* histoDalitzYieldPi0pPb 			= (TH1D*)directoryDalitzPi0pPb->Get("CorrectedYieldPi0");
	TGraphAsymmErrors* tempDalitzSystErr 		= (TGraphAsymmErrors*)directoryDalitzPi0pPb->Get("Pi0SystError");
	TGraphAsymmErrors* tempDalitzSystErrWOMaterial 	= (TGraphAsymmErrors*)directoryDalitzPi0pPb->Get("Pi0SystErrorA");
	TGraphAsymmErrors* tempDalitzSystErrA		= (TGraphAsymmErrors*)directoryDalitzPi0pPb->Get("Pi0SystErrorCatA");
	TGraphAsymmErrors* tempDalitzSystErrB		= (TGraphAsymmErrors*)directoryDalitzPi0pPb->Get("Pi0SystErrorCatB");
	TGraphAsymmErrors* tempDalitzSystErrC		= (TGraphAsymmErrors*)directoryDalitzPi0pPb->Get("Pi0SystErrorCatC");
	TGraphAsymmErrors* tempDalitzSystErrCWOMaterial = (TGraphAsymmErrors*)directoryDalitzPi0pPb->Get("Pi0SystErrorCatCWOMaterial");
	TGraphAsymmErrors* tempDalitzComplErr           = (TGraphAsymmErrors*)directoryDalitzPi0pPb->Get("Pi0ComplError");
	
	TH1D* histoPCMYieldPi0pPb =            		(TH1D*)directoryPCMPi0pPb->Get("CorrectedYieldPi0");
	TGraphAsymmErrors* tempPCMSystErr   =           (TGraphAsymmErrors*)directoryPCMPi0pPb->Get("Pi0SystError");
	TGraphAsymmErrors* tempPCMSystErrWOMaterial =  (TGraphAsymmErrors*)directoryPCMPi0pPb->Get("Pi0SystErrorA");
	TGraphAsymmErrors* tempPCMSystErrA  =		(TGraphAsymmErrors*)directoryPCMPi0pPb->Get("Pi0SystErrorCatA");
	TGraphAsymmErrors* tempPCMSystErrB  =		(TGraphAsymmErrors*)directoryPCMPi0pPb->Get("Pi0SystErrorCatB");
	TGraphAsymmErrors* tempPCMSystErrC  =		(TGraphAsymmErrors*)directoryPCMPi0pPb->Get("Pi0SystErrorCatC");
	TGraphAsymmErrors* tempPCMComplErr  =           (TGraphAsymmErrors*)directoryPCMPi0pPb->Get("Pi0ComplError");
	
	
	tempPCMSystErrC->Print();
	
	
	TH1D* histoPHOSYieldPi0pPbStat  =                    (TH1D*)fileNeutralPionPHOS->Get(nameHistoPHOS.Data());
	TH1D* histoPHOSYieldPi0pPbSyst  =                    (TH1D*)fileNeutralPionPHOS->Get(nameHistoPHOSSysErrors.Data());
	TH1D* histoPHOSYieldPi0pPbSystA =		     (TH1D*)fileNeutralPionPHOS->Get("hSysTypeA");
	TH1D* histoPHOSYieldPi0pPbSystB =		     (TH1D*)fileNeutralPionPHOS->Get("hSysTypeB");
	TH1D* histoPHOSYieldPi0pPbSystC =                    (TH1D*)fileNeutralPionPHOS->Get("hSysTypeC");
	
	TH1D* histoEMCALYieldPi0pPbStat =                    (TH1D*)directoryEMCalPi0pPb->Get(nameHistoEMCal.Data());
	TH1D* histoEMCALYieldPi0pPbSyst =                    (TH1D*)directoryEMCalPi0pPb->Get(nameHistoEMCalSysErrors.Data());
	
	
	
	TGraphAsymmErrors* temp = new TGraphAsymmErrors(histoEMCALYieldPi0pPbSyst);
	

	for(Int_t ii = 0; ii < nP_Dalitz; ii++){
	      
	      ptPi0Dalitz[ii+nDalitzOffset]		= histoDalitzYieldPi0pPb->GetBinCenter(ii + 2);
	      xPi0DalitzErr[ii+nDalitzOffset]           = xErr[ii+nDalitzOffset];
	      
	      cout<<ptPi0Dalitz[ii+nDalitzOffset]<<"  ptPi0Dalitz: "<<xPi0DalitzErr[ii+nDalitzOffset]<<endl;
	      
	      yieldPi0Dalitz[ii + nDalitzOffset]	= histoDalitzYieldPi0pPb->GetBinContent(ii + 2);
	      yPi0DalitzStatErr[ii + nDalitzOffset]	= histoDalitzYieldPi0pPb->GetBinError(ii + 2);
	      yPi0DalitzSystErr[ii + nDalitzOffset]     = tempDalitzSystErr->GetErrorYhigh(ii);
	      //Double_t materialRelErr   		= materialErr*yieldPi0Dalitz[ ii + nDalitzOffset]/100;
	      //Double_t dalitzBRRelErr   		= dalitzBRErr*yieldPi0Dalitz[ ii + nDalitzOffset]/100;
	      
	
	      yPi0DalitzSystErrMatCanceled[ii+nDalitzOffset]    = tempDalitzSystErrWOMaterial->GetErrorYhigh(ii);
	      yPi0DalitzSystErrA[ii+nDalitzOffset]   		= tempDalitzSystErrA->GetErrorYhigh(ii); 
	      yPi0DalitzSystErrB[ii+nDalitzOffset]   		= tempDalitzSystErrB->GetErrorYhigh(ii); 
	      yPi0DalitzSystErrC[ii+nDalitzOffset]		= tempDalitzSystErrC->GetErrorYhigh(ii); 
	      yPi0DalitzSystErrCMatCanceled[ii+nDalitzOffset]	= tempDalitzSystErrCWOMaterial->GetErrorYhigh(ii); 
	      yPi0DalitzTotErr[ii+nDalitzOffset]		= tempDalitzComplErr->GetErrorYhigh(ii);
	      
		      
	}
	    
	    
	    
	    for(Int_t ii = 0; ii < nP_PCM; ii++){
	      
	      ptPi0PCM[ii+nPCMOffset]		=  histoPCMYieldPi0pPb->GetBinCenter(ii + 2);
	     
	      xPi0PCMErr[ii+nPCMOffset]         =  xErr[ii+nPCMOffset];
	      
	      yieldPi0PCM[ii + nPCMOffset]	=  histoPCMYieldPi0pPb->GetBinContent(ii + 2);
	      yPi0PCMStatErr[ii + nPCMOffset]	=  histoPCMYieldPi0pPb->GetBinError(ii + 2);
	      yPi0PCMSystErr[ii + nPCMOffset]   =  tempPCMSystErr->GetErrorYhigh(ii);
	      
	      Double_t materialRelErr 		=  materialErr*yieldPi0PCM[ ii + nPCMOffset]/100;
	      Double_t twomaterialRelErr  	=  (2*materialErr)*yieldPi0PCM[ ii + nPCMOffset]/100;
	      
	    
	      
	        
	      yPi0PCMSystErrMatCanceled[ii+nPCMOffset] =  pow( (pow( yPi0PCMSystErr[ ii + nPCMOffset ],2 ) - pow(twomaterialRelErr,2) + pow(materialRelErr,2) ) , 0.5);
	  
	      yPi0PCMSystErrA[ii+nPCMOffset]   		  =  tempPCMSystErrA->GetErrorYhigh(ii);//Temporal 
	      yPi0PCMSystErrB[ii+nPCMOffset]   		  =  tempPCMSystErrB->GetErrorYhigh(ii);//Temporal 
	      yPi0PCMSystErrC[ii+nPCMOffset]		  =  tempPCMSystErrC->GetErrorYhigh(ii);//Temporal  
	      
	      cout<<"Error C: "<<yPi0PCMSystErrC[ii+nPCMOffset]/yieldPi0PCM[ii + nPCMOffset] <<endl;
	      
	      yPi0PCMSystErrCMatCanceled[ii+nPCMOffset]	  =  materialRelErr;
	      yPi0PCMTotErr[ii+nPCMOffset]		= pow( ( pow(yPi0PCMStatErr[ii+nPCMOffset],2)+ pow(yPi0PCMSystErr[ii+nPCMOffset],2) ),0.5);
	      
	  	      
	    }
	    
	    
	      
	      
	    for(Int_t ii = 0; ii < nP_PHOS; ii++){
	      
	      ptPi0PHOS[ii	  +  nPHOSOffset]  =  histoPHOSYieldPi0pPbStat->GetBinCenter(ii  + 2);
	      xPi0PHOSErr[ii	  +  nPHOSOffset]  =  xErr[ii+nPHOSOffset];
	      yieldPi0PHOS[ii 	  +  nPHOSOffset]  =  histoPHOSYieldPi0pPbStat->GetBinContent(ii + 2);
	      yPi0PHOSStatErr[ii  +  nPHOSOffset]  =  histoPHOSYieldPi0pPbStat->GetBinError(ii   + 2);
	      yPi0PHOSSystErr[ii  +  nPHOSOffset]  =  histoPHOSYieldPi0pPbSyst->GetBinError(ii   + 2);
	      yPi0PHOSSystErrA[ii +  nPHOSOffset]  =  histoPHOSYieldPi0pPbSystA->GetBinContent(ii  + 2)* yieldPi0PHOS[ii 	  +  nPHOSOffset];
	      yPi0PHOSSystErrB[ii +  nPHOSOffset]  =  histoPHOSYieldPi0pPbSystB->GetBinContent(ii  + 2)* yieldPi0PHOS[ii 	  +  nPHOSOffset];
	      yPi0PHOSSystErrC[ii +  nPHOSOffset]  =  histoPHOSYieldPi0pPbSystC->GetBinContent(ii  + 2)* yieldPi0PHOS[ii 	  +  nPHOSOffset];
	      yPi0PHOSTotErr[ii+nPHOSOffset]	   = pow( ( pow(yPi0PHOSStatErr[ii+nPHOSOffset],2)+ pow(yPi0PHOSSystErr[ii+nPHOSOffset],2) ),0.5);
	      
	      Double_t temp = pow( ( pow(yPi0PHOSSystErrA[ii+nPHOSOffset],2)+ pow(yPi0PHOSSystErrB[ii+nPHOSOffset],2) +  pow(yPi0PHOSSystErrC[ii+nPHOSOffset],2) ),0.5);
	      
		      
	    }
	    
	    for(Int_t ii = 0; ii < nP_EMCAL; ii++){
	      
	      ptPi0EMCAL[ii	  +  nEMCALOffset]  =  histoEMCALYieldPi0pPbStat->GetBinCenter(ii  + 8);
	      xPi0EMCALErr[ii	  +  nEMCALOffset]  =  xErr[ii+nEMCALOffset];
	      yieldPi0EMCAL[ii 	  +  nEMCALOffset]  =  histoEMCALYieldPi0pPbStat->GetBinContent(ii + 8);
	      yPi0EMCALStatErr[ii  +  nEMCALOffset]  =  histoEMCALYieldPi0pPbStat->GetBinError(ii   + 8);
	      yPi0EMCALSystErr[ii  +  nEMCALOffset]  =  histoEMCALYieldPi0pPbSyst->GetBinError(ii   + 8);
	      
	      Double_t energyScale = 0.03*yieldPi0EMCAL[ii + nEMCALOffset];
	      
	      yPi0EMCALSystErrA[ii +  nEMCALOffset]  = pow ( (pow( yPi0EMCALSystErr[ii + nEMCALOffset],2) - pow(energyScale,2) ),0.5);
	      yPi0EMCALSystErrB[ii +  nEMCALOffset]  =  0;
	      yPi0EMCALSystErrC[ii +  nEMCALOffset]  =  energyScale;
	      yPi0EMCALTotErr[ii+nEMCALOffset]	   = pow( ( pow(yPi0EMCALStatErr[ii+nEMCALOffset],2)+ pow(yPi0EMCALSystErr[ii+nEMCALOffset],2) ),0.5);
	      
	      Double_t temp = pow( ( pow(yPi0EMCALSystErrA[ii+nEMCALOffset],2)+ pow(yPi0EMCALSystErrB[ii+nEMCALOffset],2) +  pow(yPi0EMCALSystErrC[ii+nEMCALOffset],2) ),0.5);
	      
	      
	    }
	    
	
	 
	    TGraphAsymmErrors* graphInvYieldPi0MethodsStatErr[4];
	    TGraphAsymmErrors* graphInvYieldPi0MethodsSystErrA[4];
	    TGraphAsymmErrors* graphInvYieldPi0MethodsSystErrB[4];
	    TGraphAsymmErrors* graphInvYieldPi0MethodsSystErrC[4];
	    
	   
	    TGraphAsymmErrors* graphInvYieldPi0DalitzSystErr  			= new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0Dalitz,xPi0DalitzErr,xPi0DalitzErr,yPi0DalitzSystErr,  			yPi0DalitzSystErr);
	    TGraphAsymmErrors* graphInvYieldPi0DalitzSystErrMatCanceled         = new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0Dalitz,xPi0DalitzErr,xPi0DalitzErr,yPi0DalitzSystErrMatCanceled,        yPi0DalitzSystErrMatCanceled);
	    TGraphAsymmErrors* graphInvYieldPi0DalitzSystErrCMatCanceled 	= new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0Dalitz,xPi0DalitzErr,xPi0DalitzErr,yPi0DalitzSystErrCMatCanceled,	yPi0DalitzSystErrCMatCanceled);

	    
	    graphInvYieldPi0MethodsStatErr[kDalitz]		= new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0Dalitz,xPi0DalitzErr,xPi0DalitzErr,yPi0DalitzStatErr,			yPi0DalitzStatErr);
	    graphInvYieldPi0MethodsSystErrA[kDalitz] 		= new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0Dalitz,xPi0DalitzErr,xPi0DalitzErr,yPi0DalitzSystErrA,			yPi0DalitzSystErrA);
	    graphInvYieldPi0MethodsSystErrB[kDalitz]		= new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0Dalitz,xPi0DalitzErr,xPi0DalitzErr,yPi0DalitzSystErrB,			yPi0DalitzSystErrA);
	    graphInvYieldPi0MethodsSystErrC[kDalitz] 		= new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0Dalitz,xPi0DalitzErr,xPi0DalitzErr,yPi0DalitzSystErrC,			yPi0DalitzSystErrA);
	    
	  
	    TGraphAsymmErrors* graphInvYieldPi0PCMSystErr  			= new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0PCM,xPi0PCMErr,xPi0PCMErr,yPi0PCMSystErr,  				yPi0PCMSystErr);
	    TGraphAsymmErrors* graphInvYieldPi0PCMSystErrMatCanceled         	= new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0PCM,xPi0PCMErr,xPi0PCMErr,yPi0PCMSystErrMatCanceled,        		yPi0PCMSystErrMatCanceled);
	    TGraphAsymmErrors* graphInvYieldPi0PCMSystErrCMatCanceled 		= new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0PCM,xPi0PCMErr,xPi0PCMErr,yPi0PCMSystErrCMatCanceled,			yPi0PCMSystErrCMatCanceled);
	   
	    
	    
	    graphInvYieldPi0MethodsStatErr[kPCM]  			= new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0PCM,xPi0PCMErr,xPi0PCMErr,yPi0PCMStatErr,				yPi0PCMStatErr);
	    graphInvYieldPi0MethodsSystErrA[kPCM] 			= new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0PCM,xPi0PCMErr,xPi0PCMErr,yPi0PCMSystErrA,				yPi0PCMSystErrA);
	    graphInvYieldPi0MethodsSystErrB[kPCM] 			= new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0PCM,xPi0PCMErr,xPi0PCMErr,yPi0PCMSystErrB,				yPi0PCMSystErrB);
	    graphInvYieldPi0MethodsSystErrC[kPCM]                  	= new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0PCM,xPi0PCMErr,xPi0PCMErr,yPi0PCMSystErrC,				yPi0PCMSystErrC);
	    
	    
	    TGraphAsymmErrors* graphInvYieldPi0PHOSSystErr  				= new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0PHOS,xPi0PHOSErr,xPi0PHOSErr,yPi0PHOSSystErr,  				yPi0PHOSSystErr);
	   
	    
	    graphInvYieldPi0MethodsStatErr[kPHOS]  			= new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0PHOS,xPi0PHOSErr,xPi0PHOSErr,yPi0PHOSStatErr,				yPi0PHOSStatErr);
	    graphInvYieldPi0MethodsSystErrA[kPHOS] 			= new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0PHOS,xPi0PHOSErr,xPi0PHOSErr,yPi0PHOSSystErrA,				yPi0PHOSSystErrA);
	    graphInvYieldPi0MethodsSystErrB[kPHOS] 			= new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0PHOS,xPi0PHOSErr,xPi0PHOSErr,yPi0PHOSSystErrB,				yPi0PHOSSystErrB);
	    graphInvYieldPi0MethodsSystErrC[kPHOS] 			= new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0PHOS,xPi0PHOSErr,xPi0PHOSErr,yPi0PHOSSystErrC,				yPi0PHOSSystErrC);
	   
	    TGraphAsymmErrors* graphInvYieldPi0EMCALSystErr  				= new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0EMCAL,xPi0EMCALErr,xPi0EMCALErr,yPi0EMCALSystErr,  			yPi0EMCALSystErr);
	   
	    
	    graphInvYieldPi0MethodsStatErr[kEMCAL]  			= new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0EMCAL,xPi0EMCALErr,xPi0EMCALErr,yPi0EMCALStatErr,			yPi0EMCALStatErr);
	    graphInvYieldPi0MethodsSystErrA[kEMCAL] 			= new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0EMCAL,xPi0EMCALErr,xPi0EMCALErr,yPi0EMCALSystErrA,			yPi0EMCALSystErrA);
	    graphInvYieldPi0MethodsSystErrB[kEMCAL] 			= new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0EMCAL,xPi0EMCALErr,xPi0EMCALErr,yPi0EMCALSystErrB,			yPi0EMCALSystErrB);
	    graphInvYieldPi0MethodsSystErrC[kEMCAL] 			= new TGraphAsymmErrors(Ntotal,xValComb,yieldPi0EMCAL,xPi0EMCALErr,xPi0EMCALErr,yPi0EMCALSystErrC,			yPi0EMCALSystErrC);
	   
	    
	    TGraphAsymmErrors* ratioDalitzPCM		= (TGraphAsymmErrors*) CalculateRatio(graphInvYieldPi0MethodsStatErr[kDalitz],	graphInvYieldPi0DalitzSystErrMatCanceled,graphInvYieldPi0MethodsStatErr[kPCM],graphInvYieldPi0PCMSystErrMatCanceled,ratioABRange[0][0],ratioABRange[0][1]);
	    TGraphAsymmErrors* ratioDalitzPHOS  	= (TGraphAsymmErrors*) CalculateRatio(graphInvYieldPi0MethodsStatErr[kDalitz],	graphInvYieldPi0DalitzSystErr,graphInvYieldPi0MethodsStatErr[kPHOS],graphInvYieldPi0PHOSSystErr,ratioABRange[1][0],ratioABRange[1][1]);
	    TGraphAsymmErrors* ratioDalitzEMCAL 	= (TGraphAsymmErrors*) CalculateRatio(graphInvYieldPi0MethodsStatErr[kDalitz],	graphInvYieldPi0DalitzSystErr,graphInvYieldPi0MethodsStatErr[kEMCAL],graphInvYieldPi0EMCALSystErr,ratioABRange[2][0],ratioABRange[2][1]);
	    TGraphAsymmErrors* ratioPCMPHOS  		= (TGraphAsymmErrors*) CalculateRatio(graphInvYieldPi0MethodsStatErr[kPCM],	graphInvYieldPi0PCMSystErr,graphInvYieldPi0MethodsStatErr[kPHOS],graphInvYieldPi0PCMSystErr,ratioABRange[3][0],ratioABRange[3][1]);
	    TGraphAsymmErrors* ratioPCMEMCAL 		= (TGraphAsymmErrors*) CalculateRatio(graphInvYieldPi0MethodsStatErr[kPCM],	graphInvYieldPi0PCMSystErr,graphInvYieldPi0MethodsStatErr[kEMCAL],graphInvYieldPi0EMCALSystErr,ratioABRange[4][0],ratioABRange[4][1]);
            TGraphAsymmErrors* ratioPHOSEMCAL	 	= (TGraphAsymmErrors*) CalculateRatio(graphInvYieldPi0MethodsStatErr[kPHOS],	graphInvYieldPi0PHOSSystErr,graphInvYieldPi0MethodsStatErr[kEMCAL],graphInvYieldPi0EMCALSystErr,ratioABRange[5][0],ratioABRange[5][1]);
	  
	    TGraphAsymmErrors* graphRatioAB[6];
	    
	    graphRatioAB[0] = ratioDalitzPCM;
	    graphRatioAB[1] = ratioDalitzPHOS;
	    graphRatioAB[2] = ratioDalitzEMCAL;
	    graphRatioAB[3] = ratioPCMPHOS;
	    graphRatioAB[4] = ratioPCMEMCAL;
	    graphRatioAB[5] = ratioPHOSEMCAL;
	    
	    
	
	
	for(Int_t ii = 0; ii < 6; ii ++){
	  
	Int_t a = ratioAB[ii][0];
	Int_t b = ratioAB[ii][1];
	  
	g_methodA_stat[ii] =    (TGraphAsymmErrors*)graphInvYieldPi0MethodsStatErr[a]->Clone(Form("Pi0_%s_statErr",methodName[a].Data()));
	
        g_methodA_sysA[ii] =    (TGraphAsymmErrors*)graphInvYieldPi0MethodsSystErrA[a]->Clone(Form("Pi0_%s_statAErr",methodName[a].Data()));
	g_methodA_sysB[ii] =    (TGraphAsymmErrors*)graphInvYieldPi0MethodsSystErrB[a]->Clone(Form("Pi0_%s_sysBErr",methodName[a].Data()));
	if( ii == 0) //Material buget canceled out
        g_methodA_sysC[ii] =    (TGraphAsymmErrors*)graphInvYieldPi0DalitzSystErrCMatCanceled->Clone("Pi0_Dalitz_sysCErrMatCanceled");
	else
	g_methodA_sysC[ii] =    (TGraphAsymmErrors*)graphInvYieldPi0MethodsSystErrC[a]->Clone(Form("Pi0_%s_sysCErr",methodName[a].Data()));
	
	g_methodB_stat[ii] =    (TGraphAsymmErrors*)graphInvYieldPi0MethodsStatErr[b]->Clone(Form("Pi0_%s_statErr",methodName[a].Data()));
        g_methodB_sysA[ii] =    (TGraphAsymmErrors*)graphInvYieldPi0MethodsSystErrA[b]->Clone(Form("Pi0_%s_statAErr",methodName[a].Data()));
	g_methodB_sysB[ii] =    (TGraphAsymmErrors*)graphInvYieldPi0MethodsSystErrB[b]->Clone(Form("Pi0_%s_sysBErr",methodName[a].Data()));
	if( ii == 0) //Material buget canceled out
        g_methodB_sysC[ii] =    (TGraphAsymmErrors*)graphInvYieldPi0PCMSystErrCMatCanceled->Clone("Pi0_PCM_sysCErrMatCanceled");
	else
	g_methodB_sysC[ii] =    (TGraphAsymmErrors*)graphInvYieldPi0MethodsSystErrC[b]->Clone(Form("Pi0_%s_sysCErr",methodName[a].Data()));
	
	
	
		  
	  
	}
	
	
	
	ComputePValue(binBybin);
	
	const char *tablePvaluesbinbybin = Form("%s/tablePvaluesbinbybin.dat",outputDir.Data());
	fstream tablePvaluesbinbybinDat;
	cout << tablePvaluesbinbybin << endl;
	tablePvaluesbinbybinDat.open(tablePvaluesbinbybin, ios::out);
	
	
	
	for(Int_t iPt = 0; iPt < nPtLimits-1; iPt ++ ) {
	  
	        tablePvaluesbinbybinDat<<xValComb[iPt]<<"\t&\t";
	       
		for(Int_t ii = 0 ; ii < 6;  ii++){
		  
		  if( iPt < ratioABRange[ii][0] || iPt > ratioABRange[ii][1] ) {
		    if( ii < 5 )
		    tablePvaluesbinbybinDat<<"\t&\t";
		    else
		    tablePvaluesbinbybinDat<<"\t \\\\";
		  }
		  else {
		    if ( ii < 5 )
		    tablePvaluesbinbybinDat<<pValue[ii][iPt][0]<<"\t&\t";
		    else
		    tablePvaluesbinbybinDat<<pValue[ii][iPt][0]<<"\t\\\\";
		  }
		  
		}
		tablePvaluesbinbybinDat<<"\n";
	}  
	
	
	
	tablePvaluesbinbybinDat.close();
	
	
	
	const char *tablePvaluesANDsigma = Form("%s/tablePvaluesANDsigma.dat",outputDir.Data());
	fstream tablePvaluesANDsigmaDat;
	cout << tablePvaluesANDsigma << endl;
	tablePvaluesANDsigmaDat.open(tablePvaluesANDsigma, ios::out);
	
	
	for(Int_t ii = 0 ; ii < 6;  ii++){
	  
	  
		 //tablePvaluesANDsigmaDat << 
		 
		 
		   tablePvaluesANDsigmaDat<<methodName[ratioAB[ii][0]].Data()<<"/"<<methodName[ratioAB[ii][1]].Data()<<"\t&\t"<<xValComb[ratioABRange[ii][0]]<<" - "<< xValComb[ratioABRange[ii][1]]<<"\t&\t"<<  pValue[ii][nPtLimits][0] << "\t&\t" << pValue[ii][nPtLimits][1] <<"\\\\"<<endl;
	
		  
		   //if( iPt < ratioABRange[ii][0] || iPt > ratioABRange[ii][1] ) {
		   // if( ii < 5 )
		    //tablePvaluesbinbybinDat<<"\t&\t";
		    //else
		    //tablePvaluesbinbybinDat<<"\t \\\\";
	}
	
	tablePvaluesANDsigmaDat.close();
	
	
	
	
	
	  
            
	  ////////////////////Ploting////////////////////////////////////////////////////////////////////////////////////////////
	
	//TLatex *labelCollisionSystem = new TLatex(0.3,0.9,"p-Pb,");
	
	TLatex *labelCollisionSystempPb = new TLatex(0.16,0.9,collisionSystempPb.Data());
        SetStyleTLatex( labelCollisionSystempPb, 0.06,4);
	
	TH2F *  histo2DpPb5023GeVRatioDalitzPCM = new TH2F("histo2DpPb5023GeVRatioDalitzPCM","histo2DpPb5023GeVRatioDalitzPCM",1000,0.3,40.,1000,0.2,4.	);
	histo2DpPb5023GeVRatioDalitzPCM->GetYaxis()->SetRangeUser(0.6,2.1);
	histo2DpPb5023GeVRatioDalitzPCM->GetXaxis()->SetRangeUser(0.,20.);
	histo2DpPb5023GeVRatioDalitzPCM->GetXaxis()->SetLabelOffset(-0.015);
	SetStyleHistoTH2ForGraphs(histo2DpPb5023GeVRatioDalitzPCM, "#it{p}_{T} (GeV/#it{c})","#pi^{0} Dalitz-PCM/#pi^{0} PCM", 0.05,0.064, 0.05,0.06, 0.8,0.9, 512, 505);
	
	
	
	
	TCanvas* canvaspPb5023GeVRatioDalitzPCM = new TCanvas("canvaspPb5023GeVRatioDalitzPCM","",200,10,700,500);  // gives the page size
	DrawGammaCanvasSettings( canvaspPb5023GeVRatioDalitzPCM,  0.12, 0.02, 0.02, 0.12);
	
 	canvaspPb5023GeVRatioDalitzPCM->SetLogx();
	histo2DpPb5023GeVRatioDalitzPCM->GetXaxis()->SetRangeUser(0.55,10.5);
	histo2DpPb5023GeVRatioDalitzPCM->GetYaxis()->SetRangeUser(0.0,2.0);
	histo2DpPb5023GeVRatioDalitzPCM->DrawCopy();
	
	 TLatex *labelPValueDalitzPCM = new TLatex(0.7,0.8,Form("p-Value: %2.3f",pValue[0][nPtLimits][0]));
	 TLatex *labelPValueSigmaDalitzPCM = new TLatex(0.7,0.7,Form("#sigma: %2.3f",pValue[0][nPtLimits][1]));
	 SetStyleTLatex( labelPValueDalitzPCM, 0.05,4);
	  SetStyleTLatex( labelPValueSigmaDalitzPCM, 0.05,4);
	
	

		DrawGammaSetMarkerTGraphAsym(ratioDalitzPCM, 20, 1.5,kBlack,kBlack);
		ratioDalitzPCM->Draw("E1psame");
		labelPValueDalitzPCM->Draw("sames");
		labelPValueSigmaDalitzPCM->Draw("sames");
		labelCollisionSystempPb->Draw("sames");
		
		DrawGammaLines(0., 15. , 1, 1 ,1,kBlack);
	
		
	
	
	canvaspPb5023GeVRatioDalitzPCM->Update();
	canvaspPb5023GeVRatioDalitzPCM->Print(Form("%s/RatioDalitzPCM_pPb5023GeV.%s",outputDir.Data(),suffix.Data()));
	delete canvaspPb5023GeVRatioDalitzPCM;
	
	
	TH2F *  histo2DpPb5023GeVRatioDalitzPHOS = new TH2F("histo2DpPb5023GeVRatioDalitzPHOS","histo2DpPb5023GeVRatioDalitzPHOS",1000,0.3,40.,1000,0.2,4.	);
	histo2DpPb5023GeVRatioDalitzPHOS->GetYaxis()->SetRangeUser(0.6,2.1);
	histo2DpPb5023GeVRatioDalitzPHOS->GetXaxis()->SetRangeUser(0.,20.);
	histo2DpPb5023GeVRatioDalitzPHOS->GetXaxis()->SetLabelOffset(-0.015);
	SetStyleHistoTH2ForGraphs(histo2DpPb5023GeVRatioDalitzPHOS, "#it{p}_{T} (GeV/#it{c})","#pi^{0} Dalitz-PCM/#pi^{0} PHOS", 0.05,0.064, 0.05,0.06, 0.8,0.9, 512, 505);
	
	
	
	
	TCanvas* canvaspPb5023GeVRatioDalitzPHOS = new TCanvas("canvaspPb5023GeVRatioDalitzPHOS","",200,10,700,500);  // gives the page size
	DrawGammaCanvasSettings( canvaspPb5023GeVRatioDalitzPHOS,  0.12, 0.02, 0.02, 0.12);
	
 	canvaspPb5023GeVRatioDalitzPHOS->SetLogx();
	histo2DpPb5023GeVRatioDalitzPHOS->GetXaxis()->SetRangeUser(0.55,10.5);
	histo2DpPb5023GeVRatioDalitzPHOS->GetYaxis()->SetRangeUser(0.0,2.0);
	histo2DpPb5023GeVRatioDalitzPHOS->DrawCopy();
	
	TLatex *labelPValueDalitzPHOS = new TLatex(0.7,0.8,Form("p-Value: %2.3f",pValue[1][nPtLimits][0]));
	 TLatex *labelPValueSigmaDalitzPHOS = new TLatex(0.7,0.7,Form("#sigma: %2.3f",pValue[1][nPtLimits][1]));
	 SetStyleTLatex( labelPValueDalitzPHOS, 0.05,4);
	  SetStyleTLatex( labelPValueSigmaDalitzPHOS, 0.05,4);
	
	      
	

		DrawGammaSetMarkerTGraphAsym(ratioDalitzPHOS, 20, 1.5,kBlack,kBlack);
		ratioDalitzPHOS->Draw("E1psame");
		
		labelPValueDalitzPHOS->Draw("sames");
		labelPValueSigmaDalitzPHOS->Draw("sames");
		labelCollisionSystempPb->Draw("sames");

		
		DrawGammaLines(0., 15. , 1, 1 ,1,kBlack);
	
		
	
	
	canvaspPb5023GeVRatioDalitzPHOS->Update();
	canvaspPb5023GeVRatioDalitzPHOS->Print(Form("%s/RatioDalitzPHOS_pPb5023GeV.%s",outputDir.Data(),suffix.Data()));
	delete canvaspPb5023GeVRatioDalitzPHOS;
	
	
	TH2F *  histo2DpPb5023GeVRatioDalitzEMCAL = new TH2F("histo2DpPb5023GeVRatioDalitzEMCAL","histo2DpPb5023GeVRatioDalitzEMCAL",1000,0.3,40.,1000,0.2,4.	);
	histo2DpPb5023GeVRatioDalitzEMCAL->GetYaxis()->SetRangeUser(0.6,2.1);
	histo2DpPb5023GeVRatioDalitzEMCAL->GetXaxis()->SetRangeUser(0.,20.);
	histo2DpPb5023GeVRatioDalitzEMCAL->GetXaxis()->SetLabelOffset(-0.015);
	SetStyleHistoTH2ForGraphs(histo2DpPb5023GeVRatioDalitzEMCAL, "#it{p}_{T} (GeV/#it{c})","#pi^{0} Dalitz-PCM /#pi^{0} EMCAL", 0.05,0.064, 0.05,0.06, 0.8,0.9, 512, 505);
	
	
	
	TCanvas* canvaspPb5023GeVRatioDalitzEMCAL = new TCanvas("canvaspPb5023GeVRatioDalitzEMCAL","",200,10,700,500);  // gives the page size
	DrawGammaCanvasSettings( canvaspPb5023GeVRatioDalitzEMCAL,  0.12, 0.02, 0.02, 0.12);
	
 	canvaspPb5023GeVRatioDalitzEMCAL->SetLogx();
	histo2DpPb5023GeVRatioDalitzEMCAL->GetXaxis()->SetRangeUser(0.55,10.5);
	histo2DpPb5023GeVRatioDalitzEMCAL->GetYaxis()->SetRangeUser(0.0,2.0);
	histo2DpPb5023GeVRatioDalitzEMCAL->DrawCopy();
	
	TLatex *labelPValueDalitzEMCAL = new TLatex(0.7,0.8,Form("p-Value: %2.3f",pValue[2][nPtLimits][0]));
	 TLatex *labelPValueSigmaDalitzEMCAL = new TLatex(0.7,0.7,Form("#sigma: %2.3f",pValue[2][nPtLimits][1]));
	 SetStyleTLatex( labelPValueDalitzEMCAL, 0.05,4);
	  SetStyleTLatex( labelPValueSigmaDalitzEMCAL, 0.05,4);
	
	  
	  
	

		DrawGammaSetMarkerTGraphAsym(ratioDalitzEMCAL, 20, 1.5,kBlack,kBlack);
		ratioDalitzEMCAL->Draw("E1psame");
		labelPValueDalitzEMCAL->Draw("sames");
		labelPValueSigmaDalitzEMCAL->Draw("sames");
		  labelCollisionSystempPb->Draw("sames");
	  
	

		
		DrawGammaLines(0., 15. , 1, 1 ,1,kBlack);
	
		
	
	
	canvaspPb5023GeVRatioDalitzEMCAL->Update();
	canvaspPb5023GeVRatioDalitzEMCAL->Print(Form("%s/RatioDalitzEMCAL_pPb5023GeV.%s",outputDir.Data(),suffix.Data()));
	delete canvaspPb5023GeVRatioDalitzEMCAL;
	
	
	
	TH2F *  histo2DpPb5023GeVRatioPCMPHOS = new TH2F("histo2DpPb5023GeVRatioPCMPHOS","histo2DpPb5023GeVRatioPCMPHOS",1000,0.3,40.,1000,0.2,4.	);
	histo2DpPb5023GeVRatioPCMPHOS->GetYaxis()->SetRangeUser(0.6,2.1);
	histo2DpPb5023GeVRatioPCMPHOS->GetXaxis()->SetRangeUser(0.,20.);
	histo2DpPb5023GeVRatioPCMPHOS->GetXaxis()->SetLabelOffset(-0.015);
	SetStyleHistoTH2ForGraphs(histo2DpPb5023GeVRatioPCMPHOS, "#it{p}_{T} (GeV/#it{c})","#pi^{0} PCM /#pi^{0} PHOS", 0.05,0.064, 0.05,0.06, 0.8,0.9, 512, 505);
	
	
	
	TCanvas* canvaspPb5023GeVRatioPCMPHOS = new TCanvas("canvaspPb5023GeVRatioPCMPHOS","",200,10,700,500);  // gives the page size
	DrawGammaCanvasSettings( canvaspPb5023GeVRatioPCMPHOS,  0.12, 0.02, 0.02, 0.12);
	
 	canvaspPb5023GeVRatioPCMPHOS->SetLogx();
	histo2DpPb5023GeVRatioPCMPHOS->GetXaxis()->SetRangeUser(0.55,10.5);
	histo2DpPb5023GeVRatioPCMPHOS->GetYaxis()->SetRangeUser(0.0,2.0);
	histo2DpPb5023GeVRatioPCMPHOS->DrawCopy();
	
	TLatex *labelPValuePCMPHOS = new TLatex(0.7,0.8,Form("p-Value: %2.3f",pValue[3][nPtLimits][0]));
	 TLatex *labelPValueSigmaPCMPHOS = new TLatex(0.7,0.7,Form("#sigma: %2.3f",pValue[3][nPtLimits][1]));
	 SetStyleTLatex( labelPValuePCMPHOS, 0.05,4);
	  SetStyleTLatex( labelPValueSigmaPCMPHOS, 0.05,4);
	
	

		DrawGammaSetMarkerTGraphAsym(ratioPCMPHOS, 20, 1.5,kBlack,kBlack);
		ratioPCMPHOS->Draw("E1psame");
		
		labelPValuePCMPHOS->Draw("sames");
		labelPValueSigmaPCMPHOS->Draw("sames");
		  labelCollisionSystempPb->Draw("sames");
	  

		
		DrawGammaLines(0., 15. , 1, 1 ,1,kBlack);
	
		
	
	
	canvaspPb5023GeVRatioPCMPHOS->Update();
	canvaspPb5023GeVRatioPCMPHOS->Print(Form("%s/RatioPCMPHOS_pPb5023GeV.%s",outputDir.Data(),suffix.Data()));
	delete canvaspPb5023GeVRatioPCMPHOS;
	
	
	
	TH2F *  histo2DpPb5023GeVRatioPCMEMCAL = new TH2F("histo2DpPb5023GeVRatioPCMEMCAL","histo2DpPb5023GeVRatioPCMEMCAL",1000,0.3,40.,1000,0.2,4.	);
	histo2DpPb5023GeVRatioPCMEMCAL->GetYaxis()->SetRangeUser(0.6,2.1);
	histo2DpPb5023GeVRatioPCMEMCAL->GetXaxis()->SetRangeUser(0.,20.);
	histo2DpPb5023GeVRatioPCMEMCAL->GetXaxis()->SetLabelOffset(-0.015);
	SetStyleHistoTH2ForGraphs(histo2DpPb5023GeVRatioPCMEMCAL, "#it{p}_{T} (GeV/#it{c})","#pi^{0 }PCM /#pi^{0} EMCAL", 0.05,0.064, 0.05,0.06, 0.8,0.9, 512, 505);
	
	
	
	TCanvas* canvaspPb5023GeVRatioPCMEMCAL = new TCanvas("canvaspPb5023GeVRatioPCMEMCAL","",200,10,700,500);  // gives the page size
	DrawGammaCanvasSettings( canvaspPb5023GeVRatioPCMEMCAL,  0.12, 0.02, 0.02, 0.12);
	
 	canvaspPb5023GeVRatioPCMEMCAL->SetLogx();
	histo2DpPb5023GeVRatioPCMEMCAL->GetXaxis()->SetRangeUser(0.55,10.5);
	histo2DpPb5023GeVRatioPCMEMCAL->GetYaxis()->SetRangeUser(0.0,2.0);
	histo2DpPb5023GeVRatioPCMEMCAL->DrawCopy();
	
	TLatex *labelPValuePCMEMCAL = new TLatex(0.7,0.8,Form("p-Value: %2.3f",pValue[4][nPtLimits][0]));
	 TLatex *labelPValueSigmaPCMEMCAL = new TLatex(0.7,0.7,Form("#sigma: %2.3f",pValue[4][nPtLimits][1]));
	 SetStyleTLatex( labelPValuePCMEMCAL, 0.05,4);
	  SetStyleTLatex( labelPValueSigmaPCMEMCAL, 0.05,4);
	
	

		DrawGammaSetMarkerTGraphAsym(ratioPCMEMCAL, 20, 1.5,kBlack,kBlack);
		ratioPCMEMCAL->Draw("E1psame");
		
		labelPValuePCMEMCAL->Draw("sames");
		labelPValueSigmaPCMEMCAL->Draw("sames");
		  labelCollisionSystempPb->Draw("sames");

		
		DrawGammaLines(0., 15. , 1, 1 ,1,kBlack);
	
		
	
	
	canvaspPb5023GeVRatioPCMEMCAL->Update();
	canvaspPb5023GeVRatioPCMEMCAL->Print(Form("%s/RatioPCMEMCAL_pPb5023GeV.%s",outputDir.Data(),suffix.Data()));
	delete canvaspPb5023GeVRatioPCMEMCAL;
	
	
	
	TH2F *  histo2DpPb5023GeVRatioPHOSEMCAL = new TH2F("histo2DpPb5023GeVRatioPHOSEMCAL","histo2DpPb5023GeVRatioPHOSEMCAL",1000,0.3,40.,1000,0.2,4.	);
	histo2DpPb5023GeVRatioPHOSEMCAL->GetYaxis()->SetRangeUser(0.6,2.1);
	histo2DpPb5023GeVRatioPHOSEMCAL->GetXaxis()->SetRangeUser(0.,20.);
	histo2DpPb5023GeVRatioPHOSEMCAL->GetXaxis()->SetLabelOffset(-0.015);
	SetStyleHistoTH2ForGraphs(histo2DpPb5023GeVRatioPHOSEMCAL, "#it{p}_{T} (GeV/#it{c})","#pi^{0 }PHOS /#pi^{0} EMCAL", 0.05,0.064, 0.05,0.06, 0.8,0.9, 512, 505);
	
	
	
	TCanvas* canvaspPb5023GeVRatioPHOSEMCAL = new TCanvas("canvaspPb5023GeVRatioPHOSEMCAL","",200,10,700,500);  // gives the page size
	DrawGammaCanvasSettings( canvaspPb5023GeVRatioPHOSEMCAL,  0.12, 0.02, 0.02, 0.12);
	
 	canvaspPb5023GeVRatioPHOSEMCAL->SetLogx();
	histo2DpPb5023GeVRatioPHOSEMCAL->GetXaxis()->SetRangeUser(0.55,10.5);
	histo2DpPb5023GeVRatioPHOSEMCAL->GetYaxis()->SetRangeUser(0.0,2.0);
	histo2DpPb5023GeVRatioPHOSEMCAL->DrawCopy();
	
	
	TLatex *labelPValuePHOSEMCAL = new TLatex(0.7,0.8,Form("p-Value: %2.3f",pValue[5][nPtLimits][0]));
	TLatex *labelPValueSigmaPHOSEMCAL = new TLatex(0.7,0.7,Form("#sigma: %2.3f",pValue[5][nPtLimits][1]));
	SetStyleTLatex( labelPValuePHOSEMCAL, 0.05,4);
	SetStyleTLatex( labelPValueSigmaPHOSEMCAL, 0.05,4);
	

	DrawGammaSetMarkerTGraphAsym(ratioPHOSEMCAL, 20, 1.5,kBlack,kBlack);
	ratioPHOSEMCAL->Draw("E1psame");
		
		
	labelPValuePHOSEMCAL->Draw("sames");
	labelPValueSigmaPHOSEMCAL->Draw("sames");
	labelCollisionSystempPb->Draw("sames");


		
	DrawGammaLines(0., 15. , 1, 1 ,1,kBlack);
	
		
	
	
	canvaspPb5023GeVRatioPHOSEMCAL->Update();
	canvaspPb5023GeVRatioPHOSEMCAL->Print(Form("%s/RatioPHOSEMCAL_pPb5023GeV.%s",outputDir.Data(),suffix.Data()));
	delete canvaspPb5023GeVRatioPHOSEMCAL;
	
	
	
	TCanvas* canvasRatiosAllMethods = new TCanvas("canvasRatiosAllMethods","",3100,1080*2);  // gives the page size
	canvasRatiosAllMethods->SetTickx();
	canvasRatiosAllMethods->SetTicky();
	canvasRatiosAllMethods->SetGridx(0);
	canvasRatiosAllMethods->SetGridy(0);
	canvasRatiosAllMethods->SetLogy(1);
	canvasRatiosAllMethods->SetLeftMargin(0.00);
	canvasRatiosAllMethods->SetRightMargin(0.00);
	canvasRatiosAllMethods->SetTopMargin(0.00);
	canvasRatiosAllMethods->SetFillColor(0);
	
	  TPad* padRatiosAllMethods[6];
	
	
	  Double_t globalleftOffset  = 0.19;
	  Double_t globalrightOffset = 0.02;
	  Double_t globalbottomMarging = 0.17;
	  Double_t globaltopMarging    = 0.02;
          Double_t padWith =   ( 1.0 - (( 1.0 / 3 )*globalleftOffset ) - ( (1.0/3)*globalrightOffset ) ) / 3;
	  Double_t padHeight =  ( 1.0 - (( 1.0 / 2 )*globalbottomMarging ) - ( (1.0/2)*globaltopMarging ) ) / 2;
	   
	  cout<<"padWith: "<<padWith<<" padHeight"<<padHeight<<endl;
	 
	  
	  Double_t xMin = 0;
	  Double_t xMax = 0;
	  Double_t yMin = 0;
	  Double_t yMax = yMin + padHeight + ( ( 1.0 /2 )*globalbottomMarging) ;
	
	for ( Int_t ii = 0;  ii < 6; ii++){
	  
	  Double_t leftOffset   = 0;
	  Double_t rightOffset  = 0;
	  
	  Double_t bottomMargin =  globalbottomMarging;
	  Double_t topMargin   =  0;
	  
	  
	  xMin =  xMax;
	  xMax =  xMin + padWith; //( ( 1.0 /3 )*globalleftOffset) ;
	  
	  if ( ii >= 3 ){
	    if( ii == 3 ){
	    yMin = yMax;
	    yMax = yMin + padHeight;
	    }
	    
	    //bottomMargin =  globalbottomMarging;
	    //topMargin    =  0;
	    
	    bottomMargin = 0;
	    topMargin    = globaltopMarging;
	
	  }
 	  
 	  if( ii == 0 || ii == 3 ){
	    
	    leftOffset  = globalleftOffset;
	    
	     xMin =  0;
	     xMax =  xMin + padWith + ( ( 1.0 /3 )*globalleftOffset) ;
	    
	  }
	  
	  if( ii == 2 || ii == 5 ){
	    
	    rightOffset = globalrightOffset;
	    
	  }
	  
	  cout<<xMin<<"\t"<<yMin<<"\t"<<xMax<<"\t"<<yMax<<endl;
	  
	  
	  padRatiosAllMethods[ii] = new TPad(Form("padRatiosAllMethods_%d",ii),"",xMin,yMin,xMax,yMax,-1, -1, -2);
	  
	  padRatiosAllMethods[ii]->SetFillColor(0);
	  padRatiosAllMethods[ii]->GetFrame()->SetFillColor(0);
	  padRatiosAllMethods[ii]->SetBorderMode(0);
	  padRatiosAllMethods[ii]->SetLeftMargin(leftOffset);
	  padRatiosAllMethods[ii]->SetBottomMargin(bottomMargin);
	  padRatiosAllMethods[ii]->SetRightMargin(rightOffset);
	  padRatiosAllMethods[ii]->SetTopMargin(topMargin);
	  padRatiosAllMethods[ii]->SetLogx(1);
	  padRatiosAllMethods[ii]->Draw();
	  //padRatiosAllMethods[ii]->cd();
	  
	  //DrawGammaSetMarkerTGraphAsym(graphRatioAB[ii], 20, 1.5,kBlack,kBlack);
	  
	  //graphRatioAB[ii]->Draw("E1p");
	  
	}
	
	 TH2F *histoDummy[6];
	
	for(Int_t ii = 0; ii < 6; ii++){
	  
	  padRatiosAllMethods[ii]->cd();
	  
	  
	  
	  histoDummy[ii] = new TH2F("histoFWHMCanvas2D","histoFWHMCanvas2D",1000,0.0,15.,1000,0.0,5.);
	  histoDummy[ii]->GetYaxis()->SetRangeUser(0.6,2.1);
	  histoDummy[ii]->GetXaxis()->SetRangeUser(0.5,15.);
	  SetStyleHistoTH2ForGraphs(histoDummy[ii], "#it{p}_{T} (GeV/#it{c})","DetA/DetB", 0.05,0.064, 0.05,0.06, 0.8,0.9, 512, 505);
	   histoDummy[ii]->GetXaxis()->SetTitleOffset(1.2);
	   histoDummy[ii]->GetYaxis()->SetTitleOffset(1.2);
	   histoDummy[ii]->GetYaxis()->CenterTitle();
	  
	  histoDummy[ii]->Draw();
	  
	  DrawGammaSetMarkerTGraphAsym(graphRatioAB[ii], 20, 2.0,kBlack,kBlack,1.4,kTRUE);
	  
	  graphRatioAB[ii]->Draw("pEsame");
	  
	  DrawGammaLines(0.,15,1.0, 1.0,0.5);
	  
	   Double_t labelleftOffset = 0.05;
	   Double_t labeltopOffset = 0.88;
	   Double_t labelpValuOffset = 0;
	  
	  if( ii == 0 || ii == 3 ){
	    labelleftOffset  = 0.05+globalleftOffset;
	    labelpValuOffset = 0.6;
	  } else{
	    
	     labelpValuOffset = labelleftOffset+0.6;
	  }
	  
	  if( ii < 3 )
	    labeltopOffset =  labeltopOffset + globaltopMarging;//( 1 - globalbottomMarging )*labeltopOffset;
	    
	  
	  TLatex *labelDummy = new TLatex(labelleftOffset,labeltopOffset,collisionSystempPb.Data());
	  SetStyleTLatex( labelDummy,0.05,4);
	  labelDummy->Draw();
	  
	  TLatex *labelpValueDummy = new TLatex(labelpValuOffset,labeltopOffset-0.15,Form("p-Value: %2.3f",pValue[ii][nPtLimits][0]));
	  SetStyleTLatex( labelpValueDummy,0.05,4);
	  labelpValueDummy->Draw();
	  
	   TLatex *labelSigmaDummy = new TLatex(labelpValuOffset,labeltopOffset-0.21,Form("#sigma: %2.3f",pValue[ii][nPtLimits][1]));
	   SetStyleTLatex( labelSigmaDummy, 0.05,4);
	   labelSigmaDummy->Draw();
          
	    TLegend* legendDummy = new TLegend(labelleftOffset,labeltopOffset-0.15,labelleftOffset+0.3,labeltopOffset-.05);
	    legendDummy->SetTextSize(0.04);
	    legendDummy->SetFillColor(0);
	    legendDummy->SetLineColor(0);  
	    legendDummy->SetFillStyle(0);
	    legendDummy->Draw();
	  
	     legendDummy->AddEntry(graphRatioAB[ii],Form("%s/%s",methodName[ratioAB[ii][0]].Data(),methodName[ratioAB[ii][1]].Data()) );
	  
	  
	}
	
	canvasRatiosAllMethods->Update();
	canvasRatiosAllMethods->Print(Form("%s/RatioAllMethods.%s",outputDir.Data(),suffix.Data()));
	delete canvasRatiosAllMethods;
	
	
	
	/*
	TPad* padMassMeson01 = new TPad("padMassMeson01", "",xMin0, 0.0,xMax0, 1.,-1, -1, -2);
        padMassMeson01->SetFillColor(0);
   padMassMeson01->GetFrame()->SetFillColor(0);
   padMassMeson01->SetBorderMode(0);
   padMassMeson01->SetLeftMargin(leftOffset);
   padMassMeson01->SetBottomMargin(bottomMargin);
   padMassMeson01->SetRightMargin(0.00);
   padMassMeson01->SetTopMargin(topMargin);
   //padMassMeson01->SetLogx(1);
   padMassMeson01->Draw();
   
   
   TPad* padMassMeson02 = new TPad("padMassMeson02", "",xMin1, 0.0,xMax1, 1.,-1, -1, -2);
   padMassMeson02->SetFillColor(0);
   padMassMeson02->GetFrame()->SetFillColor(0);
   padMassMeson02->SetBorderMode(0);
   padMassMeson02->SetLeftMargin(0.0);
   padMassMeson02->SetBottomMargin(bottomMargin);
   padMassMeson02->SetRightMargin(0.00);
   padMassMeson02->SetTopMargin(topMargin);
   //padMassMeson02->SetLogx(1);
   padMassMeson02->Draw();
   
   
   TPad* padMassMeson03 = new TPad("padMassMeson03", "",xMin2, 0.0,xMax2, 1.,-1, -1, -2);
   padMassMeson03->SetFillColor(0);
   padMassMeson03->GetFrame()->SetFillColor(0);
   padMassMeson03->SetBorderMode(0);
   padMassMeson03->SetLeftMargin(0.00);
   padMassMeson03->SetBottomMargin(bottomMargin);
   padMassMeson03->SetRightMargin(rightOffset);
   padMassMeson03->SetTopMargin(topMargin);
   padMassMeson03->Draw();  
   //padMassMeson03->SetLogx(1);
   padMassMeson01->cd();
	*/
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	    
	    
	    TFile fResults(Form("ResultspPb_%s.root",dateForOutput.Data()),"RECREATE");
	   
	    ratioPCMPHOS->Write("Ratio_PCM_PHOS");
	    ratioPCMEMCAL->Write("Ratio_PCM_EMCAL");
	    ratioDalitzPHOS->Write("Ratio_Dalitz_PHOS");
	    ratioDalitzEMCAL->Write("Ratio_Dalitz_EMCAL");
	    ratioPHOSEMCAL->Write("Ratio_PHOS_EMCAL");
	    ratioDalitzPCM->Write("Ratio_Dalitz_PCM");
	    
	    
	    
		
	
	
		
}

TGraphAsymmErrors* CalculateRatio(TGraphAsymmErrors* ObjA_StatErr, TGraphAsymmErrors* ObjA_SystErr, TGraphAsymmErrors* ObjB_StatErr, TGraphAsymmErrors* ObjB_SystErr,Int_t startBin, Int_t endBin){
  
  
  //cout<<"Warnig: in this function the binning must be the same"<<endl;
  
  if( startBin < 0  || startBin >= ObjA_StatErr->GetN() || startBin >= ObjB_StatErr->GetN() || startBin > endBin || endBin >= ObjA_StatErr->GetN() || endBin >= ObjB_SystErr->GetN() ){
    
    cout<<"Error: please check the binning of the histograms"<<endl;
    
    return 0x0;
    
  }

  Int_t nPoints = (endBin - startBin) + 1;
  
  
  //Double_t yErrSystLowA[nPoints];
  //Double_t yErrSystHighA[nPoints];
  Double_t yErrTotLowA[nPoints];
  Double_t yErrTotHighA[nPoints];
  
  //Double_t yErrSystLowB[nPoints];
  //Double_t yErrSystHighB[nPoints];
  Double_t yErrTotLowB[nPoints];
  Double_t yErrTotHighB[nPoints];
    
  Double_t xVal[nPoints];
  Double_t xErrLow[nPoints];
  Double_t xErrHigh[nPoints];
  Double_t ratioY[nPoints];
  Double_t yErrSystLow[nPoints];
  Double_t yErrSystHigh[nPoints];
  Double_t yErrStatLow[nPoints];
  Double_t yErrStatHigh[nPoints];
  Double_t yErrTotLow[nPoints];
  Double_t yErrTotHigh[nPoints];
  
    
  Int_t newBin = 0;  
  
  for(Int_t i = startBin;  i <= endBin; i++){
    
        if( ObjA_StatErr->GetX()[i] != ObjB_StatErr->GetX()[i] ){
	  cout<<"Errot: In this function both graph must have the same binning in the range"<<endl;
	  return 0x0;
	}
	
	yErrTotLowA[newBin] 	= pow( ( pow( ObjA_SystErr->GetErrorYlow(i), 2 )  +  pow( ObjA_StatErr->GetErrorYlow(i), 2 ) ), 0.5 );
	yErrTotHighA[newBin]	= pow( ( pow( ObjA_SystErr->GetErrorYhigh(i),2 )  +  pow( ObjA_StatErr->GetErrorYhigh(i),2 ) ), 0.5 );
	
	yErrTotLowB[newBin] 	= pow( ( pow( ObjB_SystErr->GetErrorYlow(i), 2 )  +  pow( ObjB_StatErr->GetErrorYlow(i), 2 ) ), 0.5 );
	yErrTotHighB[newBin]	= pow( ( pow( ObjB_SystErr->GetErrorYhigh(i),2 )  +  pow( ObjB_StatErr->GetErrorYhigh(i),2 ) ), 0.5 );
	
	
	
	
	
        xVal[newBin] 		= ObjA_StatErr->GetX()[i];
	xErrLow[newBin] 	= ObjA_StatErr->GetErrorXlow(i);
	xErrHigh[newBin]	= ObjA_StatErr->GetErrorXhigh(i);
	ratioY[newBin]         =  ObjA_StatErr->GetY()[i]/ObjB_StatErr->GetY()[i];
	
	
	//yErrStatLow[newBin]     = pow( ( pow( ObjA_StatErr->GetErrorYlow(i)/ObjA_StatErr->GetY()[i] ,2) +  pow( ObjB_StatErr->GetErrorYlow(i)/ObjB_StatErr->GetY()[i] ,2) ),0.5  )*rationY[newBin];
	//yErrStatHigh[newBin]    = pow( ( pow( ObjA_StatErr->GetErrorYhigh(i)/ObjA_StatErr->GetY()[i] ,2) +  pow( ObjB_StatErr->GetErrorYhigh(i)/ObjB_StatErr->GetY()[i] ,2) ),0.5  )*rationY[newBin];
	
	
	
	
	//yErrSystLow[newBin]     = pow( ( pow( ObjA_SystErr->GetErrorYlow(i)/ObjA_SystErr->GetY()[i] ,2) +   pow( ObjB_SystErr->GetErrorYlow(i)/ObjB_SystErr->GetY()[i] ,2) ),0.5 )*rationY[newBin];
	//yErrSystHigh[newBin]    = pow( ( pow( ObjA_SystErr->GetErrorYhigh(i)/ObjA_SystErr->GetY()[i] ,2) +  pow( ObjB_SystErr->GetErrorYhigh(i)/ObjB_SystErr->GetY()[i] ,2) ),0.5  )*rationY[newBin];
	
	yErrTotLow[newBin]  = pow( ( pow( yErrTotLowA[newBin]/ObjA_SystErr->GetY()[i], 2) + pow(yErrTotLowB[newBin]/ObjB_SystErr->GetY()[i],2)  ),0.5)*ratioY[newBin];
	yErrTotHigh[newBin] = pow( ( pow( yErrTotHighA[newBin]/ObjA_SystErr->GetY()[i], 2) + pow(yErrTotHighB[newBin]/ObjB_SystErr->GetY()[i],2) ), 0.5)*ratioY[newBin];
	
	
	newBin++;
	
    
  }
    
  TGraphAsymmErrors* graphRatio = new TGraphAsymmErrors(newBin,xVal,ratioY,xErrLow,xErrHigh,yErrTotLow,yErrTotHigh);
  
  return graphRatio;
  
}

TGraphAsymmErrors* GetSpectrumWithSystErrors(TH1D* spectrum, TString fileNameSysErr,TString outputName,Int_t offsetSyst){
  

  
  ifstream fileSysErr;
 
  
  Double_t relSystErrorUp[50];
  Double_t relSystErrorDown[50];
 
  
    
	fileSysErr.open(fileNameSysErr.Data(),ios_base::in);
	cout << fileNameSysErr.Data() << endl;
	Int_t nPoints = 0;

	while(!fileSysErr.eof()){
		fileSysErr >> relSystErrorDown[nPoints] >> relSystErrorUp[nPoints];
		cout << nPoints << "\t"  << relSystErrorDown[nPoints] << "\t"  <<relSystErrorUp[nPoints] << endl;;		
		nPoints++;
	}
	fileSysErr.close();
	nPoints = nPoints-1;
		
	TGraphAsymmErrors* graphCorrectedYieldSysErr = CalculateSysErrFromRelSysHisto( spectrum, outputName.Data(), relSystErrorDown , relSystErrorUp, offsetSyst, nPoints);
	
	return graphCorrectedYieldSysErr;
	
    
}
