
/****************************************************************************************************************************
******          provided by Gamma Conversion Group, PWG-GA*****
******          Pedro Gonzalez, pedro.gonzalez.zamora@cern.ch
******          Ana Marin, marin@physi.uni-heidelberg.de    
*****           Annika Passfeld, annikapassfeld@uni-muenster,de                                                          *****
******          Friederike Bock, friederike.bock@cern.ch                                                                 *****
******          
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
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"
#include "CommonHeaders/Interpolation5023GeV.h"
#include "CommonHeaders/PlottingInterPolationRpPb5023GeV.h"
#include "TFitResultPtr.h"

  


void InterpolationpPb5023GeV(TString System="PCM", TString resultsType ="PCMPileUpCorrection", TString suffix="pdf", TString outputDir="OutputRpPb_yShift_Individual",TString FitFuncName="Tsallis", TString thesisPlots="", Int_t fixParam=1)//Bylinkin
{

	
	gROOT->SetStyle("Plain");
        TH1::AddDirectory(kFALSE);
	
	if(thesisPlots.EqualTo("thesis") ){  
	StyleSettingsThesis(suffix);
	cout<<"Entre a thesis plots"<<endl;
	}else {
	StyleSettings();
	}
	
	
	//Defining colors
	
	colorsArray[kDalitz] = kCyan+2;
	colorsArray[kPCM]    = kBlack;
	colorsArray[kPHOS]   = kRed+1;
	colorsArray[kEMCal]   = kGreen+2;
	colorsArray[kPion]   = kGreen-3;
	colorsArray[kKaon]   = kOrange-3;
	colorsArray[kProton] = kMagenta-3;
	
	Float_t pTLimit=20.;//30
	Float_t pTLimitLow=0.3;//0.8
	TString SystemLabel= CaloLabel.Data();
	if (System.CompareTo("PCM")==0){
	  //	  pTLimitLow=0.3;	
	  //	  pTLimit=14.5.;
	  SystemLabel=PCMLabel.Data();
	}
	if (System.CompareTo("Dalitz")==0 || System.CompareTo("DALITZ")==0){
	  //	  pTLimit=10.5;
	  //	  pTLimitLow=0.3;
	  SystemLabel=DalitzLabel.Data();
	}
	if (System.CompareTo("PHOS")==0 ){

	  //	  pTLimitLow=0.8;
	}
	
	if ( System.CompareTo("EMCal") == 0 || System.CompareTo("EMCAL") == 0 || FitFuncName.CompareTo("Bylinkin") == 0 ){
	
	  cout<<"WARNING: Fix parameter not applied for Bylinkin or for EMCAl system";
	  
	  fixParam = -1;
	    
	}
	  
	  
	cout << "pT Limit:  "<<pTLimit << endl;
	cout << " low pT Limit:  "<<pTLimitLow << endl;
	cout << "SytemsLabel:  "<<SystemLabel.Data() << endl;
	TString dateForOutput = ReturnDateStringForOutput();
	
	if( fixParam < 0 ) { 
	  
	  outputDir = Form("%s/%s/%s/%s/%s/%s",outputDir.Data(),dateForOutput.Data(),System.Data(),suffix.Data(),resultsType.Data(),FitFuncName.Data());
	
	} else {
	  TString funcName = Form("%sParam_n_Fixed",FitFuncName.Data());  
	  outputDir = Form("%s/%s/%s/%s/%s/%s",outputDir.Data(),dateForOutput.Data(),System.Data(),suffix.Data(),resultsType.Data(),funcName.Data());
	}
	gSystem->Exec("mkdir -p "+outputDir);
	
	Bool_t PileUpCorrection = kFALSE;
	Bool_t thesis = kFALSE;

	//Input Files GA	
	TString fileNameNeutralPionCombResultsPP	="FinalResults/CombinedResultsPP_ShiftedX_PaperRAA_16_May_2014.root";
	TString fileNameNeutralPionEMCalResultsPP	="ExternalInputpPb/EMCAL/data_EMCAL-EMCALResultsFullCorrection_PP_276_Friederike.root";
	TString fileNameNeutralPionEMCalResultsPP7TeV	="ExternalInputpPb/EMCAL/pi0Specrtum2011EMCAL_24June2015_7_Evi.root"; 
	TString fileNameNeutralPionPHOSResultsPP	="ExternalInputpPb/PHOS/CombinedResultsPP_ShiftedX_PaperRAA_ConsiderPileup7TeVPHOSData.root";
	TString fileNameNeutralPionCombResultspPb	="ExternalInputpPb/InputRpPb/ResultspPb_Tsallis_2016_06_06.root"; //low pT Cut
	TString fileNameNeutralPionPCMResultspPb        ="ExternalInputpPb/PCM/data_PCMResults_pPb_20151111_standard_CatErrors.root";
	TString fileNameNeutralPionDalitzResultspPb	="ExternalInputpPb/PCM/data_PCMResults_Dalitz_pPb_20160929.root";//data_PCMResults_Dalitz_pPb_20150806.root";data_PCMResults_Dalitz_pPb_20160601.root
	TString fileNameNeutralPionEMCalResultspPb	="ExternalInputpPb/EMCAL/data_EMCalEMCalResults_160602_newhighptbinningPi0_pPb.root";//data_EMCalEMCalResults_160215_pPb_5023_Mike.root";
	TString fileNameNeutralPionPHOSResultspPb       ="ExternalInputpPb/PHOS/20160601_Pi0InvariantSpectrum_pPb_PHOS.root";//InvariantYield_Pi0_pPb_Graph_PHOS_20160415.root";
	TString fileNameRpPbPHOS  	                ="ExternalInputpPb/PHOS/data_PHOSResults_RpPb_20160405.root";
	TString fileNamePHOSSystErrCancellation         ="ExternalInputpPb/PHOS/ComponentCancelSys.root";
	//	TFile* CommonFile=new TFile("ExternalInputpPb/Combined_pPbResults_2016_03_31.root");// Common Fit for YShift
	TFile* CommonFile=new TFile("ExternalInputpPb/InputRpPb/ResultspPb_Tsallis_2016_10_07.root");// Common Fit for YShift
	
	TString nameRebinSpectraFitsDat = Form("%s/RebinSpectraFitsParam.dat",outputDir.Data());
 	fstream  fileRebinSpectraFits;
 	fileRebinSpectraFits.open(nameRebinSpectraFitsDat.Data(), ios::out);  

	//PHOS Input

	TString NamePHOSRpPbSystErr     = "graphRpPbPHOSSystErr";
	TString NamePHOSRpPbStatErr     = "graphRpPbPHOSStatErr";

	 TFile*  fileRpPbPHOS = new TFile( fileNameRpPbPHOS.Data() );
       
	 graphRpPbPHOSSystErr = (TGraphAsymmErrors*)fileRpPbPHOS->Get(NamePHOSRpPbSystErr);
	 graphRpPbPHOSStatErr = (TGraphAsymmErrors*)fileRpPbPHOS->Get(NamePHOSRpPbStatErr);
	
	graphRpPbPHOSSystErr->RemovePoint(0);	
	graphRpPbPHOSStatErr->RemovePoint(0);	
	graphRpPbPHOSSystErr->RemovePoint(0);	
	graphRpPbPHOSStatErr->RemovePoint(0);

	//----------------Chose pPb Input-----------------------------------------
	TFile*  fileNeutralPionResultspPb;
	TDirectory* fNeutralPionResultspPbContainer;
	TH1F* histoInvYieldPi0pPb5023GeV;                   
	TH1F* histoInvYieldPi0pPb5023GeVSystErr;                   
        TGraphAsymmErrors* graphInvYieldPi0pPb5023GeVSystErr;  
	TGraphAsymmErrors* graphInvYieldPi0pPb5023GeVComplErr;
 	TGraphAsymmErrors* graphInvYieldPi0pPb5023GeVStatErr;
 	TGraphAsymmErrors* graphPHOSSystErrCancellation;  
 	TH1F*              histoPHOSSystErrCancellation;
  	TFile*             filePHOSSystErrCancellation;
	//pPb Input
	if (System.CompareTo("Comb")==0){

	  fileNeutralPionResultspPb = new TFile(fileNameNeutralPionCombResultspPb.Data());
        graphInvYieldPi0pPb5023GeVSystErr    = (TGraphAsymmErrors*)fileNeutralPionResultspPb->Get("CombinedPi0pPbSpectrumSysErr"); 
	graphInvYieldPi0pPb5023GeVComplErr   = (TGraphAsymmErrors*)fileNeutralPionResultspPb->Get("CombinedPi0pPbSpectrumTotErr");
 	graphInvYieldPi0pPb5023GeVStatErr    = (TGraphAsymmErrors*)fileNeutralPionResultspPb->Get("CombinedPi0pPbSpectrumStatErr");
	cout<<"Comb test"<<endl;
	graphInvYieldPi0pPb5023GeVStatErr->Print();
	graphInvYieldPi0pPb5023GeVSystErr->Print();
	graphInvYieldPi0pPb5023GeVComplErr->Print();
	
	//	graphInvYieldPi0pPb5023GeVStatErr->RemovePoint(0);

	cout<<" Comb"<<endl;
	graphInvYieldPi0pPb5023GeVStatErr->Print();
	graphInvYieldPi0pPb5023GeVSystErr->Print();
	graphInvYieldPi0pPb5023GeVComplErr->Print();

	}
	else if (System.CompareTo("PCM")==0){

	  fileNeutralPionResultspPb = new TFile(fileNameNeutralPionPCMResultspPb.Data());
	  fNeutralPionResultspPbContainer = (TDirectory*) fileNeutralPionResultspPb->GetDirectory("Pi0_pPb_5.023TeV_0-100%");
	if( ! fNeutralPionResultspPbContainer ) {cout<<"TList fNeutralPionResultspPbContainer does not exist: "<<endl; return;}
	
	histoInvYieldPi0pPb5023GeV              	   = (TH1F*)fNeutralPionResultspPbContainer->Get("CorrectedYieldPi0");
        graphInvYieldPi0pPb5023GeVSystErr    = (TGraphAsymmErrors*)fNeutralPionResultspPbContainer->Get("Pi0SystError"); 
	graphInvYieldPi0pPb5023GeVComplErr   = (TGraphAsymmErrors*)fNeutralPionResultspPbContainer->Get("Pi0ComplError");
 	graphInvYieldPi0pPb5023GeVStatErr    = new TGraphAsymmErrors(histoInvYieldPi0pPb5023GeV);
	cout<<"PCM test"<<endl;
	graphInvYieldPi0pPb5023GeVStatErr->Print();
	graphInvYieldPi0pPb5023GeVSystErr->Print();
	graphInvYieldPi0pPb5023GeVComplErr->Print();
	
	graphInvYieldPi0pPb5023GeVStatErr->RemovePoint(0);

	cout<<" test"<<endl;
	graphInvYieldPi0pPb5023GeVStatErr->Print();
	graphInvYieldPi0pPb5023GeVSystErr->Print();
	graphInvYieldPi0pPb5023GeVComplErr->Print();

	}else if(System.CompareTo("Dalitz")==0 ||System.CompareTo("DALITZ")==0 ) {

	  fileNeutralPionResultspPb = new TFile(fileNameNeutralPionDalitzResultspPb.Data());
	  
	  TDirectory* fNeutralPionResultspPbContainer = (TDirectory*) fileNeutralPionResultspPb->GetDirectory("Pi0_pPb_5.023TeV_0-100%");
	  
	  if( ! fNeutralPionResultspPbContainer ) {cout<<"TList fNeutralPionResultspPbContainer does not exist: "<<endl; return;}
	  
	  histoInvYieldPi0pPb5023GeV  = (TH1F*)fNeutralPionResultspPbContainer->Get("CorrectedYieldPi0");
	  graphInvYieldPi0pPb5023GeVComplErr  = (TGraphAsymmErrors*)fNeutralPionResultspPbContainer->Get("Pi0ComplError");
	  graphInvYieldPi0pPb5023GeVSystErr   = (TGraphAsymmErrors*)fNeutralPionResultspPbContainer->Get("Pi0SystError");
	  graphInvYieldPi0pPb5023GeVStatErr   = new TGraphAsymmErrors( histoInvYieldPi0pPb5023GeV );
	  graphInvYieldPi0pPb5023GeVStatErr->RemovePoint(0);
	  
	  if (! graphInvYieldPi0pPb5023GeVComplErr ){
	    
	    cout<<"graphInvYieldPi0pPb5023GeVComplErr"<<endl;
	    return ;
	  } 
	  
	  if ( !graphInvYieldPi0pPb5023GeVSystErr ){
	    
	    cout<<"graphInvYieldPi0pPb5023GeVSystErr"<<endl;
	    
	  return;
	  
	  }
	
	  if( !graphInvYieldPi0pPb5023GeVStatErr ) {
	    
	    cout<<"graphInvYieldPi0pPb5023GeVStatErr"<<endl;
	    return;
	  }
	  cout<<"Dalitz test"<<endl;
	    graphInvYieldPi0pPb5023GeVStatErr->Print();
	    graphInvYieldPi0pPb5023GeVSystErr->Print();
	    graphInvYieldPi0pPb5023GeVComplErr->Print();
	    

	} else if (System.CompareTo("EMCal")==0 ||System.CompareTo("EMCAL")==0 ){

	    fileNeutralPionResultspPb = new TFile(fileNameNeutralPionEMCalResultspPb.Data());
	    fNeutralPionResultspPbContainer = (TDirectory*) fileNeutralPionResultspPb->GetDirectory("Pi0_pPb_5.023TeV_0-100%");
	    if( ! fNeutralPionResultspPbContainer ) {cout<<"TList fNeutralPionResultspPbContainer does not exist: "<<endl; return;}
	    
	    histoInvYieldPi0pPb5023GeV              	   = (TH1F*)fNeutralPionResultspPbContainer->Get("CorrectedYieldPi0");
	    graphInvYieldPi0pPb5023GeVSystErr    = (TGraphAsymmErrors*)fNeutralPionResultspPbContainer->Get("Pi0SystError"); 
	    graphInvYieldPi0pPb5023GeVComplErr   = (TGraphAsymmErrors*)fNeutralPionResultspPbContainer->Get("Pi0ComplError");
	    graphInvYieldPi0pPb5023GeVStatErr    = new TGraphAsymmErrors(histoInvYieldPi0pPb5023GeV);
	    cout<<"EMCal test"<<endl;
	    graphInvYieldPi0pPb5023GeVStatErr->Print();
	    graphInvYieldPi0pPb5023GeVSystErr->Print();
	    graphInvYieldPi0pPb5023GeVComplErr->Print();
	    
	    graphInvYieldPi0pPb5023GeVStatErr->RemovePoint(0);
	    graphInvYieldPi0pPb5023GeVStatErr->RemovePoint(0);
	    graphInvYieldPi0pPb5023GeVStatErr->RemovePoint(0);
	    graphInvYieldPi0pPb5023GeVStatErr->RemovePoint(0);
	    graphInvYieldPi0pPb5023GeVStatErr->RemovePoint(0);
	    graphInvYieldPi0pPb5023GeVStatErr->RemovePoint(0);
	    graphInvYieldPi0pPb5023GeVStatErr->RemovePoint(0);
	    graphInvYieldPi0pPb5023GeVStatErr->RemovePoint(0);
	    graphInvYieldPi0pPb5023GeVStatErr->RemovePoint(0);
	    graphInvYieldPi0pPb5023GeVStatErr->RemovePoint(23); 
	    graphInvYieldPi0pPb5023GeVStatErr->RemovePoint(22); //remove last pT bin

	    graphInvYieldPi0pPb5023GeVSystErr->RemovePoint(24);
	    graphInvYieldPi0pPb5023GeVSystErr->RemovePoint(24);
	    graphInvYieldPi0pPb5023GeVSystErr->RemovePoint(24); 
	    graphInvYieldPi0pPb5023GeVSystErr->RemovePoint(23); 
	    graphInvYieldPi0pPb5023GeVSystErr->RemovePoint(22); //remove last pT bin

	    graphInvYieldPi0pPb5023GeVComplErr->RemovePoint(24);
	    graphInvYieldPi0pPb5023GeVComplErr->RemovePoint(24);
	    graphInvYieldPi0pPb5023GeVComplErr->RemovePoint(24);
	    graphInvYieldPi0pPb5023GeVComplErr->RemovePoint(23); 
	    graphInvYieldPi0pPb5023GeVComplErr->RemovePoint(22); //remove last pT bin
	    cout<<"EMCal test"<<endl;
	    graphInvYieldPi0pPb5023GeVStatErr->Print();
	    graphInvYieldPi0pPb5023GeVSystErr->Print();
	    graphInvYieldPi0pPb5023GeVComplErr->Print();
	} else if (System.CompareTo("PHOS")==0   ){	    
	  fileNeutralPionResultspPb   = new TFile(fileNameNeutralPionPHOSResultspPb.Data()); 
	  filePHOSSystErrCancellation = new TFile(fileNamePHOSSystErrCancellation.Data());
	  //	  graphInvYieldPi0pPb5023GeVComplErr= (TGraphAsymmErrors*) fileNeutralPionResultspPb->Get("gPi0Comp");
	  //	  graphInvYieldPi0pPb5023GeVStatErr= (TGraphAsymmErrors*) fileNeutralPionResultspPb->Get("hCor_stat");// fileNeutralPionResultspPb->Get("gPi0Stat");
	  //	  graphInvYieldPi0pPb5023GeVSystErr  = (TGraphAsymmErrors*) fileNeutralPionResultspPb->Get("hCor_syst");//fileNeutralPionResultspPb->Get("gPi0Syst");
	  histoInvYieldPi0pPb5023GeV  = (TH1F*)fileNeutralPionResultspPb->Get("hCor_stat");
	  graphInvYieldPi0pPb5023GeVStatErr   = new TGraphAsymmErrors( histoInvYieldPi0pPb5023GeV );
	  histoInvYieldPi0pPb5023GeVSystErr  = (TH1F*)fileNeutralPionResultspPb->Get("hCor_syst");
	  graphInvYieldPi0pPb5023GeVSystErr   = new TGraphAsymmErrors( histoInvYieldPi0pPb5023GeVSystErr );
	  histoPHOSSystErrCancellation  = (TH1F*)filePHOSSystErrCancellation->Get("hCancelSys");
	  graphPHOSSystErrCancellation   = new TGraphAsymmErrors( histoPHOSSystErrCancellation );
	    cout<<"PHOS test"<<endl;
	    graphInvYieldPi0pPb5023GeVStatErr->Print();
	    graphInvYieldPi0pPb5023GeVSystErr->Print();
	    //	    graphInvYieldPi0pPb5023GeVComplErr->Print();

	    //   graphInvYieldPi0pPb5023GeVStatErr->RemovePoint(0);	  
	    graphInvYieldPi0pPb5023GeVStatErr->RemovePoint(0);	  
	    graphInvYieldPi0pPb5023GeVStatErr->RemovePoint(25);	  
	    //   graphInvYieldPi0pPb5023GeVSystErr->RemovePoint(0);	  
	    graphInvYieldPi0pPb5023GeVSystErr->RemovePoint(0);	  
	    graphInvYieldPi0pPb5023GeVSystErr->RemovePoint(25);	  
	    // graphInvYieldPi0pPb5023GeVComplErr->RemovePoint(0);	  
	    // graphInvYieldPi0pPb5023GeVComplErr->RemovePoint(0);	 
	    
	    cout<<"PHOS test"<<endl;
	    graphInvYieldPi0pPb5023GeVComplErr    =  CalculateCombinedSysAndStatError( graphInvYieldPi0pPb5023GeVStatErr ,graphInvYieldPi0pPb5023GeVSystErr );
	    graphInvYieldPi0pPb5023GeVStatErr->Print();
	    graphInvYieldPi0pPb5023GeVSystErr->Print();
	    graphInvYieldPi0pPb5023GeVComplErr->Print();
	    cout<<"PHOS SystErrCancellation"<<endl;
	    graphPHOSSystErrCancellation->RemovePoint(0);	 
	    graphPHOSSystErrCancellation->RemovePoint(25);	 
	    graphPHOSSystErrCancellation->Print();
	} else {
	  
	  cout<<"Invalid System type:  "<<System.Data() <<endl;
	  return;
	  
	}


	//---------------------pp Input---------------------------------------------	
	
	if( resultsType.CompareTo("PCMPileUpCorrection") == 0 ){
	  
	       PileUpCorrection = kTRUE;
	       
	       cout<<"PileUp correction will be applied to PCM alone spetrum. "<<endl;
	       
	     	  
	}  
	
	if( thesisPlots.CompareTo("thesis") == 0 ) {
	  
	      thesis = kTRUE;
	      thesisPlotLabel = "This thesis";
	      
	} else {
	  
	  thesisPlotLabel = "ALICE";// work in progress";
	}
	  
	const char *fileNameEPS09sPi0AKK = "ExternalInputpPb/Theory/R_pi0_pPb_y0_eps09s/R_pPb5000_y0_eps09s_akk_mb.dat";
	const char *fileNameEPS09sPi0DSS = "ExternalInputpPb/Theory/R_pi0_pPb_y0_eps09s/R_pPb5000_y0_eps09s_fdss_mb.dat";
	const char *fileNameEPS09sPi0KKP = "ExternalInputpPb/Theory/R_pi0_pPb_y0_eps09s/R_pPb5000_y0_eps09s_kkp_mb.dat";
	const char *fileNameEPS09sPi0CGC = "ExternalInputpPb/Theory/ColorGlassCondensate.dat";
	
	
	
	
	
	
	TString fileNameRpPbChargePionsNew  	     	= "ExternalInputpPb/RpPb_502_PiKp_w_Kinks_12062015.root";
	TString histoNameChargedPionsRpPbSystErr     	= "hsys_RpPb_pion";
	TString histoNameChargedPionsRpPbStatErr     	= "hstat_RpPb_pion";
	TString histoNameChargedParticlesRpPbSystErr 	= "hsys_RpPb_charged";
	TString histoNameChargedParticlesRpPbStatErr 	= "hstat_RpPb_charged";
	TString histoNameChargedKaonsRpPbStatErr     	= "hstat_RpPb_kaon";
	TString histoNameChargedKaonsRpPbSystErr 	= "hsys_RpPb_kaon";
	TString histoNameChargedProtonsRpPbStatErr      = "hstat_RpPb_proton";
	TString histoNameChargedProtonsRpPbSystErr      = "hsys_RpPb_proton";
	
	
	 TFile*  fileRpPbChargedPionsNew = new TFile( fileNameRpPbChargePionsNew.Data() );
	
	
	
	//TH1D*   histo_pi_RpPb_SystErr     = (TH1D*)fileRpPbChargedPionsNew->Get(histoNameChargedPionsRpPbSystErr.Data());
	//TH1D*   histo_pi_RpPb_StatErr     = (TH1D*)fileRpPbChargedPionsNew->Get(histoNameChargedPionsRpPbStatErr.Data());
	TH1D*   histo_kaon_RpPb_StatErr   = (TH1D*)fileRpPbChargedPionsNew->Get(histoNameChargedKaonsRpPbStatErr.Data());
	TH1D*   histo_kaon_RpPb_SystErr   = (TH1D*)fileRpPbChargedPionsNew->Get(histoNameChargedKaonsRpPbSystErr.Data());
	TH1D*   histo_proton_RpPb_StatErr = (TH1D*)fileRpPbChargedPionsNew->Get(histoNameChargedProtonsRpPbStatErr.Data());
	TH1D*   histo_proton_RpPb_SystErr = (TH1D*)fileRpPbChargedPionsNew->Get(histoNameChargedProtonsRpPbSystErr.Data());
	
	TString fileNameRpPbChargePionsPaper  	     	= "ExternalInputpPb/ChargedPionRpPbPaper.root";
	TFile*  fileRpPbChargedPionsPaper = new TFile( fileNameRpPbChargePionsPaper.Data() );
	TH1D*   histo_pi_RpPb_SystErr     = (TH1D*)fileRpPbChargedPionsPaper->Get(histoNameChargedPionsRpPbSystErr.Data());
	TH1D*   histo_pi_RpPb_StatErr     = (TH1D*)fileRpPbChargedPionsPaper->Get(histoNameChargedPionsRpPbStatErr.Data());
	
	
	graphChargedPionRpPbSystErr   = new TGraphAsymmErrors(histo_pi_RpPb_SystErr);
	graphChargedPionRpPbStatErr   = new TGraphAsymmErrors(histo_pi_RpPb_StatErr);
	TGraphAsymmErrors*  graphChargedPionRpPbSystErrClone  = (TGraphAsymmErrors*) graphChargedPionRpPbSystErr->Clone();
	TGraphAsymmErrors*  graphChargedPionRpPbStatErrClone  = (TGraphAsymmErrors*) graphChargedPionRpPbStatErr->Clone();
	
	
	TGraphAsymmErrors*  graphChargedKaonRpPbStatErr   = new TGraphAsymmErrors(histo_kaon_RpPb_StatErr);
	TGraphAsymmErrors*  graphChargedKaonRpPbSystErr   = new TGraphAsymmErrors(histo_kaon_RpPb_SystErr);
	TGraphAsymmErrors*  graphChargedProtonRpPbStatErr = new TGraphAsymmErrors(histo_proton_RpPb_StatErr);
	TGraphAsymmErrors*  graphChargedProtonRpPbSystErr = new TGraphAsymmErrors(histo_proton_RpPb_SystErr);
	
	
	
	
	graphAsymmChargedParticlesRpPbSystErr = GetChargeParticlesRpPb2013("Syst");
	graphAsymmChargedParticlesRpPbStatErr = GetChargeParticlesRpPb2013("Stat");
	
	//GA pp Input	
        TFile*  fileNeutralPionCombResultsPP       = new TFile(fileNameNeutralPionCombResultsPP.Data());
	TFile*  fileNeutralPionPHOSResultsPP       = new TFile(fileNameNeutralPionPHOSResultsPP.Data());
	TFile*	fileNeutralPionEMCalResultsPP      = new TFile(fileNameNeutralPionEMCalResultsPP.Data());
	TFile*  fileNeutralPionEMCalResultsPP7TeV  = new TFile(fileNameNeutralPionEMCalResultsPP7TeV.Data());

	TString graphNameInvCrossSectionPi07TeVStatSystErr;
	TString graphNameInvCrossSectionPi07TeVStatErr;
	TString graphNameInvCrossSectionPi07TeVSysErr;
	
	TString graphNameInvCrossSectionPi02760GeVStatSystErr;
	TString graphNameInvCrossSectionPi02760GeVStatErr;
	TString graphNameInvCrossSectionPi02760GeVSysErr;
	
	TGraphAsymmErrors*  graphInvCrossSectionPi0Comb7TeVStatSystErr;
	TGraphAsymmErrors*  graphInvCrossSectionPi07TeVStatSystErr;
	TGraphAsymmErrors*  graphInvCrossSectionPi07TeVStatErr;
	TGraphAsymmErrors*  graphInvCrossSectionPi07TeVSystErr;
	
	TGraphAsymmErrors*  graphInvCrossSectionPi0Comb2760GeVStatSystErr;   	
	TGraphAsymmErrors*  graphInvCrossSectionPi02760GeVStatSystErr;   	
        TGraphAsymmErrors*  graphInvCrossSectionPi02760GeVStatErr;   	
	TGraphAsymmErrors*  graphInvCrossSectionPi02760GeVSystErr;   	
	

	TString fileNamePythia8 = "ExternalInputpPb/pythia8-5.02TeV-part-1.root";
	
	TString fitType  = "l";
	TString fitNameLabel = "Tsallis";
	
	
	
	
	
        TH1F* ppPythia = GetPPReferenceFromPythia(fileNamePythia8.Data());
	TGraphAsymmErrors*  graphppPythia = new TGraphAsymmErrors(ppPythia);
	
	
	graphppPythia->Print();
	

	
	
	
	if ( resultsType.CompareTo("Comb") == 0   ){
	  
	
	graphNameInvCrossSectionPi07TeVStatSystErr 	= Form("graphInvCrossSectionPi0%s7TeV",resultsType.Data()); 
	graphNameInvCrossSectionPi07TeVStatErr     	= Form("graphInvCrossSectionPi0%s7TeVStatErr",resultsType.Data());
	graphNameInvCrossSectionPi07TeVSysErr      	= Form("graphInvCrossSectionPi0%s7TeVSysErr",resultsType.Data());
	
	graphNameInvCrossSectionPi02760GeVStatSystErr 	= Form("graphInvCrossSectionPi0%s2760GeV",resultsType.Data());
	graphNameInvCrossSectionPi02760GeVStatErr 	= Form("graphInvCrossSectionPi0%s2760GeVStatErr",resultsType.Data());
	graphNameInvCrossSectionPi02760GeVSysErr  	= Form("graphInvCrossSectionPi0%s2760GeVSysErr",resultsType.Data());
	
	
	graphInvCrossSectionPi07TeVStatSystErr      	= (TGraphAsymmErrors*)fileNeutralPionCombResultsPP->Get(graphNameInvCrossSectionPi07TeVStatSystErr.Data());
	graphInvCrossSectionPi07TeVStatErr          	= (TGraphAsymmErrors*)fileNeutralPionCombResultsPP->Get(graphNameInvCrossSectionPi07TeVStatErr.Data());
	graphInvCrossSectionPi07TeVSystErr          	= (TGraphAsymmErrors*)fileNeutralPionCombResultsPP->Get(graphNameInvCrossSectionPi07TeVSysErr.Data());
	
	
	cout<<"************************************************Checking if they are shifted Combined**********************************************************"<<endl;
	graphInvCrossSectionPi07TeVStatErr->Print();
	
	cout<<"**************************************************"<<endl;
	
	graphInvCrossSectionPi07TeVSystErr->Print();
	
	
        graphInvCrossSectionPi02760GeVStatSystErr  	= (TGraphAsymmErrors*)fileNeutralPionCombResultsPP->Get(graphNameInvCrossSectionPi02760GeVStatSystErr.Data());	
        graphInvCrossSectionPi02760GeVStatErr      	= (TGraphAsymmErrors*)fileNeutralPionCombResultsPP->Get(graphNameInvCrossSectionPi02760GeVStatErr.Data());	
	graphInvCrossSectionPi02760GeVSystErr      	= (TGraphAsymmErrors*)fileNeutralPionCombResultsPP->Get(graphNameInvCrossSectionPi02760GeVSysErr.Data());	
	
	cout<<"Entro a Comb"<<endl;
	
	
	} else if ( resultsType.CompareTo("PHOS") == 0 ||   resultsType.CompareTo("PHOSPileUpCorrection") == 0){
	  
	
	graphNameInvCrossSectionPi07TeVStatErr     	= "graphInvCrossSectionPi0PHOSStat7TeV";
	if (resultsType.CompareTo("PHOSPileUpCorrection") == 0)   graphNameInvCrossSectionPi07TeVSysErr      	= "gInvCrossSection_Sys";// larger sys errors to account for plie-up
	else	graphNameInvCrossSectionPi07TeVSysErr      	= "graphInvCrossSectionPi0PHOSSys7TeV";//published PHOS sys errors
	
	graphNameInvCrossSectionPi02760GeVStatErr 	= "graphInvCrossSectionPi0PHOS2760GeVStatErr";
	graphNameInvCrossSectionPi02760GeVSysErr  	= "graphInvCrossSectionPi0PHOS2760GeVSysErr"; 
	
	graphInvCrossSectionPi07TeVStatErr          	= (TGraphAsymmErrors*)fileNeutralPionPHOSResultsPP->Get(graphNameInvCrossSectionPi07TeVStatErr.Data());
	graphInvCrossSectionPi07TeVSystErr          	= (TGraphAsymmErrors*)fileNeutralPionPHOSResultsPP->Get(graphNameInvCrossSectionPi07TeVSysErr.Data());
	
	
	cout<<"************************************************Checking if they are shifted Combined**********************************************************"<<endl;
	graphInvCrossSectionPi07TeVStatErr->Print();
	
	cout<<"**************************************************"<<endl;
	
	graphInvCrossSectionPi07TeVSystErr->Print();
	
	graphInvCrossSectionPi07TeVStatSystErr    =  CalculateCombinedSysAndStatError( graphInvCrossSectionPi07TeVStatErr , graphInvCrossSectionPi07TeVSystErr);

	graphInvCrossSectionPi02760GeVStatErr     = (TGraphAsymmErrors*)fileNeutralPionPHOSResultsPP->Get(graphNameInvCrossSectionPi02760GeVStatErr.Data());	
	graphInvCrossSectionPi02760GeVSystErr     = (TGraphAsymmErrors*)fileNeutralPionPHOSResultsPP->Get(graphNameInvCrossSectionPi02760GeVSysErr.Data());	
	graphInvCrossSectionPi02760GeVStatSystErr =  CalculateCombinedSysAndStatError( graphInvCrossSectionPi02760GeVStatErr , graphInvCrossSectionPi02760GeVSystErr);
	
	cout<<"************************************************Checking if they are shifted Combined**********************************************************"<<endl;
	graphInvCrossSectionPi02760GeVStatErr->Print();
	
	cout<<"**************************************************"<<endl;
	
	graphInvCrossSectionPi02760GeVSystErr->Print();
	
	cout<<"Entro a Comb"<<endl;
	
	
	} else if ( resultsType.CompareTo("EMCal") == 0 ||resultsType.CompareTo("EMCAL") == 0   ){
	  
	  cout << "hier!!!"<< endl;

	  graphNameInvCrossSectionPi07TeVStatErr     	= "pi0Stat";
	  graphNameInvCrossSectionPi07TeVSysErr      	= "pi0Syst";

	graphNameInvCrossSectionPi02760GeVStatErr 	= "graphInvCrossSectionPi0";
	graphNameInvCrossSectionPi02760GeVSysErr  	= "InvCrossSectionPi0Sys";
	
	  cout << "hier!!!"<< endl;
	  TH1D* histInvCrossSectionPi07TeVStatErr	= (TH1D*)fileNeutralPionEMCalResultsPP7TeV->Get(graphNameInvCrossSectionPi07TeVStatErr.Data());
	  TH1D* histInvCrossSectionPi07TeVSystErr	= (TH1D*)fileNeutralPionEMCalResultsPP7TeV->Get(graphNameInvCrossSectionPi07TeVSysErr.Data());
	  graphInvCrossSectionPi07TeVStatErr          	= new TGraphAsymmErrors(histInvCrossSectionPi07TeVStatErr);       graphInvCrossSectionPi07TeVStatErr->Print();
	graphInvCrossSectionPi07TeVSystErr          	= new TGraphAsymmErrors(histInvCrossSectionPi07TeVSystErr); 	  cout << "hier!!!"<< endl;graphInvCrossSectionPi07TeVStatErr->Print();
	graphInvCrossSectionPi07TeVStatSystErr    =  CalculateCombinedSysAndStatError( graphInvCrossSectionPi07TeVStatErr , graphInvCrossSectionPi07TeVSystErr);	  cout << "hier!!!"<< endl;	
	
	cout<<"************************************************Checking if they are shifted Combined**********************************************************"<<endl;
	graphInvCrossSectionPi07TeVStatErr->Print();
	
	cout<<"**************************************************"<<endl;
	
	graphInvCrossSectionPi07TeVSystErr->Print();
	
	TDirectory* fNeutralPionEMCalpPbContainer = (TDirectory*) fileNeutralPionEMCalResultsPP->GetDirectory("Pi02.76TeV");	
        graphInvCrossSectionPi02760GeVStatErr      	= (TGraphAsymmErrors*)fNeutralPionEMCalpPbContainer->Get(graphNameInvCrossSectionPi02760GeVStatErr.Data());	
	graphInvCrossSectionPi02760GeVSystErr      	= (TGraphAsymmErrors*)fNeutralPionEMCalpPbContainer->Get(graphNameInvCrossSectionPi02760GeVSysErr.Data());	
	graphInvCrossSectionPi02760GeVStatSystErr    	=  CalculateCombinedSysAndStatError( graphInvCrossSectionPi02760GeVStatErr ,graphInvCrossSectionPi02760GeVSystErr );	
	cout<<"Entro a Comb"<<endl;
	
	
	} else if(  resultsType.CompareTo("PCM") == 0 || resultsType.CompareTo("PCMPileUpCorrection") == 0 ){

	
	graphNameInvCrossSectionPi07TeVStatErr    	= "graphInvCrossSectionPi0PCMStat7TeV";
	graphNameInvCrossSectionPi07TeVSysErr     	= "graphInvCrossSectionPi0PCMSys7TeV";
	graphNameInvCrossSectionPi02760GeVStatErr 	= "graphInvCrossSectionPi0PCM2760GeVStatErr";
	graphNameInvCrossSectionPi02760GeVSysErr  	= "graphInvCrossSectionPi0PCM2760GeVSysErr";

	
	graphInvCrossSectionPi07TeVStatErr        = (TGraphAsymmErrors*)fileNeutralPionCombResultsPP->Get(graphNameInvCrossSectionPi07TeVStatErr.Data());
	graphInvCrossSectionPi07TeVSystErr        = (TGraphAsymmErrors*)fileNeutralPionCombResultsPP->Get(graphNameInvCrossSectionPi07TeVSysErr.Data());
	
	 cout<<"************************************************Checking if they are shifted **********************************************************"<<endl;
	 graphInvCrossSectionPi07TeVStatErr->Print();
	
	cout<<"**************************************************"<<endl;
	
	graphInvCrossSectionPi07TeVSystErr->Print();
	
	
	
	TGraphAsymmErrors* graphInvCrossSectionPi07TeVStatErrNoPileUpCorr = (TGraphAsymmErrors*)graphInvCrossSectionPi07TeVStatErr->Clone();
	
	if( PileUpCorrection == kTRUE ){
	  
	cout<<"entro a corregir"<<endl;  
	  
	TF1* fitCorrectionFactorsHistvsPtRatio 								= new TF1("fitCorrectionFactorsHistvsPtRatio","[0]/pow(x,[1])+[2]");
	fitCorrectionFactorsHistvsPtRatio->SetParameter(0,2.9737546081);
	fitCorrectionFactorsHistvsPtRatio->SetParameter(1,1.4795520406);
	fitCorrectionFactorsHistvsPtRatio->SetParameter(2,2.2652589579);
	
	Double_t* yStatNoPileUpCorr                             = graphInvCrossSectionPi07TeVStatErrNoPileUpCorr->GetY();

	Double_t* yStat 					= graphInvCrossSectionPi07TeVStatErr->GetY();
	Double_t* xStat 					= graphInvCrossSectionPi07TeVStatErr->GetX();
	Double_t* yStatErrUp 					= graphInvCrossSectionPi07TeVStatErr->GetEXhigh();
	Double_t* yStatErrDown 					= graphInvCrossSectionPi07TeVStatErr->GetEXlow();
	
	
	
	Double_t* ySyst 					= graphInvCrossSectionPi07TeVSystErr->GetY();
	Double_t* xSyst 					= graphInvCrossSectionPi07TeVSystErr->GetX();
	Double_t* ySystErrUp 					= graphInvCrossSectionPi07TeVSystErr->GetEXhigh();
	Double_t* ySystErrDown 					= graphInvCrossSectionPi07TeVSystErr->GetEXlow();
	
	
	const Int_t nBins = graphInvCrossSectionPi07TeVStatErr->GetN();
	
	Double_t ratioX[nBins];
	Double_t ratioY[nBins];
	Double_t errorX[nBins];
	Double_t errorY[nBins];
	
	
	cout<<"Applyin pileUp correction"<<endl;
	
	for (Int_t i = 0; i < nBins; i++){
	  
		cout << xStat[i] << "\t" << 100-fitCorrectionFactorsHistvsPtRatio->Eval(xStat[i]) << endl;
		cout << yStat[i] << "\t" << yStat[i]*(100-fitCorrectionFactorsHistvsPtRatio->Eval(xStat[i]))/100 << "\t" << yStat[i] << endl;
		yStat[i] 	   = yStat[i]*(100-fitCorrectionFactorsHistvsPtRatio->Eval(xStat[i]))/100;
		yStatErrUp[i]      = yStatErrUp[i]*(100-fitCorrectionFactorsHistvsPtRatio->Eval(xStat[i]))/100;
		yStatErrDown[i]    = yStatErrDown[i]*(100-fitCorrectionFactorsHistvsPtRatio->Eval(xStat[i]))/100;
		
		
		cout << xSyst[i] << "\t" << 100-fitCorrectionFactorsHistvsPtRatio->Eval(xSyst[i]) << endl;
		cout << ySyst[i] << "\t" << ySyst[i]*(100-fitCorrectionFactorsHistvsPtRatio->Eval(xSyst[i]))/100 << "\t" << ySyst[i] << endl;
		ySyst[i] 	   = ySyst[i]*(100-fitCorrectionFactorsHistvsPtRatio->Eval(xSyst[i]))/100;
		ySystErrUp[i]      = ySystErrUp[i]*(100-fitCorrectionFactorsHistvsPtRatio->Eval(xSyst[i]))/100;
		ySystErrDown[i]    = ySystErrDown[i]*(100-fitCorrectionFactorsHistvsPtRatio->Eval(xSyst[i]))/100;
		
		
		
		ratioX[i] = xSyst[i];
		errorX[i] = graphInvCrossSectionPi07TeVStatErr->GetErrorXhigh(i);
		ratioY[i] = yStat[i]/yStatNoPileUpCorr[i];
		errorY[i] = TMath::Sqrt(TMath::Power(graphInvCrossSectionPi07TeVStatErr->GetErrorYhigh(i)/ySyst[i],2) +TMath::Power(graphInvCrossSectionPi07TeVStatErrNoPileUpCorr->GetErrorYhigh(i)/yStatNoPileUpCorr[i],2))*ratioY[i];
		cout << "Ratio: " << ratioX[i] << "\t" <<  errorX[i] << "\t" << ratioY[i] << "\t" << errorY[i] << "\t" << errorY[i]/ratioY[i]*100 << "%"<<endl;

		
		
	}   
	cout << "*************************" << endl;
		 graphInvCrossSectionPi07TeVStatErr->Print();
	TGraphErrors* RatioNoPileUpAndPileUpCorr =  new TGraphErrors(nBins-1,ratioX,ratioY,errorX,errorY); 
	
	
	
	TCanvas* canvasRatioNoPileUpAndPileUpCorr = new TCanvas("canvasRatioNoPileUpAndPileUpCorr","Ratio between PileUp Correction And No PileUp correction spectra",200,10,1200,700);
	
	DrawGammaCanvasSettings( canvasRatioNoPileUpAndPileUpCorr,  0.1, 0.01, 0.015, 0.13);
	
	TH2F * histo2DRatioNoPileUpAndPileUpCorr = new TH2F("histo2DRatioNoPileUpAndPileUpCorr","histo2DRatioNoPileUpAndPileUpCorr",1000,0.,pTLimit,1000,0.3,2.5);
	SetStyleHistoTH2ForGraphs(histo2DRatioNoPileUpAndPileUpCorr, "#it{p}_{T} (GeV/#it{c})","PileUp Corr/No Corr", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 
	histo2DRatioNoPileUpAndPileUpCorr->GetYaxis()->SetRangeUser(0.5,1.5);
	histo2DRatioNoPileUpAndPileUpCorr->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphErr(RatioNoPileUpAndPileUpCorr,21,1.5, kBlack , kBlack);
	RatioNoPileUpAndPileUpCorr->Draw("E1psame");
	
	
		
	DrawGammaLines(0., 15.,1., 1.,2.,kBlack);
	
	TLatex *textRatioNoPileUpAndPileUpCorr = new TLatex(0.16,0.9,"pp minimum bias #sqrt{#it{s}} = 7 TeV");
	SetStyleTLatex( textRatioNoPileUpAndPileUpCorr,0.06,4); 
	textRatioNoPileUpAndPileUpCorr->Draw();
	
	
	TLegend* legendRatioNoPileUpAndPileUpCorr = new TLegend(0.18,0.15,0.9,0.21);
	legendRatioNoPileUpAndPileUpCorr->SetFillColor(0);
	legendRatioNoPileUpAndPileUpCorr->SetLineColor(0);
	legendRatioNoPileUpAndPileUpCorr->SetNColumns(2);
	legendRatioNoPileUpAndPileUpCorr->SetTextSize(0.045);
	legendRatioNoPileUpAndPileUpCorr->AddEntry(histo2DRatioNoPileUpAndPileUpCorr,"PileUp Correction / No PileUp Correction ","p");
	legendRatioNoPileUpAndPileUpCorr->Draw();
	
	
	canvasRatioNoPileUpAndPileUpCorr->SaveAs(Form("%s/Ratio_NoPileUpAndPileUpCorr.%s",outputDir.Data(),suffix.Data()));
	
	  
	  
	}
	
	
	graphInvCrossSectionPi07TeVStatSystErr    =  CalculateCombinedSysAndStatError( graphInvCrossSectionPi07TeVStatErr , graphInvCrossSectionPi07TeVSystErr);
	
	graphInvCrossSectionPi02760GeVStatErr     = (TGraphAsymmErrors*)fileNeutralPionCombResultsPP->Get(graphNameInvCrossSectionPi02760GeVStatErr.Data());	
	graphInvCrossSectionPi02760GeVSystErr     = (TGraphAsymmErrors*)fileNeutralPionCombResultsPP->Get(graphNameInvCrossSectionPi02760GeVSysErr.Data());	
	
	graphInvCrossSectionPi02760GeVStatSystErr =  CalculateCombinedSysAndStatError( graphInvCrossSectionPi02760GeVStatErr , graphInvCrossSectionPi02760GeVSystErr);
	
        cout<<"Entro a "<<endl;
	  
	  
	} else {
	  
	  cout<<"Invalid results type:  "<<resultsType.Data() <<endl;
	  return;
	  
	}
	  
     
	
	
	//**********************************************************************
	//Start R_pPb Analysis
	
	//*********************************************************************
	
	
	cout<<"LLego 1"<<endl;
	
	
	
	TGraphAsymmErrors* graphInvYieldPi0pPb5023GeVYShiftedSystErr    = (TGraphAsymmErrors*)graphInvYieldPi0pPb5023GeVSystErr->Clone();
 	TGraphAsymmErrors* graphInvYieldPi0pPb5023GeVYShiftedStatErr    = (TGraphAsymmErrors*)graphInvYieldPi0pPb5023GeVStatErr->Clone();
	TGraphAsymmErrors* graphInvYieldPi0pPb5023GeVYShiftedComplErr   = (TGraphAsymmErrors*)graphInvYieldPi0pPb5023GeVComplErr->Clone();
	
	
	cout<<"Llego 3"<<endl;
	
	

	

	
	cout<<"**********************************BinYShift ********************************************"<<endl;
	

	//Use Bylinkin fit from combined spectrum for YShift

	TF1 *fitBylinkinPi0pPb5023GeV = (TF1*)CommonFile->Get("FitCombinedPi0pPbSpectrum");	
	
	graphInvYieldPi0pPb5023GeVYShiftedComplErr  = ApplyYshiftIndividualSpectra(graphInvYieldPi0pPb5023GeVYShiftedComplErr ,fitBylinkinPi0pPb5023GeV);
	graphInvYieldPi0pPb5023GeVYShiftedSystErr   = ApplyYshiftIndividualSpectra( graphInvYieldPi0pPb5023GeVYShiftedSystErr,fitBylinkinPi0pPb5023GeV); 
	graphInvYieldPi0pPb5023GeVYShiftedStatErr   = ApplyYshiftIndividualSpectra( graphInvYieldPi0pPb5023GeVYShiftedStatErr, fitBylinkinPi0pPb5023GeV);
	

	//////////////////////////////////////////Convert to InvYield  7 TeV////////////////////////////////////////////
	
	
	TGraphAsymmErrors* graphInvYieldPi07TeVStatSystErr         = (TGraphAsymmErrors*) graphInvCrossSectionPi07TeVStatSystErr->Clone();
	TGraphAsymmErrors* graphInvYieldPi07TeVStatErr             = (TGraphAsymmErrors*) graphInvCrossSectionPi07TeVStatErr->Clone();
	TGraphAsymmErrors* graphInvYieldPi07TeVSystErr             = (TGraphAsymmErrors*) graphInvCrossSectionPi07TeVSystErr->Clone();

    

	///////////////////////////////////////////Convert to InvYield 2.76 TeV//////////////////////////////////////////
	
	TGraphAsymmErrors* graphInvYieldPi02760GeVStatSystErr      =  (TGraphAsymmErrors*)graphInvCrossSectionPi02760GeVStatSystErr->Clone();
	TGraphAsymmErrors* graphInvYieldPi02760GeVStatErr          =  (TGraphAsymmErrors*)graphInvCrossSectionPi02760GeVStatErr->Clone();
	TGraphAsymmErrors* graphInvYieldPi02760GeVSystErr          =  (TGraphAsymmErrors*)graphInvCrossSectionPi02760GeVSystErr->Clone();


	
	
	
	/////////////////////////////////////////////////Rebinning ////////////////////////////////////////////////////
	Double_t *Parameters2760GeV = NULL;
	Double_t *Parameters7TeV    = NULL;
	
	//This part obtain the "n" parameter from the fit of the combined spectrum///
	graphInvCrossSectionPi0Comb7TeVStatSystErr    = (TGraphAsymmErrors*)fileNeutralPionCombResultsPP->Get("graphInvCrossSectionPi0Comb7TeV");
	graphInvCrossSectionPi0Comb2760GeVStatSystErr = (TGraphAsymmErrors*)fileNeutralPionCombResultsPP->Get("graphInvCrossSectionPi0Comb2760GeV");
	
	
		
	
	if( FitFuncName.CompareTo("Tsallis") == 0 ){
	  
	  fitType = "l";
	  Parameters2760GeV = new Double_t[3];// = {2.4e+10,6.88,0.139}; //2.,5.,0.18
	  Parameters2760GeV[0] = 2.4e+10;
	  Parameters2760GeV[1] = 6.88;
	  Parameters2760GeV[2] = 0.139;
	  
	  Parameters7TeV    = new Double_t[3];// = {13.7e+09,7.,0.13};
	  Parameters7TeV[0] = 1.45613e+11;
	  Parameters7TeV[1] = 7;
	  Parameters7TeV[2] = 0.15;
 
          if( resultsType.CompareTo("PHOS") == 0 || resultsType.CompareTo("PHOSPileUpCorrection") == 0  ){  
	    Parameters7TeV[0] = 2.26842e+11;
	    Parameters7TeV[1] = 6.70106e+00;
	    Parameters7TeV[2] = 1.24671e-01;
	  }	  
	  fitNameLabel = "Tsallis";
	  
	  if( fixParam == 1  ){ //No applied for EMCal
	    
               Float_t minPt2760GeV = 0.4;
	       if (System.CompareTo("PHOS")== 0 ) minPt2760GeV = 0.8;
	       Float_t maxPt2760GeV = graphInvYieldPi02760GeVStatSystErr->GetXaxis()->GetBinUpEdge(graphInvYieldPi02760GeVStatSystErr->GetXaxis()->GetNbins());
	       cout<<"maxPt2760GeV: "<<maxPt2760GeV<<endl;
	       
	       TF1* FitToCombined2760GeV = FitObject("l","FitToCombined2760GeV","Pi0");
	       FitToCombined2760GeV->SetRange(minPt2760GeV,maxPt2760GeV);
 	       FitToCombined2760GeV->SetParameters(Parameters2760GeV[0],Parameters2760GeV[1],Parameters2760GeV[2]); 
	       graphInvCrossSectionPi0Comb2760GeVStatSystErr->Fit(FitToCombined2760GeV,"SQNRME+","",minPt2760GeV,maxPt2760GeV); 
	       
	       TString FitToCombined2760GeVParamToFile =  WriteParameterToFileLatexTable(FitToCombined2760GeV,kTRUE);
	       fileRebinSpectraFits << FitToCombined2760GeVParamToFile << endl;
	       
	       Parameters2760GeV[1] = FitToCombined2760GeV->GetParameter(1);

	       Float_t minPt7TeV = 0.3;
	       if (System.CompareTo("PHOS")== 0 ) minPt7TeV = 0.8;
	      
	       Float_t maxPt7TeV = graphInvYieldPi07TeVStatSystErr->GetXaxis()->GetBinUpEdge(graphInvYieldPi07TeVStatSystErr->GetXaxis()->GetNbins());
	       cout<<"maxPt7TeV: "<<maxPt7TeV<<endl;
	       TF1* FitToCombined7TeV = FitObject("l","FitToCombined7TeV","Pi0");
	      
	       FitToCombined7TeV->SetRange(minPt7TeV,maxPt7TeV);
 	       FitToCombined7TeV->SetParameters(Parameters7TeV[0],Parameters7TeV[1],Parameters7TeV[2]); 
	       graphInvCrossSectionPi0Comb7TeVStatSystErr->Fit(FitToCombined7TeV,"SQNRME+","",minPt7TeV,maxPt7TeV); 
	       TString FitToCombined7TeVParamToFile =  WriteParameterToFileLatexTable(FitToCombined7TeV,kTRUE);
	       fileRebinSpectraFits << FitToCombined7TeVParamToFile << endl;
	       Parameters7TeV[1] = FitToCombined7TeV->GetParameter(1);

	       
	             
	  }
	  

	  
	  
	} else if ( FitFuncName.CompareTo("Bylinkin") == 0){
	    
	  fitType = "tcm";
	  
	  
	  Parameters2760GeV = new Double_t[5];// = { 2.4e+10, 0.3, 1e+10,0.3,8};
	  // Parameters2760GeV[0] = 2.4e+10;
	  // Parameters2760GeV[1] = 0.3;
	  // Parameters2760GeV[2] = 1e+10;
	  // Parameters2760GeV[3] = 0.3;
	  // Parameters2760GeV[4] = 3.8;
	  // Parameters2760GeV[0] =  18.2676*xSection2760GeVppINEL* recalcBarn ;//7.4e+1;
	  // Parameters2760GeV[1] = 0.164972;
	  // Parameters2760GeV[2] =1.28029*xSection2760GeVppINEL* recalcBarn ;
	  // Parameters2760GeV[3] =0.702793;
	  // Parameters2760GeV[4] =3.16247 ;

	  Parameters2760GeV[0] =   4.31440e+09 ;//7.4e+1;
	  Parameters2760GeV[1] =  4.48647e-01 ;
	  Parameters2760GeV[2] =   1.31062e+11  ;
	  Parameters2760GeV[3]=    3.49106e-01 ;
	  Parameters2760GeV[4] =  2.78576e+00;	  
	  Parameters7TeV    = new Double_t[5];// = { 7.4e+10, 0.3, 1e+09,0.3,8};
	  // Parameters7TeV[0] =  18.2676*xSection7TeVINEL* recalcBarn ;
	  // Parameters7TeV[1] =  0.164972;
	  // Parameters7TeV[2] = 1.28029*xSection7TeVINEL* recalcBarn ;
	  // Parameters7TeV[3] =0.702793;
	  // Parameters7TeV[4] =3.16247 ;
	  Parameters7TeV[0] = 1.68474e+09   ;
	  Parameters7TeV[1] =5.57627e-01 ;
	  Parameters7TeV[2] = 1.51358e+11  ;
	  Parameters7TeV[3] = 3.98235e-01;
	  Parameters7TeV[4] =  2.80009e+00;

          if( resultsType.CompareTo("PHOS") == 0 || resultsType.CompareTo("PHOSPileUpCorrection") == 0   ){  
	    Parameters2760GeV[0] = 8.0e+11;
	    Parameters2760GeV[1] = 1.5e-01;
	    Parameters2760GeV[2] = 1.2e+10;
	    Parameters2760GeV[3] = 9.8e-01;
	    Parameters2760GeV[4] = 2.1e+00;
	  }
          if( resultsType.CompareTo("PHOS") == 0 || resultsType.CompareTo("PHOSPileUpCorrection") == 0 ||  resultsType.CompareTo("Comb") == 0 ){  
	    Parameters7TeV[0] = 8.04497e+11;
	    Parameters7TeV[1] = 1.57053e-01;
	    Parameters7TeV[2] = 1.19825e+10;
	    Parameters7TeV[3] = 9.85130e-01;
	    Parameters7TeV[4] = 2.13631e+00;
	  }	
	  
         fitNameLabel = "Bylinkin-Rostovtsev";
    
	}
	
	
	Float_t minPt = 0.4;
	if (System.CompareTo("EMCal")==0 ||System.CompareTo("EMCAL")==0 ) minPt=1.0; 
	if (System.CompareTo("PHOS")==0  ) minPt=0.8;
		Float_t maxPt = graphInvYieldPi02760GeVStatSystErr->GetXaxis()->GetBinUpEdge(graphInvYieldPi02760GeVStatSystErr->GetXaxis()->GetNbins());
		//Float_t maxPt = 8.;
		
	TGraphErrors* graphInvYieldPi02760GeVBinStatSystErr;
	TGraphErrors* graphInvYieldPi02760GeVBinStatErr;
	
	TGraphAsymmErrors* graphAInvYieldPi02760GeVBinStatSystErr;
	TGraphAsymmErrors* graphAInvYieldPi02760GeVBinSystErr;
	TGraphAsymmErrors* graphAInvYieldPi02760GeVBinStatErr;
	

	
	cout<<"Rebinning graphInvYieldPi02760GeVStatSystErr to pPb"<<endl;

	
	TF1* fitPi02760GeVBinStatSystErr    = RebinWithFitToTGraphWithMeanErr(graphInvYieldPi02760GeVStatSystErr,&graphAInvYieldPi02760GeVBinStatSystErr,graphInvYieldPi0pPb5023GeVYShiftedComplErr,0,fitType.Data(),minPt,maxPt,Parameters2760GeV,fixParam);
	fitPi02760GeVBinStatSystErr->SetName("Fit2760GeVRebin");	
	cout<<"Rebinning graphInvYieldPi02760GeVSystErr to pPb"<<endl;
	
	TString fitPi02760GeVBinStatSystErrParamToFile =  WriteParameterToFileLatexTable(fitPi02760GeVBinStatSystErr,kTRUE);
 	fileRebinSpectraFits << fitPi02760GeVBinStatSystErrParamToFile << endl;

		
	TF1* fitPi02760GeVBinSystErr        = RebinWithFitToTGraphWithMeanErr(graphInvYieldPi02760GeVSystErr,  &graphAInvYieldPi02760GeVBinSystErr,graphAInvYieldPi02760GeVBinStatSystErr,fitPi02760GeVBinStatSystErr,fitType.Data(),minPt,maxPt,Parameters2760GeV,fixParam);
	TF1* fitPi02760GeVBinStatErr	    = FillTGraphEYWithFitErr(graphInvYieldPi02760GeVStatErr,&graphAInvYieldPi02760GeVBinStatErr,graphAInvYieldPi02760GeVBinStatSystErr,fitType.Data(),minPt,maxPt,Parameters2760GeV);
	

	TGraphAsymmErrors* TGraphAymmmRatioToFitPi0PP2760GeV = (TGraphAsymmErrors*) CalculateGraphErrRatioToFit(graphInvYieldPi02760GeVStatSystErr,fitPi02760GeVBinStatSystErr);
	
	
	TCanvas* canvasRatioToFitPi0PP2760GeV = new TCanvas("canvasRatioToFitPi0PP2760GeV","Ratio to fit pp 2.76TeV",200,10,1200,700);
	
	DrawGammaCanvasSettings( canvasRatioToFitPi0PP2760GeV,  0.1, 0.01, 0.015, 0.13);
	
	TH2F * histo2DRatioToFitPi0PP2760GeV = new TH2F("histo2DRatioToFitPi0PP2760GeV","histo2DRatioToFitPi0PP2760GeV",1000,0.,pTLimit,1000,0.3,2.5);
	SetStyleHistoTH2ForGraphs(histo2DRatioToFitPi0PP2760GeV, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 
		histo2DRatioToFitPi0PP2760GeV->GetYaxis()->SetRangeUser(0.5,1.5);
	//	histo2DRatioToFitPi0PP2760GeV->GetYaxis()->SetRangeUser(0.4,1.6);
	histo2DRatioToFitPi0PP2760GeV->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(TGraphAymmmRatioToFitPi0PP2760GeV,20,1, kRed , kRed);
	TGraphAymmmRatioToFitPi0PP2760GeV->Draw("p,same,e");
		
	DrawGammaLines(0., 15.,1., 1.,2.,kBlack);
	
	TLatex *textRatioToFitPi0PP2760GeV = new TLatex(0.16,0.9,"pp minimum bias #sqrt{#it{s}} = 2.76 TeV");
	SetStyleTLatex( textRatioToFitPi0PP2760GeV,0.06,4); 
	textRatioToFitPi0PP2760GeV->Draw();
	
	TLatex *textRatioToFitPi0PP2760GeVFitFuncName = new TLatex(0.16,0.8,Form("%s Fit",fitNameLabel.Data()));
	SetStyleTLatex( textRatioToFitPi0PP2760GeVFitFuncName,0.06,4); 
	textRatioToFitPi0PP2760GeVFitFuncName->Draw();

	
	canvasRatioToFitPi0PP2760GeV->SaveAs(Form("%s/RatiotoFit_Pi0PP2760_%s_%s.%s",outputDir.Data(),resultsType.Data(),FitFuncName.Data(),suffix.Data()));
	
	
	
	
	
	
	
	minPt = 0.30;
	if (System.CompareTo("EMCal")==0 ||System.CompareTo("EMCAL")==0 ) minPt=0.60;
	if (System.CompareTo("PHOS")==0) minPt=0.8;
	 maxPt = graphInvYieldPi07TeVStatSystErr->GetXaxis()->GetBinUpEdge(graphInvYieldPi07TeVStatSystErr->GetXaxis()->GetNbins());
	
	 TGraphErrors*      graphInvYieldPi07TeVBinStatSystErr;
	 TGraphErrors*      graphInvYieldPi07TeVBinStatErr;

	 TGraphAsymmErrors* graphAInvYieldPi07TeVBinStatSystErr;
	 TGraphAsymmErrors* graphAInvYieldPi07TeVBinSystErr;
	 TGraphAsymmErrors* graphAInvYieldPi07TeVBinStatErr;
	
	
	 cout<<"Rebinning graphInvYieldPi07TeVStatSystErr to "<<endl;
	
	
	 TF1* fitPi07TeVBinStatSystErr    = RebinWithFitToTGraphWithMeanErr(graphInvYieldPi07TeVStatSystErr,&graphAInvYieldPi07TeVBinStatSystErr,   graphInvYieldPi0pPb5023GeVYShiftedComplErr,0,fitType.Data(),minPt,maxPt,Parameters7TeV,fixParam);
	 fitPi07TeVBinStatSystErr->SetName("Fit7TeVRebin");	
	 TString fitPi07TeVBinSystErrParamToFile =  WriteParameterToFileLatexTable(fitPi07TeVBinStatSystErr,kTRUE);
 	 fileRebinSpectraFits << fitPi07TeVBinSystErrParamToFile << endl;
	


	cout<<"Rebinning graphInvYieldPi07TeVSystErr to "<<endl;
	
	
	 TF1* fitPi07TeVBinSystErr        = RebinWithFitToTGraphWithMeanErr(graphInvYieldPi07TeVSystErr,&graphAInvYieldPi07TeVBinSystErr,   graphAInvYieldPi07TeVBinStatSystErr,fitPi07TeVBinStatSystErr,fitType.Data(),minPt,maxPt,Parameters7TeV,fixParam);
	 TF1* fitPi07TeVBinStatErr	  = FillTGraphEYWithFitErr(graphInvYieldPi07TeVStatErr,&graphAInvYieldPi07TeVBinStatErr,graphAInvYieldPi07TeVBinStatSystErr,fitType.Data(),minPt,maxPt,Parameters7TeV);
	
	 fitPi07TeVBinSystErr->Draw("");
	 graphInvYieldPi07TeVStatSystErr	->Draw("same");
	
	
	 cout<<"Start to computing"<<endl;
	graphInvYieldPi07TeVStatSystErr->Print();
	// TF1*fitLow7=FitObject("tcm1","fitInvCrossSectionPi07Low","Pi0");	
	// TF1*fitHigh7=FitObject("tcm2","fitInvCrossSectionPi07High","Pi0");	
	// fitLow7->SetParameters(fitPi07TeVBinStatSystErr->GetParameter(0),fitPi07TeVBinStatSystErr->GetParameter(1));
	// fitHigh7->SetParameters(fitPi07TeVBinStatSystErr->GetParameter(2),fitPi07TeVBinStatSystErr->GetParameter(3),fitPi07TeVBinStatSystErr->GetParameter(4));
	// fitLow7->SetRange(minPt,maxPt); 
	// fitHigh7->SetRange(minPt,maxPt);
	 TGraphAsymmErrors* TGraphAsymmRatioToFitPi0PP7TeV = (TGraphAsymmErrors*) CalculateGraphErrRatioToFit(graphInvYieldPi07TeVStatSystErr,fitPi07TeVBinStatSystErr);
	
	
	
	 TCanvas* canvasRatioToFitPi0PP7TeV = new TCanvas("canvasRatioToFitPi0PP7TeV","Ratio between RatioToFit pp 7TeV",200,10,1200,700);
	
	 DrawGammaCanvasSettings( canvasRatioToFitPi0PP7TeV,  0.1, 0.01, 0.015, 0.13);
	
	 TH2F * histo2DRatioToFitPi0PP7TeV = new TH2F("histo2DRatioToFitPi0PP7TeV","histo2DRatioToFitPi0PP7TeV",1000,0.,pTLimit,1000,0.3,2.5);
	SetStyleHistoTH2ForGraphs(histo2DRatioToFitPi0PP7TeV, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 
	histo2DRatioToFitPi0PP7TeV->GetYaxis()->SetRangeUser(0.5,1.5);
	//	histo2DRatioToFitPi0PP7TeV->GetYaxis()->SetRangeUser(0.4,1.6);
	histo2DRatioToFitPi0PP7TeV->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(TGraphAsymmRatioToFitPi0PP7TeV,20,1, kRed , kRed);
	//HistoRatioToFitPi0PP7TeV->Draw("E1psame");
	//HistoRatioToFitPi0PP7TeV->Draw("same");
	TGraphAsymmRatioToFitPi0PP7TeV->Draw("p,same,e");
	
	
		
	DrawGammaLines(0., 15.,1., 1.,2.,kBlack);
	
	TLatex *textRatioToFitPi0PP7TeV = new TLatex(0.16,0.9,"pp minimum bias #sqrt{#it{s}} = 7 TeV");
	SetStyleTLatex( textRatioToFitPi0PP7TeV,0.06,4); 
	textRatioToFitPi0PP7TeV->Draw();
	
	TLatex *textRatioToFitPi0PP7TeVFitFuncName = new TLatex(0.16,0.8,Form("%s Fit",fitNameLabel.Data()));
	SetStyleTLatex( textRatioToFitPi0PP7TeVFitFuncName,0.06,4); 
	textRatioToFitPi0PP7TeVFitFuncName->Draw();
	
	
	
	
	canvasRatioToFitPi0PP7TeV->SaveAs(Form("%s/RatiotoFit_Pi0PP7TeV_%s_%s.%s",outputDir.Data(),resultsType.Data(),FitFuncName.Data(),suffix.Data()));
	
	
	
	///////////////////////////////////////////
	cout<<"Just for test Checking Statistic errors"<<endl;
	 
	TGraphAsymmErrors* graphInvYieldPi07TeVStatErrClone      = ProduceTGraphAsymmToPlotErrors(graphInvYieldPi07TeVStatErr);
	TGraphAsymmErrors* graphAInvYieldPi07TeVBinStatErrClone  = ProduceTGraphAsymmToPlotErrors(graphAInvYieldPi07TeVBinStatErr);
	
	 
        TCanvas* canvasCrossCheckStatisticalErrPP7TeV = new TCanvas("canvasCrossCheckStatisticalErrPP7TeV","Cross-check statistical errors of pp 7 TeV",200,10,1200,700);
	 
	
	DrawGammaCanvasSettings( canvasCrossCheckStatisticalErrPP7TeV,  0.1, 0.01, 0.015, 0.13);
	canvasCrossCheckStatisticalErrPP7TeV->SetLogx();
	TH2F * histo2DStatisticTestPi0PP7TeV = new TH2F("histo2DStatisticTestPi0PP7TeV","histo2DStatisticTestPi0PP7TeV",1000,0.3,pTLimit,1000,0.0,50);
	SetStyleHistoTH2ForGraphs(histo2DStatisticTestPi0PP7TeV, "#it{p}_{T} (GeV/#it{c})","Statistical error %", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 
	histo2DStatisticTestPi0PP7TeV->GetYaxis()->SetRangeUser(0.001,30.0);
	histo2DStatisticTestPi0PP7TeV->GetXaxis()->SetRangeUser(0.3,20.0);
	histo2DStatisticTestPi0PP7TeV->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(    graphInvYieldPi07TeVStatErrClone,20,1, kRed , kRed,1,kTRUE,0);
	graphInvYieldPi07TeVStatErrClone->Draw("p,same,e");
	
	DrawGammaSetMarkerTGraphAsym(graphAInvYieldPi07TeVBinStatErrClone,20,1, kBlue,kBlue,1,kTRUE,0);
	graphAInvYieldPi07TeVBinStatErrClone->Draw("p,same,e");
	
	TLegend* legendStatisticalErrPP7TeV = new TLegend(0.2,0.75,0.6,0.95);
	legendStatisticalErrPP7TeV->SetFillColor(0);
	legendStatisticalErrPP7TeV->SetLineColor(0);
	legendStatisticalErrPP7TeV->SetNColumns(1);
	legendStatisticalErrPP7TeV->SetTextSize(0.03);
	legendStatisticalErrPP7TeV->AddEntry(graphInvYieldPi07TeVStatErrClone,"Measured pp #sqrt{s} = 7 TeV","pef");
	legendStatisticalErrPP7TeV->AddEntry(graphAInvYieldPi07TeVBinStatErrClone,"Calculated pp #sqrt{s} = 7 TeV","pef");
	legendStatisticalErrPP7TeV->Draw();
	
	canvasCrossCheckStatisticalErrPP7TeV->SaveAs(Form("%s/StatisticalErr_Pi0PP7TeV_%s_%s.%s",outputDir.Data(),resultsType.Data(),FitFuncName.Data(),suffix.Data()));
	
	
	
	TGraphAsymmErrors* graphInvYieldPi02760GeVStatErrClone      = ProduceTGraphAsymmToPlotErrors(graphInvYieldPi02760GeVStatErr);
	TGraphAsymmErrors* graphAInvYieldPi02760GeVBinStatErrClone  = ProduceTGraphAsymmToPlotErrors(graphAInvYieldPi02760GeVBinStatErr);
	
	 
        TCanvas* canvasCrossCheckStatisticalErrPP2760GeV = new TCanvas("canvasCrossCheckStatisticalErrPP2760GeV","Cross-check statistical errors of pp 7 TeV",200,10,1200,700);
	 
	
	DrawGammaCanvasSettings( canvasCrossCheckStatisticalErrPP2760GeV,  0.1, 0.01, 0.015, 0.13);
	canvasCrossCheckStatisticalErrPP2760GeV->SetLogx();
	TH2F * histo2DStatisticTestPi0PP2760GeV = new TH2F("histo2DStatisticTestPi0PP2760GeV","histo2DStatisticTestPi0PP2760GeV",1000,0.3,pTLimit,1000,0.0,50);
	SetStyleHistoTH2ForGraphs(histo2DStatisticTestPi0PP2760GeV, "#it{p}_{T} (GeV/#it{c})","Statistical error %", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 
	histo2DStatisticTestPi0PP2760GeV->GetYaxis()->SetRangeUser(0.001,30.0);
	histo2DStatisticTestPi0PP2760GeV->GetXaxis()->SetRangeUser(0.3,20.0);
	histo2DStatisticTestPi0PP2760GeV->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphInvYieldPi02760GeVStatErrClone,20,1, kRed , kRed,1,kTRUE,0);
	graphInvYieldPi02760GeVStatErrClone->Draw("p,same,e");
	
	DrawGammaSetMarkerTGraphAsym(graphAInvYieldPi02760GeVBinStatErrClone,20,1, kBlue,kBlue,1,kTRUE,0);
	
	graphAInvYieldPi02760GeVBinStatErrClone->Draw("p,same,e");
	
	
	TLegend* legendStatisticalErrPP2760GeV = new TLegend(0.2,0.75,0.6,0.95);
	legendStatisticalErrPP2760GeV->SetFillColor(0);
	legendStatisticalErrPP2760GeV->SetLineColor(0);
	legendStatisticalErrPP2760GeV->SetNColumns(1);
	legendStatisticalErrPP2760GeV->SetTextSize(0.03);
	legendStatisticalErrPP2760GeV->AddEntry(graphInvYieldPi02760GeVStatErrClone,"Measured pp #sqrt{s} = 2.76 TeV","pef");
	legendStatisticalErrPP2760GeV->AddEntry(graphAInvYieldPi02760GeVBinStatErrClone,"Calculated pp #sqrt{s} = 2.76 TeV","pef");
	
	legendStatisticalErrPP2760GeV->Draw();
	
	
	
	canvasCrossCheckStatisticalErrPP2760GeV->SaveAs(Form("%s/StatisticalErr_Pi0PP2760GeV_%s_%s.%s",outputDir.Data(),resultsType.Data(),FitFuncName.Data(),suffix.Data()));
	
	//////////////////////////////////////////////////////////////////////////
	
	
	graphInvYieldPi02760GeVBinStatSystErr    	= ConvertTGraphAsymmErrorstoTGraphErrors(graphAInvYieldPi02760GeVBinStatSystErr);
	graphInvYieldPi07TeVBinStatSystErr    		= ConvertTGraphAsymmErrorstoTGraphErrors(graphAInvYieldPi07TeVBinStatSystErr);

	//Compute interporlation with statistics only to get the statistics errors  
	graphInvYieldPi02760GeVBinStatErr    		= ConvertTGraphAsymmErrorstoTGraphErrors(graphAInvYieldPi02760GeVBinStatErr);
	graphInvYieldPi07TeVBinStatErr    		= ConvertTGraphAsymmErrorstoTGraphErrors(graphAInvYieldPi07TeVBinStatErr);



	TGraphErrors*  graphErrosInterPolation5023GeVBin    	   = GetInterpolSpectrum2D(graphInvYieldPi02760GeVBinStatSystErr,    graphInvYieldPi07TeVBinStatSystErr,    2760,7000,5023,"SystStat");
	TGraphErrors*  graphErrosInterPolation5023GeVBinStatErr    = GetInterpolSpectrum2D(graphInvYieldPi02760GeVBinStatErr,        graphInvYieldPi07TeVBinStatErr,        2760,7000,5023,"Stat");
	
	

	TString	  namePtvsSqrtsPlot = 	Form("%s/Pt_vs_Sqrts_%s.%s",outputDir.Data(),System.Data(),suffix.Data());
	
	if (System.CompareTo("Dalitz")==0 || System.CompareTo("DALITZ")==0){
	                  PlotInterpolationPtBins(graphPtvsSqrts,gPtvsEnergiesSystem,fPowerlawSystem,graphErrosInterPolation5023GeVBin,5,4,namePtvsSqrtsPlot);
	}	else	  PlotInterpolationPtBins(graphPtvsSqrts,gPtvsEnergiesSystem,fPowerlawSystem,graphErrosInterPolation5023GeVBin,7,5,namePtvsSqrtsPlot);

	// if (System.CompareTo("PCM")==0){
	//   PlotAlphavsPt(graphAlpha,    "PCM",    thesisPlotLabel.Data(),  Form("%s/Alpha_vs_Pt_%s.%s",   outputDir.Data(),System.Data(),suffix.Data()));
	// } else if (System.CompareTo("Dalitz")==0 || System.CompareTo("DALITZ")==0){	
	//   PlotAlphavsPt(graphAlpha,    "Dalitz",    thesisPlotLabel.Data(),  Form("%s/Alpha_vs_Pt_%s.%s",   outputDir.Data(),System.Data(),suffix.Data()));
	// }else
	PlotAlphavsPt(graphAlpha, System.Data()  ,    thesisPlotLabel.Data(),  Form("%s/Alpha_vs_Pt_%s.%s",   outputDir.Data(),System.Data(),suffix.Data()));
	
	
	
	/////////////////////////////////////////////////////////////////////
	
	
		
	cout<<"Statistics for PCM"<<endl;  //The statistical errors are taken from the interpolated pp reference with statistics
	
	Int_t     nPoints    =   graphErrosInterPolation5023GeVBin->GetN();
	
	Double_t *xValue     =   graphErrosInterPolation5023GeVBin->GetX();
	Double_t *xStatErr   =   graphInvYieldPi07TeVBinStatSystErr->GetEX(); //X errors are taken from the rebbined 7 spectrum
	
	Double_t *yValue     =   graphErrosInterPolation5023GeVBin->GetY();
	Double_t *yStatErr   =   graphErrosInterPolation5023GeVBinStatErr->GetEY(); //Taken with statistics
	
	
	graphAErrosInterPolation5023GeVBinStatErr = new TGraphAsymmErrors(nPoints,xValue,yValue,xStatErr,xStatErr,yStatErr,yStatErr);
	
	
	
	
	cout<<"Systematics for PCM"<<endl;
	graphAErrosInterPolation5023GeVBinSystErr    = CalculateSystErrors(graphErrosInterPolation5023GeVBin,graphAInvYieldPi02760GeVBinSystErr,graphAInvYieldPi07TeVBinSystErr);
	
	
	cout<<"Binning Interpolation    "<<endl;
	
	graphAErrosInterPolation5023GeVBinStatErr->Print();
	

	
	
	
	if( resultsType.CompareTo("PCM") == 0 || resultsType.CompareTo("PCMPileUpCorrection") == 0 ) {
	  if( System.CompareTo("Dalitz") == 0 || System.CompareTo("DALITZ") == 0 ) { 
	    graphErrosInterPolation5023GeVBinSystWOMatErr     = CancelOutMaterialError(graphAErrosInterPolation5023GeVBinSystErr, "PCMDalitz",graphPHOSSystErrCancellation);
	    graphInvYieldPi0pPb5023GeVYShiftedSystWOMatErr    = CancelOutMaterialError(graphInvYieldPi0pPb5023GeVYShiftedSystErr,"Dalitz",graphPHOSSystErrCancellation);
	  }else{
	    graphErrosInterPolation5023GeVBinSystWOMatErr     = CancelOutMaterialError(graphAErrosInterPolation5023GeVBinSystErr, "PCMPCM",graphPHOSSystErrCancellation);
	    graphInvYieldPi0pPb5023GeVYShiftedSystWOMatErr    = CancelOutMaterialError(graphInvYieldPi0pPb5023GeVYShiftedSystErr,"PCMPCM",graphPHOSSystErrCancellation);
	  }
	
	  CalcRpPb(graphErrosInterPolation5023GeVBinSystWOMatErr,graphAErrosInterPolation5023GeVBinStatErr,graphInvYieldPi0pPb5023GeVYShiftedSystWOMatErr, graphInvYieldPi0pPb5023GeVYShiftedStatErr,&graphRpPbSystErr,&graphRpPbStatErr);
	
	} else if (( resultsType.CompareTo("PHOS") == 0 || resultsType.CompareTo("PHOSPileUpCorrection") == 0 ) && System.CompareTo("PHOS") == 0 ){
	  graphErrosInterPolation5023GeVBinSystWOMatErr     = CancelOutMaterialError(graphAErrosInterPolation5023GeVBinSystErr, "PHOSPHOS",graphPHOSSystErrCancellation);
	    graphInvYieldPi0pPb5023GeVYShiftedSystWOMatErr    = CancelOutMaterialError(graphInvYieldPi0pPb5023GeVYShiftedSystErr,"PHOSPHOS",graphPHOSSystErrCancellation);

	    CalcRpPb(graphErrosInterPolation5023GeVBinSystWOMatErr,graphAErrosInterPolation5023GeVBinStatErr,graphInvYieldPi0pPb5023GeVYShiftedSystWOMatErr, graphInvYieldPi0pPb5023GeVYShiftedStatErr,&graphRpPbSystErr,&graphRpPbStatErr);
	
	}else if(System.CompareTo("Comb") == 0 ) {
	  
	  CalcRpPb(graphAErrosInterPolation5023GeVBinSystErr,graphAErrosInterPolation5023GeVBinStatErr,graphInvYieldPi0pPb5023GeVYShiftedSystErr, graphInvYieldPi0pPb5023GeVYShiftedStatErr,&graphRpPbSystErr,&graphRpPbStatErr,kFALSE);
	
	} else{
	  CalcRpPb(graphAErrosInterPolation5023GeVBinSystErr,graphAErrosInterPolation5023GeVBinStatErr,graphInvYieldPi0pPb5023GeVYShiftedSystErr, graphInvYieldPi0pPb5023GeVYShiftedStatErr,&graphRpPbSystErr,&graphRpPbStatErr);

	}
	
	///////////////////////////////////Quality plots //////////////////////////
	
	
	TGraphErrors* graphInvYieldPi0pPb5023GeVStatOnlyErrUpDown;
	TGraphErrors* graphInterPolation5023GeVBinStatOnlyErrUpDown;
	TGraphErrors* graphAInvYieldPi07TeVBinOnlyStatErrUpDown;
	TGraphErrors* graphAInvYieldPi02760GeVBinOnlyStatErrUpDown;
	
	TGraphErrors* graphInvYieldPi0pPb5023GeVSystOnlyErrUpDown;
	TGraphErrors* graphInterPolation5023GeVBinSystOnlyErrUpDown;
	TGraphErrors* graphAInvYieldPi07TeVBinOnlySystErrUpDown;
	TGraphErrors* graphAInvYieldPi02760GeVBinOnlySystErrUpDown;

	
	GetTGraphErrorsUpDown(graphInvYieldPi0pPb5023GeVYShiftedStatErr,    &graphInvYieldPi0pPb5023GeVStatOnlyErrUpDown);
	GetTGraphErrorsUpDown(graphAErrosInterPolation5023GeVBinStatErr,    &graphInterPolation5023GeVBinStatOnlyErrUpDown);
	GetTGraphErrorsUpDown(graphAInvYieldPi07TeVBinStatErr,              &graphAInvYieldPi07TeVBinOnlyStatErrUpDown);
        GetTGraphErrorsUpDown(graphAInvYieldPi02760GeVBinStatErr,           &graphAInvYieldPi02760GeVBinOnlyStatErrUpDown);
	
	GetTGraphErrorsUpDown(graphInvYieldPi0pPb5023GeVYShiftedSystErr,    &graphInvYieldPi0pPb5023GeVSystOnlyErrUpDown);
	GetTGraphErrorsUpDown(graphAErrosInterPolation5023GeVBinSystErr,    &graphInterPolation5023GeVBinSystOnlyErrUpDown);
	GetTGraphErrorsUpDown(graphAInvYieldPi07TeVBinSystErr,              &graphAInvYieldPi07TeVBinOnlySystErrUpDown);
        GetTGraphErrorsUpDown(graphAInvYieldPi02760GeVBinSystErr,           &graphAInvYieldPi02760GeVBinOnlySystErrUpDown);
	
	PlotErrors(graphInterPolation5023GeVBinStatOnlyErrUpDown,  graphInterPolation5023GeVBinSystOnlyErrUpDown,
		   graphAInvYieldPi07TeVBinOnlyStatErrUpDown,      graphAInvYieldPi07TeVBinOnlySystErrUpDown,
	           graphAInvYieldPi02760GeVBinOnlyStatErrUpDown,   graphAInvYieldPi02760GeVBinOnlySystErrUpDown,
		   System.Data(), thesisPlotLabel.Data(),outputDir.Data(),suffix.Data());
	

	
	/////////////Drawing InvYield pp@2.76 GeV//////////////
	

	
	TCanvas* cInvYield2760GeVBin = new TCanvas("cInvYield2760GeVBin","Invariant Yield pp@2.76 TeV",200,10,700,550);
    	DrawGammaCanvasSettings( cInvYield2760GeVBin,  0.15, 0.02, 0.03, 0.1);
	cInvYield2760GeVBin->SetLogx();
	cInvYield2760GeVBin->SetLogy();
	
	TPad* padComparisonInvYield2760GeVBin = new TPad("padComparisonInvYield2760GeVBin", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padComparisonInvYield2760GeVBin, 0.15, 0.02, 0.03, 0.1);
	padComparisonInvYield2760GeVBin->Draw();
	
	
	TH2F * histoInvYield2760GeVBin;
	histoInvYield2760GeVBin = new TH2F("histoInvYield2760GeVBin","histoInvYield2760GeVBin",1000,pTLimitLow,pTLimit,1000,2e0,10e12);
	SetStyleHistoTH2ForGraphs(histoInvYield2760GeVBin, "#it{p}_{T} (GeV/#it{c})", "#it{E}#frac{d^{2}#sigma}{d#it{p}^{3}}(pb GeV^{-2} #it{c}^{3})", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoInvYield2760GeVBin->DrawCopy(); 
	
        graphInvYieldPi02760GeVStatSystErr->SetFillColor(0);
	graphAInvYieldPi02760GeVBinStatSystErr->SetFillColor(0);
	
	
	DrawGammaSetMarkerTGraphAsym(graphAInvYieldPi02760GeVBinStatSystErr, 20,1,kBlue,kBlue);
	DrawGammaSetMarkerTGraphAsym(graphInvYieldPi02760GeVStatSystErr, 20,1,kRed,kRed);
	
	TLegend* legendInvYields2760GeVBin = new TLegend(0.2,0.15,0.6,0.35);
	legendInvYields2760GeVBin->SetFillColor(0);
	legendInvYields2760GeVBin->SetLineColor(0);
	legendInvYields2760GeVBin->SetNColumns(1);
	legendInvYields2760GeVBin->SetTextSize(0.03);
	legendInvYields2760GeVBin->AddEntry(graphInvYieldPi02760GeVStatSystErr,"pp at #sqrt{#it{s}} = 2.76 TeV","pef");
	legendInvYields2760GeVBin->AddEntry(graphAInvYieldPi02760GeVBinStatSystErr,"pp at #sqrt{#it{s}} = 2.76 TeV rebinned","pef");
	legendInvYields2760GeVBin->AddEntry(fitPi02760GeVBinStatSystErr,Form("%s fit",fitNameLabel.Data()),"l");
	
	graphAInvYieldPi02760GeVBinStatSystErr->Draw("p,same,e1");
	graphInvYieldPi02760GeVStatSystErr->Draw("p,same,e1");
	fitPi02760GeVBinStatSystErr->Draw("p,same");
	legendInvYields2760GeVBin->Draw("same");
	// fitLow->SetLineColor(kGreen);	
	// fitLow->Draw("p,same");	
	// fitHigh->SetLineColor(kRed);	
	// fitHigh->Draw("p,same");
	
	
	cInvYield2760GeVBin->Update();
	cInvYield2760GeVBin->SaveAs(Form("%s/InvYield_PP_Bin%s_2760GeV.%s",outputDir.Data(),System.Data(),suffix.Data()));
	
	
	
	
	TCanvas* cInvYield7TeVBin = new TCanvas("cInvYield7TeVBin","Invariant Yield pp@2.76 TeV",200,10,700,550);
    	DrawGammaCanvasSettings( cInvYield7TeVBin,  0.15, 0.02, 0.03, 0.1);
	cInvYield7TeVBin->SetLogx();
	cInvYield7TeVBin->SetLogy();
	
	TPad* padComparisonInvYield7TeVBin = new TPad("padComparisonInvYield7TeVBin", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padComparisonInvYield7TeVBin, 0.15, 0.02, 0.03, 0.1);
	padComparisonInvYield7TeVBin->Draw();
	
	
	TH2F * histoInvYield7TeVBin;
	histoInvYield7TeVBin = new TH2F("histoInvYield7TeVBin","histoInvYield7TeVBin",1000,pTLimitLow,pTLimit,1000,2e0,10e12);
	SetStyleHistoTH2ForGraphs(histoInvYield7TeVBin, "#it{p}_{T} (GeV/#it{c})", "#it{E}#frac{d^{2}#sigma}{d#it{p}^{3}}(pb GeV^{-2} #it{c}^{3})", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoInvYield7TeVBin->DrawCopy(); 
	
        graphInvYieldPi07TeVStatSystErr->SetFillColor(0);
	graphAInvYieldPi07TeVBinStatSystErr->SetFillColor(0);
	
	
	
	DrawGammaSetMarkerTGraphAsym(graphAInvYieldPi07TeVBinStatSystErr, 20,1,kBlue,kBlue);
	DrawGammaSetMarkerTGraphAsym(graphInvYieldPi07TeVStatSystErr, 20,1,kRed,kRed);
	
	
	TLegend* legendInvYields7TeVBin = new TLegend(0.2,0.15,0.6,0.35);
	legendInvYields7TeVBin->SetFillColor(0);
	legendInvYields7TeVBin->SetLineColor(0);
	legendInvYields7TeVBin->SetNColumns(1);
	legendInvYields7TeVBin->SetTextSize(0.03);
	legendInvYields7TeVBin->AddEntry(graphInvYieldPi07TeVStatSystErr,"pp at #sqrt{#it{s}} = 7 TeV","pef");
	legendInvYields7TeVBin->AddEntry(graphAInvYieldPi07TeVBinStatSystErr,"pp at #sqrt{#it{s}} = 7 TeV rebinned","pef");
	legendInvYields7TeVBin->AddEntry(fitPi07TeVBinStatSystErr,Form("%s fit",fitNameLabel.Data()),"l");
	
	
	graphAInvYieldPi07TeVBinStatSystErr->Draw("p,same,e1");
	graphInvYieldPi07TeVStatSystErr->Draw("p,same,e1");
	fitPi07TeVBinStatSystErr->Draw("same");
	legendInvYields7TeVBin->Draw("same");
	// fitLow7->SetLineColor(kGreen);
	// fitLow7->Draw("same");	
	// fitHigh7->SetLineColor(kRed);
	// fitHigh7->Draw("same");	

	cInvYield7TeVBin->Update();
	cInvYield7TeVBin->SaveAs(Form("%s/InvYield_PP_Bin%s_7TeV.%s",outputDir.Data(),System.Data(),suffix.Data()));
	
	
	
	
	////////////////////////////Scaling graphsErrors to be painted////////////////////
	
	
	
	TGraphErrors* graphInvYieldPi07TeVBinStatSystErrScaled       = new TGraphErrors( *graphInvYieldPi07TeVBinStatSystErr );
	TGraphErrors* graphInvYieldPi02760GeVBinScaled  	        = new TGraphErrors( *graphInvYieldPi02760GeVBinStatSystErr);
	TGraphErrors* graphErrosInterPolation5023GeVBinScaled        = new TGraphErrors( *graphErrosInterPolation5023GeVBin );
	
	TGraphErrors* graphErrosInterPolation5023GeVBinTpPb                  = new TGraphErrors( *graphErrosInterPolation5023GeVBin );
	
	
	graphInvYieldPi07TeVBinStatSystErrScaled        = ScaleGraph(graphInvYieldPi07TeVBinStatSystErrScaled,4);
	graphErrosInterPolation5023GeVBinScaled         = ScaleGraph(graphErrosInterPolation5023GeVBinScaled,2);
	graphInvYieldPi02760GeVBinScaled                = ScaleGraph(graphInvYieldPi02760GeVBinScaled,1);
	graphErrosInterPolation5023GeVBinTpPb           = ScaleGraph(graphErrosInterPolation5023GeVBinTpPb,fTpPb);
	
	
	
	
	TCanvas* cppAndpPbComarison = new TCanvas("cppAndpPbComarison","Comparison Invariant Yields",550,700);
		
	DrawGammaCanvasSettings( cppAndpPbComarison,  0.18, 0.02, 0.03, 0.1);
	cppAndpPbComarison->SetLogx();
	cppAndpPbComarison->SetLogy();
	
	TPad* padppAndpPbComarison = new TPad("padppAndpPbComarison", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padppAndpPbComarison, 0.24, 0.02, 0.03, 0.1);
	padppAndpPbComarison->Draw();
	
	
	TH2F * histoppAndpPbComparison;
	histoppAndpPbComparison = new TH2F("histoppAndpPbComparison","histoppAndpPbComparison",1000,pTLimitLow,pTLimit,1000,1e-10,1e2);
	SetStyleHistoTH2ForGraphs(histoppAndpPbComparison, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (c/GeV)^{2}", 0.032,0.04, 0.04,0.04, 1,1.9);
	//	histoppAndpPbComparison->GetXaxis()->SetRangeUser(0.6,10);
	//	histoppAndpPbComparison->GetYaxis()->SetRangeUser(1e-6,5);
	histoppAndpPbComparison->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0pPb5023GeVYShiftedComplErr,20,1,  kBlack,   kBlack);
	
	
	graphInvYieldPi0pPb5023GeVYShiftedComplErr->Draw("p,e1");
	graphErrosInterPolation5023GeVBinTpPb->Draw("XL,same");
	
	TLatex *LabelppAndpPb = new TLatex(0.45,0.9,thesisPlotLabel.Data());
	SetStyleTLatex( LabelppAndpPb, 0.04,4);
	LabelppAndpPb->Draw();
	TLatex *labelSpectraPi0pPb;
	labelSpectraPi0pPb= new TLatex(0.45,0.85,SystemLabel.Data());
	SetStyleTLatex( labelSpectraPi0pPb, 0.04,4);
	labelSpectraPi0pPb->Draw();


	graphInvYieldPi0pPb5023GeVYShiftedComplErr->SetFillColor(0);
	graphErrosInterPolation5023GeVBinTpPb->SetFillColor(0);
	TLegend* legendppAndpPbComarison = new TLegend(0.2,0.15,0.6,0.25);
	legendppAndpPbComarison->SetFillColor(0);
	legendppAndpPbComarison->SetLineColor(0);
	legendppAndpPbComarison->SetNColumns(1);
	legendppAndpPbComarison->SetTextSize(0.04);
	legendppAndpPbComarison->AddEntry(graphInvYieldPi0pPb5023GeVYShiftedComplErr,Form("%s p-Pb, MB, #sqrt{#it{s}_{NN}} = 5.02 TeV",System.Data()));
	legendppAndpPbComarison->AddEntry(graphErrosInterPolation5023GeVBinTpPb,"pp reference x #LT #it{T}_{pPb} #GT");
	
	legendppAndpPbComarison->Draw();

		
	
	cout<<"TGraphErrors pp@5.023 TeV"<<endl;
	graphErrosInterPolation5023GeVBinScaled->Print();
	
	cppAndpPbComarison->Update();
	cppAndpPbComarison->SaveAs(Form("%s/InvYield_Comparison_PP_pPb_%s.%s",outputDir.Data(),System.Data(),suffix.Data()));
	
	
		
	
	
	
	
	TCanvas* cInvYieldComparison = new TCanvas("cInvYieldComparison","Comparison Invariant Yields",200,10,700,500);
		
	DrawGammaCanvasSettings( cInvYieldComparison,  0.15, 0.02, 0.03, 0.1);
	cInvYieldComparison->SetLogx();
	cInvYieldComparison->SetLogy();
	
	TPad* padInvYieldComparison = new TPad("padInvYieldComparison", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padInvYieldComparison, 0.15, 0.02, 0.03, 0.1);
	padInvYieldComparison->Draw();
	
	
	TH2F * histoInvYieldComparison;
	histoInvYieldComparison = new TH2F("histoInvYieldComparison","histoInvYieldComparison",1000,pTLimitLow,pTLimit,1000,2e2,10e11);
	SetStyleHistoTH2ForGraphs(histoInvYieldComparison, "#it{p}_{T} (GeV/#it{c})", "#it{E}#frac{d^{2}#sigma}{d#it{p}^{3}}(pb GeV^{-2} #it{c}^{3})", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoInvYieldComparison->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphErr(graphInvYieldPi07TeVBinStatSystErrScaled,     	20, 1, kRed+1,   kRed+1);
	DrawGammaSetMarkerTGraphErr(graphInvYieldPi02760GeVBinScaled, 	20, 1, kGreen+2, kGreen+2);
	DrawGammaSetMarkerTGraphErr(graphErrosInterPolation5023GeVBinScaled,  20, 1, kBlue+1,  kBlue+1);
	
	
	
	
	
	graphInvYieldPi07TeVBinStatSystErrScaled->Draw("p,same,e1");
	graphInvYieldPi02760GeVBinScaled->Draw("p,same,e1");
	graphErrosInterPolation5023GeVBinScaled->Draw("p,same,e1");
	
	
	TLatex *LabelpPb = new TLatex(0.65,0.9,"ALICE");// work in progress");
	SetStyleTLatex( LabelpPb, 0.03,4);
	LabelpPb->Draw();
	TLatex *labelSpectraPi0LabelpPb2;
	labelSpectraPi0LabelpPb2= new TLatex(0.65,0.85,SystemLabel.Data());
       	SetStyleTLatex( labelSpectraPi0LabelpPb2, 0.03,4);
	labelSpectraPi0LabelpPb2->Draw();

	TLatex *labelSpectraPi0LabelpPb;
	if (System.CompareTo("PHOS")==0)	labelSpectraPi0LabelpPb= new TLatex(0.65,0.8,"|#it{y}_{#pi^{0},lab}| < 0.13");
	else if (System.CompareTo("EMCal")==0)	labelSpectraPi0LabelpPb= new TLatex(0.65,0.8,"|#it{y}_{#pi^{0},lab}| < 0.5");
	else	labelSpectraPi0LabelpPb= new TLatex(0.65,0.8,"|#it{y}_{#pi^{0},lab}| < 0.8");
	SetStyleTLatex( labelSpectraPi0LabelpPb, 0.03,4);
	labelSpectraPi0LabelpPb->Draw();  

	graphInvYieldPi07TeVBinStatSystErrScaled->SetFillColor(0);
	graphInvYieldPi02760GeVBinScaled->SetFillColor(0);
	graphErrosInterPolation5023GeVBinScaled->SetFillColor(0);
	TLegend* legendInvYields = new TLegend(0.2,0.15,0.6,0.35);
	legendInvYields->SetFillColor(0);
	legendInvYields->SetLineColor(0);
	legendInvYields->SetNColumns(1);
	legendInvYields->SetTextSize(0.03);
	legendInvYields->AddEntry(graphInvYieldPi07TeVBinStatSystErrScaled,"pp@7 TeV (x4)","pef");
	legendInvYields->AddEntry(graphErrosInterPolation5023GeVBinScaled,"pp@5.02 TeV - Interpolation (x2)","pef");
	legendInvYields->AddEntry(graphInvYieldPi02760GeVBinScaled,"pp@2.76 TeV (x1)","pef");
	
	
	legendInvYields->Draw();

		
	
	cInvYieldComparison->Update();
	cInvYieldComparison->SaveAs(Form("%s/InvYield_Comparison_PP_%s_2760GeV_5023GeV_7TeV.%s",outputDir.Data(),System.Data(),suffix.Data()));
	
	//cout<<"Check where the something strange is"<<endl; 
	
	
	TGraphErrors* graphRebinnedAStat = NULL; 
	TGraphErrors* graphRebinnedASyst = NULL; 
	TGraphErrors* graphRebinnedBStat = NULL; 
	TGraphErrors* graphRebinnedBSyst = NULL;
	
	TGraphErrors* graphRatioPPReference_Pythia    = CalculateRatioBetweenSpectraWithDifferentBinning(graphppPythia,graphppPythia, graphAErrosInterPolation5023GeVBinStatErr, graphAErrosInterPolation5023GeVBinSystErr,       kTRUE,  kTRUE,&graphRebinnedAStat,&graphRebinnedASyst,&graphRebinnedBStat,&graphRebinnedBSyst);
	
	//cout<<"Check where the something strange is"<<endl; 
	
  
	
	TCanvas* cInvYieldComparison_Pyhtia = new TCanvas("cInvYieldComparison_Pyhtia","Comparison Invariant Yields",200,10,700,500);
		
	DrawGammaCanvasSettings( cInvYieldComparison_Pyhtia,  0.15, 0.02, 0.03, 0.1);
	cInvYieldComparison_Pyhtia->SetLogx();
	cInvYieldComparison_Pyhtia->SetLogy();
	
	TPad* padInvYieldComparison_Pyhtia = new TPad("padInvYieldComparison_Pyhtia", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padInvYieldComparison_Pyhtia, 0.15, 0.02, 0.03, 0.1);
	padInvYieldComparison_Pyhtia->Draw();
	
	
	
	TH2F * histoInvYieldComparison_Pyhtia;
	histoInvYieldComparison_Pyhtia = new TH2F("histoInvYieldComparison_Pyhtia","histoInvYieldComparison_Pyhtia",1000,pTLimitLow,pTLimit,1000,2e0,10e11);
	SetStyleHistoTH2ForGraphs(histoInvYieldComparison_Pyhtia, "#it{p}_{T} (GeV/#it{c})", "#it{E}#frac{d^{2}#sigma}{d#it{p}^{3}}(pb GeV^{-2} #it{c}^{3})", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoInvYieldComparison_Pyhtia->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphppPythia,20, 1, kRed+1,   kRed+1);
	DrawGammaSetMarkerTGraphErr(graphErrosInterPolation5023GeVBin,  20, 1, kBlue+1,  kBlue+1);

	graphErrosInterPolation5023GeVBin->Draw("p,same,e1");
	graphppPythia->Draw("p,same,e1");
	
	
	TLatex *LabelpPb_Pyhtia = new TLatex(0.65,0.9,"ALICE");// work in progress");
	SetStyleTLatex( LabelpPb_Pyhtia, 0.03,4);
	LabelpPb_Pyhtia->Draw();
	TLatex *labelSpectraPi0LabelpPb2_Pythia;
	labelSpectraPi0LabelpPb2_Pythia= new TLatex(0.65,0.85,SystemLabel.Data());
	SetStyleTLatex( labelSpectraPi0LabelpPb2_Pythia, 0.03,4);
	labelSpectraPi0LabelpPb2_Pythia->Draw();

	TLatex *labelSpectraPi0LabelpPb_Pythia = new TLatex(0.65,0.8,"|#it{y}_{#pi^{0},lab}| < 0.8");
	SetStyleTLatex( labelSpectraPi0LabelpPb_Pythia, 0.03,4);
	labelSpectraPi0LabelpPb_Pythia->Draw();  

	graphErrosInterPolation5023GeVBin->SetFillColor(0);
	TLegend* legendInvYields_Pythia = new TLegend(0.2,0.15,0.6,0.35);
	legendInvYields_Pythia->SetFillColor(0);
	legendInvYields_Pythia->SetLineColor(0);
	legendInvYields_Pythia->SetNColumns(1);
	legendInvYields_Pythia->SetTextSize(0.03);
	legendInvYields_Pythia->AddEntry(graphErrosInterPolation5023GeVBin,"pp@5.02 TeV - Interpolation","pef");
	legendInvYields_Pythia->AddEntry(graphppPythia,"pp@5.02 TeV  Pythia","pef");	
	legendInvYields_Pythia->Draw();
	cInvYieldComparison_Pyhtia->Update();
	cInvYieldComparison_Pyhtia->SaveAs(Form("%s/InvYield_Comparison_PP_%s_Pythia.%s",outputDir.Data(),System.Data(),suffix.Data()));
	
	
	
	
	TCanvas* canvasRatio_Pythia = new TCanvas("canvasRatio_Pythia","PP reference  compared to Pythia",200,10,1200,700);
	
	DrawGammaCanvasSettings( canvasRatio_Pythia,  0.1, 0.01, 0.015, 0.13);
	
	TH2F * histo2DRatio_Pythia = new TH2F("histo2DRatio_Pythia","histo2DRatio_Pythia",1000,pTLimitLow,pTLimit,1000,0.3,2.5);
	SetStyleHistoTH2ForGraphs(histo2DRatio_Pythia, "#it{p}_{T} (GeV/#it{c})","#it{R}_{pPb}", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 
	histo2DRatio_Pythia->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphErr(graphRatioPPReference_Pythia,21,1.5, kBlack , kBlack);
	graphRatioPPReference_Pythia->Draw("E1psame");
	
	
		
	DrawGammaLines(0., 15.,1., 1.,2.,kBlack);
	
	TLatex *textRatio_Pythia = new TLatex(0.16,0.9,"p-Pb minimum bias #sqrt{#it{s}_{NN}} = 5.02 TeV");
	SetStyleTLatex( textRatio_Pythia,0.06,4); 
	textRatio_Pythia->Draw();
	
	
	TLegend* legendRatio_Pythia = new TLegend(0.18,0.15,0.9,0.21);
	legendRatio_Pythia->SetFillColor(0);
	legendRatio_Pythia->SetLineColor(0);
	legendRatio_Pythia->SetNColumns(2);
	legendRatio_Pythia->SetTextSize(0.045);
	legendRatio_Pythia->AddEntry(graphRatioPPReference_Pythia,Form("pp reference Pythia / pp reference data (%s)",System.Data()),"p");
	legendRatio_Pythia->Draw();
	
	
	canvasRatio_Pythia->SaveAs(Form("%s/Ratio_%s_Pythia.%s",outputDir.Data(),System.Data(),suffix.Data()));
	

	//////////////////////////////////////////
	
	
	
	TCanvas* cInvYieldPi0pPb5023GeVYShifted = new TCanvas("cInvYieldPi0pPb5023GeVYShifted","Pi0  pPb at 5.23 TeV X shifted",200,10,700,500);
		
	DrawGammaCanvasSettings( cInvYieldPi0pPb5023GeVYShifted,  0.15, 0.02, 0.03, 0.1);
	cInvYieldPi0pPb5023GeVYShifted->SetLogx();
	cInvYieldPi0pPb5023GeVYShifted->SetLogy();
	
	TPad* padInvYieldPi0pPb5023GeVYShifted = new TPad("padInvYieldPi0pPb5023GeVYShifted", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padInvYieldPi0pPb5023GeVYShifted, 0.15, 0.02, 0.03, 0.1);
	padInvYieldPi0pPb5023GeVYShifted->Draw();
	
	
	TH2F * histoInvYieldPi0pPb5023GeVYShifted;
	histoInvYieldPi0pPb5023GeVYShifted = new TH2F("histoInvYieldPi0pPb5023GeVYShifted","histoInvYieldPi0pPb5023GeVYShifted",1000,pTLimitLow,pTLimit,1000,2e-8,5e1);
	SetStyleHistoTH2ForGraphs(histoInvYieldPi0pPb5023GeVYShifted, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (c/GeV)^{2}", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoInvYieldPi0pPb5023GeVYShifted->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0pPb5023GeVYShiftedComplErr, 	20, 1, kRed,   kRed);
	DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0pPb5023GeVComplErr,  	        20, 1, kBlue,  kBlue);
	
	
	graphInvYieldPi0pPb5023GeVYShiftedComplErr->SetFillColor(0);
	graphInvYieldPi0pPb5023GeVComplErr->SetFillColor(0);
	
	
	graphInvYieldPi0pPb5023GeVComplErr->Draw("p,e1");
	graphInvYieldPi0pPb5023GeVYShiftedComplErr->Draw("p,same,e1");
	
	TLegend* legendInvYieldPi0pPb5023GeVYShifted = new TLegend(0.35,0.80,0.90,0.90);
	legendInvYieldPi0pPb5023GeVYShifted->SetFillColor(0);
	legendInvYieldPi0pPb5023GeVYShifted->SetLineColor(0);
	legendInvYieldPi0pPb5023GeVYShifted->SetTextSize(0.03);
	legendInvYieldPi0pPb5023GeVYShifted->AddEntry(graphInvYieldPi0pPb5023GeVComplErr,Form("#pi^{0} p-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV (%s)",System.Data()));
	legendInvYieldPi0pPb5023GeVYShifted->AddEntry(graphInvYieldPi0pPb5023GeVYShiftedComplErr,Form("#pi^{0} p-Pb at #sqrt{#it{s}_{NN}} = 5.02 #it{y}-shifted (%s)",System.Data()));
	legendInvYieldPi0pPb5023GeVYShifted->Draw("same");
	
		
	
	cInvYieldPi0pPb5023GeVYShifted->Update();
	cInvYieldPi0pPb5023GeVYShifted->SaveAs(Form("%s/InvYield_YShifted_%s_pPb_5023GeV.%s",outputDir.Data(),System.Data(),suffix.Data()));
	
	
	TCanvas* canvasRpPb = new TCanvas("canvasRpPb","R_{pPb} ",200,10,1000,700);
	
	DrawGammaCanvasSettings( canvasRpPb,  0.1, 0.01, 0.015, 0.11);
	TH2F * histo2DRpPb = new TH2F("histo2DRpPb","histo2DRpPb",1000,0.,pTLimit,1000,0.0,5.0);
	SetStyleHistoTH2ForGraphs(histo2DRpPb,pTLabel.Data(),RpPbLabel.Data(), 0.04,0.04, 0.04,0.04, 1.2,1.0, 512, 512); 
	
	
	//	histo2DRpPb->GetXaxis()->SetRangeUser(0.0,14.5);
	//	if(System.CompareTo("EMCal")==0){
	  histo2DRpPb->GetYaxis()->SetRangeUser(0.39,2.01);
	// } else 	  histo2DRpPb->GetYaxis()->SetRangeUser(0.19,1.81);
	// histo2DRpPb->GetYaxis()->CenterTitle(kTRUE);
	
	histo2DRpPb->Draw(); 

	if(System.CompareTo("Dalitz")==0 ||System.CompareTo("DALITZ")==0 ){
	  DrawGammaSetMarkerTGraphAsym(graphRpPbSystErr,20,1.5, colorsArray[kDalitz], colorsArray[kDalitz], 1, kTRUE);
	  DrawGammaSetMarkerTGraphAsym(graphRpPbStatErr,20,1.5, colorsArray[kDalitz], colorsArray[kDalitz]);
	  
	}else if(System.CompareTo("EMCal")==0){
	  DrawGammaSetMarkerTGraphAsym(graphRpPbSystErr,20,1.5, colorsArray[kEMCal], colorsArray[kEMCal], 1, kTRUE);
	  DrawGammaSetMarkerTGraphAsym(graphRpPbStatErr,20,1.5, colorsArray[kEMCal],   colorsArray[kEMCal]);
	}else if(System.CompareTo("PHOS")==0){
	  DrawGammaSetMarkerTGraphAsym(graphRpPbSystErr,20,1.5, colorsArray[kPHOS], colorsArray[kPHOS], 1, kTRUE);
	  DrawGammaSetMarkerTGraphAsym(graphRpPbStatErr,20,1.5, colorsArray[kPHOS],   colorsArray[kPHOS]);
	}else{
	  DrawGammaSetMarkerTGraphAsym(graphRpPbSystErr,20,1.5, colorsArray[kPCM], colorsArray[kPCM], 1, kTRUE);
	  DrawGammaSetMarkerTGraphAsym(graphRpPbStatErr,20,1.5, colorsArray[kPCM],   colorsArray[kPCM]);
	}


	graphRpPbSystErr->Draw("2same");
 	graphRpPbStatErr->Draw("p,same,e");
	graphRpPbStatErr->SetFillColor(0);
	
	
	
	
	DrawGammaLines(0., pTLimit,1., 1.,1.,kBlack,2);
	
	TLatex *labelpPb = new TLatex(0.15,0.92,thesisPlotLabel.Data());
	SetStyleTLatex( labelpPb, 0.04,4);
	labelpPb->Draw();
	
	TLatex *labelpPbEnergy =  new TLatex(0.15,0.87, energyLabel.Data());
	SetStyleTLatex( labelpPbEnergy, 0.04,4);
	labelpPbEnergy->Draw();
	
	
	TLegend* legendRpPb = new TLegend(0.15,0.79,0.50,0.84);
	legendRpPb->SetFillColor(0);
	legendRpPb->SetLineColor(0);
	legendRpPb->SetTextFont(42);
	legendRpPb->SetNColumns(3);
	legendRpPb->SetTextSize(0.04);
	if (System.CompareTo("PHOS")==0)	legendRpPb->AddEntry(graphRpPbStatErr,Form("%s, |#it{y}| < 0.13",SystemLabel.Data()),"p");
	else if (System.CompareTo("EMCal")==0)	legendRpPb->AddEntry(graphRpPbStatErr,Form("%s, |#it{y}| < 0.5",SystemLabel.Data()),"p");
	else	legendRpPb->AddEntry(graphRpPbStatErr,Form("%s, |#it{y}| < 0.8",SystemLabel.Data()),"p");
	legendRpPb->Draw();
	canvasRpPb->Update();
	
	

	canvasRpPb->Update(); 
	
	canvasRpPb->SaveAs(Form("%s/RpPb_%s.%s",outputDir.Data(),System.Data(),suffix.Data()));
	//--------------------Comparison to old PHOS result--------------------
	TCanvas* canvasRpPbPHOS = new TCanvas("canvasRpPbPHOS","R_{pPb} ",200,10,1000,700);
	
	DrawGammaCanvasSettings( canvasRpPbPHOS,  0.1, 0.01, 0.015, 0.11);
	TH2F * histo2DRpPbPHOS = new TH2F("histo2DRpPbPHOS","histo2DRpPbPHOS",1000,0.,pTLimit,1000,0.0,5.0);
	SetStyleHistoTH2ForGraphs(histo2DRpPbPHOS,pTLabel.Data(),RpPbLabel.Data(), 0.04,0.04, 0.04,0.04, 1.2,1.0, 512, 512); 
	
	
	//	histo2DRpPbPHOS->GetXaxis()->SetRangeUser(0.0,14.5);
	histo2DRpPbPHOS->GetYaxis()->SetRangeUser(0.19,1.81);
	histo2DRpPbPHOS->GetYaxis()->CenterTitle(kTRUE);
	
	histo2DRpPbPHOS->Draw(); 
	
	graphRpPbSystErr->Draw("2same");
 	graphRpPbStatErr->Draw("p,same,e");

	 DrawGammaSetMarkerTGraphAsym(graphRpPbPHOSSystErr,24,1.5, kRed+2,kRed+2,1.0,kTRUE);
	 DrawGammaSetMarkerTGraphAsym(graphRpPbPHOSStatErr,24,1.5, kRed+2 ,kRed+2);	
	graphRpPbPHOSSystErr->Draw("2same");
	graphRpPbPHOSStatErr->Draw("p,same,e");	
	
	DrawGammaLines(0., 14.5,1., 1.,1.,kBlack,2);
	
	TLatex *labelpPbPHOS = new TLatex(0.15,0.92,thesisPlotLabel.Data());
	SetStyleTLatex( labelpPbPHOS, 0.04,4);
	labelpPbPHOS->Draw();
	
	TLatex *labelpPbEnergyPHOS =  new TLatex(0.15,0.87, energyLabel.Data());
	SetStyleTLatex( labelpPbEnergyPHOS, 0.04,4);
	labelpPbEnergyPHOS->Draw();
	
	
	TLegend* legendRpPbPHOS = new TLegend(0.15,0.75,0.50,0.84);
	legendRpPbPHOS->SetFillColor(0);
	legendRpPbPHOS->SetLineColor(0);
	legendRpPbPHOS->SetTextFont(42);
	legendRpPbPHOS->SetNColumns(1);
	legendRpPbPHOS->SetTextSize(0.04);
	legendRpPbPHOS->AddEntry(graphRpPbPHOSStatErr,"PHOS Tsubasa","pef");
	if (System.CompareTo("PHOS")==0)	legendRpPbPHOS->AddEntry(graphRpPbStatErr,Form("%s, |#it{y}| < 0.13",SystemLabel.Data()),"p");
	else if (System.CompareTo("EMCal")==0)	legendRpPbPHOS->AddEntry(graphRpPbStatErr,Form("%s, |#it{y}| < 0.5",SystemLabel.Data()),"p");
	else	legendRpPbPHOS->AddEntry(graphRpPbStatErr,Form("%s, |#it{y}| < 0.8",SystemLabel.Data()),"p");
	legendRpPbPHOS->Draw();
	canvasRpPbPHOS->Update();
	
	

	canvasRpPbPHOS->Update(); 
	
	canvasRpPbPHOS->SaveAs(Form("%s/RpPb_CompPHOS_%s.%s",outputDir.Data(),System.Data(),suffix.Data()));
		
	//////////////////////////////Comparing to Charged Particles ///////////////////////////////////////////////
	
	
	
		
	
	DrawGammaSetMarkerTGraphAsym(graphAsymmChargedParticlesRpPbSystErr,20,1.5, kOrange+2, kOrange+2, 1, kTRUE,kOrange-9);

	
	
	DrawGammaSetMarkerTGraphAsym(graphAsymmChargedParticlesRpPbStatErr,20,1.5, kOrange+2, kOrange+2);
	graphAsymmChargedParticlesRpPbStatErr->SetFillColor(0);
	
	
	
	
	
	TCanvas* canvasRpPbChargedPar = new TCanvas("canvasRpPbChargedParCombined","Combined R_{pPb} with charged particles",200,10,1000,700);
	
	DrawGammaCanvasSettings( canvasRpPbChargedPar,  0.1, 0.01, 0.015, 0.11);
	
	
	
	
	TH2F * histo2DRpPbChargedPar = new TH2F("histo2DRpPbChargedPar","histo2DRpPbChargedPar",1000,0.,pTLimit,1000,0.0,5.0);
	SetStyleHistoTH2ForGraphs(histo2DRpPbChargedPar,pTLabel.Data(),RpPbLabel.Data(), 0.04,0.04, 0.04,0.04, 1.2,1.0, 512, 512); 
	
	
	//	histo2DRpPbChargedPar->GetXaxis()->SetRangeUser(0.0,14.5);
	histo2DRpPbChargedPar->GetYaxis()->SetRangeUser(0.19,1.81);
	histo2DRpPbChargedPar->GetYaxis()->CenterTitle();
	
	histo2DRpPbChargedPar->Draw();
	
	
	
	
	graphAsymmChargedParticlesRpPbSystErr->Draw("2p,same,e");
	graphAsymmChargedParticlesRpPbStatErr->Draw("p,same,e");
	
	
	graphRpPbSystErr->Draw("2p,same,e");
	graphRpPbStatErr->Draw("p,same,e");
	
	
	DrawGammaLines(0., 14.5,1., 1.,1.,kBlack,2);
	
	TLatex *textRpPb2 = new TLatex(0.15,0.87,energyLabel.Data());
	SetStyleTLatex( textRpPb2,0.04,4); 
	textRpPb2->Draw();
	
	TLatex *LabelpPbChargedPar = new TLatex(0.15,0.92,thesisPlotLabel.Data());
	SetStyleTLatex( LabelpPbChargedPar, 0.04,4);
	LabelpPbChargedPar->Draw();
	
	
	 TLegend* legendRpPb2 = new TLegend(0.2,0.15,0.65,0.3);
	legendRpPb2->SetFillColor(0);
	legendRpPb2->SetLineColor(0);
	legendRpPb2->SetTextSize(0.03);
	legendRpPb2->SetTextFont(42);
	legendRpPb2->SetTextSize(0.04);
	if (System.CompareTo("PHOS")==0)	legendRpPb2->AddEntry(graphRpPbStatErr,Form("%s, |#it{y}| < 0.13",SystemLabel.Data()),"p");
	else if (System.CompareTo("EMCal")==0)	legendRpPb2->AddEntry(graphRpPbStatErr,Form("%s, |#it{y}| < 0.5",SystemLabel.Data()),"p");
	else if (System.CompareTo("Comb")==0)	legendRpPb2->AddEntry(graphRpPbStatErr,"#pi^{0} combined","p");
	else	legendRpPb2->AddEntry(graphRpPbStatErr,Form("%s, |#it{y}| < 0.8",SystemLabel.Data()),"p");
	legendRpPb2->AddEntry(graphAsymmChargedParticlesRpPbStatErr,"ALICE, NSD, charged particles, |#eta_{cms}| < 0.3","p");
	
	legendRpPb2->Draw();
	
	
	
	
	
	
	
	
	canvasRpPbChargedPar->SaveAs(Form("%s/RpPb_Comparison_%s_ChargedParticles.%s",outputDir.Data(),System.Data(),suffix.Data()));

	cout<<"Llego aqui antes theory"<<endl;
	
	//return;
	
	//****************************** extracting EPS09s predictions**************************************
	
	ifstream 		inDSS;
	ifstream                inAKK;
	ifstream                inKKP;
	ifstream 		inCGC;	
	
	
	Int_t nlinesEPSsPi0fDSS = 0;
	
	inDSS.open(fileNameEPS09sPi0DSS,ios_base::in);
	
	Double_t xEPSsPi0fDSS[100],yEPSsPi0fDSS[100];
	
	Double_t xUpErrorEPSsPi0DSS[100],xDownErrorEPSsPi0DSS[100];
	Double_t yUpErrorEPSsPi0DSS[100],yDownErrorEPSsPi0DSS[100];
	
	
	while(!inDSS.eof()){
			nlinesEPSsPi0fDSS++;
			inDSS >> xEPSsPi0fDSS[nlinesEPSsPi0fDSS]  >> yEPSsPi0fDSS[nlinesEPSsPi0fDSS] >> yUpErrorEPSsPi0DSS[nlinesEPSsPi0fDSS]>>yDownErrorEPSsPi0DSS[nlinesEPSsPi0fDSS];
			yUpErrorEPSsPi0DSS[nlinesEPSsPi0fDSS] = ( yUpErrorEPSsPi0DSS[nlinesEPSsPi0fDSS] - yEPSsPi0fDSS[nlinesEPSsPi0fDSS] ) /yEPSsPi0fDSS[nlinesEPSsPi0fDSS];
			yDownErrorEPSsPi0DSS[nlinesEPSsPi0fDSS] = -1*( yDownErrorEPSsPi0DSS[nlinesEPSsPi0fDSS] - yEPSsPi0fDSS[nlinesEPSsPi0fDSS] ) /yEPSsPi0fDSS[nlinesEPSsPi0fDSS];
			xUpErrorEPSsPi0DSS[nlinesEPSsPi0fDSS]   = 0;
			xDownErrorEPSsPi0DSS[nlinesEPSsPi0fDSS] = 0;
			
			    
			cout << nlinesEPSsPi0fDSS << "         "  << xEPSsPi0fDSS[nlinesEPSsPi0fDSS] << "         "  <<yEPSsPi0fDSS[nlinesEPSsPi0fDSS]<<"         "<<yUpErrorEPSsPi0DSS[nlinesEPSsPi0fDSS]<<"          "<<yDownErrorEPSsPi0DSS[nlinesEPSsPi0fDSS]<< endl;
	
	}
	inDSS.close();
	
	cout<<"Llego aqui primer theory file"<<endl;
	
	graphPi0DSS5000 = new TGraph(nlinesEPSsPi0fDSS,xEPSsPi0fDSS,yEPSsPi0fDSS);
	graphAsymmErrorsPi0DSS5000 = new TGraphAsymmErrors(nlinesEPSsPi0fDSS,xEPSsPi0fDSS,yEPSsPi0fDSS,xDownErrorEPSsPi0DSS,xUpErrorEPSsPi0DSS,yDownErrorEPSsPi0DSS,yUpErrorEPSsPi0DSS);
	
	
	
	
	Int_t nlinesPi0AKK = 0;
		
	inAKK.open(fileNameEPS09sPi0AKK,ios_base::in);
	
	Double_t xESPsPi0AKK[100], yESPsPi0AKK[100];
	
		
	while(!inAKK.eof()){
			nlinesPi0AKK++;
			inAKK >> xESPsPi0AKK[nlinesPi0AKK]  >> yESPsPi0AKK[nlinesPi0AKK];
			cout << nlinesPi0AKK << "         "  << xESPsPi0AKK[nlinesPi0AKK] << "         "  <<xESPsPi0AKK[nlinesPi0AKK]<<endl;
	
	}
	inAKK.close();
	
	graphPi0ESP09sPi0AKK = new TGraph(nlinesPi0AKK,xESPsPi0AKK,yESPsPi0AKK);	
	
	
	Int_t nlinesPi0KKP = 0;
		
	inKKP.open(fileNameEPS09sPi0KKP,ios_base::in);
	
	Double_t xESPsPi0KKP[100], yESPsPi0KKP[100];
	
		
	while(!inKKP.eof()){
			nlinesPi0KKP++;
			inKKP >> xESPsPi0KKP[nlinesPi0KKP]  >> yESPsPi0KKP[nlinesPi0KKP];
			cout << nlinesPi0KKP << "         "  << xESPsPi0KKP[nlinesPi0KKP] << "         "  <<xESPsPi0KKP[nlinesPi0KKP]<<endl;
	
	}
	inKKP.close();
	
	
	cout<<"Final 2 "<<endl;
	//return;
	
	graphPi0ESP09sPi0KKP = new TGraph(nlinesPi0KKP,xESPsPi0KKP,yESPsPi0KKP);	
	
		
	
	graphPi0ESP09sPi0KKP->RemovePoint(0);
	graphPi0ESP09sPi0AKK->RemovePoint(0);
	graphPi0DSS5000->RemovePoint(0);
	graphAsymmErrorsPi0DSS5000->RemovePoint(0);
	
	graphPi0ESP09sPi0KKP->SetLineColor(kRed);
	graphPi0ESP09sPi0KKP->SetLineWidth(3);
	graphPi0ESP09sPi0KKP->SetLineStyle(10);
	
	
	//graphPi0ESP09sPi0KKP->RemovePoint(0);
	graphPi0ESP09sPi0AKK->SetLineColor(kBlue);
	graphPi0ESP09sPi0AKK->SetLineWidth(3);
	graphPi0ESP09sPi0AKK->SetLineStyle(7);
	
	
	
	graphPi0DSS5000->SetLineColor(8);
	graphPi0DSS5000->SetLineWidth(3);
	
	

	graphAsymmErrorsPi0DSS5000->SetFillColor(kYellow-7);
	graphAsymmErrorsPi0DSS5000->SetFillStyle(1001);
		
	
        


	Int_t nlinesPi0CGC = 0;	
	inCGC.open(fileNameEPS09sPi0CGC,ios_base::in);
	
	Double_t xESPsPi0CGC[100], yESPsPi0CGC[100];
	
		
	while(!inCGC.eof()){
			nlinesPi0CGC++;
			inCGC >> xESPsPi0CGC[nlinesPi0CGC]  >> yESPsPi0CGC[nlinesPi0CGC];
			cout << nlinesPi0CGC << "         "  << xESPsPi0CGC[nlinesPi0CGC] << "         "  <<yESPsPi0CGC[nlinesPi0CGC]<<endl;
	
	}
	inCGC.close();
	

	
	
	
	
	cout<<"Final 3"<<endl;
	
	//return;
	
	
	
	
	
	
	
	graphPi0CGC = new TGraph(nlinesPi0CGC,xESPsPi0CGC,yESPsPi0CGC);	
	
	TCanvas* canvasRpPbESP09sCGC = new TCanvas("canvasRpPbESP09sCGC"," Combined predictions for R^{#pi0}_{pPb}",200,10,1000,700);
	
	
	DrawGammaCanvasSettings( canvasRpPbESP09sCGC,  0.1, 0.01, 0.015, 0.13);

	
	
	TH2F * histo2DRpPbESP09s = new TH2F("histo2DRpPbESP09s","histo2DRpPbESP09s",1000,0.,pTLimit,1000,0.0,5.0);
	SetStyleHistoTH2ForGraphs(histo2DRpPbESP09s,pTLabel.Data(),RpPbLabel.Data(), 0.04,0.04, 0.04,0.04, 1.2,1.0, 512, 512); 
	
	
	//	histo2DRpPbESP09s->GetXaxis()->SetRangeUser(0.0,14.5);
	histo2DRpPbESP09s->GetYaxis()->SetRangeUser(0.19,1.81);
	histo2DRpPbESP09s->GetYaxis()->CenterTitle();
	
	
	histo2DRpPbESP09s->DrawCopy(); 
	
	
	//TH2F * histo2DRpPbESP09sCGC = new TH2F("histo2DRpPbESP09sCGC","histo2DRpPbESP09sCGC",1500,0.5,15.,1000,0.0,1.5);
	//SetStyleHistoTH2ForGraphs(histo2DRpPbESP09sCGC, "p_{T} (GeV/c)","R^{#pi^{0}}_{pPb}", 0.05,0.064, 0.05,0.06, 0.8,0.7, 512, 505); 
	//histo2DRpPbESP09sCGC->DrawCopy(); 
	
	
	DrawGammaLines(0., 14.5,1., 1.,2.,kBlack);
		
	
        Float_t offsetYMaxLegend = -0.15;
	Float_t xMinTheory = 0.30;
	
	TLatex *LabelRpPbESP09sCGC = new TLatex(0.15,0.92,thesisPlotLabel.Data());
	SetStyleTLatex( LabelRpPbESP09sCGC, 0.04,4);
	
	
        TLatex *textRpPbCompared = new TLatex(0.15,0.87,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
	SetStyleTLatex( textRpPbCompared,0.04,4); 
	
	TLatex *textRpPbEPS09s = new TLatex(xMinTheory,0.55 + offsetYMaxLegend,"p-Pb, #sqrt{#it{s}_{NN}} = 5.0 TeV");
	SetStyleTLatex( textRpPbEPS09s,0.04,4); 
	
	TLatex *textRpPbEPS09s2 = new TLatex(xMinTheory,0.5 + offsetYMaxLegend,"pQCD calculations");
	SetStyleTLatex( textRpPbEPS09s2,0.04,4); 
	
	
	TLatex *textRpPbESP09sCGC2= new TLatex(xMinTheory,0.23,"CGC calculations");
	SetStyleTLatex( textRpPbESP09sCGC2,0.04,4); 	
	
	TLegend* legendRpPbESP09sCGC = new TLegend(0.15,0.79,0.45,0.84);
	legendRpPbESP09sCGC->SetFillColor(0);
	legendRpPbESP09sCGC->SetLineColor(0);
	legendRpPbESP09sCGC->SetTextSize(0.04);
	legendRpPbESP09sCGC->SetTextFont(42);

	if (System.CompareTo("PHOS")==0)	legendRpPbESP09sCGC->AddEntry(graphRpPbStatErr,Form("%s, |#it{y}| < 0.13",SystemLabel.Data()),"p");
	else if (System.CompareTo("EMCal")==0)	legendRpPbESP09sCGC->AddEntry(graphRpPbStatErr,Form("%s, |#it{y}| < 0.5",SystemLabel.Data()),"p");
	else if (System.CompareTo("Comb")==0)	legendRpPbESP09sCGC->AddEntry(graphRpPbStatErr,"#pi^{0} combined","p");
	else	legendRpPbESP09sCGC->AddEntry(graphRpPbStatErr,Form("%s, |#it{y}| < 0.8",SystemLabel.Data()),"p");
	
	TLegend* legendRpPbESP09s = new TLegend(xMinTheory,0.3+offsetYMaxLegend,0.55,0.48+offsetYMaxLegend);
	legendRpPbESP09s->SetFillColor(0);
	legendRpPbESP09s->SetLineColor(0);
	legendRpPbESP09s->SetTextSize(0.04);
	legendRpPbESP09s->SetTextFont(42);
	legendRpPbESP09s->AddEntry(graphPi0ESP09sPi0KKP,"EPS09s KKP NLO");
	legendRpPbESP09s->AddEntry(graphPi0ESP09sPi0AKK,"EPS09s AKK NLO");
	legendRpPbESP09s->AddEntry(graphPi0DSS5000,"EPS09s fDSS NLO");
	legendRpPbESP09s->AddEntry(graphAsymmErrorsPi0DSS5000,"EPS09s fDSS errors","pef");
	
	
	TLegend* legendRpPbCGC = new TLegend(xMinTheory+.3,0.33+(offsetYMaxLegend/4),0.85,0.48+offsetYMaxLegend);
	legendRpPbCGC->SetFillColor(0);
	legendRpPbCGC->SetLineColor(0);
	legendRpPbCGC->SetTextSize(0.04);
	legendRpPbCGC->SetTextFont(42);
	legendRpPbCGC->AddEntry(graphPi0CGC,"Color Glas Condensate","p");

	graphPi0ESP09sPi0KKP->SetFillColor(0);
	graphPi0ESP09sPi0AKK->SetFillColor(0);
	graphPi0DSS5000->SetFillColor(0);

	graphPi0CGC->SetFillColor(0);
	graphPi0CGC->SetMarkerColor(kMagenta-2);
	graphPi0CGC->SetMarkerStyle(21);
	graphPi0CGC->SetMarkerSize(1.5);
	graphPi0CGC->SetLineColor(kGreen+3);
	graphPi0CGC->SetLineWidth(3);
	graphPi0CGC->SetLineStyle(1);

	graphPi0ESP09sPi0KKP->SetFillColor(0);
	graphPi0ESP09sPi0AKK->SetFillColor(0);
	graphPi0DSS5000->SetFillColor(0);

	
	
	histo2DRpPbESP09s->DrawCopy();
	legendRpPbESP09sCGC->Draw("same");		
	textRpPbCompared->Draw("same");
	textRpPbEPS09s->Draw("same");
 	textRpPbEPS09s2->Draw("same");
	LabelRpPbESP09sCGC->Draw("same");
	
	legendRpPbESP09s->Draw("same");
	legendRpPbCGC->Draw("same");
	

	graphAsymmErrorsPi0DSS5000->Draw("same,E3");
	graphPi0CGC->Draw("same,p");
	graphPi0DSS5000->Draw("same,p,l");
	graphPi0ESP09sPi0AKK->Draw("same,p,l");
	graphPi0ESP09sPi0KKP->Draw("same,p,l");
	graphRpPbSystErr->Draw("2p,same,e");
	graphRpPbStatErr->Draw("p,same,e");
	
	DrawGammaLines(0., 10.5,1., 1.,1.,kBlack,2);
	
	

	canvasRpPbESP09sCGC->Update(); 
	
	canvasRpPbESP09sCGC->SaveAs(Form("%s/RpPb_Comparison_%s_Models.%s",outputDir.Data(),System.Data(),suffix.Data()));
	

	
	
	//////////////////////////Comparison to Charged Pions/////////////////////////////////

	
	DrawGammaSetMarkerTGraphAsym(graphChargedPionRpPbSystErr,20,1.5, kRed, kRed, 1, kTRUE,kRed-10);

	
	DrawGammaSetMarkerTGraphAsym(graphChargedPionRpPbStatErr,20,1.5, kRed, kRed);
	graphAsymmChargedParticlesRpPbStatErr->SetFillColor(0);
	
	TCanvas* canvasRpPbChargedPions = new TCanvas("canvasRpPbChargedPions","Dalitz R_{pPb} with charged particles",200,10,1000,700);
	
	DrawGammaCanvasSettings( canvasRpPbChargedPions,  0.1, 0.01, 0.015, 0.11);
	
	
	TH2F * histo2DRpPbChargedPions = new TH2F("histo2DRpPbChargedPions","histo2DRpPbChargedPions",1000,0.,pTLimit,1000,0.0,5.0);
	SetStyleHistoTH2ForGraphs(histo2DRpPbChargedPions,pTLabel.Data(),RpPbLabel.Data(), 0.04,0.04, 0.04,0.04, 1.2,1.0, 512, 512); 
	
	
	//	histo2DRpPbChargedPions->GetXaxis()->SetRangeUser(0.0,14.5);
	histo2DRpPbChargedPions->GetYaxis()->SetRangeUser(0.19,1.81);
	histo2DRpPbChargedPions->GetYaxis()->CenterTitle();
	
	histo2DRpPbChargedPions->Draw(); 
	
	
	graphChargedPionRpPbSystErr->Draw("2p,same,e");
	graphChargedPionRpPbStatErr->Draw("p,same,e");
	
	
	graphRpPbSystErr->Draw("2p,same,e");
	graphRpPbStatErr->Draw("p,same,e");
	
	
	
	
	DrawGammaLines(0., 14.5,1., 1.,1.,kBlack,2);
	
	
	
	
	
	TLatex *textRpPbChargedPions = new TLatex(0.15,0.87,energyLabel.Data());
	SetStyleTLatex( textRpPbChargedPions,0.04,4); 
	textRpPbChargedPions->Draw();
	
	
	

	  
	
	TLatex *LabelpPbChargedPions = new TLatex(0.15,0.92,thesisPlotLabel.Data());
	SetStyleTLatex( LabelpPbChargedPions, 0.04,4);
	LabelpPbChargedPions->Draw();
	
	
	

	
	TLegend* legendRpPbChargedPions = new TLegend(0.2,0.150,0.65,0.3);
	legendRpPbChargedPions->SetFillColor(0);
	legendRpPbChargedPions->SetLineColor(0);
	legendRpPbChargedPions->SetTextSize(0.03);
	legendRpPbChargedPions->SetTextFont(42);
	legendRpPbChargedPions->SetTextSize(0.04);
	if (System.CompareTo("PHOS")==0)	legendRpPbChargedPions->AddEntry(graphRpPbStatErr,Form("%s, |#it{y}| < 0.13",SystemLabel.Data()),"p");
	else if (System.CompareTo("EMCal")==0)	legendRpPbChargedPions->AddEntry(graphRpPbStatErr,Form("%s, |#it{y}| < 0.5",SystemLabel.Data()),"p");
	else if (System.CompareTo("Comb")==0)	legendRpPbChargedPions->AddEntry(graphRpPbStatErr,"#pi^{0} combined","p");
	else	legendRpPbChargedPions->AddEntry(graphRpPbStatErr,Form("%s, |#it{y}| < 0.8",SystemLabel.Data()),"p");
	legendRpPbChargedPions->AddEntry(graphChargedPionRpPbStatErr,Form("%s, -0.5 < #it{y}_{cms} < 0   for #it{p}_{T} < 2.0 GeV/#it{c}",ChargedPionsLabel.Data()),"p");
	legendRpPbChargedPions->AddEntry((TObject*)0,"            -0.3 < #it{y}_{cms} < 0.3 for #it{p}_{T} > 2.0 GeV/#it{c}","");
	legendRpPbChargedPions->Draw();
	
	
	
	canvasRpPbChargedPions->SaveAs(Form("%s/RpPb_Comparison_%s_ChargedPions.%s",outputDir.Data(),System.Data(),suffix.Data()));


	if (System.CompareTo("PCM")==0){	
	/////
	TH1D* histoChargedPions =GraphAsymErrorsToHist_withErrors(graphChargedPionRpPbStatErr, "HistoChargedPionsRpPb");
	TH1D* histoNeutralPions =GraphAsymErrorsToHist_withErrors(graphRpPbStatErr, "HistoNeutralPionsRpPb");
	histoNeutralPions->Sumw2();
	histoChargedPions->Sumw2();
	TGraphErrors* graphRatioRpPb_ChargedPions    = CalculateRatioBetweenSpectraWithDifferentBinning(histoNeutralPions,graphRpPbSystErr ,histoChargedPions,graphChargedPionRpPbSystErr,  kTRUE,  kTRUE, &graphRebinnedAStat,&graphRebinnedASyst,&graphRebinnedBStat,&graphRebinnedBSyst);

	TCanvas* canvasRatio_ChargedPions = new TCanvas("canvasRatio_ChargedPions","Ratio RpPb NeutralPions ChargedPions",200,10,1200,700);
	
	DrawGammaCanvasSettings( canvasRatio_ChargedPions,  0.1, 0.01, 0.015, 0.13);
	
	TH2F * histo2DRatio_ChargedPions = new TH2F("histo2DRatio_ChargedPions","histo2DRatio_ChargedPions",1000,0.,pTLimit,1000,0.49,1.51);
	SetStyleHistoTH2ForGraphs(histo2DRatio_ChargedPions, "#it{p}_{T} (GeV/#it{c})","Ratio #it{R}_{pPb}", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 
	histo2DRatio_ChargedPions->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphErr(graphRatioRpPb_ChargedPions,21,1.5, kBlack , kBlack);
	graphRatioRpPb_ChargedPions->Draw("E1psame");
	

	DrawGammaLines(0., 15.,1., 1.,2.,kBlack);
	
	TLatex *textRatio_ChargedPions = new TLatex(0.16,0.9,"p-Pb minimum bias #sqrt{#it{s}_{NN}} = 5.02 TeV");
	SetStyleTLatex( textRatio_ChargedPions,0.06,4); 
	textRatio_ChargedPions->Draw();
	
	
	TLegend* legendRatio_ChargedPions = new TLegend(0.18,0.15,0.9,0.21);
	legendRatio_ChargedPions->SetFillColor(0);
	legendRatio_ChargedPions->SetLineColor(0);
	legendRatio_ChargedPions->SetNColumns(2);
	legendRatio_ChargedPions->SetTextSize(0.045);
	legendRatio_ChargedPions->AddEntry(graphRatioRpPb_ChargedPions,Form("Ratio #it{R}_{pPb} %s / #pi^{#pm} ",System.Data()),"p");
	legendRatio_ChargedPions->Draw();
	
	
	canvasRatio_ChargedPions->SaveAs(Form("%s/Ratio_%s_ChargedPions.%s",outputDir.Data(),System.Data(),suffix.Data()));

	}
	
	DrawGammaSetMarkerTGraphAsym(graphChargedKaonRpPbSystErr,20,1.5, colorsArray[kKaon],colorsArray[kKaon],1.0,kTRUE);
	DrawGammaSetMarkerTGraphAsym(graphChargedKaonRpPbStatErr,20,1.5, colorsArray[kKaon],colorsArray[kKaon]);
	graphChargedKaonRpPbStatErr->SetFillColor(0);
	
	
	DrawGammaSetMarkerTGraphAsym(graphChargedProtonRpPbSystErr,20,1.5, colorsArray[kProton], colorsArray[kProton], 1,kTRUE);
	DrawGammaSetMarkerTGraphAsym(graphChargedProtonRpPbStatErr,20,1.5, colorsArray[kProton], colorsArray[kProton]);
	graphChargedProtonRpPbStatErr->SetFillColor(0);
	
	DrawGammaSetMarkerTGraphAsym(graphChargedPionRpPbSystErrClone,20,1.5, colorsArray[kPion], colorsArray[kPion], 1,kTRUE);
	DrawGammaSetMarkerTGraphAsym(graphChargedPionRpPbStatErrClone,20,1.5, colorsArray[kPion], colorsArray[kPion]);
	graphChargedPionRpPbStatErrClone->SetFillColor(0);
	
	
	
	
	
	TCanvas* canvasRpPbChargedPionsKaonsProtons = new TCanvas("canvasRpPbChargedPionsKaonsProtons","R_{pPb} ",200,10,1000,800);
	DrawGammaCanvasSettings( canvasRpPbChargedPionsKaonsProtons,  0.1, 0.01, 0.015, 0.11);

	TH2F * histo2DRpPbChargedPionsKaonsProtons = new TH2F("histo2DRpPbChargedPionsKaonsProtons","histo2DRpPbChargedPionsKaonsProtons",1000,0.,pTLimit,1000,0.0,5.0);
	SetStyleHistoTH2ForGraphs(histo2DRpPbChargedPionsKaonsProtons,pTLabel.Data(),RpPbLabel.Data(), 0.04,0.04, 0.04,0.04, 1.2,1.0, 512, 512); 
	
	
	//	histo2DRpPbChargedPionsKaonsProtons->GetXaxis()->SetRangeUser(0.0,10.5);
	histo2DRpPbChargedPionsKaonsProtons->GetYaxis()->SetRangeUser(0.19,1.81);
	
	histo2DRpPbChargedPionsKaonsProtons->Draw(); 
	
	
	
	
	graphChargedKaonRpPbSystErr->Draw("2p,same,e");
	graphChargedKaonRpPbStatErr->Draw("p,same,e");
	graphChargedProtonRpPbSystErr->Draw("2p,same,e");
	graphChargedProtonRpPbStatErr->Draw("p,same,e");
	graphChargedPionRpPbSystErrClone->Draw("2p,same,e");
	graphChargedPionRpPbStatErrClone->Draw("p,same,e");
	
	
	graphRpPbSystErr->Draw("2p,same,e");
	graphRpPbStatErr->Draw("p,same,e");
	
	
	DrawGammaLines(0., 10.5,1., 1.,1.,kBlack,2);
	
	
	
	
	
	TLatex *textRpPbChargedPionsKaonsProtons = new TLatex(0.15,0.87,energyLabel.Data());
	SetStyleTLatex( textRpPbChargedPionsKaonsProtons,0.04,4); 
	textRpPbChargedPionsKaonsProtons->Draw();
	
	
	
	TLatex *LabelpPbChargedPionsKaonsProtons = new TLatex(0.15,0.92,thesisPlotLabel.Data());
	SetStyleTLatex( LabelpPbChargedPionsKaonsProtons, 0.04,4);
	LabelpPbChargedPionsKaonsProtons->Draw();

	
	TLegend* legendRpPbChargedPionsKaonsProtons = new TLegend(0.18,0.15,0.5,0.45);
	legendRpPbChargedPionsKaonsProtons->SetFillColor(0);
	legendRpPbChargedPionsKaonsProtons->SetLineColor(0);
	legendRpPbChargedPionsKaonsProtons->SetTextSize(0.03);
	legendRpPbChargedPionsKaonsProtons->SetTextFont(42);
	legendRpPbChargedPionsKaonsProtons->SetTextSize(0.03);
	if (System.CompareTo("PHOS")==0)	legendRpPbChargedPionsKaonsProtons->AddEntry(graphRpPbStatErr,Form("%s, |#it{y}| < 0.13",SystemLabel.Data()),"p");
	else if (System.CompareTo("EMCal")==0)	legendRpPbChargedPionsKaonsProtons->AddEntry(graphRpPbStatErr,Form("%s, |#it{y}| < 0.5",SystemLabel.Data()),"p");
	else if (System.CompareTo("Comb")==0)	legendRpPbChargedPionsKaonsProtons->AddEntry(graphRpPbStatErr,"#pi^{0} combined","p");
	else	legendRpPbChargedPionsKaonsProtons->AddEntry(graphRpPbStatErr,Form("%s, |#it{y}| < 0.8",SystemLabel.Data()),"p");
	legendRpPbChargedPionsKaonsProtons->AddEntry(graphChargedPionRpPbStatErrClone,"#pi^{+}+#pi^{-} -0.5 < #it{y}_{CMS} < 0 for #it{p}_{T} < 2.0 GeV/#it{c}","p");
	legendRpPbChargedPionsKaonsProtons->AddEntry((TObject*)0,"            -0.3 < #it{y}_{CMS}< 0.3 for #it{p}_{T} > 2.0 GeV/#it{c}","");
	legendRpPbChargedPionsKaonsProtons->AddEntry(graphChargedKaonRpPbSystErr,"K^{+} + K^{-} -0.5 < #it{y}_{CMS} < 0 for #it{p}_{T} < 2.0 GeV/#it{c}","p");
	legendRpPbChargedPionsKaonsProtons->AddEntry((TObject*)0,"            -0.3 < #it{y}_{CMS}< 0.3 for #it{p}_{T} > 2.0 GeV/#it{c}","");
	legendRpPbChargedPionsKaonsProtons->AddEntry(graphChargedProtonRpPbSystErr,"p + #bar{p} -0.5 <#it{y}_{CMS}<0 for #it{p}_{T} < 2.0 GeV/#it{c}","p");
	legendRpPbChargedPionsKaonsProtons->AddEntry((TObject*)0,"            -0.3<#it{y}_{CMS}<0.3 for #it{p}_{T} > 2.0 GeV/#it{c}","");
	
	legendRpPbChargedPionsKaonsProtons->Draw();
	
	
	canvasRpPbChargedPionsKaonsProtons->SaveAs(Form("%s/RpPb_Comparison_%s_ChargedPionsKaonsProtons.%s",outputDir.Data(),System.Data(),suffix.Data()));


	
	
	SavingFiles(outputDir, System);
	
	
	
	cout<<"Llego al final"<<endl;
	
	
	//////////////////////////Comparison to Combined RpPb /////////////////////////////////

	TFile* fCombined= new TFile("ExternalInputpPb/InputRpPb/ResultsRpPbpPb_2016_06_02.root");//without pile-up correction 	
	TGraphAsymmErrors* CombinePi0RpPbSystErr=(TGraphAsymmErrors*)fCombined->Get("CombinedPi0RpPbSystErr");
	TGraphAsymmErrors* CombinePi0RpPbStatErr=(TGraphAsymmErrors*)fCombined->Get("CombinedPi0RpPbStatErr");
	
	TCanvas* canvasRpPbCombinedPi0 = new TCanvas("canvasRpPbCombinedPi0","Dalitz R_{pPb} with charged particles",200,10,1000,700);
	
	DrawGammaCanvasSettings( canvasRpPbCombinedPi0,  0.1, 0.01, 0.015, 0.11);
	
	
	TH2F * histo2DRpPbCombinedPi0 = new TH2F("histo2DRpPbCombinedPi0","histo2DRpPbCombinedPi0",1000,0.,pTLimit,1000,0.0,5.0);
	SetStyleHistoTH2ForGraphs(histo2DRpPbCombinedPi0,pTLabel.Data(),RpPbLabel.Data(), 0.04,0.04, 0.04,0.04, 1.2,1.0, 512, 512); 
	
	
	//	histo2DRpPbCombinedPi0->GetXaxis()->SetRangeUser(0.0,14.5);
	histo2DRpPbCombinedPi0->GetYaxis()->SetRangeUser(0.19,1.81);
	histo2DRpPbCombinedPi0->GetYaxis()->CenterTitle();
	
	histo2DRpPbCombinedPi0->Draw(); 
	

	DrawGammaSetMarkerTGraphAsym(CombinePi0RpPbSystErr,24,1.5, kRed,kRed, 1,kTRUE);
	DrawGammaSetMarkerTGraphAsym(CombinePi0RpPbStatErr,24,1.5, kRed,kRed);


	graphRpPbSystErr->Draw("2p,same,e");
	graphRpPbStatErr->Draw("p,same,e");
	CombinePi0RpPbSystErr->Draw("E2,same");
	CombinePi0RpPbStatErr->Draw("p,same");	
	
	
	
	DrawGammaLines(0., 20.,1., 1.,1.,kBlack,2);
	
	TLatex *textRpPbCombinedPi0 = new TLatex(0.15,0.87,energyLabel.Data());
	SetStyleTLatex( textRpPbCombinedPi0,0.04,4); 
	textRpPbCombinedPi0->Draw();
	
	TLatex *LabelpPbCombinedPi0 = new TLatex(0.15,0.92,thesisPlotLabel.Data());
	SetStyleTLatex( LabelpPbCombinedPi0, 0.04,4);
	LabelpPbCombinedPi0->Draw();

	TLegend* legendRpPbCombinedPi0 = new TLegend(0.2,0.150,0.65,0.3);
	legendRpPbCombinedPi0->SetFillColor(0);
	legendRpPbCombinedPi0->SetLineColor(0);
	legendRpPbCombinedPi0->SetTextSize(0.03);
	legendRpPbCombinedPi0->SetTextFont(42);
	legendRpPbCombinedPi0->SetTextSize(0.04);
	if (System.CompareTo("PHOS")==0)	legendRpPbCombinedPi0->AddEntry(graphRpPbStatErr,Form("%s, |#it{y}| < 0.13",SystemLabel.Data()),"p");
	else if (System.CompareTo("EMCal")==0)	legendRpPbCombinedPi0->AddEntry(graphRpPbStatErr,Form("%s, |#it{y}| < 0.5",SystemLabel.Data()),"p");
	else if (System.CompareTo("Comb")==0)	legendRpPbCombinedPi0->AddEntry(graphRpPbStatErr,"#it{R}_{p-Pb} from combined #pi^{0} spectra","p");
	else	legendRpPbCombinedPi0->AddEntry(graphRpPbStatErr,Form("%s, |#it{y}| < 0.8",SystemLabel.Data()),"p");
	legendRpPbCombinedPi0->AddEntry(CombinePi0RpPbStatErr,"Combined individual #it{R}_{p-Pb}","p");

	legendRpPbCombinedPi0->Draw();
	
	
	
	canvasRpPbCombinedPi0->SaveAs(Form("%s/RpPb_Comparison_%s_CombinedPi0.%s",outputDir.Data(),System.Data(),suffix.Data()));

	


}

TH1F* GetPPReferenceFromPythia(TString fileName){
  
   TFile f1(fileName.Data());
   //TFile f1("pythia8-5.02TeV-part-1.root");
  // TFile f1("job_nd.12/pythia_pi0_non_diffractive.root");
   

   Double_t nEvt;
   Double_t xVal,yVal,yErr;
   
   TH1F* mult = (TH1F*)f1.Get("mult");
   TH1F* pi0PtVar = (TH1F*)f1.Get("pi0PtVar");
   TH1F* pi0Pt    = (TH1F*)f1.Get("pi0Pt");

   nEvt=mult->GetEntries()/1.e7;
   cout<<"Nevt="<< nEvt<<endl;
   
   
   pi0PtVar->Scale(1./(2*TMath::Pi()*nEvt));
   TH1F * OverPtpi0Pt=(TH1F*)pi0PtVar->Clone();
  //pi0Pt.Scale(1./(2*TMath::Pi()*nEvt));
  //   TH1F * OverPtpi0Pt=(TH1F*)pi0Pt->Clone();
     //Double_t sigma=6.789e+01;//48.03e-3;// cross section from Pythia8-tune 4C in mbar
     //Double_t sigma= 6.789e-04;
     //Double_t sigma=1.;


   for(Int_t ii=0;ii<pi0PtVar->GetNbinsX()+1;ii++){
   //for(Int_t ii=0;ii<pi0Pt->GetNbinsX()+1;ii++){

    xVal=pi0PtVar->GetBinCenter(ii+1);
    yVal=pi0PtVar->GetBinContent(ii+1);
    yErr=pi0Pt->GetBinError(ii+1);
    //OverPtpi0Pt->SetBinContent(ii+1,(yVal/(xVal*sigma))*recalcBarn);
    //OverPtpi0Pt->SetBinContent(ii+1,(yVal/xVal) * (sigma * recalcBarn) ); 
    OverPtpi0Pt->SetBinContent(ii+1,(yVal/xVal) * 1e09 ); 
    //    cout<< xVal <<" "<< OverPtpi0Pt.GetBinContent(ii+1)<<endl;

    // OverPtpi0Pt.SetBinError(ii+1,yErr/xVal);

    /*
    xVal=pi0Pt.GetBinCenter(ii+1);
    yVal=pi0Pt.GetBinContent(ii+1);
    yErr=pi0Pt.GetBinError(ii+1);
    OverPtpi0Pt.SetBinContent(ii+1,yVal/(xVal*sigma));
    // OverPtpi0Pt.SetBinError(ii+1,yErr/xVal);
    cout<< xVal <<" "<< OverPtpi0Pt.GetBinContent(ii+1)<<endl;
    */

  }
//OverPtpi0Pt->SetLineColor(2);
//OverPtpi0Pt->Draw("same");


 //if( Parameters2760GeV ) delete Parameters2760GeV;
 //if( Parameters7TeV )    delete Parameters7TeV;
	
  //canvasRpPbCombine->SaveAs(Form("%s/RpPb_Comparison__Dalitz.%s",outputDir.Data(),suffix.Data()));
  
	
 return OverPtpi0Pt;
}



void SavingFiles(TString outputDir, TString System){

      TString dateForOutput = ReturnDateStringForOutput();
  
      TFile fCombResults(Form("%s/Pi0RpPb_%s_%s.root",outputDir.Data(),System.Data(),dateForOutput.Data()),"RECREATE");
      graphAErrosInterPolation5023GeVBinStatErr->Write(Form("Pi0_pp_reference_%sBinning_StatErr",System.Data()));
      graphAErrosInterPolation5023GeVBinSystErr->Write(Form("Pi0_pp_reference_%sBinning_SystErr",System.Data())); 
      graphRpPbStatErr->Write(Form("Pi0_RpPb_%s_StatErr",System.Data()));
      graphRpPbSystErr->Write(Form("Pi0_RpPb_%s_SystErr",System.Data()));
      graphAsymmChargedParticlesRpPbStatErr->Write("RpPb_ChargedParticles_StatErr"); 
      graphAsymmChargedParticlesRpPbSystErr->Write("RpPb_ChargedParticles_SystErr");
      graphChargedPionRpPbStatErr->Write("RpPb_ChargedPions_StatErr"); 
      graphChargedPionRpPbSystErr->Write("RpPb_ChargedPions_SystErr");
      graphPi0ESP09sPi0KKP->Write("EPS09s_KKP_NLO");
      graphPi0ESP09sPi0AKK->Write("EPS09s_AKK_NLO");
      graphPi0DSS5000->Write("EPS09s_fDSS_NLO");
      graphAsymmErrorsPi0DSS5000->Write("EPS09s_fDSS_errors");
      graphPi0CGC->Write("ColorGlasCondensate");


      if( graphAlpha )     graphAlpha->Write(Form("Pi0_RpPb_%s_Alpha",System.Data()));
       
      if ( graphPtvsSqrts ){
      
	for ( Int_t i = 0;  i < graphAlpha->GetN(); i++)
	  if( graphPtvsSqrts[i] ) graphPtvsSqrts[i]->Write(Form("Pt_vs_Sqrts_%s_Bin_%d",System.Data(),i));
      }
      
  
}



