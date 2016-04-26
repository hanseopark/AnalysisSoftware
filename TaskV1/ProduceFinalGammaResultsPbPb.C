#include <Riostream.h>
#include <fstream>
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
#include "TGaxis.h"
#include "TMarker.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
//#include "AliHEPDataParser.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "CalculateGammaToPi0.h"

void ProduceFinalGammaResultsPbPb(TString inputFileName = "", TString cutSel = "", TString option = ""){

	// ******************************************************************
	// ******************** Set main style choices **********************
	// ******************************************************************
	StyleSettingsThesis();
	SetPlotStyle();

	// ******************************************************************
	// ***************** set ranges for plotting ************************
	// ******************************************************************
	Double_t doubleRatio[2];
	Double_t incRatio[2];
	Double_t doubleRatioX[2];
	TString cent ="";
	TString centralityCutNumber;
	Int_t maxBin		= 10;
	if(option.CompareTo("2.76TeV") == 0 || option.CompareTo("7TeV") == 0 ){
		doubleRatio[0] 		= 0.75;		doubleRatio[1] 		= 2.0;
		incRatio[0] 		= 0.0; 		incRatio[1] 		= 1.7;
		
		if (option.CompareTo("2.76TeV") == 0){
			cent = "2760GeV";
			doubleRatioX[0] 	= 0.; 		doubleRatioX[1] 	= 7.2;
			maxBin 				= 15;
		} else if (option.CompareTo("7TeV") == 0 ){
			cent = "7TeV";
			doubleRatioX[0] 	= 0.; 		doubleRatioX[1] 	= 12.2;
			maxBin 				= 22;
		}	
	} else if(option.CompareTo("PbPb_2.76TeV") == 0){
		centralityCutNumber = cutSel(GetEventSystemCutPosition(),3);
		doubleRatio[0] 		= 0.75; 	doubleRatio[1] 		= 2.0;
		incRatio[0] 		= 0.0; 		incRatio[1] 		= 1.7;
		doubleRatioX[0] 	= 0.8; 		doubleRatioX[1] 	= 14;
		cent 				= GetCentralityStringOutput(cutSel);
		maxBin 				= 18;
	} else if(option.CompareTo("pPb_5.023TeV") == 0){
		centralityCutNumber = cutSel(GetEventSystemCutPosition(),3);
		doubleRatio[0] 		= 0.75; 	doubleRatio[1] 		= 2.0;
		incRatio[0] 		= 0.0; 		incRatio[1] 		= 1.7;
		doubleRatioX[0] 	= 0.0; 		doubleRatioX[1] 	= 14;
		cent 				= "pPb";
		maxBin 				= 31;		
	}
	TString collisionSystem  = ReturnFullCollisionsSystem(option);
	cout << GetCentralityString(cutSel) << endl;
	if (GetCentralityString(cutSel).CompareTo("pp") != 0) collisionSystem = Form("%s %s", GetCentralityString(cutSel).Data(),collisionSystem.Data());
	
	// ******************************************************************
	// ********************** color definition **************************
	// ******************************************************************
	Color_t colorNLOcalc 				= kBlue-2;
	Color_t	colorCent					= GetColorDefaultColor(option, "", GetCentralityString(cutSel));
	Color_t	colorCentNotPi0Fitted		= GetColorDefaultColor(option, "", GetCentralityString(cutSel)) -7;
	Color_t colorErrA					= kGreen-8;
	Color_t colorErrB					= kRed-8;
	Color_t colorErrC					= kBlue-8;

	// ******************************************************************
	// ******************* Define output directory **********************
	// ******************************************************************
	TString outputDir = Form("FinalGammaResults/%s",option.Data());
	gSystem->Exec("mkdir -p "+outputDir);

	// ******************************************************************
	// ************************* read PCM data **************************
	// ******************************************************************
	TString nameDR 						= "DoubleRatioConversionTrueEffPurityABC" ;
	TString nameDRFit 					= "DoubleRatioConversionFitPurityABC" ;
	TString nameIRFit					= "histoIncRatioFitPurityA";
	if (option.CompareTo("pPb_5.023TeV") == 0){
		nameDR 							= "DoubleRatioConversionTrueEffPurity" ;
		nameDRFit 						= "DoubleRatioConversionFitPurity" ;
		nameIRFit						= "histoIncRatioFitPurity";
	}	
	TFile *fileInput 											= new TFile(inputFileName);	
	TH1D *histoDoubleRatio 										= (TH1D*) fileInput->Get(nameDR);
	TH1D *histoDoubleRatioPi0Fit 								= (TH1D*) fileInput->Get(nameDRFit);

	for (Int_t i = 1; i< histoDoubleRatioPi0Fit->GetNbinsX()+1; i++){
		cout << i << "\t" << histoDoubleRatioPi0Fit->GetBinCenter(i) << "\t" << histoDoubleRatioPi0Fit->GetBinContent(i) << "\t" << histoDoubleRatioPi0Fit->GetBinError(i) << endl;
	}	
	

	TH1D *histoIncRatio 										= (TH1D*) fileInput->Get("IncRatioPurity_trueEff");
	TH1D *histoIncRatioPi0Fit 									= (TH1D*) fileInput->Get(nameIRFit);
	TH1D *histoIncGammaSpec 									= (TH1D*) fileInput->Get("histoGammaSpecCorrPurity");
	TH1D *histoPi0 												= (TH1D*) fileInput->Get("CorrectedYieldTrueEff");
	
	TH1D *histoCocktailPi0Gamma									= (TH1D*) fileInput->Get("ptgPi0");
	TH1D *histoCocktailEtaGamma									= (TH1D*) fileInput->Get("ptgEta");
	TH1D *histoCocktailOmegaGamma								= (TH1D*) fileInput->Get("ptgOmega");
	TH1D *histoCocktailEtapGamma								= (TH1D*) fileInput->Get("ptgEtaprime");
	TH1D *histoCocktailPhiGamma									= (TH1D*) fileInput->Get("ptgPhi");
	TH1D *histoCocktailRhoGamma									= (TH1D*) fileInput->Get("ptgRho");
	TH1D *histoCocktailSigmaGamma								= (TH1D*) fileInput->Get("ptgSigma");
	
	TH1D *histoPi0Fit 											= (TH1D*) fileInput->Get("CorrectedYieldTrueEffPi0Fit");
	TH1D *histoCocktailAllGamma									= (TH1D*) fileInput->Get("ptg");
	Double_t paramsCocktail[6];
	GetFitParameter("qcd",GetCentralityString(cutSel),paramsCocktail);
	cout << paramsCocktail[0] << "\t" << paramsCocktail[1] << endl;
	TF1 *cocktailFitAllGammaForNLO 								= (TF1*)  FitObject("qcd","cocktailFit","Pi0",histoCocktailAllGamma,2.0,11,paramsCocktail,"QNRME+");
	cocktailFitAllGammaForNLO->SetRange(0,20);
	cout << WriteParameterToFile(cocktailFitAllGammaForNLO)<< endl;   
	
	TString fileNameSysErrDoubleRatio 							= Form("SystematicErrorsCalculated/SystematicErrorAveraged_DoubleRatio_PbPb_2.76TeV_%s_2015_07_29.dat",cent.Data());
	TString fileNameSysErrDoubleRatioPi0Fit 					= Form("SystematicErrorsCalculated/SystematicErrorAveraged_DoubleRatioPi0Fit_PbPb_2.76TeV_%s_2015_07_29.dat",cent.Data());
	TString fileNameSysErrIncGamma 								= Form("SystematicErrorsCalculated/SystematicErrorAveraged_GammaInc_PbPb_2.76TeV_%s_2015_07_29.dat",cent.Data());
	TString fileNameSysErrIncRatio 								= Form("SystematicErrorsCalculated/SystematicErrorAveraged_IncRatio_PbPb_2.76TeV_%s_2015_07_29.dat",cent.Data());
	TString fileNameSysErrIncRatioPi0Fit 						= Form("SystematicErrorsCalculated/SystematicErrorAveraged_IncRatioPi0Fit_PbPb_2.76TeV_%s_2015_07_29.dat",cent.Data());
	TString fileNameSysErrPi0 									= Form("SystematicErrorsCalculated/SystematicErrorAveraged_Pi0_PbPb_2.76TeV_%s_2015_07_29.dat",cent.Data());
	TString fileNameSysErrPi0Fit 								= Form("SystematicErrorsCalculated/SystematicErrorAveraged_Pi0Fit_PbPb_2.76TeV_%s_2015_07_29.dat",cent.Data());

	if (option.CompareTo("7TeV") == 0 || option.CompareTo("2.76TeV") == 0){
		fileNameSysErrDoubleRatio 							= Form("SystematicErrorsCalculated/SystematicErrorAveraged_DoubleRatio_%s_pp_19_Apr_2015.dat",option.Data());
		fileNameSysErrDoubleRatioPi0Fit 					= Form("SystematicErrorsCalculated/SystematicErrorAveraged_DoubleRatioPi0Fit_%s_pp_19_Apr_2015.dat",option.Data());
		fileNameSysErrIncGamma 								= Form("SystematicErrorsCalculated/SystematicErrorAveraged_GammaInc_%s_pp_19_Apr_2015.dat",option.Data());
		fileNameSysErrIncRatio 								= Form("SystematicErrorsCalculated/SystematicErrorAveraged_IncRatio_%s_pp_19_Apr_2015.dat",option.Data());
		fileNameSysErrIncRatioPi0Fit 						= Form("SystematicErrorsCalculated/SystematicErrorAveraged_IncRatioPi0Fit_%s_pp_19_Apr_2015.dat",option.Data());
		fileNameSysErrPi0 									= Form("SystematicErrorsCalculated/SystematicErrorAveraged_Pi0_%s_pp_19_Apr_2015.dat",option.Data());
		fileNameSysErrPi0Fit 								= Form("SystematicErrorsCalculated/SystematicErrorAveraged_Pi0Fit_%s_pp_19_Apr_2015.dat",option.Data());
	}	
	if (option.CompareTo("pPb_5.023TeV") == 0){
		fileNameSysErrDoubleRatio 							= "SystematicErrorsNew/SystematicErrorAveraged_DoubleRatio_pPb_5.023TeVMB_4_Apr_2014.dat";
		fileNameSysErrDoubleRatioPi0Fit 					= "SystematicErrorsNew/SystematicErrorAveraged_DoubleRatio_pPb_5.023TeVMB_4_Apr_2014.dat";
		fileNameSysErrIncGamma 								= "SystematicErrorsNew/SystematicErrorAveraged_GammaInc_pPb_5.023TeVMB_4_Apr_2014.dat"; 
		fileNameSysErrIncRatio 								= "SystematicErrorsNew/SystematicErrorAveraged_IncRatio_pPb_5.023TeVMB_4_Apr_2014.dat"; 
		fileNameSysErrIncRatioPi0Fit 						= "SystematicErrorsNew/SystematicErrorAveraged_IncRatio_pPb_5.023TeVMB_4_Apr_2014.dat"; 
		fileNameSysErrPi0 									= "SystematicErrorsNew/SystematicErrorAveraged_Pi0_pPb_5.023TeVMB_4_Apr_2014.dat"; 
		fileNameSysErrPi0Fit 								= "SystematicErrorsNew/SystematicErrorAveraged_Pi0_pPb_5.023TeVMB_4_Apr_2014.dat"; 
	}
	// ******************************************************************
	// *********** reading systematic errors for double ratio ***********
	// ******************************************************************
	ifstream 		fileSysErrDoubleRatio;
	fileSysErrDoubleRatio.open(fileNameSysErrDoubleRatio,ios_base::in);
	Double_t 		relSystErrorDoubleRatioUp[50];
	Double_t 		relSystErrorDoubleRatioDown[50];
	Double_t 		relSystErrorWOMaterialDoubleRatioUp[50];
	Double_t 		relSystErrorWOMaterialDoubleRatioDown[50];
	Double_t 		relSystErrorADoubleRatioUp[50];
	Double_t 		relSystErrorADoubleRatioDown[50];
	Double_t 		relSystErrorBDoubleRatioUp[50];
	Double_t 		relSystErrorBDoubleRatioDown[50];
	Double_t 		relSystErrorCDoubleRatioUp[50];
	Double_t 		relSystErrorCDoubleRatioDown[50];
	
	cout << fileNameSysErrDoubleRatio << endl;
	Int_t nPointsErrors=0;
	while(!fileSysErrDoubleRatio.eof() && nPointsErrors < 100){
		fileSysErrDoubleRatio 	>> relSystErrorDoubleRatioDown[nPointsErrors] 				>> relSystErrorDoubleRatioUp[nPointsErrors] 
								>> relSystErrorWOMaterialDoubleRatioDown[nPointsErrors] 	>> relSystErrorWOMaterialDoubleRatioUp[nPointsErrors] 
								>> relSystErrorADoubleRatioDown[nPointsErrors] 				>> relSystErrorADoubleRatioUp[nPointsErrors] 
								>> relSystErrorBDoubleRatioDown[nPointsErrors] 				>> relSystErrorBDoubleRatioUp[nPointsErrors] 
								>> relSystErrorCDoubleRatioDown[nPointsErrors] 				>> relSystErrorCDoubleRatioUp[nPointsErrors];
		cout << nPointsErrors << "\t"  << relSystErrorDoubleRatioDown[nPointsErrors] << "\t" << relSystErrorDoubleRatioUp[nPointsErrors] << "\t" 
			 << relSystErrorWOMaterialDoubleRatioDown[nPointsErrors] << "\t"  <<relSystErrorWOMaterialDoubleRatioUp[nPointsErrors] << endl;;		
		nPointsErrors++;
	}

	// ******************************************************************
	// ****** reading systematic errors for double ratio pi0 fit ********
	// ******************************************************************
	ifstream 		fileSysErrDoubleRatioPi0Fit;
	fileSysErrDoubleRatioPi0Fit.open(fileNameSysErrDoubleRatioPi0Fit, ios_base::in);
	Double_t 		relSystErrorDoubleRatioPi0FitUp[50];
	Double_t 		relSystErrorDoubleRatioPi0FitDown[50];
	Double_t 		relSystErrorWOMaterialDoubleRatioPi0FitUp[50];
	Double_t 		relSystErrorWOMaterialDoubleRatioPi0FitDown[50];
	Double_t 		relSystErrorADoubleRatioPi0FitUp[50];
	Double_t 		relSystErrorADoubleRatioPi0FitDown[50];
	Double_t 		relSystErrorBDoubleRatioPi0FitUp[50];
	Double_t 		relSystErrorBDoubleRatioPi0FitDown[50];
	Double_t 		relSystErrorCDoubleRatioPi0FitUp[50];
	Double_t 		relSystErrorCDoubleRatioPi0FitDown[50];
	
	cout << fileNameSysErrDoubleRatioPi0Fit << endl;
	nPointsErrors=0;
	while(!fileSysErrDoubleRatioPi0Fit.eof() && nPointsErrors < 100){
		fileSysErrDoubleRatioPi0Fit 	>> relSystErrorDoubleRatioPi0FitDown[nPointsErrors] 				>> relSystErrorDoubleRatioPi0FitUp[nPointsErrors] 
										>> relSystErrorWOMaterialDoubleRatioPi0FitDown[nPointsErrors] 		>> relSystErrorWOMaterialDoubleRatioPi0FitUp[nPointsErrors] 
										>> relSystErrorADoubleRatioPi0FitDown[nPointsErrors] 				>> relSystErrorADoubleRatioPi0FitUp[nPointsErrors] 
										>> relSystErrorBDoubleRatioPi0FitDown[nPointsErrors] 				>> relSystErrorBDoubleRatioPi0FitUp[nPointsErrors] 
										>> relSystErrorCDoubleRatioPi0FitDown[nPointsErrors] 				>> relSystErrorCDoubleRatioPi0FitUp[nPointsErrors];
		cout << nPointsErrors << "\t"  << relSystErrorDoubleRatioPi0FitDown[nPointsErrors] << "\t" << relSystErrorDoubleRatioPi0FitUp[nPointsErrors] << "\t" 
			 << relSystErrorWOMaterialDoubleRatioPi0FitDown[nPointsErrors] << "\t"  <<relSystErrorWOMaterialDoubleRatioPi0FitUp[nPointsErrors] << endl;;		
		nPointsErrors++;
	}

	// ******************************************************************
	// ** reading systematic errors for inclusive ratio with fitted pi0 *
	// ******************************************************************
	ifstream 		fileSysErrIncRatioPi0Fit;
	fileSysErrIncRatioPi0Fit.open(fileNameSysErrIncRatioPi0Fit,ios_base::in);
	Double_t 		relSystErrorIncRatioPi0FitUp[50];
	Double_t 		relSystErrorIncRatioPi0FitDown[50];
	Double_t 		relSystErrorWOMaterialIncRatioPi0FitUp[50];
	Double_t 		relSystErrorWOMaterialIncRatioPi0FitDown[50];
	Double_t 		relSystErrorAIncRatioPi0FitUp[50];
	Double_t 		relSystErrorAIncRatioPi0FitDown[50];
	Double_t 		relSystErrorBIncRatioPi0FitUp[50];
	Double_t 		relSystErrorBIncRatioPi0FitDown[50];
	Double_t 		relSystErrorCIncRatioPi0FitUp[50];
	Double_t 		relSystErrorCIncRatioPi0FitDown[50];
	
	cout << fileNameSysErrIncRatioPi0Fit << endl;
	nPointsErrors=0;
	while(!fileSysErrIncRatioPi0Fit.eof() && nPointsErrors < 100){
		fileSysErrIncRatioPi0Fit	>> relSystErrorIncRatioPi0FitDown[nPointsErrors] 				>> relSystErrorIncRatioPi0FitUp[nPointsErrors] 
									>> relSystErrorWOMaterialIncRatioPi0FitDown[nPointsErrors] 		>> relSystErrorWOMaterialIncRatioPi0FitUp[nPointsErrors] 
									>> relSystErrorAIncRatioPi0FitDown[nPointsErrors] 				>> relSystErrorAIncRatioPi0FitUp[nPointsErrors] 
									>> relSystErrorBIncRatioPi0FitDown[nPointsErrors] 				>> relSystErrorBIncRatioPi0FitUp[nPointsErrors] 
									>> relSystErrorCIncRatioPi0FitDown[nPointsErrors] 				>> relSystErrorCIncRatioPi0FitUp[nPointsErrors];
		cout << nPointsErrors << "\t"  << relSystErrorIncRatioPi0FitDown[nPointsErrors] << "\t"  << relSystErrorIncRatioPi0FitUp[nPointsErrors] << "\t" 
		     << relSystErrorWOMaterialIncRatioPi0FitDown[nPointsErrors] << "\t" << relSystErrorWOMaterialIncRatioPi0FitUp[nPointsErrors] << endl;;		
		nPointsErrors++;
	}
	
	// ******************************************************************
	// ********* reading systematic errors for inclusive ratio **********
	// ******************************************************************
	ifstream 		fileSysErrIncRatio;
	fileSysErrIncRatio.open(fileNameSysErrIncRatio,ios_base::in);
	Double_t 		relSystErrorIncRatioUp[50];
	Double_t 		relSystErrorIncRatioDown[50];
	Double_t 		relSystErrorWOMaterialIncRatioUp[50];
	Double_t 		relSystErrorWOMaterialIncRatioDown[50];
	Double_t 		relSystErrorAIncRatioUp[50];
	Double_t 		relSystErrorAIncRatioDown[50];
	Double_t 		relSystErrorBIncRatioUp[50];
	Double_t 		relSystErrorBIncRatioDown[50];
	Double_t 		relSystErrorCIncRatioUp[50];
	Double_t 		relSystErrorCIncRatioDown[50];
	
	cout << fileNameSysErrIncRatio << endl;
	nPointsErrors=0;
	while(!fileSysErrIncRatio.eof() && nPointsErrors < 100){
		fileSysErrIncRatio 		>> relSystErrorIncRatioDown[nPointsErrors] 					>> relSystErrorIncRatioUp[nPointsErrors] 
								>> relSystErrorWOMaterialIncRatioDown[nPointsErrors] 		>> relSystErrorWOMaterialIncRatioUp[nPointsErrors] 
								>> relSystErrorAIncRatioDown[nPointsErrors] 				>> relSystErrorAIncRatioUp[nPointsErrors] 
								>> relSystErrorBIncRatioDown[nPointsErrors] 				>> relSystErrorBIncRatioUp[nPointsErrors] 
								>> relSystErrorCIncRatioDown[nPointsErrors] 				>> relSystErrorCIncRatioUp[nPointsErrors];
		cout << nPointsErrors << "\t"  << relSystErrorIncRatioDown[nPointsErrors] << "\t"  << relSystErrorIncRatioUp[nPointsErrors] << "\t" 
		     << relSystErrorWOMaterialIncRatioDown[nPointsErrors] << "\t" << relSystErrorWOMaterialIncRatioUp[nPointsErrors] << endl;;		
		nPointsErrors++;
	}

	// ******************************************************************
	// ********* reading systematic errors for inclusive gamma **********
	// ******************************************************************
	ifstream 		fileSysErrIncGamma;
	fileSysErrIncGamma.open(fileNameSysErrIncGamma,ios_base::in);
	Double_t 		relSystErrorIncGammaUp[50];
	Double_t 		relSystErrorIncGammaDown[50];
	Double_t 		relSystErrorWOMaterialIncGammaUp[50];
	Double_t 		relSystErrorWOMaterialIncGammaDown[50];
	Double_t 		relSystErrorAIncGammaUp[50];
	Double_t 		relSystErrorAIncGammaDown[50];
	Double_t 		relSystErrorBIncGammaUp[50];
	Double_t 		relSystErrorBIncGammaDown[50];
	Double_t 		relSystErrorCIncGammaUp[50];
	Double_t 		relSystErrorCIncGammaDown[50];
	
	cout << "File-name: Gamma inc: "<< fileNameSysErrIncGamma << endl;
	nPointsErrors=0;
	while(!fileSysErrIncGamma.eof() && nPointsErrors < 100){
		fileSysErrIncGamma 		>> relSystErrorIncGammaDown[nPointsErrors] 					>> relSystErrorIncGammaUp[nPointsErrors] 
								>> relSystErrorWOMaterialIncGammaDown[nPointsErrors] 		>> relSystErrorWOMaterialIncGammaUp[nPointsErrors] 
								>> relSystErrorAIncGammaDown[nPointsErrors] 				>> relSystErrorAIncGammaUp[nPointsErrors] 
								>> relSystErrorBIncGammaDown[nPointsErrors] 				>> relSystErrorBIncGammaUp[nPointsErrors] 
								>> relSystErrorCIncGammaDown[nPointsErrors] 				>> relSystErrorCIncGammaUp[nPointsErrors];
		cout << nPointsErrors << "\t"  << relSystErrorIncGammaDown[nPointsErrors] << "\t"  <<relSystErrorIncGammaUp[nPointsErrors] << "\t" 
		     << relSystErrorWOMaterialIncGammaDown[nPointsErrors] << "\t" << relSystErrorWOMaterialIncGammaUp[nPointsErrors] << endl;;		
		nPointsErrors++;
	}
	nPointsErrors--;
	cout << "here" << endl;
	
	// ******************************************************************
	// *************** reading systematic errors for pi0 ****************
	// ******************************************************************
	ifstream 		fileSysErrPi0;
	fileSysErrPi0.open(fileNameSysErrPi0,ios_base::in);
	Double_t 		relSystErrorPi0Up[50];
	Double_t 		relSystErrorPi0Down[50];
	Double_t 		relSystErrorWOMaterialPi0Up[50];
	Double_t 		relSystErrorWOMaterialPi0Down[50];
	
	cout << "File-name Pi0: " <<fileNameSysErrPi0 << endl;
	nPointsErrors=0;
	while(!fileSysErrPi0.eof() && nPointsErrors < 100){
		fileSysErrPi0 			>> relSystErrorPi0Down[nPointsErrors] 				>> relSystErrorPi0Up[nPointsErrors] 
								>> relSystErrorWOMaterialPi0Down[nPointsErrors] 	>> relSystErrorWOMaterialPi0Up[nPointsErrors] ;
		cout << nPointsErrors << "\t"  << relSystErrorPi0Down[nPointsErrors] << "\t"  <<relSystErrorPi0Up[nPointsErrors] << "\t" 
		     << relSystErrorWOMaterialPi0Down[nPointsErrors] << "\t" << relSystErrorWOMaterialPi0Up[nPointsErrors] << endl;;		
		nPointsErrors++;
	}
	nPointsErrors--;
	
	// ******************************************************************
	// *********** reading systematic errors for pi0s fitted ************
	// ******************************************************************
	ifstream 		fileSysErrPi0Fit;
	fileSysErrPi0Fit.open(fileNameSysErrPi0Fit,ios_base::in);
	Double_t 		relSystErrorPi0FitUp[50];
	Double_t 		relSystErrorPi0FitDown[50];
	Double_t 		relSystErrorWOMaterialPi0FitUp[50];
	Double_t 		relSystErrorWOMaterialPi0FitDown[50];
	
	cout << "File-name Pi0 Fit: " << fileNameSysErrPi0Fit << endl;
	nPointsErrors=0;
	while(!fileSysErrPi0Fit.eof() && nPointsErrors < 100){
		fileSysErrPi0Fit 		>> relSystErrorPi0FitDown[nPointsErrors] 			>> relSystErrorPi0FitDown[nPointsErrors] 
								>> relSystErrorWOMaterialPi0FitDown[nPointsErrors] 	>> relSystErrorWOMaterialPi0FitUp[nPointsErrors] ;
		cout << nPointsErrors << "\t"  << relSystErrorPi0FitDown[nPointsErrors] << "\t"  <<relSystErrorPi0FitUp[nPointsErrors] << "\t" 
		     << relSystErrorWOMaterialPi0FitDown[nPointsErrors] << "\t" << relSystErrorWOMaterialPi0FitUp[nPointsErrors] << endl;;		
		nPointsErrors++;
	}
	nPointsErrors--;
	
	// ******************************************************************
	// ****************** reading theory graphs *************************
	// ******************************************************************
	TGraphAsymmErrors *graphInvYieldPbPbTheoryCTEQ61EPS09 		= NULL;
	TGraphAsymmErrors *graphInvYieldPPTheoryCT10BFG2_pdfErr 	= NULL;
	TGraphAsymmErrors *graphInvYieldPbPbTheoryEPS09 			= NULL;
	TGraphAsymmErrors *graphInvYieldPPTheoryCT10BFG2_scale 		= NULL;
	TGraphAsymmErrors* graphDRPbPbCT10BFG2_pdfErr 				= NULL;
	TGraphAsymmErrors* graphDRPbPbCTEQ61EPS09 					= NULL;
	graphInvYieldPbPbTheoryCTEQ61EPS09 							= (TGraphAsymmErrors*)fileInput->Get("PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield");
	graphInvYieldPPTheoryCT10BFG2_pdfErr 						= (TGraphAsymmErrors*)fileInput->Get("pp276CT10BFG2_sum_pdferr_InvYield");

	graphInvYieldPbPbTheoryEPS09 								= (TGraphAsymmErrors*)fileInput->Get("PbPb276EPS09BFG2_sum_scale_InvYield");
	graphInvYieldPPTheoryCT10BFG2_scale 						= (TGraphAsymmErrors*)fileInput->Get("pp276CT10BFG2_sum_scale_InvYield");

	Bool_t activateTheoryPbPb									= kFALSE;
	if (graphInvYieldPbPbTheoryCTEQ61EPS09 && graphInvYieldPPTheoryCT10BFG2_pdfErr && graphInvYieldPbPbTheoryEPS09 &&   graphInvYieldPPTheoryCT10BFG2_scale)
		activateTheoryPbPb										= kTRUE;
// 	graphInvYieldPbPbTheoryCTEQ61EPS09->RemovePoint(graphInvYieldPbPbTheoryCTEQ61EPS09->GetN()-1);
// 	graphInvYieldPPTheoryCT10BFG2_pdfErr->RemovePoint(graphInvYieldPPTheoryCT10BFG2_pdfErr->GetN()-1);
// 	
// 	graphInvYieldPbPbTheoryEPS09->RemovePoint(graphInvYieldPbPbTheoryEPS09->GetN()-1);
// 	graphInvYieldPPTheoryCT10BFG2_scale->RemovePoint(graphInvYieldPPTheoryCT10BFG2_scale->GetN()-1);

	if (activateTheoryPbPb){
		for(Int_t i = 0;i<graphInvYieldPbPbTheoryCTEQ61EPS09->GetN();i++){
			Double_t yerrlow1 		= graphInvYieldPbPbTheoryCTEQ61EPS09->GetErrorYlow(i);
			Double_t yerrlow2 		= 1*graphInvYieldPbPbTheoryEPS09->GetErrorYlow(i);
			Double_t yerrhigh1 		= graphInvYieldPbPbTheoryCTEQ61EPS09->GetErrorYhigh(i);
			Double_t yerrhigh2 		= 1*graphInvYieldPbPbTheoryEPS09->GetErrorYhigh(i);
			Double_t xerrlow 		= graphInvYieldPbPbTheoryEPS09->GetErrorXhigh(i);

			graphInvYieldPbPbTheoryCTEQ61EPS09->SetPointError( i,xerrlow,xerrlow,
														sqrt(yerrlow1*yerrlow1+ yerrlow2*yerrlow2),
														sqrt(yerrhigh1*yerrhigh1+ yerrhigh2*yerrhigh2));
		}
		
		cout << "======================================================================= DR EPS09 =================================================================" << endl;
		graphDRPbPbCTEQ61EPS09 									= (TGraphAsymmErrors*)graphInvYieldPbPbTheoryCTEQ61EPS09->Clone("graphDRPbPbCTEQ61EPS09");
		for (Int_t bin = 0; bin < graphDRPbPbCTEQ61EPS09->GetN(); bin++){
			Double_t cocktailIntegral 	= cocktailFitAllGammaForNLO->Integral(graphDRPbPbCTEQ61EPS09->GetX()[bin]-graphDRPbPbCTEQ61EPS09->GetErrorXlow(bin),
																			graphDRPbPbCTEQ61EPS09->GetX()[bin]+graphDRPbPbCTEQ61EPS09->GetErrorXhigh(bin)) 
										/(graphDRPbPbCTEQ61EPS09->GetErrorXlow(bin)+graphDRPbPbCTEQ61EPS09->GetErrorXhigh(bin));
			Double_t DRtheo	 			= 1 + (graphDRPbPbCTEQ61EPS09->GetY()[bin]/ cocktailIntegral);
			Double_t DRtheoErrDown		= graphDRPbPbCTEQ61EPS09->GetErrorYlow(bin)/ cocktailIntegral;
			Double_t DRtheoErrUp		= graphDRPbPbCTEQ61EPS09->GetErrorYhigh(bin)/ cocktailIntegral;
	// 		cout << graphDRPbPbCTEQ61EPS09->GetX()[bin]<< "\t" << cocktailIntegral << "\t" << DRtheo << "\t +" << DRtheoErrUp << "\t -" <<    DRtheoErrDown << endl;
			graphDRPbPbCTEQ61EPS09->SetPoint(bin, graphDRPbPbCTEQ61EPS09->GetX()[bin], DRtheo );
			graphDRPbPbCTEQ61EPS09->SetPointError(bin, graphDRPbPbCTEQ61EPS09->GetErrorXlow(bin), graphDRPbPbCTEQ61EPS09->GetErrorXhigh(bin), DRtheoErrDown, DRtheoErrUp );
		}	
		graphDRPbPbCTEQ61EPS09->Print();
			
		for(Int_t i = 0;i<graphInvYieldPPTheoryCT10BFG2_scale->GetN();i++){
			Double_t yerrlow1 		= graphInvYieldPPTheoryCT10BFG2_pdfErr->GetErrorYlow(i);
			Double_t yerrlow2 		= graphInvYieldPPTheoryCT10BFG2_scale->GetErrorYlow(i);
			Double_t yerrhigh1 		= graphInvYieldPPTheoryCT10BFG2_pdfErr->GetErrorYhigh(i);
			Double_t yerrhigh2 		= graphInvYieldPPTheoryCT10BFG2_scale->GetErrorYhigh(i);
			Double_t xerrlow 		= graphInvYieldPPTheoryCT10BFG2_scale->GetErrorXhigh(i);

			graphInvYieldPPTheoryCT10BFG2_pdfErr->SetPointError(i,xerrlow,xerrlow,
																sqrt(yerrlow1*yerrlow1+ yerrlow2*yerrlow2),
																sqrt(yerrhigh1*yerrhigh1+ yerrhigh2*yerrhigh2));
		}
		graphDRPbPbCT10BFG2_pdfErr 								= (TGraphAsymmErrors*)graphInvYieldPPTheoryCT10BFG2_pdfErr->Clone("graphDRPbPbCT10BFG2_pdfErr");
		for (Int_t bin = 0; bin < graphDRPbPbCT10BFG2_pdfErr->GetN(); bin++){
			Double_t cocktailIntegral 	= cocktailFitAllGammaForNLO->Integral(graphDRPbPbCT10BFG2_pdfErr->GetX()[bin]-graphDRPbPbCT10BFG2_pdfErr->GetErrorXlow(bin),
																			graphDRPbPbCT10BFG2_pdfErr->GetX()[bin]+graphDRPbPbCT10BFG2_pdfErr->GetErrorXhigh(bin)) 
										/(graphDRPbPbCT10BFG2_pdfErr->GetErrorXlow(bin)+graphDRPbPbCT10BFG2_pdfErr->GetErrorXhigh(bin));
			Double_t DRtheo	 			= 1 + (graphDRPbPbCT10BFG2_pdfErr->GetY()[bin]/ cocktailIntegral);
			Double_t DRtheoErrDown		= graphDRPbPbCT10BFG2_pdfErr->GetErrorYlow(bin)/ cocktailIntegral;
			Double_t DRtheoErrUp		= graphDRPbPbCT10BFG2_pdfErr->GetErrorYhigh(bin)/ cocktailIntegral;
	// 		cout << graphDRPbPbCT10BFG2_pdfErr->GetX()[bin]<< "\t" << cocktailIntegral << "\t" << DRtheo << "\t +" << DRtheoErrUp << "\t -" <<    DRtheoErrDown << endl;
			graphDRPbPbCT10BFG2_pdfErr->SetPoint(bin, graphDRPbPbCT10BFG2_pdfErr->GetX()[bin], DRtheo );
			graphDRPbPbCT10BFG2_pdfErr->SetPointError(bin, graphDRPbPbCT10BFG2_pdfErr->GetErrorXlow(bin), graphDRPbPbCT10BFG2_pdfErr->GetErrorXhigh(bin), DRtheoErrDown, DRtheoErrUp );
		}	
		graphDRPbPbCT10BFG2_pdfErr->Print();
	}
	
	// ******************************************************************
	// ***** read NL0 calculations and put them in proper format ********
	// ******************************************************************
	TGraphErrors *DirectPhotonDoubleNLOhalf 	= (TGraphErrors*) fileInput->Get("doubleRatioNLOhalf");
	TGraphErrors *DirectPhotonDoubleNLOone 		= (TGraphErrors*) fileInput->Get("doubleRatioNLOone");
	TGraphErrors *DirectPhotonDoubleNLOtwo  	= (TGraphErrors*) fileInput->Get("doubleRatioNLOtwo");
	DirectPhotonDoubleNLOhalf->RemovePoint(0);
	DirectPhotonDoubleNLOone->RemovePoint(0);
	DirectPhotonDoubleNLOtwo->RemovePoint(0);

	TGraphErrors *DirectPhotonNLOhalf 			= (TGraphErrors*) fileInput->Get("graphNLOCalcMuHalf");
	TGraphErrors *DirectPhotonNLOone 			= (TGraphErrors*) fileInput->Get("graphNLOCalcMuOne");
	TGraphErrors *DirectPhotonNLOtwo  			= (TGraphErrors*) fileInput->Get("graphNLOCalcMuTwo");
	DirectPhotonNLOhalf->RemovePoint(0);
	DirectPhotonNLOone->RemovePoint(0);
	DirectPhotonNLOtwo->RemovePoint(0);

	Double_t *errorup 							= new Double_t[DirectPhotonDoubleNLOhalf->GetN()];
	Double_t *errorlow 							= new Double_t[DirectPhotonDoubleNLOtwo->GetN()];
	Double_t *errorSpecup 						= new Double_t[DirectPhotonDoubleNLOhalf->GetN()];
	Double_t *errorSpeclow 						= new Double_t[DirectPhotonDoubleNLOtwo->GetN()];
	Double_t *yHalf 							= DirectPhotonDoubleNLOhalf->GetY();
	Double_t *yOne 								= DirectPhotonDoubleNLOone->GetY();
	Double_t *yTwo 								= DirectPhotonDoubleNLOtwo->GetY();
	Double_t *ySpecHalf 						= DirectPhotonNLOhalf->GetY();
	Double_t *ySpecOne 							= DirectPhotonNLOone->GetY();
	Double_t *ySpecTwo 							= DirectPhotonNLOtwo->GetY();

	for(Int_t i = 0;i<DirectPhotonDoubleNLOhalf->GetN(); i++){
		errorup[i]	 							= yHalf[i]-yOne[i];
		errorlow[i] 							= -yTwo[i]+yOne[i];
		errorSpecup[i] 							= ySpecHalf[i]-ySpecOne[i];
		errorSpeclow[i] 						= -ySpecTwo[i]+ySpecOne[i];
	}

	Int_t reduceBins = 33;
	TGraphAsymmErrors *NLODoubleRatio 			= new TGraphAsymmErrors( DirectPhotonDoubleNLOhalf->GetN()-reduceBins, DirectPhotonDoubleNLOone->GetX(), DirectPhotonDoubleNLOone->GetY(), 
																		 DirectPhotonDoubleNLOone->GetEX(), DirectPhotonDoubleNLOone->GetEX(), errorlow,errorup );
	TGraphAsymmErrors *NLO 						= new TGraphAsymmErrors( DirectPhotonDoubleNLOhalf->GetN()-reduceBins, DirectPhotonNLOone->GetX(), DirectPhotonNLOone->GetY(), DirectPhotonNLOone->GetEX(),
																		 DirectPhotonNLOone->GetEX(), errorSpeclow ,errorSpecup );
	
	SetStyleGammaNLOTGraphWithBand( NLODoubleRatio, 3.0, 1, colorNLOcalc, 1001, colorNLOcalc, 0);
	SetStyleGammaNLOTGraphWithBand( NLO, 3.0, 1, colorNLOcalc, 1001, colorNLOcalc, 0);

	// **************************************************************************
	// ******************** Draw Final Gamma Pictures ***************************
	// **************************************************************************
	TF1 *lineOne = new TF1("lineOne","1",0,16);
	lineOne->SetLineWidth(1.2);
	lineOne->SetLineColor(kGray+2);

	// **************************************************************************
	// ***************** Calculate systematic error graphs **********************
	// **************************************************************************
	TGraphAsymmErrors* graphDoubleRatioFitPi0SysErr 	= CalculateSysErrFromRelSysHisto( histoDoubleRatioPi0Fit, "DoubleRatioPi0FitSystError", relSystErrorDoubleRatioPi0FitDown, relSystErrorDoubleRatioPi0FitUp, 2,
																					   	  nPointsErrors);	
	TGraphAsymmErrors* graphDoubleRatioFitPi0SysErrA 	= CalculateSysErrFromRelSysHisto( histoDoubleRatioPi0Fit, "DoubleRatioPi0FitSystErrorA", relSystErrorADoubleRatioPi0FitDown, relSystErrorADoubleRatioPi0FitUp, 2,
																						  nPointsErrors);	
	TGraphAsymmErrors* graphDoubleRatioFitPi0SysErrB 	= CalculateSysErrFromRelSysHisto( histoDoubleRatioPi0Fit, "DoubleRatioPi0FitSystErrorB", relSystErrorBDoubleRatioPi0FitDown, relSystErrorBDoubleRatioPi0FitUp, 2,
																						  nPointsErrors);	
	TGraphAsymmErrors* graphDoubleRatioFitPi0SysErrC 	= CalculateSysErrFromRelSysHisto( histoDoubleRatioPi0Fit, "DoubleRatioPi0FitSystErrorC", relSystErrorCDoubleRatioPi0FitDown, relSystErrorCDoubleRatioPi0FitUp, 2,
																						  nPointsErrors);	
	TGraphAsymmErrors* graphDoubleRatioSysErr 			= CalculateSysErrFromRelSysHisto( histoDoubleRatio, "DoubleRatioSystError", relSystErrorDoubleRatioDown, relSystErrorDoubleRatioUp, 2,
																						  nPointsErrors);	
	TGraphAsymmErrors* graphDoubleRatioSysErrA 			= CalculateSysErrFromRelSysHisto( histoDoubleRatio, "DoubleRatioSystErrorA", relSystErrorADoubleRatioDown, relSystErrorADoubleRatioUp, 2,
																						  nPointsErrors);	
	TGraphAsymmErrors* graphDoubleRatioSysErrB 			= CalculateSysErrFromRelSysHisto( histoDoubleRatio, "DoubleRatioSystErrorB", relSystErrorBDoubleRatioDown, relSystErrorBDoubleRatioUp, 2,
																						  nPointsErrors);	
	TGraphAsymmErrors* graphDoubleRatioSysErrC 			= CalculateSysErrFromRelSysHisto( histoDoubleRatio, "DoubleRatioSystErrorC", relSystErrorCDoubleRatioDown, relSystErrorCDoubleRatioUp, 2,
																						  nPointsErrors);	

	//****************************************************************************
	//******************* draw double ratio with pi0 fitted **********************
	//****************************************************************************

	cout << "*************************************************************************************************************" << endl;
	cout << "*************************************************************************************************************" << endl;
	cout << "*************************************************************************************************************" << endl;
	
	TCanvas *canvasDoubleRatio = GetAndSetCanvas("canvasDoubleRatioFinal");
	if (!(option.CompareTo("7TeV") == 0 || option.CompareTo("2.76TeV") == 0 || option.CompareTo("pPb_5.023TeV")==0))canvasDoubleRatio->SetLogx();

	TH2D *dummyDR ;
	dummyDR = new TH2D("dummyDR", "dummyDR", 120, 0., 16, 1000., doubleRatio[0], doubleRatio[1]);
	SetStyleHistoTH2ForGraphs( dummyDR, "#it{p}_{T} (GeV/#it{c})", "(#it{N}_{#gamma_{inc}}/#it{N}_{#pi^{0}})/(#it{N}_{#gamma_{decay}}/#it{N}_{#pi^{0}})", 
							   0.045, 0.05, 0.045, 0.05, 0.85, 0.85);
	if (!(option.CompareTo("7TeV") == 0 || option.CompareTo("2.76TeV") == 0 || option.CompareTo("pPb_5.023TeV")==0)) dummyDR->GetXaxis()->SetLabelOffset(-0.015);
	dummyDR->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
	dummyDR->DrawCopy();

	for (Int_t i = 1; i< histoDoubleRatioPi0Fit->GetNbinsX()+1; i++){
		cout << i << "\t" << histoDoubleRatioPi0Fit->GetBinCenter(i) << "\t" << histoDoubleRatioPi0Fit->GetBinContent(i) << "\t" << histoDoubleRatioPi0Fit->GetBinError(i) << endl;
	}	
	graphDoubleRatioFitPi0SysErr->Print();
	
	DrawGammaSetMarker(histoDoubleRatioPi0Fit, 20, 2., colorCent, colorCent);                              	
	DrawGammaSetMarkerTGraphAsym(graphDoubleRatioFitPi0SysErr , 20, 2, colorCent, colorCent, 1, kTRUE);
	lineOne->Draw("same");
	NLODoubleRatio->Draw("p3lsame");
	
	graphDoubleRatioFitPi0SysErr->Draw("E2same");
	histoDoubleRatioPi0Fit->DrawCopy("same");

	TLegend* legendDoubleRatio = GetAndSetLegend(0.15,0.75,4);
	legendDoubleRatio->AddEntry(graphDoubleRatioFitPi0SysErr,Form("PCM, %s",collisionSystem.Data()),"pf");
	legendDoubleRatio->AddEntry(NLODoubleRatio,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
	legendDoubleRatio->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
	legendDoubleRatio->Draw();

	canvasDoubleRatio->Print(Form("%s/DoubleRatioPi0Fitted_%s.eps",outputDir.Data(),cent.Data()));

	//****************************************************************************
	// draw double ratio with pi0 fitted and compare to double ratio with pure pi0
	//****************************************************************************
	canvasDoubleRatio->cd();
	dummyDR->DrawCopy();

	lineOne->Draw("same");
	NLODoubleRatio->Draw("p3lsame");
	
	DrawGammaSetMarker(histoDoubleRatio, 20, 2., colorCentNotPi0Fitted, colorCentNotPi0Fitted);                              
	DrawGammaSetMarkerTGraphAsym(graphDoubleRatioSysErr , 20, 2, colorCentNotPi0Fitted, colorCentNotPi0Fitted, 1, kTRUE);
	graphDoubleRatioSysErr->Draw("E2same");
	histoDoubleRatio->DrawCopy("same");

	graphDoubleRatioFitPi0SysErr->Draw("E2same");
	histoDoubleRatioPi0Fit->DrawCopy("same");

	TLegend* legendDoubleRatio2 = GetAndSetLegend(0.15,0.75,4);
	legendDoubleRatio2->AddEntry(graphDoubleRatioFitPi0SysErr,Form("PCM, %s, #pi^{0} fitted",collisionSystem.Data()),"pf");
	legendDoubleRatio2->AddEntry(graphDoubleRatioSysErr,Form("PCM, %s",collisionSystem.Data()),"pf");
	legendDoubleRatio2->AddEntry(NLODoubleRatio,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
	legendDoubleRatio2->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
	legendDoubleRatio2->Draw();

	canvasDoubleRatio->Print(Form("%s/DoubleRatio_ComparedPi0FittedAndNot_%s.eps",outputDir.Data(),cent.Data()));
	
	
	//****************************************************************************
	//*************** Double Ratio with pi0 fitted with ABC errors ***************
	//****************************************************************************
	TCanvas *canvasDoubleRatioABC = GetAndSetCanvas("canvasDoubleRatioFinal");
	if (!(option.CompareTo("7TeV") == 0 || option.CompareTo("2.76TeV") == 0 || option.CompareTo("pPb_5.023TeV")==0)) canvasDoubleRatioABC->SetLogx();

	dummyDR->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
	dummyDR->DrawCopy();

	lineOne->Draw("same");
	NLODoubleRatio->Draw("p3lsame");
	histoDoubleRatioPi0Fit->DrawCopy("same");
	DrawGammaSetMarkerTGraphAsym(graphDoubleRatioFitPi0SysErrA , 20, 1, colorErrA, colorErrA, 1, kTRUE);
	DrawGammaSetMarkerTGraphAsym(graphDoubleRatioFitPi0SysErrB , 20, 1, colorErrB, colorErrB, 1, kTRUE);
	DrawGammaSetMarkerTGraphAsym(graphDoubleRatioFitPi0SysErrC , 20, 1, colorErrC, colorErrC, 1, kTRUE);

	graphDoubleRatioFitPi0SysErrA->Draw("E2same");
	graphDoubleRatioFitPi0SysErrB->Draw("E2same");
	graphDoubleRatioFitPi0SysErrC->Draw("E2same");
	graphDoubleRatioFitPi0SysErr->Draw("E2same");

	TLegend* legendDoubleRatioABC = GetAndSetLegend(0.15,0.65,6);
	legendDoubleRatioABC->AddEntry(graphDoubleRatioFitPi0SysErr,Form("PCM, %s",collisionSystem.Data()),"pf");
	legendDoubleRatioABC->AddEntry(graphDoubleRatioFitPi0SysErrA,"Error Type A","f");
	legendDoubleRatioABC->AddEntry(graphDoubleRatioFitPi0SysErrB,"Error Type B","f");
	legendDoubleRatioABC->AddEntry(graphDoubleRatioFitPi0SysErrC,"Error Type C","f");
	legendDoubleRatioABC->AddEntry(NLODoubleRatio,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
	legendDoubleRatioABC->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
	legendDoubleRatioABC->Draw();

	canvasDoubleRatioABC->Print(Form("%s/DoubleRatioPi0Fitted_ABCErrors_%s.eps",outputDir.Data(),cent.Data()));

	//************************************************************************
	//******************* Calculate error graph for inclusive ratio **********
	//************************************************************************
	cout << "here" << endl;
	TGraphAsymmErrors* graphIncRatioFitPi0SysErr 	= CalculateSysErrFromRelSysHisto( histoIncRatioPi0Fit, "IncRatioPi0FitSystError", relSystErrorIncRatioPi0FitDown, relSystErrorIncRatioPi0FitUp, 2,
																					  nPointsErrors);	
	cout << "here" << endl;
	TGraphAsymmErrors* graphIncRatioFitPi0SysErrA 	= CalculateSysErrFromRelSysHisto( histoIncRatioPi0Fit, "IncRatioPi0FitSystErrorA", relSystErrorAIncRatioPi0FitDown, relSystErrorAIncRatioPi0FitUp, 2,
																					  nPointsErrors);	
	TGraphAsymmErrors* graphIncRatioFitPi0SysErrB 	= CalculateSysErrFromRelSysHisto( histoIncRatioPi0Fit, "IncRatioPi0FitSystErrorB", relSystErrorBIncRatioPi0FitDown, relSystErrorBIncRatioPi0FitUp, 2,
																					  nPointsErrors);	
	TGraphAsymmErrors* graphIncRatioFitPi0SysErrC 	= CalculateSysErrFromRelSysHisto( histoIncRatioPi0Fit, "IncRatioPi0FitSystErrorC", relSystErrorCIncRatioPi0FitDown, relSystErrorCIncRatioPi0FitUp, 2,
																					  nPointsErrors);	
	TGraphAsymmErrors* graphIncRatioSysErr 			= CalculateSysErrFromRelSysHisto( histoIncRatio, "IncRatioSystError", relSystErrorIncRatioDown, relSystErrorIncRatioUp, 2,
																					  nPointsErrors);	
	TGraphAsymmErrors* graphIncRatioSysErrA 		= CalculateSysErrFromRelSysHisto( histoIncRatio, "IncRatioSystErrorA", relSystErrorAIncRatioDown, relSystErrorAIncRatioUp, 2,
																					  nPointsErrors);	
	TGraphAsymmErrors* graphIncRatioSysErrB 		= CalculateSysErrFromRelSysHisto( histoIncRatio, "IncRatioSystErrorB", relSystErrorBIncRatioDown, relSystErrorBIncRatioUp, 2,
																					  nPointsErrors);	
	TGraphAsymmErrors* graphIncRatioSysErrC 		= CalculateSysErrFromRelSysHisto( histoIncRatio, "IncRatioSystErrorC", relSystErrorCIncRatioDown, relSystErrorCIncRatioUp, 2,
																					  nPointsErrors);	
	// **************************************************************************
	// ******************** plotting inclusive ratio ****************************
	// **************************************************************************
	TCanvas *canvasIncRatio = GetAndSetCanvas("canvasIncRatioFinal");
	if (!(option.CompareTo("7TeV") == 0 || option.CompareTo("2.76TeV") == 0 || option.CompareTo("pPb_5.023TeV")==0)  ) canvasIncRatio->SetLogx(1);
	
	TH2D *dummyIncR ;
	dummyIncR = new TH2D("dummyIncR", "dummyIncR", 120, 0., 16, 1000., incRatio[0], incRatio[1]);
	SetStyleHistoTH2ForGraphs( dummyIncR, "#it{p}_{T} (GeV/#it{c})", "#gamma_{inc}/#pi^{0}",
							   0.045, 0.05, 0.045, 0.05, 0.85, 0.85);
	if (!(option.CompareTo("7TeV") == 0 || option.CompareTo("2.76TeV") == 0 || option.CompareTo("pPb_5.023TeV")==0)) dummyIncR->GetXaxis()->SetLabelOffset(-0.015);
	dummyIncR->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
	dummyIncR->DrawCopy();

	DrawGammaSetMarker(histoIncRatio, 20, 2., colorCentNotPi0Fitted, colorCentNotPi0Fitted);                
	DrawGammaSetMarkerTGraphAsym(graphIncRatioSysErr , 20, 2, colorCentNotPi0Fitted, colorCentNotPi0Fitted, 1, kTRUE);	
	graphIncRatioSysErr->Draw("E2same");
	histoIncRatio->DrawCopy("same");

	TLegend* legendIncRatio = GetAndSetLegend(0.12,0.2,1.5);
	legendIncRatio->AddEntry(graphIncRatioSysErr,Form("PCM, %s",collisionSystem.Data()),"pf");
	legendIncRatio->Draw();
	
	canvasIncRatio->Print(Form("%s/IncRatio_%s.eps",outputDir.Data(),cent.Data()));

	// **************************************************************************
	// ****************** plotting inclusive ratio with pi0 fitted **************
	// **************************************************************************
	
	TCanvas *canvasIncRatioPi0Fit = GetAndSetCanvas("canvasIncRatioPi0FitFinal");
	if (!(option.CompareTo("7TeV") == 0 || option.CompareTo("2.76TeV") == 0 || option.CompareTo("pPb_5.023TeV")==0)) canvasIncRatioPi0Fit->SetLogx(1);
	dummyIncR->DrawCopy();
	
	DrawGammaSetMarker(histoIncRatioPi0Fit, 20, 2., colorCent, colorCent);                
	DrawGammaSetMarkerTGraphAsym(graphIncRatioFitPi0SysErr , 20, 2, colorCent, colorCent, 1, kTRUE);	
	graphIncRatioFitPi0SysErr->Draw("E2same");
	histoIncRatioPi0Fit->DrawCopy("same");

	TLegend* legendIncRatioFit = GetAndSetLegend(0.12,0.2,1.5);
	legendIncRatioFit->AddEntry(graphIncRatioFitPi0SysErr,Form("PCM, %s",collisionSystem.Data()),"pf");
	legendIncRatioFit->Draw();

	canvasIncRatioPi0Fit->Print(Form("%s/IncRatioPi0Fitted_%s.eps",outputDir.Data(),cent.Data()));


	//************************************************************************
	//******************* Calculate error graph for pi0 **********************
	//************************************************************************	
	TGraphAsymmErrors* graphPi0FitSysErr 		= CalculateSysErrFromRelSysHisto( histoPi0Fit, "Pi0FitSystError", relSystErrorPi0FitDown, relSystErrorPi0FitUp, 2,
																				  nPointsErrors);	
	TGraphAsymmErrors* graphPi0SysErr 			= CalculateSysErrFromRelSysHisto( histoPi0, "Pi0SystError", relSystErrorPi0Down, relSystErrorPi0Up, 2,
																				  nPointsErrors);	
	 
	// **************************************************************************
	// ******************** plotting pi0 spectrum *******************************
	// **************************************************************************
	TCanvas *canvasPi0 = GetAndSetCanvas("canvasPi0Final");
	canvasPi0->SetLeftMargin(0.14);
	canvasPi0->SetLogy();
	canvasPi0->SetLogx();

	TH2D *dummyPi0 ;
	dummyPi0 = new TH2D("dummyPi0", "dummyPi0", 120, 0., 16, 1000., histoPi0->GetBinContent(maxBin)/100,histoPi0->GetBinContent(2)*10);
	SetStyleHistoTH2ForGraphs( dummyPi0, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}N_{#pi^{0}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
							   0.045, 0.05, 0.045, 0.05, 0.85, 1.2);
	dummyPi0->GetXaxis()->SetLabelOffset(-0.015);
	dummyPi0->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
	dummyPi0->DrawCopy();


	DrawGammaSetMarker(histoPi0, 20, 2., colorCentNotPi0Fitted, colorCentNotPi0Fitted);                
	DrawGammaSetMarkerTGraphAsym(graphPi0SysErr , 20, 2, colorCentNotPi0Fitted, colorCentNotPi0Fitted, 1, kTRUE);	
	graphPi0SysErr->Draw("E2same");
	histoPi0->DrawCopy("same");

	TLegend* legendPi0 = GetAndSetLegend(0.18,0.2,1.5);
	legendPi0->AddEntry(graphPi0SysErr,Form("#pi^{0} PCM, %s",collisionSystem.Data()),"pf");
	legendPi0->Draw();

	canvasPi0->Print(Form("%s/Pi0Spectrum_%s.eps",outputDir.Data(),cent.Data()));

	// **************************************************************************
	// ********************** plotting fiited pi0 spectrum **********************
	// **************************************************************************
	dummyPi0->DrawCopy();

	DrawGammaSetMarker(histoPi0Fit, 20, 2., colorCent, colorCent);                
	DrawGammaSetMarkerTGraphAsym(graphPi0FitSysErr , 20, 2, colorCent, colorCent, 1, kTRUE);	
	graphPi0FitSysErr->Draw("E2same");
	histoPi0Fit->DrawCopy("same");

	TLegend* legendPi0Fit = GetAndSetLegend(0.18,0.2,1.5);
	legendPi0Fit->AddEntry(graphPi0FitSysErr,Form("#pi^{0} PCM, %s",collisionSystem.Data()),"pf");

	legendPi0Fit->Draw();

	canvasPi0->Print(Form("%s/Pi0SprectrumFitted_%s.eps",outputDir.Data(),cent.Data()));

	
	// **************************************************************************
	// ***************** Calculate systematic error graphs **********************
	// **************************************************************************
	TGraphAsymmErrors* graphIncGammaSysErr 		= CalculateSysErrFromRelSysHisto( histoIncGammaSpec, "IncGammaSystError", relSystErrorIncGammaDown, relSystErrorIncGammaUp, 2,
																				  nPointsErrors);	
	TGraphAsymmErrors* graphIncGammaSysErrA 	= CalculateSysErrFromRelSysHisto( histoIncGammaSpec, "IncGammaSystErrorA", relSystErrorAIncGammaDown, relSystErrorAIncGammaUp, 2,
																				  nPointsErrors);	
	TGraphAsymmErrors* graphIncGammaSysErrB 	= CalculateSysErrFromRelSysHisto( histoIncGammaSpec, "IncGammaSystErrorB", relSystErrorBIncGammaDown, relSystErrorBIncGammaUp, 2,
																				  nPointsErrors);	
	TGraphAsymmErrors* graphIncGammaSysErrC 	= CalculateSysErrFromRelSysHisto( histoIncGammaSpec, "IncGammaSystErrorC", relSystErrorCIncGammaDown, relSystErrorCIncGammaUp, 2,
																				  nPointsErrors);	
	TGraphAsymmErrors* graphIncGammaSysErrW0Mat = CalculateSysErrFromRelSysHisto( histoIncGammaSpec, "IncGammaSystErrorC", relSystErrorWOMaterialIncGammaDown, relSystErrorWOMaterialIncGammaUp, 2,
																				  nPointsErrors);	
	
	
	// **************************************************************************
	// ******************* Plotting inclusive Gamma Spectrum ********************
	// **************************************************************************
	TCanvas *canvasIncGamma = GetAndSetCanvas("canvasIncGamma");
	canvasIncGamma->SetLeftMargin(0.14);
	canvasIncGamma->SetLogy();
	canvasIncGamma->SetLogx();

	TH2D *dummyGamma ;
	dummyGamma = new TH2D("dummyGamma", "dummyGamma", 120, 0., 16, 1000., histoIncGammaSpec->GetBinContent(maxBin)/100,histoIncGammaSpec->GetBinContent(2)*10);
	SetStyleHistoTH2ForGraphs( dummyGamma, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}N_{#gamma}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
							   0.045, 0.05, 0.045, 0.05, 0.85, 1.2);
	dummyGamma->GetXaxis()->SetLabelOffset(-0.015);
	dummyGamma->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
	dummyGamma->DrawCopy();


	DrawGammaSetMarker(histoIncGammaSpec, 20, 2., colorCent, colorCent);                
	DrawGammaSetMarkerTGraphAsym(graphIncGammaSysErr , 20, 2, colorCent, colorCent, 1, kTRUE);	
	graphIncGammaSysErr->Draw("E2same");
	histoIncGammaSpec->DrawCopy("same");

	TLegend* legendIncGamma = GetAndSetLegend(0.18,0.2,1.5);
	legendIncGamma->AddEntry(graphIncGammaSysErr,Form("#gamma_{inc} PCM, %s",collisionSystem.Data()),"pf");
	legendIncGamma->Draw();

	canvasIncGamma->Print(Form("%s/IncGammaSpectrum_%s.eps",outputDir.Data(),cent.Data()));
	
	// **************************************************************************
	// ***************** Calculate direct photon spectrum ***********************
	// **************************************************************************

	//_______________________ copy inclusive photon spectra _____________________
	TH1D *histoDirGammaSpectrumError 				= (TH1D*)histoIncGammaSpec->Clone("DirectPhotonSpectrum");
	TH1D *histoDirGammaSpectrumErrorSummed 			= (TH1D*)histoIncGammaSpec->Clone("DirectPhotonSpectrumSummed");
	TH1D *histoDirGammaSpectrumErrorSyst 			= (TH1D*)histoIncGammaSpec->Clone("DirectPhotonSpectrumSyst");
	TH1D *histoDirGammaSpectrumErrorSystA 			= (TH1D*)histoIncGammaSpec->Clone("DirectPhotonSpectrumSystA");
	TH1D *histoDirGammaSpectrumErrorSystB 			= (TH1D*)histoIncGammaSpec->Clone("DirectPhotonSpectrumSystB");
	TH1D *histoDirGammaSpectrumErrorSystC 			= (TH1D*)histoIncGammaSpec->Clone("DirectPhotonSpectrumSystC");
	TH1D *histoDirGammaSpectrumErrorSystWithoutMat 	= (TH1D*)histoIncGammaSpec->Clone("DirectPhotonSpectrumSyst");
	TH1D *histoDirGammaSpectrumErrorStat 			= (TH1D*)histoIncGammaSpec->Clone("DirectPhotonSpectrumStat");

	//_______________________ get arrays of double ratio errors __________________
	Double_t *systErrorsDoubleRatio 				= new Double_t[graphDoubleRatioFitPi0SysErr->GetN()];
	Double_t *systErrorsDoubleRatioA 				= new Double_t[graphDoubleRatioFitPi0SysErr->GetN()];
	Double_t *systErrorsDoubleRatioB 				= new Double_t[graphDoubleRatioFitPi0SysErr->GetN()];
	Double_t *systErrorsDoubleRatioC 				= new Double_t[graphDoubleRatioFitPi0SysErr->GetN()];
	Double_t *systErrorsDoubleRatioX 				= new Double_t[graphDoubleRatioFitPi0SysErr->GetN()];
	Double_t *systErrorsWithoutMatDoubleRatio 		= new Double_t[graphDoubleRatioFitPi0SysErr->GetN()];
	graphDoubleRatioFitPi0SysErr->Print();
	
	for (Int_t i = 0; i< graphDoubleRatioFitPi0SysErr->GetN(); i++){
		systErrorsDoubleRatio[i] 					= relSystErrorDoubleRatioPi0FitUp[i];
		systErrorsDoubleRatioA[i] 					= relSystErrorADoubleRatioPi0FitUp[i];
		systErrorsDoubleRatioB[i] 					= relSystErrorBDoubleRatioPi0FitUp[i];
		systErrorsDoubleRatioC[i] 					= relSystErrorCDoubleRatioPi0FitUp[i];
		systErrorsWithoutMatDoubleRatio[i] 			= relSystErrorWOMaterialDoubleRatioPi0FitUp[i];
	}
	
	systErrorsDoubleRatioX =  graphDoubleRatioFitPi0SysErr->GetX();

	//_______________________ copy inclusive photon spectra _____________________	
	TH1D* histoDoubleRatioWithSummedErrors 			= (TH1D*) histoDoubleRatioPi0Fit->Clone("DoubleRatioWithSummedErrors");
	TH1D* histoDoubleRatioWithStatErrors 			= (TH1D*) histoDoubleRatioPi0Fit->Clone("DoubleRatioWithStatErrors");
	TH1D* histoDoubleRatioWithSystErrors 			= (TH1D*) histoDoubleRatioPi0Fit->Clone("DoubleRatioWithSystErrors");
	TH1D* histoDoubleRatioWithSystErrorsA 			= (TH1D*) histoDoubleRatioPi0Fit->Clone("DoubleRatioWithSystErrorsA");
	TH1D* histoDoubleRatioWithSystErrorsB 			= (TH1D*) histoDoubleRatioPi0Fit->Clone("DoubleRatioWithSystErrorsB");
	TH1D* histoDoubleRatioWithSystErrorsC 			= (TH1D*) histoDoubleRatioPi0Fit->Clone("DoubleRatioWithSystErrorsC");
	TH1D* histoDoubleRatioWithSystErrorsWithoutMat 	= (TH1D*) histoDoubleRatioPi0Fit->Clone("DoubleRatioWithSystErrors");
	
	for(Int_t i = 1; i<histoDoubleRatioWithSummedErrors->GetNbinsX();i++){
// 		cout<<systErrorsDoubleRatioX[i-1]<<"  "<<histoDoubleRatioWithSummedErrors->GetBinCenter(i+1)<<endl;
		Double_t binErrorSummed = sqrt( pow( (histoDoubleRatioWithSummedErrors->GetBinError(i+1)/histoDoubleRatioWithSummedErrors->GetBinContent(i+1))*100,2) + pow(systErrorsDoubleRatio[i-1],2) );
		Double_t binErrorSyst = systErrorsDoubleRatio[i-1];
		Double_t binErrorSystA = systErrorsDoubleRatioA[i-1];
		Double_t binErrorSystB = systErrorsDoubleRatioB[i-1];
		Double_t binErrorSystC = systErrorsDoubleRatioC[i-1];
		Double_t binErrorSystWitoutMat = systErrorsWithoutMatDoubleRatio[i-1];
		Double_t binErrorStat = (histoDoubleRatioWithSummedErrors->GetBinError(i+1)/histoDoubleRatioWithSummedErrors->GetBinContent(i+1))*100;

		histoDoubleRatioWithSummedErrors->SetBinError(i+1,(binErrorSummed/100)*histoDoubleRatioWithSummedErrors->GetBinContent(i+1));
		histoDoubleRatioWithSystErrors->SetBinError(i+1,(binErrorSyst/100)*histoDoubleRatioWithSystErrors->GetBinContent(i+1));
		histoDoubleRatioWithSystErrorsA->SetBinError(i+1,(binErrorSystA/100)*histoDoubleRatioWithSystErrorsA->GetBinContent(i+1));
		histoDoubleRatioWithSystErrorsB->SetBinError(i+1,(binErrorSystB/100)*histoDoubleRatioWithSystErrorsB->GetBinContent(i+1));
		histoDoubleRatioWithSystErrorsC->SetBinError(i+1,(binErrorSystC/100)*histoDoubleRatioWithSystErrorsC->GetBinContent(i+1));
		histoDoubleRatioWithSystErrorsWithoutMat->SetBinError(i+1,(binErrorSystWitoutMat/100)*histoDoubleRatioWithSystErrorsWithoutMat->GetBinContent(i+1));
		histoDoubleRatioWithStatErrors->SetBinError(i+1,(binErrorStat/100)*histoDoubleRatioWithStatErrors->GetBinContent(i+1));
	}
   
	for(Int_t i = 1; i<histoIncGammaSpec->GetNbinsX();i++){
		histoDirGammaSpectrumError->SetBinContent(i+1,-1);
		histoDirGammaSpectrumErrorSummed->SetBinContent(i+1,-1);
		histoDirGammaSpectrumErrorSyst->SetBinContent(i+1,-1);
		histoDirGammaSpectrumErrorSystA->SetBinContent(i+1,-1);
		histoDirGammaSpectrumErrorSystB->SetBinContent(i+1,-1);
		histoDirGammaSpectrumErrorSystC->SetBinContent(i+1,-1);
		histoDirGammaSpectrumErrorStat->SetBinContent(i+1,-1);

		histoDirGammaSpectrumError->SetBinError(i+1,0);
		histoDirGammaSpectrumErrorSummed->SetBinError(i+1,0);
		histoDirGammaSpectrumErrorSyst->SetBinError(i+1,0);
		histoDirGammaSpectrumErrorSystA->SetBinError(i+1,0);
		histoDirGammaSpectrumErrorSystB->SetBinError(i+1,0);
		histoDirGammaSpectrumErrorSystC->SetBinError(i+1,0);
		histoDirGammaSpectrumErrorStat->SetBinError(i+1,0);
	}

	// get the binning of the direct photons from the DR
	TH1D *histoDirGammaSpectrumSyst 				= new TH1D(*histoDoubleRatioWithSystErrors);
	TH1D *histoDirGammaSpectrumStat 				= new TH1D(*histoDoubleRatioWithStatErrors);
	TH1D *histoDirGammaSpectrumSummed 				= new TH1D(*histoDoubleRatioWithSummedErrors);
	TH1D *histoDirGammaSpectrumSystWithoutMat	 	= new TH1D(*histoDoubleRatioWithSystErrorsWithoutMat);

	for(Int_t i = 1; i<histoDoubleRatioWithSystErrors->GetNbinsX(); i++){
		// obtain common quantities
		Double_t Rgamma 	= histoDoubleRatioWithSystErrors->GetBinContent(i+1);
		Double_t nIncGamma	= graphIncGammaSysErr->GetY()[i-1];
		
		// calculating systematics graph
		Double_t errRgamma	= histoDoubleRatioWithSystErrors->GetBinError(i+1);
		Double_t errNIncGam = graphIncGammaSysErr->GetEYhigh()[i-1];
		Double_t q1 		= 1 - 1/ Rgamma;
		
		Double_t q1Error 	= errRgamma/(Rgamma*Rgamma);
		Double_t content 	= nIncGamma * ( 1 - 1/ Rgamma);
		Double_t error 		= sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
		Double_t errDR		= content - error;
		histoDirGammaSpectrumSyst->SetBinError(i+1, error);
		histoDirGammaSpectrumSyst->SetBinContent(i+1, content);
		histoDirGammaSpectrumErrorSyst->SetBinContent(i+1, errDR);
		
		// calculating stat graphs
		errRgamma	= histoDoubleRatioWithStatErrors->GetBinError(i+1);
		errNIncGam 	= histoIncGammaSpec->GetBinError(i+1);
		q1 			= 1 - 1/ Rgamma;
		q1Error 	= errRgamma/(Rgamma*Rgamma);
		content 	= nIncGamma * ( 1 - 1/ Rgamma);
		error 		= sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
		errDR		= content - error;
		histoDirGammaSpectrumStat->SetBinError(i+1, error);
		histoDirGammaSpectrumStat->SetBinContent(i+1, content);
		histoDirGammaSpectrumErrorStat->SetBinContent(i+1, errDR);
		
		// calculating summed error graphs
		errRgamma	= histoDirGammaSpectrumSummed->GetBinError(i+1);
		errNIncGam 	= sqrt( pow( histoIncGammaSpec->GetBinError(i+1),2) + pow( graphIncGammaSysErr->GetEYhigh()[i-1], 2) );
		q1 			= 1 - 1/ Rgamma;
		q1Error 	= errRgamma/(Rgamma*Rgamma);
		content 	= nIncGamma * ( 1 - 1/ Rgamma);
		error 		= sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
		errDR		= content - error;
		histoDirGammaSpectrumSummed->SetBinError(i+1, error);
		histoDirGammaSpectrumSummed->SetBinContent(i+1, content);
		histoDirGammaSpectrumErrorSummed->SetBinContent(i+1, errDR);

		// calculating sys err without error graph
		errRgamma	= histoDirGammaSpectrumErrorSystWithoutMat->GetBinError(i+1);
		errNIncGam 	= graphIncGammaSysErrW0Mat->GetEYhigh()[i-1];
		q1 			= 1 - 1/ Rgamma;
		q1Error 	= errRgamma/(Rgamma*Rgamma);
		content 	= nIncGamma * ( 1 - 1/ Rgamma);
		error 		= sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
		errDR		= content - error;
		histoDirGammaSpectrumSystWithoutMat->SetBinError(i+1, error);
		histoDirGammaSpectrumSystWithoutMat->SetBinContent(i+1, content);
		histoDirGammaSpectrumErrorSystWithoutMat->SetBinContent(i+1, errDR);
	}
	
	// purely calculating points based on all systematic errors
	TGraphAsymmErrors *graphDirGammaSpectrumSyst = CalculateDirectPhotonPointsAndUpperLimits(histoDirGammaSpectrumErrorSyst,histoDirGammaSpectrumStat,0,0.5);
	if(graphDirGammaSpectrumSyst)graphDirGammaSpectrumSyst->SetName("graphDirGammaSpectrumSyst");
	// purely calculating points based on statistical errors
	TGraphAsymmErrors *graphDirGammaSpectrumStat = CalculateDirectPhotonPointsAndUpperLimits(histoDirGammaSpectrumErrorStat,histoDirGammaSpectrumStat,0,0.5);
	if(graphDirGammaSpectrumStat)graphDirGammaSpectrumStat->SetName("graphDirGammaSpectrumStat");
	// purely calculating points based on all systematic + statistical errors
	TGraphAsymmErrors *graphDirGammaSpectrumSummed = CalculateDirectPhotonPointsAndUpperLimits(histoDirGammaSpectrumErrorSummed,histoDirGammaSpectrumStat,0,0.5);
	if(graphDirGammaSpectrumSummed)graphDirGammaSpectrumSummed->SetName("graphDirGammaSpectrumSummed");
	
	// calculate points above confidence level summed errors
	TGraphAsymmErrors *graphDirGammaSpectrumSummedConfi = CalculateDirectPhotonPointsAndUpperLimits(histoDirGammaSpectrumErrorSummed,histoDirGammaSpectrumStat,2,0.5);
	if(graphDirGammaSpectrumSummedConfi)graphDirGammaSpectrumSummedConfi->SetName("graphDirGammaSpectrumSummedConfi");
	// calculate upperlimits summed errors
	TGraphAsymmErrors *graphDirGammaSpectrumSummedUL = CalculateDirectPhotonPointsAndUpperLimits(histoDirGammaSpectrumErrorSummed,histoDirGammaSpectrumStat,4,0.5);
	if(graphDirGammaSpectrumSummedUL)graphDirGammaSpectrumSummedUL->SetName("graphDirGammaSpectrumSummedUL");
	// calculate arrows for points with 0, error summed
	TGraphAsymmErrors *graphDirGammaSpectrumSummedAr = CalculateDirectPhotonPointsAndUpperLimits(histoDirGammaSpectrumErrorSummed,histoDirGammaSpectrumStat,5,0.5);
	if(graphDirGammaSpectrumSummedAr)graphDirGammaSpectrumSummedAr->SetName("graphDirGammaSpectrumSummedAr");
	// calculate points below confidence level summed errors
	TGraphAsymmErrors *graphDirGammaSpectrumSummedULConfi = CalculateDirectPhotonPointsAndUpperLimits(histoDirGammaSpectrumErrorSummed,histoDirGammaSpectrumStat,6,0.5);
	if(graphDirGammaSpectrumSummedULConfi)graphDirGammaSpectrumSummedULConfi->SetName("graphDirGammaSpectrumSummedULConfi");
	// calculate points below confidence level summed errors with arrows	
	TGraphAsymmErrors *graphDirGammaSpectrumSummedArConfi = CalculateDirectPhotonPointsAndUpperLimits(histoDirGammaSpectrumErrorSummed,histoDirGammaSpectrumStat,7,0.5);
	if(graphDirGammaSpectrumSummedArConfi)graphDirGammaSpectrumSummedArConfi->SetName("graphDirGammaSpectrumSummedAr");
	// purely calculating points based on systematic errors without material budget
	TGraphAsymmErrors *graphDirGammaSpectrumSystwoMat = CalculateDirectPhotonPointsAndUpperLimits(histoDirGammaSpectrumErrorSystWithoutMat,histoDirGammaSpectrumStat,0,0.5);
	if(graphDirGammaSpectrumSystwoMat) graphDirGammaSpectrumSystwoMat->SetName("graphDirGammaSpectrumSystWOMat");


	// style setting of direct photon graphs
	DrawGammaSetMarkerTGraphAsym(graphDirGammaSpectrumSyst , 20,2, colorCent, colorCent, 1, kTRUE);
	graphDirGammaSpectrumSyst->Draw("Z2");
	DrawGammaSetMarkerTGraphAsym(graphDirGammaSpectrumStat , 20, 2, colorCent, colorCent, 1, kTRUE);
	graphDirGammaSpectrumStat->Draw("pE1");
	DrawGammaSetMarkerTGraphAsym(graphDirGammaSpectrumSummedUL , 20, 0, colorCent, colorCent, 1, kTRUE);
	graphDirGammaSpectrumSummedUL->Draw("||");
	DrawGammaSetMarkerTGraphAsym(graphDirGammaSpectrumSummedAr , 1, 0, colorCent, colorCent, 1, kTRUE);
	graphDirGammaSpectrumSummedAr->Draw("|>");

	// **************************************************************************
	// ******************* Plotting direct Gamma Spectrum ***********************
	// **************************************************************************
	TCanvas *canvasDirGamma = GetAndSetCanvas("canvasDirGamma");
	canvasDirGamma->SetLeftMargin(0.14);
	canvasDirGamma->SetLogy();
	canvasDirGamma->SetLogx();

	TH2D *dummyDirGamma ;
	dummyDirGamma = new TH2D("dummyDirGamma", "dummyDirGamma", 120, 0., 16, 1000., histoIncGammaSpec->GetBinContent(maxBin)/1000,histoIncGammaSpec->GetBinContent(2)*10e-1);
	SetStyleHistoTH2ForGraphs( dummyDirGamma, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}N_{#gamma}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
							   0.045, 0.05, 0.045, 0.05, 0.85, 1.2);
	dummyDirGamma->GetXaxis()->SetLabelOffset(-0.015);
	dummyDirGamma->GetXaxis()->SetRangeUser(doubleRatioX[0],doubleRatioX[1]);
	dummyDirGamma->DrawCopy();

	graphDirGammaSpectrumSyst->Draw("Z2,same");
	graphDirGammaSpectrumStat->Draw("p,same,e1");

	TLegend* legendDirGamma = GetAndSetLegend(0.18,0.2,1.5);
	legendDirGamma->AddEntry(graphDirGammaSpectrumSyst,Form("#gamma_{dir} PCM, %s",collisionSystem.Data()),"pf");
	legendDirGamma->Draw();

	canvasDirGamma->Print(Form("%s/DirGammaSpectrum_%s.eps",outputDir.Data(),cent.Data()));
	
	
	// **************************************************************************
	// ****************** Write data to output file *****************************
	// **************************************************************************
	TString optionOutput = option;
	if (option.CompareTo("2.76TeV") == 0 || option.CompareTo("7TeV") == 0 ) optionOutput = "pp";
	const char* fileNameOutputComp = Form("Gamma_PCMResults_%s.root",optionOutput.Data());
	TFile* fileGammaSpectrum = new TFile(fileNameOutputComp,"UPDATE");		
		fileGammaSpectrum->mkdir(Form("Gamma_%s_%s",option.Data(),GetCentralityString(cutSel).Data()));
		fileGammaSpectrum->cd(Form("Gamma_%s_%s", option.Data(),GetCentralityString(cutSel).Data()));
		
			histoDoubleRatio->SetName("DoubleRatioStatError");
			SetHistogramm(histoDoubleRatio,"#it{p}_{T} (GeV/#it{c})", "(#it{N}_{#gamma_{inc}}/#it{N}_{#pi^{0}})/(#it{N}_{#gamma_{decay}}/#it{N}_{#pi^{0}})");
			histoDoubleRatio->Write("DoubleRatioStatError",TObject::kOverwrite);
			graphDoubleRatioSysErr->Write("DoubleRatioSystError",TObject::kOverwrite);
			graphDoubleRatioSysErrA->Write("DoubleRatioSystErrorA",TObject::kOverwrite);
			graphDoubleRatioSysErrB->Write("DoubleRatioSystErrorB",TObject::kOverwrite);
			graphDoubleRatioSysErrC->Write("DoubleRatioSystErrorC",TObject::kOverwrite);
			histoDoubleRatioPi0Fit->SetName("DoubleRatioPi0FitStatError");
			SetHistogramm(histoDoubleRatioPi0Fit,"#it{p}_{T} (GeV/#it{c})", "(#it{N}_{#gamma_{inc}}/#it{N}_{#pi^{0}})/(#it{N}_{#gamma_{decay}}/#it{N}_{#pi^{0}})");
			histoDoubleRatioPi0Fit->Write("DoubleRatioPi0FitStatError",TObject::kOverwrite);
			graphDoubleRatioFitPi0SysErr->Write("DoubleRatioPi0FitSystError",TObject::kOverwrite);
			graphDoubleRatioFitPi0SysErrA->Write("DoubleRatioPi0FitSystErrorA",TObject::kOverwrite);
			graphDoubleRatioFitPi0SysErrB->Write("DoubleRatioPi0FitSystErrorB",TObject::kOverwrite);
			graphDoubleRatioFitPi0SysErrC->Write("DoubleRatioPi0FitSystErrorC",TObject::kOverwrite);

			histoIncRatio->SetName("IncRatioStatError");
			SetHistogramm(histoIncRatio,"#it{p}_{T} (GeV/#it{c})", "#gamma_{inc}/#pi^{0}");
			histoIncRatio->Write("IncRatioStatError",TObject::kOverwrite);
			graphIncRatioSysErr->Write("IncRatioSystError",TObject::kOverwrite);
			graphIncRatioSysErrA->Write("IncRatioSystErrorA",TObject::kOverwrite);
			graphIncRatioSysErrB->Write("IncRatioSystErrorB",TObject::kOverwrite);
			graphIncRatioSysErrC->Write("IncRatioSystErrorC",TObject::kOverwrite);
			histoIncRatioPi0Fit->SetName("IncRatioPi0FitStatError");
			SetHistogramm(histoIncRatioPi0Fit,"#it{p}_{T} (GeV/#it{c})", "#gamma_{inc}/#pi^{0}");
			histoIncRatioPi0Fit->Write("IncRatioPi0FitStatError",TObject::kOverwrite);
			graphIncRatioFitPi0SysErr->Write("IncRatioPi0FitSystError",TObject::kOverwrite);
			graphIncRatioFitPi0SysErrA->Write("IncRatioPi0FitSystErrorA",TObject::kOverwrite);
			graphIncRatioFitPi0SysErrB->Write("IncRatioPi0FitSystErrorB",TObject::kOverwrite);
			graphIncRatioFitPi0SysErrC->Write("IncRatioPi0FitSystErrorC",TObject::kOverwrite);

			histoIncGammaSpec->SetName("IncGammaStatError");
			SetHistogramm(histoIncGammaSpec,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}N_{#gamma}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})");
			histoIncGammaSpec->Write("IncGammaStatError",TObject::kOverwrite);
			graphIncGammaSysErr->Write("IncGammaSystError",TObject::kOverwrite);
			graphIncGammaSysErrA->Write("IncGammaSystErrorA",TObject::kOverwrite);
			graphIncGammaSysErrB->Write("IncGammaSystErrorB",TObject::kOverwrite);
			graphIncGammaSysErrC->Write("IncGammaSystErrorC",TObject::kOverwrite);

			histoPi0->SetName("Pi0StatError");
			SetHistogramm(histoPi0,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}N_{#pi^{0}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})");
			histoPi0->Write("Pi0StatError",TObject::kOverwrite);
			graphPi0SysErr->Write("Pi0SystError",TObject::kOverwrite);
			histoPi0Fit->SetName("Pi0FitStatError");
			SetHistogramm(histoPi0Fit,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}N_{#pi^{0}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})");
			histoPi0Fit->Write("Pi0FitStatError",TObject::kOverwrite);
			graphPi0FitSysErr->Write("Pi0FitSystError",TObject::kOverwrite);

			if(graphDirGammaSpectrumSyst)			graphDirGammaSpectrumSyst->Write("graphDirGammaSpectrumSyst",TObject::kOverwrite);
			if(graphDirGammaSpectrumStat)			graphDirGammaSpectrumStat->Write("graphDirGammaSpectrumStat",TObject::kOverwrite);
			if(graphDirGammaSpectrumSummed)			graphDirGammaSpectrumSummed->Write("graphDirGammaSpectrumSummed",TObject::kOverwrite);
			if(graphDirGammaSpectrumSummedConfi)	graphDirGammaSpectrumSummedConfi->Write("graphDirGammaSpectrumSummedConfi",TObject::kOverwrite);
			if(graphDirGammaSpectrumSummedUL)		graphDirGammaSpectrumSummedUL->Write("graphDirGammaSpectrumSummedUL",TObject::kOverwrite);
			if(graphDirGammaSpectrumSummedAr)		graphDirGammaSpectrumSummedAr->Write("graphDirGammaSpectrumSummedAr",TObject::kOverwrite);
			if(graphDirGammaSpectrumSummedULConfi)	graphDirGammaSpectrumSummedULConfi->Write("graphDirGammaSpectrumSummedULConfi",TObject::kOverwrite);
			if(graphDirGammaSpectrumSummedArConfi)	graphDirGammaSpectrumSummedArConfi->Write("graphDirGammaSpectrumSummedAr",TObject::kOverwrite);
			if(graphDirGammaSpectrumSystwoMat) 		graphDirGammaSpectrumSystwoMat->Write("graphDirGammaSpectrumSystWOMat",TObject::kOverwrite);
				
			if (histoCocktailAllGamma)histoCocktailAllGamma->Write("CocktailSumGamma",TObject::kOverwrite);
			if (histoCocktailPi0Gamma)histoCocktailPi0Gamma->Write("CocktailPi0Gamma",TObject::kOverwrite);
			if (histoCocktailEtaGamma)histoCocktailEtaGamma->Write("CocktailEtaGamma",TObject::kOverwrite);
			if (histoCocktailOmegaGamma)histoCocktailOmegaGamma->Write("CocktailOmegaGamma",TObject::kOverwrite);
			if (histoCocktailEtapGamma)histoCocktailEtapGamma->Write("CocktailEtapGamma",TObject::kOverwrite);
			if (histoCocktailPhiGamma)histoCocktailPhiGamma->Write("CocktailPhiGamma",TObject::kOverwrite);
			if (histoCocktailRhoGamma)histoCocktailRhoGamma->Write("CocktailRhoGamma",TObject::kOverwrite);
			if (histoCocktailSigmaGamma)histoCocktailSigmaGamma->Write("CocktailSigmaGamma",TObject::kOverwrite);
	
			
			fileGammaSpectrum->mkdir("Theory");
			fileGammaSpectrum->cd("Theory");
				if (optionOutput.CompareTo("pp") == 0 ) NLODoubleRatio->Write(Form("NLODoubleRatio_%s_%s",option.Data(),GetCentralityString(cutSel).Data()),TObject::kOverwrite);
					else NLODoubleRatio->Write(Form("NLODoubleRatio_%s",GetCentralityString(cutSel).Data()),TObject::kOverwrite);
				if (optionOutput.CompareTo("pp") == 0 ) NLO->Write(Form("NLO_%s_%s",option.Data(),GetCentralityString(cutSel).Data()),TObject::kOverwrite);
					else NLO->Write(Form("NLO_%s",GetCentralityString(cutSel).Data()),TObject::kOverwrite);
				if (graphInvYieldPbPbTheoryCTEQ61EPS09) graphInvYieldPbPbTheoryCTEQ61EPS09->Write(Form("EPS09_%s",GetCentralityString(cutSel).Data()),TObject::kOverwrite);
				if (graphInvYieldPPTheoryCT10BFG2_pdfErr) graphInvYieldPPTheoryCT10BFG2_pdfErr->Write(Form("CT10BF_%s",GetCentralityString(cutSel).Data()),TObject::kOverwrite);
				if (graphDRPbPbCTEQ61EPS09) graphDRPbPbCTEQ61EPS09->Write(Form("EPS09DoubleRatio_%s",GetCentralityString(cutSel).Data()),TObject::kOverwrite);
				if (graphDRPbPbCT10BFG2_pdfErr) graphDRPbPbCT10BFG2_pdfErr->Write(Form("CT10BFDoubleRatio_%s",GetCentralityString(cutSel).Data()),TObject::kOverwrite);
	fileGammaSpectrum->Write();
	fileGammaSpectrum->Close();


}

