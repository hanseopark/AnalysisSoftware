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
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"

extern TRandom*	gRandom;
extern TBenchmark*	gBenchmark;
extern TSystem*	gSystem;
extern TMinuit*  	gMinuit;	

struct SysErrorConversion {
	Double_t value;
	Double_t error;
	//	TString name;
};

void CombineMesonMeasurements7TeV( 	TString fileNamePCM = "", 
									TString fileNameEMCAL = "",  
									TString suffix = "eps", 
									TString isMC= "", 
									TString thesisPlots = "", 
									TString bWCorrection="X"){	

	TString date = ReturnDateString();
	
	gROOT->Reset();	
	gROOT->SetStyle("Plain");
	
	StyleSettingsThesis();	
	SetPlotStyle();
	
	TString dateForOutput 			= ReturnDateStringForOutput();
	cout << dateForOutput.Data() << endl;
	//___________________________________ Declaration of files _____________________________________________
	TString collisionSystem7TeV 			= "pp, #sqrt{#it{s}} = 7 TeV";		
	
	TString fileNameTheory					= "ExternalInput/TheoryCompilationPP.root";	
	TString fileNamePHOS					= "ExternalInput/PHOS/7TeV/PHOS_pi0_7TeV_20111030_BWcorr.root";
	TString fileNameChargedPionPP 			= "ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_20_May_2015.root";
	TString fileNameChargedHadronPP			= "ExternalInput/UnidentifiedCharged/ChargedHadrinSpectraPP_20_May_2015.root";
	TString outputDir 						= Form("%s/%s/CombineMesonMeasurements7TeV%s",suffix.Data(),dateForOutput.Data(),bWCorrection.Data());
	TString nameFinalResDat 				= Form("%s/CombinedResults%s_FitResults.dat",outputDir.Data(),bWCorrection.Data());
	cout << outputDir.Data() << endl;
	cout << fileNamePCM.Data() << endl;

	gSystem->Exec("mkdir -p "+outputDir);
 	gSystem->Exec(Form("cp %s %s/InputPCM.root", fileNamePCM.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/InputPHOS.root", fileNamePHOS.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/InputEMCALLow.root", fileNameEMCAL.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/Theory.root", fileNameTheory.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/ChargedPionsPP.root", fileNameChargedPionPP.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/ChargedHadronsPP.root", fileNameChargedHadronPP.Data(), outputDir.Data()));
	
	Bool_t thesis 							= kFALSE;
	if(thesisPlots.CompareTo("thesis") == 0){
		thesis 								= kTRUE;
	}
	
	TString prefix2							= "";	
	if (isMC.CompareTo("kTRUE")==0){ 
		prefix2 = "MC";
	} else {	
		prefix2 = "Data";
	}
	
	Double_t mesonMassExpectPi0 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
	
	Double_t xSection7TeV 			= 62.22*1e-3;
	Double_t xSection7TeVV0AND 		= 54.31*1e-3;
	Double_t xSection7TeVErrUp 		= 2.18;
	Double_t xSection7TeVErrDown 	= 2.18;
	Double_t xSection7TeVppINEL 	= 73.2*1e9;	
	Double_t recalcBarn 			= 1e12; //NLO in pbarn!!!!

	Width_t		widthLinesBoxes		= 1.4;
	Width_t		widthCommonFit		= 2;
	
	// Definition of colors, styles and markers sizes
	Color_t		colorComb			= kBlue+2;
	Style_t		markerStyleComb		= 20;
	Size_t		markerSizeComb		= 2;
	
	Color_t 	colorCombLowPt 			= GetDefaultColorDiffDetectors("Comb", kFALSE, kFALSE, kFALSE);
	Color_t 	colorCombHighPt 		= GetDefaultColorDiffDetectors("Comb", kFALSE, kFALSE, kTRUE);
	Style_t 	markerStyleCombLowPt	= 20;
	Style_t 	markerStyleCombHighPt	= 20;	
	Size_t 		markerSizeComparison = 2;
	
	TString 	nameMeasGlobal[11] 	= {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMCal", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged", "PCMOtherDataset"};
	Color_t 	colorDet[11];
	Color_t 	colorDetMC[11];
	Style_t 	markerStyleDet[11];
	Style_t 	markerStyleDetMC[11];
	Size_t 		markerSizeDet[11];
	Size_t 		markerSizeDetMC[11];

	Style_t 	styleMarkerNLOMuHalf	= 24;
	Style_t 	styleMarkerNLOMuOne		= 27;
	Style_t 	styleMarkerNLOMuTwo		= 30;
	Style_t 	styleLineNLOMuHalf		= 8;
	Style_t 	styleLineNLOMuOne		= 7;
	Style_t 	styleLineNLOMuTwo		= 4;
	Style_t 	styleLineNLOMuTwoBKK	= 3;
	Style_t 	styleLineNLOMuTwoDSS	= 6;
	Size_t		sizeMarkerNLO			= 1;
	Width_t		widthLineNLO			= 2.;
	
	for (Int_t i = 0; i < 11; i++){
		colorDet[i]					= GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kFALSE, kFALSE, kTRUE);
		colorDetMC[i]				= GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kTRUE, kFALSE, kTRUE);
		markerStyleDet[i]			= GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kFALSE);
		markerStyleDetMC[i]			= GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kTRUE);
		markerSizeDet[i]			= GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kFALSE)*2;
		markerSizeDetMC[i]			= GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kTRUE)*2;
	}	
	
	//************************** Read data for PCM **************************************************
	TFile* filePCM 									= new TFile(fileNamePCM.Data());
	TH1D* histoPCMNumberOfEvents7TeV 				= (TH1D*)filePCM->Get("histoNumberOfEvents7TeV");
	TDirectory* directoryPCMPi07TeV 				= (TDirectory*)filePCM->Get("Pi07TeV"); 

		TF1* fitCorrectionFactorsHistvsPt 								= new TF1("fitCorrectionFactorsHistvsPt","[0]/pow(x,[1])+[2]");
		fitCorrectionFactorsHistvsPt->SetParameter(0,2.9737546081);
		fitCorrectionFactorsHistvsPt->SetParameter(1,1.4795520406);
		fitCorrectionFactorsHistvsPt->SetParameter(2,2.2652589579);

		TH1D* histoPCMPi0InvCrossSection7TeV					= (TH1D*)directoryPCMPi07TeV->Get("InvCrossSectionPi0");
		for (Int_t i = 2; i < histoPCMPi0InvCrossSection7TeV->GetNbinsX(); i++){
			histoPCMPi0InvCrossSection7TeV->SetBinContent(i,histoPCMPi0InvCrossSection7TeV->GetBinContent(i)*(100-fitCorrectionFactorsHistvsPt->Eval(histoPCMPi0InvCrossSection7TeV->GetBinCenter(i)))/100);
			histoPCMPi0InvCrossSection7TeV->SetBinError(i,histoPCMPi0InvCrossSection7TeV->GetBinError(i)*(100-fitCorrectionFactorsHistvsPt->Eval(histoPCMPi0InvCrossSection7TeV->GetBinCenter(i)))/100);
		}   
		
		TGraphAsymmErrors* graphPCMPi0InvCrossSectionSysA7TeV	= (TGraphAsymmErrors*)directoryPCMPi07TeV->Get("InvCrossSectionPi0SysA");
		TGraphAsymmErrors* graphPCMPi0InvCrossSectionSys7TeV	= (TGraphAsymmErrors*)directoryPCMPi07TeV->Get("InvCrossSectionPi0Sys");
		cout << "PCM sys" << endl;
		graphPCMPi0InvCrossSectionSys7TeV->Print();

		Double_t* yValue 										= graphPCMPi0InvCrossSectionSys7TeV->GetY();
		Double_t* yErrorLow 									= graphPCMPi0InvCrossSectionSys7TeV->GetEYlow();
		Double_t* yErrorHigh 									= graphPCMPi0InvCrossSectionSys7TeV->GetEYhigh();
		Double_t* xLow 											= graphPCMPi0InvCrossSectionSys7TeV->GetX();
		for (Int_t i = 0; i < graphPCMPi0InvCrossSectionSys7TeV->GetN(); i++){
			cout << xLow[i] << "\t" << 100-fitCorrectionFactorsHistvsPt->Eval(xLow[i]) << endl;
			cout << yValue[i] << "\t" << yValue[i]*(100-fitCorrectionFactorsHistvsPt->Eval(xLow[i]))/100 << "\t" << yValue[i] << endl;
			yValue[i] 	= yValue[i]*(100-fitCorrectionFactorsHistvsPt->Eval(xLow[i]))/100;
			yErrorLow[i] = yErrorLow[i]*(100-fitCorrectionFactorsHistvsPt->Eval(xLow[i]))/100;
			yErrorHigh[i] = yErrorHigh[i]*(100-fitCorrectionFactorsHistvsPt->Eval(xLow[i]))/100;
		}   
		graphPCMPi0InvCrossSectionSys7TeV->Print();
		
		TGraphAsymmErrors* graphPCMPi0InvCrossSectionStat7TeV 	= new TGraphAsymmErrors(histoPCMPi0InvCrossSection7TeV);
		cout << "PCM stat" << endl;
		graphPCMPi0InvCrossSectionStat7TeV->RemovePoint(0);
		graphPCMPi0InvCrossSectionStat7TeV->Print();
		
	//************************** Read data for EMCAL ****************************************************
	TFile* fileEMCAL								= new TFile(fileNameEMCAL.Data());
		TH1D*	histoEMCALPi0InvCrossSection7TeV 				= (TH1D*)fileEMCAL->Get("pi0Stat");
		TGraphAsymmErrors* graphEMCALPi0InvCrossSectionStat7TeV = new TGraphAsymmErrors(histoEMCALPi0InvCrossSection7TeV);
		cout << "EMCAL stat" <<endl;
		graphEMCALPi0InvCrossSectionStat7TeV->Print();
		
		TH1D*	histoEMCALPi0InvCrossSectionSys7TeV 			= (TH1D*)fileEMCAL->Get("pi0Syst");
		TGraphAsymmErrors* graphEMCALPi0InvCrossSectionSys7TeV 	= new TGraphAsymmErrors(histoEMCALPi0InvCrossSectionSys7TeV);
		cout << "EMCAL sys" <<endl;
		graphEMCALPi0InvCrossSectionSys7TeV->Print();
		
	//************************** Read data for PHOS *****************************************************
	TFile* filePHOS									= new TFile(fileNamePHOS);
		TH1D* histoPHOSPi0InvCrossSection7TeV 					= (TH1D*)filePHOS->Get("hPi07TeVStat");
		TH1D* histoPHOSPi0InvCrossSectionSys7TeV 				= (TH1D*)filePHOS->Get("hPi07TeVSys");
		histoPHOSPi0InvCrossSection7TeV->Scale(xSection7TeV*recalcBarn*0.95);
		histoPHOSPi0InvCrossSectionSys7TeV->Scale(xSection7TeV*recalcBarn*0.95);
		TGraphAsymmErrors* graphPHOSPi0InvCrossSectionStat7TeV 	= new TGraphAsymmErrors(histoPHOSPi0InvCrossSection7TeV);
		cout << "PHOS stat" <<endl;
		graphPHOSPi0InvCrossSectionStat7TeV->Print();
		TGraphAsymmErrors* graphPHOSPi0InvCrossSectionSys7TeV 	= new TGraphAsymmErrors(histoPHOSPi0InvCrossSectionSys7TeV);
		cout << "PHOS sys" <<endl;
		graphPHOSPi0InvCrossSectionStat7TeV->Print();
		
	// *******************************************************************************************************
	// ************************** Loading theory calculations ************************************************
	// *******************************************************************************************************		
	TFile* fileTheoryCompilation = new TFile(fileNameTheory.Data());
		TH1F* histoPythia8InvXSection 				= (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Spec7TeV");
		TH1F* histoPythia8InvXSection_VarBinning 	= (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Spec7TeVVarBinning");
		TGraph* graphNLOCalcMuHalf7TeV			= (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuHalf7TeV");
		TGraph* graphNLOCalcMuOne7TeV			= (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuOne7TeV");
		TGraph* graphNLOCalcMuTwo7TeV			= (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuTwo7TeV");
		TGraph* graphNLOCalcEtaMuHalf7TeV		= (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuHalf7TeV");
		TGraph* graphNLOCalcEtaMuOne7TeV			= (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuOne7TeV");
		TGraph* graphNLOCalcEtaMuTwo7TeV			= (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuTwo7TeV");
		TGraph* graphEtaToPi0NLOMuHalf7TeV 		= (TGraph*)fileTheoryCompilation->Get("graphNLOCalcEtaOverPi0MuHalf7TeV");
		TGraph* graphEtaToPi0NLOMuOne7TeV 		= (TGraph*)fileTheoryCompilation->Get("graphNLOCalcEtaOverPi0MuOne7TeV");
		TGraph* graphEtaToPi0NLOMuTwo7TeV 		= (TGraph*)fileTheoryCompilation->Get("graphNLOCalcEtaOverPi0MuTwo7TeV");
		TGraph* graphNLODSSCalcMuTwo7TeV			= (TGraph*)fileTheoryCompilation->Get("graphNLOCalcDSSInvSecPi0MuTwo7TeV");
		TGraph* graphNLOCGCCalcMuTwo7TeV			= (TGraph*)fileTheoryCompilation->Get("graphNLOCalcCGCInvCrossSec7TeV");
 		TGraphAsymmErrors* graphNLODSS14Calc7TeV	= (TGraphAsymmErrors*)fileTheoryCompilation->Get("graphNLOCalcDSS14InvCrossSec7TeV");
	
	// *******************************************************************************************************
	// ************************** Loading of Charged pions ************************************************
	// *******************************************************************************************************			
	TFile* fileChargedPionInputpp 				= new TFile(fileNameChargedPionPP.Data());
	TH1D* histoChargedPionSpecHighPtStatPP 		= (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecHighPtStat7TeV");
	TH1D* histoChargedPionSpecHighPtSystPP 		= (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecHighPtSyst7TeV");
	   
	TH1D*	histoChargedPionSpecLowPtStat7TeVCMS 			= (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtStat7TeVCMS");
	TH1D*	histoChargedPionSpecLowPtSys7TeVCMS 				= (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtSys7TeVCMS");
	TH1D*	histoChargedPionSpecLowPtStatPP7TeV 				= (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtStatPP7TeV");
	TH1D*	histoChargedPionSpecLowPtSysPP7TeV 				= (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtSysPP7TeV");

	TH1D*	histoChargedPionSpecLowPtStat7TeVCMS 				= (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtStat7TeVCMS");
	TH1D*	histoChargedPionSpecLowPtSys7TeVCMS 				= (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtSys7TeVCMS");
	TH1D*	histoChargedPionSpecLowPtStatPP7TeV 				= (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtStat7TeVALICE");
	TH1D*	histoChargedPionSpecLowPtSysPP7TeV 					= (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtSys7TeVALICE");
	TH1D*	histoChargedPionSpecHighPtStatPP7TeV 				= (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecHighPtStat7TeVALICE");
	TGraphAsymmErrors*	 graphChargedPionSpecHighPtSystPP7TeV 	= (TGraphAsymmErrors*)fileChargedPionInputpp->Get("graphChargedPionSpecHighPtSys7TeVALICE");
		
		
	// *******************************************************************************************************
	// ************************** Combination of different measurements **************************************
	// *******************************************************************************************************
	// REMARKS: 
	// 		- order of measurements defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
	//		- correlations are defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
	// 		- currently only PCM-EMCAL vs others fully implemeted energy independent
	// 		- extendable to other energies
	//		- offsets have to be determined manually, see cout's in shell from combination function, more can be uncommented
		
	
	// definition of array of histograms (NULL - means we have no measurement at this energy for this rec-method)
	// for statistical error and final value from respective method
	TH1D* statErrorCollection[11];
	for (Int_t i = 0; i< 11; i++){
		statErrorCollection[i] = NULL;
	}	
	statErrorCollection[0] = (TH1D*)histoPCMPi0InvCrossSection7TeV->Clone("statErrPCMPi0");
	statErrorCollection[1] = (TH1D*)histoPHOSPi0InvCrossSection7TeV->Clone("statErrPHOSPi0");

	
	// definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)	
	// for systematic error from respective method
	TGraphAsymmErrors* sysErrorCollection[11];
	for (Int_t i = 0; i< 11; i++){
		sysErrorCollection[i] = NULL;
	}	
	sysErrorCollection[0] = (TGraphAsymmErrors*)graphPCMPi0InvCrossSectionSys7TeV->Clone("sysErrPCMPi0");
	sysErrorCollection[1] = (TGraphAsymmErrors*)graphPHOSPi0InvCrossSectionSys7TeV->Clone("sysErrPHOSPi0");

	
	// Definition of final pt binning (has to be set manually)
	Double_t xPtLimits[35] =  {	0, 	 0.3, 0.4, 0.5, 0.6, 
								0.8, 1.0, 1.2, 1.4, 1.6, 
								1.8, 2.0, 2.2, 2.4, 2.6,
								2.8, 3.0, 3.2, 3.4, 3.6,
								3.8, 4.0, 4.5, 5.0, 5.5,
								6.0, 7.0, 8.0, 10., 12.,
								14., 16., 18., 20., 25.	};

	Double_t xPtLimitsExt[39] =  {	0, 	 0.3, 0.4, 0.5, 0.6, 
								0.8, 1.0, 1.2, 1.4, 1.6, 
								1.8, 2.0, 2.2, 2.4, 2.6,
								2.8, 3.0, 3.2, 3.4, 3.6,
								3.8, 4.0, 4.5, 5.0, 5.5,
								6.0, 7.0, 8.0, 9.0, 10.,
								11., 12., 13., 14., 15., 
								16., 18., 20., 25.	};
								
	// Definition of offsets for stat & sys see output of function in shell, make sure pt bins match							
	Int_t offSets[11]	= 	{	0, 5, 4, 0, 0,
								0, 0, 0, 0, 0, 0};
	Int_t offSetsSys[11]= 	{	1, 5, 4, 0, 2,
								0, 0, 0, 0, 0, 0};

	//	**********************************************************************************************************************
	//	******************************************* Recalculating published spectrum *****************************************
	//	**********************************************************************************************************************
	TGraph* graphWeightsOld[11];
	for (Int_t i = 0; i< 11; i++){
		graphWeightsOld[i] = NULL;
	}	
						
	// Declaration & calculation of combined spectrum							
	TString fileNameOutputWeightingOld							= Form("%s/WeightingOld.dat",outputDir.Data());
	TGraphAsymmErrors* graphCombPi0InvCrossSectionStat7TeVOld= NULL;
	TGraphAsymmErrors* graphCombPi0InvCrossSectionSys7TeVOld = NULL;
	TGraphAsymmErrors* graphCombPi0InvCrossSectionTot7TeVOld = CombinePtPointsSpectraFullCorrMat( 	statErrorCollection,	sysErrorCollection, 	
																										xPtLimits, 34,
																										offSets, offSetsSys,
																										graphCombPi0InvCrossSectionStat7TeVOld, graphCombPi0InvCrossSectionSys7TeVOld,
																										fileNameOutputWeightingOld,1
																									);
// 	graphCombPi0InvCrossSectionStat7TeVOld->RemovePoint(0);
// 	graphCombPi0InvCrossSectionStat7TeVOld->RemovePoint(graphCombPi0InvCrossSectionStat7TeVOld->GetN()-1);
	graphCombPi0InvCrossSectionStat7TeVOld->Print();
	
	// Reading weights from output file for plotting
	ifstream fileWeightsReadOld;
	fileWeightsReadOld.open(fileNameOutputWeightingOld,ios_base::in);
	cout << "reading" << fileNameOutputWeightingOld << endl;
	Double_t xValuesReadOld[50];
	Double_t weightsReadOld[11][50];
	Int_t availableMeasOld[11] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
	Int_t nMeasSetOld 			= 2;
	Int_t nPtBinsReadOld 		= 0;
	while(!fileWeightsReadOld.eof() && nPtBinsReadOld < 50){
		TString garbage = "";
		if (nPtBinsReadOld == 0){
			fileWeightsReadOld >> garbage ;//>> availableMeas[0] >> availableMeas[1] >> availableMeas[2] >> availableMeas[3];
			for (Int_t i = 0; i < nMeasSetOld; i++){
				fileWeightsReadOld >> availableMeasOld[i] ;
			}	
			cout << "read following measurements: "; 
			for (Int_t i = 0; i < 11; i++){
				cout << availableMeasOld[i] << "\t" ;
			}	
			cout << endl;
		} else {
			fileWeightsReadOld >> xValuesReadOld[nPtBinsReadOld-1];
			for (Int_t i = 0; i < nMeasSetOld; i++){
				fileWeightsReadOld >> weightsReadOld[availableMeasOld[i]][nPtBinsReadOld-1] ;
			}	
			cout << "read: "<<  nPtBinsReadOld << "\t"<< xValuesReadOld[nPtBinsReadOld-1] << "\t" ;
			for (Int_t i = 0; i < nMeasSetOld; i++){
				cout << weightsReadOld[availableMeasOld[i]][nPtBinsReadOld-1] << "\t";
			}
			cout << endl;
		}
		nPtBinsReadOld++;
	}
	nPtBinsReadOld = nPtBinsReadOld-2 ;
	fileWeightsReadOld.close();

	for (Int_t i = 0; i < nMeasSetOld; i++){
		graphWeightsOld[availableMeasOld[i]] = new TGraph(nPtBinsReadOld,xValuesReadOld,weightsReadOld[availableMeasOld[i]]);
		Int_t bin = 0;
		for (Int_t n = 0; n< nPtBinsReadOld; n++){
			if (graphWeightsOld[availableMeasOld[i]]->GetY()[bin] == 0) graphWeightsOld[availableMeasOld[i]]->RemovePoint(bin);
			else bin++;
		}	
	}	

	//	**********************************************************************************************************************
	//	******************************************* Plotting weights method only EMC *****************************************
	//	**********************************************************************************************************************
	Int_t textSizeLabelsPixel = 900*0.04;

	TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasWeights, 0.08, 0.02, 0.035, 0.09);
    canvasWeights->SetLogx();
   
	TH2F * histo2DWeights;
	histo2DWeights = new TH2F("histo2DWeights","histo2DWeights",11000,0.23,70.,1000,-0.5,1.1);
	SetStyleHistoTH2ForGraphs(histo2DWeights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
	histo2DWeights->GetXaxis()->SetMoreLogLabels();
	histo2DWeights->GetXaxis()->SetLabelOffset(-0.01);
// 	histo2DWeights->GetYaxis()->SetRangeUser(-10,10);
	histo2DWeights->Draw("copy");
	
		TLegend* legendAccWeightsOld = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetOld*1.35), 32);
		for (Int_t i = 0; i < nMeasSetOld; i++){
			DrawGammaSetMarkerTGraph(graphWeightsOld[availableMeasOld[i]],
									 markerStyleDet[availableMeasOld[i]], 
									 markerSizeDet[availableMeasOld[i]]*0.5, 
									 colorDet[availableMeasOld[i]] , 
									 colorDet[availableMeasOld[i]]);
			graphWeightsOld[availableMeasOld[i]]->Draw("p,same,e1");
			legendAccWeightsOld->AddEntry(graphWeightsOld[availableMeasOld[i]],nameMeasGlobal[availableMeasOld[i]],"p");
		}	
		legendAccWeightsOld->Draw();

		TLatex *labelWeightsEnergy = new TLatex(0.7,0.20,collisionSystem7TeV.Data());
		SetStyleTLatex( labelWeightsEnergy, 0.85*textSizeLabelsPixel,4);
		labelWeightsEnergy->SetTextFont(43);
		labelWeightsEnergy->Draw();
		TLatex *labelWeightsPi0 = new TLatex(0.7,0.16,"#pi^{0} #rightarrow #gamma#gamma");
		SetStyleTLatex( labelWeightsPi0, 0.85*textSizeLabelsPixel,4);
		labelWeightsPi0->SetTextFont(43);
		labelWeightsPi0->Draw();

// 		DrawGammaLines(0.23, 70. , 0.8, 0.8,0.1, kGray, 3);
		DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
		DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);
		
	canvasWeights->SaveAs(Form("%s/WeightsOldPublished.%s",outputDir.Data(),suffix.Data()));

	//	**********************************************************************************************************************
	//	******************************************* Calculation of spectrum including EMCal only *****************************
	//	**********************************************************************************************************************
	statErrorCollection[2] = (TH1D*)histoEMCALPi0InvCrossSection7TeV->Clone("statErrEMCALPi0");
	sysErrorCollection[2] = (TGraphAsymmErrors*)graphEMCALPi0InvCrossSectionSys7TeV->Clone("sysErrEMCALPi0");
			
	TH1D* statErrorRelCollection[11];
	for (Int_t i = 0; i< 11; i++){
		statErrorRelCollection[i] = NULL;
	}	
	for (Int_t i = 0; i < 11; i++){
		if (statErrorCollection[i]) statErrorRelCollection[i] = CalculateRelErrUpTH1D( statErrorCollection[i], Form("relativeStatError_%s", nameMeasGlobal[i].Data()));
		if (i == 2){
			for (Int_t b = 1; b < statErrorCollection[i]->GetNbinsX()+1; b++) {
				cout <<  statErrorCollection[i]->GetBinError(b)/statErrorCollection[i]->GetBinContent(b)*110 << "\t" << statErrorRelCollection[i]->GetBinContent(b) << endl;
			}	
		}	
	}
	
	TGraphAsymmErrors* sysErrorRelCollection[11];
	for (Int_t i = 0; i< 11; i++){
		sysErrorRelCollection[i] = NULL;
	}	
	for (Int_t i = 0; i < 11; i++){
		if (sysErrorCollection[i]) sysErrorRelCollection[i] = CalculateRelErrUpAsymmGraph( sysErrorCollection[i], Form("relativeSysError_%s", nameMeasGlobal[i].Data()));
	}

	
	TGraph* graphWeightsOEMC[11];
	for (Int_t i = 0; i< 11; i++){
		graphWeightsOEMC[i] = NULL;
	}	

	// Declaration & calculation of combined spectrum							
	TString fileNameOutputWeightingOEMC							= Form("%s/WeightingOEMC.dat",outputDir.Data());
	TGraphAsymmErrors* graphCombPi0InvCrossSectionStat7TeVOEMC= NULL;
	TGraphAsymmErrors* graphCombPi0InvCrossSectionSys7TeVOEMC = NULL;
	TGraphAsymmErrors* graphCombPi0InvCrossSectionTot7TeVOEMC = CombinePtPointsSpectraFullCorrMat( 	statErrorCollection,	sysErrorCollection, 	
																										xPtLimitsExt, 38,
																										offSets, offSetsSys,
																										graphCombPi0InvCrossSectionStat7TeVOEMC, graphCombPi0InvCrossSectionSys7TeVOEMC,
																										fileNameOutputWeightingOEMC,1
																									);
	
// 	graphCombPi0InvCrossSectionStat7TeVOEMC->RemovePoint(0);
// 	graphCombPi0InvCrossSectionSys7TeVOEMC->RemovePoint(0);
	graphCombPi0InvCrossSectionStat7TeVOEMC->Print();
	
	// Reading weights from output file for plotting
	ifstream fileWeightsReadOEMC;
	fileWeightsReadOEMC.open(fileNameOutputWeightingOEMC,ios_base::in);
	cout << "reading" << fileNameOutputWeightingOEMC << endl;
	Double_t xValuesReadOEMC[50];
	Double_t weightsReadOEMC[11][50];
	Int_t availableMeasOEMC[11] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
	Int_t nMeasSetOEMC 			= 3;
	Int_t nPtBinsReadOEMC 		= 0;
	while(!fileWeightsReadOEMC.eof() && nPtBinsReadOEMC < 50){
		TString garbage = "";
		if (nPtBinsReadOEMC == 0){
			fileWeightsReadOEMC >> garbage ;//>> availableMeas[0] >> availableMeas[1] >> availableMeas[2] >> availableMeas[3];
			for (Int_t i = 0; i < nMeasSetOEMC; i++){
				fileWeightsReadOEMC >> availableMeasOEMC[i] ;
			}	
			cout << "read following measurements: "; 
			for (Int_t i = 0; i < 11; i++){
				cout << availableMeasOEMC[i] << "\t" ;
			}	
			cout << endl;
		} else {
			fileWeightsReadOEMC >> xValuesReadOEMC[nPtBinsReadOEMC-1];
			for (Int_t i = 0; i < nMeasSetOEMC; i++){
				fileWeightsReadOEMC >> weightsReadOEMC[availableMeasOEMC[i]][nPtBinsReadOEMC-1] ;
			}	
			cout << "read: "<<  nPtBinsReadOEMC << "\t"<< xValuesReadOEMC[nPtBinsReadOEMC-1] << "\t" ;
			for (Int_t i = 0; i < nMeasSetOEMC; i++){
				cout << weightsReadOEMC[availableMeasOEMC[i]][nPtBinsReadOEMC-1] << "\t";
			}
			cout << endl;
		}
		nPtBinsReadOEMC++;
	}
	nPtBinsReadOEMC = nPtBinsReadOEMC-2 ;
	fileWeightsReadOEMC.close();

	for (Int_t i = 0; i < nMeasSetOEMC; i++){
		graphWeightsOEMC[availableMeasOEMC[i]] = new TGraph(nPtBinsReadOEMC,xValuesReadOEMC,weightsReadOEMC[availableMeasOEMC[i]]);
		Int_t bin = 0;
		for (Int_t n = 0; n< nPtBinsReadOEMC; n++){
			if (graphWeightsOEMC[availableMeasOEMC[i]]->GetY()[bin] == 0) graphWeightsOEMC[availableMeasOEMC[i]]->RemovePoint(bin);
			else bin++;
		}	
	}	

	//	**********************************************************************************************************************
	//	******************************************* Plotting weights method only EMC *****************************************
	//	**********************************************************************************************************************
    canvasWeights->cd();
	histo2DWeights->Draw("copy");
	
		TLegend* legendAccWeightsOEMC = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetOEMC*1.35), 32);
		for (Int_t i = 0; i < nMeasSetOEMC; i++){
			DrawGammaSetMarkerTGraph(graphWeightsOEMC[availableMeasOEMC[i]],
									 markerStyleDet[availableMeasOEMC[i]], 
									 markerSizeDet[availableMeasOEMC[i]]*0.5, 
									 colorDet[availableMeasOEMC[i]] , 
									 colorDet[availableMeasOEMC[i]]);
			graphWeightsOEMC[availableMeasOEMC[i]]->Draw("p,same,e1");
			legendAccWeightsOEMC->AddEntry(graphWeightsOEMC[availableMeasOEMC[i]],nameMeasGlobal[availableMeasOEMC[i]],"p");
		}	
		legendAccWeightsOEMC->Draw();

		labelWeightsEnergy->Draw();
		labelWeightsPi0->Draw();

// 		DrawGammaLines(0.23, 70. , 0.8, 0.8,0.1, kGray, 3);
		DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
		DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);
		
	canvasWeights->SaveAs(Form("%s/WeightsOnlyEMCal.%s",outputDir.Data(),suffix.Data()));

	

	//	**********************************************************************************************************************
	//	******************************************* Compare new spectrum to old average **************************************
	//	**********************************************************************************************************************
	
	TGraphAsymmErrors* graphCombPi0InvCrossSectionStat7TeVOEMCForCompToOld= (TGraphAsymmErrors*)graphCombPi0InvCrossSectionStat7TeVOEMC->Clone("graphCombPi0InvCrossSectionStat7TeVOEMCForCompToOld");
	TGraphAsymmErrors* graphCombPi0InvCrossSectionSys7TeVOEMCForCompToOld = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionSys7TeVOEMC->Clone("graphCombPi0InvCrossSectionSys7TeVOEMCForCompToOld");
	
	for (Int_t n = graphCombPi0InvCrossSectionStat7TeVOEMCForCompToOld->GetN()-1; graphCombPi0InvCrossSectionStat7TeVOEMCForCompToOld->GetX()[n] > 8; n-- ){
		graphCombPi0InvCrossSectionStat7TeVOEMCForCompToOld->RemovePoint(n);
	}	
	for (Int_t n = graphCombPi0InvCrossSectionSys7TeVOEMCForCompToOld->GetN()-1; graphCombPi0InvCrossSectionSys7TeVOEMCForCompToOld->GetX()[n] > 8; n-- ){
		graphCombPi0InvCrossSectionSys7TeVOEMCForCompToOld->RemovePoint(n);
	}	
	
	graphCombPi0InvCrossSectionStat7TeVOEMC->Print();
	graphCombPi0InvCrossSectionStat7TeVOld->Print();
	
// 	return;
	TGraphAsymmErrors*  graphRatioCombNewOEMCDivCombOldStat = CalculateGraphErrRatioToOtherTGraphErr(graphCombPi0InvCrossSectionStat7TeVOEMCForCompToOld,
																									 graphCombPi0InvCrossSectionStat7TeVOld, kTRUE);
	TGraphAsymmErrors*  graphRatioCombNewOECMDivCombOldSys 	= CalculateGraphErrRatioToOtherTGraphErr(graphCombPi0InvCrossSectionSys7TeVOEMCForCompToOld, 
																									 graphCombPi0InvCrossSectionSys7TeVOld, kTRUE);
	

	//	**********************************************************************************************************************
	//	******************************************* Ratio of Comb to Fit ****************************************
	//	**********************************************************************************************************************
	
	TCanvas* canvasRatioToOldCombined = new TCanvas("canvasRatioToOldCombined","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioToOldCombined, 0.1, 0.02, 0.035, 0.09);
    canvasRatioToOldCombined->SetLogx();
   
	TH2F * histo2DRatioToOldCombined;
	histo2DRatioToOldCombined = new TH2F("histo2DRatioToOldCombined","histo2DRatioToOldCombined",11000,0.23,10.,1000,0.75,1.25);
	SetStyleHistoTH2ForGraphs(histo2DRatioToOldCombined, "#it{p}_{T} (GeV/#it{c})","#frac{Comb}{Comb Old}",0.035,0.04, 0.035,0.04, 1.,1.,510,505);
	histo2DRatioToOldCombined->GetXaxis()->SetMoreLogLabels();
	histo2DRatioToOldCombined->GetXaxis()->SetLabelOffset(-0.01);
// 	histo2DRatioToOldCombined->GetYaxis()->SetRangeUser(-10,10);
	histo2DRatioToOldCombined->Draw("copy");

		DrawGammaSetMarkerTGraphAsym(graphRatioCombNewOECMDivCombOldSys, markerStyleComb+5, markerSizeComb, kGreen+2 , kGreen+2, widthLinesBoxes, kTRUE);
		graphRatioCombNewOECMDivCombOldSys->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphRatioCombNewOEMCDivCombOldStat, markerStyleComb+5, markerSizeComb, kGreen+2 , kGreen+2);
		graphRatioCombNewOEMCDivCombOldStat->Draw("p,same,e1");
		
		DrawGammaLines(0.23, 10. , 1., 1.,0.1, kGray+2);
		DrawGammaLines(0.23, 10. , 1.1, 1.1,0.1, kGray, 7);
		DrawGammaLines(0.23, 10. , 0.9, 0.9,0.1, kGray, 7);
		DrawGammaLines(0.23, 10. , 1.05, 1.05,0.1, kGray, 3);
		DrawGammaLines(0.23, 10. , 0.95, 0.95,0.1, kGray, 3);
		
		TLegend* legendRatioToOld = GetAndSetLegend2(0.67, 0.96-(0.035*3*1.35), 0.93, 0.96, 32);
		legendRatioToOld->AddEntry(graphRatioCombNewOECMDivCombOldSys,"PCM, PHOS, EMCal");
		legendRatioToOld->Draw();

		TLatex *labelRatioToOldEnergy = new TLatex(0.15,0.89,collisionSystem7TeV.Data());
		SetStyleTLatex( labelRatioToOldEnergy, 0.85*textSizeLabelsPixel,4);
		labelRatioToOldEnergy->SetTextFont(43);
		labelRatioToOldEnergy->Draw();
		TLatex *labelRatioToOldPi0 = new TLatex(0.15,0.85,"#pi^{0} #rightarrow #gamma#gamma");
		SetStyleTLatex( labelRatioToOldPi0, 0.85*textSizeLabelsPixel,4);
		labelRatioToOldPi0->SetTextFont(43);
		labelRatioToOldPi0->Draw();
		
	canvasRatioToOldCombined->SaveAs(Form("%s/RatioOfCombToCombOld_PP7TeV.%s",outputDir.Data(),suffix.Data()));
 
	
	// 	*********************************************************************************************************************
	// 	************************************ Visualize relative errors ******************************************************
	// 	*********************************************************************************************************************
	
	TCanvas* canvasRelSysErr = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelSysErr->SetLogx();
   
	TH2F * histo2DRelSysErr;
	histo2DRelSysErr = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,0.23,70.,1000,0,80.5);
	SetStyleHistoTH2ForGraphs(histo2DRelSysErr, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
	histo2DRelSysErr->GetXaxis()->SetMoreLogLabels();
	histo2DRelSysErr->GetXaxis()->SetLabelOffset(-0.01);
	histo2DRelSysErr->GetYaxis()->SetRangeUser(0,45.5);
	histo2DRelSysErr->Draw("copy");	

		TLegend* legendRelSysErr = GetAndSetLegend2(0.62, 0.94-(0.035*nMeasSetOEMC*1.35), 0.95, 0.94, 32);
		for (Int_t i = 0; i < nMeasSetOEMC; i++){
			DrawGammaSetMarkerTGraph(sysErrorRelCollection[availableMeasOEMC[i]], markerStyleDet[availableMeasOEMC[i]], markerSizeDet[availableMeasOEMC[i]]*0.5, colorDet[availableMeasOEMC[i]],
									 colorDet[availableMeasOEMC[i]]);
			sysErrorRelCollection[availableMeasOEMC[i]]->Draw("p,same,e1");
			legendRelSysErr->AddEntry(sysErrorRelCollection[availableMeasOEMC[i]],nameMeasGlobal[availableMeasOEMC[i]],"p");
		}	
		legendRelSysErr->Draw();

		TLatex *labelRelSysErrEnergy = new TLatex(0.15,0.89,collisionSystem7TeV.Data());
		SetStyleTLatex( labelRelSysErrEnergy, 0.85*textSizeLabelsPixel,4);
		labelRelSysErrEnergy->SetTextFont(43);
		labelRelSysErrEnergy->Draw();
		TLatex *labelRelSysErrPi0 = new TLatex(0.15,0.85,"#pi^{0} #rightarrow #gamma#gamma");
		SetStyleTLatex( labelRelSysErrPi0, 0.85*textSizeLabelsPixel,4);
		labelRelSysErrPi0->SetTextFont(43);
		labelRelSysErrPi0->Draw();
		
	canvasRelSysErr->SaveAs(Form("%s/RelSysErr.%s",outputDir.Data(),suffix.Data()));	
	
	// 	*********************************************************************************************************************
	// 	************************************ Visualize relative errors ******************************************************
	// 	*********************************************************************************************************************
	
	TCanvas* canvasRelStatErr = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelStatErr->SetLogx();
   
	TH2F * histo2DRelStatErr;
	histo2DRelStatErr = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,0.23,70.,1000,0,80.5);
	SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
	histo2DRelStatErr->GetXaxis()->SetMoreLogLabels();
	histo2DRelStatErr->GetXaxis()->SetLabelOffset(-0.01);
// 	histo2DRelStatErr->GetYaxis()->SetRangeUser(-10,10);
	histo2DRelStatErr->GetYaxis()->SetRangeUser(0,45.5);	
	histo2DRelStatErr->Draw("copy");	
	
		TLegend* legendRelStatErr = GetAndSetLegend2(0.14, 0.94-(0.035*nMeasSetOEMC*1.35), 0.45, 0.94, 32);
		for (Int_t i = 0; i < nMeasSetOEMC; i++){
			if (availableMeasOEMC[i]== 2){
				DrawGammaSetMarker(statErrorRelCollection[availableMeasOEMC[i]], markerStyleDet[availableMeasOEMC[i]], markerSizeDet[availableMeasOEMC[i]]*0.5, colorDet[availableMeasOEMC[i]] , 
							   colorDet[availableMeasOEMC[i]]);
				TGraphAsymmErrors* graphDummy 	= new TGraphAsymmErrors(statErrorRelCollection[availableMeasOEMC[i]]);
				DrawGammaSetMarkerTGraphAsym(graphDummy, markerStyleDet[availableMeasOEMC[i]], markerSizeDet[availableMeasOEMC[i]]*0.5, colorDet[availableMeasOEMC[i]],
									 colorDet[availableMeasOEMC[i]]);
				graphDummy->Draw("same,p,x0");
				legendRelStatErr->AddEntry(graphDummy,nameMeasGlobal[availableMeasOEMC[i]],"p");
			
			} else {
				DrawGammaSetMarker(statErrorRelCollection[availableMeasOEMC[i]], markerStyleDet[availableMeasOEMC[i]], markerSizeDet[availableMeasOEMC[i]]*0.5, colorDet[availableMeasOEMC[i]] , 
							   colorDet[availableMeasOEMC[i]]);
				statErrorRelCollection[availableMeasOEMC[i]]->Draw("p,same,e1");
				legendRelStatErr->AddEntry(statErrorRelCollection[availableMeasOEMC[i]],nameMeasGlobal[availableMeasOEMC[i]],"p");
	
			}	
		}	
		legendRelStatErr->Draw();

		TLatex *labelRelStatErrEnergy = new TLatex(0.75,0.89,collisionSystem7TeV.Data());
		SetStyleTLatex( labelRelStatErrEnergy, 0.85*textSizeLabelsPixel,4);
		labelRelStatErrEnergy->SetTextFont(43);
		labelRelStatErrEnergy->Draw();
		TLatex *labelRelStatErrPi0 = new TLatex(0.75,0.85,"#pi^{0} #rightarrow #gamma#gamma");
		SetStyleTLatex( labelRelStatErrPi0, 0.85*textSizeLabelsPixel,4);
		labelRelStatErrPi0->SetTextFont(43);
		labelRelStatErrPi0->Draw();
		
	canvasRelStatErr->SaveAs(Form("%s/RelStatErr.%s",outputDir.Data(),suffix.Data()));
		

	//************************************************************************************************************************
	//************************************** Comparison sys and stat for new and old combined ********************************
	//************************************************************************************************************************
	TGraphAsymmErrors* graphCombPi0InvCrossSectionRelStat7TeVOEMC 	= CalculateRelErrUpAsymmGraph( graphCombPi0InvCrossSectionStat7TeVOEMC, "relativeStatError_OEMC");
	TGraphAsymmErrors* graphCombPi0InvCrossSectionRelSys7TeVOEMC 	= CalculateRelErrUpAsymmGraph( graphCombPi0InvCrossSectionSys7TeVOEMC, "relativeSysError_OEMC");
	TGraphAsymmErrors* graphCombPi0InvCrossSectionRelStat7TeVOld 	= CalculateRelErrUpAsymmGraph( graphCombPi0InvCrossSectionStat7TeVOld, "relativeStatError_Old");
	TGraphAsymmErrors* graphCombPi0InvCrossSectionRelSys7TeVOld 	= CalculateRelErrUpAsymmGraph( graphCombPi0InvCrossSectionSys7TeVOld, "relativeSysError_Old");

	TGraphAsymmErrors* graphCombPi0InvCrossSectionRelTot7TeVOEMC 	= CalculateRelErrUpAsymmGraph( graphCombPi0InvCrossSectionTot7TeVOEMC, "relativeTotalError_OEMC");
	TGraphAsymmErrors* graphCombPi0InvCrossSectionRelTot7TeVOld 	= CalculateRelErrUpAsymmGraph( graphCombPi0InvCrossSectionTot7TeVOld, "relativeTotalError_Old");

	TCanvas* canvasRelTotErr = new TCanvas("canvasRelTotErr","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRelTotErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelTotErr->SetLogx();
   
	TH2F * histo2DRelTotErr;
	histo2DRelTotErr = new TH2F("histo2DRelTotErr","histo2DRelTotErr",11000,0.23,70.,1000,0,57.5);
	SetStyleHistoTH2ForGraphs(histo2DRelTotErr, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
	histo2DRelTotErr->GetXaxis()->SetMoreLogLabels();
	histo2DRelTotErr->GetXaxis()->SetLabelOffset(-0.01);
// 	histo2DRelTotErr->GetYaxis()->SetRangeUser(-10,10);
	histo2DRelTotErr->Draw("copy");	

		DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelTot7TeVOld, markerStyleComb+1, markerSizeComb, kBlack , kBlack);
		graphCombPi0InvCrossSectionRelTot7TeVOld->Draw("p,same,e1");
		DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelTot7TeVOEMC, markerStyleComb+5, markerSizeComb, kGreen+2 , kGreen+2);
		graphCombPi0InvCrossSectionRelTot7TeVOEMC->Draw("p,same,e1");

		TLegend* legendRelTotErr = GetAndSetLegend2(0.14, 0.94-(0.035*4*1.35), 0.45, 0.94, 32);
		legendRelTotErr->AddEntry(graphCombPi0InvCrossSectionRelTot7TeVOld,"PCM, PHOS","p");
		legendRelTotErr->AddEntry(graphCombPi0InvCrossSectionRelTot7TeVOEMC,"PCM, PHOS, EMCAL","p");
		legendRelTotErr->Draw();

		TLatex *labelRelTotErrEnergy = new TLatex(0.75,0.89,collisionSystem7TeV.Data());
		SetStyleTLatex( labelRelTotErrEnergy, 0.85*textSizeLabelsPixel,4);
		labelRelTotErrEnergy->SetTextFont(43);
		labelRelTotErrEnergy->Draw();
		TLatex *labelRelTotErrPi0 = new TLatex(0.75,0.85,"#pi^{0} #rightarrow #gamma#gamma");
		SetStyleTLatex( labelRelTotErrPi0, 0.85*textSizeLabelsPixel,4);
		labelRelTotErrPi0->SetTextFont(43);
		labelRelTotErrPi0->Draw();
		
	canvasRelTotErr->SaveAs(Form("%s/RelTotErr.%s",outputDir.Data(),suffix.Data()));

	histo2DRelSysErr->Draw("copy");	

		DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelSys7TeVOld, markerStyleComb+1, markerSizeComb, kBlack , kBlack);
		graphCombPi0InvCrossSectionRelSys7TeVOld->Draw("p,same,e1");
		DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelSys7TeVOEMC, markerStyleComb+5, markerSizeComb, kGreen+2 , kGreen+2);
		graphCombPi0InvCrossSectionRelSys7TeVOEMC->Draw("p,same,e1");

		legendRelTotErr->Draw();

		labelRelTotErrEnergy->Draw();
		labelRelTotErrPi0->Draw();
		
	canvasRelTotErr->SaveAs(Form("%s/RelSysErr_diffComb.%s",outputDir.Data(),suffix.Data()));

	histo2DRelStatErr->Draw("copy");	

		DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelStat7TeVOld, markerStyleComb+1, markerSizeComb, kBlack , kBlack);
		graphCombPi0InvCrossSectionRelStat7TeVOld->Draw("p,same,e1");
		DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelStat7TeVOEMC, markerStyleComb+5, markerSizeComb, kGreen+2 , kGreen+2);
		graphCombPi0InvCrossSectionRelStat7TeVOEMC->Draw("p,same,e1");

		legendRelTotErr->Draw();

		labelRelTotErrEnergy->Draw();
		labelRelTotErrPi0->Draw();
		
	canvasRelTotErr->SaveAs(Form("%s/RelStatErr_diffComb.%s",outputDir.Data(),suffix.Data()));
	
	histo2DRelTotErr->GetYaxis()->SetRangeUser(0,57.5);
	histo2DRelTotErr->GetYaxis()->SetTitle("Err (%)");
		histo2DRelTotErr->Draw("copy");	

		graphCombPi0InvCrossSectionRelTot7TeVOEMC->Draw("p,same,e1");
		DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelStat7TeVOEMC, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
		graphCombPi0InvCrossSectionRelStat7TeVOEMC->Draw("l,x0,same,e1");
		DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelSys7TeVOEMC, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
		graphCombPi0InvCrossSectionRelSys7TeVOEMC->SetLineStyle(7);
		graphCombPi0InvCrossSectionRelSys7TeVOEMC->Draw("l,x0,same,e1");

		TLegend* legendRelTotErr2 = GetAndSetLegend2(0.14, 0.94-(0.035*3*1.35), 0.45, 0.94, 32);
		legendRelTotErr2->AddEntry(graphCombPi0InvCrossSectionRelTot7TeVOEMC,"tot","p");
		legendRelTotErr2->AddEntry(graphCombPi0InvCrossSectionRelStat7TeVOEMC,"stat","l");
		legendRelTotErr2->AddEntry(graphCombPi0InvCrossSectionRelSys7TeVOEMC,"sys","l");
		legendRelTotErr2->Draw();

		labelRelTotErrEnergy->Draw();
		labelRelTotErrPi0->Draw();
		
	canvasRelTotErr->SaveAs(Form("%s/RelOEMCdecomp.%s",outputDir.Data(),suffix.Data()));

	histo2DRelTotErr->GetYaxis()->SetRangeUser(0,57.5);
	histo2DRelTotErr->GetYaxis()->SetTitle("Err (%)");
		histo2DRelTotErr->Draw("copy");	

		DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelTot7TeVOld, markerStyleComb, markerSizeComb, colorComb , colorComb);		
		graphCombPi0InvCrossSectionRelTot7TeVOld->Draw("p,same,e1");
		DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelStat7TeVOld, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
		graphCombPi0InvCrossSectionRelStat7TeVOld->Draw("l,x0,same,e1");
		DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelSys7TeVOld, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
		graphCombPi0InvCrossSectionRelSys7TeVOld->SetLineStyle(7);
		graphCombPi0InvCrossSectionRelSys7TeVOld->Draw("l,x0,same,e1");

		TLegend* legendRelTotErr3 = GetAndSetLegend2(0.14, 0.94-(0.035*3*1.35), 0.45, 0.94, 32);
		legendRelTotErr3->AddEntry(graphCombPi0InvCrossSectionRelTot7TeVOld,"tot","p");
		legendRelTotErr3->AddEntry(graphCombPi0InvCrossSectionRelStat7TeVOld,"stat","l");
		legendRelTotErr3->AddEntry(graphCombPi0InvCrossSectionRelSys7TeVOld,"sys","l");
		legendRelTotErr3->Draw();

		labelRelTotErrEnergy->Draw();
		labelRelTotErrPi0->Draw();
		
	canvasRelTotErr->SaveAs(Form("%s/RelOlddecomp.%s",outputDir.Data(),suffix.Data()));

	
	//	**********************************************************************************************************************
	//	************************************* Calculating bin shifted spectra & fitting **************************************
	//	**********************************************************************************************************************
	
	// Cloning spectra
	TGraphAsymmErrors* graphCombPi0InvCrossSectionTot7TeVOEMCUnShifted 	= (TGraphAsymmErrors*)graphCombPi0InvCrossSectionTot7TeVOEMC->Clone("Unshifted"); 
	TGraphAsymmErrors* graphCombPi0InvCrossSectionStat7TeVOEMCUnShifted = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionStat7TeVOEMC->Clone("UnshiftedStat"); 
	TGraphAsymmErrors* graphCombPi0InvCrossSectionSys7TeVOEMCUnShifted 	= (TGraphAsymmErrors*)graphCombPi0InvCrossSectionSys7TeVOEMC->Clone("UnshiftedSys"); 

	TGraphAsymmErrors* graphPCMPi0InvCrossSectionStat7TeVUnShifted 		= (TGraphAsymmErrors*)graphPCMPi0InvCrossSectionStat7TeV->Clone("UnshiftedStatPCM"); 
	TGraphAsymmErrors* graphPCMPi0InvCrossSectionSys7TeVUnShifted 		= (TGraphAsymmErrors*)graphPCMPi0InvCrossSectionSys7TeV->Clone("UnshiftedSysPCM"); 

	TGraphAsymmErrors* graphPHOSPi0InvCrossSectionStat7TeVUnshifted 	= (TGraphAsymmErrors*)graphPHOSPi0InvCrossSectionStat7TeV->Clone("UnshiftedStatPHOS"); 
	TGraphAsymmErrors* graphPHOSPi0InvCrossSectionSys7TeVUnshifted 		= (TGraphAsymmErrors*)graphPHOSPi0InvCrossSectionSys7TeV->Clone("UnshiftedSysPHOS"); 
	
	TGraphAsymmErrors* graphEMCALPi0InvCrossSectionStat7TeVUnshifted 	= (TGraphAsymmErrors*)graphEMCALPi0InvCrossSectionStat7TeV->Clone("UnshiftedStatEMCAL"); 
	TGraphAsymmErrors* graphEMCALPi0InvCrossSectionSys7TeVUnshifted 	= (TGraphAsymmErrors*)graphEMCALPi0InvCrossSectionSys7TeV->Clone("UnshiftedSysEMCAL"); 
	
	// Calculating binshifts
	Double_t paramGraph[3] 								= {1.0e12, 8., 0.13};	
	TF1* fitInvCrossSectionPi07TeV 						= FitObject("l","fitInvCrossSectionPi07TeV","Pi0",histoEMCALPi0InvCrossSection7TeV,0.3,25.,paramGraph,"QNRMEX0+");
	TF1* fitInvCrossSectionPi07TeVGraph 				= (TF1*)fitInvCrossSectionPi07TeV->Clone("fitInvCrossSectionPi07TeVGraph"); 
	
	if(bWCorrection.CompareTo("X")==0 ){
		TF1* fitTsallisPi07TeVPtMult 				= FitObject("tmpt","TsallisMultWithPtPi07TeV","Pi0");
		fitTsallisPi07TeVPtMult->SetParameters(paramGraph[0],paramGraph[1], paramGraph[2]) ; // standard parameter optimize if necessary
		graphCombPi0InvCrossSectionTot7TeVOEMC			= ApplyXshift(graphCombPi0InvCrossSectionTot7TeVOEMC, fitTsallisPi07TeVPtMult);
		graphCombPi0InvCrossSectionStat7TeVOEMC 		= ApplyXshiftIndividualSpectra (graphCombPi0InvCrossSectionTot7TeVOEMC, 
																						graphCombPi0InvCrossSectionStat7TeVOEMC, 
																						fitTsallisPi07TeVPtMult,
																						0, graphCombPi0InvCrossSectionStat7TeVOEMC->GetN());
		graphCombPi0InvCrossSectionSys7TeVOEMC 			= ApplyXshiftIndividualSpectra (graphCombPi0InvCrossSectionTot7TeVOEMC, 
																						graphCombPi0InvCrossSectionSys7TeVOEMC, 
																						fitTsallisPi07TeVPtMult, 
																						0, graphCombPi0InvCrossSectionSys7TeVOEMC->GetN());
		graphPCMPi0InvCrossSectionStat7TeV				= ApplyXshiftIndividualSpectra( graphCombPi0InvCrossSectionTot7TeVOEMC,
																					    graphPCMPi0InvCrossSectionStat7TeV,
																					    fitTsallisPi07TeVPtMult, 
																						0, 25);
		graphPCMPi0InvCrossSectionSys7TeV				= ApplyXshiftIndividualSpectra( graphCombPi0InvCrossSectionTot7TeVOEMC, 
																						graphPCMPi0InvCrossSectionSys7TeV, 
																						fitTsallisPi07TeVPtMult, 
																						0, 25);
		graphPHOSPi0InvCrossSectionStat7TeV 			= ApplyXshiftIndividualSpectra( graphCombPi0InvCrossSectionTot7TeVOEMC, 
																						graphPHOSPi0InvCrossSectionStat7TeV, 
																						fitTsallisPi07TeVPtMult,
																						4, 25);
		graphPHOSPi0InvCrossSectionSys7TeV 				= ApplyXshiftIndividualSpectra( graphCombPi0InvCrossSectionTot7TeVOEMC, 
																						graphPHOSPi0InvCrossSectionSys7TeV, 
																						fitTsallisPi07TeVPtMult,
																						4, 25);
		graphEMCALPi0InvCrossSectionStat7TeV 			= ApplyXshiftIndividualSpectra( graphCombPi0InvCrossSectionTot7TeVOEMC, 
																						graphEMCALPi0InvCrossSectionStat7TeV, 
																						fitTsallisPi07TeVPtMult,
																						3, graphCombPi0InvCrossSectionTot7TeVOEMC->GetN());		
		graphEMCALPi0InvCrossSectionSys7TeV 			= ApplyXshiftIndividualSpectra( graphCombPi0InvCrossSectionTot7TeVOEMC, 
																						graphEMCALPi0InvCrossSectionSys7TeV, 
																						fitTsallisPi07TeVPtMult,
																						3, graphCombPi0InvCrossSectionTot7TeVOEMC->GetN());		

				
		TCanvas* canvasDummy2 = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
		DrawGammaCanvasSettings( canvasDummy2,  0.13, 0.01, 0.015, 0.08);
		canvasDummy2->SetLogy();
		canvasDummy2->SetLogx();
		TH2F * histo2DDummy2;
		histo2DDummy2 = new TH2F("histo2DDummy2","histo2DDummy2",1000,0.23,70.,1000,1e-1,10e11);
		SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
		histo2DDummy2->DrawCopy(); 
	
		DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionStat7TeVOEMC, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
		graphCombPi0InvCrossSectionStat7TeVOEMC->Draw("pEsame");
		DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionTot7TeVOEMC, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
		graphCombPi0InvCrossSectionTot7TeVOEMC->Draw("pEsame");

		fitInvCrossSectionPi07TeV->SetLineColor(kBlue+2);
		fitInvCrossSectionPi07TeV->Draw("same");
		
		canvasDummy2->Update();
		canvasDummy2->Print(Form("%s/ComparisonShiftedPi0_7TeV.%s",outputDir.Data(),suffix.Data()));
	}
		
	graphCombPi0InvCrossSectionTot7TeVOEMC->Fit(fitInvCrossSectionPi07TeV,"QNRMEX0+","",0.4,50.);
	
	fitInvCrossSectionPi07TeV = FitObject("l","fitInvCrossSectionPi07TeV","Pi0",graphCombPi0InvCrossSectionTot7TeVOEMC,0.4,50.,paramGraph,"QNRMEX0+");
	fitInvCrossSectionPi07TeV = FitObject("l","fitInvCrossSectionPi07TeV","Pi0",graphCombPi0InvCrossSectionTot7TeVOEMC,0.4,50. ,paramGraph,"QNRMEX0+");

	cout << WriteParameterToFile(fitInvCrossSectionPi07TeV)<< endl;
	
	Double_t paramTCM[5] = {graphCombPi0InvCrossSectionTot7TeVOEMC->GetY()[0],0.3,graphCombPi0InvCrossSectionTot7TeVOEMC->GetY()[0]/10000,0.8,3};
	TF1* fitTCMInvCrossSectionPi07TeV = FitObject("tcm","fitTCMInvCrossSectionPi07TeV","Pi0",graphCombPi0InvCrossSectionTot7TeVOEMC,0.4,50.,paramTCM,"QNRMEX0+");
	fitTCMInvCrossSectionPi07TeV = FitObject("tcm","fitTCMInvCrossSectionPi07TeV","Pi0",graphCombPi0InvCrossSectionTot7TeVOEMC,0.4,50. ,paramTCM,"QNRMEX0+");
// 	TF1* fitTCMDecomposedL 				= new TF1("twoCompModel_DecLow",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1])",mesonMassExpectPi0,mesonMassExpectPi0,mesonMassExpectPi0));
// 	fitTCMDecomposedL->SetParameters(fitTCMInvCrossSectionPi07TeV->GetParameter(0),fitTCMInvCrossSectionPi07TeV->GetParameter(1));
// 	TF1 *fitTCMDecomposedH 				= new TF1("twoCompModel_DecH","[0]/(1 + x*x/TMath::Power(1+x*x/([1]*[1]*[2]),-[2]) )");
// // 		graphCombPi0InvCrossSectionTot7TeVOEMC->Fit(fitTCMDecomposedH,"QNRMEX0+","",5,20);
// 	fitTCMDecomposedH->SetParameters(fitTCMInvCrossSectionPi07TeV->GetParameter(2),fitTCMInvCrossSectionPi07TeV->GetParameter(3), fitTCMInvCrossSectionPi07TeV->GetParameter(4));
	cout << WriteParameterToFile(fitTCMInvCrossSectionPi07TeV)<< endl;
	
	TString forOutput= WriteParameterToFile(fitInvCrossSectionPi07TeV);
	cout<< forOutput.Data()<< endl;
	
// 	TH1D* histoRatioPythia8ToFit7TeV 					= (TH1D*) histoPythia8InvXSection->Clone();     
// 	histoRatioPythia8ToFit7TeV 							= CalculateHistoRatioToFit (histoRatioPythia8ToFit7TeV, fitTCMInvCrossSectionPi07TeV); 
// 	TH1D* histoRatioPythia8VarBinningToFit7TeV 			= (TH1D*) histoPythia8InvXSection_VarBinning->Clone();     
// 	histoRatioPythia8VarBinningToFit7TeV 				= CalculateHistoRatioToFit (histoRatioPythia8VarBinningToFit7TeV, fitTCMInvCrossSectionPi07TeV); 

// 	TGraph* graphRatioCombNLOPi07TeVMuHalf				= (TGraph*)graphNLOCalcMuHalf7TeV->Clone();
// 	TGraph* graphRatioCombNLOPi07TeVMuOne				= (TGraph*)graphNLOCalcMuOne7TeV->Clone();
// 	TGraph* graphRatioCombNLOPi07TeVMuTwo				= (TGraph*)graphNLOCalcMuTwo7TeV->Clone();	
// 	TGraphAsymmErrors* graphRatioCombNLODSS14Pi07TeV		= (TGraphAsymmErrors*)graphNLODSS14Calc7TeV->Clone();	
// 	graphRatioCombNLOPi07TeVMuHalf 						= CalculateGraphRatioToFit (graphRatioCombNLOPi07TeVMuHalf, fitTCMInvCrossSectionPi07TeV); 
// 	graphRatioCombNLOPi07TeVMuOne 						= CalculateGraphRatioToFit (graphRatioCombNLOPi07TeVMuOne, fitTCMInvCrossSectionPi07TeV); 
// 	graphRatioCombNLOPi07TeVMuTwo 						= CalculateGraphRatioToFit (graphRatioCombNLOPi07TeVMuTwo, fitTCMInvCrossSectionPi07TeV); 
// 	graphRatioCombNLODSS14Pi07TeV 						= CalculateGraphErrRatioToFit(graphRatioCombNLODSS14Pi07TeV, fitTCMInvCrossSectionPi07TeV); 
	
	TGraphAsymmErrors* graphRatioCombCombFitTot7TeVOEMC 	= (TGraphAsymmErrors*)graphCombPi0InvCrossSectionTot7TeVOEMC->Clone();
	graphRatioCombCombFitTot7TeVOEMC 						= CalculateGraphErrRatioToFit(graphRatioCombCombFitTot7TeVOEMC, fitTCMInvCrossSectionPi07TeV); 
	TGraphAsymmErrors* graphRatioCombCombFitStat7TeVOEMC 	= (TGraphAsymmErrors*)graphCombPi0InvCrossSectionStat7TeVOEMC->Clone();
	graphRatioCombCombFitStat7TeVOEMC 						= CalculateGraphErrRatioToFit(graphRatioCombCombFitStat7TeVOEMC, fitTCMInvCrossSectionPi07TeV); 
	TGraphAsymmErrors* graphRatioCombCombFitSys7TeVOEMC 	= (TGraphAsymmErrors*)graphCombPi0InvCrossSectionSys7TeVOEMC->Clone();
	graphRatioCombCombFitSys7TeVOEMC 						= CalculateGraphErrRatioToFit(graphRatioCombCombFitSys7TeVOEMC, fitTCMInvCrossSectionPi07TeV); 

	TGraphAsymmErrors* graphRatioPCMCombFitStat7TeV 		= (TGraphAsymmErrors*)graphPCMPi0InvCrossSectionStat7TeV->Clone();
	graphRatioPCMCombFitStat7TeV 						= CalculateGraphErrRatioToFit(graphRatioPCMCombFitStat7TeV, fitTCMInvCrossSectionPi07TeV); 
	TGraphAsymmErrors* graphRatioPCMCombFitSys7TeV 		= (TGraphAsymmErrors*)graphPCMPi0InvCrossSectionSys7TeV->Clone();
	graphRatioPCMCombFitSys7TeV 							= CalculateGraphErrRatioToFit(graphRatioPCMCombFitSys7TeV, fitTCMInvCrossSectionPi07TeV); 
	TGraphAsymmErrors* graphRatioPHOSCombFitStat7TeV 	= (TGraphAsymmErrors*)graphPHOSPi0InvCrossSectionStat7TeV->Clone();
	graphRatioPHOSCombFitStat7TeV 						= CalculateGraphErrRatioToFit(graphRatioPHOSCombFitStat7TeV, fitTCMInvCrossSectionPi07TeV); 
	TGraphAsymmErrors* graphRatioPHOSCombFitSys7TeV 		= (TGraphAsymmErrors*)graphPHOSPi0InvCrossSectionSys7TeV->Clone();
	graphRatioPHOSCombFitSys7TeV 						= CalculateGraphErrRatioToFit(graphRatioPHOSCombFitSys7TeV, fitTCMInvCrossSectionPi07TeV); 
	TGraphAsymmErrors* graphRatioEMCALCombFitStat7TeV 	= (TGraphAsymmErrors*)graphEMCALPi0InvCrossSectionStat7TeV->Clone();
	graphRatioEMCALCombFitStat7TeV 						= CalculateGraphErrRatioToFit(graphRatioEMCALCombFitStat7TeV, fitTCMInvCrossSectionPi07TeV); 
	TGraphAsymmErrors* graphRatioEMCALCombFitSys7TeV 	= (TGraphAsymmErrors*)graphEMCALPi0InvCrossSectionSys7TeV->Clone();
	graphRatioEMCALCombFitSys7TeV 						= CalculateGraphErrRatioToFit(graphRatioEMCALCombFitSys7TeV, fitTCMInvCrossSectionPi07TeV); 
	

	//	**********************************************************************************************************************
	//	******************************************* Ratio of Comb to Fit ****************************************
	//	**********************************************************************************************************************
	textSizeLabelsPixel = 48;
	TCanvas* canvasRatioToCombFit = new TCanvas("canvasRatioToCombFit","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioToCombFit, 0.12, 0.01, 0.01, 0.11);
    canvasRatioToCombFit->SetLogx();

		Double_t textsizeLabelsPP = 0;
		Double_t textsizeFacPP= 0;
		if (canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) <canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1()) ){
			textsizeLabelsPP = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
			textsizeFacPP = (Double_t)1./canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
		} else {
			textsizeLabelsPP = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
			textsizeFacPP = (Double_t)1./canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
		}
		cout << textsizeLabelsPP << endl;
	
	TH2F * histo2DRatioToCombFit;
	histo2DRatioToCombFit = new TH2F("histo2DRatioToCombFit","histo2DRatioToCombFit",1000,0.23,70.,1000,0.2,4.	);
	SetStyleHistoTH2ForGraphs(histo2DRatioToCombFit, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{Comb Fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP, 
							  0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
	histo2DRatioToCombFit->GetXaxis()->SetMoreLogLabels();
	histo2DRatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
// 	histo2DRatioToCombFit->GetYaxis()->SetRangeUser(-10,10);
	histo2DRatioToCombFit->GetYaxis()->SetRangeUser(0.05,2.45);
	histo2DRatioToCombFit->Draw("copy");

		DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys7TeVOEMC, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
		graphRatioCombCombFitSys7TeVOEMC->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat7TeVOEMC, markerStyleComb, markerSizeComb, colorComb , colorComb);
		graphRatioCombCombFitStat7TeVOEMC->Draw("p,same,e1");

		DrawGammaLines(0.23, 70. , 1., 1.,0.1, kGray+2);
		DrawGammaLines(0.23, 70. , 1.1, 1.1,0.1, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.9, 0.9,0.1, kGray, 7);

		TLatex *labelRatioToFitEnergy = new TLatex(0.73,0.92,collisionSystem7TeV.Data());
		SetStyleTLatex( labelRatioToFitEnergy, 0.85*textSizeLabelsPixel,4);
		labelRatioToFitEnergy->SetTextFont(43);
		labelRatioToFitEnergy->Draw();
		TLatex *labelRatioToFitPi0 = new TLatex(0.73,0.87,"#pi^{0} #rightarrow #gamma#gamma");
		SetStyleTLatex( labelRatioToFitPi0, 0.85*textSizeLabelsPixel,4);
		labelRatioToFitPi0->SetTextFont(43);
		labelRatioToFitPi0->Draw();

		
	canvasRatioToCombFit->SaveAs(Form("%s/RatioOfCombToCombFit_PP7TeV.%s",outputDir.Data(),suffix.Data()));
	
	//	**********************************************************************************************************************
	//	******************************************* Ratio of Individual meas to Fit ******************************************
	//	**********************************************************************************************************************
    
	canvasRatioToCombFit->cd();
	histo2DRatioToCombFit->Draw("copy");
	
		DrawGammaSetMarkerTGraphAsym(graphRatioPCMCombFitSys7TeV, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphRatioPCMCombFitStat7TeV, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);
		DrawGammaSetMarkerTGraphAsym(graphRatioPHOSCombFitSys7TeV, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1], widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphRatioPHOSCombFitStat7TeV, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1]);
		DrawGammaSetMarkerTGraphAsym(graphRatioEMCALCombFitSys7TeV, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphRatioEMCALCombFitStat7TeV, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);
		
		graphRatioPCMCombFitSys7TeV->Draw("E2same");
		graphRatioPHOSCombFitSys7TeV->Draw("E2same");
		graphRatioEMCALCombFitSys7TeV->Draw("E2same");
		
		graphRatioPCMCombFitStat7TeV->Draw("p,same,e");
		graphRatioPHOSCombFitStat7TeV->Draw("p,same,e");
		graphRatioEMCALCombFitStat7TeV->Draw("p,same,e");

		DrawGammaLines(0.23, 70. , 1., 1.,0.5, kGray+2);
		DrawGammaLines(0.23, 70. , 1.1, 1.1,0.5, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.9, 0.9,0.5, kGray, 7);
		
		labelRatioToFitEnergy->Draw();
		labelRatioToFitPi0->Draw();
	
		//****************************** Definition of the Legend ******************************************
		//**************** Row def ************************
		Double_t rowsLegendOnlyPi0Ratio[5] 		= {0.92,0.88,0.84,0.80,0.76};
		Double_t rowsLegendOnlyPi0RatioAbs[5] 	= {0.91,2.2,2.1,2.0,1.9};
		Double_t columnsLegendOnlyPi0Ratio[3] 	= {0.15,0.32, 0.38};
		Double_t columnsLegendOnlyPi0RatioAbs[3]= {0.15,1.04, 1.37};
		Double_t lengthBox						= 0.2/2;
		Double_t heightBox						= 0.08/2;
		//****************** first Column **************************************************
		TLatex *textPCMOnlyRatioPi0 = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[1],"PCM");
		SetStyleTLatex( textPCMOnlyRatioPi0, 0.85*textSizeLabelsPixel,4);
		textPCMOnlyRatioPi0->SetTextFont(43);
		textPCMOnlyRatioPi0->Draw();
		TLatex *textPHOSOnlyRatioPi0 = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[2],"PHOS");
		SetStyleTLatex( textPHOSOnlyRatioPi0,  0.85*textSizeLabelsPixel,4);
		textPHOSOnlyRatioPi0->SetTextFont(43);
		textPHOSOnlyRatioPi0->Draw();
		TLatex *textEMCALOnlyRatioPi0 = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[3],"EMCal");
		SetStyleTLatex( textEMCALOnlyRatioPi0,  0.85*textSizeLabelsPixel,4);
		textEMCALOnlyRatioPi0->SetTextFont(43);
		textEMCALOnlyRatioPi0->Draw();
		
		//****************** second Column *************************************************
		TLatex *textStatOnlyRatioPi0 = new TLatex(columnsLegendOnlyPi0Ratio[1],rowsLegendOnlyPi0Ratio[0] ,"stat");
		SetStyleTLatex( textStatOnlyRatioPi0, 0.85*textSizeLabelsPixel,4);
		textStatOnlyRatioPi0->SetTextFont(43);
		textStatOnlyRatioPi0->Draw();
		TLatex *textSysOnlyRatioPi0 = new TLatex(columnsLegendOnlyPi0Ratio[2] ,rowsLegendOnlyPi0Ratio[0],"syst");
		SetStyleTLatex( textSysOnlyRatioPi0, 0.85*textSizeLabelsPixel,4);
		textSysOnlyRatioPi0->SetTextFont(43);
		textSysOnlyRatioPi0->Draw();
		TMarker* markerPCMPi0OnlyRatioPi0 = CreateMarkerFromGraph(graphRatioPCMCombFitSys7TeV,columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[1],1);
		markerPCMPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[1]);
		TMarker* markerPHOSPi0OnlyRatioPi0 = CreateMarkerFromGraph(graphRatioPHOSCombFitSys7TeV, columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[2],1);
		markerPHOSPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[2]);
		TMarker* markerEMCALPi0OnlyRatioPi0 = CreateMarkerFromGraph(graphRatioEMCALCombFitSys7TeV, columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[3],1);
		markerEMCALPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[3]);
		
		TBox* boxPCMPi0OnlyRatioPi0 = CreateBoxFromGraph(graphRatioPCMCombFitSys7TeV, columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[1]- heightBox,
														 columnsLegendOnlyPi0RatioAbs[2]+ 3*lengthBox, rowsLegendOnlyPi0RatioAbs[1]+ heightBox);
		boxPCMPi0OnlyRatioPi0->Draw("l");
		TBox* boxPHOSPi0OnlyRatioPi0 = CreateBoxFromGraph(graphRatioPHOSCombFitSys7TeV, columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[2]- heightBox,
										 				  columnsLegendOnlyPi0RatioAbs[2]+ 3*lengthBox, rowsLegendOnlyPi0RatioAbs[2]+ heightBox);
		boxPHOSPi0OnlyRatioPi0->Draw("l");
		TBox* boxEMCALPi0OnlyRatioPi0 = CreateBoxFromGraph(graphRatioEMCALCombFitSys7TeV, columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[3]- heightBox,
														   columnsLegendOnlyPi0RatioAbs[2]+ 3*lengthBox, rowsLegendOnlyPi0RatioAbs[3]+ heightBox);
		boxEMCALPi0OnlyRatioPi0->Draw("l");

	canvasRatioToCombFit->SaveAs(Form("%s/RatioOfIndividualMeasToCombFit_PP7TeV.%s",outputDir.Data(),suffix.Data()));

	
		//*************************************************************************************************************
	//***************************** Comparison to Charged pions ***************************************************
	//*************************************************************************************************************

	TGraphAsymmErrors* graphCombPi0InvYieldStat7TeVOld 	= (TGraphAsymmErrors*) graphCombPi0InvCrossSectionStat7TeVOld->Clone("graphCombPi0InvYieldStat7TeVOld");
	graphCombPi0InvYieldStat7TeVOld 						= ScaleGraph(graphCombPi0InvYieldStat7TeVOld,1./xSection7TeVppINEL);
	TGraphAsymmErrors* graphCombPi0InvYieldSys7TeVOld 		= (TGraphAsymmErrors*) graphCombPi0InvCrossSectionSys7TeVOld->Clone("graphCombPi0InvYieldSys7TeVOld");
	graphCombPi0InvYieldSys7TeVOld 						= ScaleGraph(graphCombPi0InvYieldSys7TeVOld,1./xSection7TeVppINEL);

   	cout << "combined Spectrum - high Pt" << endl;
	TGraphErrors* graphChPiInvYieldHighPtStatPPHighPtCombUp 		= NULL;
	TGraphErrors* graphChPiInvYieldHighPtSystPPHighPtCombUp 		= NULL;
	TGraphErrors* graphCombPi0InvYieldStat7TeVRebinnedHighPtComb = NULL;
	TGraphErrors* graphCombPi0InvYieldSys7TeVRebinnedHighPtComb 	= NULL;
	TGraphErrors* graphRatioHighPtChPisComb7TeVOld = CalculateRatioBetweenSpectraWithDifferentBinning( graphCombPi0InvYieldStat7TeVOld, graphCombPi0InvYieldSys7TeVOld, 
																									histoChPiInvYieldHighPtStatPP, histoChPiInvYieldHighPtSystPP,  
																									kTRUE,  kTRUE, 
																									&graphCombPi0InvYieldStat7TeVRebinnedHighPtComb, &graphCombPi0InvYieldSys7TeVRebinnedHighPtComb, 
																									&graphChPiInvYieldHighPtStatPPHighPtCombUp, &graphChPiInvYieldHighPtSystPPHighPtCombUp )	;
   
	cout << "combined Spectrum - low Pt" << endl;
	TGraphErrors* graphChPiInvYieldLowPtStatPPLowPtCombUp 			= NULL;
	TGraphErrors* graphChPiInvYieldLowPtSystPPLowPtCombUp 			= NULL;
	TGraphErrors* graphCombPi0InvYieldStat7TeVRebinnedLowPtComb 	= NULL;
	TGraphErrors* graphCombPi0InvYieldSys7TeVRebinnedLowPtComb 	= NULL;
	TGraphErrors* graphRatioLowPtChPisComb7TeVOld = CalculateRatioBetweenSpectraWithDifferentBinning( graphCombPi0InvYieldStat7TeVOld, graphCombPi0InvYieldSys7TeVOld, 
																								   histoChPiInvYieldLowPtStatPP7TeV, histoChPiInvYieldLowPtSysPP7TeV,  
																								   kTRUE,  kTRUE, 
																								   &graphCombPi0InvYieldStat7TeVRebinnedLowPtComb, &graphCombPi0InvYieldSys7TeVRebinnedLowPtComb,
																								   &graphChPiInvYieldLowPtStatPPLowPtCombUp, &graphChPiInvYieldLowPtSystPPLowPtCombUp )	;
	
	cout << "combined Spectrum - low Pt CMS" << endl;
   
	TGraphErrors* graphChPiInvYieldCMSStatPPCMSCombUp 				= NULL;
	TGraphErrors* graphChPiInvYieldCMSSystPPCMSCombUp 				= NULL;
	TGraphErrors* graphCombPi0InvYieldStat7TeVRebinnedCMSComb 	= NULL;
	TGraphErrors* graphCombPi0InvYieldSys7TeVRebinnedCMSComb 	= NULL;
	TGraphErrors* graphRatioCMSChPisComb7TeVOld  = CalculateRatioBetweenSpectraWithDifferentBinning( graphCombPi0InvYieldStat7TeVOld, graphCombPi0InvYieldSys7TeVOld, 
																									  histoChPiInvYieldLowPtStat7TeVCMS, histoChPiInvYieldLowPtSys7TeVCMS,  
																								      kTRUE,  kTRUE, 
																								      &graphCombPi0InvYieldStat7TeVRebinnedCMSComb, &graphCombPi0InvYieldSys7TeVRebinnedCMSComb,
																								      &graphChPiInvYieldCMSStatPPCMSCombUp, &graphChPiInvYieldCMSSystPPCMSCombUp ) ;

   	cout << "combined Spectrum - published" << endl;
	TGraphErrors* graphChPiInvYieldPubStatPPPubCombUp 		= NULL;
	TGraphErrors* graphChPiInvYieldPubSystPPPubCombUp 		= NULL;
	TGraphErrors* graphCombPi0InvYieldStat7TeVRebinnedPubComb = NULL;
	TGraphErrors* graphCombPi0InvYieldSys7TeVRebinnedPubComb 	= NULL;
	TGraphErrors* graphRatioPubChPisComb7TeVOld = CalculateRatioBetweenSpectraWithDifferentBinning( graphCombPi0InvYieldStat7TeVOld, graphCombPi0InvYieldSys7TeVOld, 
																									histoChPiInvYieldPubStatPP, histoChPiInvYieldPubSystPP,  
																									kTRUE,  kTRUE, 
																									&graphCombPi0InvYieldStat7TeVRebinnedPubComb, &graphCombPi0InvYieldSys7TeVRebinnedPubComb, 
																									&graphChPiInvYieldPubStatPPPubCombUp, &graphChPiInvYieldPubSystPPPubCombUp )	;

   	cout << "combined Spectrum - published" << endl;
	TGraphErrors* graphChHadInvYieldPubStatPPHadCombUp 		= NULL;
	TGraphErrors* graphChHadInvYieldPubSystPPHadCombUp 		= NULL;
	TGraphErrors* graphCombPi0InvYieldStat7TeVRebinnedHadComb = NULL;
	TGraphErrors* graphCombPi0InvYieldSys7TeVRebinnedHadComb 	= NULL;
	TGraphErrors* graphRatioPubChHadsComb7TeVOld = CalculateRatioBetweenSpectraWithDifferentBinning( graphCombPi0InvYieldStat7TeVOld, graphCombPi0InvYieldSys7TeVOld, 
																									graphChargedHadronsStatPP7TeV, graphChargedHadronsSysPP7TeV,  
																									kTRUE,  kTRUE, 
																									&graphCombPi0InvYieldStat7TeVRebinnedHadComb, &graphCombPi0InvYieldSys7TeVRebinnedHadComb, 
																									&graphChHadInvYieldPubStatPPHadCombUp, &graphChHadInvYieldPubSystPPHadCombUp )	;
																									
	// ***************************************************************************************************************
	// ************************** Comparison pi0/pi+-, pi0 updated comb pp 2.76TeV ***********************
	// ***************************************************************************************************************	
	textSizeLabelsPixel = 48;
	TCanvas* canvasCompYieldPPInd = new TCanvas("canvasCompYieldPPInd","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasCompYieldPPInd,   0.12, 0.01, 0.01, 0.11);
	canvasCompYieldPPInd->SetLogx();

	TH2F * histo2DCompCombinedRatio2;
	histo2DCompCombinedRatio2 = new TH2F("histo2DCompCombinedRatio2","histo2DCompCombinedRatio2",1000,0.23,70.,1000,0.2,4.	);
	SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatio2, "#it{p}_{T} (GeV/#it{c})","#pi^{0}/#pi^{#pm}", 0.85*textsizeLabelsPP, textsizeLabelsPP, 
							  0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
	histo2DCompCombinedRatio2->GetYaxis()->SetNoExponent(kTRUE);
// 	histo2DCompCombinedRatio2->GetXaxis()->SetLabelOffset(-0.01);
	histo2DCompCombinedRatio2->GetXaxis()->SetMoreLogLabels(kTRUE);
	histo2DCompCombinedRatio2->GetXaxis()->SetNoExponent(kTRUE);
	histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.05,2.45);
	histo2DCompCombinedRatio2->DrawCopy();

		DrawGammaSetMarkerTGraphErr(graphRatioHighPtChPisComb7TeVOld, markerStyleCombHighPt, markerSizeComparison, colorCombHighPt , colorCombHighPt);
		graphRatioHighPtChPisComb7TeVOld->Draw("E1psame");
		DrawGammaSetMarkerTGraphErr(graphRatioLowPtChPisComb7TeVOld, markerStyleCombLowPt, markerSizeComparison, colorCombLowPt , colorCombLowPt);
		graphRatioLowPtChPisComb7TeVOld->Draw("E1psame");
		DrawGammaSetMarkerTGraphErr(graphRatioCMSChPisComb7TeVOld, 22, markerSizeComparison, kRed+1 , kRed+1);
		graphRatioCMSChPisComb7TeVOld->Draw("E1psame");

		TLegend* legendPi0CompChargedPionsPP = GetAndSetLegend2(0.15, 0.8, 0.9, 0.94, 0.85* textSizeLabelsPixel);
		legendPi0CompChargedPionsPP->SetNColumns(2);
		legendPi0CompChargedPionsPP->SetMargin(0.12);
		legendPi0CompChargedPionsPP->AddEntry(graphRatioLowPtChPisComb7TeVOld,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (ALICE)","p");
		legendPi0CompChargedPionsPP->AddEntry(graphRatioHighPtChPisComb7TeVOld,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (ALICE)","p");
		legendPi0CompChargedPionsPP->AddEntry(graphRatioCMSChPisComb7TeVOld,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (#pi^{#pm} from CMS)","p");
		legendPi0CompChargedPionsPP->Draw();
	
		labelRatioTheoryPP7TeV->Draw();
		
		DrawGammaLines(0., 70 , 1, 1 ,1, kGray, 1);   
   
	canvasCompYieldPPInd->Update();
	canvasCompYieldPPInd->Print(Form("%s/ComparisonChargedToNeutralOldCharged_PP7TeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));

	canvasCompYieldPPInd->cd();

	
}
	
