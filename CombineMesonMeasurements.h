/****************************************************************************************************************************
****** 	provided by Gamma Conversion Group, PWG4, 														*****
******		Ana Marin, marin@physi.uni-heidelberg.de													*****
******	   	Kathrin Koch, kkoch@physi.uni-heidelberg.de 													*****
******		Friederike Bock, friederike.bock@cern.ch													*****
*****************************************************************************************************************************/

//*****************************************************************************************************
//*********************** 7TeV variables **************************************************************
//*****************************************************************************************************
TFile* 			fileConversions;
TFile* 			fileConversionsPrelim;
TDirectory *		directoryPi07TeV;
TDirectory *		directoryEta7TeV;
Double_t 			nEvt7TeV;

TString 			collisionSystem7TeV;
const char* 		fileNameCaloPhos7TeV;
const char* 		fileNameChargedSpectra7TeV;
const char* 		fileNameChargedExpectation7TeV;
const char*		fileNameChargedSpectra2760GeV;
const char* 		fileNameSysErrEta7TeV;
const char* 		fileNamePHOSMassData7TeV;
const char* 		fileNamePHOSMassMC7TeV;
const char*		fileNameCaloPhos7TeVEta;
const char*		fileNameEMCALMassData7TeVPi0 ;
const char*		fileNameEMCALMassData7TeVEta;
const char*		fileNameEMCALMassMC7TeVPi0 ;
const char*		fileNameEMCALMassMC7TeVEta;
const char*		fileNameCaloEmcal7TeV;
const char*		fileNamePHOSEtaToPi0;
const char* 		fileNameCaloEmcalEtaToPi07TeV;
const char* 		fileNameChargedPions7TeV;
const char* 		fileNameChargedPions2760GeV;
const char* 		fileNameChargedPions900GeV;
const char* 		fileNameChargedHadrons7TeV;
const char* 	fileNameFractionPions;
const char* 		fileNameChargedHadrons2760GeV;
const char* 		fileNameChargedHadrons900GeV;


TString 			fileNameSysErrPi07TeV;
TFile* 			filePhos7TeV;
TH1D*			histoPi0Phos7TeV;
TH1D* 			histoPi0PhosSys7TeV;
TFile* 			filePhos7TeVMassData;
TFile* 			filePhos7TeVMassMC;
TFile*			filePhosEtaToPi0;
TH1D* 			histoEtaToPi0PHOS;
TH1D* 			histoEtaToPi0PHOSSys;

TH1D*			histoMassMesonPi0PHOS;
TH1D*			histoFWHMMesonPi0PHOS;
TH1D*			histoTrueMassMesonPi0PHOS;
TH1D*			histoTrueFWHMMesonPi0PHOS;
TH1D* 			histoMassMesonPi0PHOSMinusExp;
TH1D* 			histoTrueMassMesonPi0PHOSMinusExp;
TFile*			filePhos7TeVEta;
TH1D*			histoEtaPhos7TeV;
TH1D*			histoEtaPhosSys7TeV;
Double_t 			nEvtPHOS7TeV;

TFile*			fileEMCAL7TeV;
TFile*			fileEMCALEtaToPi07TeV;
TH1D*			histoPi0EMCAL7TeV;
TGraphErrors*		graphEtaEMCAL7TeV;
TGraphErrors*		graphEtaToPi0EMCAL7TeV;
TGraphErrors*		graphEtaToPi0EMCALSys7TeV;

TFile*			fileEMCAL7TeVMassData;
TFile*			fileEMCAL7TeVMassMC;
TH1D*			histoMassMesonPi0EMCAL;
TH1D*			histoFWHMMesonPi0EMCAL;
TH1D*			histoTrueMassMesonPi0EMCAL;
TH1D*			histoTrueFWHMMesonPi0EMCAL;
TFile*			fileEMCAL7TeVMassDataEta;
TFile*			fileEMCAL7TeVMassMCEta;
TH1D*			histoMassMesonEtaEMCAL;
TH1D*			histoFWHMMesonEtaEMCAL;
TH1D*			histoTrueMassMesonEtaEMCAL;
TH1D*			histoTrueFWHMMesonEtaEMCAL;
TH1D*			histoMassMesonPi0EMCALMinusExp;
TH1D*			histoTrueMassMesonPi0EMCALMinusExp;
TH1D*			histoMassMesonEtaEMCALMinusExp;
TH1D*			histoTrueMassMesonEtaEMCALMinusExp;
TH1D*			histoPi0CorrYieldMtBinShifted;
TH1D*			histoPi0CorrYieldXtBinShifted;
TH1D*			histoEtaCorrYieldMtBinShifted;
TH1D*			histoEtaCorrYieldXtBinShifted;
TH1D*			histoEtaCorrYieldMtBinShifted2760GeV;
TH1D*                   histoPi0CorrYieldMtBinShifted2760GeV;


TH1D*			histoRatioConvCaloEmcal;

TFile*			fileChargedPions7TeV;
TList*			listChargedPions7TeV;
TH1D *			histChargedPions7TeV;
TF1 *			fitChargedPions7TeV;


TFile*			fileFractionPions;
TH1D *			histFractionPions;
TFile*			fileChargedHadrons7TeV;
TH1D *			histChargedPions7TeVHighPtWideRange;
TH1D *			histChargedHadrons7TeVWideRange;
TH1D *			histChargedPions7TeVHighPt;
TF1 *			fitChargedPions7TeVHighPt;


TFile*			fileChargedPions2760GeV;
TFile*			fileChargedHadrons2760GeV;
TH1D *			histChargedPions2760GeV;
TH1D *			histChargedPions2760GeVWideRange;
TH1D *			histChargedHadrons2760GeVWideRange;
TF1 *			fitChargedPions2760GeV;
TF1 *			fitCharged2760GeVTwice;

TFile*			fileChargedPions900GeV;
TFile*			fileChargedHadrons900GeV;
TList*			listChargedPions900GeV;
TGraphAsymmErrors* 	graphChargedPions900GeV;
TGraphAsymmErrors* 	graphChargedPosPions900GeV;
TGraphAsymmErrors* 	graphChargedNegPions900GeV;
TF1 *			fitChargedPions900GeV;
TF1 *			fitCharged900GeVTwice;
TH1D *			histChargedHadrons900GeVWideRange;
TGraphAsymmErrors* 	graphChargedSpectrum900GeVWideRange;
TGraphAsymmErrors* 	graphFitChargedSpectrum900GeVWideRange;

TGraphAsymmErrors* 	graphCorrectedYieldEtaSysErrMtBinShifted;
TGraphAsymmErrors* 	graphCorrectedYieldEtaSysErrMtBinShifted2760GeV;




Double_t 			ptCharged[100];
Double_t 			chargedValue[100];
Double_t 			chargedStatError[100];
Double_t 			chargedSystError[100];		
Double_t 			chargedErrorY[100];		
Double_t 			chargedErrorX[100];
TGraphAsymmErrors* 	graphFitChargedSpect;
TGraphAsymmErrors* 	graphChargedSpectrum;
TGraphAsymmErrors* 	graphChargedSpectrumWideRange;
TGraphAsymmErrors* 	graphFitChargedSpectrumWideRange;
ifstream 			fileCharged;
TF1 *			fitCharged;
TF1 *			fitChargedTwice;
TGraphAsymmErrors* 	graphPtChargedSpectrumMC;

Double_t 			ptChargedMC[400];
Double_t 			chargedValueMC[400];
Double_t 			chargedErrorMC[400];
Int_t 			nLinesChargedMC = 		0;		
TH1D *			histoRatioConversion;
TH1D *			histoRatioChargedConvPhos;
Int_t 			npPHOS;
Double_t 			yValueCharged;
Double_t 			formerYValuePHOS = 		0;
Double_t 			xValueIntPHOS = 		0;	
Double_t 			formerYValueConv = 		0;
Double_t 			xValueConv = 			0;	
TF1 *			fitConv;
TH1D* 			histoFittingConvComp;
TH1D *			histoRatioConvCaloPhos;
TH1D *			histoRatioConvCaloPhosSys;
TH1D *			histoRatioConvConvFit;
Double_t 			xValueFitConv = 		0;	
TF1*				fitLevyConvStatErrPi07TeV;
TF1*				fitLevyConvSysErrPi07TeV;
TF1*				fitLevyConvStatErrEta7TeV;
TF1*				fitLevyConvSysErrEta7TeV;

Double_t			fitResultsPi07TeV[50];
Double_t			fitResultsEta7TeV[50];
/*Double_t			fitResultsPi07TeV[9] = {	2.724 , 	0.051 , 	0.177, 
									6.877 , 	0.058 , 	0.147,
									0.142 ,	0.002 , 	0.006};
Double_t 			fitResultsEta7TeV[9] = {	0.228 , 	0.040 , 	0.025 ,
									7.199 , 	0.597 , 	0.320 , 
									0.244 , 	0.024 , 	0.012};*/
//************** PHOS variables 7TeV ********************************
Double_t 			xValuePHOS[100];
Double_t 			xErrorPHOS[100];
Double_t 			yValuePHOS[100];
Double_t 			yValueRatioPHOSConv[100];
Double_t 			yErrorSysPHOS[100];
Double_t 			yErrorSysRatioPHOSConv[100];
const TH1D* 		histoForPHOSGraph1;
const TH1D* 		histoForPHOSGraph2;
TGraphErrors* 		graphSysErrPi0PHOS;
TGraphErrors* 		graphSysErrPi0RatioPHOSConv;

TH1D*			histoFittingConvCompEta;
TF1*				fitConvEta;
TH1D*			histoRatioConvCaloPhosEta;
TH1D*			histoRatioConvCaloPhosEtaSys;
TH1D*			histoRatioConvConvFitEta;
TH1D*			histoRatioConvCaloEmcalEta;
TGraphErrors*		graphSysErrEtaPHOS;
TGraphErrors*		graphSysErrEtaRatioPHOSConv;
TGraphErrors*		graphStaErrInvCrossSectionPi0PHOS;

TH1D*			histoInvCrossSectionPi0PHOS;
TH1D*			histoInvCrossSectionPi0SysPHOS;
TGraphAsymmErrors*		graphSysErrInvCrossSectionPi0PHOS;
TGraphAsymmErrors*		graphSysErrInvCrossSectionPi0PHOSUnscaled;
TH1D*			histoInvCrossSectionPi0EMCAL;

TH1D*			histoInvCrossSectionEtaPHOS;
TH1D*			histoInvCrossSectionEtaSysPHOS;
TGraphAsymmErrors*		graphSysErrInvCrossSectionEtaPHOS;
TGraphAsymmErrors*		graphSysErrInvCrossSectionEtaPHOSUnscaled;
TH1D*			histoInvCrossSectionEtaEMCAL;

TGraphAsymmErrors* 	graphInvCrossSectionPi0Comb7TeV;
TGraphAsymmErrors* 	graphInvCrossSectionPi0Comb7TeVStatErr;
TGraphAsymmErrors* 	graphInvCrossSectionPi0Comb7TeVSysErr;
TGraphAsymmErrors* 	graphInvCrossSectionPi0Comb7TeVUnscaled;
TGraphAsymmErrors* 	graphInvCrossSectionPi0Comb7TeVStatErrUnscaled;
TGraphAsymmErrors* 	graphInvCrossSectionPi0Comb7TeVSysErrUnscaled;

TGraphAsymmErrors* 	graphInvCrossSectionEtaComb7TeV;
TGraphAsymmErrors* 	graphInvCrossSectionEtaComb7TeVStatErr;
TGraphAsymmErrors* 	graphInvCrossSectionEtaComb7TeVSysErr;
TGraphAsymmErrors* 	graphInvCrossSectionEtaComb7TeVUnscaled;
TGraphAsymmErrors* 	graphInvCrossSectionEtaComb7TeVStatErrUnscaled;
TGraphAsymmErrors* 	graphInvCrossSectionEtaComb7TeVSysErrUnscaled;

TGraph*			graphNLOMuHalfPi07TeV;
TGraph*			graphNLOMuHalf7TeV;
TGraphErrors*		graphFitNLOMuHalfPi07TeV;
TF1*				fitNLOMuHalfPi07TeV;
TGraph*			graphNLOMuOnePi07TeV;
TGraphErrors*		graphFitNLOMuOnePi07TeV;
TF1*				fitNLOMuOnePi07TeV;
TGraph*			graphNLOMuTwoPi07TeV;
TGraphErrors*		graphFitNLOMuTwoPi07TeV;
TF1*				fitNLOMuTwoPi07TeV;
TGraph*			graphNLOMuHalfEta7TeV;
TGraphErrors*		graphFitNLOMuHalfEta7TeV;
TF1*				fitNLOMuHalfEta7TeV;
TGraph*			graphNLOMuOneEta7TeV;
TGraphErrors*		graphFitNLOMuOneEta7TeV;
TF1*				fitNLOMuOneEta7TeV;
TGraph*			graphNLOMuTwoEta7TeV;
TGraphErrors*		graphFitNLOMuTwoEta7TeV;
TF1*				fitNLOMuTwoEta7TeV;
TF1* 			fitInvCrossSectionEta7TeV;
TF1* 			fitInvCrossSectionPi0;

TGraphAsymmErrors*			graphRatioCombPHOSPi0;
TGraphAsymmErrors*			graphRatioCombPHOSPi0Sys;
TGraphAsymmErrors*			graphRatioCombConvPi0;
TH1D*			histoRatioCombEMCALPi0;
TGraphErrors*		graphSysErrRatioCombPHOSPi0;
TGraphAsymmErrors*	graphSysErrRatioCombConvPi0;

TH1D*			histoRatioCombPHOSEta7TeV;
TH1D*			histoRatioCombPHOSEta7TeVSys;
TH1D*			histoRatioCombConvEta7TeV;
TGraphErrors*		graphRatioCombEMCALEta7TeV;
TGraphErrors*		graphSysErrRatioCombPHOSEta7TeV;
TGraphAsymmErrors*	graphSysErrRatioCombConvEta7TeV;

TH1D*			histoFitInvCrossSectionEta7TeV;

TH1D*			histoNLOMuHalfPi07TeV;
TH1D*			histoNLOMuOnePi07TeV;
TH1D*			histoNLOMuTwoPi07TeV;
TH1D*			histoNLOMuHalfEta7TeV;
TH1D*			histoNLOMuOneEta7TeV;
TH1D*			histoNLOMuTwoEta7TeV;
TGraph*			graphRatioCombNLOPi07TeVMuHalf;
TGraph*			graphRatioCombNLOPi07TeVMuOne;
TGraph*			graphRatioCombNLOPi07TeVMuTwo;
TGraph*			graphRatioCombNLOBKKPi07TeVMuTwo;
TGraph*			graphRatioCombNLODSSPi07TeVMuTwo;
TGraph*			graphRatioCombNLOEta7TeVMuHalf;
TGraph*			graphRatioCombNLOEta7TeVMuOne;
TGraph*			graphRatioCombNLOEta7TeVMuTwo;

TGraphAsymmErrors*	graphRatioCombCombFitEta7TeV;
TGraphAsymmErrors*	graphRatioCombCombFitEta7TeVStat;
TGraphAsymmErrors*	graphRatioCombCombFitEta7TeVSys;
TGraphAsymmErrors*	graphRatioCombCombFit;
TGraphAsymmErrors*	graphRatioCombCombFitStat;
TGraphAsymmErrors*	graphRatioCombCombFitSys;
TGraphAsymmErrors*	graphCombinedEtaToPi0;
TGraphAsymmErrors*	graphCombinedEtaToPi0WOXErrors;
TGraphAsymmErrors*	graphStatErrCombinedEtaToPi0;
TGraphAsymmErrors*	graphSysErrCombinedEtaToPi0;
Double_t*			relSystErrorPi07TeVDown;
Double_t*			relSystErrorPi07TeVUp;
Double_t*			relSystErrorEta7TeVDown;
Double_t*			relSystErrorEta7TeVUp;
Double_t*			relSystErrorEtaPi07TeVDown;
Double_t*			relSystErrorEtaPi07TeVUp;

TH1D*			histoPi0ToChargedPhojet7TeV;
TH1D*			histoPi0ToChargedPythia7TeV;
Int_t 			nPointsEtaPi07TeV;

Double_t xPtLimits7TeV[37] =  {0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,7.0,8.0,9.0,10.,11.,12.,13.,14.,16.,18.,20.,25.};	
Double_t xPtLimits7TeVNewPHOS[34] =  {0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,7.0,8.0,10.,12.,14.,16.,18.,20.,25.};	
Double_t xPtLimitsEta7TeV[14]       =  {0.4,0.7,1.0,1.4,1.8,2.2,2.6,3.0,3.5,4.0,6.0,8.0,10.,15.};
//*****************************************************************************************************
//********************** 2.76 TeV variables ***********************************************************
//*****************************************************************************************************
const char* 		fileNameCaloPhos2760GeV;
TDirectory*       directoryPHOSPi02760GeV;
TDirectory* 		directoryPi02760GeV; 
TDirectory* 		directoryEta2760GeV; 

TGraph*     graphEtaToPi0NLOMuHalf7TeV;
TGraph*     graphEtaToPi0NLOMuOne7TeV;
TGraph*     graphEtaToPi0NLOMuTwo7TeV;
TGraph*     graphEtaToPi0NLOMuHalf2760GeV;
TGraph*     graphEtaToPi0NLOMuOne2760GeV;
TGraph*     graphEtaToPi0NLOMuTwo2760GeV;



TFile* 			filePhos2760GeV;
TH1D*			histoPi0Phos2760GeV ;
TH1D*			histoPi0PhosSys2760GeV;

TH1D*			histoNumberOfEvents2760GeV;
TH1D* 			histoMassMesonPi02760GeV;
TH1D* 			histoMassMesonPi0MinusExp2760GeV;
TH1D* 			histoFWHMMesonPi0MeV2760GeV;
TH1D* 			histoTrueMassMesonPi0MinusExp2760GeV;
TH1D* 			histoTrueFWHMMesonPi0MeV2760GeV;
TH1D* 			histoTrueMassMesonPi02760GeV;
TH1D* 			histoAccPi02760GeV;
TH1D* 			histoTrueEffPtPi02760GeV;
TH1D* 			histoPi0CorrYieldBinShifted2760GeV;
TGraphAsymmErrors* 	graphCorrectedYieldPi0SysErr2760GeV;
TGraphAsymmErrors* 	graphCorrectedYieldPi0SysErrBinShifted2760GeV;
TH1D*			histoInvCrossSectionPi02760GeV;
TGraphAsymmErrors*	graphInvCrossSectionSysPi02760GeV;
TGraphAsymmErrors*	graphInvCrossSectionSysPi02760GeVUnscaled;
TGraphAsymmErrors*	graphInvCrossSectionSysAPi02760GeV;
TGraphErrors* 		graphInvCrossSectionStaPi02760GeV;
TGraph*			graphNLOCalcMuHalf2760GeV;
TGraph*			graphNLOCalcMuOne2760GeV;
TGraph*			graphNLOCalcMuTwo2760GeV;
TGraph*        graphNLOBKKCalcMuTwo7TeV;
TGraph*        graphNLODSSCalcMuTwo7TeV;
TGraph*			graphNLOBKKCalcMuTwo2760GeV;
TGraph*			graphNLODSSCalcMuTwo2760GeV;
Double_t 			nEvt2760GeV;
Double_t			maxPtPi02760GeV;
Double_t			maxPtEta2760GeV;

TH1D* 			histoCorrectedYieldEta2760GeV;
TH1D* 			histoUnCorrectedYieldEta2760GeV;
TH1D* 			histoMassMesonEta2760GeV;
TH1D* 			histoMassMesonEtaMinusExp2760GeV;
TH1D* 			histoFWHMMesonEtaMeV2760GeV;
TH1D* 			histoTrueMassMesonEta2760GeV;
TH1D* 			histoTrueMassMesonEtaMinusExp2760GeV;
TH1D* 			histoTrueFWHMMesonEtaMeV2760GeV;
TH1D* 			histoEtaCorrYieldBinShifted2760GeV;
TH1D* 			histoInvCrossSectionEta2760GeV;
TGraphAsymmErrors*	graphInvCrossSectionSysAEta2760GeV;
TGraphAsymmErrors*	graphInvCrossSectionSysEta2760GeV;
TGraphAsymmErrors*	graphInvCrossSectionSysEta2760GeVUnscaled;
TGraphAsymmErrors*	graphCorrectedYieldEtaSysErr2760GeV;
TGraphAsymmErrors*	graphCorrectedYieldEtaSysErrBinShifted2760GeV;
TH1D* 			histoRatioEtaPi02760GeV;
TGraphAsymmErrors*	graphSystErrRatio2760GeV;
TGraph* 			graphNLOCalcEtaMuHalf2760GeV;
TGraph* 			graphNLOCalcEtaMuOne2760GeV;
TGraph* 			graphNLOCalcEtaMuTwo2760GeV;

Double_t*			relSystErrorPi02760GeVDown;
Double_t*			relSystErrorPi02760GeVUp;
Double_t*			relSystErrorEtaPi02760GeVDown;
Double_t*			relSystErrorEtaPi02760GeVUp;
Double_t*			relSystErrorEta2760GeVDown;
Double_t*			relSystErrorEta2760GeVUp;
Int_t 			nPointsPi02760GeV;
Int_t 			nPointsEta2760GeV;
TGraphAsymmErrors* 	graphRatioEtaPi0ComplErr2760GeV;
TGraphAsymmErrors* 	graphRatioEtaPi0StatErr2760GeV;
TGraphAsymmErrors* 	graphRatioEtaPi0SysErr2760GeV;
TString 			collisionSystem2760GeV;
TH1D*			histoPi0ToChargedPhojet2760GeV;
TH1D*			histoPi0ToChargedPythia2760GeV;
TH1D*			histoEtaToPi0Phojet2760GeV;
TH1D*			histoEtaToPi0Pythia2760GeV;

TGraph*			graphNLOMuHalfPi02760GeV;
TGraph*			graphNLOMuOnePi02760GeV;
TGraph*			graphNLOMuTwoPi02760GeV;
TGraphAsymmErrors* 	graphInvCrossSectionPi0Comb2760GeV;
TGraphAsymmErrors* 	graphInvCrossSectionPi0Comb2760GeVUnscaled;
TGraphAsymmErrors* 	graphInvCrossSectionPi0Comb2760GeVStatErr;
TGraphAsymmErrors* 	graphInvCrossSectionPi0Comb2760GeVStatErrUnscaled;
TGraphAsymmErrors* 	graphInvCrossSectionPi0Comb2760GeVSysErr;
TGraphAsymmErrors* 	graphInvCrossSectionPi0Comb2760GeVSysErrUnscaled;

TF1* 			fitInvCrossSectionPi02760GeV;
TH1D*			histoFitInvCrossSectionPi02760GeV;
TH1D*			histoFitInvCrossSectionFinalPi02760GeV;
TH1D*			histoRatioCombConvPi02760GeV;
TGraphAsymmErrors*		graphSysErrEtaToPi0PHOS2760GeV;
TGraphAsymmErrors*	graphSysErrRatioCombConvPi02760GeV;

TH1D*			histoNLOMuHalfPi02760GeV;
TH1D*			histoNLOMuOnePi02760GeV;
TH1D*			histoNLOMuTwoPi02760GeV;
TGraph*			graphRatioCombNLOPi02760GeVMuHalf;
TGraph*			graphRatioCombNLOPi02760GeVMuOne;
TGraph*			graphRatioCombNLOPi02760GeVMuTwo;
TGraphAsymmErrors*	graphRatioCombCombFit2760GeV;
TGraphAsymmErrors*	graphRatioCombCombFit2760GeVStat;
TGraphAsymmErrors*	graphRatioCombCombFit2760GeVSys;

TFile*			fileCharged2760GeV;
TH1D*			histoCharged2760GeV;
TF1*				fitCharged2760GeV;
TH1D*			histoRatioConversion2760GeV;
TH1D*			histoRatioChargedConvPhos2760GeV;
Int_t 			npPHOS2760GeV;

TGraph*			graphNLOMuHalfEta2760GeV;
TGraph*			graphNLOMuOneEta2760GeV;
TGraph*			graphNLOMuTwoEta2760GeV;
TGraphAsymmErrors* 	graphInvCrossSectionEtaComb2760GeV;
TGraphAsymmErrors* 	graphInvCrossSectionEtaComb2760GeVStatErr;
TGraphAsymmErrors* 	graphInvCrossSectionEtaComb2760GeVSysErr;
TGraphAsymmErrors* 	graphInvCrossSectionEtaComb2760GeVUnscaled;
TGraphAsymmErrors* 	graphInvCrossSectionEtaComb2760GeVStatErrUnscaled;
TGraphAsymmErrors* 	graphInvCrossSectionEtaComb2760GeVSysErrUnscaled;
TF1* 			fitInvCrossSectionEta2760GeV;
TH1D*			histoFitInvCrossSectionEta2760GeV;
TH1D*			histoFitFinalInvCrossSectionEta2760GeV;
TH1D*			histoRatioCombPHOSEta2760GeV;
TH1D*			histoRatioCombPHOSEta2760GeVSys;
TH1D*			histoRatioCombConvEta2760GeV;
TGraphErrors*		graphSysErrRatioCombPHOSEta2760GeV;
TGraphAsymmErrors*	graphSysErrRatioCombConvEta2760GeV;
TGraphAsymmErrors*	graphInvCrossSectionEtaComb2760GeVScaled;
TH1D*			histoNLOMuHalfEta2760GeV;
TH1D*			histoNLOMuOneEta2760GeV;
TH1D*			histoNLOMuTwoEta2760GeV;
TGraph*			graphRatioCombNLOEta2760GeVMuHalf;
TGraph*			graphRatioCombNLOEta2760GeVMuOne;
TGraph*			graphRatioCombNLOEta2760GeVMuTwo;
TGraphAsymmErrors*	graphRatioCombCombFitEta2760GeV;
TGraphAsymmErrors*	graphRatioCombCombFitEta2760GeVStat;
TGraphAsymmErrors*	graphRatioCombCombFitEta2760GeVSys;

TH1D* 			histoInvCrossSectionPi0PHOS2760GeV;
TH1D*			histoInvCrossSectionPi0SysPHOS2760GeV;
TGraphAsymmErrors*		graphSysErrInvCrossSectionPi0PHOS2760GeV;
TGraphAsymmErrors*		graphSysErrInvCrossSectionPi0PHOS2760GeVUnscaled;
TH1D*			histoRatioCombPHOSPi02760GeV;
TGraphErrors*		graphSysErrRatioCombPHOSPi02760GeV;
TH1D*			histoRatioCombPHOSPi0Sys2760GeV;
TF1*				fitLevyConvStatErrPi02760GeV;
TF1*				fitLevyConvSysErrPi02760GeV;
TF1*				fitLevyConvStatErrEta2760GeV;
TF1*				fitLevyConvSysErrEta2760GeV;
Double_t			fitResultsPi02760GeV[50];
Double_t			fitResultsEta2760GeV[50];

/*Double_t			fitResultsPi02760GeV[9] = {	2.003 , 	0.097 , 	0.158 , 
										7.267 , 	0.171 , 	0.240,
										0.139 , 	0.005 , 	0.007};
Double_t 			fitResultsEta2760GeV[9] = {	0.226 , 	0.116 , 	0.084 ,
										6.324 , 	1.297 , 	0.381 , 
										0.188 , 	0.053 , 	0.008 };*/
									
Double_t xPtLimits2760GeV[20] =  {0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,3.0,3.5,4.0,5.,6.,8.,11.,15.} ; 
Double_t xPtLimits2760GeVNew[21] =  {0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,3.0,3.5,4.0,5.,6.,8.,10,12.,15.} ; 

//*****************************************************************************************************
//********************** 900 GeV variables ************************************************************
//*****************************************************************************************************

TDirectory* 		directoryPi0900GeV; 
TDirectory* 		directoryEta900GeV; 

TH1D*			histoNumberOfEvents900GeV;
TH1D* 			histoMassMesonPi0900GeV;
TH1D* 			histoMassMesonPi0MinusExp900GeV;
TH1D* 			histoFWHMMesonPi0MeV900GeV;
TH1D* 			histoTrueMassMesonPi0MinusExp900GeV;
TH1D* 			histoTrueFWHMMesonPi0MeV900GeV;
TH1D* 			histoTrueMassMesonPi0900GeV;
TH1D* 			histoAccPi0900GeV;
TH1D* 			histoTrueEffPtPi0900GeV;
TH1D* 			histoPi0CorrYieldBinShifted900GeV;
TGraphAsymmErrors* 	graphCorrectedYieldPi0SysErr900GeV;
TGraphAsymmErrors* 	graphCorrectedYieldPi0SysErrBinShifted900GeV;
TH1D*			histoInvCrossSectionPi0900GeV;
TGraphAsymmErrors*	graphInvCrossSectionSysPi0900GeV;
TGraphAsymmErrors*	graphInvCrossSectionSysPi0900GeVUnscaled;
TGraphAsymmErrors*	graphInvCrossSectionSysPi0Unscaled;
TGraphAsymmErrors*	graphInvCrossSectionSysAPi0900GeV;
TGraphErrors* 		graphInvCrossSectionStaPi0900GeV;
TGraph*			graphNLOCalcMuHalf900GeV;
TGraph*			graphNLOCalcMuOne900GeV;
TGraph*			graphNLOCalcMuTwo900GeV;
TGraph*			graphNLOBKKCalcMuTwo900GeV;
TGraph*			graphNLODSSCalcMuTwo900GeV;
Double_t 			nEvt900GeV;

TH1D* 			histoCorrectedYieldEta900GeV;
TH1D* 			histoUnCorrectedYieldEta900GeV;
TH1D* 			histoMassMesonEta900GeV;
TH1D* 			histoMassMesonEtaMinusExp900GeV;
TH1D* 			histoFWHMMesonEtaMeV900GeV;
TH1D* 			histoTrueMassMesonEta900GeV;
TH1D* 			histoTrueMassMesonEtaMinusExp900GeV;
TH1D* 			histoTrueFWHMMesonEtaMeV900GeV;
TH1D* 			histoEtaCorrYieldBinShifted900GeV;
TH1D* 			histoInvCrossSectionEta900GeV;
TGraphAsymmErrors*	graphInvCrossSectionSysAEta900GeV;
TGraphAsymmErrors*	graphInvCrossSectionSysEta900GeV;
TGraphAsymmErrors*	graphInvCrossSectionSysEta900GeVUnscaled;
TGraphAsymmErrors*	graphInvCrossSectionSysEtaUnscaled;
TGraphAsymmErrors*	graphCorrectedYieldEtaSysErr900GeV;
TGraphAsymmErrors*	graphCorrectedYieldEtaSysErrBinShifted900GeV;
TH1D* 			histoRatioEtaPi0900GeV;
TGraphAsymmErrors*	graphSystErrRatio900GeV;
TGraph* 			graphNLOCalcEtaMuHalf900GeV;
TGraph* 			graphNLOCalcEtaMuOne900GeV;
TGraph* 			graphNLOCalcEtaMuTwo900GeV;

TGraph*        graphNLOCalcEtaMuHalf7TeV;
TGraph*        graphNLOCalcEtaMuOne7TeV;
TGraph*        graphNLOCalcEtaMuTwo7TeV;
TGraph*        graphNLOCalcMuHalf7TeV;
TGraph*        graphNLOCalcMuOne7TeV;
TGraph*        graphNLOCalcMuTwo7TeV;


TString 			collisionSystem900GeV;

const char* 		fileNameCaloPhos900GeV;
const char* 		fileNameChargedSpectra900GeV;
const char* 		fileNameChargedExpectation900GeV;
const char* 		fileNameSysErrEta900GeV;

TString 			fileNameSysErrPi0900GeV;
TFile* 			filePhos900GeV;
TH1D*			histoPi0Phos900GeV;
TH1D* 			histoPi0PhosSys900GeV;
Double_t 			nEvtPHOS900GeV;

double			p7853_d4x1y1_xval[] = 		{ 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 
										0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.05, 1.15, 
										1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.1, 2.3, 
										2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 3.9, 4.25, 4.75, 
										5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5 };
double 			p7853_d4x1y1_xerrminus[] = 	{ 0.024999999999999994, 0.024999999999999994, 0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.024999999999999967,
										0.024999999999999967, 0.025000000000000022, 0.02499999999999991, 
										0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.02499999999999991, 0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.050000000000000044, 0.04999999999999982, 
										0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.04999999999999982, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.10000000000000009, 0.09999999999999964, 
										0.10000000000000009, 0.10000000000000009, 0.10000000000000009, 0.10000000000000009, 0.09999999999999964, 0.10000000000000009, 0.10000000000000009, 0.10000000000000009, 0.25, 0.25, 
0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5 };
double 			p7853_d4x1y1_xerrplus[] = 	{ 0.025000000000000022, 0.024999999999999994, 0.024999999999999967, 0.024999999999999967, 0.025000000000000022, 0.025000000000000022,
										0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 
										0.025000000000000022, 0.02499999999999991, 0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.025000000000000022, 0.02499999999999991, 0.025000000000000022, 0.050000000000000044, 0.050000000000000044, 
										0.050000000000000044, 0.04999999999999982, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.050000000000000044, 0.04999999999999982, 0.050000000000000044, 0.10000000000000009, 0.10000000000000009, 
										0.10000000000000009, 0.09999999999999964, 0.10000000000000009, 0.10000000000000009, 0.10000000000000009, 0.10000000000000009, 0.09999999999999964, 0.10000000000000009, 0.25, 0.25, 
										0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5 };
double 			p7853_d4x1y1_yval[] = 		{ 4.976, 3.928, 3.072, 2.362, 1.828, 1.395, 1.075, 0.8468, 0.6641, 
										0.5259, 0.4136, 0.3342, 0.271, 0.2165, 0.1771, 0.1444, 0.1213, 0.09083, 0.06274, 
										0.0449, 0.03195, 0.02352, 0.01726, 0.0127, 0.009184, 0.007037, 0.005498, 0.003728, 0.002363, 
										0.001397, 9.557E-4, 6.457E-4, 3.951E-4, 2.868E-4, 1.961E-4, 1.166E-4, 8.715E-5, 5.308E-5, 2.786E-5, 
										1.336E-5, 5.917E-6, 4.014E-6, 1.828E-6, 1.182E-6, 4.854E-7, 2.656E-7 };
double 			p7853_d4x1y1_yerrminus[] =	{ 0.3545095203235027, 0.17456230979223436, 0.12847178678604887, 0.08945389874119516, 0.06140032573203501, 0.046389654018972805,
										0.03337663853655728, 0.02660244349679179, 0.020878697277368626, 
										0.01626068879229905, 0.01282653499585917, 0.010316006979447038, 0.008297590011563598, 0.006694027188471825, 0.005481788029466298, 0.004464302857109943, 0.003764306044943742, 0.0027965872058636045, 0.0019190101615155664, 
										0.0013947401191619893, 9.99599919967984E-4, 7.446475676452586E-4, 5.909314681077663E-4, 4.4721359549995795E-4, 3.326695056659086E-4, 2.619255619446105E-4, 2.105445321066306E-4, 1.444749113168096E-4, 9.625487000666512E-5, 
										6.184658438426491E-5, 4.4921375758095394E-5, 3.299484808269315E-5, 2.2901964981197576E-5, 1.8253766734567417E-5, 1.4006427096158393E-5, 1.0031948963187562E-5, 8.349874250550123E-6, 4.087542048713383E-6, 2.693120866207085E-6, 
										1.722556240010758E-6, 1.0535463919543363E-6, 8.343248767716326E-7, 5.337312057581044E-7, 2.9076450952618E-7, 1.731973440904912E-7, 1.1994569604616916E-7 };
double 			p7853_d4x1y1_yerrplus[] = 	{ 0.3545095203235027, 0.17456230979223436, 0.12847178678604887, 0.08945389874119516, 0.06140032573203501, 0.046389654018972805,
										0.03337663853655728, 0.02660244349679179, 0.020878697277368626, 
										0.01626068879229905, 0.01282653499585917, 0.010316006979447038, 0.008297590011563598, 0.006694027188471825, 0.005481788029466298, 0.004464302857109943, 0.003764306044943742, 0.0027965872058636045, 0.0019190101615155664, 
										0.0013947401191619893, 9.99599919967984E-4, 7.446475676452586E-4, 5.909314681077663E-4, 4.4721359549995795E-4, 3.326695056659086E-4, 2.619255619446105E-4, 2.105445321066306E-4, 1.444749113168096E-4, 9.625487000666512E-5, 
										6.184658438426491E-5, 4.4921375758095394E-5, 3.299484808269315E-5, 2.2901964981197576E-5, 1.8253766734567417E-5, 1.4006427096158393E-5, 1.0031948963187562E-5, 8.349874250550123E-6, 4.087542048713383E-6, 2.693120866207085E-6, 
										1.722556240010758E-6, 1.0535463919543363E-6, 8.343248767716326E-7, 5.337312057581044E-7, 2.9076450952618E-7, 1.731973440904912E-7, 1.1994569604616916E-7 };
int 				p7853_d4x1y1_numpoints = 	46;

TGraphAsymmErrors* 	graphFitChargedSpect900GeV;
TGraphAsymmErrors* 	graphChargedSpectrum900GeV;
TF1* 			fitCharged900GeV;
Double_t 			maxPtPi0900GeV;
TH1D* 			histoRatioConversion900GeV;
TH1D*			histoRatioChargedConvPhos900GeV;
TH1D*  			histoFittingConvComp900GeV;
TF1*				fitConv900GeV;
TH1D*			histoRatioConvCaloPhos900GeV;
TH1D*			histoRatioConvCaloPhos900GeVSys;
TH1D*			histoRatioConvConvFit900GeV;
TH1D*			histoForPHOSGraph900GeV1;
TH1D*			histoForPHOSGraph900GeV2;
TGraphErrors*		graphSysErrPi0PHOS900GeV;
TGraphErrors*		graphSysErrPi0RatioPHOSConv900GeV;
Int_t 			nPointsPi0900GeV;
TGraphAsymmErrors*	graphCorrectedYieldPi0SysRatioFitConv900GeV;

TH1D*			histoInvCrossSectionPi0PHOS900GeV;
TH1D*			histoInvCrossSectionPi0SysPHOS900GeV;
TGraphAsymmErrors*		graphSysErrInvCrossSectionPi0PHOS900GeV;
TGraphAsymmErrors*		graphSysErrInvCrossSectionPi0PHOS900GeVUnscaled;

TGraphAsymmErrors* 	graphInvCrossSectionPi0Comb900GeV;
TGraphAsymmErrors* 	graphInvCrossSectionPi0Comb900GeVUnscaled;
TGraphAsymmErrors* 	graphInvCrossSectionPi0Comb900GeVStatErr;
TGraphAsymmErrors* 	graphInvCrossSectionPi0Comb900GeVStatErrUnscaled;
TGraphAsymmErrors* 	graphInvCrossSectionPi0Comb900GeVSysErr;
TGraphAsymmErrors* 	graphInvCrossSectionPi0Comb900GeVSysErrUnscaled;

Int_t 			nlinesNLO900GeV = 0;
TGraph*			graphNLOMuHalfPi0900GeV;
TGraph*			graphNLOMuHalf900GeV;
TGraphErrors*		graphFitNLOMuHalfPi0900GeV;
TF1*				fitNLOMuHalfPi0900GeV;
TGraph*			graphNLOMuOnePi0900GeV;
TGraphErrors*		graphFitNLOMuOnePi0900GeV;
TF1*				fitNLOMuOnePi0900GeV;
TGraph*			graphNLOMuTwoPi0900GeV;
TGraphErrors*		graphFitNLOMuTwoPi0900GeV;
TF1*				fitNLOMuTwoPi0900GeV;

TF1* 			fitInvCrossSectionPi0900GeV;
TH1D*			histoFitInvCrossSectionPi0900GeV;
TH1D*			histoFitInvCrossSectionFinalPi0900GeV;
TGraphAsymmErrors*			graphRatioCombPHOSPi0900GeV;
TGraphAsymmErrors*			graphRatioCombPHOSPi0900GeVSys;
TGraphAsymmErrors*			graphRatioCombConvPi0900GeV;
TGraphErrors*		graphSysErrRatioCombPHOSPi0900GeV;
TGraphAsymmErrors*		graphSysErrEtaToPi0PHOS;
TGraphAsymmErrors*	graphSysErrRatioCombConvPi0900GeV;
TH1D*			histoNLOMuHalfPi0900GeV;
TH1D*			histoNLOMuOnePi0900GeV;
TH1D*			histoNLOMuTwoPi0900GeV;
TGraph*			graphRatioCombNLOPi0900GeVMuHalf;
TGraph*			graphRatioCombNLOPi0900GeVMuOne;
TGraph*			graphRatioCombNLOPi0900GeVMuTwo;
TGraph*			graphRatioCombNLOBKKPi0900GeVMuTwo;
TGraph*			graphRatioCombNLODSSPi0900GeVMuTwo;
TGraphAsymmErrors*	graphRatioCombCombFit900GeV;
TGraphAsymmErrors*	graphRatioCombCombFit900GeVStat;
TGraphAsymmErrors*	graphRatioCombCombFit900GeVSys;

Double_t*	 		relSystErrorWOMaterialPi0Down900GeV;
Double_t*			relSystErrorWOMaterialPi0Up900GeV;
Double_t*	 		relSystErrorPi0900GeVDown;
Double_t*			relSystErrorPi0900GeVUp;
Double_t*	 		relSystErrorEta900GeVDown;
Double_t*			relSystErrorEta900GeVUp;
Double_t*	 		relSystErrorEtaPi0900GeVDown;
Double_t*			relSystErrorEtaPi0900GeVUp;
TGraphAsymmErrors* 	graphRatioEtaPi0ComplErr900GeV;
TGraphAsymmErrors* 	graphRatioEtaPi0StatErr900GeV;
TGraphAsymmErrors* 	graphRatioEtaPi0SysErr900GeV;
Int_t 			nPointsEta900GeV;
TH1D*			histoPi0ToChargedPhojet900GeV;
TH1D*			histoPi0ToChargedPythia900GeV;

TGraph*			graphNLOMuHalfEta900GeV;
TGraph*			graphNLOMuOneEta900GeV;
TGraph*			graphNLOMuTwoEta900GeV;
TGraphAsymmErrors* 	graphInvCrossSectionEtaComb900GeV;
TGraphAsymmErrors* 	graphInvCrossSectionEtaComb900GeVUnscaled;
TGraphAsymmErrors* 	graphInvCrossSectionEtaComb900GeVStatErr;
TGraphAsymmErrors* 	graphInvCrossSectionEtaComb900GeVStatErrUnscaled;
TGraphAsymmErrors* 	graphInvCrossSectionEtaComb900GeVSysErr;
TGraphAsymmErrors* 	graphInvCrossSectionEtaComb900GeVSysErrUnscaled;

TF1* 			fitInvCrossSectionEta900GeV;
TH1D*			histoFitInvCrossSectionEta900GeV;
TH1D*			histoFitInvCrossSectionFinalEta900GeV;
TH1D*			histoRatioCombPHOSEta900GeV;
TH1D*			histoRatioCombPHOSEta900GeVSys;
TH1D*			histoRatioCombConvEta900GeV;
TGraphErrors*		graphSysErrRatioCombPHOSEta900GeV;
TGraphAsymmErrors*	graphSysErrRatioCombConvEta900GeV;
TGraphAsymmErrors*	graphInvCrossSectionEtaComb900GeVScaled;
TH1D*			histoNLOMuHalfEta900GeV;
TH1D*			histoNLOMuOneEta900GeV;
TH1D*			histoNLOMuTwoEta900GeV;
TGraph*			graphRatioCombNLOEta900GeVMuHalf;
TGraph*			graphRatioCombNLOEta900GeVMuOne;
TGraph*			graphRatioCombNLOEta900GeVMuTwo;
TGraphAsymmErrors*	graphRatioCombCombFitEta900GeV;
TGraphAsymmErrors*	graphRatioCombCombFitEta900GeVStat;
TGraphAsymmErrors*	graphRatioCombCombFitEta900GeVSys;

TF1*				fitLevyConvStatErrPi0900GeV;
TF1*				fitLevyConvSysErrPi0900GeV;
Double_t			fitResultsPi0900GeV[50];

// Double_t			fitResultsPi0900GeV[9] = {	2.164 , 	0.387 , 	0.370 , 	
// 										7.298 , 	0.643 , 	0.588 , 
// 										0.114 , 	0.015 , 	0.014};


Double_t xPtLimits900GeV[14] =  {0.4,0.6,0.8,1.0,1.2,1.4,1.6,2.0,2.5,3.0,3.5,4.0,5.0,7.0} ;
Double_t xPtLimits900GeVGoa[17] =  {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,3.0,3.5,4.0,5.,7.};

//***************************************************************************************************************************
//*********************************** World data Eta to Pi0 *****************************************************************
//***************************************************************************************************************************
TFile* 		fileWorldDataPi0Eta;
TGraphErrors *	graphDonaldson100GeV;
TGraphErrors *	graphDonaldson200GeV;
TGraphErrors *	graphAntille87pp;
TGraphErrors *	graphAntille87ppbar;
TGraphErrors *	graphAguilar400GeV;
TGraphErrors *	graphKourkou79pp;
TGraphErrors *	graphApana530GeV;
TGraphErrors *	graphApana800GeV;
TGraphErrors *	graphKourkou79pp52;
TGraphErrors *	graphAkesson53GeVppbar;
TGraphErrors *	graphAkesson53GeVpp;
TGraphErrors *	graphKourkou79pp62;
TGraphErrors *	graphPhenix200GeV;
TGraphErrors *	graphBanner540GeV;

//***************************************************************************************************************************
//********************************** Color definitions **********************************************************************
//***************************************************************************************************************************

//************************* Variables for styling combined plots inv Yield + Mass *******************************************
Color_t 	colorPi0900GeV 			= kRed;
Color_t 	colorPi02760GeV 			= kMagenta+1;
Color_t	colorPi07TeV				= kBlue;
Color_t 	colorEta7TeV				= kMagenta;
Style_t 	markerStyleConv		 	= 20 ;
Style_t 	markerStyleMCEtaToPi0		= 24 ;
Style_t 	markerStylePHOS		 	= 21 ;
Style_t 	markerStyleEMCAL		 	= 33 ;
Style_t 	markerStyleConvMC			= 24 ;
Style_t 	markerStylePHOSMC		 	= 25 ;
Style_t 	markerStyleEMCALMC		 	= 27 ;
Color_t 	colorConv 				= kBlack;
Color_t 	colorConvMC 				= kGray+1;
Color_t	colorEMCAL				= kGreen+2;
Color_t	colorEMCALMC				= kGreen-6;
Color_t 	colorPHOS					= kBlue+1;
Color_t 	colorPHOSMC				= kBlue-7;
Color_t 	colorPHOSMass					= kRed+1;
Color_t 	colorPHOSMCMass				= kRed-7;
Color_t  colorEMCALMass              = kGreen+2;
Color_t  colorEMCALMCMass            = kGreen-6;
Color_t 	colorConvSyst				= kGray;
Color_t 	colorPHOSSyst				= kBlue-8;
Color_t 	colorEMCALSyst				= kGreen-6;
Style_t 	fillStyleConv			 	= 0 ;
Style_t 	fillStylePHOS			 	= 0;
Style_t 	fillStyleEMCAL			 	= 0;
Color_t 	colorPhojetEtaToPi0			= kRed+2;
Color_t 	colorPhojetEtaToPi0900GeV	= kRed-6;
Color_t 	colorPythiaEtaToPi0			= kBlue+2	;
Color_t 	colorPythiaEtaToPi0900GeV	= kBlue-6	;
Size_t 	markerSizeInvYield			= 1.5;
Size_t 	markerSizeMass				= 3.2;
Size_t 	markerSizeMassMC			= 1.5;

//************************ Charged Spectrum *********************************************************************************
Color_t 	colorCharged 				= kGray;

//*********************** Variables for styling NLO - Calculations **********************************************************
Color_t 	colorNLOPi07TeVMuHalf 		= colorPi07TeV +2;
Color_t 	colorNLOPi07TeVMuOne		= colorPi07TeV;
Color_t 	colorNLOPi07TeVMuTwo		= colorPi07TeV -6 ;
Color_t 	colorNLOBKKPi07TeVMuTwo		= kCyan+3;
Color_t 	colorNLODSSPi07TeVMuTwo		= kBlack;
Color_t 	colorNLOPi0900GeVMuHalf		= colorPi0900GeV +2;
Color_t 	colorNLOPi0900GeVMuOne		= colorPi0900GeV;	
Color_t 	colorNLOPi0900GeVMuTwo		= colorPi0900GeV -6;
Color_t 	colorNLOPi02760GeVMuHalf		= colorPi02760GeV +2;
Color_t 	colorNLOPi02760GeVMuOne		= colorPi02760GeV;	
Color_t 	colorNLOPi02760GeVMuTwo		= colorPi02760GeV -6;
Color_t 	colorNLODSSPi02760GeVMuTwo	= kBlack;
Color_t 	colorNLOBKKPi0900GeVMuTwo	= kOrange;
Color_t 	colorNLODSSPi0900GeVMuTwo	= kBlack;
Color_t 	colorNLOEta7TeVMuHalf		= colorEta7TeV +2;
Color_t 	colorNLOEta7TeVMuOne		= colorEta7TeV;
Color_t 	colorNLOEta7TeVMuTwo		= colorEta7TeV -6;
Color_t 	colorNLOEta900GeVMuHalf		= colorPi0900GeV +2;
Color_t 	colorNLOEta900GeVMuOne		= colorPi0900GeV;
Color_t 	colorNLOEta900GeVMuTwo		= colorPi0900GeV -6;
Color_t 	colorNLOEta2760GeVMuHalf		= colorPi02760GeV +2;
Color_t 	colorNLOEta2760GeVMuOne		= colorPi02760GeV;
Color_t 	colorNLOEta2760GeVMuTwo		= colorPi02760GeV -6;

Style_t 	styleMarkerNLOMuHalf		= 24;
Style_t 	styleMarkerNLOMuOne			= 27;
Style_t 	styleMarkerNLOMuTwo			= 30;
Style_t 	styleLineNLOMuHalf			= 8;
Style_t 	styleLineNLOMuOne			= 7;
Style_t 	styleLineNLOMuTwo			= 4;
Style_t 	styleLineNLOMuTwoBKK		= 3;
Style_t 	styleLineNLOMuTwoDSS		= 6;
Size_t	sizeMarkerNLO				= 1;
Width_t	widthLineNLO				= 2.;

//********************* Variables Combmon Spectrum ************************************************************************
Style_t 	markerStyleCommmonSpectrum 	= 20 ;
Style_t 	markerStyleCommmonSpectrum7TeV 	= 20 ;
Style_t 	markerStyleCommmonSpectrum900GeV = 21 ;
Style_t 	markerStyleCommmonSpectrum2760GeV = 29 ;
Style_t 	markerStyleCommmonSpectrumPi07TeV = 20 ;
Style_t 	markerStyleCommmonSpectrumPi0900GeV 	= 21 ;
Style_t 	markerStyleCommmonSpectrumEta7TeV 	= 29 ;
Size_t 	markerSizeCommonSpectrum 	= 1.;
Size_t 	markerSizeCommonSpectrumPi07TeV 	= 1.8;
Size_t 	markerSizeCommonSpectrumPi0900GeV = 1.8;
Size_t 	markerSizeCommonSpectrumEta7TeV 	= 2.2;
Color_t 	colorCommonSpectrumPi0900GeV 	= colorPi0900GeV+2;
Color_t 	colorCommonSpectrumPi0900GeVBox = colorPi0900GeV-8;
Color_t 	colorCommonSpectrumPi02760GeV = colorPi02760GeV+2;
Color_t 	colorCommonSpectrumPi02760GeVBox = colorPi02760GeV-8;
Color_t 	colorCommonSpectrumPi07TeV 	= colorPi07TeV+3;
Color_t 	colorCommonSpectrumPi07TeVBox = colorPi07TeV-8;
Color_t 	colorCommonSpectrumEta7TeV 	= colorEta7TeV+2;
Color_t 	colorCommonSpectrumEta7TeVBox = colorEta7TeV-8;
Size_t 	markerSizeSpectrum 			= 1.;
Color_t 	colorEMCALPi07TeV			= colorPi07TeV-7 ;
Color_t 	colorEMCALPi0900GeV			= colorPi0900GeV-4;
Color_t 	colorEMCALEta7TeV			= colorEta7TeV-7 ;
Color_t 	colorPHOSPi07TeV			= colorPi07TeV;
Color_t 	colorPHOSPi0900GeV			= colorPi0900GeV;
Color_t 	colorPHOSEta7TeV			= colorEta7TeV;
Color_t 	colorConvPi07TeV			= colorPi07TeV+2;
Color_t 	colorConvPi0900GeV			= colorPi0900GeV+2;
Color_t 	colorConvEta7TeV			= colorEta7TeV+2;
Style_t 	styleFitCommonSpectrum 		= 1;
Width_t 	widthLinesBoxes;
Width_t 	widthStatErrBars;
Width_t 	widthCommonFit;
Width_t 	widthCommonSpectrumBoxes;
Width_t 	widthCommonErrors;
Color_t 	colorFitIndivid			= kRed+2;
Style_t	styleFitIndivid			= 1;

TString 			collisionSystemCombined;
TString 			collisionSystemCombinedReallyAll;

Double_t normalizationInvX3En[5] 		= {0.25,0.32,0.09,0.03,0.023};
Double_t normalizationInvX3EnPi0Eta[5] = {0.2,0.25,0.09,0.03,0.023};
Double_t normalizationInvXOnlySpec[5] 	= {0.65,0.78,0.09,0.03,0.026};
Double_t normalizationInvX1MesonALLEn[5]	= {0.25,0.37,0.09,0.03,0.023};
Double_t normalizationInvXSingleEn[5]	= {0.25,0.41,0.09,0.03,0.028};
Double_t normalizationInvXCTSALLEn[5] 	= {0.18,0.25,0.09,0.03,0.028};
Double_t binshiftingRightMiddle[3]		= {0.75,0.56,0.028};  


//functions for shifting in Y
TF1* fitTsallisPi07TeVPtMult;
TGraphAsymmErrors* graphConversionXSectionPi07TeV;
TGraphAsymmErrors* graphConversionXSectionPi07TeVUnscaled;
TGraphAsymmErrors* graphPHOSXSectionPi07TeV;
TGraphAsymmErrors* graphPHOSXSectionPi07TeVUnscaled;
TGraphAsymmErrors* graphEMCALXSectionPi07TeV;

TF1* fitTsallisEta7TeVPtMult;
TGraphAsymmErrors* graphConversionXSectionEta7TeV;
TGraphAsymmErrors* graphConversionXSectionEta7TeVUnscaled;
TGraphAsymmErrors* graphPHOSXSectionEta7TeV;
TGraphAsymmErrors* graphPHOSXSectionEta7TeVUnscaled;

TF1* fitTsallisPi0900GeVPtMult ;
TGraphAsymmErrors* graphConversionXSectionPi0900GeV;
TGraphAsymmErrors* graphConversionXSectionPi0900GeVUnscaled;
TGraphAsymmErrors* graphPHOSXSectionPi0900GeV;
TGraphAsymmErrors* graphPHOSXSectionPi0900GeVUnscaled;
TGraphAsymmErrors* graphRatioCombConvPi0Sys;
TGraphAsymmErrors* graphRatioCombPHOSEta7TeV;
TGraphAsymmErrors* graphRatioCombPHOSEta7TeVSys;
TGraphAsymmErrors* graphRatioCombConvEta7TeV;
TGraphAsymmErrors* graphRatioCombConvEta7TeVSys;
TGraphAsymmErrors* graphRatioCombConvPi0900GeVSys;

TGraphAsymmErrors* graphConversionXSectionPi02760GeV;
TGraphAsymmErrors* graphConversionXSectionPi02760GeVUnscaled;
TGraphAsymmErrors* graphPHOSXSectionPi02760GeV;
TGraphAsymmErrors* graphPHOSXSectionPi02760GeVUnscaled;
TGraphAsymmErrors* graphConversionXSectionPi02760GeVYShifted;
TGraphAsymmErrors* graphConversionXSectionPi02760GeVYShiftedSys;
TGraphAsymmErrors* graphPhosXSectionPi02760GeVYShifted;
TGraphAsymmErrors* graphPhosXSectionPi02760GeVYShiftedSys;
TGraphAsymmErrors* graphInvCrossSectionPi0Comb2760GeVYShifted;
TGraphAsymmErrors* graphInvCrossSectionPi0Comb2760GeVYShiftedStat;
TGraphAsymmErrors* graphInvCrossSectionPi0Comb2760GeVYShiftedSys;
TH1D* histoPi0PhosSysRAA2760GeV;
TH1D* histoInvCrossSectionPi0SysRAAPHOS2760GeV;
TGraphAsymmErrors* graphSysErrRAAInvCrossSectionPi0PHOS2760GeV;
TGraphAsymmErrors* graphConversionXSectionPi02760GeVYShiftedSysRAA;
TGraphAsymmErrors* graphPhosXSectionPi02760GeVYShiftedSysRAA;

TString fileNameCaloPhosOmega7TeV;
TFile* filePhosOmega7TeV;
TGraphAsymmErrors*	graphOmegaPhos7TeV;
TGraphAsymmErrors* graphOmegaPhosSys7TeV;
TGraphAsymmErrors* graphOmegaPhosComb7TeV;