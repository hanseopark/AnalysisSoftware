/****************************************************************************************************************************
********    provided by Gamma Conversion Group, PWGGA,                                                                  *****
********    Lucia Leardini, leardini@cern.ch                                                                            *****
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
#include "CombineMesonMeasurementsPbPbLHC11hV1.h"


extern TRandom*	gRandom;
extern TBenchmark*	gBenchmark;
extern TSystem*	gSystem;
extern TMinuit*  	gMinuit;	

struct SysErrorConversion {
	Double_t value;
	Double_t error;
	//	TString name;
};

void CombineMesonMeasurementsPbPbLHC11hV1(TString meson = "Eta",
										TString fileNamePCM = "", 
										TString fileNameEMCalFull = "",  
										TString suffix = "pdf", 
										Bool_t PaperPi0 = kTRUE,
										Bool_t noXerrorBars = kTRUE,
										TString isMC= "", 
										TString thesisPlots = "", 
										TString bWCorrection=""){	
	
	gROOT->Reset();	
	gROOT->SetStyle("Plain");
	
	StyleSettingsThesis();	
	SetPlotStyle();
	
	TString date = ReturnDateString();
	TString dateForOutput 			= ReturnDateStringForOutput();
	cout << dateForOutput.Data() << endl;
	
	
    //___________________________________ Labels definition _____________________________________________
	
	TLatex *labelPreliminary = new TLatex(0.62,0.92,"ALICE Preliminary");
	SetStyleTLatex( labelPreliminary, FontSize,4);
	
    TLatex *labelFactorPi0104 = new TLatex(0.85,0.405,"x 4#upoint10^{2}");
    SetStyleTLatex( labelFactorPi0104,FontSize,4,colorCombo0010);
    TLatex *labelFactorPi0100 = new TLatex(0.85,0.34,"x 10^{2}"); 
    SetStyleTLatex( labelFactorPi0100,FontSize,4,colorCombo2050);
    TLatex *labelFactorEta4 = new TLatex(0.85,0.247,"x 4");
    SetStyleTLatex( labelFactorEta4,FontSize,4,colorCombo0010);

    
	TLatex *labelFactorLower;
	if(meson.CompareTo("Pi0")==0){	   labelFactorLower = new TLatex(0.8,0.5,"x 4");}
	else if(meson.CompareTo("Eta")==0){labelFactorLower = new TLatex(0.805,0.35,"x 4");}
	SetStyleTLatex( labelFactorLower, FontSize,4,colorCombo0010);

	TLatex *labelFactorLowerOnlyPbPb;
	if(meson.CompareTo("Pi0")==0){	   labelFactorLowerOnlyPbPb = new TLatex(0.8,0.37,"x 4");} 
	else if(meson.CompareTo("Eta")==0){labelFactorLowerOnlyPbPb = new TLatex(0.805,0.35,"x 4");}
	SetStyleTLatex( labelFactorLowerOnlyPbPb,FontSize,4,colorCombo0010);

	TLatex *labelFactorUpper;
	if(meson.CompareTo("Pi0")==0){	   labelFactorUpper = new TLatex(0.415,0.89,"x 4");} 
	else if(meson.CompareTo("Eta")==0){labelFactorUpper = new TLatex(0.43,0.835,"x 4");}
	SetStyleTLatex( labelFactorUpper,FontSize,4,colorCombo0010);

	TLatex *labelSystOnlyPbPb;
	if(meson.CompareTo("Pi0")==0){	   labelSystOnlyPbPb= new TLatex(0.75,0.93,"#pi^{0} #rightarrow #gamma#gamma");} 
	else if(meson.CompareTo("Eta")==0){labelSystOnlyPbPb= new TLatex(0.75,0.93,"#eta #rightarrow #gamma#gamma");}
	SetStyleTLatex(labelSystOnlyPbPb,FontSize,4);
	
	TLatex *labelSyst;
	if(meson.CompareTo("Pi0")==0){	   labelSyst= new TLatex(0.6,0.93,"#pi^{0} #rightarrow #gamma#gamma");} 
	else if(meson.CompareTo("Eta")==0){labelSyst= new TLatex(0.6,0.93,"#eta #rightarrow #gamma#gamma");}
	SetStyleTLatex(labelSyst,FontSize,4);
	
	TLatex *labelSystRaa;
	if(meson.CompareTo("Pi0")==0){	   labelSystRaa= new TLatex(0.6,0.89,"#pi^{0} #rightarrow #gamma#gamma");} 
	else if(meson.CompareTo("Eta")==0){labelSystRaa= new TLatex(0.6,0.89,"#eta #rightarrow #gamma#gamma");}
	SetStyleTLatex( labelSystRaa,FontSize,4);
	
	TLatex *labelEnergyInvYieldSectionPi0LHC11h = new TLatex(0.6,0.88,collisionSystem2760GeV.Data());
	SetStyleTLatex( labelEnergyInvYieldSectionPi0LHC11h,FontSize,4);
	TLatex *labelDetSysInvYieldSectionPi0LHC11h;
	if(meson.CompareTo("Pi0")==0){
		labelDetSysInvYieldSectionPi0LHC11h= new TLatex(0.6,0.84,"#pi^{0} #rightarrow #gamma#gamma");
	} else if(meson.CompareTo("Eta")==0){
		labelDetSysInvYieldSectionPi0LHC11h= new TLatex(0.6,0.84,"#eta #rightarrow #gamma#gamma");
	}
	SetStyleTLatex( labelDetSysInvYieldSectionPi0LHC11h, FontSize,4);

	TLatex *labelEnergyInvYieldSectionPi0LHC11hnoPrelim = new TLatex(0.6,0.92,collisionSystem2760GeV.Data());
	SetStyleTLatex( labelEnergyInvYieldSectionPi0LHC11hnoPrelim, 0.035,4);
	TLatex *labelDetSysInvYieldSectionPi0LHC11hnoPrelim;
	if(meson.CompareTo("Pi0")==0){
		labelDetSysInvYieldSectionPi0LHC11hnoPrelim= new TLatex(0.6,0.88,"#pi^{0} #rightarrow #gamma#gamma");
	} else if(meson.CompareTo("Eta")==0){
		labelDetSysInvYieldSectionPi0LHC11hnoPrelim= new TLatex(0.6,0.88,"#eta #rightarrow #gamma#gamma");
	}
	SetStyleTLatex( labelDetSysInvYieldSectionPi0LHC11hnoPrelim, 0.035,4);

	TLatex *labelDetSysInvYieldSectionPi0LHC11hwithPP;
	if(meson.CompareTo("Pi0")==0){
		labelDetSysInvYieldSectionPi0LHC11hwithPP= new TLatex(0.62,0.88,"#pi^{0} #rightarrow #gamma#gamma");
	} else if(meson.CompareTo("Eta")==0){
		labelDetSysInvYieldSectionPi0LHC11hwithPP= new TLatex(0.62,0.88,"#eta #rightarrow #gamma#gamma");
	}
	SetStyleTLatex( labelDetSysInvYieldSectionPi0LHC11hwithPP, 0.035,4);


	for (Int_t i = 0; i < 11; i++){
		colorDet[i]					= GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kFALSE, kFALSE, kTRUE);
		colorDetMC[i]				= GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kTRUE, kFALSE, kTRUE);
		markerStyleDet[i]			= GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kFALSE);
		markerStyleDetMC[i]			= GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kTRUE);
		markerSizeDet[i]			= GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kFALSE)*2;
		markerSizeDetMC[i]			= GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kTRUE)*2;
	}	
	
    Double_t mesonMassExpectPi0 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    Double_t mesonMassExpectEta = TDatabasePDG::Instance()->GetParticle(221)->Mass();
	
	Double_t normErr0010 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr0010/tAA0010),2)+pow((commonCentralityErr0010/100),2),0.5);
	Double_t normErr2040 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr2040/tAA2040),2)+pow((commonCentralityErr2040/100),2),0.5);
	Double_t normErr2050 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr2050/tAA2050),2)+pow((commonCentralityErr2040/100),2),0.5);

	TBox* boxErrorNorm0010_Single = CreateBoxConv(colorComb0005Box, 0.2, 1.-normErr0010 , 0.5, 1.+normErr0010);
	TBox* boxErrorNorm2040_Single = CreateBoxConv(colorComb2040Box, 0.2, 1.-normErr2040 , 0.5, 1.+normErr2040);
	TBox* boxErrorNorm2050_Single = CreateBoxConv(colorCombo2050, 0.2, 1.-normErr2050 , 0.5, 1.+normErr2050);

	TBox* boxErrorNorm0010 = CreateBoxConv(kRed-7, 0.25, 1.-normErr0010 , 0.5, 1.+normErr0010);
	TBox* boxErrorNorm2050 = CreateBoxConv(kAzure-4, 0.5, 1.-normErr2050 , 0.75, 1.+normErr2050);

	TBox* boxErrorNorm0010Only = CreateBoxConv(kRed-7, 0.25, 1.-normErr0010 , 0.5, 1.+normErr0010);
	TBox* boxErrorNorm2050Only = CreateBoxConv(kAzure-4, 0.25, 1.-normErr2050 , 0.5, 1.+normErr2050);

	
	//___________________________________ Declaration of files _____________________________________________
	
	TString fileNameTheory					= "ExternalInputPbPb/Theory/TheoryCompilationPbPbforLHC11h.root";	

	TString filePPpublished					= "CombinedResultsPaperPP2760GeV_2016_01_08.root";//"CombinedResultsPaperX_18_Feb_2014.root";
    TString fileFinalResultsPi0PPforRAA           = "CombinedResultsPaperX_18_Feb_2014.root";
    TString fileFinalResultsEtaPPforRAA           = "CombinedResultsEta2760X_2015_09_01.root";    // 	TString filePPEtaPCM					= "LHC11hExternalInputs/CombinedResultsEta2760X_2015_09_01.root"; 
	TString fileNameChargedRatios			= "ExternalInputPbPb/IdentifiedCharged/JIRA_PWGLF-258/PbPb276.fullpT.RATIOS.20140329.root";
	TString fileNameChargedPionRAA			= "ExternalInputPbPb/IdentifiedCharged/JIRA_PWGLF-258/RAA_Pion_08052014.root";
	TString fileNameChargedKaonRAA			= "ExternalInputPbPb/IdentifiedCharged/JIRA_PWGLF-258/RAA_Kaon_08052014.root";
	TString	fileNameDataOtherEnergyInput 	= "ExternalInputPbPb/OtherExperiments/DataCompilationFromOtherEnergiesPbPbWithEta.root";
	
	TString outputDir 						= Form("%s/%s/CombineMesonMeasurements2760GeVV3%s",suffix.Data(),dateForOutput.Data(),bWCorrection.Data());
	TString nameFinalResDat 				= Form("%s/CombinedResults%s_FitResults.dat",outputDir.Data(),bWCorrection.Data());
	cout << outputDir.Data() << endl;
	cout << fileNamePCM.Data() << endl;	
	cout << fileNameEMCalFull.Data() << endl;

	gSystem->Exec("mkdir -p "+outputDir);
 	gSystem->Exec(Form("cp %s %s/InputPCM.root", fileNamePCM.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp LHC11hExternalInputs/%s %s/InputEMCalFull.root", fileNameEMCalFull.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/Theory.root", fileNameTheory.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/ChargedRatios.root", fileNameChargedRatios.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/ChargedPionRAA.root", fileNameChargedPionRAA.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/ChargedKaonRAA.root", fileNameChargedKaonRAA.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/OtherExperimentsInput.root", fileNameDataOtherEnergyInput.Data(), outputDir.Data()));
	
	TString paperPlots = Form("%s/PaperPlotsforPreview_%s",outputDir.Data(),dateForOutput.Data());
	gSystem->Exec("mkdir -p "+paperPlots);
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

	//*********************************************************************************************************//
	//*************************************        Theory inputs      *****************************************//
	//*********************************************************************************************************//
	TFile* fileTheoryGraphs = new TFile(fileNameTheory);
		
	//**************************** pQCD calculation ****************************//
    //****************************     Pi0 DSS14    ****************************//
    TGraphAsymmErrors *graphNLOCalcmuHalfDSS14InvSecPi02760GeV = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcmuHalfDSS14InvSecPi02760GeV");
	TGraphAsymmErrors *graphNLOCalcmuHalfDSS14InvYieldPi02760GeV = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcmuHalfDSS14InvSecPi02760GeV");
	TGraphAsymmErrors *graphNLOCalcmuHalfDSS14InvSecPi02760GeVNcoll0010 = (TGraphAsymmErrors*)graphNLOCalcmuHalfDSS14InvSecPi02760GeV->Clone("graphNLOCalcmuHalfDSS14InvSecPi02760GeVNcoll0010");
      graphNLOCalcmuHalfDSS14InvSecPi02760GeVNcoll0010->RemovePoint(0);
      graphNLOCalcmuHalfDSS14InvSecPi02760GeVNcoll0010 = ScaleGraph(graphNLOCalcmuHalfDSS14InvSecPi02760GeVNcoll0010,nColl0010*4/xSection2760GeVppINEL);
	TGraphAsymmErrors* graphNLOCalcmuHalfDSS14InvSecPi02760GeVNcoll2050 = (TGraphAsymmErrors*)graphNLOCalcmuHalfDSS14InvSecPi02760GeV->Clone("graphNLOCalcmuHalfDSS14InvSecPi02760GeVNcoll2050");
      graphNLOCalcmuHalfDSS14InvSecPi02760GeVNcoll2050->RemovePoint(0);
      graphNLOCalcmuHalfDSS14InvSecPi02760GeVNcoll2050 = ScaleGraph(graphNLOCalcmuHalfDSS14InvSecPi02760GeVNcoll2050,nColl2050/xSection2760GeVppINEL);	
      
	TGraphAsymmErrors *graphNLOCalcmuTwoDSS14InvSecPi02760GeV = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcmuTwoDSS14InvSecPi02760GeV");
	TGraphAsymmErrors *graphNLOCalcmuTwoDSS14InvYieldPi02760GeV = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcmuTwoDSS14InvSecPi02760GeV");
	TGraphAsymmErrors *graphNLOCalcmuTwoDSS14InvSecPi02760GeVNcoll0010 = (TGraphAsymmErrors*)graphNLOCalcmuTwoDSS14InvSecPi02760GeV->Clone("graphNLOCalcmuTwoDSS14InvSecPi02760GeVNcoll0010");
      graphNLOCalcmuTwoDSS14InvSecPi02760GeVNcoll0010->RemovePoint(0);
      graphNLOCalcmuTwoDSS14InvSecPi02760GeVNcoll0010 = ScaleGraph(graphNLOCalcmuTwoDSS14InvSecPi02760GeVNcoll0010,nColl0010*4/xSection2760GeVppINEL);
	TGraphAsymmErrors* graphNLOCalcmuTwoDSS14InvSecPi02760GeVNcoll2050 = (TGraphAsymmErrors*)graphNLOCalcmuTwoDSS14InvSecPi02760GeV->Clone("graphNLOCalcmuTwoDSS14InvSecPi02760GeVNcoll2050");
      graphNLOCalcmuTwoDSS14InvSecPi02760GeVNcoll2050->RemovePoint(0);
      graphNLOCalcmuTwoDSS14InvSecPi02760GeVNcoll2050 = ScaleGraph(graphNLOCalcmuTwoDSS14InvSecPi02760GeVNcoll2050,nColl2050/xSection2760GeVppINEL);
	
	TGraphAsymmErrors *graphNLOCalcDSS14InvSecPi02760GeV = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcDSS14InvCrossSecPi02760GeV");
	TGraphAsymmErrors *graphNLOCalcDSS14InvYieldPi02760GeV = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcDSS14InvYieldPi02760GeV");
	TGraphAsymmErrors *graphNLOCalcDSS14InvYieldPi02760GeVScaledNcoll0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcDSS14InvYieldPi02760GeV");
      graphNLOCalcDSS14InvYieldPi02760GeVScaledNcoll0010->RemovePoint(0);
      graphNLOCalcDSS14InvYieldPi02760GeVScaledNcoll0010 = ScaleGraph(graphNLOCalcDSS14InvYieldPi02760GeVScaledNcoll0010,nColl0010*4);
	TGraphAsymmErrors *graphNLOCalcDSS14InvYieldPi02760GeVScaledNcoll2050 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcDSS14InvYieldPi02760GeV");
      graphNLOCalcDSS14InvYieldPi02760GeVScaledNcoll2050->RemovePoint(0);
      graphNLOCalcDSS14InvYieldPi02760GeVScaledNcoll2050 = ScaleGraph(graphNLOCalcDSS14InvYieldPi02760GeVScaledNcoll2050,nColl2050);

	TGraphAsymmErrors* graphNLOCalcDSS14InvYieldPi0Band = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcDSS14InvCrossSecPi02760GeV");
	Double_t *yPi0PointDSS14muHalf = graphNLOCalcmuHalfDSS14InvSecPi02760GeV->GetY();
	Double_t *yPi0PointDSS14muOne = graphNLOCalcDSS14InvSecPi02760GeV->GetY();
	Double_t *yPi0PointDSS14muTwo  = graphNLOCalcmuTwoDSS14InvSecPi02760GeV->GetY();
	for(Int_t i =0; i<graphNLOCalcDSS14InvYieldPi0Band->GetN(); i++){	
		graphNLOCalcDSS14InvYieldPi0Band->SetPointEYhigh(i,yPi0PointDSS14muTwo[i]-yPi0PointDSS14muOne[i]);
		graphNLOCalcDSS14InvYieldPi0Band->SetPointEYlow(i,yPi0PointDSS14muOne[i]-yPi0PointDSS14muHalf[i]);
	}
	graphNLOCalcDSS14InvYieldPi0Band->RemovePoint(0);
	graphNLOCalcDSS14InvYieldPi0Band = ScaleGraph(graphNLOCalcDSS14InvYieldPi0Band,1./xSection2760GeVppINEL);
	
	
    //****************************     Eta DSS07    ****************************//
	TGraphAsymmErrors* graphNLOCalcmuHalfDSS07InvSecEta2760GeV = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcmuHalfDSS07InvSecEta2760GeV");
	TGraphAsymmErrors* graphNLOCalcmuHalfDSS07InvYieldEta2760GeV = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcmuHalfDSS07InvSecEta2760GeV");
	TGraphAsymmErrors* graphNLOCalcmuHalfDSS07InvSecEta2760GeVNcoll0010 = (TGraphAsymmErrors*)graphNLOCalcmuHalfDSS07InvSecEta2760GeV->Clone("graphNLOCalcmuHalfDSS07InvSecEta2760GeVNcoll0010");
	graphNLOCalcmuHalfDSS07InvSecEta2760GeVNcoll0010->RemovePoint(0);
	graphNLOCalcmuHalfDSS07InvSecEta2760GeVNcoll0010 = ScaleGraph(graphNLOCalcmuHalfDSS07InvSecEta2760GeVNcoll0010,nColl0010*4/xSection2760GeVppINEL);
	TGraphAsymmErrors* graphNLOCalcmuHalfDSS07InvSecEta2760GeVNcoll2050 = (TGraphAsymmErrors*)graphNLOCalcmuHalfDSS07InvSecEta2760GeV->Clone("graphNLOCalcmuHalfDSS07InvSecEta2760GeVNcoll2050");
	graphNLOCalcmuHalfDSS07InvSecEta2760GeVNcoll2050->RemovePoint(0);
	graphNLOCalcmuHalfDSS07InvSecEta2760GeVNcoll2050 = ScaleGraph(graphNLOCalcmuHalfDSS07InvSecEta2760GeVNcoll2050,nColl2050/xSection2760GeVppINEL);
	
	TGraphAsymmErrors* graphNLOCalcmuTwoDSS07InvSecEta2760GeV = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcmuTwoDSS07InvSecEta2760GeV");
	TGraphAsymmErrors* graphNLOCalcmuTwoDSS07InvYieldEta2760GeV = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcmuTwoDSS07InvSecEta2760GeV");
	TGraphAsymmErrors* graphNLOCalcmuTwoDSS07InvSecEta2760GeVNcoll0010 = (TGraphAsymmErrors*)graphNLOCalcmuTwoDSS07InvSecEta2760GeV->Clone("graphNLOCalcmuTwoDSS07InvSecEta2760GeVNcoll0010");
	graphNLOCalcmuTwoDSS07InvSecEta2760GeVNcoll0010->RemovePoint(0);
	graphNLOCalcmuTwoDSS07InvSecEta2760GeVNcoll0010 = ScaleGraph(graphNLOCalcmuTwoDSS07InvSecEta2760GeVNcoll0010,nColl0010*4/xSection2760GeVppINEL);
	TGraphAsymmErrors* graphNLOCalcmuTwoDSS07InvSecEta2760GeVNcoll2050 = (TGraphAsymmErrors*)graphNLOCalcmuTwoDSS07InvSecEta2760GeV->Clone("graphNLOCalcmuTwoDSS07InvSecEta2760GeVNcoll2050");
	graphNLOCalcmuTwoDSS07InvSecEta2760GeVNcoll2050->RemovePoint(0);
	graphNLOCalcmuTwoDSS07InvSecEta2760GeVNcoll2050 = ScaleGraph(graphNLOCalcmuTwoDSS07InvSecEta2760GeVNcoll2050,nColl2050/xSection2760GeVppINEL);
	
	TGraphAsymmErrors* graphNLOCalcDSS07InvSecEta2760GeV = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcDSS07InvCrossSecEta2760GeV");
	TGraphAsymmErrors* graphNLOCalcDSS07InvYieldEta2760GeV = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcDSS07InvYieldEta2760GeV");	
	TGraphAsymmErrors* graphNLOCalcDSS07InvYieldEta2760GeVScaledNcoll0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcDSS07InvYieldEta2760GeV");
	graphNLOCalcDSS07InvYieldEta2760GeVScaledNcoll0010 = ScaleGraph(graphNLOCalcDSS07InvYieldEta2760GeVScaledNcoll0010,nColl0010*4);
	graphNLOCalcDSS07InvYieldEta2760GeVScaledNcoll0010->RemovePoint(0);
	TGraphAsymmErrors* graphNLOCalcDSS07InvYieldEta2760GeVScaledNcoll2050 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcDSS07InvYieldEta2760GeV");
	graphNLOCalcDSS07InvYieldEta2760GeVScaledNcoll2050 = ScaleGraph(graphNLOCalcDSS07InvYieldEta2760GeVScaledNcoll2050,nColl2050);
	graphNLOCalcDSS07InvYieldEta2760GeVScaledNcoll2050->RemovePoint(0);
	
	TGraphAsymmErrors* graphNLOCalcDSS07InvYieldEtaBand = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphNLOCalcDSS07InvCrossSecEta2760GeV");
	Double_t *yEtaPointDSS07muHalf = graphNLOCalcmuHalfDSS07InvSecEta2760GeV->GetY();
	Double_t *yEtaPointDSS07muOne = graphNLOCalcDSS07InvSecEta2760GeV->GetY();
	Double_t *yEtaPointDSS07muTwo  = graphNLOCalcmuTwoDSS07InvSecEta2760GeV->GetY();
	for(Int_t i =0; i<graphNLOCalcDSS07InvYieldEtaBand->GetN(); i++){	
		graphNLOCalcDSS07InvYieldEtaBand->SetPointEYhigh(i,yEtaPointDSS07muTwo[i]-yEtaPointDSS07muOne[i]);
		graphNLOCalcDSS07InvYieldEtaBand->SetPointEYlow(i,yEtaPointDSS07muOne[i]-yEtaPointDSS07muHalf[i]);
	}
	graphNLOCalcDSS07InvYieldEtaBand->RemovePoint(0);
	graphNLOCalcDSS07InvYieldEtaBand = ScaleGraph(graphNLOCalcDSS07InvYieldEtaBand,1./xSection2760GeVppINEL);
	
	

    //**************************** EPOS ****************************//
    TFile* filePredictionEpos = new TFile("ExternalInputPbPb/Theory/EPOS/pi0_eta_forPbPb11h_NoUrQMD.root");

    TGraphErrors *graphEPOSPi0_0010 = (TGraphErrors*)fileTheoryGraphs->Get("epos_pi0_pt_cent0_10");
    TGraphErrors *graphEPOSEta_0010 = (TGraphErrors*)fileTheoryGraphs->Get("epos_eta_pt_cent0_10");
    
    TH1D *histoEPOSPi0_0010 = (TH1D*)GraphToHist(graphEPOSPi0_0010,14,"histoEPOSPi0_0010");
    TH1D *histoEPOSEta_0010 = (TH1D*)GraphToHist(graphEPOSEta_0010,14,"histoEPOSEta_0010");
    TH1D *histoEPOSPi0_2050 = (TH1D*)filePredictionEpos->Get("rb_pi0_pt_cent20_50");
    TH1D *histoEPOSEta_2050 = (TH1D*)filePredictionEpos->Get("rb_eta_pt_cent20_50");

    
    TH1D *histoEPOSpred_0010;
    TH1D *histoEPOSpred_2050;
    if(meson.CompareTo("Pi0")==0){
        histoEPOSpred_0010 = (TH1D*)histoEPOSPi0_0010->Clone("histoEPOSpred_0010");
        histoEPOSpred_2050 = (TH1D*)histoEPOSPi0_2050->Clone("histoEPOSpred_2050");
                    
    } else if(meson.CompareTo("Eta")==0){
        histoEPOSpred_0010 = (TH1D*)histoEPOSEta_0010->Clone("histoEPOSpred_0010");
        histoEPOSpred_2050 = (TH1D*)histoEPOSEta_2050->Clone("histoEPOSpred_2050");
        
    }

//from here down the rebin for the 0-10 only
    Int_t nbins = histoEPOSpred_0010->GetNbinsX();
    Double_t xRebin[nbins];
    for(Int_t b = 0; b<nbins; b++){
        if(histoEPOSpred_0010->GetBinCenter(b+1) < 4.){
            xRebin[b] = histoEPOSpred_0010->GetBinLowEdge(b+1);
        } else {
            xRebin[b] = xRebin[b-1] + 1.;
            cout << "new low edge: " << xRebin[b] << endl;
        }
    }
    Int_t rebinFactor = histoEPOSpred_0010->FindBin(xRebin[0]+1.) - histoEPOSpred_0010->FindBin(xRebin[0]);
    cout << "rebin factor: " << rebinFactor << endl;
    TH1D* histoEPOSRebin_0010 = (TH1D*)histoEPOSpred_0010->Rebin(nbins-1,"rebinEPOSPi0_0010",xRebin);
    TH1D* DeltaPt  = new TH1D("DeltaPt","",nbins-1,xRebin);
    for(Int_t b = 0; b<nbins; b++){
        if(histoEPOSRebin_0010->GetBinCenter(b+1) < 4.){
            DeltaPt->SetBinContent(b+1,1);
        } else {
            DeltaPt->SetBinContent(b+1,rebinFactor);
        }
    }
    DeltaPt->Sumw2();
    histoEPOSRebin_0010->Divide(DeltaPt);
// until here
    
    for (Int_t i = 1; i < histoEPOSRebin_0010->GetNbinsX()+1 ; i++){
        Double_t newBinContent = histoEPOSRebin_0010->GetBinContent(i)/(2*TMath::Pi()*histoEPOSRebin_0010->GetBinCenter(i));
        Double_t newBinError = histoEPOSRebin_0010->GetBinError(i)/(2*TMath::Pi()*histoEPOSRebin_0010->GetBinCenter(i));
        histoEPOSRebin_0010->SetBinContent(i,newBinContent);
        histoEPOSRebin_0010->SetBinError(i,newBinError);
    }
    TH1D* histoEPOSRebin_2050 = (TH1D*)histoEPOSpred_2050->Clone("histoEPOSpred_2050");
    for (Int_t i = 1; i < histoEPOSRebin_2050->GetNbinsX()+1 ; i++){
        Double_t newBinContent = histoEPOSRebin_2050->GetBinContent(i)/(2*TMath::Pi()*histoEPOSRebin_2050->GetBinCenter(i));
        Double_t newBinError = histoEPOSRebin_2050->GetBinError(i)/(2*TMath::Pi()*histoEPOSRebin_2050->GetBinCenter(i));
        histoEPOSRebin_2050->SetBinContent(i,newBinContent);
        histoEPOSRebin_2050->SetBinError(i,newBinError);
    }

    TH1D* histoEtaforRatioEPOS_0010 = (TH1D*)histoEPOSEta_0010->Clone("histoEtaforRatioEPOS_0010");
    TH1D* histoEtaforRatioEPOS_2050 = (TH1D*)histoEPOSEta_2050->Clone("histoEtaforRatioEPOS_2050");
    TH1D* histoPi0forRatioEPOS_0010 = (TH1D*)histoEPOSPi0_0010->Clone("histoPi0forRatioEPOS_0010");
    TH1D* histoPi0forRatioEPOS_2050 = (TH1D*)histoEPOSPi0_2050->Clone("histoPi0forRatioEPOS_2050");
    
    TH1D *histoEtaToPi0EPOSRebin_0010 = (TH1D*)histoEtaforRatioEPOS_0010->Clone("histoEtaToPi0EPOSRebin_0010");
      histoEtaToPi0EPOSRebin_0010->Divide(histoEtaToPi0EPOSRebin_0010,histoPi0forRatioEPOS_0010,1.,1.,"");
    TH1D *histoEtaToPi0EPOSRebin_2050 = (TH1D*)histoEtaforRatioEPOS_2050->Clone("histoEtaToPi0EPOSRebin_2050");
      histoEtaToPi0EPOSRebin_2050->Divide(histoEtaToPi0EPOSRebin_2050,histoPi0forRatioEPOS_2050,1.,1.,"");
    
	
	//**************************** Cracow model ****************************//	
	TGraphAsymmErrors *TheoryCracowPi0LowPt_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("TheoryCracowPi0LowPt_0010");
	TGraphAsymmErrors *TheoryCracowPi0LowPt_2050 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("TheoryCracowPi0LowPt_2050");

	TGraphAsymmErrors *TheoryCracowEtaLowPt_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("TheoryCracowEtaLowPt_0010");
	TGraphAsymmErrors *TheoryCracowEtaLowPt_2050 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("TheoryCracowEtaLowPt_2050");

	TGraphAsymmErrors* TheoryCracowEtaToPi0LowPt_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("TheoryCracowEtaToPi0LowPt_0010");
	TGraphAsymmErrors* TheoryCracowEtaToPi0LowPt_2050 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("TheoryCracowEtaToPi0LowPt_2050");

    TGraphAsymmErrors* TheoryCracowLowPt_0010 = NULL;
	TGraphAsymmErrors* TheoryCracowLowPt_2050 = NULL;
	if(meson.CompareTo("Pi0")==0){		
		TheoryCracowLowPt_0010 = (TGraphAsymmErrors*)TheoryCracowPi0LowPt_0010->Clone("TheoryCracowLowPt_0010"); 
		TheoryCracowLowPt_2050 = (TGraphAsymmErrors*)TheoryCracowPi0LowPt_2050->Clone("TheoryCracowLowPt_2050");
	} else if(meson.CompareTo("Eta")==0){
		TheoryCracowLowPt_0010 = (TGraphAsymmErrors*)TheoryCracowEtaLowPt_0010->Clone("TheoryCracowLowPt_0010"); 
		TheoryCracowLowPt_2050 = (TGraphAsymmErrors*)TheoryCracowEtaLowPt_2050->Clone("TheoryCracowLowPt_2050"); 
	}
	TheoryCracowLowPt_0010 = ScaleGraph(TheoryCracowLowPt_0010,4);
	
	
	//**************************** Jet Quenching NLO arXiv:1506.00838 ****************************//
	graphPi0JetQuenching18 = (TGraph*)fileTheoryGraphs->Get("graphPi0JetQuenching18_0010");
	graphPi0JetQuenching22 = (TGraph*)fileTheoryGraphs->Get("graphPi0JetQuenching22_0010");
	graphPi0JetQuenching26 = (TGraph*)fileTheoryGraphs->Get("graphPi0JetQuenching26_0010");
	
	graphEtaJetQuenching18 = (TGraph*)fileTheoryGraphs->Get("graphEtaJetQuenching18_0010");
	graphEtaJetQuenching22 = (TGraph*)fileTheoryGraphs->Get("graphEtaJetQuenching22_0010");
	graphEtaJetQuenching26 = (TGraph*)fileTheoryGraphs->Get("graphEtaJetQuenching26_0010");
	
	graphEtatoPi0RatioJetQuenching18 = (TGraph*)fileTheoryGraphs->Get("graphEtatoPi0RatioJetQuenching18_0010");
	graphEtatoPi0RatioJetQuenching22 = (TGraph*)fileTheoryGraphs->Get("graphEtatoPi0RatioJetQuenching22_0010");
	graphEtatoPi0RatioJetQuenching26 = (TGraph*)fileTheoryGraphs->Get("graphEtatoPi0RatioJetQuenching26_0010");
	
	
	TGraphAsymmErrors* graphPi0RAAJetQuenching_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphPi0RAAJetQuenching_0010");
	TGraphAsymmErrors* graphEtaRAAJetQuenching_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphEtaRAAJetQuenching_0010");
	TGraphAsymmErrors* graphEtaToPi0JetQuenching_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphEtaToPi0JetQuenching_0010");
	
	
	//**************************** WHDG ****************************//
	gWHDG_Raa_0510 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGRAA0510");
	gWHDG_Raa_0005 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGRAA0005");
	gWHDG_Raa_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGRAA0010");
	gWHDG_Raa_2050 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGRAA2050");

	gWHDG_Eta_Raa_0010 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGetaRAA0010");
	gWHDG_Eta_Raa_2050 = (TGraphAsymmErrors*)fileTheoryGraphs->Get("graphWHDGetaRAA2050");
    
	
	
	//*********************************************************************************************************//
	//********************************        Other experiment inputs     *************************************//
	//*********************************************************************************************************//
	
	TFile* fileDataOtherEnergies = new TFile(fileNameDataOtherEnergyInput);

	graphWA98_17_3GeVRAA_0013= (TGraphErrors*)fileDataOtherEnergies->Get("graphWA98RAA_0013");

	graphPHENIX200GeVRAA_0010= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVRAA_0010");
	graphPHENIX200GeVRAA_2040= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVRAA_2040");
	graphPHENIX39GeVRAA_0010= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX39GeVRAA_0010");
	graphPHENIX39GeVRAA_2040= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX39GeVRAA_2040");
	graphPHENIX62GeVRAA_0010= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX62GeVRAA_0010");
	graphPHENIX62GeVRAA_2040= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX62GeVRAA_2040");

	graphPHENIX200GeVEtaRAA_0010= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVEtaRAA_0010");
	graphPHENIX200GeVEtaRAA_2040= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVEtaRAA_2040");
	graphPHENIX200GeVEtaRAA_2060= (TGraphErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVEtaRAA_2060");

	TGraphAsymmErrors* graphEtaRAAPhenix0010 = (TGraphAsymmErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVEtaRAA_0010"); 
	TGraphAsymmErrors* graphEtaRAAPhenix2040 = (TGraphAsymmErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVEtaRAA_2040"); 
	TGraphAsymmErrors* graphEtaRAAPhenix2060 = (TGraphAsymmErrors*)fileDataOtherEnergies->Get("graphPHENIX200GeVEtaRAA_2060"); 
	
	
	
	//*********************************************************************************************************//
	//**********************************     Charged pions and kaons     **************************************//
	//*********************************************************************************************************//
	
	//**************************** Charged ratios ****************************//
	TFile* fileDataALICEChargedRatioKaonToPion	= new TFile(fileNameChargedRatios);	
	TH1D *histoStatChargedRatioKaonToPion0005 = (TH1D*)fileDataALICEChargedRatioKaonToPion->Get("hstat_PbPb276_0005_kaon_to_pion_sum");
	TGraphAsymmErrors* graphChargedRatioKaonToPion0005 = new TGraphAsymmErrors(histoStatChargedRatioKaonToPion0005); 
	TH1D *histoSystChargedRatioKaonToPion0005 = (TH1D*)fileDataALICEChargedRatioKaonToPion->Get("hsys_PbPb276_0005_kaon_to_pion_sum");
	TGraphAsymmErrors* graphChargedRatioKaonToPionSys0005 = new TGraphAsymmErrors(histoSystChargedRatioKaonToPion0005); 
	
	TH1D *histoStatChargedRatioKaonToPion0510 = (TH1D*)fileDataALICEChargedRatioKaonToPion->Get("hstat_PbPb276_0510_kaon_to_pion_sum");
	TGraphAsymmErrors* graphChargedRatioKaonToPion0510 = new TGraphAsymmErrors(histoStatChargedRatioKaonToPion0510); 
	TH1D *histoSystChargedRatioKaonToPion0510 = (TH1D*)fileDataALICEChargedRatioKaonToPion->Get("hsys_PbPb276_0510_kaon_to_pion_sum");
	TGraphAsymmErrors* graphChargedRatioKaonToPionSys0510 = new TGraphAsymmErrors(histoSystChargedRatioKaonToPion0510); 

	TH1D* histoStatChargedRatioKaonToPion0010 = (TH1D*)histoStatChargedRatioKaonToPion0510->Clone("histoStatChargedRatioKaonToPion0010");
	TH1D* histoSystChargedRatioKaonToPion0010 = (TH1D*)histoSystChargedRatioKaonToPion0510->Clone("histoSystChargedRatioKaonToPion0010");
	histoStatChargedRatioKaonToPion0010->Add(histoStatChargedRatioKaonToPion0005);
	histoSystChargedRatioKaonToPion0010->Add(histoSystChargedRatioKaonToPion0005);
	histoStatChargedRatioKaonToPion0010->Scale(0.5);
	histoSystChargedRatioKaonToPion0010->Scale(0.5);
	for (Int_t i = 1; i < histoSystChargedRatioKaonToPion0010->GetNbinsX()+1; i++){
		Double_t relErrLowerCent = 0;
		if (histoSystChargedRatioKaonToPion0005->GetBinContent(i) != 0){
			relErrLowerCent= histoSystChargedRatioKaonToPion0005->GetBinError(i)/histoSystChargedRatioKaonToPion0005->GetBinContent(i)*100 ;
		}
		Double_t relErrHigherCent = 0;
		if (histoSystChargedRatioKaonToPion0510->GetBinContent(i) != 0){
			relErrHigherCent = histoSystChargedRatioKaonToPion0510->GetBinError(i)/histoSystChargedRatioKaonToPion0510->GetBinContent(i)*100 ;
		}
		
		if (relErrHigherCent > relErrLowerCent){
			histoSystChargedRatioKaonToPion0010->SetBinError(i, histoSystChargedRatioKaonToPion0010->GetBinContent(i)*relErrHigherCent/100);
		} else {
			histoSystChargedRatioKaonToPion0010->SetBinError(i, histoSystChargedRatioKaonToPion0010->GetBinContent(i)*relErrLowerCent/100);
		}         
	}   
	TGraphAsymmErrors* graphChargedRatioKaonToPion0010 = new TGraphAsymmErrors(histoStatChargedRatioKaonToPion0010); 
	TGraphAsymmErrors* graphChargedRatioKaonToPionSys0010 = new TGraphAsymmErrors(histoSystChargedRatioKaonToPion0010); 
	
	TH1D *histoStatChargedRatioKaonToPion2040 = (TH1D*)fileDataALICEChargedRatioKaonToPion->Get("hstat_PbPb276_2040_kaon_to_pion_sum");
	TGraphAsymmErrors* graphChargedRatioKaonToPion2040 = new TGraphAsymmErrors(histoStatChargedRatioKaonToPion2040); 
	TH1D *histoSystChargedRatioKaonToPion2040 = (TH1D*)fileDataALICEChargedRatioKaonToPion->Get("hsys_PbPb276_2040_kaon_to_pion_sum");
	TGraphAsymmErrors* graphChargedRatioKaonToPionSys2040 = new TGraphAsymmErrors(histoSystChargedRatioKaonToPion2040); 


	//**************************** Charged pion ****************************//
	TFile* fileDataALICEChargedPionRAA = new TFile(fileNameChargedPionRAA);
	
	TH1D *histoStatChargedPion0005 = (TH1D*)fileDataALICEChargedPionRAA->Get("RAAPion_Stat_0_5");
	TGraphAsymmErrors* graphChargedPionRAA0005 = new TGraphAsymmErrors(histoStatChargedPion0005); 
	TH1D *histoSystChargedPion0005 = (TH1D*)fileDataALICEChargedPionRAA->Get("RAAPion_Syst_0_5");
	TGraphAsymmErrors* graphChargedPionRAASys0005 = new TGraphAsymmErrors(histoSystChargedPion0005); 
	
	TH1D *histoStatChargedPion0510 = (TH1D*)fileDataALICEChargedPionRAA->Get("RAAPion_Stat_5_10");
	TGraphAsymmErrors* graphChargedPionRAA0510 = new TGraphAsymmErrors(histoStatChargedPion0510); 
	TH1D *histoSystChargedPion0510 = (TH1D*)fileDataALICEChargedPionRAA->Get("RAAPion_Syst_5_10");
	TGraphAsymmErrors* graphChargedPionRAASys0510 = new TGraphAsymmErrors(histoSystChargedPion0510); 

	TH1D* histoStatChargedPion0010 = (TH1D*)histoStatChargedPion0510->Clone("histoStatChargedPion0010");
	TH1D* histoSystChargedPion0010 = (TH1D*)histoSystChargedPion0510->Clone("histoSystChargedPion0010");
	histoStatChargedPion0010->Add(histoStatChargedPion0005);
	histoSystChargedPion0010->Add(histoSystChargedPion0005);
	histoStatChargedPion0010->Scale(0.5);
	histoSystChargedPion0010->Scale(0.5);
	for (Int_t i = 1; i < histoSystChargedPion0010->GetNbinsX()+1; i++){
		Double_t relErrLowerCent = 0;
		if (histoSystChargedPion0005->GetBinContent(i) != 0){
			relErrLowerCent= histoSystChargedPion0005->GetBinError(i)/histoSystChargedPion0005->GetBinContent(i)*100 ;
		}
		Double_t relErrHigherCent = 0;
		if (histoSystChargedPion0510->GetBinContent(i) != 0){
			relErrHigherCent = histoSystChargedPion0510->GetBinError(i)/histoSystChargedPion0510->GetBinContent(i)*100 ;
		}
		
		if (relErrHigherCent > relErrLowerCent){
			histoSystChargedPion0010->SetBinError(i, histoSystChargedPion0010->GetBinContent(i)*relErrHigherCent/100);
		} else {
			histoSystChargedPion0010->SetBinError(i, histoSystChargedPion0010->GetBinContent(i)*relErrLowerCent/100);
		}         
	}   
	TGraphAsymmErrors* graphChargedPionRAA0010 = new TGraphAsymmErrors(histoStatChargedPion0010); 
	TGraphAsymmErrors* graphChargedPionRAASys0010 = new TGraphAsymmErrors(histoSystChargedPion0010); 
	
	TH1D *histoStatChargedPion2040 = (TH1D*)fileDataALICEChargedPionRAA->Get("RAAPion_Stat_20_40");
	TGraphAsymmErrors* graphChargedPionRAA2040 = new TGraphAsymmErrors(histoStatChargedPion2040); 
	TH1D *histoSystChargedPion2040 = (TH1D*)fileDataALICEChargedPionRAA->Get("RAAPion_Syst_20_40");
	TGraphAsymmErrors* graphChargedPionRAASys2040 = new TGraphAsymmErrors(histoSystChargedPion2040); 


	//**************************** Charged kaon ****************************//
	TFile* fileDataALICEChargedKaonRAA = new TFile(fileNameChargedKaonRAA);

	TH1D *histoStatChargedKaon0005 = (TH1D*)fileDataALICEChargedKaonRAA->Get("RAAKaon_Stat_0_5");
	TGraphAsymmErrors* graphChargedKaonRAA0005 = new TGraphAsymmErrors(histoStatChargedKaon0005); 
	TH1D *histoSystChargedKaon0005 = (TH1D*)fileDataALICEChargedKaonRAA->Get("RAAKaon_Syst_0_5");
	TGraphAsymmErrors* graphChargedKaonRAASys0005 = new TGraphAsymmErrors(histoSystChargedKaon0005); 
	
	TH1D *histoStatChargedKaon0510 = (TH1D*)fileDataALICEChargedKaonRAA->Get("RAAKaon_Stat_5_10");
	TGraphAsymmErrors* graphChargedKaonRAA0510 = new TGraphAsymmErrors(histoStatChargedKaon0510); 
	TH1D *histoSystChargedKaon0510 = (TH1D*)fileDataALICEChargedKaonRAA->Get("RAAKaon_Syst_5_10");
	TGraphAsymmErrors* graphChargedKaonRAASys0510 = new TGraphAsymmErrors(histoSystChargedKaon0510); 

	TH1D* histoStatChargedKaon0010 = (TH1D*)histoStatChargedKaon0510->Clone("histoStatChargedKaon0010");
	TH1D* histoSystChargedKaon0010 = (TH1D*)histoSystChargedKaon0510->Clone("histoSystChargedKaon0010");
	histoStatChargedKaon0010->Add(histoStatChargedKaon0005);
	histoSystChargedKaon0010->Add(histoSystChargedKaon0005);
	histoStatChargedKaon0010->Scale(0.5);
	histoSystChargedKaon0010->Scale(0.5);
	for (Int_t i = 1; i < histoSystChargedKaon0010->GetNbinsX()+1; i++){
		Double_t relErrLowerCent = 0;
		if (histoSystChargedKaon0005->GetBinContent(i) != 0){
			relErrLowerCent= histoSystChargedKaon0005->GetBinError(i)/histoSystChargedKaon0005->GetBinContent(i)*100 ;
		}
		Double_t relErrHigherCent = 0;
		if (histoSystChargedKaon0510->GetBinContent(i) != 0){
			relErrHigherCent = histoSystChargedKaon0510->GetBinError(i)/histoSystChargedKaon0510->GetBinContent(i)*100 ;
		}
		
		if (relErrHigherCent > relErrLowerCent){
			histoSystChargedKaon0010->SetBinError(i, histoSystChargedKaon0010->GetBinContent(i)*relErrHigherCent/100);
		} else {
			histoSystChargedKaon0010->SetBinError(i, histoSystChargedKaon0010->GetBinContent(i)*relErrLowerCent/100);
		}         
	}   
	TGraphAsymmErrors* graphChargedKaonRAA0010 = new TGraphAsymmErrors(histoStatChargedKaon0010); 
	TGraphAsymmErrors* graphChargedKaonRAASys0010 = new TGraphAsymmErrors(histoSystChargedKaon0010); 

	TH1D *histoStatChargedKaon2040 = (TH1D*)fileDataALICEChargedKaonRAA->Get("RAAKaon_Stat_20_40");
	TGraphAsymmErrors* graphChargedKaonRAA2040 = new TGraphAsymmErrors(histoStatChargedKaon2040); 
	TH1D *histoSystChargedKaon2040 = (TH1D*)fileDataALICEChargedKaonRAA->Get("RAAKaon_Syst_20_40");
	TGraphAsymmErrors* graphChargedKaonRAASys2040 = new TGraphAsymmErrors(histoSystChargedKaon2040); 
	
	
	
	//*********************************************************************************************************//
	//***************************************     LHC data input     ******************************************//
	//*********************************************************************************************************//
	
	//**************************** PP 2.76TeV ****************************//
	fileFinalResultsPP = 		new TFile(filePPpublished.Data());
    cout << "For the Pi0 in PP (combined analysis) " << endl; // combo of PCM, EMCal, PCM+EMCal, EMCal merged, PHOS
    TDirectoryFile* directoryPi0PP2760GeV         = (TDirectoryFile*)fileFinalResultsPP->Get("Pi02.76TeV"); 

	graphInvSectionCombStatPi02760GeV = (TGraphAsymmErrors*)directoryPi0PP2760GeV->Get("graphInvCrossSectionPi0Comb2760GeVAStatErr"); 
	graphInvSectionCombStatPi02760GeV = ScaleGraph(graphInvSectionCombStatPi02760GeV,1./xSection2760GeVppINEL);
	graphInvSectionCombSysPi02760GeV = (TGraphAsymmErrors*)directoryPi0PP2760GeV->Get("graphInvCrossSectionPi0Comb2760GeVASysErr");
	graphInvSectionCombSysPi02760GeV = ScaleGraph(graphInvSectionCombSysPi02760GeV,1./xSection2760GeVppINEL);
    
    graphInvSectionCombStatPi02760GeVPlot = (TGraphAsymmErrors*)graphInvSectionCombStatPi02760GeV->Clone("graphInvSectionCombStatPi02760GeVPlot"); 
    graphInvSectionCombSysPi02760GeVPlot = (TGraphAsymmErrors*)graphInvSectionCombSysPi02760GeV->Clone("graphInvSectionCombSysPi02760GeVPlot");

    
    cout << "Pi0 in pp 2.76TeV" << endl;    
    TFile *filePi0PPforRAA =      new TFile(fileFinalResultsPi0PPforRAA);
    TF1 *fitInvCrossSectionPi0Comb2760GeV = (TF1*)filePi0PPforRAA->Get("fitInvCrossSectionPi0Comb2760GeV_YShift");
    fitInvCrossSectionPi0Comb2760GeV->SetParameter(0, fitInvCrossSectionPi0Comb2760GeV->GetParameter(0)*1./(xSection2760GeVpp*recalcBarn)*factorToInel);
    TGraphAsymmErrors* graphInvSectionCombPi02760GeVforRAA = (TGraphAsymmErrors*)filePi0PPforRAA->Get("graphInvCrossSectionPi0Comb2760GeV_YShifted");
    graphInvSectionCombPi02760GeVforRAA = ScaleGraph(graphInvSectionCombPi02760GeVforRAA,1./(xSection2760GeVpp*recalcBarn)*factorToInel);
    TGraphAsymmErrors* graphInvSectionPCMPi02760GeVforRAA =       (TGraphAsymmErrors*)filePi0PPforRAA->Get("graphInvCrossSectionPi0PCMStat2760GeV_YShifted");
    graphInvSectionPCMPi02760GeVforRAA = ScaleGraph(graphInvSectionPCMPi02760GeVforRAA,1./(xSection2760GeVpp*recalcBarn)*factorToInel);
    TGraphAsymmErrors* graphInvSectionPCMSysPi02760GeVforRAA =        (TGraphAsymmErrors*)filePi0PPforRAA->Get("graphInvCrossSectionPi0PCMSysForRAA2760GeV_YShifted");
    graphInvSectionPCMSysPi02760GeVforRAA = ScaleGraph(graphInvSectionPCMSysPi02760GeVforRAA,1./(xSection2760GeVpp*recalcBarn)*factorToInel);
    graphInvSectionPCMSysPi02760GeVforRAA->RemovePoint(graphInvSectionPCMSysPi02760GeVforRAA->GetN()-1);
    
// 
//     cout << "For the Pi0 in PP (EMCal) " << endl;
//     graphInvSectionEMCalStatPi02760GeV = (TGraphAsymmErrors*)directoryPi0PP2760GeV->Get("graphInvCrossSectionPi0EMCal2760GeVAStatErr"); 
//     graphInvSectionEMCalStatPi02760GeV = ScaleGraph(graphInvSectionEMCalStatPi02760GeVPlot,1./xSection2760GeVppINEL);
//     graphInvSectionEMCalSysPi02760GeV = (TGraphAsymmErrors*)directoryPi0PP2760GeV->Get("graphInvCrossSectionPi0EMCal2760GeVASysErr");
//     graphInvSectionEMCalSysPi02760GeV = ScaleGraph(graphInvSectionEMCalSysPi02760GeVPlot,1./xSection2760GeVppINEL);
//     
//     graphInvSectionEMCalStatPi02760GeVPlot = (TGraphAsymmErrors*)graphInvSectionEMCalStatPi02760GeV->Clone("graphInvSectionEMCalStatPi02760GeVPlot"); 
//     graphInvSectionEMCalSysPi02760GeVPlot = (TGraphAsymmErrors*)graphInvSectionEMCalSysPi02760GeV->Clone("graphInvSectionEMCalSysPi02760GeVPlot");

    cout << "For the Eta in PP (combined analysis) " << endl; // combo of PCM, EMCal, PCM+EMCal
    TDirectoryFile* directoryEtaPP2760GeV         = (TDirectoryFile*)fileFinalResultsPP->Get("Eta2.76TeV"); 

    graphInvSectionCombStatEta2760GeV = (TGraphAsymmErrors*)directoryEtaPP2760GeV->Get("graphInvCrossSectionEtaComb2760GeVAStatErr"); 
    graphInvSectionCombStatEta2760GeV = ScaleGraph(graphInvSectionCombStatEta2760GeV,1./xSection2760GeVppINEL);
    graphInvSectionCombSysEta2760GeV = (TGraphAsymmErrors*)directoryEtaPP2760GeV->Get("graphInvCrossSectionEtaComb2760GeVASysErr");
    graphInvSectionCombSysEta2760GeV = ScaleGraph(graphInvSectionCombSysEta2760GeV,1./xSection2760GeVppINEL);

    graphInvSectionCombStatEta2760GeVPlot = (TGraphAsymmErrors*)graphInvSectionCombStatEta2760GeV->Clone("graphInvSectionCombStatEta2760GeVPlot"); 
    graphInvSectionCombSysEta2760GeVPlot = (TGraphAsymmErrors*)graphInvSectionCombSysEta2760GeV->Clone("graphInvSectionCombSysEta2760GeVPlot");
    
    cout << "Eta in pp 2.76TeV" << endl;
    TFile *fileEtaPPforRAA =      new TFile(fileFinalResultsEtaPPforRAA); 
    TF1 *fitInvCrossSectionEtaComb2760GeV = (TF1*)fileEtaPPforRAA->Get("fitInvCrossSectionEtaComb2760GeV_YShift");
    fitInvCrossSectionEtaComb2760GeV->SetParameter(0, fitInvCrossSectionEtaComb2760GeV->GetParameter(0)*1./(xSection2760GeVpp*recalcBarn)*factorToInel);
    TGraphAsymmErrors* graphInvSectionCombEta2760GeVforRAA = (TGraphAsymmErrors*)fileEtaPPforRAA->Get("graphInvCrossSectionEtaComb2760GeV_YShifted");
    graphInvSectionCombEta2760GeVforRAA = ScaleGraph(graphInvSectionCombEta2760GeVforRAA,1./(xSection2760GeVpp*recalcBarn)*factorToInel);
    TGraphAsymmErrors* graphInvSectionPCMEta2760GeVforRAA =       (TGraphAsymmErrors*)fileEtaPPforRAA->Get("graphInvCrossSectionEtaPCMStat2760GeV_YShifted");
    graphInvSectionPCMEta2760GeVforRAA = ScaleGraph(graphInvSectionPCMEta2760GeVforRAA,1./(xSection2760GeVpp*recalcBarn)*factorToInel);
    TGraphAsymmErrors* graphInvSectionPCMSysEta2760GeVforRAA =        (TGraphAsymmErrors*)fileEtaPPforRAA->Get("graphInvCrossSectionEtaPCMSysForRAA2760GeV_YShifted");
    graphInvSectionPCMSysEta2760GeVforRAA = ScaleGraph(graphInvSectionPCMSysEta2760GeVforRAA,1./(xSection2760GeVpp*recalcBarn)*factorToInel);

// 
//     cout << "For the Eta in PP (EMCal) " << endl;
//     graphInvSectionEMCalStatEta2760GeV = (TGraphAsymmErrors*)directoryEtaPP2760GeV->Get("graphInvCrossSectionEtaEMCal2760GeVAStatErr"); 
//     graphInvSectionEMCalStatEta2760GeV = ScaleGraph(graphInvSectionEMCalStatEta2760GeVPlot,1./xSection2760GeVppINEL);
//     graphInvSectionEMCalSysEta2760GeV = (TGraphAsymmErrors*)directoryEtaPP2760GeV->Get("graphInvCrossSectionEtaEMCal2760GeVASysErr");
//     graphInvSectionEMCalSysEta2760GeV = ScaleGraph(graphInvSectionEMCalSysEta2760GeVPlot,1./xSection2760GeVppINEL);
//     
//     graphInvSectionEMCalStatEta2760GeVPlot = (TGraphAsymmErrors*)graphInvSectionEMCalStatEta2760GeV->Clone("graphInvSectionEMCalStatEta2760GeVPlot"); 
//     graphInvSectionEMCalSysEta2760GeVPlot = (TGraphAsymmErrors*)graphInvSectionEMCalSysEta2760GeV->Clone("graphInvSectionEMCalSysEta2760GeVPlot");

    
	//**************************** PCM Pb-Pb data 2011 ****************************//
	cout << "Loading PCM histos" << endl;
	TFile* filePCM = new TFile(fileNamePCM.Data());
	
	TH1D* histoPCMNumberOfEventsPbPb2760GeV_0005 			= (TH1D*)filePCM->Get("histoNumberOfEventsPbPb_2.76TeV0-5%");
 	TH1D* histoPCMNumberOfEventsPbPb2760GeV_0010 			= (TH1D*)filePCM->Get("histoNumberOfEventsPbPb_2.76TeV0-10%");
	TH1D* histoPCMNumberOfEventsPbPb2760GeV_2040 			= (TH1D*)filePCM->Get("histoNumberOfEventsPbPb_2.76TeV20-40%");
	TH1D* histoPCMNumberOfEventsPbPb2760GeV_2050 			= (TH1D*)filePCM->Get("histoNumberOfEventsPbPb_2.76TeV20-50%");

	cout << "For the Pi0 in 0-10% " << endl;
	TDirectoryFile* directoryPCMPi0PbPb2760GeV_0010 		= (TDirectoryFile*)filePCM->Get("Pi0_PbPb_2.76TeV_0-10%"); 
		TH1D* histoPCMPi0MassPbPb2760GeV_0010 				= (TH1D*)directoryPCMPi0PbPb2760GeV_0010->Get("MassPi0");
		TH1D* histoPCMPi0FWHMMeVPbPb2760GeV_0010 			= (TH1D*)directoryPCMPi0PbPb2760GeV_0010->Get("FWHMPi0MeV");
		TH1D* histoPCMPi0TrueMassPbPb2760GeV_0010 			= (TH1D*)directoryPCMPi0PbPb2760GeV_0010->Get("TrueMassPi0");
		TH1D* histoPCMPi0TrueFWHMMeVPbPb2760GeV_0010 		= (TH1D*)directoryPCMPi0PbPb2760GeV_0010->Get("TrueFWHMPi0MeV");
		TH1D* histoPCMPi0AccPbPb2760GeV_0010 				= (TH1D*)directoryPCMPi0PbPb2760GeV_0010->Get("Pi0_Acceptance");
		TH1D* histoPCMPi0TrueEffPtPbPb2760GeV_0010 			= (TH1D*)directoryPCMPi0PbPb2760GeV_0010->Get("Pi0_Efficiency");
		
		TH1D* histoPCMPi0InvYieldPbPb2760GeV_0010 					= (TH1D*)directoryPCMPi0PbPb2760GeV_0010->Get("CorrectedYieldPi0");   
		TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeV_0010 	= new TGraphAsymmErrors(histoPCMPi0InvYieldPbPb2760GeV_0010);
			graphPCMPi0InvYieldStatPbPb2760GeV_0010->RemovePoint(0);

        TH1D* histoPCMPi0InvYieldPbPb2760GeVYshifted_0010                   = (TH1D*)directoryPCMPi0PbPb2760GeV_0010->Get("CorrectedYieldBinShiftedPi0");   
        TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_0010  = new TGraphAsymmErrors(histoPCMPi0InvYieldPbPb2760GeVYshifted_0010);
            graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_0010->RemovePoint(0);

		TGraphAsymmErrors* graphPCMPi0InvYieldSysA2760GeV_0010			= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_0010->Get("Pi0SystErrorA");
			graphPCMPi0InvYieldSysA2760GeV_0010->RemovePoint(graphPCMPi0InvYieldSysA2760GeV_0010->GetN()-1);
		
		TGraphAsymmErrors* graphPCMPi0CorrYieldSysErrPbPb2760GeV_0010 	= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_0010->Get("Pi0SystError"); 
			graphPCMPi0CorrYieldSysErrPbPb2760GeV_0010->RemovePoint(graphPCMPi0CorrYieldSysErrPbPb2760GeV_0010->GetN()-1);

		TGraphAsymmErrors* graphPCMPi0InvYieldSysPbPb2760GeV_0010		= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_0010->Get("Pi0SystError");
			graphPCMPi0InvYieldSysPbPb2760GeV_0010->RemovePoint(graphPCMPi0InvYieldSysPbPb2760GeV_0010->GetN()-1);

		TGraphAsymmErrors* graphPCMPi0RAAStatPbPb2760GeV_0010	= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_0010->Get("Pi0RAA");
		TH1D* histoPCMPi0RAAStatPbPb2760GeV_0010 = (TH1D*)GraphAsymErrorsToHist(graphPCMPi0RAAStatPbPb2760GeV_0010,14,"histoPCMPi0RAAStatPbPb2760GeV_0010");
		TGraphAsymmErrors* graphPCMPi0RAASysPbPb2760GeV_0010	= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_0010->Get("Pi0RAASys");

		TH1D* ratioPCMPi0MassMCDiffData_0010 = (TH1D*)histoPCMPi0TrueMassPbPb2760GeV_0010->Clone("ratioPCMPi0MassMCDiffData_0010");
			ratioPCMPi0MassMCDiffData_0010->Sumw2();
			ratioPCMPi0MassMCDiffData_0010->Add(histoPCMPi0MassPbPb2760GeV_0010,-1);
			ratioPCMPi0MassMCDiffData_0010->Divide(ratioPCMPi0MassMCDiffData_0010, histoPCMPi0MassPbPb2760GeV_0010,1.,1.,"");
			ratioPCMPi0MassMCDiffData_0010->Scale(1./mesonMassExpectPi0);

 		TH1D *histoPCMPi0AccTimesEffPbPb2760GeV_0010 = (TH1D*)histoPCMPi0TrueEffPtPbPb2760GeV_0010->Clone("histoPCMPi0AccTimesEffPbPb2760GeV_0010");
			histoPCMPi0AccTimesEffPbPb2760GeV_0010->Sumw2();
			histoPCMPi0AccTimesEffPbPb2760GeV_0010->Multiply(histoPCMPi0AccPbPb2760GeV_0010);
		
			
	cout << "For the Pi0 in 0-5% " << endl;
	TDirectoryFile* directoryPCMPi0PbPb2760GeV_0005 				= (TDirectoryFile*)filePCM->Get("Pi0_PbPb_2.76TeV_0-5%"); 
		TH1D* histoPCMPi0MassPbPb2760GeV_0005 						= (TH1D*)directoryPCMPi0PbPb2760GeV_0005->Get("MassPi0");
		TH1D* histoPCMPi0FWHMMeVPbPb2760GeV_0005 					= (TH1D*)directoryPCMPi0PbPb2760GeV_0005->Get("FWHMPi0MeV");
		TH1D* histoPCMPi0TrueMassPbPb2760GeV_0005 					= (TH1D*)directoryPCMPi0PbPb2760GeV_0005->Get("TrueMassPi0");
		TH1D* histoPCMPi0TrueFWHMMeVPbPb2760GeV_0005 				= (TH1D*)directoryPCMPi0PbPb2760GeV_0005->Get("TrueFWHMPi0MeV");
		TH1D* histoPCMPi0AccPbPb2760GeV_0005 						= (TH1D*)directoryPCMPi0PbPb2760GeV_0005->Get("Pi0_Acceptance");
		TH1D* histoPCMPi0TrueEffPtPbPb2760GeV_0005 					= (TH1D*)directoryPCMPi0PbPb2760GeV_0005->Get("Pi0_Efficiency");
		
		TH1D* histoPCMPi0InvYieldPbPb2760GeV_0005 					= (TH1D*)directoryPCMPi0PbPb2760GeV_0005->Get("CorrectedYieldPi0");   
		TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeV_0005 	= new TGraphAsymmErrors(histoPCMPi0InvYieldPbPb2760GeV_0005);
			graphPCMPi0InvYieldStatPbPb2760GeV_0005->RemovePoint(0);
			
		TGraphAsymmErrors* graphPCMPi0InvYieldSysA2760GeV_0005			= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_0005->Get("Pi0SystErrorA");
			graphPCMPi0InvYieldSysA2760GeV_0005->RemovePoint(graphPCMPi0InvYieldSysA2760GeV_0005->GetN()-1);

		TGraphAsymmErrors* graphPCMPi0CorrYieldSysErrPbPb2760GeV_0005 	= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_0005->Get("Pi0SystError"); 
			graphPCMPi0CorrYieldSysErrPbPb2760GeV_0005->RemovePoint(graphPCMPi0CorrYieldSysErrPbPb2760GeV_0005->GetN()-1);

		TGraphAsymmErrors* graphPCMPi0InvYieldSysPbPb2760GeV_0005		= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_0005->Get("Pi0SystError");
			graphPCMPi0InvYieldSysPbPb2760GeV_0005->RemovePoint(graphPCMPi0InvYieldSysPbPb2760GeV_0005->GetN()-1);

		TGraphAsymmErrors* graphPCMPi0RAAStatPbPb2760GeV_0005	= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_0005->Get("Pi0RAA");
		TGraphAsymmErrors* graphPCMPi0RAASysPbPb2760GeV_0005	= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_0005->Get("Pi0RAASys");

		TH1D* ratioPCMPi0MassMCDiffData_0005 = (TH1D*)histoPCMPi0TrueMassPbPb2760GeV_0005->Clone("ratioPCMPi0MassMCDiffData_0005");
			ratioPCMPi0MassMCDiffData_0005->Sumw2();
			ratioPCMPi0MassMCDiffData_0005->Add(histoPCMPi0MassPbPb2760GeV_0005,-1);
			ratioPCMPi0MassMCDiffData_0005->Divide(ratioPCMPi0MassMCDiffData_0005, histoPCMPi0MassPbPb2760GeV_0005,1.,1.,"");
			ratioPCMPi0MassMCDiffData_0005->Scale(1./mesonMassExpectPi0);

 		TH1D *histoPCMPi0AccTimesEffPbPb2760GeV_0005 = (TH1D*)histoPCMPi0TrueEffPtPbPb2760GeV_0005->Clone("histoPCMPi0AccTimesEffPbPb2760GeV_0005");
			histoPCMPi0AccTimesEffPbPb2760GeV_0005->Sumw2();
			histoPCMPi0AccTimesEffPbPb2760GeV_0005->Multiply(histoPCMPi0AccPbPb2760GeV_0005);

	cout << "For the Pi0 in 20-40% " << endl;
	TDirectoryFile* directoryPCMPi0PbPb2760GeV_2040 				= (TDirectoryFile*)filePCM->Get("Pi0_PbPb_2.76TeV_20-40%"); 
		TH1D* histoPCMPi0MassPbPb2760GeV_2040 						= (TH1D*)directoryPCMPi0PbPb2760GeV_2040->Get("MassPi0");
		TH1D* histoPCMPi0FWHMMeVPbPb2760GeV_2040 					= (TH1D*)directoryPCMPi0PbPb2760GeV_2040->Get("FWHMPi0MeV");
		TH1D* histoPCMPi0TrueMassPbPb2760GeV_2040 					= (TH1D*)directoryPCMPi0PbPb2760GeV_2040->Get("TrueMassPi0");
		TH1D* histoPCMPi0TrueFWHMMeVPbPb2760GeV_2040 				= (TH1D*)directoryPCMPi0PbPb2760GeV_2040->Get("TrueFWHMPi0MeV");
		TH1D* histoPCMPi0AccPbPb2760GeV_2040 						= (TH1D*)directoryPCMPi0PbPb2760GeV_2040->Get("Pi0_Acceptance");
		TH1D* histoPCMPi0TrueEffPtPbPb2760GeV_2040 					= (TH1D*)directoryPCMPi0PbPb2760GeV_2040->Get("Pi0_Efficiency");
		
		TH1D* histoPCMPi0InvYieldPbPb2760GeV_2040 					= (TH1D*)directoryPCMPi0PbPb2760GeV_2040->Get("CorrectedYieldPi0");   
		TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeV_2040 	= new TGraphAsymmErrors(histoPCMPi0InvYieldPbPb2760GeV_2040);
			graphPCMPi0InvYieldStatPbPb2760GeV_2040->RemovePoint(0);
			
		TGraphAsymmErrors* graphPCMPi0InvYieldSysA2760GeV_2040			= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2040->Get("Pi0SystErrorA");
			graphPCMPi0InvYieldSysA2760GeV_2040->RemovePoint(graphPCMPi0InvYieldSysA2760GeV_2040->GetN()-1);

		TGraphAsymmErrors* graphPCMPi0CorrYieldSysErrPbPb2760GeV_2040 	= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2040->Get("Pi0SystError"); 
			graphPCMPi0CorrYieldSysErrPbPb2760GeV_2040->RemovePoint(graphPCMPi0CorrYieldSysErrPbPb2760GeV_2040->GetN()-1);

		TGraphAsymmErrors* graphPCMPi0InvYieldSysPbPb2760GeV_2040		= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2040->Get("Pi0SystError");
			graphPCMPi0InvYieldSysPbPb2760GeV_2040->RemovePoint(graphPCMPi0InvYieldSysPbPb2760GeV_2040->GetN()-1);

		TGraphAsymmErrors* graphPCMPi0RAAStatPbPb2760GeV_2040	= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2040->Get("Pi0RAA");
		TH1D* histoPCMPi0RAAStatPbPb2760GeV_2040 = (TH1D*)GraphAsymErrorsToHist(graphPCMPi0RAAStatPbPb2760GeV_2040,14,"histoPCMPi0RAAStatPbPb2760GeV_2040");
		TGraphAsymmErrors* graphPCMPi0RAASysPbPb2760GeV_2040	= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2040->Get("Pi0RAASys");

		TH1D* ratioPCMPi0MassMCDiffData_2040 = (TH1D*)histoPCMPi0TrueMassPbPb2760GeV_2040->Clone("ratioPCMPi0MassMCDiffData_2040");
			ratioPCMPi0MassMCDiffData_2040->Sumw2();
			ratioPCMPi0MassMCDiffData_2040->Add(histoPCMPi0MassPbPb2760GeV_2040,-1);
			ratioPCMPi0MassMCDiffData_2040->Divide(ratioPCMPi0MassMCDiffData_2040, histoPCMPi0MassPbPb2760GeV_2040,1.,1.,"");
			ratioPCMPi0MassMCDiffData_2040->Scale(1./mesonMassExpectPi0);

 		TH1D *histoPCMPi0AccTimesEffPbPb2760GeV_2040 = (TH1D*)histoPCMPi0TrueEffPtPbPb2760GeV_2040->Clone("histoPCMPi0AccTimesEffPbPb2760GeV_2040");
			histoPCMPi0AccTimesEffPbPb2760GeV_2040->Sumw2();
			histoPCMPi0AccTimesEffPbPb2760GeV_2040->Multiply(histoPCMPi0AccPbPb2760GeV_2040);			
			

	cout << "For the Pi0 in 20-50% " << endl;
	TDirectory* directoryPCMPi0PbPb2760GeV_2050 					= (TDirectory*)filePCM->Get("Pi0_PbPb_2.76TeV_20-50%"); 
		TH1D* histoPCMPi0MassPbPb2760GeV_2050 						= (TH1D*)directoryPCMPi0PbPb2760GeV_2050->Get("MassPi0");
		TH1D* histoPCMPi0FWHMMeVPbPb2760GeV_2050 					= (TH1D*)directoryPCMPi0PbPb2760GeV_2050->Get("FWHMPi0MeV");
		TH1D* histoPCMPi0TrueMassPbPb2760GeV_2050 					= (TH1D*)directoryPCMPi0PbPb2760GeV_2050->Get("TrueMassPi0");
		TH1D* histoPCMPi0TrueFWHMMeVPbPb2760GeV_2050 				= (TH1D*)directoryPCMPi0PbPb2760GeV_2050->Get("TrueFWHMPi0MeV");
		TH1D* histoPCMPi0AccPbPb2760GeV_2050 						= (TH1D*)directoryPCMPi0PbPb2760GeV_2050->Get("Pi0_Acceptance");
		TH1D* histoPCMPi0TrueEffPtPbPb2760GeV_2050 					= (TH1D*)directoryPCMPi0PbPb2760GeV_2050->Get("Pi0_Efficiency");
		
	    TH1D* histoPCMPi0InvYieldPbPb2760GeV_2050 					= (TH1D*)directoryPCMPi0PbPb2760GeV_2050->Get("CorrectedYieldPi0");   
		TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeV_2050 	= new TGraphAsymmErrors(histoPCMPi0InvYieldPbPb2760GeV_2050);
			graphPCMPi0InvYieldStatPbPb2760GeV_2050->RemovePoint(0);
						
		TGraphAsymmErrors* graphPCMPi0InvYieldSysA2760GeV_2050			= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2050->Get("Pi0SystErrorA");
			graphPCMPi0InvYieldSysA2760GeV_2050->RemovePoint(graphPCMPi0InvYieldSysA2760GeV_2050->GetN()-1);
		
        TH1D* histoPCMPi0InvYieldPbPb2760GeVYshifted_2050                  = (TH1D*)directoryPCMPi0PbPb2760GeV_2050->Get("CorrectedYieldBinShiftedPi0");   
        TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_2050  = new TGraphAsymmErrors(histoPCMPi0InvYieldPbPb2760GeVYshifted_2050);
            graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_2050->RemovePoint(0);

		TGraphAsymmErrors* graphPCMPi0CorrYieldSysErrPbPb2760GeV_2050 	= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2050->Get("Pi0SystError"); 
			graphPCMPi0CorrYieldSysErrPbPb2760GeV_2050->RemovePoint(graphPCMPi0CorrYieldSysErrPbPb2760GeV_2050->GetN()-1);
		
		TGraphAsymmErrors* graphPCMPi0InvYieldSysPbPb2760GeV_2050	 	= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2050->Get("Pi0SystError");
			graphPCMPi0InvYieldSysPbPb2760GeV_2050->RemovePoint(graphPCMPi0InvYieldSysPbPb2760GeV_2050->GetN()-1);
				
		TGraphAsymmErrors* graphPCMPi0RAAStatPbPb2760GeV_2050	= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2050->Get("Pi0RAA");
		TH1D* histoPCMPi0RAAStatPbPb2760GeV_2050 = (TH1D*)GraphAsymErrorsToHist(graphPCMPi0RAAStatPbPb2760GeV_2050,14,"histoPCMPi0RAAStatPbPb2760GeV_2050");
		TGraphAsymmErrors* graphPCMPi0RAASysPbPb2760GeV_2050	= (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2050->Get("Pi0RAASys");

		TH1D* ratioPCMPi0MassMCDiffData_2050 = (TH1D*)histoPCMPi0TrueMassPbPb2760GeV_2050->Clone("ratioPCMPi0MassMCDiffData_2050");
			ratioPCMPi0MassMCDiffData_2050->Sumw2();
			ratioPCMPi0MassMCDiffData_2050->Add(histoPCMPi0MassPbPb2760GeV_2050,-1);
			ratioPCMPi0MassMCDiffData_2050->Divide(ratioPCMPi0MassMCDiffData_2050, histoPCMPi0MassPbPb2760GeV_2050,1.,1.,"");
			ratioPCMPi0MassMCDiffData_2050->Scale(1./mesonMassExpectPi0);
		
 		TH1D* histoPCMPi0AccTimesEffPbPb2760GeV_2050 = (TH1D*)histoPCMPi0TrueEffPtPbPb2760GeV_2050->Clone("histoPCMPi0AccTimesEffPbPb2760GeV_2050");
			histoPCMPi0AccTimesEffPbPb2760GeV_2050->Sumw2();
			histoPCMPi0AccTimesEffPbPb2760GeV_2050->Multiply(histoPCMPi0AccPbPb2760GeV_2050);
			
		TGraphAsymmErrors* graphPCMPi0RCPStat2760GeV	= (TGraphAsymmErrors*)filePCM->Get("Pi0RCP");
		TGraphAsymmErrors* graphPCMPi0RCPSys2760GeV	= (TGraphAsymmErrors*)filePCM->Get("Pi0RCPsys");
		

	cout << "For the Eta in 0-10% " << endl;
	TDirectory* directoryPCMEtaPbPb2760GeV_0010 					= (TDirectory*)filePCM->Get("Eta_PbPb_2.76TeV_0-10%"); 
		TH1D* histoPCMEtaMassPbPb2760GeV_0010 						= (TH1D*)directoryPCMEtaPbPb2760GeV_0010->Get("MassEta");
		TH1D* histoPCMEtaFWHMMeVPbPb2760GeV_0010					= (TH1D*)directoryPCMEtaPbPb2760GeV_0010->Get("FWHMEtaMeV");
		TH1D* histoPCMEtaTrueMassPbPb2760GeV_0010 					= (TH1D*)directoryPCMEtaPbPb2760GeV_0010->Get("TrueMassEta");
		TH1D* histoPCMEtaTrueFWHMMeVPbPb2760GeV_0010 				= (TH1D*)directoryPCMEtaPbPb2760GeV_0010->Get("TrueFWHMEtaMeV");
		TH1D* histoPCMEtaAccPbPb2760GeV_0010 						= (TH1D*)directoryPCMEtaPbPb2760GeV_0010->Get("Eta_Acceptance");
		TH1D* histoPCMEtaTrueEffPtPbPb2760GeV_0010 					= (TH1D*)directoryPCMEtaPbPb2760GeV_0010->Get("Eta_Efficiency");
		
	    TH1D* histoPCMEtaInvYieldPbPb2760GeV_0010 					= (TH1D*)directoryPCMEtaPbPb2760GeV_0010->Get("CorrectedYieldEta");   
		TGraphAsymmErrors* graphPCMEtaInvYieldStatPbPb2760GeV_0010 	= new TGraphAsymmErrors(histoPCMEtaInvYieldPbPb2760GeV_0010);
			graphPCMEtaInvYieldStatPbPb2760GeV_0010->RemovePoint(0);
			graphPCMEtaInvYieldStatPbPb2760GeV_0010->RemovePoint(0);
            
        TH1D* histoPCMEtaInvYieldPbPb2760GeVYshifted_0010                  = (TH1D*)directoryPCMEtaPbPb2760GeV_0010->Get("CorrectedYieldBinShiftedEta");   
        TGraphAsymmErrors* graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_0010  = new TGraphAsymmErrors(histoPCMEtaInvYieldPbPb2760GeVYshifted_0010);
            graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_0010->RemovePoint(0);
            graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_0010->RemovePoint(0);

        TGraphAsymmErrors* graphPCMEtaInvYieldSysA2760GeV_0010			= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_0010->Get("EtaSystErrorA");
			graphPCMEtaInvYieldSysA2760GeV_0010->RemovePoint(graphPCMEtaInvYieldSysA2760GeV_0010->GetN()-1);
		
		TGraphAsymmErrors* graphPCMEtaCorrYieldSysErrPbPb2760GeV_0010 	= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_0010->Get("EtaSystError"); 
            graphPCMEtaCorrYieldSysErrPbPb2760GeV_0010->RemovePoint(graphPCMEtaCorrYieldSysErrPbPb2760GeV_0010->GetN()-1);

		TGraphAsymmErrors* graphPCMEtaInvYieldSysPbPb2760GeV_0010		= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_0010->Get("EtaSystError");
            graphPCMEtaInvYieldSysPbPb2760GeV_0010->RemovePoint(graphPCMEtaInvYieldSysPbPb2760GeV_0010->GetN()-1);
			
		TGraphAsymmErrors* graphPCMEtaRAAStatPbPb2760GeV_0010	= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_0010->Get("EtaRAA");
		TH1D* histoPCMEtaRAAStatPbPb2760GeV_0010 = (TH1D*)GraphAsymErrorsToHist(graphPCMEtaRAAStatPbPb2760GeV_0010,10,"histoPCMEtaRAAStatPbPb2760GeV_0010");
		TGraphAsymmErrors* graphPCMEtaRAASysPbPb2760GeV_0010	= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_0010->Get("EtaRAASys");

		TH1D* histoPCMEtatoPi0Stat2760GeV_0010 					= (TH1D*)directoryPCMEtaPbPb2760GeV_0010->Get("EtatoPi0Ratio"); 
		TGraphAsymmErrors* graphPCMEtatoPi0Stat2760GeV_0010	= new TGraphAsymmErrors(histoPCMEtatoPi0Stat2760GeV_0010);
		graphPCMEtatoPi0Stat2760GeV_0010->RemovePoint(0);
		graphPCMEtatoPi0Stat2760GeV_0010->RemovePoint(0);

		TGraphAsymmErrors* graphPCMEtatoPi0Sys2760GeV_0010	= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_0010->Get("EtatoPi0RatioSys");
		graphPCMEtatoPi0Sys2760GeV_0010->RemovePoint(graphPCMEtatoPi0Sys2760GeV_0010->GetN()-1);

        TH1D* ratioPCMEtaMassMCDiffData_0010						= (TH1D*)histoPCMEtaTrueMassPbPb2760GeV_0010->Clone("ratioPCMEtaMassMCDiffData_0010");
			ratioPCMEtaMassMCDiffData_0010->Sumw2();
			ratioPCMEtaMassMCDiffData_0010->Add(histoPCMEtaMassPbPb2760GeV_0010,-1);
			ratioPCMEtaMassMCDiffData_0010->Divide(ratioPCMEtaMassMCDiffData_0010, histoPCMEtaMassPbPb2760GeV_0010,1.,1.,"");
			ratioPCMEtaMassMCDiffData_0010->Scale(1./mesonMassExpectEta);
		
 		TH1D* histoPCMEtaAccTimesEffPbPb2760GeV_0010				= (TH1D*)histoPCMEtaTrueEffPtPbPb2760GeV_0010->Clone("histoPCMEtaAccTimesEffPbPb2760GeV_0010");
			histoPCMEtaAccTimesEffPbPb2760GeV_0010->Sumw2();
			histoPCMEtaAccTimesEffPbPb2760GeV_0010->Multiply(histoPCMEtaAccPbPb2760GeV_0010);
		
		
	cout << "For the Eta in 20-50% " << endl;
	TDirectory* directoryPCMEtaPbPb2760GeV_2050 					= (TDirectory*)filePCM->Get("Eta_PbPb_2.76TeV_20-50%"); 
		TH1D* histoPCMEtaMassPbPb2760GeV_2050 						= (TH1D*)directoryPCMEtaPbPb2760GeV_2050->Get("MassEta");
		TH1D* histoPCMEtaFWHMMeVPbPb2760GeV_2050 					= (TH1D*)directoryPCMEtaPbPb2760GeV_2050->Get("FWHMEtaMeV");
		TH1D* histoPCMEtaTrueMassPbPb2760GeV_2050 					= (TH1D*)directoryPCMEtaPbPb2760GeV_2050->Get("TrueMassEta");
		TH1D* histoPCMEtaTrueFWHMMeVPbPb2760GeV_2050 				= (TH1D*)directoryPCMEtaPbPb2760GeV_2050->Get("TrueFWHMEtaMeV");
		TH1D* histoPCMEtaAccPbPb2760GeV_2050 						= (TH1D*)directoryPCMEtaPbPb2760GeV_2050->Get("Eta_Acceptance");
		TH1D* histoPCMEtaTrueEffPtPbPb2760GeV_2050 					= (TH1D*)directoryPCMEtaPbPb2760GeV_2050->Get("Eta_Efficiency");
		
	    TH1D* histoPCMEtaInvYieldPbPb2760GeV_2050 					= (TH1D*)directoryPCMEtaPbPb2760GeV_2050->Get("CorrectedYieldEta"); 
		TGraphAsymmErrors* graphPCMEtaInvYieldStatPbPb2760GeV_2050 	= new TGraphAsymmErrors(histoPCMEtaInvYieldPbPb2760GeV_2050);
			graphPCMEtaInvYieldStatPbPb2760GeV_2050->RemovePoint(0);
			graphPCMEtaInvYieldStatPbPb2760GeV_2050->RemovePoint(0);
			
        TH1D* histoPCMEtaInvYieldPbPb2760GeVYshifted_2050                  = (TH1D*)directoryPCMEtaPbPb2760GeV_2050->Get("CorrectedYieldBinShiftedEta");   
        TGraphAsymmErrors* graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_2050  = new TGraphAsymmErrors(histoPCMEtaInvYieldPbPb2760GeVYshifted_2050);
            graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_2050->RemovePoint(0);
            graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_2050->RemovePoint(0);

        TGraphAsymmErrors* graphPCMEtaInvYieldSysA2760GeV_2050	= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_2050->Get("EtaSystErrorA");
        graphPCMEtaInvYieldSysA2760GeV_2050->RemovePoint(graphPCMEtaInvYieldSysA2760GeV_2050->GetN()-1);

		TGraphAsymmErrors* graphPCMEtaCorrYieldSysErrPbPb2760GeV_2050 	= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_2050->Get("EtaSystError"); 
        graphPCMEtaCorrYieldSysErrPbPb2760GeV_2050->RemovePoint(graphPCMEtaCorrYieldSysErrPbPb2760GeV_2050->GetN()-1);
			
		TGraphAsymmErrors* graphPCMEtaInvYieldSysPbPb2760GeV_2050		= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_2050->Get("EtaSystError");
        graphPCMEtaInvYieldSysPbPb2760GeV_2050->RemovePoint(graphPCMEtaInvYieldSysPbPb2760GeV_2050->GetN()-1);

		TGraphAsymmErrors* graphPCMEtaRAAStatPbPb2760GeV_2050	= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_2050->Get("EtaRAA");
		TH1D* histoPCMEtaRAAStatPbPb2760GeV_2050 = (TH1D*)GraphAsymErrorsToHist(graphPCMEtaRAAStatPbPb2760GeV_2050,10,"histoPCMEtaRAAStatPbPb2760GeV_2050");
		TGraphAsymmErrors* graphPCMEtaRAASysPbPb2760GeV_2050	= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_2050->Get("EtaRAASys");

		TH1D* histoPCMEtatoPi0Stat2760GeV_2050 					= (TH1D*)directoryPCMEtaPbPb2760GeV_2050->Get("EtatoPi0Ratio"); 
		TGraphAsymmErrors* graphPCMEtatoPi0Stat2760GeV_2050	= new TGraphAsymmErrors(histoPCMEtatoPi0Stat2760GeV_2050);
			graphPCMEtatoPi0Stat2760GeV_2050->RemovePoint(0);
			graphPCMEtatoPi0Stat2760GeV_2050->RemovePoint(0);

		TGraphAsymmErrors* graphPCMEtatoPi0Sys2760GeV_2050	= (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_2050->Get("EtatoPi0RatioSys");
			graphPCMEtatoPi0Sys2760GeV_2050->RemovePoint(graphPCMEtatoPi0Sys2760GeV_2050->GetN()-1);

		TH1D* ratioPCMEtaMassMCDiffData_2050						= (TH1D*)histoPCMEtaTrueMassPbPb2760GeV_2050->Clone("ratioPCMEtaMassMCDiffData_2050");
			ratioPCMEtaMassMCDiffData_2050->Sumw2();
			ratioPCMEtaMassMCDiffData_2050->Add(histoPCMEtaMassPbPb2760GeV_2050,-1);
			ratioPCMEtaMassMCDiffData_2050->Divide(ratioPCMEtaMassMCDiffData_2050, histoPCMEtaMassPbPb2760GeV_2050,1.,1.,"");
			ratioPCMEtaMassMCDiffData_2050->Scale(1./mesonMassExpectEta);
		
 		TH1D* histoPCMEtaAccTimesEffPbPb2760GeV_2050				= (TH1D*)histoPCMEtaTrueEffPtPbPb2760GeV_2050->Clone("histoPCMEtaAccTimesEffPbPb2760GeV_2050");
			histoPCMEtaAccTimesEffPbPb2760GeV_2050->Sumw2();
			histoPCMEtaAccTimesEffPbPb2760GeV_2050->Multiply(histoPCMEtaAccPbPb2760GeV_2050);	
			
	
	
	//**************************** EMCal Pb-Pb data 2011 ****************************//
	cout << "Loading EMCal histos" << endl;	
	TFile* fileEMCal								= new TFile(Form("LHC11hExternalInputs/%s",fileNameEMCalFull.Data()));
	
	//for EMCal file with histos, used naming below
	TDirectory* directoryEMCalPi0PbPb2760GeV 				= (TDirectory*)fileEMCal->Get("Pi02.76TeV_PbPb");
	cout << "Pi0 0-10%" << endl;
		TH1D* histoEMCalPi0InvYieldPbPb2760GeV_0010			= (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("InvYieldPbPbStatErrPi_0010");
		TGraphAsymmErrors* graphEMCalPi0InvYieldStatPbPb2760GeV_0010 		= new TGraphAsymmErrors(histoEMCalPi0InvYieldPbPb2760GeV_0010);
			graphEMCalPi0InvYieldStatPbPb2760GeV_0010->RemovePoint(graphEMCalPi0InvYieldStatPbPb2760GeV_0010->GetN()-1);
			graphEMCalPi0InvYieldStatPbPb2760GeV_0010->RemovePoint(graphEMCalPi0InvYieldStatPbPb2760GeV_0010->GetN()-1);

		TH1D* histoEMCalPi0InvYieldSysPbPb2760GeV_0010 	= (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("InvYieldPbPbSysErrPi_0010");
		TGraphAsymmErrors* graphEMCalPi0CorrYieldSysErrPbPb2760GeV_0010 	= new TGraphAsymmErrors(histoEMCalPi0InvYieldSysPbPb2760GeV_0010); 
			graphEMCalPi0CorrYieldSysErrPbPb2760GeV_0010->RemovePoint(graphEMCalPi0CorrYieldSysErrPbPb2760GeV_0010->GetN()-1);
			graphEMCalPi0CorrYieldSysErrPbPb2760GeV_0010->RemovePoint(graphEMCalPi0CorrYieldSysErrPbPb2760GeV_0010->GetN()-1);

		TGraphAsymmErrors* graphEMCalPi0InvYieldSysPbPb2760GeV_0010 	= new TGraphAsymmErrors(histoEMCalPi0InvYieldSysPbPb2760GeV_0010);
			graphEMCalPi0InvYieldSysPbPb2760GeV_0010->RemovePoint(graphEMCalPi0InvYieldSysPbPb2760GeV_0010->GetN()-1);
			graphEMCalPi0InvYieldSysPbPb2760GeV_0010->RemovePoint(graphEMCalPi0InvYieldSysPbPb2760GeV_0010->GetN()-1);

		TH1D*	histoEMCalEtatoPi0StatPbPb2760GeV_0010 	= (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("StatErrEtatoPi0Ratio_0010");
		TGraphAsymmErrors* graphEMCalEtatoPi0Stat2760GeV_0010 	= new TGraphAsymmErrors(histoEMCalEtatoPi0StatPbPb2760GeV_0010);
			graphEMCalEtatoPi0Stat2760GeV_0010->RemovePoint(graphEMCalEtatoPi0Stat2760GeV_0010->GetN()-1);
			graphEMCalEtatoPi0Stat2760GeV_0010->RemovePoint(graphEMCalEtatoPi0Stat2760GeV_0010->GetN()-1);

		TH1D*	histoEMCalEtatoPi0SysPbPb2760GeV_0010 	= (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("SysErrEtatoPi0Ratio_0010");
		TGraphAsymmErrors* graphEMCalEtatoPi0Sys2760GeV_0010 	= new TGraphAsymmErrors(histoEMCalEtatoPi0SysPbPb2760GeV_0010);
			graphEMCalEtatoPi0Sys2760GeV_0010->RemovePoint(graphEMCalEtatoPi0Sys2760GeV_0010->GetN()-1);
			graphEMCalEtatoPi0Sys2760GeV_0010->RemovePoint(graphEMCalEtatoPi0Sys2760GeV_0010->GetN()-1);
		
	cout << "Pi0 20-50%" << endl;		
		TH1D* histoEMCalPi0InvYieldPbPb2760GeV_2050			= (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("InvYieldPbPbStatErrPi_2050");
		TGraphAsymmErrors* graphEMCalPi0InvYieldStatPbPb2760GeV_2050 		= new TGraphAsymmErrors(histoEMCalPi0InvYieldPbPb2760GeV_2050);
			graphEMCalPi0InvYieldStatPbPb2760GeV_2050->RemovePoint(graphEMCalPi0InvYieldStatPbPb2760GeV_2050->GetN()-1);
			graphEMCalPi0InvYieldStatPbPb2760GeV_2050->RemovePoint(graphEMCalPi0InvYieldStatPbPb2760GeV_2050->GetN()-1);

		TH1D* histoEMCalPi0InvYieldSysPbPb2760GeV_2050 	= (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("InvYieldPbPbSysErrPi_2050");
		TGraphAsymmErrors* graphEMCalPi0CorrYieldSysErrPbPb2760GeV_2050 	= new TGraphAsymmErrors(histoEMCalPi0InvYieldSysPbPb2760GeV_2050); 
			graphEMCalPi0CorrYieldSysErrPbPb2760GeV_2050->RemovePoint(graphEMCalPi0CorrYieldSysErrPbPb2760GeV_2050->GetN()-1);
			graphEMCalPi0CorrYieldSysErrPbPb2760GeV_2050->RemovePoint(graphEMCalPi0CorrYieldSysErrPbPb2760GeV_2050->GetN()-1);

		TGraphAsymmErrors* graphEMCalPi0InvYieldSysPbPb2760GeV_2050 	= new TGraphAsymmErrors(histoEMCalPi0InvYieldSysPbPb2760GeV_2050);
			graphEMCalPi0InvYieldSysPbPb2760GeV_2050->RemovePoint(graphEMCalPi0InvYieldSysPbPb2760GeV_2050->GetN()-1);
			graphEMCalPi0InvYieldSysPbPb2760GeV_2050->RemovePoint(graphEMCalPi0InvYieldSysPbPb2760GeV_2050->GetN()-1);

		TH1D*	histoEMCalEtatoPi0StatPbPb2760GeV_2050 	= (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("StatErrEtatoPi0Ratio_2050");
		TGraphAsymmErrors* graphEMCalEtatoPi0Stat2760GeV_2050 	= new TGraphAsymmErrors(histoEMCalEtatoPi0StatPbPb2760GeV_2050);
			graphEMCalEtatoPi0Stat2760GeV_2050->RemovePoint(graphEMCalEtatoPi0Stat2760GeV_2050->GetN()-1);
			graphEMCalEtatoPi0Stat2760GeV_2050->RemovePoint(graphEMCalEtatoPi0Stat2760GeV_2050->GetN()-1);

		TH1D*	histoEMCalEtatoPi0SysPbPb2760GeV_2050 	= (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("SysErrEtatoPi0Ratio_2050");
		TGraphAsymmErrors* graphEMCalEtatoPi0Sys2760GeV_2050 	= new TGraphAsymmErrors(histoEMCalEtatoPi0SysPbPb2760GeV_2050);
			graphEMCalEtatoPi0Sys2760GeV_2050->RemovePoint(graphEMCalEtatoPi0Sys2760GeV_2050->GetN()-1);
			graphEMCalEtatoPi0Sys2760GeV_2050->RemovePoint(graphEMCalEtatoPi0Sys2760GeV_2050->GetN()-1);

		
	TDirectory* directoryEMCalEtaPbPb2760GeV 				= (TDirectory*)fileEMCal->Get("Eta2.76TeV_PbPb");
	cout << "Eta 0-10%" << endl;
		TH1D* histoEMCalEtaInvYieldPbPb2760GeV_0010			= (TH1D*)directoryEMCalEtaPbPb2760GeV->Get("InvYieldPbPbStatErrEta_0010");
		TGraphAsymmErrors* graphEMCalEtaInvYieldStatPbPb2760GeV_0010 		= new TGraphAsymmErrors(histoEMCalEtaInvYieldPbPb2760GeV_0010);
			graphEMCalEtaInvYieldStatPbPb2760GeV_0010->RemovePoint(graphEMCalEtaInvYieldStatPbPb2760GeV_0010->GetN()-1);
			graphEMCalEtaInvYieldStatPbPb2760GeV_0010->RemovePoint(graphEMCalEtaInvYieldStatPbPb2760GeV_0010->GetN()-1);

		TH1D*	histoEMCalEtaInvYieldSysPbPb2760GeV_0010 	= (TH1D*)directoryEMCalEtaPbPb2760GeV->Get("InvYieldPbPbSysErrEta_0010");
		TGraphAsymmErrors* graphEMCalEtaCorrYieldSysErrPbPb2760GeV_0010 	= new TGraphAsymmErrors(histoEMCalEtaInvYieldSysPbPb2760GeV_0010);
			graphEMCalEtaCorrYieldSysErrPbPb2760GeV_0010->RemovePoint(graphEMCalEtaCorrYieldSysErrPbPb2760GeV_0010->GetN()-1);
			graphEMCalEtaCorrYieldSysErrPbPb2760GeV_0010->RemovePoint(graphEMCalEtaCorrYieldSysErrPbPb2760GeV_0010->GetN()-1);

		TGraphAsymmErrors* graphEMCalEtaInvYieldSysPbPb2760GeV_0010 	= new TGraphAsymmErrors(histoEMCalEtaInvYieldSysPbPb2760GeV_0010);
			graphEMCalEtaInvYieldSysPbPb2760GeV_0010->RemovePoint(graphEMCalEtaInvYieldSysPbPb2760GeV_0010->GetN()-1);
			graphEMCalEtaInvYieldSysPbPb2760GeV_0010->RemovePoint(graphEMCalEtaInvYieldSysPbPb2760GeV_0010->GetN()-1);

	cout << "Eta 20-50%" << endl;		
		TH1D* histoEMCalEtaInvYieldPbPb2760GeV_2050			= (TH1D*)directoryEMCalEtaPbPb2760GeV->Get("InvYieldPbPbStatErrEta_2050");
		TGraphAsymmErrors* graphEMCalEtaInvYieldStatPbPb2760GeV_2050 		= new TGraphAsymmErrors(histoEMCalEtaInvYieldPbPb2760GeV_2050);
			graphEMCalEtaInvYieldStatPbPb2760GeV_2050->RemovePoint(graphEMCalEtaInvYieldStatPbPb2760GeV_2050->GetN()-1);
			graphEMCalEtaInvYieldStatPbPb2760GeV_2050->RemovePoint(graphEMCalEtaInvYieldStatPbPb2760GeV_2050->GetN()-1);

		TH1D*	histoEMCalEtaInvYieldSysPbPb2760GeV_2050 	= (TH1D*)directoryEMCalEtaPbPb2760GeV->Get("InvYieldPbPbSysErrEta_2050");
		TGraphAsymmErrors* graphEMCalEtaCorrYieldSysErrPbPb2760GeV_2050 	= new TGraphAsymmErrors(histoEMCalEtaInvYieldSysPbPb2760GeV_2050);
			graphEMCalEtaCorrYieldSysErrPbPb2760GeV_2050->RemovePoint(graphEMCalEtaCorrYieldSysErrPbPb2760GeV_2050->GetN()-1);
			graphEMCalEtaCorrYieldSysErrPbPb2760GeV_2050->RemovePoint(graphEMCalEtaCorrYieldSysErrPbPb2760GeV_2050->GetN()-1);

		TGraphAsymmErrors* graphEMCalEtaInvYieldSysPbPb2760GeV_2050 	= new TGraphAsymmErrors(histoEMCalEtaInvYieldSysPbPb2760GeV_2050);
			graphEMCalEtaInvYieldSysPbPb2760GeV_2050->RemovePoint(graphEMCalEtaInvYieldSysPbPb2760GeV_2050->GetN()-1);
			graphEMCalEtaInvYieldSysPbPb2760GeV_2050->RemovePoint(graphEMCalEtaInvYieldSysPbPb2760GeV_2050->GetN()-1);
		
		
	//Raa files from EMCal <----- with mt scaled spectra
	TFile* fileEMCalRaa								= new TFile("LHC11hExternalInputs/ppReferenceSpectraEMCal27August2015.root");
		TH1D*	histoEMCalPi0RAAStatPbPb2760GeV_0010 	= (TH1D*)fileEMCalRaa->Get("RAAStatPion010");
		TGraphAsymmErrors* graphEMCalPi0RAAStatPbPb2760GeV_0010 	= new TGraphAsymmErrors(histoEMCalPi0RAAStatPbPb2760GeV_0010);
			
		TH1D*	histoEMCalPi0RAASystPbPb2760GeV_0010 	= (TH1D*)fileEMCalRaa->Get("RAASysPion010");
		TGraphAsymmErrors* graphEMCalPi0RAASysPbPb2760GeV_0010 	= new TGraphAsymmErrors(histoEMCalPi0RAASystPbPb2760GeV_0010);

		TH1D*	histoEMCalPi0RAAStatPbPb2760GeV_2050 	= (TH1D*)fileEMCalRaa->Get("RAAStatPion2050");
		TGraphAsymmErrors* graphEMCalPi0RAAStatPbPb2760GeV_2050 	= new TGraphAsymmErrors(histoEMCalPi0RAAStatPbPb2760GeV_2050);
			
		TH1D*	histoEMCalPi0RAASystPbPb2760GeV_2050 	= (TH1D*)fileEMCalRaa->Get("RAASysPion2050");
		TGraphAsymmErrors* graphEMCalPi0RAASysPbPb2760GeV_2050 	= new TGraphAsymmErrors(histoEMCalPi0RAASystPbPb2760GeV_2050);


		TH1D*	histoEMCalEtaRAAStatPbPb2760GeV_0010 	= (TH1D*)fileEMCalRaa->Get("RAAStatEta010");
		TGraphAsymmErrors* graphEMCalEtaRAAStatPbPb2760GeV_0010 	= new TGraphAsymmErrors(histoEMCalEtaRAAStatPbPb2760GeV_0010);
			
		TH1D*	histoEMCalEtaRAASystPbPb2760GeV_0010 	= (TH1D*)fileEMCalRaa->Get("RAASysEta010");
		TGraphAsymmErrors* graphEMCalEtaRAASysPbPb2760GeV_0010 	= new TGraphAsymmErrors(histoEMCalEtaRAASystPbPb2760GeV_0010);

		TH1D*	histoEMCalEtaRAAStatPbPb2760GeV_2050 	= (TH1D*)fileEMCalRaa->Get("RAAStatEta2050");
		TGraphAsymmErrors* graphEMCalEtaRAAStatPbPb2760GeV_2050 	= new TGraphAsymmErrors(histoEMCalEtaRAAStatPbPb2760GeV_2050);
			
		TH1D*	histoEMCalEtaRAASystPbPb2760GeV_2050 	= (TH1D*)fileEMCalRaa->Get("RAASysEta2050");
		TGraphAsymmErrors* graphEMCalEtaRAASysPbPb2760GeV_2050 	= new TGraphAsymmErrors(histoEMCalEtaRAASystPbPb2760GeV_2050);
		
		
		
		
	// *******************************************************************************************************
	// ************************** Combination of different measurements **************************************
	// *******************************************************************************************************
	// REMARKS: 
	// 		- order of measurements defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
	//		- correlations are defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
	// 		- currently only PCM - EMCal vs others fully implemeted energy independent
	// 		- extendable to other energies
	//		- offsets have to be determined manually, see cout's in shell from combination function, more can be uncommented
		
	
	// definition of array of histograms (NULL - means we have no measurement at this energy for this rec-method)
	// for statistical error and final value from respective method

	for (Int_t i = 0; i< 11; i++){
		statErrorCollectionLHC11h_0010[i] = NULL;
		statErrorCollectionLHC11h_2050[i] = NULL;
		statErrorCollectionEtatoPi0LHC11h_0010[i] = NULL;
		statErrorCollectionEtatoPi0LHC11h_2050[i] = NULL;
		statErrorCollectionRaaLHC11h_0010[i] = NULL;
		statErrorCollectionRaaLHC11h_2050[i] = NULL;
	}	
	if(meson.CompareTo("Pi0")==0){
		// for combination PCM - EMCal 2011
		statErrorCollectionLHC11h_0010[0] = (TH1D*)histoPCMPi0InvYieldPbPb2760GeV_0010->Clone("statErrPCMPi0_0010");
		statErrorCollectionLHC11h_0010[2] = (TH1D*)histoEMCalPi0InvYieldPbPb2760GeV_0010->Clone("statErrEMCalPi0_0010");

		statErrorCollectionLHC11h_2050[0] = (TH1D*)histoPCMPi0InvYieldPbPb2760GeV_2050->Clone("statErrPCMPi0_2050");
		statErrorCollectionLHC11h_2050[2] = (TH1D*)histoEMCalPi0InvYieldPbPb2760GeV_2050->Clone("statErrEMCalPi0_2050");

		statErrorCollectionRaaLHC11h_0010[0] = (TH1D*)histoPCMPi0RAAStatPbPb2760GeV_0010->Clone("statErrPCMPi0RAA_0010");
		statErrorCollectionRaaLHC11h_0010[2] = (TH1D*)histoEMCalPi0RAAStatPbPb2760GeV_0010->Clone("statErrEMCalPi0RAA_0010");		

		statErrorCollectionRaaLHC11h_2050[0] = (TH1D*)histoPCMPi0RAAStatPbPb2760GeV_2050->Clone("statErrPCMPi0RAA_2050");
		statErrorCollectionRaaLHC11h_2050[2] = (TH1D*)histoEMCalPi0RAAStatPbPb2760GeV_2050->Clone("statErrEMCalPi0RAA_2050");		
 
	} else if(meson.CompareTo("Eta")==0) {
		
		statErrorCollectionLHC11h_0010[0] = (TH1D*)histoPCMEtaInvYieldPbPb2760GeV_0010->Clone("statErrPCMEta_0010");
		statErrorCollectionLHC11h_0010[2] = (TH1D*)histoEMCalEtaInvYieldPbPb2760GeV_0010->Clone("statErrEMCalEta_0010");
		
		statErrorCollectionLHC11h_2050[0] = (TH1D*)histoPCMEtaInvYieldPbPb2760GeV_2050->Clone("statErrPCMEta_2050");
		statErrorCollectionLHC11h_2050[2] = (TH1D*)histoEMCalEtaInvYieldPbPb2760GeV_2050->Clone("statErrEMCalEta_2050");

		statErrorCollectionEtatoPi0LHC11h_0010[0] = (TH1D*)histoPCMEtatoPi0Stat2760GeV_0010->Clone("statErrPCMEtatoPi0_0010");
		statErrorCollectionEtatoPi0LHC11h_0010[2] = (TH1D*)histoEMCalEtatoPi0StatPbPb2760GeV_0010->Clone("statErrEMCalEtatoPi0_0010");		

		statErrorCollectionEtatoPi0LHC11h_2050[0] = (TH1D*)histoPCMEtatoPi0Stat2760GeV_2050->Clone("statErrPCMEtatoPi0_2050");
		statErrorCollectionEtatoPi0LHC11h_2050[2] = (TH1D*)histoEMCalEtatoPi0StatPbPb2760GeV_2050->Clone("statErrEMCalEtatoPi0_2050");	
		
		statErrorCollectionRaaLHC11h_0010[0] = (TH1D*)histoPCMEtaRAAStatPbPb2760GeV_0010->Clone("statErrPCMEtaRAA_0010");
		statErrorCollectionRaaLHC11h_0010[2] = (TH1D*)histoEMCalEtaRAAStatPbPb2760GeV_0010->Clone("statErrEMCalEtaRAA_0010");		

		statErrorCollectionRaaLHC11h_2050[0] = (TH1D*)histoPCMEtaRAAStatPbPb2760GeV_2050->Clone("statErrPCMEtaRAA_2050");
		statErrorCollectionRaaLHC11h_2050[2] = (TH1D*)histoEMCalEtaRAAStatPbPb2760GeV_2050->Clone("statErrEMCalEtaRAA_2050");		

		
	}
	
	// definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)	
	// for systematic error from respective method
	
	for (Int_t i = 0; i< 11; i++){
		sysErrorCollectionLHC11h_0010[i] = NULL;
		sysErrorCollectionLHC11h_2050[i] = NULL;

		sysErrorCollectionEtatoPi0LHC11h_0010[i] = NULL;
		sysErrorCollectionEtatoPi0LHC11h_2050[i] = NULL;

		sysErrorCollectionRaaLHC11h_0010[i] = NULL;
		sysErrorCollectionRaaLHC11h_2050[i] = NULL;

	}	
	if(meson.CompareTo("Pi0")==0){
		// for combination PCM - EMCal 2011
		sysErrorCollectionLHC11h_0010[0] = (TGraphAsymmErrors*)graphPCMPi0CorrYieldSysErrPbPb2760GeV_0010->Clone("sysErrPCMPi0_0010");
		sysErrorCollectionLHC11h_0010[2] = (TGraphAsymmErrors*)graphEMCalPi0CorrYieldSysErrPbPb2760GeV_0010->Clone("sysErrEMCalPi0_0010");

		sysErrorCollectionLHC11h_2050[0] = (TGraphAsymmErrors*)graphPCMPi0CorrYieldSysErrPbPb2760GeV_2050->Clone("sysErrPCMPi0_2050");
		sysErrorCollectionLHC11h_2050[2] = (TGraphAsymmErrors*)graphEMCalPi0CorrYieldSysErrPbPb2760GeV_2050->Clone("sysErrEMCalPi0_2050");

		sysErrorCollectionRaaLHC11h_0010[0] = (TGraphAsymmErrors*)graphPCMPi0RAASysPbPb2760GeV_0010->Clone("sysErrPCMPi0Raa_0010");
		sysErrorCollectionRaaLHC11h_0010[2] = (TGraphAsymmErrors*)graphEMCalPi0RAASysPbPb2760GeV_0010->Clone("sysErrEMCalPi0Raa_0010");	

		sysErrorCollectionRaaLHC11h_2050[0] = (TGraphAsymmErrors*)graphPCMPi0RAASysPbPb2760GeV_2050->Clone("sysErrPCMPi0Raa_2050");
		sysErrorCollectionRaaLHC11h_2050[2] = (TGraphAsymmErrors*)graphEMCalPi0RAASysPbPb2760GeV_2050->Clone("sysErrEMCalPi0Raa_2050");	
		
	} else if(meson.CompareTo("Eta")==0) {
		sysErrorCollectionLHC11h_0010[0] = (TGraphAsymmErrors*)graphPCMEtaCorrYieldSysErrPbPb2760GeV_0010->Clone("sysErrPCMEta_0010");
		sysErrorCollectionLHC11h_0010[2] = (TGraphAsymmErrors*)graphEMCalEtaCorrYieldSysErrPbPb2760GeV_0010->Clone("sysErrEMCalEta_0010");	
		
		sysErrorCollectionLHC11h_2050[0] = (TGraphAsymmErrors*)graphPCMEtaCorrYieldSysErrPbPb2760GeV_2050->Clone("sysErrPCMEta_2050");
		sysErrorCollectionLHC11h_2050[2] = (TGraphAsymmErrors*)graphEMCalEtaCorrYieldSysErrPbPb2760GeV_2050->Clone("sysErrEMCalEta_2050");

		sysErrorCollectionEtatoPi0LHC11h_0010[0] = (TGraphAsymmErrors*)graphPCMEtatoPi0Sys2760GeV_0010->Clone("sysErrPCMEtatoPi0_0010");
		sysErrorCollectionEtatoPi0LHC11h_0010[2] = (TGraphAsymmErrors*)graphEMCalEtatoPi0Sys2760GeV_0010->Clone("sysErrEMCalEtatoPi0_0010");	

		sysErrorCollectionEtatoPi0LHC11h_2050[0] = (TGraphAsymmErrors*)graphPCMEtatoPi0Sys2760GeV_2050->Clone("sysErrPCMEtatoPi0_2050");
		sysErrorCollectionEtatoPi0LHC11h_2050[2] = (TGraphAsymmErrors*)graphEMCalEtatoPi0Sys2760GeV_2050->Clone("sysErrEMCalEtatoPi0_2050");	

		sysErrorCollectionRaaLHC11h_0010[0] = (TGraphAsymmErrors*)graphPCMEtaRAASysPbPb2760GeV_0010->Clone("sysErrPCMEtaRaa_0010");
		sysErrorCollectionRaaLHC11h_0010[2] = (TGraphAsymmErrors*)graphEMCalEtaRAASysPbPb2760GeV_0010->Clone("sysErrEMCalEtaRaa_0010");	

		sysErrorCollectionRaaLHC11h_2050[0] = (TGraphAsymmErrors*)graphPCMEtaRAASysPbPb2760GeV_2050->Clone("sysErrPCMEtaRaa_2050");
		sysErrorCollectionRaaLHC11h_2050[2] = (TGraphAsymmErrors*)graphEMCalEtaRAASysPbPb2760GeV_2050->Clone("sysErrEMCalEtaRaa_2050");	

	}
	
								  
	Int_t textSizeLabelsPixel = 900*0.04;
	TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasWeights, 0.08, 0.02, 0.035, 0.09);
	canvasWeights->SetLogx();

	TH2F * histo2DWeights = new TH2F("histo2DWeights","histo2DWeights",11000,0.23,70.,1000,-0.5,1.1);
	SetStyleHistoTH2ForGraphs(histo2DWeights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
	histo2DWeights->GetXaxis()->SetMoreLogLabels();
	histo2DWeights->GetXaxis()->SetLabelOffset(-0.01);
	histo2DWeights->Draw("copy");

	TLatex *labelWeightsEnergy = new TLatex(0.7,0.20,collisionSystem2760GeV.Data());
	SetStyleTLatex( labelWeightsEnergy, 0.85*textSizeLabelsPixel,4);
	labelWeightsEnergy->SetTextFont(43);
	
	TLatex *labelWeightsPi0;
	if(meson.CompareTo("Pi0")==0){
		labelWeightsPi0 = new TLatex(0.7,0.16,"#pi^{0} #rightarrow #gamma#gamma");
	} else 	if(meson.CompareTo("Eta")==0){
		labelWeightsPi0 = new TLatex(0.7,0.16,"#eta #rightarrow #gamma#gamma");
	}
	SetStyleTLatex( labelWeightsPi0, 0.85*textSizeLabelsPixel,4);
	labelWeightsPi0->SetTextFont(43);
		
	//**********************************************************************************************************************//
	//******************************************* Assuming maximal correlation *********************************************//
	//**********************************************************************************************************************//
			
	for (Int_t i = 0; i< 11; i++){
			graphWeightsALHC11h_0010[i] = NULL;
			graphWeightsALHC11h_2050[i] = NULL;
	}       

	TString fileNameOutputWeightingALHC11h_0010                                           = Form("%s/0010LHC11h_WeightingMethodA%s.dat",outputDir.Data(),meson.Data());
	TString fileNameOutputWeightingALHC11h_2050                                           = Form("%s/2050LHC11h_WeightingMethodA%s.dat",outputDir.Data(),meson.Data());
	
	if(meson.CompareTo("Pi0")==0){
		// Declaration & calculation of combined spectrum

		graphCombInvYieldTot2760GeVALHC11h_0010       = CombinePtPointsSpectraFullCorrMat( statErrorCollectionLHC11h_0010, sysErrorCollectionLHC11h_0010,
																						   xPtLimitsPi0, 23, offSetsPi0, offSetsPi0Sys,
																						   graphCombInvYieldStat2760GeVALHC11h_0010, graphCombInvYieldSys2760GeVALHC11h_0010,
																						   fileNameOutputWeightingALHC11h_0010,1 );

		graphCombInvYieldTot2760GeVALHC11h_2050       = CombinePtPointsSpectraFullCorrMat( statErrorCollectionLHC11h_2050, sysErrorCollectionLHC11h_2050,
																						   xPtLimitsPi0, 23, offSetsPi0, offSetsPi0Sys,
																						   graphCombInvYieldStat2760GeVALHC11h_2050, graphCombInvYieldSys2760GeVALHC11h_2050,
																						   fileNameOutputWeightingALHC11h_2050,1 );
		graphCombInvYieldStat2760GeVALHC11h_0010->RemovePoint(0);
		graphCombInvYieldStat2760GeVALHC11h_2050->RemovePoint(0);

		graphCombInvYieldSys2760GeVALHC11h_0010->RemovePoint(0);
		graphCombInvYieldSys2760GeVALHC11h_2050->RemovePoint(0);

		graphCombInvYieldTot2760GeVALHC11h_0010->RemovePoint(0);
		graphCombInvYieldTot2760GeVALHC11h_2050->RemovePoint(0);


	} else  if(meson.CompareTo("Eta")==0){
		// Declaration & calculation of combined spectrum

		graphCombInvYieldTot2760GeVALHC11h_0010       = CombinePtPointsSpectraFullCorrMat(  statErrorCollectionLHC11h_0010,       sysErrorCollectionLHC11h_0010,        
																							xPtLimitsEta, /*17*/13,
																							offSetsEta, offSetsEtaSys,
																							graphCombInvYieldStat2760GeVALHC11h_0010, graphCombInvYieldSys2760GeVALHC11h_0010,
																							fileNameOutputWeightingALHC11h_0010,1 );

		graphCombInvYieldTot2760GeVALHC11h_2050       = CombinePtPointsSpectraFullCorrMat(  statErrorCollectionLHC11h_2050,       sysErrorCollectionLHC11h_2050,        
																							xPtLimitsEta, /*17*/13,
																							offSetsEta, offSetsEtaSys,
																							graphCombInvYieldStat2760GeVALHC11h_2050, graphCombInvYieldSys2760GeVALHC11h_2050,
																							fileNameOutputWeightingALHC11h_2050,1 );
		graphCombInvYieldStat2760GeVALHC11h_0010->RemovePoint(0);
		graphCombInvYieldStat2760GeVALHC11h_2050->RemovePoint(0);
		graphCombInvYieldStat2760GeVALHC11h_0010->RemovePoint(0);
		graphCombInvYieldStat2760GeVALHC11h_2050->RemovePoint(0);

		graphCombInvYieldSys2760GeVALHC11h_0010->RemovePoint(0);
		graphCombInvYieldSys2760GeVALHC11h_0010->RemovePoint(0);
		graphCombInvYieldSys2760GeVALHC11h_2050->RemovePoint(0);
		graphCombInvYieldSys2760GeVALHC11h_2050->RemovePoint(0);

		graphCombInvYieldTot2760GeVALHC11h_0010->RemovePoint(0);
		graphCombInvYieldTot2760GeVALHC11h_2050->RemovePoint(0);
		graphCombInvYieldTot2760GeVALHC11h_0010->RemovePoint(0);
		graphCombInvYieldTot2760GeVALHC11h_2050->RemovePoint(0);

	}

	// Reading weights from output file for plotting
	ifstream fileWeightsReadALHC11h_0010;
	fileWeightsReadALHC11h_0010.open(fileNameOutputWeightingALHC11h_0010,ios_base::in);
	ifstream fileWeightsReadALHC11h_2050;
	fileWeightsReadALHC11h_2050.open(fileNameOutputWeightingALHC11h_2050,ios_base::in);
	cout << "reading" << fileNameOutputWeightingALHC11h_0010 << " and " << fileNameOutputWeightingALHC11h_2050 << endl;
	Double_t xValuesReadALHC11h_0010[50];
	Double_t weightsReadALHC11h_0010[11][50];
	Int_t availableMeasALHC11h_0010[11] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
	Double_t xValuesReadALHC11h_2050[50];
	Double_t weightsReadALHC11h_2050[11][50];
	Int_t availableMeasALHC11h_2050[11] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

	Int_t nMeasSetALHC11h                         = 2;
	Int_t nPtBinsReadALHC11h              = 0;
	while(!fileWeightsReadALHC11h_0010.eof() && nPtBinsReadALHC11h < 50){
			TString garbage = "";
			if (nPtBinsReadALHC11h == 0){
					fileWeightsReadALHC11h_0010 >> garbage ;//>> availableMeas[0] >> availableMeas[1] >> availableMeas[2] >> availableMeas[3];
					fileWeightsReadALHC11h_2050 >> garbage ;//>> availableMeas[0] >> availableMeas[1] >> availableMeas[2] >> availableMeas[3];
					for (Int_t i = 0; i < nMeasSetALHC11h; i++){
							fileWeightsReadALHC11h_0010 >> availableMeasALHC11h_0010[i] ;
							fileWeightsReadALHC11h_2050 >> availableMeasALHC11h_2050[i] ;
					}       
					cout << "read following measurements: "; 
					for (Int_t i = 0; i < 11; i++){
							cout << availableMeasALHC11h_0010[i] << "\t" ;
							cout << availableMeasALHC11h_2050[i] << "\t" ;
					}       
					cout << endl;
			} else {
					fileWeightsReadALHC11h_0010 >> xValuesReadALHC11h_0010[nPtBinsReadALHC11h-1];
					fileWeightsReadALHC11h_2050 >> xValuesReadALHC11h_2050[nPtBinsReadALHC11h-1];
					for (Int_t i = 0; i < nMeasSetALHC11h; i++){
							fileWeightsReadALHC11h_0010 >> weightsReadALHC11h_0010[availableMeasALHC11h_0010[i]][nPtBinsReadALHC11h-1] ;
							fileWeightsReadALHC11h_2050 >> weightsReadALHC11h_2050[availableMeasALHC11h_2050[i]][nPtBinsReadALHC11h-1] ;
					}       
					cout << "read: "<<  nPtBinsReadALHC11h << " xValuesReadALHC11h_0010 "<< xValuesReadALHC11h_0010[nPtBinsReadALHC11h-1] << "\t" ;
					cout << "read: "<<  nPtBinsReadALHC11h << " xValuesReadALHC11h_2050 "<< xValuesReadALHC11h_2050[nPtBinsReadALHC11h-1] << "\t" ;
					for (Int_t i = 0; i < nMeasSetALHC11h; i++){
							cout << weightsReadALHC11h_0010[availableMeasALHC11h_0010[i]][nPtBinsReadALHC11h-1] << "\t";
							cout << weightsReadALHC11h_2050[availableMeasALHC11h_2050[i]][nPtBinsReadALHC11h-1] << "\t";
					}
					cout << endl;
			}
			nPtBinsReadALHC11h++;
	}
	nPtBinsReadALHC11h = nPtBinsReadALHC11h-2 ;
	
	fileWeightsReadALHC11h_0010.close();
	fileWeightsReadALHC11h_2050.close();

	for (Int_t i = 0; i < nMeasSetALHC11h; i++){
			graphWeightsALHC11h_0010[availableMeasALHC11h_0010[i]] = new TGraph(nPtBinsReadALHC11h,xValuesReadALHC11h_0010,weightsReadALHC11h_0010[availableMeasALHC11h_0010[i]]);
			graphWeightsALHC11h_2050[availableMeasALHC11h_2050[i]] = new TGraph(nPtBinsReadALHC11h,xValuesReadALHC11h_2050,weightsReadALHC11h_2050[availableMeasALHC11h_2050[i]]);
	graphWeightsALHC11h_0010[availableMeasALHC11h_0010[i]]->Print();
	graphWeightsALHC11h_2050[availableMeasALHC11h_2050[i]]->Print();

			Int_t bin = 0;
			for (Int_t n = 0; n< nPtBinsReadALHC11h; n++){
					if (graphWeightsALHC11h_0010[availableMeasALHC11h_0010[i]]->GetY()[bin] == 0 && graphWeightsALHC11h_2050[availableMeasALHC11h_2050[i]]->GetY()[bin] == 0){
							graphWeightsALHC11h_0010[availableMeasALHC11h_0010[i]]->RemovePoint(bin);
							graphWeightsALHC11h_2050[availableMeasALHC11h_2050[i]]->RemovePoint(bin);
					} else bin++;
			}       
	}       
	
	
	//**********************************************************************************************************************//
	//******************************************* Calculation of spectrum PCM - EMCal LHC11h *******************************//
	//**********************************************************************************************************************//
							
	for (Int_t i = 0; i< 11; i++){
		graphWeightsLHC11h_0010[i] = NULL;
		graphWeightsLHC11h_2050[i] = NULL;
	}	

	// Declaration & calculation of combined spectrum
	TString fileNameOutputWeightingLHC11h_0010					= Form("%s/0010LHC11h_WeightingEMCal%s.dat",outputDir.Data(),meson.Data());
	TString fileNameOutputWeightingLHC11h_2050					= Form("%s/2050LHC11h_WeightingEMCal%s.dat",outputDir.Data(),meson.Data());
	
	if(meson.CompareTo("Pi0")==0){
		cout << "******************************************* calculating the 0-10% combined spectra *******************************************\n\n\n" << endl;
		graphCombInvYieldTot2760GeVLHC11h_0010  = CombinePtPointsSpectraFullCorrMat( statErrorCollectionLHC11h_0010,	sysErrorCollectionLHC11h_0010, 	
																						xPtLimitsPi0, 23,
																						offSetsPi0, offSetsPi0Sys,
																						graphCombInvYieldStat2760GeVLHC11h_0010, graphCombInvYieldSys2760GeVLHC11h_0010,
																						fileNameOutputWeightingLHC11h_0010, 1);
		cout << "******************************************* calculating the 20-50% combined spectra *******************************************\n\n\n" << endl;
		graphCombInvYieldTot2760GeVLHC11h_2050 = CombinePtPointsSpectraFullCorrMat(  statErrorCollectionLHC11h_2050,	sysErrorCollectionLHC11h_2050, 	
																						xPtLimitsPi0, 23,
																						offSetsPi0, offSetsPi0Sys,
																						graphCombInvYieldStat2760GeVLHC11h_2050, graphCombInvYieldSys2760GeVLHC11h_2050,
																						fileNameOutputWeightingLHC11h_2050, 1);		
		graphCombInvYieldStat2760GeVLHC11h_0010->RemovePoint(0);
		graphCombInvYieldStat2760GeVLHC11h_2050->RemovePoint(0);

		graphCombInvYieldSys2760GeVLHC11h_0010->RemovePoint(0);
		graphCombInvYieldSys2760GeVLHC11h_2050->RemovePoint(0);

		graphCombInvYieldTot2760GeVLHC11h_0010->RemovePoint(0);
		graphCombInvYieldTot2760GeVLHC11h_2050->RemovePoint(0);

	} else 	if(meson.CompareTo("Eta")==0){
		cout << "******************************************* calculating the 0-10% combined spectra *******************************************\n\n\n" << endl;
		graphCombInvYieldTot2760GeVLHC11h_0010  = CombinePtPointsSpectraFullCorrMat( statErrorCollectionLHC11h_0010,	sysErrorCollectionLHC11h_0010, 	
																						xPtLimitsEta, /*17*/13,
																						offSetsEta, offSetsEtaSys,
																						graphCombInvYieldStat2760GeVLHC11h_0010, graphCombInvYieldSys2760GeVLHC11h_0010,
																						fileNameOutputWeightingLHC11h_0010,1 );
		cout << "******************************************* calculating the 20-50% combined spectra *******************************************\n\n\n" << endl;
		graphCombInvYieldTot2760GeVLHC11h_2050 = CombinePtPointsSpectraFullCorrMat(  statErrorCollectionLHC11h_2050,	sysErrorCollectionLHC11h_2050, 	
																						xPtLimitsEta, /*17*/13,
																						offSetsEta, offSetsEtaSys,
																						graphCombInvYieldStat2760GeVLHC11h_2050, graphCombInvYieldSys2760GeVLHC11h_2050,
																						fileNameOutputWeightingLHC11h_2050,1 );
		
		graphCombInvYieldStat2760GeVLHC11h_0010->RemovePoint(0);
		graphCombInvYieldStat2760GeVLHC11h_2050->RemovePoint(0);
		graphCombInvYieldStat2760GeVLHC11h_0010->RemovePoint(0);
		graphCombInvYieldStat2760GeVLHC11h_2050->RemovePoint(0);

		graphCombInvYieldSys2760GeVLHC11h_0010->RemovePoint(0);
		graphCombInvYieldSys2760GeVLHC11h_0010->RemovePoint(0);
		graphCombInvYieldSys2760GeVLHC11h_2050->RemovePoint(0);
		graphCombInvYieldSys2760GeVLHC11h_2050->RemovePoint(0);

		graphCombInvYieldTot2760GeVLHC11h_0010->RemovePoint(0);
		graphCombInvYieldTot2760GeVLHC11h_0010->RemovePoint(0);
		graphCombInvYieldTot2760GeVLHC11h_2050->RemovePoint(0);
		graphCombInvYieldTot2760GeVLHC11h_2050->RemovePoint(0);

	}

	// Reading weights from output file for plotting
	ifstream fileWeightsReadLHC11h_0010;
	fileWeightsReadLHC11h_0010.open(fileNameOutputWeightingLHC11h_0010,ios_base::in);
	ifstream fileWeightsReadLHC11h_2050;
	fileWeightsReadLHC11h_2050.open(fileNameOutputWeightingLHC11h_2050,ios_base::in);

	cout << "reading " << fileNameOutputWeightingLHC11h_0010 << " and " << fileNameOutputWeightingLHC11h_2050 << endl;
	Double_t xValuesReadLHC11h_0010[50];
	Double_t xValuesReadLHC11h_2050[50];
	Double_t weightsReadLHC11h_0010[11][50];
	Double_t weightsReadLHC11h_2050[11][50];
	Int_t availableMeasLHC11h_0010[11] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
	Int_t availableMeasLHC11h_2050[11] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
		
	Int_t nMeasSetLHC11h 			= 2;
	Int_t nPtBinsReadLHC11h 		= 0;
	while(!fileWeightsReadLHC11h_0010.eof() && nPtBinsReadLHC11h < 50){
		TString garbage = "";
		if (nPtBinsReadLHC11h == 0){
			fileWeightsReadLHC11h_0010 >> garbage;
			fileWeightsReadLHC11h_2050 >> garbage;
			for (Int_t i = 0; i < nMeasSetLHC11h; i++){
				fileWeightsReadLHC11h_0010 >> availableMeasLHC11h_0010[i];
				fileWeightsReadLHC11h_2050 >> availableMeasLHC11h_2050[i] ;
			}	
			cout << "read following measurements: "; 
			for (Int_t i = 0; i < 11; i++){
				cout << availableMeasLHC11h_0010[i] << "\t" ;
				cout << availableMeasLHC11h_2050[i] << "\t" ;
			}	
			cout << endl;
		} else {
			fileWeightsReadLHC11h_0010 >> xValuesReadLHC11h_0010[nPtBinsReadLHC11h-1];
			fileWeightsReadLHC11h_2050 >> xValuesReadLHC11h_2050[nPtBinsReadLHC11h-1];
			for (Int_t i = 0; i < nMeasSetLHC11h; i++){
				fileWeightsReadLHC11h_0010 >> weightsReadLHC11h_0010[availableMeasLHC11h_0010[i]][nPtBinsReadLHC11h-1] ;
				fileWeightsReadLHC11h_2050 >> weightsReadLHC11h_2050[availableMeasLHC11h_2050[i]][nPtBinsReadLHC11h-1] ;
			}	
			cout << "read: "<<  nPtBinsReadLHC11h << "\t"<< xValuesReadLHC11h_0010[nPtBinsReadLHC11h-1] << "\t" ;
			cout << "read: "<<  nPtBinsReadLHC11h << "\t"<< xValuesReadLHC11h_2050[nPtBinsReadLHC11h-1] << "\t" ;
			for (Int_t i = 0; i < nMeasSetLHC11h; i++){
				cout << weightsReadLHC11h_0010[availableMeasLHC11h_0010[i]][nPtBinsReadLHC11h-1] << "\t";
				cout << weightsReadLHC11h_2050[availableMeasLHC11h_2050[i]][nPtBinsReadLHC11h-1] << "\t";
			}
			cout << endl;
		}
		nPtBinsReadLHC11h++;
	}
	nPtBinsReadLHC11h = nPtBinsReadLHC11h-2 ;
	fileWeightsReadLHC11h_0010.close();
	fileWeightsReadLHC11h_2050.close();

	for (Int_t i = 0; i < nMeasSetLHC11h; i++){
		graphWeightsLHC11h_0010[availableMeasLHC11h_0010[i]] = new TGraph(nPtBinsReadLHC11h,xValuesReadLHC11h_0010,weightsReadLHC11h_0010[availableMeasLHC11h_0010[i]]);
		graphWeightsLHC11h_2050[availableMeasLHC11h_2050[i]] = new TGraph(nPtBinsReadLHC11h,xValuesReadLHC11h_2050,weightsReadLHC11h_2050[availableMeasLHC11h_2050[i]]);
		Int_t bin = 0;
		for (Int_t n = 0; n< nPtBinsReadLHC11h; n++){
			if (graphWeightsLHC11h_0010[availableMeasLHC11h_0010[i]]->GetY()[bin] == 0 && graphWeightsLHC11h_2050[availableMeasLHC11h_2050[i]]->GetY()[bin] == 0){
					graphWeightsLHC11h_0010[availableMeasLHC11h_0010[i]]->RemovePoint(bin);
					graphWeightsLHC11h_2050[availableMeasLHC11h_2050[i]]->RemovePoint(bin);
			} else bin++;

		}	
	}	

	//**********************************************************************************************************************//
	//******************************************* Plotting weights method only EMC *****************************************//
	//**********************************************************************************************************************//
	
	canvasWeights->cd();
	histo2DWeights->Draw("copy");

	TLegend* legendAccWeightsLHC11h = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.06*nMeasSetLHC11h*1.35), 32);
	for (Int_t i = 0; i < nMeasSetLHC11h; i++){
		DrawGammaSetMarkerTGraph(graphWeightsLHC11h_0010[availableMeasLHC11h_0010[i]],
								markerStyleDet[availableMeasLHC11h_0010[i]], 
								markerSizeDet[availableMeasLHC11h_0010[i]]*0.5, 
								colorDet[availableMeasLHC11h_0010[i]] , 
								colorDet[availableMeasLHC11h_0010[i]]);
		graphWeightsLHC11h_0010[availableMeasLHC11h_0010[i]]->Draw("p,same,e0");
		legendAccWeightsLHC11h->AddEntry(graphWeightsLHC11h_0010[availableMeasLHC11h_0010[i]],nameMeasGlobal[availableMeasLHC11h_0010[i]],"p");
		
		DrawGammaSetMarkerTGraph(graphWeightsLHC11h_2050[availableMeasLHC11h_2050[i]],
								markerStyleDet[availableMeasLHC11h_2050[i]], 
								markerSizeDet[availableMeasLHC11h_2050[i]]*0.5, 
								colorDet[availableMeasLHC11h_2050[i]]+1, 
								colorDet[availableMeasLHC11h_2050[i]]+1);
		graphWeightsLHC11h_2050[availableMeasLHC11h_2050[i]]->Draw("p,same,e0");
		legendAccWeightsLHC11h->AddEntry(graphWeightsLHC11h_2050[availableMeasLHC11h_2050[i]],nameMeasGlobal[availableMeasLHC11h_2050[i]],"p");
	}	
	legendAccWeightsLHC11h->Draw();
	labelWeightsEnergy->Draw();
	labelWeightsPi0->Draw();

//	DrawGammaLines(0.23, 70. , 0.8, 0.8,0.1, kGray, 3);
	DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
	DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
	DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
	DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);
		
	canvasWeights->SaveAs(Form("%s/%s_WeightsPCM-EMCal_LHC11h.%s",outputDir.Data(),meson.Data(),suffix.Data()));

		
	TH1D* statErrorRelCollectionLHC11h_0010[11];
	TH1D* statErrorRelCollectionLHC11h_2050[11];

	for (Int_t i = 0; i< 11; i++){
		statErrorRelCollectionLHC11h_0010[i] = NULL;
		statErrorRelCollectionLHC11h_2050[i] = NULL;
	}	
	for (Int_t i = 0; i < 11; i++){
		if (statErrorCollectionLHC11h_0010[i]) statErrorRelCollectionLHC11h_0010[i] = CalculateRelErrUpTH1D( statErrorCollectionLHC11h_0010[i], Form("relativeStatErrorLHC11h_%s_0010", nameMeasGlobal[i].Data()));
		
		if (statErrorCollectionLHC11h_2050[i]) statErrorRelCollectionLHC11h_2050[i] = CalculateRelErrUpTH1D( statErrorCollectionLHC11h_2050[i], Form("relativeStatErrorLHC11h_%s_2050", nameMeasGlobal[i].Data()));
	}
		
	TGraphAsymmErrors* sysErrorRelCollectionLHC11h_0010[11];
	TGraphAsymmErrors* sysErrorRelCollectionLHC11h_2050[11];
	for (Int_t i = 0; i< 11; i++){
		sysErrorRelCollectionLHC11h_0010[i] = NULL;
		sysErrorRelCollectionLHC11h_2050[i] = NULL;
	}	
	for (Int_t i = 0; i < 11; i++){
		if (sysErrorCollectionLHC11h_0010[i]) sysErrorRelCollectionLHC11h_0010[i] = CalculateRelErrUpAsymmGraph( sysErrorCollectionLHC11h_0010[i], Form("relativeSysErrorLHC11h_%s_0010", nameMeasGlobal[i].Data()));
		
		if (sysErrorCollectionLHC11h_2050[i]) sysErrorRelCollectionLHC11h_2050[i] = CalculateRelErrUpAsymmGraph( sysErrorCollectionLHC11h_2050[i], Form("relativeSysErrorLHC11h_%s_2050", nameMeasGlobal[i].Data()));				
	}
		


	//*********************************************************************************************************************//
	//******************************* Visualize relative errors - PCM - EMCal 2011 ****************************************//
	//*********************************************************************************************************************//
	
	TCanvas* canvasRelSysErrLHC11h = new TCanvas("canvasRelSysErrLHC11h","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRelSysErrLHC11h, 0.08, 0.02, 0.035, 0.09);
	canvasRelSysErrLHC11h->SetLogx();

		TH2F * histo2DRelSysErrLHC11h;
		histo2DRelSysErrLHC11h = new TH2F("histo2DRelSysErrLHC11h","histo2DRelSysErrLHC11h",11000,0.23,70.,1000,0,80.5);
		SetStyleHistoTH2ForGraphs(histo2DRelSysErrLHC11h, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
		histo2DRelSysErrLHC11h->GetXaxis()->SetMoreLogLabels();
		histo2DRelSysErrLHC11h->GetXaxis()->SetLabelOffset(-0.01);
		histo2DRelSysErrLHC11h->Draw("copy");	
		
		TLegend* legendRelSysErrLHC11h = GetAndSetLegend2(0.62, 0.94-(0.06*nMeasSetLHC11h*1.35), 0.95, 0.94, 32);
		for (Int_t i = 0; i < nMeasSetLHC11h; i++){
			DrawGammaSetMarkerTGraph(sysErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[i]], markerStyleDet[availableMeasLHC11h_0010[i]], 
										markerSizeDet[availableMeasLHC11h_0010[i]]*0.5, 	
										colorDet[availableMeasLHC11h_0010[i]], colorDet[availableMeasLHC11h_0010[i]]);
			sysErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[i]]->Draw("p,same,e0");
			legendRelSysErrLHC11h->AddEntry(sysErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[i]],nameMeasGlobal[availableMeasLHC11h_0010[i]],"p");

			DrawGammaSetMarkerTGraph(sysErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[i]], markerStyleDet[availableMeasLHC11h_2050[i]],
										markerSizeDet[availableMeasLHC11h_2050[i]]*0.5, 
										colorDet[availableMeasLHC11h_2050[i]]+1, colorDet[availableMeasLHC11h_2050[i]]+1);	
			sysErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[i]]->Draw("p,same,e0");
			legendRelSysErrLHC11h->AddEntry(sysErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[i]],nameMeasGlobal[availableMeasLHC11h_2050[i]],"p");
		}
		legendRelSysErrLHC11h->Draw();


		TLatex *labelRelSysErrEnergy = new TLatex(0.7,0.89,collisionSystem2760GeV.Data());
		SetStyleTLatex( labelRelSysErrEnergy, 0.85*textSizeLabelsPixel,4);
		labelRelSysErrEnergy->SetTextFont(43);
		labelRelSysErrEnergy->Draw();
			
		TLatex *labelRelSysErrPi0;
		if(meson.CompareTo("Pi0")==0){
			labelRelSysErrPi0= new TLatex(0.75,0.85,"#pi^{0} #rightarrow #gamma#gamma");
		} else if(meson.CompareTo("Eta")==0){
			labelRelSysErrPi0= new TLatex(0.75,0.85,"#eta #rightarrow #gamma#gamma");
		}
		SetStyleTLatex( labelRelSysErrPi0, 0.85*textSizeLabelsPixel,4);
		labelRelSysErrPi0->SetTextFont(43);
		labelRelSysErrPi0->Draw();
	
	canvasRelSysErrLHC11h->SaveAs(Form("%s/%s_RelSysErrLHC11h.%s",outputDir.Data(),meson.Data(),suffix.Data()));
		
		histo2DRelSysErrLHC11h->GetYaxis()->SetRangeUser(0,30.5);
		histo2DRelSysErrLHC11h->Draw("copy");	

		for (Int_t i = 0; i < nMeasSetLHC11h; i++){
			sysErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[i]]->Draw("p,same,e0");
			sysErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[i]]->Draw("p,same,e0");								
			
		}	
		legendRelSysErrLHC11h->Draw();
		labelRelSysErrEnergy->Draw();
		labelRelSysErrPi0->Draw();

	canvasRelSysErrLHC11h->SaveAs(Form("%s/%s_RelSysErrZoomedLHC11h.%s",outputDir.Data(),meson.Data(), suffix.Data()));
		
        histo2DRelSysErrLHC11h->GetYaxis()->SetRangeUser(0.,60.5);
        histo2DRelSysErrLHC11h->GetXaxis()->SetRangeUser(0.3,25.);
        if(meson.CompareTo("Eta")==0){
          histo2DRelSysErrLHC11h->GetYaxis()->SetRangeUser(0.,60.5);
          histo2DRelSysErrLHC11h->GetXaxis()->SetRangeUser(0.8,15.);
        }          
        histo2DRelSysErrLHC11h->Draw("copy"); 

          DrawGammaSetMarkerTGraph(sysErrorRelCollectionLHC11h_0010[0], markerStyleDet[0], markerSizeDet[0]*0.5, colorCombo0010 ,colorCombo0010);
          DrawGammaSetMarkerTGraph(sysErrorRelCollectionLHC11h_2050[0], markerStyleDet[0], markerSizeDet[0]*0.5, colorCombo2050, colorCombo2050);

          sysErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[0]]->Draw("p,same,e0");
          sysErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[0]]->Draw("p,same,e0");

          TLegend* legendRelSysErrLHC11hforAN = GetAndSetLegend2(0.14, 0.94-(0.06*nMeasSetLHC11h*1.35), 0.45, 0.94, 32);
          legendRelSysErrLHC11hforAN->AddEntry(sysErrorRelCollectionLHC11h_0010[0],"PCM -   0-10%","p");
          legendRelSysErrLHC11hforAN->AddEntry(sysErrorRelCollectionLHC11h_2050[0],"PCM - 20-50%","p");
          legendRelSysErrLHC11hforAN->Draw();
          labelRelSysErrEnergy->Draw();
          labelRelSysErrPi0->Draw();
            
    canvasRelSysErrLHC11h->SaveAs(Form("%s/%s_RelSysErrZoomedLHC11hforAN.%s",outputDir.Data(),meson.Data(),suffix.Data()));
        

	//*********************************************************************************************************************//
	//******************************* Visualize relative errors PCM - EMCal 2011 ******************************************//
	//*********************************************************************************************************************//
		
	TCanvas* canvasRelStatErrLHC11h = new TCanvas("canvasRelStatErrLHC11h","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasRelStatErrLHC11h, 0.08, 0.02, 0.035, 0.09);
		canvasRelStatErrLHC11h->SetLogx();
		
		TH2F * histo2DRelStatErrLHC11h;
		histo2DRelStatErrLHC11h = new TH2F("histo2DRelStatErrLHC11h","histo2DRelStatErrLHC11h",11000,0.23,70.,1000,0,80.5);
		SetStyleHistoTH2ForGraphs(histo2DRelStatErrLHC11h, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
		histo2DRelStatErrLHC11h->GetXaxis()->SetMoreLogLabels();
		histo2DRelStatErrLHC11h->GetXaxis()->SetLabelOffset(-0.01);
		histo2DRelStatErrLHC11h->Draw("copy");	

		TLegend* legendRelStatErrLHC11h = GetAndSetLegend2(0.14, 0.94-(0.06*nMeasSetLHC11h*1.35), 0.45, 0.94, 32);
		for (Int_t i = 0; i < nMeasSetLHC11h; i++){
			DrawGammaSetMarker(statErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[i]], 
							markerStyleDet[availableMeasLHC11h_0010[i]], markerSizeDet[availableMeasLHC11h_0010[i]]*0.5, 
							colorDet[availableMeasLHC11h_0010[i]] , colorDet[availableMeasLHC11h_0010[i]]);
			statErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[i]]->Draw("p,same,e0");
			legendRelStatErrLHC11h->AddEntry(statErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[i]],nameMeasGlobal[availableMeasLHC11h_0010[i]],"p");
			
			DrawGammaSetMarker(statErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[i]], 
							markerStyleDet[availableMeasLHC11h_2050[i]], markerSizeDet[availableMeasLHC11h_2050[i]]*0.5, 
							colorDet[availableMeasLHC11h_2050[i]]+1, colorDet[availableMeasLHC11h_2050[i]]+1);
			statErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[i]]->Draw("p,same,e0");
			legendRelStatErrLHC11h->AddEntry(statErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[i]],nameMeasGlobal[availableMeasLHC11h_2050[i]],"p");		
		}	
		legendRelStatErrLHC11h->Draw();

		TLatex *labelRelStatErrEnergy = new TLatex(0.7,0.89,collisionSystem2760GeV.Data());
		SetStyleTLatex( labelRelStatErrEnergy, 0.85*textSizeLabelsPixel,4);
		labelRelStatErrEnergy->SetTextFont(43);
		labelRelStatErrEnergy->Draw();
		TLatex *labelRelStatErrPi0;
		if(meson.CompareTo("Pi0")==0){
			labelRelStatErrPi0= new TLatex(0.75,0.85,"#pi^{0} #rightarrow #gamma#gamma");
		} else if(meson.CompareTo("Eta")==0){
			labelRelStatErrPi0= new TLatex(0.75,0.85,"#eta #rightarrow #gamma#gamma");
		}
		SetStyleTLatex( labelRelStatErrPi0, 0.85*textSizeLabelsPixel,4);
		labelRelStatErrPi0->SetTextFont(43);
		labelRelStatErrPi0->Draw();
			
	canvasRelStatErrLHC11h->SaveAs(Form("%s/%s_RelStatErrLHC11h.%s",outputDir.Data(),meson.Data(),suffix.Data()));

		histo2DRelStatErrLHC11h->GetYaxis()->SetRangeUser(0,30.5);
		histo2DRelStatErrLHC11h->Draw("copy");	
			for (Int_t i = 0; i < nMeasSetLHC11h; i++){
				statErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[i]]->Draw("p,same,e0");
				statErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[i]]->Draw("p,same,e0");
			}	
			legendRelStatErrLHC11h->Draw();

			labelRelStatErrEnergy->Draw();
			labelRelStatErrPi0->Draw();
			
	canvasRelStatErrLHC11h->SaveAs(Form("%s/%s_RelStatErrZoomedLHC11h.%s",outputDir.Data(),meson.Data(),suffix.Data()));
		
        histo2DRelStatErrLHC11h->GetYaxis()->SetRangeUser(0.,60.5);
        histo2DRelStatErrLHC11h->GetXaxis()->SetRangeUser(0.3,25.);
        if(meson.CompareTo("Eta")==0){
          histo2DRelStatErrLHC11h->GetYaxis()->SetRangeUser(0.,60.5);
          histo2DRelStatErrLHC11h->GetXaxis()->SetRangeUser(0.8,15.);
        }          
        histo2DRelStatErrLHC11h->Draw("copy"); 

          DrawGammaSetMarker(statErrorRelCollectionLHC11h_0010[0], markerStyleDet[0], markerSizeDet[0]*0.5, colorCombo0010 ,colorCombo0010);
          DrawGammaSetMarker(statErrorRelCollectionLHC11h_2050[0], markerStyleDet[0], markerSizeDet[0]*0.5, colorCombo2050, colorCombo2050);

          statErrorRelCollectionLHC11h_0010[availableMeasLHC11h_0010[0]]->Draw("p,same,e0");
          statErrorRelCollectionLHC11h_2050[availableMeasLHC11h_2050[0]]->Draw("p,same,e0");

          TLegend* legendRelStatErrLHC11hforAN = GetAndSetLegend2(0.14, 0.94-(0.06*nMeasSetLHC11h*1.35), 0.45, 0.94, 32);
          legendRelStatErrLHC11hforAN->AddEntry(statErrorRelCollectionLHC11h_0010[0],"PCM -   0-10%","p");
          legendRelStatErrLHC11hforAN->AddEntry(statErrorRelCollectionLHC11h_2050[0],"PCM - 20-50%","p");
          legendRelStatErrLHC11hforAN->Draw();
          labelRelStatErrEnergy->Draw();
          labelRelStatErrPi0->Draw();
            
    canvasRelStatErrLHC11h->SaveAs(Form("%s/%s_RelStatErrZoomedLHC11hforAN.%s",outputDir.Data(),meson.Data(),suffix.Data()));
		
	//*********************************************************************************************************************//
	//**************************************** Visualize relative total errors ********************************************//
	//*********************************************************************************************************************//
	
	TGraphAsymmErrors* graphCombInvYieldRelStat2760GeVALHC11h_0010        = CalculateRelErrUpAsymmGraph( graphCombInvYieldStat2760GeVALHC11h_0010, "relativeStatErrorLHC11h_MethodA_0010");
	TGraphAsymmErrors* graphCombInvYieldRelSys2760GeVALHC11h_0010         = CalculateRelErrUpAsymmGraph( graphCombInvYieldSys2760GeVALHC11h_0010, "relativeSysErrorLHC11h_MethodA_0010");
	TGraphAsymmErrors* graphCombInvYieldRelTot2760GeVALHC11h_0010         = CalculateRelErrUpAsymmGraph( graphCombInvYieldTot2760GeVALHC11h_0010, "relativeTotalErrorLHC11h_MethodA_0010");
	
	TGraphAsymmErrors* graphCombInvYieldRelTot2760GeVLHC11h_0010 = CalculateRelErrUpAsymmGraph( graphCombInvYieldTot2760GeVLHC11h_0010, "relativeTotalError_OEMC_0010");

	TGraphAsymmErrors* graphCombInvYieldRelStat2760GeVALHC11h_2050        = CalculateRelErrUpAsymmGraph( graphCombInvYieldStat2760GeVALHC11h_2050, "relativeStatErrorLHC11h_MethodA_2050");
	TGraphAsymmErrors* graphCombInvYieldRelSys2760GeVALHC11h_2050         = CalculateRelErrUpAsymmGraph( graphCombInvYieldSys2760GeVALHC11h_2050, "relativeSysErrorLHC11h_MethodA_2050");
	TGraphAsymmErrors* graphCombInvYieldRelTot2760GeVALHC11h_2050         = CalculateRelErrUpAsymmGraph( graphCombInvYieldTot2760GeVALHC11h_2050, "relativeTotalErrorLHC11h_MethodA_2050");
 
	
	//**********************************************************************************************************************//
	//************************************* Calculating bin shifted spectra & fitting **************************************//
	//**********************************************************************************************************************//
		
	// Cloning spectra
	TGraphAsymmErrors* graphCombInvYieldTot2760GeVALHC11hUnShifted_0010 		= (TGraphAsymmErrors*)graphCombInvYieldTot2760GeVALHC11h_0010->Clone("Unshifted_0010"); 
	TGraphAsymmErrors* graphCombInvYieldStat2760GeVALHC11hUnShifted_0010 		= (TGraphAsymmErrors*)graphCombInvYieldStat2760GeVALHC11h_0010->Clone("UnshiftedStat_0010"); 
	TGraphAsymmErrors* graphCombInvYieldSys2760GeVALHC11hUnShifted_0010 		= (TGraphAsymmErrors*)graphCombInvYieldSys2760GeVALHC11h_0010->Clone("UnshiftedSys_0010"); 

	TGraphAsymmErrors* graphCombInvYieldTot2760GeVALHC11hUnShifted_2050 		= (TGraphAsymmErrors*)graphCombInvYieldTot2760GeVALHC11h_2050->Clone("Unshifted_2050"); 
	TGraphAsymmErrors* graphCombInvYieldStat2760GeVALHC11hUnShifted_2050 		= (TGraphAsymmErrors*)graphCombInvYieldStat2760GeVALHC11h_2050->Clone("UnshiftedStat_2050"); 
	TGraphAsymmErrors* graphCombInvYieldSys2760GeVALHC11hUnShifted_2050  		= (TGraphAsymmErrors*)graphCombInvYieldSys2760GeVALHC11h_2050->Clone("UnshiftedSys_2050"); 

	if(meson.CompareTo("Pi0")==0){
		graphPCMInvYieldStatPbPb2760GeVUnShifted_0010 		= (TGraphAsymmErrors*)graphPCMPi0InvYieldStatPbPb2760GeV_0010->Clone("UnshiftedStatPCMPi0_0010"); 
		graphPCMInvYieldSysPbPb2760GeVUnShifted_0010 		= (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_0010->Clone("UnshiftedSysPCMPi0_0010"); 

		graphEMCalInvYieldStatPbPb2760GeVUnshifted_0010 	= (TGraphAsymmErrors*)graphEMCalPi0InvYieldStatPbPb2760GeV_0010->Clone("UnshiftedStatEMCalPi0_0010"); 
		graphEMCalInvYieldSysPbPb2760GeVUnshifted_0010 	= (TGraphAsymmErrors*)graphEMCalPi0InvYieldSysPbPb2760GeV_0010->Clone("UnshiftedSysEMCalPi0_0010"); 
		
		graphPCMInvYieldStatPbPb2760GeVUnShifted_2050 		= (TGraphAsymmErrors*)graphPCMPi0InvYieldStatPbPb2760GeV_2050->Clone("UnshiftedStatPCMPi0_2050"); 
		graphPCMInvYieldSysPbPb2760GeVUnShifted_2050 		= (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_2050->Clone("UnshiftedSysPCMPi0_2050"); 

		graphEMCalInvYieldStatPbPb2760GeVUnshifted_2050 	= (TGraphAsymmErrors*)graphEMCalPi0InvYieldStatPbPb2760GeV_2050->Clone("UnshiftedStatEMCalPi0_2050"); 
		graphEMCalInvYieldSysPbPb2760GeVUnshifted_2050 	= (TGraphAsymmErrors*)graphEMCalPi0InvYieldSysPbPb2760GeV_2050->Clone("UnshiftedSysEMCalPi0_2050"); 
		
	} else if(meson.CompareTo("Eta")==0){
		graphPCMInvYieldStatPbPb2760GeVUnShifted_0010 		= (TGraphAsymmErrors*)graphPCMEtaInvYieldStatPbPb2760GeV_0010->Clone("UnshiftedStatPCMEta_0010"); 
		graphPCMInvYieldSysPbPb2760GeVUnShifted_0010 		= (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_0010->Clone("UnshiftedSysPCMEta_0010"); 

		graphEMCalInvYieldStatPbPb2760GeVUnshifted_0010 	= (TGraphAsymmErrors*)graphEMCalEtaInvYieldStatPbPb2760GeV_0010->Clone("UnshiftedStatEMCalEta_0010"); 
		graphEMCalInvYieldSysPbPb2760GeVUnshifted_0010 	= (TGraphAsymmErrors*)graphEMCalEtaInvYieldSysPbPb2760GeV_0010->Clone("UnshiftedSysEMCalEta_0010"); 		

		graphPCMInvYieldStatPbPb2760GeVUnShifted_2050 		= (TGraphAsymmErrors*)graphPCMEtaInvYieldStatPbPb2760GeV_2050->Clone("UnshiftedStatPCMEta_2050"); 
		graphPCMInvYieldSysPbPb2760GeVUnShifted_2050 		= (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_2050->Clone("UnshiftedSysPCMEta_2050"); 

		graphEMCalInvYieldStatPbPb2760GeVUnshifted_2050 	= (TGraphAsymmErrors*)graphEMCalEtaInvYieldStatPbPb2760GeV_2050->Clone("UnshiftedStatEMCalEta_2050"); 
		graphEMCalInvYieldSysPbPb2760GeVUnshifted_2050 	= (TGraphAsymmErrors*)graphEMCalEtaInvYieldSysPbPb2760GeV_2050->Clone("UnshiftedSysEMCalEta_2050"); 		
		
	}

	// Calculating binshifts
	TF1* fitBylinkinPbPb2760GeVPtLHC11h_0010;
	TF1* fitBylinkinPbPb2760GeVPtLHC11h_2050;
    TF1* fitBylinkinPbPb2760GeVPtLHC11hYshift_0010;
    TF1* fitBylinkinPbPb2760GeVPtLHC11hYshift_2050;
	if(meson.CompareTo("Pi0")==0){
        Double_t paramTCM_0010[5] = {graphCombInvYieldTot2760GeVALHC11hUnShifted_0010->GetY()[0],0.3,graphCombInvYieldTot2760GeVALHC11hUnShifted_0010->GetY()[0],0.3,8};
		fitBylinkinPbPb2760GeVPtLHC11h_0010 = FitObject("tcm","BylinkinFitPi00010","Pi0",graphCombInvYieldTot2760GeVALHC11hUnShifted_0010,0.4,20.,paramTCM_0010);
		graphCombInvYieldTot2760GeVALHC11hUnShifted_0010->Fit(fitBylinkinPbPb2760GeVPtLHC11h_0010,"NRMEX0+","",0.4,20.);

        Double_t paramTCM_2050[5] = {graphCombInvYieldTot2760GeVALHC11hUnShifted_2050->GetY()[0],0.3,graphCombInvYieldTot2760GeVALHC11hUnShifted_2050->GetY()[0],0.3,8};        
		fitBylinkinPbPb2760GeVPtLHC11h_2050 = FitObject("tcm","BylinkinFitPi02050","Pi0",graphCombInvYieldTot2760GeVALHC11hUnShifted_2050,0.4,20.,paramTCM_2050);
		graphCombInvYieldTot2760GeVALHC11hUnShifted_2050->Fit(fitBylinkinPbPb2760GeVPtLHC11h_2050,"NRMEX0+","",0.4,20.);

        cout << WriteParameterToFile(fitBylinkinPbPb2760GeVPtLHC11h_0010)<< endl << endl;
        cout << WriteParameterToFile(fitBylinkinPbPb2760GeVPtLHC11h_2050)<< endl << endl;
        
		if(bWCorrection.CompareTo("X")==0 ){
			TF1* fitBylinkinPbPb2760GeVPtMultLHC11h_0010 = FitObject("tcmpt","BylinkinFitPi0","Pi0",graphCombInvYieldTot2760GeVALHC11h_0010);	
			
			graphCombInvYieldTot2760GeVALHC11h_0010			= ApplyXshift(graphCombInvYieldTot2760GeVALHC11h_0010, fitBylinkinPbPb2760GeVPtMultLHC11h_0010,"Pi0");
			
			graphCombInvYieldStat2760GeVALHC11h_0010 		= ApplyXshiftIndividualSpectra (graphCombInvYieldTot2760GeVALHC11h_0010, 
																							graphCombInvYieldStat2760GeVALHC11h_0010, 
																							fitBylinkinPbPb2760GeVPtMultLHC11h_0010,
																							0, 21,"Pi0");

			graphCombInvYieldSys2760GeVALHC11h_0010 		= ApplyXshiftIndividualSpectra (graphCombInvYieldTot2760GeVALHC11h_0010, 
																							graphCombInvYieldSys2760GeVALHC11h_0010, 
																							fitBylinkinPbPb2760GeVPtMultLHC11h_0010, 
																							0, 21,"Pi0");
			
			graphPCMPi0InvYieldStatPbPb2760GeV_0010			= ApplyXshiftIndividualSpectra( graphCombInvYieldTot2760GeVALHC11h_0010,
																							graphPCMPi0InvYieldStatPbPb2760GeV_0010,
																							fitBylinkinPbPb2760GeVPtMultLHC11h_0010, 
																							0, 18,"Pi0");
			
			graphPCMPi0InvYieldSysPbPb2760GeV_0010			= ApplyXshiftIndividualSpectra( graphCombInvYieldTot2760GeVALHC11h_0010, 
																							graphPCMPi0InvYieldSysPbPb2760GeV_0010, 
																							fitBylinkinPbPb2760GeVPtMultLHC11h_0010, 
																							0, 18,"Pi0");
			
			graphEMCalPi0InvYieldStatPbPb2760GeV_0010 		= ApplyXshiftIndividualSpectra( graphCombInvYieldTot2760GeVALHC11h_0010, 
																							graphEMCalPi0InvYieldStatPbPb2760GeV_0010, 
																							fitBylinkinPbPb2760GeVPtMultLHC11h_0010,
																							14, 7,"Pi0");	
			
			graphEMCalPi0InvYieldSysPbPb2760GeV_0010 		= ApplyXshiftIndividualSpectra( graphCombInvYieldTot2760GeVALHC11h_0010, 
																							graphEMCalPi0InvYieldSysPbPb2760GeV_0010, 
																							fitBylinkinPbPb2760GeVPtMultLHC11h_0010,
																							14, 7,"Pi0");	

			TF1* fitBylinkinPbPb2760GeVPtMultLHC11h_2050 = FitObject("tcmpt","BylinkinFitPi0","Pi0",graphCombInvYieldTot2760GeVALHC11h_2050);		
			graphCombInvYieldTot2760GeVALHC11h_2050			= ApplyXshift(graphCombInvYieldTot2760GeVALHC11h_2050, fitBylinkinPbPb2760GeVPtMultLHC11h_2050,"Pi0");
			
			graphCombInvYieldStat2760GeVALHC11h_2050 		= ApplyXshiftIndividualSpectra (graphCombInvYieldTot2760GeVALHC11h_2050, 
																							graphCombInvYieldStat2760GeVALHC11h_2050, 
																							fitBylinkinPbPb2760GeVPtMultLHC11h_2050,
																							0, 21,"Pi0");
			
			graphCombInvYieldSys2760GeVALHC11h_2050 		= ApplyXshiftIndividualSpectra (graphCombInvYieldTot2760GeVALHC11h_2050, 
																							graphCombInvYieldSys2760GeVALHC11h_2050, 
																							fitBylinkinPbPb2760GeVPtMultLHC11h_2050, 
																							0, 21,"Pi0");
			
			graphPCMPi0InvYieldStatPbPb2760GeV_2050			= ApplyXshiftIndividualSpectra( graphCombInvYieldTot2760GeVALHC11h_2050,
																							graphPCMPi0InvYieldStatPbPb2760GeV_2050,
																							fitBylinkinPbPb2760GeVPtMultLHC11h_2050, 
																							0, 18,"Pi0");
			
			graphPCMPi0InvYieldSysPbPb2760GeV_2050			= ApplyXshiftIndividualSpectra( graphCombInvYieldTot2760GeVALHC11h_2050, 
																							graphPCMPi0InvYieldSysPbPb2760GeV_2050, 
																							fitBylinkinPbPb2760GeVPtMultLHC11h_2050, 
																							0, 18,"Pi0");
			
			graphEMCalPi0InvYieldStatPbPb2760GeV_2050 		= ApplyXshiftIndividualSpectra( graphCombInvYieldTot2760GeVALHC11h_2050, 
																							graphEMCalPi0InvYieldStatPbPb2760GeV_2050, 
																							fitBylinkinPbPb2760GeVPtMultLHC11h_2050,
																							14, 7,"Pi0");	
			
			graphEMCalPi0InvYieldSysPbPb2760GeV_2050 		= ApplyXshiftIndividualSpectra( graphCombInvYieldTot2760GeVALHC11h_2050, 
																							graphEMCalPi0InvYieldSysPbPb2760GeV_2050, 
																							fitBylinkinPbPb2760GeVPtMultLHC11h_2050,
																							14, 7,"Pi0");	
						
            TCanvas* canvasDummy2 = new TCanvas("canvasDummy2","",200,10,1350,1350*1.15);  // gives the page size
            DrawGammaCanvasSettings( canvasDummy2, 0.16, 0.02, 0.02, 0.09);
            canvasDummy2->SetLogx();
            canvasDummy2->SetLogy();

            TH2F * histo2DDummy2 = new TH2F("histo2DDummy2","histo2DDummy2",11000,minPtRange-0.2,maxPtRange,1000,minYaxisYields,maxYaxisYields); 
                SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#pi^{0}}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",0.035,0.04, 0.035,0.04, 1.,1.6);
            histo2DDummy2->GetXaxis()->SetMoreLogLabels();
            histo2DDummy2->GetXaxis()->SetLabelOffset(-0.01);
            histo2DDummy2->Draw("copy");
		
			DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStat2760GeVALHC11hUnShifted_0010, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
			graphCombInvYieldStat2760GeVALHC11hUnShifted_0010->Draw("pEsame");
			DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStat2760GeVALHC11h_0010, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
			graphCombInvYieldStat2760GeVALHC11h_0010->Draw("pEsame");
            
			DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStat2760GeVALHC11hUnShifted_2050, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
			graphCombInvYieldStat2760GeVALHC11hUnShifted_2050->Draw("pEsame");
			DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStat2760GeVALHC11h_2050, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
			graphCombInvYieldStat2760GeVALHC11h_2050->Draw("pEsame");
			
			fitBylinkinPbPb2760GeVPtLHC11h_0010->SetLineColor(kBlue+2);
			fitBylinkinPbPb2760GeVPtLHC11h_0010->Draw("same");
			fitBylinkinPbPb2760GeVPtLHC11h_2050->SetLineColor(kBlue+2);
			fitBylinkinPbPb2760GeVPtLHC11h_2050->Draw("same");

			TLegend* legendXdummy = new TLegend(0.6,0.8,0.9,0.95);
			legendXdummy->SetFillColor(0);
			legendXdummy->SetLineColor(0);
			legendXdummy->SetTextFont(42);
			legendXdummy->SetTextSize(FontSize);
			legendXdummy->AddEntry(graphCombInvYieldStat2760GeVALHC11hUnShifted_0010,"combined unshifted","lp");
			legendXdummy->AddEntry(graphCombInvYieldStat2760GeVALHC11h_0010,"combined shifted","lp");
			legendXdummy->Draw();

			canvasDummy2->Update();
			canvasDummy2->Print(Form("%s/%s_ComparisonShiftedPi0_PbPb2760GeVLHC11h.%s",outputDir.Data(),meson.Data(),suffix.Data()));
			delete canvasDummy2;
            
            
            
            //shifting in Y for both PCM and EMCal
            fitBylinkinPbPb2760GeVPtLHC11hYshift_0010 = FitObject("tcm","BylinkinFitPi00010","Pi0",graphCombInvYieldTot2760GeVALHC11hUnShifted_0010,0.4,20.,paramTCM_0010);
            fitBylinkinPbPb2760GeVPtLHC11hYshift_0010 = ApplyYShift(graphCombInvYieldTot2760GeVALHC11hUnShifted_0010,&graphCombInvYieldTot2760GeVALHC11hYShifted_0010, "tcm", "",0.4 , paramTCM_0010,0.00001,kFALSE,"Pi0");
    //         Double_t parametersBinShiftForQCDfit_0010[5];
    //         GetFitParameter("qcd","0-10%",parametersBinShiftForQCDfit_0010);
    //         for( Int_t i = 0; i < 5; i++){
    //             cout << "parameter " << i << "\t" << parametersBinShiftForQCDfit_0010[i] << endl;
    //         }
    //         fitBylinkinPbPb2760GeVPtLHC11hYshift_0010 = FitObject("qcd","QCDfitPi00010","Pi0",graphCombInvYieldTot2760GeVALHC11hUnShifted_0010, 0.4,20.,parametersBinShiftForQCDfit_0010);
    //         fitBylinkinPbPb2760GeVPtLHC11hYshift_0010 = ApplyYShift(graphCombInvYieldTot2760GeVALHC11hUnShifted_0010,&graphCombInvYieldTot2760GeVALHC11hYShifted_0010, "qcd", "",0.4,parametersBinShiftForQCDfit_0010,0.00001,kFALSE,"Pi0");
          
            cout << "combined binshift Y" << endl;
            graphCombInvYieldSys2760GeVALHC11hYShifted_0010 = (TGraphAsymmErrors*)graphCombInvYieldSys2760GeVALHC11hUnShifted_0010->Clone("YShiftedCombSys0010");
            graphCombInvYieldSys2760GeVALHC11hYShifted_0010= ApplyYshiftIndividualSpectra( graphCombInvYieldSys2760GeVALHC11hYShifted_0010, fitBylinkinPbPb2760GeVPtLHC11hYshift_0010);
            graphCombInvYieldStat2760GeVALHC11hYShifted_0010 = (TGraphAsymmErrors*)graphCombInvYieldStat2760GeVALHC11hUnShifted_0010->Clone("YShiftedCombStat0010");
            graphCombInvYieldStat2760GeVALHC11hYShifted_0010= ApplyYshiftIndividualSpectra( graphCombInvYieldStat2760GeVALHC11hYShifted_0010, fitBylinkinPbPb2760GeVPtLHC11hYshift_0010);      

            cout << "PCM binshift Y" << endl;
            graphPCMInvYieldSysPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphPCMInvYieldSysPbPb2760GeVUnShifted_0010->Clone("YShiftedPCMSys0010");
            graphPCMInvYieldSysPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphPCMInvYieldSysPbPb2760GeVYShifted_0010, fitBylinkinPbPb2760GeVPtLHC11hYshift_0010);
            graphPCMInvYieldStatPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphPCMInvYieldStatPbPb2760GeVUnShifted_0010->Clone("YShiftedPCMStat0010");
            graphPCMInvYieldStatPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphPCMInvYieldStatPbPb2760GeVYShifted_0010, fitBylinkinPbPb2760GeVPtLHC11hYshift_0010);

//             graphPCMPi0RAAStatPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphPCMPi0RAAStatPbPb2760GeV_0010->Clone("YShiftedPCMSysRAA0010");
//             graphPCMPi0RAAStatPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphPCMPi0RAAStatPbPb2760GeVYShifted_0010, fitBylinkinPbPb2760GeVPtLHC11hYshift_0010);
//             graphPCMPi0RAASysPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphPCMPi0RAASysPbPb2760GeV_0010->Clone("YShiftedPCMSysRAA0010");
//             graphPCMPi0RAASysPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphPCMPi0RAASysPbPb2760GeVYShifted_0010, fitBylinkinPbPb2760GeVPtLHC11hYshift_0010);

            cout << "EMCal binshift Y" << endl;
            graphEMCalInvYieldSysPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphEMCalInvYieldSysPbPb2760GeVUnshifted_0010->Clone("YShiftedEMCalSys0010");
            graphEMCalInvYieldSysPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphEMCalInvYieldSysPbPb2760GeVYShifted_0010, fitBylinkinPbPb2760GeVPtLHC11hYshift_0010);
            graphEMCalInvYieldStatPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphEMCalInvYieldStatPbPb2760GeVUnshifted_0010->Clone("YShiftedEMCalStat0010");
            graphEMCalInvYieldStatPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphEMCalInvYieldStatPbPb2760GeVYShifted_0010, fitBylinkinPbPb2760GeVPtLHC11hYshift_0010);

//             graphEMCalPi0RAAStatPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphEMCalPi0RAAStatPbPb2760GeV_0010->Clone("YShiftedEMCalSysRAA0010");
//             graphEMCalPi0RAAStatPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphEMCalPi0RAAStatPbPb2760GeVYShifted_0010, fitBylinkinPbPb2760GeVPtLHC11hYshift_0010);
//             graphEMCalPi0RAASysPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphEMCalPi0RAASysPbPb2760GeV_0010->Clone("YShiftedEMCalSysRAA0010");
//             graphEMCalPi0RAASysPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphEMCalPi0RAASysPbPb2760GeVYShifted_0010, fitBylinkinPbPb2760GeVPtLHC11hYshift_0010);

            
            
            fitBylinkinPbPb2760GeVPtLHC11hYshift_2050 = FitObject("tcm","BylinkinFitPi02050","Pi0",graphCombInvYieldTot2760GeVALHC11hUnShifted_2050,0.4,20.,paramTCM_2050);
            fitBylinkinPbPb2760GeVPtLHC11hYshift_2050 = ApplyYShift(graphCombInvYieldTot2760GeVALHC11hUnShifted_2050,&graphCombInvYieldTot2760GeVALHC11hYShifted_2050, "tcm", "",0.4, paramTCM_2050,0.00001,kFALSE,"Pi0");
    //         Double_t parametersBinShiftForQCDfit_2050[5];
    //         GetFitParameter("qcd","20-50%",parametersBinShiftForQCDfit_2050);
    //         for( Int_t i = 0; i < 5; i++){
    //             cout << "parameter " << i << "\t" << parametersBinShiftForQCDfit_2050[i] << endl;
    //         }
    //         fitBylinkinPbPb2760GeVPtLHC11hYshift_2050 = FitObject("qcd","QCDfitPi02050","Pi0",graphCombInvYieldTot2760GeVALHC11hUnShifted_2050, 0.4,20.,parametersBinShiftForQCDfit_2050);
    //         fitBylinkinPbPb2760GeVPtLHC11hYshift_2050 = ApplyYShift(graphCombInvYieldTot2760GeVALHC11hUnShifted_2050,&graphCombInvYieldTot2760GeVALHC11hYShifted_2050, "qcd", "",0.4,parametersBinShiftForQCDfit_2050,0.00001,kFALSE,"Pi0");

            
            cout << "combined binshift Y" << endl;
            graphCombInvYieldSys2760GeVALHC11hYShifted_2050 = (TGraphAsymmErrors*)graphCombInvYieldSys2760GeVALHC11hUnShifted_2050->Clone("YShiftedCombSys2050");
            graphCombInvYieldSys2760GeVALHC11hYShifted_2050= ApplyYshiftIndividualSpectra( graphCombInvYieldSys2760GeVALHC11hYShifted_2050, fitBylinkinPbPb2760GeVPtLHC11hYshift_2050);
            graphCombInvYieldStat2760GeVALHC11hYShifted_2050 = (TGraphAsymmErrors*)graphCombInvYieldStat2760GeVALHC11hUnShifted_2050->Clone("YShiftedCombStat2050");
            graphCombInvYieldStat2760GeVALHC11hYShifted_2050= ApplyYshiftIndividualSpectra( graphCombInvYieldStat2760GeVALHC11hYShifted_2050, fitBylinkinPbPb2760GeVPtLHC11hYshift_2050);      

            cout << "PCM binshift Y" << endl;
            graphPCMInvYieldSysPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphPCMInvYieldSysPbPb2760GeVUnShifted_2050->Clone("YShiftedPCMSys2050");
            graphPCMInvYieldSysPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphPCMInvYieldSysPbPb2760GeVYShifted_2050, fitBylinkinPbPb2760GeVPtLHC11hYshift_2050);
            graphPCMInvYieldStatPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphPCMInvYieldStatPbPb2760GeVUnShifted_2050->Clone("YShiftedPCMStat2050");
            graphPCMInvYieldStatPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphPCMInvYieldStatPbPb2760GeVYShifted_2050, fitBylinkinPbPb2760GeVPtLHC11hYshift_2050);

//             graphPCMPi0RAAStatPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphPCMPi0RAAStatPbPb2760GeV_2050->Clone("YShiftedPCMSysRAA2050");
//             graphPCMPi0RAAStatPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphPCMPi0RAAStatPbPb2760GeVYShifted_2050, fitBylinkinPbPb2760GeVPtLHC11hYshift_2050);
//             graphPCMPi0RAASysPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphPCMPi0RAASysPbPb2760GeV_2050->Clone("YShiftedPCMSysRAA2050");
//             graphPCMPi0RAASysPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphPCMPi0RAASysPbPb2760GeVYShifted_2050, fitBylinkinPbPb2760GeVPtLHC11hYshift_2050);

            cout << "EMCal binshift Y" << endl;
            graphEMCalInvYieldSysPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphEMCalInvYieldSysPbPb2760GeVUnshifted_2050->Clone("YShiftedEMCalSys2050");
            graphEMCalInvYieldSysPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphEMCalInvYieldSysPbPb2760GeVYShifted_2050, fitBylinkinPbPb2760GeVPtLHC11hYshift_2050);
            graphEMCalInvYieldStatPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphEMCalInvYieldStatPbPb2760GeVUnshifted_2050->Clone("YShiftedEMCalStat2050");
            graphEMCalInvYieldStatPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphEMCalInvYieldStatPbPb2760GeVYShifted_2050, fitBylinkinPbPb2760GeVPtLHC11hYshift_2050);

//             graphEMCalPi0RAAStatPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphEMCalPi0RAAStatPbPb2760GeV_2050->Clone("YShiftedEMCalSysRAA2050");
//             graphEMCalPi0RAAStatPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphEMCalPi0RAAStatPbPb2760GeVYShifted_2050, fitBylinkinPbPb2760GeVPtLHC11hYshift_2050);
//             graphEMCalPi0RAASysPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphEMCalPi0RAASysPbPb2760GeV_2050->Clone("YShiftedEMCalSysRAA2050");
//             graphEMCalPi0RAASysPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphEMCalPi0RAASysPbPb2760GeVYShifted_2050, fitBylinkinPbPb2760GeVPtLHC11hYshift_2050);

            cout << WriteParameterToFile(fitBylinkinPbPb2760GeVPtLHC11hYshift_0010)<< endl << endl;
            cout << WriteParameterToFile(fitBylinkinPbPb2760GeVPtLHC11hYshift_2050)<< endl << endl;

            TCanvas* canvasDummy3 = new TCanvas("canvasDummy3","",200,10,1350,1350*1.15);  // gives the page size
            DrawGammaCanvasSettings( canvasDummy3, 0.16, 0.02, 0.02, 0.09);
            canvasDummy3->SetLogx();
            canvasDummy3->SetLogy();

            TH2F * histo2DDummy3 = new TH2F("histo2DDummy3","histo2DDummy3",11000,minPtRange-0.2,maxPtRange,1000,minYaxisYields,maxYaxisYields); 
                SetStyleHistoTH2ForGraphs(histo2DDummy3, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#pi^{0}}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",0.035,0.04, 0.035,0.04, 1.,1.6);
            histo2DDummy3->GetXaxis()->SetMoreLogLabels();
            histo2DDummy3->GetXaxis()->SetLabelOffset(-0.01);
            histo2DDummy3->Draw("copy");
        
            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStat2760GeVALHC11hUnShifted_0010, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
            graphCombInvYieldStat2760GeVALHC11hUnShifted_0010->Draw("pEsame");
            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStat2760GeVALHC11hYShifted_0010, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
            graphCombInvYieldStat2760GeVALHC11hYShifted_0010->Draw("pEsame");
            
            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStat2760GeVALHC11hUnShifted_2050, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
            graphCombInvYieldStat2760GeVALHC11hUnShifted_2050->Draw("pEsame");
            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStat2760GeVALHC11hYShifted_2050, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
            graphCombInvYieldStat2760GeVALHC11hYShifted_2050->Draw("pEsame");
            
    //         DrawGammaSetMarkerTGraphAsym(graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_0010, 25, 1.5, kBlue, kBlue, widthLinesBoxes, kTRUE);
    //         graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_0010->Draw("pEsame");
    //         DrawGammaSetMarkerTGraphAsym(graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_2050, 25, 1.5, kBlue, kBlue, widthLinesBoxes, kTRUE);
    //         graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_2050->Draw("pEsame");

            fitBylinkinPbPb2760GeVPtLHC11hYshift_0010->SetLineColor(kBlue+2);
            fitBylinkinPbPb2760GeVPtLHC11hYshift_0010->Draw("same");
            fitBylinkinPbPb2760GeVPtLHC11hYshift_2050->SetLineColor(kBlue+2);
            fitBylinkinPbPb2760GeVPtLHC11hYshift_2050->Draw("same");

            TLegend* legendYdummy = new TLegend(0.6,0.8,0.9,0.95);
            legendYdummy->SetFillColor(0);
            legendYdummy->SetLineColor(0);
            legendYdummy->SetTextFont(42);
            legendYdummy->SetTextSize(FontSize);
            legendYdummy->AddEntry(graphCombInvYieldStat2760GeVALHC11hUnShifted_0010,"combined unshifted","lp");
            legendYdummy->AddEntry(graphCombInvYieldStat2760GeVALHC11hYShifted_0010,"combined shifted","lp");
            legendYdummy->Draw();

            canvasDummy3->Update();
            canvasDummy3->Print(Form("%s/%s_ComparisonYShifted.%s",outputDir.Data(),meson.Data(),suffix.Data()));
            delete canvasDummy3;
            
		}
		
        
	} else if(meson.CompareTo("Eta")==0){
        Double_t paramTCM_0010[5] = {graphCombInvYieldTot2760GeVALHC11hUnShifted_0010->GetY()[0],0.3,graphCombInvYieldTot2760GeVALHC11hUnShifted_0010->GetY()[0],0.3,8};
		fitBylinkinPbPb2760GeVPtLHC11h_0010 = FitObject("tcm","BylinkinFitEta0010","Eta",graphCombInvYieldTot2760GeVALHC11hUnShifted_0010,1.,20., paramTCM_0010);
		graphCombInvYieldTot2760GeVALHC11hUnShifted_0010->Fit(fitBylinkinPbPb2760GeVPtLHC11h_0010,"NRMEX0+","",1.,20.);

        Double_t paramTCM_2050[5] = {graphCombInvYieldTot2760GeVALHC11hUnShifted_2050->GetY()[0],0.3,graphCombInvYieldTot2760GeVALHC11hUnShifted_2050->GetY()[0],0.3,8};        
		fitBylinkinPbPb2760GeVPtLHC11h_2050 = FitObject("tcm","BylinkinFitEta2050","Eta",graphCombInvYieldTot2760GeVALHC11hUnShifted_2050,1.,20., paramTCM_2050);
		graphCombInvYieldTot2760GeVALHC11hUnShifted_2050->Fit(fitBylinkinPbPb2760GeVPtLHC11h_2050,"NRMEX0+","",1.,20.);

        cout << WriteParameterToFile(fitBylinkinPbPb2760GeVPtLHC11h_0010)<< endl << endl;
        cout << WriteParameterToFile(fitBylinkinPbPb2760GeVPtLHC11h_2050)<< endl << endl;

        if(bWCorrection.CompareTo("X")==0 ){
			TF1* fitBylinkinPbPb2760GeVPtMultLHC11h_0010 = FitObject("tcmpt","BylinkinFitEta0010","Eta",graphCombInvYieldTot2760GeVALHC11h_0010);
			
			graphCombInvYieldTot2760GeVALHC11h_0010			= ApplyXshift(graphCombInvYieldTot2760GeVALHC11h_0010, fitBylinkinPbPb2760GeVPtMultLHC11h_0010,"Eta");
			
			graphCombInvYieldStat2760GeVALHC11h_0010 		= ApplyXshiftIndividualSpectra (graphCombInvYieldTot2760GeVALHC11h_0010, 
																							graphCombInvYieldStat2760GeVALHC11h_0010, 
																							fitBylinkinPbPb2760GeVPtMultLHC11h_0010,
																							0, 10,"Eta");
			graphCombInvYieldSys2760GeVALHC11h_0010 		= ApplyXshiftIndividualSpectra (graphCombInvYieldTot2760GeVALHC11h_0010, 
																							graphCombInvYieldSys2760GeVALHC11h_0010, 
																							fitBylinkinPbPb2760GeVPtMultLHC11h_0010, 
																							0, 10,"Eta");
			graphPCMEtaInvYieldStatPbPb2760GeV_0010			= ApplyXshiftIndividualSpectra( graphCombInvYieldTot2760GeVALHC11h_0010,
																							graphPCMEtaInvYieldStatPbPb2760GeV_0010,
																							fitBylinkinPbPb2760GeVPtMultLHC11h_0010, 
																							0, 7,"Eta");
			graphPCMEtaInvYieldSysPbPb2760GeV_0010			= ApplyXshiftIndividualSpectra( graphCombInvYieldTot2760GeVALHC11h_0010, 
																							graphPCMEtaInvYieldSysPbPb2760GeV_0010, 
																							fitBylinkinPbPb2760GeVPtMultLHC11h_0010, 
																							0, 7,"Eta");
			
			graphEMCalEtaInvYieldStatPbPb2760GeV_0010 		= ApplyXshiftIndividualSpectra( graphCombInvYieldTot2760GeVALHC11h_0010, 
																							graphEMCalEtaInvYieldStatPbPb2760GeV_0010, 
																							fitBylinkinPbPb2760GeVPtMultLHC11h_0010,
																							4, 7,"Eta");	
			
			graphEMCalEtaInvYieldSysPbPb2760GeV_0010 		= ApplyXshiftIndividualSpectra( graphCombInvYieldTot2760GeVALHC11h_0010, 
																							graphEMCalEtaInvYieldSysPbPb2760GeV_0010, 
																							fitBylinkinPbPb2760GeVPtMultLHC11h_0010,
																							4, 7,"Eta");	

			TF1* fitBylinkinPbPb2760GeVPtMultLHC11h_2050 = FitObject("tcmpt","BylinkinFitEta2050","Eta",graphCombInvYieldTot2760GeVALHC11h_2050);
			
			graphCombInvYieldTot2760GeVALHC11h_2050			= ApplyXshift(graphCombInvYieldTot2760GeVALHC11h_2050, fitBylinkinPbPb2760GeVPtMultLHC11h_2050,"Eta");
			
			graphCombInvYieldStat2760GeVALHC11h_2050 		= ApplyXshiftIndividualSpectra (graphCombInvYieldTot2760GeVALHC11h_2050, 
																							graphCombInvYieldStat2760GeVALHC11h_2050, 
																							fitBylinkinPbPb2760GeVPtMultLHC11h_2050,
																							0, 10,"Eta");
			
			graphCombInvYieldSys2760GeVALHC11h_2050 		= ApplyXshiftIndividualSpectra (graphCombInvYieldTot2760GeVALHC11h_2050, 
																							graphCombInvYieldSys2760GeVALHC11h_2050, 
																							fitBylinkinPbPb2760GeVPtMultLHC11h_2050, 
																							0, 10,"Eta");
			
			graphPCMEtaInvYieldStatPbPb2760GeV_2050			= ApplyXshiftIndividualSpectra( graphCombInvYieldTot2760GeVALHC11h_2050,
																							graphPCMEtaInvYieldStatPbPb2760GeV_2050,
																							fitBylinkinPbPb2760GeVPtMultLHC11h_2050, 
																							0, 7,"Eta");
			
			graphPCMEtaInvYieldSysPbPb2760GeV_2050			= ApplyXshiftIndividualSpectra( graphCombInvYieldTot2760GeVALHC11h_2050, 
																							graphPCMEtaInvYieldSysPbPb2760GeV_2050, 
																							fitBylinkinPbPb2760GeVPtMultLHC11h_2050, 
																							0, 7,"Eta");
			
			graphEMCalEtaInvYieldStatPbPb2760GeV_2050 		= ApplyXshiftIndividualSpectra( graphCombInvYieldTot2760GeVALHC11h_2050, 
																							graphEMCalEtaInvYieldStatPbPb2760GeV_2050, 
																							fitBylinkinPbPb2760GeVPtMultLHC11h_2050,
																							4, 7,"Eta");	
			
			graphEMCalEtaInvYieldSysPbPb2760GeV_2050 		= ApplyXshiftIndividualSpectra( graphCombInvYieldTot2760GeVALHC11h_2050, 
																							graphEMCalEtaInvYieldSysPbPb2760GeV_2050, 
																							fitBylinkinPbPb2760GeVPtMultLHC11h_2050,
																							4, 7,"Eta");	
			
					
            TCanvas* canvasDummy2 = new TCanvas("canvasDummy2","",200,10,1350,1350*1.15);  // gives the page size
            DrawGammaCanvasSettings( canvasDummy2, 0.16, 0.02, 0.02, 0.09);
            canvasDummy2->SetLogx();
            canvasDummy2->SetLogy();
            
            TH2F * histo2DDummy2 = new TH2F("histo2DDummy2","histo2DDummy2",11000,minPtRange-0.2,maxPtRange,1000,minYaxisYields,maxYaxisYields);
                SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#eta}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",0.035,0.04, 0.035,0.04, 1.,1.6);            
            histo2DDummy2->GetXaxis()->SetMoreLogLabels();
            histo2DDummy2->GetXaxis()->SetLabelOffset(-0.01);
            histo2DDummy2->Draw("copy");
		
			DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStat2760GeVALHC11hUnShifted_0010, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
			graphCombInvYieldStat2760GeVALHC11hUnShifted_0010->Draw("pEsame");
			DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStat2760GeVALHC11h_0010, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
			graphCombInvYieldStat2760GeVALHC11h_0010->Draw("pEsame");
			fitBylinkinPbPb2760GeVPtLHC11h_0010->SetLineColor(kBlue+2);
			fitBylinkinPbPb2760GeVPtLHC11h_0010->Draw("same");

			DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStat2760GeVALHC11hUnShifted_2050, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
			graphCombInvYieldStat2760GeVALHC11hUnShifted_2050->Draw("pEsame");
			DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStat2760GeVALHC11h_2050, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
			graphCombInvYieldStat2760GeVALHC11h_2050->Draw("pEsame");	
			
			fitBylinkinPbPb2760GeVPtLHC11h_2050->SetLineColor(kBlue+2);
			fitBylinkinPbPb2760GeVPtLHC11h_2050->Draw("same");

			TLegend* legendXdummy = new TLegend(0.6,0.8,0.9,0.95);
			legendXdummy->SetFillColor(0);
			legendXdummy->SetLineColor(0);
			legendXdummy->SetTextFont(42);
			legendXdummy->SetTextSize(FontSize);
			legendXdummy->AddEntry(graphCombInvYieldStat2760GeVALHC11hUnShifted_0010,"combined unshifted","lp");
			legendXdummy->AddEntry(graphCombInvYieldStat2760GeVALHC11h_0010,"combined shifted","lp");
			legendXdummy->Draw();
			
			canvasDummy2->Update();
			canvasDummy2->Print(Form("%s/%s_ComparisonShiftedEta_PbPb2760GeVLHC11h.%s",outputDir.Data(),meson.Data(),suffix.Data()));
			delete canvasDummy2;
            
            
            
            //shifting in Y for both PCM and EMCal
            
            fitBylinkinPbPb2760GeVPtLHC11hYshift_0010 = FitObject("tcm","BylinkinFitEta0010","Eta",graphCombInvYieldTot2760GeVALHC11hUnShifted_0010,1.,20.,paramTCM_0010);
            fitBylinkinPbPb2760GeVPtLHC11hYshift_0010 = ApplyYShift(graphCombInvYieldTot2760GeVALHC11hUnShifted_0010,&graphCombInvYieldTot2760GeVALHC11hYShifted_0010, "tcm", "",1. , paramTCM_0010,0.00001,kFALSE,"Eta");
    //         Double_t parametersBinShiftForQCDfit_0010[5];
    //         GetFitParameter("qcd","0-10%",parametersBinShiftForQCDfit_0010);
    //         for( Int_t i = 0; i < 5; i++){
    //             cout << "parameter " << i << "\t" << parametersBinShiftForQCDfit_0010[i] << endl;
    //         }
    //         fitBylinkinPbPb2760GeVPtLHC11hYshift_0010 = FitObject("qcd","QCDfitEta0010","Eta",graphCombInvYieldTot2760GeVALHC11hUnShifted_0010, 1.,20.,parametersBinShiftForQCDfit_0010);
    //         fitBylinkinPbPb2760GeVPtLHC11hYshift_0010 = ApplyYShift(graphCombInvYieldTot2760GeVALHC11hUnShifted_0010,&graphCombInvYieldTot2760GeVALHC11hYShifted_0010, "qcd", "",1.,parametersBinShiftForQCDfit_0010,0.00001,kFALSE,"Eta");
          
            cout << "combined binshift Y" << endl;
            graphCombInvYieldSys2760GeVALHC11hYShifted_0010 = (TGraphAsymmErrors*)graphCombInvYieldSys2760GeVALHC11hUnShifted_0010->Clone("YShiftedCombSys0010");
            graphCombInvYieldSys2760GeVALHC11hYShifted_0010= ApplyYshiftIndividualSpectra( graphCombInvYieldSys2760GeVALHC11hYShifted_0010, fitBylinkinPbPb2760GeVPtLHC11hYshift_0010);
            graphCombInvYieldStat2760GeVALHC11hYShifted_0010 = (TGraphAsymmErrors*)graphCombInvYieldStat2760GeVALHC11hUnShifted_0010->Clone("YShiftedCombStat0010");
            graphCombInvYieldStat2760GeVALHC11hYShifted_0010= ApplyYshiftIndividualSpectra( graphCombInvYieldStat2760GeVALHC11hYShifted_0010, fitBylinkinPbPb2760GeVPtLHC11hYshift_0010);      

            cout << "PCM binshift Y" << endl;
            graphPCMInvYieldSysPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphPCMInvYieldSysPbPb2760GeVUnShifted_0010->Clone("YShiftedPCMSys0010");
            graphPCMInvYieldSysPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphPCMInvYieldSysPbPb2760GeVYShifted_0010, fitBylinkinPbPb2760GeVPtLHC11hYshift_0010);
            graphPCMInvYieldStatPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphPCMInvYieldStatPbPb2760GeVUnShifted_0010->Clone("YShiftedPCMStat0010");
            graphPCMInvYieldStatPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphPCMInvYieldStatPbPb2760GeVYShifted_0010, fitBylinkinPbPb2760GeVPtLHC11hYshift_0010);

//             graphPCMEtaRAAStatPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphPCMEtaRAAStatPbPb2760GeV_0010->Clone("YShiftedPCMSysRAA0010");
//             graphPCMEtaRAAStatPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphPCMEtaRAAStatPbPb2760GeVYShifted_0010, fitBylinkinPbPb2760GeVPtLHC11hYshift_0010);
//             graphPCMEtaRAASysPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphPCMEtaRAASysPbPb2760GeV_0010->Clone("YShiftedPCMSysRAA0010");
//             graphPCMEtaRAASysPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphPCMEtaRAASysPbPb2760GeVYShifted_0010, fitBylinkinPbPb2760GeVPtLHC11hYshift_0010);
            
            cout << "EMCal binshift Y" << endl;
            graphEMCalInvYieldSysPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphEMCalInvYieldSysPbPb2760GeVUnshifted_0010->Clone("YShiftedPHOSSys0010");
            graphEMCalInvYieldSysPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphEMCalInvYieldSysPbPb2760GeVYShifted_0010, fitBylinkinPbPb2760GeVPtLHC11hYshift_0010);
            graphEMCalInvYieldStatPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphEMCalInvYieldStatPbPb2760GeVUnshifted_0010->Clone("YShiftedPHOSStat0010");
            graphEMCalInvYieldStatPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphEMCalInvYieldStatPbPb2760GeVYShifted_0010, fitBylinkinPbPb2760GeVPtLHC11hYshift_0010);

//             graphEMCalEtaRAAStatPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphEMCalEtaRAAStatPbPb2760GeV_0010->Clone("YShiftedEMCalSysRAA0010");
//             graphEMCalEtaRAAStatPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphEMCalEtaRAAStatPbPb2760GeVYShifted_0010, fitBylinkinPbPb2760GeVPtLHC11hYshift_0010);
//             graphEMCalEtaRAASysPbPb2760GeVYShifted_0010 = (TGraphAsymmErrors*)graphEMCalEtaRAASysPbPb2760GeV_0010->Clone("YShiftedEMCalSysRAA0010");
//             graphEMCalEtaRAASysPbPb2760GeVYShifted_0010 = ApplyYshiftIndividualSpectra( graphEMCalEtaRAASysPbPb2760GeVYShifted_0010, fitBylinkinPbPb2760GeVPtLHC11hYshift_0010);

            
            
            
            fitBylinkinPbPb2760GeVPtLHC11hYshift_2050 = FitObject("tcm","BylinkinFitEta2050","Eta",graphCombInvYieldTot2760GeVALHC11hUnShifted_2050,1.,20.,paramTCM_2050);
            fitBylinkinPbPb2760GeVPtLHC11hYshift_2050 = ApplyYShift(graphCombInvYieldTot2760GeVALHC11hUnShifted_2050,&graphCombInvYieldTot2760GeVALHC11hYShifted_2050, "tcm", "",1., paramTCM_2050,0.00001,kFALSE,"Eta");
    //         Double_t parametersBinShiftForQCDfit_2050[5];
    //         GetFitParameter("qcd","20-50%",parametersBinShiftForQCDfit_2050);
    //         for( Int_t i = 0; i < 5; i++){
    //             cout << "parameter " << i << "\t" << parametersBinShiftForQCDfit_2050[i] << endl;
    //         }
    //         fitBylinkinPbPb2760GeVPtLHC11hYshift_2050 = FitObject("qcd","QCDfitEta2050","Eta",graphCombInvYieldTot2760GeVALHC11hUnShifted_2050, 1.,20.,parametersBinShiftForQCDfit_2050);
    //         fitBylinkinPbPb2760GeVPtLHC11hYshift_2050 = ApplyYShift(graphCombInvYieldTot2760GeVALHC11hUnShifted_2050,&graphCombInvYieldTot2760GeVALHC11hYShifted_2050, "qcd", "",1.,parametersBinShiftForQCDfit_2050,0.00001,kFALSE,"Eta");

            
            cout << "combined binshift Y" << endl;
            graphCombInvYieldSys2760GeVALHC11hYShifted_2050 = (TGraphAsymmErrors*)graphCombInvYieldSys2760GeVALHC11hUnShifted_2050->Clone("YShiftedCombSys2050");
            graphCombInvYieldSys2760GeVALHC11hYShifted_2050= ApplyYshiftIndividualSpectra( graphCombInvYieldSys2760GeVALHC11hYShifted_2050, fitBylinkinPbPb2760GeVPtLHC11hYshift_2050);
            graphCombInvYieldStat2760GeVALHC11hYShifted_2050 = (TGraphAsymmErrors*)graphCombInvYieldStat2760GeVALHC11hUnShifted_2050->Clone("YShiftedCombStat2050");
            graphCombInvYieldStat2760GeVALHC11hYShifted_2050= ApplyYshiftIndividualSpectra( graphCombInvYieldStat2760GeVALHC11hYShifted_2050, fitBylinkinPbPb2760GeVPtLHC11hYshift_2050);      

            cout << "PCM binshift Y" << endl;
            graphPCMInvYieldSysPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphPCMInvYieldSysPbPb2760GeVUnShifted_2050->Clone("YShiftedPCMSys2050");
            graphPCMInvYieldSysPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphPCMInvYieldSysPbPb2760GeVYShifted_2050, fitBylinkinPbPb2760GeVPtLHC11hYshift_2050);
            graphPCMInvYieldStatPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphPCMInvYieldStatPbPb2760GeVUnShifted_2050->Clone("YShiftedPCMStat2050");
            graphPCMInvYieldStatPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphPCMInvYieldStatPbPb2760GeVYShifted_2050, fitBylinkinPbPb2760GeVPtLHC11hYshift_2050);

//             graphPCMEtaRAAStatPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphPCMEtaRAAStatPbPb2760GeV_2050->Clone("YShiftedPCMSysRAA2050");
//             graphPCMEtaRAAStatPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphPCMEtaRAAStatPbPb2760GeVYShifted_2050, fitBylinkinPbPb2760GeVPtLHC11hYshift_2050);
//             graphPCMEtaRAASysPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphPCMEtaRAASysPbPb2760GeV_2050->Clone("YShiftedPCMSysRAA2050");
//             graphPCMEtaRAASysPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphPCMEtaRAASysPbPb2760GeVYShifted_2050, fitBylinkinPbPb2760GeVPtLHC11hYshift_2050);
            
            cout << "EMCal binshift Y" << endl;
            graphEMCalInvYieldSysPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphEMCalInvYieldSysPbPb2760GeVUnshifted_2050->Clone("YShiftedPHOSSys2050");
            graphEMCalInvYieldSysPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphEMCalInvYieldSysPbPb2760GeVYShifted_2050, fitBylinkinPbPb2760GeVPtLHC11hYshift_2050);
            graphEMCalInvYieldStatPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphEMCalInvYieldStatPbPb2760GeVUnshifted_2050->Clone("YShiftedPHOSStat2050");
            graphEMCalInvYieldStatPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphEMCalInvYieldStatPbPb2760GeVYShifted_2050, fitBylinkinPbPb2760GeVPtLHC11hYshift_2050);

//             graphEMCalEtaRAAStatPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphEMCalEtaRAAStatPbPb2760GeV_2050->Clone("YShiftedEMCalSysRAA2050");
//             graphEMCalEtaRAAStatPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphEMCalEtaRAAStatPbPb2760GeVYShifted_2050, fitBylinkinPbPb2760GeVPtLHC11hYshift_2050);
//             graphEMCalEtaRAASysPbPb2760GeVYShifted_2050 = (TGraphAsymmErrors*)graphEMCalEtaRAASysPbPb2760GeV_2050->Clone("YShiftedEMCalSysRAA2050");
//             graphEMCalEtaRAASysPbPb2760GeVYShifted_2050 = ApplyYshiftIndividualSpectra( graphEMCalEtaRAASysPbPb2760GeVYShifted_2050, fitBylinkinPbPb2760GeVPtLHC11hYshift_2050);

            cout << WriteParameterToFile(fitBylinkinPbPb2760GeVPtLHC11hYshift_0010)<< endl << endl;
            cout << WriteParameterToFile(fitBylinkinPbPb2760GeVPtLHC11hYshift_2050)<< endl << endl;

            TCanvas* canvasDummy3 = new TCanvas("canvasDummy3","",200,10,1350,1350*1.15);  // gives the page size
            DrawGammaCanvasSettings( canvasDummy3, 0.16, 0.02, 0.02, 0.09);
            canvasDummy3->SetLogx();
            canvasDummy3->SetLogy();

            TH2F * histo2DDummy3 = new TH2F("histo2DDummy3","histo2DDummy3",11000,minPtRange-0.2,maxPtRange,1000,minYaxisYields,maxYaxisYields); 
                SetStyleHistoTH2ForGraphs(histo2DDummy3, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#eta}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",0.035,0.04, 0.035,0.04, 1.,1.6);            
            histo2DDummy3->GetXaxis()->SetMoreLogLabels();
            histo2DDummy3->GetXaxis()->SetLabelOffset(-0.01);
            histo2DDummy3->Draw("copy");
        
            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStat2760GeVALHC11hUnShifted_0010, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
            graphCombInvYieldStat2760GeVALHC11hUnShifted_0010->Draw("pEsame");
            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStat2760GeVALHC11hYShifted_0010, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
            graphCombInvYieldStat2760GeVALHC11hYShifted_0010->Draw("pEsame");
            
            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStat2760GeVALHC11hUnShifted_2050, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
            graphCombInvYieldStat2760GeVALHC11hUnShifted_2050->Draw("pEsame");
            DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStat2760GeVALHC11hYShifted_2050, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
            graphCombInvYieldStat2760GeVALHC11hYShifted_2050->Draw("pEsame");
            
//             DrawGammaSetMarkerTGraphAsym(graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_0010, 25, 1.5, kBlue, kBlue, widthLinesBoxes, kTRUE);
//             graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_0010->Draw("pEsame");
//             DrawGammaSetMarkerTGraphAsym(graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_2050, 25, 1.5, kBlue, kBlue, widthLinesBoxes, kTRUE);
//             graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_2050->Draw("pEsame");

            fitBylinkinPbPb2760GeVPtLHC11hYshift_0010->SetLineColor(kBlue+2);
            fitBylinkinPbPb2760GeVPtLHC11hYshift_0010->Draw("same");
            fitBylinkinPbPb2760GeVPtLHC11hYshift_2050->SetLineColor(kBlue+2);
            fitBylinkinPbPb2760GeVPtLHC11hYshift_2050->Draw("same");

            TLegend* legendYdummy = new TLegend(0.6,0.8,0.9,0.95);
            legendYdummy->SetFillColor(0);
            legendYdummy->SetLineColor(0);
            legendYdummy->SetTextFont(42);
            legendYdummy->SetTextSize(FontSize);
            legendYdummy->AddEntry(graphCombInvYieldStat2760GeVALHC11hUnShifted_0010,"combined unshifted","lp");
            legendYdummy->AddEntry(graphCombInvYieldStat2760GeVALHC11hYShifted_0010,"combined shifted","lp");
            legendYdummy->Draw();

            canvasDummy3->Update();
            canvasDummy3->Print(Form("%s/%s_ComparisonYShifted.%s",outputDir.Data(),meson.Data(),suffix.Data()));
            delete canvasDummy3;
            
		}
		
	}

	
	
	
	TF1* fitTCMInvYield2760GeVLHC11h_0010;
	TF1* fitTCMInvYield2760GeVLHC11h_2050;
    TF1* fitTCMInvYield2760GeVas10h_0010;
    TF1* fitTCMInvYield2760GeVas10h_2050;
    Double_t minfitPt = 0.4;
	if(meson.CompareTo("Pi0")==0){
        if(PaperPi0) minfitPt = 1.0;
        else minfitPt = 0.4;
		cout << "Fit to the 0-10% spectrum" << endl;
		Double_t paramTCM_0010[5] = {graphCombInvYieldTot2760GeVALHC11h_0010->GetY()[0],0.3,graphCombInvYieldTot2760GeVALHC11h_0010->GetY()[0],0.3,8};
		fitTCMInvYield2760GeVLHC11h_0010 = FitObject("tcm","fitBylinkin0010","Pi0",graphCombInvYieldTot2760GeVALHC11h_0010,minfitPt,20.,paramTCM_0010,"NRMEX0+");
        graphCombInvYieldTot2760GeVALHC11h_0010->Fit(fitTCMInvYield2760GeVLHC11h_0010,"NRMEX0+","",minfitPt,20.);

//         fitTCMInvYield2760GeVLHC11h_0010 = FitObject("qcd","fitQCD","Pi0",graphCombInvYieldTot2760GeVALHC11h_0010,0.4,20.);
//         graphCombInvYieldTot2760GeVALHC11h_0010->Fit(fitTCMInvYield2760GeVLHC11h_0010,"NRMEX0+","",0.4,20.);
        
        fitTCMInvYield2760GeVas10h_0010 = FitObject("qcd","fitQCD","Pi0",graphCombInvYieldTot2760GeVALHC11h_0010,0.4,20.);
        graphCombInvYieldTot2760GeVALHC11h_0010->Fit(fitTCMInvYield2760GeVas10h_0010,"NRMEX0+","",0.4,20.);

        cout << "Fit to the 20-50% spectrum" << endl;
		Double_t paramTCM_2050[5] = {graphCombInvYieldTot2760GeVALHC11h_2050->GetY()[0],0.3,graphCombInvYieldTot2760GeVALHC11h_2050->GetY()[0],0.3,8};
        fitTCMInvYield2760GeVLHC11h_2050 = FitObject("tcm","fitBylinkin2050","Pi0",graphCombInvYieldTot2760GeVALHC11h_2050,minfitPt,20.,paramTCM_2050,"NRMEX0+");
        graphCombInvYieldTot2760GeVALHC11h_2050->Fit(fitTCMInvYield2760GeVLHC11h_2050,"NRMEX0+","",minfitPt,20.);     

//         fitTCMInvYield2760GeVLHC11h_2050 = FitObject("qcd","fitQCD","Pi0",graphCombInvYieldTot2760GeVALHC11h_2050,0.4,20.);
//         graphCombInvYieldTot2760GeVALHC11h_2050->Fit(fitTCMInvYield2760GeVLHC11h_2050,"NRMEX0+","",0.4,20.);     

        fitTCMInvYield2760GeVas10h_2050 = FitObject("qcd","fitQCD","Pi0",graphCombInvYieldTot2760GeVALHC11h_2050,0.4,20.);
        graphCombInvYieldTot2760GeVALHC11h_2050->Fit(fitTCMInvYield2760GeVas10h_2050,"NRMEX0+","",0.4,20.);     
        
	} else if(meson.CompareTo("Eta")==0){
		cout << "Fit to the 0-10% spectrum" << endl;
		Double_t paramTCM_0010[5] = {graphCombInvYieldTot2760GeVALHC11h_0010->GetY()[0],0.3,graphCombInvYieldTot2760GeVALHC11h_0010->GetY()[0],0.3,8};
		fitTCMInvYield2760GeVLHC11h_0010 = FitObject("tcm","fitBylinkin0010","Eta",graphCombInvYieldTot2760GeVALHC11h_0010,1.,20.,paramTCM_0010,"NRMEX0+");
        graphCombInvYieldTot2760GeVALHC11h_0010->Fit(fitTCMInvYield2760GeVLHC11h_0010,"NRMEX0+","",1.,20.);
        
        fitTCMInvYield2760GeVas10h_0010 = FitObject("qcd","fitQCD","Eta",graphCombInvYieldTot2760GeVALHC11h_0010,1.,20.);
        graphCombInvYieldTot2760GeVALHC11h_0010->Fit(fitTCMInvYield2760GeVas10h_0010,"NRMEX0+","",1.,20.);

		cout << "Fit to the 20-50% spectrum" << endl;
// 		Double_t paramTCM_2050[5] = {1.,0.,1.,0.3,8.};
        Double_t paramTCM_2050[5] = {graphCombInvYieldTot2760GeVALHC11h_2050->GetY()[0],0.3,graphCombInvYieldTot2760GeVALHC11h_2050->GetY()[0],0.3,8};
        fitTCMInvYield2760GeVLHC11h_2050 = FitObject("tcm","fitBylinkin2050","Eta",graphCombInvYieldStat2760GeVALHC11h_2050,1.,20.,paramTCM_2050,"NRMEX0+");
        graphCombInvYieldStat2760GeVALHC11h_2050->Fit(fitTCMInvYield2760GeVLHC11h_2050,"NRMEX0+","",1.,20.);     
        
//         fitTCMInvYield2760GeVLHC11h_2050 = FitObject("qcd","fitQCD","Eta",graphCombInvYieldTot2760GeVALHC11h_2050,1.,20.);
//         graphCombInvYieldTot2760GeVALHC11h_2050->Fit(fitTCMInvYield2760GeVLHC11h_2050,"NRMEX0+","",1.,20.);     

        fitTCMInvYield2760GeVas10h_2050 = FitObject("qcd","fitQCD","Eta",graphCombInvYieldTot2760GeVALHC11h_2050,1.,20.);
        graphCombInvYieldTot2760GeVALHC11h_2050->Fit(fitTCMInvYield2760GeVas10h_2050,"NRMEX0+","",1.,20.);     

    }
	cout << WriteParameterToFile(fitTCMInvYield2760GeVLHC11h_0010)<< endl << endl;
	cout << WriteParameterToFile(fitTCMInvYield2760GeVLHC11h_2050)<< endl << endl;

//     cout << WriteParameterToFile(fitTCMInvYield2760GeVas10h_2050)<< endl << endl;
//     return;

	TGraphAsymmErrors* graphRatioCombCombFitTot2760GeVLHC11h_0010 	= (TGraphAsymmErrors*)graphCombInvYieldTot2760GeVALHC11h_0010->Clone();
	graphRatioCombCombFitTot2760GeVLHC11h_0010 						= CalculateGraphErrRatioToFit(graphRatioCombCombFitTot2760GeVLHC11h_0010, fitTCMInvYield2760GeVLHC11h_0010); 
	TGraphAsymmErrors* graphRatioCombCombFitStat2760GeVLHC11h_0010 	= (TGraphAsymmErrors*)graphCombInvYieldStat2760GeVALHC11h_0010->Clone();
	graphRatioCombCombFitStat2760GeVLHC11h_0010 						= CalculateGraphErrRatioToFit(graphRatioCombCombFitStat2760GeVLHC11h_0010, fitTCMInvYield2760GeVLHC11h_0010); 
	TGraphAsymmErrors* graphRatioCombCombFitSys2760GeVLHC11h_0010 	= (TGraphAsymmErrors*)graphCombInvYieldSys2760GeVALHC11h_0010->Clone();
	graphRatioCombCombFitSys2760GeVLHC11h_0010 						= CalculateGraphErrRatioToFit(graphRatioCombCombFitSys2760GeVLHC11h_0010, fitTCMInvYield2760GeVLHC11h_0010); 

	TGraphAsymmErrors* graphRatioCombCombFitTot2760GeVLHC11h_2050 	= (TGraphAsymmErrors*)graphCombInvYieldTot2760GeVALHC11h_2050->Clone();
	graphRatioCombCombFitTot2760GeVLHC11h_2050 						= CalculateGraphErrRatioToFit(graphRatioCombCombFitTot2760GeVLHC11h_2050, fitTCMInvYield2760GeVLHC11h_2050); 
	TGraphAsymmErrors* graphRatioCombCombFitStat2760GeVLHC11h_2050 	= (TGraphAsymmErrors*)graphCombInvYieldStat2760GeVALHC11h_2050->Clone();
	graphRatioCombCombFitStat2760GeVLHC11h_2050 						= CalculateGraphErrRatioToFit(graphRatioCombCombFitStat2760GeVLHC11h_2050, fitTCMInvYield2760GeVLHC11h_2050); 
	TGraphAsymmErrors* graphRatioCombCombFitSys2760GeVLHC11h_2050 	= (TGraphAsymmErrors*)graphCombInvYieldSys2760GeVALHC11h_2050->Clone();
	graphRatioCombCombFitSys2760GeVLHC11h_2050 						= CalculateGraphErrRatioToFit(graphRatioCombCombFitSys2760GeVLHC11h_2050, fitTCMInvYield2760GeVLHC11h_2050); 

	TGraphAsymmErrors* graphRatioPCMCombFitStat2760GeV_0010;
	TGraphAsymmErrors* graphRatioPCMCombFitSys2760GeV_0010;
	TGraphAsymmErrors* graphRatioEMCalCombFitStat2760GeV_0010;
	TGraphAsymmErrors* graphRatioEMCalCombFitSys2760GeV_0010;
	
	TGraphAsymmErrors* graphRatioPCMCombFitStat2760GeV_2050;
	TGraphAsymmErrors* graphRatioPCMCombFitSys2760GeV_2050;
	TGraphAsymmErrors* graphRatioEMCalCombFitStat2760GeV_2050;
	TGraphAsymmErrors* graphRatioEMCalCombFitSys2760GeV_2050;
	
	if(meson.CompareTo("Pi0")==0){
		
		graphRatioPCMCombFitStat2760GeV_0010		= (TGraphAsymmErrors*)graphPCMPi0InvYieldStatPbPb2760GeV_0010->Clone();
		graphRatioPCMCombFitStat2760GeV_0010						= CalculateGraphErrRatioToFit(graphRatioPCMCombFitStat2760GeV_0010, fitTCMInvYield2760GeVLHC11h_0010); 
		graphRatioPCMCombFitSys2760GeV_0010 		= (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_0010->Clone();
		graphRatioPCMCombFitSys2760GeV_0010 						= CalculateGraphErrRatioToFit(graphRatioPCMCombFitSys2760GeV_0010, fitTCMInvYield2760GeVLHC11h_0010); 
		graphRatioEMCalCombFitStat2760GeV_0010 	= (TGraphAsymmErrors*)graphEMCalPi0InvYieldStatPbPb2760GeV_0010->Clone();
		graphRatioEMCalCombFitStat2760GeV_0010 						= CalculateGraphErrRatioToFit(graphRatioEMCalCombFitStat2760GeV_0010, fitTCMInvYield2760GeVLHC11h_0010); 
		graphRatioEMCalCombFitSys2760GeV_0010 	= (TGraphAsymmErrors*)graphEMCalPi0InvYieldSysPbPb2760GeV_0010->Clone();
		graphRatioEMCalCombFitSys2760GeV_0010 						= CalculateGraphErrRatioToFit(graphRatioEMCalCombFitSys2760GeV_0010, fitTCMInvYield2760GeVLHC11h_0010); 
		
		graphRatioPCMCombFitStat2760GeV_2050		= (TGraphAsymmErrors*)graphPCMPi0InvYieldStatPbPb2760GeV_2050->Clone();
		graphRatioPCMCombFitStat2760GeV_2050						= CalculateGraphErrRatioToFit(graphRatioPCMCombFitStat2760GeV_2050, fitTCMInvYield2760GeVLHC11h_2050); 
		graphRatioPCMCombFitSys2760GeV_2050 		= (TGraphAsymmErrors*)graphPCMPi0InvYieldSysPbPb2760GeV_2050->Clone();
		graphRatioPCMCombFitSys2760GeV_2050 						= CalculateGraphErrRatioToFit(graphRatioPCMCombFitSys2760GeV_2050, fitTCMInvYield2760GeVLHC11h_2050); 
		graphRatioEMCalCombFitStat2760GeV_2050 	= (TGraphAsymmErrors*)graphEMCalPi0InvYieldStatPbPb2760GeV_2050->Clone();
		graphRatioEMCalCombFitStat2760GeV_2050 						= CalculateGraphErrRatioToFit(graphRatioEMCalCombFitStat2760GeV_2050, fitTCMInvYield2760GeVLHC11h_2050); 
		graphRatioEMCalCombFitSys2760GeV_2050 	= (TGraphAsymmErrors*)graphEMCalPi0InvYieldSysPbPb2760GeV_2050->Clone();
		graphRatioEMCalCombFitSys2760GeV_2050 						= CalculateGraphErrRatioToFit(graphRatioEMCalCombFitSys2760GeV_2050, fitTCMInvYield2760GeVLHC11h_2050); 
		
	} else if(meson.CompareTo("Eta")==0){
		
		graphRatioPCMCombFitStat2760GeV_0010		= (TGraphAsymmErrors*)graphPCMEtaInvYieldStatPbPb2760GeV_0010->Clone();
		graphRatioPCMCombFitStat2760GeV_0010						= CalculateGraphErrRatioToFit(graphRatioPCMCombFitStat2760GeV_0010, fitTCMInvYield2760GeVLHC11h_0010); 
		graphRatioPCMCombFitSys2760GeV_0010 		= (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_0010->Clone();
		graphRatioPCMCombFitSys2760GeV_0010 						= CalculateGraphErrRatioToFit(graphRatioPCMCombFitSys2760GeV_0010, fitTCMInvYield2760GeVLHC11h_0010); 
		graphRatioEMCalCombFitStat2760GeV_0010 	= (TGraphAsymmErrors*)graphEMCalEtaInvYieldStatPbPb2760GeV_0010->Clone();
		graphRatioEMCalCombFitStat2760GeV_0010 						= CalculateGraphErrRatioToFit(graphRatioEMCalCombFitStat2760GeV_0010, fitTCMInvYield2760GeVLHC11h_0010); 
		graphRatioEMCalCombFitSys2760GeV_0010 	= (TGraphAsymmErrors*)graphEMCalEtaInvYieldSysPbPb2760GeV_0010->Clone();
		graphRatioEMCalCombFitSys2760GeV_0010 						= CalculateGraphErrRatioToFit(graphRatioEMCalCombFitSys2760GeV_0010, fitTCMInvYield2760GeVLHC11h_0010); 
		
		graphRatioPCMCombFitStat2760GeV_2050		= (TGraphAsymmErrors*)graphPCMEtaInvYieldStatPbPb2760GeV_2050->Clone();
		graphRatioPCMCombFitStat2760GeV_2050						= CalculateGraphErrRatioToFit(graphRatioPCMCombFitStat2760GeV_2050, fitTCMInvYield2760GeVLHC11h_2050); 
		graphRatioPCMCombFitSys2760GeV_2050 		= (TGraphAsymmErrors*)graphPCMEtaInvYieldSysPbPb2760GeV_2050->Clone();
		graphRatioPCMCombFitSys2760GeV_2050 						= CalculateGraphErrRatioToFit(graphRatioPCMCombFitSys2760GeV_2050, fitTCMInvYield2760GeVLHC11h_2050); 
		graphRatioEMCalCombFitStat2760GeV_2050 	= (TGraphAsymmErrors*)graphEMCalEtaInvYieldStatPbPb2760GeV_2050->Clone();
		graphRatioEMCalCombFitStat2760GeV_2050 						= CalculateGraphErrRatioToFit(graphRatioEMCalCombFitStat2760GeV_2050, fitTCMInvYield2760GeVLHC11h_2050); 
		graphRatioEMCalCombFitSys2760GeV_2050 	= (TGraphAsymmErrors*)graphEMCalEtaInvYieldSysPbPb2760GeV_2050->Clone();
		graphRatioEMCalCombFitSys2760GeV_2050 						= CalculateGraphErrRatioToFit(graphRatioEMCalCombFitSys2760GeV_2050, fitTCMInvYield2760GeVLHC11h_2050); 
		
	}
		
	if(PaperPi0 && meson.CompareTo("Pi0")==0){
		
		graphRatioCombCombFitStat2760GeVLHC11h_0010->RemovePoint(0);
		graphRatioCombCombFitStat2760GeVLHC11h_0010->RemovePoint(0);
		graphRatioCombCombFitStat2760GeVLHC11h_0010->RemovePoint(0);
		
		graphRatioCombCombFitSys2760GeVLHC11h_0010->RemovePoint(0);
		graphRatioCombCombFitSys2760GeVLHC11h_0010->RemovePoint(0);
		graphRatioCombCombFitSys2760GeVLHC11h_0010->RemovePoint(0);

		graphRatioCombCombFitStat2760GeVLHC11h_2050->RemovePoint(0);
		graphRatioCombCombFitStat2760GeVLHC11h_2050->RemovePoint(0);
		graphRatioCombCombFitStat2760GeVLHC11h_2050->RemovePoint(0);

		graphRatioCombCombFitSys2760GeVLHC11h_2050->RemovePoint(0);
		graphRatioCombCombFitSys2760GeVLHC11h_2050->RemovePoint(0);
		graphRatioCombCombFitSys2760GeVLHC11h_2050->RemovePoint(0);
	}
	
	
	//**********************************************************************************************************************//
	//*********************************************     Ratio of Comb to Fit    ********************************************//
	//**********************************************************************************************************************//
	
	textSizeLabelsPixel = 48;
	TCanvas* canvasRatioToCombFitLHC11h = new TCanvas("canvasRatioToCombFitLHC11h","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioToCombFitLHC11h, 0.12, 0.01, 0.01, 0.11);
	canvasRatioToCombFitLHC11h->SetLogx();

		Double_t textsizeLabelsPP = 0;
		Double_t textsizeFacPP= 0;
		if (canvasRatioToCombFitLHC11h->XtoPixel(canvasRatioToCombFitLHC11h->GetX2()) <canvasRatioToCombFitLHC11h->YtoPixel(canvasRatioToCombFitLHC11h->GetY1()) ){
			textsizeLabelsPP = (Double_t)textSizeLabelsPixel/canvasRatioToCombFitLHC11h->XtoPixel(canvasRatioToCombFitLHC11h->GetX2()) ;
			textsizeFacPP = (Double_t)1./canvasRatioToCombFitLHC11h->XtoPixel(canvasRatioToCombFitLHC11h->GetX2()) ;
		} else {
			textsizeLabelsPP = (Double_t)textSizeLabelsPixel/canvasRatioToCombFitLHC11h->YtoPixel(canvasRatioToCombFitLHC11h->GetY1());
			textsizeFacPP = (Double_t)1./canvasRatioToCombFitLHC11h->YtoPixel(canvasRatioToCombFitLHC11h->GetY1());
		}
		cout << textsizeLabelsPP << endl;
	
		TH2F * histo2DRatioToCombFitLHC11h;
		histo2DRatioToCombFitLHC11h = new TH2F("histo2DRatioToCombFitLHC11h","histo2DRatioToCombFitLHC11h",1000,0.23,70.,1000,0.2,4.	);
		SetStyleHistoTH2ForGraphs(histo2DRatioToCombFitLHC11h, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{Comb Fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP, 
								0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
		histo2DRatioToCombFitLHC11h->GetXaxis()->SetMoreLogLabels();
		histo2DRatioToCombFitLHC11h->GetXaxis()->SetLabelOffset(-0.01);
	// 	histo2DRatioToCombFitLHC11h->GetYaxis()->SetRangeUser(-10,10);
		histo2DRatioToCombFitLHC11h->GetYaxis()->SetRangeUser(0.05,2.45);
		histo2DRatioToCombFitLHC11h->Draw("copy");

		DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys2760GeVLHC11h_0010, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
		graphRatioCombCombFitSys2760GeVLHC11h_0010->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat2760GeVLHC11h_0010, markerStyleComb, markerSizeComb, colorComb , colorComb);
		graphRatioCombCombFitStat2760GeVLHC11h_0010->Draw("p,same,e0");

		DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys2760GeVLHC11h_2050, markerStyleComb, markerSizeComb, colorComb+1, colorComb+1, widthLinesBoxes, kTRUE);
		graphRatioCombCombFitSys2760GeVLHC11h_2050->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat2760GeVLHC11h_2050, markerStyleComb, markerSizeComb, colorComb+1, colorComb+1);
		graphRatioCombCombFitStat2760GeVLHC11h_2050->Draw("p,same,e0");

		TLegend* legendCombToCombFit = new TLegend(0.2,0.75,0.5,0.94);
		legendCombToCombFit->SetFillColor(0);
		legendCombToCombFit->SetLineColor(0);
		legendCombToCombFit->SetTextFont(42);
		legendCombToCombFit->SetTextSize(0.04);
		legendCombToCombFit->AddEntry(graphRatioCombCombFitSys2760GeVLHC11h_0010,"0-10%","fp");
		legendCombToCombFit->AddEntry(graphRatioCombCombFitSys2760GeVLHC11h_2050,"20-50%","fp");
		legendCombToCombFit->Draw();

		DrawGammaLines(0.23, 70. , 1., 1.,0.1, kGray+2);
		DrawGammaLines(0.23, 70. , 1.1, 1.1,0.1, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.9, 0.9,0.1, kGray, 7);

		TLatex *labelRatioToFitEnergy = new TLatex(0.65,0.92,collisionSystem2760GeV.Data());
		SetStyleTLatex( labelRatioToFitEnergy, 0.85*textSizeLabelsPixel,4);
		labelRatioToFitEnergy->SetTextFont(43);
		labelRatioToFitEnergy->Draw();
		TLatex *labelRatioToFitPi0;
		if(meson.CompareTo("Pi0")==0){
			labelRatioToFitPi0= new TLatex(0.73,0.87,"#pi^{0} #rightarrow #gamma#gamma");
		} else if(meson.CompareTo("Eta")==0){
			labelRatioToFitPi0= new TLatex(0.73,0.87,"#eta #rightarrow #gamma#gamma");
		}
		SetStyleTLatex( labelRatioToFitPi0, 0.85*textSizeLabelsPixel,4);
		labelRatioToFitPi0->SetTextFont(43);
		labelRatioToFitPi0->Draw();

	canvasRatioToCombFitLHC11h->SaveAs(Form("%s/%s_RatioOfCombToCombFit_PbPb2760GeVLHC11h.%s",outputDir.Data(),meson.Data(),suffix.Data()));
// 	canvasRatioToCombFitLHC11h->SaveAs(Form("%s/%s_RatioOfCombToCombFit_PbPb2760GeVLHC11h.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
	
	//**********************************************************************************************************************//
	//******************************************* Ratio of Individual meas to Fit ******************************************//
	//**********************************************************************************************************************//
	
	canvasRatioToCombFitLHC11h->cd();
	histo2DRatioToCombFitLHC11h->Draw("copy");

		DrawGammaSetMarkerTGraphAsym(graphRatioPCMCombFitSys2760GeV_0010, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorPCM0010, colorPCM0010, widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphRatioPCMCombFitStat2760GeV_0010, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorPCM0010, colorPCM0010);
		DrawGammaSetMarkerTGraphAsym(graphRatioEMCalCombFitSys2760GeV_0010, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorEMCal0010, colorEMCal0010, widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphRatioEMCalCombFitStat2760GeV_0010, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorEMCal0010, colorEMCal0010);
		
		graphRatioPCMCombFitSys2760GeV_0010->Draw("E2same");
		graphRatioEMCalCombFitSys2760GeV_0010->Draw("E2same");
		
		graphRatioPCMCombFitStat2760GeV_0010->Draw("p,same,e");
		graphRatioEMCalCombFitStat2760GeV_0010->Draw("p,same,e");

		DrawGammaSetMarkerTGraphAsym(graphRatioPCMCombFitSys2760GeV_2050, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorPCM2050, colorPCM2050, widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphRatioPCMCombFitStat2760GeV_2050, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorPCM2050, colorPCM2050);
		DrawGammaSetMarkerTGraphAsym(graphRatioEMCalCombFitSys2760GeV_2050, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorEMCal2050, colorEMCal2050, widthLinesBoxes, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphRatioEMCalCombFitStat2760GeV_2050, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorEMCal2050, colorEMCal2050);
		
		graphRatioPCMCombFitSys2760GeV_2050->Draw("E2same");
		graphRatioEMCalCombFitSys2760GeV_2050->Draw("E2same");
		
		graphRatioPCMCombFitStat2760GeV_2050->Draw("p,same,e");
		graphRatioEMCalCombFitStat2760GeV_2050->Draw("p,same,e");

		DrawGammaLines(0.23, 70. , 1., 1.,0.5, kGray+2);
		DrawGammaLines(0.23, 70. , 1.1, 1.1,0.5, kGray, 7);
		DrawGammaLines(0.23, 70. , 0.9, 0.9,0.5, kGray, 7);
		
		labelRatioToFitEnergy->Draw();
		labelRatioToFitPi0->Draw();
	
		// ****************************** Definition of the Legend ******************************************
		// **************** Row def ************************
		Double_t rowsLegendOnlyPi0Ratio[5] 		= {0.92,0.88,0.84,0.80,0.76};
		Double_t rowsLegendOnlyPi0RatioAbs[5] 	= {0.91,2.2,2.1,2.0,1.9};
		Double_t columnsLegendOnlyPi0Ratio[3] 	= {0.15,0.32, 0.38};
		Double_t columnsLegendOnlyPi0RatioAbs[3]= {0.15,1.04, 1.37};
		Double_t lengthBox						= 0.2/2;
		Double_t heightBox						= 0.08/2;
		// ****************** first Column **************************************************
		TLatex *textPCMOnlyRatioPi0LHC11h = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[1],"PCM");
		SetStyleTLatex( textPCMOnlyRatioPi0LHC11h, 0.85*textSizeLabelsPixel,4);
		textPCMOnlyRatioPi0LHC11h->SetTextFont(43);
		textPCMOnlyRatioPi0LHC11h->Draw();
		TLatex *textEMCalOnlyRatioPi0LHC11h = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[3],"EMCal");
		SetStyleTLatex( textEMCalOnlyRatioPi0LHC11h,  0.85*textSizeLabelsPixel,4);
		textEMCalOnlyRatioPi0LHC11h->SetTextFont(43);
		textEMCalOnlyRatioPi0LHC11h->Draw();
		
		// ****************** second Column *************************************************
		TLatex *textStatOnlyRatioPi0LHC11h = new TLatex(columnsLegendOnlyPi0Ratio[1],rowsLegendOnlyPi0Ratio[0] ,"stat");
		SetStyleTLatex( textStatOnlyRatioPi0LHC11h, 0.85*textSizeLabelsPixel,4);
		textStatOnlyRatioPi0LHC11h->SetTextFont(43);
		textStatOnlyRatioPi0LHC11h->Draw();
		TLatex *textSysOnlyRatioPi0LHC11h = new TLatex(columnsLegendOnlyPi0Ratio[2] ,rowsLegendOnlyPi0Ratio[0],"syst");
		SetStyleTLatex( textSysOnlyRatioPi0LHC11h, 0.85*textSizeLabelsPixel,4);
		textSysOnlyRatioPi0LHC11h->SetTextFont(43);
		textSysOnlyRatioPi0LHC11h->Draw();
		
		TMarker* markerPCMPi0OnlyRatioPi0LHC11h = CreateMarkerFromGraph(graphRatioPCMCombFitSys2760GeV_0010,columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[1],1);
		markerPCMPi0OnlyRatioPi0LHC11h->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[1]);
		TMarker* markerEMCalPi0OnlyRatioPi0LHC11h = CreateMarkerFromGraph(graphRatioEMCalCombFitSys2760GeV_0010, columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[3],1);
		markerEMCalPi0OnlyRatioPi0LHC11h->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[3]);

		TBox* boxPCMPi0OnlyRatioPi0 = CreateBoxFromGraph(graphRatioPCMCombFitSys2760GeV_0010, columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[1]- heightBox,
														columnsLegendOnlyPi0RatioAbs[2]+ 3*lengthBox, rowsLegendOnlyPi0RatioAbs[1]+ heightBox);
		boxPCMPi0OnlyRatioPi0->Draw("l");
		TBox* boxEMCalPi0OnlyRatioPi0 = CreateBoxFromGraph(graphRatioEMCalCombFitSys2760GeV_0010, columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[3]- heightBox,
														columnsLegendOnlyPi0RatioAbs[2]+ 3*lengthBox, rowsLegendOnlyPi0RatioAbs[3]+ heightBox);
		
		TMarker* markerPCMPi0OnlyRatioPi0LHC11h_2050 = CreateMarkerFromGraph(graphRatioPCMCombFitSys2760GeV_2050,columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[1],1);
		markerPCMPi0OnlyRatioPi0LHC11h_2050->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[1]);
		TMarker* markerEMCalPi0OnlyRatioPi0LHC11h_2050 = CreateMarkerFromGraph(graphRatioEMCalCombFitSys2760GeV_2050, columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[3],1);
		markerEMCalPi0OnlyRatioPi0LHC11h_2050->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[3]);

		TBox* boxPCMPi0OnlyRatioPi0LHC11h_2050 = CreateBoxFromGraph(graphRatioPCMCombFitSys2760GeV_2050, columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[1]- heightBox,
														columnsLegendOnlyPi0RatioAbs[2]+ 3*lengthBox, rowsLegendOnlyPi0RatioAbs[1]+ heightBox);
		boxPCMPi0OnlyRatioPi0LHC11h_2050->Draw("l");
		TBox* boxEMCalPi0OnlyRatioPi0LHC11h_2050 = CreateBoxFromGraph(graphRatioEMCalCombFitSys2760GeV_2050, columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[3]- heightBox,
														columnsLegendOnlyPi0RatioAbs[2]+ 3*lengthBox, rowsLegendOnlyPi0RatioAbs[3]+ heightBox);

		boxEMCalPi0OnlyRatioPi0LHC11h_2050->Draw("l");
	
	canvasRatioToCombFitLHC11h->SaveAs(Form("%s/%s_RatioOfIndividualMeasToCombFit_PbPb2760GeVLHC11h.%s",outputDir.Data(),meson.Data(),suffix.Data()));

	
	
	//===================================================================================================================================================//
	//                                                           Drawing main plots                                                                      //
	//===================================================================================================================================================//
	
	//**********************************************************************************************************************//
	//***************************************** Plotting Combined Invariant Yields *****************************************//
	//**********************************************************************************************************************//
	
	TCanvas* canvasInvYieldSectionPi0LHC11h = new TCanvas("canvasInvYieldSectionPi0LHC11h","",200,10,1350,1350*1.15);  // gives the page size
	DrawGammaCanvasSettings( canvasInvYieldSectionPi0LHC11h, 0.16, 0.02, 0.02, 0.09);
	canvasInvYieldSectionPi0LHC11h->SetLogx();
	canvasInvYieldSectionPi0LHC11h->SetLogy();
	
		TH2F * histo2DInvYieldSectionPi0LHC11h;
		if(meson.CompareTo("Pi0")==0){
			histo2DInvYieldSectionPi0LHC11h = new TH2F("histo2DInvYieldSectionPi0LHC11h","histo2DInvYieldSectionPi0LHC11h",11000,minPtYields,maxPtYields,1000,minYaxisYields,maxYaxisYields); 
			SetStyleHistoTH2ForGraphs(histo2DInvYieldSectionPi0LHC11h, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#pi^{0}}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",0.035,0.04, 0.035,0.04, 1.,1.5);
			
		} else if(meson.CompareTo("Eta")==0){
			histo2DInvYieldSectionPi0LHC11h = new TH2F("histo2DInvYieldSectionPi0LHC11h","histo2DInvYieldSectionPi0LHC11h",11000,minPtYields,maxPtYields,1000,minYaxisYields,maxYaxisYields);
			SetStyleHistoTH2ForGraphs(histo2DInvYieldSectionPi0LHC11h, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#eta}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",0.035,0.04, 0.035,0.04, 1.,1.6);			
		}
		histo2DInvYieldSectionPi0LHC11h->GetXaxis()->SetMoreLogLabels();
		histo2DInvYieldSectionPi0LHC11h->GetXaxis()->SetLabelOffset(-0.01);
		histo2DInvYieldSectionPi0LHC11h->Draw("copy");
// 		histo2DInvYieldSectionPi0LHC11h->GetXaxis()->SetRangeUser(minPtRange,maxPtRange);
		
		if(PaperPi0 && meson.CompareTo("Pi0")==0){
			graphCombInvYieldSys2760GeVALHC11h_0010->RemovePoint(0);
			graphCombInvYieldStat2760GeVALHC11h_0010->RemovePoint(0);

			graphCombInvYieldSys2760GeVALHC11h_0010->RemovePoint(0);
			graphCombInvYieldStat2760GeVALHC11h_0010->RemovePoint(0);

			graphCombInvYieldSys2760GeVALHC11h_0010->RemovePoint(0);
			graphCombInvYieldStat2760GeVALHC11h_0010->RemovePoint(0);

			graphCombInvYieldSys2760GeVALHC11h_2050->RemovePoint(0);
			graphCombInvYieldStat2760GeVALHC11h_2050->RemovePoint(0);

			graphCombInvYieldSys2760GeVALHC11h_2050->RemovePoint(0);
			graphCombInvYieldStat2760GeVALHC11h_2050->RemovePoint(0);
			
			graphCombInvYieldSys2760GeVALHC11h_2050->RemovePoint(0);
			graphCombInvYieldStat2760GeVALHC11h_2050->RemovePoint(0);

		}
		
// 		DrawGammaSetMarkerTGraphAsym(graphCombInvYieldSys2760GeVALHC11h_0010, markerStyle0010, markerSizeComb, colorCombo0010 , colorCombo0010, widthLinesBoxes, kTRUE);
// 		DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStat2760GeVALHC11h_0010, markerStyle0010, markerSizeComb, colorCombo0010 , colorCombo0010);		
		TGraphAsymmErrors *graphCombInvYieldSys2760GeVALHC11hScaled_0010 = (TGraphAsymmErrors*)graphCombInvYieldSys2760GeVALHC11h_0010->Clone("graphCombInvYieldSys2760GeVALHC11hScaled_0010");
		graphCombInvYieldSys2760GeVALHC11hScaled_0010 = ScaleGraph(graphCombInvYieldSys2760GeVALHC11hScaled_0010, 4); 
		TGraphAsymmErrors *graphCombInvYieldStat2760GeVALHC11hScaled_0010 = (TGraphAsymmErrors*)graphCombInvYieldStat2760GeVALHC11h_0010->Clone("graphCombInvYieldStat2760GeVALHC11hScaled_0010");
		graphCombInvYieldStat2760GeVALHC11hScaled_0010 = ScaleGraph(graphCombInvYieldStat2760GeVALHC11hScaled_0010,4);
        
		DrawGammaSetMarkerTGraphAsym(graphCombInvYieldSys2760GeVALHC11hScaled_0010, markerStyle0010, markerSizeComb, colorCombo0010 , colorCombo0010, 2, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStat2760GeVALHC11hScaled_0010, markerStyle0010, markerSizeComb, colorCombo0010, colorCombo0010);
		graphCombInvYieldSys2760GeVALHC11hScaled_0010->Draw("E2same");
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombInvYieldStat2760GeVALHC11hScaled_0010);
		graphCombInvYieldStat2760GeVALHC11hScaled_0010->Draw("p,same,e0");	
        
		DrawGammaSetMarkerTGraphAsym(graphCombInvYieldSys2760GeVALHC11h_2050, markerStyle2050, markerSizeComb, colorCombo2050 ,colorCombo2050, 2, kTRUE);
		DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStat2760GeVALHC11h_2050, markerStyle2050, markerSizeComb,  colorCombo2050 ,colorCombo2050);
		graphCombInvYieldSys2760GeVALHC11h_2050->Draw("E2same");
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombInvYieldStat2760GeVALHC11h_2050);
		graphCombInvYieldStat2760GeVALHC11h_2050->Draw("p,same,e0");	
			
		TLegend* legendInvYieldSectionPi0LHC11h_onlyPbPb = new TLegend(0.2,0.13,0.53,0.23);
		legendInvYieldSectionPi0LHC11h_onlyPbPb->SetMargin(0.17);
		legendInvYieldSectionPi0LHC11h_onlyPbPb->SetFillColor(0);
		legendInvYieldSectionPi0LHC11h_onlyPbPb->SetLineColor(0);
		legendInvYieldSectionPi0LHC11h_onlyPbPb->SetTextFont(42);
		legendInvYieldSectionPi0LHC11h_onlyPbPb->SetTextSize(FontSize);
		legendInvYieldSectionPi0LHC11h_onlyPbPb->SetHeader(collisionSystem2760GeV.Data());
		legendInvYieldSectionPi0LHC11h_onlyPbPb->AddEntry(graphCombInvYieldSys2760GeVALHC11hScaled_0010,Form("  %s",cent0010.Data()),"fp");
		legendInvYieldSectionPi0LHC11h_onlyPbPb->AddEntry(graphCombInvYieldSys2760GeVALHC11h_2050,Form("%s",cent2050.Data()),"fp");
		legendInvYieldSectionPi0LHC11h_onlyPbPb->Draw();
		
		labelSystOnlyPbPb->Draw();
		labelFactorLowerOnlyPbPb->Draw();

	canvasInvYieldSectionPi0LHC11h->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataOnly.%s",outputDir.Data(),meson.Data(),suffix.Data()));
	canvasInvYieldSectionPi0LHC11h->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataOnly.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
	
	canvasInvYieldSectionPi0LHC11h->cd();
      TH2F * histo2DInvYieldSectionLHC11hwithPP;
      if(meson.CompareTo("Pi0")==0){
          histo2DInvYieldSectionLHC11hwithPP = new TH2F("histo2DInvYieldSectionLHC11hwithPP","histo2DInvYieldSectionLHC11hwithPP",11000,0.25,70.,1000,1e-12,1e3); 
          SetStyleHistoTH2ForGraphs(histo2DInvYieldSectionLHC11hwithPP, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#pi^{0}}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2} ",0.035,0.04, 0.035,0.04, 1.,1.65);
          
      } else if(meson.CompareTo("Eta")==0){
          histo2DInvYieldSectionLHC11hwithPP = new TH2F("histo2DInvYieldSectionLHC11hwithPP","histo2DInvYieldSectionLHC11hwithPP",11000,0.3,70.,1000,2e-10,1e3);
          SetStyleHistoTH2ForGraphs(histo2DInvYieldSectionLHC11hwithPP, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#eta}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2} ",0.035,0.04, 0.035,0.04, 1.,1.6);            
      }
      histo2DInvYieldSectionLHC11hwithPP->GetXaxis()->SetMoreLogLabels();
      histo2DInvYieldSectionLHC11hwithPP->GetXaxis()->SetLabelOffset(-0.01);
      histo2DInvYieldSectionLHC11hwithPP->Draw("copy");

		TLegend* legendSpectraPP = new TLegend(0.2,0.15,0.54,0.22);
		legendSpectraPP->SetFillColor(0);
		legendSpectraPP->SetLineColor(0);
		legendSpectraPP->SetTextFont(42);
		legendSpectraPP->SetMargin(0.17);
		legendSpectraPP->SetTextSize(FontSize);

		if(meson.CompareTo("Pi0")==0){
			
			DrawGammaSetMarkerTGraphAsym(graphInvSectionCombStatPi02760GeVPlot, markerStylepp,markerSizepp, kBlack , kBlack);
			if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphInvSectionCombStatPi02760GeVPlot);
			graphInvSectionCombStatPi02760GeVPlot->Draw("p,same,e0");
			DrawGammaSetMarkerTGraphAsym(graphInvSectionCombSysPi02760GeVPlot, markerStylepp,markerSizepp, kBlack , kBlack, widthLinesBoxes, kTRUE);
			graphInvSectionCombSysPi02760GeVPlot->Draw("E2same");
			
			legendSpectraPP->SetHeader(collisionSystemPP2760GeV.Data());
			legendSpectraPP->AddEntry(graphInvSectionCombSysPi02760GeVPlot,"arXiv:XXXX.XXXX","pf");

		} else if(meson.CompareTo("Eta")==0){

			DrawGammaSetMarkerTGraphAsym(graphInvSectionCombStatEta2760GeVPlot, markerStylepp,markerSizepp, kBlack , kBlack);
			if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphInvSectionCombStatEta2760GeVPlot);
			graphInvSectionCombStatEta2760GeVPlot->Draw("p,same,e0");
			DrawGammaSetMarkerTGraphAsym(graphInvSectionCombSysEta2760GeVPlot, markerStylepp,markerSizepp, kBlack , kBlack, widthLinesBoxes, kTRUE);
			graphInvSectionCombSysEta2760GeVPlot->Draw("E2same");
			
			legendSpectraPP->SetHeader(collisionSystemPP2760GeV.Data());
			legendSpectraPP->AddEntry(graphInvSectionCombSysEta2760GeVPlot,"arXiv:XXXX.XXXX","pf"); 

		}
		legendSpectraPP->Draw();
		
		TLegend* legendInvYieldSectionPi0LHC11h_WitPP = new TLegend(0.595,0.79,0.91,0.92); //0.17,0.13,0.5,0.24);
		legendInvYieldSectionPi0LHC11h_WitPP->SetFillColor(0);
		legendInvYieldSectionPi0LHC11h_WitPP->SetMargin(0.17);
		legendInvYieldSectionPi0LHC11h_WitPP->SetLineColor(0);
		legendInvYieldSectionPi0LHC11h_WitPP->SetTextFont(42);
		legendInvYieldSectionPi0LHC11h_WitPP->SetTextSize(FontSize);
		legendInvYieldSectionPi0LHC11h_WitPP->SetHeader(collisionSystem2760GeV.Data());
		legendInvYieldSectionPi0LHC11h_WitPP->AddEntry(graphCombInvYieldSys2760GeVALHC11hScaled_0010,Form("  %s",cent0010.Data()),"fp");
		legendInvYieldSectionPi0LHC11h_WitPP->AddEntry(graphCombInvYieldSys2760GeVALHC11h_2050,Form("%s",cent2050.Data()),"fp");
        legendInvYieldSectionPi0LHC11h_WitPP->AddEntry(fitTCMInvYield2760GeVLHC11h_0010,"fits to Pb#font[122]{-}Pb","l");
//         legendInvYieldSectionPi0LHC11h_WitPP->AddEntry((TObject*)0,"#it{A}_{e} exp(-#it{E}_{T, kin}/#it{T}_{e}) + #it{A}/#(){1 + #frac{#it{p}_{T}^{2}}{#it{T}^{2}#upoint n}}^{-n}","");
		legendInvYieldSectionPi0LHC11h_WitPP->Draw();

        fitTCMInvYield2760GeVLHC11h_0010->SetParameter(0,fitTCMInvYield2760GeVLHC11h_0010->GetParameter(0)*4);
        fitTCMInvYield2760GeVLHC11h_0010->SetParameter(2,fitTCMInvYield2760GeVLHC11h_0010->GetParameter(2)*4);
        fitTCMInvYield2760GeVLHC11h_0010->SetLineStyle(4);
        fitTCMInvYield2760GeVLHC11h_0010->Draw("same");     
        fitTCMInvYield2760GeVLHC11h_2050->SetLineStyle(4);
        fitTCMInvYield2760GeVLHC11h_2050->Draw("same");     


        fitTCMInvYield2760GeVas10h_0010->SetLineColor(kBlue);
        fitTCMInvYield2760GeVas10h_0010->SetParameter(0,fitTCMInvYield2760GeVas10h_0010->GetParameter(0)*4);
//         fitTCMInvYield2760GeVas10h_0010->SetParameter(2,fitTCMInvYield2760GeVLHC11h_0010->GetParameter(2)*4);
        fitTCMInvYield2760GeVas10h_0010->SetLineStyle(2);
//         fitTCMInvYield2760GeVas10h_0010->Draw("same");  
        fitTCMInvYield2760GeVas10h_2050->SetLineColor(kBlue);
        fitTCMInvYield2760GeVas10h_2050->SetLineStyle(2);
//         fitTCMInvYield2760GeVas10h_2050->Draw("same");     

		labelSyst->Draw();
		labelFactorLower->Draw();
// 		labelPreliminary->Draw();

		graphCombInvYieldSys2760GeVALHC11hScaled_0010->Draw("E2same");
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombInvYieldStat2760GeVALHC11hScaled_0010);
		graphCombInvYieldStat2760GeVALHC11hScaled_0010->Draw("p,same,e0");	
	
	
		graphCombInvYieldSys2760GeVALHC11h_2050->Draw("E2same");
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombInvYieldStat2760GeVALHC11h_2050);
		graphCombInvYieldStat2760GeVALHC11h_2050->Draw("p,same,e0");	
	
		
	canvasInvYieldSectionPi0LHC11h->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataOnly_withPP.%s",outputDir.Data(),meson.Data(),suffix.Data()));
	canvasInvYieldSectionPi0LHC11h->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataOnly_withPP.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
	


	
	canvasInvYieldSectionPi0LHC11h->cd();

		TH2F * histo2DInvYieldSectionPi0LHC11hcopy = new TH2F("histo2DInvYieldSectionPi0LHC11hcopy","histo2DInvYieldSectionPi0LHC11hcopy",11000,0.6,40.,1000,1e-8,1e5); 
		SetStyleHistoTH2ForGraphs(histo2DInvYieldSectionPi0LHC11hcopy, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#pi^{0}, #eta}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",0.035,0.04, 0.035,0.04, 1.,1.5);
		histo2DInvYieldSectionPi0LHC11hcopy->DrawCopy();
		
		TFile *MesonInput = new TFile(Form("%s/CombinedResultsPaperPbPb2760GeVLHC11h_%s.root", outputDir.Data(),dateForOutput.Data()));
		if(MesonInput){
	
			TGraphAsymmErrors *graphInvYieldPi0StatBothMeson_0010 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldPi0Comb2760GeVLHC11hStatErr_0010");
			TGraphAsymmErrors *graphInvYieldPi0SysBothMeson_0010= (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldPi0Comb2760GeVLHC11hSysErr_0010");

			TGraphAsymmErrors *graphInvYieldPi0StatBothMeson_2050 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldPi0Comb2760GeVLHC11hStatErr_2050");
			TGraphAsymmErrors *graphInvYieldPi0SysBothMeson_2050 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldPi0Comb2760GeVLHC11hSysErr_2050");

			TGraphAsymmErrors *graphInvYieldEtaStatBothMeson_0010 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldEtaComb2760GeVLHC11hStatErr_0010");
			TGraphAsymmErrors *graphInvYieldEtaSysBothMeson_0010= (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldEtaComb2760GeVLHC11hSysErr_0010");

			TGraphAsymmErrors *graphInvYieldEtaStatBothMeson_2050 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldEtaComb2760GeVLHC11hStatErr_2050");
			TGraphAsymmErrors *graphInvYieldEtaSysBothMeson_2050 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldEtaComb2760GeVLHC11hSysErr_2050");

			TH1D *histoEPOSPi0Rebin_0010 = (TH1D*)MesonInput->Get("EPOSpredictionPi0_0010");
			TH1D *histoEPOSPi0Rebin_2050 = (TH1D*)MesonInput->Get("EPOSpredictionPi0_2050");
			TH1D *histoEPOSEtaRebin_0010 = (TH1D*)MesonInput->Get("EPOSpredictionEta_0010");
			TH1D *histoEPOSEtaRebin_2050 = (TH1D*)MesonInput->Get("EPOSpredictionEta_2050");

			if(graphInvYieldPi0StatBothMeson_0010 && graphInvYieldEtaStatBothMeson_0010){

				TGraphAsymmErrors *TheoryCracowEtaLowPt_scaledx4_0010 = (TGraphAsymmErrors*)TheoryCracowEtaLowPt_0010->Clone();
				TheoryCracowEtaLowPt_scaledx4_0010 = ScaleGraph(TheoryCracowEtaLowPt_scaledx4_0010,4); 
				DrawGammaSetMarkerTGraphAsym(TheoryCracowEtaLowPt_scaledx4_0010, 0, 0, colorCracow0010,colorCracow0010, 5, kTRUE, colorCracow0010);
				TheoryCracowEtaLowPt_scaledx4_0010->Draw("l,same");

				DrawGammaSetMarkerTGraphAsym(TheoryCracowEtaLowPt_2050, 0, 0, colorCracow2050,colorCracow2050, 5, kTRUE,colorCracow2050);
				TheoryCracowEtaLowPt_2050->Draw("l,same");

				TGraphAsymmErrors *TheoryCracowPi0LowPt_scaledx104_0010 = (TGraphAsymmErrors*)TheoryCracowPi0LowPt_0010->Clone();
				TheoryCracowPi0LowPt_scaledx104_0010 = ScaleGraph(TheoryCracowPi0LowPt_scaledx104_0010,400); 
				DrawGammaSetMarkerTGraphAsym(TheoryCracowPi0LowPt_scaledx104_0010, 0, 0, colorCracow0010,colorCracow0010, 5, kTRUE, colorCracow0010);
				TheoryCracowPi0LowPt_scaledx104_0010->Draw("l,same");

				TGraphAsymmErrors *TheoryCracowPi0LowPt_scaledx100_2050 = (TGraphAsymmErrors*)TheoryCracowPi0LowPt_2050->Clone();
				TheoryCracowPi0LowPt_scaledx100_2050 = ScaleGraph(TheoryCracowPi0LowPt_scaledx100_2050,100); 
				DrawGammaSetMarkerTGraphAsym(TheoryCracowPi0LowPt_scaledx100_2050, 0, 0, colorCracow2050,colorCracow2050, 5, kTRUE,colorCracow2050);
				TheoryCracowPi0LowPt_scaledx100_2050->Draw("l,same");
				
				histoEPOSPi0Rebin_0010->Scale(400);
				DrawGammaSetMarker(histoEPOSPi0Rebin_0010,5,2,colorEPOS0010,colorEPOS0010);
				histoEPOSPi0Rebin_0010->SetLineStyle(6);
				histoEPOSPi0Rebin_0010->SetLineWidth(2);
				histoEPOSPi0Rebin_0010->GetXaxis()->SetRangeUser(0.,13.);
				histoEPOSPi0Rebin_0010->Draw("same,c,histo");

				histoEPOSPi0Rebin_2050->Scale(100);
				DrawGammaSetMarker(histoEPOSPi0Rebin_2050,5,2,colorEPOS2050,colorEPOS2050);
				histoEPOSPi0Rebin_2050->SetLineStyle(6);
				histoEPOSPi0Rebin_2050->SetLineWidth(2);
				histoEPOSPi0Rebin_2050->GetXaxis()->SetRangeUser(0.,13.);
				histoEPOSPi0Rebin_2050->Draw("same,c,histo");

				histoEPOSEtaRebin_0010->Scale(4);
				DrawGammaSetMarker(histoEPOSEtaRebin_0010,5,2,colorEPOS0010,colorEPOS0010);
				histoEPOSEtaRebin_0010->SetLineStyle(6);
				histoEPOSEtaRebin_0010->SetLineWidth(2);
				histoEPOSEtaRebin_0010->GetXaxis()->SetRangeUser(0.,13.);
				histoEPOSEtaRebin_0010->Draw("same,c,histo");

				DrawGammaSetMarker(histoEPOSEtaRebin_2050,5,2,colorEPOS2050,colorEPOS2050);
				histoEPOSEtaRebin_2050->SetLineStyle(6);
				histoEPOSEtaRebin_2050->SetLineWidth(2);
				histoEPOSEtaRebin_2050->GetXaxis()->SetRangeUser(0.,13.);
				histoEPOSEtaRebin_2050->Draw("same,c,histo");
				
		
				TGraphAsymmErrors *graphInvYieldPi0SysBothMeson_scaledx104_0010 = (TGraphAsymmErrors*)graphInvYieldPi0SysBothMeson_0010->Clone("");
				graphInvYieldPi0SysBothMeson_scaledx104_0010 = ScaleGraph(graphInvYieldPi0SysBothMeson_scaledx104_0010,400); 
				DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0SysBothMeson_scaledx104_0010, 24, markerSizeComb, colorCombo0010 , colorCombo0010, 2, kTRUE);
				graphInvYieldPi0SysBothMeson_scaledx104_0010->Draw("E2same");
				
				TGraphAsymmErrors *graphInvYieldPi0StatBothMeson_scaledx104_0010 = (TGraphAsymmErrors*)graphInvYieldPi0StatBothMeson_0010->Clone("");
				graphInvYieldPi0StatBothMeson_scaledx104_0010 = ScaleGraph(graphInvYieldPi0StatBothMeson_scaledx104_0010,400);
				DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0StatBothMeson_scaledx104_0010, 24, markerSizeComb, colorCombo0010, colorCombo0010);
				if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphInvYieldPi0StatBothMeson_scaledx104_0010);
				graphInvYieldPi0StatBothMeson_scaledx104_0010->Draw("p,same,e0");	
				
				TGraphAsymmErrors *graphInvYieldPi0SysBothMeson_scaledx100_2050 = (TGraphAsymmErrors*)graphInvYieldPi0SysBothMeson_2050->Clone("");
				graphInvYieldPi0SysBothMeson_scaledx100_2050 = ScaleGraph(graphInvYieldPi0SysBothMeson_scaledx100_2050,100); 
				DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0SysBothMeson_scaledx100_2050, 25, markerSizeComb, colorCombo2050 ,colorCombo2050, 2, kTRUE);
				graphInvYieldPi0SysBothMeson_scaledx100_2050->Draw("E2same");
				TGraphAsymmErrors *graphInvYieldPi0StatBothMeson_scaledx100_2050 = (TGraphAsymmErrors*)graphInvYieldPi0StatBothMeson_2050->Clone("");
				graphInvYieldPi0StatBothMeson_scaledx100_2050 = ScaleGraph(graphInvYieldPi0StatBothMeson_scaledx100_2050,100);
				DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0StatBothMeson_scaledx100_2050, 25, markerSizeComb,  colorCombo2050 ,colorCombo2050);
				if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphInvYieldPi0StatBothMeson_scaledx100_2050);
				graphInvYieldPi0StatBothMeson_scaledx100_2050->Draw("p,same,e0");	

				TGraphAsymmErrors *graphInvYieldEtaSysBothMeson_scaledx4_0010 = (TGraphAsymmErrors*)graphInvYieldEtaSysBothMeson_0010->Clone("");
				graphInvYieldEtaSysBothMeson_scaledx4_0010 = ScaleGraph(graphInvYieldEtaSysBothMeson_scaledx4_0010,4); 
				DrawGammaSetMarkerTGraphAsym(graphInvYieldEtaSysBothMeson_scaledx4_0010, markerStyle0010, markerSizeComb, colorCombo0010 , colorCombo0010, 2, kTRUE);
				graphInvYieldEtaSysBothMeson_scaledx4_0010->Draw("E2same");
				
				TGraphAsymmErrors *graphInvYieldEtaStatBothMeson_scaledx4_0010 = (TGraphAsymmErrors*)graphInvYieldEtaStatBothMeson_0010->Clone("");
				graphInvYieldEtaStatBothMeson_scaledx4_0010 = ScaleGraph(graphInvYieldEtaStatBothMeson_scaledx4_0010,4);
				DrawGammaSetMarkerTGraphAsym(graphInvYieldEtaStatBothMeson_scaledx4_0010, markerStyle0010, markerSizeComb, colorCombo0010, colorCombo0010);
				if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphInvYieldEtaStatBothMeson_scaledx4_0010);
				graphInvYieldEtaStatBothMeson_scaledx4_0010->Draw("p,same,e0");	
				
				DrawGammaSetMarkerTGraphAsym(graphInvYieldEtaSysBothMeson_2050, markerStyle2050, markerSizeComb, colorCombo2050 , colorCombo2050, 2, kTRUE);
				graphInvYieldEtaSysBothMeson_2050->Draw("E2same");
				
				DrawGammaSetMarkerTGraphAsym(graphInvYieldEtaStatBothMeson_2050, markerStyle2050, markerSizeComb, colorCombo2050, colorCombo2050);
				if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphInvYieldEtaStatBothMeson_2050);
				graphInvYieldEtaStatBothMeson_2050->Draw("p,same,e0");	
				
				
				TLegend* legtheory = new TLegend(0.6,0.85,0.8,0.95);
				legtheory->SetFillColor(0);
				legtheory->SetLineColor(0);
				legtheory->SetTextFont(42);
				legtheory->SetTextSize(FontSize);
				legtheory->SetHeader("PRC 90, 014906 (2014)");
				legtheory->AddEntry(TheoryCracowEtaLowPt_scaledx4_0010,"NEQ SHM 0-10%", "l");
				legtheory->AddEntry(TheoryCracowEtaLowPt_2050,"NEQ SHM 20-50%", "l");
				legtheory->Draw("same");
				
				TLegend* legtheory2 = new TLegend(0.6,0.73,0.8,0.83);
				legtheory2->SetFillColor(0);
				legtheory2->SetLineColor(0);
				legtheory2->SetTextFont(42);
				legtheory2->SetTextSize(FontSize);
				legtheory2->SetHeader("PRC 89, 064903 (2014)");
				legtheory2->AddEntry(histoEPOSEtaRebin_0010,"EPOS 3 0-10%", "l");
				legtheory2->AddEntry(histoEPOSEtaRebin_2050,"EPOS 2 20-50%", "l");
				legtheory2->Draw("same");
				
				TLegend* legendInvYieldSectionPi0LHC11h_onlyPbPb = new TLegend(0.2,0.13,0.7,0.23);
				legendInvYieldSectionPi0LHC11h_onlyPbPb->SetMargin(0.17);
				legendInvYieldSectionPi0LHC11h_onlyPbPb->SetFillColor(0);
				legendInvYieldSectionPi0LHC11h_onlyPbPb->SetLineColor(0);
				legendInvYieldSectionPi0LHC11h_onlyPbPb->SetTextFont(42);
				legendInvYieldSectionPi0LHC11h_onlyPbPb->SetTextSize(FontSize);
                legendInvYieldSectionPi0LHC11h_onlyPbPb->SetNColumns(2);
				legendInvYieldSectionPi0LHC11h_onlyPbPb->SetHeader(collisionSystem2760GeV.Data());
                legendInvYieldSectionPi0LHC11h_onlyPbPb->AddEntry(graphInvYieldEtaSysBothMeson_scaledx4_0010,Form("#eta,   %s",cent0010.Data()),"fp");
				legendInvYieldSectionPi0LHC11h_onlyPbPb->AddEntry(graphInvYieldPi0SysBothMeson_scaledx104_0010,Form("#pi^{0},   %s",cent0010.Data()),"fp");
                legendInvYieldSectionPi0LHC11h_onlyPbPb->AddEntry(graphInvYieldPi0SysBothMeson_2050,Form("#eta, %s",cent2050.Data()),"fp");
				legendInvYieldSectionPi0LHC11h_onlyPbPb->AddEntry(graphInvYieldPi0SysBothMeson_scaledx100_2050,Form("#pi^{0}, %s",cent2050.Data()),"fp");
				legendInvYieldSectionPi0LHC11h_onlyPbPb->Draw();
                
                labelFactorPi0104->Draw();
                labelFactorPi0100->Draw();
                labelFactorEta4->Draw();
			}
		}

// 		TGraphAsymmErrors *graphInvSectionCombStatPi02760GeVPlot_scaledx001 = (TGraphAsymmErrors*)graphInvSectionCombStatPi02760GeVPlot->Clone("");
// 		graphInvSectionCombStatPi02760GeVPlot_scaledx001 = ScaleGraph(graphInvSectionCombStatPi02760GeVPlot_scaledx001,0.1); 
// 		DrawGammaSetMarkerTGraphAsym(graphInvSectionCombStatPi02760GeVPlot_scaledx001, 30,markerSizepp, kBlack , kBlack);
// 		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphInvSectionCombStatPi02760GeVPlot_scaledx001);
// 		graphInvSectionCombStatPi02760GeVPlot_scaledx001->Draw("p,same,e0");
// 
// 		TGraphAsymmErrors *graphInvSectionCombSysPi02760GeVPlot_scaledx001 = (TGraphAsymmErrors*)graphInvSectionCombSysPi02760GeVPlot->Clone("");
// 		graphInvSectionCombSysPi02760GeVPlot_scaledx001 = ScaleGraph(graphInvSectionCombSysPi02760GeVPlot_scaledx001,0.1); 		
// 		DrawGammaSetMarkerTGraphAsym(graphInvSectionCombSysPi02760GeVPlot_scaledx001, 30,markerSizepp, kBlack , kBlack, widthLinesBoxes, kTRUE);
// 		graphInvSectionCombSysPi02760GeVPlot_scaledx001->Draw("E2same");
// 
// 
// 		TGraphAsymmErrors *graphInvSectionCombStatEta2760GeVPlot_scaledx001 = (TGraphAsymmErrors*)graphInvSectionCombStatEta2760GeVPlot->Clone("");
// 		graphInvSectionCombStatEta2760GeVPlot_scaledx001 = ScaleGraph(graphInvSectionCombStatEta2760GeVPlot_scaledx001,0.01); 
// 		DrawGammaSetMarkerTGraphAsym(graphInvSectionCombStatEta2760GeVPlot_scaledx001, markerStylepp,markerSizepp, kBlack , kBlack);
// 		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphInvSectionCombStatEta2760GeVPlot_scaledx001);
// 		graphInvSectionCombStatEta2760GeVPlot_scaledx001->Draw("p,same,e0");
// 
// 		TGraphAsymmErrors *graphInvSectionCombSysEta2760GeVPlot_scaledx001 = (TGraphAsymmErrors*)graphInvSectionCombSysEta2760GeVPlot->Clone("");
// 		graphInvSectionCombSysEta2760GeVPlot_scaledx001 = ScaleGraph(graphInvSectionCombSysEta2760GeVPlot_scaledx001,0.01); 		
// 		DrawGammaSetMarkerTGraphAsym(graphInvSectionCombSysEta2760GeVPlot_scaledx001, markerStylepp,markerSizepp, kBlack , kBlack, widthLinesBoxes, kTRUE);
// 		graphInvSectionCombSysEta2760GeVPlot_scaledx001->Draw("E2same");		
		
	canvasInvYieldSectionPi0LHC11h->SaveAs(Form("%s/YieldCombinedLHC11h_BothMesonswithPP.%s",outputDir.Data(),suffix.Data()));
	canvasInvYieldSectionPi0LHC11h->SaveAs(Form("%s/YieldCombinedLHC11h_BothMesonswithPP.%s",paperPlots.Data(),suffix.Data()));
	
	
		
	canvasInvYieldSectionPi0LHC11h->cd();
    histo2DInvYieldSectionLHC11hwithPP->Draw("copy");
// 	histo2DInvYieldSectionPi0LHC11h->GetXaxis()->SetRangeUser(minPtRange-0.5,maxPtRange);

		TGraphAsymmErrors *graphNLOforPi0PP = (TGraphAsymmErrors*)graphNLOCalcDSS14InvYieldPi0Band->Clone("graphNLOforPi0PP");
		TGraphAsymmErrors *graphNLOforPi00010 = (TGraphAsymmErrors*)graphNLOCalcDSS14InvYieldPi0Band->Clone("graphNLOforPi00010");		
		graphNLOforPi00010 = ScaleGraph(graphNLOforPi00010,nColl0010*4);
		TGraphAsymmErrors *graphNLOforPi02050 = (TGraphAsymmErrors*)graphNLOCalcDSS14InvYieldPi0Band->Clone("graphNLOforPi02050");
		graphNLOforPi02050 = ScaleGraph(graphNLOforPi02050,nColl2050);

		TGraphAsymmErrors *graphNLOforEtaPP = (TGraphAsymmErrors*)graphNLOCalcDSS07InvYieldEtaBand->Clone("graphNLOforEtaPP");
		TGraphAsymmErrors *graphNLOforEta0010 = (TGraphAsymmErrors*)graphNLOCalcDSS07InvYieldEtaBand->Clone("graphNLOforEta0010");		
		graphNLOforEta0010 = ScaleGraph(graphNLOforEta0010,nColl0010*4);
		TGraphAsymmErrors *graphNLOforEta2050 = (TGraphAsymmErrors*)graphNLOCalcDSS07InvYieldEtaBand->Clone("graphNLOforEta2050");
		graphNLOforEta2050 = ScaleGraph(graphNLOforEta2050,nColl2050);
		
		TLatex *labelNLO;
		labelNLO= new TLatex(0.61,0.75,"PDF: MSTW, FF: DSS14");
		SetStyleTLatex( labelNLO, 0.035,4);

		TLegend* legendSpectraPPnlo = new TLegend(0.2,0.12,0.52,0.34);//0.6,0.7,0.96,0.87);
		legendSpectraPPnlo->SetFillColor(0);
		legendSpectraPPnlo->SetLineColor(0);
		legendSpectraPPnlo->SetMargin(0.17);
		legendSpectraPPnlo->SetTextFont(42);
		legendSpectraPPnlo->SetTextSize(FontSize);

		if(meson.CompareTo("Pi0")==0){
						
			legendSpectraPPnlo->SetHeader(collisionSystemPP2760GeV.Data());
			legendSpectraPPnlo->AddEntry(graphInvSectionCombSysPi02760GeVPlot,"EPJC 74 (2014) 3108","pf");
 			legendSpectraPPnlo->AddEntry(graphNLOforPi0PP,"pQCD NLO #mu = #it{p}_{T}","l");
			legendSpectraPPnlo->AddEntry((TObject*)0,"PDF: MSTW, FF: DSS14","");
			legendSpectraPPnlo->AddEntry((TObject*)0,"PRD 91 014035",""); 
			legendSpectraPPnlo->AddEntry((TObject*)0,"scaled by <N_{coll}>","");
			legendSpectraPPnlo->Draw();
			
			DrawGammaSetMarkerTGraphAsym(graphNLOforPi0PP, 0, 0, colorNLO, colorNLO, 3, kTRUE, colorNLO);
			graphNLOforPi0PP->Draw("3,same");
			DrawGammaSetMarkerTGraphAsym(graphNLOforPi00010, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
			graphNLOforPi00010->Draw("3,same");
			DrawGammaSetMarkerTGraphAsym(graphNLOforPi02050, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
			graphNLOforPi02050->Draw("3,same");
			
			DrawGammaSetMarkerTGraphAsym(graphInvSectionCombStatPi02760GeVPlot, markerStylepp,markerSizepp, kBlack , kBlack);
			graphInvSectionCombStatPi02760GeVPlot->Draw("p,same,e0");
			DrawGammaSetMarkerTGraphAsym(graphInvSectionCombSysPi02760GeVPlot, markerStylepp,markerSizepp, kBlack , kBlack, widthLinesBoxes, kTRUE);//, colorComb1020-5);
			graphInvSectionCombSysPi02760GeVPlot->Draw("E2same");


		} else if(meson.CompareTo("Eta")==0){

			legendSpectraPPnlo->SetHeader(collisionSystemPP2760GeV.Data());
			legendSpectraPPnlo->AddEntry(graphInvSectionCombSysEta2760GeVPlot,"JPG 38 (2011) 124076","pf"); 
 			legendSpectraPPnlo->AddEntry(graphNLOforEtaPP,"pQCD NLO #mu = #it{p}_{T}","l");
			legendSpectraPPnlo->AddEntry((TObject*)0,"PDF: MSTW, FF: DSS07","");
			legendSpectraPPnlo->AddEntry((TObject*)0,"PRD 91 014035",""); 
			legendSpectraPPnlo->AddEntry((TObject*)0,"scaled by <N_{coll}>","");
			legendSpectraPPnlo->Draw();

			DrawGammaSetMarkerTGraphAsym(graphNLOforEtaPP, 0, 0, colorNLO, colorNLO, 3, kTRUE, colorNLO);
			graphNLOforEtaPP->Draw("3,same");
			DrawGammaSetMarkerTGraphAsym(graphNLOforEta0010, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
			graphNLOforEta0010->Draw("3,same");
			DrawGammaSetMarkerTGraphAsym(graphNLOforEta2050, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
			graphNLOforEta2050->Draw("3,same");

			DrawGammaSetMarkerTGraphAsym(graphInvSectionCombStatEta2760GeVPlot, markerStylepp,markerSizepp, kBlack , kBlack);
			graphInvSectionCombStatEta2760GeVPlot->Draw("p,same,e0");
			DrawGammaSetMarkerTGraphAsym(graphInvSectionCombSysEta2760GeVPlot, markerStylepp,markerSizepp, kBlack , kBlack, widthLinesBoxes, kTRUE);//, colorComb1020-5);
			graphInvSectionCombSysEta2760GeVPlot->Draw("E2same");
		}
		
		legendInvYieldSectionPi0LHC11h_WitPP->Draw();
		labelFactorUpper->Draw();
// 		labelPreliminary->Draw();
		labelSyst->Draw();
		
		graphCombInvYieldSys2760GeVALHC11hScaled_0010->Draw("E2same");
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombInvYieldStat2760GeVALHC11hScaled_0010);
		graphCombInvYieldStat2760GeVALHC11hScaled_0010->Draw("p,same,e0");	
	
		graphCombInvYieldSys2760GeVALHC11h_2050->Draw("E2same");
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombInvYieldStat2760GeVALHC11h_2050);
		graphCombInvYieldStat2760GeVALHC11h_2050->Draw("p,same,e0");	
	
	canvasInvYieldSectionPi0LHC11h->Update();
	canvasInvYieldSectionPi0LHC11h->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataOnly_withPPandNLO.%s",outputDir.Data(),meson.Data(),suffix.Data()));
	canvasInvYieldSectionPi0LHC11h->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataOnly_withPPandNLO.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
	
	
	canvasInvYieldSectionPi0LHC11h->cd();
// 	histo2DInvYieldSectionPi0LHC11h->GetXaxis()->SetRangeUser(minPtRange,maxPtRange);
	histo2DInvYieldSectionPi0LHC11h->DrawCopy();

		DrawGammaSetMarkerTGraphAsym(TheoryCracowLowPt_0010, 0, 0, colorCracow0010,colorCracow0010, 5, kTRUE, colorCracow0010);
		TheoryCracowLowPt_0010->Draw("l,same");
		DrawGammaSetMarkerTGraphAsym(TheoryCracowLowPt_2050, 0, 0, colorCracow2050,colorCracow2050, 5, kTRUE,colorCracow2050);
		TheoryCracowLowPt_2050->Draw("l,same");

		histoEPOSRebin_0010->Scale(4);
		DrawGammaSetMarker(histoEPOSRebin_0010,1,2,colorEPOS0010,colorEPOS0010);
		histoEPOSRebin_0010->SetLineWidth(2);
		histoEPOSRebin_0010->GetXaxis()->SetRangeUser(0.,13.);
		histoEPOSRebin_0010->Draw("same,c,histo");
		
		DrawGammaSetMarker(histoEPOSRebin_2050,1,2,colorEPOS2050,colorEPOS2050);
		histoEPOSRebin_2050->SetLineWidth(2);
		histoEPOSRebin_2050->GetXaxis()->SetRangeUser(0.,13.);
		histoEPOSRebin_2050->Draw("same,c,histo");
		
		graphCombInvYieldSys2760GeVALHC11hScaled_0010->Draw("E2same");
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombInvYieldStat2760GeVALHC11hScaled_0010);
		graphCombInvYieldStat2760GeVALHC11hScaled_0010->Draw("p,same,e0");	
	
		graphCombInvYieldSys2760GeVALHC11h_2050->Draw("E2same");
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombInvYieldStat2760GeVALHC11h_2050);
		graphCombInvYieldStat2760GeVALHC11h_2050->Draw("p,same,e0");	
		
		TLegend* legtheory = new TLegend(0.2,0.13,0.52,0.24);
		legtheory->SetFillColor(0);
		legtheory->SetLineColor(0);
		legtheory->SetTextFont(42);
		legtheory->SetTextSize(FontSize);
		legtheory->SetHeader("PRC 90, 014906 (2014)");
		legtheory->AddEntry(TheoryCracowLowPt_0010,"NEQ SHM 0-10%", "l");
		legtheory->AddEntry(TheoryCracowLowPt_2050,"NEQ SHM 20-50%", "l");
		legtheory->Draw("same");
		
		TLegend* legtheory2 = new TLegend(0.2,0.24,0.52,0.35);
		legtheory2->SetFillColor(0);
		legtheory2->SetLineColor(0);
		legtheory2->SetTextFont(42);
		legtheory2->SetTextSize(FontSize);
		legtheory2->SetHeader("PRC 89, 064903 (2014)");
		legtheory2->AddEntry(histoEPOSRebin_0010,"EPOS 3 0-10%", "l");
		legtheory2->AddEntry(histoEPOSRebin_2050,"EPOS 2 20-50%", "l");
		legtheory2->Draw("same");
		
		legendInvYieldSectionPi0LHC11h_WitPP->Draw();
		labelFactorLowerOnlyPbPb->Draw();
		labelSyst->Draw();
// 		labelPreliminary->Draw();

	canvasInvYieldSectionPi0LHC11h->Update();
	canvasInvYieldSectionPi0LHC11h->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModels.%s",outputDir.Data(),meson.Data(),suffix.Data()));
	canvasInvYieldSectionPi0LHC11h->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModels.%s",paperPlots.Data(),meson.Data(),suffix.Data()));

	
	canvasInvYieldSectionPi0LHC11h->cd();
// 	histo2DInvYieldSectionPi0LHC11h->GetXaxis()->SetRangeUser(minPtRange-0.5,maxPtRange);
    histo2DInvYieldSectionLHC11hwithPP->Draw("copy");
		

		TLegend* legendSpectraPP2 = new TLegend(0.595,0.74,0.91,0.82);
		legendSpectraPP2->SetFillColor(0);
		legendSpectraPP2->SetLineColor(0);
		legendSpectraPP2->SetTextFont(42);
		legendSpectraPP2->SetMargin(0.17);
		legendSpectraPP2->SetTextSize(FontSize);
		legendSpectraPP2->SetHeader(collisionSystemPP2760GeV.Data());

		if(meson.CompareTo("Pi0")==0){
			
			DrawGammaSetMarkerTGraphAsym(graphInvSectionCombStatPi02760GeVPlot, markerStylepp,markerSizepp, kBlack , kBlack);
			graphInvSectionCombStatPi02760GeVPlot->Draw("p,same,e0");
			DrawGammaSetMarkerTGraphAsym(graphInvSectionCombSysPi02760GeVPlot, markerStylepp,markerSizepp, kBlack , kBlack, widthLinesBoxes, kTRUE);
			graphInvSectionCombSysPi02760GeVPlot->Draw("E2same");
			
			legendSpectraPP2->AddEntry(graphInvSectionCombSysPi02760GeVPlot,"EPJC 74 (2014) 3108","pf");

		} else if(meson.CompareTo("Eta")==0){

			DrawGammaSetMarkerTGraphAsym(graphInvSectionCombStatEta2760GeVPlot, markerStylepp,markerSizepp, kBlack , kBlack);
			graphInvSectionCombStatEta2760GeVPlot->Draw("p,same,e0");
			DrawGammaSetMarkerTGraphAsym(graphInvSectionCombSysEta2760GeVPlot, markerStylepp,markerSizepp, kBlack , kBlack, widthLinesBoxes, kTRUE);
			graphInvSectionCombSysEta2760GeVPlot->Draw("E2same");
			
			legendSpectraPP2->AddEntry(graphInvSectionCombSysEta2760GeVPlot,"JPG 38 (2011) 124076","pf");

		}
		legendSpectraPP2->Draw();
		
		legendInvYieldSectionPi0LHC11h_WitPP->Draw();
// 		labelPreliminary->Draw();
		labelSyst->Draw();
		labelFactorLower->Draw();
		
		graphCombInvYieldSys2760GeVALHC11hScaled_0010->Draw("E2same");
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombInvYieldStat2760GeVALHC11hScaled_0010);
		graphCombInvYieldStat2760GeVALHC11hScaled_0010->Draw("p,same,e0");	
	
		graphCombInvYieldSys2760GeVALHC11h_2050->Draw("E2same");
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombInvYieldStat2760GeVALHC11h_2050);
		graphCombInvYieldStat2760GeVALHC11h_2050->Draw("p,same,e0");	
	
		
		TheoryCracowLowPt_0010->Draw("l,same");
		TheoryCracowLowPt_2050->Draw("l,same");

		histoEPOSRebin_0010->Draw("same,c,histo");
		histoEPOSRebin_2050->Draw("same,c,histo");

		legtheory->Draw("same");
		legtheory2->Draw("same");

	canvasInvYieldSectionPi0LHC11h->Update();
	canvasInvYieldSectionPi0LHC11h->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModels_withPP.%s",outputDir.Data(),meson.Data(),suffix.Data()));
	canvasInvYieldSectionPi0LHC11h->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModels_withPP.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
	
	
	
	//*************************************************************************************************************
	//***************************** Plot with data + models and ratios ***********************************************
	//*************************************************************************************************************

	Double_t arrayBoundariesX1_XSec[2];
	Double_t arrayBoundariesY1_XSec[6];
	Double_t relativeMarginsXXSec[3];
	Double_t relativeMarginsYXSec[3];
	textSizeLabelsPixel = 90;
	ReturnCorrectValuesForCanvasScaling(2500,4000, 1, 5,0.13, 0.025, 0.003,0.05,arrayBoundariesX1_XSec,arrayBoundariesY1_XSec,relativeMarginsXXSec,relativeMarginsYXSec);
	
	TCanvas* canvasInvYieldSectionRatios = new TCanvas("canvasInvYieldSectionRatios","",0,0,2500,4000);  // gives the page size
	DrawGammaCanvasSettings( canvasInvYieldSectionRatios,  0.13, 0.02, 0.03, 0.06);

		TPad* padInvSectionSpec = new TPad("padInvSectionSpec", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[3], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[0],-1, -1, -2);
		DrawGammaPadSettings( padInvSectionSpec, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[0], relativeMarginsYXSec[1]);
		padInvSectionSpec->Draw();
		Double_t marginXSec = relativeMarginsXXSec[0]*2500;
		Double_t textsizeLabelsXSecUp = 0;
		Double_t textsizeFacXSecUp = 0;
		if (padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) < padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1())){
			textsizeLabelsXSecUp = (Double_t)textSizeLabelsPixel/padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) ;
			textsizeFacXSecUp = (Double_t)1./padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) ;
		} else {
			textsizeLabelsXSecUp = (Double_t)textSizeLabelsPixel/padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1());
			textsizeFacXSecUp = (Double_t)1./padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1());
		}
			
		TPad* padInvSectionLowPtRatio = new TPad("padInvSectionLowPtRatio", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[4], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[3],-1, -1, -2);
		DrawGammaPadSettings( padInvSectionLowPtRatio, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[1]);
		padInvSectionLowPtRatio->Draw();
		Double_t textsizeLabelsXSecMiddle = 0;
		Double_t textsizeFacXSecMiddle = 0;
		if (padInvSectionLowPtRatio->XtoPixel(padInvSectionLowPtRatio->GetX2()) < padInvSectionLowPtRatio->YtoPixel(padInvSectionLowPtRatio->GetY1())){
			textsizeLabelsXSecMiddle = (Double_t)textSizeLabelsPixel/padInvSectionLowPtRatio->XtoPixel(padInvSectionLowPtRatio->GetX2()) ;
			textsizeFacXSecMiddle = (Double_t)1./padInvSectionLowPtRatio->XtoPixel(padInvSectionLowPtRatio->GetX2()) ;
		} else {
			textsizeLabelsXSecMiddle = (Double_t)textSizeLabelsPixel/padInvSectionLowPtRatio->YtoPixel(padInvSectionLowPtRatio->GetY1());
			textsizeFacXSecMiddle = (Double_t)1./padInvSectionLowPtRatio->YtoPixel(padInvSectionLowPtRatio->GetY1());
		}
	
		TPad* padInvSectionEPOSRatio = new TPad("padInvSectionEPOSRatio", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[5], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[4],-1, -1, -2);
		DrawGammaPadSettings( padInvSectionEPOSRatio, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[2]);
		padInvSectionEPOSRatio->Draw();
		Double_t textsizeLabelsXSecDown = 0;
		Double_t textsizeFacXSecDown = 0;
		if (padInvSectionEPOSRatio->XtoPixel(padInvSectionEPOSRatio->GetX2()) < padInvSectionEPOSRatio->YtoPixel(padInvSectionEPOSRatio->GetY1())){
			textsizeLabelsXSecDown = (Double_t)textSizeLabelsPixel/padInvSectionEPOSRatio->XtoPixel(padInvSectionEPOSRatio->GetX2()) ;
			textsizeFacXSecDown = (Double_t)1./padInvSectionEPOSRatio->XtoPixel(padInvSectionEPOSRatio->GetX2()) ;
		} else {
			textsizeLabelsXSecDown = (Double_t)textSizeLabelsPixel/padInvSectionEPOSRatio->YtoPixel(padInvSectionEPOSRatio->GetY1());
			textsizeFacXSecDown = (Double_t)1./padInvSectionEPOSRatio->YtoPixel(padInvSectionEPOSRatio->GetY1());
		}
	
		
		padInvSectionSpec->cd();
		padInvSectionSpec->SetLogy(1);
		padInvSectionSpec->SetLogx(1);
		histo2DInvYieldSectionPi0LHC11h->GetXaxis()->SetRangeUser(minPtRange,maxPtRange);
		histo2DInvYieldSectionPi0LHC11h->Draw("copy");

		DrawGammaSetMarkerTGraphAsym(graphCombInvYieldSys2760GeVALHC11h_0010, markerStyleComb, markerSizeComb, colorCombo0010 , colorCombo0010, 2, kTRUE);
		graphCombInvYieldSys2760GeVALHC11h_0010->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphCombInvYieldStat2760GeVALHC11h_0010, markerStyleComb, markerSizeComb, colorCombo0010, colorCombo0010);
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombInvYieldStat2760GeVALHC11h_0010);
		graphCombInvYieldStat2760GeVALHC11h_0010->Draw("p,same,e0");	
	
		TGraphAsymmErrors *TheoryCracowLowPtforRatio_0010 = (TGraphAsymmErrors*)TheoryCracowLowPt_0010->Clone("TheoryCracowLowPtforRatio_0010");
		TheoryCracowLowPtforRatio_0010 = ScaleGraph(TheoryCracowLowPtforRatio_0010,1./4);
		DrawGammaSetMarkerTGraphAsym(TheoryCracowLowPtforRatio_0010, 0, 0, colorCracowRatio, colorCracowRatio, 5, kTRUE, colorCracowRatio);
		TheoryCracowLowPtforRatio_0010->Draw("l,same");

		histoEPOSRebin_0010->Scale(1./4);
		histoEPOSRebin_0010->SetLineColor(colorEPOSRatio);
		histoEPOSRebin_0010->Draw("same,c,histo");

        fitTCMInvYield2760GeVLHC11h_0010->SetParameter(0,fitTCMInvYield2760GeVLHC11h_0010->GetParameter(0)/4);
        fitTCMInvYield2760GeVLHC11h_0010->SetParameter(2,fitTCMInvYield2760GeVLHC11h_0010->GetParameter(2)/4);
		fitTCMInvYield2760GeVLHC11h_0010->SetLineStyle(4);
		fitTCMInvYield2760GeVLHC11h_0010->Draw("same");		
		
		TLatex *labelDetSysRatioToModelsLHC11h;
		if(meson.CompareTo("Pi0")==0){
			labelDetSysRatioToModelsLHC11h= new TLatex(0.5,0.93,"#pi^{0} #rightarrow #gamma#gamma");
		} else if(meson.CompareTo("Eta")==0){
			labelDetSysRatioToModelsLHC11h= new TLatex(0.5,0.93,"#eta #rightarrow #gamma#gamma");
		}
		SetStyleTLatex( labelDetSysRatioToModelsLHC11h, 0.035,4);
		labelDetSysRatioToModelsLHC11h->Draw();
		
		TLatex *labelEnergyRatioToModelsLHC11h = new TLatex(0.5,0.88,collisionSystemPbPb0010.Data());
		SetStyleTLatex( labelEnergyRatioToModelsLHC11h, 0.035,4);
		labelEnergyRatioToModelsLHC11h->Draw();
// 		labelPreliminary->Draw();
		
		TLegend* legendXsectionPaper = GetAndSetLegend2(0.17, 0.3-4*0.06, 0.5, 0.3, 0.85* textSizeLabelsPixel);
		legendXsectionPaper->SetNColumns(1);
		legendXsectionPaper->SetMargin(0.2);
		if(meson.CompareTo("Pi0")==0) legendXsectionPaper->AddEntry(graphCombInvYieldSys2760GeVALHC11h_0010,"#pi^{0} ALICE","pf");
		else if(meson.CompareTo("Eta")==0) legendXsectionPaper->AddEntry(graphCombInvYieldSys2760GeVALHC11h_0010,"#eta ALICE","pf");
		legendXsectionPaper->AddEntry(TheoryCracowLowPtforRatio_0010,"NEQ SHM","l");
		legendXsectionPaper->AddEntry(histoEPOSRebin_0010,"EPOS 3","l");
		legendXsectionPaper->AddEntry(fitTCMInvYield2760GeVLHC11h_0010,"#it{A}_{e} exp(-#it{E}_{T, kin}/#it{T}_{e}) + #it{A}/#(){1 + #frac{#it{p}_{T}^{2}}{#it{T}^{2}#upoint n}}^{-n}","l");
		legendXsectionPaper->Draw();		
		
	padInvSectionLowPtRatio->cd();
	padInvSectionLowPtRatio->SetLogx(1);
		TH2F * ratio2DLowPt = new TH2F("ratio2DLowPt","ratio2DLowPt",1000,0.8,30.,1000,0.2,2.1);
		SetStyleHistoTH2ForGraphs(ratio2DLowPt, "#it{p}_{T} (GeV/#it{c})","#frac{NEQ SHM, Data}{fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle, 
								  0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.2/(textsizeFacXSecMiddle*marginXSec), 510, 505);
		ratio2DLowPt->GetXaxis()->SetRangeUser(minPtRange,maxPtRange);
		ratio2DLowPt->GetYaxis()->SetNdivisions(505);
		ratio2DLowPt->GetXaxis()->SetMoreLogLabels(kTRUE);
		ratio2DLowPt->GetXaxis()->SetNoExponent(kTRUE);
		ratio2DLowPt->GetYaxis()->SetLabelFont(42);
		ratio2DLowPt->GetYaxis()->SetLabelOffset(+0.01);
// 		ratio2DLowPt->GetXaxis()->SetTickLength(0.07);
		ratio2DLowPt->DrawCopy();

		TGraphAsymmErrors* graphRatioCombLowPt2760GeVLHC11h_0010 = (TGraphAsymmErrors*)TheoryCracowLowPtforRatio_0010->Clone("graphRatioCombLowPt2760GeVLHC11h_0010");	
		graphRatioCombLowPt2760GeVLHC11h_0010 = CalculateGraphErrRatioToFit(graphRatioCombLowPt2760GeVLHC11h_0010, fitTCMInvYield2760GeVLHC11h_0010); 
		DrawGammaSetMarkerTGraphAsym(graphRatioCombLowPt2760GeVLHC11h_0010, 0, 0, colorCracowRatio, colorCracowRatio, 5, kTRUE, colorCracowRatio);
		graphRatioCombLowPt2760GeVLHC11h_0010->Draw("l,same");

		DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat2760GeVLHC11h_0010, markerStyleComb, markerSizeComb*2, kBlack, kBlack, widthLinesBoxes, kFALSE);
		graphRatioCombCombFitStat2760GeVLHC11h_0010->SetLineWidth(widthLinesBoxes);
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRatioCombCombFitStat2760GeVLHC11h_0010);
		graphRatioCombCombFitStat2760GeVLHC11h_0010->Draw("p,same,e0");

		DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys2760GeVLHC11h_0010, markerStyleComb, markerSizeComb*2, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
		graphRatioCombCombFitSys2760GeVLHC11h_0010->SetLineWidth(0);
		graphRatioCombCombFitSys2760GeVLHC11h_0010->Draw("2,same");
		
		DrawGammaLines(0.8, 30.,1., 1.,0.1,kGray);
		
	padInvSectionEPOSRatio->cd();
	padInvSectionEPOSRatio->SetLogx(1);

		TH2F * ratio2DPythia =  new TH2F("ratio2DPythia","ratio2DPythia",1000,0.8,30.,1000,0.2,2.1);
		SetStyleHistoTH2ForGraphs(ratio2DPythia, "#it{p}_{T} (GeV/#it{c})","#frac{EPOS, Data}{fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown, 
								  0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 0.9,0.2/(textsizeFacXSecDown*marginXSec), 510, 505);
		ratio2DPythia->GetXaxis()->SetRangeUser(minPtRange,maxPtRange);
		ratio2DPythia->GetYaxis()->SetNdivisions(505);
		ratio2DPythia->GetXaxis()->SetMoreLogLabels(kTRUE);
		ratio2DPythia->GetXaxis()->SetNoExponent(kTRUE);
		ratio2DPythia->GetYaxis()->SetLabelFont(42);
		ratio2DPythia->GetYaxis()->SetLabelOffset(+0.01);
// 		ratio2DPythia->GetXaxis()->SetTickLength(0.07);
		ratio2DPythia->DrawCopy();

		TH1D* histoRatioEPOSToFit2760GeVLHC11h_0010 = (TH1D*) histoEPOSRebin_0010->Clone();     
		histoRatioEPOSToFit2760GeVLHC11h_0010 = CalculateHistoRatioToFit (histoRatioEPOSToFit2760GeVLHC11h_0010, fitTCMInvYield2760GeVLHC11h_0010); 
		DrawGammaSetMarker(histoRatioEPOSToFit2760GeVLHC11h_0010, 24, 1.5,colorEPOSRatio,colorEPOSRatio);  
		histoRatioEPOSToFit2760GeVLHC11h_0010->SetLineWidth(widthCommonFit);
		histoRatioEPOSToFit2760GeVLHC11h_0010->Draw("same,hist,c");
		
		graphRatioCombCombFitStat2760GeVLHC11h_0010->Draw("p,same,e0");
		graphRatioCombCombFitSys2760GeVLHC11h_0010->Draw("2,same");

		DrawGammaLines(0.8, 30.,1., 1.,0.1,kGray);
		
	canvasInvYieldSectionRatios->Update();
// 	canvasInvYieldSectionRatios->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModelsAndRatios_0010.%s",outputDir.Data(),meson.Data(),suffix.Data()));
// 	canvasInvYieldSectionRatios->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModelsAndRatios_0010.%s",paperPlots.Data(),meson.Data(),suffix.Data()));

	
	TCanvas* canvasInvYieldSectionRatios_2050 = new TCanvas("canvasInvYieldSectionRatios_2050","",0,0,2500,4000);  // gives the page size
	DrawGammaCanvasSettings( canvasInvYieldSectionRatios_2050,  0.13, 0.02, 0.03, 0.06);
	
		TPad* padInvSectionSpec_2050 = new TPad("padInvSectionSpec_2050", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[3], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[0],-1, -1, -2);
		DrawGammaPadSettings( padInvSectionSpec_2050, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[0], relativeMarginsYXSec[1]);
		padInvSectionSpec_2050->Draw();
		if (padInvSectionSpec_2050->XtoPixel(padInvSectionSpec_2050->GetX2()) < padInvSectionSpec_2050->YtoPixel(padInvSectionSpec_2050->GetY1())){
			textsizeLabelsXSecUp = (Double_t)textSizeLabelsPixel/padInvSectionSpec_2050->XtoPixel(padInvSectionSpec_2050->GetX2()) ;
			textsizeFacXSecUp = (Double_t)1./padInvSectionSpec_2050->XtoPixel(padInvSectionSpec_2050->GetX2()) ;
		} else {
			textsizeLabelsXSecUp = (Double_t)textSizeLabelsPixel/padInvSectionSpec_2050->YtoPixel(padInvSectionSpec_2050->GetY1());
			textsizeFacXSecUp = (Double_t)1./padInvSectionSpec_2050->YtoPixel(padInvSectionSpec_2050->GetY1());
		}
			
		TPad* padInvSectionLowPtRatio_2050 = new TPad("padInvSectionLowPtRatio_2050", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[4], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[3],-1, -1, -2);
		DrawGammaPadSettings( padInvSectionLowPtRatio_2050, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[1]);
		padInvSectionLowPtRatio_2050->Draw();
		if (padInvSectionLowPtRatio_2050->XtoPixel(padInvSectionLowPtRatio_2050->GetX2()) < padInvSectionLowPtRatio_2050->YtoPixel(padInvSectionLowPtRatio_2050->GetY1())){
			textsizeLabelsXSecMiddle = (Double_t)textSizeLabelsPixel/padInvSectionLowPtRatio_2050->XtoPixel(padInvSectionLowPtRatio_2050->GetX2()) ;
			textsizeFacXSecMiddle = (Double_t)1./padInvSectionLowPtRatio_2050->XtoPixel(padInvSectionLowPtRatio_2050->GetX2()) ;
		} else {
			textsizeLabelsXSecMiddle = (Double_t)textSizeLabelsPixel/padInvSectionLowPtRatio_2050->YtoPixel(padInvSectionLowPtRatio_2050->GetY1());
			textsizeFacXSecMiddle = (Double_t)1./padInvSectionLowPtRatio_2050->YtoPixel(padInvSectionLowPtRatio_2050->GetY1());
		}
	
		TPad* padInvSectionEPOSRatio_2050 = new TPad("padInvSectionEPOSRatio_2050", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[5], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[4],-1, -1, -2);
		DrawGammaPadSettings( padInvSectionEPOSRatio_2050, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[2]);
		padInvSectionEPOSRatio_2050->Draw();
		if (padInvSectionEPOSRatio_2050->XtoPixel(padInvSectionEPOSRatio_2050->GetX2()) < padInvSectionEPOSRatio_2050->YtoPixel(padInvSectionEPOSRatio_2050->GetY1())){
			textsizeLabelsXSecDown = (Double_t)textSizeLabelsPixel/padInvSectionEPOSRatio_2050->XtoPixel(padInvSectionEPOSRatio_2050->GetX2()) ;
			textsizeFacXSecDown = (Double_t)1./padInvSectionEPOSRatio_2050->XtoPixel(padInvSectionEPOSRatio_2050->GetX2()) ;
		} else {
			textsizeLabelsXSecDown = (Double_t)textSizeLabelsPixel/padInvSectionEPOSRatio_2050->YtoPixel(padInvSectionEPOSRatio_2050->GetY1());
			textsizeFacXSecDown = (Double_t)1./padInvSectionEPOSRatio_2050->YtoPixel(padInvSectionEPOSRatio_2050->GetY1());
		}

		padInvSectionSpec_2050->cd();
		padInvSectionSpec_2050->SetLogy(1);
		padInvSectionSpec_2050->SetLogx(1);
		histo2DInvYieldSectionPi0LHC11h->DrawCopy();
	
		TGraphAsymmErrors *TheoryCracowLowPtforRatio_2050 = (TGraphAsymmErrors*)TheoryCracowLowPt_2050->Clone("TheoryCracowLowPtforRatio_2050");
		DrawGammaSetMarkerTGraphAsym(TheoryCracowLowPtforRatio_2050, 0, 0, colorCracowRatio,colorCracowRatio, 5, kTRUE,colorCracowRatio);
		TheoryCracowLowPtforRatio_2050->Draw("l,same");

		histoEPOSRebin_2050->SetLineColor(colorEPOSRatio);
		histoEPOSRebin_2050->Draw("same,c,histo");

		graphCombInvYieldSys2760GeVALHC11h_2050->Draw("E2same");
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombInvYieldStat2760GeVALHC11h_2050);
		graphCombInvYieldStat2760GeVALHC11h_2050->Draw("p,same,e0");	

		fitTCMInvYield2760GeVLHC11h_2050->SetLineStyle(4);
		fitTCMInvYield2760GeVLHC11h_2050->Draw("same");		

		TLatex *labelEnergyRatioToModelsLHC11h2 = new TLatex(0.5,0.88,collisionSystemPbPb2050.Data());
		SetStyleTLatex( labelEnergyRatioToModelsLHC11h2, 0.035,4);
		labelEnergyRatioToModelsLHC11h2->Draw();
		labelDetSysRatioToModelsLHC11h->Draw();
		
		TLegend* legendXsectionPaper2 = GetAndSetLegend2(0.17, 0.3-4*0.06, 0.5, 0.3, 0.85* textSizeLabelsPixel);
		legendXsectionPaper2->SetNColumns(1);
		legendXsectionPaper2->SetMargin(0.2);
		if(meson.CompareTo("Pi0")==0) legendXsectionPaper2->AddEntry(graphCombInvYieldSys2760GeVALHC11h_2050,"#pi^{0} ALICE","pf");
		else if(meson.CompareTo("Eta")==0) legendXsectionPaper2->AddEntry(graphCombInvYieldSys2760GeVALHC11h_2050,"#eta ALICE","pf");
		legendXsectionPaper2->AddEntry(TheoryCracowLowPtforRatio_2050,"NEQ SHM","l");
		legendXsectionPaper2->AddEntry(histoEPOSRebin_2050,"EPOS 2","l");
		legendXsectionPaper2->AddEntry(fitTCMInvYield2760GeVLHC11h_2050,"#it{A}_{e} exp(-#it{E}_{T, kin}/#it{T}_{e}) + #it{A}/#(){1 + #frac{#it{p}_{T}^{2}}{#it{T}^{2}#upoint n}}^{-n}","l");
		legendXsectionPaper2->Draw();		
		legendXsectionPaper2->Draw();		

		
		padInvSectionLowPtRatio_2050->cd();
		padInvSectionLowPtRatio_2050->SetLogx(1);
		ratio2DLowPt->DrawCopy();

		TGraphAsymmErrors* graphRatioCombLowPt2760GeVLHC11h_2050 = (TGraphAsymmErrors*)TheoryCracowLowPt_2050->Clone("graphRatioCombLowPt2760GeVLHC11h_2050");	
		graphRatioCombLowPt2760GeVLHC11h_2050 = CalculateGraphErrRatioToFit(graphRatioCombLowPt2760GeVLHC11h_2050, fitTCMInvYield2760GeVLHC11h_2050); 

// 		DrawGammaSetMarkerTGraphAsym(graphRatioCombLowPt2760GeVLHC11h_2050, 0, 0, colorCracowRatio,colorCracowRatio, widthLinesBoxes, kTRUE,colorCracowRatio);
		DrawGammaSetMarkerTGraphAsym(graphRatioCombLowPt2760GeVLHC11h_2050, 0, 0, colorCracowRatio,colorCracowRatio, 5, kTRUE,colorCracowRatio);
		graphRatioCombLowPt2760GeVLHC11h_2050->Draw("l,same");

		DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat2760GeVLHC11h_2050, markerStyleComb, markerSizeComb*2, kBlack, kBlack, widthLinesBoxes, kFALSE);
		graphRatioCombCombFitStat2760GeVLHC11h_2050->SetLineWidth(widthLinesBoxes);
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRatioCombCombFitStat2760GeVLHC11h_2050);
		graphRatioCombCombFitStat2760GeVLHC11h_2050->Draw("p,same,e0");

		DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys2760GeVLHC11h_2050, markerStyleComb, markerSizeComb*2, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
		graphRatioCombCombFitSys2760GeVLHC11h_2050->SetLineWidth(0);
		graphRatioCombCombFitSys2760GeVLHC11h_2050->Draw("2,same");
		
		DrawGammaLines(0.8, 30.,1., 1.,0.1,kGray);
		
		
		padInvSectionEPOSRatio_2050->cd();
		padInvSectionEPOSRatio_2050->SetLogx(1);
		ratio2DPythia->DrawCopy();

		TH1D* histoRatioEPOSToFit2760GeVLHC11h_2050 = (TH1D*) histoEPOSRebin_2050->Clone();     
		histoRatioEPOSToFit2760GeVLHC11h_2050 = CalculateHistoRatioToFit (histoRatioEPOSToFit2760GeVLHC11h_2050, fitTCMInvYield2760GeVLHC11h_2050); 

		DrawGammaSetMarker(histoRatioEPOSToFit2760GeVLHC11h_2050, 24, 1.5, colorEPOSRatio,colorEPOSRatio);  
		histoRatioEPOSToFit2760GeVLHC11h_2050->SetLineWidth(widthCommonFit);
// 		histoRatioEPOSToFit2760GeVLHC11h_2050->GetXaxis()->SetRangeUser(0.5,14);
		histoRatioEPOSToFit2760GeVLHC11h_2050->Draw("same,hist,c");
		
		graphRatioCombCombFitSys2760GeVLHC11h_2050->Draw("2,same");
		graphRatioCombCombFitStat2760GeVLHC11h_2050->Draw("p,same,e0");

		DrawGammaLines(0.8, 30.,1., 1.,0.1,kGray);
		
	canvasInvYieldSectionRatios_2050->Update();
// 	canvasInvYieldSectionRatios_2050->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModelsAndRatios_2050.%s",outputDir.Data(),meson.Data(),suffix.Data()));
// 	canvasInvYieldSectionRatios_2050->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModelsAndRatios_2050.%s",paperPlots.Data(),meson.Data(),suffix.Data()));	

    
    
    
    Double_t arrayBoundariesY1_XSec_ratios[3];
    Double_t relativeMarginsXXSec_ratios[3];
    Double_t relativeMarginsYXSec_ratios[3];
    textSizeLabelsPixel = 90;
    ReturnCorrectValuesForCanvasScaling(2500,2000, 1, 2,0.13, 0.025, 0.003,0.1,arrayBoundariesX1_XSec,arrayBoundariesY1_XSec_ratios,relativeMarginsXXSec_ratios,relativeMarginsYXSec_ratios);
    
    TCanvas* canvasInvYieldSectionOnlyRatios = new TCanvas("canvasInvYieldSectionOnlyRatios","",0,0,2500,2000);  //0,0,2500,4000);  // gives the page size
    DrawGammaCanvasSettings( canvasInvYieldSectionOnlyRatios,  0.13, 0.02, 0.03, 0.06);
            
        TPad* padInvSectionUpperRatio = new TPad("padInvSectionUpperRatio", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec_ratios[1], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec_ratios[0],-1, -1, -2);
        DrawGammaPadSettings( padInvSectionUpperRatio, relativeMarginsXXSec_ratios[0], relativeMarginsXXSec_ratios[2], relativeMarginsYXSec_ratios[1], relativeMarginsYXSec_ratios[1]);
        padInvSectionUpperRatio->Draw();
        if (padInvSectionUpperRatio->XtoPixel(padInvSectionUpperRatio->GetX2()) < padInvSectionUpperRatio->YtoPixel(padInvSectionUpperRatio->GetY1())){
            textsizeLabelsXSecMiddle = (Double_t)textSizeLabelsPixel/padInvSectionUpperRatio->XtoPixel(padInvSectionUpperRatio->GetX2()) ;
            textsizeFacXSecMiddle = (Double_t)1./padInvSectionUpperRatio->XtoPixel(padInvSectionUpperRatio->GetX2()) ;
        } else {
            textsizeLabelsXSecMiddle = (Double_t)textSizeLabelsPixel/padInvSectionUpperRatio->YtoPixel(padInvSectionUpperRatio->GetY1());
            textsizeFacXSecMiddle = (Double_t)1./padInvSectionUpperRatio->YtoPixel(padInvSectionUpperRatio->GetY1());
        }
    
        TPad* padInvSectionLowerRatio = new TPad("padInvSectionLowerRatio", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec_ratios[2], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec_ratios[1],-1, -1, -2);
        DrawGammaPadSettings( padInvSectionLowerRatio, relativeMarginsXXSec_ratios[0], relativeMarginsXXSec_ratios[2], relativeMarginsYXSec_ratios[1], relativeMarginsYXSec_ratios[2]);
        padInvSectionLowerRatio->Draw();
        if (padInvSectionLowerRatio->XtoPixel(padInvSectionLowerRatio->GetX2()) < padInvSectionLowerRatio->YtoPixel(padInvSectionLowerRatio->GetY1())){
            textsizeLabelsXSecDown = (Double_t)textSizeLabelsPixel/padInvSectionLowerRatio->XtoPixel(padInvSectionLowerRatio->GetX2()) ;
            textsizeFacXSecDown = (Double_t)1./padInvSectionLowerRatio->XtoPixel(padInvSectionLowerRatio->GetX2()) ;
        } else {
            textsizeLabelsXSecDown = (Double_t)textSizeLabelsPixel/padInvSectionLowerRatio->YtoPixel(padInvSectionLowerRatio->GetY1());
            textsizeFacXSecDown = (Double_t)1./padInvSectionLowerRatio->YtoPixel(padInvSectionLowerRatio->GetY1());
        }


    padInvSectionUpperRatio->cd();
    padInvSectionUpperRatio->SetLogx(1);
        TH2F * ratio2DUpperPadOnlyRatio = new TH2F("ratio2DUpperPadOnlyRatio","ratio2DUpperPadOnlyRatio",1000,0.8,30.,1000,0.01,2.1);
        SetStyleHistoTH2ForGraphs(ratio2DUpperPadOnlyRatio, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle, 
                                  0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.2/(textsizeFacXSecMiddle*marginXSec), 510, 505);
        ratio2DUpperPadOnlyRatio->GetXaxis()->SetRangeUser(minPtRange,maxPtRange);
        ratio2DUpperPadOnlyRatio->GetYaxis()->SetNdivisions(505);
        ratio2DUpperPadOnlyRatio->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DUpperPadOnlyRatio->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DUpperPadOnlyRatio->GetYaxis()->SetLabelFont(42);
        ratio2DUpperPadOnlyRatio->GetYaxis()->SetLabelOffset(+0.01);
//      ratio2DUpperPadOnlyRatio->GetXaxis()->SetTickLength(0.07);
        ratio2DUpperPadOnlyRatio->DrawCopy();

        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRatioCombCombFitStat2760GeVLHC11h_0010);
        graphRatioCombCombFitStat2760GeVLHC11h_0010->Draw("p,same,e0");
        graphRatioCombCombFitSys2760GeVLHC11h_0010->Draw("2,same");
        
        graphRatioCombLowPt2760GeVLHC11h_0010->Draw("l,same");
        histoRatioEPOSToFit2760GeVLHC11h_0010->Draw("same,hist,c");
        
       TLatex *labelcent0010 = new TLatex(0.8,0.88," 0-10%");
       SetStyleTLatex( labelcent0010, textsizeLabelsXSecMiddle,4);
       labelcent0010->Draw();
       
        TLegend* legendXsectionPaperOnlyRatios;
        if(meson.CompareTo("Pi0")==0)legendXsectionPaperOnlyRatios = GetAndSetLegend2(0.17, 0.45-4*0.065, 0.5, 0.45, 0.85* textSizeLabelsPixel);
//         else if(meson.CompareTo("Eta")==0)legendXsectionPaperOnlyRatios = GetAndSetLegend2(0.17, 0.95-4*0.07, 0.5, 0.95, 0.85* textSizeLabelsPixel);
        else if(meson.CompareTo("Eta")==0)legendXsectionPaperOnlyRatios = GetAndSetLegend2(0.17, 0.95-4*0.065, 0.5, 0.95, 0.85* textSizeLabelsPixel);
        legendXsectionPaperOnlyRatios->SetNColumns(1);
        legendXsectionPaperOnlyRatios->SetMargin(0.2);
        legendXsectionPaperOnlyRatios->SetHeader(collisionSystem2760GeV.Data());
        if(meson.CompareTo("Pi0")==0) legendXsectionPaperOnlyRatios->AddEntry(graphRatioCombCombFitSys2760GeVLHC11h_0010,"#pi^{0} ALICE","pf");
        else if(meson.CompareTo("Eta")==0) legendXsectionPaperOnlyRatios->AddEntry(graphRatioCombCombFitSys2760GeVLHC11h_0010,"#eta ALICE","pf");
        legendXsectionPaperOnlyRatios->AddEntry(graphRatioCombLowPt2760GeVLHC11h_0010,"NEQ SHM","l");
        legendXsectionPaperOnlyRatios->AddEntry(histoRatioEPOSToFit2760GeVLHC11h_0010,"EPOS","l");
//         if(meson.CompareTo("Eta")==0) legendXsectionPaperOnlyRatios->Draw();        


        DrawGammaLines(0.8, 30.,1., 1.,0.1,kGray);
        
    padInvSectionLowerRatio->cd();
    padInvSectionLowerRatio->SetLogx(1);


        TH2F * ratio2DLowerPadOnlyRatio =  new TH2F("ratio2DLowerPadOnlyRatio","ratio2DLowerPadOnlyRatio",1000,0.8,30.,1000,0.01,2.1);
        SetStyleHistoTH2ForGraphs(ratio2DLowerPadOnlyRatio, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown, 
                                  0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 0.9,0.2/(textsizeFacXSecDown*marginXSec), 510, 505);
        ratio2DLowerPadOnlyRatio->GetXaxis()->SetRangeUser(minPtRange,maxPtRange);
        ratio2DLowerPadOnlyRatio->GetYaxis()->SetNdivisions(505);
        ratio2DLowerPadOnlyRatio->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DLowerPadOnlyRatio->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DLowerPadOnlyRatio->GetYaxis()->SetLabelFont(42);
        ratio2DLowerPadOnlyRatio->GetYaxis()->SetLabelOffset(+0.01);
//      ratio2DLowerPadOnlyRatio->GetXaxis()->SetTickLength(0.07);
        ratio2DLowerPadOnlyRatio->DrawCopy();

// 		DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat2760GeVLHC11h_2050, markerStyle2050, markerSizeComb*2, kBlack, kBlack, widthLinesBoxes, kFALSE);
// 		DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys2760GeVLHC11h_2050, markerStyle2050, markerSizeComb*2, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);

        if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRatioCombCombFitStat2760GeVLHC11h_2050);
        graphRatioCombCombFitStat2760GeVLHC11h_2050->Draw("p,same,e0");
        graphRatioCombCombFitSys2760GeVLHC11h_2050->Draw("2,same");
        
        graphRatioCombLowPt2760GeVLHC11h_2050->Draw("l,same");
        histoRatioEPOSToFit2760GeVLHC11h_2050->Draw("same,hist,c");
        
        TLatex *labelcent2050 = new TLatex(0.8,0.88,"20-50%");
        SetStyleTLatex( labelcent2050, textsizeLabelsXSecDown,4);
        labelcent2050->Draw();

        if(meson.CompareTo("Eta")==0) legendXsectionPaperOnlyRatios->Draw();        
        if(meson.CompareTo("Pi0")==0)legendXsectionPaperOnlyRatios->Draw();        

        DrawGammaLines(0.8, 30.,1., 1.,0.1,kGray);

    canvasInvYieldSectionOnlyRatios->Update();
    canvasInvYieldSectionOnlyRatios->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModelsOnlyRatios.%s",outputDir.Data(),meson.Data(),suffix.Data()));
    canvasInvYieldSectionOnlyRatios->SaveAs(Form("%s/%s_YieldCombinedLHC11h_DataWithModelsOnlyRatios.%s",paperPlots.Data(),meson.Data(),suffix.Data()));   

	
	Double_t arrayX14padratios[3];
    Double_t arrayY14padratios[3];
    Double_t relX4padratios[3];
    Double_t relY4padratios[3];
	Int_t textSizeLabelsPixel4pad = 50;
    ReturnCorrectValuesForCanvasScaling(2300,1300, 2, 2,0.08, 0.025, 0.003,0.1,arrayX14padratios,arrayY14padratios,relX4padratios,relY4padratios);
    
    TCanvas* canvasInvYieldSectionOnlyRatios4pad = new TCanvas("canvasInvYieldSectionOnlyRatios4pad","",0,0,2300,1300);  //0,0,2500,4000);  // gives the page size
    DrawGammaCanvasSettings( canvasInvYieldSectionOnlyRatios4pad,  0., 0., 0., 0.);
            
      textsizeLabelsXSecMiddle = 0;
      textsizeFacXSecMiddle = 0;
      textsizeLabelsXSecDown = 0;
      textsizeFacXSecDown = 0;

      TPad* padInvYieldUpperRatio_pad1 = new TPad("padInvYieldUpperRatio_pad1", "", arrayX14padratios[0], arrayY14padratios[1], arrayX14padratios[1], arrayY14padratios[0],-1, -1, -2);
        DrawGammaPadSettings( padInvYieldUpperRatio_pad1, relX4padratios[0], relX4padratios[1], relY4padratios[0], relY4padratios[1]);
        padInvYieldUpperRatio_pad1->Draw();
        
        TPad* padInvYieldUpperRatio_pad2 = new TPad("padInvYieldUpperRatio_pad2", "", arrayX14padratios[1], arrayY14padratios[1], arrayX14padratios[2], arrayY14padratios[0],-1, -1, -2);
        DrawGammaPadSettings( padInvYieldUpperRatio_pad2, relX4padratios[1], relX4padratios[2], relY4padratios[0], relY4padratios[1]);
        padInvYieldUpperRatio_pad2->Draw();
        
        TPad* padInvYieldLowerRatio_pad1 = new TPad("padInvYieldLowerRatio_pad1", "", arrayX14padratios[0], arrayY14padratios[2], arrayX14padratios[1], arrayY14padratios[1],-1, -1, -2);
        DrawGammaPadSettings( padInvYieldLowerRatio_pad1, relX4padratios[0], relX4padratios[1], relY4padratios[1], relY4padratios[2]);
        padInvYieldLowerRatio_pad1->Draw();

        TPad* padInvYieldLowerRatio_pad2 = new TPad("padInvYieldLowerRatio_pad2", "", arrayX14padratios[1], arrayY14padratios[2], arrayX14padratios[2], arrayY14padratios[1],-1, -1, -2);
        DrawGammaPadSettings( padInvYieldLowerRatio_pad2, relX4padratios[1], relX4padratios[2], relY4padratios[1], relY4padratios[2]);
        padInvYieldLowerRatio_pad2->Draw();

        if (padInvYieldLowerRatio_pad1->XtoPixel(padInvYieldLowerRatio_pad1->GetX2()) < padInvYieldLowerRatio_pad1->YtoPixel(padInvYieldLowerRatio_pad1->GetY1())){
            textsizeLabelsXSecDown = (Double_t)textSizeLabelsPixel4pad/padInvYieldLowerRatio_pad1->XtoPixel(padInvYieldLowerRatio_pad1->GetX2()) ;
            textsizeFacXSecDown = (Double_t)1./padInvYieldLowerRatio_pad1->XtoPixel(padInvYieldLowerRatio_pad1->GetX2()) ;
        } else {
            textsizeLabelsXSecDown = (Double_t)textSizeLabelsPixel4pad/padInvYieldLowerRatio_pad1->YtoPixel(padInvYieldLowerRatio_pad1->GetY1());
            textsizeFacXSecDown = (Double_t)1./padInvYieldLowerRatio_pad1->YtoPixel(padInvYieldLowerRatio_pad1->GetY1());
        }

        if (padInvYieldUpperRatio_pad1->XtoPixel(padInvYieldUpperRatio_pad1->GetX2()) < padInvYieldUpperRatio_pad1->YtoPixel(padInvYieldUpperRatio_pad1->GetY1())){
            textsizeLabelsXSecMiddle = (Double_t)textSizeLabelsPixel4pad/padInvYieldUpperRatio_pad1->XtoPixel(padInvYieldUpperRatio_pad1->GetX2()) ;
            textsizeFacXSecMiddle = (Double_t)1./padInvYieldUpperRatio_pad1->XtoPixel(padInvYieldUpperRatio_pad1->GetX2()) ;
        } else {
            textsizeLabelsXSecMiddle = (Double_t)textSizeLabelsPixel4pad/padInvYieldUpperRatio_pad1->YtoPixel(padInvYieldUpperRatio_pad1->GetY1());
            textsizeFacXSecMiddle = (Double_t)1./padInvYieldUpperRatio_pad1->YtoPixel(padInvYieldUpperRatio_pad1->GetY1());
        }

        
		if(MesonInput){
	
			TGraphAsymmErrors *graphInvYieldPi0Stat_0010 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldPi0Comb2760GeVLHC11hStatErr_0010");
			TGraphAsymmErrors *graphInvYieldPi0Syst_0010= (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldPi0Comb2760GeVLHC11hSysErr_0010");

			TGraphAsymmErrors *graphInvYieldPi0Stat_2050 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldPi0Comb2760GeVLHC11hStatErr_2050");
			TGraphAsymmErrors *graphInvYieldPi0Syst_2050 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldPi0Comb2760GeVLHC11hSysErr_2050");

			TGraphAsymmErrors *graphInvYieldEtaStat_0010 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldEtaComb2760GeVLHC11hStatErr_0010");
			TGraphAsymmErrors *graphInvYieldEtaSyst_0010= (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldEtaComb2760GeVLHC11hSysErr_0010");

			TGraphAsymmErrors *graphInvYieldEtaStat_2050 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldEtaComb2760GeVLHC11hStatErr_2050");
			TGraphAsymmErrors *graphInvYieldEtaSyst_2050 = (TGraphAsymmErrors*)MesonInput->Get("graphInvYieldEtaComb2760GeVLHC11hSysErr_2050");

			TF1 *fitPi0InvYield2760GeVLHC11h_0010 = (TF1*)MesonInput->Get("FitToYieldPi0_0010");
			TF1 *fitPi0InvYield2760GeVLHC11h_2050 = (TF1*)MesonInput->Get("FitToYieldPi0_2050");
			TF1 *fitEtaInvYield2760GeVLHC11h_0010 = (TF1*)MesonInput->Get("FitToYieldEta_0010");
			TF1 *fitEtaInvYield2760GeVLHC11h_2050 = (TF1*)MesonInput->Get("FitToYieldEta_2050");
			
			TH1D *histoEPOSPi0Rebin_0010 = (TH1D*)MesonInput->Get("EPOSpredictionPi0_0010");
			TH1D *histoEPOSPi0Rebin_2050 = (TH1D*)MesonInput->Get("EPOSpredictionPi0_2050");
			TH1D *histoEPOSEtaRebin_0010 = (TH1D*)MesonInput->Get("EPOSpredictionEta_0010");
			TH1D *histoEPOSEtaRebin_2050 = (TH1D*)MesonInput->Get("EPOSpredictionEta_2050");

			
			if(graphInvYieldPi0Stat_0010 && graphInvYieldEtaStat_0010 && fitPi0InvYield2760GeVLHC11h_0010 && fitEtaInvYield2760GeVLHC11h_0010){

				TGraphAsymmErrors* graphRatioPi0InvYieldFitStat2760GeVLHC11h_0010 	= (TGraphAsymmErrors*)graphInvYieldPi0Stat_0010->Clone();
				graphRatioPi0InvYieldFitStat2760GeVLHC11h_0010 						= CalculateGraphErrRatioToFit(graphRatioPi0InvYieldFitStat2760GeVLHC11h_0010, fitPi0InvYield2760GeVLHC11h_0010); 
				TGraphAsymmErrors* graphRatioPi0InvYieldFitSyst2760GeVLHC11h_0010 	= (TGraphAsymmErrors*)graphInvYieldPi0Syst_0010->Clone();
				graphRatioPi0InvYieldFitSyst2760GeVLHC11h_0010 						= CalculateGraphErrRatioToFit(graphRatioPi0InvYieldFitSyst2760GeVLHC11h_0010, fitPi0InvYield2760GeVLHC11h_0010); 

				DrawGammaSetMarkerTGraphAsym(graphRatioPi0InvYieldFitSyst2760GeVLHC11h_0010,markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
				DrawGammaSetMarkerTGraphAsym(graphRatioPi0InvYieldFitStat2760GeVLHC11h_0010, markerStyleComb, markerSizeComb, kBlack, kBlack);

				TGraphAsymmErrors* graphRatioPi0InvYieldFitStat2760GeVLHC11h_2050 	= (TGraphAsymmErrors*)graphInvYieldPi0Stat_2050->Clone();
				graphRatioPi0InvYieldFitStat2760GeVLHC11h_2050 						= CalculateGraphErrRatioToFit(graphRatioPi0InvYieldFitStat2760GeVLHC11h_2050, fitPi0InvYield2760GeVLHC11h_2050); 
				TGraphAsymmErrors* graphRatioPi0InvYieldFitSyst2760GeVLHC11h_2050 	= (TGraphAsymmErrors*)graphInvYieldPi0Syst_2050->Clone();
				graphRatioPi0InvYieldFitSyst2760GeVLHC11h_2050 						= CalculateGraphErrRatioToFit(graphRatioPi0InvYieldFitSyst2760GeVLHC11h_2050, fitPi0InvYield2760GeVLHC11h_2050); 

				DrawGammaSetMarkerTGraphAsym(graphRatioPi0InvYieldFitSyst2760GeVLHC11h_2050, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
				DrawGammaSetMarkerTGraphAsym(graphRatioPi0InvYieldFitStat2760GeVLHC11h_2050, markerStyleComb, markerSizeComb, kBlack, kBlack);

				TGraphAsymmErrors* graphRatioEtaInvYieldFitStat2760GeVLHC11h_0010 	= (TGraphAsymmErrors*)graphInvYieldEtaStat_0010->Clone();
				graphRatioEtaInvYieldFitStat2760GeVLHC11h_0010 						= CalculateGraphErrRatioToFit(graphRatioEtaInvYieldFitStat2760GeVLHC11h_0010, fitEtaInvYield2760GeVLHC11h_0010); 
				TGraphAsymmErrors* graphRatioEtaInvYieldFitSyst2760GeVLHC11h_0010 	= (TGraphAsymmErrors*)graphInvYieldEtaSyst_0010->Clone();
				graphRatioEtaInvYieldFitSyst2760GeVLHC11h_0010 						= CalculateGraphErrRatioToFit(graphRatioEtaInvYieldFitSyst2760GeVLHC11h_0010, fitEtaInvYield2760GeVLHC11h_0010); 

				DrawGammaSetMarkerTGraphAsym(graphRatioEtaInvYieldFitSyst2760GeVLHC11h_0010, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
				DrawGammaSetMarkerTGraphAsym(graphRatioEtaInvYieldFitStat2760GeVLHC11h_0010, markerStyleComb, markerSizeComb, kBlack, kBlack);

				TGraphAsymmErrors* graphRatioEtaInvYieldFitStat2760GeVLHC11h_2050 	= (TGraphAsymmErrors*)graphInvYieldEtaStat_2050->Clone();
				graphRatioEtaInvYieldFitStat2760GeVLHC11h_2050 						= CalculateGraphErrRatioToFit(graphRatioEtaInvYieldFitStat2760GeVLHC11h_2050, fitEtaInvYield2760GeVLHC11h_2050); 
				TGraphAsymmErrors* graphRatioEtaInvYieldFitSyst2760GeVLHC11h_2050 	= (TGraphAsymmErrors*)graphInvYieldEtaSyst_2050->Clone();
				graphRatioEtaInvYieldFitSyst2760GeVLHC11h_2050 						= CalculateGraphErrRatioToFit(graphRatioEtaInvYieldFitSyst2760GeVLHC11h_2050, fitEtaInvYield2760GeVLHC11h_2050); 

				DrawGammaSetMarkerTGraphAsym(graphRatioEtaInvYieldFitSyst2760GeVLHC11h_2050, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
				DrawGammaSetMarkerTGraphAsym(graphRatioEtaInvYieldFitStat2760GeVLHC11h_2050, markerStyleComb, markerSizeComb, kBlack, kBlack);

				
				TGraphAsymmErrors* graphRatioLowPtPi0_0010 = (TGraphAsymmErrors*)TheoryCracowPi0LowPt_0010->Clone();	
				graphRatioLowPtPi0_0010 = CalculateGraphErrRatioToFit(graphRatioLowPtPi0_0010, fitPi0InvYield2760GeVLHC11h_0010); 
				DrawGammaSetMarkerTGraphAsym(graphRatioLowPtPi0_0010, 0, 0, colorCracowRatio, colorCracowRatio, 5, kTRUE, colorCracowRatio);

				TH1D* histoRatioEPOSToFitPi0_0010 = (TH1D*) histoEPOSPi0Rebin_0010->Clone();     
				histoRatioEPOSToFitPi0_0010 = CalculateHistoRatioToFit (histoRatioEPOSToFitPi0_0010, fitPi0InvYield2760GeVLHC11h_0010); 
				DrawGammaSetMarker(histoRatioEPOSToFitPi0_0010, 24, 1.5,colorEPOSRatio,colorEPOSRatio);  
				histoRatioEPOSToFitPi0_0010->SetLineWidth(widthCommonFit);

				TGraphAsymmErrors* graphRatioLowPtPi0_2050 = (TGraphAsymmErrors*)TheoryCracowPi0LowPt_2050->Clone();	
				graphRatioLowPtPi0_2050 = CalculateGraphErrRatioToFit(graphRatioLowPtPi0_2050, fitPi0InvYield2760GeVLHC11h_2050); 
				DrawGammaSetMarkerTGraphAsym(graphRatioLowPtPi0_2050, 0, 0, colorCracowRatio, colorCracowRatio, 5, kTRUE, colorCracowRatio);

				TH1D* histoRatioEPOSToFitPi0_2050 = (TH1D*) histoEPOSPi0Rebin_2050->Clone();     
				histoRatioEPOSToFitPi0_2050 = CalculateHistoRatioToFit (histoRatioEPOSToFitPi0_2050, fitPi0InvYield2760GeVLHC11h_2050); 
				DrawGammaSetMarker(histoRatioEPOSToFitPi0_2050, 24, 1.5,colorEPOSRatio,colorEPOSRatio);  
				histoRatioEPOSToFitPi0_2050->SetLineWidth(widthCommonFit);
				
				TGraphAsymmErrors* graphRatioLowPtEta_0010 = (TGraphAsymmErrors*)TheoryCracowEtaLowPt_0010->Clone();	
				graphRatioLowPtEta_0010 = CalculateGraphErrRatioToFit(graphRatioLowPtEta_0010, fitEtaInvYield2760GeVLHC11h_0010); 
				DrawGammaSetMarkerTGraphAsym(graphRatioLowPtEta_0010, 0, 0, colorCracowRatio, colorCracowRatio, 5, kTRUE, colorCracowRatio);

				TH1D* histoRatioEPOSToFitEta_0010 = (TH1D*) histoEPOSEtaRebin_0010->Clone();     
				histoRatioEPOSToFitEta_0010 = CalculateHistoRatioToFit (histoRatioEPOSToFitEta_0010, fitEtaInvYield2760GeVLHC11h_0010); 
				DrawGammaSetMarker(histoRatioEPOSToFitEta_0010, 24, 1.5,colorEPOSRatio,colorEPOSRatio);  
				histoRatioEPOSToFitEta_0010->SetLineWidth(widthCommonFit);

				TGraphAsymmErrors* graphRatioLowPtEta_2050 = (TGraphAsymmErrors*)TheoryCracowEtaLowPt_2050->Clone();	
				graphRatioLowPtEta_2050 = CalculateGraphErrRatioToFit(graphRatioLowPtEta_2050, fitEtaInvYield2760GeVLHC11h_2050); 
				DrawGammaSetMarkerTGraphAsym(graphRatioLowPtEta_2050, 0, 0, colorCracowRatio, colorCracowRatio, 5, kTRUE, colorCracowRatio);

				TH1D* histoRatioEPOSToFitEta_2050 = (TH1D*) histoEPOSEtaRebin_2050->Clone();     
				histoRatioEPOSToFitEta_2050 = CalculateHistoRatioToFit (histoRatioEPOSToFitEta_2050, fitEtaInvYield2760GeVLHC11h_2050); 
				DrawGammaSetMarker(histoRatioEPOSToFitEta_2050, 24, 1.5,colorEPOSRatio,colorEPOSRatio);  
				histoRatioEPOSToFitEta_2050->SetLineWidth(widthCommonFit);

				padInvYieldUpperRatio_pad1->cd();
				padInvYieldUpperRatio_pad1->SetLogx(1);
		        TH2F * ratio2DUpPad1OnlyRatio = new TH2F("ratio2DUpPad1OnlyRatio","ratio2DUpPad1OnlyRatio",1000,0.4,40.,1000,0.01,2.1);
				SetStyleHistoTH2ForGraphs(ratio2DUpPad1OnlyRatio, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle, 
										  0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.71, 510, 505);
				ratio2DUpPad1OnlyRatio->GetXaxis()->SetRangeUser(0.75,29.);
				ratio2DUpPad1OnlyRatio->GetXaxis()->SetMoreLogLabels(kTRUE);
				ratio2DUpPad1OnlyRatio->GetXaxis()->SetNoExponent(kTRUE);
				ratio2DUpPad1OnlyRatio->GetYaxis()->SetLabelFont(42);
				ratio2DUpPad1OnlyRatio->GetYaxis()->SetLabelOffset(+0.01);
				ratio2DUpPad1OnlyRatio->DrawCopy();

					DrawGammaLines(0.75,29.,1., 1.,0.1,kGray);
					
					for(Int_t a=0; a<4; a++){
						graphRatioLowPtPi0_0010->RemovePoint(0);
					}
					histoRatioEPOSToFitPi0_0010->GetXaxis()->SetRangeUser(1.,13.);
					graphRatioLowPtPi0_0010->Draw("l,same");
					histoRatioEPOSToFitPi0_0010->Draw("same,hist,c");
					
					if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRatioPi0InvYieldFitStat2760GeVLHC11h_0010);
					graphRatioPi0InvYieldFitStat2760GeVLHC11h_0010->Draw("p,same,e0");
					graphRatioPi0InvYieldFitSyst2760GeVLHC11h_0010->Draw("E2same");
					
					
// 					TLatex *labelcent0010 = new TLatex(0.85,0.88," 0-10%");
// 					SetStyleTLatex( labelcent0010,0.85* textsizeLabelsXSecMiddle,4);
// 					labelcent0010->Draw();
					
					TLegend* legen4PadOnlyRatios_pad1 = GetAndSetLegend2(0.52, 0.22-2*0.08, 0.9, 0.22, 0.85*textSizeLabelsPixel4pad);
					legen4PadOnlyRatios_pad1->SetMargin(0.17);
					legen4PadOnlyRatios_pad1->SetHeader(collisionSystemPbPb0010.Data());
					legen4PadOnlyRatios_pad1->AddEntry(graphRatioPi0InvYieldFitSyst2760GeVLHC11h_0010,"#pi^{0} ALICE","pf");
					legen4PadOnlyRatios_pad1->Draw();
					
					TLegend* legen4PadOnlyRatiosTheory = GetAndSetLegend2(0.2, 0.2-2*0.08, 0.5, 0.2, 0.85* textSizeLabelsPixel4pad);
					legen4PadOnlyRatiosTheory->AddEntry(graphRatioCombLowPt2760GeVLHC11h_0010,"NEQ SHM","l");
					legen4PadOnlyRatiosTheory->AddEntry(histoRatioEPOSToFit2760GeVLHC11h_0010,"EPOS","l");
					legen4PadOnlyRatiosTheory->Draw();

					
					
				padInvYieldLowerRatio_pad1->cd();
				padInvYieldLowerRatio_pad1->SetLogx(1);
		        TH2F * ratio2DLowerPad1OnlyRatio = new TH2F("ratio2DLowerPad1OnlyRatio","ratio2DLowerPad1OnlyRatio",1000,0.4,40.,1000,0.01,2.1);
				SetStyleHistoTH2ForGraphs(ratio2DLowerPad1OnlyRatio, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown, 
										  0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 1,0.85, 510, 505);
				ratio2DLowerPad1OnlyRatio->GetXaxis()->SetRangeUser(0.75,29);
				ratio2DLowerPad1OnlyRatio->GetXaxis()->SetMoreLogLabels(kTRUE);
				ratio2DLowerPad1OnlyRatio->GetXaxis()->SetNoExponent(kTRUE);
				ratio2DLowerPad1OnlyRatio->GetYaxis()->SetLabelFont(42);
				ratio2DLowerPad1OnlyRatio->GetYaxis()->SetLabelOffset(+0.01);
				ratio2DLowerPad1OnlyRatio->DrawCopy();

					DrawGammaLines(0.75,29.,1., 1.,0.1,kGray);
					
					graphRatioLowPtEta_0010->RemovePoint(0);
					graphRatioLowPtEta_0010->RemovePoint(0);
					graphRatioLowPtEta_0010->RemovePoint(0);
					graphRatioLowPtEta_0010->RemovePoint(0);
					graphRatioLowPtEta_0010->Draw("l,same");
					histoRatioEPOSToFitEta_0010->GetXaxis()->SetRangeUser(1.,13.);					
					histoRatioEPOSToFitEta_0010->Draw("same,hist,c");
					
					if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRatioEtaInvYieldFitStat2760GeVLHC11h_0010);
					graphRatioEtaInvYieldFitStat2760GeVLHC11h_0010->Draw("p,same,e0");
					graphRatioEtaInvYieldFitSyst2760GeVLHC11h_0010->Draw("E2same");
					
					TLatex *labelcent0010a = new TLatex(0.77,0.88," 0-10%");
					SetStyleTLatex( labelcent0010a, 0.85*textsizeLabelsXSecDown,4);
// 					labelcent0010a->Draw();

					TLegend* legen4PadOnlyRatios_pad2 = GetAndSetLegend2(0.52, 0.38-2*0.08, 0.9, 0.38, 0.85* textSizeLabelsPixel4pad);
					legen4PadOnlyRatios_pad2->SetMargin(0.17);
					legen4PadOnlyRatios_pad2->SetHeader(collisionSystemPbPb0010.Data());
					legen4PadOnlyRatios_pad2->AddEntry(graphRatioEtaInvYieldFitSyst2760GeVLHC11h_0010,"#eta ALICE","pf"); // \\sqrt{s_\\mathrm{nn}}
					legen4PadOnlyRatios_pad2->Draw();


				padInvYieldUpperRatio_pad2->cd();
				padInvYieldUpperRatio_pad2->SetLogx(1);
				ratio2DUpPad1OnlyRatio->DrawCopy();

				// 		DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat2760GeVLHC11h_2050, markerStyle2050, markerSizeComb*2, kBlack, kBlack, widthLinesBoxes, kFALSE);
				// 		DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys2760GeVLHC11h_2050, markerStyle2050, markerSizeComb*2, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
					DrawGammaLines(0.75,29,1., 1.,0.1,kGray);

					graphRatioLowPtPi0_2050->RemovePoint(0);
					graphRatioLowPtPi0_2050->RemovePoint(0);
					graphRatioLowPtPi0_2050->RemovePoint(0);
					graphRatioLowPtPi0_2050->RemovePoint(0);
					graphRatioLowPtPi0_2050->Draw("l,same");
					histoRatioEPOSToFitPi0_2050->GetXaxis()->SetRangeUser(1.,13.);
					histoRatioEPOSToFitPi0_2050->Draw("same,hist,c");
					
					if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRatioPi0InvYieldFitStat2760GeVLHC11h_2050);
					graphRatioPi0InvYieldFitStat2760GeVLHC11h_2050->Draw("p,same,e0");
					graphRatioPi0InvYieldFitSyst2760GeVLHC11h_2050->Draw("E2same");
					
					TLatex *labelcent2050 = new TLatex(0.85,0.88,"20-50%");
					SetStyleTLatex( labelcent2050, 0.85*textsizeLabelsXSecMiddle,4);
// 					labelcent2050->Draw();

					TLegend* legen4PadOnlyRatios_pad3 = GetAndSetLegend2(0.39, 0.22-2*0.08, 0.82, 0.22, 0.85* textSizeLabelsPixel4pad);
					legen4PadOnlyRatios_pad3->SetMargin(0.17);
					legen4PadOnlyRatios_pad3->SetHeader(collisionSystemPbPb2050.Data());
					legen4PadOnlyRatios_pad3->AddEntry(graphRatioPi0InvYieldFitSyst2760GeVLHC11h_2050,"#pi^{0} ALICE","pf");
					legen4PadOnlyRatios_pad3->Draw();


				padInvYieldLowerRatio_pad2->cd();
				padInvYieldLowerRatio_pad2->SetLogx(1);
				ratio2DLowerPad1OnlyRatio->DrawCopy();
				
					DrawGammaLines(0.75,29.,1., 1.,0.1,kGray);

					graphRatioLowPtEta_2050->RemovePoint(0);
					graphRatioLowPtEta_2050->RemovePoint(0);
					graphRatioLowPtEta_2050->RemovePoint(0);
					graphRatioLowPtEta_2050->RemovePoint(0);
					graphRatioLowPtEta_2050->Draw("l,same");
					histoRatioEPOSToFitEta_2050->GetXaxis()->SetRangeUser(1.,13.);
					histoRatioEPOSToFitEta_2050->Draw("same,hist,c");
				
					if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphRatioEtaInvYieldFitStat2760GeVLHC11h_2050);
					graphRatioEtaInvYieldFitStat2760GeVLHC11h_2050->Draw("p,same,e0");
					graphRatioEtaInvYieldFitSyst2760GeVLHC11h_2050->Draw("E2same");
					
					TLatex *labelcent2050a = new TLatex(0.77,0.88,"20-50%");
					SetStyleTLatex( labelcent2050a, 0.85*textsizeLabelsXSecDown,4);
// 					labelcent2050a->Draw();

					TLegend* legen4PadOnlyRatios_pad4 = GetAndSetLegend2(0.39, 0.38-2*0.08, 0.82, 0.38, 0.85* textSizeLabelsPixel4pad);
					legen4PadOnlyRatios_pad4->SetMargin(0.17);
					legen4PadOnlyRatios_pad4->SetHeader(collisionSystemPbPb2050.Data());
					legen4PadOnlyRatios_pad4->AddEntry(graphRatioEtaInvYieldFitSyst2760GeVLHC11h_2050,"#eta ALICE","pf");
					legen4PadOnlyRatios_pad4->Draw();

			}
		}
    canvasInvYieldSectionOnlyRatios4pad->Update();
    canvasInvYieldSectionOnlyRatios4pad->SaveAs(Form("%s/YieldCombinedLHC11h_DataWithModelsOnlyRatios4Pad.%s",outputDir.Data(),suffix.Data()));
    canvasInvYieldSectionOnlyRatios4pad->SaveAs(Form("%s/YieldCombinedLHC11h_DataWithModelsOnlyRatios4Pad.%s",paperPlots.Data(),suffix.Data()));   
    
	
    
	//*********************************************************************************************************************//
	//*************************	 Combination of RAA 	*******************************************************************//
	
//     if(bWCorrection.CompareTo("X")==0 ){
//       cout << " \n\nCalculating RAA for " << meson.Data() << endl;		
//       
//       TF1* fitInvCrossSectionPi0Comb2760GeV = FitObject("l","fitInvCrossSectionPi0Comb2760GeV","Pi0");
//       fitInvCrossSectionPi0Comb2760GeV->SetParameter(1,7.5);
//       fitInvCrossSectionPi0Comb2760GeV->SetParameter(2,0.152);
//       fitInvCrossSectionPi0Comb2760GeV->SetRange(0,15);
// 
//       graphInvSectionComb->Fit(fitInvCrossSectionPi0Comb2760GeV,"QNRMEX0+","",0.4,8.);
//       cout << "Initial Fit result" << endl;
//       cout << WriteParameterToFile(fitInvCrossSectionPi0Comb2760GeV)<< endl;

//     TGraphAsymmErrors* graphPi0RAAPCM0010;
//     TGraphAsymmErrors* graphPi0RAASysPCM0010;
// 
//     TCanvas* canvasDummy5 = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
//     DrawGammaCanvasSettings( canvasDummy5,  0.1, 0.01, 0.015, 0.08);
//     canvasDummy5->SetLogy();
//     canvasDummy5->SetLogx();
//     TH2F * histo2DDummy5;
//     histo2DDummy5 = new TH2F("histo2DDummy5","histo2DDummy5",1000,0.23,30.,1000,1e-8,10);
//     SetStyleHistoTH2ForGraphs(histo2DDummy5, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
//     histo2DDummy5->DrawCopy(); 
// 
//     DrawGammaSetMarkerTGraphAsym(graphInvSectionCombPi02760GeV, 20,markerSizeCommonSpectrumPi0900GeV, kRed, kRed, widthLinesBoxes, kTRUE);
//     graphInvSectionCombPi02760GeV->Draw("pEsame");
// 
//     fitInvCrossSectionPi0Comb2760GeV->SetLineColor(kBlue+2);
//     fitInvCrossSectionPi0Comb2760GeV->Draw("same");
// 
//     canvasDummy5->Update();
//     canvasDummy5->Print(Form("%s/PPDataPlusFit.%s",outputDir.Data(),suffix.Data()));
      
//       cout << "____________ PCM ____________" << endl ;
// 
//       CalcRaa(    graphInvSectionCombStatPi02760GeVPlot, graphInvSectionCombSysEta2760GeVPlot,graphInvSectionCombStatPi02760GeVPlot, fitInvCrossSectionPi0Comb2760GeV,
//                     graphPCMInvYieldStatPbPb2760GeVYShifted_0010, graphPCMInvYieldSysPbPb2760GeVYShifted_0010,  //PbPb Yields
//                     &graphPi0RAAPCM0010, &graphPi0RAASysPCM0010,
//                     nColl0010, nCollErr0010,8.,0);
// 
// //     //  return;
// //     cout << "EMCal*********************************************" << endl ;
// //     cout << "***********************************************" << endl<< endl ;
// // 
// //     TGraphAsymmErrors* graphRAAEMCal0010;
// //     TGraphAsymmErrors* graphRAASysEMCal0010;
// //     CalcRaa(    graphInvSectionEMCalStatPi02760GeVRed, graphInvSectionEMCalSysPi02760GeVRed,graphInvSectionCombPi02760GeVOnlyStat, fitInvCrossSectionPi0Comb2760GeV,
// //                     graphYieldPi0EMCalPbPb0010StatErrYShifted, graphYieldPi0EMCalPbPb0010SysRAAErrYShifted,  //PbPb Yields
// //                     &graphRAAEMCal0010, &graphRAASysEMCal0010,
// //                     nColl0010, nCollErr0010, 8.,1);
// // 
// //     //  return; 
// //     cout << "Combined*****************************************" << endl ;
// //     cout << "***********************************************" << endl<< endl ;
// // 
// //     TGraphAsymmErrors* graphRAACombInd0010;
// //     TGraphAsymmErrors* graphRAASysCombInd0010;
// //     TGraphAsymmErrors* graphRAAPi0CombPbPb0010_2 = CombinePtPointsRAA(  graphRAAPCM0010,        graphRAASysPCM0010,
// //                                                 graphRAAEMCal0010,   graphRAASysEMCal0010,
// //                                                 graphRAACombInd0010, graphRAASysCombInd0010,
// //                                                 xPtLimitsPbPbEMCal2, 19, 0, 0,2);
// 
//       
//     }
    
	TGraphAsymmErrors* graphCombRAAStatPbPb2760GeVLHC11h_0010 = NULL;
	TGraphAsymmErrors* graphCombRAASysPbPb2760GeVLHC11h_0010 = NULL;
	TGraphAsymmErrors* graphCombRAATot2760GeVLHC11h_0010 = NULL;

	TGraphAsymmErrors* graphCombRAAStatPbPb2760GeVLHC11h_2050 = NULL;
	TGraphAsymmErrors* graphCombRAASysPbPb2760GeVLHC11h_2050 = NULL;
	TGraphAsymmErrors* graphCombRAATot2760GeVLHC11h_2050 = NULL;
	
    cout << " \n\nCombining RAA for " << meson.Data() << endl;        
//     if(bWCorrection.CompareTo("X")==0 ){
    if(kFALSE ){
      if(meson.CompareTo("Pi0")==0){
          graphCombRAATot2760GeVLHC11h_0010 = CombinePtPointsRAA(graphPCMPi0RAAStatPbPb2760GeVYShifted_0010,     graphPCMPi0RAASysPbPb2760GeVYShifted_0010,     
                                                                  graphEMCalPi0RAAStatPbPb2760GeVYShifted_0010, graphEMCalPi0RAASysPbPb2760GeVYShifted_0010,                                graphCombRAAStatPbPb2760GeVLHC11h_0010, graphCombRAASysPbPb2760GeVLHC11h_0010,
                                                                  xPtLimitsPi0, 24, 0, 1, 16, kFALSE);

          graphCombRAATot2760GeVLHC11h_2050 = CombinePtPointsRAA(graphPCMPi0RAAStatPbPb2760GeVYShifted_2050,     graphPCMPi0RAASysPbPb2760GeVYShifted_2050,     
                                                                  graphEMCalPi0RAAStatPbPb2760GeVYShifted_2050, graphEMCalPi0RAASysPbPb2760GeVYShifted_2050,                                graphCombRAAStatPbPb2760GeVLHC11h_2050, graphCombRAASysPbPb2760GeVLHC11h_2050,
                                                                  xPtLimitsPi0, 24, 0, 1, 16, kFALSE);
      
      } else if(meson.CompareTo("Eta")==0){
          
          graphCombRAATot2760GeVLHC11h_0010 = CombinePtPointsRAA(graphPCMEtaRAAStatPbPb2760GeVYShifted_0010,     graphPCMEtaRAASysPbPb2760GeVYShifted_0010,     
                                                                  graphEMCalEtaRAAStatPbPb2760GeVYShifted_0010, graphEMCalEtaRAASysPbPb2760GeVYShifted_0010,                                graphCombRAAStatPbPb2760GeVLHC11h_0010, graphCombRAASysPbPb2760GeVLHC11h_0010,
                                                                  xPtLimitsEta, 14, 0, 2, 6, kFALSE);

          graphCombRAATot2760GeVLHC11h_2050 = CombinePtPointsRAA(graphPCMEtaRAAStatPbPb2760GeVYShifted_2050,     graphPCMEtaRAASysPbPb2760GeVYShifted_2050,     
                                                                  graphEMCalEtaRAAStatPbPb2760GeVYShifted_2050, graphEMCalEtaRAASysPbPb2760GeVYShifted_2050,                                graphCombRAAStatPbPb2760GeVLHC11h_2050, graphCombRAASysPbPb2760GeVLHC11h_2050,
                                                                  xPtLimitsEta, 14, 0, 2, 6, kFALSE);
          
      }
    } else {
      
      if(meson.CompareTo("Pi0")==0){
          graphCombRAATot2760GeVLHC11h_0010 = CombinePtPointsRAA(graphPCMPi0RAAStatPbPb2760GeV_0010,     graphPCMPi0RAASysPbPb2760GeV_0010,     
                                                                  graphEMCalPi0RAAStatPbPb2760GeV_0010, graphEMCalPi0RAASysPbPb2760GeV_0010,                                graphCombRAAStatPbPb2760GeVLHC11h_0010, graphCombRAASysPbPb2760GeVLHC11h_0010,
                                                                  xPtLimitsPi0, 24, 0, 1, 16, kFALSE);

          graphCombRAATot2760GeVLHC11h_2050 = CombinePtPointsRAA(graphPCMPi0RAAStatPbPb2760GeV_2050,     graphPCMPi0RAASysPbPb2760GeV_2050,     
                                                                  graphEMCalPi0RAAStatPbPb2760GeV_2050, graphEMCalPi0RAASysPbPb2760GeV_2050,                                graphCombRAAStatPbPb2760GeVLHC11h_2050, graphCombRAASysPbPb2760GeVLHC11h_2050,
                                                                  xPtLimitsPi0, 24, 0, 1, 16, kFALSE);
      
      } else if(meson.CompareTo("Eta")==0){
          
          graphCombRAATot2760GeVLHC11h_0010 = CombinePtPointsRAA(graphPCMEtaRAAStatPbPb2760GeV_0010,     graphPCMEtaRAASysPbPb2760GeV_0010,     
                                                                  graphEMCalEtaRAAStatPbPb2760GeV_0010, graphEMCalEtaRAASysPbPb2760GeV_0010,                                graphCombRAAStatPbPb2760GeVLHC11h_0010, graphCombRAASysPbPb2760GeVLHC11h_0010,
                                                                  xPtLimitsEta, 14, 0, 2, 6, kFALSE);

          graphCombRAATot2760GeVLHC11h_2050 = CombinePtPointsRAA(graphPCMEtaRAAStatPbPb2760GeV_2050,     graphPCMEtaRAASysPbPb2760GeV_2050,     
                                                                  graphEMCalEtaRAAStatPbPb2760GeV_2050, graphEMCalEtaRAASysPbPb2760GeV_2050,                                graphCombRAAStatPbPb2760GeVLHC11h_2050, graphCombRAASysPbPb2760GeVLHC11h_2050,
                                                                  xPtLimitsEta, 14, 0, 2, 6, kFALSE);
          
      }
    }      
      
      
		if(PaperPi0 && meson.CompareTo("Pi0")==0){
			graphCombRAAStatPbPb2760GeVLHC11h_0010->RemovePoint(0);
			graphCombRAASysPbPb2760GeVLHC11h_0010->RemovePoint(0);

			graphCombRAAStatPbPb2760GeVLHC11h_0010->RemovePoint(0);
			graphCombRAASysPbPb2760GeVLHC11h_0010->RemovePoint(0);

			graphCombRAAStatPbPb2760GeVLHC11h_0010->RemovePoint(0);
			graphCombRAASysPbPb2760GeVLHC11h_0010->RemovePoint(0);

			graphCombRAASysPbPb2760GeVLHC11h_2050->RemovePoint(0);
			graphCombRAAStatPbPb2760GeVLHC11h_2050->RemovePoint(0);

			graphCombRAASysPbPb2760GeVLHC11h_2050->RemovePoint(0);
			graphCombRAAStatPbPb2760GeVLHC11h_2050->RemovePoint(0);
			
			graphCombRAASysPbPb2760GeVLHC11h_2050->RemovePoint(0);
			graphCombRAAStatPbPb2760GeVLHC11h_2050->RemovePoint(0);

		}

	TCanvas* canvasRAAcomboPi0andEta = new TCanvas("canvasRAAcomboPi0andEta","",200,10,1200,1100);  //200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasRAAcomboPi0andEta, 0.08, 0.02, 0.035, 0.09);
// 		canvasRAAcomboPi0andEta->SetLogx();
	
		TH2F * histo2DRAAcomboPi0andEta;
		histo2DRAAcomboPi0andEta = new TH2F("histo2DRAAcomboPi0andEta","histo2DRAAcomboPi0andEta",11000,0.23,70.,1000,-0.5,2.);
		SetStyleHistoTH2ForGraphs(histo2DRAAcomboPi0andEta, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}",0.035,0.04, 0.035,0.04, 1.,0.9);
		histo2DRAAcomboPi0andEta->GetYaxis()->SetRangeUser(0.,1.4);
		histo2DRAAcomboPi0andEta->GetXaxis()->SetRangeUser(0.,21);
		histo2DRAAcomboPi0andEta->Draw("copy");


		TGraphAsymmErrors* graphCombRAAStatPbPb2760GeVLHC11hPi0_0010 = NULL;
		TGraphAsymmErrors* graphCombRAASysPbPb2760GeVLHC11hPi0_0010 = NULL;
		TGraphAsymmErrors* graphCombRAATot2760GeVLHC11hPi0_0010 = NULL;

		TGraphAsymmErrors* graphCombRAAStatPbPb2760GeVLHC11hPi0_2050 = NULL;
		TGraphAsymmErrors* graphCombRAASysPbPb2760GeVLHC11hPi0_2050 = NULL;
		TGraphAsymmErrors* graphCombRAATot2760GeVLHC11hPi0_2050 = NULL;

		TGraphAsymmErrors* graphCombRAAStatPbPb2760GeVLHC11hEta_0010 = NULL;
		TGraphAsymmErrors* graphCombRAASysPbPb2760GeVLHC11hEta_0010 = NULL;
		TGraphAsymmErrors* graphCombRAATot2760GeVLHC11hEta_0010 = NULL;

		TGraphAsymmErrors* graphCombRAAStatPbPb2760GeVLHC11hEta_2050 = NULL;
		TGraphAsymmErrors* graphCombRAASysPbPb2760GeVLHC11hEta_2050 = NULL;
		TGraphAsymmErrors* graphCombRAATot2760GeVLHC11hEta_2050 = NULL;
		
//     if(bWCorrection.CompareTo("X")==0 ){
    if(kFALSE ){

          graphCombRAATot2760GeVLHC11hPi0_0010 = CombinePtPointsRAA(graphPCMPi0RAAStatPbPb2760GeVYShifted_0010,	 graphPCMPi0RAASysPbPb2760GeVYShifted_0010, 	
                                                                      graphEMCalPi0RAAStatPbPb2760GeVYShifted_0010, graphEMCalPi0RAASysPbPb2760GeVYShifted_0010,								graphCombRAAStatPbPb2760GeVLHC11hPi0_0010, graphCombRAASysPbPb2760GeVLHC11hPi0_0010,
                                                                      xPtLimitsPi0, 24, 0, 1, 16, kFALSE);

          graphCombRAATot2760GeVLHC11hPi0_2050 = CombinePtPointsRAA(graphPCMPi0RAAStatPbPb2760GeVYShifted_2050,	 graphPCMPi0RAASysPbPb2760GeVYShifted_2050, 	
                                                                      graphEMCalPi0RAAStatPbPb2760GeVYShifted_2050, graphEMCalPi0RAASysPbPb2760GeVYShifted_2050,								graphCombRAAStatPbPb2760GeVLHC11hPi0_2050, graphCombRAASysPbPb2760GeVLHC11hPi0_2050,
                                                                      xPtLimitsPi0, 24, 0, 1, 16, kFALSE);
          
          graphCombRAATot2760GeVLHC11hEta_0010 = CombinePtPointsRAA(graphPCMEtaRAAStatPbPb2760GeVYShifted_0010,	 graphPCMEtaRAASysPbPb2760GeVYShifted_0010, 	
                                                                      graphEMCalEtaRAAStatPbPb2760GeVYShifted_0010, graphEMCalEtaRAASysPbPb2760GeVYShifted_0010,								graphCombRAAStatPbPb2760GeVLHC11hEta_0010, graphCombRAASysPbPb2760GeVLHC11hEta_0010,
                                                                      xPtLimitsEta, 14, 0, 2, 6, kFALSE);

          graphCombRAATot2760GeVLHC11hEta_2050 = CombinePtPointsRAA(graphPCMEtaRAAStatPbPb2760GeVYShifted_2050,	 graphPCMEtaRAASysPbPb2760GeVYShifted_2050, 	
                                                                      graphEMCalEtaRAAStatPbPb2760GeVYShifted_2050, graphEMCalEtaRAASysPbPb2760GeVYShifted_2050,								graphCombRAAStatPbPb2760GeVLHC11hEta_2050, graphCombRAASysPbPb2760GeVLHC11hEta_2050,
                                                                      xPtLimitsEta, 14, 0, 2, 6, kFALSE);
    } else {
         
          graphCombRAATot2760GeVLHC11hPi0_0010 = CombinePtPointsRAA(graphPCMPi0RAAStatPbPb2760GeV_0010,  graphPCMPi0RAASysPbPb2760GeV_0010,     
                                                                      graphEMCalPi0RAAStatPbPb2760GeV_0010, graphEMCalPi0RAASysPbPb2760GeV_0010,                                graphCombRAAStatPbPb2760GeVLHC11hPi0_0010, graphCombRAASysPbPb2760GeVLHC11hPi0_0010,
                                                                      xPtLimitsPi0, 24, 0, 1, 16, kFALSE);

          graphCombRAATot2760GeVLHC11hPi0_2050 = CombinePtPointsRAA(graphPCMPi0RAAStatPbPb2760GeV_2050,  graphPCMPi0RAASysPbPb2760GeV_2050,     
                                                                      graphEMCalPi0RAAStatPbPb2760GeV_2050, graphEMCalPi0RAASysPbPb2760GeV_2050,                                graphCombRAAStatPbPb2760GeVLHC11hPi0_2050, graphCombRAASysPbPb2760GeVLHC11hPi0_2050,
                                                                      xPtLimitsPi0, 24, 0, 1, 16, kFALSE);
          
          graphCombRAATot2760GeVLHC11hEta_0010 = CombinePtPointsRAA(graphPCMEtaRAAStatPbPb2760GeV_0010,  graphPCMEtaRAASysPbPb2760GeV_0010,     
                                                                      graphEMCalEtaRAAStatPbPb2760GeV_0010, graphEMCalEtaRAASysPbPb2760GeV_0010,                                graphCombRAAStatPbPb2760GeVLHC11hEta_0010, graphCombRAASysPbPb2760GeVLHC11hEta_0010,
                                                                      xPtLimitsEta, 14, 0, 2, 6, kFALSE);

          graphCombRAATot2760GeVLHC11hEta_2050 = CombinePtPointsRAA(graphPCMEtaRAAStatPbPb2760GeV_2050,  graphPCMEtaRAASysPbPb2760GeV_2050,     
                                                                      graphEMCalEtaRAAStatPbPb2760GeV_2050, graphEMCalEtaRAASysPbPb2760GeV_2050,                                graphCombRAAStatPbPb2760GeVLHC11hEta_2050, graphCombRAASysPbPb2760GeVLHC11hEta_2050,
                                                                      xPtLimitsEta, 14, 0, 2, 6, kFALSE);
          
    }
          
		if(PaperPi0 && meson.CompareTo("Pi0")==0){
			graphCombRAAStatPbPb2760GeVLHC11hPi0_0010->RemovePoint(0);
			graphCombRAASysPbPb2760GeVLHC11hPi0_0010->RemovePoint(0);

			graphCombRAAStatPbPb2760GeVLHC11hPi0_0010->RemovePoint(0);
			graphCombRAASysPbPb2760GeVLHC11hPi0_0010->RemovePoint(0);

			graphCombRAAStatPbPb2760GeVLHC11hPi0_0010->RemovePoint(0);
			graphCombRAASysPbPb2760GeVLHC11hPi0_0010->RemovePoint(0);

			graphCombRAAStatPbPb2760GeVLHC11hPi0_2050->RemovePoint(0);
			graphCombRAASysPbPb2760GeVLHC11hPi0_2050->RemovePoint(0);

			graphCombRAAStatPbPb2760GeVLHC11hPi0_2050->RemovePoint(0);
			graphCombRAASysPbPb2760GeVLHC11hPi0_2050->RemovePoint(0);
			
			graphCombRAAStatPbPb2760GeVLHC11hPi0_2050->RemovePoint(0);
			graphCombRAASysPbPb2760GeVLHC11hPi0_2050->RemovePoint(0);

		}
// 		TGraphAsymmErrors *graphChargedKaonRAAStatandSyst0010 = CalculateCombinedSysAndStatError(graphChargedKaonRAA0010,graphChargedKaonRAASys0010);
// 		DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAAStatandSyst0010, 25,1.6, colorCharged,colorCharged, 0.1, kFALSE);
		DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAA0010, 25,1.5, colorCharged,colorCharged, 0.1, kFALSE);
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphChargedKaonRAA0010);
		DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAASys0010,25,1.5, colorCharged,colorCharged, 1, kTRUE);
		graphChargedKaonRAASys0010->Draw("2same");
		graphChargedKaonRAA0010->Draw("p,same"); 
// 		graphChargedKaonRAAStatandSyst0010->Draw("p,same"); 
		
		//=============Combined Pi0
		DrawGammaSetMarkerTGraphAsym(graphCombRAASysPbPb2760GeVLHC11hPi0_0010, 24, markerSizeComb, colorPCM0010, colorPCM0010, widthLinesBoxes, kTRUE);
		graphCombRAASysPbPb2760GeVLHC11hPi0_0010->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphCombRAAStatPbPb2760GeVLHC11hPi0_0010, 24, markerSizeComb, colorPCM0010, colorPCM0010);
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPb2760GeVLHC11hPi0_0010);
		graphCombRAAStatPbPb2760GeVLHC11hPi0_0010->Draw("p,same,e0");

		//=============Combined Eta
		DrawGammaSetMarkerTGraphAsym(graphCombRAASysPbPb2760GeVLHC11hEta_0010, markerStyle0010, markerSizeComb, colorPCM0010, colorPCM0010, widthLinesBoxes, kTRUE);
		graphCombRAASysPbPb2760GeVLHC11hEta_0010->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphCombRAAStatPbPb2760GeVLHC11hEta_0010, markerStyle0010, markerSizeComb, colorPCM0010 , colorPCM0010);
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPb2760GeVLHC11hEta_0010);
		graphCombRAAStatPbPb2760GeVLHC11hEta_0010->Draw("p,same,e0");


// 		TLegend* legendRAAcomboPi0andEta0010 = new TLegend(0.33,0.78,0.8,0.9);
		TLegend* legendRAAcomboPi0andEta0010 = new TLegend(0.33,0.7,0.75,0.9);
		legendRAAcomboPi0andEta0010->SetFillColor(0);
		legendRAAcomboPi0andEta0010->SetLineColor(0);
		legendRAAcomboPi0andEta0010->SetTextFont(42);
// 		legendRAAcomboPi0andEta0010->SetMargin(0.17);
		legendRAAcomboPi0andEta0010->SetNColumns(2);
		legendRAAcomboPi0andEta0010->SetTextSize(FontSize);
		legendRAAcomboPi0andEta0010->SetHeader(Form("%s",collisionSystemPbPb0010.Data()));
		legendRAAcomboPi0andEta0010->AddEntry(graphCombRAASysPbPb2760GeVLHC11hPi0_0010,"#pi^{0}","pf");
		legendRAAcomboPi0andEta0010->AddEntry(graphCombRAASysPbPb2760GeVLHC11hEta_0010,"#eta","fp");
		legendRAAcomboPi0andEta0010->AddEntry(graphChargedKaonRAA0010,"K^{#pm}","fp");
		legendRAAcomboPi0andEta0010->AddEntry((TObject*)0,"","");
		legendRAAcomboPi0andEta0010->AddEntry((TObject*)0,"PLB 736 (2014)","");
		legendRAAcomboPi0andEta0010->Draw();

		boxErrorNorm0010Only->Draw();
	canvasRAAcomboPi0andEta->SaveAs(Form("%s/RAA_combinedPi0andEta_0010.%s",outputDir.Data(),suffix.Data()));
	canvasRAAcomboPi0andEta->SaveAs(Form("%s/RAA_combinedPi0andEta_0010.%s",paperPlots.Data(),suffix.Data()));

	
	canvasRAAcomboPi0andEta->cd();
	histo2DRAAcomboPi0andEta->Draw("copy");
	
		DrawGammaSetMarkerTGraphAsym(graphCombRAASysPbPb2760GeVLHC11hPi0_2050, 25, markerSizeComb, colorEMCal2050, colorEMCal2050, widthLinesBoxes, kTRUE);
		graphCombRAASysPbPb2760GeVLHC11hPi0_2050->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphCombRAAStatPbPb2760GeVLHC11hPi0_2050, 25, markerSizeComb, colorEMCal2050, colorEMCal2050);
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPb2760GeVLHC11hPi0_2050);
		graphCombRAAStatPbPb2760GeVLHC11hPi0_2050->Draw("p,same,e0");

		DrawGammaSetMarkerTGraphAsym(graphCombRAASysPbPb2760GeVLHC11hEta_2050, markerStyle2050, markerSizeComb, colorEMCal2050, colorEMCal2050, widthLinesBoxes, kTRUE);
		graphCombRAASysPbPb2760GeVLHC11hEta_2050->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphCombRAAStatPbPb2760GeVLHC11hEta_2050, markerStyle2050, markerSizeComb, colorEMCal2050 , colorEMCal2050);
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPb2760GeVLHC11hEta_2050);
		graphCombRAAStatPbPb2760GeVLHC11hEta_2050->Draw("p,same,e0");
		
		
		TLegend* legendRAAcomboPi0andEta2050 = new TLegend(0.33,0.8,0.75,0.9);
		legendRAAcomboPi0andEta2050->SetFillColor(0);
		legendRAAcomboPi0andEta2050->SetLineColor(0);
		legendRAAcomboPi0andEta2050->SetNColumns(2);
		legendRAAcomboPi0andEta2050->SetTextFont(42);
		legendRAAcomboPi0andEta2050->SetTextSize(FontSize);
		legendRAAcomboPi0andEta2050->SetHeader(Form("%s",collisionSystemPbPb2050.Data()));
		legendRAAcomboPi0andEta2050->AddEntry(graphCombRAASysPbPb2760GeVLHC11hPi0_2050,"#pi^{0}","pf");
		legendRAAcomboPi0andEta2050->AddEntry(graphCombRAASysPbPb2760GeVLHC11hEta_2050,"#eta","fp");
		legendRAAcomboPi0andEta2050->Draw();

		boxErrorNorm2050Only->Draw();
	canvasRAAcomboPi0andEta->SaveAs(Form("%s/RAA_combinedPi0andEta_2050.%s",outputDir.Data(),suffix.Data()));
	canvasRAAcomboPi0andEta->SaveAs(Form("%s/RAA_combinedPi0andEta_2050.%s",paperPlots.Data(),suffix.Data()));

	
	
	TCanvas* canvasRAAcombo = new TCanvas("canvasRAAcombo","",200,10,1200,1100);  //200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasRAAcombo, 0.08, 0.02, 0.035, 0.09);
	
		TH2F * histo2DRAAcombo;
		histo2DRAAcombo = new TH2F("histo2DRAAcombo","histo2DRAAcombo",11000,0.23,70.,1000,-0.5,2.);
		SetStyleHistoTH2ForGraphs(histo2DRAAcombo, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}",0.035,0.04, 0.035,0.04, 1.,.9);
		histo2DRAAcombo->GetYaxis()->SetRangeUser(0.,1.4);
		histo2DRAAcombo->GetXaxis()->SetRangeUser(0.,21);
		histo2DRAAcombo->Draw("copy");
			
		DrawGammaSetMarkerTGraphAsym(graphCombRAASysPbPb2760GeVLHC11h_0010, markerStyle0010, markerSizeComb, colorCombo0010 , colorCombo0010, 2, kTRUE);
		graphCombRAASysPbPb2760GeVLHC11h_0010->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphCombRAAStatPbPb2760GeVLHC11h_0010, markerStyle0010, markerSizeComb, colorCombo0010 , colorCombo0010);
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPb2760GeVLHC11h_0010);
		graphCombRAAStatPbPb2760GeVLHC11h_0010->Draw("p,same,e0");

		DrawGammaSetMarkerTGraphAsym(graphCombRAASysPbPb2760GeVLHC11h_2050, markerStyle2050, markerSizeComb, colorCombo2050, colorCombo2050, 2, kTRUE);
		graphCombRAASysPbPb2760GeVLHC11h_2050->Draw("E2same");
		DrawGammaSetMarkerTGraphAsym(graphCombRAAStatPbPb2760GeVLHC11h_2050, markerStyle2050, markerSizeComb, colorCombo2050 , colorCombo2050);
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPb2760GeVLHC11h_2050);
		graphCombRAAStatPbPb2760GeVLHC11h_2050->Draw("p,same,e0");

		labelSystRaa->Draw();
		
		TLegend* legendRAAcombo2 = new TLegend(0.595,0.73,0.91,0.88);
		legendRAAcombo2->SetFillColor(0);
		legendRAAcombo2->SetLineColor(0);
		legendRAAcombo2->SetTextFont(42);
		legendRAAcombo2->SetMargin(0.17);
		legendRAAcombo2->SetTextSize(FontSize);
		legendRAAcombo2->SetHeader(collisionSystem2760GeV.Data());
		legendRAAcombo2->AddEntry(graphCombRAASysPbPb2760GeVLHC11h_0010,Form("%s ",cent0010.Data()),"fp");
		legendRAAcombo2->AddEntry(graphCombRAASysPbPb2760GeVLHC11h_2050,Form("%s ",cent2050.Data()),"fp");
		legendRAAcombo2->Draw();

		boxErrorNorm0010->Draw();		
		boxErrorNorm2050->Draw();

	canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combined.%s",outputDir.Data(),meson.Data(),suffix.Data()));
	canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combined.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
	
	
	canvasRAAcombo->cd();
	histo2DRAAcombo->DrawCopy();

		labelSystRaa->Draw();
		
		legendRAAcombo2->Draw();

		TLegend* legendRAATheoryPbPb = new TLegend(0.595,0.54,0.91,0.68);
		legendRAATheoryPbPb->SetFillColor(0);
		legendRAATheoryPbPb->SetLineColor(0);
		legendRAATheoryPbPb->SetTextFont(42);
		legendRAATheoryPbPb->SetTextSize(FontSize);
		legendRAATheoryPbPb->SetMargin(0.17);
		
		if(meson.CompareTo("Pi0")==0){ 
			
			DrawGammaSetMarkerTGraphAsym(graphPi0RAAJetQuenching_0010, 0, 0, colorNLO0010, colorNLO0010, widthLinesBoxes, kTRUE, colorNLO0010);
			graphPi0RAAJetQuenching_0010->Draw("3,same");

			DrawGammaSetMarkerTGraphAsym(gWHDG_Raa_0010, markerStyleCommmonSpectrum0010,markerSizeSpectrum, colorWHDG0005, colorWHDG0005,widthLinesBoxes, kTRUE);
			gWHDG_Raa_0010->SetFillStyle(fillStyleWHDG);
			gWHDG_Raa_0010->SetFillColor(colorWHDG0005);
			gWHDG_Raa_0010->Draw("3 same");

			DrawGammaSetMarkerTGraphAsym(gWHDG_Raa_2050, markerStyleCommmonSpectrum2040,markerSizeSpectrum, colorWHDG4060, colorWHDG4060,widthLinesBoxes, kTRUE);
			gWHDG_Raa_2050->SetFillStyle(fillStyleWHDG);
			gWHDG_Raa_2050->SetFillColor(colorWHDG4060);
			gWHDG_Raa_2050->Draw("3 same");
			
			legendRAATheoryPbPb->AddEntry(graphPi0RAAJetQuenching_0010,"NLO DCZW - 0-10%","f");
			legendRAATheoryPbPb->AddEntry(gWHDG_Raa_0010,"WHDG - 0-10%","f");
			legendRAATheoryPbPb->AddEntry(gWHDG_Raa_2050,"WHDG - 20-40%","f");
		
		} else if(meson.CompareTo("Eta")==0){

			DrawGammaSetMarkerTGraphAsym(graphEtaRAAJetQuenching_0010, 0, 0, colorNLO0010, colorNLO0010, widthLinesBoxes, kTRUE, colorNLO0010);
			graphEtaRAAJetQuenching_0010->Draw("3,same");

			DrawGammaSetMarkerTGraphAsym(gWHDG_Eta_Raa_0010, markerStyleCommmonSpectrum0010,markerSizeSpectrum, colorWHDG0005, colorWHDG0005,widthLinesBoxes, kTRUE);
			gWHDG_Eta_Raa_0010->SetFillStyle(fillStyleWHDG);
			gWHDG_Eta_Raa_0010->SetFillColor(colorWHDG0005);
			gWHDG_Eta_Raa_0010->Draw("3 same");

			DrawGammaSetMarkerTGraphAsym(gWHDG_Eta_Raa_2050, markerStyleCommmonSpectrum2040,markerSizeSpectrum, colorWHDG4060, colorWHDG4060,widthLinesBoxes, kTRUE);
			gWHDG_Eta_Raa_2050->SetFillStyle(fillStyleWHDG);
			gWHDG_Eta_Raa_2050->SetFillColor(colorWHDG4060);
			gWHDG_Eta_Raa_2050->Draw("3 same");

			legendRAATheoryPbPb->AddEntry(graphEtaRAAJetQuenching_0010,"NLO DCZW - 0-10%","f");
			legendRAATheoryPbPb->AddEntry(gWHDG_Eta_Raa_0010,"WHDG - 0-10%","f");
			legendRAATheoryPbPb->AddEntry(gWHDG_Eta_Raa_2050,"WHDG - 20-50%","f");

		}

		legendRAATheoryPbPb->Draw();
		boxErrorNorm0010->Draw();		
		boxErrorNorm2050->Draw();

		graphCombRAASysPbPb2760GeVLHC11h_0010->Draw("E2same");
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPb2760GeVLHC11h_0010);
		graphCombRAAStatPbPb2760GeVLHC11h_0010->Draw("p,same,e0");
		graphCombRAASysPbPb2760GeVLHC11h_2050->Draw("E2same");
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPb2760GeVLHC11h_2050);
		graphCombRAAStatPbPb2760GeVLHC11h_2050->Draw("p,same,e0");

	canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedWithTheoryModels.%s",outputDir.Data(),meson.Data(),suffix.Data()));
	canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedWithTheoryModels.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
	
    
    
    
    
    Double_t arrayBoundariesX1_4[3];
    Double_t arrayBoundariesY1_4[3];
    Double_t relativeMarginsX[3];
    Double_t relativeMarginsY[3];
    ReturnCorrectValuesForCanvasScaling(2400,2000, 2, 2,0.05, 0.005, 0.005,0.06,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);

    TCanvas* canvasRaacomparisonModels = new TCanvas("canvasRaacomparisonModels","",0,0,2400,2000);  // gives the page size
//     DrawGammaCanvasSettings( canvasRaacomparisonModels,  0.13, 0.02, 0.03, 0.06);

    TPad* padRaacomparisonModels1 = new TPad("padRaacomparisonModels1", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
    DrawGammaPadSettings( padRaacomparisonModels1, relativeMarginsX[0], relativeMarginsX[1], relativeMarginsY[1], relativeMarginsY[2]);
    padRaacomparisonModels1->Draw();

    TPad* padRaacomparisonModels2 = new TPad("padRaacomparisonModels2", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padRaacomparisonModels2, relativeMarginsX[0], relativeMarginsX[1], relativeMarginsY[0], relativeMarginsY[1]);
    padRaacomparisonModels2->Draw();
    
    TPad* padRaacomparisonModels3 = new TPad("padRaacomparisonModels3", "", arrayBoundariesX1_4[1], arrayBoundariesY1_4[2], arrayBoundariesX1_4[2], arrayBoundariesY1_4[1],-1, -1, -2);
    DrawGammaPadSettings( padRaacomparisonModels3, relativeMarginsX[1], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[2]);
    padRaacomparisonModels3->Draw();

    TPad* padRaacomparisonModels4 = new TPad("padRaacomparisonModels4", "", arrayBoundariesX1_4[1], arrayBoundariesY1_4[1], arrayBoundariesX1_4[2], arrayBoundariesY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padRaacomparisonModels4, relativeMarginsX[1], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[1]);
    padRaacomparisonModels4->Draw();

    Int_t textSizeLabelsPixelRatio = 50;
    Double_t marginRatio = 0.16*2400;
    Double_t textsizeLabelsRatioUp = 0;
    Double_t textsizeFacRatioUp = 0;
    Double_t textsizeLabelsRatioDown = 0;
    Double_t textsizeFacRatioDown = 0;

    if (padRaacomparisonModels2->XtoPixel(padRaacomparisonModels2->GetX2()) < padRaacomparisonModels2->YtoPixel(padRaacomparisonModels2->GetY1())){
        textsizeLabelsRatioUp = (Double_t)textSizeLabelsPixelRatio/padRaacomparisonModels2->XtoPixel(padRaacomparisonModels2->GetX2()) ;
        textsizeFacRatioUp = (Double_t)1./padRaacomparisonModels2->XtoPixel(padRaacomparisonModels2->GetX2()) ;
    } else {
        textsizeLabelsRatioUp = (Double_t)textSizeLabelsPixelRatio/padRaacomparisonModels2->YtoPixel(padRaacomparisonModels2->GetY1());
        textsizeFacRatioUp = (Double_t)1./padRaacomparisonModels2->YtoPixel(padRaacomparisonModels2->GetY1());
    }
    if (padRaacomparisonModels1->XtoPixel(padRaacomparisonModels1->GetX2()) < padRaacomparisonModels1->YtoPixel(padRaacomparisonModels1->GetY1())){
        textsizeLabelsRatioDown = (Double_t)textSizeLabelsPixelRatio/padRaacomparisonModels1->XtoPixel(padRaacomparisonModels1->GetX2()) ;
        textsizeFacRatioDown = (Double_t)1./padRaacomparisonModels1->XtoPixel(padRaacomparisonModels1->GetX2()) ;
    } else {
        textsizeLabelsRatioDown = (Double_t)textSizeLabelsPixelRatio/padRaacomparisonModels1->YtoPixel(padRaacomparisonModels1->GetY1());
        textsizeFacRatioDown = (Double_t)1./padRaacomparisonModels1->YtoPixel(padRaacomparisonModels1->GetY1());
    }

    TH2F * histo2DRAAAll3Up;
    histo2DRAAAll3Up = new TH2F("histo2DRAAAll3Up","histo2DRAAAll3Up",1000,0.01,21.,1000,-0.05,10.);
    histo2DRAAAll3Up->GetYaxis()->SetRangeUser(0.01,1.4);
    SetStyleHistoTH2ForGraphs(histo2DRAAAll3Up, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}", 0.85*textsizeLabelsRatioUp,textsizeLabelsRatioUp, 0.85*textsizeLabelsRatioUp,textsizeLabelsRatioUp, 1,0.3/(textsizeFacRatioUp*marginRatio), 512, 505); 
    histo2DRAAAll3Up->GetXaxis()->SetLabelFont(42);
    histo2DRAAAll3Up->GetYaxis()->SetLabelFont(42);
    
    
    TH2F * histo2DRAAAll3Down;
    histo2DRAAAll3Down = new TH2F("histo2DRAAAll3Down","histo2DRAAAll3Down",1000,0.01,21.,1000,-0.05,10.);
    histo2DRAAAll3Down->GetYaxis()->SetRangeUser(0.01,1.4);
    SetStyleHistoTH2ForGraphs(histo2DRAAAll3Down, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}", 0.85*textsizeLabelsRatioDown, textsizeLabelsRatioDown,  0.85*textsizeLabelsRatioDown, textsizeLabelsRatioDown, 1,0.3/(textsizeFacRatioDown*marginRatio), 512, 505); 
    histo2DRAAAll3Down->GetXaxis()->SetLabelFont(42);
    histo2DRAAAll3Down->GetYaxis()->SetLabelFont(42);

    
//     TH2F * histo2DRAAwithModels = new TH2F("histo2DRAAwithModels","histo2DRAAwithModels",11000,0.,21.,1000,0.01,1.41);
//     SetStyleHistoTH2ForGraphs(histo2DRAAwithModels, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA} ",0.035,0.04, 0.035,0.04, 1.,.9);

    padRaacomparisonModels1->cd();
    histo2DRAAAll3Down->DrawCopy();
 
          DrawGammaSetMarkerTGraphAsym(graphEtaRAAJetQuenching_0010, 0, 0, colorNLO0010, colorNLO0010, widthLinesBoxes, kTRUE, colorNLO0010);
      graphEtaRAAJetQuenching_0010->Draw("3,same");

      DrawGammaSetMarkerTGraphAsym(gWHDG_Eta_Raa_0010, markerStyleCommmonSpectrum0010,markerSizeSpectrum, colorWHDG0005, colorWHDG0005,widthLinesBoxes, kTRUE);
      gWHDG_Eta_Raa_0010->SetFillStyle(fillStyleWHDG);
      gWHDG_Eta_Raa_0010->SetFillColor(colorWHDG0005);
      gWHDG_Eta_Raa_0010->Draw("3 same");

      //=============Combined Eta
      DrawGammaSetMarkerTGraphAsym(graphCombRAASysPbPb2760GeVLHC11hEta_0010, markerStyle0010, markerSizeComb, colorCombo0010, colorCombo0010, widthLinesBoxes, kTRUE);
      graphCombRAASysPbPb2760GeVLHC11hEta_0010->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphCombRAAStatPbPb2760GeVLHC11hEta_0010, markerStyle0010, markerSizeComb, colorCombo0010 , colorCombo0010);
      if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPb2760GeVLHC11hEta_0010);
      graphCombRAAStatPbPb2760GeVLHC11hEta_0010->Draw("p,same,e0");

      TLegend* legendRAAmodels2 = new TLegend(0.18,0.82,0.55,0.95);
      legendRAAmodels2->SetFillColor(0);
      legendRAAmodels2->SetLineColor(0);
      legendRAAmodels2->SetTextFont(42);
      legendRAAmodels2->SetMargin(0.17);
      legendRAAmodels2->SetTextSize(0.85*textsizeLabelsRatioDown);
      legendRAAmodels2->SetHeader(collisionSystemPbPb0010.Data());
      legendRAAmodels2->AddEntry(graphCombRAASysPbPb2760GeVLHC11h_0010,"#eta ALICE","fp");
      legendRAAmodels2->Draw();
  
	    boxErrorNorm0010_Single->Draw();

      DrawGammaLines(0., 20.5 , 1, 1 ,1,kGray, 2);

    histo2DRAAAll3Down->Draw("axis,same");
    padRaacomparisonModels1->Update();
    padRaacomparisonModels2->cd();
    histo2DRAAAll3Up->DrawCopy();

	  DrawGammaSetMarkerTGraphAsym(graphPi0RAAJetQuenching_0010, 0, 0, colorNLO0010, colorNLO0010, widthLinesBoxes, kTRUE, colorNLO0010);
      graphPi0RAAJetQuenching_0010->Draw("3,same");

      DrawGammaSetMarkerTGraphAsym(gWHDG_Raa_0010, markerStyleCommmonSpectrum0010,markerSizeSpectrum, colorWHDG0005, colorWHDG0005,widthLinesBoxes, kTRUE);
      gWHDG_Raa_0010->SetFillStyle(fillStyleWHDG);
      gWHDG_Raa_0010->SetFillColor(colorWHDG0005);
      gWHDG_Raa_0010->Draw("3 same");

      TGraphAsymmErrors *copyWHDG = (TGraphAsymmErrors*) gWHDG_Raa_0010->Clone();
      DrawGammaSetMarkerTGraphAsym(copyWHDG, markerStyleCommmonSpectrum0010,markerSizeSpectrum, kBlack, kBlack,widthLinesBoxes, kTRUE);
      copyWHDG->SetFillStyle(fillStyleWHDG);
      copyWHDG->SetFillColor(kBlack);
      
      //=============Combined Pi0
      DrawGammaSetMarkerTGraphAsym(graphCombRAASysPbPb2760GeVLHC11hPi0_0010, markerStyle0010, markerSizeComb, colorCombo0010, colorCombo0010, widthLinesBoxes, kTRUE);
      graphCombRAASysPbPb2760GeVLHC11hPi0_0010->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphCombRAAStatPbPb2760GeVLHC11hPi0_0010, markerStyle0010, markerSizeComb, colorCombo0010, colorCombo0010);
      if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPb2760GeVLHC11hPi0_0010);
      graphCombRAAStatPbPb2760GeVLHC11hPi0_0010->Draw("p,same,e0");

      TLegend* legendRAAmodels1 = new TLegend(0.19,0.73,0.89,0.95);
      legendRAAmodels1->SetFillColor(0);
      legendRAAmodels1->SetLineColor(0);
      legendRAAmodels1->SetTextFont(42);
      legendRAAmodels1->SetMargin(0.17);
      legendRAAmodels1->SetNColumns(2);
      legendRAAmodels1->SetTextSize(0.85*textsizeLabelsRatioUp);
      legendRAAmodels1->SetHeader(collisionSystemPbPb0010.Data());
      legendRAAmodels1->AddEntry(graphCombRAASysPbPb2760GeVLHC11h_0010,"#pi^{0} ALICE","fp");
      legendRAAmodels1->AddEntry(graphPi0RAAJetQuenching_0010,"NLO DCZW - 0-10%","f");
      legendRAAmodels1->AddEntry((TObject*)0,"","");
      legendRAAmodels1->AddEntry(copyWHDG,"WHDG","f");
      legendRAAmodels1->Draw();

	    boxErrorNorm0010_Single->Draw();

      DrawGammaLines(0., 20.5 , 1, 1 ,1,kGray, 2);
    
    histo2DRAAAll3Up->Draw("axis,same");
    padRaacomparisonModels2->Update();
    padRaacomparisonModels3->cd();
    histo2DRAAAll3Down->DrawCopy();
            
	DrawGammaSetMarkerTGraphAsym(gWHDG_Eta_Raa_2050, markerStyleCommmonSpectrum2040,markerSizeSpectrum, colorWHDG4060, colorWHDG4060,widthLinesBoxes, kTRUE);
      gWHDG_Eta_Raa_2050->SetFillStyle(fillStyleWHDG);
      gWHDG_Eta_Raa_2050->SetFillColor(colorWHDG4060);
      gWHDG_Eta_Raa_2050->Draw("3 same");

      DrawGammaSetMarkerTGraphAsym(graphCombRAASysPbPb2760GeVLHC11hEta_2050, markerStyle2050, markerSizeComb, colorCombo2050, colorCombo2050, widthLinesBoxes, kTRUE);
      graphCombRAASysPbPb2760GeVLHC11hEta_2050->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphCombRAAStatPbPb2760GeVLHC11hEta_2050, markerStyle2050, markerSizeComb, colorCombo2050 , colorCombo2050);
      if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPb2760GeVLHC11hEta_2050);
      graphCombRAAStatPbPb2760GeVLHC11hEta_2050->Draw("p,same,e0");

      TLegend* legendRAAmodels4 = new TLegend(0.1,0.82,0.47,0.95);
      legendRAAmodels4->SetFillColor(0);
      legendRAAmodels4->SetLineColor(0);
      legendRAAmodels4->SetTextFont(42);
      legendRAAmodels4->SetMargin(0.17);
      legendRAAmodels4->SetTextSize(0.85*textsizeLabelsRatioDown);
      legendRAAmodels4->SetHeader(collisionSystemPbPb2050.Data());
      legendRAAmodels4->AddEntry(graphCombRAASysPbPb2760GeVLHC11hEta_2050,"#eta ALICE","fp");
      legendRAAmodels4->Draw();

    boxErrorNorm2050_Single->Draw();

      DrawGammaLines(0., 20.5 , 1, 1 ,1,kGray, 2);

    histo2DRAAAll3Down->Draw("axis,same");
    padRaacomparisonModels3->Update();
    padRaacomparisonModels4->cd();
    histo2DRAAAll3Up->DrawCopy();
    

	  DrawGammaSetMarkerTGraphAsym(gWHDG_Raa_2050, markerStyleCommmonSpectrum2040,markerSizeSpectrum, colorWHDG4060, colorWHDG4060,widthLinesBoxes, kTRUE);
      gWHDG_Raa_2050->SetFillStyle(fillStyleWHDG);
      gWHDG_Raa_2050->SetFillColor(colorWHDG4060);
      gWHDG_Raa_2050->Draw("3 same");
      
      //=============Combined Pi0
      DrawGammaSetMarkerTGraphAsym(graphCombRAASysPbPb2760GeVLHC11hPi0_2050, markerStyle2050, markerSizeComb, colorCombo2050, colorCombo2050, widthLinesBoxes, kTRUE);
      graphCombRAASysPbPb2760GeVLHC11hPi0_2050->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphCombRAAStatPbPb2760GeVLHC11hPi0_2050, markerStyle2050, markerSizeComb, colorCombo2050, colorCombo2050);
      if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPb2760GeVLHC11hPi0_2050);
      graphCombRAAStatPbPb2760GeVLHC11hPi0_2050->Draw("p,same,e0");
        
      
      TLegend* legendRAAmodels3 = new TLegend(0.1,0.82,0.47,0.95);
      legendRAAmodels3->SetFillColor(0);
      legendRAAmodels3->SetLineColor(0);
      legendRAAmodels3->SetTextFont(42);
      legendRAAmodels3->SetMargin(0.17);
      legendRAAmodels3->SetTextSize(0.85*textsizeLabelsRatioUp);
      legendRAAmodels3->SetHeader(collisionSystemPbPb2050.Data());
      legendRAAmodels3->AddEntry(graphCombRAASysPbPb2760GeVLHC11hPi0_2050,"#pi^{0} ALICE","fp");
      legendRAAmodels3->Draw();

      boxErrorNorm2050_Single->Draw();
      DrawGammaLines(0., 19.5 , 1, 1 ,1, kGray, 2);   
        
    histo2DRAAAll3Up->Draw("axis,same");
    padRaacomparisonModels4->Update();

    canvasRaacomparisonModels->Update();    
    canvasRaacomparisonModels->SaveAs(Form("%s/RAA_TheoryModels.%s",outputDir.Data(),suffix.Data()));
    canvasRaacomparisonModels->SaveAs(Form("%s/RAA_TheoryModels.%s",paperPlots.Data(),suffix.Data()));
    delete padRaacomparisonModels1; 
    delete padRaacomparisonModels3; 
    delete padRaacomparisonModels2; 
    delete padRaacomparisonModels4; 
    delete canvasRaacomparisonModels;

    
  	canvasRAAcombo->cd();
	histo2DRAAcombo->DrawCopy();
	
		labelSystRaa->Draw();

		TLegend* legendRAAcomboCharged = new TLegend(0.595,0.7,0.91,0.88);
		legendRAAcomboCharged->SetFillColor(0);
		legendRAAcomboCharged->SetLineColor(0);
		legendRAAcomboCharged->SetTextFont(42);
		legendRAAcomboCharged->SetMargin(0.17);
		legendRAAcomboCharged->SetTextSize(FontSize);
		legendRAAcomboCharged->SetHeader(collisionSystem2760GeV.Data());
		if(meson.CompareTo("Pi0")==0){
			legendRAAcomboCharged->AddEntry(graphCombRAASysPbPb2760GeVLHC11h_0010,Form("#pi^{0} - %s ",cent0010.Data()),"fp");
		} else if(meson.CompareTo("Eta")==0){
			legendRAAcomboCharged->AddEntry(graphCombRAASysPbPb2760GeVLHC11h_0010,Form("#eta - %s ",cent0010.Data()),"fp");

		}
	
		if(meson.CompareTo("Pi0")==0){ 
			
			DrawGammaSetMarkerTGraphAsym(graphChargedPionRAA0010, 21,2, colorCharged,colorCharged, 0.1, kFALSE);
			if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphChargedPionRAA0010);
			DrawGammaSetMarkerTGraphAsym(graphChargedPionRAASys0010,21,2, colorCharged,colorCharged, 1, kTRUE);
			graphChargedPionRAASys0010->Draw("2same");
			graphChargedPionRAA0010->Draw("p,same"); 
			legendRAAcomboCharged->AddEntry(graphChargedPionRAA0010,Form("#pi^{#pm} - %s ",cent0010.Data()),"fp");
			legendRAAcomboCharged->AddEntry((TObject*)0,"PLB 736 (2014)","");
			
		} else if(meson.CompareTo("Eta")==0){

			DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAA0010, 21,2, colorCharged,colorCharged, 0.1, kFALSE);
			if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphChargedKaonRAA0010);
			DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAASys0010,21,2, colorCharged,colorCharged, 1, kTRUE);
			graphChargedKaonRAASys0010->Draw("2same");
			graphChargedKaonRAA0010->Draw("p,same"); 
			legendRAAcomboCharged->AddEntry(graphChargedKaonRAA0010,Form("K^{#pm} - %s ",cent0010.Data()),"fp");
			legendRAAcomboCharged->AddEntry((TObject*)0,"PLB 736 (2014)","");
			
		}
		
		graphCombRAASysPbPb2760GeVLHC11h_0010->Draw("E2same");
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPb2760GeVLHC11h_0010);
		graphCombRAAStatPbPb2760GeVLHC11h_0010->Draw("p,same,e0");

		boxErrorNorm0010Only->Draw();		

		legendRAAcomboCharged->Draw();
	
	canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedWithCharged_0010.%s",outputDir.Data(),meson.Data(),suffix.Data()));
	canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedWithCharged_0010.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
	

	canvasRAAcombo->cd();
	histo2DRAAcombo->DrawCopy();
		
		labelSystRaa->Draw();
		
		TLegend* legendRAAcomboCharged2050 = new TLegend(0.595,0.7,0.91,0.88);
		legendRAAcomboCharged2050->SetFillColor(0);
		legendRAAcomboCharged2050->SetLineColor(0);
		legendRAAcomboCharged2050->SetTextFont(42);
		legendRAAcomboCharged2050->SetMargin(0.17);
		legendRAAcomboCharged2050->SetTextSize(FontSize);
		legendRAAcomboCharged2050->SetHeader(collisionSystem2760GeV.Data());
		if(meson.CompareTo("Pi0")==0){
			legendRAAcomboCharged2050->AddEntry(graphCombRAASysPbPb2760GeVLHC11h_2050,Form("#pi^{0} - %s ",cent2050.Data()),"fp");
		} else if(meson.CompareTo("Eta")==0){
			legendRAAcomboCharged2050->AddEntry(graphCombRAASysPbPb2760GeVLHC11h_2050,Form("#eta - %s ",cent2050.Data()),"fp");

		}
	
		if(meson.CompareTo("Pi0")==0){ 
			
			DrawGammaSetMarkerTGraphAsym(graphChargedPionRAA2040, 21,2, colorCharged,colorCharged, 0.1, kFALSE);
			if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphChargedPionRAA2040);
			DrawGammaSetMarkerTGraphAsym(graphChargedPionRAASys2040,21,2, colorCharged,colorCharged, 1, kTRUE);
			graphChargedPionRAASys2040->Draw("2same");
			graphChargedPionRAA2040->Draw("p,same"); 
			legendRAAcomboCharged2050->AddEntry(graphChargedPionRAA2040,"#pi^{#pm} - 20-40%","fp");
			legendRAAcomboCharged2050->AddEntry((TObject*)0,"PLB 736 (2014)","");
			
		} else if(meson.CompareTo("Eta")==0){

			DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAA2040, 21,2,  colorCharged,colorCharged, 0.1, kFALSE);
			if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphChargedKaonRAA2040);
			DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAASys2040,21,2, colorCharged,colorCharged, 1, kTRUE);
			graphChargedKaonRAASys2040->Draw("2same");
			graphChargedKaonRAA2040->Draw("p,same"); 
			legendRAAcomboCharged2050->AddEntry(graphChargedKaonRAA2040,"K^{#pm} - 20-40%","fp");
			legendRAAcomboCharged2050->AddEntry((TObject*)0,"PLB 736 (2014)","");			
			
		}
		
		graphCombRAASysPbPb2760GeVLHC11h_2050->Draw("E2same");
		if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPb2760GeVLHC11h_2050);
		graphCombRAAStatPbPb2760GeVLHC11h_2050->Draw("p,same,e0");

		boxErrorNorm2050Only->Draw();

		legendRAAcomboCharged2050->Draw();
	
	canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedWithCharged_2040.%s",outputDir.Data(),meson.Data(),suffix.Data()));
	canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedWithCharged_2040.%s",paperPlots.Data(),meson.Data(),suffix.Data()));


	Int_t textSizeLabelsPixelRAA = 50;
	Double_t marginRAA = 0.14*1200;
	Double_t textsizeLabelsRAA = 0;
	Double_t textsizeFacRAA = 0;

	if (canvasRAAcombo->XtoPixel(canvasRAAcombo->GetX2()) < canvasRAAcombo->YtoPixel(canvasRAAcombo->GetY1())){
		textsizeLabelsRAA = (Double_t)textSizeLabelsPixelRAA/canvasRAAcombo->XtoPixel(canvasRAAcombo->GetX2()) ;
		textsizeFacRAA = (Double_t)1./canvasRAAcombo->XtoPixel(canvasRAAcombo->GetX2()) ;
	} else {
		textsizeLabelsRAA = (Double_t)textSizeLabelsPixelRAA/canvasRAAcombo->YtoPixel(canvasRAAcombo->GetY1());
		textsizeFacRAA = (Double_t)1./canvasRAAcombo->YtoPixel(canvasRAAcombo->GetY1());
	}

	if(meson.CompareTo("Pi0")==0){
		canvasRAAcombo->cd();
			TH2F * histo2DRAAcomboPHENIX;
			histo2DRAAcomboPHENIX = new TH2F("histo2DRAAcomboPHENIX","histo2DRAAcomboPHENIX",11000,0.23,70.,1000,-0.5,2.);
			SetStyleHistoTH2ForGraphs(histo2DRAAcomboPHENIX, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}",0.035,0.04, 0.035,0.04, 1.,.92);
	// 		histo2DRAAcombo->GetXaxis()->SetMoreLogLabels();	
	// 		histo2DRAAcombo->GetXaxis()->SetLabelOffset(-0.01);
			histo2DRAAcomboPHENIX->GetYaxis()->SetRangeUser(0.,2.);
			histo2DRAAcomboPHENIX->GetXaxis()->SetRangeUser(0.,20.01);
			histo2DRAAcomboPHENIX->Draw("copy");
			
			DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVRAA_0010, markerStylePHENIX200GeV,markerSizePHENIX200GeV, kBlack , kBlack);
			DrawGammaSetMarkerTGraphErr(graphPHENIX62GeVRAA_0010, markerStylePHENIX62GeV,markerSizePHENIX62GeV, kBlack, kBlack);
			DrawGammaSetMarkerTGraphErr(graphPHENIX39GeVRAA_0010, markerStylePHENIX39GeV,markerSizePHENIX39GeV, kBlack , kBlack);
			DrawGammaSetMarkerTGraphErr(graphWA98_17_3GeVRAA_0013, markerStyleWA98,markerSizeWA98, kGray+2 , kGray+2);

			graphPHENIX200GeVRAA_0010->Draw("p,same,e0");	
			graphPHENIX39GeVRAA_0010->Draw("p,same,e0");	
			graphPHENIX62GeVRAA_0010->Draw("p,same,e0");	
			graphWA98_17_3GeVRAA_0013->Draw("p,same,e0");	

			graphCombRAASysPbPb2760GeVLHC11h_0010->Draw("E2same");
			if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPb2760GeVLHC11h_0010);
		graphCombRAAStatPbPb2760GeVLHC11h_0010->Draw("p,same,e0");

			TLatex *labelRAAALICEPbPb0010 = new TLatex(0.35,0.9,"#pi^{0} ALICE 0#font[122]{-}10% Pb#font[122]{-}Pb (2011)");
			SetStyleTLatex( labelRAAALICEPbPb0010, 0.85*textsizeLabelsRAA,4);
			labelRAAALICEPbPb0010->Draw();
			TLegend* legendRAASinglePbPb0010 = new TLegend(0.35,0.84,0.65,0.885);
			legendRAASinglePbPb0010->SetFillColor(0);
			legendRAASinglePbPb0010->SetLineColor(0);
			legendRAASinglePbPb0010->SetNColumns(1);
			legendRAASinglePbPb0010->SetTextFont(42);
			legendRAASinglePbPb0010->SetTextSize(0.85*textsizeLabelsRAA);
			legendRAASinglePbPb0010->SetMargin(0.17);
			legendRAASinglePbPb0010->AddEntry(graphCombRAASysPbPb2760GeVLHC11h_0010,"0#font[122]{-}10% #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","pf");
			legendRAASinglePbPb0010->Draw();

			TLatex *labelRAAPHENIXPbPb0010 = new TLatex(0.35,0.79,"#pi^{0} PHENIX 0#font[122]{-}10% Au#font[122]{-}Au");
			SetStyleTLatex( labelRAAPHENIXPbPb0010, 0.85*textsizeLabelsRAA,4);
			labelRAAPHENIXPbPb0010->Draw();

			TLegend* legendRAARHICPbPb0010 = new TLegend(0.35,0.66,0.95,0.78);
			legendRAARHICPbPb0010->SetFillColor(0);
			legendRAARHICPbPb0010->SetLineColor(0);
			legendRAARHICPbPb0010->SetNColumns(2);
			legendRAARHICPbPb0010->SetTextFont(42);
			legendRAARHICPbPb0010->SetMargin(0.17);
			legendRAARHICPbPb0010->SetTextSize(0.85*textsizeLabelsRAA);
			legendRAARHICPbPb0010->AddEntry(graphPHENIX200GeVRAA_0010,"#sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
			legendRAARHICPbPb0010->AddEntry(graphPHENIX62GeVRAA_0010,"#sqrt{#it{s}_{_{NN}}} = 62.4 GeV","p");
			legendRAARHICPbPb0010->AddEntry(graphPHENIX39GeVRAA_0010,"#sqrt{#it{s}_{_{NN}}} = 39 GeV","p");
			legendRAARHICPbPb0010->Draw();

			TLatex *labelRAAWA98PbPb0010 = new TLatex(0.35,0.61,"#pi^{0} WA98     0#font[122]{-}13% Pb#font[122]{-}Pb");
			SetStyleTLatex( labelRAAWA98PbPb0010, 0.85*textsizeLabelsRAA,4);
			labelRAAWA98PbPb0010->Draw();

			TLegend* legendRAASPSPbPb0010 = new TLegend(0.35,0.55,0.95,0.59);
			legendRAASPSPbPb0010->SetFillColor(0);
			legendRAASPSPbPb0010->SetLineColor(0);
			legendRAASPSPbPb0010->SetNColumns(2);
			legendRAASPSPbPb0010->SetTextFont(42);
			legendRAASPSPbPb0010->SetTextSize(0.85*textsizeLabelsRAA);
			legendRAASPSPbPb0010->SetMargin(0.17);
			legendRAASPSPbPb0010->AddEntry(graphWA98_17_3GeVRAA_0013,"#sqrt{#it{s}_{_{NN}}} = 17.3 GeV","p");
			legendRAASPSPbPb0010->Draw();

			boxErrorNorm0010_Single->Draw();
			DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
		
		
		canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedwithPHENIX_asPaper_0010.%s",outputDir.Data(),meson.Data(),suffix.Data()));
		canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedwithPHENIX_asPaper_0010.%s",paperPlots.Data(),meson.Data(),suffix.Data()));

		
		canvasRAAcombo->cd();
		histo2DRAAcomboPHENIX->DrawCopy();

			DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVRAA_2040, markerStylePHENIX200GeV,markerSizePHENIX200GeV, kBlack , kBlack);
			DrawGammaSetMarkerTGraphErr(graphPHENIX62GeVRAA_2040, markerStylePHENIX62GeV,markerSizePHENIX62GeV, kBlack, kBlack);
			DrawGammaSetMarkerTGraphErr(graphPHENIX39GeVRAA_2040, markerStylePHENIX39GeV,markerSizePHENIX39GeV, kBlack , kBlack);

			graphPHENIX200GeVRAA_2040->Draw("p,same,e0");	
			graphPHENIX39GeVRAA_2040->Draw("p,same,e0");	
			graphPHENIX62GeVRAA_2040->Draw("p,same,e0");	

			graphCombRAASysPbPb2760GeVLHC11h_2050->Draw("E2same");
			if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPb2760GeVLHC11h_2050);
			graphCombRAAStatPbPb2760GeVLHC11h_2050->Draw("p,same,e0");

			TLatex *labelRAAALICEPbPb2040 = new TLatex(0.5,0.87,"#pi^{0} ALICE Pb#font[122]{-}Pb (2011)");
			SetStyleTLatex( labelRAAALICEPbPb2040, 0.85*textsizeLabelsRAA,4);
			labelRAAALICEPbPb2040->Draw();
			TLegend* legendRAASinglePbPb2040 = new TLegend(0.5,0.81,0.83,0.85);
			legendRAASinglePbPb2040->SetFillColor(0);
			legendRAASinglePbPb2040->SetLineColor(0);
			legendRAASinglePbPb2040->SetNColumns(1);
			legendRAASinglePbPb2040->SetTextFont(42);
			legendRAASinglePbPb2040->SetTextSize(0.85*textsizeLabelsRAA);
			legendRAASinglePbPb2040->SetMargin(0.17);
			legendRAASinglePbPb2040->AddEntry(graphCombRAASysPbPb2760GeVLHC11h_2050,"20#font[122]{-}50% #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","pf");
			legendRAASinglePbPb2040->Draw();

			TLatex *labelRAAPHENIXPbPb2040 = new TLatex(0.5,0.75,"#pi^{0} PHENIX 20-40% Au#font[122]{-}Au");
			SetStyleTLatex( labelRAAPHENIXPbPb2040, 0.85*textsizeLabelsRAA,4);
			labelRAAPHENIXPbPb2040->Draw();

			TLegend* legendRAARHICPbPb2040 = new TLegend(0.5,0.55,0.84,0.73);
			legendRAARHICPbPb2040->SetFillColor(0);
			legendRAARHICPbPb2040->SetLineColor(0);
		// 	legendRAARHICPbPb2040->SetNColumns(2);
			legendRAARHICPbPb2040->SetTextFont(42);
			legendRAARHICPbPb2040->SetTextSize(0.85*textsizeLabelsRAA);
			legendRAARHICPbPb2040->SetMargin(0.17);
			legendRAARHICPbPb2040->AddEntry(graphPHENIX200GeVRAA_2040,"#sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
			legendRAARHICPbPb2040->AddEntry(graphPHENIX62GeVRAA_2040,"#sqrt{#it{s}_{_{NN}}} = 62.4 GeV","p");
			legendRAARHICPbPb2040->AddEntry(graphPHENIX39GeVRAA_2040,"#sqrt{#it{s}_{_{NN}}} = 39 GeV","p");
			legendRAARHICPbPb2040->Draw();

	// 		boxErrorNorm2040_Single->Draw();
			boxErrorNorm2050_Single->Draw();
			DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
			
		
		canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedwithPHENIX_asPaper_2050.%s",outputDir.Data(),meson.Data(),suffix.Data()));
		canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedwithPHENIX_asPaper_2050.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
		
	}

	
	if(meson.CompareTo("Eta")==0){ 	

		canvasRAAcombo->cd();
		histo2DRAAcombo->DrawCopy();
			
			graphCombRAASysPbPb2760GeVLHC11h_0010->Draw("E2same");
			if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPb2760GeVLHC11h_0010);
			graphCombRAAStatPbPb2760GeVLHC11h_0010->Draw("p,same,e0");

			DrawGammaSetMarkerTGraphAsym(graphEtaRAAPhenix0010, markerStylePHENIX200GeV,markerSizePHENIX200GeV, kBlack , kBlack);
			graphEtaRAAPhenix0010->Draw("p,same"); 
		
			TLegend* legendEtaRAAcomp =new TLegend(0.5,0.83,0.8,0.93);
			legendEtaRAAcomp->SetFillColor(0);
			legendEtaRAAcomp->SetLineColor(0);
			legendEtaRAAcomp->SetTextFont(42);
			legendEtaRAAcomp->SetMargin(0.17);
			legendEtaRAAcomp->SetTextSize(0.85*textsizeLabelsRAA);
			legendEtaRAAcomp->SetHeader("#eta ALICE 0#font[122]{-}10% Pb#font[122]{-}Pb (2011)");
			legendEtaRAAcomp->AddEntry(graphCombRAASysPbPb2760GeVLHC11h_0010,"0#font[122]{-}10% #sqrt{#it{s}_{_{NN}}} = 2.76 TeV", "pf");
			legendEtaRAAcomp->Draw();
			TLegend* legendEtaRAAcompPH =new TLegend(0.5,0.72,0.8,0.82);
			legendEtaRAAcompPH->SetFillColor(0);
			legendEtaRAAcompPH->SetLineColor(0);
			legendEtaRAAcompPH->SetTextFont(42);
			legendEtaRAAcompPH->SetMargin(0.17);
			legendEtaRAAcompPH->SetTextSize(0.85*textsizeLabelsRAA);
			legendEtaRAAcompPH->SetHeader("#eta PHENIX 0#font[122]{-}10% Au#font[122]{-}Au");
			legendEtaRAAcompPH->AddEntry(graphEtaRAAPhenix0010,"#sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
			legendEtaRAAcompPH->Draw();
			
			boxErrorNorm0010_Single->Draw();
			DrawGammaLines(0., 21 , 1, 1 ,1,kGray);
					
		canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedwithPHENIX_0010.%s",outputDir.Data(),meson.Data(),suffix.Data()));
		canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedwithPHENIX_0010.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
		
		canvasRAAcombo->cd();
		histo2DRAAcombo->Draw("copy");

			DrawGammaSetMarkerTGraphAsym(graphEtaRAAPhenix2040, 25,2, colorPhenix,colorPhenix, 0.1, kFALSE);
			DrawGammaSetMarkerTGraphAsym(graphEtaRAAPhenix2040,24,2, colorPhenix,colorPhenix, 1, kTRUE,colorPhenix);
			graphEtaRAAPhenix2040->SetFillStyle(3003);
			graphEtaRAAPhenix2040->Draw("2same");
			graphEtaRAAPhenix2040->Draw("p,same"); 

			DrawGammaSetMarkerTGraphAsym(graphEtaRAAPhenix2060, 25,2, colorPhenix,colorPhenix, 0.1, kFALSE);
			DrawGammaSetMarkerTGraphAsym(graphEtaRAAPhenix2060,25,2, colorPhenix,colorPhenix, 1, kTRUE, colorPhenix);
			graphEtaRAAPhenix2060->SetFillStyle(3002);
			graphEtaRAAPhenix2060->Draw("2same");
			graphEtaRAAPhenix2060->Draw("p,same"); 
		
			graphCombRAASysPbPb2760GeVLHC11h_2050->Draw("E2same");
			if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombRAAStatPbPb2760GeVLHC11h_2050);
			graphCombRAAStatPbPb2760GeVLHC11h_2050->Draw("p,same,e0");

			TLegend* legendEtaRAAcomp1 = new TLegend(0.48,0.73,0.8,0.93);
			legendEtaRAAcomp1->SetFillColor(0);
			legendEtaRAAcomp1->SetLineColor(0);
			legendEtaRAAcomp1->SetTextFont(42);
			legendEtaRAAcomp1->SetMargin(0.17);
			legendEtaRAAcomp1->SetTextSize(FontSize);
			legendEtaRAAcomp1->SetHeader("#eta - semicentral");
			legendEtaRAAcomp1->AddEntry(graphCombRAASysPbPb2760GeVLHC11h_2050,"20-50%, Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV", "pf");
			legendEtaRAAcomp1->AddEntry(graphEtaRAAPhenix2040,"20-40%, Au-Au #sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
			legendEtaRAAcomp1->AddEntry(graphEtaRAAPhenix2060,"20-60%, Au-Au #sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
			legendEtaRAAcomp1->Draw();
			
			boxErrorNorm2050Only->Draw();
		
		canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedwithPHENIX_2050.%s",outputDir.Data(),meson.Data(),suffix.Data()));
		canvasRAAcombo->SaveAs(Form("%s/%s_RAA_combinedwithPHENIX_2050.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
		
	}
	
			
	if(meson.CompareTo("Eta")==0){
		//*********************************************************************************************************************//
		//************************************	 Combination of EtatoPi0 	***************************************************//
		cout << "\n\n\n *************************	 Combination of Eta to Pi0 ratio 	************************* \n\n\n" << endl;
		
        if(bWCorrection.CompareTo("X")==0 ){
          
          
          
        }
        
		TLatex *labelEtaToPi0Energy = new TLatex(0.12,0.87,collisionSystem2760GeV.Data());
		SetStyleTLatex( labelEtaToPi0Energy,0.04,4);
// 		labelEtaToPi0Energy->SetTextFont(63);

		TLatex *labelPreliminaryEtaToPi0 = new TLatex(0.65,0.15,"ALICE Preliminary");
		SetStyleTLatex( labelPreliminaryEtaToPi0, 0.04,4);

		//===================== loading pp inputs =============================
		TFile* fEtatoPi0input = new TFile("EtaToPi0InputsForCombination.root");
		TH1D *histoPCMEtaToPi0RatioPbPb0010 = (TH1D*)fEtatoPi0input->Get("histoPCMEtaToPi0RatioPbPb0010");
		TH1D *histoPCMEtaToPi0RatioPbPb2050 = (TH1D*)fEtatoPi0input->Get("histoPCMEtaToPi0RatioPbPb2050");
		TH1D *histoPCMEtaToPi0RatiopPb = (TH1D*)fEtatoPi0input->Get("histoPCMEtaToPi0RatiopPb");
		TGraphAsymmErrors *graphPCMEtaToPi0RatioSysErrpPb = (TGraphAsymmErrors*)fEtatoPi0input->Get("graphPCMEtaToPi0RatioSysErrpPb");
			
		TGraphAsymmErrors *graphCombEtaToPi0RatioSysErrpp7TeV = (TGraphAsymmErrors*)fEtatoPi0input->Get("graphCombEtaToPi0RatioSysErrpp7TeV");
		DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RatioSysErrpp7TeV, markerStylepp, markerSizepp, kBlack, kBlack, 1, kTRUE);
		TGraphAsymmErrors *graphCombEtaToPi0Ratiopp7TeVNoXErrors = (TGraphAsymmErrors*)fEtatoPi0input->Get("graphCombEtaToPi0Ratiopp7TeVNoXErrors");
		DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Ratiopp7TeVNoXErrors, markerStylepp, markerSizepp, kBlack, kBlack, 1, kTRUE);

		TH1D *cocktailEtaToPi0Ratio7TeVRebined = (TH1D*)fEtatoPi0input->Get("cocktailEtaToPi0Ratio7TeVRebined");
		TH1D *cocktailEtaToPi0Ratio_MtScaledRebinned = (TH1D*)fEtatoPi0input->Get("cocktailEtaToPi0Ratio_MtScaledRebinned");
		TH1D *cocktailEtaToPi0Ratio_K0ScaledRebinned = (TH1D*)fEtatoPi0input->Get("cocktailEtaToPi0Ratio_K0ScaledRebinned");

		cocktailEtaToPi0Ratio7TeVRebined->GetXaxis()->SetRangeUser(0.1,16.);
		cocktailEtaToPi0Ratio_K0ScaledRebinned->GetXaxis()->SetRangeUser(0.1,16.);
		cocktailEtaToPi0Ratio_MtScaledRebinned->GetXaxis()->SetRangeUser(0.1,16.);
		cocktailEtaToPi0Ratio7TeVRebined->SetLineStyle(5);
		cocktailEtaToPi0Ratio7TeVRebined->SetLineWidth(2.5);
		cocktailEtaToPi0Ratio_K0ScaledRebinned->SetLineStyle(6);
		cocktailEtaToPi0Ratio_K0ScaledRebinned->SetLineWidth(2.5);
		cocktailEtaToPi0Ratio_MtScaledRebinned->SetLineStyle(7);
		cocktailEtaToPi0Ratio_MtScaledRebinned->SetLineWidth(2.5);
		//======================================================================
		
		TGraphAsymmErrors* graphCombEtatoPi0Stat2760GeVLHC11h_0010 = NULL;
		TGraphAsymmErrors* graphCombEtatoPi0Sys2760GeVLHC11h_0010 = NULL;
		TGraphAsymmErrors* graphCombEtatoPi0Tot2760GeVLHC11h_0010 = CombinePtPointsSpectraFullCorrMat( statErrorCollectionEtatoPi0LHC11h_0010,	sysErrorCollectionEtatoPi0LHC11h_0010, 	
																									xPtLimitsEta, /*17*/13, offSetsEta, offSetsEtaSys,
																									graphCombEtatoPi0Stat2760GeVLHC11h_0010, graphCombEtatoPi0Sys2760GeVLHC11h_0010, "weightEtatoPi0_0010.dat",1 );
		TGraphAsymmErrors* graphCombEtatoPi0Stat2760GeVLHC11h_2050 = NULL;
		TGraphAsymmErrors* graphCombEtatoPi0Sys2760GeVLHC11h_2050 = NULL;
		TGraphAsymmErrors* graphCombEtatoPi0Tot2760GeVLHC11h_2050 = CombinePtPointsSpectraFullCorrMat( statErrorCollectionEtatoPi0LHC11h_2050,	sysErrorCollectionEtatoPi0LHC11h_2050, 	
																									xPtLimitsEta, /*17*/13, offSetsEta, offSetsEtaSys,
																									graphCombEtatoPi0Stat2760GeVLHC11h_2050, graphCombEtatoPi0Sys2760GeVLHC11h_2050, "weightEtatoPi0_2050.dat",1 );

		graphCombEtatoPi0Tot2760GeVLHC11h_0010->RemovePoint(0);
		graphCombEtatoPi0Tot2760GeVLHC11h_0010->RemovePoint(0);
		graphCombEtatoPi0Tot2760GeVLHC11h_2050->RemovePoint(0);
		graphCombEtatoPi0Tot2760GeVLHC11h_2050->RemovePoint(0);
		
		graphCombEtatoPi0Stat2760GeVLHC11h_0010->RemovePoint(0);
		graphCombEtatoPi0Stat2760GeVLHC11h_0010->RemovePoint(0);
		
		graphCombEtatoPi0Sys2760GeVLHC11h_0010->RemovePoint(0);
		graphCombEtatoPi0Sys2760GeVLHC11h_0010->RemovePoint(0);

		graphCombEtatoPi0Stat2760GeVLHC11h_2050->RemovePoint(0);
		graphCombEtatoPi0Stat2760GeVLHC11h_2050->RemovePoint(0);

		graphCombEtatoPi0Sys2760GeVLHC11h_2050->RemovePoint(0);
		graphCombEtatoPi0Sys2760GeVLHC11h_2050->RemovePoint(0);


		TCanvas* canvasEtatoPi0combo = new TCanvas("canvasEtatoPi0combo","",200,10,1200,1100);  //200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasEtatoPi0combo, 0.09, 0.03, 0.035, 0.1);
		canvasEtatoPi0combo->SetLogx();
		
			TH2F * histo2DEtatoPi0combo = new TH2F("histo2DEtatoPi0combo","histo2DEtatoPi0combo",11000,minPtRange,maxPtRange,1000,0.01,1.2);
			SetStyleHistoTH2ForGraphs(histo2DEtatoPi0combo, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}",0.035,0.04, 0.035,0.04, 1.2,1.);
			histo2DEtatoPi0combo->GetXaxis()->SetMoreLogLabels();
			histo2DEtatoPi0combo->GetXaxis()->SetLabelOffset(-0.01);
            histo2DEtatoPi0combo->GetYaxis()->SetRangeUser(0.01,maxYEtatoPi0);
			histo2DEtatoPi0combo->Draw("copy");
            
			DrawGammaSetMarkerTGraphAsym(graphCombEtatoPi0Sys2760GeVLHC11h_0010, markerStyle0010, markerSizeComb, colorCombo0010 , colorCombo0010, 2, kTRUE);
			graphCombEtatoPi0Sys2760GeVLHC11h_0010->Draw("E2same");
			DrawGammaSetMarkerTGraphAsym(graphCombEtatoPi0Stat2760GeVLHC11h_0010, markerStyle0010, markerSizeComb, colorCombo0010 , colorCombo0010);
			if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombEtatoPi0Stat2760GeVLHC11h_0010);
			graphCombEtatoPi0Stat2760GeVLHC11h_0010->Draw("p,same,e0");

			DrawGammaSetMarkerTGraphAsym(graphCombEtatoPi0Sys2760GeVLHC11h_2050, markerStyle2050, markerSizeComb, colorCombo2050 , colorCombo2050, 2, kTRUE);
			graphCombEtatoPi0Sys2760GeVLHC11h_2050->Draw("E2same");
			DrawGammaSetMarkerTGraphAsym(graphCombEtatoPi0Stat2760GeVLHC11h_2050, markerStyle2050, markerSizeComb, colorCombo2050 , colorCombo2050);
			if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombEtatoPi0Stat2760GeVLHC11h_2050);
			graphCombEtatoPi0Stat2760GeVLHC11h_2050->Draw("p,same,e0");
			
			TLegend* legendEtatoPi0combo_onlyPbPb = new TLegend(0.12,0.76,0.53,0.92);
//             TLegend* legendEtatoPi0combo_onlyPbPb = GetAndSetLegend2(0.12, 0.76-3*0.06, 0.53, 0.76, 0.85* textSizeLabelsPixelEtaPi0Ratio);
			legendEtatoPi0combo_onlyPbPb->SetFillColor(0);
			legendEtatoPi0combo_onlyPbPb->SetLineColor(0);
			legendEtatoPi0combo_onlyPbPb->SetTextFont(42);
			legendEtatoPi0combo_onlyPbPb->SetTextSize(0.037);
			legendEtatoPi0combo_onlyPbPb->SetMargin(0.17);
            legendEtatoPi0combo_onlyPbPb->SetHeader(collisionSystem2760GeV.Data());
// 			legendEtatoPi0combo_onlyPbPb->AddEntry(graphCombEtatoPi0Sys2760GeVLHC11h_0010,Form("%s %s",cent0010.Data(),collisionSystem2760GeV.Data()),"fp");
// 			legendEtatoPi0combo_onlyPbPb->AddEntry(graphCombEtatoPi0Sys2760GeVLHC11h_2050,Form("%s %s",cent2050.Data(),collisionSystem2760GeV.Data()),"fp");
			legendEtatoPi0combo_onlyPbPb->AddEntry(graphCombEtatoPi0Sys2760GeVLHC11h_0010,Form("  %s",cent0010.Data()),"fp");
			legendEtatoPi0combo_onlyPbPb->AddEntry(graphCombEtatoPi0Sys2760GeVLHC11h_2050,Form("%s",cent2050.Data()),"fp");
			legendEtatoPi0combo_onlyPbPb->Draw();
			
		canvasEtatoPi0combo->SaveAs(Form("%s/%s_EtatoPi0Ratio_combined_DataOnly.%s",outputDir.Data(),meson.Data(),suffix.Data()));
		canvasEtatoPi0combo->SaveAs(Form("%s/%s_EtatoPi0Ratio_combined_DataOnly.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
		
		
		
		canvasEtatoPi0combo->cd();
			histo2DEtatoPi0combo->Draw("copy");
         
			graphCombEtatoPi0Sys2760GeVLHC11h_0010->Draw("E2same");
			if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombEtatoPi0Stat2760GeVLHC11h_0010);
			graphCombEtatoPi0Stat2760GeVLHC11h_0010->Draw("p,same,e0");

			graphCombEtatoPi0Sys2760GeVLHC11h_2050->Draw("E2same");
			if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombEtatoPi0Stat2760GeVLHC11h_2050);
			graphCombEtatoPi0Stat2760GeVLHC11h_2050->Draw("p,same,e0");

			graphCombEtaToPi0RatioSysErrpp7TeV->Draw("same,pE2");
			graphCombEtaToPi0Ratiopp7TeVNoXErrors->Draw("same,pe");
            		
// 			TLegend* legendEtatoPi0combo_onlyPbPbwithPP = new TLegend(0.12,0.73,0.53,0.92);
// 			legendEtatoPi0combo_onlyPbPbwithPP->SetFillColor(0);
// 			legendEtatoPi0combo_onlyPbPbwithPP->SetLineColor(0);
// 			legendEtatoPi0combo_onlyPbPbwithPP->SetTextFont(42);
// 			legendEtatoPi0combo_onlyPbPbwithPP->SetTextSize(0.04);
// 			legendEtatoPi0combo_onlyPbPbwithPP->SetMargin(0.17);
// 			legendEtatoPi0combo_onlyPbPbwithPP->SetHeader(collisionSystem2760GeV.Data());
// 			legendEtatoPi0combo_onlyPbPbwithPP->AddEntry(graphCombEtatoPi0Sys2760GeVLHC11h_0010,Form("  %s",cent0010.Data()),"fp");
// 			legendEtatoPi0combo_onlyPbPbwithPP->AddEntry(graphCombEtatoPi0Sys2760GeVLHC11h_2050,Form("%s",cent2050.Data()),"fp");
// 			legendEtatoPi0combo_onlyPbPbwithPP->Draw();
            legendEtatoPi0combo_onlyPbPb->Draw();

			
			TLegend* legendEtatoPi0combo_withPP = new TLegend(0.55,0.15,0.95,0.26);
			legendEtatoPi0combo_withPP->SetFillColor(0);
			legendEtatoPi0combo_withPP->SetLineColor(0);
			legendEtatoPi0combo_withPP->SetTextFont(42);
			legendEtatoPi0combo_withPP->SetTextSize(0.037);
			legendEtatoPi0combo_withPP->SetMargin(0.17);
			legendEtatoPi0combo_withPP->SetHeader(Form("#eta/#pi^{0} %s",collisionSystemPP7TeV.Data()));
			legendEtatoPi0combo_withPP->AddEntry(graphCombEtaToPi0RatioSysErrpp7TeV,"PLB 717 (2012) 162","fp");//"Phys. Lett. B 717 (2012) 162-172","fp");
			legendEtatoPi0combo_withPP->Draw();

			
		canvasEtatoPi0combo->SaveAs(Form("%s/%s_EtatoPi0Ratio_combined_DataOnlyWithPP.%s",outputDir.Data(),meson.Data(),suffix.Data()));
		canvasEtatoPi0combo->SaveAs(Form("%s/%s_EtatoPi0Ratio_combined_DataOnlyWithPP.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
		
		
		canvasEtatoPi0combo->cd();
			histo2DEtatoPi0combo->Draw("copy");
			
			DrawGammaSetMarkerTGraphAsym(graphEtaToPi0JetQuenching_0010, 0, 0, colorNLO0010,colorNLO0010, widthLinesBoxes, kTRUE, colorNLO0010);
			graphEtaToPi0JetQuenching_0010->Draw("3,same");

			DrawGammaSetMarker(histoEtaToPi0EPOSRebin_0010,1,2,colorEPOS0010,colorEPOS0010);
			DrawGammaSetMarker(histoEtaToPi0EPOSRebin_2050,1,2,colorEPOS2050,colorEPOS2050);
			histoEtaToPi0EPOSRebin_0010->SetLineWidth(2);
			histoEtaToPi0EPOSRebin_2050->SetLineWidth(2);
			
			histoEtaToPi0EPOSRebin_0010->GetXaxis()->SetRangeUser(0.,13.9);
// 			histoEtaToPi0EPOSRebin_2050->GetXaxis()->SetRangeUser(0.,20.);
// 			histoEtaToPi0EPOSRebin_0010->Draw("same,c,histo");
// 			histoEtaToPi0EPOSRebin_2050->Draw("same,c,histo");
// 
			DrawGammaSetMarkerTGraphAsym(TheoryCracowEtaToPi0LowPt_0010, 0, 0, colorCracow0010,colorCracow0010, 5, kTRUE,colorCracow0010);
// 			TheoryCracowEtaToPi0LowPt_0010->Draw("l,same");
			DrawGammaSetMarkerTGraphAsym(TheoryCracowEtaToPi0LowPt_2050, 0, 0, colorCracow2050,colorCracow2050, 5, kTRUE,colorCracow2050);
// 			TheoryCracowEtaToPi0LowPt_2050->Draw("l,same");

			graphCombEtatoPi0Sys2760GeVLHC11h_0010->Draw("E2same");
			if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombEtatoPi0Stat2760GeVLHC11h_0010);
			graphCombEtatoPi0Stat2760GeVLHC11h_0010->Draw("p,same,e0");

// 			graphCombEtatoPi0Sys2760GeVLHC11h_2050->Draw("E2same");
// 			if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombEtatoPi0Stat2760GeVLHC11h_2050);
// 			graphCombEtatoPi0Stat2760GeVLHC11h_2050->Draw("p,same,e0");
			
			labelEtaToPi0Energy->Draw();
            
            
            TLegend* legendEtatoPi0combo_onlyPbPb2 = new TLegend(0.12,0.8,0.53,0.92);
            legendEtatoPi0combo_onlyPbPb2->SetFillColor(0);
            legendEtatoPi0combo_onlyPbPb2->SetLineColor(0);
            legendEtatoPi0combo_onlyPbPb2->SetTextFont(42);
            legendEtatoPi0combo_onlyPbPb2->SetTextSize(0.037);
            legendEtatoPi0combo_onlyPbPb2->SetMargin(0.17);
            legendEtatoPi0combo_onlyPbPb2->SetHeader(collisionSystem2760GeV.Data());
//          legendEtatoPi0combo_onlyPbPb->AddEntry(graphCombEtatoPi0Sys2760GeVLHC11h_0010,Form("%s %s",cent0010.Data(),collisionSystem2760GeV.Data()),"fp");
//          legendEtatoPi0combo_onlyPbPb->AddEntry(graphCombEtatoPi0Sys2760GeVLHC11h_2050,Form("%s %s",cent2050.Data(),collisionSystem2760GeV.Data()),"fp");
            legendEtatoPi0combo_onlyPbPb2->AddEntry(graphCombEtatoPi0Sys2760GeVLHC11h_0010,Form("  %s",cent0010.Data()),"fp");
			legendEtatoPi0combo_onlyPbPb2->Draw();

			TLegend* legendRatioALICE2 = new TLegend(0.54,0.75,0.92,0.91); //0.12,0.62,0.5,0.71);
			legendRatioALICE2->SetFillColor(0);
			legendRatioALICE2->SetLineColor(0);
			legendRatioALICE2->SetTextFont(42);
			legendRatioALICE2->SetTextSize(0.037);
			legendRatioALICE2->SetHeader("NLO DCZW (#tau_{0} = 0.6 fm)");
			legendRatioALICE2->AddEntry(graphEtaToPi0JetQuenching_0010,"0-10%","f");
			legendRatioALICE2->AddEntry((TObject*)0,"arXiv:1506.00838","");
			legendRatioALICE2->Draw();

			
			TLegend* legendRatioALICE3A = new TLegend(0.41,0.78,0.76,0.9);
			legendRatioALICE3A->SetFillColor(0);
			legendRatioALICE3A->SetLineColor(0);
			legendRatioALICE3A->SetTextFont(42);
			legendRatioALICE3A->SetTextSize(0.04);
			legendRatioALICE3A->SetHeader("Phys.Rev. C89, 064903 (2014)");
			legendRatioALICE3A->AddEntry(histoEPOSRebin_0010,"EPOS 0-10%", "l");
			legendRatioALICE3A->AddEntry(histoEPOSRebin_2050,"EPOS 20-50%", "l");
// 			legendRatioALICE3A->Draw();
			TLegend* legendRatioALICE3B = new TLegend(0.41,0.66,0.76,0.78);
			legendRatioALICE3B->SetFillColor(0);
			legendRatioALICE3B->SetLineColor(0);
			legendRatioALICE3B->SetTextFont(42);
			legendRatioALICE3B->SetTextSize(0.04);
			legendRatioALICE3B->SetHeader("Phys.Rev. C90, 014906 (2014)");
			legendRatioALICE3B->AddEntry(TheoryCracowEtaToPi0LowPt_0010,"NEQ SHM 0-10%", "l");
			legendRatioALICE3B->AddEntry(TheoryCracowEtaToPi0LowPt_2050,"NEQ SHM 20-50%", "l");
// 			legendRatioALICE3B->Draw();
			
		canvasEtatoPi0combo->SaveAs(Form("%s/%s_EtatoPi0Ratio_combined.%s",outputDir.Data(),meson.Data(),suffix.Data()));
		canvasEtatoPi0combo->SaveAs(Form("%s/%s_EtatoPi0Ratio_combined.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
		
				
		canvasEtatoPi0combo->cd();
			histo2DEtatoPi0combo->Draw("copy");
			histo2DEtatoPi0combo->GetYaxis()->SetTitleOffset(1.1);
			histo2DEtatoPi0combo->GetYaxis()->SetTitle("Particle ratio");
				
			DrawGammaSetMarkerTGraphAsym(graphChargedRatioKaonToPion0010, 24,markerSizeComb, colorCharged, colorCharged);
			DrawGammaSetMarkerTGraphAsym(graphChargedRatioKaonToPionSys0010,24,markerSizeComb, colorCharged, colorCharged, 1, kTRUE);
			graphChargedRatioKaonToPionSys0010->Draw("2same");
			graphChargedRatioKaonToPion0010->Draw("p,same"); 
			
			TLegend* legendChargedRatio = new TLegend(0.12,0.76,0.95,0.92); 
			legendChargedRatio->SetFillColor(0);
			legendChargedRatio->SetLineColor(0);
			legendChargedRatio->SetTextFont(42);
			legendChargedRatio->SetTextSize(0.037);
			legendChargedRatio->SetMargin(0.17);
			legendChargedRatio->SetHeader(collisionSystem2760GeV.Data());
			legendChargedRatio->SetNColumns(2);
			legendChargedRatio->AddEntry(graphCombEtatoPi0Sys2760GeVLHC11h_0010,Form("#eta/#pi^{0},   %s",cent0010.Data()),"fp");
			legendChargedRatio->AddEntry(graphChargedRatioKaonToPionSys0010,Form("K^{#pm}/#pi^{#pm}, %s",cent0010.Data()),"fp");
			legendChargedRatio->AddEntry(graphCombEtatoPi0Sys2760GeVLHC11h_2050,Form("#eta/#pi^{0}, %s",cent2050.Data()),"fp");
			legendChargedRatio->AddEntry((TObject*)0,"PLB 736 (2014) 196",""); 
			legendChargedRatio->Draw();
			
			graphCombEtatoPi0Sys2760GeVLHC11h_0010->Draw("E2same");
			if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombEtatoPi0Stat2760GeVLHC11h_0010);
			graphCombEtatoPi0Stat2760GeVLHC11h_0010->Draw("p,same,e0");

			graphCombEtatoPi0Sys2760GeVLHC11h_2050->Draw("E2same");
			if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombEtatoPi0Stat2760GeVLHC11h_2050);
			graphCombEtatoPi0Stat2760GeVLHC11h_2050->Draw("p,same,e0");

		canvasEtatoPi0combo->SaveAs(Form("%s/%s_EtatoPi0Ratio_combined_withKaonsToPions.%s",outputDir.Data(),meson.Data(),suffix.Data()));
		canvasEtatoPi0combo->SaveAs(Form("%s/%s_EtatoPi0Ratio_combined_withKaonsToPions.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
	
        
        
        canvasEtatoPi0combo->cd();
            histo2DEtatoPi0combo->Draw("copy");
            histo2DEtatoPi0combo->GetYaxis()->SetRangeUser(0.,1.15);

            graphChargedRatioKaonToPionSys0010->Draw("2same");
            graphChargedRatioKaonToPion0010->Draw("p,same"); 
            
            graphCombEtaToPi0RatioSysErrpp7TeV->Draw("same,pE2");
            graphCombEtaToPi0Ratiopp7TeVNoXErrors->Draw("same,pe");

            TLegend* legendChargedRatio2 = new TLegend(0.12,0.7,0.48,0.92); 
            legendChargedRatio2->SetFillColor(0);
            legendChargedRatio2->SetLineColor(0);
            legendChargedRatio2->SetTextFont(42);
            legendChargedRatio2->SetTextSize(0.037);
            legendChargedRatio2->SetMargin(0.17);
            legendChargedRatio2->SetHeader(collisionSystemPbPb0010.Data());
            legendChargedRatio2->AddEntry(graphCombEtatoPi0Sys2760GeVLHC11h_0010,"#eta/#pi^{0}","fp");
//             legendChargedRatio2->AddEntry((TObject*)0,"","");
            legendChargedRatio2->AddEntry(graphChargedRatioKaonToPionSys0010,"K^{#pm}/#pi^{#pm}","fp");
            legendChargedRatio2->AddEntry((TObject*)0,"PLB 736 (2014) 196",""); 
            legendChargedRatio2->Draw();
            
            legendEtatoPi0combo_withPP->Draw();

            
            graphCombEtatoPi0Sys2760GeVLHC11h_0010->Draw("E2same");
            if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombEtatoPi0Stat2760GeVLHC11h_0010);
            graphCombEtatoPi0Stat2760GeVLHC11h_0010->Draw("p,same,e0");

//             graphCombEtatoPi0Sys2760GeVLHC11h_2050->Draw("E2same");
//             if(noXerrorBars) ProduceGraphAsymmWithoutXErrors(graphCombEtatoPi0Stat2760GeVLHC11h_2050);
//             graphCombEtatoPi0Stat2760GeVLHC11h_2050->Draw("p,same,e0");

        canvasEtatoPi0combo->SaveAs(Form("%s/%s_EtatoPi0Ratio_withPPandKaonsToPions.%s",outputDir.Data(),meson.Data(),suffix.Data()));
        canvasEtatoPi0combo->SaveAs(Form("%s/%s_EtatoPi0Ratio_withPPandKaonsToPions.%s",paperPlots.Data(),meson.Data(),suffix.Data()));
        
	}
	

	
	
	
	
	
	
	
	
//     Double_t arrayBoundariesX1_MassWidth[3];
//     Double_t arrayBoundariesY1_MassWidth[3];
//     Double_t relativeMarginsX_MassWidth[3];
//     Double_t relativeMarginsY_MassWidth[3];
//     textSizeLabelsPixel = 90;
//     ReturnCorrectValuesForCanvasScaling(2400,1600, 2, 2,0.25, 0.2, 0.003,0.2,arrayBoundariesX1_MassWidth,arrayBoundariesY1_MassWidth,relativeMarginsX_MassWidth,relativeMarginsY_MassWidth);
// 	
//     TCanvas * canvas4PartMassWidth = new TCanvas("canvas4PartMassWidth","",0,0,2400,1600);  // gives the page size        
// //     DrawGammaCanvasSettings( canvas4PartMassWidth, 0.13, 0.0, 0.02, 0.1);
//     canvas4PartMassWidth->cd();
// 
//     TPad* pad4PartMassWidth1 = new TPad("pad4PartMassWidth1", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
//     DrawGammaPadSettings( pad4PartMassWidth1, relativeMarginsX[0], relativeMarginsX[1], relativeMarginsY[0], relativeMarginsY[1]);
//     pad4PartMassWidth1->Draw();
// 
//     TPad* pad4PartMassWidth2 = new TPad("pad4PartMassWidth2", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
//     DrawGammaPadSettings( pad4PartMassWidth2, relativeMarginsX[0], relativeMarginsX[1], relativeMarginsY[1], relativeMarginsY[2]);
//     pad4PartMassWidth2->Draw();
//     
//     TPad* pad4PartMassWidth3 = new TPad("pad4PartMassWidth3", "", arrayBoundariesX1_4[1], arrayBoundariesY1_4[1], arrayBoundariesX1_4[2], arrayBoundariesY1_4[0],-1, -1, -2);
//     DrawGammaPadSettings( pad4PartMassWidth3, relativeMarginsX[1], relativeMarginsX[1], relativeMarginsY[0], relativeMarginsY[1]);
//     pad4PartMassWidth3->Draw();
// 
//     TPad* pad4PartMassWidth4 = new TPad("pad4PartMassWidth4", "", arrayBoundariesX1_4[1], arrayBoundariesY1_4[2], arrayBoundariesX1_4[2], arrayBoundariesY1_4[1],-1, -1, -2);
//     DrawGammaPadSettings( pad4PartMassWidth4, relativeMarginsX[1], relativeMarginsX[1], relativeMarginsY[1], relativeMarginsY[2]);
//     pad4PartMassWidth4->Draw();
//  
// 
// 
//     Int_t textSizeLabelsPixelMass = 50;
//     Double_t marginMass = 0.13*2400;
//     Double_t textsizeLabelsMass = 0;
//     Double_t textsizeFacMass = 0;
//     Double_t textsizeLabelsWidth = 0;
//     Double_t textsizeFacWidth = 0;
// 
//     if (pad4PartMassWidth1->XtoPixel(pad4PartMassWidth1->GetX2()) < pad4PartMassWidth1->YtoPixel(pad4PartMassWidth1->GetY1())){
//         textsizeLabelsWidth = (Double_t)textSizeLabelsPixelMass/pad4PartMassWidth1->XtoPixel(pad4PartMassWidth1->GetX2()) ;
//         textsizeFacWidth = (Double_t)1./pad4PartMassWidth1->XtoPixel(pad4PartMassWidth1->GetX2()) ;
//     } else {
//         textsizeLabelsWidth = (Double_t)textSizeLabelsPixelMass/pad4PartMassWidth1->YtoPixel(pad4PartMassWidth1->GetY1());
//         textsizeFacWidth = (Double_t)1./pad4PartMassWidth1->YtoPixel(pad4PartMassWidth1->GetY1());
//     }
//     if (pad4PartMassWidth2->XtoPixel(pad4PartMassWidth2->GetX2()) < pad4PartMassWidth2->YtoPixel(pad4PartMassWidth2->GetY1())){
//         textsizeLabelsMass = (Double_t)textSizeLabelsPixelMass/pad4PartMassWidth2->XtoPixel(pad4PartMassWidth2->GetX2()) ;
//         textsizeFacMass = (Double_t)1./pad4PartMassWidth2->XtoPixel(pad4PartMassWidth2->GetX2()) ;
//     } else {
//         textsizeLabelsMass = (Double_t)textSizeLabelsPixelMass/pad4PartMassWidth2->YtoPixel(pad4PartMassWidth2->GetY1());
//         textsizeFacMass = (Double_t)1./pad4PartMassWidth2->YtoPixel(pad4PartMassWidth2->GetY1());
//     }
// 
//     cout << textsizeLabelsMass << endl;
// 
//     TPad* padFWHMLegend1 = new TPad("padFWHMLegend1", "", 0.07, 0.815, 0.35, 0.93,-1, -1, -2);
//     DrawGammaPadSettings( padFWHMLegend1, 0., 0., 0., 0.);
//     padFWHMLegend1->Draw();
// 
// 
//     TH2D *histo2DPi0FWHM;
//     histo2DPi0FWHM = new TH2D("histo2DPi0FWHM", "histo2DPi0FWHM", 20,0.35,20. ,1000.,-30,40);
//     SetStyleHistoTH2ForGraphs(histo2DPi0FWHM, "#it{p}_{T} (GeV/#it{c})","peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth,textsizeLabelsWidth, 0.85*textsizeLabelsWidth,textsizeLabelsWidth, 1,0.3/(textsizeFacWidth*marginMass), 515, 504); 
//     histo2DPi0FWHM->GetYaxis()->SetRangeUser(-1.,18);
//     histo2DPi0FWHM->GetYaxis()->SetLabelOffset(0.01);
//     histo2DPi0FWHM->GetYaxis()->SetLabelFont(42);
//     histo2DPi0FWHM->GetXaxis()->SetLabelFont(42);
// 
//     TH2D *histo2DPi0Mass;
//     histo2DPi0Mass = new TH2D("histo2DPi0Mass", "histo2DPi0Mass", 20,0.35,20. ,1000.,125.,150);
//     SetStyleHistoTH2ForGraphs(histo2DPi0Mass, "#it{p}_{T} (GeV/#it{c})","peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass,textsizeLabelsMass, 0.85*textsizeLabelsMass,textsizeLabelsMass, 0.9,0.3/(textsizeFacMass*marginMass), 515, 510); 
//     histo2DPi0Mass->GetYaxis()->SetRangeUser(128.,143.5);
//     histo2DPi0Mass->GetXaxis()->SetLabelOffset(-0.02);
//     histo2DPi0Mass->GetYaxis()->SetLabelOffset(0.01);
//     histo2DPi0Mass->GetYaxis()->SetLabelFont(42);
//     histo2DPi0Mass->GetXaxis()->SetLabelFont(42);
// 
//     pad4PartMassWidth1->cd();
//     pad4PartMassWidth1->SetLogx();
//     histo2DPi0FWHM->DrawCopy();
//         
// //     DrawGammaSetMarker(histoPHOSWidthDataPP, markerStylePHOS, markerSizeMass, colorPHOSMass, colorPHOSMass);
// //     histoPHOSWidthDataPP->DrawCopy("same,e0,p"); 
// //     DrawGammaSetMarker(histoPHOSWidthMCPP, markerStylePHOSMC, markerSizeMass, colorPHOSMCMass , colorPHOSMCMass);
// //     histoPHOSWidthMCPP->DrawCopy("same,e0,p"); 
// // 
// //     DrawGammaSetMarker(histoPCMWidthDataPP, markerStyleConv, markerSizeMass, colorConv, colorConv);
// //     histoPCMWidthDataPP->DrawCopy("same,e0,p"); 
// //     DrawGammaSetMarker(histoPCMWidthMCPP, markerStyleConvMC, markerSizeMass, colorConvMC, colorConvMC);
// //     histoPCMWidthMCPP->DrawCopy("same,e0,p"); 
// // 
// //     DrawGammaSetMarker(histoPHOSWidthData0010, markerStylePHOS, markerSizeMass, colorPHOSMass, colorPHOSMass);
// //     DrawGammaSetMarker(histoPHOSWidthMC0010, markerStylePHOSMC, markerSizeMass, colorPHOSMCMass , colorPHOSMCMass);
// // 
// //     TLatex *labelMassPi0PP = new TLatex(0.2,0.88,collisionSystemPP.Data());
// //     SetStyleTLatex( labelMassPi0PP, 0.85*textsizeLabelsWidth,4);
// //     labelMassPi0PP->Draw();
// //     TLatex *labelLegendAMass = new TLatex(0.92,0.88,"a)");
// //     SetStyleTLatex( labelLegendAMass,0.85*textsizeLabelsWidth,4);
// //     labelLegendAMass->Draw();
// // 
// //     //********************************** Defintion of the Legend ************************************************** 
// //     Double_t columnsLegendFWHM[4]   = {0.,0.2,0.37,0.55};
// //     Double_t rowsLegendFWHM[3]      = {0.66,0.33,0.0};
// //     //******************* Text sizes *******************
// //     Size_t textSizeLeftColumnFWHM   = 0.301;
// //     Size_t textSizeTopRowFWHM   = 0.301; 
// //     Size_t textSizeSecondRowFWHM    = 0.301;
// //     //******************* Offsets ***********************
// //     Double_t offsetMarkerXFWHM  = 0.07;
// //     Double_t offsetMarkerYFWHM  = 0.07;
// //     //****************** Scale factors ******************
// //     Double_t scaleMarkerFWHM        = 1.;
// // 
// //     padFWHMLegend1->cd();
// //     //****************** first Column **************************************************
// //     TLatex *textFWHMCTS = new TLatex(columnsLegendFWHM[0],rowsLegendFWHM[1],"PCM");
// //     SetStyleTLatex( textFWHMCTS, textSizeLeftColumnFWHM,4);
// //     textFWHMCTS->Draw();
// //     TLatex *textFWHMPHOS = new TLatex(columnsLegendFWHM[0],rowsLegendFWHM[2],"PHOS");
// //     SetStyleTLatex( textFWHMPHOS, textSizeLeftColumnFWHM,4);
// //     textFWHMPHOS->Draw();
// // 
// //     //****************** second Column *************************************************
// //     TLatex *textFWHMData2 = new TLatex(columnsLegendFWHM[1],rowsLegendFWHM[0] ,"Data");
// //     SetStyleTLatex( textFWHMData2, textSizeTopRowFWHM ,4);
// //     textFWHMData2->Draw();
// //     TLatex *textFWHMMC2 = new TLatex(columnsLegendFWHM[2] ,rowsLegendFWHM[0],"MC");
// //     SetStyleTLatex( textFWHMMC2, textSizeTopRowFWHM,4);
// //     textFWHMMC2->Draw();
// // 
// //     TMarker* markerCTSPi0FWHM = CreateMarkerFromHisto(histoPCMWidthDataPP,columnsLegendFWHM[1]+ offsetMarkerXFWHM ,rowsLegendFWHM[1]+ offsetMarkerYFWHM ,scaleMarkerFWHM);
// //     markerCTSPi0FWHM->DrawMarker(columnsLegendFWHM[1]+ offsetMarkerXFWHM ,rowsLegendFWHM[1]+ offsetMarkerYFWHM);
// //     TMarker* markerPHOSPi0FWHM = CreateMarkerFromHisto(histoPHOSWidthData0010,columnsLegendFWHM[1]+ offsetMarkerXFWHM ,rowsLegendFWHM[2]+ offsetMarkerYFWHM ,scaleMarkerFWHM);
// //     markerPHOSPi0FWHM->DrawMarker(columnsLegendFWHM[1]+ offsetMarkerXFWHM ,rowsLegendFWHM[2]+ offsetMarkerYFWHM);
// // 
// //     TMarker* markerCTSPi0FWHMMC = CreateMarkerFromHisto(histoPCMWidthMCPP,columnsLegendFWHM[2]+ offsetMarkerXFWHM ,rowsLegendFWHM[1]+ offsetMarkerYFWHM ,scaleMarkerFWHM);
// //     markerCTSPi0FWHMMC->DrawMarker(columnsLegendFWHM[2]+ offsetMarkerXFWHM ,rowsLegendFWHM[1]+ offsetMarkerYFWHM);
// //     TMarker* markerPHOSPi0FWHMMC = CreateMarkerFromHisto(histoPHOSWidthMC0010,columnsLegendFWHM[2]+ offsetMarkerXFWHM ,rowsLegendFWHM[2]+ offsetMarkerYFWHM ,scaleMarkerFWHM);
// //     markerPHOSPi0FWHMMC->DrawMarker(columnsLegendFWHM[2]+ offsetMarkerXFWHM ,rowsLegendFWHM[2]+ offsetMarkerYFWHM);
// // 
// //     TLatex *textWidthConv2 = new TLatex(columnsLegendFWHM[3],rowsLegendFWHM[1] ,"FWHM/2.35");
// //     SetStyleTLatex( textWidthConv2, textSizeSecondRowFWHM,4);
// //     textWidthConv2->Draw();
// //     TLatex *textWidthPHOS2 = new TLatex(columnsLegendFWHM[3] ,rowsLegendFWHM[2],"#sigma");
// //     SetStyleTLatex( textWidthPHOS2, textSizeSecondRowFWHM,4);
// //     textWidthPHOS2->Draw();
// 
// 
//     pad4PartMassWidth1->cd();
//     pad4PartMassWidth1->SetLogx();
//     histo2DPi0FWHM->DrawCopy();
//         
//     DrawGammaSetMarker(histoPCMPi0FWHMMeVPbPb2760GeV_0010, markerStyleConv, markerSizeMass, colorConv, colorConv);
//     histoPCMPi0FWHMMeVPbPb2760GeV_0010->DrawCopy("same,e0,p"); 
//     DrawGammaSetMarker(histoPCMPi0TrueFWHMMeVPbPb2760GeV_0010, markerStyleConvMC, markerSizeMass, colorConvMC, colorConvMC);
//     histoPCMPi0TrueFWHMMeVPbPb2760GeV_0010->DrawCopy("same,e0,p"); 
// 
//     TLatex *labelMassPi0PbPb0005 = new TLatex(0.05,0.88,collisionSystemPbPb0010.Data());
//     SetStyleTLatex( labelMassPi0PbPb0005, 0.85*textsizeLabelsWidth,4);
//     labelMassPi0PbPb0005->Draw();
//     TLatex *labelLegendBMass = new TLatex(0.89,0.88,"c)");
//     SetStyleTLatex( labelLegendBMass, 0.85*textsizeLabelsWidth,4);
//     labelLegendBMass->Draw();
// 
//     pad4PartMassWidth1->Update();
//     pad4PartMassWidth2->cd();
//     pad4PartMassWidth2->SetLogx();
//     histo2DPi0Mass->DrawCopy();
//         
//     DrawGammaSetMarker(histoPCMPi0MassPbPb2760GeV_0010 , markerStyleConv, markerSizeMass, colorConv, colorConv);                    
//     histoPCMPi0MassPbPb2760GeV_0010->DrawCopy("same,e0,p"); 
//     DrawGammaSetMarker(histoPCMPi0TrueMassPbPb2760GeV_0010 , markerStyleConvMC , markerSizeMass, colorConvMC, colorConvMC);                   
//     histoPCMPi0TrueMassPbPb2760GeV_0010->DrawCopy("same,e0,p"); 
//     
//     DrawGammaLines(0.35, 20. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,1.,colorConv);
//     TLatex *labelLegendFMass = new TLatex(0.89,0.9,"f)");
//     SetStyleTLatex( labelLegendFMass, 0.85*textsizeLabelsMass,4);
//     labelLegendFMass->Draw();
// 
//     pad4PartMassWidth2->Update();
// 
//     pad4PartMassWidth3->cd();
//     pad4PartMassWidth3->SetLogx();
//     histo2DPi0FWHM->DrawCopy();
//         
//     DrawGammaSetMarker(histoPCMEtaFWHMMeVPbPb2760GeV_0010, markerStyleConv, markerSizeMass, colorConv, colorConv);
//     histoPCMEtaFWHMMeVPbPb2760GeV_0010->DrawCopy("same,e0,p"); 
//     DrawGammaSetMarker(histoPCMEtaTrueFWHMMeVPbPb2760GeV_0010, markerStyleConvMC, markerSizeMass, colorConvMC, colorConvMC);
//     histoPCMEtaTrueFWHMMeVPbPb2760GeV_0010->DrawCopy("same,e0,p"); 
// 
//     TLatex *labelMassPi0PbPb6080 = new TLatex(0.05,0.88,collisionSystemCent6080.Data());
//     SetStyleTLatex( labelMassPi0PbPb6080, 0.85*textsizeLabelsWidth,4);
//     labelMassPi0PbPb6080->Draw();
//     TLatex *labelLegendCMass = new TLatex(0.91,0.88,"b)");
//     SetStyleTLatex( labelLegendCMass, 0.85*textsizeLabelsWidth,4);
//     labelLegendCMass->Draw();
// 
//     pad4PartMassWidth3->Update();
//     pad4PartMassWidth4->cd();
//     pad4PartMassWidth4->SetLogx();
//     histo2DPi0Mass->DrawCopy();
// 
//     DrawGammaSetMarker(histoPCMPi0MassPbPb2760GeV_0010 , markerStyleConv, markerSizeMass, colorConv, colorConv);                    
//     histoPCMPi0MassPbPb2760GeV_0010->DrawCopy("same,e0,p"); 
//     DrawGammaSetMarker(histoPCMPi0TrueMassPbPb2760GeV_0010 , markerStyleConvMC , markerSizeMass, colorConvMC, colorConvMC);                   
//     histoPCMPi0TrueMassPbPb2760GeV_0010->DrawCopy("same,e0,p"); 
// 
//     DrawGammaLines(0.35, 20. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,1.,colorConv);
//     TLatex *labelLegendEMass = new TLatex(0.91,0.9,"e)");
//     SetStyleTLatex( labelLegendEMass, 0.85*textsizeLabelsMass,4);
//     labelLegendEMass->Draw();
// 
//     pad4PartMassWidth4->Update();
// 
//     canvas4PartMassWidth->Update(); 
//     canvas4PartMassWidth->SaveAs(Form("%s/MassWidth_4Parted_Paper.%s",paperPlots.Data(),suffix.Data()));
//     canvas4PartMassWidth->SaveAs(Form("%s/MassWidth_4Parted_Paper.%s",outputDir.Data(),suffix.Data()));
//     delete pad4PartMassWidth1;  
//     delete pad4PartMassWidth2;  
//     delete pad4PartMassWidth3;  
//     delete pad4PartMassWidth4;  
//     delete canvas4PartMassWidth;    
    
    
    cout << WriteParameterToFile(fitTCMInvYield2760GeVLHC11h_0010)<< endl << endl;
    cout << WriteParameterToFile(fitTCMInvYield2760GeVLHC11h_2050)<< endl << endl;

	
		
	//******************************************************************************//
	//************************* Saving of final results ****************************//
	//******************************************************************************//
	
	TFile fCombResults(Form("%s/CombinedResultsPaperPbPb2760GeVLHC11h_%s.root", outputDir.Data(),dateForOutput.Data()), "UPDATE");

		graphCombInvYieldTot2760GeVALHC11h_0010->Write(Form("graphInvYield%sComb2760GeVLHC11h_0010",meson.Data()));
		graphCombInvYieldStat2760GeVALHC11h_0010->Write(Form("graphInvYield%sComb2760GeVLHC11hStatErr_0010",meson.Data()));
		graphCombInvYieldSys2760GeVALHC11h_0010->Write(Form("graphInvYield%sComb2760GeVLHC11hSysErr_0010",meson.Data()));    
		
		graphCombInvYieldTot2760GeVALHC11h_2050->Write(Form("graphInvYield%sComb2760GeVLHC11h_2050",meson.Data()));
		graphCombInvYieldStat2760GeVALHC11h_2050->Write(Form("graphInvYield%sComb2760GeVLHC11hStatErr_2050",meson.Data()));
		graphCombInvYieldSys2760GeVALHC11h_2050->Write(Form("graphInvYield%sComb2760GeVLHC11hSysErr_2050",meson.Data()));      

		fitTCMInvYield2760GeVLHC11h_0010->Write(Form("FitToYield%s_0010",meson.Data()));
		fitTCMInvYield2760GeVLHC11h_2050->Write(Form("FitToYield%s_2050",meson.Data()));
		
        histoEPOSRebin_0010->Write(Form("EPOSprediction%s_0010",meson.Data()));
        histoEPOSRebin_2050->Write(Form("EPOSprediction%s_2050",meson.Data()));
        
        TheoryCracowLowPtforRatio_0010->Write(Form("Cracowprediction%s_2050",meson.Data()));
        TheoryCracowLowPtforRatio_2050->Write(Form("Cracowprediction%s_2050",meson.Data()));

                    
        graphPCMPi0CorrYieldSysErrPbPb2760GeV_0010->Write("graphPCMPi0CorrYieldSysErrPbPb2760GeV_0010");
        graphPCMPi0InvYieldStatPbPb2760GeV_0010->Write("graphPCMPi0InvYieldStatPbPb2760GeV_0010");
        
        graphPCMPi0CorrYieldSysErrPbPb2760GeV_2040->Write("graphPCMPi0CorrYieldSysErrPbPb2760GeV_2040");
        graphPCMPi0InvYieldStatPbPb2760GeV_2040->Write("graphPCMPi0InvYieldStatPbPb2760GeV_2040");

        graphEMCalPi0CorrYieldSysErrPbPb2760GeV_0010->Write("graphEMCalPi0CorrYieldSysErrPbPb2760GeV_0010");
        graphEMCalPi0InvYieldStatPbPb2760GeV_0010->Write("graphEMCalPi0InvYieldStatPbPb2760GeV_0010");

// 		if(meson.CompareTo("Eta")==0){
// 		
// 			graphCombEtatoPi0Sys2760GeVLHC11h_0010->Write("graphCombEtatoPi0Sys2760GeVLHC11h_0010");
// 			graphCombEtatoPi0Stat2760GeVLHC11h_0010->Write("graphCombEtatoPi0Stat2760GeVLHC11h_0010");
// 
// 			graphCombEtatoPi0Sys2760GeVLHC11h_2050->Write("graphCombEtatoPi0Sys2760GeVLHC11h_2050");
// 			graphCombEtatoPi0Stat2760GeVLHC11h_2050->Write("graphCombEtatoPi0Stat2760GeVLHC11h_2050");
// 
// 		}
		
	fCombResults.Close();
		




		
		
}
	
