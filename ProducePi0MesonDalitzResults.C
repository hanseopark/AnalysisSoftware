/****************************************************************************************************************************
****** 		provided by Gamma Conversion Group, PWG4,												*****
******		Pedro Gonzalez, pedro.gonzalez.zamora@cern.ch
*****		Ana Marin, marin@physi.uni-heidelberg.de													*****												*****													*****
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
#include "ProducePi0MesonDalitzResults.h"

Int_t ProducePi0MesonDalitzResults(TString fileNameppDalitz 	= "ExternalInput/Dalitz/data_PCMDalitzResultsFullCorrection_PP_NoBinShifting_2016-03-19.root", 
				  TString fileNamepPbDalitz 	= "ExternalInputpPb/PCM/data_PCMResults_Dalitz_pPb_20150806.root",  
				  TString suffix 		= "pdf", 
				  TString thesisPlots 		= "", 
				  TString bWCorrection="X",
				  TString isMC = "kTRUE"
 				){	

	//data_PCMDalitzResultsFullCorrection_PP_NoBinShifting_2016-03-19.root
	TString date = ReturnDateString();
	
	
	gROOT->Reset();	
	gROOT->SetStyle("Plain");
	
	StyleSettingsThesis(suffix);	
	SetPlotStyle();
	
	TString dateForOutput 			= ReturnDateStringForOutput();
	cout << dateForOutput.Data() << endl;
	//___________________________________ Declaration of files _____________________________________________
	TString collisionSystem7TeV 			= "pp, #sqrt{#it{s}} = 7 TeV";		
	TString collisionSystem2760GeV			= "pp, #sqrt{#it{s}} = 2.76 TeV";
	TString collisionSystem5023GeV			= "p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV";
	TString invYieldLabel 				= "#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}N_{#pi^{0}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})";
	TString pTLabel 				= "#it{p}_{T} (GeV/#it{c})";
	TString RpPbLabel				= "#it{R}_{pPb}";
	TString energyLabel				= "p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV";
	TString DalitzLabel 				= "#pi^{0} #rightarrow e^{+}e^{-}#gamma";
	TString ChargedPionsLabel 			= "#pi^{+} + #pi^{-}";
	TString thesisPlotLabel 			= "This thesis";
	
	
	
	Width_t						widthLinesBoxes		= 1.4;
	
		
	TString outputDir 				= Form("Pi0MesonDalitzResults%s/%s/%s/",bWCorrection.Data(),dateForOutput.Data(),suffix.Data());
	
	
	cout << outputDir.Data() << endl;
	cout << fileNameppDalitz.Data() << endl;
	
	Int_t color[20] = {kGreen+1,kBlue,kRed,620,880,591,432,422,435,420,407,416,830,404,403,802,634,608,920,923};	

	gSystem->Exec("mkdir -p "+outputDir);
 	
	gSystem->Exec(Form("cp %s %s/InputppDalitz.root",  fileNameppDalitz.Data(), outputDir.Data()));
	gSystem->Exec(Form("cp %s %s/InputpPbDalitz.root", fileNamepPbDalitz.Data(), outputDir.Data()));
	
	
	
	Bool_t thesis						= kFALSE;
	
	if(thesisPlots.CompareTo("thesis") == 0){
	  
		thesis 						= kTRUE;
		
	}
	
	TString prefix2						= "";
	
	
	if (isMC.CompareTo("kTRUE")==0){ 
		prefix2 = "MC";
	} else {	
		prefix2 = "Data";
	}
	
	TFile* fileppDalitz						= new TFile(fileNameppDalitz.Data());
	
	if( !fileppDalitz ){
	  
	  cout<<"ERROR: The file "<<fileNameppDalitz.Data()<<" does not exist"<<endl;
	  return -1;
	}
	
	TFile* filepPbDalitz						= new TFile( fileNamepPbDalitz.Data());
	
	if( !filepPbDalitz ){
	  
	  cout<<"ERROR: The file "<<fileNamepPbDalitz.Data()<<" does not exist"<<endl;
	  return -1;
	}
	
	if (suffix.CompareTo("eps")==0){
		widthLinesBoxes				= 1.4;
		widthCommonFit					= 2.;
		widthStatErrBars				= 1.5;
		widthCommonErrors				= 1.1;
		widthCommonSpectrumBoxes			= 0.99;
	} else {
		widthLinesBoxes					= 2.3;
		widthCommonFit					= 5.6;
		widthStatErrBars				= 2.6;
		widthCommonErrors				= 2.;
		widthCommonSpectrumBoxes			= 2.3;
	}
	
	
	TGraphAsymmErrors* dummyGraphSystErr;
	
	TString fileNameNLOPi0   	  	= "ExternalInput/Theory/ALICENLOcalcPi0Vogelsang7Tev.data";
	TString fileNamePi0NLODSS	  	= "ExternalInput/Theory/ALICE-DPT-DETA-PI0-7000.RES";
	TString fileNamePi0NLOBKK	  	= "ExternalInput/Theory/lhc_7000_CTEQ5M_BKK_20.dat";
	TString fileNameTheoryCompilation 	= "ExternalInput/TheoryCompilationPP.root";
	TString fileNameTheoryCompilationPbPb	= "ExternalInputPbPb/Theory/TheoryCompilationPbPb.root";
	
	
	TLatex *labelRatioNLOPi07TeV = new TLatex(0.18,0.75,"#pi^{0}, #sqrt{#it{s}} = 7 TeV");
	SetStyleTLatex( labelRatioNLOPi07TeV, 0.17,4);
	
	TLatex *labelRatioNLOPi02760GeV = new TLatex(0.18,0.75,"#pi^{0}, #sqrt{#it{s}} = 2.76 TeV");
	SetStyleTLatex( labelRatioNLOPi02760GeV, 0.13,4);
	
	
		
	TBox* boxErrorSigmaPi02760GeVRatio = CreateBoxConv(colorCommonSpectrumPi02760GeVBox, 0.25, 1.-(xSection2760GeVErr/(xSection2760GeV*1e3) ), 0.29, 1.+(xSection2760GeVErr/(xSection2760GeV*1e3)));
	TBox* boxErrorSigmaPi07TeVRatio    = CreateBoxConv(colorCommonSpectrumPi07TeVBox, 0.25, 1.-(xSection7TeVErrDown /(xSection7TeV*1e3) ), 0.29, 1.+(xSection7TeVErrUp /(xSection7TeV*1e3)));
	
	
  
	 //////////////Loading Results of 7 TeV//////////////////////////////////////////////////
	
	
	TH1D* histoDalitzNumberOfEvents7TeV 		= (TH1D*) fileppDalitz->Get("histoNumberOfEvents7TeV");
	TH1D* histoDalitzNumberOfEvents2760GeV 		= (TH1D*) fileppDalitz->Get("histoNumberOfEvents2.76TeV");
	TH1D* histoDalitzNumberOfEventspPb5023TeV	= (TH1D*) filepPbDalitz->Get("histoNumberOfEventspPb_5.023TeV0-100%");
	
	TDirectory* directoryPi0Dalitz7TeV      	= (TDirectory*) fileppDalitz->Get("Pi0Dalitz7TeV");
	TDirectory* directoryPi0Dalitz2760GeV		= (TDirectory*) fileppDalitz->Get("Pi0Dalitz2.76TeV");
	TDirectory* directoryPi0DalitzpPb5023GeV	= (TDirectory*) filepPbDalitz->Get("Pi0_pPb_5.023TeV_0-100%");
	
	
	graphDalitzPi0InvYieldSysA7TeV	= (TGraphAsymmErrors*)directoryPi0Dalitz7TeV->Get("Pi0SystErrorA");
	graphDalitzPi0InvYieldSys7TeV	= (TGraphAsymmErrors*)directoryPi0Dalitz7TeV->Get("Pi0SystError");
	graphDalitzPi0InvYieldCompl7TeV	= (TGraphAsymmErrors*)directoryPi0Dalitz7TeV->Get("Pi0ComplError");
	TH1D* histoCorretedYieldPi07TeV				= (TH1D*)directoryPi0Dalitz7TeV->Get("CorrectedYieldPi0");
	graphDalitzPi0InvYieldStat7TeV       = new TGraphAsymmErrors(histoCorretedYieldPi07TeV);
	graphDalitzPi0InvYieldStat7TeV->RemovePoint(0);
	
	
	graphDalitzPi0InvYieldSysA7TeVUnShifted  = (TGraphAsymmErrors*)graphDalitzPi0InvYieldSysA7TeV->Clone("graphDalitzPi0InvYieldSysA7TeVUnShifted");
	graphDalitzPi0InvYieldSys7TeVUnShifted   = (TGraphAsymmErrors*)graphDalitzPi0InvYieldSys7TeV->Clone("graphDalitzPi0InvYieldSys7TeVUnShifted");
	graphDalitzPi0InvYieldStat7TeVUnShifted  = (TGraphAsymmErrors*)graphDalitzPi0InvYieldStat7TeV->Clone("graphDalitzPi0InvYieldStat7TeVUnShifted");
	graphDalitzPi0InvYieldCompl7TeVUnShifted = (TGraphAsymmErrors*)graphDalitzPi0InvYieldCompl7TeV->Clone("graphDalitzPi0InvYieldCompl7TeVUnShifted");
	
	
	
	
	graphDalitzPi0InvCrossSectionStat7TeVUnShifted  		= (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldStat7TeVUnShifted,	xSection7TeV*recalcBarn);
	graphDalitzPi0InvCrossSectionSys7TeVUnShifted   		= (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldSys7TeVUnShifted,	xSection7TeV*recalcBarn);
	graphDalitzPi0InvCrossSectionSysA7TeVUnShifted  		= (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldSysA7TeVUnShifted,	xSection7TeV*recalcBarn);
	graphDalitzPi0InvCrossSectionCompl7TeVUnShifted 		= (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldCompl7TeVUnShifted,	xSection7TeV*recalcBarn);
	       
	
	
		
		
	graphDalitzPi0InvYieldSysA2760GeV	= (TGraphAsymmErrors*)directoryPi0Dalitz2760GeV->Get("Pi0SystErrorA");
	graphDalitzPi0InvYieldSys2760GeV	= (TGraphAsymmErrors*)directoryPi0Dalitz2760GeV->Get("Pi0SystError");
	graphDalitzPi0InvYieldCompl2760GeV	= (TGraphAsymmErrors*)directoryPi0Dalitz2760GeV->Get("Pi0ComplError");
	TH1D* histoCorretedYieldPi02760GeV 	= (TH1D*)directoryPi0Dalitz2760GeV->Get("CorrectedYieldPi0");
	graphDalitzPi0InvYieldStat2760GeV  	= new TGraphAsymmErrors(histoCorretedYieldPi02760GeV);
	graphDalitzPi0InvYieldStat2760GeV->RemovePoint(0);
	
	
	graphDalitzPi0InvYieldStat2760GeVUnShifted  			= (TGraphAsymmErrors*)graphDalitzPi0InvYieldStat2760GeV->Clone("graphDalitzPi0InvYieldStat2760GeVUnShifted");
	graphDalitzPi0InvYieldSysA2760GeVUnShifted  			= (TGraphAsymmErrors*)graphDalitzPi0InvYieldSysA2760GeV->Clone("graphDalitzPi0InvYieldSysA2760GeVUnShifted");
	graphDalitzPi0InvYieldSys2760GeVUnShifted   			= (TGraphAsymmErrors*)graphDalitzPi0InvYieldSys2760GeV->Clone("graphDalitzPi0InvYieldSys2760GeVUnShifted");
	graphDalitzPi0InvYieldCompl2760GeVUnShifted 			= (TGraphAsymmErrors*)graphDalitzPi0InvYieldCompl2760GeV->Clone("graphDalitzPi0InvYieldCompl2760GeVUnShifted");
	
	graphDalitzPi0InvCrossSectionStat2760GeVUnShifted  		= (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldStat2760GeVUnShifted,	xSection2760GeV*recalcBarn);
	graphDalitzPi0InvCrossSectionSys2760GeVUnShifted   		= (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldSys2760GeVUnShifted,	xSection2760GeV*recalcBarn);
	graphDalitzPi0InvCrossSectionSysA2760GeVUnShifted  		= (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldSysA2760GeVUnShifted,	xSection2760GeV*recalcBarn);
	graphDalitzPi0InvCrossSectionCompl2760GeVUnShifted 		= (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldCompl2760GeVUnShifted,	xSection2760GeV*recalcBarn);
	       
	
	TGraphAsymmErrors* graphDalitzPi0InvYieldSyspPb5023GeV   	= (TGraphAsymmErrors*)directoryPi0DalitzpPb5023GeV->Get("Pi0SystError");
	TGraphAsymmErrors* graphDaliizPi0InvYieldSysApPb5023GeV		= (TGraphAsymmErrors*)directoryPi0DalitzpPb5023GeV->Get("Pi0SystErrorA");
	TGraphAsymmErrors* graphDalitzPi0InvYieldComplpPb5023GeV	= (TGraphAsymmErrors*)directoryPi0DalitzpPb5023GeV->Get("Pi0ComplError");
	TH1D* histoCorretedYieldPi05023GeV				= (TH1D*)directoryPi0DalitzpPb5023GeV->Get("CorrectedYieldPi0");
	TGraphAsymmErrors* graphDalitzPi0InvYieldStatpPb5023GeV      	= new TGraphAsymmErrors(histoCorretedYieldPi05023GeV);
	graphDalitzPi0InvYieldStatpPb5023GeV->RemovePoint(0);
	
	
		
	if( bWCorrection.CompareTo("X") == 0 ){
	  
	  
	  
	       //------------------------------Invariant yield pi0-----------------------------------------------
	  
	       Double_t paramFitBinShiftXPP7TeV[3]	= {2.1,6.3,0.14};
	       
	       
	       TF1* fitTsallisPi07TeVPtMult = FitObject("tmpt","tmptInvYieldPP7TeV","Pi0");
	       
	       
	       fitTsallisPi07TeVPtMult->SetParameters(paramFitBinShiftXPP7TeV[0],paramFitBinShiftXPP7TeV[1], paramFitBinShiftXPP7TeV[2]) ; // standard parameter optimize if necessary
	       
	       
	       graphDalitzPi0InvYieldCompl7TeV	= ApplyXshift( graphDalitzPi0InvYieldCompl7TeV, fitTsallisPi07TeVPtMult,"Pi0");
	       graphDalitzPi0InvYieldSys7TeV    = ApplyXshiftIndividualSpectra (graphDalitzPi0InvYieldCompl7TeV, graphDalitzPi0InvYieldSys7TeV, fitTsallisPi07TeVPtMult,  0, graphDalitzPi0InvYieldSys7TeV->GetN());
	       graphDalitzPi0InvYieldStat7TeV   = ApplyXshiftIndividualSpectra (graphDalitzPi0InvYieldCompl7TeV, graphDalitzPi0InvYieldStat7TeV, fitTsallisPi07TeVPtMult, 0, graphDalitzPi0InvYieldStat7TeV->GetN());
	       graphDalitzPi0InvYieldSysA7TeV   = ApplyXshiftIndividualSpectra (graphDalitzPi0InvYieldCompl7TeV, graphDalitzPi0InvYieldSysA7TeV, fitTsallisPi07TeVPtMult, 0, graphDalitzPi0InvYieldSysA7TeV->GetN());
	       
	       
	       Double_t paramFitBinShiftXPP2760GeV[3]  = {1.5,8.0,0.13};
	       TF1* fitTsallisPi02760GeVPtMult = FitObject("tmpt","tmptInvYieldPP2760GeV","Pi0");
	       fitTsallisPi02760GeVPtMult->SetParameters(paramFitBinShiftXPP2760GeV[0],paramFitBinShiftXPP2760GeV[1], paramFitBinShiftXPP2760GeV[2]) ; // standard parameter optimize if necessary
	       	
		
	       graphDalitzPi0InvYieldCompl2760GeV    = ApplyXshift( graphDalitzPi0InvYieldCompl2760GeV, fitTsallisPi02760GeVPtMult,"Pi0");
	       graphDalitzPi0InvYieldSys2760GeV      = ApplyXshiftIndividualSpectra (graphDalitzPi0InvYieldCompl2760GeV, graphDalitzPi0InvYieldSys2760GeV, fitTsallisPi02760GeVPtMult,  0, graphDalitzPi0InvYieldSys2760GeV->GetN());
	       graphDalitzPi0InvYieldStat2760GeV     = ApplyXshiftIndividualSpectra (graphDalitzPi0InvYieldCompl2760GeV, graphDalitzPi0InvYieldStat2760GeV, fitTsallisPi02760GeVPtMult, 0, graphDalitzPi0InvYieldStat2760GeV->GetN());
	       graphDalitzPi0InvYieldSysA2760GeV     = ApplyXshiftIndividualSpectra (graphDalitzPi0InvYieldCompl2760GeV, graphDalitzPi0InvYieldSysA2760GeV, fitTsallisPi02760GeVPtMult, 0, graphDalitzPi0InvYieldSysA2760GeV->GetN());
	       
	     
	       
	       Double_t paramFitBinShiftXpPb5023GeV[3]  = {7.0,8.0,0.16};
	       TF1* fitTsallisPi0pPb5023GeVPtMult = FitObject("tmpt","tmpt","Pi0");
	       fitTsallisPi0pPb5023GeVPtMult->SetParameters(paramFitBinShiftXpPb5023GeV[0],paramFitBinShiftXpPb5023GeV[1], paramFitBinShiftXpPb5023GeV[2]) ; // standard parameter optimize if necessary
	              
	       
	       graphDalitzPi0InvYieldComplpPb5023GeV    = ApplyXshift( graphDalitzPi0InvYieldComplpPb5023GeV, fitTsallisPi0pPb5023GeVPtMult,"Pi0");
	       graphDalitzPi0InvYieldSyspPb5023GeV   	= ApplyXshiftIndividualSpectra (graphDalitzPi0InvYieldComplpPb5023GeV, graphDalitzPi0InvYieldSyspPb5023GeV, fitTsallisPi0pPb5023GeVPtMult,  0, graphDalitzPi0InvYieldSyspPb5023GeV->GetN());
	       graphDalitzPi0InvYieldStatpPb5023GeV  	= ApplyXshiftIndividualSpectra (graphDalitzPi0InvYieldComplpPb5023GeV, graphDalitzPi0InvYieldStatpPb5023GeV, fitTsallisPi0pPb5023GeVPtMult, 0, graphDalitzPi0InvYieldStatpPb5023GeV->GetN());
	       
	       
	       
	       graphDalitzPi0InvCrossSectionStat7TeV  		= (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldStat7TeV,	xSection7TeV*recalcBarn);
	       graphDalitzPi0InvCrossSectionSys7TeV   		= (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldSys7TeV,		xSection7TeV*recalcBarn);
	       graphDalitzPi0InvCrossSectionSysA7TeV  		= (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldSysA7TeV,	xSection7TeV*recalcBarn);
	       graphDalitzPi0InvCrossSectionCompl7TeV 		= (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldCompl7TeV,	xSection7TeV*recalcBarn);
	       
	       graphDalitzPi0InvCrossSectionStat2760GeV    	= (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldStat2760GeV,	xSection2760GeV*recalcBarn);
	       graphDalitzPi0InvCrossSectionSys2760GeV      	= (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldSys2760GeV, 	xSection2760GeV*recalcBarn);
	       graphDalitzPi0InvCrossSectionSysA2760GeV		= (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldSysA2760GeV,	xSection2760GeV*recalcBarn);
	       graphDalitzPi0InvCrossSectionCompl2760GeV	= (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldCompl2760GeV,	xSection2760GeV*recalcBarn);
	   
	  
	}
	
	
	
	
	
	
	
	TFile* fileTheoryCompilation = new TFile(fileNameTheoryCompilation.Data());
	
	if( ! fileTheoryCompilation ) {
	  
	  cout<<"ERROR: Please the file "<<fileNameTheoryCompilation.Data()<<" was not found"<<endl;
	  
	}
	
	TH1F* histoPythia8InvXSection                       		= (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Spec2760GeV");
	
	if( !histoPythia8InvXSection ) return 0;
	
	TGraph* graphNLOCalcMuHalf2760GeV				=    	(TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuHalf2760GeV");
	TGraph* graphNLOCalcMuOne2760GeV				=     	(TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuOne2760GeV");
	TGraph* graphNLOCalcMuTwo2760GeV				=     	(TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuTwo2760GeV");
	
	cout<<"Pedrito NLOCal"<<endl;
	TGraph* graphNLOCalcMuHalf7TeV					=    	(TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuHalf7000GeV");
	TGraph* graphNLOCalcMuOne7TeV					=     	(TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuOne7000GeV");
	TGraph* graphNLOCalcMuTwo7TeV					=     	(TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuTwo7000GeV");
	TGraph* graphNLOBKKCalcMuTwo7TeV				=  	(TGraph*)fileTheoryCompilation->Get("graphNLOCalcBKKInvSecPi0MuTwo7000GeV");
	TGraphAsymmErrors* graphNLODSS14CalcMuOne7TeV			= 	(TGraphAsymmErrors*)fileTheoryCompilation->Get("graphNLOCalcDSS14InvCrossSec7000GeV");
	
	for (Int_t   iPoint = 0; iPoint < 30; iPoint++) {
	graphNLODSS14CalcMuOne7TeV->RemovePoint(graphNLODSS14CalcMuOne7TeV->GetN()-1);
	}
	graphNLODSS14CalcMuOne7TeV->RemovePoint(0);
	
	graphNLODSS14CalcMuOne7TeV->Print();
	
	
	
	TFile* fileTheoryCompilationPbPb = new TFile(fileNameTheoryCompilationPbPb.Data());
	
	if( ! fileTheoryCompilationPbPb ) {
	  
	  cout<<"ERROR: Please the file "<<fileNameTheoryCompilationPbPb.Data()<<" was not found"<<endl;
	  
	}
	
	
      
	TGraphAsymmErrors*	graphNLOCalcmuHalfDSS14InvSecPi02760GeV =	(TGraphAsymmErrors*)fileTheoryCompilationPbPb->Get("graphNLOCalcmuHalfDSS14InvSecPi02760GeV");
	TGraphAsymmErrors* 	graphNLOCalcmuOneDSS14InvSecPi02760GeV	=	(TGraphAsymmErrors*)fileTheoryCompilationPbPb->Get("graphNLOCalcDSS14InvCrossSecPi02760GeV");
	TGraphAsymmErrors*	graphNLOCalcmuTwoDSS14InvSecPi02760GeV	=	(TGraphAsymmErrors*)fileTheoryCompilationPbPb->Get("graphNLOCalcmuTwoDSS14InvSecPi02760GeV");
	
	for (Int_t   iPoint = 0; iPoint < 30; iPoint++) {
	graphNLOCalcmuHalfDSS14InvSecPi02760GeV->RemovePoint(graphNLOCalcmuHalfDSS14InvSecPi02760GeV->GetN()-1);
	graphNLOCalcmuOneDSS14InvSecPi02760GeV->RemovePoint(graphNLOCalcmuOneDSS14InvSecPi02760GeV->GetN()-1);
	graphNLOCalcmuTwoDSS14InvSecPi02760GeV->RemovePoint(graphNLOCalcmuTwoDSS14InvSecPi02760GeV->GetN()-1);
	}
	graphNLOCalcmuHalfDSS14InvSecPi02760GeV->RemovePoint(0);
	graphNLOCalcmuOneDSS14InvSecPi02760GeV->RemovePoint(0);
	graphNLOCalcmuTwoDSS14InvSecPi02760GeV->RemovePoint(0);
	
	graphNLOCalcmuHalfDSS14InvSecPi02760GeV->Print();
	
	
	
	TString  nameFinalResDat = Form("%s/FitResultsData.dat",outputDir.Data());
        TString parametersOutput;
	fstream  fileFitFinalResults;
	fileFitFinalResults.open(nameFinalResDat.Data(), ios::out);
			  
	
	
	Double_t parametersTsallis2760GeV[3] 		= {2.,5.,0.18};//{2.18,6.88,0.139};  //2.,5.,0.18
	Double_t parametersTsallisInvCross2760GeV[3]    = {1.0e12,6.88,0.139};
		
	Double_t parametersTsallis7TeV[3]    		= {2.,5.,0.18}; //{2.3,6.3,0.14};
	Double_t parametersTsallisInvCross7TeV[3]	= {1.0e12,6.3,0.14};
	
	Double_t parametersTsalisspPb5023GeV[3]		= {7.0,8.0,0.19};
	
	Double_t parametersBylinkin2760GeV[5] 		= {0.06,0.3,6,0.3,3.8};
	Double_t parametersBylinkinInvCross2760GeV[5] 	= {graphDalitzPi0InvCrossSectionCompl2760GeV->GetY()[0],0.3,graphDalitzPi0InvCrossSectionCompl2760GeV->GetY()[0],0.3,8};
	
	
	Double_t parametersBylinkin7TeV[5]		= {0.08,0.4,.6,0.3,1.6};
	Double_t parametersBylinkinInvCross7TeV[5]      = {graphDalitzPi0InvCrossSectionCompl7TeV->GetY()[0],0.4,graphDalitzPi0InvCrossSectionCompl7TeV->GetY()[0],0.3,1.6};
	Double_t parametersBylinkinpPb5023GeV[5]	= {10,0.3,1,0.3,8.0};
	
	// Tsallis fits  
	
	TF1 *fitTsallisPi02760GeV = FitObject("l","fitTsallisPi02760GeV","Pi0");
	fitTsallisPi02760GeV->SetRange(0.6, 10.0);
	fitTsallisPi02760GeV->SetParameters(parametersTsallis2760GeV[0],parametersTsallis2760GeV[1],parametersTsallis2760GeV[2]);
	graphDalitzPi0InvYieldCompl2760GeV->Fit(fitTsallisPi02760GeV,"QNRMEX0+","",0.6,10.0);
	//graphDalitzPi0InvYieldCompl2760GeV->Fit(fitTsallisPi02760GeV,"QNRMEX0+","",0.6,10.0);
	cout<<"Fit prob Tsallis 2760 GeV:  "<<fitTsallisPi02760GeV->GetProb() <<endl;
	parametersOutput= WriteParameterToFile(fitTsallisPi02760GeV);
	
	fileFitFinalResults << parametersOutput.Data();
	TH1F* fitHistogramTsallisPi02760GeV = (TH1F*)fitTsallisPi02760GeV->GetHistogram();
	
	
	TF1 *fitTsallisInvCrossSectionPi02760GeV = FitObject("l","fitTsallisInvCrossSectionPi02760GeV","Pi0");
	fitTsallisInvCrossSectionPi02760GeV->SetRange(minPtPi0, maxPtPi0);
	fitTsallisInvCrossSectionPi02760GeV->SetParameters(parametersTsallisInvCross2760GeV[0],parametersTsallisInvCross2760GeV[1],parametersTsallisInvCross2760GeV[2]);
	
	graphDalitzPi0InvCrossSectionCompl2760GeV->Fit(fitTsallisInvCrossSectionPi02760GeV,"QNRMEX0+","",minPtPi0,maxPtPi0);
	//graphDalitzPi0InvCrossSectionCompl2760GeV->Fit(fitTsallisInvCrossSectionPi02760GeV,"QNRMEX0+","",minPtPi0,maxPtPi0);
	
	parametersOutput= WriteParameterToFile(fitTsallisInvCrossSectionPi02760GeV);
	fileFitFinalResults << parametersOutput.Data();
	TH1F* fitHistogramTsallisInvCrossSectionPi02760GeV = (TH1F*)fitTsallisInvCrossSectionPi02760GeV->GetHistogram();
	
	
	
	
	TF1 *fitTsallisPi07TeV = FitObject("l","fitTsallisPi07TeV","Pi0");
	fitTsallisPi07TeV->SetRange(0.6,10.0);
	fitTsallisPi07TeV->SetParameters(parametersTsallis7TeV[0],parametersTsallis7TeV[1],parametersTsallis7TeV[2]);
	graphDalitzPi0InvYieldCompl7TeV->Fit(fitTsallisPi07TeV,"QNRMEX0+","",0.6,10.0);
	graphDalitzPi0InvYieldCompl7TeV->Fit(fitTsallisPi07TeV,"QNRMEX0+","",0.6,10.0);
	
	cout<<"Fit prob Tsallis 7 TeV:  "<<fitTsallisPi07TeV->GetProb() <<endl;
	parametersOutput= WriteParameterToFile(fitTsallisPi07TeV);
	fileFitFinalResults << parametersOutput.Data();
	
	TH1F* fitHistogramTsallisPi07TeV = (TH1F*)fitTsallisPi07TeV->GetHistogram();
	
	
	
	TF1 *fitTsallisInvCrossSectionPi07TeV = FitObject("l","fitTsallisInvCrossSectionPi07TeV","Pi0");
	fitTsallisInvCrossSectionPi07TeV->SetRange(minPtPi0, maxPtPi0);
	fitTsallisInvCrossSectionPi07TeV->SetParameters(parametersTsallisInvCross7TeV[0],parametersTsallisInvCross7TeV[1],parametersTsallisInvCross7TeV[2]);
	
	graphDalitzPi0InvCrossSectionCompl7TeV->Fit(fitTsallisInvCrossSectionPi07TeV,"QNRMEX0+","",minPtPi0,maxPtPi0);
	graphDalitzPi0InvCrossSectionCompl7TeV->Fit(fitTsallisInvCrossSectionPi07TeV,"QNRMEX0+","",minPtPi0,maxPtPi0);
	
	parametersOutput= WriteParameterToFile(fitTsallisInvCrossSectionPi07TeV);
	fileFitFinalResults << parametersOutput.Data();
	TH1F* fitHistogramTsallisInvCrossSectionPi07TeV = (TH1F*)fitTsallisInvCrossSectionPi07TeV->GetHistogram();
	
	
	
	
	
	TF1 *fitTsallisPi0pPb5023GeV = FitObject("l","fitTsallisPi0pPb5023GeV","Pi0");
	fitTsallisPi0pPb5023GeV->SetRange(0.6,10.0);
	fitTsallisPi0pPb5023GeV->SetParameters(parametersTsalisspPb5023GeV);
	graphDalitzPi0InvYieldComplpPb5023GeV->Fit(fitTsallisPi0pPb5023GeV,"QNRMEX0+","",0.6,10.0);
	graphDalitzPi0InvYieldComplpPb5023GeV->Fit(fitTsallisPi0pPb5023GeV,"QNRMEX0+","",0.6,10.0);
	
	cout<<"Fit prob Tsallis pPb 5023 TeV:  "<<fitTsallisPi0pPb5023GeV->GetProb() <<endl;
	parametersOutput= WriteParameterToFile(fitTsallisPi0pPb5023GeV);
	fileFitFinalResults << parametersOutput.Data();
	
	TH1F* fitHistogramTsallisPi0pPb5023GeV = (TH1F*)fitTsallisPi0pPb5023GeV->GetHistogram();
	TH1F* fitHistogramTsallisPi0pPb5023GeVCopy  = (TH1F*) fitHistogramTsallisPi0pPb5023GeV->Clone();
	
	
	
	TGraphAsymmErrors* graphRatioToTsallisFitPi0PP2760GeV 	= (TGraphAsymmErrors*) CalculateGraphErrRatioToFit(graphDalitzPi0InvYieldCompl2760GeV,fitTsallisPi02760GeV);
	TGraphAsymmErrors* graphRatioToTsallisFitPi0SysPP2760GeV 	= (TGraphAsymmErrors*) CalculateGraphErrRatioToFit(graphDalitzPi0InvYieldSys2760GeV,fitTsallisPi02760GeV);
	TGraphAsymmErrors* graphRatioToTsallisFitPi0StatPP2760GeV 	= (TGraphAsymmErrors*) CalculateGraphErrRatioToFit(graphDalitzPi0InvYieldStat2760GeV,fitTsallisPi02760GeV);
	
	TGraphAsymmErrors* graphRatioToTsallisFitPi0PP7TeV 	= (TGraphAsymmErrors*) CalculateGraphErrRatioToFit(graphDalitzPi0InvYieldCompl7TeV,fitTsallisPi07TeV);
	TGraphAsymmErrors* graphRatioToTsallisFitPi0SysPP7TeV 	= (TGraphAsymmErrors*) CalculateGraphErrRatioToFit(graphDalitzPi0InvYieldSys7TeV,fitTsallisPi07TeV);
	TGraphAsymmErrors* graphRatioToTsallisFitPi0StatPP7TeV 	= (TGraphAsymmErrors*) CalculateGraphErrRatioToFit(graphDalitzPi0InvYieldStat7TeV,fitTsallisPi07TeV);
	
	TGraphAsymmErrors* graphRatioToTsallisFitPi0pPb5023TeV  = (TGraphAsymmErrors*) CalculateGraphErrRatioToFit(graphDalitzPi0InvYieldComplpPb5023GeV,fitTsallisPi0pPb5023GeV);
	TGraphAsymmErrors* graphRatioToTsallisFitPi0SyspPb5023TeV  = (TGraphAsymmErrors*) CalculateGraphErrRatioToFit(graphDalitzPi0InvYieldSyspPb5023GeV,fitTsallisPi0pPb5023GeV);
	TGraphAsymmErrors* graphRatioToTsallisFitPi0StatpPb5023TeV  = (TGraphAsymmErrors*) CalculateGraphErrRatioToFit(graphDalitzPi0InvYieldStatpPb5023GeV,fitTsallisPi0pPb5023GeV);
	
	
	
	
	
	// Bylinkin fits
	
	
	TF1 *fitBylinkinPi02760GeV = FitObject("tcm","fitBylinkinPi02760GeV","Pi0");
	fitBylinkinPi02760GeV->SetRange(0.6, 10.0);
	fitBylinkinPi02760GeV->SetParameters(parametersBylinkin2760GeV[0],parametersBylinkin2760GeV[1],parametersBylinkin2760GeV[2],parametersBylinkin2760GeV[3],parametersBylinkin2760GeV[4]);
	graphDalitzPi0InvYieldCompl2760GeV->Fit(fitBylinkinPi02760GeV,"QNRME+","",0.6,10.0);
	

	TH1F* fitHistogramBylinkinPi02760GeV = (TH1F*)fitBylinkinPi02760GeV->GetHistogram();
	
	parametersOutput= WriteParameterToFile(fitBylinkinPi02760GeV);
	fileFitFinalResults << parametersOutput.Data();
	
	
	TF1 *fitBylinkinDalitzPi0InvCrossSection2760GeV = FitObject("tcm","fitBylinkinDalitzPi0InvCrossSection2760GeV","Pi0");
	fitBylinkinDalitzPi0InvCrossSection2760GeV->SetRange(0.6, 10.0);
	fitBylinkinDalitzPi0InvCrossSection2760GeV->SetParameters(parametersBylinkinInvCross2760GeV[0],parametersBylinkinInvCross2760GeV[1],parametersBylinkinInvCross2760GeV[2],parametersBylinkinInvCross2760GeV[3],parametersBylinkinInvCross2760GeV[4]);
	graphDalitzPi0InvCrossSectionCompl2760GeV->Fit(fitBylinkinDalitzPi0InvCrossSection2760GeV,"QNRME+","",0.6,10.0);
	
	TH1F* fitHistogramBylinkinInvCrossSection2760GeV = (TH1F*)fitBylinkinDalitzPi0InvCrossSection2760GeV->GetHistogram();
	
	parametersOutput= WriteParameterToFile(fitBylinkinDalitzPi0InvCrossSection2760GeV);
	fileFitFinalResults << parametersOutput.Data();
	
	
	
	
	TF1 *fitBylinkinPi07TeV = FitObject("tcm","fitBylinkinPi07TeV","Pi0");
	fitBylinkinPi07TeV->SetRange(0.6, 10.0);
	fitBylinkinPi07TeV->SetParameters(parametersBylinkin7TeV[0],parametersBylinkin7TeV[1],parametersBylinkin7TeV[2],parametersBylinkin7TeV[3],parametersBylinkin7TeV[4]);
	graphDalitzPi0InvYieldCompl7TeV->Fit(fitBylinkinPi07TeV,"QNRME+","",0.6,10.0);
	graphDalitzPi0InvYieldCompl7TeV->Fit(fitBylinkinPi07TeV,"QNRME+","",0.6,10.0);
	//graphDalitzPi0InvYieldCompl7TeV->Fit(fitBylinkinPi07TeV,"SQNRME+","",0.6,10.0);
	
	parametersOutput= WriteParameterToFile(fitBylinkinPi07TeV);
	fileFitFinalResults << parametersOutput.Data();
	
	TH1F* fitHistogramBylinkinPi07TeV = (TH1F*)fitBylinkinPi07TeV->GetHistogram();
	
	
	TF1 *fitBylinkinDalitzPi0InvCrossSection7TeV = FitObject("tcm","fitBylinkinDalitzPi0InvCrossSection7TeV","Pi0");
	fitBylinkinDalitzPi0InvCrossSection7TeV->SetRange(0.6, 10.0);
	fitBylinkinDalitzPi0InvCrossSection7TeV->SetParameters(parametersBylinkinInvCross7TeV[0],parametersBylinkinInvCross7TeV[1],parametersBylinkinInvCross7TeV[2],parametersBylinkinInvCross7TeV[3],parametersBylinkinInvCross7TeV[4]);
	graphDalitzPi0InvCrossSectionCompl7TeV->Fit(fitBylinkinDalitzPi0InvCrossSection7TeV,"QNRME+","",0.6,10.0);
	
	

	TH1F* fitHistogramBylinkinDalitzPi0InvCrossSection7TeV = (TH1F*)fitBylinkinDalitzPi0InvCrossSection7TeV->GetHistogram();
	
	parametersOutput= WriteParameterToFile(fitBylinkinDalitzPi0InvCrossSection7TeV);
	fileFitFinalResults << parametersOutput.Data();
	
	
	
	TF1 *fitBylinkinPi0pPb5023GeV = FitObject("tcm","fitBylinkinPi0pPb5023GeV","Pi0");
	fitBylinkinPi0pPb5023GeV->SetRange(0.6, 10.0);
	fitBylinkinPi0pPb5023GeV->SetParameters(parametersBylinkinpPb5023GeV[0],parametersBylinkinpPb5023GeV[1],parametersBylinkinpPb5023GeV[2],parametersBylinkinpPb5023GeV[3],parametersBylinkinpPb5023GeV[4]);
	
	graphDalitzPi0InvYieldComplpPb5023GeV->Fit(fitBylinkinPi0pPb5023GeV,"QNRME+","",0.6,10.0);
	graphDalitzPi0InvYieldComplpPb5023GeV->Fit(fitBylinkinPi0pPb5023GeV,"SQNRME+","",0.6,10.0);
	graphDalitzPi0InvYieldComplpPb5023GeV->Fit(fitBylinkinPi0pPb5023GeV,"SQNRME+","",0.6,10.0);
	
	
	parametersOutput= WriteParameterToFile(fitBylinkinPi0pPb5023GeV);
	fileFitFinalResults << parametersOutput.Data();
	
	TH1F* fitHistogramBylinkinPi0pPb5023GeV = (TH1F*)fitBylinkinPi0pPb5023GeV->GetHistogram();
	
	
	
	
	TGraphAsymmErrors* graphRatioToBylinkinFitPi0PP2760GeV 	= (TGraphAsymmErrors*) CalculateGraphErrRatioToFit(graphDalitzPi0InvYieldCompl2760GeV,fitBylinkinPi02760GeV);
	TGraphAsymmErrors* graphRatioToBylinkinFitPi0SysPP2760GeV 	= (TGraphAsymmErrors*) CalculateGraphErrRatioToFit(graphDalitzPi0InvYieldSys2760GeV,fitBylinkinPi02760GeV);
	TGraphAsymmErrors* graphRatioToBylinkinFitPi0StatPP2760GeV 	= (TGraphAsymmErrors*) CalculateGraphErrRatioToFit(graphDalitzPi0InvYieldStat2760GeV,fitBylinkinPi02760GeV);
	
	TGraphAsymmErrors* graphRatioToBylinkinFitPi0PP7TeV 	= (TGraphAsymmErrors*) CalculateGraphErrRatioToFit(graphDalitzPi0InvYieldCompl7TeV,fitBylinkinPi07TeV);
	TGraphAsymmErrors* graphRatioToBylinkinFitPi0SysPP7TeV 	= (TGraphAsymmErrors*) CalculateGraphErrRatioToFit(graphDalitzPi0InvYieldSys7TeV,fitBylinkinPi07TeV);
	TGraphAsymmErrors* graphRatioToBylinkinFitPi0StatPP7TeV 	= (TGraphAsymmErrors*) CalculateGraphErrRatioToFit(graphDalitzPi0InvYieldStat7TeV,fitBylinkinPi07TeV);
	
	TGraphAsymmErrors* graphRatioToBylinkinFitPi0pPb5023TeV  = (TGraphAsymmErrors*) CalculateGraphErrRatioToFit(graphDalitzPi0InvYieldComplpPb5023GeV,fitBylinkinPi0pPb5023GeV);
	TGraphAsymmErrors* graphRatioToBylinkinFitPi0SyspPb5023TeV  = (TGraphAsymmErrors*) CalculateGraphErrRatioToFit(graphDalitzPi0InvYieldSyspPb5023GeV,fitBylinkinPi0pPb5023GeV);
	TGraphAsymmErrors* graphRatioToBylinkinFitPi0StatpPb5023TeV  = (TGraphAsymmErrors*) CalculateGraphErrRatioToFit(graphDalitzPi0InvYieldStatpPb5023GeV,fitBylinkinPi0pPb5023GeV);
	
	
	
	
	fitHistogramTsallisPi02760GeV->Scale(1e-1);
	fitHistogramBylinkinPi02760GeV->Scale(1e-1);
	
	fitHistogramTsallisPi07TeV->Scale(1e0);
	fitHistogramBylinkinPi07TeV->Scale(1e0);
	
	fitHistogramTsallisPi0pPb5023GeV->Scale(1e1);
	fitHistogramBylinkinPi0pPb5023GeV->Scale(1e1);
	
		
	TGraphAsymmErrors* graphDalitzPi0InvYieldCompl2760GeVScale 	= (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldCompl2760GeV,1e-1);
	TGraphAsymmErrors* graphDalitzPi0InvYieldCompl7TeVScale 	= (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldCompl7TeV,1e-0);
	
	
	TGraphAsymmErrors* graphDalitzPi0InvYieldSys2760GeVScale 	= (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldSys2760GeV,1e-1);
	TGraphAsymmErrors* graphDalitzPi0InvYieldStat2760GeVScale	= (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldStat2760GeV,1e-1);
	TGraphAsymmErrors* graphDalitzPi0InvYieldSys7TeVScale    	= (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldSys7TeV,1e0);
	TGraphAsymmErrors* graphDalitzPi0InvYieldStat7TeVScale   	= (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldStat7TeV,1e0);
	TGraphAsymmErrors* graphDalitzPi0InvYieldSyspPb5023GeVScale     = (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldSyspPb5023GeV,1e1);
	TGraphAsymmErrors* graphDalitzPi0InvYieldStatpPb5023GeVScale    = (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvYieldStatpPb5023GeV,1e1);
	
	
	////////////////////////////////////////Setting histograms//////////////////////////////////
	
	//void SetStyleTLatex( TLatex* text, Size_t textSize, Width_t lineWidth, Color_t textColor = 1, Font_t textFont = 42, Bool_t kNDC = kTRUE){

	TLatex *labelRatioFitsPi0pPb5023GeV = new TLatex(0.18,0.75,"#pi^{0}, #sqrt{#it{s_{NN}}} = 5.02 TeV");
	SetStyleTLatex( labelRatioFitsPi0pPb5023GeV, 36,4,1,43);
	
	
	TLatex *labelRatioFitsPi07TeV = new TLatex(0.18,0.75,"#pi^{0}, #sqrt{#it{s}} = 7 TeV");
	SetStyleTLatex( labelRatioFitsPi07TeV,36,4,1,43);
	
	TLatex *labelRatioFitsPi02760GeV = new TLatex(0.18,0.75,"#pi^{0}, #sqrt{#it{s}} = 2.76 TeV");
	SetStyleTLatex( labelRatioFitsPi02760GeV,36,4,1,43);
	
	DrawGammaSetMarkerTGraphAsym(graphDalitzPi0InvYieldStat2760GeVScale, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kFALSE);
	graphDalitzPi0InvYieldStat2760GeVScale->SetLineWidth(widthCommonErrors);
	DrawGammaSetMarkerTGraphAsym(graphDalitzPi0InvYieldSys2760GeVScale, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kTRUE, colorCommonSpectrumPi02760GeVBox);
	graphDalitzPi0InvYieldSys2760GeVScale->SetLineWidth(0);
	
	
	DrawGammaSetMarkerTGraphAsym(graphDalitzPi0InvYieldStat7TeV, markerStyleCommmonSpectrumPi07TeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kFALSE);
	graphDalitzPi0InvYieldStat7TeV->SetLineWidth(widthCommonErrors);	
	DrawGammaSetMarkerTGraphAsym(graphDalitzPi0InvYieldSys7TeV, markerStyleCommmonSpectrumPi07TeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kTRUE, colorCommonSpectrumPi07TeVBox);
	graphDalitzPi0InvYieldSys7TeV->SetLineWidth(widthCommonErrors);	
	
     
	
	DrawGammaSetMarkerTGraphAsym(graphDalitzPi0InvYieldStatpPb5023GeVScale, markerStyleCommmonSpectrumPi0pPb5023GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi0pPb5023GeV, colorCommonSpectrumPi0pPb5023GeV, widthCommonSpectrumBoxes, kFALSE);
	graphDalitzPi0InvYieldStatpPb5023GeVScale->SetLineWidth(widthCommonErrors);
	DrawGammaSetMarkerTGraphAsym(graphDalitzPi0InvYieldSyspPb5023GeVScale, markerStyleCommmonSpectrumPi0pPb5023GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi0pPb5023GeV, colorCommonSpectrumPi0pPb5023GeV, widthCommonSpectrumBoxes, kTRUE, colorCommonSpectrumPi0pPb5023GeVBox);
	graphDalitzPi0InvYieldSyspPb5023GeVScale->SetLineWidth(0);
	
	
	TH1F*  fitTsallisLineStyle 	=  (TH1F*)fitHistogramTsallisPi02760GeV->Clone("fitTsallisLineStyle");
	TH1F*  fitBylinkinLineStyle 	=  (TH1F*)fitHistogramBylinkinPi02760GeV->Clone("fitBylinkinLineStyle");
	
	SetStyleHisto(fitTsallisLineStyle, 			widthCommonFit, styleFitCommonSpectrum, 		kBlack);
	SetStyleHisto(fitBylinkinLineStyle,			widthCommonFit, styleBylinkinFitCommonSpectrum, 	kBlack);
	
	
	SetStyleHisto(fitHistogramTsallisPi02760GeV, 		widthCommonFit, styleFitCommonSpectrum, 		colorCommonSpectrumPi02760GeV);
	SetStyleHisto(fitHistogramBylinkinPi02760GeV,		widthCommonFit, styleBylinkinFitCommonSpectrum, 	colorCommonSpectrumPi02760GeV);
	
	SetStyleHisto(fitHistogramTsallisPi07TeV, 		widthCommonFit, styleFitCommonSpectrum, 		colorCommonSpectrumPi07TeV);
	SetStyleHisto(fitHistogramBylinkinPi07TeV, 		widthCommonFit, styleBylinkinFitCommonSpectrum, 	colorCommonSpectrumPi07TeV);
	
	
	SetStyleHisto(fitHistogramTsallisPi0pPb5023GeV, 	widthCommonFit, styleFitCommonSpectrum, 		colorCommonSpectrumPi0pPb5023GeV);
	SetStyleHisto(fitHistogramBylinkinPi0pPb5023GeV, 	widthCommonFit, styleBylinkinFitCommonSpectrum, 	colorCommonSpectrumPi0pPb5023GeV);
	
	
	
	DrawGammaSetMarkerTGraphAsym(graphRatioToTsallisFitPi0pPb5023TeV, markerStyleTsallisFitRatioPi0pPb5023GeV, markerSizeCommonSpectrum*1.6, colorCommonSpectrumPi0pPb5023GeV, colorCommonSpectrumPi0pPb5023GeV, widthCommonSpectrumBoxes, kTRUE,colorCommonSpectrumPi0pPb5023GeVBox);
	graphRatioToTsallisFitPi0pPb5023TeV->SetLineWidth(widthCommonErrors);
	
	DrawGammaSetMarkerTGraphAsym(graphRatioToBylinkinFitPi0pPb5023TeV, markerStyleBylinkinFitRatioPi0pPb5023GeV, markerSizeCommonSpectrum*1.6, colorCommonSpectrumPi0pPb5023GeV, colorCommonSpectrumPi0pPb5023GeV, widthCommonSpectrumBoxes, kTRUE);
	graphRatioToBylinkinFitPi0pPb5023TeV->SetLineWidth(widthCommonErrors);
	
	
	
	DrawGammaSetMarkerTGraphAsym(graphRatioToTsallisFitPi0PP7TeV, markerStyleTsallisFitRatioPi07TeV, markerSizeCommonSpectrum*1.6, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kTRUE,colorCommonSpectrumPi07TeVBox);
	graphRatioToTsallisFitPi0PP7TeV->SetLineWidth(widthCommonErrors);
	
	DrawGammaSetMarkerTGraphAsym(graphRatioToBylinkinFitPi0PP7TeV, markerStyleBylinkinFitRatioPi07TeV, markerSizeCommonSpectrum*1.6, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kTRUE);
	graphRatioToBylinkinFitPi0PP7TeV->SetLineWidth(widthCommonErrors);
	
	
	DrawGammaSetMarkerTGraphAsym(graphRatioToTsallisFitPi0PP2760GeV, markerStyleTsallisFitRatioPi02760GeV, markerSizeCommonSpectrum*1.6, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kTRUE,colorCommonSpectrumPi02760GeVBox);
	graphRatioToTsallisFitPi0PP2760GeV->SetLineWidth(widthCommonErrors);
	
	DrawGammaSetMarkerTGraphAsym(graphRatioToBylinkinFitPi0PP2760GeV, markerStyleBylinkinFitRatioPi02760GeV, markerSizeCommonSpectrum*1.6, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kTRUE);
	graphRatioToBylinkinFitPi0PP2760GeV->SetLineWidth(widthCommonErrors);
	
	
	Double_t scaleMarkerFits = 1.0;
 	Double_t scaleLineWidthFits = 1.;
	
	//*************** Offsets *************************
	
	Double_t offsetFitsLegendPrelAndFinalMarkerX = 0.09;
	Double_t offsetFitsLegendPrelAndFinalMarkerY = 0.03;
	Double_t offsetFitsLegendPrelAndFinalBox = 0.05;
	Double_t offsetFitsLegendPrelAndFinalLine = 0.06;
		
	Double_t rowsFitsLegendPrelAndFinalOnlyRatioPi0[8] = {0.88,0.77,0.58,0.43,0.36,0.34,0.21,0.09}; // {0.88,0.75,0.62,0.49,0.36,0.25,0.12}; //{0.9,0.79,0.69,0.54,0.41,0.26,0.12,0.0};
	Double_t columnsFitsLegendPrelAndFinalOnlyRatio[4] = {0.,0.25,0.51,0.76};
	
	Double_t textFitsSizeLeftLabelsPrelAndFinalRatio = 0.12; //;0.105;
	Double_t textFitsSizeTopLablesPrelAndFinalRatio = 0.12; //0.11;
 	Double_t textFitsSizeTopLowerLablesPrelAndFinalRatio =  0.12; // 0.105;
	
	
	//*************** first Column **********************************************************
	TLatex *textFitsTsallisSpectrumOnlyRatioPi0 = new TLatex(columnsFitsLegendPrelAndFinalOnlyRatio[0],rowsFitsLegendPrelAndFinalOnlyRatioPi0[2],"#pi^{0} Tsallis fit");
	SetStyleTLatex( textFitsTsallisSpectrumOnlyRatioPi0, textFitsSizeLeftLabelsPrelAndFinalRatio,4);
	
	TLatex *textFitsBylinkinSpectrumOnlyRatioPi0 = new TLatex(columnsFitsLegendPrelAndFinalOnlyRatio[0],rowsFitsLegendPrelAndFinalOnlyRatioPi0[3],"#pi^{0} Bylinkin fit");
	SetStyleTLatex( textFitsBylinkinSpectrumOnlyRatioPi0, textFitsSizeLeftLabelsPrelAndFinalRatio,4);
	
	
	//*************** second Column **********************************************************
	TLatex *textPi0pPb5023GeVFitsOnlyRatioPi0 = new TLatex(columnsFitsLegendPrelAndFinalOnlyRatio[1],rowsFitsLegendPrelAndFinalOnlyRatioPi0[0],"#pi^{0}, #sqrt{#it{s_{NN}}} = 5.02 TeV");
	SetStyleTLatex( textPi0pPb5023GeVFitsOnlyRatioPi0, textFitsSizeTopLablesPrelAndFinalRatio,4);
	TLatex *textPi0pPb5023GeVFitssysOnlyRatioPi0 = new TLatex(columnsFitsLegendPrelAndFinalOnlyRatio[1],rowsFitsLegendPrelAndFinalOnlyRatioPi0[1],"syst. + stat.");
	SetStyleTLatex( textPi0pPb5023GeVFitssysOnlyRatioPi0, textFitsSizeTopLowerLablesPrelAndFinalRatio,4);
	
	
	TBox* boxFitsTsallisPi0pPb5023GeVOnlyRatioPi0 	 = CreateBoxFromGraph(   	graphDalitzPi0InvYieldSyspPb5023GeVScale,columnsFitsLegendPrelAndFinalOnlyRatio[1]+offsetFitsLegendPrelAndFinalMarkerX-offsetFitsLegendPrelAndFinalLine, rowsFitsLegendPrelAndFinalOnlyRatioPi0[2]+ offsetFitsLegendPrelAndFinalMarkerY- offsetFitsLegendPrelAndFinalBox, columnsFitsLegendPrelAndFinalOnlyRatio[1]+offsetFitsLegendPrelAndFinalMarkerX+offsetFitsLegendPrelAndFinalLine, rowsFitsLegendPrelAndFinalOnlyRatioPi0[2]+ offsetFitsLegendPrelAndFinalMarkerY+offsetFitsLegendPrelAndFinalBox);
	boxFitsTsallisPi0pPb5023GeVOnlyRatioPi0->SetFillColor(colorCommonSpectrumPi0pPb5023GeVBox);
	boxFitsTsallisPi0pPb5023GeVOnlyRatioPi0->SetFillStyle(1001);
	TMarker* markerFitsTsallisPi0pPb5023GeVOnlyRatioPi0 = CreateMarkerFromGraph(	graphDalitzPi0InvYieldSyspPb5023GeVScale,columnsFitsLegendPrelAndFinalOnlyRatio[1]+offsetFitsLegendPrelAndFinalMarkerX,rowsFitsLegendPrelAndFinalOnlyRatioPi0[2]+offsetFitsLegendPrelAndFinalMarkerY ,scaleMarkerFits);
	markerFitsTsallisPi0pPb5023GeVOnlyRatioPi0->SetMarkerStyle(markerStyleTsallisFitRatioPi0pPb5023GeV);
	
	
	TBox* boxFitsBylinkinPi0pPb5023GeVOnlyRatioPi0 	 	= CreateBoxFromGraph(   	graphDalitzPi0InvYieldSyspPb5023GeVScale,columnsFitsLegendPrelAndFinalOnlyRatio[1]+offsetFitsLegendPrelAndFinalMarkerX-offsetFitsLegendPrelAndFinalLine, rowsFitsLegendPrelAndFinalOnlyRatioPi0[3]+ offsetFitsLegendPrelAndFinalMarkerY- offsetFitsLegendPrelAndFinalBox, columnsFitsLegendPrelAndFinalOnlyRatio[1]+offsetFitsLegendPrelAndFinalMarkerX+offsetFitsLegendPrelAndFinalLine, rowsFitsLegendPrelAndFinalOnlyRatioPi0[3]+ offsetFitsLegendPrelAndFinalMarkerY+offsetFitsLegendPrelAndFinalBox);
	TMarker* markerFitsBylinkinPi0pPb5023GeVOnlyRatioPi0 	= CreateMarkerFromGraph(	graphDalitzPi0InvYieldSyspPb5023GeVScale,columnsFitsLegendPrelAndFinalOnlyRatio[1]+offsetFitsLegendPrelAndFinalMarkerX,rowsFitsLegendPrelAndFinalOnlyRatioPi0[3]+offsetFitsLegendPrelAndFinalMarkerY ,scaleMarkerFits);
	markerFitsBylinkinPi0pPb5023GeVOnlyRatioPi0->SetMarkerStyle(markerStyleBylinkinFitRatioPi0pPb5023GeV);
	
	
	//////////////////////third column //////////////////////////////////////7
	
	
	TLatex *textPi07TeVFitsOnlyRatioPi0 = new TLatex(columnsFitsLegendPrelAndFinalOnlyRatio[2],rowsFitsLegendPrelAndFinalOnlyRatioPi0[0],"#pi^{0}, #sqrt{#it{s}} = 7 TeV");
	SetStyleTLatex( textPi07TeVFitsOnlyRatioPi0, textFitsSizeTopLablesPrelAndFinalRatio,4);
	TLatex *textPi07TeVFitssysOnlyRatioPi0 = new TLatex(columnsFitsLegendPrelAndFinalOnlyRatio[2],rowsFitsLegendPrelAndFinalOnlyRatioPi0[1],"syst. + stat.");
	SetStyleTLatex( textPi07TeVFitssysOnlyRatioPi0, textFitsSizeTopLowerLablesPrelAndFinalRatio,4);
	
	TBox* boxFitsTsallisPi07TeVOnlyRatioPi0 	 	= CreateBoxFromGraph(   	graphDalitzPi0InvYieldSys7TeV,columnsFitsLegendPrelAndFinalOnlyRatio[2]+offsetFitsLegendPrelAndFinalMarkerX-offsetFitsLegendPrelAndFinalLine, rowsFitsLegendPrelAndFinalOnlyRatioPi0[2]+ offsetFitsLegendPrelAndFinalMarkerY- offsetFitsLegendPrelAndFinalBox, columnsFitsLegendPrelAndFinalOnlyRatio[2]+offsetFitsLegendPrelAndFinalMarkerX+offsetFitsLegendPrelAndFinalLine, rowsFitsLegendPrelAndFinalOnlyRatioPi0[2]+ offsetFitsLegendPrelAndFinalMarkerY+offsetFitsLegendPrelAndFinalBox);
	boxFitsTsallisPi07TeVOnlyRatioPi0->SetFillColor(colorCommonSpectrumPi07TeVBox);
	boxFitsTsallisPi07TeVOnlyRatioPi0->SetFillStyle(1001);
	
	TMarker* markerFitsTsallisPi07TeVOnlyRatioPi0 		= CreateMarkerFromGraph(	graphDalitzPi0InvYieldSys7TeV,columnsFitsLegendPrelAndFinalOnlyRatio[2]+offsetFitsLegendPrelAndFinalMarkerX,rowsFitsLegendPrelAndFinalOnlyRatioPi0[2]+offsetFitsLegendPrelAndFinalMarkerY ,scaleMarkerFits);
	markerFitsTsallisPi07TeVOnlyRatioPi0->SetMarkerStyle(markerStyleTsallisFitRatioPi07TeV);
	
	
	TBox* boxFitsBylinkinPi07TeVOnlyRatioPi0 	 	= CreateBoxFromGraph(   	graphDalitzPi0InvYieldSys7TeV,columnsFitsLegendPrelAndFinalOnlyRatio[2]+offsetFitsLegendPrelAndFinalMarkerX-offsetFitsLegendPrelAndFinalLine, rowsFitsLegendPrelAndFinalOnlyRatioPi0[3]+ offsetFitsLegendPrelAndFinalMarkerY- offsetFitsLegendPrelAndFinalBox, columnsFitsLegendPrelAndFinalOnlyRatio[2]+offsetFitsLegendPrelAndFinalMarkerX+offsetFitsLegendPrelAndFinalLine, rowsFitsLegendPrelAndFinalOnlyRatioPi0[3]+ offsetFitsLegendPrelAndFinalMarkerY+offsetFitsLegendPrelAndFinalBox);
	TMarker* markerFitsBylinkinPi07TeVOnlyRatioPi0 		= CreateMarkerFromGraph(	graphDalitzPi0InvYieldSys7TeV,columnsFitsLegendPrelAndFinalOnlyRatio[2]+offsetFitsLegendPrelAndFinalMarkerX,rowsFitsLegendPrelAndFinalOnlyRatioPi0[3]+offsetFitsLegendPrelAndFinalMarkerY ,scaleMarkerFits);
	markerFitsBylinkinPi07TeVOnlyRatioPi0->SetMarkerStyle(markerStyleBylinkinFitRatioPi07TeV);
	
	
	/////////////fourth column/////////////////////////////
	
	
	TLatex *textPi02760GeVFitsOnlyRatioPi0 = new TLatex(columnsFitsLegendPrelAndFinalOnlyRatio[3],rowsFitsLegendPrelAndFinalOnlyRatioPi0[0],"#pi^{0}, #sqrt{#it{s}} = 2.76 TeV");
	SetStyleTLatex( textPi02760GeVFitsOnlyRatioPi0, textFitsSizeTopLablesPrelAndFinalRatio,4);
	TLatex *textPi02760GeVFitssysOnlyRatioPi0 = new TLatex(columnsFitsLegendPrelAndFinalOnlyRatio[3],rowsFitsLegendPrelAndFinalOnlyRatioPi0[1],"syst. + stat.");
	SetStyleTLatex( textPi02760GeVFitssysOnlyRatioPi0, textFitsSizeTopLowerLablesPrelAndFinalRatio,4);
	
	TBox* boxFitsTsallisPi02760GeVOnlyRatioPi0 	 	= CreateBoxFromGraph(   	graphDalitzPi0InvYieldSys2760GeVScale,columnsFitsLegendPrelAndFinalOnlyRatio[3]+offsetFitsLegendPrelAndFinalMarkerX-offsetFitsLegendPrelAndFinalLine, rowsFitsLegendPrelAndFinalOnlyRatioPi0[2]+ offsetFitsLegendPrelAndFinalMarkerY- offsetFitsLegendPrelAndFinalBox, columnsFitsLegendPrelAndFinalOnlyRatio[3]+offsetFitsLegendPrelAndFinalMarkerX+offsetFitsLegendPrelAndFinalLine, rowsFitsLegendPrelAndFinalOnlyRatioPi0[2]+ offsetFitsLegendPrelAndFinalMarkerY+offsetFitsLegendPrelAndFinalBox);
	boxFitsTsallisPi02760GeVOnlyRatioPi0->SetFillColor(colorCommonSpectrumPi02760GeVBox);
	boxFitsTsallisPi02760GeVOnlyRatioPi0->SetFillStyle(1001);
	
	TMarker* markerFitsTsallisPi02760GeVOnlyRatioPi0 	= CreateMarkerFromGraph(	graphDalitzPi0InvYieldSys2760GeVScale,columnsFitsLegendPrelAndFinalOnlyRatio[3]+offsetFitsLegendPrelAndFinalMarkerX,rowsFitsLegendPrelAndFinalOnlyRatioPi0[2]+offsetFitsLegendPrelAndFinalMarkerY ,scaleMarkerFits);
	markerFitsTsallisPi02760GeVOnlyRatioPi0->SetMarkerStyle(markerStyleTsallisFitRatioPi02760GeV);
	
	TBox* boxFitsBylinkinPi02760GeVOnlyRatioPi0 	 	= CreateBoxFromGraph(   	graphDalitzPi0InvYieldSys2760GeVScale,columnsFitsLegendPrelAndFinalOnlyRatio[3]+offsetFitsLegendPrelAndFinalMarkerX-offsetFitsLegendPrelAndFinalLine, rowsFitsLegendPrelAndFinalOnlyRatioPi0[3]+ offsetFitsLegendPrelAndFinalMarkerY- offsetFitsLegendPrelAndFinalBox, columnsFitsLegendPrelAndFinalOnlyRatio[3]+offsetFitsLegendPrelAndFinalMarkerX+offsetFitsLegendPrelAndFinalLine, rowsFitsLegendPrelAndFinalOnlyRatioPi0[3]+ offsetFitsLegendPrelAndFinalMarkerY+offsetFitsLegendPrelAndFinalBox);
	TMarker* markerFitsBylinkinPi02760GeVOnlyRatioPi0 	= CreateMarkerFromGraph(	graphDalitzPi0InvYieldSys2760GeVScale,columnsFitsLegendPrelAndFinalOnlyRatio[3]+offsetFitsLegendPrelAndFinalMarkerX,rowsFitsLegendPrelAndFinalOnlyRatioPi0[3]+offsetFitsLegendPrelAndFinalMarkerY ,scaleMarkerFits);
	markerFitsBylinkinPi02760GeVOnlyRatioPi0->SetMarkerStyle(markerStyleBylinkinFitRatioPi02760GeV);
	
	
	
	
	
	
	
	
	//TBox* boxFitTsallisPi07TeVOnlyRatioPi0 = CreateBoxFromGraph(graphDalitzPi0InvCrossSectionCompl7TeV,columnsNLOLegendPrelAndFinalOnlyRatio[1]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+ offsetNLOLegendPrelAndFinalMarkerY- offsetNLOLegendPrelAndFinalBox, columnsNLOLegendPrelAndFinalOnlyRatio[1]+offsetNLOLegendPrelAndFinalMarkerX+offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+ offsetNLOLegendPrelAndFinalMarkerY+offsetNLOLegendPrelAndFinalBox);
	//TMarker* markerCombinedPi07TeVOnlyRatioPi0 = CreateMarkerFromGraph(graphDalitzPi0InvCrossSectionCompl7TeV,columnsNLOLegendPrelAndFinalOnlyRatio[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY ,scaleMarkerNLO);
	
	//TBox* boxFitBylinkinPi07TeVOnlyRatioPi0 = CreateBoxFromGraph(graphDalitzPi0InvCrossSectionCompl7TeV,columnsNLOLegendPrelAndFinalOnlyRatio[1]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[3]+ offsetNLOLegendPrelAndFinalMarkerY- offsetNLOLegendPrelAndFinalBox, columnsNLOLegendPrelAndFinalOnlyRatio[1]+offsetNLOLegendPrelAndFinalMarkerX+offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[3]+ offsetNLOLegendPrelAndFinalMarkerY+offsetNLOLegendPrelAndFinalBox);
	
		
	//TLine * lineNLOPi07TeVMuHalfOnlyRatioPi0 = CreateLineFromGraph(graphNLOMuHalfPi07TeV, 	columnsNLOLegendPrelAndFinalOnlyRatio[1]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[3]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinalOnlyRatio[1]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[3]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO*1.5); 
	//TLine * lineNLOPi07TeVMuOneOnlyRatioPi0  = CreateLineFromGraph(graphNLOMuOnePi07TeV,	columnsNLOLegendPrelAndFinalOnlyRatio[1]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[4]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinalOnlyRatio[1]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[4]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO*1.5); 
	//TLine * lineNLOPi07TeVMuTwoOnlyRatioPi0  = CreateLineFromGraph(graphNLOMuTwoPi07TeV, 	columnsNLOLegendPrelAndFinalOnlyRatio[1]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[5]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinalOnlyRatio[1]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[5]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO*1.5); 
	
	
	
	
		
	
	
	
	//////////////////////////////////////////////////////////////////////////////////////////////
	
	
	
	
	  
        TCanvas* canvasDummy001 = new TCanvas("canvasDummy001","",1150,1000);  // gives the page size
	DrawGammaCanvasSettings( canvasDummy001,  0.18, 0.02, 0.03, 0.1);
	canvasDummy001->SetLogy(1);
	canvasDummy001->SetLogx(1);
	
	
	TLegend* legendCorrectedYieldMeson001 = new TLegend(0.2,0.18,0.38,0.40);
	legendCorrectedYieldMeson001->SetTextSize(0.040);
	legendCorrectedYieldMeson001->SetFillColor(0);
	legendCorrectedYieldMeson001->SetLineColor(0);
	legendCorrectedYieldMeson001->SetTextFont(42);
	
	TLegend* legendPlotLabel = new TLegend(0.71,0.85,0.95,0.95);
	legendPlotLabel->SetTextSize(0.040);
	legendPlotLabel->SetFillColor(0);
	legendPlotLabel->SetLineColor(0);
	legendPlotLabel->SetTextFont(42);
	legendPlotLabel->AddEntry( (TObject*)0,thesisPlotLabel.Data(),"");
	legendPlotLabel->AddEntry( (TObject*)0,DalitzLabel.Data(),"");
	
	
	
	
	
	TH2F * histo2DDummy001;
	histo2DDummy001 = new TH2F("histo2DDummy001","histo2DDummy001",1000,0.2,15.,10000,1e-9,1e2);
	//histo2DDummy001->GetXaxis()->SetRangeUser(0.5,15);
	SetStyleHistoTH2ForGraphs(histo2DDummy001,pTLabel.Data(), invYieldLabel.Data(), 0.04,0.04, 0.04,0.04, 1.2,1.8);
	
	histo2DDummy001->GetXaxis()->SetRangeUser(0.5,12);
	histo2DDummy001->GetYaxis()->SetRangeUser(8e-9,50);
	
	//SetStyleHistoTH2ForGraphs(histo2DDummy001,invYieldLabel.Data(), 0.04,0.04, 0.04,0.04, 1,1.55);
	histo2DDummy001->DrawCopy(); 
	       
	
	       
	    
	graphDalitzPi0InvYieldSys2760GeVScale->Draw("2same");    
	graphDalitzPi0InvYieldStat2760GeVScale->Draw("p,same,e");
	fitHistogramTsallisPi02760GeV->DrawCopy("lsame");
	      
	graphDalitzPi0InvYieldSys7TeV->Draw("2same");       
	graphDalitzPi0InvYieldStat7TeV->Draw("p,same,e");
	fitHistogramTsallisPi07TeV->DrawCopy("lsame");
	      
	      
	graphDalitzPi0InvYieldSyspPb5023GeVScale->Draw("2same");
	graphDalitzPi0InvYieldStatpPb5023GeVScale->Draw("p,same,e");
	fitHistogramTsallisPi0pPb5023GeV->DrawCopy("lsame");
	    
	      
	      
	     
	legendCorrectedYieldMeson001->AddEntry(graphDalitzPi0InvYieldSyspPb5023GeVScale,Form("%s x 10^{1}",collisionSystem5023GeV.Data()));
	legendCorrectedYieldMeson001->AddEntry(graphDalitzPi0InvYieldSys7TeV,collisionSystem7TeV.Data());
	legendCorrectedYieldMeson001->AddEntry(graphDalitzPi0InvYieldSys2760GeVScale,Form("%s x 10^{-1}",collisionSystem2760GeV.Data()));
	legendCorrectedYieldMeson001->AddEntry(fitTsallisLineStyle,"Tsallis Fit");
	      
	legendCorrectedYieldMeson001->Draw();
	legendPlotLabel->Draw();
    
	       
	
	canvasDummy001->Update();
	canvasDummy001->Print(Form("%s/InvariantYieldPPandpPb_Tsallis.%s",outputDir.Data(),suffix.Data()));
	
	
	TH2F * histo2DRatioDummy001;
	histo2DRatioDummy001 = new TH2F("histo2DRatioDummy001","histo2DRatioDummy001",1000,0.23,70.,1000,0.001,3.6);
	histo2DRatioDummy001->GetXaxis()->SetRangeUser(0.5,15);
	histo2DRatioDummy001->GetYaxis()->SetRangeUser(0.001,3.0);
	SetStyleHistoTH2ForGraphs(histo2DRatioDummy001,"#it{p}_{T} (GeV/#it{c})","#frac{spectrum}{fit}", 0.032,0.04, 0.04,0.04, 1,1.55);
	histo2DRatioDummy001->DrawCopy(); 
	  
	
	  
        TCanvas* canvasDummy002 = new TCanvas("canvasDummy002","",1150,1000);  // gives the page size
		
	DrawGammaCanvasSettings( canvasDummy002,  0.18, 0.02, 0.03, 0.1);
	canvasDummy002->SetLogy(1);
	canvasDummy002->SetLogx(1);
	
	
	TLegend* legendCorrectedYieldMeson002 = new TLegend(0.2,0.18,0.38,0.40);
	legendCorrectedYieldMeson002->SetTextSize(0.04);
	legendCorrectedYieldMeson002->SetFillColor(0);
	legendCorrectedYieldMeson002->SetLineColor(0);
	legendCorrectedYieldMeson002->SetTextFont(42);
	
	
	
	TH2F * histo2DDummy002;
	       
	histo2DDummy002 = new TH2F("histo2DDummy002","histo2DDummy002",1000,0.2,15.,10000,1e-9,1e2);
	SetStyleHistoTH2ForGraphs(histo2DDummy002,pTLabel.Data(), invYieldLabel.Data(), 0.04,0.04, 0.04,0.04, 1.2,1.8);
	
	histo2DDummy002->GetXaxis()->SetRangeUser(0.5,12);
	histo2DDummy002->GetYaxis()->SetRangeUser(8e-9,50);
	
	histo2DDummy002->DrawCopy(); 
	
	graphDalitzPi0InvYieldSys2760GeVScale->Draw("2same");
	graphDalitzPi0InvYieldStat2760GeVScale->Draw("p,same,e");
	fitHistogramBylinkinPi02760GeV->DrawCopy("lsame");
	
	
	graphDalitzPi0InvYieldSys7TeV->Draw("2same");   
	graphDalitzPi0InvYieldStat7TeV->Draw("p,same,e");
	fitHistogramBylinkinPi07TeV->DrawCopy("lsame");
	
	graphDalitzPi0InvYieldSyspPb5023GeVScale->Draw("2same");
	graphDalitzPi0InvYieldStatpPb5023GeVScale->Draw("p,same,e");
	      
	fitHistogramBylinkinPi0pPb5023GeV->DrawCopy("lsame");
	      
	      
	      
	      
	legendCorrectedYieldMeson002->AddEntry(graphDalitzPi0InvYieldSyspPb5023GeVScale,collisionSystem5023GeV.Data());
	legendCorrectedYieldMeson002->AddEntry(graphDalitzPi0InvYieldSys7TeV,collisionSystem7TeV.Data());
	legendCorrectedYieldMeson002->AddEntry(graphDalitzPi0InvYieldSys2760GeVScale,Form("%s x 10^{-1}",collisionSystem2760GeV.Data()));
	legendCorrectedYieldMeson002->AddEntry(fitBylinkinLineStyle,"Bylinkin Fit");
	      
	legendCorrectedYieldMeson002->Draw();
	legendPlotLabel->Draw();
    
	       
	
	canvasDummy002->Update();
	canvasDummy002->Print(Form("%s/InvariantYieldPPandpPb_Bylinkin.%s",outputDir.Data(),suffix.Data()));
	
	
	
	TCanvas* canvasInvYieldFitsOnlySpectraPi0 = new TCanvas("canvasInvYieldFitsOnlySpectraPi0","",200,10,1200,1200);  // gives the page size
	DrawGammaCanvasSettings( canvasInvYieldFitsOnlySpectraPi0,  0.15, 0.02, 0.03, 0.12);
	
	canvasInvYieldFitsOnlySpectraPi0->SetLogy(1);
	canvasInvYieldFitsOnlySpectraPi0->SetLogx(1);
	
	TH2F * histo2DInvYieldFitsOnlySpectraPi0 = new TH2F("histo2DInvYieldFitsOnlySpectraPi0","histo2DInvYieldFitsOnlySpectraPi0",1000,0.2,15.,10000,1e-9,1e2);
	SetStyleHistoTH2ForGraphs(histo2DInvYieldFitsOnlySpectraPi0,pTLabel.Data(), invYieldLabel.Data(), 36,36,36,36,1.5,2.0);
	
	histo2DInvYieldFitsOnlySpectraPi0->GetXaxis()->SetLabelFont(43);
	histo2DInvYieldFitsOnlySpectraPi0->GetYaxis()->SetLabelFont(43); 
	histo2DInvYieldFitsOnlySpectraPi0->GetXaxis()->SetTitleFont(63);
	histo2DInvYieldFitsOnlySpectraPi0->GetYaxis()->SetTitleFont(63);
	
	histo2DInvYieldFitsOnlySpectraPi0->GetXaxis()->SetRangeUser(0.5,12);
	histo2DInvYieldFitsOnlySpectraPi0->GetYaxis()->SetRangeUser(8e-9,50);
	histo2DInvYieldFitsOnlySpectraPi0->DrawCopy(); 
	       
	
	
        graphDalitzPi0InvYieldSys2760GeVScale->Draw("2,same");
	graphDalitzPi0InvYieldStat2760GeVScale->Draw("p,same");
	
	
	      
             
        fitHistogramTsallisPi02760GeV->DrawCopy("lsame");
	fitHistogramBylinkinPi02760GeV->DrawCopy("lsame");
	      
	 	     
	graphDalitzPi0InvYieldSys7TeV->Draw("2,same");      
	graphDalitzPi0InvYieldStat7TeV->Draw("p,same");
	
	
	fitHistogramTsallisPi07TeV->DrawCopy("lsame");
	      
	fitHistogramBylinkinPi07TeV->DrawCopy("lsame");
	      
	graphDalitzPi0InvYieldSyspPb5023GeVScale->Draw("2,same");
	graphDalitzPi0InvYieldStatpPb5023GeVScale->Draw("p,same");
	   
	      
	fitHistogramTsallisPi0pPb5023GeV->DrawCopy("lsame");
	fitHistogramBylinkinPi0pPb5023GeV->DrawCopy("lsame");
	      
	      
	TLegend* legendCorrectedYieldMeson003 = new TLegend(0.21,0.16,0.39,0.38);
	legendCorrectedYieldMeson003->SetTextSize(36);
	legendCorrectedYieldMeson003->SetFillColor(0);
	legendCorrectedYieldMeson003->SetLineColor(0);
	legendCorrectedYieldMeson003->SetTextFont(43);
	     
	legendCorrectedYieldMeson003->AddEntry(graphDalitzPi0InvYieldSyspPb5023GeVScale,Form("%s x 10^{1}",collisionSystem5023GeV.Data()));
	legendCorrectedYieldMeson003->AddEntry(graphDalitzPi0InvYieldSys7TeV,collisionSystem7TeV.Data());
	legendCorrectedYieldMeson003->AddEntry(graphDalitzPi0InvYieldSys2760GeVScale,Form("%s x 10^{-1}",collisionSystem2760GeV.Data()));
	legendCorrectedYieldMeson003->AddEntry(fitTsallisLineStyle,"Tsallis Fit");
	legendCorrectedYieldMeson003->AddEntry(fitBylinkinLineStyle,"Bylinkin Fit");
	legendCorrectedYieldMeson003->Draw();
	
	legendPlotLabel->Draw();
    
	       
	
	canvasInvYieldFitsOnlySpectraPi0->Update();
	canvasInvYieldFitsOnlySpectraPi0->Print(Form("%s/InvariantYieldPPandpPb_Tsallis_Bylinkin.%s",outputDir.Data(),suffix.Data()));
	
	
	
	
	TCanvas* canvasInvYieldFitsOnlyRatioPi0 	= new TCanvas("canvasInvYieldFitsOnlyRatioPi0","",200,10,1200,1200);  // gives the page size
	DrawGammaCanvasSettings( canvasInvYieldFitsOnlyRatioPi0,  0.15, 0.02, 0.03, 0.06);  
	
	TPad* padInvYieldFitsLegendOnlyRatioPi0 		= new TPad("padInvYieldFitsLegendOnlyRatioPi0", 	   "", 0.15, 0.8, 0.98,  1.,-1, -1, -2);
 	DrawGammaPadSettings( padInvYieldFitsLegendOnlyRatioPi0, 0.15, 0.02, 0.01, 0.0);
 	padInvYieldFitsLegendOnlyRatioPi0->Draw();
	
	TPad*  padInvYieldFitsLegendOnlyRatioPi0pPb5023GeV 	= new TPad("padInvYieldFitsLegendOnlyRatioPi0pPb5023GeV",  "", 0.00, 0.565, 1., 0.8,-1, -1, -2);
 	DrawGammaPadSettings( padInvYieldFitsLegendOnlyRatioPi0pPb5023GeV, 0.15, 0.02, 0.02, 0.);
 	padInvYieldFitsLegendOnlyRatioPi0pPb5023GeV->Draw();

 	TPad* padInvYieldFitsLegendOnlyRatioPi07TeV 	= new TPad("padInvYieldFitsLegendOnlyRatioPi07TeV",    "", 0., 0.33, 1., 0.565,-1, -1, -2);
 	DrawGammaPadSettings( padInvYieldFitsLegendOnlyRatioPi07TeV, 0.15, 0.02, 0.00, 0.);
 	padInvYieldFitsLegendOnlyRatioPi07TeV->Draw();
 	
  	TPad* padInvYieldFitsLegendOnlyRatioPi02760GeV	= new TPad("padInvYieldFitsLegendOnlyRatioPi02760GeV", "", 0.,0.0, 1., 0.33,-1, -1, -2);
  	DrawGammaPadSettings( padInvYieldFitsLegendOnlyRatioPi02760GeV, 0.15, 0.02, 0., 0.28);
  	padInvYieldFitsLegendOnlyRatioPi02760GeV->Draw();
	
	
	
	


  	padInvYieldFitsLegendOnlyRatioPi0pPb5023GeV->cd();
  	padInvYieldFitsLegendOnlyRatioPi0pPb5023GeV->SetLogx();
  	
  	TH2F * ratio2DInvYieldsFitsOnlyRatioPi0pPb5023GeV = new TH2F("ratio2DInvYieldsFitsOnlyRatioPi0pPb5023GeV","ratio2DInvYieldsFitsOnlyRatioPi0pPb5023GeV",1000,0.23,30.,1000,0.4,3.55);
  	SetStyleHistoTH2ForGraphs(ratio2DInvYieldsFitsOnlyRatioPi0pPb5023GeV, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{fit}", 36,36,36,36,3.5,1.2, 512, 505);

	//Jus for differen size of Tpads
  	ratio2DInvYieldsFitsOnlyRatioPi0pPb5023GeV->GetXaxis()->SetLabelFont(43);
  	ratio2DInvYieldsFitsOnlyRatioPi0pPb5023GeV->GetYaxis()->SetLabelFont(43); 
  	ratio2DInvYieldsFitsOnlyRatioPi0pPb5023GeV->GetXaxis()->SetTitleFont(63);
  	ratio2DInvYieldsFitsOnlyRatioPi0pPb5023GeV->GetYaxis()->SetTitleFont(63);
	ratio2DInvYieldsFitsOnlyRatioPi0pPb5023GeV->GetYaxis()->CenterTitle();
  	
  	ratio2DInvYieldsFitsOnlyRatioPi0pPb5023GeV->DrawCopy(); 
  	
 	graphRatioToTsallisFitPi0pPb5023TeV->Draw("p,E2same");
	graphRatioToBylinkinFitPi0pPb5023TeV->Draw("p,E2same");
	
 	
  	labelRatioFitsPi0pPb5023GeV->Draw();
  	
  	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	padInvYieldFitsLegendOnlyRatioPi07TeV->cd();
  	padInvYieldFitsLegendOnlyRatioPi07TeV->SetLogx();
  	
  	TH2F * ratio2DInvYieldsFitsOnlyRatioPi07TeV = new TH2F("ratio2DInvYieldsFitsOnlyRatioPi07TeV","ratio2DInvYieldsFitsOnlyRatioPi07TeV",1000,0.23,30.,1000,0.4,3.55);
  	SetStyleHistoTH2ForGraphs(ratio2DInvYieldsFitsOnlyRatioPi07TeV, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{fit}", 36,36,36,36,3.5,1.2, 512, 505);

	//Jus for differen size of Tpads
  	ratio2DInvYieldsFitsOnlyRatioPi07TeV->GetXaxis()->SetLabelFont(43);
  	ratio2DInvYieldsFitsOnlyRatioPi07TeV->GetYaxis()->SetLabelFont(43); 
  	ratio2DInvYieldsFitsOnlyRatioPi07TeV->GetXaxis()->SetTitleFont(63);
  	ratio2DInvYieldsFitsOnlyRatioPi07TeV->GetYaxis()->SetTitleFont(63);
	ratio2DInvYieldsFitsOnlyRatioPi07TeV->GetYaxis()->CenterTitle();
  	ratio2DInvYieldsFitsOnlyRatioPi07TeV->DrawCopy(); 
 
	graphRatioToTsallisFitPi0PP7TeV->Draw("p,E2same");
	graphRatioToBylinkinFitPi0PP7TeV->Draw("p,E2same");
	
 	
  	labelRatioFitsPi07TeV->Draw();
  	
  	DrawGammaLines(0., 30.,1., 1.,0.1);
	
		
  	padInvYieldFitsLegendOnlyRatioPi02760GeV->cd();
  	padInvYieldFitsLegendOnlyRatioPi02760GeV->SetLogx();
	
	
	TH2F * ratio2DInvYieldsFitsOnlyRatioPi02760GeV = new TH2F("ratio2DInvYieldsFitsOnlyRatioPi02760GeV","ratio2DInvYieldsFitsOnlyRatioPi02760GeV",1000,0.23,30.,1000,0.4,3.55);
  	SetStyleHistoTH2ForGraphs(ratio2DInvYieldsFitsOnlyRatioPi02760GeV, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{fit}", 36,36,36,36,3.5,1.2, 512, 505);

	//Jus for differen size of Tpads
  	ratio2DInvYieldsFitsOnlyRatioPi02760GeV->GetXaxis()->SetLabelFont(43);
  	ratio2DInvYieldsFitsOnlyRatioPi02760GeV->GetYaxis()->SetLabelFont(43); 
  	ratio2DInvYieldsFitsOnlyRatioPi02760GeV->GetXaxis()->SetTitleFont(63);
  	ratio2DInvYieldsFitsOnlyRatioPi02760GeV->GetYaxis()->SetTitleFont(63);
	ratio2DInvYieldsFitsOnlyRatioPi02760GeV->GetYaxis()->CenterTitle();
  	
  	ratio2DInvYieldsFitsOnlyRatioPi02760GeV->DrawCopy(); 
  	
 	graphRatioToTsallisFitPi0PP2760GeV->Draw("p,E2same");
	graphRatioToBylinkinFitPi0PP2760GeV->Draw("p,E2same");
	
 	
  	labelRatioFitsPi02760GeV->Draw();
  	
  	DrawGammaLines(0., 30.,1., 1.,0.1);
	

 			
 	padInvYieldFitsLegendOnlyRatioPi0->cd();
 	DrawGammaPadSettings( padInvYieldFitsLegendOnlyRatioPi0, 0., 0., 0.0, 0.0);
 	padInvYieldFitsLegendOnlyRatioPi0->SetBorderMode(-1);
 	padInvYieldFitsLegendOnlyRatioPi0->SetBorderSize(3);
 	padInvYieldFitsLegendOnlyRatioPi0->Draw();
 	padInvYieldFitsLegendOnlyRatioPi0->cd();
 	  	//*************** first Column **********************************************************
  	textFitsTsallisSpectrumOnlyRatioPi0->Draw();
	textFitsBylinkinSpectrumOnlyRatioPi0->Draw();

	 	//*************** second Column **********************************************************
  	textPi0pPb5023GeVFitsOnlyRatioPi0->Draw();
	textPi0pPb5023GeVFitssysOnlyRatioPi0->Draw();
	
	
	boxFitsTsallisPi0pPb5023GeVOnlyRatioPi0->Draw();
	markerFitsTsallisPi0pPb5023GeVOnlyRatioPi0->DrawMarker(columnsFitsLegendPrelAndFinalOnlyRatio[1]+offsetFitsLegendPrelAndFinalMarkerX,rowsFitsLegendPrelAndFinalOnlyRatioPi0[2]+offsetFitsLegendPrelAndFinalMarkerY);
	
	
	boxFitsBylinkinPi0pPb5023GeVOnlyRatioPi0->Draw();
	markerFitsBylinkinPi0pPb5023GeVOnlyRatioPi0->DrawMarker(columnsFitsLegendPrelAndFinalOnlyRatio[1]+offsetFitsLegendPrelAndFinalMarkerX,rowsFitsLegendPrelAndFinalOnlyRatioPi0[3]+offsetFitsLegendPrelAndFinalMarkerY);
 	
 	
//  	
//  	lineNLOPi07TeVMuHalfOnlyRatioPi0->Draw("same");
//  	lineNLOPi07TeVMuOneOnlyRatioPi0->Draw("same");
//  	lineNLOPi07TeVMuTwoOnlyRatioPi0->Draw("same");
//  	
// 	
//  	//*************** third Column **********************************************************
	
	
	textPi07TeVFitsOnlyRatioPi0->Draw();
	textPi07TeVFitssysOnlyRatioPi0->Draw();
	
	
	boxFitsTsallisPi07TeVOnlyRatioPi0->Draw();
	markerFitsTsallisPi07TeVOnlyRatioPi0->DrawMarker(columnsFitsLegendPrelAndFinalOnlyRatio[2]+offsetFitsLegendPrelAndFinalMarkerX,rowsFitsLegendPrelAndFinalOnlyRatioPi0[2]+offsetFitsLegendPrelAndFinalMarkerY);
	
	
	boxFitsBylinkinPi07TeVOnlyRatioPi0->Draw();
	markerFitsBylinkinPi07TeVOnlyRatioPi0->DrawMarker(columnsFitsLegendPrelAndFinalOnlyRatio[2]+offsetFitsLegendPrelAndFinalMarkerX,rowsFitsLegendPrelAndFinalOnlyRatioPi0[3]+offsetFitsLegendPrelAndFinalMarkerY);
	
	//****************Fourth column**********************************
	
	//textPi02670GeVFitsOnlyRatioPi0
	textPi02760GeVFitsOnlyRatioPi0->Draw();
	textPi02760GeVFitssysOnlyRatioPi0->Draw();
	
	
	boxFitsTsallisPi02760GeVOnlyRatioPi0->Draw();
	markerFitsTsallisPi02760GeVOnlyRatioPi0->DrawMarker(columnsFitsLegendPrelAndFinalOnlyRatio[3]+offsetFitsLegendPrelAndFinalMarkerX,rowsFitsLegendPrelAndFinalOnlyRatioPi0[2]+offsetFitsLegendPrelAndFinalMarkerY);
	
	
	boxFitsBylinkinPi02760GeVOnlyRatioPi0->Draw();
	markerFitsBylinkinPi02760GeVOnlyRatioPi0->DrawMarker(columnsFitsLegendPrelAndFinalOnlyRatio[3]+offsetFitsLegendPrelAndFinalMarkerX,rowsFitsLegendPrelAndFinalOnlyRatioPi0[3]+offsetFitsLegendPrelAndFinalMarkerY);
 	
  		
	canvasInvYieldFitsOnlyRatioPi0->Update();
	canvasInvYieldFitsOnlyRatioPi0->Print(Form("%s/InvariantYieldPPandpPb_Tsallis_Bylinkin_Ratio.%s",outputDir.Data(),suffix.Data()));
	
	
	
	TCanvas* canvasRatiotoFitBylinkin2760GeV = new TCanvas("canvasRatiotoFitBylinkin2760GeV","",200,10,1500,900);  // gives the page size
	canvasRatiotoFitBylinkin2760GeV->SetTickx();
	canvasRatiotoFitBylinkin2760GeV->SetTicky();
	canvasRatiotoFitBylinkin2760GeV->SetGridx(0);
	canvasRatiotoFitBylinkin2760GeV->SetGridy(0);
	canvasRatiotoFitBylinkin2760GeV->SetLogx(1);
	canvasRatiotoFitBylinkin2760GeV->SetLeftMargin(0.15);
	canvasRatiotoFitBylinkin2760GeV->SetRightMargin(0.02);
	canvasRatiotoFitBylinkin2760GeV->SetTopMargin(0.02);
	canvasRatiotoFitBylinkin2760GeV->SetBottomMargin(0.1);
	canvasRatiotoFitBylinkin2760GeV->SetFillColor(0);	
	
	
	TLegend* legendRatioToFitBylinkin2760GeV = new TLegend(0.20,0.18,0.45,0.41);
	legendRatioToFitBylinkin2760GeV->SetTextSize(0.05);
	legendRatioToFitBylinkin2760GeV->SetFillColor(0);
	legendRatioToFitBylinkin2760GeV->SetLineColor(0);
	
	histo2DRatioDummy001->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphRatioToBylinkinFitPi0SysPP2760GeV, 20, 1.5, color[2], color[2], 1, kTRUE);
	graphRatioToBylinkinFitPi0SysPP2760GeV->Draw("2same");
	
	DrawGammaSetMarkerTGraphAsym(graphRatioToBylinkinFitPi0StatPP2760GeV, 20, 1.5, color[2], color[2]);
	graphRatioToBylinkinFitPi0StatPP2760GeV->Draw("p,same,e");
	
	
	DrawGammaLines(0.23, 70.,1., 1.,0.1,kGray);
	    
       
	
	canvasRatiotoFitBylinkin2760GeV->Update();
	canvasRatiotoFitBylinkin2760GeV->Print(Form("%s/RatioToFit_PP2760GeV_Bylinkin.%s",outputDir.Data(),suffix.Data()));
	
	
	TCanvas* canvasRatiotoFitBylinkin7TeV = new TCanvas("canvasRatiotoFitBylinkin7TeV","",200,10,1500,900);  // gives the page size
	canvasRatiotoFitBylinkin7TeV->SetTickx();
	canvasRatiotoFitBylinkin7TeV->SetTicky();
	canvasRatiotoFitBylinkin7TeV->SetGridx(0);
	canvasRatiotoFitBylinkin7TeV->SetGridy(0);
	canvasRatiotoFitBylinkin7TeV->SetLogx(1);
	canvasRatiotoFitBylinkin7TeV->SetLeftMargin(0.15);
	canvasRatiotoFitBylinkin7TeV->SetRightMargin(0.02);
	canvasRatiotoFitBylinkin7TeV->SetTopMargin(0.02);
	canvasRatiotoFitBylinkin7TeV->SetBottomMargin(0.1);
	canvasRatiotoFitBylinkin7TeV->SetFillColor(0);	
	
	
	TLegend* legendRatioToFitBylinkin7TeV = new TLegend(0.20,0.18,0.45,0.41);
	legendRatioToFitBylinkin7TeV->SetTextSize(0.05);
	legendRatioToFitBylinkin7TeV->SetFillColor(0);
	legendRatioToFitBylinkin7TeV->SetLineColor(0);
	
	histo2DRatioDummy001->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphRatioToBylinkinFitPi0SysPP7TeV, 20, 1.5, color[2], color[2], 1, kTRUE);
	graphRatioToBylinkinFitPi0SysPP7TeV->Draw("2same");
	
	DrawGammaSetMarkerTGraphAsym(graphRatioToBylinkinFitPi0StatPP7TeV, 20, 1.5, color[2], color[2]);
	graphRatioToBylinkinFitPi0StatPP7TeV->Draw("p,same,e");
	
	
	DrawGammaLines(0.23, 70.,1., 1.,0.1,kGray);
	    
       
	
	canvasRatiotoFitBylinkin7TeV->Update();
	canvasRatiotoFitBylinkin7TeV->Print(Form("%s/RatioToFit_PP7TeV_Bylinkin.%s",outputDir.Data(),suffix.Data()));
	
	
		
	
	
	TCanvas* canvasRatiotoFitBylinkinpPb5023GeV = new TCanvas("canvasRatiotoFitBylinkinpPb5023GeV","",200,10,1500,900);  // gives the page size
	canvasRatiotoFitBylinkinpPb5023GeV->SetTickx();
	canvasRatiotoFitBylinkinpPb5023GeV->SetTicky();
	canvasRatiotoFitBylinkinpPb5023GeV->SetGridx(0);
	canvasRatiotoFitBylinkinpPb5023GeV->SetGridy(0);
	canvasRatiotoFitBylinkinpPb5023GeV->SetLogx(1);
	canvasRatiotoFitBylinkinpPb5023GeV->SetLeftMargin(0.15);
	canvasRatiotoFitBylinkinpPb5023GeV->SetRightMargin(0.02);
	canvasRatiotoFitBylinkinpPb5023GeV->SetTopMargin(0.02);
	canvasRatiotoFitBylinkinpPb5023GeV->SetBottomMargin(0.1);
	canvasRatiotoFitBylinkinpPb5023GeV->SetFillColor(0);	
	
	
	TLegend* legendRatioToFitBylinkinpPb5023GeV = new TLegend(0.21,0.18,0.45,0.41);
	legendRatioToFitBylinkinpPb5023GeV->SetTextSize(0.05);
	legendRatioToFitBylinkinpPb5023GeV->SetFillColor(0);
	legendRatioToFitBylinkinpPb5023GeV->SetLineColor(0);
	
	histo2DRatioDummy001->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphRatioToBylinkinFitPi0SyspPb5023TeV, 20, 1.5, color[2], color[2], 1, kTRUE);
	graphRatioToBylinkinFitPi0SyspPb5023TeV->Draw("2same");
	
	DrawGammaSetMarkerTGraphAsym(graphRatioToBylinkinFitPi0StatpPb5023TeV, 20, 1.5, color[2], color[2]);
	graphRatioToBylinkinFitPi0StatpPb5023TeV->Draw("p,same,e");
	
	
	DrawGammaLines(0.23, 70.,1., 1.,0.1,kGray);
	    
       
	
	canvasRatiotoFitBylinkinpPb5023GeV->Update();
	canvasRatiotoFitBylinkinpPb5023GeV->Print(Form("%s/RatioToFit_pPb5023GeV_Bylinkin.%s",outputDir.Data(),suffix.Data()));
	

	
	
	//Set histograms///
	TGraphAsymmErrors* graphDalitzPi0InvCrossSectionCompl2760GeVScale = (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvCrossSectionCompl2760GeV,1e-1);
	
	TGraphAsymmErrors* graphDalitzPi0InvCrossSectionSys2760GeVScale   = (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvCrossSectionSys2760GeV,1e-1);
	TGraphAsymmErrors* graphDalitzPi0InvCrossSectionStat2760GeVScale   = (TGraphAsymmErrors*)ScaleGraph(graphDalitzPi0InvCrossSectionStat2760GeV,1e-1);
	
	
	DrawGammaSetMarkerTGraphAsym(graphDalitzPi0InvCrossSectionCompl2760GeVScale, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kTRUE);
	graphDalitzPi0InvCrossSectionCompl2760GeVScale->SetLineWidth(widthCommonErrors);
	
	
	DrawGammaSetMarkerTGraphAsym(graphDalitzPi0InvCrossSectionStat7TeV, markerStyleCommmonSpectrumPi07TeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kFALSE);
	graphDalitzPi0InvCrossSectionStat7TeV->SetLineWidth(widthCommonErrors);	
	DrawGammaSetMarkerTGraphAsym(graphDalitzPi0InvCrossSectionSys7TeV, markerStyleCommmonSpectrumPi07TeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kTRUE, colorCommonSpectrumPi07TeVBox);
	graphDalitzPi0InvCrossSectionSys7TeV->SetLineWidth(widthCommonErrors);	
	
	
	DrawGammaSetMarkerTGraphAsym(graphDalitzPi0InvCrossSectionStat2760GeVScale, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kFALSE);
	graphDalitzPi0InvCrossSectionStat2760GeVScale->SetLineWidth(widthCommonErrors);
	
	DrawGammaSetMarkerTGraphAsym(graphDalitzPi0InvCrossSectionSys2760GeVScale, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kTRUE, colorCommonSpectrumPi02760GeVBox);
	graphDalitzPi0InvCrossSectionSys2760GeVScale->SetLineWidth(widthCommonErrors);
	
	
	
	
	//-------------------------Pi0 2760 GeV NLO mu = 0.5pt---------------------------------------------
	
	
	
						    
	TGraph* graphNLOMuHalfPi02760GeV = (TGraph*)graphNLOCalcMuHalf2760GeV->Clone();
	DrawGammaSetMarkerTGraph(graphNLOMuHalfPi02760GeV, styleMarkerNLOMuHalf, sizeMarkerNLO, colorNLOPi02760GeVMuHalf, colorNLOPi02760GeVMuHalf );
	
	graphNLOMuHalfPi02760GeV= ScaleGraph(graphNLOMuHalfPi02760GeV,1e-1);
	DrawGammaNLOTGraph( graphNLOMuHalfPi02760GeV, widthCommonFit, styleLineNLOMuHalf, colorNLOPi02760GeVMuHalf);
	
	TGraphAsymmErrors* graphNLODSS14MuHalfPi02760GeV = (TGraphAsymmErrors*)graphNLOCalcmuHalfDSS14InvSecPi02760GeV->Clone();
	graphNLODSS14MuHalfPi02760GeV= ScaleGraph(graphNLODSS14MuHalfPi02760GeV,1e-1);
	graphNLODSS14MuHalfPi02760GeV->SetLineWidth(widthCommonFit);
	graphNLODSS14MuHalfPi02760GeV->SetLineColor(colorNLODSS14Pi02760GeVMuHalf);
	graphNLODSS14MuHalfPi02760GeV->SetLineStyle(styleLineNLODSS14MuHalf);
	
	
	
	
	
	//------------------------- Pi0 2760 GeV NLO mu = pt ----------------------------
	TGraph* graphNLOMuOnePi02760GeV = (TGraph*)graphNLOCalcMuOne2760GeV->Clone();
	DrawGammaSetMarkerTGraph(graphNLOMuOnePi02760GeV, styleMarkerNLOMuOne, sizeMarkerNLO, colorNLOPi02760GeVMuOne, colorNLOPi02760GeVMuOne);
	graphNLOMuOnePi02760GeV= ScaleGraph(graphNLOMuOnePi02760GeV,1e-1);
	DrawGammaNLOTGraph( graphNLOMuOnePi02760GeV, widthCommonFit, styleLineNLOMuOne, colorNLOPi02760GeVMuOne);
	
	TGraphAsymmErrors* graphNLODSS14MuOnePi02760GeV = (TGraphAsymmErrors*)graphNLOCalcmuOneDSS14InvSecPi02760GeV->Clone();
	graphNLODSS14MuOnePi02760GeV= ScaleGraph(graphNLODSS14MuOnePi02760GeV,1e-1);
	graphNLODSS14MuOnePi02760GeV->SetLineWidth(widthCommonFit);
	graphNLODSS14MuOnePi02760GeV->SetLineColor(colorNLODSS14Pi02760GeVMuOne);
	graphNLODSS14MuOnePi02760GeV->SetLineStyle(styleLineNLODSS14MuOne);
	
	
	
	//------------------------- Pi0 2760 GeV NLO mu = 2pt -----------------------------
	TGraph* graphNLOMuTwoPi02760GeV = (TGraph*)graphNLOCalcMuTwo2760GeV->Clone();
	DrawGammaSetMarkerTGraph(graphNLOMuTwoPi02760GeV, styleMarkerNLOMuTwo, sizeMarkerNLO, colorNLOPi02760GeVMuTwo, colorNLOPi02760GeVMuTwo);
	graphNLOMuTwoPi02760GeV = ScaleGraph(graphNLOMuTwoPi02760GeV,1e-1);
	DrawGammaNLOTGraph( graphNLOMuTwoPi02760GeV, widthCommonFit, styleLineNLOMuTwo, colorNLOPi02760GeVMuTwo);
	
	
	
	
	TGraphAsymmErrors* graphNLODSS14MuTwoPi02760GeV = (TGraphAsymmErrors*)graphNLOCalcmuTwoDSS14InvSecPi02760GeV->Clone();
	graphNLODSS14MuTwoPi02760GeV= ScaleGraph(graphNLODSS14MuTwoPi02760GeV,1e-1);
	graphNLODSS14MuTwoPi02760GeV->SetLineWidth(widthCommonFit);
	graphNLODSS14MuTwoPi02760GeV->SetLineColor(colorNLODSS14Pi02760GeVMuTwo);
	graphNLODSS14MuTwoPi02760GeV->SetLineStyle(styleLineNLODSS14MuTwo);
	
	
	TH1D* histoPythia8Pi02760GeVInvXSection = (TH1D*)histoPythia8InvXSection->Clone();
	
	DrawGammaSetMarker(histoPythia8Pi02760GeVInvXSection, 24, 1.5, kRed+2 , kRed+2);  
        histoPythia8Pi02760GeVInvXSection->SetLineWidth(widthCommonFit);
	
	
	DrawGammaSetMarkerTF1( fitTsallisInvCrossSectionPi02760GeV, 		styleFitCommonSpectrum, widthCommonFit, colorCommonSpectrumPi02760GeV);
	DrawGammaSetMarkerTF1( fitBylinkinDalitzPi0InvCrossSection2760GeV, 	styleFitCommonSpectrum, widthCommonFit, colorCommonSpectrumPi02760GeV);
	

	
	//----------------------- Pi0 7 TeV NLO mu = pt/2 ------------------------------
	TGraph* graphNLOMuHalfPi07TeV = (TGraph*)graphNLOCalcMuHalf7TeV->Clone();
	DrawGammaSetMarkerTGraph(graphNLOMuHalfPi07TeV, styleMarkerNLOMuHalf, sizeMarkerNLO, colorNLOPi07TeVMuHalf, colorNLOPi07TeVMuHalf );
	DrawGammaNLOTGraph( graphNLOMuHalfPi07TeV, widthCommonFit, styleLineNLOMuHalf, colorNLOPi07TeVMuHalf);
	
	//------------------------- Pi0 7 TeV NLO mu = pt ----------------------------
	TGraph* graphNLOMuOnePi07TeV = (TGraph*)graphNLOCalcMuOne7TeV->Clone();
	DrawGammaSetMarkerTGraph(graphNLOMuOnePi07TeV, styleMarkerNLOMuOne, sizeMarkerNLO, colorNLOPi07TeVMuOne, colorNLOPi07TeVMuOne);
	DrawGammaNLOTGraph( graphNLOMuOnePi07TeV, 0.5*widthCommonFit, styleLineNLOMuOne, colorNLOPi07TeVMuOne);
	
	//------------------------- Pi0 7 TeV NLO mu = 2pt -----------------------------
	TGraph* graphNLOMuTwoPi07TeV = (TGraph*)graphNLOCalcMuTwo7TeV->Clone();
	graphNLOMuTwoPi07TeV->SetLineWidth(widthCommonFit*1.5);
	DrawGammaSetMarkerTGraph(graphNLOMuTwoPi07TeV, styleMarkerNLOMuTwo, sizeMarkerNLO, colorNLOPi07TeVMuTwo, colorNLOPi07TeVMuTwo);
	DrawGammaNLOTGraph( graphNLOMuTwoPi07TeV, widthCommonFit, styleLineNLOMuTwo, colorNLOPi07TeVMuTwo);
	
		
	DrawGammaSetMarkerTGraph(graphNLOBKKCalcMuTwo7TeV, styleMarkerNLOMuTwo, sizeMarkerNLO, colorNLOBKKPi07TeVMuTwo, colorNLOBKKPi07TeVMuTwo);
	
	
	DrawGammaSetMarkerTF1( fitTsallisInvCrossSectionPi07TeV,        styleFitCommonSpectrum, widthCommonFit, colorCommonSpectrumPi07TeV);
	DrawGammaSetMarkerTF1( fitBylinkinDalitzPi0InvCrossSection7TeV, styleFitCommonSpectrum, widthCommonFit, colorCommonSpectrumPi07TeV);
	
	//-------------------------Pi0 7 TeV NLO  mu = pt --------------------------------
	
	//TGraph*	graphNLODSS14MuOne7TeV	=	(TGraph*)graphNLODSS14CalcMuOne7TeV->Clone();
	//DrawGammaSetMarkerTGraph(graphNLODSS14MuOne7TeV, styleMarkerNLODSS14MuOne, sizeMarkerNLO, colorNLODSS14Pi07TeVMuOne, colorNLODSS14Pi07TeVMuOne );
	//DrawGammaNLOTGraph( graphNLODSS14MuOne7TeV, widthCommonFit, styleLineNLODSS14MuOne, colorNLODSS14Pi07TeVMuOne);
	
	TGraphAsymmErrors* graphNLODSS14MuOne7TeV = (TGraphAsymmErrors*)graphNLODSS14CalcMuOne7TeV->Clone();
	graphNLODSS14MuOne7TeV->SetLineWidth(widthCommonFit);
	graphNLODSS14MuOne7TeV->SetLineColor(colorNLODSS14Pi02760GeVMuOne);
	graphNLODSS14MuOne7TeV->SetLineStyle(styleLineNLODSS14MuOne);
	
	
	
	
	
	
	fitHistogramTsallisInvCrossSectionPi02760GeV->Scale(1e-1);
	SetStyleHisto(fitHistogramTsallisInvCrossSectionPi02760GeV, widthCommonFit, styleFitCommonSpectrum, colorCommonSpectrumPi02760GeV);
	
	fitHistogramBylinkinInvCrossSection2760GeV->Scale(1e-1);
	SetStyleHisto(fitHistogramBylinkinInvCrossSection2760GeV, widthCommonFit, styleFitCommonSpectrum, colorCommonSpectrumPi02760GeV);
	
	histoPythia8Pi02760GeVInvXSection->Scale(1e-1);
	
	
	
	
	
	
	//****************************** Definition of the Legend Ratio ******************************************
	//**************** Row def ************************
	//**************** Row def ************************
	
	Double_t rowsNLOLegendPrelAndFinal[8] = {0.9,0.81,0.70,0.59,0.48,0.36,0.23,0.11};
	
	//*************** Label sizes *********************
	
 	Double_t textSizeLeftLabelsPrelAndFinal = 0.11;
	Double_t textSizeTopLablesPrelAndFinal = 0.115;
 	Double_t textSizeTopLowerLablesPrelAndFinal = 0.11;
	
	//*************** Column def ***********************
	
	Double_t columnsNLOLegendPrelAndFinal[4] = {0.,0.35,0.61,0.89};//{0.,0.23,0.46,0.74};
	
	//*************** Size factors ********************
	
 	Double_t scaleMarkerNLO = 1.0;
 	Double_t scaleLineWidthNLO = 1.;
	
	//*************** Offsets *************************
	
	Double_t offsetNLOLegendPrelAndFinalMarkerX = 0.09;
	Double_t offsetNLOLegendPrelAndFinalMarkerY = 0.03;
	Double_t offsetNLOLegendPrelAndFinalBox = 0.05;
	Double_t offsetNLOLegendPrelAndFinalLine = 0.06;
		
	Double_t rowsNLOLegendPrelAndFinalOnlyRatioPi0[8] = {0.9,0.81,0.70,0.57,0.46,0.34,0.21,0.09}; // {0.88,0.75,0.62,0.49,0.36,0.25,0.12}; //{0.9,0.79,0.69,0.54,0.41,0.26,0.12,0.0};
	Double_t columnsNLOLegendPrelAndFinalOnlyRatio[4] = {0.,0.35,0.61,0.89};
	
	Double_t textSizeLeftLabelsPrelAndFinalRatio = 0.105;
	Double_t textSizeTopLablesPrelAndFinalRatio = 0.11;
 	Double_t textSizeTopLowerLablesPrelAndFinalRatio = 0.105;
	
	
		
	//*************** first Column **********************************************************
	TLatex *textSpectrumALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[0],rowsNLOLegendPrelAndFinal[2],"#pi^{0}");
	SetStyleTLatex( textSpectrumALLEnergies, textSizeLeftLabelsPrelAndFinal,4);
	
	TLatex *textTsallisFitALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[0],rowsNLOLegendPrelAndFinal[3],"Tsallis fit");
	SetStyleTLatex( textTsallisFitALLEnergies, textSizeLeftLabelsPrelAndFinal,4);
	
	TLatex *textBylinkinFitALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[0],rowsNLOLegendPrelAndFinal[3],"Bylinkin-Rostovtset fit");
	SetStyleTLatex( textBylinkinFitALLEnergies, textSizeLeftLabelsPrelAndFinal,4);
	
	//TLatex *textBylinkinFitALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[0],rowsNLOLegendPrelAndFinal[4],"Bylinkin fit");
	//SetStyleTLatex( textBylinkinFitALLEnergies, textSizeLeftLabelsPrelAndFinal,4);
	
	TLatex *textNLOMuHalfALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[0],rowsNLOLegendPrelAndFinal[4],"NLO #mu= 0.5 #it{p}_{T}");
	SetStyleTLatex( textNLOMuHalfALLEnergies, textSizeLeftLabelsPrelAndFinal,4);
	TLatex *textNLOMuOneALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[0],rowsNLOLegendPrelAndFinal[5],"NLO #mu= #it{p}_{T}");
	SetStyleTLatex( textNLOMuOneALLEnergies, textSizeLeftLabelsPrelAndFinal,4);
	TLatex *textNLOMuTwoALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[0],rowsNLOLegendPrelAndFinal[6],"NLO #mu= 2 #it{p}_{T}");
	SetStyleTLatex( textNLOMuTwoALLEnergies, textSizeLeftLabelsPrelAndFinal,4);
	//TLatex *textNLOMuTwoBKKALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[0],rowsNLOLegendPrelAndFinal[7],"NLO #mu= 2 #it{p}_{T} (BKK)");
	//SetStyleTLatex( textNLOMuTwoBKKALLEnergies, textSizeLeftLabelsPrelAndFinal,4);
	
	//*************** second Column **********************************************************
	TLatex *textPi07TeVNLOALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[1],rowsNLOLegendPrelAndFinal[0],"#pi^{0}, #sqrt{#it{s}} = 7 TeV (*)");
	SetStyleTLatex( textPi07TeVNLOALLEnergies, textSizeTopLablesPrelAndFinal,4);
	
	TLatex *textPi07TeVNLOsysALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[1],rowsNLOLegendPrelAndFinal[1],"syst. + stat.");
	SetStyleTLatex( textPi07TeVNLOsysALLEnergies, textSizeTopLowerLablesPrelAndFinal,4);
	
	TBox* boxCombinedPi07TeVALLEnergies = CreateBoxFromGraph(graphDalitzPi0InvCrossSectionSys7TeV,columnsNLOLegendPrelAndFinal[1]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerY- offsetNLOLegendPrelAndFinalBox, columnsNLOLegendPrelAndFinal[1]+offsetNLOLegendPrelAndFinalMarkerX+offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerY+offsetNLOLegendPrelAndFinalBox);
	
	TMarker* markerCombinedPi07TeVALLEnergies = CreateMarkerFromGraph(graphDalitzPi0InvCrossSectionSys7TeV,columnsNLOLegendPrelAndFinal[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY ,scaleMarkerNLO);
	
	
	TLine * lineTsallisFit7TeVNLOALLEnergies 	= CreateLineFromFit(fitTsallisInvCrossSectionPi07TeV, 		columnsNLOLegendPrelAndFinal[1]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[1]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO);
	//TLine * lineBylinkinFit7TeVNLOALLEnergies 	= CreateLineFromFit(fitBylinkinDalitzPi0InvCrossSection7TeV, 	columnsNLOLegendPrelAndFinal[1]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[4]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[1]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinal[4]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO);
	TLine * lineNLOPi07TeVMuHalfALLEnergies 	= CreateLineFromGraph(graphNLOMuHalfPi07TeV, 			columnsNLOLegendPrelAndFinal[1]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[4]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[1]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinal[4]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO); 
	TLine * lineNLOPi07TeVMuOneALLEnergies 		= CreateLineFromGraph(graphNLOMuOnePi07TeV,			columnsNLOLegendPrelAndFinal[1]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[5]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[1]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinal[5]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO); 
	TLine * lineNLOPi07TeVMuTwoALLEnergies 		= CreateLineFromGraph(graphNLOMuTwoPi07TeV, 			columnsNLOLegendPrelAndFinal[1]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[6]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[1]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinal[6]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO); 
	
	//*************** third Column **********************************************************
	TLatex *textPi02760GeVNLOALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[2],rowsNLOLegendPrelAndFinal[0],"#pi^{0}, #sqrt{#it{s}} = 2.76 TeV (**)");
	SetStyleTLatex( textPi02760GeVNLOALLEnergies, textSizeTopLablesPrelAndFinal,4);
	
	TLatex *textPi02760GeVNLOsysALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[2],rowsNLOLegendPrelAndFinal[1],"syst. + stat.");
	SetStyleTLatex( textPi02760GeVNLOsysALLEnergies, textSizeTopLowerLablesPrelAndFinal,4);
	
	TBox* boxCombinedPi02760GeVALLEnergies = CreateBoxFromGraph(graphDalitzPi0InvCrossSectionSys2760GeVScale,columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerY- offsetNLOLegendPrelAndFinalBox, columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX+offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerY+offsetNLOLegendPrelAndFinalBox);
	TMarker* markerCombinedPi02760GeVALLEnergies = CreateMarkerFromGraph(graphDalitzPi0InvCrossSectionSys2760GeVScale,columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY ,scaleMarkerNLO);
	
	TLine * lineTsallisFitPi02760GeVNLOALLEnergies 	=	CreateLineFromFit(fitTsallisInvCrossSectionPi02760GeV, 		columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinal[3]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO);
	TLine * lineNLOPi02760GeVMuHalfALLEnergies 	= 	CreateLineFromGraph(graphNLOMuHalfPi02760GeV,  			columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[4]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinal[4]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO); 
	TLine * lineNLOPi02760GeVMuOneALLEnergies 	= 	CreateLineFromGraph(graphNLOMuOnePi02760GeV, 			columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[5]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinal[5]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO); 
	TLine * lineNLOPi02760GeVMuTwoALLEnergies 	= 	CreateLineFromGraph(graphNLOMuTwoPi02760GeV, 			columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[6]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinal[6]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO); 
	TLine * linePythia8ALLEnergies 			= 	CreateLineFromHisto(histoPythia8Pi02760GeVInvXSection, 		columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinal[7]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinal[2]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinal[7]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO); 
	
	
	
	//*************** first Column **********************************************************
	TLatex *textFitTsallisSpectrumOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[0],rowsNLOLegendPrelAndFinalOnlyRatioPi0[2],"#pi^{0} Tsallis fit");
	SetStyleTLatex( textFitTsallisSpectrumOnlyRatioPi0, textSizeLeftLabelsPrelAndFinalRatio,4);
	
	TLatex *textFitBylinkinSpectrumOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[0],rowsNLOLegendPrelAndFinalOnlyRatioPi0[2],"#pi^{0} Bylinkin-Rostovtset fit");
	SetStyleTLatex( textFitBylinkinSpectrumOnlyRatioPi0, textSizeLeftLabelsPrelAndFinalRatio,4);
	
			
	TLatex *textNLOMuHalfOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[0],rowsNLOLegendPrelAndFinalOnlyRatioPi0[3],"NLO #mu= 0.5 #it{p}_{T}");
	SetStyleTLatex( textNLOMuHalfOnlyRatioPi0, textSizeLeftLabelsPrelAndFinalRatio,4);
	
	TLatex *textNLOMuOneOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[0],rowsNLOLegendPrelAndFinalOnlyRatioPi0[4],"NLO #mu= #it{p}_{T}");
	SetStyleTLatex( textNLOMuOneOnlyRatioPi0, textSizeLeftLabelsPrelAndFinalRatio,4);
	
	TLatex *textNLOMuTwoOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[0],rowsNLOLegendPrelAndFinalOnlyRatioPi0[5],"NLO #mu= 2 #it{p}_{T}");
	SetStyleTLatex( textNLOMuTwoOnlyRatioPi0, textSizeLeftLabelsPrelAndFinalRatio,4);
	
	TLatex *textNLOBKKMuTwoOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[0],rowsNLOLegendPrelAndFinalOnlyRatioPi0[6],"NLO #mu= 2 #it{p}_{T} (BKK)");
	SetStyleTLatex( textNLOBKKMuTwoOnlyRatioPi0, textSizeLeftLabelsPrelAndFinalRatio,4);
	
	
	//*************** second Column **********************************************************
	TLatex *textPi07TeVNLOOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[1],rowsNLOLegendPrelAndFinalOnlyRatioPi0[0],"#pi^{0}, #sqrt{#it{s}} = 7 TeV (*)");
	SetStyleTLatex( textPi07TeVNLOOnlyRatioPi0, textSizeTopLablesPrelAndFinalRatio,4);
	TLatex *textPi07TeVNLOsysOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[1],rowsNLOLegendPrelAndFinalOnlyRatioPi0[1],"syst. + stat.");
	SetStyleTLatex( textPi07TeVNLOsysOnlyRatioPi0, textSizeTopLowerLablesPrelAndFinalRatio,4);
	
	
	TBox* boxFitTsallisPi07TeVOnlyRatioPi0 = CreateBoxFromGraph(graphDalitzPi0InvCrossSectionSys7TeV,columnsNLOLegendPrelAndFinalOnlyRatio[1]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+ offsetNLOLegendPrelAndFinalMarkerY- offsetNLOLegendPrelAndFinalBox, columnsNLOLegendPrelAndFinalOnlyRatio[1]+offsetNLOLegendPrelAndFinalMarkerX+offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+ offsetNLOLegendPrelAndFinalMarkerY+offsetNLOLegendPrelAndFinalBox);
	TMarker* markerCombinedPi07TeVOnlyRatioPi0 = CreateMarkerFromGraph(graphDalitzPi0InvCrossSectionSys7TeV,columnsNLOLegendPrelAndFinalOnlyRatio[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY ,scaleMarkerNLO);
	
	//TBox* boxFitBylinkinPi07TeVOnlyRatioPi0 = CreateBoxFromGraph(graphDalitzPi0InvCrossSectionCompl7TeV,columnsNLOLegendPrelAndFinalOnlyRatio[1]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[3]+ offsetNLOLegendPrelAndFinalMarkerY- offsetNLOLegendPrelAndFinalBox, columnsNLOLegendPrelAndFinalOnlyRatio[1]+offsetNLOLegendPrelAndFinalMarkerX+offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[3]+ offsetNLOLegendPrelAndFinalMarkerY+offsetNLOLegendPrelAndFinalBox);
	
		
	TLine * lineNLOPi07TeVMuHalfOnlyRatioPi0 = CreateLineFromGraph(graphNLOMuHalfPi07TeV, 	columnsNLOLegendPrelAndFinalOnlyRatio[1]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[3]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinalOnlyRatio[1]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[3]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO*1.5); 
	TLine * lineNLOPi07TeVMuOneOnlyRatioPi0  = CreateLineFromGraph(graphNLOMuOnePi07TeV,	columnsNLOLegendPrelAndFinalOnlyRatio[1]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[4]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinalOnlyRatio[1]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[4]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO*1.5); 
	TLine * lineNLOPi07TeVMuTwoOnlyRatioPi0  = CreateLineFromGraph(graphNLOMuTwoPi07TeV, 	columnsNLOLegendPrelAndFinalOnlyRatio[1]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[5]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinalOnlyRatio[1]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[5]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO*1.5); 
	
	
	//*************** third Column **********************************************************
	TLatex *textPi02760GeVNLOOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[2],rowsNLOLegendPrelAndFinalOnlyRatioPi0[0],"#pi^{0}, #sqrt{#it{s}} = 2.76 TeV (**)");
	SetStyleTLatex( textPi02760GeVNLOOnlyRatioPi0, textSizeTopLablesPrelAndFinalRatio,4);
	TLatex *textPi02760GeVNLOsysOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[2],rowsNLOLegendPrelAndFinalOnlyRatioPi0[1],"syst. + stat.");
	SetStyleTLatex( textPi02760GeVNLOsysOnlyRatioPi0, textSizeTopLowerLablesPrelAndFinalRatio,4);
	
	TBox* boxTSallisFitPi02760GeVOnlyRatioPi0 = CreateBoxFromGraph(graphDalitzPi0InvCrossSectionSys2760GeVScale,columnsNLOLegendPrelAndFinalOnlyRatio[2]+offsetNLOLegendPrelAndFinalMarkerX-offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+ offsetNLOLegendPrelAndFinalMarkerY- offsetNLOLegendPrelAndFinalBox, columnsNLOLegendPrelAndFinalOnlyRatio[2]+offsetNLOLegendPrelAndFinalMarkerX+offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+ offsetNLOLegendPrelAndFinalMarkerY+offsetNLOLegendPrelAndFinalBox);
	TMarker* markerCombinedPi02760GeVOnlyRatioPi0 = CreateMarkerFromGraph(graphDalitzPi0InvCrossSectionSys2760GeVScale,columnsNLOLegendPrelAndFinalOnlyRatio[2]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY ,scaleMarkerNLO);
	
	
	
	
	
	TLine * lineNLOPi02760GeVMuHalfOnlyRatioPi0 = CreateLineFromGraph(graphNLOMuHalfPi02760GeV,  	columnsNLOLegendPrelAndFinalOnlyRatio[2]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[3]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinalOnlyRatio[2]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinalOnlyRatioPi0[3]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO*1.5); 
	TLine * lineNLOPi02760GeVMuOneOnlyRatioPi0  = CreateLineFromGraph(graphNLOMuOnePi02760GeV,  	columnsNLOLegendPrelAndFinalOnlyRatio[2]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[4]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinalOnlyRatio[2]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinalOnlyRatioPi0[4]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO*1.5); 
	TLine * lineNLOPi02760GeVMuTwoOnlyRatioPi0  = CreateLineFromGraph(graphNLOMuTwoPi02760GeV,	columnsNLOLegendPrelAndFinalOnlyRatio[2]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[5]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinalOnlyRatio[2]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinalOnlyRatioPi0[5]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO*1.5); 
	TLine * linePythi82760GeVOnlyRationPi0	    = CreateLineFromHisto(histoPythia8Pi02760GeVInvXSection      ,columnsNLOLegendPrelAndFinalOnlyRatio[2]+ offsetNLOLegendPrelAndFinalMarkerX- offsetNLOLegendPrelAndFinalLine, rowsNLOLegendPrelAndFinalOnlyRatioPi0[6]+offsetNLOLegendPrelAndFinalMarkerY, columnsNLOLegendPrelAndFinalOnlyRatio[2]+ offsetNLOLegendPrelAndFinalMarkerX+ offsetNLOLegendPrelAndFinalLine,rowsNLOLegendPrelAndFinalOnlyRatioPi0[6]+offsetNLOLegendPrelAndFinalMarkerY, scaleLineWidthNLO*1.5);  
	
	
	///////////////////////////////////////DSS14 label////////////////////////////////////7
	
	//*************** first Column **********************************************************
	
	TLatex *textNLODSS14MuHalfALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[0],rowsNLOLegendPrelAndFinal[4],"NLO DSS14 #mu= 0.5 #it{p}_{T}");
	SetStyleTLatex( textNLODSS14MuHalfALLEnergies, textSizeLeftLabelsPrelAndFinal,4);
	TLatex *textNLODSS14MuOneALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[0],rowsNLOLegendPrelAndFinal[5],"NLO DSS14 #mu= #it{p}_{T}");
	SetStyleTLatex( textNLODSS14MuOneALLEnergies, textSizeLeftLabelsPrelAndFinal,4);
	TLatex *textNLODSS14MuTwoALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[0],rowsNLOLegendPrelAndFinal[6],"NLO DSS14 #mu= 2 #it{p}_{T}");
	SetStyleTLatex( textNLODSS14MuTwoALLEnergies, textSizeLeftLabelsPrelAndFinal,4);
	TLatex *textPythia8ALLEnergies = new TLatex(columnsNLOLegendPrelAndFinal[0],rowsNLOLegendPrelAndFinal[7],"Pythia 8, Tune 4C");
	SetStyleTLatex( textPythia8ALLEnergies, textSizeLeftLabelsPrelAndFinal,4);
	
	
	
	
	
	

	//************************************************* End Legend ***************************************************
	
	
	
	//**************************** pQCD calculation ****************************//
	
	
	
	TCanvas* canvasInvXSectionNLOOnlySpectraPi0 = new TCanvas("canvasInvXSectionNLOOnlySpectraPi0","",200,10,1200,1200);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionNLOOnlySpectraPi0,  0.15, 0.02, 0.03, 0.06);
	
	
	TPad* padComparisonXSectionNLOOnlySpectraPi0 = new TPad("padComparisonXSectionNLOOnlySpectraPi0", "", 0., 0., 1., 0.77,-1, -1, -2);
	DrawGammaPadSettings( padComparisonXSectionNLOOnlySpectraPi0, 0.15, 0.02, 0.02, 0.09);
	padComparisonXSectionNLOOnlySpectraPi0->Draw();
	
	TPad* padXSectionNLOLegendOnlySpectraPi0 = new TPad("padXSectionNLOLegendOnlySpectraPi0", "", 0.15, 0.77, 0.97, 1.,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOLegendOnlySpectraPi0, 0.15, 0.02, 0.03, 0.);
	padXSectionNLOLegendOnlySpectraPi0->Draw();
	
	padComparisonXSectionNLOOnlySpectraPi0->cd();
	padComparisonXSectionNLOOnlySpectraPi0->SetLogy();		
	padComparisonXSectionNLOOnlySpectraPi0->SetLogx();		
	
	
	
	TH2F * histo2DInvXSectionNLOOnlySpectraPi0 = new TH2F("histo2DInvXSectionNLOOnlySpectraPi0","histo2DInvXSectionNLOOnlySpectraPi0",1000,0.23,30.,1000,2e0,10e12);
	SetStyleHistoTH2ForGraphs(histo2DInvXSectionNLOOnlySpectraPi0, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
	histo2DInvXSectionNLOOnlySpectraPi0->DrawCopy(); 
	
	graphDalitzPi0InvCrossSectionSys7TeV->Draw("2,same");
	graphDalitzPi0InvCrossSectionStat7TeV->Draw("p,same");
	
	graphDalitzPi0InvCrossSectionSys2760GeVScale->Draw("2,same");
	graphDalitzPi0InvCrossSectionStat2760GeVScale->Draw("p,same");
	graphNLOMuHalfPi02760GeV->Draw("same,l");
	graphNLOMuOnePi02760GeV->Draw("same,c");
	graphNLOMuTwoPi02760GeV->Draw("same,c");
	
	graphNLOMuHalfPi07TeV->Draw("same,c");
	graphNLOMuOnePi07TeV->Draw("same,c");
	graphNLOMuTwoPi07TeV->Draw("same,c");
	
	
	fitTsallisInvCrossSectionPi07TeV->Draw("same");
	fitHistogramTsallisInvCrossSectionPi02760GeV->Draw("same,c");
	
	
	TLatex *labelScalingPi07TeV = new TLatex(0.47,2E10,"x 1");
	SetStyleTLatex( labelScalingPi07TeV, 0.025,4,graphDalitzPi0InvCrossSectionSys7TeV->GetLineColor(),62,kFALSE);
	labelScalingPi07TeV->Draw();
	
	labelScalingPi07TeV->Draw();
	TLatex *labelScalingPi02760GeV = new TLatex(0.47,2E9,"x 10^{-1}");
	SetStyleTLatex( labelScalingPi02760GeV, 0.025,4,graphDalitzPi0InvCrossSectionSys2760GeVScale->GetLineColor(),62,kFALSE);
	labelScalingPi02760GeV->Draw();
	
	legendPlotLabel->Draw();
	
	
	padXSectionNLOLegendOnlySpectraPi0->cd();
	DrawGammaPadSettings( padXSectionNLOLegendOnlySpectraPi0, 0., 0., 0., 0.);
	padXSectionNLOLegendOnlySpectraPi0->SetBorderMode(-1);
	padXSectionNLOLegendOnlySpectraPi0->SetBorderSize(3);
	padXSectionNLOLegendOnlySpectraPi0->Draw();
	padXSectionNLOLegendOnlySpectraPi0->cd();

	//*************** first Column **********************************************************	
	textSpectrumALLEnergies->Draw();
	textTsallisFitALLEnergies->Draw();
	//textBylinkinFitALLEnergies->Draw();
	textNLOMuHalfALLEnergies->Draw();
	textNLOMuOneALLEnergies->Draw();
	textNLOMuTwoALLEnergies->Draw();
	
	//*************** second Column **********************************************************
	textPi07TeVNLOALLEnergies->Draw();
	textPi07TeVNLOsysALLEnergies->Draw();
	boxCombinedPi07TeVALLEnergies->Draw("l");
	markerCombinedPi07TeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineTsallisFit7TeVNLOALLEnergies->Draw("same");
	lineNLOPi07TeVMuHalfALLEnergies->Draw("same");
	lineNLOPi07TeVMuOneALLEnergies->Draw("same");
	lineNLOPi07TeVMuTwoALLEnergies->Draw("same");
	//lineNLOPi07TeVMuTwoBKKALLEnergies->Draw("same");
	
	//*************** third Column **********************************************************
	textPi02760GeVNLOALLEnergies->Draw();
	textPi02760GeVNLOsysALLEnergies->Draw();
	boxCombinedPi02760GeVALLEnergies->Draw("l");
	markerCombinedPi02760GeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
	lineTsallisFitPi02760GeVNLOALLEnergies->Draw("same");
	//lineBylinkinFitPi02760GeVNLOALLEnergies->Draw("same");
	lineNLOPi02760GeVMuHalfALLEnergies->Draw("same");
	lineNLOPi02760GeVMuOneALLEnergies->Draw("same");
	lineNLOPi02760GeVMuTwoALLEnergies->Draw("same");
	linePythia8ALLEnergies->Draw("same");
		
    
	       
	
	canvasInvXSectionNLOOnlySpectraPi0->Update();
	canvasInvXSectionNLOOnlySpectraPi0->Print(Form("%s/InvariantCrossSecPP_pQCD_Old.%s",outputDir.Data(),suffix.Data()));
	
	graphRatioInvCrossSectionPi07TeVTsallisFit = (TGraphAsymmErrors*)graphDalitzPi0InvCrossSectionCompl7TeV->Clone();
	graphRatioInvCrossSectionPi07TeVTsallisFit = CalculateGraphErrRatioToFit(graphRatioInvCrossSectionPi07TeVTsallisFit, fitTsallisInvCrossSectionPi07TeV); 
	
	
	graphRatioInvCrossSectionPi07TeVBylinkinFit 	=	(TGraphAsymmErrors*)graphDalitzPi0InvCrossSectionCompl7TeV->Clone();
	graphRatioInvCrossSectionPi07TeVBylinkinFit	= 	CalculateGraphErrRatioToFit(graphRatioInvCrossSectionPi07TeVBylinkinFit,fitBylinkinDalitzPi0InvCrossSection7TeV);
	
	
	DrawGammaSetMarkerTGraphAsym(graphRatioInvCrossSectionPi07TeVTsallisFit, markerStyleCommmonSpectrumPi07TeV,markerSizeCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes, kTRUE);
	graphRatioInvCrossSectionPi07TeVTsallisFit->SetLineWidth(widthCommonErrors);
	
	DrawGammaSetMarkerTGraphAsym(graphRatioInvCrossSectionPi07TeVBylinkinFit, markerStyleCommmonSpectrumPi07TeV,markerSizeCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, colorCommonSpectrumPi07TeV, widthCommonSpectrumBoxes,kTRUE);
	graphRatioInvCrossSectionPi07TeVBylinkinFit->SetLineWidth(widthCommonErrors);
	
	graphRatioInvCrossSectionPi02760GeVTsallisFit 	= (TGraphAsymmErrors*)graphDalitzPi0InvCrossSectionCompl2760GeV->Clone();
	graphRatioInvCrossSectionPi02760GeVTsallisFit 	= CalculateGraphErrRatioToFit(graphRatioInvCrossSectionPi02760GeVTsallisFit, fitTsallisInvCrossSectionPi02760GeV); 
	
	graphRatioInvCrossSectionPi02760GeVTBylinkinFit = 	(TGraphAsymmErrors*)graphDalitzPi0InvCrossSectionCompl2760GeV->Clone();
	graphRatioInvCrossSectionPi02760GeVTBylinkinFit	= 	CalculateGraphErrRatioToFit(graphRatioInvCrossSectionPi02760GeVTBylinkinFit,fitBylinkinDalitzPi0InvCrossSection2760GeV);
	
	
	DrawGammaSetMarkerTGraphAsym(graphRatioInvCrossSectionPi02760GeVTsallisFit, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes, kTRUE);
	graphRatioInvCrossSectionPi02760GeVTsallisFit->SetLineWidth(widthCommonErrors);
	
	DrawGammaSetMarkerTGraphAsym(graphRatioInvCrossSectionPi02760GeVTBylinkinFit, markerStyleCommmonSpectrum2760GeV,1.6*markerSizeCommonSpectrum, colorCommonSpectrumPi02760GeV, colorCommonSpectrumPi02760GeV, widthCommonSpectrumBoxes);
	graphRatioInvCrossSectionPi02760GeVTBylinkinFit->SetLineWidth(widthCommonErrors);
	
		
			
	graphRatioNLODalitzPi02760GeVMuHalf   		= CalculateGraphRatioToFit (graphNLOCalcMuHalf2760GeV, 	fitTsallisInvCrossSectionPi02760GeV); 
	graphRatioNLODalitzPi02760GeVMuOne    		= CalculateGraphRatioToFit (graphNLOCalcMuOne2760GeV, 	fitTsallisInvCrossSectionPi02760GeV); 
	graphRatioNLODalitzPi02760GeVMuTwo    		= CalculateGraphRatioToFit (graphNLOCalcMuTwo2760GeV, 	fitTsallisInvCrossSectionPi02760GeV);
	
	 TH1D* histoRatioPythia8ToFitTsallis                     = (TH1D*) histoPythia8InvXSection->Clone(); 
      
	histoRatioPythia8ToFitTsallis                          = CalculateHistoRatioToFit (histoRatioPythia8ToFitTsallis, fitTsallisInvCrossSectionPi02760GeV); 
   	graphRatioNLODSS14DalitzPi02760GeVMuHalf	= CalculateGraphErrRatioToFit (graphNLOCalcmuHalfDSS14InvSecPi02760GeV, fitTsallisInvCrossSectionPi02760GeV);
	graphRatioNLODSS14DalitzPi02760GeVMuOne		= CalculateGraphErrRatioToFit (graphNLOCalcmuOneDSS14InvSecPi02760GeV, fitTsallisInvCrossSectionPi02760GeV);
	graphRatioNLODSS14DalitzPi02760GeVMuTwo		= CalculateGraphErrRatioToFit (graphNLOCalcmuTwoDSS14InvSecPi02760GeV, fitTsallisInvCrossSectionPi02760GeV);
	
	
	TH1D* histoRatioPythia8ToFitBylinkin                     = (TH1D*) histoPythia8InvXSection->Clone(); 
      
	histoRatioPythia8ToFitBylinkin                         = CalculateHistoRatioToFit (histoRatioPythia8ToFitBylinkin, fitBylinkinDalitzPi0InvCrossSection2760GeV); 
   	
	
	graphRatioBylinkinFitNLODSS14DalitzPi02760GeVMuHalf		= CalculateGraphErrRatioToFit (graphNLOCalcmuHalfDSS14InvSecPi02760GeV, fitBylinkinDalitzPi0InvCrossSection2760GeV);
	graphRatioBylinkinFitNLODSS14DalitzPi02760GeVMuOne		= CalculateGraphErrRatioToFit (graphNLOCalcmuOneDSS14InvSecPi02760GeV,  fitBylinkinDalitzPi0InvCrossSection2760GeV);
	graphRatioBylinkinFitNLODSS14DalitzPi02760GeVMuTwo		= CalculateGraphErrRatioToFit (graphNLOCalcmuTwoDSS14InvSecPi02760GeV,  fitBylinkinDalitzPi0InvCrossSection2760GeV);
	
	
	
	
	
	
	
	
	
	graphRatioNLODalitzPi07TeVMuHalf   	= CalculateGraphRatioToFit (graphNLOCalcMuHalf7TeV, 	fitTsallisInvCrossSectionPi07TeV); 
	graphRatioNLODalitzPi07TeVMuOne    	= CalculateGraphRatioToFit (graphNLOCalcMuOne7TeV, 	fitTsallisInvCrossSectionPi07TeV); 
	graphRatioNLODalitzPi07TeVMuTwo    	= CalculateGraphRatioToFit (graphNLOCalcMuTwo7TeV, 	fitTsallisInvCrossSectionPi07TeV); 
	graphRatioNLODSS14DalitzPi07TeVMuOne	= CalculateGraphErrRatioToFit (graphNLODSS14MuOne7TeV,     fitTsallisInvCrossSectionPi07TeV);
	
	
	DrawGammaNLOTGraph( graphRatioNLODalitzPi02760GeVMuHalf, 0.5*widthCommonFit, styleLineNLOMuHalf, 	colorNLOPi02760GeVMuHalf);
	DrawGammaNLOTGraph( graphRatioNLODalitzPi02760GeVMuOne,  0.5*widthCommonFit, styleLineNLOMuOne, 	colorNLOPi02760GeVMuOne);
	DrawGammaNLOTGraph( graphRatioNLODalitzPi02760GeVMuTwo,  0.5*widthCommonFit, styleLineNLOMuTwo, 	colorNLOPi02760GeVMuTwo);
	
	//////////////Might be temporal
	
	DrawGammaNLOTGraph( graphRatioNLODSS14DalitzPi02760GeVMuHalf, 0.5*widthCommonFit, styleLineNLODSS14MuHalf, 	colorNLODSS14Pi02760GeVMuHalf);
	DrawGammaNLOTGraph( graphRatioNLODSS14DalitzPi02760GeVMuOne,  0.5*widthCommonFit, styleLineNLODSS14MuOne, 	colorNLODSS14Pi02760GeVMuOne);
	DrawGammaNLOTGraph( graphRatioNLODSS14DalitzPi02760GeVMuTwo,  0.5*widthCommonFit, styleLineNLODSS14MuTwo, 	colorNLODSS14Pi02760GeVMuTwo);
	
	
	
	DrawGammaNLOTGraph( graphRatioBylinkinFitNLODSS14DalitzPi02760GeVMuHalf, 0.5*widthCommonFit, styleLineNLODSS14MuHalf, 	colorNLODSS14Pi02760GeVMuHalf);
	DrawGammaNLOTGraph( graphRatioBylinkinFitNLODSS14DalitzPi02760GeVMuOne,  0.5*widthCommonFit, styleLineNLODSS14MuOne, 	colorNLODSS14Pi02760GeVMuOne);
	DrawGammaNLOTGraph( graphRatioBylinkinFitNLODSS14DalitzPi02760GeVMuTwo,  0.5*widthCommonFit, styleLineNLODSS14MuTwo, 	colorNLODSS14Pi02760GeVMuTwo);
	
	
	
	
	
	
	DrawGammaNLOTGraph( graphRatioNLODalitzPi07TeVMuHalf, 		0.5*widthCommonFit, styleLineNLOMuHalf, 	colorNLOPi07TeVMuHalf);
	DrawGammaNLOTGraph( graphRatioNLODalitzPi07TeVMuOne,		0.5*widthCommonFit, styleLineNLOMuOne, 		colorNLOPi07TeVMuOne);
	DrawGammaNLOTGraph( graphRatioNLODalitzPi07TeVMuTwo,  		0.5*widthCommonFit, styleLineNLOMuTwo, 		colorNLOPi07TeVMuTwo);
	DrawGammaNLOTGraph( graphRatioNLODSS14DalitzPi07TeVMuOne,	0.5*widthCommonFit, styleLineNLODSS14MuOne,     colorNLODSS14Pi07TeVMuOne);
	
	
	
	TCanvas* canvasInvXSectionNLOOnlyRatioPi0 	= new TCanvas("canvasInvXSectionNLOOnlyRatioPi0","",200,10,1200,834);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionNLOOnlyRatioPi0,  0.15, 0.02, 0.03, 0.06);  
	
	TPad* padXSectionNLOLegendOnlyRatioPi0 		= new TPad("padXSectionNLOLegendOnlyRatioPi0", 	   "", 0.15, 0.66, 0.957,  1.,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOLegendOnlyRatioPi0, 0.15, 0.02, 0.01, 0.0);
	padXSectionNLOLegendOnlyRatioPi0->Draw();

	TPad* padXSectionNLOOnlyRatioPi0Pi07TeV 	= new TPad("padXSectionNLOOnlyRatioPi0Pi07TeV",    "", 0., 0.366, 1., 0.66,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOOnlyRatioPi0Pi07TeV, 0.15, 0.02, 0.02, 0.);
	padXSectionNLOOnlyRatioPi0Pi07TeV->Draw();
	
	TPad* padXSectionNLOOnlyRatioPi0Pi02760GeV	= new TPad("padXSectionNLOOnlyRatioPi0Pi02760GeV", "", 0.,0.0, 1., 0.366,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLOOnlyRatioPi0Pi02760GeV, 0.15, 0.02, 0., 0.28);
	padXSectionNLOOnlyRatioPi0Pi02760GeV->Draw();
	
	
	padXSectionNLOOnlyRatioPi0Pi07TeV->cd();
	padXSectionNLOOnlyRatioPi0Pi07TeV->SetLogx();
	
	TH2F * ratio2DInvXSectionNLOPi0 = new TH2F("ratio2DInvXSectionNLOPi0","ratio2DInvXSectionNLOPi0",1000,0.23,30.,1000,0.4,3.55);
	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionNLOPi0, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 36,36,36,36,2.7,0.8, 512, 505);
	//Jus for differen size of Tpads
	ratio2DInvXSectionNLOPi0->GetXaxis()->SetLabelFont(43);
	ratio2DInvXSectionNLOPi0->GetYaxis()->SetLabelFont(43); 
	ratio2DInvXSectionNLOPi0->GetXaxis()->SetTitleFont(63);
	ratio2DInvXSectionNLOPi0->GetYaxis()->SetTitleFont(63);
	ratio2DInvXSectionNLOPi0->DrawCopy(); 
	
	graphRatioNLODalitzPi07TeVMuHalf->Draw("same,c");
	graphRatioNLODalitzPi07TeVMuOne->Draw("same,c");
	graphRatioNLODalitzPi07TeVMuTwo->Draw("same,c");
	boxErrorSigmaPi07TeVRatio->Draw();
	graphRatioInvCrossSectionPi07TeVTsallisFit->Draw("p,E2same");
	
	labelRatioNLOPi07TeV->Draw();
	
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
	padXSectionNLOOnlyRatioPi0Pi02760GeV->cd();
	padXSectionNLOOnlyRatioPi0Pi02760GeV->SetLogx();
	
	TH2F *  ratio2DInvXSectionNLOPi02760GeV;
	 ratio2DInvXSectionNLOPi02760GeV = new TH2F(" ratio2DInvXSectionNLOPi02760GeV"," ratio2DInvXSectionNLOPi02760GeV",1000,0.23,30.,1000,0.4,3.55);
	SetStyleHistoTH2ForGraphs( ratio2DInvXSectionNLOPi02760GeV, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}",36,36,36,36, 2.7,0.80, 512, 505);
	ratio2DInvXSectionNLOPi02760GeV->GetXaxis()->SetLabelFont(43);
	ratio2DInvXSectionNLOPi02760GeV->GetYaxis()->SetLabelFont(43); 
	ratio2DInvXSectionNLOPi02760GeV->GetXaxis()->SetTitleFont(63);
	ratio2DInvXSectionNLOPi02760GeV->GetYaxis()->SetTitleFont(63);
	ratio2DInvXSectionNLOPi02760GeV->DrawCopy(); 
	
	graphRatioNLODalitzPi02760GeVMuHalf->Draw("same,c");
	graphRatioNLODalitzPi02760GeVMuOne->Draw("same,c");
	graphRatioNLODalitzPi02760GeVMuTwo->Draw("same,c");
	
	boxErrorSigmaPi02760GeVRatio->Draw();
	graphRatioInvCrossSectionPi02760GeVTsallisFit->Draw("p,E2same");
	labelRatioNLOPi02760GeV->Draw();
	
	
	
	DrawGammaLines(0., 30.,1., 1.,0.1);
	
			
	padXSectionNLOLegendOnlyRatioPi0->cd();
	DrawGammaPadSettings( padXSectionNLOLegendOnlyRatioPi0, 0., 0., 0.0, 0.0);
	padXSectionNLOLegendOnlyRatioPi0->SetBorderMode(-1);
	padXSectionNLOLegendOnlyRatioPi0->SetBorderSize(3);
	padXSectionNLOLegendOnlyRatioPi0->Draw();
	padXSectionNLOLegendOnlyRatioPi0->cd();
	
	//*************** first Column **********************************************************

	textFitTsallisSpectrumOnlyRatioPi0->Draw();
	//textFitBylinkinSpectrumOnlyRatioPi0->Draw();
	textNLOMuHalfOnlyRatioPi0->Draw();
	textNLOMuOneOnlyRatioPi0->Draw();
	textNLOMuTwoOnlyRatioPi0->Draw();
	textNLOBKKMuTwoOnlyRatioPi0->Draw();
	
	
//	*************** second Column **********************************************************
	textPi07TeVNLOOnlyRatioPi0->Draw();
	textPi07TeVNLOsysOnlyRatioPi0->Draw();
	boxFitTsallisPi07TeVOnlyRatioPi0->Draw("l");
	markerCombinedPi07TeVOnlyRatioPi0->DrawMarker(columnsNLOLegendPrelAndFinalOnlyRatio[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY);
	//boxFitBylinkinPi07TeVOnlyRatioPi0->Draw("f");
	
	
	//lineNLOPi07TeVMuHalfOnlyRatioPi0->Draw("same");
	lineNLOPi07TeVMuOneOnlyRatioPi0->Draw("same");
	//lineNLOPi07TeVMuTwoOnlyRatioPi0->Draw("same");
	
	
	//*************** third Column **********************************************************
	
	textPi02760GeVNLOOnlyRatioPi0->Draw();
	textPi02760GeVNLOsysOnlyRatioPi0->Draw();
	
	boxTSallisFitPi02760GeVOnlyRatioPi0->Draw("l");
	markerCombinedPi02760GeVOnlyRatioPi0->DrawMarker(columnsNLOLegendPrelAndFinalOnlyRatio[2]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY);
	
	//boxBylinkinFitPi02760GeVOnlyRatioPi0->Draw("f");
	
	
	lineNLOPi02760GeVMuHalfOnlyRatioPi0->Draw("same");
	lineNLOPi02760GeVMuOneOnlyRatioPi0->Draw("same");
	lineNLOPi02760GeVMuTwoOnlyRatioPi0->Draw("same"); 
	linePythi82760GeVOnlyRationPi0->Draw("same");

	//************************************************* End Legend ***************************************************
	
	
	canvasInvXSectionNLOOnlyRatioPi0->Update();
	canvasInvXSectionNLOOnlyRatioPi0->Print(Form("%s/Pi0_InvXSectionNLO_OnlyRatio_Paper.%s",outputDir.Data(),suffix.Data()));
	
	
	//////////////////////////////////////////////////////////DSS14////////////////////////////////////////////
	
	
	//*************** first Ratio Column **********************************************************
	
	TLatex *textNLODSS14MuHalfOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[0],rowsNLOLegendPrelAndFinalOnlyRatioPi0[3],"NLO DSS14 #mu= 0.5 #it{p}_{T}");
	SetStyleTLatex( textNLODSS14MuHalfOnlyRatioPi0, textSizeLeftLabelsPrelAndFinalRatio,4);
	
	TLatex *textNLODSS14MuOneOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[0],rowsNLOLegendPrelAndFinalOnlyRatioPi0[4],"NLO DSS14 #mu= #it{p}_{T}");
	SetStyleTLatex( textNLODSS14MuOneOnlyRatioPi0, textSizeLeftLabelsPrelAndFinalRatio,4);
	
	TLatex *textNLODSS14MuTwoOnlyRatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[0],rowsNLOLegendPrelAndFinalOnlyRatioPi0[5],"NLO DSS14 #mu= 2 #it{p}_{T}");
	SetStyleTLatex( textNLODSS14MuTwoOnlyRatioPi0, textSizeLeftLabelsPrelAndFinalRatio,4);
	
	TLatex *textPythia8RatioPi0 = new TLatex(columnsNLOLegendPrelAndFinalOnlyRatio[0],rowsNLOLegendPrelAndFinalOnlyRatioPi0[6],"Pythia 8, Tune 4C");
	SetStyleTLatex( textPythia8RatioPi0, textSizeLeftLabelsPrelAndFinalRatio,4);
	
	
	
	
	
        TCanvas* canvasInvXSectionNLODSS14OnlySpectraPi0 = new TCanvas("canvasInvXSectionNLODSS14OnlySpectraPi0","",200,10,1200,1200);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionNLODSS14OnlySpectraPi0,  0.15, 0.02, 0.03, 0.06);
	
	
	TPad* padComparisonXSectionNLODSS14OnlySpectraPi0 = new TPad("padComparisonXSectionNLODSS14OnlySpectraPi0", "", 0., 0., 1., 0.77,-1, -1, -2);
	DrawGammaPadSettings( padComparisonXSectionNLODSS14OnlySpectraPi0, 0.15, 0.02, 0.02, 0.09);
	padComparisonXSectionNLODSS14OnlySpectraPi0->Draw();
	
	TPad* padXSectionNLODSS14LegendOnlySpectraPi0 = new TPad("padXSectionNLODSS14LegendOnlySpectraPi0", "", 0.15, 0.77, 0.97, 1.,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLODSS14LegendOnlySpectraPi0, 0.15, 0.02, 0.03, 0.);
	padXSectionNLODSS14LegendOnlySpectraPi0->Draw();
	
	padComparisonXSectionNLODSS14OnlySpectraPi0->cd();
	padComparisonXSectionNLODSS14OnlySpectraPi0->SetLogy();		
	padComparisonXSectionNLODSS14OnlySpectraPi0->SetLogx();		
	
	
	
	TH2F * histo2DInvXSectionNLODSS14OnlySpectraPi0 = new TH2F("histo2DInvXSectionNLODSS14OnlySpectraPi0","histo2DInvXSectionNLODSS14OnlySpectraPi0",1000,0.23,30.,1000,1e2,10e10);
	SetStyleHistoTH2ForGraphs(histo2DInvXSectionNLODSS14OnlySpectraPi0, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
	histo2DInvXSectionNLODSS14OnlySpectraPi0->GetXaxis()->SetRangeUser(0.4,12.0);
	histo2DInvXSectionNLODSS14OnlySpectraPi0->DrawCopy(); 
	
	graphDalitzPi0InvCrossSectionSys7TeV->Draw("2,same");
	graphDalitzPi0InvCrossSectionStat7TeV->Draw("p,same");
	
	graphDalitzPi0InvCrossSectionSys2760GeVScale->Draw("2,same");
	graphDalitzPi0InvCrossSectionStat2760GeVScale->Draw("p,same");
	
	
	
	//graphNLODSS14MuHalfPi02760GeV->Draw("same,c");
	//graphNLODSS14MuOnePi02760GeV->Draw("same,c");
	//graphNLODSS14MuTwoPi02760GeV->Draw("same,c");
	
	graphNLODSS14MuOne7TeV->GetXaxis()->SetRangeUser(0.5,10.5);

	graphNLODSS14MuOne7TeV->Draw("same,c");
	
	graphNLODSS14MuHalfPi02760GeV->GetXaxis()->SetRangeUser(0.5,10.5);
	graphNLODSS14MuOnePi02760GeV->GetXaxis()->SetRangeUser(0.5,10.5);
	graphNLODSS14MuTwoPi02760GeV->GetXaxis()->SetRangeUser(0.5,10.5);
	
	graphNLODSS14MuHalfPi02760GeV->Draw("same,c");
	graphNLODSS14MuOnePi02760GeV->Draw("same,c");
	graphNLODSS14MuTwoPi02760GeV->Draw("same,c");
	histoPythia8Pi02760GeVInvXSection->Draw("same,hist,c");

	

 	fitTsallisInvCrossSectionPi07TeV->Draw("same");
 	fitHistogramTsallisInvCrossSectionPi02760GeV->Draw("same,c");
	
 	

 	TLatex *labelScalingNLODSS14Pi07TeV = new TLatex(0.47,2E10,"x 1");
 	SetStyleTLatex( labelScalingNLODSS14Pi07TeV, 0.025,4,graphDalitzPi0InvCrossSectionSys7TeV->GetLineColor(),62,kFALSE);
 	labelScalingNLODSS14Pi07TeV->Draw();
 	
 	
 	TLatex *labelScalingNLODSS14Pi02760GeV = new TLatex(0.47,2E9,"x 10^{-1}");
 	SetStyleTLatex( labelScalingNLODSS14Pi02760GeV, 0.025,4,graphDalitzPi0InvCrossSectionSys2760GeVScale->GetLineColor(),62,kFALSE);
 	labelScalingNLODSS14Pi02760GeV->Draw();
	
	legendPlotLabel->Draw();
// 	
// 	legendPlotLabel->Draw();
 	
	padComparisonXSectionNLODSS14OnlySpectraPi0->Update();
 	padXSectionNLODSS14LegendOnlySpectraPi0->cd();
 	DrawGammaPadSettings( padXSectionNLODSS14LegendOnlySpectraPi0, 0., 0., 0., 0.);
 	padXSectionNLODSS14LegendOnlySpectraPi0->SetBorderMode(-1);
 	padXSectionNLODSS14LegendOnlySpectraPi0->SetBorderSize(3);
 	padXSectionNLODSS14LegendOnlySpectraPi0->Draw();
 	padXSectionNLODSS14LegendOnlySpectraPi0->cd();
// 
 	//*************** first Column **********************************************************	
 	textSpectrumALLEnergies->Draw();
 	textTsallisFitALLEnergies->Draw();
 	textNLODSS14MuHalfALLEnergies->Draw();
 	textNLODSS14MuOneALLEnergies->Draw();
 	textNLODSS14MuTwoALLEnergies->Draw();
	textPythia8ALLEnergies->Draw();
 	
// 	//*************** second Column **********************************************************
	textPi07TeVNLOALLEnergies->Draw();
 	textPi07TeVNLOsysALLEnergies->Draw();
 	boxCombinedPi07TeVALLEnergies->Draw("l");
 	markerCombinedPi07TeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
 	lineTsallisFit7TeVNLOALLEnergies->Draw("same");
 	//lineNLOPi07TeVMuHalfALLEnergies->Draw("same");
 	lineNLOPi07TeVMuOneALLEnergies->Draw("same");
 	//lineNLOPi07TeVMuTwoALLEnergies->Draw("same");
// 	
// 	//*************** third Column **********************************************************
 	textPi02760GeVNLOALLEnergies->Draw();
 	textPi02760GeVNLOsysALLEnergies->Draw();
 	boxCombinedPi02760GeVALLEnergies->Draw("l");
 	markerCombinedPi02760GeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
 	lineTsallisFitPi02760GeVNLOALLEnergies->Draw("same");
 	lineNLOPi02760GeVMuHalfALLEnergies->Draw("same");
 	lineNLOPi02760GeVMuOneALLEnergies->Draw("same");
 	lineNLOPi02760GeVMuTwoALLEnergies->Draw("same");
	linePythia8ALLEnergies->Draw("same");
 	
	canvasInvXSectionNLODSS14OnlySpectraPi0->Update();
	canvasInvXSectionNLODSS14OnlySpectraPi0->Print(Form("%s/InvariantCrossSecPP_pQCD_TsallisFit_DSS14.%s",outputDir.Data(),suffix.Data()));

	
	
	
	
	TCanvas* canvasInvXSectionBylinkinNLODSS14OnlySpectraPi0 = new TCanvas("canvasInvXSectionBylinkinNLODSS14OnlySpectraPi0","",200,10,1200,1200);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionBylinkinNLODSS14OnlySpectraPi0,  0.15, 0.02, 0.03, 0.06);
	
	
	TPad* padComparisonXSectionBylinkinNLODSS14OnlySpectraPi0 = new TPad("padComparisonXSectionBylinkinNLODSS14OnlySpectraPi0", "", 0., 0., 1., 0.77,-1, -1, -2);
	DrawGammaPadSettings( padComparisonXSectionBylinkinNLODSS14OnlySpectraPi0, 0.15, 0.02, 0.02, 0.09);
	padComparisonXSectionBylinkinNLODSS14OnlySpectraPi0->Draw();
	
	TPad* padXSectionBylinkinNLODSS14LegendOnlySpectraPi0 = new TPad("padXSectionBylinkinNLODSS14LegendOnlySpectraPi0", "", 0.15, 0.77, 0.97, 1.,-1, -1, -2);
	DrawGammaPadSettings( padXSectionBylinkinNLODSS14LegendOnlySpectraPi0, 0.15, 0.02, 0.03, 0.);
	padXSectionBylinkinNLODSS14LegendOnlySpectraPi0->Draw();
	
	padComparisonXSectionBylinkinNLODSS14OnlySpectraPi0->cd();
	padComparisonXSectionBylinkinNLODSS14OnlySpectraPi0->SetLogy();		
	padComparisonXSectionBylinkinNLODSS14OnlySpectraPi0->SetLogx();		
	
	
	
	TH2F * histo2DInvXSectionBylinkinNLODSS14OnlySpectraPi0 = new TH2F("histo2DInvXSectionBylinkinNLODSS14OnlySpectraPi0","histo2DInvXSectionBylinkinNLODSS14OnlySpectraPi0",1000,0.23,30.,1000,1e2,10e10);
	SetStyleHistoTH2ForGraphs(histo2DInvXSectionBylinkinNLODSS14OnlySpectraPi0, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
	histo2DInvXSectionBylinkinNLODSS14OnlySpectraPi0->GetXaxis()->SetRangeUser(0.4,12.0);
	histo2DInvXSectionBylinkinNLODSS14OnlySpectraPi0->DrawCopy(); 
	
	graphDalitzPi0InvCrossSectionSys7TeV->Draw("2,same");
	graphDalitzPi0InvCrossSectionStat7TeV->Draw("p,same");
	
	graphDalitzPi0InvCrossSectionSys2760GeVScale->Draw("2,same");
	graphDalitzPi0InvCrossSectionStat2760GeVScale->Draw("p,same");
	
	
	
	graphNLODSS14MuHalfPi02760GeV->Draw("same,c");
	graphNLODSS14MuOnePi02760GeV->Draw("same,c");
	graphNLODSS14MuTwoPi02760GeV->Draw("same,c");

	graphNLODSS14MuOne7TeV->Draw("same,c");
	
	graphNLODSS14MuHalfPi02760GeV->GetXaxis()->SetRangeUser(0.5,10.5);
	graphNLODSS14MuOnePi02760GeV->GetXaxis()->SetRangeUser(0.5,10.5);
	graphNLODSS14MuTwoPi02760GeV->GetXaxis()->SetRangeUser(0.5,10.5);
	
	graphNLODSS14MuHalfPi02760GeV->Draw("same,c");
	graphNLODSS14MuOnePi02760GeV->Draw("same,c");
	graphNLODSS14MuTwoPi02760GeV->Draw("same,c");

	
        fitBylinkinDalitzPi0InvCrossSection7TeV->Draw("same");
 	fitHistogramBylinkinInvCrossSection2760GeV->Draw("same,c");
	histoPythia8Pi02760GeVInvXSection->Draw("same,hist,c");
	
 	//fitBylinkinDalitzPi0InvCrossSection2760GeV->Draw("same");
 	//fitHistogramBylinkinDalitzPi0InvCrossSection2760GeV->Draw("same,c");
	
 	

 	TLatex *labelScalingBylinkinNLODSS14Pi07TeV = new TLatex(0.47,2E10,"x 1");
 	SetStyleTLatex( labelScalingBylinkinNLODSS14Pi07TeV, 0.025,4,graphDalitzPi0InvCrossSectionSys7TeV->GetLineColor(),62,kFALSE);
 	labelScalingBylinkinNLODSS14Pi07TeV->Draw();
 	
 	
 	TLatex *labelScalingBylinkinNLODSS14Pi02760GeV = new TLatex(0.47,2E9,"x 10^{-1}");
 	SetStyleTLatex( labelScalingBylinkinNLODSS14Pi02760GeV, 0.025,4,graphDalitzPi0InvCrossSectionSys2760GeVScale->GetLineColor(),62,kFALSE);
 	labelScalingBylinkinNLODSS14Pi02760GeV->Draw();
	
	legendPlotLabel->Draw();
	
	padComparisonXSectionBylinkinNLODSS14OnlySpectraPi0->Update();
 	padXSectionBylinkinNLODSS14LegendOnlySpectraPi0->cd();
 	DrawGammaPadSettings( padXSectionBylinkinNLODSS14LegendOnlySpectraPi0, 0., 0., 0., 0.);
 	padXSectionBylinkinNLODSS14LegendOnlySpectraPi0->SetBorderMode(-1);
 	padXSectionBylinkinNLODSS14LegendOnlySpectraPi0->SetBorderSize(3);
 	padXSectionBylinkinNLODSS14LegendOnlySpectraPi0->Draw();
 	padXSectionBylinkinNLODSS14LegendOnlySpectraPi0->cd();
// 
 	//*************** first Column **********************************************************	
 	textSpectrumALLEnergies->Draw();
 	textBylinkinFitALLEnergies->Draw();
 	textNLODSS14MuHalfALLEnergies->Draw();
 	textNLODSS14MuOneALLEnergies->Draw();
 	textNLODSS14MuTwoALLEnergies->Draw();
	textPythia8ALLEnergies->Draw();
 	
// 	//*************** second Column **********************************************************
	textPi07TeVNLOALLEnergies->Draw();
 	textPi07TeVNLOsysALLEnergies->Draw();
 	boxCombinedPi07TeVALLEnergies->Draw("l");
 	markerCombinedPi07TeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
 	lineTsallisFit7TeVNLOALLEnergies->Draw("same");
 	lineNLOPi07TeVMuHalfALLEnergies->Draw("same");
 	lineNLOPi07TeVMuOneALLEnergies->Draw("same");
 	lineNLOPi07TeVMuTwoALLEnergies->Draw("same");
// 	
// 	//*************** third Column **********************************************************
 	textPi02760GeVNLOALLEnergies->Draw();
 	textPi02760GeVNLOsysALLEnergies->Draw();
 	boxCombinedPi02760GeVALLEnergies->Draw("l");
 	markerCombinedPi02760GeVALLEnergies->DrawMarker(columnsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinal[2]+offsetNLOLegendPrelAndFinalMarkerY);
 	lineTsallisFitPi02760GeVNLOALLEnergies->Draw("same");
 	lineNLOPi02760GeVMuHalfALLEnergies->Draw("same");
 	lineNLOPi02760GeVMuOneALLEnergies->Draw("same");
 	lineNLOPi02760GeVMuTwoALLEnergies->Draw("same");
	linePythia8ALLEnergies->Draw("same");
 	
	canvasInvXSectionBylinkinNLODSS14OnlySpectraPi0->Update();
	canvasInvXSectionBylinkinNLODSS14OnlySpectraPi0->Print(Form("%s/InvariantCrossSecPP_pQCD_BylinkinFit_DSS14.%s",outputDir.Data(),suffix.Data()));

	
	
	
	
	
	TCanvas* canvasInvXSectionNLODSS14OnlyRatioPi0 	= new TCanvas("canvasInvXSectionNLODSS14OnlyRatioPi0","",200,10,1200,834);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionNLODSS14OnlyRatioPi0,  0.15, 0.02, 0.03, 0.06);  
	
 	TPad* padXSectionNLODSS14LegendOnlyRatioPi0 		= new TPad("padXSectionNLODSS14LegendOnlyRatioPi0", 	   "", 0.15, 0.66, 0.957,  1.,-1, -1, -2);
 	DrawGammaPadSettings( padXSectionNLODSS14LegendOnlyRatioPi0, 0.15, 0.02, 0.01, 0.0);
 	padXSectionNLODSS14LegendOnlyRatioPi0->Draw();

	TPad* padXSectionNLODSS14OnlyRatioPi0Pi07TeV 	= new TPad("padXSectionNLODSS14OnlyRatioPi0Pi07TeV",    "", 0., 0.366, 1., 0.66,-1, -1, -2);
	DrawGammaPadSettings( padXSectionNLODSS14OnlyRatioPi0Pi07TeV, 0.15, 0.02, 0.02, 0.);
	padXSectionNLODSS14OnlyRatioPi0Pi07TeV->Draw();
	
 	TPad* padXSectionNLODSS14OnlyRatioPi0Pi02760GeV	= new TPad("padXSectionNLODSS14OnlyRatioPi0Pi02760GeV", "", 0.,0.0, 1., 0.366,-1, -1, -2);
 	DrawGammaPadSettings( padXSectionNLODSS14OnlyRatioPi0Pi02760GeV, 0.15, 0.02, 0., 0.28);
 	padXSectionNLODSS14OnlyRatioPi0Pi02760GeV->Draw();
 	
 	
 	padXSectionNLODSS14OnlyRatioPi0Pi07TeV->cd();
 	padXSectionNLODSS14OnlyRatioPi0Pi07TeV->SetLogx();
 	
 	TH2F * ratio2DInvXSectionNLODSS14Pi0 = new TH2F("ratio2DInvXSectionNLODSS14Pi0","ratio2DInvXSectionNLODSS14Pi0",1000,0.23,30.,1000,0.4,3.55);
 	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionNLODSS14Pi0, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 36,36,36,36,2.7,0.8, 512, 505);
 	//Jus for differen size of Tpads
 	ratio2DInvXSectionNLODSS14Pi0->GetXaxis()->SetLabelFont(43);
 	ratio2DInvXSectionNLODSS14Pi0->GetYaxis()->SetLabelFont(43); 
 	ratio2DInvXSectionNLODSS14Pi0->GetXaxis()->SetTitleFont(63);
 	ratio2DInvXSectionNLODSS14Pi0->GetYaxis()->SetTitleFont(63);
 	ratio2DInvXSectionNLODSS14Pi0->DrawCopy(); 
 	
	
	graphRatioNLODSS14DalitzPi07TeVMuOne->Draw("same,c");
 	boxErrorSigmaPi07TeVRatio->Draw();
 	graphRatioInvCrossSectionPi07TeVTsallisFit->Draw("p,E2same");
 	
 	labelRatioNLOPi07TeV->Draw();
 	
 	DrawGammaLines(0., 30.,1., 1.,0.1);
 	
 	padXSectionNLODSS14OnlyRatioPi0Pi02760GeV->cd();
 	padXSectionNLODSS14OnlyRatioPi0Pi02760GeV->SetLogx();
 	
 	TH2F *  ratio2DInvXSectionNLODSS14Pi02760GeV;
	ratio2DInvXSectionNLODSS14Pi02760GeV = new TH2F(" ratio2DInvXSectionNLODSS14Pi02760GeV"," ratio2DInvXSectionNLODSS14Pi02760GeV",1000,0.23,30.,1000,0.4,3.55);
 	SetStyleHistoTH2ForGraphs( ratio2DInvXSectionNLODSS14Pi02760GeV, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}",36,36,36,36, 2.7,0.80, 512, 505);
 	ratio2DInvXSectionNLODSS14Pi02760GeV->GetXaxis()->SetLabelFont(43);
 	ratio2DInvXSectionNLODSS14Pi02760GeV->GetYaxis()->SetLabelFont(43); 
 	ratio2DInvXSectionNLODSS14Pi02760GeV->GetXaxis()->SetTitleFont(63);
 	ratio2DInvXSectionNLODSS14Pi02760GeV->GetYaxis()->SetTitleFont(63);
 	ratio2DInvXSectionNLODSS14Pi02760GeV->DrawCopy(); 
	
	DrawGammaSetMarker(histoRatioPythia8ToFitTsallis, 24, 1.5, kRed+2 , kRed+2);  
        histoRatioPythia8ToFitTsallis->SetLineWidth(widthCommonFit);
        histoRatioPythia8ToFitTsallis->Draw("same,hist,c");
 	
 	graphRatioNLODSS14DalitzPi02760GeVMuHalf->Draw("same,c");
 	graphRatioNLODSS14DalitzPi02760GeVMuOne->Draw("same,c");
 	graphRatioNLODSS14DalitzPi02760GeVMuTwo->Draw("same,c");
 	
 	boxErrorSigmaPi02760GeVRatio->Draw();
 	graphRatioInvCrossSectionPi02760GeVTsallisFit->Draw("p,E2same");
 	labelRatioNLOPi02760GeV->Draw();
	
 	
 	
 	
 	DrawGammaLines(0., 30.,1., 1.,0.1);
 	
 			
 	padXSectionNLODSS14LegendOnlyRatioPi0->cd();
 	DrawGammaPadSettings( padXSectionNLODSS14LegendOnlyRatioPi0, 0., 0., 0.0, 0.0);
 	padXSectionNLODSS14LegendOnlyRatioPi0->SetBorderMode(-1);
 	padXSectionNLODSS14LegendOnlyRatioPi0->SetBorderSize(3);
 	padXSectionNLODSS14LegendOnlyRatioPi0->Draw();
 	padXSectionNLODSS14LegendOnlyRatioPi0->cd();
 	
 	//*************** first Column **********************************************************
 	textFitTsallisSpectrumOnlyRatioPi0->Draw();
 	textNLODSS14MuHalfOnlyRatioPi0->Draw();
 	textNLODSS14MuOneOnlyRatioPi0->Draw();
 	textNLODSS14MuTwoOnlyRatioPi0->Draw();
	textPythia8RatioPi0->Draw();
 	

 	
 	//*************** second Column **********************************************************
 	textPi07TeVNLOOnlyRatioPi0->Draw();
 	textPi07TeVNLOsysOnlyRatioPi0->Draw();
 	boxFitTsallisPi07TeVOnlyRatioPi0->Draw("l");
 	markerCombinedPi07TeVOnlyRatioPi0->DrawMarker(columnsNLOLegendPrelAndFinalOnlyRatio[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY);
 	
 	
 	//lineNLOPi07TeVMuHalfOnlyRatioPi0->Draw("same");
 	lineNLOPi07TeVMuOneOnlyRatioPi0->Draw("same");
 	//lineNLOPi07TeVMuTwoOnlyRatioPi0->Draw("same");
 	
	
 	//*************** third Column **********************************************************
 	
	textPi02760GeVNLOOnlyRatioPi0->Draw();
 	textPi02760GeVNLOsysOnlyRatioPi0->Draw();
 	
 	boxTSallisFitPi02760GeVOnlyRatioPi0->Draw("l");
 	markerCombinedPi02760GeVOnlyRatioPi0->DrawMarker(columnsNLOLegendPrelAndFinalOnlyRatio[2]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY);
 	
 	lineNLOPi02760GeVMuHalfOnlyRatioPi0->Draw("same");
 	lineNLOPi02760GeVMuOneOnlyRatioPi0->Draw("same");
 	lineNLOPi02760GeVMuTwoOnlyRatioPi0->Draw("same"); 
	linePythi82760GeVOnlyRationPi0->Draw("same");
 
	//************************************************* End Legend ***************************************************
 		
	canvasInvXSectionNLODSS14OnlyRatioPi0->Update();
	canvasInvXSectionNLODSS14OnlyRatioPi0->Print(Form("%s/Pi0_InvXSectionNLO_Tsallis_OnlyRatio_DSS14.%s",outputDir.Data(),suffix.Data()));
	
	
		
	
	
	TCanvas* canvasInvXSectionBylinkinNLODSS14OnlyRatioPi0 	= new TCanvas("canvasInvXSectionBylinkinNLODSS14OnlyRatioPi0","",200,10,1200,834);  // gives the page size
	DrawGammaCanvasSettings( canvasInvXSectionBylinkinNLODSS14OnlyRatioPi0,  0.15, 0.02, 0.03, 0.06);  
	
 	TPad* padXSectionBylinkinNLODSS14LegendOnlyRatioPi0 		= new TPad("padXSectionBylinkinNLODSS14LegendOnlyRatioPi0", 	   "", 0.15, 0.66, 0.957,  1.,-1, -1, -2);
 	DrawGammaPadSettings( padXSectionBylinkinNLODSS14LegendOnlyRatioPi0, 0.15, 0.02, 0.01, 0.0);
 	padXSectionBylinkinNLODSS14LegendOnlyRatioPi0->Draw();

	TPad* padXSectionBylinkinNLODSS14OnlyRatioPi0Pi07TeV 	= new TPad("padXSectionBylinkinNLODSS14OnlyRatioPi0Pi07TeV",    "", 0., 0.366, 1., 0.66,-1, -1, -2);
	DrawGammaPadSettings( padXSectionBylinkinNLODSS14OnlyRatioPi0Pi07TeV, 0.15, 0.02, 0.02, 0.);
	padXSectionBylinkinNLODSS14OnlyRatioPi0Pi07TeV->Draw();
	
 	TPad* padXSectionBylinkinNLODSS14OnlyRatioPi0Pi02760GeV	= new TPad("padXSectionBylinkinNLODSS14OnlyRatioPi0Pi02760GeV", "", 0.,0.0, 1., 0.366,-1, -1, -2);
 	DrawGammaPadSettings( padXSectionBylinkinNLODSS14OnlyRatioPi0Pi02760GeV, 0.15, 0.02, 0., 0.28);
 	padXSectionBylinkinNLODSS14OnlyRatioPi0Pi02760GeV->Draw();
 	
 	
 	padXSectionBylinkinNLODSS14OnlyRatioPi0Pi07TeV->cd();
 	padXSectionBylinkinNLODSS14OnlyRatioPi0Pi07TeV->SetLogx();
 	
 	TH2F * ratio2DInvXSectionBylinkinNLODSS14Pi0 = new TH2F("ratio2DInvXSectionBylinkinNLODSS14Pi0","ratio2DInvXSectionBylinkinNLODSS14Pi0",1000,0.23,30.,1000,0.4,3.55);
 	SetStyleHistoTH2ForGraphs(ratio2DInvXSectionBylinkinNLODSS14Pi0, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}", 36,36,36,36,2.7,0.8, 512, 505);
 	//Jus for differen size of Tpads
 	ratio2DInvXSectionBylinkinNLODSS14Pi0->GetXaxis()->SetLabelFont(43);
 	ratio2DInvXSectionBylinkinNLODSS14Pi0->GetYaxis()->SetLabelFont(43); 
 	ratio2DInvXSectionBylinkinNLODSS14Pi0->GetXaxis()->SetTitleFont(63);
 	ratio2DInvXSectionBylinkinNLODSS14Pi0->GetYaxis()->SetTitleFont(63);
 	ratio2DInvXSectionBylinkinNLODSS14Pi0->DrawCopy(); 
 	
	
	graphRatioNLODSS14DalitzPi07TeVMuOne->Draw("same,c");
 	boxErrorSigmaPi07TeVRatio->Draw();
 	graphRatioInvCrossSectionPi07TeVBylinkinFit->Draw("p,E2same");
 	
 	labelRatioNLOPi07TeV->Draw();
 	
 	DrawGammaLines(0., 30.,1., 1.,0.1);
 	
 	padXSectionBylinkinNLODSS14OnlyRatioPi0Pi02760GeV->cd();
 	padXSectionBylinkinNLODSS14OnlyRatioPi0Pi02760GeV->SetLogx();
 	
 	TH2F *  ratio2DInvXSectionBylinkinNLODSS14Pi02760GeV;
	ratio2DInvXSectionBylinkinNLODSS14Pi02760GeV = new TH2F(" ratio2DInvXSectionBylinkinNLODSS14Pi02760GeV"," ratio2DInvXSectionBylinkinNLODSS14Pi02760GeV",1000,0.23,30.,1000,0.4,3.55);
 	SetStyleHistoTH2ForGraphs( ratio2DInvXSectionBylinkinNLODSS14Pi02760GeV, "#it{p}_{T} (GeV/#it{c})","#frac{NLO}{fit}",36,36,36,36, 2.7,0.80, 512, 505);
 	ratio2DInvXSectionBylinkinNLODSS14Pi02760GeV->GetXaxis()->SetLabelFont(43);
 	ratio2DInvXSectionBylinkinNLODSS14Pi02760GeV->GetYaxis()->SetLabelFont(43); 
 	ratio2DInvXSectionBylinkinNLODSS14Pi02760GeV->GetXaxis()->SetTitleFont(63);
 	ratio2DInvXSectionBylinkinNLODSS14Pi02760GeV->GetYaxis()->SetTitleFont(63);
 	ratio2DInvXSectionBylinkinNLODSS14Pi02760GeV->DrawCopy(); 
 	
	DrawGammaSetMarker(histoRatioPythia8ToFitBylinkin, 24, 1.5, kRed+2 , kRed+2);  
        histoRatioPythia8ToFitBylinkin->SetLineWidth(widthCommonFit);
        histoRatioPythia8ToFitBylinkin->Draw("same,hist,c");
 	
	
 	graphRatioBylinkinFitNLODSS14DalitzPi02760GeVMuHalf->Draw("same,c");
 	graphRatioBylinkinFitNLODSS14DalitzPi02760GeVMuOne->Draw("same,c");
 	graphRatioBylinkinFitNLODSS14DalitzPi02760GeVMuTwo->Draw("same,c");
 	
 	boxErrorSigmaPi02760GeVRatio->Draw();
 	graphRatioInvCrossSectionPi02760GeVTsallisFit->Draw("p,E2same");
 	labelRatioNLOPi02760GeV->Draw();
	
 	
 	
 	
 	DrawGammaLines(0., 30.,1., 1.,0.1);
 	
 			
 	padXSectionBylinkinNLODSS14LegendOnlyRatioPi0->cd();
 	DrawGammaPadSettings( padXSectionBylinkinNLODSS14LegendOnlyRatioPi0, 0., 0., 0.0, 0.0);
 	padXSectionBylinkinNLODSS14LegendOnlyRatioPi0->SetBorderMode(-1);
 	padXSectionBylinkinNLODSS14LegendOnlyRatioPi0->SetBorderSize(3);
 	padXSectionBylinkinNLODSS14LegendOnlyRatioPi0->Draw();
 	padXSectionBylinkinNLODSS14LegendOnlyRatioPi0->cd();
 	
 	//*************** first Column **********************************************************
 	textFitBylinkinSpectrumOnlyRatioPi0->Draw();
 	textNLODSS14MuHalfOnlyRatioPi0->Draw();
 	textNLODSS14MuOneOnlyRatioPi0->Draw();
 	textNLODSS14MuTwoOnlyRatioPi0->Draw();
	textPythia8RatioPi0->Draw();
 	

 	
 	//*************** second Column **********************************************************
 	textPi07TeVNLOOnlyRatioPi0->Draw();
 	textPi07TeVNLOsysOnlyRatioPi0->Draw();
 	boxFitTsallisPi07TeVOnlyRatioPi0->Draw("l");
 	markerCombinedPi07TeVOnlyRatioPi0->DrawMarker(columnsNLOLegendPrelAndFinalOnlyRatio[1]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY);
 	
 	
 	//lineNLOPi07TeVMuHalfOnlyRatioPi0->Draw("same");
 	lineNLOPi07TeVMuOneOnlyRatioPi0->Draw("same");
 	//lineNLOPi07TeVMuTwoOnlyRatioPi0->Draw("same");
 	
	
 	//*************** third Column **********************************************************
 	
	textPi02760GeVNLOOnlyRatioPi0->Draw();
 	textPi02760GeVNLOsysOnlyRatioPi0->Draw();
 	
 	boxTSallisFitPi02760GeVOnlyRatioPi0->Draw("l");
 	markerCombinedPi02760GeVOnlyRatioPi0->DrawMarker(columnsNLOLegendPrelAndFinalOnlyRatio[2]+offsetNLOLegendPrelAndFinalMarkerX,rowsNLOLegendPrelAndFinalOnlyRatioPi0[2]+offsetNLOLegendPrelAndFinalMarkerY);
 	
 	lineNLOPi02760GeVMuHalfOnlyRatioPi0->Draw("same");
 	lineNLOPi02760GeVMuOneOnlyRatioPi0->Draw("same");
 	lineNLOPi02760GeVMuTwoOnlyRatioPi0->Draw("same"); 
	linePythi82760GeVOnlyRationPi0->Draw("same");
 
	//************************************************* End Legend ***************************************************
 		
	canvasInvXSectionBylinkinNLODSS14OnlyRatioPi0->Update();
	canvasInvXSectionBylinkinNLODSS14OnlyRatioPi0->Print(Form("%s/Pi0_InvXSectionNLO_FitBylinkin_OnlyRatio_DSS14.%s",outputDir.Data(),suffix.Data()));
	
	
	
	
}
	
