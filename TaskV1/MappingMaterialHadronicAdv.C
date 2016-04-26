/*********************************************************************************************************
 ******  provided by Gamma Conversion Group, PWG4:                                                 ********
 ******                 Kathrin Koch, kkoch@physi.uni-heidelberg.de                                                                ********
 ******         Friederike Bock, fbock@physi.uni-heidelberg.de                                                     ********
 **********************************************************************************************************
 ***     This Macro can be used for the material analysis of the ALICE ITS and TPC up to an R of 180 cm   ***
 ***     To use this macro the GammaConversionTaks needs to be run first.                                                  ***
 ***     The macro can be run with 2 input files, the second one should be a Montecarlo file, if it is a  ***
 ***     data-file, modifications are necessary, taking out the histograms, which are only produced if    *** 
 ***     the ConversionTask is run on MC. They usually have MCtruth somewhere in its name.                ***
 **********************************************************************************************************/

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
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TMarker.h"
#include "TGraphAsymmErrors.h" 
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"

TString textGenerator;
TString collisionSystem;
TString textPeriod;
TString textDate;
using namespace std;

Double_t CalculateErrors(	fstream 	&fileDatOutput, 
						Int_t 	Bin, 
						Double_t 	entriesBinData, 
						Double_t 	entriesSumData, 
						Double_t 	entriesEventsData, 
						Double_t 	entriesBinMC, 
						Double_t 	entriesSumMC, 
						Double_t 	entriesEventsMC ) 
{
        if (Bin == 1) fileDatOutput << endl;
        Double_t entriesBinDataWithoutNorm 		= entriesBinData * entriesEventsData;
        Double_t entriesBinMCWithoutNorm 		= entriesBinMC * entriesEventsMC;
        Double_t statErrData 					= TMath::Sqrt(entriesBinDataWithoutNorm) /entriesEventsData;
        Double_t statErrMC	 				= TMath::Sqrt(entriesBinMCWithoutNorm) /entriesEventsMC;
        Double_t statErrComple 				= TMath::Sqrt(statErrMC*statErrMC + statErrData*statErrData);
	   
        fileDatOutput <<"\t" << Bin << "\t" << entriesBinData <<"\t"<< statErrData <<"\t" << entriesBinData/entriesSumData * 100 <<"\t" << entriesBinMC << "\t" << statErrMC << "\t" << entriesBinMC/entriesSumMC * 100 << "\t" << entriesBinData- entriesBinMC << "\t" << statErrComple << "\t" << (entriesBinData- entriesBinMC)/entriesSumMC *100  << "\t" << statErrComple/entriesSumMC *100 << endl;
        
	   return (entriesBinData -entriesBinMC);
}

void PlotStandard2D( TH2* histo2D, TString nameOutput, TString title, TString xTitle, TString yTitle, Bool_t kRangeY, Double_t startY, Double_t endY, Bool_t kRangeX, Double_t startX, Double_t endX, Bool_t kRangeZ, Double_t startZ, Double_t endZ, Int_t logX, Int_t logZ, Float_t* floatLogo, Int_t canvasSizeX = 500, Int_t canvasSizeY = 500, Double_t rightMargin = 0.12, Double_t leftMargin = 0.12, Double_t bottomMargin = 0.1, Double_t topMargin = 0.04, TString drawLogo ="Perf",Float_t titleOffsetX=1.4, Float_t titleOffsetY=1.2, Float_t labelOffsetZ=0.005){
	TCanvas * canvasStandard = new TCanvas("canvasStandard","",10,10,canvasSizeX,canvasSizeY);  // gives the page size		
	canvasStandard->SetLogx(logX);
	canvasStandard->SetLogz(logZ);
	canvasStandard->SetRightMargin(rightMargin); 		
	canvasStandard->SetLeftMargin(leftMargin); 		
	canvasStandard->SetBottomMargin(bottomMargin); 		
	canvasStandard->SetTopMargin(topMargin); 		
	canvasStandard->cd();
	if (kRangeZ){
		histo2D->GetZaxis()->SetRangeUser(startZ,endZ);
	}
	histo2D->GetZaxis()->SetLabelOffset(labelOffsetZ);
	DrawAutoGammaHisto2D(	histo2D,
									title.Data(), xTitle.Data(), yTitle.Data(),"",kRangeY, startY, endY, kRangeX, startX, endX,titleOffsetX, titleOffsetY);
	if (drawLogo.CompareTo("Data")==0){
		DrawLabelsEvents(floatLogo[0],floatLogo[1],floatLogo[2], 0.00, collisionSystem, "Data", textPeriod);
	} else {
		DrawLabelsEvents(floatLogo[0],floatLogo[1],floatLogo[2], 0.00, collisionSystem, textGenerator, textPeriod);
	}
// 	if (drawLogo.CompareTo("Perf")==0) DrawAliceLogoPerformance(floatLogo[0],floatLogo[1],floatLogo[2],floatLogo[3],0.00, textDate,collisionSystem, "","",canvasSizeX,canvasSizeY);	
// 	if (drawLogo.CompareTo("Wip")==0)DrawAliceLogo(floatLogo[0],floatLogo[1],floatLogo[2],floatLogo[3],canvasSizeX,canvasSizeY);	
// 	
	canvasStandard->Update();
	canvasStandard->SaveAs(nameOutput.Data());
	delete canvasStandard;
}


void  MappingMaterialHadronicAdv(const char *data = "myOutput", TString MCfile = "", TString cutSel = "", TString suffix = "gif", TString optEnergy="", TString optMCGenerator="", TString optPeriod=""){        
        
	gROOT->Reset();
	gROOT->SetStyle("Plain");
 
	
	if(optEnergy.CompareTo("7TeV") == 0){
			collisionSystem = 			"pp, #sqrt{#it{s}} = 7 TeV";          
	} else if( optEnergy.CompareTo("2.76TeV") == 0) {
		collisionSystem = 			"pp, #sqrt{#it{s}} = 2.76 TeV";
	} else if( optEnergy.CompareTo("900GeV") == 0) {
			collisionSystem = 			"pp, #sqrt{#it{s}} = 0.9 TeV";
	} else if( optEnergy.CompareTo("HI") == 0) {
			collisionSystem =			"PbPb, #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";
	} else {
			cout << "No correct collision system specification, has been given" << endl;
			return;         
	}
   TString suffixData = "_Data";
	TString suffixMC = "_MC";
   Bool_t kPureData = kFALSE;
	if(optMCGenerator.CompareTo("") ==0){
		textGenerator = "";
	} else if(optMCGenerator.CompareTo("Pass2/Pass4") ==0){
		textGenerator = "";
		kPureData = kTRUE;
		suffixData = "_DataPass2";
		suffixMC = "_DataPass4";
	} else {
		textGenerator = optMCGenerator;
	}
	TString outputDirectory;
	TString outputDirectory1;
	TString outputDirectory2;
	if(optPeriod.CompareTo("") ==0 || optPeriod.CompareTo("All") ==0){
		textPeriod = "";
		if (optMCGenerator.CompareTo("")!=0){
			outputDirectory =	 		Form("%s/%s/%s/%s/MappingMaterialHadronic",optMCGenerator.Data(),cutSel.Data(),optEnergy.Data(),suffix.Data());
		} else {
			outputDirectory =	 		Form("%s/%s/%s/MappingMaterialHadronic",cutSel.Data(),optEnergy.Data(),suffix.Data());
		}
		gSystem->Exec("mkdir -p "+outputDirectory);
	} else {
		textPeriod = optPeriod;
		if (optMCGenerator.CompareTo("")!=0){
			outputDirectory =	 		Form("%s/%s/%s/%s/%s/MappingMaterialHadronic",optMCGenerator.Data(),cutSel.Data(),optEnergy.Data(),suffix.Data(),optPeriod.Data());
		} else {
			outputDirectory =	 		Form("%s/%s/%s/%s/MappingMaterialHadronic",cutSel.Data(),optEnergy.Data(),suffix.Data(),optPeriod.Data());
		}
		gSystem->Exec("mkdir -p "+outputDirectory);
	}
	
	Int_t iRStart;
	TString rCutNumber =cutSel(15,1);
	if (rCutNumber.CompareTo("0") == 0 ){
		iRStart   	= 0;
	} else {
		iRStart 	= 1;
	}       
	StyleSettingsThesis();
	SetPlotStyle();

	// Which file do you want to analyse
	TString 	nameDataFile = 		(Form("%s",data));
		
	//rebinning
	Int_t	rebin = 				5;
	
		
	//How big should the right margin in 2D plots be? 
	Float_t 	rightMargin = 			0.12;
	Float_t 	leftMargin = 			0.11;
	
	//How many slices do you have?
	const int 	nBinsR = 				13;  
	const int 	nBinsZ = 				23;
	const int 	nBinsZHotZone = 		6;
	
	// How many raws and colums you want to have in your output?
	Int_t 		columnR = 			2;
	Int_t 		rowR = 				7;
	Int_t 		columnZ = 			3;
	Int_t 		rowZ = 				6;
	
	// Array for ZinR-Ranges
	Float_t	arrayRangeZInR[nBinsR] = {300.,300.,300.,300.,300.,300.,300.,300.,300.,300.,300.,300.,300.}; 
	
	// Array of Rbins
	Float_t 	arrayRBins[14] = 		{0.,3.5,5.75,9.5,13.,21.,27.5,35.,42.,55.,72.,79.5,90.,180.};
	TString 	arrayNamesRBins[13]=	{"Beam Pipe", 
									"SPD 1", 
									"SPD 2",
									"Thermal shield/Support between SPD/SDD", 
									"SDD 1", 
									"SDD 2", 
									"Thermal shield/Support between SDD/SSD", 
									"SSD 1", 
									"SSD 2", 
									"Air + TPC in. cont. vessel + C0_{2} (+ITS serv. 0.9 < |#eta| < 1.4)", 
									"CO_{2} + TPC in. field cage vessel (+ITS serv. 0.9 < |#eta| < 1.4)", 
									"TPC rods + Ne: CO_{2}: N_{2}", 
									"Ne: CO_{2}: N_{2}"};
	
	// Array of Zbins
	Float_t 	arrayZBins[24] = 		{-150.,-105.,-90.,-70.,-60.,-50.,-40.,-30.,-20.,-10., 						
												0.,10.,20.,30.,40.,50.,60.,70.,90.,105.,140.,180.,220.,300.};
	TString 	arrayNamesZBins[23]=	{
									"-150 < Z < -105", 
									"-105 < Z < -90", 
									"-90 < Z < -70", 
									"-70 < Z < -60", 
									"-60 < Z < -50", 
									"-50 < Z < -40", 
									"-40 < Z < -30", 
									"-30 < Z < -20", 
									"-20 < Z < -10", 
									"-10 < Z < 0", 
									"0 < Z < 10", 
									"10 < Z < 20", 
									"20 < Z < 30", 
									"30 < Z < 40", 
									"40 < Z < 50", 
									"50 < Z < 60", 
									"60 < Z < 70", 
									"70 < Z < 90", 
									"90 < Z < 105", 
									"105 < Z < 140", 
									"140 < Z < 180", 
									"180 < Z < 220", 
									"220 < Z < 300"};
	Float_t 	arrayZBinsHotZone[7] = 		{45.,50.,55.,60.,65.,70.,75.};
	TString 	arrayNamesZBinsHotZone[6]=	
									{"45 < Z < 50", 
									"50 < Z < 55", 
									"55 < Z < 60", 
									"60 < Z < 65", 
									"65 < Z < 70", 
									"70 < Z < 75"};
	
									
	// Array defintion for printing Logo in right upper corner
	Float_t 	floatLocationRightUp[4]= {0.8,0.605,0.15, 0.02};
	Float_t floatLocationRightUp2D[4] = {0.6,0.95,0.035, 0.02};
	Float_t 	floatLocationRightUp2[4]= {0.8,0.76,0.15, 0.02};
	Float_t 	floatLocationRightDown[4]=	{0.7,0.23,0.15, 0.02};
// 	Float_t 	floatLocationRightUp2D[4]=	{0.68,0.775,0.11, 0.02};
	Float_t 	floatLocationRightUp2DXY[4]=	{0.69,0.845,0.115, 0.025};
	Float_t 	floatLocationLeftLow2DZR[4]=	{0.1,0.25,0.04, 0.025};
	// Array defintion for printing Logo in left upper corner
	Float_t 	floatLocationLeftUp[4]=		{0.18,0.76, 0.15, 0.02};
	// Array defintion for printing text in right upper corner
	Float_t 	floatLocationRightUpText[4]=	{0.16,0.92,0.15, 0.04};
	Float_t 	floatLocationRightUpText2[4]=	{0.7,0.8,0.15, 0.04};
	// Get the histos
	TFile 		fileData(nameDataFile);  
	
	// choice of dateset
	TString nameGammaDirectory;
	TString nameGammaList;
	if(cutSel.CompareTo("") != 0){
			nameGammaDirectory = 		Form("SecHadInt_%s",  cutSel.Data());
			cout << nameGammaDirectory.Data() << endl;
	} else {
		cout << "no cutselection set"<< endl;
		return;
	}
	// labeling
	TString 		textStandardYAxis = 	"#gamma/ event scaled by multiplicity";
	TString 		textStandardYAxis2 = 	"#frac{dN_{sec. had}}{N_{sec. had}}";
	
	TDatime 	now;
	int 		iDate = now.GetDate();
	int 		iYear=iDate/10000;
	int 		iMonth=(iDate%10000)/100;
	int 		iDay=iDate%100;
	char* 	cMonth[12]={"Jan","Feb","Mar","Apr","May","Jun",
		    "Jul","Aug","Sep","Oct","Nov","Dec"};
	TString 	textDayth;
	if (iDay%10 == 1){
		textDayth = "st";
	} else if (iDay%10 == 2){
		textDayth = "nd";
	} else if (iDay%10 == 3){
		textDayth = "rd";
	} else {
		textDayth = "th";
	}
	textDate = Form("%i^{%s} %s %i",iDay, textDayth.Data(),cMonth[iMonth-1], iYear);

	TString 		textYAxisRHisto = "#frac{dN_{sec. had}}{N_{ch}dR}";
	TString 		textYAxisZHisto = "#frac{dN_{sec. had}}{N_{ch}dZ}";
	TString 		textYAxisEtaHisto = "#frac{dN_{sec. had}}{N_{ch}d#eta}";
	TString 		textYAxisSlicesHisto = "#frac{dN_{sec. had}}{N_{ch}dZ d#Phi dR}";

	TString etaCutNumber = cutSel(1,1);
	Bool_t boolRatioLog = kTRUE;
	Double_t rmin = 2.0;
	Double_t minYValueRatio = 0.2;
	Double_t maxYValueRatio = 10.;
	Double_t minYValueRestRRatio = 0.05;
	Double_t maxYValueRestRRatio = 10.;
	Int_t rebinIntPlots = 2;
	Double_t minYValueZPlot = 1e-4;
	Double_t minYValueRPlot = 1e-4;
	Double_t minYValueRPlotSPD = 1e-3;
	Double_t minYValueRPlotMinPt = 1e-7;
	Double_t minYValueZPlotMinPt = 1e-7;
	cout << "here" << endl;
	
	//------------------------------- Reading FILES ----------------------------------------------------------------------
	TDirectory 	*directorySecHadIntData = 	new TDirectory(); // definition of first folder / list  

	directorySecHadIntData = 			(TDirectory*)fileData.Get(nameGammaDirectory); 
	
	TH1F* histoEventQualityData=					 	(TH1F*)directorySecHadIntData->Get("ESD_EventQuality"); 
	TH2F *histoHadIntVsZRData=			(TH2F*)directorySecHadIntData->Get(Form("ESD_HadIntMap_ZR"));    
	TH2F *histoHadIntVsXYData=			(TH2F*)directorySecHadIntData->Get(Form("ESD_HadIntMap_XY"));    
	TH2F *histoHadIntVsXYInnerBPData=			(TH2F*)directorySecHadIntData->Get(Form("ESD_HadIntMapInnerBeampipe_XY"));    
	TH1D *histoHadIntVsZData=			(TH1D*)histoHadIntVsZRData->ProjectionX("histoHadIntVsZData");  
	ConvGammaRebinWithBinCorrection(histoHadIntVsZData,rebinIntPlots);
	TH1D *histoHadIntVsRData=			(TH1D*)histoHadIntVsZRData->ProjectionY("histoHadIntVsRData");        
 	ConvGammaRebinWithBinCorrection(histoHadIntVsRData,rebinIntPlots);
	TH1D* 	histoHadIntVsRCentralData=			(TH1D*)histoHadIntVsZRData->ProjectionY("histoHadIntVsRCentralData", histoHadIntVsZRData->GetXaxis()->FindBin(-10.),histoHadIntVsZRData->GetXaxis()->FindBin(10.));

	Double_t integralBPData = histoHadIntVsRData->Integral(histoHadIntVsRData->GetXaxis()->FindBin(2.5),histoHadIntVsRData->GetXaxis()->FindBin(3.5));
	Double_t integralGasData = histoHadIntVsRData->Integral(histoHadIntVsRData->GetXaxis()->FindBin(90),histoHadIntVsRData->GetXaxis()->FindBin(100));
	
	cout <<"Number of events in BP:" << integralBPData << endl;
	cout <<"Number of events in Gas:" << integralBPData << endl;
	
	TH2F *histoHadIntVsZRData2 =	(TH2F*) histoHadIntVsZRData->Clone("histoHadIntVsZRData2");
	TH1F* histoHadIntQualChi2NDOFData=	(TH1F*)directorySecHadIntData->Get("ESD_HadIntQual_Chi2PerDOF"); 
	TH1F* histoHadIntQualErrXData=	(TH1F*)directorySecHadIntData->Get("ESD_HadIntQual_ErrX"); 
	TH1F* histoHadIntQualErrYData=	(TH1F*)directorySecHadIntData->Get("ESD_HadIntQual_ErrY"); 
	TH1F* histoHadIntQualErrZData=	(TH1F*)directorySecHadIntData->Get("ESD_HadIntQual_ErrZ"); 
	TH1F* histoHadIntQualErr2DData=	(TH1F*)directorySecHadIntData->Get("ESD_HadIntQual_Err2D"); 
	TH1F* histoHadIntQualErr3DData=	(TH1F*)directorySecHadIntData->Get("ESD_HadIntQual_Err3D"); 
	TH1F* histoHadIntQualNTracksData=	(TH1F*)directorySecHadIntData->Get("ESD_HadIntQual_nTracksSecVtx"); 
	TH1F* 	histoHadIntQualChi2NDOFAftCutsData=	(TH1F*)directorySecHadIntData->Get(Form("ESD_HadIntQualAftCuts_Chi2PerDOF"));   
	TH1F* 	histoHadIntQualErrXAftCutsData=	(TH1F*)directorySecHadIntData->Get(Form("ESD_HadIntQualAftCuts_ErrX"));  
	TH1F* 	histoHadIntQualErrYAftCutsData=	(TH1F*)directorySecHadIntData->Get(Form("ESD_HadIntQualAftCuts_ErrY"));  
	TH1F* 	histoHadIntQualErrZAftCutsData=	(TH1F*)directorySecHadIntData->Get(Form("ESD_HadIntQualAftCuts_ErrZ"));   
	TH1F* 	histoHadIntQualErr2DAftCutsData=	(TH1F*)directorySecHadIntData->Get(Form("ESD_HadIntQualAftCuts_Err2D"));   
	TH1F* 	histoHadIntQualErr3DAftCutsData=	(TH1F*)directorySecHadIntData->Get(Form("ESD_HadIntQualAftCuts_Err3D"));   
	TH1F* 	histoHadIntMapEtaData=	(TH1F*)directorySecHadIntData->Get(Form("ESD_HadIntMap_Eta"));   
	
	Double_t binWidthR = histoHadIntVsRData->GetXaxis()->GetBinWidth(3);
	Double_t binWidthZ = histoHadIntVsZData->GetXaxis()->GetBinWidth(3);
	
	for (Int_t iX = 1 ; iX < histoHadIntVsXYInnerBPData->GetNbinsX(); iX ++){
		for (Int_t iY = 1; iY < histoHadIntVsXYInnerBPData->GetNbinsY(); iY ++) {
			Double_t rcalc = TMath::Sqrt(histoHadIntVsXYInnerBPData->GetXaxis()->GetBinCenter(iX)*histoHadIntVsXYInnerBPData->GetXaxis()->GetBinCenter(iX) + histoHadIntVsXYInnerBPData->GetYaxis()->GetBinCenter(iY)*histoHadIntVsXYInnerBPData->GetYaxis()->GetBinCenter(iY));
			if (rcalc < rmin){
				histoHadIntVsXYInnerBPData->SetBinContent(iX,iY,NULL);
			}
		}
	}

	Float_t	numberGoodEventsData = 					histoEventQualityData->GetEntries()-histoEventQualityData->GetBinContent(6)-histoEventQualityData->GetBinContent(4);

	Float_t numberGoodTriggerData = 				histoEventQualityData->GetEntries();
	Float_t numberReconstGammaData =				histoHadIntVsRData->GetEntries();
	cout<< data << "    Number of events::   " << numberGoodEventsData << "    Number of triggers::   " << numberGoodTriggerData << "    Number reconstructed vertices::    "<< numberReconstGammaData <<endl;
	Double_t meanMultiplicityData = 1;
	Float_t normFactorReconstData=				1./integralBPData;//*1./meanMultiplicityData;
	Double_t integralNTracksData = 				1./histoHadIntQualNTracksData->Integral();
	Double_t integralChi2Data = 					1./histoHadIntQualChi2NDOFData->Integral();
	Double_t integralErr2DData = 					1./histoHadIntQualErr2DData->Integral();
	Double_t integralErr3DData = 					1./histoHadIntQualErr3DData->Integral();
	Double_t integralErrXData = 					1./histoHadIntQualErrXData->Integral();
	Double_t integralErrYData = 					1./histoHadIntQualErrYData->Integral();
	Double_t integralErrZData = 					1./histoHadIntQualErrZData->Integral();
	
	//Scaling reconstr.
	GammaScalingHistogramm(histoHadIntVsRData,normFactorReconstData);
	GammaScalingHistogramm(histoHadIntMapEtaData,normFactorReconstData);
	GammaScalingHistogramm(histoHadIntVsRCentralData,normFactorReconstData);
	GammaScalingHistogramm(histoHadIntVsZData,normFactorReconstData);
	GammaScalingHistogramm(histoHadIntVsZRData2,1./numberGoodEventsData);
 	GammaScalingHistogramm(histoHadIntQualNTracksData,integralNTracksData);
// 	GammaScalingHistogramm(histoHadIntVsXYData,normFactorReconstData);	
	GammaScalingHistogramm(histoHadIntQualChi2NDOFData,integralChi2Data);
	GammaScalingHistogramm(histoHadIntQualErr2DData,integralErr2DData);
	GammaScalingHistogramm(histoHadIntQualErr3DData,integralErr3DData);
	GammaScalingHistogramm(histoHadIntQualErrXData,integralErrXData);
	GammaScalingHistogramm(histoHadIntQualErrYData,integralErrYData);
	GammaScalingHistogramm(histoHadIntQualErrZData,integralErrZData);
	
	GammaScalingHistogramm(histoHadIntQualChi2NDOFAftCutsData,1./histoHadIntQualChi2NDOFAftCutsData->Integral());
	GammaScalingHistogramm(histoHadIntQualErr2DAftCutsData,1./histoHadIntQualErr2DAftCutsData->Integral());
	GammaScalingHistogramm(histoHadIntQualErr3DAftCutsData,1./histoHadIntQualErr3DAftCutsData->Integral());
	GammaScalingHistogramm(histoHadIntQualErrXAftCutsData,1./histoHadIntQualErrXAftCutsData->Integral());
	GammaScalingHistogramm(histoHadIntQualErrYAftCutsData,1./histoHadIntQualErrZAftCutsData->Integral());
	GammaScalingHistogramm(histoHadIntQualErrZAftCutsData,1./histoHadIntQualErrZAftCutsData->Integral());
	
	
	TH1D*	histoMappingHadIntBPProjectionX0Data;
	histoMappingHadIntBPProjectionX0Data=	histoHadIntVsXYInnerBPData->ProjectionY( "histoMappingHadIntBPProjectionX0Data", histoHadIntVsXYInnerBPData->GetXaxis()->FindBin(0.)-1, histoHadIntVsXYInnerBPData->GetXaxis()->FindBin(0.)+1);
	ConvGammaRebinWithBinCorrection(histoMappingHadIntBPProjectionX0Data,1);
	GammaScalingHistogramm(histoMappingHadIntBPProjectionX0Data,1./numberGoodEventsData);

	TH1D*	histoMappingHadIntBPProjectionY0Data;
	histoMappingHadIntBPProjectionY0Data=		histoHadIntVsXYInnerBPData->ProjectionX("histoMappingHadIntBPProjectionY0Data", histoHadIntVsXYInnerBPData->GetYaxis()->FindBin(0.)-1, histoHadIntVsXYInnerBPData->GetYaxis()->FindBin(0.)+1);
	ConvGammaRebinWithBinCorrection(histoMappingHadIntBPProjectionY0Data,1);
	GammaScalingHistogramm(histoMappingHadIntBPProjectionY0Data,1./numberGoodEventsData);
	
	GammaScalingHistogramm(histoHadIntVsXYInnerBPData,1./numberGoodEventsData);
	
	// ------------------------------------- second file-----------------------------------------------------------
	TDirectory* directorySecHadIntMonteCarlo = 	new TDirectory(); // definition of first folder / list    
	TList*	listHistosSecHadIntMonteCarlo = 		new TList(); // definition of first folder / list
	TList* 	listHadIntContainerMonteCarlo = 			new TList();  // definition of following folder / list	
	TFile* montecarlo = new TFile(MCfile); 
	
	directorySecHadIntMonteCarlo = 			(TDirectory*)montecarlo->Get(nameGammaDirectory); 
	TH1F* 	histoEventQualityMonteCarlo=					(TH1F*)directorySecHadIntMonteCarlo->Get("ESD_EventQuality");
	TH2F* 	histoHadIntVsZRMonteCarlo=			(TH2F*)directorySecHadIntMonteCarlo->Get(Form("ESD_HadIntMap_ZR"));  
	TH2F* 	histoHadIntVsXYMonteCarlo=			(TH2F*)directorySecHadIntMonteCarlo->Get(Form("ESD_HadIntMap_XY"));  
	TH2F*		histoHadIntVsXYInnerBPMonteCarlo=			(TH2F*)directorySecHadIntMonteCarlo->Get(Form("ESD_HadIntMapInnerBeampipe_XY")); 
	TH1D* 	histoHadIntVsZMonteCarlo=			(TH1D*)histoHadIntVsZRMonteCarlo->ProjectionX("histoHadIntVsZMonteCarlo");
	ConvGammaRebinWithBinCorrection(histoHadIntVsZMonteCarlo,rebinIntPlots);
	TH1D* 	histoHadIntVsRMonteCarlo=			(TH1D*)histoHadIntVsZRMonteCarlo->ProjectionY("histoHadIntVsRMonteCarlo");
	ConvGammaRebinWithBinCorrection(histoHadIntVsRMonteCarlo,rebinIntPlots);
	TH1D* 	histoHadIntVsRCentralMonteCarlo=			(TH1D*)histoHadIntVsZRMonteCarlo->ProjectionY("histoHadIntVsRCentralMonteCarlo", histoHadIntVsZRMonteCarlo->GetXaxis()->FindBin(-10.),histoHadIntVsZRMonteCarlo->GetXaxis()->FindBin(10.));
// 	ConvGammaRebinWithBinCorrection(histoHadIntVsRMonteCarlo,rebinIntPlots);
	Double_t integralBPMC = histoHadIntVsRMonteCarlo->Integral(histoHadIntVsRMonteCarlo->GetXaxis()->FindBin(2.5),histoHadIntVsRMonteCarlo->GetXaxis()->FindBin(3.5));
	Double_t integralGasMC = histoHadIntVsRMonteCarlo->Integral(histoHadIntVsRMonteCarlo->GetXaxis()->FindBin(90),histoHadIntVsRMonteCarlo->GetXaxis()->FindBin(180));
	
	cout <<"Number of events in PP:" << integralBPMC << endl;
	cout <<"Number of events in Gas:" << integralGasMC << endl;
	
	TH2F* 	histoHadIntVsZRMonteCarlo2 =	(TH2F*) histoHadIntVsZRMonteCarlo->Clone("histoHadIntVsZRMonteCarlo2");
	TH1F* 	histoHadIntQualChi2NDOFMonteCarlo=	(TH1F*)directorySecHadIntMonteCarlo->Get("ESD_HadIntQual_Chi2PerDOF"); 
	TH1F* 	histoHadIntQualErrXMonteCarlo=	(TH1F*)directorySecHadIntMonteCarlo->Get("ESD_HadIntQual_ErrX"); 
	TH1F* 	histoHadIntQualErrYMonteCarlo=	(TH1F*)directorySecHadIntMonteCarlo->Get("ESD_HadIntQual_ErrY"); 
	TH1F* 	histoHadIntQualErrZMonteCarlo=	(TH1F*)directorySecHadIntMonteCarlo->Get("ESD_HadIntQual_ErrZ"); 
	TH1F* 	histoHadIntQualErr2DMonteCarlo=	(TH1F*)directorySecHadIntMonteCarlo->Get("ESD_HadIntQual_Err2D"); 
	TH1F* 	histoHadIntQualErr3DMonteCarlo=	(TH1F*)directorySecHadIntMonteCarlo->Get("ESD_HadIntQual_Err3D"); 
	TH1F* 	histoHadIntQualNTracksMonteCarlo=	(TH1F*)directorySecHadIntMonteCarlo->Get("ESD_HadIntQual_nTracksSecVtx"); 
	TH1F* 	histoHadIntQualChi2NDOFAftCutsMonteCarlo=	(TH1F*)directorySecHadIntMonteCarlo->Get(Form("ESD_HadIntQualAftCuts_Chi2PerDOF"));   
	TH1F* 	histoHadIntQualErrXAftCutsMonteCarlo=	(TH1F*)directorySecHadIntMonteCarlo->Get(Form("ESD_HadIntQualAftCuts_ErrX"));   
	TH1F* 	histoHadIntQualErrYAftCutsMonteCarlo=	(TH1F*)directorySecHadIntMonteCarlo->Get(Form("ESD_HadIntQualAftCuts_ErrY"));   
	TH1F* 	histoHadIntQualErrZAftCutsMonteCarlo=	(TH1F*)directorySecHadIntMonteCarlo->Get(Form("ESD_HadIntQualAftCuts_ErrZ"));   
	TH1F* 	histoHadIntQualErr2DAftCutsMonteCarlo=	(TH1F*)directorySecHadIntMonteCarlo->Get(Form("ESD_HadIntQualAftCuts_Err2D"));   
	TH1F* 	histoHadIntQualErr3DAftCutsMonteCarlo=	(TH1F*)directorySecHadIntMonteCarlo->Get(Form("ESD_HadIntQualAftCuts_Err3D"));  
	TH1F* 	histoHadIntMapEtaMonteCarlo=	(TH1F*)directorySecHadIntMonteCarlo->Get(Form("ESD_HadIntMap_Eta"));   
		
	for (Int_t iX = 1 ; iX < histoHadIntVsXYInnerBPMonteCarlo->GetNbinsX(); iX ++){
		for (Int_t iY = 1; iY < histoHadIntVsXYInnerBPMonteCarlo->GetNbinsY(); iY ++) {
			Double_t rcalc = TMath::Sqrt(histoHadIntVsXYInnerBPMonteCarlo->GetXaxis()->GetBinCenter(iX)*histoHadIntVsXYInnerBPMonteCarlo->GetXaxis()->GetBinCenter(iX) + histoHadIntVsXYInnerBPMonteCarlo->GetYaxis()->GetBinCenter(iY)*histoHadIntVsXYInnerBPMonteCarlo->GetYaxis()->GetBinCenter(iY));
			if (rcalc < rmin){
				histoHadIntVsXYInnerBPMonteCarlo->SetBinContent(iX,iY,NULL);
			}
		}
	}

	
	Float_t numberGoodEventsMonteCarlo = 						histoEventQualityMonteCarlo->GetEntries()-histoEventQualityMonteCarlo->GetBinContent(6)-histoEventQualityMonteCarlo->GetBinContent(4);
	Float_t numberGoodTriggerMonteCarlo = 					histoEventQualityMonteCarlo->GetEntries();
	Float_t numberReconstGammaMonteCarlo = 					histoHadIntVsRMonteCarlo->GetEntries();
	cout<< MCfile << "    Number of events::   " << numberGoodEventsMonteCarlo << "    Number of triggers::   " << numberGoodTriggerMonteCarlo << "    Number reconstructed vertices::    "<< numberReconstGammaMonteCarlo <<endl;
	
	Double_t meanMultiplitcityMonteCarlo = 1;
	Float_t normFactorReconstMonteCarlo=		1./integralGasMC;// * 1./meanMultiplitcityMonteCarlo;
	Double_t integralNTracksMonteCarlo = 1./histoHadIntQualNTracksMonteCarlo->Integral();
	Double_t integralChi2MonteCarlo = 1./histoHadIntQualChi2NDOFMonteCarlo->Integral();
	Double_t integralErr2DMonteCarlo = 1./histoHadIntQualErr2DMonteCarlo->Integral();
	Double_t integralErr3DMonteCarlo = 1./histoHadIntQualErr3DMonteCarlo->Integral();
	Double_t integralErrXMonteCarlo = 1./histoHadIntQualErrXMonteCarlo->Integral();
	Double_t integralErrYMonteCarlo = 1./histoHadIntQualErrYMonteCarlo->Integral();
	Double_t integralErrZMonteCarlo = 1./histoHadIntQualErrZMonteCarlo->Integral();
	//Scaling reconstr.
	GammaScalingHistogramm(histoHadIntVsRMonteCarlo,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoHadIntMapEtaMonteCarlo,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoHadIntVsZMonteCarlo,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoHadIntVsRCentralMonteCarlo,normFactorReconstMonteCarlo);
 	GammaScalingHistogramm(histoHadIntVsZRMonteCarlo2,1./numberGoodEventsMonteCarlo);
 	GammaScalingHistogramm(histoHadIntQualNTracksMonteCarlo,integralNTracksMonteCarlo);
// 	GammaScalingHistogramm(histoHadIntVsXYMonteCarlo,normFactorReconstMonteCarlo);
// 	GammaScalingHistogramm(histoHadIntVsXYHotZoneMonteCarlo,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoHadIntQualChi2NDOFMonteCarlo,integralChi2MonteCarlo);
	GammaScalingHistogramm(histoHadIntQualErr2DMonteCarlo,integralErr2DMonteCarlo);
	GammaScalingHistogramm(histoHadIntQualErr3DMonteCarlo,integralErr3DMonteCarlo);
	GammaScalingHistogramm(histoHadIntQualErrXMonteCarlo,integralErrXMonteCarlo);
	GammaScalingHistogramm(histoHadIntQualErrYMonteCarlo,integralErrYMonteCarlo);
	GammaScalingHistogramm(histoHadIntQualErrZMonteCarlo,integralErrZMonteCarlo);
	
	GammaScalingHistogramm(histoHadIntQualChi2NDOFAftCutsMonteCarlo,1./histoHadIntQualChi2NDOFAftCutsMonteCarlo->Integral());
	GammaScalingHistogramm(histoHadIntQualErr2DAftCutsMonteCarlo,1./histoHadIntQualErr2DAftCutsMonteCarlo->Integral());
	GammaScalingHistogramm(histoHadIntQualErr3DAftCutsMonteCarlo,1./histoHadIntQualErr3DAftCutsMonteCarlo->Integral());
	GammaScalingHistogramm(histoHadIntQualErrXAftCutsMonteCarlo,1./histoHadIntQualErrXAftCutsMonteCarlo->Integral());
	GammaScalingHistogramm(histoHadIntQualErrYAftCutsMonteCarlo,1./histoHadIntQualErrYAftCutsMonteCarlo->Integral());
	GammaScalingHistogramm(histoHadIntQualErrZAftCutsMonteCarlo,1./histoHadIntQualErrZAftCutsMonteCarlo->Integral());
	
	TH1D*	histoMappingHadIntBPProjectionX0MonteCarlo;
	histoMappingHadIntBPProjectionX0MonteCarlo=		histoHadIntVsXYInnerBPMonteCarlo->ProjectionY( "histoMappingHadIntBPProjectionX0MonteCarlo", histoHadIntVsXYInnerBPMonteCarlo->GetXaxis()->FindBin(0.)-1, histoHadIntVsXYInnerBPMonteCarlo->GetXaxis()->FindBin(0.)+1);
	ConvGammaRebinWithBinCorrection(histoMappingHadIntBPProjectionX0MonteCarlo,1.);
	GammaScalingHistogramm(histoMappingHadIntBPProjectionX0MonteCarlo,1./numberGoodEventsMonteCarlo);

	TH1D*	histoMappingHadIntBPProjectionY0MonteCarlo;
	histoMappingHadIntBPProjectionY0MonteCarlo=		histoHadIntVsXYInnerBPMonteCarlo->ProjectionX("histoMappingHadIntBPProjectionY0MonteCarlo", histoHadIntVsXYInnerBPMonteCarlo->GetYaxis()->FindBin(0.)-1, histoHadIntVsXYInnerBPMonteCarlo->GetYaxis()->FindBin(0.)+1);
	ConvGammaRebinWithBinCorrection(histoMappingHadIntBPProjectionY0MonteCarlo,1.);
	GammaScalingHistogramm(histoMappingHadIntBPProjectionY0MonteCarlo,1./numberGoodEventsMonteCarlo);
	
	GammaScalingHistogramm(histoHadIntVsXYInnerBPMonteCarlo,1./numberGoodEventsMonteCarlo);
	
	
	
	// -----------------page 1---------------------------------------------------------------------------
	TLine * 	linePhi =	 		new TLine (-3.2,1,3.2,1);
	TLine * 	lineZ = 			new TLine (-300,1,300,1);
	TLine * 	lineR = 			new TLine (0,1,200,1);
	linePhi->SetLineColor(2);
	lineZ->SetLineColor(2);
	lineR->SetLineColor(2);

	if (suffix.CompareTo("png")==0){
		PlotStandard2D( histoHadIntVsXYData , Form("%s/XY_distribution%s.%s",outputDirectory.Data(),suffixData.Data(),suffix.Data()), "", "X (cm)", "Y (cm)", kTRUE, -100., 100., kTRUE, -100., 100.,kTRUE, 1, histoHadIntVsXYData->GetMaximum(), 0, 1, floatLocationRightUp2D,3000,2640,rightMargin, leftMargin,0.08, 0.02,"Data");
	} else {
		PlotStandard2D( histoHadIntVsXYData , Form("%s/XY_distribution%s.%s",outputDirectory.Data(),suffixData.Data(),suffix.Data()), "", "X (cm)", "Y (cm)", kTRUE, -100., 100., kTRUE, -100., 100.,kTRUE, 1, histoHadIntVsXYData->GetMaximum(), 0, 1, floatLocationRightUp2D,1000,880,rightMargin, leftMargin,0.08, 0.02,"Data");
	}
	TCanvas * canvas2DimensionXYMC = new TCanvas("canvas2DimensionXYMC","",1000,880);  // gives the page size
	canvas2DimensionXYMC->SetLogz(1);
	canvas2DimensionXYMC->SetTopMargin(0.02);                                              
	canvas2DimensionXYMC->SetRightMargin(rightMargin);                                            
	canvas2DimensionXYMC->SetLeftMargin(leftMargin);
	canvas2DimensionXYMC->SetBottomMargin(0.08);
	canvas2DimensionXYMC->cd();
	histoHadIntVsXYMonteCarlo->GetZaxis()->SetRangeUser(1., histoHadIntVsXYMonteCarlo->GetMaximum());
	DrawAutoGammaHisto2D(   histoHadIntVsXYMonteCarlo,
									"", "X (cm)", "Y (cm)", "",
									kTRUE, -100., 100.,
									kTRUE, -100., 100.);
// 	DrawStructure();
// 	DrawAliceLogoPerformance2D(floatLocationRightUp2DXY[0],floatLocationRightUp2DXY[1],floatLocationRightUp2DXY[2],floatLocationRightUp2DXY[3],0.02,textDate,collisionSystem, "", textPeriod,1000,880);
	canvas2DimensionXYMC->Update();
	canvas2DimensionXYMC->SaveAs(Form("%s/XY_distribution%s.%s",outputDirectory.Data(),suffixMC.Data(),suffix.Data()));
	delete canvas2DimensionXYMC;
	delete histoHadIntVsXYData;
	delete histoHadIntVsXYMonteCarlo;

		// -----------------------  2 dim Plots -----------------------------------
// 	histoHadIntVsZRData2->GetZaxis()->SetRangeUser(1,histoHadIntVsZRData2->GetMaximum());
	if (suffix.CompareTo("png")==0){
		cout << histoHadIntVsZRData2->GetZaxis()->GetLabelOffset()<< endl;
		PlotStandard2D( histoHadIntVsZRData2 , Form("%s/ZR_distribution%s.%s",outputDirectory.Data(),suffixData.Data(),suffix.Data()), "",  "Z (cm)", "R (cm)",kTRUE,0., 90., kFALSE, -100., 100., kTRUE, histoHadIntVsZRData2->GetMinimum()*10, histoHadIntVsZRData2->GetMaximum()*2, 0, 1, floatLocationLeftLow2DZR,3000,1320,0.08, 0.05,0.08, 0.02,"Data",1.2,0.6,-0.004);

			// -----------------------  2 dim Plots -----------------------------------
	// 	histoHadIntVsZRMonteCarlo2->GetZaxis()->SetRangeUser(1,histoHadIntVsZRMonteCarlo2->GetMaximum());
		PlotStandard2D( histoHadIntVsZRMonteCarlo2 , Form("%s/ZR_distribution%s.%s",outputDirectory.Data(),suffixMC.Data(),suffix.Data()), "",  "Z (cm)", "R (cm)",kTRUE,0., 90., kFALSE, -100., 100., kTRUE,histoHadIntVsZRData2->GetMinimum()*10, histoHadIntVsZRData2->GetMaximum(), 0, 1, floatLocationLeftLow2DZR,3000,1320,0.08, 0.05,0.08, 0.02,"",1.2,0.6,-0.004);
	} else {
		cout << histoHadIntVsZRData2->GetZaxis()->GetLabelOffset()<< endl;
		PlotStandard2D( histoHadIntVsZRData2 , Form("%s/ZR_distribution%s.%s",outputDirectory.Data(),suffixData.Data(),suffix.Data()), "",  "Z (cm)", "R (cm)",kTRUE,0., 90., kFALSE, -100., 100., kTRUE, histoHadIntVsZRData2->GetMinimum()*10, histoHadIntVsZRData2->GetMaximum()*2, 0, 1, floatLocationLeftLow2DZR,1000,440,0.08, 0.05,0.08, 0.02,"Data",1.2,0.6,-0.004);

			// -----------------------  2 dim Plots -----------------------------------
	// 	histoHadIntVsZRMonteCarlo2->GetZaxis()->SetRangeUser(1,histoHadIntVsZRMonteCarlo2->GetMaximum());
		PlotStandard2D( histoHadIntVsZRMonteCarlo2 , Form("%s/ZR_distribution%s.%s",outputDirectory.Data(),suffixMC.Data(),suffix.Data()), "",  "Z (cm)", "R (cm)",kTRUE,0., 90., kFALSE, -100., 100., kTRUE,histoHadIntVsZRData2->GetMinimum()*10, histoHadIntVsZRData2->GetMaximum(), 0, 1, floatLocationLeftLow2DZR,1000,440,0.08, 0.05,0.08, 0.02,"",1.2,0.6,-0.004);
	}
	
	delete histoHadIntVsZRData2;
	delete histoHadIntVsZRMonteCarlo2;

	if (suffix.CompareTo("png")==0){
		PlotStandard2D( histoHadIntVsXYInnerBPData , Form("%s/XY_distribution_BP_ZSmaller30%s.%s",outputDirectory.Data(),suffixData.Data(),suffix.Data()), "", "X (cm)", "Y (cm)", kTRUE, -15., 15., kTRUE, -15., 15.,kTRUE, histoHadIntVsXYInnerBPData->GetMinimum()*10, histoHadIntVsXYInnerBPData->GetMaximum(), 0, 1, floatLocationRightUp2D,3000,2640,rightMargin, leftMargin,0.08, 0.02,"Data");
		PlotStandard2D( histoHadIntVsXYInnerBPMonteCarlo , Form("%s/XY_distribution_BP_ZSmaller30%s.%s",outputDirectory.Data(),suffixMC.Data(),suffix.Data()), "", "X (cm)", "Y (cm)", kTRUE, -15., 15., kTRUE, -15., 15.,kTRUE,   histoHadIntVsXYInnerBPData->GetMinimum()*10, histoHadIntVsXYInnerBPData->GetMaximum(), 0, 1, floatLocationRightUp2D,3000,2640,rightMargin, leftMargin,0.08, 0.02,"");
	} else {
		PlotStandard2D( histoHadIntVsXYInnerBPData , Form("%s/XY_distribution_BP_ZSmaller30%s.%s",outputDirectory.Data(),suffixData.Data(),suffix.Data()), "", "X (cm)", "Y (cm)", kTRUE, -15., 15., kTRUE, -15., 15.,kTRUE, histoHadIntVsXYInnerBPData->GetMinimum()*10, histoHadIntVsXYInnerBPData->GetMaximum(), 0, 1, floatLocationRightUp2D,1000,880,rightMargin, leftMargin,0.08, 0.02,"Data");
		PlotStandard2D( histoHadIntVsXYInnerBPMonteCarlo , Form("%s/XY_distribution_BP_ZSmaller30%s.%s",outputDirectory.Data(),suffixMC.Data(),suffix.Data()), "", "X (cm)", "Y (cm)", kTRUE, -15., 15., kTRUE, -15., 15.,kTRUE,   histoHadIntVsXYInnerBPData->GetMinimum()*10, histoHadIntVsXYInnerBPData->GetMaximum(), 0, 1, floatLocationRightUp2D,1000,880,rightMargin, leftMargin,0.08, 0.02,"");
	}
	delete histoHadIntVsXYInnerBPData;
	delete histoHadIntVsXYInnerBPMonteCarlo;

// --------------------------- single plots  ----------------------------
	
	// ---------------------------- R- Distribution -----------------------
	TCanvas * canvasRLin = new TCanvas("canvasRLin","",1200,1000);  // gives the page size  
	canvasRLin->SetLogy(0);
	canvasRLin->SetTopMargin(0.05);            
	canvasRLin->SetLeftMargin(0.13);            
	canvasRLin->cd();
	
	DrawAutoGammaHistosWOLeg( histoHadIntVsRData, 
							histoHadIntVsRMonteCarlo, 
							"","R (cm)",textYAxisRHisto,
							kTRUE, 1.1,minYValueRPlot,
							kFALSE,0. ,0.,
							kTRUE, 0.,180.);
	histoHadIntVsRData->Draw("same,hist,e");		
	histoHadIntVsRMonteCarlo->Draw("same,hist,e");
	histoHadIntVsRData->Draw("same,axis");								    

	TLegend* legendRPlot;
	legendRPlot = new TLegend( 0.68,0.88,0.96,0.94);
	legendRPlot->SetTextSize(0.02);                        
	legendRPlot->SetFillColor(0);
	legendRPlot->SetBorderSize(0);
	if (!kPureData){
		legendRPlot->AddEntry(histoHadIntVsRData,("Data"),"l");
		legendRPlot->AddEntry(histoHadIntVsRMonteCarlo,("MC"),"l");
	} else {
		legendRPlot->AddEntry(histoHadIntVsRData,("pass2"),"l");
		legendRPlot->AddEntry(histoHadIntVsRMonteCarlo,("pass4"),"l");
	}
	legendRPlot->Draw();
	
	Double_t vectorLinPlot[9] = {7.e-4, 7.e-4,4.e-4,4.e-4,5.4e-4,5.4e-4,4.8e-4,4.2e-4,0.7e-4};
	
// 	if(etaCutNumber.CompareTo("5") != 0  && etaCutNumber.CompareTo("6") != 0) DrawIndividualTextSlicesR (vectorLinPlot , 0.025,"lin");
	
// 	DrawAliceLogoPerformance(floatLocationRightUp[0],floatLocationRightUp[1],floatLocationRightUp[2],floatLocationRightUp[3],0.0001, textDate,collisionSystem, textGenerator, textPeriod,1200,1000); 
	canvasRLin->Update();
	canvasRLin->SaveAs(Form("%s/R_distribution_lin.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasRLin;          
	
	TCanvas * canvasRLog = new TCanvas("canvasRLog","",1200,1000);  // gives the page size
	canvasRLog->cd();
	canvasRLog->SetLeftMargin(0.13);
	//canvasRLog->SetTopMargin(0.05);            
	canvasRLog->SetLogy(1);
	
	DrawAutoGammaHistosWOLeg( histoHadIntVsRData, 
							histoHadIntVsRMonteCarlo, 
							"","R (cm)",textYAxisRHisto,
							kTRUE, 1.5,minYValueRPlot,
							kFALSE,0. ,0.,
							kTRUE, 0.,180.);		
	histoHadIntVsRData->Draw("same,hist,e");		
	histoHadIntVsRMonteCarlo->Draw("same,hist,e");
	histoHadIntVsRData->Draw("same,axis");		
	
	Double_t vectorLogPlot[9] = {1.8e-5,2e-5,3e-6,3e-6,1.5e-6,1.5e-6,1.5e-6,1.5e-6,1.5e-6};
// 	if(etaCutNumber.CompareTo("5") != 0  && etaCutNumber.CompareTo("6") != 0) DrawIndividualTextSlicesR (vectorLogPlot , 0.025,"log");
	
	legendRPlot->Draw();
// 	DrawAliceLogoPerformance(floatLocationRightUp[0],floatLocationRightUp[1],floatLocationRightUp[2],floatLocationRightUp[3],0.0001, textDate,collisionSystem, textGenerator, textPeriod,1200,1000); 
	
	canvasRLog->Update();
	canvasRLog->SaveAs(Form("%s/R_distribution_log.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasRLog;
	
	TCanvas * canvasRLogCentral = new TCanvas("canvasRLogCentral","",1200,1000);  // gives the page size
	canvasRLogCentral->cd();
	canvasRLogCentral->SetLeftMargin(0.13);
	//canvasRLogCentral->SetTopMargin(0.05);            
	canvasRLogCentral->SetLogy(1);
	
	DrawAutoGammaHistosWOLeg( histoHadIntVsRCentralData, 
							histoHadIntVsRCentralMonteCarlo, 
							"","R (cm)",textYAxisRHisto,
							kTRUE, 1.5,minYValueRPlot,
							kFALSE,0. ,0.,
							kTRUE, 0.,180.);		
// 	histoHadIntVsRCentralData->Draw("same,hist,e");		
// 	histoHadIntVsRCentralMonteCarlo->Draw("same,hist,e");
// 	histoHadIntVsRCentralData->Draw("same,axis");		
// 	
// 	if(etaCutNumber.CompareTo("5") != 0  && etaCutNumber.CompareTo("6") != 0) DrawIndividualTextSlicesR (vectorLogPlot , 0.025,"log");
	
	legendRPlot->Draw();
// 	DrawAliceLogoPerformance(floatLocationRightUp[0],floatLocationRightUp[1],floatLocationRightUp[2],floatLocationRightUp[3],0.0001, textDate,collisionSystem, textGenerator, textPeriod,1200,1000); 
	
	canvasRLogCentral->Update();
	canvasRLogCentral->SaveAs(Form("%s/R_distributionCentral_log.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasRLogCentral;
	
	TCanvas * canvasREnlargedSDDSSDLog = new TCanvas("canvasREnlargedSDDSSDLog","",1200,1000);  // gives the page size
	canvasREnlargedSDDSSDLog->cd();
	canvasREnlargedSDDSSDLog->SetLeftMargin(0.13);
//	canvasREnlargedSDDSSDLog->SetTopMargin(0.05);            
	canvasREnlargedSDDSSDLog->SetLogy(1);
	
	DrawAutoGammaHistosWOLeg( histoHadIntVsRData, 
							histoHadIntVsRMonteCarlo, 
							"","R (cm)",textYAxisRHisto,
							kTRUE, 1.8,minYValueRPlotSPD,
							kFALSE,1.e-2 ,minYValueRPlotSPD,
							kTRUE, 13.,60.);
	histoHadIntVsRData->Draw("hist");
	histoHadIntVsRData->Draw("same,hist,e");		
	histoHadIntVsRMonteCarlo->Draw("same,hist,e");
	histoHadIntVsRData->Draw("same,axis");											
	legendRPlot->Draw();
// 	if(etaCutNumber.CompareTo("5") == 0 || etaCutNumber.CompareTo("6") == 0 ){
// 		DrawAliceLogoPerformance(floatLocationLeftUp[0],floatLocationLeftUp[1],floatLocationLeftUp[2],floatLocationLeftUp[3],0.0001, textDate,collisionSystem, textGenerator, textPeriod,1200,1000);     
// 	} else {
// 		DrawAliceLogoPerformance(floatLocationRightUp[0],floatLocationRightUp[1],floatLocationRightUp[2],floatLocationRightUp[3],0.0001, textDate,collisionSystem, textGenerator, textPeriod,1200,1000); 
// 	}
	
	canvasREnlargedSDDSSDLog->Update();
	canvasREnlargedSDDSSDLog->SaveAs(Form("%s/R_distributionSDDSSDEnlarged_log.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasREnlargedSDDSSDLog;
	
	
	TCanvas * canvasREnlargedTPCLog = new TCanvas("canvasREnlargedTPCLog","",1200,1000);  // gives the page size
	canvasREnlargedTPCLog->cd();
	canvasREnlargedTPCLog->SetLeftMargin(0.13);
//	canvasREnlargedTPCLog->SetTopMargin(0.05);            
	canvasREnlargedTPCLog->SetLogy(1);
	
	DrawAutoGammaHistosWOLeg( histoHadIntVsRData, 
							histoHadIntVsRMonteCarlo, 
							"","R (cm)",textYAxisRHisto,
							kTRUE, 1.,minYValueRPlot,
							kFALSE,0. ,0.,
							kTRUE, 60.,180.);
	histoHadIntVsRData->Draw("same,hist,e");		
	histoHadIntVsRMonteCarlo->Draw("same,hist,e");
	histoHadIntVsRData->Draw("same,axis");											
	
	legendRPlot->Draw();
// 	DrawAliceLogoPerformance(floatLocationRightUp[0],floatLocationRightUp[1],floatLocationRightUp[2],floatLocationRightUp[3],0.0001, textDate,collisionSystem, textGenerator, textPeriod,1200,1000); 
	
	canvasREnlargedTPCLog->Update();
	canvasREnlargedTPCLog->SaveAs(Form("%s/R_distributionTPCEnlarged_log.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasREnlargedTPCLog;
	
	TCanvas * canvasREnlargedSPDLog = new TCanvas("canvasREnlargedSPDLog","",1200,1000);  // gives the page size
	canvasREnlargedSPDLog->cd();
	canvasREnlargedSPDLog->SetLeftMargin(0.13);
//	canvasREnlargedSPDLog->SetTopMargin(0.05);            
	canvasREnlargedSPDLog->SetLogy(1);
	
	DrawAutoGammaHistosWOLeg( histoHadIntVsRData, 
							histoHadIntVsRMonteCarlo, 
							"","R (cm)",textYAxisRHisto,
							kTRUE, 3.,minYValueRPlotSPD,
							kFALSE,3.e-4 ,minYValueRPlotSPD,
							kTRUE, 0.,13.);
	histoHadIntVsRData->Draw("hist");
	histoHadIntVsRData->Draw("same,hist,e");		
	histoHadIntVsRMonteCarlo->Draw("same,hist,e");
	histoHadIntVsRData->Draw("same,axis");		
	legendRPlot->Draw();
// 	DrawAliceLogoPerformance(floatLocationLeftUp[0],floatLocationLeftUp[1],floatLocationLeftUp[2],floatLocationLeftUp[3],0.0001, textDate,collisionSystem, textGenerator, textPeriod,1200,1000);     
	
	canvasREnlargedSPDLog->Update();
	canvasREnlargedSPDLog->SaveAs(Form("%s/R_distributionSPDEnlarged_log.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasREnlargedSPDLog;
	delete legendRPlot;

	//-------------------- BP --------------------------------------------

	TCanvas * canvasBeampipeX0 = new TCanvas("canvasBeampipeX0","",1200,1000);  // gives the page size
	canvasBeampipeX0->cd();
	canvasBeampipeX0->SetLeftMargin(0.13);
	canvasBeampipeX0->SetTopMargin(0.04);            
// 	canvasBeampipeX0->SetLogy(1);
	
	DrawAutoGammaHistosMaterialP( histoMappingHadIntBPProjectionX0Data, 
							histoMappingHadIntBPProjectionX0MonteCarlo, 
							"","Y (cm)","#frac{dN_{sec. had}}{N_{ch}dY}",
							kTRUE, 1.3,2e-4,
							kFALSE,0. ,0.,
							kTRUE, -4,4.);		
// 	DrawGammaSetMarker(histoMappingHadIntBPProjectionX0Data, 20, 0.8, kBlack, kBlack);
// 	histoMappingHadIntBPProjectionX0Data->Draw("p");		
// 	DrawGammaSetMarker(histoMappingHadIntBPProjectionX0MonteCarlo, 20, 0.8, kRed+2, kRed+2);
// 	histoMappingHadIntBPProjectionX0MonteCarlo->Draw("same,p");
// 	DrawAliceLogoPerformance(floatLocationRightUp[0],floatLocationRightUp[1],floatLocationRightUp[2],floatLocationRightUp[3],0.0001, textDate,collisionSystem, textGenerator, textPeriod,1200,1000); 
	
	canvasBeampipeX0->Update();
	canvasBeampipeX0->SaveAs(Form("%s/Y_distribution_XEqual0.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasBeampipeX0;

	
	TCanvas * canvasBeampipeY0 = new TCanvas("canvasBeampipeY0","",1200,1000);  // gives the page size
	canvasBeampipeY0->cd();
	canvasBeampipeY0->SetLeftMargin(0.13);
	canvasBeampipeY0->SetTopMargin(0.04);            
// 	canvasBeampipeY0->SetLogy(1);
	
	DrawAutoGammaHistosMaterialP( histoMappingHadIntBPProjectionY0Data, 
							histoMappingHadIntBPProjectionY0MonteCarlo, 
							"","X (cm)","#frac{dN_{sec. had}}{N_{ch}dX}",
							kTRUE, 1.3,2e-4,
							kFALSE,0. ,0.,
							kTRUE, -4.,4.);		
// 	DrawGammaSetMarker(histoMappingHadIntBPProjectionY0Data, 20, 0.8, kBlack, kBlack);
// 	histoMappingHadIntBPProjectionY0Data->Draw("p");		
// 	DrawGammaSetMarker(histoMappingHadIntBPProjectionY0MonteCarlo, 20, 0.8, kRed+2, kRed+2);
// 	histoMappingHadIntBPProjectionY0MonteCarlo->Draw("same,p");
// 	DrawAliceLogoPerformance(floatLocationRightUp[0],floatLocationRightUp[1],floatLocationRightUp[2],floatLocationRightUp[3],0.0001, textDate,collisionSystem, textGenerator, textPeriod,1200,1000); 
	
	canvasBeampipeY0->Update();
	canvasBeampipeY0->SaveAs(Form("%s/X_distribution_YEqual0.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasBeampipeY0;

		
	//-------------------- Z - Distribution ------------------------------
	
	
	TCanvas * canvasZLin = new TCanvas("canvasZLin","",1200,1000);  // gives the page size  
	canvasZLin->SetLogy(0);
	canvasZLin->SetTopMargin(0.04);            
	canvasZLin->SetLeftMargin(0.13);            
	canvasZLin->cd();

	
	DrawAutoGammaHistosWOLeg( histoHadIntVsZData, 
							histoHadIntVsZMonteCarlo, 
							"","Z (cm)",textYAxisZHisto,
							kTRUE, 1.2,minYValueZPlot,
							kFALSE,0. ,0.,
							kTRUE,  -301.,301.);
	histoHadIntVsZData->Draw("same,hist,e");		
	histoHadIntVsZMonteCarlo->Draw("same,hist,e");
	histoHadIntVsZData->Draw("same,axis");								    
	
	TLegend* legendZPlot;
	legendZPlot = new TLegend( 0.145,0.88,0.33,0.94);
	legendZPlot->SetTextSize(0.02);                        
	legendZPlot->SetFillColor(0);
	legendZPlot->SetBorderSize(0);
	if (!kPureData){
		legendZPlot->AddEntry(histoHadIntVsZData,("Data"),"l");
		legendZPlot->AddEntry(histoHadIntVsZMonteCarlo,("MC"),"l");
	} else {
		legendZPlot->AddEntry(histoHadIntVsZData,("pass2"),"l");
		legendZPlot->AddEntry(histoHadIntVsZMonteCarlo,("pass4"),"l");
	}
	legendZPlot->Draw();
// 	DrawAliceLogoPerformance(floatLocationRightUp2[0],floatLocationRightUp2[1],floatLocationRightUp2[2],floatLocationRightUp2[3],0.0001,textDate,collisionSystem, textGenerator, textPeriod,1200,1000);   
	canvasZLin->Update();
	canvasZLin->SaveAs(Form("%s/Z_distribution_lin.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasZLin;          
	

	TCanvas * canvasZLog = new TCanvas("canvasZLog","",1200,1000);  // gives the page size
	canvasZLog->cd();
	canvasZLog->SetLeftMargin(0.13);
	//canvasZLog->SetTopMargin(0.05);            
	canvasZLog->SetLogy(1);

	DrawAutoGammaHistosWOLeg( histoHadIntVsZData, 
							histoHadIntVsZMonteCarlo, 
							"","Z (cm)",textYAxisZHisto,
							kTRUE, 2.,minYValueZPlot,
							kFALSE,0. ,0.,
							kTRUE,  -301.,301.);		
	histoHadIntVsZData->Draw("hist");
	histoHadIntVsZData->Draw("same,hist,e");		
	histoHadIntVsZMonteCarlo->Draw("same,hist,e");
	histoHadIntVsZData->Draw("same,axis");		

	legendZPlot->Draw();
// 	DrawAliceLogoPerformance(floatLocationRightUp2[0],floatLocationRightUp2[1],floatLocationRightUp2[2],floatLocationRightUp2[3],0.0001, textDate,collisionSystem, textGenerator, textPeriod,1200,1000); 
	canvasZLog->Update();
	canvasZLog->SaveAs(Form("%s/Z_distribution_log.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasZLog;		

	TCanvas * canvasEtaLin = new TCanvas("canvasEtaLin","",1200,1000);  // gives the page size  
	canvasEtaLin->SetLogy(0);
	canvasEtaLin->SetTopMargin(0.04);            
	canvasEtaLin->SetLeftMargin(0.13);            
	canvasEtaLin->cd();

	
	DrawAutoGammaHistosWOLeg( histoHadIntMapEtaData, 
							histoHadIntMapEtaMonteCarlo, 
							"","#eta",textYAxisEtaHisto,
							kTRUE, 1.2,minYValueZPlot,
							kFALSE,0. ,0.,
							kTRUE,  -7.,7.);
	histoHadIntMapEtaData->Draw("same,hist,e");		
	histoHadIntMapEtaMonteCarlo->Draw("same,hist,e");
	histoHadIntMapEtaData->Draw("same,axis");								    
	
	TLegend* legendEtaPlot;
	legendEtaPlot = new TLegend( 0.145,0.88,0.33,0.94);
	legendEtaPlot->SetTextSize(0.02);                        
	legendEtaPlot->SetFillColor(0);
	legendEtaPlot->SetBorderSize(0);
	if (!kPureData){
		legendEtaPlot->AddEntry(histoHadIntMapEtaData,("Data"),"l");
		legendEtaPlot->AddEntry(histoHadIntMapEtaMonteCarlo,("MC"),"l");
	} else {
		legendEtaPlot->AddEntry(histoHadIntMapEtaData,("pass2"),"l");
		legendEtaPlot->AddEntry(histoHadIntMapEtaMonteCarlo,("pass4"),"l");
	}
	legendEtaPlot->Draw();
// 	DrawAliceLogoPerformance(floatLocationRightUp2[0],floatLocationRightUp2[1],floatLocationRightUp2[2],floatLocationRightUp2[3],0.0001,textDate,collisionSystem, textGenerator, textPeriod,1200,1000);   
	canvasEtaLin->Update();
	canvasEtaLin->SaveAs(Form("%s/Eta_distribution_lin.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasEtaLin;          

		TCanvas * canvasEtaLog = new TCanvas("canvasEtaLog","",1200,1000);  // gives the page size
	canvasEtaLog->cd();
	canvasEtaLog->SetLeftMargin(0.13);
	//canvasEtaLog->SetTopMargin(0.05);            
	canvasEtaLog->SetLogy(1);

	DrawAutoGammaHistosWOLeg( histoHadIntMapEtaData, 
							histoHadIntMapEtaMonteCarlo, 
							"","#eta",textYAxisEtaHisto,
							kTRUE, 2.,minYValueZPlot,
							kFALSE,0. ,0.,
							kTRUE,  -7.,7.);		
	histoHadIntMapEtaData->Draw("hist");
	histoHadIntMapEtaData->Draw("same,hist,e");		
	histoHadIntMapEtaMonteCarlo->Draw("same,hist,e");
	histoHadIntMapEtaData->Draw("same,axis");		

	legendEtaPlot->Draw();
// 	DrawAliceLogoPerformance(floatLocationRightUp2[0],floatLocationRightUp2[1],floatLocationRightUp2[2],floatLocationRightUp2[3],0.0001, textDate,collisionSystem, textGenerator, textPeriod,1200,1000); 
	canvasEtaLog->Update();
	canvasEtaLog->SaveAs(Form("%s/Eta_distribution_log.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasEtaLog;		

	if (kPureData){
		return;
	}
	//__________________________________________________________________________________________________________________
	// ****************************************** Definition of the projection plots ***********************************
	//__________________________________________________________________________________________________________________
	
	cout<< "here I Am" << endl;
	TH1D*	histoMappingHadIntPhiInRData[nBinsR];
	TH2F* 	histoHadIntRPhiData = (TH2F*)directorySecHadIntData->Get(Form("ESD_HadIntMap_RPhi"));  
	for(Int_t iR = 0; iR < nBinsR; iR++){
		histoMappingHadIntPhiInRData[iR]=		histoHadIntRPhiData->ProjectionY( Form("histoMappingHadIntPhiInRData_%i",iR), histoHadIntRPhiData->GetXaxis()->FindBin(arrayRBins[iR]), histoHadIntRPhiData->GetXaxis()->FindBin(arrayRBins[iR+1]));
		ConvGammaRebinWithBinCorrection(histoMappingHadIntPhiInRData[iR],rebin);
		GammaScalingHistogramm(histoMappingHadIntPhiInRData[iR],normFactorReconstData*1/600*1/(TMath::Abs(histoHadIntRPhiData->GetXaxis()->GetBinUpEdge(histoHadIntRPhiData->GetXaxis()->FindBin(arrayRBins[iR+1]))-histoHadIntRPhiData->GetXaxis()->GetBinLowEdge(histoHadIntRPhiData->GetXaxis()->FindBin(arrayRBins[iR])))));
	}
	cout<< "here I Am" << endl;
	Int_t minimumBinZNormal = 0;
	TH1D*	histoMappingHadIntPhiInZData[nBinsZ];
	TH2F* 	histoHadIntZPhiData = (TH2F*)directorySecHadIntData->Get(Form("ESD_HadIntMap_ZPhi"));  
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingHadIntPhiInZData[iZ]=		histoHadIntZPhiData->ProjectionY( Form("histoMappingHadIntPhiInZData_%i",iZ), histoHadIntZPhiData->GetXaxis()->FindBin(arrayZBins[iZ]), histoHadIntZPhiData->GetXaxis()->FindBin(arrayZBins[iZ+1]));
		if (minimumBinZNormal == 0){
			if (histoMappingHadIntPhiInZData[iZ]->GetEntries() > 0){
				minimumBinZNormal = iZ;
			}
		}
		ConvGammaRebinWithBinCorrection(histoMappingHadIntPhiInZData[iZ],rebin*2);
		GammaScalingHistogramm(histoMappingHadIntPhiInZData[iZ],normFactorReconstData*1/180*1/(TMath::Abs(histoHadIntZPhiData->GetXaxis()->GetBinUpEdge(histoHadIntZPhiData->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoHadIntZPhiData->GetXaxis()->GetBinLowEdge(histoHadIntZPhiData->GetXaxis()->FindBin(arrayZBins[iZ])))));
		
	}
	cout<< "here I Am" << endl;
	TH1D*	histoMappingHadIntRInZData[nBinsZ];
	TH2F* 	histoHadIntZRData = (TH2F*)histoHadIntVsZRData->Clone("histoHadIntVsZRMappingData");
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingHadIntRInZData[iZ]=		histoHadIntZRData->ProjectionY( Form("histoMappingHadIntRInZData_%i",iZ), histoHadIntZRData->GetXaxis()->FindBin(arrayZBins[iZ]), histoHadIntZRData->GetXaxis()->FindBin(arrayZBins[iZ+1]));
		ConvGammaRebinWithBinCorrection(histoMappingHadIntRInZData[iZ],rebin);
		GammaScalingHistogramm(histoMappingHadIntRInZData[iZ],normFactorReconstData*1/(2*TMath::Pi())*1/(TMath::Abs(histoHadIntZRData->GetXaxis()->GetBinUpEdge(histoHadIntZRData->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoHadIntZRData->GetXaxis()->GetBinLowEdge(histoHadIntZRData->GetXaxis()->FindBin(arrayZBins[iZ])))));
	}
	cout<< "here I Am" << endl;
	TH1D*	histoMappingHadIntZInRData[nBinsR];
	for(Int_t iR = 0; iR < nBinsR; iR++){
		histoMappingHadIntZInRData[iR]=		histoHadIntZRData->ProjectionX( Form("histoMappingHadIntZInRData_%i",iR), histoHadIntZRData->GetYaxis()->FindBin(arrayRBins[iR]), histoHadIntZRData->GetYaxis()->FindBin(arrayRBins[iR+1]));
		ConvGammaRebinWithBinCorrection(histoMappingHadIntZInRData[iR],rebin);
		GammaScalingHistogramm(histoMappingHadIntZInRData[iR],normFactorReconstData*1/(2*TMath::Pi())*1/(TMath::Abs(histoHadIntZRData->GetYaxis()->GetBinUpEdge(histoHadIntZRData->GetYaxis()->FindBin(arrayRBins[iR+1]))-histoHadIntZRData->GetYaxis()->GetBinLowEdge(histoHadIntZRData->GetYaxis()->FindBin(arrayRBins[iR])))));
	}
	cout<< "here I Am" << endl;
	Int_t minimumBinZSPD = 0;
	TH1D*		histoMappingHadIntSPDPhiInZData[nBinsZ];
	TH2F* 	histoHadIntSPDZPhiData = (TH2F*)directorySecHadIntData->Get(Form("ESD_HadIntMap_SPD_ZPhi"));  
	TString nameHistoSPDPhiInZData;
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingHadIntSPDPhiInZData[iZ]=		histoHadIntSPDZPhiData->ProjectionY( Form("histoMappingHadIntSPDPhiInZData_%i",iZ), histoHadIntSPDZPhiData->GetXaxis()->FindBin(arrayZBins[iZ]), histoHadIntSPDZPhiData->GetXaxis()->FindBin(arrayZBins[iZ+1]));
		if (minimumBinZSPD == 0){
			if (histoMappingHadIntSPDPhiInZData[iZ]->GetEntries() > 0){
				minimumBinZSPD = iZ;
			}
		}
		ConvGammaRebinWithBinCorrection(histoMappingHadIntSPDPhiInZData[iZ],2*rebin);
		GammaScalingHistogramm(histoMappingHadIntSPDPhiInZData[iZ],normFactorReconstData*1/25.*1/(TMath::Abs(histoHadIntSPDZPhiData->GetXaxis()->GetBinUpEdge(histoHadIntSPDZPhiData->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoHadIntSPDZPhiData->GetXaxis()->GetBinLowEdge(histoHadIntSPDZPhiData->GetXaxis()->FindBin(arrayZBins[iZ])))));		
	}

	cout<< "here I Am" << endl;
	Int_t minimumBinZSDD = 0;
	TH1D*		histMappingHadIntSDDPhiInZData[nBinsZ];
	TH2F* 	histoHadIntSDDZPhiData = (TH2F*)directorySecHadIntData->Get(Form("ESD_HadIntMap_SDD_ZPhi"));  
	TString nameHistoSDDPhiInZData;
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histMappingHadIntSDDPhiInZData[iZ]=		histoHadIntSDDZPhiData->ProjectionY( Form("histMappingHadIntSDDPhiInZData_%i",iZ), histoHadIntSDDZPhiData->GetXaxis()->FindBin(arrayZBins[iZ]), histoHadIntSDDZPhiData->GetXaxis()->FindBin(arrayZBins[iZ+1]));
		if (minimumBinZSDD == 0){
			if (histMappingHadIntSDDPhiInZData[iZ]->GetEntries() > 0){
				minimumBinZSDD = iZ;
			}
		}
		ConvGammaRebinWithBinCorrection(histMappingHadIntSDDPhiInZData[iZ],2*rebin);
		GammaScalingHistogramm(histMappingHadIntSDDPhiInZData[iZ],normFactorReconstData*1/25.*1/(TMath::Abs(histoHadIntSDDZPhiData->GetXaxis()->GetBinUpEdge(histoHadIntSDDZPhiData->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoHadIntSDDZPhiData->GetXaxis()->GetBinLowEdge(histoHadIntSDDZPhiData->GetXaxis()->FindBin(arrayZBins[iZ])))));
		
	}
	
	cout<< "here I Am" << endl;
	Int_t minimumBinZHotZone = 0;
	TH1D*		histoMappingHadIntHotZonePhiInZData[nBinsZHotZone];
	TH2F* 	histoHadIntHotZonePhiZData = (TH2F*)directorySecHadIntData->Get(Form("ESD_HadIntMap_HotZone_ZPhi"));  
	for(Int_t iZ = 0; iZ < nBinsZHotZone; iZ++){
		histoMappingHadIntHotZonePhiInZData[iZ]=		histoHadIntHotZonePhiZData->ProjectionY( Form("histoMappingHadIntHotZonePhiInZData_%i",iZ), histoHadIntHotZonePhiZData->GetXaxis()->FindBin(arrayZBinsHotZone[iZ]), histoHadIntHotZonePhiZData->GetXaxis()->FindBin(arrayZBinsHotZone[iZ+1]));
// 		if (minimumBinZHotZone == 0){
// 			if (histoMappingHadIntHotZonePhiInZData[iZ]->GetEntries() > 0){
// 				minimumBinZHotZone = iZ;
// 			}
// 		}
		ConvGammaRebinWithBinCorrection(histoMappingHadIntHotZonePhiInZData[iZ],2*rebin);
		GammaScalingHistogramm(histoMappingHadIntHotZonePhiInZData[iZ],normFactorReconstData*1/30.*1/(TMath::Abs(histoHadIntHotZonePhiZData->GetXaxis()->GetBinUpEdge(histoHadIntHotZonePhiZData->GetXaxis()->FindBin(arrayZBinsHotZone[iZ+1]))-histoHadIntHotZonePhiZData->GetXaxis()->GetBinLowEdge(histoHadIntHotZonePhiZData->GetXaxis()->FindBin(arrayZBinsHotZone[iZ])))));
	}
	
	cout<< "here I Am" << endl;
	Int_t minimumBinZITSTPC = 0;
	TH1D*		histoMappingHadIntITSTPCPhiInZData[nBinsZ];
	TH2F* 	histoHadIntITSTPCPhiZData = (TH2F*)directorySecHadIntData->Get(Form("ESD_HadIntMap_ITSTPC_ZPhi"));  
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingHadIntITSTPCPhiInZData[iZ]=		histoHadIntITSTPCPhiZData->ProjectionY( Form("histoMappingHadIntITSTPCPhiInZData_%i",iZ), histoHadIntITSTPCPhiZData->GetXaxis()->FindBin(arrayZBins[iZ]), histoHadIntITSTPCPhiZData->GetXaxis()->FindBin(arrayZBins[iZ+1]));
		if (minimumBinZITSTPC == 0){
			if (histoMappingHadIntITSTPCPhiInZData[iZ]->GetEntries() > 0){
				minimumBinZITSTPC = iZ;
			}
		}
		ConvGammaRebinWithBinCorrection(histoMappingHadIntITSTPCPhiInZData[iZ],2*rebin);
		GammaScalingHistogramm(histoMappingHadIntITSTPCPhiInZData[iZ],normFactorReconstData*1/30.*1/(TMath::Abs(histoHadIntITSTPCPhiZData->GetXaxis()->GetBinUpEdge(histoHadIntITSTPCPhiZData->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoHadIntITSTPCPhiZData->GetXaxis()->GetBinLowEdge(histoHadIntITSTPCPhiZData->GetXaxis()->FindBin(arrayZBins[iZ])))));
	}
	
	Int_t nRowsZProjections = 0;
	Int_t nColumnsZProjections = 0;
	Int_t nRowsZSPD = 0;
	Int_t nColumnsZSPD = 0;
	Int_t nRowsZSDD = 0;
	Int_t nColumnsZSDD = 0;
	Int_t nRowsZITSTPC = 0;
	Int_t nColumnsZITSTPC = 0;
	Int_t nRowsZHotZone = 0;
	Int_t nColumnsZHotZone = 0;
	
	if ( (nBinsZ - 2*minimumBinZNormal) < 15){
		nColumnsZProjections = 2;
	} else {
		nColumnsZProjections = 3;
	}
	if ( ((nBinsZ - 2*minimumBinZNormal)%nColumnsZProjections) == 0){
		nRowsZProjections = (nBinsZ - 2*minimumBinZNormal)/nColumnsZProjections ;
	} else {
		nRowsZProjections = (nBinsZ - 2*minimumBinZNormal)/nColumnsZProjections + 1;
	}
	cout << "normal minZBin: " << minimumBinZNormal << "\t col: " << nColumnsZProjections  << "\t row: " << nRowsZProjections << endl;
	
	if ( (nBinsZ - 2*minimumBinZSPD ) < 15){
		nColumnsZSPD = 2;
	} else {
		nColumnsZSPD = 3;
	}
	if ( ((nBinsZ - 2*minimumBinZSPD)%nColumnsZSPD) == 0){
		nRowsZSPD = (nBinsZ - 2*minimumBinZSPD)/nColumnsZSPD ;
	} else {
		nRowsZSPD = (nBinsZ - 2*minimumBinZSPD)/nColumnsZSPD + 1;
	}
	cout << "SPD minZBin: " << minimumBinZSPD << "\t col: " << nColumnsZSPD  << "\t row: " << nRowsZSPD << endl;

	if ( (nBinsZ - 2*minimumBinZSDD ) < 15){
		nColumnsZSDD = 2;
	} else {
		nColumnsZSDD = 3;
	}
	if ( ((nBinsZ - 2*minimumBinZSDD)%nColumnsZSDD) == 0){
		nRowsZSDD = (nBinsZ - 2*minimumBinZSDD)/nColumnsZSDD ;
	} else {
		nRowsZSDD = (nBinsZ - 2*minimumBinZSDD)/nColumnsZSDD + 1;
	}
	cout << "SDD minZBin: " << minimumBinZSDD << "\t col: " << nColumnsZSDD  << "\t row: " << nRowsZSDD << endl;
	
	if ( (nBinsZ - 2*minimumBinZITSTPC ) < 15){
		nColumnsZITSTPC = 2;
	} else {
		nColumnsZITSTPC = 3;
	}
	if ( ((nBinsZ - 2*minimumBinZITSTPC)%nColumnsZITSTPC) == 0){
		nRowsZITSTPC = (nBinsZ - 2*minimumBinZITSTPC)/nColumnsZITSTPC ;
	} else {
		nRowsZITSTPC = (nBinsZ - 2*minimumBinZITSTPC)/nColumnsZITSTPC + 1;
	}
	cout << "ITSTPC minZBin: " << minimumBinZITSTPC << "\t col: " << nColumnsZITSTPC  << "\t row: " << nRowsZITSTPC << endl;

	nColumnsZHotZone = 2;
	
	if ( ((nBinsZHotZone - 2*minimumBinZHotZone)%nColumnsZHotZone) == 0){
		nRowsZHotZone = (nBinsZ - 2*minimumBinZHotZone)/nColumnsZHotZone ;
	} else {
		nRowsZHotZone = (nBinsZ - 2*minimumBinZHotZone)/nColumnsZHotZone + 1;
	}
	cout << "HotZone minZBin: " << minimumBinZHotZone << "\t col: " << nColumnsZHotZone  << "\t row: " << nRowsZHotZone << endl;

	TH1D* 	histoMappingHadIntPhiInRMonteCarlo[nBinsR];
	TH1D* 	histoMappingHadIntZInRMonteCarlo[nBinsR];
	TH1D* 	histoMappingHadIntPhiInZMonteCarlo[nBinsZ];
	TH1D* 	histoMappingHadIntRInZMonteCarlo[nBinsZ];		
	TH1D* 	histoMappingHadIntSPDPhiInZMonteCarlo[nBinsZ];
	TH1D* 	histoMappingHadIntSDDPhiInZMonteCarlo[nBinsZ];
	TH1D* 	histoMappingHadIntITSTPCPhiInZMonteCarlo[nBinsZ];
	TH1D*		histoMappingHadIntHotZonePhiInZMonteCarlo[nBinsZHotZone];
	
		
	TString 	nameHistoRatioPhiInRMonteCarlo;
	TH1F*		histoMappingRatioPhiInR[nBinsR];
	TString 	nameHistoRatioZInRMonteCarlo;
	TH1F* 	histoMappingRatioZInR[nBinsR];
	TString	nameHistoRatioPhiInZMonteCarlo;
	TH1F* 	histoMappingRatioPhiInZ[nBinsZ];
	TString 	nameHistoSPDPhiInZMonteCarlo;
	TString 	nameHistoSDDPhiInZMonteCarlo;
	TString	nameHistoRatioSPDPhiInZMonteCarlo;
	TString	nameHistoRatioSDDPhiInZMonteCarlo;
	TH1F * 	histoMappingRatioSPDPhiInZ[nBinsZ];
	TH1F * 	histoMappingRatioSDDPhiInZ[nBinsZ];
	TString 	nameHistoITSTPCPhiInZMonteCarlo;
	TString	nameHistoRatioITSTPCPhiInZMonteCarlo;
	TH1F * 	histoMappingRatioITSTPCPhiInZ[nBinsZ];
	TString 	nameHistoHotZonePhiInZMonteCarlo;
	TString	nameHistoRatioHotZonePhiInZMonteCarlo;
	TH1F * 	histoMappingRatioHotZonePhiInZ[nBinsZHotZone];
	TString 	nameHistoRatioRInZMonteCarlo;
	TH1F* 	histoMappingRatioRInZ[nBinsZ];
	
	TH2F* 	histoHadIntRPhiMonteCarlo = (TH2F*)directorySecHadIntMonteCarlo->Get(Form("ESD_HadIntMap_RPhi"));  
	for(Int_t iR = 0; iR < nBinsR; iR++){
		histoMappingHadIntPhiInRMonteCarlo[iR]=		histoHadIntRPhiMonteCarlo->ProjectionY( Form("histoMappingHadIntPhiInRMonteCarlo_%i",iR), histoHadIntRPhiMonteCarlo->GetXaxis()->FindBin(arrayRBins[iR]), histoHadIntRPhiMonteCarlo->GetXaxis()->FindBin(arrayRBins[iR+1]));
		ConvGammaRebinWithBinCorrection(histoMappingHadIntPhiInRMonteCarlo[iR],rebin);
		GammaScalingHistogramm(histoMappingHadIntPhiInRMonteCarlo[iR],normFactorReconstMonteCarlo*1/600*1/(TMath::Abs(histoHadIntRPhiMonteCarlo->GetXaxis()->GetBinUpEdge(histoHadIntRPhiMonteCarlo->GetXaxis()->FindBin(arrayRBins[iR+1]))-histoHadIntRPhiMonteCarlo->GetXaxis()->GetBinLowEdge(histoHadIntRPhiMonteCarlo->GetXaxis()->FindBin(arrayRBins[iR])))));
		nameHistoRatioPhiInRMonteCarlo=Form("histoMappingRatioPhiInR_%02d",iR);
		histoMappingRatioPhiInR[iR]= 			(TH1F*)histoMappingHadIntPhiInRData[iR]->Clone();
		histoMappingRatioPhiInR[iR]->SetName(nameHistoRatioPhiInRMonteCarlo); 
		histoMappingRatioPhiInR[iR]->Divide(histoMappingRatioPhiInR[iR],histoMappingHadIntPhiInRMonteCarlo[iR]);
		
	}

	TH2F* 	histoHadIntZPhiMonteCarlo = (TH2F*)directorySecHadIntMonteCarlo->Get(Form("ESD_HadIntMap_ZPhi"));  
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingHadIntPhiInZMonteCarlo[iZ]=		histoHadIntZPhiMonteCarlo->ProjectionY( Form("histoMappingHadIntPhiInZMonteCarlo_%i",iZ), histoHadIntZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ]), histoHadIntZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]));
		ConvGammaRebinWithBinCorrection(histoMappingHadIntPhiInZMonteCarlo[iZ],2*rebin);
		GammaScalingHistogramm(histoMappingHadIntPhiInZMonteCarlo[iZ],normFactorReconstMonteCarlo*1/180*1/(TMath::Abs(histoHadIntZPhiMonteCarlo->GetXaxis()->GetBinUpEdge(histoHadIntZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoHadIntZPhiMonteCarlo->GetXaxis()->GetBinLowEdge(histoHadIntZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ])))));
		nameHistoRatioPhiInZMonteCarlo=Form("histoMappingRatioPhiInZ_%02d",iZ);
		histoMappingRatioPhiInZ[iZ]= 			(TH1F*)histoMappingHadIntPhiInZData[iZ]->Clone();
		histoMappingRatioPhiInZ[iZ]->SetName(nameHistoRatioPhiInZMonteCarlo); 
		histoMappingRatioPhiInZ[iZ]->Divide(histoMappingRatioPhiInZ[iZ],histoMappingHadIntPhiInZMonteCarlo[iZ]);
		
	}

	TH2F* 	histoHadIntZRMonteCarlo = 	(TH2F*) histoHadIntVsZRMonteCarlo->Clone("histoHadIntZRMonteCarlo");
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingHadIntRInZMonteCarlo[iZ]=		histoHadIntZRMonteCarlo->ProjectionY( Form("histoMappingHadIntRInZMonteCarlo_%i",iZ), histoHadIntZRMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ]), histoHadIntZRMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]));
		ConvGammaRebinWithBinCorrection(histoMappingHadIntRInZMonteCarlo[iZ],rebin);
		GammaScalingHistogramm(histoMappingHadIntRInZMonteCarlo[iZ],normFactorReconstMonteCarlo*1/(2*TMath::Pi())*1/(TMath::Abs(histoHadIntZRMonteCarlo->GetXaxis()->GetBinUpEdge(histoHadIntZRMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoHadIntZRMonteCarlo->GetXaxis()->GetBinLowEdge(histoHadIntZRMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ])))));
		nameHistoRatioRInZMonteCarlo=Form("histoMappingRatioRInZ_%02d",iZ);
		histoMappingRatioRInZ[iZ]= 			(TH1F*)histoMappingHadIntRInZData[iZ]->Clone();
		histoMappingRatioRInZ[iZ]->SetName(nameHistoRatioRInZMonteCarlo); 
		histoMappingRatioRInZ[iZ]->Divide(histoMappingRatioRInZ[iZ],histoMappingHadIntRInZMonteCarlo[iZ]);
	}

	for(Int_t iR = 0; iR < nBinsR; iR++){
		histoMappingHadIntZInRMonteCarlo[iR]=		histoHadIntZRMonteCarlo->ProjectionX( Form("histoMappingHadIntZInRMonteCarlo_%i",iR), histoHadIntZRMonteCarlo->GetYaxis()->FindBin(arrayRBins[iR]), histoHadIntZRMonteCarlo->GetYaxis()->FindBin(arrayRBins[iR+1]));
		ConvGammaRebinWithBinCorrection(histoMappingHadIntZInRMonteCarlo[iR],rebin);
		GammaScalingHistogramm(histoMappingHadIntZInRMonteCarlo[iR],normFactorReconstMonteCarlo*1/(2*TMath::Pi())*1/(TMath::Abs(histoHadIntZRMonteCarlo->GetYaxis()->GetBinUpEdge(histoHadIntZRMonteCarlo->GetYaxis()->FindBin(arrayRBins[iR+1]))-histoHadIntZRMonteCarlo->GetYaxis()->GetBinLowEdge(histoHadIntZRMonteCarlo->GetYaxis()->FindBin(arrayRBins[iR])))));
		nameHistoRatioZInRMonteCarlo=Form("histoMappingRatioZInR_%02d",iR);
		histoMappingRatioZInR[iR]= (TH1F*)histoMappingHadIntZInRData[iR]->Clone();
		histoMappingRatioZInR[iR]->SetName(nameHistoRatioZInRMonteCarlo); 
		histoMappingRatioZInR[iR]->Divide(histoMappingRatioZInR[iR],histoMappingHadIntZInRMonteCarlo[iR]);
	}

	TH2F* 	histoHadIntSPDZPhiMonteCarlo = (TH2F*)directorySecHadIntMonteCarlo->Get(Form("ESD_HadIntMap_SPD_ZPhi"));  
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingHadIntSPDPhiInZMonteCarlo[iZ]=		histoHadIntSPDZPhiMonteCarlo->ProjectionY( Form("histoMappingHadIntSPDPhiInZMonteCarlo_%i",iZ), histoHadIntSPDZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ]), histoHadIntSPDZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]));
		ConvGammaRebinWithBinCorrection(histoMappingHadIntSPDPhiInZMonteCarlo[iZ],2*rebin);
		GammaScalingHistogramm(histoMappingHadIntSPDPhiInZMonteCarlo[iZ],normFactorReconstMonteCarlo*1/25*1/(TMath::Abs(histoHadIntSPDZPhiMonteCarlo->GetXaxis()->GetBinUpEdge(histoHadIntSPDZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoHadIntSPDZPhiMonteCarlo->GetXaxis()->GetBinLowEdge(histoHadIntSPDZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ])))));
		nameHistoRatioSPDPhiInZMonteCarlo=Form("histoMappingRatioSPDPhiInZ_%02d",iZ);
		histoMappingRatioSPDPhiInZ[iZ]= 		(TH1F*)histoMappingHadIntSPDPhiInZData[iZ]->Clone();
		histoMappingRatioSPDPhiInZ[iZ]->SetName(nameHistoRatioSPDPhiInZMonteCarlo); 
		histoMappingRatioSPDPhiInZ[iZ]->Divide(histoMappingRatioSPDPhiInZ[iZ],histoMappingHadIntSPDPhiInZMonteCarlo[iZ]);
		
	}
	
	TH2F* 	histoHadIntSDDZPhiMonteCarlo = (TH2F*)directorySecHadIntMonteCarlo->Get(Form("ESD_HadIntMap_SDD_ZPhi"));  
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingHadIntSDDPhiInZMonteCarlo[iZ]=		histoHadIntSDDZPhiMonteCarlo->ProjectionY( Form("histoMappingHadIntSDDPhiInZMonteCarlo_%i",iZ), histoHadIntSDDZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ]), histoHadIntSDDZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]));
		ConvGammaRebinWithBinCorrection(histoMappingHadIntSDDPhiInZMonteCarlo[iZ],2*rebin);
		GammaScalingHistogramm(histoMappingHadIntSDDPhiInZMonteCarlo[iZ],normFactorReconstMonteCarlo*1/25*1/(TMath::Abs(histoHadIntSDDZPhiMonteCarlo->GetXaxis()->GetBinUpEdge(histoHadIntSDDZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoHadIntSDDZPhiMonteCarlo->GetXaxis()->GetBinLowEdge(histoHadIntSDDZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ])))));
		nameHistoRatioSDDPhiInZMonteCarlo=Form("histoMappingRatioSDDPhiInZ_%02d",iZ);
		histoMappingRatioSDDPhiInZ[iZ]= 		(TH1F*)histMappingHadIntSDDPhiInZData[iZ]->Clone();
		histoMappingRatioSDDPhiInZ[iZ]->SetName(nameHistoRatioSDDPhiInZMonteCarlo); 
		histoMappingRatioSDDPhiInZ[iZ]->Divide(histoMappingRatioSDDPhiInZ[iZ],histoMappingHadIntSDDPhiInZMonteCarlo[iZ]);
		
	}
	
	TH2F* 	histoHadIntHotZoneZPhiMonteCarlo = (TH2F*)directorySecHadIntMonteCarlo->Get(Form("ESD_HadIntMap_HotZone_ZPhi"));  
	for(Int_t iZ = 0; iZ < nBinsZHotZone; iZ++){
		histoMappingHadIntHotZonePhiInZMonteCarlo[iZ]=		histoHadIntHotZoneZPhiMonteCarlo->ProjectionY( Form("histoMappingHadIntHotZonePhiInZMonteCarlo_%i",iZ), histoHadIntHotZoneZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBinsHotZone[iZ]), histoHadIntHotZoneZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBinsHotZone[iZ+1]));
		ConvGammaRebinWithBinCorrection(histoMappingHadIntHotZonePhiInZMonteCarlo[iZ],2*rebin);
		GammaScalingHistogramm(histoMappingHadIntHotZonePhiInZMonteCarlo[iZ],normFactorReconstMonteCarlo*1/30*1/(TMath::Abs(histoHadIntHotZoneZPhiMonteCarlo->GetXaxis()->GetBinUpEdge(histoHadIntHotZoneZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBinsHotZone[iZ+1]))-histoHadIntHotZoneZPhiMonteCarlo->GetXaxis()->GetBinLowEdge(histoHadIntHotZoneZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBinsHotZone[iZ])))));
		nameHistoRatioHotZonePhiInZMonteCarlo=Form("histoMappingRatioHotZonePhiInZ_%02d",iZ);
		histoMappingRatioHotZonePhiInZ[iZ]= 		(TH1F*)histoMappingHadIntHotZonePhiInZData[iZ]->Clone();
		histoMappingRatioHotZonePhiInZ[iZ]->SetName(nameHistoRatioHotZonePhiInZMonteCarlo); 
		histoMappingRatioHotZonePhiInZ[iZ]->Divide(histoMappingRatioHotZonePhiInZ[iZ],histoMappingHadIntHotZonePhiInZMonteCarlo[iZ]);
	}			       

	
	TH2F* 	histoHadIntITSTPCZPhiMonteCarlo = (TH2F*)directorySecHadIntMonteCarlo->Get(Form("ESD_HadIntMap_ITSTPC_ZPhi"));  
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingHadIntITSTPCPhiInZMonteCarlo[iZ]=		histoHadIntITSTPCZPhiMonteCarlo->ProjectionY( Form("histoMappingHadIntITSTPCPhiInZMonteCarlo_%i",iZ), histoHadIntITSTPCZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ]), histoHadIntITSTPCZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]));
		ConvGammaRebinWithBinCorrection(histoMappingHadIntITSTPCPhiInZMonteCarlo[iZ],2*rebin);
		GammaScalingHistogramm(histoMappingHadIntITSTPCPhiInZMonteCarlo[iZ],normFactorReconstMonteCarlo*1/30*1/(TMath::Abs(histoHadIntITSTPCZPhiMonteCarlo->GetXaxis()->GetBinUpEdge(histoHadIntITSTPCZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoHadIntITSTPCZPhiMonteCarlo->GetXaxis()->GetBinLowEdge(histoHadIntITSTPCZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ])))));
		nameHistoRatioITSTPCPhiInZMonteCarlo=Form("histoMappingRatioITSTPCPhiInZ_%02d",iZ);
		histoMappingRatioITSTPCPhiInZ[iZ]= 		(TH1F*)histoMappingHadIntITSTPCPhiInZData[iZ]->Clone();
		histoMappingRatioITSTPCPhiInZ[iZ]->SetName(nameHistoRatioITSTPCPhiInZMonteCarlo); 
		histoMappingRatioITSTPCPhiInZ[iZ]->Divide(histoMappingRatioITSTPCPhiInZ[iZ],histoMappingHadIntITSTPCPhiInZMonteCarlo[iZ]);
	}			       
	
	Float_t minimumPhiInRBins = 0.;
	Float_t minimumZInRBins = 0.1*normFactorReconstData;
	Float_t minimumRInZBins = 0.1*normFactorReconstData;
	Float_t minimumPhiInZBins = 0;
	Float_t minimumPhiInZITSTPCBins =0.;
	Float_t minimumPhiInZSPDBins = 0.;
	Float_t minimumPhiInZSDDBins = 0.;
	Float_t minimumPhiInZHotZoneBins = 0.;
	
	//------------------- Giving the Pad Phi in R in better resolution
	
	TCanvas* canvasSinglePhiInR = new TCanvas("canvasSinglePhiInR","",400,20,2000,1400);  // gives the page size
	TPad* padSinglePhiInR = new TPad("padSinglePhiInR","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
	padSinglePhiInR->SetFillColor(0);
	padSinglePhiInR->GetFrame()->SetFillColor(0);
	padSinglePhiInR->SetBorderMode(0);
	padSinglePhiInR->Divide(rowR-3,columnR+2);
	padSinglePhiInR->Draw();
	
	for(Int_t iR = iRStart; iR < nBinsR; iR++){
		Int_t place = iR + 1;
		
		padSinglePhiInR->cd(place);
		padSinglePhiInR->cd(place)->SetRightMargin(0.01);
		padSinglePhiInR->cd(place)->SetLogy(1);
		TString nameHistoPhiInRData=Form("#Phi in R:  %s",arrayNamesRBins[iR].Data());
		
		DrawAutoGammaHistosMaterial( histoMappingHadIntPhiInRData[iR], 
										histoMappingHadIntPhiInRMonteCarlo[iR], 
										nameHistoPhiInRData,"#Phi",textYAxisSlicesHisto,
										kTRUE,1.3 ,minimumPhiInRBins,
										kFALSE,0. ,0.,
										kFALSE, 0.,180.);
		if (iR == iRStart){DrawAliceText(floatLocationRightUpText2[0],floatLocationRightUpText2[1], floatLocationRightUpText2[3]);}
	}
	padSinglePhiInR->Update();
	canvasSinglePhiInR->Print(Form("%s/pad_Phi_in_R.%s",outputDirectory.Data(),suffix.Data()));   
	delete  	padSinglePhiInR;
	delete 	canvasSinglePhiInR;               
	
	
	// --------------- Giving Phi in R for several bins -------------------------------------------------------------
	for(Int_t iR = iRStart; iR < nBinsR; iR++){
		if (histoMappingHadIntPhiInRData[iR]->GetEntries()>0){
			TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
// 			canvasSingleBin->SetLogy(1);
			canvasSingleBin->SetTopMargin(0.04);
			canvasSingleBin->cd();
			TString nameHistoPhiInRData=Form("#phi in R:  %s",arrayNamesRBins[iR].Data());
			TLatex *latexBinning = new TLatex(0.16,0.76,nameHistoPhiInRData); // Bo: this was modified
			latexBinning->SetNDC();
			latexBinning->SetTextColor(1);
			latexBinning->SetTextFont(62);
			latexBinning->SetTextSize(0.04);
			latexBinning->SetLineWidth(2);       
			histoMappingHadIntPhiInRData[iR]->GetYaxis()->SetMoreLogLabels(kTRUE);
			DrawAutoGammaHistosMaterial( histoMappingHadIntPhiInRData[iR], 
											histoMappingHadIntPhiInRMonteCarlo[iR],
											"","#phi",textYAxisSlicesHisto,
											kFALSE,1.5 ,1e-6,
											kFALSE,0. ,0.,
											kFALSE, 0.,180.);
			latexBinning->Draw("same");  
// 			DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
			
			canvasSingleBin->Update();
			canvasSingleBin->SaveAs(Form("%s/Phi_in_R_%i.%s",outputDirectory.Data(),iR,suffix.Data()));
			delete canvasSingleBin;
		}
	}               
	
	//----------- Giving Z in R as single Pad ----------------
	TCanvas* canvasSingleZInR = new TCanvas("canvasSingleZInR","",400,20,2000,1400);  // gives the page size
	
	TPad* padSingleZInR = new TPad("padSingleZInR","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
	padSingleZInR->SetFillColor(0);
	padSingleZInR->GetFrame()->SetFillColor(0);
	padSingleZInR->SetBorderMode(0);
	padSingleZInR->Divide(rowR-3,columnR+2);
	padSingleZInR->Draw();
	
	for(Int_t iR = iRStart; iR < nBinsR; iR++){
		Float_t rangeZ = arrayRangeZInR[iR];
		Int_t place = iR + 1;
		padSingleZInR->cd(place);
		padSingleZInR->cd(place)->SetRightMargin(0.01);
		padSingleZInR->cd(place)->SetLogy(1);
		TString nameHistoZInRData=Form("Z in R:  %s",arrayNamesRBins[iR].Data());
		DrawAutoGammaHistosMaterial( histoMappingHadIntZInRData[iR], 
										histoMappingHadIntZInRMonteCarlo[iR], 
										nameHistoZInRData,"Z (cm)",textYAxisSlicesHisto,
										kTRUE,2 ,minimumZInRBins,
										kFALSE,0. ,0.,
										kTRUE, -rangeZ,rangeZ);
// 		if (iR == iRStart){DrawAliceText(floatLocationRightUpText2[0],floatLocationRightUpText2[1], floatLocationRightUpText2[3]);}
	}
	
	padSingleZInR->Update();
	canvasSingleZInR->SaveAs(Form("%s/pad_Z_in_R.%s",outputDirectory.Data(),suffix.Data()));
	delete padSingleZInR;
	delete canvasSingleZInR;
	
	// -------------- Giving Z in R for all bins ----------------------------------------------
	
	for(Int_t iR = iRStart; iR < nBinsR; iR++){
		if (histoMappingHadIntZInRData[iR]->GetEntries()>0){
			TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
			canvasSingleBin->SetLogy(1);
			canvasSingleBin->cd();
			Float_t rangeZ = arrayRangeZInR[iR];
			TString nameHistoZInRData=Form("Z in R:  %s",arrayNamesRBins[iR].Data());
			TLatex *latexBinning = new TLatex(0.16,0.76,nameHistoZInRData); // Bo: this was modified
			latexBinning->SetNDC();
			latexBinning->SetTextColor(1);
			latexBinning->SetTextFont(62);
			latexBinning->SetTextSize(0.04);
			latexBinning->SetLineWidth(2);       
// 			histoMappingHadIntZInRData[iR]->GetYaxis()->SetMoreLogLabels(kTRUE);
			DrawAutoGammaHistosMaterial( histoMappingHadIntZInRData[iR], 
											histoMappingHadIntZInRMonteCarlo[iR], 
											"","Z",textYAxisSlicesHisto,
											kFALSE,2 ,minimumZInRBins,
											kFALSE,0. ,0.,
											kTRUE, -rangeZ,rangeZ);
			latexBinning->Draw("same");  
// 			DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
			
			canvasSingleBin->Update();
			canvasSingleBin->SaveAs(Form("%s/Z_in_R_%i.%s",outputDirectory.Data(),iR,suffix.Data()));
			delete canvasSingleBin;
		}
	}               
	
	
	// ------------- Giving phi in Z in singleplot
	TCanvas* canvasSinglePhiInZ = new TCanvas("canvasSinglePhiInZ","",200,10,2000,1400);  // gives the page size
	
	TPad* padSinglePhiInZ = new TPad("padSinglePhiInZ","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
	padSinglePhiInZ->SetFillColor(0);
	padSinglePhiInZ->GetFrame()->SetFillColor(0);
	padSinglePhiInZ->SetBorderMode(0);
	padSinglePhiInZ->Divide(nRowsZProjections-1 ,nColumnsZProjections+1);
	padSinglePhiInZ->Draw();
	
	
	for(Int_t iZ = minimumBinZNormal ; iZ < (nBinsZ-minimumBinZNormal); iZ++){
		Int_t place = iZ -minimumBinZNormal+ 1;
		padSinglePhiInZ->cd(place);
// 		padSinglePhiInZ->cd(place)->SetLogy(1);
		padSinglePhiInZ->cd(place)->SetRightMargin(0.01);
		TString nameHistoPhiInZData=Form("#phi in Z:  %s",arrayNamesZBins[iZ].Data());
		DrawAutoGammaHistosMaterial( histoMappingHadIntPhiInZData[iZ], 
										histoMappingHadIntPhiInZMonteCarlo[iZ], 
										nameHistoPhiInZData,"#Phi",textYAxisSlicesHisto,
										kTRUE,1.3 ,minimumPhiInZBins,
										kFALSE,0. ,0.,
										kFALSE, 0.,180.);
// 		if (iZ == 0){DrawAliceText(floatLocationRightUpText2[0],floatLocationRightUpText2[1], floatLocationRightUpText2[3]);}
	}
	
	padSinglePhiInZ->Update();
	canvasSinglePhiInZ->SaveAs(Form("%s/pad_Phi_in_Z.%s",outputDirectory.Data(),suffix.Data()));
	delete padSinglePhiInZ ;
	delete canvasSinglePhiInZ;
	
	// ************************************* Phi in Z for all bins *****************************************
	for(Int_t iZ = minimumBinZNormal ; iZ < (nBinsZ-minimumBinZNormal); iZ++){
		TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
// 		canvasSingleBin->SetLogy(1);
		canvasSingleBin->SetTopMargin(0.04);            
		canvasSingleBin->cd();
		TString nameHistoPhiInZData=Form("#phi in Z:  %s",arrayNamesZBins[iZ].Data());
		TLatex *latexBinning = new TLatex(0.16,0.76,nameHistoPhiInZData); // Bo: this was modified
		latexBinning->SetNDC();
		latexBinning->SetTextColor(1);
		latexBinning->SetTextFont(62);
		latexBinning->SetTextSize(0.04);
		latexBinning->SetLineWidth(2);       
		histoMappingHadIntPhiInZData[iZ]->GetYaxis()->SetMoreLogLabels(kTRUE);
		DrawAutoGammaHistosMaterial(histoMappingHadIntPhiInZData[iZ], 
										histoMappingHadIntPhiInZMonteCarlo[iZ],
										"","#phi",textYAxisSlicesHisto,
										kFALSE,2 ,1e-9,
										kFALSE,0. ,0.,
										kFALSE, -0,0);
		latexBinning->Draw("same");  
// 		DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
		
		canvasSingleBin->Update();
		canvasSingleBin->SaveAs(Form("%s/Phi_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
		delete canvasSingleBin;
	}               
	
	// ------------- Giving SPD Z in phi in singleplot
	TCanvas* canvasSPDPhiInZSP = new TCanvas("canvasSPDPhiInZSP","",200,10,2000,1400);  // gives the page size
	
	TPad* padSPDPhiInZSP = new TPad("padSPDPhiInZSP","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
	padSPDPhiInZSP->SetFillColor(0);
	padSPDPhiInZSP->GetFrame()->SetFillColor(0);
	padSPDPhiInZSP->SetBorderMode(0);
	padSPDPhiInZSP->Divide(nRowsZSPD-1 ,nColumnsZSPD+1);
	padSPDPhiInZSP->Draw();
	
	
	for(Int_t iZ = minimumBinZSPD ; iZ < (nBinsZ-minimumBinZSPD); iZ++){
		Int_t place = iZ -minimumBinZSPD+ 1;
		padSPDPhiInZSP->cd(place);
		padSPDPhiInZSP->cd(place)->SetLogy(1);
		padSPDPhiInZSP->cd(place)->SetRightMargin(0.01);
		nameHistoSPDPhiInZData=Form("0 cm < R < 13 cm: #phi in Z:  %s",arrayNamesZBins[iZ].Data());
		DrawAutoGammaHistosMaterial( histoMappingHadIntSPDPhiInZData[iZ], 
										histoMappingHadIntSPDPhiInZMonteCarlo[iZ], 
										nameHistoSPDPhiInZData,"#Phi",textYAxisSlicesHisto,
										kTRUE,1.3 ,minimumPhiInZSPDBins,
										kFALSE,0. ,0.,
										kFALSE, 0.,180.);
// 		if (iZ == 0){DrawAliceText(floatLocationRightUpText2[0],floatLocationRightUpText2[1], floatLocationRightUpText2[3]);}
	}
	
	padSPDPhiInZSP->Update();
	canvasSPDPhiInZSP->SaveAs(Form("%s/pad_SPD_Phi_in_Z.%s",outputDirectory.Data(),suffix.Data()));
	delete padSPDPhiInZSP ;
	delete canvasSPDPhiInZSP;
	
	// ************************************* Phi in Z for all bins *****************************************
	for(Int_t iZ = minimumBinZSPD ; iZ < (nBinsZ-minimumBinZSPD); iZ++){
		TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
// 		canvasSingleBin->SetLogy(1);
		canvasSingleBin->SetTopMargin(0.04);
		canvasSingleBin->cd();
		nameHistoSPDPhiInZData=Form("0 cm < R < 13 cm: #phi in Z:  %s",arrayNamesZBins[iZ].Data());
		TLatex *latexBinning = new TLatex(0.16,0.76,nameHistoSPDPhiInZData); // Bo: this was modified
		latexBinning->SetNDC();
		latexBinning->SetTextColor(1);
		latexBinning->SetTextFont(62);
		latexBinning->SetTextSize(0.04);
		latexBinning->SetLineWidth(2);       
		histoMappingHadIntSPDPhiInZData[iZ]->GetYaxis()->SetMoreLogLabels(kTRUE);
		DrawAutoGammaHistosMaterial(histoMappingHadIntSPDPhiInZData[iZ], 
										histoMappingHadIntSPDPhiInZMonteCarlo[iZ],
										"","#phi",textYAxisSlicesHisto,
										kFALSE,2 ,0.0000001,
										kFALSE,0. ,0.,
										kFALSE, -0,0);
		latexBinning->Draw("same");  
// 		DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
		
		canvasSingleBin->Update();
		canvasSingleBin->SaveAs(Form("%s/SPD_Phi_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
		delete canvasSingleBin;
	}               

	// ------------- Giving SPD Z in phi in singleplot
	TCanvas* canvasSDDPhiInZSP = new TCanvas("canvasSDDPhiInZSP","",200,10,2000,1400);  // gives the page size
	
	TPad* padSDDPhiInZSP = new TPad("padSDDPhiInZSP","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
	padSDDPhiInZSP->SetFillColor(0);
	padSDDPhiInZSP->GetFrame()->SetFillColor(0);
	padSDDPhiInZSP->SetBorderMode(0);
	padSDDPhiInZSP->Divide(nRowsZSDD-1 ,nColumnsZSDD+1);
	padSDDPhiInZSP->Draw();
		
	for(Int_t iZ = minimumBinZSDD ; iZ < (nBinsZ-minimumBinZSDD); iZ++){
		Int_t place = iZ -minimumBinZSDD+ 1;
		cout << "iZ " << iZ << "\t place" << place << endl;
		padSDDPhiInZSP->cd(place);
		padSDDPhiInZSP->cd(place)->SetLogy(1);
		padSDDPhiInZSP->cd(place)->SetRightMargin(0.01);
		nameHistoSDDPhiInZData=Form("13 cm < R < 25 cm: #phi in Z:  %s",arrayNamesZBins[iZ].Data());
		DrawAutoGammaHistosMaterial( histMappingHadIntSDDPhiInZData[iZ], 
										histoMappingHadIntSDDPhiInZMonteCarlo[iZ], 
										nameHistoSDDPhiInZData,"#Phi",textYAxisSlicesHisto,
										kTRUE,1.3 ,minimumPhiInZSDDBins,
										kFALSE,0. ,0.,
										kFALSE, 0.,180.);
// 		if (iZ == 0){DrawAliceText(floatLocationRightUpText2[0],floatLocationRightUpText2[1], floatLocationRightUpText2[3]);}
	}
	
	padSDDPhiInZSP->Update();
	canvasSDDPhiInZSP->SaveAs(Form("%s/pad_SDD_Phi_in_Z.%s",outputDirectory.Data(),suffix.Data()));
	delete padSDDPhiInZSP ;
	delete canvasSDDPhiInZSP;
	
	// ************************************* Phi in Z for all bins *****************************************
	for(Int_t iZ = minimumBinZSDD ; iZ < (nBinsZ-minimumBinZSDD); iZ++){
		TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
// 		canvasSingleBin->SetLogy(1);
		canvasSingleBin->SetTopMargin(0.04);
		canvasSingleBin->cd();
		nameHistoSDDPhiInZData=Form("13 cm < R < 25 cm: #phi in Z:  %s",arrayNamesZBins[iZ].Data());
		TLatex *latexBinning = new TLatex(0.16,0.76,nameHistoSDDPhiInZData); // Bo: this was modified
		latexBinning->SetNDC();
		latexBinning->SetTextColor(1);
		latexBinning->SetTextFont(62);
		latexBinning->SetTextSize(0.04);
		latexBinning->SetLineWidth(2);       
		histMappingHadIntSDDPhiInZData[iZ]->GetYaxis()->SetMoreLogLabels(kTRUE);	
		DrawAutoGammaHistosMaterial(histMappingHadIntSDDPhiInZData[iZ], 
										histoMappingHadIntSDDPhiInZMonteCarlo[iZ],
										"","#phi",textYAxisSlicesHisto,
										kFALSE,2 ,0.0000001,
										kFALSE,0. ,0.,
										kFALSE, -0,0);
		latexBinning->Draw("same");  
// 		DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
		
		canvasSingleBin->Update();
		canvasSingleBin->SaveAs(Form("%s/SDD_Phi_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
		delete canvasSingleBin;
	}               

	// ------------- Giving ITSTPC Z in phi in singleplot
	TCanvas* canvasITSTPCPhiInZSP = new TCanvas("canvasITSTPCPhiInZSP","",200,10,2000,1400);  // gives the page size
	
	TPad* padITSTPCPhiInZSP = new TPad("padITSTPCPhiInZSP","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
	padITSTPCPhiInZSP->SetFillColor(0);
	padITSTPCPhiInZSP->GetFrame()->SetFillColor(0);
	padITSTPCPhiInZSP->SetBorderMode(0);
	padITSTPCPhiInZSP->Divide(nRowsZITSTPC-2 ,nColumnsZITSTPC+2);
	padITSTPCPhiInZSP->Draw();
	
	
	for(Int_t iZ = minimumBinZITSTPC ; iZ < (nBinsZ-minimumBinZITSTPC); iZ++){
		Int_t place = iZ -minimumBinZITSTPC+ 1;
		padITSTPCPhiInZSP->cd(place);
		padITSTPCPhiInZSP->cd(place)->SetLogy(1);
		padITSTPCPhiInZSP->cd(place)->SetRightMargin(0.01);
		TString nameHistoITSTPCPhiInZData=Form("50 cm < R < 80 cm: #phi in Z:  %s",arrayNamesZBins[iZ].Data());
		DrawAutoGammaHistosMaterial( histoMappingHadIntITSTPCPhiInZData[iZ], 
						histoMappingHadIntITSTPCPhiInZMonteCarlo[iZ], 
						nameHistoITSTPCPhiInZData,"#Phi",textYAxisSlicesHisto,
						kTRUE,1.3 ,minimumPhiInZITSTPCBins,
						kFALSE,0. ,0.,
						kFALSE, 0.,180.);
// 		if (iZ == 0){DrawAliceText(floatLocationRightUpText2[0],floatLocationRightUpText2[1], floatLocationRightUpText2[3]);}
	}
	
	padITSTPCPhiInZSP->Update();
	canvasITSTPCPhiInZSP->SaveAs(Form("%s/pad_ITSTPC_Phi_in_Z.%s",outputDirectory.Data(),suffix.Data()));
	delete padITSTPCPhiInZSP ;
	delete canvasITSTPCPhiInZSP;
	
	// ************************************* ITSTPC Phi in Z for all bins *****************************************
	for(Int_t iZ = minimumBinZITSTPC ; iZ < (nBinsZ-minimumBinZITSTPC); iZ++){
		TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
// 		canvasSingleBin->SetLogy(1);
		canvasSingleBin->SetTopMargin(0.04);
		canvasSingleBin->cd();
		TString nameHistoITSTPCPhiInZData=Form("50 cm < R < 80 cm: #phi in Z:  %s",arrayNamesZBins[iZ].Data());
		TLatex *latexBinning = new TLatex(0.16,0.76,nameHistoITSTPCPhiInZData); // Bo: this was modified
		latexBinning->SetNDC();
		latexBinning->SetTextColor(1);
		latexBinning->SetTextFont(62);
		latexBinning->SetTextSize(0.04);
		latexBinning->SetLineWidth(2);       
		histoMappingHadIntITSTPCPhiInZData[iZ]->GetYaxis()->SetMoreLogLabels(kTRUE);	
		DrawAutoGammaHistosMaterial(histoMappingHadIntITSTPCPhiInZData[iZ], 
						histoMappingHadIntITSTPCPhiInZMonteCarlo[iZ],
						"","#phi",textYAxisSlicesHisto,
						kFALSE,1 ,1e-7,
						kFALSE,0. ,0.,
						kFALSE, -0,0);
		latexBinning->Draw("same");  
// 		DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
		
		canvasSingleBin->Update();
		canvasSingleBin->SaveAs(Form("%s/ITSTPC_Phi_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
		delete canvasSingleBin;
	}               
	
	
	// ------------- Giving HotZone Z in phi in singleplot
	TCanvas* canvasHotZonePhiInZSP = new TCanvas("canvasHotZonePhiInZSP","",200,10,2000,1400);  // gives the page size
	
	TPad* padHotZonePhiInZSP = new TPad("padHotZonePhiInZSP","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
	padHotZonePhiInZSP->SetFillColor(0);
	padHotZonePhiInZSP->GetFrame()->SetFillColor(0);
	padHotZonePhiInZSP->SetBorderMode(0);
	padHotZonePhiInZSP->Divide(nRowsZHotZone-2 ,nColumnsZHotZone+2);
	padHotZonePhiInZSP->Draw();
	
	
	for(Int_t iZ = minimumBinZHotZone ; iZ < (nBinsZHotZone-minimumBinZHotZone); iZ++){
		Int_t place = iZ -minimumBinZHotZone+ 1;
		padHotZonePhiInZSP->cd(place);
		padHotZonePhiInZSP->cd(place)->SetLogy(1);
		padHotZonePhiInZSP->cd(place)->SetRightMargin(0.01);
		TString nameHistoHotZonePhiInZData=Form("5.7 cm < R < 50 cm: #phi in Z:  %s",arrayNamesZBinsHotZone[iZ].Data());
		DrawAutoGammaHistosMaterial( histoMappingHadIntHotZonePhiInZData[iZ], 
						histoMappingHadIntHotZonePhiInZMonteCarlo[iZ], 
						nameHistoHotZonePhiInZData,"#Phi",textYAxisSlicesHisto,
						kTRUE,1.3 ,minimumPhiInZHotZoneBins,
						kFALSE,0. ,0.,
						kFALSE, 0.,180.);
// 		if (iZ == 0){DrawAliceText(floatLocationRightUpText2[0],floatLocationRightUpText2[1], floatLocationRightUpText2[3]);}
	}
	
	padHotZonePhiInZSP->Update();
	canvasHotZonePhiInZSP->SaveAs(Form("%s/pad_HotZone_Phi_in_Z.%s",outputDirectory.Data(),suffix.Data()));
	delete padHotZonePhiInZSP ;
	delete canvasHotZonePhiInZSP;
	
	// ************************************* HotZone Phi in Z for all bins *****************************************
	for(Int_t iZ = minimumBinZHotZone ; iZ < (nBinsZHotZone-minimumBinZHotZone); iZ++){
		TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
// 		canvasSingleBin->SetLogy(1);
		canvasSingleBin->SetTopMargin(0.04);
		canvasSingleBin->cd();
		TString nameHistoHotZonePhiInZData=Form("5.7 cm < R < 50 cm: #phi in Z:  %s",arrayNamesZBinsHotZone[iZ].Data());
		TLatex *latexBinning = new TLatex(0.16,0.76,nameHistoHotZonePhiInZData); // Bo: this was modified
		latexBinning->SetNDC();
		latexBinning->SetTextColor(1);
		latexBinning->SetTextFont(62);
		latexBinning->SetTextSize(0.04);
		latexBinning->SetLineWidth(2);       
		histoMappingHadIntHotZonePhiInZData[iZ]->GetYaxis()->SetMoreLogLabels(kTRUE);	
		DrawAutoGammaHistosMaterial(histoMappingHadIntHotZonePhiInZData[iZ], 
						histoMappingHadIntHotZonePhiInZMonteCarlo[iZ],
						"","#phi",textYAxisSlicesHisto,
						kFALSE,1 ,1e-7,
						kFALSE,0. ,0.,
						kFALSE, -0,0);
		latexBinning->Draw("same");  
// 		DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
		
		canvasSingleBin->Update();
		canvasSingleBin->SaveAs(Form("%s/HotZone_Phi_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
		delete canvasSingleBin;
	}               
	

	// ----- Giving R in Z in SinglePlot    
	TCanvas* canvasSingleRInZ = new TCanvas("canvasSingleRInZ","",400,20,2000,1400);  // gives the page size
	
	TPad* padSingleRInZ = new TPad("padSingleRInZ","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
	padSingleRInZ->SetFillColor(0);
	padSingleRInZ->GetFrame()->SetFillColor(0);
	padSingleRInZ->SetBorderMode(0);
	padSingleRInZ->Divide(nRowsZProjections-1 ,nColumnsZProjections+1 );
	padSingleRInZ->Draw();
	
	
	for(Int_t iZ = minimumBinZNormal ; iZ < (nBinsZ-minimumBinZNormal); iZ++){
		Int_t place = iZ -minimumBinZNormal+ 1;
		padSingleRInZ->cd(place);
		padSingleRInZ->cd(place)->SetLogy(1);
		padSingleRInZ->cd(place)->SetRightMargin(0.01);
		TString nameHistoRInZData=Form("R in Z:  %s",arrayNamesZBins[iZ].Data());
		DrawAutoGammaHistosMaterial( histoMappingHadIntRInZData[iZ], 
										histoMappingHadIntRInZMonteCarlo[iZ], 
										nameHistoRInZData,"R (cm)",textYAxisSlicesHisto,
										kTRUE,2. ,minimumRInZBins,
										kFALSE,0. ,0.,
										kFALSE, 0.,180.);
// 		if (iZ == 0){DrawAliceText(floatLocationRightUpText2[0],floatLocationRightUpText2[1], floatLocationRightUpText2[3]);}
	}
	
	padSingleRInZ->Update();
	canvasSingleRInZ->SaveAs(Form("%s/pad_R_in_Z.%s",outputDirectory.Data(),suffix.Data()));
	delete padSingleRInZ ;
	delete canvasSingleRInZ;
	
	for(Int_t iZ = minimumBinZNormal ; iZ < (nBinsZ-minimumBinZNormal); iZ++){
		TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
		canvasSingleBin->SetLogy(1);
		canvasSingleBin->cd();
		TString nameHistoRInZData=Form("R in Z:  %s",arrayNamesZBins[iZ].Data());
		TLatex *latexBinning = new TLatex(0.16,0.76,nameHistoRInZData); // Bo: this was modified
		latexBinning->SetNDC();
		latexBinning->SetTextColor(1);
		latexBinning->SetTextFont(62);
		latexBinning->SetTextSize(0.04);
		latexBinning->SetLineWidth(2);       
		
		DrawAutoGammaHistosMaterial(histoMappingHadIntRInZData[iZ], 
										histoMappingHadIntRInZMonteCarlo[iZ],
										"","R",textYAxisSlicesHisto,
										kFALSE,2 ,0.5*normFactorReconstData,
										kFALSE,0. ,0.,
										kTRUE, 0,180);
		latexBinning->Draw("same");  
// 		DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
		
		canvasSingleBin->Update();
		canvasSingleBin->SaveAs(Form("%s/R_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
		delete canvasSingleBin;
	}               
	
	
	
	//---------------------------- Ratio Pads ----------------------------------
			
	TCanvas* canvasRatioPhInRSP = new TCanvas("canvasRatioPhInRSP","",400,20,2000,1400);  // gives the page size
	
	TPad* padRatioPhInRSP = new TPad("padRatioPhInRSP","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
	padRatioPhInRSP->SetFillColor(0);
	padRatioPhInRSP->GetFrame()->SetFillColor(0);
	padRatioPhInRSP->SetBorderMode(0);
	padRatioPhInRSP->Divide(rowR-3,columnR+2);
	padRatioPhInRSP->Draw();
	
	for(Int_t iR = iRStart; iR < nBinsR; iR++){
			Int_t place = iR + 1;
			
			padRatioPhInRSP->cd(place);
			if (boolRatioLog) padRatioPhInRSP->cd(place)->SetLogy(1);
			padRatioPhInRSP->cd(place)->SetRightMargin(0.01);
			nameHistoRatioPhiInRMonteCarlo=Form("#phi in R:  %s",arrayNamesRBins[iR].Data());
			DrawRatioGammaHisto(    histoMappingRatioPhiInR[iR], 
											nameHistoRatioPhiInRMonteCarlo,"#Phi","norm Data/norm MC",
											kFALSE,3 ,0.000001,
											kTRUE, minYValueRatio , maxYValueRatio,
											kFALSE, 0.,180.);
			linePhi->Draw("same");
// 			if (iR == iRStart){DrawAliceText(floatLocationRightUpText2[0],floatLocationRightUpText2[1], floatLocationRightUpText2[3]);}
	}
	
	padRatioPhInRSP->Update();
	canvasRatioPhInRSP->SaveAs(Form("%s/Ratio_Phi_in_R.%s",outputDirectory.Data(),suffix.Data()));
	delete padRatioPhInRSP ;
	delete  canvasRatioPhInRSP;
	
	for(Int_t iR = iRStart; iR < nBinsR; iR++){
		if (histoMappingRatioPhiInR[iR]->GetEntries()>0){
			TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
			if (boolRatioLog) canvasSingleBin->SetLogy(1);
			canvasSingleBin->cd();
			TString nameHistoPhiInRData=Form("#phi in R:  %s",arrayNamesRBins[iR].Data());
			TLatex *latexBinning = new TLatex(0.16,0.76,nameHistoPhiInRData); // Bo: this was modified
			latexBinning->SetNDC();
			latexBinning->SetTextColor(1);
			latexBinning->SetTextFont(62);
			latexBinning->SetTextSize(0.04);
			latexBinning->SetLineWidth(2);       
			
			DrawRatioGammaHisto(    histoMappingRatioPhiInR[iR], 
											"","#Phi","norm Data/norm MC",
											kFALSE,3 ,0.000001,
											kTRUE, minYValueRatio , maxYValueRatio,
											kFALSE, 0.,180.);
			DrawGammaLines(-3.2,3.2,1,1,0.02);              
			
			latexBinning->Draw("same");  
// 			DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
			
			canvasSingleBin->Update();
			canvasSingleBin->SaveAs(Form("%s/Ratio_Phi_in_R_%i.%s",outputDirectory.Data(),iR,suffix.Data()));
			delete canvasSingleBin;
		}
	}               
		
	TCanvas* canvasRatioZInRSP = new TCanvas("canvasRatioZInRSP","",400,20,2000,1400);  // gives the page size
	
	TPad* padRatioZInRSP = new TPad("padRatioZInRSP","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
	padRatioZInRSP->SetFillColor(0);
	padRatioZInRSP->GetFrame()->SetFillColor(0);
	padRatioZInRSP->SetBorderMode(0);
	padRatioZInRSP->Divide(rowR-3,columnR+2);
	padRatioZInRSP->Draw();
	
	for(Int_t iR = iRStart; iR < nBinsR; iR++){
			Int_t place = iR + 1;
			
			padRatioZInRSP->cd(place);
			if (boolRatioLog) padRatioZInRSP->cd(place)->SetLogy(1);
			padRatioZInRSP->cd(place)->SetRightMargin(0.01);
			nameHistoRatioZInRMonteCarlo=Form("Z in R:  %s",arrayNamesRBins[iR].Data());
			Float_t rangeZ = arrayRangeZInR[iR];
			DrawRatioGammaHisto( histoMappingRatioZInR[iR], 
											nameHistoRatioZInRMonteCarlo,"Z","norm Data/norm MC",
											kFALSE,3 ,0.000001,
											kTRUE, minYValueRatio , maxYValueRatio,
											kTRUE, -rangeZ,rangeZ);
			DrawGammaLines(-rangeZ,rangeZ,1,1,0.02);
// 			if (iR == iRStart){DrawAliceText(floatLocationRightUpText2[0],floatLocationRightUpText2[1], floatLocationRightUpText2[3]);}
	}
	padRatioZInRSP->Update();
	canvasRatioZInRSP->SaveAs(Form("%s/Ratio_Z_in_R.%s",outputDirectory.Data(),suffix.Data()));
	delete padRatioZInRSP ;
	delete canvasRatioZInRSP;
	
	for(Int_t iR = iRStart; iR < nBinsR; iR++){
		if (histoMappingRatioZInR[iR]->GetEntries()>0){	
			TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
			if (boolRatioLog) canvasSingleBin->SetLogy(1);
			canvasSingleBin->cd();
			Float_t rangeZ = arrayRangeZInR[iR];
			TString nameHistoPhiInRData=Form("Z in R:  %s",arrayNamesRBins[iR].Data());
			TLatex *latexBinning = new TLatex(0.16,0.76,nameHistoPhiInRData); // Bo: this was modified
			latexBinning->SetNDC();
			latexBinning->SetTextColor(1);
			latexBinning->SetTextFont(62);
			latexBinning->SetTextSize(0.04);
			latexBinning->SetLineWidth(2);       
			
			DrawRatioGammaHisto( histoMappingRatioZInR[iR], 
											"","Z","norm Data/norm MC",
											kFALSE,3 ,0.000001,
											kTRUE, minYValueRatio , maxYValueRatio,
											kTRUE, -rangeZ,rangeZ);
			DrawGammaLines(-rangeZ,rangeZ,1,1,0.02);                
			latexBinning->Draw("same");  
// 			DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
			
			canvasSingleBin->Update();
			canvasSingleBin->SaveAs(Form("%s/Ratio_Z_in_R_%i.%s",outputDirectory.Data(),iR,suffix.Data()));
			delete canvasSingleBin;
		}
	}
	
	TCanvas* canvasRatioPhiInZSP = new TCanvas("canvasRatioPhiInZSP","",400,20,2000,1400);  // gives the page size
	
	TPad* padRatioPhiInZSP = new TPad("padRatioPhiInZSP","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
	padRatioPhiInZSP->SetFillColor(0);
	padRatioPhiInZSP->GetFrame()->SetFillColor(0);
	padRatioPhiInZSP->SetBorderMode(0);
	padRatioPhiInZSP->Divide(nRowsZProjections-1 ,nColumnsZProjections+1 );
	padRatioPhiInZSP->Draw();
	
	for(Int_t iZ = minimumBinZNormal ; iZ < (nBinsZ-minimumBinZNormal); iZ++){
		Int_t place = iZ -minimumBinZNormal+ 1;
		padRatioPhiInZSP->cd(place);
		if (boolRatioLog) padRatioPhiInZSP->cd(place)->SetLogy(1);
		padRatioPhiInZSP->cd(place)->SetRightMargin(0.01);
		nameHistoRatioPhiInZMonteCarlo=Form("#Phi in Z:  %s",arrayNamesZBins[iZ].Data());
		DrawRatioGammaHisto( histoMappingRatioPhiInZ[iZ], 
										nameHistoRatioPhiInZMonteCarlo,"#Phi","norm Data/norm MC",
										kFALSE,3 ,0.000001,
										kTRUE, minYValueRatio , maxYValueRatio,
										kFALSE, 0,0);
		linePhi->Draw("same");
// 		if (iZ == 0){DrawAliceText(floatLocationRightUpText2[0],floatLocationRightUpText2[1], floatLocationRightUpText2[3]);}
	}
	padRatioPhiInZSP->Update();
	canvasRatioPhiInZSP->SaveAs(Form("%s/Ratio_Phi_in_Z.%s",outputDirectory.Data(),suffix.Data()));
	delete padRatioPhiInZSP ;
	delete canvasRatioPhiInZSP;
	
	for(Int_t iZ = minimumBinZNormal ; iZ < (nBinsZ-minimumBinZNormal); iZ++){
		TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
		if (boolRatioLog) canvasSingleBin->SetLogy(1);
		canvasSingleBin->cd();
		TString nameHistoPhiInRData=Form("#Phi in Z:  %s",arrayNamesZBins[iZ].Data());
		TLatex *latexBinning = new TLatex(0.16,0.76,nameHistoPhiInRData); // Bo: this was modified
		latexBinning->SetNDC();
		latexBinning->SetTextColor(1);
		latexBinning->SetTextFont(62);
		latexBinning->SetTextSize(0.04);
		latexBinning->SetLineWidth(2);       
		
		DrawRatioGammaHisto( histoMappingRatioPhiInZ[iZ], 
										"","#Phi","norm Data/norm MC",
										kFALSE,3 ,0.000001,
										kTRUE, minYValueRatio , maxYValueRatio,
										kFALSE, 0,0);
		DrawGammaLines(-3.2,3.2,1,1,0.02);              
		latexBinning->Draw("same");  
// 		DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
		canvasSingleBin->Update();
		canvasSingleBin->SaveAs(Form("%s/Ratio_Phi_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
		delete canvasSingleBin;
	}               
		
	TCanvas* canvasSPDRatioPhiInZSP = new TCanvas("canvasSPDRatioPhiInZSP","",400,20,2000,1400);  // gives the page size
	
	TPad* padSPDRatioPhiInZSP = new TPad("padSPDRatioPhiInZSP","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
	padSPDRatioPhiInZSP->SetFillColor(0);
	padSPDRatioPhiInZSP->GetFrame()->SetFillColor(0);
	padSPDRatioPhiInZSP->SetBorderMode(0);
	padSPDRatioPhiInZSP->Divide(nRowsZSPD-1 ,nColumnsZSPD+1);
	padSPDRatioPhiInZSP->Draw();
	
	for(Int_t iZ = minimumBinZSPD ; iZ < (nBinsZ-minimumBinZSPD); iZ++){
		Int_t place = iZ -minimumBinZSPD + 1;
		padSPDRatioPhiInZSP->cd(place);
		if (boolRatioLog) padSPDRatioPhiInZSP->cd(place)->SetLogy(1);
		padSPDRatioPhiInZSP->cd(place)->SetRightMargin(0.01);
		nameHistoRatioSPDPhiInZMonteCarlo=Form("0 cm < R < 25 cm: #phi in Z:  %s",arrayNamesZBins[iZ].Data());
		DrawRatioGammaHisto( histoMappingRatioSPDPhiInZ[iZ], 
										nameHistoRatioSPDPhiInZMonteCarlo,"#Phi","norm Data/norm MC",
										kFALSE,3 ,0.000001,
										kTRUE, minYValueRestRRatio , maxYValueRestRRatio,
										kFALSE, 0,0);
		linePhi->Draw("same");
// 		if (iZ == 0){DrawAliceText(floatLocationRightUpText2[0],floatLocationRightUpText2[1], floatLocationRightUpText2[3]);}
	}
	
	padSPDRatioPhiInZSP->Update();
	canvasSPDRatioPhiInZSP->SaveAs(Form("%s/Ratio_SPD_Phi_in_Z.%s",outputDirectory.Data(),suffix.Data()));
	delete padSPDRatioPhiInZSP ;
	delete canvasSPDRatioPhiInZSP;
	
	for(Int_t iZ = minimumBinZSPD ; iZ < (nBinsZ-minimumBinZSPD); iZ++){
		TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
		if (boolRatioLog) canvasSingleBin->SetLogy(1);
		canvasSingleBin->cd();
		TString nameHistoPhiInRData=Form("0 cm < R < 25 cm: #phi in Z:  %s",arrayNamesZBins[iZ].Data());
		TLatex *latexBinning = new TLatex(0.16,0.76,nameHistoPhiInRData); // Bo: this was modified
		latexBinning->SetNDC();
		latexBinning->SetTextColor(1);
		latexBinning->SetTextFont(62);
		latexBinning->SetTextSize(0.04);
		latexBinning->SetLineWidth(2);       
		DrawRatioGammaHisto( histoMappingRatioSPDPhiInZ[iZ], 
										"","#Phi","norm Data/norm MC",
										kFALSE,3 ,0.000001,
										kTRUE, minYValueRestRRatio , maxYValueRestRRatio,
										kFALSE, 0,0);
		DrawGammaLines(-3.2,3.2,1,1,0.02);              
		latexBinning->Draw("same");  
// 		DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
		canvasSingleBin->Update();
		canvasSingleBin->SaveAs(Form("%s/Ratio_SPD_Phi_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
		delete canvasSingleBin;
	}               

	TCanvas* canvasSDDRatioPhiInZSP = new TCanvas("canvasSDDRatioPhiInZSP","",400,20,2000,1400);  // gives the page size
	
	TPad* padSDDRatioPhiInZSP = new TPad("padSDDRatioPhiInZSP","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
	padSDDRatioPhiInZSP->SetFillColor(0);
	padSDDRatioPhiInZSP->GetFrame()->SetFillColor(0);
	padSDDRatioPhiInZSP->SetBorderMode(0);
	padSDDRatioPhiInZSP->Divide(nRowsZSDD-1 ,nColumnsZSDD+1);
	padSDDRatioPhiInZSP->Draw();
	
	for(Int_t iZ = minimumBinZSDD ; iZ < (nBinsZ-minimumBinZSDD); iZ++){
		Int_t place = iZ -minimumBinZSDD + 1;
		padSDDRatioPhiInZSP->cd(place);
		if (boolRatioLog) padSDDRatioPhiInZSP->cd(place)->SetLogy(1);
		padSDDRatioPhiInZSP->cd(place)->SetRightMargin(0.01);
		nameHistoRatioSDDPhiInZMonteCarlo=Form("25 cm < R < 50 cm: #phi in Z:  %s",arrayNamesZBins[iZ].Data());
		DrawRatioGammaHisto( histoMappingRatioSDDPhiInZ[iZ], 
										nameHistoRatioSDDPhiInZMonteCarlo,"#Phi","norm Data/norm MC",
										kFALSE,3 ,0.000001,
										kTRUE, minYValueRestRRatio , maxYValueRestRRatio,
										kFALSE, 0,0);
		linePhi->Draw("same");
// 		if (iZ == 0){DrawAliceText(floatLocationRightUpText2[0],floatLocationRightUpText2[1], floatLocationRightUpText2[3]);}
	}
	
	padSDDRatioPhiInZSP->Update();
	canvasSDDRatioPhiInZSP->SaveAs(Form("%s/Ratio_SDD_Phi_in_Z.%s",outputDirectory.Data(),suffix.Data()));
	delete padSDDRatioPhiInZSP ;
	delete canvasSDDRatioPhiInZSP;
	
	for(Int_t iZ = minimumBinZSDD ; iZ < (nBinsZ-minimumBinZSDD); iZ++){
		TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
		if (boolRatioLog) canvasSingleBin->SetLogy(1);
		canvasSingleBin->cd();
		TString nameHistoPhiInRData=Form("25 cm < R < 50 cm: #phi in Z:  %s",arrayNamesZBins[iZ].Data());
		TLatex *latexBinning = new TLatex(0.16,0.76,nameHistoPhiInRData); // Bo: this was modified
		latexBinning->SetNDC();
		latexBinning->SetTextColor(1);
		latexBinning->SetTextFont(62);
		latexBinning->SetTextSize(0.04);
		latexBinning->SetLineWidth(2);       
		DrawRatioGammaHisto( histoMappingRatioSDDPhiInZ[iZ], 
										"","#Phi","norm Data/norm MC",
										kFALSE,3 ,0.000001,
										kTRUE, minYValueRestRRatio , maxYValueRestRRatio,
										kFALSE, 0,0);
		DrawGammaLines(-3.2,3.2,1,1,0.02);              
		latexBinning->Draw("same");  
// 		DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
		canvasSingleBin->Update();
		canvasSingleBin->SaveAs(Form("%s/Ratio_SDD_Phi_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
		delete canvasSingleBin;
	}               

			
	TCanvas* canvasITSTPCRatioPhiInZSP = new TCanvas("canvasITSTPCRatioPhiInZSP","",400,20,2000,1400);  // gives the page size
	
	TPad* padITSTPCRatioPhiInZSP = new TPad("padITSTPCRatioPhiInZSP","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
	padITSTPCRatioPhiInZSP->SetFillColor(0);
	padITSTPCRatioPhiInZSP->GetFrame()->SetFillColor(0);
	padITSTPCRatioPhiInZSP->SetBorderMode(0);
	padITSTPCRatioPhiInZSP->Divide(nRowsZITSTPC-2 ,nColumnsZITSTPC+2);
	padITSTPCRatioPhiInZSP->Draw();
	
	
	for(Int_t iZ = minimumBinZITSTPC ; iZ < (nBinsZ-minimumBinZITSTPC); iZ++){
		Int_t place = iZ -minimumBinZITSTPC + 1;
		padITSTPCRatioPhiInZSP->cd(place);
		if (boolRatioLog) padITSTPCRatioPhiInZSP->cd(place)->SetLogy(1);
		padITSTPCRatioPhiInZSP->cd(place)->SetRightMargin(0.01);
		nameHistoRatioITSTPCPhiInZMonteCarlo=Form("50 cm < R < 80 cm: #phi in Z:  %s",arrayNamesZBins[iZ].Data());
		DrawRatioGammaHisto( histoMappingRatioITSTPCPhiInZ[iZ], 
						nameHistoRatioITSTPCPhiInZMonteCarlo,"#Phi","norm Data/norm MC",
						kFALSE,3 ,0.000001,
						kTRUE, minYValueRestRRatio , maxYValueRestRRatio,
						kFALSE, 0,0);
						linePhi->Draw("same");
// 						if (iZ == 0){DrawAliceText(floatLocationRightUpText2[0],floatLocationRightUpText2[1], floatLocationRightUpText2[3]);}
	}
	
	padITSTPCRatioPhiInZSP->Update();
	canvasITSTPCRatioPhiInZSP->SaveAs(Form("%s/Ratio_ITSTPC_Phi_in_Z.%s",outputDirectory.Data(),suffix.Data()));
	delete padITSTPCRatioPhiInZSP ;
	delete canvasITSTPCRatioPhiInZSP;
	
	for(Int_t iZ = minimumBinZITSTPC ; iZ < (nBinsZ-minimumBinZITSTPC); iZ++){
		TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
		if (boolRatioLog) canvasSingleBin->SetLogy(1);
		canvasSingleBin->cd();
		TString nameHistoPhiInRData=Form("50 cm < R < 80 cm: #phi in Z:  %s",arrayNamesZBins[iZ].Data());
		TLatex *latexBinning = new TLatex(0.16,0.76,nameHistoPhiInRData); // Bo: this was modified
		latexBinning->SetNDC();
		latexBinning->SetTextColor(1);
		latexBinning->SetTextFont(62);
		latexBinning->SetTextSize(0.04);
		latexBinning->SetLineWidth(2);       
		DrawRatioGammaHisto( histoMappingRatioITSTPCPhiInZ[iZ], 
						"","#Phi","norm Data/norm MC",
						kFALSE,3 ,0.000001,
						kTRUE, minYValueRestRRatio , maxYValueRestRRatio,
						kFALSE, 0,0);
		DrawGammaLines(-3.2,3.2,1,1,0.02);              
		latexBinning->Draw("same");  
// 		DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
		canvasSingleBin->Update();
		canvasSingleBin->SaveAs(Form("%s/Ratio_ITSTPC_Phi_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
		delete canvasSingleBin;
	}               

	TCanvas* canvasHotZoneRatioPhiInZSP = new TCanvas("canvasHotZoneRatioPhiInZSP","",400,20,2000,1400);  // gives the page size
	
	TPad* padHotZoneRatioPhiInZSP = new TPad("padHotZoneRatioPhiInZSP","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
	padHotZoneRatioPhiInZSP->SetFillColor(0);
	padHotZoneRatioPhiInZSP->GetFrame()->SetFillColor(0);
	padHotZoneRatioPhiInZSP->SetBorderMode(0);
	padHotZoneRatioPhiInZSP->Divide(nRowsZHotZone-2 ,nColumnsZHotZone+2);
	padHotZoneRatioPhiInZSP->Draw();
	
	
	for(Int_t iZ = minimumBinZHotZone ; iZ < (nBinsZHotZone-minimumBinZHotZone); iZ++){
		Int_t place = iZ -minimumBinZHotZone + 1;
		padHotZoneRatioPhiInZSP->cd(place);
		if (boolRatioLog) padHotZoneRatioPhiInZSP->cd(place)->SetLogy(1);
		padHotZoneRatioPhiInZSP->cd(place)->SetRightMargin(0.01);
		nameHistoRatioHotZonePhiInZMonteCarlo=Form("5.7 cm < R < 50 cm: #phi in Z:  %s",arrayNamesZBinsHotZone[iZ].Data());
		DrawRatioGammaHisto( histoMappingRatioHotZonePhiInZ[iZ], 
						nameHistoRatioHotZonePhiInZMonteCarlo,"#Phi","norm Data/norm MC",
						kFALSE,3 ,0.000001,
						kTRUE, minYValueRestRRatio , maxYValueRestRRatio,
						kFALSE, 0,0);
						linePhi->Draw("same");
// 						if (iZ == 0){DrawAliceText(floatLocationRightUpText2[0],floatLocationRightUpText2[1], floatLocationRightUpText2[3]);}
	}
	
	padHotZoneRatioPhiInZSP->Update();
	canvasHotZoneRatioPhiInZSP->SaveAs(Form("%s/Ratio_HotZone_Phi_in_Z.%s",outputDirectory.Data(),suffix.Data()));
	delete padHotZoneRatioPhiInZSP ;
	delete canvasHotZoneRatioPhiInZSP;
	
	for(Int_t iZ = minimumBinZHotZone ; iZ < (nBinsZHotZone-minimumBinZHotZone); iZ++){
		TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
		if (boolRatioLog) canvasSingleBin->SetLogy(1);
		canvasSingleBin->cd();
		TString nameHistoPhiInRData=Form("5.7 cm < R < 50 cm: #phi in Z:  %s",arrayNamesZBinsHotZone[iZ].Data());
		TLatex *latexBinning = new TLatex(0.16,0.76,nameHistoPhiInRData); // Bo: this was modified
		latexBinning->SetNDC();
		latexBinning->SetTextColor(1);
		latexBinning->SetTextFont(62);
		latexBinning->SetTextSize(0.04);
		latexBinning->SetLineWidth(2);       
		DrawRatioGammaHisto( histoMappingRatioHotZonePhiInZ[iZ], 
						"","#Phi","norm Data/norm MC",
						kFALSE,3 ,0.000001,
						kTRUE, minYValueRestRRatio , maxYValueRestRRatio,
						kFALSE, 0,0);
		DrawGammaLines(-3.2,3.2,1,1,0.02);              
		latexBinning->Draw("same");  
// 		DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
		canvasSingleBin->Update();
		canvasSingleBin->SaveAs(Form("%s/Ratio_HotZone_Phi_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
		delete canvasSingleBin;
	}               

		
	TCanvas* canvasRatioRInZSP = new TCanvas("canvasRatioRInZSP","",400,20,2000,1400);  // gives the page size
	
	TPad* padRatioRInZSP = new TPad("padRatioRInZSP","",0.01,0.01,0.99,0.99,0);   // gives the size of the histo areas 
	padRatioRInZSP->SetFillColor(0);
	padRatioRInZSP->GetFrame()->SetFillColor(0);
	padRatioRInZSP->SetBorderMode(0);
	padRatioRInZSP->Divide(nRowsZProjections-1 ,nColumnsZProjections+1 );
	padRatioRInZSP->Draw();
	
	for(Int_t iZ = minimumBinZNormal ; iZ < (nBinsZ-minimumBinZNormal); iZ++){
		Int_t place = iZ -minimumBinZNormal+ 1;
			padRatioRInZSP->cd(place);
			if (boolRatioLog) padRatioRInZSP->cd(place)->SetLogy(1);
			padRatioRInZSP->cd(place)->SetRightMargin(0.01);
			nameHistoRatioRInZMonteCarlo=Form("R in Z:  %s",arrayNamesZBins[iZ].Data());
			DrawRatioGammaHisto( histoMappingRatioRInZ[iZ], 
											nameHistoRatioRInZMonteCarlo,"R (cm)","norm Data/norm MC",
											kFALSE,3 ,0.000001,
											kTRUE, minYValueRatio , maxYValueRatio,
											kFALSE, 0,0);
			lineR->Draw("same");
// 			if (iZ == 0){DrawAliceText(floatLocationRightUpText2[0],floatLocationRightUpText2[1], floatLocationRightUpText2[3]);}
	}
	padRatioRInZSP->Update();
	canvasRatioRInZSP->SaveAs(Form("%s/Ratio_R_in_Z.%s",outputDirectory.Data(),suffix.Data()));
	delete padRatioRInZSP ;
	delete canvasRatioRInZSP;
	
	for(Int_t iZ = minimumBinZNormal ; iZ < (nBinsZ-minimumBinZNormal); iZ++){
		TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
		if (boolRatioLog) canvasSingleBin->SetLogy(1);
		canvasSingleBin->cd();
		TString nameHistoPhiInRData=Form("R in Z:  %s",arrayNamesZBins[iZ].Data());
		TLatex *latexBinning = new TLatex(0.16,0.76,nameHistoPhiInRData); // Bo: this was modified
		latexBinning->SetNDC();
		latexBinning->SetTextColor(1);
		latexBinning->SetTextFont(62);
		latexBinning->SetTextSize(0.04);
		latexBinning->SetLineWidth(2);       
		
		DrawRatioGammaHisto( histoMappingRatioRInZ[iZ], 
										"","R (cm)","norm Data/norm MC",
										kFALSE,3 ,0.000001,
										kTRUE, minYValueRatio , maxYValueRatio,
										kFALSE, 0,0);
		DrawGammaLines(0,200,1,1,0.02);         
		latexBinning->Draw("same");  
// 		DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
		canvasSingleBin->Update();
		canvasSingleBin->SaveAs(Form("%s/Ratio_R_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
		delete canvasSingleBin;
	}               
		
	TCanvas * canvasChi2NDOF = new TCanvas("canvasChi2NDOF","",10,10,500,500);  // gives the page size
	canvasChi2NDOF->SetLogy(1);
	canvasChi2NDOF->SetLeftMargin(0.16);
	canvasChi2NDOF->cd();
	DrawAutoGammaHistos( histoHadIntQualChi2NDOFData,
						histoHadIntQualChi2NDOFMonteCarlo, 
						"", "#chi^{2}/n_{dof}",textStandardYAxis2,
						kTRUE, 3,1e-8,
						kFALSE,0. ,0., 
						kFALSE, 0.1,40.);
	
	canvasChi2NDOF->Update();
	canvasChi2NDOF->SaveAs(Form("%s/HadInt_Chi2NDOF.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasChi2NDOF;


	TCanvas * canvasChi2NDOFAftCuts = new TCanvas("canvasChi2NDOFAftCuts","",10,10,500,500);  // gives the page size
	canvasChi2NDOFAftCuts->SetLogy(1);
	canvasChi2NDOFAftCuts->SetLeftMargin(0.16);
	canvasChi2NDOFAftCuts->cd();
	DrawAutoGammaHistos( histoHadIntQualChi2NDOFAftCutsData,
						histoHadIntQualChi2NDOFAftCutsMonteCarlo, 
						"", "#chi^{2}/n_{dof}",textStandardYAxis2,
						kTRUE, 3,1e-2,
						kFALSE,0. ,0., 
						kTRUE, 0.,5.);
	DrawLabelsEvents(floatLocationRightUp2D[0]+0.09,floatLocationRightUp2D[1]-0.1,floatLocationRightUp2D[2], 0.00, collisionSystem, textGenerator,textPeriod);
	canvasChi2NDOFAftCuts->Update();
	canvasChi2NDOFAftCuts->SaveAs(Form("%s/HadInt_Chi2NDOF_AftCuts.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasChi2NDOFAftCuts;

	
	TCanvas * canvasErr2D = new TCanvas("canvasErr2D","",10,10,500,500);  // gives the page size
	canvasErr2D->SetLogy(1);
	canvasErr2D->SetLeftMargin(0.16);
	canvasErr2D->cd();
	DrawAutoGammaHistos( histoHadIntQualErr2DData,
						histoHadIntQualErr2DMonteCarlo, 
						"", "#Delta_{2D}^{sec. vtx}",textStandardYAxis2,
						kTRUE, 3,5e-6,
						kFALSE,0. ,0., 
						kTRUE, 0.,4.);
	
	canvasErr2D->Update();
	canvasErr2D->SaveAs(Form("%s/HadInt_Err2D.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasErr2D;


	TCanvas * canvasErr2DAftCuts = new TCanvas("canvasErr2DAftCuts","",10,10,500,500);  // gives the page size
	canvasErr2DAftCuts->SetLogy(1);
	canvasErr2DAftCuts->SetLeftMargin(0.16);
	canvasErr2DAftCuts->cd();
	DrawAutoGammaHistos( histoHadIntQualErr2DAftCutsData,
						histoHadIntQualErr2DAftCutsMonteCarlo, 
						"", "#Delta_{2D}^{sec. vtx}",textStandardYAxis2,
						kTRUE, 3,8e-4,
						kFALSE,0. ,0., 
						kTRUE, 0.,2.);
	DrawLabelsEvents(floatLocationRightUp2D[0]+0.09,floatLocationRightUp2D[1]-0.1,floatLocationRightUp2D[2], 0.00, collisionSystem, textGenerator,textPeriod);
	canvasErr2DAftCuts->Update();
	canvasErr2DAftCuts->SaveAs(Form("%s/HadInt_Err2D_AftCuts.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasErr2DAftCuts;

	TCanvas * canvasErr3D = new TCanvas("canvasErr3D","",10,10,500,500);  // gives the page size
	canvasErr3D->SetLogy(1);
	canvasErr3D->SetLeftMargin(0.16);
	canvasErr3D->cd();
	DrawAutoGammaHistos( histoHadIntQualErr3DData,
						histoHadIntQualErr3DMonteCarlo, 
						"", "#Delta_{3D}^{sec. vtx}",textStandardYAxis2,
						kTRUE, 3,5e-6,
						kFALSE,0. ,0., 
						kTRUE, 0.,4.);
	
	canvasErr3D->Update();
	canvasErr3D->SaveAs(Form("%s/HadInt_Err3D.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasErr3D;

	TCanvas * canvasErr3DAftCuts = new TCanvas("canvasErr3DAftCuts","",10,10,500,500);  // gives the page size
	canvasErr3DAftCuts->SetLogy(1);
	canvasErr3DAftCuts->SetLeftMargin(0.16);
	canvasErr3DAftCuts->cd();
	DrawAutoGammaHistos( histoHadIntQualErr3DAftCutsData,
						histoHadIntQualErr3DAftCutsMonteCarlo, 
						"", "#Delta_{3D}^{sec. vtx}",textStandardYAxis2,
						kTRUE, 3,8e-4,
						kFALSE,0. ,0., 
						kTRUE, 0.,2.);
	DrawLabelsEvents(floatLocationRightUp2D[0]+0.09,floatLocationRightUp2D[1]-0.1,floatLocationRightUp2D[2], 0.00, collisionSystem, textGenerator,textPeriod);
	canvasErr3DAftCuts->Update();
	canvasErr3DAftCuts->SaveAs(Form("%s/HadInt_Err3D_AftCuts.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasErr3DAftCuts;

	
	TCanvas * canvasErrX = new TCanvas("canvasErrX","",10,10,500,500);  // gives the page size
	canvasErrX->SetLogy(1);
	canvasErrX->SetLeftMargin(0.16);
	canvasErrX->cd();
	DrawAutoGammaHistos( histoHadIntQualErrXData,
						histoHadIntQualErrXMonteCarlo, 
						"", "#Delta_{X}^{sec. vtx}",textStandardYAxis2,
						kTRUE, 3,5e-6,
						kFALSE,0. ,0., 
						kTRUE, 0.,4.);
	
	canvasErrX->Update();
	canvasErrX->SaveAs(Form("%s/HadInt_ErrX.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasErrX;

	TCanvas * canvasErrXAftCuts = new TCanvas("canvasErrXAftCuts","",10,10,500,500);  // gives the page size
	canvasErrXAftCuts->SetLogy(1);
	canvasErrXAftCuts->SetLeftMargin(0.16);
	canvasErrXAftCuts->cd();
	DrawAutoGammaHistos( histoHadIntQualErrXAftCutsData,
						histoHadIntQualErrXAftCutsMonteCarlo, 
						"", "#Delta_{X}^{sec. vtx}",textStandardYAxis2,
						kTRUE, 3,5e-5,
						kFALSE,0. ,0., 
						kTRUE, 0.,2.);
	DrawLabelsEvents(floatLocationRightUp2D[0]+0.09,floatLocationRightUp2D[1]-0.1,floatLocationRightUp2D[2], 0.00, collisionSystem, textGenerator,textPeriod);
	canvasErrXAftCuts->Update();
	canvasErrXAftCuts->SaveAs(Form("%s/HadInt_ErrX_AftCuts.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasErrXAftCuts;

	TCanvas * canvasErrY = new TCanvas("canvasErrY","",10,10,500,500);  // gives the page size
	canvasErrY->SetLogy(1);
	canvasErrY->SetLeftMargin(0.16);
	canvasErrY->cd();
	DrawAutoGammaHistos( histoHadIntQualErrYData,
						histoHadIntQualErrYMonteCarlo, 
						"", "#Delta_{Y}^{sec. vtx}",textStandardYAxis2,
						kTRUE, 3,5e-6,
						kFALSE,0. ,0., 
						kTRUE, 0.,4.);
	
	canvasErrY->Update();
	canvasErrY->SaveAs(Form("%s/HadInt_ErrY.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasErrY;

	TCanvas * canvasErrYAftCuts = new TCanvas("canvasErrYAftCuts","",10,10,500,500);  // gives the page size
	canvasErrYAftCuts->SetLogy(1);
	canvasErrYAftCuts->SetLeftMargin(0.16);
	canvasErrYAftCuts->cd();
	DrawAutoGammaHistos( histoHadIntQualErrYAftCutsData,
						histoHadIntQualErrYAftCutsMonteCarlo, 
						"", "#Delta_{Y}^{sec. vtx}",textStandardYAxis2,
						kTRUE, 3,5e-5,
						kFALSE,0. ,0., 
						kTRUE, 0.,2.);
	DrawLabelsEvents(floatLocationRightUp2D[0]+0.09,floatLocationRightUp2D[1]-0.1,floatLocationRightUp2D[2], 0.00, collisionSystem, textGenerator,textPeriod);
	canvasErrYAftCuts->Update();
	canvasErrYAftCuts->SaveAs(Form("%s/HadInt_ErrY_AftCuts.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasErrYAftCuts;

	
	TCanvas * canvasErrZ = new TCanvas("canvasErrZ","",10,10,500,500);  // gives the page size
	canvasErrZ->SetLogy(1);
	canvasErrZ->SetLeftMargin(0.16);
	canvasErrZ->cd();
	DrawAutoGammaHistos( histoHadIntQualErrZData,
						histoHadIntQualErrZMonteCarlo, 
						"", "#Delta_{Z}^{sec. vtx}",textStandardYAxis2,
						kTRUE, 3,5e-6,
						kFALSE,0. ,0., 
						kTRUE, 0.,4.);
	
	canvasErrZ->Update();
	canvasErrZ->SaveAs(Form("%s/HadInt_ErrZ.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasErrZ;

	TCanvas * canvasErrZAftCuts = new TCanvas("canvasErrZAftCuts","",10,10,500,500);  // gives the page size
	canvasErrZAftCuts->SetLogy(1);
	canvasErrZAftCuts->SetLeftMargin(0.16);
	canvasErrZAftCuts->cd();
	DrawAutoGammaHistos( histoHadIntQualErrZAftCutsData,
						histoHadIntQualErrZAftCutsMonteCarlo, 
						"", "#Delta_{Z}^{sec. vtx}",textStandardYAxis2,
						kTRUE, 3,5e-5,
						kFALSE,0. ,0., 
						kTRUE, 0.,2.);
	DrawLabelsEvents(floatLocationRightUp2D[0]+0.09,floatLocationRightUp2D[1]-0.1,floatLocationRightUp2D[2], 0.00, collisionSystem, textGenerator, textPeriod);
	canvasErrZAftCuts->Update();
	canvasErrZAftCuts->SaveAs(Form("%s/HadInt_ErrZ_AftCuts.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasErrZAftCuts;

	TCanvas * canvasNTracks = new TCanvas("canvasNTracks","",10,10,500,500);  // gives the page size
	canvasNTracks->SetLogy(1);
	canvasNTracks->SetLeftMargin(0.16);
	canvasNTracks->cd();
	DrawAutoGammaHistos( histoHadIntQualNTracksData,
						histoHadIntQualNTracksMonteCarlo, 
						"", "N_{Tracks per Vtx}",textStandardYAxis2,
						kTRUE, 3,2.5e-7,
						kFALSE,0. ,0., 
						kFALSE, 0.,2.);
	DrawLabelsEvents(floatLocationRightUp2D[0]+0.09,floatLocationRightUp2D[1]-0.1,floatLocationRightUp2D[2], 0.00, collisionSystem, textGenerator,textPeriod);
	canvasNTracks->Update();
	canvasNTracks->SaveAs(Form("%s/HadInt_nTracksPerVtx.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasNTracks;

	
// TFile* fileMappingOutputData = new TFile(Form("../MappingOutputData_%s.root",cutSel.Data()),"UPDATE");
// 	if (textGenerator.CompareTo("Phojet")==0){
// 		histoHadIntVsRData->GetXaxis()->SetRangeUser(0.,180.);
// 		histoHadIntVsRData->Write(Form("RDistribution_Data_%s",textPeriod.Data()),TObject::kOverwrite);
// 		histoHadIntVsZData->Write(Form("ZDistribution_Data_%s",textPeriod.Data()),TObject::kOverwrite);
// 	}
// 	histoHadIntVsRMonteCarlo->GetXaxis()->SetRangeUser(0.,180.);
// 	histoHadIntVsRMonteCarlo->Write(Form("RDistribution_MC_%s_%s",textGenerator.Data(),textPeriod.Data()),TObject::kOverwrite);
// 	histoHadIntVsZMonteCarlo->Write(Form("ZDistribution_MC_%s_%s",textGenerator.Data(),textPeriod.Data()),TObject::kOverwrite);
// 	fileMappingOutputData->Write();
// 	fileMappingOutputData->Close();
// 
// 	TFile* fileMappingDetailed = new TFile(Form("%s/%s/MappingHadronicDetailed.root",cutSel.Data(),optEnergy.Data()),"UPDATE");
// 	fileMappingDetailed->mkdir(textPeriod.Data());
// 	fileMappingDetailed->cd(textPeriod.Data());
// 	for(Int_t iR = iRStart; iR < nBinsR; iR++){
// 		histoMappingHadIntPhiInRData[iR]->Write(Form("Data_Conversion_Mapping_Phi_in_R_%02d",iR),TObject::kOverwrite);
// 		histoMappingHadIntPhiInRMonteCarlo[iR]->Write(Form("MC_Conversion_Mapping_Phi_in_R_%02d",iR),TObject::kOverwrite);
// 		histoMappingHadIntZInRData[iR]->Write(Form("Data_Conversion_Mapping_Z_in_R_%02d",iR),TObject::kOverwrite);
// 		histoMappingHadIntZInRMonteCarlo[iR]->Write(Form("MC_Conversion_Mapping_Z_in_R_%02d",iR),TObject::kOverwrite);
// 	}
// 	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
// 		histoMappingHadIntPhiInZData[iZ]->Write(Form("Data_Conversion_Mapping_Phi_in_Z_%02d",iZ),TObject::kOverwrite);
// 		histoMappingHadIntPhiInZMonteCarlo[iZ]->Write(Form("MC_Conversion_Mapping_Phi_in_Z_%02d",iZ),TObject::kOverwrite);
// 		histoMappingHadIntSPDPhiInZData[iZ]->Write(Form("Data_Conversion_Mapping_SPD_Phi_in_Z_%02d",iZ),TObject::kOverwrite);
// 		histoMappingHadIntSPDPhiInZMonteCarlo[iZ]->Write(Form("MC_Conversion_Mapping_SPD_Phi_in_Z_%02d",iZ),TObject::kOverwrite);
// 		histMappingHadIntSDDPhiInZData[iZ]->Write(Form("Data_Conversion_Mapping_SDD_Phi_in_Z_%02d",iZ),TObject::kOverwrite);
// 		histoMappingHadIntSDDPhiInZMonteCarlo[iZ]->Write(Form("MC_Conversion_Mapping_SDD_Phi_in_Z_%02d",iZ),TObject::kOverwrite);
// 		histoMappingHadIntITSTPCPhiInZData[iZ]->Write(Form("Data_Conversion_Mapping_ITSTPC_Phi_in_Z_%02d",iZ),TObject::kOverwrite);
// 		histoMappingHadIntITSTPCPhiInZMonteCarlo[iZ]->Write(Form("MC_Conversion_Mapping_ITSTPC_Phi_in_Z_%02d",iZ),TObject::kOverwrite);
// 		histoMappingHadIntRInZData[iZ]->Write(Form("Data_Conversion_Mapping_R_in_Z_%02d",iZ),TObject::kOverwrite);
// 		histoMappingHadIntRInZMonteCarlo[iZ]->Write(Form("MC_Conversion_Mapping_R_in_Z_%02d",iZ),TObject::kOverwrite);
// 	}
// 	histoHadIntVsRData->GetXaxis()->SetRangeUser(0.,180.);
// 	histoHadIntVsRData->Write("RDistribution_Data",TObject::kOverwrite);
// 	histoHadIntVsRMonteCarlo->GetXaxis()->SetRangeUser(0.,180.);
// 	histoHadIntVsRMonteCarlo->Write(Form("RDistribution_MC_%s",textGenerator.Data()),TObject::kOverwrite);
// 	histoHadIntVsZData->Write("ZDistribution_Data",TObject::kOverwrite);
// 	histoHadIntVsZMonteCarlo->Write(Form("ZDistribution_MC_%s",textGenerator.Data()),TObject::kOverwrite);
// 	fileMappingDetailed->Write();
// 	fileMappingDetailed->Close();
		
	
}
