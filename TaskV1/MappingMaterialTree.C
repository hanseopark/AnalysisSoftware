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
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"

TString textGenerator;
TString collisionSystem;
TString textPeriod;
TString textDate;
Bool_t thesis = kFALSE;

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

Bool_t CheckIfHistoHasEntriesInPlotRange( TH1D* histo, Double_t xStart, Double_t xEnd, Double_t yStart){
	for (Int_t i = histo->GetXaxis()->FindBin(xStart); i < histo->GetXaxis()->FindBin(xEnd); i++){
		if (histo->GetBinContent(i) > yStart ){
			return kTRUE;
		}
	}
	return kFALSE;
}	

Bool_t CheckIfHistoHasEntriesInPlotRange( TH1F* histo, Double_t xStart, Double_t xEnd, Double_t yStart){
	for (Int_t i = histo->GetXaxis()->FindBin(xStart); i < histo->GetXaxis()->FindBin(xEnd); i++){
		if (histo->GetBinContent(i) > yStart ){
			return kTRUE;
		}
	}
	return kFALSE;
}	

void PlotStandard2D( TH2* histo2D, TString nameOutput, TString title, TString xTitle, TString yTitle, Bool_t kRangeY, Double_t startY, Double_t endY, Bool_t kRangeX, Double_t startX, Double_t endX, Bool_t kRangeZ, Double_t startZ, Double_t endZ, Int_t logX, Int_t logZ, Float_t* floatLogo, Int_t canvasSizeX = 500, Int_t canvasSizeY = 500, Double_t rightMargin = 0.12, Double_t leftMargin = 0.12, Double_t bottomMargin = 0.1, Double_t topMargin = 0.04, TString drawLogo ="Perf"){
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
	DrawAutoGammaHisto2D(	histo2D,
									title.Data(), xTitle.Data(), yTitle.Data(),"",kRangeY, startY, endY, kRangeX, startX, endX);
	if(!thesis){
		if (drawLogo.CompareTo("Perf")==0) DrawAliceLogoPerformance2D(floatLogo[0],floatLogo[1],floatLogo[2],floatLogo[3],0.00, textDate,collisionSystem, "","",canvasSizeX,canvasSizeY);	
		if (drawLogo.CompareTo("Wip")==0) DrawAliceLogo(floatLogo[0],floatLogo[1],floatLogo[2],floatLogo[3],canvasSizeX,canvasSizeY);	
	}
	canvasStandard->Update();
	canvasStandard->SaveAs(nameOutput.Data());
	delete canvasStandard;
}


void  MappingMaterialTree(const char *data = "myOutput", TString MCfile = "", TString cutSel = "", TString suffix = "gif", TString optEnergy="", TString optMCGenerator="", TString optPeriod="",TString optThesis =""){        
        
	gROOT->Reset();
	gROOT->SetStyle("Plain");
 	
	collisionSystem = ReturnFullCollisionsSystem(optEnergy);
   if (collisionSystem.CompareTo("") == 0){
      cout << "No correct collision system specification, has been given" << endl;
      return;
   }
   
	if(optMCGenerator.CompareTo("") ==0){
		textGenerator = "";
	} else {
		textGenerator = optMCGenerator;
	}
	TString outputDirectory;
	TString outputDirectory1;
	if(optPeriod.CompareTo("") ==0 || optPeriod.CompareTo("All") ==0){
		if (textGenerator.CompareTo("")!=0){
			outputDirectory =	 		Form("%s/%s/%s/%s/MappingMaterial",textGenerator.Data(),cutSel.Data(),optEnergy.Data(),suffix.Data());
			outputDirectory1 = Form("%s/%s",textGenerator.Data(),cutSel.Data());
		} else {
			textPeriod = "";
			outputDirectory =	 		Form("%s/%s/%s/MappingMaterial",cutSel.Data(),optEnergy.Data(),suffix.Data());
			outputDirectory1 = Form("%s",cutSel.Data());
		}
		gSystem->Exec("mkdir -p "+outputDirectory);
	} else if (textGenerator.CompareTo("")!=0){
		textPeriod = optPeriod;
		outputDirectory =	 		Form("%s/%s/%s/%s/%s/MappingMaterial",textGenerator.Data(),cutSel.Data(),optEnergy.Data(),suffix.Data(),optPeriod.Data());
		outputDirectory1 = Form("%s/%s",textGenerator.Data(),cutSel.Data());
		gSystem->Exec("mkdir -p "+outputDirectory);
	} else {
		textPeriod = optPeriod;
		outputDirectory =	 		Form("%s/%s/%s/%s/MappingMaterial",cutSel.Data(),optEnergy.Data(),suffix.Data(),optPeriod.Data());
		outputDirectory1 = Form("%s/%s",textGenerator.Data(),cutSel.Data());
		gSystem->Exec("mkdir "+outputDirectory);
	}
	
     
	if(optThesis.CompareTo("thesis") == 0){// means we want to plot values for the pi0
		thesis = kTRUE;
	}

	Int_t iRStart = 0;
	StyleSettingsThesis();
	SetPlotStyle();

	// Which file do you want to analyse
	TString 	nameDataFile = 		(Form("%s",data));
	
	//rebinning
	Int_t	rebin = 				4; //4
// 	if(optEnergy.CompareTo("7TeV") == 0){
// 		rebin = 4;
// 	} else {
// 		rebin = 8;
// 	}
		
	//How big should the right margin in 2D plots be? 
	Float_t 	rightMargin = 			0.12;
	Float_t 	leftMargin = 			0.11;
	
	//How many slices do you have?
	const int 	nBinsR = 				13;  
	const int 	nBinsZ = 				18;
	
	// How many raws and colums you want to have in your output?
	Int_t 		columnR = 			2;
	Int_t 		rowR = 				7;
	Int_t 		columnZ = 			3;
	Int_t 		rowZ = 				6;
	
	// Array for ZinR-Ranges
	Float_t	arrayRangeZInR[nBinsR] = {30.,30.,40.,50.,60.,70.,80.,100.,120.,140.,160.,180.,240.}; 
	
	// Array of Rbins
	Float_t 	arrayRBins[14] = 		{0.,3.5,5.75,9.5,13.,21.,27.5,35.,42.,55.,72.,79.5,90.,180.};
	TString 	arrayNamesRangesRBins[13]=	{"0 cm < R < 3.5 cm", 
									"3.5 cm < R < 5.75 cm", 
									"5.75 cm < R < 9.5 cm",
									"9.5 cm < R < 13 cm", 
									"13 cm < R < 21 cm", 
									"21 cm < R < 27.5 cm", 
									"27.5 cm < R < 35 cm", 
									"35 cm < R < 42 cm", 
									"42 cm < R < 55 cm", 
									"55 cm < R < 72 cm", 
									"72 cm < R < 79.5 cm", 
									"79.5 cm < R < 90 cm", 
									"90 cm < R < 180 cm"};								
	TString 	arrayNamesRBins[13]=	{"Beam Pipe", 
									"SPD 1", 
									"SPD 2",
									"Thermal shield/Support between SPD/SDD", 
									"SDD 1", 
									"SDD 2", 
									"Thermal shield/Support between SDD/SSD", 
									"SSD 1", 
									"SSD 2", 
									"Air + TPC in. cont. vessel + CO_{2}", 
									"CO_{2} + TPC in. field cage vessel", 
									"TPC rods + Ne: CO_{2}: N_{2}", 
									"Ne: CO_{2}: N_{2}"};
	TString 	arrayNamesAddRBins[13]=	{"", 
									"", 
									"",
									"", 
									"", 
									"", 
									"", 
									"", 
									"", 
									"+ ITS & TPC services/support ", 
									"+ ITS & TPC services/support", 
									"", 
									""};
	
	// Array of Zbins
	Float_t 	arrayZBins[19] = 		{-300.,-200.,-150.,-100.,-75.,-50.,-30.,-20.,-10.,0.,10.,20.,30.,50.,75.,100.,150.,200.,300.};
	TString 	arrayNamesZBins[18]=	{"-300 cm < Z < -200 cm", 
									"-200 cm < Z < -150 cm", 
									"-150 cm < Z < -100 cm", 
									"-100 cm < Z < -75 cm", 
									"-75 cm < Z < -50 cm", 
									"-50 cm < Z < -30 cm", 
									"-30 cm < Z < -20 cm", 
									"-20 cm < Z < -10 cm", 
									"-10 cm < Z < 0 cm", 
									"0 cm < Z < 10 cm", 
									"10 cm < Z < 20 cm", 
									"20 cm < Z < 30 cm", 
									"30 cm < Z < 50 cm", 
									"50 cm < Z < 75 cm", 
									"75 cm < Z < 100 cm", 
									"100 cm < Z < 150 cm", 
									"150 cm < Z < 200 cm", 
									"200 cm < Z < 300 cm"};
	
	// Array defintion for printing Logo in right upper corner
	Float_t 	floatLocationRightUp[4]= {0.8,0.605,0.09, 0.02};
	Float_t 	floatLocationRightUpR[4]= {0.8,0.55,0.09, 0.02};
	Float_t 	floatLocationRightUp2[4]= {0.8,0.76,0.09, 0.02};
	Float_t 	floatLocationRightDown[4]=	{0.7,0.23,0.09, 0.02};
	Float_t 	floatLocationRightUp2D[4]=	{0.68,0.775,0.10, 0.02};
	Float_t 	floatLocationRightUp2DXY[4]=	{0.72,0.97,0.10, 0.025};
	Float_t 	floatLocationRightUp2DZR[4]=	{0.72,0.4,0.10, 0.025};
	// Array defintion for printing Logo in left upper corner
	Float_t 	floatLocationRight[4] = 	{0.57,0.20,0.04, 0.02};
	Float_t 	floatLocationLeftUp[4]=		{0.18,0.76, 0.09, 0.02};
	Float_t 	floatLocationLeftUpEta[4]=		{0.18,0.8, 0.09, 0.02};
	// Array defintion for printing text in right upper corner
	Float_t 	floatLocationRightUpText[4]=	{0.16,0.92,0.09, 0.04};
	Float_t 	floatLocationRightUpText2[4]=	{0.7,0.8,0.09, 0.04};
	// Get the histos
	TFile 		fileData(nameDataFile);  
	
	// choice of dateset
	TString nameGammaDirectory;
	TString cutSelMat = cutSel;
   cout << cutSelMat.Data() << endl;
	if(cutSelMat.CompareTo("") != 0){
		nameGammaDirectory = 		Form("GammaConv_%s",  cutSelMat.Data());
		cout << nameGammaDirectory.Data() << endl;
	}       
// 	return;
	// labeling
	TString 		textStandardYAxis = 	"#gamma/ event scaled by multiplicity";
	TString 		textStandardYAxis2 = 	"#frac{N_{#gamma}}{N_{ch}}";
	
	textDate = ReturnDateString();

	TString 		textYAxisRHisto = "#frac{1}{N_{ch}} #frac{dN_{#gamma}}{dR} (cm^{-1}) ";
	TString 		textYAxisEtaHisto = "#frac{1}{N_{ch}} #frac{dN_{#gamma}}{d#eta}";
	TString 		textYAxisCEHisto = "#frac{1}{N_{ch}} #frac{dN_{#gamma}}{dR dZ} (cm^{-1}) ";
	TString 		textYAxisZHisto = "#frac{1}{N_{ch}} #frac{dN_{#gamma}}{dZ} (cm^{-1}) ";
	TString 		textYAxisSlicesHisto = "#frac{dN_{#gamma}}{N_{ch}dZ d#varphi dR} (cm^{-2}) ";

	TString etaCutNumber = cutSel(1,1);
	cout << "eta cutnumber = " << etaCutNumber.Data() << endl;
	TLatex *latexEtaRange;
	Bool_t boolRatioLog = kFALSE;
	Double_t minYValueRatio = 0.4;
	Double_t maxYValueRatio = 2.;
	Double_t minYValueRestRRatio = 0.4;
	Double_t maxYValueRestRRatio = 2.;
	Int_t rebinIntPlots = 1;
	Int_t rebinPhiPlots = 2;
	Double_t minYValueZPlot = 1e-6;
	Double_t minYValueRPlot = 1e-6;
	Double_t maxYValueRPlot = 1e-1;
	Double_t minYValueEtaPlot = 1e-6;
	Double_t minYValueRPlotSPD = 1e-6;
	Double_t maxYFactorZPlot = 2.;
	Double_t maxYFactorZPlotlin = 0.8;
	Bool_t outer = kFALSE;
	Double_t rmin = 1.8;
	
	Double_t  floatPtValues7TeVPi0[34] =   {0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,
														2.2,2.4,2.6,2.8,3.,3.2,3.4,3.6,3.8, 4.,4.5,5.,5.5,6.,7.,8.,10.,12.,16.,20.,25};

	Double_t vectorLinPlot[9] = {7.e-4, 7.e-4,4.e-4,4.e-4,5.4e-4,5.4e-4,4.8e-4,4.2e-4,0.7e-4};
	Double_t vectorLogPlot[9] = {1.8e-5,2e-5,3e-6,3e-6,1.5e-6,1.5e-6,1.5e-6,1.5e-6,1.5e-6};
	if (optEnergy.CompareTo("2.76TeV") == 0 && etaCutNumber.CompareTo("0") == 0){
		vectorLinPlot[0] = 7.e-4;
		vectorLinPlot[1] = 7.e-4;
		vectorLinPlot[2] = 6.e-4;
		vectorLinPlot[3] = 6.e-4;
		vectorLinPlot[4] = 8e-4;
		vectorLinPlot[5] = 8e-4;
		vectorLinPlot[6] = 8e-4;
		vectorLinPlot[7] = 6e-4;
		vectorLinPlot[8] = 0.7e-4;		
		
		vectorLogPlot[0] = 1.8e-5;
		vectorLogPlot[1] = 2e-5;
		vectorLogPlot[2] = 3e-6;
		vectorLogPlot[3] = 3e-6;
		vectorLogPlot[4] = 1.5e-6;
		vectorLogPlot[5] = 1.5e-6;
		vectorLogPlot[6] = 1.5e-6;
		vectorLogPlot[7] = 1.5e-6;
		vectorLogPlot[8] = 1.5e-6;
	} else if (optEnergy.CompareTo("2.76TeV") == 0 && etaCutNumber.CompareTo("4") == 0){
		vectorLinPlot[0] = 8.e-4;
		vectorLinPlot[1] = 8.e-4;
		vectorLinPlot[2] = 6.e-4;
		vectorLinPlot[3] = 7.e-4;
		vectorLinPlot[4] = 9e-4;
		vectorLinPlot[5] = 12e-4;
		vectorLinPlot[6] = 1e-4;
		vectorLinPlot[7] = 6e-4;
		vectorLinPlot[8] = 0.7e-4;
		
		vectorLogPlot[0] = 1.8e-5;
		vectorLogPlot[1] = 2e-5;
		vectorLogPlot[2] = 1.5e-6;
		vectorLogPlot[3] = 1.5e-6;
		vectorLogPlot[4] = 1.5e-6;
		vectorLogPlot[5] = 1.5e-6;
		vectorLogPlot[6] = 4e-7;
		vectorLogPlot[7] = 4e-7;
		vectorLogPlot[8] = 4e-7;
	} else if (optEnergy.CompareTo("7TeV") == 0 && etaCutNumber.CompareTo("0") == 0){
		vectorLinPlot[0] = 7.e-4;
		vectorLinPlot[1] = 8.e-4;
		vectorLinPlot[2] = 5.e-4;
		vectorLinPlot[3] = 5.e-4;
		vectorLinPlot[4] = 5.6e-4;
		vectorLinPlot[5] = 5.6e-4;
		vectorLinPlot[6] = 5.6e-4;
		vectorLinPlot[7] = 5e-4;
		vectorLinPlot[8] = 0.7e-4;		
		
		vectorLogPlot[0] = 1.8e-5;
		vectorLogPlot[1] = 2e-5;
		vectorLogPlot[2] = 3e-6;
		vectorLogPlot[3] = 3e-6;
		vectorLogPlot[4] = 1.5e-6;
		vectorLogPlot[5] = 1.5e-6;
		vectorLogPlot[6] = 1.5e-6;
		vectorLogPlot[7] = 1.5e-6;
		vectorLogPlot[8] = 1.5e-6;
	} else if (optEnergy.CompareTo("7TeV") == 0 && etaCutNumber.CompareTo("4") == 0){
		vectorLinPlot[0] = 8.e-4;
		vectorLinPlot[1] = 8.e-4;
		vectorLinPlot[2] = 6.e-4;
		vectorLinPlot[3] = 7.e-4;
		vectorLinPlot[4] = 9e-4;
		vectorLinPlot[5] = 12e-4;
		vectorLinPlot[6] = 1e-4;
		vectorLinPlot[7] = 6e-4;
		vectorLinPlot[8] = 0.7e-4;
		
		vectorLogPlot[0] = 1.8e-5;
		vectorLogPlot[1] = 2e-5;
		vectorLogPlot[2] = 1.5e-6;
		vectorLogPlot[3] = 1.5e-6;
		vectorLogPlot[4] = 1.5e-6;
		vectorLogPlot[5] = 1.5e-6;
		vectorLogPlot[6] = 4e-7;
		vectorLogPlot[7] = 4e-7;
		vectorLogPlot[8] = 4e-7;
	} else if (optEnergy.CompareTo("pPb_5.023TeV") == 0 && etaCutNumber.CompareTo("0") == 0){

		vectorLinPlot[0] = 7.e-4;
		vectorLinPlot[1] = 7.e-4;
		vectorLinPlot[2] = 4.e-4;
		vectorLinPlot[3] = 4.e-4;
		vectorLinPlot[4] = 5.4e-4;
		vectorLinPlot[5] = 5.4e-4;
		vectorLinPlot[6] = 4.8e-4;
		vectorLinPlot[7] = 4.2e-4;
		vectorLinPlot[8] = 0.7e-4;
		
		vectorLogPlot[0] = 1.8e-6;
		vectorLogPlot[1] = 2e-6;
		vectorLogPlot[2] = 3e-7;
		vectorLogPlot[3] = 3e-7;
		vectorLogPlot[4] = 1.5e-7;
		vectorLogPlot[5] = 1.5e-7;
		vectorLogPlot[6] = 1.5e-7;
		vectorLogPlot[7] = 1.5e-7;
		vectorLogPlot[8] = 1.5e-7;	
	}
	
	Double_t maxEta = 0.;
	TString multCutNumbers = cutSel(4,2);
	cout << multCutNumbers.Data() << endl;
	TString nameNTrackHisto = "ESD_NumberOfGoodESDTracks09Vtx";
 	if(etaCutNumber.CompareTo("0") == 0){
			latexEtaRange = 	new TLatex(0.15,0.92,"|#eta| < 0.9 "); 
			latexEtaRange->SetNDC();
			latexEtaRange->SetTextFont(62);
			latexEtaRange->SetTextSize(0.04);
			latexEtaRange->SetLineWidth(6);      
			rebinIntPlots = 2;
			maxEta = 0.9;
			arrayRangeZInR[0] = 15.; arrayRangeZInR[1] = 15.; arrayRangeZInR[2] = 18.; arrayRangeZInR[3] = 20.; arrayRangeZInR[4] = 30.;arrayRangeZInR[5] = 40.;
			arrayRangeZInR[6] = 45.; arrayRangeZInR[7] = 55.; arrayRangeZInR[8] = 60.; arrayRangeZInR[9] = 80.; arrayRangeZInR[10] = 100.; arrayRangeZInR[11] = 100.;
			arrayRangeZInR[12] = 180.; 
			nameNTrackHisto = "ESD_NumberOfGoodESDTracks09Vtx";
				cout << multCutNumbers.Data() << endl;
			if (optEnergy.CompareTo("pPb_5.023TeV") == 0) {
				if (multCutNumbers.CompareTo("00") != 0){
					minYValueRPlot = 1e-7;
					maxYValueRPlot = 3e-3;
				}	
			}	
			
	} else if(etaCutNumber.CompareTo("1") == 0){
			latexEtaRange = 	new TLatex(0.15,0.92,"|#eta| < 0.1 "); 
			latexEtaRange->SetNDC();
			latexEtaRange->SetTextFont(62);
			latexEtaRange->SetTextSize(0.04);
			latexEtaRange->SetLineWidth(6);      
			rebinIntPlots = 2;
			maxEta = 0.1;
	} else if(etaCutNumber.CompareTo("2") == 0){
			latexEtaRange = 	new TLatex(0.15,0.92,"0.9 < |#eta| < 1.4 "); 
			latexEtaRange->SetNDC();
			latexEtaRange->SetTextFont(62);
			latexEtaRange->SetTextSize(0.04);
			latexEtaRange->SetLineWidth(6);
			rebinIntPlots = 4;
			minYValueZPlot = 1e-7;
			minYValueRPlot = 1e-7;
			minYValueRPlotSPD = 1e-7;
			boolRatioLog = kTRUE;
			minYValueRatio = 0.2;
			maxYValueRatio = 10.;
			minYValueRestRRatio = 0.05;
			maxYValueRestRRatio = 10.;
			maxYFactorZPlot = 2. ;
			maxYFactorZPlotlin = 0.8;
			outer = kTRUE;
			maxEta = 1.4;
			nameNTrackHisto = "ESD_NumberOfGoodESDTracks0914Vtx";
	} else if(etaCutNumber.CompareTo("3") == 0){
			latexEtaRange = 	new TLatex(0.15,0.92,"1.4 < |#eta| < 1.8 "); 
			latexEtaRange->SetNDC();
			latexEtaRange->SetTextFont(62);
			latexEtaRange->SetTextSize(0.04);
			latexEtaRange->SetLineWidth(6);
			rebinIntPlots = 4;
			minYValueZPlot = 1e-7;
			minYValueRPlot = 1e-7;
			minYValueRPlotSPD = 1e-7;
			boolRatioLog = kTRUE;
			minYValueRatio = 0.2;
			maxYValueRatio = 10.;
			minYValueRestRRatio = 0.05;
			maxYValueRestRRatio = 10.;
			maxYFactorZPlot = 2. ;
			maxYFactorZPlotlin = 0.8;
			outer = kTRUE;
			maxEta = 1.8;
	} else if(etaCutNumber.CompareTo("4") == 0){
			latexEtaRange = 	new TLatex(0.15,0.92,"|#eta| < 1.4 "); 
			latexEtaRange->SetNDC();
			latexEtaRange->SetTextFont(62);
			latexEtaRange->SetTextSize(0.04);
			latexEtaRange->SetLineWidth(6);
			rebinIntPlots = 4;
			minYValueZPlot = 1e-7;
			minYValueRPlot = 1e-7;
			minYValueRPlotSPD = 1e-7;
			boolRatioLog = kTRUE;
			minYValueRatio = 0.2;
			maxYValueRatio = 10.;
			minYValueEtaPlot = 5e-7;
			minYValueRestRRatio = 0.05;
			maxYValueRestRRatio = 10.;
			maxYFactorZPlot = 2. ;
			maxYFactorZPlotlin = 0.8;
			outer = kTRUE;
			maxEta = 1.4;
			nameNTrackHisto = "ESD_NumberOfGoodESDTracks14Vtx";
	} else if(etaCutNumber.CompareTo("5") == 0){
			latexEtaRange = 	new TLatex(0.15,0.92,"2.5 < |#eta| < 10.0 "); 
			latexEtaRange->SetNDC();
			latexEtaRange->SetTextFont(62);
			latexEtaRange->SetTextSize(0.04);
			latexEtaRange->SetLineWidth(6);
			rebinIntPlots = 4;
			minYValueZPlot = 1e-7;
			minYValueRPlot = 1e-7;
			minYValueRPlotSPD = 1e-7;
			boolRatioLog = kTRUE;
			minYValueRatio = 0.2;
			maxYValueRatio = 10.;
			minYValueRestRRatio = 0.05;
			maxYValueRestRRatio = 10.;
			maxYFactorZPlot = 2. ;
			maxYFactorZPlotlin = 0.8;
			outer = kTRUE;
			maxEta = 10.;
	} else if(etaCutNumber.CompareTo("6") == 0){
			latexEtaRange = 	new TLatex(0.15,0.92,"|#eta| < 10.0 "); 
			latexEtaRange->SetNDC();
			latexEtaRange->SetTextFont(62);
			latexEtaRange->SetTextSize(0.04);
			latexEtaRange->SetLineWidth(6);      
			rebinIntPlots = 2;
			maxEta = 10.;    
	} else if(etaCutNumber.CompareTo("7") == 0){
         latexEtaRange =   new TLatex(0.15,0.92,"|#eta| < 0.65 "); 
         latexEtaRange->SetNDC();
         latexEtaRange->SetTextFont(62);
         latexEtaRange->SetTextSize(0.04);
         latexEtaRange->SetLineWidth(6);      
         rebinIntPlots = 2;
         maxEta = 0.65;
   }
	Double_t doubleLatexNamingBinsX = 0.16;
	Double_t doubleLatexNamingBinsRatioX = 0.11;
	Double_t doubleLatexNamingBinsX2 = 0.24;
	Double_t doubleLatexNamingBinsY = 0.9;
	Double_t doubleLatexNamingBinsY2 = 0.86;
	Double_t doubleLatexNamingBinsY3 = 0.82;
	Size_t sizeTextNameBins = 0.035;
	
	//------------------------------- Reading FILES ----------------------------------------------------------------------
	TDirectory 	*directoryGammaConvData = 	new TDirectory(); 

	directoryGammaConvData = 			(TDirectory*)fileData.Get(nameGammaDirectory); 
	TH1F* histoEventQualityData=					 	(TH1F*)directoryGammaConvData->Get("ESD_EventQuality");
	TH1F* histoVertexZData=					 	(TH1F*)directoryGammaConvData->Get("Z_Vertex_distribution");
	TH1F* histoConversionPointVsEtaData=			(TH1F*)directoryGammaConvData->Get("ESD_ConversionMapping_Eta");
	TH2F *histoConversionPointVsZRData=				(TH2F*)directoryGammaConvData->Get("ESD_ConversionMapping_ZR");
	TH2F *histoConversionPointVsZPhiData2=				(TH2F*)directoryGammaConvData->Get("ESD_ConversionMapping_ZPhi");
	TH2F *histoConversionPointVsRPhiData2=				(TH2F*)directoryGammaConvData->Get("ESD_ConversionMapping_RPhi");
	TH2F *histoConversionPointVsXYData=				(TH2F*)directoryGammaConvData->Get("ESD_ConversionMapping_XY");
	TH2F *histoConversionPointVsXYBPData=			(TH2F*)directoryGammaConvData->Get("ESD_ConversionMappingInnerBeampipe_XY");
	TH1D *histoConversionPointVsRData=				(TH1D*)histoConversionPointVsZRData->ProjectionY("histoConversionPointVsRData");
	ConvGammaRebinWithBinCorrection(histoConversionPointVsRData,rebinIntPlots);
	TH1D *histoConversionPointVsZData=				(TH1D*)histoConversionPointVsZRData->ProjectionX("histoConversionPointVsZData");
	ConvGammaRebinWithBinCorrection(histoConversionPointVsZData,rebinIntPlots);
	TH1D *histoConversionPointVsZCentralElectrodeData=(TH1D*)histoConversionPointVsZRData->ProjectionY("ESD_ConversionZR_CEY",histoConversionPointVsZRData->GetXaxis()->FindBin(-1),histoConversionPointVsZRData->GetXaxis()->FindBin(1));
	TH1F *histoNumberOfGoodESDTracksData=			(TH1F*)directoryGammaConvData->Get(nameNTrackHisto.Data());
	TH2F *histoConversionPointVsZRData2 =			(TH2F*) histoConversionPointVsZRData->Clone("histoConversionPointVsZRData2");
	TH2F *histoConversionPointVsXYData2 =			(TH2F*) histoConversionPointVsXYData->Clone("histoConversionPointVsXYData2");
	Double_t binWidthR = histoConversionPointVsRData->GetXaxis()->GetBinWidth(3);
	Double_t binWidthZ = histoConversionPointVsZData->GetXaxis()->GetBinWidth(3);

	TH2F *histoConversionPointVsZRPlotData =  (TH2F*) histoConversionPointVsZRData->Clone("ESD_ConversionMapping_ZRPlot");
	/*for (Int_t iZ = 1 ; iZ < histoConversionPointVsZRData->GetNbinsX(); iZ ++){
		for (Int_t iR = 1; iR < histoConversionPointVsZRData->GetNbinsY(); iR ++) {
			histoConversionPointVsZRPlotData->Fill(histoConversionPointVsZRData->GetXaxis()->GetBinCenter(iZ), histoConversionPointVsZRData->GetYaxis()->GetBinCenter(iR),  histoConversionPointVsZRData->GetBinContent(iZ,iR));
		}
	}
	*/
	TH2F *histoConversionPointVsXYPlotData =  new TH2F("ESD_ConversionMapping_XYPlot","",2400,-120.,120.,2400,-120.,120.);
	for (Int_t iX = 1 ; iX < histoConversionPointVsXYData->GetNbinsX(); iX ++){
		for (Int_t iY = 1; iY < histoConversionPointVsXYData->GetNbinsY(); iY ++) {
			Double_t radius = TMath::Sqrt(histoConversionPointVsXYData->GetXaxis()->GetBinCenter(iX)* histoConversionPointVsXYData->GetXaxis()->GetBinCenter(iX) + histoConversionPointVsXYData->GetYaxis()->GetBinCenter(iY)* histoConversionPointVsXYData->GetYaxis()->GetBinCenter(iY));
// 			if (radius > rmin) {
				histoConversionPointVsXYPlotData->Fill(histoConversionPointVsXYData->GetXaxis()->GetBinCenter(iX), histoConversionPointVsXYData->GetYaxis()->GetBinCenter(iY),  histoConversionPointVsXYData->GetBinContent(iX,iY));
// 			} 
		}
	}
	Float_t	numberGoodEventsData = 		histoEventQualityData->GetBinContent(1);
	Float_t numberGoodTriggerData = 		histoEventQualityData->GetEntries();
	Float_t numberReconstGammaData =		histoConversionPointVsRData->GetEntries();
   
	cout<< data << "    Number of events::   " << numberGoodEventsData << "    Number of triggers::   " << numberGoodTriggerData << "    Number reconstructed gammas::    "<< numberReconstGammaData <<endl;
	cout << "bu" << endl;
	Double_t meanMultiplitcityData = 				histoNumberOfGoodESDTracksData->GetMean();
	Float_t normFactorReconstData=					1./numberGoodEventsData*1./meanMultiplitcityData;
   Float_t numberReconstGammaAbove100Data =    histoConversionPointVsRData->Integral(histoConversionPointVsRData->GetXaxis()->FindBin(100),histoConversionPointVsRData->GetNbinsX());
// 	if (optEnergy.CompareTo("PbPb_2.76TeV")==0) normFactorReconstData = 1./numberReconstGammaAbove100Data;
	//Scaling reconstr.
	histoConversionPointVsZCentralElectrodeData->Scale(1/2.);
	GammaScalingHistogramm(histoVertexZData,1./numberGoodEventsData);
	GammaScalingHistogramm(histoConversionPointVsRData,normFactorReconstData);
	GammaScalingHistogramm(histoConversionPointVsZCentralElectrodeData,normFactorReconstData);
	ConvGammaRebinWithBinCorrection(histoConversionPointVsZCentralElectrodeData,2*rebinIntPlots);
	GammaScalingHistogramm(histoConversionPointVsZData,normFactorReconstData);
	GammaScalingHistogramm(histoConversionPointVsZRData2,normFactorReconstData);
	GammaScalingHistogramm(histoConversionPointVsXYData2,normFactorReconstData);
	GammaScalingHistogramm(histoConversionPointVsXYBPData,normFactorReconstData);
	GammaScalingHistogramm(histoConversionPointVsEtaData,normFactorReconstData);
	
	TH1D*	histoMappingPhiInRData[nBinsR];
	TH2F* 	histoMappingRPhiData = (TH2F*)directoryGammaConvData->Get("ESD_ConversionMapping_RPhi");
	for(Int_t iR = 0; iR < nBinsR; iR++){
		histoMappingPhiInRData[iR]=		histoMappingRPhiData->ProjectionY( Form("histoMappingPhiInRData_%i",iR), histoMappingRPhiData->GetXaxis()->FindBin(arrayRBins[iR]), histoMappingRPhiData->GetXaxis()->FindBin(arrayRBins[iR+1]));
 		ConvGammaRebinWithBinCorrection(histoMappingPhiInRData[iR],rebinPhiPlots);
		GammaScalingHistogramm(histoMappingPhiInRData[iR],normFactorReconstData*1/600*1/(TMath::Abs(histoMappingRPhiData->GetXaxis()->GetBinUpEdge(histoMappingRPhiData->GetXaxis()->FindBin(arrayRBins[iR+1]))-histoMappingRPhiData->GetXaxis()->GetBinLowEdge(histoMappingRPhiData->GetXaxis()->FindBin(arrayRBins[iR])))));
	}

	Int_t minimumBinZNormal = 0;
	TH1D*	histoMappingPhiInZData[nBinsZ];
	TH2F* 	histoMappingZPhiData = (TH2F*)directoryGammaConvData->Get("ESD_ConversionMapping_ZPhi");
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingPhiInZData[iZ]=		histoMappingZPhiData->ProjectionY( Form("histoMappingPhiInZData_%i",iZ), histoMappingZPhiData->GetXaxis()->FindBin(arrayZBins[iZ]), histoMappingZPhiData->GetXaxis()->FindBin(arrayZBins[iZ+1]));
		if (minimumBinZNormal == 0){
			if (histoMappingPhiInZData[iZ]->GetEntries() > 0){
				minimumBinZNormal = iZ;
			}
		}
 		ConvGammaRebinWithBinCorrection(histoMappingPhiInZData[iZ],rebinPhiPlots);
		GammaScalingHistogramm(histoMappingPhiInZData[iZ],normFactorReconstData*1/180*1/(TMath::Abs(histoMappingZPhiData->GetXaxis()->GetBinUpEdge(histoMappingZPhiData->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoMappingZPhiData->GetXaxis()->GetBinLowEdge(histoMappingZPhiData->GetXaxis()->FindBin(arrayZBins[iZ])))));
		
	}

	TH1D*	histoMappingRInZData[nBinsZ];
	TH2F* 	histoMappingZRData = (TH2F*)directoryGammaConvData->Get("ESD_ConversionMapping_ZR");
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingRInZData[iZ]=		histoMappingZRData->ProjectionY( Form("histoMappingRInZData_%i",iZ), histoMappingZRData->GetXaxis()->FindBin(arrayZBins[iZ]), histoMappingZRData->GetXaxis()->FindBin(arrayZBins[iZ+1]));
 		ConvGammaRebinWithBinCorrection(histoMappingRInZData[iZ],rebin);
		GammaScalingHistogramm(histoMappingRInZData[iZ],normFactorReconstData*1/(2*TMath::Pi())*1/(TMath::Abs(histoMappingZRData->GetXaxis()->GetBinUpEdge(histoMappingZRData->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoMappingZRData->GetXaxis()->GetBinLowEdge(histoMappingZRData->GetXaxis()->FindBin(arrayZBins[iZ])))));
	}

	TH1D*	histoMappingZInRData[nBinsR];
	for(Int_t iR = 0; iR < nBinsR; iR++){
		histoMappingZInRData[iR]=		histoMappingZRData->ProjectionX( Form("histoMappingZInRData_%i",iR), histoMappingZRData->GetYaxis()->FindBin(arrayRBins[iR]), histoMappingZRData->GetYaxis()->FindBin(arrayRBins[iR+1]));
		ConvGammaRebinWithBinCorrection(histoMappingZInRData[iR],rebin);
		GammaScalingHistogramm(histoMappingZInRData[iR],normFactorReconstData*1/(2*TMath::Pi())*1/(TMath::Abs(histoMappingZRData->GetYaxis()->GetBinUpEdge(histoMappingZRData->GetYaxis()->FindBin(arrayRBins[iR+1]))-histoMappingZRData->GetYaxis()->GetBinLowEdge(histoMappingZRData->GetYaxis()->FindBin(arrayRBins[iR])))));
	}
	
	TH1D*	histoMappingMidPtPhiInRData[nBinsR];
	TH2F* 	histoMappingMiRPhiData = (TH2F*)directoryGammaConvData->Get("ESD_ConversionMappingMidPt_RPhi");
	for(Int_t iR = 0; iR < nBinsR; iR++){
		histoMappingMidPtPhiInRData[iR]=		histoMappingMiRPhiData->ProjectionY( Form("histoMappingMidPtPhiInRData_%i",iR), histoMappingMiRPhiData->GetXaxis()->FindBin(arrayRBins[iR]), histoMappingMiRPhiData->GetXaxis()->FindBin(arrayRBins[iR+1]));
 		ConvGammaRebinWithBinCorrection(histoMappingMidPtPhiInRData[iR],rebinPhiPlots);
		GammaScalingHistogramm(histoMappingMidPtPhiInRData[iR],normFactorReconstData*1/600*1/(TMath::Abs(histoMappingMiRPhiData->GetXaxis()->GetBinUpEdge(histoMappingMiRPhiData->GetXaxis()->FindBin(arrayRBins[iR+1]))-histoMappingMiRPhiData->GetXaxis()->GetBinLowEdge(histoMappingMiRPhiData->GetXaxis()->FindBin(arrayRBins[iR])))));
	}

	TH1D*	histoMappingMidPtPhiInZData[nBinsZ];
	TH2F* 	histoMappingMidPtZPhiData = (TH2F*)directoryGammaConvData->Get("ESD_ConversionMappingMidPt_ZPhi");
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingMidPtPhiInZData[iZ]=		histoMappingMidPtZPhiData->ProjectionY( Form("histoMappingMidPtPhiInZData_%i",iZ), histoMappingMidPtZPhiData->GetXaxis()->FindBin(arrayZBins[iZ]), histoMappingMidPtZPhiData->GetXaxis()->FindBin(arrayZBins[iZ+1]));
		ConvGammaRebinWithBinCorrection(histoMappingMidPtPhiInZData[iZ],rebinPhiPlots);
		GammaScalingHistogramm(histoMappingMidPtPhiInZData[iZ],normFactorReconstData*1/180*1/(TMath::Abs(histoMappingMidPtZPhiData->GetXaxis()->GetBinUpEdge(histoMappingMidPtZPhiData->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoMappingMidPtZPhiData->GetXaxis()->GetBinLowEdge(histoMappingMidPtZPhiData->GetXaxis()->FindBin(arrayZBins[iZ])))));
		
	}

	TH1D*	histoMappingMidPtRInZData[nBinsZ];
	TH2F* histoMappingMidPtZRData = (TH2F*)directoryGammaConvData->Get("ESD_ConversionMappingMidPt_ZR");
	TH1D * histoConversionPointVsRMidPtData=	(TH1D*)histoMappingMidPtZRData->ProjectionY("histoConversionPointVsRMidPtData");
	ConvGammaRebinWithBinCorrection(histoConversionPointVsRMidPtData,rebinIntPlots);
	GammaScalingHistogramm(histoConversionPointVsRMidPtData, normFactorReconstData);
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingMidPtRInZData[iZ]=		histoMappingMidPtZRData->ProjectionY( Form("histoMappingMidPtRInZData_%i",iZ), histoMappingMidPtZRData->GetXaxis()->FindBin(arrayZBins[iZ]), histoMappingMidPtZRData->GetXaxis()->FindBin(arrayZBins[iZ+1]));
 		ConvGammaRebinWithBinCorrection(histoMappingMidPtRInZData[iZ],rebin);
		GammaScalingHistogramm(histoMappingMidPtRInZData[iZ],normFactorReconstData*1/(2*TMath::Pi())*1/(TMath::Abs(histoMappingMidPtZRData->GetXaxis()->GetBinUpEdge(histoMappingMidPtZRData->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoMappingMidPtZRData->GetXaxis()->GetBinLowEdge(histoMappingMidPtZRData->GetXaxis()->FindBin(arrayZBins[iZ])))));
	}

	TH1D*	histoMappingMidPtZInRData[nBinsR];
	for(Int_t iR = 0; iR < nBinsR; iR++){
		histoMappingMidPtZInRData[iR]=		histoMappingMidPtZRData->ProjectionX( Form("histoMappingMidPtZInRData_%i",iR), histoMappingMidPtZRData->GetYaxis()->FindBin(arrayRBins[iR]), histoMappingMidPtZRData->GetYaxis()->FindBin(arrayRBins[iR+1]));
		ConvGammaRebinWithBinCorrection(histoMappingMidPtZInRData[iR],rebin);
		GammaScalingHistogramm(histoMappingMidPtZInRData[iR],normFactorReconstData*1/(2*TMath::Pi())*1/(TMath::Abs(histoMappingMidPtZRData->GetYaxis()->GetBinUpEdge(histoMappingMidPtZRData->GetYaxis()->FindBin(arrayRBins[iR+1]))-histoMappingMidPtZRData->GetYaxis()->GetBinLowEdge(histoMappingMidPtZRData->GetYaxis()->FindBin(arrayRBins[iR])))));
	}
	
	
	Int_t minimumBinZSPD = 0;
	TH1D*		histoMappingSPDPhiInZData[nBinsZ];
	TH2F* 	histoMappingSPDZPhiData = (TH2F*)directoryGammaConvData->Get("ESD_ConversionMapping_SPD_ZPhi");
	TString nameHistoSPDPhiInZData;
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingSPDPhiInZData[iZ]=		histoMappingSPDZPhiData->ProjectionY( Form("histoMappingPhiInZ_SPD_Data_%i",iZ), histoMappingSPDZPhiData->GetXaxis()->FindBin(arrayZBins[iZ]), histoMappingSPDZPhiData->GetXaxis()->FindBin(arrayZBins[iZ+1]));
		if (minimumBinZSPD == 0){
			if (histoMappingSPDPhiInZData[iZ]->GetEntries() > 0){
				minimumBinZSPD = iZ;
			}
		}
 		ConvGammaRebinWithBinCorrection(histoMappingSPDPhiInZData[iZ],rebinPhiPlots);
		GammaScalingHistogramm(histoMappingSPDPhiInZData[iZ],normFactorReconstData*1/9.5*1/(TMath::Abs(histoMappingSPDZPhiData->GetXaxis()->GetBinUpEdge(histoMappingSPDZPhiData->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoMappingSPDZPhiData->GetXaxis()->GetBinLowEdge(histoMappingSPDZPhiData->GetXaxis()->FindBin(arrayZBins[iZ])))));
		
	}

	Int_t minimumBinZSPDTh = 0;
	TH1D*		histoMappingSPDThPhiInZData[nBinsZ];
	TH2F* 	histoMappingSPDThZPhiData = (TH2F*)directoryGammaConvData->Get("ESD_ConversionMapping_SPDTh_ZPhi");
	TString nameHistoSPDThPhiInZData;
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingSPDThPhiInZData[iZ]=		histoMappingSPDThZPhiData->ProjectionY( Form("histoMappingPhiInZ_SPDTh_Data_%i",iZ), histoMappingSPDThZPhiData->GetXaxis()->FindBin(arrayZBins[iZ]), histoMappingSPDThZPhiData->GetXaxis()->FindBin(arrayZBins[iZ+1]));
		if (minimumBinZSPDTh == 0){
			if (histoMappingSPDThPhiInZData[iZ]->GetEntries() > 0){
				minimumBinZSPDTh = iZ;
			}
		}
 		ConvGammaRebinWithBinCorrection(histoMappingSPDThPhiInZData[iZ],rebinPhiPlots);
		GammaScalingHistogramm(histoMappingSPDThPhiInZData[iZ],normFactorReconstData*1/3.5*1/(TMath::Abs(histoMappingSPDThZPhiData->GetXaxis()->GetBinUpEdge(histoMappingSPDThZPhiData->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoMappingSPDThZPhiData->GetXaxis()->GetBinLowEdge(histoMappingSPDThZPhiData->GetXaxis()->FindBin(arrayZBins[iZ])))));
		
	}
	
	Int_t minimumBinZSDD = 0;
	TH1D*		histoMappingSDDPhiInZData[nBinsZ];
	TH2F* 	histoMappingSDDZPhiData = (TH2F*)directoryGammaConvData->Get("ESD_ConversionMapping_SDD_ZPhi");
	TString nameHistoSDDPhiInZData;
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingSDDPhiInZData[iZ]=		histoMappingSDDZPhiData->ProjectionY( Form("histoMappingPhiInZ_SDD_Data_%i",iZ), histoMappingSDDZPhiData->GetXaxis()->FindBin(arrayZBins[iZ]), histoMappingSDDZPhiData->GetXaxis()->FindBin(arrayZBins[iZ+1]));
		if (minimumBinZSDD == 0){
			if (histoMappingSDDPhiInZData[iZ]->GetEntries() > 0){
				minimumBinZSDD = iZ;
			}
		}
 		ConvGammaRebinWithBinCorrection(histoMappingSDDPhiInZData[iZ],rebinPhiPlots);
		GammaScalingHistogramm(histoMappingSDDPhiInZData[iZ],normFactorReconstData*1/22.*1/(TMath::Abs(histoMappingSDDZPhiData->GetXaxis()->GetBinUpEdge(histoMappingSDDZPhiData->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoMappingSDDZPhiData->GetXaxis()->GetBinLowEdge(histoMappingSDDZPhiData->GetXaxis()->FindBin(arrayZBins[iZ])))));
		
	}

	Int_t minimumBinZSSD = 0;
	TH1D*		histoMappingSSDPhiInZData[nBinsZ];
	TH2F* 	histoMappingSSDZPhiData = (TH2F*)directoryGammaConvData->Get("ESD_ConversionMapping_SSD_ZPhi");
	TString nameHistoSSDPhiInZData;
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingSSDPhiInZData[iZ]=		histoMappingSSDZPhiData->ProjectionY( Form("histoMappingPhiInZ_SSD_Data_%i",iZ), histoMappingSSDZPhiData->GetXaxis()->FindBin(arrayZBins[iZ]), histoMappingSSDZPhiData->GetXaxis()->FindBin(arrayZBins[iZ+1]));
		if (minimumBinZSSD == 0){
			if (histoMappingSSDPhiInZData[iZ]->GetEntries() > 0){
				minimumBinZSSD = iZ;
			}
		}
 		ConvGammaRebinWithBinCorrection(histoMappingSSDPhiInZData[iZ],rebinPhiPlots);
		GammaScalingHistogramm(histoMappingSSDPhiInZData[iZ],normFactorReconstData*1/20.*1/(TMath::Abs(histoMappingSSDZPhiData->GetXaxis()->GetBinUpEdge(histoMappingSSDZPhiData->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoMappingSSDZPhiData->GetXaxis()->GetBinLowEdge(histoMappingSSDZPhiData->GetXaxis()->FindBin(arrayZBins[iZ])))));
		
	}

	
	Int_t minimumBinZITSTPC = 0;
	TH1D*		histoMappingITSTPCPhiInZData[nBinsZ];
	TH2F* 	histoMappingITSTPCPhiZData = (TH2F*)directoryGammaConvData->Get("ESD_ConversionMapping_ITSTPC_ZPhi");
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingITSTPCPhiInZData[iZ]=		histoMappingITSTPCPhiZData->ProjectionY( Form("histoMappingPhiInZ_ITSTPC_Data_%i",iZ), histoMappingITSTPCPhiZData->GetXaxis()->FindBin(arrayZBins[iZ]), histoMappingITSTPCPhiZData->GetXaxis()->FindBin(arrayZBins[iZ+1]));
		if (minimumBinZITSTPC == 0){
			if (histoMappingITSTPCPhiInZData[iZ]->GetEntries() > 0){
				minimumBinZITSTPC = iZ;
			}
		}
 		ConvGammaRebinWithBinCorrection(histoMappingITSTPCPhiInZData[iZ],rebinPhiPlots);
		GammaScalingHistogramm(histoMappingITSTPCPhiInZData[iZ],normFactorReconstData*1/17.5*1/(TMath::Abs(histoMappingITSTPCPhiZData->GetXaxis()->GetBinUpEdge(histoMappingITSTPCPhiZData->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoMappingITSTPCPhiZData->GetXaxis()->GetBinLowEdge(histoMappingITSTPCPhiZData->GetXaxis()->FindBin(arrayZBins[iZ])))));
	}
	
	Int_t nRowsZProjections = 0;
	Int_t nColumnsZProjections = 0;
	Int_t nRowsZSPD = 0;
	Int_t nColumnsZSPD = 0;
	Int_t nRowsZSPDTh = 0;
	Int_t nColumnsZSPDTh = 0;
	Int_t nRowsZITSTPC = 0;
	Int_t nColumnsZITSTPC = 0;
	
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

	if ( (nBinsZ - 2*minimumBinZSPDTh ) < 15){
		nColumnsZSPDTh = 2;
	} else {
		nColumnsZSPDTh = 3;
	}
	if ( ((nBinsZ - 2*minimumBinZSPDTh)%nColumnsZSPDTh) == 0){
		nRowsZSPDTh = (nBinsZ - 2*minimumBinZSPDTh)/nColumnsZSPDTh ;
	} else {
		nRowsZSPDTh = (nBinsZ - 2*minimumBinZSPDTh)/nColumnsZSPDTh + 1;
	}
	cout << "SPDTh minZBin: " << minimumBinZSPDTh << "\t col: " << nColumnsZSPDTh  << "\t row: " << nRowsZSPDTh << endl;
	
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

	// ------------------------------------- second file-----------------------------------------------------------
	TH1D* 	histoMappingPhiInRMonteCarlo[nBinsR];
	TH1D* 	histoMappingPhiInRMCConv[nBinsR];
	TH1D* 	histoMappingZInRMonteCarlo[nBinsR];
	TH1D* 	histoMappingZInRMCConv[nBinsR];
	TH1D* 	histoMappingPhiInZMonteCarlo[nBinsZ];
	TH1D* 	histoMappingRInZMonteCarlo[nBinsZ];
	TH1D* 	histoMappingMidPtPhiInRMonteCarlo[nBinsR];
	TH1D* 	histoMappingMidPtZInRMonteCarlo[nBinsR];
	TH1D* 	histoMappingMidPtPhiInZMonteCarlo[nBinsZ];
	TH1D* 	histoMappingMidPtRInZMonteCarlo[nBinsZ];		
	TH1D* 	histoMappingSPDPhiInZMonteCarlo[nBinsZ];
	TH1D* 	histoMappingSPDThPhiInZMonteCarlo[nBinsZ];
	TH1D* 	histoMappingITSTPCPhiInZMonteCarlo[nBinsZ];
	
	TDirectory* directoryGammaConvMonteCarlo = 	new TDirectory(); // definition of first folder / list    

	TString 	nameHistoRatioPhiInRMonteCarlo;
	TH1F*		histoMappingRatioPhiInR[nBinsR];
	TString 	nameHistoRatioZInRMonteCarlo;
	TH1F* 	histoMappingRatioZInR[nBinsR];
	TString	nameHistoRatioPhiInZMonteCarlo;
	TH1F* 	histoMappingRatioPhiInZ[nBinsZ];
	TString 	nameHistoSPDPhiInZMonteCarlo;
	TString 	nameHistoSDDPhiInZMonteCarlo;
	TString 	nameHistoSSDPhiInZMonteCarlo;
	TString 	nameHistoSPDThPhiInZMonteCarlo;
	TString	nameHistoRatioSPDPhiInZMonteCarlo;
	TString	nameHistoRatioSSDPhiInZMonteCarlo;
	TString	nameHistoRatioSDDPhiInZMonteCarlo;
	TString	nameHistoRatioSPDThPhiInZMonteCarlo;
	TH1F * 	histoMappingRatioSPDPhiInZ[nBinsZ];
	TH1F * 	histoMappingRatioSSDPhiInZ[nBinsZ];
	TH1F * 	histoMappingRatioSDDPhiInZ[nBinsZ];
	TH1F * 	histoMappingRatioSPDThPhiInZ[nBinsZ];
	TString 	nameHistoITSTPCPhiInZMonteCarlo;
	TString	nameHistoRatioITSTPCPhiInZMonteCarlo;
	TH1F * 	histoMappingRatioITSTPCPhiInZ[nBinsZ];
	TString 	nameHistoRatioRInZMonteCarlo;
	TH1F* 	histoMappingRatioRInZ[nBinsZ];
	
	TFile* montecarlo = new TFile(MCfile); 
	
	directoryGammaConvMonteCarlo = 			(TDirectory*)montecarlo->Get(nameGammaDirectory); 
	
	TH1F* 	histoEventQualityMonteCarlo=					(TH1F*)directoryGammaConvMonteCarlo->Get("ESD_EventQuality");
	TH1F* histoVertexZMonteCarlo=					 	(TH1F*)directoryGammaConvMonteCarlo->Get("Z_Vertex_distribution");
	TH1F* 	histoConversionPointVsEtaMonteCarlo=		(TH1F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMapping_Eta");
	TH2F* 	histoConversionPointVsZRMonteCarlo=			(TH2F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMapping_ZR");
	TH2F* 	histoConversionPointVsXYMonteCarlo=			(TH2F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMapping_XY");
	TH2F* 	histoConversionPointVsXYBPMonteCarlo=		(TH2F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMappingInnerBeampipe_XY");
	TH1D* 	histoConversionPointVsRMonteCarlo=			(TH1D*)histoConversionPointVsZRMonteCarlo->ProjectionY("histoConversionPointVsRMonteCarlo");
	ConvGammaRebinWithBinCorrection(histoConversionPointVsRMonteCarlo,rebinIntPlots);
	TH1D* 	histoConversionPointVsZMonteCarlo=			(TH1D*)histoConversionPointVsZRMonteCarlo->ProjectionX("histoConversionPointVsZMonteCarlo");
	ConvGammaRebinWithBinCorrection(histoConversionPointVsZMonteCarlo,rebinIntPlots);

// 	TH2F* 	histoMCConversionPointVsZRMonteCarlo=		(TH2F*)directoryGammaConvMonteCarlo->Get("MC_ConversionMapping_ZR");
// 	TH2F* 	histoMCConversionPointVsXYMonteCarlo=		(TH2F*)directoryGammaConvMonteCarlo->Get("MC_ConversionMapping_XY");
// 	TH1D* 	histoMCConversionVsR=			(TH1D*)histoMCConversionPointVsZRMonteCarlo->ProjectionY("MC_Conversion_R");
// 	Double_t nGammasMC = histoMCConversionVsR->GetEntries();
// 	ConvGammaRebinWithBinCorrection(histoMCConversionVsR,rebinIntPlots);
// 	TH1D* 	histoMCConversionVsZ=			(TH1D*)histoMCConversionPointVsZRMonteCarlo->ProjectionX("MC_Conversion_Z");
// 	ConvGammaRebinWithBinCorrection(histoMCConversionVsZ,rebinIntPlots);
		
	TH1D* histoDeltaPt =			 	new TH1D("histoDeltaPt","",32,floatPtValues7TeVPi0);
	for(Int_t iPt=1;iPt<32+1;iPt++){
		histoDeltaPt->SetBinContent(iPt,floatPtValues7TeVPi0[iPt]-floatPtValues7TeVPi0[iPt-1]);
		histoDeltaPt->SetBinError(iPt,0);
	}


	
	TH1F* 	histoTrueConversionPointVsEta=		(TH1F*)directoryGammaConvMonteCarlo->Get("ESD_TrueConversionMapping_Eta");
	TH2F* 	histoTrueConversionPointVsZR=			(TH2F*)directoryGammaConvMonteCarlo->Get("ESD_TrueConversionMapping_ZR");	
	TH1D* 	histoTrueConversionPointVSR=			(TH1D*)histoTrueConversionPointVsZR->ProjectionY("ESD_TrueConversion_R");
	ConvGammaRebinWithBinCorrection(histoTrueConversionPointVSR,rebinIntPlots);
	TH1D* 	histoTrueConversionPointVSZ=			(TH1D*)histoTrueConversionPointVsZR->ProjectionX("ESD_TrueConversion_Z");
	ConvGammaRebinWithBinCorrection(histoTrueConversionPointVSZ,rebinIntPlots);

	TH1F* 	histoTrueDalitzPi0VsEta=		(TH1F*)directoryGammaConvMonteCarlo->Get("ESD_TrueDalPi0Mapping_Eta");
	TH2F* 	histoTrueDalitzPi0VsZR=			(TH2F*)directoryGammaConvMonteCarlo->Get("ESD_TrueDalPi0Mapping_ZR");	
	TH1D* 	histoTrueDalitzPi0VsR=			(TH1D*)histoTrueDalitzPi0VsZR->ProjectionY("ESD_TrueDalitzPi0_R");
	ConvGammaRebinWithBinCorrection(histoTrueDalitzPi0VsR,rebinIntPlots);
	TH1D* 	histoTrueDalitzPi0VsZ=			(TH1D*)histoTrueDalitzPi0VsZR->ProjectionX("ESD_TrueDalitzPi0_Z");
	ConvGammaRebinWithBinCorrection(histoTrueDalitzPi0VsZ,rebinIntPlots);

	TH1F* 	histoTrueDalitzEtaVsEta=		(TH1F*)directoryGammaConvMonteCarlo->Get("ESD_TrueDalEtaMapping_Eta");
	TH2F* 	histoTrueDalitzEtaVsZR=			(TH2F*)directoryGammaConvMonteCarlo->Get("ESD_TrueDalEtaMapping_ZR");	
	TH1D* 	histoTrueDalitzEtaVsR=			(TH1D*)histoTrueDalitzEtaVsZR->ProjectionY("ESD_TrueDalitzEta_R");
	ConvGammaRebinWithBinCorrection(histoTrueDalitzEtaVsR,rebinIntPlots);
	TH1D* 	histoTrueDalitzEtaVsZ=			(TH1D*)histoTrueDalitzEtaVsZR->ProjectionX("ESD_TrueDalitzEta_Z");
	ConvGammaRebinWithBinCorrection(histoTrueDalitzEtaVsZ,rebinIntPlots);

	
	TH1F* 	histoTrueCombinatoricsVsEta=		(TH1F*)directoryGammaConvMonteCarlo->Get("ESD_TrueCombMapping_Eta");
	TH2F* 	histoTrueCombinatoricsVsZR=		(TH2F*)directoryGammaConvMonteCarlo->Get("ESD_TrueCombMapping_ZR");	
	TH1D* 	histoTrueCombinatoricsVsR=			(TH1D*)histoTrueCombinatoricsVsZR->ProjectionY("ESD_TrueConvCombinatorial_R");
	ConvGammaRebinWithBinCorrection(histoTrueCombinatoricsVsR,rebinIntPlots);
	TH1D* 	histoTrueCombinatoricsVsZ=			(TH1D*)histoTrueCombinatoricsVsZR->ProjectionX("ESD_TrueConvCombinatorial_Z");
	ConvGammaRebinWithBinCorrection(histoTrueCombinatoricsVsZ,rebinIntPlots);

	TH1F* 	histoTrueHadVsEta=		(TH1F*)directoryGammaConvMonteCarlo->Get("ESD_TrueHadMapping_Eta");
	TH2F* 	histoTrueHadVsZR=					(TH2F*)directoryGammaConvMonteCarlo->Get("ESD_TrueHadMapping_ZR");	
	TH1D* 	histoTrueHadVsR=			(TH1D*)histoTrueHadVsZR->ProjectionY("ESD_TrueConvHadronicBck_R");
	ConvGammaRebinWithBinCorrection(histoTrueHadVsR,rebinIntPlots);
	TH1D* 	histoTrueHadVsZ=			(TH1D*)histoTrueHadVsZR->ProjectionX("ESD_TrueConvHadronicBck_Z");
	ConvGammaRebinWithBinCorrection(histoTrueHadVsZ,rebinIntPlots);

	TH1F* 	histoTrueConversionPointPrimVsEta=		(TH1F*)directoryGammaConvMonteCarlo->Get("ESD_TruePrimMapping_Eta");
	TH2F* 	histoTrueConversionPointPrimVsZR=		(TH2F*)directoryGammaConvMonteCarlo->Get("ESD_TruePrimMapping_ZR");	
	TH1D* 	histoTrueConversionPointPrimVSR=			(TH1D*)histoTrueConversionPointPrimVsZR->ProjectionY("ESD_TrueConvPrimaryGamma_R");
	ConvGammaRebinWithBinCorrection(histoTrueConversionPointPrimVSR,rebinIntPlots);
	TH1D* 	histoTrueConversionPointPrimVSZ=			(TH1D*)histoTrueConversionPointPrimVsZR->ProjectionX("ESD_TrueConvPrimaryGamma_Z");
	ConvGammaRebinWithBinCorrection(histoTrueConversionPointPrimVSZ,rebinIntPlots);

	TH1F* 	histoTrueConversionPointSecVSEta = (TH1F*) histoTrueConversionPointVsEta->Clone("histoTrueConversionPointSecVSEta");
	histoTrueConversionPointSecVSEta->Sumw2();
	histoTrueConversionPointSecVSEta->Add(histoTrueConversionPointPrimVsEta,-1.);
	TH1F* 	histoTrueConversionPointSecVSR = (TH1F*) histoTrueConversionPointVSR->Clone("histoTrueConversionPointSecVSR");
	histoTrueConversionPointSecVSR->Sumw2();
	histoTrueConversionPointSecVSR->Add(histoTrueConversionPointPrimVSR,-1.);
	TH1F* 	histoTrueConversionPointSecVSZ = (TH1F*) histoTrueConversionPointVSZ->Clone("histoTrueConversionPointSecVSZ");
	histoTrueConversionPointSecVSZ->Sumw2();
	histoTrueConversionPointSecVSZ->Add(histoTrueConversionPointPrimVSZ,-1.);
	
	TH1D* 	histoConversionPointVsZCentralElectrodeMonteCarlo=	(TH1D*)histoConversionPointVsZRMonteCarlo->ProjectionY("ESD_ConversionZR_CEY",histoConversionPointVsZRMonteCarlo->GetXaxis()->FindBin(-1),histoConversionPointVsZRMonteCarlo->GetXaxis()->FindBin(1));
// 	TH1D* 	histoMCConversionPointVsZCentralElectrodeMonteCarlo=	(TH1D*)histoMCConversionPointVsZRMonteCarlo->ProjectionY("MC_ConversionZR_CEY",histoMCConversionPointVsZRMonteCarlo->GetXaxis()->FindBin(-1),histoMCConversionPointVsZRMonteCarlo->GetXaxis()->FindBin(1));
	
	TH1F* 	histoNumberOfGoodESDTracksMonteCarlo=				(TH1F*)directoryGammaConvMonteCarlo->Get(nameNTrackHisto.Data());                       
		
	TH2F* 	histoConversionPointVsZRMonteCarlo2 =	(TH2F*) histoConversionPointVsZRMonteCarlo->Clone("histoConversionPointVsZRMonteCarlo2");
	TH2F* 	histoConversionPointVsXYMonteCarlo2 =	(TH2F*) histoConversionPointVsXYMonteCarlo->Clone("histoConversionPointVsXYMonteCarlo2");
		
	Float_t numberGoodEventsMonteCarlo = 					histoEventQualityMonteCarlo->GetBinContent(1);
	Float_t numberGoodTriggerMonteCarlo = 					histoEventQualityMonteCarlo->GetEntries();
	Float_t numberReconstGammaMonteCarlo = 					histoConversionPointVsRMonteCarlo->GetEntries();
	cout<< MCfile << "    Number of events::   " << numberGoodEventsMonteCarlo << "    Number of triggers::   " << numberGoodTriggerMonteCarlo << "    Number reconstructed gammas::    "<< numberReconstGammaMonteCarlo <<endl;
		
	Double_t meanMultiplitcityMonteCarlo = 					histoNumberOfGoodESDTracksMonteCarlo->GetMean();
	Float_t normFactorReconstMonteCarlo=						1./numberGoodEventsMonteCarlo * 1./meanMultiplitcityMonteCarlo;
	Float_t numberReconstGammaAbove100MonteCarlo =    histoConversionPointVsRMonteCarlo->Integral(histoConversionPointVsRMonteCarlo->GetXaxis()->FindBin(100),histoConversionPointVsRMonteCarlo->GetNbinsX());
//    if (optEnergy.CompareTo("PbPb_2.76TeV")==0) normFactorReconstMonteCarlo = 1./numberReconstGammaAbove100MonteCarlo;
   
	
// 	histoMCConversionPointVsZRMonteCarlo->Scale(1./histoMCConversionPointVsZRMonteCarlo->GetEntries());
// 	histoMCConversionPointVsXYMonteCarlo->Scale(1./histoMCConversionPointVsXYMonteCarlo->GetEntries());
	histoConversionPointVsZCentralElectrodeMonteCarlo->Scale(1/2.);
	ConvGammaRebinWithBinCorrection(histoConversionPointVsZCentralElectrodeMonteCarlo,2*rebinIntPlots);
// 	histoMCConversionPointVsZCentralElectrodeMonteCarlo->Scale(1/2.);
// 	ConvGammaRebinWithBinCorrection(histoMCConversionPointVsZCentralElectrodeMonteCarlo,2*rebinIntPlots);
	GammaScalingHistogramm(histoConversionPointVsRMonteCarlo,normFactorReconstMonteCarlo);
	
	GammaScalingHistogramm(histoVertexZMonteCarlo,1./numberGoodEventsMonteCarlo);
// 	GammaScalingHistogramm(histoMCConversionVsR,1./nGammasMC * 1./meanMultiplitcityMonteCarlo);
	GammaScalingHistogramm(histoTrueConversionPointVSR,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoConversionPointVsXYBPMonteCarlo,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoTrueDalitzPi0VsR,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoTrueDalitzEtaVsR,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoTrueCombinatoricsVsR,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoTrueHadVsR,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoConversionPointVsZCentralElectrodeMonteCarlo,normFactorReconstMonteCarlo);
// 	GammaScalingHistogramm(histoMCConversionPointVsZCentralElectrodeMonteCarlo,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoConversionPointVsZMonteCarlo,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoConversionPointVsZRMonteCarlo2,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoConversionPointVsXYMonteCarlo2,normFactorReconstMonteCarlo);
// 	GammaScalingHistogramm(histoMCConversionVsR,1./histoMCConversionVsR->GetEntries());
// 	GammaScalingHistogramm(histoMCConversionVsZ,1./histoMCConversionVsZ->GetEntries());
	GammaScalingHistogramm(histoTrueConversionPointVSZ,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoTrueDalitzPi0VsZ,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoTrueDalitzEtaVsZ,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoTrueCombinatoricsVsZ,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoTrueConversionPointPrimVSR,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoTrueConversionPointPrimVSZ,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoTrueConversionPointSecVSR,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoTrueConversionPointSecVSZ,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoTrueHadVsZ,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoConversionPointVsEtaMonteCarlo,normFactorReconstMonteCarlo);	
	GammaScalingHistogramm(histoTrueDalitzPi0VsEta,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoTrueDalitzEtaVsEta,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoTrueCombinatoricsVsEta,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoTrueHadVsEta,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoTrueConversionPointVsEta,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoTrueConversionPointPrimVsEta,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoTrueConversionPointSecVSEta,normFactorReconstMonteCarlo);

	TH2F* 	histoMappingRPhiMonteCarlo = (TH2F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMapping_RPhi");
	for(Int_t iR = 0; iR < nBinsR; iR++){
		histoMappingPhiInRMonteCarlo[iR]=		histoMappingRPhiMonteCarlo->ProjectionY( Form("histoMappingPhiInRMonteCarlo_%i",iR), histoMappingRPhiMonteCarlo->GetXaxis()->FindBin(arrayRBins[iR]), histoMappingRPhiMonteCarlo->GetXaxis()->FindBin(arrayRBins[iR+1]));
 		ConvGammaRebinWithBinCorrection(histoMappingPhiInRMonteCarlo[iR],rebinPhiPlots);
		GammaScalingHistogramm(histoMappingPhiInRMonteCarlo[iR],normFactorReconstMonteCarlo*1/600*1/(TMath::Abs(histoMappingRPhiMonteCarlo->GetXaxis()->GetBinUpEdge(histoMappingRPhiMonteCarlo->GetXaxis()->FindBin(arrayRBins[iR+1]))-histoMappingRPhiMonteCarlo->GetXaxis()->GetBinLowEdge(histoMappingRPhiMonteCarlo->GetXaxis()->FindBin(arrayRBins[iR])))));
		nameHistoRatioPhiInRMonteCarlo=Form("histoMappingRatioPhiInR_%02d",iR);
		histoMappingRatioPhiInR[iR]= 			(TH1F*)histoMappingPhiInRData[iR]->Clone();
		histoMappingRatioPhiInR[iR]->SetName(nameHistoRatioPhiInRMonteCarlo); 
		histoMappingRatioPhiInR[iR]->Divide(histoMappingRatioPhiInR[iR],histoMappingPhiInRMonteCarlo[iR]);
		
	}

// 	TH2F* 	histoMappingRPhiMCConv = (TH2F*)directoryGammaConvMonteCarlo->Get("MC_ConversionMapping_RPhi");
// 	for(Int_t iR = 0; iR < nBinsR; iR++){
// 		histoMappingPhiInRMCConv[iR]=		histoMappingRPhiMCConv->ProjectionY( Form("histoMappingPhiInRMCConv_%i",iR), histoMappingRPhiMCConv->GetXaxis()->FindBin(arrayRBins[iR]), histoMappingRPhiMCConv->GetXaxis()->FindBin(arrayRBins[iR+1]));
//  		ConvGammaRebinWithBinCorrection(histoMappingPhiInRMCConv[iR],rebinPhiPlots);
// 		GammaScalingHistogramm(histoMappingPhiInRMCConv[iR],normFactorReconstMonteCarlo*1/600*1/(TMath::Abs(histoMappingRPhiMCConv->GetXaxis()->GetBinUpEdge(histoMappingRPhiMCConv->GetXaxis()->FindBin(arrayRBins[iR+1]))-histoMappingRPhiMCConv->GetXaxis()->GetBinLowEdge(histoMappingRPhiMCConv->GetXaxis()->FindBin(arrayRBins[iR])))));
// 	}

	
	TH2F* 	histoMappingZPhiMonteCarlo = (TH2F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMapping_ZPhi");
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingPhiInZMonteCarlo[iZ]=		histoMappingZPhiMonteCarlo->ProjectionY( Form("histoMappingPhiInZMonteCarlo_%i",iZ), histoMappingZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ]), histoMappingZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]));
		ConvGammaRebinWithBinCorrection(histoMappingPhiInZMonteCarlo[iZ],rebinPhiPlots);
		GammaScalingHistogramm(histoMappingPhiInZMonteCarlo[iZ],normFactorReconstMonteCarlo*1/180*1/(TMath::Abs(histoMappingZPhiMonteCarlo->GetXaxis()->GetBinUpEdge(histoMappingZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoMappingZPhiMonteCarlo->GetXaxis()->GetBinLowEdge(histoMappingZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ])))));
		nameHistoRatioPhiInZMonteCarlo=Form("histoMappingRatioPhiInZ_%02d",iZ);
		histoMappingRatioPhiInZ[iZ]= 			(TH1F*)histoMappingPhiInZData[iZ]->Clone();
		histoMappingRatioPhiInZ[iZ]->SetName(nameHistoRatioPhiInZMonteCarlo); 
		histoMappingRatioPhiInZ[iZ]->Divide(histoMappingRatioPhiInZ[iZ],histoMappingPhiInZMonteCarlo[iZ]);
		
	}

	TH2F* 	histoMappingZRMonteCarlo = (TH2F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMapping_ZR");
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingRInZMonteCarlo[iZ]=		histoMappingZRMonteCarlo->ProjectionY( Form("histoMappingRInZMonteCarlo_%i",iZ), histoMappingZRMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ]), histoMappingZRMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]));
 		ConvGammaRebinWithBinCorrection(histoMappingRInZMonteCarlo[iZ],rebin);
		GammaScalingHistogramm(histoMappingRInZMonteCarlo[iZ],normFactorReconstMonteCarlo*1/(2*TMath::Pi())*1/(TMath::Abs(histoMappingZRMonteCarlo->GetXaxis()->GetBinUpEdge(histoMappingZRMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoMappingZRMonteCarlo->GetXaxis()->GetBinLowEdge(histoMappingZRMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ])))));
		nameHistoRatioRInZMonteCarlo=Form("histoMappingRatioRInZ_%02d",iZ);
		histoMappingRatioRInZ[iZ]= 			(TH1F*)histoMappingRInZData[iZ]->Clone();
		histoMappingRatioRInZ[iZ]->SetName(nameHistoRatioRInZMonteCarlo); 
		histoMappingRatioRInZ[iZ]->Divide(histoMappingRatioRInZ[iZ],histoMappingRInZMonteCarlo[iZ]);
	}
	for(Int_t iR = 0; iR < nBinsR; iR++){
		histoMappingZInRMonteCarlo[iR]=		histoMappingZRMonteCarlo->ProjectionX( Form("histoMappingZInRMonteCarlo_%i",iR), histoMappingZRMonteCarlo->GetYaxis()->FindBin(arrayRBins[iR]), histoMappingZRMonteCarlo->GetYaxis()->FindBin(arrayRBins[iR+1]));
		ConvGammaRebinWithBinCorrection(histoMappingZInRMonteCarlo[iR],rebin);
		GammaScalingHistogramm(histoMappingZInRMonteCarlo[iR],normFactorReconstMonteCarlo*1/(2*TMath::Pi())*1/(TMath::Abs(histoMappingZRMonteCarlo->GetYaxis()->GetBinUpEdge(histoMappingZRMonteCarlo->GetYaxis()->FindBin(arrayRBins[iR+1]))-histoMappingZRMonteCarlo->GetYaxis()->GetBinLowEdge(histoMappingZRMonteCarlo->GetYaxis()->FindBin(arrayRBins[iR])))));
		nameHistoRatioZInRMonteCarlo=Form("histoMappingRatioZInR_%02d",iR);
		histoMappingRatioZInR[iR]= (TH1F*)histoMappingZInRData[iR]->Clone();
		histoMappingRatioZInR[iR]->SetName(nameHistoRatioZInRMonteCarlo); 
		histoMappingRatioZInR[iR]->Divide(histoMappingRatioZInR[iR],histoMappingZInRMonteCarlo[iR]);
	}

// 	TH2F* 	histoMappingZRMCConv = (TH2F*)directoryGammaConvMonteCarlo->Get("MC_ConversionMapping_ZR");
// 	for(Int_t iR = 0; iR < nBinsR; iR++){
// 		histoMappingZInRMCConv[iR]=		histoMappingZRMCConv->ProjectionX( Form("histoMappingZInRMCConv_%i",iR), histoMappingZRMCConv->GetYaxis()->FindBin(arrayRBins[iR]), histoMappingZRMCConv->GetYaxis()->FindBin(arrayRBins[iR+1]));
// 		ConvGammaRebinWithBinCorrection(histoMappingZInRMCConv[iR],rebin);
// 		GammaScalingHistogramm(histoMappingZInRMCConv[iR],normFactorReconstMonteCarlo*1/(2*TMath::Pi())*1/(TMath::Abs(histoMappingZRMCConv->GetYaxis()->GetBinUpEdge(histoMappingZRMCConv->GetYaxis()->FindBin(arrayRBins[iR+1]))-histoMappingZRMCConv->GetYaxis()->GetBinLowEdge(histoMappingZRMCConv->GetYaxis()->FindBin(arrayRBins[iR])))));
// 	}


	TH2F* 	histoMappingMidPtRPhiMonteCarlo = (TH2F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMappingMidPt_RPhi");
	for(Int_t iR = 0; iR < nBinsR; iR++){
		histoMappingMidPtPhiInRMonteCarlo[iR]=		histoMappingMidPtRPhiMonteCarlo->ProjectionY( Form("histoMappingMidPtPhiInRMonteCarlo_%i",iR), histoMappingMidPtRPhiMonteCarlo->GetXaxis()->FindBin(arrayRBins[iR]), histoMappingMidPtRPhiMonteCarlo->GetXaxis()->FindBin(arrayRBins[iR+1]));
 		ConvGammaRebinWithBinCorrection(histoMappingMidPtPhiInRMonteCarlo[iR],rebinPhiPlots);
		GammaScalingHistogramm(histoMappingMidPtPhiInRMonteCarlo[iR],normFactorReconstMonteCarlo*1/600*1/(TMath::Abs(histoMappingMidPtRPhiMonteCarlo->GetXaxis()->GetBinUpEdge(histoMappingMidPtRPhiMonteCarlo->GetXaxis()->FindBin(arrayRBins[iR+1]))-histoMappingMidPtRPhiMonteCarlo->GetXaxis()->GetBinLowEdge(histoMappingMidPtRPhiMonteCarlo->GetXaxis()->FindBin(arrayRBins[iR])))));		
	}

	TH2F* 	histoMappingMidPtZPhiMonteCarlo = (TH2F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMappingMidPt_ZPhi");
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingMidPtPhiInZMonteCarlo[iZ]=		histoMappingMidPtZPhiMonteCarlo->ProjectionY( Form("histoMappingMidPtPhiInZMonteCarlo_%i",iZ), histoMappingMidPtZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ]), histoMappingMidPtZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]));
		ConvGammaRebinWithBinCorrection(histoMappingMidPtPhiInZMonteCarlo[iZ],rebinPhiPlots);
		GammaScalingHistogramm(histoMappingMidPtPhiInZMonteCarlo[iZ],normFactorReconstMonteCarlo*1/180*1/(TMath::Abs(histoMappingMidPtZPhiMonteCarlo->GetXaxis()->GetBinUpEdge(histoMappingMidPtZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoMappingMidPtZPhiMonteCarlo->GetXaxis()->GetBinLowEdge(histoMappingMidPtZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ])))));
	}

	TH2F* 	histoMappingMidPtZRMonteCarlo = (TH2F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMappingMidPt_ZR");
	TH1D * histoConversionPointVsRMidPtMonteCarlo=	(TH1D*)histoMappingMidPtZRMonteCarlo->ProjectionY("histoConversionPointVsRMidPtMonteCarlo");
	ConvGammaRebinWithBinCorrection(histoConversionPointVsRMidPtMonteCarlo,rebinIntPlots);
	GammaScalingHistogramm(histoConversionPointVsRMidPtMonteCarlo, normFactorReconstMonteCarlo);
	
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingMidPtRInZMonteCarlo[iZ]=		histoMappingMidPtZRMonteCarlo->ProjectionY( Form("histoMappingMidPtRInZMonteCarlo_%i",iZ), histoMappingMidPtZRMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ]), histoMappingMidPtZRMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]));
 		ConvGammaRebinWithBinCorrection(histoMappingMidPtRInZMonteCarlo[iZ],rebin);
		GammaScalingHistogramm(histoMappingMidPtRInZMonteCarlo[iZ],normFactorReconstMonteCarlo*1/(2*TMath::Pi())*1/(TMath::Abs(histoMappingMidPtZRMonteCarlo->GetXaxis()->GetBinUpEdge(histoMappingMidPtZRMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoMappingMidPtZRMonteCarlo->GetXaxis()->GetBinLowEdge(histoMappingMidPtZRMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ])))));
	}

	for(Int_t iR = 0; iR < nBinsR; iR++){
		histoMappingMidPtZInRMonteCarlo[iR]=		histoMappingMidPtZRMonteCarlo->ProjectionX( Form("histoMappingMidPtZInRMonteCarlo_%i",iR), histoMappingMidPtZRMonteCarlo->GetYaxis()->FindBin(arrayRBins[iR]), histoMappingMidPtZRMonteCarlo->GetYaxis()->FindBin(arrayRBins[iR+1]));
		ConvGammaRebinWithBinCorrection(histoMappingMidPtZInRMonteCarlo[iR],rebin);
		GammaScalingHistogramm(histoMappingMidPtZInRMonteCarlo[iR],normFactorReconstMonteCarlo*1/(2*TMath::Pi())*1/(TMath::Abs(histoMappingMidPtZRMonteCarlo->GetYaxis()->GetBinUpEdge(histoMappingMidPtZRMonteCarlo->GetYaxis()->FindBin(arrayRBins[iR+1]))-histoMappingMidPtZRMonteCarlo->GetYaxis()->GetBinLowEdge(histoMappingMidPtZRMonteCarlo->GetYaxis()->FindBin(arrayRBins[iR])))));
	}

	
	TH2F* 	histoTrueMappingRPhiMonteCarlo = (TH2F*)directoryGammaConvMonteCarlo->Get("ESD_TrueConversionMapping_RPhi");
	TH1D* 	histoTrueMappingPhiInRMonteCarlo[nBinsR];
	for(Int_t iR = 0; iR < nBinsR; iR++){
		histoTrueMappingPhiInRMonteCarlo[iR]=		histoTrueMappingRPhiMonteCarlo->ProjectionY( Form("histoTrueMappingPhiInRMonteCarlo_%i",iR), histoTrueMappingRPhiMonteCarlo->GetXaxis()->FindBin(arrayRBins[iR]), histoTrueMappingRPhiMonteCarlo->GetXaxis()->FindBin(arrayRBins[iR+1]));
 		ConvGammaRebinWithBinCorrection(histoTrueMappingPhiInRMonteCarlo[iR],rebinPhiPlots);
		GammaScalingHistogramm(histoTrueMappingPhiInRMonteCarlo[iR],normFactorReconstMonteCarlo*1/600*1/(TMath::Abs(histoTrueMappingRPhiMonteCarlo->GetXaxis()->GetBinUpEdge(histoTrueMappingRPhiMonteCarlo->GetXaxis()->FindBin(arrayRBins[iR+1]))-histoTrueMappingRPhiMonteCarlo->GetXaxis()->GetBinLowEdge(histoTrueMappingRPhiMonteCarlo->GetXaxis()->FindBin(arrayRBins[iR])))));
	}

	TH2F* 	histoTrueMappingZPhiMonteCarlo = (TH2F*)directoryGammaConvMonteCarlo->Get("ESD_TrueConversionMapping_ZPhi");
	TH1D* 	histoTrueMappingPhiInZMonteCarlo[nBinsZ];
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoTrueMappingPhiInZMonteCarlo[iZ]=		histoTrueMappingZPhiMonteCarlo->ProjectionY( Form("histoTrueMappingPhiInZMonteCarlo_%i",iZ), histoTrueMappingZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ]), histoTrueMappingZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]));
		ConvGammaRebinWithBinCorrection(histoTrueMappingPhiInZMonteCarlo[iZ],rebinPhiPlots);
		GammaScalingHistogramm(histoTrueMappingPhiInZMonteCarlo[iZ],normFactorReconstMonteCarlo*1/180*1/(TMath::Abs(histoTrueMappingZPhiMonteCarlo->GetXaxis()->GetBinUpEdge(histoTrueMappingZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoTrueMappingZPhiMonteCarlo->GetXaxis()->GetBinLowEdge(histoTrueMappingZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ])))));
	}

	TH2F* 	histoTrueMappingZRMonteCarlo = (TH2F*)directoryGammaConvMonteCarlo->Get("ESD_TrueConversionMapping_ZR");
	TH1D* 	histoTrueMappingRInZMonteCarlo[nBinsZ];		
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoTrueMappingRInZMonteCarlo[iZ]=		histoTrueMappingZRMonteCarlo->ProjectionY( Form("histoTrueMappingRInZMonteCarlo_%i",iZ), histoTrueMappingZRMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ]), histoTrueMappingZRMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]));
 		ConvGammaRebinWithBinCorrection(histoTrueMappingRInZMonteCarlo[iZ],rebin);
		GammaScalingHistogramm(histoTrueMappingRInZMonteCarlo[iZ],normFactorReconstMonteCarlo*1/(2*TMath::Pi())*1/(TMath::Abs(histoTrueMappingZRMonteCarlo->GetXaxis()->GetBinUpEdge(histoTrueMappingZRMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoTrueMappingZRMonteCarlo->GetXaxis()->GetBinLowEdge(histoTrueMappingZRMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ])))));
	}
	TH1D* 	histoTrueMappingZInRMonteCarlo[nBinsR];
	for(Int_t iR = 0; iR < nBinsR; iR++){
		histoTrueMappingZInRMonteCarlo[iR]=		histoTrueMappingZRMonteCarlo->ProjectionX( Form("histoTrueMappingZInRMonteCarlo_%i",iR), histoTrueMappingZRMonteCarlo->GetYaxis()->FindBin(arrayRBins[iR]), histoTrueMappingZRMonteCarlo->GetYaxis()->FindBin(arrayRBins[iR+1]));
		ConvGammaRebinWithBinCorrection(histoTrueMappingZInRMonteCarlo[iR],rebin);
		GammaScalingHistogramm(histoTrueMappingZInRMonteCarlo[iR],normFactorReconstMonteCarlo*1/(2*TMath::Pi())*1/(TMath::Abs(histoTrueMappingZRMonteCarlo->GetYaxis()->GetBinUpEdge(histoTrueMappingZRMonteCarlo->GetYaxis()->FindBin(arrayRBins[iR+1]))-histoTrueMappingZRMonteCarlo->GetYaxis()->GetBinLowEdge(histoTrueMappingZRMonteCarlo->GetYaxis()->FindBin(arrayRBins[iR])))));
	}

	TH2F* 	histoMappingSPDZPhiMonteCarlo = (TH2F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMapping_SPD_ZPhi");
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingSPDPhiInZMonteCarlo[iZ]=		histoMappingSPDZPhiMonteCarlo->ProjectionY( Form("histoMappingPhiInZMonteCarlo_SPD_%i",iZ), histoMappingSPDZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ]), histoMappingSPDZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]));
 		ConvGammaRebinWithBinCorrection(histoMappingSPDPhiInZMonteCarlo[iZ],rebinPhiPlots);
		GammaScalingHistogramm(histoMappingSPDPhiInZMonteCarlo[iZ],normFactorReconstMonteCarlo*1/9.5*1/(TMath::Abs(histoMappingSPDZPhiMonteCarlo->GetXaxis()->GetBinUpEdge(histoMappingSPDZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoMappingSPDZPhiMonteCarlo->GetXaxis()->GetBinLowEdge(histoMappingSPDZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ])))));
		nameHistoRatioSPDPhiInZMonteCarlo=Form("histoMappingRatioSPDPhiInZ_%02d",iZ);
		histoMappingRatioSPDPhiInZ[iZ]= 		(TH1F*)histoMappingSPDPhiInZData[iZ]->Clone();
		histoMappingRatioSPDPhiInZ[iZ]->SetName(nameHistoRatioSPDPhiInZMonteCarlo); 
		histoMappingRatioSPDPhiInZ[iZ]->Divide(histoMappingRatioSPDPhiInZ[iZ],histoMappingSPDPhiInZMonteCarlo[iZ]);
		
	}
	
	TH2F* 	histoMappingSPDThZPhiMonteCarlo = (TH2F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMapping_SPDTh_ZPhi");
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingSPDThPhiInZMonteCarlo[iZ]=		histoMappingSPDThZPhiMonteCarlo->ProjectionY( Form("histoMappingPhiInZMonteCarlo_SPDTh_%i",iZ), histoMappingSPDThZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ]), histoMappingSPDThZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]));
 		ConvGammaRebinWithBinCorrection(histoMappingSPDThPhiInZMonteCarlo[iZ],rebinPhiPlots);
		GammaScalingHistogramm(histoMappingSPDThPhiInZMonteCarlo[iZ],normFactorReconstMonteCarlo*1/3.5*1/(TMath::Abs(histoMappingSPDThZPhiMonteCarlo->GetXaxis()->GetBinUpEdge(histoMappingSPDThZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoMappingSPDThZPhiMonteCarlo->GetXaxis()->GetBinLowEdge(histoMappingSPDThZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ])))));
		nameHistoRatioSPDThPhiInZMonteCarlo=Form("histoMappingRatioSPDThPhiInZ_%02d",iZ);
		histoMappingRatioSPDThPhiInZ[iZ]= 		(TH1F*)histoMappingSPDThPhiInZData[iZ]->Clone();
		histoMappingRatioSPDThPhiInZ[iZ]->SetName(nameHistoRatioSPDThPhiInZMonteCarlo); 
		histoMappingRatioSPDThPhiInZ[iZ]->Divide(histoMappingRatioSPDThPhiInZ[iZ],histoMappingSPDThPhiInZMonteCarlo[iZ]);
		
	}
	
	TH1D*		histoMappingSDDPhiInZMonteCarlo[nBinsZ];
	TH2F* 	histoMappingSDDZPhiMonteCarlo = (TH2F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMapping_SDD_ZPhi");
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingSDDPhiInZMonteCarlo[iZ]=		histoMappingSDDZPhiMonteCarlo->ProjectionY( Form("histoMappingPhiInZ_SDD_MonteCarlo_%i",iZ), histoMappingSDDZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ]), histoMappingSDDZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]));
		ConvGammaRebinWithBinCorrection(histoMappingSDDPhiInZMonteCarlo[iZ],rebinPhiPlots);
		GammaScalingHistogramm(histoMappingSDDPhiInZMonteCarlo[iZ],normFactorReconstMonteCarlo*1/22.*1/(TMath::Abs(histoMappingSDDZPhiMonteCarlo->GetXaxis()->GetBinUpEdge(histoMappingSDDZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoMappingSDDZPhiMonteCarlo->GetXaxis()->GetBinLowEdge(histoMappingSDDZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ])))));
		nameHistoRatioSDDPhiInZMonteCarlo=Form("histoMappingRatioSDDPhiInZ_%02d",iZ);
		histoMappingRatioSDDPhiInZ[iZ]= 		(TH1F*)histoMappingSDDPhiInZData[iZ]->Clone();
		histoMappingRatioSDDPhiInZ[iZ]->SetName(nameHistoRatioSDDPhiInZMonteCarlo); 
		histoMappingRatioSDDPhiInZ[iZ]->Divide(histoMappingRatioSDDPhiInZ[iZ],histoMappingSDDPhiInZMonteCarlo[iZ]);
		
	}

	TH1D*		histoMappingSSDPhiInZMonteCarlo[nBinsZ];
	TH2F* 	histoMappingSSDZPhiMonteCarlo = (TH2F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMapping_SSD_ZPhi");
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingSSDPhiInZMonteCarlo[iZ]=		histoMappingSSDZPhiMonteCarlo->ProjectionY( Form("histoMappingPhiInZ_SSD_MonteCarlo_%i",iZ), histoMappingSSDZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ]), histoMappingSSDZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]));
 		ConvGammaRebinWithBinCorrection(histoMappingSSDPhiInZMonteCarlo[iZ],rebinPhiPlots);
		GammaScalingHistogramm(histoMappingSSDPhiInZMonteCarlo[iZ],normFactorReconstMonteCarlo*1/20.*1/(TMath::Abs(histoMappingSSDZPhiMonteCarlo->GetXaxis()->GetBinUpEdge(histoMappingSSDZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoMappingSSDZPhiMonteCarlo->GetXaxis()->GetBinLowEdge(histoMappingSSDZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ])))));
		nameHistoRatioSSDPhiInZMonteCarlo=Form("histoMappingRatioSSDPhiInZ_%02d",iZ);
		histoMappingRatioSSDPhiInZ[iZ]= 		(TH1F*)histoMappingSSDPhiInZData[iZ]->Clone();
		histoMappingRatioSSDPhiInZ[iZ]->SetName(nameHistoRatioSSDPhiInZMonteCarlo); 
		histoMappingRatioSSDPhiInZ[iZ]->Divide(histoMappingRatioSSDPhiInZ[iZ],histoMappingSSDPhiInZMonteCarlo[iZ]);
		
	}

	
	TH2F* 	histoMappingITSTPCZPhiMonteCarlo = (TH2F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMapping_ITSTPC_ZPhi");
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingITSTPCPhiInZMonteCarlo[iZ]=		histoMappingITSTPCZPhiMonteCarlo->ProjectionY( Form("histoMappingPhiInZMonteCarlo_ITSTPC_%i",iZ), histoMappingITSTPCZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ]), histoMappingITSTPCZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1])); 		ConvGammaRebinWithBinCorrection(histoMappingITSTPCPhiInZMonteCarlo[iZ],rebinPhiPlots);
		GammaScalingHistogramm(histoMappingITSTPCPhiInZMonteCarlo[iZ],normFactorReconstMonteCarlo*1/17.5*1/(TMath::Abs(histoMappingITSTPCZPhiMonteCarlo->GetXaxis()->GetBinUpEdge(histoMappingITSTPCZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ+1]))-histoMappingITSTPCZPhiMonteCarlo->GetXaxis()->GetBinLowEdge(histoMappingITSTPCZPhiMonteCarlo->GetXaxis()->FindBin(arrayZBins[iZ])))));
		nameHistoRatioITSTPCPhiInZMonteCarlo=Form("histoMappingRatioITSTPCPhiInZ_%02d",iZ);
		histoMappingRatioITSTPCPhiInZ[iZ]= 		(TH1F*)histoMappingITSTPCPhiInZData[iZ]->Clone();
		histoMappingRatioITSTPCPhiInZ[iZ]->SetName(nameHistoRatioITSTPCPhiInZMonteCarlo); 
		histoMappingRatioITSTPCPhiInZ[iZ]->Divide(histoMappingRatioITSTPCPhiInZ[iZ],histoMappingITSTPCPhiInZMonteCarlo[iZ]);
	}			       
	
	
	for (Int_t iX = 1 ; iX < histoConversionPointVsXYBPData->GetNbinsX(); iX ++){
		for (Int_t iY = 1; iY < histoConversionPointVsXYBPData->GetNbinsY(); iY ++) {
			Double_t rcalc = TMath::Sqrt(histoConversionPointVsXYBPData->GetXaxis()->GetBinCenter(iX)*histoConversionPointVsXYBPData->GetXaxis()->GetBinCenter(iX) + histoConversionPointVsXYBPData->GetYaxis()->GetBinCenter(iY)*histoConversionPointVsXYBPData->GetYaxis()->GetBinCenter(iY));
// 			if (rcalc < rmin){
// 				histoConversionPointVsXYBPData->SetBinContent(iX,iY,NULL);
// 			}
		}
	}

	for (Int_t iX = 1 ; iX < histoConversionPointVsXYBPMonteCarlo->GetNbinsX(); iX ++){
		for (Int_t iY = 1; iY < histoConversionPointVsXYBPMonteCarlo->GetNbinsY(); iY ++) {
			Double_t rcalc = TMath::Sqrt(histoConversionPointVsXYBPMonteCarlo->GetXaxis()->GetBinCenter(iX)*histoConversionPointVsXYBPMonteCarlo->GetXaxis()->GetBinCenter(iX) + histoConversionPointVsXYBPMonteCarlo->GetYaxis()->GetBinCenter(iY)*histoConversionPointVsXYBPMonteCarlo->GetYaxis()->GetBinCenter(iY));
// 			if (rcalc < rmin){
// 				histoConversionPointVsXYBPMonteCarlo->SetBinContent(iX,iY,NULL);
// 			}
		}
	}

	
	PlotStandard2D( histoConversionPointVsXYBPData , Form("%s/XY_distribution_BP_Data.%s",outputDirectory.Data(),suffix.Data()), "", "X (cm)", "Y (cm)", kTRUE, -15., 15., kTRUE, -15., 15.,kFALSE,  4.e-6, 1.e-3, 0, 1, floatLocationRightUp2DXY,1000,880,rightMargin, leftMargin,0.08, 0.02,"Perf");
	delete histoConversionPointVsXYBPData;
	PlotStandard2D( histoConversionPointVsXYBPMonteCarlo , Form("%s/XY_distribution_BP_MC.%s",outputDirectory.Data(),suffix.Data()), "", "X (cm)", "Y (cm)", kTRUE, -15., 15., kTRUE, -15., 15.,kFALSE,  4.e-6, 1.e-3, 0, 1, floatLocationRightUp2DXY,1000,880,rightMargin, leftMargin,0.08, 0.02,"Perf");
	delete histoConversionPointVsXYBPMonteCarlo;
	
	
	TLine * 	linePhi =	 		new TLine (-3.2,1,3.2,1);
	TLine * 	lineZ = 			new TLine (-300,1,300,1);
	TLine * 	lineR = 			new TLine (0,1,200,1);
	linePhi->SetLineColor(2);
	lineZ->SetLineColor(2);
	lineR->SetLineColor(2);
	
	Double_t integratedRData[nBinsR];       
	Double_t integratedRMonteCarlo[nBinsR];                 
	Double_t integratedMidPtRData[nBinsR];       
	Double_t integratedMidPtRMonteCarlo[nBinsR];                 
	
	
	for(Int_t iR = iRStart; iR < nBinsR; iR++){
		integratedRData[iR]=histoMappingPhiInRData[iR]->Integral() ;                 
		integratedRMonteCarlo[iR]=histoMappingPhiInRMonteCarlo[iR]->Integral() ;     
		integratedMidPtRData[iR]=histoMappingMidPtPhiInRData[iR]->Integral() ;                 
		integratedMidPtRMonteCarlo[iR]=histoMappingMidPtPhiInRMonteCarlo[iR]->Integral() ;     
	}
	
	Double_t integratedZInRData[nBinsR];       
	Double_t integratedZInRMonteCarlo[nBinsR];                 
	Double_t integratedMidPtZInRData[nBinsR];       
	Double_t integratedMidPtZInRMonteCarlo[nBinsR];                 
	
	for(Int_t iR = iRStart; iR < nBinsR; iR++){
		integratedZInRData[iR]=histoMappingZInRData[iR]->Integral() ;                        
		integratedZInRMonteCarlo[iR]=histoMappingZInRMonteCarlo[iR]->Integral() ;    
		integratedMidPtZInRData[iR]=histoMappingMidPtZInRData[iR]->Integral() ;                        
		integratedMidPtZInRMonteCarlo[iR]=histoMappingMidPtZInRMonteCarlo[iR]->Integral() ;    
	}
	
	Double_t integratedZData[nBinsZ];  
	Double_t integratedZMonteCarlo[nBinsZ];                    
	Double_t integratedMidPtZData[nBinsZ];  
	Double_t integratedMidPtZMonteCarlo[nBinsZ];                    
	
	for(Int_t iZ = minimumBinZNormal ; iZ < (nBinsZ-minimumBinZNormal); iZ++){
		integratedZData[iZ]=histoMappingPhiInZData[iZ]->Integral() ;                 
		integratedZMonteCarlo[iZ]=histoMappingPhiInZMonteCarlo[iZ]->Integral() ;     
		integratedMidPtZData[iZ]=histoMappingMidPtPhiInZData[iZ]->Integral() ;                 
		integratedMidPtZMonteCarlo[iZ]=histoMappingMidPtPhiInZMonteCarlo[iZ]->Integral() ;     
	}       
	
	Double_t integratedRInZData[nBinsZ];       
	Double_t integratedRInZMonteCarlo[nBinsZ];                 
	Double_t integratedMidPtRInZData[nBinsZ];       
	Double_t integratedMidPtRInZMonteCarlo[nBinsZ];                 
	
	for(Int_t iZ = minimumBinZNormal ; iZ < (nBinsZ-minimumBinZNormal); iZ++){
		integratedRInZData[iZ]=histoMappingRInZData[iZ]->Integral() ;                        
		integratedRInZMonteCarlo[iZ]=histoMappingRInZMonteCarlo[iZ]->Integral() ;            
		integratedMidPtRInZData[iZ]=histoMappingMidPtRInZData[iZ]->Integral() ;                        
		integratedMidPtRInZMonteCarlo[iZ]=histoMappingMidPtRInZMonteCarlo[iZ]->Integral() ;            

	}
				
	
	Double_t integratedBeforeFirstLayerData = histoConversionPointVsRData->Integral(histoConversionPointVsRData->GetXaxis()->FindBin(0.), histoConversionPointVsRData->GetXaxis()->FindBin(3.9));
	Double_t integratedAfterFirstLayerData = histoConversionPointVsRData->Integral(histoConversionPointVsRData->GetXaxis()->FindBin(0.), histoConversionPointVsRData->GetXaxis()->FindBin(4.5));
	Double_t integratedBeforeFirstLayerMonteCarlo = histoConversionPointVsRMonteCarlo->Integral(histoConversionPointVsRData->GetXaxis()->FindBin(0.), histoConversionPointVsRData->GetXaxis()->FindBin(3.9));
	Double_t integratedAfterFirstLayerMonteCarlo = histoConversionPointVsRMonteCarlo->Integral(histoConversionPointVsRData->GetXaxis()->FindBin(0.), histoConversionPointVsRData->GetXaxis()->FindBin(4.5));
	Double_t integratedBeforeFirstLayerMidPtData = histoConversionPointVsRMidPtData->Integral(histoConversionPointVsRMidPtData->GetXaxis()->FindBin(0.), histoConversionPointVsRMidPtData->GetXaxis()->FindBin(3.9));
	Double_t integratedAfterFirstLayerMidPtData = histoConversionPointVsRMidPtData->Integral(histoConversionPointVsRMidPtData->GetXaxis()->FindBin(0.), histoConversionPointVsRMidPtData->GetXaxis()->FindBin(4.5));
	Double_t integratedBeforeFirstLayerMidPtMonteCarlo = histoConversionPointVsRMidPtMonteCarlo->Integral(histoConversionPointVsRMidPtData->GetXaxis()->FindBin(0.), histoConversionPointVsRMidPtData->GetXaxis()->FindBin(3.9));
	Double_t integratedAfterFirstLayerMidPtMonteCarlo = histoConversionPointVsRMidPtMonteCarlo->Integral(histoConversionPointVsRMidPtData->GetXaxis()->FindBin(0.), histoConversionPointVsRMidPtData->GetXaxis()->FindBin(4.5));
	
	
	Double_t sumPhiInRData = 0., 		sumPhiInZData =0.;
	Double_t sumPhiInRMonteCarlo = 0., sumPhiInZMonteCarlo =0.;
	Double_t sumZinRData = 0., 		sumRinZData =0.;
	Double_t sumZinRMonteCarlo = 0., 	sumRinZMonteCarlo =0.;
	
	for(Int_t iR = iRStart; iR < nBinsR; iR++){
			sumPhiInRData = 			sumPhiInRData + integratedRData[iR];
			sumPhiInRMonteCarlo = 		sumPhiInRMonteCarlo + integratedRMonteCarlo[iR];
			sumZinRData = 				sumZinRData + integratedZInRData[iR];
			sumZinRMonteCarlo = 		sumZinRMonteCarlo + integratedZInRMonteCarlo[iR];
	}
	for(Int_t iZ = 1; iZ < nBinsZ-1; iZ++){
			sumPhiInZData = 			sumPhiInZData + integratedZData[iZ];
			sumPhiInZMonteCarlo = 		sumPhiInZMonteCarlo + integratedZMonteCarlo[iZ];
			sumRinZData = 				sumRinZData + integratedRInZData[iZ];
			sumRinZMonteCarlo = 		sumRinZMonteCarlo + integratedRInZMonteCarlo[iZ];
	}
	
	Double_t sumMidPtPhiInRData = 0., 		sumMidPtPhiInZData =0.;
	Double_t sumMidPtPhiInRMonteCarlo = 0., sumMidPtPhiInZMonteCarlo =0.;
	Double_t sumMidPtZinRData = 0., 		sumMidPtRinZData =0.;
	Double_t sumMidPtZinRMonteCarlo = 0., 	sumMidPtRinZMonteCarlo =0.;
	
	for(Int_t iR = iRStart; iR < nBinsR; iR++){
			sumMidPtPhiInRData = 			sumMidPtPhiInRData + integratedMidPtRData[iR];
			sumMidPtPhiInRMonteCarlo = 		sumMidPtPhiInRMonteCarlo + integratedMidPtRMonteCarlo[iR];
			sumMidPtZinRData = 				sumMidPtZinRData + integratedMidPtZInRData[iR];
			sumMidPtZinRMonteCarlo = 		sumMidPtZinRMonteCarlo + integratedMidPtZInRMonteCarlo[iR];
	}
	for(Int_t iZ = 1; iZ < nBinsZ-1; iZ++){
			sumMidPtPhiInZData = 			sumMidPtPhiInZData + integratedMidPtZData[iZ];
			sumMidPtPhiInZMonteCarlo = 		sumMidPtPhiInZMonteCarlo + integratedMidPtZMonteCarlo[iZ];
			sumMidPtRinZData = 				sumMidPtRinZData + integratedMidPtRInZData[iZ];
			sumMidPtRinZMonteCarlo = 		sumMidPtRinZMonteCarlo + integratedMidPtRInZMonteCarlo[iZ];
	}
	
	
	const char *nameDataOutputFile = Form("%s/SystematicErrorMaterialBudget%s_%s_%s.dat",outputDirectory1.Data(),optEnergy.Data(),textGenerator.Data(),textPeriod.Data());
	fstream fileDatOutput;
	fileDatOutput.open(nameDataOutputFile, ios::out);
	fileDatOutput << "#Calculating Integrals" << endl;       
	fileDatOutput << "------------------------------------------------------------------------------------------" << endl;
	fileDatOutput << "# This file is created to display the Integrals of the different bins of the different diagrams. The forth column displays the Integral over the data in that bin divided by the sum over all bins of data, the same does the sixth column for montecarlo. The eigth column displays the difference of data - montecarlo divided by montecarlo of that bin. Therefore you should put the data in the first input and the Montecarlo in the second." << endl;
	fileDatOutput << "------------------------------------------------------------------------------------------" << endl;
	
	fileDatOutput << "data :\t" 			<< data 	 	<< endl;
	fileDatOutput << "montecarlo :\t" 		<< MCfile 	<< endl;
	fileDatOutput << "------------------------------------------------------------------------------------------" << endl;
	fileDatOutput << endl;
	fileDatOutput << "\t data \t montecarlo " 		<< endl;
	fileDatOutput << "Number of events" 			<< "\t" 	<< numberGoodEventsData 		<< "\t" 	<< numberGoodEventsMonteCarlo 	<< endl;
	fileDatOutput << "Number of triggers" 			<< "\t" 	<< numberGoodTriggerData 	<< "\t" 	<< numberGoodTriggerMonteCarlo 	<< endl;
	fileDatOutput << "Number reconstructed gammas"	<< "\t" 	<< numberReconstGammaData 	<< "\t" 	<< numberReconstGammaMonteCarlo 	<<endl;
	fileDatOutput << "Mean Multiplicity" 			<< "\t" 	<< meanMultiplitcityData 	<< "\t" 	<< meanMultiplitcityMonteCarlo 	<< endl;
	fileDatOutput << endl;
	fileDatOutput << endl;
	
	fileDatOutput << "------------------------------------------------------------------------------------------" << endl;           
	fileDatOutput << "graph \t bin \t data \t stat err data \t % \t montecarlo \t stat err MC\t % \t data-montecarlo\t stat err D-MC \t % \t stat err %" <<endl;
	
	Double_t 	positiveErrorAllpt = 0;
	Double_t 	negativeErrorAllpt = 0;
	Double_t 	positiveErrorMidpt = 0;
	Double_t 	negativeErrorMidpt = 0;
	
	TGraph* graphIntegratedRData =new TGraph(nBinsR-iRStart);
	TGraph* graphIntegratedRMonteCarlo =new TGraph(nBinsR-iRStart);
	TGraph* graphIntegratedMidPtRData =new TGraph(nBinsR-iRStart);
	TGraph* graphIntegratedMidPtRMonteCarlo =new TGraph(nBinsR-iRStart);
// 	TGraph* graphIntegratedZData =0x00;
// 	TGraph* graphIntegratedZMonteCarlo =0x00;
// 	TGraph* graphIntegratedMidPtZData =0x00;
// 	TGraph* graphIntegratedMidPtZMonteCarlo =0x00;
	TGraph* graphRelativeDevR = new TGraph(nBinsR-iRStart);
// 	TGraph* graphRelativeDevZ = 0x00;
	TGraph* graphRelativeDevMidPtR = new TGraph(nBinsR-iRStart);
// 	TGraph* graphRelativeDevMitPtZ = 0x00;
	TGraph* graphTotalErrors = new TGraph(10);
	
	
	fileDatOutput << "Phi in R" <<endl;
//                      fileDatOutput <<"\t0\t" << integratedRData[0] << "\t" << integratedRData[0]/sumPhiInRData * 100 <<"\t" << integratedRMonteCarlo[0] << "\t" << integratedRMonteCarlo[0]/sumPhiInRMonteCarlo * 100 << "\t" << integratedRData[0]- integratedRMonteCarlo[0] << "\t" << (integratedRData[0]- integratedRMonteCarlo[0])/integratedRMonteCarlo[0] *100 << endl;   

	Double_t deviation= 0;
	for(Int_t iR = iRStart; iR < nBinsR; iR++){
			deviation = CalculateErrors(fileDatOutput, iR, integratedRData[iR], sumPhiInRData, numberGoodEventsData, integratedRMonteCarlo[iR], sumPhiInRMonteCarlo, numberGoodEventsMonteCarlo);
			if(deviation < 0){ 
				negativeErrorAllpt = negativeErrorAllpt + deviation;
			} else {
				positiveErrorAllpt = positiveErrorAllpt + deviation;
			}  
	}
	fileDatOutput<< endl;
	fileDatOutput << "Relative Errors" << endl;
	for(Int_t iR = iRStart; iR < nBinsR; iR++){
		fileDatOutput 	<< "\t" 	<< iR 	
					<< "\t" 	<< integratedRData[iR] 						<< "\t" << integratedRData[iR]/sumPhiInRData * 100 
					<< "\t" 	<< integratedRMonteCarlo[iR] 					<< "\t" << integratedRMonteCarlo[iR]/sumPhiInRMonteCarlo * 100 
					<< "\t" 	<< integratedRData[iR]- integratedRMonteCarlo[iR] << "\t" << (integratedRData[iR]- integratedRMonteCarlo[iR])/integratedRMonteCarlo[iR]* 100 << endl;
		graphIntegratedRData->SetPoint(iR-iRStart, (arrayRBins[iR+1] + arrayRBins[iR])/2., integratedRData[iR]);
		graphIntegratedRMonteCarlo->SetPoint(iR-iRStart, (arrayRBins[iR+1] + arrayRBins[iR])/2., integratedRMonteCarlo[iR]);			
		graphRelativeDevR->SetPoint(iR-iRStart, (arrayRBins[iR+1] + arrayRBins[iR])/2., (integratedRData[iR]- integratedRMonteCarlo[iR])/integratedRMonteCarlo[iR]* 100);
		cout << arrayRBins[iR] << "\t" << arrayRBins[iR+1] << "\t" <<(arrayRBins[iR+1] + arrayRBins[iR])/2. << endl;
	}
	fileDatOutput 	<< "\t" << "sum" 
				<< "\t" << sumPhiInRData	 					<< "\t" << "100" 
				<< "\t" << sumPhiInRMonteCarlo 				<< "\t" << "100" 
				<< "\t" << sumPhiInRData - sumPhiInRMonteCarlo 	<< "\t" << (sumPhiInRData - sumPhiInRMonteCarlo)/sumPhiInRMonteCarlo *100 	<< endl;
	fileDatOutput 	<< "\t" << "Error neg \t" 	<<  negativeErrorAllpt/sumPhiInRMonteCarlo *100 	
				<< "\t" << "Error pos \t" 	<<  positiveErrorAllpt/sumPhiInRMonteCarlo *100 << endl;
				
	graphTotalErrors->SetPoint(0,0.,(sumPhiInRData - sumPhiInRMonteCarlo)/sumPhiInRMonteCarlo *100);
	graphTotalErrors->SetPoint(1,1.,positiveErrorAllpt/sumPhiInRMonteCarlo *100);
	graphTotalErrors->SetPoint(2,-1.,negativeErrorAllpt/sumPhiInRMonteCarlo *100);
	fileDatOutput << endl;
	fileDatOutput << "------------------------------------------------------------------------------------------" << endl;           
	
	fileDatOutput 	<< "\t" << "before 1st layer" 
				<< "\t" << integratedBeforeFirstLayerData	 					<< "\t" << integratedBeforeFirstLayerData/histoConversionPointVsRData->Integral() * 100
				<< "\t" << integratedBeforeFirstLayerMonteCarlo 				<< "\t" <<  integratedBeforeFirstLayerMonteCarlo/histoConversionPointVsRMonteCarlo->Integral() * 100 
				<< "\t" << integratedBeforeFirstLayerData - integratedBeforeFirstLayerMonteCarlo 	<< "\t" << (integratedBeforeFirstLayerData - integratedBeforeFirstLayerMonteCarlo)/integratedBeforeFirstLayerMonteCarlo *100 	<< endl;
	fileDatOutput 	<< "\t" << "after 1st layer" 
				<< "\t" << integratedAfterFirstLayerData	 					<< "\t" << integratedAfterFirstLayerData/ histoConversionPointVsRData->Integral() * 100
				<< "\t" << integratedAfterFirstLayerMonteCarlo 				<< "\t" <<  integratedAfterFirstLayerMonteCarlo/histoConversionPointVsRMonteCarlo->Integral() * 100 
				<< "\t" << integratedAfterFirstLayerData - integratedAfterFirstLayerMonteCarlo 	<< "\t" << (integratedAfterFirstLayerData - integratedAfterFirstLayerMonteCarlo)/integratedAfterFirstLayerMonteCarlo *100 	<< endl;
	fileDatOutput << "------------------------------------------------------------------------------------------" << endl;           
	graphTotalErrors->SetPoint(6,2.,(integratedBeforeFirstLayerData - integratedBeforeFirstLayerMonteCarlo)/integratedBeforeFirstLayerMonteCarlo *100);
	graphTotalErrors->SetPoint(7,3.,(integratedAfterFirstLayerData - integratedAfterFirstLayerMonteCarlo)/integratedAfterFirstLayerMonteCarlo *100);
	
	
	
	fileDatOutput << "Phi in R - Mid Pt" <<endl;
//                      fileDatOutput <<"\t0\t" << integratedRData[0] << "\t" << integratedRData[0]/sumPhiInRData * 100 <<"\t" << integratedRMonteCarlo[0] << "\t" << integratedRMonteCarlo[0]/sumPhiInRMonteCarlo * 100 << "\t" << integratedRData[0]- integratedRMonteCarlo[0] << "\t" << (integratedRData[0]- integratedRMonteCarlo[0])/integratedRMonteCarlo[0] *100 << endl;   

	Double_t deviationMidPt= 0;
	for(Int_t iR = iRStart; iR < nBinsR; iR++){
			deviationMidPt = CalculateErrors(fileDatOutput, iR, integratedMidPtRData[iR], sumMidPtPhiInRData, numberGoodEventsData, integratedMidPtRMonteCarlo[iR], sumMidPtPhiInRMonteCarlo, numberGoodEventsMonteCarlo);
			if(deviationMidPt < 0){ 
				negativeErrorMidpt = negativeErrorMidpt + deviationMidPt;
			} else {
				positiveErrorMidpt = positiveErrorMidpt + deviationMidPt;
			}  
	}
	fileDatOutput<< endl;
	fileDatOutput << "Relative Errors" << endl;
	for(Int_t iR = iRStart; iR < nBinsR; iR++){
		fileDatOutput 	<< "\t" 	<< iR 	
					<< "\t" 	<< integratedMidPtRData[iR] 						<< "\t" << integratedMidPtRData[iR]/sumMidPtPhiInRData * 100 
					<< "\t" 	<< integratedMidPtRMonteCarlo[iR] 					<< "\t" << integratedMidPtRMonteCarlo[iR]/sumMidPtPhiInRMonteCarlo * 100 
					<< "\t" 	<< integratedMidPtRData[iR]- integratedMidPtRMonteCarlo[iR] << "\t" << (integratedMidPtRData[iR]- integratedMidPtRMonteCarlo[iR])/integratedMidPtRMonteCarlo[iR]* 100 << endl;
		graphIntegratedMidPtRData->SetPoint(iR-iRStart, (arrayRBins[iR+1] + arrayRBins[iR])/2., integratedMidPtRData[iR]);
		graphIntegratedMidPtRMonteCarlo->SetPoint(iR-iRStart, (arrayRBins[iR+1] + arrayRBins[iR])/2., integratedMidPtRMonteCarlo[iR]);			
		graphRelativeDevMidPtR->SetPoint(iR-iRStart, (arrayRBins[iR+1] + arrayRBins[iR])/2., (integratedMidPtRData[iR]- integratedMidPtRMonteCarlo[iR])/integratedMidPtRMonteCarlo[iR]* 100);
	}
	fileDatOutput 	<< "\t" << "sum" 
				<< "\t" << sumMidPtPhiInRData	 					<< "\t" << "100" 
				<< "\t" << sumMidPtPhiInRMonteCarlo 				<< "\t" << "100" 
				<< "\t" << sumMidPtPhiInRData - sumMidPtPhiInRMonteCarlo 	<< "\t" << (sumMidPtPhiInRData - sumMidPtPhiInRMonteCarlo)/sumMidPtPhiInRMonteCarlo *100 	<< endl;
	fileDatOutput 	<< "\t" << "Error neg \t" 	<<  negativeErrorMidpt/sumMidPtPhiInRMonteCarlo *100 	
				<< "\t" << "Error pos \t" 	<<  positiveErrorMidpt/sumMidPtPhiInRMonteCarlo *100 << endl;
	
	graphTotalErrors->SetPoint(3,0.5,(sumMidPtPhiInRData - sumMidPtPhiInRMonteCarlo)/sumMidPtPhiInRMonteCarlo *100);
	graphTotalErrors->SetPoint(4,1.5,positiveErrorMidpt/sumMidPtPhiInRMonteCarlo *100);
	graphTotalErrors->SetPoint(5,-1.5,negativeErrorMidpt/sumMidPtPhiInRMonteCarlo *100);

	fileDatOutput << "------------------------------------------------------------------------------------------" << endl;           
	
	fileDatOutput 	<< "\t" << "before 1st layer" 
				<< "\t" << integratedBeforeFirstLayerMidPtData	 					<< "\t" << integratedBeforeFirstLayerMidPtData/histoConversionPointVsRMidPtData->Integral() * 100
				<< "\t" << integratedBeforeFirstLayerMidPtMonteCarlo 				<< "\t" <<  integratedBeforeFirstLayerMonteCarlo/histoConversionPointVsRMidPtMonteCarlo->Integral() * 100 
				<< "\t" << integratedBeforeFirstLayerMidPtData - integratedBeforeFirstLayerMidPtMonteCarlo 	<< "\t" << (integratedBeforeFirstLayerMidPtData - integratedBeforeFirstLayerMidPtMonteCarlo)/integratedBeforeFirstLayerMidPtMonteCarlo *100 	<< endl;
	fileDatOutput 	<< "\t" << "after 1st layer" 
				<< "\t" << integratedAfterFirstLayerMidPtData	 					<< "\t" << integratedAfterFirstLayerMidPtData/ histoConversionPointVsRMidPtData->Integral() * 100
				<< "\t" << integratedAfterFirstLayerMidPtMonteCarlo 				<< "\t" <<  integratedAfterFirstLayerMidPtMonteCarlo/histoConversionPointVsRMidPtMonteCarlo->Integral() * 100 
				<< "\t" << integratedAfterFirstLayerMidPtData - integratedAfterFirstLayerMidPtMonteCarlo 	<< "\t" << (integratedAfterFirstLayerMidPtData - integratedAfterFirstLayerMidPtMonteCarlo)/integratedAfterFirstLayerMidPtMonteCarlo *100 	<< endl;
	
	graphTotalErrors->SetPoint(8,2.,(integratedBeforeFirstLayerMidPtData - integratedBeforeFirstLayerMidPtMonteCarlo)/integratedBeforeFirstLayerMidPtMonteCarlo *100);
	graphTotalErrors->SetPoint(9,3.,(integratedAfterFirstLayerMidPtData - integratedAfterFirstLayerMidPtMonteCarlo)/integratedAfterFirstLayerMidPtMonteCarlo *100);
	
	fileDatOutput << endl;
	fileDatOutput << "------------------------------------------------------------------------------------------" << endl;           
	
	
	fileDatOutput << "Phi in Z" <<endl;
	for(Int_t iZ = 1; iZ < nBinsZ-1; iZ++){
		fileDatOutput 	<< "\t" 	<< iZ
					<< "\t" 	<< integratedZData[iZ] 						<< "\t" << integratedZData[iZ]/sumPhiInZData * 100 
					<< "\t" 	<< integratedZMonteCarlo[iZ] 					<< "\t" << integratedZMonteCarlo[iZ]/sumPhiInZMonteCarlo * 100 
					<< "\t" 	<< integratedZData[iZ]- integratedZMonteCarlo[iZ] << "\t" << (integratedZData[iZ]- integratedZMonteCarlo[iZ])/integratedZMonteCarlo[iZ]*100 << endl;
	}
	fileDatOutput 	<< "\t" << "sum" 
				<< "\t" << sumPhiInZData 					<< "\t" << "100" 
				<< "\t" << sumPhiInZMonteCarlo 				<< "\t" << "100" 
				<< "\t" << sumPhiInZData - sumPhiInZMonteCarlo 	<< "\t" << (sumPhiInZData-sumPhiInZMonteCarlo)/sumPhiInZMonteCarlo *100 << endl;
	fileDatOutput << endl;
			
	fileDatOutput << "Z in R" <<endl;
	for(Int_t iR = iRStart; iR < nBinsR; iR++){
		fileDatOutput 	<< "\t" 	<< iR 
					<< "\t" 	<< integratedZInRData[iR] 							<< "\t" << integratedZInRData[iR]/sumZinRData * 100 
					<< "\t" 	<< integratedZInRMonteCarlo[iR] 						<< "\t" << integratedZInRMonteCarlo[iR]/sumZinRMonteCarlo * 100 
					<< "\t" 	<< integratedZInRData[iR]- integratedZInRMonteCarlo[iR] 	<< "\t" << (integratedZInRData[iR]- integratedZInRMonteCarlo[iR])/ integratedZInRMonteCarlo[iR]*100 << endl;
	}
	fileDatOutput 	<< "\t" << "sum" 
				<< "\t" << sumZinRData 					<< "\t" << "100" 
				<< "\t" << sumZinRMonteCarlo 				<< "\t" << "100" 
				<< "\t" << sumZinRData - sumZinRMonteCarlo 	<< "\t" << (sumZinRData-sumZinRMonteCarlo)/sumZinRMonteCarlo *100 << endl;
	fileDatOutput << endl;

	
	fileDatOutput << "R in Z" <<endl;
	for(Int_t iZ = 1; iZ < nBinsZ-1; iZ++){
		fileDatOutput 	<< "\t" << iZ	
					<< "\t" << integratedRInZData[iZ] 							<< "\t" << integratedRInZData[iZ]/sumRinZData * 100 
					<< "\t" << integratedRInZMonteCarlo[iZ] 					<< "\t" << integratedRInZMonteCarlo[iZ]/sumRinZMonteCarlo * 100 
					<< "\t" << integratedRInZData[iZ]- integratedRInZMonteCarlo[iZ] 	<< "\t" << (integratedRInZData[iZ]- integratedRInZMonteCarlo[iZ])/ integratedRInZMonteCarlo[iZ]*100 << endl;
	}
	fileDatOutput 	<< "\t" << "sum" 
				<< "\t" << sumRinZData 					<< "\t" << "100" 
				<< "\t" << sumRinZMonteCarlo 				<< "\t" << "100" 
				<< "\t" << sumRinZData - sumRinZMonteCarlo 	<< "\t" << (sumRinZData-sumRinZMonteCarlo)/sumRinZMonteCarlo *100 << endl;
	fileDatOutput << endl;
	fileDatOutput << "------------------------------------------------------------------------------------------" << endl;           
	// hier ende
	fileDatOutput.close();    


// --------------------------- single plots  ----------------------------
	
	// ---------------------------- R- Distribution -----------------------
	TCanvas * canvasRLin = new TCanvas("canvasRLin","",1200,1000);  // gives the page size  
	canvasRLin->SetLogy(0);
	canvasRLin->SetTopMargin(0.05);            
	canvasRLin->SetLeftMargin(0.13);            
	canvasRLin->cd();
	
	DrawAutoGammaHistosWOLeg( histoConversionPointVsRData, 
							histoConversionPointVsRMonteCarlo, 
							"","R (cm)",textYAxisRHisto,
							kTRUE, 1.2,minYValueRPlot,
							kFALSE,0. ,0.,
							kTRUE, 0.,180.);
	histoTrueConversionPointVSR->SetLineColor(kYellow-7);
	histoTrueConversionPointVSR->SetFillColor(kYellow-7);
	histoTrueConversionPointVSR->Draw("same,hist");
	histoTrueCombinatoricsVsR->SetLineColor(kOrange-5);
	histoTrueCombinatoricsVsR->SetFillColor(kOrange-5);
	histoTrueCombinatoricsVsR->SetFillStyle(3344);
	histoTrueCombinatoricsVsR->Draw("same,hist");
	histoTrueHadVsR->SetLineColor(kGray+3);
	histoTrueHadVsR->SetFillColor(kGray+3);
	histoTrueHadVsR->SetFillStyle(3245);
	histoTrueHadVsR->Draw("same,hist");
	histoTrueDalitzPi0VsR->SetLineColor(kBlue-9);
	histoTrueDalitzPi0VsR->SetFillColor(kBlue-9);
	histoTrueDalitzPi0VsR->SetFillStyle(3244);
	histoTrueDalitzPi0VsR->Draw("same,hist");
	histoTrueDalitzEtaVsR->SetLineColor(kBlue-3);
	histoTrueDalitzEtaVsR->SetFillColor(kBlue-3);
	histoTrueDalitzEtaVsR->SetFillStyle(3002);
	histoTrueDalitzEtaVsR->Draw("same,hist");		
	histoConversionPointVsRData->Draw("same,hist,e");		
	histoConversionPointVsRMonteCarlo->Draw("same,hist,e");
	
	histoConversionPointVsRData->Draw("same,axis");								    

	Bool_t kHadFullR = CheckIfHistoHasEntriesInPlotRange(histoTrueHadVsR, 0,180, minYValueRPlot);
	Bool_t kDalitzEtaFullR = CheckIfHistoHasEntriesInPlotRange(histoTrueDalitzEtaVsR, 0,180, minYValueRPlot);
	Bool_t kDalitzPi0FullR = CheckIfHistoHasEntriesInPlotRange(histoTrueDalitzPi0VsR, 0,180, minYValueRPlot);
	Bool_t kCombFullR =CheckIfHistoHasEntriesInPlotRange(histoTrueCombinatoricsVsR, 0,180, minYValueRPlot);
	Bool_t kSecFullR = CheckIfHistoHasEntriesInPlotRange(histoTrueConversionPointSecVSR, 0,180, minYValueRPlot);
	Int_t nLegendFullR = 3;
	if (kHadFullR) nLegendFullR++;
	if (kCombFullR) nLegendFullR++;
	if (kSecFullR) nLegendFullR++;
	if (kDalitzEtaFullR) nLegendFullR++;
	if (kDalitzPi0FullR) nLegendFullR++;
	Double_t yLegendRFull = 0.027*(nLegendFullR-1);
	
	TLegend* legendRPlot;
	legendRPlot = new TLegend( 0.6,0.94-yLegendRFull,0.8,0.94);
	legendRPlot->SetTextSize(0.025);                        
	legendRPlot->SetFillColor(0);
	legendRPlot->SetBorderSize(0);
	legendRPlot->SetMargin(0.2);
	legendRPlot->AddEntry(histoConversionPointVsRData,("Data"),"l,p");
	legendRPlot->AddEntry(histoConversionPointVsRMonteCarlo,("MC conversion candidates"),"l");
	legendRPlot->AddEntry(histoTrueConversionPointVSR,("MC true conversion"),"f");
	if (kDalitzPi0FullR)legendRPlot->AddEntry(histoTrueDalitzPi0VsR,("MC true #pi^{0} Dalitz"),"f");
	if (kDalitzEtaFullR)legendRPlot->AddEntry(histoTrueDalitzEtaVsR,("MC true #eta Dalitz"),"f");
	if (kCombFullR)legendRPlot->AddEntry(histoTrueCombinatoricsVsR,("MC true combinatorics"),"f");
	if (kHadFullR)legendRPlot->AddEntry(histoTrueHadVsR,("MC true hadronic background"),"f");
	legendRPlot->Draw();
	
	
	
	if(etaCutNumber.CompareTo("2") != 0  && etaCutNumber.CompareTo("3") != 0) DrawIndividualTextSlicesR (vectorLinPlot , 0.025,"lin");
	
	if(!thesis)DrawAliceLogoPerformance(floatLocationRightUpR[0],floatLocationRightUpR[1],floatLocationRightUpR[2],floatLocationRightUpR[3],0.0001, textDate,collisionSystem, textGenerator, textPeriod,1200,1000); 
	canvasRLin->Update();
	canvasRLin->SaveAs(Form("%s/R_distribution_lin.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasRLin;          
	
	
	TCanvas * canvasRLog = new TCanvas("canvasRLog","",1200,1000);  // gives the page size
	canvasRLog->cd();
	canvasRLog->SetLeftMargin(0.13);
	//canvasRLog->SetTopMargin(0.05);            
	canvasRLog->SetLogy(1);
	
	DrawAutoGammaHistosWOLeg( histoConversionPointVsRData, 
							histoConversionPointVsRMonteCarlo, 
							"","R (cm)",textYAxisRHisto,
							kTRUE, 1.5,minYValueRPlot,
							kFALSE,0. ,0.,
							kTRUE, 0.,180.);		
	histoTrueConversionPointVSR->Draw("same,hist");
	histoTrueCombinatoricsVsR->Draw("same,hist");
	histoTrueHadVsR->Draw("same,hist");
	histoTrueDalitzPi0VsR->Draw("same,hist");
	histoTrueDalitzEtaVsR->Draw("same,hist");		
	histoConversionPointVsRData->Draw("same,hist,e");		
	histoConversionPointVsRMonteCarlo->Draw("same,hist,e");
	histoConversionPointVsRData->Draw("same,axis");		
	
	if(etaCutNumber.CompareTo("2") != 0  && etaCutNumber.CompareTo("3") != 0) DrawIndividualTextSlicesR (vectorLogPlot , 0.025,"log");
	
	legendRPlot->Draw();
	if(!thesis)DrawAliceLogoPerformance(floatLocationRightUpR[0],floatLocationRightUpR[1],floatLocationRightUpR[2],floatLocationRightUpR[3],0.0001, textDate,collisionSystem, textGenerator, textPeriod,1200,1000); 
	
	canvasRLog->Update();
	canvasRLog->SaveAs(Form("%s/R_distribution_log.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasRLog;
		
	TCanvas * canvasREnlargedSDDSSDLog = new TCanvas("canvasREnlargedSDDSSDLog","",1200,1000);  // gives the page size
	canvasREnlargedSDDSSDLog->cd();
	canvasREnlargedSDDSSDLog->SetLeftMargin(0.13);
//	canvasREnlargedSDDSSDLog->SetTopMargin(0.05);            
	canvasREnlargedSDDSSDLog->SetLogy(1);
	
	DrawAutoGammaHistosWOLeg( histoConversionPointVsRData, 
							histoConversionPointVsRMonteCarlo, 
							"","R (cm)",textYAxisRHisto,
							kTRUE, 2.4,minYValueRPlot,
							kFALSE,1.e-2 ,minYValueRPlotSPD,
							kTRUE, 13.,60.);
	if( etaCutNumber.CompareTo("2") == 0 || etaCutNumber.CompareTo("3") == 0 ){
		histoConversionPointVsRData->GetYaxis()->SetRangeUser(minYValueRPlotSPD,6e-3);
	} else {
		histoConversionPointVsRData->GetYaxis()->SetRangeUser(minYValueRPlotSPD,2e-3);
	}
	histoConversionPointVsRData->Draw("hist");
	histoTrueConversionPointVSR->Draw("same,hist");
	histoTrueCombinatoricsVsR->Draw("same,hist");
	histoTrueHadVsR->Draw("same,hist");
	histoTrueDalitzPi0VsR->Draw("same,hist");
	histoTrueDalitzEtaVsR->Draw("same,hist");		
	histoConversionPointVsRData->Draw("same,hist,e");		
	histoConversionPointVsRMonteCarlo->Draw("same,hist,e");
	histoConversionPointVsRData->Draw("same,axis");											

	Bool_t kHadSDDR = CheckIfHistoHasEntriesInPlotRange(histoTrueHadVsR, 13.,60., minYValueRPlotSPD);
	Bool_t kDalitzEtaSDDR = CheckIfHistoHasEntriesInPlotRange(histoTrueDalitzEtaVsR, 13.,60.0, minYValueRPlotSPD);
	Bool_t kDalitzPi0SDDR = CheckIfHistoHasEntriesInPlotRange(histoTrueDalitzPi0VsR, 13.,60., minYValueRPlotSPD);
	Bool_t kCombSDDR =CheckIfHistoHasEntriesInPlotRange(histoTrueCombinatoricsVsR, 13.,60., minYValueRPlotSPD);
	Bool_t kSecSDDR = CheckIfHistoHasEntriesInPlotRange(histoTrueConversionPointSecVSR, 13.,60., minYValueRPlotSPD);
	Int_t nLegendSDDR = 3;
	if (kHadSDDR) nLegendSDDR++;
	if (kCombSDDR) nLegendSDDR++;
	if (kSecSDDR) nLegendSDDR++;
	if (kDalitzEtaSDDR) nLegendSDDR++;
	if (kDalitzPi0SDDR) nLegendSDDR++;
	Double_t yLegendRSDD = 0.027*(nLegendSDDR-1);
	
	TLegend* legendRPlotSDD;
	if (etaCutNumber.CompareTo("2") == 0 || etaCutNumber.CompareTo("3") == 0 ){
		legendRPlotSDD = new TLegend( 0.35,0.96-yLegendRSDD,0.55,0.96);
	} else {
		legendRPlotSDD = new TLegend( 0.16,0.96-yLegendRSDD,0.36,0.96);
	}
	legendRPlotSDD->SetTextSize(0.025);                        
	legendRPlotSDD->SetFillColor(0);
	legendRPlotSDD->SetBorderSize(0);
	legendRPlotSDD->SetMargin(0.2);
	legendRPlotSDD->AddEntry(histoConversionPointVsRData,("Data"),"l,p");
	legendRPlotSDD->AddEntry(histoConversionPointVsRMonteCarlo,("MC conversion candidates"),"l");
	legendRPlotSDD->AddEntry(histoTrueConversionPointVSR,("MC true conversion"),"f");
	if (kDalitzPi0SDDR)legendRPlotSDD->AddEntry(histoTrueDalitzPi0VsR,("MC true #pi^{0} Dalitz"),"f");
	if (kDalitzEtaSDDR)legendRPlotSDD->AddEntry(histoTrueDalitzEtaVsR,("MC true #eta Dalitz"),"f");
	if (kCombSDDR)legendRPlotSDD->AddEntry(histoTrueCombinatoricsVsR,("MC true combinatorics"),"f");
	if (kHadSDDR)legendRPlotSDD->AddEntry(histoTrueHadVsR,("MC true hadronic background"),"f");
	legendRPlotSDD->Draw();

	if(etaCutNumber.CompareTo("2") == 0 || etaCutNumber.CompareTo("3") == 0 ){
		if(!thesis)DrawAliceLogoPerformance(floatLocationLeftUp[0],floatLocationLeftUp[1]+0.05,floatLocationLeftUp[2],floatLocationLeftUp[3],0.0001, textDate,collisionSystem, textGenerator, textPeriod,1200,1000);     
	} else {
		if(!thesis)DrawAliceLogoPerformance(floatLocationRightUpR[0],floatLocationRightUpR[1]+0.25,floatLocationRightUpR[2],floatLocationRightUpR[3],0.0001, textDate,collisionSystem, textGenerator, textPeriod,1200,1000); 
	}
	
	canvasREnlargedSDDSSDLog->Update();
	canvasREnlargedSDDSSDLog->SaveAs(Form("%s/R_distributionSDDSSDEnlarged_log.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasREnlargedSDDSSDLog;
	
	
	TCanvas * canvasREnlargedTPCLog = new TCanvas("canvasREnlargedTPCLog","",1200,1000);  // gives the page size
	canvasREnlargedTPCLog->cd();
	canvasREnlargedTPCLog->SetLeftMargin(0.13);
//	canvasREnlargedTPCLog->SetTopMargin(0.05);            
	canvasREnlargedTPCLog->SetLogy(1);
	
	DrawAutoGammaHistosWOLeg( histoConversionPointVsRData, 
							histoConversionPointVsRMonteCarlo, 
							"","R (cm)",textYAxisRHisto,
							kTRUE, 1.,minYValueRPlot,
							kFALSE,0. ,0.,
							kTRUE, 60.,180.);
	histoTrueConversionPointVSR->Draw("same,hist");
	histoTrueCombinatoricsVsR->Draw("same,hist");
	histoTrueHadVsR->Draw("same,hist");
	histoTrueDalitzPi0VsR->Draw("same,hist");
	histoTrueDalitzEtaVsR->Draw("same,hist");		
	histoConversionPointVsRData->Draw("same,hist,e");		
	histoConversionPointVsRMonteCarlo->Draw("same,hist,e");
	histoConversionPointVsRData->Draw("same,axis");											
	
	Bool_t kHadTPCR = CheckIfHistoHasEntriesInPlotRange(histoTrueHadVsR, 60.,180., minYValueRPlot);
	Bool_t kDalitzEtaTPCR = CheckIfHistoHasEntriesInPlotRange(histoTrueDalitzEtaVsR, 60.,180., minYValueRPlot);
	Bool_t kDalitzPi0TPCR = CheckIfHistoHasEntriesInPlotRange(histoTrueDalitzPi0VsR, 60.,180., minYValueRPlot);
	Bool_t kCombTPCR =CheckIfHistoHasEntriesInPlotRange(histoTrueCombinatoricsVsR, 60.,180., minYValueRPlot);
	Bool_t kSecTPCR = CheckIfHistoHasEntriesInPlotRange(histoTrueConversionPointSecVSR, 60.,180., minYValueRPlot);
	Int_t nLegendTPCR = 3;
	if (kHadTPCR) nLegendTPCR++;
	if (kCombTPCR) nLegendTPCR++;
	if (kSecTPCR) nLegendTPCR++;
	if (kDalitzEtaTPCR) nLegendTPCR++;
	if (kDalitzPi0TPCR) nLegendTPCR++;
	Double_t yLegendRTPC = 0.027*(nLegendTPCR-1);
	
	TLegend* legendRPlotTPC;
	legendRPlotTPC = new TLegend( 0.6,0.96-yLegendRTPC,0.8,0.96);
	legendRPlotTPC->SetTextSize(0.025);                        
	legendRPlotTPC->SetFillColor(0);
	legendRPlotTPC->SetBorderSize(0);
	legendRPlotTPC->SetMargin(0.2);
	legendRPlotTPC->AddEntry(histoConversionPointVsRData,("Data"),"l,p");
	legendRPlotTPC->AddEntry(histoConversionPointVsRMonteCarlo,("MC conversion candidates"),"l");
	legendRPlotTPC->AddEntry(histoTrueConversionPointVSR,("MC true conversion"),"f");
	if (kDalitzPi0TPCR)legendRPlotTPC->AddEntry(histoTrueDalitzPi0VsR,("MC true #pi^{0} Dalitz"),"f");
	if (kDalitzEtaTPCR)legendRPlotTPC->AddEntry(histoTrueDalitzEtaVsR,("MC true #eta Dalitz"),"f");
	if (kCombTPCR)legendRPlotTPC->AddEntry(histoTrueCombinatoricsVsR,("MC true combinatorics"),"f");
	if (kHadTPCR)legendRPlotTPC->AddEntry(histoTrueHadVsR,("MC true hadronic background"),"f");
	legendRPlotTPC->Draw();

	
	if(!thesis)DrawAliceLogoPerformance(floatLocationRightUpR[0],floatLocationRightUpR[1],floatLocationRightUpR[2],floatLocationRightUpR[3],0.0001, textDate,collisionSystem, textGenerator, textPeriod,1200,1000); 
	
	canvasREnlargedTPCLog->Update();
	canvasREnlargedTPCLog->SaveAs(Form("%s/R_distributionTPCEnlarged_log.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasREnlargedTPCLog;

	
	TCanvas * canvasREnlargedSPDLog = new TCanvas("canvasREnlargedSPDLog","",1200,1000);  // gives the page size
	canvasREnlargedSPDLog->cd();
	canvasREnlargedSPDLog->SetLeftMargin(0.13);
//	canvasREnlargedSPDLog->SetTopMargin(0.05);            
	canvasREnlargedSPDLog->SetLogy(1);
	
	DrawAutoGammaHistosWOLeg( histoConversionPointVsRData, 
							histoConversionPointVsRMonteCarlo, 
							"","R (cm)",textYAxisRHisto,
							kTRUE, 1.5,minYValueRPlotSPD,
							kFALSE,3.e-3 ,minYValueRPlotSPD,
							kTRUE, 0.,13.);
	if( etaCutNumber.CompareTo("2") == 0 || etaCutNumber.CompareTo("3") == 0 ){
		histoConversionPointVsRData->GetYaxis()->SetRangeUser(minYValueRPlotSPD,5e-3);
	} else {
		histoConversionPointVsRData->GetYaxis()->SetRangeUser(minYValueRPlotSPD,1e-2);
	}
	histoConversionPointVsRData->Draw("hist");
	histoTrueConversionPointVSR->Draw("same,hist");
	histoTrueCombinatoricsVsR->Draw("same,hist");
	histoTrueHadVsR->Draw("same,hist");
	histoTrueDalitzPi0VsR->Draw("same,hist");
	histoTrueDalitzEtaVsR->Draw("same,hist");		
	histoConversionPointVsRData->Draw("same,hist,e");		
	histoConversionPointVsRMonteCarlo->Draw("same,hist,e");
	histoConversionPointVsRData->Draw("same,axis");		
	
	Bool_t kHadSPDR = CheckIfHistoHasEntriesInPlotRange(histoTrueHadVsR, 0.,13., minYValueRPlotSPD);
	Bool_t kDalitzEtaSPDR = CheckIfHistoHasEntriesInPlotRange(histoTrueDalitzEtaVsR, 0.,13., minYValueRPlotSPD);
	Bool_t kDalitzPi0SPDR = CheckIfHistoHasEntriesInPlotRange(histoTrueDalitzPi0VsR,0.,13., minYValueRPlotSPD);
	Bool_t kCombSPDR =CheckIfHistoHasEntriesInPlotRange(histoTrueCombinatoricsVsR, 0.,13., minYValueRPlotSPD);
	Bool_t kSecSPDR = CheckIfHistoHasEntriesInPlotRange(histoTrueConversionPointSecVSR, 0.,13., minYValueRPlotSPD);
	Int_t nLegendSPDR = 3;
	if (kHadSPDR) nLegendSPDR++;
	if (kCombSPDR) nLegendSPDR++;
	if (kSecSPDR) nLegendSPDR++;
	if (kDalitzEtaSPDR) nLegendSPDR++;
	if (kDalitzPi0SPDR) nLegendSPDR++;
	Double_t yLegendRSPD = 0.027*(nLegendSPDR-1);
	
	TLegend* legendRPlotSPD;
	legendRPlotSPD = new TLegend( 0.6,0.96-yLegendRSPD,0.8,0.96);
	legendRPlotSPD->SetTextSize(0.025);                        
	legendRPlotSPD->SetFillColor(0);
	legendRPlotSPD->SetBorderSize(0);
	legendRPlotSPD->SetMargin(0.2);
	legendRPlotSPD->AddEntry(histoConversionPointVsRData,("Data"),"l,p");
	legendRPlotSPD->AddEntry(histoConversionPointVsRMonteCarlo,("MC conversion candidates"),"l");
	legendRPlotSPD->AddEntry(histoTrueConversionPointVSR,("MC true conversion"),"f");
	if (kDalitzPi0SPDR)legendRPlotSPD->AddEntry(histoTrueDalitzPi0VsR,("MC true #pi^{0} Dalitz"),"f");
	if (kDalitzEtaSPDR)legendRPlotSPD->AddEntry(histoTrueDalitzEtaVsR,("MC true #eta Dalitz"),"f");
	if (kCombSPDR)legendRPlotSPD->AddEntry(histoTrueCombinatoricsVsR,("MC true combinatorics"),"f");
	if (kHadSPDR)legendRPlotSPD->AddEntry(histoTrueHadVsR,("MC true hadronic background"),"f");
	legendRPlotSPD->Draw();

	if(!thesis)DrawAliceLogoPerformance(floatLocationLeftUp[0],floatLocationLeftUp[1],floatLocationLeftUp[2],floatLocationLeftUp[3],0.0001, textDate,collisionSystem, textGenerator, textPeriod,1200,1000);     
	
	canvasREnlargedSPDLog->Update();
	canvasREnlargedSPDLog->SaveAs(Form("%s/R_distributionSPDEnlarged_log.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasREnlargedSPDLog;
	delete legendRPlot;
	
	//---------------------------------Eta distribution -------------------------------------------------------
	TCanvas * canvasEtaLog = new TCanvas("canvasEtaLog","",1200,1000);  // gives the page size
	canvasEtaLog->cd();
	canvasEtaLog->SetLeftMargin(0.13);
	//canvasEtaLog->SetTopMargin(0.05);            
	canvasEtaLog->SetLogy(1);
	
	DrawAutoGammaHistosWOLeg( histoConversionPointVsEtaData, 
							histoConversionPointVsEtaMonteCarlo, 
							"","#eta",textYAxisEtaHisto,
							kTRUE, 3.,minYValueEtaPlot,
							kFALSE,0. ,0.,
							kTRUE,-maxEta-0.2 ,maxEta+0.2);		
	histoTrueConversionPointPrimVsEta->SetLineColor(kYellow-7);
	histoTrueConversionPointPrimVsEta->SetFillColor(kYellow-7);
	histoTrueConversionPointPrimVsEta->Draw("same,hist");
	histoTrueConversionPointSecVSEta->SetLineColor(kSpring+2);
	histoTrueConversionPointSecVSEta->SetFillColor(kSpring+2);
	histoTrueConversionPointSecVSEta->Draw("same,hist");
	histoTrueCombinatoricsVsEta->SetLineColor(kOrange-5);
	histoTrueCombinatoricsVsEta->SetFillColor(kOrange-5);
	histoTrueCombinatoricsVsEta->SetFillStyle(3344);
	histoTrueCombinatoricsVsEta->Draw("same,hist");
	histoTrueHadVsEta->SetLineColor(kGray+3);
	histoTrueHadVsEta->SetFillColor(kGray+3);
	histoTrueHadVsEta->SetFillStyle(3245);
	histoTrueHadVsEta->Draw("same,hist");
	histoTrueDalitzPi0VsEta->SetLineColor(kBlue-9);
	histoTrueDalitzPi0VsEta->SetFillColor(kBlue-9);
	histoTrueDalitzPi0VsEta->SetFillStyle(3244);
	histoTrueDalitzPi0VsEta->Draw("same,hist");
	histoTrueDalitzEtaVsEta->SetLineColor(kBlue-3);
	histoTrueDalitzEtaVsEta->SetFillColor(kBlue-3);
	histoTrueDalitzEtaVsEta->SetFillStyle(3002);
	histoTrueDalitzEtaVsEta->Draw("same,hist");		
	histoConversionPointVsEtaData->Draw("same,hist,e");		
	histoConversionPointVsEtaMonteCarlo->Draw("same,hist,e");
	histoConversionPointVsEtaData->Draw("same,axis");		
	
	
	Bool_t kHadEta = CheckIfHistoHasEntriesInPlotRange(histoTrueHadVsEta, -maxEta-0.05 ,maxEta+0.05, minYValueEtaPlot);
	Bool_t kDalitzEtaEta = CheckIfHistoHasEntriesInPlotRange(histoTrueDalitzEtaVsEta, -maxEta-0.05 ,maxEta+0.05, minYValueEtaPlot);
	Bool_t kDalitzPi0Eta = CheckIfHistoHasEntriesInPlotRange(histoTrueDalitzPi0VsEta,-maxEta-0.05 ,maxEta+0.05, minYValueEtaPlot);
	Bool_t kCombEta =CheckIfHistoHasEntriesInPlotRange(histoTrueCombinatoricsVsEta, -maxEta-0.05 ,maxEta+0.05, minYValueEtaPlot);
	Bool_t kSecEta = CheckIfHistoHasEntriesInPlotRange(histoTrueConversionPointSecVSEta, -maxEta-0.05 ,maxEta+0.05, minYValueEtaPlot);
	Int_t nLegendEta = 3;
	if (kHadEta) nLegendEta++;
	if (kCombEta) nLegendEta++;
	if (kSecEta) nLegendEta++;
	if (kDalitzEtaEta) nLegendEta++;
	if (kDalitzPi0Eta) nLegendEta++;
	Double_t yLegendEta;
	if (nLegendEta%2 == 1){
		yLegendEta = (0.027*(nLegendEta+1))/2;
	} else {
		yLegendEta = 0.027*nLegendEta/2;
	}
	
	TLegend* legendEtaPlot;
	legendEtaPlot = new TLegend( 0.3,0.96-yLegendEta,0.96,0.96);
	legendEtaPlot->SetTextSize(0.025);                        
	legendEtaPlot->SetFillColor(0);
	legendEtaPlot->SetBorderSize(0);
	legendEtaPlot->SetMargin(0.13);
	legendEtaPlot->SetNColumns(2);
	legendEtaPlot->AddEntry(histoConversionPointVsEtaData,("Data"),"l,p");
	legendEtaPlot->AddEntry(histoConversionPointVsEtaMonteCarlo,("MC conversion candidates"),"l");
	legendEtaPlot->AddEntry(histoTrueConversionPointPrimVsEta,("MC true primary conversion"),"f");
	if (kSecEta)legendEtaPlot->AddEntry(histoTrueConversionPointSecVSEta,("MC true secondary conversion"),"f");
	if (kDalitzPi0Eta)legendEtaPlot->AddEntry(histoTrueDalitzPi0VsEta,("MC true #pi^{0} Dalitz"),"f");
	if (kDalitzEtaEta)legendEtaPlot->AddEntry(histoTrueDalitzEtaVsEta,("MC true #eta Dalitz"),"f");
	if (kCombEta)legendEtaPlot->AddEntry(histoTrueCombinatoricsVsEta,("MC true combinatorics"),"f");
	if (kHadEta)legendEtaPlot->AddEntry(histoTrueHadVsEta,("MC true hadronic background"),"f");
	legendEtaPlot->Draw();

	if(!thesis)DrawAliceLogoPerformance(floatLocationLeftUpEta[0],floatLocationLeftUpEta[1],floatLocationLeftUpEta[2],floatLocationLeftUpEta[3],0.0001, textDate,collisionSystem, textGenerator, textPeriod,1200,1000); 
	
	canvasEtaLog->Update();
	canvasEtaLog->SaveAs(Form("%s/Eta_distribution_log.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasEtaLog;

	
	// -----------------------  2 dim Plots -----------------------------------
	
// 	PlotStandard2D( histoMCConversionPointVsZRMonteCarlo , Form("%s/ZR_MCdistribution.%s",outputDirectory.Data(),suffix.Data()), "", "Z (cm)", "R (cm)", kTRUE, 0, 100., kTRUE, -150., 150.,kTRUE,  5.e-8, 1.3e-4, 0, 1, floatLocationRightUp2DXY,500,440,rightMargin-0.01, leftMargin,0.08, 0.02,"");
// 	delete histoMCConversionPointVsZRMonteCarlo;
// 	PlotStandard2D( histoMCConversionPointVsXYMonteCarlo , Form("%s/XY_MCdistribution.%s",outputDirectory.Data(),suffix.Data()), "", "X (cm)", "Y (cm)", kTRUE, -100, 100., kTRUE, -100., 100.,kTRUE,  5.e-8, 1.3e-4, 0, 1, floatLocationRightUp2DXY,500,440,rightMargin-0.01, leftMargin,0.08, 0.02,"");
// 	delete histoMCConversionPointVsXYMonteCarlo;	
	
	
	
	TCanvas * canvas2DimensionZR = new TCanvas("canvas2DimensionZR","",500,440);  // gives the page size
	canvas2DimensionZR->SetLogz(1);     
	canvas2DimensionZR->SetTopMargin(0.02);                                              
	canvas2DimensionZR->SetRightMargin(rightMargin-0.02);                                            
	canvas2DimensionZR->SetLeftMargin(leftMargin);
	canvas2DimensionZR->SetBottomMargin(0.08);
	canvas2DimensionZR->cd();
	histoConversionPointVsZRPlotData->GetZaxis()->SetRangeUser(1,histoConversionPointVsZRPlotData->GetMaximum());
	DrawAutoGammaHisto2D(  histoConversionPointVsZRPlotData,
									"", "Z (cm)", "R (cm)", "",
									kTRUE, 0., 100.,
									kTRUE, -150., 150.);
	latexEtaRange->Draw();
	if(etaCutNumber.CompareTo("2") != 0  && etaCutNumber.CompareTo("3") != 0) DrawStructureZRNew();
	DrawLabelsEvents(floatLocationRight[0],floatLocationRight[1],floatLocationRight[2], 0.00, collisionSystem, "", textPeriod);
// 	if(!thesis) DrawAliceLogoPerformance2D(floatLocationRightUp2DZR[0], floatLocationRightUp2DZR[1], floatLocationRightUp2DZR[2], floatLocationRightUp2DZR[3], 0.02, textDate, collisionSystem, "", textPeriod, 500,440);
	
	canvas2DimensionZR->Update();
	canvas2DimensionZR->SaveAs(Form("%s/ZR_distribution.%s",outputDirectory.Data(),suffix.Data()));

	if (suffix.CompareTo("png")==0){
		TCanvas * canvas3dPlotZR = new TCanvas("canvas3dPlotZR","",800,720);  // gives the page size
		canvas3dPlotZR->SetLogz(1);     		
		canvas3dPlotZR->cd();
		canvas3dPlotZR->SetTopMargin(0.00);                                              
		canvas3dPlotZR->SetRightMargin(0);                                            
		canvas3dPlotZR->SetLeftMargin(0);
		canvas3dPlotZR->SetBottomMargin(0.0);

		histoConversionPointVsZRPlotData->GetZaxis()->SetRangeUser(1,histoConversionPointVsRPhiData2->GetMaximum());
		DrawAutoGammaHisto2D(  histoConversionPointVsZRPlotData,
										"", "Z (cm)", "R (cm)", "",
										kTRUE, 0., 180.,
										kTRUE, -200., 200.);
		if (suffix.CompareTo("png")==0) histoConversionPointVsZRPlotData->Draw("col");
		canvas3dPlotZR->Update();
		canvas3dPlotZR->SaveAs(Form("%s/ZR_distribution_plain.%s",outputDirectory.Data(),suffix.Data()));
	}
	delete canvas2DimensionZR;
	
	TCanvas * canvas2DimensionRPhi = NULL;
	if (suffix.CompareTo("png")==0) canvas2DimensionRPhi =new TCanvas("canvas2DimensionRPhi","",695,800);  // gives the page size
	else canvas2DimensionRPhi =new TCanvas("canvas2DimensionRPhi","",500,440);  // gives the page size
	canvas2DimensionRPhi->SetLogz(1);     
	if (suffix.CompareTo("png")==0){
		canvas2DimensionRPhi->SetTopMargin(0.0);                                              
		canvas2DimensionRPhi->SetRightMargin(0.0);                                            
		canvas2DimensionRPhi->SetLeftMargin(0.0);
		canvas2DimensionRPhi->SetBottomMargin(0.0);
	} else {
		canvas2DimensionRPhi->SetTopMargin(0.02);                                              
		canvas2DimensionRPhi->SetRightMargin(rightMargin-0.02);                                            
		canvas2DimensionRPhi->SetLeftMargin(leftMargin);
		canvas2DimensionRPhi->SetBottomMargin(0.08);		
	}	
	canvas2DimensionRPhi->cd();
	histoConversionPointVsRPhiData2->GetZaxis()->SetRangeUser(1,histoConversionPointVsRPhiData2->GetMaximum());
	DrawAutoGammaHisto2D(  histoConversionPointVsRPhiData2,
									"", "R (cm)", "#varphi (cm)", "",
									kFALSE, 0., 100.,
									kTRUE, 0., 180.);
	if (suffix.CompareTo("png")==0) histoConversionPointVsRPhiData2->Draw("col");
	canvas2DimensionRPhi->Update();
	canvas2DimensionRPhi->SaveAs(Form("%s/RPhi_distribution_plain.%s",outputDirectory.Data(),suffix.Data()));

	TCanvas * canvas2DimensionZPhi = NULL;  // gives the page size
	if (suffix.CompareTo("png")==0) canvas2DimensionZPhi =new TCanvas("canvas2DimensionZPhi","",800,800);  // gives the page size
	else canvas2DimensionZPhi =new TCanvas("canvas2DimensionZPhi","",500,440);  // gives the page size

	canvas2DimensionZPhi->SetLogz(1);     
	if (suffix.CompareTo("png")==0){
		canvas2DimensionZPhi->SetTopMargin(0.0);                                              
		canvas2DimensionZPhi->SetRightMargin(0.0);                                            
		canvas2DimensionZPhi->SetLeftMargin(0.0);
		canvas2DimensionZPhi->SetBottomMargin(0.0);
	} else {
		canvas2DimensionZPhi->SetTopMargin(0.02);                                              
		canvas2DimensionZPhi->SetRightMargin(rightMargin-0.02);                                            
		canvas2DimensionZPhi->SetLeftMargin(leftMargin);
		canvas2DimensionZPhi->SetBottomMargin(0.08);		
	}	
	canvas2DimensionZPhi->cd();
	histoConversionPointVsZPhiData2->GetZaxis()->SetRangeUser(1,histoConversionPointVsRPhiData2->GetMaximum());
	DrawAutoGammaHisto2D(  histoConversionPointVsZPhiData2,
									"", "Z (cm)", "#varphi (cm)", "",
									kFALSE, 0., 100.,
									kTRUE, -200., 200.);
	if (suffix.CompareTo("png")==0) histoConversionPointVsZPhiData2->Draw("col");
	canvas2DimensionZPhi->Update();
	canvas2DimensionZPhi->SaveAs(Form("%s/ZPhi_distribution_plain.%s",outputDirectory.Data(),suffix.Data()));
	
	
	TCanvas * canvas2DimensionXY = new TCanvas("canvas2DimensionXY","",500,440);  // gives the page size
	canvas2DimensionXY->SetLogz(1);
	canvas2DimensionXY->SetTopMargin(0.02);                                              
	canvas2DimensionXY->SetRightMargin(rightMargin);                                            
	canvas2DimensionXY->SetLeftMargin(leftMargin);
	canvas2DimensionXY->SetBottomMargin(0.08);
	canvas2DimensionXY->cd();
	histoConversionPointVsXYPlotData->GetZaxis()->SetRangeUser(1,histoConversionPointVsXYPlotData->GetMaximum());
	DrawAutoGammaHisto2D(   histoConversionPointVsXYPlotData,
									"", "X (cm)", "Y (cm)", "",
									kTRUE, -250., 250.,
									kTRUE, -250., 250.);
	if(etaCutNumber.CompareTo("2") != 0  && etaCutNumber.CompareTo("3") != 0) DrawStructureNew();
	latexEtaRange->SetTextColor(kWhite);
	latexEtaRange->Draw();

	if(!thesis) 	DrawAliceLogoPerformance2D(floatLocationRightUp2DXY[0], floatLocationRightUp2DXY[1], floatLocationRightUp2DXY[2], floatLocationRightUp2DXY[3], 0.02, textDate, collisionSystem, "", textPeriod, 500,440);
	canvas2DimensionXY->Update();
	canvas2DimensionXY->SaveAs(Form("%s/XY_distribution.%s",outputDirectory.Data(),suffix.Data()));
	delete canvas2DimensionXY;
						    
	TH1D* histoDiffXYDistribution;
	TH1D* histoDiffZRDistribution;
	TCanvas * canvas2DimensionDiffXY = new TCanvas("canvas2DimensionDiffXY","",500,440);  // gives the page size
	canvas2DimensionDiffXY->SetLogz(1);
	canvas2DimensionDiffXY->SetTopMargin(0.02);                                              
	canvas2DimensionDiffXY->SetRightMargin(rightMargin);                                            
	canvas2DimensionDiffXY->SetLeftMargin(leftMargin);
	canvas2DimensionDiffXY->SetBottomMargin(0.08);
	canvas2DimensionDiffXY->cd();
	histoDiffXYDistribution = (TH1D*)histoConversionPointVsXYData2->Clone();
	histoDiffXYDistribution->Divide(histoConversionPointVsXYData2,histoConversionPointVsXYMonteCarlo2,1.,1.,"B");
	histoDiffXYDistribution->GetZaxis()->SetRangeUser(0.05,15.);
	histoDiffXYDistribution->Draw("colz");     
	if(!thesis)DrawAliceLogo(floatLocationRightUp2D [0],floatLocationRightUp2D[1],floatLocationRightUp2D[2],floatLocationRightUp2D[3],500,440);
	canvas2DimensionDiffXY->Update();
	canvas2DimensionDiffXY->SaveAs(Form("%s/DiffXYDistributions.%s",outputDirectory.Data(),suffix.Data()));
	delete canvas2DimensionDiffXY;

	TCanvas * canvas2DimensionDiffZR = new TCanvas("canvas2DimensionDiffZR","",500,440);  // gives the page size
	canvas2DimensionDiffZR->SetLogz(1);     
	canvas2DimensionDiffZR->SetTopMargin(0.02);                                              
	canvas2DimensionDiffZR->SetRightMargin(rightMargin);                                            
	canvas2DimensionDiffZR->SetLeftMargin(leftMargin);
	canvas2DimensionDiffZR->SetBottomMargin(0.08);
	canvas2DimensionDiffZR->cd();

	histoDiffZRDistribution = (TH1D*)histoConversionPointVsZRData2->Clone();
	histoDiffZRDistribution->Divide(histoConversionPointVsZRData2,histoConversionPointVsZRMonteCarlo2,1.,1.,"B");
	histoDiffZRDistribution->GetZaxis()->SetRangeUser(0.05,15.);
	histoDiffZRDistribution->Draw("colz");     
	if(!thesis)DrawAliceLogo(floatLocationRightUp2D [0],floatLocationRightUp2D[1],floatLocationRightUp2D[2],floatLocationRightUp2D[3],500,440);
	canvas2DimensionDiffZR->Update();
	canvas2DimensionDiffZR->SaveAs(Form("%s/DiffZRDistributions.%s",outputDirectory.Data(),suffix.Data()));
	delete canvas2DimensionDiffZR;

	
	//-------------------- Z - Distribution ------------------------------
	
	
	TCanvas * canvasZLin = new TCanvas("canvasZLin","",1200,1000);  // gives the page size  
	canvasZLin->SetLogy(0);
	canvasZLin->SetTopMargin(0.04);            
	canvasZLin->SetLeftMargin(0.13);            
	canvasZLin->cd();

	
	DrawAutoGammaHistosWOLeg( histoConversionPointVsZData, 
							histoConversionPointVsZMonteCarlo, 
							"","Z (cm)",textYAxisZHisto,
							kTRUE, maxYFactorZPlotlin*1.5 ,minYValueZPlot,
							kFALSE,0. ,0.,
							kTRUE,  -201.,201.);
	histoTrueConversionPointVSZ->SetLineColor(kYellow-7);
	histoTrueConversionPointVSZ->SetFillColor(kYellow-7);
	histoTrueConversionPointVSZ->Draw("same,hist");
	histoTrueCombinatoricsVsZ->SetLineColor(kOrange-5);
	histoTrueCombinatoricsVsZ->SetFillColor(kOrange-5);
	histoTrueCombinatoricsVsZ->Draw("same,hist");
	histoTrueCombinatoricsVsZ->SetFillStyle(3344);
	histoTrueHadVsZ->SetLineColor(kGray+3);
	histoTrueHadVsZ->SetFillColor(kGray+3);
	histoTrueHadVsZ->SetFillStyle(3245);
	histoTrueHadVsZ->Draw("same,hist");
	histoTrueDalitzPi0VsZ->SetLineColor(kBlue-9);
	histoTrueDalitzPi0VsZ->SetFillColor(kBlue-9);
	histoTrueDalitzPi0VsZ->SetFillStyle(3244);
	histoTrueDalitzPi0VsZ->Draw("same,hist");
	histoTrueDalitzEtaVsZ->SetLineColor(kBlue-3);
	histoTrueDalitzEtaVsZ->SetFillColor(kBlue-3);
	histoTrueDalitzEtaVsZ->SetFillStyle(3002);
	histoTrueDalitzEtaVsZ->Draw("same,hist");		
	histoConversionPointVsZData->Draw("same,hist,e");		
	histoConversionPointVsZMonteCarlo->Draw("same,hist,e");
	histoConversionPointVsZData->Draw("same,axis");								    
	
	Bool_t kHadZ = CheckIfHistoHasEntriesInPlotRange(histoTrueHadVsZ, -201.,201., minYValueZPlot);
	Bool_t kDalitzEtaZ = CheckIfHistoHasEntriesInPlotRange(histoTrueDalitzEtaVsZ, -201.,201., minYValueZPlot);
	Bool_t kDalitzPi0Z = CheckIfHistoHasEntriesInPlotRange(histoTrueDalitzPi0VsZ,-201.,201., minYValueZPlot);
	Bool_t kCombZ =CheckIfHistoHasEntriesInPlotRange(histoTrueCombinatoricsVsZ,-201.,201., minYValueZPlot);
	Bool_t kSecZ = CheckIfHistoHasEntriesInPlotRange(histoTrueConversionPointSecVSZ, -201.,201., minYValueZPlot);
	Int_t nLegendZ = 3;
	if (kHadZ) nLegendZ++;
	if (kCombZ) nLegendZ++;
	if (kSecZ) nLegendZ++;
	if (kDalitzEtaZ) nLegendZ++;
	if (kDalitzPi0Z) nLegendZ++;
	Double_t yLegendZ = 0.027*(nLegendZ-1);
	
	TLegend* legendZPlot;
	legendZPlot = new TLegend( 0.145,0.94-yLegendZ,0.345,0.94);
	legendZPlot->SetTextSize(0.025);                        
	legendZPlot->SetFillColor(0);
	legendZPlot->SetBorderSize(0);
	legendZPlot->SetMargin(0.2);
	legendZPlot->AddEntry(histoConversionPointVsZData,("Data"),"l,p");
	legendZPlot->AddEntry(histoConversionPointVsZMonteCarlo,("MC conversion candidates"),"l");
	legendZPlot->AddEntry(histoTrueConversionPointVSZ,("MC true conversion"),"f");
	if (kDalitzPi0Z)legendZPlot->AddEntry(histoTrueDalitzPi0VsZ,("MC true #pi^{0} Dalitz"),"f");
	if (kDalitzEtaZ)legendZPlot->AddEntry(histoTrueDalitzEtaVsZ,("MC true #eta Dalitz"),"f");
	if (kCombZ)legendZPlot->AddEntry(histoTrueCombinatoricsVsZ,("MC true combinatorics"),"f");
	if (kHadZ)legendZPlot->AddEntry(histoTrueHadVsZ,("MC true hadronic background"),"f");
	legendZPlot->Draw();

	if(!thesis)DrawAliceLogoPerformance(floatLocationRightUp2[0], floatLocationRightUp2[1], floatLocationRightUp2[2], floatLocationRightUp2[3], 0.0001, textDate, collisionSystem, textGenerator, textPeriod,1200, 1000);   
	canvasZLin->Update();
	canvasZLin->SaveAs(Form("%s/Z_distribution_lin.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasZLin;          
	

	TCanvas * canvasZLog = new TCanvas("canvasZLog","",1200,1000);  // gives the page size
	canvasZLog->cd();
	canvasZLog->SetLeftMargin(0.13);
	//canvasZLog->SetTopMargin(0.05);            
	canvasZLog->SetLogy(1);

	DrawAutoGammaHistosWOLeg( histoConversionPointVsZData, 
							histoConversionPointVsZMonteCarlo, 
							"","Z (cm)",textYAxisZHisto,
							kTRUE, maxYFactorZPlot,minYValueZPlot,
							kFALSE,0. ,0.,
							kTRUE,  -201.,201.);		
 	if( etaCutNumber.CompareTo("3") == 0 ){
 		histoConversionPointVsZData->GetYaxis()->SetRangeUser(minYValueZPlot,4e-3);
 	} else if (etaCutNumber.CompareTo("2") == 0 || etaCutNumber.CompareTo("4") == 0){
		histoConversionPointVsZData->GetYaxis()->SetRangeUser(minYValueZPlot,10e-3);
	} else {
 		histoConversionPointVsZData->GetYaxis()->SetRangeUser(minYValueZPlot,1.1e-3);
 	}
	histoConversionPointVsZData->Draw("hist");
	histoTrueConversionPointVSZ->Draw("same,hist");
	histoTrueCombinatoricsVsZ->Draw("same,hist");
	histoTrueHadVsZ->Draw("same,hist");
	histoTrueDalitzPi0VsZ->Draw("same,hist");
	histoTrueDalitzEtaVsZ->Draw("same,hist");		
	histoConversionPointVsZData->Draw("same,hist,e");		
	histoConversionPointVsZMonteCarlo->Draw("same,hist,e");
	histoConversionPointVsZData->Draw("same,axis");		

	legendZPlot->Draw();
	if(!thesis)DrawAliceLogoPerformance(floatLocationRightUp2[0], floatLocationRightUp2[1], floatLocationRightUp2[2], floatLocationRightUp2[3], 0.0001, textDate, collisionSystem, textGenerator, textPeriod,1200,1000); 
	canvasZLog->Update();
	canvasZLog->SaveAs(Form("%s/Z_distribution_log.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasZLog;		
	
	TCanvas * canvasRLogSec = new TCanvas("canvasRLogSec","",1200,1000);  // gives the page size
	canvasRLogSec->cd();
	canvasRLogSec->SetLeftMargin(0.13);
	//canvasRLogSec->SetTopMargin(0.05);            
	canvasRLogSec->SetLogy(1);
	
	DrawAutoGammaHistosWOLeg( histoConversionPointVsRData, 
							histoConversionPointVsRMonteCarlo, 
							"","R (cm)",textYAxisRHisto,
							kTRUE, 0.5,minYValueRPlot,
							kFALSE,0. ,0.,
							kTRUE, 0.,180.);		
	if( etaCutNumber.CompareTo("2") == 0 || etaCutNumber.CompareTo("3") == 0 ){
		histoConversionPointVsRData->GetYaxis()->SetRangeUser(minYValueRPlot,2e-2);
	} else if (etaCutNumber.CompareTo("4") == 0){
		cout << "here" << endl;
		histoConversionPointVsRData->GetYaxis()->SetRangeUser(minYValueRPlot,4e-3);
	} else {
		histoConversionPointVsRData->GetYaxis()->SetRangeUser(minYValueRPlot,2.1e-3);
	}
	if ( optEnergy.CompareTo("pPb_5.023TeV") == 0 && multCutNumbers.CompareTo("00") != 0) histoConversionPointVsRData->GetYaxis()->SetRangeUser(minYValueRPlot,maxYValueRPlot);
	
	histoConversionPointVsRData->Draw("hist");
	histoTrueConversionPointPrimVSR->SetLineColor(kYellow-7);
	histoTrueConversionPointPrimVSR->SetFillColor(kYellow-7);
	histoTrueConversionPointPrimVSR->Draw("same,hist");
	histoTrueConversionPointSecVSR->SetLineColor(kSpring+2);
	histoTrueConversionPointSecVSR->SetFillColor(kSpring+2);
	histoTrueConversionPointSecVSR->Draw("same,hist");
	histoTrueCombinatoricsVsR->Draw("same,hist");
	histoTrueHadVsR->Draw("same,hist");
	histoTrueDalitzPi0VsR->Draw("same,hist");
	histoTrueDalitzEtaVsR->Draw("same,hist");		
	histoConversionPointVsRData->Draw("same,hist,e");		
	histoConversionPointVsRMonteCarlo->Draw("same,hist,e");
	histoConversionPointVsRData->Draw("same,axis");						

	if(etaCutNumber.CompareTo("2") != 0  && etaCutNumber.CompareTo("3") != 0) DrawIndividualTextSlicesR (vectorLogPlot , 0.025,"log");
	
	Double_t yLegendFullRSec= 0.027*(nLegendFullR);
	
	TLegend* legendRPlotSec;
	legendRPlotSec = new TLegend( 0.6,0.95-yLegendFullRSec,0.8,0.95);
	legendRPlotSec->SetTextSize(0.025);                        
	legendRPlotSec->SetFillColor(0);
	legendRPlotSec->SetBorderSize(0);
	legendRPlotSec->SetMargin(0.2);
	legendRPlotSec->AddEntry(histoConversionPointVsRData,("Data"),"l,p");
	legendRPlotSec->AddEntry(histoConversionPointVsRMonteCarlo,("MC conversion candidates"),"l");
	legendRPlotSec->AddEntry(histoTrueConversionPointPrimVSR,("MC true primary conversion"),"f");
	if (kSecFullR)legendRPlotSec->AddEntry(histoTrueConversionPointSecVSR,("MC true secondary conversion"),"f");
	if (kDalitzPi0FullR)legendRPlotSec->AddEntry(histoTrueDalitzPi0VsR,("MC true #pi^{0} Dalitz"),"f");
	if (kDalitzEtaFullR)legendRPlotSec->AddEntry(histoTrueDalitzEtaVsR,("MC true #eta Dalitz"),"f");
	if (kCombFullR)legendRPlotSec->AddEntry(histoTrueCombinatoricsVsR,("MC true combinatorics"),"f");
	if (kHadFullR)legendRPlotSec->AddEntry(histoTrueHadVsR,("MC true hadronic background"),"f");
	legendRPlotSec->Draw();

	if(!thesis)DrawAliceLogoPerformance(floatLocationRightUpR[0], floatLocationRightUpR[1], floatLocationRightUpR[2], floatLocationRightUpR[3], 0.0001, textDate,collisionSystem, textGenerator, textPeriod,1200,1000); 
	
	canvasRLogSec->Update();
	canvasRLogSec->SaveAs(Form("%s/R_distributionSec_log.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasRLogSec;

	TCanvas * canvasZLogSec = new TCanvas("canvasZLogSec","",1200,1000);  // gives the page size
	canvasZLogSec->cd();
	canvasZLogSec->SetLeftMargin(0.13);
	//canvasZLogSec->SetTopMargin(0.05);            
	canvasZLogSec->SetLogy(1);
	
	DrawAutoGammaHistosWOLeg( histoConversionPointVsZData, 
							histoConversionPointVsZMonteCarlo, 
							"","Z (cm)",textYAxisZHisto,
							kTRUE, 2.,minYValueZPlot,
							kFALSE,0. ,0.,
							kTRUE,  -201.,201.);		
	if( etaCutNumber.CompareTo("3") == 0 ){
		histoConversionPointVsZData->GetYaxis()->SetRangeUser(minYValueZPlot,2e-3);
	} else if (etaCutNumber.CompareTo("2") == 0 ){
		histoConversionPointVsZData->GetYaxis()->SetRangeUser(minYValueZPlot,3e-2);
	} else if (etaCutNumber.CompareTo("4") == 0){
		histoConversionPointVsZData->GetYaxis()->SetRangeUser(minYValueZPlot,10e-3);
	} else {
		histoConversionPointVsZData->GetYaxis()->SetRangeUser(minYValueZPlot,1.1e-3);
	}
	histoConversionPointVsZData->Draw("hist");
	histoTrueConversionPointPrimVSZ->SetLineColor(kYellow-7);
	histoTrueConversionPointPrimVSZ->SetFillColor(kYellow-7);
	histoTrueConversionPointPrimVSZ->Draw("same,hist");
	histoTrueConversionPointSecVSZ->SetLineColor(kSpring+2);
	histoTrueConversionPointSecVSZ->SetFillColor(kSpring+2);
	histoTrueConversionPointSecVSZ->Draw("same,hist");
	histoTrueCombinatoricsVsZ->Draw("same,hist");
	histoTrueHadVsZ->Draw("same,hist");
	histoTrueDalitzPi0VsZ->Draw("same,hist");
	histoTrueDalitzEtaVsZ->Draw("same,hist");		
	histoConversionPointVsZData->Draw("same,hist,e");		
	histoConversionPointVsZMonteCarlo->Draw("same,hist,e");
	histoConversionPointVsZData->Draw("same,axis");						

	if (nLegendZ > 4){
		Double_t yLegendZSec= 0.027*(nLegendZ/2);
		TLegend* legendZPlotSec;
		legendZPlotSec = new TLegend( 0.145,0.94-yLegendZSec,0.8,0.94);
		legendZPlotSec->SetTextSize(0.025);                        
		legendZPlotSec->SetFillColor(0);
		legendZPlotSec->SetBorderSize(0);
		legendZPlotSec->SetNColumns(2);
		legendZPlotSec->SetMargin(0.14);
		legendZPlotSec->AddEntry(histoConversionPointVsZData,("Data"),"l,p");
		legendZPlotSec->AddEntry(histoConversionPointVsZMonteCarlo,("MC conversion candidates"),"l");
		legendZPlotSec->AddEntry(histoTrueConversionPointPrimVSZ,("MC true primary conversion"),"f");
		if (kSecZ)legendZPlotSec->AddEntry(histoTrueConversionPointSecVSZ,("MC true secondary conversion"),"f");
		if (kDalitzPi0Z)legendZPlotSec->AddEntry(histoTrueDalitzPi0VsZ,("MC true #pi^{0} Dalitz"),"f");
		if (kDalitzEtaZ)legendZPlotSec->AddEntry(histoTrueDalitzEtaVsZ,("MC true #eta Dalitz"),"f");
		if (kCombZ)legendZPlotSec->AddEntry(histoTrueCombinatoricsVsZ,("MC true combinatorics"),"f");
		if (kHadZ)legendZPlotSec->AddEntry(histoTrueHadVsZ,("MC true hadronic background"),"f");
		legendZPlotSec->Draw();
	} else {
		TLegend* legendZPlotSec;
		Double_t yLegendZSec= 0.027*(nLegendZ);
		legendZPlotSec = new TLegend( 0.145,0.94-yLegendZSec,0.345,0.94);
		legendZPlotSec->SetTextSize(0.025);                        
		legendZPlotSec->SetFillColor(0);
		legendZPlotSec->SetBorderSize(0);
		legendZPlotSec->SetMargin(0.2);
		legendZPlotSec->AddEntry(histoConversionPointVsZData,("Data"),"l,p");
		legendZPlotSec->AddEntry(histoConversionPointVsZMonteCarlo,("MC conversion candidates"),"l");
		legendZPlotSec->AddEntry(histoTrueConversionPointPrimVSZ,("MC true primary conversion"),"f");
		if (kSecZ)legendZPlotSec->AddEntry(histoTrueConversionPointSecVSZ,("MC true secondary conversion"),"f");
		if (kDalitzPi0Z)legendZPlotSec->AddEntry(histoTrueDalitzPi0VsZ,("MC true #pi^{0} Dalitz"),"f");
		if (kDalitzEtaZ)legendZPlotSec->AddEntry(histoTrueDalitzEtaVsZ,("MC true #eta Dalitz"),"f");
		if (kCombZ)legendZPlotSec->AddEntry(histoTrueCombinatoricsVsZ,("MC true combinatorics"),"f");
		if (kHadZ)legendZPlotSec->AddEntry(histoTrueHadVsZ,("MC true hadronic background"),"f");
		legendZPlotSec->Draw();
	}
	
	
	if(!thesis)DrawAliceLogoPerformance(floatLocationRightUp2[0], floatLocationRightUp2[1], floatLocationRightUp2[2], floatLocationRightUp2[3], 0.0001, textDate,collisionSystem, textGenerator, textPeriod,1200,1000); 
	
	canvasZLogSec->Update();
	canvasZLogSec->SaveAs(Form("%s/Z_distributionSec_log.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasZLogSec;
	
	// --------------- Giving Phi in R for several bins -------------------------------------------------------------
	for(Int_t iR = iRStart; iR < nBinsR; iR++){
		if (histoMappingPhiInRData[iR]->GetEntries()>0 && histoMappingPhiInRMonteCarlo[iR]->GetEntries()>0){
			TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
			if (iR == 11)	canvasSingleBin->SetLogy(1);
			canvasSingleBin->SetTopMargin(0.04);
			canvasSingleBin->cd();
			TString nameHistoPhiInRData=Form("#varphi in R:  %s",arrayNamesRangesRBins[iR].Data());
			TLatex *latexBinning = new TLatex(doubleLatexNamingBinsX, doubleLatexNamingBinsY,nameHistoPhiInRData); 
			SetStyleTLatex( latexBinning, sizeTextNameBins,2);
			TLatex *latexBinning2 = new TLatex(doubleLatexNamingBinsX2, doubleLatexNamingBinsY2,arrayNamesRBins[iR]); 
			SetStyleTLatex( latexBinning2, sizeTextNameBins,2);
			
			DrawAutoGamma3Histos( histoMappingPhiInRData[iR], 
											histoMappingPhiInRMonteCarlo[iR],
											histoTrueMappingPhiInRMonteCarlo[iR],
											"","#varphi",textYAxisSlicesHisto,
											kTRUE,1.2 ,1e-10,
											kFALSE,0. ,0.,
											kFALSE, 0.,180.);
			latexBinning->Draw("same");  
			latexBinning2->Draw("same");  
			if (outer) {
				TLatex *latexBinning3 = new TLatex(doubleLatexNamingBinsX2, doubleLatexNamingBinsY3,arrayNamesAddRBins[iR]);
				SetStyleTLatex( latexBinning3, sizeTextNameBins,2);
				latexBinning3->Draw("same"); 
			}
			
// 			if(!thesis)DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
			
			canvasSingleBin->Update();
			canvasSingleBin->SaveAs(Form("%s/Phi_in_R_%i.%s",outputDirectory.Data(),iR,suffix.Data()));
			delete canvasSingleBin;
		}
	}               
	
	// -------------- Giving Z in R for all bins ----------------------------------------------
	
	for(Int_t iR = iRStart; iR < nBinsR; iR++){
		if (histoMappingZInRData[iR]->GetEntries()>0 && histoMappingZInRMonteCarlo[iR]->GetEntries()>0){
			TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
			canvasSingleBin->SetLogy(1);
			canvasSingleBin->cd();
			Float_t rangeZ = arrayRangeZInR[iR];
			TString nameHistoPhiInRData=Form("Z in R:  %s",arrayNamesRangesRBins[iR].Data());
			TLatex *latexBinning = new TLatex(doubleLatexNamingBinsX, doubleLatexNamingBinsY,nameHistoPhiInRData); 
			SetStyleTLatex( latexBinning, sizeTextNameBins,2);
			TLatex *latexBinning2 = new TLatex(doubleLatexNamingBinsX2, doubleLatexNamingBinsY2,arrayNamesRBins[iR]); 
			SetStyleTLatex( latexBinning2, sizeTextNameBins,2);
			if (iR == nBinsR -2 || iR == nBinsR -1){
				histoMappingZInRData[iR]->GetYaxis()->SetRangeUser(histoMappingZInRData[iR]->GetMaximum()*1e-2, 3*histoMappingZInRData[iR]->GetMaximum());
			} else {
				histoMappingZInRData[iR]->GetYaxis()->SetRangeUser(histoMappingZInRData[iR]->GetMaximum()*1e-1, 3*histoMappingZInRData[iR]->GetMaximum());
			}
			DrawAutoGamma3Histos( histoMappingZInRData[iR], 
											histoMappingZInRMonteCarlo[iR], 
											histoTrueMappingZInRMonteCarlo[iR],
											"","Z",textYAxisSlicesHisto,
											kFALSE,2 ,1e-9,
											kFALSE,0. ,0.,
											kTRUE, -rangeZ,rangeZ);	
			latexBinning->Draw("same");  
			latexBinning2->Draw("same");  
			if (outer) {
				TLatex *latexBinning3 = new TLatex(doubleLatexNamingBinsX2, doubleLatexNamingBinsY3,arrayNamesAddRBins[iR]);
				SetStyleTLatex( latexBinning3, sizeTextNameBins,2);
				latexBinning3->Draw("same"); 
			}
// 			if(!thesis)DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
			
			canvasSingleBin->Update();
			canvasSingleBin->SaveAs(Form("%s/Z_in_R_%i.%s",outputDirectory.Data(),iR,suffix.Data()));
			delete canvasSingleBin;
		}
	}               
	
		
	// ************************************* Phi in Z for all bins *****************************************
	for(Int_t iZ = minimumBinZNormal ; iZ < (nBinsZ-minimumBinZNormal); iZ++){
		if (histoMappingPhiInZData[iZ]->GetEntries() > 0 && histoMappingPhiInZMonteCarlo[iZ]->GetEntries() > 0 ){
			TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
	// 		canvasSingleBin->SetLogy(1);
			canvasSingleBin->SetTopMargin(0.04);
			canvasSingleBin->cd();
			TString nameHistoPhiInZData=Form("#varphi in Z:  %s",arrayNamesZBins[iZ].Data());
			TLatex *latexBinning = new TLatex(doubleLatexNamingBinsX, doubleLatexNamingBinsY,nameHistoPhiInZData); // Bo: this was modified
			latexBinning->SetNDC();
			latexBinning->SetTextColor(1);
			latexBinning->SetTextFont(62);
			latexBinning->SetTextSize(sizeTextNameBins);
			latexBinning->SetLineWidth(2);       
			
			DrawAutoGamma3Histos(histoMappingPhiInZData[iZ], 
											histoMappingPhiInZMonteCarlo[iZ],
											histoTrueMappingPhiInZMonteCarlo[iZ],
											"","#varphi",textYAxisSlicesHisto,
											kTRUE,1.2 ,1e-9,
											kFALSE,0. ,0.,
											kFALSE, -0,0);
			latexBinning->Draw("same");  
	// 		if(!thesis)DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
			
			canvasSingleBin->Update();
			canvasSingleBin->SaveAs(Form("%s/Phi_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
			delete canvasSingleBin;
		}
	}               
	
	// ------------- Giving SPD Z in phi in singleplot
	for(Int_t iZ = minimumBinZSPD ; iZ < (nBinsZ-minimumBinZSPD); iZ++){
		if (histoMappingSPDPhiInZData[iZ]->GetEntries() > 0 && histoMappingSPDPhiInZMonteCarlo[iZ]->GetEntries() > 0 ){
			TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
	// 		canvasSingleBin->SetLogy(1);
			canvasSingleBin->SetTopMargin(0.04);
			canvasSingleBin->cd();
			nameHistoSPDPhiInZData=Form("0 cm < R < 9.5 cm: #varphi in Z:  %s",arrayNamesZBins[iZ].Data());
			TLatex *latexBinning = new TLatex(doubleLatexNamingBinsX, doubleLatexNamingBinsY,nameHistoSPDPhiInZData); // Bo: this was modified
			latexBinning->SetNDC();
			latexBinning->SetTextColor(1);
			latexBinning->SetTextFont(62);
			latexBinning->SetTextSize(sizeTextNameBins);
			latexBinning->SetLineWidth(2);       
			
			DrawAutoGammaHistosMaterial(histoMappingSPDPhiInZData[iZ], 
											histoMappingSPDPhiInZMonteCarlo[iZ],
											"","#varphi",textYAxisSlicesHisto,
											kTRUE,1.2 ,5e-10,
											kFALSE,0. ,0.,
											kFALSE, -0,0);
			latexBinning->Draw("same");  
	// 		if(!thesis)DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
			
			canvasSingleBin->Update();
			canvasSingleBin->SaveAs(Form("%s/SPD_Phi_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
			delete canvasSingleBin;
		}
	}               
	
	for(Int_t iZ = minimumBinZSPDTh ; iZ < (nBinsZ-minimumBinZSPDTh); iZ++){
		if (histoMappingSPDThPhiInZData[iZ]->GetEntries() > 0 && histoMappingSPDThPhiInZMonteCarlo[iZ]->GetEntries() > 0 ){
			TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
	// 		canvasSingleBin->SetLogy(1);
			canvasSingleBin->SetTopMargin(0.04);
			canvasSingleBin->cd();
			nameHistoSPDThPhiInZData=Form("9.5 cm < R < 13 cm: #varphi in Z:  %s",arrayNamesZBins[iZ].Data());
			TLatex *latexBinning = new TLatex(doubleLatexNamingBinsX, doubleLatexNamingBinsY,nameHistoSPDThPhiInZData); // Bo: this was modified
			latexBinning->SetNDC();
			latexBinning->SetTextColor(1);
			latexBinning->SetTextFont(62);
			latexBinning->SetTextSize(sizeTextNameBins);
			latexBinning->SetLineWidth(2);       
			
			DrawAutoGammaHistosMaterial(histoMappingSPDThPhiInZData[iZ], 
											histoMappingSPDThPhiInZMonteCarlo[iZ],
											"","#varphi",textYAxisSlicesHisto,
											kTRUE,1.2 ,5e-10,
											kFALSE,0. ,0.,
											kFALSE, -0,0);
			latexBinning->Draw("same");  
	// 		if(!thesis)DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
			
			canvasSingleBin->Update();
			canvasSingleBin->SaveAs(Form("%s/SPDTh_Phi_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
			delete canvasSingleBin;
		}
	}               

	for(Int_t iZ = minimumBinZSDD ; iZ < (nBinsZ-minimumBinZSDD); iZ++){
		if (histoMappingSDDPhiInZData[iZ]->GetEntries() > 0 && histoMappingSDDPhiInZMonteCarlo[iZ]->GetEntries() > 0 ){
			TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
	// 		canvasSingleBin->SetLogy(1);
			canvasSingleBin->SetTopMargin(0.04);
			canvasSingleBin->cd();
			nameHistoSDDPhiInZData=Form("13 cm < R < 35 cm: #varphi in Z:  %s",arrayNamesZBins[iZ].Data());
			TLatex *latexBinning = new TLatex(doubleLatexNamingBinsX, doubleLatexNamingBinsY,nameHistoSDDPhiInZData); // Bo: this was modified
			latexBinning->SetNDC();
			latexBinning->SetTextColor(1);
			latexBinning->SetTextFont(62);
			latexBinning->SetTextSize(sizeTextNameBins);
			latexBinning->SetLineWidth(2);       
			
			DrawAutoGammaHistosMaterial(histoMappingSDDPhiInZData[iZ], 
											histoMappingSDDPhiInZMonteCarlo[iZ],
											"","#varphi",textYAxisSlicesHisto,
											kTRUE,1.2 ,5e-10,
											kFALSE,0. ,0.,
											kFALSE, -0,0);
			latexBinning->Draw("same");  
	// 		if(!thesis)DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
			
			canvasSingleBin->Update();
			canvasSingleBin->SaveAs(Form("%s/SDD_Phi_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
			delete canvasSingleBin;
		}
	}               
	
	for(Int_t iZ = minimumBinZSSD ; iZ < (nBinsZ-minimumBinZSSD); iZ++){
		if (histoMappingSSDPhiInZData[iZ]->GetEntries() > 0 && histoMappingSSDPhiInZMonteCarlo[iZ]->GetEntries() > 0 ){
			TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
	// 		canvasSingleBin->SetLogy(1);
			canvasSingleBin->SetTopMargin(0.04);
			canvasSingleBin->cd();
			nameHistoSSDPhiInZData=Form("35 cm < R < 55 cm: #varphi in Z:  %s",arrayNamesZBins[iZ].Data());
			TLatex *latexBinning = new TLatex(doubleLatexNamingBinsX, doubleLatexNamingBinsY,nameHistoSSDPhiInZData); // Bo: this was modified
			latexBinning->SetNDC();
			latexBinning->SetTextColor(1);
			latexBinning->SetTextFont(62);
			latexBinning->SetTextSize(sizeTextNameBins);
			latexBinning->SetLineWidth(2);       
			
			DrawAutoGammaHistosMaterial(histoMappingSSDPhiInZData[iZ], 
											histoMappingSSDPhiInZMonteCarlo[iZ],
											"","#varphi",textYAxisSlicesHisto,
											kTRUE,1.2 ,5e-10,
											kFALSE,0. ,0.,
											kFALSE, -0,0);
			latexBinning->Draw("same");  
	// 		if(!thesis)DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
			
			canvasSingleBin->Update();
			canvasSingleBin->SaveAs(Form("%s/SSD_Phi_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
			delete canvasSingleBin;
		}
	}               
	
	
	// ************************************* ITSTPC Phi in Z for all bins *****************************************
	for(Int_t iZ = minimumBinZITSTPC ; iZ < (nBinsZ-minimumBinZITSTPC); iZ++){
		if (histoMappingITSTPCPhiInZData[iZ]->GetEntries() > 0 && histoMappingITSTPCPhiInZMonteCarlo[iZ]->GetEntries() > 0 ){
			TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
	// 		canvasSingleBin->SetLogy(1);
			canvasSingleBin->SetTopMargin(0.04);
			canvasSingleBin->cd();
			TString nameHistoITSTPCPhiInZData=Form("55 cm < R < 72.5 cm: #varphi in Z:  %s",arrayNamesZBins[iZ].Data());
			TLatex *latexBinning = new TLatex(doubleLatexNamingBinsX, doubleLatexNamingBinsY,nameHistoITSTPCPhiInZData); // Bo: this was modified
			latexBinning->SetNDC();
			latexBinning->SetTextColor(1);
			latexBinning->SetTextFont(62);
			latexBinning->SetTextSize(sizeTextNameBins);
			latexBinning->SetLineWidth(2);       
			
			DrawAutoGammaHistosMaterial(histoMappingITSTPCPhiInZData[iZ], 
							histoMappingITSTPCPhiInZMonteCarlo[iZ],
							"","#varphi",textYAxisSlicesHisto,
							kTRUE,1.2 ,5e-10,
							kFALSE,0. ,0.,
							kFALSE, -0,0);
			latexBinning->Draw("same");  
	// 		if(!thesis)DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
			
			canvasSingleBin->Update();
			canvasSingleBin->SaveAs(Form("%s/ITSTPC_Phi_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
			delete canvasSingleBin;
		}
	}               
	

	// ----- Giving R in Z in SinglePlot    
	for(Int_t iZ = minimumBinZNormal ; iZ < (nBinsZ-minimumBinZNormal); iZ++){
		if (histoMappingRInZData[iZ]->GetEntries() > 0 && histoMappingRInZMonteCarlo[iZ]->GetEntries() > 0 ){
			TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
			canvasSingleBin->SetLogy(1);
			canvasSingleBin->cd();
			TString nameHistoRInZData=Form("R in Z:  %s",arrayNamesZBins[iZ].Data());
			TLatex *latexBinning = new TLatex(doubleLatexNamingBinsX, doubleLatexNamingBinsY,nameHistoRInZData); // Bo: this was modified
			latexBinning->SetNDC();
			latexBinning->SetTextColor(1);
			latexBinning->SetTextFont(62);
			latexBinning->SetTextSize(sizeTextNameBins);
			latexBinning->SetLineWidth(2);       
			
			DrawAutoGamma3Histos(histoMappingRInZData[iZ], 
											histoMappingRInZMonteCarlo[iZ],
											histoTrueMappingRInZMonteCarlo[iZ],
											"","R",textYAxisSlicesHisto,
											kTRUE,8 ,5e-10,
											kFALSE,0. ,0.,
											kTRUE, 0,180);
			latexBinning->Draw("same");  
	// 		if(!thesis)DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
			
			canvasSingleBin->Update();
			canvasSingleBin->SaveAs(Form("%s/R_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
			delete canvasSingleBin;
		}
	}               
	
	
		
	for(Int_t iR = iRStart; iR < nBinsR; iR++){
		if (histoMappingRatioPhiInR[iR]->GetEntries()>0){
			TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
			DrawGammaCanvasSettings( canvasSingleBin, 0.07, 0.015, 0.035, 0.075);
			if (boolRatioLog) canvasSingleBin->SetLogy(1);
			canvasSingleBin->cd();
			TString nameHistoPhiInRData=Form("#varphi in R:  %s",arrayNamesRBins[iR].Data());
			TLatex *latexBinning = new TLatex(doubleLatexNamingBinsRatioX, doubleLatexNamingBinsY,nameHistoPhiInRData); // Bo: this was modified
			latexBinning->SetNDC();
			latexBinning->SetTextColor(1);
			latexBinning->SetTextFont(62);
			latexBinning->SetTextSize(sizeTextNameBins);
			latexBinning->SetLineWidth(2);       
			
			DrawRatioGammaHisto(    histoMappingRatioPhiInR[iR], 
											"","#varphi","norm Data/norm MC",
											kFALSE,3 ,0.000001,
											kTRUE, minYValueRatio , maxYValueRatio,
											kFALSE, 0.,180.);
			DrawGammaLines(0,2*TMath::Pi(),1,1,0.02);              
			
			latexBinning->Draw("same");  
// 			if(!thesis)DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
			
			canvasSingleBin->Update();
			canvasSingleBin->SaveAs(Form("%s/Ratio_Phi_in_R_%i.%s",outputDirectory.Data(),iR,suffix.Data()));
			delete canvasSingleBin;
		}
	}               
		
	for(Int_t iR = iRStart; iR < nBinsR; iR++){
		if (histoMappingRatioZInR[iR]->GetEntries()>0){	
			TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
			DrawGammaCanvasSettings( canvasSingleBin, 0.07, 0.015, 0.035, 0.075);
			if (boolRatioLog) canvasSingleBin->SetLogy(1);
			canvasSingleBin->cd();
			Float_t rangeZ = arrayRangeZInR[iR];
			TString nameHistoPhiInRData=Form("Z in R:  %s",arrayNamesRBins[iR].Data());
			TLatex *latexBinning = new TLatex(doubleLatexNamingBinsRatioX, doubleLatexNamingBinsY,nameHistoPhiInRData); // Bo: this was modified
			latexBinning->SetNDC();
			latexBinning->SetTextColor(1);
			latexBinning->SetTextFont(62);
			latexBinning->SetTextSize(sizeTextNameBins);
			latexBinning->SetLineWidth(2);       
			
			DrawRatioGammaHisto( histoMappingRatioZInR[iR], 
											"","Z","norm Data/norm MC",
											kFALSE,3 ,0.000001,
											kTRUE, minYValueRatio , maxYValueRatio,
											kTRUE, -rangeZ,rangeZ);
			DrawGammaLines(-rangeZ,rangeZ,1,1,0.02);                
			latexBinning->Draw("same");  
// 			if(!thesis)DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
			
			canvasSingleBin->Update();
			canvasSingleBin->SaveAs(Form("%s/Ratio_Z_in_R_%i.%s",outputDirectory.Data(),iR,suffix.Data()));
			delete canvasSingleBin;
		}
	}
	
	for(Int_t iZ = minimumBinZNormal ; iZ < (nBinsZ-minimumBinZNormal); iZ++){
		if (histoMappingRatioPhiInZ[iZ]->GetEntries()>0){	
			TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
			DrawGammaCanvasSettings( canvasSingleBin, 0.07, 0.015, 0.035, 0.075);
			if (boolRatioLog) canvasSingleBin->SetLogy(1);
			canvasSingleBin->cd();
			TString nameHistoPhiInRData=Form("#varphi in Z:  %s",arrayNamesZBins[iZ].Data());
			TLatex *latexBinning = new TLatex(doubleLatexNamingBinsRatioX, doubleLatexNamingBinsY,nameHistoPhiInRData); // Bo: this was modified
			latexBinning->SetNDC();
			latexBinning->SetTextColor(1);
			latexBinning->SetTextFont(62);
			latexBinning->SetTextSize(sizeTextNameBins);
			latexBinning->SetLineWidth(2);       
			
			DrawRatioGammaHisto( histoMappingRatioPhiInZ[iZ], 
											"","#varphi","norm Data/norm MC",
											kFALSE,3 ,0.000001,
											kTRUE, minYValueRatio , maxYValueRatio,
											kFALSE, 0,0);
			DrawGammaLines(0,2*TMath::Pi(),1,1,0.02);              
			latexBinning->Draw("same");  
	// 		if(!thesis)DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
			canvasSingleBin->Update();
			canvasSingleBin->SaveAs(Form("%s/Ratio_Phi_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
			delete canvasSingleBin;
		}
	}               
	
	for(Int_t iZ = minimumBinZSPD ; iZ < (nBinsZ-minimumBinZSPD); iZ++){
		if (histoMappingRatioSPDPhiInZ[iZ]->GetEntries()>0){	
			TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
			DrawGammaCanvasSettings( canvasSingleBin, 0.07, 0.015, 0.035, 0.075);
			if (boolRatioLog) canvasSingleBin->SetLogy(1);
			canvasSingleBin->cd();
			TString nameHistoPhiInRData=Form("0 cm < R < 9.5 cm: #varphi in Z:  %s",arrayNamesZBins[iZ].Data());
			TLatex *latexBinning = new TLatex(doubleLatexNamingBinsRatioX, doubleLatexNamingBinsY,nameHistoPhiInRData); // Bo: this was modified
			latexBinning->SetNDC();
			latexBinning->SetTextColor(1);
			latexBinning->SetTextFont(62);
			latexBinning->SetTextSize(sizeTextNameBins);
			latexBinning->SetLineWidth(2);       
			DrawRatioGammaHisto( histoMappingRatioSPDPhiInZ[iZ], 
											"","#varphi","norm Data/norm MC",
											kFALSE,3 ,0.000001,
											kTRUE, minYValueRestRRatio , maxYValueRestRRatio,
											kFALSE, 0,0);
			DrawGammaLines(0,2*TMath::Pi(),1,1,0.02);              
			latexBinning->Draw("same");  
	// 		if(!thesis)DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
			canvasSingleBin->Update();
			canvasSingleBin->SaveAs(Form("%s/Ratio_SPD_Phi_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
			delete canvasSingleBin;
		}
	}               
		
	for(Int_t iZ = minimumBinZSPDTh ; iZ < (nBinsZ-minimumBinZSPDTh); iZ++){
		if (histoMappingRatioSPDThPhiInZ[iZ]->GetEntries()>0){	
			TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
			DrawGammaCanvasSettings( canvasSingleBin, 0.07, 0.015, 0.035, 0.075);
			if (boolRatioLog) canvasSingleBin->SetLogy(1);
			canvasSingleBin->cd();
			TString nameHistoPhiInRData=Form("9.5 cm < R < 13 cm: #varphi in Z:  %s",arrayNamesZBins[iZ].Data());
			TLatex *latexBinning = new TLatex(doubleLatexNamingBinsRatioX, doubleLatexNamingBinsY,nameHistoPhiInRData); // Bo: this was modified
			latexBinning->SetNDC();
			latexBinning->SetTextColor(1);
			latexBinning->SetTextFont(62);
			latexBinning->SetTextSize(sizeTextNameBins);
			latexBinning->SetLineWidth(2);       
			DrawRatioGammaHisto( histoMappingRatioSPDThPhiInZ[iZ], 
											"","#varphi","norm Data/norm MC",
											kFALSE,3 ,0.000001,
											kTRUE, minYValueRestRRatio , maxYValueRestRRatio,
											kFALSE, 0,0);
			DrawGammaLines(0,2*TMath::Pi(),1,1,0.02);              
			latexBinning->Draw("same");  
	// 		if(!thesis)DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
			canvasSingleBin->Update();
			canvasSingleBin->SaveAs(Form("%s/Ratio_SPDTh_Phi_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
			delete canvasSingleBin;
		}
	}               

	for(Int_t iZ = minimumBinZSDD ; iZ < (nBinsZ-minimumBinZSDD); iZ++){
		if (histoMappingRatioSDDPhiInZ[iZ]->GetEntries()>0){	
			TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
			DrawGammaCanvasSettings( canvasSingleBin, 0.07, 0.015, 0.035, 0.075);
			if (boolRatioLog) canvasSingleBin->SetLogy(1);
			canvasSingleBin->cd();
			TString nameHistoPhiInRData=Form("13 cm < R < 35 cm: #varphi in Z:  %s",arrayNamesZBins[iZ].Data());
			TLatex *latexBinning = new TLatex(doubleLatexNamingBinsRatioX, doubleLatexNamingBinsY,nameHistoPhiInRData); // Bo: this was modified
			latexBinning->SetNDC();
			latexBinning->SetTextColor(1);
			latexBinning->SetTextFont(62);
			latexBinning->SetTextSize(sizeTextNameBins);
			latexBinning->SetLineWidth(2);       
			DrawRatioGammaHisto( histoMappingRatioSDDPhiInZ[iZ], 
											"","#varphi","norm Data/norm MC",
											kFALSE,3 ,0.000001,
											kTRUE, minYValueRestRRatio , maxYValueRestRRatio,
											kFALSE, 0,0);
			DrawGammaLines(0,2*TMath::Pi(),1,1,0.02);              
			latexBinning->Draw("same");  
	// 		if(!thesis)DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
			canvasSingleBin->Update();
			canvasSingleBin->SaveAs(Form("%s/Ratio_SDD_Phi_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
			delete canvasSingleBin;
		}
	}               

	for(Int_t iZ = minimumBinZSSD ; iZ < (nBinsZ-minimumBinZSSD); iZ++){
		if (histoMappingRatioSSDPhiInZ[iZ]->GetEntries()>0){	
			TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
			DrawGammaCanvasSettings( canvasSingleBin, 0.07, 0.015, 0.035, 0.075);
			if (boolRatioLog) canvasSingleBin->SetLogy(1);
			canvasSingleBin->cd();
			TString nameHistoPhiInRData=Form("35 cm < R < 55 cm: #varphi in Z:  %s",arrayNamesZBins[iZ].Data());
			TLatex *latexBinning = new TLatex(doubleLatexNamingBinsRatioX, doubleLatexNamingBinsY,nameHistoPhiInRData); // Bo: this was modified
			latexBinning->SetNDC();
			latexBinning->SetTextColor(1);
			latexBinning->SetTextFont(62);
			latexBinning->SetTextSize(sizeTextNameBins);
			latexBinning->SetLineWidth(2);       
			DrawRatioGammaHisto( histoMappingRatioSSDPhiInZ[iZ], 
											"","#varphi","norm Data/norm MC",
											kFALSE,3 ,0.000001,
											kTRUE, minYValueRestRRatio , maxYValueRestRRatio,
											kFALSE, 0,0);
			DrawGammaLines(0,2*TMath::Pi(),1,1,0.02);              
			latexBinning->Draw("same");  
	// 		if(!thesis)DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
			canvasSingleBin->Update();
			canvasSingleBin->SaveAs(Form("%s/Ratio_SSD_Phi_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
			delete canvasSingleBin;
		}
	}               

	
	for(Int_t iZ = minimumBinZITSTPC ; iZ < (nBinsZ-minimumBinZITSTPC); iZ++){
		if (histoMappingRatioITSTPCPhiInZ[iZ]->GetEntries()>0){	
			TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
			DrawGammaCanvasSettings( canvasSingleBin, 0.07, 0.015, 0.035, 0.075);
			if (boolRatioLog) canvasSingleBin->SetLogy(1);
			canvasSingleBin->cd();
			TString nameHistoPhiInRData=Form("55 cm < R < 72.5 cm: #varphi in Z:  %s",arrayNamesZBins[iZ].Data());
			TLatex *latexBinning = new TLatex(doubleLatexNamingBinsRatioX, doubleLatexNamingBinsY,nameHistoPhiInRData); // Bo: this was modified
			latexBinning->SetNDC();
			latexBinning->SetTextColor(1);
			latexBinning->SetTextFont(62);
			latexBinning->SetTextSize(sizeTextNameBins);
			latexBinning->SetLineWidth(2);       
			DrawRatioGammaHisto( histoMappingRatioITSTPCPhiInZ[iZ], 
							"","#varphi","norm Data/norm MC",
							kFALSE,3 ,0.000001,
							kTRUE, minYValueRestRRatio , maxYValueRestRRatio,
							kFALSE, 0,0);
			DrawGammaLines(0,2*TMath::Pi(),1,1,0.02);              
			latexBinning->Draw("same");  
	// 		if(!thesis)DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
			canvasSingleBin->Update();
			canvasSingleBin->SaveAs(Form("%s/Ratio_ITSTPC_Phi_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
			delete canvasSingleBin;
		}
	}               
			
	for(Int_t iZ = minimumBinZNormal ; iZ < (nBinsZ-minimumBinZNormal); iZ++){
		if (histoMappingRatioRInZ[iZ]->GetEntries()>0){	
			TCanvas * canvasSingleBin = new TCanvas("canvasSingleBin","",500,350);  // gives the page size
			DrawGammaCanvasSettings( canvasSingleBin, 0.07, 0.015, 0.035, 0.075);
			if (boolRatioLog) canvasSingleBin->SetLogy(1);
			canvasSingleBin->cd();
			TString nameHistoPhiInRData=Form("R in Z:  %s",arrayNamesZBins[iZ].Data());
			TLatex *latexBinning = new TLatex(doubleLatexNamingBinsRatioX, doubleLatexNamingBinsY,nameHistoPhiInRData); // Bo: this was modified
			latexBinning->SetNDC();
			latexBinning->SetTextColor(1);
			latexBinning->SetTextFont(62);
			latexBinning->SetTextSize(sizeTextNameBins);
			latexBinning->SetLineWidth(2);       
			
			DrawRatioGammaHisto( histoMappingRatioRInZ[iZ], 
											"","R (cm)","norm Data/norm MC",
											kFALSE,3 ,0.000001,
											kTRUE, minYValueRatio , maxYValueRatio,
											kFALSE, 0,0);
			DrawGammaLines(0,180,1,1,0.02);         
			latexBinning->Draw("same");  
	// 		if(!thesis)DrawAliceText(floatLocationRightUpText[0],floatLocationRightUpText[1], floatLocationRightUpText[3]);             
			canvasSingleBin->Update();
			canvasSingleBin->SaveAs(Form("%s/Ratio_R_in_Z_%i.%s",outputDirectory.Data(),iZ,suffix.Data()));
			delete canvasSingleBin;
		}
	}               
		
	if (histoConversionPointVsZCentralElectrodeData->GetEntries() > 0 && histoConversionPointVsZCentralElectrodeMonteCarlo->GetEntries() > 0){
		TCanvas * canvasCentralElectrode = new TCanvas("canvasCentralElectrode","",1200,1000);  // gives the page size
		canvasCentralElectrode->SetLogy(1);
		canvasCentralElectrode->cd();
		canvasCentralElectrode->SetLeftMargin(0.15);
         DrawAutoGammaHistos( histoConversionPointVsZCentralElectrodeData, 
                              histoConversionPointVsZCentralElectrodeMonteCarlo, 
                              "","R (cm) ",textYAxisCEHisto,
                              kTRUE, 2.5,2.e-8,
                              kFALSE,0. ,0.,
                              kTRUE, 68.,140.);
		
         TLatex *latexBinning = new TLatex(doubleLatexNamingBinsX+0.05, doubleLatexNamingBinsY,"Central Electrode"); 
         SetStyleTLatex( latexBinning, sizeTextNameBins,2);
         TLatex *latexBinning2 = new TLatex(doubleLatexNamingBinsX+0.05, doubleLatexNamingBinsY2,"-1 cm < Z < 1 cm"); 
         SetStyleTLatex( latexBinning2, sizeTextNameBins,2);
         latexBinning->Draw("same");	
         latexBinning2->Draw("same");	
         if(!thesis)DrawAliceLogo1D(floatLocationRightUp [0],floatLocationRightUp[1],floatLocationRightUp[2],floatLocationRightUp[3],1200,1000);                              
		canvasCentralElectrode->Update();
		canvasCentralElectrode->SaveAs(Form("%s/RCentralElectrode.%s",outputDirectory.Data(),suffix.Data()));
		delete canvasCentralElectrode;
	}
	
	TCanvas * canvasZVertexPos = new TCanvas("canvasZVertexPos","",1200,1000);  // gives the page size
	DrawGammaCanvasSettings( canvasZVertexPos, 0.12, 0.015, 0.035, 0.09);
// 		canvasZVertexPos->SetLogy(1);
	canvasZVertexPos->cd();
				DrawAutoGammaHistos( histoVertexZData, 
												histoVertexZMonteCarlo, 
												"","Z_{vtx} (cm) ","#frac{1}{N_{ev}} #frac{dN_{ev}}{dZ_{vtx}} (cm^{-1})",
												kTRUE, 1.1,2.e-8,
												kFALSE,0. ,0.,
												kTRUE, -20.,20.);
		histoVertexZData->GetYaxis()->SetTitleOffset(1.2);
		histoVertexZData->Draw("pe");
		histoVertexZMonteCarlo->Draw("same,histe");
   canvasZVertexPos->Update();
   canvasZVertexPos->SaveAs(Form("%s/ZVertexDistribution.%s",outputDirectory.Data(),suffix.Data()));
// 		delete canvasZVertexPos;

   
   DrawGammaCanvasSettings( canvasZVertexPos, 0.09, 0.015, 0.01, 0.09);
   canvasZVertexPos->cd();
   canvasZVertexPos->SetLogy(1);
      GammaScalingHistogramm(histoNumberOfGoodESDTracksData,1./numberGoodEventsData);
      GammaScalingHistogramm(histoNumberOfGoodESDTracksMonteCarlo,1./numberGoodEventsMonteCarlo);
            DrawAutoGammaHistos( histoNumberOfGoodESDTracksData, 
                                    histoNumberOfGoodESDTracksMonteCarlo, 
                                    "","# good tracks #eta < 0.9","counts normalized",
                                    kTRUE, 1.5,2.e-8,
                                    kFALSE,0. ,0.,
                                    kFALSE, histoNumberOfGoodESDTracksData->GetMean()-3* histoNumberOfGoodESDTracksData->GetRMS(),histoNumberOfGoodESDTracksData->GetMean()+3* histoNumberOfGoodESDTracksData->GetRMS());
			
			histoNumberOfGoodESDTracksData->GetYaxis()->SetTitleOffset(1.1);
			histoNumberOfGoodESDTracksData->Draw("pe");
			histoNumberOfGoodESDTracksMonteCarlo->Draw("same,histe");

            cout << histoNumberOfGoodESDTracksData->GetMean() << endl;
            cout << histoNumberOfGoodESDTracksData->GetRMS() << endl;
   canvasZVertexPos->Update();
   canvasZVertexPos->SaveAs(Form("%s/ESDTrackDistribution.%s",outputDirectory.Data(),suffix.Data()));
   delete canvasZVertexPos;
  
      
      
	TFile* fileMappingOutputData = new TFile(Form("MappingOutputData_%s.root",cutSel.Data()),"UPDATE");
	if (textGenerator.CompareTo("Phojet")==0){
		histoConversionPointVsRData->GetXaxis()->SetRangeUser(0.,180.);
		histoConversionPointVsRData->Write(Form("RDistribution_Data_%s",textPeriod.Data()),TObject::kOverwrite);
		histoConversionPointVsZData->Write(Form("ZDistribution_Data_%s",textPeriod.Data()),TObject::kOverwrite);
	}
	histoConversionPointVsRMonteCarlo->GetXaxis()->SetRangeUser(0.,180.);
// 	histoMCConversionVsR->Write(Form("MC_ConversionMapping_R_%s_%s",textGenerator.Data(),textPeriod.Data()),TObject::kOverwrite);
// 	histoMCConversionVsZ->Write(Form("MC_ConversionMapping_Z_%s_%s",textGenerator.Data(),textPeriod.Data()),TObject::kOverwrite);
	histoConversionPointVsRMonteCarlo->Write(Form("RDistribution_MC_%s_%s",textGenerator.Data(),textPeriod.Data()),TObject::kOverwrite);
	cout << "Here" << endl;
	histoConversionPointVsZMonteCarlo->Write(Form("ZDistribution_MC_%s_%s",textGenerator.Data(),textPeriod.Data()),TObject::kOverwrite);
	cout << "Here" << endl;
	graphIntegratedRData->Write(Form("graphIntegratedR_Data_%s_%s",textGenerator.Data(),textPeriod.Data()),TObject::kOverwrite);
	cout << "Here" << endl;
	graphIntegratedRMonteCarlo->Write(Form("graphIntegratedR_MC_%s_%s",textGenerator.Data(),textPeriod.Data()),TObject::kOverwrite);
	cout << "Here" << endl;
	graphIntegratedMidPtRData->Write(Form("graphIntegratedMidPtR_Data_%s_%s",textGenerator.Data(),textPeriod.Data()),TObject::kOverwrite);
	cout << "Here" << endl;
	graphIntegratedMidPtRMonteCarlo->Write(Form("graphIntegratedMidPtR_MC_%s_%s",textGenerator.Data(),textPeriod.Data()),TObject::kOverwrite);
	cout << "Here" << endl;
	graphRelativeDevR->Write(Form("graphRelativeDevR_%s_%s",textGenerator.Data(),textPeriod.Data()),TObject::kOverwrite);
	cout << "Here" << endl;
	graphRelativeDevMidPtR->Write(Form("graphRelativeDevMidPtR_%s_%s",textGenerator.Data(),textPeriod.Data()),TObject::kOverwrite);
	cout << "Here" << endl;
	graphTotalErrors->Write(Form("graphTotalErrors_%s_%s",textGenerator.Data(),textPeriod.Data()),TObject::kOverwrite);
	cout << "Here" << endl;
	
	fileMappingOutputData->Write();
	fileMappingOutputData->Close();

	delete graphIntegratedRData;
	delete graphIntegratedRMonteCarlo;
	delete graphIntegratedMidPtRData;
	delete graphIntegratedMidPtRMonteCarlo;
	delete graphRelativeDevR;
	delete graphRelativeDevMidPtR;
	delete graphTotalErrors;

	TH2F *histoConversionPointVsZRSDD1Data=	(TH2F*)directoryGammaConvData->Get("ESD_ConversionMapping_SDD1_ZR");
	TH2F *histoConversionPointVsZPhiSDD1Data=	(TH2F*)directoryGammaConvData->Get("ESD_ConversionMapping_SDD1_ZPhi");
	TH2F *histoConversionPointVsRPhiSDD1Data=	(TH2F*)directoryGammaConvData->Get("ESD_ConversionMapping_SDD1_RPhi");
	TH2F *histoConversionPointVsZRSDD1MC=	(TH2F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMapping_SDD1_ZR");
	TH2F *histoConversionPointVsZPhiSDD1MC=	(TH2F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMapping_SDD1_ZPhi");
	TH2F *histoConversionPointVsRPhiSDD1MC=	(TH2F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMapping_SDD1_RPhi");
	TH1D* histoMappingPhiSDD1_ZSmaller0_Data= histoConversionPointVsZPhiSDD1Data->ProjectionY( "histoMappingPhiSDD1_ZSmaller0_Data", 1, histoConversionPointVsZPhiSDD1Data->GetXaxis()->FindBin(0.));
	GammaScalingHistogramm(histoMappingPhiSDD1_ZSmaller0_Data,normFactorReconstData);
	TH1D* histoMappingPhiSDD1_ZLarger0_Data= histoConversionPointVsZPhiSDD1Data->ProjectionY( "histoMappingPhiSDD1_ZLarger0_Data", histoConversionPointVsZPhiSDD1Data->GetXaxis()->FindBin(0.), histoConversionPointVsZPhiSDD1Data->GetNbinsX() );
	GammaScalingHistogramm(histoMappingPhiSDD1_ZLarger0_Data,normFactorReconstData);
	TH1D* histoMappingPhiSDD1_ZSmaller0_MC= histoConversionPointVsZPhiSDD1MC->ProjectionY( "histoMappingPhiSDD1_ZSmaller0_MC", 1, histoConversionPointVsZPhiSDD1MC->GetXaxis()->FindBin(0.));
	GammaScalingHistogramm(histoMappingPhiSDD1_ZSmaller0_MC,normFactorReconstMonteCarlo);
	TH1D* histoMappingPhiSDD1_ZLarger0_MC= histoConversionPointVsZPhiSDD1MC->ProjectionY( "histoMappingPhiSDD1_ZLarger0_MC", histoConversionPointVsZPhiSDD1MC->GetXaxis()->FindBin(0.), histoConversionPointVsZPhiSDD1MC->GetNbinsX() );
	GammaScalingHistogramm(histoMappingPhiSDD1_ZLarger0_MC,normFactorReconstMonteCarlo);
	
	TH1D* histoMappingZLowestLadderSDD1_Data= histoConversionPointVsZPhiSDD1Data->ProjectionX( "histoMappingZLowestLadderSDD1_Data", 1, histoConversionPointVsZPhiSDD1Data->GetYaxis()->FindBin(0.4));
	GammaScalingHistogramm(histoMappingZLowestLadderSDD1_Data,normFactorReconstData);
	TH1D* histoMappingZLowestLadderSDD1_MC= histoConversionPointVsZPhiSDD1MC->ProjectionX( "histoMappingZLowestLadderSDD1_MC", 1, histoConversionPointVsZPhiSDD1MC->GetYaxis()->FindBin(0.4));
	GammaScalingHistogramm(histoMappingZLowestLadderSDD1_MC,normFactorReconstMonteCarlo);


	GammaScalingHistogramm(histoConversionPointVsZRSDD1Data,normFactorReconstData);
	GammaScalingHistogramm(histoConversionPointVsZPhiSDD1Data,normFactorReconstData);
	GammaScalingHistogramm(histoConversionPointVsRPhiSDD1Data,normFactorReconstData);
	GammaScalingHistogramm(histoConversionPointVsRPhiSDD1MC,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoConversionPointVsZRSDD1MC,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoConversionPointVsZPhiSDD1MC,normFactorReconstMonteCarlo);
	
	Double_t minimumPlots2D = normFactorReconstData;
	if (normFactorReconstData > normFactorReconstMonteCarlo) minimumPlots2D = normFactorReconstMonteCarlo;
	
	cout << normFactorReconstData << "\t"<< normFactorReconstMonteCarlo << endl;
	
	TLatex* latexRRangeSDD1_2 = 	new TLatex(0.15,0.88,"SDD1: 12 < R (cm) < 20"); 
	latexRRangeSDD1_2->SetNDC();
	latexRRangeSDD1_2->SetTextFont(62);
	latexRRangeSDD1_2->SetTextSize(0.04);
	latexRRangeSDD1_2->SetLineWidth(6);      
	TLatex* latexZRangeSDD1_2 = 	new TLatex(0.25,0.84,"-#infty < Z (cm) < 0"); 
	latexZRangeSDD1_2->SetNDC();
	latexZRangeSDD1_2->SetTextFont(62);
	latexZRangeSDD1_2->SetTextSize(0.04);
	latexZRangeSDD1_2->SetLineWidth(6);      

	TLatex* latexRRangeSDD1_1 = 	new TLatex(0.15,0.88,"SDD1: 12 < R (cm) < 20");
	latexRRangeSDD1_1->SetNDC();
	latexRRangeSDD1_1->SetTextFont(62);
	latexRRangeSDD1_1->SetTextSize(0.04);
	latexRRangeSDD1_1->SetLineWidth(6);      
	TLatex* latexZRangeSDD1_1 = 	new TLatex(0.25,0.84,"0 < Z (cm) < #infty");
	latexZRangeSDD1_1->SetNDC();
	latexZRangeSDD1_1->SetTextFont(62);
	latexZRangeSDD1_1->SetTextSize(0.04);
	latexZRangeSDD1_1->SetLineWidth(6);      
	TLatex* latexPhiRangeSDD1_1 = 	new TLatex(0.25,0.84,"0 < #varphi (rad) < 0.4");
	latexPhiRangeSDD1_1->SetNDC();
	latexPhiRangeSDD1_1->SetTextFont(62);
	latexPhiRangeSDD1_1->SetTextSize(0.04);
	latexPhiRangeSDD1_1->SetLineWidth(6);      
	
	TCanvas * canvasPhiSDD1 = new TCanvas("canvasPhiSDD1","",1200,1000);  // gives the page size
	canvasPhiSDD1->SetLogy(1);
	canvasPhiSDD1->cd();
	canvasPhiSDD1->SetLeftMargin(0.12);
		DrawAutoGammaHistosMaterial( histoMappingPhiSDD1_ZSmaller0_Data, 
							histoMappingPhiSDD1_ZSmaller0_MC, 
							"","#varphi (rad) ","Normalized counts",
							kFALSE, 2.5,2.e-10,
							kFALSE,0. ,0.,
							kFALSE, 0.,0.);
		latexEtaRange->SetTextColor(kBlack);
		latexEtaRange->Draw();
		latexRRangeSDD1_2->Draw();
		latexZRangeSDD1_2->Draw();
	canvasPhiSDD1->Update();
	canvasPhiSDD1->SaveAs(Form("%s/PhiSDD1_ZSmaller0.%s",outputDirectory.Data(),suffix.Data()));
		DrawAutoGammaHistosMaterial( histoMappingPhiSDD1_ZLarger0_Data, 
							histoMappingPhiSDD1_ZLarger0_MC, 
							"","#varphi (rad) ","Normalized counts",
							kFALSE, 2.5,2.e-10,
							kFALSE,0. ,0.,
							kFALSE, 0.,0.);
		latexEtaRange->Draw();
		latexRRangeSDD1_1->Draw();
		latexZRangeSDD1_1->Draw();
	canvasPhiSDD1->Update();
	canvasPhiSDD1->SaveAs(Form("%s/PhiSDD1_ZLarger0.%s",outputDirectory.Data(),suffix.Data()));
		DrawAutoGammaHistosMaterial( histoMappingZLowestLadderSDD1_Data, 
							histoMappingZLowestLadderSDD1_MC, 
							"","R (cm) ","Normalized counts",
							kTRUE, 10,2.e-9,
							kFALSE,0. ,0.,
							kFALSE, 0.,0.);
		latexEtaRange->Draw();
		latexRRangeSDD1_1->Draw();
		latexPhiRangeSDD1_1->Draw();
	canvasPhiSDD1->Update();
	canvasPhiSDD1->SaveAs(Form("%s/ZSDD1LowestLadder.%s",outputDirectory.Data(),suffix.Data()));


	delete canvasPhiSDD1;

	
	TLatex* latexRRangeSDD1 = 	new TLatex(0.15,0.88,"SDD1: 12 < R (cm) < 20"); 
	latexRRangeSDD1->SetNDC();
	latexRRangeSDD1->SetTextFont(62);
	latexRRangeSDD1->SetTextSize(0.04);
	latexRRangeSDD1->SetLineWidth(6);      

	TCanvas * canvas2DimensionZRSDD1 = new TCanvas("canvas2DimensionZRSDD1","",500,440);  // gives the page size
	canvas2DimensionZRSDD1->SetLogz(1);     
	canvas2DimensionZRSDD1->SetTopMargin(0.02);                                              
	canvas2DimensionZRSDD1->SetRightMargin(rightMargin);
	canvas2DimensionZRSDD1->SetLeftMargin(leftMargin);
	canvas2DimensionZRSDD1->SetBottomMargin(0.08);
	canvas2DimensionZRSDD1->cd();
	histoConversionPointVsZRSDD1Data->GetZaxis()->SetRangeUser(minimumPlots2D*10,histoConversionPointVsZRSDD1Data->GetMaximum());
	DrawAutoGammaHisto2D(  histoConversionPointVsZRSDD1Data,
									"", "Z (cm)", "R (cm)", "",
									kFALSE, 0., 100.,
									kFALSE, -150., 150.);	
	latexEtaRange->SetTextColor(kBlack);
	latexEtaRange->Draw();
	latexRRangeSDD1->Draw();
	DrawLabelsEvents(floatLocationRight[0],floatLocationRight[1],floatLocationRight[2], 0.00, collisionSystem, "Data", textPeriod);
	
	canvas2DimensionZRSDD1->Update();
	canvas2DimensionZRSDD1->SaveAs(Form("%s/ZRSDD1_Data_distribution.%s",outputDirectory.Data(),suffix.Data()));

	histoConversionPointVsZRSDD1MC->GetZaxis()->SetRangeUser(minimumPlots2D*10,histoConversionPointVsZRSDD1Data->GetMaximum());
	DrawAutoGammaHisto2D(  histoConversionPointVsZRSDD1MC,
									"", "Z (cm)", "R (cm)", "",
									kFALSE, 0., 100.,
									kFALSE, -150., 150.);	
	latexEtaRange->SetTextColor(kBlack);
	latexEtaRange->Draw();
	latexRRangeSDD1->Draw();
	DrawLabelsEvents(floatLocationRight[0],floatLocationRight[1],floatLocationRight[2], 0.00, collisionSystem, textGenerator, textPeriod);
	
	canvas2DimensionZRSDD1->Update();
	canvas2DimensionZRSDD1->SaveAs(Form("%s/ZRSDD1_MC_distribution.%s",outputDirectory.Data(),suffix.Data()));
	delete canvas2DimensionZRSDD1;
	
	TCanvas * canvas2DimensionRPhiSDD1 = NULL;
	canvas2DimensionRPhiSDD1 =new TCanvas("canvas2DimensionRPhiSDD1","",500,440);  // gives the page size
	canvas2DimensionRPhiSDD1->SetLogz(1);     
	canvas2DimensionRPhiSDD1->SetTopMargin(0.02);                                              
	canvas2DimensionRPhiSDD1->SetRightMargin(rightMargin);                                            
	canvas2DimensionRPhiSDD1->SetLeftMargin(leftMargin);
	canvas2DimensionRPhiSDD1->SetBottomMargin(0.08);		
	canvas2DimensionRPhiSDD1->cd();
	histoConversionPointVsRPhiSDD1Data->GetZaxis()->SetRangeUser(minimumPlots2D,histoConversionPointVsRPhiSDD1Data->GetMaximum());
	DrawAutoGammaHisto2D(  histoConversionPointVsRPhiSDD1Data,
									"", "R (cm)", "#varphi (cm)", "",
									kFALSE, 0., 100.,
									kFALSE, 0., 180.);
	latexEtaRange->Draw();
	latexRRangeSDD1->Draw();
	DrawLabelsEvents(0.15,floatLocationRight[1],floatLocationRight[2], 0.00, collisionSystem, "Data", textPeriod);
	canvas2DimensionRPhiSDD1->Update();
	canvas2DimensionRPhiSDD1->SaveAs(Form("%s/RPhiSDD1_Data_distribution.%s",outputDirectory.Data(),suffix.Data()));

	histoConversionPointVsRPhiSDD1MC->GetZaxis()->SetRangeUser(minimumPlots2D,histoConversionPointVsRPhiSDD1Data->GetMaximum());
	DrawAutoGammaHisto2D(  histoConversionPointVsRPhiSDD1MC,
									"", "R (cm)", "#varphi (cm)", "",
									kFALSE, 0., 100.,
									kFALSE, 0., 180.);
	latexEtaRange->Draw();
	latexRRangeSDD1->Draw();
	DrawLabelsEvents(0.15,floatLocationRight[1],floatLocationRight[2], 0.00, collisionSystem, textGenerator, textPeriod);
	canvas2DimensionRPhiSDD1->Update();
	canvas2DimensionRPhiSDD1->SaveAs(Form("%s/RPhiSDD1_MC_distribution.%s",outputDirectory.Data(),suffix.Data()));

	delete canvas2DimensionRPhiSDD1;
	
	
	TCanvas * canvas2DimensionZPhiSDD1 = NULL;  // gives the page size
	canvas2DimensionZPhiSDD1 =new TCanvas("canvas2DimensionZPhiSDD1","",500,440); 
	canvas2DimensionZPhiSDD1->SetLogz(1);     
	canvas2DimensionZPhiSDD1->SetTopMargin(0.02);                                              
	canvas2DimensionZPhiSDD1->SetRightMargin(rightMargin);                                            
	canvas2DimensionZPhiSDD1->SetLeftMargin(leftMargin);
	canvas2DimensionZPhiSDD1->SetBottomMargin(0.08);		
	canvas2DimensionZPhiSDD1->cd();
	histoConversionPointVsZPhiSDD1Data->GetZaxis()->SetRangeUser(minimumPlots2D,histoConversionPointVsZPhiSDD1Data->GetMaximum());
	DrawAutoGammaHisto2D(  histoConversionPointVsZPhiSDD1Data,
									"", "Z (cm)", "#varphi (cm)", "",
									kFALSE, 0., 100.,
									kFALSE, -200., 200.);
	latexEtaRange->Draw();
	latexRRangeSDD1->Draw();
	DrawLabelsEvents(0.15,floatLocationRight[1],floatLocationRight[2], 0.00, collisionSystem, "Data", textPeriod);
	
	canvas2DimensionZPhiSDD1->Update();
	canvas2DimensionZPhiSDD1->SaveAs(Form("%s/ZPhiSDD1_Data_distribution.%s",outputDirectory.Data(),suffix.Data()));

	histoConversionPointVsZPhiSDD1MC->GetZaxis()->SetRangeUser(minimumPlots2D,histoConversionPointVsZPhiSDD1Data->GetMaximum());
	DrawAutoGammaHisto2D(  histoConversionPointVsZPhiSDD1MC,
									"", "Z (cm)", "#varphi (cm)", "",
									kFALSE, 0., 100.,
									kFALSE, -200., 200.);
	latexEtaRange->Draw();
	latexRRangeSDD1->Draw();
	DrawLabelsEvents(0.15,floatLocationRight[1],floatLocationRight[2], 0.00, collisionSystem, textGenerator, textPeriod);
	
	canvas2DimensionZPhiSDD1->Update();
	canvas2DimensionZPhiSDD1->SaveAs(Form("%s/ZPhiSDD1_MC_distribution.%s",outputDirectory.Data(),suffix.Data()));
	delete canvas2DimensionZPhiSDD1;
	
	TH2F *histoConversionPointVsZRSDD2Data=	(TH2F*)directoryGammaConvData->Get("ESD_ConversionMapping_SDD2_ZR");
	TH2F *histoConversionPointVsZPhiSDD2Data=	(TH2F*)directoryGammaConvData->Get("ESD_ConversionMapping_SDD2_ZPhi");
	TH2F *histoConversionPointVsRPhiSDD2Data=	(TH2F*)directoryGammaConvData->Get("ESD_ConversionMapping_SDD2_RPhi");
	TH2F *histoConversionPointVsZRSDD2MC=	(TH2F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMapping_SDD2_ZR");
	TH2F *histoConversionPointVsZPhiSDD2MC=	(TH2F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMapping_SDD2_ZPhi");
	TH2F *histoConversionPointVsRPhiSDD2MC=	(TH2F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMapping_SDD2_RPhi");
	TH1D* histoMappingPhiSDD2_ZSmaller0_Data= histoConversionPointVsZPhiSDD2Data->ProjectionY( "histoMappingPhiSDD2_ZSmaller0_Data", 1, histoConversionPointVsZPhiSDD2Data->GetXaxis()->FindBin(0.));
	GammaScalingHistogramm(histoMappingPhiSDD2_ZSmaller0_Data,normFactorReconstData);
	TH1D* histoMappingPhiSDD2_ZLarger0_Data= histoConversionPointVsZPhiSDD2Data->ProjectionY( "histoMappingPhiSDD2_ZLarger0_Data", histoConversionPointVsZPhiSDD2Data->GetXaxis()->FindBin(0.), histoConversionPointVsZPhiSDD2Data->GetNbinsX() );
	GammaScalingHistogramm(histoMappingPhiSDD2_ZLarger0_Data,normFactorReconstData);
	TH1D* histoMappingPhiSDD2_ZSmaller0_MC= histoConversionPointVsZPhiSDD2MC->ProjectionY( "histoMappingPhiSDD2_ZSmaller0_MC", 1, histoConversionPointVsZPhiSDD2MC->GetXaxis()->FindBin(0.));
	GammaScalingHistogramm(histoMappingPhiSDD2_ZSmaller0_MC,normFactorReconstMonteCarlo);
	TH1D* histoMappingPhiSDD2_ZLarger0_MC= histoConversionPointVsZPhiSDD2MC->ProjectionY( "histoMappingPhiSDD2_ZLarger0_MC", histoConversionPointVsZPhiSDD2MC->GetXaxis()->FindBin(0.), histoConversionPointVsZPhiSDD2MC->GetNbinsX() );
	GammaScalingHistogramm(histoMappingPhiSDD2_ZLarger0_MC,normFactorReconstMonteCarlo);
	
	TH1D* histoMappingZLowestLadderSDD2_Data= histoConversionPointVsZPhiSDD2Data->ProjectionX( "histoMappingZLowestLadderSDD2_Data", 1, histoConversionPointVsZPhiSDD2Data->GetYaxis()->FindBin(0.4));
	GammaScalingHistogramm(histoMappingZLowestLadderSDD2_Data,normFactorReconstData);
	TH1D* histoMappingZLowestLadderSDD2_MC= histoConversionPointVsZPhiSDD2MC->ProjectionX( "histoMappingZLowestLadderSDD2_MC", 1, histoConversionPointVsZPhiSDD2MC->GetYaxis()->FindBin(0.4));
	GammaScalingHistogramm(histoMappingZLowestLadderSDD2_MC,normFactorReconstMonteCarlo);

	GammaScalingHistogramm(histoConversionPointVsZRSDD2Data,normFactorReconstData);
	GammaScalingHistogramm(histoConversionPointVsZPhiSDD2Data,normFactorReconstData);
	GammaScalingHistogramm(histoConversionPointVsRPhiSDD2Data,normFactorReconstData);
	GammaScalingHistogramm(histoConversionPointVsRPhiSDD2MC,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoConversionPointVsZRSDD2MC,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoConversionPointVsZPhiSDD2MC,normFactorReconstMonteCarlo);
	

	TLatex* latexRRangeSDD2_2 = 	new TLatex(0.15,0.88,"SDD2: 20 < R (cm) < 28"); 
	latexRRangeSDD2_2->SetNDC();
	latexRRangeSDD2_2->SetTextFont(62);
	latexRRangeSDD2_2->SetTextSize(0.04);
	latexRRangeSDD2_2->SetLineWidth(6);      
	TLatex* latexZRangeSDD2_2 = 	new TLatex(0.25,0.84,"-#infty < Z (cm) < 0"); 
	latexZRangeSDD2_2->SetNDC();
	latexZRangeSDD2_2->SetTextFont(62);
	latexZRangeSDD2_2->SetTextSize(0.04);
	latexZRangeSDD2_2->SetLineWidth(6);      

	TLatex* latexRRangeSDD2_1 = 	new TLatex(0.15,0.88,"SDD2: 20 < R (cm) < 28");
	latexRRangeSDD2_1->SetNDC();
	latexRRangeSDD2_1->SetTextFont(62);
	latexRRangeSDD2_1->SetTextSize(0.04);
	latexRRangeSDD2_1->SetLineWidth(6);      
	TLatex* latexZRangeSDD2_1 = 	new TLatex(0.25,0.84,"0 < Z (cm) < #infty");
	latexZRangeSDD2_1->SetNDC();
	latexZRangeSDD2_1->SetTextFont(62);
	latexZRangeSDD2_1->SetTextSize(0.04);
	latexZRangeSDD2_1->SetLineWidth(6);      
	TLatex* latexPhiRangeSDD2_1 = 	new TLatex(0.25,0.84,"0 < #varphi (rad) < 0.4");
	latexPhiRangeSDD2_1->SetNDC();
	latexPhiRangeSDD2_1->SetTextFont(62);
	latexPhiRangeSDD2_1->SetTextSize(0.04);
	latexPhiRangeSDD2_1->SetLineWidth(6);      
	
	TCanvas * canvasPhiSDD2 = new TCanvas("canvasPhiSDD2","",1200,1000);  // gives the page size
	canvasPhiSDD2->SetLogy(1);
	canvasPhiSDD2->cd();
	canvasPhiSDD2->SetLeftMargin(0.12);
		DrawAutoGammaHistosMaterial( histoMappingPhiSDD2_ZSmaller0_Data, 
							histoMappingPhiSDD2_ZSmaller0_MC, 
							"","#varphi (rad) ","Normalized counts",
							kFALSE, 2.5,2.e-10,
							kFALSE,0. ,0.,
							kFALSE, 0.,0.);
		latexEtaRange->SetTextColor(kBlack);
		latexEtaRange->Draw();
		latexRRangeSDD2_2->Draw();
		latexZRangeSDD2_2->Draw();
	canvasPhiSDD2->Update();
	canvasPhiSDD2->SaveAs(Form("%s/PhiSDD2_ZSmaller0.%s",outputDirectory.Data(),suffix.Data()));
		DrawAutoGammaHistosMaterial( histoMappingPhiSDD2_ZLarger0_Data, 
							histoMappingPhiSDD2_ZLarger0_MC, 
							"","#varphi (rad) ","Normalized counts",
							kFALSE, 2.5,2.e-10,
							kFALSE,0. ,0.,
							kFALSE, 0.,0.);
		latexEtaRange->Draw();
		latexRRangeSDD2_1->Draw();
		latexZRangeSDD2_1->Draw();
	canvasPhiSDD2->Update();
	canvasPhiSDD2->SaveAs(Form("%s/PhiSDD2_ZLarger0.%s",outputDirectory.Data(),suffix.Data()));
		DrawAutoGammaHistosMaterial( histoMappingZLowestLadderSDD2_Data, 
							histoMappingZLowestLadderSDD2_MC, 
							"","R (cm) ","Normalized counts",
							kTRUE, 10,2.e-9,
							kFALSE,0. ,0.,
							kFALSE, 0.,0.);
		latexEtaRange->Draw();
		latexRRangeSDD2_1->Draw();
		latexPhiRangeSDD2_1->Draw();
	canvasPhiSDD2->Update();
	canvasPhiSDD2->SaveAs(Form("%s/ZSDD2LowestLadder.%s",outputDirectory.Data(),suffix.Data()));


	delete canvasPhiSDD2;
	
	TLatex* latexRRangeSDD2 = 	new TLatex(0.15,0.88,"SDD2: 20 < R (cm) < 28"); 
	latexRRangeSDD2->SetNDC();
	latexRRangeSDD2->SetTextFont(62);
	latexRRangeSDD2->SetTextSize(0.04);
	latexRRangeSDD2->SetLineWidth(6);      

	TCanvas * canvas2DimensionZRSDD2 = new TCanvas("canvas2DimensionZRSDD2","",500,440);  // gives the page size
	canvas2DimensionZRSDD2->SetLogz(1);     
	canvas2DimensionZRSDD2->SetTopMargin(0.02);                                              
	canvas2DimensionZRSDD2->SetRightMargin(rightMargin);
	canvas2DimensionZRSDD2->SetLeftMargin(leftMargin);
	canvas2DimensionZRSDD2->SetBottomMargin(0.08);
	canvas2DimensionZRSDD2->cd();
	histoConversionPointVsZRSDD2Data->GetZaxis()->SetRangeUser(minimumPlots2D*10,histoConversionPointVsZRSDD2Data->GetMaximum());
	DrawAutoGammaHisto2D(  histoConversionPointVsZRSDD2Data,
									"", "Z (cm)", "R (cm)", "",
									kFALSE, 0., 100.,
									kFALSE, -150., 150.);	
	latexEtaRange->SetTextColor(kBlack);
	latexEtaRange->Draw();
	latexRRangeSDD2->Draw();
	DrawLabelsEvents(floatLocationRight[0],floatLocationRight[1],floatLocationRight[2], 0.00, collisionSystem, "Data", textPeriod);
	
	canvas2DimensionZRSDD2->Update();
	canvas2DimensionZRSDD2->SaveAs(Form("%s/ZRSDD2_Data_distribution.%s",outputDirectory.Data(),suffix.Data()));

	histoConversionPointVsZRSDD2MC->GetZaxis()->SetRangeUser(minimumPlots2D*10,histoConversionPointVsZRSDD2Data->GetMaximum());
	DrawAutoGammaHisto2D(  histoConversionPointVsZRSDD2MC,
									"", "Z (cm)", "R (cm)", "",
									kFALSE, 0., 100.,
									kFALSE, -150., 150.);	
	latexEtaRange->SetTextColor(kBlack);
	latexEtaRange->Draw();
	latexRRangeSDD2->Draw();
	DrawLabelsEvents(floatLocationRight[0],floatLocationRight[1],floatLocationRight[2], 0.00, collisionSystem, textGenerator, textPeriod);
	
	canvas2DimensionZRSDD2->Update();
	canvas2DimensionZRSDD2->SaveAs(Form("%s/ZRSDD2_MC_distribution.%s",outputDirectory.Data(),suffix.Data()));
	delete canvas2DimensionZRSDD2;
	
	TCanvas * canvas2DimensionRPhiSDD2 = NULL;
	canvas2DimensionRPhiSDD2 =new TCanvas("canvas2DimensionRPhiSDD2","",500,440);  // gives the page size
	canvas2DimensionRPhiSDD2->SetLogz(1);     
	canvas2DimensionRPhiSDD2->SetTopMargin(0.02);                                              
	canvas2DimensionRPhiSDD2->SetRightMargin(rightMargin);                                            
	canvas2DimensionRPhiSDD2->SetLeftMargin(leftMargin);
	canvas2DimensionRPhiSDD2->SetBottomMargin(0.08);		
	canvas2DimensionRPhiSDD2->cd();
	histoConversionPointVsRPhiSDD2Data->GetZaxis()->SetRangeUser(minimumPlots2D,histoConversionPointVsRPhiSDD2Data->GetMaximum());
	DrawAutoGammaHisto2D(  histoConversionPointVsRPhiSDD2Data,
									"", "R (cm)", "#varphi (cm)", "",
									kFALSE, 0., 100.,
									kFALSE, 0., 180.);
	latexEtaRange->Draw();
	latexRRangeSDD2->Draw();
	DrawLabelsEvents(0.15,floatLocationRight[1],floatLocationRight[2], 0.00, collisionSystem, "Data", textPeriod);
	canvas2DimensionRPhiSDD2->Update();
	canvas2DimensionRPhiSDD2->SaveAs(Form("%s/RPhiSDD2_Data_distribution.%s",outputDirectory.Data(),suffix.Data()));

	histoConversionPointVsRPhiSDD2MC->GetZaxis()->SetRangeUser(minimumPlots2D,histoConversionPointVsRPhiSDD2Data->GetMaximum());
	DrawAutoGammaHisto2D(  histoConversionPointVsRPhiSDD2MC,
									"", "R (cm)", "#varphi (cm)", "",
									kFALSE, 0., 100.,
									kFALSE, 0., 180.);
	latexEtaRange->Draw();
	latexRRangeSDD2->Draw();
	DrawLabelsEvents(0.15,floatLocationRight[1],floatLocationRight[2], 0.00, collisionSystem, textGenerator, textPeriod);
	canvas2DimensionRPhiSDD2->Update();
	canvas2DimensionRPhiSDD2->SaveAs(Form("%s/RPhiSDD2_MC_distribution.%s",outputDirectory.Data(),suffix.Data()));

	delete canvas2DimensionRPhiSDD2;
	
	
	TCanvas * canvas2DimensionZPhiSDD2 = NULL;  // gives the page size
	canvas2DimensionZPhiSDD2 =new TCanvas("canvas2DimensionZPhiSDD2","",500,440); 
	canvas2DimensionZPhiSDD2->SetLogz(1);     
	canvas2DimensionZPhiSDD2->SetTopMargin(0.02);                                              
	canvas2DimensionZPhiSDD2->SetRightMargin(rightMargin);                                            
	canvas2DimensionZPhiSDD2->SetLeftMargin(leftMargin);
	canvas2DimensionZPhiSDD2->SetBottomMargin(0.08);		
	canvas2DimensionZPhiSDD2->cd();
	histoConversionPointVsZPhiSDD2Data->GetZaxis()->SetRangeUser(minimumPlots2D,histoConversionPointVsZPhiSDD2Data->GetMaximum());
	DrawAutoGammaHisto2D(  histoConversionPointVsZPhiSDD2Data,
									"", "Z (cm)", "#varphi (cm)", "",
									kFALSE, 0., 100.,
									kFALSE, -200., 200.);
	latexEtaRange->Draw();
	latexRRangeSDD2->Draw();
	DrawLabelsEvents(0.15,floatLocationRight[1],floatLocationRight[2], 0.00, collisionSystem, "Data", textPeriod);
	
	canvas2DimensionZPhiSDD2->Update();
	canvas2DimensionZPhiSDD2->SaveAs(Form("%s/ZPhiSDD2_Data_distribution.%s",outputDirectory.Data(),suffix.Data()));

	histoConversionPointVsZPhiSDD2MC->GetZaxis()->SetRangeUser(minimumPlots2D,histoConversionPointVsZPhiSDD2Data->GetMaximum());
	DrawAutoGammaHisto2D(  histoConversionPointVsZPhiSDD2MC,
									"", "Z (cm)", "#varphi (cm)", "",
									kFALSE, 0., 100.,
									kFALSE, -200., 200.);
	latexEtaRange->Draw();
	latexRRangeSDD2->Draw();
	DrawLabelsEvents(0.15,floatLocationRight[1],floatLocationRight[2], 0.00, collisionSystem, textGenerator, textPeriod);
	
	canvas2DimensionZPhiSDD2->Update();
	canvas2DimensionZPhiSDD2->SaveAs(Form("%s/ZPhiSDD2_MC_distribution.%s",outputDirectory.Data(),suffix.Data()));
	delete canvas2DimensionZPhiSDD2;
	
	TH2F *histoConversionPointVsZRSDDThermalData=	(TH2F*)directoryGammaConvData->Get("ESD_ConversionMapping_SDDThermal_ZR");
	TH2F *histoConversionPointVsZPhiSDDThermalData=	(TH2F*)directoryGammaConvData->Get("ESD_ConversionMapping_SDDThermal_ZPhi");
	TH2F *histoConversionPointVsRPhiSDDThermalData=	(TH2F*)directoryGammaConvData->Get("ESD_ConversionMapping_SDDThermal_RPhi");
	TH2F *histoConversionPointVsZRSDDThermalMC=	(TH2F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMapping_SDDThermal_ZR");
	TH2F *histoConversionPointVsZPhiSDDThermalMC=	(TH2F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMapping_SDDThermal_ZPhi");
	TH2F *histoConversionPointVsRPhiSDDThermalMC=	(TH2F*)directoryGammaConvMonteCarlo->Get("ESD_ConversionMapping_SDDThermal_RPhi");
	
	
	GammaScalingHistogramm(histoConversionPointVsZRSDDThermalData,normFactorReconstData);
	GammaScalingHistogramm(histoConversionPointVsZPhiSDDThermalData,normFactorReconstData);
	GammaScalingHistogramm(histoConversionPointVsRPhiSDDThermalData,normFactorReconstData);
	GammaScalingHistogramm(histoConversionPointVsRPhiSDDThermalMC,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoConversionPointVsZRSDDThermalMC,normFactorReconstMonteCarlo);
	GammaScalingHistogramm(histoConversionPointVsZPhiSDDThermalMC,normFactorReconstMonteCarlo);
	
	TLatex* latexRRangeSDDThermal = 	new TLatex(0.15,0.88,"SDD Thermal Shield: 28 < R (cm) < 32"); 
	latexRRangeSDDThermal->SetNDC();
	latexRRangeSDDThermal->SetTextFont(62);
	latexRRangeSDDThermal->SetTextSize(0.04);
	latexRRangeSDDThermal->SetLineWidth(6);      

	TCanvas * canvas2DimensionZRSDDThermal = new TCanvas("canvas2DimensionZRSDDThermal","",500,440);  // gives the page size
	canvas2DimensionZRSDDThermal->SetLogz(1);     
	canvas2DimensionZRSDDThermal->SetTopMargin(0.02);                                              
	canvas2DimensionZRSDDThermal->SetRightMargin(rightMargin);
	canvas2DimensionZRSDDThermal->SetLeftMargin(leftMargin);
	canvas2DimensionZRSDDThermal->SetBottomMargin(0.08);
	canvas2DimensionZRSDDThermal->cd();
	histoConversionPointVsZRSDDThermalData->GetZaxis()->SetRangeUser(minimumPlots2D*10,histoConversionPointVsZRSDDThermalData->GetMaximum());
	DrawAutoGammaHisto2D(  histoConversionPointVsZRSDDThermalData,
									"", "Z (cm)", "R (cm)", "",
									kFALSE, 0., 100.,
									kFALSE, -150., 150.);	
	latexEtaRange->SetTextColor(kBlack);
	latexEtaRange->Draw();
	latexRRangeSDDThermal->Draw();
	DrawLabelsEvents(floatLocationRight[0],floatLocationRight[1],floatLocationRight[2], 0.00, collisionSystem, "Data", textPeriod);
	
	canvas2DimensionZRSDDThermal->Update();
	canvas2DimensionZRSDDThermal->SaveAs(Form("%s/ZRSDDThermal_Data_distribution.%s",outputDirectory.Data(),suffix.Data()));

	histoConversionPointVsZRSDDThermalMC->GetZaxis()->SetRangeUser(minimumPlots2D*10,histoConversionPointVsZRSDDThermalData->GetMaximum());
	DrawAutoGammaHisto2D(  histoConversionPointVsZRSDDThermalMC,
									"", "Z (cm)", "R (cm)", "",
									kFALSE, 0., 100.,
									kFALSE, -150., 150.);	
	latexEtaRange->SetTextColor(kBlack);
	latexEtaRange->Draw();
	latexRRangeSDDThermal->Draw();
	DrawLabelsEvents(floatLocationRight[0],floatLocationRight[1],floatLocationRight[2], 0.00, collisionSystem, textGenerator, textPeriod);
	
	canvas2DimensionZRSDDThermal->Update();
	canvas2DimensionZRSDDThermal->SaveAs(Form("%s/ZRSDDThermal_MC_distribution.%s",outputDirectory.Data(),suffix.Data()));
	delete canvas2DimensionZRSDDThermal;
	
	TCanvas * canvas2DimensionRPhiSDDThermal = NULL;
	canvas2DimensionRPhiSDDThermal =new TCanvas("canvas2DimensionRPhiSDDThermal","",500,440);  // gives the page size
	canvas2DimensionRPhiSDDThermal->SetLogz(1);     
	canvas2DimensionRPhiSDDThermal->SetTopMargin(0.02);                                              
	canvas2DimensionRPhiSDDThermal->SetRightMargin(rightMargin);                                            
	canvas2DimensionRPhiSDDThermal->SetLeftMargin(leftMargin);
	canvas2DimensionRPhiSDDThermal->SetBottomMargin(0.08);		
	canvas2DimensionRPhiSDDThermal->cd();
	histoConversionPointVsRPhiSDDThermalData->GetZaxis()->SetRangeUser(minimumPlots2D,histoConversionPointVsRPhiSDDThermalData->GetMaximum());
	DrawAutoGammaHisto2D(  histoConversionPointVsRPhiSDDThermalData,
									"", "R (cm)", "#varphi (cm)", "",
									kFALSE, 0., 100.,
									kFALSE, 0., 180.);
	latexEtaRange->Draw();
	latexRRangeSDDThermal->Draw();
	DrawLabelsEvents(0.15,floatLocationRight[1],floatLocationRight[2], 0.00, collisionSystem, "Data", textPeriod);
	canvas2DimensionRPhiSDDThermal->Update();
	canvas2DimensionRPhiSDDThermal->SaveAs(Form("%s/RPhiSDDThermal_Data_distribution.%s",outputDirectory.Data(),suffix.Data()));

	histoConversionPointVsRPhiSDDThermalMC->GetZaxis()->SetRangeUser(minimumPlots2D,histoConversionPointVsRPhiSDDThermalData->GetMaximum());
	DrawAutoGammaHisto2D(  histoConversionPointVsRPhiSDDThermalMC,
									"", "R (cm)", "#varphi (cm)", "",
									kFALSE, 0., 100.,
									kFALSE, 0., 180.);
	latexEtaRange->Draw();
	latexRRangeSDDThermal->Draw();
	DrawLabelsEvents(0.15,floatLocationRight[1],floatLocationRight[2], 0.00, collisionSystem, textGenerator, textPeriod);
	canvas2DimensionRPhiSDDThermal->Update();
	canvas2DimensionRPhiSDDThermal->SaveAs(Form("%s/RPhiSDDThermal_MC_distribution.%s",outputDirectory.Data(),suffix.Data()));

	delete canvas2DimensionRPhiSDDThermal;
	
	
	TCanvas * canvas2DimensionZPhiSDDThermal = NULL;  // gives the page size
	canvas2DimensionZPhiSDDThermal =new TCanvas("canvas2DimensionZPhiSDDThermal","",500,440); 
	canvas2DimensionZPhiSDDThermal->SetLogz(1);     
	canvas2DimensionZPhiSDDThermal->SetTopMargin(0.02);                                              
	canvas2DimensionZPhiSDDThermal->SetRightMargin(rightMargin);                                            
	canvas2DimensionZPhiSDDThermal->SetLeftMargin(leftMargin);
	canvas2DimensionZPhiSDDThermal->SetBottomMargin(0.08);		
	canvas2DimensionZPhiSDDThermal->cd();
	histoConversionPointVsZPhiSDDThermalData->GetZaxis()->SetRangeUser(minimumPlots2D,histoConversionPointVsZPhiSDDThermalData->GetMaximum());
	DrawAutoGammaHisto2D(  histoConversionPointVsZPhiSDDThermalData,
									"", "Z (cm)", "#varphi (cm)", "",
									kFALSE, 0., 100.,
									kFALSE, -200., 200.);
	latexEtaRange->Draw();
	latexRRangeSDDThermal->Draw();
	DrawLabelsEvents(0.15,floatLocationRight[1],floatLocationRight[2], 0.00, collisionSystem, "Data", textPeriod);
	
	canvas2DimensionZPhiSDDThermal->Update();
	canvas2DimensionZPhiSDDThermal->SaveAs(Form("%s/ZPhiSDDThermal_Data_distribution.%s",outputDirectory.Data(),suffix.Data()));

	histoConversionPointVsZPhiSDDThermalMC->GetZaxis()->SetRangeUser(minimumPlots2D,histoConversionPointVsZPhiSDDThermalData->GetMaximum());
	DrawAutoGammaHisto2D(  histoConversionPointVsZPhiSDDThermalMC,
									"", "Z (cm)", "#varphi (cm)", "",
									kFALSE, 0., 100.,
									kFALSE, -200., 200.);
	latexEtaRange->Draw();
	latexRRangeSDDThermal->Draw();
	DrawLabelsEvents(0.15,floatLocationRight[1],floatLocationRight[2], 0.00, collisionSystem, textGenerator, textPeriod);
	
	canvas2DimensionZPhiSDDThermal->Update();
	canvas2DimensionZPhiSDDThermal->SaveAs(Form("%s/ZPhiSDDThermal_MC_distribution.%s",outputDirectory.Data(),suffix.Data()));
	delete canvas2DimensionZPhiSDDThermal;
	
	TFile* fileMappingDetailed = new TFile(Form("%s/%s/MappingDetailed.root",outputDirectory1.Data(),optEnergy.Data()),"UPDATE");
	fileMappingDetailed->mkdir(textPeriod.Data());
	fileMappingDetailed->cd(textPeriod.Data());
	cout << "Here" << endl;
	for(Int_t iR = iRStart; iR < nBinsR; iR++){
		histoMappingPhiInRData[iR]->Write(Form("Data_Conversion_Mapping_Phi_in_R_%02d",iR),TObject::kOverwrite);
		histoMappingPhiInRMonteCarlo[iR]->Write(Form("MC_Conversion_Mapping_Phi_in_R_%02d",iR),TObject::kOverwrite);
		histoMappingZInRData[iR]->Write(Form("Data_Conversion_Mapping_Z_in_R_%02d",iR),TObject::kOverwrite);
		histoMappingZInRMonteCarlo[iR]->Write(Form("MC_Conversion_Mapping_Z_in_R_%02d",iR),TObject::kOverwrite);
	}
	for(Int_t iZ = 0; iZ < nBinsZ; iZ++){
		histoMappingPhiInZData[iZ]->Write(Form("Data_Conversion_Mapping_Phi_in_Z_%02d",iZ),TObject::kOverwrite);
		histoMappingPhiInZMonteCarlo[iZ]->Write(Form("MC_Conversion_Mapping_Phi_in_Z_%02d",iZ),TObject::kOverwrite);
		histoMappingSPDPhiInZData[iZ]->Write(Form("Data_Conversion_Mapping_SPD_Phi_in_Z_%02d",iZ),TObject::kOverwrite);
		histoMappingSPDPhiInZMonteCarlo[iZ]->Write(Form("MC_Conversion_Mapping_SPD_Phi_in_Z_%02d",iZ),TObject::kOverwrite);
		histoMappingSPDThPhiInZData[iZ]->Write(Form("Data_Conversion_Mapping_SPDTh_Phi_in_Z_%02d",iZ),TObject::kOverwrite);
		histoMappingSPDThPhiInZMonteCarlo[iZ]->Write(Form("MC_Conversion_Mapping_SPDTh_Phi_in_Z_%02d",iZ),TObject::kOverwrite);
		histoMappingITSTPCPhiInZData[iZ]->Write(Form("Data_Conversion_Mapping_ITSTPC_Phi_in_Z_%02d",iZ),TObject::kOverwrite);
		histoMappingITSTPCPhiInZMonteCarlo[iZ]->Write(Form("MC_Conversion_Mapping_ITSTPC_Phi_in_Z_%02d",iZ),TObject::kOverwrite);
		histoMappingRInZData[iZ]->Write(Form("Data_Conversion_Mapping_R_in_Z_%02d",iZ),TObject::kOverwrite);
		histoMappingRInZMonteCarlo[iZ]->Write(Form("MC_Conversion_Mapping_R_in_Z_%02d",iZ),TObject::kOverwrite);
	}
	histoConversionPointVsRData->GetXaxis()->SetRangeUser(0.,180.);
	cout << "Here" << endl;
	histoConversionPointVsRData->Write("RDistribution_Data",TObject::kOverwrite);
	cout << "Here" << endl;
	histoConversionPointVsRMonteCarlo->GetXaxis()->SetRangeUser(0.,180.);
	histoConversionPointVsRMonteCarlo->Write(Form("RDistribution_MC_%s",textGenerator.Data()),TObject::kOverwrite);
	histoConversionPointVsZData->Write("ZDistribution_Data",TObject::kOverwrite);
	histoConversionPointVsZMonteCarlo->Write(Form("ZDistribution_MC_%s",textGenerator.Data()),TObject::kOverwrite);
	histoConversionPointVsXYData->Write("XYDistribution_Data",TObject::kOverwrite);
	histoConversionPointVsXYMonteCarlo->Write("XYDistribution_MonteCarlo",TObject::kOverwrite);
	histoConversionPointVsZRData->Write("ZRDistribution_Data",TObject::kOverwrite);
	histoConversionPointVsZRMonteCarlo->Write("ZRDistribution_MonteCarlo",TObject::kOverwrite);
	cout << "Here" << endl;
// 	histoMCConversionVsR->Write("MC_ConversionMapping_R",TObject::kOverwrite);
// 	histoMCConversionVsZ->Write("MC_ConversionMapping_Z",TObject::kOverwrite);
	histoDiffXYDistribution->Write("RatioXYDistribution",TObject::kOverwrite);
	cout << "Here" << endl;
	histoDiffZRDistribution->Write("RatioZRDistribution",TObject::kOverwrite);
	fileMappingDetailed->Write();
	fileMappingDetailed->Close();
}
