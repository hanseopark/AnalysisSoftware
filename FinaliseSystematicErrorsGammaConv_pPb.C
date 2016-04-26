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
#include "TGraphAsymmErrors.h" 
#include "TGaxis.h"
#include "TMarker.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"

void FinalyseSystematicErrorsGammapPb(const char* nameDataFileErrors ="", TString energy="", Int_t nBinsDoubleRatio =1, Int_t nBinsIncRatio =1 , Int_t nBinsPi0 =1 , Int_t nBinsGamma =1, 	Int_t offSetBeginning = 0, TString additionalName = "",TString additionalNameOutput = "",TString suffix = "eps"){
	
	StyleSettingsThesis();	
	SetPlotStyle();
	
	TString date = ReturnDateString();
	TString dateForOutput = ReturnDateStringForOutput();
	
	Int_t color[20] = {860,894,807,880,418,403,802,923,634,432,404,435,420,407,416,830,404,608,920,1};
	
	Int_t 	numberOfEntriesPos = 0;
	Int_t 	numberOfEntriesNeg = 0;
		
	const Int_t ErrorDoubleRatio = 12;
	const Int_t ErrorIncRatio = 9;
	const Int_t ErrorPi0 = 9;
	const Int_t ErrorGammaSpec = 7;

	TString nameCutVariationsDoubleRatio[ErrorDoubleRatio] = {"Chi2","PsiPair","Cluster","dEdxE","dEdxPi","Qt","SinglePt","Alpha","IntRange","CocktailEta","CocktailParam","CocktailEtaNorm"};//,"Fit"};
	TString nameCutVariationsIncRatio[ErrorIncRatio] = {"Chi2","PsiPair","Cluster","dEdxE","dEdxPi","Qt","SinglePt","Alpha","IntRange"};//,"Fit"};
	TString nameCutVariationsPi0[ErrorPi0] = {"Chi2","PsiPair","Cluster","dEdxE","dEdxPi","Qt","SinglePt","Alpha","IntRange"};//,"Fit"};
	TString nameCutVariationsGamma[ErrorGammaSpec] = {"Chi2","PsiPair","Cluster","dEdxE","dEdxPi","Qt","SinglePt"};
	TString nameCutVariationsGammaLegend[ErrorDoubleRatio] = {"#chi^{2} Cut","#psi_{pair} Cut","Findable Cls. TPC","e^{#pm} dEdx","#pi^{#pm} dEdx","q_{T} Cut","min. e^{#pm} #it{p}_{T}","#pi^{0} #alpha_{asym} Cut","#pi^{0} Int. Range","Cocktail #eta Shape","Cocktail #pi^{0} Shape","Cocktail #eta Norm" };

	TFile *fileSystematics = new TFile(nameDataFileErrors);

	Int_t ErrorDoubleRatioToRun = 9;
	Int_t ErrorIncRatioToRun = 9;
	Int_t ErrorPi0ToRun = 9;
	Int_t ErrorGammaSpecToRun = 7;
	
	Double_t* ptBinsDoubleRatio;
	Double_t* ptBinsDoubleRatioErr;
	Double_t* ptBinsIncRatio;
	Double_t* ptBinsIncRatioErr;
	Double_t* ptBinsPi0;
	Double_t* ptBinsPi0Err;
	Double_t* ptBinsGamma;
	Double_t* ptBinsGammaErr;

	const Int_t ConstnBinsDoubleRatio = nBinsDoubleRatio;
	const Int_t ConstnBinsIncRatio = nBinsIncRatio;
	const Int_t ConstnBinsPi0 = nBinsPi0;
	const Int_t ConstnBinsGamma = nBinsGamma;


	TGraphAsymmErrors  **graphPosErrorsDoubleRatio = new TGraphAsymmErrors*[ErrorDoubleRatio];
	TGraphAsymmErrors  **graphNegErrorsDoubleRatio = new TGraphAsymmErrors*[ErrorDoubleRatio];

	TGraphAsymmErrors  **graphPosErrorsDoubleRatioPi0Fit = new TGraphAsymmErrors*[ErrorDoubleRatio];
	TGraphAsymmErrors  **graphNegErrorsDoubleRatioPi0Fit = new TGraphAsymmErrors*[ErrorDoubleRatio];

	Double_t* errorsNegDoubleRatio[ErrorDoubleRatio];
	Double_t errorsNegCorrDoubleRatio[ErrorDoubleRatio][ConstnBinsDoubleRatio];
	Double_t errorsNegSummedDoubleRatio[ConstnBinsDoubleRatio];
	Double_t errorsNegCorrSummedDoubleRatio[ConstnBinsDoubleRatio];
	Double_t errorsNegCorrMatSummedDoubleRatio[ConstnBinsDoubleRatio];

	Double_t* errorsNegErrDoubleRatio[ErrorDoubleRatio];
	Double_t errorsNegErrCorrDoubleRatio[ErrorDoubleRatio][ConstnBinsDoubleRatio];
	Double_t errorsNegErrSummedDoubleRatio[ConstnBinsDoubleRatio];
	Double_t errorsNegErrCorrSummedDoubleRatio[ConstnBinsDoubleRatio];

	Double_t* errorsPosDoubleRatio[ErrorDoubleRatio];
	Double_t errorsPosCorrDoubleRatio[ErrorDoubleRatio][ConstnBinsDoubleRatio];
	Double_t errorsPosSummedDoubleRatio[ConstnBinsDoubleRatio];
	Double_t errorsPosCorrSummedDoubleRatio[ConstnBinsDoubleRatio];
	Double_t errorsPosCorrMatSummedDoubleRatio[ConstnBinsDoubleRatio];

	Double_t* errorsPosErrDoubleRatio[ErrorDoubleRatio];
	Double_t errorsPosErrSummedDoubleRatio[ConstnBinsDoubleRatio];
	Double_t errorsPosErrCorrDoubleRatio[ErrorDoubleRatio][ConstnBinsDoubleRatio];
	Double_t errorsPosErrCorrSummedDoubleRatio[ConstnBinsDoubleRatio];

	Double_t* errorsNegDoubleRatioPi0Fit[ErrorDoubleRatio];
	Double_t errorsNegCorrDoubleRatioPi0Fit[ErrorDoubleRatio][ConstnBinsDoubleRatio];
	Double_t errorsNegSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];
	Double_t errorsNegCorrSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];
	Double_t errorsNegCorrMatSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];

	Double_t* errorsNegErrDoubleRatioPi0Fit[ErrorDoubleRatio];
	Double_t errorsNegErrCorrDoubleRatioPi0Fit[ErrorDoubleRatio][ConstnBinsDoubleRatio];
	Double_t errorsNegErrSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];
	Double_t errorsNegErrCorrSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];

	Double_t* errorsPosDoubleRatioPi0Fit[ErrorDoubleRatio];
	Double_t errorsPosCorrDoubleRatioPi0Fit[ErrorDoubleRatio][ConstnBinsDoubleRatio];
	Double_t errorsPosSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];
	Double_t errorsPosCorrSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];
	Double_t errorsPosCorrMatSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];

	Double_t* errorsPosErrDoubleRatioPi0Fit[ErrorDoubleRatio];
	Double_t errorsPosErrSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];
	Double_t errorsPosErrCorrDoubleRatioPi0Fit[ErrorDoubleRatio][ConstnBinsDoubleRatio];
	Double_t errorsPosErrCorrSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];

	Double_t errorsNegSummedDoubleRatioA[ConstnBinsDoubleRatio];
	Double_t errorsNegCorrSummedDoubleRatioA[ConstnBinsDoubleRatio];
	Double_t errorsNegCorrMatSummedDoubleRatioA[ConstnBinsDoubleRatio];

	Double_t errorsNegErrSummedDoubleRatioA[ConstnBinsDoubleRatio];
	Double_t errorsNegErrCorrSummedDoubleRatioA[ConstnBinsDoubleRatio];

	Double_t errorsPosSummedDoubleRatioA[ConstnBinsDoubleRatio];
	Double_t errorsPosCorrSummedDoubleRatioA[ConstnBinsDoubleRatio];
	Double_t errorsPosCorrMatSummedDoubleRatioA[ConstnBinsDoubleRatio];

	Double_t errorsPosErrSummedDoubleRatioA[ConstnBinsDoubleRatio];
	Double_t errorsPosErrCorrSummedDoubleRatioA[ConstnBinsDoubleRatio];

	Double_t errorsNegSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];
	Double_t errorsNegCorrSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];
	Double_t errorsNegCorrMatSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];

	Double_t errorsNegErrSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];
	Double_t errorsNegErrCorrSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];

	Double_t errorsPosSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];
	Double_t errorsPosCorrSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];
	Double_t errorsPosCorrMatSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];

	Double_t errorsPosErrSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];
	Double_t errorsPosErrCorrSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];

	Double_t errorsNegSummedDoubleRatioB[ConstnBinsDoubleRatio];
	Double_t errorsNegCorrSummedDoubleRatioB[ConstnBinsDoubleRatio];
	Double_t errorsNegCorrMatSummedDoubleRatioB[ConstnBinsDoubleRatio];

	Double_t errorsNegErrSummedDoubleRatioB[ConstnBinsDoubleRatio];
	Double_t errorsNegErrCorrSummedDoubleRatioB[ConstnBinsDoubleRatio];

	Double_t errorsPosSummedDoubleRatioB[ConstnBinsDoubleRatio];
	Double_t errorsPosCorrSummedDoubleRatioB[ConstnBinsDoubleRatio];
	Double_t errorsPosCorrMatSummedDoubleRatioB[ConstnBinsDoubleRatio];

	Double_t errorsPosErrSummedDoubleRatioB[ConstnBinsDoubleRatio];
	Double_t errorsPosErrCorrSummedDoubleRatioB[ConstnBinsDoubleRatio];

	Double_t errorsNegSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];
	Double_t errorsNegCorrSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];
	Double_t errorsNegCorrMatSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];

	Double_t errorsNegErrSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];
	Double_t errorsNegErrCorrSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];

	Double_t errorsPosSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];
	Double_t errorsPosCorrSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];
	Double_t errorsPosCorrMatSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];

	Double_t errorsPosErrSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];
	Double_t errorsPosErrCorrSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];

	Double_t errorsNegSummedDoubleRatioC[ConstnBinsDoubleRatio];
	Double_t errorsNegCorrSummedDoubleRatioC[ConstnBinsDoubleRatio];
	Double_t errorsNegCorrMatSummedDoubleRatioC[ConstnBinsDoubleRatio];

	Double_t errorsNegErrSummedDoubleRatioC[ConstnBinsDoubleRatio];
	Double_t errorsNegErrCorrSummedDoubleRatioC[ConstnBinsDoubleRatio];

	Double_t errorsPosSummedDoubleRatioC[ConstnBinsDoubleRatio];
	Double_t errorsPosCorrSummedDoubleRatioC[ConstnBinsDoubleRatio];
	Double_t errorsPosCorrMatSummedDoubleRatioC[ConstnBinsDoubleRatio];

	Double_t errorsPosErrSummedDoubleRatioC[ConstnBinsDoubleRatio];
	Double_t errorsPosErrCorrSummedDoubleRatioC[ConstnBinsDoubleRatio];

	Double_t errorsNegSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];
	Double_t errorsNegCorrSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];
	Double_t errorsNegCorrMatSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];

	Double_t errorsNegErrSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];
	Double_t errorsNegErrCorrSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];

	Double_t errorsPosSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];
	Double_t errorsPosCorrSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];
	Double_t errorsPosCorrMatSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];

	Double_t errorsPosErrSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];
	Double_t errorsPosErrCorrSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];

	Double_t errorsMeanDoubleRatio[ErrorDoubleRatio][ConstnBinsDoubleRatio];
	Double_t errorsMeanCorrDoubleRatio[ErrorDoubleRatio][ConstnBinsDoubleRatio];
	Double_t errorsMeanSummedDoubleRatio[ConstnBinsDoubleRatio];
	Double_t errorsMeanCorrSummedDoubleRatio[ConstnBinsDoubleRatio];
	Double_t errorsMeanCorrMatSummedDoubleRatio[ConstnBinsDoubleRatio];

	Double_t errorsMeanErrDoubleRatio[ErrorDoubleRatio][ConstnBinsDoubleRatio];
	Double_t errorsMeanErrCorrDoubleRatio[ErrorDoubleRatio][ConstnBinsDoubleRatio];
	Double_t errorsMeanErrSummedDoubleRatio[ConstnBinsDoubleRatio];
	Double_t errorsMeanErrCorrSummedDoubleRatio[ConstnBinsDoubleRatio];
	Double_t errorsMeanErrCorrMatSummedDoubleRatio[ConstnBinsDoubleRatio];

	Double_t errorsMeanDoubleRatioPi0Fit[ErrorDoubleRatio][ConstnBinsDoubleRatio];
	Double_t errorsMeanCorrDoubleRatioPi0Fit[ErrorDoubleRatio][ConstnBinsDoubleRatio];
	Double_t errorsMeanSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];
	Double_t errorsMeanCorrSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];
	Double_t errorsMeanCorrMatSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];

	Double_t errorsMeanErrDoubleRatioPi0Fit[ErrorDoubleRatio][ConstnBinsDoubleRatio];
	Double_t errorsMeanErrCorrDoubleRatioPi0Fit[ErrorDoubleRatio][ConstnBinsDoubleRatio];
	Double_t errorsMeanErrSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];
	Double_t errorsMeanErrCorrSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];
	Double_t errorsMeanErrCorrMatSummedDoubleRatioPi0Fit[ConstnBinsDoubleRatio];

	Double_t errorsMeanSummedDoubleRatioA[ConstnBinsDoubleRatio];
	Double_t errorsMeanCorrSummedDoubleRatioA[ConstnBinsDoubleRatio];
	Double_t errorsMeanCorrMatSummedDoubleRatioA[ConstnBinsDoubleRatio];

	Double_t errorsMeanErrSummedDoubleRatioA[ConstnBinsDoubleRatio];
	Double_t errorsMeanErrCorrSummedDoubleRatioA[ConstnBinsDoubleRatio];
	Double_t errorsMeanErrCorrMatSummedDoubleRatioA[ConstnBinsDoubleRatio];

	Double_t errorsMeanSummedDoubleRatioB[ConstnBinsDoubleRatio];
	Double_t errorsMeanCorrSummedDoubleRatioB[ConstnBinsDoubleRatio];
	Double_t errorsMeanCorrMatSummedDoubleRatioB[ConstnBinsDoubleRatio];

	Double_t errorsMeanErrSummedDoubleRatioB[ConstnBinsDoubleRatio];
	Double_t errorsMeanErrCorrSummedDoubleRatioB[ConstnBinsDoubleRatio];
	Double_t errorsMeanErrCorrMatSummedDoubleRatioB[ConstnBinsDoubleRatio];

	Double_t errorsMeanSummedDoubleRatioC[ConstnBinsDoubleRatio];
	Double_t errorsMeanCorrSummedDoubleRatioC[ConstnBinsDoubleRatio];
	Double_t errorsMeanCorrMatSummedDoubleRatioC[ConstnBinsDoubleRatio];

	Double_t errorsMeanErrSummedDoubleRatioC[ConstnBinsDoubleRatio];
	Double_t errorsMeanErrCorrSummedDoubleRatioC[ConstnBinsDoubleRatio];
	Double_t errorsMeanErrCorrMatSummedDoubleRatioC[ConstnBinsDoubleRatio];


	Double_t errorsMeanSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];
	Double_t errorsMeanCorrSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];
	Double_t errorsMeanCorrMatSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];

	Double_t errorsMeanErrSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];
	Double_t errorsMeanErrCorrSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];
	Double_t errorsMeanErrCorrMatSummedDoubleRatioPi0FitA[ConstnBinsDoubleRatio];

	Double_t errorsMeanSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];
	Double_t errorsMeanCorrSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];
	Double_t errorsMeanCorrMatSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];

	Double_t errorsMeanErrSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];
	Double_t errorsMeanErrCorrSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];
	Double_t errorsMeanErrCorrMatSummedDoubleRatioPi0FitB[ConstnBinsDoubleRatio];

	Double_t errorsMeanSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];
	Double_t errorsMeanCorrSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];
	Double_t errorsMeanCorrMatSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];

	Double_t errorsMeanErrSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];
	Double_t errorsMeanErrCorrSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];
	Double_t errorsMeanErrCorrMatSummedDoubleRatioPi0FitC[ConstnBinsDoubleRatio];

	Double_t errorsMeanCorrMatSummedDoubleRatioPi0FitWithoutMat[ConstnBinsDoubleRatio];
	Double_t errorsMaterialBudgetDoubleRatioPi0Fit[ConstnBinsDoubleRatio];
	Double_t errorsMeanErrCorrMatSummedDoubleRatioPi0FitWithoutMat[ConstnBinsDoubleRatio];
	Double_t errorsErrMaterialBudgetDoubleRatioPi0Fit[ConstnBinsDoubleRatio];

	for (Int_t l = 0; l < ConstnBinsDoubleRatio; l++){
		errorsPosSummedDoubleRatio[l] = 0.;
		errorsNegSummedDoubleRatio[l] = 0.;
		errorsMeanSummedDoubleRatio[l] = 0.;
		errorsPosCorrSummedDoubleRatio[l] = 0.;
		errorsNegCorrSummedDoubleRatio[l] = 0.;
		errorsMeanCorrSummedDoubleRatio[l] = 0.;

		errorsPosSummedDoubleRatioPi0Fit[l] = 0.;
		errorsNegSummedDoubleRatioPi0Fit[l] = 0.;
		errorsMeanSummedDoubleRatioPi0Fit[l] = 0.;
		errorsPosCorrSummedDoubleRatioPi0Fit[l] = 0.;
		errorsNegCorrSummedDoubleRatioPi0Fit[l] = 0.;
		errorsMeanCorrSummedDoubleRatioPi0Fit[l] = 0.;

		errorsPosSummedDoubleRatioA[l] = 0.;
		errorsNegSummedDoubleRatioA[l] = 0.;
		errorsMeanSummedDoubleRatioA[l] = 0.;
		errorsPosCorrSummedDoubleRatioA[l] = 0.;
		errorsNegCorrSummedDoubleRatioA[l] = 0.;
		errorsMeanCorrSummedDoubleRatioA[l] = 0.;

		errorsPosSummedDoubleRatioPi0FitA[l] = 0.;
		errorsNegSummedDoubleRatioPi0FitA[l] = 0.;
		errorsMeanSummedDoubleRatioPi0FitA[l] = 0.;
		errorsPosCorrSummedDoubleRatioPi0FitA[l] = 0.;
		errorsNegCorrSummedDoubleRatioPi0FitA[l] = 0.;
		errorsMeanCorrSummedDoubleRatioPi0FitA[l] = 0.;

		errorsPosSummedDoubleRatioB[l] = 0.;
		errorsNegSummedDoubleRatioB[l] = 0.;
		errorsMeanSummedDoubleRatioB[l] = 0.;
		errorsPosCorrSummedDoubleRatioB[l] = 0.;
		errorsNegCorrSummedDoubleRatioB[l] = 0.;
		errorsMeanCorrSummedDoubleRatioB[l] = 0.;

		errorsPosSummedDoubleRatioPi0FitB[l] = 0.;
		errorsNegSummedDoubleRatioPi0FitB[l] = 0.;
		errorsMeanSummedDoubleRatioPi0FitB[l] = 0.;
		errorsPosCorrSummedDoubleRatioPi0FitB[l] = 0.;
		errorsNegCorrSummedDoubleRatioPi0FitB[l] = 0.;
		errorsMeanCorrSummedDoubleRatioPi0FitB[l] = 0.;

		errorsPosSummedDoubleRatioC[l] = 0.;
		errorsNegSummedDoubleRatioC[l] = 0.;
		errorsMeanSummedDoubleRatioC[l] = 0.;
		errorsPosCorrSummedDoubleRatioC[l] = 0.;
		errorsNegCorrSummedDoubleRatioC[l] = 0.;
		errorsMeanCorrSummedDoubleRatioC[l] = 0.;

		errorsPosSummedDoubleRatioPi0FitC[l] = 0.;
		errorsNegSummedDoubleRatioPi0FitC[l] = 0.;
		errorsMeanSummedDoubleRatioPi0FitC[l] = 0.;
		errorsPosCorrSummedDoubleRatioPi0FitC[l] = 0.;
		errorsNegCorrSummedDoubleRatioPi0FitC[l] = 0.;
		errorsMeanCorrSummedDoubleRatioPi0FitC[l] = 0.;
	}


	
	TGraphAsymmErrors  **graphPosErrorsIncRatio = new TGraphAsymmErrors*[ErrorIncRatio];
	TGraphAsymmErrors  **graphNegErrorsIncRatio = new TGraphAsymmErrors*[ErrorIncRatio];

	TGraphAsymmErrors  **graphPosErrorsIncRatioPi0Fit = new TGraphAsymmErrors*[ErrorIncRatio];
	TGraphAsymmErrors  **graphNegErrorsIncRatioPi0Fit = new TGraphAsymmErrors*[ErrorIncRatio];

	Double_t* errorsNegIncRatio[ErrorIncRatio];
	Double_t errorsNegCorrIncRatio[ErrorIncRatio][ConstnBinsIncRatio];
	Double_t errorsNegSummedIncRatio[ConstnBinsIncRatio];
	Double_t errorsNegCorrSummedIncRatio[ConstnBinsIncRatio];
	Double_t errorsNegCorrMatSummedIncRatio[ConstnBinsIncRatio];

	Double_t* errorsNegErrIncRatio[ErrorIncRatio];
	Double_t errorsNegErrCorrIncRatio[ErrorIncRatio][ConstnBinsIncRatio];
	Double_t errorsNegErrSummedIncRatio[ConstnBinsIncRatio];
	Double_t errorsNegErrCorrSummedIncRatio[ConstnBinsIncRatio];

	Double_t* errorsPosIncRatio[ErrorIncRatio];
	Double_t errorsPosCorrIncRatio[ErrorIncRatio][ConstnBinsIncRatio];
	Double_t errorsPosSummedIncRatio[ConstnBinsIncRatio];
	Double_t errorsPosCorrSummedIncRatio[ConstnBinsIncRatio];
	Double_t errorsPosCorrMatSummedIncRatio[ConstnBinsIncRatio];

	Double_t* errorsPosErrIncRatio[ErrorIncRatio];
	Double_t errorsPosErrSummedIncRatio[ConstnBinsIncRatio];
	Double_t errorsPosErrCorrIncRatio[ErrorIncRatio][ConstnBinsIncRatio];
	Double_t errorsPosErrCorrSummedIncRatio[ConstnBinsIncRatio];

	Double_t* errorsNegIncRatioPi0Fit[ErrorIncRatio];
	Double_t errorsNegCorrIncRatioPi0Fit[ErrorIncRatio][ConstnBinsIncRatio];
	Double_t errorsNegSummedIncRatioPi0Fit[ConstnBinsIncRatio];
	Double_t errorsNegCorrSummedIncRatioPi0Fit[ConstnBinsIncRatio];
	Double_t errorsNegCorrMatSummedIncRatioPi0Fit[ConstnBinsIncRatio];

	Double_t* errorsNegErrIncRatioPi0Fit[ErrorIncRatio];
	Double_t errorsNegErrCorrIncRatioPi0Fit[ErrorIncRatio][ConstnBinsIncRatio];
	Double_t errorsNegErrSummedIncRatioPi0Fit[ConstnBinsIncRatio];
	Double_t errorsNegErrCorrSummedIncRatioPi0Fit[ConstnBinsIncRatio];

	Double_t* errorsPosIncRatioPi0Fit[ErrorIncRatio];
	Double_t errorsPosCorrIncRatioPi0Fit[ErrorIncRatio][ConstnBinsIncRatio];
	Double_t errorsPosSummedIncRatioPi0Fit[ConstnBinsIncRatio];
	Double_t errorsPosCorrSummedIncRatioPi0Fit[ConstnBinsIncRatio];
	Double_t errorsPosCorrMatSummedIncRatioPi0Fit[ConstnBinsIncRatio];

	Double_t* errorsPosErrIncRatioPi0Fit[ErrorIncRatio];
	Double_t errorsPosErrSummedIncRatioPi0Fit[ConstnBinsIncRatio];
	Double_t errorsPosErrCorrIncRatioPi0Fit[ErrorIncRatio][ConstnBinsIncRatio];
	Double_t errorsPosErrCorrSummedIncRatioPi0Fit[ConstnBinsIncRatio];

	Double_t errorsNegSummedIncRatioA[ConstnBinsIncRatio];
	Double_t errorsNegCorrSummedIncRatioA[ConstnBinsIncRatio];
	Double_t errorsNegCorrMatSummedIncRatioA[ConstnBinsIncRatio];

	Double_t errorsNegErrSummedIncRatioA[ConstnBinsIncRatio];
	Double_t errorsNegErrCorrSummedIncRatioA[ConstnBinsIncRatio];

	Double_t errorsPosSummedIncRatioA[ConstnBinsIncRatio];
	Double_t errorsPosCorrSummedIncRatioA[ConstnBinsIncRatio];
	Double_t errorsPosCorrMatSummedIncRatioA[ConstnBinsIncRatio];

	Double_t errorsPosErrSummedIncRatioA[ConstnBinsIncRatio];
	Double_t errorsPosErrCorrSummedIncRatioA[ConstnBinsIncRatio];

	Double_t errorsNegSummedIncRatioPi0FitA[ConstnBinsIncRatio];
	Double_t errorsNegCorrSummedIncRatioPi0FitA[ConstnBinsIncRatio];
	Double_t errorsNegCorrMatSummedIncRatioPi0FitA[ConstnBinsIncRatio];

	Double_t errorsNegErrSummedIncRatioPi0FitA[ConstnBinsIncRatio];
	Double_t errorsNegErrCorrSummedIncRatioPi0FitA[ConstnBinsIncRatio];

	Double_t errorsPosSummedIncRatioPi0FitA[ConstnBinsIncRatio];
	Double_t errorsPosCorrSummedIncRatioPi0FitA[ConstnBinsIncRatio];
	Double_t errorsPosCorrMatSummedIncRatioPi0FitA[ConstnBinsIncRatio];

	Double_t errorsPosErrSummedIncRatioPi0FitA[ConstnBinsIncRatio];
	Double_t errorsPosErrCorrSummedIncRatioPi0FitA[ConstnBinsIncRatio];

	Double_t errorsNegSummedIncRatioB[ConstnBinsIncRatio];
	Double_t errorsNegCorrSummedIncRatioB[ConstnBinsIncRatio];
	Double_t errorsNegCorrMatSummedIncRatioB[ConstnBinsIncRatio];

	Double_t errorsNegErrSummedIncRatioB[ConstnBinsIncRatio];
	Double_t errorsNegErrCorrSummedIncRatioB[ConstnBinsIncRatio];

	Double_t errorsPosSummedIncRatioB[ConstnBinsIncRatio];
	Double_t errorsPosCorrSummedIncRatioB[ConstnBinsIncRatio];
	Double_t errorsPosCorrMatSummedIncRatioB[ConstnBinsIncRatio];

	Double_t errorsPosErrSummedIncRatioB[ConstnBinsIncRatio];
	Double_t errorsPosErrCorrSummedIncRatioB[ConstnBinsIncRatio];

	Double_t errorsNegSummedIncRatioPi0FitB[ConstnBinsIncRatio];
	Double_t errorsNegCorrSummedIncRatioPi0FitB[ConstnBinsIncRatio];
	Double_t errorsNegCorrMatSummedIncRatioPi0FitB[ConstnBinsIncRatio];

	Double_t errorsNegErrSummedIncRatioPi0FitB[ConstnBinsIncRatio];
	Double_t errorsNegErrCorrSummedIncRatioPi0FitB[ConstnBinsIncRatio];

	Double_t errorsPosSummedIncRatioPi0FitB[ConstnBinsIncRatio];
	Double_t errorsPosCorrSummedIncRatioPi0FitB[ConstnBinsIncRatio];
	Double_t errorsPosCorrMatSummedIncRatioPi0FitB[ConstnBinsIncRatio];

	Double_t errorsPosErrSummedIncRatioPi0FitB[ConstnBinsIncRatio];
	Double_t errorsPosErrCorrSummedIncRatioPi0FitB[ConstnBinsIncRatio];

	Double_t errorsNegSummedIncRatioC[ConstnBinsIncRatio];
	Double_t errorsNegCorrSummedIncRatioC[ConstnBinsIncRatio];
	Double_t errorsNegCorrMatSummedIncRatioC[ConstnBinsIncRatio];

	Double_t errorsNegErrSummedIncRatioC[ConstnBinsIncRatio];
	Double_t errorsNegErrCorrSummedIncRatioC[ConstnBinsIncRatio];

	Double_t errorsPosSummedIncRatioC[ConstnBinsIncRatio];
	Double_t errorsPosCorrSummedIncRatioC[ConstnBinsIncRatio];
	Double_t errorsPosCorrMatSummedIncRatioC[ConstnBinsIncRatio];

	Double_t errorsPosErrSummedIncRatioC[ConstnBinsIncRatio];
	Double_t errorsPosErrCorrSummedIncRatioC[ConstnBinsIncRatio];

	Double_t errorsNegSummedIncRatioPi0FitC[ConstnBinsIncRatio];
	Double_t errorsNegCorrSummedIncRatioPi0FitC[ConstnBinsIncRatio];
	Double_t errorsNegCorrMatSummedIncRatioPi0FitC[ConstnBinsIncRatio];

	Double_t errorsNegErrSummedIncRatioPi0FitC[ConstnBinsIncRatio];
	Double_t errorsNegErrCorrSummedIncRatioPi0FitC[ConstnBinsIncRatio];

	Double_t errorsPosSummedIncRatioPi0FitC[ConstnBinsIncRatio];
	Double_t errorsPosCorrSummedIncRatioPi0FitC[ConstnBinsIncRatio];
	Double_t errorsPosCorrMatSummedIncRatioPi0FitC[ConstnBinsIncRatio];

	Double_t errorsPosErrSummedIncRatioPi0FitC[ConstnBinsIncRatio];
	Double_t errorsPosErrCorrSummedIncRatioPi0FitC[ConstnBinsIncRatio];

	Double_t errorsMeanIncRatio[ErrorIncRatio][ConstnBinsIncRatio];
	Double_t errorsMeanCorrIncRatio[ErrorIncRatio][ConstnBinsIncRatio];
	Double_t errorsMeanSummedIncRatio[ConstnBinsIncRatio];
	Double_t errorsMeanCorrSummedIncRatio[ConstnBinsIncRatio];
	Double_t errorsMeanCorrMatSummedIncRatio[ConstnBinsIncRatio];

	Double_t errorsMeanErrIncRatio[ErrorIncRatio][ConstnBinsIncRatio];
	Double_t errorsMeanErrCorrIncRatio[ErrorIncRatio][ConstnBinsIncRatio];
	Double_t errorsMeanErrSummedIncRatio[ConstnBinsIncRatio];
	Double_t errorsMeanErrCorrSummedIncRatio[ConstnBinsIncRatio];
	Double_t errorsMeanErrCorrMatSummedIncRatio[ConstnBinsIncRatio];

	Double_t errorsMeanIncRatioPi0Fit[ErrorIncRatio][ConstnBinsIncRatio];
	Double_t errorsMeanCorrIncRatioPi0Fit[ErrorIncRatio][ConstnBinsIncRatio];
	Double_t errorsMeanSummedIncRatioPi0Fit[ConstnBinsIncRatio];
	Double_t errorsMeanCorrSummedIncRatioPi0Fit[ConstnBinsIncRatio];
	Double_t errorsMeanCorrMatSummedIncRatioPi0Fit[ConstnBinsIncRatio];

	Double_t errorsMeanErrIncRatioPi0Fit[ErrorIncRatio][ConstnBinsIncRatio];
	Double_t errorsMeanErrCorrIncRatioPi0Fit[ErrorIncRatio][ConstnBinsIncRatio];
	Double_t errorsMeanErrSummedIncRatioPi0Fit[ConstnBinsIncRatio];
	Double_t errorsMeanErrCorrSummedIncRatioPi0Fit[ConstnBinsIncRatio];
	Double_t errorsMeanErrCorrMatSummedIncRatioPi0Fit[ConstnBinsIncRatio];

	Double_t errorsMeanSummedIncRatioA[ConstnBinsIncRatio];
	Double_t errorsMeanCorrSummedIncRatioA[ConstnBinsIncRatio];
	Double_t errorsMeanCorrMatSummedIncRatioA[ConstnBinsIncRatio];

	Double_t errorsMeanErrSummedIncRatioA[ConstnBinsIncRatio];
	Double_t errorsMeanErrCorrSummedIncRatioA[ConstnBinsIncRatio];
	Double_t errorsMeanErrCorrMatSummedIncRatioA[ConstnBinsIncRatio];

	Double_t errorsMeanSummedIncRatioB[ConstnBinsIncRatio];
	Double_t errorsMeanCorrSummedIncRatioB[ConstnBinsIncRatio];
	Double_t errorsMeanCorrMatSummedIncRatioB[ConstnBinsIncRatio];

	Double_t errorsMeanErrSummedIncRatioB[ConstnBinsIncRatio];
	Double_t errorsMeanErrCorrSummedIncRatioB[ConstnBinsIncRatio];
	Double_t errorsMeanErrCorrMatSummedIncRatioB[ConstnBinsIncRatio];

	Double_t errorsMeanSummedIncRatioC[ConstnBinsIncRatio];
	Double_t errorsMeanCorrSummedIncRatioC[ConstnBinsIncRatio];
	Double_t errorsMeanCorrMatSummedIncRatioC[ConstnBinsIncRatio];

	Double_t errorsMeanErrSummedIncRatioC[ConstnBinsIncRatio];
	Double_t errorsMeanErrCorrSummedIncRatioC[ConstnBinsIncRatio];
	Double_t errorsMeanErrCorrMatSummedIncRatioC[ConstnBinsIncRatio];

	Double_t errorsMeanSummedIncRatioPi0FitA[ConstnBinsIncRatio];
	Double_t errorsMeanCorrSummedIncRatioPi0FitA[ConstnBinsIncRatio];
	Double_t errorsMeanCorrMatSummedIncRatioPi0FitA[ConstnBinsIncRatio];

	Double_t errorsMeanErrSummedIncRatioPi0FitA[ConstnBinsIncRatio];
	Double_t errorsMeanErrCorrSummedIncRatioPi0FitA[ConstnBinsIncRatio];
	Double_t errorsMeanErrCorrMatSummedIncRatioPi0FitA[ConstnBinsIncRatio];

	Double_t errorsMeanSummedIncRatioPi0FitB[ConstnBinsIncRatio];
	Double_t errorsMeanCorrSummedIncRatioPi0FitB[ConstnBinsIncRatio];
	Double_t errorsMeanCorrMatSummedIncRatioPi0FitB[ConstnBinsIncRatio];

	Double_t errorsMeanErrSummedIncRatioPi0FitB[ConstnBinsIncRatio];
	Double_t errorsMeanErrCorrSummedIncRatioPi0FitB[ConstnBinsIncRatio];
	Double_t errorsMeanErrCorrMatSummedIncRatioPi0FitB[ConstnBinsIncRatio];

	Double_t errorsMeanSummedIncRatioPi0FitC[ConstnBinsIncRatio];
	Double_t errorsMeanCorrSummedIncRatioPi0FitC[ConstnBinsIncRatio];
	Double_t errorsMeanCorrMatSummedIncRatioPi0FitC[ConstnBinsIncRatio];

	Double_t errorsMeanErrSummedIncRatioPi0FitC[ConstnBinsIncRatio];
	Double_t errorsMeanErrCorrSummedIncRatioPi0FitC[ConstnBinsIncRatio];
	Double_t errorsMeanErrCorrMatSummedIncRatioPi0FitC[ConstnBinsIncRatio];	


	for (Int_t l = 0; l < ConstnBinsIncRatio; l++){
		errorsPosSummedIncRatio[l] = 0.;
		errorsNegSummedIncRatio[l] = 0.;
		errorsMeanSummedIncRatio[l] = 0.;
		errorsPosCorrSummedIncRatio[l] = 0.;
		errorsNegCorrSummedIncRatio[l] = 0.;
		errorsMeanCorrSummedIncRatio[l] = 0.;

		errorsPosSummedIncRatioPi0Fit[l] = 0.;
		errorsNegSummedIncRatioPi0Fit[l] = 0.;
		errorsMeanSummedIncRatioPi0Fit[l] = 0.;
		errorsPosCorrSummedIncRatioPi0Fit[l] = 0.;
		errorsNegCorrSummedIncRatioPi0Fit[l] = 0.;
		errorsMeanCorrSummedIncRatioPi0Fit[l] = 0.;

		errorsPosSummedIncRatioA[l] = 0.;
		errorsNegSummedIncRatioA[l] = 0.;
		errorsMeanSummedIncRatioA[l] = 0.;
		errorsPosCorrSummedIncRatioA[l] = 0.;
		errorsNegCorrSummedIncRatioA[l] = 0.;
		errorsMeanCorrSummedIncRatioA[l] = 0.;

		errorsPosSummedIncRatioPi0FitA[l] = 0.;
		errorsNegSummedIncRatioPi0FitA[l] = 0.;
		errorsMeanSummedIncRatioPi0FitA[l] = 0.;
		errorsPosCorrSummedIncRatioPi0FitA[l] = 0.;
		errorsNegCorrSummedIncRatioPi0FitA[l] = 0.;
		errorsMeanCorrSummedIncRatioPi0FitA[l] = 0.;

		errorsPosSummedIncRatioB[l] = 0.;
		errorsNegSummedIncRatioB[l] = 0.;
		errorsMeanSummedIncRatioB[l] = 0.;
		errorsPosCorrSummedIncRatioB[l] = 0.;
		errorsNegCorrSummedIncRatioB[l] = 0.;
		errorsMeanCorrSummedIncRatioB[l] = 0.;

		errorsPosSummedIncRatioPi0FitB[l] = 0.;
		errorsNegSummedIncRatioPi0FitB[l] = 0.;
		errorsMeanSummedIncRatioPi0FitB[l] = 0.;
		errorsPosCorrSummedIncRatioPi0FitB[l] = 0.;
		errorsNegCorrSummedIncRatioPi0FitB[l] = 0.;
		errorsMeanCorrSummedIncRatioPi0FitB[l] = 0.;

		errorsPosSummedIncRatioC[l] = 0.;
		errorsNegSummedIncRatioC[l] = 0.;
		errorsMeanSummedIncRatioC[l] = 0.;
		errorsPosCorrSummedIncRatioC[l] = 0.;
		errorsNegCorrSummedIncRatioC[l] = 0.;
		errorsMeanCorrSummedIncRatioC[l] = 0.;

		errorsPosSummedIncRatioPi0FitC[l] = 0.;
		errorsNegSummedIncRatioPi0FitC[l] = 0.;
		errorsMeanSummedIncRatioPi0FitC[l] = 0.;
		errorsPosCorrSummedIncRatioPi0FitC[l] = 0.;
		errorsNegCorrSummedIncRatioPi0FitC[l] = 0.;
		errorsMeanCorrSummedIncRatioPi0FitC[l] = 0.;

	}

	
	TGraphAsymmErrors  **graphPosErrorsPi0 = new TGraphAsymmErrors*[ErrorPi0];
	TGraphAsymmErrors  **graphNegErrorsPi0 = new TGraphAsymmErrors*[ErrorPi0];

	TGraphAsymmErrors  **graphPosErrorsPi0Fit = new TGraphAsymmErrors*[ErrorPi0];
	TGraphAsymmErrors  **graphNegErrorsPi0Fit = new TGraphAsymmErrors*[ErrorPi0];

	Double_t* errorsNegPi0[ErrorPi0];
	Double_t errorsNegCorrPi0[ErrorPi0][ConstnBinsPi0];
	Double_t errorsNegSummedPi0[ConstnBinsPi0];
	Double_t errorsNegCorrSummedPi0[ConstnBinsPi0];
	Double_t errorsNegCorrMatSummedPi0[ConstnBinsPi0];

	Double_t* errorsNegErrPi0[ErrorPi0];
	Double_t errorsNegErrCorrPi0[ErrorPi0][ConstnBinsPi0];
	Double_t errorsNegErrSummedPi0[ConstnBinsPi0];
	Double_t errorsNegErrCorrSummedPi0[ConstnBinsPi0];

	Double_t* errorsPosPi0[ErrorPi0];
	Double_t errorsPosCorrPi0[ErrorPi0][ConstnBinsPi0];
	Double_t errorsPosSummedPi0[ConstnBinsPi0];
	Double_t errorsPosCorrSummedPi0[ConstnBinsPi0];
	Double_t errorsPosCorrMatSummedPi0[ConstnBinsPi0];

	Double_t* errorsPosErrPi0[ErrorPi0];
	Double_t errorsPosErrSummedPi0[ConstnBinsPi0];
	Double_t errorsPosErrCorrPi0[ErrorPi0][ConstnBinsPi0];
	Double_t errorsPosErrCorrSummedPi0[ConstnBinsPi0];

	Double_t* errorsNegPi0Fit[ErrorPi0];
	Double_t errorsNegCorrPi0Fit[ErrorPi0][ConstnBinsPi0];
	Double_t errorsNegSummedPi0Fit[ConstnBinsPi0];
	Double_t errorsNegCorrSummedPi0Fit[ConstnBinsPi0];
	Double_t errorsNegCorrMatSummedPi0Fit[ConstnBinsPi0];

	Double_t* errorsNegErrPi0Fit[ErrorPi0];
	Double_t errorsNegErrCorrPi0Fit[ErrorPi0][ConstnBinsPi0];
	Double_t errorsNegErrSummedPi0Fit[ConstnBinsPi0];
	Double_t errorsNegErrCorrSummedPi0Fit[ConstnBinsPi0];

	Double_t* errorsPosPi0Fit[ErrorPi0];
	Double_t errorsPosCorrPi0Fit[ErrorPi0][ConstnBinsPi0];
	Double_t errorsPosSummedPi0Fit[ConstnBinsPi0];
	Double_t errorsPosCorrSummedPi0Fit[ConstnBinsPi0];
	Double_t errorsPosCorrMatSummedPi0Fit[ConstnBinsPi0];

	Double_t* errorsPosErrPi0Fit[ErrorPi0];
	Double_t errorsPosErrSummedPi0Fit[ConstnBinsPi0];
	Double_t errorsPosErrCorrPi0Fit[ErrorPi0][ConstnBinsPi0];
	Double_t errorsPosErrCorrSummedPi0Fit[ConstnBinsPi0];

	Double_t errorsMeanPi0[ErrorPi0][ConstnBinsPi0];
	Double_t errorsMeanCorrPi0[ErrorPi0][ConstnBinsPi0];
	Double_t errorsMeanSummedPi0[ConstnBinsPi0];
	Double_t errorsMeanCorrSummedPi0[ConstnBinsPi0];
	Double_t errorsMeanCorrMatSummedPi0[ConstnBinsPi0];
	Double_t errorsMeanCorrMatSummedPi0WithoutMaterial[ConstnBinsPi0];

	Double_t errorsMeanErrPi0[ErrorPi0][ConstnBinsPi0];
	Double_t errorsMeanErrCorrPi0[ErrorPi0][ConstnBinsPi0];
	Double_t errorsMeanErrSummedPi0[ConstnBinsPi0];
	Double_t errorsMeanErrCorrSummedPi0[ConstnBinsPi0];
	Double_t errorsMeanErrCorrMatSummedPi0[ConstnBinsPi0];

	Double_t errorsMeanPi0Fit[ErrorPi0][ConstnBinsPi0];
	Double_t errorsMeanCorrPi0Fit[ErrorPi0][ConstnBinsPi0];
	Double_t errorsMeanSummedPi0Fit[ConstnBinsPi0];
	Double_t errorsMeanCorrSummedPi0Fit[ConstnBinsPi0];
	Double_t errorsMeanCorrMatSummedPi0Fit[ConstnBinsPi0];
	Double_t errorsMeanCorrMatSummedPi0FitWithoutMaterial[ConstnBinsPi0];

	Double_t errorsMeanErrPi0Fit[ErrorPi0][ConstnBinsPi0];
	Double_t errorsMeanErrCorrPi0Fit[ErrorPi0][ConstnBinsPi0];
	Double_t errorsMeanErrSummedPi0Fit[ConstnBinsPi0];
	Double_t errorsMeanErrCorrSummedPi0Fit[ConstnBinsPi0];
	Double_t errorsMeanErrCorrMatSummedPi0Fit[ConstnBinsPi0];

	
	for (Int_t l = 0; l < ConstnBinsPi0; l++){
		errorsPosSummedPi0[l] = 0.;
		errorsNegSummedPi0[l] = 0.;
		errorsMeanSummedPi0[l] = 0.;
		errorsPosCorrSummedPi0[l] = 0.;
		errorsNegCorrSummedPi0[l] = 0.;
		errorsMeanCorrSummedPi0[l] = 0.;

		errorsPosSummedPi0Fit[l] = 0.;
		errorsNegSummedPi0Fit[l] = 0.;
		errorsMeanSummedPi0Fit[l] = 0.;
		errorsPosCorrSummedPi0Fit[l] = 0.;
		errorsNegCorrSummedPi0Fit[l] = 0.;
		errorsMeanCorrSummedPi0Fit[l] = 0.;
	}

	
	TGraphAsymmErrors  **graphPosErrorsGamma = new TGraphAsymmErrors*[ErrorGammaSpec];
	TGraphAsymmErrors  **graphNegErrorsGamma = new TGraphAsymmErrors*[ErrorGammaSpec];

	Double_t* errorsNegGamma[ErrorGammaSpec];
	Double_t errorsNegCorrGamma[ErrorGammaSpec][ConstnBinsGamma];
	Double_t errorsNegSummedGamma[ConstnBinsGamma];
	Double_t errorsNegCorrSummedGamma[ConstnBinsGamma];
	Double_t errorsNegCorrMatSummedGamma[ConstnBinsGamma];

	Double_t* errorsNegErrGamma[ErrorGammaSpec];
	Double_t errorsNegErrCorrGamma[ErrorGammaSpec][ConstnBinsGamma];
	Double_t errorsNegErrSummedGamma[ConstnBinsGamma];
	Double_t errorsNegErrCorrSummedGamma[ConstnBinsGamma];

	Double_t* errorsPosGamma[ErrorGammaSpec];
	Double_t errorsPosCorrGamma[ErrorGammaSpec][ConstnBinsGamma];
	Double_t errorsPosSummedGamma[ConstnBinsGamma];
	Double_t errorsPosCorrSummedGamma[ConstnBinsGamma];
	Double_t errorsPosCorrMatSummedGamma[ConstnBinsGamma];

	Double_t* errorsPosErrGamma[ErrorGammaSpec];
	Double_t errorsPosErrSummedGamma[ConstnBinsGamma];
	Double_t errorsPosErrCorrGamma[ErrorGammaSpec][ConstnBinsGamma];
	Double_t errorsPosErrCorrSummedGamma[ConstnBinsGamma];

	Double_t errorsNegSummedGammaA[ConstnBinsGamma];
	Double_t errorsNegCorrSummedGammaA[ConstnBinsGamma];
	Double_t errorsNegCorrMatSummedGammaA[ConstnBinsGamma];

	Double_t errorsNegErrSummedGammaA[ConstnBinsGamma];
	Double_t errorsNegErrCorrSummedGammaA[ConstnBinsGamma];

	Double_t errorsPosSummedGammaA[ConstnBinsGamma];
	Double_t errorsPosCorrSummedGammaA[ConstnBinsGamma];
	Double_t errorsPosCorrMatSummedGammaA[ConstnBinsGamma];

	Double_t errorsPosErrSummedGammaA[ConstnBinsGamma];
	Double_t errorsPosErrCorrSummedGammaA[ConstnBinsGamma];

	Double_t errorsNegSummedGammaB[ConstnBinsGamma];
	Double_t errorsNegCorrSummedGammaB[ConstnBinsGamma];
	Double_t errorsNegCorrMatSummedGammaB[ConstnBinsGamma];

	Double_t errorsNegErrSummedGammaB[ConstnBinsGamma];
	Double_t errorsNegErrCorrSummedGammaB[ConstnBinsGamma];

	Double_t errorsPosSummedGammaB[ConstnBinsGamma];
	Double_t errorsPosCorrSummedGammaB[ConstnBinsGamma];
	Double_t errorsPosCorrMatSummedGammaB[ConstnBinsGamma];

	Double_t errorsPosErrSummedGammaB[ConstnBinsGamma];
	Double_t errorsPosErrCorrSummedGammaB[ConstnBinsGamma];

	Double_t errorsNegSummedGammaC[ConstnBinsGamma];
	Double_t errorsNegCorrSummedGammaC[ConstnBinsGamma];
	Double_t errorsNegCorrMatSummedGammaC[ConstnBinsGamma];

	Double_t errorsNegErrSummedGammaC[ConstnBinsGamma];
	Double_t errorsNegErrCorrSummedGammaC[ConstnBinsGamma];

	Double_t errorsPosSummedGammaC[ConstnBinsGamma];
	Double_t errorsPosCorrSummedGammaC[ConstnBinsGamma];
	Double_t errorsPosCorrMatSummedGammaC[ConstnBinsGamma];

	Double_t errorsPosErrSummedGammaC[ConstnBinsGamma];
	Double_t errorsPosErrCorrSummedGammaC[ConstnBinsGamma];


	Double_t errorsMeanGamma[ErrorGammaSpec][ConstnBinsGamma];
	Double_t errorsMeanCorrGamma[ErrorGammaSpec][ConstnBinsGamma];
	Double_t errorsMeanSummedGamma[ConstnBinsGamma];
	Double_t errorsMeanCorrSummedGamma[ConstnBinsGamma];
	Double_t errorsMeanCorrMatSummedGamma[ConstnBinsGamma];
	Double_t errorsMeanCorrMatSummedGammaWithoutMat[ConstnBinsGamma];
	Double_t errorsMaterialBudget[ConstnBinsGamma];

	Double_t errorsMeanErrGamma[ErrorGammaSpec][ConstnBinsGamma];
	Double_t errorsMeanErrCorrGamma[ErrorGammaSpec][ConstnBinsGamma];
	Double_t errorsMeanErrSummedGamma[ConstnBinsGamma];
	Double_t errorsMeanErrCorrSummedGamma[ConstnBinsGamma];
	Double_t errorsMeanErrCorrMatSummedGamma[ConstnBinsGamma];
	Double_t errorsMeanErrCorrMatSummedGammaWithoutMat[ConstnBinsGamma];
	Double_t errorsErrMaterialBudget[ConstnBinsGamma];

	Double_t errorsMeanSummedGammaA[ConstnBinsGamma];
	Double_t errorsMeanCorrSummedGammaA[ConstnBinsGamma];
	Double_t errorsMeanCorrMatSummedGammaA[ConstnBinsGamma];

	Double_t errorsMeanErrSummedGammaA[ConstnBinsGamma];
	Double_t errorsMeanErrCorrSummedGammaA[ConstnBinsGamma];
	Double_t errorsMeanErrCorrMatSummedGammaA[ConstnBinsGamma];

	Double_t errorsMeanSummedGammaB[ConstnBinsGamma];
	Double_t errorsMeanCorrSummedGammaB[ConstnBinsGamma];
	Double_t errorsMeanCorrMatSummedGammaB[ConstnBinsGamma];

	Double_t errorsMeanErrSummedGammaB[ConstnBinsGamma];
	Double_t errorsMeanErrCorrSummedGammaB[ConstnBinsGamma];
	Double_t errorsMeanErrCorrMatSummedGammaB[ConstnBinsGamma];

	Double_t errorsMeanSummedGammaC[ConstnBinsGamma];
	Double_t errorsMeanCorrSummedGammaC[ConstnBinsGamma];
	Double_t errorsMeanCorrMatSummedGammaC[ConstnBinsGamma];

	Double_t errorsMeanErrSummedGammaC[ConstnBinsGamma];
	Double_t errorsMeanErrCorrSummedGammaC[ConstnBinsGamma];
	Double_t errorsMeanErrCorrMatSummedGammaC[ConstnBinsGamma];

	for (Int_t l = 0; l < ConstnBinsGamma; l++){
		errorsPosSummedGamma[l] = 0.;
		errorsNegSummedGamma[l] = 0.;
		errorsMeanSummedGamma[l] = 0.;
		errorsPosCorrSummedGamma[l] = 0.;
		errorsNegCorrSummedGamma[l] = 0.;
		errorsMeanCorrSummedGamma[l] = 0.;

		errorsPosSummedGammaA[l] = 0.;
		errorsNegSummedGammaA[l] = 0.;
		errorsMeanSummedGammaA[l] = 0.;
		errorsPosCorrSummedGammaA[l] = 0.;
		errorsNegCorrSummedGammaA[l] = 0.;
		errorsMeanCorrSummedGammaA[l] = 0.;

		errorsPosSummedGammaB[l] = 0.;
		errorsNegSummedGammaB[l] = 0.;
		errorsMeanSummedGammaB[l] = 0.;
		errorsPosCorrSummedGammaB[l] = 0.;
		errorsNegCorrSummedGammaB[l] = 0.;
		errorsMeanCorrSummedGammaB[l] = 0.;

		errorsPosSummedGammaC[l] = 0.;
		errorsNegSummedGammaC[l] = 0.;
		errorsMeanSummedGammaC[l] = 0.;
		errorsPosCorrSummedGammaC[l] = 0.;
		errorsNegCorrSummedGammaC[l] = 0.;
		errorsMeanCorrSummedGammaC[l] = 0.;
	}

	TGraphErrors* negativeErrorsDoubleRatio[ErrorDoubleRatio];
	TGraphErrors* negativeErrorsSummedDoubleRatio;
	TGraphErrors* positiveErrorsDoubleRatio[ErrorDoubleRatio];
	TGraphErrors* positiveErrorsSummedDoubleRatio;
	TGraphErrors* negativeErrorsCorrDoubleRatio[ErrorDoubleRatio];
	TGraphErrors* negativeErrorsCorrSummedDoubleRatio;
	TGraphErrors* positiveErrorsCorrDoubleRatio[ErrorDoubleRatio];
	TGraphErrors* positiveErrorsCorrSummedDoubleRatio;
	TGraphErrors* meanErrorsDoubleRatio[ErrorDoubleRatio];
	TGraphErrors* meanErrorsSummedDoubleRatio;
	TGraphErrors* meanErrorsCorrDoubleRatio[ErrorDoubleRatio];
	TGraphErrors* meanErrorsCorrSummedDoubleRatio;
	TGraphErrors* meanErrorsCorrSummedIncMatDoubleRatio;

	TGraphErrors* negativeErrorsDoubleRatioPi0Fit[ErrorDoubleRatio];
	TGraphErrors* negativeErrorsSummedDoubleRatioPi0Fit;
	TGraphErrors* positiveErrorsDoubleRatioPi0Fit[ErrorDoubleRatio];
	TGraphErrors* positiveErrorsSummedDoubleRatioPi0Fit;
	TGraphErrors* negativeErrorsCorrDoubleRatioPi0Fit[ErrorDoubleRatio];
	TGraphErrors* negativeErrorsCorrSummedDoubleRatioPi0Fit;
	TGraphErrors* positiveErrorsCorrDoubleRatioPi0Fit[ErrorDoubleRatio];
	TGraphErrors* positiveErrorsCorrSummedDoubleRatioPi0Fit;
	TGraphErrors* meanErrorsDoubleRatioPi0Fit[ErrorDoubleRatio];
	TGraphErrors* meanErrorsSummedDoubleRatioPi0Fit;
	TGraphErrors* meanErrorsCorrDoubleRatioPi0Fit[ErrorDoubleRatio];
	TGraphErrors* meanErrorsCorrSummedDoubleRatioPi0Fit;
	TGraphErrors* meanErrorsCorrSummedIncMatDoubleRatioPi0Fit;
	TGraphErrors* meanErrorsCorrSummedWithoutMatDoubleRatioPi0Fit;
	TGraphErrors* meanErrorsCorrSummedMaterialDoubleRatioPi0Fit;
	TGraphErrors* meanErrorsCorrSummedMaterialDoubleRatioPi0FitCocktail;

	TGraphErrors* negativeErrorsSummedDoubleRatioA;
	TGraphErrors* positiveErrorsSummedDoubleRatioA;
	TGraphErrors* negativeErrorsCorrSummedDoubleRatioA;
	TGraphErrors* positiveErrorsCorrSummedDoubleRatioA;
	TGraphErrors* meanErrorsSummedDoubleRatioA;
	TGraphErrors* meanErrorsCorrSummedDoubleRatioA;
	TGraphErrors* meanErrorsCorrSummedIncMatDoubleRatioA;

	TGraphErrors* negativeErrorsSummedDoubleRatioPi0FitA;
	TGraphErrors* positiveErrorsSummedDoubleRatioPi0FitA;
	TGraphErrors* negativeErrorsCorrSummedDoubleRatioPi0FitA;
	TGraphErrors* positiveErrorsCorrSummedDoubleRatioPi0FitA;
	TGraphErrors* meanErrorsSummedDoubleRatioPi0FitA;
	TGraphErrors* meanErrorsCorrSummedDoubleRatioPi0FitA;
	TGraphErrors* meanErrorsCorrSummedIncMatDoubleRatioPi0FitA;
	
	TGraphErrors* negativeErrorsSummedDoubleRatioB;
	TGraphErrors* positiveErrorsSummedDoubleRatioB;
	TGraphErrors* negativeErrorsCorrSummedDoubleRatioB;
	TGraphErrors* positiveErrorsCorrSummedDoubleRatioB;
	TGraphErrors* meanErrorsSummedDoubleRatioB;
	TGraphErrors* meanErrorsCorrSummedDoubleRatioB;
	TGraphErrors* meanErrorsCorrSummedIncMatDoubleRatioB;

	TGraphErrors* negativeErrorsSummedDoubleRatioPi0FitB;
	TGraphErrors* positiveErrorsSummedDoubleRatioPi0FitB;
	TGraphErrors* negativeErrorsCorrSummedDoubleRatioPi0FitB;
	TGraphErrors* positiveErrorsCorrSummedDoubleRatioPi0FitB;
	TGraphErrors* meanErrorsSummedDoubleRatioPi0FitB;
	TGraphErrors* meanErrorsCorrSummedDoubleRatioPi0FitB;
	TGraphErrors* meanErrorsCorrSummedIncMatDoubleRatioPi0FitB;
	
	TGraphErrors* negativeErrorsSummedDoubleRatioC;
	TGraphErrors* positiveErrorsSummedDoubleRatioC;
	TGraphErrors* negativeErrorsCorrSummedDoubleRatioC;
	TGraphErrors* positiveErrorsCorrSummedDoubleRatioC;
	TGraphErrors* meanErrorsSummedDoubleRatioC;
	TGraphErrors* meanErrorsCorrSummedDoubleRatioC;
	TGraphErrors* meanErrorsCorrSummedIncMatDoubleRatioC;

	TGraphErrors* negativeErrorsSummedDoubleRatioPi0FitC;
	TGraphErrors* positiveErrorsSummedDoubleRatioPi0FitC;
	TGraphErrors* negativeErrorsCorrSummedDoubleRatioPi0FitC;
	TGraphErrors* positiveErrorsCorrSummedDoubleRatioPi0FitC;
	TGraphErrors* meanErrorsSummedDoubleRatioPi0FitC;
	TGraphErrors* meanErrorsCorrSummedDoubleRatioPi0FitC;
	TGraphErrors* meanErrorsCorrSummedIncMatDoubleRatioPi0FitC;


	TGraphErrors* negativeErrorsIncRatio[ErrorIncRatio];
	TGraphErrors* negativeErrorsSummedIncRatio;
	TGraphErrors* positiveErrorsIncRatio[ErrorIncRatio];
	TGraphErrors* positiveErrorsSummedIncRatio;
	TGraphErrors* negativeErrorsCorrIncRatio[ErrorIncRatio];
	TGraphErrors* negativeErrorsCorrSummedIncRatio;
	TGraphErrors* positiveErrorsCorrIncRatio[ErrorIncRatio];
	TGraphErrors* positiveErrorsCorrSummedIncRatio;
	TGraphErrors* meanErrorsIncRatio[ErrorIncRatio];
	TGraphErrors* meanErrorsSummedIncRatio;
	TGraphErrors* meanErrorsCorrIncRatio[ErrorIncRatio];
	TGraphErrors* meanErrorsCorrSummedIncRatio;
	TGraphErrors* meanErrorsCorrSummedIncMatIncRatio;

	TGraphErrors* negativeErrorsIncRatioPi0Fit[ErrorIncRatio];
	TGraphErrors* negativeErrorsSummedIncRatioPi0Fit;
	TGraphErrors* positiveErrorsIncRatioPi0Fit[ErrorIncRatio];
	TGraphErrors* positiveErrorsSummedIncRatioPi0Fit;
	TGraphErrors* negativeErrorsCorrIncRatioPi0Fit[ErrorIncRatio];
	TGraphErrors* negativeErrorsCorrSummedIncRatioPi0Fit;
	TGraphErrors* positiveErrorsCorrIncRatioPi0Fit[ErrorIncRatio];
	TGraphErrors* positiveErrorsCorrSummedIncRatioPi0Fit;
	TGraphErrors* meanErrorsIncRatioPi0Fit[ErrorIncRatio];
	TGraphErrors* meanErrorsSummedIncRatioPi0Fit;
	TGraphErrors* meanErrorsCorrIncRatioPi0Fit[ErrorIncRatio];
	TGraphErrors* meanErrorsCorrSummedIncRatioPi0Fit;
	TGraphErrors* meanErrorsCorrSummedIncMatIncRatioPi0Fit;

	TGraphErrors* negativeErrorsSummedIncRatioA;
	TGraphErrors* positiveErrorsSummedIncRatioA;
	TGraphErrors* negativeErrorsCorrSummedIncRatioA;
	TGraphErrors* positiveErrorsCorrSummedIncRatioA;
	TGraphErrors* meanErrorsSummedIncRatioA;
	TGraphErrors* meanErrorsCorrSummedIncRatioA;
	TGraphErrors* meanErrorsCorrSummedIncMatIncRatioA;

	TGraphErrors* negativeErrorsSummedIncRatioPi0FitA;
	TGraphErrors* positiveErrorsSummedIncRatioPi0FitA;
	TGraphErrors* negativeErrorsCorrSummedIncRatioPi0FitA;
	TGraphErrors* positiveErrorsCorrSummedIncRatioPi0FitA;
	TGraphErrors* meanErrorsSummedIncRatioPi0FitA;
	TGraphErrors* meanErrorsCorrSummedIncRatioPi0FitA;
	TGraphErrors* meanErrorsCorrSummedIncMatIncRatioPi0FitA;

	TGraphErrors* negativeErrorsSummedIncRatioB;
	TGraphErrors* positiveErrorsSummedIncRatioB;
	TGraphErrors* negativeErrorsCorrSummedIncRatioB;
	TGraphErrors* positiveErrorsCorrSummedIncRatioB;
	TGraphErrors* meanErrorsSummedIncRatioB;
	TGraphErrors* meanErrorsCorrSummedIncRatioB;
	TGraphErrors* meanErrorsCorrSummedIncMatIncRatioB;

	TGraphErrors* negativeErrorsSummedIncRatioPi0FitB;
	TGraphErrors* positiveErrorsSummedIncRatioPi0FitB;
	TGraphErrors* negativeErrorsCorrSummedIncRatioPi0FitB;
	TGraphErrors* positiveErrorsCorrSummedIncRatioPi0FitB;
	TGraphErrors* meanErrorsSummedIncRatioPi0FitB;
	TGraphErrors* meanErrorsCorrSummedIncRatioPi0FitB;
	TGraphErrors* meanErrorsCorrSummedIncMatIncRatioPi0FitB;

	TGraphErrors* negativeErrorsSummedIncRatioC;
	TGraphErrors* positiveErrorsSummedIncRatioC;
	TGraphErrors* negativeErrorsCorrSummedIncRatioC;
	TGraphErrors* positiveErrorsCorrSummedIncRatioC;
	TGraphErrors* meanErrorsSummedIncRatioC;
	TGraphErrors* meanErrorsCorrSummedIncRatioC;
	TGraphErrors* meanErrorsCorrSummedIncMatIncRatioC;

	TGraphErrors* negativeErrorsSummedIncRatioPi0FitC;
	TGraphErrors* positiveErrorsSummedIncRatioPi0FitC;
	TGraphErrors* negativeErrorsCorrSummedIncRatioPi0FitC;
	TGraphErrors* positiveErrorsCorrSummedIncRatioPi0FitC;
	TGraphErrors* meanErrorsSummedIncRatioPi0FitC;
	TGraphErrors* meanErrorsCorrSummedIncRatioPi0FitC;
	TGraphErrors* meanErrorsCorrSummedIncMatIncRatioPi0FitC;

	TGraphErrors* negativeErrorsPi0[ErrorPi0];
	TGraphErrors* negativeErrorsSummedPi0;
	TGraphErrors* positiveErrorsPi0[ErrorPi0];
	TGraphErrors* positiveErrorsSummedPi0;
	TGraphErrors* negativeErrorsCorrPi0[ErrorPi0];
	TGraphErrors* negativeErrorsCorrSummedPi0;
	TGraphErrors* positiveErrorsCorrPi0[ErrorPi0];
	TGraphErrors* positiveErrorsCorrSummedPi0;
	TGraphErrors* meanErrorsPi0[ErrorPi0];
	TGraphErrors* meanErrorsSummedPi0;
	TGraphErrors* meanErrorsCorrPi0[ErrorPi0];
	TGraphErrors* meanErrorsCorrSummedPi0;
	TGraphErrors* meanErrorsCorrSummedIncMatPi0;
	TGraphErrors* meanErrorsCorrSummedIncMatPi0WithoutMat;

	TGraphErrors* negativeErrorsPi0Fit[ErrorPi0];
	TGraphErrors* negativeErrorsSummedPi0Fit;
	TGraphErrors* positiveErrorsPi0Fit[ErrorPi0];
	TGraphErrors* positiveErrorsSummedPi0Fit;
	TGraphErrors* negativeErrorsCorrPi0Fit[ErrorPi0];
	TGraphErrors* negativeErrorsCorrSummedPi0Fit;
	TGraphErrors* positiveErrorsCorrPi0Fit[ErrorPi0];
	TGraphErrors* positiveErrorsCorrSummedPi0Fit;
	TGraphErrors* meanErrorsPi0Fit[ErrorPi0];
	TGraphErrors* meanErrorsSummedPi0Fit;
	TGraphErrors* meanErrorsCorrPi0Fit[ErrorPi0];
	TGraphErrors* meanErrorsCorrSummedPi0Fit;
	TGraphErrors* meanErrorsCorrSummedIncMatPi0Fit;
	TGraphErrors* meanErrorsCorrSummedIncMatPi0FitWithoutMat;

	TGraphErrors* negativeErrorsGamma[ErrorGammaSpec];
	TGraphErrors* negativeErrorsSummedGamma;
	TGraphErrors* positiveErrorsGamma[ErrorGammaSpec];
	TGraphErrors* positiveErrorsSummedGamma;
	TGraphErrors* negativeErrorsCorrGamma[ErrorGammaSpec];
	TGraphErrors* negativeErrorsCorrSummedGamma;
	TGraphErrors* positiveErrorsCorrGamma[ErrorGammaSpec];
	TGraphErrors* positiveErrorsCorrSummedGamma;
	TGraphErrors* meanErrorsGamma[ErrorGammaSpec];
	TGraphErrors* meanErrorsSummedGamma;
	TGraphErrors* meanErrorsCorrGamma[ErrorGammaSpec];
	TGraphErrors* meanErrorsCorrSummedGamma;
	TGraphErrors* meanErrorsCorrSummedIncMatGamma;
	TGraphErrors* meanErrorsCorrSummedWithoutMatGamma;
	TGraphErrors* meanErrorsCorrSummedMaterialGamma;

	TGraphErrors* negativeErrorsSummedGammaA;
	TGraphErrors* positiveErrorsSummedGammaA;
	TGraphErrors* negativeErrorsCorrSummedGammaA;
	TGraphErrors* positiveErrorsCorrSummedGammaA;
	TGraphErrors* meanErrorsSummedGammaA;
	TGraphErrors* meanErrorsCorrSummedGammaA;
	TGraphErrors* meanErrorsCorrSummedIncMatGammaA;

	TGraphErrors* negativeErrorsSummedGammaB;
	TGraphErrors* positiveErrorsSummedGammaB;
	TGraphErrors* negativeErrorsCorrSummedGammaB;
	TGraphErrors* positiveErrorsCorrSummedGammaB;
	TGraphErrors* meanErrorsSummedGammaB;
	TGraphErrors* meanErrorsCorrSummedGammaB;
	TGraphErrors* meanErrorsCorrSummedIncMatGammaB;

	TGraphErrors* negativeErrorsSummedGammaC;
	TGraphErrors* positiveErrorsSummedGammaC;
	TGraphErrors* negativeErrorsCorrSummedGammaC;
	TGraphErrors* positiveErrorsCorrSummedGammaC;
	TGraphErrors* meanErrorsSummedGammaC;
	TGraphErrors* meanErrorsCorrSummedGammaC;
	TGraphErrors* meanErrorsCorrSummedIncMatGammaC;



	Int_t nBinsGraphDoubeRatio = 0;
	Int_t nBinsGraphGamma = 0;
	for(Int_t i = 0; i<ErrorDoubleRatioToRun; i++){
		cout<<"------------------------------------------------>  "<<nameCutVariationsDoubleRatio[i]<<endl;

		graphPosErrorsDoubleRatio[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatio_SystErrorRelPos_%s%s_%s",nameCutVariationsDoubleRatio[i].Data(),additionalNameOutput.Data(),additionalName.Data()));
		cout << "here" << endl;
		graphNegErrorsDoubleRatio[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatio_SystErrorRelNeg_%s%s_%s",nameCutVariationsDoubleRatio[i].Data(),additionalNameOutput.Data(),additionalName.Data()));
		cout << "here" << endl;
		if(!nameCutVariationsDoubleRatio[i].CompareTo("CocktailEtaNorm")){
			graphPosErrorsDoubleRatio[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatio_SystErrorRelPos_%s%s_%s",nameCutVariationsDoubleRatio[i].Data(),additionalNameOutput.Data(),additionalName.Data()));
			graphNegErrorsDoubleRatio[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatio_SystErrorRelNeg_%s%s_%s",nameCutVariationsDoubleRatio[i].Data(),additionalNameOutput.Data(),additionalName.Data()));
			for(Int_t jj = 0;jj<graphPosErrorsDoubleRatio[i]->GetN();jj++){
				graphPosErrorsDoubleRatio[i]->SetPoint(jj, graphPosErrorsDoubleRatio[i]->GetX()[jj], 1.4);
				graphNegErrorsDoubleRatio[i]->SetPoint(jj, graphNegErrorsDoubleRatio[i]->GetX()[jj], 1.4);
			}
		}
		cout << "here" << endl;
// 		graphPosErrorsDoubleRatio[i]->RemovePoint(0);
// 		graphNegErrorsDoubleRatio[i]->RemovePoint(0);
// 		graphPosErrorsDoubleRatio[i]->RemovePoint(0);
// 		graphNegErrorsDoubleRatio[i]->RemovePoint(0);
		cout << "here" << endl;
		
		graphPosErrorsDoubleRatioPi0Fit[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatioFit_SystErrorRelPos_%s%s_%s",nameCutVariationsDoubleRatio[i].Data(),additionalNameOutput.Data(),additionalName.Data()));
		cout << "here" << endl;
		graphNegErrorsDoubleRatioPi0Fit[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatioFit_SystErrorRelNeg_%s%s_%s",nameCutVariationsDoubleRatio[i].Data(),additionalNameOutput.Data(),additionalName.Data()));
		cout << "here" << endl;
		if(!nameCutVariationsDoubleRatio[i].CompareTo("CocktailEtaNorm")){
			graphPosErrorsDoubleRatioPi0Fit[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatioFit_SystErrorRelPos_%s%s_%s",nameCutVariationsDoubleRatio[i].Data(),additionalNameOutput.Data(),additionalName.Data()));
			graphNegErrorsDoubleRatioPi0Fit[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("DoubleRatioFit_SystErrorRelNeg_%s%s_%s",nameCutVariationsDoubleRatio[i].Data(),additionalNameOutput.Data(),additionalName.Data()));
			for(Int_t jj = 0;jj<graphPosErrorsDoubleRatioPi0Fit[i]->GetN();jj++){
				graphPosErrorsDoubleRatioPi0Fit[i]->SetPoint(jj, graphPosErrorsDoubleRatioPi0Fit[i]->GetX()[jj], 1.4);
				graphNegErrorsDoubleRatioPi0Fit[i]->SetPoint(jj, graphNegErrorsDoubleRatioPi0Fit[i]->GetX()[jj], 1.4);
			}
		}

// 		graphPosErrorsDoubleRatioPi0Fit[i]->RemovePoint(0);
// 		graphNegErrorsDoubleRatioPi0Fit[i]->RemovePoint(0);
// 		graphPosErrorsDoubleRatioPi0Fit[i]->RemovePoint(0);
// 		graphNegErrorsDoubleRatioPi0Fit[i]->RemovePoint(0);
		cout << "here" << endl;
		Int_t type = 0;
		if(!nameCutVariationsDoubleRatio[i].CompareTo("Chi2") || !nameCutVariationsDoubleRatio[i].CompareTo("dEdxE") || !nameCutVariationsDoubleRatio[i].CompareTo("PsiPair")  || !nameCutVariationsDoubleRatio[i].CompareTo("Cluster") || !nameCutVariationsDoubleRatio[i].CompareTo("CocktailParam") || !nameCutVariationsDoubleRatio[i].CompareTo("CocktailEta") || !nameCutVariationsDoubleRatio[i].CompareTo("TOF") || !nameCutVariationsDoubleRatio[i].CompareTo("SinglePt") || !nameCutVariationsDoubleRatio[i].CompareTo("Alpha") ){ // A Errors
			type = 1;
		}
		if(!nameCutVariationsDoubleRatio[i].CompareTo("qT") || !nameCutVariationsDoubleRatio[i].CompareTo("dEdxPi") || !nameCutVariationsDoubleRatio[i].CompareTo("IntRange")){ // B Errors
			type = 2;
		}
		if(!nameCutVariationsDoubleRatio[i].CompareTo("CocktailEtaNorm")){ // C Errors
			type = 3;
		}


		if(i==0){
			ptBinsDoubleRatio = graphPosErrorsDoubleRatio[i]->GetX();
			nBinsGraphDoubeRatio = graphPosErrorsDoubleRatio[i]->GetN();
			ptBinsDoubleRatioErr = graphPosErrorsDoubleRatio[i]->GetEXhigh();
		}

		errorsNegDoubleRatio[i] = graphNegErrorsDoubleRatio[i]->GetY();
		errorsNegErrDoubleRatio[i] = graphNegErrorsDoubleRatio[i]->GetEYhigh();
		errorsPosDoubleRatio[i] = graphPosErrorsDoubleRatio[i]->GetY();
		errorsPosErrDoubleRatio[i] = graphPosErrorsDoubleRatio[i]->GetEYhigh();

		errorsPosDoubleRatioPi0Fit[i] = graphPosErrorsDoubleRatioPi0Fit[i]->GetY();
		errorsPosErrDoubleRatioPi0Fit[i] = graphPosErrorsDoubleRatioPi0Fit[i]->GetEYhigh();
		errorsNegDoubleRatioPi0Fit[i] = graphNegErrorsDoubleRatioPi0Fit[i]->GetY();
		errorsNegErrDoubleRatioPi0Fit[i] = graphNegErrorsDoubleRatioPi0Fit[i]->GetEYhigh();

		CalculateMeanSysErr(errorsMeanDoubleRatio[i], errorsMeanErrDoubleRatio[i], errorsPosDoubleRatio[i], errorsNegDoubleRatio[i], ConstnBinsDoubleRatio);
		CalculateMeanSysErr(errorsMeanDoubleRatioPi0Fit[i], errorsMeanErrDoubleRatioPi0Fit[i], errorsPosDoubleRatioPi0Fit[i], errorsNegDoubleRatioPi0Fit[i], ConstnBinsDoubleRatio);

		CorrectSystematicErrorsWithMean(errorsPosDoubleRatio[i],errorsPosErrDoubleRatio[i], errorsPosCorrDoubleRatio[i], errorsPosErrCorrDoubleRatio[i], ConstnBinsDoubleRatio);
		CorrectSystematicErrorsWithMean(errorsNegDoubleRatio[i],errorsNegErrDoubleRatio[i], errorsNegCorrDoubleRatio[i], errorsNegErrCorrDoubleRatio[i], ConstnBinsDoubleRatio);
		CorrectSystematicErrorsWithMean(errorsMeanDoubleRatio[i], errorsMeanErrDoubleRatio[i], errorsMeanCorrDoubleRatio[i], errorsMeanErrCorrDoubleRatio[i], ConstnBinsDoubleRatio);

		CorrectSystematicErrorsWithMean(errorsPosDoubleRatioPi0Fit[i],errorsPosErrDoubleRatioPi0Fit[i], errorsPosCorrDoubleRatioPi0Fit[i], errorsPosErrCorrDoubleRatioPi0Fit[i], ConstnBinsDoubleRatio);
		CorrectSystematicErrorsWithMean(errorsNegDoubleRatioPi0Fit[i],errorsNegErrDoubleRatioPi0Fit[i], errorsNegCorrDoubleRatioPi0Fit[i], errorsNegErrCorrDoubleRatioPi0Fit[i], ConstnBinsDoubleRatio);
		CorrectSystematicErrorsWithMean(errorsMeanDoubleRatioPi0Fit[i], errorsMeanErrDoubleRatioPi0Fit[i], errorsMeanCorrDoubleRatioPi0Fit[i], errorsMeanErrCorrDoubleRatioPi0Fit[i], ConstnBinsDoubleRatio);


		negativeErrorsDoubleRatio[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegDoubleRatio[i] ,ptBinsDoubleRatioErr ,errorsNegErrDoubleRatio[i] );
		meanErrorsDoubleRatio[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanDoubleRatio[i] ,ptBinsDoubleRatioErr ,errorsMeanErrDoubleRatio[i] );
		positiveErrorsDoubleRatio[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosDoubleRatio[i] ,ptBinsDoubleRatioErr ,errorsPosErrDoubleRatio[i] );
		negativeErrorsCorrDoubleRatio[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrDoubleRatio[i] ,ptBinsDoubleRatioErr ,errorsNegErrCorrDoubleRatio[i] );
		meanErrorsCorrDoubleRatio[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrDoubleRatio[i] ,ptBinsDoubleRatioErr ,errorsMeanErrCorrDoubleRatio[i] );
		positiveErrorsCorrDoubleRatio[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrDoubleRatio[i] ,ptBinsDoubleRatioErr ,errorsPosErrCorrDoubleRatio[i] );

		negativeErrorsDoubleRatioPi0Fit[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegDoubleRatioPi0Fit[i] ,ptBinsDoubleRatioErr ,errorsNegErrDoubleRatioPi0Fit[i] );
		meanErrorsDoubleRatioPi0Fit[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanDoubleRatioPi0Fit[i] ,ptBinsDoubleRatioErr ,errorsMeanErrDoubleRatioPi0Fit[i] );
		positiveErrorsDoubleRatioPi0Fit[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosDoubleRatioPi0Fit[i] ,ptBinsDoubleRatioErr ,errorsPosErrDoubleRatioPi0Fit[i] );
		negativeErrorsCorrDoubleRatioPi0Fit[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrDoubleRatioPi0Fit[i] ,ptBinsDoubleRatioErr ,errorsNegErrCorrDoubleRatioPi0Fit[i] );
		meanErrorsCorrDoubleRatioPi0Fit[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrDoubleRatioPi0Fit[i] ,ptBinsDoubleRatioErr ,errorsMeanErrCorrDoubleRatioPi0Fit[i] );
		positiveErrorsCorrDoubleRatioPi0Fit[i] = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrDoubleRatioPi0Fit[i] ,ptBinsDoubleRatioErr ,errorsPosErrCorrDoubleRatioPi0Fit[i] );


		if(i<ErrorIncRatioToRun){
			cout << "second here" << endl;
			graphPosErrorsIncRatio[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("IncRatio_SystErrorRelPos_%s%s_%s",nameCutVariationsDoubleRatio[i].Data(),additionalNameOutput.Data(),additionalName.Data()));
			cout << Form("IncRatio_SystErrorRelPos_%s%s_%s",nameCutVariationsDoubleRatio[i].Data(),additionalNameOutput.Data(),additionalName.Data()) << endl;
			graphNegErrorsIncRatio[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("IncRatio_SystErrorRelNeg_%s%s_%s",nameCutVariationsDoubleRatio[i].Data(),additionalNameOutput.Data(),additionalName.Data()));
			cout << Form("IncRatio_SystErrorRelNeg_%s%s_%s",nameCutVariationsDoubleRatio[i].Data(),additionalNameOutput.Data(),additionalName.Data()) << endl;
			
// 			graphPosErrorsIncRatio[i]->RemovePoint(0);
// 			graphNegErrorsIncRatio[i]->RemovePoint(0);
// 			graphPosErrorsIncRatio[i]->RemovePoint(0);
// 			graphNegErrorsIncRatio[i]->RemovePoint(0);
			cout << "second here" << endl;
			
			graphPosErrorsIncRatioPi0Fit[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("IncRatioFit_SystErrorRelPos_%s%s_%s",nameCutVariationsDoubleRatio[i].Data(),additionalNameOutput.Data(),additionalName.Data()));
			graphNegErrorsIncRatioPi0Fit[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("IncRatioFit_SystErrorRelNeg_%s%s_%s",nameCutVariationsDoubleRatio[i].Data(),additionalNameOutput.Data(),additionalName.Data()));

// 			graphPosErrorsIncRatioPi0Fit[i]->RemovePoint(0);
// 			graphNegErrorsIncRatioPi0Fit[i]->RemovePoint(0);
// 			graphPosErrorsIncRatioPi0Fit[i]->RemovePoint(0);
// 			graphNegErrorsIncRatioPi0Fit[i]->RemovePoint(0);

			if(i==0){
				ptBinsIncRatio = graphPosErrorsIncRatio[i]->GetX();
				nBinsGraphDoubeRatio = graphPosErrorsIncRatio[i]->GetN();
				ptBinsIncRatioErr = graphPosErrorsIncRatio[i]->GetEXhigh();
			}

			errorsNegIncRatio[i] = graphNegErrorsIncRatio[i]->GetY();
			errorsNegErrIncRatio[i] = graphNegErrorsIncRatio[i]->GetEYhigh();
			errorsPosIncRatio[i] = graphPosErrorsIncRatio[i]->GetY();
			errorsPosErrIncRatio[i] = graphPosErrorsIncRatio[i]->GetEYhigh();

			errorsPosIncRatioPi0Fit[i] = graphPosErrorsIncRatioPi0Fit[i]->GetY();
			errorsPosErrIncRatioPi0Fit[i] = graphPosErrorsIncRatioPi0Fit[i]->GetEYhigh();
			errorsNegIncRatioPi0Fit[i] = graphNegErrorsIncRatioPi0Fit[i]->GetY();
			errorsNegErrIncRatioPi0Fit[i] = graphNegErrorsIncRatioPi0Fit[i]->GetEYhigh();


			CalculateMeanSysErr(errorsMeanIncRatio[i], errorsMeanErrIncRatio[i], errorsPosIncRatio[i], errorsNegIncRatio[i], ConstnBinsIncRatio);
			CalculateMeanSysErr(errorsMeanIncRatioPi0Fit[i], errorsMeanErrIncRatioPi0Fit[i], errorsPosIncRatioPi0Fit[i], errorsNegIncRatioPi0Fit[i], ConstnBinsIncRatio);

			CorrectSystematicErrorsWithMean(errorsPosIncRatio[i],errorsPosErrIncRatio[i], errorsPosCorrIncRatio[i], errorsPosErrCorrIncRatio[i], ConstnBinsIncRatio);
			CorrectSystematicErrorsWithMean(errorsNegIncRatio[i],errorsNegErrIncRatio[i], errorsNegCorrIncRatio[i], errorsNegErrCorrIncRatio[i], ConstnBinsIncRatio);
			CorrectSystematicErrorsWithMean(errorsMeanIncRatio[i], errorsMeanErrIncRatio[i], errorsMeanCorrIncRatio[i], errorsMeanErrCorrIncRatio[i], ConstnBinsIncRatio);

			CorrectSystematicErrorsWithMean(errorsPosIncRatioPi0Fit[i],errorsPosErrIncRatioPi0Fit[i], errorsPosCorrIncRatioPi0Fit[i], errorsPosErrCorrIncRatioPi0Fit[i], ConstnBinsIncRatio);
			CorrectSystematicErrorsWithMean(errorsNegIncRatioPi0Fit[i],errorsNegErrIncRatioPi0Fit[i], errorsNegCorrIncRatioPi0Fit[i], errorsNegErrCorrIncRatioPi0Fit[i], ConstnBinsIncRatio);
			CorrectSystematicErrorsWithMean(errorsMeanIncRatioPi0Fit[i], errorsMeanErrIncRatioPi0Fit[i], errorsMeanCorrIncRatioPi0Fit[i], errorsMeanErrCorrIncRatioPi0Fit[i], ConstnBinsIncRatio);

			negativeErrorsIncRatio[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegIncRatio[i] ,ptBinsIncRatioErr ,errorsNegErrIncRatio[i] );
			meanErrorsIncRatio[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanIncRatio[i] ,ptBinsIncRatioErr ,errorsMeanErrIncRatio[i] );
			positiveErrorsIncRatio[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosIncRatio[i] ,ptBinsIncRatioErr ,errorsPosErrIncRatio[i] );
			negativeErrorsCorrIncRatio[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrIncRatio[i] ,ptBinsIncRatioErr ,errorsNegErrCorrIncRatio[i] );
			meanErrorsCorrIncRatio[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrIncRatio[i] ,ptBinsIncRatioErr ,errorsMeanErrCorrIncRatio[i] );
			positiveErrorsCorrIncRatio[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrIncRatio[i] ,ptBinsIncRatioErr ,errorsPosErrCorrIncRatio[i] );

			negativeErrorsIncRatioPi0Fit[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegIncRatioPi0Fit[i] ,ptBinsIncRatioErr ,errorsNegErrIncRatioPi0Fit[i] );
			meanErrorsIncRatioPi0Fit[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanIncRatioPi0Fit[i] ,ptBinsIncRatioErr ,errorsMeanErrIncRatioPi0Fit[i] );
			positiveErrorsIncRatioPi0Fit[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosIncRatioPi0Fit[i] ,ptBinsIncRatioErr ,errorsPosErrIncRatioPi0Fit[i] );
			negativeErrorsCorrIncRatioPi0Fit[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrIncRatioPi0Fit[i] ,ptBinsIncRatioErr ,errorsNegErrCorrIncRatioPi0Fit[i] );
			meanErrorsCorrIncRatioPi0Fit[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrIncRatioPi0Fit[i] ,ptBinsIncRatioErr ,errorsMeanErrCorrIncRatioPi0Fit[i] );
			positiveErrorsCorrIncRatioPi0Fit[i] = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrIncRatioPi0Fit[i] ,ptBinsIncRatioErr ,errorsPosErrCorrIncRatioPi0Fit[i] );


			for (Int_t l = 0; l < ConstnBinsIncRatio; l++){
				errorsPosSummedIncRatio[l] = errorsPosSummedIncRatio[l]+pow(errorsPosIncRatio[i][l],2);
				errorsNegSummedIncRatio[l] = errorsNegSummedIncRatio[l]+ pow(errorsNegIncRatio[i][l],2);
				errorsMeanSummedIncRatio[l] = errorsMeanSummedIncRatio[l]+ pow(errorsMeanIncRatio[i][l],2);
				errorsPosCorrSummedIncRatio[l] = errorsPosCorrSummedIncRatio[l]+pow(errorsPosCorrIncRatio[i][l],2);
				errorsNegCorrSummedIncRatio[l] = errorsNegCorrSummedIncRatio[l] +pow(errorsNegCorrIncRatio[i][l],2);
				errorsMeanCorrSummedIncRatio[l] =errorsMeanCorrSummedIncRatio[l]+ pow(errorsMeanCorrIncRatio[i][l],2);

				errorsPosSummedIncRatioPi0Fit[l] = errorsPosSummedIncRatioPi0Fit[l]+pow(errorsPosIncRatioPi0Fit[i][l],2);
				errorsNegSummedIncRatioPi0Fit[l] = errorsNegSummedIncRatioPi0Fit[l]+ pow(errorsNegIncRatioPi0Fit[i][l],2);
				errorsMeanSummedIncRatioPi0Fit[l] = errorsMeanSummedIncRatioPi0Fit[l]+ pow(errorsMeanIncRatioPi0Fit[i][l],2);
				errorsPosCorrSummedIncRatioPi0Fit[l] = errorsPosCorrSummedIncRatioPi0Fit[l]+pow(errorsPosCorrIncRatioPi0Fit[i][l],2);
				errorsNegCorrSummedIncRatioPi0Fit[l] = errorsNegCorrSummedIncRatioPi0Fit[l] +pow(errorsNegCorrIncRatioPi0Fit[i][l],2);
				errorsMeanCorrSummedIncRatioPi0Fit[l] =errorsMeanCorrSummedIncRatioPi0Fit[l]+ pow(errorsMeanCorrIncRatioPi0Fit[i][l],2);
				if(type == 1){
					errorsPosSummedIncRatioA[l] = errorsPosSummedIncRatioA[l]+pow(errorsPosIncRatio[i][l],2);
					errorsNegSummedIncRatioA[l] = errorsNegSummedIncRatioA[l]+ pow(errorsNegIncRatio[i][l],2);
					errorsMeanSummedIncRatioA[l] = errorsMeanSummedIncRatioA[l]+ pow(errorsMeanIncRatio[i][l],2);
					errorsPosCorrSummedIncRatioA[l] = errorsPosCorrSummedIncRatioA[l]+pow(errorsPosCorrIncRatio[i][l],2);
					errorsNegCorrSummedIncRatioA[l] = errorsNegCorrSummedIncRatioA[l] +pow(errorsNegCorrIncRatio[i][l],2);
					errorsMeanCorrSummedIncRatioA[l] =errorsMeanCorrSummedIncRatioA[l]+ pow(errorsMeanCorrIncRatio[i][l],2);

					errorsPosSummedIncRatioPi0FitA[l] = errorsPosSummedIncRatioPi0FitA[l]+pow(errorsPosIncRatioPi0Fit[i][l],2);
					errorsNegSummedIncRatioPi0FitA[l] = errorsNegSummedIncRatioPi0FitA[l]+ pow(errorsNegIncRatioPi0Fit[i][l],2);
					errorsMeanSummedIncRatioPi0FitA[l] = errorsMeanSummedIncRatioPi0FitA[l]+ pow(errorsMeanIncRatioPi0Fit[i][l],2);
					errorsPosCorrSummedIncRatioPi0FitA[l] = errorsPosCorrSummedIncRatioPi0FitA[l]+pow(errorsPosCorrIncRatioPi0Fit[i][l],2);
					errorsNegCorrSummedIncRatioPi0FitA[l] = errorsNegCorrSummedIncRatioPi0FitA[l] +pow(errorsNegCorrIncRatioPi0Fit[i][l],2);
					errorsMeanCorrSummedIncRatioPi0FitA[l] =errorsMeanCorrSummedIncRatioPi0FitA[l]+ pow(errorsMeanCorrIncRatioPi0Fit[i][l],2);
				}
				if(type == 2){
					errorsPosSummedIncRatioB[l] = errorsPosSummedIncRatioB[l]+pow(errorsPosIncRatio[i][l],2);
					errorsNegSummedIncRatioB[l] = errorsNegSummedIncRatioB[l]+ pow(errorsNegIncRatio[i][l],2);
					errorsMeanSummedIncRatioB[l] = errorsMeanSummedIncRatioB[l]+ pow(errorsMeanIncRatio[i][l],2);
					errorsPosCorrSummedIncRatioB[l] = errorsPosCorrSummedIncRatioB[l]+pow(errorsPosCorrIncRatio[i][l],2);
					errorsNegCorrSummedIncRatioB[l] = errorsNegCorrSummedIncRatioB[l] +pow(errorsNegCorrIncRatio[i][l],2);
					errorsMeanCorrSummedIncRatioB[l] =errorsMeanCorrSummedIncRatioB[l]+ pow(errorsMeanCorrIncRatio[i][l],2);

					errorsPosSummedIncRatioPi0FitB[l] = errorsPosSummedIncRatioPi0FitB[l]+pow(errorsPosIncRatioPi0Fit[i][l],2);
					errorsNegSummedIncRatioPi0FitB[l] = errorsNegSummedIncRatioPi0FitB[l]+ pow(errorsNegIncRatioPi0Fit[i][l],2);
					errorsMeanSummedIncRatioPi0FitB[l] = errorsMeanSummedIncRatioPi0FitB[l]+ pow(errorsMeanIncRatioPi0Fit[i][l],2);
					errorsPosCorrSummedIncRatioPi0FitB[l] = errorsPosCorrSummedIncRatioPi0FitB[l]+pow(errorsPosCorrIncRatioPi0Fit[i][l],2);
					errorsNegCorrSummedIncRatioPi0FitB[l] = errorsNegCorrSummedIncRatioPi0FitB[l] +pow(errorsNegCorrIncRatioPi0Fit[i][l],2);
					errorsMeanCorrSummedIncRatioPi0FitB[l] =errorsMeanCorrSummedIncRatioPi0FitB[l]+ pow(errorsMeanCorrIncRatioPi0Fit[i][l],2);
				}
				if(type == 3){
					errorsPosSummedIncRatioC[l] = errorsPosSummedIncRatioC[l]+pow(errorsPosIncRatio[i][l],2);
					errorsNegSummedIncRatioC[l] = errorsNegSummedIncRatioC[l]+ pow(errorsNegIncRatio[i][l],2);
					errorsMeanSummedIncRatioC[l] = errorsMeanSummedIncRatioC[l]+ pow(errorsMeanIncRatio[i][l],2);
					errorsPosCorrSummedIncRatioC[l] = errorsPosCorrSummedIncRatioC[l]+pow(errorsPosCorrIncRatio[i][l],2);
					errorsNegCorrSummedIncRatioC[l] = errorsNegCorrSummedIncRatioC[l] +pow(errorsNegCorrIncRatio[i][l],2);
					errorsMeanCorrSummedIncRatioC[l] =errorsMeanCorrSummedIncRatioC[l]+ pow(errorsMeanCorrIncRatio[i][l],2);

					errorsPosSummedIncRatioPi0FitC[l] = errorsPosSummedIncRatioPi0FitC[l]+pow(errorsPosIncRatioPi0Fit[i][l],2);
					errorsNegSummedIncRatioPi0FitC[l] = errorsNegSummedIncRatioPi0FitC[l]+ pow(errorsNegIncRatioPi0Fit[i][l],2);
					errorsMeanSummedIncRatioPi0FitC[l] = errorsMeanSummedIncRatioPi0FitC[l]+ pow(errorsMeanIncRatioPi0Fit[i][l],2);
					errorsPosCorrSummedIncRatioPi0FitC[l] = errorsPosCorrSummedIncRatioPi0FitC[l]+pow(errorsPosCorrIncRatioPi0Fit[i][l],2);
					errorsNegCorrSummedIncRatioPi0FitC[l] = errorsNegCorrSummedIncRatioPi0FitC[l] +pow(errorsNegCorrIncRatioPi0Fit[i][l],2);
					errorsMeanCorrSummedIncRatioPi0FitC[l] =errorsMeanCorrSummedIncRatioPi0FitC[l]+ pow(errorsMeanCorrIncRatioPi0Fit[i][l],2);
				}
			}
		}

		if(i<ErrorPi0ToRun){

			graphPosErrorsPi0[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("Pi0_SystErrorRelPos_%s%s_%s",nameCutVariationsDoubleRatio[i].Data(),additionalNameOutput.Data(),additionalName.Data()));
			graphNegErrorsPi0[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("Pi0_SystErrorRelNeg_%s%s_%s",nameCutVariationsDoubleRatio[i].Data(),additionalNameOutput.Data(),additionalName.Data()));

			graphPosErrorsPi0Fit[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("Pi0Fit_SystErrorRelPos_%s%s_%s",nameCutVariationsDoubleRatio[i].Data(),additionalNameOutput.Data(),additionalName.Data()));
			graphNegErrorsPi0Fit[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("Pi0Fit_SystErrorRelNeg_%s%s_%s",nameCutVariationsDoubleRatio[i].Data(),additionalNameOutput.Data(),additionalName.Data()));

// 			graphPosErrorsPi0[i]->RemovePoint(0);
// 			graphNegErrorsPi0[i]->RemovePoint(0);
// 			graphPosErrorsPi0[i]->RemovePoint(0);
// 			graphNegErrorsPi0[i]->RemovePoint(0);
// 
// 			graphPosErrorsPi0Fit[i]->RemovePoint(0);
// 			graphNegErrorsPi0Fit[i]->RemovePoint(0);
// 			graphPosErrorsPi0Fit[i]->RemovePoint(0);
// 			graphNegErrorsPi0Fit[i]->RemovePoint(0);

			if(i==0){
				ptBinsPi0 = graphPosErrorsPi0[i]->GetX();
				nBinsGraphDoubeRatio = graphPosErrorsPi0[i]->GetN();
				ptBinsPi0Err = graphPosErrorsPi0[i]->GetEXhigh();
			}



			errorsNegPi0[i] = graphNegErrorsPi0[i]->GetY();
			errorsNegErrPi0[i] = graphNegErrorsPi0[i]->GetEYhigh();
			errorsPosPi0[i] = graphPosErrorsPi0[i]->GetY();
			errorsPosErrPi0[i] = graphPosErrorsPi0[i]->GetEYhigh();


			errorsPosPi0Fit[i] = graphPosErrorsPi0Fit[i]->GetY();
			errorsPosErrPi0Fit[i] = graphPosErrorsPi0Fit[i]->GetEYhigh();
			errorsNegPi0Fit[i] = graphNegErrorsPi0Fit[i]->GetY();
			errorsNegErrPi0Fit[i] = graphNegErrorsPi0Fit[i]->GetEYhigh();


			CalculateMeanSysErr(errorsMeanPi0[i], errorsMeanErrPi0[i], errorsPosPi0[i], errorsNegPi0[i], ConstnBinsPi0);
			CalculateMeanSysErr(errorsMeanPi0Fit[i], errorsMeanErrPi0Fit[i], errorsPosPi0Fit[i], errorsNegPi0Fit[i], ConstnBinsPi0);

			CorrectSystematicErrorsWithMean(errorsPosPi0[i],errorsPosErrPi0[i], errorsPosCorrPi0[i], errorsPosErrCorrPi0[i], ConstnBinsPi0);
			CorrectSystematicErrorsWithMean(errorsNegPi0[i],errorsNegErrPi0[i], errorsNegCorrPi0[i], errorsNegErrCorrPi0[i], ConstnBinsPi0);
			CorrectSystematicErrorsWithMean(errorsMeanPi0[i], errorsMeanErrPi0[i], errorsMeanCorrPi0[i], errorsMeanErrCorrPi0[i], ConstnBinsPi0);

			CorrectSystematicErrorsWithMean(errorsPosPi0Fit[i],errorsPosErrPi0Fit[i], errorsPosCorrPi0Fit[i], errorsPosErrCorrPi0Fit[i], ConstnBinsPi0);
			CorrectSystematicErrorsWithMean(errorsNegPi0Fit[i],errorsNegErrPi0Fit[i], errorsNegCorrPi0Fit[i], errorsNegErrCorrPi0Fit[i], ConstnBinsPi0);
			CorrectSystematicErrorsWithMean(errorsMeanPi0Fit[i], errorsMeanErrPi0Fit[i], errorsMeanCorrPi0Fit[i], errorsMeanErrCorrPi0Fit[i], ConstnBinsPi0);

			negativeErrorsPi0[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsNegPi0[i] ,ptBinsPi0Err ,errorsNegErrPi0[i] );
			meanErrorsPi0[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanPi0[i] ,ptBinsPi0Err ,errorsMeanErrPi0[i] );
			positiveErrorsPi0[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsPosPi0[i] ,ptBinsPi0Err ,errorsPosErrPi0[i] );
			negativeErrorsCorrPi0[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsNegCorrPi0[i] ,ptBinsPi0Err ,errorsNegErrCorrPi0[i] );
			meanErrorsCorrPi0[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanCorrPi0[i] ,ptBinsPi0Err ,errorsMeanErrCorrPi0[i] );
			positiveErrorsCorrPi0[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsPosCorrPi0[i] ,ptBinsPi0Err ,errorsPosErrCorrPi0[i] );

			negativeErrorsPi0Fit[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsNegPi0Fit[i] ,ptBinsPi0Err ,errorsNegErrPi0Fit[i] );
			meanErrorsPi0Fit[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanPi0Fit[i] ,ptBinsPi0Err ,errorsMeanErrPi0Fit[i] );
			positiveErrorsPi0Fit[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsPosPi0Fit[i] ,ptBinsPi0Err ,errorsPosErrPi0Fit[i] );
			negativeErrorsCorrPi0Fit[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsNegCorrPi0Fit[i] ,ptBinsPi0Err ,errorsNegErrCorrPi0Fit[i] );
			meanErrorsCorrPi0Fit[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanCorrPi0Fit[i] ,ptBinsPi0Err ,errorsMeanErrCorrPi0Fit[i] );
			positiveErrorsCorrPi0Fit[i] = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsPosCorrPi0Fit[i] ,ptBinsPi0Err ,errorsPosErrCorrPi0Fit[i] );

			for (Int_t l = 0; l < ConstnBinsPi0; l++){

				errorsPosSummedPi0[l] = errorsPosSummedPi0[l]+pow(errorsPosPi0[i][l],2);
				errorsNegSummedPi0[l] = errorsNegSummedPi0[l]+ pow(errorsNegPi0[i][l],2);
				errorsMeanSummedPi0[l] = errorsMeanSummedPi0[l]+ pow(errorsMeanPi0[i][l],2);
				errorsPosCorrSummedPi0[l] = errorsPosCorrSummedPi0[l]+pow(errorsPosCorrPi0[i][l],2);
				errorsNegCorrSummedPi0[l] = errorsNegCorrSummedPi0[l] +pow(errorsNegCorrPi0[i][l],2);
				errorsMeanCorrSummedPi0[l] =errorsMeanCorrSummedPi0[l]+ pow(errorsMeanCorrPi0[i][l],2);

				errorsPosSummedPi0Fit[l] = errorsPosSummedPi0Fit[l]+pow(errorsPosPi0Fit[i][l],2);
				errorsNegSummedPi0Fit[l] = errorsNegSummedPi0Fit[l]+ pow(errorsNegPi0Fit[i][l],2);
				errorsMeanSummedPi0Fit[l] = errorsMeanSummedPi0Fit[l]+ pow(errorsMeanPi0Fit[i][l],2);
				errorsPosCorrSummedPi0Fit[l] = errorsPosCorrSummedPi0Fit[l]+pow(errorsPosCorrPi0Fit[i][l],2);
				errorsNegCorrSummedPi0Fit[l] = errorsNegCorrSummedPi0Fit[l] +pow(errorsNegCorrPi0Fit[i][l],2);
				errorsMeanCorrSummedPi0Fit[l] =errorsMeanCorrSummedPi0Fit[l]+ pow(errorsMeanCorrPi0Fit[i][l],2);
			}
		}


		if(i<ErrorGammaSpecToRun){
			graphPosErrorsGamma[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("Gamma_SystErrorRelPos_%s%s_%s",nameCutVariationsDoubleRatio[i].Data(),additionalNameOutput.Data(),additionalName.Data()));
			graphNegErrorsGamma[i] = (TGraphAsymmErrors*)fileSystematics->Get(Form("Gamma_SystErrorRelNeg_%s%s_%s",nameCutVariationsDoubleRatio[i].Data(),additionalNameOutput.Data(),additionalName.Data()));

// 			graphPosErrorsGamma[i]->RemovePoint(0);
// 			graphNegErrorsGamma[i]->RemovePoint(0);
// 			graphPosErrorsGamma[i]->RemovePoint(0);
// 			graphNegErrorsGamma[i]->RemovePoint(0);

			if(i==0){
				ptBinsGamma = graphPosErrorsGamma[i]->GetX();
				nBinsGraphGamma = graphPosErrorsGamma[i]->GetN();
				ptBinsGammaErr = graphPosErrorsGamma[i]->GetEXhigh();
			}

			errorsPosGamma[i] = graphPosErrorsGamma[i]->GetY();
			errorsPosErrGamma[i] = graphPosErrorsGamma[i]->GetEYhigh();
			errorsNegGamma[i] = graphNegErrorsGamma[i]->GetY();
			errorsNegErrGamma[i] = graphNegErrorsGamma[i]->GetEYhigh();

			CalculateMeanSysErr(errorsMeanGamma[i], errorsMeanErrGamma[i], errorsPosGamma[i], errorsNegGamma[i], ConstnBinsGamma);
			CorrectSystematicErrorsWithMean(errorsPosGamma[i],errorsPosErrGamma[i], errorsPosCorrGamma[i], errorsPosErrCorrGamma[i], ConstnBinsGamma);
			CorrectSystematicErrorsWithMean(errorsNegGamma[i],errorsNegErrGamma[i], errorsNegCorrGamma[i], errorsNegErrCorrGamma[i], ConstnBinsGamma);
			CorrectSystematicErrorsWithMean(errorsMeanGamma[i], errorsMeanErrGamma[i], errorsMeanCorrGamma[i], errorsMeanErrCorrGamma[i], ConstnBinsGamma);

			negativeErrorsGamma[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegGamma[i] ,ptBinsGammaErr ,errorsNegErrGamma[i] );
			meanErrorsGamma[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanGamma[i] ,ptBinsGammaErr ,errorsMeanErrGamma[i] );
			positiveErrorsGamma[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosGamma[i] ,ptBinsGammaErr ,errorsPosErrGamma[i] );
			negativeErrorsCorrGamma[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegCorrGamma[i] ,ptBinsGammaErr ,errorsNegErrCorrGamma[i] );
			meanErrorsCorrGamma[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanCorrGamma[i] ,ptBinsGammaErr ,errorsMeanErrCorrGamma[i] );
			positiveErrorsCorrGamma[i] = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosCorrGamma[i] ,ptBinsGammaErr ,errorsPosErrCorrGamma[i] );



			for (Int_t l = 0; l < ConstnBinsGamma; l++){
				errorsPosSummedGamma[l] = errorsPosSummedGamma[l]+pow(errorsPosGamma[i][l],2);
				errorsNegSummedGamma[l] = errorsNegSummedGamma[l]+ pow(errorsNegGamma[i][l],2);
				errorsMeanSummedGamma[l] = errorsMeanSummedGamma[l]+ pow(errorsMeanGamma[i][l],2);
				errorsPosCorrSummedGamma[l] = errorsPosCorrSummedGamma[l]+pow(errorsPosCorrGamma[i][l],2);
				errorsNegCorrSummedGamma[l] = errorsNegCorrSummedGamma[l] +pow(errorsNegCorrGamma[i][l],2);
				errorsMeanCorrSummedGamma[l] =errorsMeanCorrSummedGamma[l]+ pow(errorsMeanCorrGamma[i][l],2);
				if(type == 1){
					errorsPosSummedGammaA[l] = errorsPosSummedGammaA[l]+pow(errorsPosGamma[i][l],2);
					errorsNegSummedGammaA[l] = errorsNegSummedGammaA[l]+ pow(errorsNegGamma[i][l],2);
					errorsMeanSummedGammaA[l] = errorsMeanSummedGammaA[l]+ pow(errorsMeanGamma[i][l],2);
					errorsPosCorrSummedGammaA[l] = errorsPosCorrSummedGammaA[l]+pow(errorsPosCorrGamma[i][l],2);
					errorsNegCorrSummedGammaA[l] = errorsNegCorrSummedGammaA[l] +pow(errorsNegCorrGamma[i][l],2);
					errorsMeanCorrSummedGammaA[l] =errorsMeanCorrSummedGammaA[l]+ pow(errorsMeanCorrGamma[i][l],2);
				}
				if(type == 2){
					errorsPosSummedGammaB[l] = errorsPosSummedGammaB[l]+pow(errorsPosGamma[i][l],2);
					errorsNegSummedGammaB[l] = errorsNegSummedGammaB[l]+ pow(errorsNegGamma[i][l],2);
					errorsMeanSummedGammaB[l] = errorsMeanSummedGammaB[l]+ pow(errorsMeanGamma[i][l],2);
					errorsPosCorrSummedGammaB[l] = errorsPosCorrSummedGammaB[l]+pow(errorsPosCorrGamma[i][l],2);
					errorsNegCorrSummedGammaB[l] = errorsNegCorrSummedGammaB[l] +pow(errorsNegCorrGamma[i][l],2);
					errorsMeanCorrSummedGammaB[l] =errorsMeanCorrSummedGammaB[l]+ pow(errorsMeanCorrGamma[i][l],2);
				}
				if(type == 3){
					errorsPosSummedGammaC[l] = errorsPosSummedGammaC[l]+pow(errorsPosGamma[i][l],2);
					errorsNegSummedGammaC[l] = errorsNegSummedGammaC[l]+ pow(errorsNegGamma[i][l],2);
					errorsMeanSummedGammaC[l] = errorsMeanSummedGammaC[l]+ pow(errorsMeanGamma[i][l],2);
					errorsPosCorrSummedGammaC[l] = errorsPosCorrSummedGammaC[l]+pow(errorsPosCorrGamma[i][l],2);
					errorsNegCorrSummedGammaC[l] = errorsNegCorrSummedGammaC[l] +pow(errorsNegCorrGamma[i][l],2);
					errorsMeanCorrSummedGammaC[l] =errorsMeanCorrSummedGammaC[l]+ pow(errorsMeanCorrGamma[i][l],2);
				}
			}
		}

		for (Int_t l = 0; l < ConstnBinsDoubleRatio; l++){

			errorsPosSummedDoubleRatio[l] = errorsPosSummedDoubleRatio[l]+pow(errorsPosDoubleRatio[i][l],2);
			errorsNegSummedDoubleRatio[l] = errorsNegSummedDoubleRatio[l]+ pow(errorsNegDoubleRatio[i][l],2);
			errorsMeanSummedDoubleRatio[l] = errorsMeanSummedDoubleRatio[l]+ pow(errorsMeanDoubleRatio[i][l],2);
			errorsPosCorrSummedDoubleRatio[l] = errorsPosCorrSummedDoubleRatio[l]+pow(errorsPosCorrDoubleRatio[i][l],2);
			errorsNegCorrSummedDoubleRatio[l] = errorsNegCorrSummedDoubleRatio[l] +pow(errorsNegCorrDoubleRatio[i][l],2);
			errorsMeanCorrSummedDoubleRatio[l] =errorsMeanCorrSummedDoubleRatio[l]+ pow(errorsMeanCorrDoubleRatio[i][l],2);

			errorsPosSummedDoubleRatioPi0Fit[l] = errorsPosSummedDoubleRatioPi0Fit[l]+pow(errorsPosDoubleRatioPi0Fit[i][l],2);
			errorsNegSummedDoubleRatioPi0Fit[l] = errorsNegSummedDoubleRatioPi0Fit[l]+ pow(errorsNegDoubleRatioPi0Fit[i][l],2);
			errorsMeanSummedDoubleRatioPi0Fit[l] = errorsMeanSummedDoubleRatioPi0Fit[l]+ pow(errorsMeanDoubleRatioPi0Fit[i][l],2);
			errorsPosCorrSummedDoubleRatioPi0Fit[l] = errorsPosCorrSummedDoubleRatioPi0Fit[l]+pow(errorsPosCorrDoubleRatioPi0Fit[i][l],2);
			errorsNegCorrSummedDoubleRatioPi0Fit[l] = errorsNegCorrSummedDoubleRatioPi0Fit[l] +pow(errorsNegCorrDoubleRatioPi0Fit[i][l],2);
			errorsMeanCorrSummedDoubleRatioPi0Fit[l] = errorsMeanCorrSummedDoubleRatioPi0Fit[l]+ pow(errorsMeanCorrDoubleRatioPi0Fit[i][l],2);

			if(type == 1){
				errorsPosSummedDoubleRatioA[l] = errorsPosSummedDoubleRatioA[l]+pow(errorsPosDoubleRatio[i][l],2);
				errorsNegSummedDoubleRatioA[l] = errorsNegSummedDoubleRatioA[l]+ pow(errorsNegDoubleRatio[i][l],2);
				errorsMeanSummedDoubleRatioA[l] = errorsMeanSummedDoubleRatioA[l]+ pow(errorsMeanDoubleRatio[i][l],2);
				errorsPosCorrSummedDoubleRatioA[l] = errorsPosCorrSummedDoubleRatioA[l]+pow(errorsPosCorrDoubleRatio[i][l],2);
				errorsNegCorrSummedDoubleRatioA[l] = errorsNegCorrSummedDoubleRatioA[l] +pow(errorsNegCorrDoubleRatio[i][l],2);
				errorsMeanCorrSummedDoubleRatioA[l] =errorsMeanCorrSummedDoubleRatioA[l]+ pow(errorsMeanCorrDoubleRatio[i][l],2);

				errorsPosSummedDoubleRatioPi0FitA[l] = errorsPosSummedDoubleRatioPi0FitA[l]+pow(errorsPosDoubleRatioPi0Fit[i][l],2);
				errorsNegSummedDoubleRatioPi0FitA[l] = errorsNegSummedDoubleRatioPi0FitA[l]+ pow(errorsNegDoubleRatioPi0Fit[i][l],2);
				errorsMeanSummedDoubleRatioPi0FitA[l] = errorsMeanSummedDoubleRatioPi0FitA[l]+ pow(errorsMeanDoubleRatioPi0Fit[i][l],2);
				errorsPosCorrSummedDoubleRatioPi0FitA[l] = errorsPosCorrSummedDoubleRatioPi0FitA[l]+pow(errorsPosCorrDoubleRatioPi0Fit[i][l],2);
				errorsNegCorrSummedDoubleRatioPi0FitA[l] = errorsNegCorrSummedDoubleRatioPi0FitA[l] +pow(errorsNegCorrDoubleRatioPi0Fit[i][l],2);
				errorsMeanCorrSummedDoubleRatioPi0FitA[l] = errorsMeanCorrSummedDoubleRatioPi0FitA[l]+ pow(errorsMeanCorrDoubleRatioPi0Fit[i][l],2);
			}
			if(type == 2){
				errorsPosSummedDoubleRatioB[l] = errorsPosSummedDoubleRatioB[l]+pow(errorsPosDoubleRatio[i][l],2);
				errorsNegSummedDoubleRatioB[l] = errorsNegSummedDoubleRatioB[l]+ pow(errorsNegDoubleRatio[i][l],2);
				errorsMeanSummedDoubleRatioB[l] = errorsMeanSummedDoubleRatioB[l]+ pow(errorsMeanDoubleRatio[i][l],2);
				errorsPosCorrSummedDoubleRatioB[l] = errorsPosCorrSummedDoubleRatioB[l]+pow(errorsPosCorrDoubleRatio[i][l],2);
				errorsNegCorrSummedDoubleRatioB[l] = errorsNegCorrSummedDoubleRatioB[l] +pow(errorsNegCorrDoubleRatio[i][l],2);
				errorsMeanCorrSummedDoubleRatioB[l] =errorsMeanCorrSummedDoubleRatioB[l]+ pow(errorsMeanCorrDoubleRatio[i][l],2);

				errorsPosSummedDoubleRatioPi0FitB[l] = errorsPosSummedDoubleRatioPi0FitB[l]+pow(errorsPosDoubleRatioPi0Fit[i][l],2);
				errorsNegSummedDoubleRatioPi0FitB[l] = errorsNegSummedDoubleRatioPi0FitB[l]+ pow(errorsNegDoubleRatioPi0Fit[i][l],2);
				errorsMeanSummedDoubleRatioPi0FitB[l] = errorsMeanSummedDoubleRatioPi0FitB[l]+ pow(errorsMeanDoubleRatioPi0Fit[i][l],2);
				errorsPosCorrSummedDoubleRatioPi0FitB[l] = errorsPosCorrSummedDoubleRatioPi0FitB[l]+pow(errorsPosCorrDoubleRatioPi0Fit[i][l],2);
				errorsNegCorrSummedDoubleRatioPi0FitB[l] = errorsNegCorrSummedDoubleRatioPi0FitB[l] +pow(errorsNegCorrDoubleRatioPi0Fit[i][l],2);
				errorsMeanCorrSummedDoubleRatioPi0FitB[l] = errorsMeanCorrSummedDoubleRatioPi0FitB[l]+ pow(errorsMeanCorrDoubleRatioPi0Fit[i][l],2);
			}
			if(type == 3){
				errorsPosSummedDoubleRatioC[l] = errorsPosSummedDoubleRatioC[l]+pow(errorsPosDoubleRatio[i][l],2);
				errorsNegSummedDoubleRatioC[l] = errorsNegSummedDoubleRatioC[l]+ pow(errorsNegDoubleRatio[i][l],2);
				errorsMeanSummedDoubleRatioC[l] = errorsMeanSummedDoubleRatioC[l]+ pow(errorsMeanDoubleRatio[i][l],2);
				errorsPosCorrSummedDoubleRatioC[l] = errorsPosCorrSummedDoubleRatioC[l]+pow(errorsPosCorrDoubleRatio[i][l],2);
				errorsNegCorrSummedDoubleRatioC[l] = errorsNegCorrSummedDoubleRatioC[l] +pow(errorsNegCorrDoubleRatio[i][l],2);
				errorsMeanCorrSummedDoubleRatioC[l] =errorsMeanCorrSummedDoubleRatioC[l]+ pow(errorsMeanCorrDoubleRatio[i][l],2);

				errorsPosSummedDoubleRatioPi0FitC[l] = errorsPosSummedDoubleRatioPi0FitC[l]+pow(errorsPosDoubleRatioPi0Fit[i][l],2);
				errorsNegSummedDoubleRatioPi0FitC[l] = errorsNegSummedDoubleRatioPi0FitC[l]+ pow(errorsNegDoubleRatioPi0Fit[i][l],2);
				errorsMeanSummedDoubleRatioPi0FitC[l] = errorsMeanSummedDoubleRatioPi0FitC[l]+ pow(errorsMeanDoubleRatioPi0Fit[i][l],2);
				errorsPosCorrSummedDoubleRatioPi0FitC[l] = errorsPosCorrSummedDoubleRatioPi0FitC[l]+pow(errorsPosCorrDoubleRatioPi0Fit[i][l],2);
				errorsNegCorrSummedDoubleRatioPi0FitC[l] = errorsNegCorrSummedDoubleRatioPi0FitC[l] +pow(errorsNegCorrDoubleRatioPi0Fit[i][l],2);
				errorsMeanCorrSummedDoubleRatioPi0FitC[l] = errorsMeanCorrSummedDoubleRatioPi0FitC[l]+ pow(errorsMeanCorrDoubleRatioPi0Fit[i][l],2);
			}
		}
	}


	for(Int_t l = 0; l <ConstnBinsDoubleRatio; l++){


		errorsPosSummedDoubleRatio[l] = pow(errorsPosSummedDoubleRatio[l],0.5);
		errorsMeanSummedDoubleRatio[l] = pow(errorsMeanSummedDoubleRatio[l],0.5);
		errorsPosErrSummedDoubleRatio[l] = errorsPosSummedDoubleRatio[l]*0.001;
		errorsMeanErrSummedDoubleRatio[l] = errorsMeanSummedDoubleRatio[l]*0.001;
		errorsNegSummedDoubleRatio[l] = -pow(errorsNegSummedDoubleRatio[l],0.5);
		errorsNegErrSummedDoubleRatio[l] = errorsNegSummedDoubleRatio[l]*0.001;

		errorsPosCorrMatSummedDoubleRatio[l] = pow(errorsPosCorrSummedDoubleRatio[l]+ pow(4.50 ,2.),0.5);
		errorsMeanCorrMatSummedDoubleRatio[l] = pow(errorsMeanCorrSummedDoubleRatio[l]+ pow(4.50 ,2.),0.5);
		errorsNegCorrMatSummedDoubleRatio[l] = -pow(errorsNegCorrSummedDoubleRatio[l]+ pow(4.50 ,2.),0.5);

		errorsPosCorrSummedDoubleRatio[l] = pow(errorsPosCorrSummedDoubleRatio[l],0.5);
		errorsMeanCorrSummedDoubleRatio[l] = pow(errorsMeanCorrSummedDoubleRatio[l],0.5);
		errorsPosErrCorrSummedDoubleRatio[l] = errorsPosCorrSummedDoubleRatio[l]*0.001;
		errorsMeanErrCorrSummedDoubleRatio[l] = errorsMeanCorrSummedDoubleRatio[l]*0.001;
		errorsMeanErrCorrMatSummedDoubleRatio[l] = errorsMeanCorrMatSummedDoubleRatio[l]*0.001;
		errorsNegCorrSummedDoubleRatio[l] = -pow(errorsNegCorrSummedDoubleRatio[l],0.5);
		errorsNegErrCorrSummedDoubleRatio[l] = errorsNegCorrSummedDoubleRatio[l]*0.001;

		errorsPosSummedDoubleRatioA[l] = pow(errorsPosSummedDoubleRatioA[l],0.5);
		errorsMeanSummedDoubleRatioA[l] = pow(errorsMeanSummedDoubleRatioA[l],0.5);
		errorsPosErrSummedDoubleRatioA[l] = errorsPosSummedDoubleRatioA[l]*0.001;
		errorsMeanErrSummedDoubleRatioA[l] = errorsMeanSummedDoubleRatioA[l]*0.001;
		errorsNegSummedDoubleRatioA[l] = -pow(errorsNegSummedDoubleRatioA[l],0.5);
		errorsNegErrSummedDoubleRatioA[l] = errorsNegSummedDoubleRatioA[l]*0.001;

		errorsPosCorrMatSummedDoubleRatioA[l] = pow(errorsPosCorrSummedDoubleRatioA[l],0.5);
		errorsMeanCorrMatSummedDoubleRatioA[l] = pow(errorsMeanCorrSummedDoubleRatioA[l],0.5);
		errorsNegCorrMatSummedDoubleRatioA[l] = -pow(errorsNegCorrSummedDoubleRatioA[l],0.5);

		errorsPosCorrSummedDoubleRatioA[l] = pow(errorsPosCorrSummedDoubleRatioA[l],0.5);
		errorsMeanCorrSummedDoubleRatioA[l] = pow(errorsMeanCorrSummedDoubleRatioA[l],0.5);
		errorsPosErrCorrSummedDoubleRatioA[l] = errorsPosCorrSummedDoubleRatioA[l]*0.001;
		errorsMeanErrCorrSummedDoubleRatioA[l] = errorsMeanCorrSummedDoubleRatioA[l]*0.001;
		errorsMeanErrCorrMatSummedDoubleRatioA[l] = errorsMeanCorrMatSummedDoubleRatioA[l]*0.001;
		errorsNegCorrSummedDoubleRatioA[l] = -pow(errorsNegCorrSummedDoubleRatioA[l],0.5);
		errorsNegErrCorrSummedDoubleRatioA[l] = errorsNegCorrSummedDoubleRatioA[l]*0.001;

		errorsPosSummedDoubleRatioB[l] = pow(errorsPosSummedDoubleRatioB[l],0.5);
		errorsMeanSummedDoubleRatioB[l] = pow(errorsMeanSummedDoubleRatioB[l],0.5);
		errorsPosErrSummedDoubleRatioB[l] = errorsPosSummedDoubleRatioB[l]*0.001;
		errorsMeanErrSummedDoubleRatioB[l] = errorsMeanSummedDoubleRatioB[l]*0.001;
		errorsNegSummedDoubleRatioB[l] = -pow(errorsNegSummedDoubleRatioB[l],0.5);
		errorsNegErrSummedDoubleRatioB[l] = errorsNegSummedDoubleRatioB[l]*0.001;

		errorsPosCorrMatSummedDoubleRatioB[l] = pow(errorsPosCorrSummedDoubleRatioB[l],0.5);
		errorsMeanCorrMatSummedDoubleRatioB[l] = pow(errorsMeanCorrSummedDoubleRatioB[l],0.5);
		errorsNegCorrMatSummedDoubleRatioB[l] = -pow(errorsNegCorrSummedDoubleRatioB[l],0.5);

		errorsPosCorrSummedDoubleRatioB[l] = pow(errorsPosCorrSummedDoubleRatioB[l],0.5);
		errorsMeanCorrSummedDoubleRatioB[l] = pow(errorsMeanCorrSummedDoubleRatioB[l],0.5);
		errorsPosErrCorrSummedDoubleRatioB[l] = errorsPosCorrSummedDoubleRatioB[l]*0.001;
		errorsMeanErrCorrSummedDoubleRatioB[l] = errorsMeanCorrSummedDoubleRatioB[l]*0.001;
		errorsMeanErrCorrMatSummedDoubleRatioB[l] = errorsMeanCorrMatSummedDoubleRatioB[l]*0.001;
		errorsNegCorrSummedDoubleRatioB[l] = -pow(errorsNegCorrSummedDoubleRatioB[l],0.5);
		errorsNegErrCorrSummedDoubleRatioB[l] = errorsNegCorrSummedDoubleRatioB[l]*0.001;

		errorsPosSummedDoubleRatioC[l] = pow(errorsPosSummedDoubleRatioC[l],0.5);
		errorsMeanSummedDoubleRatioC[l] = pow(errorsMeanSummedDoubleRatioC[l],0.5);
		errorsPosErrSummedDoubleRatioC[l] = errorsPosSummedDoubleRatioC[l]*0.001;
		errorsMeanErrSummedDoubleRatioC[l] = errorsMeanSummedDoubleRatioC[l]*0.001;
		errorsNegSummedDoubleRatioC[l] = -pow(errorsNegSummedDoubleRatioC[l],0.5);
		errorsNegErrSummedDoubleRatioC[l] = errorsNegSummedDoubleRatioC[l]*0.001;

		errorsPosCorrMatSummedDoubleRatioC[l] = pow(errorsPosCorrSummedDoubleRatioC[l]+ pow(4.50 ,2.),0.5);
		errorsMeanCorrMatSummedDoubleRatioC[l] = pow(errorsMeanCorrSummedDoubleRatioC[l]+ pow(4.50 ,2.),0.5);
		errorsNegCorrMatSummedDoubleRatioC[l] = -pow(errorsNegCorrSummedDoubleRatioC[l]+ pow(4.50 ,2.),0.5);

		errorsPosCorrSummedDoubleRatioC[l] = pow(errorsPosCorrSummedDoubleRatioC[l],0.5);
		errorsMeanCorrSummedDoubleRatioC[l] = pow(errorsMeanCorrSummedDoubleRatioC[l],0.5);
		errorsPosErrCorrSummedDoubleRatioC[l] = errorsPosCorrSummedDoubleRatioC[l]*0.001;
		errorsMeanErrCorrSummedDoubleRatioC[l] = errorsMeanCorrSummedDoubleRatioC[l]*0.001;
		errorsMeanErrCorrMatSummedDoubleRatioC[l] = errorsMeanCorrMatSummedDoubleRatioC[l]*0.001;
		errorsNegCorrSummedDoubleRatioC[l] = -pow(errorsNegCorrSummedDoubleRatioC[l],0.5);
		errorsNegErrCorrSummedDoubleRatioC[l] = errorsNegCorrSummedDoubleRatioC[l]*0.001;


		errorsPosSummedDoubleRatioPi0Fit[l] = pow(errorsPosSummedDoubleRatioPi0Fit[l],0.5);
		errorsMeanSummedDoubleRatioPi0Fit[l] = pow(errorsMeanSummedDoubleRatioPi0Fit[l],0.5);
		errorsPosErrSummedDoubleRatioPi0Fit[l] = errorsPosSummedDoubleRatioPi0Fit[l]*0.001;
		errorsMeanErrSummedDoubleRatioPi0Fit[l] = errorsMeanSummedDoubleRatioPi0Fit[l]*0.001;
		errorsNegSummedDoubleRatioPi0Fit[l] = -pow(errorsNegSummedDoubleRatioPi0Fit[l],0.5);
		errorsNegErrSummedDoubleRatioPi0Fit[l] = errorsNegSummedDoubleRatioPi0Fit[l]*0.001;

		errorsPosCorrMatSummedDoubleRatioPi0Fit[l] = pow(errorsPosCorrSummedDoubleRatioPi0Fit[l]+ pow(4.50 ,2.),0.5);
		errorsMeanCorrMatSummedDoubleRatioPi0Fit[l] = pow(errorsMeanCorrSummedDoubleRatioPi0Fit[l]+ pow(4.50 ,2.),0.5);
		errorsNegCorrMatSummedDoubleRatioPi0Fit[l] = -pow(errorsNegCorrSummedDoubleRatioPi0Fit[l]+ pow(4.50 ,2.),0.5);

		errorsPosCorrSummedDoubleRatioPi0Fit[l] = pow(errorsPosCorrSummedDoubleRatioPi0Fit[l],0.5);
		errorsMeanCorrSummedDoubleRatioPi0Fit[l] = pow(errorsMeanCorrSummedDoubleRatioPi0Fit[l],0.5);
		errorsPosErrCorrSummedDoubleRatioPi0Fit[l] = errorsPosCorrSummedDoubleRatioPi0Fit[l]*0.001;
		errorsMeanErrCorrSummedDoubleRatioPi0Fit[l] = errorsMeanCorrSummedDoubleRatioPi0Fit[l]*0.001;
		errorsMeanErrCorrMatSummedDoubleRatioPi0Fit[l] = errorsMeanCorrMatSummedDoubleRatioPi0Fit[l]*0.001;
		errorsNegCorrSummedDoubleRatioPi0Fit[l] = -pow(errorsNegCorrSummedDoubleRatioPi0Fit[l],0.5);
		errorsNegErrCorrSummedDoubleRatioPi0Fit[l] = errorsNegCorrSummedDoubleRatioPi0Fit[l]*0.001;

		errorsMeanErrCorrMatSummedDoubleRatioPi0FitWithoutMat[l] = errorsMeanCorrMatSummedDoubleRatioPi0FitWithoutMat[l]*0.001;
		errorsErrMaterialBudgetDoubleRatioPi0Fit[l] = errorsMaterialBudgetDoubleRatioPi0Fit[l]*0.001;


		errorsPosSummedDoubleRatioPi0FitA[l] = pow(errorsPosSummedDoubleRatioPi0FitA[l],0.5);
		errorsMeanSummedDoubleRatioPi0FitA[l] = pow(errorsMeanSummedDoubleRatioPi0FitA[l],0.5);
		errorsPosErrSummedDoubleRatioPi0FitA[l] = errorsPosSummedDoubleRatioPi0FitA[l]*0.001;
		errorsMeanErrSummedDoubleRatioPi0FitA[l] = errorsMeanSummedDoubleRatioPi0FitA[l]*0.001;
		errorsNegSummedDoubleRatioPi0FitA[l] = -pow(errorsNegSummedDoubleRatioPi0FitA[l],0.5);
		errorsNegErrSummedDoubleRatioPi0FitA[l] = errorsNegSummedDoubleRatioPi0FitA[l]*0.001;

		errorsPosCorrMatSummedDoubleRatioPi0FitA[l] = pow(errorsPosCorrSummedDoubleRatioPi0FitA[l],0.5);
		errorsMeanCorrMatSummedDoubleRatioPi0FitA[l] = pow(errorsMeanCorrSummedDoubleRatioPi0FitA[l],0.5);
		errorsNegCorrMatSummedDoubleRatioPi0FitA[l] = -pow(errorsNegCorrSummedDoubleRatioPi0FitA[l],0.5);

		errorsPosCorrSummedDoubleRatioPi0FitA[l] = pow(errorsPosCorrSummedDoubleRatioPi0FitA[l],0.5);
		errorsMeanCorrSummedDoubleRatioPi0FitA[l] = pow(errorsMeanCorrSummedDoubleRatioPi0FitA[l],0.5);
		errorsPosErrCorrSummedDoubleRatioPi0FitA[l] = errorsPosCorrSummedDoubleRatioPi0FitA[l]*0.001;
		errorsMeanErrCorrSummedDoubleRatioPi0FitA[l] = errorsMeanCorrSummedDoubleRatioPi0FitA[l]*0.001;
		errorsMeanErrCorrMatSummedDoubleRatioPi0FitA[l] = errorsMeanCorrMatSummedDoubleRatioPi0FitA[l]*0.001;
		errorsNegCorrSummedDoubleRatioPi0FitA[l] = -pow(errorsNegCorrSummedDoubleRatioPi0FitA[l],0.5);
		errorsNegErrCorrSummedDoubleRatioPi0FitA[l] = errorsNegCorrSummedDoubleRatioPi0FitA[l]*0.001;

		errorsPosSummedDoubleRatioPi0FitB[l] = pow(errorsPosSummedDoubleRatioPi0FitB[l],0.5);
		errorsMeanSummedDoubleRatioPi0FitB[l] = pow(errorsMeanSummedDoubleRatioPi0FitB[l],0.5);
		errorsPosErrSummedDoubleRatioPi0FitB[l] = errorsPosSummedDoubleRatioPi0FitB[l]*0.001;
		errorsMeanErrSummedDoubleRatioPi0FitB[l] = errorsMeanSummedDoubleRatioPi0FitB[l]*0.001;
		errorsNegSummedDoubleRatioPi0FitB[l] = -pow(errorsNegSummedDoubleRatioPi0FitB[l],0.5);
		errorsNegErrSummedDoubleRatioPi0FitB[l] = errorsNegSummedDoubleRatioPi0FitB[l]*0.001;

		errorsPosCorrMatSummedDoubleRatioPi0FitB[l] = pow(errorsPosCorrSummedDoubleRatioPi0FitB[l],0.5);
		errorsMeanCorrMatSummedDoubleRatioPi0FitB[l] = pow(errorsMeanCorrSummedDoubleRatioPi0FitB[l],0.5);
		errorsNegCorrMatSummedDoubleRatioPi0FitB[l] = -pow(errorsNegCorrSummedDoubleRatioPi0FitB[l],0.5);

		errorsPosCorrSummedDoubleRatioPi0FitB[l] = pow(errorsPosCorrSummedDoubleRatioPi0FitB[l],0.5);
		errorsMeanCorrSummedDoubleRatioPi0FitB[l] = pow(errorsMeanCorrSummedDoubleRatioPi0FitB[l],0.5);
		errorsPosErrCorrSummedDoubleRatioPi0FitB[l] = errorsPosCorrSummedDoubleRatioPi0FitB[l]*0.001;
		errorsMeanErrCorrSummedDoubleRatioPi0FitB[l] = errorsMeanCorrSummedDoubleRatioPi0FitB[l]*0.001;
		errorsMeanErrCorrMatSummedDoubleRatioPi0FitB[l] = errorsMeanCorrMatSummedDoubleRatioPi0FitB[l]*0.001;
		errorsNegCorrSummedDoubleRatioPi0FitB[l] = -pow(errorsNegCorrSummedDoubleRatioPi0FitB[l],0.5);
		errorsNegErrCorrSummedDoubleRatioPi0FitB[l] = errorsNegCorrSummedDoubleRatioPi0FitB[l]*0.001;

		errorsPosSummedDoubleRatioPi0FitC[l] = pow(errorsPosSummedDoubleRatioPi0FitC[l],0.5);
		errorsMeanSummedDoubleRatioPi0FitC[l] = pow(errorsMeanSummedDoubleRatioPi0FitC[l],0.5);
		errorsPosErrSummedDoubleRatioPi0FitC[l] = errorsPosSummedDoubleRatioPi0FitC[l]*0.001;
		errorsMeanErrSummedDoubleRatioPi0FitC[l] = errorsMeanSummedDoubleRatioPi0FitC[l]*0.001;
		errorsNegSummedDoubleRatioPi0FitC[l] = -pow(errorsNegSummedDoubleRatioPi0FitC[l],0.5);
		errorsNegErrSummedDoubleRatioPi0FitC[l] = errorsNegSummedDoubleRatioPi0FitC[l]*0.001;

		errorsPosCorrMatSummedDoubleRatioPi0FitC[l] = pow(errorsPosCorrSummedDoubleRatioPi0FitC[l]+ pow(4.50 ,2.),0.5);
		errorsMeanCorrMatSummedDoubleRatioPi0FitC[l] = pow(errorsMeanCorrSummedDoubleRatioPi0FitC[l]+ pow(4.50 ,2.),0.5);
		errorsNegCorrMatSummedDoubleRatioPi0FitC[l] = -pow(errorsNegCorrSummedDoubleRatioPi0FitC[l]+ pow(4.50 ,2.),0.5);

		errorsPosCorrSummedDoubleRatioPi0FitC[l] = pow(errorsPosCorrSummedDoubleRatioPi0FitC[l],0.5);
		errorsMeanCorrSummedDoubleRatioPi0FitC[l] = pow(errorsMeanCorrSummedDoubleRatioPi0FitC[l],0.5);
		errorsPosErrCorrSummedDoubleRatioPi0FitC[l] = errorsPosCorrSummedDoubleRatioPi0FitC[l]*0.001;
		errorsMeanErrCorrSummedDoubleRatioPi0FitC[l] = errorsMeanCorrSummedDoubleRatioPi0FitC[l]*0.001;
		errorsMeanErrCorrMatSummedDoubleRatioPi0FitC[l] = errorsMeanCorrMatSummedDoubleRatioPi0FitC[l]*0.001;
		errorsNegCorrSummedDoubleRatioPi0FitC[l] = -pow(errorsNegCorrSummedDoubleRatioPi0FitC[l],0.5);
		errorsNegErrCorrSummedDoubleRatioPi0FitC[l] = errorsNegCorrSummedDoubleRatioPi0FitC[l]*0.001;

	}

	for (Int_t l = 0; l <ConstnBinsIncRatio; l++){

		errorsPosSummedIncRatio[l] = pow(errorsPosSummedIncRatio[l],0.5);
		errorsMeanSummedIncRatio[l] = pow(errorsMeanSummedIncRatio[l],0.5);
		errorsPosErrSummedIncRatio[l] = errorsPosSummedIncRatio[l]*0.001;
		errorsMeanErrSummedIncRatio[l] = errorsMeanSummedIncRatio[l]*0.001;
		errorsNegSummedIncRatio[l] = -pow(errorsNegSummedIncRatio[l],0.5);
		errorsNegErrSummedIncRatio[l] = errorsNegSummedIncRatio[l]*0.001;
		errorsPosCorrMatSummedIncRatio[l] = pow(errorsPosCorrSummedIncRatio[l]+ pow(4.50 ,2.),0.5);
		errorsMeanCorrMatSummedIncRatio[l] = pow(errorsMeanCorrSummedIncRatio[l]+ pow(4.50 ,2.),0.5);
		errorsNegCorrMatSummedIncRatio[l] = -pow(errorsNegCorrSummedIncRatio[l]+ pow(4.50 ,2.),0.5);

		errorsPosCorrSummedIncRatio[l] = pow(errorsPosCorrSummedIncRatio[l],0.5);
		errorsMeanCorrSummedIncRatio[l] = pow(errorsMeanCorrSummedIncRatio[l],0.5);
		errorsPosErrCorrSummedIncRatio[l] = errorsPosCorrSummedIncRatio[l]*0.001;
		errorsMeanErrCorrSummedIncRatio[l] = errorsMeanCorrSummedIncRatio[l]*0.001;
		errorsMeanErrCorrMatSummedIncRatio[l] = errorsMeanCorrMatSummedIncRatio[l]*0.001;
		errorsNegCorrSummedIncRatio[l] = -pow(errorsNegCorrSummedIncRatio[l],0.5);
		errorsNegErrCorrSummedIncRatio[l] = errorsNegCorrSummedIncRatio[l]*0.001;


		errorsPosSummedIncRatioA[l] = pow(errorsPosSummedIncRatioA[l],0.5);
		errorsMeanSummedIncRatioA[l] = pow(errorsMeanSummedIncRatioA[l],0.5);
		errorsPosErrSummedIncRatioA[l] = errorsPosSummedIncRatioA[l]*0.001;
		errorsMeanErrSummedIncRatioA[l] = errorsMeanSummedIncRatioA[l]*0.001;
		errorsNegSummedIncRatioA[l] = -pow(errorsNegSummedIncRatioA[l],0.5);
		errorsNegErrSummedIncRatioA[l] = errorsNegSummedIncRatioA[l]*0.001;
		errorsPosCorrMatSummedIncRatioA[l] = pow(errorsPosCorrSummedIncRatioA[l],0.5);
		errorsMeanCorrMatSummedIncRatioA[l] = pow(errorsMeanCorrSummedIncRatioA[l],0.5);
		errorsNegCorrMatSummedIncRatioA[l] = -pow(errorsNegCorrSummedIncRatioA[l],0.5);

		errorsPosCorrSummedIncRatioA[l] = pow(errorsPosCorrSummedIncRatioA[l],0.5);
		errorsMeanCorrSummedIncRatioA[l] = pow(errorsMeanCorrSummedIncRatioA[l],0.5);
		errorsPosErrCorrSummedIncRatioA[l] = errorsPosCorrSummedIncRatioA[l]*0.001;
		errorsMeanErrCorrSummedIncRatioA[l] = errorsMeanCorrSummedIncRatioA[l]*0.001;
		errorsMeanErrCorrMatSummedIncRatioA[l] = errorsMeanCorrMatSummedIncRatioA[l]*0.001;
		errorsNegCorrSummedIncRatioA[l] = -pow(errorsNegCorrSummedIncRatioA[l],0.5);
		errorsNegErrCorrSummedIncRatioA[l] = errorsNegCorrSummedIncRatioA[l]*0.001;

		errorsPosSummedIncRatioB[l] = pow(errorsPosSummedIncRatioB[l],0.5);
		errorsMeanSummedIncRatioB[l] = pow(errorsMeanSummedIncRatioB[l],0.5);
		errorsPosErrSummedIncRatioB[l] = errorsPosSummedIncRatioB[l]*0.001;
		errorsMeanErrSummedIncRatioB[l] = errorsMeanSummedIncRatioB[l]*0.001;
		errorsNegSummedIncRatioB[l] = -pow(errorsNegSummedIncRatioB[l],0.5);
		errorsNegErrSummedIncRatioB[l] = errorsNegSummedIncRatioB[l]*0.001;
		errorsPosCorrMatSummedIncRatioB[l] = pow(errorsPosCorrSummedIncRatioB[l],0.5);
		errorsMeanCorrMatSummedIncRatioB[l] = pow(errorsMeanCorrSummedIncRatioB[l],0.5);
		errorsNegCorrMatSummedIncRatioB[l] = -pow(errorsNegCorrSummedIncRatioB[l],0.5);

		errorsPosCorrSummedIncRatioB[l] = pow(errorsPosCorrSummedIncRatioB[l],0.5);
		errorsMeanCorrSummedIncRatioB[l] = pow(errorsMeanCorrSummedIncRatioB[l],0.5);
		errorsPosErrCorrSummedIncRatioB[l] = errorsPosCorrSummedIncRatioB[l]*0.001;
		errorsMeanErrCorrSummedIncRatioB[l] = errorsMeanCorrSummedIncRatioB[l]*0.001;
		errorsMeanErrCorrMatSummedIncRatioB[l] = errorsMeanCorrMatSummedIncRatioB[l]*0.001;
		errorsNegCorrSummedIncRatioB[l] = -pow(errorsNegCorrSummedIncRatioB[l],0.5);
		errorsNegErrCorrSummedIncRatioB[l] = errorsNegCorrSummedIncRatioB[l]*0.001;

		errorsPosSummedIncRatioC[l] = pow(errorsPosSummedIncRatioC[l],0.5);
		errorsMeanSummedIncRatioC[l] = pow(errorsMeanSummedIncRatioC[l],0.5);
		errorsPosErrSummedIncRatioC[l] = errorsPosSummedIncRatioC[l]*0.001;
		errorsMeanErrSummedIncRatioC[l] = errorsMeanSummedIncRatioC[l]*0.001;
		errorsNegSummedIncRatioC[l] = -pow(errorsNegSummedIncRatioC[l],0.5);
		errorsNegErrSummedIncRatioC[l] = errorsNegSummedIncRatioC[l]*0.001;
		errorsPosCorrMatSummedIncRatioC[l] = pow(errorsPosCorrSummedIncRatioC[l]+ pow(4.50 ,2.),0.5);
		errorsMeanCorrMatSummedIncRatioC[l] = pow(errorsMeanCorrSummedIncRatioC[l]+ pow(4.50 ,2.),0.5);
		errorsNegCorrMatSummedIncRatioC[l] = -pow(errorsNegCorrSummedIncRatioC[l]+ pow(4.50 ,2.),0.5);

		errorsPosCorrSummedIncRatioC[l] = pow(errorsPosCorrSummedIncRatioC[l],0.5);
		errorsMeanCorrSummedIncRatioC[l] = pow(errorsMeanCorrSummedIncRatioC[l],0.5);
		errorsPosErrCorrSummedIncRatioC[l] = errorsPosCorrSummedIncRatioC[l]*0.001;
		errorsMeanErrCorrSummedIncRatioC[l] = errorsMeanCorrSummedIncRatioC[l]*0.001;
		errorsMeanErrCorrMatSummedIncRatioC[l] = errorsMeanCorrMatSummedIncRatioC[l]*0.001;
		errorsNegCorrSummedIncRatioC[l] = -pow(errorsNegCorrSummedIncRatioC[l],0.5);
		errorsNegErrCorrSummedIncRatioC[l] = errorsNegCorrSummedIncRatioC[l]*0.001;

		errorsPosSummedIncRatioPi0Fit[l] = pow(errorsPosSummedIncRatioPi0Fit[l],0.5);
		errorsMeanSummedIncRatioPi0Fit[l] = pow(errorsMeanSummedIncRatioPi0Fit[l],0.5);
		errorsPosErrSummedIncRatioPi0Fit[l] = errorsPosSummedIncRatioPi0Fit[l]*0.001;
		errorsMeanErrSummedIncRatioPi0Fit[l] = errorsMeanSummedIncRatioPi0Fit[l]*0.001;
		errorsNegSummedIncRatioPi0Fit[l] = -pow(errorsNegSummedIncRatioPi0Fit[l],0.5);
		errorsNegErrSummedIncRatioPi0Fit[l] = errorsNegSummedIncRatioPi0Fit[l]*0.001;
		errorsPosCorrMatSummedIncRatioPi0Fit[l] = pow(errorsPosCorrSummedIncRatioPi0Fit[l]+ pow(4.50 ,2.),0.5);
		errorsMeanCorrMatSummedIncRatioPi0Fit[l] = pow(errorsMeanCorrSummedIncRatioPi0Fit[l]+ pow(4.50 ,2.),0.5);
		errorsNegCorrMatSummedIncRatioPi0Fit[l] = -pow(errorsNegCorrSummedIncRatioPi0Fit[l]+ pow(4.50 ,2.),0.5);

		errorsPosCorrSummedIncRatioPi0Fit[l] = pow(errorsPosCorrSummedIncRatioPi0Fit[l],0.5);
		errorsMeanCorrSummedIncRatioPi0Fit[l] = pow(errorsMeanCorrSummedIncRatioPi0Fit[l],0.5);
		errorsPosErrCorrSummedIncRatioPi0Fit[l] = errorsPosCorrSummedIncRatioPi0Fit[l]*0.001;
		errorsMeanErrCorrSummedIncRatioPi0Fit[l] = errorsMeanCorrSummedIncRatioPi0Fit[l]*0.001;
		errorsMeanErrCorrMatSummedIncRatioPi0Fit[l] = errorsMeanCorrMatSummedIncRatioPi0Fit[l]*0.001;
		errorsNegCorrSummedIncRatioPi0Fit[l] = -pow(errorsNegCorrSummedIncRatioPi0Fit[l],0.5);
		errorsNegErrCorrSummedIncRatioPi0Fit[l] = errorsNegCorrSummedIncRatioPi0Fit[l]*0.001;


		errorsPosSummedIncRatioPi0FitA[l] = pow(errorsPosSummedIncRatioPi0FitA[l],0.5);
		errorsMeanSummedIncRatioPi0FitA[l] = pow(errorsMeanSummedIncRatioPi0FitA[l],0.5);
		errorsPosErrSummedIncRatioPi0FitA[l] = errorsPosSummedIncRatioPi0FitA[l]*0.001;
		errorsMeanErrSummedIncRatioPi0FitA[l] = errorsMeanSummedIncRatioPi0FitA[l]*0.001;
		errorsNegSummedIncRatioPi0FitA[l] = -pow(errorsNegSummedIncRatioPi0FitA[l],0.5);
		errorsNegErrSummedIncRatioPi0FitA[l] = errorsNegSummedIncRatioPi0FitA[l]*0.001;
		errorsPosCorrMatSummedIncRatioPi0FitA[l] = pow(errorsPosCorrSummedIncRatioPi0FitA[l],0.5);
		errorsMeanCorrMatSummedIncRatioPi0FitA[l] = pow(errorsMeanCorrSummedIncRatioPi0FitA[l],0.5);
		errorsNegCorrMatSummedIncRatioPi0FitA[l] = -pow(errorsNegCorrSummedIncRatioPi0FitA[l],0.5);

		errorsPosCorrSummedIncRatioPi0FitA[l] = pow(errorsPosCorrSummedIncRatioPi0FitA[l],0.5);
		errorsMeanCorrSummedIncRatioPi0FitA[l] = pow(errorsMeanCorrSummedIncRatioPi0FitA[l],0.5);
		errorsPosErrCorrSummedIncRatioPi0FitA[l] = errorsPosCorrSummedIncRatioPi0FitA[l]*0.001;
		errorsMeanErrCorrSummedIncRatioPi0FitA[l] = errorsMeanCorrSummedIncRatioPi0FitA[l]*0.001;
		errorsMeanErrCorrMatSummedIncRatioPi0FitA[l] = errorsMeanCorrMatSummedIncRatioPi0FitA[l]*0.001;
		errorsNegCorrSummedIncRatioPi0FitA[l] = -pow(errorsNegCorrSummedIncRatioPi0FitA[l],0.5);
		errorsNegErrCorrSummedIncRatioPi0FitA[l] = errorsNegCorrSummedIncRatioPi0FitA[l]*0.001;

		errorsPosSummedIncRatioPi0FitB[l] = pow(errorsPosSummedIncRatioPi0FitB[l],0.5);
		errorsMeanSummedIncRatioPi0FitB[l] = pow(errorsMeanSummedIncRatioPi0FitB[l],0.5);
		errorsPosErrSummedIncRatioPi0FitB[l] = errorsPosSummedIncRatioPi0FitB[l]*0.001;
		errorsMeanErrSummedIncRatioPi0FitB[l] = errorsMeanSummedIncRatioPi0FitB[l]*0.001;
		errorsNegSummedIncRatioPi0FitB[l] = -pow(errorsNegSummedIncRatioPi0FitB[l],0.5);
		errorsNegErrSummedIncRatioPi0FitB[l] = errorsNegSummedIncRatioPi0FitB[l]*0.001;
		errorsPosCorrMatSummedIncRatioPi0FitB[l] = pow(errorsPosCorrSummedIncRatioPi0FitB[l],0.5);
		errorsMeanCorrMatSummedIncRatioPi0FitB[l] = pow(errorsMeanCorrSummedIncRatioPi0FitB[l],0.5);
		errorsNegCorrMatSummedIncRatioPi0FitB[l] = -pow(errorsNegCorrSummedIncRatioPi0FitB[l],0.5);

		errorsPosCorrSummedIncRatioPi0FitB[l] = pow(errorsPosCorrSummedIncRatioPi0FitB[l],0.5);
		errorsMeanCorrSummedIncRatioPi0FitB[l] = pow(errorsMeanCorrSummedIncRatioPi0FitB[l],0.5);
		errorsPosErrCorrSummedIncRatioPi0FitB[l] = errorsPosCorrSummedIncRatioPi0FitB[l]*0.001;
		errorsMeanErrCorrSummedIncRatioPi0FitB[l] = errorsMeanCorrSummedIncRatioPi0FitB[l]*0.001;
		errorsMeanErrCorrMatSummedIncRatioPi0FitB[l] = errorsMeanCorrMatSummedIncRatioPi0FitB[l]*0.001;
		errorsNegCorrSummedIncRatioPi0FitB[l] = -pow(errorsNegCorrSummedIncRatioPi0FitB[l],0.5);
		errorsNegErrCorrSummedIncRatioPi0FitB[l] = errorsNegCorrSummedIncRatioPi0FitB[l]*0.001;

		errorsPosSummedIncRatioPi0FitC[l] = pow(errorsPosSummedIncRatioPi0FitC[l],0.5);
		errorsMeanSummedIncRatioPi0FitC[l] = pow(errorsMeanSummedIncRatioPi0FitC[l],0.5);
		errorsPosErrSummedIncRatioPi0FitC[l] = errorsPosSummedIncRatioPi0FitC[l]*0.001;
		errorsMeanErrSummedIncRatioPi0FitC[l] = errorsMeanSummedIncRatioPi0FitC[l]*0.001;
		errorsNegSummedIncRatioPi0FitC[l] = -pow(errorsNegSummedIncRatioPi0FitC[l],0.5);
		errorsNegErrSummedIncRatioPi0FitC[l] = errorsNegSummedIncRatioPi0FitC[l]*0.001;
		errorsPosCorrMatSummedIncRatioPi0FitC[l] = pow(errorsPosCorrSummedIncRatioPi0FitC[l]+ pow(4.50 ,2.),0.5);
		errorsMeanCorrMatSummedIncRatioPi0FitC[l] = pow(errorsMeanCorrSummedIncRatioPi0FitC[l]+ pow(4.50 ,2.),0.5);
		errorsNegCorrMatSummedIncRatioPi0FitC[l] = -pow(errorsNegCorrSummedIncRatioPi0FitC[l]+ pow(4.50 ,2.),0.5);

		errorsPosCorrSummedIncRatioPi0FitC[l] = pow(errorsPosCorrSummedIncRatioPi0FitC[l],0.5);
		errorsMeanCorrSummedIncRatioPi0FitC[l] = pow(errorsMeanCorrSummedIncRatioPi0FitC[l],0.5);
		errorsPosErrCorrSummedIncRatioPi0FitC[l] = errorsPosCorrSummedIncRatioPi0FitC[l]*0.001;
		errorsMeanErrCorrSummedIncRatioPi0FitC[l] = errorsMeanCorrSummedIncRatioPi0FitC[l]*0.001;
		errorsMeanErrCorrMatSummedIncRatioPi0FitC[l] = errorsMeanCorrMatSummedIncRatioPi0FitC[l]*0.001;
		errorsNegCorrSummedIncRatioPi0FitC[l] = -pow(errorsNegCorrSummedIncRatioPi0FitC[l],0.5);
		errorsNegErrCorrSummedIncRatioPi0FitC[l] = errorsNegCorrSummedIncRatioPi0FitC[l]*0.001;



	}

	for (Int_t l = 0; l <ConstnBinsPi0; l++){

		errorsPosSummedPi0[l] = pow(errorsPosSummedPi0[l],0.5);
		errorsMeanSummedPi0[l] = pow(errorsMeanSummedPi0[l],0.5);
		errorsPosErrSummedPi0[l] = errorsPosSummedPi0[l]*0.001;
		errorsMeanErrSummedPi0[l] = errorsMeanSummedPi0[l]*0.001;
		errorsNegSummedPi0[l] = -pow(errorsNegSummedPi0[l],0.5);
		errorsNegErrSummedPi0[l] = errorsNegSummedPi0[l]*0.001;
		errorsPosCorrMatSummedPi0[l] = pow(errorsPosCorrSummedPi0[l]+ pow(9.0 ,2.),0.5);
		errorsMeanCorrMatSummedPi0[l] = pow(errorsMeanCorrSummedPi0[l]+ pow(9.0 ,2.),0.5);
		errorsNegCorrMatSummedPi0[l] = -pow(errorsNegCorrSummedPi0[l]+ pow(9.0 ,2.),0.5);
		errorsMeanCorrMatSummedPi0WithoutMaterial[l] = pow(errorsMeanCorrSummedPi0[l],0.5);


		errorsPosCorrSummedPi0[l] = pow(errorsPosCorrSummedPi0[l],0.5);
		errorsMeanCorrSummedPi0[l] = pow(errorsMeanCorrSummedPi0[l],0.5);
		errorsPosErrCorrSummedPi0[l] = errorsPosCorrSummedPi0[l]*0.001;
		errorsMeanErrCorrSummedPi0[l] = errorsMeanCorrSummedPi0[l]*0.001;
		errorsMeanErrCorrMatSummedPi0[l] = errorsMeanCorrMatSummedPi0[l]*0.001;
		errorsNegCorrSummedPi0[l] = -pow(errorsNegCorrSummedPi0[l],0.5);
		errorsNegErrCorrSummedPi0[l] = errorsNegCorrSummedPi0[l]*0.001;

		errorsPosSummedPi0Fit[l] = pow(errorsPosSummedPi0Fit[l],0.5);
		errorsMeanSummedPi0Fit[l] = pow(errorsMeanSummedPi0Fit[l],0.5);
		errorsPosErrSummedPi0Fit[l] = errorsPosSummedPi0Fit[l]*0.001;
		errorsMeanErrSummedPi0Fit[l] = errorsMeanSummedPi0Fit[l]*0.001;
		errorsNegSummedPi0Fit[l] = -pow(errorsNegSummedPi0Fit[l],0.5);
		errorsNegErrSummedPi0Fit[l] = errorsNegSummedPi0Fit[l]*0.001;
		errorsPosCorrMatSummedPi0Fit[l] = pow(errorsPosCorrSummedPi0Fit[l]+ pow(9.0 ,2.),0.5);
		errorsMeanCorrMatSummedPi0Fit[l] = pow(errorsMeanCorrSummedPi0Fit[l]+ pow(9.0 ,2.),0.5);
		errorsNegCorrMatSummedPi0Fit[l] = -pow(errorsNegCorrSummedPi0Fit[l]+ pow(9.0 ,2.),0.5);
		errorsMeanCorrMatSummedPi0FitWithoutMaterial[l] = pow(errorsMeanCorrSummedPi0Fit[l],0.5);


		errorsPosCorrSummedPi0Fit[l] = pow(errorsPosCorrSummedPi0Fit[l],0.5);
		errorsMeanCorrSummedPi0Fit[l] = pow(errorsMeanCorrSummedPi0Fit[l],0.5);
		errorsPosErrCorrSummedPi0Fit[l] = errorsPosCorrSummedPi0Fit[l]*0.001;
		errorsMeanErrCorrSummedPi0Fit[l] = errorsMeanCorrSummedPi0Fit[l]*0.001;
		errorsMeanErrCorrMatSummedPi0Fit[l] = errorsMeanCorrMatSummedPi0Fit[l]*0.001;
		errorsNegCorrSummedPi0Fit[l] = -pow(errorsNegCorrSummedPi0Fit[l],0.5);
		errorsNegErrCorrSummedPi0Fit[l] = errorsNegCorrSummedPi0Fit[l]*0.001;
	}


	for (Int_t l = 0; l <ConstnBinsGamma; l++){
		errorsPosSummedGamma[l] = pow(errorsPosSummedGamma[l],0.5);
		errorsMeanSummedGamma[l] = pow(errorsMeanSummedGamma[l],0.5);
		errorsPosErrSummedGamma[l] = errorsPosSummedGamma[l]*0.001;
		errorsMeanErrSummedGamma[l] = errorsMeanSummedGamma[l]*0.001;
		errorsNegSummedGamma[l] = -pow(errorsNegSummedGamma[l],0.5);
		errorsNegErrSummedGamma[l] = errorsNegSummedGamma[l]*0.001;
		errorsPosCorrMatSummedGamma[l] = pow(errorsPosCorrSummedGamma[l]+ pow(4.50 ,2.),0.5);
		errorsMeanCorrMatSummedGamma[l] = pow(errorsMeanCorrSummedGamma[l]+ pow(4.50 ,2.),0.5);
		errorsNegCorrMatSummedGamma[l] = -pow(errorsNegCorrSummedGamma[l]+ pow(4.50 ,2.),0.5);

		errorsMeanCorrMatSummedGammaWithoutMat[l] =  pow(errorsMeanCorrSummedGamma[l],0.5);
		errorsMaterialBudget[l] = 4.5;

		errorsPosCorrSummedGamma[l] = pow(errorsPosCorrSummedGamma[l],0.5);
		errorsMeanCorrSummedGamma[l] = pow(errorsMeanCorrSummedGamma[l],0.5);
		errorsPosErrCorrSummedGamma[l] = errorsPosCorrSummedGamma[l]*0.001;
		errorsMeanErrCorrSummedGamma[l] = errorsMeanCorrSummedGamma[l]*0.001;
		errorsMeanErrCorrMatSummedGamma[l] = errorsMeanCorrMatSummedGamma[l]*0.001;
		errorsNegCorrSummedGamma[l] = -pow(errorsNegCorrSummedGamma[l],0.5);
		errorsNegErrCorrSummedGamma[l] = errorsNegCorrSummedGamma[l]*0.001;

		errorsMeanErrCorrMatSummedGammaWithoutMat[l] = errorsMeanCorrMatSummedGammaWithoutMat[l]*0.001;
		errorsErrMaterialBudget[l] = errorsMaterialBudget[l]*0.001;

		errorsPosSummedGammaA[l] = pow(errorsPosSummedGammaA[l],0.5);
		errorsMeanSummedGammaA[l] = pow(errorsMeanSummedGammaA[l],0.5);
		errorsPosErrSummedGammaA[l] = errorsPosSummedGammaA[l]*0.001;
		errorsMeanErrSummedGammaA[l] = errorsMeanSummedGammaA[l]*0.001;
		errorsNegSummedGammaA[l] = -pow(errorsNegSummedGammaA[l],0.5);
		errorsNegErrSummedGammaA[l] = errorsNegSummedGammaA[l]*0.001;
		errorsPosCorrMatSummedGammaA[l] = pow(errorsPosCorrSummedGammaA[l],0.5);
		errorsMeanCorrMatSummedGammaA[l] = pow(errorsMeanCorrSummedGammaA[l],0.5);
		errorsNegCorrMatSummedGammaA[l] = -pow(errorsNegCorrSummedGammaA[l],0.5);

		errorsPosCorrSummedGammaA[l] = pow(errorsPosCorrSummedGammaA[l],0.5);
		errorsMeanCorrSummedGammaA[l] = pow(errorsMeanCorrSummedGammaA[l],0.5);
		errorsPosErrCorrSummedGammaA[l] = errorsPosCorrSummedGammaA[l]*0.001;
		errorsMeanErrCorrSummedGammaA[l] = errorsMeanCorrSummedGammaA[l]*0.001;
		errorsMeanErrCorrMatSummedGammaA[l] = errorsMeanCorrMatSummedGammaA[l]*0.001;
		errorsNegCorrSummedGammaA[l] = -pow(errorsNegCorrSummedGammaA[l],0.5);
		errorsNegErrCorrSummedGammaA[l] = errorsNegCorrSummedGammaA[l]*0.001;

		errorsPosSummedGammaB[l] = pow(errorsPosSummedGammaB[l],0.5);
		errorsMeanSummedGammaB[l] = pow(errorsMeanSummedGammaB[l],0.5);
		errorsPosErrSummedGammaB[l] = errorsPosSummedGammaB[l]*0.001;
		errorsMeanErrSummedGammaB[l] = errorsMeanSummedGammaB[l]*0.001;
		errorsNegSummedGammaB[l] = -pow(errorsNegSummedGammaB[l],0.5);
		errorsNegErrSummedGammaB[l] = errorsNegSummedGammaB[l]*0.001;
		errorsPosCorrMatSummedGammaB[l] = pow(errorsPosCorrSummedGammaB[l],0.5);
		errorsMeanCorrMatSummedGammaB[l] = pow(errorsMeanCorrSummedGammaB[l],0.5);
		errorsNegCorrMatSummedGammaB[l] = -pow(errorsNegCorrSummedGammaB[l],0.5);

		errorsPosCorrSummedGammaB[l] = pow(errorsPosCorrSummedGammaB[l],0.5);
		errorsMeanCorrSummedGammaB[l] = pow(errorsMeanCorrSummedGammaB[l],0.5);
		errorsPosErrCorrSummedGammaB[l] = errorsPosCorrSummedGammaB[l]*0.001;
		errorsMeanErrCorrSummedGammaB[l] = errorsMeanCorrSummedGammaB[l]*0.001;
		errorsMeanErrCorrMatSummedGammaB[l] = errorsMeanCorrMatSummedGammaB[l]*0.001;
		errorsNegCorrSummedGammaB[l] = -pow(errorsNegCorrSummedGammaB[l],0.5);
		errorsNegErrCorrSummedGammaB[l] = errorsNegCorrSummedGammaB[l]*0.001;

		errorsPosSummedGammaC[l] = pow(errorsPosSummedGammaC[l],0.5);
		errorsMeanSummedGammaC[l] = pow(errorsMeanSummedGammaC[l],0.5);
		errorsPosErrSummedGammaC[l] = errorsPosSummedGammaC[l]*0.001;
		errorsMeanErrSummedGammaC[l] = errorsMeanSummedGammaC[l]*0.001;
		errorsNegSummedGammaC[l] = -pow(errorsNegSummedGammaC[l],0.5);
		errorsNegErrSummedGammaC[l] = errorsNegSummedGammaC[l]*0.001;
		errorsPosCorrMatSummedGammaC[l] = pow(errorsPosCorrSummedGammaC[l]+ pow(4.50 ,2.),0.5);
		errorsMeanCorrMatSummedGammaC[l] = pow(errorsMeanCorrSummedGammaC[l]+ pow(4.50 ,2.),0.5);
		errorsNegCorrMatSummedGammaC[l] = -pow(errorsNegCorrSummedGammaC[l]+ pow(4.50 ,2.),0.5);

		errorsPosCorrSummedGammaC[l] = pow(errorsPosCorrSummedGammaC[l],0.5);
		errorsMeanCorrSummedGammaC[l] = pow(errorsMeanCorrSummedGammaC[l],0.5);
		errorsPosErrCorrSummedGammaC[l] = errorsPosCorrSummedGammaC[l]*0.001;
		errorsMeanErrCorrSummedGammaC[l] = errorsMeanCorrSummedGammaC[l]*0.001;
		errorsMeanErrCorrMatSummedGammaC[l] = errorsMeanCorrMatSummedGammaC[l]*0.001;
		errorsNegCorrSummedGammaC[l] = -pow(errorsNegCorrSummedGammaC[l],0.5);
		errorsNegErrCorrSummedGammaC[l] = errorsNegCorrSummedGammaC[l]*0.001;
	}

	negativeErrorsSummedDoubleRatio = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegSummedDoubleRatio ,ptBinsDoubleRatioErr ,errorsNegErrSummedDoubleRatio );
	negativeErrorsCorrSummedDoubleRatio = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrSummedDoubleRatio ,ptBinsDoubleRatioErr ,errorsNegErrCorrSummedDoubleRatio );
	positiveErrorsSummedDoubleRatio = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosSummedDoubleRatio ,ptBinsDoubleRatioErr ,errorsPosErrSummedDoubleRatio );
	positiveErrorsCorrSummedDoubleRatio = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrSummedDoubleRatio ,ptBinsDoubleRatioErr ,errorsPosErrCorrSummedDoubleRatio );
	meanErrorsSummedDoubleRatio = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanSummedDoubleRatio ,ptBinsDoubleRatioErr ,errorsMeanErrSummedDoubleRatio );
	meanErrorsCorrSummedDoubleRatio = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrSummedDoubleRatio ,ptBinsDoubleRatioErr ,errorsMeanErrCorrSummedDoubleRatio );
	meanErrorsCorrSummedIncMatDoubleRatio = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrMatSummedDoubleRatio ,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedDoubleRatio );

	negativeErrorsSummedDoubleRatioPi0Fit = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegSummedDoubleRatioPi0Fit ,ptBinsDoubleRatioErr ,errorsNegErrSummedDoubleRatioPi0Fit );
	negativeErrorsCorrSummedDoubleRatioPi0Fit = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrSummedDoubleRatioPi0Fit ,ptBinsDoubleRatioErr ,errorsNegErrCorrSummedDoubleRatioPi0Fit );
	positiveErrorsSummedDoubleRatioPi0Fit = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosSummedDoubleRatioPi0Fit ,ptBinsDoubleRatioErr ,errorsPosErrSummedDoubleRatioPi0Fit );
	positiveErrorsCorrSummedDoubleRatioPi0Fit = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrSummedDoubleRatioPi0Fit ,ptBinsDoubleRatioErr ,errorsPosErrCorrSummedDoubleRatioPi0Fit );
	meanErrorsSummedDoubleRatioPi0Fit = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanSummedDoubleRatioPi0Fit ,ptBinsDoubleRatioErr ,errorsMeanErrSummedDoubleRatioPi0Fit );
	meanErrorsCorrSummedDoubleRatioPi0Fit = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrSummedDoubleRatioPi0Fit ,ptBinsDoubleRatioErr ,errorsMeanErrCorrSummedDoubleRatioPi0Fit );
	meanErrorsCorrSummedIncMatDoubleRatioPi0Fit = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrMatSummedDoubleRatioPi0Fit ,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedDoubleRatioPi0Fit );

	meanErrorsCorrSummedMaterialDoubleRatioPi0FitCocktail = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrMatSummedDoubleRatioPi0Fit ,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedDoubleRatioPi0Fit );
	for(Int_t i = 0;i<meanErrorsCorrSummedMaterialDoubleRatioPi0FitCocktail->GetN();i++){
		Double_t x,y;
		meanErrorsCorrSummedMaterialDoubleRatioPi0FitCocktail->GetPoint(i,x,y);
		meanErrorsCorrSummedMaterialDoubleRatioPi0FitCocktail->SetPoint(i,x,0);
	}

	negativeErrorsSummedDoubleRatioA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegSummedDoubleRatioA ,ptBinsDoubleRatioErr ,errorsNegErrSummedDoubleRatioA );
	negativeErrorsCorrSummedDoubleRatioA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrSummedDoubleRatioA ,ptBinsDoubleRatioErr ,errorsNegErrCorrSummedDoubleRatioA );
	positiveErrorsSummedDoubleRatioA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosSummedDoubleRatioA ,ptBinsDoubleRatioErr ,errorsPosErrSummedDoubleRatioA );
	positiveErrorsCorrSummedDoubleRatioA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrSummedDoubleRatioA ,ptBinsDoubleRatioErr ,errorsPosErrCorrSummedDoubleRatioA );
	meanErrorsSummedDoubleRatioA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanSummedDoubleRatioA ,ptBinsDoubleRatioErr ,errorsMeanErrSummedDoubleRatioA );
	meanErrorsCorrSummedDoubleRatioA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrSummedDoubleRatioA ,ptBinsDoubleRatioErr ,errorsMeanErrCorrSummedDoubleRatioA );
	meanErrorsCorrSummedIncMatDoubleRatioA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrMatSummedDoubleRatioA ,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedDoubleRatioA );

	negativeErrorsSummedDoubleRatioPi0FitA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegSummedDoubleRatioPi0FitA ,ptBinsDoubleRatioErr ,errorsNegErrSummedDoubleRatioPi0FitA );
	negativeErrorsCorrSummedDoubleRatioPi0FitA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrSummedDoubleRatioPi0FitA ,ptBinsDoubleRatioErr ,errorsNegErrCorrSummedDoubleRatioPi0FitA );
	positiveErrorsSummedDoubleRatioPi0FitA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosSummedDoubleRatioPi0FitA ,ptBinsDoubleRatioErr ,errorsPosErrSummedDoubleRatioPi0FitA );
	positiveErrorsCorrSummedDoubleRatioPi0FitA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrSummedDoubleRatioPi0FitA ,ptBinsDoubleRatioErr ,errorsPosErrCorrSummedDoubleRatioPi0FitA );
	meanErrorsSummedDoubleRatioPi0FitA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanSummedDoubleRatioPi0FitA ,ptBinsDoubleRatioErr ,errorsMeanErrSummedDoubleRatioPi0FitA );
	meanErrorsCorrSummedDoubleRatioPi0FitA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrSummedDoubleRatioPi0FitA ,ptBinsDoubleRatioErr ,errorsMeanErrCorrSummedDoubleRatioPi0FitA );
	meanErrorsCorrSummedIncMatDoubleRatioPi0FitA = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrMatSummedDoubleRatioPi0FitA ,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedDoubleRatioPi0FitA );

	negativeErrorsSummedDoubleRatioB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegSummedDoubleRatioB ,ptBinsDoubleRatioErr ,errorsNegErrSummedDoubleRatioB );
	negativeErrorsCorrSummedDoubleRatioB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrSummedDoubleRatioB ,ptBinsDoubleRatioErr ,errorsNegErrCorrSummedDoubleRatioB );
	positiveErrorsSummedDoubleRatioB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosSummedDoubleRatioB ,ptBinsDoubleRatioErr ,errorsPosErrSummedDoubleRatioB );
	positiveErrorsCorrSummedDoubleRatioB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrSummedDoubleRatioB ,ptBinsDoubleRatioErr ,errorsPosErrCorrSummedDoubleRatioB );
	meanErrorsSummedDoubleRatioB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanSummedDoubleRatioB ,ptBinsDoubleRatioErr ,errorsMeanErrSummedDoubleRatioB );
	meanErrorsCorrSummedDoubleRatioB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrSummedDoubleRatioB ,ptBinsDoubleRatioErr ,errorsMeanErrCorrSummedDoubleRatioB );
	meanErrorsCorrSummedIncMatDoubleRatioB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrMatSummedDoubleRatioB ,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedDoubleRatioB );

	negativeErrorsSummedDoubleRatioPi0FitB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegSummedDoubleRatioPi0FitB ,ptBinsDoubleRatioErr ,errorsNegErrSummedDoubleRatioPi0FitB );
	negativeErrorsCorrSummedDoubleRatioPi0FitB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrSummedDoubleRatioPi0FitB ,ptBinsDoubleRatioErr ,errorsNegErrCorrSummedDoubleRatioPi0FitB );
	positiveErrorsSummedDoubleRatioPi0FitB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosSummedDoubleRatioPi0FitB ,ptBinsDoubleRatioErr ,errorsPosErrSummedDoubleRatioPi0FitB );
	positiveErrorsCorrSummedDoubleRatioPi0FitB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrSummedDoubleRatioPi0FitB ,ptBinsDoubleRatioErr ,errorsPosErrCorrSummedDoubleRatioPi0FitB );
	meanErrorsSummedDoubleRatioPi0FitB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanSummedDoubleRatioPi0FitB ,ptBinsDoubleRatioErr ,errorsMeanErrSummedDoubleRatioPi0FitB );
	meanErrorsCorrSummedDoubleRatioPi0FitB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrSummedDoubleRatioPi0FitB ,ptBinsDoubleRatioErr ,errorsMeanErrCorrSummedDoubleRatioPi0FitB );
	meanErrorsCorrSummedIncMatDoubleRatioPi0FitB = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrMatSummedDoubleRatioPi0FitB ,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedDoubleRatioPi0FitB );

	negativeErrorsSummedDoubleRatioC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegSummedDoubleRatioC ,ptBinsDoubleRatioErr ,errorsNegErrSummedDoubleRatioC );
	negativeErrorsCorrSummedDoubleRatioC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrSummedDoubleRatioC ,ptBinsDoubleRatioErr ,errorsNegErrCorrSummedDoubleRatioC );
	positiveErrorsSummedDoubleRatioC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosSummedDoubleRatioC ,ptBinsDoubleRatioErr ,errorsPosErrSummedDoubleRatioC );
	positiveErrorsCorrSummedDoubleRatioC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrSummedDoubleRatioC ,ptBinsDoubleRatioErr ,errorsPosErrCorrSummedDoubleRatioC );
	meanErrorsSummedDoubleRatioC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanSummedDoubleRatioC ,ptBinsDoubleRatioErr ,errorsMeanErrSummedDoubleRatioC );
	meanErrorsCorrSummedDoubleRatioC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrSummedDoubleRatioC ,ptBinsDoubleRatioErr ,errorsMeanErrCorrSummedDoubleRatioC );
	meanErrorsCorrSummedIncMatDoubleRatioC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrMatSummedDoubleRatioC ,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedDoubleRatioC );

	negativeErrorsSummedDoubleRatioPi0FitC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegSummedDoubleRatioPi0FitC ,ptBinsDoubleRatioErr ,errorsNegErrSummedDoubleRatioPi0FitC );
	negativeErrorsCorrSummedDoubleRatioPi0FitC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsNegCorrSummedDoubleRatioPi0FitC ,ptBinsDoubleRatioErr ,errorsNegErrCorrSummedDoubleRatioPi0FitC );
	positiveErrorsSummedDoubleRatioPi0FitC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosSummedDoubleRatioPi0FitC ,ptBinsDoubleRatioErr ,errorsPosErrSummedDoubleRatioPi0FitC );
	positiveErrorsCorrSummedDoubleRatioPi0FitC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsPosCorrSummedDoubleRatioPi0FitC ,ptBinsDoubleRatioErr ,errorsPosErrCorrSummedDoubleRatioPi0FitC );
	meanErrorsSummedDoubleRatioPi0FitC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanSummedDoubleRatioPi0FitC ,ptBinsDoubleRatioErr ,errorsMeanErrSummedDoubleRatioPi0FitC );
	meanErrorsCorrSummedDoubleRatioPi0FitC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrSummedDoubleRatioPi0FitC ,ptBinsDoubleRatioErr ,errorsMeanErrCorrSummedDoubleRatioPi0FitC );
	meanErrorsCorrSummedIncMatDoubleRatioPi0FitC = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio ,errorsMeanCorrMatSummedDoubleRatioPi0FitC ,ptBinsDoubleRatioErr ,errorsMeanErrCorrMatSummedDoubleRatioPi0FitC );


	meanErrorsCorrSummedWithoutMatDoubleRatioPi0Fit = new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio , errorsMeanCorrMatSummedDoubleRatioPi0FitWithoutMat,ptBinsDoubleRatioErr , errorsMeanErrCorrMatSummedDoubleRatioPi0FitWithoutMat);
	meanErrorsCorrSummedMaterialDoubleRatioPi0Fit =   new TGraphErrors(ConstnBinsDoubleRatio,ptBinsDoubleRatio , errorsMaterialBudgetDoubleRatioPi0Fit,ptBinsDoubleRatioErr , errorsErrMaterialBudgetDoubleRatioPi0Fit);

	negativeErrorsSummedIncRatio = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegSummedIncRatio ,ptBinsIncRatioErr ,errorsNegErrSummedIncRatio );
	negativeErrorsCorrSummedIncRatio = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrSummedIncRatio ,ptBinsIncRatioErr ,errorsNegErrCorrSummedIncRatio );
	positiveErrorsSummedIncRatio = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosSummedIncRatio ,ptBinsIncRatioErr ,errorsPosErrSummedIncRatio );
	positiveErrorsCorrSummedIncRatio = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrSummedIncRatio ,ptBinsIncRatioErr ,errorsPosErrCorrSummedIncRatio );
	meanErrorsSummedIncRatio = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanSummedIncRatio ,ptBinsIncRatioErr ,errorsMeanErrSummedIncRatio );
	meanErrorsCorrSummedIncRatio = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrSummedIncRatio ,ptBinsIncRatioErr ,errorsMeanErrCorrSummedIncRatio );
	meanErrorsCorrSummedIncMatIncRatio = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrMatSummedIncRatio ,ptBinsIncRatioErr ,errorsMeanErrCorrMatSummedIncRatio );

	negativeErrorsSummedIncRatioPi0Fit = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegSummedIncRatioPi0Fit ,ptBinsIncRatioErr ,errorsNegErrSummedIncRatioPi0Fit );
	negativeErrorsCorrSummedIncRatioPi0Fit = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrSummedIncRatioPi0Fit ,ptBinsIncRatioErr ,errorsNegErrCorrSummedIncRatioPi0Fit );
	positiveErrorsSummedIncRatioPi0Fit = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosSummedIncRatioPi0Fit ,ptBinsIncRatioErr ,errorsPosErrSummedIncRatioPi0Fit );
	positiveErrorsCorrSummedIncRatioPi0Fit = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrSummedIncRatioPi0Fit ,ptBinsIncRatioErr ,errorsPosErrCorrSummedIncRatioPi0Fit );
	meanErrorsSummedIncRatioPi0Fit = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanSummedIncRatioPi0Fit ,ptBinsIncRatioErr ,errorsMeanErrSummedIncRatioPi0Fit );
	meanErrorsCorrSummedIncRatioPi0Fit = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrSummedIncRatioPi0Fit ,ptBinsIncRatioErr ,errorsMeanErrCorrSummedIncRatioPi0Fit );
	meanErrorsCorrSummedIncMatIncRatioPi0Fit = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrMatSummedIncRatioPi0Fit ,ptBinsIncRatioErr ,errorsMeanErrCorrMatSummedIncRatioPi0Fit );

	negativeErrorsSummedIncRatioA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegSummedIncRatioA ,ptBinsIncRatioErr ,errorsNegErrSummedIncRatioA );
	negativeErrorsCorrSummedIncRatioA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrSummedIncRatioA ,ptBinsIncRatioErr ,errorsNegErrCorrSummedIncRatioA );
	positiveErrorsSummedIncRatioA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosSummedIncRatioA ,ptBinsIncRatioErr ,errorsPosErrSummedIncRatioA );
	positiveErrorsCorrSummedIncRatioA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrSummedIncRatioA ,ptBinsIncRatioErr ,errorsPosErrCorrSummedIncRatioA );
	meanErrorsSummedIncRatioA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanSummedIncRatioA ,ptBinsIncRatioErr ,errorsMeanErrSummedIncRatioA );
	meanErrorsCorrSummedIncRatioA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrSummedIncRatioA ,ptBinsIncRatioErr ,errorsMeanErrCorrSummedIncRatioA );
	meanErrorsCorrSummedIncMatIncRatioA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrMatSummedIncRatioA ,ptBinsIncRatioErr ,errorsMeanErrCorrMatSummedIncRatioA );

	negativeErrorsSummedIncRatioPi0FitA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegSummedIncRatioPi0FitA ,ptBinsIncRatioErr ,errorsNegErrSummedIncRatioPi0FitA );
	negativeErrorsCorrSummedIncRatioPi0FitA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrSummedIncRatioPi0FitA ,ptBinsIncRatioErr ,errorsNegErrCorrSummedIncRatioPi0FitA );
	positiveErrorsSummedIncRatioPi0FitA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosSummedIncRatioPi0FitA ,ptBinsIncRatioErr ,errorsPosErrSummedIncRatioPi0FitA );
	positiveErrorsCorrSummedIncRatioPi0FitA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrSummedIncRatioPi0FitA ,ptBinsIncRatioErr ,errorsPosErrCorrSummedIncRatioPi0FitA );
	meanErrorsSummedIncRatioPi0FitA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanSummedIncRatioPi0FitA ,ptBinsIncRatioErr ,errorsMeanErrSummedIncRatioPi0FitA );
	meanErrorsCorrSummedIncRatioPi0FitA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrSummedIncRatioPi0FitA ,ptBinsIncRatioErr ,errorsMeanErrCorrSummedIncRatioPi0FitA );
	meanErrorsCorrSummedIncMatIncRatioPi0FitA = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrMatSummedIncRatioPi0FitA ,ptBinsIncRatioErr ,errorsMeanErrCorrMatSummedIncRatioPi0FitA );

	negativeErrorsSummedIncRatioB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegSummedIncRatioB ,ptBinsIncRatioErr ,errorsNegErrSummedIncRatioB );
	negativeErrorsCorrSummedIncRatioB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrSummedIncRatioB ,ptBinsIncRatioErr ,errorsNegErrCorrSummedIncRatioB );
	positiveErrorsSummedIncRatioB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosSummedIncRatioB ,ptBinsIncRatioErr ,errorsPosErrSummedIncRatioB );
	positiveErrorsCorrSummedIncRatioB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrSummedIncRatioB ,ptBinsIncRatioErr ,errorsPosErrCorrSummedIncRatioB );
	meanErrorsSummedIncRatioB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanSummedIncRatioB ,ptBinsIncRatioErr ,errorsMeanErrSummedIncRatioB );
	meanErrorsCorrSummedIncRatioB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrSummedIncRatioB ,ptBinsIncRatioErr ,errorsMeanErrCorrSummedIncRatioB );
	meanErrorsCorrSummedIncMatIncRatioB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrMatSummedIncRatioB ,ptBinsIncRatioErr ,errorsMeanErrCorrMatSummedIncRatioB );

	negativeErrorsSummedIncRatioPi0FitB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegSummedIncRatioPi0FitB ,ptBinsIncRatioErr ,errorsNegErrSummedIncRatioPi0FitB );
	negativeErrorsCorrSummedIncRatioPi0FitB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrSummedIncRatioPi0FitB ,ptBinsIncRatioErr ,errorsNegErrCorrSummedIncRatioPi0FitB );
	positiveErrorsSummedIncRatioPi0FitB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosSummedIncRatioPi0FitB ,ptBinsIncRatioErr ,errorsPosErrSummedIncRatioPi0FitB );
	positiveErrorsCorrSummedIncRatioPi0FitB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrSummedIncRatioPi0FitB ,ptBinsIncRatioErr ,errorsPosErrCorrSummedIncRatioPi0FitB );
	meanErrorsSummedIncRatioPi0FitB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanSummedIncRatioPi0FitB ,ptBinsIncRatioErr ,errorsMeanErrSummedIncRatioPi0FitB );
	meanErrorsCorrSummedIncRatioPi0FitB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrSummedIncRatioPi0FitB ,ptBinsIncRatioErr ,errorsMeanErrCorrSummedIncRatioPi0FitB );
	meanErrorsCorrSummedIncMatIncRatioPi0FitB = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrMatSummedIncRatioPi0FitB ,ptBinsIncRatioErr ,errorsMeanErrCorrMatSummedIncRatioPi0FitB );

	negativeErrorsSummedIncRatioC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegSummedIncRatioC ,ptBinsIncRatioErr ,errorsNegErrSummedIncRatioC );
	negativeErrorsCorrSummedIncRatioC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrSummedIncRatioC ,ptBinsIncRatioErr ,errorsNegErrCorrSummedIncRatioC );
	positiveErrorsSummedIncRatioC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosSummedIncRatioC ,ptBinsIncRatioErr ,errorsPosErrSummedIncRatioC );
	positiveErrorsCorrSummedIncRatioC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrSummedIncRatioC ,ptBinsIncRatioErr ,errorsPosErrCorrSummedIncRatioC );
	meanErrorsSummedIncRatioC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanSummedIncRatioC ,ptBinsIncRatioErr ,errorsMeanErrSummedIncRatioC );
	meanErrorsCorrSummedIncRatioC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrSummedIncRatioC ,ptBinsIncRatioErr ,errorsMeanErrCorrSummedIncRatioC );
	meanErrorsCorrSummedIncMatIncRatioC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrMatSummedIncRatioC ,ptBinsIncRatioErr ,errorsMeanErrCorrMatSummedIncRatioC );

	negativeErrorsSummedIncRatioPi0FitC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegSummedIncRatioPi0FitC ,ptBinsIncRatioErr ,errorsNegErrSummedIncRatioPi0FitC );
	negativeErrorsCorrSummedIncRatioPi0FitC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsNegCorrSummedIncRatioPi0FitC ,ptBinsIncRatioErr ,errorsNegErrCorrSummedIncRatioPi0FitC );
	positiveErrorsSummedIncRatioPi0FitC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosSummedIncRatioPi0FitC ,ptBinsIncRatioErr ,errorsPosErrSummedIncRatioPi0FitC );
	positiveErrorsCorrSummedIncRatioPi0FitC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsPosCorrSummedIncRatioPi0FitC ,ptBinsIncRatioErr ,errorsPosErrCorrSummedIncRatioPi0FitC );
	meanErrorsSummedIncRatioPi0FitC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanSummedIncRatioPi0FitC ,ptBinsIncRatioErr ,errorsMeanErrSummedIncRatioPi0FitC );
	meanErrorsCorrSummedIncRatioPi0FitC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrSummedIncRatioPi0FitC ,ptBinsIncRatioErr ,errorsMeanErrCorrSummedIncRatioPi0FitC );
	meanErrorsCorrSummedIncMatIncRatioPi0FitC = new TGraphErrors(ConstnBinsIncRatio,ptBinsIncRatio ,errorsMeanCorrMatSummedIncRatioPi0FitC ,ptBinsIncRatioErr ,errorsMeanErrCorrMatSummedIncRatioPi0FitC );

	negativeErrorsSummedPi0 = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsNegSummedPi0 ,ptBinsPi0Err ,errorsNegErrSummedPi0 );
	negativeErrorsCorrSummedPi0 = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsNegCorrSummedPi0 ,ptBinsPi0Err ,errorsNegErrCorrSummedPi0 );
	positiveErrorsSummedPi0 = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsPosSummedPi0 ,ptBinsPi0Err ,errorsPosErrSummedPi0 );
	positiveErrorsCorrSummedPi0 = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsPosCorrSummedPi0 ,ptBinsPi0Err ,errorsPosErrCorrSummedPi0 );
	meanErrorsSummedPi0 = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanSummedPi0 ,ptBinsPi0Err ,errorsMeanErrSummedPi0 );
	meanErrorsCorrSummedPi0 = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanCorrSummedPi0 ,ptBinsPi0Err ,errorsMeanErrCorrSummedPi0 );
	meanErrorsCorrSummedIncMatPi0 = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanCorrMatSummedPi0 ,ptBinsPi0Err ,errorsMeanErrCorrMatSummedPi0 );
	meanErrorsCorrSummedIncMatPi0WithoutMat = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanCorrMatSummedPi0WithoutMaterial ,ptBinsPi0Err ,errorsMeanErrCorrMatSummedPi0 );

	negativeErrorsSummedPi0Fit = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsNegSummedPi0Fit ,ptBinsPi0Err ,errorsNegErrSummedPi0Fit );
	negativeErrorsCorrSummedPi0Fit = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsNegCorrSummedPi0Fit ,ptBinsPi0Err ,errorsNegErrCorrSummedPi0Fit );
	positiveErrorsSummedPi0Fit = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsPosSummedPi0Fit ,ptBinsPi0Err ,errorsPosErrSummedPi0Fit );
	positiveErrorsCorrSummedPi0Fit = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsPosCorrSummedPi0Fit ,ptBinsPi0Err ,errorsPosErrCorrSummedPi0Fit );
	meanErrorsSummedPi0Fit = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanSummedPi0Fit ,ptBinsPi0Err ,errorsMeanErrSummedPi0Fit );
	meanErrorsCorrSummedPi0Fit = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanCorrSummedPi0Fit ,ptBinsPi0Err ,errorsMeanErrCorrSummedPi0Fit );
	meanErrorsCorrSummedIncMatPi0Fit = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanCorrMatSummedPi0Fit ,ptBinsPi0Err ,errorsMeanErrCorrMatSummedPi0Fit );
	meanErrorsCorrSummedIncMatPi0FitWithoutMat = new TGraphErrors(ConstnBinsPi0,ptBinsPi0 ,errorsMeanCorrMatSummedPi0FitWithoutMaterial ,ptBinsPi0Err ,errorsMeanErrCorrMatSummedPi0Fit );

	negativeErrorsSummedGamma = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegSummedGamma ,ptBinsGammaErr ,errorsNegErrSummedGamma );
	negativeErrorsCorrSummedGamma = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegCorrSummedGamma ,ptBinsGammaErr ,errorsNegErrCorrSummedGamma );
	positiveErrorsSummedGamma = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosSummedGamma ,ptBinsGammaErr ,errorsPosErrSummedGamma );
	positiveErrorsCorrSummedGamma = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosCorrSummedGamma ,ptBinsGammaErr ,errorsPosErrCorrSummedGamma );
	meanErrorsSummedGamma = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanSummedGamma ,ptBinsGammaErr ,errorsMeanErrSummedGamma );
	meanErrorsCorrSummedGamma = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanCorrSummedGamma ,ptBinsGammaErr ,errorsMeanErrCorrSummedGamma );
	meanErrorsCorrSummedIncMatGamma = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanCorrMatSummedGamma ,ptBinsGammaErr ,errorsMeanErrCorrMatSummedGamma );

	meanErrorsCorrSummedWithoutMatGamma = new TGraphErrors(ConstnBinsGamma,ptBinsGamma , errorsMeanCorrMatSummedGammaWithoutMat,ptBinsGammaErr , errorsMeanErrCorrMatSummedGammaWithoutMat);
	meanErrorsCorrSummedMaterialGamma = new TGraphErrors(ConstnBinsGamma,ptBinsGamma , errorsMaterialBudget,ptBinsGammaErr , errorsErrMaterialBudget);

	negativeErrorsSummedGammaA = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegSummedGammaA ,ptBinsGammaErr ,errorsNegErrSummedGammaA );
	negativeErrorsCorrSummedGammaA = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegCorrSummedGammaA ,ptBinsGammaErr ,errorsNegErrCorrSummedGammaA );
	positiveErrorsSummedGammaA = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosSummedGammaA ,ptBinsGammaErr ,errorsPosErrSummedGammaA );
	positiveErrorsCorrSummedGammaA = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosCorrSummedGammaA ,ptBinsGammaErr ,errorsPosErrCorrSummedGammaA );
	meanErrorsSummedGammaA = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanSummedGammaA ,ptBinsGammaErr ,errorsMeanErrSummedGammaA );
	meanErrorsCorrSummedGammaA = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanCorrSummedGammaA ,ptBinsGammaErr ,errorsMeanErrCorrSummedGammaA );
	meanErrorsCorrSummedIncMatGammaA = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanCorrMatSummedGammaA ,ptBinsGammaErr ,errorsMeanErrCorrMatSummedGammaA );

	negativeErrorsSummedGammaB = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegSummedGammaB ,ptBinsGammaErr ,errorsNegErrSummedGammaB );
	negativeErrorsCorrSummedGammaB = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegCorrSummedGammaB ,ptBinsGammaErr ,errorsNegErrCorrSummedGammaB );
	positiveErrorsSummedGammaB = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosSummedGammaB ,ptBinsGammaErr ,errorsPosErrSummedGammaB );
	positiveErrorsCorrSummedGammaB = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosCorrSummedGammaB ,ptBinsGammaErr ,errorsPosErrCorrSummedGammaB );
	meanErrorsSummedGammaB = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanSummedGammaB ,ptBinsGammaErr ,errorsMeanErrSummedGammaB );
	meanErrorsCorrSummedGammaB = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanCorrSummedGammaB ,ptBinsGammaErr ,errorsMeanErrCorrSummedGammaB );
	meanErrorsCorrSummedIncMatGammaB = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanCorrMatSummedGammaB ,ptBinsGammaErr ,errorsMeanErrCorrMatSummedGammaB );

	negativeErrorsSummedGammaC = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegSummedGammaC ,ptBinsGammaErr ,errorsNegErrSummedGammaC );
	negativeErrorsCorrSummedGammaC = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsNegCorrSummedGammaC ,ptBinsGammaErr ,errorsNegErrCorrSummedGammaC );
	positiveErrorsSummedGammaC = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosSummedGammaC ,ptBinsGammaErr ,errorsPosErrSummedGammaC );
	positiveErrorsCorrSummedGammaC = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsPosCorrSummedGammaC ,ptBinsGammaErr ,errorsPosErrCorrSummedGammaC );
	meanErrorsSummedGammaC = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanSummedGammaC ,ptBinsGammaErr ,errorsMeanErrSummedGammaC );
	meanErrorsCorrSummedGammaC = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanCorrSummedGammaC ,ptBinsGammaErr ,errorsMeanErrCorrSummedGammaC );
	meanErrorsCorrSummedIncMatGammaC = new TGraphErrors(ConstnBinsGamma,ptBinsGamma ,errorsMeanCorrMatSummedGammaC ,ptBinsGammaErr ,errorsMeanErrCorrMatSummedGammaC );

	cout << "here" << endl;
	Int_t i=0;
	while(ptBinsDoubleRatio[i]<0.3){
		for (Int_t j = 0; j < ErrorDoubleRatioToRun; j++)
			meanErrorsCorrDoubleRatioPi0Fit[j]->RemovePoint(0);
		meanErrorsCorrSummedDoubleRatioPi0Fit->RemovePoint(0);
		meanErrorsCorrSummedIncMatDoubleRatioPi0Fit->RemovePoint(0);
		meanErrorsCorrSummedDoubleRatioPi0FitA->RemovePoint(0);
		meanErrorsCorrSummedDoubleRatioPi0FitB->RemovePoint(0);
		meanErrorsCorrSummedIncMatDoubleRatioPi0FitC->RemovePoint(0);
		i++;
	}		
	
	TCanvas* canvasNewSysErrMean = new TCanvas("canvasNewSysErrMean","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasNewSysErrMean, 0.08, 0.01, 0.015, 0.09);
	
	TH2D *histo2DNewSysErrMean ;
// 	if (meson.CompareTo("Pi0")==0 ){
		histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "histo2DNewSysErrMean", 20,0.,ptBinsDoubleRatio[ConstnBinsDoubleRatio-1]+1,1000.,-0.5,30.);
// 	} else { 
// 		histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "histo2DNewSysErrMean", 20,0.,ptBins[ConstnBinsDoubleRatio-1]+1,1000.,-0.5,50.);
// 	}
	histo2DNewSysErrMean->SetYTitle("mean smoothed systematic Err %");
	histo2DNewSysErrMean->SetXTitle("p_{T} (GeV/c)");
	histo2DNewSysErrMean->GetYaxis()->SetNdivisions(510,kTRUE);
	histo2DNewSysErrMean->GetYaxis()->SetLabelSize(0.03);
	histo2DNewSysErrMean->GetYaxis()->SetTitleSize(0.04);
	histo2DNewSysErrMean->GetYaxis()->SetDecimals();
	histo2DNewSysErrMean->GetYaxis()->SetTitleOffset(0.9);
	histo2DNewSysErrMean->GetXaxis()->SetTitleOffset(1.);
	histo2DNewSysErrMean->GetXaxis()->SetLabelSize(0.03);
	histo2DNewSysErrMean->GetXaxis()->SetTitleSize(0.04);
	histo2DNewSysErrMean->SetTitle("");
	histo2DNewSysErrMean->Draw();
	
	TLegend* legendDoubleRatio;
	legendDoubleRatio= new TLegend(0.15,0.6,0.57,0.95);
	legendDoubleRatio->SetTextSize(0.035);
	legendDoubleRatio->SetFillColor(0);
	legendDoubleRatio->SetBorderSize(0);
	if (ErrorDoubleRatioToRun> 8) legendDoubleRatio->SetNColumns(2);
	for(Int_t i = 0; i< ErrorDoubleRatioToRun ; i++){
		DrawGammaSetMarkerTGraphErr(meanErrorsCorrDoubleRatioPi0Fit[i], 20+i, 1.,color[i],color[i]);
		meanErrorsCorrDoubleRatioPi0Fit[i]->Draw("p,csame");
		legendDoubleRatio->AddEntry(meanErrorsCorrDoubleRatioPi0Fit[i],nameCutVariationsGammaLegend[i].Data(),"p");	
	}
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedDoubleRatioPi0Fit, 20, 1.,1,1);
	meanErrorsCorrSummedDoubleRatioPi0Fit->Draw("p,csame");
	legendDoubleRatio->AddEntry(meanErrorsCorrSummedDoubleRatioPi0Fit,"quad. sum.","p");
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatDoubleRatioPi0Fit, 20, 1.,kRed,kRed);
	meanErrorsCorrSummedIncMatDoubleRatioPi0Fit->Draw("p,csame");
	legendDoubleRatio->AddEntry(meanErrorsCorrSummedIncMatDoubleRatioPi0Fit,"quad. sum., inc. mat.","p");
	legendDoubleRatio->Draw();
	
	TLatex *labelMeasurement;
	labelMeasurement= new TLatex(0.65,0.89,Form("(#gamma_{inc}/#pi^{0})/(#gamma_{decay}/#pi^{0})"));
	SetStyleTLatex( labelMeasurement, 0.038,4);
	labelMeasurement->Draw();

	TLatex *labelCentrality = new TLatex(0.65,0.93,Form("%s p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV"	,additionalName.Data()	));
	SetStyleTLatex( labelCentrality, 0.038,4);
	labelCentrality->Draw();

	canvasNewSysErrMean->Update();
	canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculated/DoubleRatio_SysMeanNewWithMean_%s%s_%s.%s", energy.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));

	histo2DNewSysErrMean->Draw();

	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedDoubleRatioPi0FitA, 21, 1.,kRed-8,kRed-8);
	meanErrorsCorrSummedDoubleRatioPi0FitA->Draw("p,csame");
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedDoubleRatioPi0FitB, 22, 1.,kGreen-8,kGreen-8);
	meanErrorsCorrSummedDoubleRatioPi0FitB->Draw("p,csame");
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatDoubleRatioPi0FitC, 29, 1.,kBlue-8,kBlue-8);
	meanErrorsCorrSummedIncMatDoubleRatioPi0FitC->Draw("p,csame");
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatDoubleRatioPi0Fit, 20, 1.,kRed,kRed);
	meanErrorsCorrSummedIncMatDoubleRatioPi0Fit->Draw("p,csame");
	
	TLegend* legendABCDoubleRatio;
	legendABCDoubleRatio= new TLegend(0.15,0.6,0.57,0.95);
	legendABCDoubleRatio->SetTextSize(0.035);
	legendABCDoubleRatio->SetFillColor(0);
	legendABCDoubleRatio->SetBorderSize(0);
	legendABCDoubleRatio->AddEntry(meanErrorsCorrSummedDoubleRatioPi0FitA,"Type A","p");
	legendABCDoubleRatio->AddEntry(meanErrorsCorrSummedDoubleRatioPi0FitB,"Type B","p");
	legendABCDoubleRatio->AddEntry(meanErrorsCorrSummedIncMatDoubleRatioPi0FitC,"Type C","p");
	legendABCDoubleRatio->AddEntry(meanErrorsCorrSummedIncMatDoubleRatioPi0Fit,"quad. sum.","p");
	legendABCDoubleRatio->Draw();
	
	labelMeasurement->Draw();
	labelCentrality->Draw();

	canvasNewSysErrMean->Update();
	canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculated/DoubleRatio_SysMeanNewWithMeanABC_%s%s_%s.%s", energy.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));

	i=0;
	while(ptBinsIncRatio[i]<0.3){
		for (Int_t j = 0; j < ErrorIncRatioToRun; j++)
			meanErrorsCorrIncRatioPi0Fit[j]->RemovePoint(0);
		meanErrorsCorrSummedIncRatioPi0Fit->RemovePoint(0);
		meanErrorsCorrSummedIncMatIncRatioPi0Fit->RemovePoint(0);
		meanErrorsCorrSummedIncRatioPi0FitA->RemovePoint(0);
		meanErrorsCorrSummedIncRatioPi0FitB->RemovePoint(0);
		meanErrorsCorrSummedIncMatIncRatioPi0FitC->RemovePoint(0);
		i++;
	}		
	
	histo2DNewSysErrMean->Draw();
	
	TLegend* legendIncRatio;
	legendIncRatio= new TLegend(0.15,0.6,0.57,0.95);
	legendIncRatio->SetTextSize(0.035);
	legendIncRatio->SetFillColor(0);
	legendIncRatio->SetBorderSize(0);
	if (ErrorIncRatioToRun> 8) legendIncRatio->SetNColumns(2);
	for(Int_t i = 0; i< ErrorIncRatioToRun ; i++){
		DrawGammaSetMarkerTGraphErr(meanErrorsCorrIncRatioPi0Fit[i], 20+i, 1.,color[i],color[i]);
		meanErrorsCorrIncRatioPi0Fit[i]->Draw("p,csame");
		legendIncRatio->AddEntry(meanErrorsCorrIncRatioPi0Fit[i],nameCutVariationsGammaLegend[i].Data(),"p");	
	}
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncRatioPi0Fit, 20, 1.,1,1);
	meanErrorsCorrSummedIncRatioPi0Fit->Draw("p,csame");
	legendIncRatio->AddEntry(meanErrorsCorrSummedIncRatioPi0Fit,"quad. sum.","p");
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatIncRatioPi0Fit, 20, 1.,kRed,kRed);
	meanErrorsCorrSummedIncMatIncRatioPi0Fit->Draw("p,csame");
	legendIncRatio->AddEntry(meanErrorsCorrSummedIncMatIncRatioPi0Fit,"quad. sum., inc. mat.","p");
	legendIncRatio->Draw();
	
	TLatex *labelMeasurement2;
	labelMeasurement2= new TLatex(0.65,0.89,Form("#gamma_{inc}/#pi^{0}"));
	SetStyleTLatex( labelMeasurement2, 0.038,4);
	labelMeasurement2->Draw();

	labelCentrality->Draw();

	canvasNewSysErrMean->Update();
	canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculated/IncRatio_SysMeanNewWithMean_%s%s_%s.%s", energy.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));

	histo2DNewSysErrMean->Draw();

	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncRatioPi0FitA, 21, 1.,kRed-8,kRed-8);
	meanErrorsCorrSummedIncRatioPi0FitA->Draw("p,csame");
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncRatioPi0FitB, 22, 1.,kGreen-8,kGreen-8);
	meanErrorsCorrSummedIncRatioPi0FitB->Draw("p,csame");
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatIncRatioPi0FitC, 29, 1.,kBlue-8,kBlue-8);
	meanErrorsCorrSummedIncMatIncRatioPi0FitC->Draw("p,csame");
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatIncRatioPi0Fit, 20, 1.,kRed,kRed);
	meanErrorsCorrSummedIncMatIncRatioPi0Fit->Draw("p,csame");
	
	TLegend* legendABCIncRatio;
	legendABCIncRatio= new TLegend(0.15,0.6,0.57,0.95);
	legendABCIncRatio->SetTextSize(0.035);
	legendABCIncRatio->SetFillColor(0);
	legendABCIncRatio->SetBorderSize(0);
	legendABCIncRatio->AddEntry(meanErrorsCorrSummedIncRatioPi0FitA,"Type A","p");
	legendABCIncRatio->AddEntry(meanErrorsCorrSummedIncRatioPi0FitB,"Type B","p");
	legendABCIncRatio->AddEntry(meanErrorsCorrSummedIncMatIncRatioPi0FitC,"Type C","p");
	legendABCIncRatio->AddEntry(meanErrorsCorrSummedIncMatIncRatioPi0Fit,"quad. sum.","p");
	legendABCIncRatio->Draw();
	
	labelMeasurement2->Draw();
	labelCentrality->Draw();

	canvasNewSysErrMean->Update();
	canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculated/IncRatio_SysMeanNewWithMeanABC_%s%s_%s.%s", energy.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));
	
	i=0;
	while(ptBinsGamma[i]<0.3){
		for (Int_t j = 0; j < ErrorGammaSpecToRun; j++)
			meanErrorsCorrGamma[j]->RemovePoint(0);
		meanErrorsCorrSummedGamma->RemovePoint(0);
		meanErrorsCorrSummedIncMatGamma->RemovePoint(0);
		meanErrorsCorrSummedGammaA->RemovePoint(0);
		meanErrorsCorrSummedGammaB->RemovePoint(0);
		meanErrorsCorrSummedIncMatGammaC->RemovePoint(0);
		i++;
	}		
	
	histo2DNewSysErrMean->Draw();
	
	TLegend* legendGamma;
	legendGamma= new TLegend(0.15,0.6,0.57,0.95);
	legendGamma->SetTextSize(0.035);
	legendGamma->SetFillColor(0);
	legendGamma->SetBorderSize(0);
	if (ErrorGammaSpecToRun> 8) legendGamma->SetNColumns(2);
	for(Int_t i = 0; i< ErrorGammaSpecToRun ; i++){
		DrawGammaSetMarkerTGraphErr(meanErrorsCorrGamma[i], 20+i, 1.,color[i],color[i]);
		meanErrorsCorrGamma[i]->Draw("p,csame");
		legendGamma->AddEntry(meanErrorsCorrGamma[i],nameCutVariationsGammaLegend[i].Data(),"p");	
	}
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedGamma, 20, 1.,1,1);
	meanErrorsCorrSummedGamma->Draw("p,csame");
	legendGamma->AddEntry(meanErrorsCorrSummedGamma,"quad. sum.","p");
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatGamma, 20, 1.,kRed,kRed);
	meanErrorsCorrSummedIncMatGamma->Draw("p,csame");
	legendGamma->AddEntry(meanErrorsCorrSummedIncMatGamma,"quad. sum., inc. mat.","p");
	legendGamma->Draw();
	
	TLatex *labelMeasurement3;
	labelMeasurement3= new TLatex(0.65,0.89,Form("#gamma_{inc}"));
	SetStyleTLatex( labelMeasurement3, 0.038,4);
	labelMeasurement3->Draw();

	labelCentrality->Draw();

	canvasNewSysErrMean->Update();
	canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculated/IncGamma_SysMeanNewWithMean_%s%s_%s.%s", energy.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));

	histo2DNewSysErrMean->Draw();

	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedGammaA, 21, 1.,kRed-8,kRed-8);
	meanErrorsCorrSummedGammaA->Draw("p,csame");
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedGammaB, 22, 1.,kGreen-8,kGreen-8);
	meanErrorsCorrSummedGammaB->Draw("p,csame");
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatGammaC, 29, 1.,kBlue-8,kBlue-8);
	meanErrorsCorrSummedIncMatGammaC->Draw("p,csame");
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMatGamma, 20, 1.,kRed,kRed);
	meanErrorsCorrSummedIncMatGamma->Draw("p,csame");
	
	TLegend* legendABCGamma;
	legendABCGamma= new TLegend(0.15,0.6,0.57,0.95);
	legendABCGamma->SetTextSize(0.035);
	legendABCGamma->SetFillColor(0);
	legendABCGamma->SetBorderSize(0);
	legendABCGamma->AddEntry(meanErrorsCorrSummedGammaA,"Type A","p");
	legendABCGamma->AddEntry(meanErrorsCorrSummedGammaB,"Type B","p");
	legendABCGamma->AddEntry(meanErrorsCorrSummedIncMatGammaC,"Type C","p");
	legendABCGamma->AddEntry(meanErrorsCorrSummedIncMatGamma,"quad. sum.","p");
	legendABCGamma->Draw();
	
	labelMeasurement3->Draw();
	labelCentrality->Draw();

	canvasNewSysErrMean->Update();
	canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculated/IncGamma_SysMeanNewWithMeanABC_%s%s_%s.%s", energy.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));

	
	const char *SysErrDatnameDoubleRatio = Form("SystematicErrorsCalculated/SystematicErrorAveraged_DoubleRatio_%s%s_%s.dat",energy.Data(),additionalNameOutput.Data(),dateForOutput.Data());
	fstream SysErrDatAverDoubleRatio;
	cout << SysErrDatnameDoubleRatio << endl;
	SysErrDatAverDoubleRatio.open(SysErrDatnameDoubleRatio, ios::out);
	for (Int_t l=0; l< ConstnBinsDoubleRatio; l++){
		SysErrDatAverDoubleRatio << "-"<< errorsMeanCorrMatSummedDoubleRatioPi0Fit[l] << "\t" <<errorsMeanCorrMatSummedDoubleRatioPi0Fit[l] << "\t"  << "-"<< errorsMeanCorrSummedDoubleRatioPi0Fit[l] << "\t" <<errorsMeanCorrSummedDoubleRatioPi0Fit[l]  << "\t"  << "-" << errorsMeanCorrSummedDoubleRatioPi0FitA[l] << "\t" <<errorsMeanCorrSummedDoubleRatioPi0FitA[l] << "\t"  << "-" << errorsMeanCorrSummedDoubleRatioPi0FitB[l] << "\t" <<errorsMeanCorrSummedDoubleRatioPi0FitB[l] << "\t"  << "-" << errorsMeanCorrMatSummedDoubleRatioPi0FitC[l] << "\t" <<errorsMeanCorrMatSummedDoubleRatioPi0FitC[l]  <<endl;
	}

	SysErrDatAverDoubleRatio.close();

	const char *SysErrDatnameIncRatio = Form("SystematicErrorsCalculated/SystematicErrorAveraged_IncRatio_%s%s_%s.dat",energy.Data(),additionalNameOutput.Data(),dateForOutput.Data());
	fstream SysErrDatAverIncRatio;
	cout << SysErrDatnameIncRatio << endl;
	SysErrDatAverIncRatio.open(SysErrDatnameIncRatio, ios::out);
	for (Int_t l=0; l< ConstnBinsIncRatio; l++){
		SysErrDatAverIncRatio << "-"<< errorsMeanCorrMatSummedIncRatioPi0Fit[l] << "\t" <<errorsMeanCorrMatSummedIncRatioPi0Fit[l] << "\t"  << "-"<< errorsMeanCorrSummedIncRatioPi0Fit[l] << "\t" <<errorsMeanCorrSummedIncRatioPi0Fit[l]  << "\t"  << "-" << errorsMeanCorrSummedIncRatioPi0FitA[l] << "\t" <<errorsMeanCorrSummedIncRatioPi0FitA[l] << "\t"  << "-" << errorsMeanCorrSummedIncRatioPi0FitB[l] << "\t" <<errorsMeanCorrSummedIncRatioPi0FitB[l] << "\t"  << "-" << errorsMeanCorrMatSummedIncRatioPi0FitC[l] << "\t" <<errorsMeanCorrMatSummedIncRatioPi0FitC[l]  <<endl;
	}

	SysErrDatAverIncRatio.close();

	const char *SysErrDatnameGamma = Form("SystematicErrorsCalculated/SystematicErrorAveraged_GammaInc_%s%s_%s.dat",energy.Data(),additionalNameOutput.Data(),dateForOutput.Data());
	fstream SysErrDatAverGamma;
	cout << SysErrDatnameGamma << endl;
	SysErrDatAverGamma.open(SysErrDatnameGamma, ios::out);
	for (Int_t l=0; l< ConstnBinsIncRatio; l++){
		SysErrDatAverGamma << "-"<< errorsMeanCorrMatSummedGamma[l] << "\t" <<errorsMeanCorrMatSummedGamma[l] << "\t"  << "-"<< errorsMeanCorrSummedGamma[l] << "\t" <<errorsMeanCorrSummedGamma[l]  << "\t"  << "-" << errorsMeanCorrSummedGammaA[l] << "\t" <<errorsMeanCorrSummedGammaA[l] << "\t"  << "-" << errorsMeanCorrSummedGammaB[l] << "\t" <<errorsMeanCorrSummedGammaB[l] << "\t"  << "-" << errorsMeanCorrMatSummedGammaC[l] << "\t" <<errorsMeanCorrMatSummedGammaC[l]  <<endl;
	}

	SysErrDatAverGamma.close();

	
}