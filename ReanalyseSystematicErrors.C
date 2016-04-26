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

void ReanalyseSystematicErrors(const char* dataFileNegative ="", const char* dataFilePositive ="", TString energy="", TString meson = "", Int_t numberOfPtBins =1 ,Int_t numberCutStudies=1, Double_t decisionBoundary=1.){
	
	StyleSettingsThesis();	
	SetPlotStyle();
	
	TDatime now;
	int iDate = now.GetDate();
	int iYear=iDate/10000;
	int iMonth=(iDate%10000)/100;
	int iDay=iDate%100;
	char* cMonth[12]={"Jan","Feb","Mar","Apr","May","Jun",
		    "Jul","Aug","Sep","Oct","Nov","Dec"};
	char cStamp1[25],cStamp2[25];
	sprintf(cStamp1,"%i_%s_%i",iDay, cMonth[iMonth-1], iYear);
	sprintf(cStamp2,"%i/%.2d/%i",iDay, iMonth, iYear);
	
	Int_t color[20] = {860,894,620,880,591,403,802,923,634,432,422,435,420,407,416,830,404,608,920,1};
	
	Double_t readoutPos[1000];
	Double_t readoutNeg[1000];
	Int_t 	numberOfEntriesPos = 0;
	Int_t 	numberOfEntriesNeg = 0;
	
	ifstream 		fileSysErrUp;
	fileSysErrUp.open(dataFilePositive,ios_base::in);
	cout << dataFilePositive << endl;
	while(!fileSysErrUp.eof()){
		fileSysErrUp >> readoutPos[numberOfEntriesPos];
		cout << readoutPos[numberOfEntriesPos] << endl;
		numberOfEntriesPos++;
	}
	fileSysErrUp.close();

	ifstream 		fileSysErrDown;
	fileSysErrDown.open(dataFileNegative,ios_base::in);
	cout << dataFileNegative << endl;
	while(!fileSysErrDown.eof()){
		fileSysErrDown >> readoutNeg[numberOfEntriesNeg];
		cout << readoutNeg[numberOfEntriesNeg] << endl;
		numberOfEntriesNeg++;
	}
	fileSysErrDown.close();
	
	const Int_t nPtBins = numberOfPtBins;
	const Int_t nCutStudies = numberCutStudies;
	Double_t ptBins[nPtBins];
	Double_t ptBinsErr[nPtBins];
	
	TString nameCutVariation7TeV[10] = {"Cluster variation", "Single E/P p_{t} variation","dE/dx variation", "#chi^{2} #gamma variation", "#alpha variation","diff Periods","Yield extraction","BG variation","",""};  
	TString nameCutForOutputVariation7TeV[10] = {"Cluster", "SingleEPPt","dEdx", "Chi2Gamma", "Alpha","Periods","YieldExtraction","BG","",""};  

	TString nameCutVariation900GeV[10] = {"Cluster variation", "Single E/P p_{t} variation","dE/dx variation", "#chi^{2} #gamma variation", "#alpha variation","Yield extraction","BG variation","","",""};  

	TString nameCutVariation[10];
	
	Int_t offsetAtEnd;
	if (energy.CompareTo("7TeV") == 0) {
		offsetAtEnd = 3;
		for (Int_t i = 0; i < 10; i++){
			nameCutVariation[i] = nameCutVariation7TeV[i];
		}
	} else {
		offsetAtEnd = 2;
		for (Int_t i = 0; i < 10; i++){
			nameCutVariation[i] = nameCutVariation900GeV[i];
		}
	}
	
	cout<< "pt Bins" << endl;
	for (Int_t i = 0; i < nPtBins; i++){
		ptBins[i] = readoutPos[i*(numberCutStudies*3+offsetAtEnd+1)];
		ptBinsErr[i] = 0;
	}
	
	Double_t errorsNegCutSt1[nPtBins];
	Double_t errorsNegCutSt2[nPtBins];
	Double_t errorsNegCutSt3[nPtBins];
	Double_t errorsNegCutSt4[nPtBins];
	Double_t errorsNegCutSt5[nPtBins];
	Double_t errorsNegCutSt6[nPtBins];
	Double_t errorsNegCutSt7[nPtBins];
	Double_t errorsNegCutSt8[nPtBins];
	Double_t errorsNegCutSt9[nPtBins];
	Double_t errorsNegCutSt10[nPtBins];
	Double_t errorsNegCorrCutSt1[nPtBins];
	Double_t errorsNegCorrCutSt2[nPtBins];
	Double_t errorsNegCorrCutSt3[nPtBins];
	Double_t errorsNegCorrCutSt4[nPtBins];
	Double_t errorsNegCorrCutSt5[nPtBins];
	Double_t errorsNegCorrCutSt6[nPtBins];
	Double_t errorsNegCorrCutSt7[nPtBins];
	Double_t errorsNegCorrCutSt8[nPtBins];
	Double_t errorsNegCorrCutSt9[nPtBins];
	Double_t errorsNegCorrCutSt10[nPtBins];
	Double_t errorsNegSummed[nPtBins];
	Double_t errorsNegCorrSummed[nPtBins];
	Double_t errorsNegCorrMatSummed[nPtBins];
	
	Double_t errorsNegErrCutSt1[nPtBins];
	Double_t errorsNegErrCutSt2[nPtBins];
	Double_t errorsNegErrCutSt3[nPtBins];
	Double_t errorsNegErrCutSt4[nPtBins];
	Double_t errorsNegErrCutSt5[nPtBins];
	Double_t errorsNegErrCutSt6[nPtBins];
	Double_t errorsNegErrCutSt7[nPtBins];
	Double_t errorsNegErrCutSt8[nPtBins];
	Double_t errorsNegErrCutSt9[nPtBins];
	Double_t errorsNegErrCutSt10[nPtBins];
	Double_t errorsNegErrCorrCutSt1[nPtBins];
	Double_t errorsNegErrCorrCutSt2[nPtBins];
	Double_t errorsNegErrCorrCutSt3[nPtBins];
	Double_t errorsNegErrCorrCutSt4[nPtBins];
	Double_t errorsNegErrCorrCutSt5[nPtBins];
	Double_t errorsNegErrCorrCutSt6[nPtBins];
	Double_t errorsNegErrCorrCutSt7[nPtBins];
	Double_t errorsNegErrCorrCutSt8[nPtBins];
	Double_t errorsNegErrCorrCutSt9[nPtBins];
	Double_t errorsNegErrCorrCutSt10[nPtBins];
	Double_t errorsNegErrSummed[nPtBins];
	Double_t errorsNegErrCorrSummed[nPtBins];
	
	Double_t errorsPosCutSt1[nPtBins];
	Double_t errorsPosCutSt2[nPtBins];
	Double_t errorsPosCutSt3[nPtBins];
	Double_t errorsPosCutSt4[nPtBins];
	Double_t errorsPosCutSt5[nPtBins];
	Double_t errorsPosCutSt6[nPtBins];
	Double_t errorsPosCutSt7[nPtBins];
	Double_t errorsPosCutSt8[nPtBins];
	Double_t errorsPosCutSt9[nPtBins];
	Double_t errorsPosCutSt10[nPtBins];
	Double_t errorsPosCorrCutSt1[nPtBins];
	Double_t errorsPosCorrCutSt2[nPtBins];
	Double_t errorsPosCorrCutSt3[nPtBins];
	Double_t errorsPosCorrCutSt4[nPtBins];
	Double_t errorsPosCorrCutSt5[nPtBins];
	Double_t errorsPosCorrCutSt6[nPtBins];
	Double_t errorsPosCorrCutSt7[nPtBins];
	Double_t errorsPosCorrCutSt8[nPtBins];
	Double_t errorsPosCorrCutSt9[nPtBins];
	Double_t errorsPosCorrCutSt10[nPtBins];
	Double_t errorsPosSummed[nPtBins];
	Double_t errorsPosCorrSummed[nPtBins];
	Double_t errorsPosCorrMatSummed[nPtBins];
	
	Double_t errorsPosErrCutSt1[nPtBins];
	Double_t errorsPosErrCutSt2[nPtBins];
	Double_t errorsPosErrCutSt3[nPtBins];
	Double_t errorsPosErrCutSt4[nPtBins];
	Double_t errorsPosErrCutSt5[nPtBins];
	Double_t errorsPosErrCutSt6[nPtBins];
	Double_t errorsPosErrCutSt7[nPtBins];
	Double_t errorsPosErrCutSt8[nPtBins];
	Double_t errorsPosErrCutSt9[nPtBins];
	Double_t errorsPosErrCutSt10[nPtBins];
	Double_t errorsPosErrSummed[nPtBins];
	Double_t errorsPosErrCorrCutSt1[nPtBins];
	Double_t errorsPosErrCorrCutSt2[nPtBins];
	Double_t errorsPosErrCorrCutSt3[nPtBins];
	Double_t errorsPosErrCorrCutSt4[nPtBins];
	Double_t errorsPosErrCorrCutSt5[nPtBins];
	Double_t errorsPosErrCorrCutSt6[nPtBins];
	Double_t errorsPosErrCorrCutSt7[nPtBins];
	Double_t errorsPosErrCorrCutSt8[nPtBins];
	Double_t errorsPosErrCorrCutSt9[nPtBins];
	Double_t errorsPosErrCorrCutSt10[nPtBins];
	Double_t errorsPosErrCorrSummed[nPtBins];
	
	Double_t errorsMeanCutSt1[nPtBins];
	Double_t errorsMeanCutSt2[nPtBins];
	Double_t errorsMeanCutSt3[nPtBins];
	Double_t errorsMeanCutSt4[nPtBins];
	Double_t errorsMeanCutSt5[nPtBins];
	Double_t errorsMeanCutSt6[nPtBins];
	Double_t errorsMeanCutSt7[nPtBins];
	Double_t errorsMeanCutSt8[nPtBins];
	Double_t errorsMeanCutSt9[nPtBins];
	Double_t errorsMeanCutSt10[nPtBins];
	Double_t errorsMeanCorrCutSt1[nPtBins];
	Double_t errorsMeanCorrCutSt2[nPtBins];
	Double_t errorsMeanCorrCutSt3[nPtBins];
	Double_t errorsMeanCorrCutSt4[nPtBins];
	Double_t errorsMeanCorrCutSt5[nPtBins];
	Double_t errorsMeanCorrCutSt6[nPtBins];
	Double_t errorsMeanCorrCutSt7[nPtBins];
	Double_t errorsMeanCorrCutSt8[nPtBins];
	Double_t errorsMeanCorrCutSt9[nPtBins];
	Double_t errorsMeanCorrCutSt10[nPtBins];
	Double_t errorsMeanSummed[nPtBins];
	Double_t errorsMeanCorrSummed[nPtBins];
	Double_t errorsMeanCorrMatSummed[nPtBins];
	
	Double_t errorsMeanErrCutSt1[nPtBins];
	Double_t errorsMeanErrCutSt2[nPtBins];
	Double_t errorsMeanErrCutSt3[nPtBins];
	Double_t errorsMeanErrCutSt4[nPtBins];
	Double_t errorsMeanErrCutSt5[nPtBins];
	Double_t errorsMeanErrCutSt6[nPtBins];
	Double_t errorsMeanErrCutSt7[nPtBins];
	Double_t errorsMeanErrCutSt8[nPtBins];
	Double_t errorsMeanErrCutSt9[nPtBins];
	Double_t errorsMeanErrCutSt10[nPtBins];
	Double_t errorsMeanErrCorrCutSt1[nPtBins];
	Double_t errorsMeanErrCorrCutSt2[nPtBins];
	Double_t errorsMeanErrCorrCutSt3[nPtBins];
	Double_t errorsMeanErrCorrCutSt4[nPtBins];
	Double_t errorsMeanErrCorrCutSt5[nPtBins];
	Double_t errorsMeanErrCorrCutSt6[nPtBins];
	Double_t errorsMeanErrCorrCutSt7[nPtBins];
	Double_t errorsMeanErrCorrCutSt8[nPtBins];
	Double_t errorsMeanErrCorrCutSt9[nPtBins];
	Double_t errorsMeanErrCorrCutSt10[nPtBins];
	Double_t errorsMeanErrSummed[nPtBins];
	Double_t errorsMeanErrCorrSummed[nPtBins];
	Double_t errorsMeanErrCorrMatSummed[nPtBins];
	
	TGraphErrors* negativeErrors[nCutStudies];
	TGraphErrors* negativeErrorsSummed;
	TGraphErrors* positiveErrors[nCutStudies];
	TGraphErrors* positiveErrorsSummed;
	TGraphErrors* negativeErrorsCorr[nCutStudies];
	TGraphErrors* negativeErrorsCorrSummed;
	TGraphErrors* positiveErrorsCorr[nCutStudies];
	TGraphErrors* positiveErrorsCorrSummed;
	TGraphErrors* meanErrors[nCutStudies];
	TGraphErrors* meanErrorsSummed;
	TGraphErrors* meanErrorsCorr[nCutStudies];
	TGraphErrors* meanErrorsCorrSummed;
	TGraphErrors* meanErrorsCorrSummedIncMat;
	
	if (numberCutStudies >= 1){
		cout << "Cutstudies 1" << endl;
		ReadOutFromSysErrVector(readoutPos, errorsPosCutSt1,errorsPosErrCutSt1,1, numberCutStudies ,offsetAtEnd ,nPtBins,decisionBoundary);
		ReadOutFromSysErrVector(readoutNeg, errorsNegCutSt1,errorsNegErrCutSt1,1, numberCutStudies ,offsetAtEnd ,nPtBins,decisionBoundary);
		
		CalculateMeanSysErr(errorsMeanCutSt1, errorsMeanErrCutSt1, errorsPosCutSt1, errorsNegCutSt1, nPtBins);
		
		CorrectSystematicErrorsWithMean(errorsPosCutSt1,errorsPosErrCutSt1, errorsPosCorrCutSt1, errorsPosErrCorrCutSt1, nPtBins);
		CorrectSystematicErrorsWithMean(errorsNegCutSt1,errorsNegErrCutSt1, errorsNegCorrCutSt1, errorsNegErrCorrCutSt1, nPtBins);
		CorrectSystematicErrorsWithMean(errorsMeanCutSt1, errorsMeanErrCutSt1, errorsMeanCorrCutSt1, errorsMeanErrCorrCutSt1, nPtBins);
		
		negativeErrors[0] = new TGraphErrors(nPtBins,ptBins ,errorsNegCutSt1 ,ptBinsErr ,errorsNegErrCutSt1 );
		meanErrors[0] = new TGraphErrors(nPtBins,ptBins ,errorsMeanCutSt1 ,ptBinsErr ,errorsMeanErrCutSt1 );
		positiveErrors[0] = new TGraphErrors(nPtBins,ptBins ,errorsPosCutSt1 ,ptBinsErr ,errorsPosErrCutSt1 );
		negativeErrorsCorr[0] = new TGraphErrors(nPtBins,ptBins ,errorsNegCorrCutSt1 ,ptBinsErr ,errorsNegErrCorrCutSt1 );
		meanErrorsCorr[0] = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrCutSt1 ,ptBinsErr ,errorsMeanErrCorrCutSt1 );
		positiveErrorsCorr[0] = new TGraphErrors(nPtBins,ptBins ,errorsPosCorrCutSt1 ,ptBinsErr ,errorsPosErrCorrCutSt1 );
		for (Int_t l = 0; l < nPtBins; l++){
			errorsPosSummed[l] = pow(errorsPosCutSt1[l],2);
			errorsNegSummed[l] = pow(errorsNegCutSt1[l],2);
			errorsMeanSummed[l] = pow(errorsMeanCutSt1[l],2);
			errorsPosCorrSummed[l] = pow(errorsPosCorrCutSt1[l],2);
			errorsNegCorrSummed[l] = pow(errorsNegCorrCutSt1[l],2);
			errorsMeanCorrSummed[l] = pow(errorsMeanCorrCutSt1[l],2);
			cout << l <<errorsPosCorrSummed[l] << "\t" << errorsNegCorrSummed[l] << "\t" << errorsMeanCorrSummed[l]<<endl;
		} 
	}
	if (numberCutStudies >= 2){
		cout << "Cutstudies 2" << endl;
		ReadOutFromSysErrVector(readoutPos, errorsPosCutSt2,errorsPosErrCutSt2,2, numberCutStudies ,offsetAtEnd ,nPtBins,decisionBoundary);
		ReadOutFromSysErrVector(readoutNeg, errorsNegCutSt2,errorsNegErrCutSt2,2, numberCutStudies ,offsetAtEnd ,nPtBins,decisionBoundary);
		CalculateMeanSysErr(errorsMeanCutSt2, errorsMeanErrCutSt2, errorsPosCutSt2, errorsNegCutSt2, nPtBins);
		
		CorrectSystematicErrorsWithMean(errorsPosCutSt2,errorsPosErrCutSt2, errorsPosCorrCutSt2, errorsPosErrCorrCutSt2, nPtBins);
		CorrectSystematicErrorsWithMean(errorsNegCutSt2,errorsNegErrCutSt2, errorsNegCorrCutSt2, errorsNegErrCorrCutSt2, nPtBins);
		CorrectSystematicErrorsWithMean(errorsMeanCutSt2, errorsMeanErrCutSt2, errorsMeanCorrCutSt2, errorsMeanErrCorrCutSt2, nPtBins);
		
		negativeErrors[1] = new TGraphErrors(nPtBins,ptBins ,errorsNegCutSt2 ,ptBinsErr ,errorsNegErrCutSt2 );
		meanErrors[1] = new TGraphErrors(nPtBins,ptBins ,errorsMeanCutSt2 ,ptBinsErr ,errorsMeanErrCutSt2 );
		positiveErrors[1] = new TGraphErrors(nPtBins,ptBins ,errorsPosCutSt2 ,ptBinsErr ,errorsPosErrCutSt2 );
		negativeErrorsCorr[1] = new TGraphErrors(nPtBins,ptBins ,errorsNegCorrCutSt2 ,ptBinsErr ,errorsNegErrCorrCutSt2 );
		meanErrorsCorr[1] = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrCutSt2 ,ptBinsErr ,errorsMeanErrCorrCutSt2 );
		positiveErrorsCorr[1] = new TGraphErrors(nPtBins,ptBins ,errorsPosCorrCutSt2 ,ptBinsErr ,errorsPosErrCorrCutSt2 );

		for (Int_t l = 0; l < nPtBins; l++){
			errorsPosSummed[l] = errorsPosSummed[l] + pow(errorsPosCutSt2[l],2);
			errorsNegSummed[l] = errorsNegSummed[l] + pow(errorsNegCutSt2[l],2);
			errorsMeanSummed[l] = errorsMeanSummed[l] + pow(errorsMeanCutSt2[l],2);
			errorsPosCorrSummed[l] = errorsPosCorrSummed[l] + pow(errorsPosCorrCutSt2[l],2);
			errorsNegCorrSummed[l] = errorsNegCorrSummed[l] + pow(errorsNegCorrCutSt2[l],2);
			errorsMeanCorrSummed[l] = errorsMeanCorrSummed[l] + pow(errorsMeanCorrCutSt2[l],2);
			
		} 
	}
	if (numberCutStudies >= 3){	
		cout << "Cutstudies 3" << endl;
		ReadOutFromSysErrVector(readoutPos, errorsPosCutSt3,errorsPosErrCutSt3,3, numberCutStudies ,offsetAtEnd ,nPtBins,decisionBoundary);
		ReadOutFromSysErrVector(readoutNeg, errorsNegCutSt3,errorsNegErrCutSt3,3, numberCutStudies ,offsetAtEnd ,nPtBins,decisionBoundary);
		CalculateMeanSysErr(errorsMeanCutSt3, errorsMeanErrCutSt3, errorsPosCutSt3, errorsNegCutSt3, nPtBins);
		
		CorrectSystematicErrorsWithMean(errorsPosCutSt3,errorsPosErrCutSt3, errorsPosCorrCutSt3, errorsPosErrCorrCutSt3, nPtBins);
		CorrectSystematicErrorsWithMean(errorsNegCutSt3,errorsNegErrCutSt3, errorsNegCorrCutSt3, errorsNegErrCorrCutSt3, nPtBins);
		CorrectSystematicErrorsWithMean(errorsMeanCutSt3, errorsMeanErrCutSt3, errorsMeanCorrCutSt3, errorsMeanErrCorrCutSt3, nPtBins);
			
		negativeErrors[2] = new TGraphErrors(nPtBins,ptBins ,errorsNegCutSt3 ,ptBinsErr ,errorsNegErrCutSt3 );
		meanErrors[2] = new TGraphErrors(nPtBins,ptBins ,errorsMeanCutSt3 ,ptBinsErr ,errorsMeanErrCutSt3 );
		positiveErrors[2] = new TGraphErrors(nPtBins,ptBins ,errorsPosCutSt3 ,ptBinsErr ,errorsPosErrCutSt3 );
		negativeErrorsCorr[2] = new TGraphErrors(nPtBins,ptBins ,errorsNegCorrCutSt3 ,ptBinsErr ,errorsNegErrCorrCutSt3 );
		meanErrorsCorr[2] = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrCutSt3 ,ptBinsErr ,errorsMeanErrCorrCutSt3 );
		positiveErrorsCorr[2] = new TGraphErrors(nPtBins,ptBins ,errorsPosCorrCutSt3 ,ptBinsErr ,errorsPosErrCorrCutSt3 );

		for (Int_t l = 0; l < nPtBins; l++){
			errorsPosSummed[l] = errorsPosSummed[l] + pow(errorsPosCutSt3[l],2);
			errorsNegSummed[l] = errorsNegSummed[l] + pow(errorsNegCutSt3[l],2);
			errorsMeanSummed[l] = errorsMeanSummed[l] + pow(errorsMeanCutSt3[l],2);
			errorsPosCorrSummed[l] = errorsPosCorrSummed[l] + pow(errorsPosCorrCutSt3[l],2);
			errorsNegCorrSummed[l] = errorsNegCorrSummed[l] + pow(errorsNegCorrCutSt3[l],2);
			errorsMeanCorrSummed[l] = errorsMeanCorrSummed[l] + pow(errorsMeanCorrCutSt3[l],2);
		}
	}
	if (numberCutStudies >= 4){
		cout << "Cutstudies 4" << endl;
		ReadOutFromSysErrVector(readoutPos, errorsPosCutSt4,errorsPosErrCutSt4,4, numberCutStudies ,offsetAtEnd ,nPtBins,decisionBoundary);
		ReadOutFromSysErrVector(readoutNeg, errorsNegCutSt4,errorsNegErrCutSt4,4, numberCutStudies ,offsetAtEnd ,nPtBins,decisionBoundary);
		CalculateMeanSysErr(errorsMeanCutSt4, errorsMeanErrCutSt4, errorsPosCutSt4, errorsNegCutSt4, nPtBins);
		
		CorrectSystematicErrorsWithMean(errorsPosCutSt4,errorsPosErrCutSt4, errorsPosCorrCutSt4, errorsPosErrCorrCutSt4, nPtBins);
		CorrectSystematicErrorsWithMean(errorsNegCutSt4,errorsNegErrCutSt4, errorsNegCorrCutSt4, errorsNegErrCorrCutSt4, nPtBins);
		CorrectSystematicErrorsWithMean(errorsMeanCutSt4, errorsMeanErrCutSt4, errorsMeanCorrCutSt4, errorsMeanErrCorrCutSt4, nPtBins);
			
		negativeErrors[3] = new TGraphErrors(nPtBins,ptBins ,errorsNegCutSt4 ,ptBinsErr ,errorsNegErrCutSt4 );
		meanErrors[3] = new TGraphErrors(nPtBins,ptBins ,errorsMeanCutSt4 ,ptBinsErr ,errorsMeanErrCutSt4 );
		positiveErrors[3] = new TGraphErrors(nPtBins,ptBins ,errorsPosCutSt4 ,ptBinsErr ,errorsPosErrCutSt4 );
		negativeErrorsCorr[3] = new TGraphErrors(nPtBins,ptBins ,errorsNegCorrCutSt4 ,ptBinsErr ,errorsNegErrCorrCutSt4 );
		meanErrorsCorr[3] = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrCutSt4 ,ptBinsErr ,errorsMeanErrCorrCutSt4 );
		positiveErrorsCorr[3] = new TGraphErrors(nPtBins,ptBins ,errorsPosCorrCutSt4 ,ptBinsErr ,errorsPosErrCorrCutSt4 );
		
		for (Int_t l = 0; l < nPtBins; l++){
			errorsPosSummed[l] = errorsPosSummed[l] + pow(errorsPosCutSt4[l],2);
			errorsNegSummed[l] = errorsNegSummed[l] + pow(errorsNegCutSt4[l],2);
			errorsMeanSummed[l] = errorsMeanSummed[l] + pow(errorsMeanCutSt4[l],2);
			errorsPosCorrSummed[l] = errorsPosCorrSummed[l] + pow(errorsPosCorrCutSt4[l],2);
			errorsNegCorrSummed[l] = errorsNegCorrSummed[l] + pow(errorsNegCorrCutSt4[l],2);
			errorsMeanCorrSummed[l] = errorsMeanCorrSummed[l] + pow(errorsMeanCorrCutSt4[l],2);
		}
	}
	if (numberCutStudies >= 5){	
		cout << "Cutstudies 5" << endl;
		ReadOutFromSysErrVector(readoutPos, errorsPosCutSt5,errorsPosErrCutSt5,5, numberCutStudies ,offsetAtEnd ,nPtBins,decisionBoundary);
		ReadOutFromSysErrVector(readoutNeg, errorsNegCutSt5,errorsNegErrCutSt5,5, numberCutStudies ,offsetAtEnd ,nPtBins,decisionBoundary);
		CalculateMeanSysErr(errorsMeanCutSt5, errorsMeanErrCutSt5, errorsPosCutSt5, errorsNegCutSt5, nPtBins);
		
		CorrectSystematicErrorsWithMean(errorsPosCutSt5,errorsPosErrCutSt5, errorsPosCorrCutSt5, errorsPosErrCorrCutSt5, nPtBins);
		CorrectSystematicErrorsWithMean(errorsNegCutSt5,errorsNegErrCutSt5, errorsNegCorrCutSt5, errorsNegErrCorrCutSt5, nPtBins);
		CorrectSystematicErrorsWithMean(errorsMeanCutSt5, errorsMeanErrCutSt5, errorsMeanCorrCutSt5, errorsMeanErrCorrCutSt5, nPtBins);
		
		negativeErrors[4] = new TGraphErrors(nPtBins,ptBins ,errorsNegCutSt5 ,ptBinsErr ,errorsNegErrCutSt5 );
		meanErrors[4] = new TGraphErrors(nPtBins,ptBins ,errorsMeanCutSt5 ,ptBinsErr ,errorsMeanErrCutSt5 );
		positiveErrors[4] = new TGraphErrors(nPtBins,ptBins ,errorsPosCutSt5 ,ptBinsErr ,errorsPosErrCutSt5 );
		negativeErrorsCorr[4] = new TGraphErrors(nPtBins,ptBins ,errorsNegCorrCutSt5 ,ptBinsErr ,errorsNegErrCorrCutSt5 );
		meanErrorsCorr[4] = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrCutSt5 ,ptBinsErr ,errorsMeanErrCorrCutSt5 );
		positiveErrorsCorr[4] = new TGraphErrors(nPtBins,ptBins ,errorsPosCorrCutSt5 ,ptBinsErr ,errorsPosErrCorrCutSt5 );
		
		for (Int_t l = 0; l < nPtBins; l++){
			errorsPosSummed[l] = errorsPosSummed[l] + pow(errorsPosCutSt5[l],2);
			errorsNegSummed[l] = errorsNegSummed[l] + pow(errorsNegCutSt5[l],2);
			errorsMeanSummed[l] = errorsMeanSummed[l] + pow(errorsMeanCutSt5[l],2);
			errorsPosCorrSummed[l] = errorsPosCorrSummed[l] + pow(errorsPosCorrCutSt5[l],2);
			errorsNegCorrSummed[l] = errorsNegCorrSummed[l] + pow(errorsNegCorrCutSt5[l],2);
			errorsMeanCorrSummed[l] = errorsMeanCorrSummed[l] + pow(errorsMeanCorrCutSt5[l],2);
		}
	}
	if (numberCutStudies >= 6){	
		cout << "Cutstudies 6" << endl;
		ReadOutFromSysErrVector(readoutPos, errorsPosCutSt6,errorsPosErrCutSt6,6, numberCutStudies ,offsetAtEnd ,nPtBins,decisionBoundary);
		ReadOutFromSysErrVector(readoutNeg, errorsNegCutSt6,errorsNegErrCutSt6,6, numberCutStudies ,offsetAtEnd ,nPtBins,decisionBoundary);
		CalculateMeanSysErr(errorsMeanCutSt6, errorsMeanErrCutSt6, errorsPosCutSt6, errorsNegCutSt6, nPtBins);
			
		CorrectSystematicErrorsWithMean(errorsPosCutSt6,errorsPosErrCutSt6, errorsPosCorrCutSt6, errorsPosErrCorrCutSt6, nPtBins);
		CorrectSystematicErrorsWithMean(errorsNegCutSt6,errorsNegErrCutSt6, errorsNegCorrCutSt6, errorsNegErrCorrCutSt6, nPtBins);
		CorrectSystematicErrorsWithMean(errorsMeanCutSt6, errorsMeanErrCutSt6, errorsMeanCorrCutSt6, errorsMeanErrCorrCutSt6, nPtBins);
		
		negativeErrors[5] = new TGraphErrors(nPtBins,ptBins ,errorsNegCutSt6 ,ptBinsErr ,errorsNegErrCutSt6 );
		meanErrors[5] = new TGraphErrors(nPtBins,ptBins ,errorsMeanCutSt6 ,ptBinsErr ,errorsMeanErrCutSt6 );
		positiveErrors[5] = new TGraphErrors(nPtBins,ptBins ,errorsPosCutSt6 ,ptBinsErr ,errorsPosErrCutSt6 );
		negativeErrorsCorr[5] = new TGraphErrors(nPtBins,ptBins ,errorsNegCorrCutSt6 ,ptBinsErr ,errorsNegErrCorrCutSt6 );
		meanErrorsCorr[5] = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrCutSt6 ,ptBinsErr ,errorsMeanErrCorrCutSt6 );
		positiveErrorsCorr[5] = new TGraphErrors(nPtBins,ptBins ,errorsPosCorrCutSt6 ,ptBinsErr ,errorsPosErrCorrCutSt6 );
		
		if (energy.CompareTo("7TeV")!=0){
			for (Int_t l = 0; l < nPtBins; l++){
				errorsPosSummed[l] = errorsPosSummed[l] + pow(errorsPosCutSt6[l],2);
				errorsNegSummed[l] = errorsNegSummed[l] + pow(errorsNegCutSt6[l],2);
				errorsMeanSummed[l] = errorsMeanSummed[l] + pow(errorsMeanCutSt6[l],2);
				errorsPosCorrSummed[l] = errorsPosCorrSummed[l] + pow(errorsPosCorrCutSt6[l],2);
				errorsNegCorrSummed[l] = errorsNegCorrSummed[l] + pow(errorsNegCorrCutSt6[l],2);
				errorsMeanCorrSummed[l] = errorsMeanCorrSummed[l] + pow(errorsMeanCorrCutSt6[l],2);
			}
		}
	}
	if (numberCutStudies >= 7){	
		cout << "Cutstudies 7" << endl;
		ReadOutFromSysErrVector(readoutPos, errorsPosCutSt7,errorsPosErrCutSt7,7, numberCutStudies ,offsetAtEnd ,nPtBins,decisionBoundary);
		ReadOutFromSysErrVector(readoutNeg, errorsNegCutSt7,errorsNegErrCutSt7,7, numberCutStudies ,offsetAtEnd ,nPtBins,decisionBoundary);
		CalculateMeanSysErr(errorsMeanCutSt7, errorsMeanErrCutSt7, errorsPosCutSt7, errorsNegCutSt7, nPtBins);
		
		CorrectSystematicErrorsWithMean(errorsPosCutSt7,errorsPosErrCutSt7, errorsPosCorrCutSt7, errorsPosErrCorrCutSt7, nPtBins);
		CorrectSystematicErrorsWithMean(errorsNegCutSt7,errorsNegErrCutSt7, errorsNegCorrCutSt7, errorsNegErrCorrCutSt7, nPtBins);
		CorrectSystematicErrorsWithMean(errorsMeanCutSt7, errorsMeanErrCutSt7, errorsMeanCorrCutSt7, errorsMeanErrCorrCutSt7, nPtBins);
		
		negativeErrors[6] = new TGraphErrors(nPtBins,ptBins ,errorsNegCutSt7 ,ptBinsErr ,errorsNegErrCutSt7 );
		meanErrors[6] = new TGraphErrors(nPtBins,ptBins ,errorsMeanCutSt7 ,ptBinsErr ,errorsMeanErrCutSt7 );
		positiveErrors[6] = new TGraphErrors(nPtBins,ptBins ,errorsPosCutSt7 ,ptBinsErr ,errorsPosErrCutSt7 );
		negativeErrorsCorr[6] = new TGraphErrors(nPtBins,ptBins ,errorsNegCorrCutSt7 ,ptBinsErr ,errorsNegErrCorrCutSt7 );
		meanErrorsCorr[6] = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrCutSt7 ,ptBinsErr ,errorsMeanErrCorrCutSt7 );
		positiveErrorsCorr[6] = new TGraphErrors(nPtBins,ptBins ,errorsPosCorrCutSt7 ,ptBinsErr ,errorsPosErrCorrCutSt7 );
		
		for (Int_t l = 0; l < nPtBins; l++){
			errorsPosSummed[l] = errorsPosSummed[l] + pow(errorsPosCutSt7[l],2);
			errorsNegSummed[l] = errorsNegSummed[l] + pow(errorsNegCutSt7[l],2);
			errorsMeanSummed[l] = errorsMeanSummed[l] + pow(errorsMeanCutSt7[l],2);
			errorsPosCorrSummed[l] = errorsPosCorrSummed[l] + pow(errorsPosCorrCutSt7[l],2);
			errorsNegCorrSummed[l] = errorsNegCorrSummed[l] + pow(errorsNegCorrCutSt7[l],2);
			errorsMeanCorrSummed[l] = errorsMeanCorrSummed[l] + pow(errorsMeanCorrCutSt7[l],2);
		}
	}
	if (numberCutStudies >= 8){	
		cout << "Cutstudies 8" << endl;
		ReadOutFromSysErrVector(readoutPos, errorsPosCutSt8,errorsPosErrCutSt8,8, numberCutStudies ,offsetAtEnd ,nPtBins,decisionBoundary);
		ReadOutFromSysErrVector(readoutNeg, errorsNegCutSt8,errorsNegErrCutSt8,8, numberCutStudies ,offsetAtEnd ,nPtBins,decisionBoundary);
		CalculateMeanSysErr(errorsMeanCutSt8, errorsMeanErrCutSt8, errorsPosCutSt8, errorsNegCutSt8, nPtBins);
		
		CorrectSystematicErrorsWithMean(errorsPosCutSt8,errorsPosErrCutSt8, errorsPosCorrCutSt8, errorsPosErrCorrCutSt8, nPtBins);
		CorrectSystematicErrorsWithMean(errorsNegCutSt8,errorsNegErrCutSt8, errorsNegCorrCutSt8, errorsNegErrCorrCutSt8, nPtBins);
		CorrectSystematicErrorsWithMean(errorsMeanCutSt8, errorsMeanErrCutSt8, errorsMeanCorrCutSt8, errorsMeanErrCorrCutSt8, nPtBins);
		
		negativeErrors[7] = new TGraphErrors(nPtBins,ptBins ,errorsNegCutSt8 ,ptBinsErr ,errorsNegErrCutSt8 );
		meanErrors[7] = new TGraphErrors(nPtBins,ptBins ,errorsMeanCutSt8 ,ptBinsErr ,errorsMeanErrCutSt8 );
		positiveErrors[7] = new TGraphErrors(nPtBins,ptBins ,errorsPosCutSt8 ,ptBinsErr ,errorsPosErrCutSt8 );
		negativeErrorsCorr[7] = new TGraphErrors(nPtBins,ptBins ,errorsNegCorrCutSt8 ,ptBinsErr ,errorsNegErrCorrCutSt8 );
		meanErrorsCorr[7] = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrCutSt8 ,ptBinsErr ,errorsMeanErrCorrCutSt8 );
		positiveErrorsCorr[7] = new TGraphErrors(nPtBins,ptBins ,errorsPosCorrCutSt8 ,ptBinsErr ,errorsPosErrCorrCutSt8 );
		
		for (Int_t l = 0; l < nPtBins; l++){
			errorsPosSummed[l] = errorsPosSummed[l] + pow(errorsPosCutSt8[l],2);
			errorsNegSummed[l] = errorsNegSummed[l] + pow(errorsNegCutSt8[l],2);
			errorsMeanSummed[l] = errorsMeanSummed[l] + pow(errorsMeanCutSt8[l],2);
			errorsPosCorrSummed[l] = errorsPosCorrSummed[l] + pow(errorsPosCorrCutSt8[l],2);
			errorsNegCorrSummed[l] = errorsNegCorrSummed[l] + pow(errorsNegCorrCutSt8[l],2);
			errorsMeanCorrSummed[l] = errorsMeanCorrSummed[l] + pow(errorsMeanCorrCutSt8[l],2);
		}
	}
	if (numberCutStudies >= 9){
		cout << "Cutstudies 9" << endl;
		ReadOutFromSysErrVector(readoutPos, errorsPosCutSt9,errorsPosErrCutSt9,9, numberCutStudies ,offsetAtEnd ,nPtBins,decisionBoundary);	
		ReadOutFromSysErrVector(readoutNeg, errorsNegCutSt9,errorsNegErrCutSt9,9, numberCutStudies ,offsetAtEnd ,nPtBins,decisionBoundary);
		CalculateMeanSysErr(errorsMeanCutSt9, errorsMeanErrCutSt9, errorsPosCutSt9, errorsNegCutSt9, nPtBins);
		
		CorrectSystematicErrorsWithMean(errorsPosCutSt9,errorsPosErrCutSt9, errorsPosCorrCutSt9, errorsPosErrCorrCutSt9, nPtBins);
		CorrectSystematicErrorsWithMean(errorsNegCutSt9,errorsNegErrCutSt9, errorsNegCorrCutSt9, errorsNegErrCorrCutSt9, nPtBins);
		CorrectSystematicErrorsWithMean(errorsMeanCutSt9, errorsMeanErrCutSt9, errorsMeanCorrCutSt9, errorsMeanErrCorrCutSt9, nPtBins);
		
		negativeErrors[8] = new TGraphErrors(nPtBins,ptBins ,errorsNegCutSt9 ,ptBinsErr ,errorsNegErrCutSt9 );
		meanErrors[8] = new TGraphErrors(nPtBins,ptBins ,errorsMeanCutSt9 ,ptBinsErr ,errorsMeanErrCutSt9 );
		positiveErrors[8] = new TGraphErrors(nPtBins,ptBins ,errorsPosCutSt9 ,ptBinsErr ,errorsPosErrCutSt9 );
		negativeErrorsCorr[8] = new TGraphErrors(nPtBins,ptBins ,errorsNegCorrCutSt9 ,ptBinsErr ,errorsNegErrCorrCutSt9 );
		meanErrorsCorr[8] = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrCutSt9 ,ptBinsErr ,errorsMeanErrCorrCutSt9 );
		positiveErrorsCorr[8] = new TGraphErrors(nPtBins,ptBins ,errorsPosCorrCutSt9 ,ptBinsErr ,errorsPosErrCorrCutSt9 );
		
		for (Int_t l = 0; l < nPtBins; l++){
			errorsPosSummed[l] = errorsPosSummed[l] + pow(errorsPosCutSt9[l],2);
			errorsNegSummed[l] = errorsNegSummed[l] + pow(errorsNegCutSt9[l],2);
			errorsMeanSummed[l] = errorsMeanSummed[l] + pow(errorsMeanCutSt9[l],2);
			errorsPosCorrSummed[l] = errorsPosCorrSummed[l] + pow(errorsPosCorrCutSt9[l],2);
			errorsNegCorrSummed[l] = errorsNegCorrSummed[l] + pow(errorsNegCorrCutSt9[l],2);
			errorsMeanCorrSummed[l] = errorsMeanCorrSummed[l] + pow(errorsMeanCorrCutSt9[l],2);
		}
	}
	if (numberCutStudies >= 10){	
		cout << "Cutstudies 10" << endl;
		ReadOutFromSysErrVector(readoutPos, errorsPosCutSt10,errorsPosErrCutSt10,10, numberCutStudies ,offsetAtEnd ,nPtBins,decisionBoundary);
		ReadOutFromSysErrVector(readoutNeg, errorsNegCutSt10,errorsNegErrCutSt10,10, numberCutStudies ,offsetAtEnd ,nPtBins,decisionBoundary);		
		CalculateMeanSysErr(errorsMeanCutSt10, errorsMeanErrCutSt10, errorsPosCutSt10, errorsNegCutSt10, nPtBins);
			
		CorrectSystematicErrorsWithMean(errorsPosCutSt10,errorsPosErrCutSt10, errorsPosCorrCutSt10, errorsPosErrCorrCutSt10, nPtBins);
		CorrectSystematicErrorsWithMean(errorsNegCutSt10,errorsNegErrCutSt10, errorsNegCorrCutSt10, errorsNegErrCorrCutSt10, nPtBins);
		CorrectSystematicErrorsWithMean(errorsMeanCutSt10, errorsMeanErrCutSt10, errorsMeanCorrCutSt10, errorsMeanErrCorrCutSt10, nPtBins);
		
		negativeErrors[9] = new TGraphErrors(nPtBins,ptBins ,errorsNegCutSt10 ,ptBinsErr ,errorsNegErrCutSt10 );
		meanErrors[9] = new TGraphErrors(nPtBins,ptBins ,errorsMeanCutSt10 ,ptBinsErr ,errorsMeanErrCutSt10 );
		positiveErrors[9] = new TGraphErrors(nPtBins,ptBins ,errorsPosCutSt10 ,ptBinsErr ,errorsPosErrCutSt10 );
		negativeErrorsCorr[9] = new TGraphErrors(nPtBins,ptBins ,errorsNegCorrCutSt10 ,ptBinsErr ,errorsNegErrCorrCutSt10 );
		meanErrorsCorr[9] = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrCutSt10 ,ptBinsErr ,errorsMeanErrCorrCutSt10 );
		positiveErrorsCorr[9] = new TGraphErrors(nPtBins,ptBins ,errorsPosCorrCutSt10 ,ptBinsErr ,errorsPosErrCorrCutSt10 );

		for (Int_t l = 0; l < nPtBins; l++){
			errorsPosSummed[l] = errorsPosSummed[l] + pow(errorsPosCutSt10[l],2);
			errorsNegSummed[l] = errorsNegSummed[l] + pow(errorsNegCutSt10[l],2);
			errorsMeanSummed[l] = errorsMeanSummed[l] + pow(errorsMeanCutSt10[l],2);
			errorsPosCorrSummed[l] = errorsPosCorrSummed[l] + pow(errorsPosCorrCutSt10[l],2);
			errorsNegCorrSummed[l] = errorsNegCorrSummed[l] + pow(errorsNegCorrCutSt10[l],2);
			errorsMeanCorrSummed[l] = errorsMeanCorrSummed[l] + pow(errorsMeanCorrCutSt10[l],2);
		}
	}	
	
	
	for (Int_t l = 0; l < nPtBins; l++){
			errorsPosSummed[l] = pow(errorsPosSummed[l],0.5);
			errorsMeanSummed[l] = pow(errorsMeanSummed[l],0.5);
			errorsPosErrSummed[l] = errorsPosSummed[l]*0.001;
			errorsMeanErrSummed[l] = errorsMeanSummed[l]*0.001;
			errorsNegSummed[l] = -pow(errorsNegSummed[l],0.5);
			errorsNegErrSummed[l] = errorsNegSummed[l]*0.001;
// 			errorsPosCorrMatSummed[l] = pow(errorsPosCorrSummed[l]+ pow(TMath::Sqrt(2)*6.3 ,2.),0.5);
// 			errorsMeanCorrMatSummed[l] = pow(errorsMeanCorrSummed[l]+ pow(TMath::Sqrt(2)*6.3 ,2.),0.5);
// 			errorsNegCorrMatSummed[l] = -pow(errorsNegCorrSummed[l]+ pow(TMath::Sqrt(2)*3.4 ,2.),0.5);
			errorsPosCorrMatSummed[l] = pow(errorsPosCorrSummed[l]+ pow(2*4.50 ,2.),0.5);
			errorsMeanCorrMatSummed[l] = pow(errorsMeanCorrSummed[l]+ pow(2*4.50 ,2.),0.5);
			errorsNegCorrMatSummed[l] = -pow(errorsNegCorrSummed[l]+ pow(2*4.50 ,2.),0.5);
			
			errorsPosCorrSummed[l] = pow(errorsPosCorrSummed[l],0.5);
			errorsMeanCorrSummed[l] = pow(errorsMeanCorrSummed[l],0.5);
			errorsPosErrCorrSummed[l] = errorsPosCorrSummed[l]*0.001;
			errorsMeanErrCorrSummed[l] = errorsMeanCorrSummed[l]*0.001;
			errorsMeanErrCorrMatSummed[l] = errorsMeanCorrMatSummed[l]*0.001;
			errorsNegCorrSummed[l] = -pow(errorsNegCorrSummed[l],0.5);
			errorsNegErrCorrSummed[l] = errorsNegCorrSummed[l]*0.001;
			cout << l <<errorsPosCorrSummed[l] << "\t" << errorsNegCorrSummed[l] << endl;	
	}
	negativeErrorsSummed = new TGraphErrors(nPtBins,ptBins ,errorsNegSummed ,ptBinsErr ,errorsNegErrSummed );
	negativeErrorsCorrSummed = new TGraphErrors(nPtBins,ptBins ,errorsNegCorrSummed ,ptBinsErr ,errorsNegErrCorrSummed );
	positiveErrorsSummed = new TGraphErrors(nPtBins,ptBins ,errorsPosSummed ,ptBinsErr ,errorsPosErrSummed );
	positiveErrorsCorrSummed = new TGraphErrors(nPtBins,ptBins ,errorsPosCorrSummed ,ptBinsErr ,errorsPosErrCorrSummed );
	meanErrorsSummed = new TGraphErrors(nPtBins,ptBins ,errorsMeanSummed ,ptBinsErr ,errorsMeanErrSummed );
	meanErrorsCorrSummed = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSummed ,ptBinsErr ,errorsMeanErrCorrSummed );
	meanErrorsCorrSummedIncMat = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrMatSummed ,ptBinsErr ,errorsMeanErrCorrMatSummed );
		
	negativeErrorsCorrSummed->Print();
	positiveErrorsCorrSummed->Print();
	meanErrorsCorrSummed->Print();
		//*********************** Comparison of Measuremnts in ALICE ***************************************************
	TCanvas* canvasSysErrNeg = new TCanvas("canvasSysErrNeg","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasSysErrNeg, 0.08, 0.01, 0.015, 0.09);
	TH2D *histo2DSysErrNeg ;
	
	
	if (energy.CompareTo("7TeV")==0 && meson.CompareTo("Pi0EtaBinning") == 0 ) {
		histo2DSysErrNeg = new TH2D("histo2DSysErrNeg", "histo2DSysErrNeg", 20,0.,ptBins[nPtBins-1]+1,1000.,-15.,0.5);
	} else if (energy.CompareTo("7TeV")==0 && meson.CompareTo("Pi0") == 0 ) {
		histo2DSysErrNeg = new TH2D("histo2DSysErrNeg", "histo2DSysErrNeg", 20,0.,ptBins[nPtBins-1]+1,1000.,-45.,0.5);
	} else if (energy.CompareTo("900GeV")==0 && meson.CompareTo("Pi0") == 0 ) { 
		histo2DSysErrNeg = new TH2D("histo2DSysErrNeg", "histo2DSysErrNeg", 20,0.,ptBins[nPtBins-1]+1,1000.,-15.,0.5);
	} else if (energy.CompareTo("900GeV")==0 && meson.CompareTo("Pi0EtaBinning") == 0 ) { 
		histo2DSysErrNeg = new TH2D("histo2DSysErrNeg", "histo2DSysErrNeg", 20,0.,ptBins[nPtBins-1]+1,1000.,-6.,0.5);
	} else if (energy.CompareTo("2.76TeV")==0 && meson.CompareTo("Pi0EtaBinning") == 0 ) { 
		histo2DSysErrNeg = new TH2D("histo2DSysErrNeg", "histo2DSysErrNeg", 20,0.,ptBins[nPtBins-1]+1,1000.,-6.,0.5);
	} else {	
		histo2DSysErrNeg = new TH2D("histo2DSysErrNeg", "histo2DSysErrNeg", 20,0.,ptBins[nPtBins-1]+1,1000.,-25.,0.5);
	}
	histo2DSysErrNeg->SetYTitle("negative systematic Err %");
	histo2DSysErrNeg->SetXTitle("p_{T} (GeV/c)");
	histo2DSysErrNeg->GetYaxis()->SetNdivisions(510,kTRUE);
	histo2DSysErrNeg->GetYaxis()->SetLabelSize(0.03);
	histo2DSysErrNeg->GetYaxis()->SetTitleSize(0.04);
	histo2DSysErrNeg->GetYaxis()->SetDecimals();
	histo2DSysErrNeg->GetYaxis()->SetTitleOffset(0.9);
	histo2DSysErrNeg->GetXaxis()->SetTitleOffset(1.);
	histo2DSysErrNeg->GetXaxis()->SetLabelSize(0.03);
	histo2DSysErrNeg->GetXaxis()->SetTitleSize(0.04);
	histo2DSysErrNeg->SetTitle("");
	histo2DSysErrNeg->Draw();
	
	TLegend* legendNegative = new TLegend(0.45,0.15,0.75,0.43);
	legendNegative->SetTextSize(0.035);
	legendNegative->SetFillColor(0);
	legendNegative->SetBorderSize(0);
	for(Int_t i = 0; i< numberCutStudies ; i++){
		if (energy.CompareTo("7TeV")==0 && i == 5){ 
			cout << "skipped Period contribution" << endl;	
		} else {
			DrawGammaSetMarkerTGraphErr(negativeErrors[i], 20+i, 1.,color[i],color[i]);
			negativeErrors[i]->Draw("p,csame");
			legendNegative->AddEntry(negativeErrors[i],nameCutVariation[i].Data(),"p");
		}
	}
	DrawGammaSetMarkerTGraphErr(negativeErrorsSummed, 20, 1.,1,1);
	negativeErrorsSummed->Draw("p,csame");
	legendNegative->AddEntry(negativeErrorsSummed,"quadratically summed","p");
	legendNegative->Draw();
	canvasSysErrNeg->Update();
	canvasSysErrNeg->SaveAs(Form("SystematicErrors/SysNeg_%s_%s.eps",meson.Data(), energy.Data()));
	
	delete canvasSysErrNeg;

	TCanvas* canvasSysErrPos = new TCanvas("canvasSysErrPos","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasSysErrPos, 0.08, 0.01, 0.015, 0.09);
	TH2D *histo2DSysErrPos ;
	
	
	if (energy.CompareTo("7TeV")==0 && meson.CompareTo("Pi0EtaBinning") == 0 ) {
		histo2DSysErrPos = new TH2D("histo2DSysErrPos", "histo2DSysErrPos", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,10.);
	} else if (energy.CompareTo("7TeV")==0 && meson.CompareTo("Pi0") == 0 ) {
		histo2DSysErrPos = new TH2D("histo2DSysErrPos", "histo2DSysErrPos", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,35.);
	} else if (energy.CompareTo("7TeV")==0 && meson.CompareTo("Eta") == 0 ) {
		histo2DSysErrPos = new TH2D("histo2DSysErrPos", "histo2DSysErrPos", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,30.);
	} else if (energy.CompareTo("900GeV")==0 && meson.CompareTo("Pi0") == 0 ) { 
		histo2DSysErrPos = new TH2D("histo2DSysErrPos", "histo2DSysErrPos", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,8.);
	} else if (energy.CompareTo("900GeV")==0 && meson.CompareTo("Pi0EtaBinning") == 0 ) { 
		histo2DSysErrPos = new TH2D("histo2DSysErrPos", "histo2DSysErrPos", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,8.);
	} else if (energy.CompareTo("2.76TeV")==0 && meson.CompareTo("Pi0EtaBinning") == 0 ) { 
		histo2DSysErrPos = new TH2D("histo2DSysErrPos", "histo2DSysErrPos", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,10.);
	} else if (energy.CompareTo("900GeV")==0 && meson.CompareTo("Eta") == 0 ) { 
		histo2DSysErrPos = new TH2D("histo2DSysErrPos", "histo2DSysErrPos", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,40.);
	} else {
		histo2DSysErrPos = new TH2D("histo2DSysErrPos", "histo2DSysErrPos", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,25.);
	}
	histo2DSysErrPos->SetYTitle("positive systematic Err %");
	histo2DSysErrPos->SetXTitle("p_{T} (GeV/c)");
	histo2DSysErrPos->GetYaxis()->SetNdivisions(510,kTRUE);
	histo2DSysErrPos->GetYaxis()->SetLabelSize(0.03);
	histo2DSysErrPos->GetYaxis()->SetTitleSize(0.04);
	histo2DSysErrPos->GetYaxis()->SetDecimals();
	histo2DSysErrPos->GetYaxis()->SetTitleOffset(0.9);
	histo2DSysErrPos->GetXaxis()->SetTitleOffset(1.);
	histo2DSysErrPos->GetXaxis()->SetLabelSize(0.03);
	histo2DSysErrPos->GetXaxis()->SetTitleSize(0.04);
	histo2DSysErrPos->SetTitle("");
	histo2DSysErrPos->Draw();
	
	TLegend* legendPositive = new TLegend(0.6,0.65,0.95,0.95);
	legendPositive->SetTextSize(0.035);
	legendPositive->SetFillColor(0);
	legendPositive->SetBorderSize(0);
	for(Int_t i = 0; i< numberCutStudies ; i++){
		if (energy.CompareTo("7TeV")==0 && i == 5){ 
			cout << "skipped Period contribution" << endl;
		} else {
			DrawGammaSetMarkerTGraphErr(positiveErrors[i], 20+i, 1.,color[i],color[i]);
			positiveErrors[i]->Draw("p,csame");
			legendPositive->AddEntry(positiveErrors[i],nameCutVariation[i].Data(),"p");
		}	
	}
	DrawGammaSetMarkerTGraphErr(positiveErrorsSummed, 20, 1.,1,1);
	positiveErrorsSummed->Draw("p,csame");
	legendPositive->AddEntry(positiveErrorsSummed,"quadratically summed","p");
	legendPositive->Draw();
	canvasSysErrPos->Update();
	canvasSysErrPos->SaveAs(Form("SystematicErrors/SysPos_%s_%s.eps",meson.Data(), energy.Data()));
	
	delete canvasSysErrPos;

	
	TCanvas* canvasNewSysErrPos = new TCanvas("canvasNewSysErrPos","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasNewSysErrPos, 0.08, 0.01, 0.015, 0.09);
	TH2D *histo2DNewSysErrPos ;
	
	
	if (energy.CompareTo("7TeV")==0 && meson.CompareTo("Pi0EtaBinning") == 0 ) {
		histo2DNewSysErrPos = new TH2D("histo2DNewSysErrPos", "histo2DNewSysErrPos", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,10.);
	} else if (energy.CompareTo("7TeV")==0 && meson.CompareTo("Pi0") == 0 ) {
		histo2DNewSysErrPos = new TH2D("histo2DNewSysErrPos", "histo2DNewSysErrPos", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,50.);
	} else if (energy.CompareTo("7TeV")==0 && meson.CompareTo("Eta") == 0 ) {
		histo2DNewSysErrPos = new TH2D("histo2DNewSysErrPos", "histo2DNewSysErrPos", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,30.);
	} else if (energy.CompareTo("900GeV")==0 && meson.CompareTo("Pi0") == 0 ) { 
		histo2DNewSysErrPos = new TH2D("histo2DNewSysErrPos", "histo2DNewSysErrPos", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,8.);
	} else if (energy.CompareTo("900GeV")==0 && meson.CompareTo("Pi0EtaBinning") == 0 ) { 
		histo2DNewSysErrPos = new TH2D("histo2DNewSysErrPos", "histo2DNewSysErrPos", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,8.);
	} else if (energy.CompareTo("2.76TeV")==0 && meson.CompareTo("Pi0EtaBinning") == 0 ) { 
		histo2DNewSysErrPos = new TH2D("histo2DNewSysErrPos", "histo2DNewSysErrPos", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,10.);
	} else if (energy.CompareTo("900GeV")==0 && meson.CompareTo("Eta") == 0 ) { 
		histo2DNewSysErrPos = new TH2D("histo2DNewSysErrPos", "histo2DNewSysErrPos", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,40.);
	} else {
		histo2DNewSysErrPos = new TH2D("histo2DNewSysErrPos", "histo2DNewSysErrPos", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,25.);
	}
	histo2DNewSysErrPos->SetYTitle("positive systematic Err %");
	histo2DNewSysErrPos->SetXTitle("p_{T} (GeV/c)");
	histo2DNewSysErrPos->GetYaxis()->SetNdivisions(510,kTRUE);
	histo2DNewSysErrPos->GetYaxis()->SetLabelSize(0.03);
	histo2DNewSysErrPos->GetYaxis()->SetTitleSize(0.04);
	histo2DNewSysErrPos->GetYaxis()->SetDecimals();
	histo2DNewSysErrPos->GetYaxis()->SetTitleOffset(0.9);
	histo2DNewSysErrPos->GetXaxis()->SetTitleOffset(1.);
	histo2DNewSysErrPos->GetXaxis()->SetLabelSize(0.03);
	histo2DNewSysErrPos->GetXaxis()->SetTitleSize(0.04);
	histo2DNewSysErrPos->SetTitle("");
	histo2DNewSysErrPos->Draw();
	
	TLegend* legendPositiveNew = new TLegend(0.6,0.65,0.95,0.95);
	legendPositiveNew->SetTextSize(0.035);
	legendPositiveNew->SetFillColor(0);
	legendPositiveNew->SetBorderSize(0);
	for(Int_t i = 0; i< numberCutStudies ; i++){
		if (energy.CompareTo("7TeV")==0 && i == 5){ 
			cout << "skipped Period contribution" << endl;
		} else {
			DrawGammaSetMarkerTGraphErr(positiveErrorsCorr[i], 20+i, 1.,color[i],color[i]);
			positiveErrorsCorr[i]->Draw("p,csame");
			legendPositiveNew->AddEntry(positiveErrorsCorr[i],nameCutVariation[i].Data(),"p");
		}	
	}
	DrawGammaSetMarkerTGraphErr(positiveErrorsCorrSummed, 20, 1.,1,1);
	positiveErrorsCorrSummed->Draw("p,csame");
	legendPositiveNew->AddEntry(positiveErrorsCorrSummed,"quadratically summed","p");
	legendPositiveNew->Draw();
	canvasNewSysErrPos->Update();
	canvasNewSysErrPos->SaveAs(Form("SystematicErrors/SysPosNewWithMean_%s_%s.eps",meson.Data(), energy.Data()));
	
	delete canvasNewSysErrPos;

	TCanvas* canvasNewSysErrNeg = new TCanvas("canvasNewSysErrNeg","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasNewSysErrNeg, 0.08, 0.01, 0.015, 0.09);
	TH2D *histo2DNewSysErrNeg ;
	
	
	if (energy.CompareTo("7TeV")==0 && meson.CompareTo("Pi0EtaBinning") == 0 ) {
		histo2DNewSysErrNeg = new TH2D("histo2DNewSysErrNeg", "histo2DNewSysErrNeg", 20,0.,ptBins[nPtBins-1]+1,1000.,-15.,0.5);
	} else if (energy.CompareTo("7TeV")==0 && meson.CompareTo("Pi0") == 0 ) {
		histo2DNewSysErrNeg = new TH2D("histo2DNewSysErrNeg", "histo2DNewSysErrNeg", 20,0.,ptBins[nPtBins-1]+1,1000.,-45.,0.5);
	} else if (energy.CompareTo("900GeV")==0 && meson.CompareTo("Pi0") == 0 ) { 
		histo2DNewSysErrNeg = new TH2D("histo2DNewSysErrNeg", "histo2DNewSysErrNeg", 20,0.,ptBins[nPtBins-1]+1,1000.,-15.,0.5);
	} else if (energy.CompareTo("900GeV")==0 && meson.CompareTo("Pi0EtaBinning") == 0 ) { 
		histo2DNewSysErrNeg = new TH2D("histo2DNewSysErrNeg", "histo2DNewSysErrNeg", 20,0.,ptBins[nPtBins-1]+1,1000.,-6.,0.5);
	} else if (energy.CompareTo("2.76TeV")==0 && meson.CompareTo("Pi0EtaBinning") == 0 ) { 
		histo2DNewSysErrNeg = new TH2D("histo2DNewSysErrNeg", "histo2DNewSysErrNeg", 20,0.,ptBins[nPtBins-1]+1,1000.,-6.,0.5);
	} else {	
		histo2DNewSysErrNeg = new TH2D("histo2DNewSysErrNeg", "histo2DNewSysErrNeg", 20,0.,ptBins[nPtBins-1]+1,1000.,-25.,0.5);
	}
	histo2DNewSysErrNeg->SetYTitle("negative systematic Err %");
	histo2DNewSysErrNeg->SetXTitle("p_{T} (GeV/c)");
	histo2DNewSysErrNeg->GetYaxis()->SetNdivisions(510,kTRUE);
	histo2DNewSysErrNeg->GetYaxis()->SetLabelSize(0.03);
	histo2DNewSysErrNeg->GetYaxis()->SetTitleSize(0.04);
	histo2DNewSysErrNeg->GetYaxis()->SetDecimals();
	histo2DNewSysErrNeg->GetYaxis()->SetTitleOffset(0.9);
	histo2DNewSysErrNeg->GetXaxis()->SetTitleOffset(1.);
	histo2DNewSysErrNeg->GetXaxis()->SetLabelSize(0.03);
	histo2DNewSysErrNeg->GetXaxis()->SetTitleSize(0.04);
	histo2DNewSysErrNeg->SetTitle("");
	histo2DNewSysErrNeg->Draw();
	
	TLegend* legendNegativeNew = new TLegend(0.45,0.15,0.75,0.43);
	legendNegativeNew->SetTextSize(0.035);
	legendNegativeNew->SetFillColor(0);
	legendNegativeNew->SetBorderSize(0);
	for(Int_t i = 0; i< numberCutStudies ; i++){
		if (energy.CompareTo("7TeV")==0 && i == 5){ 
			cout << "skipped Period contribution" << endl;	
		} else {
			DrawGammaSetMarkerTGraphErr(negativeErrorsCorr[i], 20+i, 1.,color[i],color[i]);
			negativeErrorsCorr[i]->Draw("p,csame");
			legendNegativeNew->AddEntry(negativeErrorsCorr[i],nameCutVariation[i].Data(),"p");
		}
	}
	DrawGammaSetMarkerTGraphErr(negativeErrorsCorrSummed, 20, 1.,1,1);
	negativeErrorsCorrSummed->Draw("p,csame");
	legendNegativeNew->AddEntry(negativeErrorsCorrSummed,"quadratically summed","p");
	legendNegativeNew->Draw();
	canvasNewSysErrNeg->Update();
	canvasNewSysErrNeg->SaveAs(Form("SystematicErrors/SysNegNewWithMean_%s_%s.eps",meson.Data(), energy.Data()));
	
	delete canvasNewSysErrNeg;

	
	TCanvas* canvasSysErrMean = new TCanvas("canvasSysErrMean","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasSysErrMean, 0.08, 0.01, 0.015, 0.09);
	TH2D *histo2DSysErrMean ;
	
	
	if (energy.CompareTo("7TeV")==0 && meson.CompareTo("Pi0EtaBinning") == 0 ) {
		histo2DSysErrMean = new TH2D("histo2DSysErrMean", "histo2DSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,0.,10.);
	} else if (energy.CompareTo("7TeV")==0 && meson.CompareTo("Pi0") == 0 ) {
		histo2DSysErrMean = new TH2D("histo2DSysErrMean", "histo2DSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,0.,25.);
	} else if (energy.CompareTo("7TeV")==0 && meson.CompareTo("Eta") == 0 ) {
		histo2DSysErrMean = new TH2D("histo2DSysErrMean", "histo2DSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,0.,30.);
	} else if (energy.CompareTo("900GeV")==0 && meson.CompareTo("Pi0") == 0 ) { 
		histo2DSysErrMean = new TH2D("histo2DSysErrMean", "histo2DSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,0.,8.);
	} else if (energy.CompareTo("900GeV")==0 && meson.CompareTo("Pi0EtaBinning") == 0 ) { 
		histo2DSysErrMean = new TH2D("histo2DSysErrMean", "histo2DSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,0.,8.);
	} else if (energy.CompareTo("2.76TeV")==0 && meson.CompareTo("Pi0EtaBinning") == 0 ) { 
		histo2DSysErrMean = new TH2D("histo2DSysErrMean", "histo2DSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,0.,10.);
	} else if (energy.CompareTo("900GeV")==0 && meson.CompareTo("Eta") == 0 ) { 
		histo2DSysErrMean = new TH2D("histo2DSysErrMean", "histo2DSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,0.,40.);
	} else {
		histo2DSysErrMean = new TH2D("histo2DSysErrMean", "histo2DSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,0.,25.);
	}
	histo2DSysErrMean->SetYTitle("mean systematic Err %");
	histo2DSysErrMean->SetXTitle("p_{T} (GeV/c)");
	histo2DSysErrMean->GetYaxis()->SetNdivisions(510,kTRUE);
	histo2DSysErrMean->GetYaxis()->SetLabelSize(0.03);
	histo2DSysErrMean->GetYaxis()->SetTitleSize(0.04);
	histo2DSysErrMean->GetYaxis()->SetDecimals();
	histo2DSysErrMean->GetYaxis()->SetTitleOffset(0.9);
	histo2DSysErrMean->GetXaxis()->SetTitleOffset(1.);
	histo2DSysErrMean->GetXaxis()->SetLabelSize(0.03);
	histo2DSysErrMean->GetXaxis()->SetTitleSize(0.04);
	histo2DSysErrMean->SetTitle("");
	histo2DSysErrMean->Draw();
	
	TLegend* legendMean = new TLegend(0.13,0.65,0.35,0.95);
	legendMean->SetTextSize(0.035);
	legendMean->SetFillColor(0);
	legendMean->SetBorderSize(0);
	for(Int_t i = 0; i< numberCutStudies ; i++){
		if (energy.CompareTo("7TeV")==0 && i == 5){ 
			cout << "skipped Period contribution" << endl;
		} else {
			DrawGammaSetMarkerTGraphErr(meanErrors[i], 20+i, 1.,color[i],color[i]);
			meanErrors[i]->Draw("p,csame");
			legendMean->AddEntry(meanErrors[i],nameCutVariation[i].Data(),"p");
		}	
	}
	DrawGammaSetMarkerTGraphErr(meanErrorsSummed, 20, 1.,1,1);
	meanErrorsSummed->Draw("p,csame");
	legendMean->AddEntry(meanErrorsSummed,"quadratically summed","p");
	legendMean->Draw();
	canvasSysErrMean->Update();
	canvasSysErrMean->SaveAs(Form("SystematicErrors/SysMean_%s_%s.eps",meson.Data(), energy.Data()));
	
	delete canvasSysErrMean;

	if (energy.CompareTo("7TeV")==0 && meson.CompareTo("Pi0") == 0 ) {
		TCanvas* canvasSysErrMeanSingle = new TCanvas("canvasSysErrMeanSingle","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasSysErrMeanSingle, 0.08, 0.01, 0.015, 0.09);
		

		for(Int_t i = 0; i< numberCutStudies ; i++){
			if ( i == 5){ 
				cout << "skipped Period contribution" << endl;
			} else {
				histo2DSysErrMean->SetMaximum(20);
				histo2DSysErrMean->Draw();
				meanErrors[i]->Draw("p,csame");
				canvasSysErrMeanSingle->Update();
				canvasSysErrMeanSingle->SaveAs(Form("SystematicErrors/SysMean_%sCutVariation_%s_%s.eps",nameCutForOutputVariation7TeV[i].Data(),meson.Data(), energy.Data()));
			}	
		}
		delete canvasSysErrMean;
	}
	
	TCanvas* canvasNewSysErrMean = new TCanvas("canvasNewSysErrMean","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasNewSysErrMean, 0.08, 0.01, 0.015, 0.09);
	TH2D *histo2DNewSysErrMean ;
	
	
	if (energy.CompareTo("7TeV")==0 && meson.CompareTo("Pi0EtaBinning") == 0 ) {
		histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "histo2DNewSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,10.);
	} else if (energy.CompareTo("7TeV")==0 && meson.CompareTo("Pi0") == 0 ) {
		histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "histo2DNewSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,25.);
	} else if (energy.CompareTo("7TeV")==0 && meson.CompareTo("Eta") == 0 ) {
		histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "histo2DNewSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,30.);
	} else if (energy.CompareTo("900GeV")==0 && meson.CompareTo("Pi0") == 0 ) { 
		histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "histo2DNewSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,25.);
	} else if (energy.CompareTo("900GeV")==0 && meson.CompareTo("Pi0EtaBinning") == 0 ) { 
		histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "histo2DNewSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,25.);
	} else if (energy.CompareTo("2.76TeV")==0 && meson.CompareTo("Pi0EtaBinning") == 0 ) { 
		histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "histo2DNewSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,10.);
	} else if (energy.CompareTo("900GeV")==0 && meson.CompareTo("Eta") == 0 ) { 
		histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "histo2DNewSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,40.);
	} else {
		histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "histo2DNewSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,25.);
	}
	histo2DNewSysErrMean->SetYTitle("positive systematic Err %");
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
	
	TLegend* legendMeanNew = new TLegend(0.15,0.65,0.35,0.95);
	legendMeanNew->SetTextSize(0.035);
	legendMeanNew->SetFillColor(0);
	legendMeanNew->SetBorderSize(0);
	for(Int_t i = 0; i< numberCutStudies ; i++){
		if (energy.CompareTo("7TeV")==0 && i == 5){ 
			cout << "skipped Period contribution" << endl;
		} else {
			DrawGammaSetMarkerTGraphErr(meanErrorsCorr[i], 20+i, 1.,color[i],color[i]);
			meanErrorsCorr[i]->Draw("p,csame");
			legendMeanNew->AddEntry(meanErrorsCorr[i],nameCutVariation[i].Data(),"p");
		}	
	}
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummed, 20, 1.,1,1);
	meanErrorsCorrSummed->Draw("p,csame");
	legendMeanNew->AddEntry(meanErrorsCorrSummed,"quadratically summed","p");
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kRed,kRed);
	meanErrorsCorrSummedIncMat->Draw("p,csame");
	legendMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"quadratically summed, inc mat","p");
	legendMeanNew->Draw();
	
	canvasNewSysErrMean->Update();
	canvasNewSysErrMean->SaveAs(Form("SystematicErrors/SysMeanNewWithMean_%s_%s.eps",meson.Data(), energy.Data()));
	
	delete canvasNewSysErrMean;

	if (energy.CompareTo("7TeV")==0 && meson.CompareTo("Pi0") == 0 ) {
		TCanvas* canvasSysErrMeanSingle = new TCanvas("canvasSysErrMeanSingle","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasSysErrMeanSingle, 0.08, 0.01, 0.015, 0.09);
		

		for(Int_t i = 0; i< numberCutStudies ; i++){
			if ( i == 5){ 
				cout << "skipped Period contribution" << endl;
			} else {
				histo2DSysErrMean->Draw();
				meanErrorsCorr[i]->Draw("p,csame");
				canvasSysErrMeanSingle->Update();
				canvasSysErrMeanSingle->SaveAs(Form("SystematicErrors/SysMeanSmoothed_%sCutVariation_%s_%s.eps",nameCutForOutputVariation7TeV[i].Data(),meson.Data(), energy.Data()));
			}	
		}
		delete canvasSysErrMean;
	}
	
		const char *SysErrDatname = Form("SystematicErrors/SystematicError_%s_%s_%s.dat",meson.Data(),energy.Data(),cStamp1);
		fstream SysErrDat;
		cout << SysErrDatname << endl;
		SysErrDat.open(SysErrDatname, ios::out);
		for (Int_t l=0; l< nPtBins; l++){
			SysErrDat <<errorsNegCorrMatSummed[l] << "\t" <<errorsPosCorrMatSummed[l] << "\t"  <<errorsNegCorrSummed[l] << "\t" <<errorsPosCorrSummed[l]  << endl;
		}

		SysErrDat.close();

		const char *SysErrDatnameMean = Form("SystematicErrors/SystematicErrorAveraged_%s_%s_%s.dat",meson.Data(),energy.Data(),cStamp1);
		fstream SysErrDatAver;
		cout << SysErrDatnameMean << endl;
		SysErrDatAver.open(SysErrDatnameMean, ios::out);
		for (Int_t l=0; l< nPtBins; l++){
			SysErrDatAver << "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
		}

		SysErrDatAver.close();

		const char *SysErrDatnameMeanPaper = Form("SystematicErrors/SystematicErrorAveragedPaper_%s_%s_%s.dat",meson.Data(),energy.Data(),cStamp1);
		fstream SysErrDatAverPaper;
		cout << SysErrDatnameMeanPaper << endl;
		SysErrDatAverPaper.open(SysErrDatnameMeanPaper, ios::out);
		for (Int_t l=0; l< nPtBins; l++){
// "Cluster variation"-1, "Single E/P p_{t} variation"-2,"dE/dx variation"-3, "#chi^{2} #gamma variation"-4, "#alpha variation"-5,"Yield extraction"-6,"BG variation"-6,"","",""
			SysErrDatAverPaper << "p_{T}" << "\t" << ptBins[l] << "\t" << errorsMeanCorrCutSt7[l]<<"+"<< errorsMeanCorrCutSt8[l]<<"="<< TMath::Sqrt(errorsMeanCorrCutSt7[l]*errorsMeanCorrCutSt7[l]+ errorsMeanCorrCutSt8[l]*errorsMeanCorrCutSt8[l]) << "\t" <<errorsMeanCorrCutSt3[l] << "\t" <<errorsMeanCorrCutSt4[l] << "\t" <<errorsMeanCorrCutSt1[l]<<"+"<< errorsMeanCorrCutSt2[l]<<"+" <<errorsMeanCorrCutSt5[l]<< "="<< TMath::Sqrt(errorsMeanCorrCutSt1[l]*errorsMeanCorrCutSt1[l]+errorsMeanCorrCutSt2[l]*errorsMeanCorrCutSt2[l]+errorsMeanCorrCutSt5[l]*errorsMeanCorrCutSt5[l])<<"\t" << endl;
		}

		SysErrDatAverPaper.close();

		
		
		const char *SysErrDatname2 = Form("SystematicErrors/SystematicErrorAllComp_%s_%s_%s.dat",meson.Data(),energy.Data(),cStamp1);
		fstream SysErrDat2;
		cout << SysErrDatname2 << endl;
		SysErrDat2.open(SysErrDatname2, ios::out);
		SysErrDat2 << "p_{T}" << "\t" 		 ;
		if (numberCutStudies >= 1){
			SysErrDat2 << nameCutVariation[0].Data() << "\t" 		 ;
		}
		if (numberCutStudies >= 2){
			SysErrDat2 << nameCutVariation[1].Data() << "\t"  ;
		}
		if (numberCutStudies >= 3){	
			SysErrDat2 << nameCutVariation[2].Data() << "\t" ;
		}
		if (numberCutStudies >= 4){
			SysErrDat2 << nameCutVariation[3].Data() << "\t" ;
		}
		if (numberCutStudies >= 5){	
			SysErrDat2 << nameCutVariation[4].Data() << "\t" ;
		}
		if (numberCutStudies >= 6){	
			if (energy.CompareTo("7TeV")==0){ 
				cout << "skipped Period contribution" << endl;	
			} else {
			SysErrDat2 << nameCutVariation[5].Data() << "\t" ;
			}
		}
		if (numberCutStudies >= 7){	
			SysErrDat2 << nameCutVariation[6].Data() << "\t" ;
		}
		if (numberCutStudies >= 8){	
			SysErrDat2 << nameCutVariation[7].Data() << "\t" ;
		}
		if (numberCutStudies >= 9){
			SysErrDat2 << nameCutVariation[8].Data() << "\t" ;
		}
		if (numberCutStudies >= 10){	
			SysErrDat2 << nameCutVariation[9].Data() << "\t" ;
		}	
		SysErrDat2 << "MatErr" <<"\t"<< "Summed Err w Mat"<<  "\t" <<"Summed Err w/o Mat"  << endl;

		SysErrDat2 << "Negative Table detailed results" << endl;
		for (Int_t l=0; l< nPtBins; l++){
			SysErrDat2 << ptBins[l] << "\t" ;
			if (numberCutStudies >= 1){
				SysErrDat2 << errorsNegCorrCutSt1[l] << "\t" 	;	 
			}
			if (numberCutStudies >= 2){
				SysErrDat2 << errorsNegCorrCutSt2[l] << "\t"  ;
			}
			if (numberCutStudies >= 3){	
				SysErrDat2 << errorsNegCorrCutSt3[l] << "\t" ;
			}
			if (numberCutStudies >= 4){
				SysErrDat2 << errorsNegCorrCutSt4[l] << "\t" ;
			}
			if (numberCutStudies >= 5){	
				SysErrDat2 << errorsNegCorrCutSt5[l] << "\t" ;
			}
			if (numberCutStudies >= 6){	
				if (energy.CompareTo("7TeV")==0){ 
					cout << "skipped Period contribution" << endl;	
				} else {
					SysErrDat2 << errorsNegCorrCutSt6[l] << "\t" ;
				}
			}
			if (numberCutStudies >= 7){	
				SysErrDat2 << errorsNegCorrCutSt7[l] << "\t" ;
			}
			if (numberCutStudies >= 8){	
				SysErrDat2 << errorsNegCorrCutSt8[l] << "\t" ;
			}
			if (numberCutStudies >= 9){
				SysErrDat2 << errorsNegCorrCutSt9[l] << "\t" ;
			}
			if (numberCutStudies >= 10){	
				SysErrDat2 << errorsNegCorrCutSt10[l] << "\t" ;
			}	
			SysErrDat2 << TMath::Sqrt(2)*6.2  <<"\t"<<  errorsNegCorrMatSummed[l] <<  "\t" <<errorsNegCorrSummed[l]  << endl;
		}
		SysErrDat2 << endl;
		SysErrDat2 << endl;
		SysErrDat2 << "Positive Table detailed results" << endl;
		for (Int_t l=0; l< nPtBins; l++){
			SysErrDat2 << ptBins[l] << "\t" ;
			if (numberCutStudies >= 1){
				SysErrDat2 << errorsPosCorrCutSt1[l] << "\t" ;		 
			}
			if (numberCutStudies >= 2){
				SysErrDat2 << errorsPosCorrCutSt2[l] << "\t"  ;
			}
			if (numberCutStudies >= 3){	
				SysErrDat2 << errorsPosCorrCutSt3[l] << "\t" ;
			}
			if (numberCutStudies >= 4){
				SysErrDat2 << errorsPosCorrCutSt4[l] << "\t" ;
			}
			if (numberCutStudies >= 5){	
				SysErrDat2 << errorsPosCorrCutSt5[l] << "\t" ;
			}
			if (numberCutStudies >= 6){	
				if (energy.CompareTo("7TeV")==0){ 
					cout << "skipped Period contribution" << endl;	
				} else {
					SysErrDat2 << errorsPosCorrCutSt6[l] << "\t" ;
				}
			}
			if (numberCutStudies >= 7){	
				SysErrDat2 << errorsPosCorrCutSt7[l] << "\t" ;
			}
			if (numberCutStudies >= 8){	
				SysErrDat2 << errorsPosCorrCutSt8[l] << "\t" ;
			}
			if (numberCutStudies >= 9){
				SysErrDat2 << errorsPosCorrCutSt9[l] << "\t" ;
			}
			if (numberCutStudies >= 10){	
				SysErrDat2 << errorsPosCorrCutSt10[l] << "\t" ;
			}	
			SysErrDat2 <<TMath::Sqrt(2)*6.2 << "\t" <<errorsPosCorrMatSummed[l] <<  "\t" <<errorsPosCorrSummed[l]  << endl;
		}

		SysErrDat2.close();
	

}