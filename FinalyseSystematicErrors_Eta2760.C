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

void FinalyseSystematicErrors_Eta2760(const char* nameDataFileErrors ="", TString energy="", TString meson = "", Int_t numberOfPtBins =1 ,Int_t numberCutStudies=1, Double_t decisionBoundary=1., Int_t offSetBeginning = 0, TString additionalName = "",TString additionalNameOutput = ""){
	
	StyleSettingsThesis();	
	SetPlotStyle();
	
	TString date = ReturnDateString();
	TString dateForOutput = ReturnDateStringForOutput();
	
	Int_t color[20] = {860,894,807,880,418,403,802,923,634,432,404,435,420,407,416,830,404,608,920,1};
	
	Double_t readoutPos[1000];
	Double_t readoutNeg[1000];
	Int_t 	numberOfEntriesPos = 0;
	Int_t 	numberOfEntriesNeg = 0;
	
	TFile* fileErrorInput= new TFile(nameDataFileErrors);
	const Int_t nPtBins = numberOfPtBins;
	const Int_t nCuts = numberCutStudies;
	Double_t* ptBins;
	Double_t* ptBinsErr;
	
	TString nameCutVariation7TeV[14] = {"Yield extraction","BG variation","dE/dx e-line","dE/dx #pi-line","TPC cluster","TOF variation","Single e^{#pm} p_{t}",  "#chi^{2} #gamma","q_{t,max}","R", "V0-finder","#alpha","#psi_{pair}", "#eta"};  
	TString nameCutForOutputVariation7TeV[14] = {"Yield extraction","BG variation","dE/dx e-line","dE/dx #pi-line","TPC cluster","TOF variation","Single e^{#pm} p_{t}",  "#chi^{2} #gamma","q_{t,max}", "V0-finder","R","#alpha","#psi_{pair}", "#eta"};  

	TString nameCutVariation900GeV[14] = {"Yield extraction","dE/dx e-line","dE/dx #pi-line","TPC cluster","TOF variation","Single e^{#pm} p_{t}",  "#chi^{2} #gamma","q_{t,max}","#alpha","R","#psi_{pair}", "V0-finder", "#eta","BG variation"};  
	TString nameCutVariationPbPb2760GeV[15] = {"Yield extraction","dE/dx e-line","dE/dx #pi-line","TPC cluster","TOF variation","Single e^{#pm} p_{t}",  "#chi^{2} #gamma","q_{t,max}","#alpha", "#psi_{pair}", "#eta","Pileup Estimate" ,"V0-finder","R", "BG variation"};  
	TString nameCutVariation2760GeV[9] = {"Yield extraction","BG variation", "TPC cluster","Single e^{#pm} p_{t}","#chi^{2} #gamma","q_{t,max}","#alpha","BGEstimate","Pileup Estimate"};  

	TString nameCutVariation[16];
	TString nameCutVariationSC[16];
	
	TString nameCutVariationSCPbPb2760GeV[15] = {"YieldExtraction","dEdxE","dEdxPi", "Cluster", "TOF", "SinglePt", "Chi2", "Qt", "Alpha", "PsiPair", "Eta","BGEstimate","V0Offline", "R","BG"};
	
	TString nameCutVariationSC2760GeV[8] = {"YieldExtraction","BG", "TPCCluster", "SinglePt", "Chi2", "Qt", "Alpha","BGEstimate"};
	
	
// 	if (energy.CompareTo("7TeV") == 0) {
// 		for (Int_t i = 0; i < numberCutStudies; i++){
// 			nameCutVariation[i] = nameCutVariation7TeV[i];
// 			nameCutVariationSC[i] = nameCutVariationSC2760GeV[i];
// 		}
// 	} 
	if (energy.CompareTo("2.76TeV") == 0) { 
		for (Int_t i = 0; i < numberCutStudies; i++){
			nameCutVariation[i] = nameCutVariation2760GeV[i];
			nameCutVariationSC[i] = nameCutVariationSC2760GeV[i];
			cout<< "Defining names"<<  i<< "\t" <<nameCutVariation[i].Data() << "\t" << nameCutVariationSC[i].Data()<<  endl;
		}	
	} else if (energy.CompareTo("PbPb_2.76TeV") == 0) {
		for (Int_t i = 0; i < numberCutStudies; i++){
			nameCutVariation[i] = nameCutVariationPbPb2760GeV[i];
			nameCutVariationSC[i] = nameCutVariationSCPbPb2760GeV[i];
		}
	}
// 	else { 
// 		for (Int_t i = 0; i < numberCutStudies; i++){
// 			nameCutVariation[i] = nameCutVariation900GeV[i];
// 			nameCutVariationSC[i] = nameCutVariationSC2760GeV[i];
// 		}
// 	}
	
	Double_t* errorsNeg[nCuts];
	Double_t errorsNegCorr[nCuts][nPtBins];
	Double_t errorsNegSummed[nPtBins];
	Double_t errorsNegCorrSummed[nPtBins];
	Double_t errorsNegCorrMatSummed[nPtBins];
	
	Double_t* errorsNegErr[nCuts];
	Double_t errorsNegErrCorr[nCuts][nPtBins];
	Double_t errorsNegErrSummed[nPtBins];
	Double_t errorsNegErrCorrSummed[nPtBins];
	
	Double_t* errorsPos[nCuts];
	Double_t errorsPosCorr[nCuts][nPtBins];
	Double_t errorsPosSummed[nPtBins];
	Double_t errorsPosCorrSummed[nPtBins];
	Double_t errorsPosCorrMatSummed[nPtBins];
	
	Double_t* errorsPosErr[nCuts];
	Double_t errorsPosErrSummed[nPtBins];
	Double_t errorsPosErrCorr[nCuts][nPtBins];
	Double_t errorsPosErrCorrSummed[nPtBins];
	
	Double_t errorsMean[nCuts][nPtBins];
	Double_t errorsMeanCorr[nCuts][nPtBins];
	Double_t errorsMeanSummed[nPtBins];
	Double_t errorsMeanCorrSummed[nPtBins];
	Double_t errorsMeanCorrMatSummed[nPtBins];
	
	Double_t errorsMeanErr[nCuts][nPtBins];
	Double_t errorsMeanErrCorr[nCuts][nPtBins];
	Double_t errorsMeanErrSummed[nPtBins];
	Double_t errorsMeanErrCorrSummed[nPtBins];
	Double_t errorsMeanErrCorrMatSummed[nPtBins];
	
	TGraphErrors* negativeErrors[nCuts];
	TGraphErrors* negativeErrorsSummed;
	TGraphErrors* positiveErrors[nCuts];
	TGraphErrors* positiveErrorsSummed;
	TGraphErrors* negativeErrorsCorr[nCuts];
	TGraphErrors* negativeErrorsCorrSummed;
	TGraphErrors* positiveErrorsCorr[nCuts];
	TGraphErrors* positiveErrorsCorrSummed;
	TGraphErrors* meanErrors[nCuts];
	TGraphErrors* meanErrorsSummed;
	TGraphErrors* meanErrorsCorr[nCuts];
	TGraphErrors* meanErrorsCorrSummed;
	TGraphErrors* meanErrorsCorrSummedIncMat;
	
	for (Int_t l = 0; l < nPtBins; l++){
			errorsPosSummed[l] = 0.;
			errorsNegSummed[l] = 0.;
			errorsMeanSummed[l] = 0.;
			errorsPosCorrSummed[l] = 0.;
			errorsNegCorrSummed[l] = 0.;
			errorsMeanCorrSummed[l] = 0.;
	} 
	
	for (Int_t i = 0; i < nCuts; i++){
	  TGraphAsymmErrors* graphPosErrors;
	  TGraphAsymmErrors* graphNegErrors;
	  if (i == 0){
	    TString nameGraphPos;
	    TString nameGraphNeg;
	    if (energy.CompareTo("PbPb_2.76TeV") == 0) {
	      nameGraphPos = Form("%s_SystErrorRelPos_%s_%s",meson.Data(),nameCutVariationSC[i].Data(),additionalName.Data()	);
	      nameGraphNeg = Form("%s_SystErrorRelNeg_%s_%s",meson.Data(),nameCutVariationSC[i].Data(),additionalName.Data()	);
	    } else {	
	      nameGraphPos = Form("%s_SystErrorRelPos_%s_pp",meson.Data(),nameCutVariationSC[i].Data());
	      nameGraphNeg = Form("%s_SystErrorRelNeg_%s_pp",meson.Data(),nameCutVariationSC[i].Data());
	    }
	    cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
	    
	    graphPosErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
	    graphNegErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
	  } else if (i == 10 && energy.CompareTo("2.76TeV") == 0)  {
	    TString nameGraphPos;
	    TString nameGraphNeg;
	    nameGraphPos = Form("%s_SystErrorRel_%s_pp",meson.Data(),nameCutVariationSC[i].Data(),additionalName.Data() );
	    nameGraphNeg = Form("%s_SystErrorRel_%s_pp",meson.Data(),nameCutVariationSC[i].Data(),additionalName.Data() );
	    cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
	    graphPosErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
	    graphNegErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
	  } else if ( energy.CompareTo("PbPb_2.76TeV") == 0 && nameCutVariationSC[i].CompareTo("BGEstimate")==0  ){
	    cout << "including pileup for PbPb" << endl;
	    TString nameGraphPos;
	    TString nameGraphNeg;
	    nameGraphPos = Form("%s_SystErrorRelPos_%s%s",meson.Data(),"dEdxE",additionalNameOutput.Data()  );
	    nameGraphNeg = Form("%s_SystErrorRelNeg_%s%s",meson.Data(),"dEdxE",additionalNameOutput.Data()  );
	    graphPosErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
	    graphNegErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
	  } else if ( energy.CompareTo("2.76TeV") == 0 && nameCutVariationSC[i].CompareTo("BGEstimate")==0  ){
	    
	    TString nameGraphPos;
	    TString nameGraphNeg;
	    nameGraphPos = Form("%s_SystErrorRel_%s_pp",meson.Data(),nameCutVariationSC[i].Data()  );
	    nameGraphNeg = Form("%s_SystErrorRel_%s_pp",meson.Data(),nameCutVariationSC[i].Data()  );
	    cout << "including pileup for pp2760" << nameGraphPos.Data()<< " "<< nameGraphNeg.Data()<< endl;
	    graphPosErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
	    graphNegErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());

	  } else {
	    TString nameGraphPos;
	    TString nameGraphNeg;
	    if (energy.CompareTo("PbPb_2.76TeV") == 0) {
	      if (nameCutVariationSC[i].CompareTo("Alpha")==0 && 
		  (additionalNameOutput.CompareTo("0005") || additionalNameOutput.CompareTo("0510") || additionalNameOutput.CompareTo("0010"))){
		nameGraphPos = Form("%s_SystErrorRelPos_%s%s",meson.Data(),nameCutVariationSC[i].Data(),"1020");
		nameGraphNeg = Form("%s_SystErrorRelNeg_%s%s",meson.Data(),nameCutVariationSC[i].Data(),"1020");   
	      } else {
		nameGraphPos = Form("%s_SystErrorRelPos_%s%s",meson.Data(),nameCutVariationSC[i].Data(),additionalNameOutput.Data()  );
		nameGraphNeg = Form("%s_SystErrorRelNeg_%s%s",meson.Data(),nameCutVariationSC[i].Data(),additionalNameOutput.Data()  );
	      }
	    } else {	
	      nameGraphPos = Form("%s_SystErrorRelPos_%s",meson.Data(),nameCutVariationSC[i].Data());
	      nameGraphNeg = Form("%s_SystErrorRelNeg_%s",meson.Data(),nameCutVariationSC[i].Data());
	    }
	    cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
	    graphPosErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
	    graphNegErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
	  }
	  for (Int_t j = 0; j < offSetBeginning; j++){
	    graphPosErrors->RemovePoint(0);
	    graphNegErrors->RemovePoint(0);
	  }
	  if (i == 0){
	    ptBins = graphNegErrors->GetX();
	    ptBinsErr = graphNegErrors->GetEXhigh();
	  }
	  errorsNeg[i] = graphNegErrors->GetY();
	  errorsNegErr[i] = graphNegErrors->GetEYhigh();
	  errorsPos[i] = graphPosErrors->GetY();
	  errorsPosErr[i] = graphPosErrors->GetEYhigh();
	  
	  cout << nameCutVariationSC[i].Data() << endl;
	  Float_t	constant ;
	  Float_t	constantErr;
	  if (energy.CompareTo("PbPb_2.76TeV") == 0 && nameCutVariationSC[i].CompareTo("V0Offline") == 0) {
	    // 			cout << Form("CutStudies/PbPb_2.76TeV/%s_RatioRawYields.root",meson.Data()) << endl;
	    TFile *InputFileComparisonOffline = new TFile(Form("CutStudies/PbPb_2.76TeV/%s_RatioRawYields.root",meson.Data()));
	    // 			cout <<  Form("%s_fitRatioConstCorrectedYieldCutV0Offline%s_data", meson.Data(), additionalNameOutput.Data()) << endl;
	    TF1 *	fitRatioConstCorrectedYieldCut =(TF1*)InputFileComparisonOffline->Get( Form("%s_fitRatioConstCorrectedYieldCutV0Offline%s_data", meson.Data(), additionalNameOutput.Data()));
	    cout << WriteParameterToFile(fitRatioConstCorrectedYieldCut)<< endl;
	    constant = (1-fitRatioConstCorrectedYieldCut->GetParameter(0))/2*100;
	    constantErr = fitRatioConstCorrectedYieldCut->GetParError(0)/2*100;
	    for (Int_t l = 0; l < nPtBins; l++){
	      errorsNeg[i][l] = constant;
	      errorsMean[i][l] = constant;
	      errorsPos[i][l] = constant;
	      errorsNegCorr[i][l] =constant;
	      errorsMeanCorr[i][l] = constant;
	      errorsPosCorr[i][l] =constant;
	      errorsNegErr[i][l] = constantErr;
	      errorsMeanErr[i][l] = constantErr;
	      errorsPosErr[i][l] = constantErr;
	      errorsNegErrCorr[i][l] =constantErr;
	      errorsMeanErrCorr[i][l] = constantErr;
	      errorsPosErrCorr[i][l] =constantErr;	
	    }
	  } else if (energy.CompareTo("PbPb_2.76TeV") == 0 && nameCutVariationSC[i].CompareTo("PsiPair") == 0) {   
	    if (additionalNameOutput.CompareTo("0005") == 0 || additionalNameOutput.CompareTo("0510") == 0 || additionalNameOutput.CompareTo("0010") == 0){
	      constant = 4;
	      constantErr = 0.04;
	    } else if (additionalNameOutput.CompareTo("1020") == 0 || additionalNameOutput.CompareTo("2040") == 0 ){
	      constant = 3;
	      constantErr = 0.03;
	    } else if (additionalNameOutput.CompareTo("4060") == 0 || additionalNameOutput.CompareTo("6080") == 0 ){
	      constant = 2.5;
	      constantErr = 0.025;
	    }   
	    for (Int_t l = 0; l < nPtBins; l++){
	      errorsNeg[i][l] = constant;
	      errorsMean[i][l] = constant;
	      errorsPos[i][l] = constant;
	      errorsNegCorr[i][l] =constant;
	      errorsMeanCorr[i][l] = constant;
	      errorsPosCorr[i][l] =constant;
	      errorsNegErr[i][l] = constantErr;
	      errorsMeanErr[i][l] = constantErr;
	      errorsPosErr[i][l] = constantErr;
	      errorsNegErrCorr[i][l] =constantErr;
	      errorsMeanErrCorr[i][l] = constantErr;
	      errorsPosErrCorr[i][l] =constantErr;   
	    }	
	  } else if (energy.CompareTo("PbPb_2.76TeV") == 0 && nameCutVariationSC[i].CompareTo("BGEstimate") == 0) {	
	    if (additionalNameOutput.CompareTo("0005") == 0 || additionalNameOutput.CompareTo("0510") == 0 || additionalNameOutput.CompareTo("0010") == 0 || additionalNameOutput.CompareTo("1020") == 0 || additionalNameOutput.CompareTo("2040") == 0){
	      constant = 0;
	      constantErr = 0.0001;
	    } else if ( additionalNameOutput.CompareTo("4060") == 0 ){
	      constant = 0.5;
	      constantErr = 0.005;
	    } else if (additionalNameOutput.CompareTo("6080") == 0 ){
	      constant = 2.;
	      constantErr = 0.02;
	    }   
	    for (Int_t l = 0; l < nPtBins; l++){
	      errorsNeg[i][l] = constant;
	      errorsMean[i][l] = constant;
	      errorsPos[i][l] = constant;
	      errorsNegCorr[i][l] =constant;
	      errorsMeanCorr[i][l] = constant;
	      errorsPosCorr[i][l] =constant;
	      errorsNegErr[i][l] = constantErr;
	      errorsMeanErr[i][l] = constantErr;
	      errorsPosErr[i][l] = constantErr;
	      errorsNegErrCorr[i][l] =constantErr;
	      errorsMeanErrCorr[i][l] = constantErr;
	      errorsPosErrCorr[i][l] =constantErr;   
         }
	  } else if (energy.CompareTo("PbPb_2.76TeV") == 0 && nameCutVariationSC[i].CompareTo("Eta") == 0) {   
	    if (additionalNameOutput.CompareTo("0005") == 0 || additionalNameOutput.CompareTo("0510") == 0 || additionalNameOutput.CompareTo("0010") == 0 || additionalNameOutput.CompareTo("1020") == 0 ){
	      constant = 2;
	      constantErr = 0.02;
	    } else if (additionalNameOutput.CompareTo("2040") == 0 || additionalNameOutput.CompareTo("4060") == 0 || additionalNameOutput.CompareTo("6080") == 0 ){
	      constant = 1.8;
	      constantErr = 0.018;
	    }   
	    for (Int_t l = 0; l < nPtBins; l++){
	      errorsNeg[i][l] = constant;
	      errorsMean[i][l] = constant;
	      errorsPos[i][l] = constant;
	      errorsNegCorr[i][l] =constant;
	      errorsMeanCorr[i][l] = constant;
	      errorsPosCorr[i][l] =constant;
	      errorsNegErr[i][l] = constantErr;
	      errorsMeanErr[i][l] = constantErr;
	      errorsPosErr[i][l] = constantErr;
	      errorsNegErrCorr[i][l] =constantErr;
	      errorsMeanErrCorr[i][l] = constantErr;
	      errorsPosErrCorr[i][l] =constantErr;   
	    }   
	  } else if (energy.CompareTo("PbPb_2.76TeV") == 0 && nameCutVariationSC[i].CompareTo("Alpha") == 0) {   
	    if (additionalNameOutput.CompareTo("0005") == 0 || additionalNameOutput.CompareTo("0510") == 0 || additionalNameOutput.CompareTo("0010") == 0 ){
	      constant = 2;
	      constantErr = 0.02;
	    } else if (additionalNameOutput.CompareTo("1020") == 0 ){
	      constant = 1.5;
	      constantErr = 0.015;
	    } else if (additionalNameOutput.CompareTo("2040") == 0 || additionalNameOutput.CompareTo("4060") == 0 || additionalNameOutput.CompareTo("6080") == 0 ){
	      constant = 1.;
	      constantErr = 0.01;
	    }   
	    for (Int_t l = 0; l < nPtBins; l++){
	      errorsNeg[i][l] = constant;
	      errorsMean[i][l] = constant;
	      errorsPos[i][l] = constant;
	      errorsNegCorr[i][l] =constant;
	      errorsMeanCorr[i][l] = constant;
	      errorsPosCorr[i][l] =constant;
	      errorsNegErr[i][l] = constantErr;
	      errorsMeanErr[i][l] = constantErr;
	      errorsPosErr[i][l] = constantErr;
	      errorsNegErrCorr[i][l] =constantErr;
	      errorsMeanErrCorr[i][l] = constantErr;
	      errorsPosErrCorr[i][l] =constantErr;   
	    }      
	  } else {
	    cout<< " Calculating Syst. Errors with Mean::" << i<<endl;
	    CalculateMeanSysErr(errorsMean[i], errorsMeanErr[i], errorsPos[i], errorsNeg[i], nPtBins);	
	    CorrectSystematicErrorsWithMean(errorsPos[i],errorsPosErr[i], errorsPosCorr[i], errorsPosErrCorr[i], nPtBins);
	    CorrectSystematicErrorsWithMean(errorsNeg[i],errorsNegErr[i], errorsNegCorr[i], errorsNegErrCorr[i], nPtBins);
	    CorrectSystematicErrorsWithMean(errorsMean[i], errorsMeanErr[i], errorsMeanCorr[i], errorsMeanErrCorr[i], nPtBins);
	  }
	  for (Int_t l = 0; l < nPtBins; l++){
	    errorsPosSummed[l] = errorsPosSummed[l]+pow(errorsPos[i][l],2);
	    errorsNegSummed[l] = errorsNegSummed[l]+ pow(errorsNeg[i][l],2);
	    errorsMeanSummed[l] = errorsMeanSummed[l]+ pow(errorsMean[i][l],2);
	    errorsPosCorrSummed[l] = errorsPosCorrSummed[l]+pow(errorsPosCorr[i][l],2);
	    errorsNegCorrSummed[l] = errorsNegCorrSummed[l] +pow(errorsNegCorr[i][l],2);
	    errorsMeanCorrSummed[l] =errorsMeanCorrSummed[l]+ pow(errorsMeanCorr[i][l],2);
	  }
	  negativeErrors[i] = new TGraphErrors(nPtBins,ptBins ,errorsNeg[i] ,ptBinsErr ,errorsNegErr[i] );
	  meanErrors[i] = new TGraphErrors(nPtBins,ptBins ,errorsMean[i] ,ptBinsErr ,errorsMeanErr[i] );
	  positiveErrors[i] = new TGraphErrors(nPtBins,ptBins ,errorsPos[i] ,ptBinsErr ,errorsPosErr[i] );
	  negativeErrorsCorr[i] = new TGraphErrors(nPtBins,ptBins ,errorsNegCorr[i] ,ptBinsErr ,errorsNegErrCorr[i] );
	  meanErrorsCorr[i] = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorr[i] ,ptBinsErr ,errorsMeanErrCorr[i] );
	  positiveErrorsCorr[i] = new TGraphErrors(nPtBins,ptBins ,errorsPosCorr[i] ,ptBinsErr ,errorsPosErrCorr[i] );
	  
	}
	
	Double_t errorMaterial;
	if (energy.CompareTo("PbPb_2.76TeV") == 0) {
	  errorMaterial = 4.86;
	} else {
	  errorMaterial = 4.50;
	}
	
	for (Int_t l = 0; l < nPtBins; l++){
	  errorsPosSummed[l] = pow(errorsPosSummed[l],0.5);
	  errorsMeanSummed[l] = pow(errorsMeanSummed[l],0.5);
	  errorsPosErrSummed[l] = errorsPosSummed[l]*0.001;
	  errorsMeanErrSummed[l] = errorsMeanSummed[l]*0.001;
	  errorsNegSummed[l] = -pow(errorsNegSummed[l],0.5);
	  errorsNegErrSummed[l] = errorsNegSummed[l]*0.001;
	  errorsPosCorrMatSummed[l] = pow(errorsPosCorrSummed[l]+ pow(2*errorMaterial ,2.),0.5);
	  errorsMeanCorrMatSummed[l] = pow(errorsMeanCorrSummed[l]+ pow(2*errorMaterial ,2.),0.5);
	  errorsNegCorrMatSummed[l] = -pow(errorsNegCorrSummed[l]+ pow(2*errorMaterial ,2.),0.5);
	  
	  errorsPosCorrSummed[l] = pow(errorsPosCorrSummed[l],0.5);
	  errorsMeanCorrSummed[l] = pow(errorsMeanCorrSummed[l],0.5);
	  errorsPosErrCorrSummed[l] = errorsPosCorrSummed[l]*0.001;
	  errorsMeanErrCorrSummed[l] = errorsMeanCorrSummed[l]*0.001;
	  errorsMeanErrCorrMatSummed[l] = errorsMeanCorrMatSummed[l]*0.001;
	  errorsNegCorrSummed[l] = -pow(errorsNegCorrSummed[l],0.5);
	  errorsNegErrCorrSummed[l] = errorsNegCorrSummed[l]*0.001;
	  // 		cout << errorsMeanSummed[l] << "\t" << errorsMeanCorrMatSummed[l]<< endl;
	}
	
	Double_t errorsMat[nPtBins];
	for (Int_t l = 0; l < nPtBins; l++){
	  errorsMat[l] = 2*errorMaterial;
	  
	}
	TGraphErrors* graphMaterialError = new TGraphErrors(nPtBins,ptBins ,errorsMat ,ptBinsErr ,errorsMeanErrSummed );
	
	cout << "here" << endl;
	negativeErrorsSummed = new TGraphErrors(nPtBins,ptBins ,errorsNegSummed ,ptBinsErr ,errorsNegErrSummed );
	negativeErrorsCorrSummed = new TGraphErrors(nPtBins,ptBins ,errorsNegCorrSummed ,ptBinsErr ,errorsNegErrCorrSummed );
	positiveErrorsSummed = new TGraphErrors(nPtBins,ptBins ,errorsPosSummed ,ptBinsErr ,errorsPosErrSummed );
	positiveErrorsCorrSummed = new TGraphErrors(nPtBins,ptBins ,errorsPosCorrSummed ,ptBinsErr ,errorsPosErrCorrSummed );
	meanErrorsSummed = new TGraphErrors(nPtBins,ptBins ,errorsMeanSummed ,ptBinsErr ,errorsMeanErrSummed );
	meanErrorsCorrSummed = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSummed ,ptBinsErr ,errorsMeanErrCorrSummed );
	meanErrorsCorrSummedIncMat = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrMatSummed ,ptBinsErr ,errorsMeanErrCorrMatSummed );
	
	
	cout << "here" << endl;
	
	TCanvas* canvasSysErrMean = new TCanvas("canvasSysErrMean","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasSysErrMean, 0.08, 0.01, 0.015, 0.09);
	TH2D *histo2DSysErrMean ;
	
	
	if (energy.CompareTo("PbPb_2.76TeV")==0 && additionalName.CompareTo("0-20%")==0 ) {
	  histo2DSysErrMean = new TH2D("histo2DSysErrMean", "histo2DSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,0.,65.);
	} else	if (energy.CompareTo("PbPb_2.76TeV")==0 && additionalName.CompareTo("40-60%")==0 ) {
	  histo2DSysErrMean = new TH2D("histo2DSysErrMean", "histo2DSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,0.,80.);
	} else if (energy.CompareTo("PbPb_2.76TeV")==0 && additionalName.CompareTo("20-40%")==0 ) {
	  histo2DSysErrMean = new TH2D("histo2DSysErrMean", "histo2DSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,0.,60.);		
	} else if (energy.CompareTo("PbPb_2.76TeV")==0 && additionalName.CompareTo("60-80%")==0 ) {
	  histo2DSysErrMean = new TH2D("histo2DSysErrMean", "histo2DSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,0.,70.);		
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
	
	TLegend* legendMean = new TLegend(0.13,0.6,0.5,0.95);
	legendMean->SetTextSize(0.035);
	legendMean->SetFillColor(0);
	legendMean->SetBorderSize(0);
	if (numberCutStudies> 9) legendMean->SetNColumns(2);
	for(Int_t i = 0; i< numberCutStudies ; i++){
	  DrawGammaSetMarkerTGraphErr(meanErrors[i], 20+i, 1.,color[i],color[i]);
	  meanErrors[i]->Draw("p,csame");
	  legendMean->AddEntry(meanErrors[i],nameCutVariation[i].Data(),"p");
	}
	DrawGammaSetMarkerTGraphErr(meanErrorsSummed, 20, 1.,1,1);
	meanErrorsSummed->Draw("p,csame");
	legendMean->AddEntry(meanErrorsSummed,"quad. sum.","p");
	legendMean->Draw();
	canvasSysErrMean->Update();
	canvasSysErrMean->SaveAs(Form("SystematicErrorsNew/SysMean_%s_%s%s_%s.eps",meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data()));
	
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
	      canvasSysErrMeanSingle->SaveAs(Form("SystematicErrorsNew/SysMean_%sCutVariation_%s_%s%s_%s.eps",nameCutForOutputVariation7TeV[i].Data(),meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data()));
	    }	
	  }
	  delete canvasSysErrMean;
	}
	
	TCanvas* canvasNewSysErrMean = new TCanvas("canvasNewSysErrMean","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasNewSysErrMean, 0.08, 0.01, 0.015, 0.09);
	TH2D *histo2DNewSysErrMean ;
	
	
	if (energy.CompareTo("PbPb_2.76TeV")==0 && additionalName.CompareTo("0-20%")==0 ) {
	  histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "histo2DNewSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,60.);
	} else	if (energy.CompareTo("PbPb_2.76TeV")==0 && additionalName.CompareTo("40-60%")==0 ) {
	  histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "histo2DNewSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,55.);
	} else	if (energy.CompareTo("PbPb_2.76TeV")==0 && additionalName.CompareTo("20-40%")==0 ) {
	  histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "histo2DNewSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,50.);
	} else	if (energy.CompareTo("PbPb_2.76TeV")==0 && additionalName.CompareTo("60-80%")==0 ) {
	  histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "histo2DNewSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,50.);
	} else 	if (energy.CompareTo("PbPb_2.76TeV")==0) {
	  histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "histo2DNewSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,50.);
	} else {
	  histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "histo2DNewSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,25.);
	}
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
	
	TLegend* legendMeanNew;
	if (energy.CompareTo("PbPb_2.76TeV")==0 && additionalName.CompareTo("60-80%")==0 ){
	  legendMeanNew= new TLegend(0.13,0.6,0.55,0.95);
	} else {
	  legendMeanNew= new TLegend(0.2,0.6,0.62,0.95);
	}
	
	legendMeanNew->SetTextSize(0.035);
	legendMeanNew->SetFillColor(0);
	legendMeanNew->SetBorderSize(0);
	if (numberCutStudies> 9) legendMeanNew->SetNColumns(2);
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
	legendMeanNew->AddEntry(meanErrorsCorrSummed,"quad. sum.","p");
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kRed,kRed);
	meanErrorsCorrSummedIncMat->Draw("p,csame");
	legendMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"quad. sum., inc. mat.","p");
	legendMeanNew->Draw();
	
	if (energy.CompareTo("PbPb_2.76TeV")==0){
	  TLatex *labelCentrality = new TLatex(0.65,0.93,Form("%s Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV"	,additionalName.Data()	));
	  SetStyleTLatex( labelCentrality, 0.038,4);
	  labelCentrality->Draw();
	} else {
	  if (meson.Contains("Pi0")){
	    TLatex *labelCollisionSystem = new TLatex(0.75,0.93,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}, pp, #sqrt{#it{s}} = 2.76 TeV");
	    SetStyleTLatex( labelCollisionSystem, 0.038,4);
	    labelCollisionSystem->Draw();
	  } else {
	    TLatex *labelCollisionSystem = new TLatex(0.75,0.93,"#eta #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}, pp, #sqrt{#it{s}} = 2.76 TeV");
	    SetStyleTLatex( labelCollisionSystem, 0.038,4);
	    labelCollisionSystem->Draw();
	  }
	}
	
	canvasNewSysErrMean->Update();
	canvasNewSysErrMean->SaveAs(Form("SystematicErrorsNew/SysMeanNewWithMean_%s_%s%s_%s.eps",meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data()));
	
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
	      canvasSysErrMeanSingle->SaveAs(Form("SystematicErrorsNew/SysMeanSmoothed_%sCutVariation_%s_%s%s_%s.eps", nameCutForOutputVariation7TeV[i].Data(), meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data()));
	      
	    }	
	  }
	  delete canvasSysErrMean;
	}
	
	const char *SysErrDatname = Form("SystematicErrorsNew/SystematicError_%s_%s%s_%s.dat",meson.Data(),energy.Data(),additionalNameOutput.Data(),dateForOutput.Data());
	fstream SysErrDat;
	cout << SysErrDatname << endl;
	SysErrDat.open(SysErrDatname, ios::out);
	for (Int_t l=0; l< nPtBins; l++){
	  SysErrDat <<errorsNegCorrMatSummed[l] << "\t" <<errorsPosCorrMatSummed[l] << "\t"  <<errorsNegCorrSummed[l] << "\t" <<errorsPosCorrSummed[l]  << endl;
	}
	
	SysErrDat.close();
	
	const char *SysErrDatnameMean = Form("SystematicErrorsNew/SystematicErrorAveraged_%s_%s%s_%s.dat",meson.Data(),energy.Data(),additionalNameOutput.Data(),dateForOutput.Data());
	fstream SysErrDatAver;
	cout << SysErrDatnameMean << endl;
	SysErrDatAver.open(SysErrDatnameMean, ios::out);
	for (Int_t l=0; l< nPtBins; l++){
	  SysErrDatAver << "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
	}
	
	SysErrDatAver.close();
	
	Double_t errorsMeanCorrPID[nPtBins];	
	Double_t errorsMeanCorrSignalExtraction[nPtBins];	
	Double_t errorsMeanCorrTrackReco[nPtBins];	
	Double_t errorsMeanCorrPhotonReco[nPtBins];	
	
	for (Int_t l=0; l< nPtBins; l++){
	  if (energy.CompareTo("PbPb_2.76TeV")==0){
	    // 				"YieldExtraction-0","dEdxE-1","dEdxPi-2", "Cluster-3", "TOF-4", "SinglePt-5", "Chi2-6", "Qt-7", "Alpha-8", "PsiPair-9", "Eta-10", "BGEstimate-11","V0Offline-12", "R-13", "BG-14",
	    if (numberCutStudies>8){
	      errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]+ errorsMeanCorr[1][l]*errorsMeanCorr[1][l]+errorsMeanCorr[8][l]*errorsMeanCorr[8][l]);
	    } else if (numberCutStudies>13){   
	      errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]+ errorsMeanCorr[13][l]*errorsMeanCorr[13][l]+errorsMeanCorr[8][l]*errorsMeanCorr[8][l]);
	    } else{
	      errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]);
	    }
	    errorsMeanCorrPID[l] =TMath::Sqrt(errorsMeanCorr[1][l]*errorsMeanCorr[1][l]+ errorsMeanCorr[2][l]*errorsMeanCorr[2][l]+ errorsMeanCorr[4][l]*errorsMeanCorr[4][l]);
	    if (numberCutStudies>10){
	      errorsMeanCorrPhotonReco[l] =TMath::Sqrt( errorsMeanCorr[6][l]* errorsMeanCorr[6][l]+errorsMeanCorr[7][l]*errorsMeanCorr[7][l]+ errorsMeanCorr[9][l]*errorsMeanCorr[9][l]+ errorsMeanCorr[10][l]*errorsMeanCorr[10][l]);
	    } else if (numberCutStudies>9){
	      errorsMeanCorrPhotonReco[l] =TMath::Sqrt( errorsMeanCorr[6][l]* errorsMeanCorr[6][l]+errorsMeanCorr[7][l]*errorsMeanCorr[7][l]+ errorsMeanCorr[9][l]*errorsMeanCorr[9][l]);
	    }   
	    else {
	      errorsMeanCorrPhotonReco[l] =TMath::Sqrt( errorsMeanCorr[6][l]* errorsMeanCorr[6][l]+errorsMeanCorr[7][l]*errorsMeanCorr[7][l]);
	    }
	    errorsMeanCorrTrackReco[l] = TMath::Sqrt(errorsMeanCorr[3][l]*errorsMeanCorr[3][l]+errorsMeanCorr[5][l]*errorsMeanCorr[5][l]);
	  } else if (energy.CompareTo("2.76TeV")==0 && meson.CompareTo("Eta") == 0  ) {
	    // {"YieldExtraction-0","BG-1", "TPCCluster-2", "SinglePt-3", "Chi2-4", "Qt-5", "Alpha-6"};
	    errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]+errorsMeanCorr[1][l]*errorsMeanCorr[1][l]);
	    errorsMeanCorrTrackReco[l] = TMath::Sqrt(errorsMeanCorr[2][l]*errorsMeanCorr[2][l]+errorsMeanCorr[3][l]*errorsMeanCorr[3][l]);
	    errorsMeanCorrPhotonReco[l] = TMath::Sqrt( errorsMeanCorr[4][l]* errorsMeanCorr[4][l]+errorsMeanCorr[5][l]*errorsMeanCorr[5][l]);
	    errorsMeanCorrPID[l] = 0;
	  } else {
	    // 				"YieldExtraction-0","BG-1", "BGN-2", "dEdxE-3","dEdxPi-4", "Cluster-5", "SinglePt-6", "Chi2-7","Qt-8", "Alpha-9", "TrackShare-10",  "V0Offline-11","R-12" ,"PsiPair-13", "Eta-14"
	    if (numberCutStudies>9){
	      errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]+ errorsMeanCorr[1][l]*errorsMeanCorr[1][l]+ errorsMeanCorr[2][l]*errorsMeanCorr[2][l]+errorsMeanCorr[9][l]*errorsMeanCorr[9][l]);
	    } else {
	      errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]+ errorsMeanCorr[1][l]*errorsMeanCorr[1][l]+ errorsMeanCorr[2][l]*errorsMeanCorr[2][l]);
	    }
	    errorsMeanCorrPID[l] =TMath::Sqrt(errorsMeanCorr[3][l]*errorsMeanCorr[3][l]+ errorsMeanCorr[4][l]*errorsMeanCorr[4][l]);
	    errorsMeanCorrPhotonReco[l] = errorsMeanCorr[7][l] ;
	    if (numberCutStudies>12){
	      errorsMeanCorrPhotonReco[l] =TMath::Sqrt( errorsMeanCorr[7][l]* errorsMeanCorr[7][l]+errorsMeanCorr[8][l]*errorsMeanCorr[8][l]+errorsMeanCorr[11][l]*errorsMeanCorr[11][l]+errorsMeanCorr[12][l]*errorsMeanCorr[12][l]);
	    } else if (numberCutStudies>11){
	      errorsMeanCorrPhotonReco[l] =TMath::Sqrt( errorsMeanCorr[7][l]* errorsMeanCorr[7][l]+errorsMeanCorr[8][l]*errorsMeanCorr[8][l] +errorsMeanCorr[11][l]*errorsMeanCorr[11][l]);
	    } else {
	      errorsMeanCorrPhotonReco[l] =TMath::Sqrt( errorsMeanCorr[7][l]* errorsMeanCorr[7][l]+errorsMeanCorr[8][l]*errorsMeanCorr[8][l]);
	    }
	    errorsMeanCorrTrackReco[l] = TMath::Sqrt(errorsMeanCorr[5][l]*errorsMeanCorr[5][l]+errorsMeanCorr[6][l]*errorsMeanCorr[6][l]);
	    
	  }
	}
	TGraphErrors* meanErrorsPID = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrPID ,ptBinsErr ,errorsMeanErrCorrSummed );
	TGraphErrors* meanErrorsPhotonReco = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrPhotonReco ,ptBinsErr ,errorsMeanErrCorrSummed );
	TGraphErrors* meanErrorsSignalExtraction = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSignalExtraction ,ptBinsErr ,errorsMeanErrCorrSummed );
	TGraphErrors* meanErrorsTrackReco = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrTrackReco ,ptBinsErr ,errorsMeanErrCorrSummed );
	
	cout << "here" << endl;
   
   
	TCanvas* canvasSummedErrMean = new TCanvas("canvasSummedErrMean","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasSummedErrMean, 0.08, 0.01, 0.015, 0.09);
	TH2D *histo2DSummedErrMean ;
	
	
	if (energy.CompareTo("PbPb_2.76TeV")==0 && meson.CompareTo("Eta") == 0 ) {
	  histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "histo2DSummedErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,100.);
	} else	if (energy.CompareTo("PbPb_2.76TeV")==0 ) {
	  histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "histo2DSummedErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,50.);
	} else 	if (energy.CompareTo("2.76TeV")==0 && meson.CompareTo("Eta") == 0  ){
	  histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "histo2DSummedErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,40.);
	} else {
	  histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "histo2DSummedErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,25.);
	}
	histo2DSummedErrMean->SetYTitle("mean smoothed systematic Err %");
	histo2DSummedErrMean->SetXTitle("p_{T} (GeV/c)");
	histo2DSummedErrMean->GetYaxis()->SetNdivisions(510,kTRUE);
	histo2DSummedErrMean->GetYaxis()->SetLabelSize(0.03);
	histo2DSummedErrMean->GetYaxis()->SetTitleSize(0.04);
	histo2DSummedErrMean->GetYaxis()->SetDecimals();
	histo2DSummedErrMean->GetYaxis()->SetTitleOffset(0.9);
	histo2DSummedErrMean->GetXaxis()->SetTitleOffset(1.);
	histo2DSummedErrMean->GetXaxis()->SetLabelSize(0.03);
	histo2DSummedErrMean->GetXaxis()->SetTitleSize(0.04);
	histo2DSummedErrMean->SetTitle("");
	histo2DSummedErrMean->Draw();
	
	TLegend* legendSummedMeanNew;
	if (energy.CompareTo("PbPb_2.76TeV")==0 && additionalName.CompareTo("60-80%")==0 ){
	  legendSummedMeanNew= new TLegend(0.13,0.7,0.55,0.95);
	} else {
	  legendSummedMeanNew= new TLegend(0.15,0.7,0.45,0.95);
	}
	legendSummedMeanNew->SetTextSize(0.035);
	legendSummedMeanNew->SetFillColor(0);
	legendSummedMeanNew->SetBorderSize(0);
	
	
	DrawGammaSetMarkerTGraphErr(meanErrorsSignalExtraction, 20, 1.,color[0],color[0]);
	meanErrorsSignalExtraction->Draw("p,csame");
   cout << "here" << endl;
	DrawGammaSetMarkerTGraphErr(meanErrorsPID, 21, 1.,color[1],color[1]);
	meanErrorsPID->Draw("p,csame");
   cout << "here" << endl;
	DrawGammaSetMarkerTGraphErr(meanErrorsTrackReco, 22, 1.,color[2],color[2]);
	meanErrorsTrackReco->Draw("p,csame");
   cout << "here" << endl;
	DrawGammaSetMarkerTGraphErr(meanErrorsPhotonReco, 23, 1.,color[3],color[3]);
	meanErrorsPhotonReco->Draw("p,csame");
   cout << "here" << endl;
   if (energy.CompareTo("2.76TeV")==0 ){ 
      DrawGammaSetMarkerTGraphErr(meanErrorsCorr[10], 25, 1.,color[5],color[5]);
      meanErrorsCorr[10]->Draw("p,csame");
      cout << "here" << endl;
   }
   if (energy.CompareTo("PbPb_2.76TeV")==0 && (additionalName.CompareTo("60-80%")==0 || additionalName.CompareTo("40-60%")==0)){ 
      DrawGammaSetMarkerTGraphErr(meanErrorsCorr[11], 25, 1.,color[5],color[5]);
      meanErrorsCorr[11]->Draw("p,csame");
      cout << "here" << endl;
   }
	DrawGammaSetMarkerTGraphErr(graphMaterialError, 24, 1.,color[4],color[4]);
	graphMaterialError->Draw("p,csame");
   cout << "here" << endl;
	
   legendSummedMeanNew->AddEntry(meanErrorsSignalExtraction,"Signal Extraction","p");
   legendSummedMeanNew->AddEntry(meanErrorsPID,"Electron PID","p");
   legendSummedMeanNew->AddEntry(meanErrorsTrackReco,"Track Reconstruction","p");
   legendSummedMeanNew->AddEntry(meanErrorsPhotonReco,"Photon Reconstruction","p");
   if (energy.CompareTo("2.76TeV")==0 )legendSummedMeanNew->AddEntry(meanErrorsCorr[10],"Pileup Estimate","p");
   if (energy.CompareTo("PbPb_2.76TeV")==0 && (additionalName.CompareTo("60-80%")==0 || additionalName.CompareTo("40-60%")==0)) legendSummedMeanNew->AddEntry(meanErrorsCorr[11],"Pileup Estimate","p");
   legendSummedMeanNew->AddEntry(graphMaterialError,"Material","p");
	
// 	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummed, 20, 1.,1,1);
// 	meanErrorsCorrSummed->Draw("p,csame");
// 	legendSummedMeanNew->AddEntry(meanErrorsCorrSummed,"quad. sum.","p");
	
   DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kBlack,kBlack);
   meanErrorsCorrSummedIncMat->Draw("p,csame");
   legendSummedMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"quad. sum.","p");
   legendSummedMeanNew->Draw();
   
   if (energy.CompareTo("PbPb_2.76TeV")==0){
     TLatex *labelCentrality = new TLatex(0.65,0.93,Form("%s Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV"	,additionalName.Data()	));
     SetStyleTLatex( labelCentrality, 0.038,4);
     labelCentrality->Draw();
   } else {
     if (meson.Contains("Pi0")){
       TLatex *labelCollisionSystem = new TLatex(0.55,0.93,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}, pp, #sqrt{#it{s}} = 2.76 TeV");
       SetStyleTLatex( labelCollisionSystem, 0.038,4);
       labelCollisionSystem->Draw();
     } else {
       TLatex *labelCollisionSystem = new TLatex(0.55,0.93,"#eta #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}, pp, #sqrt{#it{s}} = 2.76 TeV");
       SetStyleTLatex( labelCollisionSystem, 0.038,4);
       labelCollisionSystem->Draw();
     }
   }
   
   canvasSummedErrMean->Update();
   canvasSummedErrMean->SaveAs(Form("SystematicErrorsNew/SysErrorSummedVisu_%s_%s%s_%s.eps",meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data()));
   
   delete canvasSummedErrMean;
   
   
   
   const char *SysErrDatnameMeanPaper = Form("SystematicErrorsNew/SystematicErrorAveragedPaper_%s_%s%s_%s.dat",meson.Data(),energy.Data(),additionalNameOutput.Data(),dateForOutput.Data());
   fstream SysErrDatAverPaper;
   cout << SysErrDatnameMeanPaper << endl;
   SysErrDatAverPaper.open(SysErrDatnameMeanPaper, ios::out);
   SysErrDatAverPaper  << "p_{T}" << "\t Material \t Yield Extraction \t PID \t photon reco \t track recon \t summed" <<  endl;
   for (Int_t l=0; l< nPtBins; l++){
     SysErrDatAverPaper << ptBins[l] <<"\t" << errorMaterial*2 << "\t" << errorsMeanCorrSignalExtraction[l] << "\t" << errorsMeanCorrPID[l]<< "\t" << errorsMeanCorrPhotonReco[l]<< "\t" <<errorsMeanCorrTrackReco[l] <<"\t" << errorsMeanCorrMatSummed[l]<< endl;
   }
   
   SysErrDatAverPaper.close();
}
