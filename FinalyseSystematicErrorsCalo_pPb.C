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

void FinalyseSystematicErrorsCalo_pPb(	const char* nameDataFileErrors ="", 
											TString energy="", 
											TString meson = "", 
											Int_t numberOfPtBins =1, 
											Int_t numberCutStudies=1, 
											Int_t offSetBeginning = 0, 
											TString additionalName = "pp", 
											TString additionalNameOutput = "", 
											TString suffix = "eps",
											Int_t mode = 4
										){
	
	StyleSettingsThesis();	
	SetPlotStyle();
	
	TString date = ReturnDateString();
	TString dateForOutput = ReturnDateStringForOutput();
	TString collisionSystem= ReturnFullCollisionsSystem(energy);
	
	Color_t color[20] = {860,kCyan-2,807,880,418,kBlue+2,kViolet+2,kRed+2,kTeal+3,kGray+2,kOrange+2,435,kBlue-7,kSpring+4,416,830,404,608,kCyan-2,1};
	
	Int_t 	numberOfEntriesPos = 0;
	Int_t 	numberOfEntriesNeg = 0;
	
	TFile* fileErrorInput= new TFile(nameDataFileErrors);
	const Int_t nPtBins = numberOfPtBins;
	const Int_t nCuts = numberCutStudies;
	Double_t* ptBins;
	Double_t* ptBinsErr;
	
	TString nameCutVariation2760GeV[14] = {"Yield extraction", "min E_{cluster}", "min # cells", "clu. energy calibration", 
										  "tr. match. to cl.", "M_{02}", "mat. infront of EMCal","clusterizer","clu. timing", 
										  "feed down correction", "energy scale", "material #pi^{0}", "fit function for yield", 
										  "trigger normalization"};  
	
	TString nameCutVariationSC2760GeV[14] = {"YieldExtraction", "ClusterMinEnergy", "ClusterNCells", "NonLinearity", 
											"TrackMatching", "ClusterM02","ClusterMaterialTRD","Clusterizer","ClusterTiming", 
											"FeedDownCorrection","EnergyScale","MatPi0", "FitFunc", "Trigger"};
											
	TString nameCutVariation5023GeV[14] = {"Yield extraction", "min E_{cluster}", "min # cells", "clu. energy calibration", 
										  "M_{02}", "mat. infront of EMCal","Rapidity", "Alpha", 
										  "energy scale", "material #pi^{0}", "fit function for yield", 
										  "trigger normalization", 
										  "tr. match. to cl.","clusterizer"};  
	
	TString nameCutVariationSC5023GeV[14] = {"YieldExtraction", "MinEClus", "MinNCell", "NonLinCor", 
											"M02","TRD","Rapidity","Alpha",
											"EnergyScale","MatPi0", "FitFunc", "Trigger", 
											"TrackMatching","Clusterizer"};
											
	TString nameCutVariation[14];
	TString nameCutVariationSC[14];
// 	Bool_t bsmooth[13] = {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	Bool_t bsmooth[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	if (energy.CompareTo("2.76TeV") == 0) {
		for (Int_t i = 0; i < numberCutStudies; i++){
			nameCutVariation[i] = nameCutVariation2760GeV[i];
			nameCutVariationSC[i] = nameCutVariationSC2760GeV[i];
		}
	}
	else if (energy.CompareTo("5.023TeV") == 0) {
		for (Int_t i = 0; i < numberCutStudies; i++){
			nameCutVariation[i] = nameCutVariation5023GeV[i];
			nameCutVariationSC[i] = nameCutVariationSC5023GeV[i];
		}
	}
	
	gSystem->Exec("mkdir -p SystematicErrorsCalculatedCalo");
	
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
		if (i == 0 ){
			TString nameGraphPos = Form("%s_SystErrorRelPos_%s_%s",meson.Data(),nameCutVariationSC[i].Data(),additionalName.Data()	);
			TString nameGraphNeg = Form("%s_SystErrorRelNeg_%s_%s",meson.Data(),nameCutVariationSC[i].Data(),additionalName.Data()	);
			cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
			graphPosErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
			graphNegErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
		} /*else if (i == 6 || i == 9 || i == 10 || i == 11|| i == 12 || i==13){	
			TString nameGraphPos = Form("%s_SystErrorRelPos_%s%s",meson.Data(),nameCutVariationSC[5].Data(),additionalNameOutput.Data()  );
			TString nameGraphNeg = Form("%s_SystErrorRelNeg_%s%s",meson.Data(),nameCutVariationSC[5].Data(),additionalNameOutput.Data()  );
			cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
			graphPosErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
			graphNegErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());			} */
		else {
			TString nameGraphPos = Form("%s_SystErrorRelPos_%s%s",meson.Data(),nameCutVariationSC[i].Data(),additionalNameOutput.Data()  );
			TString nameGraphNeg = Form("%s_SystErrorRelNeg_%s%s",meson.Data(),nameCutVariationSC[i].Data(),additionalNameOutput.Data()  );
			cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
			graphPosErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
			graphNegErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
		}
		
		for (Int_t j = 0; j < offSetBeginning; j++){
			graphPosErrors->RemovePoint(0);
			graphNegErrors->RemovePoint(0);
		}
		if (i == 0) {
			ptBins = graphNegErrors->GetX();
			ptBinsErr = graphNegErrors->GetEXhigh();
		}
		errorsNeg[i] = graphNegErrors->GetY();
		errorsNegErr[i] = graphNegErrors->GetEYhigh();
		errorsPos[i] = graphPosErrors->GetY();
		errorsPosErr[i] = graphPosErrors->GetEYhigh();
			
		cout << nameCutVariationSC[i].Data() << endl;
		CalculateMeanSysErr(errorsMean[i], errorsMeanErr[i], errorsPos[i], errorsNeg[i], nPtBins);	
		CorrectSystematicErrorsWithMean(errorsPos[i],errorsPosErr[i], errorsPosCorr[i], errorsPosErrCorr[i], nPtBins);
		CorrectSystematicErrorsWithMean(errorsNeg[i],errorsNegErr[i], errorsNegCorr[i], errorsNegErrCorr[i], nPtBins);
		CorrectSystematicErrorsWithMean(errorsMean[i], errorsMeanErr[i], errorsMeanCorr[i], errorsMeanErrCorr[i], nPtBins);

		// loop for adjusting errors to fit sys err
		for (Int_t k = 0; k < nPtBins; k++){
			errorsMeanErr[i][k] = 0.03;
			errorsMeanErrCorr[i][k] = 0.03;
		}	
		
		if (bsmooth[i]){ // only apply if error needs to be smoothed
			// Cut 1
			if (nameCutVariationSC[i].CompareTo("ClusterMinEnergy")==0  && additionalNameOutput.CompareTo("") == 0){
				for (Int_t k = 0; k < nPtBins; k++){
					Double_t error = 0.264+0.15*ptBins[k];
					errorsMean[i][k] = error;
					errorsMeanErr[i][k] = 0.01*error;
					errorsMeanCorr[i][k] = error;
					errorsMeanErrCorr[i][k] = 0.01*error;
				}
			} else if (nameCutVariationSC[i].CompareTo("ClusterMinEnergy")==0  && additionalNameOutput.CompareTo("EMC1") == 0){
				for (Int_t k = 0; k < nPtBins; k++){
					Double_t error = 2.57-0.54*ptBins[k]+0.034*ptBins[k]*ptBins[k];
					errorsMean[i][k] = error;
					errorsMeanErr[i][k] = 0.01*error;
					errorsMeanCorr[i][k] = error;
					errorsMeanErrCorr[i][k] = 0.01*error;
				}			
			}	

			// Cut 2
			if (nameCutVariationSC[i].CompareTo("ClusterNCells")==0 && additionalNameOutput.CompareTo("") == 0){ //&& meson.Contains("Pi0")
				for (Int_t k = 0; k < nPtBins; k++){
					Double_t error = 0.73+(-0.18)*ptBins[k]+(0.05)*ptBins[k]*ptBins[k];
					errorsMean[i][k] = error;
					errorsMeanErr[i][k] = 0.01*error;
					errorsMeanCorr[i][k] = error;
					errorsMeanErrCorr[i][k] = 0.01*error;
				}	
			} else if (nameCutVariationSC[i].CompareTo("ClusterNCells")==0 && additionalNameOutput.CompareTo("EMC1") == 0){ //&& meson.Contains("Pi0")
				for (Int_t k = 0; k < nPtBins; k++){
					Double_t error = 1.0898294285;
					errorsMean[i][k] = error;
					errorsMeanErr[i][k] = 0.01*error;
					errorsMeanCorr[i][k] = error;
					errorsMeanErrCorr[i][k] = 0.01*error;
				}	
			}	
			
			// Cut 3
			if (nameCutVariationSC[i].CompareTo("NonLinearity")==0 ){ //&& meson.Contains("Pi0")
				Int_t start = 2;
				if (additionalNameOutput.CompareTo("EMC1") == 0) start = 0;
				for (Int_t k = start; k < nPtBins; k++){
					cout << ptBins[k] << endl;
					Double_t error = 3.2785720508;
	// 				Double_t error = 4.07-1.6*ptBins[k]+0.326*ptBins[k]*ptBins[k]-0.0011*ptBins[k]*ptBins[k]*ptBins[k]*ptBins[k]; // parametrisation with No NonLinearity in 
					errorsMean[i][k] = error;
					errorsMeanErr[i][k] = error*0.01;
					errorsMeanCorr[i][k] = error;
					errorsMeanErrCorr[i][k] = error*0.01;
				}
			}

			// Cut 4
			if (nameCutVariationSC[i].CompareTo("TrackMatching")==0 && additionalNameOutput.CompareTo("") == 0){
				for (Int_t k = 1; k < nPtBins; k++){
					Double_t error = 1.62;
					errorsMean[i][k] = error;
					errorsMeanErr[i][k] = 0.01*error;
					errorsMeanCorr[i][k] = error;
					errorsMeanErrCorr[i][k] = 0.01*error;
				}	
			} else if (nameCutVariationSC[i].CompareTo("TrackMatching")==0 && additionalNameOutput.CompareTo("EMC1") == 0){	
				for (Int_t k = 1; k < nPtBins; k++){
					Double_t error = 3.3375247088;
					errorsMean[i][k] = error;
					errorsMeanErr[i][k] = 0.01*error;
					errorsMeanCorr[i][k] = error;
					errorsMeanErrCorr[i][k] = 0.01*error;
				}				
			}
			
			// Cut 5
			if (nameCutVariationSC[i].CompareTo("ClusterM02")==0 ){ //&& meson.Contains("Pi0")
				for (Int_t k = 0; k < nPtBins; k++){
					Double_t error = 1.4; // parametrisation
					errorsMean[i][k] = error;
					errorsMeanErr[i][k] = error*0.01;
					errorsMeanCorr[i][k] = error;
					errorsMeanErrCorr[i][k] = error*0.01;
				}
			} 
			// Cut 6
			if (nameCutVariationSC[i].CompareTo("ClusterMaterialTRD")==0 ){
				for (Int_t k = 0; k < nPtBins; k++){
					errorsMean[i][k] = 7.07; //(5% for TRD mat, 5% for TOF mat added in quadrature)
					errorsMeanErr[i][k] = 0.07;
					errorsMeanCorr[i][k] = 7.07;
					errorsMeanErrCorr[i][k] = 0.07;
				}	
			}

			// Cut 7
			if (nameCutVariationSC[i].CompareTo("Clusterizer")==0 ){
				Int_t start = 2;
				if (additionalNameOutput.CompareTo("EMC1")==0)  start = 2;
				for (Int_t k = 0; k < nPtBins; k++){
					Double_t error = 2.2666741720;
					errorsMean[i][k] = error;
					errorsMeanErr[i][k] = 0.01*error;
					errorsMeanCorr[i][k] = error;
					errorsMeanErrCorr[i][k] = 0.01*error;
				}	
			}	

			// Cut 8 
			if (nameCutVariationSC[i].CompareTo("ClusterTiming")==0 && additionalNameOutput.CompareTo("") == 0){
				for (Int_t k = 0; k < nPtBins; k++){
					Double_t error = 0.57;
					errorsMean[i][k] = error;
					errorsMeanErr[i][k] = 0.01*error;
					errorsMeanCorr[i][k] = error;
					errorsMeanErrCorr[i][k] = 0.01*error;
				}	
			} else if (nameCutVariationSC[i].CompareTo("ClusterTiming")==0 && additionalNameOutput.CompareTo("EMC1") == 0){
				for (Int_t k = 0; k < nPtBins; k++){
					Double_t error = 1.26;
					errorsMean[i][k] = error;
					errorsMeanErr[i][k] = 0.01*error;
					errorsMeanCorr[i][k] = error;
					errorsMeanErrCorr[i][k] = 0.01*error;
				}				
			}	
			
			// Cut 9
			if (nameCutVariationSC[i].CompareTo("FeedDownCorrection")==0 ){ //&& meson.Contains("Pi0")
				for (Int_t k = 0; k < nPtBins; k++){
					Double_t error = 2.99999999999999989e-01*( ( (TMath::Exp((2.03566999505870525e-02/ptBins[k])
																+2.79316247698967624e-01)+(-1.27582401175446813e+00))
																* (1.0-(1.0/(TMath::Exp(-((ptBins[k]-1.18295907157045721e+00)/1.00000000000841915e-01))+1))))
															+((1.0/(TMath::Exp(-((ptBins[k]-1.18295907157045721e+00)/1.00000000000841915e-01))+1))
																*(TMath::Exp((1.07513843613599516e-02/ptBins[k])+1.13403734664314237e+00)+(-3.06439046375136526e+00))))*100;
					errorsMean[i][k] = error;
					errorsMeanErr[i][k] = error*0.01;
					errorsMeanCorr[i][k] = error;
					errorsMeanErrCorr[i][k] = error*0.01;
				}
			}

			// Cut 10
			if (nameCutVariationSC[i].CompareTo("EnergyScale")==0 ){
				for (Int_t k = 0; k < nPtBins; k++){
					errorsMean[i][k] = 2.9; //(5% for TRD mat, 5% for TOF mat added in quadrature)
					errorsMeanErr[i][k] = 0.03;
					errorsMeanCorr[i][k] = 2.9;
					errorsMeanErrCorr[i][k] = 0.03;
				}	
			}
			
			// Cut 11
			if (nameCutVariationSC[i].CompareTo("MatPi0")==0 ){ //&& meson.Contains("Pi0")
				for (Int_t k = 0; k < nPtBins; k++){
					Double_t error = TMath::Exp((-2.39068509299404530e+00)+(-2.21799999999999997e+00)*ptBins[k])*100;
					errorsMean[i][k] = error;
					errorsMeanErr[i][k] = error*0.01;
					errorsMeanCorr[i][k] = error;
					errorsMeanErrCorr[i][k] = error*0.01;
				}
			}
			
			// Cut 12
			if (nameCutVariationSC[i].CompareTo("FitFunc")==0 ){ //&& meson.Contains("Pi0")
				for (Int_t k = 0; k < nPtBins; k++){
					Double_t error = (TMath::Exp((2.50000000000000000e+00)+(-6.78000000000000025e+00)*ptBins[k])+1.49999999999999994e-02)*100;
					errorsMean[i][k] = error;
					errorsMeanErr[i][k] = error*0.01;
					errorsMeanCorr[i][k] = error;
					errorsMeanErrCorr[i][k] = error*0.01;
				}
			}
			
			// Cut 13
			if (nameCutVariationSC[i].CompareTo("Trigger")==0 && additionalNameOutput.CompareTo("EMC1") == 0){
				for (Int_t k = 0; k < nPtBins; k++){
					Double_t error = 4.6;
					errorsMean[i][k] = error;
					errorsMeanErr[i][k] = 0.01*error;
					errorsMeanCorr[i][k] = error;
					errorsMeanErrCorr[i][k] = 0.01*error;
				}				
				
			}	
		}
		
		if (!nameCutVariationSC[i].CompareTo("ClusterMaterialTRD")==0){
			cout << "errors added quadratically" << endl;
			for (Int_t l = 0; l < nPtBins; l++){
				errorsPosSummed[l] = errorsPosSummed[l]+pow(errorsPos[i][l],2);
				errorsNegSummed[l] = errorsNegSummed[l]+ pow(errorsNeg[i][l],2);
				errorsMeanSummed[l] = errorsMeanSummed[l]+ pow(errorsMean[i][l],2);
				errorsPosCorrSummed[l] = errorsPosCorrSummed[l]+pow(errorsPosCorr[i][l],2);
				errorsNegCorrSummed[l] = errorsNegCorrSummed[l] +pow(errorsNegCorr[i][l],2);
				errorsMeanCorrSummed[l] =errorsMeanCorrSummed[l]+ pow(errorsMeanCorr[i][l],2);
			}
		}	
		negativeErrors[i] = new TGraphErrors(nPtBins,ptBins ,errorsNeg[i] ,ptBinsErr ,errorsNegErr[i] );
		meanErrors[i] = new TGraphErrors(nPtBins,ptBins ,errorsMean[i] ,ptBinsErr ,errorsMeanErr[i] );
		positiveErrors[i] = new TGraphErrors(nPtBins,ptBins ,errorsPos[i] ,ptBinsErr ,errorsPosErr[i] );
		negativeErrorsCorr[i] = new TGraphErrors(nPtBins,ptBins ,errorsNegCorr[i] ,ptBinsErr ,errorsNegErrCorr[i] );
		meanErrorsCorr[i] = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorr[i] ,ptBinsErr ,errorsMeanErrCorr[i] );
		positiveErrorsCorr[i] = new TGraphErrors(nPtBins,ptBins ,errorsPosCorr[i] ,ptBinsErr ,errorsPosErrCorr[i] );
	}
	
	
	for (Int_t l = 0; l < nPtBins; l++){
		errorsPosSummed[l] = pow(errorsPosSummed[l],0.5);
		errorsMeanSummed[l] = pow(errorsMeanSummed[l],0.5);
		errorsPosErrSummed[l] = errorsPosSummed[l]*0.001;
		errorsMeanErrSummed[l] = errorsMeanSummed[l]*0.001;
		errorsNegSummed[l] = -pow(errorsNegSummed[l],0.5);
		errorsNegErrSummed[l] = errorsNegSummed[l]*0.001;
// 		cout << nameCutVariationSC[6].Data() << endl;
		// add PCM & EMCal material errors
		errorsPosCorrMatSummed[l] = pow(errorsPosCorrSummed[l]+  pow(errorsPosCorr[6][l],2) ,0.5);
		errorsMeanCorrMatSummed[l] = pow(errorsMeanCorrSummed[l]+  pow(errorsMeanCorr[6][l],2),0.5);
		errorsNegCorrMatSummed[l] = -pow(errorsNegCorrSummed[l]+  pow(errorsNegCorr[6][l],2),0.5);
		
		errorsPosCorrSummed[l] = pow(errorsPosCorrSummed[l],0.5);
		errorsMeanCorrSummed[l] = pow(errorsMeanCorrSummed[l],0.5);
		errorsPosErrCorrSummed[l] = errorsPosCorrSummed[l]*0.001;
		errorsMeanErrCorrSummed[l] = errorsMeanCorrSummed[l]*0.001;
		errorsMeanErrCorrMatSummed[l] = errorsMeanCorrMatSummed[l]*0.001;
		errorsNegCorrSummed[l] = -pow(errorsNegCorrSummed[l],0.5);
		errorsNegErrCorrSummed[l] = errorsNegCorrSummed[l]*0.001;
// 		cout << errorsMeanSummed[l] << "\t" << errorsMeanCorrMatSummed[l]<< endl;
	}
	
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
	if (meson.Contains("Pi0") ){
		histo2DSysErrMean = new TH2D("histo2DSysErrMean", "histo2DSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,0.,30.);
	} else {
		histo2DSysErrMean = new TH2D("histo2DSysErrMean", "histo2DSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,0.,65.);
	}	
	histo2DSysErrMean->SetYTitle("mean systematic Err %");
	histo2DSysErrMean->SetXTitle("#it{p}_{T} (GeV/#it{c})");
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
	
	TLegend* legendMean; 
	if (meson.Contains("Pi0")){
		legendMean= new TLegend(0.15,0.6,0.57,0.95);
	} else {
		legendMean= new TLegend(0.20,0.6,0.62,0.95);
	}
	legendMean->SetTextSize(0.035);
	legendMean->SetFillColor(0);
	legendMean->SetBorderSize(0);
	if (numberCutStudies> 9) legendMean->SetNColumns(2);
	for(Int_t i = 0; i< numberCutStudies ; i++){
// 		if (additionalNameOutput.CompareTo("EMC1") == 0 && nameCutVariation[i].CompareTo("clusterizer") == 0){
// 		  cout << "jump clusterizer contribution for EMC1" << endl;	
// 		} else {
			DrawGammaSetMarkerTGraphErr(meanErrors[i], 20+i, 1.,color[i],color[i]);
			meanErrors[i]->Draw("pE0,csame");
			legendMean->AddEntry(meanErrors[i],nameCutVariation[i].Data(),"p");
// 		}	
	}
	legendMean->Draw();
	canvasSysErrMean->Update();
	canvasSysErrMean->SaveAs(Form("SystematicErrorsCalculatedCalo/SysMean_%s_%s%s_%s.%s",meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));
	
	delete canvasSysErrMean;
	
	TCanvas* canvasNewSysErrMean = new TCanvas("canvasNewSysErrMean","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasNewSysErrMean, 0.08, 0.01, 0.015, 0.09);
	
	TH2D *histo2DNewSysErrMean ;
	if (meson.Contains("Pi0")){
		histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "histo2DNewSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,30.);
	} else { 
		histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "histo2DNewSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,65.);
	}
	histo2DNewSysErrMean->SetYTitle("mean smoothed systematic Err %");
	histo2DNewSysErrMean->SetXTitle("#it{p}_{T} (GeV/#it{c})");
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
	if (meson.Contains("Pi0")){
		legendMeanNew= new TLegend(0.12,0.6,0.67,0.95);
	} else {
		legendMeanNew= new TLegend(0.23,0.65,0.72,0.95);
	}	
	legendMeanNew->SetTextSize(0.035);
	legendMeanNew->SetMargin(0.1);
	legendMeanNew->SetFillColor(0);
	legendMeanNew->SetFillStyle(0);
	legendMeanNew->SetBorderSize(0);
	if (numberCutStudies> 9) legendMeanNew->SetNColumns(2);
	for(Int_t i = 0; i< numberCutStudies ; i++){
// 		if (additionalNameOutput.CompareTo("EMC1") == 0 && nameCutVariation[i].CompareTo("clusterizer") == 0){
// 		  cout << "jump clusterizer contribution for EMC1" << endl;	
// 		} else {
			DrawGammaSetMarkerTGraphErr(meanErrorsCorr[i], 20+i, 1.,color[i],color[i]);
			meanErrorsCorr[i]->Draw("pX0,csame");
			legendMeanNew->AddEntry(meanErrorsCorr[i],nameCutVariation[i].Data(),"p");	
// 		}	
	}

	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kBlack,kBlack);
	meanErrorsCorrSummedIncMat->Draw("p,csame");
	legendMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"quad. sum.","p");
	legendMeanNew->Draw();
	
	TLatex *labelMeson;
	if (meson.Contains("Pi0")){
		labelMeson= new TLatex(0.75,0.89,Form("#pi^{0} #rightarrow #gamma_{EMC}#gamma_{EMC}"));
	} else {
		labelMeson= new TLatex(0.75,0.89,Form("#eta #rightarrow #gamma_{EMC}#gamma_{EMC}"));
	}
	SetStyleTLatex( labelMeson, 0.038,4);
	labelMeson->Draw();

	TLatex *labelTrig;
	if (additionalNameOutput.CompareTo("")==0){
		labelTrig= new TLatex(0.75,0.84,Form("MB LHC11a"));
	} else if (additionalNameOutput.CompareTo("EMC1")==0){
		labelTrig= new TLatex(0.75,0.84,Form("EMC1 LHC11a"));
	}
	SetStyleTLatex( labelTrig, 0.038,4);
	labelTrig->Draw();
	
	
	TLatex *labelCentrality = new TLatex(0.75,0.93,Form("%s",collisionSystem.Data()	));
	SetStyleTLatex( labelCentrality, 0.038,4);
	labelCentrality->Draw();

	canvasNewSysErrMean->Update();
	canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculatedCalo/SysMeanNewWithMean_%s_%s%s_%s.%s",meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));
	
	
	canvasNewSysErrMean->cd();
	histo2DNewSysErrMean->Draw();

		Int_t cut=0;
		Double_t minPt = 0.8;
		Double_t maxPt = 6;
		if (additionalNameOutput.CompareTo("EMC1") == 0){
			minPt = 2.6;
			maxPt = 16;
		}	
		TF1* pol0 = new TF1("pol0","[0]",minPt,maxPt); //

		TF1* pol1 = new TF1("pol1","[0]+[1]*x",minPt,maxPt); //

		TF1* pol2 = new TF1("pol2","[0]+[1]*x+[2]*x*x",minPt,maxPt); //

		TF1* pol4 = new TF1("pol4","[0]+[1]*x+[2]*x*x+[3]*x*x*x*x",minPt,maxPt); //
// 		pol4->SetParLimits(3,0,10);
// 		pol4->SetParLimits(2,0,10);
		
		meanErrorsCorr[cut]->Fit(pol0,"QNRMEX0+","",minPt,maxPt);
		meanErrorsCorr[cut]->Fit(pol1,"QNRMEX0+","",minPt,maxPt);
		meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",minPt,maxPt);
		meanErrorsCorr[cut]->Fit(pol4,"QNRMEX0+","",minPt,maxPt);
		
		cout << WriteParameterToFile(pol0) << endl;
		cout << WriteParameterToFile(pol1) << endl;
		cout << WriteParameterToFile(pol2) << endl;
		cout << WriteParameterToFile(pol4) << endl;
		
		pol0->SetLineColor(kBlack);
		pol1->SetLineColor(kRed+2);
		pol2->SetLineColor(kBlue+2);
		pol4->SetLineColor(kGreen+2);
		
		DrawGammaSetMarkerTGraphErr(meanErrorsCorr[cut], 20+cut, 1.,color[cut],color[cut]);
		meanErrorsCorr[cut]->Draw("p,csame");
		pol0->Draw("same");
		pol1->Draw("same");
		pol2->Draw("same");
		pol4->Draw("same");
		
	canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculatedCalo/SysMeanNewWithMeanSingle_%s_%s%s_%s_Cut%i.%s",meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data(),cut,suffix.Data()));
	
	

	const char *SysErrDatname = Form("SystematicErrorsCalculatedCalo/SystematicError_%s_%s%s_%s.dat",meson.Data(),energy.Data(),additionalNameOutput.Data(),dateForOutput.Data());
	fstream SysErrDat;
	cout << SysErrDatname << endl;
	SysErrDat.open(SysErrDatname, ios::out);
	for (Int_t l=0; l< nPtBins; l++){
		SysErrDat <<ptBins[l] << "\t" << errorsNegCorrMatSummed[l] << "\t" <<errorsPosCorrMatSummed[l] << "\t"  <<errorsNegCorrSummed[l] << "\t" <<errorsPosCorrSummed[l]  << endl;
	}

	SysErrDat.close();

	const char *SysErrDatnameMean = Form("SystematicErrorsCalculatedCalo/SystematicErrorAveraged_%s_%s%s_%s.dat",meson.Data(),energy.Data(),additionalNameOutput.Data(),dateForOutput.Data());
	fstream SysErrDatAver;
	cout << SysErrDatnameMean << endl;
	SysErrDatAver.open(SysErrDatnameMean, ios::out);
	for (Int_t l=0; l< nPtBins; l++){
		SysErrDatAver << ptBins[l] << "\t-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
	}

	SysErrDatAver.close();
	
	Double_t errorsMeanCorrPID[nPtBins];	
	Double_t errorsMeanCorrSignalExtraction[nPtBins];	
	Double_t errorsMeanCorrTrackReco[nPtBins];	
	Double_t errorsMeanCorrPhotonReco[nPtBins];	
	Double_t errorsMeanCorrClusterDescription[nPtBins];	
	
	for (Int_t l=0; l< nPtBins; l++){
			// YieldExtraction	X i=0
			// NonLinCor		X i=3
			// MinEClus		X i=1
			// M02			X i=4
			// MinNCell		X i=2
			// TRD			X i=5
			// Rapidity		X i=6
			// Alpha		X i=7
		errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]+errorsMeanCorr[12][l]*errorsMeanCorr[12][l]);
		// 1 - 2 - 4 - 5 - 
		errorsMeanCorrClusterDescription[l] = TMath::Sqrt(errorsMeanCorr[1][l]*errorsMeanCorr[1][l]+errorsMeanCorr[2][l]*errorsMeanCorr[2][l]
															+errorsMeanCorr[4][l]*errorsMeanCorr[4][l]+errorsMeanCorr[5][l]*errorsMeanCorr[5][l]);
// 		errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]+errorsMeanCorr[12][l]*errorsMeanCorr[12][l]);
		// cluster description in MC: ClusterMinEnergy 1, ClusterNCells 2, ClusterM02 5, ClusterTrackMatching 4, Cluster timing 8
// 		errorsMeanCorrClusterDescription[l] = TMath::Sqrt(errorsMeanCorr[1][l]*errorsMeanCorr[1][l]+errorsMeanCorr[2][l]*errorsMeanCorr[2][l]
// 															+errorsMeanCorr[4][l]*errorsMeanCorr[4][l]+errorsMeanCorr[5][l]*errorsMeanCorr[5][l]
// 															+errorsMeanCorr[8][l]*errorsMeanCorr[8][l]);
	}
		
		
	TGraphErrors* meanErrorsSignalExtraction = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSignalExtraction ,ptBinsErr ,errorsMeanErrCorrSummed );
	TGraphErrors* meanErrorsClusterDescrip = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrClusterDescription ,ptBinsErr ,errorsMeanErrCorrSummed );
	
    cout << "here" << endl;
   
   
	TCanvas* canvasSummedErrMean = new TCanvas("canvasSummedErrMean","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasSummedErrMean, 0.08, 0.01, 0.015, 0.09);
	TH2D *histo2DSummedErrMean ;
	
	if (meson.Contains("Pi0") ){	
		histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "histo2DSummedErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,25.);
	} else {
		histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "histo2DSummedErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,65.);
	}
	histo2DSummedErrMean->SetYTitle("mean smoothed systematic Err %");
	histo2DSummedErrMean->SetXTitle("#it{p}_{T} (GeV/#it{c})");
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
	if (meson.Contains("Pi0")){
		legendSummedMeanNew= new TLegend(0.13,0.7,0.65,0.95);
	} else {
		legendSummedMeanNew= new TLegend(0.20,0.7,0.7,0.95);
	}
	legendSummedMeanNew->SetTextSize(0.035);
	legendSummedMeanNew->SetFillColor(0);
	legendSummedMeanNew->SetFillStyle(0);	
	legendSummedMeanNew->SetBorderSize(0);
	legendSummedMeanNew->SetNColumns(2);
	legendSummedMeanNew->SetMargin(0.1);
	
	// Non linearity
	if (numberCutStudies>3){
		DrawGammaSetMarkerTGraphErr(meanErrorsCorr[3], 25, 1.,color[8],color[8]);
		meanErrorsCorr[3]->Draw("pX0,csame");
	}
	// Material infront of EMCAL
	if (numberCutStudies>6){
		DrawGammaSetMarkerTGraphErr(meanErrorsCorr[6], 21, 1.,color[7],color[7]);
		meanErrorsCorr[6]->Draw("pX0,csame");
	}
	// trigger systematics
	if ( numberCutStudies>13 && additionalNameOutput.CompareTo("EMC1") == 0 ) {
		DrawGammaSetMarkerTGraphErr(meanErrorsCorr[13], 28, 1.,color[5],color[5]);
		meanErrorsCorr[13]->Draw("pX0,csame");
	}	
	// Clusterizer
	if (numberCutStudies>7){
		DrawGammaSetMarkerTGraphErr(meanErrorsCorr[7], 24, 1.,color[1],color[1]);
		meanErrorsCorr[7]->Draw("pX0,csame");
	}	
	// Feed down correction
	if (numberCutStudies>9){
		DrawGammaSetMarkerTGraphErr(meanErrorsCorr[9], 28, 1.2,color[10],color[10]);
		meanErrorsCorr[9]->Draw("pX0,csame");
	}
	// Material pion
	if (numberCutStudies>11){
		DrawGammaSetMarkerTGraphErr(meanErrorsCorr[11], 29, 1.5,color[3],color[3]);
		meanErrorsCorr[11]->Draw("pX0,csame");
	}
	// Energy scale
	if (numberCutStudies>10){
		DrawGammaSetMarkerTGraphErr(meanErrorsCorr[10], 33, 2.,color[4],color[4]);
		meanErrorsCorr[10]->Draw("pX0,csame");
	}
	
	// Cluster description in MC
	DrawGammaSetMarkerTGraphErr(meanErrorsClusterDescrip, 22, 1.,color[2],color[2]);
	meanErrorsClusterDescrip->Draw("pX0,csame");
	// Signal extraction error
	DrawGammaSetMarkerTGraphErr(meanErrorsSignalExtraction, 20, 1.,color[0],color[0]);
	meanErrorsSignalExtraction->Draw("pX0,csame");
	
	legendSummedMeanNew->AddEntry(meanErrorsSignalExtraction,"Signal Extraction","p");
	legendSummedMeanNew->AddEntry(meanErrorsClusterDescrip,"Cluster Description","p");
	if (numberCutStudies>3)		legendSummedMeanNew->AddEntry(meanErrorsCorr[3],"Cluster Energy Calibration","p");
	if (numberCutStudies>10)	legendSummedMeanNew->AddEntry(meanErrorsCorr[10],"Energy Scale","p");
	if (numberCutStudies>7)		legendSummedMeanNew->AddEntry(meanErrorsCorr[7],"Clusterizer","p");
	if (numberCutStudies>9)		legendSummedMeanNew->AddEntry(meanErrorsCorr[9],"Feed down correction","p");
	if (numberCutStudies>11)	legendSummedMeanNew->AddEntry(meanErrorsCorr[11],"Material #pi^{0}","p");
	if (numberCutStudies>6)		legendSummedMeanNew->AddEntry(meanErrorsCorr[6],"Mat. infront of EMCal","p");
	if (numberCutStudies>13 && additionalNameOutput.CompareTo("EMC1") == 0 ) 
								legendSummedMeanNew->AddEntry( meanErrorsCorr[13],"Trigger Normalization","p");
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kBlack,kBlack);
	meanErrorsCorrSummedIncMat->Draw("p,csame");
	legendSummedMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"quad. sum.","p");
	legendSummedMeanNew->Draw();
	
	labelMeson->Draw();
	labelCentrality->Draw();
	labelTrig->Draw();
	
	canvasSummedErrMean->Update();
	canvasSummedErrMean->SaveAs(Form("SystematicErrorsCalculatedCalo/SysErrorSummedVisu_%s_%s%s_%s.%s",meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));
	
	delete canvasSummedErrMean;
	
}