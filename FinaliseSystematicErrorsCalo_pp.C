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

void FinaliseSystematicErrorsCalo_pp(const char* nameDataFileErrors ="", TString energy="", TString meson = "", Int_t numberOfPtBins =1 ,Int_t numberCutStudies=1, Int_t offSetBeginning = 0, TString additionalName = "pp",TString additionalNameOutput = "",TString suffix = "eps", Int_t mode = 2){
	
	StyleSettingsThesis();	
	SetPlotStyle();
	
	TString date = ReturnDateString();
	TString dateForOutput = ReturnDateStringForOutput();
	TString collisionSystem= ReturnFullCollisionsSystem(energy);
	
	Color_t color[20] = {860,894,807,880,418,403,802,923,634,432,404,435,420,407,416,830,404,608,kCyan-2,1};
	
	Int_t 	numberOfEntriesPos = 0;
	Int_t 	numberOfEntriesNeg = 0;
	
	TFile* fileErrorInput= new TFile(nameDataFileErrors);
	const Int_t nPtBins = numberOfPtBins;
	const Int_t nCuts = numberCutStudies;
	Double_t* ptBins;
	Double_t* ptBinsErr;
	
	TString nameCutVariation2760GeV[7] = {"Yield extraction", "min E_{cluster}", "min # cells", "clu. energy calibration", "V0 tr. match. to cl.", "M_{02}", "cluster timing"};
	TString nameCutVariation[7];
	TString nameCutVariationSC[7];
	
	TString nameCutVariationSC2760GeV[7] = {"YieldExtraction", "ClusterMinEnergy", "ClusterNCells", "NonLinearity", "TrackMatching", "ClusterM02", "Timing"};
		
	if (energy.CompareTo("2.76TeV") == 0) {
		for (Int_t i = 0; i < numberCutStudies; i++){
			nameCutVariation[i] = nameCutVariation2760GeV[i];
			nameCutVariationSC[i] = nameCutVariationSC2760GeV[i];
		}
	}
	
	gSystem->Exec("mkdir -p SystematicErrorsCalculatedConvCalo");
	
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
			TString nameGraphPos = Form("%s_SystErrorRelPos_%s_%s",meson.Data(),nameCutVariationSC[i].Data(),additionalName.Data()	);
			TString nameGraphNeg = Form("%s_SystErrorRelNeg_%s_%s",meson.Data(),nameCutVariationSC[i].Data(),additionalName.Data()	);
			cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
			graphPosErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
			graphNegErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
		} else {
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

		if (nameCutVariationSC[i].CompareTo("NonLinearity")==0 ){ //&& meson.Contains("Pi0")
			for (Int_t k = 1; k < nPtBins; k++){
				cout << ptBins[k] << endl;
				Double_t error = 2.;
// 				Double_t error = 0.79-0.21*ptBins[k]+0.17*ptBins[k]*ptBins[k]-0.00058*ptBins[k]*ptBins[k]*ptBins[k]*ptBins[k]; // parametrisation with No NonLinearity in 
				errorsMean[i][k] = error;
				errorsMeanErr[i][k] = error*0.01;
				errorsMeanCorr[i][k] = error;
				errorsMeanErrCorr[i][k] = error*0.01;
			}
		}
		if (nameCutVariationSC[i].CompareTo("ClusterM02")==0 ){ //&& meson.Contains("Pi0")
			for (Int_t k = 0; k < nPtBins; k++){
				Double_t error = 2.2;
// 				Double_t error = 0.598+(0.285)*ptBins[k]-(1.08724e-01)*ptBins[k]*ptBins[k]+(4.26047e-03)*ptBins[k]*ptBins[k]*ptBins[k]*ptBins[k]; // parametrisation
				errorsMean[i][k] = error;
				errorsMeanErr[i][k] = error*0.01;
				errorsMeanCorr[i][k] = error;
				errorsMeanErrCorr[i][k] = error*0.01;
			}
		}
		if (nameCutVariationSC[i].CompareTo("ClusterMaterialTRD")==0 ){
			for (Int_t k = 0; k < nPtBins; k++){
				errorsMean[i][k] = 4.24; //(5% for TRD mat, 5% for TOF mat added in quadrature)
				errorsMeanErr[i][k] = 0.04;
				errorsMeanCorr[i][k] = 4.24;
				errorsMeanErrCorr[i][k] = 0.04;
			}	
		}

		if (nameCutVariationSC[i].CompareTo("ClusterMinEnergy")==0  ){
			if(meson.Contains("Pi0")){
				for (Int_t k = 0; k < nPtBins; k++){
					errorsMean[i][k] = 2.;
					errorsMeanErr[i][k] = 0.015;
					errorsMeanCorr[i][k] = 2.;
					errorsMeanErrCorr[i][k] = 0.015;
				}
			} else {
				for (Int_t k = 0; k < nPtBins; k++){
					errorsMean[i][k] = 5.;
					errorsMeanErr[i][k] = 0.05;
					errorsMeanCorr[i][k] = 5.;
					errorsMeanErrCorr[i][k] = 0.05;
				}
				
			}	
		}

		if (nameCutVariationSC[i].CompareTo("ClusterNCells")==0 ){ //&& meson.Contains("Pi0")
			for (Int_t k = 0; k < nPtBins; k++){
				Double_t error = 0.91+(0.255)*ptBins[k];
				errorsMean[i][k] = error;
				errorsMeanErr[i][k] = 0.032;
				errorsMeanCorr[i][k] = error;
				errorsMeanErrCorr[i][k] = 0.032;
			}	
		}

		if (nameCutVariationSC[i].CompareTo("ConvPhi")==0 && meson.Contains("Eta") && !meson.Contains("Pi0")){
			cout << "here" << endl;
			for (Int_t k = 0; k < nPtBins; k++){
				Double_t error = 15;
				errorsMean[i][k] = error;
				errorsMeanErr[i][k] = 0.01*error;
				errorsMeanCorr[i][k] = error;
				errorsMeanErrCorr[i][k] = 0.01*error;
			}	
		}

		if (nameCutVariationSC[i].CompareTo("Alpha")==0 && meson.Contains("Pi0")){
			for (Int_t k = 2; k < nPtBins; k++){
				Double_t error = 2.67140e-01 ;//+5.61581e-02*ptBins[k]*ptBins[k]+7.93500e-04*ptBins[k]*ptBins[k]*ptBins[k]*ptBins[k];
				errorsMean[i][k] = error;
				errorsMeanErr[i][k] = 0.01*error;
				errorsMeanCorr[i][k] = error;
				errorsMeanErrCorr[i][k] = 0.01*error;
			}	
		}
		if (nameCutVariationSC[i].CompareTo("Qt")==0 && meson.Contains("Eta") && !meson.Contains("Pi0")){
			for (Int_t k = 0; k < nPtBins; k++){
				Double_t error = 6;
				errorsMean[i][k] = error;
				errorsMeanErr[i][k] = 0.01*error;
				errorsMeanCorr[i][k] = error;
				errorsMeanErrCorr[i][k] = 0.01*error;
			}	
		}
		if (nameCutVariationSC[i].CompareTo("Chi2")==0 && meson.Contains("Eta") && !meson.Contains("Pi0")){
			for (Int_t k = 0; k < nPtBins; k++){
				Double_t error = 4;
				errorsMean[i][k] = error;
				errorsMeanErr[i][k] = 0.01*error;
				errorsMeanCorr[i][k] = error;
				errorsMeanErrCorr[i][k] = 0.01*error;
			}	
		}
		if (nameCutVariationSC[i].CompareTo("dEdxE")==0 && meson.Contains("Eta") && !meson.Contains("Pi0")){
			for (Int_t k = 0; k < nPtBins; k++){
				Double_t error = 5;
				errorsMean[i][k] = error;
				errorsMeanErr[i][k] = 0.01*error;
				errorsMeanCorr[i][k] = error;
				errorsMeanErrCorr[i][k] = 0.01*error;
			}	
		}
		if (nameCutVariationSC[i].CompareTo("dEdxPi")==0 && meson.Contains("Eta") && !meson.Contains("Pi0")){
			for (Int_t k = 0; k < nPtBins; k++){
				Double_t error = 5;
				errorsMean[i][k] = error;
				errorsMeanErr[i][k] = 0.01*error;
				errorsMeanCorr[i][k] = error;
				errorsMeanErrCorr[i][k] = 0.01*error;
			}	
		}
		if (nameCutVariationSC[i].CompareTo("ClusterTrackMatching")==0 && meson.Contains("Eta") && !meson.Contains("Pi0")){
			for (Int_t k = 0; k < nPtBins; k++){
				Double_t error = 10;
				errorsMean[i][k] = error;
				errorsMeanErr[i][k] = 0.01*error;
				errorsMeanCorr[i][k] = error;
				errorsMeanErrCorr[i][k] = 0.01*error;
			}	
		}


		if (nameCutVariationSC[i].CompareTo("SinglePt")==0 && meson.Contains("Eta") && !meson.Contains("Pi0")){
			for (Int_t k = 0; k < nPtBins; k++){
				Double_t error = 10;
				errorsMean[i][k] = error;
				errorsMeanCorr[i][k] = error;
				errorsMeanErr[i][k] = 0.01*error;
				errorsMeanErrCorr[i][k] = 0.01*error;
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
	
	Double_t errorMaterial = 4.50;
	
	for (Int_t l = 0; l < nPtBins; l++){
		errorsPosSummed[l] = pow(errorsPosSummed[l],0.5);
		errorsMeanSummed[l] = pow(errorsMeanSummed[l],0.5);
		errorsPosErrSummed[l] = errorsPosSummed[l]*0.001;
		errorsMeanErrSummed[l] = errorsMeanSummed[l]*0.001;
		errorsNegSummed[l] = -pow(errorsNegSummed[l],0.5);
		errorsNegErrSummed[l] = errorsNegSummed[l]*0.001;
// 		cout << nameCutVariationSC[14].Data() << endl;
		// add PCM & EMCal material errors
		errorsPosCorrMatSummed[l] = pow(errorsPosCorrSummed[l]+ pow(errorMaterial ,2.) + pow(errorsPosCorr[14][l],2) ,0.5);
		errorsMeanCorrMatSummed[l] = pow(errorsMeanCorrSummed[l]+ pow(errorMaterial ,2.)+ pow(errorsMeanCorr[14][l],2),0.5);
		errorsNegCorrMatSummed[l] = -pow(errorsNegCorrSummed[l]+ pow(errorMaterial ,2.)+ pow(errorsNegCorr[14][l],2),0.5);
		
		errorsPosCorrSummed[l] = pow(errorsPosCorrSummed[l],0.5);
		errorsMeanCorrSummed[l] = pow(errorsMeanCorrSummed[l],0.5);
		errorsPosErrCorrSummed[l] = errorsPosCorrSummed[l]*0.001;
		errorsMeanErrCorrSummed[l] = errorsMeanCorrSummed[l]*0.001;
		errorsMeanErrCorrMatSummed[l] = errorsMeanCorrMatSummed[l]*0.001;
		errorsNegCorrSummed[l] = -pow(errorsNegCorrSummed[l],0.5);
		errorsNegErrCorrSummed[l] = errorsNegCorrSummed[l]*0.001;
// 		cout << errorsMeanSummed[l] << "\t" << errorsMeanCorrMatSummed[l]<< endl;
	}
	Double_t ptBinsMaterial[nPtBins];
	Double_t errorsMat[nPtBins];
	for (Int_t l = 0; l < nPtBins; l++){
		errorsMat[l] = errorMaterial;
		ptBinsMaterial[l] = ptBins[l]+0.05;
		
	}
	TGraphErrors* graphMaterialError = new TGraphErrors(nPtBins,ptBinsMaterial ,errorsMat ,ptBinsErr ,errorsMeanErrSummed );
	
	cout << __LINE__ << " here" <<  endl;
	negativeErrorsSummed = new TGraphErrors(nPtBins,ptBins ,errorsNegSummed ,ptBinsErr ,errorsNegErrSummed );
	negativeErrorsCorrSummed = new TGraphErrors(nPtBins,ptBins ,errorsNegCorrSummed ,ptBinsErr ,errorsNegErrCorrSummed );
	positiveErrorsSummed = new TGraphErrors(nPtBins,ptBins ,errorsPosSummed ,ptBinsErr ,errorsPosErrSummed );
	positiveErrorsCorrSummed = new TGraphErrors(nPtBins,ptBins ,errorsPosCorrSummed ,ptBinsErr ,errorsPosErrCorrSummed );
	meanErrorsSummed = new TGraphErrors(nPtBins,ptBins ,errorsMeanSummed ,ptBinsErr ,errorsMeanErrSummed );
	meanErrorsCorrSummed = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSummed ,ptBinsErr ,errorsMeanErrCorrSummed );
	meanErrorsCorrSummedIncMat = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrMatSummed ,ptBinsErr ,errorsMeanErrCorrMatSummed );

	
	
	cout << __LINE__ << " here" << endl;
	
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
		DrawGammaSetMarkerTGraphErr(meanErrors[i], 20+i, 1.,color[i],color[i]);
		meanErrors[i]->Draw("pE0,csame");
		legendMean->AddEntry(meanErrors[i],nameCutVariation[i].Data(),"p");
	}
	DrawGammaSetMarkerTGraphErr(graphMaterialError, 24, 1.,color[10],color[10]);
	graphMaterialError->Draw("pX0,csame");
	legendMean->AddEntry(graphMaterialError,"Inner Material","p");

/*	
	DrawGammaSetMarkerTGraphErr(meanErrorsSummed, 20, 1.,1,1);
	meanErrorsSummed->Draw("pE0,csame");
	legendMean->AddEntry(meanErrorsSummed,"quad. sum.","p");*/
	legendMean->Draw();
	canvasSysErrMean->Update();
	canvasSysErrMean->SaveAs(Form("SystematicErrorsCalculatedConvCalo/SysMean_%s_%s%s_%s.%s",meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));
	
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
		DrawGammaSetMarkerTGraphErr(meanErrorsCorr[i], 20+i, 1.,color[i],color[i]);
		meanErrorsCorr[i]->Draw("pX0,csame");
		legendMeanNew->AddEntry(meanErrorsCorr[i],nameCutVariation[i].Data(),"p");	
	}

	DrawGammaSetMarkerTGraphErr(graphMaterialError, 24, 1.,color[10],color[10]);
	graphMaterialError->Draw("pX0,csame");
	legendMeanNew->AddEntry(graphMaterialError,"Inner Material","p");
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kBlack,kBlack);
	meanErrorsCorrSummedIncMat->Draw("p,csame");
	legendMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"quad. sum.","p");
	legendMeanNew->Draw();
	
	TLatex *labelMeson;
	if (meson.Contains("Pi0")){
		labelMeson= new TLatex(0.75,0.89,Form("#pi^{0} #rightarrow #gamma_{conv}#gamma_{calo}"));
	} else {
		labelMeson= new TLatex(0.75,0.89,Form("#eta #rightarrow #gamma_{conv}#gamma_{calo}"));
	}
	SetStyleTLatex( labelMeson, 0.038,4);
	labelMeson->Draw();

	TLatex *labelCentrality = new TLatex(0.75,0.93,Form("%s",collisionSystem.Data()	));
	SetStyleTLatex( labelCentrality, 0.038,4);
	labelCentrality->Draw();

	canvasNewSysErrMean->Update();
	canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculatedConvCalo/SysMeanNewWithMean_%s_%s%s_%s.%s",meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));
	
	
	canvasNewSysErrMean->cd();
	histo2DNewSysErrMean->Draw();

		Int_t cut=4;
		Double_t minPt = 0.6;
		Double_t maxPt = ptBins[nPtBins-2]+1;
// 		Double_t maxPt = 3.6;
// 		TF1* pol4 = new TF1("pol2","[0]+[1]*x+[2]*x*x",minPt,maxPt); //
		TF1* pol4 = new TF1("pol4","[0]+[1]*x+[2]*x*x+[3]*x*x*x*x",minPt,maxPt); //
		pol4->SetParLimits(3,0,10);
// 		pol4->FixParameter(2,0);
// 		pol4->SetParLimits(1,0,1);
// 		pol4->FixParameter(0,0.8);
		
		meanErrorsCorr[cut]->Fit(pol4,"NRMEX0+","",minPt,maxPt);
// 		pol4->SetParameter(0,0.598);
// 		pol4->SetParameter(1,0.285);
// 		pol4->SetParameter(2,0.112);
// 		pol4->SetParameter(3,-0.0004);
		pol4->SetLineColor(kRed+2);
		
		
		
		DrawGammaSetMarkerTGraphErr(meanErrorsCorr[cut], 20+cut, 1.,color[cut],color[cut]);
		meanErrorsCorr[cut]->Draw("p,csame");
		pol4->Draw("same");
		
	canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculatedConvCalo/SysMeanNewWithMeanSingle_%s_%s%s_%s.%s",meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));
	
	
// 	delete canvasNewSysErrMean;

		const char *SysErrDatname = Form("SystematicErrorsCalculatedConvCalo/SystematicError_%s_%s%s_%s.dat",meson.Data(),energy.Data(),additionalNameOutput.Data(),dateForOutput.Data());
		fstream SysErrDat;
		cout << SysErrDatname << endl;
		SysErrDat.open(SysErrDatname, ios::out);
		for (Int_t l=0; l< nPtBins; l++){
			SysErrDat <<errorsNegCorrMatSummed[l] << "\t" <<errorsPosCorrMatSummed[l] << "\t"  <<errorsNegCorrSummed[l] << "\t" <<errorsPosCorrSummed[l]  << endl;
		}

		SysErrDat.close();

		const char *SysErrDatnameMean = Form("SystematicErrorsCalculatedConvCalo/SystematicErrorAveraged_%s_%s%s_%s.dat",meson.Data(),energy.Data(),additionalNameOutput.Data(),dateForOutput.Data());
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
		Double_t errorsMeanCorrClusterDescription[nPtBins];	
		
		for (Int_t l=0; l< nPtBins; l++){
// 			"YieldExtraction"-0,"dEdxE"-1,"dEdxPi"-2, "TPCCluster"-3, "SinglePt"-4, "Chi2"-5, "Qt"-6, "Alpha"-7, "ConvPhi"-8, "ClusterMinEnergy"-9, "ClusterNCells"-10, "NonLinearity"-11, "ClusterTrackMatching" -12, "ClusterM02" -13, "ClusterMaterialTRD" -14
			// grouping:
			// Signal extraction: Yield extraction 0, Alpha 7 
			if (numberCutStudies>8){
				errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]+errorsMeanCorr[7][l]*errorsMeanCorr[7][l]+2.5*2.5);	
			} else {
				errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]+2.5*2.5);
			}
			// PID: dEdxE 1, dEdxPi 2
			errorsMeanCorrPID[l] =TMath::Sqrt(errorsMeanCorr[1][l]*errorsMeanCorr[1][l]+ errorsMeanCorr[2][l]*errorsMeanCorr[2][l]);
			// photon reco: Chi2 5, Qt 6
			errorsMeanCorrPhotonReco[l] =TMath::Sqrt( errorsMeanCorr[5][l]* errorsMeanCorr[5][l]+errorsMeanCorr[6][l]*errorsMeanCorr[6][l]);
			// track reconstruction: TPCCluster 3, Single pt 4, ConvPhi 8
			errorsMeanCorrTrackReco[l] = TMath::Sqrt(errorsMeanCorr[3][l]*errorsMeanCorr[3][l]+errorsMeanCorr[4][l]*errorsMeanCorr[4][l]+errorsMeanCorr[8][l]*errorsMeanCorr[8][l]);				
			// cluster description in MC: ClusterMinEnergy 9, ClusterNCells 10, ClusterM02 13
			errorsMeanCorrClusterDescription[l] = TMath::Sqrt(errorsMeanCorr[9][l]*errorsMeanCorr[9][l]+errorsMeanCorr[10][l]*errorsMeanCorr[10][l]+errorsMeanCorr[13][l]*errorsMeanCorr[13][l]);
		}
		
		
	TGraphErrors* meanErrorsPID = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrPID ,ptBinsErr ,errorsMeanErrCorrSummed );
	TGraphErrors* meanErrorsPhotonReco = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrPhotonReco ,ptBinsErr ,errorsMeanErrCorrSummed );
	TGraphErrors* meanErrorsSignalExtraction = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSignalExtraction ,ptBinsErr ,errorsMeanErrCorrSummed );
	TGraphErrors* meanErrorsTrackReco = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrTrackReco ,ptBinsErr ,errorsMeanErrCorrSummed );
	TGraphErrors* meanErrorsClusterDescrip = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrClusterDescription ,ptBinsErr ,errorsMeanErrCorrSummed );
	
    cout << __LINE__ << " here" << endl;
   
   
	TCanvas* canvasSummedErrMean = new TCanvas("canvasSummedErrMean","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasSummedErrMean, 0.08, 0.01, 0.015, 0.09);
	TH2D *histo2DSummedErrMean ;
	
	if (meson.Contains("Pi0") ){	
		histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "histo2DSummedErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,30.);
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
	
	// Signal extraction error
	DrawGammaSetMarkerTGraphErr(meanErrorsSignalExtraction, 20, 1.,color[0],color[0]);
	meanErrorsSignalExtraction->Draw("pX0,csame");
	// electron PID error
	DrawGammaSetMarkerTGraphErr(meanErrorsPID, 21, 1.,color[1],color[1]);
	meanErrorsPID->Draw("pX0,csame");
	// track reconstruction error
	DrawGammaSetMarkerTGraphErr(meanErrorsTrackReco, 22, 1.,color[2],color[2]);
	meanErrorsTrackReco->Draw("pX0,csame");
	// Photon reco error
	DrawGammaSetMarkerTGraphErr(meanErrorsPhotonReco, 23, 1.,color[3],color[3]);
	meanErrorsPhotonReco->Draw("pX0,csame");
	// Non linearity
	DrawGammaSetMarkerTGraphErr(meanErrorsCorr[11], 25, 1.,color[5],color[5]);
	meanErrorsCorr[11]->Draw("pX0,csame");
	// Material infront of EMCAL
	DrawGammaSetMarkerTGraphErr(meanErrorsCorr[14], 21, 1.,color[7],color[7]);
	meanErrorsCorr[14]->Draw("pX0,csame");
	// PCM Material budget 
	DrawGammaSetMarkerTGraphErr(graphMaterialError, 24, 1.,color[4],color[4]);
	graphMaterialError->Draw("pX0,csame");
	// Track matching V0 to EMCAL
	DrawGammaSetMarkerTGraphErr(meanErrorsCorr[12], 20, 1.,color[18],color[18]);
	meanErrorsCorr[12]->Draw("pX0,csame");
	// Cluster description in MC
	DrawGammaSetMarkerTGraphErr(meanErrorsClusterDescrip, 22, 1.,color[8],color[8]);
	meanErrorsClusterDescrip->Draw("pX0,csame");
	
	legendSummedMeanNew->AddEntry(meanErrorsSignalExtraction,"Signal Extraction","p");
	legendSummedMeanNew->AddEntry(meanErrorsPID,"Electron PID","p");
	legendSummedMeanNew->AddEntry(meanErrorsTrackReco,"Track Reconstruction","p");
	legendSummedMeanNew->AddEntry(meanErrorsPhotonReco,"Photon Reconstruction","p");
    legendSummedMeanNew->AddEntry(graphMaterialError,"Inner Material","p");
	legendSummedMeanNew->AddEntry(meanErrorsClusterDescrip,"Cluster Description","p");
	legendSummedMeanNew->AddEntry(meanErrorsCorr[11],"Cluster Energy Calibration","p");
	legendSummedMeanNew->AddEntry(meanErrorsCorr[12],"V0 tr. match. to cluster","p");
	legendSummedMeanNew->AddEntry(meanErrorsCorr[14],"Mat. infront of EMCal","p");
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kBlack,kBlack);
	meanErrorsCorrSummedIncMat->Draw("p,csame");
	legendSummedMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"quad. sum.","p");
	legendSummedMeanNew->Draw();
	
	labelMeson->Draw();
	labelCentrality->Draw();
	
	canvasSummedErrMean->Update();
	canvasSummedErrMean->SaveAs(Form("SystematicErrorsCalculatedConvCalo/SysErrorSummedVisu_%s_%s%s_%s.%s",meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));
	
	delete canvasSummedErrMean;
// 
// 	
// 	
// 	const char *SysErrDatnameMeanPaper = Form("SystematicErrorsCalculatedConvCalo/SystematicErrorAveragedPaper_%s_%s%s_%s.dat",meson.Data(),energy.Data(),additionalNameOutput.Data(),dateForOutput.Data());
// 	fstream SysErrDatAverPaper;
// 	cout << SysErrDatnameMeanPaper << endl;
// 	SysErrDatAverPaper.open(SysErrDatnameMeanPaper, ios::out);
// 	SysErrDatAverPaper  << "#it{p}_{T}" << "\t Material \t Yield Extraction \t PID \t photon reco \t track recon \t summed" <<  endl;
// 	for (Int_t l=0; l< nPtBins; l++){
// 		SysErrDatAverPaper << ptBins[l] <<"\t" << errorMaterial*2 << "\t" << errorsMeanCorrSignalExtraction[l] << "\t" << errorsMeanCorrPID[l]<< "\t" << errorsMeanCorrPhotonReco[l]<< "\t" <<errorsMeanCorrTrackReco[l] <<"\t" << errorsMeanCorrMatSummed[l]<< endl;
// 	}
// 
// 	SysErrDatAverPaper.close();
	
}