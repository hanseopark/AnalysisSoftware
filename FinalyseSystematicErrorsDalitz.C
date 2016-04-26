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
#include "CommonHeaders/CombinationFunctions.h"

void FinalyseSystematicErrorsDalitz(const char* nameDataFileErrors ="", TString energy="", TString meson = "", Int_t numberOfPtBins =1 ,Int_t numberCutStudies=1, Double_t decisionBoundary=1., Int_t offSetBeginning = 0, TString additionalName = "",TString additionalNameOutput = ""){
	
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


	
	
        TString nameCutVariation[18];
        TString nameCutVariationSC[18];
      
        


	TString nameCutVariation900GeV[14] = {"Yield extraction","dE/dx e-line","dE/dx #pi-line","TPC cluster","TOF variation","Single e^{#pm} p_{t}",  "#chi^{2} #gamma","q_{t,max}","#alpha","R","#psi_{pair}", "V0-finder", "#eta","BG variation"};  
	TString nameCutVariationPbPb_2760GeV[14] = {"Yield extraction","BG variation","dE/dx e-line","dE/dx #pi-line","TPC cluster","TOF variation","Single e^{#pm} p_{t}",  "#chi^{2} #gamma","q_{t,max}","#alpha", "V0-finder","R","#psi_{pair}", "#eta"};  

	
	
	TString nameCutVariationSCPbPb_2760GeV[14] = {"YieldExtraction","BG","dEdxE","dEdxPi", "Cluster", "TOF", "SinglePt", "Chi2", "Qt", "Alpha", "V0Offline", "R", "PsiPair", "Eta"};
	
	
        TString nameCutVariation2760GeV[18]    = {"YieldExtraction","dEdx e-line sec","dEdx #pi-line sec","dEdx e-line prim","Single e^{#pm} p_{t} sec","TPCcls sec","DCAxy","Single e^{#pm} p_{t} prim","TPCcls prim","Chi2","Qt","MassCut","Eta","Alpha","PsiPair + ITS","Background","dEdx #pi-line prim","V0-finder"};
        TString nameCutVariationSC2760GeV[18]  = {"YieldExtraction","dEdxEG","dEdxPiG","dEdxEE","SinglePtG","TPCclsG","DCAxy","SinglePtE","TPCclsE","Chi2","Qt","MassCut","Eta","Alpha","PsiPairITS","Background","dEdxPiE","V0-finder"};
	                                        // {"YieldExtraction","dEdxEG-1","dEdxPiG-2","dEdxEE-3","SinglePtG-4","TPCclsG-5","DCAxy-6","SinglePtE-7","TPCclsE-8","Chi2-9","Qt-10","MassCut-11","Eta-12","Alpha-13","PsiPairITS-14","Background-15","dEdxPiE-16","V0-finder-17"};
						
						
        TString nameCutVariation7TeV[18]           = {"YieldExtraction","dEdx e-line sec","dEdx #pi-line sec","dEdx e-line prim","Single e^{#pm} p_{t} sec","TPCcls sec","DCAxy","Single e^{#pm} p_{t} prim","TPCcls prim","Chi2","Qt","MassCut","Eta","Alpha","PsiPair + ITS","Background","dEdx #pi-line prim","V0-finder"};
        TString nameCutForOutputVariation7TeV[18]  = {"YieldExtraction","dEdxEG","dEdxPiG","dEdxEE","SinglePtG","TPCclsG","DCAxy","SinglePtE","TPCclsE","Chi2","Qt","MassCut","Eta","Alpha","PsiPairITS","Background","dEdxPiE","V0-finder"};

	
	//TString nameCutVariationpPb_5023GeV[18]    = {"YieldExtraction","dEdx e-line sec","dEdx #pi-line sec","dEdx e-line prim","Single e^{#pm} p_{t} sec","TPCcls sec","DCAxy","Single e^{#pm} p_{t} prim","TPCcls prim","Chi2","Qt","MassCut","Eta","Alpha","PsiPair + ITS","Background","dEdx #pi-line prim","V0-finder"};
        //TString nameCutVariationSCpPb_5023GeV[18]  = {"YieldExtraction","dEdxEG","dEdxPiG","dEdxEE","SinglePtG","TPCClsG","DCAxy","SinglePtE","TPCClsE","Chi2","Qt","MassCut","Eta","Alpha","PsiPairITS","Background","dEdxPiE","V0-finder"};
	
        TString nameCutVariationpPb_5023GeV[18]    = {"YieldExtraction","dEdx e-line sec","dEdx #pi-line sec","dEdx e-line prim","Single e^{#pm} p_{t} sec","TPCcls sec","DCAxy","Single e^{#pm} p_{t} prim","TPCcls prim","Chi2","Qt","MassCut","Alpha","PsiPair + ITS","Background","dEdx #pi-line prim","Eta","V0-finder"};
        TString nameCutVariationSCpPb_5023GeV[18]  = {"YieldExtraction","dEdxEG","dEdxPiG","dEdxEE","SinglePtG","TPCClsG","DCAxy","SinglePtE","TPCClsE","Chi2","Qt","MassCut","Alpha","PsiPairITS","Background","dEdxPiE","Eta","V0-finder"};
	enum{kYieldExtraction=0,kdEdxEG,kdEdxPiG,kdEdxEE,kSinglePtG,kTPCClsG,kDCAxy,kSinglePtE,kTPCClsE,kChi2,kQt,kMassCut,kAlpha,kPsiPairITS,kBackground,kdEdxPiE,kEta,kV0finder};                 

                          
	
	
	if (energy.CompareTo("2.76TeV") == 0) {
		for (Int_t i = 0; i < numberCutStudies; i++){
			nameCutVariation[i] = nameCutVariation2760GeV[i];
			nameCutVariationSC[i] = nameCutVariationSC2760GeV[i];
		}	
	} else if (energy.CompareTo("PbPb_2.76TeV") == 0) {
		for (Int_t i = 0; i < numberCutStudies; i++){
			nameCutVariation[i]   = nameCutVariationPbPb_2760GeV[i];
			nameCutVariationSC[i] = nameCutVariationSCPbPb_2760GeV[i];
		}
	}
	else if (energy.CompareTo("7TeV") == 0) {
	  	for (Int_t i = 0; i < numberCutStudies; i++){
			nameCutVariation[i]   = nameCutVariation7TeV[i];
			nameCutVariationSC[i] = nameCutForOutputVariation7TeV[i];
		}
	}
	else if ( energy.CompareTo("pPb_5.023TeV") == 0 ) {
	  
		for (Int_t i = 0; i < numberCutStudies; i++){
			nameCutVariation[i]   = nameCutVariationpPb_5023GeV[i];
			nameCutVariationSC[i] = nameCutVariationSCpPb_5023GeV[i];
		}
	}
	
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
			if (energy.CompareTo("PbPb_2.76TeV") == 0 || energy.CompareTo("pPb_5.023TeV") == 0 ) {
				nameGraphPos = Form("%s_SystErrorRelPos_%s_%s",meson.Data(),nameCutVariationSC[i].Data(),additionalName.Data()	);
				nameGraphNeg = Form("%s_SystErrorRelNeg_%s_%s",meson.Data(),nameCutVariationSC[i].Data(),additionalName.Data()	);
			} else {	
				nameGraphPos = Form("%s_SystErrorRelPos_%s_pp",meson.Data(),nameCutVariationSC[i].Data());
				nameGraphNeg = Form("%s_SystErrorRelNeg_%s_pp",meson.Data(),nameCutVariationSC[i].Data());
			}
			cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
			
			graphPosErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
			graphNegErrors = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
		} else {
			TString nameGraphPos;
			TString nameGraphNeg;
			if (energy.CompareTo("PbPb_2.76TeV") == 0) {
				nameGraphPos = Form("%s_SystErrorRelPos_%s%s",meson.Data(),nameCutVariationSC[i].Data(),additionalNameOutput.Data()	);
				nameGraphNeg = Form("%s_SystErrorRelNeg_%s%s",meson.Data(),nameCutVariationSC[i].Data(),additionalNameOutput.Data()	);
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
			ptBins    = graphNegErrors->GetX();
			ptBinsErr = graphNegErrors->GetEXhigh();
		}

		errorsNeg[i]    = graphNegErrors->GetY();
		errorsNegErr[i] = graphNegErrors->GetEYhigh();
		errorsPos[i]    = graphPosErrors->GetY();
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
		} else {
			CalculateMeanSysErr(errorsMean[i], errorsMeanErr[i], errorsPos[i], errorsNeg[i], nPtBins);	
			CorrectSystematicErrorsWithMean(errorsPos[i],errorsPosErr[i], errorsPosCorr[i], errorsPosErrCorr[i], nPtBins);
			CorrectSystematicErrorsWithMean(errorsNeg[i],errorsNegErr[i], errorsNegCorr[i], errorsNegErrCorr[i], nPtBins);
			CorrectSystematicErrorsWithMean(errorsMean[i], errorsMeanErr[i], errorsMeanCorr[i], errorsMeanErrCorr[i], nPtBins);
		        cout<<"Entro a la media"<<endl;	
		
		}
		for (Int_t l = 0; l < nPtBins; l++){
			errorsPosSummed[l] = errorsPosSummed[l]+pow(errorsPos[i][l],2);
			errorsNegSummed[l] = errorsNegSummed[l]+ pow(errorsNeg[i][l],2);
			errorsMeanSummed[l] = errorsMeanSummed[l]+ pow(errorsMean[i][l],2);
			cout<<"errorsMeanSummed: "<< l <<"  "<<errorsMeanSummed[l]<<endl;
			errorsPosCorrSummed[l] = errorsPosCorrSummed[l]+pow(errorsPosCorr[i][l],2);
			errorsNegCorrSummed[l] = errorsNegCorrSummed[l] +pow(errorsNegCorr[i][l],2);
			errorsMeanCorrSummed[l] =errorsMeanCorrSummed[l]+ pow(errorsMeanCorr[i][l],2);
			cout<<"errorsMeanCorrSummed: "<< l << " "<< errorsMeanCorrSummed[l]<<endl;
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
	
	Double_t errorBRDalitz;
	
	errorBRDalitz = 2.98;
	
	
	for (Int_t l = 0; l < nPtBins; l++){
	  
		errorsPosSummed[l] = pow(errorsPosSummed[l]   + pow( errorBRDalitz,2.) ,0.5);
		errorsMeanSummed[l] = pow(errorsMeanSummed[l] + pow( errorBRDalitz,2.) ,0.5);
		
		errorsPosErrSummed[l] = errorsPosSummed[l]*0.001;
		errorsMeanErrSummed[l] = errorsMeanSummed[l]*0.001;
		
		errorsNegSummed[l] = -pow(errorsNegSummed[l] + pow( errorBRDalitz,2.) ,0.5);
		errorsNegErrSummed[l] = errorsNegSummed[l]*0.001;

		errorsPosCorrMatSummed[l] =  pow(errorsPosCorrSummed[l] + pow( errorBRDalitz,2.) + pow(errorMaterial,2.),0.5);
		errorsMeanCorrMatSummed[l] = pow(errorsMeanCorrSummed[l]+ pow( errorBRDalitz,2.) + pow(errorMaterial,2.),0.5);
		errorsNegCorrMatSummed[l] = -pow(errorsNegCorrSummed[l] + pow( errorBRDalitz,2.) + pow(errorMaterial,2.),0.5);
		
		errorsPosCorrSummed[l] = pow(errorsPosCorrSummed[l]   + pow( errorBRDalitz,2.) ,0.5);
		errorsMeanCorrSummed[l] = pow(errorsMeanCorrSummed[l] + pow( errorBRDalitz,2.) ,0.5);
		errorsPosErrCorrSummed[l] = errorsPosCorrSummed[l]*0.001;
		errorsMeanErrCorrSummed[l] = errorsMeanCorrSummed[l]*0.001;
		errorsMeanErrCorrMatSummed[l] = errorsMeanCorrMatSummed[l]*0.001;
		errorsNegCorrSummed[l] = -pow(errorsNegCorrSummed[l] + pow( errorBRDalitz,2.) ,0.5);
		errorsNegErrCorrSummed[l] = errorsNegCorrSummed[l]*0.001;
// 		cout << errorsMeanSummed[l] << "\t" << errorsMeanCorrMatSummed[l]<< endl;
	}
	
	Double_t errorsMat[nPtBins];
	for (Int_t l = 0; l < nPtBins; l++){
		errorsMat[l] = errorMaterial;
		
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
		histo2DSysErrMean = new TH2D("histo2DSysErrMean", "histo2DSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,0.,50.);
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
	
	TLegend* legendMean = new TLegend(0.1,0.6,0.36,0.95);
	legendMean->SetTextSize(0.030);
	legendMean->SetFillColor(0);
	legendMean->SetBorderSize(0);

	if (numberCutStudies> 9) legendMean->SetNColumns(2);

	for(Int_t i = 0; i< numberCutStudies ; i++){
		DrawGammaSetMarkerTGraphErr(meanErrors[i], 20+i, 1.,color[i],color[i]);
		meanErrors[i]->Draw("p,csame");
		legendMean->AddEntry(meanErrors[i],nameCutVariation[i].Data(),"p");
                cout<<nameCutVariation[i].Data()<<endl;
	}

	DrawGammaSetMarkerTGraphErr(meanErrorsSummed, 20, 1.,1,1);
	meanErrorsSummed->Draw("p,csame");
	legendMean->AddEntry(meanErrorsSummed,"quad. sum. incl. BR dalitz","p");
	legendMean->Draw();
	canvasSysErrMean->Update();
	canvasSysErrMean->SaveAs(Form("SystematicErrorsNew/SysMean_%s_%s_%s_%s.pdf",meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data()));
	
	delete canvasSysErrMean;

	if (energy.CompareTo("7TeV")==0 && meson.CompareTo("Pi0") == 0 ) {
		TCanvas* canvasSysErrMeanSingle = new TCanvas("canvasSysErrMeanSingle","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasSysErrMeanSingle, 0.08, 0.01, 0.015, 0.09);
		

		for(Int_t i = 0; i< numberCutStudies ; i++){

				histo2DSysErrMean->SetMaximum(20);
				histo2DSysErrMean->Draw();
				meanErrors[i]->Draw("p,csame");
				canvasSysErrMeanSingle->Update();
				canvasSysErrMeanSingle->SaveAs(Form("SystematicErrorsNew/SysMean_%sCutVariation_%s_%s_%s_%s.pdf",nameCutForOutputVariation7TeV[i].Data(),meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data()));
			
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
		histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "histo2DNewSysErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,50.);
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
		legendMeanNew= new TLegend(0.1,0.6,0.38,0.95);
	}

	legendMeanNew->SetTextSize(0.025);
	legendMeanNew->SetFillColor(0);
	legendMeanNew->SetBorderSize(0);

	if (numberCutStudies> 9) legendMeanNew->SetNColumns(2);

	for(Int_t i = 0; i< numberCutStudies ; i++){

			DrawGammaSetMarkerTGraphErr(meanErrorsCorr[i], 20+i, 1.,color[i],color[i]);
			meanErrorsCorr[i]->Draw("p,csame");
			legendMeanNew->AddEntry(meanErrorsCorr[i],nameCutVariation[i].Data(),"p");
			
	}
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummed, 20, 1.,1,1);
	meanErrorsCorrSummed->Draw("p,csame");
	legendMeanNew->AddEntry(meanErrorsCorrSummed,"quad. sum. inc. BR dalitz","p");
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kRed,kRed);
	meanErrorsCorrSummedIncMat->Draw("p,csame");
	legendMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"quad. sum., inc. BR and mat.","p");
	legendMeanNew->Draw();
	
	if (energy.CompareTo("PbPb_2.76TeV")==0){
		TLatex *labelCentrality = new TLatex(0.65,0.93,Form("%s Pb-Pb #sqrt{s_{NN}} = 2.76 TeV"	,additionalName.Data()	));
		SetStyleTLatex( labelCentrality, 0.038,4);
		labelCentrality->Draw();
	} else if ( energy.CompareTo("pPb_5.023TeV") == 0 ) {
		TLatex *labelCentrality = new TLatex(0.65,0.93,Form("%s p-Pb #sqrt{s_{NN}} = 5.02 TeV"	,additionalName.Data()	));
		SetStyleTLatex( labelCentrality, 0.038,4);
		labelCentrality->Draw();
	}
	else {
		if (meson.Contains("Pi0")){
			TLatex *labelCollisionSystem = new TLatex(0.65,0.93,Form("#pi^{0} #rightarrow #gamma e^{+}e^{-}, pp, #sqrt{s} = %s",energy.Data()));
			SetStyleTLatex( labelCollisionSystem, 0.038,4);
			labelCollisionSystem->Draw();
		} else {
			TLatex *labelCollisionSystem = new TLatex(0.65,0.93,Form("#eta #rightarrow #gamma e^{+}e^{-}, pp, #sqrt{s} = %s",energy.Data()));
			SetStyleTLatex( labelCollisionSystem, 0.038,4);
			labelCollisionSystem->Draw();
		}
	}
	
	canvasNewSysErrMean->Update();
	canvasNewSysErrMean->SaveAs(Form("SystematicErrorsNew/SysMeanNewWithMean_%s_%s_%s_%s.pdf",meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data()));
	
	delete canvasNewSysErrMean;


	if (energy.CompareTo("7TeV")==0 && meson.CompareTo("Pi0") == 0 ) {

		TCanvas* canvasSysErrMeanSingle = new TCanvas("canvasSysErrMeanSingle","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasSysErrMeanSingle, 0.08, 0.01, 0.015, 0.09);
		

		for(Int_t i = 0; i< numberCutStudies ; i++){
				histo2DSysErrMean->Draw();
				meanErrorsCorr[i]->Draw("p,csame");
				canvasSysErrMeanSingle->Update();
				canvasSysErrMeanSingle->SaveAs(Form("SystematicErrorsNew/SysMeanSmoothed_%sCutVariation_%s_%s_%s_%s.pdf", nameCutForOutputVariation7TeV[i].Data(), meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data()));
				
		}
		delete canvasSysErrMean;

	}

		const char *SysErrDatname = Form("SystematicErrorsNew/SystematicError_%s_%s_%s_%s.dat",meson.Data(),energy.Data(),additionalNameOutput.Data(),dateForOutput.Data());
		fstream SysErrDat;
		cout << SysErrDatname << endl;
		SysErrDat.open(SysErrDatname, ios::out);
		for (Int_t l=0; l< nPtBins; l++){
			SysErrDat <<errorsNegCorrMatSummed[l] << "\t" <<errorsPosCorrMatSummed[l] << "\t"  <<errorsNegCorrSummed[l] << "\t" <<errorsPosCorrSummed[l]  << endl;
		}

		SysErrDat.close();

		const char *SysErrDatnameMean = Form("SystematicErrorsNew/SystematicErrorAveraged_%s_%s_%s_%s.dat",meson.Data(),energy.Data(),additionalNameOutput.Data(),dateForOutput.Data());
		fstream SysErrDatAver;
		cout << SysErrDatnameMean << endl;
		SysErrDatAver.open(SysErrDatnameMean, ios::out);
		for (Int_t l=0; l< nPtBins; l++){
			SysErrDatAver << "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
		}

		SysErrDatAver.close();
		
		


                Double_t errorsMeanCorrPID[nPtBins];    
                //Double_t errorsMeanCorrSignalExtraction[nPtBins];       
                //Double_t errorsMeanCorrTrackReco[nPtBins];      
                Double_t errorsMeanCorrPhotonReco[nPtBins];

                Double_t errorsMeanCorrSignalExtraction[nPtBins];
                Double_t errorsMeanCorrElectronSelection[nPtBins];
                Double_t errorsMeanCorrTrackReco[nPtBins];
                Double_t errorsMeanCorrPhotonSelection[nPtBins];
                Double_t errorsMeanCorrRecEfficiency[nPtBins];
                Double_t errorsMeanCorrRejPi0GGChannel[nPtBins];
		Double_t errorsMeanCorrBackground[nPtBins];
		
		/////////////Define error clasiffication/////////////////////
		
		Double_t errorsMeanCorrTypeA[nPtBins];
                Double_t errorsMeanCorrTypeB[nPtBins];
                Double_t errorsMeanCorrTypeC[nPtBins];
		Double_t errorsMeanCorrTypeCNoMat[nPtBins];
               
               
	
		
		for (Int_t l=0; l< nPtBins; l++){
		  
			if (energy.CompareTo("PbPb_2.76TeV")==0){
// 				"YieldExtraction-0","BG-1","dEdxE-2","dEdxPi-3", "Cluster-4", "TOF-5", "SinglePt-6", "Chi2-7", "Qt-8", "Alpha-9", "V0Offline-10", "R-11", "PsiPair-12", "Eta-13"
				if (numberCutStudies>9){
					errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]+ errorsMeanCorr[1][l]*errorsMeanCorr[1][l]+errorsMeanCorr[9][l]*errorsMeanCorr[9][l]);
				} else{
					errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]+ errorsMeanCorr[1][l]*errorsMeanCorr[1][l]);
				}
				errorsMeanCorrPID[l] =TMath::Sqrt(errorsMeanCorr[2][l]*errorsMeanCorr[2][l]+ errorsMeanCorr[3][l]*errorsMeanCorr[3][l]+ errorsMeanCorr[5][l]*errorsMeanCorr[5][l]);
				if (numberCutStudies>10){
					errorsMeanCorrPhotonReco[l] =TMath::Sqrt( errorsMeanCorr[7][l]* errorsMeanCorr[7][l]+errorsMeanCorr[8][l]*errorsMeanCorr[8][l]+ errorsMeanCorr[10][l]*errorsMeanCorr[10][l]);
				} else {
					errorsMeanCorrPhotonReco[l] =TMath::Sqrt( errorsMeanCorr[7][l]* errorsMeanCorr[7][l]+errorsMeanCorr[8][l]*errorsMeanCorr[8][l]);
				}
				errorsMeanCorrTrackReco[l] = TMath::Sqrt(errorsMeanCorr[4][l]*errorsMeanCorr[4][l]+errorsMeanCorr[6][l]*errorsMeanCorr[6][l]);

			} else {
                                // {"YieldExtraction","dEdxEG-1","dEdxPiG-2","dEdxEE-3","dEdxPiE-4","SinglePtG-5","TPCclsG-6","DCAxy-7","SinglePtE-8","TPCclsE-9","Chi2-10","Qt-11","MassCut-12","Eta-13","Alpha-14","PsiPairITS-15","Background-16","dEdxPiE-17"};
                                 errorsMeanCorrSignalExtraction[l]  = TMath::Sqrt(  errorsMeanCorr[0][l]  * errorsMeanCorr[0][l] );
                                if ( numberCutStudies > 15 ) {
                                     errorsMeanCorrElectronSelection[l] = TMath::Sqrt(  errorsMeanCorr[1][l]  * errorsMeanCorr[1][l]  + errorsMeanCorr[2][l] * errorsMeanCorr[2][l]  + errorsMeanCorr[3][l] * errorsMeanCorr[3][l] +  errorsMeanCorr[15][l] * errorsMeanCorr[15][l] );
                                } else{
                                     errorsMeanCorrElectronSelection[l] = TMath::Sqrt(  errorsMeanCorr[1][l]  * errorsMeanCorr[1][l]  + errorsMeanCorr[2][l] * errorsMeanCorr[2][l]  + errorsMeanCorr[3][l] * errorsMeanCorr[3][l] );
                                }
                                                            
                                errorsMeanCorrTrackReco[l]         = TMath::Sqrt(  errorsMeanCorr[4][l]  * errorsMeanCorr[4][l]  + errorsMeanCorr[5][l] * errorsMeanCorr[5][l]  + errorsMeanCorr[6][l] * errorsMeanCorr[6][l] + errorsMeanCorr[7][l]*errorsMeanCorr[7][l] + errorsMeanCorr[8][l]*errorsMeanCorr[8][l] );
                                errorsMeanCorrPhotonSelection[l]   = TMath::Sqrt(  errorsMeanCorr[9][l]  * errorsMeanCorr[9][l]  + errorsMeanCorr[10][l]* errorsMeanCorr[10][l] + errorsMeanCorr[11][l]* errorsMeanCorr[11][l] );
                                errorsMeanCorrRecEfficiency[l]     = TMath::Sqrt(  errorsMeanCorr[12][l] * errorsMeanCorr[12][l] + errorsMeanCorr[16][l] * errorsMeanCorr[16][l] );
                                errorsMeanCorrRejPi0GGChannel[l]   = TMath::Sqrt(  errorsMeanCorr[13][l] * errorsMeanCorr[13][l] );
				errorsMeanCorrBackground[l]        = TMath::Sqrt(  errorsMeanCorr[14][l] * errorsMeanCorr[14][l] );
				Double_t sumTemp = pow( errorsMeanCorrSignalExtraction[l],2) + pow( errorsMeanCorrElectronSelection[l],2) + pow( errorsMeanCorrTrackReco[l],2) + pow( errorsMeanCorrPhotonSelection[l],2) + pow( errorsMeanCorrRecEfficiency[l],2) + pow( errorsMeanCorrRejPi0GGChannel[l],2 ) + pow( errorsMeanCorrBackground[l],2);
				Double_t sumTempMatBR = pow( ( sumTemp + errorMaterial*errorMaterial + errorBRDalitz*errorBRDalitz ),0.5);
				//cout<<l<<"  sumTemp: "<<sumTemp<< " sqrt "<<pow(sumTemp,0.5)<<" BRMat "<<sumTempMatBR<<endl;
				
				errorsMeanCorrTypeA[l] =  TMath::Sqrt( pow ( errorsMeanCorr[kYieldExtraction][l], 2) + pow( errorsMeanCorr[kdEdxEE][l],2) + pow( errorsMeanCorr[kdEdxEG][l],2) +  pow( errorsMeanCorr[kSinglePtE][l],2 ) + pow( errorsMeanCorr[kSinglePtG][l],2) + pow( errorsMeanCorr[kTPCClsE][l],2) + pow( errorsMeanCorr[kTPCClsG][l],2) + pow(errorsMeanCorr[kDCAxy][l],2)+ pow(errorsMeanCorr[kBackground][l],2)); 
				errorsMeanCorrTypeB[l] =  TMath::Sqrt ( pow ( errorsMeanCorr[kdEdxPiE][l], 2 ) + pow( errorsMeanCorr[kdEdxPiG][l],2 ) + pow( errorsMeanCorr[kChi2][l],2) + pow( errorsMeanCorr[kQt][l],2) + pow( errorsMeanCorr[kMassCut][l],2 ) + pow(  errorsMeanCorr[kEta][l],2 ) + pow( errorsMeanCorr[kAlpha][l], 2 ) + pow( errorsMeanCorr[kPsiPairITS][l],2 ) );
				errorsMeanCorrTypeC[l] =  TMath::Sqrt(  pow (  errorMaterial,2 ) + pow( errorBRDalitz,2) );
				errorsMeanCorrTypeCNoMat[l] = errorBRDalitz;
				
				//cout<<"sumTemp: "<<pow( sumTemp, 0.5 )<<"   SumTempAC: "<< TMath::Sqrt( errorsMeanCorrTypeA[l] + errorsMeanCorrTypeC[l] )<<endl;
				
			}
		}


        TGraphErrors* meanErrorsSignalExtraction     =  new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSignalExtraction,  ptBinsErr ,errorsMeanErrCorrSummed );
        TGraphErrors* meanErrorsElectronSelection    =  new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrElectronSelection, ptBinsErr ,errorsMeanErrCorrSummed );
        TGraphErrors* meanErrorsTrackReco            =  new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrTrackReco,         ptBinsErr ,errorsMeanErrCorrSummed );
        TGraphErrors* meanErrorsPhotonSelection      =  new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrPhotonSelection,   ptBinsErr ,errorsMeanErrCorrSummed );
        TGraphErrors* meanErorsRecEfficiency         =  new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrRecEfficiency,     ptBinsErr ,errorsMeanErrCorrSummed );
        TGraphErrors* meanErrorsRejPi0GGChannel      =  new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrRejPi0GGChannel,   ptBinsErr ,errorsMeanErrCorrSummed );
        TGraphErrors* meanErrorsBackground           =  new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrBackground,        ptBinsErr ,errorsMeanErrCorrSummed );
	
	
	TGraphErrors* meanErrorsTypeA    	=  new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrTypeA,   	ptBinsErr ,errorsMeanErrCorrSummed );
        TGraphErrors* meanErrorsTypeB           =  new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrTypeB,   	ptBinsErr ,errorsMeanErrCorrSummed );
	TGraphErrors* meanErrorsTypeC		=  new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrTypeC,   	ptBinsErr ,errorsMeanErrCorrSummed );
	TGraphErrors* meanErrorsTypeCNoMat	=  new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrTypeCNoMat,   ptBinsErr ,errorsMeanErrCorrSummed );
	
	




        TCanvas* canvasSummedErrMean = new TCanvas("canvasSummedErrMean","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasSummedErrMean, 0.08, 0.01, 0.015, 0.09);
	TH2D *histo2DSummedErrMean ;
	
	
	if (energy.CompareTo("PbPb_2.76TeV")==0 && meson.CompareTo("Eta") == 0 ) {
		histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "histo2DSummedErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,100.);
	} else	if (energy.CompareTo("PbPb_2.76TeV")==0 ) {
		histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "histo2DSummedErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,50.);
	} else 	if (energy.CompareTo("2.76TeV")==0 && meson.CompareTo("Eta") == 0  ){
		histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "histo2DSummedErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,40.);
	} else  if(  energy.CompareTo("pPb_5.023TeV")==0 ) {
		histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "histo2DSummedErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,50.);	
	}else {
		histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "histo2DSummedErrMean", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,50.);
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
		legendSummedMeanNew= new TLegend(0.1,0.7,0.38,0.95);

	}
	legendSummedMeanNew->SetTextSize(0.030);
	legendSummedMeanNew->SetFillColor(0);
	legendSummedMeanNew->SetBorderSize(0);
	
	DrawGammaSetMarkerTGraphErr(meanErrorsSignalExtraction, 20, 1.,color[0],color[0]);
	meanErrorsSignalExtraction->Draw("p,csame");

	DrawGammaSetMarkerTGraphErr(meanErrorsElectronSelection, 21, 1.,color[1],color[1]);
	meanErrorsElectronSelection->Draw("p,csame");
	DrawGammaSetMarkerTGraphErr(meanErrorsTrackReco, 22, 1.,color[2],color[2]);
	meanErrorsTrackReco->Draw("p,csame");
	DrawGammaSetMarkerTGraphErr(meanErrorsPhotonSelection, 23, 1.,color[3],color[3]);
	meanErrorsPhotonSelection->Draw("p,csame");
        
        DrawGammaSetMarkerTGraphErr(meanErorsRecEfficiency, 24, 1.,color[4],color[4]);
        meanErorsRecEfficiency->Draw("p,csame");

        DrawGammaSetMarkerTGraphErr(meanErrorsRejPi0GGChannel, 25, 1.,color[5],color[5]);
        meanErrorsRejPi0GGChannel->Draw("p,csame");
	
	DrawGammaSetMarkerTGraphErr(meanErrorsBackground, 26, 1.,color[6],color[6]);
        meanErrorsBackground->Draw("p,csame");
	
	//meanErrorsBackground

	DrawGammaSetMarkerTGraphErr(graphMaterialError, 27, 1.,color[7],color[7]);
	graphMaterialError->Draw("p,csame");
	

	legendSummedMeanNew->AddEntry(meanErrorsSignalExtraction,"Signal Extraction","p");
        legendSummedMeanNew->AddEntry(meanErrorsElectronSelection,"Electron PID","p");
	legendSummedMeanNew->AddEntry(meanErrorsTrackReco,"Track Reconstruction","p");
	legendSummedMeanNew->AddEntry(meanErrorsPhotonSelection,"Photon Selection","p");
	legendSummedMeanNew->AddEntry(meanErorsRecEfficiency,"Reconstruction efficiency","p");
        legendSummedMeanNew->AddEntry(meanErrorsRejPi0GGChannel,"#pi^{0} #rightarrow #gamma#gamma rejection","p");
	legendSummedMeanNew->AddEntry(meanErrorsBackground,"Background","p");
        legendSummedMeanNew->AddEntry(graphMaterialError,"Material","p");
	
	
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kBlack,kBlack);
	meanErrorsCorrSummedIncMat->Draw("p,csame");
	legendSummedMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"quad. sum. + BR dalitz","p");
	legendSummedMeanNew->Draw();
	
	if (energy.CompareTo("PbPb_2.76TeV")==0){

		TLatex *labelCentrality = new TLatex(0.65,0.93,Form("%s Pb-Pb #sqrt{s_{NN}} = 2.76 TeV"	,additionalName.Data()	));
		SetStyleTLatex( labelCentrality, 0.038,4);
		labelCentrality->Draw();

	} else if ( energy.CompareTo("pPb_5.023TeV")==0) {
		TLatex *labelCentrality = new TLatex(0.65,0.93,Form("%s p-Pb #sqrt{s_{NN}} = 5.02 TeV"	,additionalName.Data()	));
		SetStyleTLatex( labelCentrality, 0.038,4);
		labelCentrality->Draw();
	}else {

		if (meson.Contains("Pi0")){

			TLatex *labelCollisionSystem = new TLatex(0.55,0.93,Form("#pi^{0} #rightarrow #gamma e^{+}e^{-}, pp, #sqrt{s} = %s",energy.Data()));
			SetStyleTLatex( labelCollisionSystem, 0.038,4);
			labelCollisionSystem->Draw();

		} else {

			TLatex *labelCollisionSystem = new TLatex(0.55,0.93,Form("#eta #rightarrow #gamma e^{+}e^{-}, pp, #sqrt{s} = %s",energy.Data()));
			SetStyleTLatex( labelCollisionSystem, 0.038,4);
			labelCollisionSystem->Draw();

		}
	}
	
	canvasSummedErrMean->Update();
	canvasSummedErrMean->SaveAs(Form("SystematicErrorsNew/SysErrorSummedVisu_%s_%s_%s_%s.pdf",meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data()));
	
	delete canvasSummedErrMean;


		
		
		const char *SysErrDatnameMeanPaper = Form("SystematicErrorsNew/SystematicErrorAveragedPaper_%s_%s_%s_%s.dat",meson.Data(),energy.Data(),additionalNameOutput.Data(),dateForOutput.Data());
		fstream SysErrDatAverPaper;
		cout << SysErrDatnameMeanPaper << endl;
		SysErrDatAverPaper.open(SysErrDatnameMeanPaper, ios::out);
		SysErrDatAverPaper  << "p_{T}" <<"\t DalitzBR "<<"\t Material \t Yield Extraction \t PID \t track recon \t photon recon \t RecEfficiency \t Pi0GGRejection \t Background \t summed" <<  endl;

		for (Int_t l=0; l< nPtBins; l++){

			SysErrDatAverPaper << ptBins[l]<<"\t"<<errorBRDalitz<<"\t"<< errorMaterial <<"\t"<< errorsMeanCorrSignalExtraction[l]<<"\t"<< errorsMeanCorrElectronSelection[l]<<"\t"<< errorsMeanCorrTrackReco[l]<<"\t"<<errorsMeanCorrPhotonSelection[l] <<"\t"<< errorsMeanCorrRecEfficiency[l]<<"\t"<<errorsMeanCorrRejPi0GGChannel[l]<<"\t"<<errorsMeanCorrBackground[l]<<"\t"<< errorsMeanCorrMatSummed[l]<< endl;

		}

 		SysErrDatAverPaper.close();
		
		
		
		
	TCanvas* canvasErrMeanClassified = new TCanvas("canvasErrMeanClassified","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasErrMeanClassified, 0.08, 0.01, 0.015, 0.09);
	TH2D *histo2DErrMeanClassified ;
	
	if (energy.CompareTo("PbPb_2.76TeV")==0 && meson.CompareTo("Eta") == 0 ) {
		histo2DErrMeanClassified = new TH2D("histo2DErrMeanClassified", "histo2DErrMeanClassified", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,100.);
	} else	if (energy.CompareTo("PbPb_2.76TeV")==0 ) {
		histo2DErrMeanClassified = new TH2D("histo2DErrMeanClassified", "histo2DErrMeanClassified", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,50.);
	} else 	if (energy.CompareTo("2.76TeV")==0 && meson.CompareTo("Eta") == 0  ){
		histo2DErrMeanClassified = new TH2D("histo2DErrMeanClassified", "histo2DErrMeanClassified", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,40.);
	} else  if(  energy.CompareTo("pPb_5.023TeV")==0 ) {
		histo2DErrMeanClassified = new TH2D("histo2DErrMeanClassified", "histo2DErrMeanClassified", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,50.);	
	}else {
		histo2DErrMeanClassified = new TH2D("histo2DErrMeanClassified", "histo2DErrMeanClassified", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,50.);
	}
	
	histo2DErrMeanClassified->SetYTitle("mean smoothed systematic Err %");
	histo2DErrMeanClassified->SetXTitle("p_{T} (GeV/c)");
	histo2DErrMeanClassified->GetYaxis()->SetNdivisions(510,kTRUE);
	histo2DErrMeanClassified->GetYaxis()->SetLabelSize(0.03);
	histo2DErrMeanClassified->GetYaxis()->SetTitleSize(0.04);
	histo2DErrMeanClassified->GetYaxis()->SetDecimals();
	histo2DErrMeanClassified->GetYaxis()->SetTitleOffset(0.9);
	histo2DErrMeanClassified->GetXaxis()->SetTitleOffset(1.);
	histo2DErrMeanClassified->GetXaxis()->SetLabelSize(0.03);
	histo2DErrMeanClassified->GetXaxis()->SetTitleSize(0.04);
	histo2DErrMeanClassified->SetTitle("");
	histo2DErrMeanClassified->Draw();
	
	
	TLegend* legendErrMeanClassified;

	if (energy.CompareTo("PbPb_2.76TeV")==0 && additionalName.CompareTo("60-80%")==0 ){
		legendErrMeanClassified= new TLegend(0.13,0.7,0.55,0.95);
	} else {
		legendErrMeanClassified= new TLegend(0.1,0.7,0.38,0.95);

	}
	legendErrMeanClassified->SetTextSize(0.030);
	legendErrMeanClassified->SetFillColor(0);
	legendErrMeanClassified->SetBorderSize(0);
	
	DrawGammaSetMarkerTGraphErr(meanErrorsTypeA, 20, 1.,color[0],color[0]);
	meanErrorsTypeA->Draw("p,csame");
	DrawGammaSetMarkerTGraphErr(meanErrorsTypeB, 21, 1.,color[1],color[1]);
	meanErrorsTypeB->Draw("p,csame");
	DrawGammaSetMarkerTGraphErr(meanErrorsTypeC, 22, 1.,color[2],color[2]);
	meanErrorsTypeC->Draw("p,csame");
	
	DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.5,kBlack,kBlack,1.4,kTRUE);
	meanErrorsCorrSummedIncMat->Draw("p,csame");
	
	//DrawGammaSetMarkerTGraphErr(meanErrorsTypeCNoMat, 23, 1.,color[3],color[3]);
	//meanErrorsTypeCNoMat->Draw("p,csame");
	
	
	legendErrMeanClassified->AddEntry(meanErrorsTypeA,"Errors type A","p");
        legendErrMeanClassified->AddEntry(meanErrorsTypeB,"Errors type B","p");
	legendErrMeanClassified->AddEntry(meanErrorsTypeC,"Errors type C","p");
	//legendErrMeanClassified->AddEntry(meanErrorsTypeCNoMat,"Errors type C No Material","p");
	legendErrMeanClassified->AddEntry(meanErrorsCorrSummedIncMat,"quad. sum ","p");
	
	
	
	legendErrMeanClassified->Draw();
	
	canvasErrMeanClassified->Update();
	canvasErrMeanClassified->SaveAs(Form("SystematicErrorsNew/SysErrorAverageClassified_%s_%s_%s_%s.pdf",meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data()));
	
	delete canvasErrMeanClassified;
	
	
	
	const char *SysErrMeanClassified = Form("SystematicErrorsNew/SystematicErrorAveragedClassified_%s_%s_%s_%s.dat",meson.Data(),energy.Data(),additionalNameOutput.Data(),dateForOutput.Data());
	fstream SysErrMeanClassifiedDat;
	cout << SysErrMeanClassified << endl;
	SysErrMeanClassifiedDat.open(SysErrMeanClassified, ios::out);
	

	for (Int_t l=0; l< nPtBins; l++){

			SysErrMeanClassifiedDat << errorsMeanCorrTypeA[l]<<"\t"<<errorsMeanCorrTypeB[l]<<"\t"<<errorsMeanCorrTypeC[l]<<"\t"<<errorsMeanCorrTypeCNoMat[l]<<endl;

	}

 	SysErrMeanClassifiedDat.close();
	
	//Print type Errors
	
	cout<<"pT"<<"\t"<<nameCutVariationpPb_5023GeV[kYieldExtraction]<<"\t"<< nameCutVariationpPb_5023GeV[kdEdxEE]<<"\t"<<nameCutVariationpPb_5023GeV[kdEdxEG]<<"\t"<<nameCutVariationpPb_5023GeV[kSinglePtE]<<"\t"<<nameCutVariationpPb_5023GeV[kSinglePtG]<<"\t"<<nameCutVariationpPb_5023GeV[kTPCClsE]<<"\t"<<nameCutVariationpPb_5023GeV[kTPCClsG]<<"\t"<<nameCutVariationpPb_5023GeV[kDCAxy]<<"\t"<<nameCutVariationpPb_5023GeV[kBackground]<<endl;
	
	for( Int_t l=0; l < nPtBins; l++){
	  
	 cout<<ptBins[l]<<"\t"<<errorsMeanCorr[kYieldExtraction][l]<<"\t"<< errorsMeanCorr[kdEdxEE][l]<<"\t"<<errorsMeanCorr[kdEdxEG][l]<<"\t"<<errorsMeanCorr[kSinglePtE][l]<<"\t"<<errorsMeanCorr[kSinglePtG][l]<<"\t"<<errorsMeanCorr[kTPCClsE][l]<<"\t"<<errorsMeanCorr[kTPCClsG][l]<<"\t"<<errorsMeanCorr[kDCAxy][l]<<"\t"<<errorsMeanCorr[kBackground][l]<<endl;
	 
	}
	
	cout<<"pT"<<"\t"<<nameCutVariationpPb_5023GeV[kdEdxPiE]<<"\t"<<nameCutVariationpPb_5023GeV[kdEdxPiG]<<"\t"<<nameCutVariationpPb_5023GeV[kChi2]<<"\t"<<nameCutVariationpPb_5023GeV[kQt]<<"\t"<<nameCutVariationpPb_5023GeV[kMassCut]<<"\t"<<nameCutVariationpPb_5023GeV[kEta]<<"\t"<<nameCutVariationpPb_5023GeV[kAlpha]<<"\t"<<nameCutVariationpPb_5023GeV[kPsiPairITS]<<endl;
	
	
	for( Int_t l = 0;  l < nPtBins; l++){
	  
	cout<<ptBins[l]<<"\t"<<errorsMeanCorr[kdEdxPiE][l]<<"\t"<<errorsMeanCorr[kdEdxPiG][l]<<"\t"<<errorsMeanCorr[kChi2][l]<<"\t"<<errorsMeanCorr[kQt][l]<<"\t"<<errorsMeanCorr[kMassCut][l]<<"\t"<<errorsMeanCorr[kEta][l]<<"\t"<<errorsMeanCorr[kAlpha][l]<<"\t"<<errorsMeanCorr[kPsiPairITS][l]<<endl;
	
	}
				
	
	
	//errorsMeanCorrTypeB[l] =  TMath::Sqrt ( pow ( errorsMeanCorr[kdEdxPiE][l], 2 ) + pow( errorsMeanCorr[kdEdxPiG][l],2 ) + pow( errorsMeanCorr[kChi2][l],2) + pow( errorsMeanCorr[kQt][l],2) + pow( errorsMeanCorr[kMassCut][l],2 ) + pow(  errorsMeanCorr[kEta][l],2 ) + pow( errorsMeanCorr[kAlpha][l], 2 ) + pow( errorsMeanCorr[kPsiPairITS][l],2 ) );
				//errorsMeanCorrTypeC[l] =  TMath::Sqrt(  pow (  errorMaterial,2 ) + pow( errorBRDalitz,2) );
				
	////cout<<
	
	
	
		

	

		
		
}
