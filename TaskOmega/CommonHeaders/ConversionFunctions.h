#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TObjString.h"
// #include "../CommonHeaders/Chi2Fit.h"

struct AliConvDataObject {
	Double_t valueX;
	Double_t errorXLow;
	Double_t errorXHigh;
	Double_t valueY;
	Double_t errorYStatLow;
	Double_t errorYStatHigh;
	Double_t errorYSystLow;
	Double_t errorYSystHigh;
	Double_t errorYTotLow;
	Double_t errorYTotHigh;
};

//********************* global Variables *********************************************************************
Double_t ncoll0005 = 1684.4;
Double_t ncoll0510 = 1316;
Double_t ncoll0010 = 1500;
Double_t ncoll1020 = 921.2;
Double_t ncoll0020 = 1210.6;
Double_t ncoll2040 = 438.4;
Double_t ncoll4060 = 127.7;
Double_t ncoll4050 = 171.25;
Double_t ncoll5060 = 84.28;
Double_t ncoll6080 = 26.71;
Double_t ncoll6070 = 37.855;
Double_t ncoll7080 = 15.575;
Double_t ncoll8090 = 6.293;
Double_t ncoll7590 = 8.219;

Double_t nCollErr0005 = 190;
Double_t nCollErr0510 = 140;
Double_t nCollErr0010 = 165;
Double_t nCollErr1020 = 96;
Double_t nCollErr0020 = 130.5;
Double_t nCollErr2040 = 42.;
Double_t nCollErr4060 = 11;
Double_t nCollErr4050 = 16;
Double_t nCollErr5060 = 6.95;
Double_t nCollErr6080 = 2;
Double_t nCollErr6070 = 2.85;
Double_t nCollErr7080 = 1.035;
Double_t nCollErr8090 = 0.325;
Double_t nCollErr7590 = 0.473;

Double_t tAA0005 = 26.32;
Double_t tAA0510 = 20.56;
Double_t tAA0010 = 23.44;
Double_t tAA1020 = 14.39;
Double_t tAA0020 = 18.915;
Double_t tAA2040 = 6.85;
Double_t tAA4060 = 1.996;
Double_t tAA6080 = 0.4174;

Double_t tAAErr0005 = 0.84224;
Double_t tAAErr0510 = 0.65792;
Double_t tAAErr0010 = 0.75008;
Double_t tAAErr1020 = 0.44609;
Double_t tAAErr0020 = 0.5958225;
Double_t tAAErr2040 = 0.22605;
Double_t tAAErr4060 = 0.097804;
Double_t tAAErr6080 = 0.026296;



//*********************** declaration of functions defined in this header ***********************************
Float_t CalculateMeanPt(const TF1* );
void CalculateMeanPtWithError(const TF1* , Float_t& , Float_t& );
TH1D* CalculateHistoRatioToFit (TH1D*, TF1*);
TH1F* CalculateHistoRatioToFit (TH1F*, TF1*);
TH1D* CorrectHistoToBinCenter (TH1D*);
TGraphErrors* CalculateGraphErrRatioToFit (TGraphErrors* , TF1* );
TGraphAsymmErrors* CalculateGraphErrRatioToFit (TGraphAsymmErrors* , TF1* );
TGraph* CalculateGraphRatioToFit (TGraph* , TF1* );
TH1D* CalculateHistoRatioToFitNLO (TH1D* , TF1* , Double_t );
TGraphAsymmErrors* CalculateSysErrFromRelSysHisto( TH1D* , TString , Double_t* , Double_t* , Int_t , const Int_t  );
TGraphAsymmErrors* CalculateSysErrAFromRelSysHisto( TH1D* , TString , Double_t* , Double_t* , Int_t , const Int_t  );
TGraphAsymmErrors* CalculateSysErrFromRelSysHistoComplete( TH1D* , TString , Double_t* , Double_t* , Int_t , const Int_t  );
TGraphAsymmErrors* CalculateCombinedSysAndStatError( TGraphAsymmErrors* , TGraphAsymmErrors* );
TH1D* CalculateShiftedSpectrumInXtOrMt(TH1D* , TH1D* , TString );
TH1D* CalculateMassMinusExpectedMass(TH1D* , Double_t );
void RemoveScalingWithPt(TH1D* );
void RemoveScalingWithPtGraph(TGraphAsymmErrors* );
void ScaleWithPtGraph(TGraphAsymmErrors* );
TGraph* ScaleGraph (TGraph* , Double_t );
TGraphAsymmErrors* ScaleGraph (TGraphAsymmErrors* , Double_t );
TGraphErrors* ScaleGraph (TGraphErrors* , Double_t );
Double_t* ExtractRelErrDownAsymmGraph(TGraphAsymmErrors* );
Double_t* ExtractRelErrUpAsymmGraph(TGraphAsymmErrors* );
// TGraphAsymmErrors* CalculateGraphAsymErrRatioToGraphErr (TGraphAsymmErrors* , TGraphAsymmErrors* );
TGraphAsymmErrors* RebinCombPi0Graph(TGraphAsymmErrors* , TH1D* );
TGraphErrors* RebinNLOGraph(TGraphErrors* );
// TH1D *GraphToHist(TGraphErrors *graph,Int_t maxPt = 50, TString name = "");
// TH1D* RebinTH1D(TH1D* , TH1D* , Bool_t );
TH1D* CalculateFitToFitRatio(TF1 *, TF1 *);
TGraphAsymmErrors* CombinePtPointsEta(TH1D* ,          TGraphAsymmErrors* ,
                                      TH1D* ,             TGraphAsymmErrors* ,
//                                       TGraphAsymmErrors* ,    TGraphAsymmErrors* ,  
                                      Int_t  );
TGraphErrors* CombinePtPointsEtaCalcStat(TGraphAsymmErrors* , TH1D* , TH1D* , Int_t   );
TH1D* CalculateWeightedAveragePCM( TH1D* , TH1D* );
// TGraphAsymmErrors* CombinePtPointsSpectra( TH1D* ,    TGraphAsymmErrors* , TH1D* ,     TGraphAsymmErrors* , TGraphAsymmErrors* &,  TGraphAsymmErrors* &,  Double_t* , Int_t , Int_t ,        Int_t ,  Int_t , Bool_t , Int_t );
// TGraphAsymmErrors* CombinePtPointsRAA( TGraphAsymmErrors* ,      TGraphAsymmErrors* , TGraphAsymmErrors* ,      TGraphAsymmErrors* , TGraphAsymmErrors* &,  TGraphAsymmErrors* &,  Double_t* , Int_t , Int_t ,        Int_t ,  Int_t , Bool_t , Int_t );
// void CalculateFitResults(TF1* , TF1* , Double_t*  ,TString , Double_t );
void ReadOutFromSysErrVector(Double_t* , Double_t* , Double_t* , Int_t , Int_t ,Int_t , Int_t ,Double_t );
void CalculateMeanSysErr(Double_t* , Double_t* , Double_t* , Double_t* , Int_t );   
void CorrectSystematicErrorsWithMean(Double_t* ,Double_t* , Double_t* , Double_t* , Int_t );
void ProduceGraphAsymmWithoutXErrors(TGraphAsymmErrors* );
void ProduceGraphPartialXErrors(TGraphErrors* , Double_t );
void ProduceGraphAsymmPartialXErrors(TGraphAsymmErrors* , Double_t );
void ProduceGraphFixedXErrors(TGraphErrors* , Double_t );
void ProduceGraphAsymmFixedXErrors(TGraphAsymmErrors* , Double_t );
void ProduceGraphAsymmWithoutXYErrors(TGraphAsymmErrors* );
void ProduceGraphErrWithoutXErrors(TGraphErrors* );
void ProduceGraphErrDisplacedX(TGraphErrors* ,Double_t );
TGraphAsymmErrors* Add2TGraphAsymmErrorsSameBinning(TGraphAsymmErrors* ,TGraphAsymmErrors* );
TGraphAsymmErrors* CorrectTGraphAsymmErrorsToBinCenter(TGraphAsymmErrors* );
TH1D* ShortChargedHadronHisto(TH1D* );
TH1D* ConvertChargedHadronHisto(TH1D* , TH1D* );
// TF1 *BinShiftTH1D(TH1D *, TH1D **, TString , TString ,Double_t , Double_t* );
// TF1 *ApplyYShift(TGraphAsymmErrors *, TGraphAsymmErrors **, TString , TString ,Double_t , Double_t* , Double_t , Bool_t );
TGraphAsymmErrors* ApplyYshiftIndividualSpectra(TGraphAsymmErrors * , TF1 *);
Double_t bin_shift_x(TF1 *, Double_t , Double_t , Double_t );
TGraphAsymmErrors *ApplyXshift(TGraphAsymmErrors *, TF1 *);
// TGraphAsymmErrors* ApplyXshiftIndividualSpectra(TGraphAsymmErrors *, TGraphAsymmErrors * , TF1 *, Int_t , Int_t );
TString ReturnDateString();
TString ReturnDateStringForOutput();
Double_t ReturnRapidityStringAndDouble(TString , TString& );
Double_t ReturnDeltaEta(TString);
// void ReturnSeparatedCutNumber(TString , TString& ,TString& ,TString& , Bool_t );
Float_t ReturnBackgroundMult(TString );
Double_t GetNCollFromCutNumber (TString );
Double_t GetScalingFactorSecCorrection(TString );
Double_t GetNCollFromName (TString );
Double_t GetTAAFromName (TString );
Double_t GetNCollErrFromCutNumber (TString );
Double_t GetNCollErrFromName (TString );
Double_t GetTAAErrFromName (TString );
void ReturnParameterSetFittingPbPb(TString , Double_t* );
void ReturnParameterSetFittingPbPbFromString(TString , Double_t* );
TString GetCentralityString(TString );
// void CalcRaa(  TGraphAsymmErrors* , TGraphAsymmErrors* , TGraphAsymmErrors* , TF1* , TGraphAsymmErrors* , TGraphAsymmErrors* , TGraphAsymmErrors** , TGraphAsymmErrors** , Double_t , Double_t , Double_t , Int_t );
Int_t GetBinning(TObject *, Double_t* );
// Int_t CompareBinning(Int_t , Double_t* , Int_t , Double_t* , Double_t* , Int_t* , Int_t* , TString , Double_t );
// Int_t FindFirstCommonBin(Double_t* , Double_t* , Double_t);
// void RebinObjects(TObject* ,TObject* , Double_t* , Int_t* , Int_t , AliConvDataObject* , TString , TString  ,Bool_t );
// TGraphErrors* CalculateRatioBetweenSpectraWithDifferentBinning(TObject* , TObject* , TObject* , TObject* , Bool_t scaleByBinCenterA=kTRUE,  Bool_t scaleByBinCenterB=kTRUE);
Int_t GetNEvents (TH1D* );
Int_t GetNEvents (TH1F* );
TString AnalyseTPCClusterCut(Int_t );
TString AnalyseAlphaMesonCut(Int_t );
TString AnalyseTPCdEdxCutElectronLine(Int_t);
TString AnalyseTPCdEdxCutPionLine(TString);
TString AnalyseChi2GammaCut(Int_t , Int_t);
TString AnalyseChi2MesonCut(Int_t );
TString AnalyseCosPointCut(Int_t );
TString AnalyseChi2PsiPair(Int_t );
TString AnalyseQtMaxCut(Int_t );  
TString AnalyseSinglePtCut(Int_t);
TString AnalyseBackgroundScheme(TString );
TString AnalyseEtaCut(Int_t );
TString AnalyseRCut(Int_t );
TString AnalyseRapidityMesonCut(Int_t );
TString AnalysePsiPair(Int_t, Int_t );
Double_t AnalyseDCAZPhotonCutValue(Int_t dcaZPhoton);
TString AnalyseDCAZPhotonCut(Int_t dcaZPhoton);
  
//*********************** definition of functions defined in this header ***********************************
Float_t CalculateMeanPt(const TF1* fit){

	// fit is the invariant cross section (or proportinal to it):
	// fit ~ 1/pT dN/dpT

	TF1* fs = (TF1*) fit->Clone();
	fs->SetName("fs");

	TF1* f1 = new TF1("f1","x*fs",0.,16.);
	TF1* f2 = new TF1("f2","x^2*fs",0.,16.);

	Float_t i1 = f1->Integral(0.,16.);
	Float_t i2 = f2->Integral(0.,16.);
	

	delete fs;
	delete f1;
	delete f2;

	return i2/i1;
}

void CalculateMeanPtWithError(const TF1* fSigma, Float_t& meanpt, Float_t& meanpt_error){

	meanpt = CalculateMeanPt(fSigma);

	// quadratic summ of errors
	Float_t sum_err2 = 0.;

	// loop over the parameters of the function fSigma
	for (Int_t i=0; i<fSigma->GetNpar(); i++) {
		
		Float_t par = fSigma->GetParameter(i);
		Float_t sigma = fSigma->GetParError(i);

		// mean pt for a function with par(i) = par(i) + sigma(par(i))
		TF1* fmod_up = (TF1*) fSigma->Clone();      
		fmod_up->SetName("fmod_up");
		fmod_up->SetParameter(i, par+sigma);
		Float_t delta_meanpt_up = fabs(meanpt - CalculateMeanPt(fmod_up)); 

		// mean pt for a function with par(i) = par(i) - sigma(par(i))
		TF1* fmod_down = (TF1*) fSigma->Clone();
		fmod_up->SetName("fmod_down");
		fmod_down->SetParameter(i, par-sigma);
		Float_t delta_meanpt_down = fabs(meanpt - CalculateMeanPt(fmod_down)); 

		// take the larger one of the two deltas
		Float_t delta_meanpt = delta_meanpt_up;
		if (delta_meanpt_up < delta_meanpt_down) delta_meanpt = delta_meanpt_down;
		
		sum_err2 += delta_meanpt*delta_meanpt;

		delete fmod_up;
		delete fmod_down;

	}

	meanpt_error = sqrt(sum_err2);
}

TH1D* InvertHisto (TH1D* histo){
	TH1D* histo2 = (TH1D*)histo->Clone("Dummy");
	for( Int_t I = 1; I < histo->GetNbinsX() +1 ;I++){
		Double_t yValue = 1./histo2->GetBinContent(I);
		Double_t yError = 1./histo2->GetBinContent(I)*histo2->GetBinError(I)/histo2->GetBinContent(I);
		histo2->SetBinContent(I,yValue);
		histo2->SetBinError(I,yError);
	}
	return histo2;
}


TH1D* CalculateHistoRatioToFit (TH1D* histo, TF1* fit){
	TH1D* histo2 = (TH1D*)histo->Clone("Dummy");
	for( Int_t I = 1; I < histo->GetNbinsX() +1 ;I++){
		Double_t xValue = histo2->GetBinCenter(I);
		Double_t yValue = fit->Eval(xValue);
		Double_t formerYValue = histo2->GetBinContent(I);
		if (yValue != 0){
			histo2->SetBinContent(I,formerYValue/yValue);
			histo2->SetBinError(I,histo2->GetBinError(I)/yValue);
		}
	}
	return histo2;
}

TH1F* CalculateHistoRatioToFit (TH1F* histo, TF1* fit){
	TH1F* histo2 = (TH1F*)histo->Clone("Dummy");
	for( Int_t I = 1; I < histo->GetNbinsX() +1 ;I++){
		Double_t xValue = histo2->GetBinCenter(I);
		Double_t yValue = fit->Eval(xValue);
		Double_t formerYValue = histo2->GetBinContent(I);
		if (yValue != 0){
			histo2->SetBinContent(I,formerYValue/yValue);
			histo2->SetBinError(I,histo2->GetBinError(I)/yValue);
		}
	}
	return histo2;
}

TH1D* CorrectHistoToBinCenter (TH1D* histo){
	histo->Sumw2();
	for( Int_t I = 1; I < histo->GetNbinsX() +1 ;I++){
		Double_t xValue = histo->GetBinCenter(I);
		Double_t yValue = histo->GetBinContent(I);
		Double_t yValueError = histo->GetBinError(I);
		histo->SetBinContent(I,yValue/xValue);
		histo->SetBinError(I,yValueError/xValue);
	}
	return histo;
}

TGraphErrors* CalculateGraphErrRatioToFit (TGraphErrors* graph_Org, TF1* fit){
	TGraphErrors* graph = (TGraphErrors*)graph_Org->Clone("Dummy");
	Double_t * xValue = graph->GetX(); 
	Double_t * yValue = graph->GetY();
	Double_t* xError = graph->GetEX();
	Double_t* yError = graph->GetEY();
	Int_t nPoints = graph->GetN();
	for (Int_t i = 0; i < nPoints; i++){
		yValue[i] = yValue[i]/fit->Eval(xValue[i]);
		yError[i] = yError[i]/fit->Eval(xValue[i]);
	}
	TGraphErrors* returnGraph =  new TGraphErrors(nPoints,xValue,yValue,xError,yError); 
	return returnGraph;
}

TGraphAsymmErrors* CalculateGraphErrRatioToFit (TGraphAsymmErrors* graph_Org, TF1* fit){
	TGraphAsymmErrors* graph = (TGraphAsymmErrors*)graph_Org->Clone("Dummy");
	Double_t * xValue = graph->GetX(); 
	Double_t * yValue = graph->GetY();
	Double_t* xErrorLow = graph->GetEXlow();
	Double_t* xErrorHigh = graph->GetEXhigh();
	Double_t* yErrorLow = graph->GetEYlow();
	Double_t* yErrorHigh = graph->GetEYhigh();
	Int_t nPoints = graph->GetN();
	for (Int_t i = 0; i < nPoints; i++){
		yValue[i] = yValue[i]/fit->Eval(xValue[i]);
		yErrorLow[i] = yErrorLow[i]/fit->Eval(xValue[i]);
		yErrorHigh[i] = yErrorHigh[i]/fit->Eval(xValue[i]);
	}
	TGraphAsymmErrors* returnGraph =  new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh); 
	return returnGraph;
}

TGraph* CalculateGraphRatioToFit (TGraph* graph, TF1* fit){
	Double_t * xValue = graph->GetX(); 
	Double_t * yValue = graph->GetY();
	Int_t nPoints = graph->GetN();
	for (Int_t i = 0; i < nPoints; i++){
		yValue[i] = yValue[i]/fit->Eval(xValue[i]);
	}
	TGraph* returnGraph =  new TGraph(nPoints,xValue,yValue); 
	return returnGraph;
}

TH1D* CalculateHistoRatioToFitNLO (TH1D* histo, TF1* fit, Double_t startX){
	Double_t beginBin = histo->GetXaxis()->FindBin(startX);
	for (Int_t L = 0; L < beginBin; L++){
		histo->SetBinContent(L+1,-1);
	}
	for (Int_t L = beginBin; L < histo->GetNbinsX()+1; L++){
		Double_t xValue = histo->GetBinCenter(L);
		Double_t yValue = fit->Eval(xValue);
		Double_t formerYValue = histo->GetBinContent(L);
		if (formerYValue > 0){
			histo->SetBinContent(L,formerYValue/yValue);
                        histo->SetBinError(L,0);//histo->GetBinError(L)/yValue);
		}
	}	
	return histo;
}

TGraphAsymmErrors* CalculateSysErrFromRelSysHisto( TH1D* histo, TString nameGraph, Double_t* relSystErrorDown, Double_t* relSystErrorUp, Int_t offset, const Int_t nPoints ){
	Double_t xValueCorr[nPoints];
	Double_t yValueCorr[nPoints];
	//	Double_t yErrorCorr[nPoints];
	Double_t xErrorCorr[nPoints];
	Double_t systErrorDown[nPoints];
	Double_t systErrorUp[nPoints];
	for (Int_t i = 0; i < nPoints; i++){
		xValueCorr[i] = histo->GetBinCenter(i+offset);
		yValueCorr[i] = histo->GetBinContent(i+offset);
		xErrorCorr[i] = histo->GetBinWidth(i+offset)/2.;
		systErrorDown[i] = relSystErrorDown[i]*yValueCorr[i]/100.*(-1);
		systErrorUp[i] = relSystErrorUp[i]*yValueCorr[i]/100.;
		cout << i << "\t" << i+offset <<"\t" <<xValueCorr[i] << "\t" << yValueCorr[i] << "\t" << relSystErrorDown[i] << endl;
	}
	cout << "create graph" << endl;
	TGraphAsymmErrors* graphSys = new TGraphAsymmErrors(nPoints,xValueCorr,yValueCorr,xErrorCorr,xErrorCorr,systErrorDown,systErrorUp);
	graphSys->SetName(nameGraph);
	cout << "graph could be generated" << endl;
	return graphSys;
}

TGraphAsymmErrors* CalculateSysErrAFromRelSysHisto( TH1D* histo, TString nameGraph, Double_t* relSystErrorDown, Double_t* relSystErrorUp, Int_t offset, const Int_t nPoints ){
	Double_t xValueCorr[nPoints];
	Double_t yValueCorr[nPoints];
	//	Double_t yErrorCorr[nPoints];
	Double_t xErrorCorr[nPoints];
	Double_t systErrorADown[nPoints];
	Double_t systErrorAUp[nPoints];
	
	Double_t relSystErrorADown[nPoints];
	Double_t relSystErrorAUp[nPoints];
	Double_t relSysMatDown=2.*6.2;
	Double_t relSysMatUp  =2.*3.4;
	for (Int_t i = 0; i < nPoints; i++){
		if( TMath::Abs(relSystErrorDown[i])>=relSysMatDown &&  TMath::Abs(relSystErrorUp[i])>=relSysMatUp){
			relSystErrorADown[i] = -TMath::Sqrt(relSystErrorDown[i]*relSystErrorDown[i]-relSysMatDown*relSysMatDown);
			relSystErrorAUp[i] = TMath::Sqrt(relSystErrorUp[i]*relSystErrorUp[i]-relSysMatUp*relSysMatUp);
		}else{
			relSystErrorADown[i] = -TMath::Sqrt(relSystErrorDown[i]*relSystErrorDown[i]);
			relSystErrorAUp[i] = TMath::Sqrt(relSystErrorUp[i]*relSystErrorUp[i]);
		}
		xValueCorr[i] = histo->GetBinCenter(i+offset);
		yValueCorr[i] = histo->GetBinContent(i+offset);
		xErrorCorr[i] = histo->GetBinWidth(i+offset)/2.;
		systErrorADown[i] = relSystErrorADown[i]*yValueCorr[i]/100.*(-1);
		systErrorAUp[i] = relSystErrorAUp[i]*yValueCorr[i]/100.;
		cout << xValueCorr[i] << "\t" << yValueCorr[i] << endl;
		
	}
	
	TGraphAsymmErrors* graphSys = new TGraphAsymmErrors(nPoints,xValueCorr,yValueCorr,xErrorCorr,xErrorCorr,systErrorADown,systErrorAUp);
	graphSys->SetName(nameGraph);
	return graphSys;
}

TGraphAsymmErrors* CalculateSysErrFromRelSysHistoComplete( TH1D* histo, TString nameGraph, Double_t* relSystErrorDown, Double_t* relSystErrorUp, Int_t offset, const Int_t nPoints ){
	Double_t xValueCorr[nPoints];
	Double_t yValueCorr[nPoints];
	Double_t yErrorCorr[nPoints];
	Double_t xErrorCorr[nPoints];
	Double_t systErrorDown[nPoints];
	Double_t systErrorUp[nPoints];
	for (Int_t i = 0; i < nPoints; i++){
		xValueCorr[i] = histo->GetBinCenter(i+offset);
		xErrorCorr[i] = histo->GetBinWidth(i+offset)/2.;
		yValueCorr[i] = histo->GetBinContent(i+offset);
		yErrorCorr[i] = histo->GetBinError(i+offset);
		systErrorDown[i] = TMath::Sqrt((TMath::Power(yErrorCorr[i],2) + TMath::Power(relSystErrorDown[i]*yValueCorr[i]/100.,2)));
		systErrorUp[i] = TMath::Sqrt(TMath::Power(yErrorCorr[i],2) + TMath::Power(relSystErrorUp[i]*yValueCorr[i]/100.,2));
		cout << xValueCorr[i] << "\t" << yValueCorr[i] << endl;
	}
	
	TGraphAsymmErrors* graphSys = new TGraphAsymmErrors(nPoints,xValueCorr,yValueCorr,xErrorCorr,xErrorCorr,systErrorDown,systErrorUp);
	
// 	graphSys->Print();
	graphSys->SetName(nameGraph);
	return graphSys;
}

TGraphAsymmErrors* CalculateCombinedSysAndStatError( TGraphAsymmErrors* graphStat, TGraphAsymmErrors* graphSys){
   TGraphAsymmErrors* graphStatCopy = (TGraphAsymmErrors*)graphStat->Clone("graphStatCopy");
   TGraphAsymmErrors* graphSysCopy = (TGraphAsymmErrors*)graphSys->Clone("graphSysCopy");
   Double_t* xValue = graphStatCopy->GetX();
   Double_t* xErrorLow = graphStatCopy->GetEXlow();
   Double_t* xErrorHigh = graphStatCopy->GetEXhigh();
   Double_t* yValueStat = graphStatCopy->GetY();
   Double_t* yErrorLowStat = graphStatCopy->GetEYlow();
   Double_t* yErrorHighStat = graphStatCopy->GetEYhigh();
   Double_t* yErrorLowSys = graphSysCopy->GetEYlow();
   Double_t* yErrorHighSys = graphSysCopy->GetEYhigh();
   Int_t nPoints = graphStatCopy->GetN();
   Double_t yErrorLowComb[nPoints];
   Double_t yErrorHighComb[nPoints];
   for (Int_t i = 0; i < nPoints; i++){
      yErrorLowComb[i] = TMath::Sqrt((TMath::Power(yErrorLowStat[i],2) + TMath::Power(yErrorLowSys[i],2)));
      yErrorHighComb[i] = TMath::Sqrt((TMath::Power(yErrorHighStat[i],2) + TMath::Power(yErrorHighSys[i],2)));
      
   }
   
   TGraphAsymmErrors* graphSysAndComb = new TGraphAsymmErrors(nPoints,xValue,yValueStat,xErrorLow,xErrorHigh,yErrorLowComb,yErrorHighComb);
   
//    graphSysAndComb->Print();
   graphSysAndComb->SetName(graphSys->GetName());
   return graphSysAndComb;
}


TH1D* CalculateShiftedSpectrumInXtOrMt(TH1D* ptSpectrum, TH1D* spectrum, TString nameNewSpectrum){
	Double_t dummyValue = 0;
	Double_t dummyError = 0;
	TH1D* shiftedSpectrum = (TH1D*)spectrum->Clone(nameNewSpectrum);
	for (Int_t i = 0; i < shiftedSpectrum->GetNbinsX()+1; i++){
		dummyValue = ptSpectrum->GetBinContent(i);
		dummyError = ptSpectrum->GetBinError(i);
		shiftedSpectrum->SetBinContent(i,dummyValue);
		shiftedSpectrum->SetBinError(i,dummyError);
	}
	return shiftedSpectrum;
}

TH1D* CalculateMassMinusExpectedMass(TH1D* histo, Double_t mass){
	TH1D* histoNew= (TH1D*) histo->Clone();
	for ( Int_t l=0; l < histoNew->GetNbinsX()+1; l++){
		Double_t intermediateValue = histoNew->GetBinContent(l);
		Double_t intermediateError = histoNew->GetBinError(l);
		if (intermediateValue != 0) {
			histoNew->SetBinContent(l,(intermediateValue - mass)*1000);
			histoNew->SetBinError(l,intermediateError*1000);
		}
	}	
	return histoNew;
}

void RemoveScalingWithPt(TH1D* histo){
	for (Int_t i = 0; i < histo->GetNbinsX()+1;i++){
		histo->SetBinContent(i, histo->GetBinContent(i)*histo->GetBinCenter(i));
		histo->SetBinError(i, histo->GetBinError(i)*histo->GetBinCenter(i));
	}
	return;
}

void RemoveScalingWithPtGraph(TGraphAsymmErrors* graph){
	Double_t * xValue = graph->GetX(); 
	Double_t * yValue = graph->GetY();
	Double_t* yErrorHigh = graph->GetEYhigh();
	Double_t* yErrorLow = graph->GetEYlow();
	Int_t nPoints = graph->GetN();
	for (Int_t i = 0; i < nPoints; i++){
		yValue[i] = yValue[i]*xValue[i];
		yErrorHigh[i] = yErrorHigh[i]*xValue[i];
		yErrorLow[i] = yErrorLow[i]*xValue[i];
	}
	return;
}

void ScaleWithPtGraph(TGraphAsymmErrors* graph){
	Double_t * xValue = graph->GetX(); 
	Double_t * yValue = graph->GetY();
	Double_t* yErrorHigh = graph->GetEYhigh();
	Double_t* yErrorLow = graph->GetEYlow();
	Int_t nPoints = graph->GetN();
	for (Int_t i = 0; i < nPoints; i++){
		yValue[i] = yValue[i]/xValue[i];
		yErrorHigh[i] = yErrorHigh[i]/xValue[i];
		yErrorLow[i] = yErrorLow[i]/xValue[i];
	}
	return;
}

TGraph* ScaleGraph (TGraph* graph, Double_t scaleFac){
   TGraph* dummyGraph = (TGraph*)graph->Clone(Form("%s_Scaled",graph->GetName()));
	Double_t * xValue = dummyGraph->GetX();
	Double_t * yValue = dummyGraph->GetY();
	
	Int_t nPoints = dummyGraph->GetN();
	for (Int_t i = 0; i < nPoints; i++){
		yValue[i] = yValue[i]*scaleFac;
	}
	TGraph* returnGraph = new TGraph(nPoints,xValue,yValue);
	return returnGraph;
}

TGraphAsymmErrors* ScaleGraph (TGraphAsymmErrors* graph, Double_t scaleFac){
   TGraphAsymmErrors* dummyGraph = (TGraphAsymmErrors*)graph->Clone(Form("%s_Scaled",graph->GetName()));
   
	Double_t * xValue = dummyGraph->GetX(); 
	Double_t * yValue = dummyGraph->GetY();
	Double_t* xErrorLow = dummyGraph->GetEXlow();
	Double_t* xErrorHigh = dummyGraph->GetEXhigh();
	Double_t* yErrorLow = dummyGraph->GetEYlow();
	Double_t* yErrorHigh = dummyGraph->GetEYhigh();
	Int_t nPoints = dummyGraph->GetN();
	for (Int_t i = 0; i < nPoints; i++){
		yValue[i] = yValue[i]*scaleFac;
		yErrorLow[i] = yErrorLow[i]*scaleFac;
		yErrorHigh[i] = yErrorHigh[i]*scaleFac;
	}
	TGraphAsymmErrors* returnGraph =  new TGraphAsymmErrors(nPoints,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh); 
	return returnGraph;
}

TGraphErrors* ScaleGraph (TGraphErrors* graph, Double_t scaleFac){
	TGraphErrors* dummyGraph = (TGraphErrors*)graph->Clone(Form("%s_Scaled",graph->GetName()));
	Double_t * xValue = dummyGraph->GetX(); 
	Double_t * yValue = dummyGraph->GetY();
	Double_t* xError = dummyGraph->GetEX();
	Double_t* yError = dummyGraph->GetEY();
	Int_t nPoints = dummyGraph->GetN();
	for (Int_t i = 0; i < nPoints; i++){
		yValue[i] = yValue[i]*scaleFac;
		yError[i] = yError[i]*scaleFac;
	}
	TGraphErrors* returnGraph =  new TGraphErrors(nPoints,xValue,yValue,xError,yError); 

	return returnGraph;
}

Double_t* ExtractRelErrDownAsymmGraph(TGraphAsymmErrors* graph){
	Double_t * yValue = graph->GetY();
	Double_t* yErrorLow = graph->GetEYlow();
	Double_t* yErrorLow2 = yErrorLow;
	Int_t nPoints = graph->GetN();
	for (Int_t i = 0; i < nPoints; i++){
		yErrorLow2[i] = -yErrorLow2[i]/yValue[i]*100;
	}
	return yErrorLow2;
}

Double_t* ExtractRelErrUpAsymmGraph(TGraphAsymmErrors* graph){
	Double_t * yValue = graph->GetY();
	Double_t* yErrorHigh = graph->GetEYhigh();
	Double_t* yErrorHigh2 = yErrorHigh;
	Int_t nPoints = graph->GetN();
	
	for (Int_t i = 0; i < nPoints; i++){
		yErrorHigh2[i] = yErrorHigh2[i]/yValue[i]*100;
	}
	return yErrorHigh2;
}

TGraphAsymmErrors* CalculateGraphAsymErrRatioToGraphErr (TGraphAsymmErrors* graphAsymA_v1, TGraphAsymmErrors* graphAsymB_v1,Bool_t binomial = kFALSE, Double_t commonPercentageADown = 0., Double_t commonPercentageAUp = 0., Double_t commonPercentageBDown = 0., Double_t commonPercentageBUp = 0.){
	TGraphAsymmErrors *graphAsymA = (TGraphAsymmErrors*) graphAsymA_v1->Clone("graphAsymA");

	Double_t * xValueAsymA = graphAsymA->GetX(); 
	Double_t * yValueAsymA = graphAsymA->GetY();
	Double_t* xErrorLowAsymA = graphAsymA->GetEXlow();
	Double_t* xErrorHighAsymA = graphAsymA->GetEXhigh();
	Double_t* yErrorLowAsymA = graphAsymA->GetEYlow();
	Double_t* yErrorHighAsymA = graphAsymA->GetEYhigh();
	
	TGraphAsymmErrors *graphAsymB = (TGraphAsymmErrors*) graphAsymB_v1->Clone("graphAsymB");
	//	Double_t * xValueAsymB = graphAsymB->GetX(); 
	Double_t * yValueAsymB = graphAsymB->GetY();
	Double_t* yErrorLowAsymB = graphAsymB->GetEYlow();
	
	//	Double_t* xErrorHighAsymB = graphAsymB->GetEXhigh();
	//      Double_t* yErrorLowAsymB = graphAsymB->GetEYlow();
	Double_t* yErrorHighAsymB = graphAsymB->GetEYhigh();
	
	Int_t nPointsB = graphAsymB->GetN();
	Int_t nPoints = graphAsymA->GetN();
	if (nPoints > nPointsB) nPoints = nPointsB;
	Double_t yValueRatio[nPoints];
	for (Int_t i = 0; i < nPoints; i++){
		yValueRatio[i] = yValueAsymA[i]/yValueAsymB[i];
	//       cout << "before low: "<<yErrorLowAsymA[i]/yValueAsymA[i] << "\t" << commonPercentageADown << "\t"  <<  yErrorLowAsymB[i]/yValueAsymB[i] << "\t"  << commonPercentageBDown << endl;
	//       cout << "before up: "<<yErrorHighAsymA[i]/yValueAsymA[i] << "\t" << commonPercentageAUp << "\t"  <<  yErrorHighAsymB[i]/yValueAsymB[i] << "\t"  <<  commonPercentageBUp << endl;
		yErrorLowAsymA[i] = TMath::Sqrt(pow(yErrorLowAsymA[i],2) - pow(commonPercentageADown*yValueAsymA[i]/100,2));
		yErrorHighAsymA[i] = TMath::Sqrt(pow(yErrorHighAsymA[i],2) - pow(commonPercentageAUp*yValueAsymA[i]/100,2));
		yErrorLowAsymB[i] = TMath::Sqrt(pow(yErrorLowAsymB[i],2) - pow(commonPercentageBDown*yValueAsymB[i]/100,2));
		yErrorHighAsymB[i] = TMath::Sqrt(pow(yErrorHighAsymB[i],2) - pow(commonPercentageBUp*yValueAsymB[i]/100,2));
	//       cout << "after low: "<< yErrorLowAsymA[i]/yValueAsymA[i] << "\t"  <<  yErrorLowAsymB[i]/yValueAsymB[i] << endl;
	//       cout << "after up: "<< yErrorHighAsymA[i]/yValueAsymA[i] << "\t"  <<  yErrorHighAsymB[i]/yValueAsymB[i] << endl;
		if (binomial){
			Double_t w = yValueAsymA[i]/yValueAsymB[i];
			yErrorLowAsymA[i] = TMath::Abs( ( (1.-2.*w)*yErrorLowAsymA[i]*yErrorLowAsymA[i] + w*w*yErrorLowAsymB[i]*yErrorLowAsymB[i] )/(yValueAsymB[i]*yValueAsymB[i]) );
			yErrorHighAsymA[i] = TMath::Abs( ( (1.-2.*w)*yErrorHighAsymA[i]*yErrorHighAsymA[i] + w*w*yErrorHighAsymB[i]*yErrorHighAsymB[i] )/(yValueAsymB[i]*yValueAsymB[i]) );
		} else {
			
			yErrorLowAsymA[i] =  TMath::Sqrt( pow(yErrorLowAsymA[i]/yValueAsymB[i],2)  + pow( yErrorLowAsymB[i]*yValueAsymA[i]/pow(yValueAsymB[i],2),2) );
			yErrorHighAsymA[i] = TMath::Sqrt( pow(yErrorHighAsymA[i]/yValueAsymB[i],2) + pow( yErrorHighAsymB[i]*yValueAsymA[i]/pow(yValueAsymB[i],2),2) );
		}
	}
	TGraphAsymmErrors* returnGraph =  new TGraphAsymmErrors(nPoints,xValueAsymA,yValueRatio,xErrorLowAsymA,xErrorHighAsymA,yErrorLowAsymA,yErrorHighAsymA); 
	return returnGraph;
}

TGraphAsymmErrors* RebinCombPi0Graph(TGraphAsymmErrors* graphAsymV1, TH1D* histoGraph){
	TGraphAsymmErrors *graphAsym = (TGraphAsymmErrors*) graphAsymV1->Clone("graphAsym");

	Double_t * xValueAsym = graphAsym->GetX(); 
	Double_t * yValueAsym = graphAsym->GetY();
	//Double_t* xErrorLowAsym = graphAsym->GetEXlow();
	//Double_t* xErrorHighAsym = graphAsym->GetEXhigh();
	Double_t* yErrorLowAsym = graphAsym->GetEYlow();
	Double_t* yErrorHighAsym = graphAsym->GetEYhigh();
	Int_t nPoints = graphAsym->GetN();
	
	Int_t nBins = histoGraph->GetNbinsX();
	Double_t *newBinningY = new Double_t[nBins];
	Double_t *newBinningX = new Double_t[nBins];
	Double_t *newBinningErrorLowY = new Double_t[nBins];
	Double_t *newBinningErrorLowX = new Double_t[nBins];
	Double_t *newBinningErrorHighY = new Double_t[nBins];
	Double_t *newBinningErrorHighX = new Double_t[nBins];
	
	Int_t notusedBins = 0;

	for(Int_t bin = 1; bin <nBins+1; bin++){
		Double_t yValue = 0;
		Double_t yErrorLow = 0;
		Double_t yErrorHigh = 0;
		Int_t nMergedBins = 0;
		Double_t binCenter = histoGraph->GetBinCenter(bin);
		Double_t binWidth = histoGraph->GetBinWidth(bin);
	
		for(Int_t graphBin = 0; graphBin<nPoints; graphBin++){
		
		if(histoGraph->FindBin(xValueAsym[graphBin]) == bin){
			yValue = yValue+yValueAsym[graphBin];
			yErrorLow = yErrorLow+yErrorLowAsym[graphBin];
			yErrorHigh = yErrorHigh+yErrorHighAsym[graphBin];
			nMergedBins++;
		}
		
		}
		if(nMergedBins == 0){
		notusedBins++;
		continue;
		}
		
		newBinningY[bin-1-notusedBins] = yValue/nMergedBins;
		newBinningX[bin-1-notusedBins] = binCenter;
		newBinningErrorLowX[bin-1-notusedBins] = binWidth/2;
		newBinningErrorHighX[bin-1-notusedBins] = binWidth/2;
		newBinningErrorLowY[bin-1-notusedBins] = yErrorLow/nMergedBins;
		newBinningErrorHighY[bin-1-notusedBins] = yErrorHigh/nMergedBins;
		
		histoGraph->SetBinContent(bin,newBinningY[bin-1-notusedBins]);
		if(newBinningErrorLowY[bin-1-notusedBins]>=newBinningErrorHighY[bin-1-notusedBins]) histoGraph->SetBinError(bin,newBinningErrorLowY[bin-1-notusedBins]);
		else histoGraph->SetBinError(bin,newBinningErrorHighY[bin-1-notusedBins]);
	}

	Double_t *newnewBinningY = new Double_t[nBins-notusedBins];
	Double_t *newnewBinningX = new Double_t[nBins-notusedBins];
	Double_t *newnewBinningErrorLowY = new Double_t[nBins-notusedBins];
	Double_t *newnewBinningErrorLowX = new Double_t[nBins-notusedBins];
	Double_t *newnewBinningErrorHighY = new Double_t[nBins-notusedBins];
	Double_t *newnewBinningErrorHighX = new Double_t[nBins-notusedBins];

	
	for(Int_t graphBin = 0; graphBin<nBins-notusedBins; graphBin++){
		newnewBinningY[graphBin] = newBinningY[graphBin];
		newnewBinningX[graphBin] = newBinningX[graphBin];
		newnewBinningErrorLowY[graphBin] = newBinningErrorLowY[graphBin];
		newnewBinningErrorLowX[graphBin] = newBinningErrorLowX[graphBin];
		newnewBinningErrorHighY[graphBin] = newBinningErrorHighY[graphBin];
		newnewBinningErrorHighX[graphBin] = newBinningErrorHighX[graphBin];
	}

	TGraphAsymmErrors* returnGraph =  new TGraphAsymmErrors(nBins-notusedBins,newnewBinningX,newnewBinningY,newnewBinningErrorLowX,newnewBinningErrorHighX,newnewBinningErrorLowY,newnewBinningErrorHighY);

	delete graphAsym;
	delete[] newBinningY;
	delete[] newBinningX;
	delete[] newBinningErrorLowY;
	delete[] newBinningErrorLowX;
	delete[] newBinningErrorHighY;
	delete[] newBinningErrorHighX;
	delete[] newnewBinningY;
	delete[] newnewBinningX;
	delete[] newnewBinningErrorLowY;
	delete[] newnewBinningErrorLowX;
	delete[] newnewBinningErrorHighY;
	delete[] newnewBinningErrorHighX;
	return returnGraph;
}

TGraphErrors* RebinNLOGraph(TGraphErrors* graphV1){
	TGraphErrors *graph = (TGraphErrors*) graphV1->Clone("graph");

	Double_t * xValue = graph->GetX(); 
	Double_t * yValue = graph->GetY();
	Int_t nPoints = graph->GetN();
	
	Int_t newBinning = nPoints*0.5;
	Double_t *newBinningY = new Double_t[newBinning];
	Double_t *newBinningX = new Double_t[newBinning];
	Double_t *newBinningErrorY = new Double_t[newBinning];
	Double_t *newBinningErrorX = new Double_t[newBinning];
	Int_t newBin = 0;

	for(Int_t i = 2; i<nPoints; i=i+2){
		Double_t yValueOne = 0;
		Double_t yValueTwo = 0;
		Double_t yValueThree = 0;
		Double_t yValueNew = 0;
		Double_t xValueOne = 0;
		Double_t xValueTwo = 0;
		Double_t xValueThree = 0;
		Double_t xValueNew = 0;

		
		yValueOne = yValue[i-1];
		yValueTwo = yValue[i];
		yValueThree = yValue[i+1];

		xValueOne = xValue[i-1];
		xValueTwo = xValue[i];
		xValueThree = xValue[i+1];
		

		yValueNew = (yValueOne + yValueTwo+yValueThree)/3;
		xValueNew = (xValueOne + xValueTwo+xValueThree)/3;

		//   cout<<i<<"    "<<xValueOne<<"  "<<xValueTwo<<"   "<<xValueThree<<"   "<<  xValueNew<<endl;
		newBinningY[newBin] = yValueNew;
		newBinningX[newBin] = xValueNew;
		newBinningErrorY[newBin] = 1.;
		newBinningErrorX[newBin] = xValueTwo-xValueOne;
		//cout<<-(xValueOne-xValueTwo)*0.5<<endl;

		newBin++;
	}
	
	
	TGraphErrors* returnGraph =  new TGraphErrors(newBinning,newBinningX,newBinningY,newBinningErrorX,newBinningErrorY);

	delete graph;
	return returnGraph;
}

TH1D *GraphToHist(TGraphErrors *graph,Int_t maxPt = 50, TString name = ""){
	Double_t *xValue =  graph->GetX(); 
	Double_t *yValue = graph->GetY();
	Double_t * Ex = graph->GetEX();
	Int_t  nPoints = graph->GetN();
	Int_t maxPoints = 0;

	for(Int_t i = 0; i<nPoints; i++){
		if(xValue[i]<=maxPt) maxPoints++;
	}

	Double_t *newBinningX = new Double_t[maxPoints];
	for(Int_t i = 0;i<maxPoints;i++) newBinningX[i] = xValue[i]-Ex[i];

	TH1D *hist = new TH1D(name,"",maxPoints-1,newBinningX);

	for(Int_t i = 1;i<maxPoints;i++) hist->SetBinContent(i,yValue[i-1]);

	return hist;
}

TH1D* RebinTH1D(TH1D* histoV2, TH1D* histoBinningV2, Bool_t deltaPt = kFALSE){

	Double_t binArray[60];

	TH1D *histo = (TH1D*) histoV2->Clone(histoV2->GetName());
	TH1D *histoBinning = (TH1D*) histoBinningV2->Clone(histoBinningV2->GetName());

	Int_t nBinsPi0 = histoBinning->GetNbinsX();
	Bool_t sameStartBin = kFALSE;
	Int_t missingStartBins = 0;

	for(Int_t bin = 0; bin<nBinsPi0; bin++){

		if(!sameStartBin && (histoBinning->GetBinLowEdge(bin+1) < histo->GetBinLowEdge(1))){
			missingStartBins++;
			continue;
		}
		else if(!sameStartBin && histoBinning->GetBinLowEdge(bin+1) >= histo->GetBinLowEdge(1)){
			sameStartBin = kTRUE;
		}

		binArray[bin-missingStartBins] = histoBinning->GetBinLowEdge(bin+1);
		if(bin == (nBinsPi0-1)) binArray[bin+1-missingStartBins] = histoBinning->GetXaxis()->GetBinUpEdge(bin+1);
	}
	
	nBinsPi0 = nBinsPi0 - missingStartBins;
	
	TH1D *Rebin = new TH1D(Form("Rebin_%s",histoV2->GetName()),"",nBinsPi0,binArray);
	Rebin->Sumw2();
	TH1D *DeltaPt = new TH1D(Form("DeltaPt_%s",histoV2->GetName()),"",nBinsPi0,binArray);
	DeltaPt->Sumw2();


	for(Int_t iPt = 0; iPt < nBinsPi0; iPt++){
		Int_t startBin = histo->GetXaxis()->FindBin(binArray[iPt]+0.001);
		Int_t endBin = histo->GetXaxis()->FindBin(binArray[iPt+1]-0.001);
		Int_t nMergedBins = 1+endBin - startBin;
		Rebin->SetBinContent(iPt+1,nMergedBins);
		Rebin->SetBinError(iPt+1,0);
		Double_t diffPt = binArray[iPt+1]-binArray[iPt];
		DeltaPt->SetBinContent(iPt+1,diffPt);
		DeltaPt->SetBinError(iPt+1,0);
	}


	
	histo = (TH1D*) histo->Rebin(nBinsPi0,histo->GetName(),binArray);
	if(!deltaPt) histo->Divide(Rebin);
	if(deltaPt) histo->Divide(DeltaPt);

	delete Rebin;
	delete DeltaPt;
	delete histoBinning;

	return histo;
}

TH1D* CalculateFitToFitRatio(TF1 *funcA, TF1 *funcB){

	Double_t maxA = funcA->GetXmax();
	Double_t maxB = funcB->GetXmax();
	Double_t newMax = 0;
	if(maxA>maxB) newMax = maxB;
	else newMax = maxA;
	
	Int_t nBinsX = newMax/0.1;
	
	TH1D *returnHist = new TH1D(Form("%s/%s",funcA->GetName(),funcB->GetName()),"",nBinsX,0,newMax);

	for(Int_t i = 1; i<nBinsX+1; i++){
		Double_t binCenter = returnHist->GetBinCenter(i);
		Double_t newBinContent = funcA->Eval(binCenter)/funcB->Eval(binCenter);
		returnHist->SetBinContent(i,newBinContent);
		returnHist->SetBinError(i,0);
	}
	return returnHist;
}

TGraphAsymmErrors* CombinePtPointsEta(TH1D* histoConv,          TGraphAsymmErrors* graphSystConv,
												  TH1D* histoPHOS,             TGraphAsymmErrors* graphSystPHOS,
// 												  TGraphAsymmErrors* graphStatComb, 	TGraphAsymmErrors* graphSystComb,  
												  Int_t  bin0PHOS){
	
	cout << "combining eta points"<< endl;  
	TGraphErrors*   graphStatErrPHOS = new TGraphErrors(histoPHOS);
	TGraphErrors*   graphStatErrConv = new TGraphErrors(histoConv);
	
	Int_t nPtLimitsEta=14;
	Double_t xPtLimitsEta[14]       =  {0.4,0.7,1.0,1.4,1.8,2.2,2.6,3.0,3.5,4.0,6.0,8.0,10.,15.};
	
	Int_t nPtBinsCombEta=13;  
	Double_t xCombEta[13];
	Double_t xECombEta[13];
	Double_t xSectionCombEta[13];
	Double_t xSectionCombErrEta[13];
	Double_t xSectionCombErrLEta[13];
	Double_t xSectionCombErrHEta[13];
	
	Int_t nPHOSEta = graphSystPHOS->GetN();
	Double_t * xPHOSEta =  graphSystPHOS->GetX();
	Double_t * yPHOSEta =  graphSystPHOS->GetY();
	Double_t * eySysPHOSEta =  graphSystPHOS->GetEYlow();
	Double_t * exSysPHOSEta =  graphSystPHOS->GetEXlow();
	Double_t * eyStaPHOSEta =  graphStatErrPHOS->GetEY();
	Double_t * eTot2PHOSEta;
	eTot2PHOSEta = new Double_t[nPHOSEta];
	Double_t * eTotPHOSEta;
	eTotPHOSEta = new Double_t [nPHOSEta];
	
	for(Int_t i=0;i<nPHOSEta;i++){
		eTot2PHOSEta[i]=(eyStaPHOSEta[i]*eyStaPHOSEta[i]+eySysPHOSEta[i]*eySysPHOSEta[i]);
		eTotPHOSEta[i]=TMath::Sqrt( eTot2PHOSEta[i]);
		cout<< "PHOS::"<< xPHOSEta[i]<< " "<<yPHOSEta[i]<< " " <<  eTotPHOSEta[i]<< endl;
	}
	
	Int_t  offset=1;
	Int_t nPCMEta;
	Double_t * xPCMEta;
	Double_t * yPCMEta;
	Double_t * exSysPCMEta;
	Double_t * eySysLPCMEta;
	Double_t * eySysHPCMEta;
	
	nPCMEta= graphSystConv->GetN();
	xPCMEta = graphSystConv->GetX();
	yPCMEta= graphSystConv->GetY();
	exSysPCMEta = graphSystConv->GetEXlow();
	eySysLPCMEta = graphSystConv->GetEYlow();
	eySysHPCMEta = graphSystConv->GetEYhigh();
	
	Double_t * eyStaPCMEta = graphStatErrConv->GetEY();
	
	Double_t eTotL2PCMEta[nPCMEta];
	Double_t eTotH2PCMEta[nPCMEta];
	Double_t eTotLPCMEta[nPCMEta];
	Double_t eTotHPCMEta[nPCMEta];
	for(Int_t i=0;i<nPCMEta;i++){
		eTotH2PCMEta[i]=(eyStaPCMEta[i+offset]*eyStaPCMEta[i+offset]+eySysHPCMEta[i]*eySysHPCMEta[i]);
		eTotHPCMEta[i]=TMath::Sqrt( eTotH2PCMEta[i]);
		eTotL2PCMEta[i]=(eyStaPCMEta[i+offset]*eyStaPCMEta[i+offset]+eySysLPCMEta[i]*eySysLPCMEta[i]);
		eTotLPCMEta[i]=TMath::Sqrt( eTotL2PCMEta[i]);
		cout<< "PCM::"<< xPCMEta[i]<< " "<<yPCMEta[i]<< " " <<  eTotLPCMEta[i]<< " "<< eyStaPCMEta[i+offset]<<" "<< eySysHPCMEta[i]<< endl;
	}
	cout<<endl;  
	
	Bool_t okPHOS,okPCM;
	Int_t  bin0PCM=0;
	
	for (Int_t i=0;i<nPtLimitsEta-1;i++){
		Double_t xCenterEta = 0.5*(xPtLimitsEta[i+1]+xPtLimitsEta[i]);
		cout<< "xCenterEta::"<< i<< " " <<xCenterEta<< endl;
		okPHOS=kFALSE;
		okPCM=kFALSE;
						
		xCombEta[i]=xCenterEta;
	
	
		if((i-bin0PHOS)>=0){
			if ( xPHOSEta[i-bin0PHOS] == xCenterEta){
				okPHOS=kTRUE;
			}
		}
	
		if( (i-bin0PCM) >= 0){
			if ( xPCMEta[i-bin0PCM] == xCenterEta){
				okPCM=kTRUE;
			}
		}
	
	
		if ( okPHOS && okPCM ){
			xECombEta[i]=exSysPCMEta[i-bin0PCM];
			if( eTotL2PCMEta[i-bin0PCM]!=0. &&  eTotH2PCMEta[i-bin0PCM]!=0. && eTot2PHOSEta[i-bin0PHOS] !=0.){
				if(yPCMEta[i-bin0PCM]> yPHOSEta[i-bin0PHOS]){
					xSectionCombEta[i] = (yPCMEta[i-bin0PCM]/eTotL2PCMEta[i-bin0PCM] +  yPHOSEta[i-bin0PHOS]/eTot2PHOSEta[i-bin0PHOS])/
					(1.0/ eTotL2PCMEta[i-bin0PCM] +  1.0/eTot2PHOSEta[i-bin0PHOS]);
					xSectionCombErrEta[i] = pow((1.0/ eTotL2PCMEta[i-bin0PCM] +  1.0/eTot2PHOSEta[i-bin0PHOS]),-0.5);
				
					cout<< " PHOS,PCM_L::"<< xSectionCombEta[i]<< " " << xSectionCombErrEta[i] << " "  
					<< yPCMEta[i-bin0PCM]<< " "<< eTotLPCMEta[i-bin0PCM] << " "
					<< yPHOSEta[i-bin0PHOS]<<" "<< eTotPHOSEta[i-bin0PHOS]<< endl;
				}else{
					xSectionCombEta[i] = (yPCMEta[i-bin0PCM]/eTotH2PCMEta[i-bin0PCM] +  yPHOSEta[i-bin0PHOS]/eTot2PHOSEta[i-bin0PHOS])/
					(1.0/ eTotH2PCMEta[i-bin0PCM] +  1.0/eTot2PHOSEta[i-bin0PHOS]);
					xSectionCombErrEta[i] = pow((1.0/ eTotH2PCMEta[i-bin0PCM] +  1.0/eTot2PHOSEta[i-bin0PHOS]),-0.5);
				
					cout<< " PHOS,PCM_H::"<< xSectionCombEta[i]<< " " << xSectionCombErrEta[i] << " "
					<< yPCMEta[i-bin0PCM]<< " "<< eTotHPCMEta[i-bin0PCM] << " "
					<< yPHOSEta[i-bin0PHOS]<<" "<< eTotPHOSEta[i-bin0PHOS] <<endl;
				}
			}
			xSectionCombErrLEta[i]= xSectionCombErrEta[i];
			xSectionCombErrHEta[i]= xSectionCombErrEta[i];
		}
	
		if ( okPHOS && !okPCM ){
			xECombEta[i]=exSysPHOSEta[i-bin0PHOS];
			if( eTot2PHOSEta[i-bin0PHOS] !=0.){
				xSectionCombEta[i] = yPHOSEta[i-bin0PHOS];
				xSectionCombErrEta[i] = pow((eTot2PHOSEta[i-bin0PHOS]),0.5);
				cout<< " PHOS, OK::"<< xSectionCombEta[i]<< " " << xSectionCombErrEta[i] << " "  
				<< yPHOSEta[i-bin0PHOS]<<" "<< eTotPHOSEta[i-bin0PHOS] <<endl;
			}
			xSectionCombErrLEta[i]= xSectionCombErrEta[i];
			xSectionCombErrHEta[i]= xSectionCombErrEta[i];
		}
	
		if ( !okPHOS && okPCM ){
			xECombEta[i]=exSysPCMEta[i-bin0PCM];
			if( eTotL2PCMEta[i-bin0PCM] !=0. && eTotH2PCMEta[i-bin0PCM]){
				xSectionCombEta[i] = yPCMEta[i-bin0PCM];
				xSectionCombErrLEta[i] = pow((eTotL2PCMEta[i-bin0PCM]),0.5);      // Asymmetric errors needed
				xSectionCombErrHEta[i] = pow((eTotH2PCMEta[i-bin0PCM]),0.5);      // Asymmetric errors needed
				cout<< " PCM_OK::"<< xSectionCombEta[i]<< " " << xSectionCombErrEta[i] << " "  
				<< yPCMEta[i-bin0PCM]<< " "<< eTotLPCMEta[i-bin0PCM] << " "<< endl;
			}
		}
	} //  end of pt point
	
	TGraphAsymmErrors* returnGraph = new TGraphAsymmErrors(nPtBinsCombEta,xCombEta,xSectionCombEta,xECombEta,xECombEta,xSectionCombErrLEta,xSectionCombErrHEta);
// 	returnGraph->Print();
	return returnGraph;
}

TGraphErrors* CombinePtPointsEtaCalcStat(TGraphAsymmErrors* graphSystErr, TH1D* histoConv, TH1D* histoPHOS, Int_t  bin0PHOS ){
       
	cout << "combining eta points"<< endl;  
	TGraphErrors*   graphStatErrPHOS = new TGraphErrors(histoPHOS);
	TGraphErrors*   graphStatErrConv = new TGraphErrors(histoConv);
	graphStatErrConv->RemovePoint(0);
	
	Int_t nPtLimitsEta=14;
	Double_t xPtLimitsEta[14]       =  {0.4,0.7,1.0,1.4,1.8,2.2,2.6,3.0,3.5,4.0,6.0,8.0,10.,15.};
	
	Int_t nPtBinsCombEta=14;  
//        Double_t xCombEta[14];
	Double_t xSectionCombErrEtaStat[14];
	
	Double_t * xPHOSEta =  graphStatErrPHOS->GetX();
	Double_t * eyStaPHOSEta =  graphStatErrPHOS->GetEY();
	
	Double_t * xValueSystError = graphSystErr->GetX();
	Double_t * xErrorSystError = graphSystErr->GetEXlow();
	Double_t * yValueSystError = graphSystErr->GetY();
	
	Double_t * xPCMEta;
	
	xPCMEta = graphStatErrConv->GetX();    
	Double_t * eyStaPCMEta = graphStatErrConv->GetEY();
	
	Bool_t okPHOS,okPCM;
	Int_t  bin0PCM=0;
	
	for (Int_t i=0;i<nPtLimitsEta-1;i++){
		Double_t xCenterEta = 0.5*(xPtLimitsEta[i+1]+xPtLimitsEta[i]);
		cout<< "xCenterEta::"<< i<< " " <<xCenterEta<< endl;
		okPHOS=kFALSE;
		okPCM=kFALSE;
					
// xCombEta[i]=xCenterEta;
	
		if((i-bin0PHOS)>=0){
			if ( xPHOSEta[i-bin0PHOS] == xCenterEta){
				okPHOS=kTRUE;
			}
		}
	
		if( (i-bin0PCM) >= 0){
			if ( xPCMEta[i-bin0PCM] == xCenterEta){
				okPCM=kTRUE;
			}
		}
	
		if ( okPHOS && okPCM ){
			xSectionCombErrEtaStat[i]= pow((1.0/ TMath::Power(eyStaPCMEta[i-bin0PCM],2.) +  1.0/TMath::Power(eyStaPHOSEta[i-bin0PHOS],2.)),-0.5);
		}
		if ( okPHOS && !okPCM ){
			xSectionCombErrEtaStat[i]=eyStaPHOSEta[i-bin0PHOS];
		}
		if ( !okPHOS && okPCM ){
			xSectionCombErrEtaStat[i]=eyStaPCMEta[i-bin0PCM];
		}
	} //  end of pt point
	
	TGraphErrors* returnGraph = new TGraphErrors(nPtBinsCombEta,xValueSystError,yValueSystError,xErrorSystError,xSectionCombErrEtaStat);
// 	returnGraph->Print();
	return returnGraph;
}

TH1D* CalculateWeightedAveragePCM( TH1D* histoConvOnfly,	TH1D* histoConvOffline){
   
	TH1D* returnHisto = (TH1D*)histoConvOnfly->Clone("returnHisto");
	TGraphErrors* graphStatErrOffline = new TGraphErrors(histoConvOffline);  
	TGraphErrors* graphStatErrOnfly = new TGraphErrors(histoConvOnfly);
	
	Int_t nOffline = graphStatErrOffline->GetN();
	Double_t * xOffline =  graphStatErrOffline->GetX();
	Double_t * yOffline =  graphStatErrOffline->GetY();
	Double_t * eyStaOffline =  graphStatErrOffline->GetEY();
	
	for(Int_t i=0;i<nOffline;i++){
		cout<< "Offline:: "<< i<< "\t" <<xOffline[i]<< " "<<yOffline[i]<< " " <<  eyStaOffline[i]<< endl;
	}
	
	Int_t nOnfly = graphStatErrOnfly->GetN();
	Double_t * xOnfly =  graphStatErrOnfly->GetX();
	Double_t * yOnfly =  graphStatErrOnfly->GetY();
	Double_t * eyStaOnfly =  graphStatErrOnfly->GetEY();
	
	for(Int_t i=0;i<nOnfly;i++){
		cout<< "Onfly::"<< i<< "\t" << xOnfly[i]<< " "<<yOnfly[i]<< " " <<  eyStaOnfly[i]<< endl;
	}
	
	cout<<endl;  
	
	Bool_t okOffline,okOnfly;     
	for (Int_t i=0;i<nOnfly;i++){
		okOffline=kFALSE;
		okOnfly=kFALSE;
						
		if (  yOffline[i]!= 0.){
			okOffline=kTRUE;
		}
	
		if (  yOnfly[i]!= 0.){
			okOnfly=kTRUE;
		}
	
		if ( okOffline && okOnfly ){
			if (eyStaOffline[i]== 0. && eyStaOnfly[i]>0. ){
				returnHisto->SetBinContent(i+1, 	yOnfly[i]);
				returnHisto->SetBinError(i+1, eyStaOnfly[i]);
				cout<< " Combined::"<< i << "\t empty Offline " << yOnfly[i] << "\t" << eyStaOnfly[i]<<  endl;
			} else if( eyStaOnfly[i]!=0. &&  eyStaOffline[i]!=0. && eyStaOnfly[i]>0. && eyStaOffline[i]>0.){
				Double_t wOffline = 1./eyStaOffline[i];
				Double_t wOnfly = 1./eyStaOnfly[i];
				Double_t wSum = wOnfly+wOffline;
				Double_t xSectionComb = (wOnfly*yOnfly[i] +  wOffline*yOffline[i])/ wSum;
				Double_t xSectionCombErr= pow(1./2.* wOnfly/wSum* eyStaOnfly[i]*eyStaOnfly[i] + 1./2.* wOffline/wSum* eyStaOffline[i]*eyStaOffline[i],0.5);
				cout<< " Combined::"<< i << "\t" <<xSectionComb<< " " << xSectionCombErr << " " << yOnfly[i]<< " "<< eyStaOnfly[i] << " "
				<< yOffline[i]<<" "<< eyStaOffline[i]<< endl;
				returnHisto->SetBinContent(i+1, 	xSectionComb);
				returnHisto->SetBinError(i+1, xSectionCombErr);
			} else {
				returnHisto->SetBinContent(i+1, 0.);
				returnHisto->SetBinError(i+1, 0.);
				cout<< " Combined::"<< i << "\t empty" << endl;
				cout << "here" << endl;
			}
		} else {
			returnHisto->SetBinContent(i+1, 0.);
			returnHisto->SetBinError(i+1, 0.);
			cout<< " Combined::"<< i << "\t empty" << endl;
		}
	}
	return returnHisto;
}

TGraphAsymmErrors* CombinePtPointsSpectra( TH1D* histoConv,		TGraphAsymmErrors* graphSystConv,
													TH1D* histoPHOS,		TGraphAsymmErrors* graphSystPHOS,
 													TGraphAsymmErrors* &graphStatComb, 	TGraphAsymmErrors* &graphSystComb,  
													Double_t* xPtLimits,	Int_t nPtLimits,
													Int_t offset,			Int_t bin0PCM,	 Int_t bin0PHOS, Bool_t kRemoveLastPCMPoint=kFALSE, Int_t nBinsPCMRem = 1){
   
	TGraphErrors* graphStatErrPHOS = new TGraphErrors(histoPHOS);  
	TGraphErrors* graphStatErrConv = new TGraphErrors(histoConv);
	TGraphAsymmErrors* graphSystConvClone = (TGraphAsymmErrors*)graphSystConv->Clone("DummyConv");  
	TGraphAsymmErrors* graphSystPHOSClone = (TGraphAsymmErrors*)graphSystPHOS->Clone("DummyConv");  
	
	Double_t xComb[50];
	Double_t xEComb[50];
	Double_t xSectionComb[50];
	Double_t xSectionCombErr[50];
	Double_t xSectionCombErrL[50];
	Double_t xSectionCombErrH[50];
	Double_t xSectionCombStatErr[50];
	Double_t xSectionCombSysErr[50];
	
	cout << "********************************************************************************" << endl;
	cout << "************************** PHOS ************************************************" << endl;
	cout << "********************************************************************************" << endl;
	graphStatErrPHOS->Print();
	
	cout << "********************************************************************************" << endl;
	cout << "************************** PCM ************************************************" << endl;
	cout << "********************************************************************************" << endl;
	graphStatErrConv->Print();
   
	Int_t nPHOS = graphSystPHOSClone->GetN();
	Double_t * xPHOS =  graphStatErrPHOS->GetX();
	Double_t * yPHOS =  graphStatErrPHOS->GetY();
	Double_t * eySysPHOS =  graphSystPHOSClone->GetEYlow();
	Double_t * exSysPHOS =  graphSystPHOSClone->GetEXlow();
	Double_t * eyStaPHOS =  graphStatErrPHOS->GetEY();
	Double_t * eTot2PHOS;
	eTot2PHOS = new Double_t[nPHOS];
	Double_t * eTotPHOS;
	eTotPHOS=new Double_t [nPHOS];
	
	for(Int_t i=0;i<nPHOS;i++){
		eTot2PHOS[i]=(eyStaPHOS[i]*eyStaPHOS[i]+eySysPHOS[i]*eySysPHOS[i]);
		eTotPHOS[i]=TMath::Sqrt( eTot2PHOS[i]);
 		cout<< "PHOS::"<< xPHOS[i]<< " "<<yPHOS[i]<< " " <<  eTotPHOS[i]<< endl;
	}
	
	Int_t nPCM;
	Double_t * xPCM;
	Double_t * yPCM;
	Double_t * exSysPCM;
	Double_t * eySysLPCM;
	Double_t * eySysHPCM;
	
	nPCM= graphSystConvClone->GetN();
	if (kRemoveLastPCMPoint) {
		nPCM= nPCM-nBinsPCMRem;
	}
	xPCM = graphSystConvClone->GetX();
	yPCM= graphSystConvClone->GetY();
	exSysPCM = graphSystConvClone->GetEXlow();
	eySysLPCM = graphSystConvClone->GetEYlow();
	eySysHPCM = graphSystConvClone->GetEYhigh();

	Double_t * eyStaPCM = graphStatErrConv->GetEY();
// 	graphSystConvClone->Print();
// 	graphStatErrConv->Print();
	Double_t eTotL2PCM[nPCM];
	Double_t eTotH2PCM[nPCM];
	Double_t eTotLPCM[nPCM];
// 	Double_t eTotHPCM[nPCM];
	for(Int_t i=0;i<nPCM;i++){
// 		arrayPCMSpec[i].valueX
// 		arrayPCMSpec[i].errorXHigh = histoCorrectedYieldCut[j]->GetBinError(i);

		eTotH2PCM[i]=(eyStaPCM[i+offset]*eyStaPCM[i+offset]+eySysHPCM[i]*eySysHPCM[i]);
// 		eTotHPCM[i]=TMath::Sqrt( eTotH2PCM[i]);
		eTotL2PCM[i]=(eyStaPCM[i+offset]*eyStaPCM[i+offset]+eySysLPCM[i]*eySysLPCM[i]);
		eTotLPCM[i]=TMath::Sqrt( eTotL2PCM[i]);
 		cout<< "PCM::"<< xPCM[i]<< " "<<yPCM[i]<< " " <<  eTotLPCM[i]<< " "<< eyStaPCM[i+offset]<<" "<< eySysHPCM[i]<< endl;
	}
	
	
	cout<<endl;  
	
	Bool_t okPHOS,okPCM;     
	for (Int_t i=0;i<nPtLimits-1;i++){
		Double_t xCenter = 0.5*(xPtLimits[i+1]+xPtLimits[i]);
// 		cout<< "xCenter::"<< i<< " " <<xCenter<< endl;
		okPHOS=kFALSE;
		okPCM=kFALSE;
				
		xComb[i]=xCenter;

		xSectionComb[i] = 0;
		xSectionCombErrL[i]= 0;
		xSectionCombErrH[i]= 0;	

		
		if((i-bin0PHOS)>=0){
 			cout << xPHOS[i-bin0PHOS] << "\t" << yPHOS[i-bin0PHOS] << endl;
			if ( xPHOS[i-bin0PHOS] == xCenter && yPHOS[i-bin0PHOS]!= 0.){
				okPHOS=kTRUE;
			}
// 			cout << "PHOS ok " << okPHOS << endl;
		}
	
		if( (i-bin0PCM) >= 0){
			cout << xPCM[i-bin0PCM] <<"\t" <<   yPCM[i-bin0PCM]  << endl;
			if (i-bin0PCM < nPCM){
				if ( xPCM[i-bin0PCM] == xCenter && yPCM[i-bin0PCM] !=0){
					okPCM=kTRUE;
				}
			}
// 			cout << "PCM ok " << okPCM << endl;
		}
	
		if ( okPHOS && okPCM ){
			xEComb[i]=exSysPCM[i-bin0PCM];
			if( eTotL2PCM[i-bin0PCM]!=0. &&  eTotH2PCM[i-bin0PCM]!=0. && eTot2PHOS[i-bin0PHOS] !=0.){
				if(yPCM[i-bin0PCM]> yPHOS[i-bin0PHOS]){
					Double_t wPHOS = 1./eTot2PHOS[i-bin0PHOS];
					Double_t wPCM = 1./eTotL2PCM[i-bin0PCM];
					Double_t wSum = wPCM+wPHOS;
					xSectionComb[i] = (wPCM*yPCM[i-bin0PCM] +  wPHOS*yPHOS[i-bin0PHOS])/ wSum;
					xSectionCombErr[i] = pow((wPCM +  wPHOS),-0.5);    
					xSectionCombStatErr[i] = pow(1./2.* wPCM/wSum* eyStaPCM[i-bin0PCM+offset]*eyStaPCM[i-bin0PCM+offset] + 1./2.* wPHOS/wSum* eyStaPHOS[i-bin0PHOS]*eyStaPHOS[i-bin0PHOS],0.5);
					xSectionCombSysErr[i] = pow(1./2.* wPCM/wSum* eySysHPCM[i-bin0PCM]*eySysHPCM[i-bin0PCM] + 1./2.* wPHOS/wSum* eySysPHOS[i-bin0PHOS]*eySysPHOS[i-bin0PHOS],0.5);		
					cout<< " PHOS,PCM_L::"<< xSectionComb[i]<< " " << xSectionCombErr[i] << " " << xSectionCombStatErr[i] << " " << xSectionCombSysErr[i] << endl;
				}else{
					Double_t wPHOS = 1./eTot2PHOS[i-bin0PHOS];
					Double_t wPCM = 1./eTotH2PCM[i-bin0PCM];
					Double_t wSum = wPCM+wPHOS;
					xSectionComb[i] = (wPCM*yPCM[i-bin0PCM] +  wPHOS*yPHOS[i-bin0PHOS])/ wSum;
					xSectionCombErr[i] = pow((wPCM +  wPHOS),-0.5);    
					xSectionCombStatErr[i] = pow(1./2.* wPCM/wSum* eyStaPCM[i-bin0PCM+offset]*eyStaPCM[i-bin0PCM+offset] + 1./2.* wPHOS/wSum* eyStaPHOS[i-bin0PHOS]*eyStaPHOS[i-bin0PHOS],0.5);
					xSectionCombSysErr[i] = pow(1./2.* wPCM/wSum* eySysHPCM[i-bin0PCM]*eySysHPCM[i-bin0PCM] + 1./2.* wPHOS/wSum* eySysPHOS[i-bin0PHOS]*eySysPHOS[i-bin0PHOS],0.5);		
					cout<< " PHOS,PCM_H::"<< xSectionComb[i]<< " " << xSectionCombErr[i] << " " << xSectionCombStatErr[i] << " " << xSectionCombSysErr[i] << endl;
				}
			}
			xSectionCombErrL[i]= xSectionCombErr[i];
			xSectionCombErrH[i]= xSectionCombErr[i];
		}
	
		if ( okPHOS && !okPCM ){
			xEComb[i]=exSysPHOS[i-bin0PHOS];
			if( eTot2PHOS[i-bin0PHOS] !=0.){
				xSectionComb[i] = yPHOS[i-bin0PHOS];
				xSectionCombErr[i] = pow((eTot2PHOS[i-bin0PHOS]),0.5);
				xSectionCombStatErr[i] = eyStaPHOS[i-bin0PHOS];
				xSectionCombSysErr[i] = eySysPHOS[i-bin0PHOS];
				cout<< " PHOS_OK::"<< xSectionComb[i]<< " " << xSectionCombErr[i] << " "  << xSectionCombStatErr[i] << " " << xSectionCombSysErr[i] << endl;
			}
			xSectionCombErrL[i]= xSectionCombErr[i];
			xSectionCombErrH[i]= xSectionCombErr[i];
		}
	
		if ( !okPHOS && okPCM ){
			xEComb[i]=exSysPCM[i-bin0PCM];
			if( eTotL2PCM[i-bin0PCM] !=0. && eTotH2PCM[i-bin0PCM]){
				xSectionComb[i] = yPCM[i-bin0PCM];
				xSectionCombErr[i] = pow((eTotL2PCM[i-bin0PCM]),0.5);      // Asymmetric errors needed
				xSectionCombStatErr[i] = eyStaPCM[i-bin0PCM+offset];
				xSectionCombSysErr[i] = eySysHPCM[i-bin0PCM];
				cout<< " PCM_OK::"<< xSectionComb[i]<< " " << xSectionCombErr[i] << " "  << xSectionCombStatErr[i] << " " << xSectionCombSysErr[i] << endl;
			}
			xSectionCombErrL[i]= xSectionCombErr[i];
			xSectionCombErrH[i]= xSectionCombErr[i];
		}
	}
	
	graphStatComb = new TGraphAsymmErrors(nPtLimits-1,xComb,xSectionComb,xEComb,xEComb,xSectionCombStatErr,xSectionCombStatErr);
// 	cout << "statistical errors only" << endl;
// 	graphStatComb->Print();
	graphSystComb = new TGraphAsymmErrors(nPtLimits-1,xComb,xSectionComb,xEComb,xEComb,xSectionCombSysErr,xSectionCombSysErr);   
// 	cout << "systematic errors only" << endl;
// 	graphSystComb->Print();
// 	cout << "stat+sys errors" << endl;
	TGraphAsymmErrors* returnGraph = new TGraphAsymmErrors(nPtLimits-1,xComb,xSectionComb,xEComb,xEComb,xSectionCombErrL,xSectionCombErrH);
// 	returnGraph->Print();
	Int_t b = 0;
	while (xSectionComb[b] == 0){
		graphStatComb->RemovePoint(0);
		graphSystComb->RemovePoint(0);
		returnGraph->RemovePoint(0);
		b++;
	}
// 	cout << "statistical errors only" << endl;
// 	graphStatComb->Print();
// 	cout << "systematic errors only" << endl;
// 	graphSystComb->Print();
// 	cout << "stat+sys errors" << endl;
// 	returnGraph->Print();

	return returnGraph;
}

TGraphAsymmErrors* CombinePtPointsRAA( TGraphAsymmErrors* graphStatErrConv,		TGraphAsymmErrors* graphSystConv,
													TGraphAsymmErrors* graphStatErrPHOS,		TGraphAsymmErrors* graphSystPHOS,
 													TGraphAsymmErrors* &graphStatComb, 	TGraphAsymmErrors* &graphSystComb,  
													Double_t* xPtLimits,	Int_t nPtLimits,
													Int_t offset,			Int_t bin0PCM,	 Int_t bin0PHOS, Bool_t kRemoveLastPCMPoint=kFALSE, Int_t nBinsPCMRem = 1){
   
 	cout << "PHOS histogramm" << endl;
	graphStatErrPHOS->Print();
 	cout << "PCM histogramm" << endl;
	graphStatErrConv->Print();
	Double_t xComb[50];
	Double_t xEComb[50];
	Double_t xSectionComb[50];
	Double_t xSectionCombErr[50];
	Double_t xSectionCombErrL[50];
	Double_t xSectionCombErrH[50];
	Double_t xSectionCombStatErr[50];
	Double_t xSectionCombSysErr[50];
	
	Int_t nPHOS = graphSystPHOS->GetN();
	Double_t * xPHOS =  graphSystPHOS->GetX();
	Double_t * yPHOS =  graphSystPHOS->GetY();
	Double_t * eySysPHOS =  graphSystPHOS->GetEYlow();
	Double_t * exSysPHOS =  graphSystPHOS->GetEXlow();
	Double_t * eyStaPHOS =  graphStatErrPHOS->GetEYlow();
	Double_t * eTot2PHOS;
	eTot2PHOS = new Double_t[nPHOS];
	Double_t * eTotPHOS;
	eTotPHOS=new Double_t [nPHOS];
	
	for(Int_t i=0;i<nPHOS;i++){
		eTot2PHOS[i]=(eyStaPHOS[i]*eyStaPHOS[i]+eySysPHOS[i]*eySysPHOS[i]);
		eTotPHOS[i]=TMath::Sqrt( eTot2PHOS[i]);
 		cout<< "PHOS::"<< xPHOS[i]<< " "<<yPHOS[i]<< " " <<  eTotPHOS[i]<< endl;
	}
	
	Int_t nPCM;
	Double_t * xPCM;
	Double_t * yPCM;
	Double_t * exSysPCM;
	Double_t * eySysLPCM;
	Double_t * eySysHPCM;
	
	nPCM= graphSystConv->GetN();
	if (kRemoveLastPCMPoint) {
		nPCM= nPCM-nBinsPCMRem;
	}
	xPCM = graphSystConv->GetX();
	yPCM= graphSystConv->GetY();
	exSysPCM = graphSystConv->GetEXlow();
	eySysLPCM = graphSystConv->GetEYlow();
	eySysHPCM = graphSystConv->GetEYhigh();

	Double_t * eyStaPCM = graphStatErrConv->GetEYlow();

	 	cout << "PHOS sys" << endl;
	graphSystPHOS->Print();
 	cout << "PCM sys" << endl;
	graphSystConv->Print();

	
	Double_t eTotL2PCM[nPCM];
	Double_t eTotH2PCM[nPCM];
	Double_t eTotLPCM[nPCM];
	Double_t eTotHPCM[nPCM];
	for(Int_t i=0;i<nPCM;i++){
		eTotH2PCM[i]=(eyStaPCM[i+offset]*eyStaPCM[i+offset]+eySysHPCM[i]*eySysHPCM[i]);
		eTotHPCM[i]=TMath::Sqrt( eTotH2PCM[i]);
		eTotL2PCM[i]=(eyStaPCM[i+offset]*eyStaPCM[i+offset]+eySysLPCM[i]*eySysLPCM[i]);
		eTotLPCM[i]=TMath::Sqrt( eTotL2PCM[i]);
 		cout<< "PCM::"<< xPCM[i]<< " "<<yPCM[i]<< " " <<  eTotLPCM[i]<< " "<< eyStaPCM[i+offset]<<" "<< eySysHPCM[i]<< endl;
	}
	cout<<endl;  
	
	Bool_t okPHOS,okPCM;     
	for (Int_t i=0;i<nPtLimits-1;i++){
		Double_t xCenter = 0.5*(xPtLimits[i+1]+xPtLimits[i]);
		okPHOS=kFALSE;
		okPCM=kFALSE;
				
		xComb[i]=xCenter;

		xSectionComb[i] = 0;
		xSectionCombErrL[i]= 0;
		xSectionCombErrH[i]= 0;	

		Double_t decisionBoundary = 0.0000001;
		if((i-bin0PHOS)>=0){
 			cout << xPHOS[i-bin0PHOS] << "\t" << xPCM[i-bin0PCM] << "  xCenter::"<< i<< " " <<xCenter<< endl;
			if ( TMath::Abs(xPHOS[i-bin0PHOS] - xCenter) < decisionBoundary && yPHOS[i-bin0PHOS]!= 0.){
				okPHOS=kTRUE;
			}
		}
	
		if( (i-bin0PCM) >= 0){
			if (i-bin0PCM < nPCM){
				if ( TMath::Abs(xPCM[i-bin0PCM] - xCenter) < decisionBoundary && yPCM[i-bin0PCM]!= 0.){
					okPCM=kTRUE;
				}
			}
		}
	
		if ( okPHOS && okPCM ){
			xEComb[i]=exSysPCM[i-bin0PCM];
			if( eTotL2PCM[i-bin0PCM]!=0. &&  eTotH2PCM[i-bin0PCM]!=0. && eTot2PHOS[i-bin0PHOS] !=0.){
				if(yPCM[i-bin0PCM]> yPHOS[i-bin0PHOS]){
					Double_t wPHOS = 1./eTot2PHOS[i-bin0PHOS];
					Double_t wPCM = 1./eTotL2PCM[i-bin0PCM];
					Double_t wSum = wPCM+wPHOS;
					xSectionComb[i] = (wPCM*yPCM[i-bin0PCM] +  wPHOS*yPHOS[i-bin0PHOS])/ wSum;
					xSectionCombErr[i] = pow((wPCM +  wPHOS),-0.5);    
					xSectionCombStatErr[i] = pow(1./2.* wPCM/wSum* eyStaPCM[i-bin0PCM+offset]*eyStaPCM[i-bin0PCM+offset] + 1./2.* wPHOS/wSum* eyStaPHOS[i-bin0PHOS]*eyStaPHOS[i-bin0PHOS],0.5);
					xSectionCombSysErr[i] = pow(1./2.* wPCM/wSum* eySysHPCM[i-bin0PCM]*eySysHPCM[i-bin0PCM] + 1./2.* wPHOS/wSum* eySysPHOS[i-bin0PHOS]*eySysPHOS[i-bin0PHOS],0.5);		
					cout<< " PHOS,PCM_L::"<< xSectionComb[i]<< " " << xSectionCombErr[i] << " " << yPCM[i-bin0PCM]<< " "<< eTotLPCM[i-bin0PCM] << " "
					<< yPHOS[i-bin0PHOS]<<" "<< eTotPHOS[i-bin0PHOS]<< endl;
					cout << "stat & sys separated ::" << xSectionCombStatErr[i] << " " << xSectionCombSysErr[i] << endl;
				}else{
					Double_t wPHOS = 1./eTot2PHOS[i-bin0PHOS];
					Double_t wPCM = 1./eTotH2PCM[i-bin0PCM];
					Double_t wSum = wPCM+wPHOS;
					xSectionComb[i] = (wPCM*yPCM[i-bin0PCM] +  wPHOS*yPHOS[i-bin0PHOS])/ wSum;
					xSectionCombErr[i] = pow((wPCM +  wPHOS),-0.5);    
					xSectionCombStatErr[i] = pow(1./2.* wPCM/wSum* eyStaPCM[i-bin0PCM+offset]*eyStaPCM[i-bin0PCM+offset] + 1./2.* wPHOS/wSum* eyStaPHOS[i-bin0PHOS]*eyStaPHOS[i-bin0PHOS],0.5);
					xSectionCombSysErr[i] = pow(1./2.* wPCM/wSum* eySysHPCM[i-bin0PCM]*eySysHPCM[i-bin0PCM] + 1./2.* wPHOS/wSum* eySysPHOS[i-bin0PHOS]*eySysPHOS[i-bin0PHOS],0.5);		
					cout<< " PHOS,PCM_H::"<< xSectionComb[i]<< " " << xSectionCombErr[i] << " " << yPCM[i-bin0PCM]<< " "<< eTotHPCM[i-bin0PCM] << " "
					<< yPHOS[i-bin0PHOS]<<" "<< eTotPHOS[i-bin0PHOS] <<endl;
					cout << "stat & sys separated ::" << xSectionCombStatErr[i] << " " << xSectionCombSysErr[i] << endl;
				}
			}
			xSectionCombErrL[i]= xSectionCombErr[i];
			xSectionCombErrH[i]= xSectionCombErr[i];
		}
	
		if ( okPHOS && !okPCM ){
			xEComb[i]=exSysPHOS[i-bin0PHOS];
			if( eTot2PHOS[i-bin0PHOS] !=0.){
				xSectionComb[i] = yPHOS[i-bin0PHOS];
				xSectionCombErr[i] = pow((eTot2PHOS[i-bin0PHOS]),0.5);
				xSectionCombStatErr[i] = eyStaPHOS[i-bin0PHOS];
				xSectionCombSysErr[i] = eySysPHOS[i-bin0PHOS];
				cout<< " PHOS_OK::"<< xSectionComb[i]<< " " << xSectionCombErr[i] << " "  << yPHOS[i-bin0PHOS]<<" "<< eTotPHOS[i-bin0PHOS] <<endl;
				cout << "stat & sys separated ::" << xSectionCombStatErr[i] << " " << xSectionCombSysErr[i] << endl;
			}
			xSectionCombErrL[i]= xSectionCombErr[i];
			xSectionCombErrH[i]= xSectionCombErr[i];
		}
	
		if ( !okPHOS && okPCM ){
			xEComb[i]=exSysPCM[i-bin0PCM];
			if( eTotL2PCM[i-bin0PCM] !=0. && eTotH2PCM[i-bin0PCM]){
				xSectionComb[i] = yPCM[i-bin0PCM];
				xSectionCombErr[i] = pow((eTotL2PCM[i-bin0PCM]),0.5);      // Asymmetric errors needed
				xSectionCombStatErr[i] = eyStaPCM[i-bin0PCM+offset];
				xSectionCombSysErr[i] = eySysHPCM[i-bin0PCM];
				cout<< " PCM_OK::"<< xSectionComb[i]<< " " << xSectionCombErr[i] << " "  << yPCM[i-bin0PCM]<< " "<< eTotLPCM[i-bin0PCM] << " "<< endl;
				cout << "stat & sys separated ::" << xSectionCombStatErr[i] << " " << xSectionCombSysErr[i] << endl;
			}
			xSectionCombErrL[i]= xSectionCombErr[i];
			xSectionCombErrH[i]= xSectionCombErr[i];
		}
	}
	
	graphStatComb = new TGraphAsymmErrors(nPtLimits-1,xComb,xSectionComb,xEComb,xEComb,xSectionCombStatErr,xSectionCombStatErr);
	graphSystComb = new TGraphAsymmErrors(nPtLimits-1,xComb,xSectionComb,xEComb,xEComb,xSectionCombSysErr,xSectionCombSysErr);   
	TGraphAsymmErrors* returnGraph = new TGraphAsymmErrors(nPtLimits-1,xComb,xSectionComb,xEComb,xEComb,xSectionCombErrL,xSectionCombErrH);

	Int_t b = 0;
	while (xSectionComb[b] == 0){
		graphStatComb->RemovePoint(0);
		graphSystComb->RemovePoint(0);
		returnGraph->RemovePoint(0);
		b++;
	}

	cout << "statistical errors only" << endl;
	graphStatComb->Print();
	cout << "systematic errors only" << endl;
	graphSystComb->Print();
	cout << "stat+sys errors" << endl;
 	returnGraph->Print();
	
	return returnGraph;
}


void CalculateFitResults(TF1* fitstat, TF1* fitsys, Double_t* fitresults ,TString fitName = "", Double_t sigma = 1 ){
	Int_t nPar = fitstat->GetNpar();
	cout << nPar << endl;
	for (Int_t i =0; i < nPar; i++){
		if (fitName.CompareTo("Levy") == 0 && i == 0){
			fitresults[i*3] = fitsys->GetParameter(i)/sigma;
			fitresults[i*3+1] = fitstat->GetParError(i)/sigma;
			fitresults[i*3+2] = TMath::Sqrt(TMath::Abs(pow((fitsys->GetParError(i)/sigma),2)-pow((fitstat->GetParError(i)/sigma),2)));	
		} else {
			fitresults[i*3] = fitsys->GetParameter(i);
			fitresults[i*3+1] = fitstat->GetParError(i);
			fitresults[i*3+2] = TMath::Sqrt(TMath::Abs(pow((fitsys->GetParError(i)),2)-pow((fitstat->GetParError(i)),2)));	
		}
		cout << fitresults[i*3] << "\t" << fitresults[i*3+1] << "\t" << fitresults[i*3+2] << endl;
	}
	return;
}

void ReadOutFromSysErrVector(Double_t* readoutVector, Double_t* fillInVector, Double_t* fillInVectorError, Int_t cutNr, Int_t nrCutStudies,Int_t offsetAtEnd, Int_t numberOfLines,Double_t decisionBoundary){
	Double_t dummyValue;
	Double_t dummyError;
	Int_t l = cutNr*3-2;
	Int_t totalLineEntries = 1+ nrCutStudies*3+offsetAtEnd;
	for (Int_t i = 0; i < numberOfLines; i++){
		dummyValue = readoutVector[l];	
		dummyError = readoutVector[l+1];
		cout << "readout \t " << l << "\t" << dummyValue << "\t" << dummyError << endl;
		if (dummyError != 0 &&TMath::Abs(dummyValue)/dummyError > decisionBoundary){
			fillInVector[i] = dummyValue/TMath::Sqrt(2);	
			fillInVectorError[i] = fillInVector[i]*0.001;
		} else {
			fillInVector[i] = 0;	
			fillInVectorError[i] = 0;
		}
		l = l + totalLineEntries;
	}
	return;
}

void CalculateMeanSysErr(Double_t* fillInVector, Double_t* fillInVectorError, Double_t* posError, Double_t* negError, Int_t numberOfLines){
	Double_t dummyValue;
	Double_t dummyError;
	for (Int_t i = 0; i < numberOfLines; i++){
		dummyValue = (posError[i] + TMath::Abs(negError[i]))/2;	
		cout << posError[i] << "\t" << negError[i] << "\t" << dummyValue << endl;
		dummyError = dummyValue*0.001;
		fillInVector[i] = dummyValue;
		fillInVectorError[i] = dummyError;
	}
	return;
}

void CorrectSystematicErrorsWithMean(Double_t* oldErrorVector,Double_t* oldErrorVectorsErrors, Double_t* newErrorVector, Double_t* newErrorVectorError, Int_t numberOfPtBins){
	if (numberOfPtBins > 2) {	
		for (Int_t i = 0; i < numberOfPtBins; i++){
			if (i == 0){ // bin 1
				if (oldErrorVector[i] == 0 && oldErrorVector[i+1] == 0){
					newErrorVector[i] = oldErrorVector[i+2];
					newErrorVectorError[i] = oldErrorVectorsErrors[i+2];
				} else if (oldErrorVector[i] == 0){
					newErrorVector[i] = oldErrorVector[i+1];
					newErrorVectorError[i] = oldErrorVectorsErrors[i+1];
				} else {
					newErrorVector[i] = oldErrorVector[i];
					newErrorVectorError[i] = newErrorVector[i]*0.001;
				}
			} else if (i == 1){ // bin 2
				if (oldErrorVector[i] == 0 && oldErrorVector[i+1] != 0){
					newErrorVector[i] = oldErrorVector[i-1]*0.2 + oldErrorVector[i+1]*0.8;
					newErrorVectorError[i] = newErrorVector[i]*0.001;
					if (oldErrorVector[i-1] == 0 && newErrorVector[i-1] != 0){
						newErrorVector[i] = newErrorVector[i-1]*0.2 + oldErrorVector[i+1]*0.8 ;
						newErrorVectorError[i] = newErrorVector[i]*0.001;
					}
				} else if (oldErrorVector[i] == 0 && oldErrorVector[i+1] == 0){
					newErrorVector[i] = oldErrorVector[i-1]*0.2 + oldErrorVector[i+2]*0.8;
					newErrorVectorError[i] = newErrorVector[i]*0.001;
					if (oldErrorVector[i-1] == 0 && newErrorVector[i-1] != 0){
						newErrorVector[i] = (newErrorVector[i-1] + oldErrorVector[i+1])/2. ;
						newErrorVectorError[i] = newErrorVector[i]*0.001;
					}
				} else {
					newErrorVector[i] = oldErrorVector[i];
					newErrorVectorError[i] = newErrorVector[i]*0.001;
				}
			}else if (i == numberOfPtBins - 2 ){ // previous to the last bin
				if (oldErrorVector[i] == 0 && oldErrorVector[i+1] != 0){
					newErrorVector[i] = oldErrorVector[i-1]*0.8 + oldErrorVector[i+1]*0.2;
					newErrorVectorError[i] = newErrorVector[i]*0.001;
					if (oldErrorVector[i-1] == 0 && newErrorVector[i-1] != 0){
						newErrorVector[i] = newErrorVector[i-1]*0.8 + oldErrorVector[i+1]*0.2 ;
						newErrorVectorError[i] = newErrorVector[i]*0.001;
					}
				} else if (oldErrorVector[i] == 0 && oldErrorVector[i+1] == 0){
					newErrorVector[i] = (oldErrorVector[i-1] + oldErrorVector[i+1])/2. ;
					newErrorVectorError[i] = newErrorVector[i]*0.001;
					if (oldErrorVector[i-1] == 0 && newErrorVector[i-1] != 0){
						newErrorVector[i] = (newErrorVector[i-1] + oldErrorVector[i+1])/2. ;
						newErrorVectorError[i] = newErrorVector[i]*0.001;
					}
				} else {
					newErrorVector[i] = oldErrorVector[i];
					newErrorVectorError[i] = newErrorVector[i]*0.001;
				}
			} else if (i == numberOfPtBins - 3){ //3rd form end
				if (oldErrorVector[i] == 0 && oldErrorVector[i+1] != 0){
					if (TMath::Abs(oldErrorVector[i-1]) < TMath::Abs(oldErrorVector[i+1])){
						newErrorVector[i] = oldErrorVector[i-1];
						newErrorVectorError[i] = newErrorVector[i]*0.001;
					} else {
						newErrorVector[i] = oldErrorVector[i+1];
						newErrorVectorError[i] = newErrorVector[i]*0.001;
					}
				} else if (oldErrorVector[i] == 0 && oldErrorVector[i+1] == 0 && oldErrorVector[i+2] != 0){
					if (TMath::Abs(oldErrorVector[i-1]) < TMath::Abs(oldErrorVector[i+2])){
						newErrorVector[i] = oldErrorVector[i-1];
						newErrorVectorError[i] = newErrorVector[i]*0.001;
					} else {
						newErrorVector[i] = oldErrorVector[i+2];
						newErrorVectorError[i] = newErrorVector[i]*0.001;
					}
				} else {
					newErrorVector[i] = oldErrorVector[i];
					newErrorVectorError[i] = oldErrorVectorsErrors[i];
				}
			} else if (i == numberOfPtBins - 1 ){ //last
				if (oldErrorVector[i] == 0 && oldErrorVector[i-1] != 0){
					newErrorVector[i] = oldErrorVector[i-1];
					newErrorVectorError[i] = oldErrorVectorsErrors[i-1];				
				} else if (oldErrorVector[i] == 0 && oldErrorVector[i-1] == 0 && newErrorVector[i-1] != 0){
					newErrorVector[i] = newErrorVector[i-1];
					newErrorVectorError[i] = newErrorVectorError[i-1];				
				} else {
					newErrorVector[i] = oldErrorVector[i];
					newErrorVectorError[i] = newErrorVector[i]*0.001;
				}
			} else {	 //every where else
				if (oldErrorVector[i] == 0 && oldErrorVector[i+1] == 0 && i+2 != numberOfPtBins){
					newErrorVector[i] = (oldErrorVector[i-1] + oldErrorVector[i+2])/2. ;
					newErrorVectorError[i] = newErrorVector[i]*0.001;
					if (oldErrorVector[i-1] == 0 && newErrorVector[i-1] != 0){
						newErrorVector[i] = (newErrorVector[i-1] + oldErrorVector[i+2])/2. ;
						newErrorVectorError[i] = newErrorVector[i]*0.001;
					}
				} else if (oldErrorVector[i] == 0){
					newErrorVector[i] = (oldErrorVector[i-1] + oldErrorVector[i+1])/2. ;
					newErrorVectorError[i] = newErrorVector[i]*0.001;
					if (oldErrorVector[i-1] == 0 && newErrorVector[i-1] != 0){
						newErrorVector[i] = (newErrorVector[i-1] + oldErrorVector[i+1])/2. ;
						newErrorVectorError[i] = newErrorVector[i]*0.001;
					}
				} else {
					newErrorVector[i] = oldErrorVector[i];
					newErrorVectorError[i] = newErrorVector[i]*0.001;
				}	
			}		
		}
	} else {
		newErrorVector[0] = (oldErrorVector[0] + oldErrorVector[1])/2.;
		newErrorVector[1] = newErrorVector[0];
		newErrorVectorError[0] = newErrorVector[0]*0.001;
		newErrorVectorError[1] = newErrorVectorError[0];
	}
}

void ProduceGraphAsymmWithoutXErrors(TGraphAsymmErrors* inputgraph){
	Int_t n = inputgraph->GetN();
	Double_t* xValue = inputgraph->GetX();
	Double_t* xErrorHigh = inputgraph->GetEXhigh();
	Double_t* xErrorLow = inputgraph->GetEXlow();
	Double_t* yValue = inputgraph->GetY();
	Double_t* yErrorLow = inputgraph->GetEYlow();
	Double_t* yErrorHigh = inputgraph->GetEYhigh();
	for (Int_t i= 0; i < n; i++){
		xErrorHigh[i] = 0.;
		xErrorLow[i] = 0.;
	}
	inputgraph = new TGraphAsymmErrors(n,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
// 	inputgraph->Print();
}

void ProduceGraphPartialXErrors(TGraphErrors* inputgraph, Double_t part){  // kk
	Int_t n = inputgraph->GetN();
	Double_t* xValue = inputgraph->GetX();
	Double_t* xError = inputgraph->GetEX();
	Double_t* yValue = inputgraph->GetY();
	Double_t* yError = inputgraph->GetEY();
	for (Int_t i= 0; i < n; i++){
		xError[i] = part*xError[i];
		xError[i] = part*xError[i];
	}
	inputgraph = new TGraphErrors(n,xValue,yValue,xError,yError);
// 	inputgraph->Print();
}

void ProduceGraphAsymmPartialXErrors(TGraphAsymmErrors* inputgraph, Double_t part){  // kk
	Int_t n = inputgraph->GetN();
	Double_t* xValue = inputgraph->GetX();
	Double_t* xErrorHigh = inputgraph->GetEXhigh();
	Double_t* xErrorLow = inputgraph->GetEXlow();
	Double_t* yValue = inputgraph->GetY();
	Double_t* yErrorLow = inputgraph->GetEYlow();
	Double_t* yErrorHigh = inputgraph->GetEYhigh();
	for (Int_t i= 0; i < n; i++){
		xErrorHigh[i] = part*xErrorHigh[i];
		xErrorLow[i] = part*xErrorLow[i];


	}
	inputgraph = new TGraphAsymmErrors(n,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
// 	inputgraph->Print();
}

void ProduceGraphFixedXErrors(TGraphErrors* inputgraph, Double_t widthXError){ 
	Int_t n = inputgraph->GetN();
	Double_t* xValue = inputgraph->GetX();
	Double_t* xError = inputgraph->GetEX();
	Double_t* yValue = inputgraph->GetY();
	Double_t* yError = inputgraph->GetEY();
	for (Int_t i= 0; i < n; i++){
		xError[i] = widthXError/2;
		xError[i] = widthXError/2;
	}
	inputgraph = new TGraphErrors(n,xValue,yValue,xError,yError);
// 	inputgraph->Print();
}

void ProduceGraphAsymmFixedXErrors(TGraphAsymmErrors* inputgraph, Double_t widthXError){ 
	Int_t n = inputgraph->GetN();
	Double_t* xValue = inputgraph->GetX();
	Double_t* xErrorHigh = inputgraph->GetEXhigh();
	Double_t* xErrorLow = inputgraph->GetEXlow();
	Double_t* yValue = inputgraph->GetY();
	Double_t* yErrorLow = inputgraph->GetEYlow();
	Double_t* yErrorHigh = inputgraph->GetEYhigh();
	for (Int_t i= 0; i < n; i++){
		xErrorHigh[i] = widthXError/2;
		xErrorLow[i] = widthXError/2;


	}
	inputgraph = new TGraphAsymmErrors(n,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
// 	inputgraph->Print();
}

void ProduceGraphAsymmWithoutXYErrors(TGraphAsymmErrors* inputgraph){  // kk
	Int_t n = inputgraph->GetN();
	Double_t* xValue = inputgraph->GetX();
	Double_t* xErrorHigh = inputgraph->GetEXhigh();
	Double_t* xErrorLow = inputgraph->GetEXlow();
	Double_t* yValue = inputgraph->GetY();
	Double_t* yErrorLow = inputgraph->GetEYlow();
	Double_t* yErrorHigh = inputgraph->GetEYhigh();
	for (Int_t i= 0; i < n; i++){
		xErrorHigh[i] = 0.;
		xErrorLow[i] = 0.;
		yErrorHigh[i] = 0.;
		yErrorLow[i] = 0.;

	}
	inputgraph = new TGraphAsymmErrors(n,xValue,yValue,xErrorLow,xErrorHigh,yErrorLow,yErrorHigh);
// 	inputgraph->Print();
}

void ProduceGraphErrWithoutXErrors(TGraphErrors* inputgraph){
	Int_t n = inputgraph->GetN();
	Double_t* xValue = inputgraph->GetX();
	Double_t* xErrorHigh = inputgraph->GetEX();
	Double_t* yValue = inputgraph->GetY();
	Double_t* yErrorHigh = inputgraph->GetEY();
	for (Int_t i= 0; i < n; i++){
		xErrorHigh[i] = 0.;
	}
	inputgraph = new TGraphErrors(n,xValue,yValue,xErrorHigh,yErrorHigh);
// 	inputgraph->Print();
}

void ProduceGraphErrDisplacedX(TGraphErrors* inputgraph,Double_t displaceX){
	Int_t n = inputgraph->GetN();
	Double_t* xValue = inputgraph->GetX();
	Double_t* xErrorHigh = inputgraph->GetEX();
	Double_t* yValue = inputgraph->GetY();
	Double_t* yErrorHigh = inputgraph->GetEY();
	for (Int_t i= 0; i < n; i++){
	       xValue[i] = xValue[i]+displaceX;
		xErrorHigh[i] = 0.;

	}
	inputgraph = new TGraphErrors(n,xValue,yValue,xErrorHigh,yErrorHigh);
// 	inputgraph->Print();
}

TGraphAsymmErrors* Add2TGraphAsymmErrorsSameBinning(TGraphAsymmErrors* inputgraph1,TGraphAsymmErrors* inputgraph2){
	Double_t * xValue1 = inputgraph1->GetX(); 
	Double_t * yValue1 = inputgraph1->GetY();
	Double_t* xErrorLow1 = inputgraph1->GetEXlow();
	Double_t* xErrorHigh1 = inputgraph1->GetEXhigh();
	Double_t* yErrorLow1 = inputgraph1->GetEYlow();
	Double_t* yErrorHigh1 = inputgraph1->GetEYhigh();

	Double_t * yValue2 = inputgraph2->GetY();
	Double_t* yErrorLow2 = inputgraph2->GetEYlow();
	Double_t* yErrorHigh2 = inputgraph2->GetEYhigh();

	Double_t * yValue = inputgraph2->GetY();
	Double_t* yErrorLow = inputgraph2->GetEYlow();
	Double_t* yErrorHigh =inputgraph2->GetEYhigh();

	Int_t nPoints = inputgraph1->GetN();
	for (Int_t i = 0; i < nPoints; i++){
		yValue[i] = yValue1[i]+ yValue2[i];
		yErrorLow[i] = TMath::Sqrt(TMath::Power(yErrorLow1[i],2) + TMath::Power(yErrorLow2[i],2));
		yErrorHigh[i] = TMath::Sqrt(TMath::Power(yErrorHigh1[i],2) + TMath::Power(yErrorHigh2[i],2));
	}
	TGraphAsymmErrors* returnGraph =  new TGraphAsymmErrors(nPoints,xValue1,yValue,xErrorLow1,xErrorHigh1,yErrorLow,yErrorHigh); 
	return returnGraph;
}

TGraphAsymmErrors* CorrectTGraphAsymmErrorsToBinCenter(TGraphAsymmErrors* inputgraph1){
	Double_t * xValue1 = inputgraph1->GetX(); 
	Double_t * yValue1 = inputgraph1->GetY();
	Double_t* xErrorLow1 = inputgraph1->GetEXlow();
	Double_t* xErrorHigh1 = inputgraph1->GetEXhigh();
	Double_t* yErrorLow1 = inputgraph1->GetEYlow();
	Double_t* yErrorHigh1 = inputgraph1->GetEYhigh();

	Double_t * yValueNew = inputgraph1->GetY();
	Double_t* yErrorLowNew = inputgraph1->GetEYlow();
	Double_t* yErrorHighNew = inputgraph1->GetEYhigh();

	Int_t nPoints = inputgraph1->GetN();
	for (Int_t i = 0; i < nPoints; i++){
		yValueNew[i] = yValue1[i]/xValue1[i];
		yErrorLowNew[i] = yErrorLow1[i]/xValue1[i];
		yErrorHighNew[i] = yErrorHigh1[i]/xValue1[i];
	}
	TGraphAsymmErrors* returnGraph =  new TGraphAsymmErrors(nPoints,xValue1,yValueNew,xErrorLow1,xErrorHigh1,yErrorLowNew,yErrorHighNew); 
	return returnGraph;
}

TH1D* ShortChargedHadronHisto(TH1D* histIn) {
	const Int_t nPtBins = 57;
	Double_t xBins[nPtBins+1] = 
		{0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
		0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
		1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
		2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
		4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
		11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0};

	
	TH1D* hist = new TH1D("hist", "", nPtBins, xBins);
	hist->SetTitle(histIn->GetTitle());
	hist->GetXaxis()->SetTitle(histIn->GetXaxis()->GetTitle());
	hist->GetYaxis()->SetTitle(histIn->GetYaxis()->GetTitle());

	const Double_t deltapt = 0.0001;

	for(Int_t bin = 1; bin <= nPtBins; bin++) {

		// check bin size
		if(TMath::Abs(hist->GetXaxis()->GetBinLowEdge(bin) -
			histIn->GetXaxis()->GetBinLowEdge(bin)) > deltapt) {
		cout << "pt edge low does not agree!!!!!!!!!!!!!!!" << endl;
		return 0;
		}
		if(TMath::Abs(hist->GetXaxis()->GetBinUpEdge(bin) -
			histIn->GetXaxis()->GetBinUpEdge(bin)) > deltapt) {
		cout << "pt edge high does not agree!!!!!!!!!!!!!!!" << endl;
		return 0;
		}
		if(TMath::Abs(hist->GetXaxis()->GetBinCenter(bin) -
			histIn->GetXaxis()->GetBinCenter(bin)) > deltapt) {
		cout << "pt center does not agree!!!!!!!!!!!!!!!" << endl;
		return 0;
		}

		hist->SetBinContent(bin, histIn->GetBinContent(bin));
		cout << "\n converting n\n" << endl;
		cout << "pt in   " << histIn->GetXaxis()->GetBinCenter(bin) << "    bin content in  " << histIn->GetBinContent(bin) << endl;
		cout << " pt out   " << hist->GetXaxis()->GetBinCenter(bin) << "    bin content out   " << hist->GetBinContent(bin) << endl;


		hist->SetBinError(bin, histIn->GetBinError(bin));
	}

	return hist;
}

TH1D* ConvertChargedHadronHisto(TH1D* histIn, TH1D* fraction ) {
	const Int_t nPtBins = 57;
	Double_t xBins[nPtBins+1] = 
		{0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
		0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
		1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
		2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
		4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
		11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0};

	
	TH1D* hist = new TH1D("hist", "", nPtBins, xBins);
	hist->SetTitle(histIn->GetTitle());
	hist->GetXaxis()->SetTitle(histIn->GetXaxis()->GetTitle());
	hist->GetYaxis()->SetTitle(histIn->GetYaxis()->GetTitle());

	const Double_t deltapt = 0.0001;

	for(Int_t bin = 1; bin <= nPtBins; bin++) {

		// check bin size
		if(TMath::Abs(hist->GetXaxis()->GetBinLowEdge(bin) -
			histIn->GetXaxis()->GetBinLowEdge(bin)) > deltapt) {
		cout << "pt edge low does not agree!!!!!!!!!!!!!!!" << endl;
		return 0;
		}
		if(TMath::Abs(hist->GetXaxis()->GetBinUpEdge(bin) -
			histIn->GetXaxis()->GetBinUpEdge(bin)) > deltapt) {
		cout << "pt edge high does not agree!!!!!!!!!!!!!!!" << endl;
		return 0;
		}
		if(TMath::Abs(hist->GetXaxis()->GetBinCenter(bin) -
			histIn->GetXaxis()->GetBinCenter(bin)) > deltapt) {
		cout << "pt center does not agree!!!!!!!!!!!!!!!" << endl;
		return 0;
		}

		if(TMath::Abs(hist->GetXaxis()->GetBinLowEdge(bin) -
			fraction->GetXaxis()->GetBinLowEdge(bin)) > deltapt) {
		cout << "pt edge fraction low does not agree!!!!!!!!!!!!!!!" << endl;
		return 0;
		}
		if(TMath::Abs(hist->GetXaxis()->GetBinUpEdge(bin) -
			fraction->GetXaxis()->GetBinUpEdge(bin)) > deltapt) {
		cout << "pt edge fraction high does not agree!!!!!!!!!!!!!!!" << endl;
		return 0;
		}
		if(TMath::Abs(hist->GetXaxis()->GetBinCenter(bin) -
			fraction->GetXaxis()->GetBinCenter(bin)) > deltapt) {
		cout << "pt center fraction does not agree!!!!!!!!!!!!!!!" << endl;
		return 0;
		}

		Double_t content = histIn->GetBinContent(bin)*fraction->GetBinContent(bin);
		Double_t error = TMath::Sqrt(pow(histIn->GetBinError(bin)*fraction->GetBinContent(bin),2)+ pow(histIn->GetBinContent(bin)*fraction->GetBinError(bin),2) );     

		hist->SetBinContent(bin, content);
		cout << "\n converting n\n" << endl;
		cout << "pt in   " << histIn->GetXaxis()->GetBinCenter(bin) << "    bin content in  " << histIn->GetBinContent(bin) << endl;
		cout << "pt frac   " << fraction->GetXaxis()->GetBinCenter(bin) << "    bin content fraction  " << fraction->GetBinContent(bin) << endl;
		cout << " pt out   " << hist->GetXaxis()->GetBinCenter(bin) << "    bin content out   " << hist->GetBinContent(bin) << endl;
		hist->SetBinError(bin, error);
	}

	return hist;
}

TF1 *BinShiftTH1D(TH1D *UnshiftedYield, TH1D **ShiftedYield, TString mesonType =  "Pi0" ,TString BinShiftType="L", TString NameBinShiftHist = "DummyFit",Double_t minPtForFitsDummy = 0.6, Double_t* parameters = 0x00){
	TF1 *CurrentFit = new TF1();
	TH1D *CurrentHist  = (TH1D*)UnshiftedYield->Clone(""); 
	(*ShiftedYield) = (TH1D*)CurrentHist->Clone(NameBinShiftHist);
		
	Int_t binNumber = CurrentHist->GetXaxis()->GetNbins();
	Double_t globalRatioBS = 0;
	Double_t testGlobalRatioBS = 1;
	Int_t colorBS = 1;
	Double_t maxPtPi0BS = CurrentHist->GetXaxis()->GetBinUpEdge(CurrentHist->GetNbinsX());
	
	if(BinShiftType.BeginsWith("h") || BinShiftType.BeginsWith("H")){
		CurrentFit = FitObject("h","fitBinShiftingPi0",mesonType.Data(),CurrentHist,minPtForFitsDummy,maxPtPi0BS,parameters,"QNRME+");
	}
	else if(BinShiftType.BeginsWith("l") || BinShiftType.BeginsWith("L")){
		CurrentFit = FitObject("l","fitBinShiftingPi0",mesonType.Data(),CurrentHist,minPtForFitsDummy,maxPtPi0BS,parameters,"QNRME+");
	}
	else if(BinShiftType.BeginsWith("p") || BinShiftType.BeginsWith("P")){
		CurrentFit = FitObject("p","fitBinShiftingPi0",mesonType.Data(),CurrentHist,minPtForFitsDummy,maxPtPi0BS,parameters,"QNRME+");
	}
	else if(BinShiftType.BeginsWith("rad") || BinShiftType.BeginsWith("RAD")){
		CurrentFit = FitObject("rad","fitBinShiftingPi0",mesonType.Data(),CurrentHist,minPtForFitsDummy,maxPtPi0BS,parameters,"QNRME+");
	}
	else if(BinShiftType.BeginsWith("qcd") || BinShiftType.BeginsWith("QCD")){
		CurrentFit = FitObject("qcd","fitBinShiftingPi0",mesonType.Data(),CurrentHist,minPtForFitsDummy,maxPtPi0BS,parameters,"QNRME+");
	}
	else if(BinShiftType.BeginsWith("xqcd") || BinShiftType.BeginsWith("XQCD")){
		CurrentFit = FitObject("xqcd","fitBinShiftingPi0",mesonType.Data(),CurrentHist,minPtForFitsDummy,maxPtPi0BS,parameters,"QNRME+");
	}
	CurrentFit->SetRange(minPtForFitsDummy,maxPtPi0BS);
	
	while(globalRatioBS != testGlobalRatioBS){	
		if(colorBS == 200) break;
		testGlobalRatioBS = globalRatioBS;
		globalRatioBS = 0.;
		colorBS++;
		// fit acc+eff corrected yield with Hagedorn function
		(*ShiftedYield)->Fit(CurrentFit,"QRM0");
	
		// apply bin shift correction
		for (Int_t ib=2; ib<=binNumber; ib++) {
			Double_t ptMin = CurrentHist->GetBinLowEdge(ib);
			Double_t ptMax = ptMin + CurrentHist->GetBinWidth(ib);
			// the bin shift affected value of the fit function in current bin
			Double_t shiftedValue = CurrentFit->Integral(ptMin,ptMax) / (ptMax - ptMin);
			// the correct value at the bin center
			Double_t trueValue = CurrentFit->Eval((ptMax + ptMin)/2.);
				
		
			// the bin shift correction factor
			Double_t ratio =  shiftedValue / trueValue;
				
			(*ShiftedYield)->SetBinContent(ib,  CurrentHist->GetBinContent(ib) / ratio);
			(*ShiftedYield)->SetBinError(ib,  CurrentHist->GetBinError(ib) / ratio);
			globalRatioBS = globalRatioBS + ratio;
		}
		
		globalRatioBS = globalRatioBS/(binNumber-1);
		
		(*ShiftedYield)->SetMarkerStyle(24);
		(*ShiftedYield)->SetMarkerSize(0.9);
		(*ShiftedYield)->SetMarkerColor(kBlack);
		
		DrawGammaSetMarkerTF1( CurrentFit, 1, 0.4, kBlue-4);
			
		Double_t parameter[10];
		CurrentFit->GetParameters(parameter);
		CurrentFit->SetParameter(0,parameter[0]);
		CurrentFit->SetParameter(1,parameter[1]);
		CurrentFit->SetParameter(2,parameter[2]);
		CurrentFit->SetParameter(3,parameter[3]);
		CurrentFit->SetParameter(4,parameter[4]);
		CurrentFit->SetParameter(5,parameter[5]);
		CurrentFit->SetParameter(6,parameter[6]);
		CurrentFit->SetParameter(7,parameter[7]);
		CurrentFit->SetParameter(8,parameter[8]);
		CurrentFit->SetParameter(9,parameter[9]);
		
		cout<<colorBS<<" ";
	}
// 	cout << WriteParameterToFile(CurrentFit)<< endl;	
	return CurrentFit;
	
}

TF1 *ApplyYShift(TGraphAsymmErrors *UnshiftedYield, TGraphAsymmErrors **ShiftedYield, TString BinShiftType, TString NameBinShiftHist,Double_t minPtForFitsDummy = 0.6, Double_t* parameters = 0x00, Double_t accuracy = 0.0001, Bool_t incParlimits = kFALSE){
// 	cout << "entered Y shift" << endl;
	TF1 *CurrentFit = new TF1();
	TGraphAsymmErrors *CurrentGraph  = (TGraphAsymmErrors*)UnshiftedYield->Clone(""); 
	(*ShiftedYield) = (TGraphAsymmErrors*)CurrentGraph->Clone(NameBinShiftHist);	
		
	Int_t binNumber = CurrentGraph->GetN();
	Double_t* xvalue = UnshiftedYield->GetX();
	Double_t* xvalueErr = UnshiftedYield->GetEXlow();
	Double_t *yValue     = UnshiftedYield->GetY();
	Double_t *yValueErrlow  = UnshiftedYield->GetEYlow();
// 	Double_t *yValueErrhigh = UnshiftedYield->GetEYhigh();	
// 	Double_t xbins[binNumber+1];
	
	Double_t globalRatioBS = 0;
	Double_t testGlobalRatioBS = 1;
	Int_t colorBS = 1;
	Double_t maxPtPi0BS = xvalue[binNumber] + xvalueErr[binNumber];
	
	
	if(BinShiftType.BeginsWith("h") || BinShiftType.BeginsWith("H")){
		CurrentFit = FitObject("h","fitBinShiftingPi0","Pi0",CurrentGraph,minPtForFitsDummy,maxPtPi0BS,parameters);
	}
	else if(BinShiftType.BeginsWith("l") || BinShiftType.BeginsWith("L")){
		CurrentFit = FitObject("l","fitBinShiftingPi0","Pi0",CurrentGraph,minPtForFitsDummy,maxPtPi0BS,parameters);
	}
	else if(BinShiftType.BeginsWith("p") || BinShiftType.BeginsWith("P")){
		CurrentFit = FitObject("p","fitBinShiftingPi0","Pi0",CurrentGraph,minPtForFitsDummy,maxPtPi0BS,parameters);
	}
	else if(BinShiftType.BeginsWith("rad") || BinShiftType.BeginsWith("RAD")){
		if (!incParlimits){
			CurrentFit = FitObject("rad","fitBinShiftingPi0","Pi0",CurrentGraph,minPtForFitsDummy,maxPtPi0BS,parameters);
		} else {
			CurrentFit = FitObject("rad","fitBinShiftingPi0","Pi0",CurrentGraph,minPtForFitsDummy,maxPtPi0BS);
			SetParametersLimitsForFit(CurrentFit, 5, parameters);
		}	
	}
	else if(BinShiftType.BeginsWith("qcd") || BinShiftType.BeginsWith("QCD")){
		CurrentFit = FitObject("qcd","fitBinShiftingPi0","Pi0",CurrentGraph,minPtForFitsDummy,maxPtPi0BS,parameters);
	}
	else if(BinShiftType.BeginsWith("xqcd") || BinShiftType.BeginsWith("XQCD")){
		CurrentFit = FitObject("xqcd","fitBinShiftingPi0","Pi0",CurrentGraph,minPtForFitsDummy,maxPtPi0BS,parameters);
	}

	CurrentFit->SetRange(minPtForFitsDummy,maxPtPi0BS);
	
	
	
	Double_t ratio[binNumber+1];	
	while(TMath::Abs(globalRatioBS - testGlobalRatioBS) > accuracy){
		if(colorBS == 200) break;
		testGlobalRatioBS = globalRatioBS;
		globalRatioBS = 0.;
		colorBS++;
		(*ShiftedYield)->Fit(CurrentFit,"QRM0");
		
// 		cout << WriteParameterToFile(CurrentFit) << endl;
		Double_t *xPoint     = (*ShiftedYield)->GetX();
		Double_t *yPoint     = (*ShiftedYield)->GetY();
		
		Double_t *errorXlow  = (*ShiftedYield)->GetEXlow();
// 		Double_t *errorXhigh = (*ShiftedYield)->GetEXhigh();
		Double_t *errorYlow  = (*ShiftedYield)->GetEYlow();
// 		Double_t *errorYhigh = (*ShiftedYield)->GetEYhigh();	
	
		// apply bin shift correction
		for (Int_t ib=0; ib<binNumber; ib++) {
			Double_t ptMin = xPoint[ib]-errorXlow[ib];
			Double_t ptMax = xPoint[ib]+errorXlow[ib];
			
			// the bin shift affected value of the fit function in current bin
			Double_t shiftedValue = CurrentFit->Integral(ptMin,ptMax) / (ptMax - ptMin);
			// the correct value at the bin center
			Double_t trueValue = CurrentFit->Eval((ptMax + ptMin)/2.);
			// the bin shift correction factor
			ratio[ib] =  shiftedValue / trueValue;
			errorYlow[ib] = 	yValueErrlow[ib]/ ratio[ib];
			errorYlow[ib] = 	yValueErrlow[ib]/ ratio[ib];
			yPoint[ib] = yValue[ib]/ratio[ib];
			(*ShiftedYield)->SetMarkerStyle(24);
			(*ShiftedYield)->SetMarkerSize(0.9);
			(*ShiftedYield)->SetMarkerColor(kBlack);
			globalRatioBS = globalRatioBS + ratio[ib];
		}
		
		globalRatioBS = globalRatioBS/(binNumber);
		
		DrawGammaSetMarkerTF1( CurrentFit, 1, 0.4, kBlue-4);
			
		Double_t parameter[10];
		CurrentFit->GetParameters(parameter);
		
		CurrentFit->SetParameter(0,parameter[0]);
		CurrentFit->SetParameter(1,parameter[1]);
		CurrentFit->SetParameter(2,parameter[2]);
		CurrentFit->SetParameter(3,parameter[3]);
		CurrentFit->SetParameter(4,parameter[4]);
		CurrentFit->SetParameter(5,parameter[5]);
		CurrentFit->SetParameter(6,parameter[6]);
		CurrentFit->SetParameter(7,parameter[7]);
		CurrentFit->SetParameter(8,parameter[8]);
		CurrentFit->SetParameter(9,parameter[9]);
// 		cout<<colorBS<<" " << globalRatioBS << "  " << testGlobalRatioBS << endl;
		
	}
// 	for (Int_t ib=0; ib<=binNumber; ib++) {
// 		cout << ratio[ib] << "\t" << endl;
// 	}
	return CurrentFit;
	
}

//_________________________________________________________________________________________________
TGraphAsymmErrors* ApplyYshiftIndividualSpectra(TGraphAsymmErrors *IndividualSpectum , TF1 *commonFit) {
// 	cout << "entered bin shifting individual spectra y  "<<IndividualSpectum->GetName()  << endl;
	TGraphAsymmErrors*	IndividualSpectumShifted = (TGraphAsymmErrors*)IndividualSpectum->Clone();
   TGraphAsymmErrors*   dummySpec = (TGraphAsymmErrors*)IndividualSpectum->Clone();
   
	Double_t* xvalueI = IndividualSpectumShifted->GetX();
	Double_t* xvalueIErrlow = IndividualSpectumShifted->GetEXlow();
	Double_t* xvalueIErrhigh = IndividualSpectumShifted->GetEXhigh();
	Double_t* yvalueI = IndividualSpectumShifted->GetY();
	Double_t* yvalueIErrlow = dummySpec->GetEYlow();
	Double_t* yvalueIErrhigh = dummySpec->GetEYhigh();
	Int_t numberPointsI = IndividualSpectumShifted->GetN();
   
	for (Int_t ip = 0; ip < numberPointsI; ip++){
//       cout << "before:  " << yvalueIErrlow[ip]/yvalueI[ip] << "\t" << yvalueIErrhigh[ip]/yvalueI[ip] << endl;
		Double_t ptMin = xvalueI[ip]-xvalueIErrlow[ip];
		Double_t ptMax = xvalueI[ip]+xvalueIErrhigh[ip];
      yvalueIErrlow[ip]  = yvalueIErrlow[ip]/yvalueI[ip];
		yvalueIErrhigh[ip]  = yvalueIErrhigh[ip]/yvalueI[ip];	
		// the bin shift affected value of the fit function in current bin
		Double_t shiftedValue = commonFit->Integral(ptMin,ptMax) / (ptMax - ptMin);
		// the correct value at the bin center
		Double_t trueValue = commonFit->Eval((ptMax + ptMin)/2.);
		// the bin shift correction factor
		Double_t ratio =  shiftedValue / trueValue;
      Double_t output = yvalueI[ip]/ratio;
		IndividualSpectumShifted->SetPoint(ip,xvalueI[ip],yvalueI[ip]/ratio);
		IndividualSpectumShifted->SetPointEYlow(ip,yvalueIErrlow[ip]*output);
		IndividualSpectumShifted->SetPointEYhigh(ip,yvalueIErrhigh[ip]*output);		
	}
	
// 	IndividualSpectumShifted->Print();
 	return IndividualSpectumShifted;
}

//_______________________________________________________________________________________
Double_t bin_shift_x(TF1 *fYield, Double_t ptMin, Double_t ptMax, Double_t yValue){

	//
	// sample code for bin shift in x direction
	// Klaus Reygers
	//

	// define function describing yield dN/dpT 
	// parameter [3] needed for root finding (f(pTlw) - meas = 0)

	// define the bin limits
	Float_t pT1 = ptMin; // lower bin limit in GeV/c
	Float_t pT2 = ptMax; // upper bin limit in GeV/c

	// integrate yield of that bin
	Float_t meas = 1./(pT2 - pT1) * fYield->Integral(pT1,pT2);

	// Set offset to integrated yield (variable "meas") 
	// in order to solve the eq. f(pTlw) - meas = 0
	fYield->SetParameter(3, meas);
	
	// prepare root finding with root's BrentRootFinder
	ROOT::Math::WrappedTF1 wfYield(*fYield); // create wrapper function
	ROOT::Math::BrentRootFinder brf;         // create root finder
	brf.SetFunction(wfYield, pT1, pT2);       // set range for function wfYield
	
	// solve the equation f(pTlw) - meas = 0
	brf.Solve(); 
	Double_t pt0 = brf.Root();

	// correct position of the point according to 
	// Lafferty, Wyatt., Nucl. Instr. and Meth. A 355, 541, 1995
	// printf("%4.1f < pt < %4.1f GeV/c, pt0 = %6.3f GeV/c, y=%g, int=%g\n",ptMin,ptMax,pt0,yValue,meas);
	return pt0;
	
   if(yValue == 0) {}
}

//__________________________ApplyXshift___________________________________________________
TGraphAsymmErrors *ApplyXshift(TGraphAsymmErrors *spectrum, TF1 *tsallis){
	//----------------------------------------------------------------------
	// This function takes a spectrum, fits it by a function tsallis
	// and calculates pt-shift iteratively.
	// Return value: a new spectrum with pt-shifts for each point.
	// Yuri Kharlov. 19.12.2011
	// modified by Friederike Bock : 19.12.2011
	//----------------------------------------------------------------------

// 	cout << "entered bin shifting x" << endl;
	TGraphAsymmErrors *spectrumShift = (TGraphAsymmErrors*)spectrum->Clone();
	spectrumShift->SetName(Form("%s_xshift",spectrum->GetName()));
	Int_t nIter = 10;
// 	cout << "loaded spectrum" << endl;
	RemoveScalingWithPtGraph(spectrumShift);
  
	// Define a function Int(f)-f0

//   	Double_t* xvalue = spectrum->GetX();
// 	Double_t* xvalueErr = spectrum->GetEXlow();
	Int_t numberPoints = spectrum->GetN();
// 	cout << "number of bin in x" << numberPoints<< endl;
// 	cout << "red x values from spectrum" << endl;
		
	TString formula = Form("%s - [3]",tsallis->GetName());
	TF1 * fYield = new TF1("dTsalis", formula, 0.2,25.) ;

	for (Int_t iter=0; iter<nIter; iter++) {
		Double_t *xPoint     = spectrumShift->GetX();
		Double_t *yPoint     = spectrumShift->GetY();
		
		Double_t *errorXlow  = spectrumShift->GetEXlow();
		Double_t *errorXhigh = spectrumShift->GetEXhigh();
		Double_t *errorYlow  = spectrumShift->GetEYlow();
		Double_t *errorYhigh = spectrumShift->GetEYhigh();

		
		spectrumShift->Fit(tsallis, "NRMEXQ0","", xPoint[0]-errorXlow[0],xPoint[numberPoints-1]+errorXhigh[numberPoints-1]);
		fYield->SetParameter(0,tsallis->GetParameter(0));
		fYield->SetParameter(1,tsallis->GetParameter(1));
		fYield->SetParameter(2,tsallis->GetParameter(2));
		
// 		cout << "iteration " <<iter << endl;
		for (Int_t i=0; i<numberPoints; i++) {
			Double_t ptMin = xPoint[i]-errorXlow[i];
			Double_t ptMax = xPoint[i]+errorXhigh[i];
			fYield->SetParameter(3,0.);
			Double_t pt0 = bin_shift_x(fYield,ptMin,ptMax,yPoint[i]);
			Double_t dX = xPoint[i] - pt0;
			spectrumShift->SetPoint(i,pt0,yPoint[i]);
			spectrumShift->SetPointEXlow (i,errorXlow[i] -dX);
			spectrumShift->SetPointEXhigh(i,errorXhigh[i]+dX);
			if (iter == nIter-1) {
				cout << Form("%4.1f<pt<%4.1f GeV/c, pt0=%6.3f GeV/c, Ed3s/dp3 = %.3g + %.3g - %.3g pb/GeV2",
				ptMin,ptMax,pt0,yPoint[i],errorYhigh[i],errorYlow[i]) << endl;
			}
		}
	}
	
// 	cout<<"%Integral measured::"<<100.*tsallis->Integral(xvalue[0]-xvalueErr[0],xvalue[numberPoints-1]+xvalueErr[numberPoints-1])/tsallis->Integral(0.,100.)<<endl;
// 	cout<<"%Integral NOT measured::"<<100.*(1.-tsallis->Integral(xvalue[0]-xvalueErr[0],xvalue[numberPoints-1]+xvalueErr[numberPoints-1])/tsallis->Integral(0.,100.)) <<endl;

	
	ScaleWithPtGraph(spectrumShift);
// 	spectrumShift->Print();
	return spectrumShift;
}

//_________________________________________________________________________________________________
TGraphAsymmErrors* ApplyXshiftIndividualSpectra(TGraphAsymmErrors *CombinedSpectrum, TGraphAsymmErrors *IndividualSpectum , TF1 *tsallis, Int_t startBin =0, Int_t numberOfCommonBins = 0){
// 	CombinedSpectrum->Print();
// 	cout << "entered bin shifting individual spectra x" << endl;
	TGraphAsymmErrors *spectrumShifted = (TGraphAsymmErrors*)CombinedSpectrum->Clone();
  	Double_t* xvalue = spectrumShifted->GetX();
// 	Double_t* yvalue = spectrumShifted->GetY();
	Double_t* xvalueErrUp = spectrumShifted->GetEXhigh();
	Double_t* xvalueErrLow = spectrumShifted->GetEXlow();
		
	TGraphAsymmErrors*	IndividualSpectumShifted = (TGraphAsymmErrors*)IndividualSpectum->Clone();
// 	IndividualSpectum->Print();
	Double_t* xvalueI = IndividualSpectumShifted->GetX();
	
	Double_t* xvalueIErrlow = IndividualSpectumShifted->GetEXlow();
	Double_t* xvalueIErrhigh = IndividualSpectumShifted->GetEXhigh();
	Double_t* yvalueI = IndividualSpectumShifted->GetY();
	Double_t* yvalueIErrlow = IndividualSpectumShifted->GetEYlow();
	Double_t* yvalueIErrhigh = IndividualSpectumShifted->GetEYhigh();
	Int_t numberPointsI = IndividualSpectumShifted->GetN();
		
// 	IndividualSpectumShifted->Print();
	Int_t iBinsInd = 0;
	for (Int_t iBins = startBin ; iBins < numberOfCommonBins; iBins++){
// 		cout << iBinsInd << "\t"<<xvalue[iBins] << "\t"<< xvalueI[iBinsInd]<< endl;
// 		cout << "Before Shift" << yvalue[iBinsInd] << "\t" << yvalueI[iBinsInd] << endl;
		yvalueI[iBinsInd] = yvalueI[iBinsInd]*xvalueI[iBinsInd]/xvalue[iBins];
// 		cout << "after Shift" << yvalue[iBinsInd] << "\t" << yvalueI[iBinsInd] << endl;
		yvalueIErrhigh[iBinsInd] = yvalueIErrhigh[iBinsInd]*xvalueI[iBinsInd]/xvalue[iBins];
		yvalueIErrlow[iBinsInd] = yvalueIErrlow[iBinsInd]*xvalueI[iBinsInd]/xvalue[iBins];
		xvalueI[iBinsInd] = xvalue[iBins];
		xvalueIErrhigh[iBinsInd] = xvalueErrUp[iBins];
		xvalueIErrlow[iBinsInd] = xvalueErrLow[iBins];
		iBinsInd++;
	}
// 	cout << iBinsInd << "\t" << numberPointsI << endl;
	
	TString formula = Form("%s - [3]",tsallis->GetName());
	TF1 * fYield = new TF1("dTsalis", formula, 0.2,25.) ;
	fYield->SetParameter(0,tsallis->GetParameter(0));
	fYield->SetParameter(1,tsallis->GetParameter(1));
	fYield->SetParameter(2,tsallis->GetParameter(2));

	for (Int_t ip = (iBinsInd-1); ip < numberPointsI; ip++){
		Double_t ptMin = xvalueI[ip]-xvalueIErrlow[ip];
		Double_t ptMax = xvalueI[ip]+xvalueIErrhigh[ip];
		fYield->SetParameter(3,0.);
		Double_t pt0 = bin_shift_x(fYield,ptMin,ptMax,yvalueI[ip]*xvalueI[ip]);
		Double_t dX = xvalueI[ip] - pt0;
		IndividualSpectumShifted->SetPoint(ip,pt0,yvalueI[ip]*xvalueI[ip]/pt0);
		IndividualSpectumShifted->SetPointEXlow(ip,xvalueIErrlow[ip] -dX);
		IndividualSpectumShifted->SetPointEXhigh(ip,xvalueIErrhigh[ip]+dX);		
	}
	
// 	IndividualSpectumShifted->Print();
 	return IndividualSpectumShifted;
}

TString ReturnDateString(){
	TDatime 	today;
	int 		iDate = today.GetDate();
	int 		iYear=iDate/10000;
	int 		iMonth=(iDate%10000)/100;
	int 		iDay=iDate%100;
	TString 	cMonth[12]={"Jan","Feb","Mar","Apr","May","Jun",
		    "Jul","Aug","Sep","Oct","Nov","Dec"};
	TString 	textDayth;
	if (iDay== 11){
		textDayth = "st";
	} else if  (iDay== 12){
		textDayth = "th";
	} else if  (iDay== 13){
		textDayth = "th";
	} else if  (iDay%10 == 1){
		textDayth = "st";
	} else if (iDay%10 == 2){
		textDayth = "nd";
	} else if (iDay%10 == 3){
		textDayth = "rd";
	} else {
		textDayth = "th";
	}
	return Form("%i^{%s} %s %i",iDay, textDayth.Data(),cMonth[iMonth-1].Data(), iYear);
}

TString ReturnDateStringForOutput(){
	TDatime 	today;
	int 		iDate = today.GetDate();
	int 		iYear=iDate/10000;
	int 		iMonth=(iDate%10000)/100;
	int 		iDay=iDate%100;
	TString 	cMonth[12]={"Jan","Feb","Mar","Apr","May","Jun",
		    "Jul","Aug","Sep","Oct","Nov","Dec"};
	TString 	textDayth;
	if (iDay== 11){
		textDayth = "st";
	} else if  (iDay== 12){
		textDayth = "th";
	} else if  (iDay== 13){
		textDayth = "th";
	} else if  (iDay%10 == 1){
		textDayth = "st";
	} else if (iDay%10 == 2){
		textDayth = "nd";
	} else if (iDay%10 == 3){
		textDayth = "rd";
	} else {
		textDayth = "th";
	}
	return Form("%i_%s_%i",iDay, cMonth[iMonth-1].Data(), iYear);

}

Double_t ReturnRapidityStringAndDouble(TString cutSel, TString& rapidityRangeDummy){
	TString rapitdityCutNumberDummy = cutSel(4,1);
	if (rapitdityCutNumberDummy.CompareTo("0") == 0){
		cout << "using rapidity of 0.9" << endl;
		rapidityRangeDummy 	= "1.35";
		return 1.35*2;
	}
	if (rapitdityCutNumberDummy.CompareTo("1") == 0){
		cout << "using rapidity of 0.8" << endl;
		rapidityRangeDummy 	= "0.8";
		return 1.6;
	}
	if (rapitdityCutNumberDummy.CompareTo("2") == 0){
		cout << "using rapidity of 0.7" << endl;
		rapidityRangeDummy 	= "0.7";
		return 1.4;
	}
	if (rapitdityCutNumberDummy.CompareTo("3") == 0){
		cout << "using rapidity of 0.6" << endl;
		rapidityRangeDummy 	= "0.6";
		return 1.2;
	}
	if (rapitdityCutNumberDummy.CompareTo("4") == 0){
		cout << "using rapidity of 0.5" << endl;
		rapidityRangeDummy 	= "0.5";
		return 1.;
	}
	if (rapitdityCutNumberDummy.CompareTo("5") == 0){
		cout << "using rapidity of 0.85" << endl;
		rapidityRangeDummy 	= "0.85";
		return 0.85*2;
	}
	if (rapitdityCutNumberDummy.CompareTo("6") == 0){
		cout << "using rapidity of 0.75" << endl;
		rapidityRangeDummy 	= "0.75";
		return 0.75*2;
	}
	if (rapitdityCutNumberDummy.CompareTo("7") == 0){
		cout << "using rapidity of 0.3" << endl;
		rapidityRangeDummy   = "0.3";
		return 0.6;
	}
	if (rapitdityCutNumberDummy.CompareTo("8") == 0){
		cout << "using rapidity of 0.35" << endl;
		rapidityRangeDummy   = "0.35";
		return 0.7;
	}
	if (rapitdityCutNumberDummy.CompareTo("9") == 0){
		cout << "using rapidity of 0.4" << endl;
		rapidityRangeDummy   = "0.4";
		return 0.8;   
	} else {
		cout <<  " no rapidity Range selected" << endl;
		return 1.;
	}
}

Double_t ReturnDeltaEta(TString cutSel){
	
	TString etaCutNumber(cutSel(1,1));
	if (etaCutNumber.CompareTo("0")==0){
		cout << "using eta for gammas of 0.9" << endl;
		return  1.8;
	}
	else if (etaCutNumber.CompareTo("1")==0){
		cout << "using eta for gammas of 0.6" << endl;
		return  1.2;
	}
	else if (etaCutNumber.CompareTo("2")==0){
		cout << "using eta for gammas of 1.4" << endl;
		return  2.8;
	}
	else if (etaCutNumber.CompareTo("3")==0){
		cout << "using eta for gammas of 0.8" << endl;
		return 1.6;
	}
	else if (etaCutNumber.CompareTo("4")==0){
		cout << "using eta for gammas of 0.75" << endl;
		return  1.5;
	}
	else if (etaCutNumber.CompareTo("5")==0){
		cout << "using eta for gammas of 0.5" << endl;
		return  1.0;
	}
	else if (etaCutNumber.CompareTo("6")==0){
		cout << "using eta for gammas of 5.0" << endl;
		return  10.0;
	}
	else if (etaCutNumber.CompareTo("7")==0){
		cout << "using eta for gammas of 0.3" << endl;
		return  0.6;
	}
	else if (etaCutNumber.CompareTo("8")==0){
		cout << "using eta for gammas of 0.4" << endl;
		return 0.8;
	}
	
	cout << "Eta Value NOT found!!! using eta for gammas of 0.9" << endl;
	return 1.8;	
}

void ReturnSeparatedCutNumber(TString cutSel, TString& gammaCutNumber,TString& electronCutNumber,TString& mesonCutNumber, Bool_t kDalitz=kFALSE){
	
	TObjArray *arr;
	arr = cutSel.Tokenize("_");
	TObjString* objstrGamma;
	TObjString* objstrElectron;
	TObjString* objstrMeson;
	
	if (kDalitz){
		objstrGamma = (TObjString*)arr->At(0);
		objstrElectron = (TObjString*)arr->At(1);
		objstrMeson = (TObjString*)arr->At(2);
		
		gammaCutNumber= objstrGamma->GetString();
		electronCutNumber = objstrElectron->GetString();
		mesonCutNumber = objstrMeson->GetString();
	}
	else {
		objstrGamma = (TObjString*)arr->At(0);
		objstrMeson = (TObjString*)arr->At(1);
		
		gammaCutNumber= objstrGamma->GetString();
		mesonCutNumber = objstrMeson->GetString();
	}
	cout << cutSel.Data() << "\t" << gammaCutNumber.Data() << "\t" << electronCutNumber.Data() << "\t" << mesonCutNumber.Data() << endl;
	return;
}

void ReturnSeparatedCutNumberAdvanced(TString cutSel, TString& eventCutNumber, TString& gammaCutNumber, TString& clusterCutNumber, TString& electronCutNumber,TString& mesonCutNumber, Int_t type=0){
	
	TObjArray *arr;
	arr = cutSel.Tokenize("_");
	TObjString* objstrEvent;
	TObjString* objstrGamma;
	TObjString* objstrCluster;
	TObjString* objstrElectron;
	TObjString* objstrMeson;
	
	if (type == 0){ // PCM-PCM
		objstrEvent = (TObjString*)arr->At(0);
		objstrGamma = (TObjString*)arr->At(1);
		objstrMeson = (TObjString*)arr->At(2);
		
		eventCutNumber = objstrEvent->GetString();
		gammaCutNumber = objstrGamma->GetString();
		mesonCutNumber = objstrMeson->GetString();

	} else if (type == 1){ //PCM dalitz
		objstrEvent = (TObjString*)arr->At(0);
		objstrGamma = (TObjString*)arr->At(1);
		objstrElectron = (TObjString*)arr->At(2);
		objstrMeson = (TObjString*)arr->At(3);
		
		eventCutNumber= objstrEvent->GetString();
		gammaCutNumber= objstrGamma->GetString();
		electronCutNumber = objstrElectron->GetString();
		mesonCutNumber = objstrMeson->GetString();
	} else if (type == 2){ //PCM-EMCAL
		objstrEvent = (TObjString*)arr->At(0);
		objstrGamma = (TObjString*)arr->At(1);
		objstrCluster = (TObjString*)arr->At(2);
		objstrMeson = (TObjString*)arr->At(3);
		
		eventCutNumber= objstrEvent->GetString();
		gammaCutNumber= objstrGamma->GetString();
		clusterCutNumber = objstrCluster->GetString();
		mesonCutNumber = objstrMeson->GetString();
	} else if (type == 3){ //PCM-PHOS
		objstrEvent = (TObjString*)arr->At(0);
		objstrGamma = (TObjString*)arr->At(1);
		objstrCluster = (TObjString*)arr->At(2);
		objstrMeson = (TObjString*)arr->At(3);
		
		eventCutNumber= objstrEvent->GetString();
		gammaCutNumber= objstrGamma->GetString();
		clusterCutNumber = objstrCluster->GetString();
		mesonCutNumber = objstrMeson->GetString();
	}  else if (type == 4){ //EMCAL-EMCAL
		objstrEvent = (TObjString*)arr->At(0);
		objstrCluster = (TObjString*)arr->At(1);
		objstrMeson = (TObjString*)arr->At(2);
		
		eventCutNumber= objstrEvent->GetString();
		clusterCutNumber = objstrCluster->GetString();
		mesonCutNumber = objstrMeson->GetString();
	}  else if (type == 5){ //PHOS-PHOS
		objstrEvent = (TObjString*)arr->At(0);
		objstrCluster = (TObjString*)arr->At(1);
		objstrMeson = (TObjString*)arr->At(2);
		
		eventCutNumber= objstrEvent->GetString();
		clusterCutNumber = objstrCluster->GetString();
		mesonCutNumber = objstrMeson->GetString();
	}
	
	cout << cutSel.Data() << "\t" << eventCutNumber.Data() << "\t" << gammaCutNumber.Data() << "\t" <<  clusterCutNumber.Data() << "\t" <<electronCutNumber.Data() << "\t" << mesonCutNumber.Data() << endl;
	return;
}



Int_t ReturnSeparatedCutNumberPiPlPiMiPiZero(TString cutSel, TString& typeCutNumber, TString& eventCutNumber, TString& gammaCutNumber, TString& clusterCutNumber, TString& pionCutNumber, TString& neutralPionCutNumber, TString& mesonCutNumber){
	
	TObjArray *arr;
	arr = cutSel.Tokenize("_");
	TObjString* objstrType;
	TObjString* objstrEvent;
	TObjString* objstrGamma;
	TObjString* objstrCluster;
	TObjString* objstrPion;
	TObjString* objstrNeutralPion;
	TObjString* objstrMeson;
	Int_t mode = -1;
	objstrType  = (TObjString*)arr->At(0);
	typeCutNumber  = objstrType->GetString();
	
	if (typeCutNumber.CompareTo("0") == 0){ // PCM-PCM
		objstrType  = (TObjString*)arr->At(0);
		objstrEvent = (TObjString*)arr->At(1);
		objstrGamma = (TObjString*)arr->At(2);
		objstrPion = (TObjString*)arr->At(3);
		objstrNeutralPion = (TObjString*)arr->At(4);
		objstrMeson = (TObjString*)arr->At(5);
		
		typeCutNumber  = objstrType->GetString();
		eventCutNumber= objstrEvent->GetString();
		gammaCutNumber= objstrGamma->GetString();
		pionCutNumber = objstrPion->GetString();
		neutralPionCutNumber = objstrNeutralPion->GetString();
		mesonCutNumber = objstrMeson->GetString();
		mode = 0;
	} else { 
	if (typeCutNumber.CompareTo("1") == 0)//PCM-calo
	{
		objstrType  = (TObjString*)arr->At(0);
		objstrEvent = (TObjString*)arr->At(1);
		objstrGamma = (TObjString*)arr->At(2);
		objstrCluster = (TObjString*)arr->At(3);
		objstrPion = (TObjString*)arr->At(4);
		objstrNeutralPion = (TObjString*)arr->At(5);
		objstrMeson = (TObjString*)arr->At(6);
		
		typeCutNumber  = objstrType->GetString();
		eventCutNumber= objstrEvent->GetString();
		gammaCutNumber= objstrGamma->GetString();
		clusterCutNumber = objstrCluster->GetString();
		pionCutNumber = objstrPion->GetString();
		neutralPionCutNumber = objstrNeutralPion->GetString();
		mesonCutNumber = objstrMeson->GetString();
		
		TString firstLetter(clusterCutNumber(0,1));
		if (firstLetter.CompareTo("1") == 0) mode = 2; else mode = 3;
	}
	else
	if (typeCutNumber.CompareTo("2") == 0)//calo-calo
	{
		objstrType  = (TObjString*)arr->At(0);
		objstrEvent = (TObjString*)arr->At(1);
		objstrCluster = (TObjString*)arr->At(2);
		objstrPion = (TObjString*)arr->At(3);
		objstrNeutralPion = (TObjString*)arr->At(4);
		objstrMeson = (TObjString*)arr->At(5);
		
		typeCutNumber  = objstrType->GetString();
		eventCutNumber= objstrEvent->GetString();
		clusterCutNumber = objstrCluster->GetString();
		pionCutNumber = objstrPion->GetString();
		neutralPionCutNumber = objstrNeutralPion->GetString();
		mesonCutNumber = objstrMeson->GetString();
		
		TString firstLetter(clusterCutNumber(0,1));
		if (firstLetter.CompareTo("1") == 0) mode = 4; else mode = 5;
	}
		
	} 
	
	cout << cutSel.Data() << "\t" << typeCutNumber.Data() << "\t" << eventCutNumber.Data() << "\t" << gammaCutNumber.Data() << "\t" <<  clusterCutNumber.Data() << "\t" <<pionCutNumber.Data() << "\t" <<neutralPionCutNumber.Data() << "\t" << mesonCutNumber.Data() << endl;
	return mode;
}


Float_t ReturnBackgroundMult(TString cutSel){
	TString fBackgroundMultCutNumberDummy = cutSel(2,1);
	if (fBackgroundMultCutNumberDummy.CompareTo("0") == 0){
		cout << "using number of events for BG 5" << endl;
		return 5;
	}
	if (fBackgroundMultCutNumberDummy.CompareTo("1") == 0){
		cout << "using number of events for BG 10" << endl;
		return 10;
	}
	if (fBackgroundMultCutNumberDummy.CompareTo("2") == 0){
		cout << "using number of events for BG 15" << endl;
		return 15;
	}
	if (fBackgroundMultCutNumberDummy.CompareTo("3") == 0){
		cout << "using number of events for BG 20" << endl;
		return 20;
	}
	if (fBackgroundMultCutNumberDummy.CompareTo("4") == 0){
		cout << "using number of events for BG 2" << endl;
		return 2;
	}
	if (fBackgroundMultCutNumberDummy.CompareTo("5") == 0){
		cout << "using number of events for BG 50" << endl;
		return 50;
	}
	if (fBackgroundMultCutNumberDummy.CompareTo("6") == 0){
		cout << "using number of events for BG 80" << endl;
		return 80;
	}
	if (fBackgroundMultCutNumberDummy.CompareTo("7") == 0){
		cout << "using number of events for BG 100" << endl;
		return 100;
	}
	return 0;
}   

Double_t GetNCollFromCutNumber (TString cutNumber){
	TString systemCutNumber = cutNumber(0,1);
	TString centralityCutNumber = cutNumber(1,2);
	if (systemCutNumber.CompareTo("1") == 0 || systemCutNumber.CompareTo("2") == 0 || systemCutNumber.CompareTo("5") == 0){
		if (centralityCutNumber.CompareTo("02") == 0){ //0-20%
			return ncoll0020;
		} else if (centralityCutNumber.CompareTo("01") == 0){ //0-10%
			return ncoll0010;
		} else if (centralityCutNumber.CompareTo("12") == 0){ //10-20%
			return ncoll1020;
		} else if (centralityCutNumber.CompareTo("24") == 0){ //20-40%
			return ncoll2040;
		} else if (centralityCutNumber.CompareTo("46") == 0){ //40-60%
			return ncoll4060;
		} else if (centralityCutNumber.CompareTo("45") == 0){ //40-50%
			return ncoll4050;   
		} else if (centralityCutNumber.CompareTo("56") == 0){ //50-60%
			return ncoll5060;   
		} else if (centralityCutNumber.CompareTo("68") == 0){ //60-80%
			return ncoll6080;
		} else if (centralityCutNumber.CompareTo("67") == 0){ //60-70%
			return ncoll6070;   
		} else if (centralityCutNumber.CompareTo("78") == 0){ //70-80%
			return ncoll7080;   
		} else if (centralityCutNumber.CompareTo("89") == 0){ //80-90%
			return ncoll8090;   
		} else if (centralityCutNumber.CompareTo("48") == 0){ //40-80%
			return 77.1;   
		} else if (centralityCutNumber.CompareTo("04") == 0){ //0-40%
			return 740.;   
		}      
	} else if (systemCutNumber.CompareTo("3") == 0 || systemCutNumber.CompareTo("6") == 0){   
		if (centralityCutNumber.CompareTo("01") == 0){ //0-5%
			return ncoll0005;
		} else if (centralityCutNumber.CompareTo("12") == 0){ //0-5%
			return ncoll0510;
		} 
	} else if (systemCutNumber.CompareTo("4") == 0 || systemCutNumber.CompareTo("7") == 0){
		if (centralityCutNumber.CompareTo("69") == 0){ //75-90%
			return ncoll7590;
		}
	} else return 1.;
	return 1.;
}

Double_t GetScalingFactorSecCorrection(TString cutNumber){
	
	TString systemCutNumber = cutNumber(0,1);
	TString centralityCutNumber = cutNumber(1,2);
	if (systemCutNumber.CompareTo("1") == 0 || systemCutNumber.CompareTo("2") == 0 || systemCutNumber.CompareTo("5") == 0){
		if (centralityCutNumber.CompareTo("02") == 0){ //0-20%
			return 1./0.302 -1.; 
		} else if (centralityCutNumber.CompareTo("01") == 0){ //0-10%
			return 1./0.303 -1.;
		} else if (centralityCutNumber.CompareTo("12") == 0){ //10-20%
			return 1./0.2989 -1.;
		} else if (centralityCutNumber.CompareTo("24") == 0){ //20-40%
			return 1./0.308 -1.;
		} else if (centralityCutNumber.CompareTo("46") == 0){ //40-60%
			return 1./0.344 -1.;
		} else if (centralityCutNumber.CompareTo("45") == 0){ //40-50%
			return 1./0.342 -1.;   
		} else if (centralityCutNumber.CompareTo("56") == 0){ //50-60%
			return 1./0.346 -1.;   
		} else if (centralityCutNumber.CompareTo("68") == 0){ //60-80%
			return 1./0.3979 -1.;
		} else if (centralityCutNumber.CompareTo("67") == 0){ //60-70%
			return 1./0.39 -1.; 
		} else if (centralityCutNumber.CompareTo("78") == 0){ //70-80%
			return 1./0.41 -1.;
		} else if (centralityCutNumber.CompareTo("89") == 0){ //80-90%
			return 1./0.42 -1.;
		} else if (centralityCutNumber.CompareTo("48") == 0){ //40-80%
			return 1./0.36 -1.;
		} else if (centralityCutNumber.CompareTo("04") == 0){ //0-40%
			return 1./0.303 -1.;
		} else if (centralityCutNumber.CompareTo("08") == 0){ //0-40%
			return 1./0.3 -1.;
		}
	} else if (systemCutNumber.CompareTo("3") == 0 || systemCutNumber.CompareTo("6") == 0){   
		if (centralityCutNumber.CompareTo("01") == 0){ //0-5%
			return 1./0.306  -1.;
		} else if (centralityCutNumber.CompareTo("12") == 0){ //0-5%
			return 1./0.301  -1.;
		} 
	} else if (systemCutNumber.CompareTo("4") == 0 || systemCutNumber.CompareTo("7") == 0){
		if (centralityCutNumber.CompareTo("69") == 0){ //75-90%
		return 1./0.42 -1.;
		}
	} else if (systemCutNumber.CompareTo("8") == 0 || systemCutNumber.CompareTo("9") == 0){
		if (centralityCutNumber.CompareTo("00") == 0){ //MB
			return 1./0.52 -1.; 
		} else if (centralityCutNumber.CompareTo("02") == 0){ //00-20%
			return 1./0.24 -1.;
		} else if (centralityCutNumber.CompareTo("24") == 0){ //20-40%
			return 1./0.39 -1.;
		} else if (centralityCutNumber.CompareTo("46") == 0){ //40-60%
			return 1./0.61 -1.;
		} else if (centralityCutNumber.CompareTo("68") == 0){ //60-80%
			return 1./1.11 -1.;
		} else if (centralityCutNumber.CompareTo("60") == 0){ //60-100%
			return 1./1.62 -1.; 
		} else if (centralityCutNumber.CompareTo("80") == 0){ //80-100%
			return 1./3.00 -1.;
		} 
	} else return 0.;
	return 0.;
}

Double_t GetNCollFromName (TString name){
	if (name.CompareTo("0020") == 0){ //0-20%
		return ncoll0020;
	} else if (name.CompareTo("0005") == 0){ //0-5%
		return ncoll0005;
	} else if (name.CompareTo("0510") == 0){ //0-5%
		return ncoll0510;	
	} else if (name.CompareTo("0010") == 0){ //0-10%
		return ncoll0010;
	} else if (name.CompareTo("1020") == 0){ //10-20%
		return ncoll1020;	
	} else if (name.CompareTo("2040") == 0){ //20-40%
		return ncoll2040;
	} else if (name.CompareTo("4060") == 0){ //40-60%
		return ncoll4060;
	} else if (name.CompareTo("4050") == 0){ //40-50%
		return ncoll4050;
	} else if (name.CompareTo("5060") == 0){ //40-60%
		return ncoll5060;   
	} else if (name.CompareTo("6080") == 0){ //60-80%
		return ncoll6080;
	} else if (name.CompareTo("6070") == 0){ //60-80%
		return ncoll6070;   
	} else if (name.CompareTo("7080") == 0){ //60-80%
		return ncoll7080;      
	} else if (name.CompareTo("8090") == 0){ //60-80%
		return ncoll8090;
	} else if (name.CompareTo("7590") == 0){ //60-80%
		return ncoll7590;            
	} 
	else return 1.;
}

Double_t GetTAAFromName (TString name){
	if (name.CompareTo("0020") == 0){ //0-20%
		return tAA0020;
	} else if (name.CompareTo("0005") == 0){ //0-5%
		return tAA0005;
	} else if (name.CompareTo("0510") == 0){ //0-5%
		return tAA0510;
	} else if (name.CompareTo("0010") == 0){ //0-10%
		return tAA0010;
	} else if (name.CompareTo("1020") == 0){ //10-20%
		return tAA1020;
	} else if (name.CompareTo("2040") == 0){ //20-40%
		return tAA2040;
	} else if (name.CompareTo("4060") == 0){ //40-60%
		return tAA4060;
	} else if (name.CompareTo("6080") == 0){ //60-80%
		return tAA6080;
	} 
	else return 1.;
}

Double_t GetNCollErrFromCutNumber (TString cutNumber){
	TString systemCutNumber = cutNumber(0,1);
	TString centralityCutNumber = cutNumber(1,2);
	if (systemCutNumber.CompareTo("1") == 0 || systemCutNumber.CompareTo("2") == 0 || systemCutNumber.CompareTo("5") == 0){
		if (centralityCutNumber.CompareTo("02") == 0){ //0-20%
			return nCollErr0020;
		} else if (centralityCutNumber.CompareTo("01") == 0){ //0-10%
			return nCollErr0010;
		} else if (centralityCutNumber.CompareTo("12") == 0){ //10-20%
			return nCollErr1020;
		} else if (centralityCutNumber.CompareTo("24") == 0){ //20-40%
			return nCollErr2040;
		} else if (centralityCutNumber.CompareTo("46") == 0){ //40-60%
			return nCollErr4060;
		} else if (centralityCutNumber.CompareTo("45") == 0){ //40-50%
			return nCollErr4050;   
		} else if (centralityCutNumber.CompareTo("56") == 0){ //50-60%
			return nCollErr5060;   
		} else if (centralityCutNumber.CompareTo("68") == 0){ //60-80%
			return nCollErr6080;
		} else if (centralityCutNumber.CompareTo("67") == 0){ //60-70%
			return nCollErr6070;   
		} else if (centralityCutNumber.CompareTo("78") == 0){ //70-80%
			return nCollErr7080;   
		} else if (centralityCutNumber.CompareTo("89") == 0){ //80-90%
			return nCollErr8090;   
		}      
	} else if (systemCutNumber.CompareTo("3") == 0 || systemCutNumber.CompareTo("6") == 0){   
		if (centralityCutNumber.CompareTo("01") == 0){ //0-5%
			return nCollErr0005;
		} else if (centralityCutNumber.CompareTo("12") == 0){ //0-5%
			return nCollErr0510;
		} 
	} else if (systemCutNumber.CompareTo("4") == 0 || systemCutNumber.CompareTo("7") == 0){
		if (centralityCutNumber.CompareTo("69") == 0){ //75-90%
			return nCollErr7590;
		}
	} else return 1.;
	return 1;
}

Double_t GetNCollErrFromName (TString name){
	if (name.CompareTo("0020") == 0){ //0-20%
		return nCollErr0020;
	} else if (name.CompareTo("0005") == 0){ //0-5%
		return nCollErr0005;
	} else if (name.CompareTo("0510") == 0){ //0-5%
		return nCollErr0510;	
	} else if (name.CompareTo("0010") == 0){ //0-10%
		return nCollErr0010;
	} else if (name.CompareTo("1020") == 0){ //10-20%
		return nCollErr1020;	
	} else if (name.CompareTo("2040") == 0){ //20-40%
		return nCollErr2040;
	} else if (name.CompareTo("4060") == 0){ //40-60%
		return nCollErr4060;
	} else if (name.CompareTo("4050") == 0){ //40-60%
		return nCollErr4050;
	} else if (name.CompareTo("5060") == 0){ //40-60%
		return nCollErr5060;   
	} else if (name.CompareTo("6080") == 0){ //60-80%
		return nCollErr6080;
	} else if (name.CompareTo("6070") == 0){ //60-70%
		return nCollErr6070;
	} else if (name.CompareTo("7080") == 0){ //70-80%
		return nCollErr7080;   
	} else if (name.CompareTo("8090") == 0){ //80-90%
		return nCollErr8090;      
	} else if (name.CompareTo("7590") == 0){ //75-90%
		return nCollErr7590;      
	} 
	else return 1.;
	return 1.;
}

Double_t GetTAAErrFromName (TString name){
	if (name.CompareTo("0020") == 0){ //0-20%
		return tAAErr0020;
	} else if (name.CompareTo("0005") == 0){ //0-5%
		return tAAErr0005;
	} else if (name.CompareTo("0510") == 0){ //5-10%
		return tAAErr0510;
	} else if (name.CompareTo("0010") == 0){ //0-10%
		return tAAErr0010;
	} else if (name.CompareTo("1020") == 0){ //10-20%
		return tAAErr1020;
	} else if (name.CompareTo("2040") == 0){ //20-40%
		return tAAErr2040;
	} else if (name.CompareTo("4060") == 0){ //40-60%
		return tAAErr4060;
	} else if (name.CompareTo("6080") == 0){ //60-80%
		return tAAErr6080;
	} 
	else return 1.;
   return 1.;
}

void ReturnParameterSetFittingPbPb(TString cutsel, Double_t* parameters){
	TString centralityCutNumber = cutsel(0,3);
	cout << "bla here" << endl;
	if (centralityCutNumber.CompareTo("102") == 0 || centralityCutNumber.CompareTo("101") == 0 || centralityCutNumber.CompareTo("112") == 0  || centralityCutNumber.CompareTo("502") == 0 || centralityCutNumber.CompareTo("501") == 0 || centralityCutNumber.CompareTo("512") == 0 ){ //0-20%
		parameters[0] = 100.;
		parameters[1] = 1500.;
		parameters[2] = 4.;
		parameters[3] = 40.;
		parameters[4] = 0.07;
		parameters[5] = 0.7;
		parameters[6] = 3.;
		parameters[7] = 4.5;
		parameters[8] = 2.;
		parameters[9] = 18.;
	} else if (centralityCutNumber.CompareTo("301") == 0 || centralityCutNumber.CompareTo("601") == 0){ //0-5%
		parameters[0] = 50.;
		parameters[1] = 1500.;
		parameters[2] = 4.;
		parameters[3] = 40.;
		parameters[4] = 0.07;
		parameters[5] = 0.7;
		parameters[6] = 3.;
		parameters[7] = 6.;
		parameters[8] = 2.;
		parameters[9] = 18.;
	} else if (centralityCutNumber.CompareTo("312") == 0 || centralityCutNumber.CompareTo("612") == 0 ){ //5-10%
		parameters[0] = 50.;
		parameters[1] = 1500.;
		parameters[2] = 4.;
		parameters[3] = 40.;
		parameters[4] = 0.07;
		parameters[5] = 0.7;
		parameters[6] = 3.;
		parameters[7] = 6.;
		parameters[8] = 2.;
		parameters[9] = 18.;
	} else if (centralityCutNumber.CompareTo("124") == 0 || centralityCutNumber.CompareTo("524") == 0){ //20-40%
		parameters[0] = 50.;
		parameters[1] = 900.;
		parameters[2] = 3.;
		parameters[3] = 30.;
		parameters[4] = 0.06;
		parameters[5] = 0.5;
		parameters[6] = 3.;
		parameters[7] = 4.;
		parameters[8] = 2.;
		parameters[9] = 18.;
	} else if (centralityCutNumber.CompareTo("146") == 0 || centralityCutNumber.CompareTo("546") == 0){ //40-60%
		parameters[0] = 10.;
		parameters[1] = 200.;
		parameters[2] = 4.;
		parameters[3] = 40.;
		parameters[4] = 0.07;
		parameters[5] = 0.7;
		parameters[6] = 4.2;
		parameters[7] = 4.5;
		parameters[8] = 2.;
		parameters[9] = 18.;
	} else if (centralityCutNumber.CompareTo("168") == 0 || centralityCutNumber.CompareTo("568") == 0){ //60-80%
		parameters[0] = 1.;
		parameters[1] = 80.;
		parameters[2] = 2.5;
		parameters[3] = 30.;
		parameters[4] = 0.05;
		parameters[5] = 0.5;
		parameters[6] = 2.;
		parameters[7] = 3.3;
		parameters[8] = 2.;
		parameters[9] = 18.;
	} 	else {
		parameters[0] = 1.;
		parameters[1] = 80.;
		parameters[2] = 2.5;
		parameters[3] = 30.;
		parameters[4] = 0.05;
		parameters[5] = 0.5;
		parameters[6] = 2.;
		parameters[7] = 3.3;
		parameters[8] = 2.;
		parameters[9] = 18.;
	} 
	return;
}

void ReturnParameterSetFittingPbPbFromString(TString centralityCutNumber, Double_t* parameters){
	if (centralityCutNumber.CompareTo("0020") == 0 || centralityCutNumber.CompareTo("0010") == 0 || centralityCutNumber.CompareTo("1020") == 0  ){ //0-20%
		parameters[0] = 100.;
		parameters[1] = 1500.;
		parameters[2] = 4.;
		parameters[3] = 40.;
		parameters[4] = 0.07;
		parameters[5] = 0.7;
		parameters[6] = 3.;
		parameters[7] = 4.5;
		parameters[8] = 2.;
		parameters[9] = 18.;
		parameters[10] = 150.;
		parameters[11] = 11.5;
		parameters[12] = 0.135;
		parameters[13] = 4.3;
		parameters[14] = 7.6;
	} else if (centralityCutNumber.CompareTo("0005") == 0){ //0-5%
		parameters[0] = 50.;
		parameters[1] = 1500.;
		parameters[2] = 4.;
		parameters[3] = 40.;
		parameters[4] = 0.07;
		parameters[5] = 0.7;
		parameters[6] = 3.;
		parameters[7] = 6.;
		parameters[8] = 2.;
		parameters[9] = 18.;
		parameters[10] = 150.;
		parameters[11] = 11.5;
		parameters[12] = 0.135;
		parameters[13] = 4.3;
		parameters[14] = 7.6;
	} else if (centralityCutNumber.CompareTo("0510") == 0 ){ //5-10%
		parameters[0] = 50.;
		parameters[1] = 1500.;
		parameters[2] = 4.;
		parameters[3] = 40.;
		parameters[4] = 0.07;
		parameters[5] = 0.7;
		parameters[6] = 3.;
		parameters[7] = 6.;
		parameters[8] = 2.;
		parameters[9] = 18.;
		parameters[10] = 150.;
		parameters[11] = 11.5;
		parameters[12] = 0.135;
		parameters[13] = 4.3;
		parameters[14] = 7.6;
	} else if (centralityCutNumber.CompareTo("2040") == 0){ //20-40%
		parameters[0] = 50.;
		parameters[1] = 900.;
		parameters[2] = 3.;
		parameters[3] = 30.;
		parameters[4] = 0.06;
		parameters[5] = 0.5;
		parameters[6] = 3.;
		parameters[7] = 4.;
		parameters[8] = 2.;
		parameters[9] = 18.;
		parameters[10] = 90.;
		parameters[11] = 10.;
		parameters[12] = 0.135;
		parameters[13] = 3.5;
		parameters[14] = 7.7;
	} else if (centralityCutNumber.CompareTo("4060") == 0){ //40-60%
		parameters[0] = 10.;
		parameters[1] = 200.;
		parameters[2] = 4.;
		parameters[3] = 40.;
		parameters[4] = 0.07;
		parameters[5] = 0.7;
		parameters[6] = 4.2;
		parameters[7] = 4.5;
		parameters[8] = 2.;
		parameters[9] = 18.;
		parameters[10] = 27.;
		parameters[11] = 10.;
		parameters[12] = 0.135;
		parameters[13] = 4.4;
		parameters[14] = 7.5;
	} else if (centralityCutNumber.CompareTo("6080") == 0){ //60-80%
		parameters[0] = 1.;
		parameters[1] = 80.;
		parameters[2] = 2.5;
		parameters[3] = 30.;
		parameters[4] = 0.05;
		parameters[5] = 0.5;
		parameters[6] = 2.;
		parameters[7] = 3.3;
		parameters[8] = 2.;
		parameters[9] = 18.;
		parameters[10] = 10.;
		parameters[11] = 10.;
		parameters[12] = 0.135;
		parameters[13] = 3.;
		parameters[14] = 7.;
		cout << "bla here" << endl;
	} 	else {
// 		Double_t parameter6080[10] = {10.,80.,2.5,30.,0.05,.5,3.,3.3,2.,18.};
		parameters[0] = 1.;
		parameters[1] = 80.;
		parameters[2] = 2.5;
		parameters[3] = 30.;
		parameters[4] = 0.05;
		parameters[5] = 0.5;
		parameters[6] = 2.;
		parameters[7] = 3.3;
		parameters[8] = 2.;
		parameters[9] = 18.;
		parameters[10] = 10.;
		parameters[11] = 10.;
		parameters[12] = 0.135;
		parameters[13] = 3.;
		parameters[14] = 7.;
	} 
	return;
}

TString GetCentralityString(TString cutNumber){
	TString centralityCutNumberStart = cutNumber(1,1);
	TString centralityCutNumberEnd = cutNumber(2,1);
	TString ppCutNumber = cutNumber(0,1);
	if (ppCutNumber.CompareTo("0") ==0){
		return "pp"; 
	} else if ( ppCutNumber.CompareTo("1") ==0 || ppCutNumber.CompareTo("2") ==0 || ppCutNumber.CompareTo("5") ==0 || ppCutNumber.CompareTo("8") ==0 || ppCutNumber.CompareTo("9") ==0){       
		if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
			return "0-100%"; 
		} else if( centralityCutNumberEnd.CompareTo("0")!=0){
			return Form("%i-%i%s", centralityCutNumberStart.Atoi()*10,centralityCutNumberEnd.Atoi()*10,"%");
		} else {
			return Form("%i-%i%s", centralityCutNumberStart.Atoi()*10,centralityCutNumberEnd.Atoi()*10,"%");
		}
	} else if (ppCutNumber.CompareTo("3") ==0 || ppCutNumber.CompareTo("6") ==0){
		if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
			return "0-45%"; 
		} else {
			return Form("%i-%i%s", centralityCutNumberStart.Atoi()*5,centralityCutNumberEnd.Atoi()*5,"%");
		}
	} else if (ppCutNumber.CompareTo("4") ==0 || ppCutNumber.CompareTo("7") ==0){
		if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
			return "45-95%"; 
		} else {
			return Form("%i-%i%s",45+centralityCutNumberStart.Atoi()*5,45+centralityCutNumberEnd.Atoi()*5,"%");
		}
	} else return "";	
}	

TString GetCentralityStringWoPer(TString cutNumber){
	TString centralityCutNumberStart = cutNumber(1,1);
	TString centralityCutNumberEnd = cutNumber(2,1);
	TString ppCutNumber = cutNumber(0,1);
	if (ppCutNumber.CompareTo("0") ==0){
		return "pp"; 
	} else if ( ppCutNumber.CompareTo("1") ==0 || ppCutNumber.CompareTo("2") ==0 || ppCutNumber.CompareTo("5") ==0 || ppCutNumber.CompareTo("8") ==0 || ppCutNumber.CompareTo("9") ==0){       
		if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
			return "0-100"; 
		} else {
			return Form("%i-%i", centralityCutNumberStart.Atoi()*10,centralityCutNumberEnd.Atoi()*10);
		}
	} else if (ppCutNumber.CompareTo("3") ==0 || ppCutNumber.CompareTo("6") ==0){
		if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
			return "0-45"; 
		} else {
			return Form("%i-%i", centralityCutNumberStart.Atoi()*5,centralityCutNumberEnd.Atoi()*5);
		}
	} else if (ppCutNumber.CompareTo("4") ==0 || ppCutNumber.CompareTo("7") ==0){
		if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
			return "45-95"; 
		} else {
			return Form("%i-%i",45+centralityCutNumberStart.Atoi()*5,45+centralityCutNumberEnd.Atoi()*5);
		}
	} else return ""; 
}  


void CalcRaa( 	TGraphAsymmErrors* graphPPSpectrum, TGraphAsymmErrors* graphPPSpectrumSystNoMat, TGraphAsymmErrors* graphPPCombinedSpectrum, TF1* fitPP,
					TGraphAsymmErrors* graphPbPbSpectrum, TGraphAsymmErrors* graphPbPbSpectrumSysNoMat,  //PbPb Yields
					TGraphAsymmErrors** graphRAA, TGraphAsymmErrors** graphRAASys,
					Double_t fNcoll, Double_t fNcollError, Double_t maxPtPP = 0., Int_t nSubEnd = 0){
	//--------------------------------------- extrapolating the pp CONV yield in pT --------------------------------------------------
//  	graphPPSpectrum->Print();
// 	cout << endl;
//  	graphPbPbSpectrum->Print();
	TGraphAsymmErrors* dummyPPSpectrum =(TGraphAsymmErrors*)graphPPCombinedSpectrum->Clone("dummyPPSpectrum");
	Double_t* xBinsPPFit = dummyPPSpectrum->GetX();
	Double_t* xBinsPPErrFitLow = dummyPPSpectrum->GetEXlow();
	Double_t* xBinsPPErrFitHigh = dummyPPSpectrum->GetEXhigh();
	Int_t nBinsPPFit = dummyPPSpectrum->GetN();
	for (Int_t i = 0; i < nBinsPPFit; i++){
// 		xBinsPPErrFitLow[i] = 0.5*2*(xBinsPPErrFitLow[i]);
// 		xBinsPPErrFitHigh[i] = 0.5*2*(xBinsPPErrFitHigh[i]);
// 		xBinsPPErrFitLow[i] = 0.001*(xBinsPPFit[i]);
// 		xBinsPPErrFitHigh[i] = 0.001*(xBinsPPFit[i]);		
		cout << xBinsPPFit[i] << "\t" << xBinsPPErrFitLow[i] << "\t" <<  xBinsPPErrFitHigh[i] << endl;
	}
	dummyPPSpectrum->Print();
// 	TString fittingOptions ="QSEM0";
	TString fittingOptions ="QSEM0EX0";
	
	Double_t* xBinsPP =  graphPPSpectrum->GetX();
	Double_t* xBinsErrPP = graphPPSpectrum->GetEXlow();
	Int_t nBinsPP = graphPPSpectrum->GetN();
	Double_t* yPP =  graphPPSpectrum->GetY();
	Double_t* yErrLowPP = graphPPSpectrum->GetEYlow();
	Double_t* yErrHighPP = graphPPSpectrum->GetEYhigh();
	Double_t* yErrLowPPSys = graphPPSpectrumSystNoMat->GetEYlow();
	Double_t* yErrHighPPSys = graphPPSpectrumSystNoMat->GetEYhigh();
	
	Double_t maxPtPPForSpec;
	if (maxPtPP !=0.){
		maxPtPPForSpec = maxPtPP;
	} else {
		maxPtPPForSpec = xBinsPP[nBinsPP-1];
	}
	
	Double_t* xBinsPbPb =  graphPbPbSpectrum->GetX();
	Double_t* xBinsErrPbPb =  graphPbPbSpectrum->GetEXlow();
	Double_t* yPbPb =  graphPbPbSpectrum->GetY();
	
	Double_t* xBinsPPComb =  graphPPCombinedSpectrum->GetX();
   //Int_t nBinsPPComb =  graphPPSpectrum->GetN();
	Int_t firstBinPbPb=0; 
	Double_t decisionBoundary = 0.0000001;
	while (TMath::Abs(xBinsPP[firstBinPbPb] - xBinsPbPb[0]) >decisionBoundary){
		cout << xBinsPP[firstBinPbPb] << "\t" << xBinsPbPb[0] << endl;
		firstBinPbPb++;
	}
	Double_t fitBeginpp = xBinsPPComb[0];
// 	Double_t fitEndpp = xBinsPPComb[nBinsPPComb-1];
	if (!fitPP) fitPP = FitObject("l","fitRAARefTsallis","Pi0");
// 	fitPP->SetParameter(0,1.6);
// 	fitPP->SetParameter(1,7.3);
// 	fitPP->SetParameter(1,0.14);
	TFitResultPtr resultPP = dummyPPSpectrum->Fit(fitPP,fittingOptions.Data() , "", fitBeginpp, maxPtPPForSpec);
// 	cout << "Fit result with TFitResultPtr" << endl;
	cout << WriteParameterToFile(fitPP)<< endl;
	
	TF1* fitPPPowerlaw = FitObject("powPure","fitRAARefPowerlaw","Pi0");
	TFitResultPtr resultPPPowerlaw = dummyPPSpectrum->Fit(fitPPPowerlaw,fittingOptions.Data() , "", 3., maxPtPPForSpec);
	cout << WriteParameterToFile(fitPPPowerlaw)<< endl;
	

	
// 	TCanvas* canvasDummy6 = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
// 	DrawGammaCanvasSettings( canvasDummy6,  0.1, 0.01, 0.015, 0.08);
// 	canvasDummy6->SetLogy();
// 	canvasDummy6->SetLogx();
// 	TH2F * histo2DDummy6;
// 	histo2DDummy6 = new TH2F("histo2DDummy6","histo2DDummy5",1000,0.23,30.,1000,1e-8,10);
// 	SetStyleHistoTH2ForGraphs(histo2DDummy6, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
// 	histo2DDummy6->DrawCopy(); 
// 
// 	DrawGammaSetMarkerTGraphAsym(graphPPCombinedSpectrum, 20,1, kRed, kRed, 1, kTRUE);
// 	graphPPCombinedSpectrum->Draw("pEsame");
// 	
// 	fitPP->SetLineColor(kBlue+2);
// 	fitPP->Draw("same");
// 	
// 	canvasDummy6->Update();
// 	canvasDummy6->Print("debugSpectrum.eps");
// 
	
	TGraphAsymmErrors* graphPPSpectrumExtended = (TGraphAsymmErrors*) graphPbPbSpectrum->Clone("graphPPSpectrumExtended");
	TGraphAsymmErrors* graphPPSpectrumExtendedSys = (TGraphAsymmErrors*) graphPbPbSpectrumSysNoMat->Clone("graphPPSpectrumExtendedSys"); 
	
	(*graphRAA) = new TGraphAsymmErrors(graphPbPbSpectrum->GetN()-nSubEnd);
	(*graphRAASys) = new TGraphAsymmErrors(graphPbPbSpectrum->GetN()-nSubEnd);
	
	for (Int_t i=0;i<graphPbPbSpectrum->GetN()-nSubEnd ;i++){
// 		Double_t xCenter = 0.5*(xPtLimits[i+1]+xPtLimits[i]);
// 		cout <<  xBinsPP[i+firstBinPbPb] << "\t" << xBinsPbPb[i] << "\t" << xBinsErrPP[i+firstBinPbPb] << "\t" << xBinsErrPbPb[i] << endl;
		if (TMath::Abs(xBinsPP[i+firstBinPbPb] - xBinsPbPb[i]) < decisionBoundary && TMath::Abs(xBinsErrPP[i+firstBinPbPb] - xBinsErrPbPb[i]) < decisionBoundary && xBinsPP[i] < maxPtPP ){
			cout<< "xCenter::"<< xBinsPP[i+firstBinPbPb]<< " " << xBinsErrPP[i+firstBinPbPb] << " " << yPP[i+firstBinPbPb]<< "  "<<xBinsPbPb[i]<< " " <<xBinsErrPbPb[i] <<"  " << yPbPb[i] / fNcoll<< " " << yPbPb[i] /(fNcoll*yPP[i+firstBinPbPb])<< endl; 
			graphPPSpectrumExtended->SetPoint(i,xBinsPbPb[i],yPP[i+firstBinPbPb]);
			graphPPSpectrumExtendedSys->SetPoint(i,xBinsPbPb[i],yPP[i+firstBinPbPb]);
			graphPPSpectrumExtended->SetPointError(i, xBinsErrPbPb[i], xBinsErrPbPb[i], yErrLowPP[i+firstBinPbPb], yErrHighPP[i+firstBinPbPb]);
			graphPPSpectrumExtendedSys->SetPointError(i, xBinsErrPbPb[i], xBinsErrPbPb[i], yErrLowPPSys[i+firstBinPbPb], yErrHighPPSys[i+firstBinPbPb]);			
		} else {
			Double_t ptStart = xBinsPbPb[i] - xBinsErrPbPb[i];
			Double_t ptEnd = xBinsPbPb[i] + xBinsErrPbPb[i];
			Double_t binWidth = ptEnd-ptStart;
			Double_t yieldPP = fitPP->Integral(ptStart, ptEnd, resultPP->GetParams()) / binWidth;
			Double_t errorYieldPP = fitPP->IntegralError(ptStart, ptEnd, resultPP->GetParams(), resultPP->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
			cout<< "xCenter::"<< ptStart<< " " << ptEnd << " " << yieldPP<< "+-" << errorYieldPP <<" "  << yPbPb[i] / fNcoll<< " " << yPbPb[i] /(fNcoll*yieldPP)<< endl; 
			Double_t yieldPPPowerlaw = fitPPPowerlaw->Integral(ptStart, ptEnd, resultPPPowerlaw->GetParams()) / binWidth;
			Double_t errorYieldPPPowerlaw = fitPPPowerlaw->IntegralError(ptStart, ptEnd, resultPPPowerlaw->GetParams(), resultPPPowerlaw->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
			cout<< "xCenter::"<< ptStart<< " " << ptEnd << " " << yieldPPPowerlaw<< "+-" << errorYieldPPPowerlaw << " " << yPbPb[i] / fNcoll<< " " << yPbPb[i] /(fNcoll*yieldPP)<< endl; 
			cout << abs(yieldPPPowerlaw-yieldPP)/yieldPP *100 << endl;
			graphPPSpectrumExtended->SetPoint(i,xBinsPbPb[i],yieldPP);
			graphPPSpectrumExtended->SetPointError(i, xBinsErrPbPb[i], xBinsErrPbPb[i], errorYieldPP, errorYieldPP);
			graphPPSpectrumExtendedSys->SetPoint(i,xBinsPbPb[i],yieldPP);
			graphPPSpectrumExtendedSys->SetPointError(i, xBinsErrPbPb[i], xBinsErrPbPb[i],yieldPP* yErrLowPPSys[nBinsPP-1]/yPP[nBinsPP-1], yieldPP*yErrHighPPSys[nBinsPP-1]/yPP[nBinsPP-1]);
		}
	} 
	graphPPSpectrumExtended->Print();
// 	graphPbPbSpectrumSysNoMat->Print();
	Double_t* yPPExtended =  graphPPSpectrumExtended->GetY();
	for(Int_t i=0; i<(*graphRAA)->GetN(); i++){
		(*graphRAA)->SetPoint(i,xBinsPbPb[i],yPbPb[i] /(fNcoll* yPPExtended[i]));
		Double_t errYStat = pow( pow(graphPbPbSpectrum->GetErrorYlow(i)/yPbPb[i],2.) + 
										pow(graphPPSpectrumExtended->GetErrorYlow(i)/yPPExtended[i],2.), 0.5)*yPbPb[i] /(fNcoll* yPPExtended[i]);
		
		(*graphRAA)->SetPointError(i, xBinsErrPbPb[i], xBinsErrPbPb[i], errYStat, errYStat);
		(*graphRAASys)->SetPoint(i,xBinsPbPb[i],yPbPb[i] /(fNcoll* yPPExtended[i]));
		Double_t errYlow = pow( pow(graphPbPbSpectrumSysNoMat->GetErrorYlow(i)/yPbPb[i],2.) + 
										pow(graphPPSpectrumSystNoMat->GetErrorYlow(i)/yPPExtended[i],2.) 
							  , 0.5)*yPbPb[i] /(fNcoll* yPPExtended[i]); //+ pow(fNcollError/fNcoll,2.) fNcollError not taking into account 
		Double_t errYhigh = pow( pow(graphPbPbSpectrumSysNoMat->GetErrorYhigh(i)/yPbPb[i],2.) + 
										 pow(graphPPSpectrumSystNoMat->GetErrorYhigh(i)/yPPExtended[i],2.) 
							  , 0.5)*yPbPb[i] /(fNcoll* yPPExtended[i]); //+ pow(fNcollError/fNcoll,2.) fNcollError not taking into account 
                       
      cout << i << "\t" << yPbPb[i]<< "\t" << yPPExtended[i]<< "\t" << graphPbPbSpectrumSysNoMat->GetErrorYlow(i) <<"\t" << errYlow << "\t" << errYhigh << endl;
		(*graphRAASys)->SetPointError(i, xBinsErrPbPb[i], xBinsErrPbPb[i], errYlow, errYhigh);
	}
	Int_t b = 0;
	while (yPbPb[b] == 0){
		(*graphRAA)->RemovePoint(0);
		(*graphRAASys)->RemovePoint(0);
		b++;
	}	
	
	delete dummyPPSpectrum;
	
	if(fNcollError){}

}

Int_t GetBinning(TObject *Obj_Dummy, Double_t* doubleBinningX){
	TString ClassName = Obj_Dummy->ClassName();
	if(ClassName.BeginsWith("TH1")){
		TH1D *histo = (TH1D*)Obj_Dummy;
		Int_t bin = 0;
		for (Int_t i = 1; i < histo->GetNbinsX()+1; i++){
			if (histo->GetBinContent(i) != 0){
				doubleBinningX[bin] = histo->GetBinLowEdge(i);
				doubleBinningX[bin+1] = histo->GetXaxis()->GetBinUpEdge(i);
				bin++;
			}
		}
// 		for (Int_t i = 0; i < bin+1; i++){
// 			cout << doubleBinningX[i] << "\t," ;
// 		}
// 		cout << endl;
		return bin+1;
	} else if(ClassName.CompareTo("TGraphErrors")==0){
		TGraphErrors *graph = (TGraphErrors*)Obj_Dummy;
		Double_t* binCenter = graph->GetX();
		Double_t* binContent = graph->GetY();
		Double_t* binError = graph->GetEX();
		Int_t nBins = graph->GetN();
		Int_t bin = 0;
		for (Int_t i = 0; i < nBins; i++){
			if (binContent[i] != 0){
				doubleBinningX[bin] = binCenter[i]-binError[i];
				doubleBinningX[bin+1] = binCenter[i]+binError[i];
				bin++;
			}
		}
// 		for (Int_t i = 0; i < nBins+1; i++){
// 			cout << doubleBinningX[i] << "\t," ;
// 		}
// 		cout << endl;
		return bin+1;
	} else if(ClassName.CompareTo("TGraphAsymmErrors")==0){
		TGraphAsymmErrors *graph = (TGraphAsymmErrors*)Obj_Dummy;
		Double_t* binCenter = graph->GetX();
		Double_t* binContent = graph->GetY();
		Double_t* binErrorDown = graph->GetEXlow();
		Double_t* binErrorUp = graph->GetEXhigh();
		Int_t nBins = graph->GetN();
		Int_t bin = 0;
		for (Int_t i = 0; i < nBins; i++){
			if (binContent[i] != 0){
				doubleBinningX[bin] = binCenter[i]-binErrorDown[i];
				doubleBinningX[bin+1] = binCenter[i]+binErrorUp[i];
				bin++;
			}
		}
// 		for (Int_t i = 0; i < nBins+1; i++){
// 			cout << doubleBinningX[i] << "\t," ;
// 		}
// 		cout << endl;
		return bin+1;
	} else {
		cout << " class not defined" << endl;
		return 0;
	}
	
}

Int_t CompareBinning(Int_t nBinsA, Double_t* binningA, Int_t nBinsB, Double_t* binningB, Double_t* newBinning, Int_t* nBinsToBeCombinedA, Int_t* nBinsToBeCombinedB, TString returnStr = "", Double_t decisionBoundary = 0.0000001){
	Int_t startingBin = 0;
	//Finding startBin
	cout << "Binning A" << endl;
		for (Int_t i = 0; i < nBinsA ; i++){
			cout << binningA[i] << "\t" ;
		}
	cout << endl;
	cout << "Binning B" << endl;
		for (Int_t i = 0; i < nBinsB ; i++){
			cout << binningB[i] << "\t" ;
		}
	cout << endl;	
	if (binningA[0] < binningB[0]){
		while ((binningB[0] - binningA[startingBin])>decisionBoundary && startingBin <nBinsA){
			startingBin++;
		}
		cout << "Binning A" << endl;
		for (Int_t i = startingBin; i < nBinsA ; i++){
			cout << binningA[i] << "\t" ;
		}
		cout << endl;
		cout << "Binning B" << endl;
		for (Int_t i = 0; i < nBinsB ; i++){
			cout << binningB[i] << "\t" ;
		}
		cout << endl;
	
		Int_t c= 1;
		Int_t startBinNewBins = 1;
		Int_t binsToBeMergedB = 1;
		Int_t binsToBeMergedA = 1;
		cout << "Binning A starts earlier, combined binning will start at " << startingBin << " with " << binningA[startingBin] << endl;
		newBinning[0] = binningA[startingBin];
		for (Int_t i = startingBin+1; i < nBinsA; i++){
			if (c > nBinsB-1) return startBinNewBins;
			newBinning[startBinNewBins] = binningA[i];
			binsToBeMergedB = 1;
			binsToBeMergedA = 1;
//  			cout << "entering loop " << binningA[i] << "\t" << binningB[c] << endl;
			while( (binningA[i] - binningB[c])>decisionBoundary){
				cout << binningA[i] << "\t" << binningB[c] << endl;
				c++;
				binsToBeMergedB++; 
				if (!((binningA[i] - binningB[c])>decisionBoundary)){
//  					cout << "nächstes bin ist größer in B" << endl;
					if (TMath::Abs(binningA[i]- binningB[c]) < decisionBoundary) {
//  						cout << " alles super" << endl;
					} else {
//  						cout << " müssen einen hoch" << endl;
						i++;
						newBinning[startBinNewBins] = binningA[i];
						binsToBeMergedA++;
					}
				}
				if (c > nBinsB-1) return startBinNewBins;
			}
			nBinsToBeCombinedB[startBinNewBins-1] = binsToBeMergedB;
			nBinsToBeCombinedA[startBinNewBins-1] = binsToBeMergedA;
			if (c > nBinsB-1) return startBinNewBins;
			cout << "exiting loop " << binningA[i] << "\t" << binningB[c] << "\t bins needed to merged A :"<<nBinsToBeCombinedA[startBinNewBins-1] <<" B :"<<nBinsToBeCombinedB[startBinNewBins-1] <<endl;
			startBinNewBins++;
			c++;
		}			
		return startBinNewBins;
	} else  if (binningB[0] < binningA[0]){
		while (!(TMath::Abs(binningA[0]- binningB[startingBin]) < decisionBoundary) && startingBin< nBinsB){
			cout << "deviation to first bin in A for startBin: "<< startingBin << "\t" <<TMath::Abs(binningA[0]- binningB[startingBin]) << endl;
			startingBin++;
		} 
		Int_t check2 = 0;
		while (startingBin == nBinsB){
			cout << "Failed to evalute starting point in attempt " << check2	<< endl;
			check2++;
			startingBin=0;
			while (!(TMath::Abs(binningA[check2]- binningB[startingBin]) < decisionBoundary) && startingBin< nBinsB){
				cout << TMath::Abs(binningA[check2]- binningB[startingBin]) << endl;
				startingBin++;	
			}
		}
		cout << "Binning A" << endl;
		for (Int_t i = 0; i < nBinsA ; i++){
			cout << binningA[i] << "\t" ;
		}
		cout << endl;
		cout << "Binning B" << endl;
		for (Int_t i = startingBin; i < nBinsB ; i++){
			cout << binningB[i] << "\t" ;
		}
		cout << endl;
		cout << "Binning B starts earlier, combined binning will start at " << startingBin << " with " << binningB[startingBin] << endl;
		
		Int_t c= startingBin+1;
		Int_t startBinNewBins = 1;
		Int_t binsToBeMergedB = 1;
		Int_t binsToBeMergedA = 1;
		newBinning[0] = binningA[0];
		for (Int_t i = 1; i < nBinsA; i++){
			if (c > nBinsB-1) return startBinNewBins;
			newBinning[startBinNewBins] = binningA[i];
			binsToBeMergedB = 1;
			binsToBeMergedA = 1;
// 			cout << "entering loop " << binningA[i] << "\t" << binningB[c] << endl;
			while( (binningA[i] - binningB[c])>decisionBoundary){
// 				cout << binningA[i] << "\t" << binningB[c] << endl;
				if (!((binningA[i] - binningB[c+1])>decisionBoundary)){
// 					cout << "nächstes bin ist größer in B" << endl;
					if (TMath::Abs(binningA[i]- binningB[c+1]) < decisionBoundary) {
// 						cout << " alles super" << endl;
					} else {
// 						cout << " müssen einen hoch" << endl;
						i++;
						newBinning[startBinNewBins] = binningA[i];
						binsToBeMergedA++;
					}
				}
				c++;
				binsToBeMergedB++; 
				if (c > nBinsB-1) return startBinNewBins-1;
			}
			nBinsToBeCombinedB[startBinNewBins-1] = binsToBeMergedB;
			nBinsToBeCombinedA[startBinNewBins-1] = binsToBeMergedA;
 			cout << "exiting loop " << binningA[i] << "\t" << binningB[c] << "\t bins needed to merged A :"<<nBinsToBeCombinedA[startBinNewBins-1] <<" B :"<<nBinsToBeCombinedB[startBinNewBins-1] <<endl;
			startBinNewBins++;
			c++;
		}			
		return startBinNewBins;

	} else {
		cout << "Binning A" << endl;
		for (Int_t i = 0; i < nBinsA ; i++){
			cout << binningA[i] << "\t" ;
		}
		cout << endl;
		cout << "Binning B" << endl;
		for (Int_t i = startingBin; i < nBinsB ; i++){
			cout << binningB[i] << "\t" ;
		}
		cout << endl;

		cout << "Both start at the same value " << binningA[0] << endl;
		
		Int_t c= 1;
		Int_t startBinNewBins = 1;
		Int_t binsToBeMergedB = 1;
		Int_t binsToBeMergedA = 1;
		newBinning[0] = binningA[0];
		for (Int_t i = 1; i < nBinsA; i++){
			if (c > nBinsB-1) return startBinNewBins-1;
			newBinning[startBinNewBins] = binningA[i];
			binsToBeMergedB = 1;
			binsToBeMergedA = 1;
// 			cout << "entering loop " << binningA[i] << "\t" << binningB[c] << endl;
			while( (binningA[i] - binningB[c])>decisionBoundary){
// 				cout << binningA[i] << "\t" << binningB[c] << endl;
				if (!((binningA[i] - binningB[c+1])>decisionBoundary)){
// 					cout << "nächstes bin ist größer in B" << endl;
					if (TMath::Abs(binningA[i]- binningB[c+1]) < decisionBoundary) {
// 						cout << " alles super" << endl;
					} else {
// 						cout << " müssen einen hoch" << endl;
						i++;
						newBinning[startBinNewBins] = binningA[i];
						binsToBeMergedA++;
					}
				}
				c++;
				binsToBeMergedB++; 
				if (c > nBinsB-1) return startBinNewBins-1;
			}
			nBinsToBeCombinedB[startBinNewBins-1] = binsToBeMergedB;
			nBinsToBeCombinedA[startBinNewBins-1] = binsToBeMergedA;
// 			cout << "exiting loop " << binningA[i] << "\t" << binningB[c] << "\t bins needed to merged A :"<<nBinsToBeCombinedA[startBinNewBins-1] <<" B :"<<nBinsToBeCombinedB[startBinNewBins-1] <<endl;
			startBinNewBins++;
			c++;
		}			
		return startBinNewBins;

	}
	
	
	return 0;

   if(returnStr){}


}

Int_t FindFirstCommonBin(Double_t* vectorNewBinning, Double_t* oldBinning, Double_t decisionBoundary = 0.0000001){
	Int_t startingBin = 0;
	//Finding startBin
	
	while (!(TMath::Abs(vectorNewBinning[0] - oldBinning[startingBin])<decisionBoundary)){ 
		startingBin++;
		cout << startingBin << "\t" <<!(TMath::Abs(vectorNewBinning[0] - oldBinning[startingBin])<decisionBoundary) << endl;
	}
	return startingBin;
}

void RebinObjects(TObject* Obj_DummyStat,TObject* Obj_DummySyst, Double_t* vectorNewBinning, Int_t* vectorRebinFactors, Int_t nCommonBins, AliConvDataObject* outputObject, TString ClassNameStat, TString ClassNameSyst ,Bool_t scaleByBinCenter=kFALSE){ //
	Double_t binningOldStat[200] ;
	Double_t binningOldSyst[200] ;
	Int_t validBinsOldStat= GetBinning(Obj_DummyStat, binningOldStat);
	Int_t validBinsOldSys= GetBinning(Obj_DummySyst, binningOldSyst);  
   
	Int_t firstBinStat = FindFirstCommonBin(vectorNewBinning, binningOldStat);
	Int_t firstBinSyst = FindFirstCommonBin(vectorNewBinning, binningOldSyst);
	
	cout << "FirstBin stat " << firstBinStat<< "\t" << binningOldStat[firstBinStat]	<< endl;
	cout << "FirstBin sys " << firstBinSyst	<< "\t" << binningOldSyst[firstBinSyst] << endl;
// 
	cout << "statistical Errors" <<endl;
	if(ClassNameStat.BeginsWith("TH1")){
		TH1D *histo = (TH1D*)Obj_DummyStat;
		Int_t indBin = 0;
		
		while (vectorNewBinning[0] > histo->GetBinLowEdge(indBin))	{ indBin++;}
		for (Int_t commonBin = 0; commonBin < nCommonBins-1; commonBin++){
			outputObject[commonBin].valueY= 0;
			outputObject[commonBin].errorYStatLow = 0;
			outputObject[commonBin].valueX = (vectorNewBinning[commonBin] + vectorNewBinning[commonBin+1])/2;
			outputObject[commonBin].errorXLow = (vectorNewBinning[commonBin] + vectorNewBinning[commonBin+1])/2 -vectorNewBinning[commonBin];
			outputObject[commonBin].errorXHigh = (vectorNewBinning[commonBin] + vectorNewBinning[commonBin+1])/2 -vectorNewBinning[commonBin];
			for (Int_t add = 1; add < vectorRebinFactors[commonBin]+1;add++){
				if (scaleByBinCenter){
					cout << indBin << "\t" << histo->GetBinCenter(indBin) << "\t"<< (histo->GetBinContent(indBin)*histo->GetBinCenter(indBin)*histo->GetBinWidth(indBin)*2*TMath::Pi()) << endl;
					outputObject[commonBin].valueY = outputObject[commonBin].valueY+(histo->GetBinContent(indBin)*histo->GetBinCenter(indBin)*histo->GetBinWidth(indBin)*2*TMath::Pi());
					outputObject[commonBin].errorYStatLow = outputObject[commonBin].errorYStatLow+TMath::Power((histo->GetBinError(indBin)*histo->GetBinCenter(indBin)*histo->GetBinWidth(indBin)*2*TMath::Pi()),2);
				}	else {
					cout << indBin << "\t" << histo->GetBinCenter(indBin) << "\t"<<(histo->GetBinContent(indBin)*histo->GetBinWidth(indBin)) << endl;
					outputObject[commonBin].valueY = outputObject[commonBin].valueY+(histo->GetBinContent(indBin)*histo->GetBinWidth(indBin));
					outputObject[commonBin].errorYStatLow = outputObject[commonBin].errorYStatLow+TMath::Power((histo->GetBinContent(indBin)*histo->GetBinWidth(indBin)),2);
				}
				indBin++;
			}
			outputObject[commonBin].errorYStatLow = TMath::Sqrt(outputObject[commonBin].errorYStatLow);
			outputObject[commonBin].errorYStatHigh = outputObject[commonBin].errorYStatLow;
			cout << commonBin<< "\t"  <<outputObject[commonBin].valueX << "\t +- " << outputObject[commonBin].errorXLow << " \t " << outputObject[commonBin].valueY<< "\t+" << outputObject[commonBin].errorYStatHigh<<"\t-"<< outputObject[commonBin].errorYStatLow<<endl;
		}
	} else if(ClassNameStat.CompareTo("TGraphErrors")==0){
		TGraphErrors *graph = (TGraphErrors*)Obj_DummyStat;
		TGraphErrors *graphCopy = (TGraphErrors*)graph->Clone("GraphCopy");
		Double_t* valueX = graphCopy->GetX();
		Double_t* valueY = graphCopy->GetY();
		Double_t* errorX = graphCopy->GetEX();
		Double_t* errorY = graphCopy->GetEY();
		for (Int_t i = 0; i < graphCopy->GetN();i++){
			if (scaleByBinCenter){
				valueY[i] = valueY[i]*valueX[i]*(errorX[i]+errorX[i])*2*TMath::Pi();
				errorY[i] = errorY[i]*valueX[i]*(errorX[i]+errorX[i])*2*TMath::Pi();
			}	else {
				errorY[i] = errorY[i]*(errorX[i]+errorX[i]);
			}
		}
		Int_t indBin = firstBinStat;
		for (Int_t commonBin = 0; commonBin < nCommonBins-1; commonBin++){
			outputObject[commonBin].valueY= 0;
			outputObject[commonBin].errorYStatLow = 0;
			outputObject[commonBin].valueX = (vectorNewBinning[commonBin] + vectorNewBinning[commonBin+1])/2;
			outputObject[commonBin].errorXLow = (vectorNewBinning[commonBin] + vectorNewBinning[commonBin+1])/2 -vectorNewBinning[commonBin];
			outputObject[commonBin].errorXHigh = (vectorNewBinning[commonBin] + vectorNewBinning[commonBin+1])/2 -vectorNewBinning[commonBin];
			for (Int_t add = 1; add < vectorRebinFactors[commonBin]+1;add++){
				cout << indBin << "\t" << valueY[indBin] << endl;
				outputObject[commonBin].valueY = outputObject[commonBin].valueY+valueY[indBin];
				outputObject[commonBin].errorYStatLow = outputObject[commonBin].errorYStatLow+TMath::Power(errorY[indBin],2);
				indBin++;
			}
			outputObject[commonBin].errorYStatLow = TMath::Sqrt(outputObject[commonBin].errorYStatLow);
			outputObject[commonBin].errorYStatHigh = outputObject[commonBin].errorYStatLow;
			cout << commonBin<< "\t"  <<outputObject[commonBin].valueX << "\t +- " << outputObject[commonBin].errorXLow << " \t " << outputObject[commonBin].valueY<< "\t+" << outputObject[commonBin].errorYStatHigh<<"\t-"<< outputObject[commonBin].errorYStatLow<<endl;
		}
	} else if(ClassNameStat.CompareTo("TGraphAsymmErrors")==0){
		TGraphAsymmErrors *graph = (TGraphAsymmErrors*)Obj_DummyStat;
		TGraphAsymmErrors *graphCopy = (TGraphAsymmErrors*)graph->Clone("GraphCopy");
		Double_t* valueX = graphCopy->GetX();
		Double_t* valueY = graphCopy->GetY();
		Double_t* errorXlow = graphCopy->GetEXlow();
		Double_t* errorXhigh = graphCopy->GetEXhigh();
		Double_t* errorYlow = graphCopy->GetEYlow();
		Double_t* errorYhigh = graphCopy->GetEYhigh();
		for (Int_t i = 0; i < graphCopy->GetN();i++){
			if (scaleByBinCenter){
				valueY[i] = valueY[i]*valueX[i]*(errorXlow[i]+errorXhigh[i])*2*TMath::Pi();
				errorYlow[i] = errorYlow[i]*valueX[i]*(errorXlow[i]+errorXhigh[i])*2*TMath::Pi();
				errorYhigh[i] = errorYhigh[i]*valueX[i]*(errorXlow[i]+errorXhigh[i])*2*TMath::Pi();
			}	else {
				valueY[i] = valueY[i]*(errorXlow[i]+errorXhigh[i]);
				errorYlow[i] = errorYlow[i]*(errorXlow[i]+errorXhigh[i]);
				errorYhigh[i] = errorYhigh[i]*(errorXlow[i]+errorXhigh[i]);
			}
		}
		cout<< "after modification" << endl;
 		graphCopy->Print();
		Int_t indBin = firstBinStat;
		for (Int_t commonBin = 0; commonBin < nCommonBins-1; commonBin++){
			outputObject[commonBin].valueY = 0;
			outputObject[commonBin].errorYStatLow = 0;
			outputObject[commonBin].errorYStatHigh = 0;
			outputObject[commonBin].valueX = (vectorNewBinning[commonBin] + vectorNewBinning[commonBin+1])/2;
			outputObject[commonBin].errorXLow = (vectorNewBinning[commonBin] + vectorNewBinning[commonBin+1])/2 -vectorNewBinning[commonBin];
			outputObject[commonBin].errorXHigh = (vectorNewBinning[commonBin] + vectorNewBinning[commonBin+1])/2 -vectorNewBinning[commonBin];

			for (Int_t add = 1; add < vectorRebinFactors[commonBin]+1;add++){
				cout << indBin << "\t" << valueY[indBin] << endl;
				outputObject[commonBin].valueY = outputObject[commonBin].valueY+valueY[indBin];
				outputObject[commonBin].errorYStatHigh = outputObject[commonBin].errorYStatHigh + TMath::Power(errorYhigh[indBin],2);
				outputObject[commonBin].errorYStatLow = outputObject[commonBin].errorYStatLow + TMath::Power(errorYlow[indBin],2);
				indBin++;
			}
			outputObject[commonBin].errorYStatLow = TMath::Sqrt(outputObject[commonBin].errorYStatLow);
			outputObject[commonBin].errorYStatHigh = TMath::Sqrt(outputObject[commonBin].errorYStatHigh);
			cout << commonBin<< "\t"  <<outputObject[commonBin].valueX << "\t +- " << outputObject[commonBin].errorXLow << " \t " << outputObject[commonBin].valueY<< "\t+" << outputObject[commonBin].errorYStatHigh<<"\t-"<< outputObject[commonBin].errorYStatLow<<endl;
		}
	} else {
		cout << " class for Stat not defined" << endl;
		return;
	}
	cout << "systematic errors" << endl;
	if(ClassNameSyst.BeginsWith("TH1")){
		cout << "\n \n TH1" << endl << endl;
		TH1D *histo = (TH1D*)Obj_DummySyst;
		Int_t indBin = 0;
		while (vectorNewBinning[0] > histo->GetBinLowEdge(indBin))	{ indBin++;}
		for (Int_t commonBin = 0; commonBin < nCommonBins-1; commonBin++){
			outputObject[commonBin].errorYSystLow = 0;
			for (Int_t add = 1; add < vectorRebinFactors[commonBin]+1;add++){
// 				if (scaleByBinCenter){
// 					cout << indBin << "\t" << histo->GetBinCenter(indBin) << "\t"<< histo->GetBinWidth(indBin) <<"\t" <<histo->GetBinError(indBin)/histo->GetBinContent(indBin)*100 <<"\t"<< histo->GetBinError(indBin)/histo->GetBinContent(indBin)*histo->GetBinWidth(indBin)<<endl; //
				outputObject[commonBin].errorYSystLow = outputObject[commonBin].errorYSystLow + histo->GetBinError(indBin)/histo->GetBinContent(indBin)*histo->GetBinWidth(indBin);
// 					cout << indBin << "\t" << outputObject[commonBin].errorYSystLow<<endl; //
// 				}	else {
// 					cout << indBin << "\t" << histo->GetBinCenter(indBin) << "\t"<<(histo->GetBinContent(indBin)*histo->GetBinWidth(indBin)) << endl;
// 					outputObject[commonBin].errorYSystLow = outputObject[commonBin].errorYSystLow+(histo->GetBinContent(indBin)*histo->GetBinWidth(indBin))*histo->GetBinContent(indBin);
    
// 				}
				indBin++;
			}
			outputObject[commonBin].errorYSystLow = outputObject[commonBin].errorYSystLow/(outputObject[commonBin].errorXLow*2)*outputObject[commonBin].valueY;
			outputObject[commonBin].errorYSystHigh = outputObject[commonBin].errorYSystLow;
			outputObject[commonBin].errorYTotHigh = TMath::Sqrt(TMath::Power(outputObject[commonBin].errorYSystHigh,2) + TMath::Power(outputObject[commonBin].errorYStatHigh,2));
			outputObject[commonBin].errorYTotLow = TMath::Sqrt(TMath::Power(outputObject[commonBin].errorYSystLow,2) + TMath::Power(outputObject[commonBin].errorYStatLow,2));
			
 			cout << commonBin<< "\t"  <<outputObject[commonBin].valueX << "\t +- " << outputObject[commonBin].errorXLow << " \t " << outputObject[commonBin].valueY<< "\t+" << outputObject[commonBin].errorYSystHigh<<"\t-"<< outputObject[commonBin].errorYSystLow<< "\t+" << outputObject[commonBin].errorYTotHigh<<"\t-"<< outputObject[commonBin].errorYTotLow<<"\t syst \t"<< outputObject[commonBin].errorYSystLow/outputObject[commonBin].valueY*100 <<"%"<<endl;
		}
	} else if(ClassNameSyst.CompareTo("TGraphErrors")==0){
		cout << "\n \n TGraphErrors" << endl << endl;
		TGraphErrors *graph = (TGraphErrors*)Obj_DummySyst;
		TGraphErrors *graphCopy = (TGraphErrors*)graph->Clone("GraphCopy");
		Double_t* valueX = graphCopy->GetX();
		Double_t* valueY = graphCopy->GetY();
		Double_t* errorX = graphCopy->GetEX();
		Double_t* errorY = graphCopy->GetEY();
		for (Int_t i = 0; i < graphCopy->GetN();i++){
			if (scaleByBinCenter){
     			  	cout << i <<"\t"<<"before: " <<valueY[i] << "\t" << errorY[i] << "\t" << errorY[i]/valueY[i]*100<<"\t" << (errorX[i]+errorX[i])<< endl;
				errorY[i] = errorY[i]*(errorX[i]+errorX[i])/valueY[i];
				valueY[i] = valueY[i]*valueX[i]*(errorX[i]+errorX[i])*2*TMath::Pi();
			}	else {
				errorY[i] = errorY[i]*(errorX[i]+errorX[i])/valueY[i];
				errorY[i] = errorY[i]*(errorX[i]+errorX[i])*valueY[i] ;
			}
		}
		Int_t indBin = firstBinSyst;
		for (Int_t commonBin = 0; commonBin < nCommonBins-1; commonBin++){
			outputObject[commonBin].errorYSystLow = 0;
			for (Int_t add = 1; add < vectorRebinFactors[commonBin]+1;add++){
// 				cout << indBin << "\t" << valueY[indBin] << endl;
				outputObject[commonBin].errorYSystLow = outputObject[commonBin].errorYSystLow+errorY[indBin];
				cout << indBin << "\t" << valueY[indBin] << "\t"<< outputObject[commonBin].errorYSystLow<<endl;
				indBin++;
			}
			outputObject[commonBin].errorYSystLow = outputObject[commonBin].errorYSystLow*outputObject[commonBin].valueY/ (outputObject[commonBin].errorXLow+outputObject[commonBin].errorXHigh);
			outputObject[commonBin].errorYSystHigh = outputObject[commonBin].errorYSystLow;
			outputObject[commonBin].errorYTotHigh = TMath::Sqrt(TMath::Power(outputObject[commonBin].errorYSystHigh,2) + TMath::Power(outputObject[commonBin].errorYStatHigh,2));
			outputObject[commonBin].errorYTotLow = TMath::Sqrt(TMath::Power(outputObject[commonBin].errorYSystLow,2) + TMath::Power(outputObject[commonBin].errorYStatLow,2));
			
			cout << commonBin<< "\t"  <<outputObject[commonBin].valueX << "\t +- " << outputObject[commonBin].errorXLow << " \t " << outputObject[commonBin].valueY<< "\t+" << outputObject[commonBin].errorYSystHigh<<"\t-"<< outputObject[commonBin].errorYSystLow<< "\t+" << outputObject[commonBin].errorYTotHigh<<"\t-"<< outputObject[commonBin].errorYTotLow<<endl;
		}
	} else if(ClassNameSyst.CompareTo("TGraphAsymmErrors")==0){
		cout << "\n \n TGraphAsymmErrors" << endl << endl;
		TGraphAsymmErrors *graph = (TGraphAsymmErrors*)Obj_DummySyst;
		TGraphAsymmErrors *graphCopy = (TGraphAsymmErrors*)graph->Clone("GraphCopy");
		graphCopy->Print();
		Double_t* valueX = graphCopy->GetX();
		Double_t* valueY = graphCopy->GetY();
		Double_t* errorXlow = graphCopy->GetEXlow();
		Double_t* errorXhigh = graphCopy->GetEXhigh();
		Double_t* errorYlow = graphCopy->GetEYlow();
		Double_t* errorYhigh = graphCopy->GetEYhigh();
		for (Int_t i = 0; i < graphCopy->GetN();i++){
			if (scaleByBinCenter){
 			  	cout << i <<"\t"<<"before: " <<valueY[i] << "\t" << errorYlow[i] << "\t" << errorYlow[i]/valueY[i]*100<<"\t" << (errorXlow[i]+errorXhigh[i])<< endl;
				errorYlow[i] = errorYlow[i]*(errorXlow[i]+errorXhigh[i])/valueY[i];
				errorYhigh[i] = errorYhigh[i]*(errorXlow[i]+errorXhigh[i])/valueY[i];
				valueY[i] = valueY[i]*valueX[i]*(errorXlow[i]+errorXhigh[i])*2*TMath::Pi();
			}	else {
				errorYlow[i] = errorYlow[i]*(errorXlow[i]+errorXhigh[i])*valueY[i];
				errorYhigh[i] = errorYhigh[i]*(errorXlow[i]+errorXhigh[i])*valueY[i];
				valueY[i] = valueY[i]*(errorXlow[i]+errorXhigh[i]);
				
			}
		}
// 		graphCopy->Print();
		Int_t indBin = firstBinSyst;
		cout << firstBinSyst << endl;
		for (Int_t commonBin = 0; commonBin < nCommonBins-1; commonBin++){
			outputObject[commonBin].errorYSystLow = 0;
			outputObject[commonBin].errorYSystHigh = 0;
			for (Int_t add = 1; add < vectorRebinFactors[commonBin]+1;add++){
 				
				outputObject[commonBin].errorYSystHigh = outputObject[commonBin].errorYSystHigh + errorYhigh[indBin];
				outputObject[commonBin].errorYSystLow = outputObject[commonBin].errorYSystLow + errorYlow[indBin];
				cout << indBin << "\t" << valueY[indBin] << "\t"<< outputObject[commonBin].errorYSystHigh<<endl;
				indBin++;
			}
			outputObject[commonBin].errorYSystLow = outputObject[commonBin].errorYSystLow*outputObject[commonBin].valueY/ (outputObject[commonBin].errorXLow+outputObject[commonBin].errorXHigh);
			outputObject[commonBin].errorYSystHigh = outputObject[commonBin].errorYSystHigh*outputObject[commonBin].valueY/ (outputObject[commonBin].errorXLow+outputObject[commonBin].errorXHigh);
			outputObject[commonBin].errorYTotHigh = TMath::Sqrt(TMath::Power(outputObject[commonBin].errorYSystHigh,2) + TMath::Power(outputObject[commonBin].errorYStatHigh,2));
			outputObject[commonBin].errorYTotLow = TMath::Sqrt(TMath::Power(outputObject[commonBin].errorYSystLow,2) + TMath::Power(outputObject[commonBin].errorYStatLow,2));
			
 			cout << commonBin<< "\t"  <<outputObject[commonBin].valueX << "\t +- " << outputObject[commonBin].errorXLow << " \t " << outputObject[commonBin].valueY<< "\t+" << outputObject[commonBin].errorYSystHigh<<"\t-"<< outputObject[commonBin].errorYSystLow<< "\t+" << outputObject[commonBin].errorYTotHigh<<"\t-"<< outputObject[commonBin].errorYTotLow<<"\t syst \t"<< outputObject[commonBin].errorYSystLow/outputObject[commonBin].valueY*100 <<"%"<<endl;
		}
	} else {
		cout << " class for Syst not defined" << endl;
		return;
	}
	return;
   if (validBinsOldStat || validBinsOldSys){}
}

//*******************************************************************************************
//** Function which compares 2 Spectra with each other and returns the ratio               **
//**  statistical and systematic errors have to be handed over                             **
//**   - Graphs should be handed over as Clones of the original object, otherwise          **
//**     they will be modified                                                             **
//*******************************************************************************************
TGraphErrors* CalculateRatioBetweenSpectraWithDifferentBinning(TObject* Obj_DummyAStat, TObject* Obj_DummyASyst, TObject* Obj_DummyBStat, TObject* Obj_DummyBSyst,  Bool_t scaleByBinCenterA=kTRUE,  Bool_t scaleByBinCenterB=kTRUE, TGraphErrors** graphRebinnedAStat = NULL, TGraphErrors** graphRebinnedASyst = NULL, TGraphErrors** graphRebinnedBStat = NULL, TGraphErrors** graphRebinnedBSyst = NULL ){
	
	cout << "Reading from Object A " << endl;
	TString ClassNameA = Obj_DummyAStat->ClassName();
	Int_t nBinsXA = 0;
	if(ClassNameA.BeginsWith("TH1")){
		TH1D *histo = (TH1D*)Obj_DummyAStat;
		nBinsXA = histo->GetNbinsX()+1;
	} else if(ClassNameA.BeginsWith("TGraph")){
		TGraphErrors *graph = (TGraphErrors*)Obj_DummyAStat;
		nBinsXA = graph->GetN()+1;
	}
	Double_t binningXA[nBinsXA];
	Int_t validBinsA = GetBinning(Obj_DummyAStat, binningXA);
	
	cout << "Reading from Object B " << endl;
	TString ClassNameB = Obj_DummyBStat->ClassName();
	Int_t nBinsXB = 0 ;
	if(ClassNameB.BeginsWith("TH1")){
		TH1D *histo = (TH1D*)Obj_DummyBStat;
		nBinsXB = histo->GetNbinsX()+1;
	} else if(ClassNameB.BeginsWith("TGraph")){
		TGraphErrors *graph = (TGraphErrors*)Obj_DummyBStat;
		nBinsXB = graph->GetN()+1;
	}
	Double_t binningXB[nBinsXB];
	Int_t validBinsB= GetBinning(Obj_DummyBStat, binningXB);
	
	Int_t nBinsComb;
	if (nBinsXB < nBinsXA){
		nBinsComb= nBinsXA;
	} else {
		nBinsComb= nBinsXB;
	}
// 	for (Int_t i = 0; i < nBinsXB; i++){
// 		cout << binningXB[i] << "\t," ;
// 	}
	Double_t binningCombined[nBinsComb];
	Int_t binningCombinedBinsToBeMergedA[nBinsComb];
	Int_t binningCombinedBinsToBeMergedB[nBinsComb];
	Int_t nBinsNew = CompareBinning( validBinsA, binningXA,validBinsB, binningXB, binningCombined, binningCombinedBinsToBeMergedA,binningCombinedBinsToBeMergedB, "A");
	
	cout << "Object A"  << endl;
	AliConvDataObject rebinnedA[nBinsComb];
	RebinObjects(Obj_DummyAStat,Obj_DummyASyst, binningCombined, binningCombinedBinsToBeMergedA, nBinsNew,rebinnedA, Obj_DummyAStat->ClassName(), Obj_DummyASyst->ClassName(),scaleByBinCenterA);

	cout << "Object B"  << endl;	
	AliConvDataObject rebinnedB[nBinsComb];
	RebinObjects(Obj_DummyBStat,Obj_DummyBSyst, binningCombined, binningCombinedBinsToBeMergedB, nBinsNew,rebinnedB, Obj_DummyBStat->ClassName(), Obj_DummyBSyst->ClassName(), scaleByBinCenterB);
	
	Double_t ratioX[nBinsComb];
	Double_t errorX[nBinsComb];
	Double_t ratioY[nBinsComb];
	Double_t errorY[nBinsComb];
	cout << nBinsNew-1 << "\t" <<nBinsComb << endl;
	for (Int_t commonBin = 0; commonBin < nBinsNew-1; commonBin++){
		cout << "A "<<commonBin<< "\t"  <<rebinnedA[commonBin].valueX << "\t +- " << rebinnedA[commonBin].errorXLow << " \t " << rebinnedA[commonBin].valueY<< "\t+" << rebinnedA[commonBin].errorYStatHigh<<"\t-"<< rebinnedA[commonBin].errorYStatLow<< "\t+" << rebinnedA[commonBin].errorYSystHigh<<"\t-"<< rebinnedA[commonBin].errorYSystLow<< "\t+" << rebinnedA[commonBin].errorYTotHigh<<"\t-"<< rebinnedA[commonBin].errorYTotLow<< "\t" << rebinnedA[commonBin].errorYTotHigh/rebinnedA[commonBin].valueY*100 << "%"<<endl;
		cout << "B " <<commonBin<< "\t"  <<rebinnedB[commonBin].valueX << "\t +- " << rebinnedB[commonBin].errorXLow << " \t " << rebinnedB[commonBin].valueY<< "\t+" << rebinnedB[commonBin].errorYStatHigh<<"\t-"<< rebinnedB[commonBin].errorYStatLow<< "\t+" << rebinnedB[commonBin].errorYSystHigh<<"\t-"<< rebinnedB[commonBin].errorYSystLow<< "\t+" << rebinnedB[commonBin].errorYTotHigh<<"\t-"<< rebinnedB[commonBin].errorYTotLow<< "\t" << rebinnedB[commonBin].errorYTotHigh/rebinnedB[commonBin].valueY*100 << "%"<<endl;
		ratioX[commonBin] = rebinnedA[commonBin].valueX;
		errorX[commonBin] = rebinnedA[commonBin].errorXHigh;
		ratioY[commonBin] = rebinnedA[commonBin].valueY/rebinnedB[commonBin].valueY;
		errorY[commonBin] = TMath::Sqrt(TMath::Power(rebinnedA[commonBin].errorYTotHigh/rebinnedA[commonBin].valueY,2) +TMath::Power(rebinnedB[commonBin].errorYTotHigh/rebinnedB[commonBin].valueY,2))*ratioY[commonBin];
		cout << "Ratio: " << ratioX[commonBin] << "\t" <<  errorX[commonBin] << "\t" << ratioY[commonBin] << "\t" << errorY[commonBin] << "\t" << errorY[commonBin]/ratioY[commonBin]*100 << "%"<<endl;
	}
	
    
	Double_t rebinnedSpectrumAY[nBinsNew-1];
	Double_t rebinnedSpectrumAYStatErr[nBinsNew-1];
	Double_t rebinnedSpectrumAYSysErr[nBinsNew-1];
	Double_t rebinnedSpectrumBY[nBinsNew-1];
	Double_t rebinnedSpectrumBYStatErr[nBinsNew-1];
	Double_t rebinnedSpectrumBYSysErr[nBinsNew-1];
	for (Int_t commonBin = 0; commonBin < nBinsNew-1; commonBin++){
		rebinnedSpectrumAY[commonBin] = rebinnedA[commonBin].valueY;
		rebinnedSpectrumBY[commonBin] = rebinnedB[commonBin].valueY;
		rebinnedSpectrumAYStatErr[commonBin] = rebinnedA[commonBin].errorYStatHigh;
		rebinnedSpectrumBYStatErr[commonBin] = rebinnedB[commonBin].errorYStatHigh;
		rebinnedSpectrumAYSysErr[commonBin] = TMath::Sqrt(rebinnedA[commonBin].errorYSystHigh*rebinnedA[commonBin].errorYSystHigh - 0.045*rebinnedA[commonBin].valueY*0.045*rebinnedA[commonBin].valueY) ;
		rebinnedSpectrumBYSysErr[commonBin] = rebinnedB[commonBin].errorYSystHigh;
	}   
		
	(*graphRebinnedAStat) =  new TGraphErrors(nBinsNew-1,ratioX,rebinnedSpectrumAY,errorX,rebinnedSpectrumAYStatErr); 
	(*graphRebinnedASyst) =  new TGraphErrors(nBinsNew-1,ratioX,rebinnedSpectrumAY,errorX,rebinnedSpectrumAYSysErr); 
	(*graphRebinnedBStat) =  new TGraphErrors(nBinsNew-1,ratioX,rebinnedSpectrumBY,errorX,rebinnedSpectrumBYStatErr); 
	(*graphRebinnedBSyst) =  new TGraphErrors(nBinsNew-1,ratioX,rebinnedSpectrumBY,errorX,rebinnedSpectrumBYSysErr); 
		
	TGraphErrors* returnGraph =  new TGraphErrors(nBinsNew-1,ratioX,ratioY,errorX,errorY); 
	return returnGraph;
	
}

Int_t GetNEvents (TH1D* histo){
	Int_t nEvents = histo->GetBinContent(1)+(histo->GetBinContent(1)/(histo->GetBinContent(1)+histo->GetBinContent(5)))*histo->GetBinContent(6);
	// BinContent 1 - good events
	// BinContent 2 - centrality not selected
	// BinContent 3 - MC corrupt
	// BinContent 4 - no Trigger Bit
	// BinContent 5 - Zvertex-position, 
	// BinContent 6 - no Contributors to vtx
	// BinContent 7 - PileUp
	cout <<"nEvents new: "<< nEvents << "\t nEvents old: "<< histo->GetEntries()-histo->GetBinContent(5)-histo->GetBinContent(7)-histo->GetBinContent(4)<< endl; 
	return nEvents;
}

Int_t GetNEvents (TH1F* histo){
	Int_t nEvents = histo->GetBinContent(1)+(histo->GetBinContent(1)/(histo->GetBinContent(1)+histo->GetBinContent(5)))*histo->GetBinContent(6);
	// BinContent 1 - good events
	// BinContent 2 - centrality not selected
	// BinContent 3 - MC corrupt
	// BinContent 4 - no Trigger Bit
	// BinContent 5 - Zvertex-position, 
	// BinContent 6 - no Contributors to vtx
	// BinContent 7 - PileUp
	cout <<"nEvents new: "<< nEvents << "\t nEvents old: "<< histo->GetEntries()-histo->GetBinContent(5)-histo->GetBinContent(7)-histo->GetBinContent(4)<< endl; 
	return nEvents;
}

TString AnalyseTPCClusterCut(Int_t clsTPCCut){   
	switch(clsTPCCut){
	case 0: // 0
		return "min. TPC cluster: 0";
	case 1:  // 60
		return "min. TPC cluster: 60";
	case 2:  // 80
		return "min. TPC cluster: 80";
	case 3:  // 100
		return "min. TPC cluster: 100";     
	case 4:  // 95% of findable clusters
		return "TPC cluster/uncorr findable Clusters: 0.95";     
	case 5:  // 0% of findable clusters
		return "TPC cluster/corr findable Clusters: 0.";     
	case 6:  // 80% of findable clusters
		return "TPC cluster/corr findable Clusters: 0.8";     
	case 7:  // 0% of findable clusters
		return "TPC cluster/uncorr findable Clusters: 0.35";          
	case 8:
		return "TPC cluster/corr findable Clusters: 0.35";          
	case 9:
		return "TPC cluster/corr findable Clusters: 0.6";          
	default:
		return "no cluster cut defined";
	}  
}

TString AnalyseTPCdEdxCutElectronLine(Int_t ededxSigmaCut){
	switch(ededxSigmaCut){
	case 0: // -10,10
		return "-10 < #sigma_{e} < 10";
	case 1: // -5,5
		return "-5 < #sigma_{e} < 5";
	case 2: // -3,5
		return "-3 < #sigma_{e} < 5";
	case 3: // -4,5
		return "-4 < #sigma_{e} < 5";
	case 4: // -6,7
		return "-6 < #sigma_{e} < 7";
	case 5: // -4,4
		return "-4 < #sigma_{e} < 4";
	case 6: // -2.5,4
		return "-2.5 < #sigma_{e} < 4";
	case 7: // -2,3.5
		return "-2 < #sigma_{e} < 3.5";
	default:
		return "no dEdx cut defined";
	}
	return kTRUE;
}

///________________________________________________________________________
TString AnalyseAlphaMesonCut(Int_t alphaMesonCut){ 
	switch(alphaMesonCut){
	case 0:	// 0- 0.7
		return "|#alpha_{meson}| < 0.7";
	case 1:	// 0-0.5
		return "|#alpha_{meson}| < 0.5";
	case 2:	// 0.5-1
		return "0.5 < |#alpha_{meson}| < 1";
	case 3:	// 0.0-1
		return "|#alpha_{meson}| < 1";
	case 4:	// 0-0.65
		return "|#alpha_{meson}| < 0.65";
	case 5:	// 0-0.75
		return "|#alpha_{meson}| < 0.75";
	case 6:	// 0-0.8
		return "|#alpha_{meson}| < 0.8";
	case 7:	// 0.0-0.85
		return "|#alpha_{meson}| < 0.85";
	case 8:	// 0.0-0.6
		return "|#alpha_{meson}| < 0.6";
	case 9:	// 0.0-0.3
		return "|#alpha_{meson}| < 0.3";
	default:
		return "no alpha cut defined";
	}
}

TString AnalyseMCSmearingCut(Int_t mcSmearingCut){ 
	switch(mcSmearingCut){
	case 0:
		return "no additional smearing" ;
	case 1:
		return "#sqrt{0.011^{2} + 0.007^{2} #times P_{#gamma}^{2}} #times Rand.Gaus(0,1)" ;
	case 2:
		return "#sqrt{0.022^{2} + 0.007^{2} #times P_{#gamma}^{2}} #times Rand.Gaus(0,1)" ;
	case 3:
		return "#sqrt{0.044^{2} + 0.007^{2} #times P_{#gamma}^{2}} #times Rand.Gaus(0,1)" ;
	case 7:
		return "#sqrt{0.011^{2} + 0.014^{2} #times P_{#gamma}^{2}} #times Rand.Gaus(0,1)" ;
	case 8:
		return "#sqrt{0.011^{2} + 0.0035^{2} #times P_{#gamma}^{2}} #times Rand.Gaus(0,1)" ;
	case 9:
		return "#sqrt{0.011^{2} + 0.028^{2} #times P_{#gamma}^{2}} #times Rand.Gaus(0,1)" ;
	default:
		return "smearing cut not defined";
	}
}

TString AnalyseTPCdEdxCutPionLine(TString sPionCut){
	cout << sPionCut << endl;
	TString sPidedxSigmaCut = sPionCut(0,1);
	Int_t pidedxSigmaCut = sPidedxSigmaCut.Atoi();
	cout << "pidedxSigmaCut: " << pidedxSigmaCut << endl;
	Double_t fPIDnSigmaAbovePionLine = 0;
	Double_t fPIDnSigmaAbovePionLineHighPt=0;
	switch(pidedxSigmaCut){
	case 0:  // -10
		fPIDnSigmaAbovePionLine=-10;
		fPIDnSigmaAbovePionLineHighPt=-10;
		break;
	case 1:   // 0
		fPIDnSigmaAbovePionLine=0;
		fPIDnSigmaAbovePionLineHighPt=-10;
		break;
	case 2:  // 1
		fPIDnSigmaAbovePionLine=1;
		fPIDnSigmaAbovePionLineHighPt=-10;
		break;
	case 3:  // 1
		fPIDnSigmaAbovePionLine=2.5;
		fPIDnSigmaAbovePionLineHighPt=-10;
		break;
	case 4:  // 1
		fPIDnSigmaAbovePionLine=0.5;
		fPIDnSigmaAbovePionLineHighPt=-10;
		break;
	case 5:  // 1
		fPIDnSigmaAbovePionLine=2.;
		fPIDnSigmaAbovePionLineHighPt=-10;
		break;
	case 6:  // 1
		fPIDnSigmaAbovePionLine=2.;
		fPIDnSigmaAbovePionLineHighPt=0.5;
		break;
	case 7:  // 1
		fPIDnSigmaAbovePionLine=3.5;
		fPIDnSigmaAbovePionLineHighPt=-10;
		break;
	case 8:  // 1
		fPIDnSigmaAbovePionLine=2.;
		fPIDnSigmaAbovePionLineHighPt=1.;
		break;
	case 9:
		fPIDnSigmaAbovePionLine=3.0; // We need a bit less tight cut on dE/dx
		fPIDnSigmaAbovePionLineHighPt=-10;
		break;
	default:
		cout << "pion line cut unknown" << endl;
	}

	TString sPiMomdedxSigmaCut = sPionCut(1,1);
	Int_t piMomdedxSigmaCut = sPiMomdedxSigmaCut.Atoi();
	Double_t fPIDMinPnSigmaAbovePionLine = 0;
	cout << "piMomdedxSigmaCut: " << piMomdedxSigmaCut << endl;
	switch(piMomdedxSigmaCut){
	case 0:  // 0.5 GeV
		fPIDMinPnSigmaAbovePionLine=0.5;
		break;
	case 1:  // 1. GeV
		fPIDMinPnSigmaAbovePionLine=1.;
		break;
	case 2:  // 1.5 GeV
		fPIDMinPnSigmaAbovePionLine=1.5;
		break;
	case 3:  // 20.0 GeV
		fPIDMinPnSigmaAbovePionLine=20.;
		break;
	case 4:  // 50.0 GeV
		fPIDMinPnSigmaAbovePionLine=50.;
		break;
	case 5:  // 0.3 GeV
		fPIDMinPnSigmaAbovePionLine=0.3;
		break;
	case 6:  // 0.25 GeV
		fPIDMinPnSigmaAbovePionLine=0.25;
		break;
	case 7:  // 0.4 GeV
		fPIDMinPnSigmaAbovePionLine=0.4;
		break;
	case 8:  // 0.2 GeV
		fPIDMinPnSigmaAbovePionLine=0.2;
		break;
	default:
		cout << "pion line minimum pt cut unknown" << endl;
		
	}

	TString sPiMaxMomdedxSigmaCut = sPionCut(2,1);
	Int_t piMaxMomdedxSigmaCut = sPiMaxMomdedxSigmaCut.Atoi();
	Double_t fPIDMaxPnSigmaAbovePionLine = 0;
	cout << "piMaxMomdedxSigmaCut: " << piMaxMomdedxSigmaCut << endl;
	switch(piMaxMomdedxSigmaCut){
	case 0:  // 100. GeV
		fPIDMaxPnSigmaAbovePionLine=100.;
		break;
	case 1:  // 5. GeV
		fPIDMaxPnSigmaAbovePionLine=5.;
		break;
	case 2:  // 4. GeV
		fPIDMaxPnSigmaAbovePionLine=4.;
		break;
	case 3:  // 3.5 GeV
		fPIDMaxPnSigmaAbovePionLine=3.5;
		break;
	case 4:  // 3. GeV
		fPIDMaxPnSigmaAbovePionLine=3.;
		break;
	case 5:  // 7. GeV
		fPIDMaxPnSigmaAbovePionLine=7.;
		break;
	default:
		cout << "pion line minimum pt cut unknown" << endl;
	}
	return Form("rejected #pi for #sigma_{#pi} < %.2f (%.2f GeV/c < p_{T}_{#pi} < %.2f GeV/c), #sigma_{#pi} < %.2f (p_{T}_{#pi} > %.2f GeV/c)",fPIDnSigmaAbovePionLine, fPIDMinPnSigmaAbovePionLine,fPIDMaxPnSigmaAbovePionLine,fPIDnSigmaAbovePionLineHighPt,fPIDMaxPnSigmaAbovePionLine);
}

TString AnalyseChi2GammaCut(Int_t chi2GammaCut, Int_t psiPairCut){   // Set Cut

	TString psiPairCutString = "";
	Bool_t k2DPsiPairChi2 = kFALSE;
	switch(psiPairCut) {
		case 0:
			psiPairCutString = "|#Psi_{Pair}| < 10000";
			break;
		case 1:
			psiPairCutString = "|#Psi_{Pair}| < 0.1";
			break;
		case 2:
			psiPairCutString = "|#Psi_{Pair}| < 0.05";
			break;
		case 3:
			psiPairCutString = "|#Psi_{Pair}| < 0.035";
			break;
		case 4:
			psiPairCutString = "|#Psi_{Pair}| < 0.2";
			break;
		case 5:
			psiPairCutString = "|#Psi_{Pair}| < 0.1";
			k2DPsiPairChi2= kTRUE;
			break;
		case 6:
			psiPairCutString = "|#Psi_{Pair}| < 0.05";
			k2DPsiPairChi2= kTRUE;
			break;
		case 7:
			psiPairCutString = "|#Psi_{Pair}| < 0.035";
			k2DPsiPairChi2= kTRUE;
			break;
		case 8:
			psiPairCutString =  "|#Psi_{Pair}| < 0.2";
			k2DPsiPairChi2= kTRUE;
			break;
		case 9:
			psiPairCutString =  "|#Psi_{Pair}| < 0.5";
			break;
		default:
			psiPairCutString =  "#Psi_{Pair} cut not defined";
			break;
	}
	
	TString chi2CutString = "";
	switch(chi2GammaCut){
	case 0: // 100
		chi2CutString = "#chi_{#gamma}^{2} < 100";
		break;
	case 1:  // 50
		chi2CutString = "#chi_{#gamma}^{2} < 50";
		break;
	case 2:  // 30
		chi2CutString = "#chi_{#gamma}^{2} < 30";
		break;
	case 3:
		chi2CutString = "#chi_{#gamma}^{2} < 200";
		break;
	case 4:
		chi2CutString = "#chi_{#gamma}^{2} < 500";
		break;
	case 5:
		chi2CutString = "#chi_{#gamma}^{2} < 1000000";
		break;
	case 6:
		chi2CutString = "#chi_{#gamma}^{2} < 5";
		break;
	case 7:
		chi2CutString = "#chi_{#gamma}^{2} < 10";
		break;
	case 8:
		chi2CutString = "#chi_{#gamma}^{2} < 20";
		break;
	case 9:
		chi2CutString = "#chi_{#gamma}^{2} < 15";
		break;
	default:
		chi2CutString = "#chi_{#gamma}^{2} cut unknown";
		break;
	}
	if (k2DPsiPairChi2){
		return Form("2D cut: %s, %s", chi2CutString.Data(), psiPairCutString.Data());
	} else {
		return Form("1D cut: %s", chi2CutString.Data());
	}
	return "";
}

TString AnalyseQtMaxCut(Int_t QtMaxCut){   
		
	switch(QtMaxCut){
	case 0: //
		return "no q_{T}_{#gamma} cut applied";
	case 1:
		return "q_{T}_{#gamma} < 0.1 GeV/c";
	case 2:
		return "q_{T}_{#gamma} < 0.07 GeV/c";
	case 3:
		return "q_{T}_{#gamma} < 0.05 GeV/c";
	case 4:
		return "q_{T}_{#gamma} < 0.03 GeV/c";
	case 5:
		return "q_{T}_{#gamma} < 0.02 GeV/c";
	case 6:
		return "2D ellipse q_{T}_{#gamma} < 0.02 GeV/c, #alpha < 0.95";
	case 7:
		return "q_{T}_{#gamma} < 0.15 GeV/c";  
	case 8:
		return "2D ellipse q_{T}_{#gamma} < 0.05 GeV/c, #alpha < 0.95";
	case 9:
		return "2D ellipse q_{T}_{#gamma} < 0.03 GeV/c, #alpha < 0.95";
	default:
		return "no q_{T}_{#gamma} cut defined";
	}
}

TString AnalyseSinglePtCut(Int_t singlePtCut){
	switch(singlePtCut){
	case 0: // 0.050 GeV
		return "p_{T}_{e^{#pm}} > 0.050 GeV/c";
	case 1:  // 0.100 GeV
		return "p_{T}_{e^{#pm}} > 0.100 GeV/c";
	case 2:  // 0.150 GeV
		return "p_{T}_{e^{#pm}} > 0.150 GeV/c";
	case 3:  // 0.200 GeV
		return "p_{T}_{e^{#pm}} > 0.200 GeV/c";
	case 4:  // 0.075 GeV
		return "p_{T}_{e^{#pm}} > 0.075 GeV/c";
	case 5:  // 0.125 GeV
		return "p_{T}_{e^{#pm}} > 0.125 GeV/c";
	case 6:  // 0.04 GeV
		return "p_{T}_{e^{#pm}} > 0.040 GeV/c";
	case 7:  // 0.0 GeV
		return "p_{T}_{e^{#pm}} > 0.0 GeV/c";
	default:
		return "p_{T}_{e^{#pm}} cut not defined";
	}
}

TString AnalyseDCAZPhotonCut(Int_t dcaZPhoton){
	switch(dcaZPhoton){
		case 0:  //
			return "|dca_{Z}| < 1000 cm"; 
		case 1:  //
			return "|dca_{Z}| < 10 cm"; 
		case 2:  //
			return "|dca_{Z}| < 5 cm"; 
		case 3:  //
			return "|dca_{Z}| < 4 cm"; 
		case 4:  //
			return "|dca_{Z}| < 3 cm"; 
		case 5:  //
			return "|dca_{Z}| < 2.5 cm"; 
		case 6:  //
			return "|dca_{Z}| < 2 cm"; 
		case 7:  //
			return "|dca_{Z}| < 1.5 cm"; 
		case 8:  //
			return "|dca_{Z}| < 1 cm"; 
		case 9:  //
			return "|dca_{Z}| < 0.5 cm"; 
		default:
			return "|dca_{Z}| cut not defined";
	}
}

Double_t AnalyseDCAZPhotonCutValue(Int_t dcaZPhoton){
	switch(dcaZPhoton){   
		case 0:  //
			return 1000.; 
		case 1:  //
			return 10.; 
		case 2:  //
			return 5.; 
		case 3:  //
			return 4.; 
		case 4:  //
			return 3.; 
		case 5:  //
			return 2.5; 
		case 6:  //
			return 2.; 
		case 7:  //
			return 1.5; 
		case 8:  //
			return 1.; 
		case 9:  //
			return 0.5; 
		default:
			return 1000;
	}
}

TString AnalyseCosPointCut(Int_t cosPoint){
	switch(cosPoint){
		case 0:  //
			return "cos(#Theta_{point}) > -1"; 
		case 1:  //
			return "cos(#Theta_{point}) > 0"; 
		case 2:  //
			return "cos(#Theta_{point}) > 0.5"; 
		case 3:  //
			return "cos(#Theta_{point}) > 0.75"; 
		case 4:  //
			return "cos(#Theta_{point}) > 0.85"; 
		case 5:  //
			return "cos(#Theta_{point}) > 0.88"; 
		case 6:  //
			return "cos(#Theta_{point}) > 0.9"; 
		case 7:  //
			return "cos(#Theta_{point}) > 0.95"; 
		default:
			return "cos(#Theta_{point}) cut not defined";
	}
}


///________________________________________________________________________
TString AnalyseBackgroundScheme(TString sBackgroundScheme){
	TString sBackgroundSchemeB = sBackgroundScheme(0,1);
	Int_t BackgroundScheme = sBackgroundSchemeB.Atoi();
	TString bGScheme = "";
	cout << "BackgroundScheme: " << BackgroundScheme << endl;

	
	switch(BackgroundScheme){
	case 0: //Rotation
		bGScheme="Rotation";
		break;
	case 1: // mixed event with V0 multiplicity
		bGScheme="mixed event, V0 mult";
		break;
	case 2: // mixed event with track multiplicity
		bGScheme="mixed event, track mult";
		break;
	case 3: //Rotation
		bGScheme="Rotation, with prob.";
		break;
	case 4: //No BG calculation
		bGScheme="no BG calculated";
		break;
	case 5: //Rotation
		bGScheme="Rotation, new BG handler";
		break;
	case 6: // mixed event with V0 multiplicity
		bGScheme="mixed event, new BG handler, V0 mult";
		break;
	case 7: // mixed event with track multiplicity
		bGScheme="mixed event, new BG handler, track mult";
		break;
	case 8: //Rotation
		bGScheme="Rotation, new BG handler, with prob. ";
		break;
	default:
		bGScheme="no BG method selected";
		
	}
	
	TString sNumberOfBGEvents = sBackgroundScheme(1,1);
	Int_t NumberOfBGEvents = sNumberOfBGEvents.Atoi();
	Int_t fNumberOfBGEvents = 0;
	cout << "NumberOfBGEvents: " << NumberOfBGEvents << endl;

	switch(NumberOfBGEvents){
	case 0:
		fNumberOfBGEvents = 5;
		break;
	case 1:
		fNumberOfBGEvents = 10;
		break;
	case 2:
		fNumberOfBGEvents = 15;
		break;
	case 3:
		fNumberOfBGEvents = 20;
		break;
	case 4:
		fNumberOfBGEvents = 2;
		break;
	case 5:
		fNumberOfBGEvents = 50;
		break;
	case 6:
		fNumberOfBGEvents = 80;
		break;
	case 7:
		fNumberOfBGEvents = 100;
		break;
	default:
		cout<<"Warning: NumberOfBGEvents not defined "<<NumberOfBGEvents<<endl;
	}

	TString sDegreesForRotationMethod = sBackgroundScheme(1,2);
	Int_t DegreesForRotationMethod = sDegreesForRotationMethod.Atoi();
	Int_t fnDegreeRotationPMForBG = 0;
	cout << "DegreesForRotationMethod: " << DegreesForRotationMethod << endl;
	
	switch(DegreesForRotationMethod){
	case 0:
		fnDegreeRotationPMForBG = 5;
		break;
	case 1:
		fnDegreeRotationPMForBG = 10;
		break;
	case 2:
		fnDegreeRotationPMForBG = 15;
		break;
	case 3:
		fnDegreeRotationPMForBG = 20;
		break;
	default:
		cout<<"Warning: DegreesForRotationMethod not defined "<<DegreesForRotationMethod<<endl;
	
	}

	if (BackgroundScheme == 0 || BackgroundScheme == 3 || BackgroundScheme == 5 || BackgroundScheme == 8){
		return Form("%s within #pm %i°, %i photons per pool",bGScheme.Data(), fnDegreeRotationPMForBG, fNumberOfBGEvents)   ;
	} else {
		return Form("%s with %i photons per pool",bGScheme.Data(), fNumberOfBGEvents)   ;
	}
	return "";
}

TString AnalyseEtaCut(Int_t etaCut){ 

	switch(etaCut){
	case 0: // 0.9
		return "#eta_{#gamma,e^{#pm}} < 0.9";
	case 1:  // 0.6
		return "#eta_{#gamma,e^{#pm}} < 0.6";
	case 2:  // 1.4
		return "#eta_{#gamma,e^{#pm}} < 1.4";
	case 3: // 0.65
		return "#eta_{#gamma,e^{#pm}} < 0.65";
	case 4: // 0.75
		return "#eta_{#gamma,e^{#pm}} < 0.75";
	case 5: // 0.5
		return "#eta_{#gamma,e^{#pm}} < 0.5";
	case 6: // 5.
		return "#eta_{#gamma,e^{#pm}} < 5";
	case 7: // 0.1 - 0.8
		return "#eta_{#gamma,e^{#pm}} < 0.7";
	case 8: // 0.1 - 0.8
		return "#eta_{#gamma,e^{#pm}} < 0.4";
	case 9: // 10
		return "#eta_{#gamma,e^{#pm}} < 10";
	default:
		return "no #eta_{#gamma,e^{#pm}} cut defined";
	}
}

TString AnalyseEtaCutpPb(Int_t etaCut){ 

	switch(etaCut){
	case 0: // 0.9
		return "#eta_{#gamma,e^{#pm}} < 0.9";
	case 1:  // 1.2
		return "#eta_{#gamma,e^{#pm}} < 0.6";
	case 2:  // 1.4
		return "#eta_{#gamma,e^{#pm}} < 1.4";
	case 3: // 0.8
		return "#eta_{#gamma,e^{#pm}} < 0.8";
	case 4: // 0.75
		return "#eta_{#gamma,e^{#pm}} < 0.75";
	case 5: // 0.9 - 1.4
		return "#eta_{#gamma,e^{#pm}} < 0.5";
	case 6: // 5.
		return "#eta_{#gamma,e^{#pm}} < 5";
	case 7: // 0.1 - 0.8
		return "#eta_{#gamma,e^{#pm}} < 0.3";
	case 8: // 0.1 - 0.8
		return "#eta_{#gamma,e^{#pm}} < 0.4";
	case 9: // 10
		return "#eta_{#gamma,e^{#pm}} < 10";
	default:
		return "no #eta_{#gamma,e^{#pm}} cut defined";
	}
}


TString AnalyseRCut(Int_t RCut){
	// Set Cut
	switch(RCut){
	case 0:
		return "0 cm < R_{conv, #gamma} < 180 cm";
	case 1:
		return "2.8 cm < R_{conv, #gamma} < 180 cm";
	case 2:
		return "5 cm < R_{conv, #gamma} < 180 cm";
	case 3:
		return "10 cm < R_{conv, #gamma} < 70 cm";
	case 4:
		return "5 cm < R_{conv, #gamma} < 70 cm";
	case 5:
		return "10 cm < R_{conv, #gamma} < 180 cm";
	case 6:
		return "20 cm < R_{conv, #gamma} < 180 cm";
	case 7:
		return "35 cm < R_{conv, #gamma} < 180 cm"; // before 26 (19.4.2014)qdel 
	case 8:
		return "12.5 cm < R_{conv, #gamma} < 180 cm";
	case 9:
		return "7.5 cm < R_{conv, #gamma} < 180 cm";
	default:
		return "R cut not defined";
	}
}

TString AnalyseRCutAndQuality(Int_t RCut, Int_t photonQualitCut){
	// Set Cut
	TString stringRCut = "";
	switch(RCut){
	case 0:
		stringRCut= "0 cm < R_{conv, #gamma} < 180 cm";
		break;
	case 1:
		stringRCut= "2.8 cm < R_{conv, #gamma} < 180 cm";
		break;
	case 2:
		stringRCut=  "5 cm < R_{conv, #gamma} < 180 cm";
		break;
	case 3:
		stringRCut= "10 cm < R_{conv, #gamma} < 70 cm";
		break;
	case 4:
		stringRCut= "5 cm < R_{conv, #gamma} < 70 cm";
		break;
	case 5:
		stringRCut= "10 cm < R_{conv, #gamma} < 180 cm";
		break;
	case 6:
		stringRCut= "20 cm < R_{conv, #gamma} < 180 cm";
		break;
	case 7:
		stringRCut= "35 cm < R_{conv, #gamma} < 180 cm"; // before 26 (19.4.2014)qdel 
		break;
	case 8:
		stringRCut= "12.5 cm < R_{conv, #gamma} < 180 cm";
		break;
	case 9:
		stringRCut= "7.5 cm < R_{conv, #gamma} < 180 cm";
		break;
	default:
		stringRCut= "R cut not defined";
		break;
	}
	
	TString stringPhotonQuality = "";
	switch(photonQualitCut) {
		case 0:
			stringPhotonQuality =  "photon Quality: 1,2,3";
			break;
		case 2:
			stringPhotonQuality =  "photon Quality: 1";
			break;
		case 3:
			stringPhotonQuality =  "photon Quality: 2";
			break;
		case 4:
			stringPhotonQuality =  "photon Quality: 3";
			break;
		default:
			stringPhotonQuality =  "photon Quality cut not defined";
			break;
	}
	return Form("%s, %s", stringRCut.Data(), stringPhotonQuality.Data());
	
}


TString AnalysePsiPair(Int_t PsiPairCut, Int_t chi2GammaCut){
	
		
	TString psiPairCutString = "";
	Bool_t k2DPsiPairChi2 = kFALSE;
	switch(PsiPairCut) {
		case 0:
			psiPairCutString = "|#Psi_{Pair}| < 10000";
			break;
		case 1:
			psiPairCutString = "|#Psi_{Pair}| < 0.1";
			break;
		case 2:
			psiPairCutString = "|#Psi_{Pair}| < 0.05";
			break;
		case 3:
			psiPairCutString = "|#Psi_{Pair}| < 0.035";
			break;
		case 4:
			psiPairCutString = "|#Psi_{Pair}| < 0.2";
			break;
		case 5:
			psiPairCutString = "|#Psi_{Pair}| < 0.1";
			k2DPsiPairChi2= kTRUE;
			break;
		case 6:
			psiPairCutString = "|#Psi_{Pair}| < 0.05";
			k2DPsiPairChi2= kTRUE;
			break;
		case 7:
			psiPairCutString = "|#Psi_{Pair}| < 0.035";
			k2DPsiPairChi2= kTRUE;
			break;
		case 8:
			psiPairCutString =  "|#Psi_{Pair}| < 0.2";
			k2DPsiPairChi2= kTRUE;
			break;
		case 9:
			psiPairCutString =  "|#Psi_{Pair}| < 0.5";
			break;
		default:
			psiPairCutString =  "#Psi_{Pair} cut not defined";
			break;
	}
	
	TString chi2CutString = "";
	switch(chi2GammaCut){
	case 0: // 100
		chi2CutString = "#chi_{#gamma}^{2} < 100";
		break;
	case 1:  // 50
		chi2CutString = "#chi_{#gamma}^{2} < 50";
		break;
	case 2:  // 30
		chi2CutString = "#chi_{#gamma}^{2} < 30";
		break;
	case 3:
		chi2CutString = "#chi_{#gamma}^{2} < 200";
		break;
	case 4:
		chi2CutString = "#chi_{#gamma}^{2} < 500";
		break;
	case 5:
		chi2CutString = "#chi_{#gamma}^{2} < 1000000";
		break;
	case 6:
		chi2CutString = "#chi_{#gamma}^{2} < 5";
		break;
	case 7:
		chi2CutString = "#chi_{#gamma}^{2} < 10";
		break;
	case 8:
		chi2CutString = "#chi_{#gamma}^{2} < 20";
		break;
	case 9:
		chi2CutString = "#chi_{#gamma}^{2} < 15";
		break;
	default:
		chi2CutString = "#chi_{#gamma}^{2} cut unknown";
		break;
	}
	if (k2DPsiPairChi2){
		return Form("2D cut: %s, %s", psiPairCutString.Data(), chi2CutString.Data());
	} else {
		return Form("1D cut: %s", psiPairCutString.Data());
	}
	return "";
}   

TString AnalysePsiPairAndR(Int_t PsiPairCut, Int_t RCut ){
	TString psiPairCut = "";
	switch(PsiPairCut) {
	case 0:
		psiPairCut = "|#Psi_{Pair}| < 10000";
		break;
	case 1:
		psiPairCut = "|#Psi_{Pair}| < 0.1";
		break;
	case 2:
		psiPairCut = "|#Psi_{Pair}| < 0.05";
		break;
	case 3:
		psiPairCut = "|#Psi_{Pair}| < 0.035";
		break;
	case 4:
		psiPairCut = "|#Psi_{Pair}| < 0.15";
		break;
	case 5:
		psiPairCut = "#Psi_{Pair} < (0.1 - 0.1/1 * #Delta #Phi )";
		break;
	case 6:
		psiPairCut = "#Psi_{Pair} < (0.05 - 0.05/1 * #Delta #Phi )";
		break;
	case 7:
		psiPairCut = "#Psi_{Pair} < (0.035 - 0.035/1 * #Delta #Phi )";
		break;
	case 8:
		psiPairCut = "#Psi_{Pair} < (0.2 - 0.2/1 * #Delta #Phi )";
		break;
	case 9:
		psiPairCut = "#Psi_{Pair} < 0.5";
		break;
	default:
		psiPairCut = "#Psi_{Pair} cut not defined";
		break;
	}
	TString RCutString = "";
	switch(RCut){
	case 0:
		RCutString =  "0 cm < R_{conv, #gamma} < 180 cm";
		break;
	case 1:
		RCutString =  "2.8 cm < R_{conv, #gamma} < 180 cm";
		break;
	case 2:
		RCutString =  "5 cm < R_{conv, #gamma} < 180 cm";
		break;
	case 3:
		RCutString =  "10 cm < R_{conv, #gamma} < 70 cm";
		break;
	case 4:
		RCutString =  "5 cm < R_{conv, #gamma} < 70 cm";
		break;
	case 5:
		RCutString =  "10 cm < R_{conv, #gamma} < 180 cm";
		break;
	case 6:
		RCutString =  "20 cm < R_{conv, #gamma} < 180 cm";
		break;
	case 7:
		RCutString =  "26 cm < R_{conv, #gamma} < 180 cm";
		break;
	case 8:
		RCutString =  "35 cm < R_{conv, #gamma} < 180 cm";
		break;
	case 9:
		RCutString =  "5 cm < R_{conv, #gamma} < 35 cm";
		break;
	default:
		RCutString =  "R cut not defined";
		break;
	}
	return Form("%s, %s", psiPairCut.Data(), RCutString.Data());
}   


TString AnalyseRapidityMesonCut(Int_t RapidityMesonCut){ 
   // Set Cut
	switch(RapidityMesonCut){
	case 0:  //
		return "|y_{meson}| < 0.9";
	case 1:  //
		return "|y_{meson}| < 0.8";
	case 2:  //
		return "|y_{meson}| < 0.7";
	case 3:  //
		return "|y_{meson}| < 0.6";
	case 4:  //
		return "|y_{meson}| < 0.5";
	case 5:  //
		return "|y_{meson}| < 0.65";
	case 6:  //
		return "|y_{meson}| < 0.75";
	case 7:  //
		return "|y_{meson}| < 0.3";
	case 8:  //
		return "|y_{meson}| < 0.35";
	case 9:  //
		return "|y_{meson}| < 0.4";
	
	default:
		return "rapidity cut not defined";
	}
}


TString AnalyseRapidityMesonCutpPb(Int_t RapidityMesonCut){ 
	// Set Cut
	switch(RapidityMesonCut){
	case 0:  //
		return "|y_{meson}| < 0.9";
	case 1:  //
		return "|y_{meson}| < 0.8";
	case 2:  //
		return "|y_{meson}| < 0.7";
	case 3:  //
		return "|y_{meson}| < 0.6";
	case 4:  //
		return "|y_{meson}| < 0.5";
	case 5:  //
		return "|y_{meson}| < 0.65";
	case 6:  //
		return "|y_{meson}| < 0.75";
	case 7:  //
		return "0.165 < y_{c.m.s, meson} < 0.765";
	case 8:  //
		return "|y_{meson}| < 0.35";
	case 9:  //
		return "|y_{meson}| < 0.4";
	
	default:
		return "rapidity cut not defined";
	}
}

TString AnalyseV0ReaderCut(Int_t V0ReaderCut){ 
	// Set Cut
	switch(V0ReaderCut){
	case 0:  //
		return "Onfly V0finder";
	case 1:  //
		return "Offline V0finder";
	default:
		return "V0finder cut not defined";
	}
}

TString AnalyseSpecialTriggerCut(Int_t SpecialTrigger){ 
	// Set Cut
	switch(SpecialTrigger){
	case 0:
		return "without SDD, V0OR";
	case 1:
		return "without SDD, V0AND";
	case 2:
		return "with SDD, V0OR";
	case 3:
		return "with SDD, V0AND";
	default:
		return "special Trigger cut not defined";
	}
}

TString ReturnTextReconstructionProcess(Int_t mode){
	switch (mode){
		case 0:
			return "PCM";
		case 1: 
			return "PCM, Dalitz";
		case 2:
			return "PCM, EMCAL";
		case 3: 
			return "PCM, PHOS";
		case 4:
			return "EMCAL, EMCAL";
		case 5: 
			return "PHOS, PHOS";
		default:
			return "not known";
	}
	
}	


TString ReturnFullTextReconstructionProcess(Int_t mode){ 
	switch (mode){
		case 0:
			return "#gamma's rec. with PCM";
		case 1: 
			return "#gamma's rec. with PCM, Dalitz";
		case 2:
			return "#gamma's rec. with PCM, EMCAL";
		case 3: 
			return "#gamma's rec. with PCM, PHOS";
		case 4:
			return "#gamma's rec. with EMCAL, EMCAL";
		case 5: 
			return "#gamma's rec. with PHOS, PHOS";
		default:
			return "not known";
	}
}

TString ReturnFullTextMeson(TString fEnergyFlagOpt, TString textProcessOpt){ 

	if(fEnergyFlagOpt.CompareTo("7TeV") == 0){
		return Form("pp #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 7 TeV ",textProcessOpt.Data());
	} else if(fEnergyFlagOpt.CompareTo("8TeV") == 0){
		return Form("pp #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 8 TeV ",textProcessOpt.Data());
    } else if( fEnergyFlagOpt.CompareTo("900GeV") == 0) {
		return  Form("pp #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 900 GeV ",textProcessOpt.Data());
	} else if( fEnergyFlagOpt.CompareTo("2.76TeV") == 0) {
		return Form("pp #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 2.76 TeV ",textProcessOpt.Data());
	} else if( fEnergyFlagOpt.CompareTo("PbPb_2.76TeV") == 0) {
		return  Form("PbPb #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 2.76 TeV ",textProcessOpt.Data());
	} else if( fEnergyFlagOpt.CompareTo("pPb_5.023TeV") == 0) {
		return  Form("pPb #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 5.02 TeV ",textProcessOpt.Data());
	} else {
		cout << "No correct collision system specification, has been given" << endl;
		return "";
	}
}

TString ReturnFullTextMesonPiPlPiMiGamma(TString fEnergyFlagOpt, TString textProcessOpt){ 

	if(fEnergyFlagOpt.CompareTo("7TeV") == 0){
		return Form("pp #rightarrow %s (#rightarrow #pi^{+}#pi^{-}#gamma) + X @ 7 TeV ",textProcessOpt.Data());
	} else if( fEnergyFlagOpt.CompareTo("8TeV") == 0) {
		return  Form("pp #rightarrow %s (#rightarrow #pi^{+}#pi^{-}#gamma) + X @ 8 TeV ",textProcessOpt.Data());
	} else if( fEnergyFlagOpt.CompareTo("900GeV") == 0) {
		return  Form("pp #rightarrow %s (#rightarrow #pi^{+}#pi^{-}#gamma) + X @ 900 GeV ",textProcessOpt.Data());
	} else if( fEnergyFlagOpt.CompareTo("2.76TeV") == 0) {
		return Form("pp #rightarrow %s (#rightarrow #pi^{+}#pi^{-}#gamma) + X @ 2.76 TeV ",textProcessOpt.Data());
	} else if( fEnergyFlagOpt.CompareTo("PbPb_2.76TeV") == 0) {
		return  Form("PbPb #rightarrow %s (#rightarrow #pi^{+}#pi^{-}#gamma) + X @ 2.76 TeV ",textProcessOpt.Data());
	} else if( fEnergyFlagOpt.CompareTo("pPb_5.023TeV") == 0) {
		return  Form("pPb #rightarrow %s (#rightarrow #pi^{+}#pi^{-}#gamma) + X @ 5.02 TeV ",textProcessOpt.Data());
	} else {
		cout << "No correct collision system specification, has been given" << endl;
		return "";
	}

}


TString ReturnFullCollisionsSystem(TString fEnergyFlagOpt){ 
	if(fEnergyFlagOpt.CompareTo("7TeV") == 0){
		return  "pp, #sqrt{#it{s}} = 7 TeV";
	} else if( fEnergyFlagOpt.CompareTo("8TeV") == 0) {
		return  "pp, #sqrt{#it{s}} = 8 TeV";	
	} else if( fEnergyFlagOpt.CompareTo("900GeV") == 0) {
		return  "pp, #sqrt{#it{s}} = 900 GeV";
	} else if( fEnergyFlagOpt.CompareTo("2.76TeV") == 0) {
		return  "pp, #sqrt{#it{s}} = 2.76 TeV";
	} else if( (fEnergyFlagOpt.CompareTo("PbPb_2.76TeV") == 0) || (fEnergyFlagOpt.CompareTo("HI") == 0) ) {
		return "Pb-Pb, #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";
	} else if( fEnergyFlagOpt.CompareTo("pPb_5.023TeV") == 0) {
		return "p-Pb, #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";
	} else {
		cout << "No correct collision system specification, has been given" << endl;
		return "";
	}
}

Double_t ReturnCollisionEnergy(TString fEnergyFlagOpt){ 
	if(fEnergyFlagOpt.CompareTo("7TeV") == 0){
		return  7000;
	} else if( fEnergyFlagOpt.CompareTo("8TeV") == 0) {
		return 8000; 	
	} else if( fEnergyFlagOpt.CompareTo("2.76TeV") == 0) {
		return 2760; 
	} else if( fEnergyFlagOpt.CompareTo("900GeV") == 0) {
		return 900;
	} else if( (fEnergyFlagOpt.CompareTo("PbPb_2.76TeV") == 0) || (fEnergyFlagOpt.CompareTo("HI") == 0) ) {
		return 2760;
	} else if( fEnergyFlagOpt.CompareTo("pPb_5.023TeV") == 0) {
		return 5023;
	} else {
		cout << "No correct collision system energy specification, has been given" << endl;
		return 1;     
	}
}

Double_t ReturnCorrectK0ScalingFactor(TString fEnergyFlagOpt, TString cutNr){
	if(fEnergyFlagOpt.CompareTo("7TeV") == 0){
		return 1./0.75 -1.;
	} else if( fEnergyFlagOpt.CompareTo("8TeV") == 0) {
		return  1./0.75 -1.;		
	} else if( fEnergyFlagOpt.CompareTo("2.76TeV") == 0) {
		return  1./0.685 -1.;
	} else if( fEnergyFlagOpt.CompareTo("900GeV") == 0) {
		return 1./0.6 -1.;
	} else if( (fEnergyFlagOpt.CompareTo("PbPb_2.76TeV") == 0) || (fEnergyFlagOpt.CompareTo("HI") == 0) ) {
		return GetScalingFactorSecCorrection(cutNr.Data());
	} else if( fEnergyFlagOpt.CompareTo("pPb_5.023TeV") == 0) {
		// return 0.;
		return  GetScalingFactorSecCorrection(cutNr.Data());
	} else {
		cout << "No correct collision system specification, has been given" << endl;
		return 0.;     
	}
}

TString AnalyseChi2MesonCut(Int_t chi2GammaCut){   // Set Cut
	switch(chi2GammaCut){
		case 0: // 100
			return "#chi_{meson}^{2} < 100";
		case 1:  // 50
			return "#chi_{meson}^{2} < 50";
		case 2:  // 30
			return "#chi_{meson}^{2} < 30";
		case 3:
			return "#chi_{meson}^{2} < 200";
		case 4:
			return "#chi_{meson}^{2} < 500";
		case 5:
			return "#chi_{meson}^{2} < 1000";  
		default:
			return "#chi_{meson}^{2} cut unknown";
	}
}

TString AnalyseChi2PsiPair(Int_t PsiPairCut ){
	switch(PsiPairCut) {
		case 26:
			return "2D: #chi_{#gamma}^{2} = 30 |#Psi_{Pair}| < 0.05";
		case 22:
			return "1D: #chi_{#gamma}^{2} = 30 |#Psi_{Pair}| < 0.05";
		case 16:
			return "2D: #chi_{#gamma}^{2} = 50 |#Psi_{Pair}| < 0.05";
		case 86:
			return "2D: #chi_{#gamma}^{2} = 20|#Psi_{Pair}| < 0.05";
		case 25:
			return "2D: #chi_{#gamma}^{2} = 30 |#Psi_{Pair}| < 0.1";
		case 27:
			return "2D: #chi_{#gamma}^{2} = 30 |#Psi_{Pair}| < 0.035";
		default:
			return " #chi_{#gamma}^{2} #Psi_{Pair} cut not defined";
	}
}  

TString AnalysePhotonQuality(Int_t photonQualitCut ){
	switch(photonQualitCut) {
		case 0:
			return "photon Quality: 1,2,3";
		case 2:
			return "photon Quality: 1 (TPC only photons)";
		case 3:
			return "photon Quality: 2 (1 leg with #geq 2 ITS hits)";
		case 4:
			return "photon Quality: 3 (both legs with #geq 2 ITS hits)";
		default:
			return "photon Quality cut not defined";
	}
}  


TString ReturnMesonString (TString mesonName){
	if(mesonName.CompareTo("Pi0") == 0 || mesonName.CompareTo("Pi0EtaBinning") == 0){
		return "#pi^{0}";
	} else if(mesonName.CompareTo("Eta") == 0)  {
		return "#eta";
	} else if(mesonName.CompareTo("Omega") == 0)  {
		return "#omega";
	} else if(mesonName.CompareTo("EtaPrim") == 0)  {
		return "#eta'";
	} else {
		cout << "No correct meson has been selected" << endl;
		return "";
	}
}

Bool_t ReturnMesonOption (TString mesonName){
	if(mesonName.CompareTo("Pi0") == 0 || mesonName.CompareTo("Pi0EtaBinning") == 0){
		return kTRUE;
	} else if(mesonName.CompareTo("Eta") == 0)  {
		return kFALSE;
	} else {
		return kFALSE;
	}
}

Double_t BWdndptPi0(Double_t *x, Double_t *par) {
	
	Double_t pt = x[0]*1000.; //GeV->MeV
	Double_t t=par[1];   
	Double_t beta=abs(par[2]);
	Double_t yt=0.5*TMath::Log((1+beta)/(1-beta)); // atanh(beta)
	Double_t m0=135.; //pi0
	Double_t mt=TMath::Sqrt(m0*m0+pt*pt);
	
	Double_t f = par[0]*pt*mt*TMath::BesselI(0,pt*sinh(yt)/t)*TMath::BesselK(1,TMath::Abs(mt*cosh(yt)/t));
	
	return f;
}

Double_t BWdndptEta(Double_t *x, Double_t *par) {

	Double_t pt = x[0]*1000.; //GeV->MeV
	Double_t t=par[1];
	Double_t beta=abs(par[2]);
	Double_t yt=0.5*TMath::Log((1+beta)/(1-beta)); // atanh(beta)
	Double_t m0=548.; //eta
	Double_t mt=TMath::Sqrt(m0*m0+pt*pt);
	
	Double_t f = par[0]*pt*mt*TMath::BesselI(0,pt*sinh(yt)/t)*TMath::BesselK(1,TMath::Abs(mt*cosh(yt)/t));
	
	return f;
}

void CorrectGammaDataUnfold(TH1D* histoGammaSpecCorr,TH1D* histoConvProb, TH1D* histoRecoEff, Double_t deltaEtaDummy, Double_t scalingDummy, Double_t nEvtDummy){
	histoGammaSpecCorr->Divide(histoGammaSpecCorr,histoConvProb,1.,1.,"");
	histoGammaSpecCorr->Divide(histoGammaSpecCorr,histoRecoEff,1.,1.,"");
	histoGammaSpecCorr->Scale(1./deltaEtaDummy);
	histoGammaSpecCorr->Scale(scalingDummy);
	histoGammaSpecCorr->Scale(1./nEvtDummy);
	for (Int_t i = 1; i < histoGammaSpecCorr->GetNbinsX()+1 ; i++){
		Double_t newBinContent = histoGammaSpecCorr->GetBinContent(i)/histoGammaSpecCorr->GetBinCenter(i);
		Double_t newBinError = histoGammaSpecCorr->GetBinError(i)/histoGammaSpecCorr->GetBinCenter(i);
		histoGammaSpecCorr->SetBinContent(i,newBinContent);
		histoGammaSpecCorr->SetBinError(i,newBinError);
	}
}

TH1D *FixEfficiency(TH1D *FixedEff, TH1D* Eff, TString option, TString centString){
	if(option.CompareTo("2.76TeV") == 0){
		FixedEff->SetBinContent(8,(Eff->GetBinContent(7)+Eff->GetBinContent(9))/2);
		FixedEff->SetBinContent(13,(Eff->GetBinContent(12)+Eff->GetBinContent(14))/2);
	}
	
	if(centString.CompareTo("0-5%") == 0){
		FixedEff->SetBinContent(7,(Eff->GetBinContent(6)+Eff->GetBinContent(8))/2);
		FixedEff->SetBinContent(10,(Eff->GetBinContent(9)+Eff->GetBinContent(11))/2);
		FixedEff->SetBinContent(12,(Eff->GetBinContent(11)+Eff->GetBinContent(13))/2);
		FixedEff->SetBinContent(15,(Eff->GetBinContent(14)+Eff->GetBinContent(16))/2);
	}
	if(centString.CompareTo("5-10%") == 0){
		FixedEff->SetBinContent(15,(Eff->GetBinContent(14)+Eff->GetBinContent(16))/2);
		FixedEff->SetBinContent(18,(Eff->GetBinContent(17)+Eff->GetBinContent(19))/2);
	}
	if(centString.CompareTo("0-10%") == 0){
		FixedEff->SetBinContent(10,(Eff->GetBinContent(9)+Eff->GetBinContent(11))/2);
		FixedEff->SetBinContent(13,(Eff->GetBinContent(12)+Eff->GetBinContent(14))/2);
		FixedEff->SetBinContent(16,(Eff->GetBinContent(15)+Eff->GetBinContent(17))/2);
	}
	if(centString.CompareTo("0-20%") == 0){
		FixedEff->SetBinContent(10,(Eff->GetBinContent(9)+Eff->GetBinContent(11))/2);
		FixedEff->SetBinContent(15,(Eff->GetBinContent(14)+Eff->GetBinContent(16))/2);
	}
	if(centString.CompareTo("0-80%") == 0){
		FixedEff->SetBinContent(11,(Eff->GetBinContent(10)+Eff->GetBinContent(12))/2);
		FixedEff->SetBinContent(14,(Eff->GetBinContent(13)+Eff->GetBinContent(15))/2);
	}
	if(centString.CompareTo("10-20%") == 0){
		FixedEff->SetBinContent(10,(Eff->GetBinContent(9)+Eff->GetBinContent(11))/2);
		FixedEff->SetBinContent(13,(Eff->GetBinContent(12)+Eff->GetBinContent(16))/2);
		FixedEff->SetBinContent(14,(FixedEff->GetBinContent(13)+Eff->GetBinContent(16))/2);
		FixedEff->SetBinContent(15,(FixedEff->GetBinContent(14)+Eff->GetBinContent(16))/2);
		FixedEff->SetBinContent(17,(Eff->GetBinContent(16)+Eff->GetBinContent(18))/2);
	}
	if(centString.CompareTo("20-30%") == 0){
		FixedEff->SetBinContent(14,(Eff->GetBinContent(13)+Eff->GetBinContent(15))/2);
		FixedEff->SetBinContent(16,(Eff->GetBinContent(15)+Eff->GetBinContent(17))/2);
	}
	if(centString.CompareTo("20-40%") == 0){
		FixedEff->SetBinContent(12,(Eff->GetBinContent(11)+Eff->GetBinContent(13))/2);
		FixedEff->SetBinContent(14,(Eff->GetBinContent(13)+Eff->GetBinContent(15))/2);
		FixedEff->SetBinContent(16,(Eff->GetBinContent(15)+Eff->GetBinContent(17))/2);
	}
	if(centString.CompareTo("20-50%") == 0){
		FixedEff->SetBinContent(12,(Eff->GetBinContent(11)+Eff->GetBinContent(13))/2);
		FixedEff->SetBinContent(14,(Eff->GetBinContent(13)+Eff->GetBinContent(15))/2);
		FixedEff->SetBinContent(16,(Eff->GetBinContent(15)+Eff->GetBinContent(17))/2);
	}
	if(centString.CompareTo("30-50%") == 0){
		FixedEff->SetBinContent(11,(Eff->GetBinContent(10)+Eff->GetBinContent(12))/2);
		FixedEff->SetBinContent(16,(Eff->GetBinContent(15)+Eff->GetBinContent(17))/2);
	}
	if(centString.CompareTo("40-60%") == 0){
		FixedEff->SetBinContent(8,(Eff->GetBinContent(7)+Eff->GetBinContent(9))/2);
		FixedEff->SetBinContent(12,(FixedEff->GetBinContent(13)+Eff->GetBinContent(16))/2);
		FixedEff->SetBinContent(13,(FixedEff->GetBinContent(12)+Eff->GetBinContent(16))/2);
		FixedEff->SetBinContent(14,(FixedEff->GetBinContent(13)+Eff->GetBinContent(16))/2);
		FixedEff->SetBinContent(15,(FixedEff->GetBinContent(14)+Eff->GetBinContent(16))/2);
		FixedEff->SetBinContent(17,(Eff->GetBinContent(16)+Eff->GetBinContent(18))/2);
	}
	if(centString.CompareTo("50-80%") == 0){
		FixedEff->SetBinContent(11,(Eff->GetBinContent(10)+Eff->GetBinContent(12))/2);
		FixedEff->SetBinContent(10,(Eff->GetBinContent(9)+FixedEff->GetBinContent(11))/2);
		FixedEff->SetBinContent(11,(FixedEff->GetBinContent(10)+Eff->GetBinContent(12))/2);
		FixedEff->SetBinContent(15,(Eff->GetBinContent(14)+Eff->GetBinContent(16))/2);
	}
	if(centString.CompareTo("60-80%") == 0){
		FixedEff->SetBinContent(11,(Eff->GetBinContent(10)+Eff->GetBinContent(12))/2);
		FixedEff->SetBinContent(10,(Eff->GetBinContent(9)+FixedEff->GetBinContent(11))/2);
		FixedEff->SetBinContent(11,(FixedEff->GetBinContent(10)+Eff->GetBinContent(12))/2);
		FixedEff->SetBinContent(10,(Eff->GetBinContent(9)+FixedEff->GetBinContent(11))/2);
	}
	if(centString.CompareTo("40-80%") == 0){
		FixedEff->SetBinContent(11,(Eff->GetBinContent(10)+Eff->GetBinContent(12))/2);
		FixedEff->SetBinContent(14,(Eff->GetBinContent(13)+Eff->GetBinContent(15))/2);
	}

	return FixedEff;


}


Int_t CalcDeltaPtOverPt(const TGraphAsymmErrors* gScaledPPStatErr, const TGraphAsymmErrors* gScaledPPSysErr, 
         const TGraphAsymmErrors* gPbPbStatErr, const TGraphAsymmErrors* gPbPbSysErr, 
         TGraphAsymmErrors* gResultStatErr, TGraphAsymmErrors* gResultSysErr, 
         bool bResultAtProtonProtonPt = true) {

	//
	// Klaus Reygers, 15-Sep-2013
	//

	//
	// Calculate relative pT shift between the Ncoll scaled pp invariant neutral pion yield and the Pb+Pb inv. yield.
	// For each Pb+Pb data points the pT,pp is calculated at which the saled pp inv. yield has the same value
	// as the Pb+Pb yield. The result of this routine is the normalized difference (pT - pT,pp)/pT,pp.
	// which cann be seen a measure of the fractional energy loss of a parton.
	// By default this result is given at pT,pp. In case of bResultAtProtonProtonPt == false it is given at the A+A pT.
	// See arXiv:1208.2254 for details.
	//
	// Klaus Reygers, 15-Sep-2013
	//
		Int_t iStatus = 0;

	// define Tsallis function to parameterize scaled pp invariant yield
	Double_t mass = 0.135; // pi0 mass
	TF1* fTsallis = new TF1("Tsallis", Form("[0] * pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1])",
					mass, mass, mass),0.,30.);
	fTsallis->SetParNames("A","n","T_{Tsallis} (GeV/c)");
	fTsallis->SetParameters(1.,5.,0.18);

	//
	// fit scaled pp invariant yield
	//

	// copy scaled pp spectrum as fit would modify it (not allows due to 'const')
	TGraphAsymmErrors* gScaledPPStatErrCopy = dynamic_cast<TGraphAsymmErrors*>(gScaledPPStatErr->Clone());

	gScaledPPStatErrCopy->Fit(fTsallis, "0qN");
	Double_t prob = fTsallis->GetProb();
	Double_t p_A = fTsallis->GetParameter(0);
	Double_t p_n = fTsallis->GetParameter(1);
	Double_t p_T = fTsallis->GetParameter(2);
	delete gScaledPPStatErrCopy;

	// check fit status
	if (prob < 0.05) { 
		cout << "CalcDeltaPtOverPt: Warning: bad fit of pp spectrum (p-Value = " << prob << ")" << endl;
		iStatus = 1;
	}

	// function used to determine value of pTppStar at which the pp yield (pT * dN/dpT = pT^2 * inv. yield) 
	// takes a given value yieldPbPb
	TF1* fdNdPtScaledPP = new TF1("fdNdPtScaledPP", 
				Form("[0] * x * x * pow(1.+(sqrt(x*x+%.10f*%.10f)-%.10f)/([1]*[2]), -[1]) - [3]",
					mass, mass, mass),0.,30.);
	fdNdPtScaledPP->SetParNames("A","n","T_{Tsallis} (GeV/c)","yieldPbPb");


	// create root finder
	ROOT::Math::BrentRootFinder brf;

	//
	// loop over Pb+Pb data point and calculate pTppStar - pT for each Pb+Pb point with a given pT 
	//
	for (Int_t i=0; i<gPbPbStatErr->GetN(); i++) {

		Double_t pTPbPb = gPbPbStatErr->GetX()[i];
		Double_t sfAA = pTPbPb * pTPbPb;
		Double_t yPbPb = sfAA * gPbPbStatErr->GetY()[i];  // yield pT * dN/dpT = pT^2 * inv. yield
		Double_t yPbPbStatErr = sfAA * gPbPbStatErr->GetEYhigh()[i];
		Double_t yPbPbSysErr = sfAA * gPbPbSysErr->GetEYhigh()[i];

		// determine relative errors of pp yields at the same pT
		Double_t pTpp = gScaledPPStatErr->GetX()[i];
		Double_t sfpp = pTpp * pTpp;
		Double_t ypp = sfpp * gScaledPPStatErr->GetY()[i]; // yield pT * dN/dpT = pT^2 * inv. yield
		Double_t yppStatErr = sfpp * gScaledPPStatErr->GetEYhigh()[i];
		Double_t yppRelStatErr = yppStatErr/ypp;
		Double_t yppSysErr = sfpp * gScaledPPSysErr->GetEYhigh()[i];
		Double_t yppRelSysErr = yppSysErr/ypp;

		// combine the error of the Pb+Pb yield and the pp yield at pT,pp assuming that one can neglect the
		// the change in the relative error of the pp yield between pT and pT,pp
		Double_t yCombinedStatErr = sqrt(TMath::Power(yPbPbStatErr, 2) + TMath::Power(yPbPb * yppRelStatErr, 2));
		Double_t yCombinedSysErr = sqrt(TMath::Power(yPbPbSysErr, 2) + TMath::Power(yPbPb * yppRelSysErr, 2));

		//
		// calculate the pTppStar at which the scaled pp spectrum has the value of the PbPb yield
		//

		// create wrapper function
		fdNdPtScaledPP->SetParameters(p_A, p_n, p_T, yPbPb);
		ROOT::Math::WrappedTF1 wf(*fdNdPtScaledPP);

		brf.SetFunction(wf, pTPbPb, 20.);
		brf.SetNpx(1000);
		brf.Solve();
		Double_t pTppStar = brf.Root(); // pT at which scaled pp yield takes value of Pb+Pb yield
		
		Double_t DeltaPt = pTppStar - pTPbPb;

		Double_t inverseSlope = -1./fdNdPtScaledPP->Derivative(pTppStar);

		Double_t pTppStatErr = inverseSlope * yCombinedStatErr;
		Double_t pTppSysErr = inverseSlope * yCombinedSysErr;

		// sanity check
		if (fabs(fdNdPtScaledPP->Eval(pTppStar)) > 1e-5) {
		cout << "CalcDeltaPtOverPt: Warning: root calculation incorrect: " << fabs(fdNdPtScaledPP->Eval(pTppStar)) << endl;
		iStatus = 1;
		}

		// determine at which pT the result is shown
		Double_t pTPlot = pTppStar;  // show result at p+p pT
		if (!bResultAtProtonProtonPt) pTPlot = pTPbPb;  // show result at A+A pT

		// result with statistical error
		gResultStatErr->GetX()[i] = pTPlot;
		gResultStatErr->GetY()[i] = DeltaPt/pTppStar;
		gResultStatErr->GetEXhigh()[i] = 0;
		gResultStatErr->GetEXlow()[i] = 0;
		gResultStatErr->GetEYhigh()[i] = pTPbPb/(pTppStar*pTppStar) * pTppStatErr;  // Gaussian error propagation
		gResultStatErr->GetEYlow()[i] = pTPbPb/(pTppStar*pTppStar) * pTppStatErr;   // Gaussian error propagation

		// result with systematic error
		gResultSysErr->GetX()[i] = pTPlot;
		gResultSysErr->GetY()[i] = DeltaPt/pTppStar;
		gResultSysErr->GetEXhigh()[i] = 0.2;  // to make error boxes visible
		gResultSysErr->GetEXlow()[i] = 0.2;  // to make error boxes visible
		gResultSysErr->GetEYhigh()[i] = pTPbPb/(pTppStar*pTppStar) * pTppSysErr;  // Gaussian error propagation
		gResultSysErr->GetEYlow()[i] = pTPbPb/(pTppStar*pTppStar) * pTppSysErr;   // Gaussian error propagation

	}

	delete fTsallis;
	delete fdNdPtScaledPP;

	return iStatus;

}
