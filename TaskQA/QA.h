/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Daniel Muehlheim, d.muehlheim@cern.ch                          *****
 *******************************************************************************/

#ifndef QA_H
#define QA_H

//**************************************************************************************************************************
//************************************************* QA header **************************************************************
//**************************************************************************************************************************
//**************************************************************************************************************************

using namespace std;

#include <Riostream.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <math.h>
#include <queue>
#include <algorithm>

#include <TApplication.h>
#include <TArrow.h>
#include <TASImage.h>
#include <TAxis.h>
#include <TAttAxis.h>
#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TFrame.h>
#include <TF1.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <THnSparse.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TList.h>
#include <TMath.h>
#include <TMarker.h>
#include <TMinuit.h>
#include <TMultiGraph.h>
#include <TObject.h>
#include <TPaletteAxis.h>
#include <TPaveLabel.h>
#include <TPostScript.h>
#include <TProfile2D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TString.h>
#include <TSystem.h>
#include <TVirtualFitter.h>

#include "../CommonHeaders/PlottingMeson.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/AdjustHistRange.h"
// #include "../CommonHeaders/ExtractSignalBinning.h"
// #include "../CommonHeaders/ExtractSignalPlotting.h"

//**************************************************************************************************************************
Double_t 	fFWHMFunc;
Double_t 	fFWHMFuncError;

struct CellQAObj{
	CellQAObj(){
		cellIDsEnergy.clear();
		cellIDsTime.clear();
		cellIDsHotCells1D.clear();
		cellIDsHotCellsTime1D.clear();
		cellIDsHotCells2D.clear();
		cellIDsMissing.clear();
		goodCells.clear();
	}
	Double_t EnergyMean[2];
	Double_t EnergySigma[2];
	Double_t TimeMean[2];
	Double_t TimeSigma[2];
	Double_t HotCells1D[2];
	Double_t HotCellsTime1D[2];
	Double_t HotCells2D[9][2];
	std::vector<Int_t> cellIDsEnergy;
	std::vector<Int_t> cellIDsTime;
	std::vector<Int_t> cellIDsHotCells1D;
	std::vector<Int_t> cellIDsHotCellsTime1D;
	std::vector<Int_t> cellIDsHotCells2D;
	std::vector<Int_t> cellIDsMissing;
	std::vector<Int_t> goodCells;
};

void CalculateFWHM(TF1 * fFunc, Double_t startMass, Double_t endMass);

class MesonFit{

public:
	MesonFit(){
		ResetPointer();
		Init();
		return;
	}

	~MesonFit(){
		Init();
		return;
	}

	TH1D* GetGammaGamma(){return fGammaGamma;}
	TH1D* GetBck(){return fBck;}
	TH1D* GetBckNormPi0(){return fBckNormPi;}
	TH1D* GetBckNormEta(){return fBckNormEta;}
	TH1D* GetSignalPi0(){return fSignalPi;}
	TH1D* GetSignalEta(){return fSignalEta;}
	TF1* GetFitPi0(){return fFitRecoPi;}
	TF1* GetFitEta(){return fFitRecoEta;}

	void GetMeson(Bool_t isPi0, Double_t &width,Double_t &widthErr,Double_t &mass, Double_t &massErr){
		if(isPi0){
            width = widthPi;
            widthErr = widthPiErr;
			mass = massPi;
			massErr = massPiErr;
		}else{
            width = widthEta;
            widthErr = widthEtaErr;
			mass = massEta;
			massErr = massEtaErr;
		}
		return;
	}

	void GetMesonRatios(Bool_t isPi0, Double_t &ratio,Double_t &ratioErr){
		if(isPi0){
			ratio = ratioPi0;
			ratioErr = ratioPi0Err;
		}else{
			ratio = ratioEta;
			ratioErr = ratioEtaErr;
		}
		return;
	}

    Bool_t DoFitting(TH2D* fInvMassMesonPt, TH2D* fInvMassBGPt, Double_t nEventsBin1, Int_t mode, TString saveCanvas, TString name, Bool_t doPrint, Bool_t doLog, fstream& fLog){
		Init();

        //Fitting starts at pT = 1 GeV/c
        fGammaGamma = (TH1D*)fInvMassMesonPt->ProjectionX("ESD_Mother_InvMass",10,fInvMassMesonPt->GetNbinsY());
        fBck = (TH1D*)fInvMassBGPt->ProjectionX("ESD_BG_InvMass",10,fInvMassBGPt->GetNbinsY());

		for (Int_t binx= 0; binx < fGammaGamma->GetNbinsX()+1; binx++){
			if(fGammaGamma->GetBinContent(binx) == 0){
				fGammaGamma->SetBinError(binx,1.);
				fGammaGamma->SetBinContent(binx,0.);
			}
		}
		fBckNormPi = (TH1D*)fBck->Clone("fBckNormPi");
		fBckNormEta = (TH1D*)fBck->Clone("fBckNormEta");
		fGammaGamma->Sumw2();
		fBck->Sumw2();
		fBckNormPi->Sumw2();
		fBckNormEta->Sumw2();

		Double_t    rPi= fGammaGamma->Integral(fGammaGamma->GetXaxis()->FindBin(0.17),fGammaGamma->GetXaxis()->FindBin(0.3));
		Double_t    bPi= fBck->Integral(fBck->GetXaxis()->FindBin(0.17),fBck->GetXaxis()->FindBin(0.3));
		Double_t    normPi = 1;

		Double_t    rEta= fGammaGamma->Integral(fGammaGamma->GetXaxis()->FindBin(0.58),fGammaGamma->GetXaxis()->FindBin(0.8));
		Double_t    bEta= fBck->Integral(fBck->GetXaxis()->FindBin(0.58),fBck->GetXaxis()->FindBin(0.8));
		Double_t    normEta = 1;

		if(bPi != 0) normPi = rPi/bPi;
		fBckNormPi->Scale(normPi);
		if(bEta != 0) normEta = rEta/bEta;
		fBckNormEta->Scale(normEta);

		Int_t numberOfZeros = 0;
		for (Int_t i = 1; i < fBckNormPi->GetNbinsX()+1; i++){
			if (fBckNormPi->GetBinContent(i) == 0){
				numberOfZeros++;
				if (normPi > 1.){
					fBckNormPi->SetBinError(i,1.);
					fBckNormPi->SetBinContent(i,0.);
				}
			}
        }
		numberOfZeros = 0;
		for (Int_t i = 1; i < fBckNormEta->GetNbinsX()+1; i++){
			if (fBckNormEta->GetBinContent(i) == 0){
				numberOfZeros++;
				if (normEta > 1.){
					fBckNormEta->SetBinError(i,1.);
					fBckNormEta->SetBinContent(i,0.);
				}
			}
		}

		TCanvas* canvasMass = new TCanvas("canvasMass","",200,10,1350,900);  // gives the page size
		fSignalPi = (TH1D*)fGammaGamma->Clone("fSignalPi");
		fSignalPi->SetTitle("");
		fSignalPi->Sumw2();
		fSignalPi->Add(fBckNormPi,-1.);
        if(fSignalPi->Integral(fSignalPi->FindBin(0.1),fSignalPi->FindBin(0.17))<=50){
          cout << "Total integral <= 10 for pi0 in " << name.Data() << ", skipping fit and returning..." << endl;
          return kFALSE;
        }

		if (mode == 0){
			fSignalPi->Rebin(2);
		} else {
			fSignalPi->Rebin(4);
		}
        fSignalPi->GetXaxis()->SetRangeUser( 0.,0.3);

		Double_t mesonAmplitude =fSignalPi->GetMaximum();
		Double_t mesonAmplitudeMinPi  = mesonAmplitude*80./100.;
		Double_t mesonAmplitudeMaxPi = mesonAmplitude*115./100.;
        fFitRecoPi = new TF1("GaussExpLinear","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)",0.05,0.25);
		fFitRecoPi->SetParameter(0,mesonAmplitude);
		fFitRecoPi->SetParameter(1,0.135);
		fFitRecoPi->SetParLimits(0,mesonAmplitudeMinPi,mesonAmplitudeMaxPi);
        fFitRecoPi->SetParLimits(1,0.1,0.2);

        if(mode==2){
          fFitRecoPi->SetParameter(2,0.008);
          fFitRecoPi->SetParameter(3,0.015);
          fFitRecoPi->SetParLimits(2,0.005,0.02);
          fFitRecoPi->SetParLimits(3,0.010,0.04);
        }else if(mode==4){
          fFitRecoPi->SetParameter(2,0.012);
          fFitRecoPi->SetParameter(3,0.020);
          fFitRecoPi->SetParLimits(2,0.01,0.03);
          fFitRecoPi->SetParLimits(3,0.01,0.04);
        }else{
          fFitRecoPi->SetParameter(2,0.003);
          fFitRecoPi->SetParameter(3,0.020);
          fFitRecoPi->SetParLimits(2,0.001,0.01);
          fFitRecoPi->SetParLimits(3,0.001,0.05);
        }
        fSignalPi->Draw();
        fSignalPi->Fit(fFitRecoPi,"SINRMQEC+","",0.05,0.25);
        TFitResultPtr result = fSignalPi->Fit(fFitRecoPi,"SINRMQEC+","",0.05,0.25);

		fFitRecoPi->SetLineColor(3);
		fFitRecoPi->SetLineWidth(1);
		fFitRecoPi->SetLineStyle(1);
		fFitRecoPi->Draw("same");
		canvasMass->SaveAs(Form("%s/Pi0_%s.eps",saveCanvas.Data(),name.Data()));

		Double_t integral = 0;
		Double_t integralErr = 0;
		Double_t integralBG = 0;

		if(doPrint) cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		if(doLog) fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		if(!(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0) ){
			if(doPrint) cout << "Pi0 Fitting failed in with status " << gMinuit->fCstatu.Data() <<endl << endl;
			if(doLog) fLog << "Pi0 Fitting failed in with status " << gMinuit->fCstatu.Data() <<endl << endl;
			TString status = gMinuit->fCstatu.Data();
			if(doPrint) cout << status.Data() << endl;
			if (status.Contains("PROBLEMS")){
                CalculateFWHM(fFitRecoPi,0.1,0.15);
				widthPi = fFWHMFunc;
				widthPiErr = fFWHMFuncError;
				massPi = fFitRecoPi->GetParameter(1);
				massPiErr = fFitRecoPi->GetParError(1);
                integral = fFitRecoPi->Integral(0.1, 0.15, result->GetParams()) / fSignalPi->GetBinWidth(10);
                integralErr = fFitRecoPi->IntegralError(0.1, 0.15, result->GetParams(), result->GetCovarianceMatrix().GetMatrixArray() ) / fSignalPi->GetBinWidth(10);
                integralBG = fFitRecoPi->GetParameter(4)*0.15 + fFitRecoPi->GetParameter(5)/2 *0.15*0.15- (fFitRecoPi->GetParameter(4)*0.1 + fFitRecoPi->GetParameter(5)/2 *0.1*0.1);
                if(doPrint) cout << "Pi0 full width: "  << widthPi << "\t +-" << widthPiErr << "\t Mass: "<< massPi << "\t+-" << massPiErr  << endl;
				if(doPrint) cout << "integral Pi: " << integral << "\t +-" << integralErr << "\t integral BG : "<< integralBG<< endl;
                if(doLog) fLog << "Pi0 full width: "  << widthPi << "\t +-" << widthPiErr << "\t Mass: "<< massPi << "\t+-" << massPiErr  << endl;
				if(doLog) fLog << "integral Pi: " << integral << "\t +-" << integralErr << "\t integral BG : "<< integralBG<< endl;
			} else {
				if(doPrint) cout << "here" << endl;
			}
		} else {
            CalculateFWHM(fFitRecoPi,0.1,0.17);
			widthPi = fFWHMFunc;
			widthPiErr = fFWHMFuncError;
			massPi = fFitRecoPi->GetParameter(1);
			massPiErr = fFitRecoPi->GetParError(1);
            integral = fFitRecoPi->Integral(0.1, 0.15, result->GetParams()) / fSignalPi->GetBinWidth(10);
            integralErr = fFitRecoPi->IntegralError(0.1, 0.15, result->GetParams(), result->GetCovarianceMatrix().GetMatrixArray() ) / fSignalPi->GetBinWidth(10);
            integralBG = fFitRecoPi->GetParameter(4)*0.15 + fFitRecoPi->GetParameter(5)/2 *0.15*0.15- (fFitRecoPi->GetParameter(4)*0.1 + fFitRecoPi->GetParameter(5)/2 *0.1*0.1);
            if(doPrint) cout << "Pi0 full width: "  << widthPi << "\t +-" << widthPiErr << "\t Mass: "<< massPi << "\t+-" << massPiErr  << endl;
			if(doPrint) cout << "integral Pi: " << integral << "\t +-" << integralErr << "\t integral BG : "<< integralBG<< endl;
            if(doLog) fLog << "Pi0 full width: "  << widthPi << "\t +-" << widthPiErr << "\t Mass: "<< massPi << "\t+-" << massPiErr  << endl;
			if(doLog) fLog << "integral Pi: " << integral << "\t +-" << integralErr << "\t integral BG : "<< integralBG<< endl;
		}

		fSignalEta = (TH1D*)fGammaGamma->Clone("fSignalEta");
		fSignalEta->SetTitle("");
		fSignalEta->Sumw2();
		fSignalEta->Add(fBckNormEta,-1.);
        if(fSignalEta->Integral(fSignalEta->FindBin(0.5),fSignalEta->FindBin(0.57))<=25){
          cout << "Total integral <= 5 for eta in " << name.Data() << ", skipping fit and returning..." << endl;
          return kFALSE;
        }


		fSignalEta->Rebin(4);
        fSignalEta->GetXaxis()->SetRangeUser(0.4,0.7);


        fFitRecoEta = new TF1("GaussExpLinearEta","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)",0.45,0.65);
		Double_t mesonAmplitudeEta =fSignalEta->GetMaximum();
		Double_t mesonAmplitudeMinEta  = mesonAmplitudeEta*80./110.;
		Double_t mesonAmplitudeMaxEta = mesonAmplitudeEta*115./100.;

		fFitRecoEta->SetParameter(0,mesonAmplitudeEta);
        fFitRecoEta->SetParameter(1,0.540);
		fFitRecoEta->SetParLimits(0,mesonAmplitudeMinEta,mesonAmplitudeMaxEta);
        fFitRecoEta->SetParLimits(1,0.45,0.60);

        if(mode==2){
          fFitRecoEta->SetParameter(2,0.020);
          fFitRecoEta->SetParameter(3,0.020);
          fFitRecoEta->SetParLimits(2,0.010,0.040);
          fFitRecoEta->SetParLimits(3,0.010,0.030);
        }else if(mode==4){
          fFitRecoEta->SetParameter(2,0.030);
          fFitRecoEta->SetParameter(3,0.025);
          fFitRecoEta->SetParLimits(2,0.015,0.050);
          fFitRecoEta->SetParLimits(3,0.010,0.040);
        }else{
          fFitRecoEta->SetParameter(2,0.010);
          fFitRecoEta->SetParameter(3,0.007);
          fFitRecoEta->SetParLimits(2,0.002,0.050);
          fFitRecoEta->SetParLimits(3,0.005,0.026);
        }

		fSignalEta->Draw("");
        fSignalEta->Fit(fFitRecoEta,"SINRMQEC+","",0.45,0.65);
        TFitResultPtr resultEta = fSignalEta->Fit(fFitRecoEta,"SINRQMEC+","",0.45,0.65);

		fFitRecoEta->SetLineColor(3);
		fFitRecoEta->SetLineWidth(1);
		fFitRecoEta->SetLineStyle(1);
		fFitRecoEta->Draw("same");
		canvasMass->SaveAs(Form("%s/Eta_%s.eps",saveCanvas.Data(),name.Data()));

		Double_t integralEta = 0;
		Double_t integralEtaErr = 0;
		Double_t integralEtaBG = 0;

		if(doPrint) cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		if(doLog) fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		if(!(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0) ){
			if(doPrint) cout << "Fitting failed in with status " << gMinuit->fCstatu.Data() <<endl << endl;
			if(doLog) fLog << "Fitting failed in with status " << gMinuit->fCstatu.Data() <<endl << endl;
		} else {
			CalculateFWHM(fFitRecoEta,0.5,0.57);
			widthEta = fFWHMFunc;
			widthEtaErr = fFWHMFuncError;
			massEta = fFitRecoEta->GetParameter(1);
			massEtaErr = fFitRecoEta->GetParError(1);
            if(doPrint) cout << "Eta full width: "  << widthEta << "\t +-" << widthEtaErr << "\t Mass: "<< massEta << "\t+-" << massEtaErr  << endl;
            if(doLog) fLog << "Eta full width: "  << widthEta << "\t +-" << widthEtaErr << "\t Mass: "<< massEta << "\t+-" << massEtaErr  << endl;

			integralEta = fFitRecoEta->Integral(0.5,0.57, resultEta->GetParams()) / fSignalEta->GetBinWidth(10);
			integralEtaErr = fFitRecoEta->IntegralError(0.5,0.57, resultEta->GetParams(), resultEta->GetCovarianceMatrix().GetMatrixArray() ) / fSignalEta->GetBinWidth(10);
			integralEtaBG = fFitRecoEta->GetParameter(4)*0.57 + fFitRecoEta->GetParameter(5)/2 *0.57*0.57- (fFitRecoEta->GetParameter(4)*0.5 + fFitRecoEta->GetParameter(5)/2 *0.5*0.5);
			if(doPrint) cout << "integral Eta: " << integralEta << "\t +-" << integralEtaErr << "\t integral BG : "<< integralEtaBG<< endl;
			if(doLog) fLog << "integral Eta: " << integralEta << "\t +-" << integralEtaErr << "\t integral BG : "<< integralEtaBG<< endl;
		}

		ratioPi0 = (Double_t)	(integral-integralBG)/nEventsBin1;
		ratioPi0Err = sqrt( pow(integralErr/nEventsBin1,2)  + pow( sqrt(nEventsBin1)*(integral-integralBG)/pow(nEventsBin1,2),2) );
		ratioEta = (Double_t)	(integralEta-integralEtaBG)/nEventsBin1;
		ratioEtaErr = sqrt( pow(integralEtaErr/nEventsBin1,2)  + pow( sqrt(nEventsBin1)*(integralEta-integralEtaBG)/pow(nEventsBin1,2),2) );

		delete canvasMass;

        return kTRUE;
	}

private:
	TH1D* fGammaGamma;
	TH1D* fBck;
	TH1D* fBckNormPi;
	TH1D* fBckNormEta;

	TH1D* fSignalPi;
	TF1* fFitRecoPi;
	TH1D* fSignalEta;
	TF1* fFitRecoEta;

	Double_t widthPi;
	Double_t widthPiErr;
	Double_t massPi;
	Double_t massPiErr;

	Double_t widthEta;
	Double_t widthEtaErr;
	Double_t massEta;
	Double_t massEtaErr;

	Double_t ratioPi0;
	Double_t ratioPi0Err;
	Double_t ratioEta;
	Double_t ratioEtaErr;

	void ResetPointer(){
		fGammaGamma = 0x0;
		fBck = 0x0;
		fBckNormPi = 0x0;
		fBckNormEta = 0x0;
		fSignalPi = 0x0;
		fFitRecoPi = 0x0;
		fSignalEta = 0x0;
		fFitRecoEta = 0x0;
		return;
	}

	void Init(){
		if(fGammaGamma) delete fGammaGamma;
		if(fBck) delete fBck;
		if(fBckNormPi) delete fBckNormPi;
		if(fBckNormEta) delete fBckNormEta;
		if(fSignalPi) delete fSignalPi;
		if(fFitRecoPi) delete fFitRecoPi;
		if(fSignalEta) delete fSignalEta;
		if(fFitRecoEta) delete fFitRecoEta;

		ResetPointer();

		widthPi = 0;
		widthPiErr = 0;
		massPi = 0;
		massPiErr = 0;

		widthEta = 0;
		widthEtaErr = 0;
		massEta = 0;
		massEtaErr = 0;

		ratioPi0 = 0;
		ratioPi0Err = 0;
		ratioEta = 0;
		ratioEtaErr = 0;

		return;
	}
};

//**************************************************************************************************************************
//*********************** declaration of functions defined in this header ***********************************
Bool_t CheckForData8TeV(TString set);
Bool_t CheckForTriggerData8TeV(TString set);
Bool_t CheckGoodCell(CellQAObj* obj, Int_t cellID);
void SaveCanvas(TCanvas* cvs, TString f, Bool_t logx = kFALSE, Bool_t logy = kFALSE, Bool_t logz = kFALSE);
void SaveWriteCanvas(TCanvas* cvs, TString f, Bool_t logx = kFALSE, Bool_t logy = kFALSE, Bool_t logz = kFALSE);
void EditTH1(std::vector<TString> &globalRuns, Bool_t doEquiDistantXaxis, TH1* hist, Style_t markerStyle, Size_t markerSize, Color_t markerColor, Color_t lineColor);
void OnlyEditTH1(TH1* hist, Style_t markerStyle, Size_t markerSize, Color_t markerColor, Color_t lineColor);
void EditTH1NoRunwise(TH1* hist, Style_t markerStyle, Size_t markerSize, Color_t markerColor, Color_t lineColor, Float_t xOffset=1, Float_t yOffset=0.8);
void EditTH2(TH2 *hist, Float_t titleOffsetX, Float_t titleOffsetY);
void EditRunwiseHists(TH1D *hist, Int_t nHist, TString title, Double_t Ymin = -1, Double_t Ymax = -1);
Bool_t readin(TString fileRuns, std::vector<TString> &vec, Bool_t output = kTRUE, Bool_t badCells = kFALSE);
void DrawAutoGammaCompare3H(TH1* histo1,
					 TH1* histo2,
					 TH1* histo3,
					 TString Title, TString XTitle, TString YTitle,
					 Bool_t YRangeMax, Double_t YMaxFactor, Double_t YMinimum,
					 Bool_t YRange, Double_t YMin , Double_t YMax,
					 Bool_t XRange, Double_t XMin, Double_t XMax,
					 Float_t xOffset=1., Float_t yOffset=1.7, TString data="Data", TString mc1="MC1", TString mc2="MC2");
void DrawAutoGammaHistMatch3H(TH1* histo1,
					 TH1* histo2,
					 TH1* histo3,
					 TString Title, TString XTitle, TString YTitle,
					 Bool_t YRangeMax, Double_t YMaxFactor, Double_t YMinimum,
					 Bool_t YRange, Double_t YMin , Double_t YMax,
					 Bool_t XRange, Double_t XMin, Double_t XMax,
					 Float_t xOffset=1., Float_t yOffset=1.7, TString legS1="Sum", TString legS2="Accepted", TString legS3="Rejected", Bool_t fillStyle=kTRUE);
void CalculateFractionMatches(TH1D* cFraction, TH2F* cESD_Mother, TH2F* cESD_Mother_Matched, Int_t cBin, Double_t ptMin, Double_t ptMax);
void DrawAutoGammaHistoPaper( TH1* histo1,
					TString Title, TString XTitle, TString YTitle,
					Bool_t YRangeMax, Double_t YMaxFactor, Double_t YMinimum,
					Bool_t YRange, Double_t YMin ,Double_t YMax,
                    Bool_t XRange, Double_t XMin, Double_t XMax, Double_t yOffset=1., Double_t yTitleSize = 0.06);
void DrawAutoGammaHistoPaper2D( TH2* histo1,
					TString Title, TString XTitle, TString YTitle,
					Bool_t YRangeMax, Double_t YMaxFactor, Double_t YMinimum,
					Bool_t YRange, Double_t YMin ,Double_t YMax,
					Bool_t XRange, Double_t XMin, Double_t XMax, Double_t xOffset=1, Double_t yOffset=1.);
void PlotCaloQAModule(TH2 *src, Int_t nCaloModules, TString xLabel, TString yLabel, Double_t factorLow = 1, Double_t factorHigh = 1, Int_t fixRange = 0, Double_t range = 0, Bool_t doXRange = kFALSE, Int_t minXRange = 0, Int_t maxXRange = 0, Bool_t kEnergyPlot = kFALSE,
                      Double_t titleOffsetX = 1, Double_t titleOffsetY = 1);
void PlotCellMeanVsSigma(CellQAObj* obj, Int_t nCaloCells, TH2* hist, TString xLabel, TString yLabel, Bool_t XRange, Float_t XMin, Float_t XMax, Bool_t YRange, Float_t YMin, Float_t YMax, Bool_t kEnergy, Bool_t kMC,
                         Double_t titleOffsetX = 1, Double_t titleOffsetY = 1);
TH2D** PlotCellMeanVsSigmaForRunwise(Int_t nCaloCells, TH2* histEnergy, TH2* histTime, TString xLabel, TString yLabel, TString xLabelT, TString yLabelT, Bool_t kMC, Double_t titleOffsetX = 1, Double_t titleOffsetY = 1);
void PlotHotCells(CellQAObj* obj, Int_t iSw, Int_t nCaloCells, TH2* hist, TString xLabel, TString yLabel, Bool_t XRange, Float_t XMin, Float_t XMax, Bool_t YRange, Float_t YMin, Float_t YMax, Bool_t kMC,
                         Double_t titleOffsetX = 1, Double_t titleOffsetY = 1);
void CheckCellsDataMC(CellQAObj* obj, TH2* fHistDataCell, TH2* fHistMCCell, TString xLabel, TString yLabel, Int_t nCaloCells, TString plotData, TString plotMC);
Double_t GetHistogramIntegral(TH1D* hist, Float_t lowX, Float_t highX);
Double_t GetHistogramIntegralError(TH1D* hist, Float_t lowX, Float_t highX);
//**********************************************************************************************************

//*********************** definition of functions **********************************************************
Bool_t CheckForData8TeV(TString set)
{
    if(   set.BeginsWith("LHC12") ||
          set.BeginsWith("LHC15h1") ||
          set.BeginsWith("LHC15h2")
      ) return true;
	else return false;
}

Bool_t CheckForData7TeV(TString set, Bool_t k10=kFALSE)
{
    if(		(set.CompareTo("LHC10_pass4")==0 && k10) ||
            set.CompareTo("LHC10b_pass4")==0 ||
            set.CompareTo("LHC10c_pass4")==0 ||
            set.CompareTo("LHC10d_pass4")==0 ||
            set.CompareTo("LHC10e_pass4")==0 ||
            set.CompareTo("LHC10f_pass4")==0 ) return true;
    else return false;
}

Bool_t CheckForTriggerData8TeV(TString set)
{
	if(set.BeginsWith("LHC12") && set.Contains("-kEMC")) return true;
	else return false;
}

void CheckNEntries(TH2* hist){
  if(hist->GetEntries()<=1.) cout << "WARNING: Histogram: '" << hist->GetName() << "' has <= 1 entries. Check if running over correct file and if settings of previous macros has been correctly chosen!" << endl;
  return;
}
void CheckNEntries(TH1* hist){
  if(hist->GetEntries()<=1.) cout << "WARNING: Histogram: '" << hist->GetName() << "' has <= 1 entries. Check if running over correct file and if settings of previous macros has been correctly chosen!" << endl;
  return;
}

void CheckForNegativeEntries(TH2* hist){
  for(Int_t iX=0; iX<hist->GetNbinsX(); iX++){
    for(Int_t iY=0; iY<hist->GetNbinsY(); iY++){
      if(hist->GetBinContent(iX,iY) <= 0.){
        hist->SetBinContent(iX,iY,0.);
        hist->SetBinError(iX,iY,0.);
      }
    }
  }
  return;
}

void CheckForNegativeEntries(TH1* hist){
  for(Int_t iX=0; iX<hist->GetNbinsX(); iX++){
    if(hist->GetBinContent(iX) <= 0.){
      hist->SetBinContent(iX,0.);
      hist->SetBinError(iX,0.);
    }
  }
  return;
}

Bool_t CheckForCaloHist(TString temp){
  if(    temp.CompareTo("hCaloNClusters")     == 0
      || temp.CompareTo("hCaloNClustersQA")   == 0 ) return kTRUE;
  else return kFALSE;
}

template<class ForwardIt>
void SetLogBinningXTH(ForwardIt* histoRebin){
	TAxis *axisafter = histoRebin->GetXaxis();
	Int_t bins = axisafter->GetNbins();
	Double_t from = axisafter->GetXmin();
	Double_t to = axisafter->GetXmax();
	Double_t *newbins = new Double_t[bins+1];
	newbins[0] = from;
	Double_t factor = TMath::Power(to/from, 1./bins);
	for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
	axisafter->Set(bins, newbins);
	delete [] newbins;
	return;
}

inline void WriteHistogram(TObject* hist){if(hist) hist->Write();return;}

void SaveCanvasAndWriteHistogram(TCanvas* canvas, TObject* hist, TString output){
	if(hist) hist->Write();
	canvas->Update();
	canvas->SaveAs(output.Data());
	canvas->Clear();
	return;
}

void SaveCanvasOnly(TCanvas* canvas, TString output){
	canvas->Update();
	canvas->SaveAs(output.Data());
	canvas->Clear();
	return;
}

void WriteHistogramTH2DVec(TFile* file, std::vector<TH2D*> vec, TString dir){
	file->cd();
	TDirectory *tmpDir = file->mkdir(dir.Data());
	tmpDir->cd();
	for(Int_t j=0; j<(Int_t) vec.size(); j++)
	{
		WriteHistogram(vec.at(j));
        delete vec.at(j);
	}
	return;
}
void WriteHistogramTH1DVec(TFile* file, std::vector<TH1D*> vec, TString dir){
	file->cd();
	TDirectory *tmpDir = file->mkdir(dir.Data());
	tmpDir->cd();
	for(Int_t j=0; j<(Int_t) vec.size(); j++)
	{
		WriteHistogram(vec.at(j));
		delete vec.at(j);
	}
	return;
}
void WriteHistogramTProfileVec(TFile* file, std::vector<TProfile*> vec, TString dir){
	file->cd();
	TDirectory *tmpDir = file->mkdir(dir.Data());
	tmpDir->cd();
	for(Int_t j=0; j<(Int_t) vec.size(); j++)
	{
		WriteHistogram(vec.at(j));
		delete vec.at(j);
	}
	return;
}
void WriteHistogramTF1Vec(TFile* file, std::vector<TF1*> vec, TString dir){
	file->cd();
	TDirectory *tmpDir = file->mkdir(dir.Data());
	tmpDir->cd();
	for(Int_t j=0; j<(Int_t) vec.size(); j++)
	{
		WriteHistogram(vec.at(j));
		delete vec.at(j);
	}
	return;
}

void DeleteVecTH1D(std::vector<TH1D*> vec){
	for(Int_t i=0; i<(Int_t) vec.size(); i++){
        delete vec.at(i);
	}
}
void DeleteVecTH2D(std::vector<TH2D*> vec){
	for(Int_t i=0; i<(Int_t) vec.size(); i++){
		delete vec.at(i);
	}
}
void DeleteVecTString(std::vector<TString> vec){
	for(Int_t i=0; i<(Int_t) vec.size(); i++){
		delete vec.at(i);
	}
}

void SaveCanvas(TCanvas* cvs, TString f, Bool_t logx, Bool_t logy, Bool_t logz)
{
	cvs->SetLogx(logx);
	cvs->SetLogy(logy);
	cvs->SetLogz(logz);
	cvs->Update();
	cvs->SaveAs(f.Data());
	cvs->Clear();
	return;
}

void SaveWriteCanvas(TCanvas* cvs, TString f, Bool_t logx, Bool_t logy, Bool_t logz)
{
	cvs->SetLogx(logx);
	cvs->SetLogy(logy);
	cvs->SetLogz(logz);
	cvs->Update();
	cvs->SaveAs(f.Data());
	cvs->Write();
	cvs->Clear();
	return;
}

void SetXRange(TH2* hist, Int_t min, Int_t max){

    if( max < min ) cout << "ERROR: SetXRangeTH2, max (" << max << ") < min (" << min << ")" << endl;
    if( min < 1 ) min = 1;
    if( max > hist->GetNbinsX() ) max = hist->GetNbinsX();

    hist->GetXaxis()->SetRange(min,max);

    return;
}

void SetYRange(TH2* hist, Int_t min, Int_t max){

    if( max < min ) cout << "ERROR: SetYRangeTH2, max (" << max << ") < min (" << min << ")" << endl;
    if( min < 1 ) min = 1;
    if( max > hist->GetNbinsY() ) max = hist->GetNbinsY();

    hist->GetYaxis()->SetRange(min,max);

    return;
}

void SetXRange(TH1* hist, Int_t min, Int_t max){

    if( max < min ) cout << "ERROR: SetXRangeTH1, max (" << max << ") < min (" << min << ")" << endl;
    if( min < 1 ) min = 1;
    if( max > hist->GetNbinsX() ) max = hist->GetNbinsX();

    hist->GetXaxis()->SetRange(min,max);

    return;
}

void SetZMinMaxTH2(TH2* hist, Int_t minX, Int_t maxX, Int_t minY, Int_t maxY, Bool_t setZero = kFALSE){
	if(hist->GetEntries()<2) return;

    if(setZero){
      for(Int_t iX=minX; iX<=maxX; iX++){
        for(Int_t iY=minY; iY<=maxY; iY++){
          Double_t temp = hist->GetBinContent(iX,iY);
          if(temp<=0.) hist->SetBinContent(iX,iY,0.);
        }
      }
    }

    Double_t min = 0;
    Double_t max = 0;
    Bool_t bSet = kTRUE;

    for(Int_t iX=minX; iX<=maxX; iX++){
        for(Int_t iY=minY; iY<=maxY; iY++){
            Double_t temp = hist->GetBinContent(iX,iY);
            if(temp!=0.){
              if(bSet){
                min = temp;
                max = temp;
                bSet = kFALSE;
              }
              if(temp > max) max = temp;
              if(temp < min) min = temp;
            }
        }
    }

    hist->GetZaxis()->SetRangeUser(min,max);

	return;
}

void SetTH1RangeForLogY(TH1* hist, Int_t minX, Int_t maxX){
    for(Int_t iX=minX; iX<=maxX; iX++){
        Double_t temp = hist->GetBinContent(iX);
        if(temp<=0.) hist->SetBinContent(iX,0.);
    }
    return;
}

void GetMinMaxBin(TH1* hist, Int_t &iMin, Int_t &iMax){
    if(!hist) {cout << "INFO: NULL pointer given to GetMinMaxBin!'" << endl; return;}
	iMin = 1;
	iMax = hist->GetNbinsX();
	Double_t min = hist->GetBinContent(iMin);
	Double_t max = hist->GetBinContent(iMax);

    while(min==0.){
		min = hist->GetBinContent(++iMin);
		if(iMin==iMax){
			iMin=1;
			break;
		}
	}

    while(max==0.){
		max = hist->GetBinContent(--iMax);
		if(iMax==1){
			iMax = hist->GetNbinsX();
            cout << "INFO: Empty histogram given to GetMinMaxBin: '" << hist->GetName() << "'" << endl;
			break;
		}
	}
	return;
}

void GetMinMaxBin(TH2* hist2D, Int_t &iMin, Int_t &iMax){
    if(!hist2D){cout << "INFO: NULL pointer given to GetMinMaxBin!'" << endl; return;}
    //if(hist2D->IsA()!=TH2::Class()){cout << Form("INFO: No TH2 given to GetMinMaxBin(TH2*): %s!'",hist2D->GetName()) << endl; return;}
    TH1* hist = (TH1*) hist2D->ProjectionX(Form("%s_project",hist2D->GetName()),1,hist2D->GetNbinsY());
	GetMinMaxBin(hist, iMin, iMax);
	delete hist;
	return;
}

void GetMinMaxBinY(TH2* hist2D, Int_t &iMin, Int_t &iMax){
    if(!hist2D){cout << "INFO: NULL pointer given to GetMinMaxBin!'" << endl; return;}
    //if(hist2D->IsA()!=TH2::Class()){cout << Form("INFO: No TH2 given to GetMinMaxBin(TH2*): %s!'",hist2D->GetName()) << endl; return;}
    TH1* hist = (TH1*) hist2D->ProjectionY(Form("%s_project",hist2D->GetName()),1,hist2D->GetNbinsX());
	GetMinMaxBin(hist, iMin, iMax);
	delete hist;
	return;
}

void GetMinMaxBin(std::vector<TH1D*> &vec, Int_t &iMin, Int_t &iMax){

    if(vec.size()<1) { cout << "WARNING: in GetMinMaxBin(std::vector<TH1D*>,Int_t,Int_t - vector size <1!" << endl; return;}
	TH1* temp = (TH1*) vec.at(0);
	GetMinMaxBin(temp,iMin,iMax);

    Int_t iMinTemp = 0;
	Int_t iMaxTemp = 0;
	for(Int_t i=1; i<(Int_t) vec.size(); i++){
		temp = (TH1*) vec.at(i);
		if(temp->Integral()==0) continue;
		GetMinMaxBin(temp,iMinTemp,iMaxTemp);
		if(iMinTemp<iMin) iMin = iMinTemp;
		if(iMaxTemp>iMax) iMax = iMaxTemp;
	}

	return;
}

void GetMinMaxBin(std::vector<TH2D*> &vec, Int_t &iMin, Int_t &iMax){

    if(vec.size()<1) {cout << "WARNING: in GetMinMaxBin(std::vector<TH2D*>,Int_t,Int_t - vector size <1!" << endl; return;}
	TH2* tempTH2 = (TH2*) vec.at(0);
	TH1* temp = (TH1*) tempTH2->ProjectionX("proj",1,tempTH2->GetNbinsY());
	GetMinMaxBin(temp,iMin,iMax);
	delete temp;

	Int_t iMinTemp = 0;
	Int_t iMaxTemp = 0;
	for(Int_t i=1; i<(Int_t) vec.size(); i++){
		tempTH2 = (TH2*) vec.at(i);
		temp = (TH1*) tempTH2->ProjectionX("proj",1,tempTH2->GetNbinsY());
		if(temp->Integral()==0) continue;
		GetMinMaxBin(temp,iMinTemp,iMaxTemp);
		if(iMinTemp<iMin) iMin = iMinTemp;
		if(iMaxTemp>iMax) iMax = iMaxTemp;
		delete temp;
	}
	return;
}

Double_t GetMaximumBinValueTH1(TH1* hist){
	Double_t value = hist->GetBinContent(1);
	Double_t binX = hist->GetBinCenter(1);
	for(Int_t i=2; i<hist->GetNbinsX()+1; i++){
		Double_t temp = hist->GetBinContent(i);
		if(temp>value){
			value = temp;
			binX = hist->GetBinCenter(i);
		}
	}
	return binX;
}

void EditTH1(std::vector<TString> &globalRuns, Bool_t doEquiDistantXaxis, TH1* hist, Style_t markerStyle, Size_t markerSize, Color_t markerColor, Color_t lineColor )
{
	hist->SetMarkerStyle(markerStyle);
	hist->SetMarkerSize(markerSize);
	hist->SetMarkerColor(markerColor);
	hist->SetLineColor(lineColor);

	hist->GetXaxis()->SetLabelFont(42);
	hist->GetXaxis()->SetTitleFont(62);
    hist->GetXaxis()->SetLabelSize(0.03);
	hist->GetXaxis()->SetTitleSize(0.035);
    hist->GetXaxis()->SetTitleOffset(0.9);
	hist->GetXaxis()->SetNoExponent();

	hist->GetYaxis()->SetLabelFont(42);
	hist->GetYaxis()->SetTitleFont(62);
    hist->GetYaxis()->SetLabelSize(0.03);
	hist->GetYaxis()->SetTitleSize(0.035);
	hist->GetYaxis()->SetTitleOffset(1.1);
	hist->GetYaxis()->SetDecimals();

	if( doEquiDistantXaxis )
	{
		for(Int_t i=0; i<(Int_t)globalRuns.size(); i++)
		{
			hist->GetXaxis()->SetBinLabel(i+1, globalRuns.at(i).Data());
			hist->GetXaxis()->LabelsOption("v");
		}

		if((Int_t)globalRuns.size()<80) hist->GetXaxis()->SetLabelSize(0.04);
		else hist->GetXaxis()->SetLabelSize(0.02);
		hist->GetXaxis()->SetTitleOffset(-0.3);
	}

	return;
}

void OnlyEditTH1(TH1* hist, Style_t markerStyle, Size_t markerSize, Color_t markerColor, Color_t lineColor )
{
	hist->SetMarkerStyle(markerStyle);
	hist->SetMarkerSize(markerSize);
	hist->SetMarkerColor(markerColor);
	hist->SetLineColor(lineColor);

	hist->GetXaxis()->SetLabelFont(42);
	hist->GetXaxis()->SetTitleFont(62);
	hist->GetXaxis()->SetLabelSize(0.02);
	hist->GetXaxis()->SetTitleSize(0.035);
	hist->GetXaxis()->SetTitleOffset(1.3);
    hist->GetXaxis()->SetNoExponent();

	hist->GetYaxis()->SetLabelFont(42);
	hist->GetYaxis()->SetTitleFont(62);
	hist->GetYaxis()->SetLabelSize(0.02);
	hist->GetYaxis()->SetTitleSize(0.035);
	hist->GetYaxis()->SetTitleOffset(1.1);
	hist->GetYaxis()->SetDecimals();

	return;
}

void EditTH1NoRunwise(TH1* hist, Style_t markerStyle, Size_t markerSize, Color_t markerColor, Color_t lineColor, Float_t xOffset, Float_t yOffset)
{
	hist->SetMarkerStyle(markerStyle);
	hist->SetMarkerSize(markerSize);
	hist->SetMarkerColor(markerColor);
	hist->SetLineColor(lineColor);

	hist->GetXaxis()->SetLabelFont(42);
	hist->GetXaxis()->SetTitleFont(62);
    hist->GetXaxis()->SetTitleSize(0.04);
	hist->GetXaxis()->SetLabelSize(0.035);
	hist->GetXaxis()->SetTitleOffset(xOffset);

	hist->GetYaxis()->SetLabelFont(42);
	hist->GetYaxis()->SetTitleFont(62);
	hist->GetYaxis()->SetLabelSize(0.035);
    hist->GetYaxis()->SetTitleSize(0.04);
	hist->GetYaxis()->SetNoExponent();
	hist->GetYaxis()->SetTitleOffset(yOffset);

	return;
}


void EditTH2(TH2* hist, Float_t titleOffsetX, Float_t titleOffsetY)
{
	hist->GetYaxis()->SetLabelFont(42);
	hist->GetXaxis()->SetLabelFont(42);
	hist->GetYaxis()->SetTitleFont(62);
	hist->GetXaxis()->SetTitleFont(62);

	hist->GetYaxis()->SetTitleSize(0.08);
	hist->GetXaxis()->SetTitleSize(0.08);
	hist->GetYaxis()->SetLabelSize(0.07);
	hist->GetXaxis()->SetLabelSize(0.07);
	hist->GetZaxis()->SetLabelSize(0.07);
	hist->GetYaxis()->SetDecimals();
	hist->GetYaxis()->SetTitleOffset(titleOffsetY);
	hist->GetXaxis()->SetTitleOffset(titleOffsetX);
	hist->SetOption("colz");

	return;
}

void EditTH1ForRunwise(TH1* hist, Float_t titleOffsetX, Float_t titleOffsetY)
{
	hist->GetYaxis()->SetLabelFont(42);
	hist->GetXaxis()->SetLabelFont(42);
	hist->GetYaxis()->SetTitleFont(62);
	hist->GetXaxis()->SetTitleFont(62);

	hist->GetYaxis()->SetTitleSize(0.08);
	hist->GetXaxis()->SetTitleSize(0.08);
	hist->GetYaxis()->SetLabelSize(0.07);
	hist->GetXaxis()->SetLabelSize(0.07);
	hist->GetZaxis()->SetLabelSize(0.07);
	hist->GetYaxis()->SetDecimals();
	hist->GetYaxis()->SetTitleOffset(titleOffsetY);
	hist->GetXaxis()->SetTitleOffset(titleOffsetX);

	return;
}

void EditRunwiseHists(TH1D* hist, Int_t nHist, TString title, Double_t Ymin, Double_t Ymax)
{
	//2,4,5 20-30
	Int_t markerStyles[14]={2,4,5,20,21,22,23,24,25,26,27,28,29,30};
	Int_t nMarker = nHist % 14;

	Int_t nColor = (nHist + 1) % 45;
	if(nColor == 10){nColor = 45;}
	if(nColor == 18){nColor = 46;}
	if(nColor == 19){nColor = 47;}
	if(nColor == 0){nColor = 48;}

	//cout << "nHist: " << nHist << "nMarker: " << nMarker << "nColor: " << nColor << endl;
	hist->SetMarkerStyle(markerStyles[nMarker]);
	hist->SetMarkerColor(nColor);
	hist->SetLineColor(nColor);
	hist->SetLineWidth(0.8);

	if(Ymin!=-1 && Ymax!=-1) hist->GetYaxis()->SetRangeUser(Ymin,Ymax);
	hist->SetTitle(title.Data());

	hist->GetXaxis()->SetLabelFont(42);
	hist->GetXaxis()->SetTitleFont(62);
	hist->GetXaxis()->SetLabelSize(0.03);
	hist->GetXaxis()->SetTitleSize(0.035);
	hist->GetXaxis()->SetTitleOffset(1.1);
	hist->GetXaxis()->SetNoExponent();

	hist->GetYaxis()->SetLabelFont(42);
	hist->GetYaxis()->SetTitleFont(62);
	hist->GetYaxis()->SetLabelSize(0.03);
	hist->GetYaxis()->SetTitleSize(0.035);
	hist->GetYaxis()->SetTitleOffset(1.1);
	hist->GetYaxis()->SetDecimals();

	return;
}

Bool_t readin(TString fileRuns, std::vector<TString> &vec, Bool_t output, Bool_t badCells)
{
	if(output) cout << Form("\nReading from %s...", fileRuns.Data()) << endl;
	fstream file;
	TString fVar;
	Int_t totalN=0;
	file.open(fileRuns.Data(), ios::in);
	if(file.good())
	{
		file.seekg(0L, ios::beg);
		if(output && !badCells) cout << "Processing Runs: \"";
		if(output && badCells) cout << "Read in bad cells: \"";
		while(!file.eof())
		{
			file >> fVar;
			if(fVar.Sizeof()>1)
			{
				if(output) cout << fVar.Data() << ", ";
				vec.push_back(fVar);
				totalN++;
			}
		}
		if(output) cout << "\"" << endl;
	}
	file.close();
	if(output && !badCells) cout << "...done!\n\nIn total " << totalN << " Runs will be processed!\n" << endl;
	if(output && badCells) cout << "...done!\n\nIn total " << totalN << " bad cells were read in!\n" << endl;
	if(!badCells){
		if(totalN > 0) return kTRUE;
		else return kFALSE;
	}
	else return kTRUE;
}

void DrawFit(TH1D* tempHist, Int_t i, Double_t fittedValue, Float_t* runRanges, TString DataSets, TString plotDataSets, Double_t xD, Double_t yD, Double_t sD, Double_t width){

  TString labelling="";
  if(((TString)tempHist->GetName()).CompareTo(Form("hPi0Mass_%s",DataSets.Data()))==0) labelling = Form("Fit to %s: m_{#pi^{0}} = %.02f MeV",plotDataSets.Data(),fittedValue*1E3);
  else if(((TString)tempHist->GetName()).CompareTo(Form("hEtaMass_%s",DataSets.Data()))==0) labelling = Form("Fit to %s: m_{#eta} = %.02f MeV",plotDataSets.Data(),fittedValue*1E3);
  else if(((TString)tempHist->GetName()).CompareTo(Form("hPi0Width_%s",DataSets.Data()))==0) labelling = Form("Fit to %s: #sigma_{#pi^{0}} = %.02f MeV",plotDataSets.Data(),fittedValue*1E3);
  else if(((TString)tempHist->GetName()).CompareTo(Form("hEtaWidth_%s",DataSets.Data()))==0) labelling = Form("Fit to %s: #sigma_{#eta} = %.02f MeV",plotDataSets.Data(),fittedValue*1E3);

  if(labelling.IsNull()) return;
  else{
    if(runRanges[0]==0 && runRanges[1]==0){
      runRanges[0]=tempHist->GetXaxis()->GetBinLowEdge(1);
      runRanges[1]=tempHist->GetXaxis()->GetBinUpEdge(tempHist->GetNbinsX());
    }
    TLatex* drawMesonMass = new TLatex(xD,yD-(sD*1.25*(i+1)),labelling.Data());
    drawMesonMass->SetTextSize(sD);
    drawMesonMass->SetNDC();
    drawMesonMass->Draw("SAME");
    TF1* drawMesonMassFit = new TF1(Form("drawMesonFit_%s",tempHist->GetName()),Form("%.05f",fittedValue),runRanges[0],runRanges[1]);
    drawMesonMassFit->SetLineColor(tempHist->GetLineColor());
    drawMesonMassFit->SetLineWidth(width);
    drawMesonMassFit->Draw("SAME");
  }

  return;
}

void DrawPeriodQAHistoTH1(TCanvas* canvas,Double_t leftMargin,Double_t rightMargin,Double_t topMargin,Double_t bottomMargin,Bool_t logX, Bool_t logY, Bool_t logZ,
					TH1* fHist, TString title,TString xTitle,TString yTitle, Double_t xOffset, Double_t yOffset)
{
	canvas->cd();
	canvas->SetLeftMargin(leftMargin);
	canvas->SetRightMargin(rightMargin);
	canvas->SetTopMargin(topMargin);
	canvas->SetBottomMargin(bottomMargin);
	//x-axis
	fHist->GetXaxis()->SetTitleOffset(xOffset);
	fHist->SetXTitle(xTitle.Data());
	fHist->GetXaxis()->SetLabelFont(42);
	fHist->GetXaxis()->SetTitleFont(62);
    fHist->GetXaxis()->SetTitleSize(0.04);
    fHist->GetXaxis()->SetLabelSize(0.035);
	//y-axis
	fHist->GetYaxis()->SetTitleOffset(yOffset);
	fHist->SetYTitle(yTitle.Data());
	fHist->GetYaxis()->SetLabelFont(42);
	fHist->GetYaxis()->SetTitleFont(62);
    fHist->GetYaxis()->SetTitleSize(0.04);
    fHist->GetYaxis()->SetLabelSize(0.035);
	fHist->GetYaxis()->SetDecimals();
	//draw
	fHist->SetTitle(title.Data());
	fHist->DrawCopy("e,hist");

	canvas->SetLogx(logX); canvas->SetLogy(logY); canvas->SetLogz(logZ);
	return;
}

void DrawPeriodQAHistoTH1(TCanvas* canvas,Double_t leftMargin,Double_t rightMargin,Double_t topMargin,Double_t bottomMargin,Bool_t logX, Bool_t logY, Bool_t logZ,
					TH1* fHist, TString title,TString xTitle,TString yTitle, Double_t xOffset, Double_t yOffset,
					Double_t p1,Double_t p2,Double_t p3, TString fCollisionSystem, TString plotDataSets, TString fClusters, Int_t textAlign = 11)
{
	DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,logX,logY,logZ,
						fHist,title,xTitle,yTitle,xOffset,yOffset);
	PutProcessLabelAndEnergyOnPlot(p1,p2,p3, fCollisionSystem.Data(), plotDataSets.Data(), fClusters.Data(), 62, 0.03, "", 1, 1.25, textAlign);

	return;
}

void DrawPeriodQAHistoTH2(TCanvas* canvas,Double_t leftMargin,Double_t rightMargin,Double_t topMargin,Double_t bottomMargin,Bool_t logX, Bool_t logY, Bool_t logZ,
					TH2* fHist, TString title,TString xTitle,TString yTitle, Double_t xOffset, Double_t yOffset)
{
	canvas->cd();
	canvas->SetLeftMargin(leftMargin);
	canvas->SetRightMargin(rightMargin);
	canvas->SetTopMargin(topMargin);
	canvas->SetBottomMargin(bottomMargin);
	//x-axis
	fHist->GetXaxis()->SetTitleOffset(xOffset);
	fHist->SetXTitle(xTitle.Data());
	fHist->GetXaxis()->SetLabelFont(42);
	fHist->GetXaxis()->SetTitleFont(62);
    fHist->GetXaxis()->SetTitleSize(0.04);
	fHist->GetXaxis()->SetLabelSize(0.035);
	//y-axis
	fHist->GetYaxis()->SetTitleOffset(yOffset);
	fHist->SetYTitle(yTitle.Data());
	fHist->GetYaxis()->SetLabelFont(42);
	fHist->GetYaxis()->SetTitleFont(62);
    fHist->GetYaxis()->SetTitleSize(0.04);
	fHist->GetYaxis()->SetLabelSize(0.035);
	fHist->GetYaxis()->SetDecimals();
    //z-axis
    fHist->GetZaxis()->SetLabelFont(42);
    fHist->GetZaxis()->SetTitleFont(62);
    fHist->GetZaxis()->SetTitleSize(0.04);
    fHist->GetZaxis()->SetLabelSize(0.035);
	//draw
    if(title.IsNull()) fHist->SetTitle("");
    else{
      if(topMargin<0.06) canvas->SetTopMargin(0.06);
      fHist->SetTitle(title.Data());
    }
	fHist->DrawCopy("colz");

	canvas->SetLogx(logX); canvas->SetLogy(logY); canvas->SetLogz(logZ);
	return;
}

void DrawPeriodQAHistoTH2(TCanvas* canvas,Double_t leftMargin,Double_t rightMargin,Double_t topMargin,Double_t bottomMargin,Bool_t logX, Bool_t logY, Bool_t logZ,
					TH2* fHist, TString title,TString xTitle,TString yTitle, Double_t xOffset, Double_t yOffset,
					Double_t p1,Double_t p2,Double_t p3, TString fCollisionSystem, TString plotDataSets, TString fClusters, Int_t textAlign = 11)
{
	DrawPeriodQAHistoTH2(canvas,leftMargin,rightMargin,topMargin,bottomMargin,logX,logY,logZ,
						fHist,title,xTitle,yTitle,xOffset,yOffset);
	PutProcessLabelAndEnergyOnPlot(p1,p2,p3, fCollisionSystem.Data(), plotDataSets.Data(), fClusters.Data(), 62, 0.03, "", 1, 1.25, textAlign);
	return;
}

void DrawPeriodQACompareHistoTH1(TCanvas* canvas,Double_t leftMargin,Double_t rightMargin,Double_t topMargin,Double_t bottomMargin,Bool_t logX, Bool_t logY, Bool_t logZ,
					std::vector<TH1D*> vecHist, TString title,TString xTitle,TString yTitle, Double_t xOffset, Double_t yOffset,
					TString labelData, Color_t *color, Bool_t adjustRange, Double_t rangeMin, Double_t rangeMax, Bool_t useErrors,
					Double_t p1,Double_t p2,Double_t p3, TString fCollisionSystem, TString* plotDataSets, TString fClusters, Int_t textAlign = 11)
{
    if((Int_t) vecHist.size()<1) {
        cout << "ERROR: DrawPeriodQACompareHistoTH1, vector size <1! Returning...";
		return;
	}

	canvas->cd();
	canvas->SetLeftMargin(leftMargin);
	canvas->SetRightMargin(rightMargin);
	canvas->SetTopMargin(topMargin);
	canvas->SetBottomMargin(bottomMargin);

	//TLegend* leg1 = new TLegend( 0.14,0.93-(0.04*((Int_t) vecHist.size())),0.41,0.95);
//	TLegend* leg1 = new TLegend( 0.77,0.83-(0.04*((Int_t) vecHist.size())),0.97,0.85);
//	leg1->SetTextSize(0.03);
//	leg1->SetFillColor(0);

	TLegend *leg1 = new TLegend(0.16,0.96,0.95,0.99);
	leg1->SetNColumns((Int_t) vecHist.size());
	leg1->SetFillColor(0);
	leg1->SetLineColor(0);
	leg1->SetTextSize(0.03);
	leg1->SetTextFont(42);

	if(adjustRange) AdjustHistRange(vecHist,rangeMin,rangeMax,useErrors,0,0,logY);

//	Double_t leg1X1 = 0.77;
	Int_t iEnd = (Int_t) vecHist.size();
	for(Int_t i=0; i<iEnd; i++){
		TH1D* fHist = (TH1D*) vecHist.at(i);
//		Double_t tempD;
		if(i==0){
            DrawGammaSetMarker(fHist, 24, 0.5, color[i], color[i]);
			leg1->AddEntry(fHist,labelData.Data(),"p");
//			tempD = 1.08-(((Double_t) labelData.Length() * 2.)/100.);
		}else{
            DrawGammaSetMarker(fHist, 1, 0.5, color[i], color[i]);
			leg1->AddEntry(fHist,plotDataSets[i].Data());
//			tempD = 1.08-(((Double_t) plotDataSets[i].Length() * 2.)/100.);
		}
//		if(tempD<leg1X1) leg1X1 = tempD;
	}

	Int_t iStart = (Int_t) vecHist.size()-1;
	for(Int_t i=iStart; i>=0; i--){
		TH1D* fHist = (TH1D*) vecHist.at(i);
		//x-axis
		fHist->GetXaxis()->SetTitleOffset(xOffset);
		fHist->SetXTitle(xTitle.Data());
		fHist->GetXaxis()->SetLabelFont(42);
		fHist->GetXaxis()->SetTitleFont(62);
        fHist->GetXaxis()->SetTitleSize(0.04);
		fHist->GetXaxis()->SetLabelSize(0.035);
		//y-axis
		fHist->GetYaxis()->SetTitleOffset(yOffset);
		fHist->SetYTitle(yTitle.Data());
		fHist->GetYaxis()->SetLabelFont(42);
		fHist->GetYaxis()->SetTitleFont(62);
        fHist->GetYaxis()->SetTitleSize(0.04);
		fHist->GetYaxis()->SetLabelSize(0.035);
		fHist->GetYaxis()->SetDecimals();
		//draw
		fHist->SetTitle(title.Data());

        if(iStart==0) fHist->DrawCopy("x0,e,p");
        else if(i==0) fHist->DrawCopy("x0,e,p,same");
        else if(i==iStart) fHist->DrawCopy("e,hist");
		else fHist->DrawCopy("e,hist,same");
	}
//	leg1->SetX1(leg1X1);
	leg1->Draw();

	PutProcessLabelAndEnergyOnPlot(p1,p2,p3, fCollisionSystem.Data(), fClusters.Data(), "", 62, 0.03, "", 1, 1.25, textAlign);
	canvas->SetLogx(logX); canvas->SetLogy(logY); canvas->SetLogz(logZ);

	return;
}

void DrawPeriodQACompareHistoRatioTH1(TCanvas* canvas,Double_t leftMargin,Double_t rightMargin,Double_t topMargin,Double_t bottomMargin,Bool_t logX, Bool_t logY, Bool_t logZ,
					std::vector<TH1D*> vecHist, TString title,TString xTitle,TString yTitle, Double_t xOffset, Double_t yOffset,
					TString labelData, Color_t *color, Bool_t adjustRange, Double_t rangeMin, Double_t rangeMax, Bool_t useErrors,
					Double_t p1,Double_t p2,Double_t p3, TString fCollisionSystem, TString* plotDataSets, TString fClusters, Int_t textAlign = 11)
{
	if((Int_t) vecHist.size()<2) {
        cout << "INFO: DrawPeriodQACompareHistoRatioTH1, vector size <2! Histogram name at(0): " << endl;
        if((Int_t) vecHist.size()==1) cout << "\t" << vecHist.at(0)->GetName() << "...skipping CompareHistoRatio!";
        else cout << "\nERROR: size of vector is zero!\n";
        cout << " Returning..." << endl;
		return;
	}

	canvas->cd();
	canvas->SetLeftMargin(0.05);
	canvas->SetRightMargin(rightMargin);
	canvas->SetTopMargin(topMargin);
	canvas->SetBottomMargin(bottomMargin);

//	TLegend* leg1 = new TLegend( 0.08,0.93-(0.04*((Int_t) vecHist.size())),0.25,0.95);
//	leg1->SetTextSize(0.03);
//	leg1->SetFillColor(0);

	TLegend *leg1 = new TLegend(0.15,0.96,0.95,0.99);
	leg1->SetNColumns((Int_t) vecHist.size()-1);
	leg1->SetFillColor(0);
	leg1->SetLineColor(0);
	leg1->SetTextSize(0.03);
	leg1->SetTextFont(42);

	std::vector<TH1D*> vecTemp;
    Int_t tempMin = vecHist.at(0)->GetXaxis()->GetFirst();
    Int_t tempMax = vecHist.at(0)->GetXaxis()->GetLast();
    for(Int_t iMC=1; iMC<(Int_t) vecHist.size(); iMC++){
      tempMin = (tempMin > vecHist.at(iMC)->GetXaxis()->GetFirst()) ? tempMin : vecHist.at(iMC)->GetXaxis()->GetFirst();
      tempMax = (tempMax < vecHist.at(iMC)->GetXaxis()->GetLast()) ? tempMax : vecHist.at(iMC)->GetXaxis()->GetLast();
    }
    SetXRange(vecHist.at(0), tempMin, tempMax);

	for(Int_t iMC=1; iMC<(Int_t) vecHist.size(); iMC++){
		vecTemp.push_back((TH1D*) vecHist.at(iMC));
		Int_t iVec = iMC-1;
        if(vecHist.at(0)->GetNbinsX() != vecTemp.at(iVec)->GetNbinsX()){
          cout << "*******" << endl;
          cout << "WARNING: in DrawPeriodQACompareHistoRatioTH1, histogram: " << vecHist.at(0)->GetName() << endl;
          cout << "Different number of bins: " << vecHist.at(0)->GetNbinsX() << " & iVec: " << iVec << ": " << vecTemp.at(iVec)->GetNbinsX() << endl;
          cout << "Correcting by cloning from data hist...";
          TH1D* temp = (TH1D*) vecTemp.at(iVec)->Clone("duplicateHist");
          TH1D* temp2 = (TH1D*) vecHist.at(0)->Clone("newHist");
          for(Int_t iX=temp->GetXaxis()->FindBin(temp->GetXaxis()->GetXmin()); iX<=temp->GetXaxis()->FindBin(temp->GetXaxis()->GetXmax()) ; iX++)
             temp2->SetBinContent(iX,temp->GetBinContent(iX));
          vecTemp.at(iVec) = temp2;
          delete temp;
          cout << "done!" << endl;
          cout << "*******" << endl;
        }

        SetXRange(vecTemp.at(iVec), tempMin, tempMax);
//                cout << vecHist.at(0)->GetName() << endl;
//                cout << vecHist.at(0)->GetXaxis()->GetXmax() << endl;
//                cout << vecTemp.at(iVec)->GetXaxis()->GetXmax() << endl;cout << "*******" << endl;
//                cout << vecHist.at(0)->GetXaxis()->GetFirst() << endl;
//                cout << vecTemp.at(iVec)->GetXaxis()->GetFirst() << endl;cout << "*******" << endl;
//                cout << vecHist.at(0)->GetXaxis()->GetLast() << endl;
//                cout << vecTemp.at(iVec)->GetXaxis()->GetLast() << endl;cout << "*******" << endl;
		vecTemp.at(iVec)->Divide(vecHist.at(0));
		//x-axis
		vecTemp.at(iVec)->GetXaxis()->SetTitleOffset(xOffset);
		vecTemp.at(iVec)->SetXTitle(xTitle.Data());
		vecTemp.at(iVec)->GetXaxis()->SetLabelFont(42);
		vecTemp.at(iVec)->GetXaxis()->SetTitleFont(62);
        vecTemp.at(iVec)->GetXaxis()->SetTitleSize(0.04);
		vecTemp.at(iVec)->GetXaxis()->SetLabelSize(0.035);
		//y-axis
		vecTemp.at(iVec)->GetYaxis()->SetTitleOffset(yOffset);
        vecTemp.at(iVec)->SetYTitle("");
        vecTemp.at(iVec)->GetYaxis()->SetRangeUser(0,2);
		vecTemp.at(iVec)->GetYaxis()->SetLabelFont(42);
		vecTemp.at(iVec)->GetYaxis()->SetTitleFont(62);
        vecTemp.at(iVec)->GetYaxis()->SetTitleSize(0.04);
		vecTemp.at(iVec)->GetYaxis()->SetLabelSize(0.035);
		vecTemp.at(iVec)->GetYaxis()->SetDecimals();
		//draw
		vecTemp.at(iVec)->SetTitle(title.Data());
	}

	//if(adjustRange) AdjustHistRange(vecTemp,rangeMin,rangeMax,useErrors);

//	Double_t leg1X2 = 0.25;
	for(Int_t iMC=0; iMC<(Int_t) vecTemp.size(); iMC++){
        DrawGammaSetMarker(vecTemp.at(iMC), 24, 0.5, color[iMC+1], color[iMC+1]);
		TString tempStr = Form("%s / %s",plotDataSets[iMC+1].Data(), labelData.Data());
//		Double_t tempD = (((Double_t) tempStr.Length())/100.)+0.15;
//		if(tempD>leg1X2) leg1X2 = tempD;
		leg1->AddEntry(vecTemp.at(iMC), tempStr.Data(),"p");

		if(iMC==0) vecTemp.at(iMC)->DrawCopy("x0,e,p");
		else vecTemp.at(iMC)->DrawCopy("x0,e,p,same");
	}
	TLine *line = new TLine();
	line->SetLineStyle(2);
	line->SetLineWidth(2);
	line->SetLineColor(1);
	line->DrawLine(vecHist.at(0)->GetXaxis()->GetBinLowEdge(vecHist.at(0)->GetXaxis()->GetFirst()),1,
				   vecHist.at(0)->GetXaxis()->GetBinUpEdge(vecHist.at(0)->GetXaxis()->GetLast()),1);
	//leg1->SetX2(leg1X2);
	leg1->Draw();
	vecTemp.clear();

	PutProcessLabelAndEnergyOnPlot(p1,p2,p3, fCollisionSystem.Data(), fClusters.Data(), "", 62, 0.03, "", 1, 1.25, textAlign);
	canvas->SetLogx(kFALSE); canvas->SetLogy(kFALSE); canvas->SetLogz(kFALSE);

	return;
}

void DrawVectorOverviewTH2D(TCanvas* canvas, std::vector<TH2D*> vec, TString saveString, TString outputDirDataSet, TString suffix,
							Double_t padL, Double_t padR, Double_t padT, Double_t padB, Double_t offsetX, Double_t offsetY,
							Double_t latX, Double_t latY, TBox* box, Bool_t logY, Bool_t logZ){
	Int_t iD = 1;
	Int_t varDraw = 0;
	canvas->Divide(4,4);
	for(Int_t j=0; j<(Int_t) vec.size(); j++, iD++)
	{
		if( iD > 16 )
		{
			SaveCanvas(canvas, Form("%s/%s_p%i.%s", outputDirDataSet.Data(), saveString.Data(), varDraw++, suffix.Data()));
			canvas->Divide(4,4);
			iD = 1;
		}
		TPad* currPad = (TPad*) canvas->cd(iD);
		currPad->SetLeftMargin(padL);
		currPad->SetRightMargin(padR);
		currPad->SetTopMargin(padT);
		currPad->SetBottomMargin(padB);
		currPad->SetLogy(logY);
		currPad->SetLogz(logZ);

		TString title = ((TH2D*) vec.at(j))->GetTitle();
		EditTH2(((TH2D*) vec.at(j)), offsetX, offsetY);
        ((TH2D*) vec.at(j))->SetTitle("");
		((TH2D*) vec.at(j))->Draw();

		if(box) box->Draw();
		TLatex *alice = new TLatex(latX,latY,title.Data());
		alice->SetNDC();
		alice->SetTextColor(1);
		alice->SetTextSize(0.08);
		alice->Draw();
	}
	SaveCanvas(canvas, Form("%s/%s_p%i.%s", outputDirDataSet.Data(), saveString.Data(), varDraw, suffix.Data()));

	return;
}

void DrawVectorOverviewTH1D(TCanvas* canvas, std::vector<TH1D*> vec, TString saveString, TString outputDirDataSet, TString suffix,
							Double_t padL, Double_t padR, Double_t padT, Double_t padB, Double_t offsetX, Double_t offsetY,
							Double_t latX, Double_t latY, TBox* box, Bool_t logY, Bool_t logZ){
	Int_t iD = 1;
	Int_t varDraw = 0;
	canvas->Divide(4,4);
	for(Int_t j=0; j<(Int_t) vec.size(); j++, iD++)
	{
		if( iD > 16 )
		{
			SaveCanvas(canvas, Form("%s/%s_p%i.%s", outputDirDataSet.Data(), saveString.Data(), varDraw++, suffix.Data()));
			canvas->Divide(4,4);
			iD = 1;
		}
		TPad* currPad = (TPad*) canvas->cd(iD);
		currPad->SetLeftMargin(padL);
		currPad->SetRightMargin(padR);
		currPad->SetTopMargin(padT);
		currPad->SetBottomMargin(padB);
		currPad->SetLogy(logY);
		currPad->SetLogz(logZ);

		TString title = ((TH1D*) vec.at(j))->GetTitle();
		EditTH1ForRunwise(((TH1D*) vec.at(j)), offsetX, offsetY);
        ((TH1D*) vec.at(j))->SetTitle("");
		((TH1D*) vec.at(j))->Draw();

		if(box) box->Draw();
		TLatex *alice = new TLatex(latX,latY,title.Data());
		alice->SetNDC();
		alice->SetTextColor(1);
		alice->SetTextSize(0.08);
		alice->Draw();
	}
	SaveCanvas(canvas, Form("%s/%s_p%i.%s", outputDirDataSet.Data(), saveString.Data(), varDraw, suffix.Data()));

	return;
}

void DrawVectorOverviewMissingCells(TCanvas* canvas, std::vector<TH2D*> vec, std::vector<TH2D*> vecClusterFired, std::vector<TH2D*> &vecMissingIDs,
									TString saveString, TString outputDirDataSet, TString suffix, TString DataSets){
	Int_t iD = 1;
	Int_t varDraw = 0;
	canvas->Divide(4,4);
	for(Int_t j=0; j<(Int_t) vec.size(); j++, iD++)
	{
		if( iD > 16 ){
			SaveCanvas(canvas, Form("%s/%s_p%i.%s", outputDirDataSet.Data(), saveString.Data(),varDraw++,suffix.Data()));
			canvas->Divide(4,4);
			iD = 1;
		}
		TPad* currPad = (TPad*) canvas->cd(iD);
		currPad->SetLeftMargin(0.13);
		currPad->SetRightMargin(0.15);
		currPad->SetTopMargin(0.1);
		currPad->SetBottomMargin(0.14);

		TString title = ((TH2D*) vec.at(j))->GetTitle();
		EditTH1ForRunwise(((TH2D*) vec.at(j)), 0.8, 0.9);
        ((TH2D*) vec.at(j))->SetTitle("");

		TH2D* tempClusterMissing = (TH2D*) ((TH2D*) vec.at(j))->Clone(Form("%s-%s",title.Data(),DataSets.Data()));
		tempClusterMissing->Add(((TH2D*) vecClusterFired.at(j)),-1);
		tempClusterMissing->GetZaxis()->SetRangeUser(0,2);
		tempClusterMissing->DrawCopy("COLZ");
		tempClusterMissing->SetTitle(title.Data());
		TString tempTitle = title;
		tempTitle.Remove(6,19);
		tempClusterMissing->SetName(Form("%s-%s",tempTitle.Data(),DataSets.Data()));
		vecMissingIDs.push_back(tempClusterMissing);
		((TH2D*) vec.at(j))->SetTitle(title.Data());

		TLatex *alice = new TLatex(0.25,0.93,title.Data());
		alice->SetNDC();
		alice->SetTextColor(1);
		alice->SetTextSize(0.08);
		alice->Draw();
	}
	SaveCanvas(canvas, Form("%s/%s_p%i.%s", outputDirDataSet.Data(), saveString.Data(),varDraw, suffix.Data()));

	return;
}

void DrawVectorRunwiseTH1D(TCanvas *canvas, TLegend *legendRuns, std::vector<TH1D*> vec, std::vector<TString> vecRuns,
						   Double_t lowerAdjust, Double_t higherAdjust, Bool_t adjustIncludeError, Double_t addRight,
						   Double_t putLabel1, Double_t putLabel2, Double_t putLabel3, Double_t putLabel4, Double_t putLabel5, Double_t putLabel6,
						   Bool_t doTrigger, TString fTrigger, Bool_t data, TString outputDirDataSet, TString saveString, TString plotDataSets, Bool_t doXaxisExp,
						   TString fCollisionSystem, TString calo, TString suffix, Bool_t logX, Bool_t logY, Bool_t logZ, Bool_t doXRange = kTRUE, Bool_t doYRange = kTRUE){

    if(vec.size()==0) { cout << "WARNING: in DrawVectorRunwiseTH1D - vector size 0!" << endl; return;}
    if(vec.size()!=vecRuns.size()) { cout << "WARNING: in DrawVectorRunwiseTH1D - vec.size()!=vecRuns.size()" << endl; return;}
	Int_t min = 0;
	Int_t max = 0;
	if(doXRange) GetMinMaxBin(vec,min,max);

	if(doYRange) AdjustHistRange(vec,lowerAdjust,higherAdjust,adjustIncludeError);
	for(Int_t h=0; h<(Int_t) vec.size(); h++)
	{
		TString draw = (h==0)?"p":"p, same";
		EditRunwiseHists(((TH1D*) vec.at(h)), h, "");
		if(doXaxisExp) ((TH1D*) vec.at(h))->GetXaxis()->SetNoExponent(kFALSE);
        if(doXRange) SetXRange(((TH1D*) vec.at(h)), min, max);
		((TH1D*) vec.at(h))->Draw(draw.Data());
		legendRuns->AddEntry(((TH1D*) vec.at(h)),Form("%s",vecRuns.at(h).Data()),"p");
	}
	legendRuns->Draw();

    //if(calo.IsNull()) putLabel5+=0.04;
    if(doTrigger && data){
      PutProcessLabelAndEnergyOnPlot(putLabel1-addRight, putLabel2, putLabel3, fCollisionSystem.Data(), plotDataSets.Data(), fTrigger.Data());
      PutProcessLabelAndEnergyOnPlot(putLabel4-addRight, putLabel5, putLabel6, calo.Data(), "", "");
    }else{
      PutProcessLabelAndEnergyOnPlot(putLabel1-addRight, putLabel2, putLabel3, fCollisionSystem.Data(), plotDataSets.Data(), calo.Data());
    }
	SaveCanvas(canvas, Form("%s/%s_%s.%s", outputDirDataSet.Data(), saveString.Data(), plotDataSets.Data(),suffix.Data()),logX,logY,logZ);
	legendRuns->Clear();

	return;
}

void DrawVectorRunwiseBadCells(TCanvas *canvas, TLegend *legendRuns, std::vector<TH1D*> vec, Int_t badCellsSize, Int_t iBad, std::vector<TString> vecRuns,
						   Double_t lowerAdjust, Double_t higherAdjust, Bool_t adjustIncludeError, Double_t addRight,
						   Double_t putLabel1, Double_t putLabel2, Double_t putLabel3, Double_t putLabel4, Double_t putLabel5, Double_t putLabel6,
						   Bool_t doTrigger, TString fTrigger, Bool_t data, TString outputDirDataSet, TString saveString, TString plotDataSets, Bool_t doXaxisExp,
						   TString fCollisionSystem, TString calo, TString suffix, Bool_t logX, Bool_t logY, Bool_t logZ){

	Int_t hC = 0;
	AdjustHistRange(vec,lowerAdjust,higherAdjust,adjustIncludeError);
	for(Int_t h=0; h<(Int_t) vec.size(); h++){
		if(h%badCellsSize!=iBad) continue;
		TString draw = (hC==0)?"p":"p, same";
		EditRunwiseHists(((TH1D*) vec.at(h)), hC, "");
		if(doXaxisExp) ((TH1D*) vec.at(h))->GetXaxis()->SetNoExponent(kFALSE);
		((TH1D*) vec.at(h))->Draw(draw.Data());
		legendRuns->AddEntry(((TH1D*) vec.at(h)),Form("%s",vecRuns.at(hC).Data()),"p");
		hC++;
	}
	legendRuns->Draw();

	PutProcessLabelAndEnergyOnPlot(putLabel1-addRight, putLabel2, putLabel3, fCollisionSystem.Data(), calo.Data(), plotDataSets.Data());
	if( doTrigger && data) PutProcessLabelAndEnergyOnPlot(putLabel4-addRight, putLabel5, putLabel6, fTrigger.Data(), "", "");
	SaveCanvas(canvas, Form("%s/%s_%s.%s",outputDirDataSet.Data(),saveString.Data(),plotDataSets.Data(),suffix.Data()),logX,logY,logZ);
	legendRuns->Clear();

	return;
}

void DrawAutoGammaCompare3H( TH1* histo1,
                     TH1* histo2,
                     TH1* histo3,
                     TString Title, TString XTitle, TString YTitle,
                     Bool_t YRangeMax, Double_t YMaxFactor, Double_t YMinimum,
                     Bool_t YRange, Double_t YMin ,Double_t YMax,
                     Bool_t XRange, Double_t XMin, Double_t XMax,
                     Float_t xOffset, Float_t yOffset, TString data, TString mc1, TString mc2) {
    if (YRangeMax && !XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        if(maxRangeR < histo2->GetMaximum()){
            maxRangeR = histo2->GetMaximum();
        }
        Double_t minRangeR = histo1->GetMinimum();
        if(minRangeR > histo2->GetMinimum()){
            minRangeR = histo2->GetMinimum();
        }
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);
    }
    if (YRangeMax && XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        if(maxRangeR < histo2->GetMaximum()){
            maxRangeR = histo2->GetMaximum();
        }
        Double_t minRangeR = histo1->GetMinimum();
        if(minRangeR > histo2->GetMinimum()){
            minRangeR = histo2->GetMinimum();
        }
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);
    }
    if (YRange && XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);
    }
    if (!YRangeMax && !YRange && XRange){
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);
    }
    if (YRange && !XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);
    }

    if(Title.CompareTo("") != 0){
        histo1->SetTitle(Title.Data());
    }else{histo1->SetTitle("");
    histo2->SetTitle("");}
    if(XTitle.CompareTo("") != 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.CompareTo("") != 0){
        histo1->SetYTitle(YTitle.Data());
    }

    DrawGammaSetMarker(histo1, 1, 0.5, kRed+1, kRed+1);
    DrawGammaSetMarker(histo2, 1, 0.5, kMagenta+2, kMagenta+2);
    DrawGammaSetMarker(histo3, 24, 0.5, kBlack, kBlack);
    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);

    histo1->GetYaxis()->SetLabelSize(0.035);
    histo1->GetYaxis()->SetTitleSize(0.04);
    histo1->GetYaxis()->SetDecimals();
    histo1->GetYaxis()->SetTitleOffset(yOffset);
    histo1->GetXaxis()->SetTitleSize(0.04);
    histo1->GetXaxis()->SetLabelSize(0.035);
    histo1->GetXaxis()->SetTitleOffset(xOffset);
    histo1->DrawCopy("e,hist");
    histo2->DrawCopy("e,hist,same");
    histo3->DrawCopy("x0,e,p,same");

    TLegend* leg1 = new TLegend( 0.77,0.81,0.97,0.95);
    leg1->SetTextSize(0.04);
    leg1->SetFillColor(0);
    leg1->AddEntry(histo3,data.Data(),"p");
    leg1->AddEntry(histo1,mc1.Data());
    leg1->AddEntry(histo2,mc2.Data());
    leg1->Draw();
}

void DrawAutoGammaHistMatch3H( TH1* histo1,
					 TH1* histo2,
					 TH1* histo3,
					 TString Title, TString XTitle, TString YTitle,
					 Bool_t YRangeMax, Double_t YMaxFactor, Double_t YMinimum,
					 Bool_t YRange, Double_t YMin ,Double_t YMax,
					 Bool_t XRange, Double_t XMin, Double_t XMax,
					 Float_t xOffset, Float_t yOffset, TString legS1, TString legS2, TString legS3, Bool_t fillStyle) {
	if (YRangeMax && !XRange){
		YRange = kFALSE;
		Double_t maxRangeR = histo1->GetMaximum();
		if(maxRangeR < histo2->GetMaximum()){
			maxRangeR = histo2->GetMaximum();
		}
		Double_t minRangeR = histo1->GetMinimum();
		if(minRangeR > histo2->GetMinimum()){
			minRangeR = histo2->GetMinimum();
		}
		if(YMinimum > minRangeR){minRangeR = YMinimum;}
		histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);
	}
	if (YRangeMax && XRange){
		YRange = kFALSE;
		Double_t maxRangeR = histo1->GetMaximum();
		if(maxRangeR < histo2->GetMaximum()){
			maxRangeR = histo2->GetMaximum();
		}
		Double_t minRangeR = histo1->GetMinimum();
		if(minRangeR > histo2->GetMinimum()){
			minRangeR = histo2->GetMinimum();
		}
		if(YMinimum > minRangeR){minRangeR = YMinimum;}
		histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);
	}
	if (YRange && XRange){
		histo1->GetYaxis()->SetRangeUser(YMin, YMax);
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);
	}
	if (!YRangeMax && !YRange && XRange){
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);
	}
	if (YRange && !XRange){
		histo1->GetYaxis()->SetRangeUser(YMin, YMax);
	}

	if(Title.CompareTo("") != 0){
		histo1->SetTitle(Title.Data());
	}else{histo1->SetTitle("");
	histo2->SetTitle("");}
	if(XTitle.CompareTo("") != 0){
		histo1->SetXTitle(XTitle.Data());
	}
	if(YTitle.CompareTo("") != 0){
		histo1->SetYTitle(YTitle.Data());
	}

    DrawGammaSetMarker(histo1, 24, 0.5, kBlack, kBlack);
    DrawGammaSetMarker(histo2, 1, 0.5, kRed+1, kRed+1);
    DrawGammaSetMarker(histo3, 1, 0.5, kMagenta+2, kMagenta+2);

	histo1->GetYaxis()->SetLabelFont(42);
	histo1->GetXaxis()->SetLabelFont(42);
	histo1->GetYaxis()->SetTitleFont(62);
	histo1->GetXaxis()->SetTitleFont(62);

	histo1->GetYaxis()->SetLabelSize(0.035);
	histo1->GetYaxis()->SetTitleSize(0.04);
	histo1->GetYaxis()->SetDecimals();
	histo1->GetYaxis()->SetTitleOffset(yOffset);
	histo1->GetXaxis()->SetTitleSize(0.04);
	histo1->GetXaxis()->SetLabelSize(0.035);
	histo1->GetXaxis()->SetTitleOffset(xOffset);
	histo1->DrawCopy("AXIS");
	if(fillStyle)
	{
		histo2->SetFillColor(kRed+1);
		histo2->SetFillStyle(3003);
	}
	histo2->DrawCopy("hist,same");
	histo3->DrawCopy("hist,same");
	histo1->DrawCopy("p,hist,same");
	TLegend* leg1;
	if(legS2.BeginsWith("TrueMC") && legS3.BeginsWith("TrueMC"))
	{
		leg1 = new TLegend( 0.73,0.73,0.97,0.93);
	}else
	{
		leg1 = new TLegend( 0.77,0.79,0.97,0.93);
	}
	leg1->SetTextSize(0.04);
	leg1->SetFillColor(0);
	leg1->AddEntry(histo1,legS1.Data(),"p");
	leg1->AddEntry(histo2,legS2.Data());
	leg1->AddEntry(histo3,legS3.Data());
	leg1->Draw();
}

void PlotCutHistoReasons(TCanvas *canvas, Double_t leftMargin, Double_t rightMargin, Double_t topMargin, Double_t bottomMargin,
                         TH2D* src, TString xLabel, TString yLabel, Double_t factorLow, Double_t factorHigh,
                         Int_t fixRange, Double_t range){
  const Int_t fNBinsPt 			= 35;
  Double_t fBinsPt[fNBinsPt+1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8,
                                       2.0, 2.4, 2.8, 3.2, 3.6,
                                       4.0, 4.4, 4.8, 5.2, 5.6,
                                       6.0, 7.0, 8.0,
                                       9.0, 10, 12, 14, 16, 18,
                                       20, 25, 30, 40, 50, 75, 100};

  TH1D* fDeltaPt = new TH1D("deltaPt","",fNBinsPt,fBinsPt);
  for(Int_t iPt=1;iPt<fNBinsPt+1;iPt++){
      fDeltaPt->SetBinContent(iPt,fBinsPt[iPt]-fBinsPt[iPt-1]);
      fDeltaPt->SetBinError(iPt,0);
  }

  canvas->cd();
  canvas->SetLeftMargin(leftMargin);
  canvas->SetRightMargin(rightMargin);
  canvas->SetTopMargin(topMargin);
  if(topMargin<0.06) canvas->SetTopMargin(0.06);
  canvas->SetBottomMargin(bottomMargin);

  Int_t nColumns=0;
  for(Int_t i=1; i<=src->GetNbinsX(); i++)
  {
    TString tempStr = src->GetXaxis()->GetBinLabel(i);
    TH1D* hist = (TH1D*) src->ProjectionY(Form("projectSrc_%i",i),i,i);
    if(tempStr.CompareTo("") && tempStr.CompareTo("out") && hist->GetEntries()>0) nColumns++;
    delete hist;
  }

  TLegend *legend = new TLegend(0.15,0.95,0.9,0.98);
  legend->SetNColumns(nColumns);
  legend->SetFillColor(0);
  legend->SetLineColor(0);
  legend->SetTextSize(0.03);
  legend->SetTextFont(42);

  std::vector<TH1D*> vecHists;
  for(Int_t i=1; i<=src->GetNbinsX(); i++)
  {
    TString tempStr = src->GetXaxis()->GetBinLabel(i);
    if(tempStr.CompareTo("") && tempStr.CompareTo("out")){
      TH1D* hist = (TH1D*) src->ProjectionY(Form("projectSrc_%i",i),i,i);
      if(hist->GetEntries()<1.){
        delete hist;
        continue;
      }

      hist->Sumw2();
      TH1D* histRebin = (TH1D*) hist->Rebin(fNBinsPt, Form("rebinProjectSrc_%i",i), fBinsPt);
      histRebin->Divide(fDeltaPt);
      delete hist;

      Int_t nMarker = 20 + (i % 15);
      Int_t nColor = (i + 1) % 45;
      if(nColor == 10){nColor = 45;} if(nColor == 18){nColor = 46;} if(nColor == 19){nColor = 47;}
      DrawGammaSetMarker(histRebin, nMarker, 0.8, nColor, nColor);

      histRebin->SetTitle("");
      histRebin->SetXTitle(xLabel.Data());
      histRebin->SetYTitle(yLabel.Data());

      histRebin->GetYaxis()->SetLabelFont(42);
      histRebin->GetXaxis()->SetLabelFont(42);
      histRebin->GetYaxis()->SetTitleFont(62);
      histRebin->GetXaxis()->SetTitleFont(62);

      histRebin->GetYaxis()->SetLabelSize(0.035);
      histRebin->GetYaxis()->SetTitleSize(0.04);
      histRebin->GetYaxis()->SetTitleOffset(1.0);
      histRebin->GetYaxis()->SetDecimals();

      histRebin->GetXaxis()->SetLabelSize(0.035);
      histRebin->GetXaxis()->SetTitleSize(0.04);
      histRebin->GetXaxis()->SetTitleOffset(1.1);

      legend->AddEntry(histRebin,tempStr.Data(),"p");
      vecHists.push_back(histRebin);
    }
  }

  Int_t minB=0;
  Int_t maxB=1;
  GetMinMaxBin(vecHists,minB,maxB);
  AdjustHistRange(vecHists,factorLow,factorHigh,kFALSE,fixRange,range);

  for(Int_t i=0; i<(Int_t)vecHists.size(); i++)
  {
      SetXRange(vecHists.at(i),minB,maxB);
      TString draw = (i==0)?"p":"p, same";
      vecHists.at(i)->DrawCopy(draw.Data());
  }
  legend->Draw("same");

  vecHists.clear();
  delete fDeltaPt;

  return;
}

void CalculateFractionMatches(TH1D* cFraction, TH2F* cESD_Mother, TH2F* cESD_Mother_Matched, Int_t cBin, Double_t ptMin, Double_t ptMax)
{
	cESD_Mother->GetXaxis()->SetRangeUser(0,0.8);
	cESD_Mother->GetYaxis()->SetRangeUser(ptMin,ptMax);
	cESD_Mother_Matched->GetXaxis()->SetRangeUser(0,0.8);
	cESD_Mother_Matched->GetYaxis()->SetRangeUser(ptMin,ptMax);
	Double_t nESD_Mother = cESD_Mother->GetEffectiveEntries();
	Double_t nESD_Mother_Matched = cESD_Mother_Matched->GetEffectiveEntries();
	if(nESD_Mother_Matched > 0){
		Double_t nFractionMatches = nESD_Mother_Matched / (nESD_Mother + nESD_Mother_Matched);
		cFraction->SetBinContent(cBin, nFractionMatches);
		Double_t A = nESD_Mother_Matched; Double_t dA = sqrt(A);
		Double_t B = nESD_Mother; Double_t dB = sqrt(B);
		Double_t AB = A+B; Double_t AB2 = pow(AB,2);
		cFraction->SetBinError(cBin, sqrt( pow(A*dB / AB2,2) + pow(dA*((1/AB)-(A/AB2)),2) ));
	}else{
        cout << "WARNING: for bin " << cBin << ", range: " << ptMin << " - " << ptMax << ", nESD_Mother_Matched equals zero!" << endl;
		cFraction->SetBinContent(cBin, 0.0001);
		cFraction->SetBinError(cBin, 0);
	}

	return;
}

void DrawAutoGammaHistoPaper( TH1* histo1,
					TString Title, TString XTitle, TString YTitle,
					Bool_t YRangeMax, Double_t YMaxFactor, Double_t YMinimum,
					Bool_t YRange, Double_t YMin ,Double_t YMax,
                    Bool_t XRange, Double_t XMin, Double_t XMax, Double_t yOffset, Double_t yTitleSize) {
	if (YRangeMax && !XRange){
		YRange = kFALSE;
		Double_t maxRangeR = histo1->GetMaximum();
		Double_t minRangeR = histo1->GetMinimum();
		if(YMinimum > minRangeR){minRangeR = YMinimum;}
		histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);
	}
	if (YRangeMax && XRange){
		YRange = kFALSE;
		Double_t maxRangeR = histo1->GetMaximum();
		Double_t minRangeR = histo1->GetMinimum();
		if(YMinimum > minRangeR){minRangeR = YMinimum;}
		histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);
	}
	if (YRange && XRange){
		histo1->GetYaxis()->SetRangeUser(YMin, YMax);
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);
	}
	if (!YRangeMax && !YRange && XRange){
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);
	}

	if (YRange && !XRange){
		histo1->GetYaxis()->SetRangeUser(YMin, YMax);
	}


	histo1->SetTitle(Title.Data());

	if(XTitle.CompareTo("") != 0){
		histo1->SetXTitle(XTitle.Data());
	}
	if(YTitle.CompareTo("") != 0){
		histo1->SetYTitle(YTitle.Data());
	}

	histo1->GetYaxis()->SetLabelFont(42);
	histo1->GetXaxis()->SetLabelFont(42);
	histo1->GetYaxis()->SetTitleFont(62);
	histo1->GetXaxis()->SetTitleFont(62);
	histo1->GetYaxis()->SetLabelSize(0.051);
    histo1->GetYaxis()->SetTitleSize(yTitleSize);
	histo1->GetYaxis()->SetDecimals();
	histo1->GetXaxis()->SetTitleOffset(1);
	histo1->GetYaxis()->SetTitleOffset(yOffset);
	histo1->GetXaxis()->SetTitleSize(0.06);
	histo1->GetXaxis()->SetLabelSize(0.051);
	histo1->DrawCopy("e,hist");
}

void DrawAutoGammaHistoPaper2D( TH2* histo1,
					TString Title, TString XTitle, TString YTitle,
					Bool_t YRangeMax, Double_t YMaxFactor, Double_t YMinimum,
					Bool_t YRange, Double_t YMin ,Double_t YMax,
					Bool_t XRange, Double_t XMin, Double_t XMax, Double_t xOffset, Double_t yOffset) {
	if (YRangeMax && !XRange){
		YRange = kFALSE;
		Double_t maxRangeR = histo1->GetMaximum();
		Double_t minRangeR = histo1->GetMinimum();
		if(YMinimum > minRangeR){minRangeR = YMinimum;}
		histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);
	}
	if (YRangeMax && XRange){
		YRange = kFALSE;
		Double_t maxRangeR = histo1->GetMaximum();
		Double_t minRangeR = histo1->GetMinimum();
		if(YMinimum > minRangeR){minRangeR = YMinimum;}
		histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);
	}
	if (YRange && XRange){
		histo1->GetYaxis()->SetRangeUser(YMin, YMax);
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);
	}
	if (!YRangeMax && !YRange && XRange){
		histo1->GetXaxis()->SetRangeUser(XMin, XMax);
	}

	if (YRange && !XRange){
		histo1->GetYaxis()->SetRangeUser(YMin, YMax);
	}


	histo1->SetTitle(Title.Data());

	if(XTitle.CompareTo("") != 0){
		histo1->SetXTitle(XTitle.Data());
	}
	if(YTitle.CompareTo("") != 0){
		histo1->SetYTitle(YTitle.Data());
	}

	histo1->GetYaxis()->SetLabelFont(42);
	histo1->GetXaxis()->SetLabelFont(42);
	histo1->GetYaxis()->SetTitleFont(62);
	histo1->GetXaxis()->SetTitleFont(62);
	histo1->GetYaxis()->SetLabelSize(0.051);
	histo1->GetYaxis()->SetTitleSize(0.06);
	histo1->GetYaxis()->SetDecimals();
	histo1->GetXaxis()->SetTitleOffset(xOffset);
	histo1->GetYaxis()->SetTitleOffset(yOffset);
	histo1->GetXaxis()->SetTitleSize(0.06);
	histo1->GetXaxis()->SetLabelSize(0.051);
}

void DrawArmenterosLines(Bool_t drawLines = kTRUE){
	TLatex *k0s = new TLatex(0.48,0.73,"K^{0}_{s}");
	k0s->SetNDC();
	k0s->SetTextFont(62);
	k0s->SetTextSize(0.05);
	k0s->SetLineWidth(4);
	k0s->Draw();

    TLatex *lambda = new TLatex(0.71,0.47,"#Lambda");
	lambda->SetNDC();
	lambda->SetTextFont(62);
	lambda->SetTextSize(0.05);
	lambda->SetLineWidth(4);
	lambda->Draw();

    TLatex *lambdabar = new TLatex(0.26,0.47,"#bar{#Lambda}");
	lambdabar->SetNDC();
	lambdabar->SetTextFont(62);
	lambdabar->SetTextSize(0.05);
	lambdabar->SetLineWidth(4);
	lambdabar->Draw();

    TLatex *gamma = new TLatex(0.495,0.16,"#gamma");
	gamma->SetNDC();
	gamma->SetTextFont(62);
	gamma->SetTextSize(0.05);
	gamma->SetLineWidth(4);
	gamma->Draw();

    if(drawLines){
      TF1 *felipseGamma = new TF1 ("felipseGamma", "TMath::Sqrt(1-(x*x)/([0]*[0]))*[1]", -0.95, 0.95);
      felipseGamma->SetParameter(0,0.95);
      felipseGamma->SetParameter(1,0.03);
      felipseGamma->SetLineColor(kBlack);
      felipseGamma->Draw("same");
      TF1 *felipseK0s = new TF1 ("felipseK0s", "TMath::Sqrt(1-(x*x)/([0]*[0]))*[1]", -0.828, 0.828);
      felipseK0s->SetParameter(0,0.828);
      felipseK0s->SetParameter(1,0.206);
      felipseK0s->SetLineColor(kViolet+2);
      felipseK0s->Draw("same");
      TF1 *felipseLambda = new TF1 ("felipseLambda", "TMath::Sqrt(1-((x-[2])*(x-[2]))/([0]*[0]))*[1]", 0.511, 0.871);
      felipseLambda->SetParameter(0,0.18);
      felipseLambda->SetParameter(1,0.101);
      felipseLambda->SetParameter(2,0.691);
      felipseLambda->SetLineColor(kGreen+2);
      felipseLambda->Draw("same");
      TF1 *felipseLambdaBar = new TF1 ("felipseLambdaBar", "TMath::Sqrt(1-((x-[2])*(x-[2]))/([0]*[0]))*[1]", -0.871, -0.511);
      felipseLambdaBar->SetParameter(0,0.18);
      felipseLambdaBar->SetParameter(1,0.101);
      felipseLambdaBar->SetParameter(2,-0.691);
      felipseLambdaBar->SetLineColor(kGreen+2);
      felipseLambdaBar->Draw("same");
    }
	return;
}

void PlotCaloQAModule(TH2* src, Int_t nCaloModules, TString xLabel, TString yLabel, Double_t factorLow, Double_t factorHigh, Int_t fixRange, Double_t range, Bool_t doXRange, Int_t minXRange, Int_t maxXRange, Bool_t kEnergyPlot,
					  Double_t titleOffsetX, Double_t titleOffsetY){
	const Int_t fNBinsClusterPt 			= 60;
	Double_t fBinsClusterPt[fNBinsClusterPt+1] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
										 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
										 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8,
										 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8,
										 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.4, 7.8, 8.2, 8.6,
										 9.0, 9.5, 10, 11, 12, 14, 16, 18, 20, 25, 30};

	TH1D* fDeltaPt = 0;
	if(kEnergyPlot){
		fDeltaPt = new TH1D("deltaPt","",fNBinsClusterPt,fBinsClusterPt);
		for(Int_t iPt=1;iPt<fNBinsClusterPt+1;iPt++){
			fDeltaPt->SetBinContent(iPt,fBinsClusterPt[iPt]-fBinsClusterPt[iPt-1]);
			fDeltaPt->SetBinError(iPt,0);
		}
	}

	TLegend *legend = new TLegend(0.12,0.95,0.9,0.98);
	legend->SetNColumns(nCaloModules);
	legend->SetFillColor(0);
	legend->SetLineColor(0);
	legend->SetTextSize(0.03);
	legend->SetTextFont(42);

	std::vector<TH1D*> vecHists;
	for(Int_t i=0; i<nCaloModules; i++)
	{
		TH1D* hist = (TH1D*) src->ProjectionX(Form("projectSrc_%i",i),i+1,i+1);

		if(kEnergyPlot){
			hist = (TH1D*) hist->Rebin(fNBinsClusterPt, Form("rebinProjectSrc_%i",i), fBinsClusterPt);
			hist->Divide(fDeltaPt);
		}

		Int_t nMarker = 20 + (i % 15);
		Int_t nColor = (i + 1) % 45;
		if(nColor == 10){nColor = 45;} if(nColor == 18){nColor = 46;} if(nColor == 19){nColor = 47;}
		DrawGammaSetMarker(hist, nMarker, 0.8, nColor, nColor);

        hist->SetTitle("");
        if(doXRange) SetXRange(hist,minXRange,maxXRange);
		hist->SetXTitle(xLabel.Data());
		hist->SetYTitle(yLabel.Data());

		hist->GetYaxis()->SetLabelFont(42);
		hist->GetXaxis()->SetLabelFont(42);
		hist->GetYaxis()->SetTitleFont(62);
		hist->GetXaxis()->SetTitleFont(62);

		hist->GetYaxis()->SetLabelSize(0.035);
        hist->GetYaxis()->SetTitleSize(0.04);
		hist->GetYaxis()->SetTitleOffset(titleOffsetY);
		hist->GetYaxis()->SetDecimals();

		hist->GetXaxis()->SetLabelSize(0.035);
        hist->GetXaxis()->SetTitleSize(0.04);
		hist->GetXaxis()->SetTitleOffset(titleOffsetX);

		hist->Sumw2();
		legend->AddEntry(hist,Form("SM: %i",i),"p");
		vecHists.push_back(hist);
	}

	AdjustHistRange(vecHists,factorLow,factorHigh,kTRUE,fixRange,range);

	for(Int_t i=0; i<(Int_t)vecHists.size(); i++)
	{
		TString draw = (i==0)?"p":"p, same";
		vecHists.at(i)->DrawCopy(draw.Data());
	}
	legend->Draw("same");

	vecHists.clear();
	delete fDeltaPt;

	return;
}

void PlotCellMeanVsSigma(CellQAObj* obj, Int_t nCaloCells, TH2* hist, TString xLabel, TString yLabel, Bool_t XRange, Float_t XMin, Float_t XMax, Bool_t YRange, Float_t YMin, Float_t YMax, Bool_t kEnergy, Bool_t kMC, Double_t titleOffsetX, Double_t titleOffsetY)
{
	TH2D* outHist;
	if(kEnergy){
		if(obj) outHist = new TH2D("CellEnergy","CellEnergy",(obj->EnergyMean[1]+0.15)*200,0,obj->EnergyMean[1]+0.15,(obj->EnergySigma[1]+0.15)*200,0,obj->EnergySigma[1]+0.15);
		else outHist = new TH2D("CellEnergy","CellEnergy",100,0,0.5,80,0,0.4);
	}else{
		if(kMC) outHist = new TH2D("CellTime","CellTime", 80, 6E-7,6.4E-7, 40, 0, 0.2E-7);
		else{
            if(obj){
              Double_t binX = 1E9*((obj->TimeMean[1]+0.02E-6)-(obj->TimeMean[0]-0.02E-6));
              Double_t binY = 4*1E8*((obj->TimeSigma[1]+0.05E-6)-(obj->TimeSigma[0]-0.05E-6));
              if(binX>1E3){cout << "WARNING: PlotCellMeanVsSigma: CellTime, number of bins for x-axis'" << binX << "', limiting to 1000" << endl; binX = 1E3;}
              if(binY>1E3){cout << "WARNING: PlotCellMeanVsSigma: CellTime, number of bins for y-axis'" << binY << "', limiting to 1000" << endl; binY = 1E3;}
              outHist = new TH2D("CellTime","CellTime", binX, (obj->TimeMean[0]-0.02E-6),(obj->TimeMean[1]+0.02E-6), binY, (obj->TimeSigma[0]-0.05E-6), (obj->TimeSigma[1]+0.05E-6));
            }
			else outHist = new TH2D("CellTime","CellTime", 100, -0.04E-6,0.04E-6, 50, 0, 0.1E-6);
		}
	}

	if(XRange) outHist->GetXaxis()->SetRangeUser(XMin, XMax);
	if(YRange) outHist->GetYaxis()->SetRangeUser(YMin, YMax);
    outHist->SetTitle("");
	outHist->SetXTitle(xLabel.Data());
	outHist->SetYTitle(yLabel.Data());

	outHist->GetYaxis()->SetLabelFont(42);
	outHist->GetXaxis()->SetLabelFont(42);
	outHist->GetYaxis()->SetTitleFont(62);
	outHist->GetXaxis()->SetTitleFont(62);

	outHist->GetYaxis()->SetLabelSize(0.035);
    outHist->GetYaxis()->SetTitleSize(0.04);
	outHist->GetYaxis()->SetTitleOffset(titleOffsetY);
	outHist->GetYaxis()->SetDecimals();

	outHist->GetXaxis()->SetLabelSize(0.035);
    outHist->GetXaxis()->SetTitleSize(0.04);
	outHist->GetXaxis()->SetTitleOffset(titleOffsetX);

	for(Int_t iY=1; iY<nCaloCells+1; iY++)
	{
		TH1D* temp;
		if(kEnergy) temp = (TH1D*) hist->ProjectionX("energy",iY,iY);
		else temp = (TH1D*) hist->ProjectionX("time",iY,iY);
		Double_t mean = temp->GetMean();
		if(mean>=outHist->GetXaxis()->GetBinUpEdge(outHist->GetNbinsX())){ mean = outHist->GetXaxis()->GetBinCenter(outHist->GetNbinsX());}
		Double_t rms = temp->GetRMS();
		if(rms>=outHist->GetYaxis()->GetBinUpEdge(outHist->GetNbinsY())){ rms = outHist->GetYaxis()->GetBinCenter(outHist->GetNbinsY());}
		if(temp->Integral(1,temp->GetNbinsX())>0){
			outHist->Fill(mean,rms);
			if(obj){
				if(kEnergy && ( mean < obj->EnergyMean[0] || mean > obj->EnergyMean[1] || rms < obj->EnergySigma[0] || rms > obj->EnergySigma[1] )){
					if(!CheckGoodCell(obj,iY-1)) obj->cellIDsEnergy.push_back(iY-1);
				}
				if(!kEnergy && ( mean < obj->TimeMean[0] || mean > obj->TimeMean[1] || rms < obj->TimeSigma[0] || rms > obj->TimeSigma[1]) ){
					if(!CheckGoodCell(obj,iY-1)) obj->cellIDsTime.push_back(iY-1);
				}
			}
		}
		delete temp;
	}

	outHist->DrawCopy("COLZ");
	delete outHist;
	return;
}

TH2D** PlotCellMeanVsSigmaForRunwise(Int_t nCaloCells, TH2* histEnergy, TH2* histTime, TString xLabel, TString yLabel, TString xLabelT, TString yLabelT, Bool_t kMC, Double_t titleOffsetX, Double_t titleOffsetY)
{
	TH2D** outHist = new TH2D*[2];
	Double_t xMinMaxEnergy[2], xMinMaxTime[2];
	Double_t yMinMaxEnergy[2], yMinMaxTime[2];

	xMinMaxEnergy[0]=0; xMinMaxEnergy[1]=0.5;
	yMinMaxEnergy[0]=0; yMinMaxEnergy[1]=0.5;
	outHist[0] = new TH2D("CellEnergy","CellEnergy",100,xMinMaxEnergy[0],xMinMaxEnergy[1],100,yMinMaxEnergy[0],yMinMaxEnergy[1]);

	if(kMC){
		xMinMaxTime[0]=6E-7; xMinMaxTime[1]=6.4E-7;
		yMinMaxTime[0]=0; yMinMaxTime[1]=1E-7;
		outHist[1] = new TH2D("CellTime","CellTime", 100,xMinMaxTime[0],xMinMaxTime[1],100,yMinMaxTime[0],yMinMaxTime[1]);
	}
	else{
		xMinMaxTime[0]=-0.2E-6; xMinMaxTime[1]=0.2E-6;
		yMinMaxTime[0]=0; yMinMaxTime[1]=0.2E-6;
		outHist[1] = new TH2D("CellTime","CellTime", 100,xMinMaxTime[0],xMinMaxTime[1],100,yMinMaxTime[0],yMinMaxTime[1]);
	}

	outHist[0]->SetXTitle(xLabel.Data());
	outHist[0]->SetYTitle(yLabel.Data());
	outHist[1]->SetXTitle(xLabelT.Data());
	outHist[1]->SetYTitle(yLabelT.Data());

	for(Int_t i=0; i<2; i++){
        outHist[i]->SetTitle("");
		outHist[i]->GetYaxis()->SetLabelFont(42);
		outHist[i]->GetXaxis()->SetLabelFont(42);
		outHist[i]->GetYaxis()->SetTitleFont(62);
		outHist[i]->GetXaxis()->SetTitleFont(62);

		outHist[i]->GetYaxis()->SetLabelSize(0.035);
        outHist[i]->GetYaxis()->SetTitleSize(0.04);
		outHist[i]->GetYaxis()->SetTitleOffset(titleOffsetY);
		outHist[i]->GetYaxis()->SetDecimals();

		outHist[i]->GetXaxis()->SetLabelSize(0.035);
        outHist[i]->GetXaxis()->SetTitleSize(0.04);
		outHist[i]->GetXaxis()->SetTitleOffset(titleOffsetX);
	}

	for(Int_t iY=1; iY<nCaloCells+1; iY++)
	{
		TH1D* temp[2];
		temp[0] = (TH1D*) histEnergy->ProjectionX("energy",iY,iY);
		temp[1] = (TH1D*) histTime->ProjectionX("time",iY,iY);

		Double_t meanE = temp[0]->GetMean();
		if(meanE>=xMinMaxEnergy[1]){ meanE = outHist[0]->GetXaxis()->GetBinCenter(outHist[0]->GetNbinsX());}
		Double_t rmsE = temp[0]->GetRMS();
		if(rmsE>=yMinMaxEnergy[1]){ rmsE = outHist[0]->GetYaxis()->GetBinCenter(outHist[0]->GetNbinsY());}

		Double_t meanT = temp[1]->GetMean();
		if(meanT>=xMinMaxTime[1]){ meanT = outHist[1]->GetXaxis()->GetBinCenter(outHist[1]->GetNbinsX());}
		Double_t rmsT = temp[1]->GetRMS();
		if(rmsT>=yMinMaxTime[1]){ rmsT = outHist[1]->GetYaxis()->GetBinCenter(outHist[1]->GetNbinsY());}

		if(temp[0]->Integral(1,temp[0]->GetNbinsX())>0) outHist[0]->Fill(meanE,rmsE);
		if(temp[1]->Integral(1,temp[1]->GetNbinsX())>0) outHist[1]->Fill(meanT,rmsT);

		delete temp[0];
		delete temp[1];
	}
	outHist[0]->SetOption("COLZ");
	outHist[1]->SetOption("COLZ");
	return outHist;
}



void PlotHotCells(CellQAObj* obj, Int_t iSw, Int_t nCaloCells, TH2* hist, TString xLabel, TString yLabel, Bool_t XRange, Float_t XMin, Float_t XMax, Bool_t YRange, Float_t YMin, Float_t YMax, Bool_t kMC, Double_t titleOffsetX, Double_t titleOffsetY)
{
	if(iSw==0){
		Double_t min = 0;
		Double_t max = 0;
		TH1D* outHist;
		if(kMC) {
			min = -0.05E6;
			max = 0.05E6;
			outHist = new TH1D("CellHotCells","", 100,min,max);
		}
		else{
			if(obj){
				min = obj->HotCells1D[0];
				max = obj->HotCells1D[1];
				if(min-(max/10) < -10) outHist = new TH1D("CellHotCells","", 1000,-10,max+(max/10));
				else outHist = new TH1D("CellHotCells","", 1000,min-(max/10),max+(max/10));
			}else{
				min = -0.5E6;
				max = 5E6;
				outHist = new TH1D("CellHotCells","", 1000,min,max);
			}
		}
		//SetLogBinningXTH(outHist);
		if(XRange) outHist->GetXaxis()->SetRangeUser(XMin, XMax);
		if(YRange) outHist->GetYaxis()->SetRangeUser(YMin, YMax);
        outHist->SetTitle("");
		outHist->SetXTitle(xLabel.Data());
		outHist->SetYTitle(yLabel.Data());

		outHist->GetYaxis()->SetLabelFont(42);
		outHist->GetXaxis()->SetLabelFont(42);
		outHist->GetYaxis()->SetTitleFont(62);
		outHist->GetXaxis()->SetTitleFont(62);

		outHist->GetYaxis()->SetLabelSize(0.035);
        outHist->GetYaxis()->SetTitleSize(0.04);
		outHist->GetYaxis()->SetTitleOffset(titleOffsetY);
		outHist->GetYaxis()->SetDecimals();

		outHist->GetXaxis()->SetLabelSize(0.035);
        outHist->GetXaxis()->SetTitleSize(0.04);
		outHist->GetXaxis()->SetTitleOffset(titleOffsetX);

//		Double_t averageTemp=0;
//		Double_t sigmaTemp=0;
//		Double_t nCells=0;
//		for(Int_t iY=1; iY<nCaloCells+1; iY++)
//		{
//			TH1D* temp = (TH1D*) hist->ProjectionX("number",iY,iY);
//			if(temp->Integral(1,temp->GetNbinsX())>10) {
//				averageTemp+=temp->Integral(1,temp->GetNbinsX());
//				sigmaTemp+=TMath::Power(temp->Integral(1,temp->GetNbinsX()),2);
//				nCells++;
//			}
//			delete temp;
//		}
//		Double_t average=averageTemp/nCells;
//		Double_t sigma=sqrt(1/(nCells-1)*(sigmaTemp-1/nCells*TMath::Power(averageTemp,2)));

		for(Int_t iY=1; iY<nCaloCells+1; iY++)
		{
			TH1D* temp = (TH1D*) hist->ProjectionX("numberOfEntries",iY,iY);
			if(temp->Integral(1,temp->GetNbinsX())>0) {
				Double_t nFired = temp->Integral(1,temp->GetXaxis()->GetNbins());
				if(!kMC && obj && nFired>=max) nFired = max;
				else if(nFired>=(max-1)) nFired = max-1;
				outHist->Fill(nFired);
				if(obj && ( nFired <= obj->HotCells1D[0] || nFired >= obj->HotCells1D[1] )){
					if(!CheckGoodCell(obj,iY-1)) obj->cellIDsHotCells1D.push_back(iY-1);
				}
			}
			delete temp;
			//if(temp->GetMean()>0.4) cout << iY << endl;
		}

		outHist->DrawCopy();
		delete outHist;
		return;
	}else if(iSw==1){
		TH2D* outHist2D = new TH2D("CellHotCells2D","", 100,0.00009,110, 9,0.1,1.);
		SetLogBinningXTH(outHist2D);
		if(XRange) outHist2D->GetXaxis()->SetRangeUser(XMin, XMax);
		if(YRange) outHist2D->GetYaxis()->SetRangeUser(YMin, YMax);
        outHist2D->SetTitle("");
		outHist2D->SetXTitle(xLabel.Data());
		outHist2D->SetYTitle(yLabel.Data());

		outHist2D->GetYaxis()->SetLabelFont(42);
		outHist2D->GetXaxis()->SetLabelFont(42);
		outHist2D->GetYaxis()->SetTitleFont(62);
		outHist2D->GetXaxis()->SetTitleFont(62);

		outHist2D->GetYaxis()->SetLabelSize(0.035);
        outHist2D->GetYaxis()->SetTitleSize(0.04);
		outHist2D->GetYaxis()->SetTitleOffset(titleOffsetY);
		outHist2D->GetYaxis()->SetDecimals();

		outHist2D->GetXaxis()->SetLabelSize(0.035);
        outHist2D->GetXaxis()->SetTitleSize(0.04);
		outHist2D->GetXaxis()->SetTitleOffset(titleOffsetX);

		for(Int_t iY=1; iY<nCaloCells+1; iY++)
		{
			TH1D* temp = (TH1D*) hist->ProjectionX("numberOfEntriesAbove",iY,iY);
			for(Int_t iBin=1; iBin<10; iBin++)
			{
				if(temp->Integral(1,temp->GetNbinsX())>0){
					//Double_t nTotal = temp->Integral(1,temp->GetXaxis()->GetNbins());
					Double_t nCell = temp->Integral(1,iBin);
					Double_t nCellAbove = temp->Integral(iBin+1,temp->GetXaxis()->GetNbins());
					Double_t fraction = 100;
					if(nCell>0) fraction = nCellAbove / nCell;
					if(fraction>100) fraction = 100;
					if(fraction<0.0001) fraction = 0.0001;
					outHist2D->Fill(fraction,((Double_t)iBin+0.05)/10.);
					if(obj && (fraction < obj->HotCells2D[iBin-1][0] || fraction > obj->HotCells2D[iBin-1][1]) ){
						if(!CheckGoodCell(obj,iY-1)) obj->cellIDsHotCells2D.push_back(iY-1);
					}
					//if(fraction>10){ cout << "iBin:" << iBin << ", cell" << iY << endl;}
				}
			}
			delete temp;
		}

		outHist2D->DrawCopy("COLZ");
		delete outHist2D;
		return;
	}else if(iSw==2){
		TH1D* outHist1D;
		if(obj && (obj->HotCellsTime1D[1]+0.7 > 1)) outHist1D = new TH1D("CellTime1D","", ((obj->HotCellsTime1D[1]+0.7)-1)*100,1,(obj->HotCellsTime1D[1]+0.7));
		else outHist1D = new TH1D("CellTime1D","", 110,1,2.1);
		//SetLogBinningXTH(outHist1D);
		if(XRange) outHist1D->GetXaxis()->SetRangeUser(XMin, XMax);
		if(YRange) outHist1D->GetYaxis()->SetRangeUser(YMin, YMax);
        outHist1D->SetTitle("");
		outHist1D->SetXTitle(xLabel.Data());
		outHist1D->SetYTitle(yLabel.Data());

		outHist1D->GetYaxis()->SetLabelFont(42);
		outHist1D->GetXaxis()->SetLabelFont(42);
		outHist1D->GetYaxis()->SetTitleFont(62);
		outHist1D->GetXaxis()->SetTitleFont(62);

		outHist1D->GetYaxis()->SetLabelSize(0.035);
        outHist1D->GetYaxis()->SetTitleSize(0.04);
		outHist1D->GetYaxis()->SetTitleOffset(titleOffsetY);
		outHist1D->GetYaxis()->SetDecimals();

		outHist1D->GetXaxis()->SetLabelSize(0.035);
        outHist1D->GetXaxis()->SetTitleSize(0.04);
		outHist1D->GetXaxis()->SetTitleOffset(titleOffsetX);

		Double_t fracLimit = 0;
		if(obj) fracLimit = obj->HotCellsTime1D[1]+0.6;
		else fracLimit = 2;

		for(Int_t iY=1; iY<nCaloCells+1; iY++)
		{
			TH1D* temp = (TH1D*) hist->ProjectionX("numberOfEntriesAbove",iY,iY);
			if(temp->Integral(1,temp->GetNbinsX())>0){
				Double_t nTotal = temp->Integral(1,temp->GetXaxis()->GetNbins());
				Double_t nCell = 0;
				Double_t localMaximum = temp->GetBinCenter(temp->GetMaximumBin());
				nCell = temp->Integral(temp->GetXaxis()->FindBin(localMaximum-0.02E-6),temp->GetXaxis()->FindBin(localMaximum+0.02E-6));
				Double_t fraction = fracLimit;
				if(nCell>0) fraction = nTotal / nCell;
				if(fraction>fracLimit) fraction = fracLimit;
				outHist1D->Fill(fraction);
				if(obj && ( fraction < obj->HotCellsTime1D[0] || fraction > obj->HotCellsTime1D[1] )){
					if(!CheckGoodCell(obj,iY-1)) obj->cellIDsHotCellsTime1D.push_back(iY-1);
				}
				//cout << nTotal << ", " << nCell << endl;
			}
			delete temp;
		}

		outHist1D->DrawCopy();
		delete outHist1D;
		return;
	}
	else{cout << "ERROR: iSw not defined: " << iSw << endl; return;}
}

void CheckCellsDataMC(CellQAObj* obj, TH2* fHistDataCell, TH2* fHistMCCell, TString xLabel, TString yLabel, Int_t nCaloCells, TString plotData, TString plotMC){
	TH2I* compare = new TH2I(Form("%sVs%s",plotData.Data(),plotMC.Data()),"",3,-1,2,nCaloCells+1000,0,nCaloCells+1000);
	compare->GetXaxis()->SetBinLabel(1,"onlyMC");
	compare->GetXaxis()->SetBinLabel(2,"ok");
	compare->GetXaxis()->SetBinLabel(3,"onlyData");
	compare->SetXTitle(xLabel.Data());
	compare->SetYTitle(yLabel.Data());

	compare->GetYaxis()->SetLabelFont(42);
	compare->GetXaxis()->SetLabelFont(42);
	compare->GetYaxis()->SetTitleFont(62);
	compare->GetXaxis()->SetTitleFont(62);

	compare->GetYaxis()->SetLabelSize(0.035);
    compare->GetYaxis()->SetTitleSize(0.04);
	compare->GetYaxis()->SetTitleOffset(1);
	compare->GetYaxis()->SetDecimals();

	compare->GetXaxis()->SetLabelSize(0.035);
    compare->GetXaxis()->SetTitleSize(0.04);
	compare->GetXaxis()->SetTitleOffset(1);

	for(Int_t iY=1; iY<nCaloCells+1; iY++){
		TH1D* temp = (TH1D*) fHistDataCell->ProjectionX("checkData",iY,iY);
		TH1D* tempMC = (TH1D*) fHistMCCell->ProjectionX("checkMC",iY,iY);
		if(temp->Integral(1,temp->GetNbinsX()) > 0 && tempMC->Integral(1,tempMC->GetNbinsX()) == 0) {
			compare->SetBinContent(3,iY,2);
			if(obj){
				if(!CheckGoodCell(obj,iY-1)) obj->cellIDsMissing.push_back(iY-1);
			}
		}
		else if(temp->Integral(1,temp->GetNbinsX()) == 0 && tempMC->Integral(1,tempMC->GetNbinsX()) > 0){
			compare->SetBinContent(1,iY,2);
			if(obj){
				if(!CheckGoodCell(obj,iY-1)) obj->cellIDsMissing.push_back(iY-1);
			}
		}
		else compare->SetBinContent(2,iY,1);

		delete temp;
		delete tempMC;
	}

	compare->GetZaxis()->SetRangeUser(0,2);
	compare->DrawCopy("COLZ");
	delete compare;
	return;
}

TH2D* CompareDeadCellsRunwise(TH1D* fHistCell, Int_t nCaloCells, TString plotRun){
	TH2D* compare = new TH2D(Form("DeadCells_Run-%s",plotRun.Data()),Form("%s",plotRun.Data()),2,0,2,nCaloCells,0,nCaloCells);
	compare->GetXaxis()->SetBinLabel(1,"only in Data");
	compare->GetXaxis()->SetBinLabel(2,"only in MC");

	compare->GetYaxis()->SetLabelFont(42);
	compare->GetXaxis()->SetLabelFont(42);
	compare->GetYaxis()->SetTitleFont(62);
	compare->GetXaxis()->SetTitleFont(62);

	compare->GetYaxis()->SetLabelSize(0.035);
    compare->GetYaxis()->SetTitleSize(0.04);
	compare->GetYaxis()->SetTitleOffset(1);
	compare->GetYaxis()->SetDecimals();

	compare->GetXaxis()->SetLabelSize(0.035);
    compare->GetXaxis()->SetTitleSize(0.04);
	compare->GetXaxis()->SetTitleOffset(1);

	for(Int_t iY=1; iY<nCaloCells+1; iY++){
		Double_t nFired = (Double_t) fHistCell->GetBinContent(iY);
		if(nFired > 1) {
			compare->SetBinContent(1,iY,1);
			compare->SetBinContent(2,iY,-2);
		}
	}
	//compare->GetZaxis()->SetRangeUser(0,2);
	return compare;
}

void CheckBadCellCandidates(TH2* fHistDataCell, TH2* fHistMCPytCell, TH2* fHistMCPhoCell, std::vector<Int_t> &vec, TString str)
{
	if(str.CompareTo("HotCellsTime1D")==0){
		std::vector<Int_t> erase;
		for(Int_t i=0; i<(Int_t)vec.size(); i++)
		{
			Int_t cellID = vec.at(i);
			TH1D* tempData = (TH1D*) fHistDataCell->ProjectionX("numberOfEntriesAboveData",cellID+1, cellID+1);
			TH1D* tempMCPyt= (TH1D*) fHistMCPytCell->ProjectionX("numberOfEntriesAboveMCPyt",cellID+1, cellID+1);
			TH1D* tempMCPho = (TH1D*) fHistMCPhoCell->ProjectionX("numberOfEntriesAboveMCPho",cellID+1, cellID+1);
			Double_t nData = tempData->Integral(tempData->GetXaxis()->FindBin(-0.02E-6),tempData->GetXaxis()->FindBin(0.02E-6));
			Double_t nMCPyt = tempMCPyt->Integral(tempMCPyt->GetXaxis()->FindBin(0.5E-6),tempMCPyt->GetXaxis()->FindBin(0.7E-6));
			Double_t nMCPho = tempMCPho->Integral(tempMCPho->GetXaxis()->FindBin(0.5E-6),tempMCPho->GetXaxis()->FindBin(0.7E-6));
			if(nMCPyt < 10 && nMCPho < 10) erase.push_back(i);

			delete tempData;
			delete tempMCPyt;
			delete tempMCPho;
		}
		for(Int_t i=0; i<(Int_t)erase.size(); i++) vec.erase(vec.begin()+(erase.at(i)-i));
	}
	return;
}


void CheckBadCellCandidatesVec(std::vector<TH2D*> &vecHistMC, std::vector<Int_t> &vec, TString str)
{
	if(str.CompareTo("HotCellsTime1D")==0){
		std::vector<Int_t> erase;
		for(Int_t i=0; i<(Int_t)vec.size(); i++)
		{
			Int_t cellID = vec.at(i);
			Bool_t boolErase = kTRUE;
			for(Int_t j=1; j<(Int_t)vecHistMC.size(); j++)
			{
				TH1D* tempMC = (TH1D*) vecHistMC.at(j)->ProjectionX("numberOfEntriesAboveMCPyt",cellID+1, cellID+1);
				Double_t nMC = tempMC->Integral(tempMC->GetXaxis()->FindBin(0.5E-6),tempMC->GetXaxis()->FindBin(0.7E-6));
				if(nMC >= 10) boolErase = kFALSE;
				delete tempMC;
			}
			if(boolErase) erase.push_back(i);
		}
		for(Int_t i=0; i<(Int_t)erase.size(); i++) vec.erase(vec.begin()+(erase.at(i)-i));
	}
	return;
}

void CollectAndPlotBadCellCandidates(fstream &str, std::vector<Int_t> &allCell, std::vector<Int_t> &vec, TString where){
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "Suspicious cells from: " << where.Data() << endl;
	cout << "Total number of " << vec.size()<< " cells" << endl;
	str << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	str << "Suspicious cells from: " << where.Data() << endl;
	str << "Total number of " << vec.size()<< " cells" << endl;
	for(Int_t iC=0; iC<(Int_t)vec.size(); iC++){
		Int_t currentCell = vec.at(iC);
		allCell.push_back(currentCell);
		cout << currentCell << ", ";
		str << currentCell << ", ";
	}
	cout << endl;
	str << endl;
	return;
}

void PlotBadCellReasons(CellQAObj* obj, std::vector<Int_t> allCells, TCanvas* canvas, TString outputDir, TString suffix, TString fClusters, TString fPlot, TString outputPlot, TString fCollisionSystem){
	Int_t atPlotting = 50;
	Int_t iCount = 0;
	do{
		atPlotting = (iCount+1)*50;
		if((Int_t)allCells.size()<=atPlotting) atPlotting = (Int_t)allCells.size();
		TH2D* fHistBadCellCand = new TH2D("BadCellReasons","BadCellCandidates",6,0.5,6.5,atPlotting-iCount*50,0+iCount*50,atPlotting);
		fHistBadCellCand->GetXaxis()->SetBinLabel(1,"Energy #mu/#sigma");
		fHistBadCellCand->GetXaxis()->SetBinLabel(2,"Time #mu/#sigma");
		fHistBadCellCand->GetXaxis()->SetBinLabel(3,"HotCells1D");
		fHistBadCellCand->GetXaxis()->SetBinLabel(4,"HotCellsTime1D");
		fHistBadCellCand->GetXaxis()->SetBinLabel(5,"HotCells2D");
		fHistBadCellCand->GetXaxis()->SetBinLabel(6,"MissingMC");
		fHistBadCellCand->GetZaxis()->SetRangeUser(0,2);
		for(Int_t iC=iCount*50; iC<atPlotting; iC++){
			fHistBadCellCand->GetYaxis()->SetBinLabel(iC+1-iCount*50,Form("%i",allCells.at(iC)));
			std::vector<Int_t>::iterator it;
			it = find (obj->cellIDsEnergy.begin(), obj->cellIDsEnergy.end(), allCells.at(iC));
			if (it != obj->cellIDsEnergy.end()) fHistBadCellCand->SetBinContent(1,iC+1-iCount*50,1);
			it = find (obj->cellIDsTime.begin(), obj->cellIDsTime.end(), allCells.at(iC));
			if (it != obj->cellIDsTime.end()) fHistBadCellCand->SetBinContent(2,iC+1-iCount*50,1);
			it = find (obj->cellIDsHotCells1D.begin(), obj->cellIDsHotCells1D.end(), allCells.at(iC));
			if (it != obj->cellIDsHotCells1D.end()) fHistBadCellCand->SetBinContent(3,iC+1-iCount*50,1);
			it = find (obj->cellIDsHotCellsTime1D.begin(), obj->cellIDsHotCellsTime1D.end(), allCells.at(iC));
			if (it != obj->cellIDsHotCellsTime1D.end()) fHistBadCellCand->SetBinContent(4,iC+1-iCount*50,1);
			it = find (obj->cellIDsHotCells2D.begin(), obj->cellIDsHotCells2D.end(), allCells.at(iC));
			if (it != obj->cellIDsHotCells2D.end()) fHistBadCellCand->SetBinContent(5,iC+1-iCount*50,1);
			it = find (obj->cellIDsMissing.begin(), obj->cellIDsMissing.end(), allCells.at(iC));
			if (it != obj->cellIDsMissing.end()) fHistBadCellCand->SetBinContent(6,iC+1-iCount*50,1);
		}
		DrawAutoGammaHisto2D(fHistBadCellCand,
							 Form("%s - %s - %s",fCollisionSystem.Data(), fPlot.Data(), fClusters.Data()),
							 "Bad Cell Candidate found by...",
							 "Cell ID",
							 "",
							 0,0,0,
							 0,0,0,
							 1,1);
		canvas->SetLogx(0); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
		canvas->SaveAs(Form("%s/Cells/ReasonsBadCell_%s_%i.%s", outputDir.Data(), outputPlot.Data(),iCount++,suffix.Data()));
		canvas->Clear();
		delete fHistBadCellCand;
	}while(atPlotting != (Int_t)allCells.size());

	return;
}


void PlotBadCellOverview(Bool_t bEnergy, Bool_t bMC, TH2* fHistCell, std::vector<Int_t> allCells, TCanvas* canvas, TString outputDir, TString suffix, TString fClusters, TString fPlot, TString outputPlot, TString fCollisionSystem){
	if(bEnergy){
		Int_t atPlotting = 0;
		Int_t iCount = 0;
		do{
			atPlotting = (iCount+1)*50;
			if((Int_t)allCells.size()<=atPlotting) atPlotting = (Int_t)allCells.size();
			TH2D* fHistBadCellCand = new TH2D("BadCellCandidates","BadCellCandidates",fHistCell->GetXaxis()->GetNbins(),fHistCell->GetXaxis()->GetXmin(),fHistCell->GetXaxis()->GetXmax(),atPlotting-iCount*50,0+iCount*50,atPlotting);
			for(Int_t iC=iCount*50; iC<atPlotting; iC++){
				fHistBadCellCand->GetYaxis()->SetBinLabel(iC+1-iCount*50,Form("%i",allCells.at(iC)));
				for(Int_t iX=1; iX<fHistCell->GetXaxis()->GetNbins()+1; iX++){
					fHistBadCellCand->SetBinContent(iX,iC+1-iCount*50,fHistCell->GetBinContent(iX,allCells.at(iC)+1));
				}
			}
			DrawAutoGammaHisto2D(fHistBadCellCand,
								 Form("%s - %s - %s",fCollisionSystem.Data(), fPlot.Data(), fClusters.Data()),
								 "Cell Energy (GeV)",
								 "Cell ID of Bad Cell Candidate",
								 "",
								 0,0,0,
								 1,0,10,
								 1,1);
			canvas->SetLogx(0); canvas->SetLogy(0); canvas->SetLogz(1); canvas->Update();
			canvas->SaveAs(Form("%s/Cells/CellEnergyVsCellID_%s_%i.%s", outputDir.Data(), outputPlot.Data(),iCount++,suffix.Data()));
			canvas->Clear();
			delete fHistBadCellCand;
		}while(atPlotting != (Int_t)allCells.size());
	}else{
		Int_t atPlotting = 0;
		Int_t iCount = 0;
		Double_t ranges[2];
		if(bMC){
			ranges[0]=0.5E-6; ranges[1]=0.7E-6;
		}else{
			ranges[0]=-0.2E-6; ranges[1]=0.2E-6;
		}
		do{
			atPlotting = (iCount+1)*50;
			if((Int_t)allCells.size()<=atPlotting) atPlotting = (Int_t)allCells.size();
			TH2D* fHistBadCellCand = new TH2D("BadCellCandidates","BadCellCandidates",fHistCell->GetXaxis()->GetNbins(),fHistCell->GetXaxis()->GetXmin(),fHistCell->GetXaxis()->GetXmax(),atPlotting-iCount*50,0+iCount*50,atPlotting);
			for(Int_t iC=iCount*50; iC<atPlotting; iC++){
				fHistBadCellCand->GetYaxis()->SetBinLabel(iC+1-iCount*50,Form("%i",allCells.at(iC)));
				for(Int_t iX=1; iX<fHistCell->GetXaxis()->GetNbins()+1; iX++){
					fHistBadCellCand->SetBinContent(iX,iC+1-iCount*50,fHistCell->GetBinContent(iX,allCells.at(iC)+1));
				}
			}
			DrawAutoGammaHisto2D(fHistBadCellCand,
								 Form("%s - %s - %s",fCollisionSystem.Data(), fPlot.Data(), fClusters.Data()),
								 "Cell Time (#mus)",
								 "Cell ID of Bad Cell Candidate",
								 "",
								 0,0,0,
								 0,0,0,
								 1,1);
			canvas->SetLogx(0); canvas->SetLogy(0); canvas->SetLogz(1); canvas->Update();
			canvas->SaveAs(Form("%s/Cells/CellTimeVsCellID_NoZoom_%s_%i.%s", outputDir.Data(), outputPlot.Data(),iCount,suffix.Data()));
			canvas->Clear();
			DrawAutoGammaHisto2D(fHistBadCellCand,
								 Form("%s - %s - %s",fCollisionSystem.Data(), fPlot.Data(), fClusters.Data()),
								 "Cell Time (#mus)",
								 "Cell ID of Bad Cell Candidate",
								 "",
								 0,0,0,
								 1,ranges[0],ranges[1],
								 1,1);
			canvas->SetLogx(0); canvas->SetLogy(0); canvas->SetLogz(1); canvas->Update();
			canvas->SaveAs(Form("%s/Cells/CellTimeVsCellID_%s_%i.%s", outputDir.Data(), outputPlot.Data(),iCount++,suffix.Data()));
			canvas->Clear();
			delete fHistBadCellCand;
		}while(atPlotting != (Int_t)allCells.size());
	}
	return;
}

void SaveBadCellCandidates(std::vector<Int_t> vec,TString str){
	if((Int_t)vec.size()>0){
		TH1D* hist = new TH1D(str.Data(),str.Data(),(Int_t)vec.size(),0,(Int_t)vec.size());
		for(Int_t ii=0; ii<(Int_t)vec.size(); ii++){
			hist->GetXaxis()->SetBinLabel(ii+1,Form("%i",vec.at(ii)));
		}
		WriteHistogram(hist);
	}
	return;
}

void PlotBadCellComparison(TH2* fHistBadCellCandData, TH2* fHistBadCellCandMCPyt, TH2* fHistBadCellCandMCPho, std::vector<Int_t> allCells, TCanvas* canvas, TString outputDir, TString suffix, TString calo, TString fPlotData, TString fPlotMCPyt, TString fPlotMCPho, TString fCollisionSystem){
    SetXRange(fHistBadCellCandData,1,fHistBadCellCandData->GetNbinsX());
    SetXRange(fHistBadCellCandMCPyt,1,fHistBadCellCandMCPyt->GetNbinsX());
    SetXRange(fHistBadCellCandMCPho,1,fHistBadCellCandMCPho->GetNbinsX());

	for(Int_t iCell=0; iCell<(Int_t)allCells.size(); iCell++)
	{
		TH1D* fHistDataProject = (TH1D*) fHistBadCellCandData->ProjectionX(Form("dataProject_%i",allCells.at(iCell)),allCells.at(iCell)+1,allCells.at(iCell)+1);
		TH1D* fHistMCPytProject = (TH1D*) fHistBadCellCandMCPyt->ProjectionX(Form("mcPytProject_%i",allCells.at(iCell)),allCells.at(iCell)+1,allCells.at(iCell)+1);
		TH1D* fHistMCPhoProject = (TH1D*) fHistBadCellCandMCPho->ProjectionX(Form("mcPhoProject_%i",allCells.at(iCell)),allCells.at(iCell)+1,allCells.at(iCell)+1);

		DrawAutoGammaCompare3H(fHistMCPytProject, fHistMCPhoProject, fHistDataProject,
                            "",
							"Cell Energy (GeV)",
							Form("CellID %i",allCells.at(iCell)),
							0,0,0,
							1,0.5,1E5,
							1,0,10,
							1,1.3,fPlotData.Data(),fPlotMCPyt.Data(),fPlotMCPho.Data());
		PutProcessLabelAndEnergyOnPlot(0.8, 0.8, 0.03, fCollisionSystem.Data(), Form("%s clusters", calo.Data()), "");
		canvas->SetLogx(0); canvas->SetLogy(1); canvas->SetLogz(0); canvas->Update();
		canvas->SaveAs(Form("%s/Cells/Detailed/Cell%i_EnergyComparison.%s", outputDir.Data(), allCells.at(iCell), suffix.Data()));
		canvas->Clear();

		delete fHistDataProject;
		delete fHistMCPytProject;
		delete fHistMCPhoProject;
	}
	return;
}

void PlotBadCellComparisonVec( std::vector<TH2D*> DataMCHists, Color_t *color, std::vector<Int_t> allCells, TCanvas* canvas, TString outputDir, TString suffix, TString calo, TString* plotDataSets, TString fCollisionSystem)
{
	std::vector<TH1D*> fVecDataMC;
	for(Int_t i=0; i<(Int_t)DataMCHists.size(); i++){
      SetXRange(DataMCHists.at(i),1,DataMCHists.at(i)->GetNbinsX());
	}

	for(Int_t iCell=0; iCell<(Int_t)allCells.size(); iCell++)
	{
		TLegend* leg1 = new TLegend( 0.6,0.81,0.97,0.95);
		leg1->SetTextSize(0.04);
		leg1->SetFillColor(0);

		for(Int_t i=0; i<(Int_t)DataMCHists.size(); i++){
			if(i==0) fVecDataMC.push_back((TH1D*) DataMCHists.at(i)->ProjectionX(Form("dataProject_%i",allCells.at(iCell)),allCells.at(iCell)+1,allCells.at(iCell)+1));
			else fVecDataMC.push_back((TH1D*) DataMCHists.at(i)->ProjectionX(Form("mcProject_%i",allCells.at(iCell)),allCells.at(iCell)+1,allCells.at(iCell)+1));
		}

		AdjustHistRange(fVecDataMC,5,5,kTRUE);
		Int_t iStart = (Int_t) fVecDataMC.size()-1;
		Int_t minB=0; Int_t maxB=0;
		for(Int_t i=iStart; i>=0; i--){
			TH1D* fHistDataMC = fVecDataMC.at(i);
			GetMinMaxBin(fHistDataMC,minB,maxB);
            SetXRange(fHistDataMC,minB,maxB+10);
            fHistDataMC->SetTitle("");
			fHistDataMC->SetXTitle("Cell Energy (GeV)");
			fHistDataMC->SetYTitle(Form("Cell ID %i",allCells.at(iCell)));

			if(i==0){
                DrawGammaSetMarker(fHistDataMC, 24, 0.5, color[i], color[i]);
				leg1->AddEntry(fHistDataMC,plotDataSets[i].Data(),"p");
			}else{
                DrawGammaSetMarker(fHistDataMC, 1, 0.5, color[i], color[i]);
				leg1->AddEntry(fHistDataMC,plotDataSets[i].Data());
			}

			fHistDataMC->GetYaxis()->SetLabelFont(42);
			fHistDataMC->GetXaxis()->SetLabelFont(42);
			fHistDataMC->GetYaxis()->SetTitleFont(62);
			fHistDataMC->GetXaxis()->SetTitleFont(62);

			fHistDataMC->GetYaxis()->SetLabelSize(0.035);
			fHistDataMC->GetYaxis()->SetTitleSize(0.04);
			fHistDataMC->GetYaxis()->SetDecimals();
			fHistDataMC->GetYaxis()->SetTitleOffset(1.3);
			fHistDataMC->GetXaxis()->SetTitleSize(0.04);
			fHistDataMC->GetXaxis()->SetLabelSize(0.035);
			fHistDataMC->GetXaxis()->SetTitleOffset(1.0);
			if(i==0) fHistDataMC->DrawCopy("x0,e,p,same");
			else if(i==iStart) fHistDataMC->DrawCopy("e,hist");
			else fHistDataMC->DrawCopy("e,hist,same");
		}
		leg1->Draw("same");

		PutProcessLabelAndEnergyOnPlot(0.7, 0.8, 0.03, fCollisionSystem.Data(), Form("%s clusters", calo.Data()), "");
		SaveCanvas(canvas,Form("%s/Cells/Detailed/Cell%i_EnergyComparison.%s", outputDir.Data(), allCells.at(iCell), suffix.Data()),0,1,0);
		delete leg1;

		for(Int_t i=0; i<(Int_t)fVecDataMC.size(); i++) delete fVecDataMC.at(i);
		fVecDataMC.clear();
	}

	return;
}

void CheckHotAndColdCellsEFracRunwise(fstream &fLogRunwiseHotCells, fstream &fLogRunwiseColdCells, TH1D* tempClusEFraccCellBefore, TProfile* BadCells, Int_t nCaloCells, Bool_t simpleOutput){
	if(tempClusEFraccCellBefore->Integral()<1E4){
		fLogRunwiseHotCells << "NotEnoughCellsFired-" << tempClusEFraccCellBefore->Integral() << "-ToObtainNoisyChannels!" << endl;
		return;
	}
	Bool_t doColdCells = kTRUE;
	if(tempClusEFraccCellBefore->Integral()<1E5){
		fLogRunwiseColdCells << "NotEnoughCellsFired-" << tempClusEFraccCellBefore->Integral() << "-ToObtainColdChannels!" << endl;
		doColdCells = kFALSE;
	}

	Int_t noisyCells = 0;
	Int_t coldCells = 0;
	Double_t nCellsBuffer = 100;
	Double_t mean = 0;
	Double_t hmean = 0;
	Int_t cell=1;
	Int_t stopCell = 1;
	do{
		Double_t nCurrentEFrac = (Double_t) tempClusEFraccCellBefore->GetBinContent(cell++);
		if(nCurrentEFrac>0){
			hmean+=(nCurrentEFrac/nCellsBuffer);
			stopCell++;
		}
	}while(stopCell<nCellsBuffer);

	queue <Double_t> lastNFired;
	for(Int_t iCell = 1; iCell<=stopCell; iCell++){
		Double_t nCurrentEFrac = (Double_t) tempClusEFraccCellBefore->GetBinContent(iCell);
		if((nCurrentEFrac>2*hmean && nCurrentEFrac>10)||(nCurrentEFrac==0)) stopCell++;
		else {
			mean+=nCurrentEFrac/nCellsBuffer;
			lastNFired.push(nCurrentEFrac);
		}
	}

	for(Int_t iCell = 1; iCell<nCaloCells+1; iCell++){
		if(BadCells && BadCells->GetBinContent(iCell)>0) continue;
		Double_t nCurrentEFrac = (Double_t) tempClusEFraccCellBefore->GetBinContent(iCell);
		if( (nCurrentEFrac>2*mean && nCurrentEFrac>80) || (nCurrentEFrac>3*mean && nCurrentEFrac>20) || (nCurrentEFrac>4*mean && nCurrentEFrac>8) || (nCurrentEFrac>5*mean && nCurrentEFrac>5)){
			if(simpleOutput) fLogRunwiseHotCells << iCell-1 << endl;
			else fLogRunwiseHotCells << iCell-1 << "-CellEFraction:" << nCurrentEFrac << "-CurrentMean:" << mean << endl;
			//fLogRunwiseBadCells << iCell << endl;
			//tempClusIncludedCellBefore->SetBinContent(iCell,0);
			noisyCells++;
		}
		else if( doColdCells && ((nCurrentEFrac<mean/3 && mean>=80) || (nCurrentEFrac<mean/5 && mean>=40 && mean<80) || (nCurrentEFrac<mean/8 && mean>=10 && mean<40) || (nCurrentEFrac<mean/10 && mean<10)) ){
			if(simpleOutput) fLogRunwiseColdCells << iCell-1 << endl;
			else fLogRunwiseColdCells << iCell-1 << "-CellEFraction:" << nCurrentEFrac << "-CurrentMean:" << mean << endl;
			coldCells++;
		}
		else if(iCell>stopCell){
            if(!doColdCells && nCurrentEFrac==0) continue;
			mean -= lastNFired.front()/nCellsBuffer;
			mean += nCurrentEFrac/nCellsBuffer;
			lastNFired.pop();
			lastNFired.push(nCurrentEFrac);
		}
	}
	if(noisyCells == 0) fLogRunwiseHotCells << "NoNoisyCellsFound!" << endl;
	if(coldCells == 0) fLogRunwiseColdCells << "NoColdCellsFound!" << endl;
	return;
}

template<class ForwardIt>
void selection_sort(ForwardIt begin, ForwardIt end)
{
	for (ForwardIt i = begin; i != end; ++i)
		std::iter_swap(i, std::min_element(i, end));
}

Float_t GetMedianTH2(TH2* hist)
{
	std::vector<Float_t> vectorMedian;
	for(Int_t x=1; x<=hist->GetXaxis()->GetNbins(); x++)
	{
		for(Int_t y=1; y<=hist->GetYaxis()->GetNbins(); y++)
		{
			if(!hist->GetBinContent(x,y)==0) vectorMedian.push_back(hist->GetBinContent(x,y));
		}
	}
	selection_sort(vectorMedian.begin(), vectorMedian.end());
	Int_t median = (Int_t) vectorMedian.size()/2;
	return vectorMedian.at(median);
}

Float_t GetMeanTH2(TH2* hist)
{
	Float_t Mean = 0;
	Float_t nZeros = 0;
	for(Int_t x=1; x<=hist->GetXaxis()->GetNbins(); x++)
	{
		for(Int_t y=1; y<=hist->GetYaxis()->GetNbins(); y++)
		{
			if(hist->GetBinContent(x,y)==0) nZeros++;
			else Mean+=hist->GetBinContent(x,y);
		}
	}
	Mean/=(((Float_t) hist->GetXaxis()->GetNbins()) * ((Float_t) hist->GetYaxis()->GetNbins())) - nZeros;
	return Mean;
}

Bool_t CheckGoodCell(CellQAObj* obj, Int_t cellID){
	std::vector<Int_t>::iterator it;
	it = find (obj->goodCells.begin(), obj->goodCells.end(), cellID);
	if (it == obj->goodCells.end()) return false;
	return true;
}
void FillGoodCells(CellQAObj* cellQAData, Int_t nGoodCells, Int_t* goodCells){
	for(Int_t i=0; i<nGoodCells; i++) cellQAData->goodCells.push_back(goodCells[i]);
	return;
}

Double_t GetHistogramIntegral(TH1D* hist, Float_t lowX, Float_t highX)
{
	Double_t Integral = 0;
	Double_t bgerror = 0;
	TAxis *axis = hist->GetXaxis();
	int bmin = axis->FindBin(lowX); //in your case xmin=-1.5
	int bmax = axis->FindBin(highX); //in your case xmax=0.8
	Integral = hist->IntegralAndError(bmin,bmax,bgerror,"");
	return Integral;
}

Double_t GetHistogramIntegralError(TH1D* hist, Float_t lowX, Float_t highX)
{
	Double_t Integral = 0;
	Double_t bgerror = 0;
	TAxis *axis = hist->GetXaxis();
	int bmin = axis->FindBin(lowX); //in your case xmin=-1.5
	int bmax = axis->FindBin(highX); //in your case xmax=0.8
	Integral = hist->IntegralAndError(bmin,bmax,bgerror,"");
	return bgerror;
}

void setQAEnergy(CellQAObj* obj, Double_t min, Double_t max, Double_t minSigma, Double_t maxSigma){
	obj->EnergyMean[0]=min;obj->EnergyMean[1]=max;
	obj->EnergySigma[0]=minSigma;obj->EnergySigma[1]=maxSigma;
	return;
}
void setQATime(CellQAObj* obj, Double_t min, Double_t max, Double_t minSigma, Double_t maxSigma){
	obj->TimeMean[0]=min;obj->TimeMean[1]=max;
	obj->TimeSigma[0]=minSigma;obj->TimeSigma[1]=maxSigma;
	return;
}
void setQAHotCells1D(CellQAObj* obj, Double_t min, Double_t max, Double_t minT, Double_t maxT){
	obj->HotCells1D[0]=min;obj->HotCells1D[1]=max;
	obj->HotCellsTime1D[0]=minT;obj->HotCellsTime1D[1]=maxT;
	return;
}
void setQAHotCells2D(CellQAObj* obj, Int_t size, Double_t* min, Double_t* max){
	if(size!=9){
		cout << "Wrong size of array given to setQAHotCells2D, returning..." << endl;
		return;
	}
	for(Int_t iFill=0; iFill<size; iFill++)
	{
		obj->HotCells2D[iFill][0]=min[iFill];
		obj->HotCells2D[iFill][1]=max[iFill];
	}
	return;
}

void CalculateFWHM(TF1 * fFunc, Double_t startMass, Double_t endMass)
{
   // Default function
   TF1* fFunc_def;
   fFunc_def = new TF1("fFunc_def","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",startMass,endMass);
   fFunc_def->SetParameter(0,fFunc->GetParameter(0));
   fFunc_def->SetParameter(1,fFunc->GetParameter(1));
   fFunc_def->SetParameter(2,fFunc->GetParameter(2));
   fFunc_def->SetParameter(3,fFunc->GetParameter(3));



   //FWHM
   fFWHMFunc = fFunc->GetX(fFunc_def->GetParameter(0)*0.5,fFunc_def->GetParameter(1), endMass) - fFunc_def->GetX(fFunc_def->GetParameter(0)*0.5,startMass,fFunc_def->GetParameter(1));

   //FWHM error +
   TF1* fFunc_plus;
   // fFunc_plus = fFunc;
   fFunc_plus = new TF1("fFunc_plus","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",startMass,endMass);

   fFunc_plus->SetParameter(0,fFunc->GetParameter(0) + fFunc->GetParError(0));
   fFunc_plus->SetParameter(1,fFunc->GetParameter(1) + fFunc->GetParError(1));
   fFunc_plus->SetParameter(2,fFunc->GetParameter(2) + fFunc->GetParError(2));
   fFunc_plus->SetParameter(3,fFunc->GetParameter(3) + fFunc->GetParError(3));
   Double_t FWHM_plus = fFunc_plus->GetX(fFunc_plus->GetParameter(0)*0.5,fFunc_plus->GetParameter(1), endMass) - fFunc_plus->GetX(fFunc_plus->GetParameter(0)*0.5,startMass,fFunc_plus->GetParameter(1));

   //FWHM error -
   TF1* fFunc_minus;
   // fFunc_minus = fFunc;
   fFunc_minus = new TF1("fFunc_minus","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",startMass,endMass);
   fFunc_minus->SetParameter(0,fFunc->GetParameter(0) - fFunc->GetParError(0));
   fFunc_minus->SetParameter(1,fFunc->GetParameter(1) - fFunc->GetParError(1));
   fFunc_minus->SetParameter(2,fFunc->GetParameter(2) - fFunc->GetParError(2));
   fFunc_minus->SetParameter(3,fFunc->GetParameter(3) - fFunc->GetParError(3));
   Double_t FWHM_minus =  fFunc_minus->GetX(fFunc_minus->GetParameter(0)*0.5,fFunc_minus->GetParameter(1), endMass) -fFunc_minus->GetX(fFunc_minus->GetParameter(0)*0.5,startMass,fFunc_minus->GetParameter(1));

   Double_t Error1 = TMath::Abs(fFWHMFunc-FWHM_plus);
   Double_t Error2 = TMath::Abs(fFWHMFunc-FWHM_minus);

   if(Error1>=Error2) fFWHMFuncError = Error1;
   if(Error1<Error2) fFWHMFuncError = Error2;

   delete fFunc_plus;
   delete fFunc_minus;
   delete fFunc_def;

   return;
}


#endif // QA_H
