#include <Riostream.h>
#include <fstream>
#include "TMath.h"
#include <stdio.h>
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

struct SysErrorConversion {
   Double_t value;
   Double_t error;
};


void GammaCutStudiesV2(TString cutFile = "CombineCuts.dat",TString energy="",TString cutVariationName ="",TString suffix = "eps", Int_t mode = 9){
	
	gStyle->SetOptStat(0);   //gStyle->SetOptFit(1111);
	
	TString outputDir = Form("CutStudies/%s/%s",energy.Data(),cutVariationName.Data());
	TString outputFileDir = Form("CutStudies/%s",energy.Data());
	gSystem->Exec("mkdir -p "+outputDir);

	TString collisionSystem= ReturnFullCollisionsSystem(energy);   
	if (collisionSystem.CompareTo("") == 0){
		cout << "No correct collision system specification, has been given" << endl;
		return;     
	}

	cout<<endl;cout<<endl;cout<<endl;
	cout<<"Processing '"<<cutVariationName<<"' Study for "<<collisionSystem<<endl;
	cout<<endl;
	TString eventCutSelection[50];
	TString gammaCutSelection[50];
	TString clusterCutSelection[50];
	TString electronCutSelection[50];
	TString mesonCutSelection[50];
	TString cutSelection[50];
	
	ifstream in(cutFile);
	cout<<"=========================="<<endl;
	cout<<"Available Cuts:"<<endl;
	string currentCutNumber;
	Int_t number = 0;

	while(getline(in, currentCutNumber)){
		cutSelection[number] = currentCutNumber;
		cout<<"---> "<<currentCutNumber<<endl;
		
		if (mode == 9){
			ReturnSeparatedCutNumber(cutSelection[number], gammaCutSelection[number], electronCutSelection[number],mesonCutSelection[number]);
			eventCutSelection[number] = gammaCutSelection[number](0,7);
			gammaCutSelection[number] = gammaCutSelection[number](7,gammaCutSelection[number].Length()-7);
			cout << eventCutSelection[number].Data() << "\t" << gammaCutSelection[number].Data() << endl;
		} else {
			ReturnSeparatedCutNumberAdvanced(cutSelection[number],eventCutSelection[number], gammaCutSelection[number], clusterCutSelection[number], electronCutSelection[number], mesonCutSelection[number], mode);
		}	
	}

	cout<<"=========================="<<endl;

	TF1 *One = new TF1("One","1",0,25);
	One->SetLineWidth(1.2);
	One->SetLineColor(1);

	TString currentCorrectedFile = "";
	TH1D **InclusiveGammaSpectrum = new TH1D*[number];
	TH1D **InclusiveGammaSpectrumRatio = new TH1D*[number];
	TH1D **Pi0Spectrum = new TH1D*[number];
	TH1D **Pi0SpectrumRatio = new TH1D*[number];
	TH1D **Pi0SpectrumFit = new TH1D*[number];
	TH1D **Pi0SpectrumFitRatio = new TH1D*[number];
	TH1D **InclusiveGammaToPi0Ratio = new TH1D*[number];
	TH1D **InclusiveGammaToPi0RatioRatio = new TH1D*[number];
	TH1D **InclusiveGammaToPi0RatioFit = new TH1D*[number];
	TH1D **InclusiveGammaToPi0RatioFitRatio = new TH1D*[number];
	TH1D **DirectPhotonDoubleRatio = new TH1D*[number];
	TH1D **DirectPhotonDoubleRatioRatio = new TH1D*[number];
	TH1D **DirectPhotonDoubleRatioFit = new TH1D*[number];
	TH1D **DirectPhotonDoubleRatioFitRatio = new TH1D*[number];
	// TH1D **DirectPhotonSpectrum = new TH1D*[number];
	// TH1D **DirectPhotonSpectrumRatio = new TH1D*[number];

	TH1D **Purity = new TH1D*[number];
	TH1D **GammaEff = new TH1D*[number];
	TH1D **RawGamma = new TH1D*[number];
	TH1D **PurityRatio = new TH1D*[number];
	TH1D **GammaEffRatio = new TH1D*[number];
	TH1D **RawGammaRatio = new TH1D*[number];

	TFile **currentFinalFile = new TFile*[number];
	TFile **currentCorrectionFile = new TFile*[number];

	TCanvas *GammaSpectrumCanvas = GetAndSetCanvas("GammaSpectra");GammaSpectrumCanvas->SetLogy();
	TCanvas *Pi0SpectrumCanvas = GetAndSetCanvas("Pi0Spectra");Pi0SpectrumCanvas->SetLogy();
	TCanvas *Pi0SpectrumFitCanvas = GetAndSetCanvas("Pi0FitSpectra");Pi0SpectrumFitCanvas->SetLogy();
	TCanvas *InclusiveGammaToPi0RatioCanvas = GetAndSetCanvas("InclusiveRatios");
	TCanvas *InclusiveGammaToPi0RatioFitCanvas = GetAndSetCanvas("InclusiveRatiosFit");
	TCanvas *DirectPhotonDoubleRatioCanvas= GetAndSetCanvas("DoubleRatios");
	TCanvas *DirectPhotonDoubleRatioFitCanvas= GetAndSetCanvas("DoubleRatiosFit");
	//TCanvas *DirectPhotonSpectrumCanvas= GetAndSetCanvas("DirectPhotonSpectrum");DirectPhotonSpectrumCanvas->SetLogy();

	TCanvas *GammaSpectrumRatioCanvas = GetAndSetCanvas("GammaSpectraRatio");
	TCanvas *Pi0SpectrumRatioCanvas = GetAndSetCanvas("Pi0SpectraRatio");
	TCanvas *Pi0SpectrumFitRatioCanvas = GetAndSetCanvas("Pi0FitSpectraRatio");
	TCanvas *InclusiveGammaToPi0RatioRatioCanvas = GetAndSetCanvas("InclusiveRatiosRatio");
	TCanvas *InclusiveGammaToPi0RatioFitRatioCanvas = GetAndSetCanvas("InclusiveRatiosFitRatio");
	TCanvas *DirectPhotonDoubleRatioRatioCanvas= GetAndSetCanvas("DoubleRatiosRatio");
	TCanvas *DirectPhotonDoubleRatioFitRatioCanvas= GetAndSetCanvas("DoubleRatiosFitRatio");
	//TCanvas *DirectPhotonSpectrumRatioCanvas= GetAndSetCanvas("DirectPhotonSpectrumRatio");
	
	TCanvas *RawGammaCanvas = GetAndSetCanvas("RawGamma");RawGammaCanvas->SetLogy();
	TCanvas *PurityCanvas = GetAndSetCanvas("Purity");
	TCanvas *GammaEffCanvas = GetAndSetCanvas("GammaEff");

	TCanvas *RawGammaRatioCanvas = GetAndSetCanvas("RawGammaRatio");
	TCanvas *PurityRatioCanvas = GetAndSetCanvas("PurityRatio");
	TCanvas *GammaEffRatioCanvas = GetAndSetCanvas("GammaEffRatio");

	TLegend *GammaSpectrumLegend = GetAndSetLegend(0.3,0.7,number);
	TLegend *GammaSpectrumRatioLegend = GetAndSetLegend(0.15,0.75,number);
	TLegend *Pi0SpectrumLegend = GetAndSetLegend(0.3,0.7,number);
	TLegend *Pi0SpectrumRatioLegend = GetAndSetLegend(0.15,0.75,number);
	TLegend *Pi0SpectrumFitLegend = GetAndSetLegend(0.3,0.7,number);
	TLegend *Pi0SpectrumFitRatioLegend = GetAndSetLegend(0.15,0.75,number);
	TLegend *PurityLegend = GetAndSetLegend(0.2,0.3,number);
	TLegend *PurityRatioLegend = GetAndSetLegend(0.15,0.15,number);
	TLegend *GammaEffLegend = GetAndSetLegend(0.3,0.3,number);
	TLegend *GammaEffRatioLegend = GetAndSetLegend(0.15,0.75,number);
	TLegend *InclusiveGammaToPi0RatioLegend = GetAndSetLegend(0.3,0.7,number);
	TLegend *InclusiveGammaToPi0RatioRatioLegend = GetAndSetLegend(0.15,0.75,number);
	TLegend *InclusiveGammaToPi0RatioFitLegend = GetAndSetLegend(0.3,0.7,number);
	TLegend *InclusiveGammaToPi0RatioFitRatioLegend = GetAndSetLegend(0.15,0.75,number);
	TLegend *DirectPhotonDoubleRatioLegend = GetAndSetLegend(0.3,0.7,number);
	TLegend *DirectPhotonDoubleRatioRatioLegend = GetAndSetLegend(0.15,0.75,number);
	TLegend *DirectPhotonDoubleRatioFitLegend = GetAndSetLegend(0.3,0.7,number);
	TLegend *DirectPhotonDoubleRatioFitRatioLegend = GetAndSetLegend(0.15,0.75,number);
	// TLegend *DirectPhotonSpectrumLegend = GetAndSetLegend(0.3,0.7,number);
	// TLegend *DirectPhotonSpectrumRatioLegend = GetAndSetLegend(0.15,0.15,number);
	TLegend *RawGammaLegend = GetAndSetLegend(0.2,0.3,number);
	TLegend *RawGammaRatioLegend = GetAndSetLegend(0.15,0.15,number);


	Color_t color[10] = {kBlack,kRed,kBlue,kGreen,kMagenta,kYellow,kCyan,kGreen-2,kAzure,kGray};
	Int_t markerType = 24;
	TString cutStringsName[50];
	
	
	for(Int_t i = 0; i<number; i++){
	
		TString fEventCutSelection;
		TString fGammaCutSelection;
		TString fElectronCutSelection;
		TString fMesonCutSelection;    
		TString fClusterCutSelection;
		if (mode == 9){
			ReturnSeparatedCutNumber(cutSelection[i].Data(), fGammaCutSelection, fElectronCutSelection,fMesonCutSelection);
			fEventCutSelection = fGammaCutSelection(0,7);
			fGammaCutSelection = fGammaCutSelection(7,fGammaCutSelection.Length()-7);
			cout << fEventCutSelection.Data() << "\t" << fGammaCutSelection.Data() << endl;
		} else {
			ReturnSeparatedCutNumberAdvanced(cutSelection[i].Data(),fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);
		}	

		if (cutVariationName.Contains("SpecialTrigg")){
			TString fTrigger = fEventCutSelection(GetEventSelectSpecialTriggerCutPosition(),1);
			cutStringsName[i] = AnalyseSpecialTriggerCut(fTrigger.Atoi());      
		} else if (cutVariationName.Contains("V0Reader")){
			TString fV0Reader = fGammaCutSelection(GetPhotonV0FinderCutPosition(fGammaCutSelection),1);
			cutStringsName[i] = AnalyseV0ReaderCut(fV0Reader.Atoi());      
		} else if (cutVariationName.Contains("Eta")){
			TString fEtaCut = fGammaCutSelection(GetPhotonEtaCutPosition(fGammaCutSelection),1);	
			cout << fGammaCutSelection.Data() << "\t"<<fEtaCut.Data() << endl;
			cutStringsName[i] = AnalyseEtaCut(fEtaCut.Atoi());
		} else if (cutVariationName.Contains("RCut")){
			TString fRCut = fGammaCutSelection(GetPhotonMinRCutPosition(fGammaCutSelection),1);
			cutStringsName[i] = AnalyseRCut(fRCut.Atoi());
		} else if (cutVariationName.Contains("SinglePt")){
			TString fSinglePtCut = fGammaCutSelection(GetPhotonSinglePtCutPosition(fGammaCutSelection),1);
			cutStringsName[i] = AnalyseSinglePtCut(fSinglePtCut.Atoi());
		} else if (cutVariationName.Contains("Cluster")){
			TString fClusterCut = fGammaCutSelection(GetPhotonClsTPCCutPosition(fGammaCutSelection),1);
			cutStringsName[i] = AnalyseTPCClusterCut(fClusterCut.Atoi());
			cout << i << "\t" << fClusterCut.Data() << "\t" << cutStringsName[i].Data()<< endl;
		} else if (cutVariationName.Contains("dEdxE")){	 
			TString fdEdxCut = fGammaCutSelection(GetPhotonEDedxSigmaCutPosition(fGammaCutSelection),1);
			cutStringsName[i] = AnalyseTPCdEdxCutElectronLine(fdEdxCut.Atoi());
		} else if (cutVariationName.Contains("dEdxPi")){    
			TString fdEdxCut = fGammaCutSelection(GetPhotonPiDedxSigmaCutPosition(fGammaCutSelection),3);
			cutStringsName[i] = AnalyseTPCdEdxCutPionLine(fdEdxCut.Data());      
		} else if (cutVariationName.Contains("Qt")){
			TString fQtCut = fGammaCutSelection(GetPhotonQtMaxCutPosition(fGammaCutSelection),1);
			cutStringsName[i] = AnalyseQtMaxCut(fQtCut.Atoi());
		} else if (cutVariationName.Contains("Chi2")){
			TString fChi2Cut = fGammaCutSelection(GetPhotonChi2GammaCutPosition(fGammaCutSelection),1);
			TString fPsiPairCut = fGammaCutSelection(GetPhotonPsiPairCutPosition(fGammaCutSelection),1);
			cutStringsName[i] = AnalyseChi2GammaCut(fChi2Cut.Atoi(),fPsiPairCut.Atoi());
		} else if (cutVariationName.Contains("PsiPairAndR")){
			TString fPsiPairCut = fGammaCutSelection(GetPhotonPsiPairCutPosition(fGammaCutSelection),1);
			TString fRCut = fGammaCutSelection(GetPhotonMinRCutPosition(fGammaCutSelection),1);
			cutStringsName[i] = AnalysePsiPairAndR(fPsiPairCut.Atoi(),fRCut.Atoi());      
		} else if (cutVariationName.Contains("PsiPair")){
			TString fPsiPairCut = fGammaCutSelection(GetPhotonPsiPairCutPosition(fGammaCutSelection),1);
			TString fChi2Cut = fGammaCutSelection(GetPhotonChi2GammaCutPosition(fGammaCutSelection),1);
			cutStringsName[i] = AnalysePsiPair(fPsiPairCut.Atoi(),fChi2Cut.Atoi());   
		} else if (cutVariationName.Contains("DCAZPhoton")){   
			TString fDCAZCut = fGammaCutSelection(GetPhotonDcaZPrimVtxCutPosition(fGammaCutSelection),1);
			cutStringsName[i] = AnalyseDCAZPhotonCut(fDCAZCut.Atoi());
		} else if (cutVariationName.Contains("CosPoint")){
			TString fCosPoint = fGammaCutSelection(GetPhotonCosinePointingAngleCutPosition(fGammaCutSelection),1);
			cutStringsName[i] = AnalyseCosPointCut(fCosPoint.Atoi());      	 	  	
		} else if (cutVariationName.Contains("PhotonQuality")){
			TString fPhotonQuality = fGammaCutSelection(GetPhotonSharedElectronCutPosition(fGammaCutSelection),1);
			cutStringsName[i] = AnalysePhotonQuality(fPhotonQuality.Atoi());      	 	  
		} else if (cutVariationName.Contains("BG")){
			TString fBGCut = fMesonCutSelection(GetMesonBGSchemeCutPosition(),3);
			cutStringsName[i] = AnalyseBackgroundScheme(fBGCut.Data());   
		} else if (cutVariationName.Contains("Rapidity")){
			TString fRapidityCut = fMesonCutSelection(GetMesonRapidityCutPosition(),1);
			cutStringsName[i] = AnalyseRapidityMesonCut(fRapidityCut.Atoi());      
		} else if (cutVariationName.Contains("Alpha")){
			TString fAlphaCut = fMesonCutSelection(GetMesonAlphaCutPosition(),1);
			cutStringsName[i] = AnalyseAlphaMesonCut(fAlphaCut.Atoi());
		} else if (cutVariationName.Contains("Cent")){
			cutStringsName[i] = GetCentralityString(fGammaCutSelection.Data());
		} else if (cutVariationName.Contains("DiffRapWindow")){
			TString fRapidityCut = fMesonCutSelection(GetMesonRapidityCutPosition(),1);
			cutStringsName[i] = AnalyseRapidityMesonCutpPb(fRapidityCut.Atoi());      
		} else if (cutVariationName.Contains("MCSmearing")){
			TString fMCSmearing = fMesonCutSelection(GetMesonUseMCPSmearingCutPosition(),1);
			cutStringsName[i] = AnalyseMCSmearingCut(fMCSmearing.Atoi());      
		} else if (cutVariationName.Contains("IntRange")){
			if (i==0) cutStringsName[i] = "standard #pi^{0} integration range";
			if (i==1) cutStringsName[i] = "wide #pi^{0} integration range";
			if (i==2) cutStringsName[i] = "narrow #pi^{0} integration range";
		} else {
			cutStringsName[i] = cutSelection[i].Data();
		}

		
		currentCorrectedFile = Form("%s/%s/Gamma_Pi0_data_GammaConvV1_InclusiveRatio_0-100.root",cutSelection[i].Data(),energy.Data());
		currentFinalFile[i] = new TFile(currentCorrectedFile);
	
		InclusiveGammaSpectrum[i] = (TH1D*) currentFinalFile[i]->Get("histoGammaSpecCorrPurity");
		InclusiveGammaSpectrum[i]->SetTitle("");
		DrawGammaSetMarker(InclusiveGammaSpectrum[i], markerType, 2.0, color[i], color[i]);
		InclusiveGammaSpectrumRatio[i] = (TH1D*) InclusiveGammaSpectrum[i]->Clone(Form("InclusiveGammaSpectrumRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
		InclusiveGammaSpectrumRatio[i]->Divide(InclusiveGammaSpectrumRatio[i],InclusiveGammaSpectrum[0],1,1,"b");
		SetHistogramm(InclusiveGammaSpectrumRatio[i],"p_{T} (GeV/c)","Ratios of #gamma Spectra",0.0,2);
		
		InclusiveGammaToPi0Ratio[i] = (TH1D*) currentFinalFile[i]->Get("IncRatioPurity_trueEff");
		if(i == 1 && cutVariationName.Contains("IntRange")) InclusiveGammaToPi0Ratio[i] = (TH1D*) currentFinalFile[i]->Get("IncRatioPurity_trueEffWide");
		if(i == 2 && cutVariationName.Contains("IntRange")) InclusiveGammaToPi0Ratio[i] = (TH1D*) currentFinalFile[i]->Get("IncRatioPurity_trueEffNarrow");
		InclusiveGammaToPi0Ratio[i]->SetTitle("");
		DrawGammaSetMarker(InclusiveGammaToPi0Ratio[i], markerType, 2.0, color[i], color[i]);
		InclusiveGammaToPi0RatioRatio[i] = (TH1D*) InclusiveGammaToPi0Ratio[i]->Clone(Form("InclusiveGammaToPi0RatioRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
		InclusiveGammaToPi0RatioRatio[i]->Divide(InclusiveGammaToPi0RatioRatio[i],InclusiveGammaToPi0Ratio[0],1,1,"b");
		SetHistogramm(InclusiveGammaToPi0RatioRatio[i],"p_{T} (GeV/c)","Ratios of #gamma/#pi^{0} Ratios",0,2);
		
		Pi0Spectrum[i] = (TH1D*) currentFinalFile[i]->Get("CorrectedYieldTrueEff");
		if(i == 1 && cutVariationName.Contains("IntRange")) Pi0Spectrum[i] = (TH1D*) currentFinalFile[i]->Get("CorrectedYieldTrueEffWide");
		if(i == 2 && cutVariationName.Contains("IntRange")) Pi0Spectrum[i] = (TH1D*) currentFinalFile[i]->Get("CorrectedYieldTrueEffNarrow");
		Pi0Spectrum[i]->SetTitle("");
		DrawGammaSetMarker(Pi0Spectrum[i], markerType, 2.0, color[i], color[i]);
		Pi0SpectrumRatio[i] = (TH1D*) Pi0Spectrum[i]->Clone(Form("Pi0SpectrumRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
		Pi0SpectrumRatio[i]->Divide(Pi0SpectrumRatio[i],Pi0Spectrum[0],1,1,"b");
		SetHistogramm(Pi0SpectrumRatio[i],"p_{T} (GeV/c)","Ratios of #pi^{0} Spectra",0.0,2.0);

		Pi0SpectrumFit[i] = (TH1D*) currentFinalFile[i]->Get("CorrectedYieldTrueEffPi0Fit");
		if(i == 1 && cutVariationName.Contains("IntRange")) Pi0SpectrumFit[i] = (TH1D*) currentFinalFile[i]->Get("CorrectedYieldTrueEffPi0FitWide");
		if(i == 2 && cutVariationName.Contains("IntRange")) Pi0SpectrumFit[i] = (TH1D*) currentFinalFile[i]->Get("CorrectedYieldTrueEffPi0FitNarrow");
		Pi0SpectrumFit[i]->SetTitle("");
		DrawGammaSetMarker(Pi0SpectrumFit[i], markerType, 2.0, color[i], color[i]);
		Pi0SpectrumFitRatio[i] = (TH1D*) Pi0SpectrumFit[i]->Clone(Form("Pi0SpectrumFitRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
		Pi0SpectrumFitRatio[i]->Divide(Pi0SpectrumFitRatio[i],Pi0SpectrumFit[0],1,1,"b");
		SetHistogramm(Pi0SpectrumFitRatio[i],"p_{T} (GeV/c)","Ratios of #pi^{0} Spectra",0.0,2.0);
		
		InclusiveGammaToPi0Ratio[i] = (TH1D*) currentFinalFile[i]->Get("IncRatioPurity_trueEff");
		if(i == 1 && cutVariationName.Contains("IntRange")) InclusiveGammaToPi0Ratio[i] = (TH1D*) currentFinalFile[i]->Get("IncRatioPurity_trueEffWide");
		if(i == 2 && cutVariationName.Contains("IntRange")) InclusiveGammaToPi0Ratio[i] = (TH1D*) currentFinalFile[i]->Get("IncRatioPurity_trueEffNarrow");
		InclusiveGammaToPi0Ratio[i]->SetTitle("");
		DrawGammaSetMarker(InclusiveGammaToPi0Ratio[i], markerType, 2.0, color[i], color[i]);
		InclusiveGammaToPi0RatioRatio[i] = (TH1D*) InclusiveGammaToPi0Ratio[i]->Clone(Form("InclusiveGammaToPi0RatioRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
		InclusiveGammaToPi0RatioRatio[i]->Divide(InclusiveGammaToPi0RatioRatio[i],InclusiveGammaToPi0Ratio[0],1,1,"b");
		SetHistogramm(InclusiveGammaToPi0RatioRatio[i],"p_{T} (GeV/c)","Ratios of #gamma/#pi^{0} Ratios",0,2);

		InclusiveGammaToPi0RatioFit[i] = (TH1D*) currentFinalFile[i]->Get("histoIncRatioFitPurity");
		if(i == 1 && cutVariationName.Contains("IntRange")) InclusiveGammaToPi0RatioFit[i] = (TH1D*) currentFinalFile[i]->Get("histoIncRatioFitPurityWide");
		if(i == 2 && cutVariationName.Contains("IntRange")) InclusiveGammaToPi0RatioFit[i] = (TH1D*) currentFinalFile[i]->Get("histoIncRatioFitPurityNarrow");
		if(i == 1 && cutVariationName.Contains("Fit")) InclusiveGammaToPi0RatioFit[i] = (TH1D*) currentFinalFile[i]->Get("histoIncRatioLowFitPurity");
		if(i == 2 && cutVariationName.Contains("Fit")) InclusiveGammaToPi0RatioFit[i] = (TH1D*) currentFinalFile[i]->Get("histoIncRatioHighFitPurity");

		InclusiveGammaToPi0RatioFit[i]->SetTitle("");
		DrawGammaSetMarker(InclusiveGammaToPi0RatioFit[i], markerType, 2.0, color[i], color[i]);
		InclusiveGammaToPi0RatioFitRatio[i] = (TH1D*) InclusiveGammaToPi0RatioFit[i]->Clone(Form("InclusiveGammaToPi0RatioFitRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
		InclusiveGammaToPi0RatioFitRatio[i]->Divide(InclusiveGammaToPi0RatioFitRatio[i],InclusiveGammaToPi0RatioFit[0],1,1,"b");
		SetHistogramm(InclusiveGammaToPi0RatioFitRatio[i],"p_{T} (GeV/c)","Ratios of #gamma/#pi^{0}_{Fit} Ratios",0,2);

		DirectPhotonDoubleRatio[i] = (TH1D*) currentFinalFile[i]->Get("DoubleRatioConversionTrueEffPurity");
		if(i == 1 && cutVariationName.Contains("IntRange")) DirectPhotonDoubleRatio[i] = (TH1D*) currentFinalFile[i]->Get("DoubleRatioConversionTrueEffPurityWide");
		if(i == 2 && cutVariationName.Contains("IntRange")) DirectPhotonDoubleRatio[i] = (TH1D*) currentFinalFile[i]->Get("DoubleRatioConversionTrueEffPurityNarrow");
		if(i == 1 && cutVariationName.CompareTo("CocktailEta") == 0) DirectPhotonDoubleRatio[i] = (TH1D*) currentFinalFile[i]->Get("DoubleRatioConversionTrueEffPurityHigh");
		if(i == 2 && cutVariationName.CompareTo("CocktailEta") == 0) DirectPhotonDoubleRatio[i] = (TH1D*) currentFinalFile[i]->Get("DoubleRatioConversionTrueEffPurityLow");
		if(i == 1 && cutVariationName.CompareTo("CocktailEtaNorm") == 0) DirectPhotonDoubleRatio[i] = (TH1D*) currentFinalFile[i]->Get("DoubleRatioConversionTrueEffPurityEtaHigh");
		if(i == 2 && cutVariationName.CompareTo("CocktailEtaNorm") == 0) DirectPhotonDoubleRatio[i] = (TH1D*) currentFinalFile[i]->Get("DoubleRatioConversionTrueEffPurityEtaLow");
		if(i == 1 && cutVariationName.Contains("CocktailParam")) DirectPhotonDoubleRatio[i] = (TH1D*) currentFinalFile[i]->Get("DoubleRatioConversionTrueEffPurityModA");
		if(i == 2 && cutVariationName.Contains("CocktailParam")) DirectPhotonDoubleRatio[i] = (TH1D*) currentFinalFile[i]->Get("DoubleRatioConversionTrueEffPurityModB");
		DirectPhotonDoubleRatio[i]->SetTitle("");
		DrawGammaSetMarker(DirectPhotonDoubleRatio[i], markerType, 2.0, color[i], color[i]);
		DirectPhotonDoubleRatioRatio[i] = (TH1D*) DirectPhotonDoubleRatio[i]->Clone(Form("DirectPhotonDoubleRatioRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
		DirectPhotonDoubleRatioRatio[i]->Divide(DirectPhotonDoubleRatioRatio[i],DirectPhotonDoubleRatio[0],1,1,"b");
		SetHistogramm(DirectPhotonDoubleRatioRatio[i],"p_{T} (GeV/c)","Ratios of Double Ratios",0,2);

		DirectPhotonDoubleRatioFit[i] = (TH1D*) currentFinalFile[i]->Get("DoubleRatioConversionFitPurity");
		if(i == 1 && cutVariationName.Contains("IntRange")) DirectPhotonDoubleRatioFit[i] = (TH1D*) currentFinalFile[i]->Get("DoubleRatioConversionFitPurityWide");
		if(i == 2 && cutVariationName.Contains("IntRange")) DirectPhotonDoubleRatioFit[i] = (TH1D*) currentFinalFile[i]->Get("DoubleRatioConversionFitPurityNarrow");
		if(i == 1 && cutVariationName.Contains("Fit")) DirectPhotonDoubleRatioFit[i] = (TH1D*) currentFinalFile[i]->Get("DoubleRatioConversionLowFitPurity");
		if(i == 2 && cutVariationName.Contains("Fit")) DirectPhotonDoubleRatioFit[i] = (TH1D*) currentFinalFile[i]->Get("DoubleRatioConversionHighFitPurity");
		if(i == 1 && cutVariationName.CompareTo("CocktailEta") == 0) DirectPhotonDoubleRatioFit[i] = (TH1D*) currentFinalFile[i]->Get("DoubleRatioConversionFitPurityHigh");
		if(i == 2 && cutVariationName.CompareTo("CocktailEta") == 0) DirectPhotonDoubleRatioFit[i] = (TH1D*) currentFinalFile[i]->Get("DoubleRatioConversionFitPurityLow");
		if(i == 1 && cutVariationName.CompareTo("CocktailEtaNorm") == 0) DirectPhotonDoubleRatioFit[i] = (TH1D*) currentFinalFile[i]->Get("DoubleRatioConversionFitPurityEtaHigh");
		if(i == 2 && cutVariationName.CompareTo("CocktailEtaNorm") == 0) DirectPhotonDoubleRatioFit[i] = (TH1D*) currentFinalFile[i]->Get("DoubleRatioConversionFitPurityEtaLow");
		if(i == 1 && cutVariationName.Contains("CocktailParam")) DirectPhotonDoubleRatioFit[i] = (TH1D*) currentFinalFile[i]->Get("DoubleRatioConversionFitPurityModA");
		if(i == 2 && cutVariationName.Contains("CocktailParam")) DirectPhotonDoubleRatioFit[i] = (TH1D*) currentFinalFile[i]->Get("DoubleRatioConversionFitPurityModB");
		DirectPhotonDoubleRatioFit[i]->SetTitle("");
		DrawGammaSetMarker(DirectPhotonDoubleRatioFit[i], markerType, 2.0, color[i], color[i]);
		DirectPhotonDoubleRatioFitRatio[i] = (TH1D*) DirectPhotonDoubleRatioFit[i]->Clone(Form("DirectPhotonDoubleRatioFitRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
		DirectPhotonDoubleRatioFitRatio[i]->Divide(DirectPhotonDoubleRatioFitRatio[i],DirectPhotonDoubleRatioFit[0],1,1,"b");
		SetHistogramm(DirectPhotonDoubleRatioFitRatio[i],"p_{T} (GeV/c)","Ratios of Double Ratios",0,2);

		// DirectPhotonSpectrum[i] = (TH1D*) currentFinalFile[i]->Get("histoDirectPhotonSpectrum");
		// DirectPhotonSpectrum[i]->SetTitle("");
		// DrawGammaSetMarker(DirectPhotonSpectrum[i], markerType, 2.0, color[i], color[i]);
		// DirectPhotonSpectrumRatio[i] = (TH1D*) DirectPhotonSpectrum[i]->Clone(Form("DirectPhotonSpectrumRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
		// DirectPhotonSpectrumRatio[i]->Divide(DirectPhotonSpectrumRatio[i],DirectPhotonSpectrum[0],1,1,"b");
		// SetHistogramm(DirectPhotonSpectrumRatio[i],"p_{T} (GeV/c)","Ratios of Direct Photon Spectra",0,2.5);


		PlotCanvas(i,number,GammaSpectrumCanvas,InclusiveGammaSpectrum[i],GammaSpectrumLegend,GammaSpectrumRatioCanvas,InclusiveGammaSpectrumRatio[i],GammaSpectrumRatioLegend,cutStringsName[i],One);
		PlotCanvas(i,number,Pi0SpectrumCanvas,Pi0Spectrum[i],Pi0SpectrumLegend,Pi0SpectrumRatioCanvas,Pi0SpectrumRatio[i],Pi0SpectrumRatioLegend,cutStringsName[i],One);
		PlotCanvas(i,number,Pi0SpectrumFitCanvas,Pi0SpectrumFit[i],Pi0SpectrumFitLegend,Pi0SpectrumFitRatioCanvas,Pi0SpectrumFitRatio[i],Pi0SpectrumFitRatioLegend,cutStringsName[i],One);
		PlotCanvas(i,number,InclusiveGammaToPi0RatioCanvas,InclusiveGammaToPi0Ratio[i],InclusiveGammaToPi0RatioLegend,InclusiveGammaToPi0RatioRatioCanvas,InclusiveGammaToPi0RatioRatio[i],InclusiveGammaToPi0RatioRatioLegend,cutStringsName[i],One);
		PlotCanvas(i,number,InclusiveGammaToPi0RatioFitCanvas,InclusiveGammaToPi0RatioFit[i],InclusiveGammaToPi0RatioFitLegend,InclusiveGammaToPi0RatioFitRatioCanvas,InclusiveGammaToPi0RatioFitRatio[i],InclusiveGammaToPi0RatioFitRatioLegend,cutStringsName[i],One);
		PlotCanvas(i,number,DirectPhotonDoubleRatioCanvas,DirectPhotonDoubleRatio[i],DirectPhotonDoubleRatioLegend,DirectPhotonDoubleRatioRatioCanvas,DirectPhotonDoubleRatioRatio[i],DirectPhotonDoubleRatioRatioLegend,cutStringsName[i],One,One);
		PlotCanvas(i,number,DirectPhotonDoubleRatioFitCanvas,DirectPhotonDoubleRatioFit[i],DirectPhotonDoubleRatioFitLegend,DirectPhotonDoubleRatioFitRatioCanvas,DirectPhotonDoubleRatioFitRatio[i],DirectPhotonDoubleRatioFitRatioLegend,cutStringsName[i],One,One);
		//PlotCanvas(i,number,DirectPhotonSpectrumCanvas,DirectPhotonSpectrum[i],DirectPhotonSpectrumLegend,DirectPhotonSpectrumRatioCanvas,DirectPhotonSpectrumRatio[i],DirectPhotonSpectrumRatioLegend,cutStringsName[i],One);
		
		currentCorrectedFile = Form("%s/%s/Gamma_Pi0_data_GammaConvV1Correction_%s.root",cutSelection[i].Data(),energy.Data(),cutSelection[i].Data());
		currentCorrectionFile[i] = new TFile(currentCorrectedFile);

		RawGamma[i] = (TH1D*) currentCorrectionFile[i]->Get("RawGammaSpectrum");
		RawGamma[i]->SetTitle("");
		DrawGammaSetMarker(RawGamma[i], markerType, 2.0, color[i], color[i]);
		RawGammaRatio[i] = (TH1D*) RawGamma[i]->Clone(Form("RawGammaRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
		RawGammaRatio[i]->Divide(RawGammaRatio[i],RawGamma[0],1,1,"b");
		SetHistogramm(RawGammaRatio[i],"p_{T} (GeV/c)","Ratios of Raw #gamma",0,2);

		Purity[i] = (TH1D*) currentCorrectionFile[i]->Get("MCGammaTruePurity");
		Purity[i]->SetTitle("");
		DrawGammaSetMarker(Purity[i], markerType, 2.0, color[i], color[i]);
		PurityRatio[i] = (TH1D*) Purity[i]->Clone(Form("PurityRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
		PurityRatio[i]->Divide(PurityRatio[i],Purity[0],1,1,"b");
		SetHistogramm(PurityRatio[i],"p_{T} (GeV/c)","Ratios of Purities",0.8,1.2);
			
		GammaEff[i] = (TH1D*) currentCorrectionFile[i]->Get("MCGammaPrimaryRecoEffMCPt");
		GammaEff[i]->SetTitle("");
		DrawGammaSetMarker(GammaEff[i], markerType, 2.0, color[i], color[i]);
		GammaEffRatio[i] = (TH1D*) GammaEff[i]->Clone(Form("GammaEffRatio_%s/%s",cutSelection[i].Data(),cutSelection[0].Data()));
		GammaEffRatio[i]->Divide(GammaEffRatio[i],GammaEff[0],1,1,"b");
		SetHistogramm(GammaEffRatio[i],"p_{T} (GeV/c)","Ratios of Purities",0,2);

		PlotCanvas(i,number,RawGammaCanvas,RawGamma[i],RawGammaLegend,RawGammaRatioCanvas,RawGammaRatio[i],RawGammaRatioLegend,cutStringsName[i],One);
		PlotCanvas(i,number,PurityCanvas,Purity[i],PurityLegend,PurityRatioCanvas,PurityRatio[i],PurityRatioLegend,cutStringsName[i],One);
		PlotCanvas(i,number,GammaEffCanvas,GammaEff[i],GammaEffLegend,GammaEffRatioCanvas,GammaEffRatio[i],GammaEffRatioLegend,cutStringsName[i],One);
		
	}

	GammaSpectrumCanvas->SaveAs(Form("%s/%s_GammaSpectra_%s.eps",outputDir.Data(),cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()));
	GammaSpectrumRatioCanvas->SaveAs(Form("%s/%s_GammaSpectraRatios_%s.eps",outputDir.Data(),cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()));
	Pi0SpectrumCanvas->SaveAs(Form("%s/%s_Pi0Spectra_%s.eps",outputDir.Data(),cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()));
	Pi0SpectrumRatioCanvas->SaveAs(Form("%s/%s_Pi0SpectraRatios_%s.eps",outputDir.Data(),cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()));
	Pi0SpectrumFitCanvas->SaveAs(Form("%s/%s_Pi0SpectraFit_%s.eps",outputDir.Data(),cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()));
	Pi0SpectrumFitRatioCanvas->SaveAs(Form("%s/%s_Pi0SpectraFitRatios_%s.eps",outputDir.Data(),cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()));
	InclusiveGammaToPi0RatioCanvas->SaveAs(Form("%s/%s_InclusiveRatios_%s.eps",outputDir.Data(),cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()));
	InclusiveGammaToPi0RatioRatioCanvas->SaveAs(Form("%s/%s_InclusiveRatiosRatios_%s.eps",outputDir.Data(),cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()));
	InclusiveGammaToPi0RatioFitCanvas->SaveAs(Form("%s/%s_InclusiveRatiosFit_%s.eps",outputDir.Data(),cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()));
	InclusiveGammaToPi0RatioFitRatioCanvas->SaveAs(Form("%s/%s_InclusiveRatiosFitRatios_%s.eps",outputDir.Data(),cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()));
	DirectPhotonDoubleRatioCanvas->SaveAs(Form("%s/%s_DoubleRatios_%s.eps",outputDir.Data(),cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()));
	DirectPhotonDoubleRatioRatioCanvas->SaveAs(Form("%s/%s_DoubleRatiosRatios_%s.eps",outputDir.Data(),cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()));
	DirectPhotonDoubleRatioFitCanvas->SaveAs(Form("%s/%s_DoubleRatiosFit_%s.eps",outputDir.Data(),cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()));
	DirectPhotonDoubleRatioFitRatioCanvas->SaveAs(Form("%s/%s_DoubleRatiosFitRatios_%s.eps",outputDir.Data(),cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()));
	// DirectPhotonSpectrumCanvas->SaveAs(Form("%s/%s_DirectPhotonSpectrum_%s.eps",outputDir.Data(),cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()));
	// DirectPhotonSpectrumRatioCanvas->SaveAs(Form("%s/%s_DirectPhotonSpectrumRatios_%s.eps",outputDir.Data(),cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()));
	RawGammaCanvas->SaveAs(Form("%s/%s_RawGamma_%s.eps",outputDir.Data(),cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()));
	RawGammaRatioCanvas->SaveAs(Form("%s/%s_RawGammaRatios_%s.eps",outputDir.Data(),cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()));
	PurityCanvas->SaveAs(Form("%s/%s_Purity_%s.eps",outputDir.Data(),cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()));
	PurityRatioCanvas->SaveAs(Form("%s/%s_PurityRatios_%s.eps",outputDir.Data(),cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()));
	GammaEffCanvas->SaveAs(Form("%s/%s_GammaEff_%s.eps",outputDir.Data(),cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()));
	GammaEffRatioCanvas->SaveAs(Form("%s/%s_GammaEffRatios_%s.eps",outputDir.Data(),cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()));
	


	// Start Systematic Error Calculation

	Int_t NBinsPt = InclusiveGammaSpectrum[0]->GetNbinsX();
	const Int_t NBinstPtConst = NBinsPt+1;

	const Int_t ConstNumberOfCuts = number;

	Double_t BinValueGamma[NBinstPtConst];
	BinValueGamma[0]=0.;

	Double_t  BinsXCenter[NBinstPtConst];
	Double_t  BinsXWidth[NBinstPtConst];
	Double_t BinValue[NBinstPtConst];
	BinsXCenter[0] = 0;
	BinsXWidth[0]=0.;
	BinValue[0]=0.;
	for (Int_t i = 1; i < NBinsPt +1; i++){
		BinsXCenter[i] = InclusiveGammaSpectrum[0]->GetBinCenter(i);
		BinsXWidth[i]= InclusiveGammaSpectrum[0]->GetBinWidth(i)/2.;
	}

	Double_t BinValueGamma[NBinstPtConst];
	BinValueGamma[0]=0.;

	SysErrorConversion SysErrCutGamma[ConstNumberOfCuts][NBinstPtConst];
	SysErrorConversion SysErrCutGammaRaw[ConstNumberOfCuts][NBinstPtConst];

	for (Int_t j = 0; j < number; j++){
		for (Int_t i = 1; i < NBinsPt +1; i++){
			BinValueGamma[i]= InclusiveGammaSpectrum[0]->GetBinContent(i);
			SysErrCutGamma[j][i].value = InclusiveGammaSpectrum[j]->GetBinContent(i);
			SysErrCutGamma[j][i].error = InclusiveGammaSpectrum[j]->GetBinError(i);
			SysErrCutGammaRaw[j][i].value = RawGamma[j]->GetBinContent(i);
			SysErrCutGammaRaw[j][i].error = RawGamma[j]->GetBinError(i);
		}
	}

	Double_t DifferenceCutGamma[ConstNumberOfCuts][NBinstPtConst];
	Double_t DifferenceErrorCutGamma[ConstNumberOfCuts][NBinstPtConst];

	Double_t LargestDiffGammaNeg[NBinstPtConst];
	Double_t LargestDiffGammaPos[NBinstPtConst];
	Double_t LargestDiffGammaErrorNeg[NBinstPtConst];
	Double_t LargestDiffGammaErrorPos[NBinstPtConst];

	Double_t LargestDiffGammaRelNeg[NBinstPtConst];
	Double_t LargestDiffGammaRelPos[NBinstPtConst];
	Double_t LargestDiffGammaRelErrorNeg[NBinstPtConst];
	Double_t LargestDiffGammaRelErrorPos[NBinstPtConst];

	Double_t RelDifferenceCutGamma[ConstNumberOfCuts][NBinstPtConst];
	Double_t RelDifferenceErrorCutGamma[ConstNumberOfCuts][NBinstPtConst];
	Double_t RelDifferenceRawCut[ConstNumberOfCuts][NBinstPtConst];
			
	for (Int_t j = 1; j < number; j++){
		for ( Int_t i = 1; i < NBinstPtConst; i++) {
			DifferenceCutGamma[j][i]=0.;
			DifferenceErrorCutGamma[j][i]=0.;
			LargestDiffGammaNeg[i]=0.;
			LargestDiffGammaPos[i]=0.;
			LargestDiffGammaErrorNeg[i]=0.;
			LargestDiffGammaErrorPos[i]=0.;
			RelDifferenceCutGamma[j][i]=0.;
			RelDifferenceRawCut[j][i]=0.;
			RelDifferenceErrorCutGamma[j][i]=0.;
		}
	}

	for(Int_t j = 1; j < number; j++){
		for (Int_t i = 1; i < NBinsPt +1; i++){
			//Calculate differences
			DifferenceCutGamma[j][i] = SysErrCutGamma[j][i].value - SysErrCutGamma[0][i].value;
			DifferenceErrorCutGamma[j][i] = TMath::Sqrt(TMath::Abs(TMath::Power(SysErrCutGamma[j][i].error,2)-TMath::Power(SysErrCutGamma[0][i].error,2)));
			if(SysErrCutGamma[0][i].value != 0){
				RelDifferenceCutGamma[j][i] = DifferenceCutGamma[j][i]/SysErrCutGamma[0][i].value*100. ;
				RelDifferenceErrorCutGamma[j][i] = DifferenceErrorCutGamma[j][i]/SysErrCutGamma[0][i].value*100. ;
			} else {
				RelDifferenceCutGamma[j][i] = -10000.;
				RelDifferenceErrorCutGamma[j][i] = 100. ;
			}
			if(SysErrCutGammaRaw[0][i].value != 0){
				RelDifferenceRawCut[j][i] = (SysErrCutGammaRaw[j][i].value - SysErrCutGammaRaw[0][i].value)/SysErrCutGammaRaw[0][i].value*100. ;
			} else {
				RelDifferenceRawCut[j][i] = -10000.;
			}
					
			if(DifferenceCutGamma[j][i] < 0){
				if (TMath::Abs(LargestDiffGammaNeg[i]) < TMath::Abs(DifferenceCutGamma[j][i]) && RelDifferenceRawCut[j][i] > -75.){
				LargestDiffGammaNeg[i] = DifferenceCutGamma[j][i];
				LargestDiffGammaErrorNeg[i] = DifferenceErrorCutGamma[j][i];
				}
			}else{
				if (TMath::Abs(LargestDiffGammaPos[i]) < TMath::Abs(DifferenceCutGamma[j][i]) && RelDifferenceRawCut[j][i] > -75.){
				LargestDiffGammaPos[i] = DifferenceCutGamma[j][i];
				LargestDiffGammaErrorPos[i] = DifferenceErrorCutGamma[j][i];
				}
			}
		}
	}

	cout << "done filling" << endl;
	const char *SysErrDatnameGamma = Form("%s/Gamma_SystematicErrorCutStudies.dat",outputDir.Data());
	fstream SysErrDatGamma;
	SysErrDatGamma.open(SysErrDatnameGamma, ios::out);
	SysErrDatGamma << "Calculation of the systematic error due to the yield cuts" << endl;

	cout << "works" << endl;
	for (Int_t l=0; l< number; l++){
		if (l == 0) {
			SysErrDatGamma << endl <<"Bin" << "\t" << cutSelection[l] << "\t" <<endl;
			for(Int_t i = 1; i < (NBinsPt +1); i++){
				SysErrDatGamma << BinsXCenter[i] << "\t" << SysErrCutGamma[l][i].value << "\t" << SysErrCutGamma[l][i].error << endl;	
			}
		} else{
			SysErrDatGamma << endl <<"Bin" << "\t" << cutSelection[l] << "\t" << "Error " << "\t Dif to Cut1" << endl;
			for(Int_t i = 1; i < (NBinsPt +1); i++){
				if (RelDifferenceRawCut[l][i] > -75.){
				SysErrDatGamma << BinsXCenter[i] << "\t" << SysErrCutGamma[l][i].value << "\t" << SysErrCutGamma[l][i].error << "\t" <<  DifferenceCutGamma[l][i] << "\t"<< DifferenceErrorCutGamma[l][i] << "\t"<< RelDifferenceCutGamma[l][i] <<  "\t" << RelDifferenceErrorCutGamma[l][i] <<"\t" << RelDifferenceRawCut[l][i]<< endl;
				} else {
				SysErrDatGamma << BinsXCenter[i] << "\t" << SysErrCutGamma[l][i].value << "\t" << SysErrCutGamma[l][i].error << "\t" <<  DifferenceCutGamma[l][i] << "\t"<< DifferenceErrorCutGamma[l][i] << "\t"<< RelDifferenceCutGamma[l][i] <<  "\t" << RelDifferenceErrorCutGamma[l][i] <<"\t" << RelDifferenceRawCut[l][i]  <<"\t not considered in largest dev" <<endl;
				}
			}
		}
	}


	SysErrDatGamma << endl;
	SysErrDatGamma << endl;
	SysErrDatGamma << "Bin" << "\t" << "Largest Dev Neg" << "\t" << "Largest Dev Pos"  << endl;
	for(Int_t i = 1; i < (NBinsPt +1); i++){
		SysErrDatGamma << BinsXCenter[i]  << "\t" << LargestDiffGammaNeg[i] << "\t" <<LargestDiffGammaErrorNeg[i]<< "\t" << LargestDiffGammaPos[i] << "\t" << LargestDiffGammaErrorPos[i]<<endl;
	}
	SysErrDatGamma << endl << endl <<"Bin" << "\t" << "Largest Dev Neg rel" << "\t" << "Largest Dev Pos rel"  << endl;
	for(Int_t i = 0; i < (NBinsPt +1); i++){
		if ( SysErrCutGamma[0][i].value != 0.){
			LargestDiffGammaRelNeg[i] = - LargestDiffGammaNeg[i]/SysErrCutGamma[0][i].value*100.;
			LargestDiffGammaRelPos[i] = LargestDiffGammaPos[i]/SysErrCutGamma[0][i].value*100.;
			LargestDiffGammaRelErrorNeg[i] = - LargestDiffGammaErrorNeg[i]/SysErrCutGamma[0][i].value*100.;
			LargestDiffGammaRelErrorPos[i] = LargestDiffGammaErrorPos[i]/SysErrCutGamma[0][i].value*100.;
			if (i > 0) SysErrDatGamma << BinsXCenter[i] << "\t" << LargestDiffGammaNeg[i]/SysErrCutGamma[0][i].value*100. << "\t" << LargestDiffGammaErrorNeg[i]/SysErrCutGamma[0][i].value*100. << "\t" << LargestDiffGammaPos[i]/SysErrCutGamma[0][i].value*100. << "\t" << LargestDiffGammaErrorPos[i]/SysErrCutGamma[0][i].value*100.<<endl;
		} else {
			LargestDiffGammaRelNeg[i] = 0.;
			LargestDiffGammaRelPos[i] = 0.;
			LargestDiffGammaRelErrorNeg[i] = 0.;
			LargestDiffGammaRelErrorPos[i] = 0.;
		}
	}

	SysErrDatGamma.close();
				
	TGraphAsymmErrors* SystErrGraphNegGamma = new TGraphAsymmErrors(NBinsPt+1, BinsXCenter, LargestDiffGammaRelNeg, BinsXWidth, BinsXWidth, LargestDiffGammaRelErrorNeg, LargestDiffGammaRelErrorNeg);
	SystErrGraphNegGamma->SetName(Form("Gamma_SystErrorRelNeg_%s",cutVariationName.Data()));
	TGraphAsymmErrors* SystErrGraphPosGamma = new TGraphAsymmErrors(NBinsPt+1, BinsXCenter, LargestDiffGammaRelPos, BinsXWidth, BinsXWidth, LargestDiffGammaRelErrorPos, LargestDiffGammaRelErrorPos);
	SystErrGraphPosGamma->SetName(Form("Gamma_SystErrorRelPos_%s",cutVariationName.Data()));
	const char* Outputname = Form("%s/Gamma_SystematicErrorCuts.root",outputFileDir.Data());
	TFile* SystematicErrorFile = new TFile(Outputname,"UPDATE");
	SystErrGraphPosGamma->Write(Form("Gamma_SystErrorRelPos_%s_%s",cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()),TObject::kOverwrite);
	SystErrGraphNegGamma->Write(Form("Gamma_SystErrorRelNeg_%s_%s",cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()),TObject::kOverwrite);
	SystematicErrorFile->Write();
	SystematicErrorFile->Close();
		


	//*************************************************************************************************
	//******************** Output of the systematic Error due to Signal extraction for DoubleRatio ****
	//*************************************************************************************************

	Int_t NumberOfCutsDoubleRatio = number;
	if(cutVariationName.CompareTo("CocktailEtaNorm") == 0 || cutVariationName.CompareTo("CocktailEta") == 0 || cutVariationName.Contains("CocktailParam") || cutVariationName.CompareTo("Yield") == 0 || cutVariationName.Contains("Fit") || cutVariationName.CompareTo("Purity") == 0){
		NumberOfCutsDoubleRatio = 3;
	}
	if(cutVariationName.CompareTo("Charged") == 0){
		NumberOfCutsDoubleRatio = 2;
	}

	const Int_t ConstNumberOfCutsDoubleRatio = NumberOfCutsDoubleRatio;
		
	Int_t NBinsPtdoubleR = DirectPhotonDoubleRatio[0]->GetNbinsX();
		
	const Int_t NBinstPtConstdoubleR = NBinsPtdoubleR+1;
		
	Double_t BinValueDoubleRatio[NBinstPtConstdoubleR];
	BinValueDoubleRatio[0]=0.;
		
	Double_t  BinsXCenterDoubleR[NBinstPtConstdoubleR];
	Double_t  BinsXWidthDoubleR[NBinstPtConstdoubleR];
	BinsXCenterDoubleR[0] = 0;
	BinsXWidthDoubleR[0]=0.;
	BinValue[0]=0.;
		
		
	for (Int_t i = 1; i <  NBinsPtdoubleR+1; i++){
		BinsXCenterDoubleR[i] = DirectPhotonDoubleRatio[0]->GetBinCenter(i);
		BinsXWidthDoubleR[i]= DirectPhotonDoubleRatio[0]->GetBinWidth(i)/2.;
	}
		
	SysErrorConversion SysErrCutDoubleRatio[ConstNumberOfCutsDoubleRatio][NBinstPtConstdoubleR];

	for (Int_t j = 0; j < NumberOfCutsDoubleRatio; j++){
		for (Int_t i = 1; i < NBinsPtdoubleR+1; i++){
			BinValueDoubleRatio[i]=DirectPhotonDoubleRatio[0]->GetBinContent(i);
			SysErrCutDoubleRatio[j][i].value =DirectPhotonDoubleRatio[j]->GetBinContent(i);
			SysErrCutDoubleRatio[j][i].error =DirectPhotonDoubleRatio[j]->GetBinError(i);
		}
	}
	Double_t DifferenceCutDoubleRatio[ConstNumberOfCutsDoubleRatio][NBinstPtConstdoubleR];
	Double_t DifferenceErrorCutDoubleRatio[ConstNumberOfCutsDoubleRatio][NBinstPtConstdoubleR];

	Double_t LargestDiffDoubleRatioNeg[NBinstPtConstdoubleR];
	Double_t LargestDiffDoubleRatioPos[NBinstPtConstdoubleR];
	Double_t LargestDiffDoubleRatioErrorNeg[NBinstPtConstdoubleR];
	Double_t LargestDiffDoubleRatioErrorPos[NBinstPtConstdoubleR];

	Double_t LargestDiffDoubleRatioRelNeg[NBinstPtConstdoubleR];
	Double_t LargestDiffDoubleRatioRelPos[NBinstPtConstdoubleR];
	Double_t LargestDiffDoubleRatioRelErrorNeg[NBinstPtConstdoubleR];
	Double_t LargestDiffDoubleRatioRelErrorPos[NBinstPtConstdoubleR];

	Double_t RelDifferenceCutDoubleRatio[ConstNumberOfCutsDoubleRatio][NBinstPtConstdoubleR];
	Double_t RelDifferenceErrorCutDoubleRatio[ConstNumberOfCutsDoubleRatio][NBinstPtConstdoubleR];

	for (Int_t j = 1; j < NumberOfCutsDoubleRatio; j++){
		for ( Int_t i = 1; i < NBinstPtConstdoubleR; i++) {
			DifferenceCutDoubleRatio[j][i]=0.;
			DifferenceErrorCutDoubleRatio[j][i]=0.;
			LargestDiffDoubleRatioNeg[i]=0.;
			LargestDiffDoubleRatioPos[i]=0.;
			LargestDiffDoubleRatioErrorNeg[i]=0.;
			LargestDiffDoubleRatioErrorPos[i]=0.;
			RelDifferenceCutDoubleRatio[j][i]=0.;
			RelDifferenceErrorCutDoubleRatio[j][i]=0.;
		}
	}
		
	for(Int_t j = 1; j < NumberOfCutsDoubleRatio; j++){
		for (Int_t i = 1; i < NBinsPtdoubleR +1; i++){
			//Calculate differences
			DifferenceCutDoubleRatio[j][i] = SysErrCutDoubleRatio[j][i].value - SysErrCutDoubleRatio[0][i].value;
			DifferenceErrorCutDoubleRatio[j][i] = TMath::Sqrt(TMath::Abs(TMath::Power(SysErrCutDoubleRatio[j][i].error,2)-TMath::Power(SysErrCutDoubleRatio[0][i].error,2)));
			if(SysErrCutDoubleRatio[0][i].value != 0){
				RelDifferenceCutDoubleRatio[j][i] = DifferenceCutDoubleRatio[j][i]/SysErrCutDoubleRatio[0][i].value*100. ;
				RelDifferenceErrorCutDoubleRatio[j][i] = DifferenceErrorCutDoubleRatio[j][i]/SysErrCutDoubleRatio[0][i].value*100. ;
			} else {
				RelDifferenceCutDoubleRatio[j][i] = -10000.;
				RelDifferenceErrorCutDoubleRatio[j][i] = 100. ;
			}
						
			if(DifferenceCutDoubleRatio[j][i] < 0){
				if (TMath::Abs(LargestDiffDoubleRatioNeg[i]) < TMath::Abs(DifferenceCutDoubleRatio[j][i])){
				LargestDiffDoubleRatioNeg[i] = DifferenceCutDoubleRatio[j][i];
				LargestDiffDoubleRatioErrorNeg[i] = DifferenceErrorCutDoubleRatio[j][i];
				}
			}else{
				if (TMath::Abs(LargestDiffDoubleRatioPos[i]) < TMath::Abs(DifferenceCutDoubleRatio[j][i])){
				LargestDiffDoubleRatioPos[i] = DifferenceCutDoubleRatio[j][i];
				LargestDiffDoubleRatioErrorPos[i] = DifferenceErrorCutDoubleRatio[j][i];
				}
			}
		}
	}
		
	cout << "done filling" << endl;
	const char *SysErrDatnameDoubleRatio = Form("%s/DoubleRatio_SystematicErrorCutStudies.dat",outputDir.Data());
	fstream SysErrDatDoubleRatio;
	SysErrDatDoubleRatio.open(SysErrDatnameDoubleRatio, ios::out);
	SysErrDatDoubleRatio << "Calculation of the systematic error due to the yield cuts" << endl;

	cout << "works" << endl;
	for (Int_t l=0; l< NumberOfCutsDoubleRatio; l++){
		if (l == 0) {
			SysErrDatDoubleRatio << endl <<"Bin" << "\t" << cutSelection[l] << "\t" <<endl;
			for(Int_t i = 1; i < (NBinsPtdoubleR +1); i++){
				SysErrDatDoubleRatio << BinsXCenterDoubleR[i] << "\t" << SysErrCutDoubleRatio[l][i].value << "\t" << SysErrCutDoubleRatio[l][i].error << endl;	
			}
		} else{
			SysErrDatDoubleRatio << endl <<"Bin" << "\t" << cutSelection[l] << "\t" << "Error " << "\t Dif to Cut1" << endl;
			for(Int_t i = 1; i < (NBinsPtdoubleR +1); i++){
				SysErrDatDoubleRatio << BinsXCenterDoubleR[i] << "\t" << SysErrCutDoubleRatio[l][i].value << "\t" << SysErrCutDoubleRatio[l][i].error << "\t" <<  DifferenceCutDoubleRatio[l][i] << "\t"<< DifferenceErrorCutDoubleRatio[l][i] << "\t"<< RelDifferenceCutDoubleRatio[l][i] <<  "\t" << RelDifferenceErrorCutDoubleRatio[l][i] << endl;
			}
		}
	}


	SysErrDatDoubleRatio << endl;
	SysErrDatDoubleRatio << endl;
	SysErrDatDoubleRatio << "Bin" << "\t" << "Largest Dev Neg" << "\t" << "Largest Dev Pos"  << endl;
	for(Int_t i = 1; i < (NBinsPtdoubleR +1); i++){
		SysErrDatDoubleRatio << BinsXCenterDoubleR[i]  << "\t" << LargestDiffDoubleRatioNeg[i] << "\t" <<LargestDiffDoubleRatioErrorNeg[i]<< "\t" << LargestDiffDoubleRatioPos[i] << "\t" << LargestDiffDoubleRatioErrorPos[i]<<endl;
	}
	SysErrDatDoubleRatio << endl << endl <<"Bin" << "\t" << "Largest Dev Neg rel" << "\t" << "Largest Dev Pos rel"  << endl;
	for(Int_t i = 0; i < (NBinsPtdoubleR +1); i++){
		if ( SysErrCutDoubleRatio[0][i].value != 0.){
			LargestDiffDoubleRatioRelNeg[i] = - LargestDiffDoubleRatioNeg[i]/SysErrCutDoubleRatio[0][i].value*100.;
			LargestDiffDoubleRatioRelPos[i] = LargestDiffDoubleRatioPos[i]/SysErrCutDoubleRatio[0][i].value*100.;
			LargestDiffDoubleRatioRelErrorNeg[i] = - LargestDiffDoubleRatioErrorNeg[i]/SysErrCutDoubleRatio[0][i].value*100.;
			LargestDiffDoubleRatioRelErrorPos[i] = LargestDiffDoubleRatioErrorPos[i]/SysErrCutDoubleRatio[0][i].value*100.;
			if (i > 0) SysErrDatDoubleRatio << BinsXCenterDoubleR[i] << "\t" << LargestDiffDoubleRatioNeg[i]/SysErrCutDoubleRatio[0][i].value*100. << "\t" << LargestDiffDoubleRatioErrorNeg[i]/SysErrCutDoubleRatio[0][i].value*100. << "\t" << LargestDiffDoubleRatioPos[i]/SysErrCutDoubleRatio[0][i].value*100. << "\t" << LargestDiffDoubleRatioErrorPos[i]/SysErrCutDoubleRatio[0][i].value*100.<<endl;
		} else {
			LargestDiffDoubleRatioRelNeg[i] = 0.;
			LargestDiffDoubleRatioRelPos[i] = 0.;
			LargestDiffDoubleRatioRelErrorNeg[i] = 0.;
			LargestDiffDoubleRatioRelErrorPos[i] = 0.;
		}
	}

	SysErrDatDoubleRatio.close();
				
	TGraphAsymmErrors* SystErrGraphNegDoubleRatio = new TGraphAsymmErrors(NBinsPtdoubleR+1, BinsXCenterDoubleR, LargestDiffDoubleRatioRelNeg, BinsXWidthDoubleR, BinsXWidthDoubleR, LargestDiffDoubleRatioRelErrorNeg, LargestDiffDoubleRatioRelErrorNeg);
	SystErrGraphNegDoubleRatio->SetName(Form("DoubleRatio_SystErrorRelNeg_%s",cutVariationName.Data()));
	TGraphAsymmErrors* SystErrGraphPosDoubleRatio = new TGraphAsymmErrors(NBinsPtdoubleR+1, BinsXCenterDoubleR, LargestDiffDoubleRatioRelPos, BinsXWidthDoubleR, BinsXWidthDoubleR, LargestDiffDoubleRatioRelErrorPos, LargestDiffDoubleRatioRelErrorPos);
	SystErrGraphPosDoubleRatio->SetName(Form("DoubleRatio_SystErrorRelPos_%s",cutVariationName.Data()));
	const char* Outputname = Form("%s/Gamma_SystematicErrorCuts.root",outputFileDir.Data());
	TFile* SystematicErrorFile = new TFile(Outputname,"UPDATE");
	SystErrGraphPosDoubleRatio->Write(Form("DoubleRatio_SystErrorRelPos_%s_%s",cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()),TObject::kOverwrite);
	SystErrGraphNegDoubleRatio->Write(Form("DoubleRatio_SystErrorRelNeg_%s_%s",cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()),TObject::kOverwrite);
	SystematicErrorFile->Write();
	SystematicErrorFile->Close();
		
	
	//*************************************************************************************************
	//***************** Output of the systematic Error due to Signal extraction for DoubleRatioFit ****
	//*************************************************************************************************

	Int_t NBinsPtdoubleR = DirectPhotonDoubleRatioFit[0]->GetNbinsX();
	const Int_t NBinstPtConstdoubleR = NBinsPtdoubleR+1;
		
	Double_t BinValueDoubleRatioFit[NBinstPtConstdoubleR];
	BinValueDoubleRatioFit[0]=0.;

	Double_t  BinsXCenterDoubleR[NBinstPtConstdoubleR];
	Double_t  BinsXWidthDoubleR[NBinstPtConstdoubleR];
	BinsXCenterDoubleR[0] = 0;
	BinsXWidthDoubleR[0]=0.;
	BinValue[0]=0.;
	for (Int_t i = 1; i <  NBinsPtdoubleR+1; i++){
		BinsXCenterDoubleR[i] = DirectPhotonDoubleRatioFit[0]->GetBinCenter(i);
		BinsXWidthDoubleR[i]= DirectPhotonDoubleRatioFit[0]->GetBinWidth(i)/2.;
	}

	SysErrorConversion SysErrCutDoubleRatioFit[ConstNumberOfCutsDoubleRatio][NBinstPtConstdoubleR];

	for (Int_t j = 0; j < NumberOfCutsDoubleRatio; j++){
		for (Int_t i = 1; i <  NBinsPtdoubleR+1; i++){
			BinValueDoubleRatioFit[i]=DirectPhotonDoubleRatioFit[0]->GetBinContent(i);
			SysErrCutDoubleRatioFit[j][i].value =DirectPhotonDoubleRatioFit[j]->GetBinContent(i);
			SysErrCutDoubleRatioFit[j][i].error =DirectPhotonDoubleRatioFit[j]->GetBinError(i);
		}
	}

	Double_t DifferenceCutDoubleRatioFit[ConstNumberOfCutsDoubleRatio][NBinstPtConstdoubleR];
	Double_t DifferenceErrorCutDoubleRatioFit[ConstNumberOfCutsDoubleRatio][NBinstPtConstdoubleR];

	Double_t LargestDiffDoubleRatioFitNeg[NBinstPtConstdoubleR];
	Double_t LargestDiffDoubleRatioFitPos[NBinstPtConstdoubleR];
	Double_t LargestDiffDoubleRatioFitErrorNeg[NBinstPtConstdoubleR];
	Double_t LargestDiffDoubleRatioFitErrorPos[NBinstPtConstdoubleR];

	Double_t LargestDiffDoubleRatioFitRelNeg[NBinstPtConstdoubleR];
	Double_t LargestDiffDoubleRatioFitRelPos[NBinstPtConstdoubleR];
	Double_t LargestDiffDoubleRatioFitRelErrorNeg[NBinstPtConstdoubleR];
	Double_t LargestDiffDoubleRatioFitRelErrorPos[NBinstPtConstdoubleR];

	Double_t RelDifferenceCutDoubleRatioFit[ConstNumberOfCutsDoubleRatio][NBinstPtConstdoubleR];
	Double_t RelDifferenceErrorCutDoubleRatioFit[ConstNumberOfCutsDoubleRatio][NBinstPtConstdoubleR];
				
	for (Int_t j = 1; j < NumberOfCutsDoubleRatio; j++){
		for ( Int_t i = 1; i < NBinstPtConstdoubleR; i++) {
			DifferenceCutDoubleRatioFit[j][i]=0.;
			DifferenceErrorCutDoubleRatioFit[j][i]=0.;
			LargestDiffDoubleRatioFitNeg[i]=0.;
			LargestDiffDoubleRatioFitPos[i]=0.;
			LargestDiffDoubleRatioFitErrorNeg[i]=0.;
			LargestDiffDoubleRatioFitErrorPos[i]=0.;
			RelDifferenceCutDoubleRatioFit[j][i]=0.;
			RelDifferenceErrorCutDoubleRatioFit[j][i]=0.;
		}
	}
		
	for(Int_t j = 1; j < NumberOfCutsDoubleRatio; j++){
		for (Int_t i = 1; i < NBinsPtdoubleR +1; i++){
			//Calculate differences
			DifferenceCutDoubleRatioFit[j][i] = SysErrCutDoubleRatioFit[j][i].value - SysErrCutDoubleRatioFit[0][i].value;
			DifferenceErrorCutDoubleRatioFit[j][i] = TMath::Sqrt(TMath::Abs(TMath::Power(SysErrCutDoubleRatioFit[j][i].error,2)-TMath::Power(SysErrCutDoubleRatioFit[0][i].error,2)));
			if(SysErrCutDoubleRatioFit[0][i].value != 0){
				RelDifferenceCutDoubleRatioFit[j][i] = DifferenceCutDoubleRatioFit[j][i]/SysErrCutDoubleRatioFit[0][i].value*100. ;
				RelDifferenceErrorCutDoubleRatioFit[j][i] = DifferenceErrorCutDoubleRatioFit[j][i]/SysErrCutDoubleRatioFit[0][i].value*100. ;
			} else {
				RelDifferenceCutDoubleRatioFit[j][i] = -10000.;
				RelDifferenceErrorCutDoubleRatioFit[j][i] = 100. ;
			}
						
			if(DifferenceCutDoubleRatioFit[j][i] < 0){
				if (TMath::Abs(LargestDiffDoubleRatioFitNeg[i]) < TMath::Abs(DifferenceCutDoubleRatioFit[j][i])){
				LargestDiffDoubleRatioFitNeg[i] = DifferenceCutDoubleRatioFit[j][i];
				LargestDiffDoubleRatioFitErrorNeg[i] = DifferenceErrorCutDoubleRatioFit[j][i];
				}
			}else{
				if (TMath::Abs(LargestDiffDoubleRatioFitPos[i]) < TMath::Abs(DifferenceCutDoubleRatioFit[j][i])){
				LargestDiffDoubleRatioFitPos[i] = DifferenceCutDoubleRatioFit[j][i];
				LargestDiffDoubleRatioFitErrorPos[i] = DifferenceErrorCutDoubleRatioFit[j][i];
				}
			}
		}
	}
		
	cout << "done filling" << endl;
	const char *SysErrDatnameDoubleRatioFit = Form("%s/DoubleRatioFit_SystematicErrorCutStudies.dat",outputDir.Data());
	fstream SysErrDatDoubleRatioFit;
	SysErrDatDoubleRatioFit.open(SysErrDatnameDoubleRatioFit, ios::out);
	SysErrDatDoubleRatioFit << "Calculation of the systematic error due to the yield cuts" << endl;

	cout << "works" << endl;
	for (Int_t l=0; l< NumberOfCutsDoubleRatio; l++){
		if (l == 0) {
			SysErrDatDoubleRatioFit << endl <<"Bin" << "\t" << cutSelection[l] << "\t" <<endl;
			for(Int_t i = 1; i < (NBinsPtdoubleR +1); i++){
				SysErrDatDoubleRatioFit << BinsXCenterDoubleR[i] << "\t" << SysErrCutDoubleRatioFit[l][i].value << "\t" << SysErrCutDoubleRatioFit[l][i].error << endl;	
			}
		} else{
			SysErrDatDoubleRatioFit << endl <<"Bin" << "\t" << cutSelection[l] << "\t" << "Error " << "\t Dif to Cut1" << endl;
			for(Int_t i = 1; i < (NBinsPtdoubleR +1); i++){
				SysErrDatDoubleRatioFit << BinsXCenterDoubleR[i] << "\t" << SysErrCutDoubleRatioFit[l][i].value << "\t" << SysErrCutDoubleRatioFit[l][i].error << "\t" <<  DifferenceCutDoubleRatioFit[l][i] << "\t"<< DifferenceErrorCutDoubleRatioFit[l][i] << "\t"<< RelDifferenceCutDoubleRatioFit[l][i] <<  "\t" << RelDifferenceErrorCutDoubleRatioFit[l][i] << endl;
			}
		}
	}


	SysErrDatDoubleRatioFit << endl;
	SysErrDatDoubleRatioFit << endl;
	SysErrDatDoubleRatioFit << "Bin" << "\t" << "Largest Dev Neg" << "\t" << "Largest Dev Pos"  << endl;
	for(Int_t i = 1; i < (NBinsPtdoubleR +1); i++){
		SysErrDatDoubleRatioFit << BinsXCenterDoubleR[i]  << "\t" << LargestDiffDoubleRatioFitNeg[i] << "\t" <<LargestDiffDoubleRatioFitErrorNeg[i]<< "\t" << LargestDiffDoubleRatioFitPos[i] << "\t" << LargestDiffDoubleRatioFitErrorPos[i]<<endl;
	}
	SysErrDatDoubleRatioFit << endl << endl <<"Bin" << "\t" << "Largest Dev Neg rel" << "\t" << "Largest Dev Pos rel"  << endl;
	for(Int_t i = 0; i < (NBinsPtdoubleR +1); i++){
		if ( SysErrCutDoubleRatioFit[0][i].value != 0.){
			LargestDiffDoubleRatioFitRelNeg[i] = - LargestDiffDoubleRatioFitNeg[i]/SysErrCutDoubleRatioFit[0][i].value*100.;
			LargestDiffDoubleRatioFitRelPos[i] = LargestDiffDoubleRatioFitPos[i]/SysErrCutDoubleRatioFit[0][i].value*100.;
			LargestDiffDoubleRatioFitRelErrorNeg[i] = - LargestDiffDoubleRatioFitErrorNeg[i]/SysErrCutDoubleRatioFit[0][i].value*100.;
			LargestDiffDoubleRatioFitRelErrorPos[i] = LargestDiffDoubleRatioFitErrorPos[i]/SysErrCutDoubleRatioFit[0][i].value*100.;
			if (i > 0) SysErrDatDoubleRatioFit << BinsXCenterDoubleR[i] << "\t" << LargestDiffDoubleRatioFitNeg[i]/SysErrCutDoubleRatioFit[0][i].value*100. << "\t" << LargestDiffDoubleRatioFitErrorNeg[i]/SysErrCutDoubleRatioFit[0][i].value*100. << "\t" << LargestDiffDoubleRatioFitPos[i]/SysErrCutDoubleRatioFit[0][i].value*100. << "\t" << LargestDiffDoubleRatioFitErrorPos[i]/SysErrCutDoubleRatioFit[0][i].value*100.<<endl;
		} else {
			LargestDiffDoubleRatioFitRelNeg[i] = 0.;
			LargestDiffDoubleRatioFitRelPos[i] = 0.;
			LargestDiffDoubleRatioFitRelErrorNeg[i] = 0.;
			LargestDiffDoubleRatioFitRelErrorPos[i] = 0.;
		}
	}

	SysErrDatDoubleRatioFit.close();
				
	TGraphAsymmErrors* SystErrGraphNegDoubleRatioFit = new TGraphAsymmErrors(NBinsPtdoubleR+1, BinsXCenterDoubleR, LargestDiffDoubleRatioFitRelNeg, BinsXWidthDoubleR, BinsXWidthDoubleR, LargestDiffDoubleRatioFitRelErrorNeg, LargestDiffDoubleRatioFitRelErrorNeg);
	SystErrGraphNegDoubleRatioFit->SetName(Form("DoubleRatioFit_SystErrorRelNeg_%s",cutVariationName.Data()));
	TGraphAsymmErrors* SystErrGraphPosDoubleRatioFit = new TGraphAsymmErrors(NBinsPtdoubleR+1, BinsXCenterDoubleR, LargestDiffDoubleRatioFitRelPos, BinsXWidthDoubleR, BinsXWidthDoubleR, LargestDiffDoubleRatioFitRelErrorPos, LargestDiffDoubleRatioFitRelErrorPos);
	SystErrGraphPosDoubleRatioFit->SetName(Form("DoubleRatioFit_SystErrorRelPos_%s",cutVariationName.Data()));
	const char* Outputname = Form("%s/Gamma_SystematicErrorCuts.root",outputFileDir.Data());
	TFile* SystematicErrorFile = new TFile(Outputname,"UPDATE");
	SystErrGraphPosDoubleRatioFit->Write(Form("DoubleRatioFit_SystErrorRelPos_%s_%s",cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()),TObject::kOverwrite);
	SystErrGraphNegDoubleRatioFit->Write(Form("DoubleRatioFit_SystErrorRelNeg_%s_%s",cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()),TObject::kOverwrite);
	SystematicErrorFile->Write();
	SystematicErrorFile->Close();


	//*************************************************************************************************
	//******************** Output of the systematic Error due to Signal extraction for IncRatio ****
	//*************************************************************************************************

	Int_t NumberOfCutsIncRatio = number;
	if(cutVariationName.CompareTo("Yield") == 0 || cutVariationName.Contains("Fit") || cutVariationName.CompareTo("Purity") == 0){
		NumberOfCutsIncRatio = 3;
	}

	const Int_t ConstNumberOfCutsIncRatio = NumberOfCutsIncRatio;
		
	Int_t NBinsPtdoubleR = InclusiveGammaToPi0Ratio[0]->GetNbinsX();
		
	const Int_t NBinstPtConstdoubleR = NBinsPtdoubleR+1;
		
	Double_t BinValueIncRatio[NBinstPtConstdoubleR];
	BinValueIncRatio[0]=0.;
		
	Double_t  BinsXCenterIncR[NBinstPtConstdoubleR];
	Double_t  BinsXWidthIncR[NBinstPtConstdoubleR];
	BinsXCenterIncR[0] = 0;
	BinsXWidthIncR[0]=0.;
	BinValue[0]=0.;
		
		
	for (Int_t i = 1; i <  NBinsPtdoubleR+1; i++){
		BinsXCenterIncR[i] = InclusiveGammaToPi0Ratio[0]->GetBinCenter(i);
		BinsXWidthIncR[i]= InclusiveGammaToPi0Ratio[0]->GetBinWidth(i)/2.;
	}
		
	SysErrorConversion SysErrCutIncRatio[ConstNumberOfCutsIncRatio][NBinstPtConstdoubleR];

	for (Int_t j = 0; j < NumberOfCutsIncRatio; j++){
		for (Int_t i = 1; i < NBinsPtdoubleR+1; i++){
			BinValueIncRatio[i]=InclusiveGammaToPi0Ratio[0]->GetBinContent(i);
			SysErrCutIncRatio[j][i].value =InclusiveGammaToPi0Ratio[j]->GetBinContent(i);
			SysErrCutIncRatio[j][i].error =InclusiveGammaToPi0Ratio[j]->GetBinError(i);
		}
	}
	Double_t DifferenceCutIncRatio[ConstNumberOfCutsIncRatio][NBinstPtConstdoubleR];
	Double_t DifferenceErrorCutIncRatio[ConstNumberOfCutsIncRatio][NBinstPtConstdoubleR];

	Double_t LargestDiffIncRatioNeg[NBinstPtConstdoubleR];
	Double_t LargestDiffIncRatioPos[NBinstPtConstdoubleR];
	Double_t LargestDiffIncRatioErrorNeg[NBinstPtConstdoubleR];
	Double_t LargestDiffIncRatioErrorPos[NBinstPtConstdoubleR];

	Double_t LargestDiffIncRatioRelNeg[NBinstPtConstdoubleR];
	Double_t LargestDiffIncRatioRelPos[NBinstPtConstdoubleR];
	Double_t LargestDiffIncRatioRelErrorNeg[NBinstPtConstdoubleR];
	Double_t LargestDiffIncRatioRelErrorPos[NBinstPtConstdoubleR];

	Double_t RelDifferenceCutIncRatio[ConstNumberOfCutsIncRatio][NBinstPtConstdoubleR];
	Double_t RelDifferenceErrorCutIncRatio[ConstNumberOfCutsIncRatio][NBinstPtConstdoubleR];

	for (Int_t j = 1; j < NumberOfCutsIncRatio; j++){
		for ( Int_t i = 1; i < NBinstPtConstdoubleR; i++) {
			DifferenceCutIncRatio[j][i]=0.;
			DifferenceErrorCutIncRatio[j][i]=0.;
			LargestDiffIncRatioNeg[i]=0.;
			LargestDiffIncRatioPos[i]=0.;
			LargestDiffIncRatioErrorNeg[i]=0.;
			LargestDiffIncRatioErrorPos[i]=0.;
			RelDifferenceCutIncRatio[j][i]=0.;
			RelDifferenceErrorCutIncRatio[j][i]=0.;
		}
	}
		
	for(Int_t j = 1; j < NumberOfCutsIncRatio; j++){
		for (Int_t i = 1; i < NBinsPtdoubleR +1; i++){
			//Calculate differences
			DifferenceCutIncRatio[j][i] = SysErrCutIncRatio[j][i].value - SysErrCutIncRatio[0][i].value;
			DifferenceErrorCutIncRatio[j][i] = TMath::Sqrt(TMath::Abs(TMath::Power(SysErrCutIncRatio[j][i].error,2)-TMath::Power(SysErrCutIncRatio[0][i].error,2)));
			if(SysErrCutIncRatio[0][i].value != 0){
				RelDifferenceCutIncRatio[j][i] = DifferenceCutIncRatio[j][i]/SysErrCutIncRatio[0][i].value*100. ;
				RelDifferenceErrorCutIncRatio[j][i] = DifferenceErrorCutIncRatio[j][i]/SysErrCutIncRatio[0][i].value*100. ;
			} else {
				RelDifferenceCutIncRatio[j][i] = -10000.;
				RelDifferenceErrorCutIncRatio[j][i] = 100. ;
			}
						
			if(DifferenceCutIncRatio[j][i] < 0){
				if (TMath::Abs(LargestDiffIncRatioNeg[i]) < TMath::Abs(DifferenceCutIncRatio[j][i])){
				LargestDiffIncRatioNeg[i] = DifferenceCutIncRatio[j][i];
				LargestDiffIncRatioErrorNeg[i] = DifferenceErrorCutIncRatio[j][i];
				}
			}else{
				if (TMath::Abs(LargestDiffIncRatioPos[i]) < TMath::Abs(DifferenceCutIncRatio[j][i])){
				LargestDiffIncRatioPos[i] = DifferenceCutIncRatio[j][i];
				LargestDiffIncRatioErrorPos[i] = DifferenceErrorCutIncRatio[j][i];
				}
			}
		}
	}
		
	cout << "done filling" << endl;
	const char *SysErrDatnameIncRatio = Form("%s/IncRatio_SystematicErrorCutStudies.dat",outputDir.Data());
	fstream SysErrDatIncRatio;
	SysErrDatIncRatio.open(SysErrDatnameIncRatio, ios::out);
	SysErrDatIncRatio << "Calculation of the systematic error due to the yield cuts" << endl;

	cout << "works" << endl;
	for (Int_t l=0; l< NumberOfCutsIncRatio; l++){
		if (l == 0) {
			SysErrDatIncRatio << endl <<"Bin" << "\t" << cutSelection[l] << "\t" <<endl;
			for(Int_t i = 1; i < (NBinsPtdoubleR +1); i++){
				SysErrDatIncRatio << BinsXCenterIncR[i] << "\t" << SysErrCutIncRatio[l][i].value << "\t" << SysErrCutIncRatio[l][i].error << endl;	
			}
		} else{
			SysErrDatIncRatio << endl <<"Bin" << "\t" << cutSelection[l] << "\t" << "Error " << "\t Dif to Cut1" << endl;
			for(Int_t i = 1; i < (NBinsPtdoubleR +1); i++){
				SysErrDatIncRatio << BinsXCenterIncR[i] << "\t" << SysErrCutIncRatio[l][i].value << "\t" << SysErrCutIncRatio[l][i].error << "\t" <<  DifferenceCutIncRatio[l][i] << "\t"<< DifferenceErrorCutIncRatio[l][i] << "\t"<< RelDifferenceCutIncRatio[l][i] <<  "\t" << RelDifferenceErrorCutIncRatio[l][i] << endl;
			}
		}
	}


	SysErrDatIncRatio << endl;
	SysErrDatIncRatio << endl;
	SysErrDatIncRatio << "Bin" << "\t" << "Largest Dev Neg" << "\t" << "Largest Dev Pos"  << endl;
	for(Int_t i = 1; i < (NBinsPtdoubleR +1); i++){
		SysErrDatIncRatio << BinsXCenterIncR[i]  << "\t" << LargestDiffIncRatioNeg[i] << "\t" <<LargestDiffIncRatioErrorNeg[i]<< "\t" << LargestDiffIncRatioPos[i] << "\t" << LargestDiffIncRatioErrorPos[i]<<endl;
	}
	SysErrDatIncRatio << endl << endl <<"Bin" << "\t" << "Largest Dev Neg rel" << "\t" << "Largest Dev Pos rel"  << endl;
	for(Int_t i = 0; i < (NBinsPtdoubleR +1); i++){
		if ( SysErrCutIncRatio[0][i].value != 0.){
			LargestDiffIncRatioRelNeg[i] = - LargestDiffIncRatioNeg[i]/SysErrCutIncRatio[0][i].value*100.;
			LargestDiffIncRatioRelPos[i] = LargestDiffIncRatioPos[i]/SysErrCutIncRatio[0][i].value*100.;
			LargestDiffIncRatioRelErrorNeg[i] = - LargestDiffIncRatioErrorNeg[i]/SysErrCutIncRatio[0][i].value*100.;
			LargestDiffIncRatioRelErrorPos[i] = LargestDiffIncRatioErrorPos[i]/SysErrCutIncRatio[0][i].value*100.;
			if (i > 0) SysErrDatIncRatio << BinsXCenterIncR[i] << "\t" << LargestDiffIncRatioNeg[i]/SysErrCutIncRatio[0][i].value*100. << "\t" << LargestDiffIncRatioErrorNeg[i]/SysErrCutIncRatio[0][i].value*100. << "\t" << LargestDiffIncRatioPos[i]/SysErrCutIncRatio[0][i].value*100. << "\t" << LargestDiffIncRatioErrorPos[i]/SysErrCutIncRatio[0][i].value*100.<<endl;
		} else {
			LargestDiffIncRatioRelNeg[i] = 0.;
			LargestDiffIncRatioRelPos[i] = 0.;
			LargestDiffIncRatioRelErrorNeg[i] = 0.;
			LargestDiffIncRatioRelErrorPos[i] = 0.;
		}
	}

	SysErrDatIncRatio.close();
				
	TGraphAsymmErrors* SystErrGraphNegIncRatio = new TGraphAsymmErrors(NBinsPtdoubleR+1, BinsXCenterIncR, LargestDiffIncRatioRelNeg, BinsXWidthIncR, BinsXWidthIncR, LargestDiffIncRatioRelErrorNeg, LargestDiffIncRatioRelErrorNeg);
	SystErrGraphNegIncRatio->SetName(Form("IncRatio_SystErrorRelNeg_%s",cutVariationName.Data()));
	TGraphAsymmErrors* SystErrGraphPosIncRatio = new TGraphAsymmErrors(NBinsPtdoubleR+1, BinsXCenterIncR, LargestDiffIncRatioRelPos, BinsXWidthIncR, BinsXWidthIncR, LargestDiffIncRatioRelErrorPos, LargestDiffIncRatioRelErrorPos);
	SystErrGraphPosIncRatio->SetName(Form("IncRatio_SystErrorRelPos_%s",cutVariationName.Data()));
	const char* Outputname = Form("%s/Gamma_SystematicErrorCuts.root",outputFileDir.Data());
	TFile* SystematicErrorFile = new TFile(Outputname,"UPDATE");
	SystErrGraphPosIncRatio->Write(Form("IncRatio_SystErrorRelPos_%s_%s",cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()),TObject::kOverwrite);
	SystErrGraphNegIncRatio->Write(Form("IncRatio_SystErrorRelNeg_%s_%s",cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()),TObject::kOverwrite);
	if(cutVariationName.CompareTo("Purity") == 0){
		SystErrGraphPosIncRatio->Write(Form("Gamma_SystErrorRelPos_%s_%s",cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()),TObject::kOverwrite);
		SystErrGraphNegIncRatio->Write(Form("Gamma_SystErrorRelNeg_%s_%s",cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()),TObject::kOverwrite);
	}

	SystematicErrorFile->Write();
	SystematicErrorFile->Close();
		
	
	//*************************************************************************************************
	//******************** Output of the systematic Error due to Signal extraction for IncRatioFit ****
	//*************************************************************************************************

	Int_t NBinsPtdoubleR = InclusiveGammaToPi0RatioFit[0]->GetNbinsX();
	const Int_t NBinstPtConstdoubleR = NBinsPtdoubleR+1;
		
	Double_t BinValueIncRatioFit[NBinstPtConstdoubleR];
	BinValueIncRatioFit[0]=0.;

	Double_t  BinsXCenterIncR[NBinstPtConstdoubleR];
	Double_t  BinsXWidthIncR[NBinstPtConstdoubleR];
	BinsXCenterIncR[0] = 0;
	BinsXWidthIncR[0]=0.;
	BinValue[0]=0.;
	for (Int_t i = 1; i <  NBinsPtdoubleR+1; i++){
		BinsXCenterIncR[i] = InclusiveGammaToPi0RatioFit[0]->GetBinCenter(i);
		BinsXWidthIncR[i]= InclusiveGammaToPi0RatioFit[0]->GetBinWidth(i)/2.;
	}

	SysErrorConversion SysErrCutIncRatioFit[ConstNumberOfCutsIncRatio][NBinstPtConstdoubleR];

	for (Int_t j = 0; j < NumberOfCutsIncRatio; j++){
		for (Int_t i = 1; i <  NBinsPtdoubleR+1; i++){
			BinValueIncRatioFit[i]=InclusiveGammaToPi0RatioFit[0]->GetBinContent(i);
			SysErrCutIncRatioFit[j][i].value =InclusiveGammaToPi0RatioFit[j]->GetBinContent(i);
			SysErrCutIncRatioFit[j][i].error =InclusiveGammaToPi0RatioFit[j]->GetBinError(i);
		}
	}

	Double_t DifferenceCutIncRatioFit[ConstNumberOfCutsIncRatio][NBinstPtConstdoubleR];
	Double_t DifferenceErrorCutIncRatioFit[ConstNumberOfCutsIncRatio][NBinstPtConstdoubleR];

	Double_t LargestDiffIncRatioFitNeg[NBinstPtConstdoubleR];
	Double_t LargestDiffIncRatioFitPos[NBinstPtConstdoubleR];
	Double_t LargestDiffIncRatioFitErrorNeg[NBinstPtConstdoubleR];
	Double_t LargestDiffIncRatioFitErrorPos[NBinstPtConstdoubleR];

	Double_t LargestDiffIncRatioFitRelNeg[NBinstPtConstdoubleR];
	Double_t LargestDiffIncRatioFitRelPos[NBinstPtConstdoubleR];
	Double_t LargestDiffIncRatioFitRelErrorNeg[NBinstPtConstdoubleR];
	Double_t LargestDiffIncRatioFitRelErrorPos[NBinstPtConstdoubleR];

	Double_t RelDifferenceCutIncRatioFit[ConstNumberOfCutsIncRatio][NBinstPtConstdoubleR];
	Double_t RelDifferenceErrorCutIncRatioFit[ConstNumberOfCutsIncRatio][NBinstPtConstdoubleR];
				
	for (Int_t j = 1; j < NumberOfCutsIncRatio; j++){
		for ( Int_t i = 1; i < NBinstPtConstdoubleR; i++) {
			DifferenceCutIncRatioFit[j][i]=0.;
			DifferenceErrorCutIncRatioFit[j][i]=0.;
			LargestDiffIncRatioFitNeg[i]=0.;
			LargestDiffIncRatioFitPos[i]=0.;
			LargestDiffIncRatioFitErrorNeg[i]=0.;
			LargestDiffIncRatioFitErrorPos[i]=0.;
			RelDifferenceCutIncRatioFit[j][i]=0.;
			RelDifferenceErrorCutIncRatioFit[j][i]=0.;
		}
	}
		
	for(Int_t j = 1; j < NumberOfCutsIncRatio; j++){
		for (Int_t i = 1; i < NBinsPtdoubleR +1; i++){
			//Calculate differences
			DifferenceCutIncRatioFit[j][i] = SysErrCutIncRatioFit[j][i].value - SysErrCutIncRatioFit[0][i].value;
			DifferenceErrorCutIncRatioFit[j][i] = TMath::Sqrt(TMath::Abs(TMath::Power(SysErrCutIncRatioFit[j][i].error,2)-TMath::Power(SysErrCutIncRatioFit[0][i].error,2)));
			if(SysErrCutIncRatioFit[0][i].value != 0){
				RelDifferenceCutIncRatioFit[j][i] = DifferenceCutIncRatioFit[j][i]/SysErrCutIncRatioFit[0][i].value*100. ;
				RelDifferenceErrorCutIncRatioFit[j][i] = DifferenceErrorCutIncRatioFit[j][i]/SysErrCutIncRatioFit[0][i].value*100. ;
			} else {
				RelDifferenceCutIncRatioFit[j][i] = -10000.;
				RelDifferenceErrorCutIncRatioFit[j][i] = 100. ;
			}
						
			if(DifferenceCutIncRatioFit[j][i] < 0){
				if (TMath::Abs(LargestDiffIncRatioFitNeg[i]) < TMath::Abs(DifferenceCutIncRatioFit[j][i])){
				LargestDiffIncRatioFitNeg[i] = DifferenceCutIncRatioFit[j][i];
				LargestDiffIncRatioFitErrorNeg[i] = DifferenceErrorCutIncRatioFit[j][i];
				}
			}else{
				if (TMath::Abs(LargestDiffIncRatioFitPos[i]) < TMath::Abs(DifferenceCutIncRatioFit[j][i])){
				LargestDiffIncRatioFitPos[i] = DifferenceCutIncRatioFit[j][i];
				LargestDiffIncRatioFitErrorPos[i] = DifferenceErrorCutIncRatioFit[j][i];
				}
			}
		}
	}
		
	cout << "done filling" << endl;
	const char *SysErrDatnameIncRatioFit = Form("%s/IncRatioFit_SystematicErrorCutStudies.dat",outputDir.Data());
	fstream SysErrDatIncRatioFit;
	SysErrDatIncRatioFit.open(SysErrDatnameIncRatioFit, ios::out);
	SysErrDatIncRatioFit << "Calculation of the systematic error due to the yield cuts" << endl;

	cout << "works" << endl;
	for (Int_t l=0; l< NumberOfCutsIncRatio; l++){
		if (l == 0) {
			SysErrDatIncRatioFit << endl <<"Bin" << "\t" << cutSelection[l] << "\t" <<endl;
			for(Int_t i = 1; i < (NBinsPtdoubleR +1); i++){
				SysErrDatIncRatioFit << BinsXCenterIncR[i] << "\t" << SysErrCutIncRatioFit[l][i].value << "\t" << SysErrCutIncRatioFit[l][i].error << endl;	
			}
		} else{
			SysErrDatIncRatioFit << endl <<"Bin" << "\t" << cutSelection[l] << "\t" << "Error " << "\t Dif to Cut1" << endl;
			for(Int_t i = 1; i < (NBinsPtdoubleR +1); i++){
				SysErrDatIncRatioFit << BinsXCenterIncR[i] << "\t" << SysErrCutIncRatioFit[l][i].value << "\t" << SysErrCutIncRatioFit[l][i].error << "\t" <<  DifferenceCutIncRatioFit[l][i] << "\t"<< DifferenceErrorCutIncRatioFit[l][i] << "\t"<< RelDifferenceCutIncRatioFit[l][i] <<  "\t" << RelDifferenceErrorCutIncRatioFit[l][i] << endl;
			}
		}
	}


	SysErrDatIncRatioFit << endl;
	SysErrDatIncRatioFit << endl;
	SysErrDatIncRatioFit << "Bin" << "\t" << "Largest Dev Neg" << "\t" << "Largest Dev Pos"  << endl;
	for(Int_t i = 1; i < (NBinsPtdoubleR +1); i++){
		SysErrDatIncRatioFit << BinsXCenterIncR[i]  << "\t" << LargestDiffIncRatioFitNeg[i] << "\t" <<LargestDiffIncRatioFitErrorNeg[i]<< "\t" << LargestDiffIncRatioFitPos[i] << "\t" << LargestDiffIncRatioFitErrorPos[i]<<endl;
	}
	SysErrDatIncRatioFit << endl << endl <<"Bin" << "\t" << "Largest Dev Neg rel" << "\t" << "Largest Dev Pos rel"  << endl;
	for(Int_t i = 0; i < (NBinsPtdoubleR +1); i++){
		if ( SysErrCutIncRatioFit[0][i].value != 0.){
			LargestDiffIncRatioFitRelNeg[i] = - LargestDiffIncRatioFitNeg[i]/SysErrCutIncRatioFit[0][i].value*100.;
			LargestDiffIncRatioFitRelPos[i] = LargestDiffIncRatioFitPos[i]/SysErrCutIncRatioFit[0][i].value*100.;
			LargestDiffIncRatioFitRelErrorNeg[i] = - LargestDiffIncRatioFitErrorNeg[i]/SysErrCutIncRatioFit[0][i].value*100.;
			LargestDiffIncRatioFitRelErrorPos[i] = LargestDiffIncRatioFitErrorPos[i]/SysErrCutIncRatioFit[0][i].value*100.;
			if (i > 0) SysErrDatIncRatioFit << BinsXCenterIncR[i] << "\t" << LargestDiffIncRatioFitNeg[i]/SysErrCutIncRatioFit[0][i].value*100. << "\t" << LargestDiffIncRatioFitErrorNeg[i]/SysErrCutIncRatioFit[0][i].value*100. << "\t" << LargestDiffIncRatioFitPos[i]/SysErrCutIncRatioFit[0][i].value*100. << "\t" << LargestDiffIncRatioFitErrorPos[i]/SysErrCutIncRatioFit[0][i].value*100.<<endl;
		} else {
			LargestDiffIncRatioFitRelNeg[i] = 0.;
			LargestDiffIncRatioFitRelPos[i] = 0.;
			LargestDiffIncRatioFitRelErrorNeg[i] = 0.;
			LargestDiffIncRatioFitRelErrorPos[i] = 0.;
		}
	}

	SysErrDatIncRatioFit.close();
				
	TGraphAsymmErrors* SystErrGraphNegIncRatioFit = new TGraphAsymmErrors(NBinsPtdoubleR+1, BinsXCenterIncR, LargestDiffIncRatioFitRelNeg, BinsXWidthIncR, BinsXWidthIncR, LargestDiffIncRatioFitRelErrorNeg, LargestDiffIncRatioFitRelErrorNeg);
	SystErrGraphNegIncRatioFit->SetName(Form("IncRatioFit_SystErrorRelNeg_%s",cutVariationName.Data()));
	TGraphAsymmErrors* SystErrGraphPosIncRatioFit = new TGraphAsymmErrors(NBinsPtdoubleR+1, BinsXCenterIncR, LargestDiffIncRatioFitRelPos, BinsXWidthIncR, BinsXWidthIncR, LargestDiffIncRatioFitRelErrorPos, LargestDiffIncRatioFitRelErrorPos);
	SystErrGraphPosIncRatioFit->SetName(Form("IncRatioFit_SystErrorRelPos_%s",cutVariationName.Data()));
	const char* Outputname = Form("%s/Gamma_SystematicErrorCuts.root",outputFileDir.Data());
	TFile* SystematicErrorFile = new TFile(Outputname,"UPDATE");
	SystErrGraphPosIncRatioFit->Write(Form("IncRatioFit_SystErrorRelPos_%s_%s",cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()),TObject::kOverwrite);
	SystErrGraphNegIncRatioFit->Write(Form("IncRatioFit_SystErrorRelNeg_%s_%s",cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()),TObject::kOverwrite);
	SystematicErrorFile->Write();
	SystematicErrorFile->Close();



	// -------------------------------------------------------------------------

	//*************************************************************************************************
	//******************** Output of the systematic Error due to Signal extraction for Pi0 ****
	//*************************************************************************************************

	Int_t NumberOfCutsPi0 = number;
	if(cutVariationName.CompareTo("Yield") == 0 || cutVariationName.Contains("Fit") || cutVariationName.CompareTo("Purity") == 0){
		NumberOfCutsPi0 = 3;
	}

	const Int_t ConstNumberOfCutsPi0 = NumberOfCutsPi0;
		
	Int_t NBinsPtdoubleR = InclusiveGammaToPi0Ratio[0]->GetNbinsX();
		
	const Int_t NBinstPtConstdoubleR = NBinsPtdoubleR+1;
		
	Double_t BinValuePi0[NBinstPtConstdoubleR];
	BinValuePi0[0]=0.;
		
	Double_t  BinsXCenterIncR[NBinstPtConstdoubleR];
	Double_t  BinsXWidthIncR[NBinstPtConstdoubleR];
	BinsXCenterIncR[0] = 0;
	BinsXWidthIncR[0]=0.;
	BinValue[0]=0.;
		
		
	for (Int_t i = 1; i <  NBinsPtdoubleR+1; i++){
		BinsXCenterIncR[i] = InclusiveGammaToPi0Ratio[0]->GetBinCenter(i);
		BinsXWidthIncR[i]= InclusiveGammaToPi0Ratio[0]->GetBinWidth(i)/2.;
	}
		
	SysErrorConversion SysErrCutPi0[ConstNumberOfCutsPi0][NBinstPtConstdoubleR];

	for (Int_t j = 0; j < NumberOfCutsPi0; j++){
		for (Int_t i = 1; i < NBinsPtdoubleR+1; i++){
			BinValuePi0[i]=InclusiveGammaToPi0Ratio[0]->GetBinContent(i);
			SysErrCutPi0[j][i].value =InclusiveGammaToPi0Ratio[j]->GetBinContent(i);
			SysErrCutPi0[j][i].error =InclusiveGammaToPi0Ratio[j]->GetBinError(i);
		}
	}
	Double_t DifferenceCutPi0[ConstNumberOfCutsPi0][NBinstPtConstdoubleR];
	Double_t DifferenceErrorCutPi0[ConstNumberOfCutsPi0][NBinstPtConstdoubleR];

	Double_t LargestDiffPi0Neg[NBinstPtConstdoubleR];
	Double_t LargestDiffPi0Pos[NBinstPtConstdoubleR];
	Double_t LargestDiffPi0ErrorNeg[NBinstPtConstdoubleR];
	Double_t LargestDiffPi0ErrorPos[NBinstPtConstdoubleR];

	Double_t LargestDiffPi0RelNeg[NBinstPtConstdoubleR];
	Double_t LargestDiffPi0RelPos[NBinstPtConstdoubleR];
	Double_t LargestDiffPi0RelErrorNeg[NBinstPtConstdoubleR];
	Double_t LargestDiffPi0RelErrorPos[NBinstPtConstdoubleR];

	Double_t RelDifferenceCutPi0[ConstNumberOfCutsPi0][NBinstPtConstdoubleR];
	Double_t RelDifferenceErrorCutPi0[ConstNumberOfCutsPi0][NBinstPtConstdoubleR];

	for (Int_t j = 1; j < NumberOfCutsPi0; j++){
		for ( Int_t i = 1; i < NBinstPtConstdoubleR; i++) {
			DifferenceCutPi0[j][i]=0.;
			DifferenceErrorCutPi0[j][i]=0.;
			LargestDiffPi0Neg[i]=0.;
			LargestDiffPi0Pos[i]=0.;
			LargestDiffPi0ErrorNeg[i]=0.;
			LargestDiffPi0ErrorPos[i]=0.;
			RelDifferenceCutPi0[j][i]=0.;
			RelDifferenceErrorCutPi0[j][i]=0.;
		}
	}
		
	for(Int_t j = 1; j < NumberOfCutsPi0; j++){
		for (Int_t i = 1; i < NBinsPtdoubleR +1; i++){
			//Calculate differences
			DifferenceCutPi0[j][i] = SysErrCutPi0[j][i].value - SysErrCutPi0[0][i].value;
			DifferenceErrorCutPi0[j][i] = TMath::Sqrt(TMath::Abs(TMath::Power(SysErrCutPi0[j][i].error,2)-TMath::Power(SysErrCutPi0[0][i].error,2)));
			if(SysErrCutPi0[0][i].value != 0){
				RelDifferenceCutPi0[j][i] = DifferenceCutPi0[j][i]/SysErrCutPi0[0][i].value*100. ;
				RelDifferenceErrorCutPi0[j][i] = DifferenceErrorCutPi0[j][i]/SysErrCutPi0[0][i].value*100. ;
			} else {
				RelDifferenceCutPi0[j][i] = -10000.;
				RelDifferenceErrorCutPi0[j][i] = 100. ;
			}
						
			if(DifferenceCutPi0[j][i] < 0){
				if (TMath::Abs(LargestDiffPi0Neg[i]) < TMath::Abs(DifferenceCutPi0[j][i])){
				LargestDiffPi0Neg[i] = DifferenceCutPi0[j][i];
				LargestDiffPi0ErrorNeg[i] = DifferenceErrorCutPi0[j][i];
				}
			}else{
				if (TMath::Abs(LargestDiffPi0Pos[i]) < TMath::Abs(DifferenceCutPi0[j][i])){
				LargestDiffPi0Pos[i] = DifferenceCutPi0[j][i];
				LargestDiffPi0ErrorPos[i] = DifferenceErrorCutPi0[j][i];
				}
			}
		}
	}
		
	cout << "done filling" << endl;
	const char *SysErrDatnamePi0 = Form("%s/Pi0_SystematicErrorCutStudies.dat",outputDir.Data());
	fstream SysErrDatPi0;
	SysErrDatPi0.open(SysErrDatnamePi0, ios::out);
	SysErrDatPi0 << "Calculation of the systematic error due to the yield cuts" << endl;

	cout << "works" << endl;
	for (Int_t l=0; l< NumberOfCutsPi0; l++){
		if (l == 0) {
			SysErrDatPi0 << endl <<"Bin" << "\t" << cutSelection[l] << "\t" <<endl;
			for(Int_t i = 1; i < (NBinsPtdoubleR +1); i++){
				SysErrDatPi0 << BinsXCenterIncR[i] << "\t" << SysErrCutPi0[l][i].value << "\t" << SysErrCutPi0[l][i].error << endl;	
			}
		} else{
			SysErrDatPi0 << endl <<"Bin" << "\t" << cutSelection[l] << "\t" << "Error " << "\t Dif to Cut1" << endl;
			for(Int_t i = 1; i < (NBinsPtdoubleR +1); i++){
				SysErrDatPi0 << BinsXCenterIncR[i] << "\t" << SysErrCutPi0[l][i].value << "\t" << SysErrCutPi0[l][i].error << "\t" <<  DifferenceCutPi0[l][i] << "\t"<< DifferenceErrorCutPi0[l][i] << "\t"<< RelDifferenceCutPi0[l][i] <<  "\t" << RelDifferenceErrorCutPi0[l][i] << endl;
			}
		}
	}


	SysErrDatPi0 << endl;
	SysErrDatPi0 << endl;
	SysErrDatPi0 << "Bin" << "\t" << "Largest Dev Neg" << "\t" << "Largest Dev Pos"  << endl;
	for(Int_t i = 1; i < (NBinsPtdoubleR +1); i++){
		SysErrDatPi0 << BinsXCenterIncR[i]  << "\t" << LargestDiffPi0Neg[i] << "\t" <<LargestDiffPi0ErrorNeg[i]<< "\t" << LargestDiffPi0Pos[i] << "\t" << LargestDiffPi0ErrorPos[i]<<endl;
	}
	SysErrDatPi0 << endl << endl <<"Bin" << "\t" << "Largest Dev Neg rel" << "\t" << "Largest Dev Pos rel"  << endl;
	for(Int_t i = 0; i < (NBinsPtdoubleR +1); i++){
		if ( SysErrCutPi0[0][i].value != 0.){
			LargestDiffPi0RelNeg[i] = - LargestDiffPi0Neg[i]/SysErrCutPi0[0][i].value*100.;
			LargestDiffPi0RelPos[i] = LargestDiffPi0Pos[i]/SysErrCutPi0[0][i].value*100.;
			LargestDiffPi0RelErrorNeg[i] = - LargestDiffPi0ErrorNeg[i]/SysErrCutPi0[0][i].value*100.;
			LargestDiffPi0RelErrorPos[i] = LargestDiffPi0ErrorPos[i]/SysErrCutPi0[0][i].value*100.;
			if (i > 0) SysErrDatPi0 << BinsXCenterIncR[i] << "\t" << LargestDiffPi0Neg[i]/SysErrCutPi0[0][i].value*100. << "\t" << LargestDiffPi0ErrorNeg[i]/SysErrCutPi0[0][i].value*100. << "\t" << LargestDiffPi0Pos[i]/SysErrCutPi0[0][i].value*100. << "\t" << LargestDiffPi0ErrorPos[i]/SysErrCutPi0[0][i].value*100.<<endl;
		} else {
			LargestDiffPi0RelNeg[i] = 0.;
			LargestDiffPi0RelPos[i] = 0.;
			LargestDiffPi0RelErrorNeg[i] = 0.;
			LargestDiffPi0RelErrorPos[i] = 0.;
		}
	}

	SysErrDatPi0.close();
				
	TGraphAsymmErrors* SystErrGraphNegPi0 = new TGraphAsymmErrors(NBinsPtdoubleR+1, BinsXCenterIncR, LargestDiffPi0RelNeg, BinsXWidthIncR, BinsXWidthIncR, LargestDiffPi0RelErrorNeg, LargestDiffPi0RelErrorNeg);
	SystErrGraphNegPi0->SetName(Form("Pi0_SystErrorRelNeg_%s",cutVariationName.Data()));
	TGraphAsymmErrors* SystErrGraphPosPi0 = new TGraphAsymmErrors(NBinsPtdoubleR+1, BinsXCenterIncR, LargestDiffPi0RelPos, BinsXWidthIncR, BinsXWidthIncR, LargestDiffPi0RelErrorPos, LargestDiffPi0RelErrorPos);
	SystErrGraphPosPi0->SetName(Form("Pi0_SystErrorRelPos_%s",cutVariationName.Data()));
	const char* Outputname = Form("%s/Gamma_SystematicErrorCuts.root",outputFileDir.Data());
	TFile* SystematicErrorFile = new TFile(Outputname,"UPDATE");
	SystErrGraphPosPi0->Write(Form("Pi0_SystErrorRelPos_%s_%s",cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()),TObject::kOverwrite);
	SystErrGraphNegPi0->Write(Form("Pi0_SystErrorRelNeg_%s_%s",cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()),TObject::kOverwrite);
	if(cutVariationName.CompareTo("Purity") == 0){
		SystErrGraphPosPi0->Write(Form("Gamma_SystErrorRelPos_%s_%s",cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()),TObject::kOverwrite);
		SystErrGraphNegPi0->Write(Form("Gamma_SystErrorRelNeg_%s_%s",cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()),TObject::kOverwrite);
	}

	SystematicErrorFile->Write();
	SystematicErrorFile->Close();
		
	
	//*************************************************************************************************
	//******************** Output of the systematic Error due to Signal extraction for Pi0Fit ****
	//*************************************************************************************************

	Int_t NBinsPtdoubleR = InclusiveGammaToPi0RatioFit[0]->GetNbinsX();
	const Int_t NBinstPtConstdoubleR = NBinsPtdoubleR+1;
		
	Double_t BinValuePi0Fit[NBinstPtConstdoubleR];
	BinValuePi0Fit[0]=0.;

	Double_t  BinsXCenterIncR[NBinstPtConstdoubleR];
	Double_t  BinsXWidthIncR[NBinstPtConstdoubleR];
	BinsXCenterIncR[0] = 0;
	BinsXWidthIncR[0]=0.;
	BinValue[0]=0.;
	for (Int_t i = 1; i <  NBinsPtdoubleR+1; i++){
		BinsXCenterIncR[i] = InclusiveGammaToPi0RatioFit[0]->GetBinCenter(i);
		BinsXWidthIncR[i]= InclusiveGammaToPi0RatioFit[0]->GetBinWidth(i)/2.;
	}

	SysErrorConversion SysErrCutPi0Fit[ConstNumberOfCutsPi0][NBinstPtConstdoubleR];

	for (Int_t j = 0; j < NumberOfCutsPi0; j++){
		for (Int_t i = 1; i <  NBinsPtdoubleR+1; i++){
			BinValuePi0Fit[i]=InclusiveGammaToPi0RatioFit[0]->GetBinContent(i);
			SysErrCutPi0Fit[j][i].value =InclusiveGammaToPi0RatioFit[j]->GetBinContent(i);
			SysErrCutPi0Fit[j][i].error =InclusiveGammaToPi0RatioFit[j]->GetBinError(i);
		}
	}

	Double_t DifferenceCutPi0Fit[ConstNumberOfCutsPi0][NBinstPtConstdoubleR];
	Double_t DifferenceErrorCutPi0Fit[ConstNumberOfCutsPi0][NBinstPtConstdoubleR];

	Double_t LargestDiffPi0FitNeg[NBinstPtConstdoubleR];
	Double_t LargestDiffPi0FitPos[NBinstPtConstdoubleR];
	Double_t LargestDiffPi0FitErrorNeg[NBinstPtConstdoubleR];
	Double_t LargestDiffPi0FitErrorPos[NBinstPtConstdoubleR];

	Double_t LargestDiffPi0FitRelNeg[NBinstPtConstdoubleR];
	Double_t LargestDiffPi0FitRelPos[NBinstPtConstdoubleR];
	Double_t LargestDiffPi0FitRelErrorNeg[NBinstPtConstdoubleR];
	Double_t LargestDiffPi0FitRelErrorPos[NBinstPtConstdoubleR];

	Double_t RelDifferenceCutPi0Fit[ConstNumberOfCutsPi0][NBinstPtConstdoubleR];
	Double_t RelDifferenceErrorCutPi0Fit[ConstNumberOfCutsPi0][NBinstPtConstdoubleR];
				
	for (Int_t j = 1; j < NumberOfCutsPi0; j++){
		for ( Int_t i = 1; i < NBinstPtConstdoubleR; i++) {
			DifferenceCutPi0Fit[j][i]=0.;
			DifferenceErrorCutPi0Fit[j][i]=0.;
			LargestDiffPi0FitNeg[i]=0.;
			LargestDiffPi0FitPos[i]=0.;
			LargestDiffPi0FitErrorNeg[i]=0.;
			LargestDiffPi0FitErrorPos[i]=0.;
			RelDifferenceCutPi0Fit[j][i]=0.;
			RelDifferenceErrorCutPi0Fit[j][i]=0.;
		}
	}
		
	for(Int_t j = 1; j < NumberOfCutsPi0; j++){
		for (Int_t i = 1; i < NBinsPtdoubleR +1; i++){
			//Calculate differences
			DifferenceCutPi0Fit[j][i] = SysErrCutPi0Fit[j][i].value - SysErrCutPi0Fit[0][i].value;
			DifferenceErrorCutPi0Fit[j][i] = TMath::Sqrt(TMath::Abs(TMath::Power(SysErrCutPi0Fit[j][i].error,2)-TMath::Power(SysErrCutPi0Fit[0][i].error,2)));
			if(SysErrCutPi0Fit[0][i].value != 0){
				RelDifferenceCutPi0Fit[j][i] = DifferenceCutPi0Fit[j][i]/SysErrCutPi0Fit[0][i].value*100. ;
				RelDifferenceErrorCutPi0Fit[j][i] = DifferenceErrorCutPi0Fit[j][i]/SysErrCutPi0Fit[0][i].value*100. ;
			} else {
				RelDifferenceCutPi0Fit[j][i] = -10000.;
				RelDifferenceErrorCutPi0Fit[j][i] = 100. ;
			}
						
			if(DifferenceCutPi0Fit[j][i] < 0){
				if (TMath::Abs(LargestDiffPi0FitNeg[i]) < TMath::Abs(DifferenceCutPi0Fit[j][i])){
				LargestDiffPi0FitNeg[i] = DifferenceCutPi0Fit[j][i];
				LargestDiffPi0FitErrorNeg[i] = DifferenceErrorCutPi0Fit[j][i];
				}
			}else{
				if (TMath::Abs(LargestDiffPi0FitPos[i]) < TMath::Abs(DifferenceCutPi0Fit[j][i])){
				LargestDiffPi0FitPos[i] = DifferenceCutPi0Fit[j][i];
				LargestDiffPi0FitErrorPos[i] = DifferenceErrorCutPi0Fit[j][i];
				}
			}
		}
	}
		
	cout << "done filling" << endl;
	const char *SysErrDatnamePi0Fit = Form("%s/Pi0Fit_SystematicErrorCutStudies.dat",outputDir.Data());
	fstream SysErrDatPi0Fit;
	SysErrDatPi0Fit.open(SysErrDatnamePi0Fit, ios::out);
	SysErrDatPi0Fit << "Calculation of the systematic error due to the yield cuts" << endl;

	cout << "works" << endl;
	for (Int_t l=0; l< NumberOfCutsPi0; l++){
		if (l == 0) {
			SysErrDatPi0Fit << endl <<"Bin" << "\t" << cutSelection[l] << "\t" <<endl;
			for(Int_t i = 1; i < (NBinsPtdoubleR +1); i++){
				SysErrDatPi0Fit << BinsXCenterIncR[i] << "\t" << SysErrCutPi0Fit[l][i].value << "\t" << SysErrCutPi0Fit[l][i].error << endl;	
			}
		} else{
			SysErrDatPi0Fit << endl <<"Bin" << "\t" << cutSelection[l] << "\t" << "Error " << "\t Dif to Cut1" << endl;
			for(Int_t i = 1; i < (NBinsPtdoubleR +1); i++){
				SysErrDatPi0Fit << BinsXCenterIncR[i] << "\t" << SysErrCutPi0Fit[l][i].value << "\t" << SysErrCutPi0Fit[l][i].error << "\t" <<  DifferenceCutPi0Fit[l][i] << "\t"<< DifferenceErrorCutPi0Fit[l][i] << "\t"<< RelDifferenceCutPi0Fit[l][i] <<  "\t" << RelDifferenceErrorCutPi0Fit[l][i] << endl;
			}
		}
	}


	SysErrDatPi0Fit << endl;
	SysErrDatPi0Fit << endl;
	SysErrDatPi0Fit << "Bin" << "\t" << "Largest Dev Neg" << "\t" << "Largest Dev Pos"  << endl;
	for(Int_t i = 1; i < (NBinsPtdoubleR +1); i++){
		SysErrDatPi0Fit << BinsXCenterIncR[i]  << "\t" << LargestDiffPi0FitNeg[i] << "\t" <<LargestDiffPi0FitErrorNeg[i]<< "\t" << LargestDiffPi0FitPos[i] << "\t" << LargestDiffPi0FitErrorPos[i]<<endl;
	}
	SysErrDatPi0Fit << endl << endl <<"Bin" << "\t" << "Largest Dev Neg rel" << "\t" << "Largest Dev Pos rel"  << endl;
	for(Int_t i = 0; i < (NBinsPtdoubleR +1); i++){
		if ( SysErrCutPi0Fit[0][i].value != 0.){
			LargestDiffPi0FitRelNeg[i] = - LargestDiffPi0FitNeg[i]/SysErrCutPi0Fit[0][i].value*100.;
			LargestDiffPi0FitRelPos[i] = LargestDiffPi0FitPos[i]/SysErrCutPi0Fit[0][i].value*100.;
			LargestDiffPi0FitRelErrorNeg[i] = - LargestDiffPi0FitErrorNeg[i]/SysErrCutPi0Fit[0][i].value*100.;
			LargestDiffPi0FitRelErrorPos[i] = LargestDiffPi0FitErrorPos[i]/SysErrCutPi0Fit[0][i].value*100.;
			if (i > 0) SysErrDatPi0Fit << BinsXCenterIncR[i] << "\t" << LargestDiffPi0FitNeg[i]/SysErrCutPi0Fit[0][i].value*100. << "\t" << LargestDiffPi0FitErrorNeg[i]/SysErrCutPi0Fit[0][i].value*100. << "\t" << LargestDiffPi0FitPos[i]/SysErrCutPi0Fit[0][i].value*100. << "\t" << LargestDiffPi0FitErrorPos[i]/SysErrCutPi0Fit[0][i].value*100.<<endl;
		} else {
			LargestDiffPi0FitRelNeg[i] = 0.;
			LargestDiffPi0FitRelPos[i] = 0.;
			LargestDiffPi0FitRelErrorNeg[i] = 0.;
			LargestDiffPi0FitRelErrorPos[i] = 0.;
		}
	}

	SysErrDatPi0Fit.close();
				
	TGraphAsymmErrors* SystErrGraphNegPi0Fit = new TGraphAsymmErrors(NBinsPtdoubleR+1, BinsXCenterIncR, LargestDiffPi0FitRelNeg, BinsXWidthIncR, BinsXWidthIncR, LargestDiffPi0FitRelErrorNeg, LargestDiffPi0FitRelErrorNeg);
	SystErrGraphNegPi0Fit->SetName(Form("Pi0Fit_SystErrorRelNeg_%s",cutVariationName.Data()));
	TGraphAsymmErrors* SystErrGraphPosPi0Fit = new TGraphAsymmErrors(NBinsPtdoubleR+1, BinsXCenterIncR, LargestDiffPi0FitRelPos, BinsXWidthIncR, BinsXWidthIncR, LargestDiffPi0FitRelErrorPos, LargestDiffPi0FitRelErrorPos);
	SystErrGraphPosPi0Fit->SetName(Form("Pi0Fit_SystErrorRelPos_%s",cutVariationName.Data()));
	const char* Outputname = Form("%s/Gamma_SystematicErrorCuts.root",outputFileDir.Data());
	TFile* SystematicErrorFile = new TFile(Outputname,"UPDATE");
	SystErrGraphPosPi0Fit->Write(Form("Pi0Fit_SystErrorRelPos_%s_%s",cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()),TObject::kOverwrite);
	SystErrGraphNegPi0Fit->Write(Form("Pi0Fit_SystErrorRelNeg_%s_%s",cutVariationName.Data(),(GetCentralityString(cutSelection[0])).Data()),TObject::kOverwrite);
	SystematicErrorFile->Write();
	SystematicErrorFile->Close();





}


void PlotCanvas(Int_t i, Int_t number, TCanvas *canvas, TH1D *spectrum, TLegend *legend, TCanvas *ratiocanvas, TH1D *ratio, TLegend *ratiolegend, TString cutSelection, TF1 *OneA = NULL, TF1 *OneB = NULL){
   

	if(spectrum->GetMinimum()<=0.0)
		spectrum->SetMinimum(spectrum->GetBinContent(spectrum->GetNbinsX())/10.);

	canvas->cd();
	spectrum->GetYaxis()->SetTitleSize(0.035);
	spectrum->GetXaxis()->SetTitleSize(0.035);
	spectrum->GetYaxis()->SetLabelSize(0.03);
	spectrum->GetXaxis()->SetLabelSize(0.03);
	spectrum->GetYaxis()->SetTitleOffset(1.);
	spectrum->GetXaxis()->SetTitleOffset(1.);
	if(!i){
		spectrum->DrawCopy();
		if(OneB)OneB->Draw("same");
	}
	else spectrum->DrawCopy("same");
	legend->AddEntry(spectrum,cutSelection,"p");
	if(i==number-1) legend->Draw();
	ratiocanvas->cd();
	//    ratiocanvas->SetGridx();
	//    ratiocanvas->SetGridy();
	if(!i){
		ratio->DrawCopy("e1][");
		OneA->Draw("same");
	}
	else ratio->DrawCopy("e1same][");
	ratiolegend->AddEntry(ratio,cutSelection,"p");
	if(i==number-1) ratiolegend->Draw();
   
}

