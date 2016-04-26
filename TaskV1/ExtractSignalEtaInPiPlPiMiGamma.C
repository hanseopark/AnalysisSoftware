// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion

//This file is not supposed to be run on outputfiles of the GammaConv-Software before the 30th Sept 2010.

#include <stdlib.h>
#include <iostream>
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
#include "TProfile2D.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "../CommonHeaders/PlottingMeson.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TMath.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "ExtractSignalEtaInPiPlPiMiGamma.h"
#include "../CommonHeaders/ExtractSignalBinning.h"
#include "../CommonHeaders/ExtractSignalPlotting.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "THnSparse.h"




// Main Function
void ExtractSignalEtaInPiPlPiMiGamma(TString file="", TString cutSelection="", TString Suffix="", TString optionMC="", TString optionEnergy="",  TString optionUseMinBiasEff="", TString optionPeriod="", TString optionAdvancedMesonQA="",Int_t numberOfBins = 30, Bool_t addSig = kFALSE) {
	gROOT->Reset();

	cout << "THIS macro is not up to date concerning the CUTSELECTION please don't run it in this state!" << endl;
	return;
	
	if(optionAdvancedMesonQA.Contains("AdvancedMesonQA")){fAdvancedMesonQA = kTRUE;}

	TString fGammaCutSelection;
	TString fPionCutSelection;
	TString fMesonCutSelection;
	fCutSelection = cutSelection;
	TString fCutSelectionRead = cutSelection;
	ReturnSeparatedCutNumber(cutSelection, fGammaCutSelection, fPionCutSelection,fMesonCutSelection,kTRUE);
	TString fGammaCutSelectionRead = fGammaCutSelection.Data();
	TString fMesonCutSelectionRead = fMesonCutSelection.Data();
	if (addSig) {
		cout << "running added Signal" << endl;
		cout << fGammaCutSelection.Data() << endl;
		fGammaCutSelection.Replace(GetEventRejectExtraSignalsCutPosition(),1,"2");
		cout << fGammaCutSelection.Data() << endl;
		fGammaCutSelectionRead = fGammaCutSelection;
		fMesonCutSelectionRead = fMesonCutSelection;
		fCutSelectionRead = Form("%s_%s", fGammaCutSelection.Data(), fMesonCutSelection.Data());
		cout << fCutSelectionRead.Data() << endl;
	}
	if(optionUseMinBiasEff.CompareTo("MinBiasEffOnly")==0 && optionMC.CompareTo("kTRUE") == 0){
		cout << "calculating MinBias Eff" << endl;
		cout << fGammaCutSelection.Data() << endl;
		fGammaCutSelection.Replace(GetEventCentralityMinCutPosition(),2,"00");
		fGammaCutSelectionRead = fGammaCutSelection;
		fMesonCutSelectionRead = fMesonCutSelection;      
		cout << fGammaCutSelection.Data() << endl;
		fCutSelectionRead = Form("%s_%s", fGammaCutSelection.Data(), fMesonCutSelection.Data());
		cout << fCutSelectionRead.Data() << endl;
	}
	
	StyleSettingsThesis();
	SetPlotStyle();
		
	fEnergyFlag = optionEnergy;
	fPrefix="Eta";
	fPeriodFlag = optionPeriod;

	TString outputDir = Form("%s/%s/%s/ExtractSignal",cutSelection.Data(),optionEnergy.Data(),Suffix.Data());
	gSystem->Exec("mkdir -p "+outputDir);
	
	cout << "Pictures are saved as " << Suffix.Data()<< endl;
	fdate = ReturnDateString();
	
	//****************************** Specification of collision system ************************************************
	TString textProcess = ReturnMesonString (fPrefix);
	if(textProcess.CompareTo("") == 0 ){
		cout << "Meson unknown" << endl;
		return ;
	}
	
	fTextMeasurement = ReturnFullTextMesonPiPlPiMiGamma(fEnergyFlag, textProcess);
	fCollisionSystem = ReturnFullCollisionsSystem(fEnergyFlag);
	if (fCollisionSystem.CompareTo("") == 0){
		cout << "No correct collision system specification, has been given" << endl;
		return;
	}

	//****************************** Choice of Fitting procedure ******************************************************
	cout << "Gaussian fit chosen ..." << endl;
	
	if(cutSelection.Length() == 0){
		cout<<"ERROR: Cut selection is not set, please do!"<<endl;
		return;
	}

	//***************************** Specification Data/MC ************************************************************
	if(optionMC.CompareTo("kTRUE") == 0){
		fIsMC = 1;
		fPrefix2 = "MC";
	} else {
		fIsMC = 0;
		fPrefix2 = "data";
	}
	
	//**************************** Determine Centrality *************************************************************
	centralityString = GetCentralityString(fGammaCutSelection);
	if (centralityString.CompareTo("pp")==0){
		fTextCent = "MinBias";  
	} else {
		fTextCent = Form("%s central", centralityString.Data());
	}
	
	//***************************** Initialization of variables according to meson type ******************************
	Initialize("Eta",numberOfBins);

	
	//************************* Start of Main routine ***************************************************************
	const char* fFileErrLogDatname = Form("%s/%s/%s_%s_FileErrLog%s_%s.dat",cutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelectionRead.Data());
	fFileErrLog.open(fFileErrLogDatname, ios::out);

	TFile f(file.Data());
	
	TList *TopDir =(TList*)f.Get("GammaConvEtaPiPlPiMiGamma_1");
	if(TopDir == NULL){
		cout<<"ERROR: TopDir not Found"<<endl;
		return;
	}
	TList *HistosGammaConversion = (TList*)TopDir->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
	if(HistosGammaConversion == NULL){
		cout<<"ERROR: " << Form("Cut Number %s",fCutSelectionRead.Data()) << " not Found in File"<<endl;
		return;
	}
	TList *ESDContainer = (TList*) HistosGammaConversion->FindObject(Form("%s ESD histograms",fCutSelectionRead.Data()));
	TList *BackgroundContainer = (TList*) HistosGammaConversion->FindObject(Form("%s Back histograms",fCutSelectionRead.Data()));
	TList *MotherContainer = (TList*) HistosGammaConversion->FindObject(Form("%s Mother histograms",fCutSelectionRead.Data()));
	cout << fMesonCutSelectionRead.Data() << endl;
	cout << fGammaCutSelectionRead.Data() << endl;   
	fNumberOfGoodESDTracks = (TH1D*)ESDContainer->FindObject("GoodESDTracks");
	fEventQuality = (TH1D*)ESDContainer->FindObject("NEvents");
		
	TString rapidityRange;
	fYMaxMeson =  ReturnRapidityStringAndDouble(fMesonCutSelection, rapidityRange);
	fBackgroundMultNumber = ReturnBackgroundMult(fMesonCutSelection);
		
	TString ObjectNameESD		= "ESD_Mother_InvMass_Pt";
	TString ObjectNameBck		= "ESD_Background_InvMass_Pt";
	
	SetCorrectMCHistogrammNames(optionPeriod.Data(), optionEnergy.Data());
	
	fGammaGammaInvMassVSPt= (TH2D*)ESDContainer->FindObject(ObjectNameESD.Data());
	fBckInvMassVSPt = (TH2D*)ESDContainer->FindObject(ObjectNameBck.Data());

	const char* FileDataLogname = Form("%s/%s/%s_%s_EffiCheck_RAWDATA%s_%s.dat",cutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelectionRead.Data());
	fFileDataLog.open(FileDataLogname, ios::out);

	ProduceBckProperWeighting(BackgroundContainer,MotherContainer);
	
	
	if(fIsMC){
		TList *MCContainer = (TList*)HistosGammaConversion->FindObject(Form("%s MC histograms",fCutSelectionRead.Data()));
		TList *TrueConversionContainer = (TList*)HistosGammaConversion->FindObject(Form("%s True histograms",fCutSelectionRead.Data()));

		
		if( fMesonId == 221){
			fHistoMCMesonPtWithinAcceptance = (TH1D*)MCContainer->FindObject(ObjectNameMCEtaAcc.Data());
			fHistoMCMesonPt = (TH1D*)MCContainer->FindObject(ObjectNameMCEta.Data());   // Not the best; better having a 2D Pt_vs_Rapid in case we change limits
			fHistoMCMesonPtWOWeights =(TH1D*)MCContainer->FindObject(ObjectNameMCEtaWOWeights.Data());
		}
		fHistoMCMesonPt->Sumw2();
		fHistoMCMesonPtWithinAcceptance->Sumw2();
		if (fHistoMCMesonPtWOWeights){
		fHistoMCMesonPtWeights = (TH1D*)fHistoMCMesonPtWOWeights->Clone("WeightsMeson");
		fHistoMCMesonPtWeights->Divide(fHistoMCMesonPt,fHistoMCMesonPtWOWeights, 1.,1.,"B");
			
		}   
		fHistoTrueMesonInvMassVSPt = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrue.Data());
		cout << "here" << endl;
		FillMassMCTrueMesonHistosArray(fHistoTrueMesonInvMassVSPt);
		cout << "here" << endl;
		fHistoTrueGGMesonInvMassVSPt = (TH2D*)TrueConversionContainer->FindObject(ObjectNameContaminationGG.Data());
		FillMassMCTrueGGMesonHistosArray(fHistoTrueGGMesonInvMassVSPt);
		fHistoTrueDalitzMesonInvMassVSPt = (TH2D*)TrueConversionContainer->FindObject(ObjectNameContaminationDalitz.Data());
		FillMassMCTrueDalitzMesonHistosArray(fHistoTrueDalitzMesonInvMassVSPt);   
		
	}
	
	fMesonMassExpect = TDatabasePDG::Instance()->GetParticle(fMesonId)->Mass();
	if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 || fEnergyFlag.CompareTo("pPb_5.023TeV") == 0){
		fNEvents = fEventQuality->GetBinContent(1);
	} else {
		fNEvents =  GetNEvents(fEventQuality);
	}
	
	TH1D *fBck = (TH1D*)fBckInvMassVSPt->ProjectionX("ESD_Background_InvMass");
	TH1D *fGammaGamma = (TH1D*)fGammaGammaInvMassVSPt->ProjectionX("ESD_Mother_InvMass");

	cout<< "The mass of the meson is: "<< fMesonMassExpect<< " Events analysed: "<< fNEvents<< endl;
	
	// Process the 1D invariant mass histos
	fGammaGamma->SetTitle(Form("%s %s",fGammaGamma->GetTitle(),fCutSelection.Data()));
	fGammaGamma->Rebin(fNRebinGlobal);
	//fGammaGamma->Scale(1./fNRebinGlobal);
	fBck->Rebin(fNRebinGlobal);
	//fBck->Scale(1./fNRebinGlobal);
	ProcessEM( fGammaGamma , fBck, fBGFitRange);
	fHistoMappingBackNormInvMass = fBckNorm;
	fHistoMappingSignalInvMass = fSignal;
	
	fGammaGamma->DrawCopy();
	fHistoMappingBackNormInvMass->DrawCopy("same");
	fHistoMappingSignalInvMass->DrawCopy("same");

	// Function to Project the 2D histos InvariantMass VS Pt into Invariant Mass spectrum
	FillMassHistosArray(fGammaGammaInvMassVSPt);
	
	
	ProcessEM( fMesonFullPtSignal, fMesonFullPtBackground, fBGFitRange);
	fMesonFullPtBackNorm = fBckNorm;
	
	ProcessEM( fFittingHistMidPtSignal, fFittingHistMidPtBackground, fBGFitRange);
	fFittingHistMidPtSignalSub = fSignal;
	fFileErrLog << "Using exp fit"<<endl;
	FitSubtractedInvMassInPtBins(fFittingHistMidPtSignalSub, fMesonIntRange,200,kTRUE);
	fFitSignalInvMassMidPt = fFitReco;
	
	
	TString fDecayChannel = "#pi^{+}#pi^{-}#gamma";

	delete fMidPt;
	
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){ // BEGIN ANALYSIS for each Pt bin
		//       cout << "Begin Analysis Pt Bin " << iPt <<endl;
		// Function to subtract GG minus Bck
		ProcessBckFitSubtraction(fHistoMappingPiPlPiMiGammaInvMassPtBin[iPt],iPt,fPeakRange,fFitRange);
		
		ProcessEM( fHistoMappingPiPlPiMiGammaInvMassPtBin[iPt], fHistoMappingBackInvMassPtBin[iPt], fBGFitRange);
		fHistoMappingSignalInvMassPtBin[iPt] = fSignal;
		fHistoMappingBackNormInvMassPtBin[iPt] = fBckNorm;

		FitWithPol2ForBG(fHistoMappingPiPlPiMiGammaInvMassPtBin[iPt], fMesonIntRange,iPt,kFALSE);
		fFitWithPol2ForBG[iPt] = fFitReco;
		
		//       cout<< "iPt"<< iPt<< " "<< "standard range"<<endl;
		// Fitting the subtracted spectra
		fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "normal range/right normalization" << endl;

		fFitSignalInvMassPtBin[iPt]=0x00;
		fFitSignalInvMassBackFitPtBin[iPt]=0x00;
		fFileErrLog << "Using exp fit"<<endl;
		FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntRange,iPt,kTRUE);
		fFitSignalInvMassPtBin[iPt] = fFitReco;
		fFitSignalPeakPosInvMassPtBin[iPt] = fFitGausExp;
		fFitBckInvMassPtBin[iPt] = fFitLinearBck;
		fMesonYieldsResidualBckFunc[iPt] = fIntLinearBck;
		fMesonYieldsResidualBckFuncError[iPt] = fIntLinearBckError;
		FitSubtractedInvMassInPtBins(fHistoMappingGGInvMassBackFitPtBin[iPt], fMesonIntRange,iPt,kTRUE);
		fFitSignalInvMassBackFitPtBin[iPt] = fFitReco;
		fFitSignalPeakPosInvMassBackFitPtBin[iPt] = fFitGausExp;
		fFitBckInvMassBackFitPtBin[iPt] = fFitLinearBck;
		fMesonYieldsResidualBckFuncBackFit[iPt] = fIntLinearBck;
		fMesonYieldsResidualBckFuncBackFitError[iPt] = fIntLinearBckError;
				//GetFWHM
		CalculateFWHM( fFitSignalInvMassPtBin[iPt]);
		fMesonFWHM[iPt] = fFWHMFunc;
		fMesonFWHMError[iPt] = fFWHMFuncError;

		
		if (fFitSignalInvMassPtBin[iPt] !=0x00){
			fMesonMass[iPt] = fFitSignalInvMassPtBin[iPt]->GetParameter(1);
			fMesonMassError[iPt] = fFitSignalInvMassPtBin[iPt]->GetParError(1);
			fMesonWidth[iPt] = fFitSignalInvMassPtBin[iPt]->GetParameter(2);
			fMesonWidthError[iPt] = fFitSignalInvMassPtBin[iPt]->GetParError(2);
			fMesonCurIntRange[0] = fMesonIntRange[0] - (fMesonMassExpect-fMesonMass[iPt]);
			fMesonCurIntRangeWide[0] = fMesonIntRangeWide[0] - (fMesonMassExpect-fMesonMass[iPt]);
			fMesonCurIntRangeNarrow[0] = fMesonIntRangeNarrow[0] - (fMesonMassExpect-fMesonMass[iPt]);
			fMesonCurIntRange[1] = fMesonIntRange[1] - (fMesonMassExpect-fMesonMass[iPt]);
			fMesonCurIntRangeWide[1] = fMesonIntRangeWide[1] - (fMesonMassExpect-fMesonMass[iPt]);
			fMesonCurIntRangeNarrow[1] = fMesonIntRangeNarrow[1] - (fMesonMassExpect-fMesonMass[iPt]);
			
		} else {
			fMesonMass[iPt] = 0.;
			fMesonMassError[iPt] = 0.;
			fMesonWidth[iPt] = 0.;
			fMesonWidthError[iPt] = 0.;
			fMesonCurIntRange[0] = fMesonIntRange[0];
			fMesonCurIntRangeWide[0] = fMesonIntRangeWide[0];
			fMesonCurIntRangeNarrow[0] = fMesonIntRangeNarrow[0];
			fMesonCurIntRange[1] = fMesonIntRange[1];
			fMesonCurIntRangeWide[1] = fMesonIntRangeWide[1];
			fMesonCurIntRangeNarrow[1] = fMesonIntRangeNarrow[1];
		}

		if (fFitSignalInvMassBackFitPtBin[iPt] !=0x00){
			cout<<fFitSignalInvMassBackFitPtBin[iPt]->GetParameter(1)<<endl;
			fMesonMassBackFit[iPt] = fFitSignalInvMassBackFitPtBin[iPt]->GetParameter(1);
			fMesonMassBackFitError[iPt] = fFitSignalInvMassBackFitPtBin[iPt]->GetParError(1);
			fMesonWidthBackFit[iPt] = fFitSignalInvMassBackFitPtBin[iPt]->GetParameter(2);
			fMesonWidthBackFitError[iPt] = fFitSignalInvMassBackFitPtBin[iPt]->GetParError(2);
			fMesonCurIntRangeBackFit[0] = fMesonIntRange[0] - (fMesonMassExpect-fMesonMassBackFit[iPt]);
			fMesonCurIntRangeBackFit[1] = fMesonIntRange[1] - (fMesonMassExpect-fMesonMassBackFit[iPt]);
		}
		else{
			fMesonMassBackFit[iPt] = 0;
			fMesonMassBackFitError[iPt] = 0;
			fMesonWidthBackFit[iPt] = 0;
			fMesonWidthBackFitError[iPt] = 0;
			fMesonCurIntRangeBackFit[0] = fMesonIntRange[0];
			fMesonCurIntRangeBackFit[1] = fMesonIntRange[1];
		}


		IntegrateHistoInvMass( fHistoMappingPiPlPiMiGammaInvMassPtBin[iPt], fMesonCurIntRange);
		fPiPlPiMiGammaYields[iPt] = fYields;
		fPiPlPiMiGammaYieldsError[iPt] = fYieldsError;
		IntegrateHistoInvMass( fHistoMappingPiPlPiMiGammaInvMassPtBin[iPt], fMesonCurIntRangeWide);
		fPiPlPiMiGammaYieldsWide[iPt] = fYields;
		fPiPlPiMiGammaYieldsWideError[iPt] = fYieldsError;
		IntegrateHistoInvMass( fHistoMappingPiPlPiMiGammaInvMassPtBin[iPt], fMesonCurIntRangeNarrow);
		fPiPlPiMiGammaYieldsNarrow[iPt] = fYields;
		fPiPlPiMiGammaYieldsNarrowError[iPt] = fYieldsError;

		// Integrate the bck histo
		IntegrateHistoInvMass( fHistoMappingBackNormInvMassPtBin[iPt], fMesonCurIntRange);
		fBckYields[iPt] = fYields;
		fBckYieldsError[iPt] = fYieldsError;
		IntegrateHistoInvMass( fHistoMappingBackNormInvMassPtBin[iPt], fMesonCurIntRangeWide);
		fBckYieldsWide[iPt] = fYields;
		fBckYieldsWideError[iPt] = fYieldsError;
		IntegrateHistoInvMass( fHistoMappingBackNormInvMassPtBin[iPt], fMesonCurIntRangeNarrow);
		fBckYieldsNarrow[iPt] = fYields;
		fBckYieldsNarrowError[iPt] = fYieldsError;

		// Integrate the signal histo
		fFileDataLog<< endl <<"Signal histo normal range, right norm "<< fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
		IntegrateHistoInvMassStream( fHistoMappingGGInvMassBackFitPtBin[iPt], fMesonCurIntRangeBackFit);
		fMesonYieldsBackFit[iPt] = fYields;
		fMesonYieldsBackFitError[iPt] = fYieldsError;
		IntegrateHistoInvMassStream( fHistoMappingSignalInvMassPtBin[iPt], fMesonCurIntRange);
		fMesonYields[iPt] = fYields;
		fMesonYieldsError[iPt] = fYieldsError;
		fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
		fFileDataLog<< endl <<"Signal histo wide range, right norm " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
		IntegrateHistoInvMassStream( fHistoMappingSignalInvMassPtBin[iPt], fMesonCurIntRangeWide);
		fMesonYieldsWide[iPt] = fYields;
		fMesonYieldsWideError[iPt] = fYieldsError;
		fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
		fFileDataLog<< endl <<"Signal histo narrow range, right norm" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
		IntegrateHistoInvMassStream( fHistoMappingSignalInvMassPtBin[iPt], fMesonCurIntRangeNarrow);
		fMesonYieldsNarrow[iPt] = fYields;
		fMesonYieldsNarrowError[iPt] = fYieldsError;
		fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
		if(fIsMC){
			fFileDataLog<< endl <<"True histo normal range" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			fFitTrueSignalInvMassPtBin[iPt]=0x00;		
			fFileErrLog << "Using exp fit"<<endl;
			FitTrueInvMassInPtBins(fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonIntRange,iPt);
			
			//	FitSubtractedInvMassInPtBins(fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonIntRange,iPt,kFALSE);
			if (fHistoMappingTrueMesonInvMassPtBins[iPt]->GetEntries() !=0){
				fFitTrueSignalInvMassPtBin[iPt] = fFitReco;
				if (fFitTrueSignalInvMassPtBin[iPt] != 0x00){
					fMesonTrueMass[iPt] = fFitTrueSignalInvMassPtBin[iPt]->GetParameter(1);
					fMesonTrueMassError[iPt] = fFitTrueSignalInvMassPtBin[iPt]->GetParError(1);
					CalculateFWHM(fFitTrueSignalInvMassPtBin[iPt]);
					fMesonTrueFWHM[iPt] = fFWHMFunc;
					fMesonTrueFWHMError[iPt] = fFWHMFuncError;
					fFileDataLog << "TrueFWHM \t" << fMesonTrueFWHM[iPt] << "\t +-" << fMesonTrueFWHMError[iPt] << endl;
					fMesonTrueIntRange[0] = fMesonIntRange[0] - (fMesonMassExpect-fMesonTrueMass[iPt]);
					fMesonTrueIntRangeWide[0] = fMesonIntRangeWide[0] - (fMesonMassExpect-fMesonTrueMass[iPt]);
					fMesonTrueIntRangeNarrow[0] = fMesonIntRangeNarrow[0] - (fMesonMassExpect-fMesonTrueMass[iPt]);
					fMesonTrueIntRange[1] = fMesonIntRange[1] - (fMesonMassExpect-fMesonTrueMass[iPt]);
					fMesonTrueIntRangeWide[1] = fMesonIntRangeWide[1] - (fMesonMassExpect-fMesonTrueMass[iPt]);
					fMesonTrueIntRangeNarrow[1] = fMesonIntRangeNarrow[1] - (fMesonMassExpect-fMesonTrueMass[iPt]);
				} else {
					fMesonTrueMass[iPt] = 0.;
					fMesonTrueMassError[iPt] = 1.;
					fMesonTrueFWHM[iPt] = 0.;
					fMesonTrueFWHMError[iPt] = 0.;
					fMesonTrueIntRange[0] = fMesonIntRange[0];
					fMesonTrueIntRangeWide[0] = fMesonIntRangeWide[0];
					fMesonTrueIntRangeNarrow[0] = fMesonIntRangeNarrow[0];
					fMesonTrueIntRange[1] = fMesonIntRange[1];
					fMesonTrueIntRangeWide[1] = fMesonIntRangeWide[1];
					fMesonTrueIntRangeNarrow[1] = fMesonIntRangeNarrow[1];

				}

			}

			cout << "here" << endl;
			IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonTrueIntRange);
			fMesonTrueYields[iPt] = fYields;
			fMesonTrueYieldsError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

			fFileDataLog<< endl <<"True histo wide range" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonTrueIntRangeWide);
			fMesonTrueYieldsWide[iPt] = fYields;
			fMesonTrueYieldsWideError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

			fFileDataLog<< endl <<"True histo narrow range" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonTrueIntRangeNarrow);
			fMesonTrueYieldsNarrow[iPt] = fYields;
			fMesonTrueYieldsNarrowError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

			cout<< "Analyse reweighted histos" << endl;
					
			fFileDataLog<< endl <<"TrueGG cont histo" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			IntegrateHistoInvMassStream( fHistoMappingTrueGGMesonInvMassPtBins[iPt], fMesonTrueIntRange);
			fMesonTrueGGContYields[iPt] = fYields;
			fMesonTrueGGContYieldsError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;
				
			fFileDataLog<< endl <<"TrueGG cont histo wide range  " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			IntegrateHistoInvMassStream( fHistoMappingTrueGGMesonInvMassPtBins[iPt], fMesonTrueIntRangeWide);
			fMesonTrueGGContYieldsWide[iPt] = fYields;
			fMesonTrueGGContYieldsWideError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;
				
			fFileDataLog<< endl <<"TrueGG cont histo narrow range " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			IntegrateHistoInvMassStream( fHistoMappingTrueGGMesonInvMassPtBins[iPt], fMesonTrueIntRangeNarrow);
			fMesonTrueGGContYieldsNarrow[iPt] = fYields;
			fMesonTrueGGContYieldsNarrowError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

			fFileDataLog<< endl <<"TrueDalitz cont histo " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			IntegrateHistoInvMassStream( fHistoMappingTrueDalitzInvMassPtBins[iPt], fMesonTrueIntRange);
			fMesonTrueDalitzContYields[iPt] = fYields;
			fMesonTrueDalitzContYieldsError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

			fFileDataLog<< endl <<"TrueDalitz cont histo wide range " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			IntegrateHistoInvMassStream( fHistoMappingTrueDalitzInvMassPtBins[iPt], fMesonTrueIntRangeWide);
			fMesonTrueDalitzContYieldsWide[iPt] = fYields;
			fMesonTrueDalitzContYieldsWideError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

			fFileDataLog<< endl <<"TrueDalitz cont histo narrow range " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			IntegrateHistoInvMassStream( fHistoMappingTrueDalitzInvMassPtBins[iPt], fMesonTrueIntRangeNarrow);
			fMesonTrueDalitzContYieldsNarrow[iPt] = fYields;
			fMesonTrueDalitzContYieldsNarrowError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

			if( (fPiPlPiMiGammaYields[iPt] - fMesonTrueYields[iPt]) > 0) {
				fMesonTrueSB[iPt]   = fMesonTrueYields[iPt] / ( fPiPlPiMiGammaYields[iPt] - fMesonTrueYields[iPt] );
				fMesonTrueSign[iPt] = fMesonTrueYields[iPt] / pow( ( fPiPlPiMiGammaYields[iPt] - fMesonTrueYields[iPt] ) , 0.5);
				fMesonTrueSBError[iPt] = 0;
				fMesonTrueSignError[iPt] = 0;
			}
			else {
				fMesonTrueSB[iPt] = 0.;
				fMesonTrueSign[iPt] = 0.;
				fMesonTrueSBError[iPt] = 0.;
				fMesonTrueSignError[iPt] = 0.;
			}
		}
		fFileDataLog<< "Residual Background leftover norm integration/right norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFunc[iPt] << "\t +- \t" << fMesonYieldsResidualBckFuncError[iPt] << endl<< endl;
		fTotalBckYields[iPt] = fBckYields[iPt] + fMesonYieldsResidualBckFunc[iPt];
		fTotalBckYieldsError[iPt] = pow(fBckYieldsError[iPt]*fBckYieldsError[iPt] + fMesonYieldsResidualBckFuncError[iPt]*fMesonYieldsResidualBckFuncError[iPt],0.5);

		fMesonYieldsCorResidualBckFunc[iPt] = fMesonYields[iPt]- fMesonYieldsResidualBckFunc[iPt];
		fMesonYieldsCorResidualBckFuncError[iPt] =
			pow((fMesonYieldsError[iPt]*fMesonYieldsError[iPt]+
				fMesonYieldsResidualBckFuncError[iPt]*fMesonYieldsResidualBckFuncError[iPt]),0.5);
		fMesonYieldsPerEvent[iPt]= fMesonYieldsCorResidualBckFunc[iPt]/fNEvents;
		fMesonYieldsPerEventError[iPt]= fMesonYieldsCorResidualBckFuncError[iPt]/fNEvents;

		fMesonYieldsCorResidualBckFuncBackFit[iPt] = fMesonYieldsBackFit[iPt]- fMesonYieldsResidualBckFuncBackFit[iPt];
		fMesonYieldsCorResidualBckFuncBackFitError[iPt] =
			pow((fMesonYieldsBackFitError[iPt]*fMesonYieldsBackFitError[iPt]+
				fMesonYieldsResidualBckFuncBackFitError[iPt]*fMesonYieldsResidualBckFuncBackFitError[iPt]),0.5);
		fMesonYieldsPerEventBackFit[iPt]= fMesonYieldsCorResidualBckFuncBackFit[iPt]/fNEvents;
		fMesonYieldsPerEventBackFitError[iPt]= fMesonYieldsCorResidualBckFuncBackFitError[iPt]/fNEvents;


		//Integrate Fit Function
		IntegrateFitFunc( fFitSignalPeakPosInvMassPtBin[iPt], fHistoMappingSignalInvMassPtBin[iPt], fMesonCurIntRange);
		fMesonYieldsFunc[iPt]=fYieldsFunc;


			
		if( fFitBckInvMassPtBin[iPt]->Integral(fMesonMass[iPt]-fMesonFWHM[iPt], fMesonMass[iPt]+fFitSignalInvMassPtBin[iPt]->GetParameter(2))!=0){
			Double_t background = fFitBckInvMassPtBin[iPt]->Integral(fMesonMass[iPt]-fMesonFWHM[iPt], 
																	fMesonMass[iPt]+fFitSignalInvMassPtBin[iPt]->GetParameter(2));
			Double_t backgroundErr = fFitBckInvMassPtBin[iPt]->IntegralError(fMesonMass[iPt]-fMesonFWHM[iPt], 
																			fMesonMass[iPt]+fFitSignalInvMassPtBin[iPt]->GetParameter(2));
			Double_t signal = fFitSignalInvMassPtBin[iPt]->Integral(fMesonMass[iPt]-fMesonFWHM[iPt], 
																	fMesonMass[iPt]+fFitSignalInvMassPtBin[iPt]->GetParameter(2)) - background;
			Double_t signalErr =pow( pow(fFitSignalInvMassPtBin[iPt]->IntegralError(fMesonMass[iPt]-fMesonFWHM[iPt],
																					fMesonMass[iPt]+fFitSignalInvMassPtBin[iPt]->GetParameter(2)),2 )+ pow(backgroundErr,2),0.5);
			fMesonSB[iPt] = signal/ background;
			fMesonSBError[iPt] = pow( pow(signalErr/background,2.)+pow(signal/(background *background )*backgroundErr ,2.) ,0.5); 
			fMesonSign[iPt] = signal/ pow(background + signal,0.5);
			fMesonSignError[iPt] = pow(pow( (pow(background + signal ,0.5) - 0.5*signal*pow(background+signal,-0.5))/(background+signal) * signalErr ,2) + pow( 0.5*pow(signal+background,-1.5),2) ,0.5); 
		}else{
			fMesonSB[iPt] = 0.;
			fMesonSBError[iPt] = 0.;
			fMesonSign[iPt] = 0.;
			fMesonSignError[iPt] = 0.;
		}

		//       cout<< "iPt"<< iPt<< " "<< "FWHM done"<<endl;

		// Wide integration mass window
		//       cout<< "iPt"<< iPt<< " "<< "wide range"<<endl;
		fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "wide range/right normalization" << endl;

		fFileErrLog << "Using exp fit"<<endl;
		FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntRangeWide,iPt,kTRUE);
		fMesonYieldsResidualBckFuncWide[iPt] = fIntLinearBck;
		fMesonYieldsResidualBckFuncWideError[iPt] = fIntLinearBckError;

		//FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt],fMesonIntRangeWide,iPt,kFALSE);
		fFileDataLog<< "Residual Background leftover wide integration/right norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFuncWide[iPt] <<"\t +- \t" << fMesonYieldsResidualBckFuncWideError[iPt] <<endl<< endl;
		fMesonYieldsCorResidualBckFuncWide[iPt] = fMesonYieldsWide[iPt]- fMesonYieldsResidualBckFuncWide[iPt];

		fTotalBckYieldsWide[iPt] = fBckYieldsWide[iPt] + fMesonYieldsResidualBckFuncWide[iPt];
		fTotalBckYieldsWideError[iPt] = pow(fBckYieldsWideError[iPt]*fBckYieldsWideError[iPt] + fMesonYieldsResidualBckFuncWideError[iPt]*fMesonYieldsResidualBckFuncWideError[iPt],0.5);

		fMesonYieldsCorResidualBckFuncWideError[iPt] =
			pow((fMesonYieldsWideError[iPt]*fMesonYieldsWideError[iPt]+
				fMesonYieldsResidualBckFuncWideError[iPt]*fMesonYieldsResidualBckFuncWideError[iPt]),0.5);
		fMesonYieldsPerEventWide[iPt]= fMesonYieldsCorResidualBckFuncWide[iPt]/fNEvents;
		fMesonYieldsPerEventWideError[iPt]= fMesonYieldsCorResidualBckFuncWideError[iPt]/fNEvents;

		if( fTotalBckYieldsWide[iPt]!=0){
			fMesonSBWide[iPt] = fMesonYieldsCorResidualBckFuncWide[iPt]/fTotalBckYieldsWide[iPt];
			fMesonSBWideError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncWideError[iPt]/fTotalBckYieldsWide[iPt],2.)+pow(fMesonYieldsCorResidualBckFuncWide[iPt]/(fTotalBckYieldsWide[iPt]*fTotalBckYieldsWide[iPt])*fTotalBckYieldsWideError[iPt],2.) ,0.5);
			fMesonSignWide[iPt] = fMesonYieldsCorResidualBckFuncWide[iPt]/pow(fTotalBckYieldsWide[iPt],0.5);
			fMesonSignWideError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncWideError[iPt]/pow(fTotalBckYieldsWide[iPt],0.5),2.)+pow(0.5*fMesonYieldsCorResidualBckFuncWide[iPt]/pow(fTotalBckYieldsWide[iPt],1.5)*fTotalBckYieldsWideError[iPt],2.) ,0.5);

		}else{
			fMesonSBWide[iPt] = 0.;
			fMesonSBWideError[iPt] = 0.;
			fMesonSignWide[iPt] = 0.;
			fMesonSignWideError[iPt] = 0.;
		}


		// Narrow integration mass window
		fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "narrow range/right normalization" << endl; ;
		//       cout<< "iPt"<< iPt<< " "<< "narrow range"<<endl;

		fFileErrLog << "Using exp fit"<<endl;
		FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntRangeNarrow,iPt,kTRUE);
		fMesonYieldsResidualBckFuncNarrow[iPt] = fIntLinearBck;
		fMesonYieldsResidualBckFuncNarrowError[iPt] = fIntLinearBckError;

		//FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt],fMesonIntRangeNarrow,iPt,kFALSE);
		fFileDataLog<< "Residual Background leftover narrow integration/right norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFuncNarrow[iPt] <<"\t +- \t" << fMesonYieldsResidualBckFuncNarrowError[iPt]<<endl<< endl;
		fMesonYieldsCorResidualBckFuncNarrow[iPt] = fMesonYieldsNarrow[iPt]- fMesonYieldsResidualBckFuncNarrow[iPt];
		fMesonYieldsCorResidualBckFuncNarrowError[iPt] =
			pow((fMesonYieldsNarrowError[iPt]*fMesonYieldsNarrowError[iPt]+
				fMesonYieldsResidualBckFuncNarrowError[iPt]*fMesonYieldsResidualBckFuncNarrowError[iPt]),0.5);
		fMesonYieldsPerEventNarrow[iPt]= fMesonYieldsCorResidualBckFuncNarrow[iPt]/fNEvents;
		fMesonYieldsPerEventNarrowError[iPt]= fMesonYieldsCorResidualBckFuncNarrowError[iPt]/fNEvents;

		fTotalBckYieldsNarrow[iPt] = fBckYieldsNarrow[iPt] + fMesonYieldsResidualBckFuncNarrow[iPt];
		fTotalBckYieldsNarrowError[iPt] = pow(fBckYieldsNarrowError[iPt]*fBckYieldsNarrowError[iPt] + fMesonYieldsResidualBckFuncNarrowError[iPt]*fMesonYieldsResidualBckFuncNarrowError[iPt],0.5);


		if( fTotalBckYieldsNarrow[iPt]!=0){
			fMesonSBNarrow[iPt] = fMesonYieldsCorResidualBckFuncNarrow[iPt]/fTotalBckYieldsNarrow[iPt];
			fMesonSBNarrowError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncNarrowError[iPt]/fTotalBckYieldsNarrow[iPt],2.)+pow(fMesonYieldsCorResidualBckFuncNarrow[iPt]/(fTotalBckYieldsNarrow[iPt]*fTotalBckYieldsNarrow[iPt])*fTotalBckYieldsNarrowError[iPt],2.) ,0.5);
			fMesonSignNarrow[iPt] = fMesonYieldsCorResidualBckFuncNarrow[iPt]/pow(fTotalBckYieldsNarrow[iPt],0.5);
			fMesonSignNarrowError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncNarrowError[iPt]/pow(fTotalBckYieldsNarrow[iPt],0.5),2.)+pow(0.5*fMesonYieldsCorResidualBckFuncNarrow[iPt]/pow(fTotalBckYieldsNarrow[iPt],1.5)*fTotalBckYieldsNarrowError[iPt],2.) ,0.5);

		}else{
			fMesonSBNarrow[iPt] = 0.;
			fMesonSBNarrowError[iPt] = 0.;
			fMesonSignNarrow[iPt] = 0.;
			fMesonSignNarrowError[iPt] = 0.;
		}

		//////////////////////////////// Start Analysis with  Normalization at the left of the Meson Peak

		// Function to subtract GG minus Bck
		cout << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
		ProcessEM( fHistoMappingPiPlPiMiGammaInvMassPtBin[iPt], fHistoMappingBackInvMassPtBin[iPt], fBGFitRangeLeft);
		fHistoMappingSignalInvMassLeftPtBin[iPt] = fSignal;
		fHistoMappingBackNormInvMassLeftPtBin[iPt] = fBckNorm;


		//       cout<< "iPt"<< iPt<< " "<< "standard range"<<endl;
		// Fitting the subtracted spectra
		fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "normal range/left normalization" << endl;

		fFitInvMassLeftPtBin[iPt] =0x00;
		fFileErrLog << "Using exp fit"<<endl;
		FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntRange,iPt,kTRUE);
		fFitInvMassLeftPtBin[iPt] = fFitReco;
		fFitSignalPeakPosInvMassLeftPtBin[iPt] = fFitGausExp;
		fFitBckInvMassLeftPtBin[iPt] = fFitLinearBck;
		fMesonYieldsResidualBckFuncLeft[iPt] = fIntLinearBck;
		fMesonYieldsResidualBckFuncLeftError[iPt] = fIntLinearBckError;
		//FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntRange,iPt,kFALSE);
		CalculateFWHM(fFitInvMassLeftPtBin[iPt]);
		fMesonFWHMLeft[iPt] = fFWHMFunc;
		fMesonFWHMLeftError[iPt] = fFWHMFuncError;
		
		if (fFitInvMassLeftPtBin[iPt] !=0x00){
			fMesonMassLeft[iPt] = fFitInvMassLeftPtBin[iPt]->GetParameter(1);
			fMesonMassLeftError[iPt] = fFitInvMassLeftPtBin[iPt]->GetParError(1);
			fMesonWidthLeft[iPt] = fFitInvMassLeftPtBin[iPt]->GetParameter(2);
			fMesonWidthLeftError[iPt] = fFitInvMassLeftPtBin[iPt]->GetParError(2);
			fMesonCurLeftIntRange[0] = fMesonIntRange[0] - (fMesonMassExpect-fMesonMassLeft[iPt]);
			fMesonCurLeftIntRangeWide[0] = fMesonIntRangeWide[0] - (fMesonMassExpect-fMesonMassLeft[iPt]);
			fMesonCurLeftIntRangeNarrow[0] = fMesonIntRangeNarrow[0] - (fMesonMassExpect-fMesonMassLeft[iPt]);
			fMesonCurLeftIntRange[1] = fMesonIntRange[1] - (fMesonMassExpect-fMesonMassLeft[iPt]);
			fMesonCurLeftIntRangeWide[1] = fMesonIntRangeWide[1] - (fMesonMassExpect-fMesonMassLeft[iPt]);
			fMesonCurLeftIntRangeNarrow[1] = fMesonIntRangeNarrow[1] - (fMesonMassExpect-fMesonMassLeft[iPt]);   
		} else {
			fMesonMassLeft[iPt] = 0.;
			fMesonMassLeftError[iPt] = 0.;
			fMesonWidthLeft[iPt] = 0.;
			fMesonWidthLeftError[iPt] = 0.;
			fMesonCurLeftIntRange[0] = fMesonIntRange[0];
			fMesonCurLeftIntRangeWide[0] = fMesonIntRangeWide[0];
			fMesonCurLeftIntRangeNarrow[0] = fMesonIntRangeNarrow[0];
			fMesonCurLeftIntRange[1] = fMesonIntRange[1];
			fMesonCurLeftIntRangeWide[1] = fMesonIntRangeWide[1];
			fMesonCurLeftIntRangeNarrow[1] = fMesonIntRangeNarrow[1];
		}

		// Integrate the bck histo
		IntegrateHistoInvMass( fHistoMappingBackNormInvMassLeftPtBin[iPt], fMesonCurLeftIntRange);
		fBckYieldsLeft[iPt] = fYields;
		fBckYieldsLeftError[iPt] = fYieldsError;
		IntegrateHistoInvMass( fHistoMappingBackNormInvMassLeftPtBin[iPt], fMesonCurLeftIntRangeWide);
		fBckYieldsLeftWide[iPt] = fYields;
		fBckYieldsLeftWideError[iPt] = fYieldsError;
		IntegrateHistoInvMass( fHistoMappingBackNormInvMassLeftPtBin[iPt], fMesonCurLeftIntRangeNarrow);
		fBckYieldsLeftNarrow[iPt] = fYields;
		fBckYieldsLeftNarrowError[iPt] = fYieldsError;

		// Integrate the signal histo
		fFileDataLog<< endl <<"Signal histo normal range, left norm " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
		IntegrateHistoInvMassStream( fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonCurLeftIntRange);
		fMesonYieldsLeft[iPt] = fYields;
		fMesonYieldsLeftError[iPt] = fYieldsError;
		fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
		fFileDataLog<< endl <<"Signal histo wide range, left norm " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
		IntegrateHistoInvMassStream( fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonCurLeftIntRangeWide);
		fMesonYieldsLeftWide[iPt] = fYields;
		fMesonYieldsLeftWideError[iPt] = fYieldsError;
		fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
		fFileDataLog<< endl <<"Signal histo narrow range, left norm " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
		IntegrateHistoInvMassStream( fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonCurLeftIntRangeNarrow);
		fMesonYieldsLeftNarrow[iPt] = fYields;
		fMesonYieldsLeftNarrowError[iPt] = fYieldsError;
		fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;

		fFileDataLog<< "Residual Background leftover norm integration/left norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFuncLeft[iPt] <<"\t +- \t" << fMesonYieldsResidualBckFuncLeftError[iPt]<<endl<< endl;
		fMesonYieldsCorResidualBckFuncLeft[iPt] = fMesonYieldsLeft[iPt]- fMesonYieldsResidualBckFuncLeft[iPt];
		fMesonYieldsCorResidualBckFuncLeftError[iPt] =
			pow((fMesonYieldsLeftError[iPt]*fMesonYieldsLeftError[iPt]+
				fMesonYieldsResidualBckFuncLeftError[iPt]*fMesonYieldsResidualBckFuncLeftError[iPt]),0.5);
		fMesonYieldsLeftPerEvent[iPt]= fMesonYieldsCorResidualBckFuncLeft[iPt]/fNEvents;
		fMesonYieldsLeftPerEventError[iPt]= fMesonYieldsCorResidualBckFuncLeftError[iPt]/fNEvents;

		fTotalBckYieldsLeft[iPt] = fBckYieldsLeft[iPt] + fMesonYieldsResidualBckFuncLeft[iPt];
		fTotalBckYieldsLeftError[iPt] = pow(fBckYieldsLeftError[iPt]*fBckYieldsLeftError[iPt] + fMesonYieldsResidualBckFuncLeftError[iPt]*fMesonYieldsResidualBckFuncLeftError[iPt],0.5);

		//Integrate Fit Function
		IntegrateFitFunc( fFitSignalPeakPosInvMassLeftPtBin[iPt], fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonCurLeftIntRange);
		fMesonYieldsFuncLeft[iPt]=fYieldsFunc;

		//GetFWHM
		
		if( fFitBckInvMassLeftPtBin[iPt]->Integral(fMesonMassLeft[iPt]-fMesonFWHMLeft[iPt], fMesonMassLeft[iPt]+fFitInvMassLeftPtBin[iPt]->GetParameter(2))!=0){
			Double_t background = fFitBckInvMassLeftPtBin[iPt]->Integral(fMesonMassLeft[iPt]-fMesonFWHMLeft[iPt], fMesonMassLeft[iPt]+fFitInvMassLeftPtBin[iPt]->GetParameter(2));
			Double_t backgroundErr = fFitBckInvMassLeftPtBin[iPt]->IntegralError(fMesonMassLeft[iPt]-fMesonFWHMLeft[iPt], fMesonMassLeft[iPt]+fFitInvMassLeftPtBin[iPt]->GetParameter(2));
			Double_t signal = fFitInvMassLeftPtBin[iPt]->Integral(fMesonMassLeft[iPt]-fMesonFWHMLeft[iPt], fMesonMassLeft[iPt]+fFitInvMassLeftPtBin[iPt]->GetParameter(2)) - background;
			Double_t signalErr =pow( pow(fFitInvMassLeftPtBin[iPt]->IntegralError(fMesonMassLeft[iPt]-fMesonFWHMLeft[iPt], fMesonMassLeft[iPt]+fFitInvMassLeftPtBin[iPt]->GetParameter(2)),2 )+ pow(backgroundErr,2),0.5);
			fMesonSBLeft[iPt] = signal/ background;
			fMesonSBLeftError[iPt] = pow( pow(signalErr/background,2.)+pow(signal/(background *background )*backgroundErr ,2.) ,0.5); 
			fMesonSignLeft[iPt] = signal/ pow(background + signal,0.5);
			fMesonSignLeftError[iPt] = pow(pow( (pow(background + signal ,0.5) - 0.5*signal*pow(background+signal,-0.5))/(background+signal) * signalErr ,2) + pow( 0.5*pow(signal+background,-1.5),2) ,0.5); 
		}else{
			fMesonSBLeft[iPt] = 0.;
			fMesonSBLeftError[iPt] = 0.;
			fMesonSignLeft[iPt] = 0.;
			fMesonSignLeftError[iPt] = 0.;
		}

		// Wide integration mass window
		//       cout<< "iPt"<< iPt<< " "<< "wide range"<<endl;
		fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "wide range/left normalization" << endl;

		fFileErrLog << "Using exp fit"<<endl;
		FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntRangeWide,iPt,kTRUE);
		fMesonYieldsResidualBckFuncLeftWide[iPt] = fIntLinearBck;
		fMesonYieldsResidualBckFuncLeftWideError[iPt] = fIntLinearBckError;

		//FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt],fMesonIntRangeWide,iPt,kFALSE);
		fFileDataLog<< "Residual Background leftover wide integration/left norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFuncLeftWide[iPt]<<"\t +- \t" << fMesonYieldsResidualBckFuncLeftWideError[iPt] <<endl<< endl;
		fMesonYieldsCorResidualBckFuncLeftWide[iPt] = fMesonYieldsLeftWide[iPt]- fMesonYieldsResidualBckFuncLeftWide[iPt];
		fMesonYieldsCorResidualBckFuncLeftWideError[iPt] =
			pow((fMesonYieldsLeftWideError[iPt]*fMesonYieldsLeftWideError[iPt]+ fMesonYieldsResidualBckFuncLeftWideError[iPt]*fMesonYieldsResidualBckFuncLeftWideError[iPt]),0.5);
		fMesonYieldsLeftPerEventWide[iPt]= fMesonYieldsCorResidualBckFuncLeftWide[iPt]/fNEvents;
		fMesonYieldsLeftPerEventWideError[iPt]= fMesonYieldsCorResidualBckFuncLeftWideError[iPt]/fNEvents;

		fTotalBckYieldsLeftWide[iPt] = fBckYieldsLeftWide[iPt] + fMesonYieldsResidualBckFuncLeftWide[iPt];
		fTotalBckYieldsLeftWideError[iPt] = pow(fBckYieldsLeftWideError[iPt]*fBckYieldsLeftWideError[iPt] + fMesonYieldsResidualBckFuncLeftWideError[iPt]*fMesonYieldsResidualBckFuncLeftWideError[iPt],0.5);


		if( fTotalBckYieldsLeftWide[iPt]!=0){
			fMesonSBLeftWide[iPt] = fMesonYieldsCorResidualBckFuncLeftWide[iPt]/fTotalBckYieldsLeftWide[iPt];
			fMesonSBLeftWideError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncLeftWideError[iPt]/fTotalBckYieldsLeftWide[iPt],2.)+pow(fMesonYieldsCorResidualBckFuncLeftWide[iPt]/(fTotalBckYieldsLeftWide[iPt]*fTotalBckYieldsLeftWide[iPt])*fTotalBckYieldsLeftWideError[iPt],2.) ,0.5);
			fMesonSignLeftWide[iPt] = fMesonYieldsCorResidualBckFuncLeftWide[iPt]/pow(fTotalBckYieldsLeftWide[iPt],0.5);
			fMesonSignLeftWideError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncLeftWideError[iPt]/pow(fTotalBckYieldsLeftWide[iPt],0.5),2.)+pow(0.5*fMesonYieldsCorResidualBckFuncLeftWide[iPt]/pow(fTotalBckYieldsLeftWide[iPt],1.5)*fTotalBckYieldsLeftWideError[iPt],2.) ,0.5);
		}else{
			fMesonSBLeftWide[iPt] = 0.;
			fMesonSBLeftWideError[iPt] = 0.;
			fMesonSignLeftWide[iPt] = 0.;
			fMesonSignLeftWideError[iPt] = 0.;
		}


		// Narrow integration mass window
		//       cout<< "iPt"<< iPt<< " "<< "narrow range"<<endl;
		fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "narrow range/left normalization" << endl;

		fFileErrLog << "Using exp fit"<<endl;
		FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntRangeNarrow,iPt,kTRUE);
		fMesonYieldsResidualBckFuncLeftNarrow[iPt] = fIntLinearBck;
		fMesonYieldsResidualBckFuncLeftNarrowError[iPt] = fIntLinearBckError;

		//FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt],fMesonIntRangeNarrow,iPt,kFALSE);
		fFileDataLog<< "Residual Background leftover narrow integration/left norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFuncLeftNarrow[iPt]<<"\t +- \t" << fMesonYieldsResidualBckFuncLeftNarrowError[iPt] <<endl<< endl;
		fMesonYieldsCorResidualBckFuncLeftNarrow[iPt] = fMesonYieldsLeftNarrow[iPt]- fMesonYieldsResidualBckFuncLeftNarrow[iPt];
		fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt] =
			pow((fMesonYieldsLeftNarrowError[iPt]*fMesonYieldsLeftNarrowError[iPt]+
				fMesonYieldsResidualBckFuncLeftNarrowError[iPt]*fMesonYieldsResidualBckFuncLeftNarrowError[iPt]),0.5);
		fMesonYieldsLeftPerEventNarrow[iPt]= fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]/fNEvents;
		fMesonYieldsLeftPerEventNarrowError[iPt]= fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt]/fNEvents;

		fTotalBckYieldsLeftNarrow[iPt] = fBckYieldsLeftNarrow[iPt] + fMesonYieldsResidualBckFuncLeftNarrow[iPt];
		fTotalBckYieldsLeftNarrowError[iPt] = pow(fBckYieldsLeftNarrowError[iPt]*fBckYieldsLeftNarrowError[iPt] + fMesonYieldsResidualBckFuncLeftNarrowError[iPt]*fMesonYieldsResidualBckFuncLeftNarrowError[iPt],0.5);


		if( fTotalBckYieldsLeftNarrow[iPt]!=0){
			fMesonSBLeftNarrow[iPt] = fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]/fTotalBckYieldsLeftNarrow[iPt];
			fMesonSBLeftNarrowError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt]/fTotalBckYieldsLeftNarrow[iPt],2.)+pow(fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]/(fTotalBckYieldsLeftNarrow[iPt]*fTotalBckYieldsLeftNarrow[iPt])*fTotalBckYieldsLeftNarrowError[iPt],2.) ,0.5);
			fMesonSignLeftNarrow[iPt] = fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]/pow(fTotalBckYieldsLeftNarrow[iPt],0.5);
			fMesonSignLeftNarrowError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt]/pow(fTotalBckYieldsLeftNarrow[iPt],0.5),2.)+pow(0.5*fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]/pow(fTotalBckYieldsLeftNarrow[iPt],1.5)*fTotalBckYieldsLeftNarrowError[iPt],2.) ,0.5);

		}else{
			fMesonSBLeftNarrow[iPt] = 0.;
			fMesonSBLeftNarrowError[iPt] = 0.;
			fMesonSignLeftNarrow[iPt] = 0.;
			fMesonSignLeftNarrowError[iPt] = 0.;
		}

	}



	//******************** Data OUTPUTFILE ***************************************************
	const char* fileNameSysErrDat = Form("%s/%s/%s_%s_SystematicErrorYieldExtraction_RAWDATA_%s.dat",fCutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPrefix2.Data(), cutSelection.Data());
	fstream fileSysErrDat;
	fileSysErrDat.open(fileNameSysErrDat, ios::out);
	fileSysErrDat << "Calculation of the systematic error due to the yield extraction RAWDATA " << endl;
	fileSysErrDat <<  endl;
	fileSysErrDat << "fPiPlPiMiGammaYields" << endl;
	fileSysErrDat << "Bin \t Right \t Right Wide \t Right Narr" << endl;
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << fPiPlPiMiGammaYields[iPt] << "+-" << fPiPlPiMiGammaYieldsError[iPt] << "\t" <<
			fPiPlPiMiGammaYieldsWide[iPt] << "+-" << fPiPlPiMiGammaYieldsWideError[iPt] << "\t" <<
			fPiPlPiMiGammaYieldsNarrow[iPt] << "+-" << fPiPlPiMiGammaYieldsNarrowError[iPt] << endl;

	}
	fileSysErrDat <<  endl;
	fileSysErrDat << "fTotalBckYields" << endl;
	fileSysErrDat << "Bin \t Right \t Right Wide \t Right Narr \t Left \t Left Wide \t Left Narr" << endl;
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" <<
			fTotalBckYields[iPt] << "+-" << fTotalBckYieldsError[iPt] << "\t" <<
			fTotalBckYieldsWide[iPt] << "+-" << fTotalBckYieldsWideError[iPt] << "\t" <<
			fTotalBckYieldsNarrow[iPt] << "+-" << fTotalBckYieldsNarrowError[iPt] << "\t" <<
			fTotalBckYieldsLeft[iPt] << "+-" << fTotalBckYieldsLeftError[iPt]<< "\t" <<
			fTotalBckYieldsLeftWide[iPt]<< "+-" << fTotalBckYieldsLeftWideError[iPt]<< "\t" <<
			fTotalBckYieldsLeftNarrow[iPt]<< "+-" << fTotalBckYieldsLeftNarrowError[iPt] << endl;
	}
	fileSysErrDat <<  endl;
	fileSysErrDat << "fMesonYields" << endl;
	fileSysErrDat << "Bin \t Right \t Right Wide \t Right Narr \t Left \t Left Wide \t Left Narr" << endl;
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" <<
			fMesonYieldsCorResidualBckFunc[iPt] << "+-" << fMesonYieldsCorResidualBckFuncError[iPt] << "\t" <<
			fMesonYieldsCorResidualBckFuncWide[iPt] << "+-" << fMesonYieldsCorResidualBckFuncWideError[iPt] << "\t" <<
			fMesonYieldsCorResidualBckFuncNarrow[iPt] << "+-" << fMesonYieldsCorResidualBckFuncNarrowError[iPt] << "\t" <<
			fMesonYieldsCorResidualBckFuncLeft[iPt]<< "+-" << fMesonYieldsCorResidualBckFuncLeftError[iPt]<< "\t" <<
			fMesonYieldsCorResidualBckFuncLeftWide[iPt]<< "+-" << fMesonYieldsCorResidualBckFuncLeftWideError[iPt]<< "\t" <<
			fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]<< "+-" << fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt] << endl;
	}
	if(fIsMC){
		fileSysErrDat <<  endl;
		fileSysErrDat << "TrueYields" << endl;
		fileSysErrDat << "Bin \t True \t True Wide \t True Narr " << endl;
		for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
			fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" <<
				fMesonTrueYields[iPt] << "\t" <<
				fMesonTrueYieldsWide[iPt] << "\t" <<
				fMesonTrueYieldsNarrow[iPt] << endl;
		}
	}
	fileSysErrDat.close();
	//******************************** OUTPUT END ******************************************************

	TString nameMeson = Form("%s/%s_%s_MesonWithBck%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
	TString nameCanvas = "MesonWithBckCanvas";
	TString namePad = "MesonWithBckPad";
	PlotInvMassInPtBins( fHistoMappingPiPlPiMiGammaInvMassPtBin, fHistoMappingBackNormInvMassPtBin, nameMeson, nameCanvas, namePad, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fTextCent);
		
	TString nameMesonSub= Form("%s/%s_%s_MesonSubtracted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(), fPeriodFlag.Data(), fCutSelection.Data(),Suffix.Data());
	TString nameCanvasSub= "MesonCanvasSubtracted";
	TString namePadSub= "MesonPadSubtracted";
	PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalInvMassPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitSignalInvMassPtBin, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fTextCent);

	nameMesonSub= Form("%s/%s_%s_MesonSubtractedBackFit%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(), fPeriodFlag.Data(), fCutSelection.Data(),Suffix.Data());
	nameCanvasSub= "MesonCanvasSubtractedBackFit";
	namePadSub= "MesonPadSubtractedBackFit";
	PlotWithFitSubtractedInvMassInPtBins( fHistoMappingGGInvMassBackFitPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitSignalInvMassBackFitPtBin, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fTextCent);

	nameMesonSub= Form("%s/%s_%s_MesonBackFit%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(), fPeriodFlag.Data(), fCutSelection.Data(),Suffix.Data());
	nameCanvasSub= "MesonCanvasBackFit";
	namePadSub= "MesonPadBackFit";
	PlotWithFitSubtractedInvMassInPtBins( fHistoMappingPiPlPiMiGammaInvMassPtBin, fHistoMappingTrueMesonInvMassPtBins, fBackgroundFitPol, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fTextCent);

	nameMeson= Form("%s/%s_%s_MesonWithBckLeft%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
	nameCanvas = "MesonWithBckCanvasLeft";
	namePad = "MesonWithBckPadLeft";
	PlotInvMassInPtBins( fHistoMappingPiPlPiMiGammaInvMassPtBin, fHistoMappingBackNormInvMassLeftPtBin, nameMeson, nameCanvas, namePad,  fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fTextCent);

	nameMesonSub= Form("%s/%s_%s_MesonWithBGAndFit%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
	nameCanvasSub= "MesonCanvasWithBGAndFit";
	namePadSub= "MesonPadWithBGAndFit";
	PlotWithFitSubtractedInvMassInPtBins( fHistoMappingPiPlPiMiGammaInvMassPtBin, fHistoMappingTrueMesonInvMassPtBins,fFitWithPol2ForBG, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC , fDecayChannel, fTextCent);

	nameMesonSub= Form("%s/%s_%s_MesonSubtractedLeft%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
	nameCanvasSub= "MesonCanvasSubtractedLeft";
	namePadSub= "MesonPadSubtractedLeft";
	PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalInvMassLeftPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitInvMassLeftPtBin, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fTextCent);

	PlotExampleInvMassBins(fHistoMappingPiPlPiMiGammaInvMassPtBin[fExampleBin], fHistoMappingSignalInvMassPtBin[fExampleBin], fHistoMappingBackNormInvMassPtBin[fExampleBin] , fFitSignalInvMassPtBin[fExampleBin], fExampleBin, outputDir.Data(),Suffix.Data(), fMesonMassRange, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2, fThesis, fCollisionSystem, fBinsPt, fDecayChannel);
		
	if(fIsMC){
		TString nameMesonTrue= Form("%s/%s_%s_TrueMesonFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
		TString nameCanvasTrue= "TrueMesonCanvasFitted";
		TString namePadTrue= "TrueMesonPadFitted";
		PlotWithFitSubtractedInvMassInPtBins(fHistoMappingTrueMesonInvMassPtBins, fHistoMappingTrueMesonInvMassPtBins, fFitTrueSignalInvMassPtBin, nameMesonTrue, nameCanvasTrue, namePadTrue, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fTextCent);
			
		nameMesonTrue= Form("%s/%s_%s_TrueMesonContamination%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
		nameCanvasTrue= "TrueMesonCanvasSec";
		namePadTrue= "TrueMesonPadSec";
		PlotInvMassSecondaryInPtBins( fHistoMappingTrueMesonInvMassPtBins, fHistoMappingTrueGGMesonInvMassPtBins, fHistoMappingTrueDalitzInvMassPtBins, nameMesonTrue, nameCanvasTrue, namePadTrue, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel , fTextCent);
	}

	CreatePtHistos();
	FillPtHistos();

	if(fIsMC){
	
		FillHistosArrayMC(fHistoMCMesonPtWithinAcceptance,fHistoMCMesonPt,fDeltaPt,Form("%s/%s_%s_MCYieldMesonFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));

		CalculateMesonAcceptance();
		//       cout << "Calculated MesonAcceptance" << endl;

		fNameHistoEffi="MesonEffiPt";
		cout << fNameHistoEffi.Data() << endl;
		CalculateMesonEfficiency(fHistoYieldMeson,fHistoYieldTrueGGMeson,fNameHistoEffi,Form("%s/%s_%s_RecYieldMesonFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));
		fHistoMonteMesonEffiPt = fHistoMCMesonEffiPt;
		fHistoMonteMesonEffiFitPt = fHistoMCMesonEffiFitPt;
			
		fNameHistoEffi="MesonEffiBackFitPt";
		CalculateMesonEfficiency(fHistoYieldMesonBackFit,fHistoYieldTrueGGMeson,fNameHistoEffi,Form("%s/%s_%s_RecYieldMesonBackFitFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));
		fHistoMonteMesonEffiBackFitPt = fHistoMCMesonEffiPt;

		fNameHistoEffi="MesonWideEffiPt";
		cout << fNameHistoEffi.Data() << endl;
		CalculateMesonEfficiency(fHistoYieldMesonWide,fHistoYieldTrueGGMesonWide,fNameHistoEffi,Form("%s/%s_%s_RecYieldMesonWideFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));
		fHistoMonteMesonWideEffiPt = fHistoMCMesonEffiPt;
		fHistoMonteMesonWideEffiFitPt = fHistoMCMesonEffiFitPt;

		fNameHistoEffi="MesonNarrowEffiPt";
		cout << fNameHistoEffi.Data() << endl;
		CalculateMesonEfficiency(fHistoYieldMesonNarrow,fHistoYieldTrueGGMesonNarrow,fNameHistoEffi,Form("%s/%s_%s_RecYieldMesonNarrowFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));
		fHistoMonteMesonNarrowEffiPt = fHistoMCMesonEffiPt;
		fHistoMonteMesonNarrowEffiFitPt = fHistoMCMesonEffiFitPt;
			
		fNameHistoEffi="MesonLeftEffiPt";
		cout << fNameHistoEffi.Data() << endl;
		CalculateMesonEfficiency(fHistoYieldMesonLeft,fHistoYieldTrueGGMeson,fNameHistoEffi,Form("%s/%s_%s_RecYieldMesonLefFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));
		fHistoMonteMesonLeftEffiPt = fHistoMCMesonEffiPt;
		fHistoMonteMesonLeftEffiFitPt = fHistoMCMesonEffiFitPt;
			
		fNameHistoEffi="MesonLeftNarrowEffiPt";
		cout << fNameHistoEffi.Data() << endl;
		CalculateMesonEfficiency(fHistoYieldMesonLeftNarrow,fHistoYieldTrueGGMesonNarrow,fNameHistoEffi,Form("%s/%s_%s_RecYieldMesonLefNarrowFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));
		fHistoMonteMesonLeftNarrowEffiPt = fHistoMCMesonEffiPt;
		fHistoMonteMesonLeftNarrowEffiFitPt = fHistoMCMesonEffiFitPt;
			
		fNameHistoEffi="MesonLeftWideEffiPt";
		cout << fNameHistoEffi.Data() << endl;
		CalculateMesonEfficiency(fHistoYieldMesonLeftWide,fHistoYieldTrueGGMesonWide,fNameHistoEffi,Form("%s/%s_%s_RecYieldMesonLeftWideFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));
		fHistoMonteMesonLeftWideEffiPt = fHistoMCMesonEffiPt;
		fHistoMonteMesonLeftWideEffiFitPt = fHistoMCMesonEffiFitPt;
			
		// True Meson (only once case, because no normalization
		fNameHistoEffi="TrueMesonEffiPt";
		cout << fNameHistoEffi.Data() << endl;
		CalculateMesonEfficiency(fHistoYieldTrueMeson,NULL,fNameHistoEffi,Form("%s/%s_%s_TrueYieldMesonFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));
		fHistoMCTrueMesonEffiPt = fHistoMCMesonEffiPt; 
		fHistoMCTrueMesonEffiFitPt = fHistoMCMesonEffiFitPt; 
			
		fNameHistoEffi="TrueMesonWideEffiPt";
		cout << fNameHistoEffi.Data() << endl;
		CalculateMesonEfficiency(fHistoYieldTrueMesonWide,NULL,fNameHistoEffi,Form("%s/%s_%s_TrueYieldMesonWideFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));
		fHistoMCTrueMesonWideEffiPt = fHistoMCMesonEffiPt;
		fHistoMCTrueMesonWideEffiFitPt = fHistoMCMesonEffiFitPt;
			
		fNameHistoEffi="TrueMesonNarrowEffiPt";
		cout << fNameHistoEffi.Data() << endl;
		CalculateMesonEfficiency(fHistoYieldTrueMesonNarrow,NULL,fNameHistoEffi,Form("%s/%s_%s_TrueYieldMesonNarrowFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));
		fHistoMCTrueMesonNarrowEffiPt = fHistoMCMesonEffiPt;
		fHistoMCTrueMesonNarrowEffiFitPt = fHistoMCMesonEffiFitPt;
		
		fNameHistoFrac="TrueGGContFrac";
		fHistoYieldTrueGGCont= CalulateContaminationFraction(fHistoYieldTrueMeson, fHistoYieldTrueGGMeson, fNameHistoFrac);
		fNameHistoFrac="TrueDalitzContFrac";
		fHistoYieldTrueDalitzCont= CalulateContaminationFraction(fHistoYieldTrueMeson, fHistoYieldTrueDalitzMeson, fNameHistoFrac);
			
		fNameHistoFrac="TrueGGContFracNarrow";
		fHistoYieldTrueGGContNarrow= CalulateContaminationFraction(fHistoYieldTrueMesonNarrow, fHistoYieldTrueGGMesonNarrow, fNameHistoFrac);
		fNameHistoFrac="TrueDalitzContFracNarrow";
		fHistoYieldTrueDalitzContNarrow= CalulateContaminationFraction(fHistoYieldTrueMesonNarrow, fHistoYieldTrueDalitzMesonNarrow, fNameHistoFrac);
			
		fNameHistoFrac="TrueGGContFracWide";
		fHistoYieldTrueGGContWide= CalulateContaminationFraction(fHistoYieldTrueMesonWide, fHistoYieldTrueGGMesonWide, fNameHistoFrac);
		fNameHistoFrac="TrueDalitzContFracWide";
		fHistoYieldTrueDalitzContWide= CalulateContaminationFraction(fHistoYieldTrueMesonWide, fHistoYieldTrueDalitzMesonWide, fNameHistoFrac);

		SaveCorrectionHistos(fCutSelection, fPrefix2);
	}
	SaveHistos(fIsMC, fCutSelection, fPrefix2);
	cout << "here" << endl;
	fFileErrLog.close();
	cout << "here" << endl;
	fFileDataLog.close();
	cout << "here" << endl;
	Delete();
	cout << "here" << endl;
}

void ProcessBckFitSubtraction(TH1D *fGammaGamma, Int_t i, Double_t * fPeakRangeDummy, Double_t *fFitRangeDummy)
{
	
	fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i] = (TH1D*)fGammaGamma->Clone(Form("GG_WithoutSigal_%i",i));
	fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->Sumw2();

	for (Int_t binx= 0; binx < fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetNbinsX()+1; binx++){
		if(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinContent(binx) == 0){
			fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinError(binx,1.);
			fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinContent(binx,0.);
		}
		if(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinCenter(binx) > fPeakRangeDummy[0] && fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinCenter(binx) < fPeakRangeDummy[1]){
			fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinContent(binx,0.0);
			fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinError(binx,0.0);
		}
	}

	fBackgroundFitPol[i] = NULL;
	fBackgroundFitPol[i] = new TF1("bla","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x",fFitRangeDummy[0],fFitRangeDummy[1]);
	fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->Fit(fBackgroundFitPol[i],"QRME0");

	for (Int_t binx= 0; binx < fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetNbinsX()+1; binx++){
		if(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinCenter(binx) > fFitRangeDummy[0] && fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinCenter(binx) < fFitRangeDummy[1]){
			fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinContent(binx,fBackgroundFitPol[i]->Eval(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinCenter(binx)));
			fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinError(binx,fGammaGamma->GetBinError(binx));
		}
	}
	
	fHistoMappingGGInvMassBackFitPtBin[i] = (TH1D*) fGammaGamma->Clone(Form("GG_SubtractedSignal_%i",i));
	fHistoMappingGGInvMassBackFitPtBin[i]->Sumw2();
	fHistoMappingGGInvMassBackFitPtBin[i]->Add(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i],-1);

}

void ProcessEM(TH1D* fGammaGamma, TH1D* fBck, Double_t * fBGFitRangeEM)
{
	for (Int_t binx= 0; binx < fGammaGamma->GetNbinsX()+1; binx++){
		if(fGammaGamma->GetBinContent(binx) == 0){
			fGammaGamma->SetBinError(binx,1.);
			fGammaGamma->SetBinContent(binx,0.);
		}
	}
	fBckNorm = (TH1D*)fBck->Clone("fBckNorm");
	fGammaGamma->Sumw2();
	fBck->Sumw2();
	fBckNorm->Sumw2();
	
	Double_t 	r= fGammaGamma->Integral(fGammaGamma->GetXaxis()->FindBin(fBGFitRangeEM[0]),fGammaGamma->GetXaxis()->FindBin(fBGFitRangeEM[1]));
	Double_t 	b= fBck->Integral(fBck->GetXaxis()->FindBin(fBGFitRangeEM[0]),fBck->GetXaxis()->FindBin(fBGFitRangeEM[1]));
	Double_t 	norm = 1;
	
	if(b != 0) norm = r/b;
	fBckNorm->Scale(norm);
	//    cout<<"r="<<r<<" b="<<b<<" r/b="<<r/b<< " " << endl;
	
	Int_t numberOfZeros = 0;
	for (Int_t i = 1; i < fBckNorm->GetNbinsX()+1; i++){
		if (fBckNorm->GetBinContent(i) == 0){
			numberOfZeros++;
			if (norm > 1.){
				fBckNorm->SetBinError(i,1.);
				fBckNorm->SetBinContent(i,0.);
			}
		}
	}
	//    cout << "counted " << numberOfZeros << " in the normalized BG : "<< (Double_t)numberOfZeros/fBck->GetNbinsX() << endl;
	
	fSignal = (TH1D*)fGammaGamma->Clone("fSignal");
	fSignal->Sumw2();
	if ((Double_t)numberOfZeros/fBck->GetNbinsX()< 0.25) fSignal->Add(fBckNorm,-1.);
}


void ProcessRatioSignalBackground(TH1D* fGammaGamma, TH1D* fBck)
{
	for (Int_t binx= 0; binx < fGammaGamma->GetNbinsX()+1; binx++){
		if(fGammaGamma->GetBinContent(binx) == 0){
			fGammaGamma->SetBinError(binx,1.);
			fGammaGamma->SetBinContent(binx,0.);
		}
	}

	fGammaGamma->Sumw2();
	fBck->Sumw2();
	fRatioSB = (TH1D*)fGammaGamma->Clone("RatioSB");
	fRatioSB->Divide(fGammaGamma, fBck, 1.,1.,"B");
}

void ProduceBckProperWeighting(TList* fBackgroundContainer,TList* fMotherContainer){
	
	THnSparseF* fSparseMotherZM;
	THnSparseF* fSparseBckZM;
	fSparseMotherZM = (THnSparseF*)fMotherContainer->FindObject("Back_Mother_InvMass_Pt_z_m");
	fSparseBckZM = (THnSparseF*)fBackgroundContainer->FindObject("Back_Back_InvMass_Pt_z_m");
	
	
	
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fHistoWeightsBGZbinVsMbin[iPt] = new  TH2F("BGWeights", "", fSparseMotherZM->GetAxis(2)->GetNbins(),  0, fSparseMotherZM->GetAxis(2)->GetNbins(),fSparseMotherZM->GetAxis(3)->GetNbins(),  0, fSparseMotherZM->GetAxis(3)->GetNbins());
		fHistoWeightsBGZbinVsMbin[iPt]->GetYaxis()->SetTitle("M-bins");
		fHistoWeightsBGZbinVsMbin[iPt]->GetXaxis()->SetTitle("Z-bins");
		fHistoWeightsBGZbinVsMbin[iPt]->Sumw2();
		fHistoFillPerEventBGZbinVsMbin[iPt] = new  TH2F("BGPoolsFillstatus", "", fSparseMotherZM->GetAxis(2)->GetNbins(),  0, fSparseMotherZM->GetAxis(2)->GetNbins(),fSparseMotherZM->GetAxis(3)->GetNbins(),  0, fSparseMotherZM->GetAxis(3)->GetNbins());
		fHistoFillPerEventBGZbinVsMbin[iPt]->GetYaxis()->SetTitle("M-bins");
		fHistoFillPerEventBGZbinVsMbin[iPt]->GetXaxis()->SetTitle("Z-bins");
		fHistoFillPerEventBGZbinVsMbin[iPt]->Sumw2();
	
		for (Int_t z=0;z < fSparseMotherZM->GetAxis(2)->GetNbins()-1;z++){
			for (Int_t m = 0; m < fSparseMotherZM->GetAxis(3)->GetNbins(); m++){ 
				// pt
				fSparseMotherZM->GetAxis(1)->SetRange((fSparseMotherZM->GetAxis(1))->FindBin(fBinsPt[iPt]+0.001),(fSparseMotherZM->GetAxis(1))->FindBin(fBinsPt[iPt+1]-0.001));
				fSparseBckZM->GetAxis(1)->SetRange((fSparseBckZM->GetAxis(1))->FindBin(fBinsPt[iPt]+0.001),(fSparseBckZM->GetAxis(1))->FindBin(fBinsPt[iPt+1]-0.001));
				// z
				fSparseMotherZM->GetAxis(2)->SetRange(z, z);
				fSparseBckZM->GetAxis(2)->SetRange(z, z);
				// m
				fSparseMotherZM->GetAxis(3)->SetRange(m,m);
				fSparseBckZM->GetAxis(3)->SetRange(m,m);
				
				fHistoMotherZMProj = (TH1D*)fSparseMotherZM->Projection(0);
				fHistoMotherZMProj->Sumw2();
				fHistoBckZMProj = (TH1D*)fSparseBckZM->Projection(0);
				fHistoBckZMProj->Sumw2();
			
				fScalingFactorBck[z][m]= 1./fBackgroundMultNumber;
				if (m==0 && z ==0){
				if(fHistoMappingBackInvMassPtBin[iPt]!= NULL){
					delete fHistoMappingBackInvMassPtBin[iPt];
					fHistoMappingBackInvMassPtBin[iPt]=NULL;
				}
				fNameHistoBack = Form("Mapping_Back_InvMass_in_Pt_Bin%02d", iPt);
				fHistoMappingBackInvMassPtBin[iPt]= (TH1D*)fHistoBckZMProj->Clone(fNameHistoBack);
				fHistoMappingBackInvMassPtBin[iPt]->Sumw2();
				for (Int_t ii = 0; ii < fHistoBckZMProj->GetNbinsX()+1; ii++){
					fHistoMappingBackInvMassPtBin[iPt]->SetBinContent(ii,0.);
					fHistoMappingBackInvMassPtBin[iPt]->SetBinError(ii,0.);
				}
				}
				Int_t startBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[0]);
				Int_t endBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[1]);
				if (fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral) != 0) {
				fScalingFactorBck[z][m] = fHistoMotherZMProj->Integral(startBinIntegral,endBinIntegral)/fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral);
	//                cout << z << "\t" << m << "\t"  << fScalingFactorBck[z][m] << "\t" << 20./fBackgroundMultNumber << endl;
				if ( fScalingFactorBck[z][m]> (20./fBackgroundMultNumber) ){
	//                   cout << "fail safe entered" << endl;
					fScalingFactorBck[z][m]=1./fBackgroundMultNumber;
	//                   cout << "\t" <<  z << "\t" << m << "\t"  << fScalingFactorBck[z][m] << "\t" << 20./fBackgroundMultNumber << endl;
				}
				}
				fHistoMappingBackInvMassPtBin[iPt]->Add(fHistoBckZMProj,fScalingFactorBck[z][m]);
				fHistoWeightsBGZbinVsMbin[iPt]->Fill(z+0.5,m+0.5,fScalingFactorBck[z][m]);
				fHistoFillPerEventBGZbinVsMbin[iPt]->Fill(z+0.5,m+0.5,fHistoBckZMProj->GetEntries());
				fHistoMotherZMProj->Clear();
				fHistoBckZMProj->Clear();

			}
		}
		fHistoMappingBackInvMassPtBin[iPt]->Rebin(fNRebin[iPt]);
		//fHistoMappingBackInvMassPtBin[iPt]->Scale(1./fNRebin[iPt]);
		for (Int_t ii = 0; ii < fHistoMappingBackInvMassPtBin[iPt]->GetNbinsX()+1; ii++){
			if(fHistoMappingBackInvMassPtBin[iPt]->GetBinContent(ii) == 0){
				fHistoMappingBackInvMassPtBin[iPt]->SetBinContent(ii,0.);
				fHistoMappingBackInvMassPtBin[iPt]->SetBinError(ii,1.);
			}
		}
		fFileDataLog << "Scaling Background factors for Pt bin " << iPt << " z m " << endl;
		// 		cout << "Scaling Background factors for Pt bin " << iPt << " z m " << endl;
		for (Int_t z=0; z < fSparseMotherZM->GetAxis(2)->GetNbins()-1; z++){
			fFileDataLog << fScalingFactorBck[z][0] << "\t" << fScalingFactorBck[z][1] << "\t" << fScalingFactorBck[z][2] << "\t" << fScalingFactorBck[z][3] << endl;
			// 			cout << fScalingFactorBck[z][0] << "\t" << fScalingFactorBck[z][1] << "\t" << fScalingFactorBck[z][2] << "\t" << fScalingFactorBck[z][3] << endl;
		}
	}
	
	for (Int_t z=0;z < fSparseMotherZM->GetAxis(2)->GetNbins()-1;z++){
		for (Int_t m = 0; m < fSparseMotherZM->GetAxis(3)->GetNbins(); m++) {

			// pt
			fSparseMotherZM->GetAxis(1)->SetRange((fSparseMotherZM->GetAxis(1))->FindBin(fFullPt[0]+0.001),(fSparseMotherZM->GetAxis(1))->FindBin(fFullPt[1]-0.001));
			fSparseBckZM->GetAxis(1)->SetRange((fSparseBckZM->GetAxis(1))->FindBin(fFullPt[0]+0.001),(fSparseBckZM->GetAxis(1))->FindBin(fFullPt[1]-0.001));
			// z
			fSparseMotherZM->GetAxis(2)->SetRange(z+1, z+1);
			fSparseBckZM->GetAxis(2)->SetRange(z+1, z+1);
			// m
			fSparseMotherZM->GetAxis(3)->SetRange(m+1,m+1);
			fSparseBckZM->GetAxis(3)->SetRange(m+1,m+1);
			
			fHistoMotherZMProj = (TH1D*)fSparseMotherZM->Projection(0);
			fHistoMotherZMProj->Sumw2();
			fHistoBckZMProj = (TH1D*)fSparseBckZM->Projection(0);
			fHistoBckZMProj->Sumw2();
			
			fScalingFactorBck[z][m]= 1./fBackgroundMultNumber;
			if (m==0 && z ==0){
				fNameHistoBack = "Mapping_Back_InvMass_FullPt";
				fMesonFullPtBackground = (TH1D*)fHistoBckZMProj->Clone(fNameHistoBack);
				fMesonFullPtBackground->Sumw2();
				for (Int_t ii = 0; ii < fHistoBckZMProj->GetNbinsX()+1; ii++){
				fMesonFullPtBackground->SetBinContent(ii,0.);
				fMesonFullPtBackground->SetBinError(ii,0.);
				}
			}
			Int_t startBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[0]);
			Int_t endBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[1]);
			if (fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral) != 0) {
				fScalingFactorBck[z][m] = fHistoMotherZMProj->Integral(startBinIntegral,endBinIntegral)/fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral);
				if ( fScalingFactorBck[z][m]>20./fBackgroundMultNumber ){
				fScalingFactorBck[z][m]=1./fBackgroundMultNumber;
				}
			}
			fMesonFullPtBackground->Add(fHistoBckZMProj,fScalingFactorBck[z][m]);

			fHistoMotherZMProj->Clear();
			fHistoBckZMProj->Clear();
		}
	}
	fMesonFullPtBackground->Rebin(fNRebin[4]);
	//fMesonFullPtBackground->Scale(1./fNRebin[4]);
	

	for (Int_t z=0;z < fSparseMotherZM->GetAxis(2)->GetNbins()-1;z++){
		for (Int_t m = 0; m < fSparseMotherZM->GetAxis(3)->GetNbins(); m++) {
			
			// pt
			fSparseMotherZM->GetAxis(1)->SetRange((fSparseMotherZM->GetAxis(1))->FindBin(fMidPt[0]+0.001),(fSparseMotherZM->GetAxis(1))->FindBin(fMidPt[1]-0.001));
			fSparseBckZM->GetAxis(1)->SetRange((fSparseBckZM->GetAxis(1))->FindBin(fMidPt[0]+0.001),(fSparseBckZM->GetAxis(1))->FindBin(fMidPt[1]-0.001));
			// z
			fSparseMotherZM->GetAxis(2)->SetRange(z+1, z+1);
			fSparseBckZM->GetAxis(2)->SetRange(z+1, z+1);
			// m
			fSparseMotherZM->GetAxis(3)->SetRange(m+1,m+1);
			fSparseBckZM->GetAxis(3)->SetRange(m+1,m+1);
			
			fHistoMotherZMProj = (TH1D*)fSparseMotherZM->Projection(0);
			fHistoMotherZMProj->Sumw2();
			fHistoBckZMProj = (TH1D*)fSparseBckZM->Projection(0);
			fHistoBckZMProj->Sumw2();
			
			fScalingFactorBck[z][m]= 1./fBackgroundMultNumber;
			if (m==0 && z ==0){
				fNameHistoBack = "Mapping_Back_InvMass_MidPt";
				fFittingHistMidPtBackground = (TH1D*)fHistoBckZMProj->Clone(fNameHistoBack);
				fFittingHistMidPtBackground->Sumw2();
				for (Int_t ii = 0; ii < fHistoBckZMProj->GetNbinsX()+1; ii++){
				fFittingHistMidPtBackground->SetBinContent(ii,0.);
				fFittingHistMidPtBackground->SetBinError(ii,0.);
				}
			}
			Int_t startBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[0]);
			Int_t endBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[1]);
			if (fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral) != 0) {
				fScalingFactorBck[z][m] = fHistoMotherZMProj->Integral(startBinIntegral,endBinIntegral)/fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral);
				if ( fScalingFactorBck[z][m]>20./fBackgroundMultNumber ){
				fScalingFactorBck[z][m]=1./fBackgroundMultNumber;
				}
			}
			fFittingHistMidPtBackground->Add(fHistoBckZMProj,fScalingFactorBck[z][m]);
			
			fHistoMotherZMProj->Clear();
			fHistoBckZMProj->Clear();
		}
	}
	fFittingHistMidPtBackground->Rebin(fNRebin[4]);
	//fFittingHistMidPtBackground->Scale(1./fNRebin[4]);
}

void ProduceBckWithoutWeighting(TH2D *fBckInvMassVSPtDummy){
	//calculation background for midPt without weighting
	Int_t startBinMidPt = fBckInvMassVSPtDummy->GetYaxis()->FindBin(fMidPt[0]+0.001);
	Int_t endBinMidPt = fBckInvMassVSPtDummy->GetYaxis()->FindBin(fMidPt[1]-0.001);
	fFittingHistMidPtBackground = new TH1D("Mapping_Back_InvMass_MidPt","Mapping_Back_InvMass_MidPt",fBckInvMassVSPtDummy->GetNbinsX(),0.,1.);
		fFittingHistMidPtBackground->Sumw2();
	fBckInvMassVSPtDummy->ProjectionX("Mapping_Back_InvMass_MidPt",startBinMidPt,endBinMidPt);
	fFittingHistMidPtBackground=(TH1D*)gDirectory->Get("Mapping_Back_InvMass_MidPt");
	fFittingHistMidPtBackground->Rebin(fNRebin[4]);

	//calulation background for fullPt without weighting
	fMesonFullPtBackground = new TH1D("Mapping_Back_InvMass_FullPt","Mapping_Back_InvMass_FullPt",fBckInvMassVSPtDummy->GetNbinsX(),0.,1.);
	fMesonFullPtBackground->Sumw2();
	Int_t startBinFullPt = fBckInvMassVSPtDummy->GetYaxis()->FindBin(fFullPt[0]+0.001);
	Int_t endBinFullPt = fBckInvMassVSPtDummy->GetYaxis()->FindBin(fFullPt[1]-0.001);
	fBckInvMassVSPtDummy->ProjectionX("Mapping_Back_InvMass_FullPt",startBinFullPt,endBinFullPt);
	fMesonFullPtBackground=(TH1D*)gDirectory->Get("Mapping_Back_InvMass_FullPt");
	fMesonFullPtBackground->Rebin(fNRebin[4]);

	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fNameHistoBack = Form("Mapping_Back_InvMass_in_Pt_Bin%02d", iPt);
		if(fHistoMappingBackInvMassPtBin[iPt]!= NULL){
			delete fHistoMappingBackInvMassPtBin[iPt];
			fHistoMappingBackInvMassPtBin[iPt]=NULL;
		}
		fHistoMappingBackInvMassPtBin[iPt]=new TH1D(fNameHistoBack.Data(),fNameHistoBack.Data(),fBckInvMassVSPtDummy->GetNbinsX(),0.,1.);
			fHistoMappingBackInvMassPtBin[iPt]->Sumw2();
		Int_t startBin = fBckInvMassVSPtDummy->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
		Int_t endBin = fBckInvMassVSPtDummy->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

		fBckInvMassVSPtDummy->ProjectionX(fNameHistoBack.Data(),startBin,endBin);
		fHistoMappingBackInvMassPtBin[iPt]=(TH1D*)gDirectory->Get(fNameHistoBack.Data());
		if(fNRebin[iPt]>1){
			fHistoMappingBackInvMassPtBin[iPt]->Rebin(fNRebin[iPt]);
		}

	}
}


void FillMassHistosArray(TH2D* fGammaGammaInvMassVSPtDummy)
{
	fFittingHistMidPtSignal = new TH1D("Mapping_GG_InvMass_MidPt","Mapping_GG_InvMass_MidPt",fGammaGammaInvMassVSPtDummy->GetNbinsX(),0.,(double)fGammaGammaInvMassVSPtDummy->GetXaxis()->GetBinUpEdge(fGammaGammaInvMassVSPtDummy->GetNbinsX()));
	cout << fGammaGammaInvMassVSPtDummy->GetXaxis()->GetBinUpEdge(fGammaGammaInvMassVSPtDummy->GetNbinsX()) << endl;
	Int_t startBinMidPt = fGammaGammaInvMassVSPtDummy->GetYaxis()->FindBin(fMidPt[0]+0.001);
	Int_t endBinMidPt = fGammaGammaInvMassVSPtDummy->GetYaxis()->FindBin(fMidPt[1]-0.001);
		fFittingHistMidPtSignal->Sumw2();
	fGammaGammaInvMassVSPtDummy->ProjectionX("Mapping_GG_InvMass_MidPt",startBinMidPt,endBinMidPt);

	fFittingHistMidPtSignal=(TH1D*)gDirectory->Get("Mapping_GG_InvMass_MidPt");
	fFittingHistMidPtSignal->Rebin(fNRebin[4]);
	//fFittingHistMidPtSignal->Scale(1./fNRebin[4]);
	//    cout << "Mid pt geschrieben" << endl;

	fMesonFullPtSignal = new TH1D("Mapping_GG_InvMass_FullPt","Mapping_GG_InvMass_FullPt",fGammaGammaInvMassVSPtDummy->GetNbinsX(),0.,fGammaGammaInvMassVSPtDummy->GetXaxis()->GetBinUpEdge(fGammaGammaInvMassVSPtDummy->GetNbinsX()));
		fMesonFullPtSignal->Sumw2();
	Int_t startBinFullPt = fGammaGammaInvMassVSPtDummy->GetYaxis()->FindBin(fFullPt[0]+0.001);
	Int_t endBinFullPt = fGammaGammaInvMassVSPtDummy->GetYaxis()->FindBin(fFullPt[1]-0.001);

	fGammaGammaInvMassVSPtDummy->ProjectionX("Mapping_GG_InvMass_FullPt",startBinFullPt,endBinFullPt);

	fMesonFullPtSignal=(TH1D*)gDirectory->Get("Mapping_GG_InvMass_FullPt");
	fMesonFullPtSignal->Rebin(fNRebin[4]);

	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fNameHistoGG = Form("Mapping_GG_InvMass_in_Pt_Bin%02d", iPt);

		if(fHistoMappingPiPlPiMiGammaInvMassPtBin[iPt]!= NULL){
			delete fHistoMappingPiPlPiMiGammaInvMassPtBin[iPt];
			fHistoMappingPiPlPiMiGammaInvMassPtBin[iPt]=NULL;
		}
		fHistoMappingPiPlPiMiGammaInvMassPtBin[iPt]=new TH1D(fNameHistoGG.Data(),fNameHistoGG.Data(),fGammaGammaInvMassVSPtDummy->GetNbinsX(),0.,fGammaGammaInvMassVSPtDummy->GetXaxis()->GetBinUpEdge(fGammaGammaInvMassVSPtDummy->GetNbinsX()));
			fHistoMappingPiPlPiMiGammaInvMassPtBin[iPt]->Sumw2();

		Int_t startBin = fGammaGammaInvMassVSPtDummy->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
		Int_t endBin = fGammaGammaInvMassVSPtDummy->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

		//       cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
		//       cout<< "bin values::"<< fGammaGammaInvMassVSPtDummy->GetYaxis()->GetBinCenter(startBin)<< " "
		//           << fGammaGammaInvMassVSPtDummy->GetYaxis()->GetBinCenter(endBin)<< endl;

		fGammaGammaInvMassVSPtDummy->ProjectionX(fNameHistoGG.Data(),startBin,endBin);


		fHistoMappingPiPlPiMiGammaInvMassPtBin[iPt]=(TH1D*)gDirectory->Get(fNameHistoGG.Data());

		if(fNRebin[iPt]>1){
			fHistoMappingPiPlPiMiGammaInvMassPtBin[iPt]->Rebin(fNRebin[iPt]);
			//fHistoMappingPiPlPiMiGammaInvMassPtBin[iPt]->Scale(1./fNRebin[iPt]);
		}
	}
	//    cout << "each pt written" << endl;

}
void FillMassMCTrueMesonHistosArray(TH2D* fHistoTrueMesonInvMassVSPtFill)
{
	fHistoTrueMesonInvMassVSPtFill->Sumw2();
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fNameHistoTrue = Form("Mapping_TrueMeson_InvMass_in_Pt_Bin%02d", iPt);
		if(fHistoMappingTrueMesonInvMassPtBins[iPt]!= NULL){
			delete fHistoMappingTrueMesonInvMassPtBins[iPt];
			fHistoMappingTrueMesonInvMassPtBins[iPt]=NULL;
		}
		fHistoMappingTrueMesonInvMassPtBins[iPt] = new TH1D(fNameHistoTrue.Data(),fNameHistoTrue.Data(),fHistoTrueMesonInvMassVSPtFill->GetNbinsX(),0.,fHistoTrueMesonInvMassVSPtFill->GetXaxis()->GetBinUpEdge(fHistoTrueMesonInvMassVSPtFill->GetNbinsX()));
			fHistoMappingTrueMesonInvMassPtBins[iPt]->Sumw2();
		Int_t startBin = fHistoTrueMesonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
		Int_t endBin = fHistoTrueMesonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

		//       cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
		//       cout<< "bin values::"<< fHistoTrueMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(startBin)<< " "
		//           << fHistoTrueMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(endBin)<< endl;

		fHistoTrueMesonInvMassVSPtFill->ProjectionX(fNameHistoTrue.Data(),startBin,endBin,"e");
			
		fHistoMappingTrueMesonInvMassPtBins[iPt]=(TH1D*)gDirectory->Get(fNameHistoTrue.Data());
			cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonInvMassPtBins[iPt]->GetEntries() << endl;
		if(fNRebin[iPt]>1){
			fHistoMappingTrueMesonInvMassPtBins[iPt]->Rebin(fNRebin[iPt]);
		}
		fHistoMappingTrueMesonInvMassPtBins[iPt]->SetLineWidth(1);
		fHistoMappingTrueMesonInvMassPtBins[iPt]->SetLineColor(2);
	}
}


void FillMassMCTrueGGMesonHistosArray(TH2D* fHistoTrueGGMesonInvMassVSPtFill)
{
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fNameHistoTrueGG = Form("Mapping_TrueGGMeson_InvMass_in_Pt_Bin%02d", iPt);
		if(fHistoMappingTrueGGMesonInvMassPtBins[iPt]!= NULL){
			delete fHistoMappingTrueGGMesonInvMassPtBins[iPt];
			fHistoMappingTrueGGMesonInvMassPtBins[iPt]=NULL;
		}
		fHistoMappingTrueGGMesonInvMassPtBins[iPt] = new TH1D(fNameHistoTrueGG.Data(),fNameHistoTrueGG.Data(),fHistoTrueGGMesonInvMassVSPtFill->GetNbinsX(),0.,fHistoTrueGGMesonInvMassVSPtFill->GetXaxis()->GetBinUpEdge(fHistoTrueGGMesonInvMassVSPtFill->GetNbinsX()));
			fHistoMappingTrueGGMesonInvMassPtBins[iPt]->Sumw2();
		Int_t startBin = fHistoTrueGGMesonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
		Int_t endBin = fHistoTrueGGMesonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

		//       cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
		//       cout<< "bin values::"<< fHistoTrueGGMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(startBin)<< " "
		//           << fHistoTrueGGMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(endBin)<< endl;

		fHistoTrueGGMesonInvMassVSPtFill->ProjectionX(fNameHistoTrueGG.Data(),startBin,endBin);
		fHistoMappingTrueGGMesonInvMassPtBins[iPt]=(TH1D*)gDirectory->Get(fNameHistoTrueGG.Data());
		if(fNRebin[iPt]>1){
			fHistoMappingTrueGGMesonInvMassPtBins[iPt]->Rebin(fNRebin[iPt]);
			//fHistoMappingTrueGGMesonInvMassPtBins[iPt]->Scale(1./fNRebin[iPt]);

		}
		fHistoMappingTrueGGMesonInvMassPtBins[iPt]->SetLineWidth(1);
		fHistoMappingTrueGGMesonInvMassPtBins[iPt]->SetLineColor(2);
	}
}

void FillMassMCTrueDalitzMesonHistosArray(TH2D* fHistoTrueDalitzMesonInvMassVSPtFill)
{
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fNameHistoTrueDalitz = Form("Mapping_TrueDalitzMeson_InvMass_in_Pt_Bin%02d", iPt);
		if(fHistoMappingTrueDalitzInvMassPtBins[iPt]!= NULL){
			delete fHistoMappingTrueDalitzInvMassPtBins[iPt];
			fHistoMappingTrueDalitzInvMassPtBins[iPt]=NULL;
		}
		fHistoMappingTrueDalitzInvMassPtBins[iPt] = new TH1D(fNameHistoTrueDalitz.Data(),fNameHistoTrueDalitz.Data(),fHistoTrueDalitzMesonInvMassVSPtFill->GetNbinsX(),0.,fHistoTrueDalitzMesonInvMassVSPtFill->GetXaxis()->GetBinUpEdge(fHistoTrueDalitzMesonInvMassVSPtFill->GetNbinsX()));
			fHistoMappingTrueDalitzInvMassPtBins[iPt]->Sumw2();
		Int_t startBin = fHistoTrueDalitzMesonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
		Int_t endBin = fHistoTrueDalitzMesonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

		//       cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
		//       cout<< "bin values::"<< fHistoTrueDalitzMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(startBin)<< " "
		//           << fHistoTrueDalitzMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(endBin)<< endl;

		fHistoTrueDalitzMesonInvMassVSPtFill->ProjectionX(fNameHistoTrueDalitz.Data(),startBin,endBin);
		fHistoMappingTrueDalitzInvMassPtBins[iPt]=(TH1D*)gDirectory->Get(fNameHistoTrueDalitz.Data());
		if(fNRebin[iPt]>1){
			fHistoMappingTrueDalitzInvMassPtBins[iPt]->Rebin(fNRebin[iPt]);
			//fHistoMappingTrueDalitzInvMassPtBins[iPt]->Scale(1./fNRebin[iPt]);

		}
		fHistoMappingTrueDalitzInvMassPtBins[iPt]->SetLineWidth(1);
		fHistoMappingTrueDalitzInvMassPtBins[iPt]->SetLineColor(2);
	}
}

TH1D* CalulateContaminationFraction(TH1D* histoRawYield, TH1D* histoRawYieldSec, TString nameHistoFrac){
	TH1D* histoFracSec = (TH1D*)histoRawYieldSec->Clone(nameHistoFrac.Data());
	histoFracSec->Divide(histoFracSec,histoRawYield,1.,1.,"B");
	return histoFracSec;
}

void CreatePtHistos(){

	fDeltaPt =			 	new TH1D("deltaPt","",fNBinsPt,fBinsPt);

	fHistoYieldMeson = 			new TH1D("histoYieldMeson","",fNBinsPt,fBinsPt);
	fHistoYieldMesonBackFit = 			new TH1D("histoYieldMesonBackFit","",fNBinsPt,fBinsPt);
	fHistoYieldTrueMeson = 		new TH1D("histoYieldTrueMeson","",fNBinsPt,fBinsPt);
	fHistoYieldTrueGGMeson = 		new TH1D("histoYieldTrueSecMeson","",fNBinsPt,fBinsPt);
	fHistoYieldTrueDalitzMeson = 		new TH1D("histoYieldTrueSecFromK0SMeson","",fNBinsPt,fBinsPt);

	fHistoYieldMesonPerEvent = 	new TH1D("histoYieldMesonPerEvent","",fNBinsPt,fBinsPt);
	fHistoYieldMesonPerEventBackFit = 	new TH1D("histoYieldMesonPerEventBackFit","",fNBinsPt,fBinsPt);
	fHistoSignMeson = 			new TH1D("histoSignMeson","",fNBinsPt,fBinsPt);
	fHistoSBMeson = 			new TH1D("histoSBMeson","",fNBinsPt,fBinsPt);
	fHistoMassMeson = 			new TH1D("histoMassMeson","",fNBinsPt,fBinsPt);
	fHistoWidthMeson = 			new TH1D("histoWidthMeson","",fNBinsPt,fBinsPt);
	fHistoFWHMMeson = 			new TH1D("histoFWHMMeson","",fNBinsPt,fBinsPt);
	fHistoTrueMassMeson = 		new TH1D("histoTrueMassMeson","",fNBinsPt,fBinsPt);
	fHistoTrueFWHMMeson = 		new TH1D("histoTrueFWHMMeson","",fNBinsPt,fBinsPt);
	fHistoTrueSignMeson = 		new TH1D("histoTrueSignMeson","",fNBinsPt,fBinsPt);
	fHistoTrueSBMeson = 		new TH1D("histoTrueSBMeson","",fNBinsPt,fBinsPt);

	fHistoYieldMesonNarrow = 	new TH1D("histoYieldMesonNarrow","",fNBinsPt,fBinsPt);
	fHistoYieldTrueMesonNarrow = 	new TH1D("histoYieldTrueMesonNarrow","",fNBinsPt,fBinsPt);
	fHistoYieldTrueGGMesonNarrow = 		new TH1D("histoYieldTrueSecMesonNarrow","",fNBinsPt,fBinsPt);
	fHistoYieldTrueDalitzMesonNarrow = 		new TH1D("histoYieldTrueSecFromK0SMesonNarrow","",fNBinsPt,fBinsPt);
	fHistoYieldMesonPerEventNarrow = new TH1D("histoYieldMesonPerEventNarrow","",fNBinsPt,fBinsPt);
	fHistoSignMesonNarrow = 		new TH1D("histoSignMesonNarrow","",fNBinsPt,fBinsPt);
	fHistoSBMesonNarrow = 		new TH1D("histoSBMesonNarrow","",fNBinsPt,fBinsPt);

	fHistoYieldMesonWide = 		new TH1D("histoYieldMesonWide","",fNBinsPt,fBinsPt);
	fHistoYieldTrueMesonWide = 	new TH1D("histoYieldTrueMesonWide","",fNBinsPt,fBinsPt);
	fHistoYieldTrueGGMesonWide = 		new TH1D("histoYieldTrueSecMesonWide","",fNBinsPt,fBinsPt);
	fHistoYieldTrueDalitzMesonWide = 		new TH1D("histoYieldTrueSecFromK0SMesonWide","",fNBinsPt,fBinsPt);
	fHistoYieldMesonPerEventWide = new TH1D("histoYieldMesonPerEventWide","",fNBinsPt,fBinsPt);
	fHistoSignMesonWide = 		new TH1D("histoSignMesonWide","",fNBinsPt,fBinsPt);
	fHistoSBMesonWide = 		new TH1D("histoSBMesonWide","",fNBinsPt,fBinsPt);

	// Histos for normalization at the left of the peak

	fHistoYieldMesonLeft = 		new TH1D("histoYieldMesonLeft","",fNBinsPt,fBinsPt);
	fHistoYieldMesonLeftPerEvent = new TH1D("histoYieldMesonLeftPerEvent","",fNBinsPt,fBinsPt);
	fHistoSignMesonLeft = 		new TH1D("histoSignMesonLeft","",fNBinsPt,fBinsPt);
	fHistoSBMesonLeft = 		new TH1D("histoSBMesonLeft","",fNBinsPt,fBinsPt);
	fHistoMassMesonLeft = 		new TH1D("histoMassMesonLeft","",fNBinsPt,fBinsPt);
	fHistoWidthMesonLeft = 		new TH1D("histoWidthMesonLeft","",fNBinsPt,fBinsPt);
	fHistoFWHMMesonLeft = 		new TH1D("histoFWHMMesonLeft","",fNBinsPt,fBinsPt);


	fHistoYieldMesonLeftNarrow = 	new TH1D("histoYieldMesonLeftNarrow","",fNBinsPt,fBinsPt);
	fHistoYieldMesonLeftPerEventNarrow = new TH1D("histoYieldMesonLeftPerEventNarrow","",fNBinsPt,fBinsPt);
	fHistoSignMesonLeftNarrow = 	new TH1D("histoSignMesonLeftNarrow","",fNBinsPt,fBinsPt);
	fHistoSBMesonLeftNarrow = 	new TH1D("histoSBMesonLeftNarrow","",fNBinsPt,fBinsPt);


	fHistoYieldMesonLeftWide = 	new TH1D("histoYieldMesonLeftWide","",fNBinsPt,fBinsPt);
	fHistoYieldMesonLeftPerEventWide = new TH1D("histoYieldMesonLeftPerEventWide","",fNBinsPt,fBinsPt);
	fHistoSignMesonLeftWide = 	new TH1D("histoSignMesonLeftWide","",fNBinsPt,fBinsPt);
	fHistoSBMesonLeftWide = 		new TH1D("histoSBMesonLeftWide","",fNBinsPt,fBinsPt);
}

void FillPtHistos()
{
	for(Int_t iPt=fStartPtBin+1;iPt<fNBinsPt+1;iPt++){

		fDeltaPt->SetBinContent(iPt,fBinsPt[iPt]-fBinsPt[iPt-1]);
		fDeltaPt->SetBinError(iPt,0);


		fHistoMassMeson->SetBinContent(iPt,fMesonMass[iPt-1]);
		fHistoMassMeson->SetBinError(iPt,fMesonMassError[iPt-1]);
		// fHistoWidthMeson->SetBinContent(iPt);
		fHistoFWHMMeson->SetBinContent(iPt,fMesonFWHM[iPt-1]);
		fHistoFWHMMeson->SetBinError(iPt,fMesonFWHMError[iPt-1]);

		if (fIsMC) {
			fHistoTrueMassMeson->SetBinContent(iPt,fMesonTrueMass[iPt-1]);
			fHistoTrueMassMeson->SetBinError(iPt,fMesonTrueMassError[iPt-1]);
			
			fHistoTrueFWHMMeson->SetBinContent(iPt,fMesonTrueFWHM[iPt-1]);
			fHistoTrueFWHMMeson->SetBinError(iPt,fMesonTrueFWHMError[iPt-1]);
			
			fHistoTrueSignMeson->SetBinContent(iPt,fMesonTrueSign[iPt-1]);
			fHistoTrueSignMeson->SetBinError(iPt,fMesonTrueSignError[iPt-1]);
			fHistoTrueSBMeson->SetBinContent(iPt,fMesonTrueSB[iPt-1]);
			fHistoTrueSBMeson->SetBinError(iPt,fMesonTrueSBError[iPt-1]);
		}

		fHistoSignMeson->SetBinContent(iPt,fMesonSign[iPt-1]);
		fHistoSignMeson->SetBinError(iPt,fMesonSignError[iPt-1]);
		fHistoSBMeson->SetBinContent(iPt,fMesonSB[iPt-1]);
		fHistoSBMeson->SetBinError(iPt,fMesonSBError[iPt-1]);

		fHistoYieldMeson->SetBinContent(iPt,fMesonYieldsCorResidualBckFunc[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMeson->SetBinError(iPt,fMesonYieldsCorResidualBckFuncError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		fHistoYieldMesonBackFit->SetBinContent(iPt,fMesonYieldsCorResidualBckFuncBackFit[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonBackFit->SetBinError(iPt,fMesonYieldsCorResidualBckFuncBackFitError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));


		if (fIsMC) {
			fHistoYieldTrueMeson->SetBinContent(iPt,fMesonTrueYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueMeson->SetBinError(iPt,fMesonTrueYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueGGMeson->SetBinContent(iPt,fMesonTrueGGContYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueGGMeson->SetBinError(iPt,fMesonTrueGGContYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueDalitzMeson->SetBinContent(iPt,fMesonTrueDalitzContYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueDalitzMeson->SetBinError(iPt,fMesonTrueDalitzContYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		
		}
		fHistoYieldMesonPerEvent->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFunc[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonPerEvent->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		fHistoYieldMesonPerEventBackFit->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncBackFit[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonPerEventBackFit->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncBackFitError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		// Narrow integration window
		fHistoSignMesonNarrow->SetBinContent(iPt,fMesonSignNarrow[iPt-1]);
		fHistoSignMesonNarrow->SetBinError(iPt,fMesonSignNarrowError[iPt-1]);
		fHistoSBMesonNarrow->SetBinContent(iPt,fMesonSBNarrow[iPt-1]);
		fHistoSBMesonNarrow->SetBinError(iPt,fMesonSBNarrowError[iPt-1]);

		fHistoYieldMesonNarrow->SetBinContent(iPt,fMesonYieldsCorResidualBckFuncNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonNarrow->SetBinError(iPt,fMesonYieldsCorResidualBckFuncNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		if (fIsMC) {
			fHistoYieldTrueMesonNarrow->SetBinContent(iPt,fMesonTrueYieldsNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueMesonNarrow->SetBinError(iPt,fMesonTrueYieldsNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueGGMesonNarrow->SetBinContent(iPt,fMesonTrueGGContYieldsNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueGGMesonNarrow->SetBinError(iPt,fMesonTrueGGContYieldsNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueDalitzMesonNarrow->SetBinContent(iPt,fMesonTrueDalitzContYieldsNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueDalitzMesonNarrow->SetBinError(iPt,fMesonTrueDalitzContYieldsNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		
		}
		fHistoYieldMesonPerEventNarrow->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonPerEventNarrow->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		// Wide integration window
		fHistoSignMesonWide->SetBinContent(iPt,fMesonSignWide[iPt-1]);
		fHistoSignMesonWide->SetBinError(iPt,fMesonSignWideError[iPt-1]);
		fHistoSBMesonWide->SetBinContent(iPt,fMesonSBWide[iPt-1]);
		fHistoSBMesonWide->SetBinError(iPt,fMesonSBWideError[iPt-1]);

		fHistoYieldMesonWide->SetBinContent(iPt,fMesonYieldsCorResidualBckFuncWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonWide->SetBinError(iPt,fMesonYieldsCorResidualBckFuncWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		if (fIsMC) {
			fHistoYieldTrueMesonWide->SetBinContent(iPt,fMesonTrueYieldsWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueMesonWide->SetBinError(iPt,fMesonTrueYieldsWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueGGMesonWide->SetBinContent(iPt,fMesonTrueGGContYieldsWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueGGMesonWide->SetBinError(iPt,fMesonTrueGGContYieldsWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueDalitzMesonWide->SetBinContent(iPt,fMesonTrueDalitzContYieldsWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueDalitzMesonWide->SetBinError(iPt,fMesonTrueDalitzContYieldsWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		
		}
		fHistoYieldMesonPerEventWide->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonPerEventWide->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		// Histos for integration at the left of the peak
		fHistoMassMesonLeft->SetBinContent(iPt,fMesonMassLeft[iPt-1]);
		fHistoMassMesonLeft->SetBinError(iPt,fMesonMassLeftError[iPt-1]);
		fHistoFWHMMesonLeft->SetBinContent(iPt,fMesonFWHMLeft[iPt-1]);
		fHistoFWHMMesonLeft->SetBinError(iPt,fMesonFWHMLeftError[iPt-1]);

		fHistoSignMesonLeft->SetBinContent(iPt,fMesonSignLeft[iPt-1]);
		fHistoSignMesonLeft->SetBinError(iPt,fMesonSignLeftError[iPt-1]);
		fHistoSBMesonLeft->SetBinContent(iPt,fMesonSBLeft[iPt-1]);
		fHistoSBMesonLeft->SetBinError(iPt,fMesonSBLeftError[iPt-1]);

		fHistoYieldMesonLeft->SetBinContent(iPt,fMesonYieldsCorResidualBckFuncLeft[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeft->SetBinError(iPt,fMesonYieldsCorResidualBckFuncLeftError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftPerEvent->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeft[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftPerEvent->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeftError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		// Narrow integration window
		fHistoSignMesonLeftNarrow->SetBinContent(iPt,fMesonSignLeftNarrow[iPt-1]);
		fHistoSignMesonLeftNarrow->SetBinError(iPt,fMesonSignLeftNarrowError[iPt-1]);
		fHistoSBMesonLeftNarrow->SetBinContent(iPt,fMesonSBLeftNarrow[iPt-1]);
		fHistoSBMesonLeftNarrow->SetBinError(iPt,fMesonSBLeftNarrowError[iPt-1]);

		fHistoYieldMesonLeftNarrow->SetBinContent(iPt,fMesonYieldsCorResidualBckFuncLeftNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftNarrow->SetBinError(iPt,fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftPerEventNarrow->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeftNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftPerEventNarrow->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		// Wide integration window
		fHistoSignMesonLeftWide->SetBinContent(iPt,fMesonSignLeftWide[iPt-1]);
		fHistoSignMesonLeftWide->SetBinError(iPt,fMesonSignLeftWideError[iPt-1]);
		fHistoSBMesonLeftWide->SetBinContent(iPt,fMesonSBLeftWide[iPt-1]);
		fHistoSBMesonLeftWide->SetBinError(iPt,fMesonSBLeftWideError[iPt-1]);

		fHistoYieldMesonLeftWide->SetBinContent(iPt,fMesonYieldsCorResidualBckFuncLeftWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftWide->SetBinError(iPt,fMesonYieldsCorResidualBckFuncLeftWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftPerEventWide->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeftWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftPerEventWide->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeftWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
	}
}





void FitSubtractedInvMassInPtBins(TH1D* fHistoMappingSignalInvMassPtBinSingle, Double_t* fMesonIntRangeFit, Int_t ptBin, Bool_t vary)
{

	//    cout<<"Start Fitting spectra"<<endl;
	fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassRange[0],fMesonMassRange[1]);
	Double_t mesonAmplitude =fHistoMappingSignalInvMassPtBinSingle->GetMaximum();
	Double_t mesonAmplitudeMin = mesonAmplitude*50./100.;
	Double_t mesonAmplitudeMax = mesonAmplitude*115./100.;

	fFitReco= NULL;
	fFitReco = new TF1("GaussExpLinear","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)",fMesonFitRange[0],fMesonFitRange[1]);

	fFitGausExp =NULL;
	fFitGausExp = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);

	fFitLinearBck = NULL;
	fFitLinearBck = new TF1("Linear","[0]+[1]*x",fMesonFitRange[0],fMesonFitRange[1]);


	fFitReco->SetParameter(0,mesonAmplitude);
	fFitReco->SetParameter(1,fMesonMassExpect);
	fFitReco->SetParameter(2,fMesonWidthExpect);
	if(vary){
		fFitReco->SetParameter(3,fMesonLambdaTail);
	} else {
		fFitReco->FixParameter(3,fMesonLambdaTail);
	}
	fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
	fFitReco->SetParLimits(1,fMesonMassRange[0],fMesonMassRange[1]);
	fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);
	if(vary){fFitReco->SetParLimits(3,fMesonLambdaTailRange[0],fMesonLambdaTailRange[1]);}

	fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
	fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");

	fFitReco->SetLineColor(3);
	fFitReco->SetLineWidth(1);
	fFitReco->SetLineStyle(1);

	if (vary) fMesonLambdaTail = fFitReco->GetParameter(3);

	fFitGausExp->SetParameter(0,fFitReco->GetParameter(0));
	fFitGausExp->SetParameter(1,fFitReco->GetParameter(1));
	fFitGausExp->SetParameter(2,fFitReco->GetParameter(2));
	fFitGausExp->SetParameter(3,fFitReco->GetParameter(3));

	fFitGausExp->SetParError(0,fFitReco->GetParError(0));
	fFitGausExp->SetParError(1,fFitReco->GetParError(1));
	fFitGausExp->SetParError(2,fFitReco->GetParError(2));
	fFitGausExp->SetParError(3,fFitReco->GetParError(3));

	fFitLinearBck->SetParameter(0,fFitReco->GetParameter(4));
	fFitLinearBck->SetParameter(1,fFitReco->GetParameter(5));

	fFitLinearBck->SetParError(0,fFitReco->GetParError(4));
	fFitLinearBck->SetParError(1,fFitReco->GetParError(5));

	Int_t binCenterStart;
	Double_t startBinEdge;
	Int_t binCenterEnd;
	Double_t endBinEdge;

	TVirtualFitter * fitter = TVirtualFitter::GetFitter();

	fIntLinearBck = 0;
	fIntLinearBckError = 0;
	if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("PROBLEMS") == 0){
		binCenterStart = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeFit[0]-(fMesonMassExpect-fFitReco->GetParameter(1)));
		startBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterStart)- 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
		binCenterEnd = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeFit[1]-(fMesonMassExpect-fFitReco->GetParameter(1)));
		endBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterEnd)+ 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);

		Int_t nFreePar = fFitReco->GetNumberFreeParameters();
		double * covMatrix = fitter->GetCovarianceMatrix();

		Float_t intLinearBack = fFitLinearBck->GetParameter(0)*(endBinEdge-startBinEdge)+
			0.5*fFitLinearBck->GetParameter(1)*(endBinEdge*endBinEdge-startBinEdge*startBinEdge);

		Float_t errorLinearBck = pow((pow( (endBinEdge-startBinEdge)*fFitReco->GetParError(4),2)+pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fFitReco->GetParError(5),2)+2*covMatrix[nFreePar*nFreePar-2]*(endBinEdge-startBinEdge)*0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5);

		fFileDataLog << "Parameter for bin " << ptBin << endl;
		fFileDataLog << "Gausexp: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<< "\t "<< fFitReco->GetParameter(3) <<"+-" << fFitReco->GetParError(3)<<endl;
		fFileDataLog << "Linear: \t"<<fFitReco->GetParameter(4)<<"+-" << fFitReco->GetParError(4) << "\t "<<fFitReco->GetParameter(5) <<"+-" << fFitReco->GetParError(5)<< endl;

		fIntLinearBck = intLinearBack/fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
		fIntLinearBckError = errorLinearBck/fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
	} else {
		fFileErrLog << "Fitting failed in " << ptBin << " with status " << gMinuit->fCstatu.Data() <<endl << endl;
	}
	fFitReco->DrawCopy("same");
}

void FitTrueInvMassInPtBins(TH1D* fHistoMappingSignalInvMassPtBinSingle, Double_t* fMesonIntRangeFit, Int_t ptBin)
{
	//    cout<<"Start Fitting spectra"<<endl;
	fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassRange[0],fMesonMassRange[1]);
	Double_t mesonAmplitude =fHistoMappingSignalInvMassPtBinSingle->GetMaximum();
	Double_t mesonAmplitudeMin;
	Double_t mesonAmplitudeMax;
	mesonAmplitudeMin = mesonAmplitude*98./100.;
	mesonAmplitudeMax = mesonAmplitude*115./100.;

	fFitReco = NULL;
	fFitReco = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);

	fFitReco->SetParameter(0,mesonAmplitude);
	fFitReco->SetParameter(1,fMesonMassExpect);
	fFitReco->SetParameter(2,fMesonWidthExpect);
	fFitReco->SetParameter(3,fMesonLambdaTail);

	fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
	fFitReco->SetParLimits(1,fMesonMassRange[0],fMesonMassRange[1]);
	fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);
	fFitReco->SetParLimits(3,fMesonLambdaTailRange[0],fMesonLambdaTailRange[1]);

	fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
	fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");

	fFitReco->SetLineColor(3);
	fFitReco->SetLineWidth(1);
	fFitReco->SetLineStyle(1);

	//    Int_t binCenterStart = 0;
	//    Double_t startBinEdge = 0;;
	//    Int_t binCenterEnd = 0;
	//    Double_t endBinEdge = 0;

	if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 ){
	//       binCenterStart = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeFit[0]-(fMesonMassExpect-fFitReco->GetParameter(1)));
	//       startBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterStart)- 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
	//       binCenterEnd = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeFit[1]-(fMesonMassExpect-fFitReco->GetParameter(1)));
	//       endBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterEnd)+ 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);

		fFileDataLog << "Parameter for bin " << ptBin << endl;
		fFileDataLog << "Gausexp: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<< "\t "<< fFitReco->GetParameter(3) <<"+-" << fFitReco->GetParError(3)<<endl;
	} else {
		fFileErrLog << "Fitting failed in " << ptBin << " with status " << gMinuit->fCstatu.Data() <<endl << endl;
	}
	fFitReco->DrawCopy("same");
	if (fMesonIntRangeFit){}
}


void FitWithPol2ForBG(TH1D* fHistoMappingSignalInvMassPtBinSingle,Double_t * fMesonIntRangeFit, Int_t ptBin, Bool_t vary)
{
	//    cout<<"Start Fitting spectra"<<endl;
	Int_t startBinSearch = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeFit[0]);
	Int_t endBinSearch = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeFit[1]);
	Double_t mesonAmplitude = 0;
	for (Int_t i = startBinSearch; i < endBinSearch+1; i++){
		if (mesonAmplitude < fHistoMappingSignalInvMassPtBinSingle->GetBinContent(i)){
			mesonAmplitude = fHistoMappingSignalInvMassPtBinSingle->GetBinContent(i);
		}
	}
	mesonAmplitude = mesonAmplitude-fHistoMappingSignalInvMassPtBinSingle->GetBinContent(endBinSearch);

	fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassRange[0] ,fMesonMassRange[1]);
	Double_t mesonAmplitudeMin = mesonAmplitude*70./100.;
	Double_t mesonAmplitudeMax = mesonAmplitude*110./100.;
	Double_t fitRange[2] = {0.3,0.8};

	fFitReco = NULL;
	fFitReco = new TF1("GaussExpLinear","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x+[6]*x*x)+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x+[6]*x*x)",fitRange[0] ,fitRange[1]);

	fFitGausExp = NULL;
	fFitGausExp = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",0.05,0.3);

	fFitLinearBck = NULL;
	fFitLinearBck = new TF1("Linear","[0]+[1]*x+[2]*x*x",fitRange[0] ,fitRange[1]);

	fFitReco->SetParameter(0,mesonAmplitude);
	fFitReco->SetParameter(1,fMesonMassExpect);
	fFitReco->SetParameter(2,fMesonWidthExpect);
	if(vary){
		fFitReco->SetParameter(3,fMesonLambdaTail);
	} else {
		fFitReco->FixParameter(3,fMesonLambdaTail);
	}
	fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
	fFitReco->SetParLimits(1,fMesonMassExpect*95./100.,fMesonMassExpect*105./100.);
	fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);
	if(vary){fFitReco->SetParLimits(3,fMesonLambdaTailRange[0],fMesonLambdaTailRange[1]);}

	fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
	fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");

	fFitReco->SetLineColor(3);
	fFitReco->SetLineWidth(1);
	fFitReco->SetLineStyle(1);

	//if (vary) fMesonLambdaTail = fFitReco->GetParameter(3);

	fFitGausExp->SetParameter(0,fFitReco->GetParameter(0));
	fFitGausExp->SetParameter(1,fFitReco->GetParameter(1));
	fFitGausExp->SetParameter(2,fFitReco->GetParameter(2));
	fFitGausExp->SetParameter(3,fFitReco->GetParameter(3));

	fFitGausExp->SetParError(0,fFitReco->GetParError(0));
	fFitGausExp->SetParError(1,fFitReco->GetParError(1));
	fFitGausExp->SetParError(2,fFitReco->GetParError(2));
	fFitGausExp->SetParError(3,fFitReco->GetParError(3));

	fFitLinearBck->SetParameter(0,fFitReco->GetParameter(4));
	fFitLinearBck->SetParameter(1,fFitReco->GetParameter(5));
	fFitLinearBck->SetParameter(2,fFitReco->GetParameter(6));

	fFitLinearBck->SetParError(0,fFitReco->GetParError(4));
	fFitLinearBck->SetParError(1,fFitReco->GetParError(5));
	fFitLinearBck->SetParError(2,fFitReco->GetParError(6));

	if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 ){
		fFileDataLog << "Parameter for bin " << ptBin << endl;
		fFileDataLog << "Gausexp: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<< "\t "<< fFitReco->GetParameter(3) <<"+-" << fFitReco->GetParError(3)<<endl;
		fFileDataLog << "Quadratic: \t"<<fFitReco->GetParameter(4)<<"+-" << fFitReco->GetParError(4) << "\t "<<fFitReco->GetParameter(5) <<"+-" << fFitReco->GetParError(5) << "\t "<<fFitReco->GetParameter(6) <<"+-" << fFitReco->GetParError(6)<< endl;
	} else {
		fFileErrLog << "Fitting failed in " << ptBin << " with status " << gMinuit->fCstatu.Data() <<endl << endl;
	}
	fFitReco->DrawCopy("same");
}




void IntegrateHistoInvMass(TH1D * fHistoMappingSignalInvMassPtBinSingle, Double_t * fMesonIntRangeInt)
{
	Int_t binLowMassMeson = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[0]);
	Int_t binHighMassMeson = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[1]);
	fYields = fHistoMappingSignalInvMassPtBinSingle->IntegralAndError(binLowMassMeson,binHighMassMeson,fYieldsError);
}


void IntegrateHistoInvMassStream(TH1D * fHistoMappingSignalInvMassPtBinSingle, Double_t * fMesonIntRangeInt)
{
	Int_t binLowMassMeson = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[0]);
	Int_t binHighMassMeson = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[1]);
	fYields = fHistoMappingSignalInvMassPtBinSingle->IntegralAndError(binLowMassMeson,binHighMassMeson,fYieldsError);
	for ( Int_t M = binLowMassMeson; M < binHighMassMeson+1; M++){
		fFileDataLog << M << "\t" << fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(M) <<"\t" <<fHistoMappingSignalInvMassPtBinSingle->GetBinContent(M)<< "+-"<< fHistoMappingSignalInvMassPtBinSingle->GetBinError(M)<< endl;
	}
}


void IntegrateFitFunc(TF1 * fFunc, TH1D *  fHistoMappingSignalInvMassPtBinSingle,Double_t * fMesonIntRangeInt)
{
	fYieldsFunc = fFunc->Integral(fMesonIntRangeInt[0],fMesonIntRangeInt[1])/fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
}



void FillHistosArrayMC(TH1D* fHistoMCMesonPtWithinAcceptanceFill, TH1D * fHistoMCMesonPtFill, TH1D * fDeltaPtFill, TString nameCanvas)
{
//    Char_t nameHisto[100] = "fHistoMCMesonPtEtaWithinAcceptance";
	fHistoMCMesonPtWithinAcceptanceFill->Sumw2();
	fHistoMCMesonWithinAccepPt = (TH1D*)fHistoMCMesonPtWithinAcceptanceFill->Rebin(fNBinsPt,"",fBinsPt); // Proper bins in Pt
	fHistoMCMesonWithinAccepPt->Divide(fDeltaPtFill);
	fHistoMCMesonPtFill->Sumw2();
	fHistoMCMesonPt1 = (TH1D*)fHistoMCMesonPtFill->Rebin(fNBinsPt,"",fBinsPt); // Proper bins in Pt
	fHistoMCMesonPt1->Divide(fDeltaPtFill);
	if (nameCanvas){}	
}



void CalculateMesonAcceptance()
{
	fHistoMCMesonAcceptPt = new TH1D("fMCMesonAccepPt","",fNBinsPt,fBinsPt);
	fHistoMCMesonAcceptPt->Sumw2();

	fHistoMCMesonAcceptPt->Divide(fHistoMCMesonWithinAccepPt,fHistoMCMesonPt1,1.,1.,"B");
	fHistoMCMesonAcceptPt->DrawCopy();
	fFileDataLog << endl << "Calculation of the Acceptance" << endl;
	for ( Int_t i = 1; i < fHistoMCMesonAcceptPt->GetNbinsX()+1 ; i++){
		fFileDataLog << "Bin " << i << "\t"<< fHistoMCMesonAcceptPt->GetBinCenter(i)<< "\t" << fHistoMCMesonAcceptPt->GetBinContent(i) << "\t" << fHistoMCMesonAcceptPt->GetBinError(i) <<endl;
	}
}

void CalculateMesonEfficiency(TH1D* fMC_fMesonYieldsPt, TH1D* fMC_SecondaryYieldPt, TString nameEfi, TString nameCanvas )
{
	fHistoMCMesonEffiPt = new TH1D(nameEfi.Data(),"",fNBinsPt,fBinsPt);	
	fHistoMCMesonEffiPt->Sumw2();
	fHistoMCMesonEffiPt->Add(fMC_fMesonYieldsPt,1.);
	fHistoMCMesonEffiPt->Add(fMC_SecondaryYieldPt,-1.);
		
		
	TH1D* fHistoMCMesonEffiScaledByPt = (TH1D*)fHistoMCMesonEffiPt->Clone("ScaledByPt");
	for (Int_t i = 1; i < fHistoMCMesonEffiScaledByPt->GetNbinsX()+1; i++){
		fHistoMCMesonEffiScaledByPt->SetBinContent(i, fHistoMCMesonEffiScaledByPt->GetBinContent(i)/fNEvents);//fHistoMCMesonEffiScaledByPt->GetBinCenter(i)/2/TMath::Pi()
		fHistoMCMesonEffiScaledByPt->SetBinError(i, fHistoMCMesonEffiScaledByPt->GetBinError(i)/fNEvents	);//fHistoMCMesonEffiScaledByPt->GetBinCenter(i)/2/TMath::Pi()
	}
		
	fHistoMCMesonEffiPt->Divide(fHistoMCMesonEffiPt,fHistoMCMesonWithinAccepPt,1.,1.,"B");
	fFileDataLog << endl << "Calculation of the Efficiency" << nameEfi.Data()<< endl;
	for ( Int_t i = 1; i < fHistoMCMesonEffiPt->GetNbinsX()+1 ; i++){
		fFileDataLog << "Bin " << i << "\t" << fHistoMCMesonEffiPt->GetBinCenter(i)<< "\t"<< fHistoMCMesonEffiPt->GetBinContent(i) << "\t" << fHistoMCMesonEffiPt->GetBinError(i) <<endl;
	}
	fHistoMCMesonEffiFitPt = (TH1D*)fHistoMCMesonEffiPt->Clone(Form("%s_Fit",nameEfi.Data()));
	if (nameCanvas){} 
}



void SaveHistos(Int_t optionMC, TString fCutID, TString fPrefix3)
{
	const char* nameOutput = Form("%s/%s/%s_%s_GammaConvEtaInPiPlPiMiGamma%s_%s.root",fCutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPrefix3.Data(),fPeriodFlag.Data(),fCutID.Data());
	fOutput1 = new TFile(nameOutput,"RECREATE");

	cout << "Begin writing Uncorrected File" << endl;
	
	fHistoYieldMeson->Write();
	fHistoYieldMesonPerEvent->Write();
	fHistoYieldMesonBackFit->Write();
	fHistoYieldMesonPerEventBackFit->Write();
	fHistoSignMeson->Write();
	fHistoSBMeson->Write();

	fHistoYieldMesonNarrow->Write();
	fHistoYieldMesonPerEventNarrow->Write();
	fHistoSignMesonNarrow->Write();
	fHistoSBMesonNarrow->Write();

	fHistoYieldMesonWide->Write();
	fHistoYieldMesonPerEventWide->Write();
	fHistoSignMesonWide->Write();
	fHistoSBMesonWide->Write();

	fHistoMassMeson->Write();
	fHistoWidthMeson->Write();
	fHistoFWHMMeson->Write();
	fDeltaPt->Write();
	
	fHistoYieldMesonLeft->Write();
	fHistoYieldMesonLeftPerEvent->Write();
	fHistoSignMesonLeft->Write();
	fHistoSBMesonLeft->Write();

	fHistoYieldMesonLeftNarrow->Write();
	fHistoYieldMesonLeftPerEventNarrow->Write();
	fHistoSignMesonLeftNarrow->Write();
	fHistoSBMesonLeftNarrow->Write();

	fHistoYieldMesonLeftWide->Write();
	fHistoYieldMesonLeftPerEventWide->Write();
	fHistoSignMesonLeftWide->Write();
	fHistoSBMesonLeftWide->Write();

	fHistoMassMesonLeft->Write();
	fHistoWidthMesonLeft->Write();
	fHistoFWHMMesonLeft->Write();
	fMesonFullPtSignal->Write();
	fMesonFullPtBackground->Write();
	fMesonFullPtBackNorm->SetName("Mapping_BackNorm_InvMass_FullPt");
	fMesonFullPtBackNorm->Write();
	fNumberOfGoodESDTracks->Write();
	fEventQuality->Write();

	TString nameHistoSignal;
	TString nameHistoPeakPos;
	TString nameHistoBckNorm;
	TString fitnameSignal;
	TString nameHistoSignalPos;
	for(Int_t ii =fStartPtBin;ii<fNBinsPt;ii++){
		fHistoWeightsBGZbinVsMbin[ii]->Write(Form("BGWeights_%02d", ii));
		fHistoFillPerEventBGZbinVsMbin[ii]->Scale(1./fNEvents);
		fHistoFillPerEventBGZbinVsMbin[ii]->Write(Form("BGPoolsFillstatus_%02d", ii));
	
		fHistoMappingPiPlPiMiGammaInvMassPtBin[ii]->Write();
		nameHistoBckNorm = Form("Mapping_BckNorm_InvMass_in_Pt_Bin%02d", ii);
		fHistoMappingBackNormInvMassPtBin[ii]->Write(nameHistoBckNorm.Data());
		nameHistoSignal = Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d", ii);
		fHistoMappingSignalInvMassPtBin[ii]->Write(nameHistoSignal.Data());
		fitnameSignal = Form("Signal_InvMassFit_in_Pt_Bin%02d", ii);
		if(fFitSignalInvMassPtBin[ii]!=0x00) fFitSignalInvMassPtBin[ii]->Write(fitnameSignal.Data());
	}

	if(optionMC){
		fHistoTrueSignMeson->Write();
		fHistoTrueSBMeson->Write();
		fHistoMCMesonPtWithinAcceptance->Write();
		fHistoMCMesonWithinAccepPt->Write(); // Proper bins in Pt
		fHistoMCMesonPt1->Write(); // Proper bins in Pt
		fHistoYieldTrueMeson->Write();
		fHistoYieldTrueMesonWide->Write();
		fHistoYieldTrueMesonNarrow->Write();
		
		fHistoTrueMesonInvMassVSPt->Write();
		fHistoYieldTrueGGMeson->Write();
		fHistoYieldTrueDalitzMeson->Write();
		fHistoYieldTrueGGMesonWide->Write();
		fHistoYieldTrueDalitzMesonWide->Write();
		fHistoYieldTrueGGMesonNarrow->Write();
		fHistoYieldTrueDalitzMesonNarrow->Write();
		for(Int_t ii =fStartPtBin;ii<fNBinsPt;ii++){
			fHistoMappingTrueMesonInvMassPtBins[ii]->Write();
			fHistoMappingTrueGGMesonInvMassPtBins[ii]->Write();
			fHistoMappingTrueDalitzInvMassPtBins[ii]->Write();	
			if (fFitTrueSignalInvMassPtBin[ii]!=0x00) fFitTrueSignalInvMassPtBin[ii]->Write();
		}
	}

	cout << "End writing Uncorrected File" << endl;
	
	fOutput1->Write();
	fOutput1->Close();
}

void SaveCorrectionHistos(TString fCutID, TString fPrefix3)
{
	const char* nameOutput = Form("%s/%s/%s_%s_GammaConvV1CorrectionHistos%s_%s.root",fCutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPrefix3.Data(),fPeriodFlag.Data(),fCutID.Data());
	fOutput2 = new TFile(nameOutput,"RECREATE");
	cout<<"======================================================"<<endl;
	cout<<"======================================================"<<endl;
	cout<<"======================================================"<<endl;
	cout<<nameOutput<<endl;
	cout<<"======================================================"<<endl;
	cout<<"======================================================"<<endl;
	cout<<"======================================================"<<endl;

	cout << "Begin writing Correction File" << endl;
	fHistoMCMesonAcceptPt->Write();
	fHistoTrueMesonInvMassVSPt->Write();
	fHistoMonteMesonEffiPt->Write();
	fHistoMonteMesonEffiBackFitPt->Write();
	fHistoMonteMesonNarrowEffiPt->Write();
	fHistoMonteMesonWideEffiPt->Write();
	fHistoMonteMesonLeftEffiPt->Write();
	fHistoMonteMesonLeftNarrowEffiPt->Write();
	fHistoMonteMesonLeftWideEffiPt->Write();
	fHistoMCTrueMesonEffiPt->Write();
	fHistoMCTrueMesonNarrowEffiPt->Write();
	fHistoMCTrueMesonWideEffiPt->Write();
	fHistoMonteMesonEffiFitPt->Write();
	fHistoMonteMesonNarrowEffiFitPt->Write();
	fHistoMonteMesonWideEffiFitPt->Write();
	fHistoMonteMesonLeftEffiFitPt->Write();
	fHistoMonteMesonLeftNarrowEffiFitPt->Write();
	fHistoMonteMesonLeftWideEffiFitPt->Write();
	fHistoMCTrueMesonEffiFitPt->Write();
	fHistoMCTrueMesonNarrowEffiFitPt->Write();
	fHistoMCTrueMesonWideEffiFitPt->Write();
	
	fHistoTrueMassMeson->Write();
	fHistoTrueFWHMMeson->Write();
	fHistoMCMesonPt1->SetName("MC_Meson_genPt");
	fHistoMCMesonPt1->Write(); // Proper bins in Pt
	fHistoMCMesonPt->SetName("MC_Meson_genPt_oldBin");
	fHistoMCMesonPt->Scale(1./fHistoMCMesonPt->GetBinWidth(5));
	//    fHistoMCMesonPt->GetXaxis()->SetRangeUser(0.,25.);
	fHistoMCMesonPt->Write(); // Proper bins in Pt
	if (fHistoMCMesonPtWOWeights){
		fHistoMCMesonPtWOWeights->SetName("MC_Meson_genPt_WOWeights");
		fHistoMCMesonPtWOWeights->Scale(1./fHistoMCMesonPtWOWeights->GetBinWidth(5));
	//       fHistoMCMesonPtWOWeights->GetXaxis()->SetRangeUser(0.,25.);
		fHistoMCMesonPtWOWeights->Write(); // Proper bins in Pt
	}
	if (fHistoMCMesonPtWeights){
		fHistoMCMesonPtWeights->Write("MC_Meson_genPt_Weights"); // Proper bins in Pt
	}   
	fEventQuality->Write();
	fHistoYieldTrueGGMeson->Write();
	fHistoYieldTrueDalitzMeson->Write();
	fHistoYieldTrueGGMesonWide->Write();
	fHistoYieldTrueDalitzMesonWide->Write();
	fHistoYieldTrueGGMesonNarrow->Write();
	fHistoYieldTrueDalitzMesonNarrow->Write();
	fHistoYieldTrueGGCont->Write();
	fHistoYieldTrueDalitzCont->Write();
	fHistoYieldTrueGGContWide->Write();
	fHistoYieldTrueDalitzContWide->Write();
	fHistoYieldTrueGGContNarrow->Write();
	fHistoYieldTrueDalitzContNarrow->Write();

	fHistoYieldTrueMeson->Write();
	fHistoYieldTrueMesonWide->Write();
	fHistoYieldTrueMesonNarrow->Write();
	
	cout << "end writing Correction File" << endl;
	
	fOutput2->Write();
	fOutput2->Close();
}

void Initialize(TString setPi0,Int_t numberOfBins){
	if (setPi0.CompareTo("Eta") == 0){
		fNBinsPt = 		numberOfBins;
		fBinsPt= 			new Double_t[20];
		fNRebin = 		new Int_t[19];

		if (fEnergyFlag.CompareTo("7TeV") == 0) {
			fStartPtBin=	5;
			fColumn = 	4;
			fRow = 		3;
			
			if (fNBinsPt > 12) {
				cout << "You have chosen to have more than 12 bins for Eta, this is not possible, it will be reduced to 12" << endl;
				fNBinsPt = 12;
			}
			for (Int_t i = 0; i < fNBinsPt+2; i++) {
				fBinsPt[i] = fBinsEta7TeVPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsEta7TeVPtRebin[i];
			}
			fExampleBin = 6;
		} else if (fEnergyFlag.CompareTo("2.76TeV") == 0) {
			fStartPtBin=	1;
			fColumn = 	3;
			fRow = 		3;

			if (fNBinsPt > 7) {
				cout << "You have chosen to have more than 7 bins for Eta, this is not possible, it will be reduced to 7" << endl;
				fNBinsPt = 7;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] = fBinsEta2760GeVPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsEta2760GeVPtRebin[i];
			}
			fExampleBin = 4;
		} else if (fEnergyFlag.CompareTo("900GeV") == 0) {
			fStartPtBin=	1;
			fColumn = 	2;
			fRow = 		2;

			if (fNBinsPt > 3) {
				cout << "You have chosen to have more than 3 bins for Eta, this is not possible, it will be reduced to 3" << endl;
				fNBinsPt = 3;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] = fBinsEta900GeVPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsEta900GeVPtRebin[i];
			}
			fExampleBin = 2;
		} else if( fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0) { 
	//          fStartPtBin=	1;
	//          fColumn = 	2;
	//          fRow = 		2;
	// 
				TString Cent = fCutSelection(0,3);
				if (Cent.CompareTo("601")==0 ||Cent.CompareTo("612")==0 || Cent.CompareTo("512")==0  ){  //min bias binning           
				fStartPtBin=   1;
				fColumn =   4;
				fRow =      3;
				if (fNBinsPt > 12) {
					cout << "You have chosen to have more than 12 bins, this is not possible, it will be reduced to 12" << endl;
					fNBinsPt = 12;
				}
				for (Int_t i = 0; i < fNBinsPt+1; i++) {
					fBinsPt[i] = fBinsEtaHIPtLHC11h[i];
					if (i < fNBinsPt+1) fNRebin[i] = fBinsEtaHIPtRebinLHC11hFinerBinning[i];
				}
				} else {
				fStartPtBin=   1;
				fColumn =   4;
				fRow =      3;
				if (fNBinsPt > 12) {
					cout << "You have chosen to have more than 12 bins, this is not possible, it will be reduced to 12" << endl;
					fNBinsPt = 12;
				}
				for (Int_t i = 0; i < fNBinsPt+1; i++) {
					fBinsPt[i] = fBinsEtaHIPtLHC11h[i];
					if (i < fNBinsPt+1) fNRebin[i] = fBinsEtaHIPtRebinLHC11h[i];
				}
				}         
			fExampleBin = 4;
		} else if( fEnergyFlag.CompareTo("pPb_5.023TeV") == 0) { 
			fStartPtBin=   3;
			fColumn =   5;
			fRow =      3;

			if (fNBinsPt > 14) {
				cout << "You have chosen to have more than 14 bins, this is not possible, it will be reduced to 14" << endl;
				fNBinsPt = 14;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] = fBinsEtapPbPt3Body[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsEtapPbPt3BodyRebin[i];
			}

			fExampleBin = 6;
		}
		fPeakRange = 		new Double_t[2]; 	fPeakRange[0]=	0.5; 	fPeakRange[1]=	0.58; //eta 0.9
		fFitRange = 		new Double_t[2]; 	fFitRange[0]=	0.4; 	fFitRange[1]=	0.7; //eta 0.9
		fBGFitRange = 		new Double_t[2]; 	fBGFitRange[0]=	0.58; 	fBGFitRange[1]=	0.7; //eta 0.9
		fBGFitRangeLeft = 	new Double_t[2]; 	fBGFitRangeLeft[0]=	0.4; 	fBGFitRangeLeft[1]=	0.48;  // eta 09
		fMesonPlotRange = 	new Double_t[2]; 	fMesonPlotRange[0]=	0.53; 	fMesonPlotRange[1]=	0.560;
		fMesonIntRange = 	new Double_t[2]; 	fMesonIntRange[0] = 0.50; 	fMesonIntRange[1]=	0.57;
		fMesonIntRangeWide = new Double_t[2]; 	fMesonIntRangeWide[0]=0.48; 	fMesonIntRangeWide[1]=0.58;
		fMesonIntRangeNarrow = new Double_t[2]; fMesonIntRangeNarrow[0]=0.515; fMesonIntRangeNarrow[1]=0.56;
		fMesonMassRange = 	new Double_t[2]; 	fMesonMassRange[0]=	0.4; 	fMesonMassRange[1]=	0.7;;
		fMesonFitRange = 	new Double_t[2]; 	fMesonFitRange[0]=	0.4; 	fMesonFitRange[1]=	0.7;
		//		fMesonFitRangeWithoutPeak = new Double_t[2]; fMesonFitRangeWithoutPeak[0] = 0.4; fMesonFitRangeWithoutPeak[0] = 0.7;
		fMesonId=			221;
		if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0){
			fMesonWidthExpect = 0.010;
		} else {
			fMesonWidthExpect = 0.005;
		}   
		fMesonLambdaTail = 	0.007;
		fMesonWidthRange = 	new Double_t[2]; 	 fMesonWidthRange[0]=0.002;     fMesonWidthRange[1]=0.020;
		fMesonLambdaTailRange = new Double_t[2]; fMesonLambdaTailRange[0]=0.0005; fMesonLambdaTailRange[1]=0.026;
		fMidPt = 							new Double_t[2]; fMidPt[0] = 1.2; fMidPt[1] = 2.5;
	} else {
		 cout <<"ERROR no meson specified" << endl;
		return;
	}	
	fFullPt = 						new Double_t[2]; fFullPt[0] = 0.4; fFullPt[1] = 15.;

	fMesonCurIntRange = 				new Double_t[2];
	fMesonCurIntRangeBackFit = 				new Double_t[2];
	fMesonCurIntRangeWide = 				new Double_t[2];
	fMesonCurIntRangeNarrow = 			new Double_t[2];
	fMesonCurLeftIntRange = 				new Double_t[2];
	fMesonCurLeftIntRangeWide = 			new Double_t[2];
	fMesonCurLeftIntRangeNarrow = 		new Double_t[2];
	fMesonTrueIntRange = 				new Double_t[2];
	fMesonTrueIntRangeWide = 			new Double_t[2];
	fMesonTrueIntRangeNarrow = 			new Double_t[2];
	
	fPiPlPiMiGammaYields = 						new Double_t[fNBinsPt];
	fMesonYieldsBackFit =					new Double_t[fNBinsPt];
	fBckYields = 						new Double_t[fNBinsPt];
	fTotalBckYields = 					new Double_t[fNBinsPt];
	fMesonYields = 					new Double_t[fNBinsPt];
	fMesonTrueYields = 					new Double_t[fNBinsPt];
	fMesonTrueGGContYields = 					new Double_t[fNBinsPt];
	fMesonTrueDalitzContYields = 					new Double_t[fNBinsPt];
	fMesonYieldsFunc = 					new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFunc = 		new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncBackFit = 		new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFunc = 		new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncBackFit = 		new Double_t[fNBinsPt];
	fMesonYieldsPerEvent = 				new Double_t[fNBinsPt];
	fMesonYieldsPerEventBackFit = 				new Double_t[fNBinsPt];
	fMesonMass = 						new Double_t[fNBinsPt];
	fMesonMassBackFit = 						new Double_t[fNBinsPt];
	fMesonWidth = 						new Double_t[fNBinsPt];
	fMesonWidthBackFit = 						new Double_t[fNBinsPt];
	fMesonSB = 						new Double_t[fNBinsPt];
	fMesonSign = 						new Double_t[fNBinsPt];
	fMesonFWHM = 						new Double_t[fNBinsPt];
	fMesonFWHMAlpha01 =					new Double_t[fNBinsPt];
	fMesonTrueMass = 					new Double_t[fNBinsPt];
	fMesonTrueFWHM = 					new Double_t[fNBinsPt];
	fMesonTrueSB = 					new Double_t[fNBinsPt];
	fMesonTrueSign = 					new Double_t[fNBinsPt];


	// Normalization at the left of the peak
	fPiPlPiMiGammaYieldsLeft = 					new Double_t[fNBinsPt];
	fBckYieldsLeft = 					new Double_t[fNBinsPt];
	fTotalBckYieldsLeft = 				new Double_t[fNBinsPt];
	fMesonYieldsLeft = 					new Double_t[fNBinsPt];
	fMesonYieldsFuncLeft = 				new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncLeft = 		new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncLeft = 	new Double_t[fNBinsPt];
	fMesonYieldsLeftPerEvent = 			new Double_t[fNBinsPt];
	fMesonMassLeft = 					new Double_t[fNBinsPt];
	fMesonWidthLeft = 					new Double_t[fNBinsPt];
	fMesonSBLeft = 					new Double_t[fNBinsPt];
	fMesonSignLeft = 					new Double_t[fNBinsPt];
	fMesonFWHMLeft = 					new Double_t[fNBinsPt];

	// Narrow Integration Window
	fPiPlPiMiGammaYieldsNarrow = 					new Double_t[fNBinsPt];
	fBckYieldsNarrow = 					new Double_t[fNBinsPt];
	fTotalBckYieldsNarrow =	 			new Double_t[fNBinsPt];
	fMesonYieldsNarrow = 				new Double_t[fNBinsPt];
	fMesonTrueYieldsNarrow = 			new Double_t[fNBinsPt];
	fMesonTrueGGContYieldsNarrow = 					new Double_t[fNBinsPt];
	fMesonTrueDalitzContYieldsNarrow = 					new Double_t[fNBinsPt];
	fMesonYieldsFuncNarrow = 			new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncNarrow = 	new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncNarrow = 	new Double_t[fNBinsPt];
	fMesonYieldsPerEventNarrow = 			new Double_t[fNBinsPt];
	fMesonSBNarrow = 					new Double_t[fNBinsPt];
	fMesonSignNarrow = 					new Double_t[fNBinsPt];

	fPiPlPiMiGammaYieldsLeftNarrow = 				new Double_t[fNBinsPt];
	fBckYieldsLeftNarrow = 				new Double_t[fNBinsPt];
	fTotalBckYieldsLeftNarrow = 			new Double_t[fNBinsPt];
	fMesonYieldsLeftNarrow = 			new Double_t[fNBinsPt];
	fMesonYieldsFuncLeftNarrow = 			new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncLeftNarrow = new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncLeftNarrow = new Double_t[fNBinsPt];
	fMesonYieldsLeftPerEventNarrow = 		new Double_t[fNBinsPt];
	fMesonSBLeftNarrow = 				new Double_t[fNBinsPt];
	fMesonSignLeftNarrow = 				new Double_t[fNBinsPt];

	// Wide Integration Window
	fPiPlPiMiGammaYieldsWide = 					new Double_t[fNBinsPt];
	fBckYieldsWide = 					new Double_t[fNBinsPt];
	fTotalBckYieldsWide =	 			new Double_t[fNBinsPt];
	fMesonYieldsWide = 					new Double_t[fNBinsPt];
	fMesonTrueYieldsWide = 				new Double_t[fNBinsPt];
	fMesonTrueGGContYieldsWide = 					new Double_t[fNBinsPt];
	fMesonTrueDalitzContYieldsWide = 					new Double_t[fNBinsPt];
	fMesonYieldsFuncWide = 				new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncWide = 		new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncWide =	new Double_t[fNBinsPt];
	fMesonYieldsPerEventWide = 			new Double_t[fNBinsPt];
	fMesonSBWide = 					new Double_t[fNBinsPt];
	fMesonSignWide = 					new Double_t[fNBinsPt];

	fPiPlPiMiGammaYieldsLeftWide = 				new Double_t[fNBinsPt];
	fBckYieldsLeftWide = 				new Double_t[fNBinsPt];
	fTotalBckYieldsLeftWide = 			new Double_t[fNBinsPt];
	fMesonYieldsLeftWide = 				new Double_t[fNBinsPt];
	fMesonYieldsFuncLeftWide = 			new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncLeftWide =	new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncLeftWide = new Double_t[fNBinsPt];
	fMesonYieldsLeftPerEventWide = 		new Double_t[fNBinsPt];
	fMesonSBLeftWide = 					new Double_t[fNBinsPt];
	fMesonSignLeftWide = 				new Double_t[fNBinsPt];

	fPiPlPiMiGammaYieldsError = 					new Double_t[fNBinsPt];
	fMesonYieldsBackFitError =				new Double_t[fNBinsPt];
	fBckYieldsError = 					new Double_t[fNBinsPt];
	fTotalBckYieldsError = 				new Double_t[fNBinsPt];
	fMesonYieldsError = 				new Double_t[fNBinsPt];
	fMesonTrueYieldsError = 				new Double_t[fNBinsPt];
	fMesonTrueGGContYieldsError = 					new Double_t[fNBinsPt];
	fMesonTrueDalitzContYieldsError = 					new Double_t[fNBinsPt];
	fMesonYieldsFuncError = 				new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncError = 	new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncBackFitError = 	new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncError =	new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncBackFitError =new Double_t[fNBinsPt];
	fMesonYieldsPerEventError = 			new Double_t[fNBinsPt];
	fMesonYieldsPerEventBackFitError = 			new Double_t[fNBinsPt];
	fMesonMassError = 					new Double_t[fNBinsPt];
	fMesonMassBackFitError = 					new Double_t[fNBinsPt];
	fMesonWidthError = 					new Double_t[fNBinsPt];
	fMesonWidthBackFitError = 					new Double_t[fNBinsPt];
	fMesonSBError = 					new Double_t[fNBinsPt];
	fMesonSignError = 					new Double_t[fNBinsPt];
	fMesonFWHMError = 					new Double_t[fNBinsPt];
	fMesonFWHMAlpha01Error = 				new Double_t[fNBinsPt];
	fMesonTrueMassError = 				new Double_t[fNBinsPt];
	fMesonTrueFWHMError = 				new Double_t[fNBinsPt];
	fMesonTrueSBError = 				new Double_t[fNBinsPt];
	fMesonTrueSignError = 				new Double_t[fNBinsPt];

	fPiPlPiMiGammaYieldsLeftError = 				new Double_t[fNBinsPt];
	fBckYieldsLeftError =	 			new Double_t[fNBinsPt];
	fTotalBckYieldsLeftError = 			new Double_t[fNBinsPt];
	fMesonYieldsLeftError = 				new Double_t[fNBinsPt];
	fMesonYieldsFuncLeftError = 			new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncLeftError = 	new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncLeftError = new Double_t[fNBinsPt];
	fMesonYieldsLeftPerEventError = 		new Double_t[fNBinsPt];
	fMesonMassLeftError = 				new Double_t[fNBinsPt];
	fMesonWidthLeftError = 				new Double_t[fNBinsPt];
	fMesonSBLeftError = 				new Double_t[fNBinsPt];
	fMesonSignLeftError = 				new Double_t[fNBinsPt];
	fMesonFWHMLeftError = 				new Double_t[fNBinsPt];

	// Narrow integration Window
	fPiPlPiMiGammaYieldsNarrowError = 				new Double_t[fNBinsPt];
	fBckYieldsNarrowError = 				new Double_t[fNBinsPt];
	fTotalBckYieldsNarrowError = 			new Double_t[fNBinsPt];
	fMesonYieldsNarrowError = 			new Double_t[fNBinsPt];
	fMesonTrueYieldsNarrowError = 		new Double_t[fNBinsPt];
	fMesonTrueGGContYieldsNarrowError = 					new Double_t[fNBinsPt];
	fMesonTrueDalitzContYieldsNarrowError = 					new Double_t[fNBinsPt];
	fMesonYieldsFuncNarrowError = 		new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncNarrowError = new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncNarrowError = new Double_t[fNBinsPt];
	fMesonYieldsPerEventNarrowError = 		new Double_t[fNBinsPt];
	fMesonSBNarrowError = 				new Double_t[fNBinsPt];
	fMesonSignNarrowError = 				new Double_t[fNBinsPt];

	fPiPlPiMiGammaYieldsLeftNarrowError = 			new Double_t[fNBinsPt];
	fBckYieldsLeftNarrowError = 			new Double_t[fNBinsPt];
	fTotalBckYieldsLeftNarrowError = 		new Double_t[fNBinsPt];
	fMesonYieldsLeftNarrowError = 		new Double_t[fNBinsPt];
	fMesonYieldsFuncLeftNarrowError = 		new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncLeftNarrowError = new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncLeftNarrowError = new Double_t[fNBinsPt];
	fMesonYieldsLeftPerEventNarrowError = 	new Double_t[fNBinsPt];
	fMesonSBLeftNarrowError = 			new Double_t[fNBinsPt];
	fMesonSignLeftNarrowError = 			new Double_t[fNBinsPt];

	// Wide integration Window
	fPiPlPiMiGammaYieldsWideError = 				new Double_t[fNBinsPt];
	fBckYieldsWideError = 				new Double_t[fNBinsPt];
	fTotalBckYieldsWideError = 			new Double_t[fNBinsPt];
	fMesonYieldsWideError = 				new Double_t[fNBinsPt];
	fMesonTrueYieldsWideError = 			new Double_t[fNBinsPt];
	fMesonTrueGGContYieldsWideError = 					new Double_t[fNBinsPt];
	fMesonTrueDalitzContYieldsWideError = 					new Double_t[fNBinsPt];	
	fMesonYieldsFuncWideError = 			new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncWideError = 	new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncWideError = new Double_t[fNBinsPt];
	fMesonYieldsPerEventWideError = 		new Double_t[fNBinsPt];
	fMesonSBWideError = 				new Double_t[fNBinsPt];
	fMesonSignWideError = 				new Double_t[fNBinsPt];

	fPiPlPiMiGammaYieldsLeftWideError = 			new Double_t[fNBinsPt];
	fBckYieldsLeftWideError = 			new Double_t[fNBinsPt];
	fTotalBckYieldsLeftWideError = 		new Double_t[fNBinsPt];
	fMesonYieldsLeftWideError = 			new Double_t[fNBinsPt];
	fMesonYieldsFuncLeftWideError = 		new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncLeftWideError = new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncLeftWideError = new Double_t[fNBinsPt];
	fMesonYieldsLeftPerEventWideError = 	new Double_t[fNBinsPt];
	fMesonSBLeftWideError = 				new Double_t[fNBinsPt];
	fMesonSignLeftWideError = 			new Double_t[fNBinsPt];

	fHistoMappingTrueMesonInvMassPtBins = 	new TH1D*[fNBinsPt];    
	fHistoMappingTrueGGMesonInvMassPtBins = 	new TH1D*[fNBinsPt];    
	fHistoMappingTrueDalitzInvMassPtBins = 	new TH1D*[fNBinsPt];    
		fHistoWeightsBGZbinVsMbin = new TH2F*[fNBinsPt];    
	fHistoFillPerEventBGZbinVsMbin = new TH2F*[fNBinsPt];    
	
	fHistoMappingPiPlPiMiGammaInvMassPtBin = 		new TH1D*[fNBinsPt];    
	fHistoMappingGGInvMassBackFitPtBin =		new TH1D*[fNBinsPt];
	fHistoMappingGGInvMassBackFitWithoutSignalPtBin =		new TH1D*[fNBinsPt];
	fHistoMappingBackInvMassPtBin = 		new TH1D*[fNBinsPt];
	fHistoMappingBackNormInvMassPtBin = 	new TH1D*[fNBinsPt];
	fHistoMappingSignalInvMassPtBin =		new TH1D*[fNBinsPt];
	fHistoMappingRatioSBInvMassPtBin= 		new TH1D*[fNBinsPt];


	fBackgroundFitPol = 			new TF1*[fNBinsPt];
	fFitSignalInvMassPtBin = 			new TF1*[fNBinsPt];
	fFitSignalInvMassBackFitPtBin = 			new TF1*[fNBinsPt];
	fFitTrueSignalInvMassPtBin = 			new TF1*[fNBinsPt];
	fFitSignalPeakPosInvMassPtBin = 		new TF1*[fNBinsPt];
	fFitSignalPeakPosInvMassBackFitPtBin = 		new TF1*[fNBinsPt];
	fFitBckInvMassPtBin = 				new TF1*[fNBinsPt];
	fFitBckInvMassBackFitPtBin =				new TF1*[fNBinsPt];
	fFitRatioInvMassPtBin = 				new TF1*[fNBinsPt];
	// Histograms for normalization on the left of the peak
	fHistoMappingBackNormInvMassLeftPtBin = new TH1D*[fNBinsPt];
	fHistoMappingSignalInvMassLeftPtBin = 	new TH1D*[fNBinsPt];
	fHistoMappingPeakPosInvMassPtBin = 	new TH1D*[fNBinsPt];

	fFitInvMassLeftPtBin = 				new TF1*[fNBinsPt];
	fFitSignalPeakPosInvMassLeftPtBin = 	new TF1*[fNBinsPt];
	fFitBckInvMassLeftPtBin = 			new TF1*[fNBinsPt];
	fFitWithPol2ForBG = 				new TF1*[fNBinsPt];
	fFitPeakPosPtBin = 					new TF1*[fNBinsPt];
	fMesonMassPeakPos =					new Double_t[fNBinsPt];
	fMesonMassPeakPosError = 			new Double_t[fNBinsPt];

	for(Int_t i = 0;i<fNBinsPt; i++){
		fHistoMappingTrueMesonInvMassPtBins[i] = NULL;
		fHistoMappingTrueGGMesonInvMassPtBins[i] = NULL;
		fHistoMappingTrueDalitzInvMassPtBins[i] = NULL;

		fHistoMappingPiPlPiMiGammaInvMassPtBin[i] = NULL;
		fHistoMappingGGInvMassBackFitPtBin[i] = NULL;
		fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i] = NULL;
		fBackgroundFitPol[i] = NULL;
		fHistoMappingBackInvMassPtBin[i] = NULL;
		fHistoMappingBackNormInvMassPtBin[i] = NULL;
		fHistoMappingSignalInvMassPtBin[i] = NULL;
		fHistoMappingRatioSBInvMassPtBin[i] = NULL;

		fFitSignalInvMassPtBin[i] = NULL;
		fFitSignalInvMassBackFitPtBin[i] = NULL;
		fFitTrueSignalInvMassPtBin[i] = NULL;
		fFitSignalPeakPosInvMassPtBin[i] = NULL;
		fFitSignalPeakPosInvMassBackFitPtBin[i] = NULL;
		fFitBckInvMassPtBin[i] = NULL;
		fFitBckInvMassBackFitPtBin[i] = NULL;
		fFitRatioInvMassPtBin[i] = NULL;
		// Histograms for normalization on the left of the peak
		fHistoMappingBackNormInvMassLeftPtBin[i] = NULL;
		fHistoMappingSignalInvMassLeftPtBin[i] = NULL;
		fHistoMappingPeakPosInvMassPtBin[i] = NULL;
		
		fFitInvMassLeftPtBin[i] = NULL;
		fFitSignalPeakPosInvMassLeftPtBin[i] = NULL;
		fFitBckInvMassLeftPtBin[i] = NULL;
		fFitWithPol2ForBG[i] = NULL;
		fFitPeakPosPtBin[i] = NULL;
	}
}

void CalculateFWHM(TF1 * fFunc)
{
	// Default function
	TF1* fFunc_def;
	fFunc_def = new TF1("fFunc_def","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);
	fFunc_def->SetParameter(0,fFunc->GetParameter(0));
	fFunc_def->SetParameter(1,fFunc->GetParameter(1));
	fFunc_def->SetParameter(2,fFunc->GetParameter(2));
	fFunc_def->SetParameter(3,fFunc->GetParameter(3));



	//FWHM
	fFWHMFunc = fFunc->GetX(fFunc_def->GetParameter(0)*0.5,fFunc_def->GetParameter(1), fMesonFitRange[1]) - fFunc_def->GetX(fFunc_def->GetParameter(0)*0.5,fMesonFitRange[0],fFunc_def->GetParameter(1));

	//FWHM error +
	TF1* fFunc_plus;
	//	fFunc_plus = fFunc;
	fFunc_plus = new TF1("fFunc_plus","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);

	fFunc_plus->SetParameter(0,fFunc->GetParameter(0) + fFunc->GetParError(0));
	fFunc_plus->SetParameter(1,fFunc->GetParameter(1) + fFunc->GetParError(1));
	fFunc_plus->SetParameter(2,fFunc->GetParameter(2) + fFunc->GetParError(2));
	fFunc_plus->SetParameter(3,fFunc->GetParameter(3) + fFunc->GetParError(3));
	Double_t FWHM_plus = fFunc_plus->GetX(fFunc_plus->GetParameter(0)*0.5,fFunc_plus->GetParameter(1), fMesonFitRange[1]) - fFunc_plus->GetX(fFunc_plus->GetParameter(0)*0.5,fMesonFitRange[0],fFunc_plus->GetParameter(1));

	//FWHM error -
	TF1* fFunc_minus;
	//	fFunc_minus = fFunc;
	fFunc_minus = new TF1("fFunc_minus","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);
	fFunc_minus->SetParameter(0,fFunc->GetParameter(0) - fFunc->GetParError(0));
	fFunc_minus->SetParameter(1,fFunc->GetParameter(1) - fFunc->GetParError(1));
	fFunc_minus->SetParameter(2,fFunc->GetParameter(2) - fFunc->GetParError(2));
	fFunc_minus->SetParameter(3,fFunc->GetParameter(3) - fFunc->GetParError(3));
	Double_t FWHM_minus =  fFunc_minus->GetX(fFunc_minus->GetParameter(0)*0.5,fFunc_minus->GetParameter(1), fMesonFitRange[1]) -fFunc_minus->GetX(fFunc_minus->GetParameter(0)*0.5,fMesonFitRange[0],fFunc_minus->GetParameter(1));

	Double_t Error1 = TMath::Abs(fFWHMFunc-FWHM_plus);
	Double_t Error2 = TMath::Abs(fFWHMFunc-FWHM_minus);

	if(Error1>=Error2) fFWHMFuncError = Error1;
	if(Error1<Error2) fFWHMFuncError = Error2;
}

void Delete(){
	if (fBinsPt) delete fBinsPt;
	if (fPeakRange) delete fPeakRange;
	if (fFitRange) delete fFitRange;
	if (fBGFitRange) delete fBGFitRange;
	if (fBGFitRangeLeft) delete fBGFitRangeLeft;
	if (fMesonPlotRange) delete fMesonPlotRange;
	if (fMesonIntRange) delete fMesonIntRange;
	if (fMesonIntRangeWide) delete fMesonIntRangeWide;
	if (fMesonIntRangeNarrow) delete fMesonIntRangeNarrow;
	if (fMesonMassRange) delete fMesonMassRange;
	if (fMesonFitRange) delete fMesonFitRange;
	if (fMesonWidthRange) delete fMesonWidthRange;
	if (fMesonLambdaTailRange) delete fMesonLambdaTailRange;
	if (fNRebin) delete fNRebin;
	if (fPiPlPiMiGammaYields) delete fPiPlPiMiGammaYields;
	if (fMesonYieldsBackFit) delete fMesonYieldsBackFit;
	if (fBckYields) delete fBckYields;
	if (fMesonYields) delete fMesonYields;
	if (fMesonTrueYields) delete fMesonTrueYields;
	if (fMesonYieldsFunc) delete fMesonYieldsFunc;
	if (fMesonYieldsResidualBckFunc) delete fMesonYieldsResidualBckFunc;
	if (fMesonYieldsResidualBckFuncBackFit) delete fMesonYieldsResidualBckFuncBackFit;
	if (fMesonYieldsCorResidualBckFunc) delete fMesonYieldsCorResidualBckFunc;
	if (fMesonYieldsCorResidualBckFuncBackFit) delete fMesonYieldsCorResidualBckFuncBackFit;
	if (fMesonYieldsPerEvent) delete fMesonYieldsPerEvent;
	if (fMesonYieldsPerEventBackFit) delete fMesonYieldsPerEventBackFit;
	if (fMesonMass) delete fMesonMass;
	if (fMesonWidth) delete fMesonWidth;
	if (fMesonSB) delete fMesonSB;
	if (fMesonSign) delete fMesonSign;
	if (fMesonTrueSB) delete fMesonTrueSB;
	if (fMesonTrueSign) delete fMesonTrueSign;
	if (fMesonFWHM) delete fMesonFWHM;
	if (fPiPlPiMiGammaYieldsLeft) delete fPiPlPiMiGammaYieldsLeft;
	if (fBckYieldsLeft) delete fBckYieldsLeft;
	if (fMesonYieldsLeft) delete fMesonYieldsLeft;
	if (fMesonYieldsFuncLeft) delete fMesonYieldsFuncLeft;
	if (fMesonYieldsResidualBckFuncLeft) delete fMesonYieldsResidualBckFuncLeft;
	if (fMesonYieldsCorResidualBckFuncLeft) delete fMesonYieldsCorResidualBckFuncLeft;
	if (fMesonYieldsLeftPerEvent) delete fMesonYieldsLeftPerEvent;
	if (fMesonMassLeft) delete fMesonMassLeft;
	if (fMesonWidthLeft) delete fMesonWidthLeft;
	if (fMesonSBLeft) delete fMesonSBLeft;
	if (fMesonSignLeft) delete fMesonSignLeft;
	if (fMesonFWHMLeft) delete fMesonFWHMLeft;
	if (fPiPlPiMiGammaYieldsNarrow) delete fPiPlPiMiGammaYieldsNarrow;
	if (fBckYieldsNarrow) delete fBckYieldsNarrow;
	if (fMesonYieldsNarrow) delete fMesonYieldsNarrow;
	if (fMesonTrueYieldsNarrow) delete fMesonTrueYieldsNarrow;
	if (fMesonYieldsFuncNarrow) delete fMesonYieldsFuncNarrow;
	if (fMesonYieldsResidualBckFuncNarrow) delete fMesonYieldsResidualBckFuncNarrow;
	if (fMesonYieldsCorResidualBckFuncNarrow) delete fMesonYieldsCorResidualBckFuncNarrow;
	if (fMesonYieldsPerEventNarrow) delete fMesonYieldsPerEventNarrow;
	if (fMesonSBNarrow) delete fMesonSBNarrow;
	if (fMesonSignNarrow) delete fMesonSignNarrow;
	if (fPiPlPiMiGammaYieldsLeftNarrow) delete fPiPlPiMiGammaYieldsLeftNarrow;
	if (fBckYieldsLeftNarrow) delete fBckYieldsLeftNarrow;
	if (fMesonYieldsLeftNarrow) delete fMesonYieldsLeftNarrow;
	if (fMesonYieldsFuncLeftNarrow) delete fMesonYieldsFuncLeftNarrow;
	if (fMesonYieldsResidualBckFuncLeftNarrow) delete fMesonYieldsResidualBckFuncLeftNarrow;
	if (fMesonYieldsCorResidualBckFuncLeftNarrow) delete fMesonYieldsCorResidualBckFuncLeftNarrow;
	if (fMesonYieldsLeftPerEventNarrow) delete fMesonYieldsLeftPerEventNarrow;
	if (fMesonSBLeftNarrow) delete fMesonSBLeftNarrow;
	if (fMesonSignLeftNarrow) delete fMesonSignLeftNarrow;
	if (fPiPlPiMiGammaYieldsWide) delete fPiPlPiMiGammaYieldsWide;
	if (fBckYieldsWide) delete fBckYieldsWide;
	if (fMesonYieldsWide) delete fMesonYieldsWide;
	if (fMesonTrueYieldsWide) delete fMesonTrueYieldsWide;
	if (fMesonYieldsFuncWide) delete fMesonYieldsFuncWide;
	if (fMesonYieldsResidualBckFuncWide) delete fMesonYieldsResidualBckFuncWide;
	if (fMesonYieldsCorResidualBckFuncWide) delete fMesonYieldsCorResidualBckFuncWide;
	if (fMesonYieldsPerEventWide) delete fMesonYieldsPerEventWide;
	if (fMesonSBWide) delete fMesonSBWide;
	if (fMesonSignWide) delete fMesonSignWide;
	if (fPiPlPiMiGammaYieldsLeftWide) delete fPiPlPiMiGammaYieldsLeftWide;
	if (fBckYieldsLeftWide) delete fBckYieldsLeftWide;
	if (fMesonYieldsLeftWide) delete fMesonYieldsLeftWide;
	if (fMesonYieldsFuncLeftWide) delete fMesonYieldsFuncLeftWide;
	if (fMesonYieldsResidualBckFuncLeftWide) delete fMesonYieldsResidualBckFuncLeftWide;
	if (fMesonYieldsCorResidualBckFuncLeftWide) delete fMesonYieldsCorResidualBckFuncLeftWide;
	if (fMesonYieldsLeftPerEventWide) delete fMesonYieldsLeftPerEventWide;
	if (fMesonSBLeftWide) delete fMesonSBLeftWide;
	if (fMesonSignLeftWide) delete fMesonSignLeftWide;
	if (fMesonYieldsBackFitError) delete fMesonYieldsBackFitError;
	if (fBckYieldsError) delete fBckYieldsError;
	if (fMesonYieldsError) delete fMesonYieldsError;
	if (fMesonYieldsFuncError) delete fMesonYieldsFuncError;
	if (fMesonYieldsResidualBckFuncError) delete fMesonYieldsResidualBckFuncError;
	if (fMesonYieldsResidualBckFuncBackFitError) delete fMesonYieldsResidualBckFuncBackFitError;
	if (fMesonYieldsCorResidualBckFuncError) delete fMesonYieldsCorResidualBckFuncError;
	if (fMesonYieldsCorResidualBckFuncBackFitError) delete fMesonYieldsCorResidualBckFuncBackFitError;
	if (fMesonYieldsPerEventError) delete fMesonYieldsPerEventError;
	if (fMesonYieldsPerEventBackFitError) delete fMesonYieldsPerEventBackFitError;
	if (fMesonMassError) delete fMesonMassError;
	if (fMesonWidthError) delete fMesonWidthError;
	if (fMesonSBError) delete fMesonSBError;
	if (fMesonSignError) delete fMesonSignError;
	if (fMesonTrueSBError) delete fMesonTrueSBError;
	if (fMesonTrueSignError) delete fMesonTrueSignError;
	if (fMesonFWHMError) delete fMesonFWHMError;
	if (fPiPlPiMiGammaYieldsLeftError) delete fPiPlPiMiGammaYieldsLeftError;
	if (fBckYieldsLeftError) delete fBckYieldsLeftError;
	if (fMesonYieldsLeftError) delete fMesonYieldsLeftError;
	if (fMesonYieldsFuncLeftError) delete fMesonYieldsFuncLeftError;
	if (fMesonYieldsResidualBckFuncLeftError) delete fMesonYieldsResidualBckFuncLeftError;
	if (fMesonYieldsCorResidualBckFuncLeftError) delete fMesonYieldsCorResidualBckFuncLeftError;
	if (fMesonYieldsLeftPerEventError) delete fMesonYieldsLeftPerEventError;
	if (fMesonMassLeftError) delete fMesonMassLeftError;
	if (fMesonWidthLeftError) delete fMesonWidthLeftError;
	if (fMesonSBLeftError) delete fMesonSBLeftError;
	if (fMesonSignLeftError) delete fMesonSignLeftError;
	if (fMesonFWHMLeftError) delete fMesonFWHMLeftError;
	if (fPiPlPiMiGammaYieldsNarrowError) delete fPiPlPiMiGammaYieldsNarrowError;
	if (fBckYieldsNarrowError) delete fBckYieldsNarrowError;
	if (fMesonYieldsNarrowError) delete fMesonYieldsNarrowError;
	if (fMesonYieldsFuncNarrowError) delete fMesonYieldsFuncNarrowError;
	if (fMesonYieldsResidualBckFuncNarrowError) delete fMesonYieldsResidualBckFuncNarrowError;
	if (fMesonYieldsCorResidualBckFuncNarrowError) delete fMesonYieldsCorResidualBckFuncNarrowError;
	if (fMesonYieldsPerEventNarrowError) delete fMesonYieldsPerEventNarrowError;
	if (fPiPlPiMiGammaYieldsLeftNarrowError) delete fPiPlPiMiGammaYieldsLeftNarrowError;
	if (fBckYieldsLeftNarrowError) delete fBckYieldsLeftNarrowError;
	if (fMesonYieldsLeftNarrowError) delete fMesonYieldsLeftNarrowError;
	if (fMesonYieldsFuncLeftNarrowError) delete fMesonYieldsFuncLeftNarrowError;
	if (fMesonYieldsResidualBckFuncLeftNarrowError) delete fMesonYieldsResidualBckFuncLeftNarrowError;
	if (fMesonYieldsCorResidualBckFuncLeftNarrowError) delete fMesonYieldsCorResidualBckFuncLeftNarrowError;
	if (fMesonYieldsLeftPerEventNarrowError) delete fMesonYieldsLeftPerEventNarrowError;
	if (fMesonSBLeftNarrowError) delete fMesonSBLeftNarrowError;
	if (fMesonSignLeftNarrowError) delete fMesonSignLeftNarrowError;
	if (fPiPlPiMiGammaYieldsWideError) delete fPiPlPiMiGammaYieldsWideError;
	if (fBckYieldsWideError) delete fBckYieldsWideError;
	if (fMesonYieldsWideError) delete fMesonYieldsWideError;
	if (fMesonYieldsFuncWideError) delete fMesonYieldsFuncWideError;
	if (fMesonYieldsResidualBckFuncWideError) delete fMesonYieldsResidualBckFuncWideError;
	if (fMesonYieldsCorResidualBckFuncWideError) delete fMesonYieldsCorResidualBckFuncWideError;
	if (fMesonYieldsPerEventWideError) delete fMesonYieldsPerEventWideError;
	if (fMesonSBWideError) delete fMesonSBWideError;
	if (fMesonSignWideError) delete fMesonSignWideError;
	if (fPiPlPiMiGammaYieldsLeftWideError) delete fPiPlPiMiGammaYieldsLeftWideError;
	if (fBckYieldsLeftWideError) delete fBckYieldsLeftWideError;
	if (fMesonYieldsLeftWideError) delete fMesonYieldsLeftWideError;
	if (fMesonYieldsFuncLeftWideError) delete fMesonYieldsFuncLeftWideError;
	if (fMesonYieldsResidualBckFuncLeftWideError) delete fMesonYieldsResidualBckFuncLeftWideError;
	if (fMesonYieldsCorResidualBckFuncLeftWideError) delete fMesonYieldsCorResidualBckFuncLeftWideError;
	if (fMesonYieldsLeftPerEventWideError) delete fMesonYieldsLeftPerEventWideError;
	if (fMesonSBLeftWideError) delete fMesonSBLeftWideError;
	if (fMesonSignLeftWideError) delete fMesonSignLeftWideError;
	if (fHistoMappingTrueMesonInvMassPtBins) delete fHistoMappingTrueMesonInvMassPtBins;
	if (fHistoMappingPiPlPiMiGammaInvMassPtBin) delete fHistoMappingPiPlPiMiGammaInvMassPtBin;
	if (fHistoMappingBackInvMassPtBin) delete fHistoMappingBackInvMassPtBin;
	if (fHistoMappingBackNormInvMassPtBin) delete fHistoMappingBackNormInvMassPtBin;
	if (fHistoMappingSignalInvMassPtBin) delete fHistoMappingSignalInvMassPtBin;
	if (fHistoMappingRatioSBInvMassPtBin) delete fHistoMappingRatioSBInvMassPtBin;
	if (fFitSignalInvMassPtBin) delete fFitSignalInvMassPtBin;
	if (fFitSignalPeakPosInvMassPtBin) delete fFitSignalPeakPosInvMassPtBin;
	if (fFitBckInvMassPtBin) delete fFitBckInvMassPtBin;
	if (fHistoMappingBackNormInvMassLeftPtBin) delete fHistoMappingBackNormInvMassLeftPtBin;
	if (fHistoMappingSignalInvMassLeftPtBin) delete fHistoMappingSignalInvMassLeftPtBin;
	if (fFitInvMassLeftPtBin) delete fFitInvMassLeftPtBin;
	if (fFitSignalPeakPosInvMassLeftPtBin) delete fFitSignalPeakPosInvMassLeftPtBin;
	if (fFitBckInvMassLeftPtBin) delete fFitBckInvMassLeftPtBin;
	if (fFitRatioInvMassPtBin) delete fFitRatioInvMassPtBin;
	if (fFitWithPol2ForBG) delete fFitWithPol2ForBG;
	if (fMesonFWHMAlpha01) delete fMesonFWHMAlpha01;
	if (fMesonFWHMAlpha01Error) delete fMesonFWHMAlpha01Error;
	if (fHistoWeightsBGZbinVsMbin) delete fHistoWeightsBGZbinVsMbin;
	if (fHistoFillPerEventBGZbinVsMbin) delete fHistoFillPerEventBGZbinVsMbin;
}


void SetCorrectMCHistogrammNames(TString optionEnergyPeriods, TString optionEnergyEnergy){
	cout << "standard MC chosen" << endl;
	ObjectNameTrue = "ESD_TrueMotherPiPlPiMiGamma_InvMass_Pt";
	ObjectNameContaminationGG = "ESD_TrueMotherGG_InvMass_Pt";;
	ObjectNameContaminationDalitz = "ESD_TrueMotherDalitz_InvMass_Pt";
	ObjectNameMCEtaAcc = "MC_EtaInAcc_Pt";
	ObjectNameMCEta = "MC_Eta_Pt";
	ObjectNameMCEtaWOWeights = "MC_Eta_Pt";
	if (optionEnergyPeriods || optionEnergyEnergy){}
}

