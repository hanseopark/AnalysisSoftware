//  **********************************************************************************
//  ******     provided by Gamma Conversion Group, PWGGA,                        *****
//  ******     Friederike Bock, friederike.bock@cern.ch                          *****
//  **********************************************************************************

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
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TMarker.h"
#include "TGraphAsymmErrors.h" 
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TPaveText.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
// #include "CalculateGammaToPi0.h"
//#include "ExtractSignal.h" 

struct SysErrorConversion {
	Double_t value;
	Double_t error;
	// TString name;
};

Double_t Pi0DalitzBR = 0.0;
Double_t Pi0GGBR     = 0.0;

void CorrectYieldDalitz(TH1D* histoCorrectedYield,TH1D* histoRawGGYield, TH1D* histoEffiPt, TH1D* histoAcceptance, Double_t deltaRapid, Double_t scaling, Double_t nEvt, TString nameMeson ){
	histoCorrectedYield->Sumw2();
	histoCorrectedYield->Add(histoRawGGYield,-1.);
	histoCorrectedYield->Scale(1./nEvt);
	histoCorrectedYield->Divide(histoCorrectedYield,histoEffiPt,1.,1.,"");
	histoCorrectedYield->Divide(histoCorrectedYield,histoAcceptance,1.,1.,"");
	histoCorrectedYield->Scale(1./deltaRapid);
	histoCorrectedYield->Scale(scaling);
	for (Int_t i = 1; i < histoCorrectedYield->GetNbinsX()+1 ; i++){
		Double_t newBinContent = histoCorrectedYield->GetBinContent(i)/histoCorrectedYield->GetBinCenter(i);
		Double_t newBinError = histoCorrectedYield->GetBinError(i)/histoCorrectedYield->GetBinCenter(i);
		histoCorrectedYield->SetBinContent(i,newBinContent);
		histoCorrectedYield->SetBinError(i,newBinError);
	}
	if ( nameMeson.CompareTo("Pi0") == 0 || nameMeson.CompareTo("Pi0EtaBinning") == 0 ){
	
		cout<<"Scaling by Pi0DalitzBR: "<<Pi0DalitzBR<<endl;	
		histoCorrectedYield->Scale(1./Pi0DalitzBR);	 
		
	}else{
	  
		histoCorrectedYield->Scale(1./0.000068);
	}
}


void CorrectYield(TH1D* histoCorrectedYield, TH1D* histoRawSecYield, TH1D* histoRawSecYieldFromK0s, TH1D* histoEffiPt, TH1D* histoAcceptance, Double_t deltaRapid, Double_t scaling, Double_t nEvt, TString nameMeson){
	histoCorrectedYield->Sumw2();
	histoCorrectedYield->Add(histoRawSecYield,-1.);
	histoCorrectedYield->Add(histoRawSecYieldFromK0s,-1.);
	histoCorrectedYield->Scale(1./nEvt);
	histoCorrectedYield->Divide(histoCorrectedYield,histoEffiPt,1.,1.,"");
	histoCorrectedYield->Divide(histoCorrectedYield,histoAcceptance,1.,1.,"");
	histoCorrectedYield->Scale(1./deltaRapid);
	histoCorrectedYield->Scale(scaling);
	for (Int_t i = 1; i < histoCorrectedYield->GetNbinsX()+1 ; i++){
		Double_t newBinContent = histoCorrectedYield->GetBinContent(i)/histoCorrectedYield->GetBinCenter(i);
		Double_t newBinError = histoCorrectedYield->GetBinError(i)/histoCorrectedYield->GetBinCenter(i);
		histoCorrectedYield->SetBinContent(i,newBinContent);
		histoCorrectedYield->SetBinError(i,newBinError);
	}
	if (nameMeson.CompareTo("Pi0") == 0 ||nameMeson.CompareTo("Pi0EtaBinning") == 0 ){
	
		cout<<"Scaling by Pi0GGBR: "<<Pi0GGBR<<endl;
		histoCorrectedYield->Scale( 1./Pi0GGBR );
		
		
	}else{
		histoCorrectedYield->Scale(1./0.3931);
	}
	
}

void CompileFullCorrectionFactor(TH1D* histoEffiPt, TH1D* histoAcceptance, Double_t deltaRapid){
	histoEffiPt->Sumw2();
	histoEffiPt->Multiply(histoEffiPt,histoAcceptance,1.,1.,"");
	histoEffiPt->Scale(deltaRapid);
}


void ScaleMCYield(TH1D* histoCorrectedToBeScaled, Double_t deltaRapid, Double_t scaling, Double_t nEvtMC, TString nameMeson, Bool_t optionDalitz ){
	histoCorrectedToBeScaled->Sumw2();
	histoCorrectedToBeScaled->Scale(1./deltaRapid);
	histoCorrectedToBeScaled->Scale(scaling);
	histoCorrectedToBeScaled->Scale(1./nEvtMC);
	for (Int_t i = 1; i < histoCorrectedToBeScaled->GetNbinsX()+1 ; i++){
		Double_t newBinContent = histoCorrectedToBeScaled->GetBinContent(i)/histoCorrectedToBeScaled->GetBinCenter(i);
		Double_t newBinError = histoCorrectedToBeScaled->GetBinError(i)/histoCorrectedToBeScaled->GetBinCenter(i);
		histoCorrectedToBeScaled->SetBinContent(i,newBinContent);
		histoCorrectedToBeScaled->SetBinError(i,newBinError);
	}
	if (nameMeson.CompareTo("Pi0") == 0 ||nameMeson.CompareTo("Pi0EtaBinning") == 0 ){
	  
		if (!optionDalitz){
		  
		  	     cout<<"Scaling MCYield by Pi0GGBR: "<<Pi0GGBR<<endl;
			     histoCorrectedToBeScaled->Scale(1./Pi0GGBR); 
			
		} else {
			
			     cout<<"Scaling MCYield by Pi0DalitzBR: "<<Pi0DalitzBR<<endl;
			     histoCorrectedToBeScaled->Scale( 1./Pi0DalitzBR);
		}
	} else{
		if (!optionDalitz){
			histoCorrectedToBeScaled->Scale(1./0.3931);
		} else {
			histoCorrectedToBeScaled->Scale(1./6.8e-5);
		}
		
	}
}


void  CorrectSignalDalitzV2(TString fileNameUnCorrectedFile = "myOutput", TString fileNameCorrectionFile = "", TString fCutSelection ="", TString suffix = "gif", TString nameMeson ="", TString isMC ="", TString optionEnergy="", TString optionPeriod="", TString fEstimatePileup="", Bool_t optDalitz = kFALSE, TString MCSample="", Int_t mode = 9){  
	
	gROOT->Reset();   
	gROOT->SetStyle("Plain");
	
	StyleSettingsThesis();  
	
	SetPlotStyle();
	
	TString fDalitz="";
	Bool_t kDalitz = optDalitz;
	if (kDalitz){
		fDalitz = "Dalitz";
	} else {
		fDalitz = "";
	}
	
	
	
	
	Bool_t scaleGGCont = kFALSE;
	
	
	if ( isMC.CompareTo("kFALSE") == 0 ) {
	  
	   scaleGGCont 	= kTRUE;
	   cout<<"Scaled gg is set as kTRUE "<<endl;
	  
	  
	} else if ( isMC.CompareTo("kTRUE") == 0 ){
	  
	    scaleGGCont = kFALSE;
	    cout<<"Scaled gg is set as kFALSE because of MC "<<endl;
	  
	}

	
	TString outputDir = Form("%s/%s/%s/%s/CorrectSignalV2%s",fCutSelection.Data(),optionEnergy.Data(),optionPeriod.Data(),suffix.Data(), fDalitz.Data());
	gSystem->Exec("mkdir "+outputDir);
		
	
	TString date = ReturnDateString();
   
	//declaration for printing logo 
	//    Float_t  pictDrawingCoordinates[9] =      {0.55, 0.8, 0.25, 0.04, 0.7,0.5, 0.18, 0.035,0};
	//    Float_t  pictDrawingCoordinatesEff[9] =   {0.55, 0.8, 0.18, 0.04, 0.3, 0.2, 0.18, 0.035,0};
	//    Float_t  pictDrawingCoordinatesAcc[9] =   {0.63, 0.2, 0.40, 0.04, 0.7,0.33, 0.18, 0.035,0};
	Float_t  pictDrawingCoordinatesInv[9] =   {0.63, 0.8, 0.40, 0.04, 0.7, 0.43, 0.18, 0.035,0}; // kk
	//    Float_t  pictDrawingCoordinatesMass[9] =  {0.4, 0.8, 0.75, 0.04, 0.2,0.7, 0.18, 0.035,0};
	//    Float_t  pictDrawingCoordinatesFWHM[9] =  {0.4, 0.8, 0.75, 0.04, 0.2,0.68, 0.18, 0.035,0};
	//    Bool_t   pictDrawingOptions[4] =          {kFALSE, kFALSE, kFALSE, kTRUE};
	TString  prefix2="";
   
	
	
	
	TString fEventCutSelection="";
	TString fGammaCutSelection="";
	TString fClusterCutSelection="";
	TString fElectronCutSelection="";
	TString fMesonCutSelection="";
	
	if (mode == 9  ){
	  
		if (kDalitz){
			ReturnSeparatedCutNumber(fCutSelection, fGammaCutSelection, fElectronCutSelection,fMesonCutSelection,kTRUE);
		} else {
			ReturnSeparatedCutNumber(fCutSelection, fGammaCutSelection, fElectronCutSelection,fMesonCutSelection);
		}
		fEventCutSelection = fGammaCutSelection(0,7);
		fGammaCutSelection = fGammaCutSelection(7,fGammaCutSelection.Length()-7);
		cout << fEventCutSelection.Data() << "\t" << fGammaCutSelection.Data() << endl;
	} else {
		ReturnSeparatedCutNumberAdvanced(fCutSelection,fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);
	}


	TString collisionSystem= ReturnFullCollisionsSystem(optionEnergy);
	Double_t energy = ReturnCollisionEnergy( optionEnergy);
	Double_t doubleAddFactorK0s = ReturnCorrectK0ScalingFactor( optionEnergy,  fEventCutSelection);
	if (isMC.CompareTo("kTRUE") ==0){ 
		prefix2 =         "MC";
		doubleAddFactorK0s = 0.;
	//       pictDrawingOptions[1] =    kTRUE;
	} else { 
		prefix2 =         "data";
	//       pictDrawingOptions[1] =    kFALSE;
	}
	
	cout << "The additional K0 correction factor is: "  << doubleAddFactorK0s<<endl;
	TString textMeson=ReturnMesonString ( nameMeson);
	if (textMeson.CompareTo("") == 0) return;
	//    pictDrawingOptions[3] = ReturnMesonOption (nameMeson);

	
	if (collisionSystem.CompareTo("") == 0){
		cout << "No correct collision system specification, has been given" << endl;
		return;     
	}

	
	TString intermediate = GetCentralityString(fGammaCutSelection);
	TString fTextCent ="";
	if (intermediate.CompareTo("pp")==0){
		fTextCent = "MinBias";  
		intermediate = "";
	} else {
		fTextCent = Form("%s central", intermediate.Data());
	}
	if (optionPeriod.CompareTo("") != 0){
		collisionSystem = Form("%s, %s",collisionSystem.Data(),optionPeriod.Data());
	}   
		
	TString rapidityRange = "";
	Double_t deltaRapid =  ReturnRapidityStringAndDouble(fMesonCutSelection, rapidityRange);
		
	//Variable defintion
	Double_t scaling = 1./(2.*TMath::Pi());
	
	// File definitions
	TFile fileUncorrected(fileNameUnCorrectedFile.Data());  
	if (fileUncorrected.IsZombie()) return;
	TH1F *histoNumberOfGoodESDTracksVtx =  (TH1F*)fileUncorrected.Get("GoodESDTracks");
	TH1D *histoEventQuality =        (TH1D*)fileUncorrected.Get("NEvents");
	TH1D *histoUnCorrectedYield =       (TH1D*)fileUncorrected.Get("histoYieldMeson");
	TH1D *histoUnCorrectedYieldWide =      (TH1D*)fileUncorrected.Get("histoYieldMesonWide");
	TH1D *histoUnCorrectedYieldNarrow = (TH1D*)fileUncorrected.Get("histoYieldMesonNarrow");
	TH1D *histoUnCorrectedYieldLeft =      (TH1D*)fileUncorrected.Get("histoYieldMesonLeft");
	TH1D *histoUnCorrectedYieldLeftWide =  (TH1D*)fileUncorrected.Get("histoYieldMesonLeftWide");
	TH1D *histoUnCorrectedYieldLeftNarrow = (TH1D*)fileUncorrected.Get("histoYieldMesonLeftNarrow");
	TH1D *histoFWHMMeson =           (TH1D*)fileUncorrected.Get("histoFWHMMeson");
	TH1D *histoFWHMMesonLeft =          (TH1D*)fileUncorrected.Get("histoFWHMMesonLeft");
	
	TH1D *histoMassMeson =           (TH1D*)fileUncorrected.Get("histoMassMeson");
	TH1D *histoMassMesonLeft =          (TH1D*)fileUncorrected.Get("histoMassMesonLeft");
	
	TH1D *histoMesonSignalFullPtInvMass=   (TH1D*)fileUncorrected.Get("Mapping_GG_InvMass_FullPt");
	TH1D *histoMesonBckNormFullPtInvMass=  (TH1D*)fileUncorrected.Get("Mapping_BackNorm_InvMass_FullPt");
	
	TH1D *histoSignificanceMeson	=			(TH1D*)fileUncorrected.Get("histoSignMeson");
	TH1D *histoSBMeson		= 			(TH1D*)fileUncorrected.Get("histoSBMeson");

	Float_t nEvt = 0;	
	//if (	optionEnergy.CompareTo("PbPb_2.76TeV") == 0 || optionEnergy.CompareTo("pPb_5.023TeV") == 0 ){
	  if (	optionEnergy.CompareTo("PbPb_2.76TeV") == 0 ){
	  
	  
		nEvt = histoEventQuality->GetBinContent(1);
		
	} else {
	  
		nEvt = GetNEvents(histoEventQuality);
		// BinContent 5 - Zvertex-position, BinContent 4 - no Trigger Bit, BinContent 7 - PileUp 
	}
   
	//    pictDrawingCoordinates[8] =         nEvt;
	TH1D *deltaPt =               (TH1D*)fileUncorrected.Get("deltaPt");
	
	for (Int_t i = 0; i < deltaPt->GetNbinsX() +1; i++){
		deltaPt->SetBinError(i, 0);
	}  
	
	Double_t maxPtMeson= histoUnCorrectedYield->GetXaxis()->GetBinUpEdge(histoUnCorrectedYield->GetNbinsX());
	
	TFile* fileCorrectionsSecPi0 = new TFile("ExternalInput/PCM/SecondaryFractionHistogramms7TeV.root");
	if (fileCorrectionsSecPi0->IsZombie()) return;
	TH1D* histoDefaultTrueSecFracMeson =            (TH1D*)fileCorrectionsSecPi0->Get("TrueSecFrac");
	TH1D* histoDefaultTrueSecFracMesonWide =        (TH1D*)fileCorrectionsSecPi0->Get("TrueSecFracWide");
	TH1D* histoDefaultTrueSecFracMesonNarrow =      (TH1D*)fileCorrectionsSecPi0->Get("TrueSecFracNarrow");
	TH1D* histoDefaultTrueSecFracFromK0SMeson =     (TH1D*)fileCorrectionsSecPi0->Get("TrueSecFracFromK0S");
	TH1D* histoDefaultTrueSecFracFromK0SMesonWide = (TH1D*)fileCorrectionsSecPi0->Get("TrueSecFracFromK0SWide");
	TH1D* histoDefaultTrueSecFracFromK0SMesonNarrow = (TH1D*)fileCorrectionsSecPi0->Get("TrueSecFracFromK0SNarrow");
	
	TH1D* histoDefaultTrueSecFracNotFromK0sMeson = (TH1D*)histoDefaultTrueSecFracMeson->Clone("histoDefaultTrueSecFracNotFromK0sMeson");
	histoDefaultTrueSecFracNotFromK0sMeson->Sumw2();
	histoDefaultTrueSecFracNotFromK0sMeson->Add(histoDefaultTrueSecFracFromK0SMeson,-1);
	TH1D* histoDefaultTrueSecFracModeledToData = (TH1D*)histoDefaultTrueSecFracFromK0SMeson->Clone("histoDefaultTrueSecFracModeledToData");
	histoDefaultTrueSecFracModeledToData->Sumw2();
	histoDefaultTrueSecFracModeledToData->Scale(1+doubleAddFactorK0s);
	histoDefaultTrueSecFracModeledToData->Add(histoDefaultTrueSecFracNotFromK0sMeson,1);
	
	TF1* fitDefaultSecFrac = new TF1("fitDefaultSecFrac","[0]/pow(x,[1])");
	fitDefaultSecFrac->SetRange(0.3,16.);
	TFitResultPtr resultSecFrac = histoDefaultTrueSecFracMeson->Fit(fitDefaultSecFrac,"SINRME+","",0.3,16.);
	TF1* fitDefaultSecFracWide = new TF1("fitDefaultSecFrac","[0]/pow(x,[1])");
	fitDefaultSecFracWide->SetRange(0.3,16.);
	TFitResultPtr resultSecFracWide = histoDefaultTrueSecFracMesonWide->Fit(fitDefaultSecFracWide,"SINRME+","",0.3,16.);
	TF1* fitDefaultSecFracNarrow = new TF1("fitDefaultSecFrac","[0]/pow(x,[1])");
	fitDefaultSecFracNarrow->SetRange(0.3,16.);
	TFitResultPtr resultSecFracNarrow = histoDefaultTrueSecFracMesonNarrow->Fit(fitDefaultSecFracNarrow,"SINRME+","",0.3,16.);
	TF1* fitDefaultSecFracFromK0 = new TF1("fitDefaultSecFrac","[0]/pow(x,[1])");
	fitDefaultSecFracFromK0->SetRange(0.3,16.);
	TFitResultPtr resultSecFracFromK0 = histoDefaultTrueSecFracFromK0SMeson->Fit(fitDefaultSecFracFromK0,"SINRME+","",0.3,16.);
	TF1* fitDefaultSecFracFromK0Wide = new TF1("fitDefaultSecFrac","[0]/pow(x,[1])");
	fitDefaultSecFracFromK0Wide->SetRange(0.3,16.);
	TFitResultPtr resultSecFracFromK0Wide = histoDefaultTrueSecFracFromK0SMesonWide->Fit(fitDefaultSecFracFromK0Wide,"SINRME+","",0.3,16.);
	TF1* fitDefaultSecFracFromK0Narrow = new TF1("fitDefaultSecFrac","[0]/pow(x,[1])");
	fitDefaultSecFracFromK0Narrow->SetRange(0.3,16.);
	TFitResultPtr resultSecFracFromK0Narrow = histoDefaultTrueSecFracFromK0SMesonNarrow->Fit(fitDefaultSecFracFromK0Narrow,"SINRME+","",0.3,16.);
	
	TFile* fileCorrections =         new TFile(fileNameCorrectionFile.Data());
	if (fileCorrections->IsZombie()) return;
	TH1F *histoEventQualityMC =         	(TH1F*)fileCorrections->Get("NEvents");
	TH1D *histoEffiPt =  			(TH1D*)fileCorrections->Get("MesonEffiPt"); //not yet correct MesonEffiPt
	TH1D *histoEffiNarrowPt = 		(TH1D*)fileCorrections->Get("MesonNarrowEffiPt");
	TH1D *histoEffiWidePt = 		(TH1D*)fileCorrections->Get("MesonWideEffiPt");
	TH1D *histoEffiLeftPt = 		(TH1D*)fileCorrections->Get("MesonLeftEffiPt");
	TH1D *histoEffiLeftNarrowPt = 		(TH1D*)fileCorrections->Get("MesonLeftNarrowEffiPt");
	TH1D *histoEffiLeftWidePt = 		(TH1D*)fileCorrections->Get("MesonLeftWideEffiPt");
	TH1D *histoAcceptance=              	(TH1D*)fileCorrections->Get("fMCMesonAccepPt");
	
	TH1D *histoTrueEffiPt = NULL;
	TH1D *histoTrueEffiNarrowPt = NULL;
	TH1D *histoTrueEffiWidePt = NULL;
	histoTrueEffiPt =             		(TH1D*)fileCorrections->Get("TrueMesonEffiPt"); 
	histoTrueEffiNarrowPt =       		(TH1D*)fileCorrections->Get("TrueMesonNarrowEffiPt");
	histoTrueEffiWidePt =         		(TH1D*)fileCorrections->Get("TrueMesonWideEffiPt");

	TH1D* histoInputMesonPt =           	(TH1D*)fileCorrections->Get("MC_Meson_genPt");
	TH1D* histoInputMesonOldBinPt =     	(TH1D*)fileCorrections->Get("MC_Meson_genPt_oldBin");
	TH1D* histoInputMesonOldBinPtWOWeights = NULL;
	TH1D* histoInputMesonOldBinPtWeights = 	 NULL;
	histoInputMesonOldBinPtWOWeights =     	(TH1D*)fileCorrections->Get("MC_Meson_genPt_WOWeights");
	histoInputMesonOldBinPtWeights =     	(TH1D*)fileCorrections->Get("MC_Meson_genPt_Weights");
	TH1D* histoMCInputAddedSig = 		NULL;
	TH1D* histoMCInputWOWeightingAddedSig = NULL;
	TH1D* histoMCInputWeightsAddedSig = 	NULL;
	histoMCInputAddedSig =     		(TH1D*)fileCorrections->Get("MC_Meson_genPt_oldBin_AddedSig");
	histoMCInputWOWeightingAddedSig =     	(TH1D*)fileCorrections->Get("MC_Meson_genPt_WOWeights_AddedSig");
	histoMCInputWeightsAddedSig =     	(TH1D*)fileCorrections->Get("MC_Meson_genPt_Weights_AddedSig");
	
	
	TH1D* histoTrueMassMeson =          	(TH1D*)fileCorrections->Get("histoTrueMassMeson");
	TH1D* histoTrueFWHMMeson =          	(TH1D*)fileCorrections->Get("histoTrueFWHMMeson");
	
	TArrayD* fArrayBRPi0Meson = 	    	(TArrayD*)fileCorrections->Get("fArrayBRPi0Meson");
	
	
	if ( nameMeson.CompareTo("Pi0") == 0 || nameMeson.CompareTo("Pi0EtaBinning") == 0 ){
	  
	      if( isMC.CompareTo("kTRUE") == 0 ){
	    
		cout<<"Taking Branching ratios from MC"<<endl;
		Pi0GGBR     =  fArrayBRPi0Meson->GetAt(2);
		Pi0DalitzBR =  fArrayBRPi0Meson->GetAt(3);
		
	      } else {
		
		cout<<"Taking Branching ratios from DPG 2013 "<<endl;
		Pi0GGBR     =  fArrayBRPi0Meson->GetAt(0); //  0.98823 DPG
		Pi0DalitzBR =  fArrayBRPi0Meson->GetAt(1); //  0.01174 DPG
		
	      }
    
		cout<<"The Branching ratios were set as follow:"<<endl;
		cout<<"Pi0->GG    : "<<Pi0GGBR<<endl;
		cout<<"Pi0->e+e-G : "<<Pi0DalitzBR<<endl;
	     
	} 
	
	

	TString centralityCutNumber = fEventCutSelection(0,3);
	TString centralityString = GetCentralityString(fEventCutSelection);

	TH1D *histoTrueEffiPtFixed 			= (TH1D*)histoTrueEffiPt->Clone("histoTrueEffiPtFixed");
	TH1D *histoTrueEffiNarrowPtFixed 		= (TH1D*)histoTrueEffiNarrowPt->Clone("TrueMesonNarrowEffiPtFixed");
	TH1D *histoTrueEffiWidePtFixed 			= (TH1D*)histoTrueEffiWidePt->Clone("TrueMesonWideEffiPt");
	TH1D *histoEffiPtFixed 				= (TH1D*)histoEffiPt->Clone("histoEffiPtFixed");
	TH1D *histoEffiNarrowPtFixed 			= (TH1D*)histoEffiNarrowPt->Clone("MesonNarrowEffiPtFixed");
	TH1D *histoEffiWidePtFixed 			= (TH1D*)histoEffiWidePt->Clone("MesonWideEffiPt");	
	
	//TH1D *histoTrueNormEffiRatio 			= (TH1D*)histoEffiPt->Clone("histoTrueNormEffiRatio");
	
	//histoTrueNormEffiRatio->Divide(histoEffiPt,histoTrueEffiPt,1.,1.,"B");
	
	
	
	
	
	

	histoTrueEffiPtFixed = FixEfficiency(histoTrueEffiPtFixed,histoTrueEffiPt,optionEnergy,centralityString);
	histoTrueEffiNarrowPtFixed = FixEfficiency(histoTrueEffiNarrowPtFixed,histoTrueEffiNarrowPt,optionEnergy,centralityString);
	histoTrueEffiWidePtFixed = FixEfficiency(histoTrueEffiWidePtFixed,histoTrueEffiWidePt,optionEnergy,centralityString);
	histoEffiPtFixed = FixEfficiency(histoEffiPtFixed,histoEffiPt,optionEnergy,centralityString);
	histoEffiNarrowPtFixed = FixEfficiency(histoEffiNarrowPtFixed,histoEffiNarrowPt,optionEnergy,centralityString);
	histoEffiWidePtFixed = FixEfficiency(histoEffiWidePtFixed,histoEffiWidePt,optionEnergy,centralityString);

	Float_t nEvtMC = 0;
	if (optionEnergy.CompareTo("PbPb_2.76TeV") == 0){
		nEvtMC = histoEventQualityMC->GetBinContent(1);
	} else {
		nEvtMC = GetNEvents(histoEventQualityMC);
		// BinContent 5 - Zvertex-position, BinContent 4 - no Trigger Bit, BinContent 7 - PileUp 
	}
			
	TH1D *histoYieldTrueSecFracMeson = NULL;
	TH1D *histoYieldTrueSecFracMesonWide = NULL;
	TH1D *histoYieldTrueSecFracMesonNarrow = NULL;
	TH1D *histoYieldTrueSecFracFromK0SMesonNarrow = NULL;
	TH1D *histoYieldTrueSecFracFromK0SMeson = NULL;
	TH1D *histoYieldTrueSecFracFromK0SMesonWide = NULL;
	TH1D *histoYieldTrueSecFracMeson_orig = NULL;
	TH1D *histoYieldTrueSecFracFromK0SMeson_orig = NULL;
	TF1* fitpPbSecFrac = new TF1("fitpPbSecFrac","[0]/pow(x,[1])");
	TF1* fitpPbSecFracWide = new TF1("fitpPbSecFrac","[0]/pow(x,[1])");
	TF1* fitpPbSecFracNarrow = new TF1("fitpPbSecFrac","[0]/pow(x,[1])");
	TF1* fitpPbSecFracFromK0 = new TF1("fitpPbSecFrac","[0]/pow(x,[1])");
	TF1* fitpPbSecFracFromK0Wide = new TF1("fitpPbSecFrac","[0]/pow(x,[1])");
	TF1* fitpPbSecFracFromK0Narrow = new TF1("fitpPbSecFrac","[0]/pow(x,[1])");
	TH1D* histoYieldSecMesonLeft = NULL;
	TH1D* histoYieldSecMeson = NULL;
	TH1D* histoYieldSecMesonLeftNarrow = NULL;
	TH1D* histoYieldSecMesonNarrow = NULL;
	TH1D* histoYieldSecMesonLeftWide = NULL;
	TH1D* histoYieldSecMesonWide = NULL;
	TH1D* histoYieldSecFromK0SMeson = NULL;
	TH1D* histoYieldSecFromK0SMesonLeft = NULL;
	TH1D* histoYieldSecFromK0SMesonLeftNarrow = NULL;
	TH1D* histoYieldSecFromK0SMesonLeftWide = NULL;
	TH1D* histoYieldSecFromK0SMesonNarrow = NULL;
	TH1D* histoYieldSecFromK0SMesonWide = NULL;
	TH1D *histoYieldTrueGGFracMeson = NULL;
	TH1D *histoYieldTrueGGFracMesonWide = NULL;
	TH1D *histoYieldTrueGGFracMesonNarrow = NULL;
	TH1D *histoYieldTrueGGFracMesonForData = NULL;
	TH1D *histoYieldTrueGGFracMesonWideForData = NULL;
	TH1D *histoYieldTrueGGFracMesonNarrowForData = NULL;
	TH1D* histoYieldGGMesonLeft = NULL;
	TH1D* histoYieldGGMeson = NULL;
	TH1D* histoYieldGGMesonLeftNarrow = NULL;
	TH1D* histoYieldGGMesonNarrow = NULL;
	TH1D* histoYieldGGMesonLeftWide = NULL;
	TH1D* histoYieldGGMesonWide = NULL;
	
	if (!optDalitz){
		histoYieldTrueSecFracMeson =              (TH1D*)fileCorrections->Get("TrueSecFrac");
		histoYieldTrueSecFracMeson_orig =          (TH1D*)histoYieldTrueSecFracMeson->Clone("TrueSecFrac_orig");
		histoYieldTrueSecFracMesonWide =          (TH1D*)fileCorrections->Get("TrueSecFracWide");
		histoYieldTrueSecFracMesonNarrow =        (TH1D*)fileCorrections->Get("TrueSecFracNarrow");
		histoYieldTrueSecFracFromK0SMeson =       (TH1D*)fileCorrections->Get("TrueSecFracFromK0S");
		histoYieldTrueSecFracFromK0SMeson_orig =       (TH1D*)histoYieldTrueSecFracFromK0SMeson->Clone("TrueSecFracFromK0S_orig");
		histoYieldTrueSecFracFromK0SMesonWide =   (TH1D*)fileCorrections->Get("TrueSecFracFromK0SWide");
		histoYieldTrueSecFracFromK0SMesonNarrow = (TH1D*)fileCorrections->Get("TrueSecFracFromK0SNarrow");
		if ((nameMeson.CompareTo("Pi0")==0 || nameMeson.CompareTo("Pi0EtaBinning")==0) && isMC.CompareTo("kTRUE") !=0 && optionEnergy.CompareTo("pPb_5.023TeV") != 0){
			cout << isMC.Data() << endl;
			for (Int_t i = 2; i < histoYieldTrueSecFracMeson->GetNbinsX()+1; i++){
				Double_t ptStart = histoYieldTrueSecFracMeson->GetXaxis()->GetBinLowEdge(i);
				Double_t ptEnd = histoYieldTrueSecFracMeson->GetXaxis()->GetBinUpEdge(i);
				Double_t binWidth = ptEnd-ptStart;
				Double_t secFrac = fitDefaultSecFrac->Integral(ptStart, ptEnd, resultSecFrac->GetParams()) / binWidth;
				Double_t errorSecFrac = fitDefaultSecFrac->IntegralError(ptStart, ptEnd, resultSecFrac->GetParams(), resultSecFrac->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
				histoYieldTrueSecFracMeson->SetBinContent(i, secFrac);
				histoYieldTrueSecFracMeson->SetBinError(i, errorSecFrac);
				
				secFrac = fitDefaultSecFracWide->Integral(ptStart, ptEnd, resultSecFracWide->GetParams()) / binWidth;
				errorSecFrac = fitDefaultSecFracWide->IntegralError(ptStart, ptEnd, resultSecFracWide->GetParams(), resultSecFracWide->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
				histoYieldTrueSecFracMesonWide->SetBinContent(i, secFrac);
				histoYieldTrueSecFracMesonWide->SetBinError(i, errorSecFrac);
				
				secFrac = fitDefaultSecFracNarrow->Integral(ptStart, ptEnd, resultSecFracNarrow->GetParams()) / binWidth;
				errorSecFrac = fitDefaultSecFracNarrow->IntegralError(ptStart, ptEnd, resultSecFracNarrow->GetParams(), resultSecFracNarrow->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
				histoYieldTrueSecFracMesonNarrow->SetBinContent(i, secFrac);
				histoYieldTrueSecFracMesonNarrow->SetBinError(i, errorSecFrac);
				
				secFrac = fitDefaultSecFracFromK0->Integral(ptStart, ptEnd, resultSecFracFromK0->GetParams()) / binWidth;
				errorSecFrac = fitDefaultSecFracFromK0->IntegralError(ptStart, ptEnd, resultSecFracFromK0->GetParams(), resultSecFracFromK0->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
				histoYieldTrueSecFracFromK0SMeson->SetBinContent(i, secFrac);
				histoYieldTrueSecFracFromK0SMeson->SetBinError(i, errorSecFrac);
				
				secFrac = fitDefaultSecFracFromK0Wide->Integral(ptStart, ptEnd, resultSecFracFromK0Wide->GetParams()) / binWidth;
				errorSecFrac = fitDefaultSecFracFromK0Wide->IntegralError(ptStart, ptEnd, resultSecFracFromK0Wide->GetParams(), resultSecFracFromK0Wide->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
				histoYieldTrueSecFracFromK0SMesonWide->SetBinContent(i, secFrac);
				histoYieldTrueSecFracFromK0SMesonWide->SetBinError(i, errorSecFrac);
				
				secFrac = fitDefaultSecFracFromK0Narrow->Integral(ptStart, ptEnd, resultSecFracFromK0Narrow->GetParams()) / binWidth;
				errorSecFrac = fitDefaultSecFracFromK0Narrow->IntegralError(ptStart, ptEnd, resultSecFracFromK0Narrow->GetParams(), resultSecFracFromK0Narrow->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
				histoYieldTrueSecFracFromK0SMesonNarrow->SetBinContent(i, secFrac);
				histoYieldTrueSecFracFromK0SMesonNarrow->SetBinError(i, errorSecFrac);
			}	
		} else if ((nameMeson.CompareTo("Pi0")==0 || nameMeson.CompareTo("Pi0EtaBinning")==0) && isMC.CompareTo("kTRUE") !=0 && optionEnergy.CompareTo("pPb_5.023TeV") == 0){

			cout << "Fit1"<< endl;  
			fitpPbSecFrac->SetRange(0.3,14.);
			TFitResultPtr resultpPbSecFrac = histoYieldTrueSecFracMeson->Fit(fitpPbSecFrac,"SINRME+","",0.3,5.);
			cout << "Fit2"<< endl;  
			fitpPbSecFracWide->SetRange(0.3,5.);
			TFitResultPtr resultpPbSecFracWide = histoYieldTrueSecFracMesonWide->Fit(fitpPbSecFracWide,"SINRME+","",0.3,5.);

			fitpPbSecFracNarrow->SetRange(0.3,5.);
			TFitResultPtr resultpPbSecFracNarrow = histoYieldTrueSecFracMesonNarrow->Fit(fitpPbSecFracNarrow,"SINRME+","",0.3,5.);
			cout << "Fit3"<< endl;
			fitpPbSecFracFromK0->SetRange(0.3,5.);
			TFitResultPtr resultpPbSecFracFromK0 = histoYieldTrueSecFracFromK0SMeson->Fit(fitpPbSecFracFromK0,"SINRME+","",0.3,5.);

			fitpPbSecFracFromK0Wide->SetRange(0.3,5.);
			TFitResultPtr resultpPbSecFracFromK0Wide =histoYieldTrueSecFracFromK0SMesonWide->Fit(fitpPbSecFracFromK0Wide,"SINRME+","",0.3,5.);

			fitpPbSecFracFromK0Narrow->SetRange(0.3,5.);
			TFitResultPtr resultpPbSecFracFromK0Narrow = histoYieldTrueSecFracFromK0SMesonNarrow->Fit(fitpPbSecFracFromK0Narrow,"SINRME+","",0.3,5.);
			
			for (Int_t i = 2; i < histoYieldTrueSecFracMeson->GetNbinsX()+1; i++){
				Double_t ptStart = histoYieldTrueSecFracMeson->GetXaxis()->GetBinLowEdge(i);
				Double_t ptEnd = histoYieldTrueSecFracMeson->GetXaxis()->GetBinUpEdge(i);
				Double_t binWidth = ptEnd-ptStart;
				Double_t secFrac = fitpPbSecFrac->Integral(ptStart, ptEnd, resultpPbSecFrac->GetParams()) / binWidth;
				Double_t errorSecFrac = fitpPbSecFrac->IntegralError(ptStart, ptEnd, resultpPbSecFrac->GetParams(), resultpPbSecFrac->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
				histoYieldTrueSecFracMeson->SetBinContent(i, secFrac);
				histoYieldTrueSecFracMeson->SetBinError(i, errorSecFrac);
					
				secFrac = fitpPbSecFracWide->Integral(ptStart, ptEnd, resultpPbSecFracWide->GetParams()) / binWidth;
				errorSecFrac = fitpPbSecFracWide->IntegralError(ptStart, ptEnd, resultpPbSecFracWide->GetParams(), resultpPbSecFracWide->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
				histoYieldTrueSecFracMesonWide->SetBinContent(i, secFrac);
				histoYieldTrueSecFracMesonWide->SetBinError(i, errorSecFrac);
					
				secFrac = fitpPbSecFracNarrow->Integral(ptStart, ptEnd, resultpPbSecFracNarrow->GetParams()) / binWidth;
				errorSecFrac = fitpPbSecFracNarrow->IntegralError(ptStart, ptEnd, resultpPbSecFracNarrow->GetParams(), resultpPbSecFracNarrow->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
				histoYieldTrueSecFracMesonNarrow->SetBinContent(i, secFrac);
				histoYieldTrueSecFracMesonNarrow->SetBinError(i, errorSecFrac);
					
				secFrac = fitpPbSecFracFromK0->Integral(ptStart, ptEnd, resultpPbSecFracFromK0->GetParams()) / binWidth;
				errorSecFrac = fitpPbSecFracFromK0->IntegralError(ptStart, ptEnd, resultpPbSecFracFromK0->GetParams(), resultpPbSecFracFromK0->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
				histoYieldTrueSecFracFromK0SMeson->SetBinContent(i, secFrac);
				histoYieldTrueSecFracFromK0SMeson->SetBinError(i, errorSecFrac);
					
				secFrac = fitpPbSecFracFromK0Wide->Integral(ptStart, ptEnd, resultpPbSecFracFromK0Wide->GetParams()) / binWidth;
				errorSecFrac = fitpPbSecFracFromK0Wide->IntegralError(ptStart, ptEnd, resultpPbSecFracFromK0Wide->GetParams(), resultpPbSecFracFromK0Wide->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
				histoYieldTrueSecFracFromK0SMesonWide->SetBinContent(i, secFrac);
				histoYieldTrueSecFracFromK0SMesonWide->SetBinError(i, errorSecFrac);
					
				secFrac = fitpPbSecFracFromK0Narrow->Integral(ptStart, ptEnd, resultpPbSecFracFromK0Narrow->GetParams()) / binWidth;
				errorSecFrac = fitpPbSecFracFromK0Narrow->IntegralError(ptStart, ptEnd, resultpPbSecFracFromK0Narrow->GetParams(), resultpPbSecFracFromK0Narrow->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
				histoYieldTrueSecFracFromK0SMesonNarrow->SetBinContent(i, secFrac);
				histoYieldTrueSecFracFromK0SMesonNarrow->SetBinError(i, errorSecFrac);
			}
		}
		
		histoYieldSecMeson = (TH1D*)histoUnCorrectedYield->Clone("SecFracMeson");
		histoYieldSecMeson->Sumw2();
		histoYieldSecMeson->Multiply(histoYieldTrueSecFracMeson);
		histoYieldSecMesonLeft = (TH1D*)histoUnCorrectedYieldLeft->Clone("SecFracMesonLeft");
		histoYieldSecMesonLeft->Sumw2();
		histoYieldSecMesonLeft->Multiply(histoYieldTrueSecFracMeson);
		histoYieldSecMesonNarrow = (TH1D*)histoUnCorrectedYieldNarrow->Clone("SecFracMesonNarrow");
		histoYieldSecMesonNarrow->Sumw2();
		histoYieldSecMesonNarrow->Multiply(histoYieldTrueSecFracMesonNarrow);
		histoYieldSecMesonLeftNarrow = (TH1D*)histoUnCorrectedYieldLeftNarrow->Clone("SecFracMesonLeftNarrow");
		histoYieldSecMesonLeftNarrow->Sumw2();
		histoYieldSecMesonLeftNarrow->Multiply(histoYieldTrueSecFracMesonNarrow);
		histoYieldSecMesonWide = (TH1D*)histoUnCorrectedYieldWide->Clone("SecFracMesonWide");
		histoYieldSecMesonWide->Sumw2();
		histoYieldSecMesonWide->Multiply(histoYieldTrueSecFracMesonWide);
		histoYieldSecMesonLeftWide = (TH1D*)histoUnCorrectedYieldLeftWide->Clone("SecFracMesonLeftWide");
		histoYieldSecMesonLeftWide->Sumw2();
		histoYieldSecMesonLeftWide->Multiply(histoYieldTrueSecFracMesonWide);
		
		histoYieldSecFromK0SMeson = (TH1D*)histoUnCorrectedYield->Clone("SecFracFromK0SMeson");
		histoYieldSecFromK0SMeson->Sumw2();
		histoYieldSecFromK0SMeson->Multiply(histoYieldTrueSecFracFromK0SMeson);
		histoYieldSecFromK0SMeson->Scale(doubleAddFactorK0s);
		histoYieldSecFromK0SMesonLeft = (TH1D*)histoUnCorrectedYieldLeft->Clone("SecFracFromK0SMesonLeft");
		histoYieldSecFromK0SMesonLeft->Sumw2();
		histoYieldSecFromK0SMesonLeft->Multiply(histoYieldTrueSecFracFromK0SMeson);
		histoYieldSecFromK0SMesonLeft->Scale(doubleAddFactorK0s);
		histoYieldSecFromK0SMesonNarrow = (TH1D*)histoUnCorrectedYieldNarrow->Clone("SecFracFromK0SMesonNarrow");
		histoYieldSecFromK0SMesonNarrow->Sumw2();
		histoYieldSecFromK0SMesonNarrow->Multiply(histoYieldTrueSecFracFromK0SMesonNarrow);
		histoYieldSecFromK0SMesonNarrow->Scale(doubleAddFactorK0s);
		histoYieldSecFromK0SMesonLeftNarrow = (TH1D*)histoUnCorrectedYieldLeftNarrow->Clone("SecFracFromK0SMesonLeftNarrow");
		histoYieldSecFromK0SMesonLeftNarrow->Sumw2();
		histoYieldSecFromK0SMesonLeftNarrow->Multiply(histoYieldTrueSecFracFromK0SMesonNarrow);
		histoYieldSecFromK0SMesonLeftNarrow->Scale(doubleAddFactorK0s);
		histoYieldSecFromK0SMesonWide = (TH1D*)histoUnCorrectedYieldWide->Clone("SecFracFromK0SMesonWide");
		histoYieldSecFromK0SMesonWide->Sumw2();
		histoYieldSecFromK0SMesonWide->Multiply(histoYieldTrueSecFracFromK0SMesonWide);
		histoYieldSecFromK0SMesonWide->Scale(doubleAddFactorK0s);
		histoYieldSecFromK0SMesonLeftWide = (TH1D*)histoUnCorrectedYieldLeftWide->Clone("SecFracFromK0SMesonLeftWide");
		histoYieldSecFromK0SMesonLeftWide->Sumw2();
		histoYieldSecFromK0SMesonLeftWide->Multiply(histoYieldTrueSecFracFromK0SMesonWide);
		histoYieldSecFromK0SMesonLeftWide->Scale(doubleAddFactorK0s);
	} else {
	  
	  
		histoYieldTrueGGFracMeson =            (TH1D*)fileCorrections->Get("TrueGGFrac");
		histoYieldTrueGGFracMesonWide =        (TH1D*)fileCorrections->Get("TrueGGFracWide");
		histoYieldTrueGGFracMesonNarrow =      (TH1D*)fileCorrections->Get("TrueGGFracNarrow");
		
			
		
		if( scaleGGCont ) {
		  
		    histoYieldTrueGGFracMesonForData       = (TH1D*)fileCorrections->Get("TrueGGFracForData");
		    histoYieldTrueGGFracMesonWideForData   = (TH1D*)fileCorrections->Get("TrueGGFracWideForData");
		    histoYieldTrueGGFracMesonNarrowForData = (TH1D*)fileCorrections->Get("TrueGGFracNarrowForData");
		     cout<<"The gg contamination fraction will be scaled"<<endl;
		    
		} else {
		  
		  
		    /*
		    histoYieldTrueGGFracMesonForData       = (TH1D*)histoYieldTrueGGFracMeson->Clone("TrueGGFracForData");
		    histoYieldTrueGGFracMesonWideForData   = (TH1D*)histoYieldTrueGGFracMesonWide->Clone("TrueGGFracForData");
		    histoYieldTrueGGFracMesonNarrowForData = (TH1D*)histoYieldTrueGGFracMesonNarrow->Clone("TrueGGFracForData");
		    */
		    
		    histoYieldTrueGGFracMesonForData       = (TH1D*)histoYieldTrueGGFracMeson->Clone("TrueGGFracForData");
		    histoYieldTrueGGFracMesonWideForData   = (TH1D*)histoYieldTrueGGFracMesonWide->Clone("TrueGGFracWideForData");
		    histoYieldTrueGGFracMesonNarrowForData = (TH1D*)histoYieldTrueGGFracMesonNarrow->Clone("TrueGGFracNarrowForData");
		    
		    		    
		}
		
		
		
		
		
		histoYieldGGMeson = (TH1D*)histoUnCorrectedYield->Clone("GGFracMeson");
		histoYieldGGMeson->Sumw2();
		histoYieldGGMeson->Multiply(histoYieldTrueGGFracMesonForData);
		
		
		
		histoYieldGGMesonLeft = (TH1D*)histoUnCorrectedYieldLeft->Clone("GGFracMesonLeft");
		histoYieldGGMesonLeft->Sumw2();
		histoYieldGGMesonLeft->Multiply(histoYieldTrueGGFracMesonForData);
		
		
		histoYieldGGMesonNarrow = (TH1D*)histoUnCorrectedYieldNarrow->Clone("GGFracMesonNarrow");
		histoYieldGGMesonNarrow->Sumw2();
		histoYieldGGMesonNarrow->Multiply(histoYieldTrueGGFracMesonNarrowForData);
		
				
		histoYieldGGMesonLeftNarrow = (TH1D*)histoUnCorrectedYieldLeftNarrow->Clone("GGFracMesonLeftNarrow");
		histoYieldGGMesonLeftNarrow->Sumw2();
		histoYieldGGMesonLeftNarrow->Multiply(histoYieldTrueGGFracMesonNarrowForData);
		
		
		histoYieldGGMesonWide = (TH1D*)histoUnCorrectedYieldWide->Clone("GGFracMesonWide");
		histoYieldGGMesonWide->Sumw2();
		histoYieldGGMesonWide->Multiply(histoYieldTrueGGFracMesonWideForData);
		
		
		histoYieldGGMesonLeftWide = (TH1D*)histoUnCorrectedYieldLeftWide->Clone("GGFracMesonLeftWide");
		histoYieldGGMesonLeftWide->Sumw2();
		histoYieldGGMesonLeftWide->Multiply(histoYieldTrueGGFracMesonWideForData);
		
	}
	
	Double_t mesonMassExpect = 0;
	if(nameMeson.CompareTo("Pi0") == 0  || nameMeson.CompareTo("Pi0EtaBinning") == 0 ) mesonMassExpect = TDatabasePDG::Instance()->GetParticle(111)->Mass();
	if(nameMeson.CompareTo("Eta") == 0 ) mesonMassExpect = TDatabasePDG::Instance()->GetParticle(221)->Mass();
	

		
	TString fileNameDCAData = Form("%s/%s/%s_Data_GammaConvV1DCATestAnalysed%s.root",fCutSelection.Data(),optionEnergy.Data(),nameMeson.Data(),optionPeriod.Data());
	
	TFile* fileDCAAnalysisData =         new TFile(fileNameDCAData.Data());
	Bool_t kDCAFileDataExists = kTRUE;
	if (fileDCAAnalysisData->IsZombie()) kDCAFileDataExists = kFALSE;
	cout << kDCAFileDataExists << endl;
	TH1D *histoFracCatvsPt[6];
	TH1D *histoFracIntHistBGvsPt[6];
	TH1D *histoCorrectionFactorsHistCat[6];
	for (Int_t i = 0; i < 6 ; i++){
		histoFracCatvsPt[i] = NULL;
		histoFracIntHistBGvsPt[i] = NULL;
		histoCorrectionFactorsHistCat[i] = NULL;
	}   
	TH1D* histoCorrectionFactorsHistvsPt = NULL;
	TH1D* histoCorrectionFactorsFitvsPt = NULL;
	TH1D* histoCorrectionFactorsHistvsPtCatA = NULL;
	TH1D* histoCorrectionFactorsHistvsPtCatC = NULL;
	TH1D* histoCorrectionFactorsHistvsPtCatD = NULL;
	TH1D* histoDCAZUnderMesonAllCat_AllPt = NULL;
	if (kDCAFileDataExists){
		histoCorrectionFactorsHistvsPt =             (TH1D*)fileDCAAnalysisData->Get("fHistCorrectionFactorsHistAllCat_vsPt");
		histoCorrectionFactorsFitvsPt =             (TH1D*)fileDCAAnalysisData->Get("fHistCorrectionFactorsFitAllCat_vsPt");
		histoCorrectionFactorsHistvsPtCatA =             (TH1D*)fileDCAAnalysisData->Get("fHistCorrectionFactorsHistvsPt_0");
		histoCorrectionFactorsHistvsPtCatC =             (TH1D*)fileDCAAnalysisData->Get("fHistCorrectionFactorsHistvsPt_1");
		histoCorrectionFactorsHistvsPtCatD =             (TH1D*)fileDCAAnalysisData->Get("fHistCorrectionFactorsHistvsPt_2");
		histoDCAZUnderMesonAllCat_AllPt =  (TH1D*)fileDCAAnalysisData->Get("HistDCAZUnderMesonAllCat_AllPt");
		for (Int_t i = 0; i < 6 ; i++){
			histoFracCatvsPt[i] = (TH1D*)fileDCAAnalysisData->Get(Form("fHistFracCat_%i_vsPt",i+1));
			histoFracIntHistBGvsPt[i] = (TH1D*)fileDCAAnalysisData->Get(Form("fHistFracIntHistBGvsPt_Cat_%i_Variant_1",i+1));
			histoCorrectionFactorsHistCat[i] = (TH1D*)fileDCAAnalysisData->Get(Form("fHistCorrectionFactorsHistCat_%i_Variant_1_vsPt",i+1));
		}
		
	}
	
	
	TString fileNameDCAMonteCarlo = Form("%s/%s/%s_MC_GammaConvV1DCATestAnalysed%s.root",fCutSelection.Data(),optionEnergy.Data(),nameMeson.Data(),optionPeriod.Data());

	TFile* fileDCAAnalysisMonteCarlo =         new TFile(fileNameDCAMonteCarlo.Data());
	Bool_t kDCAFileMCExists = kTRUE;
	if (fileDCAAnalysisMonteCarlo->IsZombie()) kDCAFileMCExists = kFALSE;
	cout << kDCAFileMCExists << endl;
	
	TH1D *histoMCFracCatvsPt[6];
	TH1D *histoMCFracIntHistBGvsPt[6];
	TH1D *histoMCCorrectionFactorsHistCat[6];
	for (Int_t i = 0; i < 6 ; i++){
		histoMCFracCatvsPt[i] = NULL;
		histoMCFracIntHistBGvsPt[i] = NULL;
		histoMCCorrectionFactorsHistCat[i] = NULL;
	}
	//    TH1D *histoMCCorrectionFactorsHistvsPt;
	//    TH1D *histoMCCorrectionFactorsFitvsPt;
	TH1D* histoMCDCAZUnderMesonAllCat_AllPt = NULL;
	TH1D* histoMCDCAZGarbageAllCat_AllPt= NULL;
	TH1D* histoMCDCAZTrueBackgroundAllCat_AllPt= NULL;
	TH1D* histoMCDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt= NULL;
	TH1D* histoMCDCAZTrueSecondaryMesonFromEtaAllCat_AllPt= NULL;
	TH1D* histoMCDCAZTrueSecondaryMesonFromK0sAllCat_AllPt= NULL;
	TH1D* histoMCDCAZTruePrimaryMesonDalitzAllCat_AllPt= NULL;
	TH1D* histoMCDCAZTruePrimaryMesonGammaGammaAllCat_AllPt= NULL;

	if (kDCAFileMCExists){
	//       histoMCCorrectionFactorsHistvsPt =             (TH1D*)fileDCAAnalysisMonteCarlo->Get("fHistCorrectionFactorsHistAllCat_vsPt");
	//       histoMCCorrectionFactorsFitvsPt =             (TH1D*)fileDCAAnalysisMonteCarlo->Get("fHistCorrectionFactorsFitAllCat_vsPt");
		histoMCDCAZUnderMesonAllCat_AllPt =  (TH1D*)fileDCAAnalysisMonteCarlo->Get("HistDCAZUnderMesonAllCat_AllPt");
		histoMCDCAZGarbageAllCat_AllPt =  (TH1D*)fileDCAAnalysisMonteCarlo->Get("HistDCAZGarbageAllCat_AllPt");
		histoMCDCAZTrueBackgroundAllCat_AllPt =  (TH1D*)fileDCAAnalysisMonteCarlo->Get("HistDCAZTrueBackgroundAllCat_AllPt");
		histoMCDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt =  (TH1D*)fileDCAAnalysisMonteCarlo->Get("HistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt");
		histoMCDCAZTrueSecondaryMesonFromEtaAllCat_AllPt =  (TH1D*)fileDCAAnalysisMonteCarlo->Get("HistDCAZTrueSecondaryMesonFromEtaAllCat_AllPt");
		histoMCDCAZTrueSecondaryMesonFromK0sAllCat_AllPt =  (TH1D*)fileDCAAnalysisMonteCarlo->Get("HistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt");
		histoMCDCAZTruePrimaryMesonDalitzAllCat_AllPt =  (TH1D*)fileDCAAnalysisMonteCarlo->Get("HistDCAZTruePrimaryMesonDalitzAllCat_AllPt");
		histoMCDCAZTruePrimaryMesonGammaGammaAllCat_AllPt =  (TH1D*)fileDCAAnalysisMonteCarlo->Get("HistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt");
		for (Int_t i = 0; i < 6 ; i++){
			histoMCFracCatvsPt[i] = (TH1D*)fileDCAAnalysisMonteCarlo->Get(Form("fHistFracCat_%i_vsPt",i+1));
			histoMCFracIntHistBGvsPt[i] = (TH1D*)fileDCAAnalysisMonteCarlo->Get(Form("fHistFracIntHistBGvsPt_Cat_%i_Variant_1",i+1));
			histoMCCorrectionFactorsHistCat[i] = (TH1D*)fileDCAAnalysisMonteCarlo->Get(Form("fHistCorrectionFactorsHistCat_%i_Variant_1_vsPt",i+1));
		}
	}
	
	
	

	Color_t  colorCatAll[3]    = {kBlack,kRed+2,kBlue+1};
	Color_t  colorCatAllMC[3]    = {kGray+2,kRed+3,kBlue+3};
	Color_t  colorCat[6]    = { kRed+1, 807, 800, kGreen+2, kCyan+2, kBlue+1};
	Color_t  colorCatMC[6]    = { kRed+3, 807+2, 800+2, kGreen+4, kCyan+4, kBlue+3};
	Style_t  styleCat[6]    = { 20, 21, 29, 33, 20, 21};
	Style_t  styleCatMC[6]    = { 24, 25, 30, 27, 24, 25};
	
        TH1D *histoBGEstimateA    = NULL;
	TH1D *histoBGEstimateB    = NULL;
	TH1D *histoBGEstimateCatA = NULL;
	TH1D *histoBGEstimateCatC = NULL;
	TH1D *histoBGEstimateCatD = NULL;
	
	
	if( ! kDalitz ) {
	
	histoBGEstimateA = (TH1D*)histoYieldTrueSecFracMeson->Clone("histoBGEstimateA");
	histoBGEstimateB = (TH1D*)histoYieldTrueSecFracMeson->Clone("histoBGEstimateB");
	histoBGEstimateCatA = (TH1D*)histoYieldTrueSecFracMeson->Clone("histoBGEstimateCatA");
	histoBGEstimateCatC = (TH1D*)histoYieldTrueSecFracMeson->Clone("histoBGEstimateCatC");
	histoBGEstimateCatD = (TH1D*)histoYieldTrueSecFracMeson->Clone("histoBGEstimateCatD");
	
	}
	

	
	if (kDCAFileDataExists && kDCAFileMCExists){
	  
		TCanvas* canvasDCAMCComponents = new TCanvas("canvasDCAMCComponents","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasDCAMCComponents, 0.08, 0.02, 0.02, 0.09);
		canvasDCAMCComponents->SetLogy();
			if (histoDCAZUnderMesonAllCat_AllPt){
				DrawAutoGammaMesonHistos( histoDCAZUnderMesonAllCat_AllPt, 
								"","dca_{z} #gamma (cm)", "d(dca_{z})/N_{evt}", 
								kFALSE, 2.,1, kFALSE,
								kTRUE,1e-8,2*histoMCDCAZUnderMesonAllCat_AllPt->GetMaximum(), 
								kTRUE, -6., 6.);
				DrawGammaSetMarker(histoDCAZUnderMesonAllCat_AllPt, 20, 1., kBlack, kBlack);
				histoDCAZUnderMesonAllCat_AllPt->GetYaxis()->SetTitleOffset(0.8);
				histoDCAZUnderMesonAllCat_AllPt->DrawCopy("p,e1"); 

				}
				if (histoMCDCAZUnderMesonAllCat_AllPt){
				DrawGammaSetMarker(histoMCDCAZUnderMesonAllCat_AllPt, 24, 1., kGray, kGray);
				histoMCDCAZUnderMesonAllCat_AllPt->DrawCopy("same,p,e1"); 
				}   
				if (histoMCDCAZTruePrimaryMesonGammaGammaAllCat_AllPt){
				DrawGammaSetMarker(histoMCDCAZTruePrimaryMesonGammaGammaAllCat_AllPt, 20, 1., kRed+2, kRed+2);
				histoMCDCAZTruePrimaryMesonGammaGammaAllCat_AllPt->DrawCopy("same,p,e1"); 

				} 
				if (histoMCDCAZTruePrimaryMesonDalitzAllCat_AllPt){
				DrawGammaSetMarker(histoMCDCAZTruePrimaryMesonDalitzAllCat_AllPt, 20, 1., kGreen+2, kGreen+2);
				histoMCDCAZTruePrimaryMesonDalitzAllCat_AllPt->DrawCopy("same,p,e1"); 
				}
				if (histoMCDCAZTrueBackgroundAllCat_AllPt){
					DrawGammaSetMarker(histoMCDCAZTrueBackgroundAllCat_AllPt, 20, 1., kPink+2, kPink+2);
				histoMCDCAZTrueBackgroundAllCat_AllPt->DrawCopy("same,p,e1"); 
				}
				if (histoMCDCAZTrueSecondaryMesonFromK0sAllCat_AllPt){
					DrawGammaSetMarker(histoMCDCAZTrueSecondaryMesonFromK0sAllCat_AllPt, 20, 1., 807, 807);
				histoMCDCAZTrueSecondaryMesonFromK0sAllCat_AllPt->DrawCopy("same,p,e1"); 
				}
				if (histoMCDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt){
					DrawGammaSetMarker(histoMCDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt, 20, 1., kViolet+2, kViolet+2);
				histoMCDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt->DrawCopy("same,p,e1"); 
				}
				if (histoMCDCAZGarbageAllCat_AllPt){
				DrawGammaSetMarker(histoMCDCAZGarbageAllCat_AllPt, 20, 1., kCyan+2, kCyan+2);
				histoMCDCAZGarbageAllCat_AllPt->DrawCopy("same,p,e1"); 
				}
				
				TLatex *labelEnergy = new TLatex(0.11,0.9,Form("%s %s",intermediate.Data(),collisionSystem.Data()));
				SetStyleTLatex( labelEnergy, 0.04,4);
				labelEnergy->Draw();

			TLegend* legendDCAMCComponents0 = new TLegend(0.7,0.7,0.85,0.95);
			legendDCAMCComponents0->SetTextSize(0.04);         
			legendDCAMCComponents0->SetFillColor(0);
			legendDCAMCComponents0->SetLineColor(0);
			legendDCAMCComponents0->SetNColumns(1);
			if (histoDCAZUnderMesonAllCat_AllPt)legendDCAMCComponents0->AddEntry(histoDCAZUnderMesonAllCat_AllPt,"Data","p");
			if (histoMCDCAZUnderMesonAllCat_AllPt)legendDCAMCComponents0->AddEntry(histoMCDCAZUnderMesonAllCat_AllPt,"Total MC","p");
			if (histoMCDCAZTruePrimaryMesonGammaGammaAllCat_AllPt)legendDCAMCComponents0->AddEntry(histoMCDCAZTruePrimaryMesonGammaGammaAllCat_AllPt,Form("Prim"),"p");
			if (histoMCDCAZTruePrimaryMesonDalitzAllCat_AllPt)legendDCAMCComponents0->AddEntry(histoMCDCAZTruePrimaryMesonDalitzAllCat_AllPt,Form("Dalitz"),"l");
			if (histoMCDCAZTrueSecondaryMesonFromK0sAllCat_AllPt && histoMCDCAZTrueSecondaryMesonFromK0sAllCat_AllPt->GetEntries() > 0)legendDCAMCComponents0->AddEntry(histoMCDCAZTrueSecondaryMesonFromK0sAllCat_AllPt,"Sec. #pi^{0} from K^{0}_{s}","p");
			if (histoMCDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt && histoMCDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt->GetEntries() > 0)legendDCAMCComponents0->AddEntry(histoMCDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt,"Sec. #pi^{0} from X","p");
			if (histoMCDCAZTrueBackgroundAllCat_AllPt)legendDCAMCComponents0->AddEntry(histoMCDCAZTrueBackgroundAllCat_AllPt,"#gamma#gamma BG","p");
			if (histoMCDCAZGarbageAllCat_AllPt)legendDCAMCComponents0->AddEntry(histoMCDCAZGarbageAllCat_AllPt,"garbage","p");
			legendDCAMCComponents0->Draw();         
		canvasDCAMCComponents->Update(); 
		canvasDCAMCComponents->SaveAs(Form("%s/%s_MC_DCAzDecomposition.%s",outputDir.Data(),nameMeson.Data(),suffix.Data()));
			
		
		
	}   
	
	if (kDCAFileDataExists){
		TCanvas* canvasCorrFrac = new TCanvas("canvasCorrFrac","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasCorrFrac, 0.08, 0.02, 0.02, 0.09);

		canvasCorrFrac->cd();
		DrawAutoGammaMesonHistos( histoFracCatvsPt[0], 
								"", "p_{T,#pi^{0}} (GeV/c)", "N_{#pi^{0} per cat}/(N_{#pi^{0}}) (%)", 
								kFALSE, 2.,1e-8, kFALSE,
								kTRUE, 0, 100, 
								kFALSE, 0., 10.);
		DrawGammaSetMarker(histoFracCatvsPt[0], styleCat[0], 1., colorCat[0], colorCat[0]);      
		histoFracCatvsPt[0]->GetYaxis()->SetTitleOffset(0.8);
		histoFracCatvsPt[0]->DrawCopy("p,e1"); 
	
		if(histoMCFracCatvsPt[0]){
			DrawGammaSetMarker(histoMCFracCatvsPt[0], styleCatMC[0], 1., colorCatMC[0], colorCatMC[0]);      
			histoMCFracCatvsPt[0]->GetYaxis()->SetTitleOffset(0.8);
			histoMCFracCatvsPt[0]->DrawCopy("same,p,e1"); 
		}
		TLegend* legendFractionCat = new TLegend(0.65,0.7,0.95,0.95);
		legendFractionCat->SetTextSize(0.04);         
		legendFractionCat->SetFillColor(0);
		legendFractionCat->SetLineColor(0);
		if (kDCAFileMCExists) legendFractionCat->SetNColumns(3);
			else legendFractionCat->SetNColumns(2);
		legendFractionCat->AddEntry((TObject*)0,"Cat 1","");
		legendFractionCat->AddEntry(histoFracCatvsPt[0],"Data","p");
		if(histoMCFracCatvsPt[0])legendFractionCat->AddEntry(histoMCFracCatvsPt[0],"MC","p");
		
		for (Int_t i = 1; i< 6; i++){
			DrawGammaSetMarker(histoFracCatvsPt[i], styleCat[i], 1., colorCat[i], colorCat[i]);
			histoFracCatvsPt[i]->DrawCopy("same,p,e1"); 
			if(histoMCFracCatvsPt[i]){
				DrawGammaSetMarker(histoMCFracCatvsPt[i], styleCatMC[i], 1., colorCatMC[i], colorCatMC[i]);      
				histoMCFracCatvsPt[i]->DrawCopy("same,p,e1"); 
			}
			legendFractionCat->AddEntry((TObject*)0,Form("Cat %i",i+1),"");
			legendFractionCat->AddEntry(histoFracCatvsPt[i],"Data","p");
			if(histoMCFracCatvsPt[i])legendFractionCat->AddEntry(histoMCFracCatvsPt[i],"MC","p");
		}   
		legendFractionCat->Draw();
		canvasCorrFrac->Update(); 
		canvasCorrFrac->SaveAs(Form("%s/%s_FractionPerCatVsPt_ComparedToMC.%s",outputDir.Data(),nameMeson.Data(),suffix.Data()));
		
		
	}   
	
	if (kDCAFileDataExists){
		TCanvas* canvasCorrFrac = new TCanvas("canvasCorrFrac","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasCorrFrac, 0.08, 0.02, 0.02, 0.09);

		canvasCorrFrac->cd();
		DrawAutoGammaMesonHistos( histoFracIntHistBGvsPt[0], 
								"", "p_{T,#pi^{0}} (GeV/c)", "BG/Total (%)", 
								kFALSE, 2.,1e-8, kFALSE,
								kTRUE, 0, 30, 
								kFALSE, 0., 10.);
		DrawGammaSetMarker(histoFracIntHistBGvsPt[0], styleCat[0], 1., colorCat[0], colorCat[0]);      
		histoFracIntHistBGvsPt[0]->GetYaxis()->SetTitleOffset(0.8);
		histoFracIntHistBGvsPt[0]->DrawCopy("p,e1"); 
	
		if(histoMCFracIntHistBGvsPt[0]){
			DrawGammaSetMarker(histoMCFracIntHistBGvsPt[0], styleCatMC[0], 1., colorCatMC[0], colorCatMC[0]);      
			histoMCFracIntHistBGvsPt[0]->DrawCopy("same,p,e1"); 
		}
		TLegend* legendFractionCat = new TLegend(0.65,0.7,0.95,0.95);
		legendFractionCat->SetTextSize(0.04);         
		legendFractionCat->SetFillColor(0);
		legendFractionCat->SetLineColor(0);
		if (kDCAFileMCExists) legendFractionCat->SetNColumns(3);
			else legendFractionCat->SetNColumns(2);
		legendFractionCat->AddEntry((TObject*)0,"Cat 1","");
		legendFractionCat->AddEntry(histoFracIntHistBGvsPt[0],"Data","p");
		if(histoMCFracIntHistBGvsPt[0])legendFractionCat->AddEntry(histoMCFracIntHistBGvsPt[0],"MC","p");
		
		for (Int_t i = 1; i< 6; i++){
			DrawGammaSetMarker(histoFracIntHistBGvsPt[i], styleCat[i], 1., colorCat[i], colorCat[i]);
			histoFracIntHistBGvsPt[i]->DrawCopy("same,p,e1"); 
			if(histoMCFracIntHistBGvsPt[i]){
				DrawGammaSetMarker(histoMCFracIntHistBGvsPt[i], styleCatMC[i], 1., colorCatMC[i], colorCatMC[i]);      
				histoMCFracIntHistBGvsPt[i]->DrawCopy("same,p,e1"); 
			}
			legendFractionCat->AddEntry((TObject*)0,Form("Cat %i",i+1),"");
			legendFractionCat->AddEntry(histoFracIntHistBGvsPt[i],"Data","p");
			if(histoMCFracIntHistBGvsPt[i])legendFractionCat->AddEntry(histoMCFracIntHistBGvsPt[i],"MC","p");
		}   
		legendFractionCat->Draw();
		canvasCorrFrac->Update(); 
		canvasCorrFrac->SaveAs(Form("%s/%s_FracBGOverIntHist_ComparedToMC.%s",outputDir.Data(),nameMeson.Data(),suffix.Data()));
		
		
	}   
	
	Color_t colorMethod[5] = {kBlack, kCyan+2, kRed+2, kGreen+2, kBlue+2};
	Style_t styleMethod[5] = {20,24,21,29,33};
	TString nameMethod[5] = {"A","B","A","C","D"};

	if (kDCAFileDataExists){
	  
		TCanvas* canvasCorrFrac = new TCanvas("canvasCorrFrac","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasCorrFrac, 0.08, 0.02, 0.02, 0.09);

		canvasCorrFrac->cd();
		DrawAutoGammaMesonHistos( histoCorrectionFactorsHistvsPt, 
								"", "p_{T,#pi^{0}} (GeV/c)", "Contamination from Pileup (%)", 
								kFALSE, 2.,1e-8, kFALSE,
								kTRUE, 0, 30, 
								kFALSE, 0., 10.);
		DrawGammaSetMarker(histoCorrectionFactorsHistvsPt, styleMethod[0], 1.2, colorMethod[0], colorMethod[0]);      
		histoCorrectionFactorsHistvsPt->GetYaxis()->SetTitleOffset(0.8);
		histoCorrectionFactorsHistvsPt->DrawCopy("p,e1"); 
		TF1* fitCorrectionFactorsHistvsPt = new TF1("fitCorrectionFactorsHistvsPt","[0]/pow(x,[1])+[2]");
		fitCorrectionFactorsHistvsPt->SetRange(0.4, maxPtMeson);
		TFitResultPtr resultCorrectionFactorsHistvsPt = histoCorrectionFactorsHistvsPt->Fit(fitCorrectionFactorsHistvsPt,"SINRME+","",0.4, maxPtMeson);
		TString bla= WriteParameterToFile(fitCorrectionFactorsHistvsPt);
		cout << bla.Data()<< endl;
		fitCorrectionFactorsHistvsPt->SetLineColor(colorMethod[0]);
		fitCorrectionFactorsHistvsPt->Draw("same");
		
		TF1* fitCorrectionFactorsFitvsPt = NULL;
		TFitResultPtr resultCorrectionFactorsFitvsPt ;
		if (histoCorrectionFactorsFitvsPt) {
			DrawGammaSetMarker(histoCorrectionFactorsFitvsPt, styleMethod[1], 1.2, colorMethod[1], colorMethod[1]);      
			histoCorrectionFactorsFitvsPt->DrawCopy("same,p,e1"); 
			fitCorrectionFactorsFitvsPt = new TF1("fitCorrectionFactorsFitvsPt","[0]/pow(x,[1])+[2]");
			fitCorrectionFactorsFitvsPt->SetRange(0.4, maxPtMeson);
			resultCorrectionFactorsFitvsPt = histoCorrectionFactorsFitvsPt->Fit(fitCorrectionFactorsFitvsPt,"SINRME+","",0.4, maxPtMeson);
			fitCorrectionFactorsFitvsPt->SetLineColor(colorMethod[1]);
			fitCorrectionFactorsFitvsPt->Draw("same");
		}
		DrawGammaSetMarker(histoCorrectionFactorsHistvsPtCatA, styleMethod[2], 1.2, colorMethod[2], colorMethod[2]);      
		histoCorrectionFactorsHistvsPtCatA->DrawCopy("same,p,e1"); 
		TF1* fitCorrectionFactorsHistvsPtCatA = new TF1("fitCorrectionFactorsHistvsPtCatA","[0]/pow(x,[1])+[2]");
		fitCorrectionFactorsHistvsPtCatA->SetRange(0.4, maxPtMeson);
		
		TFitResultPtr resultCorrectionFactorsHistvsPtCatA = histoCorrectionFactorsHistvsPtCatA->Fit(fitCorrectionFactorsHistvsPtCatA,"SINRME+","",0.4, maxPtMeson);
		fitCorrectionFactorsHistvsPtCatA->SetLineColor(colorMethod[2]);
		fitCorrectionFactorsHistvsPtCatA->Draw("same");

		DrawGammaSetMarker(histoCorrectionFactorsHistvsPtCatC, styleMethod[3], 1.2, colorMethod[3], colorMethod[3]);      
		histoCorrectionFactorsHistvsPtCatC->DrawCopy("same,p,e1"); 
		TF1* fitCorrectionFactorsHistvsPtCatC = new TF1("fitCorrectionFactorsHistvsPtCatC","[0]/pow(x,[1])+[2]");
		fitCorrectionFactorsHistvsPtCatC->SetRange(0.4, maxPtMeson);
		TFitResultPtr resultCorrectionFactorsHistvsPtCatC = histoCorrectionFactorsHistvsPtCatC->Fit(fitCorrectionFactorsHistvsPtCatC,"SINRME+","",0.4, maxPtMeson);
		fitCorrectionFactorsHistvsPtCatC->SetLineColor(colorMethod[3]);
		fitCorrectionFactorsHistvsPtCatC->Draw("same");

		DrawGammaSetMarker(histoCorrectionFactorsHistvsPtCatD, styleMethod[4], 1.2, colorMethod[4], colorMethod[4]);      
		histoCorrectionFactorsHistvsPtCatD->DrawCopy("same,p,e1"); 
		TF1* fitCorrectionFactorsHistvsPtCatD = new TF1("fitCorrectionFactorsHistvsPtCatD","[0]/pow(x,[1])+[2]");
		fitCorrectionFactorsHistvsPtCatD->SetRange(0.4, maxPtMeson);
		TFitResultPtr resultCorrectionFactorsHistvsPtCatD = histoCorrectionFactorsHistvsPtCatD->Fit(fitCorrectionFactorsHistvsPtCatD,"SINRME+","",0.4, maxPtMeson);
		fitCorrectionFactorsHistvsPtCatD->SetLineColor(colorMethod[4]);
		fitCorrectionFactorsHistvsPtCatD->Draw("same");
				
		TLegend* legendFractionCat = new TLegend(0.5,0.8,0.95,0.95);
		legendFractionCat->SetTextSize(0.03);         
		legendFractionCat->SetFillColor(0);
		legendFractionCat->SetLineColor(0);
		legendFractionCat->SetNColumns(2);
		legendFractionCat->AddEntry((TObject*)0,"Method A","");
		legendFractionCat->AddEntry(histoCorrectionFactorsHistvsPt,"Data","p");
		if (histoCorrectionFactorsFitvsPt) legendFractionCat->AddEntry((TObject*)0,"Method B","");
		if (histoCorrectionFactorsFitvsPt) legendFractionCat->AddEntry(histoCorrectionFactorsFitvsPt,"Data","p");
		legendFractionCat->AddEntry((TObject*)0,"Method A sep Cat","");
		legendFractionCat->AddEntry(histoCorrectionFactorsHistvsPtCatA,"Data","p");
		legendFractionCat->AddEntry((TObject*)0,"Method C sep Cat","");
		legendFractionCat->AddEntry(histoCorrectionFactorsHistvsPtCatC,"Data","p");
		legendFractionCat->AddEntry((TObject*)0,"Method D sep Cat","");
		legendFractionCat->AddEntry(histoCorrectionFactorsHistvsPtCatD,"Data","p");
	
		TLatex *labelEnergy = new TLatex(0.11,0.9,Form("%s %s",intermediate.Data(),collisionSystem.Data()));
		SetStyleTLatex( labelEnergy, 0.04,4);
		labelEnergy->Draw();
		
		legendFractionCat->Draw();
		canvasCorrFrac->Update(); 
		canvasCorrFrac->SaveAs(Form("%s/%s_FinalBGEstimate_AllMethods.%s",outputDir.Data(),nameMeson.Data(),suffix.Data()));
		
		for (Int_t i = 2; i < histoBGEstimateA->GetNbinsX()+1; i++){
			Double_t ptStart = histoBGEstimateA->GetXaxis()->GetBinLowEdge(i);
			Double_t ptEnd = histoBGEstimateA->GetXaxis()->GetBinUpEdge(i);
			Double_t binWidth = ptEnd-ptStart;
			Double_t bgEstimate = (100-fitCorrectionFactorsHistvsPt->Integral(ptStart, ptEnd, resultCorrectionFactorsHistvsPt->GetParams()) / binWidth )/100.;
			Double_t errorBGEstimate = (fitCorrectionFactorsHistvsPt->IntegralError(ptStart, ptEnd, resultCorrectionFactorsHistvsPt->GetParams(), resultCorrectionFactorsHistvsPt->GetCovarianceMatrix().GetMatrixArray() ) / binWidth )/100.;
			histoBGEstimateA->SetBinContent(i, bgEstimate);
			histoBGEstimateA->SetBinError(i, errorBGEstimate);
			
			if (fitCorrectionFactorsFitvsPt){
				bgEstimate = (100-fitCorrectionFactorsFitvsPt->Integral(ptStart, ptEnd, resultCorrectionFactorsFitvsPt->GetParams()) / binWidth )/100.;
				errorBGEstimate = (fitCorrectionFactorsFitvsPt->IntegralError(ptStart, ptEnd, resultCorrectionFactorsFitvsPt->GetParams(), resultCorrectionFactorsFitvsPt->GetCovarianceMatrix().GetMatrixArray() ) / binWidth )/100.;
				histoBGEstimateB->SetBinContent(i, bgEstimate);
				histoBGEstimateB->SetBinError(i, errorBGEstimate);
			}
			
			bgEstimate = (100-fitCorrectionFactorsHistvsPtCatA->Integral(ptStart, ptEnd, resultCorrectionFactorsHistvsPtCatA->GetParams()) / binWidth )/100.;
			errorBGEstimate = (fitCorrectionFactorsHistvsPtCatA->IntegralError(ptStart, ptEnd, resultCorrectionFactorsHistvsPtCatA->GetParams(), resultCorrectionFactorsHistvsPtCatA->GetCovarianceMatrix().GetMatrixArray() ) / binWidth )/100.;
			histoBGEstimateCatA->SetBinContent(i, bgEstimate);
			histoBGEstimateCatA->SetBinError(i, errorBGEstimate);
			
			bgEstimate = (100-fitCorrectionFactorsHistvsPtCatC->Integral(ptStart, ptEnd, resultCorrectionFactorsHistvsPtCatC->GetParams()) / binWidth )/100.;
			errorBGEstimate = (fitCorrectionFactorsHistvsPtCatC->IntegralError(ptStart, ptEnd, resultCorrectionFactorsHistvsPtCatC->GetParams(), resultCorrectionFactorsHistvsPtCatC->GetCovarianceMatrix().GetMatrixArray() ) / binWidth )/100.;
			histoBGEstimateCatC->SetBinContent(i, bgEstimate);
			histoBGEstimateCatC->SetBinError(i, errorBGEstimate);
			
			bgEstimate = (100-fitCorrectionFactorsHistvsPtCatD->Integral(ptStart, ptEnd, resultCorrectionFactorsHistvsPtCatD->GetParams()) / binWidth )/100.;
			errorBGEstimate = (fitCorrectionFactorsHistvsPtCatD->IntegralError(ptStart, ptEnd, resultCorrectionFactorsHistvsPtCatD->GetParams(), resultCorrectionFactorsHistvsPtCatD->GetCovarianceMatrix().GetMatrixArray() ) / binWidth )/100.;
			histoBGEstimateCatD->SetBinContent(i, bgEstimate);
			histoBGEstimateCatD->SetBinError(i, errorBGEstimate);
		}

		canvasCorrFrac->cd();
		DrawAutoGammaMesonHistos( histoBGEstimateA, 
								"", "p_{T,#pi^{0}} (GeV/c)", "Correction factor", 
								kFALSE, 2.,1e-8, kFALSE,
								kTRUE, 0.6, 1, 
								kFALSE, 0., 10.);
		DrawGammaSetMarker(histoBGEstimateA, styleMethod[0], 1.2, colorMethod[0], colorMethod[0]);      
		histoBGEstimateA->GetYaxis()->SetTitleOffset(0.8);
		histoBGEstimateA->DrawCopy("p,e1"); 
		
		if ( fitCorrectionFactorsFitvsPt) {
			DrawGammaSetMarker(histoBGEstimateB, styleMethod[1], 1.2, colorMethod[1], colorMethod[1]);      
			histoBGEstimateB->DrawCopy("same,p,e1"); 
		}
		
		DrawGammaSetMarker(histoBGEstimateCatA, styleMethod[2], 1.2, colorMethod[2], colorMethod[2]);      
		histoBGEstimateCatA->DrawCopy("same,p,e1"); 
		
		DrawGammaSetMarker(histoBGEstimateCatC, styleMethod[3], 1.2, colorMethod[3], colorMethod[3]);      
		histoBGEstimateCatC->DrawCopy("same,p,e1"); 
		
		DrawGammaSetMarker(histoBGEstimateCatD, styleMethod[4], 1.2, colorMethod[4], colorMethod[4]);      
		histoBGEstimateCatD->DrawCopy("same,p,e1"); 
				
		TLegend* legendDiffMethods = new TLegend(0.5,0.15,0.95,0.3);
		legendDiffMethods->SetTextSize(0.03);         
		legendDiffMethods->SetFillColor(0);
		legendDiffMethods->SetLineColor(0);
		legendDiffMethods->SetNColumns(2);
		legendDiffMethods->AddEntry((TObject*)0,"Method A","");
		legendDiffMethods->AddEntry(histoBGEstimateA,"Data","p");
		if (fitCorrectionFactorsFitvsPt){
			legendDiffMethods->AddEntry((TObject*)0,"Method B","");
			legendDiffMethods->AddEntry(histoBGEstimateB,"Data","p");
		}
		legendDiffMethods->AddEntry((TObject*)0,"Method A sep Cat","");
		legendDiffMethods->AddEntry(histoBGEstimateCatA,"Data","p");
		legendDiffMethods->AddEntry((TObject*)0,"Method C sep Cat","");
		legendDiffMethods->AddEntry(histoBGEstimateCatC,"Data","p");
		legendDiffMethods->AddEntry((TObject*)0,"Method D sep Cat","");
		legendDiffMethods->AddEntry(histoBGEstimateCatD,"Data","p");   
		legendDiffMethods->Draw();
		
		TLatex *labelEnergy2 = new TLatex(0.11,0.15,Form("%s %s",intermediate.Data(),collisionSystem.Data()));
		SetStyleTLatex( labelEnergy2, 0.04,4);
		labelEnergy2->Draw();

		canvasCorrFrac->Update(); 
		canvasCorrFrac->SaveAs(Form("%s/%s_FinalCorrectionFactor_AllMethods.%s",outputDir.Data(),nameMeson.Data(),suffix.Data()));
		
		
	}   
	
	
     
	
	//**********************************************************************************
	//******************** Mass Plot *********************************************
	//**********************************************************************************
	TCanvas* canvasMass = new TCanvas("canvasMass","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasMass, 0.13, 0.02, 0.02, 0.09);
	
	if (nameMeson.CompareTo("Pi0") == 0 || nameMeson.CompareTo("Pi0EtaBinning") == 0 ){
		histoMassMeson->GetYaxis()->SetRangeUser(0.130,0.140);
	} else {
		histoMassMeson->GetYaxis()->SetRangeUser(0.54,0.56);
	}               
	histoMassMeson->GetYaxis()->SetNdivisions(510); 
	
	DrawAutoGammaMesonHistos( histoMassMeson, 
								"", "p_{T} (GeV/c)", Form("Mass for %s in |y| < %s (GeV/c^{2})",textMeson.Data(), rapidityRange.Data()), 
								kFALSE, 0., 0.7, kFALSE,
								kFALSE, 0., 0.7, 
								kFALSE, 0., 10.);
	DrawGammaSetMarker(histoMassMeson, 22, 0.8, kBlack, kBlack);                
	histoMassMeson->DrawCopy("same,e1,p"); 
	
	DrawGammaSetMarker(histoTrueMassMeson, 24, 0.8, kRed+2, kRed+2);
	histoTrueMassMeson->DrawCopy("same,e1,p"); 
	
// 	DrawGammaSetMarker(histoMassMesonLeft, 20, 0.8, kBlue, kBlue);
// 	histoMassMesonLeft->DrawCopy("same,e1,p"); 
	
	DrawGammaLines(0., maxPtMeson,mesonMassExpect, mesonMassExpect,0.1);
	
	TLegend* legendMass = new TLegend(0.15,0.12,0.5,0.25);
	legendMass->SetTextSize(0.02);         
	legendMass->SetFillColor(0);
	legendMass->SetFillStyle(0);
	legendMass->SetLineColor(0);
	legendMass->AddEntry(histoMassMeson,"standard");
	if(nameMeson.CompareTo("Pi0") == 0  || nameMeson.CompareTo("Pi0EtaBinning") == 0 ) legendMass->AddEntry(histoTrueMassMeson,"True reconstructed #pi^{0}");
	if(nameMeson.CompareTo("Eta") == 0 ) legendMass->AddEntry(histoTrueMassMeson,"True reconstructed #eta");
	
	legendMass->Draw();
	
	canvasMass->Update();
	
	canvasMass->SaveAs(Form("%s/%s_%s_Mass_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));

	canvasMass->cd();
	
	if (mode==2 || mode == 3){
		if (nameMeson.CompareTo("Pi0") == 0 || nameMeson.CompareTo("Pi0EtaBinning") == 0 ){
			histoMassMeson->GetYaxis()->SetRangeUser(0.120,0.140);
		} else {
			histoMassMeson->GetYaxis()->SetRangeUser(0.52,0.58);
		}               

		histoMassMeson->DrawCopy("e1,p"); 
		histoTrueMassMeson->DrawCopy("same,e1,p"); 
			
		TLegend* legendMass2 = new TLegend(0.55,0.12,0.95,0.25);
		legendMass2->SetTextSize(0.02);         
		legendMass2->SetFillColor(0);
		legendMass2->SetFillStyle(0);
		legendMass2->SetLineColor(0);
		legendMass2->AddEntry(histoMassMeson,"standard");
		if(nameMeson.CompareTo("Pi0") == 0  || nameMeson.CompareTo("Pi0EtaBinning") == 0 ) legendMass2->AddEntry(histoTrueMassMeson,"True reconstructed #pi^{0}");
		if(nameMeson.CompareTo("Eta") == 0 ) legendMass2->AddEntry(histoTrueMassMeson,"True reconstructed #eta");
		
		TH1D* histoTrueMassCaloPhotonMeson =          (TH1D*)fileCorrections->Get("histoTrueMassMesonCaloPhoton");
		if (histoTrueMassCaloPhotonMeson){
			DrawGammaSetMarker(histoTrueMassCaloPhotonMeson, 25, 0.8, kGreen+2, kGreen+2);
			histoTrueMassCaloPhotonMeson->DrawCopy("same,e1,p"); 
			if(nameMeson.CompareTo("Pi0") == 0  || nameMeson.CompareTo("Pi0EtaBinning") == 0 ) legendMass2->AddEntry(histoTrueMassCaloPhotonMeson,"True reconstructed #pi^{0}, cluster real #gamma");
			if(nameMeson.CompareTo("Eta") == 0 ) legendMass2->AddEntry(histoTrueMassCaloPhotonMeson,"True reconstructed #eta, cluster real #gamma");
		}
		TH1D* histoTrueMassCaloConvPhotonMeson =          (TH1D*)fileCorrections->Get("histoTrueMassMesonCaloConvPhoton");
		if (histoTrueMassCaloConvPhotonMeson){
			DrawGammaSetMarker(histoTrueMassCaloConvPhotonMeson, 25, 0.8, kCyan+2, kCyan+2);
			histoTrueMassCaloConvPhotonMeson->DrawCopy("same,e1,p"); 
			if(nameMeson.CompareTo("Pi0") == 0  || nameMeson.CompareTo("Pi0EtaBinning") == 0 ) legendMass2->AddEntry(histoTrueMassCaloConvPhotonMeson,"True reconstructed #pi^{0}, cluster conv #gamma");
			if(nameMeson.CompareTo("Eta") == 0 ) legendMass2->AddEntry(histoTrueMassCaloConvPhotonMeson,"True reconstructed #eta, cluster conv #gamma");
		}
		TH1D* histoTrueMassCaloMergedClusterMeson =          (TH1D*)fileCorrections->Get("histoTrueMassMesonCaloMergedCluster");
		if (histoTrueMassCaloMergedClusterMeson){
			DrawGammaSetMarker(histoTrueMassCaloMergedClusterMeson, 25, 0.8, kViolet+2, kViolet+2);
			histoTrueMassCaloMergedClusterMeson->DrawCopy("same,e1,p"); 
			if(nameMeson.CompareTo("Pi0") == 0  || nameMeson.CompareTo("Pi0EtaBinning") == 0 ) legendMass2->AddEntry(histoTrueMassCaloMergedClusterMeson,"True reconstructed #pi^{0}, merged cluster #gamma");
			if(nameMeson.CompareTo("Eta") == 0 ) legendMass2->AddEntry(histoTrueMassCaloMergedClusterMeson,"True reconstructed #eta, merged cluster #gamma");
		}
		
		DrawGammaLines(0., maxPtMeson,mesonMassExpect, mesonMassExpect,0.1);
		
		
		legendMass2->Draw();
		canvasMass->Update();	
		canvasMass->SaveAs(Form("%s/%s_%s_MassAddedInfos_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
	}
	delete canvasMass;
	delete legendMass;
	//**********************************************************************************
	//******************** FWHM Plot *********************************************
	//**********************************************************************************
	TCanvas* canvasFWHM = new TCanvas("canvasFWHM","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasFWHM, 0.13, 0.02, 0.04, 0.09);
	
	histoFWHMMeson->Sumw2();
	histoFWHMMeson->Scale(1./2.35);

	DrawAutoGammaMesonHistos( histoFWHMMeson, 
								"", "p_{T} (GeV/c)", Form("FWHM/2.35 for %s in |y| < %s (GeV/c^{2})",textMeson.Data(), rapidityRange.Data()), 
								kFALSE, 1.5,-20., kFALSE,
								kTRUE, -0.004, 0.020, 
								kFALSE, 0., 10.);  
	
	histoFWHMMeson->GetYaxis()->SetNdivisions(510); 
	DrawGammaSetMarker(histoFWHMMeson, 22, 0.8, kBlack, kBlack);                
	histoFWHMMeson->DrawCopy("same,e1,p"); 
	
	histoTrueFWHMMeson->Sumw2();
	histoTrueFWHMMeson->Scale(1./2.35);              
	DrawGammaSetMarker(histoTrueFWHMMeson, 24, 0.8, kRed+2, kRed+2);
	histoTrueFWHMMeson->DrawCopy("same,e1,p"); 
	
	//DrawGammaSetMarker(histoFWHMMesonLeft, 20, 0.8, kBlue, kBlue);
	//histoFWHMMesonLeft->Sumw2();
	//histoFWHMMesonLeft->Scale(1./2.35);
	//histoFWHMMesonLeft->DrawCopy("same,e1,p"); 
	
	TLegend* legendFWHM = new TLegend(0.15,0.1,0.5,0.2);
	legendFWHM->SetTextSize(0.02);         
	legendFWHM->SetFillColor(0);
	legendFWHM->AddEntry(histoFWHMMeson,"normal");
	//legendFWHM->AddEntry(histoFWHMMesonLeft,"left norm");
	if(nameMeson.CompareTo("Pi0") == 0 || nameMeson.CompareTo("Pi0EtaBinning") == 0 ) legendFWHM->AddEntry(histoTrueFWHMMeson,"True reconstructed #pi^{0}");
	if(nameMeson.CompareTo("Eta") == 0 ) legendFWHM->AddEntry(histoTrueFWHMMeson,"True reconstructed #eta");
	legendFWHM->Draw();
	canvasFWHM->Update();
	
	canvasFWHM->SaveAs(Form("%s/%s_%s_FWHM_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));

	if (mode==2 || mode == 3){
	
		histoFWHMMeson->GetYaxis()->SetRangeUser(-0.004, 0.030); 
		histoFWHMMeson->DrawCopy("e1,p"); 
		histoTrueFWHMMeson->DrawCopy("same,e1,p"); 
			
		TLegend* legendFWHM2 = new TLegend(0.55,0.12,0.95,0.25);
		legendFWHM2->SetTextSize(0.02);         
		legendFWHM2->SetFillColor(0);
		legendFWHM2->SetFillStyle(0);
		legendFWHM2->SetLineColor(0);
		legendFWHM2->AddEntry(histoFWHMMeson,"standard");
		if(nameMeson.CompareTo("Pi0") == 0  || nameMeson.CompareTo("Pi0EtaBinning") == 0 ) legendFWHM2->AddEntry(histoTrueFWHMMeson,"True reconstructed #pi^{0}");
		if(nameMeson.CompareTo("Eta") == 0 ) legendFWHM2->AddEntry(histoTrueFWHMMeson,"True reconstructed #eta");
		
		TH1D* histoTrueFWHMCaloPhotonMeson =          (TH1D*)fileCorrections->Get("histoTrueFWHMMesonCaloPhoton");
		if (histoTrueFWHMCaloPhotonMeson){
			histoTrueFWHMCaloPhotonMeson->Scale(1./2.35);
			DrawGammaSetMarker(histoTrueFWHMCaloPhotonMeson, 25, 0.8, kGreen+2, kGreen+2);
			histoTrueFWHMCaloPhotonMeson->DrawCopy("same,e1,p"); 
			if(nameMeson.CompareTo("Pi0") == 0  || nameMeson.CompareTo("Pi0EtaBinning") == 0 ) legendFWHM2->AddEntry(histoTrueFWHMCaloPhotonMeson,"True reconstructed #pi^{0}, cluster real #gamma");
			if(nameMeson.CompareTo("Eta") == 0 ) legendFWHM2->AddEntry(histoTrueFWHMCaloPhotonMeson,"True reconstructed #eta, cluster real #gamma");
		}
		TH1D* histoTrueFWHMCaloConvPhotonMeson =          (TH1D*)fileCorrections->Get("histoTrueFWHMMesonCaloConvPhoton");
		if (histoTrueFWHMCaloConvPhotonMeson){
			histoTrueFWHMCaloConvPhotonMeson->Scale(1./2.35);
			DrawGammaSetMarker(histoTrueFWHMCaloConvPhotonMeson, 25, 0.8, kCyan+2, kCyan+2);
			histoTrueFWHMCaloConvPhotonMeson->DrawCopy("same,e1,p"); 
			if(nameMeson.CompareTo("Pi0") == 0  || nameMeson.CompareTo("Pi0EtaBinning") == 0 ) legendFWHM2->AddEntry(histoTrueFWHMCaloConvPhotonMeson,"True reconstructed #pi^{0}, cluster conv #gamma");
			if(nameMeson.CompareTo("Eta") == 0 ) legendFWHM2->AddEntry(histoTrueFWHMCaloConvPhotonMeson,"True reconstructed #eta, cluster conv #gamma");
		}
		TH1D* histoTrueFWHMCaloMergedClusterMeson =          (TH1D*)fileCorrections->Get("histoTrueFWHMMesonCaloMergedCluster");
		if (histoTrueFWHMCaloMergedClusterMeson){
			histoTrueFWHMCaloMergedClusterMeson->Scale(1./2.35);
			DrawGammaSetMarker(histoTrueFWHMCaloMergedClusterMeson, 25, 0.8, kViolet+2, kViolet+2);
			histoTrueFWHMCaloMergedClusterMeson->DrawCopy("same,e1,p"); 
			if(nameMeson.CompareTo("Pi0") == 0  || nameMeson.CompareTo("Pi0EtaBinning") == 0 ) legendFWHM2->AddEntry(histoTrueFWHMCaloMergedClusterMeson,"True reconstructed #pi^{0}, merged cluster #gamma");
			if(nameMeson.CompareTo("Eta") == 0 ) legendFWHM2->AddEntry(histoTrueFWHMCaloMergedClusterMeson,"True reconstructed #eta, merged cluster #gamma");
		}
		
		
		legendFWHM2->Draw();
		canvasFWHM->Update();	
		canvasFWHM->SaveAs(Form("%s/%s_%s_FWHMAddedInfos_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
	}

	delete canvasFWHM;
	delete legendFWHM;
	
	//**************************** Fitting efficiencies ************************************************************************
	Double_t maxPtMesonEffFit;
	Double_t minPtMesonEffFit;
	Int_t offsetCorrectionHighPt;
	TF1* fitTrueEffi;
	
	if (optionEnergy.CompareTo("PbPb_2.76TeV")==0){
		maxPtMesonEffFit = maxPtMeson;
		minPtMesonEffFit = 1.2;
		offsetCorrectionHighPt= 1;
		fitTrueEffi = new TF1("EffiFitDummy","1 - [0]*exp([2]*x)+[2]");
	} else {
		maxPtMesonEffFit = 4.;
		minPtMesonEffFit = 1.2;
		if (nameMeson.Contains("Eta")){
			offsetCorrectionHighPt= -1;
		} else {
			offsetCorrectionHighPt= -2;
			cout << "doing pi0" << endl;
		}
		fitTrueEffi = new TF1("EffiFitDummy","1 - [0]*exp([2]*x)+[2]");
	}
	fitTrueEffi->SetRange(minPtMesonEffFit,maxPtMesonEffFit);
	TFitResultPtr resultEffi = histoTrueEffiPt->Fit(fitTrueEffi,"SINRME+","",minPtMesonEffFit,maxPtMesonEffFit);

	TH1D* histoTrueEffiPtFit = (TH1D*)histoTrueEffiPt->Clone("histoTrueEffiPtFit");
	
	for (Int_t i = histoTrueEffiPt->GetXaxis()->FindBin(minPtMesonEffFit)+1; i < histoTrueEffiPt->GetXaxis()->FindBin(maxPtMesonEffFit)+offsetCorrectionHighPt; i++){
		Double_t ptStart = histoTrueEffiPt->GetXaxis()->GetBinLowEdge(i);
		Double_t ptEnd = histoTrueEffiPt->GetXaxis()->GetBinUpEdge(i);
		Double_t binWidth = ptEnd-ptStart;
		Double_t effi = fitTrueEffi->Integral(ptStart, ptEnd, resultEffi->GetParams()) / binWidth;
		Double_t errorEffi = fitTrueEffi->IntegralError(ptStart, ptEnd, resultEffi->GetParams(), resultEffi->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
		histoTrueEffiPtFit->SetBinContent(i, effi);
		histoTrueEffiPtFit->SetBinError(i, errorEffi);
	}
	TH1D* histoRatioTrueEffiDivFitted = (TH1D*)histoTrueEffiPt->Clone("histoRatioTrueEffiDivFitted");
	histoRatioTrueEffiDivFitted->Divide(histoTrueEffiPt, histoTrueEffiPtFit, 1. ,1., "B");
	
	TF1* fitTrueEffiNarrow = new TF1("EffiNarrowFitDummy","1 - [0]*exp([2]*x)+[2] ");
	fitTrueEffiNarrow->SetRange(minPtMesonEffFit,maxPtMesonEffFit);
	//       fitTrueEffiNarrow->SetParameter(2,1-fitTrueEffiNarrowHighPtCut[i]->GetParameter(0)   );
	TFitResultPtr resultEffiNarrow = histoTrueEffiNarrowPt->Fit(fitTrueEffiNarrow,"SINRME+","",minPtMesonEffFit,maxPtMesonEffFit);

	TH1D* histoTrueEffiNarrowPtFit = (TH1D*)histoTrueEffiNarrowPt->Clone("histoTrueEffiNarrowPtFit");
	for (Int_t i = histoTrueEffiNarrowPt->GetXaxis()->FindBin(minPtMesonEffFit)+1; i < histoTrueEffiNarrowPt->GetXaxis()->FindBin(maxPtMesonEffFit)+offsetCorrectionHighPt; i++){
		Double_t ptStart = histoTrueEffiNarrowPt->GetXaxis()->GetBinLowEdge(i);
		Double_t ptEnd = histoTrueEffiNarrowPt->GetXaxis()->GetBinUpEdge(i);
		Double_t binWidth = ptEnd-ptStart;
		Double_t effi = fitTrueEffiNarrow->Integral(ptStart, ptEnd, resultEffiNarrow->GetParams()) / binWidth;
		Double_t errorEffiNarrow = fitTrueEffiNarrow->IntegralError(ptStart, ptEnd, resultEffiNarrow->GetParams(), resultEffiNarrow->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
		histoTrueEffiNarrowPtFit->SetBinContent(i, effi);
		histoTrueEffiNarrowPtFit->SetBinError(i, errorEffiNarrow);
	}

	TF1* fitTrueEffiWide = new TF1("EffiWideFitDummy","1 - [0]*exp([2]*x)+[2] ");
	fitTrueEffiWide->SetRange(minPtMesonEffFit,maxPtMesonEffFit);
	//       fitTrueEffiWide->SetParameter(2,1-fitTrueEffiWideHighPtCut[i]->GetParameter(0)   );
	TFitResultPtr resultEffiWide = histoTrueEffiWidePt->Fit(fitTrueEffiWide,"SINRME+","",minPtMesonEffFit,maxPtMesonEffFit);

	TH1D* histoTrueEffiWidePtFit = (TH1D*)histoTrueEffiWidePt->Clone("histoTrueEffiWidePtFit");
	for (Int_t i = histoTrueEffiWidePt->GetXaxis()->FindBin(minPtMesonEffFit)+1; i < histoTrueEffiWidePt->GetXaxis()->FindBin(maxPtMesonEffFit)+offsetCorrectionHighPt; i++){
		Double_t ptStart = histoTrueEffiWidePt->GetXaxis()->GetBinLowEdge(i);
		Double_t ptEnd = histoTrueEffiWidePt->GetXaxis()->GetBinUpEdge(i);
		Double_t binWidth = ptEnd-ptStart;
		Double_t effi = fitTrueEffiWide->Integral(ptStart, ptEnd, resultEffiWide->GetParams()) / binWidth;
		Double_t errorEffiWide = fitTrueEffiWide->IntegralError(ptStart, ptEnd, resultEffiWide->GetParams(), resultEffiWide->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
		histoTrueEffiWidePtFit->SetBinContent(i, effi);
		histoTrueEffiWidePtFit->SetBinError(i, errorEffiWide);
	}

	
	TH1D* histoTrueEffiPtUnmod 						= (TH1D*) histoTrueEffiPt->Clone("histoTrueEffiPtUnmod"); 
	TH1D* histoTrueEffiNarrowPtUnmod 				= (TH1D*) histoTrueEffiNarrowPt->Clone("histoTrueEffiNarrowPtUnmod"); 
	TH1D* histoTrueEffiWidePtUnmod 				= (TH1D*) histoTrueEffiWidePt->Clone("histoTrueEffiWidePtUnmod"); 
	
	
	//************************************Efficiencies for Dalitz-Calo *****************************************************
	  if( mode == 6 || mode == 7 ) {
	  
		
	        TH1D* histoTrueEffiPtWOWeights 			= (TH1D*)histoTrueEffiPt->Clone("histoTrueEffiPtWOWeights"); //Temporal
		TH1D* histoTrueEffiNarrowPtWOWeights		= (TH1D*)histoTrueEffiNarrowPt->Clone("histoTrueEffiNarrowPtWOWeights");
		TH1D* histoTrueEffiWidePtWOWeights		= (TH1D*)histoTrueEffiWidePt->Clone("histoTrueEffiWidePtWOWeights");
	       
		TCanvas* canvasCompEffSimple = new TCanvas("canvasCompEffSimple","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasCompEffSimple, 0.10, 0.01, 0.035, 0.09);

		TH1D* histoRatioEffWOWeightingNormalEff		= (TH1D*) histoEffiPt->Clone(); 
		histoRatioEffWOWeightingNormalEff->Divide(histoRatioEffWOWeightingNormalEff, histoTrueEffiPtWOWeights, 1., 1., "B");
		TH1D* histoRatioEffWOWeightingNormalEffNarrow	= (TH1D*) histoEffiNarrowPt->Clone(); 
		histoRatioEffWOWeightingNormalEffNarrow->Divide(histoRatioEffWOWeightingNormalEffNarrow, histoTrueEffiNarrowPtWOWeights, 1., 1., "B");
		TH1D* histoRatioEffWOWeightingNormalEffWide	= (TH1D*) histoEffiWidePt->Clone(); 
		histoRatioEffWOWeightingNormalEffWide->Divide(histoRatioEffWOWeightingNormalEffWide, histoTrueEffiWidePtWOWeights, 1., 1., "B");

		// Calculation & Plotting of correction factor for Normal integration window
		DrawAutoGammaMesonHistos( histoRatioEffWOWeightingNormalEff, 
									"", "#it{p}_{T} (GeV/#it{c})", Form("#epsilon_{eff,%s, rec}/#epsilon_{eff,%s, true wo weights} ", textMeson.Data(), textMeson.Data()), 
									kFALSE, 1.3, 3e-6, kFALSE,
									kTRUE, 0.8, 1.5, 
									kFALSE, 0., 10.);
		DrawGammaSetMarker(histoRatioEffWOWeightingNormalEff, 24, 1., 807, 807);   
		histoRatioEffWOWeightingNormalEff->Draw("e1");

		
		TF1* fitEffiBiasWOWeightsNormalPol0 		= new TF1("fitEffiBiasWOWeightsNormalPol0","[0]",0.4,maxPtMeson);
		TF1* fitEffiBiasWOWeightsNormalPol1 		= new TF1("fitEffiBiasWOWeightsNormalPol1","[0]/pow(x,[1])+[2]",0.4,maxPtMeson);
		fitEffiBiasWOWeightsNormalPol1->SetParLimits(2,0.5,1.5);
		
		histoRatioEffWOWeightingNormalEff->Fit(fitEffiBiasWOWeightsNormalPol0,"NRME+","",1.0,maxPtMeson);
		//cout << WriteParameterToFile(fitEffiBiasWOWeightsNormalPol0) << endl;
		TH1D* histoRatioEffWOWeightingNormalEffCFPol0 = (TH1D*)histoRatioEffWOWeightingNormalEff->Clone("histoRatioEffWOWeightingNormalEffCFPol0");
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(histoRatioEffWOWeightingNormalEffCFPol0);
		histoRatioEffWOWeightingNormalEffCFPol0->SetStats(kFALSE);
		histoRatioEffWOWeightingNormalEffCFPol0->SetFillColor(806);
		histoRatioEffWOWeightingNormalEffCFPol0->SetMarkerSize(0);
		histoRatioEffWOWeightingNormalEffCFPol0->Draw("e3,same");
		
		histoRatioEffWOWeightingNormalEff->Fit(fitEffiBiasWOWeightsNormalPol1,"NRME+","",1.0,maxPtMeson	);
		//cout << WriteParameterToFile(fitEffiBiasWOWeightsNormalPol1) << endl;
		TH1D* histoRatioEffWOWeightingNormalEffCFPol1 = (TH1D*)histoRatioEffWOWeightingNormalEff->Clone("histoRatioEffWOWeightingNormalEffCFPol1");
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(histoRatioEffWOWeightingNormalEffCFPol1);
		histoRatioEffWOWeightingNormalEffCFPol1->SetStats(kFALSE);
		histoRatioEffWOWeightingNormalEffCFPol1->SetFillColor(791);
		histoRatioEffWOWeightingNormalEffCFPol1->SetFillStyle(3003);
		histoRatioEffWOWeightingNormalEffCFPol1->SetMarkerSize(0);
		for (Int_t i=1; i< histoTrueEffiPt->GetNbinsX(); i++){
			if (histoTrueEffiPt->GetBinContent(i) == 0){
				histoRatioEffWOWeightingNormalEffCFPol1->SetBinContent(i,1);
				histoRatioEffWOWeightingNormalEffCFPol1->SetBinContent(i,0);
				histoRatioEffWOWeightingNormalEffCFPol0->SetBinContent(i,1);
				histoRatioEffWOWeightingNormalEffCFPol0->SetBinContent(i,0);				
			}	
		}	
		histoRatioEffWOWeightingNormalEffCFPol1->Draw("e3,same");		
		fitEffiBiasWOWeightsNormalPol0->SetLineColor(807);
		fitEffiBiasWOWeightsNormalPol0->SetLineStyle(1);
		fitEffiBiasWOWeightsNormalPol0->Draw("same");
		fitEffiBiasWOWeightsNormalPol1->SetLineColor(797);
		fitEffiBiasWOWeightsNormalPol1->SetLineStyle(7);
		fitEffiBiasWOWeightsNormalPol1->Draw("same");
		histoRatioEffWOWeightingNormalEff->Draw("e1,same");
		
		//PutProcessLabelAndEnergyOnPlot(0.72, 0.25, 28, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data(), 63, 0.03);
		
		canvasCompEffSimple->Update();
		canvasCompEffSimple->SaveAs(Form("%s/%s_EffiCompW0WeightingNormalRatio_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));		
		
		// Calculation & Plotting of correction factor for narrow integration window
		DrawGammaSetMarker(histoRatioEffWOWeightingNormalEffNarrow, 24, 1., kGreen+2, kGreen+2);   
		histoRatioEffWOWeightingNormalEffNarrow->Draw("e1,same");

		DrawAutoGammaMesonHistos( histoRatioEffWOWeightingNormalEffNarrow, 
									"", "#it{p}_{T} (GeV/#it{c})", Form("#epsilon_{eff,%s, rec}/#epsilon_{eff,%s, true wo weights} ", textMeson.Data(), textMeson.Data()), 
									kFALSE, 1.3, 3e-6, kFALSE,
									kTRUE, 0.8, 1.5, 
									kFALSE, 0., 10.);
		DrawGammaSetMarker(histoRatioEffWOWeightingNormalEffNarrow, 24, 1., kGreen+2, kGreen+2);   
		histoRatioEffWOWeightingNormalEffNarrow->Draw("e1");

		
		TF1* fitEffiBiasWOWeightsNormalPol0Nar 		= new TF1("fitEffiBiasWOWeightsNormalPol0Nar","[0]",0.4,maxPtMeson);
		TF1* fitEffiBiasWOWeightsNormalPol1Nar 		= new TF1("fitEffiBiasWOWeightsNormalPol1Nar","[0]/pow(x,[1])+[2]",0.4,maxPtMeson);
		fitEffiBiasWOWeightsNormalPol1Nar->SetParLimits(2,0.5,1.5);
		
		histoRatioEffWOWeightingNormalEffNarrow->Fit(fitEffiBiasWOWeightsNormalPol0Nar,"NRME+","",1.0,maxPtMeson);
		//cout << WriteParameterToFile(fitEffiBiasWOWeightsNormalPol0Nar) << endl;
		TH1D* histoRatioEffWOWeightingNormalEffCFPol0Nar = (TH1D*)histoRatioEffWOWeightingNormalEffNarrow->Clone("histoRatioEffWOWeightingNormalEffCFPol0Nar");
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(histoRatioEffWOWeightingNormalEffCFPol0Nar);
		histoRatioEffWOWeightingNormalEffCFPol0Nar->SetStats(kFALSE);
		histoRatioEffWOWeightingNormalEffCFPol0Nar->SetFillColor(kGreen-7);
		histoRatioEffWOWeightingNormalEffCFPol0Nar->SetMarkerSize(0);
		histoRatioEffWOWeightingNormalEffCFPol0Nar->Draw("e3,same");
		
		histoRatioEffWOWeightingNormalEffNarrow->Fit(fitEffiBiasWOWeightsNormalPol1Nar,"NRME+","",1.0,maxPtMeson	);
		//cout << WriteParameterToFile(fitEffiBiasWOWeightsNormalPol1Nar) << endl;
		TH1D* histoRatioEffWOWeightingNormalEffCFPol1Nar = (TH1D*)histoRatioEffWOWeightingNormalEffNarrow->Clone("histoRatioEffWOWeightingNormalEffCFPol1Nar");
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(histoRatioEffWOWeightingNormalEffCFPol1Nar);
		histoRatioEffWOWeightingNormalEffCFPol1Nar->SetStats(kFALSE);
		histoRatioEffWOWeightingNormalEffCFPol1Nar->SetFillColor(kGreen-6);
		histoRatioEffWOWeightingNormalEffCFPol1Nar->SetFillStyle(3003);
		histoRatioEffWOWeightingNormalEffCFPol1Nar->SetMarkerSize(0);
		for (Int_t i=1; i< histoTrueEffiPt->GetNbinsX(); i++){
			if (histoTrueEffiPt->GetBinContent(i) == 0){
				histoRatioEffWOWeightingNormalEffCFPol1Nar->SetBinContent(i,1);
				histoRatioEffWOWeightingNormalEffCFPol1Nar->SetBinContent(i,0);
				histoRatioEffWOWeightingNormalEffCFPol0Nar->SetBinContent(i,1);
				histoRatioEffWOWeightingNormalEffCFPol0Nar->SetBinContent(i,0);				
			}	
		}	
		histoRatioEffWOWeightingNormalEffCFPol1Nar->Draw("e3,same");		
		fitEffiBiasWOWeightsNormalPol0Nar->SetLineColor(kGreen+1);
		fitEffiBiasWOWeightsNormalPol0Nar->SetLineStyle(1);
		fitEffiBiasWOWeightsNormalPol0Nar->Draw("same");
		fitEffiBiasWOWeightsNormalPol1Nar->SetLineColor(kGreen+3);
		fitEffiBiasWOWeightsNormalPol1Nar->SetLineStyle(7);
		fitEffiBiasWOWeightsNormalPol1Nar->Draw("same");
		histoRatioEffWOWeightingNormalEffNarrow->Draw("e1,same");
		
		//PutProcessLabelAndEnergyOnPlot(0.72, 0.25, 28, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data(), 63, 0.03);
		
		canvasCompEffSimple->Update();
		canvasCompEffSimple->SaveAs(Form("%s/%s_EffiCompW0WeightingNormalRatioNarrow_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));		

		// Calculation & Plotting of correction factor for wide integration window
		DrawAutoGammaMesonHistos( histoRatioEffWOWeightingNormalEffWide, 
									"", "#it{p}_{T} (GeV/#it{c})", Form("#epsilon_{eff,%s, rec}/#epsilon_{eff,%s, true wo weights} ", textMeson.Data(), textMeson.Data()), 
									kFALSE, 1.3, 3e-6, kFALSE,
									kTRUE, 0.5, 1.5, 
									kFALSE, 0., 10.);
		DrawGammaSetMarker(histoRatioEffWOWeightingNormalEffWide, 24, 1., kCyan+2, kCyan+2);   
		histoRatioEffWOWeightingNormalEffWide->Draw("e1");

		
		TF1* fitEffiBiasWOWeightsNormalPol0Wi 		= new TF1("fitEffiBiasWOWeightsNormalPol0Wi","[0]",0.4,maxPtMeson);
		TF1* fitEffiBiasWOWeightsNormalPol1Wi 		= new TF1("fitEffiBiasWOWeightsNormalPol1Wi","[0]/pow(x,[1])+[2]",0.4,maxPtMeson);
		fitEffiBiasWOWeightsNormalPol1Wi->SetParLimits(2,0.5,1.5);
		
		histoRatioEffWOWeightingNormalEff->Fit(fitEffiBiasWOWeightsNormalPol0Wi,"NRME+","",0.4,maxPtMeson);
		//cout << WriteParameterToFile(fitEffiBiasWOWeightsNormalPol0Wi) << endl;
		TH1D* histoRatioEffWOWeightingNormalEffCFPol0Wi = (TH1D*)histoRatioEffWOWeightingNormalEffWide->Clone("histoRatioEffWOWeightingNormalEffCFPol0Wi");
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(histoRatioEffWOWeightingNormalEffCFPol0Wi);
		histoRatioEffWOWeightingNormalEffCFPol0Wi->SetStats(kFALSE);
		histoRatioEffWOWeightingNormalEffCFPol0Wi->SetFillColor(kCyan-7);
		histoRatioEffWOWeightingNormalEffCFPol0Wi->SetMarkerSize(0);
		histoRatioEffWOWeightingNormalEffCFPol0Wi->Draw("e3,same");
		
		histoRatioEffWOWeightingNormalEff->Fit(fitEffiBiasWOWeightsNormalPol1Wi,"NRME+","",0.4,maxPtMeson	);
		//cout << WriteParameterToFile(fitEffiBiasWOWeightsNormalPol1Wi) << endl;
		TH1D* histoRatioEffWOWeightingNormalEffCFPol1Wi = (TH1D*)histoRatioEffWOWeightingNormalEffWide->Clone("histoRatioEffWOWeightingNormalEffCFPol1Wi");
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(histoRatioEffWOWeightingNormalEffCFPol1Wi);
		histoRatioEffWOWeightingNormalEffCFPol1Wi->SetStats(kFALSE);
		histoRatioEffWOWeightingNormalEffCFPol1Wi->SetFillColor(kCyan-6);
		histoRatioEffWOWeightingNormalEffCFPol1Wi->SetFillStyle(3003);
		histoRatioEffWOWeightingNormalEffCFPol1Wi->SetMarkerSize(0);
		for (Int_t i=1; i< histoTrueEffiPt->GetNbinsX(); i++){
			if (histoTrueEffiPt->GetBinContent(i) == 0){
				histoRatioEffWOWeightingNormalEffCFPol1Wi->SetBinContent(i,1);
				histoRatioEffWOWeightingNormalEffCFPol1Wi->SetBinContent(i,0);
				histoRatioEffWOWeightingNormalEffCFPol0Wi->SetBinContent(i,1);
				histoRatioEffWOWeightingNormalEffCFPol0Wi->SetBinContent(i,0);				
			}	
		}	
		histoRatioEffWOWeightingNormalEffCFPol1Wi->Draw("e3,same");		
		fitEffiBiasWOWeightsNormalPol0Wi->SetLineColor(kCyan+1);
		fitEffiBiasWOWeightsNormalPol0Wi->SetLineStyle(1);
		fitEffiBiasWOWeightsNormalPol0Wi->Draw("same");
		fitEffiBiasWOWeightsNormalPol1Wi->SetLineColor(kCyan+3);
		fitEffiBiasWOWeightsNormalPol1Wi->SetLineStyle(7);
		fitEffiBiasWOWeightsNormalPol1Wi->Draw("same");
		histoRatioEffWOWeightingNormalEffWide->Draw("e1,same");
		
		//PutProcessLabelAndEnergyOnPlot(0.72, 0.25, 28, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data(), 63, 0.03);
		
		canvasCompEffSimple->Update();
		canvasCompEffSimple->SaveAs(Form("%s/%s_EffiCompW0WeightingNormalRatioWide_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));		
		
		
		
		if( optionEnergy.CompareTo("pPb_5.023TeV") == 0 ) {
		
		histoTrueEffiPt->Multiply(histoTrueEffiPt,histoRatioEffWOWeightingNormalEffCFPol1 );
		histoTrueEffiNarrowPt->Multiply(histoTrueEffiNarrowPt,histoRatioEffWOWeightingNormalEffCFPol1Nar );
		histoTrueEffiWidePt->Multiply(histoTrueEffiWidePt,histoRatioEffWOWeightingNormalEffCFPol1Wi );
		
		} else if( optionEnergy.CompareTo("2.76TeV") == 0 ) {
		  
		histoTrueEffiPt->Multiply(histoTrueEffiPt,histoRatioEffWOWeightingNormalEffCFPol0 );
		histoTrueEffiNarrowPt->Multiply(histoTrueEffiNarrowPt,histoRatioEffWOWeightingNormalEffCFPol0Nar );
		histoTrueEffiWidePt->Multiply(histoTrueEffiWidePt,histoRatioEffWOWeightingNormalEffCFPol0Wi );
		 
		  
		}
		
		
		
		
		// plotting of final comparison
		TH1D* histoRatioEffWOWeightingTrueEffCorr		= (TH1D*) histoEffiPt->Clone(); 
		histoRatioEffWOWeightingTrueEffCorr->Divide(histoRatioEffWOWeightingTrueEffCorr, histoTrueEffiPt, 1., 1., "B");

		histoRatioEffWOWeightingNormalEff->Draw("e1");
		DrawGammaSetMarker(histoRatioEffWOWeightingTrueEffCorr, 20, 1.5, kAzure-6, kAzure-6);   
		histoRatioEffWOWeightingTrueEffCorr->Draw("same,e1");
		
		canvasCompEffSimple->Update();
		canvasCompEffSimple->SaveAs(Form("%s/%s_EffiCompW0WeightingNormalRatioAfterFix_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));		
		
		
	  
	  
	}
	
	
	
	
	//**********************************************************************************
	//******************** Acceptance Plot *********************************************
	//**********************************************************************************
	if (isMC.CompareTo("kTRUE") ==0){
		TCanvas* canvasAcceptance = new TCanvas("canvasAcceptance","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasAcceptance, 0.13, 0.02, 0.02, 0.09);
		
		histoAcceptance->SetXTitle("p_{T} (GeV/c)");
		if (nameMeson.CompareTo("Pi0") == 0 || nameMeson.CompareTo("Pi0EtaBinning") == 0 ){
			if (optionEnergy.CompareTo("pPb_5.023TeV")==0) histoAcceptance->GetYaxis()->SetRangeUser(0.2,1.02);
			else histoAcceptance->GetYaxis()->SetRangeUser(0.7,1.02);   
		} else {
			if (optionEnergy.CompareTo("pPb_5.023TeV")==0) histoAcceptance->GetYaxis()->SetRangeUser(0.1,1.);
				else  histoAcceptance->GetYaxis()->SetRangeUser(0.5,1.);
		}  
		histoAcceptance->SetYTitle( Form("Geom. Acceptance for %s in |y| < %s",textMeson.Data(),rapidityRange.Data()));  
		histoAcceptance->GetXaxis()->SetLabelSize(0.03);
		histoAcceptance->GetYaxis()->SetLabelSize(0.03);
		histoAcceptance->GetYaxis()->SetTitleOffset(1.2);
		histoAcceptance->SetMarkerStyle(22);
		histoAcceptance->SetMarkerSize(1.);
		histoAcceptance->SetMarkerColor(kAzure-6);
		histoAcceptance->SetLineColor(kAzure-6);
		histoAcceptance->DrawCopy("e1"); 
		
		canvasAcceptance->Update();
			
		canvasAcceptance->SaveAs(Form("%s/%s_Acceptance_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
		delete canvasAcceptance;
			
		if ( (nameMeson.CompareTo("Pi0")==0 || nameMeson.CompareTo("Pi0EtaBinning")==0)){
			if (!kDalitz){
				//**********************************************************************************
				//******************** Secondary Fraction     **************************************
				//**********************************************************************************
				TCanvas* canvasSecFrac = new TCanvas("canvasSecFrac","",200,10,1350,900);  // gives the page size
				DrawGammaCanvasSettings( canvasSecFrac, 0.09, 0.02, 0.04, 0.09);
				//       canvasSecFrac->SetLogy(1); 
							
				DrawAutoGammaMesonHistos( histoYieldTrueSecFracMeson_orig, 
										"", "p_{T} (GeV/c)", "#frac{X->#pi^{0}}{#pi^{0}}", 
										kTRUE, 1.5, 0, kFALSE,
										kFALSE, 0., 0.7, 
										kFALSE, 0., 10.);
				histoYieldTrueSecFracMeson_orig->GetYaxis()->SetTitleOffset(0.9);      
				DrawGammaSetMarker(histoYieldTrueSecFracMeson_orig, 22, 1., kBlack, kBlack);                             
				DrawGammaSetMarker(histoYieldTrueSecFracFromK0SMeson_orig, 22, 1., kBlue, kBlue);                              
				histoYieldTrueSecFracMeson_orig->DrawCopy("e1");  
				histoYieldTrueSecFracFromK0SMeson_orig->DrawCopy("e1,same");  
				if (optionEnergy.CompareTo("pPb_5.023TeV")==0) {
					fitpPbSecFrac->SetLineColor(kBlack);  
					fitpPbSecFracFromK0->SetLineColor(kBlue);	
					fitpPbSecFrac->Draw("same");
					fitpPbSecFracFromK0->Draw("same");
				} else {	
					fitDefaultSecFrac->SetLineColor(kBlack);  
					fitDefaultSecFracFromK0->SetLineColor(kBlue);
					fitDefaultSecFrac->Draw("same");
					fitDefaultSecFracFromK0->Draw("same");
				}
				
				TLegend* legendSecFrac = new TLegend(0.6,0.8,0.98,0.96);
				legendSecFrac->SetTextSize(0.03);         
				legendSecFrac->SetFillColor(0);
				legendSecFrac->AddEntry(histoYieldTrueSecFracMeson,"X= All Particles");
				if (optionEnergy.CompareTo("pPb_5.023TeV")==0) legendSecFrac->AddEntry(fitpPbSecFrac,"fit on total fraction");
					else legendSecFrac->AddEntry(fitDefaultSecFrac,"fit on total fraction");
				legendSecFrac->AddEntry(histoYieldTrueSecFracFromK0SMeson,"X = K_{s}^{0}");
				if (optionEnergy.CompareTo("pPb_5.023TeV")==0) legendSecFrac->AddEntry(fitpPbSecFracFromK0,"fit on fraction from K_{s}^{0}");
					else legendSecFrac->AddEntry(fitDefaultSecFracFromK0,"fit on fraction from K_{s}^{0}");
				legendSecFrac->Draw();

				canvasSecFrac->Update();
				
				canvasSecFrac->SaveAs(Form("%s/%s_FracSecondaries_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
				delete canvasSecFrac;
			}
		}
		
		//**********************************************************************************
		//******************** Efficiency Simple Plot **************************************
		//**********************************************************************************
		TCanvas* canvasEffSimple = new TCanvas("canvasEffSimple","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasEffSimple, 0.13, 0.02, 0.03, 0.09);
		//       canvasEffSimple->SetLogy(1);  
					
		DrawAutoGammaMesonHistos( histoTrueEffiPtUnmod, 
									"", "p_{T} (GeV/c)", "True Efficiency", 
									kTRUE, 1.5, 8e-5, kFALSE,
									kFALSE, 0., 0.7, 
									kFALSE, 0., 10.);
				
		DrawGammaSetMarker(histoTrueEffiPtUnmod, 22, 1., kBlack, kBlack);   
		histoTrueEffiPtUnmod->DrawCopy("e1");   
		
		
		if ( mode == 6 || mode == 7 ){
			DrawGammaSetMarker(histoEffiPt, 25, 1., kGreen+2, kGreen+2);   
			histoEffiPt->DrawCopy("same,e1,p");
			//DrawGammaSetMarker(histoTrueEffiPtWOWeights, 24, 1., 807, 807);   
			//histoTrueEffiPtWOWeights->DrawCopy("same,e1,p");
			DrawGammaSetMarker(histoTrueEffiPt, 21, 1., kBlue+1, kBlue+1);   
			histoTrueEffiPt->DrawCopy("same,e1,p");
			
		}
		TLegend* legendEff = GetAndSetLegend2(0.25,0.13,0.45,0.24, 28);
		legendEff->SetMargin(0.15);
		
		if ( mode == 6 || mode == 7 ){
			legendEff->AddEntry(histoTrueEffiPtUnmod,"validated efficiency, w/o weights");
			//legendEff->AddEntry(histoTrueEffiPtWOWeights,"validated efficiency, w/o weights"); 
			legendEff->AddEntry(histoEffiPt,"reconstructed efficiency, as in Data"); 
			legendEff->AddEntry(histoTrueEffiPt,"corr validated efficiency"); 
			
		} else {
			legendEff->AddEntry(histoTrueEffiPtUnmod,"validated efficiency");
		}	
		legendEff->Draw();
		
		canvasEffSimple->Update();
		//PutProcessLabelAndEnergyOnPlot(0.72, 0.25, 28, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data(), 63, 0.03);
		
		
		//canvasEffSimple->Update();
		
		canvasEffSimple->SaveAs(Form("%s/%s_TrueEffSimple_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
		delete canvasEffSimple;
		

		//*********************************************************************************
		//********************** True Efficiency Plot ******************************************
		//*********************************************************************************
		
		

		TCanvas* canvasEffi = new TCanvas("canvasEffi","",1350,1200);  // gives the page size
		DrawGammaCanvasSettings( canvasEffi, 0.13, 0.02, 0.03, 0.09);

		DrawAutoGammaMesonHistos( histoTrueEffiPt, 
									"", "p_{T} (GeV/c)", "True Efficiency", 
									kFALSE, 0.75, 3e-6, kFALSE,
									kFALSE, 0., 0.7, 
									kFALSE, 0., 10.);

		DrawGammaSetMarker(histoTrueEffiPt, 22, 1., kBlack, kBlack);                               
		histoTrueEffiPt->DrawCopy("e1");    
		
		//right Side Normalization narrow
		DrawGammaSetMarker(histoTrueEffiNarrowPt, 26, 1., kGray+1, kGray+1);                               
		histoTrueEffiNarrowPt->DrawCopy("e1,same"); 
		
		//       //right Side Normalization wide
		DrawGammaSetMarker(histoTrueEffiWidePt, 26, 1., kGray+3, kGray+3);                               
		histoTrueEffiWidePt->DrawCopy("e1,same"); 
		
		if (optionEnergy.CompareTo("PbPb_2.76TeV")==0 ){         //|| optionEnergy.CompareTo("2.76TeV")==0
			fitTrueEffi->SetLineColor(kRed+2);
			fitTrueEffi->Draw("same");
			fitTrueEffiNarrow->SetLineColor(kRed-2);
			fitTrueEffiNarrow->Draw("same");
			fitTrueEffiWide->SetLineColor(kRed-4);
			fitTrueEffiWide->Draw("same");
			DrawGammaSetMarker(histoTrueEffiPtFit, 20, 1., kRed+2, kRed+2);                               
			histoTrueEffiPtFit->DrawCopy("e1,same"); 
			DrawGammaSetMarker(histoTrueEffiNarrowPtFit, 26, 1., kRed-2, kRed-2);                               
			histoTrueEffiNarrowPtFit->DrawCopy("e1,same"); 
			DrawGammaSetMarker(histoTrueEffiWidePtFit, 26, 1., kRed-4, kRed-4);                              
			histoTrueEffiWidePtFit->DrawCopy("e1,same"); 
		} 
		
		
		TLegend* legendTrueEff = new TLegend(0.6,0.13,0.98,0.24);
		legendTrueEff->SetTextSize(0.02);         
		legendTrueEff->SetFillColor(0);
		legendTrueEff->AddEntry(histoTrueEffiPt,"true normal");
		if (optionEnergy.CompareTo("PbPb_2.76TeV")==0 ) legendTrueEff->AddEntry(fitTrueEffi,"true fitted"); //|| optionEnergy.CompareTo("2.76TeV")==0 
		legendTrueEff->AddEntry(histoTrueEffiWidePt,"true wide int");
		if (optionEnergy.CompareTo("PbPb_2.76TeV")==0 )  legendTrueEff->AddEntry(fitTrueEffiWide,"true wide fitted");
		legendTrueEff->AddEntry(histoTrueEffiNarrowPt,"true narrow int");
		if (optionEnergy.CompareTo("PbPb_2.76TeV")==0 )  legendTrueEff->AddEntry(fitTrueEffiNarrow,"true narrow fitted");
		legendTrueEff->Draw();
		
		canvasEffi->Update();

		canvasEffi->SaveAs(Form("%s/%s_%s_TrueEfficiency_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
		delete legendTrueEff;
		delete canvasEffi;
	}
	
	//**********************************************************************************
	//*************************** MC Yield *********************************************
	//**********************************************************************************

	TCanvas* canvasMCYieldMeson = new TCanvas("canvasMCYieldMeson","",1350,1500);  // gives the page size
	DrawGammaCanvasSettings( canvasMCYieldMeson, 0.13, 0.02, 0.02, 0.09);
	canvasMCYieldMeson->SetLogy();

	TH1D *histoMCYieldMesonOldBin = (TH1D*)histoInputMesonOldBinPt->Clone();
	histoMCYieldMesonOldBin->SetName("MCYield_Meson_oldBin");
	ScaleMCYield(histoMCYieldMesonOldBin,  deltaRapid,  scaling,  nEvtMC,  nameMeson ,optDalitz );
	DrawAutoGammaMesonHistos( histoMCYieldMesonOldBin, 
							"", "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}",
							kFALSE, 3., 4e-10, kTRUE,
							kFALSE, 0., 0.7, 
							kFALSE, 0., 10.);
	if (histoInputMesonOldBinPtWOWeights){
		ScaleMCYield(histoInputMesonOldBinPtWOWeights,  deltaRapid,  scaling,  nEvtMC,  nameMeson ,optDalitz );
		histoInputMesonOldBinPtWOWeights->SetName("MCYield_Meson_oldBinWOWeights");
	}   
	if (histoMCInputAddedSig){
		ScaleMCYield(histoMCInputAddedSig,  deltaRapid,  scaling,  nEvtMC,  nameMeson ,optDalitz );
		histoMCInputAddedSig->SetName("MCYield_Meson_oldBin_AddedSig");
	}   
	if (histoMCInputWOWeightingAddedSig){
		ScaleMCYield(histoMCInputWOWeightingAddedSig,  deltaRapid,  scaling,  nEvtMC,  nameMeson ,optDalitz );
		histoMCInputWOWeightingAddedSig->SetName("MCYield_Meson_oldBinWOWeights_AddedSig");
	}   

	TH1D *histoMCYieldMeson = (TH1D*)histoInputMesonPt->Clone();
	histoMCYieldMeson->SetName("MCYield_Meson");
	ScaleMCYield(histoMCYieldMeson,  deltaRapid,  scaling,  nEvtMC,  nameMeson ,optDalitz );
	DrawAutoGammaMesonHistos(histoMCYieldMeson , 
							"", "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}",
							kFALSE, 3., 4e-10, kTRUE,
							kFALSE, 0., 0.7, 
							kFALSE, 0., histoUnCorrectedYield->GetXaxis()->GetBinUpEdge(histoUnCorrectedYield->GetNbinsX()));

	TF1* fitTsallisMC;
	if (nameMeson.CompareTo("Pi0")==0 || nameMeson.CompareTo("Pi0EtaBinning")==0 ){
		fitTsallisMC= FitObject("l","fitTsallisMC","Pi0",histoMCYieldMesonOldBin,0.3,histoUnCorrectedYield->GetXaxis()->GetBinUpEdge(histoUnCorrectedYield->GetNbinsX()),NULL,"QNRME+");
	} else { 
		fitTsallisMC= FitObject("l","fitTsallisMC","Eta",histoMCYieldMesonOldBin,0.3,histoUnCorrectedYield->GetXaxis()->GetBinUpEdge(histoUnCorrectedYield->GetNbinsX()),NULL,"QNRME+");
	}
	DrawGammaSetMarkerTF1(fitTsallisMC, 1, 1.5, kBlue);
	TString forOutput= WriteParameterToFile(fitTsallisMC);
	cout << forOutput.Data()<< endl;
	histoMCYieldMesonOldBin->Draw("same");
	fitTsallisMC->Draw("same");
	
	canvasMCYieldMeson->SaveAs(Form("%s/%s_histoMCYieldMeson_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
	delete canvasMCYieldMeson;             
	

	//**********************************************************************************
	//******************** InvMass Plot full pt range **********************************
	//**********************************************************************************
	TCanvas* canvasInvMassFull = new TCanvas("canvasInvMassFull","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasInvMassFull, 0.13, 0.02, 0.06, 0.09);

	Double_t widthMeV = histoMesonSignalFullPtInvMass->GetBinWidth(10)*1000;
	cout << widthMeV << endl ;
	DrawAutoGammaMesonHistos( histoMesonSignalFullPtInvMass, 
								"", "M_{#gamma#gamma} (GeV/c^{2})", Form("Counts/ %2.1e MeV/c^{2}",widthMeV), 
								kTRUE, 1.1, 0., kFALSE,
								kFALSE, 0., 0.7, 
								kFALSE, 0., 10.);

	DrawGammaSetMarker(histoMesonSignalFullPtInvMass, 20, 0.4, kBlack, kBlack);                               
	histoMesonSignalFullPtInvMass->DrawCopy("e1");  
	histoMesonBckNormFullPtInvMass->SetLineColor(kBlue);
	histoMesonBckNormFullPtInvMass->SetLineWidth(0.99);
	histoMesonBckNormFullPtInvMass->DrawCopy("same,hist");

	TLatex* textPi0InvMass;
	if (nameMeson.CompareTo("Pi0") == 0  || nameMeson.CompareTo("Pi0EtaBinning") == 0 ){
		textPi0InvMass = new TLatex(0.44,0.8+pictDrawingCoordinatesInv[3],"p_{T}^{#gamma#gamma} > 0.4 GeV/c");
	} else {
		textPi0InvMass = new TLatex(0.44,0.8+pictDrawingCoordinatesInv[3],"p_{T}^{#gamma#gamma} > 0.6 GeV/c");
	}
	textPi0InvMass->SetNDC();
	textPi0InvMass->SetTextColor(1);
	textPi0InvMass->SetTextSize(pictDrawingCoordinatesInv[7]);
	textPi0InvMass->Draw();
	canvasInvMassFull->Update();

	canvasInvMassFull->SaveAs(Form("%s/%s_%s_InvMassFullPt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
	delete canvasInvMassFull;
		
	TCanvas* canvasInvMassFullWO = new TCanvas("canvasInvMassFullWO","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasInvMassFullWO, 0.13, 0.02, 0.06, 0.09);  
	DrawAutoGammaMesonHistos( histoMesonSignalFullPtInvMass, 
								"", "M_{#gamma#gamma} (GeV/c^{2})", Form("Counts/ %2.1e MeV/c^{2}",widthMeV), 
								kTRUE, 1.1, 0., kFALSE,
								kFALSE, 0., 0.7, 
								kFALSE, 0., 10.);

	DrawGammaSetMarker(histoMesonSignalFullPtInvMass, 20, 0.4, kBlack, kBlack);                               
	histoMesonSignalFullPtInvMass->DrawCopy("e1");  
		
	if (nameMeson.CompareTo("Pi0") == 0  || nameMeson.CompareTo("Pi0EtaBinning") == 0 ){
		textPi0InvMass = new TLatex(0.44,0.8+pictDrawingCoordinatesInv[3],"p_{T}^{#gamma#gamma} > 0.4 GeV/c");
	} else {
		textPi0InvMass = new TLatex(0.44,0.8+pictDrawingCoordinatesInv[3],"p_{T}^{#gamma#gamma} > 0.6 GeV/c");
	}
	textPi0InvMass->SetNDC();
	textPi0InvMass->SetTextColor(1);
	textPi0InvMass->SetTextSize(pictDrawingCoordinatesInv[7]);
	textPi0InvMass->Draw();

	canvasInvMassFullWO->Update();

	canvasInvMassFullWO->SaveAs(Form("%s/%s_%s_InvMassFullPtWithoutBG_%s.%s",outputDir.Data(), nameMeson.Data(), prefix2.Data(), fCutSelection.Data(), suffix.Data()));
	delete canvasInvMassFullWO;



	//**********************************************************************************
	//******************** RAW Yield spectrum ******************************************
	//**********************************************************************************
		
	TCanvas* canvasRAWYield = new TCanvas("canvasRAWYield","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRAWYield, 0.13, 0.02, 0.02, 0.09); 
	canvasRAWYield->SetLogy(1);         

	TH1D* histoUnCorrectedYieldDrawing = (TH1D*)histoUnCorrectedYield->Clone();
	histoUnCorrectedYieldDrawing->Scale(1./nEvt);
	DrawAutoGammaMesonHistos( histoUnCorrectedYieldDrawing, 
								"", "p_{T} (GeV/c)", "RAW Yield/ N_{Evt}", 
								kTRUE, 3., 4e-10, kTRUE,
								kFALSE, 0., 0.7, 
								kFALSE, 0., 10.);
	histoUnCorrectedYieldDrawing->SetLineWidth(0.5);                
	DrawGammaSetMarker(histoUnCorrectedYieldDrawing, 20, 0.5, kBlack, kBlack);                             
	histoUnCorrectedYieldDrawing->DrawCopy("e1");   

	canvasRAWYield->Update();

	canvasRAWYield->SaveAs(Form("%s/%s_%s_RAWYieldPt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
	delete canvasRAWYield;

	//***********************************************************************************************
	//***************************  correction for yield with True effi ******************************
	//***********************************************************************************************
		
	
	TH1D* histoCorrectedYieldNorm = (TH1D*)histoUnCorrectedYield->Clone();
	histoCorrectedYieldNorm->SetName("CorrectedYieldNormEff");
	//    TH1D* histoCorrectedYieldNormBackFit = (TH1D*)histoUnCorrectedYieldBackFit->Clone();
	//    histoCorrectedYieldNormBackFit->SetName("CorrectedYieldNormEffBackFit");
	TH1D* histoCorrectedYieldTrue = (TH1D*)histoUnCorrectedYield->Clone();
	histoCorrectedYieldTrue->SetName("CorrectedYieldTrueEff");
	//    TH1D* histoCorrectedYieldTrueBackFit = (TH1D*)histoUnCorrectedYieldBackFit->Clone();
	//    histoCorrectedYieldTrueBackFit->SetName("CorrectedYieldTrueEffBackFit");
	TH1D* histoCorrectedYieldTrueNarrow = (TH1D*)histoUnCorrectedYieldNarrow->Clone();
	histoCorrectedYieldTrueNarrow->SetName("CorrectedYieldTrueEffNarrow");
	TH1D* histoCorrectedYieldTrueWide = (TH1D*)histoUnCorrectedYieldWide->Clone();
	histoCorrectedYieldTrueWide->SetName("CorrectedYieldTrueEffWide");
	TH1D* histoCorrectedYieldFixed = (TH1D*)histoUnCorrectedYield->Clone();
	histoCorrectedYieldFixed->SetName("CorrectedYieldEffFixed");
	TH1D* histoCorrectedYieldNarrowFixed = (TH1D*)histoUnCorrectedYieldNarrow->Clone();
	histoCorrectedYieldNarrowFixed->SetName("CorrectedYieldEffNarrowFixed");
	TH1D* histoCorrectedYieldWideFixed = (TH1D*)histoUnCorrectedYieldWide->Clone();
	histoCorrectedYieldWideFixed->SetName("CorrectedYieldEffWideFixed");
	TH1D* histoCorrectedYieldTrueFixed = (TH1D*)histoUnCorrectedYield->Clone();
	histoCorrectedYieldTrueFixed->SetName("CorrectedYieldTrueEffFixed");
	TH1D* histoCorrectedYieldTrueNarrowFixed = (TH1D*)histoUnCorrectedYieldNarrow->Clone();
	histoCorrectedYieldTrueNarrowFixed->SetName("CorrectedYieldTrueEffNarrowFixed");
	TH1D* histoCorrectedYieldTrueWideFixed = (TH1D*)histoUnCorrectedYieldWide->Clone();
	histoCorrectedYieldTrueWideFixed->SetName("CorrectedYieldTrueEffWideFixed");
	TH1D* histoCorrectedYieldTrueLeft = (TH1D*)histoUnCorrectedYieldLeft->Clone();
	histoCorrectedYieldTrueLeft->SetName("CorrectedYieldTrueEffLeft");
	TH1D* histoCorrectedYieldTrueLeftNarrow = (TH1D*)histoUnCorrectedYieldLeftNarrow->Clone();
	histoCorrectedYieldTrueLeftNarrow->SetName("CorrectedYieldTrueEffLeftNarrow");
	TH1D* histoCorrectedYieldTrueLeftWide = (TH1D*)histoUnCorrectedYieldLeftWide->Clone();
	histoCorrectedYieldTrueLeftWide->SetName("CorrectedYieldTrueEffLeftWide");
	TH1D* histoCorrectedYieldTrueFitted = (TH1D*)histoUnCorrectedYield->Clone();
	histoCorrectedYieldTrueFitted->SetName("CorrectedYieldTrueEffFitted");
	TH1D* histoCorrectedYieldTrueNarrowFitted = (TH1D*)histoUnCorrectedYieldNarrow->Clone();
	histoCorrectedYieldTrueNarrowFitted->SetName("CorrectedYieldTrueEffNarrowFitted");
	TH1D* histoCorrectedYieldTrueWideFitted = (TH1D*)histoUnCorrectedYieldWide->Clone();
	histoCorrectedYieldTrueWideFitted->SetName("CorrectedYieldTrueEffWideFitted");
	TH1D* histoCorrectedYieldTrueLeftFitted = (TH1D*)histoUnCorrectedYieldLeft->Clone();
	histoCorrectedYieldTrueLeftFitted->SetName("CorrectedYieldTrueEffLeftFitted");
	TH1D* histoCorrectedYieldTrueLeftNarrowFitted = (TH1D*)histoUnCorrectedYieldLeftNarrow->Clone();
	histoCorrectedYieldTrueLeftNarrowFitted->SetName("CorrectedYieldTrueEffLeftNarrowFitted");
	TH1D* histoCorrectedYieldTrueLeftWideFitted = (TH1D*)histoUnCorrectedYieldLeftWide->Clone();
	histoCorrectedYieldTrueLeftWideFitted->SetName("CorrectedYieldTrueEffLeftWideFitted");
	
	TH1D* histoCompleteCorr = (TH1D*)histoTrueEffiPt->Clone();

	if (!optDalitz){
		CorrectYield(histoCorrectedYieldNorm, histoYieldSecMeson, histoYieldSecFromK0SMeson, histoEffiPt, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
		
		CorrectYield(histoCorrectedYieldTrue, histoYieldSecMeson, histoYieldSecFromK0SMeson, histoTrueEffiPt, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
		CompileFullCorrectionFactor( histoCompleteCorr, histoAcceptance, deltaRapid);
	//       CorrectYield(histoCorrectedYieldNormBackFit, histoYieldSecMesonBackFit, histoYieldSecFromK0SMesonBackFit, histoEffiPtBackFit, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
	//       CorrectYield(histoCorrectedYieldTrueBackFit, histoYieldSecMesonBackFit, histoYieldSecFromK0SMesonBackFit, histoTrueEffiPt, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
		CorrectYield(histoCorrectedYieldTrueNarrow, histoYieldSecMesonNarrow, histoYieldSecFromK0SMesonNarrow , histoTrueEffiNarrowPt, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
		CorrectYield(histoCorrectedYieldTrueWide, histoYieldSecMesonWide, histoYieldSecFromK0SMesonWide, histoTrueEffiWidePt, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
		
		CorrectYield(histoCorrectedYieldFixed, histoYieldSecMeson, histoYieldSecFromK0SMeson, histoEffiPtFixed, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
		CorrectYield(histoCorrectedYieldNarrowFixed, histoYieldSecMesonNarrow, histoYieldSecFromK0SMesonNarrow , histoEffiNarrowPtFixed, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
		CorrectYield(histoCorrectedYieldWideFixed, histoYieldSecMesonWide, histoYieldSecFromK0SMesonWide, histoEffiWidePtFixed, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);

		CorrectYield(histoCorrectedYieldTrueFixed, histoYieldSecMeson, histoYieldSecFromK0SMeson, histoTrueEffiPtFixed, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
		CorrectYield(histoCorrectedYieldTrueNarrowFixed, histoYieldSecMesonNarrow, histoYieldSecFromK0SMesonNarrow , histoTrueEffiNarrowPtFixed, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
		CorrectYield(histoCorrectedYieldTrueWideFixed, histoYieldSecMesonWide, histoYieldSecFromK0SMesonWide, histoTrueEffiWidePtFixed, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
	
		CorrectYield(histoCorrectedYieldTrueLeft, histoYieldSecMesonLeft, histoYieldSecFromK0SMesonLeft, histoTrueEffiPt, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
		CorrectYield(histoCorrectedYieldTrueLeftNarrow, histoYieldSecMesonLeftNarrow, histoYieldSecFromK0SMesonLeftNarrow, histoTrueEffiNarrowPt, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
		CorrectYield(histoCorrectedYieldTrueLeftWide, histoYieldSecMesonLeftWide, histoYieldSecFromK0SMesonLeftWide, histoTrueEffiWidePt, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);

		CorrectYield(histoCorrectedYieldTrueFitted, histoYieldSecMeson, histoYieldSecFromK0SMeson, histoTrueEffiPtFit, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
		CorrectYield(histoCorrectedYieldTrueNarrowFitted, histoYieldSecMesonNarrow, histoYieldSecFromK0SMesonNarrow , histoTrueEffiNarrowPtFit, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
		CorrectYield(histoCorrectedYieldTrueWideFitted, histoYieldSecMesonWide, histoYieldSecFromK0SMesonWide, histoTrueEffiWidePtFit, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
		CorrectYield(histoCorrectedYieldTrueLeftFitted, histoYieldSecMesonLeft, histoYieldSecFromK0SMesonLeft, histoTrueEffiPtFit, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
		CorrectYield(histoCorrectedYieldTrueLeftNarrowFitted, histoYieldSecMesonLeftNarrow, histoYieldSecFromK0SMesonLeftNarrow, histoTrueEffiNarrowPtFit, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
		CorrectYield(histoCorrectedYieldTrueLeftWideFitted, histoYieldSecMesonLeftWide, histoYieldSecFromK0SMesonLeftWide, histoTrueEffiWidePtFit, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
		
		if (isMC.CompareTo("kFALSE") ==0 && kDCAFileDataExists){
	//          histoCorrectedYieldNorm->Multiply(histoBGEstimateA);
				histoCorrectedYieldTrue->Multiply(histoBGEstimateCatA);
				histoCorrectedYieldTrueNarrow->Multiply(histoBGEstimateCatA);
				histoCorrectedYieldTrueWide->Multiply(histoBGEstimateCatA);
				histoCorrectedYieldFixed->Multiply(histoBGEstimateCatA);
				histoCorrectedYieldNarrowFixed->Multiply(histoBGEstimateCatA);
				histoCorrectedYieldWideFixed->Multiply(histoBGEstimateCatA);
				histoCorrectedYieldTrueFixed->Multiply(histoBGEstimateCatA);
				histoCorrectedYieldTrueNarrowFixed->Multiply(histoBGEstimateCatA);
				histoCorrectedYieldTrueWideFixed->Multiply(histoBGEstimateCatA);
				histoCorrectedYieldTrueLeft->Multiply(histoBGEstimateCatA);
				histoCorrectedYieldTrueLeftNarrow->Multiply(histoBGEstimateCatA);
				histoCorrectedYieldTrueLeftWide->Multiply(histoBGEstimateCatA);
				histoCorrectedYieldTrueFitted->Multiply(histoBGEstimateCatA);
				histoCorrectedYieldTrueNarrowFitted->Multiply(histoBGEstimateCatA);
				histoCorrectedYieldTrueWideFitted->Multiply(histoBGEstimateCatA);
				histoCorrectedYieldTrueLeftFitted->Multiply(histoBGEstimateCatA);
				histoCorrectedYieldTrueLeftNarrowFitted->Multiply(histoBGEstimateCatA);
				histoCorrectedYieldTrueLeftWideFitted->Multiply(histoBGEstimateCatA);
		}   

	} else {
		CorrectYieldDalitz(histoCorrectedYieldNorm, histoYieldGGMeson,histoEffiPt,     histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
		CorrectYieldDalitz(histoCorrectedYieldTrue, histoYieldGGMeson,histoTrueEffiPt, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
		CorrectYieldDalitz(histoCorrectedYieldTrueNarrow, histoYieldGGMesonNarrow, histoTrueEffiNarrowPt, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
		CorrectYieldDalitz(histoCorrectedYieldTrueWide, histoYieldGGMesonWide, histoTrueEffiWidePt, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);

		CorrectYieldDalitz(histoCorrectedYieldTrueLeft, histoYieldGGMesonLeft, histoTrueEffiPt, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
		CorrectYieldDalitz(histoCorrectedYieldTrueLeftNarrow, histoYieldGGMesonLeftNarrow,histoTrueEffiNarrowPt, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
		CorrectYieldDalitz(histoCorrectedYieldTrueLeftWide, histoYieldGGMesonLeftWide, histoTrueEffiWidePt, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
	}
	
	cout << "here" << endl;
	TCanvas* canvasCorrecftedYield = new TCanvas("canvasCorrecftedYield","",1350,1500);  // gives the page size
	DrawGammaCanvasSettings( canvasCorrecftedYield, 0.13, 0.02, 0.02, 0.09);   
	canvasCorrecftedYield->SetLogy();   

	TPad* padCorrectedYieldHistos = new TPad("padCorrectedYieldHistos", "", 0., 0.25, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( padCorrectedYieldHistos, 0.12, 0.02, 0.02, 0.);
	padCorrectedYieldHistos->Draw();

	TPad* padCorrectedYieldRatios = new TPad("padCorrectedYieldRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
	DrawGammaPadSettings( padCorrectedYieldRatios, 0.12, 0.02, 0., 0.2);
	padCorrectedYieldRatios->Draw();

	padCorrectedYieldHistos->cd();
	padCorrectedYieldHistos->SetLogy();    

	cout << "here" << endl;
	DrawAutoGammaMesonHistos( histoCorrectedYieldTrue, 
								"", "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}", 
								kTRUE, 3., 4e-10, kTRUE,
								kFALSE, 0., 0.7, 
								kFALSE, 0., 10.);
	cout << "here" << endl; 
	DrawGammaSetMarker(histoCorrectedYieldTrue, 22, 1., kBlack, kBlack);                             
	histoCorrectedYieldTrue->DrawCopy("e1");  
	//right Side Normalization narrow
	cout << "here" << endl;
	DrawGammaSetMarker(histoCorrectedYieldTrueNarrow, 26, 1., kGray+1, kGray+1);                                
	histoCorrectedYieldTrueNarrow->DrawCopy("e1,same"); 
	//right Side Normalization wide
	cout << "here" << endl;
	DrawGammaSetMarker(histoCorrectedYieldTrueWide, 26, 1., kGray+3, kGray+3);                             
	histoCorrectedYieldTrueWide->DrawCopy("e1,same"); 
	cout << "here" << endl;
	//left Side Normalization 
	DrawGammaSetMarker(histoCorrectedYieldTrueLeft, 20, 1., kBlue, kBlue);                              
	histoCorrectedYieldTrueLeft->DrawCopy("e1,same"); 
	cout << "here" << endl;
	//left Side Normalization narrow
	DrawGammaSetMarker(histoCorrectedYieldTrueLeftNarrow, 24, 1., kBlue-5, kBlue-5);                      
	histoCorrectedYieldTrueLeftNarrow->DrawCopy("e1,same"); 
	cout << "here" << endl;
	//left Side Normalization wide
	DrawGammaSetMarker(histoCorrectedYieldTrueLeftWide, 24, 1., kBlue+2, kBlue+2);                               
	histoCorrectedYieldTrueLeftWide->DrawCopy("e1,same");    
	
	cout << "here" << endl; 
	TLegend* legendYield3 = new TLegend(0.15,0.03,0.66,0.19);
	legendYield3->SetTextSize(0.02);       
	legendYield3->SetFillColor(0);
	legendYield3->AddEntry(histoCorrectedYieldTrue,"corr true eff/right norm");
	legendYield3->AddEntry(histoCorrectedYieldTrueWide,"corr true eff wide int /right norm");
	legendYield3->AddEntry(histoCorrectedYieldTrueNarrow,"corr true eff narrow int /right norm");
	legendYield3->AddEntry(histoCorrectedYieldTrueLeft,Form("corr true eff /left norm"));
	legendYield3->AddEntry(histoCorrectedYieldTrueLeftWide,"corr true eff wide int /left norm");
	legendYield3->AddEntry(histoCorrectedYieldTrueLeftNarrow,"corr true eff narrow int /left norm");
	legendYield3->Draw();

	padCorrectedYieldRatios->cd();
	padCorrectedYieldRatios->SetTickx();
	padCorrectedYieldRatios->SetTicky();
	padCorrectedYieldRatios->SetLogy(0);
	TH1D *RatioTrue = (TH1D*) histoCorrectedYieldTrue->Clone(); 
	RatioTrue->Divide(RatioTrue,histoCorrectedYieldTrue,1.,1.,"");
	TH1D *RatioTrueWide = (TH1D*) histoCorrectedYieldTrue->Clone();   
	RatioTrueWide->Divide(RatioTrueWide,histoCorrectedYieldTrueWide,1.,1.,"");
	TH1D *RatioTrueNarrow = (TH1D*) histoCorrectedYieldTrue->Clone(); 
	RatioTrueNarrow->Divide(RatioTrueNarrow,histoCorrectedYieldTrueNarrow,1.,1.,"");
	TH1D *RatioTrueLeft = (TH1D*) histoCorrectedYieldTrue->Clone();   
	RatioTrueLeft->Divide(RatioTrueLeft,histoCorrectedYieldTrueLeft,1.,1.,"");
	TH1D *RatioTrueLeftWide = (TH1D*) histoCorrectedYieldTrue->Clone();  
	RatioTrueLeftWide->Divide(RatioTrueLeftWide,histoCorrectedYieldTrueLeftWide,1.,1.,"");
	TH1D *RatioTrueLeftNarrow = (TH1D*) histoCorrectedYieldTrue->Clone();   
	RatioTrueLeftNarrow->Divide(RatioTrueLeftNarrow,histoCorrectedYieldTrueLeftNarrow,1.,1.,"");
	TH1D *RatioTrueMCInput = (TH1D*) histoCorrectedYieldTrue->Clone();   
	RatioTrueMCInput->Divide(RatioTrueMCInput,histoMCYieldMeson,1.,1.,"");
	TH1D *RatioNormal = (TH1D*) histoCorrectedYieldTrue->Clone();  
	RatioNormal->Divide(RatioNormal,histoCorrectedYieldNorm,1.,1.,"");

	RatioTrue->SetYTitle("#frac{standard}{modified}"); 
	RatioTrue->GetYaxis()->SetRangeUser(0.8,1.23);
	RatioTrue->GetYaxis()->SetLabelSize(0.07);
	RatioTrue->GetYaxis()->SetNdivisions(505);
	RatioTrue->GetYaxis()->SetTitleSize(0.1); 
	RatioTrue->GetYaxis()->SetDecimals();
	RatioTrue->GetYaxis()->SetTitleOffset(0.5);
	RatioTrue->GetXaxis()->SetTitleSize(0.11);   
	RatioTrue->GetXaxis()->SetLabelSize(0.08);
	RatioTrue->SetMarkerStyle(22);
	RatioTrue->SetMarkerSize(1.);
	RatioTrue->SetMarkerColor(kBlack);
	RatioTrue->SetLineColor(kBlack);
	RatioTrue->DrawCopy("p,e1"); 

	DrawGammaSetMarker(RatioTrue, 22, 1., kBlack, kBlack);                               
	RatioTrue->DrawCopy("p,e1");  
	//right Side Normalization narrow
	DrawGammaSetMarker(RatioTrueNarrow, 26, 1., kGray+1, kGray+1);                               
	RatioTrueNarrow->DrawCopy("e1,same"); 
	//right Side Normalization wide
	DrawGammaSetMarker(RatioTrueWide, 26, 1., kGray+3, kGray+3);                               
	RatioTrueWide->DrawCopy("e1,same"); 
	//left Side Normalization 
	DrawGammaSetMarker(RatioTrueLeft, 20, 1., kBlue, kBlue);                             
	RatioTrueLeft->DrawCopy("e1,same"); 
	//left Side Normalization narrow
	DrawGammaSetMarker(RatioTrueLeftNarrow, 24, 1., kBlue-5, kBlue-5);                        
	RatioTrueLeftNarrow->DrawCopy("e1,same"); 
	//left Side Normalization wide
	DrawGammaSetMarker(RatioTrueLeftWide, 24, 1., kBlue+2, kBlue+2);                              
	RatioTrueLeftWide->DrawCopy("e1,same"); 

	canvasCorrecftedYield->Update();
	canvasCorrecftedYield->SaveAs(Form("%s/%s_%s_CorrectedYieldTrueEff_%s.%s",outputDir.Data(), nameMeson.Data(), prefix2.Data(),  fCutSelection.Data(), suffix.Data()));

	
	
	if (isMC.CompareTo("kTRUE") ==0){
		canvasCorrecftedYield->cd();   

		padCorrectedYieldHistos->cd();
		padCorrectedYieldHistos->SetLogy();    

		DrawGammaSetMarker(histoCorrectedYieldTrue, 22, 1.2, kBlack, kBlack);                             
		histoCorrectedYieldTrue->DrawCopy("e1");  
		DrawGammaSetMarker(histoMCYieldMeson, 24, 1.2, kRed+2, kRed+2);                             
		histoMCYieldMeson->DrawCopy("e1,same"); 
		DrawGammaSetMarker(histoCorrectedYieldNorm, 25, 1.2, kGreen+2, kGreen+2);                               
		histoCorrectedYieldNorm->DrawCopy("e1,same"); 
		
		
		
		cout << "here" << endl; 
		TLegend* legendYield3 = new TLegend(0.15,0.03,0.66,0.19);
		legendYield3->SetTextSize(0.02);       
		legendYield3->SetFillColor(0);
		legendYield3->AddEntry(histoCorrectedYieldTrue,"corr true eff");
		legendYield3->AddEntry(histoMCYieldMeson,"MC input");
		legendYield3->AddEntry(histoCorrectedYieldNorm,"normal eff");
		legendYield3->Draw();

		padCorrectedYieldRatios->cd();
		padCorrectedYieldRatios->SetTickx();
		padCorrectedYieldRatios->SetTicky();
		padCorrectedYieldRatios->SetLogy(0);
		
		DrawGammaSetMarker(RatioTrue, 22, 1.2, kBlack, kBlack);                               
		RatioTrue->DrawCopy("p,e1");  
		DrawGammaSetMarker(RatioTrueMCInput, 24, 1.2, kRed+2, kRed+2);                              
		RatioTrueMCInput->DrawCopy("e1,same"); 
		DrawGammaSetMarker(RatioNormal, 25, 1.2, kGreen+2, kGreen+2);                               
		RatioNormal->DrawCopy("e1,same"); 

		canvasCorrecftedYield->Update();
		canvasCorrecftedYield->SaveAs(Form("%s/%s_%s_CorrectedYield_SanityCheck_%s.%s",outputDir.Data(), nameMeson.Data(), prefix2.Data(),  fCutSelection.Data(), suffix.Data()));
			
		
	}

	delete canvasCorrecftedYield;
	delete legendYield3;

	//***********************************************************************************************
	//***************************  Secondary RAW Yield  *********************************************
	//***********************************************************************************************
	if (!optDalitz && doubleAddFactorK0s >= 0){
		TCanvas* canvasRAWYieldSec = new TCanvas("canvasRAWYieldSec","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasRAWYieldSec, 0.13, 0.02, 0.02, 0.09); 
		canvasRAWYieldSec->SetLogy(1);         

		DrawGammaSetMarker(histoUnCorrectedYieldDrawing, 20, 1., kBlack, kBlack);                              
		histoUnCorrectedYieldDrawing->Draw("e1");
		histoYieldSecMeson->Scale(1./nEvt);
		DrawGammaSetMarker(histoYieldSecMeson, 22, 1., kBlue, kBlue);                              
		histoYieldSecMeson->DrawCopy("same,e1");  
		histoYieldSecFromK0SMeson->Scale(1./nEvt);
		DrawGammaSetMarker(histoYieldSecFromK0SMeson, 21, 1., kCyan, kCyan);                             
		histoYieldSecFromK0SMeson->DrawCopy("same,e1");    

		TLegend* legendSecRAWYield = new TLegend(0.6,0.8,0.97,0.95);
		legendSecRAWYield->SetTextSize(0.03);        
		legendSecRAWYield->SetFillColor(0);
		legendSecRAWYield->SetBorderSize(0);
		legendSecRAWYield->AddEntry(histoUnCorrectedYieldDrawing,"RAW yield");
		legendSecRAWYield->AddEntry(histoYieldSecMeson,"total secondaries");
		legendSecRAWYield->AddEntry(histoYieldSecFromK0SMeson,"additional secondaries from K^{0}_{s}");
		legendSecRAWYield->Draw();

		canvasRAWYieldSec->Update();
		canvasRAWYieldSec->SaveAs(Form("%s/%s_%s_RAWYieldSecPt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
		delete canvasRAWYieldSec;
	} else if (optDalitz){
	  
	  
		if( isMC.CompareTo("kTRUE") == 0 ) {
	  
		cout<<"***************************Computing Ratio and Fitting***************************************************"<<endl;
     
     
		TH1D *RatioTrueMCInput02 = (TH1D*) histoCorrectedYieldTrue->Clone("RatioTrueMCInput02");  
		      RatioTrueMCInput02->Divide(RatioTrueMCInput02,histoMCYieldMeson,1.,1.,"");
   
		Double_t minPtRatioTrueMCInput = 0.05;
		Double_t maxPtRatioTrueMCInput = 10.0;
   
  
   		TF1* fitRatioTrueMCInput = new TF1("fitRatioTrueMCInput","pol0");
   		fitRatioTrueMCInput->SetRange(minPtRatioTrueMCInput,maxPtRatioTrueMCInput);
		TFitResultPtr resultFitRatioTrueMCInput = RatioTrueMCInput02->Fit(fitRatioTrueMCInput,"SINRME+","",minPtRatioTrueMCInput,maxPtRatioTrueMCInput);
   
		TPaveText* fitParamRatioTrueMCInput = new TPaveText(6.0,0.03,9.6,0.29);
   
		// TPaveText *pt = new TPaveText(.05,.1,.95,.8);
		fitParamRatioTrueMCInput->AddText(Form("chi2   %f",fitRatioTrueMCInput->GetChisquare()));
		fitParamRatioTrueMCInput->AddText(Form("Param0 %f",fitRatioTrueMCInput->GetParameter(0)));
		fitParamRatioTrueMCInput->AddText(Form("Param0 Err #pm %f",fitRatioTrueMCInput->GetParError(0)));
  
     
     
		TCanvas* canvasCorrectedYieldMCInput = new TCanvas("canvasCorrectedYieldMCInput","",1350,1500);  // gives the page size
		DrawGammaCanvasSettings( canvasCorrectedYieldMCInput, 0.13, 0.02, 0.02, 0.09);
		canvasCorrectedYieldMCInput->SetLogy();

		TPad* padCorrectedYieldMCInput = new TPad("padCorrectedYieldMCInput", "", 0., 0.25, 1., 1.,-1, -1, -2);
		DrawGammaPadSettings( padCorrectedYieldMCInput, 0.12, 0.02, 0.02, 0.);
		padCorrectedYieldMCInput->Draw();

		TPad* padCorrectedYieldMCInputRatios = new TPad("padCorrectedYieldMCInputRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
		DrawGammaPadSettings( padCorrectedYieldMCInputRatios, 0.12, 0.02, 0., 0.2);
		padCorrectedYieldMCInputRatios->Draw();

		padCorrectedYieldMCInput->cd();
		padCorrectedYieldMCInput->SetLogy(); 

		DrawAutoGammaMesonHistos( histoCorrectedYieldTrue, 
                             "", "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}", 
                             kTRUE, 3., 4e-10, kTRUE,
                             kFALSE, 0., 0.7, 
                             kFALSE, 0., 10.);
		DrawGammaSetMarker(histoCorrectedYieldTrue, 22, 1., kBlack, kBlack);
		histoCorrectedYieldTrue->DrawCopy("e1");


		DrawGammaSetMarker(histoMCYieldMeson, 24, 1., kRed+2, kRed+2);                                                                           
		histoMCYieldMeson->DrawCopy("e1,same"); 
   
   
        
        
        
     
		TLegend* legendYieldMCInput = new TLegend(0.15,0.03,0.66,0.19);
		legendYieldMCInput->SetTextSize(0.02);                     
		legendYieldMCInput->SetFillColor(0);
		legendYieldMCInput->AddEntry(histoCorrectedYieldTrue,"corr true eff/right norm");
		legendYieldMCInput->AddEntry(histoMCYieldMeson,"MC input");

		legendYieldMCInput->Draw();
   
		fitParamRatioTrueMCInput->Draw("sames");

		//if (fThesis.CompareTo("thesis") == 0)DrawAliceLogoPi0WorkInProgress(pictDrawingCoordinates[0], pictDrawingCoordinates[1], pictDrawingCoordinates[2], pictDrawingCoordinates[3], pictDrawingCoordinates[4], pictDrawingCoordinates[5], pictDrawingCoordinates[6], pictDrawingCoordinates[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,1125,kDalitz);

		padCorrectedYieldMCInputRatios->cd();
		padCorrectedYieldMCInputRatios->SetTickx();
		padCorrectedYieldMCInputRatios->SetTicky();
		padCorrectedYieldMCInputRatios->SetLogy(0);
  
 
		RatioTrueMCInput02->SetYTitle("#frac{standard}{modified}");   
		RatioTrueMCInput02->GetYaxis()->SetRangeUser(0.8,1.23);
		RatioTrueMCInput02->GetYaxis()->SetLabelSize(0.07);
		RatioTrueMCInput02->GetYaxis()->SetNdivisions(505);
		RatioTrueMCInput02->GetYaxis()->SetTitleSize(0.1);    
		RatioTrueMCInput02->GetYaxis()->SetDecimals();
		RatioTrueMCInput02->GetYaxis()->SetTitleOffset(0.5);
		RatioTrueMCInput02->GetXaxis()->SetTitleSize(0.11);   
		RatioTrueMCInput02->GetXaxis()->SetLabelSize(0.08);
		RatioTrueMCInput02->SetMarkerStyle(22);
		RatioTrueMCInput02->SetMarkerSize(1.);
		RatioTrueMCInput02->SetMarkerColor(kBlack);
		RatioTrueMCInput02->SetLineColor(kBlack);
		RatioTrueMCInput02->DrawCopy("p,e1"); 

		DrawGammaSetMarker(RatioTrueMCInput02, 22, 1., kBlack, kBlack);                                                                                
		RatioTrueMCInput02->DrawCopy("p,e1");    
   

		DrawGammaLines(0., maxPtMeson,1., 1.,0.1);
		fitRatioTrueMCInput->DrawCopy("same");

		canvasCorrectedYieldMCInput->Update();
		canvasCorrectedYieldMCInput->SaveAs(Form("%s/%s_%s_CorrectedYieldMCInput_%s.%s",outputDir.Data(), nameMeson.Data(), prefix2.Data(),  fCutSelection.Data(), suffix.Data()));

		delete canvasCorrectedYieldMCInput;
		delete legendYieldMCInput;
		
		}
	  
	  
	  
		TCanvas* canvasRAWYieldSec = new TCanvas("canvasRAWYieldSec","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasRAWYieldSec, 0.13, 0.02, 0.02, 0.09); 
		canvasRAWYieldSec->SetLogy(1);         

		DrawGammaSetMarker(histoUnCorrectedYieldDrawing, 20, 1., kBlack, kBlack);                              
		histoUnCorrectedYieldDrawing->Draw("e1");
		histoYieldGGMeson->Scale(1./nEvt);
		DrawGammaSetMarker(histoYieldGGMeson, 22, 1., kBlue, kBlue);                               
		histoYieldGGMeson->DrawCopy("same,e1");   
			
		TLegend* legendSecRAWYield = new TLegend(0.6,0.8,0.97,0.95);
		legendSecRAWYield->SetTextSize(0.03);        
		legendSecRAWYield->SetFillColor(0);
		legendSecRAWYield->SetBorderSize(0);
		legendSecRAWYield->AddEntry(histoUnCorrectedYieldDrawing,"RAW yield");
		legendSecRAWYield->AddEntry(histoYieldGGMeson,"Contamination from #gamma#gamma");
		legendSecRAWYield->Draw();

		canvasRAWYieldSec->Update();
		canvasRAWYieldSec->SaveAs(Form("%s/%s_%s_RAWYieldContGGPt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
		delete canvasRAWYieldSec;
		
		
		
		TCanvas* canvasFractionContGG = new TCanvas("canvasFractionContGG","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasFractionContGG, 0.13, 0.02, 0.02, 0.09);	
		
     
		DrawAutoGammaMesonHistos( histoYieldTrueGGFracMeson, 
                             "", "p_{T} (GeV/c)", "#frac{True #pi^{0} #rightarrow #gamma#gamma}{ True #pi^{0} #rightarrow #gamma#gamma + True #pi^{0} #rightarrow e^{+}e^{-}#gamma }", 
                             kTRUE, 3., 4e-10, kTRUE,
                             kFALSE, 0., 0.30, 
                             kFALSE, 0., 10.);
		histoYieldTrueGGFracMeson->SetLineWidth(0.5);	
		histoYieldTrueGGFracMeson->GetYaxis()->SetRangeUser(0.0,0.35);
		DrawGammaSetMarker(histoYieldTrueGGFracMeson, 20, 0.5, kBlack, kBlack);										 
		histoYieldTrueGGFracMeson->DrawCopy("e1"); 	

		canvasFractionContGG->Update();

		canvasFractionContGG->SaveAs(Form("%s/%s_%s_FractionContGGPt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
		delete canvasFractionContGG;
		
		
		if( isMC.CompareTo("kFALSE") == 0 ) {
		
		
		TCanvas* canvasFractionContGGForData = new TCanvas("canvasFractionContGGForData","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasFractionContGGForData, 0.13, 0.02, 0.02, 0.09);	
		
     
		DrawAutoGammaMesonHistos( histoYieldTrueGGFracMesonForData, 
                             "", "p_{T} (GeV/c)", Form("#frac{ #frac{%f}{%f} #times True #pi^{0} #rightarrow #gamma#gamma}{ #frac{%f}{%f} #times True #pi^{0} #rightarrow #gamma#gamma + #frac{%f}{%f} #times True #pi^{0} #rightarrow e^{+}e^{-}#gamma}",fArrayBRPi0Meson->GetAt(0),fArrayBRPi0Meson->GetAt(2),
						       fArrayBRPi0Meson->GetAt(0),fArrayBRPi0Meson->GetAt(2),fArrayBRPi0Meson->GetAt(1),fArrayBRPi0Meson->GetAt(3)), 
                             kTRUE, 3., 4e-10, kTRUE,
                             kFALSE, 0., 0.30, 
                             kFALSE, 0., 10.);
		histoYieldTrueGGFracMesonForData->SetLineWidth(0.5);	
		histoYieldTrueGGFracMesonForData->GetYaxis()->SetRangeUser(0.0,0.35);
		DrawGammaSetMarker(histoYieldTrueGGFracMesonForData, 20, 0.5, kBlack, kBlack);										 
		histoYieldTrueGGFracMesonForData->DrawCopy("e1"); 	

		canvasFractionContGGForData->Update();

		canvasFractionContGGForData->SaveAs(Form("%s/%s_%s_FractionContGGPtForData_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
		delete canvasFractionContGGForData;
		
		
		}
		
		
		
		////////////////Significance/////////////////////////////////////////////////////////////
          
		TCanvas* canvasSB = new TCanvas("canvasSB","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasSB, 0.13, 0.02, 0.02, 0.09);	
		canvasSB->SetLogy(1);			

		TH1D* histoSBMesonDrawing = (TH1D*)histoSBMeson->Clone();
  
	     
		DrawAutoGammaMesonHistos(histoSBMesonDrawing,
                                   "", "p_{T} (GeV/c)", "Signal/Background",
                                   kTRUE, 5., 10e-10, kTRUE,
                                   kFALSE, 0.0, 0.030,
                                   kFALSE, 0., 10.);
		histoSBMesonDrawing->SetLineWidth(0.8);					 
		DrawGammaSetMarker(histoSBMesonDrawing, 20, 0.5, kBlue, kBlue);										 
		histoSBMesonDrawing->DrawCopy("e1"); 	

		//if (fThesis.CompareTo("thesis") == 0)DrawAliceLogoPi0Performance(pictDrawingCoordinatesInv[0], pictDrawingCoordinatesInv[1], pictDrawingCoordinatesInv[2], pictDrawingCoordinatesInv[3], pictDrawingCoordinatesInv[4], pictDrawingCoordinatesInv[5], pictDrawingCoordinatesInv[6], pictDrawingCoordinatesInv[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900,date,"MinBias",kDalitz);
		canvasSB->Update();

		canvasSB->SaveAs(Form("%s/%s_%s_SBPt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
		delete canvasSB;
      
      
      
		TCanvas* canvasSignificance = new TCanvas("canvasSignificance","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasSignificance, 0.13, 0.02, 0.02, 0.09);	
		canvasSignificance->SetLogy(1);			

		TH1D* histoSignificanceMesonDrawing = (TH1D*)histoSignificanceMeson->Clone();
      
		DrawAutoGammaMesonHistos( histoSignificanceMesonDrawing, 
                             "", "p_{T} (GeV/c)", "#frac{Signal}{#sqrt{Background}}", 
                              kTRUE, 5., 10e-10, kTRUE,
                              kFALSE, 0.0, 0.030,
                              kFALSE, 0., 10.);
				      
		histoSignificanceMesonDrawing->SetLineWidth(0.8);					 
		DrawGammaSetMarker(histoSignificanceMesonDrawing, 20, 0.5, kBlue, kBlue);										 
		histoSignificanceMesonDrawing->DrawCopy("e1"); 	

		//if (fThesis.CompareTo("thesis") == 0)DrawAliceLogoPi0Performance(pictDrawingCoordinatesInv[0], pictDrawingCoordinatesInv[1], pictDrawingCoordinatesInv[2], pictDrawingCoordinatesInv[3], pictDrawingCoordinatesInv[4], pictDrawingCoordinatesInv[5], pictDrawingCoordinatesInv[6], pictDrawingCoordinatesInv[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900,date,"MinBias",kDalitz);
		canvasSignificance->Update();
  
		canvasSignificance->SaveAs(Form("%s/%s_%s_SignificancePt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
		delete canvasSignificance;
      
      
		
		
	}

	Int_t nBinsPt =   histoCorrectedYieldTrue->GetNbinsX();
	TH1D *histoCorrectedYieldTrueMt = NULL;
	TH1D *histoCorrectedYieldTrueNarrowMt = NULL;
	TH1D *histoCorrectedYieldTrueWideMt = NULL;
	TH1D *histoCorrectedYieldTrueLeftMt = NULL;
	TH1D *histoCorrectedYieldTrueLeftNarrowMt = NULL;
	TH1D *histoCorrectedYieldTrueLeftWideMt = NULL;
	TH1D *histoRatioTrueMt = NULL;
	TH1D *histoRatioTrueWideMt = NULL;
	TH1D *histoRatioTrueNarrowMt = NULL;
	TH1D *histoRatioTrueLeftMt = NULL;
	TH1D *histoRatioTrueLeftWideMt = NULL;
	TH1D *histoRatioTrueLeftNarrowMt = NULL;
	

	if (!kDalitz){
		//***********************************************************************************************
		//***************************  correction for yield in mt bins **********************************
		//***********************************************************************************************
		

		Double_t * binsMt = NULL;
		binsMt=        new Double_t[nBinsPt+1];
		binsMt[0] =       mesonMassExpect;


		for(Int_t iPt=1;iPt<nBinsPt+1;iPt++){
			binsMt[iPt]=pow((histoCorrectedYieldTrue->GetXaxis()->GetBinUpEdge(iPt)*histoCorrectedYieldTrue->GetXaxis()->GetBinUpEdge(iPt)+mesonMassExpect*mesonMassExpect),0.5);
			cout << "recalculation pt to mt:    " << iPt<<"     pt      "<< histoCorrectedYieldTrue->GetXaxis()->GetBinUpEdge(iPt)<< "      mt    " << binsMt[iPt]<< endl;
		}

		TH1F *deltaMt =   new TH1F("deltaMt","",nBinsPt,binsMt);

		histoCorrectedYieldTrueMt = new TH1D("CorrectedYieldTrueEff_Mt","",nBinsPt,binsMt);
		histoCorrectedYieldTrueMt->SetName("CorrectedYieldTrueEff_Mt");
		histoCorrectedYieldTrueNarrowMt = new TH1D("CorrectedYieldTrueEffNarrow_Mt","",nBinsPt,binsMt);
		histoCorrectedYieldTrueNarrowMt->SetName("CorrectedYieldTrueEffNarrow_Mt");
		histoCorrectedYieldTrueWideMt = new TH1D("CorrectedYieldTrueEffWide_Mt","",nBinsPt,binsMt);
		histoCorrectedYieldTrueWideMt->SetName("CorrectedYieldTrueEffWide_Mt");

		histoCorrectedYieldTrueLeftMt = new TH1D("CorrectedYieldTrueEffLeft_Mt","",nBinsPt,binsMt);
		histoCorrectedYieldTrueLeftMt->SetName("CorrectedYieldTrueEffLeft_Mt");
		histoCorrectedYieldTrueLeftNarrowMt = new TH1D("CorrectedYieldTrueEffLeftNarrow_Mt","",nBinsPt,binsMt);
		histoCorrectedYieldTrueLeftNarrowMt->SetName("CorrectedYieldTrueEffLeftNarrow_Mt");
		histoCorrectedYieldTrueLeftWideMt = new TH1D("CorrectedYieldTrueEffLeftWide_Mt","",nBinsPt,binsMt);
		histoCorrectedYieldTrueLeftWideMt->SetName("CorrectedYieldTrueEffLeftWide_Mt");

		for(Int_t iPt=1;iPt<nBinsPt+1;iPt++){
			deltaMt->SetBinContent(iPt,binsMt[iPt]-binsMt[iPt-1]);
			deltaMt->SetBinError(iPt,0);
		}

		for(Int_t iPt=1;iPt<nBinsPt+1;iPt++){
			histoCorrectedYieldTrueMt->SetBinContent(iPt,histoCorrectedYieldTrue->GetBinContent(iPt));
			histoCorrectedYieldTrueMt->SetBinError(iPt,histoCorrectedYieldTrue->GetBinError(iPt));

			histoCorrectedYieldTrueNarrowMt->SetBinContent(iPt,histoCorrectedYieldTrueNarrow->GetBinContent(iPt));
			histoCorrectedYieldTrueNarrowMt->SetBinError(iPt,histoCorrectedYieldTrueNarrow->GetBinError(iPt));

			histoCorrectedYieldTrueWideMt->SetBinContent(iPt,histoCorrectedYieldTrueWide->GetBinContent(iPt));
			histoCorrectedYieldTrueWideMt->SetBinError(iPt,histoCorrectedYieldTrueWide->GetBinError(iPt));

			histoCorrectedYieldTrueLeftMt->SetBinContent(iPt,histoCorrectedYieldTrueLeft->GetBinContent(iPt));
			histoCorrectedYieldTrueLeftMt->SetBinError(iPt,histoCorrectedYieldTrueLeft->GetBinError(iPt));

			histoCorrectedYieldTrueLeftNarrowMt->SetBinContent(iPt,histoCorrectedYieldTrueLeftNarrow->GetBinContent(iPt));
			histoCorrectedYieldTrueLeftNarrowMt->SetBinError(iPt,histoCorrectedYieldTrueLeftNarrow->GetBinError(iPt));

			histoCorrectedYieldTrueLeftWideMt->SetBinContent(iPt,histoCorrectedYieldTrueLeftWide->GetBinContent(iPt));
			histoCorrectedYieldTrueLeftWideMt->SetBinError(iPt,histoCorrectedYieldTrueLeftWide->GetBinError(iPt));
		}


		histoRatioTrueMt = (TH1D*) histoCorrectedYieldTrueMt->Clone(); 
		histoRatioTrueMt->Divide(histoRatioTrueMt,histoCorrectedYieldTrueMt,1.,1.,"");
		histoRatioTrueWideMt = (TH1D*) histoCorrectedYieldTrueMt->Clone();   
		histoRatioTrueWideMt->Divide(histoRatioTrueWideMt,histoCorrectedYieldTrueWideMt,1.,1.,"");
		histoRatioTrueNarrowMt = (TH1D*) histoCorrectedYieldTrueMt->Clone(); 
		histoRatioTrueNarrowMt->Divide(histoRatioTrueNarrowMt,histoCorrectedYieldTrueNarrowMt,1.,1.,"");
		histoRatioTrueLeftMt = (TH1D*) histoCorrectedYieldTrueMt->Clone();   
		histoRatioTrueLeftMt->Divide(histoRatioTrueMt,histoCorrectedYieldTrueLeftMt,1.,1.,"");
		histoRatioTrueLeftWideMt = (TH1D*) histoCorrectedYieldTrueMt->Clone();  
		histoRatioTrueLeftWideMt->Divide(histoRatioTrueLeftWideMt,histoCorrectedYieldTrueLeftWideMt,1.,1.,"");
		histoRatioTrueLeftNarrowMt = (TH1D*) histoCorrectedYieldTrueMt->Clone();   
		histoRatioTrueLeftNarrowMt->Divide(histoRatioTrueLeftNarrowMt,histoCorrectedYieldTrueLeftNarrowMt,1.,1.,"");

		//*************************************************************************************************
		//*********************** Plotting Corrected Yield in mt - bins normal eff ************************
		//*************************************************************************************************

		TCanvas* canvasCorrecftedYieldMt = new TCanvas("canvasCorrecftedYieldMt","",1350,1500);  // gives the page size
		DrawGammaCanvasSettings( canvasCorrecftedYieldMt, 0.13, 0.02, 0.02, 0.09); 
		canvasCorrecftedYieldMt->SetLogy(); 

		TPad* padCorrectedYieldHistosMt = new TPad("padCorrectedYieldHistosMt", "", 0., 0.25, 1., 1.,-1, -1, -2);
		DrawGammaPadSettings( padCorrectedYieldHistosMt, 0.12, 0.02, 0.02, 0.);
		padCorrectedYieldHistosMt->Draw();

		TPad* padCorrectedYieldRatiosMt = new TPad("padCorrectedYieldRatiosMt", "", 0., 0., 1., 0.25,-1, -1, -2);
		DrawGammaPadSettings( padCorrectedYieldRatiosMt, 0.12, 0.02, 0., 0.2);
		padCorrectedYieldRatiosMt->Draw();

		padCorrectedYieldHistosMt->cd();
		padCorrectedYieldHistosMt->SetLogy();     

		
		DrawAutoGammaMesonHistos( histoCorrectedYieldTrueMt, 
									"", "m_{t} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{m_{t}dm_{t}dy} (c/GeV)^{2}", 
									kTRUE, 3., 4e-10, kTRUE,
									kFALSE, 0., 0.7, 
									kFALSE, 0., 15.);
		DrawGammaSetMarker(histoCorrectedYieldTrueMt, 22, 1., kBlack, kBlack);                              
		histoCorrectedYieldTrueMt->DrawCopy("e1,same");    
		//right Side Normalization narrow
		DrawGammaSetMarker(histoCorrectedYieldTrueNarrowMt, 26, 1., kGray+1, kGray+1);                                 
		histoCorrectedYieldTrueNarrowMt->DrawCopy("e1,same"); 
		//right Side Normalization wide
		DrawGammaSetMarker(histoCorrectedYieldTrueWideMt, 26, 1., kGray+3, kGray+3);                              
		histoCorrectedYieldTrueWideMt->DrawCopy("e1,same"); 
		//left Side Normalization
		DrawGammaSetMarker(histoCorrectedYieldTrueLeftMt, 22, 1., kBlue, kBlue);                               
		histoCorrectedYieldTrueLeftMt->DrawCopy("e1,same");   
		//left Side Normalization narrow
		DrawGammaSetMarker(histoCorrectedYieldTrueLeftNarrowMt, 26, 1., kBlue-5, kBlue-5);                                
		histoCorrectedYieldTrueLeftNarrowMt->DrawCopy("e1,same"); 
		//left Side Normalization wide
		DrawGammaSetMarker(histoCorrectedYieldTrueLeftWideMt, 26, 1., kBlue+2, kBlue+2);                             
		histoCorrectedYieldTrueLeftWideMt->DrawCopy("e1,same"); 


		TLegend* legendYield_Mt = new TLegend(0.15,0.03,0.66,0.19);
		legendYield_Mt->SetTextSize(0.02);        
		legendYield_Mt->SetFillColor(0);
		legendYield_Mt->AddEntry(histoCorrectedYieldTrueMt,"corr true /right norm");
		legendYield_Mt->AddEntry(histoCorrectedYieldTrueWideMt,"corr true wide int /right norm");
		legendYield_Mt->AddEntry(histoCorrectedYieldTrueNarrowMt,"corr true narrow int /right norm");
		legendYield_Mt->AddEntry(histoCorrectedYieldTrueLeftMt,Form("corr true  /left norm"));
		legendYield_Mt->AddEntry(histoCorrectedYieldTrueLeftWideMt,"corr true wide int /left norm");
		legendYield_Mt->AddEntry(histoCorrectedYieldTrueLeftNarrowMt,"corr true narrow int /left norm");

		legendYield_Mt->Draw();

		padCorrectedYieldRatiosMt->cd();
		histoRatioTrueMt->SetYTitle("#frac{standard}{modified}");   
		histoRatioTrueMt->SetXTitle("m_{t} (GeV/c)");
		if(nameMeson.CompareTo("Pi0") == 0  || nameMeson.CompareTo("Pi0EtaBinning") == 0 ) histoRatioTrueMt->GetYaxis()->SetRangeUser(0.8,1.23);  
		if(nameMeson.CompareTo("Eta") == 0 ) histoRatioTrueMt->GetYaxis()->SetRangeUser(0.4,1.63);
		histoRatioTrueMt->GetXaxis()->SetRangeUser(0.,15.);
		histoRatioTrueMt->GetYaxis()->SetNdivisions(505);
		//histoRatioTrueMt->GetXaxis()->SetRangeUser(0.,1.5);
		histoRatioTrueMt->GetYaxis()->SetLabelSize(0.08);
		histoRatioTrueMt->GetYaxis()->SetTitleSize(0.1);   
		histoRatioTrueMt->GetYaxis()->SetDecimals();
		histoRatioTrueMt->GetYaxis()->SetTitleOffset(0.42);
		histoRatioTrueMt->GetXaxis()->SetTitleSize(0.11);  
		histoRatioTrueMt->GetXaxis()->SetLabelSize(0.08);

		DrawGammaSetMarker(histoRatioTrueMt, 22, 1., kBlack, kBlack);                              
		histoRatioTrueMt->DrawCopy("e1");   
		//right Side Normalization narrow
		DrawGammaSetMarker(histoRatioTrueNarrowMt, 26, 1., kGray+1, kGray+1);                                 
		histoRatioTrueNarrowMt->DrawCopy("e1,same"); 
		//right Side Normalization wide
		DrawGammaSetMarker(histoRatioTrueWideMt, 26, 1., kGray+3, kGray+3);                              
		histoRatioTrueWideMt->DrawCopy("e1,same"); 
		//left Side Normalization
		DrawGammaSetMarker(histoRatioTrueLeftMt, 22, 1., kBlue, kBlue);                               
		histoRatioTrueLeftMt->DrawCopy("e1,same");   
		//left Side Normalization narrow
		DrawGammaSetMarker(histoRatioTrueLeftNarrowMt, 26, 1., kBlue-5, kBlue-5);                                
		histoRatioTrueLeftNarrowMt->DrawCopy("e1,same"); 
		//left Side Normalization wide
		DrawGammaSetMarker(histoRatioTrueLeftWideMt, 26, 1., kBlue+2, kBlue+2);                             
		histoRatioTrueLeftWideMt->DrawCopy("e1,same"); 
		DrawGammaLines(0., maxPtMeson,1., 1.,0.1);

		canvasCorrecftedYieldMt->Update();

		canvasCorrecftedYieldMt->SaveAs(Form("%s/%s_%s_CorrectedYield_Mtspectra_TrueEff_%s.%s", outputDir.Data(), nameMeson.Data() ,prefix2.Data() , fCutSelection.Data(), suffix.Data()));
		delete canvasCorrecftedYieldMt;
		delete legendYield_Mt;
			
		//***********************************************************************************************
		//***************************  correction for yield in Xt bins **********************************
		//***********************************************************************************************
		Double_t * binsXt = NULL;
		binsXt=  new Double_t[nBinsPt+1];
			
		for(Int_t iPt=0;iPt<nBinsPt+1;iPt++){
			binsXt[iPt]=2*histoCorrectedYieldTrue->GetXaxis()->GetBinUpEdge(iPt)/energy;
	//          cout << "recalculation pt to xt:    " << iPt<<"     pt      "<< histoCorrectedYieldTrue->GetXaxis()->GetBinUpEdge(iPt)<< "      xt    " << binsXt[iPt]<< endl;
		}
	}
	//*************************************************************************************************
	//******************** Output of the systematic Error due to Signal extraction ********************
	//*************************************************************************************************
	Double_t  binsXCenter[50];
	Double_t  binsXWidth[50];
	//    Double_t binYValue[50];
	//    binsXCenter[0] =     0;
	binsXWidth[0]=       0.;
	//    binYValue[0]=        0.;
	for (Int_t i = 1; i < nBinsPt +1; i++){
		binsXCenter[i] =  histoCorrectedYieldTrue->GetBinCenter(i);
		binsXWidth[i]=    histoCorrectedYieldTrue->GetBinWidth(i)/2.;
	}


	SysErrorConversion sysErrNormal[50];
	SysErrorConversion sysErrLeft[50];
	SysErrorConversion sysErrWide[50];
	SysErrorConversion sysErrNarrow[50];
	SysErrorConversion sysErrLeftWide[50];
	SysErrorConversion sysErrLeftNarrow[50];

	for (Int_t i = 1; i < nBinsPt +1; i++){
		if (optionEnergy.CompareTo("PbPb_2.76TeV")==0  ){  //|| optionEnergy.CompareTo("2.76TeV")==0
	//          binYValue[i] = histoCorrectedYieldTrueFitted->GetBinContent(i);
			sysErrNormal[i].value = histoCorrectedYieldTrueFitted->GetBinContent(i);
			sysErrNormal[i].error = histoCorrectedYieldTrueFitted->GetBinError(i);
			sysErrLeft[i].value = histoCorrectedYieldTrueLeftFitted->GetBinContent(i);
			sysErrLeft[i].error = histoCorrectedYieldTrueLeftFitted->GetBinError(i);
			sysErrNarrow[i].value = histoCorrectedYieldTrueNarrowFitted->GetBinContent(i);
			sysErrNarrow[i].error = histoCorrectedYieldTrueNarrowFitted->GetBinError(i);
			sysErrLeftNarrow[i].value = histoCorrectedYieldTrueLeftNarrowFitted->GetBinContent(i);
			sysErrLeftNarrow[i].error = histoCorrectedYieldTrueLeftNarrowFitted->GetBinError(i);
			sysErrWide[i].value = histoCorrectedYieldTrueWideFitted->GetBinContent(i);
			sysErrWide[i].error = histoCorrectedYieldTrueWideFitted->GetBinError(i);
			sysErrLeftWide[i].value = histoCorrectedYieldTrueLeftWideFitted->GetBinContent(i);
			sysErrLeftWide[i].error = histoCorrectedYieldTrueLeftWideFitted->GetBinError(i); 
		} else {
	//          binYValue[i] = histoCorrectedYieldTrue->GetBinContent(i);
			sysErrNormal[i].value = histoCorrectedYieldTrue->GetBinContent(i);
			sysErrNormal[i].error = histoCorrectedYieldTrue->GetBinError(i);
			sysErrLeft[i].value = histoCorrectedYieldTrueLeft->GetBinContent(i);
			sysErrLeft[i].error = histoCorrectedYieldTrueLeft->GetBinError(i);
			sysErrNarrow[i].value = histoCorrectedYieldTrueNarrow->GetBinContent(i);
			sysErrNarrow[i].error = histoCorrectedYieldTrueNarrow->GetBinError(i);
			sysErrLeftNarrow[i].value = histoCorrectedYieldTrueLeftNarrow->GetBinContent(i);
			sysErrLeftNarrow[i].error = histoCorrectedYieldTrueLeftNarrow->GetBinError(i);
			sysErrWide[i].value = histoCorrectedYieldTrueWide->GetBinContent(i);
			sysErrWide[i].error = histoCorrectedYieldTrueWide->GetBinError(i);
			sysErrLeftWide[i].value = histoCorrectedYieldTrueLeftWide->GetBinContent(i);
			sysErrLeftWide[i].error = histoCorrectedYieldTrueLeftWide->GetBinError(i); 
		}
	}
	Double_t differenceLeft[50];
	Double_t differenceLeftWide[50];
	Double_t differenceLeftNarrow[50];
	Double_t differenceWide[50];
	Double_t differenceNarrow[50];
	Double_t largestDifferenceNeg[50];
	Double_t largestDifferencePos[50];
	Double_t differenceLeftError[50];
	Double_t differenceLeftWideError[50];
	Double_t differenceLeftNarrowError[50];
	Double_t differenceWideError[50];
	Double_t differenceNarrowError[50];
	Double_t largestDifferenceNegError[50];
	Double_t largestDifferencePosError[50];
	Double_t relDifferenceLeft[50];
	Double_t relDifferenceLeftWide[50];
	Double_t relDifferenceLeftNarrow[50];
	Double_t relDifferenceWide[50];
	Double_t relDifferenceNarrow[50];
	Double_t relLargestDifferenceNeg[50];
	Double_t relLargestDifferencePos[50];
	Double_t relDifferenceLeftError[50];
	Double_t relDifferenceLeftWideError[50];
	Double_t relDifferenceLeftNarrowError[50];
	Double_t relDifferenceWideError[50];
	Double_t relDifferenceNarrowError[50];
	Double_t relLargestDifferenceNegError[50];
	Double_t relLargestDifferencePosError[50];
	for (Int_t i = 1; i < nBinsPt +1; i++){
		largestDifferenceNeg[i] =     0;
		largestDifferencePos[i] =     0;
		largestDifferenceNegError[i] =   0;
		largestDifferencePosError[i] =   0;
		relLargestDifferenceNeg[i] =     0;
		relLargestDifferencePos[i] =     0;
		relLargestDifferenceNegError[i] =   0;
		relLargestDifferencePosError[i] =   0;
		//Calculate differences
		differenceLeft[i] =        sysErrLeft[i].value - sysErrNormal[i].value;
		differenceLeftError[i] =      TMath::Sqrt(TMath::Abs(TMath::Power(sysErrLeft[i].error,2)-TMath::Power(sysErrNormal[i].error,2)));
		differenceLeftNarrow[i] =     sysErrLeftNarrow[i].value - sysErrNormal[i].value;
		differenceLeftNarrowError[i] =   TMath::Sqrt(TMath::Abs(TMath::Power(sysErrLeftNarrow[i].error,2)-TMath::Power(sysErrNormal[i].error,2)));
		differenceLeftWide[i] =          sysErrLeftWide[i].value - sysErrNormal[i].value;
		differenceLeftWideError[i] =     TMath::Sqrt(TMath::Abs(TMath::Power(sysErrLeftWide[i].error,2)-TMath::Power(sysErrNormal[i].error,2)));
		differenceNarrow[i] =         sysErrNarrow[i].value - sysErrNormal[i].value;
		differenceNarrowError[i] =       TMath::Sqrt(TMath::Abs(TMath::Power(sysErrNarrow[i].error,2)-TMath::Power(sysErrNormal[i].error,2)));
		differenceWide[i] =        sysErrWide[i].value - sysErrNormal[i].value;
		differenceWideError[i] =      TMath::Sqrt(TMath::Abs(TMath::Power(sysErrWide[i].error,2)-TMath::Power(sysErrNormal[i].error,2)));

		if (sysErrNormal[i].value != 0){
			relDifferenceLeft[i] =     (differenceLeft[i]/sysErrNormal[i].value)*100.;
			relDifferenceLeftError[i] =   (differenceLeftError[i]/sysErrNormal[i].value)*100.;
			relDifferenceLeftNarrow[i] =  (differenceLeftNarrow[i]/sysErrNormal[i].value)*100.;
			relDifferenceLeftNarrowError[i] = (differenceLeftNarrowError[i]/sysErrNormal[i].value)*100.;
			relDifferenceLeftWide[i] =    (differenceLeftWide[i]/sysErrNormal[i].value)*100.;
			relDifferenceLeftWideError[i] = (differenceLeftWideError[i]/sysErrNormal[i].value)*100.;
			relDifferenceNarrow[i] =   (differenceNarrow[i]/sysErrNormal[i].value)*100.;
			relDifferenceNarrowError[i] = (differenceNarrowError[i]/sysErrNormal[i].value)*100.;
			relDifferenceWide[i] =     (differenceWide[i]/sysErrNormal[i].value)*100.;
			relDifferenceWideError[i] =   (differenceWideError[i]/sysErrNormal[i].value)*100.;
		} else {
			relDifferenceLeft[i] =     0.;
			relDifferenceLeftError[i] =   0.;
			relDifferenceLeftNarrow[i] =  0.;
			relDifferenceLeftNarrowError[i] = 0.;
			relDifferenceLeftWide[i] =    0.;
			relDifferenceLeftWideError[i] = 0.;
			relDifferenceNarrow[i] =      0.;
			relDifferenceNarrowError[i] = 0.;
			relDifferenceWide[i] =     0.;
			relDifferenceWideError[i] =   0.;
		}

		//Find biggest Deviation
		if (TMath::Abs(relDifferenceLeft[i]) < 75. ){
			if(differenceLeft[i] < 0){
				largestDifferenceNeg[i] =     differenceLeft[i];
				largestDifferenceNegError[i] =   differenceLeftError[i];
				relLargestDifferenceNeg[i] =     relDifferenceLeft[i];
				relLargestDifferenceNegError[i] =   relDifferenceLeftError[i];
			}else{   
				largestDifferencePos[i] =     differenceLeft[i]; 
				largestDifferencePosError[i] =   differenceLeftError[i]; 
				relLargestDifferencePos[i] =     relDifferenceLeft[i]; 
				relLargestDifferencePosError[i] =   relDifferenceLeftError[i]; 
			}  
		}
		if (TMath::Abs(relDifferenceNarrow[i]) < 75.){
			if(differenceNarrow[i] < 0){
				if(differenceNarrow[i] < largestDifferenceNeg[i]){
				largestDifferenceNeg[i] =     differenceNarrow[i]; 
				largestDifferenceNegError[i] =   differenceNarrowError[i];  
				relLargestDifferenceNeg[i] =     relDifferenceNarrow[i]; 
				relLargestDifferenceNegError[i] =   relDifferenceNarrowError[i];  
				}
			} else {
				if(differenceNarrow[i] > largestDifferencePos[i]){
				largestDifferencePos[i] =     differenceNarrow[i];
				largestDifferencePosError[i] =   differenceNarrowError[i];
				relLargestDifferencePos[i] =     relDifferenceNarrow[i];
				relLargestDifferencePosError[i] =   relDifferenceNarrowError[i];
				}
			}  
		}
		if (TMath::Abs(relDifferenceWide[i]) < 75.){
			if(differenceWide[i] < 0){
				if(differenceWide[i] < largestDifferenceNeg[i]){
				largestDifferenceNeg[i] =     differenceWide[i];
				largestDifferenceNegError[i] =   differenceWideError[i];
				relLargestDifferenceNeg[i] =     relDifferenceWide[i];
				relLargestDifferenceNegError[i] =   relDifferenceWideError[i];
				}
			} else {
				if(differenceWide[i] > largestDifferencePos[i]) {
				largestDifferencePos[i] =     differenceWide[i];   
				largestDifferencePosError[i] =   differenceWideError[i]; 
				relLargestDifferencePos[i] =     relDifferenceWide[i];   
				relLargestDifferencePosError[i] =   relDifferenceWideError[i]; 
				}
			}
		}
	}  

	cout << "done filling" << endl;
	const char *nameFileSysErrDat = Form("%s/%s/%s_%s_SystematicErrorYieldExtraction_%s.dat",fCutSelection.Data(),optionEnergy.Data() ,nameMeson.Data(),prefix2.Data(),fCutSelection.Data());
	fstream fileSysErrDat;
	fileSysErrDat.open(nameFileSysErrDat, ios::out);
	fileSysErrDat << "Calculation of the systematic error due to the yield extraction" << endl;
	fileSysErrDat << "Bin" << "\t" << "Normal value" << "\t" << "Normal error" << endl;
	for(Int_t i = 1; i < (nBinsPt +1); i++){
		fileSysErrDat << i << "\t" << sysErrNormal[i].value << "\t" << sysErrNormal[i].error << endl;
	}
	fileSysErrDat << endl;
	fileSysErrDat << "Bin" << "\t" << "Narrow value" << "\t" << "Narrow error" << "\t" << "Diff to Norm" << endl;
	for(Int_t i = 1; i < (nBinsPt +1); i++){
		fileSysErrDat << i << "\t" << sysErrNarrow[i].value << "\t" << sysErrNarrow[i].error << "\t" << differenceNarrow[i] << "\t" << differenceNarrowError[i] << "\t" << relDifferenceNarrow[i] << "\t" << relDifferenceNarrowError[i] <<endl;
	}
	fileSysErrDat << endl;
	fileSysErrDat << "Bin" << "\t" << "Wide value" << "\t" << "Wide error" << "\t" << "Diff to Norm" << endl;
	for(Int_t i = 1; i < (nBinsPt +1); i++){
		fileSysErrDat << i << "\t" << sysErrWide[i].value << "\t" << sysErrWide[i].error << "\t" << differenceWide[i] << "\t" << differenceWideError[i] << "\t" << relDifferenceWide[i] << "\t" << relDifferenceWideError[i] <<endl;
	}
	fileSysErrDat << endl;
	fileSysErrDat << "Bin" << "\t" << "Left value" << "\t" << "Left error" << "\t" << "Diff to Norm" << endl;
	for(Int_t i = 1; i < (nBinsPt +1); i++){
		fileSysErrDat << i << "\t" << sysErrLeft[i].value << "\t" << sysErrLeft[i].error << "\t" << differenceLeft[i] << "\t" << differenceLeftError[i] << "\t" << relDifferenceLeft[i] << "\t" << relDifferenceLeftError[i] <<endl;
	}
	fileSysErrDat << endl;
	fileSysErrDat << "Bin" << "\t" << "Left Narrow value" << "\t" << "Left Narrow error" << "\t" << "Diff to Norm" << endl;
	for(Int_t i = 1; i < (nBinsPt +1); i++){
		fileSysErrDat << i << "\t" << sysErrLeftNarrow[i].value << "\t" << sysErrLeftNarrow[i].error << "\t" << differenceLeftNarrow[i] <<  "\t" << differenceLeftNarrowError[i] << "\t" << relDifferenceLeftNarrow[i] << "\t" << relDifferenceLeftNarrowError[i] <<endl;
	}
	fileSysErrDat << endl;
	fileSysErrDat << "Bin" << "\t" << "Left Wide value" << "\t" << "Left Wide error" << "\t" << "Diff to Norm" << endl;
	for(Int_t i = 1; i < (nBinsPt +1); i++){
		fileSysErrDat << i << "\t" << sysErrLeftWide[i].value << "\t" << sysErrLeftWide[i].error << "\t" << differenceLeftWide[i] << "\t" << differenceLeftWideError[i] << "\t" << relDifferenceLeftWide[i] << "\t" << relDifferenceLeftWideError[i] <<endl;
	}
	fileSysErrDat << endl;

	fileSysErrDat << "Bin" << "\t" << "Largest Dev Neg" << "\t" << "Largest Dev Pos"  << endl;
	for(Int_t i = 1; i < (nBinsPt +1); i++){
		fileSysErrDat << i << "\t" << largestDifferenceNeg[i] <<  "\t" << largestDifferenceNegError[i] << "\t" << largestDifferencePos[i] << "\t" << largestDifferencePosError[i]<<endl;
	}
	fileSysErrDat << "Bin" << "\t" << "Largest Rel Dev Neg" << "\t" << "Largest Rel Dev Pos"  << endl;
	for(Int_t i = 1; i < (nBinsPt +1); i++){
		fileSysErrDat << i << "\t" << relLargestDifferenceNeg[i] <<  "\t" << relLargestDifferenceNegError[i] << "\t" << relLargestDifferencePos[i] << "\t" << relLargestDifferencePosError[i]<<endl;
	}

	fileSysErrDat.close();
	
	TGraphAsymmErrors* SystErrGraphNeg = new TGraphAsymmErrors(nBinsPt+1, binsXCenter, relLargestDifferenceNeg, binsXWidth, binsXWidth, relLargestDifferenceNegError, relLargestDifferenceNegError);
	TGraphAsymmErrors* SystErrGraphPos = new TGraphAsymmErrors(nBinsPt+1, binsXCenter, relLargestDifferencePos, binsXWidth, binsXWidth, relLargestDifferencePosError, relLargestDifferencePosError);

	Double_t relBGEstimate[50];
	Double_t relBGEstimateError[50];
	if( !kDalitz ){
	for (Int_t i = 1; i < nBinsPt +1; i++){
			relBGEstimateError[i] = 0.;
			if ( TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateCatC->GetBinContent(i)) > TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateCatD->GetBinContent(i)) && TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateA->GetBinContent(i)) > TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateA->GetBinContent(i))){
				relBGEstimate[i] = TMath::Abs(histoBGEstimateCatA->GetBinContent(i)- histoBGEstimateCatC->GetBinContent(i)) *100;
			} else if ( TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateCatD->GetBinContent(i)) > TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateCatC->GetBinContent(i)) && TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateA->GetBinContent(i)) > TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateA->GetBinContent(i))) {
			relBGEstimate[i] = TMath::Abs(histoBGEstimateCatA->GetBinContent(i)- histoBGEstimateCatC->GetBinContent(i)) *100;
			} else {
			relBGEstimate[i] = TMath::Abs(histoBGEstimateCatA->GetBinContent(i)- histoBGEstimateA->GetBinContent(i)) *100;
			}
	 }   
	}
	
	TGraphAsymmErrors* SystErrGraphBGEstimate = NULL;
	
	SystErrGraphBGEstimate = new TGraphAsymmErrors(nBinsPt+1, binsXCenter, relBGEstimate, binsXWidth, binsXWidth, relBGEstimateError, relBGEstimateError);
		
	const char* nameOutput2 = Form("%s/%s/%s_%s_GammaConv_OnlyCorrectionFactor%s_%s.root",fCutSelection.Data(),optionEnergy.Data(),nameMeson.Data(),prefix2.Data(),optionPeriod.Data(),fCutSelection.Data());
	TFile* correctedOutput2 = new TFile(nameOutput2,"RECREATE");     
		histoAcceptance->Write();
		histoTrueEffiPt->Write("TrueMesonEffiPt");
		histoCompleteCorr->Write("EffiTimesAcceptanceTimesDeltaY");
	correctedOutput2->Write();
	correctedOutput2->Close();
	
	
	const char* nameOutput = Form("%s/%s/%s_%s_GammaConvV1%sCorrection%s_%s.root",fCutSelection.Data(),optionEnergy.Data(),nameMeson.Data(),prefix2.Data(),fDalitz.Data(),optionPeriod.Data(),fCutSelection.Data());
	TFile* correctedOutput = new TFile(nameOutput,"RECREATE");     

	histoCorrectedYieldNorm->Write();
	SystErrGraphPos->Write(Form("%s_SystErrorRelPos_YieldExtraction_%s",nameMeson.Data(),centralityString.Data()),TObject::kOverwrite);
	SystErrGraphNeg->Write(Form("%s_SystErrorRelNeg_YieldExtraction_%s",nameMeson.Data(),centralityString.Data()),TObject::kOverwrite);
	if (SystErrGraphBGEstimate) SystErrGraphBGEstimate->Write(Form("%s_SystErrorRel_BGEstimate_%s",nameMeson.Data(),centralityString.Data()),TObject::kOverwrite);
	if (histoBGEstimateCatA)histoBGEstimateCatA->Write("BGEstimateFromPileup");
	if(histoCorrectionFactorsHistvsPtCatA)histoCorrectionFactorsHistvsPtCatA->Write("PileupContamination");
	histoCorrectedYieldTrue->Write();
	histoCorrectedYieldTrueNarrow->Write();
	histoCorrectedYieldTrueWide->Write();
	histoCorrectedYieldFixed->Write();
	histoCorrectedYieldNarrowFixed->Write();
	histoCorrectedYieldWideFixed->Write();
	histoCorrectedYieldTrueFixed->Write();
	histoCorrectedYieldTrueNarrowFixed->Write();
	histoCorrectedYieldTrueWideFixed->Write();
	histoCorrectedYieldTrueLeft->Write();
	histoCorrectedYieldTrueLeftWide->Write();
	histoCorrectedYieldTrueLeftNarrow->Write();
	histoCorrectedYieldTrueFitted->Write();
	histoCorrectedYieldTrueNarrowFitted->Write();
	histoCorrectedYieldTrueWideFitted->Write();
	histoCorrectedYieldTrueLeftFitted->Write();
	histoCorrectedYieldTrueLeftWideFitted->Write();
	histoCorrectedYieldTrueLeftNarrowFitted->Write();
	
	
	
	/*********************New**************************/
   
	  histoYieldTrueGGFracMeson->Write();
	  histoYieldTrueGGFracMesonWide->Write();
	  histoYieldTrueGGFracMesonNarrow->Write();
	  
	  histoYieldTrueGGFracMesonForData->Write();
	  histoYieldTrueGGFracMesonWideForData->Write();
	  histoYieldTrueGGFracMesonNarrowForData->Write();
	

	  histoYieldGGMeson->Write();
	  histoYieldGGMesonLeft->Write();
   
	  histoYieldGGMesonNarrow->Write();
	  histoYieldGGMesonLeftNarrow->Write();
   
	  histoYieldGGMesonWide->Write();   
	  histoYieldGGMesonLeftWide->Write();
      
      /*************************************************/  
	
	
	
	
	
	
	if (!kDalitz){
		histoCorrectedYieldTrueMt->Write();
		histoCorrectedYieldTrueNarrowMt->Write();
		histoCorrectedYieldTrueWideMt->Write();
		histoCorrectedYieldTrueLeftMt->Write();
		histoCorrectedYieldTrueLeftWideMt->Write();
		histoCorrectedYieldTrueLeftNarrowMt->Write();

	}
	
	histoUnCorrectedYield->Write();

	histoFWHMMeson->Write();
	histoMassMeson->Write();
	
	histoAcceptance->Write();
	histoTrueEffiPt->Write("TrueMesonEffiPt");
	//histoTrueNormEffiRatio->Write();
	//    histoCompleteCorr->Write("EffiTimesAcceptanceTimesDeltaY");
	histoTrueEffiPtFit->Write("TrueMesonEffiPtFitted");
	histoRatioTrueEffiDivFitted->Write("TrueMesonEffiPtDivFittedEffi");
	histoFWHMMeson->Write();
	histoTrueFWHMMeson->Write();
	histoTrueMassMeson->Write();
	histoMesonSignalFullPtInvMass->SetName("FullInvariantMass");
	histoMesonSignalFullPtInvMass->Write();
	histoEventQuality->Write();
	histoNumberOfGoodESDTracksVtx->Write();
	histoInputMesonPt->Write();
	histoMCYieldMeson->Write();
	cout<<"save"<<endl;
	histoMCYieldMesonOldBin->Write();
	if (histoInputMesonOldBinPtWOWeights) histoInputMesonOldBinPtWOWeights->Write();
	if (histoInputMesonOldBinPtWeights) histoInputMesonOldBinPtWeights->Write("WeightsMeson");
	if (histoMCInputAddedSig)   histoMCInputAddedSig->Write();
	if (histoMCInputWOWeightingAddedSig)   histoMCInputWOWeightingAddedSig->Write();
	if (histoMCInputWeightsAddedSig) histoMCInputWeightsAddedSig->Write("WeightsMeson_AddedSig");
	
	histoUnCorrectedYieldDrawing->SetName("histoYieldMesonPerEvent");
	histoUnCorrectedYieldDrawing->Write();
	deltaPt->Write("deltaPt");
	correctedOutput->Write();
	correctedOutput->Close();
	
	if (histoEffiNarrowPt || histoEffiWidePt || histoEffiLeftPt || histoEffiLeftNarrowPt || histoEffiLeftWidePt){}
}
