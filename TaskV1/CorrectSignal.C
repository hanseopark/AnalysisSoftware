/****************************************************************************************************************************
 ****** 	provided by Gamma Conversion Group, PWG4,								*****
 ******		Ana Marin, marin@physi.uni-heidelberg.de								*****
 ******	   	Kathrin Koch, kkoch@physi.uni-heidelberg.de 								*****
 ******		Friederike Bock, friederike.bock@cern.ch								*****
 *****************************************************************************************************************************/

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
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctions.h"
// #include "CalculateGammaToPi0.h"
//#include "ExtractSignal.h" 

struct SysErrorConversion {
   Double_t value;
   Double_t error;
   //	TString name;
};

void CorrectYieldDalitz(TH1D* histoCorrectedYield,TH1D* histoRawGGYield, TH1D* histoEffiPt, TH1D* histoAcceptance, Double_t deltaRapid, Double_t scaling, Double_t nEvt, TString nameMeson){
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
   if (nameMeson.CompareTo("Pi0") == 0 ||nameMeson.CompareTo("Pi0EtaBinning") == 0 ){
      histoCorrectedYield->Scale(1./0.01198);
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
      histoCorrectedYield->Scale(1./0.98798);
   }else{
      histoCorrectedYield->Scale(1./0.3931);
   }
}

void ScaleMCYield(TH1D* histoCorrectedToBeScaled, Double_t deltaRapid, Double_t scaling, Double_t nEvtMC, TString nameMeson, TString optionDalitz ){
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
      if (optionDalitz.CompareTo("kFALSE")==0){
         histoCorrectedToBeScaled->Scale(1./0.98798);
      } else {
         histoCorrectedToBeScaled->Scale(1./0.01198);
      }
   }else{
      if (optionDalitz.CompareTo("kFALSE")==0){
         histoCorrectedToBeScaled->Scale(1./0.3931);
      } else {
         histoCorrectedToBeScaled->Scale(1./6.8e-5);
      }
		
   }
}


void  CorrectSignal(const char *fileNameUnCorrectedFile = "myOutput", const char *fileNameCorrectionFile = "", TString cutSelection ="", const char *suffix = "gif", TString nameMeson ="", const TString isMC ="", TString optionEnergy="", TString multFlag="", TString thesis="", TString optDalitz = "kFALSE"){	
   
   gROOT->Reset();	
   gROOT->SetStyle("Plain");
	
   StyleSettingsThesis();	
	
   SetPlotStyle();
		
   TString fThesis = thesis;
	
   TString fDalitz="";
   Bool_t kDalitz = kFALSE;
   if (optDalitz.CompareTo("kTRUE")==0){
      fDalitz = "Dalitz";
      kDalitz = kTRUE;
   } else {
      fDalitz = "";
      kDalitz = kFALSE;
   }
	
   TString outputDir = Form("%s/%s/%s/CorrectSignal%s",cutSelection.Data(),optionEnergy.Data(),suffix, fDalitz.Data());
   gSystem->Exec("mkdir "+outputDir);
	
   TString collisionSystem= ReturnFullCollisionsSystem(optionEnergy);
   Double_t energy = ReturnCollisionEnergy( optionEnergy);
   Double_t doubleAddFactorK0s = ReturnCorrectK0ScalingFactor( optionEnergy,  cutSelection);
	
   if (collisionSystem.CompareTo("") == 0){
      cout << "No correct collision system specification, has been given" << endl;
      return;		
   }
  
   TString date = ReturnDateString();
   
   //declaration for printing logo 
   Float_t 	pictDrawingCoordinates[9] = 		{0.55, 0.8, 0.25, 0.04, 0.7,0.5, 0.18, 0.035,0};
   Float_t 	pictDrawingCoordinatesEff[9] = 	{0.55, 0.8, 0.18, 0.04, 0.3, 0.2, 0.18, 0.035,0};
   Float_t 	pictDrawingCoordinatesAcc[9] = 	{0.63, 0.2, 0.40, 0.04, 0.7,0.33, 0.18, 0.035,0};
   Float_t 	pictDrawingCoordinatesInv[9] = 	{0.63, 0.8, 0.40, 0.04, 0.7, 0.43, 0.18, 0.035,0}; // kk
   Float_t 	pictDrawingCoordinatesMass[9] = 	{0.4, 0.8, 0.75, 0.04, 0.2,0.7, 0.18, 0.035,0};
   Float_t 	pictDrawingCoordinatesFWHM[9] = 	{0.4, 0.8, 0.75, 0.04, 0.2,0.68, 0.18, 0.035,0};
   Bool_t 	pictDrawingOptions[4] = 			{kFALSE, kFALSE, kFALSE, kTRUE};
   TString 	prefix2="";
   TString 	finder( cutSelection(1,1));
	
   if (isMC.CompareTo("kTRUE") ==0){ 
      prefix2 = 			"MC";
      doubleAddFactorK0s = 0.;
      pictDrawingOptions[1] = 	kTRUE;
   } else {	
      prefix2 = 			"data";
      pictDrawingOptions[1] = 	kFALSE;
   }
   cout << "The additional K0 correction factor is: "  << doubleAddFactorK0s<<endl;
   TString textMeson=ReturnMesonString ( nameMeson);
   pictDrawingOptions[3] = ReturnMesonOption (nameMeson);
   
   TString fGammaCutSelection  ="";
   TString fElectronCutSelection="";
   TString fMesonCutSelection="";
	
   if (kDalitz){
      ReturnSeparatedCutNumber(cutSelection, fGammaCutSelection, fElectronCutSelection,fMesonCutSelection,kTRUE);
   } else {
      ReturnSeparatedCutNumber(cutSelection, fGammaCutSelection, fElectronCutSelection,fMesonCutSelection);
   }
	
   TString rapidityRange = "";
   Double_t deltaRapid =  ReturnRapidityStringAndDouble(fMesonCutSelection, rapidityRange);
	   
   //Variable defintion
   Double_t scaling = 1./(2.*TMath::Pi());
	
   // File definitions
   TString fileNameUnCorrectedFileChar = 	fileNameUnCorrectedFile;
   TFile fileUncorrected(fileNameUnCorrectedFileChar.Data());  
   TH1F *histoNumberOfGoodESDTracksVtx = 	(TH1F*)fileUncorrected.Get("GoodESDTracks");
   TH1D *histoEventQuality = 			(TH1D*)fileUncorrected.Get("NEvents");
   TH1D *histoUnCorrectedYield =			(TH1D*)fileUncorrected.Get("histoYieldMeson");
//    TH1D *histoUnCorrectedYieldBackFit =			(TH1D*)fileUncorrected.Get("histoYieldMesonBackFit");
   TH1D *histoUnCorrectedYieldWide =		(TH1D*)fileUncorrected.Get("histoYieldMesonWide");
   TH1D *histoUnCorrectedYieldNarrow =	(TH1D*)fileUncorrected.Get("histoYieldMesonNarrow");
   TH1D *histoUnCorrectedYieldLeft = 		(TH1D*)fileUncorrected.Get("histoYieldMesonLeft");
   TH1D *histoUnCorrectedYieldLeftWide = 	(TH1D*)fileUncorrected.Get("histoYieldMesonLeftWide");
   TH1D *histoUnCorrectedYieldLeftNarrow = (TH1D*)fileUncorrected.Get("histoYieldMesonLeftNarrow");
   TH1D *histoFWHMMeson = 				(TH1D*)fileUncorrected.Get("histoFWHMMeson");
   TH1D *histoFWHMMesonLeft = 			(TH1D*)fileUncorrected.Get("histoFWHMMesonLeft");
	
   TH1D *histoMassMeson = 				(TH1D*)fileUncorrected.Get("histoMassMeson");
   TH1D *histoMassMesonLeft = 			(TH1D*)fileUncorrected.Get("histoMassMesonLeft");
	
   TH1D *histoMesonSignalFullPtInvMass=	(TH1D*)fileUncorrected.Get("Mapping_GG_InvMass_FullPt");
   TH1D *histoMesonBckNormFullPtInvMass=	(TH1D*)fileUncorrected.Get("Mapping_BackNorm_InvMass_FullPt");

    Float_t nEvt = 0;
   if (optionEnergy.CompareTo("PbPb_2.76TeV") == 0){
      nEvt = histoEventQuality->GetBinContent(1);
   } else {
      nEvt = GetNEvents(histoEventQuality);
      // BinContent 5 - Zvertex-position, BinContent 4 - no Trigger Bit, BinContent 7 - PileUp 
   }
  
   pictDrawingCoordinates[8] = 			nEvt;
   TH1D *deltaPt =					(TH1D*)fileUncorrected.Get("deltaPt");
	
   for (Int_t i = 0; i < deltaPt->GetNbinsX() +1; i++){
      deltaPt->SetBinError(i, 0);
   }	
   Double_t maxPtMeson= histoUnCorrectedYield->GetXaxis()->GetBinUpEdge(histoUnCorrectedYield->GetNbinsX());
	
   TFile* fileCorrectionsSecPi0 = new TFile("ExternalInput/PCM/SecondaryFractionHistogramms7TeV.root");
	
   TH1D*	histoDefaultTrueSecFracMeson = 				(TH1D*)fileCorrectionsSecPi0->Get("TrueSecFrac");
   TH1D*	histoDefaultTrueSecFracMesonWide = 			(TH1D*)fileCorrectionsSecPi0->Get("TrueSecFracWide");
   TH1D*	histoDefaultTrueSecFracMesonNarrow = 		(TH1D*)fileCorrectionsSecPi0->Get("TrueSecFracNarrow");
   TH1D*	histoDefaultTrueSecFracFromK0SMeson = 		(TH1D*)fileCorrectionsSecPi0->Get("TrueSecFracFromK0S");
   TH1D*	histoDefaultTrueSecFracFromK0SMesonWide = (TH1D*)fileCorrectionsSecPi0->Get("TrueSecFracFromK0SWide");
   TH1D*	histoDefaultTrueSecFracFromK0SMesonNarrow = (TH1D*)fileCorrectionsSecPi0->Get("TrueSecFracFromK0SNarrow");
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
   
   TFile* fileCorrections = 			new TFile(fileNameCorrectionFile);
   TH1F *histoEventQualityMC = 			(TH1F*)fileCorrections->Get("NEvents");
   TH1D *histoEffiPt =	(TH1D*)fileCorrections->Get("MesonEffiPt"); //not yet correct MesonEffiPt
//    TH1D *histoEffiPtBackFit =	(TH1D*)fileCorrections->Get("MesonEffiBackFitPt"); //not yet correct MesonEffiPt
   TH1D *histoEffiNarrowPt = (TH1D*)fileCorrections->Get("MesonNarrowEffiPt");
   TH1D *histoEffiWidePt = (TH1D*)fileCorrections->Get("MesonWideEffiPt");
   TH1D *histoEffiLeftPt = (TH1D*)fileCorrections->Get("MesonLeftEffiPt");
   TH1D *histoEffiLeftNarrowPt = (TH1D*)fileCorrections->Get("MesonLeftNarrowEffiPt");
   TH1D *histoEffiLeftWidePt = (TH1D*)fileCorrections->Get("MesonLeftWideEffiPt");
   TH1D *histoAcceptance= 					(TH1D*)fileCorrections->Get("fMCMesonAccepPt");
   
   TH1D *histoTrueEffiPt = NULL;
   TH1D *histoTrueEffiNarrowPt = NULL;
   TH1D *histoTrueEffiWidePt = NULL;
//    if( optionEnergy.CompareTo("2.76TeV") == 0) {
// 		histoTrueEffiPt =					(TH1D*)fileCorrections->Get("TrueMesonEffiPt_Fit"); 
// 		histoTrueEffiNarrowPt = 		(TH1D*)fileCorrections->Get("TrueMesonNarrowEffiPt_Fit");
// 		histoTrueEffiWidePt = 			(TH1D*)fileCorrections->Get("TrueMesonWideEffiPt_Fit");
// 	} else {
		histoTrueEffiPt =					(TH1D*)fileCorrections->Get("TrueMesonEffiPt"); 
		histoTrueEffiNarrowPt = 		(TH1D*)fileCorrections->Get("TrueMesonNarrowEffiPt");
		histoTrueEffiWidePt = 			(TH1D*)fileCorrections->Get("TrueMesonWideEffiPt");
// 	}
	TH1D*	histoInputMesonPt = 				(TH1D*)fileCorrections->Get("MC_Meson_genPt");
	TH1D*	histoInputMesonOldBinPt = 		(TH1D*)fileCorrections->Get("MC_Meson_genPt_oldBin");
   TH1D* histoInputMesonOldBinPtWOWeights = NULL;
   TH1D* histoInputMesonOldBinPtWeights = NULL;
   histoInputMesonOldBinPtWOWeights =     (TH1D*)fileCorrections->Get("MC_Meson_genPt_WOWeights");
   histoInputMesonOldBinPtWeights =     (TH1D*)fileCorrections->Get("MC_Meson_genPt_Weights");
   
   
	TH1D* histoTrueMassMeson =				(TH1D*)fileCorrections->Get("histoTrueMassMeson");
	TH1D* histoTrueFWHMMeson = 			(TH1D*)fileCorrections->Get("histoTrueFWHMMeson");

   TString centralityCutNumber = cutSelection(0,3);
   TString centralityString = GetCentralityString(cutSelection);

   TH1D *histoTrueEffiPtFixed = (TH1D*)histoTrueEffiPt->Clone("histoTrueEffiPtFixed");
   TH1D *histoTrueEffiNarrowPtFixed = (TH1D*)histoTrueEffiNarrowPt->Clone("TrueMesonNarrowEffiPtFixed");
   TH1D *histoTrueEffiWidePtFixed = (TH1D*)histoTrueEffiWidePt->Clone("TrueMesonWideEffiPt");
   TH1D *histoEffiPtFixed = (TH1D*)histoEffiPt->Clone("histoEffiPtFixed");
   TH1D *histoEffiNarrowPtFixed = (TH1D*)histoEffiNarrowPt->Clone("MesonNarrowEffiPtFixed");
   TH1D *histoEffiWidePtFixed = (TH1D*)histoEffiWidePt->Clone("MesonWideEffiPt");

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
   TH1D* histoYieldSecMesonLeft = NULL;
   TH1D* histoYieldSecMeson = NULL;
//    TH1D* histoYieldSecMesonBackFit = NULL;
   TH1D* histoYieldSecMesonLeftNarrow = NULL;
   TH1D* histoYieldSecMesonNarrow = NULL;
   TH1D* histoYieldSecMesonLeftWide = NULL;
   TH1D* histoYieldSecMesonWide = NULL;
   TH1D* histoYieldSecFromK0SMeson = NULL;
//    TH1D* histoYieldSecFromK0SMesonBackFit = NULL;
   TH1D* histoYieldSecFromK0SMesonLeft = NULL;
   TH1D* histoYieldSecFromK0SMesonLeftNarrow = NULL;
   TH1D* histoYieldSecFromK0SMesonLeftWide = NULL;
   TH1D* histoYieldSecFromK0SMesonNarrow = NULL;
   TH1D* histoYieldSecFromK0SMesonWide = NULL;
   TH1D *histoYieldTrueGGFracMeson = NULL;
   TH1D *histoYieldTrueGGFracMesonWide = NULL;
   TH1D *histoYieldTrueGGFracMesonNarrow = NULL;
   TH1D* histoYieldGGMesonLeft = NULL;
   TH1D* histoYieldGGMeson = NULL;
//    TH1D* histoYieldGGMesonBackFit;
   TH1D* histoYieldGGMesonLeftNarrow = NULL;
   TH1D* histoYieldGGMesonNarrow = NULL;
   TH1D* histoYieldGGMesonLeftWide = NULL;
   TH1D* histoYieldGGMesonWide = NULL;
   if (optDalitz.CompareTo("kFALSE")==0){
      histoYieldTrueSecFracMeson = 					(TH1D*)fileCorrections->Get("TrueSecFrac");
      histoYieldTrueSecFracMesonWide = 			(TH1D*)fileCorrections->Get("TrueSecFracWide");
      histoYieldTrueSecFracMesonNarrow = 			(TH1D*)fileCorrections->Get("TrueSecFracNarrow");
      histoYieldTrueSecFracFromK0SMeson = 		(TH1D*)fileCorrections->Get("TrueSecFracFromK0S");
      histoYieldTrueSecFracFromK0SMesonWide = 	(TH1D*)fileCorrections->Get("TrueSecFracFromK0SWide");
      histoYieldTrueSecFracFromK0SMesonNarrow = (TH1D*)fileCorrections->Get("TrueSecFracFromK0SNarrow");
      if ((nameMeson.CompareTo("Pi0")==0 || nameMeson.CompareTo("Pi0EtaBinning")==0) && isMC.CompareTo("kTRUE") !=0){
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
      }
      histoYieldSecMeson = (TH1D*)histoUnCorrectedYield->Clone("SecFracMeson");
      histoYieldSecMeson->Sumw2();
      histoYieldSecMeson->Multiply(histoYieldTrueSecFracMeson);
//       histoYieldSecMesonBackFit = (TH1D*)histoUnCorrectedYieldBackFit->Clone("SecFracMesonBackFit");
//       histoYieldSecMesonBackFit->Sumw2();
//       histoYieldSecMesonBackFit->Multiply(histoYieldTrueSecFracMeson);    
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
		
//       histoYieldSecFromK0SMesonBackFit = (TH1D*)histoUnCorrectedYieldBackFit->Clone("SecFracFromK0SMesonBackFit");
//       histoYieldSecFromK0SMesonBackFit->Sumw2();
//       histoYieldSecFromK0SMesonBackFit->Multiply(histoYieldTrueSecFracFromK0SMeson);
//       histoYieldSecFromK0SMesonBackFit->Scale(doubleAddFactorK0s);
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
      histoYieldTrueGGFracMeson = 					(TH1D*)fileCorrections->Get("TrueGGFrac");
      histoYieldTrueGGFracMesonWide = 			(TH1D*)fileCorrections->Get("TrueGGFracWide");
      histoYieldTrueGGFracMesonNarrow = 			(TH1D*)fileCorrections->Get("TrueGGFracNarrow");
      histoYieldGGMeson = (TH1D*)histoUnCorrectedYield->Clone("GGFracMeson");
      histoYieldGGMeson->Sumw2();
//       histoYieldGGMesonBackFit = (TH1D*)histoUnCorrectedYieldBackFit->Clone("GGFracMesonBackFit");
//       histoYieldGGMesonBackFit->Sumw2();
      histoYieldGGMeson->Multiply(histoYieldTrueGGFracMeson);
      histoYieldGGMesonLeft = (TH1D*)histoUnCorrectedYieldLeft->Clone("GGFracMesonLeft");
      histoYieldGGMesonLeft->Sumw2();
      histoYieldGGMesonLeft->Multiply(histoYieldTrueGGFracMeson);
      histoYieldGGMesonNarrow = (TH1D*)histoUnCorrectedYieldNarrow->Clone("GGFracMesonNarrow");
      histoYieldGGMesonNarrow->Sumw2();
      histoYieldGGMesonNarrow->Multiply(histoYieldTrueGGFracMesonNarrow);
      histoYieldGGMesonLeftNarrow = (TH1D*)histoUnCorrectedYieldLeftNarrow->Clone("GGFracMesonLeftNarrow");
      histoYieldGGMesonLeftNarrow->Sumw2();
      histoYieldGGMesonLeftNarrow->Multiply(histoYieldTrueGGFracMesonNarrow);
      histoYieldGGMesonWide = (TH1D*)histoUnCorrectedYieldWide->Clone("GGFracMesonWide");
      histoYieldGGMesonWide->Sumw2();
      histoYieldGGMesonWide->Multiply(histoYieldTrueGGFracMesonWide);
      histoYieldGGMesonLeftWide = (TH1D*)histoUnCorrectedYieldLeftWide->Clone("GGFracMesonLeftWide");
      histoYieldGGMesonLeftWide->Sumw2();
      histoYieldGGMesonLeftWide->Multiply(histoYieldTrueGGFracMesonWide);
   }
	
   Double_t mesonMassExpect = 0;
   if(nameMeson.CompareTo("Pi0") == 0  || nameMeson.CompareTo("Pi0EtaBinning") == 0 ) mesonMassExpect = TDatabasePDG::Instance()->GetParticle(111)->Mass();
   if(nameMeson.CompareTo("Eta") == 0 ) mesonMassExpect = TDatabasePDG::Instance()->GetParticle(221)->Mass();
	
	
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
                             "", "#it{p}_{T} (GeV/#it{c})", Form("Mass for %s in |#it{y}| < %s (GeV/#it{c}^{2})",textMeson.Data(), rapidityRange.Data()), 
                             kFALSE, 0., 0.7, kFALSE,
                             kFALSE, 0., 0.7, 
                             kFALSE, 0., 10.);
   DrawGammaSetMarker(histoMassMeson, 22, 0.8, kBlack, kBlack);					 
   histoMassMeson->DrawCopy("same,e1,p"); 
	
   DrawGammaSetMarker(histoTrueMassMeson, 24, 0.8, kRed+2, kRed+2);
   histoTrueMassMeson->DrawCopy("same,e1,p"); 
	
   DrawGammaSetMarker(histoMassMesonLeft, 20, 0.8, kBlue, kBlue);
   histoMassMesonLeft->DrawCopy("same,e1,p"); 
	
   DrawGammaLines(0., maxPtMeson,mesonMassExpect, mesonMassExpect,0.1);
	
   TLegend* legendMass = new TLegend(0.15,0.12,0.5,0.25);
   legendMass->SetTextSize(0.02);			
   legendMass->SetFillColor(0);
   legendMass->AddEntry(histoMassMeson,Form("normal  %s", cutSelection.Data()));
   legendMass->AddEntry(histoMassMesonLeft,"left norm");
   if(nameMeson.CompareTo("Pi0") == 0  || nameMeson.CompareTo("Pi0EtaBinning") == 0 ) legendMass->AddEntry(histoTrueMassMeson,"True reconstructed #pi^{0}");
   if(nameMeson.CompareTo("Eta") == 0 ) legendMass->AddEntry(histoTrueMassMeson,"True reconstructed #eta");
	
   legendMass->Draw();
	
   if (fThesis.CompareTo("thesis") == 0)DrawAliceLogoPi0Performance(pictDrawingCoordinatesMass[0], pictDrawingCoordinatesMass[1], pictDrawingCoordinatesMass[2], pictDrawingCoordinatesMass[3], pictDrawingCoordinatesMass[4], pictDrawingCoordinatesMass[5], pictDrawingCoordinatesMass[6], pictDrawingCoordinatesMass[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900,date,"MinBias",kDalitz);
	
   canvasMass->Update();
	
   canvasMass->SaveAs(Form("%s/%s_%s_Mass_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),cutSelection.Data(),suffix));
   delete canvasMass;
   delete legendMass;
   //**********************************************************************************
   //******************** FWHM Plot *********************************************
   //**********************************************************************************
   TCanvas* canvasFWHM = new TCanvas("canvasFWHM","",200,10,1350,900);  // gives the page size
   DrawGammaCanvasSettings( canvasFWHM, 0.13, 0.02, 0.04, 0.09);
	
   histoFWHMMeson->Sumw2();
   histoFWHMMeson->Scale(1./2.36);

   DrawAutoGammaMesonHistos( histoFWHMMeson, 
                             "", "p_{T} (GeV/c)", Form("FWHM/2.36 for %s in |y| < %s (GeV/c^{2})",textMeson.Data(), rapidityRange.Data()), 
                             kFALSE, 1.5,-20., kFALSE,
                             kTRUE, -0.004, 0.020, 
                             kFALSE, 0., 10.);	
	
   histoFWHMMeson->GetYaxis()->SetNdivisions(510); 
   DrawGammaSetMarker(histoFWHMMeson, 22, 0.8, kBlack, kBlack);					 
   histoFWHMMeson->DrawCopy("same,e1,p"); 
	
   histoTrueFWHMMeson->Sumw2();
   histoTrueFWHMMeson->Scale(1./2.36);					 
   DrawGammaSetMarker(histoTrueFWHMMeson, 24, 0.8, kRed+2, kRed+2);
   histoTrueFWHMMeson->DrawCopy("same,e1,p"); 
	
   DrawGammaSetMarker(histoFWHMMesonLeft, 20, 0.8, kBlue, kBlue);
   histoFWHMMesonLeft->Sumw2();
   histoFWHMMesonLeft->Scale(1./2.36);
   histoFWHMMesonLeft->DrawCopy("same,e1,p"); 
	
   TLegend* legendFWHM = new TLegend(0.15,0.1,0.5,0.2);
   legendFWHM->SetTextSize(0.02);			
   legendFWHM->SetFillColor(0);
   legendFWHM->AddEntry(histoFWHMMeson,"normal");
   legendFWHM->AddEntry(histoFWHMMesonLeft,"left norm");
   if(nameMeson.CompareTo("Pi0") == 0 || nameMeson.CompareTo("Pi0EtaBinning") == 0 ) legendFWHM->AddEntry(histoTrueFWHMMeson,"True reconstructed #pi^{0}");
   if(nameMeson.CompareTo("Eta") == 0 ) legendFWHM->AddEntry(histoTrueFWHMMeson,"True reconstructed #eta");
   legendFWHM->Draw();
   if (fThesis.CompareTo("thesis") == 0)DrawAliceLogoPi0Performance(pictDrawingCoordinatesFWHM[0], pictDrawingCoordinatesFWHM[1], pictDrawingCoordinatesFWHM[2], pictDrawingCoordinatesFWHM[3], pictDrawingCoordinatesFWHM[4], pictDrawingCoordinatesFWHM[5], pictDrawingCoordinatesFWHM[6], pictDrawingCoordinatesFWHM[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900,date,"MinBias",kDalitz);
	
   canvasFWHM->Update();
	
   canvasFWHM->SaveAs(Form("%s/%s_%s_FWHM_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),cutSelection.Data(),suffix));
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
   // 		 "[0]*exp(-pow(x-[1],2)/(pow([2],2)))+[3]*log(pow(x,0.5)-[4])*(1-exp([5]-x))"
   fitTrueEffi->SetRange(minPtMesonEffFit,maxPtMesonEffFit);
   // 		fitTrueEffi->SetRange(0.3,maxPtMesonEffFit);
   // 		fitTrueEffi->SetParameter(2,1-fitTrueEffiHighPtCut[i]->GetParameter(0)   );
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
   // 		fitTrueEffiNarrow->SetParameter(2,1-fitTrueEffiNarrowHighPtCut[i]->GetParameter(0)   );
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
   // 		fitTrueEffiWide->SetParameter(2,1-fitTrueEffiWideHighPtCut[i]->GetParameter(0)   );
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
		
      if (fThesis.CompareTo("thesis") == 0)DrawAliceLogoPi0Performance(pictDrawingCoordinatesAcc[0], pictDrawingCoordinatesAcc[1], pictDrawingCoordinatesAcc[2], pictDrawingCoordinatesAcc[3], pictDrawingCoordinatesAcc[4], pictDrawingCoordinatesAcc[5], pictDrawingCoordinatesAcc[6], pictDrawingCoordinatesAcc[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900,date,"MinBias",kDalitz);
      canvasAcceptance->Update();
         
      canvasAcceptance->SaveAs(Form("%s/%s_Acceptance_%s.%s",outputDir.Data(),nameMeson.Data(),cutSelection.Data(),suffix));
      delete canvasAcceptance;
         
      if ( (nameMeson.CompareTo("Pi0")==0 || nameMeson.CompareTo("Pi0EtaBinning")==0) && isMC.CompareTo("kTRUE") !=0){
         if (!kDalitz){
            //**********************************************************************************
            //******************** Secondary Fraction		 **************************************
            //**********************************************************************************
            TCanvas* canvasSecFrac = new TCanvas("canvasSecFrac","",200,10,1350,900);  // gives the page size
            DrawGammaCanvasSettings( canvasSecFrac, 0.09, 0.02, 0.04, 0.09);
            //       canvasSecFrac->SetLogy(1);	
								
            DrawAutoGammaMesonHistos( histoYieldTrueSecFracMeson, 
                                      "", "p_{T} (GeV/c)", "#frac{X->#pi^{0}}{#pi^{0}}", 
                                      kTRUE, 1.5, 0, kFALSE,
                                      kFALSE, 0., 0.7, 
                                      kFALSE, 0., 10.);
            histoYieldTrueSecFracMeson->GetYaxis()->SetTitleOffset(0.9);		
            DrawGammaSetMarker(histoYieldTrueSecFracMeson, 22, 1., kBlack, kBlack);										 
            fitDefaultSecFrac->SetLineColor(kBlack);	
            DrawGammaSetMarker(histoYieldTrueSecFracFromK0SMeson, 22, 1., kBlue, kBlue);										 
            fitDefaultSecFracFromK0->SetLineColor(kBlue);
				
            histoYieldTrueSecFracMeson->DrawCopy("e1"); 	
            fitDefaultSecFrac->Draw("same");
            histoYieldTrueSecFracFromK0SMeson->DrawCopy("e1,same"); 	
            fitDefaultSecFracFromK0->Draw("same");
				
            TLegend* legendSecFrac = new TLegend(0.6,0.8,0.98,0.96);
            legendSecFrac->SetTextSize(0.03);			
            legendSecFrac->SetFillColor(0);
            legendSecFrac->AddEntry(histoYieldTrueSecFracMeson,"X= All Particles");
            legendSecFrac->AddEntry(fitDefaultSecFrac,"fit on total fraction");
            legendSecFrac->AddEntry(histoYieldTrueSecFracFromK0SMeson,"X = K_{s}^{0}");
            legendSecFrac->AddEntry(fitDefaultSecFracFromK0,"fit on fraction from K_{s}^{0}");
            legendSecFrac->Draw();

            /*     if (fThesis.CompareTo("thesis") == 0)DrawAliceLogoPi0Performance(pictDrawingCoordinatesAcc[0], pictDrawingCoordinatesAcc[1], pictDrawingCoordinatesAcc[2], pictDrawingCoordinatesAcc[3], pictDrawingCoordinatesAcc[4], pictDrawingCoordinatesAcc[5], pictDrawingCoordinatesAcc[6], pictDrawingCoordinatesAcc[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900,date,"MinBias",kDalitz);
             */     canvasSecFrac->Update();
				
            canvasSecFrac->SaveAs(Form("%s/%s_FracSecondaries_%s.%s",outputDir.Data(),nameMeson.Data(),cutSelection.Data(),suffix));
            delete canvasSecFrac;
         }
      }
		
      //**********************************************************************************
      //******************** Efficiency Simple Plot **************************************
      //**********************************************************************************
      TCanvas* canvasEffSimple = new TCanvas("canvasEffSimple","",200,10,1350,900);  // gives the page size
      DrawGammaCanvasSettings( canvasEffSimple, 0.13, 0.02, 0.03, 0.09);
      //       canvasEffSimple->SetLogy(1);	
						
      DrawAutoGammaMesonHistos( histoTrueEffiPt, 
                                "", "p_{T} (GeV/c)", "True Efficiency", 
                                kTRUE, 2., 3e-6, kFALSE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
				
      DrawGammaSetMarker(histoTrueEffiPt, 22, 1., kBlack, kBlack);	
      DrawGammaSetMarker(histoTrueEffiPtFixed, 22, 1., kRed, kRed);	
      histoTrueEffiPt->DrawCopy("e1"); 	
      histoTrueEffiPtFixed->DrawCopy("e1same"); 	
      if (fThesis.CompareTo("thesis") == 0)DrawAliceLogoPi0Performance(pictDrawingCoordinatesAcc[0], pictDrawingCoordinatesAcc[1], pictDrawingCoordinatesAcc[2], pictDrawingCoordinatesAcc[3], pictDrawingCoordinatesAcc[4], pictDrawingCoordinatesAcc[5], pictDrawingCoordinatesAcc[6], pictDrawingCoordinatesAcc[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900,date,"MinBias",kDalitz);
      canvasEffSimple->Update();
		
      canvasEffSimple->SaveAs(Form("%s/%s_TrueEffSimple_%s.%s",outputDir.Data(),nameMeson.Data(),cutSelection.Data(),suffix));
      delete canvasEffSimple;
		

      //*********************************************************************************
      //********************** True Efficiency Plot ******************************************
      //*********************************************************************************
		
		

      TCanvas* canvasEffi = new TCanvas("canvasEffi","",1350,1200);  // gives the page size
      DrawGammaCanvasSettings( canvasEffi, 0.13, 0.02, 0.03, 0.09);

      DrawAutoGammaMesonHistos( histoTrueEffiPt, 
                                "", "p_{T} (GeV/c)", "True Efficiency", 
                                kTRUE, 0.75, 3e-6, kFALSE,
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
		
      if (optionEnergy.CompareTo("PbPb_2.76TeV")==0 ){			//|| optionEnergy.CompareTo("2.76TeV")==0
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


      if (fThesis.CompareTo("thesis") == 0)DrawAliceLogoPi0Performance(pictDrawingCoordinatesEff[0], pictDrawingCoordinatesEff[1], pictDrawingCoordinatesEff[2], pictDrawingCoordinatesEff[3], pictDrawingCoordinatesEff[4], pictDrawingCoordinatesEff[5], pictDrawingCoordinatesEff[6], pictDrawingCoordinatesEff[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,1005,date,"MinBias",kDalitz);

      canvasEffi->Update();

      canvasEffi->SaveAs(Form("%s/%s_%s_TrueEfficiency_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),cutSelection.Data(),suffix));
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
   ScaleMCYield(histoMCYieldMesonOldBin,  deltaRapid,  scaling,  nEvtMC,  nameMeson ,optDalitz);
   DrawAutoGammaMesonHistos( histoMCYieldMesonOldBin, 
                           "", "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev.}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}",
                           kFALSE, 3., 4e-10, kTRUE,
                           kFALSE, 0., 0.7, 
                           kFALSE, 0., 10.);
   if (histoInputMesonOldBinPtWOWeights){
      ScaleMCYield(histoInputMesonOldBinPtWOWeights,  deltaRapid,  scaling,  nEvtMC,  nameMeson ,optDalitz);
      histoInputMesonOldBinPtWOWeights->SetName("MCYield_Meson_oldBinWOWeights");
   }   
   TH1D *histoMCYieldMeson = (TH1D*)histoInputMesonPt->Clone();
   histoMCYieldMeson->SetName("MCYield_Meson");
   ScaleMCYield(histoMCYieldMeson,  deltaRapid,  scaling,  nEvtMC,  nameMeson ,optDalitz);
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
   
   canvasMCYieldMeson->SaveAs(Form("%s/%s_histoMCYieldMeson_%s.%s",outputDir.Data(),nameMeson.Data(),cutSelection.Data(),suffix));
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

   if (fThesis.CompareTo("thesis") == 0)DrawAliceLogoPi0Performance(pictDrawingCoordinatesInv[0], pictDrawingCoordinatesInv[1], pictDrawingCoordinatesInv[2], pictDrawingCoordinatesInv[3], pictDrawingCoordinatesInv[4], pictDrawingCoordinatesInv[5], pictDrawingCoordinatesInv[6], pictDrawingCoordinatesInv[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900,date,"MinBias",kDalitz);
   canvasInvMassFull->Update();

   canvasInvMassFull->SaveAs(Form("%s/%s_%s_InvMassFullPt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),cutSelection.Data(),suffix));
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

   if (fThesis.CompareTo("thesis") == 0)DrawAliceLogoPi0Performance(pictDrawingCoordinatesInv[0], pictDrawingCoordinatesInv[1], pictDrawingCoordinatesInv[2], pictDrawingCoordinatesInv[3], pictDrawingCoordinatesInv[4], pictDrawingCoordinatesInv[5], pictDrawingCoordinatesInv[6], pictDrawingCoordinatesInv[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900,date,"MinBias",kDalitz);
   canvasInvMassFullWO->Update();

   canvasInvMassFullWO->SaveAs(Form("%s/%s_%s_InvMassFullPtWithoutBG_%s.%s",outputDir.Data(), nameMeson.Data(), prefix2.Data(), cutSelection.Data(), suffix));
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

   if (fThesis.CompareTo("thesis") == 0)DrawAliceLogoPi0Performance(pictDrawingCoordinatesInv[0], pictDrawingCoordinatesInv[1], pictDrawingCoordinatesInv[2], pictDrawingCoordinatesInv[3], pictDrawingCoordinatesInv[4], pictDrawingCoordinatesInv[5], pictDrawingCoordinatesInv[6], pictDrawingCoordinatesInv[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900,date,"MinBias",kDalitz);
   canvasRAWYield->Update();

   canvasRAWYield->SaveAs(Form("%s/%s_%s_RAWYieldPt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),cutSelection.Data(),suffix));
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
	

   if (optDalitz.CompareTo("kFALSE")==0){
      CorrectYield(histoCorrectedYieldNorm, histoYieldSecMeson, histoYieldSecFromK0SMeson, histoEffiPt, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
      CorrectYield(histoCorrectedYieldTrue, histoYieldSecMeson, histoYieldSecFromK0SMeson, histoTrueEffiPt, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
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

   } else {
      CorrectYieldDalitz(histoCorrectedYieldNorm, histoYieldGGMeson, histoEffiPt, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
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
	DrawGammaSetMarker(histoMCYieldMeson, 24, 1., kRed+2, kRed+2);										 
   histoMCYieldMeson->DrawCopy("e1,same"); 
	DrawGammaSetMarker(histoCorrectedYieldNorm, 24, 1., kGreen+2, kGreen+2);										 
   histoCorrectedYieldNorm->DrawCopy("e1,same"); 
	
	
	
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
	legendYield3->AddEntry(histoMCYieldMeson,"MC input");
	legendYield3->AddEntry(histoCorrectedYieldNorm,"normal eff");
   legendYield3->Draw();

   if (fThesis.CompareTo("thesis") == 0)DrawAliceLogoPi0WorkInProgress(pictDrawingCoordinates[0], pictDrawingCoordinates[1], pictDrawingCoordinates[2], pictDrawingCoordinates[3], pictDrawingCoordinates[4], pictDrawingCoordinates[5], pictDrawingCoordinates[6], pictDrawingCoordinates[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,1125,kDalitz);

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

	DrawGammaSetMarker(RatioTrueMCInput, 24, 1., kRed+2, kRed+2);										 
   RatioTrueMCInput->DrawCopy("e1,same"); 
	DrawGammaSetMarker(RatioNormal, 24, 1., kGreen+2, kGreen+2);										 
   RatioNormal->DrawCopy("e1,same"); 

   canvasCorrecftedYield->Update();
   canvasCorrecftedYield->SaveAs(Form("%s/%s_%s_CorrectedYieldTrueEff_%s.%s",outputDir.Data(), nameMeson.Data(), prefix2.Data(),  cutSelection.Data(), suffix));

   delete canvasCorrecftedYield;
   delete legendYield3;
   


   //***********************************************************************************************
   //***************************  Secondary RAW Yield  *********************************************
   //***********************************************************************************************
   if (optDalitz.CompareTo("kFALSE")==0 && doubleAddFactorK0s >= 0){
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

      // 		if (fThesis.CompareTo("thesis") == 0)DrawAliceLogoPi0Performance(pictDrawingCoordinatesInv[0], pictDrawingCoordinatesInv[1], pictDrawingCoordinatesInv[2], pictDrawingCoordinatesInv[3], pictDrawingCoordinatesInv[4], pictDrawingCoordinatesInv[5], pictDrawingCoordinatesInv[6], pictDrawingCoordinatesInv[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900,date,"MinBias",kDalitz);

      TLegend* legendSecRAWYield = new TLegend(0.6,0.8,0.97,0.95);
      legendSecRAWYield->SetTextSize(0.03);			
      legendSecRAWYield->SetFillColor(0);
      legendSecRAWYield->SetBorderSize(0);
      legendSecRAWYield->AddEntry(histoUnCorrectedYieldDrawing,"RAW yield");
      legendSecRAWYield->AddEntry(histoYieldSecMeson,"total secondaries");
      legendSecRAWYield->AddEntry(histoYieldSecFromK0SMeson,"additional secondaries from K^{0}_{s}");
      legendSecRAWYield->Draw();

      canvasRAWYieldSec->Update();
      canvasRAWYieldSec->SaveAs(Form("%s/%s_%s_RAWYieldSecPt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),cutSelection.Data(),suffix));
      delete canvasRAWYieldSec;
   } else if (!optDalitz.CompareTo("kFALSE")==0){
      TCanvas* canvasRAWYieldSec = new TCanvas("canvasRAWYieldSec","",200,10,1350,900);  // gives the page size
      DrawGammaCanvasSettings( canvasRAWYieldSec, 0.13, 0.02, 0.02, 0.09);	
      canvasRAWYieldSec->SetLogy(1);			

      DrawGammaSetMarker(histoUnCorrectedYieldDrawing, 20, 1., kBlack, kBlack);										 
      histoUnCorrectedYieldDrawing->Draw("e1");
      histoYieldGGMeson->Scale(1./nEvt);
      DrawGammaSetMarker(histoYieldGGMeson, 22, 1., kBlue, kBlue);										 
      histoYieldGGMeson->DrawCopy("same,e1"); 	
			
      // 		if (fThesis.CompareTo("thesis") == 0)DrawAliceLogoPi0Performance(pictDrawingCoordinatesInv[0], pictDrawingCoordinatesInv[1], pictDrawingCoordinatesInv[2], pictDrawingCoordinatesInv[3], pictDrawingCoordinatesInv[4], pictDrawingCoordinatesInv[5], pictDrawingCoordinatesInv[6], pictDrawingCoordinatesInv[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,900,date,"MinBias",kDalitz);

      TLegend* legendSecRAWYield = new TLegend(0.6,0.8,0.97,0.95);
      legendSecRAWYield->SetTextSize(0.03);			
      legendSecRAWYield->SetFillColor(0);
      legendSecRAWYield->SetBorderSize(0);
      legendSecRAWYield->AddEntry(histoUnCorrectedYieldDrawing,"RAW yield");
      legendSecRAWYield->AddEntry(histoYieldGGMeson,"Contamination from #gamma#gamma");
      legendSecRAWYield->Draw();

      canvasRAWYieldSec->Update();
      canvasRAWYieldSec->SaveAs(Form("%s/%s_%s_RAWYieldContGGPt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),cutSelection.Data(),suffix));
      delete canvasRAWYieldSec;
   }

   Int_t nBinsPt = 	histoCorrectedYieldTrue->GetNbinsX();
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
      binsMt= 			new Double_t[nBinsPt+1];
      binsMt[0] = 		mesonMassExpect;


      for(Int_t iPt=1;iPt<nBinsPt+1;iPt++){
         binsMt[iPt]=pow((histoCorrectedYieldTrue->GetXaxis()->GetBinUpEdge(iPt)*histoCorrectedYieldTrue->GetXaxis()->GetBinUpEdge(iPt)+mesonMassExpect*mesonMassExpect),0.5);
         cout << "recalculation pt to mt:    " << iPt<<"     pt      "<< histoCorrectedYieldTrue->GetXaxis()->GetBinUpEdge(iPt)<< "      mt    " << binsMt[iPt]<< endl;
      }

      TH1F *deltaMt = 	new TH1F("deltaMt","",nBinsPt,binsMt);

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

      if (fThesis.CompareTo("thesis") == 0)DrawAliceLogoPi0WorkInProgress(pictDrawingCoordinates[0], pictDrawingCoordinates[1], pictDrawingCoordinates[2], pictDrawingCoordinates[3], pictDrawingCoordinates[4], pictDrawingCoordinates[5], pictDrawingCoordinates[6], pictDrawingCoordinates[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3],1350,1125,kDalitz);

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

      canvasCorrecftedYieldMt->SaveAs(Form("%s/%s_%s_CorrectedYield_Mtspectra_TrueEff_%s.%s", outputDir.Data(), nameMeson.Data() ,prefix2.Data() , cutSelection.Data(), suffix));
      delete canvasCorrecftedYieldMt;
      delete legendYield_Mt;
			
      //***********************************************************************************************
      //***************************  correction for yield in Xt bins **********************************
      //***********************************************************************************************
      Double_t * binsXt = NULL;
      binsXt=  new Double_t[nBinsPt+1];
			
      for(Int_t iPt=0;iPt<nBinsPt+1;iPt++){
         binsXt[iPt]=2*histoCorrectedYieldTrue->GetXaxis()->GetBinUpEdge(iPt)/energy;
         cout << "recalculation pt to xt:    " << iPt<<"     pt      "<< histoCorrectedYieldTrue->GetXaxis()->GetBinUpEdge(iPt)<< "      xt    " << binsXt[iPt]<< endl;
      }
	}
   //*************************************************************************************************
   //******************** Output of the systematic Error due to Signal extraction ********************
   //*************************************************************************************************
   Double_t  binsXCenter[50];
   Double_t  binsXWidth[50];
//    Double_t binYValue[50];
//    binsXCenter[0] = 		0;
   binsXWidth[0]=			0.;
//    binYValue[0]=			0.;
   for (Int_t i = 1; i < nBinsPt +1; i++){
      binsXCenter[i] = 	histoCorrectedYieldTrue->GetBinCenter(i);
      binsXWidth[i]= 	histoCorrectedYieldTrue->GetBinWidth(i)/2.;
   }


   SysErrorConversion sysErrNormal[50];
   SysErrorConversion sysErrLeft[50];
   SysErrorConversion sysErrWide[50];
   SysErrorConversion sysErrNarrow[50];
   SysErrorConversion sysErrLeftWide[50];
   SysErrorConversion sysErrLeftNarrow[50];

   for (Int_t i = 1; i < nBinsPt +1; i++){
      if (optionEnergy.CompareTo("PbPb_2.76TeV")==0  ){	//|| optionEnergy.CompareTo("2.76TeV")==0
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
      largestDifferenceNeg[i] = 		0;
      largestDifferencePos[i] = 		0;
      largestDifferenceNegError[i] = 	0;
      largestDifferencePosError[i] = 	0;
      relLargestDifferenceNeg[i] = 		0;
      relLargestDifferencePos[i] = 		0;
      relLargestDifferenceNegError[i] = 	0;
      relLargestDifferencePosError[i] = 	0;
      //Calculate differences
      differenceLeft[i] = 			sysErrLeft[i].value - sysErrNormal[i].value;
      differenceLeftError[i] = 		TMath::Sqrt(TMath::Abs(TMath::Power(sysErrLeft[i].error,2)-TMath::Power(sysErrNormal[i].error,2)));
      differenceLeftNarrow[i] = 		sysErrLeftNarrow[i].value - sysErrNormal[i].value;
      differenceLeftNarrowError[i] = 	TMath::Sqrt(TMath::Abs(TMath::Power(sysErrLeftNarrow[i].error,2)-TMath::Power(sysErrNormal[i].error,2)));
      differenceLeftWide[i] = 			sysErrLeftWide[i].value - sysErrNormal[i].value;
      differenceLeftWideError[i] = 		TMath::Sqrt(TMath::Abs(TMath::Power(sysErrLeftWide[i].error,2)-TMath::Power(sysErrNormal[i].error,2)));
      differenceNarrow[i] = 			sysErrNarrow[i].value - sysErrNormal[i].value;
      differenceNarrowError[i] = 		TMath::Sqrt(TMath::Abs(TMath::Power(sysErrNarrow[i].error,2)-TMath::Power(sysErrNormal[i].error,2)));
      differenceWide[i] = 			sysErrWide[i].value - sysErrNormal[i].value;
      differenceWideError[i] = 		TMath::Sqrt(TMath::Abs(TMath::Power(sysErrWide[i].error,2)-TMath::Power(sysErrNormal[i].error,2)));

      if (sysErrNormal[i].value != 0){
         relDifferenceLeft[i] = 		(differenceLeft[i]/sysErrNormal[i].value)*100.;
         relDifferenceLeftError[i] = 	(differenceLeftError[i]/sysErrNormal[i].value)*100.;
         relDifferenceLeftNarrow[i] = 	(differenceLeftNarrow[i]/sysErrNormal[i].value)*100.;
         relDifferenceLeftNarrowError[i] = (differenceLeftNarrowError[i]/sysErrNormal[i].value)*100.;
         relDifferenceLeftWide[i] = 	(differenceLeftWide[i]/sysErrNormal[i].value)*100.;
         relDifferenceLeftWideError[i] = (differenceLeftWideError[i]/sysErrNormal[i].value)*100.;
         relDifferenceNarrow[i] = 	(differenceNarrow[i]/sysErrNormal[i].value)*100.;
         relDifferenceNarrowError[i] = (differenceNarrowError[i]/sysErrNormal[i].value)*100.;
         relDifferenceWide[i] = 		(differenceWide[i]/sysErrNormal[i].value)*100.;
         relDifferenceWideError[i] = 	(differenceWideError[i]/sysErrNormal[i].value)*100.;
      } else {
         relDifferenceLeft[i] = 		0.;
         relDifferenceLeftError[i] = 	0.;
         relDifferenceLeftNarrow[i] = 	0.;
         relDifferenceLeftNarrowError[i] = 0.;
         relDifferenceLeftWide[i] = 	0.;
         relDifferenceLeftWideError[i] = 0.;
         relDifferenceNarrow[i] =	 	0.;
         relDifferenceNarrowError[i] = 0.;
         relDifferenceWide[i] = 		0.;
         relDifferenceWideError[i] = 	0.;
      }

      //Find biggest Deviation
      cout << "new bin " << i << endl;
      cout << TMath::Abs(relDifferenceLeft[i]) << endl;
      if (TMath::Abs(relDifferenceLeft[i]) < 75. ){
         if(differenceLeft[i] < 0){
            largestDifferenceNeg[i] = 		differenceLeft[i];
            largestDifferenceNegError[i] = 	differenceLeftError[i];
            relLargestDifferenceNeg[i] = 		relDifferenceLeft[i];
            relLargestDifferenceNegError[i] = 	relDifferenceLeftError[i];
         }else{	
            largestDifferencePos[i] = 		differenceLeft[i]; 
            largestDifferencePosError[i] = 	differenceLeftError[i]; 
            relLargestDifferencePos[i] = 		relDifferenceLeft[i]; 
            relLargestDifferencePosError[i] = 	relDifferenceLeftError[i]; 
         }	
      }
      //       cout << TMath::Abs(relDifferenceLeftNarrow[i]) << endl;
      // 		if (TMath::Abs(relDifferenceLeftNarrow[i]) < 75. ){
      // 			if(differenceLeftNarrow[i] < 0){
      // 				if(differenceLeftNarrow[i] < largestDifferenceNeg[i]){
      // 					largestDifferenceNeg[i] = 		differenceLeftNarrow[i];	
      // 					largestDifferenceNegError[i] = 	differenceLeftNarrowError[i];	
      // 					relLargestDifferenceNeg[i] = 		relDifferenceLeftNarrow[i];	
      // 					relLargestDifferenceNegError[i] = 	relDifferenceLeftNarrowError[i];	
      // 
      // 				}
      // 			} else {
      // 				if(differenceLeftNarrow[i] > largestDifferencePos[i]){
      // 					largestDifferencePos[i] = 		differenceLeftNarrow[i];
      // 					largestDifferencePosError[i] = 	differenceLeftNarrowError[i];
      // 					relLargestDifferencePos[i] = 		relDifferenceLeftNarrow[i];
      // 					relLargestDifferencePosError[i] = 	relDifferenceLeftNarrowError[i];
      // 				}
      // 			}
      //       cout << TMath::Abs(relDifferenceLeftWide[i]) << endl;
      // 		if (TMath::Abs(relDifferenceLeftWide[i]) < 75. ){
      // 			if(differenceLeftWide[i] < 0){
      // 				if(differenceLeftWide[i] < largestDifferenceNeg[i]){
      // 					largestDifferenceNeg[i] = 		differenceLeftWide[i];	
      // 					largestDifferenceNegError[i] = 	differenceLeftWideError[i];	
      // 					relLargestDifferenceNeg[i] = 		relDifferenceLeftWide[i];
      // 					relLargestDifferenceNegError[i] = 	relDifferenceLeftWideError[i];
      // 
      // 				}
      // 			} else {
      // 				if(differenceLeftWide[i] > largestDifferencePos[i]){
      // 					largestDifferencePos[i] = 		differenceLeftWide[i];	
      // 					largestDifferencePosError[i] = 	differenceLeftWideError[i];	
      // 					relLargestDifferencePos[i] = 		relDifferenceLeftWide[i];	
      // 					relLargestDifferencePosError[i] = 	relDifferenceLeftWideError[i];	
      // 				}
      // 			}
      // 		}
      cout << TMath::Abs(relDifferenceNarrow[i]) << endl;
      if (TMath::Abs(relDifferenceNarrow[i]) < 75.){
         if(differenceNarrow[i] < 0){
            if(differenceNarrow[i] < largestDifferenceNeg[i]){
               largestDifferenceNeg[i] = 		differenceNarrow[i];	
               largestDifferenceNegError[i] = 	differenceNarrowError[i];	
               relLargestDifferenceNeg[i] = 		relDifferenceNarrow[i];	
               relLargestDifferenceNegError[i] = 	relDifferenceNarrowError[i];	
            }
         } else {
            if(differenceNarrow[i] > largestDifferencePos[i]){
               largestDifferencePos[i] = 		differenceNarrow[i];
               largestDifferencePosError[i] = 	differenceNarrowError[i];
               relLargestDifferencePos[i] = 		relDifferenceNarrow[i];
               relLargestDifferencePosError[i] = 	relDifferenceNarrowError[i];
            }
         }	
      }
      cout << TMath::Abs(relDifferenceWide[i]) << endl;
      if (TMath::Abs(relDifferenceWide[i]) < 75.){
         if(differenceWide[i] < 0){
            if(differenceWide[i] < largestDifferenceNeg[i]){
               largestDifferenceNeg[i] = 		differenceWide[i];
               largestDifferenceNegError[i] = 	differenceWideError[i];
               relLargestDifferenceNeg[i] = 		relDifferenceWide[i];
               relLargestDifferenceNegError[i] = 	relDifferenceWideError[i];
            }
         } else {
            if(differenceWide[i] > largestDifferencePos[i]) {
               largestDifferencePos[i] = 		differenceWide[i];	
               largestDifferencePosError[i] = 	differenceWideError[i];	
               relLargestDifferencePos[i] = 		relDifferenceWide[i];	
               relLargestDifferencePosError[i] = 	relDifferenceWideError[i];	
            }
         }
      }
   }	

   cout << "done filling" << endl;
   const char *nameFileSysErrDat = Form("%s/%s/%s_%s_SystematicErrorYieldExtraction_%s.dat",cutSelection.Data(),optionEnergy.Data() ,nameMeson.Data(),prefix2.Data(),cutSelection.Data());
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
      fileSysErrDat << i << "\t" << sysErrLeftNarrow[i].value << "\t" << sysErrLeftNarrow[i].error << "\t" << differenceLeftNarrow[i] << 	"\t" << differenceLeftNarrowError[i] << "\t" << relDifferenceLeftNarrow[i] << "\t" << relDifferenceLeftNarrowError[i] <<endl;
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

   const char* nameOutput = Form("%s/%s/%s_%s_GammaConvV1%sCorrection_%s.root",cutSelection.Data(),optionEnergy.Data(),nameMeson.Data(),prefix2.Data(),fDalitz.Data(),cutSelection.Data());
   TFile* correctedOutput = new TFile(nameOutput,"RECREATE");		

   histoCorrectedYieldNorm->Write();
//    histoCorrectedYieldNormBackFit->Write();

   SystErrGraphPos->Write(Form("%s_SystErrorRelPos_YieldExtraction_%s",nameMeson.Data(),centralityString.Data()),TObject::kOverwrite);
   SystErrGraphNeg->Write(Form("%s_SystErrorRelNeg_YieldExtraction_%s",nameMeson.Data(),centralityString.Data()),TObject::kOverwrite);

   histoCorrectedYieldTrue->Write();
//    histoCorrectedYieldTrueBackFit->Write();
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
   if (histoInputMesonOldBinPtWeights) histoInputMesonOldBinPtWeights->Write("WeightsPi0");
   //    histoRatioComparisonPi0MCtoRec->Write();
   //    histoRatioComparisonPi0MCtoRecTrue->Write();
   histoUnCorrectedYieldDrawing->SetName("histoYieldMesonPerEvent");
   histoUnCorrectedYieldDrawing->Write();
   deltaPt->Write("deltaPt");
   correctedOutput->Write();
   correctedOutput->Close();
   
   if (histoEffiNarrowPt || histoEffiWidePt || histoEffiLeftPt || histoEffiLeftNarrowPt || histoEffiLeftWidePt || multFlag){}
}
