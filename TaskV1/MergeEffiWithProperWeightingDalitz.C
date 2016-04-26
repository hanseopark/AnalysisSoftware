/****************************************************************************************************************************
 ****** 	provided by Gamma Conversion Group, PWG4, 														*****
 ******		Ana Marin, marin@physi.uni-heidelberg.de													*****
 ******	   	Kathrin Koch, kkoch@physi.uni-heidelberg.de 													*****
 ******		Friederike Bock, friederike.bock@cern.ch													*****
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
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "ExtractSignalDalitz.h"
#include "../CommonHeaders/ExtractSignalBinningDalitz.h"
#include "../CommonHeaders/ExtractSignalPlotting.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "THnSparse.h"



void MergeEffiWithProperWeightingDalitz(TString sourcePath, TString fCutSelection, TString mesonType, TString fSuffix, TString fEnergyFlag, TString nameFileCorrectionFileFull, TString nameFileCorrectionFileBC, TString nameFileCorrectionFileD){
  
	
   TString nameFilesSamples[6] = {Form("%s/mergedALL/GammaConvV1Data.root",sourcePath.Data()), Form("%s/mergedBC/GammaConvV1Data.root",sourcePath.Data()), Form("%s/mergedDE/GammaConvV1Data.root",sourcePath.Data()), Form("%s/mergedALL/GammaConvV1MC.root",sourcePath.Data()), Form("%s/mergedBC/GammaConvV1MC.root",sourcePath.Data()), Form("%s/mergedDE/GammaConvV1MC.root",sourcePath.Data())};
   
   
  if (fEnergyFlag.CompareTo("2.76TeV")==0 || fEnergyFlag.CompareTo("HI")==0){
		nameFilesSamples[0] = Form("%s/mergedAll/GammaConvV1Data.root",sourcePath.Data());
		nameFilesSamples[1] = Form("%s/mergedAll/GammaConvV1Data.root",sourcePath.Data());
		nameFilesSamples[2] = Form("%s/mergedAll/GammaConvV1Data.root",sourcePath.Data());
		nameFilesSamples[3] = Form("%s/mergedAll/GammaConvV1MC.root",sourcePath.Data());
		nameFilesSamples[4] = Form("%s/WithoutAddedSignals/GammaConvV1MC.root",sourcePath.Data());
		nameFilesSamples[5] = Form("%s/WithAddedSignals/GammaConvV1MC.root",sourcePath.Data());
  }
   TString outputDir = Form("%s/%s/%s/ExtractSignalDalitz",fCutSelection.Data(),fEnergyFlag.Data(),fSuffix.Data());
   gSystem->Exec("mkdir "+outputDir);
   
   //canvasEffSimple->SaveAs(Form("%s/%s_ComparisonEffiMerged_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));

   fstream 	fFileDataLog; 
   const char* fFileDataLogname = Form("%s_%s_ComparisonEffiMerged_%s.dat",outputDir.Data(),mesonType.Data(),fCutSelection.Data());
   fFileDataLog.open(fFileDataLogname, ios::out);

   
   Double_t numberEvents[6];
   
   for (Int_t nFiles = 0; nFiles < 6; nFiles ++){
     
      fileInput = new TFile(nameFilesSamples[nFiles]);
      TList *TopDir =(TList*)fileInput->Get("GammaConvDalitzV1");   

      TList *HistosGammaConversion = (TList*)TopDir->FindObject(Form("Cut Number %s",fCutSelection.Data()));
      TList *ESDContainer = (TList*) HistosGammaConversion->FindObject(Form("%s ESD histograms",fCutSelection.Data()));   
      
      TH1F* fEventQuality = (TH1F*)ESDContainer->FindObject("NEvents");
		if (fEnergyFlag.CompareTo("HI")==0){
			numberEvents[nFiles] = fEventQuality->GetBinContent(1);
		} else {
					      
			numberEvents[nFiles] =  GetNEvents(fEventQuality);
			//numberEvents[nFiles] = fEventQuality->GetEntries()-fEventQuality->GetBinContent(5)-fEventQuality->GetBinContent(7)-fEventQuality->GetBinContent(4);
		}
      delete fileInput;
      cout << numberEvents[nFiles] << endl;
   }	
   Double_t ratioData[2];
   ratioData[0] = numberEvents[1]/numberEvents[0];
   ratioData[1] = numberEvents[2]/numberEvents[0];
	
   Double_t ratioMC[2];
   ratioMC[0] = numberEvents[4]/numberEvents[3];
   ratioMC[1] = numberEvents[5]/numberEvents[3];
	
   fFileDataLog << "Ratios Data:  BC/Full =" << ratioData[0] << "  DE/Full ="  << ratioData[1] << endl;
   fFileDataLog << "Ratios MC:    BC/Full =" << ratioMC[0] << "    DE/Full ="  << ratioMC[1] << endl;
	
   Double_t weighting[2];
   weighting[0] = ratioData[0]/ratioMC[0];
   weighting[1] = ratioData[1]/ratioMC[1];
   
   fFileDataLog << "Weighting MC : BC =" << weighting[0] << "  DE = "  << weighting[1] << endl;
   fFileDataLog << "Corrected Ratio MC : BC/Full=" << ratioMC[0]*weighting[0] << "  DE/Full ="  << ratioMC[1]*weighting[1] << endl;
      
	
   TH1D *histoEffiWeighted;
   TH1D *histoEffiNarrowWeighted;
   TH1D *histoEffiWideWeighted;
   TH1D *histoEffiLeftWeighted;
   TH1D *histoEffiLeftNarrowWeighted;
   TH1D *histoEffiLeftWideWeighted;
   TH1D *histoTrueEffiWeighted;
   TH1D *histoTrueEffiNarrowWeighted;
   TH1D *histoTrueEffiWideWeighted;
   TH1D *histoMCGammaPurityWeighted;
   TH1D *histoMCGammaTruePurityWeighted;
   TH1D *histoGammaMCPrimaryRecoEffWeighted;
   TH1D *histoGammaMCPrimaryRecoEffMCPtWeighted;
   // 	TH1D *histoMCGammaBackgroundWeighted;
   TH1D *histoMCGammaConvProbWeighted;
   TH1D *histoGammaRecoEffWeighted;
   TH1D *histoSecFracMesonWeighted;
   TH1D *histoSecFracMesonWideWeighted;
   TH1D *histoSecFracMesonNarrowWeighted;
   TH1D *histoSecFracFromK0SMesonWeighted;
   TH1D *histoSecFracFromK0SMesonWideWeighted;
   TH1D *histoSecFracFromK0SMesonNarrowWeighted;
   TH1D *fHistoFracAllGammaToSecWeighted;
   TH1D *fHistoFracAllGammaToSecOldBinWeighted;
   TH1D *fHistoFracAllGammaToSecFromXFromK0sWeighted;
   TH1D *fHistoFracAllGammaToSecFromXFromK0sOldBinWeighted;
	
   file = new TFile (nameFileCorrectionFileFull);
   TH1D *histoEffiPt =					(TH1D*)file->Get("MesonEffiPt"); //not yet correct MesonEffiPt
   TH1D *histoEffiNarrowPt = 				(TH1D*)file->Get("MesonNarrowEffiPt");
   TH1D *histoEffiWidePt = 				(TH1D*)file->Get("MesonWideEffiPt");
   TH1D *histoEffiLeftPt = 				(TH1D*)file->Get("MesonLeftEffiPt");
   TH1D *histoEffiLeftNarrowPt = 			(TH1D*)file->Get("MesonLeftNarrowEffiPt");
   TH1D *histoEffiLeftWidePt = 				(TH1D*)file->Get("MesonLeftWideEffiPt");
   TH1D *histoTrueEffiPt =				(TH1D*)file->Get("TrueMesonEffiPt"); //not yet correct MesonEffiPt
   TH1D *histoTrueEffiNarrowPt = 			(TH1D*)file->Get("TrueMesonNarrowEffiPt");
   TH1D *histoTrueEffiWidePt = 				(TH1D*)file->Get("TrueMesonWideEffiPt");
   TH1D *histoYieldTrueGGFracMeson = 			(TH1D*)file->Get("TrueGGFrac");
   TH1D *histoYieldTrueGGFracMesonWide = 		(TH1D*)file->Get("TrueGGFracWide");
   TH1D *histoYieldTrueGGFracMesonNarrow = 		(TH1D*)file->Get("TrueGGFracNarrow");
   
   
   
   /*TH1D *histoMCGammaPurity =				(TH1D*)file->Get("fMCGammaPurity"); //not yet correct MesonEffiPt
   TH1D *histoMCGammaTruePurity = 			(TH1D*)file->Get("fMCGammaTruePurity"); //not yet correct MesonEffiPt;
   TH1D *histoGammaMCPrimaryRecoEff = 			(TH1D*)file->Get("fMCGammaPrimaryRecoEff"); //not yet correct MesonEffiPt;
   TH1D *histoGammaMCPrimaryRecoEffMCPt = 		(TH1D*)file->Get("fMCGammaPrimaryRecoEffMCPt"); //not yet correct MesonEffiPt;
   // 	TH1D *histoMCGammaBackground =				(TH1D*)file->Get("fMCGammaBackground"); //not yet correct MesonEffiPt
   // 	histoMCGammaBackground->Scale(1./numberEvents[4]);
   TH1D *histoMCGammaConvProb = 			(TH1D*)file->Get("fMCGammaConvProb");
   TH1D *histoGammaRecoEff = 				(TH1D*)file->Get("fMCGammaRecoEff");
   TH1D *histoYieldTrueSecFracMeson = 			(TH1D*)file->Get("TrueSecFrac");
   TH1D *histoYieldTrueSecFracMesonWide = 		(TH1D*)file->Get("TrueSecFracWide");
   TH1D *histoYieldTrueSecFracMesonNarrow = 		(TH1D*)file->Get("TrueSecFracNarrow");
   TH1D *histoYieldTrueSecFracFromK0SMeson = 		(TH1D*)file->Get("TrueSecFracFromK0S");
   TH1D *histoYieldTrueSecFracFromK0SMesonWide = 	(TH1D*)file->Get("TrueSecFracFromK0SWide");
   TH1D *histoYieldTrueSecFracFromK0SMesonNarrow = 	(TH1D*)file->Get("TrueSecFracFromK0SNarrow");
   TH1D *fHistoFracAllGammaToSec = (TH1D*)file->Get("FracAllGammaToSec");
   TH1D *fHistoFracAllGammaToSecOldBin = (TH1D*)file->Get("FracAllGammaToSecOldBin");
   TH1D *fHistoFracAllGammaToSecFromXFromK0s = (TH1D*)file->Get("FracAllGammaToSecFromXFromK0s");
   TH1D *fHistoFracAllGammaToSecFromXFromK0sOldBin = (TH1D*)file->Get("FracAllGammaToSecFromXFromK0sOldBin");*/

   fileBC = new TFile (nameFileCorrectionFileBC);
   
   TH1D *histoEffiPtBC =				(TH1D*)fileBC->Get("MesonEffiPt"); //not yet correct MesonEffiPt
   TH1D *histoEffiNarrowPtBC = 				(TH1D*)fileBC->Get("MesonNarrowEffiPt");
   TH1D *histoEffiWidePtBC = 				(TH1D*)fileBC->Get("MesonWideEffiPt");
   TH1D *histoEffiLeftPtBC = 				(TH1D*)fileBC->Get("MesonLeftEffiPt");
   TH1D *histoEffiLeftNarrowPtBC = 			(TH1D*)fileBC->Get("MesonLeftNarrowEffiPt");
   TH1D *histoEffiLeftWidePtBC = 			(TH1D*)fileBC->Get("MesonLeftWideEffiPt");
   TH1D *histoTrueEffiPtBC =				(TH1D*)fileBC->Get("TrueMesonEffiPt"); //not yet correct MesonEffiPt
   TH1D *histoTrueEffiNarrowPtBC = 			(TH1D*)fileBC->Get("TrueMesonNarrowEffiPt");
   TH1D *histoTrueEffiWidePtBC = 			(TH1D*)fileBC->Get("TrueMesonWideEffiPt");
   
   TH1D *histoYieldTrueGGFracMesonBC = 			(TH1D*)fileBC->Get("TrueGGFrac");
   TH1D *histoYieldTrueGGFracMesonWideBC = 		(TH1D*)fileBC->Get("TrueGGFracWide");
   TH1D *histoYieldTrueGGFracMesonNarrowBC = 		(TH1D*)fileBC->Get("TrueGGFracNarrow");
   
   
   /*TH1D *histoMCGammaPurityBC =				(TH1D*)fileBC->Get("fMCGammaPurity"); //not yet correct MesonEffiPt
   TH1D *histoMCGammaTruePurityBC = 			(TH1D*)fileBC->Get("fMCGammaTruePurity"); //not yet correct MesonEffiPt;
   TH1D *histoGammaMCPrimaryRecoEffBC = 			(TH1D*)fileBC->Get("fMCGammaPrimaryRecoEff"); //not yet correct MesonEffiPt;
   TH1D *histoGammaMCPrimaryRecoEffMCPtBC = 			(TH1D*)fileBC->Get("fMCGammaPrimaryRecoEffMCPt"); //not yet correct MesonEffiPt;
   // 	TH1D *histoMCGammaBackgroundBC =				(TH1D*)fileBC->Get("fMCGammaBackground"); //not yet correct MesonEffiPt
   // 	histoMCGammaBackgroundBC->Scale(1./numberEvents[4]);
   TH1D *histoMCGammaConvProbBC = 			(TH1D*)fileBC->Get("fMCGammaConvProb");
   TH1D *histoGammaRecoEffBC = 				(TH1D*)fileBC->Get("fMCGammaRecoEff");
   TH1D *histoYieldTrueSecFracMesonBC = 					(TH1D*)fileBC->Get("TrueSecFrac");
   TH1D *histoYieldTrueSecFracMesonWideBC = 			(TH1D*)fileBC->Get("TrueSecFracWide");
   TH1D *histoYieldTrueSecFracMesonNarrowBC = 			(TH1D*)fileBC->Get("TrueSecFracNarrow");
   TH1D *histoYieldTrueSecFracFromK0SMesonBC = 		(TH1D*)fileBC->Get("TrueSecFracFromK0S");
   TH1D *histoYieldTrueSecFracFromK0SMesonWideBC = 	(TH1D*)fileBC->Get("TrueSecFracFromK0SWide");
   TH1D *histoYieldTrueSecFracFromK0SMesonNarrowBC = (TH1D*)fileBC->Get("TrueSecFracFromK0SNarrow");
   TH1D *fHistoFracAllGammaToSecBC = (TH1D*)fileBC->Get("FracAllGammaToSec");
   TH1D *fHistoFracAllGammaToSecOldBinBC = (TH1D*)fileBC->Get("FracAllGammaToSecOldBin");
   TH1D *fHistoFracAllGammaToSecFromXFromK0sBC = (TH1D*)fileBC->Get("FracAllGammaToSecFromXFromK0s");
   TH1D *fHistoFracAllGammaToSecFromXFromK0sOldBinBC = (TH1D*)fileBC->Get("FracAllGammaToSecFromXFromK0sOldBin");*/
	
   fileD = new TFile (nameFileCorrectionFileD);
   TH1D *histoEffiPtD =					(TH1D*)fileD->Get("MesonEffiPt"); //not yet correct MesonEffiPt
   TH1D *histoEffiNarrowPtD = 				(TH1D*)fileD->Get("MesonNarrowEffiPt");
   TH1D *histoEffiWidePtD = 				(TH1D*)fileD->Get("MesonWideEffiPt");
   TH1D *histoEffiLeftPtD = 				(TH1D*)fileD->Get("MesonLeftEffiPt");
   TH1D *histoEffiLeftNarrowPtD = 			(TH1D*)fileD->Get("MesonLeftNarrowEffiPt");
   TH1D *histoEffiLeftWidePtD = 			(TH1D*)fileD->Get("MesonLeftWideEffiPt");
   TH1D *histoTrueEffiPtD =				(TH1D*)fileD->Get("TrueMesonEffiPt"); //not yet correct MesonEffiPt
   TH1D *histoTrueEffiNarrowPtD = 			(TH1D*)fileD->Get("TrueMesonNarrowEffiPt");
   TH1D *histoTrueEffiWidePtD = 			(TH1D*)fileD->Get("TrueMesonWideEffiPt");
   TH1D *histoYieldTrueGGFracMesonD = 			(TH1D*)fileBC->Get("TrueGGFrac");
   TH1D *histoYieldTrueGGFracMesonWideD = 		(TH1D*)fileBC->Get("TrueGGFracWide");
   TH1D *histoYieldTrueGGFracMesonNarrowD = 		(TH1D*)fileBC->Get("TrueGGFracNarrow");
   
   /*TH1D *histoMCGammaPurityD =				(TH1D*)fileD->Get("fMCGammaPurity"); //not yet correct MesonEffiPt
   TH1D *histoMCGammaTruePurityD = 			(TH1D*)fileD->Get("fMCGammaTruePurity"); //not yet correct MesonEffiPt;
   TH1D *histoGammaMCPrimaryRecoEffD = 			(TH1D*)fileD->Get("fMCGammaPrimaryRecoEff"); //not yet correct MesonEffiPt;
   TH1D *histoGammaMCPrimaryRecoEffMCPtD = 			(TH1D*)fileD->Get("fMCGammaPrimaryRecoEffMCPt"); //not yet correct MesonEffiPt;
   // 	TH1D *histoMCGammaBackgroundD =				(TH1D*)fileD->Get("fMCGammaBackground"); //not yet correct MesonEffiPt
   // 	histoMCGammaBackgroundD->Scale(1./numberEvents[5]);
   TH1D *histoMCGammaConvProbD = 			(TH1D*)fileD->Get("fMCGammaConvProb");
   TH1D *histoGammaRecoEffD = 				(TH1D*)fileD->Get("fMCGammaRecoEff");
   TH1D *histoYieldTrueSecFracMesonD = 					(TH1D*)fileD->Get("TrueSecFrac");
   TH1D *histoYieldTrueSecFracMesonWideD = 			(TH1D*)fileD->Get("TrueSecFracWide");
   TH1D *histoYieldTrueSecFracMesonNarrowD = 			(TH1D*)fileD->Get("TrueSecFracNarrow");
   TH1D *histoYieldTrueSecFracFromK0SMesonD = 		(TH1D*)fileD->Get("TrueSecFracFromK0S");
   TH1D *histoYieldTrueSecFracFromK0SMesonWideD = 	(TH1D*)fileD->Get("TrueSecFracFromK0SWide");
   TH1D *histoYieldTrueSecFracFromK0SMesonNarrowD = (TH1D*)fileD->Get("TrueSecFracFromK0SNarrow");
   TH1D *fHistoFracAllGammaToSecD = (TH1D*)fileD->Get("FracAllGammaToSec");
   TH1D *fHistoFracAllGammaToSecOldBinD = (TH1D*)fileD->Get("FracAllGammaToSecOldBin");
   TH1D *fHistoFracAllGammaToSecFromXFromK0sD = (TH1D*)fileD->Get("FracAllGammaToSecFromXFromK0s");
   TH1D *fHistoFracAllGammaToSecFromXFromK0sOldBinD = (TH1D*)fileD->Get("FracAllGammaToSecFromXFromK0sOldBin");*/

	if (fEnergyFlag.CompareTo("2.76TeV")==0 || fEnergyFlag.CompareTo("HI")==0  ){
	  
		histoEffiWeighted= (TH1D*) histoEffiPt->Clone();
		histoEffiNarrowWeighted=(TH1D*) histoEffiLeftNarrowPt->Clone();
		histoEffiWideWeighted=(TH1D*) histoEffiWidePt->Clone();
		histoEffiLeftWeighted=(TH1D*) histoEffiLeftPt->Clone();
		histoEffiLeftNarrowWeighted=(TH1D*) histoEffiLeftNarrowPt->Clone();
		histoEffiLeftWideWeighted=(TH1D*) histoEffiLeftWidePt->Clone();
		histoTrueEffiWeighted=(TH1D*) histoTrueEffiPt->Clone();
		histoTrueEffiNarrowWeighted=(TH1D*) histoTrueEffiNarrowPt->Clone();
		histoTrueEffiWideWeighted=(TH1D*) histoTrueEffiWidePt->Clone();
		
		histoGGFracMesonWeighted=(TH1D*)histoYieldTrueGGFracMeson->Clone();
		histoGGFracMesonWideWeighted=(TH1D*)histoYieldTrueGGFracMesonWide->Clone();
		histoGGFracMesonNarrowWeighted=(TH1D*)histoYieldTrueGGFracMesonNarrow->Clone();
   
		
		/*histoMCGammaPurityWeighted=(TH1D*) histoMCGammaPurity->Clone();
		histoMCGammaTruePurityWeighted=(TH1D*) histoMCGammaTruePurity->Clone();
		histoGammaMCPrimaryRecoEffWeighted=(TH1D*) histoGammaMCPrimaryRecoEff->Clone();
		histoGammaMCPrimaryRecoEffMCPtWeighted=(TH1D*) histoGammaMCPrimaryRecoEffMCPt->Clone();
		histoMCGammaConvProbWeighted=(TH1D*) histoMCGammaConvProb->Clone();
		histoGammaRecoEffWeighted=(TH1D*) histoGammaRecoEff->Clone();
		histoSecFracMesonWeighted=(TH1D*) histoYieldTrueSecFracMesonBC->Clone();
		histoSecFracMesonWideWeighted=(TH1D*) histoYieldTrueSecFracMesonWideBC->Clone();
		histoSecFracMesonNarrowWeighted=(TH1D*) histoYieldTrueSecFracMesonNarrowBC->Clone();
		histoSecFracFromK0SMesonWeighted=(TH1D*) histoYieldTrueSecFracFromK0SMesonBC->Clone();
		histoSecFracFromK0SMesonWideWeighted=(TH1D*) histoYieldTrueSecFracFromK0SMesonWideBC->Clone();
		histoSecFracFromK0SMesonNarrowWeighted=(TH1D*) histoYieldTrueSecFracFromK0SMesonNarrowBC->Clone();
		fHistoFracAllGammaToSecWeighted=(TH1D*) fHistoFracAllGammaToSecBC->Clone();
		fHistoFracAllGammaToSecOldBinWeighted=(TH1D*) fHistoFracAllGammaToSecOldBinBC->Clone();
		fHistoFracAllGammaToSecFromXFromK0sWeighted=(TH1D*) fHistoFracAllGammaToSecFromXFromK0sBC->Clone();
		fHistoFracAllGammaToSecFromXFromK0sOldBinWeighted=(TH1D*) fHistoFracAllGammaToSecFromXFromK0sOldBinBC->Clone();*/
		
		TCanvas* canvasEffSimple = new TCanvas("canvasEffSimple","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasEffSimple, 0.13, 0.02, 0.02, 0.09);
// 		canvasEffSimple->SetLogy(1);	
						
		DrawAutoGammaMesonHistos( histoTrueEffiPtBC, 
										"", "p_{t} (GeV/c)", "True Efficiency", 
										kTRUE, 2., 3e-6, kFALSE,
										kFALSE, 0., 0.7, 
										kFALSE, 0., 10.);
						
		DrawGammaSetMarker(histoTrueEffiPtBC, 22, 1., kBlue, kBlue);										 
		histoTrueEffiPtBC->DrawCopy("e1,same"); 	
		
		DrawGammaSetMarker(histoTrueEffiPtD, 22, 1., kRed, kRed);										 
		histoTrueEffiPtD->DrawCopy("e1,same"); 		
		
		TLegend* legendYield = new TLegend(0.6,0.15,0.95,0.3);
		legendYield->SetTextSize(0.02);			
		legendYield->SetFillColor(0);
		legendYield->AddEntry(histoTrueEffiPtBC,"Without Added Signals");
		legendYield->AddEntry(histoTrueEffiPtD,"With Added Signals");
		legendYield->Draw();

		canvasEffSimple->Update();
		canvasEffSimple->SaveAs(Form("%s/%s_ComparisonMesonDifferentMC_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));
		delete canvasEffSimple;
		
	} else {
	  
		histoEffiWeighted= (TH1D*) histoEffiPtBC->Clone();
		histoEffiWeighted->Scale(ratioData[0]);
		histoEffiWeighted->Add(histoEffiPtD, ratioData[1]);
		
		histoEffiNarrowWeighted=(TH1D*) histoEffiLeftNarrowPtBC->Clone();
		histoEffiNarrowWeighted->Scale(ratioData[0]);
		histoEffiNarrowWeighted->Add(histoEffiNarrowPtD, ratioData[1]);
		
		histoEffiWideWeighted=(TH1D*) histoEffiWidePtBC->Clone();
		histoEffiWideWeighted->Scale(ratioData[0]);
		histoEffiWideWeighted->Add(histoEffiWidePtD, ratioData[1]);
		
		histoEffiLeftWeighted=(TH1D*) histoEffiLeftPtBC->Clone();
		histoEffiLeftWeighted->Scale(ratioData[0]);
		histoEffiLeftWeighted->Add(histoEffiLeftPtD, ratioData[1]);
		
		histoEffiLeftNarrowWeighted=(TH1D*) histoEffiLeftNarrowPtBC->Clone();
		histoEffiLeftNarrowWeighted->Scale(ratioData[0]);
		histoEffiLeftNarrowWeighted->Add(histoEffiLeftNarrowPtD, ratioData[1]);
		
		histoEffiLeftWideWeighted=(TH1D*) histoEffiLeftWidePtBC->Clone();
		histoEffiLeftWideWeighted->Scale(ratioData[0]);
		histoEffiLeftWideWeighted->Add(histoEffiLeftWidePtD, ratioData[1]);
		
		histoTrueEffiWeighted=(TH1D*) histoTrueEffiPtBC->Clone();
		histoTrueEffiWeighted->Scale(ratioData[0]);
		histoTrueEffiWeighted->Add(histoTrueEffiPtD, ratioData[1]);
		
		histoTrueEffiNarrowWeighted=(TH1D*) histoTrueEffiNarrowPtBC->Clone();
		histoTrueEffiNarrowWeighted->Scale(ratioData[0]);
		histoTrueEffiNarrowWeighted->Add(histoTrueEffiNarrowPtD, ratioData[1]);
		
		histoTrueEffiWideWeighted=(TH1D*) histoTrueEffiWidePtBC->Clone();
		histoTrueEffiWideWeighted->Scale(ratioData[0]);
		histoTrueEffiWideWeighted->Add(histoTrueEffiWidePtD, ratioData[1]);
				
		histoGGFracMesonWeighted=(TH1D*)histoYieldTrueGGFracMesonBC->Clone();
		histoGGFracMesonWeighted->Scale(ratioData[0]);
		histoGGFracMesonWeighted->Add(histoYieldTrueGGFracMesonD,ratioData[1]);
			
		histoGGFracMesonWideWeighted=(TH1D*)histoYieldTrueGGFracMesonWideBC->Clone();
		histoGGFracMesonWideWeighted->Scale(ratioData[0]);
		histoGGFracMesonWideWeighted->Add(histoYieldTrueGGFracMesonWideD,ratioData[1]);
		
		histoGGFracMesonNarrowWeighted=(TH1D*)histoYieldTrueGGFracMesonNarrowBC->Clone();
		histoGGFracMesonNarrowWeighted->Scale(ratioData[0]);
		histoGGFracMesonNarrowWeighted->Add(histoYieldTrueGGFracMesonNarrowD,ratioData[1]);
		
		
		/*histoMCGammaPurityWeighted=(TH1D*) histoMCGammaPurityBC->Clone();
		histoMCGammaPurityWeighted->Scale(ratioData[0]);
		histoMCGammaPurityWeighted->Add(histoMCGammaPurityD, ratioData[1]);
		histoMCGammaTruePurityWeighted=(TH1D*) histoMCGammaTruePurityBC->Clone();
		histoMCGammaTruePurityWeighted->Scale(ratioData[0]);
		histoMCGammaTruePurityWeighted->Add(histoMCGammaTruePurityD, ratioData[1]);
		histoGammaMCPrimaryRecoEffWeighted=(TH1D*) histoGammaMCPrimaryRecoEffBC->Clone();
		histoGammaMCPrimaryRecoEffWeighted->Scale(ratioData[0]);
		histoGammaMCPrimaryRecoEffWeighted->Add(histoGammaMCPrimaryRecoEffD, ratioData[1]);
		histoGammaMCPrimaryRecoEffMCPtWeighted=(TH1D*) histoGammaMCPrimaryRecoEffMCPtBC->Clone();
		histoGammaMCPrimaryRecoEffMCPtWeighted->Scale(ratioData[0]);
		histoGammaMCPrimaryRecoEffMCPtWeighted->Add(histoGammaMCPrimaryRecoEffMCPtD, ratioData[1]);
		// 
		// 	histoMCGammaBackgroundWeighted=(TH1D*) histoMCGammaBackgroundBC->Clone();
		// 	histoMCGammaBackgroundWeighted->Scale(ratioData[0]);
		// 	histoMCGammaBackgroundWeighted->Add(histoMCGammaBackgroundD, ratioData[1]);

		histoMCGammaConvProbWeighted=(TH1D*) histoMCGammaConvProbBC->Clone();
		histoMCGammaConvProbWeighted->Scale(ratioData[0]);
		histoMCGammaConvProbWeighted->Add(histoMCGammaConvProbD, ratioData[1]);
		histoGammaRecoEffWeighted=(TH1D*) histoGammaRecoEffBC->Clone();
		histoGammaRecoEffWeighted->Scale(ratioData[0]);
		histoGammaRecoEffWeighted->Add(histoGammaRecoEffD, ratioData[1]);
		
		histoSecFracMesonWeighted=(TH1D*) histoYieldTrueSecFracMesonBC->Clone();
		histoSecFracMesonWeighted->Scale(ratioData[0]);
		histoSecFracMesonWeighted->Add(histoYieldTrueSecFracMesonD, ratioData[1]);
		histoSecFracMesonWideWeighted=(TH1D*) histoYieldTrueSecFracMesonWideBC->Clone();
		histoSecFracMesonWideWeighted->Scale(ratioData[0]);
		histoSecFracMesonWideWeighted->Add(histoYieldTrueSecFracMesonWideD, ratioData[1]);
		histoSecFracMesonNarrowWeighted=(TH1D*) histoYieldTrueSecFracMesonNarrowBC->Clone();
		histoSecFracMesonNarrowWeighted->Scale(ratioData[0]);
		histoSecFracMesonNarrowWeighted->Add(histoYieldTrueSecFracMesonNarrowD, ratioData[1]);
		
		histoSecFracFromK0SMesonWeighted=(TH1D*) histoYieldTrueSecFracFromK0SMesonBC->Clone();
		histoSecFracFromK0SMesonWeighted->Scale(ratioData[0]);
		histoSecFracFromK0SMesonWeighted->Add(histoYieldTrueSecFracFromK0SMesonD, ratioData[1]);
		histoSecFracFromK0SMesonWideWeighted=(TH1D*) histoYieldTrueSecFracFromK0SMesonWideBC->Clone();
		histoSecFracFromK0SMesonWideWeighted->Scale(ratioData[0]);
		histoSecFracFromK0SMesonWideWeighted->Add(histoYieldTrueSecFracFromK0SMesonWideD, ratioData[1]);
		histoSecFracFromK0SMesonNarrowWeighted=(TH1D*) histoYieldTrueSecFracFromK0SMesonNarrowBC->Clone();
		histoSecFracFromK0SMesonNarrowWeighted->Scale(ratioData[0]);
		histoSecFracFromK0SMesonNarrowWeighted->Add(histoYieldTrueSecFracFromK0SMesonNarrowD, ratioData[1]);
		
		fHistoFracAllGammaToSecWeighted=(TH1D*) fHistoFracAllGammaToSecBC->Clone();
		fHistoFracAllGammaToSecWeighted->Scale(ratioData[0]);
		fHistoFracAllGammaToSecWeighted->Add(fHistoFracAllGammaToSecD, ratioData[1]);
		fHistoFracAllGammaToSecOldBinWeighted=(TH1D*) fHistoFracAllGammaToSecOldBinBC->Clone();
		fHistoFracAllGammaToSecOldBinWeighted->Scale(ratioData[0]);
		fHistoFracAllGammaToSecOldBinWeighted->Add(fHistoFracAllGammaToSecOldBinD, ratioData[1]);
		fHistoFracAllGammaToSecFromXFromK0sWeighted=(TH1D*) fHistoFracAllGammaToSecFromXFromK0sBC->Clone();
		fHistoFracAllGammaToSecFromXFromK0sWeighted->Scale(ratioData[0]);
		fHistoFracAllGammaToSecFromXFromK0sWeighted->Add(fHistoFracAllGammaToSecFromXFromK0sD, ratioData[1]);
		fHistoFracAllGammaToSecFromXFromK0sOldBinWeighted=(TH1D*) fHistoFracAllGammaToSecFromXFromK0sOldBinBC->Clone();
		fHistoFracAllGammaToSecFromXFromK0sOldBinWeighted->Scale(ratioData[0]);
		fHistoFracAllGammaToSecFromXFromK0sOldBinWeighted->Add(fHistoFracAllGammaToSecFromXFromK0sOldBinD, ratioData[1]);*/


		//Drawing different Efficiencys
		TCanvas* canvasEffSimple = new TCanvas("canvasEffSimple","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasEffSimple, 0.13, 0.02, 0.02, 0.09);
		canvasEffSimple->SetLogy(1);	
						
		DrawAutoGammaMesonHistos( histoTrueEffiWeighted, 
										"", "p_{t} (GeV/c)", "True Efficiency", 
										kTRUE, 2., 3e-6, kFALSE,
										kFALSE, 0., 0.7, 
										kFALSE, 0., 10.);
				
		DrawGammaSetMarker(histoTrueEffiWeighted, 22, 1., kBlack, kBlack);										 
		histoTrueEffiWeighted->DrawCopy("e1"); 	
		
		DrawGammaSetMarker(histoTrueEffiPtBC, 22, 1., kBlue, kBlue);										 
		histoTrueEffiPtBC->DrawCopy("e1,same"); 	
		
		DrawGammaSetMarker(histoTrueEffiPtD, 22, 1., kRed, kRed);										 
		histoTrueEffiPtD->DrawCopy("e1,same"); 		
		
		TLegend* legendYield = new TLegend(0.6,0.15,0.95,0.3);
		legendYield->SetTextSize(0.02);			
		legendYield->SetFillColor(0);
		legendYield->AddEntry(histoTrueEffiWeighted,"merged with weighting");
		legendYield->AddEntry(histoTrueEffiPtBC,"Period B&C");
		legendYield->AddEntry(histoTrueEffiPtD,"Period D&E");
		legendYield->Draw();
			
			
		//	DrawAliceLogoPi0Performance(pictDrawingCoordinatesAcc[0], pictDrawingCoordinatesAcc[1], pictDrawingCoordinatesAcc[2], pictDrawingCoordinatesAcc[3], pictDrawingCoordinatesAcc[4], pictDrawingCoordinatesAcc[5], pictDrawingCoordinatesAcc[6], pictDrawingCoordinatesAcc[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3]);
		canvasEffSimple->Update();
		canvasEffSimple->SaveAs(Form("%s/%s_ComparisonEffiMerged_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));
		delete canvasEffSimple;

		
		/*TCanvas* canvasEffSimple = new TCanvas("canvasEffSimple","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasEffSimple, 0.13, 0.02, 0.02, 0.09);
		canvasEffSimple->SetLogy(1);	
						
		DrawAutoGammaMesonHistos( histoGammaRecoEffWeighted, 
										"", "p_{t} (GeV/c)", "#gamma reco. Efficiency", 
										kTRUE, 2., 3e-6, kFALSE,
										kFALSE, 0., 0.7, 
										kFALSE, 0., 10.);
				
		DrawGammaSetMarker(histoGammaRecoEffWeighted, 22, 1., kBlack, kBlack);										 
		histoGammaRecoEffWeighted->DrawCopy("e1"); 	
		
		DrawGammaSetMarker(histoGammaRecoEffBC, 22, 1., kBlue, kBlue);										 
		histoGammaRecoEffBC->DrawCopy("e1,same"); 	
		
		DrawGammaSetMarker(histoGammaRecoEffD, 22, 1., kRed, kRed);										 
		histoGammaRecoEffD->DrawCopy("e1,same"); 		
		
		TLegend* legendYield = new TLegend(0.6,0.15,0.95,0.3);
		legendYield->SetTextSize(0.02);			
		legendYield->SetFillColor(0);
		legendYield->AddEntry(histoGammaRecoEffWeighted,"merged with weighting");
		legendYield->AddEntry(histoGammaRecoEffBC,"Period B&C");
		legendYield->AddEntry(histoGammaRecoEffD,"Period D&E");
		legendYield->Draw();
			
			
		//	DrawAliceLogoPi0Performance(pictDrawingCoordinatesAcc[0], pictDrawingCoordinatesAcc[1], pictDrawingCoordinatesAcc[2], pictDrawingCoordinatesAcc[3], pictDrawingCoordinatesAcc[4], pictDrawingCoordinatesAcc[5], pictDrawingCoordinatesAcc[6], pictDrawingCoordinatesAcc[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3]);
		canvasEffSimple->Update();
		canvasEffSimple->SaveAs(Form("%s/%s_ComparisonGammaEffiMerged_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));
		delete canvasEffSimple;*/

		/*TCanvas* canvasEffSimple = new TCanvas("canvasEffSimple","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasEffSimple, 0.13, 0.02, 0.02, 0.09);
		canvasEffSimple->SetLogy(1);	
						
		DrawAutoGammaMesonHistos( histoGammaMCPrimaryRecoEffWeighted, 
										"", "p_{t} (GeV/c)", "#gamma reco. Efficiency", 
										kTRUE, 2., 3e-6, kFALSE,
										kFALSE, 0., 0.7, 
										kFALSE, 0., 10.);
				
		DrawGammaSetMarker(histoGammaMCPrimaryRecoEffWeighted, 22, 1., kBlack, kBlack);										 
		histoGammaMCPrimaryRecoEffWeighted->DrawCopy("e1"); 	
		
		DrawGammaSetMarker(histoGammaMCPrimaryRecoEffBC, 22, 1., kBlue, kBlue);										 
		histoGammaMCPrimaryRecoEffBC->DrawCopy("e1,same"); 	
		
		DrawGammaSetMarker(histoGammaMCPrimaryRecoEffD, 22, 1., kRed, kRed);										 
		histoGammaMCPrimaryRecoEffD->DrawCopy("e1,same"); 		
		
		TLegend* legendYield = new TLegend(0.6,0.15,0.95,0.3);
		legendYield->SetTextSize(0.02);			
		legendYield->SetFillColor(0);
		legendYield->AddEntry(histoGammaMCPrimaryRecoEffWeighted,"merged with weighting");
		legendYield->AddEntry(histoGammaMCPrimaryRecoEffBC,"Period B&C");
		legendYield->AddEntry(histoGammaMCPrimaryRecoEffD,"Period D&E");
		legendYield->Draw();
			
			
		//	DrawAliceLogoPi0Performance(pictDrawingCoordinatesAcc[0], pictDrawingCoordinatesAcc[1], pictDrawingCoordinatesAcc[2], pictDrawingCoordinatesAcc[3], pictDrawingCoordinatesAcc[4], pictDrawingCoordinatesAcc[5], pictDrawingCoordinatesAcc[6], pictDrawingCoordinatesAcc[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3]);
		canvasEffSimple->Update();
		canvasEffSimple->SaveAs(Form("%s/%s_ComparisonGammaPrimaryEffiMerged_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));
		delete canvasEffSimple;


		TCanvas* canvasEffSimple = new TCanvas("canvasEffSimple","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasEffSimple, 0.13, 0.02, 0.02, 0.09);
		canvasEffSimple->SetLogy(1);	
						
		DrawAutoGammaMesonHistos( histoGammaMCPrimaryRecoEffMCPtWeighted, 
										"", "p_{t} (GeV/c)", "#gamma reco. Efficiency", 
										kTRUE, 2., 3e-6, kFALSE,
										kFALSE, 0., 0.7, 
										kFALSE, 0., 10.);
				
		DrawGammaSetMarker(histoGammaMCPrimaryRecoEffMCPtWeighted, 22, 1., kBlack, kBlack);										 
		histoGammaMCPrimaryRecoEffMCPtWeighted->DrawCopy("e1"); 	
		
		DrawGammaSetMarker(histoGammaMCPrimaryRecoEffMCPtBC, 22, 1., kBlue, kBlue);										 
		histoGammaMCPrimaryRecoEffMCPtBC->DrawCopy("e1,same"); 	
		
		DrawGammaSetMarker(histoGammaMCPrimaryRecoEffMCPtD, 22, 1., kRed, kRed);										 
		histoGammaMCPrimaryRecoEffMCPtD->DrawCopy("e1,same"); 		
		
		TLegend* legendYield = new TLegend(0.6,0.15,0.95,0.3);
		legendYield->SetTextSize(0.02);			
		legendYield->SetFillColor(0);
		legendYield->AddEntry(histoGammaMCPrimaryRecoEffMCPtWeighted,"merged with weighting");
		legendYield->AddEntry(histoGammaMCPrimaryRecoEffMCPtBC,"Period B&C");
		legendYield->AddEntry(histoGammaMCPrimaryRecoEffMCPtD,"Period D&E");
		legendYield->Draw();
			
			
		//	DrawAliceLogoPi0Performance(pictDrawingCoordinatesAcc[0], pictDrawingCoordinatesAcc[1], pictDrawingCoordinatesAcc[2], pictDrawingCoordinatesAcc[3], pictDrawingCoordinatesAcc[4], pictDrawingCoordinatesAcc[5], pictDrawingCoordinatesAcc[6], pictDrawingCoordinatesAcc[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3]);
		canvasEffSimple->Update();
		canvasEffSimple->SaveAs(Form("%s/%s_ComparisonGammaPrimaryEffiMCPtMerged_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));
		delete canvasEffSimple;
		*/


		/*TCanvas* canvasEffSimple = new TCanvas("canvasEffSimple","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasEffSimple, 0.13, 0.02, 0.02, 0.09);	
						
		DrawAutoGammaMesonHistos( histoMCGammaConvProbWeighted, 
										"", "p_{t} (GeV/c)", "#gamma Conv Prob", 
										kFALSE, 2., 3e-6, kFALSE,
										kTRUE, 0.02, 0.14, 
										kFALSE, 0., 10.);
				
		DrawGammaSetMarker(histoMCGammaConvProbWeighted, 22, 1., kBlack, kBlack);										 
		histoMCGammaConvProbWeighted->DrawCopy("e1"); 	
		
		DrawGammaSetMarker(histoMCGammaConvProbBC, 22, 1., kBlue, kBlue);										 
		histoMCGammaConvProbBC->DrawCopy("e1,same"); 	
		
		DrawGammaSetMarker(histoMCGammaConvProbD, 22, 1., kRed, kRed);										 
		histoMCGammaConvProbD->DrawCopy("e1,same"); 		
		
		TLegend* legendYield = new TLegend(0.6,0.15,0.95,0.3);
		legendYield->SetTextSize(0.02);			
		legendYield->SetFillColor(0);
		legendYield->AddEntry(histoMCGammaConvProbWeighted,"merged with weighting");
		legendYield->AddEntry(histoMCGammaConvProbBC,"Period B&C");
		legendYield->AddEntry(histoMCGammaConvProbD,"Period D&E");
		legendYield->Draw();
			
			
		//	DrawAliceLogoPi0Performance(pictDrawingCoordinatesAcc[0], pictDrawingCoordinatesAcc[1], pictDrawingCoordinatesAcc[2], pictDrawingCoordinatesAcc[3], pictDrawingCoordinatesAcc[4], pictDrawingCoordinatesAcc[5], pictDrawingCoordinatesAcc[6], pictDrawingCoordinatesAcc[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3]);
		canvasEffSimple->Update();
		canvasEffSimple->SaveAs(Form("%s/%s_ComparisonGammaProbMerged_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));
		delete canvasEffSimple;
		*/

		/*TCanvas* canvasEffSimple = new TCanvas("canvasEffSimple","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasEffSimple, 0.13, 0.02, 0.02, 0.09);
		canvasEffSimple->SetLogy(1);	
						
		DrawAutoGammaMesonHistos( histoMCGammaPurityWeighted, 
										"", "p_{t} (GeV/c)", "#gamma Purity", 
										kFALSE, 2., 3e-6, kFALSE,
										kTRUE, 0.5, 1.05, 
										kFALSE, 0., 10.);
				
		DrawGammaSetMarker(histoMCGammaPurityWeighted, 22, 1., kBlack, kBlack);										 
		histoMCGammaPurityWeighted->DrawCopy("e1"); 	
		
		DrawGammaSetMarker(histoMCGammaPurityBC, 22, 1., kBlue, kBlue);										 
		histoMCGammaPurityBC->DrawCopy("e1,same"); 	
		
		DrawGammaSetMarker(histoMCGammaPurityD, 22, 1., kRed, kRed);										 
		histoMCGammaPurityD->DrawCopy("e1,same"); 		
		
		TLegend* legendYield = new TLegend(0.6,0.15,0.95,0.3);
		legendYield->SetTextSize(0.02);			
		legendYield->SetFillColor(0);
		legendYield->AddEntry(histoMCGammaPurityWeighted,"merged with weighting");
		legendYield->AddEntry(histoMCGammaPurityBC,"Period B&C");
		legendYield->AddEntry(histoMCGammaPurityD,"Period D&E");
		legendYield->Draw();
			
			
		//	DrawAliceLogoPi0Performance(pictDrawingCoordinatesAcc[0], pictDrawingCoordinatesAcc[1], pictDrawingCoordinatesAcc[2], pictDrawingCoordinatesAcc[3], pictDrawingCoordinatesAcc[4], pictDrawingCoordinatesAcc[5], pictDrawingCoordinatesAcc[6], pictDrawingCoordinatesAcc[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3]);
		canvasEffSimple->Update();
		canvasEffSimple->SaveAs(Form("%s/%s_ComparisonGammaPurityMerged_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));
		delete canvasEffSimple;

		TCanvas* canvasEffSimple = new TCanvas("canvasEffSimple","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasEffSimple, 0.13, 0.02, 0.02, 0.09);
		canvasEffSimple->SetLogy(1);	*/
						
		/*DrawAutoGammaMesonHistos( histoMCGammaTruePurityWeighted, 
										"", "p_{t} (GeV/c)", "#gamma TruePurity", 
										kFALSE, 2., 3e-6, kFALSE,
										kTRUE, 0.5, 1.05, 
										kFALSE, 0., 10.);
				
		DrawGammaSetMarker(histoMCGammaTruePurityWeighted, 22, 1., kBlack, kBlack);										 
		histoMCGammaTruePurityWeighted->DrawCopy("e1"); 	
		
		DrawGammaSetMarker(histoMCGammaTruePurityBC, 22, 1., kBlue, kBlue);										 
		histoMCGammaTruePurityBC->DrawCopy("e1,same"); 	
		
		DrawGammaSetMarker(histoMCGammaTruePurityD, 22, 1., kRed, kRed);										 
		histoMCGammaTruePurityD->DrawCopy("e1,same"); 		
		
		TLegend* legendYield = new TLegend(0.6,0.15,0.95,0.3);
		legendYield->SetTextSize(0.02);			
		legendYield->SetFillColor(0);
		legendYield->AddEntry(histoMCGammaTruePurityWeighted,"merged with weighting");
		legendYield->AddEntry(histoMCGammaTruePurityBC,"Period B&C");
		legendYield->AddEntry(histoMCGammaTruePurityD,"Period D&E");
		legendYield->Draw();
		
		//	DrawAliceLogoPi0Performance(pictDrawingCoordinatesAcc[0], pictDrawingCoordinatesAcc[1], pictDrawingCoordinatesAcc[2], pictDrawingCoordinatesAcc[3], pictDrawingCoordinatesAcc[4], pictDrawingCoordinatesAcc[5], pictDrawingCoordinatesAcc[6], pictDrawingCoordinatesAcc[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3]);
		canvasEffSimple->Update();
		canvasEffSimple->SaveAs(Form("%s/%s_ComparisonGammaTruePurityMerged_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));
		delete canvasEffSimple;*/
	}

   fileSum = new TFile (nameFileCorrectionFileFull,"UPDATE");
   histoEffiWeighted->Write("MesonEffiPt",TObject::kOverwrite);
   histoEffiNarrowWeighted->Write("MesonNarrowEffiPt",TObject::kOverwrite);
   histoEffiWideWeighted->Write("MesonWideEffiPt",TObject::kOverwrite);
   histoEffiLeftWeighted->Write("MesonLeftEffiPt",TObject::kOverwrite);
   histoEffiLeftNarrowWeighted->Write("MesonLeftNarrowEffiPt",TObject::kOverwrite);
   histoEffiLeftWideWeighted->Write("MesonLeftWideEffiPt",TObject::kOverwrite);
   histoTrueEffiWeighted->Write("TrueMesonEffiPt",TObject::kOverwrite);
   histoTrueEffiNarrowWeighted->Write("TrueMesonNarrowEffiPt",TObject::kOverwrite);
   histoTrueEffiWideWeighted->Write("TrueMesonWideEffiPt",TObject::kOverwrite);
   histoGGFracMesonWeighted->Write("TrueGGFrac",TObject::kOverwrite);
   histoGGFracMesonWideWeighted->Write("TrueGGFracWide",TObject::kOverwrite);
   histoGGFracMesonNarrowWeighted->Write("TrueGGFracNarrow",TObject::kOverwrite);
   
   /*histoMCGammaPurityWeighted->Write("fMCGammaPurity",TObject::kOverwrite);
   histoMCGammaTruePurityWeighted->Write("fMCGammaTruePurity",TObject::kOverwrite);
   // 	histoMCGammaBackgroundWeighted->Write("fMCGammaBackground",TObject::kOverwrite);
   histoMCGammaConvProbWeighted->Write("fMCGammaConvProb",TObject::kOverwrite);
   histoGammaRecoEffWeighted->Write("fMCGammaRecoEff",TObject::kOverwrite);
   histoGammaMCPrimaryRecoEffWeighted->Write("fMCGammaPrimaryRecoEff",TObject::kOverwrite);
   histoGammaMCPrimaryRecoEffMCPtWeighted->Write("fMCGammaPrimaryRecoEffMCPt",TObject::kOverwrite);
   histoSecFracMesonWeighted->Write("TrueSecFrac",TObject::kOverwrite);
   histoSecFracMesonNarrowWeighted->Write("TrueSecFracNarrow",TObject::kOverwrite);
   histoSecFracMesonWideWeighted->Write("TrueSecFracWide",TObject::kOverwrite);
   histoSecFracFromK0SMesonWeighted->Write("TrueSecFracFromK0S",TObject::kOverwrite);
   histoSecFracFromK0SMesonNarrowWeighted->Write("TrueSecFracFromK0SNarrow",TObject::kOverwrite);
   histoSecFracFromK0SMesonWideWeighted->Write("TrueSecFracFromK0SWide",TObject::kOverwrite);
   fHistoFracAllGammaToSecD->Write("FracAllGammaToSec",TObject::kOverwrite);
   fHistoFracAllGammaToSecOldBinD->Write("FracAllGammaToSecOldBin",TObject::kOverwrite);
   fHistoFracAllGammaToSecFromXFromK0sD->Write("FracAllGammaToSecFromXFromK0s",TObject::kOverwrite);
   fHistoFracAllGammaToSecFromXFromK0sOldBinD->Write("FracAllGammaToSecFromXFromK0sOldBin",TObject::kOverwrite);*/
        
   fileSum->Write();
   fileSum->Close();
   fFileDataLog.close();
	
	
}
