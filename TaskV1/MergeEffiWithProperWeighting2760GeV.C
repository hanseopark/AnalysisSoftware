/****************************************************************************************************************************
 ****** 	provided by Gamma Conversion Group, PWGGA, 	            													*****
 ******		Friederike Bock, friederike.bock@cern.ch				                									*****
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



void MergeEffiWithProperWeighting2760GeV(   TString fCutSelection, 
                                            TString mesonType, 
                                            TString fSuffix, 
                                            TString fEnergyFlag, 
                                            TString nameFileCorrectionFileFull, 
                                            TString nameFileMinBias, 
                                            TString nameFileAddedSignal
                                        ){
	
	gSystem->Exec("cp "+nameFileCorrectionFileFull+" "+nameFileMinBias);
	TString outputDir 						= Form("%s/%s/%s/ExtractSignal", fCutSelection.Data(), fEnergyFlag.Data(), fSuffix.Data());
	gSystem->Exec("mkdir "+outputDir);
   
	TFile* file 							= new TFile (nameFileCorrectionFileFull);
	TH1D *histoTrueEffiPt 					= (TH1D*)file->Get("TrueMesonEffiPt"); 
	TH1D *histoTrueYield 					= (TH1D*)file->Get("histoYieldTrueMeson"); 
	TH1D *histoTrueEffiNarrowPt 			= (TH1D*)file->Get("TrueMesonNarrowEffiPt");
	TH1D *histoTrueEffiWidePt 				= (TH1D*)file->Get("TrueMesonWideEffiPt");
	TH1D* histoTrueMassMeson 				= (TH1D*)file->Get("histoTrueMassMeson");
	TH1D* histoTrueFWHMMeson 				= (TH1D*)file->Get("histoTrueFWHMMeson");

	
	TFile* fileMinBias 						= new TFile (nameFileMinBias);
	TH1D *histoTrueEffiPtMinBias 			= (TH1D*)fileMinBias->Get("TrueMesonEffiPt"); 
	TH1D *histoTrueYieldMinBias 			= (TH1D*)fileMinBias->Get("histoYieldTrueMeson"); 
	TH1D *histoTrueEffiNarrowPtMinBias 		= (TH1D*)fileMinBias->Get("TrueMesonNarrowEffiPt");
	TH1D *histoTrueEffiWidePtMinBias 		= (TH1D*)fileMinBias->Get("TrueMesonWideEffiPt");
	TH1D* histoTrueMassMesonMinBias 		= (TH1D*)fileMinBias->Get("histoTrueMassMeson");
	TH1D* histoTrueFWHMMesonMinBias 		= (TH1D*)fileMinBias->Get("histoTrueFWHMMeson");
		
	TFile* fileAddedSignal 					= new TFile (nameFileAddedSignal);
	TH1F *histoEventQualityAddSig 			= NULL;
	TH1D *histoMCInputAddedSig 				= NULL;
	TH1D *histoMCInputWOWeightingAddedSig 	= NULL;
	TH1D *histoMCInputWeightsAddedSig 		= NULL;
	histoEventQualityAddSig 				= (TH1F*)fileAddedSignal->Get("NEvents");
	histoMCInputAddedSig 					= (TH1D*)fileAddedSignal->Get("MC_Meson_genPt_oldBin"); 
	histoMCInputWOWeightingAddedSig 		= (TH1D*)fileAddedSignal->Get("MC_Meson_genPt_WOWeights"); 
	histoMCInputWeightsAddedSig 			= (TH1D*)fileAddedSignal->Get("MC_Meson_genPt_Weights"); 
	TH1D *histoTrueYieldAddedSignal 		= (TH1D*)fileAddedSignal->Get("histoYieldTrueMeson"); 
	TH1D *histoTrueEffiPtAddedSignal 		= (TH1D*)fileAddedSignal->Get("TrueMesonEffiPt"); 
	TH1D *histoTrueEffiNarrowPtAddedSignal 	= (TH1D*)fileAddedSignal->Get("TrueMesonNarrowEffiPt");
	TH1D *histoTrueEffiWidePtAddedSignal 	= (TH1D*)fileAddedSignal->Get("TrueMesonWideEffiPt");
	TH1D* histoTrueMassMesonAddedSignal 	= (TH1D*)fileAddedSignal->Get("histoTrueMassMeson");
	TH1D* histoTrueFWHMMesonAddedSignal 	= (TH1D*)fileAddedSignal->Get("histoTrueFWHMMeson");
	
	TH1D *histoTrueEffiPtWeighted 			=	(TH1D*)histoTrueEffiPtAddedSignal->Clone("histoTrueEffiPtWeighted"); 
   
	for (Int_t i = 1; i < histoTrueEffiPtMinBias->GetNbinsX()+1 ; i++){
		Double_t relErrMinBias 		= histoTrueEffiPtMinBias->GetBinError(i)/histoTrueEffiPtMinBias->GetBinContent(i)*100;
		Double_t relErrAddedSignal 	= histoTrueEffiPtAddedSignal->GetBinError(i)/histoTrueEffiPtAddedSignal->GetBinContent(i)*100;
				
		Double_t weightMinBias 		= 1/TMath::Power(histoTrueEffiPtMinBias->GetBinError(i),2);
		Double_t weightAddSig 		= 1/TMath::Power(histoTrueEffiPtAddedSignal->GetBinError(i),2);
		Double_t weightSum 			= weightMinBias + weightAddSig;
		Double_t weightedEffi 		= (weightMinBias*histoTrueEffiPtMinBias->GetBinContent(i) + weightAddSig * histoTrueEffiPtAddedSignal->GetBinContent(i))/weightSum;
		Double_t weightedEffiErr 	= pow((weightMinBias +  weightAddSig),-0.5);
		
		if (isfinite(weightedEffi) && isfinite(weightedEffiErr)){
			histoTrueEffiPt->SetBinContent(i, weightedEffi);
			histoTrueEffiPt->SetBinError(i, weightedEffiErr);
		}
		cout << histoTrueEffiPtMinBias->GetBinContent(i) << "\t" << histoTrueEffiPtAddedSignal->GetBinContent(i)<<"\t" << histoTrueEffiPtWeighted->GetBinContent(i) << "\t" << histoTrueEffiPtWeighted->GetBinError(i) << endl;
		
		cout << "NORMAL: pt: "<< histoTrueEffiPtMinBias->GetBinCenter(i) <<"\t error min Bias: " << relErrMinBias << "\t error added sig: "  << relErrAddedSignal << endl;
		if (relErrAddedSignal < relErrMinBias && isfinite(relErrAddedSignal) ){
//          histoTrueEffiPt->SetBinContent(i,histoTrueEffiPtAddedSignal->GetBinContent(i));
//          histoTrueEffiPt->SetBinError(i,histoTrueEffiPtAddedSignal->GetBinError(i));
         histoTrueMassMeson->SetBinContent(i,histoTrueMassMesonAddedSignal->GetBinContent(i));
         histoTrueMassMeson->SetBinError(i,histoTrueMassMesonAddedSignal->GetBinError(i));
         histoTrueFWHMMeson->SetBinContent(i,histoTrueFWHMMesonAddedSignal->GetBinContent(i));
         histoTrueFWHMMeson->SetBinError(i,histoTrueFWHMMesonAddedSignal->GetBinError(i));
		} else {
// 			histoTrueEffiPt->SetBinContent(i,histoTrueEffiPtMinBias->GetBinContent(i));
//          histoTrueEffiPt->SetBinError(i,histoTrueEffiPtMinBias->GetBinError(i));
         histoTrueMassMeson->SetBinContent(i,histoTrueMassMesonMinBias->GetBinContent(i));
         histoTrueMassMeson->SetBinError(i,histoTrueMassMesonMinBias->GetBinError(i));
         histoTrueFWHMMeson->SetBinContent(i,histoTrueFWHMMesonMinBias->GetBinContent(i));
         histoTrueFWHMMeson->SetBinError(i,histoTrueFWHMMesonMinBias->GetBinError(i));
		}
		Double_t relErrMinBiasNarrow = histoTrueEffiNarrowPtMinBias->GetBinError(i)/histoTrueEffiNarrowPtMinBias->GetBinContent(i)*100;
		Double_t relErrAddedSignalNarrow = histoTrueEffiNarrowPtAddedSignal->GetBinError(i)/histoTrueEffiNarrowPtAddedSignal->GetBinContent(i)*100;
		cout << "NARROW: pt: "<< histoTrueEffiPtMinBias->GetBinCenter(i) <<"\t error min Bias: " << relErrMinBiasNarrow << "\t error added sig: "  << relErrAddedSignalNarrow << endl;
		
		Double_t weightMinBiasNarrow = 1/TMath::Power(histoTrueEffiNarrowPtMinBias->GetBinError(i),2);
		Double_t weightAddSigNarrow = 1/TMath::Power(histoTrueEffiNarrowPtAddedSignal->GetBinError(i),2);
		Double_t weightSumNarrow = weightMinBiasNarrow + weightAddSigNarrow;
		Double_t weightedEffiNarrow = (weightMinBiasNarrow*histoTrueEffiNarrowPtMinBias->GetBinContent(i) + weightAddSigNarrow * histoTrueEffiNarrowPtAddedSignal->GetBinContent(i))/weightSumNarrow;
		Double_t weightedEffiErrNarrow = pow((weightMinBiasNarrow +  weightAddSigNarrow),-0.5);

		if (isfinite(weightedEffiNarrow) && isfinite(weightedEffiErrNarrow)){
			histoTrueEffiNarrowPt->SetBinContent(i, weightedEffiNarrow);
			histoTrueEffiNarrowPt->SetBinError(i, weightedEffiErrNarrow);
		}
		
// 		if ( relErrAddedSignalNarrow< relErrMinBiasNarrow && isfinite(relErrAddedSignalNarrow) ){
// 			histoTrueEffiNarrowPt->SetBinContent(i,histoTrueEffiNarrowPtAddedSignal->GetBinContent(i));
//          histoTrueEffiNarrowPt->SetBinError(i,histoTrueEffiNarrowPtAddedSignal->GetBinError(i));
// 		} else {
// 			histoTrueEffiNarrowPt->SetBinContent(i,histoTrueEffiNarrowPtMinBias->GetBinContent(i));
//          histoTrueEffiNarrowPt->SetBinError(i,histoTrueEffiNarrowPtMinBias->GetBinError(i));
// 		}

		Double_t relErrMinBiasWide = histoTrueEffiWidePtMinBias->GetBinError(i)/histoTrueEffiWidePtMinBias->GetBinContent(i)*100;
		Double_t relErrAddedSignalWide = histoTrueEffiWidePtAddedSignal->GetBinError(i)/histoTrueEffiWidePtAddedSignal->GetBinContent(i)*100;
		cout << "WIDE: : pt: "<< histoTrueEffiPtMinBias->GetBinCenter(i) <<"\t error min Bias: " <<  relErrMinBiasWide << "\t error added sig: "  << relErrAddedSignalWide << endl;
		
		Double_t weightMinBiasWide = 1/TMath::Power(histoTrueEffiWidePtMinBias->GetBinError(i),2);
		Double_t weightAddSigWide = 1/TMath::Power(histoTrueEffiWidePtAddedSignal->GetBinError(i),2);
		Double_t weightSumWide = weightMinBiasWide + weightAddSigWide;
		Double_t weightedEffiWide = (weightMinBiasWide*histoTrueEffiWidePtMinBias->GetBinContent(i) + weightAddSigWide * histoTrueEffiWidePtAddedSignal->GetBinContent(i))/weightSumWide;
		Double_t weightedEffiErrWide = pow((weightMinBiasWide +  weightAddSigWide),-0.5);

		if (isfinite(weightedEffiWide) && isfinite(weightedEffiErrWide)){
			histoTrueEffiWidePt->SetBinContent(i, weightedEffiWide);
			histoTrueEffiWidePt->SetBinError(i, weightedEffiErrWide);
		}

		
// 		if (relErrAddedSignalWide < relErrMinBiasWide  && isfinite(relErrAddedSignalWide) ){
// 			histoTrueEffiWidePt->SetBinContent(i,histoTrueEffiWidePtAddedSignal->GetBinContent(i));
//          histoTrueEffiWidePt->SetBinError(i,histoTrueEffiWidePtAddedSignal->GetBinError(i));
// 		} else {
//          histoTrueEffiWidePt->SetBinContent(i,histoTrueEffiWidePtMinBias->GetBinContent(i));
//          histoTrueEffiWidePt->SetBinError(i,histoTrueEffiWidePtMinBias->GetBinError(i));
// 			
// 		}
	}
	
	TH1D* ratioMinBiasFinal = (TH1D*) histoTrueEffiPtMinBias->Clone("ratioMinBiasFinal");
	ratioMinBiasFinal->Divide(histoTrueEffiPtMinBias,histoTrueEffiPt , 1.,1.,"B");
	TH1D* ratioAddSignalFinal = (TH1D*) histoTrueEffiPtAddedSignal->Clone("ratioAddSignalFinal");
	ratioAddSignalFinal->Divide(histoTrueEffiPtAddedSignal,histoTrueEffiPt , 1.,1.,"B");
	TH1D* ratioMinBiasAddSignal = (TH1D*) histoTrueEffiPtMinBias->Clone("ratioAddSignalFinal");
	ratioMinBiasAddSignal->Divide(histoTrueEffiPtMinBias,histoTrueEffiPtAddedSignal , 1.,1.,"B");

	
	TCanvas* canvasFraction2 = new TCanvas("canvasFraction2","",1550,1200);  // gives the page size
	canvasFraction2->SetTickx();
	canvasFraction2->SetTicky();
	canvasFraction2->SetGridx(0);
	canvasFraction2->SetGridy(0);
	canvasFraction2->SetLogy(0);
	canvasFraction2->SetLeftMargin(0.13);
	canvasFraction2->SetRightMargin(0.02);
	canvasFraction2->SetTopMargin(0.02);
	canvasFraction2->SetFillColor(0);


	DrawGammaSetMarker(ratioMinBiasFinal, 24, 1., kBlue+2, kBlue+2);
	DrawAutoGammaMesonHistos( ratioMinBiasFinal,
					"", "p_{T} (GeV/c)", "effi A/ effi B",
					kFALSE, 5., 10e-10, kTRUE,
					kTRUE, 0.5, 1.5,
					kTRUE, 0., 7.9);
	
	DrawGammaSetMarker(ratioAddSignalFinal, 25, 1., kRed, kRed);
	ratioAddSignalFinal->Draw("same");
	DrawGammaSetMarker(ratioMinBiasAddSignal, 26, 1., kGreen+2, kGreen+2);
	ratioMinBiasAddSignal->Draw("same");
	canvasFraction2->Update();
	DrawGammaLines(0., 8.,1., 1.,0.1);

	TLegend* legendMultDataPP = new TLegend(0.15,0.8,0.4,0.95);
	legendMultDataPP->SetFillColor(0);
	legendMultDataPP->SetLineColor(0);
	legendMultDataPP->SetTextSize(0.04);
	legendMultDataPP->AddEntry(ratioAddSignalFinal,"Added Signal / final","p");
	legendMultDataPP->AddEntry(ratioMinBiasAddSignal,"Min Bias/ added Signal","p");
	legendMultDataPP->AddEntry(ratioMinBiasFinal, "Min Bias/ final","p");
	legendMultDataPP->Draw();

	
	canvasFraction2->SaveAs(Form("%s/%s_MC_RatioComparisonEffiMerged_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));

	
	//Drawing different Efficiencys
	TCanvas* canvasEffSimple = new TCanvas("canvasEffSimple","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasEffSimple, 0.13, 0.02, 0.02, 0.09);
	canvasEffSimple->SetLogy(1);	
					
	DrawAutoGammaMesonHistos( histoTrueEffiPt, 
									"", "p_{t} (GeV/c)", "True Efficiency", 
									kTRUE, 2., 3e-6, kFALSE,
									kFALSE, 0., 0.7, 
									kFALSE, 0., 10.);
			
	DrawGammaSetMarker(histoTrueEffiPt, 24, 1., kBlack, kBlack);										 
	histoTrueEffiPt->DrawCopy("e1"); 	
	
	DrawGammaSetMarker(histoTrueEffiPtMinBias, 25, 1., kBlue, kBlue);										 
	histoTrueEffiPtMinBias->DrawCopy("e1,same"); 	
	
	
	DrawGammaSetMarker(histoTrueEffiPtAddedSignal, 26, 1., kRed, kRed);										 
	histoTrueEffiPtAddedSignal->DrawCopy("e1,same"); 		

// 	DrawGammaSetMarker(histoTrueEffiPtWeighted, 26, 1., kGreen+2, kGreen+2);										 
// 	histoTrueEffiPtWeighted->DrawCopy("e1,same"); 	

	
	TLegend* legendYield = new TLegend(0.6,0.15,0.95,0.3);
	legendYield->SetTextSize(0.02);			
	legendYield->SetFillColor(0);
	legendYield->AddEntry(histoTrueEffiPt,"merged ");
// 	legendYield->AddEntry(histoTrueEffiPtWeighted,"weighted ");
	legendYield->AddEntry(histoTrueEffiPtMinBias,"min Bias");
	legendYield->AddEntry(histoTrueEffiPtAddedSignal,"added Signal");
	legendYield->Draw();
		
		
	//	DrawAliceLogoPi0Performance(pictDrawingCoordinatesAcc[0], pictDrawingCoordinatesAcc[1], pictDrawingCoordinatesAcc[2], pictDrawingCoordinatesAcc[3], pictDrawingCoordinatesAcc[4], pictDrawingCoordinatesAcc[5], pictDrawingCoordinatesAcc[6], pictDrawingCoordinatesAcc[7], pictDrawingCoordinates[8], collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], pictDrawingOptions[3]);
	canvasEffSimple->Update();
	canvasEffSimple->SaveAs(Form("%s/%s_MC_ComparisonEffiMerged_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));
	
		//Drawing different Efficiencys
	canvasEffSimple->SetLogy(0);	
					
	DrawAutoGammaMesonHistos( histoTrueEffiPt, 
									"", "p_{t} (GeV/c)", "True Efficiency", 
									kTRUE, 0.75, 3e-6, kFALSE,
									kFALSE, 0., 0.7, 
									kFALSE, 0., 10.);
			
	histoTrueEffiPt->DrawCopy("e1"); 	
	histoTrueEffiPtMinBias->DrawCopy("e1,same"); 	
// 	histoTrueEffiPtWeighted->DrawCopy("e1,same"); 	
	histoTrueEffiPtAddedSignal->DrawCopy("e1,same"); 			
	legendYield->Draw();
		
	canvasEffSimple->Update();
	canvasEffSimple->SaveAs(Form("%s/%s_MC_ComparisonEffiMerged_linY_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));

	delete canvasEffSimple;

	
	
   TFile* fileSum = new TFile (nameFileCorrectionFileFull,"UPDATE");
	histoTrueEffiPt->Write("TrueMesonEffiPt",TObject::kOverwrite);
	histoTrueEffiNarrowPt->Write("TrueMesonNarrowEffiPt",TObject::kOverwrite);
	histoTrueEffiWidePt->Write("TrueMesonWideEffiPt",TObject::kOverwrite);
	histoTrueMassMeson->Write("histoTrueMassMeson",TObject::kOverwrite);
	histoTrueFWHMMeson->Write("histoTrueFWHMMeson",TObject::kOverwrite);
	histoEventQualityAddSig->Write("NEvents_AddedSig",TObject::kOverwrite);
	histoMCInputAddedSig->Write("MC_Meson_genPt_oldBin_AddedSig",TObject::kOverwrite);
	histoMCInputWOWeightingAddedSig->Write("MC_Meson_genPt_WOWeights_AddedSig",TObject::kOverwrite);
	histoMCInputWeightsAddedSig->Write("MC_Meson_genPt_Weights_AddedSig",TObject::kOverwrite); 
	fileSum->Write();
   fileSum->Close();
		
}
