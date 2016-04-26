/*******************************************************************************
 ****** 	provided by Gamma Conversion Group, PWGGA,								*****
 ******		Ana Marin, marin@physi.uni-heidelberg.de								*****
 ******	  	Martin Wilde, m_wild03@uni-muenster.de 								*****
 ******		Friederike Bock, friederike.bock@cern.ch								*****
 *******************************************************************************/

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
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "ProduceFinalResultspPb.h"

void ProduceFinalResultspPb(TString fileName = "myOutput", TString cutSelection = "", TString suffix = "gif", TString textMeson = "Pi0" ,TString makeBinShiftWithFunction = "", TString optionEnergy = "", TString fileName2= "", Bool_t optNoBinShift=kFALSE, TString optDalitz = "",Int_t mode = 1){
	gROOT->Reset();	
	gROOT->SetStyle("Plain");
	
	if (optionEnergy.CompareTo("pPb_5.023TeV") != 0){
	  cout << "this macro is intended to be used on pPb only, it can not run on "<< optionEnergy.Data()<< endl;  
	  return;
	}
   
	Bool_t kDalitz = kFALSE;
	
        if (optDalitz.CompareTo("kTRUE")==0){
                kDalitz = kTRUE;
        }
   
	StyleSettingsThesis();	
	SetPlotStyle();
	
	date = ReturnDateString();

	minPtForFits=0.4;
	collisionSystem ="p-Pb, #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";
	centralityString = GetCentralityString(cutSelection.Data());
	nColl = GetNCollFromCutNumber(cutSelection.Data());
	nCollErr = GetNCollErrFromCutNumber(cutSelection.Data());
	Int_t offsetSyst = 2;
	Bool_t kMeson = kTRUE;
	if (textMeson.CompareTo("Pi0")==0){
	  
	    if( kDalitz == kFALSE ){
	  
		if (centralityString.CompareTo("0-20%") == 0 ){
			fileNameSysErr = Form("SystematicErrorsNew/SystematicErrorAveraged_%s_DummypPb.dat", textMeson.Data()); 
		} else if (centralityString.CompareTo("20-40%") == 0){
			fileNameSysErr = Form("SystematicErrorsNew/SystematicErrorAveraged_%s_DummypPb.dat", textMeson.Data()); 
		}else if (centralityString.CompareTo("40-60%") == 0){
			fileNameSysErr = Form("SystematicErrorsNew/SystematicErrorAveraged_%s_DummypPb.dat", textMeson.Data()); 
		}else if (centralityString.CompareTo("60-80%") == 0){
			fileNameSysErr = Form("SystematicErrorsNew/SystematicErrorAveraged_%s_DummypPb.dat", textMeson.Data()); 
		} else {
			fileNameSysErr = Form("SystematicErrorsNew/SystematicErrorAveraged_%s_pPb_5.023TeV_MB_29_Apr_2014.dat", textMeson.Data()); 
		}
		offsetSyst = 2;
		kMeson = kTRUE;
	    } else {
	      
	    	  if (centralityString.CompareTo("0-20%") == 0 ){
			  fileNameSysErr = Form("SystematicErrorsNew/SystematicErrorAveraged_%s_DummypPb.dat", textMeson.Data()); 
		  } else if (centralityString.CompareTo("20-40%") == 0){
			  fileNameSysErr = Form("SystematicErrorsNew/SystematicErrorAveraged_%s_DummypPb.dat", textMeson.Data()); 
		  }else if (centralityString.CompareTo("40-60%") == 0){
			  fileNameSysErr = Form("SystematicErrorsNew/SystematicErrorAveraged_%s_DummypPb.dat", textMeson.Data()); 
		  }else if (centralityString.CompareTo("60-80%") == 0){
			  fileNameSysErr = Form("SystematicErrorsNew/SystematicErrorAveraged_%s_DummypPb.dat", textMeson.Data()); 
		  } else {
			  fileNameSysErr = Form("SystematicErrorsNew/SystematicErrorAveraged_Pi0_pPb_5.023TeVMB_27_Mar_2014.dat", textMeson.Data()); 
		  }
		 offsetSyst = 3;
		 kMeson = kTRUE;
			  
	    }
		
	} else {
		if (centralityString.CompareTo("0-20%") == 0 ){
			fileNameSysErr = Form("SystematicErrorsNew/SystematicErrorAveraged_%s_DummypPb.dat", textMeson.Data()); 
		} else if (centralityString.CompareTo("20-40%") == 0){
			fileNameSysErr = Form("SystematicErrorsNew/SystematicErrorAveraged_%s_DummypPb.dat", textMeson.Data()); 
		}else if (centralityString.CompareTo("40-60%") == 0){
			fileNameSysErr = Form("SystematicErrorsNew/SystematicErrorAveraged_%s_DummypPb.dat", textMeson.Data()); 
		}else if (centralityString.CompareTo("60-80%") == 0){
			fileNameSysErr = Form("SystematicErrorsNew/SystematicErrorAveraged_%s_DummypPb.dat", textMeson.Data()); 
		} else {
			fileNameSysErr = Form("SystematicErrorsNew/SystematicErrorAveraged_Eta_pPb_5.023TeVMB_27_Mar_2014.dat"); 
		}
		offsetSyst = 4;
		kMeson = kFALSE;
	}
	cout << fileNameSysErr.Data() << endl;
	
	TString outputDir = Form("%s/%s/%s/ProduceFinalResultspPb",cutSelection.Data(),optionEnergy.Data(),suffix.Data());
	gSystem->Exec("mkdir -p "+outputDir);
	
	TString fEventCutSelection="";
	TString fGammaCutSelection="";
	TString fClusterCutSelection="";
	TString fElectronCutSelection="";
	TString fMesonCutSelection="";
	
	if (mode == 9){
		if (kDalitz){
			ReturnSeparatedCutNumber(cutSelection, fGammaCutSelection, fElectronCutSelection,fMesonCutSelection,kTRUE);
		} else {
			ReturnSeparatedCutNumber(cutSelection, fGammaCutSelection, fElectronCutSelection,fMesonCutSelection);
		}
		
		fEventCutSelection = fGammaCutSelection(0,7);
		fGammaCutSelection = fGammaCutSelection(7,fGammaCutSelection.Length()-7);
		cout << fEventCutSelection.Data() << "\t" << fGammaCutSelection.Data() << endl;
	} else if (kDalitz){
	        ReturnSeparatedCutNumberAdvanced(cutSelection,fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);
	}else {
		ReturnSeparatedCutNumberAdvanced(cutSelection,fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);
	}
   
	deltaEta =  ReturnRapidityStringAndDouble(fMesonCutSelection, rapidityRange);
				
	//declaration for printing logo 
	prefix2 = "data";
	pictDrawingOptions[1] = kFALSE;
	
	if (makeBinShiftWithFunction.CompareTo("Hagedorn")==0){ 
		kHag = kTRUE;
		kLevy = kFALSE;
	}else { 
		if(makeBinShiftWithFunction.CompareTo("Levy")==0){
			kHag = kFALSE;
			kLevy = kTRUE;
		}else {
			cout << "No known fit for Binshift defined, Levy will be taken!" <<endl;
			kHag = kFALSE;
			kLevy = kTRUE;
		}
	}
	if (textMeson.Contains("Pi0")){
		mesonMassExpect = TDatabasePDG::Instance()->GetParticle(111)->Mass();
	} else {
		mesonMassExpect = TDatabasePDG::Instance()->GetParticle(221)->Mass();
	}
	
	// File definitions
	cout << "fileName: " << fileName.Data()<< endl;
	file = 					new TFile(fileName);
	histoCorrectedYield =			(TH1D*)file->Get("CorrectedYieldTrueEff"); 
	cout << "Bins" << endl;
	for (Int_t i = 0; i < histoCorrectedYield->GetNbinsX()+1; i++){
		cout << i<<"\t" <<histoCorrectedYield->GetXaxis()->GetBinCenter(i)<< "\t" <<histoCorrectedYield->GetBinContent(i)<< endl;
	}
	
	histoFWHMMeson = 		(TH1D*)file->Get("histoFWHMMeson");    
	histoMassMeson = 		(TH1D*)file->Get("histoMassMeson");
	histoEventQualtity = 		(TH1F*)file->Get("NEvents");
	histoTrueFWHMMeson = 		(TH1D*)file->Get("histoTrueFWHMMeson");    
	histoTrueMassMeson = 		(TH1D*)file->Get("histoTrueMassMeson");
	TH1D* histoAcceptance= 		(TH1D*)file->Get("fMCMesonAccepPt");
	TH1D* histoTrueEffPt =		(TH1D*)file->Get("TrueMesonEffiPt"); 
	TH1D* histoMCInput =     	(TH1D*)file->Get("MCYield_Meson_oldBin");
	TH1D* histoMCInputWOWeighting = (TH1D*)file->Get("MCYield_Meson_oldBinWOWeights");
	TH1D* histoMCInputWeights =     (TH1D*)file->Get("WeightsMeson");
	TH1D* histoMCInputAddedSig = NULL;
	TH1D* histoMCInputWOWeightingAddedSig = NULL;
	TH1D* histoMCInputWeightsAddedSig = NULL;
	histoMCInputAddedSig =   (TH1D*)file->Get("MCYield_Meson_oldBin_AddedSig");
	histoMCInputWOWeightingAddedSig =   (TH1D*)file->Get("MCYield_Meson_oldBinWOWeights_AddedSig");
	histoMCInputWeightsAddedSig = (TH1D*)file->Get("WeightsMeson_AddedSig");
    
	TH1D* histoBGEstimateFromPileup = (TH1D*)file->Get("BGEstimateFromPileup");
	TH1D* histoUncorrectedYield = 	(TH1D*)file->Get("histoYieldMeson");
	
	histoMassMesonMinusExp= CalculateMassMinusExpectedMass(histoMassMeson,mesonMassExpect);
	Double_t maxMass = 0;
	Double_t minMass = 0;
	for (Int_t i = 1; i < histoMassMesonMinusExp->GetNbinsX(); i++){
	  if (maxMass < histoMassMesonMinusExp->GetBinContent(i)) maxMass = histoMassMesonMinusExp->GetBinContent(i);
	  if (minMass > histoMassMesonMinusExp->GetBinContent(i)) minMass = histoMassMesonMinusExp->GetBinContent(i);
	}
	if (minMass < 0) histoMassMesonMinusExp->GetYaxis()->SetRangeUser(minMass*1.5, maxMass*1.5);
	else histoMassMesonMinusExp->GetYaxis()->SetRangeUser(minMass*0.5, maxMass*1.5);
	
	    histoTrueMassMesonMinusExp = CalculateMassMinusExpectedMass(histoTrueMassMeson,mesonMassExpect);
	maxMass = 0;
	minMass = 0;
	for (Int_t i = 1; i < histoTrueMassMesonMinusExp->GetNbinsX(); i++){
	    if (maxMass < histoTrueMassMesonMinusExp->GetBinContent(i)) maxMass = histoTrueMassMesonMinusExp->GetBinContent(i);
	    if (minMass > histoTrueMassMesonMinusExp->GetBinContent(i)) minMass = histoTrueMassMesonMinusExp->GetBinContent(i);
	}
	if (minMass < 0) histoTrueMassMesonMinusExp->GetYaxis()->SetRangeUser(minMass*1.5, maxMass*1.5);
	else histoTrueMassMesonMinusExp->GetYaxis()->SetRangeUser(minMass*0.5, maxMass*1.5);
 
	histoFWHMMesonMeV = (TH1D*) histoFWHMMeson->Clone();
	histoFWHMMesonMeV->Scale(1000.);
	histoTrueFWHMMesonMeV = (TH1D*) histoTrueFWHMMeson->Clone();
	histoTrueFWHMMesonMeV->Scale(1000.);
	    
	nEvt = histoEventQualtity->GetBinContent(1);
	pictDrawingCoordinates[8] = nEvt;
		
 	fileSysErr.open(fileNameSysErr.Data(),ios_base::in);
	cout << fileNameSysErr.Data() << endl;

	while(!fileSysErr.eof()){
		fileSysErr >> relSystErrorDown[nPoints] >> relSystErrorUp[nPoints]>>	relSystErrorWOMaterialDown[nPoints] >> relSystErrorWOMaterialUp[nPoints];
		cout << nPoints << "\t"  << relSystErrorDown[nPoints] << "\t"  <<relSystErrorUp[nPoints] << "\t" << relSystErrorWOMaterialDown[nPoints] << "\t"  <<relSystErrorWOMaterialUp[nPoints] << endl;;		
		nPoints++;
	}
	fileSysErr.close();
	nPoints = nPoints-1;
	
	maxPt = histoCorrectedYield->GetXaxis()->GetBinUpEdge(histoCorrectedYield->GetNbinsX());
			
	//**********************************************************************************
    //******************** FWHM Plot *********************************************
    //**********************************************************************************
	TCanvas* canvasFWHM = new TCanvas("canvasFWHM","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasFWHM, 0.08, 0.02, 0.04, 0.09);
   
	histoFWHMMeson->GetYaxis()->SetNdivisions(510); 
	DrawGammaSetMarker(histoFWHMMeson, 22, 0.8, kBlack, kBlack);					  
	DrawAutoGammaMesonHistos( histoFWHMMeson, 
								 "", "p_{T} (GeV/c)","FWHM/2.36 (GeV/c^{2})", 
								 kFALSE, 3.,0., kFALSE,
								 kTRUE, -0.004, 0.020, 
								 kFALSE, 0., 10.);
	histoFWHMMeson->GetYaxis()->SetTitleOffset(1.);
	histoFWHMMeson->DrawCopy("e1,p"); 
	DrawGammaSetMarker(histoTrueFWHMMeson, 26, 0.8, kRed+2, kRed+2);
	histoTrueFWHMMeson->DrawCopy("same,e1,p"); 
	
	TLegend* legendFWHM = new TLegend(0.3,0.1,0.7,0.2);
	legendFWHM->SetTextSize(0.02);			
	legendFWHM->SetFillColor(0);
	legendFWHM->SetLineColor(0);
	legendFWHM->AddEntry(histoFWHMMeson,Form("#pi^{0}"));
	legendFWHM->AddEntry(histoTrueFWHMMeson,"True reconstructed #pi^{0}");
	legendFWHM->Draw();
	canvasFWHM->Update();
	DrawAliceLogoPi0Performance(pictDrawingCoordinatesMass[0], pictDrawingCoordinatesMass[1], pictDrawingCoordinatesMass[2], pictDrawingCoordinatesMass[3], pictDrawingCoordinatesMass[4], pictDrawingCoordinatesMass[5], pictDrawingCoordinatesMass[6], pictDrawingCoordinatesMass[7], pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], kMeson,1350,900,date, centralityString,kFALSE);
	
	canvasFWHM->SaveAs(Form("%s/%s_FWHM_%s_%s.%s",outputDir.Data(),prefix2.Data(),textMeson.Data(),cutSelection.Data(),suffix.Data()));
	delete canvasFWHM;
	
	
	//**********************************************************************************
	//******************** Mass_Pi0 Plot *********************************************
	//**********************************************************************************
	TCanvas* canvasMass = new TCanvas("canvasMass","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasMass, 0.08, 0.02, 0.04, 0.09);
	
	DrawGammaSetMarker(histoMassMeson, 22, 0.8, kBlack, kBlack);					 
	if (textMeson.Contains("Pi0")){      
	      DrawAutoGammaMesonHistos( histoMassMeson, 
				"", "p_{T} (GeV/c)", Form("Mass (GeV/c^{2})"), 
				kFALSE, 3.,0.,  kFALSE,
				kTRUE, 0.130,0.140, 
				kFALSE, 0., 10.);
	} else {
	  DrawAutoGammaMesonHistos( histoMassMeson, 
			  "", "p_{T} (GeV/c)", "Mass (GeV/c^{2})", 
			  kFALSE, 3.,0., kFALSE,
			  kTRUE, 0.545,0.560, 
			  kFALSE, 0., 10.);
	}

	histoMassMeson->GetYaxis()->SetNdivisions(510); 
	histoMassMeson->GetYaxis()->SetTitleOffset(1.);
	histoMassMeson->DrawCopy("e1,p"); 
	
	DrawGammaSetMarker(histoTrueMassMeson, 26, 0.8, kRed+2, kRed+2);
	histoTrueMassMeson->DrawCopy("same,e1,p"); 
		
	TLegend* legendMass = new TLegend(0.15,0.1,0.5,0.2);
	legendMass->SetTextSize(0.02);			
	legendMass->SetFillColor(0);
	legendMass->AddEntry(histoMassMeson,"Data");
	legendMass->AddEntry(histoTrueMassMeson,"MonteCarlo");
	legendMass->Draw();
	DrawAliceLogoPi0Performance(pictDrawingCoordinatesMass[0], pictDrawingCoordinatesMass[1], pictDrawingCoordinatesMass[2], pictDrawingCoordinatesMass[3], pictDrawingCoordinatesMass[4], pictDrawingCoordinatesMass[5], pictDrawingCoordinatesMass[6], pictDrawingCoordinatesMass[7], pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], kMeson,1350,900,date, centralityString,kFALSE);
	
	canvasMass->Update();
	canvasMass->SaveAs(Form("%s/%s_Mass_%s_%s.%s",outputDir.Data(),prefix2.Data(),textMeson.Data(),cutSelection.Data(),suffix.Data()));	
	delete canvasMass;
	
	nameFinalResDat = Form("%s/%s/%s_%s_FinalExtraction_%s.dat",cutSelection.Data(),optionEnergy.Data(), prefix2.Data(), textMeson.Data(), cutSelection.Data());
	fileFinalResults.open(nameFinalResDat.Data(), ios::out);
	TString applyBinShift = "";
	if (textMeson.CompareTo("Pi0")==0){
		//*****************************************************************************************************
		//**************************** Binshifted Spectra *****************************************************
		//*****************************************************************************************************
		TCanvas* canvasBinShifted = new TCanvas("canvasBinShifted","",1350,1500);  // gives the page size
		DrawGammaCanvasSettings( canvasBinShifted, 0.13, 0.02, 0.02, 0.09);
		canvasBinShifted->SetLogy();	
		
		TPad* padBinShiftedHistos = new TPad("padBinShiftedHistos", "", 0., 0.25, 1., 1.,-1, -1, -2);
		DrawGammaPadSettings( padBinShiftedHistos, 0.12, 0.02, 0.02, 0.);
		padBinShiftedHistos->Draw();
		
		TPad* padBinShiftedRatios = new TPad("padBinShiftedRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
		DrawGammaPadSettings( padBinShiftedRatios, 0.12, 0.02, 0., 0.2);
		padBinShiftedRatios->Draw();
		
		padBinShiftedHistos->cd();
		padBinShiftedHistos->SetLogy();		
		histoRecBinShift = (TH1D*) histoCorrectedYield->Clone("correctedYield");	
		DrawGammaSetMarker(histoRecBinShift, 20, 0.9, kGray+1, kGray+1);
		histoCorrYieldBinShifted= (TH1D*)histoRecBinShift->Clone();

		Double_t parametersBinShift[10];
		ReturnParameterSetFittingPbPb(cutSelection,parametersBinShift);
		for( Int_t i = 0; i < 10; i++){
			cout << "parameter " << i << "\t" << parametersBinShift[i] << endl;
      }
      
		fitBinShifting = BinShiftTH1D(histoRecBinShift, &histoCorrYieldBinShifted, textMeson.Data(), "rad", "fitBinShifting",minPtForFits,parametersBinShift);
		
		DrawGammaSetMarker(histoCorrYieldBinShifted, 24, 0.9, kBlack, kBlack);
		DrawAutoGammaMesonHistos( histoCorrYieldBinShifted, 
						"", "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}", 
						kFALSE, 3., 4e-10, kTRUE,
						kFALSE, 3e-8,10, 
						kFALSE, 0., 10.);
		DrawGammaSetMarkerTF1( fitBinShifting, 1, 0.4, kBlue-4);
		fitBinShifting->DrawCopy("same");
		histoRecBinShift->DrawCopy("same,e1"); 
			
		forOutput= WriteParameterToFile(fitBinShifting);
		fileFinalResults<< forOutput.Data()<< endl;	
				
		TLegend* legendYieldBinShifted = new TLegend(0.2,0.05,0.5,0.2);
		legendYieldBinShifted->SetFillColor(0);
		legendYieldBinShifted->SetLineColor(0);
		legendYieldBinShifted->AddEntry(histoRecBinShift,"Yield","p");
		legendYieldBinShifted->AddEntry(histoCorrYieldBinShifted,"Yield Corr (iter)","p");
		legendYieldBinShifted->AddEntry(fitBinShifting,"Fit","l");
		legendYieldBinShifted->Draw();
		
		DrawAliceLogoPi0WorkInProgress(pictDrawingCoordinates[0], pictDrawingCoordinates[1], pictDrawingCoordinates[2], pictDrawingCoordinates[3], pictDrawingCoordinates[4], pictDrawingCoordinates[5], pictDrawingCoordinates[6], pictDrawingCoordinates[7], pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], kMeson,1350,1125,kFALSE,centralityString);
		
		padBinShiftedRatios->cd();
		padBinShiftedRatios->SetLogy(0);
		
		histoRatioBinShifted = (TH1D*) histoRecBinShift->Clone();	
		histoRatioBinShifted->Divide(histoRatioBinShifted,histoCorrYieldBinShifted,1.,1.,"");
		for (Int_t i = 0; i < (histoRatioBinShifted->GetNbinsX()+1); i++){
			histoRatioBinShifted->SetBinError(i,0.);
		}
		cout << "HERE" << endl;
		DrawGammaHistoRatioLowerPanel( histoRatioBinShifted, "#frac{value}{shifted value}", 0.95, 1.42, 505, 0.08, 0.1, 0.42, 0.08, 0.11, 1.);
		DrawGammaSetMarker(histoRatioBinShifted, 20, 0.5, kBlack, kBlack);
		histoRatioBinShifted->DrawCopy("p"); 
		
		DrawGammaLines(0., maxPt,1., 1.,0.1);
		
		canvasBinShifted->Update();
		
		canvasBinShifted->SaveAs(Form("%s/%s_%s_CorrectedYieldBinShifted_%s_%s.%s",outputDir.Data(),textMeson.Data(),prefix2.Data(),makeBinShiftWithFunction.Data(),cutSelection.Data(),suffix.Data()));
		
		//***************************************************************************************************
		//*************************** Fitting  Pi0 Spectrum *************************************************
		//***************************************************************************************************
		
	// 	fileFinalResults << "Final Results for Pi0" << endl << endl;
		
		//****************************************************************************************************
		//************************** Fitting of corrected normal spectrum ************************************
		//****************************************************************************************************
		
		TCanvas* canvasFitting = new TCanvas("canvasFitting","",200,10,1350,900);  // gives the page size
		
		histoFitting = (TH1D*) histoCorrYieldBinShifted->Clone();	
		
		//****************************** Fit histCorr with Levy *****************************************
		cout << "fitting with Levy" << endl;
		fitPtLevy = FitObject("l","fitPtLevyPi0","Pi0",histoFitting,minPtForFits,maxPt,NULL,"QNRME+");
		DrawGammaSetMarkerTF1(fitPtLevy, 1, 1.5, kBlue);
		kLevySucc = kTRUE;	
		forOutput= WriteParameterToFile(fitPtLevy);
		fileFinalResults<< forOutput.Data()<< endl;
		
		//**************************** Fit histCorr with Hagedorn  **********************************
		cout << "fitting with Hagedorn" << endl;
		fitPtHagedorn = FitObject("h","fitPtHagedornPi0","Pi0",histoFitting,minPtForFits,maxPt,NULL,"QNRME+");
		DrawGammaSetMarkerTF1( fitPtHagedorn, 1, 1.5, kGreen+2);
		fitPtHagedorn->SetLineStyle(7);
		kHagSucc = kTRUE;	
		forOutput= WriteParameterToFile(fitPtHagedorn);
		fileFinalResults<< forOutput.Data()<< endl;
		
		//**************************** Fit histCorr with Radu's fuction  ****************************
		cout << "fitting with Raduslav" << endl;
		fitPtRadu = FitObject("rad","fitPtRaduPi0","Pi0",histoFitting,minPtForFits,maxPt,parametersBinShift,"QNRME+");
		DrawGammaSetMarkerTF1( fitPtRadu, 1, 1.5, kRed+2);
		fitPtRadu->SetLineStyle(7);
		kRaduSucc = kTRUE;	
		forOutput= WriteParameterToFile(fitPtRadu);
		fileFinalResults<< forOutput.Data()<< endl;
		
		//************************** Fit histCorr with Boltzmann  ***********************************

		fitPtBoltzmann = FitObject("b","fitPtBoltzmannPi0","Pi0",histoFitting,minPtForFits,2.,NULL,"QNRME+");
		DrawGammaSetMarkerTF1( fitPtBoltzmann, 1, 1.5, kMagenta+4);
		kBoltzSucc = kTRUE;	
		forOutput= WriteParameterToFile(fitPtBoltzmann);
		fileFinalResults<< forOutput.Data()<< endl;
		
		//*************************** Fit histCorr with Exponential **********************************
		fitPtExp = FitObject("e","fitPtExpPi0","Pi0",histoFitting,minPtForFits,2.,NULL,"QNRME+");
		DrawGammaSetMarkerTF1( fitPtExp, 1, 1.5, kOrange+7);
		kExpSucc = kTRUE;	
		forOutput= WriteParameterToFile(fitPtExp);
		fileFinalResults<< forOutput.Data()<< endl;
		
		//**************************** Fit histCorr with Powerlaw  *********************************
		fitPtPowerlaw = FitObject("p","fitPtPowerlawPi0","Pi0",histoFitting,1.5,maxPt,NULL,"QNRME+");
		DrawGammaSetMarkerTF1( fitPtPowerlaw, 1, 1.5, kTeal);
		kPowSucc = kTRUE;	
		forOutput= WriteParameterToFile(fitPtPowerlaw);
		fileFinalResults<< forOutput.Data()<< endl;
		
		//**************************** Fit histCorr with ModPowerlaw  *********************************
		fitPtModPowerlaw = FitObject("m","fitPtModPowerlawPi0","Pi0",histoFitting,minPtForFits,maxPt,NULL,"QNRME+");
		DrawGammaSetMarkerTF1( fitPtModPowerlaw, 1, 1.5, kMagenta+2);
		fitPtModPowerlaw->SetLineStyle(4);
		kModPowSucc = kTRUE;	
		forOutput= WriteParameterToFile(fitPtModPowerlaw);
		fileFinalResults<< forOutput.Data()<< endl;
		
		//**************************** Fit histCorr with ModPowerlaw  *********************************
		TF1* fitPtQCD = FitObject("qcd","fitPtModPowerlawPi0","Pi0",histoFitting,minPtForFits,maxPt,NULL,"QNRME+");
		DrawGammaSetMarkerTF1( fitPtQCD, 1, 1.5, kBlue-3);
		fitPtQCD->SetLineStyle(5);
		forOutput= WriteParameterToFile(fitPtQCD);
		fileFinalResults<< forOutput.Data()<< endl;
		
		//*************************** Calculating Ratios *******************************************
		histoRatioFitLevy = (TH1D*) histoCorrYieldBinShifted->Clone();	
		histoRatioFitHag = (TH1D*) histoCorrYieldBinShifted->Clone();	
		histoRatioFitRadu = (TH1D*) histoCorrYieldBinShifted->Clone();	
		histoRatioFitBoltz = (TH1D*) histoCorrYieldBinShifted->Clone();	
		histoRatioFitExp = (TH1D*) histoCorrYieldBinShifted->Clone();	
		histoRatioFitPow = (TH1D*) histoCorrYieldBinShifted->Clone();	
		histoRatioFitModPow = (TH1D*) histoCorrYieldBinShifted->Clone();	
		TH1D* histoRatioFitQCD = (TH1D*) histoCorrYieldBinShifted->Clone();	
		
		if(kLevySucc) {
			histoRatioFitLevy = CalculateHistoRatioToFit (histoRatioFitLevy, fitPtLevy); 
		}
		if(kHagSucc) {
			histoRatioFitHag = CalculateHistoRatioToFit (histoRatioFitHag, fitPtHagedorn);
		}
		if(kRaduSucc) {
			histoRatioFitRadu = CalculateHistoRatioToFit (histoRatioFitRadu, fitPtRadu);
		}
		if(kBoltzSucc) {
			histoRatioFitBoltz = CalculateHistoRatioToFit (histoRatioFitBoltz, fitPtBoltzmann);
		}
		if(kExpSucc) {
			histoRatioFitExp = CalculateHistoRatioToFit (histoRatioFitExp, fitPtExp);
		}
		if(kPowSucc) {
			histoRatioFitPow = CalculateHistoRatioToFit (histoRatioFitPow, fitPtPowerlaw);
		}	
		if(kModPowSucc) {
			histoRatioFitModPow = CalculateHistoRatioToFit (histoRatioFitModPow, fitPtModPowerlaw);
		}
		histoRatioFitQCD = CalculateHistoRatioToFit (histoRatioFitQCD, fitPtQCD);
		
		delete canvasFitting;
		//**********************************************************************************************
		//********************** Plotting of fitted spectrum *******************************************
		//**********************************************************************************************	
		TCanvas* canvasFittingSpectra = new TCanvas("canvasFittingSpectra","",1350,1500);  // gives the page size
		DrawGammaCanvasSettings( canvasFittingSpectra, 0.13, 0.02, 0.02, 0.09);
		canvasFittingSpectra->SetLogy();	
		
		TPad* padFittedSpectraHistos = new TPad("padFittedSpectraHistos", "", 0., 0.25, 1., 1.,-1, -1, -2);
		DrawGammaPadSettings( padFittedSpectraHistos, 0.12, 0.02, 0.02, 0.);
		padFittedSpectraHistos->Draw();
		
		TPad* padFittedSpectraRatios = new TPad("padFittedSpectraRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
		DrawGammaPadSettings( padFittedSpectraRatios, 0.12, 0.02, 0., 0.2);
		padFittedSpectraRatios->Draw();
		
		padFittedSpectraHistos->cd();
		padFittedSpectraHistos->SetLogy();		

		DrawGammaSetMarker(histoCorrYieldBinShifted, 22, 1., kBlack, kBlack);
		DrawAutoGammaMesonHistos( histoCorrYieldBinShifted, 
								"", "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}", 
								kTRUE, 3., 4e-10, kTRUE,
								kFALSE, 3e-8,10, 
								kFALSE, 0., 10.);		
		histoCorrectedYield->Draw("e1,x0"); 
		
		fitPtLevy->Draw("same");
		fitPtHagedorn->Draw("same");
		fitPtRadu->Draw("same");
		fitPtBoltzmann->Draw("same");
		fitPtExp->Draw("same");
		fitPtPowerlaw->Draw("same");
		fitPtModPowerlaw->Draw("same");
		fitPtQCD->Draw("same");
		
		TLegend* legendFit = new TLegend(0.15,0.02,0.55,0.15);
		legendFit->SetTextSize(0.02);			
		legendFit->SetFillColor(0);
		legendFit->AddEntry(histoCorrYieldBinShifted,"#pi^{0}");
		legendFit->AddEntry(fitPtLevy,"Levy fit");
		legendFit->AddEntry(fitPtHagedorn,"Hagedorn fit");
		legendFit->AddEntry(fitPtRadu,"Raduslav fit");
		legendFit->AddEntry(fitPtBoltzmann,"Boltzman fit");
		legendFit->AddEntry(fitPtExp,"Exponential fit");
		legendFit->AddEntry(fitPtPowerlaw,"Powerlaw fit");
		legendFit->AddEntry(fitPtModPowerlaw,"Mod Powerlaw fit");
		legendFit->AddEntry(fitPtQCD,"QCD inspired fit");
		legendFit->Draw();
		
		DrawAliceLogoPi0WorkInProgress(pictDrawingCoordinates[0], pictDrawingCoordinates[1], pictDrawingCoordinates[2], pictDrawingCoordinates[3], pictDrawingCoordinates[4], pictDrawingCoordinates[5], pictDrawingCoordinates[6], pictDrawingCoordinates[7], pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], kTRUE,1350,1125,kFALSE,centralityString);
		canvasFittingSpectra->Update();
		
		padFittedSpectraRatios->cd();
		DrawGammaHistoRatioLowerPanel( histoRatioFitLevy, "#frac{value}{fit}", 0.5, 1.55, 505, 0.08, 0.1, 0.4, 0.08, 0.11, 1.);
		
		DrawGammaSetMarker(histoRatioFitLevy, 21, 0.5, kBlue, kBlue);
		histoRatioFitLevy->Draw("e1,x0");
		
		DrawGammaSetMarker(histoRatioFitHag, 21, 0.5, kGreen+2, kGreen+2);
		if(kHagSucc)histoRatioFitHag->Draw("e1,x0,same");
		
		DrawGammaSetMarker(histoRatioFitRadu, 21, 0.5, kRed+2, kRed+2);
		if(kRaduSucc)histoRatioFitRadu->Draw("e1,x0,same");
		
		DrawGammaSetMarker(histoRatioFitBoltz, 21, 0.5, kMagenta+4, kMagenta+4);
		if(kBoltzSucc)histoRatioFitBoltz->Draw("e1,x0,same");
		
		DrawGammaSetMarker(histoRatioFitExp, 21, 0.5, kOrange+7, kOrange+7);
		if(kExpSucc)histoRatioFitExp->Draw("e1,x0,same");

		DrawGammaSetMarker(histoRatioFitPow, 21, 0.5, kTeal, kTeal);
		if(kPowSucc)histoRatioFitPow->Draw("e1,x0,same");
		
		DrawGammaSetMarker(histoRatioFitModPow, 21, 0.5, kMagenta+2, kMagenta+2);
		if(kModPowSucc)histoRatioFitModPow->Draw("e1,same");
		
		DrawGammaSetMarker(histoRatioFitQCD, 21, 0.5, kBlue-3, kBlue-3);
		histoRatioFitQCD->Draw("e1,same");
		
		DrawGammaLines(0., maxPt ,1., 1.,0.1);
		
		canvasFittingSpectra->SaveAs(Form("%s/%s_%s_CorrectedYieldFitted_%s.%s",outputDir.Data(),textMeson.Data(),prefix2.Data(),cutSelection.Data(),suffix.Data()));
		delete canvasFittingSpectra;
	}	
	
	if (textMeson.CompareTo("Eta")==0){
	  //***************************************************************************************************
	  //*************************** Fitting  Pi0 Spectrum *************************************************
	  //***************************************************************************************************
	  
	//    fileFinalResults << "Final Results for Pi0" << endl << endl;
	  
	  //****************************************************************************************************
	  //************************** Fitting of corrected normal spectrum ************************************
	  //****************************************************************************************************
	  
	  cout << "Here" << endl;
	  TCanvas* canvasFitting = new TCanvas("canvasFitting","",200,10,1350,900);  // gives the page size
	  
	  histoFitting = (TH1D*) histoCorrectedYield->Clone();   
	  
	  //****************************** Fit histCorr with Levy *****************************************
	  cout << "fitting with Levy" << endl;
	  fitPtLevy = FitObject("l","fitPtLevyPi0","Eta",histoFitting,minPtForFits,maxPt,NULL,"QNRMEI+");
	  DrawGammaSetMarkerTF1(fitPtLevy, 1, 1.5, kBlue);
	  kLevySucc = kTRUE;   
	  forOutput= WriteParameterToFile(fitPtLevy);
	  fileFinalResults<< forOutput.Data()<< endl;
	  
	  //**************************** Fit histCorr with Hagedorn  **********************************
	  cout << "fitting with Hagedorn" << endl;
	  fitPtHagedorn = FitObject("h","fitPtHagedornPi0","Eta",histoFitting,minPtForFits,maxPt,NULL,"QNRMEI+");
	  DrawGammaSetMarkerTF1( fitPtHagedorn, 1, 1.5, kGreen+2);
	  fitPtHagedorn->SetLineStyle(7);
	  kHagSucc = kTRUE; 
	  forOutput= WriteParameterToFile(fitPtHagedorn);
	  fileFinalResults<< forOutput.Data()<< endl;
	  
      /*     //**************************** Fit histCorr with Radu's fuction  ****************************
	  cout << "fitting with Raduslav" << endl;
	  fitPtRadu = FitObject("rad","fitPtRaduPi0","Eta",histoFitting,minPtForFits,maxPt,parametersBinShift,"QNRMEI+");
	  DrawGammaSetMarkerTF1( fitPtRadu, 1, 1.5, kRed+2);
	  fitPtRadu->SetLineStyle(7);
	  kRaduSucc = kTRUE;   
	  forOutput= WriteParameterToFile(fitPtRadu);
	  fileFinalResults<< forOutput.Data()<< endl;
      */     
	  //************************** Fit histCorr with Boltzmann  ***********************************

	  fitPtBoltzmann = FitObject("b","fitPtBoltzmannPi0","Eta",histoFitting,minPtForFits,2.,NULL,"QNRMEI+");
	  DrawGammaSetMarkerTF1( fitPtBoltzmann, 1, 1.5, kMagenta+4);
	  kBoltzSucc = kTRUE;  
	  forOutput= WriteParameterToFile(fitPtBoltzmann);
	  fileFinalResults<< forOutput.Data()<< endl;
	  
	  //*************************** Fit histCorr with Exponential **********************************
	  fitPtExp = FitObject("e","fitPtExpPi0","Eta",histoFitting,minPtForFits,2.,NULL,"QNRMEI+");
	  DrawGammaSetMarkerTF1( fitPtExp, 1, 1.5, kOrange+7);
	  kExpSucc = kTRUE; 
	  forOutput= WriteParameterToFile(fitPtExp);
	  fileFinalResults<< forOutput.Data()<< endl;
	  
	  //**************************** Fit histCorr with Powerlaw  *********************************
	  fitPtPowerlaw = FitObject("p","fitPtPowerlawPi0","Eta",histoFitting,1.5,maxPt,NULL,"QNRMEI+");
	  DrawGammaSetMarkerTF1( fitPtPowerlaw, 1, 1.5, kTeal);
	  kPowSucc = kTRUE; 
	  forOutput= WriteParameterToFile(fitPtPowerlaw);
	  fileFinalResults<< forOutput.Data()<< endl;
	  
	  //**************************** Fit histCorr with ModPowerlaw  *********************************
	  fitPtModPowerlaw = FitObject("m","fitPtModPowerlawPi0","Eta",histoFitting,minPtForFits,maxPt,NULL,"QNRMEI+");
	  DrawGammaSetMarkerTF1( fitPtModPowerlaw, 1, 1.5, kMagenta+2);
	  fitPtModPowerlaw->SetLineStyle(4);
	  kModPowSucc = kTRUE; 
	  forOutput= WriteParameterToFile(fitPtModPowerlaw);
	  fileFinalResults<< forOutput.Data()<< endl;
	  
      //       //**************************** Fit histCorr with ModPowerlaw  *********************************
      //       TF1* fitPtQCD = FitObject("qcd","fitPtModPowerlawPi0","Eta",histoFitting,minPtForFits,maxPt,NULL,"QNRMEI+");
      //       DrawGammaSetMarkerTF1( fitPtQCD, 1, 1.5, kBlue-3);
      //       fitPtQCD->SetLineStyle(5);
      //       forOutput= WriteParameterToFile(fitPtQCD);
      //       fileFinalResults<< forOutput.Data()<< endl;
	  
	  //*************************** Calculating Ratios *******************************************
	  histoRatioFitLevy = (TH1D*) histoCorrectedYield->Clone(); 
	  histoRatioFitHag = (TH1D*) histoCorrectedYield->Clone();  
      //       histoRatioFitRadu = (TH1D*) histoCorrectedYield->Clone(); 
	  histoRatioFitBoltz = (TH1D*) histoCorrectedYield->Clone();   
	  histoRatioFitExp = (TH1D*) histoCorrectedYield->Clone();  
	  histoRatioFitPow = (TH1D*) histoCorrectedYield->Clone();  
	  histoRatioFitModPow = (TH1D*) histoCorrectedYield->Clone();  
      //       TH1D* histoRatioFitQCD = (TH1D*) histoCorrectedYield->Clone();  
	  
	  if(kLevySucc) {
	      histoRatioFitLevy = CalculateHistoRatioToFit (histoRatioFitLevy, fitPtLevy); 
	  }
	  if(kHagSucc) {
	      histoRatioFitHag = CalculateHistoRatioToFit (histoRatioFitHag, fitPtHagedorn);
	  }
      //       if(kRaduSucc) {
      //          histoRatioFitRadu = CalculateHistoRatioToFit (histoRatioFitRadu, fitPtRadu);
      //       }
	  if(kBoltzSucc) {
	      histoRatioFitBoltz = CalculateHistoRatioToFit (histoRatioFitBoltz, fitPtBoltzmann);
	  }
	  if(kExpSucc) {
	      histoRatioFitExp = CalculateHistoRatioToFit (histoRatioFitExp, fitPtExp);
	  }
	  if(kPowSucc) {
	      histoRatioFitPow = CalculateHistoRatioToFit (histoRatioFitPow, fitPtPowerlaw);
	  }  
	  if(kModPowSucc) {
	      histoRatioFitModPow = CalculateHistoRatioToFit (histoRatioFitModPow, fitPtModPowerlaw);
	  }
      //       histoRatioFitQCD = CalculateHistoRatioToFit (histoRatioFitQCD, fitPtQCD);
	  
	  delete canvasFitting;
	  //**********************************************************************************************
	  //********************** Plotting of fitted spectrum *******************************************
	  //**********************************************************************************************   
	  TCanvas* canvasFittingSpectra = new TCanvas("canvasFittingSpectra","",1350,1500);  // gives the page size
	  DrawGammaCanvasSettings( canvasFittingSpectra, 0.13, 0.02, 0.02, 0.09);
	  canvasFittingSpectra->SetLogy(); 
	  
	  TPad* padFittedSpectraHistos = new TPad("padFittedSpectraHistos", "", 0., 0.25, 1., 1.,-1, -1, -2);
	  DrawGammaPadSettings( padFittedSpectraHistos, 0.12, 0.02, 0.02, 0.);
	  padFittedSpectraHistos->Draw();
	  
	  TPad* padFittedSpectraRatios = new TPad("padFittedSpectraRatios", "", 0., 0., 1., 0.25,-1, -1, -2);
	  DrawGammaPadSettings( padFittedSpectraRatios, 0.12, 0.02, 0., 0.2);
	  padFittedSpectraRatios->Draw();
	  
	  padFittedSpectraHistos->cd();
	  padFittedSpectraHistos->SetLogy();     

	  DrawGammaSetMarker(histoCorrectedYield, 22, 1., kBlack, kBlack);
	  DrawAutoGammaMesonHistos( histoCorrectedYield, 
			    "", "p_{T} (GeV/c)", "#frac{1}{2#pi N_{ev}} #frac{d^{2}N}{p_{T}dp_{T}dy} (c/GeV)^{2}", 
			    kTRUE, 3., 4e-10, kTRUE,
			    kFALSE, 3e-8,10, 
			    kFALSE, 0., 10.);    
	  histoCorrectedYield->Draw("e1,x0"); 
	  
	  fitPtLevy->Draw("same");
	  fitPtHagedorn->Draw("same");
      //       fitPtRadu->Draw("same");
	  fitPtBoltzmann->Draw("same");
	  fitPtExp->Draw("same");
	  fitPtPowerlaw->Draw("same");
	  fitPtModPowerlaw->Draw("same");
      //       fitPtQCD->Draw("same");
	  
	  TLegend* legendFit = new TLegend(0.15,0.02,0.55,0.15);
	  legendFit->SetTextSize(0.02);       
	  legendFit->SetFillColor(0);
	  legendFit->AddEntry(histoCorrYieldBinShifted,"#pi^{0}");
	  legendFit->AddEntry(fitPtLevy,"Levy fit");
	  legendFit->AddEntry(fitPtHagedorn,"Hagedorn fit");
      //       legendFit->AddEntry(fitPtRadu,"Raduslav fit");
	  legendFit->AddEntry(fitPtBoltzmann,"Boltzman fit");
	  legendFit->AddEntry(fitPtExp,"Exponential fit");
	  legendFit->AddEntry(fitPtPowerlaw,"Powerlaw fit");
	  legendFit->AddEntry(fitPtModPowerlaw,"Mod Powerlaw fit");
      //       legendFit->AddEntry(fitPtQCD,"QCD inspired fit");
	  legendFit->Draw();
	  
	  DrawAliceLogoPi0WorkInProgress(pictDrawingCoordinates[0], pictDrawingCoordinates[1], pictDrawingCoordinates[2], pictDrawingCoordinates[3], pictDrawingCoordinates[4], pictDrawingCoordinates[5], pictDrawingCoordinates[6], pictDrawingCoordinates[7], pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], kTRUE,1350,1125,kFALSE,centralityString);
	  canvasFittingSpectra->Update();
	  
	  padFittedSpectraRatios->cd();
	  DrawGammaHistoRatioLowerPanel( histoRatioFitLevy, "#frac{value}{fit}", 0.5, 1.55, 505, 0.08, 0.1, 0.4, 0.08, 0.11, 1.);
	  
	  DrawGammaSetMarker(histoRatioFitLevy, 21, 0.5, kBlue, kBlue);
	  histoRatioFitLevy->Draw("e1,x0");
	  
	  DrawGammaSetMarker(histoRatioFitHag, 21, 0.5, kGreen+2, kGreen+2);
	  if(kHagSucc)histoRatioFitHag->Draw("e1,x0,same");
	  
      /*     DrawGammaSetMarker(histoRatioFitRadu, 21, 0.5, kRed+2, kRed+2);
	  if(kRaduSucc)histoRatioFitRadu->Draw("e1,x0,same");
      */     
	  DrawGammaSetMarker(histoRatioFitBoltz, 21, 0.5, kMagenta+4, kMagenta+4);
	  if(kBoltzSucc)histoRatioFitBoltz->Draw("e1,x0,same");
	  
	  DrawGammaSetMarker(histoRatioFitExp, 21, 0.5, kOrange+7, kOrange+7);
	  if(kExpSucc)histoRatioFitExp->Draw("e1,x0,same");

	  DrawGammaSetMarker(histoRatioFitPow, 21, 0.5, kTeal, kTeal);
	  if(kPowSucc)histoRatioFitPow->Draw("e1,x0,same");
	  
	  DrawGammaSetMarker(histoRatioFitModPow, 21, 0.5, kMagenta+2, kMagenta+2);
	  if(kModPowSucc)histoRatioFitModPow->Draw("e1,same");
	  
      //       DrawGammaSetMarker(histoRatioFitQCD, 21, 0.5, kBlue-3, kBlue-3);
      //       histoRatioFitQCD->Draw("e1,same");
	  
	  DrawGammaLines(0., maxPt ,1., 1.,0.1);
	  
	  canvasFittingSpectra->SaveAs(Form("%s/%s_%s_CorrectedYieldFitted_%s.%s",outputDir.Data(),textMeson.Data(),prefix2.Data(),cutSelection.Data(),suffix.Data()));
	  delete canvasFittingSpectra;
	}  
   	
// 	if (textMeson.CompareTo("Eta") == 0){
// 		nPoints--;
// 	}
	//***************************** Reading systematic error and plotting Pi0 *******************
	graphCorrectedYieldSysErr      = CalculateSysErrFromRelSysHisto( histoCorrectedYield, Form("%sSystError",textMeson.Data()),relSystErrorDown , relSystErrorUp, offsetSyst, nPoints);
	graphCorrectedYieldSysErrA     = CalculateSysErrFromRelSysHisto( histoCorrectedYield, Form("%sSystErrorA",textMeson.Data()),relSystErrorWOMaterialDown , relSystErrorWOMaterialUp, offsetSyst, nPoints);
	graphCorrectedYieldStatPlusSys = CalculateSysErrFromRelSysHistoComplete( histoCorrectedYield, Form("%sComplError",textMeson.Data()),relSystErrorDown , relSystErrorUp, offsetSyst, nPoints);
	
	if (textMeson.CompareTo("Pi0")==0){
		graphCorrectedYieldSysErrBinShifted = CalculateSysErrFromRelSysHisto( histoCorrYieldBinShifted, Form("%sSystErrorBinShifted",textMeson.Data()),relSystErrorDown , relSystErrorUp, offsetSyst, nPoints);
		graphCorrectedYieldSysErrABinShifted = CalculateSysErrFromRelSysHisto( histoCorrYieldBinShifted, Form("%sSystErrorABinShifted",textMeson.Data()),relSystErrorWOMaterialDown , relSystErrorWOMaterialUp, offsetSyst, nPoints);
		graphCorrectedYieldStatPlusSysBinShifted = CalculateSysErrFromRelSysHistoComplete( histoCorrYieldBinShifted, Form("%sComplErrorBinShifted",textMeson.Data()), relSystErrorDown , relSystErrorUp, offsetSyst, nPoints);
	}
	
	
	//********************************************************************************************
	//*********************** Systematic Error ***************************************************
	//********************************************************************************************
	
	//Draw Pi0 with systematic error
	TCanvas* canvasSysErrorConversion2 = new TCanvas("canvasSysErrorConversion2","",1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasSysErrorConversion2, 0.13, 0.02, 0.02, 0.09);
	canvasSysErrorConversion2->SetLogy();
	
	TH2F * histo2DSysErrorConversion = new TH2F("histo2DSysErrorConversion","histo2DSysErrorConversion",10,0.,maxPt,1000.,4e-8,900.);
	histo2DSysErrorConversion->SetXTitle("p_{T} [GeV/c]");
	histo2DSysErrorConversion->SetYTitle("#frac{1}{2#pi N_{ev}} #frac{d^{2}N^{#pi^{0}}}{p_{T}dp_{T}dy} (c/GeV)^{2}");
	histo2DSysErrorConversion->SetTitle("");
	histo2DSysErrorConversion->GetYaxis()->SetDecimals();
	histo2DSysErrorConversion->GetYaxis()->SetLabelSize(0.03);
	histo2DSysErrorConversion->GetYaxis()->SetTitleSize(0.03);
	histo2DSysErrorConversion->GetYaxis()->SetTitleOffset(1.6);
	histo2DSysErrorConversion->GetXaxis()->SetTitleSize(0.03);
	histo2DSysErrorConversion->GetXaxis()->SetTitleOffset(1.);
	histo2DSysErrorConversion->GetXaxis()->SetNdivisions(510,kTRUE);
	histo2DSysErrorConversion->DrawCopy();  
	
	DrawGammaSetMarkerTGraph(graphCorrectedYieldSysErr, 20, 0.2, kBlack, kBlack);
	graphCorrectedYieldSysErr->SetFillColor(kGray+1);
	graphCorrectedYieldSysErr->Draw("same,2,p");
	DrawGammaSetMarker(histoCorrectedYield, 20, 0.5, kBlack, kBlack);
	histoCorrectedYield->Draw("same,e1,p");
	
	histo2DSysErrorConversion->DrawCopy("same");  
	
	DrawAliceLogoPi0WorkInProgress(pictDrawingCoordinates[0], pictDrawingCoordinates[1], pictDrawingCoordinates[2], 0.03, pictDrawingCoordinates[4], pictDrawingCoordinates[5], pictDrawingCoordinates[6], 0.04, pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], kMeson,1350,900,kFALSE,centralityString);
	
	canvasSysErrorConversion2->Update();	
	canvasSysErrorConversion2->SaveAs(Form("%s/%s_%s_CorrectedYieldwithSysErrorConversion_%s.%s",outputDir.Data(),textMeson.Data() ,prefix2.Data(),cutSelection.Data(),suffix.Data()));
	delete canvasSysErrorConversion2;

	
	TGraphAsymmErrors* graphCorrectedYieldStat;
	TGraphAsymmErrors* graphRAA;
	TGraphAsymmErrors* graphRAASys;
	TH1D* 	histoCorrectedYieldScaledNColl;
	TGraphAsymmErrors* graphCorrectedYieldSysErrBinShiftedScaledNColl ;
	
	histoNumberOfEvents =  new TH1D("histoNumberOfEvents","histoNumberOfEvents",1,0.,1.);
	histoNumberOfEvents->SetBinContent(1,nEvt);
		

	
	if (textMeson.CompareTo("Pi0") ==0){
		//****************** Pi0 corrected including systematic error with binshift ***************************************
		TCanvas* canvasSysErrorConversion3 = new TCanvas("canvasSysErrorConversion3","",1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasSysErrorConversion3, 0.13, 0.02, 0.02, 0.09);
		canvasSysErrorConversion3->SetLogy();
		
		TH2F * histo2DSysErrorConversionBinShifted = new TH2F("histo2DSysErrorConversionBinShifted","histo2DSysErrorConversionBinShifted",10,0.,maxPt,1000.,4e-8,900.);
		histo2DSysErrorConversionBinShifted->SetXTitle("p_{T} [GeV/c]");
		histo2DSysErrorConversionBinShifted->SetYTitle("#frac{1}{2#pi N_{ev}} #frac{d^{2}N^{#pi^{0}}}{p_{T}dp_{T}dy} (c/GeV)^{2}");
		histo2DSysErrorConversionBinShifted->SetTitle("");
		histo2DSysErrorConversionBinShifted->GetYaxis()->SetDecimals();
		histo2DSysErrorConversionBinShifted->GetYaxis()->SetLabelSize(0.03);
		histo2DSysErrorConversionBinShifted->GetYaxis()->SetTitleSize(0.03);
		histo2DSysErrorConversionBinShifted->GetYaxis()->SetTitleOffset(1.6);
		histo2DSysErrorConversionBinShifted->GetXaxis()->SetTitleSize(0.03);
		histo2DSysErrorConversionBinShifted->GetXaxis()->SetTitleOffset(1.);
		histo2DSysErrorConversionBinShifted->GetXaxis()->SetNdivisions(510,kTRUE);
		histo2DSysErrorConversionBinShifted->DrawCopy();  
		
		DrawGammaSetMarkerTGraph(graphCorrectedYieldSysErrBinShifted, 20, 0.2, kBlack, kBlack);
		graphCorrectedYieldSysErrBinShifted->SetFillColor(kGray+1);
		graphCorrectedYieldSysErrBinShifted->Draw("same,2,p");
		DrawGammaSetMarker(histoCorrYieldBinShifted, 20, 0.5, kBlack, kBlack);
		histoCorrYieldBinShifted->Draw("same,e1,x0,p");
		
		histo2DSysErrorConversionBinShifted->DrawCopy("same");  
				
		DrawAliceLogoPi0WorkInProgress(pictDrawingCoordinates[0], pictDrawingCoordinates[1], pictDrawingCoordinates[2], 0.03, pictDrawingCoordinates[4], pictDrawingCoordinates[5], pictDrawingCoordinates[6], 0.04, pictDrawingCoordinates[8],collisionSystem, pictDrawingOptions[1], pictDrawingOptions[2], kMeson,1350,900,kFALSE,centralityString);
		
		canvasSysErrorConversion3->Update();	
		canvasSysErrorConversion3->SaveAs(Form("%s/%s_%s_CorrectedYieldwithSysErrorConversionBinShifted_%s.%s", outputDir.Data(), textMeson.Data() ,prefix2.Data(), cutSelection.Data(), suffix.Data()));
		delete canvasSysErrorConversion3;
				
		//***************************************************************************************************
		//*************************** Fitting  Pi0 Spectrum *************************************************
		//***************************************************************************************************
		
		fileFinalResults << "Final Results for Pi0 with Systematic Errors" << endl << endl;
		
		//****************************************************************************************************
		//************************** Fitting of corrected normal spectrum ************************************
		//****************************************************************************************************
		
		TCanvas* canvasFittingSysErr = new TCanvas("canvasFittingSysErr","",200,10,1350,900);  // gives the page size
		
		graphCorrectedYieldSysErrForFit = (TGraphAsymmErrors*) graphCorrectedYieldStatPlusSysBinShifted->Clone();	
		
		//****************************** Fit histCorr with Levy *****************************************
		fitPtLevySysErr = FitObject("l","fitPtLevyPi0SysErr","Pi0",graphCorrectedYieldSysErrForFit,minPtForFits,maxPt,NULL,"QNRME+");
		DrawGammaSetMarkerTF1(fitPtLevySysErr, 1, 1.5, kBlue);
		forOutput= WriteParameterToFile(fitPtLevySysErr);
		fileFinalResults<< forOutput.Data()<< endl;
		
		//**************************** Fit histCorr with Hagedorn  **********************************
		fitPtHagedornSysErr = FitObject("h","fitPtHagedornPi0SysErr","Pi0",graphCorrectedYieldSysErrForFit,minPtForFits,maxPt,NULL,"QNRME+");
		DrawGammaSetMarkerTF1( fitPtHagedornSysErr, 1, 1.5, kGreen+2);
		fitPtHagedornSysErr->SetLineStyle(7);
		forOutput= WriteParameterToFile(fitPtHagedornSysErr);
		fileFinalResults<< forOutput.Data()<< endl;
		
		//************************** Fit histCorr with Boltzmann  ***********************************
		fitPtBoltzmannSysErr = FitObject("b","fitPtBoltzmannPi0SysErr","Pi0",graphCorrectedYieldSysErrForFit,0.3,2.,NULL,"QNRME+");
		DrawGammaSetMarkerTF1( fitPtBoltzmannSysErr, 1, 1.5, kMagenta+4);
		forOutput= WriteParameterToFile(fitPtBoltzmannSysErr);
		fileFinalResults<< forOutput.Data()<< endl;
		
		//*************************** Fit histCorr with Exponential **********************************
		fitPtExpSysErr = FitObject("e","fitPtExpPi0SysErr","Pi0",graphCorrectedYieldSysErrForFit,0.3,2.,NULL,"QNRME+");
		DrawGammaSetMarkerTF1( fitPtExpSysErr, 1, 1.5, kOrange+7);
		forOutput= WriteParameterToFile(fitPtExpSysErr);
		fileFinalResults<< forOutput.Data()<< endl;
		
		//**************************** Fit histCorr with Powerlaw  *********************************
		fitPtPowerlawSysErr = FitObject("p","fitPtPowerlawPi0SysErr","Pi0",graphCorrectedYieldSysErrForFit,1.5,maxPt,NULL,"QNRME+");
		DrawGammaSetMarkerTF1( fitPtPowerlawSysErr, 1, 1.5, kTeal);
		forOutput= WriteParameterToFile(fitPtPowerlawSysErr);
		fileFinalResults<< forOutput.Data()<< endl;
		
		//**************************** Fit histCorr with ModPowerlaw  *********************************
		fitPtModPowerlawSysErr = FitObject("m","fitPtModPowerlawPi0SysErr","Pi0",graphCorrectedYieldSysErrForFit,0.3,maxPt,NULL,"QNRME+");
		DrawGammaSetMarkerTF1( fitPtModPowerlawSysErr, 1, 1.5, kMagenta+2);
		forOutput= WriteParameterToFile(fitPtModPowerlawSysErr);
		fileFinalResults<< forOutput.Data()<< endl;
		
		delete canvasFittingSysErr;
			
		fileFinalResults.close();						

		//*********************************************************************************************************
		//********************** ComparisonFile Output ************************************************************
		//*********************************************************************************************************
		if( optNoBinShift){
			applyBinShift = "NoBinShifting";
			cout << applyBinShift << endl;
		}
		
// 		graphCorrectedYieldStat = new TGraphAsymmErrors(histoCorrYieldBinShifted);
// 		graphCorrectedYieldStat->RemovePoint(0);
// 		
// 		cout << "using nColl" << nColl << "+-" << nCollErr << endl;
// 		CalcRaa( 	graphInvSectionPCMPi02760GeV, graphInvSectionPCMSysPi02760GeV,graphInvSectionCombPi02760GeV, fitInvCrossSectionPi0Comb2760GeV,
// 						graphCorrectedYieldStat, graphCorrectedYieldSysErrABinShifted,  //PbPb Yields
// 						&graphRAA, &graphRAASys,
// 						nColl, nCollErr);
// 	// 	graphRAA->RemovePoint(graphRAA->GetN()-1);
// 	// 	graphRAASys->RemovePoint(graphRAASys->GetN()-1);
// 		histoCorrectedYieldScaledNColl = (TH1D*)histoCorrYieldBinShifted->Clone("histoCorrectedYieldScaledNColl");
// 		histoCorrectedYieldScaledNColl->Scale(1/nColl);
// 		graphCorrectedYieldSysErrBinShiftedScaledNColl = (TGraphAsymmErrors*) graphCorrectedYieldSysErrBinShifted->Clone("graphCorrectedYieldSysErrBinShiftedScaledNColl");
// 		graphCorrectedYieldSysErrBinShiftedScaledNColl= ScaleGraph(graphCorrectedYieldSysErrBinShiftedScaledNColl,1./nColl);
// 		TCanvas* canvasRAA = new TCanvas("canvasRAA","",1200,900);  // gives the page size
// 		DrawGammaCanvasSettings( canvasRAA, 0.13, 0.02, 0.02, 0.09);
// 	// 	canvasRAA->SetLogx();
// 		
// 		TH2F * histo2DRAADummy = new TH2F("histo2DRAADummy","histo2DRAADummy",1000,0.,11.,1000,0,1.5);
// 		SetStyleHistoTH2ForGraphs(histo2DRAADummy, "p_{t} (GeV/c)","RAA", 0.032,0.04, 0.04,0.04, 1,1.55);
// 		histo2DRAADummy->DrawCopy(); 
// 		
// 		DrawGammaSetMarkerTGraphAsym(graphRAA, 20,1, kRed, kRed, 0.1, kFALSE);
// 		DrawGammaSetMarkerTGraphAsym(graphRAASys,20,1, kRed, kRed, 1, kTRUE, kRed-7);
// 		graphRAASys->Draw("2same");
// 		graphRAA->Draw("p,same");
// 		
// 		canvasRAA->Update();	
// 		canvasRAA->SaveAs(Form("%s/%s_%s_RAA_%s.%s", outputDir.Data(), textMeson.Data() ,prefix2.Data(), cutSelection.Data(), suffix.Data()));
	}
	
	if (fileName2.CompareTo("")!=0 && textMeson.CompareTo("Eta")==0){
      
		if (centralityString.CompareTo("0-20%") == 0 ){
			fileNameSysErrPi0EtaBinning = "SystematicErrorsNew/SystematicErrorAveraged_Pi0EtaBinning_DummypPb.dat"; 
		} else if (centralityString.CompareTo("20-40%") == 0){
			fileNameSysErrPi0EtaBinning = "SystematicErrorsNew/SystematicErrorAveraged_Pi0EtaBinning_DummypPb.dat"; 
		}else if (centralityString.CompareTo("40-60%") == 0){
			fileNameSysErrPi0EtaBinning = "SystematicErrorsNew/SystematicErrorAveraged_Pi0EtaBinning_DummypPb.dat"; 
		}else if (centralityString.CompareTo("60-80%") == 0){
			fileNameSysErrPi0EtaBinning = "SystematicErrorsNew/SystematicErrorAveraged_Pi0EtaBinning_DummypPb.dat"; 
		} else {
			fileNameSysErrPi0EtaBinning = "SystematicErrorsNew/SystematicErrorAveraged_Pi0EtaBinning_pPb_5.023TeVMB_28_Mar_2014.dat"; 
		}

		histoCorrectedYieldEta = (TH1D*) histoCorrectedYield->Clone(); 
		
		file2 =               new TFile(fileName2);
		histoCorrectedYieldPi0EtaBinning =         (TH1D*)file2->Get("CorrectedYieldTrueEff"); 
		cout << "Bins" << endl;
		for (Int_t i = 0; i < histoCorrectedYieldPi0EtaBinning->GetNbinsX()+1; i++){
			cout << i<<"\t" <<histoCorrectedYieldPi0EtaBinning->GetXaxis()->GetBinCenter(i)<< "\t"<< histoCorrectedYieldEta->GetBinContent(i) << "\t" <<histoCorrectedYieldPi0EtaBinning->GetBinContent(i)<<  "\t" <<histoCorrectedYieldEta->GetBinContent(i)/histoCorrectedYieldPi0EtaBinning->GetBinContent(i)<< endl;
		}
		
		fileSysErr2.open(fileNameSysErrPi0EtaBinning.Data(),ios_base::in);
		cout << fileNameSysErrPi0EtaBinning.Data() << endl;
		while(!fileSysErr2.eof() && nPointsPi0EtaBinning<50){
			fileSysErr2 >> relPi0EtaBinningSystErrorDown[nPointsPi0EtaBinning] >> relPi0EtaBinningSystErrorUp[nPointsPi0EtaBinning]>> relPi0EtaBinningSystErrorWOMaterialDown[nPointsPi0EtaBinning] >> relPi0EtaBinningSystErrorWOMaterialUp[nPointsPi0EtaBinning];
			cout << nPointsPi0EtaBinning << "\t"  << relPi0EtaBinningSystErrorDown[nPointsPi0EtaBinning] << "\t"  <<relPi0EtaBinningSystErrorUp[nPointsPi0EtaBinning] << "\t" << relPi0EtaBinningSystErrorWOMaterialDown[nPointsPi0EtaBinning] << "\t"  <<relPi0EtaBinningSystErrorWOMaterialUp[nPointsPi0EtaBinning] << endl;;     
			nPointsPi0EtaBinning++;
		}
		fileSysErr2.close();
		nPointsPi0EtaBinning = nPointsPi0EtaBinning-1;
		
		
		histoRatioEtaPi0 = (TH1D*) histoCorrectedYieldEta->Clone(); 
		histoRatioEtaPi0->Divide(histoRatioEtaPi0,histoCorrectedYieldPi0EtaBinning,1.,1.,"");
			
		//***************************** Only ratio *************************************************
		TCanvas* canvasRatioEtaPi0 = new TCanvas("canvasRatioEtaPi0","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasRatioEtaPi0, 0.08, 0.01, 0.015, 0.09);
		
		
		DrawAutoGammaMesonHistos( histoRatioEtaPi0, 
				"", "p_{T} (GeV/c)", "#frac{#eta}{#pi^{0}}", 
				kFALSE, 3., 4e-10, kTRUE,
				kTRUE, 0.0,1.0,
				kFALSE, 0., 10.);
		histoRatioEtaPi0->GetYaxis()->SetLabelSize(0.03);
		histoRatioEtaPi0->GetYaxis()->SetTitleSize(0.04);
		histoRatioEtaPi0->GetYaxis()->SetDecimals();
		histoRatioEtaPi0->GetYaxis()->SetTitleOffset(0.9);
		histoRatioEtaPi0->GetXaxis()->SetTitleOffset(1.);
		histoRatioEtaPi0->GetXaxis()->SetLabelSize(0.03);
		histoRatioEtaPi0->GetXaxis()->SetTitleSize(0.04);
		histoRatioEtaPi0->GetYaxis()->SetNdivisions(510);
		DrawGammaSetMarker(histoRatioEtaPi0, 20, 1., kBlack, kBlack);
		histoRatioEtaPi0->DrawCopy("e1,x0"); 
		
		//       DrawGammaLines(0., maxPtEta ,0.45, 0.45,0.1);
		
		cout << " calculation of eta/pi0" << endl;
		for (Int_t i = 0; i < nPoints + 1; i++){
			ratioXValue[i] = histoRatioEtaPi0->GetBinCenter(i+offsetSyst);
			ratioYValue[i] = histoRatioEtaPi0->GetBinContent(i+offsetSyst);
			ratioXError[i] = histoRatioEtaPi0->GetBinWidth(i+offsetSyst)/2.;
			ratioSysUpError[i]= TMath::Sqrt(TMath::Power(relSystErrorWOMaterialUp[i]/100* histoCorrectedYieldEta->GetBinContent(i+offsetSyst)/histoCorrectedYieldPi0EtaBinning->GetBinContent(i+offsetSyst),2) + TMath::Power(histoCorrectedYieldEta->GetBinContent(i+offsetSyst)*relPi0EtaBinningSystErrorWOMaterialUp[i]/100/histoCorrectedYieldPi0EtaBinning->GetBinContent(i+offsetSyst),2));
			ratioSysDownError[i]= TMath::Sqrt(TMath::Power(relSystErrorWOMaterialDown[i]/100* histoCorrectedYieldEta->GetBinContent(i+offsetSyst)/histoCorrectedYieldPi0EtaBinning->GetBinContent(i+offsetSyst),2) + TMath::Power(histoCorrectedYieldEta->GetBinContent(i+offsetSyst)*relPi0EtaBinningSystErrorWOMaterialDown[i]/100/histoCorrectedYieldPi0EtaBinning->GetBinContent(i+offsetSyst),2));
		}
		
		graphSystErrRatio = new TGraphAsymmErrors(nPoints,ratioXValue,ratioYValue,ratioXError,ratioXError,ratioSysDownError,ratioSysUpError);
		graphSystErrRatio->SetFillColor(kGray+1);
		
		graphSystErrRatio->Draw("p,2,same");   
		histoRatioEtaPi0->DrawCopy("same,p,x0,e1"); 
		histoRatioEtaPi0->DrawCopy("same,axis"); 
		
		TLegend* legendRatio = new TLegend(0.6,0.12,0.97,0.22);
		legendRatio->SetTextSize(0.03);        
		legendRatio->SetFillColor(0);
		legendRatio->AddEntry(histoRatioEtaPi0,Form("ALICE (%s)",collisionSystem.Data()),"p");
		legendRatio->AddEntry(graphSystErrRatio,"systematic uncertainty","f");
		legendRatio->Draw();

		canvasRatioEtaPi0->Update();
		canvasRatioEtaPi0->SaveAs(Form("%s/%s_Pi0EtaRatio_%s.%s",outputDir.Data(), prefix2.Data(), cutSelection.Data(), suffix.Data()));

		delete canvasRatioEtaPi0;
		delete legendRatio;
	
   }
	
	
	const char* fileNameOutputComp = Form("%s_PCMResults_pPb%s.root",prefix2.Data(),applyBinShift.Data());
	TFile* fileOutputForComparisonFullyCorrected = new TFile(fileNameOutputComp,"UPDATE");		
		if (textMeson.CompareTo("Pi0") == 0){
			histoNumberOfEvents->Write(Form("histoNumberOfEvents%s%s",optionEnergy.Data(),centralityString.Data()),TObject::kOverwrite);
		}
		fileOutputForComparisonFullyCorrected->mkdir(Form("%s_%s_%s",textMeson.Data(),optionEnergy.Data(),centralityString.Data()));
		fileOutputForComparisonFullyCorrected->cd(Form("%s_%s_%s",textMeson.Data(), optionEnergy.Data(),centralityString.Data()));
		histoCorrectedYield->Write(Form("CorrectedYield%s",textMeson.Data()),TObject::kOverwrite);
		histoUncorrectedYield->Write(Form("RAWYieldPerEvents%s",textMeson.Data()),TObject::kOverwrite);
		histoAcceptance->Write(Form("Acceptance%s",textMeson.Data()),TObject::kOverwrite);
		histoTrueEffPt->Write(Form("Efficiency%s",textMeson.Data()),TObject::kOverwrite);
		histoMCInput->Write(Form("%sMCInput",textMeson.Data()),TObject::kOverwrite);
		if (histoMCInputWOWeighting)histoMCInputWOWeighting->Write(Form("%sMCInputWOWeights",textMeson.Data()),TObject::kOverwrite);
		if (histoMCInputWeights)histoMCInputWeights->Write(Form("%sMCInputWeights",textMeson.Data()),TObject::kOverwrite);
		if (histoMCInputAddedSig)histoMCInputAddedSig->Write(Form("%sMCInputAddedSig",textMeson.Data()),TObject::kOverwrite);
		if (histoMCInputWOWeightingAddedSig)histoMCInputWOWeightingAddedSig->Write(Form("%sMCInputWOWeightsAddedSig",textMeson.Data()),TObject::kOverwrite);
		if (histoMCInputWeightsAddedSig)histoMCInputWeightsAddedSig->Write(Form("%sMCInputWeightsAddedSig",textMeson.Data()),TObject::kOverwrite);
		
		if (textMeson.CompareTo("Pi0") == 0){
		    histoCorrYieldBinShifted->Write("CorrectedYieldPi0BinShifted",TObject::kOverwrite);
		    graphCorrectedYieldSysErrBinShifted->Write(Form("%sSystErrorBinShifted",textMeson.Data()),TObject::kOverwrite);
		    graphCorrectedYieldSysErrABinShifted->Write(Form("%sSystErrorABinShifted",textMeson.Data()),TObject::kOverwrite);
		    graphCorrectedYieldStatPlusSysBinShifted->Write(Form("%sComplErrorBinShifted",textMeson.Data()),TObject::kOverwrite);
		}
		if (fileName2.CompareTo("")!=0 && textMeson.CompareTo("Eta")==0){
		    histoRatioEtaPi0->Write("EtatoPi0Ratio",TObject::kOverwrite); 
		    graphSystErrRatio->Write("EtatoPi0RatioSys",TObject::kOverwrite);
		}
		histoFWHMMeson->Write(Form("FWHM%s",textMeson.Data()),TObject::kOverwrite);
		histoMassMeson->Write(Form("Mass%s",textMeson.Data()),TObject::kOverwrite);
		histoMassMesonMinusExp->Write(Form("Mass%sMinusExp",textMeson.Data()),TObject::kOverwrite);
		histoTrueMassMeson->Write(Form("TrueMass%s",textMeson.Data()),TObject::kOverwrite);
		histoTrueMassMesonMinusExp->Write(Form("TrueMass%sMinusExp",textMeson.Data()),TObject::kOverwrite);
		histoTrueFWHMMesonMeV->Write(Form("TrueFWHM%sMeV",textMeson.Data()),TObject::kOverwrite);
		histoFWHMMesonMeV->Write(Form("FWHM%sMeV",textMeson.Data()),TObject::kOverwrite);
		graphCorrectedYieldSysErr->Write(Form("%sSystError",textMeson.Data()),TObject::kOverwrite);
		graphCorrectedYieldSysErrA->Write(Form("%sSystErrorA",textMeson.Data()),TObject::kOverwrite);
		graphCorrectedYieldStatPlusSys->Write(Form("%sComplError",textMeson.Data()),TObject::kOverwrite);
		
		
	fileOutputForComparisonFullyCorrected->Write();
	fileOutputForComparisonFullyCorrected->Close();

}