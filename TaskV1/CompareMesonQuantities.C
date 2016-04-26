//  **********************************************************************************
//  ******     provided by Gamma Conversion Group, PWGGA,                        *****
//  ******     Friederike Bock, friederike.bock@cern.ch                          *****
//  **********************************************************************************

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
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h" 
#include "TGaxis.h"
#include "TMath.h"
#include "TMarker.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
// #include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ExtractSignalBinning.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/PlottingMeson.h"

void CompareMesonQuantities(    const char *dataFilename        = "rawSignalData", 
                                const char *mcFilename          = "rawSignalMC", 
                                TString fCutSelection           = "", 
                                TString mesonType               = "Pi0", 
                                TString fSuffix                 = "", 
                                TString energyFlag              = "" ,
                                Int_t numberOfBins              = 25,
                                Int_t mode                      = 0
                           )
{
	gROOT->Reset();
	// mode:	0 // new output PCM-PCM
	//			1 // new output PCM dalitz
	//			2 // new output PCM-Calo
	//			3 // new output Calo-Calo
    //          4 // new output EMCAL-EMCAL
    //          5 // new output PHOS-PHOS
	//			9 // old output PCM-PCM


	StyleSettingsThesis(fSuffix);	
	SetPlotStyle();
	TFile fileRawSignalData(dataFilename);
	TFile fileRawSignalMC(mcFilename);
	


	cout << dataFilename << endl;
	cout << mcFilename << endl;
	cout << fCutSelection.Data() << endl;
	cout << mesonType.Data() << endl;
	cout << fSuffix.Data()<< endl;
	cout << energyFlag.Data() << endl;
	cout << numberOfBins << endl;
	
	Double_t 	*fMesonRange = 		NULL;
	TString outputDir = 			Form("%s/%s/%s/ExtractSignal",fCutSelection.Data(),energyFlag.Data(),fSuffix.Data());
	TString nameLineShapePlot=		Form("%s/%s_MesonLineShapeCompared_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data());
	TString nameLineShapePlotLeft=		Form("%s/%s_MesonLineShapeComparedLeft_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data());
	
	TString fEventCutSelection = "";
	TString fGammaCutSelection = "";
	TString fClusterCutSelection = "";
	TString fElectronCutSelection = "";
	TString fMesonCutSelection = "";
	ReturnSeparatedCutNumberAdvanced(fCutSelection,fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);
	
	InitializeBinning(mesonType, numberOfBins, energyFlag, "", mode, fEventCutSelection, fClusterCutSelection);
	
	if (mesonType.CompareTo("Pi0") == 0 || mesonType.CompareTo("Pi0EtaBinning") == 0){
		fMesonRange 		= new Double_t[2]; 
		fMesonRange[0]		= 0.; 
		fMesonRange[1]		= 0.3;
	} else if (mesonType.CompareTo("Eta") == 0){
		fMesonRange 		= new Double_t[2]; 
		fMesonRange[0]		= 0.35; 
		fMesonRange[1]		= 0.79;
	} else if (mesonType.CompareTo("EtaPrim") == 0){
		fMesonRange 			= new Double_t[2];
		fMesonRange[0]			= 0.9; 	
		fMesonRange[1]			= 1.;;
	}

	
	cout << fStartPtBin << endl;
	
	//****************************** Specification of collision system ************************************************
	/*TString 	fTextMeasurement;
	string 	textProcess;
	if(mesonType.CompareTo("Pi0") == 0|| mesonType.CompareTo("Pi0EtaBinning") == 0){textProcess = "#pi^{0}";}
	else {textProcess = "#eta";}

	if(energyFlag.CompareTo("7TeV") == 0){
		fTextMeasurement = Form("pp #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 7 TeV ",textProcess.c_str());
	} else if(energyFlag.CompareTo("8TeV") == 0){
		fTextMeasurement = Form("pp #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 8 TeV ",textProcess.c_str());
	} else if( energyFlag.CompareTo("2.76TeV") == 0) {	
		fTextMeasurement = Form("pp #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 2.76 TeV ",textProcess.c_str());
	} else if( energyFlag.CompareTo("900GeV") == 0) {	
		fTextMeasurement = Form("pp #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 900 GeV ",textProcess.c_str());
	} else if( energyFlag.CompareTo("HI") == 0) {
		if(mesonType.CompareTo("Eta") == 0){
			cout << "No eta analysis can be carried out for Heavy Ions" << endl;
			return;
		}
		if(mesonType.CompareTo("Pi0EtaBinning") == 0){
			cout << "No Pi0 Analysis in eta binning can be carried out for Heavy Ions" << endl;
			return;
		}	
		fTextMeasurement = Form("PbPb #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 2.76 TeV ",textProcess.c_str());
	} else {
		cout << "No correct collision system specification, has been given" << endl;
		return;
		
	}*/


	//******************************* Reading histograms **************************************************************
	TH1D *  histoSignalDataInvMassPtBin[50];
	TH1D *  histoSignalMCInvMassPtBin[50];
	TH1D *  histoTrueMCInvMassPtBin[50];

	for(Int_t j=0;j<2;j++){

		TCanvas * canvasDummy = new TCanvas("canvasDummy","",2800,1800);  // gives the page size	
		canvasDummy->SetTopMargin(0.02);
		canvasDummy->SetBottomMargin(0.02);
		canvasDummy->SetRightMargin(0.02);
		canvasDummy->SetLeftMargin(0.02);
		TString histonameSignal;
		TString histonameMCTruth;
		for(Int_t iPt=fStartPtBin; iPt<fNBinsPt; iPt++){
			
			if(j==0){
				histonameSignal = Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d", iPt);
			}else{
				histonameSignal = Form("fHistoMappingSignalInvMassLeft_in_Pt_Bin%02d", iPt);
			}
			histoSignalDataInvMassPtBin[iPt]=(TH1D*)fileRawSignalData.Get(histonameSignal);
			histoSignalDataInvMassPtBin[iPt]->Scale(1./histoSignalDataInvMassPtBin[iPt]->GetMaximum());
			histoSignalMCInvMassPtBin[iPt]=(TH1D*)fileRawSignalMC.Get(histonameSignal);
			histoSignalMCInvMassPtBin[iPt]->Scale(1./histoSignalMCInvMassPtBin[iPt]->GetMaximum());
			histonameMCTruth = Form("Mapping_TrueMeson_InvMass_in_Pt_Bin%02d", iPt);
			histoTrueMCInvMassPtBin[iPt]=(TH1D*)fileRawSignalMC.Get(histonameMCTruth);
			histoTrueMCInvMassPtBin[iPt]->Scale(1./histoTrueMCInvMassPtBin[iPt]->GetMaximum());
		}	
		
		delete canvasDummy;
		
		TCanvas * canvasLineShape = new TCanvas("CanvasLineShape","",2800,1800);  // gives the page size	
		canvasLineShape->SetTopMargin(0.00);
		canvasLineShape->SetBottomMargin(0.00);
		canvasLineShape->SetRightMargin(0.00);
		canvasLineShape->SetLeftMargin(0.00);
		
		TPad * padLineShape = new TPad("PadLineShape","",0.0,0.0,1.,1.,0);   // gives the size of the histo areas 
		padLineShape->SetFillColor(0);
		padLineShape->GetFrame()->SetFillColor(0);
		padLineShape->SetBorderMode(0);
		padLineShape->SetLogy(0);
		padLineShape->Divide(fColumn,fRow,0.0,0.0);
		padLineShape->Draw();
		
		cout<<"fColumn: "<<fColumn<<" fRow: "<<fRow<<endl;
		
		Double_t relWidthLogo;
		if (mesonType.CompareTo("Pi0") == 0){
			relWidthLogo=0.5;
		} else {
			relWidthLogo=0.3;
		}
		Double_t padXWidth = 1400/fColumn; 
		Double_t padYWidth = 900/fRow;
		
		
		Int_t place = 0;
		for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
			cout<<"Pt: "<<iPt<<" of "<<fNBinsPt<<endl;
			Double_t startPt = fBinsPt[iPt];	
			Double_t endPt = fBinsPt[iPt+1];	
			
			cout << startPt << "\t" << endPt << endl;
			
			place = place + 1;						//give the right place in the page
			if(place == fColumn) {
				
				iPt--;
				padLineShape->cd(place);
				
				Double_t nPixels = 13;
				Double_t textHeight = 0.08;
				
				/*if (padLineShape->cd(place)->XtoPixel(padLineShape->cd(place)->GetX2()) < padLineShape->cd(place)->YtoPixel(padLineShape->cd(place)->GetY1())){
					textHeight = (Double_t)nPixels/padLineShape->cd(place)->XtoPixel(padLineShape->cd(place)->GetX2()) ;
				} else {
					textHeight = (Double_t)nPixels/padLineShape->cd(place)->YtoPixel(padLineShape->cd(place)->GetY1());
				}*/
				
				TString textAlice = "ALICE performance";

				TString textProcess = ReturnMesonString (mesonType);
				if(textProcess.CompareTo("") == 0 ){
					cout << "Meson unknown" << endl;
					return ;
				}
				
				TString decayChannel = Form("%s #rightarrow #gamma#gamma", textProcess.Data());
				TString energyText = ReturnFullCollisionsSystem(energyFlag);
				if (energyText.CompareTo("") == 0){
					cout << "No correct collision system specification, has been given" << endl;
					return;
				}
				TString DetectionChannel = ReturnFullTextReconstructionProcess(mode);			
				TString date=ReturnDateString();
				Double_t startTextX = 0.10;
				Double_t startTextY = 0.8;
				Double_t differenceText = textHeight*1.25;

				TLatex *alice = 		new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
				TLatex *latexDate = 	new TLatex(startTextX, (startTextY-1.25*differenceText), date.Data());
				TLatex *energy = 		new TLatex(startTextX, (startTextY-2.25*differenceText), energyText.Data());
				TLatex *process = 		new TLatex(startTextX, (startTextY-3.25*differenceText), decayChannel.Data());
				TLatex *detprocess = 	new TLatex(startTextX, (startTextY-4.25*differenceText), DetectionChannel.Data());

				alice->SetNDC();
				alice->SetTextColor(1);
				alice->SetTextSize(textHeight*1.3);
				alice->Draw();

				latexDate->SetNDC();
				latexDate->SetTextColor(1);
				latexDate->SetTextSize(textHeight);
				latexDate->Draw();
		
				energy->SetNDC();
				energy->SetTextColor(1);
				energy->SetTextSize(textHeight);
				energy->Draw();

				process->SetNDC(); 
				process->SetTextColor(1);
				process->SetTextSize(textHeight);
				process->Draw();

				detprocess->SetNDC(); 
				detprocess->SetTextColor(1);
				detprocess->SetTextSize(textHeight);
				detprocess->Draw();

				
				TLegend* legendLineShape = new TLegend(startTextX,startTextY-4.75*differenceText,1,startTextY-(4.75+2.)*differenceText);
				legendLineShape->SetTextSize(textHeight);			
				legendLineShape->SetTextFont(62);
				legendLineShape->SetFillColor(0);
				legendLineShape->SetFillStyle(0);
				legendLineShape->SetLineWidth(0);
				legendLineShape->SetLineColor(0);
				legendLineShape->SetMargin(0.15);
				Size_t markersize = histoSignalDataInvMassPtBin[fStartPtBin]->GetMarkerSize();
				histoSignalDataInvMassPtBin[fStartPtBin]->SetMarkerSize(2*markersize);
				legendLineShape->AddEntry(histoSignalDataInvMassPtBin[fStartPtBin],"Data","ep");
				Size_t markersize2 = histoSignalMCInvMassPtBin[fStartPtBin]->GetMarkerSize();
				histoSignalMCInvMassPtBin[fStartPtBin]->SetMarkerSize(2*markersize2);
				legendLineShape->AddEntry(histoSignalMCInvMassPtBin[fStartPtBin],"MC reconstructed","ep");
				Size_t linesize = histoTrueMCInvMassPtBin[fStartPtBin]->GetLineWidth();
				histoTrueMCInvMassPtBin[fStartPtBin]->SetLineWidth(linesize);
				legendLineShape->AddEntry(histoTrueMCInvMassPtBin[fStartPtBin],"MC truth" ,"l");
				legendLineShape->Draw();
				
			} else {	
				
				padLineShape->cd(place);
				padLineShape->cd(place)->SetTopMargin(0.12);
				padLineShape->cd(place)->SetBottomMargin(0.15);
				padLineShape->cd(place)->SetRightMargin(0.05);
				padLineShape->cd(place)->SetLeftMargin(0.15);
		
				if (histoTrueMCInvMassPtBin[iPt]) {
					DrawGammaHistoColored( histoTrueMCInvMassPtBin[iPt], 
							Form("%3.2f GeV/c < p_{t} < %3.2f GeV/c",startPt,endPt),
							"M_{#gamma#gamma} (GeV/c^{2})", "Counts",
							fMesonRange[0],fMesonRange[1],1,634,-1);
				}
				cout << "here" << endl;
				if (histoSignalDataInvMassPtBin[iPt]) {
					DrawGammaHistoColored( histoSignalDataInvMassPtBin[iPt], 
							Form("%3.2f GeV/c < p_{t} < %3.2f GeV/c",startPt,endPt),
							"M_{#gamma#gamma} (GeV/c^{2})", "Counts",
							fMesonRange[0],fMesonRange[1],0,1,20,0.8);
				}
				cout << "here" << endl;
				if (histoSignalMCInvMassPtBin[iPt]){ 
					DrawGammaHistoColored( histoSignalMCInvMassPtBin[iPt], 
							Form("%3.2f GeV/c < p_{t} < %3.2f GeV/c",startPt,endPt),
							"M_{#gamma#gamma} (GeV/c^{2})", "Counts",
							fMesonRange[0],fMesonRange[1],0,860,1,0.8);
				}
				cout << "here" << endl;
			}
			
		}
		cout << "saving" << endl;
		cout << nameLineShapePlot.Data() << endl;
		if(j==0) {
			canvasLineShape->SaveAs(nameLineShapePlot.Data());
		} else {
			canvasLineShape->SaveAs(nameLineShapePlotLeft.Data());
		}
		cout << "deleting" << endl;
		delete padLineShape;
		delete canvasLineShape;

	}
}



