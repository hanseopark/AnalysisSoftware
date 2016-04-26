//  **********************************************************************************
//  ******     provided by Gamma Conversion Group, PWGGA,                        *****
//  ******     Pedro Gonzalez Zamora, pedro.gonzalez.zamora@cern.c               *****
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
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "ExtractSignalDalitz.h"
#include "../CommonHeaders/ExtractSignalBinningDalitz.h"
#include "../CommonHeaders/ExtractSignalPlotting.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "THnSparse.h"

void CompareMesonQuantitiesDalitz(const char *dataFilename = "", const char *mcFilename = "", TString fCutSelection="", TString mesonType = "Pi0", TString fSuffix ="eps", TString energyFlag ="")
{

	StyleSettingsThesis();	
	SetPlotStyle();
	TFile fileRawSignalData(dataFilename);
	TFile fileRawSignalMC(mcFilename);

	cout << dataFilename << endl;
	cout << mcFilename << endl;
	cout << fCutSelection.Data() << endl;
	cout << mesonType.Data() << endl;
	cout << fSuffix.Data()<< endl;
	cout << energyFlag.Data() << endl;
	//cout << fConference.Data() << endl;
	//cout << numberOfBins << endl;
	
	Double_t 	*fMesonRange = 		NULL;
	//TString outputDir = 			Form("/home/admin1/leardini/Results/analysisCuts24Nov/cut126withAddSign/%s/%s/%s/ExtractSignalDalitz",fCutSelection.Data(),energyFlag.Data(),fSuffix.Data());
	TString outputDir = 			Form("%s/%s/%s/ExtractSignalDalitz",fCutSelection.Data(),option.Data(),fSuffix.Data());
	//TString nameLineShapePlot=		Form("%s_MesonLineShapeCompared_%s.%s",mesonType.Data(),fCutSelection.Data(),fSuffix.Data());
	TString nameLineShapePlot=		Form("%s/%s_MesonLineShapeCompared_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data());
	// **************************** Defining rows, columns, and range for plotting ********************************************
	TString fDecayChannel = "e^{+}e^{-}#gamma";
	
	
	//****************************** Specification of collision system ************************************************
	TString 	fTextMeasurement;
	

        TString textProcess = ReturnMesonString (mesonType);
	if(textProcess.CompareTo("") == 0 ){
	    cout << "Meson unknown" << endl;
	    return ;
	}
	
	fTextMeasurement = ReturnFullTextMeson(energyFlag, textProcess);
	/*fCollisionSystem = ReturnFullCollisionsSystem(fEnergyFlag);
	if (fCollisionSystem.CompareTo("") == 0){
	      cout << "No correct collision system specification, has been given" << endl;
	      return;
	}*/

	//******************************* Reading histograms **************************************************************
	TH1D *  histoSignalDataInvMassPtBin[30];
	TH1D *  histoSignalMCInvMassPtBin[30];
	TH1D *  histoTrueMCInvMassPtBin[30];
	
	//Retrieving bins///



	TArrayD *hArrayDBinsPt   	= (TArrayD*)fileRawSignalData.Get("fArrayDBinsPt");
	TArrayD *fArrayDMesonMassRange  = (TArrayD*)fileRawSignalData.Get("fArrayDMesonMassRange");
	TArrayI *hArrayIParametersBins	= (TArrayI*)fileRawSignalData.Get("fArrayIParametersBins");
	
	
	fStartPtBin = 	hArrayIParametersBins->At(0);
	fNBinsPt    =   hArrayIParametersBins->At(1);
	
	fColumn     = 	hArrayIParametersBins->At(2);
	fRow 	    = 	hArrayIParametersBins->At(3);
	
	fMesonRange = 	new Double_t[2]; 	fMesonRange[0]=	fArrayDMesonMassRange->At(0); 	fMesonRange[1]=	fArrayDMesonMassRange->At(1);
	

	TCanvas * canvasDummy = new TCanvas("canvasDummy","",2800,1800);  // gives the page size	
	canvasDummy->SetTopMargin(0.02);
	canvasDummy->SetBottomMargin(0.02);
	canvasDummy->SetRightMargin(0.02);
	canvasDummy->SetLeftMargin(0.02);
	
	TString histonameSignal;
	TString histonameMCTruth;
	for(Int_t iPt=fStartPtBin; iPt<fNBinsPt; iPt++){
		histonameSignal = Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d", iPt);
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
	canvasLineShape->SetTopMargin(0.02);
	canvasLineShape->SetBottomMargin(0.02);
	canvasLineShape->SetRightMargin(0.02);
	canvasLineShape->SetLeftMargin(0.02);
	
	TPad * padLineShape = new TPad("PadLineShape","",0.0,0.0,1.,1.,0);   // gives the size of the histo areas 
	padLineShape->SetFillColor(0);
	padLineShape->GetFrame()->SetFillColor(0);
	padLineShape->SetBorderMode(0);
	padLineShape->SetLogy(0);
	padLineShape->Divide(fColumn,fRow);
	padLineShape->Draw();
	
	cout<<"fColumn: "<<fColumn<<" fRow: "<<fRow<<endl;
	
	Double_t relWidthLogo;
	if (mesonType.CompareTo("Pi0") == 0){
		relWidthLogo=0.5;
	} else {
		relWidthLogo=0.3;
	}
	Double_t padXWidth = 2800/fColumn; 
	Double_t padYWidth = 1800/fRow;
	
	Int_t place = 0;
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		cout<<"Pt: "<<iPt<<" of "<<fNBinsPt<<endl;
		//Double_t startPt = fBinsPt[iPt];	
		//Double_t endPt = fBinsPt[iPt+1];
		
		Double_t startPt = hArrayDBinsPt->At(iPt);
		Double_t endPt   = hArrayDBinsPt->At(iPt+1);
		
		place = place + 1;						//give the right place in the page
		if (place == fColumn){
			iPt--;
			padLineShape->cd(place);
			padLineShape->cd(place)->SetTopMargin(0.12);
			padLineShape->cd(place)->SetBottomMargin(0.15);
			padLineShape->cd(place)->SetRightMargin(0.05);
			padLineShape->cd(place)->SetLeftMargin(0.15);
			
			string textAlice = "ALICE performance";
			Double_t textHeight = 0.055;
			Double_t startTextX = 0.0;
			Double_t startTextY = 0.8; //0.75;
			Double_t differenceText = textHeight*1.25;
			Double_t coordinatesStartPadX = 0.25;
			Double_t coordinatesStartPadY = 0.25;
			Double_t coordinatesEndPadX = (coordinatesStartPadX*padXWidth+relWidthLogo*padXWidth)/padXWidth;
			Double_t coordinatesEndPadY = (coordinatesStartPadY*padYWidth+relWidthLogo*padXWidth)/padYWidth;
			
	
			TLatex *alice = new TLatex(startTextX,(startTextY+(1*differenceText)),Form("%s",textAlice.c_str())); 
			TLatex *process = new TLatex(startTextX, (startTextY), fTextMeasurement);
			TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo", coordinatesStartPadX, coordinatesStartPadY, coordinatesEndPadX ,coordinatesEndPadY, -1, -1, -2);
			
			alice->SetNDC();
			alice->SetTextColor(2);
			alice->SetTextSize(textHeight);
			alice->Draw();
			
			process->SetNDC(kTRUE); 
			process->SetTextSize(textHeight);
			process->Draw();
			
			histoSignalDataInvMassPtBin[3]->SetLineColor(1);
			histoSignalMCInvMassPtBin[3]->SetLineColor(860);
			histoTrueMCInvMassPtBin[3]->SetLineColor(634);
			
			TLegend* legendLineShape = new TLegend(0.02,0.01,0.9,0.22);
			legendLineShape->SetTextSize(0.05);			
			legendLineShape->SetFillColor(0);
			legendLineShape->AddEntry(histoSignalDataInvMassPtBin[3],"Data","l");
			legendLineShape->AddEntry(histoSignalMCInvMassPtBin[3],"MC reconstructed","l");
			legendLineShape->AddEntry(histoTrueMCInvMassPtBin[3],"MC truth" ,"l");
			legendLineShape->Draw();
			
			myPadLogo->SetFillColor(0);
			myPadLogo->SetFrameFillColor(0); 
			myPadLogo->SetBorderMode(0);
			myPadLogo->SetBorderSize(0);
			myPadLogo->SetFrameBorderMode(0);
			myPadLogo->SetLeftMargin(0.0);
			myPadLogo->SetTopMargin(0.0);
			myPadLogo->SetBottomMargin(0.0);
			myPadLogo->SetRightMargin(0.0);
			
			TASImage *myAliceLogo = new TASImage("ALICE_logo.eps");
			myPadLogo->Draw();  
			myPadLogo->cd();
			myAliceLogo->Draw();	
		} else {	
			
			padLineShape->cd(place);
			padLineShape->cd(place)->SetTopMargin(0.12);
			padLineShape->cd(place)->SetBottomMargin(0.15);
			padLineShape->cd(place)->SetRightMargin(0.05);
			padLineShape->cd(place)->SetLeftMargin(0.15);
	
			DrawGammaHistoColored( histoSignalDataInvMassPtBin[iPt], 
						Form("%3.2f GeV/c < p_{t} < %3.2f GeV/c",startPt,endPt),
									Form("M_{%s} (GeV/c^{2})",fDecayChannel.Data()), Form("dN_{%s}/dM_{%s}",fDecayChannel.Data(), fDecayChannel.Data()),
									fMesonRange[0],fMesonRange[1],1,1);
									
			DrawGammaHistoColored( histoSignalMCInvMassPtBin[iPt], 
						Form("%3.2f GeV/c < p_{t} < %3.2f GeV/c",startPt,endPt),
						Form("M_{%s} (GeV/c^{2})",fDecayChannel.Data()), Form("dN_{%s}/dM_{%s}",fDecayChannel.Data(), fDecayChannel.Data()),
						fMesonRange[0],fMesonRange[1],0,860);
			
			DrawGammaHistoColored( histoTrueMCInvMassPtBin[iPt], 
								Form("%3.2f GeV/c < p_{t} < %3.2f GeV/c",startPt,endPt),
								Form("M_{%s} (GeV/c^{2})",fDecayChannel.Data()), Form("dN_{%s}/dM_{%s}",fDecayChannel.Data(), fDecayChannel.Data()),
							   fMesonRange[0],fMesonRange[1],0,634);							
		}
		
	}

	canvasLineShape->Print(nameLineShapePlot.Data());
 	delete padLineShape;
 	delete canvasLineShape;
		
}
