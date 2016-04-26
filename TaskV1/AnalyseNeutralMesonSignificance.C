/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Ana Marin, marin@physi.uni-heidelberg.de                       *****
 ******     Martin Wilde, m_wild03@uni-muenster.de                         *****
 ******     Friederike Bock, friederike.bock@cern.ch                       *****
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
#include "TTree.h"
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
#include "THnSparse.h"

TString collisionSystem;
TString date;
TString fDetectionProcess ="";

//**********************************************************************************
//******************* return minimum for 2 D histo  ********************************
//**********************************************************************************
Double_t FindSmallestEntryIn2D(TH2F* histo){
	Double_t minimum = 1;
	for (Int_t i = 1; i<histo->GetNbinsX(); i++){
		for (Int_t j = 1; j<histo->GetNbinsY(); j++){
			if (histo->GetBinContent(i,j) < minimum && histo->GetBinContent(i,j) > 0){
				minimum = histo->GetBinContent(i,j);
			}
		}
	}
	return minimum;
}

//**********************************************************************************
//******************* return minimum for 1 D histo  ********************************
//**********************************************************************************
Double_t FindSmallestEntryIn1D(TH1D* histo){
	Double_t minimum = 1;
	for (Int_t i = 1; i<histo->GetNbinsX(); i++){
		if (histo->GetBinContent(i) < minimum && histo->GetBinContent(i) > 0){
			minimum = histo->GetBinContent(i);
		}
	}
	return minimum;
}

//**********************************************************************************
//******************* return maximum for 1 D histo  ********************************
//**********************************************************************************
Double_t FindLargestEntryIn1D(TH1D* histo){
	Double_t maximum = 1;
	for (Int_t i = 1; i<histo->GetNbinsX(); i++){
		if (histo->GetBinContent(i) > maximum ){
			maximum = histo->GetBinContent(i);
		}
	}
	return maximum;
}

//**********************************************************************************
//******************* Standardized plotting of 2D plots ****************************
//**********************************************************************************
void PlotStandard2D( TH2* histo2D, TString nameOutput, TString title, TString xTitle, TString yTitle, Bool_t kRangeY, Double_t startY, Double_t endY, Bool_t kRangeX, Double_t startX, Double_t endX, Int_t logX, Int_t logY, Int_t logZ, Float_t* floatLogo, Int_t canvasSizeX = 500, Int_t canvasSizeY = 500, TString generator ="" , TString period =""){
	TCanvas * canvasStandard = new TCanvas("canvasStandard","",10,10,canvasSizeX,canvasSizeY);  // gives the page size      
	canvasStandard->SetLogx(logX);
	canvasStandard->SetLogy(logY);
	canvasStandard->SetLogz(logZ);
	canvasStandard->SetRightMargin(0.12);     
	canvasStandard->SetLeftMargin(0.12);      
	canvasStandard->SetBottomMargin(0.1);     
	canvasStandard->SetTopMargin(0.04);       
	canvasStandard->cd();
	histo2D->SetTitle("");
	DrawAutoGammaHisto2D(   histo2D,
							title.Data(), xTitle.Data(), yTitle.Data(),"",kRangeY, startY, endY, kRangeX, startX, endX);
	histo2D->GetXaxis()->SetTitleOffset(1.05);
	//    cout << histo2D->GetYaxis()->GetTitleOffset() << endl;
	histo2D->GetYaxis()->SetTitleOffset(1.35);
	if (logX==1){
	//       cout << histo2D->GetXaxis()->GetLabelOffset() << endl;
		histo2D->GetXaxis()->SetLabelOffset(0.);
	}   
		
	histo2D->Draw("colz");
	DrawLabelsEvents(floatLogo[0],floatLogo[1],floatLogo[2], 0.00, collisionSystem, generator, period);
	TLatex *detprocess = 	new TLatex(floatLogo[0], floatLogo[1] - 3.2*floatLogo[2], fDetectionProcess);
	detprocess->SetNDC(); 
	detprocess->SetTextColor(1);
	detprocess->SetTextSize(floatLogo[2]);
	detprocess->Draw();

   
   canvasStandard->Update();
   canvasStandard->SaveAs(nameOutput.Data());
   delete canvasStandard;
}

//**********************************************************************************
//************ Main function to plot addition meson distributions ******************
//**********************************************************************************
void AnalyseNeutralMesonSignificance(	TString fileNameData = "myOutput",
										TString fileNameMC = "myOutput", 
										TString cutSel = "", 
										TString suffix = "gif", 
										TString optEnergy = "", 
										TString optPeriod = "", 
										TString optGenerator = "", 
										Int_t mode = 9){

	//**********************************************************************************
	//**************************** Set global variables ********************************
	//**********************************************************************************	
	gROOT->Reset();   
	gROOT->SetStyle("Plain");
	
	StyleSettingsThesis();  
	SetPlotStyle();
	
	date 						= ReturnDateString();
	collisionSystem 			= ReturnFullCollisionsSystem(optEnergy);
	TString centralityString 	= GetCentralityString(cutSel.Data());
	fDetectionProcess 			= ReturnFullTextReconstructionProcess(mode);
	
	Float_t floatLocationRightUp2D[4] = {0.45,0.95,0.035, 0.02};
	Float_t floatLocationLeftDown2D[4] = {0.15,0.25,0.035, 0.02};
	Float_t floatLocationRightDown2D[4] = {0.45,0.25,0.035, 0.02};
	
	TString outputDirectory 	= Form("%s/%s/%s/AnalyseNeutralMesonSignificance",cutSel.Data(),optEnergy.Data(),suffix.Data());
	gSystem->Exec("mkdir -p "+outputDirectory);

	Bool_t enableExtConvCaloQA = kFALSE;
	
	//**********************************************************************************
	//**************************** Defining binning ************************************
	//**********************************************************************************
	Int_t startBinPi0 			= 0;
	if (mode == 2 )startBinPi0 	= 3;
	if (mode == 4 )startBinPi0 	= 3;

	Int_t startBinEta 			= 0;
	if (mode == 2 )startBinEta 	= 2;
	if (mode == 4 )startBinEta 	= 3;	
	
	Double_t ptBins[16] 		= {	0.2, 0.3, 0.4, 0.5, 0.6, 
									0.8, 1.0, 1.5, 2.0, 3.0, 
									4.0, 6.0, 8.0, 10.0, 15.0, 
									20.0};

	//**********************************************************************************
	//**************************** Read data file **************************************
	//**********************************************************************************
	TFile fileData(fileNameData.Data());
   	TString nameMainDir 							= "";
	if (mode == 9 || mode == 0) nameMainDir 		= "GammaConvV1";
	else if (mode == 2 || mode == 3) nameMainDir 	= "GammaConvCalo";
	else if (mode == 4 || mode == 5) nameMainDir 	= "GammaCalo";
	
	TList *TopDirData =(TList*)fileData.Get(nameMainDir.Data());
	if(TopDirData == NULL){
		cout<<"ERROR: TopDirData not Found"<<endl;
		return;
	}
	TList *HistosGammaConversionData 	= (TList*)TopDirData->FindObject(Form("Cut Number %s",cutSel.Data()));
	if(HistosGammaConversionData == NULL){
		cout<<"ERROR: " << Form("Cut Number %s",cutSel.Data()) << " not Found in File"<<endl;
		return;
	}
	TList *ESDContainerData 			= (TList*) HistosGammaConversionData->FindObject(Form("%s ESD histograms",cutSel.Data()));

	TH1D* fEventQualityData 			= (TH1D*)ESDContainerData->FindObject("NEvents");
	Double_t nEventsData 				= 0;
	if (optEnergy.CompareTo("PbPb_2.76TeV") == 0 || optEnergy.CompareTo("pPb_5.023TeV") == 0){
		nEventsData 					= fEventQualityData->GetBinContent(1);
	} else {
		nEventsData 					=  GetNEvents(fEventQualityData);
	}

	TH2F* histoPi0RapidityPtData 		= (TH2F*) ESDContainerData->FindObject("ESD_MotherPi0_Pt_Y");
	histoPi0RapidityPtData->Sumw2();
	TH1D* histoRecPi0RapidityPtDatapTBins[15];
	for (Int_t i =startBinPi0; i < 15; i++){
		histoRecPi0RapidityPtDatapTBins[i] = (TH1D*)histoPi0RapidityPtData->ProjectionY(Form("histoRecPi0RapidityPtData_%dBin",i),
																				   histoPi0RapidityPtData->GetXaxis()->FindBin(ptBins[i]),
																				   histoPi0RapidityPtData->GetXaxis()->FindBin(ptBins[i+1]));
		histoRecPi0RapidityPtDatapTBins[i]->Scale(1./nEventsData);
	}

	TH2F* histoEtaRapidityPtData 		= (TH2F*) ESDContainerData->FindObject("ESD_MotherEta_Pt_Y");
	histoEtaRapidityPtData->Sumw2();
	TH1D* histoRecEtaRapidityPtDatapTBins[15];
	for (Int_t i =startBinEta; i < 15; i++){
		histoRecEtaRapidityPtDatapTBins[i] = (TH1D*)histoEtaRapidityPtData->ProjectionY(Form("histoRecEtaRapidityPtData_%dBin",i),
																				   histoEtaRapidityPtData->GetXaxis()->FindBin(ptBins[i]),
																				   histoEtaRapidityPtData->GetXaxis()->FindBin(ptBins[i+1]));
		histoRecEtaRapidityPtDatapTBins[i]->Scale(1./nEventsData);
	}


	TH2F* histoPi0AlphaPtData 			= (TH2F*) ESDContainerData->FindObject("ESD_MotherPi0_Pt_Alpha");
	histoPi0AlphaPtData->Sumw2();
	TH2F* histoEtaAlphaPtData 			= (TH2F*) ESDContainerData->FindObject("ESD_MotherEta_Pt_Alpha");
	histoEtaAlphaPtData->Sumw2();

    TH2F* histoPi0OpenPtData           = (TH2F*) ESDContainerData->FindObject("ESD_MotherPi0_Pt_OpenAngle");
    histoPi0OpenPtData->Sumw2();
    TH2F* histoEtaOpenPtData           = (TH2F*) ESDContainerData->FindObject("ESD_MotherEta_Pt_OpenAngle");
    histoEtaOpenPtData->Sumw2();
   
	histoPi0RapidityPtData->Scale(1./nEventsData);
	histoEtaRapidityPtData->Scale(1./nEventsData);
	histoPi0AlphaPtData->Scale(1./nEventsData);
	histoEtaAlphaPtData->Scale(1./nEventsData);
    histoPi0OpenPtData->Scale(1./nEventsData);
    histoEtaOpenPtData->Scale(1./nEventsData);
	
	
	//**********************************************************************************
	//**************************** Read MC file ****************************************
	//**********************************************************************************	
	TFile* fileMC 						= new TFile(fileNameMC.Data());
	TList *TopDirMC 					= (TList*)fileMC->Get(nameMainDir.Data());
	if(TopDirMC == NULL){
		cout<<"ERROR: TopDirMC not Found"<<endl;
		return;
	}
	TList *HistosGammaConversionMC 		= (TList*)TopDirMC->FindObject(Form("Cut Number %s",cutSel.Data()));
	if(HistosGammaConversionMC == NULL){
		cout<<"ERROR: " << Form("Cut Number %s",cutSel.Data()) << " not Found in File"<<endl;
		return;
	}
	TList *ESDContainerMC 				= (TList*) HistosGammaConversionMC->FindObject(Form("%s ESD histograms",cutSel.Data()));

	TH1D* fEventQualityMC 				= (TH1D*)ESDContainerMC->FindObject("NEvents");
	Double_t nEventsMC 					= 0;
	if (optEnergy.CompareTo("PbPb_2.76TeV") == 0 || optEnergy.CompareTo("pPb_5.023TeV") == 0){
		nEventsMC 						= fEventQualityMC->GetBinContent(1);
	}else {
		nEventsMC 						=  GetNEvents(fEventQualityMC);
	}

	TH2F* histoPi0RapidityPtMC 			= (TH2F*) ESDContainerMC->FindObject("ESD_MotherPi0_Pt_Y");
	histoPi0RapidityPtMC->Sumw2();
	TH2F* histoEtaRapidityPtMC 			= (TH2F*) ESDContainerMC->FindObject("ESD_MotherEta_Pt_Y");
	histoEtaRapidityPtMC->Sumw2();
	
	TH1D* histoRecPi0RapidityPtMCpTBins[15];
	for (Int_t i =startBinPi0; i < 15; i++){
		histoRecPi0RapidityPtMCpTBins[i] = (TH1D*)histoPi0RapidityPtMC->ProjectionY(Form("histoRecPi0RapidityPtMC_%dBin",i),
																				   histoPi0RapidityPtMC->GetXaxis()->FindBin(ptBins[i]),
																				   histoPi0RapidityPtMC->GetXaxis()->FindBin(ptBins[i+1]),"e");
		histoRecPi0RapidityPtMCpTBins[i]->Scale(1./nEventsMC);
	}	
	TH1D* histoRecEtaRapidityPtMCpTBins[15];
	for (Int_t i =startBinEta; i < 15; i++){
		histoRecEtaRapidityPtMCpTBins[i] = (TH1D*)histoEtaRapidityPtMC->ProjectionY(Form("histoRecEtaRapidityPtMC_%dBin",i),
																				   histoEtaRapidityPtMC->GetXaxis()->FindBin(ptBins[i]),
																				   histoEtaRapidityPtMC->GetXaxis()->FindBin(ptBins[i+1]),"e");
		histoRecEtaRapidityPtMCpTBins[i]->Scale(1./nEventsMC);
	}

	TH2F* histoPi0AlphaPtMC 			= (TH2F*) ESDContainerMC->FindObject("ESD_MotherPi0_Pt_Alpha");
	histoPi0AlphaPtMC->Sumw2();
	TH2F* histoEtaAlphaPtMC 			= (TH2F*) ESDContainerMC->FindObject("ESD_MotherEta_Pt_Alpha");
	histoEtaAlphaPtMC->Sumw2();
    TH2F* histoPi0OpenPtMC             = (TH2F*) ESDContainerMC->FindObject("ESD_MotherPi0_Pt_OpenAngle");
    histoPi0OpenPtMC->Sumw2();
    TH2F* histoEtaOpenPtMC             = (TH2F*) ESDContainerMC->FindObject("ESD_MotherEta_Pt_OpenAngle");
    histoEtaOpenPtMC->Sumw2();

    
	TList *TrueContainerMC 				= (TList*) HistosGammaConversionMC->FindObject(Form("%s True histograms",cutSel.Data()));
	TH2F* histoTruePi0RapidityPtMC 		= (TH2F*) TrueContainerMC->FindObject("ESD_TruePi0_Pt_Y");
	histoTruePi0RapidityPtMC->Sumw2();
	TH2F* histoTrueEtaRapidityPtMC 		= (TH2F*) TrueContainerMC->FindObject("ESD_TrueEta_Pt_Y");
	histoTrueEtaRapidityPtMC->Sumw2();
	TH2F* histoTruePi0AlphaPtMC 		= (TH2F*) TrueContainerMC->FindObject("ESD_TruePi0_Pt_Alpha");
	histoTruePi0AlphaPtMC->Sumw2();
	TH2F* histoTrueEtaAlphaPtMC 		= (TH2F*) TrueContainerMC->FindObject("ESD_TrueEta_Pt_Alpha");
	histoTrueEtaAlphaPtMC->Sumw2();
    TH2F* histoTruePi0OpenPtMC         = (TH2F*) TrueContainerMC->FindObject("ESD_TruePi0_Pt_OpenAngle");
    histoTruePi0OpenPtMC->Sumw2();
    TH2F* histoTrueEtaOpenPtMC         = (TH2F*) TrueContainerMC->FindObject("ESD_TrueEta_Pt_OpenAngle");
    histoTrueEtaOpenPtMC->Sumw2();
	
	TH2F* histoTruePi0CaloConvPhotonConvRPtE 	= NULL;
	TH2F* histoTruePi0CaloConvPhotonConvRAlphaE = NULL;
	TH2F* histoTruePi0CaloConvPhotonPtEAlphaE 	= NULL;
	TH2F* histoTrueEtaCaloConvPhotonConvRPtE 	= NULL;
	TH2F* histoTrueEtaCaloConvPhotonConvRAlphaE = NULL;
	TH2F* histoTrueEtaCaloConvPhotonPtEAlphaE 	= NULL;
	
	TH1D* histoTruePi0CaloConvPhotonPtE_FullR 	= NULL;
	TH1D* histoTruePi0CaloConvPhotonPtE_RL180 	= NULL;
	TH1D* histoTruePi0CaloConvPhotonPtE_RS180 	= NULL;
	TH1D* histoTruePi0CaloConvPhotonAlphaE_FullR = NULL;
	TH1D* histoTruePi0CaloConvPhotonAlphaE_R180360 = NULL;
    TH1D* histoTruePi0CaloConvPhotonAlphaE_R360460 = NULL;
	TH1D* histoTruePi0CaloConvPhotonAlphaE_RS180 = NULL;
	TH1D* histoTruePi0CaloConvPhotonConvR		= NULL;
	TH1D* histoTruePi0CaloConvPhotonAlphaE_FullPt = NULL;
	TH1D* histoTruePi0CaloConvPhotonAlphaE_Pt0T500 = NULL;
	TH1D* histoTruePi0CaloConvPhotonAlphaE_Pt500T1 = NULL;
	TH1D* histoTruePi0CaloConvPhotonAlphaE_Pt1T2 = NULL;
	TH1D* histoTruePi0CaloConvPhotonAlphaE_PtA2 = NULL;
	
	histoTruePi0CaloConvPhotonConvRPtE 			= (TH2F*) TrueContainerMC->FindObject("ESD_TruePi0CaloConvPhoton_ConvR_PtE"); 
	if (histoTruePi0CaloConvPhotonConvRPtE) enableExtConvCaloQA = kTRUE;
	if (enableExtConvCaloQA){
		histoTruePi0CaloConvPhotonPtE_FullR 	= (TH1D*)histoTruePi0CaloConvPhotonConvRPtE->ProjectionY(	"histoTruePi0CaloConvPhotonPtE_FullR",
																											0, -1, "e");
		histoTruePi0CaloConvPhotonPtE_RL180 	= (TH1D*)histoTruePi0CaloConvPhotonConvRPtE->ProjectionY(	"histoTruePi0CaloConvPhotonPtE_RL180",
																											histoTruePi0CaloConvPhotonConvRPtE->GetXaxis()->FindBin(180),
																											histoTruePi0CaloConvPhotonConvRPtE->GetXaxis()->FindBin(460),
																											"e");
		histoTruePi0CaloConvPhotonPtE_RS180 	= (TH1D*)histoTruePi0CaloConvPhotonConvRPtE->ProjectionY(	"histoTruePi0CaloConvPhotonPtE_RS180",
																											histoTruePi0CaloConvPhotonConvRPtE->GetXaxis()->FindBin(0.),
																											histoTruePi0CaloConvPhotonConvRPtE->GetXaxis()->FindBin(180),
																											"e");
		histoTruePi0CaloConvPhotonConvR			= (TH1D*)histoTruePi0CaloConvPhotonConvRPtE->ProjectionX(	"histoTruePi0CaloConvPhotonConvR",
																											0, -1, "e");
		histoTruePi0CaloConvPhotonConvR->Rebin(2);
		histoTruePi0CaloConvPhotonConvRAlphaE 	= (TH2F*) TrueContainerMC->FindObject("ESD_TruePi0CaloConvPhoton_ConvR_AlphaE"); 
		histoTruePi0CaloConvPhotonAlphaE_FullR 	= (TH1D*)histoTruePi0CaloConvPhotonConvRAlphaE->ProjectionY("histoTruePi0CaloConvPhotonAlphaE_FullR",
																											0, -1, "e");
		histoTruePi0CaloConvPhotonAlphaE_R180360 	= (TH1D*)histoTruePi0CaloConvPhotonConvRAlphaE->ProjectionY("histoTruePi0CaloConvPhotonAlphaE_R180360",
																											histoTruePi0CaloConvPhotonConvRAlphaE->GetXaxis()->FindBin(180),
																											histoTruePi0CaloConvPhotonConvRAlphaE->GetXaxis()->FindBin(360),
																											"e");
        histoTruePi0CaloConvPhotonAlphaE_R360460    = (TH1D*)histoTruePi0CaloConvPhotonConvRAlphaE->ProjectionY("histoTruePi0CaloConvPhotonAlphaE_R360460",
                                                                                                            histoTruePi0CaloConvPhotonConvRAlphaE->GetXaxis()->FindBin(360),
                                                                                                            histoTruePi0CaloConvPhotonConvRAlphaE->GetXaxis()->FindBin(460),
                                                                                                            "e");
        histoTruePi0CaloConvPhotonAlphaE_RS180 	= (TH1D*)histoTruePi0CaloConvPhotonConvRAlphaE->ProjectionY("histoTruePi0CaloConvPhotonAlphaE_RS180",
																											histoTruePi0CaloConvPhotonConvRAlphaE->GetXaxis()->FindBin(0.),
																											histoTruePi0CaloConvPhotonConvRAlphaE->GetXaxis()->FindBin(180),
																											"e");
		
		histoTruePi0CaloConvPhotonPtEAlphaE 	= (TH2F*) TrueContainerMC->FindObject("ESD_TruePi0CaloConvPhoton_PtE_AlphaE"); 
		histoTruePi0CaloConvPhotonAlphaE_FullPt 	= (TH1D*)histoTruePi0CaloConvPhotonPtEAlphaE->ProjectionY("histoTruePi0CaloConvPhotonAlphaE_FullPt",
																											0, -1, "e");
		histoTruePi0CaloConvPhotonAlphaE_Pt0T500 	= (TH1D*)histoTruePi0CaloConvPhotonPtEAlphaE->ProjectionY("histoTruePi0CaloConvPhotonAlphaE_Pt0T500",
																											histoTruePi0CaloConvPhotonPtEAlphaE->GetXaxis()->FindBin(0.),
																											histoTruePi0CaloConvPhotonPtEAlphaE->GetXaxis()->FindBin(0.5),
																											"e");
		histoTruePi0CaloConvPhotonAlphaE_Pt500T1 	= (TH1D*)histoTruePi0CaloConvPhotonPtEAlphaE->ProjectionY("histoTruePi0CaloConvPhotonAlphaE_Pt500T1",
																											histoTruePi0CaloConvPhotonPtEAlphaE->GetXaxis()->FindBin(0.5),
																											histoTruePi0CaloConvPhotonPtEAlphaE->GetXaxis()->FindBin(1.),
																											"e");
		histoTruePi0CaloConvPhotonAlphaE_Pt1T2 		= (TH1D*)histoTruePi0CaloConvPhotonPtEAlphaE->ProjectionY("histoTruePi0CaloConvPhotonAlphaE_Pt1T2",
																											histoTruePi0CaloConvPhotonPtEAlphaE->GetXaxis()->FindBin(1.),
																											histoTruePi0CaloConvPhotonPtEAlphaE->GetXaxis()->FindBin(2.),
																											"e");
		
		histoTruePi0CaloConvPhotonAlphaE_PtA2 		= (TH1D*)histoTruePi0CaloConvPhotonPtEAlphaE->ProjectionY("histoTruePi0CaloConvPhotonAlphaE_PtA2",
																											histoTruePi0CaloConvPhotonPtEAlphaE->GetXaxis()->FindBin(2.),
																											histoTruePi0CaloConvPhotonPtEAlphaE->GetXaxis()->FindBin(35.),
																											"e");
		
		histoTrueEtaCaloConvPhotonConvRPtE 		= (TH2F*) TrueContainerMC->FindObject("ESD_TrueEtaCaloConvPhoton_ConvR_PtE"); 
		histoTrueEtaCaloConvPhotonConvRAlphaE 	= (TH2F*) TrueContainerMC->FindObject("ESD_TrueEtaCaloConvPhoton_ConvR_AlphaE"); 
		histoTrueEtaCaloConvPhotonPtEAlphaE 	= (TH2F*) TrueContainerMC->FindObject("ESD_TrueEtaCaloConvPhoton_PtE_AlphaE"); 
	}
	
	TH1D* histoTruePi0RapidityPtMCpTBins[15];
	for (Int_t i =startBinPi0; i < 15; i++){
		histoTruePi0RapidityPtMCpTBins[i] = (TH1D*)histoTruePi0RapidityPtMC->ProjectionY(Form("histoTruePi0RapidityPtMC_%dBin",i),
																				   histoTruePi0RapidityPtMC->GetXaxis()->FindBin(ptBins[i]),
																				   histoTruePi0RapidityPtMC->GetXaxis()->FindBin(ptBins[i+1]), "e");
		histoTruePi0RapidityPtMCpTBins[i]->Scale(1./nEventsMC);
	}
							
	TH1D* histoTruePi0AlphaPtMCpTBins[15];
	for (Int_t i =startBinPi0; i < 15; i++){
		histoTruePi0AlphaPtMCpTBins[i] = (TH1D*)histoTruePi0AlphaPtMC->ProjectionY(Form("histoTruePi0AlphaPtMC_%dBin",i),
																				   histoTruePi0AlphaPtMC->GetXaxis()->FindBin(ptBins[i]),
																				   histoTruePi0AlphaPtMC->GetXaxis()->FindBin(ptBins[i+1]),"e");
		histoTruePi0AlphaPtMCpTBins[i]->Scale(1./nEventsMC);
	}

	TH1D* histoTrueEtaRapidityPtMCpTBins[15];
	for (Int_t i =startBinEta; i < 15; i++){
		histoTrueEtaRapidityPtMCpTBins[i] = (TH1D*)histoTrueEtaRapidityPtMC->ProjectionY(Form("histoTrueEtaRapidityPtMC_%dBin",i),
																				   histoTrueEtaRapidityPtMC->GetXaxis()->FindBin(ptBins[i]),
																				   histoTrueEtaRapidityPtMC->GetXaxis()->FindBin(ptBins[i+1]),"e");
		histoTrueEtaRapidityPtMCpTBins[i]->Scale(1./nEventsMC);
	}

	TH1D* histoTrueEtaAlphaPtMCpTBins[15];
	for (Int_t i =startBinEta; i < 15; i++){
		histoTrueEtaAlphaPtMCpTBins[i] = (TH1D*)histoTrueEtaAlphaPtMC->ProjectionY(Form("histoTrueEtaAlphaPtMC_%dBin",i),
																				   histoTrueEtaAlphaPtMC->GetXaxis()->FindBin(ptBins[i]),
																				   histoTrueEtaAlphaPtMC->GetXaxis()->FindBin(ptBins[i+1]),"e");
		histoTrueEtaAlphaPtMCpTBins[i]->Scale(1./nEventsMC);
	}

		
	TH2F* histoBGPi0RapidityPtMC = (TH2F*) histoPi0RapidityPtMC->Clone("histoBGPi0RapidityPtMC");
	histoBGPi0RapidityPtMC->Sumw2();
	histoBGPi0RapidityPtMC->Add(histoTruePi0RapidityPtMC,-1);
	TH2F* histoBGEtaRapidityPtMC = (TH2F*) histoEtaRapidityPtMC->Clone("histoBGEtaRapidityPtMC");
	histoBGEtaRapidityPtMC->Sumw2();
	histoBGEtaRapidityPtMC->Add(histoTrueEtaRapidityPtMC,-1);
	TH2F* histoBGPi0AlphaPtMC = (TH2F*) histoPi0AlphaPtMC->Clone("histoBGPi0AlphaPtMC");
	histoBGPi0AlphaPtMC->Sumw2();
	histoBGPi0AlphaPtMC->Add(histoTruePi0AlphaPtMC,-1);
    TH2F* histoBGPi0OpenPtMC = (TH2F*) histoPi0OpenPtMC->Clone("histoBGPi0OpenPtMC");
    histoBGPi0OpenPtMC->Sumw2();
    histoBGPi0OpenPtMC->Add(histoTruePi0OpenPtMC,-1);
	
	TH1D* histoBGPi0AlphaPtMCCpTBins[15];
	for (Int_t i =startBinPi0; i < 15; i++){
		histoBGPi0AlphaPtMCCpTBins[i] = (TH1D*)histoBGPi0AlphaPtMC->ProjectionY(Form("histoBGPi0AlphaPtMC_%dBin",i),
																				   histoBGPi0AlphaPtMC->GetXaxis()->FindBin(ptBins[i]),
																				   histoBGPi0AlphaPtMC->GetXaxis()->FindBin(ptBins[i+1]),"e");
		histoBGPi0AlphaPtMCCpTBins[i]->Scale(1./nEventsMC);
	}
	
	TH2F* histoBGEtaAlphaPtMC = (TH2F*) histoEtaAlphaPtMC->Clone("histoBGEtaAlphaPtMC");
	histoBGEtaAlphaPtMC->Sumw2();
	histoBGEtaAlphaPtMC->Add(histoTrueEtaAlphaPtMC,-1);
    TH2F* histoBGEtaOpenPtMC = (TH2F*) histoEtaOpenPtMC->Clone("histoBGEtaOpenPtMC");
    histoBGEtaOpenPtMC->Sumw2();
    histoBGEtaOpenPtMC->Add(histoTrueEtaOpenPtMC,-1);

    TH1D* histoBGEtaAlphaPtMCCpTBins[15];
	for (Int_t i =startBinEta; i < 15; i++){
		histoBGEtaAlphaPtMCCpTBins[i] = (TH1D*)histoBGEtaAlphaPtMC->ProjectionY(Form("histoBGEtaAlphaPtMC_%dBin",i),
																				   histoBGEtaAlphaPtMC->GetXaxis()->FindBin(ptBins[i]),
																				   histoBGEtaAlphaPtMC->GetXaxis()->FindBin(ptBins[i+1]),"e");
		histoBGEtaAlphaPtMCCpTBins[i]->Scale(1./nEventsMC);
	}

	
	histoPi0RapidityPtMC->Scale(1./nEventsMC);
	histoEtaRapidityPtMC->Scale(1./nEventsMC);
	histoPi0AlphaPtMC->Scale(1./nEventsMC);
	histoEtaAlphaPtMC->Scale(1./nEventsMC);
    histoPi0OpenPtMC->Scale(1./nEventsMC);
    histoEtaOpenPtMC->Scale(1./nEventsMC);

    histoTruePi0RapidityPtMC->Scale(1./nEventsMC);
	histoTrueEtaRapidityPtMC->Scale(1./nEventsMC);
	histoTruePi0AlphaPtMC->Scale(1./nEventsMC);
	histoTrueEtaAlphaPtMC->Scale(1./nEventsMC);
    histoTruePi0OpenPtMC->Scale(1./nEventsMC);
    histoTrueEtaOpenPtMC->Scale(1./nEventsMC);

    histoBGPi0RapidityPtMC->Scale(1./nEventsMC);
	histoBGEtaRapidityPtMC->Scale(1./nEventsMC);
	histoBGPi0AlphaPtMC->Scale(1./nEventsMC);
	histoBGEtaAlphaPtMC->Scale(1./nEventsMC);
    histoBGPi0OpenPtMC->Scale(1./nEventsMC);
    histoBGEtaOpenPtMC->Scale(1./nEventsMC);
	
	
	//**********************************************************************************
	//********************* Define minima and maxima for all 2D plots ******************
	//**********************************************************************************	
	Double_t maximumEtaY =0;
	if (histoEtaRapidityPtData->GetMaximum() > histoEtaRapidityPtMC->GetMaximum()) maximumEtaY = histoEtaRapidityPtData->GetMaximum();
		else maximumEtaY = histoEtaRapidityPtMC->GetMaximum();
	Double_t maximumPi0Y =0;
	if (histoPi0RapidityPtData->GetMaximum() > histoPi0RapidityPtMC->GetMaximum()) maximumPi0Y = histoPi0RapidityPtData->GetMaximum();
		else maximumPi0Y = histoPi0RapidityPtMC->GetMaximum();
	Double_t maximumEtaAlpha =0;
	if (histoEtaAlphaPtData->GetMaximum() > histoEtaAlphaPtMC->GetMaximum()) maximumEtaAlpha = histoEtaAlphaPtData->GetMaximum();
		else maximumEtaAlpha = histoEtaAlphaPtMC->GetMaximum();
	Double_t maximumPi0Alpha =0;
	if (histoPi0AlphaPtData->GetMaximum() > histoPi0AlphaPtMC->GetMaximum()) maximumPi0Alpha = histoPi0AlphaPtData->GetMaximum();
		else maximumPi0Alpha = histoPi0AlphaPtMC->GetMaximum();
    Double_t maximumEtaOpen =0;
    if (histoEtaOpenPtData->GetMaximum() > histoEtaOpenPtMC->GetMaximum()) maximumEtaOpen = histoEtaOpenPtData->GetMaximum();
        else maximumEtaOpen = histoEtaOpenPtMC->GetMaximum();
    Double_t maximumPi0Open =0;
    if (histoPi0OpenPtData->GetMaximum() > histoPi0OpenPtMC->GetMaximum()) maximumPi0Open = histoPi0OpenPtData->GetMaximum();
        else maximumPi0Open = histoPi0OpenPtMC->GetMaximum();
	
	Double_t minimum = 0;
	Double_t minimumData = FindSmallestEntryIn2D(histoPi0RapidityPtData);
	Double_t minimumMC = FindSmallestEntryIn2D(histoPi0RapidityPtMC);
	if (minimumData > minimumMC ){
		if (optGenerator.Contains("JetJet")){
			minimum = minimumMC*1000;
		} else {
			minimum = minimumMC*10;
		}	
	} else minimum = minimumData*10;
	cout << minimumData << "\t" << minimumMC << "\t" << minimum << "\t"<< 1/fEventQualityMC->GetEntries()<< endl;
		
	
	//**********************************************************************************
	//**************************** Plot 2D distributions *******************************
	//**********************************************************************************	
	histoPi0RapidityPtData->GetZaxis()->SetRangeUser(minimum, maximumPi0Y);
	PlotStandard2D( histoPi0RapidityPtData , 
					Form("%s/RecPi0_Y_Pt_Data.%s",outputDirectory.Data(),suffix.Data()), 
					"", "#it{p}_{#pi^{0},T} (GeV/#it{c})", "#it{y}",  
					kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1, 
					floatLocationRightUp2D,500,500,"Data", optPeriod);
	histoEtaRapidityPtData->GetZaxis()->SetRangeUser(minimum, maximumEtaY);
	PlotStandard2D( histoEtaRapidityPtData ,
					Form("%s/RecEta_Y_Pt_Data.%s",outputDirectory.Data(),suffix.Data()),
					"", "#it{p}_{#eta,T} (GeV/#it{c})", "#it{y}",  
					kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1, 
					floatLocationRightUp2D,500,500,"Data", optPeriod);
	histoPi0AlphaPtData->GetZaxis()->SetRangeUser(minimum, maximumPi0Alpha);
	PlotStandard2D( histoPi0AlphaPtData ,
					Form("%s/RecPi0_Alpha_Pt_Data.%s",outputDirectory.Data(),suffix.Data()), 
					"", "#it{p}_{#pi^{0},T} (GeV/#it{c})", "#it{#alpha}",  
					kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1, 
					floatLocationRightUp2D,500,500,"Data", optPeriod);
	histoEtaAlphaPtData->GetZaxis()->SetRangeUser(minimum, maximumEtaAlpha);
	PlotStandard2D( histoEtaAlphaPtData , 
					Form("%s/RecEta_Alpha_Pt_Data.%s",outputDirectory.Data(),suffix.Data()), 
					"", "#it{p}_{#eta,T} (GeV/#it{c})", "#it{#alpha}",  
					kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1, 
					floatLocationRightUp2D,500,500,"Data", optPeriod);
    histoPi0OpenPtData->GetZaxis()->SetRangeUser(minimum, maximumPi0Open);
    PlotStandard2D( histoPi0OpenPtData ,
                    Form("%s/RecPi0_Open_Pt_Data.%s",outputDirectory.Data(),suffix.Data()), 
                    "", "#it{p}_{#pi^{0},T} (GeV/#it{c})", "#it{#theta}_{#pi^{0}, cand}",  
                    kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1, 
                    floatLocationRightUp2D,500,500,"Data", optPeriod);
    histoEtaOpenPtData->GetZaxis()->SetRangeUser(minimum, maximumEtaOpen);
    PlotStandard2D( histoEtaOpenPtData , 
                    Form("%s/RecEta_Open_Pt_Data.%s",outputDirectory.Data(),suffix.Data()), 
                    "", "#it{p}_{#eta,T} (GeV/#it{c})", "#it{#theta}_{#eta, cand}",  
                    kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1, 
                    floatLocationRightUp2D,500,500,"Data", optPeriod);

	histoPi0RapidityPtMC->GetZaxis()->SetRangeUser(minimum, maximumPi0Y);
	PlotStandard2D( histoPi0RapidityPtMC , 
					Form("%s/RecPi0_Y_Pt_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()), 
					"", "#it{p}_{#pi^{0},T} (GeV/#it{c})", "#it{y}",  
					kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1, 
					floatLocationRightUp2D,500,500,optGenerator, optPeriod);
	histoEtaRapidityPtMC->GetZaxis()->SetRangeUser(minimum, maximumEtaY);
	PlotStandard2D( histoEtaRapidityPtMC , 
					Form("%s/RecEta_Y_Pt_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()), 
					"", "#it{p}_{#eta,T} (GeV/#it{c})", "#it{y}",   
					kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1, 
					floatLocationRightUp2D,500,500,optGenerator, optPeriod);
	histoPi0AlphaPtMC->GetZaxis()->SetRangeUser(minimum, maximumPi0Alpha);
	PlotStandard2D( histoPi0AlphaPtMC , 
					Form("%s/RecPi0_Alpha_Pt_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()), 
					"", "#it{p}_{#pi^{0},T} (GeV/#it{c})", "#it{#alpha}",   
					kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1, 
					floatLocationRightUp2D,500,500,optGenerator, optPeriod);
	histoEtaAlphaPtMC->GetZaxis()->SetRangeUser(minimum, maximumEtaAlpha);
	PlotStandard2D( histoEtaAlphaPtMC , 
					Form("%s/RecEta_Alpha_Pt_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()), 
					"", "#it{p}_{#eta,T} (GeV/#it{c})", "#it{#alpha}",   
					kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1, 
					floatLocationRightUp2D,500,500,optGenerator, optPeriod);
    histoPi0OpenPtMC->GetZaxis()->SetRangeUser(minimum, maximumPi0Open);
    PlotStandard2D( histoPi0OpenPtMC , 
                    Form("%s/RecPi0_Open_Pt_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()), 
                    "", "#it{p}_{#pi^{0},T} (GeV/#it{c})", "#it{#theta}_{#pi^{0}, cand}",  
                    kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1, 
                    floatLocationRightUp2D,500,500,optGenerator, optPeriod);
    histoEtaOpenPtMC->GetZaxis()->SetRangeUser(minimum, maximumEtaOpen);
    PlotStandard2D( histoEtaOpenPtMC , 
                    Form("%s/RecEta_Open_Pt_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()), 
                    "", "#it{p}_{#eta,T} (GeV/#it{c})", "#it{#theta}_{#eta, cand}",  
                    kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1, 
                    floatLocationRightUp2D,500,500,optGenerator, optPeriod);
	
	histoTruePi0RapidityPtMC->GetZaxis()->SetRangeUser(minimum, maximumPi0Y);
	PlotStandard2D( histoTruePi0RapidityPtMC ,
					Form("%s/TruePi0_Y_Pt_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()), 
					"", "#it{p}_{#pi^{0},T} (GeV/#it{c})", "#it{y}",   
					kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1, 
					floatLocationRightUp2D,500,500,optGenerator, optPeriod);
	histoTrueEtaRapidityPtMC->GetZaxis()->SetRangeUser(minimum, maximumEtaY);
	PlotStandard2D( histoTrueEtaRapidityPtMC , 
					Form("%s/TrueEta_Y_Pt_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()), 
					"", "#it{p}_{#eta,T} (GeV/#it{c})", "#it{y}",  
					kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1,
					floatLocationRightUp2D,500,500,optGenerator, optPeriod);
	histoTruePi0AlphaPtMC->GetZaxis()->SetRangeUser(minimum, maximumPi0Alpha);
	PlotStandard2D( histoTruePi0AlphaPtMC , 
					Form("%s/TruePi0_Alpha_Pt_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()),
					"", "#it{p}_{#pi^{0},T} (GeV/#it{c})", "#it{#alpha}",   
					kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1, floatLocationRightUp2D,500,500,optGenerator, optPeriod);
	histoTrueEtaAlphaPtMC->GetZaxis()->SetRangeUser(minimum, maximumEtaAlpha);
	PlotStandard2D( histoTrueEtaAlphaPtMC , 
					Form("%s/TrueEta_Alpha_Pt_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()), 
					"", "#it{p}_{#eta,T} (GeV/#it{c})", "#it{#alpha}",  
					kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1, 
					floatLocationRightUp2D,500,500,optGenerator, optPeriod);
    histoTruePi0OpenPtMC->GetZaxis()->SetRangeUser(minimum, maximumPi0Open);
    PlotStandard2D( histoTruePi0OpenPtMC , 
                    Form("%s/TruePi0_Open_Pt_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()),
                    "", "#it{p}_{#pi^{0},T} (GeV/#it{c})", "#it{#theta}_{#pi^{0}}",  
                    kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1, floatLocationRightUp2D,500,500,optGenerator, optPeriod);
    histoTrueEtaOpenPtMC->GetZaxis()->SetRangeUser(minimum, maximumEtaOpen);
    PlotStandard2D( histoTrueEtaOpenPtMC , 
                    Form("%s/TrueEta_Open_Pt_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()), 
                    "", "#it{p}_{#eta,T} (GeV/#it{c})", "#it{#theta}_{#eta}",  
                    kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1, 
                    floatLocationRightUp2D,500,500,optGenerator, optPeriod);
	
	histoBGPi0RapidityPtMC->GetZaxis()->SetRangeUser(minimum, maximumPi0Y);
	PlotStandard2D( histoBGPi0RapidityPtMC , 
					Form("%s/BGPi0_Y_Pt_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()), 
					"", "#it{p}_{#pi^{0},T} (GeV/#it{c})", "#it{y}",  
					kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1, 
					floatLocationRightUp2D,500,500,optGenerator, optPeriod);
	histoBGEtaRapidityPtMC->GetZaxis()->SetRangeUser(minimum, maximumEtaY);
	PlotStandard2D( histoBGEtaRapidityPtMC , 
					Form("%s/BGEta_Y_Pt_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()), 
					"", "#it{p}_{#eta,T} (GeV/#it{c})", "#it{y}",  
					kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1, 
					floatLocationRightUp2D,500,500,optGenerator, optPeriod);
	histoBGPi0AlphaPtMC->GetZaxis()->SetRangeUser(minimum, maximumPi0Alpha);
	PlotStandard2D( histoBGPi0AlphaPtMC , 
					Form("%s/BGPi0_Alpha_Pt_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()), 
					"", "#it{p}_{#pi^{0},T} (GeV/#it{c})", "#it{#theta}_{BG #pi^{0} mass window}",  
					kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1, 
					floatLocationRightUp2D,500,500,optGenerator, optPeriod);
	histoBGEtaAlphaPtMC->GetZaxis()->SetRangeUser(minimum, maximumEtaAlpha);
	PlotStandard2D( histoBGEtaAlphaPtMC , 
					Form("%s/BGEta_Alpha_Pt_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()), 
					"", "#it{p}_{#eta,T} (GeV/#it{c})", "#it{#theta}_{BG #eta mass window}",  
					kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1, 
					floatLocationRightUp2D,500,500,optGenerator, optPeriod);
    histoBGPi0OpenPtMC->GetZaxis()->SetRangeUser(minimum, maximumPi0Open);
    PlotStandard2D( histoBGPi0OpenPtMC , 
                    Form("%s/BGPi0_Open_Pt_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()), 
                    "", "#it{p}_{#pi^{0},T} (GeV/#it{c})", "#it{#alpha}",   
                    kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1, 
                    floatLocationRightUp2D,500,500,optGenerator, optPeriod);
    histoBGEtaOpenPtMC->GetZaxis()->SetRangeUser(minimum, maximumEtaOpen);
    PlotStandard2D( histoBGEtaOpenPtMC , 
                    Form("%s/BGEta_Open_Pt_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()), 
                    "", "#it{p}_{#eta,T} (GeV/#it{c})", "#it{#alpha}",   
                    kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1, 
                    floatLocationRightUp2D,500,500,optGenerator, optPeriod);

	if (enableExtConvCaloQA){
		cout << "line " << __LINE__ << endl;
		histoTruePi0CaloConvPhotonConvRPtE->GetZaxis()->SetRangeUser(FindSmallestEntryIn2D(histoTruePi0CaloConvPhotonConvRPtE), 
																	 histoTruePi0CaloConvPhotonConvRPtE->GetMaximum());
		PlotStandard2D( histoTruePi0CaloConvPhotonConvRPtE ,
						Form("%s/TruePi0CaloConvPhoton_ConvR_PtElectron_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()), 
						"", "#it{R}_{conv, cluster}", "#it{p}_{T, e clus}",   
						kTRUE, 0.0, 10., kFALSE, 30., 180., 0, 0, 1, 
						floatLocationRightUp2D,500,500,optGenerator, optPeriod);
	
		cout << "line " << __LINE__ << endl;
		histoTrueEtaCaloConvPhotonConvRPtE->GetZaxis()->SetRangeUser(FindSmallestEntryIn2D(histoTrueEtaCaloConvPhotonConvRPtE), 
																	 histoTrueEtaCaloConvPhotonConvRPtE->GetMaximum());

		PlotStandard2D( histoTrueEtaCaloConvPhotonConvRPtE , 
						Form("%s/TrueEtaCaloConvPhoton_ConvR_PtElectron_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()), 
						"", "#it{R}_{conv, cluster}", "#it{p}_{T, e clus}",   
						kTRUE, 0.0, 10., kFALSE, 30., 180., 0, 0, 1, 
						floatLocationRightUp2D,500,500,optGenerator, optPeriod);
		cout << "line " << __LINE__ << endl;
		histoTruePi0CaloConvPhotonConvRAlphaE->GetZaxis()->SetRangeUser(FindSmallestEntryIn2D(histoTruePi0CaloConvPhotonConvRAlphaE), 
																	 histoTruePi0CaloConvPhotonConvRAlphaE->GetMaximum());

		PlotStandard2D( histoTruePi0CaloConvPhotonConvRAlphaE , 
						Form("%s/TruePi0CaloConvPhoton_ConvR_AlphaElectron_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()), 
						"", "#it{R}_{conv, cluster}", "#it{#it{#alpha}}_{e, clus} to other conversion partner",   
						kFALSE, 0.0, 10., kFALSE, 30., 180., 0, 0, 1, 
						floatLocationLeftDown2D,500,500,optGenerator, optPeriod);
		
		cout << "line " << __LINE__ << endl;
		histoTrueEtaCaloConvPhotonConvRAlphaE->GetZaxis()->SetRangeUser(FindSmallestEntryIn2D(histoTrueEtaCaloConvPhotonConvRAlphaE), 
																	 histoTrueEtaCaloConvPhotonConvRAlphaE->GetMaximum());
		PlotStandard2D( histoTrueEtaCaloConvPhotonConvRAlphaE , 
						Form("%s/TrueEtaCaloConvPhoton_ConvR_AlphaElectron_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()), 
						"", "#it{R}_{conv, cluster}", "#it{#it{#alpha}}_{e, clus} to other conversion partner",   
						kFALSE, 0.0, 10., kFALSE, 30., 180., 0, 0, 1, 
						floatLocationLeftDown2D,500,500,optGenerator, optPeriod);
		
		cout << "line " << __LINE__ << endl;
		histoTruePi0CaloConvPhotonPtEAlphaE->GetZaxis()->SetRangeUser(FindSmallestEntryIn2D(histoTruePi0CaloConvPhotonPtEAlphaE), 
																	 histoTruePi0CaloConvPhotonPtEAlphaE->GetMaximum());
		PlotStandard2D( histoTruePi0CaloConvPhotonPtEAlphaE , 
						Form("%s/TruePi0CaloConvPhoton_PtElectron_AlphaElectron_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()), 
						"", "#it{p}_{T, e clus}", "#it{#it{#alpha}}_{e, clus} to other conversion partner",   
						kFALSE, 0.0, 10., kTRUE, 0., 15., 0, 0, 1, 
						floatLocationRightDown2D,500,500,optGenerator, optPeriod);
		
		cout << "line " << __LINE__ << endl;
		histoTrueEtaCaloConvPhotonPtEAlphaE->GetZaxis()->SetRangeUser(FindSmallestEntryIn2D(histoTrueEtaCaloConvPhotonPtEAlphaE), 
																	 histoTrueEtaCaloConvPhotonPtEAlphaE->GetMaximum());
		PlotStandard2D( histoTrueEtaCaloConvPhotonPtEAlphaE , 
						Form("%s/TrueEtaCaloConvPhoton_PtElectron_AlphaElectron_%s.%s",outputDirectory.Data(),optGenerator.Data(),suffix.Data()), 
						"", "#it{p}_{T, e clus}", "#it{#it{#alpha}}_{e, clus} to other conversion partner",   
						kFALSE, 0.0, 10., kTRUE, 0., 15., 0, 0, 1, 
						floatLocationRightDown2D,500,500,optGenerator, optPeriod);
		
	}	
	
	//**********************************************************************************
	//****************** Plot 1D projections for pi0 alpha distributions ***************
	//**********************************************************************************		
	TCanvas * canvasSinglePtSlice = new TCanvas("canvasSinglePtSlice","",10,10,500,500);  // gives the page size
	canvasSinglePtSlice->SetLeftMargin(0.09);
	canvasSinglePtSlice->SetRightMargin(0.01);
	canvasSinglePtSlice->SetBottomMargin(0.06);
	canvasSinglePtSlice->SetTopMargin(0.01);
	canvasSinglePtSlice->SetLogy(1);
	canvasSinglePtSlice->cd();
	
	for (Int_t i = startBinPi0; i < 15; i++){
		DrawGammaSetMarker(histoBGPi0AlphaPtMCCpTBins[i],20,0.8, kBlack , kBlack);
		DrawAutoGammaHisto( histoBGPi0AlphaPtMCCpTBins[i],
						"", "#it{#alpha} ","dN/d#it{#alpha}",
						kFALSE, 1000.,1e5*FindSmallestEntryIn1D(histoTruePi0AlphaPtMCpTBins[i]),
						kTRUE,FindSmallestEntryIn1D(histoTruePi0AlphaPtMCpTBins[i]), FindLargestEntryIn1D(histoBGPi0AlphaPtMCCpTBins[i]) ,
						kFALSE, 0.,60.);  
		histoBGPi0AlphaPtMCCpTBins[i]->Draw("hist,e1");
		DrawGammaSetMarker(histoTruePi0AlphaPtMCpTBins[i],20,0.8, kRed+2 , kRed+2);
		histoTruePi0AlphaPtMCpTBins[i]->Draw("pe1,same");
		
		TLegend* legendAlpha = new TLegend( 0.75,0.82,0.93,0.91);
		legendAlpha->SetTextSize(0.03);         
		legendAlpha->SetFillColor(0);
		legendAlpha->SetLineColor(0);
		legendAlpha->AddEntry(histoTruePi0AlphaPtMCpTBins[i],"True #pi^{0}");
		legendAlpha->AddEntry(histoBGPi0AlphaPtMCCpTBins[i],"BG");
		legendAlpha->Draw();
		
		TLatex *labelPtrange = new TLatex(0.55,0.94,Form("%2.2f GeV/#it{c}< #it{p}_{T,#pi^{0}} < %2.2f GeV/#it{c}",ptBins[i], ptBins[i+1]));
		SetStyleTLatex( labelPtrange, 0.03,4);
		labelPtrange->Draw();

		canvasSinglePtSlice->Update();
		canvasSinglePtSlice->SaveAs(Form("%s/Pi0_Alpha_%s_pTBin%d.%s",outputDirectory.Data(), optGenerator.Data(), i,suffix.Data()));
	}
	
	//**********************************************************************************
	//****************** Plot 1D projections for eta alpha distributions ***************
	//**********************************************************************************		
	for (Int_t i = startBinEta; i < 15; i++){
		DrawGammaSetMarker(histoBGEtaAlphaPtMCCpTBins[i],20,0.8, kBlack , kBlack);
		DrawAutoGammaHisto( histoBGEtaAlphaPtMCCpTBins[i],
						"", "#it{#alpha} ","dN/d#it{#alpha}",
						kFALSE, 1000.,1e5*FindSmallestEntryIn1D(histoTrueEtaAlphaPtMCpTBins[i]),
						kTRUE,FindSmallestEntryIn1D(histoTrueEtaAlphaPtMCpTBins[i]),FindLargestEntryIn1D(histoBGEtaAlphaPtMCCpTBins[i]) , 
						kFALSE, 0.,60.);  
		histoBGEtaAlphaPtMCCpTBins[i]->Draw("hist,e1");
		DrawGammaSetMarker(histoTrueEtaAlphaPtMCpTBins[i],20,0.8, kRed+2 , kRed+2);
		histoTrueEtaAlphaPtMCpTBins[i]->Draw("pe1,same");
		
		TLegend* legendAlpha = new TLegend( 0.75,0.82,0.93,0.91);
		legendAlpha->SetTextSize(0.03);         
		legendAlpha->SetFillColor(0);
		legendAlpha->SetLineColor(0);
		legendAlpha->AddEntry(histoTrueEtaAlphaPtMCpTBins[i],"True #eta");
		legendAlpha->AddEntry(histoBGEtaAlphaPtMCCpTBins[i],"BG");
		legendAlpha->Draw();
		
		TLatex *labelPtrange = new TLatex(0.55,0.94,Form("%2.2f GeV/#it{c}< #it{p}_{T,#eta} < %2.2f GeV/#it{c}",ptBins[i], ptBins[i+1]));
		SetStyleTLatex( labelPtrange, 0.03,4);
		labelPtrange->Draw();

		canvasSinglePtSlice->Update();
		canvasSinglePtSlice->SaveAs(Form("%s/Eta_Alpha_%s_pTBin%d.%s",outputDirectory.Data(), optGenerator.Data(), i,suffix.Data()));
	}
	
	//**********************************************************************************
	//****************** Plot 1D projections for pi0 rapidity distributions ************
	//**********************************************************************************		
	for (Int_t i = startBinPi0; i < 15; i++){
		DrawGammaSetMarker(histoRecPi0RapidityPtDatapTBins[i],20,0.8, kBlack , kBlack);
		DrawAutoGammaHisto( histoRecPi0RapidityPtDatapTBins[i],
						"", "y ","dN/dy",
						kFALSE, 1000.,1e5*FindSmallestEntryIn1D(histoTruePi0AlphaPtMCpTBins[i]),
						kTRUE,FindSmallestEntryIn1D(histoTruePi0AlphaPtMCpTBins[i]), FindLargestEntryIn1D(histoRecPi0RapidityPtDatapTBins[i]) ,
						kFALSE, 0.,60.);  
		histoRecPi0RapidityPtDatapTBins[i]->Draw("hist,e1");
		DrawGammaSetMarker(histoRecPi0RapidityPtMCpTBins[i],24,0.8, kBlue+1 , kBlue+1);
		histoRecPi0RapidityPtMCpTBins[i]->Draw("hist,e1,same");

		DrawGammaSetMarker(histoTruePi0RapidityPtMCpTBins[i],20,0.8, kRed+2 , kRed+2);
		histoTruePi0RapidityPtMCpTBins[i]->Draw("pe1,same");
		
		TLegend* legendAlpha = new TLegend( 0.7,0.82,0.93,0.91);
		legendAlpha->SetTextSize(0.03);         
		legendAlpha->SetFillColor(0);
		legendAlpha->SetLineColor(0);
		legendAlpha->AddEntry(histoTruePi0RapidityPtMCpTBins[i],"True #pi^{0}");
		legendAlpha->AddEntry(histoRecPi0RapidityPtDatapTBins[i],"Rec Data");
		legendAlpha->AddEntry(histoRecPi0RapidityPtMCpTBins[i],"Rec MC");
		legendAlpha->Draw();
		
		TLatex *labelPtrange = new TLatex(0.55,0.94,Form("%2.2f GeV/#it{c}< #it{p}_{T,#pi^{0}} < %2.2f GeV/#it{c}",ptBins[i], ptBins[i+1]));
		SetStyleTLatex( labelPtrange, 0.03,4);
		labelPtrange->Draw();

		canvasSinglePtSlice->Update();
		canvasSinglePtSlice->SaveAs(Form("%s/Pi0_Rapidity_%s_pTBin_%d.%s",outputDirectory.Data(), optGenerator.Data(), i, suffix.Data()));
	}

	//**********************************************************************************
	//****************** Plot 1D projections for eta rapidity distributions ************
	//**********************************************************************************			
	for (Int_t i = startBinPi0; i < 15; i++){
		DrawGammaSetMarker(histoRecEtaRapidityPtDatapTBins[i],20,0.8, kBlack , kBlack);
		DrawAutoGammaHisto( histoRecEtaRapidityPtDatapTBins[i],
						"", "y ","dN/dy",
						kFALSE, 1000.,1e5*FindSmallestEntryIn1D(histoTrueEtaAlphaPtMCpTBins[i]),
						kTRUE,FindSmallestEntryIn1D(histoTrueEtaAlphaPtMCpTBins[i]), FindLargestEntryIn1D(histoRecEtaRapidityPtDatapTBins[i]) ,
						kFALSE, 0.,60.);  
		histoRecEtaRapidityPtDatapTBins[i]->Draw("hist,e1");
		DrawGammaSetMarker(histoRecEtaRapidityPtMCpTBins[i],24,0.8, kBlue+1 , kBlue+1);
		histoRecEtaRapidityPtMCpTBins[i]->Draw("hist,e1,same");

		DrawGammaSetMarker(histoTrueEtaRapidityPtMCpTBins[i],20,0.8, kRed+2 , kRed+2);
		histoTrueEtaRapidityPtMCpTBins[i]->Draw("pe1,same");
		
		TLegend* legendAlpha = new TLegend( 0.7,0.82,0.93,0.91);
		legendAlpha->SetTextSize(0.03);         
		legendAlpha->SetFillColor(0);
		legendAlpha->SetLineColor(0);
		legendAlpha->AddEntry(histoTrueEtaRapidityPtMCpTBins[i],"True #eta");
		legendAlpha->AddEntry(histoRecEtaRapidityPtDatapTBins[i],"Rec Data");
		legendAlpha->AddEntry(histoRecEtaRapidityPtMCpTBins[i],"Rec MC");
		legendAlpha->Draw();
		
		TLatex *labelPtrange = new TLatex(0.55,0.94,Form("%2.2f GeV/#it{c}< #it{p}_{T,#eta} < %2.2f GeV/#it{c}",ptBins[i], ptBins[i+1]));
		SetStyleTLatex( labelPtrange, 0.03,4);
		labelPtrange->Draw();

		canvasSinglePtSlice->Update();
		canvasSinglePtSlice->SaveAs(Form("%s/Eta_Rapidity_%s_pTBin%d.%s",outputDirectory.Data(), optGenerator.Data(), i,suffix.Data()));
	}

	if (enableExtConvCaloQA){
		canvasSinglePtSlice->SetLogy(0);
		canvasSinglePtSlice->SetBottomMargin(0.09);
		DrawGammaSetMarker(histoTruePi0CaloConvPhotonConvR,20,0.5, kBlack , kBlack);
		DrawAutoGammaHisto( histoTruePi0CaloConvPhotonConvR,
						"", "R_{conv, clus} ","dN/dR",
						kFALSE, 1000.,1e-5,
						kTRUE,0, 1.2*FindLargestEntryIn1D(histoTruePi0CaloConvPhotonConvR) ,
						kFALSE, 0.,60.);  
		canvasSinglePtSlice->Update();
		canvasSinglePtSlice->SaveAs(Form("%s/ConvR_ClusE_%s.%s",outputDirectory.Data(), optGenerator.Data(), suffix.Data()));

		DrawGammaSetMarker(histoTruePi0CaloConvPhotonAlphaE_FullR,20,0.5, kBlack , kBlack);
		DrawAutoGammaHisto( histoTruePi0CaloConvPhotonAlphaE_FullR,
						"", "#it{#alpha}_{e, clus} ","dN/d#it{#alpha}",
						kFALSE, 1000.,1e-5,
						kTRUE,0, 1.2*FindLargestEntryIn1D(histoTruePi0CaloConvPhotonAlphaE_FullR) ,
						kFALSE, 0.,60.);  
		
		DrawGammaSetMarker(histoTruePi0CaloConvPhotonAlphaE_R180360,24,0.5, kBlue+2 , kBlue+2);
		histoTruePi0CaloConvPhotonAlphaE_R180360->Draw("same,e1,p");
        DrawGammaSetMarker(histoTruePi0CaloConvPhotonAlphaE_R360460,25,0.5, kGreen+2 , kGreen+2);
        histoTruePi0CaloConvPhotonAlphaE_R360460->Draw("same,e1,p");
        DrawGammaSetMarker(histoTruePi0CaloConvPhotonAlphaE_RS180,20,0.5, kRed+2 , kRed+2);
		histoTruePi0CaloConvPhotonAlphaE_RS180->Draw("same,e1,p");

		TLegend* legendAlphaClus = new TLegend( 0.15,0.79,0.32,0.91);
		legendAlphaClus->SetTextSize(0.03);         
		legendAlphaClus->SetFillColor(0);
		legendAlphaClus->SetLineColor(0);
		legendAlphaClus->AddEntry(histoTruePi0CaloConvPhotonAlphaE_FullR,"R_{conv, clus} > 0 cm");
		legendAlphaClus->AddEntry(histoTruePi0CaloConvPhotonAlphaE_RS180,"R_{conv, clus} < 180 cm");
        legendAlphaClus->AddEntry(histoTruePi0CaloConvPhotonAlphaE_R180360,"180 cm < R_{conv, clus} < 360 cm");
        legendAlphaClus->AddEntry(histoTruePi0CaloConvPhotonAlphaE_R360460,"360 cm < R_{conv, clus} < 460 cm");
		legendAlphaClus->Draw();

		canvasSinglePtSlice->Update();
		canvasSinglePtSlice->SaveAs(Form("%s/Alpha_ClusE_RSlices_%s.%s",outputDirectory.Data(), optGenerator.Data(), suffix.Data()));

		DrawGammaSetMarker(histoTruePi0CaloConvPhotonAlphaE_FullPt,20,0.5, kBlack , kBlack);
		DrawAutoGammaHisto( histoTruePi0CaloConvPhotonAlphaE_FullPt,
						"", "#it{#alpha}_{e, clus} ","dN/d#it{#alpha}",
						kFALSE, 1000.,1e-5,
						kTRUE,0, 1.2*FindLargestEntryIn1D(histoTruePi0CaloConvPhotonAlphaE_FullPt) ,
						kFALSE, 0.,60.);  
		
		DrawGammaSetMarker(histoTruePi0CaloConvPhotonAlphaE_Pt0T500,24,0.5, kBlue+2 , kBlue+2);
		histoTruePi0CaloConvPhotonAlphaE_Pt0T500->Draw("same,e1,p");
		DrawGammaSetMarker(histoTruePi0CaloConvPhotonAlphaE_Pt500T1,20,0.5, kRed+2 , kRed+2);
		histoTruePi0CaloConvPhotonAlphaE_Pt500T1->Draw("same,e1,p");
		DrawGammaSetMarker(histoTruePi0CaloConvPhotonAlphaE_Pt1T2,24,0.5, kViolet+2 , kViolet+2);
		histoTruePi0CaloConvPhotonAlphaE_Pt1T2->Draw("same,e1,p");
		DrawGammaSetMarker(histoTruePi0CaloConvPhotonAlphaE_PtA2,20,0.5, kCyan+2 , kCyan+2);
		histoTruePi0CaloConvPhotonAlphaE_PtA2->Draw("same,e1,p");

		TLegend* legendAlphaClusPt = new TLegend( 0.15,0.72,0.32,0.91);
		legendAlphaClusPt->SetTextSize(0.03);         
		legendAlphaClusPt->SetFillColor(0);
		legendAlphaClusPt->SetLineColor(0);
		legendAlphaClusPt->AddEntry(histoTruePi0CaloConvPhotonAlphaE_FullPt,"0 GeV/c < p_{T, clus E} < 25 GeV/c");
		legendAlphaClusPt->AddEntry(histoTruePi0CaloConvPhotonAlphaE_Pt0T500,"0.0 GeV/c < p_{T, clus E} < 0.5 GeV/c");
		legendAlphaClusPt->AddEntry(histoTruePi0CaloConvPhotonAlphaE_Pt500T1,"0.5 GeV/c < p_{T, clus E} < 1 GeV/c");
		legendAlphaClusPt->AddEntry(histoTruePi0CaloConvPhotonAlphaE_Pt1T2,"1 GeV/c < p_{T, clus E} < 2 GeV/c");
		legendAlphaClusPt->AddEntry(histoTruePi0CaloConvPhotonAlphaE_PtA2,"2 GeV/c < p_{T, clus E} < 25 GeV/c");
		legendAlphaClusPt->Draw();

		canvasSinglePtSlice->Update();
		canvasSinglePtSlice->SaveAs(Form("%s/Alpha_ClusE_PtSlices_%s.%s",outputDirectory.Data(), optGenerator.Data(), suffix.Data()));

		canvasSinglePtSlice->SetLogy(1);
		DrawGammaSetMarker(histoTruePi0CaloConvPhotonPtE_FullR,20,0.5, kBlack , kBlack);
		DrawAutoGammaHisto( histoTruePi0CaloConvPhotonPtE_FullR,
						"", "p_{T,e, clus} ","dN/dp_{T}",
						kFALSE, 1000.,1e-5,
						kTRUE,1, 1.2*FindLargestEntryIn1D(histoTruePi0CaloConvPhotonPtE_FullR) ,
						kFALSE, 0.,60.);  
		DrawGammaSetMarker(histoTruePi0CaloConvPhotonPtE_RL180,24,0.5, kBlue+2 , kBlue+2);
		histoTruePi0CaloConvPhotonPtE_RL180->Draw("same,e1,p");
		DrawGammaSetMarker(histoTruePi0CaloConvPhotonPtE_RS180,20,0.5, kRed+2 , kRed+2);
		histoTruePi0CaloConvPhotonPtE_RS180->Draw("same,e1,p");

		TLegend* legendPtClus = new TLegend( 0.6,0.82,0.93,0.91);
		legendPtClus->SetTextSize(0.03);         
		legendPtClus->SetFillColor(0);
		legendPtClus->SetLineColor(0);
		legendPtClus->AddEntry(histoTruePi0CaloConvPhotonPtE_FullR,"R_{conv, clus} > 0 cm");
		legendPtClus->AddEntry(histoTruePi0CaloConvPhotonPtE_RL180,"R_{conv, clus} > 180 cm");
		legendPtClus->AddEntry(histoTruePi0CaloConvPhotonPtE_RS180,"R_{conv, clus} < 180 cm");
		legendPtClus->Draw();

		canvasSinglePtSlice->Update();
		canvasSinglePtSlice->SaveAs(Form("%s/PtE_ClusE_%s.%s",outputDirectory.Data(), optGenerator.Data(), suffix.Data()));
	}
	
	Bool_t enableMesonTree = kFALSE;
	cout << "here" << endl;
	TList *listTreeMCCalo 		= (TList*)HistosGammaConversionMC->FindObject(Form("%s True ClusterComb tree",cutSel.Data()));
	if(listTreeMCCalo != NULL) enableMesonTree = kTRUE;
	cout << "here" << endl;
	if (enableMesonTree){
		
		//Declaration of leaves types
		Float_t         InvMass;
		Float_t         RConv;
		Float_t         OpenAngleRPrimVtx;
		Float_t         InvMassRTOF;
		Float_t         Pt;
		UChar_t         cat;
		
		TTree *clusterTree 		= (TTree*)listTreeMCCalo->FindObject("ESD_ConvGamma_Pt_Dcaz_R_Eta");
		clusterTree->SetBranchAddress("InvMass",&InvMass);
		clusterTree->SetBranchAddress("RConv",&RConv);
		clusterTree->SetBranchAddress("OpenAngleRPrimVtx",&OpenAngleRPrimVtx);
		clusterTree->SetBranchAddress("InvMassRTOF",&InvMassRTOF);
		clusterTree->SetBranchAddress("Pt",&Pt);
		clusterTree->SetBranchAddress("cat",&cat);

		TH2F* hist2DTrueGamma_InvMassPrimVtx_InvMassRTof = new TH2F("hist2DTrueGamma_InvMassPrimVtx_InvMassRTof","", 800, 0, 0.8,  1200, 0, 1.2);
		TH2F* hist2DTruePi0Prim_InvMassPrimVtx_InvMassRTof = new TH2F("hist2DTruePi0Prim_InvMassPrimVtx_InvMassRTof","", 800, 0, 0.8,  1200, 0, 1.2);
		TH2F* hist2DTrueEtaPrim_InvMassPrimVtx_InvMassRTof = new TH2F("hist2DTrueEtaPrim_InvMassPrimVtx_InvMassRTof","", 800, 0, 0.8,  1200, 0, 1.2);
		
		Long64_t nEntriesPairs 				= clusterTree->GetEntries();
		Int_t nPairs 							= 0;
		cout << "Number of P to be processed: " << nEntriesPairs << endl;
		//    (*ptGammaAssosiatedPi0)[k]
		Long64_t nbytesPairs 					= 0;

		for (Long64_t i=0; i<nEntriesPairs;i++) {
			nbytesPairs 						+= clusterTree->GetEntry(i); 
			if (cat ==0){
				hist2DTrueGamma_InvMassPrimVtx_InvMassRTof->Fill(InvMass,InvMassRTOF);
			}
			if (cat == 1){
				hist2DTruePi0Prim_InvMassPrimVtx_InvMassRTof->Fill(InvMass,InvMassRTOF);
			}
			if (cat == 2){
				hist2DTrueEtaPrim_InvMassPrimVtx_InvMassRTof->Fill(InvMass,InvMassRTOF);
			}

			
		}
		TFile* fOutputGamma = new TFile(Form("%s/outputTreeProjections.root",outputDirectory.Data()),"RECREATE");\
			hist2DTrueGamma_InvMassPrimVtx_InvMassRTof->Write();
			hist2DTruePi0Prim_InvMassPrimVtx_InvMassRTof->Write();
			hist2DTrueEtaPrim_InvMassPrimVtx_InvMassRTof->Write();
		fOutputGamma->Write();
		fOutputGamma->Close();

		
	}
	
}
