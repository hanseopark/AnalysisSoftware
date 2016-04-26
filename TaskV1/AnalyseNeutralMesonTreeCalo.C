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
void AnalyseNeutralMesonTreeCalo(		TString fileNameMC = "myOutput", 
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
	
	TString outputDirectory 	= Form("%s/%s/%s/AnalyseNeutralMesonCaloTree",cutSel.Data(),optEnergy.Data(),suffix.Data());
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
   	TString nameMainDir 							= "";
	if (mode == 9 || mode == 0) nameMainDir 		= "GammaConvV1";
	else if (mode == 2 || mode == 3) nameMainDir 	= "GammaConvCalo";
	else if (mode == 4 || mode == 5) nameMainDir 	= "GammaCalo";
		
	
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

		TH2F* hist2DTrueGamma_InvMassPrimVtx_InvMassRTof = new TH2F("hist2DTrueGamma_InvMassPrimVtx_InvMassRTof","", 800, 0, 0.8,  1800, 0, 1.8);
		hist2DTrueGamma_InvMassPrimVtx_InvMassRTof->GetXaxis()->SetTitle("M_{clus,clus}, R_{primVtx}");
		hist2DTrueGamma_InvMassPrimVtx_InvMassRTof->GetYaxis()->SetTitle("M_{clus,clus}, R_{375 cm}");
		TH2F* hist2DTrueGamma_InvMassPrimVtx_RConv = new TH2F("hist2DTrueGamma_InvMassPrimVtx_RConv","", 800, 0, 0.8,  800, 0, 400);
		hist2DTrueGamma_InvMassPrimVtx_RConv->GetXaxis()->SetTitle("M_{clus,clus}, R_{primVtx}");
		hist2DTrueGamma_InvMassPrimVtx_RConv->GetYaxis()->SetTitle("R_{conv}");
		TH2F* hist2DTrueGamma_InvMassRTof_RConv = new TH2F("hist2DTrueGamma_InvMassRTof_RConv","",  1800, 0, 1.8,  800, 0, 400);
		hist2DTrueGamma_InvMassRTof_RConv->GetXaxis()->SetTitle("M_{clus,clus}, R_{375 cm}");
		hist2DTrueGamma_InvMassRTof_RConv->GetYaxis()->SetTitle("R_{conv}");
		TH2F* hist2DTrueGamma_InvMass_OpenAngle = new TH2F("hist2DTrueGamma_InvMass_OpenAngle","",  800, 0, 0.8,  400, 0, TMath::Pi());
		hist2DTrueGamma_InvMass_OpenAngle->GetXaxis()->SetTitle("M_{clus,clus}, R_{primVtx}");
		hist2DTrueGamma_InvMass_OpenAngle->GetYaxis()->SetTitle("#Theta_{Open}");
		TH2F* hist2DTrueGamma_InvMassRTof_OpenAngle = new TH2F("hist2DTrueGamma_InvMassRTof_OpenAngle","",  1800, 0, 1.8,  400, 0, TMath::Pi());
		hist2DTrueGamma_InvMassRTof_OpenAngle->GetXaxis()->SetTitle("M_{clus,clus}, R_{375 cm}");
		hist2DTrueGamma_InvMassRTof_OpenAngle->GetYaxis()->SetTitle("#Theta_{Open}");

		TH2F* hist2DTrueGammaRConvS250cm_InvMassPrimVtx_InvMassRTof = new TH2F("hist2DTrueGammaRConvS250cm_InvMassPrimVtx_InvMassRTof","", 800, 0, 0.8,  1800, 0, 1.8);
		hist2DTrueGammaRConvS250cm_InvMassPrimVtx_InvMassRTof->GetXaxis()->SetTitle("M_{clus,clus}, R_{primVtx}");
		hist2DTrueGammaRConvS250cm_InvMassPrimVtx_InvMassRTof->GetYaxis()->SetTitle("M_{clus,clus}, R_{375 cm}");
		TH2F* hist2DTrueGammaRConvG250cm_InvMassPrimVtx_InvMassRTof = new TH2F("hist2DTrueGammaRConvG250cm_InvMassPrimVtx_InvMassRTof","", 800, 0, 0.8,  1800, 0, 1.8);
		hist2DTrueGammaRConvG250cm_InvMassPrimVtx_InvMassRTof->GetXaxis()->SetTitle("M_{clus,clus}, R_{primVtx}");
		hist2DTrueGammaRConvG250cm_InvMassPrimVtx_InvMassRTof->GetYaxis()->SetTitle("M_{clus,clus}, R_{375 cm}");

		TH2F* hist2DTruePi0Prim_InvMassPrimVtx_InvMassRTof = new TH2F("hist2DTruePi0Prim_InvMassPrimVtx_InvMassRTof","", 800, 0, 0.8,  1800, 0, 1.8);
		hist2DTruePi0Prim_InvMassPrimVtx_InvMassRTof->GetXaxis()->SetTitle("M_{clus,clus}, R_{primVtx}");
		hist2DTruePi0Prim_InvMassPrimVtx_InvMassRTof->GetYaxis()->SetTitle("M_{clus,clus}, R_{375 cm}");
		TH2F* hist2DTruePi0Prim_InvMassPrimVtx_RConv = new TH2F("hist2DTruePi0Prim_InvMassPrimVtx_RConv","", 800, 0, 0.8,  800, 0, 400);
		hist2DTruePi0Prim_InvMassPrimVtx_RConv->GetXaxis()->SetTitle("M_{clus,clus}, R_{primVtx}");
		hist2DTruePi0Prim_InvMassPrimVtx_RConv->GetYaxis()->SetTitle("R_{conv}");
		TH2F* hist2DTruePi0Prim_InvMassRTof_RConv = new TH2F("hist2DTruePi0Prim_InvMassRTof_RConv","",  1800, 0, 1.8,  800, 0, 400);
		hist2DTruePi0Prim_InvMassRTof_RConv->GetXaxis()->SetTitle("M_{clus,clus}, R_{375 cm}");
		hist2DTruePi0Prim_InvMassRTof_RConv->GetYaxis()->SetTitle("R_{conv}");
		TH2F* hist2DTruePi0Prim_InvMass_OpenAngle = new TH2F("hist2DTruePi0Prim_InvMass_OpenAngle","",  800, 0, 0.8,  400, 0, TMath::Pi());
		hist2DTruePi0Prim_InvMass_OpenAngle->GetXaxis()->SetTitle("M_{clus,clus}, R_{primVtx}");
		hist2DTruePi0Prim_InvMass_OpenAngle->GetYaxis()->SetTitle("#Theta_{Open}");


		TH2F* hist2DTruePi0SecK0s_InvMassPrimVtx_InvMassRTof = new TH2F("hist2DTruePi0SecK0s_InvMassPrimVtx_InvMassRTof","", 800, 0, 0.8,  1800, 0, 1.8);
		hist2DTruePi0SecK0s_InvMassPrimVtx_InvMassRTof->GetXaxis()->SetTitle("M_{clus,clus}, R_{primVtx}");
		hist2DTruePi0SecK0s_InvMassPrimVtx_InvMassRTof->GetYaxis()->SetTitle("M_{clus,clus}, R_{375 cm}");
		TH2F* hist2DTruePi0SecLambda_InvMassPrimVtx_InvMassRTof = new TH2F("hist2DTruePi0SecLambda_InvMassPrimVtx_InvMassRTof","", 800, 0, 0.8,  1800, 0, 1.8);
		hist2DTruePi0SecLambda_InvMassPrimVtx_InvMassRTof->GetXaxis()->SetTitle("M_{clus,clus}, R_{primVtx}");
		hist2DTruePi0SecLambda_InvMassPrimVtx_InvMassRTof->GetYaxis()->SetTitle("M_{clus,clus}, R_{375 cm}");
		TH2F* hist2DTruePi0SecMat_InvMassPrimVtx_InvMassRTof = new TH2F("hist2DTruePi0SecMat_InvMassPrimVtx_InvMassRTof","", 800, 0, 0.8,  1800, 0, 1.8);
		hist2DTruePi0SecMat_InvMassPrimVtx_InvMassRTof->GetXaxis()->SetTitle("M_{clus,clus}, R_{primVtx}");
		hist2DTruePi0SecMat_InvMassPrimVtx_InvMassRTof->GetYaxis()->SetTitle("M_{clus,clus}, R_{375 cm}");
		TH2F* hist2DTruePi0SecMat_InvMassPrimVtx_RConv = new TH2F("hist2DTruePi0SecMat_InvMassPrimVtx_RConv","", 800, 0, 0.8,  800, 0, 400);
		hist2DTruePi0SecMat_InvMassPrimVtx_RConv->GetXaxis()->SetTitle("M_{clus,clus}, R_{primVtx}");
		hist2DTruePi0SecMat_InvMassPrimVtx_RConv->GetYaxis()->SetTitle("R_{conv}");
		TH2F* hist2DTruePi0SecMat_InvMassRTof_RConv = new TH2F("hist2DTruePi0SecMat_InvMassRTof_RConv","",  1800, 0, 1.8,  800, 0, 400);
		hist2DTruePi0SecMat_InvMassRTof_RConv->GetXaxis()->SetTitle("M_{clus,clus}, R_{375 cm}");
		hist2DTruePi0SecMat_InvMassRTof_RConv->GetYaxis()->SetTitle("R_{conv}");
		
		TH2F* hist2DTrueEtaPrim_InvMassPrimVtx_InvMassRTof = new TH2F("hist2DTrueEtaPrim_InvMassPrimVtx_InvMassRTof","", 800, 0, 0.8,  1500, 0, 3.0);
		hist2DTrueEtaPrim_InvMassPrimVtx_InvMassRTof->GetXaxis()->SetTitle("M_{clus,clus}, R_{primVtx}");
		hist2DTrueEtaPrim_InvMassPrimVtx_InvMassRTof->GetYaxis()->SetTitle("M_{clus,clus}, R_{375 cm}");
		TH2F* hist2DTrueEtaPrim_InvMassPrimVtx_RConv = new TH2F("hist2DTrueEtaPrim_InvMassPrimVtx_RConv","", 800, 0, 0.8,  800, 0, 400);
		hist2DTrueEtaPrim_InvMassPrimVtx_RConv->GetXaxis()->SetTitle("M_{clus,clus}, R_{primVtx}");
		hist2DTrueEtaPrim_InvMassPrimVtx_RConv->GetYaxis()->SetTitle("R_{conv}");
		TH2F* hist2DTrueEtaPrim_InvMassRTof_RConv = new TH2F("hist2DTrueEtaPrim_InvMassRTof_RConv","",  1500, 0, 3.0,  800, 0, 400);
		hist2DTrueEtaPrim_InvMassRTof_RConv->GetXaxis()->SetTitle("M_{clus,clus}, R_{375 cm}");
		hist2DTrueEtaPrim_InvMassRTof_RConv->GetYaxis()->SetTitle("R_{conv}");
		
		Long64_t nEntriesPairs 				= clusterTree->GetEntries();
		Int_t nPairs 							= 0;
		cout << "Number of P to be processed: " << nEntriesPairs << endl;
		//    (*ptGammaAssosiatedPi0)[k]
		Long64_t nbytesPairs 					= 0;

		for (Long64_t i=0; i<nEntriesPairs;i++) {
			nbytesPairs 						+= clusterTree->GetEntry(i); 
			if (Pt >0.5) {
				if (cat ==0){
					hist2DTrueGamma_InvMassPrimVtx_InvMassRTof->Fill(InvMass,InvMassRTOF);
					hist2DTrueGamma_InvMassPrimVtx_RConv->Fill(InvMass,RConv);
					hist2DTrueGamma_InvMassRTof_RConv->Fill(InvMassRTOF,RConv);
					if (RConv < 250){
						hist2DTrueGammaRConvS250cm_InvMassPrimVtx_InvMassRTof->Fill(InvMass,InvMassRTOF);
					} else {
						hist2DTrueGammaRConvG250cm_InvMassPrimVtx_InvMassRTof->Fill(InvMass,InvMassRTOF);
					}	
					hist2DTrueGamma_InvMass_OpenAngle->Fill(InvMass,OpenAngleRPrimVtx);
					hist2DTrueGamma_InvMassRTof_OpenAngle->Fill(InvMassRTOF,OpenAngleRPrimVtx);
				}
			}	
			if (Pt >1) {	
				if (cat == 1){
					hist2DTruePi0Prim_InvMassPrimVtx_InvMassRTof->Fill(InvMass,InvMassRTOF);
					hist2DTruePi0Prim_InvMassPrimVtx_RConv->Fill(InvMass,RConv);
					hist2DTruePi0Prim_InvMassRTof_RConv->Fill(InvMassRTOF,RConv);
					hist2DTruePi0Prim_InvMass_OpenAngle->Fill(InvMass,OpenAngleRPrimVtx);

					
				}
				if (cat == 3){
					hist2DTruePi0SecK0s_InvMassPrimVtx_InvMassRTof->Fill(InvMass,InvMassRTOF);
				}
				if (cat == 5){
					hist2DTruePi0SecLambda_InvMassPrimVtx_InvMassRTof->Fill(InvMass,InvMassRTOF);
				}
				if (cat == 6){
					hist2DTruePi0SecMat_InvMassPrimVtx_InvMassRTof->Fill(InvMass,InvMassRTOF);
					hist2DTruePi0SecMat_InvMassPrimVtx_RConv->Fill(InvMass,RConv);
					hist2DTruePi0SecMat_InvMassRTof_RConv->Fill(InvMassRTOF,RConv);
				}
				if (cat == 2){
					hist2DTrueEtaPrim_InvMassPrimVtx_InvMassRTof->Fill(InvMass,InvMassRTOF);
					hist2DTrueEtaPrim_InvMassPrimVtx_RConv->Fill(InvMass,RConv);
					hist2DTrueEtaPrim_InvMassRTof_RConv->Fill(InvMassRTOF,RConv);
				}
			}
			
		}
		TFile* fOutputGamma = new TFile(Form("%s/outputTreeProjections.root",outputDirectory.Data()),"RECREATE");\
			hist2DTrueGamma_InvMassPrimVtx_InvMassRTof->Write();
			hist2DTrueGamma_InvMassPrimVtx_RConv->Write();
			hist2DTrueGamma_InvMassRTof_RConv->Write();
			hist2DTrueGammaRConvS250cm_InvMassPrimVtx_InvMassRTof->Write();
			hist2DTrueGammaRConvG250cm_InvMassPrimVtx_InvMassRTof->Write();
			hist2DTrueGamma_InvMassRTof_OpenAngle->Write();
			hist2DTrueGamma_InvMass_OpenAngle->Write();
			hist2DTruePi0Prim_InvMassPrimVtx_InvMassRTof->Write();
			hist2DTruePi0Prim_InvMassPrimVtx_RConv->Write();
			hist2DTruePi0Prim_InvMassRTof_RConv->Write();
			hist2DTruePi0Prim_InvMass_OpenAngle->Write();

			hist2DTruePi0SecK0s_InvMassPrimVtx_InvMassRTof->Write();
			hist2DTruePi0SecLambda_InvMassPrimVtx_InvMassRTof->Write();
			hist2DTruePi0SecMat_InvMassPrimVtx_InvMassRTof->Write();
			hist2DTruePi0SecMat_InvMassPrimVtx_RConv->Write();
			hist2DTruePi0SecMat_InvMassRTof_RConv->Write();

			hist2DTrueEtaPrim_InvMassPrimVtx_InvMassRTof->Write();
			hist2DTrueEtaPrim_InvMassPrimVtx_RConv->Write();
			hist2DTrueEtaPrim_InvMassRTof_RConv->Write();

		fOutputGamma->Write();
		fOutputGamma->Close();

		
	}
	
}
