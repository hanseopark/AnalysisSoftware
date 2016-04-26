/********************************************************************************************************
****** 		provided by Gamma Conversion Group,  													*****
******		Lucia Leardini , lucia.leardini@cern.ch													*****
********************************************************************************************************/
#include <Riostream.h>
#include <fstream>
#include "TMath.h"
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TPaveLabel.h>
#include <TPaveText.h>
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
#include "/home/admin1/leardini/photonconv/AnalysisSoftware/CommonHeaders/PlottingGammaConversionHistos.h"
#include "/home/admin1/leardini/photonconv/AnalysisSoftware/CommonHeaders/PlottingGammaConversionAdditional.h"
#include "/home/admin1/leardini/photonconv/AnalysisSoftware/CommonHeaders/FittingGammaConversion.h"
#include "/home/admin1/leardini/photonconv/AnalysisSoftware/CommonHeaders/ConversionFunctions.h"

void PutGraphTogetherQA(Bool_t partialOutput = kFALSE, TString outputDir = "plotQA_27Nov",TString suffix = "pdf"){
	gROOT->Reset();	
	gROOT->SetStyle("Plain");
	
	StyleSettingsThesis();
	SetPlotStyle();
	TGaxis::SetMaxDigits(6);
		
	Color_t	 colorComb_0010_Data	= kRed+1;
	Color_t	 colorComb_1020_Data	= kOrange-2;
	Color_t	 colorComb_2040_Data	= kGreen+2;
	Color_t	 colorComb_4060_Data	= kCyan+2;
//	Color_t	 colorComb6080_Data	= kBlue+1;

	Color_t	 colorComb_0010_MC	= kRed+1;
	Color_t	 colorComb_1020_MC	= kOrange-2;
	Color_t	 colorComb_2040_MC	= kGreen+2;
	Color_t	 colorComb_4060_MC	= kCyan+2;
//	Color_t	 colorComb6080_MC	= kBlue+1;
	
	Style_t 	markerStyle_0010_Data = 20 ;
	Style_t 	markerStyle_1020_Data = 21 ;
	Style_t 	markerStyle_2040_Data = 33 ;
	Style_t 	markerStyle_4060_Data = 20 ;
//	Style_t 	markerStyle6080_Data = 21 ;

	Style_t 	markerStyle_0010_MC = 20 ;
	Style_t 	markerStyle_1020_MC = 25 ;
	Style_t 	markerStyle_2040_MC = 33 ;
	Style_t 	markerStyle_4060_MC = 20 ;
//	Style_t 	markerStyle6080_MC = 21 ;

	Style_t 	markerStyleMC_0010_Data 	= 24 ;
	Style_t 	markerStyleMC_1020_Data 	= 21 ;
	Style_t 	markerStyleMC_2040_Data 	= 27 ;
	Style_t 	markerStyleMC_4060_Data 	= 24 ;
//	Style_t 	markerStyleMC6080_Data 	= 25 ;

	Style_t 	markerStyleMC_0010_MC 	= 24 ;
	Style_t 	markerStyleMC_1020_MC 	= 25 ;
	Style_t 	markerStyleMC_2040_MC 	= 27 ;
	Style_t 	markerStyleMC_4060_MC 	= 24 ;
//	Style_t 	markerStyleMC6080_MC 	= 25 ;	


	Size_t 	markerSize_0010_Data 	= 2.;
	Size_t 	markerSize_1020_Data 	= 2.;
	Size_t 	markerSize_2040_Data 	= 2.5;
	Size_t 	markerSize_4060_Data 	= 2.;
//	Size_t 	markerSize6080_Data 	= 2.;
	
	Size_t 	markerSize_0010_MC 	= 2.;
	Size_t 	markerSize_1020_MC 	= 2.;
	Size_t 	markerSize_2040_MC 	= 2.5;
	Size_t 	markerSize_4060_MC 	= 2.;
//	Size_t 	markerSize6080_MC 	= 2.;

	TString collisionSystemPbPb_0010_Data = "0-10% Pb-Pb #sqrt{s_{NN}} = 2.76 TeV";;		
	TString collisionSystemPbPb_1020_Data = "10-20% Pb-Pb #sqrt{s_{NN}} = 2.76 TeV";		
	TString collisionSystemPbPb_2040_Data = "20-40% Pb-Pb #sqrt{s_{NN}} = 2.76 TeV";		
	TString collisionSystemPbPb_4060_Data = "40-60% Pb-Pb #sqrt{s_{NN}} = 2.76 TeV";		
//	TString collisionSystemPbPb6080_Data = "60-80% Pb-Pb #sqrt{s_{NN}} = 2.76 TeV";		

	TString collisionSystem0010 = "0-10%";		
	TString collisionSystem1020 = "10-20%";	
	TString collisionSystem2040 = "20-40%";		
	TString collisionSystem4060 = "40-60%";		
//	TString collisionSystem6080 = "60-80%";		
	
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////// 			Data 			/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	const char* InputnameData = Form("%s/OutputAnalysisQA_Data.root",outputDir.Data());
	TFile* QAOutputFileData = new TFile(InputnameData);
	cout << "Loading data histos..."<< endl;

	TGraphErrors* EventsPerRun_0010_Data = (TGraphErrors*)QAOutputFileData->Get("EventsPerRun_0010");
	TGraphErrors* EventsPerRun_1020_Data = (TGraphErrors*)QAOutputFileData->Get("EventsPerRun_1020");
	TGraphErrors* EventsPerRun_2040_Data = (TGraphErrors*)QAOutputFileData->Get("EventsPerRun_2040");
	TGraphErrors* EventsPerRun_4060_Data = (TGraphErrors*)QAOutputFileData->Get("EventsPerRun_4060");

	//electron cut out of the pion
	TGraphErrors* NSigmadEdxElectron_0010_Data = (TGraphErrors*)QAOutputFileData->Get("meanNSigmadEdxElectron_0010");
	TGraphErrors* NSigmadEdxElectron_1020_Data = (TGraphErrors*)QAOutputFileData->Get("meanNSigmadEdxElectron_1020");
	TGraphErrors* NSigmadEdxElectron_2040_Data = (TGraphErrors*)QAOutputFileData->Get("meanNSigmadEdxElectron_2040");
	TGraphErrors* NSigmadEdxElectron_4060_Data = (TGraphErrors*)QAOutputFileData->Get("meanNSigmadEdxElectron_4060");

	TGraphErrors* NSigmadEdxPositron_0010_Data = (TGraphErrors*)QAOutputFileData->Get("meanNSigmadEdxPositron_0010");
	TGraphErrors* NSigmadEdxPositron_1020_Data = (TGraphErrors*)QAOutputFileData->Get("meanNSigmadEdxPositron_1020");
	TGraphErrors* NSigmadEdxPositron_2040_Data = (TGraphErrors*)QAOutputFileData->Get("meanNSigmadEdxPositron_2040");
	TGraphErrors* NSigmadEdxPositron_4060_Data = (TGraphErrors*)QAOutputFileData->Get("meanNSigmadEdxPositron_4060");

	TGraphErrors* sigmaNSigmadEdxElectron_0010_Data = (TGraphErrors*)QAOutputFileData->Get("sigmaelNSigmaElectronCut_0010");
	TGraphErrors* sigmaNSigmadEdxElectron_1020_Data = (TGraphErrors*)QAOutputFileData->Get("sigmaelNSigmaElectronCut_1020");
	TGraphErrors* sigmaNSigmadEdxElectron_2040_Data = (TGraphErrors*)QAOutputFileData->Get("sigmaelNSigmaElectronCut_2040");
	TGraphErrors* sigmaNSigmadEdxElectron_4060_Data = (TGraphErrors*)QAOutputFileData->Get("sigmaelNSigmaElectronCut_4060");

	TGraphErrors* sigmaNSigmadEdxPositron_0010_Data = (TGraphErrors*)QAOutputFileData->Get("sigmaposNSigmaElectronCut_0010");
	TGraphErrors* sigmaNSigmadEdxPositron_1020_Data = (TGraphErrors*)QAOutputFileData->Get("sigmaposNSigmaElectronCut_1020");
	TGraphErrors* sigmaNSigmadEdxPositron_2040_Data = (TGraphErrors*)QAOutputFileData->Get("sigmaposNSigmaElectronCut_2040");
	TGraphErrors* sigmaNSigmadEdxPositron_4060_Data = (TGraphErrors*)QAOutputFileData->Get("sigmaposNSigmaElectronCut_4060");

	TGraphErrors* widthNSigmadEdxElectron_0010_Data = (TGraphErrors*)QAOutputFileData->Get("widthelNSigmaElectronCut_0010");
	TGraphErrors* widthNSigmadEdxElectron_1020_Data = (TGraphErrors*)QAOutputFileData->Get("widthelNSigmaElectronCut_1020");
	TGraphErrors* widthNSigmadEdxElectron_2040_Data = (TGraphErrors*)QAOutputFileData->Get("widthelNSigmaElectronCut_2040");
	TGraphErrors* widthNSigmadEdxElectron_4060_Data = (TGraphErrors*)QAOutputFileData->Get("widthelNSigmaElectronCut_4060");

	TGraphErrors* widthNSigmadEdxPositron_0010_Data = (TGraphErrors*)QAOutputFileData->Get("widthposNSigmaElectronCut_0010");
	TGraphErrors* widthNSigmadEdxPositron_1020_Data = (TGraphErrors*)QAOutputFileData->Get("widthposNSigmaElectronCut_1020");
	TGraphErrors* widthNSigmadEdxPositron_2040_Data = (TGraphErrors*)QAOutputFileData->Get("widthposNSigmaElectronCut_2040");
	TGraphErrors* widthNSigmadEdxPositron_4060_Data = (TGraphErrors*)QAOutputFileData->Get("widthposNSigmaElectronCut_4060");


//	cut over the proton line
	TGraphErrors* NSigmadEdxElectronProtonCut_0010_Data = (TGraphErrors*)QAOutputFileData->Get("meanelNSigmaProtonCut_0010");
	TGraphErrors* NSigmadEdxElectronProtonCut_1020_Data = (TGraphErrors*)QAOutputFileData->Get("meanelNSigmaProtonCut_1020");
	TGraphErrors* NSigmadEdxElectronProtonCut_2040_Data = (TGraphErrors*)QAOutputFileData->Get("meanelNSigmaProtonCut_2040");
	TGraphErrors* NSigmadEdxElectronProtonCut_4060_Data = (TGraphErrors*)QAOutputFileData->Get("meanelNSigmaProtonCut_4060");

	TGraphErrors* NSigmadEdxPositronProtonCut_0010_Data = (TGraphErrors*)QAOutputFileData->Get("meanposNSigmaProtonCut_0010");
	TGraphErrors* NSigmadEdxPositronProtonCut_1020_Data = (TGraphErrors*)QAOutputFileData->Get("meanposNSigmaProtonCut_1020");
	TGraphErrors* NSigmadEdxPositronProtonCut_2040_Data = (TGraphErrors*)QAOutputFileData->Get("meanposNSigmaProtonCut_2040");
	TGraphErrors* NSigmadEdxPositronProtonCut_4060_Data = (TGraphErrors*)QAOutputFileData->Get("meanposNSigmaProtonCut_4060");

	TGraphErrors* sigmaNSigmadEdxElectronProtonCut_0010_Data = (TGraphErrors*)QAOutputFileData->Get("sigmaelNSigmaProtonCut_0010");
	TGraphErrors* sigmaNSigmadEdxElectronProtonCut_1020_Data = (TGraphErrors*)QAOutputFileData->Get("sigmaelNSigmaProtonCut_1020");
	TGraphErrors* sigmaNSigmadEdxElectronProtonCut_2040_Data = (TGraphErrors*)QAOutputFileData->Get("sigmaelNSigmaProtonCut_2040");
	TGraphErrors* sigmaNSigmadEdxElectronProtonCut_4060_Data = (TGraphErrors*)QAOutputFileData->Get("sigmaelNSigmaProtonCut_4060");

	TGraphErrors* sigmaNSigmadEdxPositronProtonCut_0010_Data = (TGraphErrors*)QAOutputFileData->Get("sigmaposNSigmaProtonCut_0010");
	TGraphErrors* sigmaNSigmadEdxPositronProtonCut_1020_Data = (TGraphErrors*)QAOutputFileData->Get("sigmaposNSigmaProtonCut_1020");
	TGraphErrors* sigmaNSigmadEdxPositronProtonCut_2040_Data = (TGraphErrors*)QAOutputFileData->Get("sigmaposNSigmaProtonCut_2040");
	TGraphErrors* sigmaNSigmadEdxPositronProtonCut_4060_Data = (TGraphErrors*)QAOutputFileData->Get("sigmaposNSigmaProtonCut_4060");

	TGraphErrors* widthNSigmadEdxElectronProtonCut_0010_Data = (TGraphErrors*)QAOutputFileData->Get("widthelNSigmaProtonCut_0010");
	TGraphErrors* widthNSigmadEdxElectronProtonCut_1020_Data = (TGraphErrors*)QAOutputFileData->Get("widthelNSigmaProtonCut_1020");
	TGraphErrors* widthNSigmadEdxElectronProtonCut_2040_Data = (TGraphErrors*)QAOutputFileData->Get("widthelNSigmaProtonCut_2040");
	TGraphErrors* widthNSigmadEdxElectronProtonCut_4060_Data = (TGraphErrors*)QAOutputFileData->Get("widthelNSigmaProtonCut_4060");

	TGraphErrors* widthNSigmadEdxPositronProtonCut_0010_Data = (TGraphErrors*)QAOutputFileData->Get("widthposNSigmaProtonCut_0010");
	TGraphErrors* widthNSigmadEdxPositronProtonCut_1020_Data = (TGraphErrors*)QAOutputFileData->Get("widthposNSigmaProtonCut_1020");
	TGraphErrors* widthNSigmadEdxPositronProtonCut_2040_Data = (TGraphErrors*)QAOutputFileData->Get("widthposNSigmaProtonCut_2040");
	TGraphErrors* widthNSigmadEdxPositronProtonCut_4060_Data = (TGraphErrors*)QAOutputFileData->Get("widthposNSigmaProtonCut_4060");

	//photon pt
	TGraphErrors* PhotonPt_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonPt_0010");
	TGraphErrors* PhotonPt_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonPt_1020");
	TGraphErrors* PhotonPt_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonPt_2040");
	TGraphErrors* PhotonPt_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonPt_4060");

	//eta
	TGraphErrors* ElectronEta_0010_Data = (TGraphErrors*)QAOutputFileData->Get("ElectronEta_0010");
	TGraphErrors* ElectronEta_1020_Data = (TGraphErrors*)QAOutputFileData->Get("ElectronEta_1020");
	TGraphErrors* ElectronEta_2040_Data = (TGraphErrors*)QAOutputFileData->Get("ElectronEta_2040");
	TGraphErrors* ElectronEta_4060_Data = (TGraphErrors*)QAOutputFileData->Get("ElectronEta_4060");

	TGraphErrors* PositronEta_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PositronEta_0010");
	TGraphErrors* PositronEta_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PositronEta_1020");
	TGraphErrors* PositronEta_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PositronEta_2040");
	TGraphErrors* PositronEta_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PositronEta_4060");

	TGraphErrors* PhotonEta_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonEta_0010");
	TGraphErrors* PhotonEta_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonEta_1020");
	TGraphErrors* PhotonEta_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonEta_2040");
	TGraphErrors* PhotonEta_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonEta_4060");

	//eta neg
	TGraphErrors* ElectronEtaNeg_0010_Data = (TGraphErrors*)QAOutputFileData->Get("ElectronEtaNeg_0010");
	TGraphErrors* ElectronEtaNeg_1020_Data = (TGraphErrors*)QAOutputFileData->Get("ElectronEtaNeg_1020");
	TGraphErrors* ElectronEtaNeg_2040_Data = (TGraphErrors*)QAOutputFileData->Get("ElectronEtaNeg_2040");
	TGraphErrors* ElectronEtaNeg_4060_Data = (TGraphErrors*)QAOutputFileData->Get("ElectronEtaNeg_4060");

	TGraphErrors* PositronEtaNeg_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PositronEtaNeg_0010");
	TGraphErrors* PositronEtaNeg_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PositronEtaNeg_1020");
	TGraphErrors* PositronEtaNeg_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PositronEtaNeg_2040");
	TGraphErrors* PositronEtaNeg_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PositronEtaNeg_4060");

	TGraphErrors* PhotonEtaNeg_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonEtaNeg_0010");
	TGraphErrors* PhotonEtaNeg_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonEtaNeg_1020");
	TGraphErrors* PhotonEtaNeg_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonEtaNeg_2040");
	TGraphErrors* PhotonEtaNeg_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonEtaNeg_4060");

	//eta pos
	TGraphErrors* ElectronEtaPos_0010_Data = (TGraphErrors*)QAOutputFileData->Get("ElectronEtaPos_0010");
	TGraphErrors* ElectronEtaPos_1020_Data = (TGraphErrors*)QAOutputFileData->Get("ElectronEtaPos_1020");
	TGraphErrors* ElectronEtaPos_2040_Data = (TGraphErrors*)QAOutputFileData->Get("ElectronEtaPos_2040");
	TGraphErrors* ElectronEtaPos_4060_Data = (TGraphErrors*)QAOutputFileData->Get("ElectronEtaPos_4060");

	TGraphErrors* PositronEtaPos_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PositronEtaPos_0010");
	TGraphErrors* PositronEtaPos_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PositronEtaPos_1020");
	TGraphErrors* PositronEtaPos_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PositronEtaPos_2040");
	TGraphErrors* PositronEtaPos_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PositronEtaPos_4060");

	TGraphErrors* PhotonEtaPos_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonEtaPos_0010");
	TGraphErrors* PhotonEtaPos_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonEtaPos_1020");
	TGraphErrors* PhotonEtaPos_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonEtaPos_2040");
	TGraphErrors* PhotonEtaPos_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonEtaPos_4060");

		
	cout << "Loading data histos...TPC sectors only Nev" << endl;
	TGraphErrors* PhotonAllSectorsEtaNeg_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonAllSectorsn_0010");
	TGraphErrors* PhotonAllSectorsEtaNeg_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonAllSectorsn_1020");
	TGraphErrors* PhotonAllSectorsEtaNeg_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonAllSectorsn_2040");
	TGraphErrors* PhotonAllSectorsEtaNeg_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonAllSectorsn_4060");

	TGraphErrors* PhotonAllSectorsEtaPos_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonAllSectorsp_0010");
	TGraphErrors* PhotonAllSectorsEtaPos_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonAllSectorsp_1020");
	TGraphErrors* PhotonAllSectorsEtaPos_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonAllSectorsp_2040");
	TGraphErrors* PhotonAllSectorsEtaPos_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonAllSectorsp_4060");

	TGraphErrors* PhotonSector0EtaNeg_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector0n_0010");
	TGraphErrors* PhotonSector0EtaNeg_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector0n_1020");
	TGraphErrors* PhotonSector0EtaNeg_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector0n_2040");
	TGraphErrors* PhotonSector0EtaNeg_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector0n_4060");

	TGraphErrors* PhotonSector0EtaPos_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector0p_0010");
	TGraphErrors* PhotonSector0EtaPos_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector0p_1020");
	TGraphErrors* PhotonSector0EtaPos_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector0p_2040");
	TGraphErrors* PhotonSector0EtaPos_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector0p_4060");

	TGraphErrors* PhotonSector1EtaNeg_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector1n_0010");
	TGraphErrors* PhotonSector1EtaNeg_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector1n_1020");	
	TGraphErrors* PhotonSector1EtaNeg_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector1n_2040");
	TGraphErrors* PhotonSector1EtaNeg_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector1n_4060");

	TGraphErrors* PhotonSector1EtaPos_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector1p_0010");
	TGraphErrors* PhotonSector1EtaPos_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector1p_1020");
	TGraphErrors* PhotonSector1EtaPos_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector1p_2040");
	TGraphErrors* PhotonSector1EtaPos_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector1p_4060");
//	TGraphErrors* PhotonSector1EtaPos6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector1p_6080");

	TGraphErrors* PhotonSector2EtaNeg_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector2n_0010");
	TGraphErrors* PhotonSector2EtaNeg_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector2n_1020");
	TGraphErrors* PhotonSector2EtaNeg_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector2n_2040");
	TGraphErrors* PhotonSector2EtaNeg_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector2n_4060");
//	TGraphErrors* PhotonSector2EtaNeg6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector2n_6080");

	TGraphErrors* PhotonSector2EtaPos_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector2p_0010");
	TGraphErrors* PhotonSector2EtaPos_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector2p_1020");
	TGraphErrors* PhotonSector2EtaPos_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector2p_2040");
	TGraphErrors* PhotonSector2EtaPos_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector2p_4060");
//	TGraphErrors* PhotonSector2EtaPos6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector2p_6080");

	TGraphErrors* PhotonSector3EtaNeg_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector3n_0010");
	TGraphErrors* PhotonSector3EtaNeg_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector3n_1020");
	TGraphErrors* PhotonSector3EtaNeg_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector3n_2040");
	TGraphErrors* PhotonSector3EtaNeg_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector3n_4060");
//	TGraphErrors* PhotonSector3EtaNeg6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector3n_6080");

	TGraphErrors* PhotonSector3EtaPos_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector3p_0010");
	TGraphErrors* PhotonSector3EtaPos_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector3p_1020");
	TGraphErrors* PhotonSector3EtaPos_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector3p_2040");
	TGraphErrors* PhotonSector3EtaPos_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector3p_4060");
//	TGraphErrors* PhotonSector3EtaPos6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector3p_6080");

	TGraphErrors* PhotonSector4EtaNeg_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector4n_0010");
	TGraphErrors* PhotonSector4EtaNeg_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector4n_1020");
	TGraphErrors* PhotonSector4EtaNeg_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector4n_2040");
	TGraphErrors* PhotonSector4EtaNeg_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector4n_4060");
//	TGraphErrors* PhotonSector4EtaNeg6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector4n_6080");

	TGraphErrors* PhotonSector4EtaPos_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector4p_0010");
	TGraphErrors* PhotonSector4EtaPos_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector4p_1020");
	TGraphErrors* PhotonSector4EtaPos_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector4p_2040");
	TGraphErrors* PhotonSector4EtaPos_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector4p_4060");
//	TGraphErrors* PhotonSector4EtaPos6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector4p_6080");

	TGraphErrors* PhotonSector5EtaNeg_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector5n_0010");
	TGraphErrors* PhotonSector5EtaNeg_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector5n_1020");
	TGraphErrors* PhotonSector5EtaNeg_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector5n_2040");
	TGraphErrors* PhotonSector5EtaNeg_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector5n_4060");
//	TGraphErrors* PhotonSector5EtaNeg6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector5n_6080");

	TGraphErrors* PhotonSector5EtaPos_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector5p_0010");
	TGraphErrors* PhotonSector5EtaPos_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector5p_1020");
	TGraphErrors* PhotonSector5EtaPos_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector5p_2040");
	TGraphErrors* PhotonSector5EtaPos_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector5p_4060");
//	TGraphErrors* PhotonSector5EtaPos6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector5p_6080");

	TGraphErrors* PhotonSector6EtaNeg_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector6n_0010");
	TGraphErrors* PhotonSector6EtaNeg_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector6n_1020");
	TGraphErrors* PhotonSector6EtaNeg_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector6n_2040");
	TGraphErrors* PhotonSector6EtaNeg_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector6n_4060");
//	TGraphErrors* PhotonSector6EtaNeg6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector6n_6080");

	TGraphErrors* PhotonSector6EtaPos_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector6p_0010");
	TGraphErrors* PhotonSector6EtaPos_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector6p_1020");
	TGraphErrors* PhotonSector6EtaPos_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector6p_2040");
	TGraphErrors* PhotonSector6EtaPos_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector6p_4060");
//	TGraphErrors* PhotonSector6EtaPos6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector6p_6080");

	TGraphErrors* PhotonSector7EtaNeg_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector7n_0010");
	TGraphErrors* PhotonSector7EtaNeg_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector7n_1020");
	TGraphErrors* PhotonSector7EtaNeg_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector7n_2040");
	TGraphErrors* PhotonSector7EtaNeg_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector7n_4060");
//	TGraphErrors* PhotonSector7EtaNeg6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector7n_6080");

	TGraphErrors* PhotonSector7EtaPos_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector7p_0010");
	TGraphErrors* PhotonSector7EtaPos_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector7p_1020");
	TGraphErrors* PhotonSector7EtaPos_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector7p_2040");
	TGraphErrors* PhotonSector7EtaPos_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector7p_4060");
//	TGraphErrors* PhotonSector7EtaPos6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector7p_6080");

	TGraphErrors* PhotonSector8EtaNeg_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector8n_0010");
	TGraphErrors* PhotonSector8EtaNeg_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector8n_1020");
	TGraphErrors* PhotonSector8EtaNeg_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector8n_2040");
	TGraphErrors* PhotonSector8EtaNeg_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector8n_4060");
//	TGraphErrors* PhotonSector8EtaNeg6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector8n_6080");

	TGraphErrors* PhotonSector8EtaPos_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector8p_0010");
	TGraphErrors* PhotonSector8EtaPos_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector8p_1020");
	TGraphErrors* PhotonSector8EtaPos_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector8p_2040");
	TGraphErrors* PhotonSector8EtaPos_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector8p_4060");
//	TGraphErrors* PhotonSector8EtaPos6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector8p_6080");

	TGraphErrors* PhotonSector9EtaNeg_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector9n_0010");
	TGraphErrors* PhotonSector9EtaNeg_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector9n_1020");
	TGraphErrors* PhotonSector9EtaNeg_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector9n_2040");
	TGraphErrors* PhotonSector9EtaNeg_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector9n_4060");
//	TGraphErrors* PhotonSector9EtaNeg6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector9n_6080");

	TGraphErrors* PhotonSector9EtaPos_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector9p_0010");
	TGraphErrors* PhotonSector9EtaPos_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector9p_1020");
	TGraphErrors* PhotonSector9EtaPos_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector9p_2040");
	TGraphErrors* PhotonSector9EtaPos_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector9p_4060");
//	TGraphErrors* PhotonSector9EtaPos6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector9p_6080");

	TGraphErrors* PhotonSector10EtaNeg_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector10n_0010");
	TGraphErrors* PhotonSector10EtaNeg_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector10n_1020");
	TGraphErrors* PhotonSector10EtaNeg_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector10n_2040");
	TGraphErrors* PhotonSector10EtaNeg_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector10n_4060");
//	TGraphErrors* PhotonSector10EtaNeg6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector10n_6080");

	TGraphErrors* PhotonSector10EtaPos_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector10p_0010");
	TGraphErrors* PhotonSector10EtaPos_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector10p_1020");
	TGraphErrors* PhotonSector10EtaPos_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector10p_2040");
	TGraphErrors* PhotonSector10EtaPos_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector10p_4060");
//	TGraphErrors* PhotonSector10EtaPos6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector10p_6080");

	TGraphErrors* PhotonSector11EtaNeg_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector11n_0010");
	TGraphErrors* PhotonSector11EtaNeg_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector11n_1020");
	TGraphErrors* PhotonSector11EtaNeg_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector11n_2040");
	TGraphErrors* PhotonSector11EtaNeg_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector11n_4060");
//	TGraphErrors* PhotonSector11EtaNeg6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector11n_6080");

	TGraphErrors* PhotonSector11EtaPos_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector11p_0010");
	TGraphErrors* PhotonSector11EtaPos_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector11p_1020");
	TGraphErrors* PhotonSector11EtaPos_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector11p_2040");
	TGraphErrors* PhotonSector11EtaPos_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector11p_4060");
//	TGraphErrors* PhotonSector11EtaPos6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector11p_6080");

	TGraphErrors* PhotonSector12EtaNeg_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector12n_0010");
	TGraphErrors* PhotonSector12EtaNeg_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector12n_1020");
	TGraphErrors* PhotonSector12EtaNeg_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector12n_2040");
	TGraphErrors* PhotonSector12EtaNeg_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector12n_4060");
//	TGraphErrors* PhotonSector12EtaNeg6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector12n_6080");

	TGraphErrors* PhotonSector12EtaPos_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector12p_0010");
	TGraphErrors* PhotonSector12EtaPos_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector12p_1020");
	TGraphErrors* PhotonSector12EtaPos_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector12p_2040");
	TGraphErrors* PhotonSector12EtaPos_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector12p_4060");
//	TGraphErrors* PhotonSector12EtaPos6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector12p_6080");

	TGraphErrors* PhotonSector13EtaNeg_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector13n_0010");
	TGraphErrors* PhotonSector13EtaNeg_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector13n_1020");
	TGraphErrors* PhotonSector13EtaNeg_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector13n_2040");
	TGraphErrors* PhotonSector13EtaNeg_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector13n_4060");
//	TGraphErrors* PhotonSector13EtaNeg6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector13n_6080");

	TGraphErrors* PhotonSector13EtaPos_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector13p_0010");
	TGraphErrors* PhotonSector13EtaPos_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector13p_1020");
	TGraphErrors* PhotonSector13EtaPos_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector13p_2040");
	TGraphErrors* PhotonSector13EtaPos_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector13p_4060");
//	TGraphErrors* PhotonSector13EtaPos6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector13p_6080");

	TGraphErrors* PhotonSector14EtaNeg_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector14n_0010");
	TGraphErrors* PhotonSector14EtaNeg_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector14n_1020");
	TGraphErrors* PhotonSector14EtaNeg_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector14n_2040");
	TGraphErrors* PhotonSector14EtaNeg_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector14n_4060");
//	TGraphErrors* PhotonSector14EtaNeg6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector14n_6080");

	TGraphErrors* PhotonSector14EtaPos_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector14p_0010");
	TGraphErrors* PhotonSector14EtaPos_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector14p_1020");
	TGraphErrors* PhotonSector14EtaPos_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector14p_2040");
	TGraphErrors* PhotonSector14EtaPos_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector14p_4060");
//	TGraphErrors* PhotonSector14EtaPos6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector14p_6080");

	TGraphErrors* PhotonSector15EtaNeg_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector15n_0010");
	TGraphErrors* PhotonSector15EtaNeg_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector15n_1020");
	TGraphErrors* PhotonSector15EtaNeg_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector15n_2040");
	TGraphErrors* PhotonSector15EtaNeg_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector15n_4060");
//	TGraphErrors* PhotonSector15EtaNeg6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector15n_6080");

	TGraphErrors* PhotonSector15EtaPos_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector15p_0010");
	TGraphErrors* PhotonSector15EtaPos_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector15p_1020");
	TGraphErrors* PhotonSector15EtaPos_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector15p_2040");
	TGraphErrors* PhotonSector15EtaPos_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector15p_4060");
//	TGraphErrors* PhotonSector15EtaPos6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector15p_6080");

	TGraphErrors* PhotonSector16EtaNeg_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector16n_0010");
	TGraphErrors* PhotonSector16EtaNeg_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector16n_1020");
	TGraphErrors* PhotonSector16EtaNeg_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector16n_2040");
	TGraphErrors* PhotonSector16EtaNeg_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector16n_4060");
//	TGraphErrors* PhotonSector16EtaNeg6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector16n_6080");

	TGraphErrors* PhotonSector16EtaPos_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector16p_0010");
	TGraphErrors* PhotonSector16EtaPos_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector16p_1020");
	TGraphErrors* PhotonSector16EtaPos_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector16p_2040");
	TGraphErrors* PhotonSector16EtaPos_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector16p_4060");
//	TGraphErrors* PhotonSector16EtaPos6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector16p_6080");

	TGraphErrors* PhotonSector17EtaNeg_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector17n_0010");
	TGraphErrors* PhotonSector17EtaNeg_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector17n_1020");
	TGraphErrors* PhotonSector17EtaNeg_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector17n_2040");
	TGraphErrors* PhotonSector17EtaNeg_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector17n_4060");
//	TGraphErrors* PhotonSector17EtaNeg6080_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector17n_6080");

	TGraphErrors* PhotonSector17EtaPos_0010_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector17p_0010");
	TGraphErrors* PhotonSector17EtaPos_1020_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector17p_1020");
	TGraphErrors* PhotonSector17EtaPos_2040_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector17p_2040");
	TGraphErrors* PhotonSector17EtaPos_4060_Data = (TGraphErrors*)QAOutputFileData->Get("PhotonSector17p_4060");
	

	cout << "Drawing the data plots" << endl;

	cout << "Events" << endl;
	DrawGammaSetMarkerTGraphErr( EventsPerRun_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr( EventsPerRun_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr( EventsPerRun_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr( EventsPerRun_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);

	//Nsigma
	cout << "Nsigma plots electron cut" << endl;
	//electron cut out of the pion
	DrawGammaSetMarkerTGraphErr( NSigmadEdxElectron_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxElectron_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxElectron_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxElectron_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);

	DrawGammaSetMarkerTGraphErr( NSigmadEdxPositron_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxPositron_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxPositron_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxPositron_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);

	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxElectron_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxElectron_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxElectron_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxElectron_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);

	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxPositron_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxPositron_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxPositron_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxPositron_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);

	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxElectron_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxElectron_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxElectron_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxElectron_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);

	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxPositron_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxPositron_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxPositron_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxPositron_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);

	
	cout << "Nsigma plots proton cut" << endl;
	DrawGammaSetMarkerTGraphErr( NSigmadEdxElectronProtonCut_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxElectronProtonCut_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxElectronProtonCut_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxElectronProtonCut_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);

	DrawGammaSetMarkerTGraphErr( NSigmadEdxPositronProtonCut_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxPositronProtonCut_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxPositronProtonCut_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxPositronProtonCut_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);

	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxElectronProtonCut_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxElectronProtonCut_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxElectronProtonCut_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxElectronProtonCut_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);

	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxPositronProtonCut_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxPositronProtonCut_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxPositronProtonCut_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxPositronProtonCut_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);

	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxElectronProtonCut_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxElectronProtonCut_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxElectronProtonCut_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxElectronProtonCut_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);

	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxPositronProtonCut_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxPositronProtonCut_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxPositronProtonCut_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxPositronProtonCut_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);


	//photon pt
	cout << "Photon pt" << endl;
	DrawGammaSetMarkerTGraphErr( PhotonPt_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr( PhotonPt_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr( PhotonPt_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr( PhotonPt_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);


	//eta
	cout << "Eta plots" << endl;
	DrawGammaSetMarkerTGraphErr( ElectronEta_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr( ElectronEta_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr( ElectronEta_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr( ElectronEta_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);

	DrawGammaSetMarkerTGraphErr( PositronEta_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr( PositronEta_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr( PositronEta_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr( PositronEta_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);

	DrawGammaSetMarkerTGraphErr( PhotonEta_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr( PhotonEta_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr( PhotonEta_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr( PhotonEta_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);

	//eta neg
	DrawGammaSetMarkerTGraphErr( ElectronEtaNeg_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr( ElectronEtaNeg_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr( ElectronEtaNeg_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr( ElectronEtaNeg_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);	

	DrawGammaSetMarkerTGraphErr( PositronEtaNeg_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr( PositronEtaNeg_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr( PositronEtaNeg_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr( PositronEtaNeg_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);

	DrawGammaSetMarkerTGraphErr( PhotonEtaNeg_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr( PhotonEtaNeg_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr( PhotonEtaNeg_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr( PhotonEtaNeg_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);

	//eta pos
	DrawGammaSetMarkerTGraphErr( ElectronEtaPos_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr( ElectronEtaPos_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr( ElectronEtaPos_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr( ElectronEtaPos_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);

	DrawGammaSetMarkerTGraphErr( PositronEtaPos_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr( PositronEtaPos_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr( PositronEtaPos_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr( PositronEtaPos_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);

	DrawGammaSetMarkerTGraphErr( PhotonEtaPos_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr( PhotonEtaPos_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr( PhotonEtaPos_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr( PhotonEtaPos_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);


	//number of gammas
	cout << "Number of gammas" << endl;
	DrawGammaSetMarkerTGraphErr(PhotonAllSectorsEtaNeg_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonAllSectorsEtaNeg_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonAllSectorsEtaNeg_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonAllSectorsEtaNeg_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonAllSectorsEtaNeg6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonAllSectorsEtaPos_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonAllSectorsEtaPos_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonAllSectorsEtaPos_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonAllSectorsEtaPos_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonAllSectorsEtaPos6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector0EtaNeg_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector0EtaNeg_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector0EtaNeg_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector0EtaNeg_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector0EtaNeg6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector0EtaPos_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector0EtaPos_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector0EtaPos_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector0EtaPos_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector0EtaPos6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);  

	DrawGammaSetMarkerTGraphErr(PhotonSector1EtaNeg_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector1EtaNeg_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);	
	DrawGammaSetMarkerTGraphErr(PhotonSector1EtaNeg_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector1EtaNeg_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector1EtaNeg6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector1EtaPos_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector1EtaPos_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector1EtaPos_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector1EtaPos_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector1EtaPos6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector2EtaNeg_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector2EtaNeg_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector2EtaNeg_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector2EtaNeg_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector2EtaNeg6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector2EtaPos_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector2EtaPos_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector2EtaPos_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector2EtaPos_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector2EtaPos6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector3EtaNeg_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector3EtaNeg_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector3EtaNeg_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector3EtaNeg_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector3EtaNeg6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector3EtaPos_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector3EtaPos_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector3EtaPos_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector3EtaPos_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector3EtaPos6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector4EtaNeg_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector4EtaNeg_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector4EtaNeg_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector4EtaNeg_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector4EtaNeg6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector4EtaPos_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector4EtaPos_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector4EtaPos_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector4EtaPos_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector4EtaPos6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector5EtaNeg_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector5EtaNeg_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector5EtaNeg_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector5EtaNeg_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector5EtaNeg6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector5EtaPos_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector5EtaPos_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector5EtaPos_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector5EtaPos_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector5EtaPos6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector6EtaNeg_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector6EtaNeg_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector6EtaNeg_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector6EtaNeg_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector6EtaNeg6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector6EtaPos_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector6EtaPos_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector6EtaPos_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector6EtaPos_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector6EtaPos6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector7EtaNeg_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector7EtaNeg_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector7EtaNeg_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector7EtaNeg_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector7EtaNeg6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector7EtaPos_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector7EtaPos_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector7EtaPos_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector7EtaPos_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector7EtaPos6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector8EtaNeg_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector8EtaNeg_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector8EtaNeg_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector8EtaNeg_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector8EtaNeg6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector8EtaPos_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector8EtaPos_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector8EtaPos_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector8EtaPos_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector8EtaPos6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector9EtaNeg_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector9EtaNeg_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector9EtaNeg_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector9EtaNeg_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector9EtaNeg6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector9EtaPos_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector9EtaPos_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector9EtaPos_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector9EtaPos_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector9EtaPos6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector10EtaNeg_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector10EtaNeg_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector10EtaNeg_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector10EtaNeg_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector10EtaNeg6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector10EtaPos_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector10EtaPos_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector10EtaPos_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector10EtaPos_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector10EtaPos6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector11EtaNeg_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector11EtaNeg_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector11EtaNeg_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector11EtaNeg_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector11EtaNeg6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector11EtaPos_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector11EtaPos_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector11EtaPos_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector11EtaPos_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector11EtaPos6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector12EtaNeg_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector12EtaNeg_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector12EtaNeg_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector12EtaNeg_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector12EtaNeg6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector12EtaPos_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector12EtaPos_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector12EtaPos_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector12EtaPos_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector12EtaPos6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector13EtaNeg_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector13EtaNeg_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector13EtaNeg_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector13EtaNeg_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector13EtaNeg6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector13EtaPos_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector13EtaPos_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector13EtaPos_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector13EtaPos_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector13EtaPos6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector14EtaNeg_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector14EtaNeg_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector14EtaNeg_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector14EtaNeg_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector14EtaNeg6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector14EtaPos_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector14EtaPos_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector14EtaPos_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector14EtaPos_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector14EtaPos6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector15EtaNeg_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector15EtaNeg_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector15EtaNeg_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector15EtaNeg_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector15EtaNeg6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector15EtaPos_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector15EtaPos_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector15EtaPos_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector15EtaPos_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector15EtaPos6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector16EtaNeg_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector16EtaNeg_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector16EtaNeg_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector16EtaNeg_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector16EtaNeg6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector16EtaPos_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector16EtaPos_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector16EtaPos_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector16EtaPos_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector16EtaPos6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector17EtaNeg_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector17EtaNeg_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector17EtaNeg_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector17EtaNeg_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector17EtaNeg6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);

	DrawGammaSetMarkerTGraphErr(PhotonSector17EtaPos_0010_Data,  markerStyle_0010_Data, markerSize_0010_Data, colorComb_0010_Data, colorComb_0010_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector17EtaPos_1020_Data,  markerStyle_1020_Data, markerSize_1020_Data, colorComb_1020_Data, colorComb_1020_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector17EtaPos_2040_Data,  markerStyle_2040_Data, markerSize_2040_Data, colorComb_2040_Data, colorComb_2040_Data);
	DrawGammaSetMarkerTGraphErr(PhotonSector17EtaPos_4060_Data,  markerStyle_4060_Data, markerSize_4060_Data, colorComb_4060_Data, colorComb_4060_Data);
//	DrawGammaSetMarkerTGraphErr(PhotonSector17EtaPos6080_Data,  markerStyle6080_Data, markerSize6080_Data, colorComb6080_Data, colorComb6080_Data);



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////// 			MC 			/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	const char* InputnameMC = Form("%s/OutputAnalysisQA_MC.root",outputDir.Data());
	TFile* QAOutputFileMC = new TFile(InputnameMC);

	cout << "Loading MC histos..."<< endl;

	TGraphErrors* EventsPerRun_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("EventsPerRun_0010");
	TGraphErrors* EventsPerRun_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("EventsPerRun_1020");
	TGraphErrors* EventsPerRun_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("EventsPerRun_2040");
	TGraphErrors* EventsPerRun_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("EventsPerRun_4060");

	//nsigma

	//electron cut out of the pion
	TGraphErrors* NSigmadEdxElectron_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("meanNSigmadEdxElectron_0010");
	TGraphErrors* NSigmadEdxElectron_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("meanNSigmadEdxElectron_1020");
	TGraphErrors* NSigmadEdxElectron_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("meanNSigmadEdxElectron_2040");
	TGraphErrors* NSigmadEdxElectron_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("meanNSigmadEdxElectron_4060");

	TGraphErrors* NSigmadEdxPositron_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("meanNSigmadEdxPositron_0010");
	TGraphErrors* NSigmadEdxPositron_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("meanNSigmadEdxPositron_1020");
	TGraphErrors* NSigmadEdxPositron_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("meanNSigmadEdxPositron_2040");
	TGraphErrors* NSigmadEdxPositron_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("meanNSigmadEdxPositron_4060");

	TGraphErrors* sigmaNSigmadEdxElectron_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("sigmaelNSigmaElectronCut_0010");
	TGraphErrors* sigmaNSigmadEdxElectron_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("sigmaelNSigmaElectronCut_1020");
	TGraphErrors* sigmaNSigmadEdxElectron_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("sigmaelNSigmaElectronCut_2040");
	TGraphErrors* sigmaNSigmadEdxElectron_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("sigmaelNSigmaElectronCut_4060");

	TGraphErrors* sigmaNSigmadEdxPositron_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("sigmaposNSigmaElectronCut_0010");
	TGraphErrors* sigmaNSigmadEdxPositron_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("sigmaposNSigmaElectronCut_1020");
	TGraphErrors* sigmaNSigmadEdxPositron_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("sigmaposNSigmaElectronCut_2040");
	TGraphErrors* sigmaNSigmadEdxPositron_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("sigmaposNSigmaElectronCut_4060");

	TGraphErrors* widthNSigmadEdxElectron_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("widthelNSigmaElectronCut_0010");
	TGraphErrors* widthNSigmadEdxElectron_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("widthelNSigmaElectronCut_1020");
	TGraphErrors* widthNSigmadEdxElectron_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("widthelNSigmaElectronCut_2040");
	TGraphErrors* widthNSigmadEdxElectron_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("widthelNSigmaElectronCut_4060");

	TGraphErrors* widthNSigmadEdxPositron_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("widthposNSigmaElectronCut_0010");
	TGraphErrors* widthNSigmadEdxPositron_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("widthposNSigmaElectronCut_1020");
	TGraphErrors* widthNSigmadEdxPositron_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("widthposNSigmaElectronCut_2040");
	TGraphErrors* widthNSigmadEdxPositron_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("widthposNSigmaElectronCut_4060");


	//cut over the proton line
	TGraphErrors* NSigmadEdxElectronProtonCut_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("meanelNSigmaProtonCut_0010");
	TGraphErrors* NSigmadEdxElectronProtonCut_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("meanelNSigmaProtonCut_1020");
	TGraphErrors* NSigmadEdxElectronProtonCut_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("meanelNSigmaProtonCut_2040");
	TGraphErrors* NSigmadEdxElectronProtonCut_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("meanelNSigmaProtonCut_4060");

	TGraphErrors* NSigmadEdxPositronProtonCut_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("meanposNSigmaProtonCut_0010");
	TGraphErrors* NSigmadEdxPositronProtonCut_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("meanposNSigmaProtonCut_1020");
	TGraphErrors* NSigmadEdxPositronProtonCut_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("meanposNSigmaProtonCut_2040");
	TGraphErrors* NSigmadEdxPositronProtonCut_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("meanposNSigmaProtonCut_4060");

	TGraphErrors* sigmaNSigmadEdxElectronProtonCut_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("sigmaelNSigmaProtonCut_0010");
	TGraphErrors* sigmaNSigmadEdxElectronProtonCut_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("sigmaelNSigmaProtonCut_1020");
	TGraphErrors* sigmaNSigmadEdxElectronProtonCut_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("sigmaelNSigmaProtonCut_2040");
	TGraphErrors* sigmaNSigmadEdxElectronProtonCut_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("sigmaelNSigmaProtonCut_4060");

	TGraphErrors* sigmaNSigmadEdxPositronProtonCut_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("sigmaposNSigmaProtonCut_0010");
	TGraphErrors* sigmaNSigmadEdxPositronProtonCut_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("sigmaposNSigmaProtonCut_1020");
	TGraphErrors* sigmaNSigmadEdxPositronProtonCut_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("sigmaposNSigmaProtonCut_2040");
	TGraphErrors* sigmaNSigmadEdxPositronProtonCut_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("sigmaposNSigmaProtonCut_4060");

	TGraphErrors* widthNSigmadEdxElectronProtonCut_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("widthelNSigmaProtonCut_0010");
	TGraphErrors* widthNSigmadEdxElectronProtonCut_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("widthelNSigmaProtonCut_1020");
	TGraphErrors* widthNSigmadEdxElectronProtonCut_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("widthelNSigmaProtonCut_2040");
	TGraphErrors* widthNSigmadEdxElectronProtonCut_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("widthelNSigmaProtonCut_4060");

	TGraphErrors* widthNSigmadEdxPositronProtonCut_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("widthposNSigmaProtonCut_0010");
	TGraphErrors* widthNSigmadEdxPositronProtonCut_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("widthposNSigmaProtonCut_1020");
	TGraphErrors* widthNSigmadEdxPositronProtonCut_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("widthposNSigmaProtonCut_2040");
	TGraphErrors* widthNSigmadEdxPositronProtonCut_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("widthposNSigmaProtonCut_4060");


	//photon pt
	TGraphErrors* PhotonPt_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonPt_0010");
	TGraphErrors* PhotonPt_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonPt_1020");
	TGraphErrors* PhotonPt_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonPt_2040");
	TGraphErrors* PhotonPt_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonPt_4060");


	//eta
	TGraphErrors* ElectronEta_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("ElectronEta_0010");
	TGraphErrors* ElectronEta_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("ElectronEta_1020");
	TGraphErrors* ElectronEta_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("ElectronEta_2040");
	TGraphErrors* ElectronEta_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("ElectronEta_4060");

	TGraphErrors* PositronEta_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PositronEta_0010");
	TGraphErrors* PositronEta_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PositronEta_1020");
	TGraphErrors* PositronEta_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PositronEta_2040");
	TGraphErrors* PositronEta_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PositronEta_4060");

	TGraphErrors* PhotonEta_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonEta_0010");
	TGraphErrors* PhotonEta_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonEta_1020");
	TGraphErrors* PhotonEta_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonEta_2040");
	TGraphErrors* PhotonEta_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonEta_4060");

	//eta neg
	TGraphErrors* ElectronEtaNeg_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("ElectronEtaNeg_0010");
	TGraphErrors* ElectronEtaNeg_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("ElectronEtaNeg_1020");
	TGraphErrors* ElectronEtaNeg_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("ElectronEtaNeg_2040");
	TGraphErrors* ElectronEtaNeg_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("ElectronEtaNeg_4060");

	TGraphErrors* PositronEtaNeg_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PositronEtaNeg_0010");
	TGraphErrors* PositronEtaNeg_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PositronEtaNeg_1020");
	TGraphErrors* PositronEtaNeg_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PositronEtaNeg_2040");
	TGraphErrors* PositronEtaNeg_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PositronEtaNeg_4060");

	TGraphErrors* PhotonEtaNeg_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonEtaNeg_0010");
	TGraphErrors* PhotonEtaNeg_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonEtaNeg_1020");
	TGraphErrors* PhotonEtaNeg_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonEtaNeg_2040");
	TGraphErrors* PhotonEtaNeg_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonEtaNeg_4060");

	//eta pos
	TGraphErrors* ElectronEtaPos_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("ElectronEtaPos_0010");
	TGraphErrors* ElectronEtaPos_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("ElectronEtaPos_1020");
	TGraphErrors* ElectronEtaPos_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("ElectronEtaPos_2040");
	TGraphErrors* ElectronEtaPos_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("ElectronEtaPos_4060");

	TGraphErrors* PositronEtaPos_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PositronEtaPos_0010");
	TGraphErrors* PositronEtaPos_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PositronEtaPos_1020");
	TGraphErrors* PositronEtaPos_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PositronEtaPos_2040");
	TGraphErrors* PositronEtaPos_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PositronEtaPos_4060");

	TGraphErrors* PhotonEtaPos_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonEtaPos_0010");
	TGraphErrors* PhotonEtaPos_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonEtaPos_1020");
	TGraphErrors* PhotonEtaPos_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonEtaPos_2040");
	TGraphErrors* PhotonEtaPos_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonEtaPos_4060");


	cout << "Loading MC histos...TPC sectors only Nev"<< endl;
	TGraphErrors* PhotonAllSectorsEtaNeg_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonAllSectorsn_0010");
	TGraphErrors* PhotonAllSectorsEtaNeg_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonAllSectorsn_1020");
	TGraphErrors* PhotonAllSectorsEtaNeg_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonAllSectorsn_2040");
	TGraphErrors* PhotonAllSectorsEtaNeg_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonAllSectorsn_4060");
//	TGraphErrors* PhotonAllSectorsEtaNeg6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonAllSectorsn_6080");

	TGraphErrors* PhotonAllSectorsEtaPos_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonAllSectorsp_0010");
	TGraphErrors* PhotonAllSectorsEtaPos_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonAllSectorsp_1020");
	TGraphErrors* PhotonAllSectorsEtaPos_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonAllSectorsp_2040");
	TGraphErrors* PhotonAllSectorsEtaPos_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonAllSectorsp_4060");
//	TGraphErrors* PhotonAllSectorsEtaPos6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonAllSectorsp_6080");

	TGraphErrors* PhotonSector0EtaNeg_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector0n_0010");
	TGraphErrors* PhotonSector0EtaNeg_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector0n_1020");
	TGraphErrors* PhotonSector0EtaNeg_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector0n_2040");
	TGraphErrors* PhotonSector0EtaNeg_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector0n_4060");
//	TGraphErrors* PhotonSector0EtaNeg6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector0n_6080");

	TGraphErrors* PhotonSector0EtaPos_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector0p_0010");
	TGraphErrors* PhotonSector0EtaPos_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector0p_1020");
	TGraphErrors* PhotonSector0EtaPos_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector0p_2040");
	TGraphErrors* PhotonSector0EtaPos_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector0p_4060");
//	TGraphErrors* PhotonSector0EtaPos6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector0p_6080");

	TGraphErrors* PhotonSector1EtaNeg_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector1n_0010");
	TGraphErrors* PhotonSector1EtaNeg_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector1n_1020");	
	TGraphErrors* PhotonSector1EtaNeg_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector1n_2040");
	TGraphErrors* PhotonSector1EtaNeg_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector1n_4060");
//	TGraphErrors* PhotonSector1EtaNeg6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector1n_6080");

	TGraphErrors* PhotonSector1EtaPos_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector1p_0010");
	TGraphErrors* PhotonSector1EtaPos_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector1p_1020");
	TGraphErrors* PhotonSector1EtaPos_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector1p_2040");
	TGraphErrors* PhotonSector1EtaPos_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector1p_4060");
//	TGraphErrors* PhotonSector1EtaPos6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector1p_6080");

	TGraphErrors* PhotonSector2EtaNeg_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector2n_0010");
	TGraphErrors* PhotonSector2EtaNeg_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector2n_1020");
	TGraphErrors* PhotonSector2EtaNeg_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector2n_2040");
	TGraphErrors* PhotonSector2EtaNeg_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector2n_4060");
//	TGraphErrors* PhotonSector2EtaNeg6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector2n_6080");

	TGraphErrors* PhotonSector2EtaPos_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector2p_0010");
	TGraphErrors* PhotonSector2EtaPos_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector2p_1020");
	TGraphErrors* PhotonSector2EtaPos_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector2p_2040");
	TGraphErrors* PhotonSector2EtaPos_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector2p_4060");
//	TGraphErrors* PhotonSector2EtaPos6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector2p_6080");

	TGraphErrors* PhotonSector3EtaNeg_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector3n_0010");
	TGraphErrors* PhotonSector3EtaNeg_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector3n_1020");
	TGraphErrors* PhotonSector3EtaNeg_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector3n_2040");
	TGraphErrors* PhotonSector3EtaNeg_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector3n_4060");
//	TGraphErrors* PhotonSector3EtaNeg6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector3n_6080");

	TGraphErrors* PhotonSector3EtaPos_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector3p_0010");
	TGraphErrors* PhotonSector3EtaPos_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector3p_1020");
	TGraphErrors* PhotonSector3EtaPos_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector3p_2040");
	TGraphErrors* PhotonSector3EtaPos_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector3p_4060");
//	TGraphErrors* PhotonSector3EtaPos6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector3p_6080");

	TGraphErrors* PhotonSector4EtaNeg_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector4n_0010");
	TGraphErrors* PhotonSector4EtaNeg_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector4n_1020");
	TGraphErrors* PhotonSector4EtaNeg_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector4n_2040");
	TGraphErrors* PhotonSector4EtaNeg_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector4n_4060");
//	TGraphErrors* PhotonSector4EtaNeg6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector4n_6080");

	TGraphErrors* PhotonSector4EtaPos_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector4p_0010");
	TGraphErrors* PhotonSector4EtaPos_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector4p_1020");
	TGraphErrors* PhotonSector4EtaPos_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector4p_2040");
	TGraphErrors* PhotonSector4EtaPos_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector4p_4060");
//	TGraphErrors* PhotonSector4EtaPos6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector4p_6080");

	TGraphErrors* PhotonSector5EtaNeg_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector5n_0010");
	TGraphErrors* PhotonSector5EtaNeg_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector5n_1020");
	TGraphErrors* PhotonSector5EtaNeg_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector5n_2040");
	TGraphErrors* PhotonSector5EtaNeg_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector5n_4060");
//	TGraphErrors* PhotonSector5EtaNeg6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector5n_6080");

	TGraphErrors* PhotonSector5EtaPos_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector5p_0010");
	TGraphErrors* PhotonSector5EtaPos_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector5p_1020");
	TGraphErrors* PhotonSector5EtaPos_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector5p_2040");
	TGraphErrors* PhotonSector5EtaPos_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector5p_4060");
//	TGraphErrors* PhotonSector5EtaPos6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector5p_6080");

	TGraphErrors* PhotonSector6EtaNeg_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector6n_0010");
	TGraphErrors* PhotonSector6EtaNeg_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector6n_1020");
	TGraphErrors* PhotonSector6EtaNeg_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector6n_2040");
	TGraphErrors* PhotonSector6EtaNeg_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector6n_4060");
//	TGraphErrors* PhotonSector6EtaNeg6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector6n_6080");

	TGraphErrors* PhotonSector6EtaPos_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector6p_0010");
	TGraphErrors* PhotonSector6EtaPos_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector6p_1020");
	TGraphErrors* PhotonSector6EtaPos_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector6p_2040");
	TGraphErrors* PhotonSector6EtaPos_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector6p_4060");
//	TGraphErrors* PhotonSector6EtaPos6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector6p_6080");

	TGraphErrors* PhotonSector7EtaNeg_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector7n_0010");
	TGraphErrors* PhotonSector7EtaNeg_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector7n_1020");
	TGraphErrors* PhotonSector7EtaNeg_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector7n_2040");
	TGraphErrors* PhotonSector7EtaNeg_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector7n_4060");
//	TGraphErrors* PhotonSector7EtaNeg6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector7n_6080");

	TGraphErrors* PhotonSector7EtaPos_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector7p_0010");
	TGraphErrors* PhotonSector7EtaPos_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector7p_1020");
	TGraphErrors* PhotonSector7EtaPos_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector7p_2040");
	TGraphErrors* PhotonSector7EtaPos_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector7p_4060");
//	TGraphErrors* PhotonSector7EtaPos6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector7p_6080");

	TGraphErrors* PhotonSector8EtaNeg_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector8n_0010");
	TGraphErrors* PhotonSector8EtaNeg_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector8n_1020");
	TGraphErrors* PhotonSector8EtaNeg_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector8n_2040");
	TGraphErrors* PhotonSector8EtaNeg_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector8n_4060");
//	TGraphErrors* PhotonSector8EtaNeg6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector8n_6080");

	TGraphErrors* PhotonSector8EtaPos_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector8p_0010");
	TGraphErrors* PhotonSector8EtaPos_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector8p_1020");
	TGraphErrors* PhotonSector8EtaPos_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector8p_2040");
	TGraphErrors* PhotonSector8EtaPos_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector8p_4060");
//	TGraphErrors* PhotonSector8EtaPos6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector8p_6080");

	TGraphErrors* PhotonSector9EtaNeg_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector9n_0010");
	TGraphErrors* PhotonSector9EtaNeg_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector9n_1020");
	TGraphErrors* PhotonSector9EtaNeg_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector9n_2040");
	TGraphErrors* PhotonSector9EtaNeg_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector9n_4060");
//	TGraphErrors* PhotonSector9EtaNeg6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector9n_6080");

	TGraphErrors* PhotonSector9EtaPos_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector9p_0010");
	TGraphErrors* PhotonSector9EtaPos_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector9p_1020");
	TGraphErrors* PhotonSector9EtaPos_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector9p_2040");
	TGraphErrors* PhotonSector9EtaPos_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector9p_4060");
//	TGraphErrors* PhotonSector9EtaPos6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector9p_6080");

	TGraphErrors* PhotonSector10EtaNeg_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector10n_0010");
	TGraphErrors* PhotonSector10EtaNeg_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector10n_1020");
	TGraphErrors* PhotonSector10EtaNeg_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector10n_2040");
	TGraphErrors* PhotonSector10EtaNeg_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector10n_4060");
//	TGraphErrors* PhotonSector10EtaNeg6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector10n_6080");

	TGraphErrors* PhotonSector10EtaPos_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector10p_0010");
	TGraphErrors* PhotonSector10EtaPos_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector10p_1020");
	TGraphErrors* PhotonSector10EtaPos_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector10p_2040");
	TGraphErrors* PhotonSector10EtaPos_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector10p_4060");
//	TGraphErrors* PhotonSector10EtaPos6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector10p_6080");

	TGraphErrors* PhotonSector11EtaNeg_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector11n_0010");
	TGraphErrors* PhotonSector11EtaNeg_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector11n_1020");
	TGraphErrors* PhotonSector11EtaNeg_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector11n_2040");
	TGraphErrors* PhotonSector11EtaNeg_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector11n_4060");
//	TGraphErrors* PhotonSector11EtaNeg6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector11n_6080");

	TGraphErrors* PhotonSector11EtaPos_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector11p_0010");
	TGraphErrors* PhotonSector11EtaPos_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector11p_1020");
	TGraphErrors* PhotonSector11EtaPos_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector11p_2040");
	TGraphErrors* PhotonSector11EtaPos_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector11p_4060");
//	TGraphErrors* PhotonSector11EtaPos6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector11p_6080");

	TGraphErrors* PhotonSector12EtaNeg_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector12n_0010");
	TGraphErrors* PhotonSector12EtaNeg_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector12n_1020");
	TGraphErrors* PhotonSector12EtaNeg_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector12n_2040");
	TGraphErrors* PhotonSector12EtaNeg_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector12n_4060");
//	TGraphErrors* PhotonSector12EtaNeg6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector12n_6080");

	TGraphErrors* PhotonSector12EtaPos_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector12p_0010");
	TGraphErrors* PhotonSector12EtaPos_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector12p_1020");
	TGraphErrors* PhotonSector12EtaPos_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector12p_2040");
	TGraphErrors* PhotonSector12EtaPos_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector12p_4060");
//	TGraphErrors* PhotonSector12EtaPos6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector12p_6080");

	TGraphErrors* PhotonSector13EtaNeg_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector13n_0010");
	TGraphErrors* PhotonSector13EtaNeg_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector13n_1020");
	TGraphErrors* PhotonSector13EtaNeg_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector13n_2040");
	TGraphErrors* PhotonSector13EtaNeg_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector13n_4060");
//	TGraphErrors* PhotonSector13EtaNeg6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector13n_6080");

	TGraphErrors* PhotonSector13EtaPos_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector13p_0010");
	TGraphErrors* PhotonSector13EtaPos_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector13p_1020");
	TGraphErrors* PhotonSector13EtaPos_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector13p_2040");
	TGraphErrors* PhotonSector13EtaPos_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector13p_4060");
//	TGraphErrors* PhotonSector13EtaPos6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector13p_6080");

	TGraphErrors* PhotonSector14EtaNeg_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector14n_0010");
	TGraphErrors* PhotonSector14EtaNeg_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector14n_1020");
	TGraphErrors* PhotonSector14EtaNeg_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector14n_2040");
	TGraphErrors* PhotonSector14EtaNeg_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector14n_4060");
//	TGraphErrors* PhotonSector14EtaNeg6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector14n_6080");

	TGraphErrors* PhotonSector14EtaPos_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector14p_0010");
	TGraphErrors* PhotonSector14EtaPos_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector14p_1020");
	TGraphErrors* PhotonSector14EtaPos_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector14p_2040");
	TGraphErrors* PhotonSector14EtaPos_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector14p_4060");
//	TGraphErrors* PhotonSector14EtaPos6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector14p_6080");

	TGraphErrors* PhotonSector15EtaNeg_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector15n_0010");
	TGraphErrors* PhotonSector15EtaNeg_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector15n_1020");
	TGraphErrors* PhotonSector15EtaNeg_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector15n_2040");
	TGraphErrors* PhotonSector15EtaNeg_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector15n_4060");
//	TGraphErrors* PhotonSector15EtaNeg6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector15n_6080");

	TGraphErrors* PhotonSector15EtaPos_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector15p_0010");
	TGraphErrors* PhotonSector15EtaPos_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector15p_1020");
	TGraphErrors* PhotonSector15EtaPos_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector15p_2040");
	TGraphErrors* PhotonSector15EtaPos_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector15p_4060");
//	TGraphErrors* PhotonSector15EtaPos6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector15p_6080");

	TGraphErrors* PhotonSector16EtaNeg_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector16n_0010");
	TGraphErrors* PhotonSector16EtaNeg_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector16n_1020");
	TGraphErrors* PhotonSector16EtaNeg_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector16n_2040");
	TGraphErrors* PhotonSector16EtaNeg_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector16n_4060");
//	TGraphErrors* PhotonSector16EtaNeg6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector16n_6080");

	TGraphErrors* PhotonSector16EtaPos_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector16p_0010");
	TGraphErrors* PhotonSector16EtaPos_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector16p_1020");
	TGraphErrors* PhotonSector16EtaPos_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector16p_2040");
	TGraphErrors* PhotonSector16EtaPos_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector16p_4060");
//	TGraphErrors* PhotonSector16EtaPos6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector16p_6080");

	TGraphErrors* PhotonSector17EtaNeg_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector17n_0010");
	TGraphErrors* PhotonSector17EtaNeg_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector17n_1020");
	TGraphErrors* PhotonSector17EtaNeg_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector17n_2040");
	TGraphErrors* PhotonSector17EtaNeg_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector17n_4060");
//	TGraphErrors* PhotonSector17EtaNeg6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector17n_6080");

	TGraphErrors* PhotonSector17EtaPos_0010_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector17p_0010");
	TGraphErrors* PhotonSector17EtaPos_1020_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector17p_1020");
	TGraphErrors* PhotonSector17EtaPos_2040_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector17p_2040");
	TGraphErrors* PhotonSector17EtaPos_4060_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector17p_4060");
//	TGraphErrors* PhotonSector17EtaPos6080_MC = (TGraphErrors*)QAOutputFileMC->Get("PhotonSector17p_6080");


	cout << "Drawing MC" << endl;


	DrawGammaSetMarkerTGraphErr( EventsPerRun_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr( EventsPerRun_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr( EventsPerRun_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr( EventsPerRun_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);

	//nsigma
	DrawGammaSetMarkerTGraphErr( NSigmadEdxElectron_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxElectron_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxElectron_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxElectron_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);

	DrawGammaSetMarkerTGraphErr( NSigmadEdxPositron_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxPositron_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxPositron_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxPositron_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);

	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxElectron_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxElectron_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxElectron_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxElectron_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);

	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxPositron_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxPositron_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxPositron_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxPositron_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);

	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxElectron_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxElectron_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxElectron_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxElectron_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);

	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxPositron_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxPositron_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxPositron_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxPositron_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);


	DrawGammaSetMarkerTGraphErr( NSigmadEdxElectronProtonCut_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxElectronProtonCut_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxElectronProtonCut_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxElectronProtonCut_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);

	DrawGammaSetMarkerTGraphErr( NSigmadEdxPositronProtonCut_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxPositronProtonCut_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxPositronProtonCut_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr( NSigmadEdxPositronProtonCut_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);

	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxElectronProtonCut_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxElectronProtonCut_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxElectronProtonCut_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxElectronProtonCut_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);

	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxPositronProtonCut_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxPositronProtonCut_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxPositronProtonCut_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr( sigmaNSigmadEdxPositronProtonCut_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);

	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxElectronProtonCut_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxElectronProtonCut_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxElectronProtonCut_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxElectronProtonCut_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);

	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxPositronProtonCut_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxPositronProtonCut_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxPositronProtonCut_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr( widthNSigmadEdxPositronProtonCut_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);

	//photon pt 
	DrawGammaSetMarkerTGraphErr( PhotonPt_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr( PhotonPt_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr( PhotonPt_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr( PhotonPt_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);


	//eta
	DrawGammaSetMarkerTGraphErr( ElectronEta_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr( ElectronEta_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr( ElectronEta_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr( ElectronEta_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);

	DrawGammaSetMarkerTGraphErr( PositronEta_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr( PositronEta_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr( PositronEta_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr( PositronEta_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);

	DrawGammaSetMarkerTGraphErr( PhotonEta_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr( PhotonEta_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr( PhotonEta_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr( PhotonEta_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);

	//etaneg
	DrawGammaSetMarkerTGraphErr( ElectronEtaNeg_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr( ElectronEtaNeg_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr( ElectronEtaNeg_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr( ElectronEtaNeg_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);

	DrawGammaSetMarkerTGraphErr( PositronEtaNeg_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr( PositronEtaNeg_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr( PositronEtaNeg_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr( PositronEtaNeg_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);

	DrawGammaSetMarkerTGraphErr( PhotonEtaNeg_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr( PhotonEtaNeg_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr( PhotonEtaNeg_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr( PhotonEtaNeg_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);

	//etapos
	DrawGammaSetMarkerTGraphErr( ElectronEtaPos_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr( ElectronEtaPos_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr( ElectronEtaPos_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr( ElectronEtaPos_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);

	DrawGammaSetMarkerTGraphErr( PositronEtaPos_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr( PositronEtaPos_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr( PositronEtaPos_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr( PositronEtaPos_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);

	DrawGammaSetMarkerTGraphErr( PhotonEtaPos_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr( PhotonEtaPos_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr( PhotonEtaPos_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr( PhotonEtaPos_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);


	DrawGammaSetMarkerTGraphErr(PhotonAllSectorsEtaNeg_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonAllSectorsEtaNeg_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonAllSectorsEtaNeg_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonAllSectorsEtaNeg_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonAllSectorsEtaNeg6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonAllSectorsEtaPos_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonAllSectorsEtaPos_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonAllSectorsEtaPos_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonAllSectorsEtaPos_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonAllSectorsEtaPos6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector0EtaNeg_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector0EtaNeg_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector0EtaNeg_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector0EtaNeg_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector0EtaNeg6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector0EtaPos_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector0EtaPos_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector0EtaPos_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector0EtaPos_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector0EtaPos6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);  

	DrawGammaSetMarkerTGraphErr(PhotonSector1EtaNeg_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector1EtaNeg_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);	
	DrawGammaSetMarkerTGraphErr(PhotonSector1EtaNeg_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector1EtaNeg_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector1EtaNeg6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector1EtaPos_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector1EtaPos_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector1EtaPos_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector1EtaPos_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector1EtaPos6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector2EtaNeg_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector2EtaNeg_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector2EtaNeg_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector2EtaNeg_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector2EtaNeg6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector2EtaPos_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector2EtaPos_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector2EtaPos_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector2EtaPos_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector2EtaPos6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector3EtaNeg_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector3EtaNeg_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector3EtaNeg_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector3EtaNeg_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector3EtaNeg6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector3EtaPos_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector3EtaPos_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector3EtaPos_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector3EtaPos_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector3EtaPos6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector4EtaNeg_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector4EtaNeg_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector4EtaNeg_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector4EtaNeg_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector4EtaNeg6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector4EtaPos_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector4EtaPos_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector4EtaPos_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector4EtaPos_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector4EtaPos6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector5EtaNeg_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector5EtaNeg_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector5EtaNeg_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector5EtaNeg_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector5EtaNeg6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector5EtaPos_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector5EtaPos_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector5EtaPos_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector5EtaPos_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector5EtaPos6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector6EtaNeg_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector6EtaNeg_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector6EtaNeg_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector6EtaNeg_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector6EtaNeg6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector6EtaPos_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector6EtaPos_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector6EtaPos_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector6EtaPos_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector6EtaPos6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector7EtaNeg_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector7EtaNeg_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector7EtaNeg_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector7EtaNeg_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector7EtaNeg6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector7EtaPos_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector7EtaPos_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector7EtaPos_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector7EtaPos_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector7EtaPos6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector8EtaNeg_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector8EtaNeg_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector8EtaNeg_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector8EtaNeg_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector8EtaNeg6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector8EtaPos_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector8EtaPos_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector8EtaPos_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector8EtaPos_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector8EtaPos6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector9EtaNeg_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector9EtaNeg_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector9EtaNeg_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector9EtaNeg_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector9EtaNeg6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector9EtaPos_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector9EtaPos_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector9EtaPos_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector9EtaPos_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector9EtaPos6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector10EtaNeg_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector10EtaNeg_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector10EtaNeg_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector10EtaNeg_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector10EtaNeg6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector10EtaPos_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector10EtaPos_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector10EtaPos_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector10EtaPos_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector10EtaPos6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector11EtaNeg_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector11EtaNeg_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector11EtaNeg_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector11EtaNeg_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector11EtaNeg6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector11EtaPos_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector11EtaPos_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector11EtaPos_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector11EtaPos_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
///	DrawGammaSetMarkerTGraphErr(PhotonSector11EtaPos6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector12EtaNeg_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector12EtaNeg_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector12EtaNeg_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector12EtaNeg_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector12EtaNeg6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector12EtaPos_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector12EtaPos_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector12EtaPos_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector12EtaPos_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector12EtaPos6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector13EtaNeg_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector13EtaNeg_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector13EtaNeg_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector13EtaNeg_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector13EtaNeg6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector13EtaPos_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector13EtaPos_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector13EtaPos_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector13EtaPos_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector13EtaPos6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector14EtaNeg_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector14EtaNeg_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector14EtaNeg_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector14EtaNeg_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector14EtaNeg6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector14EtaPos_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector14EtaPos_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector14EtaPos_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector14EtaPos_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector14EtaPos6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector15EtaNeg_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector15EtaNeg_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector15EtaNeg_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector15EtaNeg_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector15EtaNeg6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector15EtaPos_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector15EtaPos_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector15EtaPos_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector15EtaPos_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector15EtaPos6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector16EtaNeg_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector16EtaNeg_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector16EtaNeg_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector16EtaNeg_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector16EtaNeg6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector16EtaPos_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector16EtaPos_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector16EtaPos_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector16EtaPos_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector16EtaPos6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector17EtaNeg_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector17EtaNeg_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector17EtaNeg_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector17EtaNeg_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector17EtaNeg6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);

	DrawGammaSetMarkerTGraphErr(PhotonSector17EtaPos_0010_MC,  markerStyleMC_0010_MC, markerSize_0010_MC, colorComb_0010_MC, colorComb_0010_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector17EtaPos_1020_MC,  markerStyleMC_1020_MC, markerSize_1020_MC, colorComb_1020_MC, colorComb_1020_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector17EtaPos_2040_MC,  markerStyleMC_2040_MC, markerSize_2040_MC, colorComb_2040_MC, colorComb_2040_MC);
	DrawGammaSetMarkerTGraphErr(PhotonSector17EtaPos_4060_MC,  markerStyleMC_4060_MC, markerSize_4060_MC, colorComb_4060_MC, colorComb_4060_MC);
//	DrawGammaSetMarkerTGraphErr(PhotonSector17EtaPos6080_MC,  markerStyleMC6080_MC, markerSize6080_MC, colorComb6080_MC, colorComb6080_MC);


   Float_t floatLocationLeft[4] = {0.7,0.9,0.04,0.02};
   TString collisionSystem = "Pb-Pb, #sqrt{s_{NN}} = 2.76 TeV";
   TString textPeriod = "Data LHC11h,";
   TString optMCGenerator = "MC Hijing";

	TPaveText *MC = new TPaveText(0.1,0.925,0.2,0.995,"NDC");
	MC->SetTextSize(0.04);
	MC->SetLineColor(0);
	MC->SetFillColor(0);
	MC->SetFillStyle(0);
	MC->SetShadowColor(0);
	MC->SetTextAlign(12);
	MC->AddText("Data");
	MC->AddText("MC");

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////// 	            Drawing part		/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if(partialOutput){

		cout << "Plot N sigma dE/dx" << endl;
		TCanvas* canvasNsigmaElectron = new TCanvas("canvasNsigmaElectron","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasNsigmaElectron, 0.1, 0.03, 0.08, 0.08);
		//canvasNsigmaElectron->SetLogy();
		TH2F * histo2DnNSigmaE;
		histo2DnNSigmaE = new TH2F("histo2DnNSigmaE","histo2DnNSigmaE",2700,167900,170600,10000,-10,10);
		histo2DnNSigmaE->GetYaxis()->SetRangeUser(-.6,1.2);
		SetStyleHistoTH2ForGraphs(histo2DnNSigmaE, "Run number","<n#sigma dE^{e^{-}}/dx>/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histo2DnNSigmaE->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		NSigmadEdxElectron_0010_Data->Draw("pesame");
		NSigmadEdxElectron_1020_Data->Draw("pesame");
		NSigmadEdxElectron_2040_Data->Draw("pesame");
		NSigmadEdxElectron_4060_Data->Draw("pesame");
	// 	NSigmadEdxElectron6080_Data->Draw("pesame");
		NSigmadEdxElectron_0010_MC->Draw("pesame");
		NSigmadEdxElectron_1020_MC->Draw("pesame");
		NSigmadEdxElectron_2040_MC->Draw("pesame");
		NSigmadEdxElectron_4060_MC->Draw("pesame");
	// 	NSigmadEdxElectron6080_MC->Draw("pesame");

		TLegend* legendNSE = new TLegend(0.22,0.93,0.99,.99);
		legendNSE->SetFillColor(0);
		legendNSE->SetLineColor(0);
		legendNSE->SetTextSize(0.04);
		legendNSE->SetNColumns(4);
	// 	legendNSE->AddEntry((TObject*)0, "Data:","");
		legendNSE->AddEntry(NSigmadEdxElectron_0010_Data,collisionSystem0010.Data(),"p");
		legendNSE->AddEntry(NSigmadEdxElectron_1020_Data,collisionSystem1020.Data(),"p");
		legendNSE->AddEntry(NSigmadEdxElectron_2040_Data,collisionSystem2040.Data(),"p");
		legendNSE->AddEntry(NSigmadEdxElectron_4060_Data,collisionSystem4060.Data(),"p");
	//	legendNSE->AddEntry(NSigmadEdxElectron6080_Data,collisionSystem6080.Data(),"p");
		legendNSE->AddEntry(NSigmadEdxElectron_0010_MC,collisionSystem0010.Data(),"p");
		legendNSE->AddEntry(NSigmadEdxElectron_1020_MC,collisionSystem1020.Data(),"p");
		legendNSE->AddEntry(NSigmadEdxElectron_2040_MC,collisionSystem2040.Data(),"p");
		legendNSE->AddEntry(NSigmadEdxElectron_4060_MC,collisionSystem4060.Data(),"p");
	//	legendNSE->AddEntry(NSigmadEdxElectron6080_MC,collisionSystem6080.Data(),"p");
		MC->Draw();
		legendNSE->Draw();
		canvasNsigmaElectron->SaveAs(Form("%s/NSigmaElectron.%s",outputDir.Data(),suffix.Data()));
		delete canvasNsigmaElectron;


		TCanvas* canvasNsigmaPositron = new TCanvas("canvasNsigmaPositron","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasNsigmaPositron, 0.1, 0.03, 0.08, 0.08);
		//canvasNsigmaPositron->SetLogy();
		TH2F * histo2DnNSigmaP;
		histo2DnNSigmaP = new TH2F("histo2DnNSigmaP","histo2DnNSigmaP",2700,167900,170600,10000,-10,10);
		histo2DnNSigmaP->GetYaxis()->SetRangeUser(-.6,1.2);
		SetStyleHistoTH2ForGraphs(histo2DnNSigmaP, "Run number","<n#sigma dE^{e^{+}}/dx>/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histo2DnNSigmaP->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		NSigmadEdxPositron_0010_Data->Draw("pesame");
		NSigmadEdxPositron_1020_Data->Draw("pesame");
		NSigmadEdxPositron_2040_Data->Draw("pesame");
		NSigmadEdxPositron_4060_Data->Draw("pesame");
	// 	NSigmadEdxPositron6080_Data->Draw("pesame");
		NSigmadEdxPositron_0010_MC->Draw("pesame");
		NSigmadEdxPositron_1020_MC->Draw("pesame");
		NSigmadEdxPositron_2040_MC->Draw("pesame");
		NSigmadEdxPositron_4060_MC->Draw("pesame");
	// 	NSigmadEdxPositron6080_MC->Draw("pesame");
		TLegend* legendNSP = new TLegend(0.22,0.93,0.99,.99);
		legendNSP->SetFillColor(0);
		legendNSP->SetLineColor(0);
		legendNSP->SetTextSize(0.04);
		legendNSP->SetNColumns(4);
	// 	legendNSP->AddEntry((TObject*)0, "Data:","");
		legendNSP->AddEntry(NSigmadEdxPositron_0010_Data,collisionSystem0010.Data(),"p");
		legendNSP->AddEntry(NSigmadEdxPositron_1020_Data,collisionSystem1020.Data(),"p");
		legendNSP->AddEntry(NSigmadEdxPositron_2040_Data,collisionSystem2040.Data(),"p");
		legendNSP->AddEntry(NSigmadEdxPositron_4060_Data,collisionSystem4060.Data(),"p");
	//	legendNSP->AddEntry(NSigmadEdxPositron6080_Data,collisionSystem6080.Data(),"p");
		legendNSP->AddEntry(NSigmadEdxPositron_0010_MC,collisionSystem0010.Data(),"p");
		legendNSP->AddEntry(NSigmadEdxPositron_1020_MC,collisionSystem1020.Data(),"p");
		legendNSP->AddEntry(NSigmadEdxPositron_2040_MC,collisionSystem2040.Data(),"p");
		legendNSP->AddEntry(NSigmadEdxPositron_4060_MC,collisionSystem4060.Data(),"p");
	//	legendNSP->AddEntry(NSigmadEdxPositron6080_MC,collisionSystem6080.Data(),"p");
		MC->Draw();
		legendNSP->Draw();
		canvasNsigmaPositron->SaveAs(Form("%s/NSigmaPositron.%s",outputDir.Data(),suffix.Data()));
		delete canvasNsigmaPositron;

		cout << "Plot sigma Nsigma dE/dx" << endl;
		TCanvas* sigmaEL = new TCanvas("sigmaEL","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(sigmaEL, 0.1, 0.03, 0.08, 0.08);
		TH2F *hsigmaEL;
		hsigmaEL = new TH2F("hsigmaEL","hsigmaEL",2700,167900,170600,10000,-10,10);
		hsigmaEL->GetYaxis()->SetRangeUser(0.6,1.9);
		SetStyleHistoTH2ForGraphs(hsigmaEL, "Run number","#sigma of n#sigma dE^{e^{-}}/dx/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		hsigmaEL->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		sigmaNSigmadEdxElectron_0010_Data->Draw("pesame");
		sigmaNSigmadEdxElectron_1020_Data->Draw("pesame");
		sigmaNSigmadEdxElectron_2040_Data->Draw("pesame");
		sigmaNSigmadEdxElectron_4060_Data->Draw("pesame");
		sigmaNSigmadEdxElectron_0010_MC->Draw("pesame");
		sigmaNSigmadEdxElectron_1020_MC->Draw("pesame");
		sigmaNSigmadEdxElectron_2040_MC->Draw("pesame");
		sigmaNSigmadEdxElectron_4060_MC->Draw("pesame");

		TLegend* sigmaLEL = new TLegend(0.22,0.93,0.99,.99);
		sigmaLEL->SetFillColor(0);
		sigmaLEL->SetLineColor(0);
		sigmaLEL->SetTextSize(0.04);
		sigmaLEL->SetNColumns(4);
		sigmaLEL->AddEntry(sigmaNSigmadEdxElectron_0010_Data,collisionSystem0010.Data(),"p");
		sigmaLEL->AddEntry(sigmaNSigmadEdxElectron_1020_Data,collisionSystem1020.Data(),"p");
		sigmaLEL->AddEntry(sigmaNSigmadEdxElectron_2040_Data,collisionSystem2040.Data(),"p");
		sigmaLEL->AddEntry(sigmaNSigmadEdxElectron_4060_Data,collisionSystem4060.Data(),"p");
		sigmaLEL->AddEntry(sigmaNSigmadEdxElectron_0010_MC,collisionSystem0010.Data(),"p");
		sigmaLEL->AddEntry(sigmaNSigmadEdxElectron_1020_MC,collisionSystem1020.Data(),"p");
		sigmaLEL->AddEntry(sigmaNSigmadEdxElectron_2040_MC,collisionSystem2040.Data(),"p");
		sigmaLEL->AddEntry(sigmaNSigmadEdxElectron_4060_MC,collisionSystem4060.Data(),"p");

		MC->Draw();
		sigmaLEL->Draw();
		sigmaEL->SaveAs(Form("%s/sigmaNSigmaElectron.%s",outputDir.Data(),suffix.Data()));
		delete sigmaEL;

		cout << "Plot sigma Nsigma dE/dx" << endl;
		TCanvas* sigmaPS = new TCanvas("sigmaPS","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(sigmaPS, 0.1, 0.03, 0.08, 0.08);
		TH2F *hsigmaPS;
		hsigmaPS = new TH2F("hsigmaPS","hsigmaPS",2700,167900,170600,10000,-10,10);
		hsigmaPS->GetYaxis()->SetRangeUser(0.6,1.9);
		SetStyleHistoTH2ForGraphs(hsigmaPS, "Run number","#sigma of n#sigma dE^{e^{-}}/dx/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		hsigmaPS->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		sigmaNSigmadEdxPositron_0010_Data->Draw("pesame");
		sigmaNSigmadEdxPositron_1020_Data->Draw("pesame");
		sigmaNSigmadEdxPositron_2040_Data->Draw("pesame");
		sigmaNSigmadEdxPositron_4060_Data->Draw("pesame");
		sigmaNSigmadEdxPositron_0010_MC->Draw("pesame");
		sigmaNSigmadEdxPositron_1020_MC->Draw("pesame");
		sigmaNSigmadEdxPositron_2040_MC->Draw("pesame");
		sigmaNSigmadEdxPositron_4060_MC->Draw("pesame");

		TLegend* sigmaLPS = new TLegend(0.22,0.93,0.99,.99);
		sigmaLPS->SetFillColor(0);
		sigmaLPS->SetLineColor(0);
		sigmaLPS->SetTextSize(0.04);
		sigmaLPS->SetNColumns(4);
		sigmaLPS->AddEntry(sigmaNSigmadEdxPositron_0010_Data,collisionSystem0010.Data(),"p");
		sigmaLPS->AddEntry(sigmaNSigmadEdxPositron_1020_Data,collisionSystem1020.Data(),"p");
		sigmaLPS->AddEntry(sigmaNSigmadEdxPositron_2040_Data,collisionSystem2040.Data(),"p");
		sigmaLPS->AddEntry(sigmaNSigmadEdxPositron_4060_Data,collisionSystem4060.Data(),"p");
		sigmaLPS->AddEntry(sigmaNSigmadEdxPositron_0010_MC,collisionSystem0010.Data(),"p");
		sigmaLPS->AddEntry(sigmaNSigmadEdxPositron_1020_MC,collisionSystem1020.Data(),"p");
		sigmaLPS->AddEntry(sigmaNSigmadEdxPositron_2040_MC,collisionSystem2040.Data(),"p");
		sigmaLPS->AddEntry(sigmaNSigmadEdxPositron_4060_MC,collisionSystem4060.Data(),"p");

		MC->Draw();
		sigmaLPS->Draw();
		sigmaPS->SaveAs(Form("%s/sigmaNSigmaPositron.%s",outputDir.Data(),suffix.Data()));
		delete sigmaPS;


		cout << "Plot width Nwidth dE/dx" << endl;
		TCanvas* widthEL = new TCanvas("widthEL","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(widthEL, 0.1, 0.03, 0.08, 0.08);
		TH2F *hwidthEL;
		hwidthEL = new TH2F("hwidthEL","hwidthEL",2700,167900,170600,10000,-10,10);
		hwidthEL->GetYaxis()->SetRangeUser(0.,5.);
		SetStyleHistoTH2ForGraphs(hwidthEL, "Run number","width of n#sigma dE^{e^{-}}/dx/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		hwidthEL->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		widthNSigmadEdxElectron_0010_Data->Draw("pesame");
		widthNSigmadEdxElectron_1020_Data->Draw("pesame");
		widthNSigmadEdxElectron_2040_Data->Draw("pesame");
		widthNSigmadEdxElectron_4060_Data->Draw("pesame");
		widthNSigmadEdxElectron_0010_MC->Draw("pesame");
		widthNSigmadEdxElectron_1020_MC->Draw("pesame");
		widthNSigmadEdxElectron_2040_MC->Draw("pesame");
		widthNSigmadEdxElectron_4060_MC->Draw("pesame");

		TLegend* widthLEL = new TLegend(0.22,0.93,0.99,.99);
		widthLEL->SetFillColor(0);
		widthLEL->SetLineColor(0);
		widthLEL->SetTextSize(0.04);
		widthLEL->SetNColumns(4);
		widthLEL->AddEntry(widthNSigmadEdxElectron_0010_Data,collisionSystem0010.Data(),"p");
		widthLEL->AddEntry(widthNSigmadEdxElectron_1020_Data,collisionSystem1020.Data(),"p");
		widthLEL->AddEntry(widthNSigmadEdxElectron_2040_Data,collisionSystem2040.Data(),"p");
		widthLEL->AddEntry(widthNSigmadEdxElectron_4060_Data,collisionSystem4060.Data(),"p");
		widthLEL->AddEntry(widthNSigmadEdxElectron_0010_MC,collisionSystem0010.Data(),"p");
		widthLEL->AddEntry(widthNSigmadEdxElectron_1020_MC,collisionSystem1020.Data(),"p");
		widthLEL->AddEntry(widthNSigmadEdxElectron_2040_MC,collisionSystem2040.Data(),"p");
		widthLEL->AddEntry(widthNSigmadEdxElectron_4060_MC,collisionSystem4060.Data(),"p");

		MC->Draw();
		widthLEL->Draw();
		widthEL->SaveAs(Form("%s/widthNSigmaElectron.%s",outputDir.Data(),suffix.Data()));
		delete widthEL;

		cout << "Plot width Nwidth dE/dx" << endl;
		TCanvas* widthPS = new TCanvas("widthPS","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(widthPS, 0.1, 0.03, 0.08, 0.08);
		TH2F *hwidthPS;
		hwidthPS = new TH2F("hwidthPS","hwidthPS",2700,167900,170600,10000,-10,10);
		hwidthPS->GetYaxis()->SetRangeUser(0.,5.);
		SetStyleHistoTH2ForGraphs(hwidthPS, "Run number","width of n#sigma dE^{e^{-}}/dx/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		hwidthPS->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		widthNSigmadEdxPositron_0010_Data->Draw("pesame");
		widthNSigmadEdxPositron_1020_Data->Draw("pesame");
		widthNSigmadEdxPositron_2040_Data->Draw("pesame");
		widthNSigmadEdxPositron_4060_Data->Draw("pesame");
		widthNSigmadEdxPositron_0010_MC->Draw("pesame");
		widthNSigmadEdxPositron_1020_MC->Draw("pesame");
		widthNSigmadEdxPositron_2040_MC->Draw("pesame");
		widthNSigmadEdxPositron_4060_MC->Draw("pesame");

		TLegend* widthLPS = new TLegend(0.22,0.93,0.99,.99);
		widthLPS->SetFillColor(0);
		widthLPS->SetLineColor(0);
		widthLPS->SetTextSize(0.04);
		widthLPS->SetNColumns(4);
		widthLPS->AddEntry(widthNSigmadEdxPositron_0010_Data,collisionSystem0010.Data(),"p");
		widthLPS->AddEntry(widthNSigmadEdxPositron_1020_Data,collisionSystem1020.Data(),"p");
		widthLPS->AddEntry(widthNSigmadEdxPositron_2040_Data,collisionSystem2040.Data(),"p");
		widthLPS->AddEntry(widthNSigmadEdxPositron_4060_Data,collisionSystem4060.Data(),"p");
		widthLPS->AddEntry(widthNSigmadEdxPositron_0010_MC,collisionSystem0010.Data(),"p");
		widthLPS->AddEntry(widthNSigmadEdxPositron_1020_MC,collisionSystem1020.Data(),"p");
		widthLPS->AddEntry(widthNSigmadEdxPositron_2040_MC,collisionSystem2040.Data(),"p");
		widthLPS->AddEntry(widthNSigmadEdxPositron_4060_MC,collisionSystem4060.Data(),"p");

		MC->Draw();
		widthLPS->Draw();
		widthPS->SaveAs(Form("%s/widthNSigmaPositron.%s",outputDir.Data(),suffix.Data()));
		delete widthPS;



		////with the proton cut
		cout << "Plot N sigma dE/dx" << endl;
		TCanvas* canvasNsigmaElectronProtonCut = new TCanvas("canvasNsigmaElectronProtonCut","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasNsigmaElectronProtonCut, 0.1, 0.03, 0.08, 0.08);
		TH2F * histo2DnNSigmaEProtonCut;
		histo2DnNSigmaEProtonCut = new TH2F("histo2DnNSigmaEProtonCut","histo2DnNSigmaEProtonCut",2700,167900,170600,10000,-10,10);
		histo2DnNSigmaEProtonCut->GetYaxis()->SetRangeUser(-.6,1.2);
		SetStyleHistoTH2ForGraphs(histo2DnNSigmaEProtonCut, "Run number","<n#sigma dE^{e^{-}}/dx>/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histo2DnNSigmaEProtonCut->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		NSigmadEdxElectronProtonCut_0010_Data->Draw("pesame");
		NSigmadEdxElectronProtonCut_1020_Data->Draw("pesame");
		NSigmadEdxElectronProtonCut_2040_Data->Draw("pesame");
		NSigmadEdxElectronProtonCut_4060_Data->Draw("pesame");
		NSigmadEdxElectronProtonCut_0010_MC->Draw("pesame");
		NSigmadEdxElectronProtonCut_1020_MC->Draw("pesame");
		NSigmadEdxElectronProtonCut_2040_MC->Draw("pesame");
		NSigmadEdxElectronProtonCut_4060_MC->Draw("pesame");
	// 	NSigmadEdxElectron6080_MC->Draw("pesame");

		TLegend* legendNSEProtonCut = new TLegend(0.22,0.93,0.99,.99);
		legendNSEProtonCut->SetFillColor(0);
		legendNSEProtonCut->SetLineColor(0);
		legendNSEProtonCut->SetTextSize(0.04);
		legendNSEProtonCut->SetNColumns(4);
	// 	legendNSEProtonCut->AddEntry((TObject*)0, "Data:","");
		legendNSEProtonCut->AddEntry(NSigmadEdxElectronProtonCut_0010_Data,collisionSystem0010.Data(),"p");
		legendNSEProtonCut->AddEntry(NSigmadEdxElectronProtonCut_1020_Data,collisionSystem1020.Data(),"p");
		legendNSEProtonCut->AddEntry(NSigmadEdxElectronProtonCut_2040_Data,collisionSystem2040.Data(),"p");
		legendNSEProtonCut->AddEntry(NSigmadEdxElectronProtonCut_4060_Data,collisionSystem4060.Data(),"p");
	//	legendNSEProtonCut->AddEntry(NSigmadEdxElectron6080_Data,collisionSystem6080.Data(),"p");
		legendNSEProtonCut->AddEntry(NSigmadEdxElectronProtonCut_0010_MC,collisionSystem0010.Data(),"p");
		legendNSEProtonCut->AddEntry(NSigmadEdxElectronProtonCut_1020_MC,collisionSystem1020.Data(),"p");
		legendNSEProtonCut->AddEntry(NSigmadEdxElectronProtonCut_2040_MC,collisionSystem2040.Data(),"p");
		legendNSEProtonCut->AddEntry(NSigmadEdxElectronProtonCut_4060_MC,collisionSystem4060.Data(),"p");
	//	legendNSEProtonCut->AddEntry(NSigmadEdxElectron6080_MC,collisionSystem6080.Data(),"p");
		MC->Draw();
		legendNSEProtonCut->Draw();
		canvasNsigmaElectronProtonCut->SaveAs(Form("%s/NSigmaElectronProtonCut.%s",outputDir.Data(),suffix.Data()));
		delete canvasNsigmaElectronProtonCut;


		TCanvas* canvasNsigmaPositronProtonCut = new TCanvas("canvasNsigmaPositronProtonCut","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasNsigmaPositronProtonCut, 0.1, 0.03, 0.08, 0.08);
		//canvasNsigmaPositronProtonCut->SetLogy();
		TH2F * histo2DnNSigmaPProtonCut;
		histo2DnNSigmaPProtonCut = new TH2F("histo2DnNSigmaPProtonCut","histo2DnNSigmaPProtonCut",2700,167900,170600,10000,-10,10);
		histo2DnNSigmaPProtonCut->GetYaxis()->SetRangeUser(-.6,1.2);
		SetStyleHistoTH2ForGraphs(histo2DnNSigmaPProtonCut, "Run number","<n#sigma dE^{e^{+}}/dx>/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histo2DnNSigmaPProtonCut->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		NSigmadEdxPositronProtonCut_0010_Data->Draw("pesame");
		NSigmadEdxPositronProtonCut_1020_Data->Draw("pesame");
		NSigmadEdxPositronProtonCut_2040_Data->Draw("pesame");
		NSigmadEdxPositronProtonCut_4060_Data->Draw("pesame");
		NSigmadEdxPositronProtonCut_0010_MC->Draw("pesame");
		NSigmadEdxPositronProtonCut_1020_MC->Draw("pesame");
		NSigmadEdxPositronProtonCut_2040_MC->Draw("pesame");
		NSigmadEdxPositronProtonCut_4060_MC->Draw("pesame");
		TLegend* legendNSPProtonCut = new TLegend(0.22,0.93,0.99,.99);
		legendNSPProtonCut->SetFillColor(0);
		legendNSPProtonCut->SetLineColor(0);
		legendNSPProtonCut->SetTextSize(0.04);
		legendNSPProtonCut->SetNColumns(4);
	// 	legendNSPProtonCut->AddEntry((TObject*)0, "Data:","");
		legendNSPProtonCut->AddEntry(NSigmadEdxPositronProtonCut_0010_Data,collisionSystem0010.Data(),"p");
		legendNSPProtonCut->AddEntry(NSigmadEdxPositronProtonCut_1020_Data,collisionSystem1020.Data(),"p");
		legendNSPProtonCut->AddEntry(NSigmadEdxPositronProtonCut_2040_Data,collisionSystem2040.Data(),"p");
		legendNSPProtonCut->AddEntry(NSigmadEdxPositronProtonCut_4060_Data,collisionSystem4060.Data(),"p");
	//	legendNSPProtonCut->AddEntry(NSigmadEdxPositron6080_Data,collisionSystem6080.Data(),"p");
		legendNSPProtonCut->AddEntry(NSigmadEdxPositronProtonCut_0010_MC,collisionSystem0010.Data(),"p");
		legendNSPProtonCut->AddEntry(NSigmadEdxPositronProtonCut_1020_MC,collisionSystem1020.Data(),"p");
		legendNSPProtonCut->AddEntry(NSigmadEdxPositronProtonCut_2040_MC,collisionSystem2040.Data(),"p");
		legendNSPProtonCut->AddEntry(NSigmadEdxPositronProtonCut_4060_MC,collisionSystem4060.Data(),"p");
	//	legendNSPProtonCut->AddEntry(NSigmadEdxPositron6080_MC,collisionSystem6080.Data(),"p");
		MC->Draw();
		legendNSPProtonCut->Draw();
		canvasNsigmaPositronProtonCut->SaveAs(Form("%s/NSigmaPositronProtonCut.%s",outputDir.Data(),suffix.Data()));
		delete canvasNsigmaPositronProtonCut;

		cout << "Plot sigma Nsigma dE/dx" << endl;
		TCanvas* sigmaELProtonCut = new TCanvas("sigmaELProtonCut","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(sigmaELProtonCut, 0.1, 0.03, 0.08, 0.08);
		TH2F *hsigmaELProtonCut;
		hsigmaELProtonCut = new TH2F("hsigmaELProtonCut","hsigmaELProtonCut",2700,167900,170600,10000,-10,10);
		hsigmaELProtonCut->GetYaxis()->SetRangeUser(0.6,1.9);
		SetStyleHistoTH2ForGraphs(hsigmaELProtonCut, "Run number","#sigma of n#sigma dE^{e^{-}}/dx/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		hsigmaELProtonCut->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		sigmaNSigmadEdxElectronProtonCut_0010_Data->Draw("pesame");
		sigmaNSigmadEdxElectronProtonCut_1020_Data->Draw("pesame");
		sigmaNSigmadEdxElectronProtonCut_2040_Data->Draw("pesame");
		sigmaNSigmadEdxElectronProtonCut_4060_Data->Draw("pesame");
		sigmaNSigmadEdxElectronProtonCut_0010_MC->Draw("pesame");
		sigmaNSigmadEdxElectronProtonCut_1020_MC->Draw("pesame");
		sigmaNSigmadEdxElectronProtonCut_2040_MC->Draw("pesame");
		sigmaNSigmadEdxElectronProtonCut_4060_MC->Draw("pesame");

		TLegend* sigmaLELProtonCut = new TLegend(0.22,0.93,0.99,.99);
		sigmaLELProtonCut->SetFillColor(0);
		sigmaLELProtonCut->SetLineColor(0);
		sigmaLELProtonCut->SetTextSize(0.04);
		sigmaLELProtonCut->SetNColumns(4);
		sigmaLELProtonCut->AddEntry(sigmaNSigmadEdxElectronProtonCut_0010_Data,collisionSystem0010.Data(),"p");
		sigmaLELProtonCut->AddEntry(sigmaNSigmadEdxElectronProtonCut_1020_Data,collisionSystem1020.Data(),"p");
		sigmaLELProtonCut->AddEntry(sigmaNSigmadEdxElectronProtonCut_2040_Data,collisionSystem2040.Data(),"p");
		sigmaLELProtonCut->AddEntry(sigmaNSigmadEdxElectronProtonCut_4060_Data,collisionSystem4060.Data(),"p");
		sigmaLELProtonCut->AddEntry(sigmaNSigmadEdxElectronProtonCut_0010_MC,collisionSystem0010.Data(),"p");
		sigmaLELProtonCut->AddEntry(sigmaNSigmadEdxElectronProtonCut_1020_MC,collisionSystem1020.Data(),"p");
		sigmaLELProtonCut->AddEntry(sigmaNSigmadEdxElectronProtonCut_2040_MC,collisionSystem2040.Data(),"p");
		sigmaLELProtonCut->AddEntry(sigmaNSigmadEdxElectronProtonCut_4060_MC,collisionSystem4060.Data(),"p");

		MC->Draw();
		sigmaLELProtonCut->Draw();
		sigmaELProtonCut->SaveAs(Form("%s/sigmaNSigmaElectronProtonCut.%s",outputDir.Data(),suffix.Data()));
		delete sigmaELProtonCut;

		cout << "Plot sigma Nsigma dE/dx" << endl;
		TCanvas* sigmaPSProtonCut = new TCanvas("sigmaPSProtonCut","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(sigmaPSProtonCut, 0.1, 0.03, 0.08, 0.08);
		TH2F *hsigmaPSProtonCut;
		hsigmaPSProtonCut = new TH2F("hsigmaPSProtonCut","hsigmaPSProtonCut",2700,167900,170600,10000,-10,10);
		hsigmaPSProtonCut->GetYaxis()->SetRangeUser(0.6,1.9);
		SetStyleHistoTH2ForGraphs(hsigmaPSProtonCut, "Run number","#sigma of n#sigma dE^{e^{-}}/dx/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		hsigmaPSProtonCut->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		sigmaNSigmadEdxPositronProtonCut_0010_Data->Draw("pesame");
		sigmaNSigmadEdxPositronProtonCut_1020_Data->Draw("pesame");
		sigmaNSigmadEdxPositronProtonCut_2040_Data->Draw("pesame");
		sigmaNSigmadEdxPositronProtonCut_4060_Data->Draw("pesame");
		sigmaNSigmadEdxPositronProtonCut_0010_MC->Draw("pesame");
		sigmaNSigmadEdxPositronProtonCut_1020_MC->Draw("pesame");
		sigmaNSigmadEdxPositronProtonCut_2040_MC->Draw("pesame");
		sigmaNSigmadEdxPositronProtonCut_4060_MC->Draw("pesame");

		TLegend* sigmaLPSProtonCut = new TLegend(0.22,0.93,0.99,.99);
		sigmaLPSProtonCut->SetFillColor(0);
		sigmaLPSProtonCut->SetLineColor(0);
		sigmaLPSProtonCut->SetTextSize(0.04);
		sigmaLPSProtonCut->SetNColumns(4);
		sigmaLPSProtonCut->AddEntry(sigmaNSigmadEdxPositronProtonCut_0010_Data,collisionSystem0010.Data(),"p");
		sigmaLPSProtonCut->AddEntry(sigmaNSigmadEdxPositronProtonCut_1020_Data,collisionSystem1020.Data(),"p");
		sigmaLPSProtonCut->AddEntry(sigmaNSigmadEdxPositronProtonCut_2040_Data,collisionSystem2040.Data(),"p");
		sigmaLPSProtonCut->AddEntry(sigmaNSigmadEdxPositronProtonCut_4060_Data,collisionSystem4060.Data(),"p");
		sigmaLPSProtonCut->AddEntry(sigmaNSigmadEdxPositronProtonCut_0010_MC,collisionSystem0010.Data(),"p");
		sigmaLPSProtonCut->AddEntry(sigmaNSigmadEdxPositronProtonCut_1020_MC,collisionSystem1020.Data(),"p");
		sigmaLPSProtonCut->AddEntry(sigmaNSigmadEdxPositronProtonCut_2040_MC,collisionSystem2040.Data(),"p");
		sigmaLPSProtonCut->AddEntry(sigmaNSigmadEdxPositronProtonCut_4060_MC,collisionSystem4060.Data(),"p");

		MC->Draw();
		sigmaLPSProtonCut->Draw();
		sigmaPSProtonCut->SaveAs(Form("%s/sigmaNSigmaPositronProtonCut.%s",outputDir.Data(),suffix.Data()));
		delete sigmaPSProtonCut;


		cout << "Plot width Nwidth dE/dx" << endl;
		TCanvas* widthELProtonCut = new TCanvas("widthELProtonCut","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(widthELProtonCut, 0.1, 0.03, 0.08, 0.08);
		TH2F *hwidthELProtonCut;
		hwidthELProtonCut = new TH2F("hwidthELProtonCut","hwidthELProtonCut",2700,167900,170600,10000,-10,10);
		hwidthELProtonCut->GetYaxis()->SetRangeUser(0.,5.);
		SetStyleHistoTH2ForGraphs(hwidthELProtonCut, "Run number","width of n#sigma dE^{e^{-}}/dx/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		hwidthELProtonCut->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		widthNSigmadEdxElectronProtonCut_0010_Data->Draw("pesame");
		widthNSigmadEdxElectronProtonCut_1020_Data->Draw("pesame");
		widthNSigmadEdxElectronProtonCut_2040_Data->Draw("pesame");
		widthNSigmadEdxElectronProtonCut_4060_Data->Draw("pesame");
		widthNSigmadEdxElectronProtonCut_0010_MC->Draw("pesame");
		widthNSigmadEdxElectronProtonCut_1020_MC->Draw("pesame");
		widthNSigmadEdxElectronProtonCut_2040_MC->Draw("pesame");
		widthNSigmadEdxElectronProtonCut_4060_MC->Draw("pesame");

		TLegend* widthLELProtonCut = new TLegend(0.22,0.93,0.99,.99);
		widthLELProtonCut->SetFillColor(0);
		widthLELProtonCut->SetLineColor(0);
		widthLELProtonCut->SetTextSize(0.04);
		widthLELProtonCut->SetNColumns(4);
		widthLELProtonCut->AddEntry(widthNSigmadEdxElectronProtonCut_0010_Data,collisionSystem0010.Data(),"p");
		widthLELProtonCut->AddEntry(widthNSigmadEdxElectronProtonCut_1020_Data,collisionSystem1020.Data(),"p");
		widthLELProtonCut->AddEntry(widthNSigmadEdxElectronProtonCut_2040_Data,collisionSystem2040.Data(),"p");
		widthLELProtonCut->AddEntry(widthNSigmadEdxElectronProtonCut_4060_Data,collisionSystem4060.Data(),"p");
		widthLELProtonCut->AddEntry(widthNSigmadEdxElectronProtonCut_0010_MC,collisionSystem0010.Data(),"p");
		widthLELProtonCut->AddEntry(widthNSigmadEdxElectronProtonCut_1020_MC,collisionSystem1020.Data(),"p");
		widthLELProtonCut->AddEntry(widthNSigmadEdxElectronProtonCut_2040_MC,collisionSystem2040.Data(),"p");
		widthLELProtonCut->AddEntry(widthNSigmadEdxElectronProtonCut_4060_MC,collisionSystem4060.Data(),"p");

		MC->Draw();
		widthLELProtonCut->Draw();
		widthELProtonCut->SaveAs(Form("%s/widthNSigmaElectronProtonCut.%s",outputDir.Data(),suffix.Data()));
		delete widthELProtonCut;

		cout << "Plot width Nwidth dE/dx" << endl;
		TCanvas* widthProtonCut = new TCanvas("widthPSProtonCut","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(widthProtonCut, 0.1, 0.03, 0.08, 0.08);
		TH2F *hwidthPSProtonCut;
		hwidthPSProtonCut = new TH2F("hwidthPSProtonCut","hwidthPSProtonCut",2700,167900,170600,10000,-10,10);
		hwidthPSProtonCut->GetYaxis()->SetRangeUser(0.,5.);
		SetStyleHistoTH2ForGraphs(hwidthPSProtonCut, "Run number","width of n#sigma dE^{e^{-}}/dx/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		hwidthPSProtonCut->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		widthNSigmadEdxPositronProtonCut_0010_Data->Draw("pesame");
		widthNSigmadEdxPositronProtonCut_1020_Data->Draw("pesame");
		widthNSigmadEdxPositronProtonCut_2040_Data->Draw("pesame");
		widthNSigmadEdxPositronProtonCut_4060_Data->Draw("pesame");
		widthNSigmadEdxPositronProtonCut_0010_MC->Draw("pesame");
		widthNSigmadEdxPositronProtonCut_1020_MC->Draw("pesame");
		widthNSigmadEdxPositronProtonCut_2040_MC->Draw("pesame");
		widthNSigmadEdxPositronProtonCut_4060_MC->Draw("pesame");

		TLegend* widthLPSProtonCut = new TLegend(0.22,0.93,0.99,.99);
		widthLPSProtonCut->SetFillColor(0);
		widthLPSProtonCut->SetLineColor(0);
		widthLPSProtonCut->SetTextSize(0.04);
		widthLPSProtonCut->SetNColumns(4);
		widthLPSProtonCut->AddEntry(widthNSigmadEdxPositronProtonCut_0010_Data,collisionSystem0010.Data(),"p");
		widthLPSProtonCut->AddEntry(widthNSigmadEdxPositronProtonCut_1020_Data,collisionSystem1020.Data(),"p");
		widthLPSProtonCut->AddEntry(widthNSigmadEdxPositronProtonCut_2040_Data,collisionSystem2040.Data(),"p");
		widthLPSProtonCut->AddEntry(widthNSigmadEdxPositronProtonCut_4060_Data,collisionSystem4060.Data(),"p");
		widthLPSProtonCut->AddEntry(widthNSigmadEdxPositronProtonCut_0010_MC,collisionSystem0010.Data(),"p");
		widthLPSProtonCut->AddEntry(widthNSigmadEdxPositronProtonCut_1020_MC,collisionSystem1020.Data(),"p");
		widthLPSProtonCut->AddEntry(widthNSigmadEdxPositronProtonCut_2040_MC,collisionSystem2040.Data(),"p");
		widthLPSProtonCut->AddEntry(widthNSigmadEdxPositronProtonCut_4060_MC,collisionSystem4060.Data(),"p");

		MC->Draw();
		widthLPSProtonCut->Draw();
		widthProtonCut->SaveAs(Form("%s/widthNSigmaPositronProtonCut.%s",outputDir.Data(),suffix.Data()));
		delete widthProtonCut;





		cout << "Plot Eta" << endl;
		TCanvas* canvaselectronEta = new TCanvas("canvaselectronEta","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvaselectronEta,0.1, 0.03, 0.08, 0.08);
		//canvaselectronEta->SetLogy();
		TH2F * histoelectronEta;
		histoelectronEta = new TH2F("histoelectronEta","histoelectronEta",2700,167900,170600,10000,-0.6,0.6);
		histoelectronEta->GetYaxis()->SetRangeUser(-0.04,0.08);
		SetStyleHistoTH2ForGraphs(histoelectronEta, "Run number","<#eta_{e^{-}}>/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,1.);
		histoelectronEta->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		ElectronEta_0010_Data->Draw("pesame");
		ElectronEta_1020_Data->Draw("pesame");
		ElectronEta_2040_Data->Draw("pesame");
		ElectronEta_4060_Data->Draw("pesame");
	// 	ElectronEta6080_Data->Draw("pesame");
		ElectronEta_0010_MC->Draw("pesame");
		ElectronEta_1020_MC->Draw("pesame");
		ElectronEta_2040_MC->Draw("pesame");
		ElectronEta_4060_MC->Draw("pesame");
	// 	ElectronEta6080_MC->Draw("pesame");
		TLegend* legendelectronEta = new TLegend(0.22,0.93,0.99,.99);
		legendelectronEta->SetFillColor(0);
		legendelectronEta->SetLineColor(0);
		legendelectronEta->SetTextSize(0.04);
		legendelectronEta->SetNColumns(4);
	// 	legendelectronEta->AddEntry((TObject*)0, "Data:","");
		legendelectronEta->AddEntry(ElectronEta_0010_Data,collisionSystem0010.Data(),"p");
		legendelectronEta->AddEntry(ElectronEta_1020_Data,collisionSystem1020.Data(),"p");
		legendelectronEta->AddEntry(ElectronEta_2040_Data,collisionSystem2040.Data(),"p");
		legendelectronEta->AddEntry(ElectronEta_4060_Data,collisionSystem4060.Data(),"p");
	//	legendelectronEta->AddEntry(ElectronEta6080_Data,collisionSystem6080.Data(),"p");
		legendelectronEta->AddEntry(ElectronEta_0010_MC,collisionSystem0010.Data(),"p");
		legendelectronEta->AddEntry(ElectronEta_1020_MC,collisionSystem1020.Data(),"p");
		legendelectronEta->AddEntry(ElectronEta_2040_MC,collisionSystem2040.Data(),"p");
		legendelectronEta->AddEntry(ElectronEta_4060_MC,collisionSystem4060.Data(),"p");
	//	legendelectronEta->AddEntry(ElectronEta6080_MC,collisionSystem6080.Data(),"p");
		MC->Draw();
		legendelectronEta->Draw();
		canvaselectronEta->SaveAs(Form("%s/ElectronEta.%s",outputDir.Data(),suffix.Data()));
		delete canvaselectronEta;

		TCanvas* canvaselectronEtaNeg = new TCanvas("canvaselectronEtaNeg","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvaselectronEtaNeg,0.1, 0.03, 0.08, 0.08);
		//canvaselectronEtaNeg->SetLogy();
		TH2F * histoelectronEtaNeg;
		histoelectronEtaNeg = new TH2F("histoelectronEtaNeg","histoelectronEtaNeg",2700,167900,170600,10000,-1.0,1.0);
		histoelectronEtaNeg->GetYaxis()->SetRangeUser(-0.53,-0.43);
		SetStyleHistoTH2ForGraphs(histoelectronEtaNeg, "Run number","<#eta_{e^{-}}>/N_{evt,run} (#eta < 0)",0.03,0.04, 0.03,0.04, 0.75,1.);
		histoelectronEtaNeg->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		ElectronEtaNeg_0010_Data->Draw("psame");
		ElectronEtaNeg_1020_Data->Draw("psame");
		ElectronEtaNeg_2040_Data->Draw("psame");
		ElectronEtaNeg_4060_Data->Draw("psame");
	// 	ElectronEtaNeg6080_Data->Draw("psame");
		ElectronEtaNeg_0010_MC->Draw("psame");
		ElectronEtaNeg_1020_MC->Draw("psame");
		ElectronEtaNeg_2040_MC->Draw("psame");
		ElectronEtaNeg_4060_MC->Draw("psame");
	// 	ElectronEtaNeg6080_MC->Draw("psame");
		TLegend* legendelectronEtaNeg = new TLegend(0.22,0.93,0.99,.99);
		legendelectronEtaNeg->SetFillColor(0);
		legendelectronEtaNeg->SetLineColor(0);
		legendelectronEtaNeg->SetTextSize(0.04);
		legendelectronEtaNeg->SetNColumns(4);
	// 	legendelectronEtaNeg->AddEntry((TObject*)0, "Data:","");
		legendelectronEtaNeg->AddEntry(ElectronEtaNeg_0010_Data,collisionSystem0010.Data(),"p");
		legendelectronEtaNeg->AddEntry(ElectronEtaNeg_1020_Data,collisionSystem1020.Data(),"p");
		legendelectronEtaNeg->AddEntry(ElectronEtaNeg_2040_Data,collisionSystem2040.Data(),"p");
		legendelectronEtaNeg->AddEntry(ElectronEtaNeg_4060_Data,collisionSystem4060.Data(),"p");
	//	legendelectronEtaNeg->AddEntry(ElectronEtaNeg6080_Data,collisionSystem6080.Data(),"p");
		legendelectronEtaNeg->AddEntry(ElectronEtaNeg_0010_MC,collisionSystem0010.Data(),"p");
		legendelectronEtaNeg->AddEntry(ElectronEtaNeg_1020_MC,collisionSystem1020.Data(),"p");
		legendelectronEtaNeg->AddEntry(ElectronEtaNeg_2040_MC,collisionSystem2040.Data(),"p");
		legendelectronEtaNeg->AddEntry(ElectronEtaNeg_4060_MC,collisionSystem4060.Data(),"p");
	//	legendelectronEtaNeg->AddEntry(ElectronEtaNeg6080_MC,collisionSystem6080.Data(),"p");
		MC->Draw();
		legendelectronEtaNeg->Draw();
		canvaselectronEtaNeg->SaveAs(Form("%s/ElectronEtaNeg.%s",outputDir.Data(),suffix.Data()));
		delete canvaselectronEtaNeg;

		TCanvas* canvaselectronEtaPos = new TCanvas("canvaselectronEtaPos","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvaselectronEtaPos,0.1, 0.03, 0.08, 0.08);
		//canvaselectronEtaPos->SetLogy();
		TH2F * histoelectronEtaPos;
		histoelectronEtaPos = new TH2F("histoelectronEtaPos","histoelectronEtaPos",2700,167900,170600,10000,-0.6,0.6);
		histoelectronEtaPos->GetYaxis()->SetRangeUser(0.43,0.53);
		SetStyleHistoTH2ForGraphs(histoelectronEtaPos, "Run number","<#eta_{e^{-}}>/N_{evt,run} (#eta > 0)",0.03,0.04, 0.03,0.04, 0.75,1.);
		histoelectronEtaPos->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		ElectronEtaPos_0010_Data->Draw("psame");
		ElectronEtaPos_1020_Data->Draw("psame");
		ElectronEtaPos_2040_Data->Draw("psame");
		ElectronEtaPos_4060_Data->Draw("psame");
	// 	ElectronEtaPos6080_Data->Draw("psame");
		ElectronEtaPos_0010_MC->Draw("psame");
		ElectronEtaPos_1020_MC->Draw("psame");
		ElectronEtaPos_2040_MC->Draw("psame");
		ElectronEtaPos_4060_MC->Draw("psame");
	// 	ElectronEtaPos6080_MC->Draw("psame");
		TLegend* legendelectronEtaPos = new TLegend(0.22,0.93,0.99,.99);
		legendelectronEtaPos->SetFillColor(0);
		legendelectronEtaPos->SetLineColor(0);
		legendelectronEtaPos->SetTextSize(0.04);
		legendelectronEtaPos->SetNColumns(4);
	// 	legendelectronEtaPos->AddEntry((TObject*)0, "Data:","");
		legendelectronEtaPos->AddEntry(ElectronEtaPos_0010_Data,collisionSystem0010.Data(),"p");
		legendelectronEtaPos->AddEntry(ElectronEtaPos_1020_Data,collisionSystem1020.Data(),"p");
		legendelectronEtaPos->AddEntry(ElectronEtaPos_2040_Data,collisionSystem2040.Data(),"p");
		legendelectronEtaPos->AddEntry(ElectronEtaPos_4060_Data,collisionSystem4060.Data(),"p");
	//	legendelectronEtaPos->AddEntry(ElectronEtaPos6080_Data,collisionSystem6080.Data(),"p");
		legendelectronEtaPos->AddEntry(ElectronEtaPos_0010_MC,collisionSystem0010.Data(),"p");
		legendelectronEtaPos->AddEntry(ElectronEtaPos_1020_MC,collisionSystem1020.Data(),"p");
		legendelectronEtaPos->AddEntry(ElectronEtaPos_2040_MC,collisionSystem2040.Data(),"p");
		legendelectronEtaPos->AddEntry(ElectronEtaPos_4060_MC,collisionSystem4060.Data(),"p");
	//	legendelectronEtaPos->AddEntry(ElectronEtaPos6080_MC,collisionSystem6080.Data(),"p");
		MC->Draw();
		legendelectronEtaPos->Draw();
		canvaselectronEtaPos->SaveAs(Form("%s/ElectronEtaPos.%s",outputDir.Data(),suffix.Data()));
		delete canvaselectronEtaPos;


		TCanvas* canvaspositronEta = new TCanvas("canvaspositronEta","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvaspositronEta, 0.1, 0.03, 0.08, 0.08);
		//canvaspositronEta->SetLogy();
		TH2F * histopositronEta;
		histopositronEta = new TH2F("histopositronEta","histopositronEta",2700,167900,170600,10000,-0.6,0.6);
		histopositronEta->GetYaxis()->SetRangeUser(-0.04,0.08);
		SetStyleHistoTH2ForGraphs(histopositronEta, "Run number","<#eta_{e^{+}}>/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,1.);
		histopositronEta->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		PositronEta_0010_Data->Draw("pesame");
		PositronEta_1020_Data->Draw("pesame");
		PositronEta_2040_Data->Draw("pesame");
		PositronEta_4060_Data->Draw("pesame");
	//	PositronEta6080_Data->Draw("pesame");
		PositronEta_0010_MC->Draw("pesame");
		PositronEta_1020_MC->Draw("pesame");
		PositronEta_2040_MC->Draw("pesame");
		PositronEta_4060_MC->Draw("pesame");
	//	PositronEta6080_MC->Draw("pesame");
		
		TLegend* legendpositronEta = new TLegend(0.22,0.93,0.99,.99);
		legendpositronEta->SetFillColor(0);
		legendpositronEta->SetLineColor(0);
		legendpositronEta->SetTextSize(0.04);
		legendpositronEta->SetNColumns(4);
	// 	legendpositronEta->AddEntry((TObject*)0, "Data:","");
		legendpositronEta->AddEntry(PositronEta_0010_Data,collisionSystem0010.Data(),"p");
		legendpositronEta->AddEntry(PositronEta_1020_Data,collisionSystem1020.Data(),"p");
		legendpositronEta->AddEntry(PositronEta_2040_Data,collisionSystem2040.Data(),"p");
		legendpositronEta->AddEntry(PositronEta_4060_Data,collisionSystem4060.Data(),"p");
	//	legendpositronEta->AddEntry(PositronEta6080_Data,collisionSystem6080.Data(),"p");
		legendpositronEta->AddEntry(PositronEta_0010_MC,collisionSystem0010.Data(),"p");
		legendpositronEta->AddEntry(PositronEta_1020_MC,collisionSystem1020.Data(),"p");
		legendpositronEta->AddEntry(PositronEta_2040_MC,collisionSystem2040.Data(),"p");
		legendpositronEta->AddEntry(PositronEta_4060_MC,collisionSystem4060.Data(),"p");
	//	legendpositronEta->AddEntry(PositronEta6080_MC,collisionSystem6080.Data(),"p");
		MC->Draw();
		legendpositronEta->Draw();
		canvaspositronEta->SaveAs(Form("%s/PositronEta.%s",outputDir.Data(),suffix.Data()));
		delete canvaspositronEta;

		TCanvas* canvaspositronEtaNeg = new TCanvas("canvaspositronEtaNeg","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvaspositronEtaNeg,0.1, 0.03, 0.08, 0.08);
		//canvaspositronEtaNeg->SetLogy();
		TH2F * histopositronEtaNeg;
		histopositronEtaNeg = new TH2F("histopositronEtaNeg","histopositronEtaNeg",2700,167900,170600,10000,-0.6,0.6);
		histopositronEtaNeg->GetYaxis()->SetRangeUser(-0.52,-0.44);
		SetStyleHistoTH2ForGraphs(histopositronEtaNeg, "Run number","<#eta_{e^{+}}>/N_{evt,run} (#eta < 0)",0.03,0.04, 0.03,0.04, 0.75,1.);
		histopositronEtaNeg->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		PositronEtaNeg_0010_Data->Draw("psame");
		PositronEtaNeg_1020_Data->Draw("psame");
		PositronEtaNeg_2040_Data->Draw("psame");
		PositronEtaNeg_4060_Data->Draw("psame");
	// 	PositronEtaNeg6080_Data->Draw("psame");
		PositronEtaNeg_0010_MC->Draw("psame");
		PositronEtaNeg_1020_MC->Draw("psame");
		PositronEtaNeg_2040_MC->Draw("psame");
		PositronEtaNeg_4060_MC->Draw("psame");
	// 	PositronEtaNeg6080_MC->Draw("psame");
		TLegend* legendpositronEtaNeg = new TLegend(0.22,0.93,0.99,.99);
		legendpositronEtaNeg->SetFillColor(0);
		legendpositronEtaNeg->SetLineColor(0);
		legendpositronEtaNeg->SetTextSize(0.04);
		legendpositronEtaNeg->SetNColumns(4);
	// 	legendpositronEtaNeg->AddEntry((TObject*)0, "Data:","");
		legendpositronEtaNeg->AddEntry(PositronEtaNeg_0010_Data,collisionSystem0010.Data(),"p");
		legendpositronEtaNeg->AddEntry(PositronEtaNeg_1020_Data,collisionSystem1020.Data(),"p");
		legendpositronEtaNeg->AddEntry(PositronEtaNeg_2040_Data,collisionSystem2040.Data(),"p");
		legendpositronEtaNeg->AddEntry(PositronEtaNeg_4060_Data,collisionSystem4060.Data(),"p");
	//	legendpositronEtaNeg->AddEntry(PositronEtaNeg6080_Data,collisionSystem6080.Data(),"p");
		legendpositronEtaNeg->AddEntry(PositronEtaNeg_0010_MC,collisionSystem0010.Data(),"p");
		legendpositronEtaNeg->AddEntry(PositronEtaNeg_1020_MC,collisionSystem1020.Data(),"p");
		legendpositronEtaNeg->AddEntry(PositronEtaNeg_2040_MC,collisionSystem2040.Data(),"p");
		legendpositronEtaNeg->AddEntry(PositronEtaNeg_4060_MC,collisionSystem4060.Data(),"p");
	//	legendpositronEtaNeg->AddEntry(PositronEtaNeg6080_MC,collisionSystem6080.Data(),"p");
		MC->Draw();
		legendpositronEtaNeg->Draw();
		canvaspositronEtaNeg->SaveAs(Form("%s/PositronEtaNeg.%s",outputDir.Data(),suffix.Data()));
		delete canvaspositronEtaNeg;

		TCanvas* canvaspositronEtaPos = new TCanvas("canvaspositronEtaPos","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvaspositronEtaPos,0.1, 0.03, 0.08, 0.08);
		//canvaspositronEtaPos->SetLogy();
		TH2F * histopositronEtaPos;
		histopositronEtaPos = new TH2F("histopositronEtaPos","histopositronEtaPos",2700,167900,170600,10000,-0.6,0.6);
		histopositronEtaPos->GetYaxis()->SetRangeUser(0.45,0.53);
		SetStyleHistoTH2ForGraphs(histopositronEtaPos, "Run number","<#eta_{e^{+}}>/N_{evt,run} (#eta > 0)",0.03,0.04, 0.03,0.04, 0.75,1.);
		histopositronEtaPos->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		PositronEtaPos_0010_Data->Draw("psame");
		PositronEtaPos_1020_Data->Draw("psame");
		PositronEtaPos_2040_Data->Draw("psame");
		PositronEtaPos_4060_Data->Draw("psame");
	// 	PositronEtaPos6080_Data->Draw("psame");
		PositronEtaPos_0010_MC->Draw("psame");
		PositronEtaPos_1020_MC->Draw("psame");
		PositronEtaPos_2040_MC->Draw("psame");
		PositronEtaPos_4060_MC->Draw("psame");
	// 	PositronEtaPos6080_MC->Draw("psame");
		TLegend* legendpositronEtaPos = new TLegend(0.22,0.93,0.99,.99);
		legendpositronEtaPos->SetFillColor(0);
		legendpositronEtaPos->SetLineColor(0);
		legendpositronEtaPos->SetTextSize(0.04);
		legendpositronEtaPos->SetNColumns(4);
	// 	legendpositronEtaPos->AddEntry((TObject*)0, "Data:","");
		legendpositronEtaPos->AddEntry(PositronEtaPos_0010_Data,collisionSystem0010.Data(),"p");
		legendpositronEtaPos->AddEntry(PositronEtaPos_1020_Data,collisionSystem1020.Data(),"p");
		legendpositronEtaPos->AddEntry(PositronEtaPos_2040_Data,collisionSystem2040.Data(),"p");
		legendpositronEtaPos->AddEntry(PositronEtaPos_4060_Data,collisionSystem4060.Data(),"p");
	//	legendpositronEtaPos->AddEntry(PositronEtaPos6080_Data,collisionSystem6080.Data(),"p");
		legendpositronEtaPos->AddEntry(PositronEtaPos_0010_MC,collisionSystem0010.Data(),"p");
		legendpositronEtaPos->AddEntry(PositronEtaPos_1020_MC,collisionSystem1020.Data(),"p");
		legendpositronEtaPos->AddEntry(PositronEtaPos_2040_MC,collisionSystem2040.Data(),"p");
		legendpositronEtaPos->AddEntry(PositronEtaPos_4060_MC,collisionSystem4060.Data(),"p");
	//	legendpositronEtaPos->AddEntry(PositronEtaPos6080_MC,collisionSystem6080.Data(),"p");
		MC->Draw();
		legendpositronEtaPos->Draw();
		canvaspositronEtaPos->SaveAs(Form("%s/PositronEtaPos.%s",outputDir.Data(),suffix.Data()));
		delete canvaspositronEtaPos;





		cout << "Plot photon" << endl;
		TCanvas* canvasphotonPt = new TCanvas("canvasphotonPt","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasphotonPt, 0.1, 0.03, 0.08, 0.08);
		//canvasphotonPt->SetLogy();
		TH2F * histophotonPt;
		histophotonPt = new TH2F("histophotonPt","histophotonPt",2700,167900,170600,10000,-1,1);
		histophotonPt->GetYaxis()->SetRangeUser(0.5,0.7);
		SetStyleHistoTH2ForGraphs(histophotonPt, "Run number","<#eta_{#gamma}>/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,1.);
		histophotonPt->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		PhotonPt_0010_Data->Draw("pesame");
		PhotonPt_1020_Data->Draw("pesame");
		PhotonPt_2040_Data->Draw("pesame");
		PhotonPt_4060_Data->Draw("pesame");
	// 	PhotonPt6080_Data->Draw("pesame");
		PhotonPt_0010_MC->Draw("pesame");
		PhotonPt_1020_MC->Draw("pesame");
		PhotonPt_2040_MC->Draw("pesame");
		PhotonPt_4060_MC->Draw("pesame");
	// 	PhotonPt6080_MC->Draw("pesame");
		TLegend* legendphotonPt = new TLegend(0.22,0.93,0.99,.99);
		legendphotonPt->SetFillColor(0);
		legendphotonPt->SetLineColor(0);
		legendphotonPt->SetTextSize(0.04);
		legendphotonPt->SetNColumns(4);
	// 	legendphotonPt->AddEntry((TObject*)0, "Data:","");
		legendphotonPt->AddEntry(PhotonPt_0010_Data,collisionSystem0010.Data(),"p");
		legendphotonPt->AddEntry(PhotonPt_1020_Data,collisionSystem1020.Data(),"p");
		legendphotonPt->AddEntry(PhotonPt_2040_Data,collisionSystem2040.Data(),"p");
		legendphotonPt->AddEntry(PhotonPt_4060_Data,collisionSystem4060.Data(),"p");
	//	legendphotonPt->AddEntry(PhotonPt6080_Data,collisionSystem6080.Data(),"p");
		legendphotonPt->AddEntry(PhotonPt_0010_MC,collisionSystem0010.Data(),"p");
		legendphotonPt->AddEntry(PhotonPt_1020_MC,collisionSystem1020.Data(),"p");
		legendphotonPt->AddEntry(PhotonPt_2040_MC,collisionSystem2040.Data(),"p");
		legendphotonPt->AddEntry(PhotonPt_4060_MC,collisionSystem4060.Data(),"p");
	//	legendphotonPt->AddEntry(PhotonPt6080_MC,collisionSystem6080.Data(),"p");
		MC->Draw();
		legendphotonPt->Draw();
		canvasphotonPt->SaveAs(Form("%s/PhotonPt.%s",outputDir.Data(),suffix.Data()));
		delete canvasphotonPt;



		TCanvas* canvasphotonEta = new TCanvas("canvasphotonEta","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasphotonEta, 0.1, 0.03, 0.08, 0.08);
		//canvasphotonEta->SetLogy();
		TH2F * histophotonEta;
		histophotonEta = new TH2F("histophotonEta","histophotonEta",2700,167900,170600,10000,-1,1);
		histophotonEta->GetYaxis()->SetRangeUser(-0.04,0.08);
		SetStyleHistoTH2ForGraphs(histophotonEta, "Run number","<#eta_{#gamma}>/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,1.);
		histophotonEta->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		PhotonEta_0010_Data->Draw("pesame");
		PhotonEta_1020_Data->Draw("pesame");
		PhotonEta_2040_Data->Draw("pesame");
		PhotonEta_4060_Data->Draw("pesame");
	// 	PhotonEta6080_Data->Draw("pesame");
		PhotonEta_0010_MC->Draw("pesame");
		PhotonEta_1020_MC->Draw("pesame");
		PhotonEta_2040_MC->Draw("pesame");
		PhotonEta_4060_MC->Draw("pesame");
	// 	PhotonEta6080_MC->Draw("pesame");
		TLegend* legendphotonEta = new TLegend(0.22,0.93,0.99,.99);
		legendphotonEta->SetFillColor(0);
		legendphotonEta->SetLineColor(0);
		legendphotonEta->SetTextSize(0.04);
		legendphotonEta->SetNColumns(4);
	// 	legendphotonEta->AddEntry((TObject*)0, "Data:","");
		legendphotonEta->AddEntry(PhotonEta_0010_Data,collisionSystem0010.Data(),"p");
		legendphotonEta->AddEntry(PhotonEta_1020_Data,collisionSystem1020.Data(),"p");
		legendphotonEta->AddEntry(PhotonEta_2040_Data,collisionSystem2040.Data(),"p");
		legendphotonEta->AddEntry(PhotonEta_4060_Data,collisionSystem4060.Data(),"p");
	//	legendphotonEta->AddEntry(PhotonEta6080_Data,collisionSystem6080.Data(),"p");
		legendphotonEta->AddEntry(PhotonEta_0010_MC,collisionSystem0010.Data(),"p");
		legendphotonEta->AddEntry(PhotonEta_1020_MC,collisionSystem1020.Data(),"p");
		legendphotonEta->AddEntry(PhotonEta_2040_MC,collisionSystem2040.Data(),"p");
		legendphotonEta->AddEntry(PhotonEta_4060_MC,collisionSystem4060.Data(),"p");
	//	legendphotonEta->AddEntry(PhotonEta6080_MC,collisionSystem6080.Data(),"p");
		MC->Draw();
		legendphotonEta->Draw();
		canvasphotonEta->SaveAs(Form("%s/PhotonEta.%s",outputDir.Data(),suffix.Data()));
		delete canvasphotonEta;

		TCanvas* canvasphotonEtaNeg = new TCanvas("canvasphotonEtaNeg","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasphotonEtaNeg, 0.1, 0.03, 0.08, 0.08);
		//canvasphotonEtaNeg->SetLogy();
		TH2F * histophotonEtaNeg;
		histophotonEtaNeg = new TH2F("histophotonEtaNeg","histophotonEtaNeg",2700,167900,170600,10000,-1,1);
		histophotonEtaNeg->GetYaxis()->SetRangeUser(-0.52,-0.44);
		SetStyleHistoTH2ForGraphs(histophotonEtaNeg, "Run number","<#eta_{#gamma}>/N_{evt,run} (#eta <0)",0.03,0.04, 0.03,0.04, 0.75,1.);
		histophotonEtaNeg->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		PhotonEtaNeg_0010_Data->Draw("psame");
		PhotonEtaNeg_1020_Data->Draw("psame");
		PhotonEtaNeg_2040_Data->Draw("psame");
		PhotonEtaNeg_4060_Data->Draw("psame");
	// 	PhotonEtaNeg6080_Data->Draw("psame");
		PhotonEtaNeg_0010_MC->Draw("psame");
		PhotonEtaNeg_1020_MC->Draw("psame");
		PhotonEtaNeg_2040_MC->Draw("psame");
		PhotonEtaNeg_4060_MC->Draw("psame");
	// 	PhotonEtaNeg6080_MC->Draw("psame");
		TLegend* legendphotonEtaNeg = new TLegend(0.22,0.93,0.99,.99);
		legendphotonEtaNeg->SetFillColor(0);
		legendphotonEtaNeg->SetLineColor(0);
		legendphotonEtaNeg->SetTextSize(0.04);
		legendphotonEtaNeg->SetNColumns(4);
	// 	legendphotonEtaNeg->AddEntry((TObject*)0, "Data:","");
		legendphotonEtaNeg->AddEntry(PhotonEtaNeg_0010_Data,collisionSystem0010.Data(),"p");
		legendphotonEtaNeg->AddEntry(PhotonEtaNeg_1020_Data,collisionSystem1020.Data(),"p");
		legendphotonEtaNeg->AddEntry(PhotonEtaNeg_2040_Data,collisionSystem2040.Data(),"p");
		legendphotonEtaNeg->AddEntry(PhotonEtaNeg_4060_Data,collisionSystem4060.Data(),"p");
	//	legendphotonEtaNeg->AddEntry(PhotonEtaNeg6080_Data,collisionSystem6080.Data(),"p");
		legendphotonEtaNeg->AddEntry(PhotonEtaNeg_0010_MC,collisionSystem0010.Data(),"p");
		legendphotonEtaNeg->AddEntry(PhotonEtaNeg_1020_MC,collisionSystem1020.Data(),"p");
		legendphotonEtaNeg->AddEntry(PhotonEtaNeg_2040_MC,collisionSystem2040.Data(),"p");
		legendphotonEtaNeg->AddEntry(PhotonEtaNeg_4060_MC,collisionSystem4060.Data(),"p");
	//	legendphotonEtaNeg->AddEntry(PhotonEtaNeg6080_MC,collisionSystem6080.Data(),"p");
		MC->Draw();
		legendphotonEtaNeg->Draw();
		canvasphotonEtaNeg->SaveAs(Form("%s/PhotonEtaNeg.%s",outputDir.Data(),suffix.Data()));
		delete canvasphotonEtaNeg;

		TCanvas* canvasphotonEtaPos = new TCanvas("canvasphotonEtaPos","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasphotonEtaPos, 0.1, 0.03, 0.08, 0.08);
		//canvasphotonEtaPos->SetLogy();
		TH2F * histophotonEtaPos;
		histophotonEtaPos = new TH2F("histophotonEtaPos","histophotonEtaPos",2700,167900,170600,10000,-1,1);
		histophotonEtaPos->GetYaxis()->SetRangeUser(0.45,0.53);
		SetStyleHistoTH2ForGraphs(histophotonEtaPos, "Run number","<#eta_{#gamma}>/N_{evt,run} (#eta <0)",0.03,0.04, 0.03,0.04, 0.75,1.);
		histophotonEtaPos->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		PhotonEtaPos_0010_Data->Draw("psame");
		PhotonEtaPos_1020_Data->Draw("psame");
		PhotonEtaPos_2040_Data->Draw("psame");
		PhotonEtaPos_4060_Data->Draw("psame");
	// 	PhotonEtaPos6080_Data->Draw("psame");
		PhotonEtaPos_0010_MC->Draw("psame");
		PhotonEtaPos_1020_MC->Draw("psame");
		PhotonEtaPos_2040_MC->Draw("psame");
		PhotonEtaPos_4060_MC->Draw("psame");
	// 	PhotonEtaPos6080_MC->Draw("psame");
		TLegend* legendphotonEtaPos = new TLegend(0.22,0.93,0.99,.99);
		legendphotonEtaPos->SetFillColor(0);
		legendphotonEtaPos->SetLineColor(0);
		legendphotonEtaPos->SetTextSize(0.04);
		legendphotonEtaPos->SetNColumns(4);
	// 	legendphotonEtaPos->AddEntry((TObject*)0, "Data:","");
		legendphotonEtaPos->AddEntry(PhotonEtaPos_0010_Data,collisionSystem0010.Data(),"p");
		legendphotonEtaPos->AddEntry(PhotonEtaPos_1020_Data,collisionSystem1020.Data(),"p");
		legendphotonEtaPos->AddEntry(PhotonEtaPos_2040_Data,collisionSystem2040.Data(),"p");
		legendphotonEtaPos->AddEntry(PhotonEtaPos_4060_Data,collisionSystem4060.Data(),"p");
	//	legendphotonEtaPos->AddEntry(PhotonEtaPos6080_Data,collisionSystem6080.Data(),"p");
		legendphotonEtaPos->AddEntry(PhotonEtaPos_0010_MC,collisionSystem0010.Data(),"p");
		legendphotonEtaPos->AddEntry(PhotonEtaPos_1020_MC,collisionSystem1020.Data(),"p");
		legendphotonEtaPos->AddEntry(PhotonEtaPos_2040_MC,collisionSystem2040.Data(),"p");
		legendphotonEtaPos->AddEntry(PhotonEtaPos_4060_MC,collisionSystem4060.Data(),"p");
	//	legendphotonEtaPos->AddEntry(PhotonEtaPos6080_MC,collisionSystem6080.Data(),"p");
		MC->Draw();
		legendphotonEtaPos->Draw();
		canvasphotonEtaPos->SaveAs(Form("%s/PhotonEtaPos.%s",outputDir.Data(),suffix.Data()));
		delete canvasphotonEtaPos;
	} 

	if(!partialOutput){
		cout << "Plot gamma sectors" << endl;
		///////////////////////////////////////////////////////// A side ///////////////////////////////////////////////////////////////////////
		TCanvas* canvasGammaAllSectorsAside = new TCanvas("canvasGammaAllSectorsAside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaAllSectorsAside,0.1, 0.03, 0.08, 0.08);
		TH2F * histoGAllSectorsAside;
		histoGAllSectorsAside = new TH2F("histoGAllSectorsAside","histoGAllSectorsAside",2700,167900,170600,10000,-100,100);
		histoGAllSectorsAside->GetYaxis()->SetRangeUser(0.,22.);
		SetStyleHistoTH2ForGraphs(histoGAllSectorsAside, "Run number","N_{#gamma}^{A} all sectors/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGAllSectorsAside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		TPaveText *sectorAll = new TPaveText(0.15,0.82,0.3,0.88,"NDC");
		sectorAll->SetTextSize(0.04);
		sectorAll->SetLineColor(0);
		sectorAll->SetFillColor(0);
		sectorAll->SetFillStyle(0);
		sectorAll->SetShadowColor(0);
		sectorAll->SetTextAlign(12);
		sectorAll->AddText("All sectors");
		sectorAll->Draw(); 

		PhotonAllSectorsEtaPos_0010_Data->Draw("psame");
		PhotonAllSectorsEtaPos_1020_Data->Draw("psame");
		PhotonAllSectorsEtaPos_2040_Data->Draw("psame");
		PhotonAllSectorsEtaPos_4060_Data->Draw("psame");
	//	PhotonAllSectorsEtaPos6080_Data->Draw("psame");
		PhotonAllSectorsEtaPos_0010_MC->Draw("psame");
		PhotonAllSectorsEtaPos_1020_MC->Draw("psame");
		PhotonAllSectorsEtaPos_2040_MC->Draw("psame");
		PhotonAllSectorsEtaPos_4060_MC->Draw("psame");
	//	PhotonAllSectorsEtaPos6080_MC->Draw("psame");
		
		TLegend* legendGAllSectors = new TLegend(0.22,0.93,0.99,.99);
		legendGAllSectors->SetFillColor(0);
		legendGAllSectors->SetLineColor(0);
		legendGAllSectors->SetTextSize(0.04);
		legendGAllSectors->SetNColumns(4);
		legendGAllSectors->AddEntry(PhotonAllSectorsEtaPos_0010_Data,collisionSystem0010.Data(),"p");
		legendGAllSectors->AddEntry(PhotonAllSectorsEtaPos_1020_Data,collisionSystem1020.Data(),"p");
		legendGAllSectors->AddEntry(PhotonAllSectorsEtaPos_2040_Data,collisionSystem2040.Data(),"p");
		legendGAllSectors->AddEntry(PhotonAllSectorsEtaPos_4060_Data,collisionSystem4060.Data(),"p");
	//	legendGAllSectors->AddEntry(PhotonAllSectorsEtaPos6080_Data,collisionSystem6080.Data(),"p");
		legendGAllSectors->AddEntry(PhotonAllSectorsEtaPos_0010_MC,collisionSystem0010.Data(),"p");
		legendGAllSectors->AddEntry(PhotonAllSectorsEtaPos_1020_MC,collisionSystem1020.Data(),"p");
		legendGAllSectors->AddEntry(PhotonAllSectorsEtaPos_2040_MC,collisionSystem2040.Data(),"p");
		legendGAllSectors->AddEntry(PhotonAllSectorsEtaPos_4060_MC,collisionSystem4060.Data(),"p");
	//	legendGAllSectors->AddEntry(PhotonAllSectorsEtaPos6080_MC,collisionSystem6080.Data(),"p");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaAllSectorsAside->SaveAs(Form("%s/GammaAllSectors_Aside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaAllSectorsAside;


		
		//Sector0
		TCanvas* canvasGammaSector0Aside = new TCanvas("canvasGammaSector0Aside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector0Aside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector0Aside;
		histoGSector0Aside = new TH2F("histoGSector0Aside","histoGSector0Aside",2700,167900,170600,10000,-100,100);
		histoGSector0Aside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector0Aside, "Run number","N_{#gamma}^{A} sector 0/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector0Aside->Draw("copy");
			DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		TPaveText *sector0 = new TPaveText(0.15,0.82,0.3,0.88,"NDC");
		sector0->SetTextSize(0.04);
		sector0->SetLineColor(0);
		sector0->SetFillStyle(0);
		sector0->SetFillColor(0);
		sector0->SetShadowColor(0);
		sector0->SetTextAlign(12);
		sector0->AddText("Sector 0");
		sector0->Draw();

		PhotonSector0EtaPos_0010_Data->Draw("psame");
		PhotonSector0EtaPos_1020_Data->Draw("psame");
		PhotonSector0EtaPos_2040_Data->Draw("psame");
		PhotonSector0EtaPos_4060_Data->Draw("psame");
	//	PhotonSector0EtaPos6080_Data->Draw("psame");
		PhotonSector0EtaPos_0010_MC->Draw("psame");
		PhotonSector0EtaPos_1020_MC->Draw("psame");
		PhotonSector0EtaPos_2040_MC->Draw("psame");
		PhotonSector0EtaPos_4060_MC->Draw("psame");
	//	PhotonSector0EtaPos6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector0Aside->SaveAs(Form("%s/GammaSector0_Aside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector0Aside;


		//Sector1
		TCanvas* canvasGammaSector1Aside = new TCanvas("canvasGammaSector1Aside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector1Aside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector1Aside;
		histoGSector1Aside = new TH2F("histoGSector1Aside","histoGSector1Aside",2700,167900,170600,10000,-100,100);
		histoGSector1Aside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector1Aside, "Run number","N_{#gamma}^{A} sector 1/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector1Aside->Draw("copy");
			DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		TPaveText *sector1 = new TPaveText(0.15,0.82,0.3,0.88,"NDC");
		sector1->SetTextSize(0.04);
		sector1->SetLineColor(0);
		sector1->SetFillStyle(0);
		sector1->SetFillColor(0);
		sector1->SetShadowColor(0);
		sector1->SetTextAlign(12);
		sector1->AddText("Sector 1");
		sector1->Draw();

		PhotonSector1EtaPos_0010_Data->Draw("psame");
		PhotonSector1EtaPos_1020_Data->Draw("psame");
		PhotonSector1EtaPos_2040_Data->Draw("psame");
		PhotonSector1EtaPos_4060_Data->Draw("psame");
	//	PhotonSector1EtaPos6080_Data->Draw("psame");
		PhotonSector1EtaPos_0010_MC->Draw("psame");
		PhotonSector1EtaPos_1020_MC->Draw("psame");
		PhotonSector1EtaPos_2040_MC->Draw("psame");
		PhotonSector1EtaPos_4060_MC->Draw("psame");
	//	PhotonSector1EtaPos6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector1Aside->SaveAs(Form("%s/GammaSector1_Aside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector1Aside;


		//Sector2
		TCanvas* canvasGammaSector2Aside = new TCanvas("canvasGammaSector2Aside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector2Aside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector2Aside;
		histoGSector2Aside = new TH2F("histoGSector2Aside","histoGSector2Aside",2700,167900,170600,10000,-100,100);
		histoGSector2Aside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector2Aside, "Run number","N_{#gamma}^{A} sector 2/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector2Aside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		TPaveText *sector2 = new TPaveText(0.15,0.82,0.3,0.88,"NDC");
		sector2->SetTextSize(0.04);
		sector2->SetLineColor(0);
		sector2->SetFillStyle(0);
		sector2->SetFillColor(0);
		sector2->SetShadowColor(0);
		sector2->SetTextAlign(12);
		sector2->AddText("Sector 2");
		sector2->Draw(); 

		PhotonSector2EtaPos_0010_Data->Draw("psame");
		PhotonSector2EtaPos_1020_Data->Draw("psame");
		PhotonSector2EtaPos_2040_Data->Draw("psame");
	//	PhotonSector2EtaPos_4060_Data->Draw("psame");
	//	PhotonSector2EtaPos6080_Data->Draw("psame");
		PhotonSector2EtaPos_0010_MC->Draw("psame");
		PhotonSector2EtaPos_1020_MC->Draw("psame");
		PhotonSector2EtaPos_2040_MC->Draw("psame");
	// 	PhotonSector2EtaPos_4060_MC->Draw("psame");
	//	PhotonSector2EtaPos6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector2Aside->SaveAs(Form("%s/GammaSector2_Aside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector2Aside;


		//Sector3
		TCanvas* canvasGammaSector3Aside = new TCanvas("canvasGammaSector3Aside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector3Aside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector3Aside;
		histoGSector3Aside = new TH2F("histoGSector3Aside","histoGSector3",2700,167900,170600,10000,-100,100);
		histoGSector3Aside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector3Aside, "Run number","N_{#gamma}^{A} sector 3/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector3Aside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		TPaveText *sector3 = new TPaveText(0.15,0.82,0.3,0.88,"NDC");
		sector3->SetTextSize(0.04);
		sector3->SetLineColor(0);
		sector3->SetFillStyle(0);
		sector3->SetFillColor(0);
		sector3->SetShadowColor(0);
		sector3->SetTextAlign(12);
		sector3->AddText("Sector 3");
		sector3->Draw(); 

		PhotonSector3EtaPos_0010_Data->Draw("psame");
		PhotonSector3EtaPos_1020_Data->Draw("psame");
		PhotonSector3EtaPos_2040_Data->Draw("psame");
		PhotonSector3EtaPos_4060_Data->Draw("psame");
	//	PhotonSector3EtaPos6080_Data->Draw("psame");
		PhotonSector3EtaPos_0010_MC->Draw("psame");
		PhotonSector3EtaPos_1020_MC->Draw("psame");
		PhotonSector3EtaPos_2040_MC->Draw("psame");
		PhotonSector3EtaPos_4060_MC->Draw("psame");
	//	PhotonSector3EtaPos6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector3Aside->SaveAs(Form("%s/GammaSector3_Aside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector3Aside;


		//Sector4
		TCanvas* canvasGammaSector4Aside = new TCanvas("canvasGammaSector4Aside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector4Aside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector4Aside;
		histoGSector4Aside = new TH2F("histoGSector4Aside","histoGSector4Aside",2700,167900,170600,10000,-100,100);
		histoGSector4Aside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector4Aside, "Run number","N_{#gamma}^{A} sector 4/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector4Aside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		TPaveText *sector4 = new TPaveText(0.15,0.82,0.3,0.88,"NDC");
		sector4->SetTextSize(0.04);
		sector4->SetLineColor(0);
		sector4->SetFillStyle(0);
		sector4->SetFillColor(0);
		sector4->SetShadowColor(0);
		sector4->SetTextAlign(12);
		sector4->AddText("Sector 4");
		sector4->Draw(); 

		PhotonSector4EtaPos_0010_Data->Draw("psame");
		PhotonSector4EtaPos_1020_Data->Draw("psame");
		PhotonSector4EtaPos_2040_Data->Draw("psame");
		PhotonSector4EtaPos_4060_Data->Draw("psame");
	//	PhotonSector4EtaPos6080_Data->Draw("psame");
		PhotonSector4EtaPos_0010_MC->Draw("psame");
		PhotonSector4EtaPos_1020_MC->Draw("psame");
		PhotonSector4EtaPos_2040_MC->Draw("psame");
		PhotonSector4EtaPos_4060_MC->Draw("psame");
	//	PhotonSector4EtaPos6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector4Aside->SaveAs(Form("%s/GammaSector4_Aside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector4Aside;


		//Sector5
		TCanvas* canvasGammaSector5Aside = new TCanvas("canvasGammaSector5Aside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector5Aside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector5Aside;
		histoGSector5Aside = new TH2F("histoGSector5Aside","histoGSector5Aside",2700,167900,170600,10000,-100,100);
		histoGSector5Aside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector5Aside, "Run number","N_{#gamma}^{A} sector 5/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector5Aside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		TPaveText *sector5 = new TPaveText(0.15,0.82,0.3,0.88,"NDC");
		sector5->SetTextSize(0.04);
		sector5->SetLineColor(0);
		sector5->SetFillStyle(0);
		sector5->SetFillColor(0);
		sector5->SetShadowColor(0);
		sector5->SetTextAlign(12);
		sector5->AddText("Sector 5");
		sector5->Draw(); 

		PhotonSector5EtaPos_0010_Data->Draw("psame");
		PhotonSector5EtaPos_1020_Data->Draw("psame");
		PhotonSector5EtaPos_2040_Data->Draw("psame");
		PhotonSector5EtaPos_4060_Data->Draw("psame");
	//	PhotonSector5EtaPos6080_Data->Draw("psame");
		PhotonSector5EtaPos_0010_MC->Draw("psame");
		PhotonSector5EtaPos_1020_MC->Draw("psame");
		PhotonSector5EtaPos_2040_MC->Draw("psame");
		PhotonSector5EtaPos_4060_MC->Draw("psame");
	//	PhotonSector5EtaPos6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector5Aside->SaveAs(Form("%s/GammaSector5_Aside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector5Aside;


		//Sector6
		TCanvas* canvasGammaSector6Aside = new TCanvas("canvasGammaSector6Aside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector6Aside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector6Aside;
		histoGSector6Aside = new TH2F("histoGSector6Aside","histoGSector6Aside",2700,167900,170600,10000,-100,100);
		histoGSector6Aside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector6Aside, "Run number","N_{#gamma}^{A} sector 6/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector6Aside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		TPaveText *sector6 = new TPaveText(0.15,0.82,0.3,0.88,"NDC");
		sector6->SetTextSize(0.04);
		sector6->SetLineColor(0);
		sector6->SetFillStyle(0);
		sector6->SetFillColor(0);
		sector6->SetShadowColor(0);
		sector6->SetTextAlign(12);
		sector6->AddText("Sector 6");
		sector6->Draw();

		PhotonSector6EtaPos_0010_Data->Draw("psame");
		PhotonSector6EtaPos_1020_Data->Draw("psame");
		PhotonSector6EtaPos_2040_Data->Draw("psame");
		PhotonSector6EtaPos_4060_Data->Draw("psame");
	//	PhotonSector6EtaPos6080_Data->Draw("psame");
		PhotonSector6EtaPos_0010_MC->Draw("psame");
		PhotonSector6EtaPos_1020_MC->Draw("psame");
		PhotonSector6EtaPos_2040_MC->Draw("psame");
		PhotonSector6EtaPos_4060_MC->Draw("psame");
	//	PhotonSector6EtaPos6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector6Aside->SaveAs(Form("%s/GammaSector6_Aside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector6Aside;


		//Sector7
		TCanvas* canvasGammaSector7Aside = new TCanvas("canvasGammaSector7Aside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector7Aside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector7Aside;
		histoGSector7Aside = new TH2F("histoGSector7Aside","histoGSector7Aside",2700,167900,170600,10000,-100,100);
		histoGSector7Aside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector7Aside, "Run number","N_{#gamma}^{A} sector 7/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector7Aside->Draw("copy");
			DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		TPaveText *sector7 = new TPaveText(0.15,0.82,0.3,0.88,"NDC");
		sector7->SetTextSize(0.04);
		sector7->SetLineColor(0);
		sector7->SetFillStyle(0);
		sector7->SetFillColor(0);
		sector7->SetShadowColor(0);
		sector7->SetTextAlign(12);
		sector7->AddText("Sector 7");
		sector7->Draw();

		PhotonSector7EtaPos_0010_Data->Draw("psame");
		PhotonSector7EtaPos_1020_Data->Draw("psame");
		PhotonSector7EtaPos_2040_Data->Draw("psame");
		PhotonSector7EtaPos_4060_Data->Draw("psame");
	//	PhotonSector7EtaPos6080_Data->Draw("psame");
		PhotonSector7EtaPos_0010_MC->Draw("psame");
		PhotonSector7EtaPos_1020_MC->Draw("psame");
		PhotonSector7EtaPos_2040_MC->Draw("psame");
		PhotonSector7EtaPos_4060_MC->Draw("psame");
	//	PhotonSector7EtaPos6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector7Aside->SaveAs(Form("%s/GammaSector7_Aside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector7Aside;


		//Sector8
		TCanvas* canvasGammaSector8Aside = new TCanvas("canvasGammaSector8Aside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector8Aside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector8Aside;
		histoGSector8Aside = new TH2F("histoGSector8Aside","histoGSector8Aside",2700,167900,170600,10000,-100,100);
		histoGSector8Aside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector8Aside, "Run number","N_{#gamma}^{A} sector 8/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector8Aside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		TPaveText *sector8 = new TPaveText(0.15,0.82,0.3,0.88,"NDC");
		sector8->SetTextSize(0.04);
		sector8->SetLineColor(0);
		sector8->SetFillStyle(0);
		sector8->SetFillColor(0);
		sector8->SetShadowColor(0);
		sector8->SetTextAlign(12);
		sector8->AddText("Sector 8");
		sector8->Draw();

		PhotonSector8EtaPos_0010_Data->Draw("psame");
		PhotonSector8EtaPos_1020_Data->Draw("psame");
		PhotonSector8EtaPos_2040_Data->Draw("psame");
		PhotonSector8EtaPos_4060_Data->Draw("psame");
	//	PhotonSector8EtaPos6080_Data->Draw("psame");
		PhotonSector8EtaPos_0010_MC->Draw("psame");
		PhotonSector8EtaPos_1020_MC->Draw("psame");
		PhotonSector8EtaPos_2040_MC->Draw("psame");
		PhotonSector8EtaPos_4060_MC->Draw("psame");
	//	PhotonSector8EtaPos6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();

		canvasGammaSector8Aside->SaveAs(Form("%s/GammaSector8_Aside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector8Aside;

		//Sector9
		TCanvas* canvasGammaSector9Aside = new TCanvas("canvasGammaSector9Aside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector9Aside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector9Aside;
		histoGSector9Aside = new TH2F("histoGSector9Aside","histoGSector9Aside",2700,167900,170600,10000,-100,100);
		histoGSector9Aside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector9Aside, "Run number","N_{#gamma}^{A} sector 9/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector9Aside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		TPaveText *sector9 = new TPaveText(0.15,0.82,0.3,0.88,"NDC");
		sector9->SetTextSize(0.04);
		sector9->SetLineColor(0);
		sector9->SetFillStyle(0);
		sector9->SetFillColor(0);
		sector9->SetShadowColor(0);
		sector9->SetTextAlign(12);
		sector9->AddText("Sector 9");
		sector9->Draw(); 

		PhotonSector9EtaPos_0010_Data->Draw("psame");
		PhotonSector9EtaPos_1020_Data->Draw("psame");
		PhotonSector9EtaPos_2040_Data->Draw("psame");
		PhotonSector9EtaPos_4060_Data->Draw("psame");
	//	PhotonSector9EtaPos6080_Data->Draw("psame");
		PhotonSector9EtaPos_0010_MC->Draw("psame");
		PhotonSector9EtaPos_1020_MC->Draw("psame");
		PhotonSector9EtaPos_2040_MC->Draw("psame");
		PhotonSector9EtaPos_4060_MC->Draw("psame");
	//	PhotonSector9EtaPos6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector9Aside->SaveAs(Form("%s/GammaSector9_Aside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector9Aside;


		//Sector10
		TCanvas* canvasGammaSector10Aside = new TCanvas("canvasGammaSector10Aside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector10Aside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector10Aside;
		histoGSector10Aside = new TH2F("histoGSector10Aside","histoGSector10Aside",2700,167900,170600,10000,-100,100);
		histoGSector10Aside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector10Aside, "Run number","N_{#gamma}^{A} sector 10/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector10Aside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		TPaveText *sector10 = new TPaveText(0.15,0.82,0.3,0.88,"NDC");
		sector10->SetTextSize(0.04);
		sector10->SetLineColor(0);
		sector10->SetFillStyle(0);
		sector10->SetFillColor(0);
		sector10->SetShadowColor(0);
		sector10->SetTextAlign(12);
		sector10->AddText("Sector 10");
		sector10->Draw();

		PhotonSector10EtaPos_0010_Data->Draw("psame");
		PhotonSector10EtaPos_1020_Data->Draw("psame");
		PhotonSector10EtaPos_2040_Data->Draw("psame");
		PhotonSector10EtaPos_4060_Data->Draw("psame");
	//	PhotonSector10EtaPos6080_Data->Draw("psame");
		PhotonSector10EtaPos_0010_MC->Draw("psame");
		PhotonSector10EtaPos_1020_MC->Draw("psame");
		PhotonSector10EtaPos_2040_MC->Draw("psame");
		PhotonSector10EtaPos_4060_MC->Draw("psame");
	//	PhotonSector10EtaPos6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector10Aside->SaveAs(Form("%s/GammaSector10_Aside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector10Aside;


		//Sector11
		TCanvas* canvasGammaSector11Aside = new TCanvas("canvasGammaSector11Aside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector11Aside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector11Aside;
		histoGSector11Aside = new TH2F("histoGSector11Aside","histoGSector11Aside",2700,167900,170600,10000,-100,100);
		histoGSector11Aside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector11Aside, "Run number","N_{#gamma}^{A} sector 11/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector11Aside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		TPaveText *sector11 = new TPaveText(0.15,0.82,0.3,0.88,"NDC");
		sector11->SetTextSize(0.04);
		sector11->SetLineColor(0);
		sector11->SetFillStyle(0);
		sector11->SetFillColor(0);
		sector11->SetShadowColor(0);
		sector11->SetTextAlign(12);
		sector11->AddText("Sector 11");
		sector11->Draw(); 

		PhotonSector11EtaPos_0010_Data->Draw("psame");
		PhotonSector11EtaPos_1020_Data->Draw("psame");
		PhotonSector11EtaPos_2040_Data->Draw("psame");
		PhotonSector11EtaPos_4060_Data->Draw("psame");
	//	PhotonSector11EtaPos6080_Data->Draw("psame");
		PhotonSector11EtaPos_0010_MC->Draw("psame");
		PhotonSector11EtaPos_1020_MC->Draw("psame");
		PhotonSector11EtaPos_2040_MC->Draw("psame");
		PhotonSector11EtaPos_4060_MC->Draw("psame");
	//	PhotonSector11EtaPos6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector11Aside->SaveAs(Form("%s/GammaSector11_Aside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector11Aside;


		//Sector12
		TCanvas* canvasGammaSector12Aside = new TCanvas("canvasGammaSector12Aside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector12Aside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector12Aside;
		histoGSector12Aside = new TH2F("histoGSector12Aside","histoGSector12Aside",2700,167900,170600,10000,-100,100);
		histoGSector12Aside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector12Aside, "Run number","N_{#gamma}^{A} sector 12/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector12Aside->Draw("copy");
			DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		TPaveText *sector12 = new TPaveText(0.15,0.82,0.3,0.88,"NDC");
		sector12->SetTextSize(0.04);
		sector12->SetLineColor(0);
		sector12->SetFillStyle(0);
		sector12->SetFillColor(0);
		sector12->SetShadowColor(0);
		sector12->SetTextAlign(12);
		sector12->AddText("Sector 12");
		sector12->Draw();

		PhotonSector12EtaPos_0010_Data->Draw("psame");
		PhotonSector12EtaPos_1020_Data->Draw("psame");
		PhotonSector12EtaPos_2040_Data->Draw("psame");
		PhotonSector12EtaPos_4060_Data->Draw("psame");
	//	PhotonSector12EtaPos6080_Data->Draw("psame");
		PhotonSector12EtaPos_0010_MC->Draw("psame");
		PhotonSector12EtaPos_1020_MC->Draw("psame");
		PhotonSector12EtaPos_2040_MC->Draw("psame");
		PhotonSector12EtaPos_4060_MC->Draw("psame");
	//	PhotonSector12EtaPos6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector12Aside->SaveAs(Form("%s/GammaSector12_Aside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector12Aside;


		//Sector13
		TCanvas* canvasGammaSector13Aside = new TCanvas("canvasGammaSector13Aside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector13Aside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector13Aside;
		histoGSector13Aside = new TH2F("histoGSector13Aside","histoGSector13Aside",2700,167900,170600,10000,-100,100);
		histoGSector13Aside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector13Aside, "Run number","N_{#gamma}^{A} sector 13/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector13Aside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		TPaveText *sector13 = new TPaveText(0.15,0.82,0.3,0.88,"NDC");
		sector13->SetTextSize(0.04);
		sector13->SetLineColor(0);
		sector13->SetFillStyle(0);
		sector13->SetFillColor(0);
		sector13->SetShadowColor(0);
		sector13->SetTextAlign(12);
		sector13->AddText("Sector 13");
		sector13->Draw(); 

		PhotonSector13EtaPos_0010_Data->Draw("psame");
		PhotonSector13EtaPos_1020_Data->Draw("psame");
		PhotonSector13EtaPos_2040_Data->Draw("psame");
		PhotonSector13EtaPos_4060_Data->Draw("psame");
	//	PhotonSector13EtaPos6080_Data->Draw("psame");
		PhotonSector13EtaPos_0010_MC->Draw("psame");
		PhotonSector13EtaPos_1020_MC->Draw("psame");
		PhotonSector13EtaPos_2040_MC->Draw("psame");
		PhotonSector13EtaPos_4060_MC->Draw("psame");
	//	PhotonSector13EtaPos6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();

		canvasGammaSector13Aside->SaveAs(Form("%s/GammaSector13_Aside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector13Aside;


		//Sector14
		TCanvas* canvasGammaSector14Aside = new TCanvas("canvasGammaSector14Aside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector14Aside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector14Aside;
		histoGSector14Aside = new TH2F("histoGSector14Aside","histoGSector14Aside",2700,167900,170600,10000,-100,100);
		histoGSector14Aside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector14Aside, "Run number","N_{#gamma}^{A} sector 14/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector14Aside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		TPaveText *sector14 = new TPaveText(0.15,0.82,0.3,0.88,"NDC");
		sector14->SetTextSize(0.04);
		sector14->SetLineColor(0);
		sector14->SetFillStyle(0);
		sector14->SetFillColor(0);
		sector14->SetShadowColor(0);
		sector14->SetTextAlign(12);
		sector14->AddText("Sector 14");
		sector14->Draw(); 

		PhotonSector14EtaPos_0010_Data->Draw("psame");
		PhotonSector14EtaPos_1020_Data->Draw("psame");
		PhotonSector14EtaPos_2040_Data->Draw("psame");
		PhotonSector14EtaPos_4060_Data->Draw("psame");
	//	PhotonSector14EtaPos6080_Data->Draw("psame");
		PhotonSector14EtaPos_0010_MC->Draw("psame");
		PhotonSector14EtaPos_1020_MC->Draw("psame");
		PhotonSector14EtaPos_2040_MC->Draw("psame");
		PhotonSector14EtaPos_4060_MC->Draw("psame");
	//	PhotonSector14EtaPos6080_MC->Draw("psame");
		
		MC->Draw();
		legendGAllSectors->Draw();

		canvasGammaSector14Aside->SaveAs(Form("%s/GammaSector14_Aside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector14Aside;

		//Sector15
		TCanvas* canvasGammaSector15Aside = new TCanvas("canvasGammaSector15Aside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector15Aside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector15Aside;
		histoGSector15Aside = new TH2F("histoGSector15Aside","histoGSector15Aside",2700,167900,170600,10000,-100,100);
		histoGSector15Aside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector15Aside, "Run number","N_{#gamma}^{A} sector 15/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector15Aside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		TPaveText *sector15 = new TPaveText(0.15,0.82,0.3,0.88,"NDC");
		sector15->SetTextSize(0.04);
		sector15->SetLineColor(0);
		sector15->SetFillStyle(0);
		sector15->SetFillColor(0);
		sector15->SetShadowColor(0);
		sector15->SetTextAlign(12);
		sector15->AddText("Sector 15");
		sector15->Draw();

		PhotonSector15EtaPos_0010_Data->Draw("psame");
		PhotonSector15EtaPos_1020_Data->Draw("psame");
		PhotonSector15EtaPos_2040_Data->Draw("psame");
		PhotonSector15EtaPos_4060_Data->Draw("psame");
	//	PhotonSector15EtaPos6080_Data->Draw("psame");
		PhotonSector15EtaPos_0010_MC->Draw("psame");
		PhotonSector15EtaPos_1020_MC->Draw("psame");
		PhotonSector15EtaPos_2040_MC->Draw("psame");
		PhotonSector15EtaPos_4060_MC->Draw("psame");
	//	PhotonSector15EtaPos6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();

		canvasGammaSector15Aside->SaveAs(Form("%s/GammaSector15_Aside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector15Aside;


		//Sector16
		TCanvas* canvasGammaSector16Aside = new TCanvas("canvasGammaSector16Aside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector16Aside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector16Aside;
		histoGSector16Aside = new TH2F("histoGSector16Aside","histoGSector16Aside",2700,167900,170600,10000,-100,100);
		histoGSector16Aside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector16Aside, "Run number","N_{#gamma}^{A} sector 16/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector16Aside->Draw("copy");
			DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		TPaveText *sector16 = new TPaveText(0.15,0.82,0.3,0.88,"NDC");
		sector16->SetTextSize(0.04);
		sector16->SetLineColor(0);
		sector16->SetFillStyle(0);
		sector16->SetFillColor(0);
		sector16->SetShadowColor(0);
		sector16->SetTextAlign(12);
		sector16->AddText("Sector 16");
		sector16->Draw();

		PhotonSector16EtaPos_0010_Data->Draw("psame");
		PhotonSector16EtaPos_1020_Data->Draw("psame");
		PhotonSector16EtaPos_2040_Data->Draw("psame");
		PhotonSector16EtaPos_4060_Data->Draw("psame");
	//	PhotonSector16EtaPos6080_Data->Draw("psame");
		PhotonSector16EtaPos_0010_MC->Draw("psame");
		PhotonSector16EtaPos_1020_MC->Draw("psame");
		PhotonSector16EtaPos_2040_MC->Draw("psame");
		PhotonSector16EtaPos_4060_MC->Draw("psame");
	//	PhotonSector16EtaPos6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();

		canvasGammaSector16Aside->SaveAs(Form("%s/GammaSector16_Aside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector16Aside;

		//Sector17
		TCanvas* canvasGammaSector17Aside = new TCanvas("canvasGammaSector17Aside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector17Aside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector17Aside;
		histoGSector17Aside = new TH2F("histoGSector17Aside","histoGSector17Aside",2700,167900,170600,10000,-100,100);
		histoGSector17Aside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector17Aside, "Run number","N_{#gamma}^{A} sector 17/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector17Aside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		TPaveText *sector17 = new TPaveText(0.15,0.82,0.3,0.88,"NDC");
		sector17->SetTextSize(0.04);
		sector17->SetLineColor(0);
		sector17->SetFillStyle(0);
		sector17->SetFillColor(0);
		sector17->SetShadowColor(0);
		sector17->SetTextAlign(12);
		sector17->AddText("Sector 17");
		sector17->Draw();

		PhotonSector17EtaPos_0010_Data->Draw("psame");
		PhotonSector17EtaPos_1020_Data->Draw("psame");
		PhotonSector17EtaPos_2040_Data->Draw("psame");
		PhotonSector17EtaPos_4060_Data->Draw("psame");
	//	PhotonSector17EtaPos6080_Data->Draw("psame"); 
		PhotonSector17EtaPos_0010_MC->Draw("psame");
		PhotonSector17EtaPos_1020_MC->Draw("psame");
		PhotonSector17EtaPos_2040_MC->Draw("psame");
		PhotonSector17EtaPos_4060_MC->Draw("psame");
	//	PhotonSector17EtaPos6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();

		canvasGammaSector17Aside->SaveAs(Form("%s/GammaSector17_Aside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector17Aside;




	///////////////////////////////////////////////////////// C side ///////////////////////////////////////////////////////////////////////
		TCanvas* canvasGammaAllSectorsCside = new TCanvas("canvasGammaAllSectorsCside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaAllSectorsCside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGAllSectorsCside;
		histoGAllSectorsCside = new TH2F("histoGAllSectorsCside","histoGAllSectorsCside",2700,167900,170600,10000,-100,100);
		histoGAllSectorsCside->GetYaxis()->SetRangeUser(0.,22.);
		SetStyleHistoTH2ForGraphs(histoGAllSectorsCside, "Run number","N_{#gamma}^{C} all sectors/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGAllSectorsCside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		sectorAll->Draw(); 

		PhotonAllSectorsEtaNeg_0010_Data->Draw("psame");
		PhotonAllSectorsEtaNeg_1020_Data->Draw("psame");
		PhotonAllSectorsEtaNeg_2040_Data->Draw("psame");
		PhotonAllSectorsEtaNeg_4060_Data->Draw("psame");
	//	PhotonAllSectorsEtaNeg6080_Data->Draw("psame");
		PhotonAllSectorsEtaNeg_0010_MC->Draw("psame");
		PhotonAllSectorsEtaNeg_1020_MC->Draw("psame");
		PhotonAllSectorsEtaNeg_2040_MC->Draw("psame");
		PhotonAllSectorsEtaNeg_4060_MC->Draw("psame");
	//	PhotonAllSectorsEtaNeg6080_MC->Draw("psame");
		
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaAllSectorsCside->SaveAs(Form("%s/GammaAllSectors_Cside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaAllSectorsCside;

		
		//Sector0
		TCanvas* canvasGammaSector0Cside = new TCanvas("canvasGammaSector0Cside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector0Cside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector0Cside;
		histoGSector0Cside = new TH2F("histoGSector0Cside","histoGSector0Cside",2700,167900,170600,10000,-100,100);
		histoGSector0Cside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector0Cside, "Run number","N_{#gamma}^{C} sector 0/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector0Cside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		sector0->Draw(); 

		PhotonSector0EtaNeg_0010_Data->Draw("psame");
		PhotonSector0EtaNeg_1020_Data->Draw("psame");
		PhotonSector0EtaNeg_2040_Data->Draw("psame");
		PhotonSector0EtaNeg_4060_Data->Draw("psame");
	//	PhotonSector0EtaNeg6080_Data->Draw("psame");
		PhotonSector0EtaNeg_0010_MC->Draw("psame");
		PhotonSector0EtaNeg_1020_MC->Draw("psame");
		PhotonSector0EtaNeg_2040_MC->Draw("psame");
		PhotonSector0EtaNeg_4060_MC->Draw("psame");
	//	PhotonSector0EtaNeg6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();

		canvasGammaSector0Cside->SaveAs(Form("%s/GammaSector0_Cside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector0Cside;

		//Sector1
		TCanvas* canvasGammaSector1Cside = new TCanvas("canvasGammaSector1Cside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector1Cside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector1Cside;
		histoGSector1Cside = new TH2F("histoGSector1Cside","histoGSector1Cside",2700,167900,170600,10000,-100,100);
		histoGSector1Cside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector1Cside, "Run number","N_{#gamma}^{C} sector 1/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector1Cside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		sector1->Draw(); 

		PhotonSector1EtaNeg_0010_Data->Draw("psame");
		PhotonSector1EtaNeg_1020_Data->Draw("psame");
		PhotonSector1EtaNeg_2040_Data->Draw("psame");
		PhotonSector1EtaNeg_4060_Data->Draw("psame");
	//	PhotonSector1EtaNeg6080_Data->Draw("psame");
		PhotonSector1EtaNeg_0010_MC->Draw("psame");
		PhotonSector1EtaNeg_1020_MC->Draw("psame");
		PhotonSector1EtaNeg_2040_MC->Draw("psame");
		PhotonSector1EtaNeg_4060_MC->Draw("psame");
	//	PhotonSector1EtaNeg6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();

		canvasGammaSector1Cside->SaveAs(Form("%s/GammaSector1_Cside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector1Cside;

		//Sector2
		TCanvas* canvasGammaSector2Cside = new TCanvas("canvasGammaSector2Cside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector2Cside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector2Cside;
		histoGSector2Cside = new TH2F("histoGSector2Cside","histoGSector2Cside",2700,167900,170600,10000,-100,100);
		histoGSector2Cside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector2Cside, "Run number","N_{#gamma}^{C} sector 2/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector2Cside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		sector2->Draw(); 

		PhotonSector2EtaNeg_0010_Data->Draw("psame");
		PhotonSector2EtaNeg_1020_Data->Draw("psame");
		PhotonSector2EtaNeg_2040_Data->Draw("psame");
		PhotonSector2EtaNeg_4060_Data->Draw("psame");
	//	PhotonSector2EtaNeg6080_Data->Draw("psame");
		PhotonSector2EtaNeg_0010_MC->Draw("psame");
		PhotonSector2EtaNeg_1020_MC->Draw("psame");
		PhotonSector2EtaNeg_2040_MC->Draw("psame");
		PhotonSector2EtaNeg_4060_MC->Draw("psame");
	//	PhotonSector2EtaNeg6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();

		canvasGammaSector2Cside->SaveAs(Form("%s/GammaSector2_Cside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector2Cside;

		//Sector3
		TCanvas* canvasGammaSector3Cside = new TCanvas("canvasGammaSector3Cside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector3Cside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector3Cside;
		histoGSector3Cside = new TH2F("histoGSector3Cside","histoGSector3",2700,167900,170600,10000,-100,100);
		histoGSector3Cside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector3Cside, "Run number","N_{#gamma}^{C} sector 3/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector3Cside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		sector3->Draw(); 

		PhotonSector3EtaNeg_0010_Data->Draw("x x psame");
		PhotonSector3EtaNeg_1020_Data->Draw("x x psame");
		PhotonSector3EtaNeg_2040_Data->Draw("x x psame");
		PhotonSector3EtaNeg_4060_Data->Draw("x x psame");
	//	PhotonSector3EtaNeg6080_Data->Draw("psame");
		PhotonSector3EtaNeg_0010_MC->Draw("x x psame");
		PhotonSector3EtaNeg_1020_MC->Draw("x x psame");
		PhotonSector3EtaNeg_2040_MC->Draw("x x psame");
		PhotonSector3EtaNeg_4060_MC->Draw("x x psame");
	//	PhotonSector3EtaNeg6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();

		canvasGammaSector3Cside->SaveAs(Form("%s/GammaSector3_Cside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector3Cside;


		//Sector4
		TCanvas* canvasGammaSector4Cside = new TCanvas("canvasGammaSector4Cside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector4Cside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector4Cside;
		histoGSector4Cside = new TH2F("histoGSector4Cside","histoGSector4Cside",2700,167900,170600,10000,-100,100);
		histoGSector4Cside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector4Cside, "Run number","N_{#gamma}^{C} sector 4/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector4Cside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		sector4->Draw(); 

		PhotonSector4EtaNeg_0010_Data->Draw("psame");
		PhotonSector4EtaNeg_1020_Data->Draw("psame");
		PhotonSector4EtaNeg_2040_Data->Draw("psame");
		PhotonSector4EtaNeg_4060_Data->Draw("psame");
	//	PhotonSector4EtaNeg6080_Data->Draw("psame");
		PhotonSector4EtaNeg_0010_MC->Draw("psame");
		PhotonSector4EtaNeg_1020_MC->Draw("psame");
		PhotonSector4EtaNeg_2040_MC->Draw("psame");
		PhotonSector4EtaNeg_4060_MC->Draw("psame");
	//	PhotonSector4EtaNeg6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector4Cside->SaveAs(Form("%s/GammaSector4_Cside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector4Cside;


		//Sector5
		TCanvas* canvasGammaSector5Cside = new TCanvas("canvasGammaSector5Cside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector5Cside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector5Cside;
		histoGSector5Cside = new TH2F("histoGSector5Cside","histoGSector5Cside",2700,167900,170600,10000,-100,100);
		histoGSector5Cside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector5Cside, "Run number","N_{#gamma}^{C} sector 5/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector5Cside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		sector5->Draw(); 

		PhotonSector5EtaNeg_0010_Data->Draw("psame");
		PhotonSector5EtaNeg_1020_Data->Draw("psame");
		PhotonSector5EtaNeg_2040_Data->Draw("psame");
		PhotonSector5EtaNeg_4060_Data->Draw("psame");
	//	PhotonSector5EtaNeg6080_Data->Draw("psame");
		PhotonSector5EtaNeg_0010_MC->Draw("psame");
		PhotonSector5EtaNeg_1020_MC->Draw("psame");
		PhotonSector5EtaNeg_2040_MC->Draw("psame");
		PhotonSector5EtaNeg_4060_MC->Draw("psame");
	//	PhotonSector5EtaNeg6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector5Cside->SaveAs(Form("%s/GammaSector5_Cside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector5Cside;


		//Sector6
		TCanvas* canvasGammaSector6Cside = new TCanvas("canvasGammaSector6Cside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector6Cside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector6Cside;
		histoGSector6Cside = new TH2F("histoGSector6Cside","histoGSector6Cside",2700,167900,170600,10000,-100,100);
		histoGSector6Cside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector6Cside, "Run number","N_{#gamma}^{C} sector 6/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector6Cside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		sector6->Draw();

		PhotonSector6EtaNeg_0010_Data->Draw("psame");
		PhotonSector6EtaNeg_1020_Data->Draw("psame");
		PhotonSector6EtaNeg_2040_Data->Draw("psame");
		PhotonSector6EtaNeg_4060_Data->Draw("psame");
	//	PhotonSector6EtaNeg6080_Data->Draw("psame");
		PhotonSector6EtaNeg_0010_MC->Draw("psame");
		PhotonSector6EtaNeg_1020_MC->Draw("psame");
		PhotonSector6EtaNeg_2040_MC->Draw("psame");
		PhotonSector6EtaNeg_4060_MC->Draw("psame");
	//	PhotonSector6EtaNeg6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector6Cside->SaveAs(Form("%s/GammaSector6_Cside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector6Cside;


		//Sector7
		TCanvas* canvasGammaSector7Cside = new TCanvas("canvasGammaSector7Cside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector7Cside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector7Cside;
		histoGSector7Cside = new TH2F("histoGSector7Cside","histoGSector7Cside",2700,167900,170600,10000,-100,100);
		histoGSector7Cside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector7Cside, "Run number","N_{#gamma}^{C} sector 7/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector7Cside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		sector7->Draw(); 

		PhotonSector7EtaNeg_0010_Data->Draw("psame");
		PhotonSector7EtaNeg_1020_Data->Draw("psame");
		PhotonSector7EtaNeg_2040_Data->Draw("psame");
		PhotonSector7EtaNeg_4060_Data->Draw("psame");
	//	PhotonSector7EtaNeg6080_Data->Draw("psame");
		PhotonSector7EtaNeg_0010_MC->Draw("psame");
		PhotonSector7EtaNeg_1020_MC->Draw("psame");
		PhotonSector7EtaNeg_2040_MC->Draw("psame");
		PhotonSector7EtaNeg_4060_MC->Draw("psame");
	//	PhotonSector7EtaNeg6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector7Cside->SaveAs(Form("%s/GammaSector7_Cside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector7Cside;


		//Sector8
		TCanvas* canvasGammaSector8Cside = new TCanvas("canvasGammaSector8Cside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector8Cside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector8Cside;
		histoGSector8Cside = new TH2F("histoGSector8Cside","histoGSector8Cside",2700,167900,170600,10000,-100,100);
		histoGSector8Cside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector8Cside, "Run number","N_{#gamma}^{C} sector 8/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector8Cside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		sector8->Draw();

		PhotonSector8EtaNeg_0010_Data->Draw("psame");
		PhotonSector8EtaNeg_1020_Data->Draw("psame");
		PhotonSector8EtaNeg_2040_Data->Draw("psame");
		PhotonSector8EtaNeg_4060_Data->Draw("psame");
	//	PhotonSector8EtaNeg6080_Data->Draw("psame");
		PhotonSector8EtaNeg_0010_MC->Draw("psame");
		PhotonSector8EtaNeg_1020_MC->Draw("psame");
		PhotonSector8EtaNeg_2040_MC->Draw("psame");
		PhotonSector8EtaNeg_4060_MC->Draw("psame");
	//	PhotonSector8EtaNeg6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector8Cside->SaveAs(Form("%s/GammaSector8_Cside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector8Cside;


		//Sector9
		TCanvas* canvasGammaSector9Cside = new TCanvas("canvasGammaSector9Cside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector9Cside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector9Cside;
		histoGSector9Cside = new TH2F("histoGSector9Cside","histoGSector9Cside",2700,167900,170600,10000,-100,100);
		histoGSector9Cside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector9Cside, "Run number","N_{#gamma}^{C} sector 9/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector9Cside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		sector9->Draw(); 

		PhotonSector9EtaNeg_0010_Data->Draw("psame");
		PhotonSector9EtaNeg_1020_Data->Draw("psame");
		PhotonSector9EtaNeg_2040_Data->Draw("psame");
		PhotonSector9EtaNeg_4060_Data->Draw("psame");
	//	PhotonSector9EtaNeg6080_Data->Draw("psame");
		PhotonSector9EtaNeg_0010_MC->Draw("psame");
		PhotonSector9EtaNeg_1020_MC->Draw("psame");
		PhotonSector9EtaNeg_2040_MC->Draw("psame");
		PhotonSector9EtaNeg_4060_MC->Draw("psame");
	//	PhotonSector9EtaNeg6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector9Cside->SaveAs(Form("%s/GammaSector9_Cside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector9Cside;


		//Sector10
		TCanvas* canvasGammaSector10Cside = new TCanvas("canvasGammaSector10Cside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector10Cside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector10Cside;
		histoGSector10Cside = new TH2F("histoGSector10Cside","histoGSector10Cside",2700,167900,170600,10000,-100,100);
		histoGSector10Cside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector10Cside, "Run number","N_{#gamma}^{C} sector 10/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector10Cside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		sector10->Draw();

		PhotonSector10EtaNeg_0010_Data->Draw("psame");
		PhotonSector10EtaNeg_1020_Data->Draw("psame");
		PhotonSector10EtaNeg_2040_Data->Draw("psame");
		PhotonSector10EtaNeg_4060_Data->Draw("psame");
	//	PhotonSector10EtaNeg6080_Data->Draw("psame");
		PhotonSector10EtaNeg_0010_MC->Draw("psame");
		PhotonSector10EtaNeg_1020_MC->Draw("psame");
		PhotonSector10EtaNeg_2040_MC->Draw("psame");
		PhotonSector10EtaNeg_4060_MC->Draw("psame");
	//	PhotonSector10EtaNeg6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector10Cside->SaveAs(Form("%s/GammaSector10_Cside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector10Cside;

		//Sector11
		TCanvas* canvasGammaSector11Cside = new TCanvas("canvasGammaSector11Cside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector11Cside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector11Cside;
		histoGSector11Cside = new TH2F("histoGSector11Cside","histoGSector11Cside",2700,167900,170600,10000,-100,100);
		histoGSector11Cside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector11Cside, "Run number","N_{#gamma}^{C} sector 11/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector11Cside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		sector11->Draw(); 

		PhotonSector11EtaNeg_0010_Data->Draw("psame");
		PhotonSector11EtaNeg_1020_Data->Draw("psame");
		PhotonSector11EtaNeg_2040_Data->Draw("psame");
		PhotonSector11EtaNeg_4060_Data->Draw("psame");
	//	PhotonSector11EtaNeg6080_Data->Draw("psame");
		PhotonSector11EtaNeg_0010_MC->Draw("psame");
		PhotonSector11EtaNeg_1020_MC->Draw("psame");
		PhotonSector11EtaNeg_2040_MC->Draw("psame");
		PhotonSector11EtaNeg_4060_MC->Draw("psame");
	//	PhotonSector11EtaNeg6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector11Cside->SaveAs(Form("%s/GammaSector11_Cside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector11Cside;

		//Sector12
		TCanvas* canvasGammaSector12Cside = new TCanvas("canvasGammaSector12Cside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector12Cside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector12Cside;
		histoGSector12Cside = new TH2F("histoGSector12Cside","histoGSector12Cside",2700,167900,170600,10000,-100,100);
		histoGSector12Cside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector12Cside, "Run number","N_{#gamma}^{C} sector 12/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector12Cside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		sector12->Draw(); 

		PhotonSector12EtaNeg_0010_Data->Draw("psame");
		PhotonSector12EtaNeg_1020_Data->Draw("psame");
		PhotonSector12EtaNeg_2040_Data->Draw("psame");
		PhotonSector12EtaNeg_4060_Data->Draw("psame");
	//	PhotonSector12EtaNeg6080_Data->Draw("psame");
		PhotonSector12EtaNeg_0010_MC->Draw("psame");
		PhotonSector12EtaNeg_1020_MC->Draw("psame");
		PhotonSector12EtaNeg_2040_MC->Draw("psame");
		PhotonSector12EtaNeg_4060_MC->Draw("psame");
	//	PhotonSector12EtaNeg6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector12Cside->SaveAs(Form("%s/GammaSector12_Cside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector12Cside;


		//Sector13
		TCanvas* canvasGammaSector13Cside = new TCanvas("canvasGammaSector13Cside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector13Cside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector13Cside;
		histoGSector13Cside = new TH2F("histoGSector13Cside","histoGSector13Cside",2700,167900,170600,10000,-100,100);
		histoGSector13Cside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector13Cside, "Run number","N_{#gamma}^{C} sector 13/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector13Cside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		sector13->Draw(); 

		PhotonSector13EtaNeg_0010_Data->Draw("psame");
		PhotonSector13EtaNeg_1020_Data->Draw("psame");
		PhotonSector13EtaNeg_2040_Data->Draw("psame");
		PhotonSector13EtaNeg_4060_Data->Draw("psame");
	//	PhotonSector13EtaNeg6080_Data->Draw("psame");
		PhotonSector13EtaNeg_0010_MC->Draw("psame");
		PhotonSector13EtaNeg_1020_MC->Draw("psame");
		PhotonSector13EtaNeg_2040_MC->Draw("psame");
		PhotonSector13EtaNeg_4060_MC->Draw("psame");
	//	PhotonSector13EtaNeg6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector13Cside->SaveAs(Form("%s/GammaSector13_Cside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector13Cside;


		//Sector14
		TCanvas* canvasGammaSector14Cside = new TCanvas("canvasGammaSector14Cside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector14Cside,0.1, 0.03, 0.08, 0.08);

		sector14->Draw(); 	

		TH2F * histoGSector14Cside;
		histoGSector14Cside = new TH2F("histoGSector14Cside","histoGSector14Cside",2700,167900,170600,10000,-100,100);
		histoGSector14Cside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector14Cside, "Run number","N_{#gamma}^{C} sector 14/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector14Cside->Draw("copy");
			DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		PhotonSector14EtaNeg_0010_Data->Draw("psame");
		PhotonSector14EtaNeg_1020_Data->Draw("psame");
		PhotonSector14EtaNeg_2040_Data->Draw("psame");
		PhotonSector14EtaNeg_4060_Data->Draw("psame");
	//	PhotonSector14EtaNeg6080_Data->Draw("psame");
		PhotonSector14EtaNeg_0010_MC->Draw("psame");
		PhotonSector14EtaNeg_1020_MC->Draw("psame");
		PhotonSector14EtaNeg_2040_MC->Draw("psame");	
		PhotonSector14EtaNeg_4060_MC->Draw("psame");
	//	PhotonSector14EtaNeg6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector14Cside->SaveAs(Form("%s/GammaSector14_Cside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector14Cside;

		//Sector15
		TCanvas* canvasGammaSector15Cside = new TCanvas("canvasGammaSector15Cside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector15Cside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector15Cside;
		histoGSector15Cside = new TH2F("histoGSector15Cside","histoGSector15Cside",2700,167900,170600,10000,-100,100);
		histoGSector15Cside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector15Cside, "Run number","N_{#gamma}^{C} sector 15/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector15Cside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		sector15->Draw();

		PhotonSector15EtaNeg_0010_Data->Draw("psame");
		PhotonSector15EtaNeg_1020_Data->Draw("psame");
		PhotonSector15EtaNeg_2040_Data->Draw("psame");
		PhotonSector15EtaNeg_4060_Data->Draw("psame");
	//	PhotonSector15EtaNeg6080_Data->Draw("psame");
		PhotonSector15EtaNeg_0010_MC->Draw("psame");
		PhotonSector15EtaNeg_1020_MC->Draw("psame");
		PhotonSector15EtaNeg_2040_MC->Draw("psame");
		PhotonSector15EtaNeg_4060_MC->Draw("psame");
	//	PhotonSector15EtaNeg6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector15Cside->SaveAs(Form("%s/GammaSector15_Cside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector15Cside;


		//Sector16
		TCanvas* canvasGammaSector16Cside = new TCanvas("canvasGammaSector16Cside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector16Cside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector16Cside;
		histoGSector16Cside = new TH2F("histoGSector16Cside","histoGSector16Cside",2700,167900,170600,10000,-100,100);
		histoGSector16Cside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector16Cside, "Run number","N_{#gamma}^{C} sector 16/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector16Cside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		sector16->Draw(); 

		PhotonSector16EtaNeg_0010_Data->Draw("psame");
		PhotonSector16EtaNeg_1020_Data->Draw("psame");
		PhotonSector16EtaNeg_2040_Data->Draw("psame");
		PhotonSector16EtaNeg_4060_Data->Draw("psame");
	//	PhotonSector16EtaNeg6080_Data->Draw("psame");
		PhotonSector16EtaNeg_0010_MC->Draw("psame");
		PhotonSector16EtaNeg_1020_MC->Draw("psame");
		PhotonSector16EtaNeg_2040_MC->Draw("psame");
		PhotonSector16EtaNeg_4060_MC->Draw("psame");
	//	PhotonSector16EtaNeg6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector16Cside->SaveAs(Form("%s/GammaSector16_Cside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector16Cside;


		//Sector17
		TCanvas* canvasGammaSector17Cside = new TCanvas("canvasGammaSector17Cside","",200,10,2000,1000);  // gives the page size
		DrawGammaCanvasSettings(canvasGammaSector17Cside,0.1, 0.03, 0.08, 0.08);
		
		TH2F * histoGSector17Cside;
		histoGSector17Cside = new TH2F("histoGSector17Cside","histoGSector17Cside",2700,167900,170600,10000,-100,100);
		histoGSector17Cside->GetYaxis()->SetRangeUser(0.,1.5);
		SetStyleHistoTH2ForGraphs(histoGSector17Cside, "Run number","N_{#gamma}^{C} sector 17/N_{evt,run}",0.03,0.04, 0.03,0.04, 0.75,0.95);
		histoGSector17Cside->Draw("copy");
		DrawLabelsEvents(floatLocationLeft[0],floatLocationLeft[1],floatLocationLeft[2], floatLocationLeft[3], collisionSystem,  textPeriod, optMCGenerator);
		sector17->Draw();

		PhotonSector17EtaNeg_0010_Data->Draw("psame");
		PhotonSector17EtaNeg_1020_Data->Draw("psame");
		PhotonSector17EtaNeg_2040_Data->Draw("psame");
		PhotonSector17EtaNeg_4060_Data->Draw("psame");
	//	PhotonSector17EtaNeg6080_Data->Draw("psame"); 
		PhotonSector17EtaNeg_0010_MC->Draw("psame");
		PhotonSector17EtaNeg_1020_MC->Draw("psame");
		PhotonSector17EtaNeg_2040_MC->Draw("psame");
		PhotonSector17EtaNeg_4060_MC->Draw("psame");
	//	PhotonSector17EtaNeg6080_MC->Draw("psame");
		MC->Draw();
		legendGAllSectors->Draw();
		canvasGammaSector17Cside->SaveAs(Form("%s/GammaSector17_Cside.%s",outputDir.Data(),suffix.Data()));
		delete canvasGammaSector17Cside;
	}


}

