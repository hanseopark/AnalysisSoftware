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
#include "TTree.h"
#include "TMinuit.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TVectorT.h"

typedef TVectorT<double> TVectorD;


void BuildHistogramsForSecHadInteractions(TString fileName = "fbock_SecHadIntTree.root", TString CutSelection="214507" , Bool_t kMC=0, Bool_t merge=0, TString nameOutput = "MappingDetailedSecHadInt"){
	//Reset ROOT and connect tree file
	TString minSecHadTrackCutNumber = CutSelection(0,1);
	TString etaCutNumber = CutSelection(1,1);
	TString minPtCutNumber = CutSelection(2,1);
	TString chi2CutNumber = CutSelection(3,1);
	TString err2DCutNumber = CutSelection(4,1);
	TString dcaTrackCutNumber = CutSelection(5,1);
	
	Int_t minSecHadTrackCut=0;
	if (minSecHadTrackCutNumber.CompareTo("2") == 0){
		minSecHadTrackCut = 2;
 	} else if (minSecHadTrackCutNumber.CompareTo("3") == 0 ){
		minSecHadTrackCut = 3;
 	}else if (minSecHadTrackCutNumber.CompareTo("4") == 0 ){
		minSecHadTrackCut = 4;
 	}
	cout << "Minimum number of tracks per sec Vtx: " << minSecHadTrackCut << endl;
	
	Double_t etaMinCut = -0.1;
	Double_t etaMaxCut = 0.;
	if (etaCutNumber.CompareTo("0") == 0){
		etaMaxCut = 0.9;
 	} else if (etaCutNumber.CompareTo("1") == 0 ){
		etaMaxCut = 0.1;
 	}else if (etaCutNumber.CompareTo("2") == 0 ){
		etaMaxCut = 1.4;
		etaMinCut = 0.9;
 	}else if (etaCutNumber.CompareTo("3") == 0 ){
		etaMaxCut = 10.;
		etaMinCut = 1.4;
 	} else if (etaCutNumber.CompareTo("4") == 0 ){
		etaMaxCut = 2.5;
		etaMinCut = 1.8;
 	} else if (etaCutNumber.CompareTo("5") == 0 ){
		etaMaxCut = 10.;
		etaMinCut = 2.5;
 	} else if (etaCutNumber.CompareTo("6") == 0 ){
		etaMaxCut = 10.;
		etaMinCut = -0.1;
 	}
 	if (etaMinCut != -0.1){
		cout << "Eta range: " <<   etaMinCut  <<  " < |eta| < " << etaMaxCut << endl;
	} else {
		cout << "Eta range: " <<  "|eta| < " << etaMaxCut << endl;
	}
	
	Double_t minPtTrackCut = 0.0;
	if (minPtCutNumber.CompareTo("0") == 0){
		minPtTrackCut = 0.05;
 	} else if (minPtCutNumber.CompareTo("1") == 0 ){
		minPtTrackCut = 0.1;
 	} else if (minPtCutNumber.CompareTo("2") == 0 ){
		minPtTrackCut = 0.15;
	} else if (minPtCutNumber.CompareTo("3") == 0 ){
		minPtTrackCut = 0.2;
 	} else if (minPtCutNumber.CompareTo("4") == 0 ){
		minPtTrackCut = 0.3;
	}
	cout<< "Min pt for tracks: " << minPtTrackCut << endl;
	
	Double_t chi2PerDOFCut = 0.0;
	if (chi2CutNumber.CompareTo("0") == 0){
		chi2PerDOFCut = 3.;
 	} else if (chi2CutNumber.CompareTo("1") == 0 ){
		chi2PerDOFCut = 5.;
 	} else if (chi2CutNumber.CompareTo("2") == 0 ){
		chi2PerDOFCut = 10.;
	} else if (chi2CutNumber.CompareTo("3") == 0 ){
		chi2PerDOFCut = 15.;
 	} else if (chi2CutNumber.CompareTo("4") == 0 ){
		chi2PerDOFCut = 20.;
	}
	cout<< "Chi2 per DOF Cut: " << chi2PerDOFCut << endl;
	
	Double_t err2DCut = 0.0;
	if (err2DCutNumber.CompareTo("0") == 0){
		err2DCut = 0.5;
 	} else if (err2DCutNumber.CompareTo("1") == 0 ){
		err2DCut = 0.8;
 	} else if (err2DCutNumber.CompareTo("2") == 0 ){
		err2DCut = 1.;
	} else if (err2DCutNumber.CompareTo("3") == 0 ){
		err2DCut = 1.5;
 	} else if (err2DCutNumber.CompareTo("4") == 0 ){
		err2DCut = 2.0;
	}
	cout<< "Err2D Cut: " << err2DCut << endl;
	
	Double_t dcaTrackCut = 0.; 
	Double_t dcaTrackAbsCut = 0.; 
	if (dcaTrackCutNumber.CompareTo("0") == 0){
		dcaTrackCut = 200.;
		dcaTrackAbsCut = 1.;
 	} else if (dcaTrackCutNumber.CompareTo("1") == 0 ){
		dcaTrackCut = 0.75;
		dcaTrackAbsCut = 1.;
 	} else if (dcaTrackCutNumber.CompareTo("2") == 0 ){
		dcaTrackCut = 0.1;
		dcaTrackAbsCut = 1.;
	} else if (dcaTrackCutNumber.CompareTo("3") == 0 ){
		dcaTrackCut = 0.125;
		dcaTrackAbsCut = 1.;
 	} else if (dcaTrackCutNumber.CompareTo("4") == 0 ){
		dcaTrackCut = 0.15;
		dcaTrackAbsCut = 1.;
	} else if (dcaTrackCutNumber.CompareTo("5") == 0 ){
		dcaTrackCut = 0.2;
		dcaTrackAbsCut = 1.;
	} else if (dcaTrackCutNumber.CompareTo("6") == 0 ){
		dcaTrackCut = 0.5;
		dcaTrackAbsCut = 1.;
	} else if (dcaTrackCutNumber.CompareTo("7") == 0 ){
		dcaTrackCut = 0.75;
		dcaTrackAbsCut = 1.;	
	} else if (dcaTrackCutNumber.CompareTo("8") == 0 ){
		dcaTrackCut = 1.;
		dcaTrackAbsCut = 1.;
	} else if (dcaTrackCutNumber.CompareTo("9") == 0 ){
		dcaTrackCut = 1.5;
		dcaTrackAbsCut = 1.;
	}
	
	cout<< "DCA Cut: " << dcaTrackCut << "\t" << dcaTrackAbsCut << endl;
	
	Double_t rSPD=9.5;
	Double_t rSPDTh=13;
	Double_t rSDD=35;
	Double_t rSSD=55;
	Double_t rTPC=72.5;
	Double_t rHotZoneMin = 5.7;
	Double_t rHotZoneMax = 50.;
	Double_t zHotZoneMin = 45.;
	Double_t zHotZoneMax = 75.;
	Double_t zBeamPipeInner = 30.;	
	
   gROOT->Reset();

   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fileName.Data());
   if (!f) {
      f = new TFile(fileName.Data());
   }
   TTree *Event = (TTree*)gDirectory->Get("Event");

	
	
	//Declaration of leaves types
   Double_t        primVtxZ;
   Char_t          SPDVtx;
   Int_t           nContrVtx;
   Int_t           fMultiplicity;

   // Set branch addresses.
   Event->SetBranchAddress("primVtxZ",&primVtxZ);
   Event->SetBranchAddress("SPDVtx",&SPDVtx);
   Event->SetBranchAddress("nContrVtx",&nContrVtx);
   Event->SetBranchAddress("fMultiplicity",&fMultiplicity);

   ULong64_t nEntriesEvent = Event->GetEntries();


	cout << "number of Events: " << nEntriesEvent<< endl;
	
   TTree *Vertices = (TTree*)gDirectory->Get("Vertices");

	//Declaration of leaves types
   Int_t           nTracksPerVtx;
   Double_t        fX;
   Double_t        fY;
   Double_t        fZ;
   Double_t        fR;
   Double_t        fPhi;
   Double_t        fEta;
   Double_t        fErrX;
   Double_t        fErrY;
   Double_t        fErrZ;
   Double_t        fErr2D;
   Double_t        fErr3D;
   Double_t        fChiPerDOF;
	Int_t				 fMultTrack;
	TVectorD *dcaTrack = new  TVectorD ;
	TVectorD *pTTrack = new  TVectorD ;
	TVectorD *primVtx = new  TVectorD ;
//    Double_t dcaTrack[3];
//    Double_t pTTrack[3];

	
   // Set branch addresses.
   Vertices->SetBranchAddress("nTrackPerVtx",&nTracksPerVtx);
   Vertices->SetBranchAddress("X",&fX);
   Vertices->SetBranchAddress("Y",&fY);
   Vertices->SetBranchAddress("Z",&fZ);
   Vertices->SetBranchAddress("R",&fR);
   Vertices->SetBranchAddress("Phi",&fPhi);
   Vertices->SetBranchAddress("Eta",&fEta);
   Vertices->SetBranchAddress("ErrX",&fErrX);
   Vertices->SetBranchAddress("ErrY",&fErrY);
   Vertices->SetBranchAddress("ErrZ",&fErrZ);
   Vertices->SetBranchAddress("Err2D",&fErr2D);
   Vertices->SetBranchAddress("Err3D",&fErr3D);
	Vertices->SetBranchAddress("Mult",&fMultTrack);
   Vertices->SetBranchAddress("Chi2PerDOF",&fChiPerDOF);
   Vertices->SetBranchAddress("dcaTrack.",&dcaTrack); 
   Vertices->SetBranchAddress("pTTrack.",&pTTrack) ;
	Vertices->SetBranchAddress("primaryVertex.",&primVtx) ;
	
//     This is the loop skeleton
//       To read only selected branches, Insert statements like:
// Vertices->SetBranchStatus("*",0);  // disable all branches
// TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

			//EventQuality-plot
		Int_t nXBinsEvtQ			= 9;
		Double_t firstXBinEvtQ	= -1.5;
		Double_t lastXBinEvtQ	= 7.5;

	//______________Hadronic Interaction _________________
		Int_t nMappingZHadInt 			= 3000;
		Double_t firstBinMappingZHadInt 	= -300.;
		Double_t lastBinMappingZHadInt 	= 300.;

		Int_t nMappingRHadInt 			= 1800;
		Double_t firstBinMappingRHadInt 	= 0.;
		Double_t lastBinMappingRHadInt 	= 180.;

		Int_t nMappingXHadInt 			= 2400;
		Double_t firstBinMappingXHadInt 	= -120.;
		Double_t lastBinMappingXHadInt 	= 120.;

		Int_t nMappingYHadInt 			= 2400;
		Double_t firstBinMappingYHadInt 	= -120.;
		Double_t lastBinMappingYHadInt 	= 120.;

		Int_t nMappingPhiHadInt 		= 400;
		Double_t firstBinMappingPhiHadInt= -TMath::Pi();
		Double_t lastBinMappingPhiHadInt = TMath::Pi();

		Int_t nMappingEtaHadInt 			= 400;
		Double_t firstBinMappingEtaHadInt 	= -10.;
		Double_t lastBinMappingEtaHadInt 	= 10.;

		Int_t nMappingXHadIntHotZone 			= 1000;
		Double_t firstBinMappingXHadIntHotZone 	= -50.;
		Double_t lastBinMappingXHadIntHotZone 	= 50.;

		Int_t nMappingYHadIntHotZone 			= 1000;
		Double_t firstBinMappingYHadIntHotZone 	= -50.;
		Double_t lastBinMappingYHadIntHotZone 	= 50.;
 
		Int_t nMappingZHadIntHotZone 			= 300;
		Double_t firstBinMappingZHadIntHotZone 	= 45.;
		Double_t lastBinMappingZHadIntHotZone 	= 75.;
 
		Int_t nMappingXBPIntHad 			= 600;
		Double_t firstBinMappingXHadIntBP	= -15.;
		Double_t lastBinMappingXHadIntBP 	= 15.;
 
		Int_t nMappingYBPIntHad 			= 600;
		Double_t firstBinMappingYHadIntBP 	= -15.;
		Double_t lastBinMappingYHadIntBP 	= 15.;

		Int_t nErr1DHadInt 			= 1000;
		Double_t firstBinErr1DHadInt 	= 0.;
		Double_t lastBinErr1DHadInt 	= 10.;

		Int_t nErr2DHadInt 			= 1000;
		Double_t firstBinErr2DHadInt 	= 0.;
		Double_t lastBinErr2DHadInt 	= 10.;

		Int_t nErr3DHadInt 			= 1000;
		Double_t firstBinErr3DHadInt 	= 0.;
		Double_t lastBinErr3DHadInt 	= 10.;

		Int_t nChi2HadInt 			= 200;
		Double_t firstBinChi2HadInt 	= 0;
		Double_t lastBinChi2HadInt	 	= 50.;

		Int_t nPtHadInt 			= 200;
		Double_t firstBinPtHadInt 	= 0;
		Double_t lastBinPtHadInt	 	= 20.;

		Int_t nESDSectrk		= 15;
		Double_t firstBinESDSectrk= -0.5;
		Double_t lastBinESDSectrk = 14.5;

	TH1F* histoEventQual ;
	TH2F* histoZR ;
	TH2F* histoXY;
 	TH2F* histoZPhi ;
	TH2F* histoRPhi ;
 	TH2F* histoZPhi_SPD;
	TH2F* histoZPhi_SPDTh;
	TH2F* histoZPhi_SSD;
	TH2F* histoZPhi_SDD ;
	TH2F* histoZPhi_ITSTPC;
	TH2F* histoZPhi_HotZone;
	TH2F* histoXY_HotZone ;
 	TH2F* histoXY_Beampipe;
	TH1F* histoEta;
	TH1F* histoErrX;
	TH1F* histoErrY ;
	TH1F* histoErrZ;
	TH1F* histoErr2D;
	TH1F* histoErr3D ;
	TH1F* histoChi2PerDOF;
	TH1F* histoDCA;
	TH1F* histoDCAXY ;
	TH1F* histoPt ;
	TH1F* histoNTrackSecVtx;
	TH1F* histoErrX_AftCuts;
	TH1F* histoErrY_AftCuts;
	TH1F* histoErrZ_AftCuts;
	TH1F* histoErr2D_AftCuts;
	TH1F* histoErr3D_AftCuts ;
	TH1F* histoChi2PerDOF_AftCuts ;
	TH1F* histoNTrackSecVtx_AftCuts;
	TH1F* histoDCAXY_AftCuts ;
	TH1F* histoDCA_AftCuts ;
	TH1F* histoPt_AftCuts ;
	
		
		
	TString fileNameOutput;
	if (kMC){
		fileNameOutput = Form("%s_MC.root",nameOutput.Data());
	} else {
		fileNameOutput = Form("%s_Data.root",nameOutput.Data());
	}
	TFile* fileMappingDetailedSecHadInt = new TFile(fileNameOutput.Data());
	TDirectory*  directorySecHad = 		(TDirectory*)fileMappingDetailedSecHadInt->Get(Form("SecHadInt_%s",  CutSelection.Data())); 	
	if (!merge || directorySecHad==0) {			
		histoEventQual = new TH1F("ESD_EventQuality","ESD_EventQuality",nXBinsEvtQ,firstXBinEvtQ,lastXBinEvtQ);
		
		histoZR =  new TH2F("ESD_HadIntMap_ZR","", nMappingZHadInt, firstBinMappingZHadInt, lastBinMappingZHadInt, nMappingRHadInt, firstBinMappingRHadInt, lastBinMappingRHadInt);
		histoXY =  new TH2F("ESD_HadIntMap_XY","", nMappingXHadInt, firstBinMappingXHadInt, lastBinMappingXHadInt, nMappingYHadInt, firstBinMappingYHadInt, lastBinMappingYHadInt);
 		histoZPhi = new 	TH2F("ESD_HadIntMap_ZPhi","", nMappingZHadInt, firstBinMappingZHadInt, lastBinMappingZHadInt, nMappingPhiHadInt, firstBinMappingPhiHadInt, lastBinMappingPhiHadInt);
		histoRPhi = new 	TH2F("ESD_HadIntMap_RPhi","",nMappingRHadInt, firstBinMappingRHadInt, lastBinMappingRHadInt, nMappingPhiHadInt, firstBinMappingPhiHadInt, lastBinMappingPhiHadInt);
 		histoZPhi_SPD = new 	TH2F(	"ESD_HadIntMap_SPD_ZPhi","",nMappingZHadInt, firstBinMappingZHadInt, lastBinMappingZHadInt, nMappingPhiHadInt, firstBinMappingPhiHadInt, lastBinMappingPhiHadInt);
		histoZPhi_SPDTh = new 	TH2F(	"ESD_HadIntMap_SPDTh_ZPhi","",nMappingZHadInt, firstBinMappingZHadInt, lastBinMappingZHadInt, nMappingPhiHadInt, firstBinMappingPhiHadInt, lastBinMappingPhiHadInt);	
		histoZPhi_SSD = new 	TH2F("ESD_HadIntMap_SSD_ZPhi","",nMappingZHadInt, firstBinMappingZHadInt, lastBinMappingZHadInt, nMappingPhiHadInt, firstBinMappingPhiHadInt, lastBinMappingPhiHadInt);
		histoZPhi_SDD = new 	TH2F("ESD_HadIntMap_SDD_ZPhi","",nMappingZHadInt, firstBinMappingZHadInt, lastBinMappingZHadInt, nMappingPhiHadInt, firstBinMappingPhiHadInt, lastBinMappingPhiHadInt);	
		histoZPhi_ITSTPC = new 	TH2F("ESD_HadIntMap_ITSTPC_ZPhi","",nMappingZHadInt, firstBinMappingZHadInt, lastBinMappingZHadInt, nMappingPhiHadInt, firstBinMappingPhiHadInt, lastBinMappingPhiHadInt);
		histoZPhi_HotZone = new 	TH2F("ESD_HadIntMap_HotZone_ZPhi","",nMappingZHadIntHotZone, firstBinMappingZHadIntHotZone, lastBinMappingZHadIntHotZone, nMappingPhiHadInt, firstBinMappingPhiHadInt, lastBinMappingPhiHadInt);
		histoXY_HotZone = new 	TH2F("ESD_HadIntMap_HotZone_XY" ,"" , nMappingXHadIntHotZone, firstBinMappingXHadIntHotZone, lastBinMappingXHadIntHotZone, nMappingYHadIntHotZone, firstBinMappingYHadIntHotZone, lastBinMappingYHadIntHotZone);
 		histoXY_Beampipe = new 	TH2F("ESD_HadIntMapInnerBeampipe_XY" ,"" , nMappingXBPIntHad, firstBinMappingXHadIntBP, lastBinMappingXHadIntBP, nMappingYBPIntHad, firstBinMappingYHadIntBP, lastBinMappingYHadIntBP);	
		histoEta = new TH1F( "ESD_HadIntMap_Eta","", nMappingEtaHadInt, firstBinMappingEtaHadInt, lastBinMappingEtaHadInt);	

		histoErrX = new TH1F( "ESD_HadIntQual_ErrX","", nErr1DHadInt, firstBinErr1DHadInt, lastBinErr1DHadInt);
		histoErrY = new TH1F( "ESD_HadIntQual_ErrY","", nErr1DHadInt, firstBinErr1DHadInt, lastBinErr1DHadInt);
		histoErrZ = new TH1F( "ESD_HadIntQual_ErrZ","",  nErr1DHadInt, firstBinErr1DHadInt, lastBinErr1DHadInt);
		histoErr2D = new TH1F( "ESD_HadIntQual_Err2D","", nErr2DHadInt, firstBinErr2DHadInt, lastBinErr2DHadInt);
		histoErr3D = new TH1F( "ESD_HadIntQual_Err3D","", nErr3DHadInt, firstBinErr3DHadInt, lastBinErr3DHadInt);
		histoChi2PerDOF = new TH1F( "ESD_HadIntQual_Chi2PerDOF","", nChi2HadInt, firstBinChi2HadInt, lastBinChi2HadInt);
		histoDCA = new TH1F( "ESD_HadIntQual_DCA","", nErr1DHadInt, firstBinErr1DHadInt, lastBinErr1DHadInt);
		histoDCAXY = new TH1F( "ESD_HadIntQual_DCAXY","", nErr1DHadInt, firstBinErr1DHadInt, lastBinErr1DHadInt);
		histoPt = new TH1F( "ESD_HadIntQual_Pt","", nPtHadInt, firstBinPtHadInt, lastBinPtHadInt);
		histoNTrackSecVtx = new TH1F("ESD_HadIntQual_nTracksSecVtx","", nESDSectrk, firstBinESDSectrk, lastBinESDSectrk);
// 		
		histoErrX_AftCuts = new TH1F( "ESD_HadIntQualAftCuts_ErrX","", nErr1DHadInt, firstBinErr1DHadInt, lastBinErr1DHadInt);
		histoErrY_AftCuts = new TH1F( "ESD_HadIntQualAftCuts_ErrY","", nErr1DHadInt, firstBinErr1DHadInt, lastBinErr1DHadInt);
		histoErrZ_AftCuts = new TH1F( "ESD_HadIntQualAftCuts_ErrZ","",  nErr1DHadInt, firstBinErr1DHadInt, lastBinErr1DHadInt);
		histoErr2D_AftCuts = new TH1F( "ESD_HadIntQualAftCuts_Err2D","", nErr2DHadInt, firstBinErr2DHadInt, lastBinErr2DHadInt);
		histoErr3D_AftCuts = new TH1F( "ESD_HadIntQualAftCuts_Err3D","", nErr3DHadInt, firstBinErr3DHadInt, lastBinErr3DHadInt);
		histoChi2PerDOF_AftCuts = new TH1F( "ESD_HadIntQualAftCuts_Chi2PerDOF","", nChi2HadInt, firstBinChi2HadInt, lastBinChi2HadInt);
		histoNTrackSecVtx_AftCuts = new TH1F("ESD_HadIntQualAftCuts_nTracksSecVtx","", nESDSectrk, firstBinESDSectrk, lastBinESDSectrk);
		histoDCAXY_AftCuts = new TH1F( "ESD_HadIntQualAftCuts_DCAXY","", nErr1DHadInt, firstBinErr1DHadInt, lastBinErr1DHadInt);
		histoDCA_AftCuts = new TH1F( "ESD_HadIntQualAftCuts_DCA","", nErr1DHadInt, firstBinErr1DHadInt, lastBinErr1DHadInt);
		histoPt_AftCuts = new TH1F( "ESD_HadIntQualAftCuts_Pt","", nPtHadInt, firstBinPtHadInt, lastBinPtHadInt);
	} else {
		histoEventQual = 	(TH1F*)directorySecHad->Get("ESD_EventQuality");
		histoNTrackSecVtx = 	(TH1F*)directorySecHad->Get("ESD_HadIntQual_nTracksSecVtx");
		histoNTrackSecVtx_AftCuts = 	(TH1F*)directorySecHad->Get("ESD_HadIntQualAftCuts_nTracksSecVtx");
	
		histoRPhi = 	(TH2F*)directorySecHad->Get("ESD_HadIntMap_RPhi");
		histoXY = 	(TH2F*)directorySecHad->Get("ESD_HadIntMap_XY");
		histoZR = 	(TH2F*)directorySecHad->Get("ESD_HadIntMap_ZR");
		histoZPhi = 	(TH2F*)directorySecHad->Get("ESD_HadIntMap_ZPhi");
		histoZPhi_SPD = 	(TH2F*)directorySecHad->Get("ESD_HadIntMap_SPD_ZPhi"); 
		histoZPhi_SPDTh = 	(TH2F*)directorySecHad->Get("ESD_HadIntMap_SPDTh_ZPhi"); 
		histoZPhi_SSD = 	(TH2F*)directorySecHad->Get("ESD_HadIntMap_SSD_ZPhi"); 
		histoZPhi_SDD = 	(TH2F*)directorySecHad->Get("ESD_HadIntMap_SDD_ZPhi"); 
		histoZPhi_ITSTPC = 	(TH2F*)directorySecHad->Get("ESD_HadIntMap_ITSTPC_ZPhi"); 
		histoZPhi_HotZone = 	(TH2F*)directorySecHad->Get("ESD_HadIntMap_HotZone_ZPhi"); 
		histoXY_HotZone = 	(TH2F*)directorySecHad->Get("ESD_HadIntMap_HotZone_XY"); 
		histoXY_Beampipe = 	(TH2F*)directorySecHad->Get("ESD_HadIntMapInnerBeampipe_XY"); 
		histoEta = 	(TH1F*)directorySecHad->Get("ESD_HadIntMap_Eta");

		histoChi2PerDOF = 	(TH1F*)directorySecHad->Get("ESD_HadIntQual_Chi2PerDOF");
		histoErr2D = 	(TH1F*)directorySecHad->Get("ESD_HadIntQual_Err2D");
		histoErr3D = 	(TH1F*)directorySecHad->Get("ESD_HadIntQual_Err3D");
		histoErrX = 	(TH1F*)directorySecHad->Get("ESD_HadIntQual_ErrX");
		histoErrY = 	(TH1F*)directorySecHad->Get("ESD_HadIntQual_ErrY");
		histoErrZ = 	(TH1F*)directorySecHad->Get("ESD_HadIntQual_ErrZ");
		histoDCA = 	(TH1F*)directorySecHad->Get("ESD_HadIntQual_DCA");
		histoDCAXY = 	(TH1F*)directorySecHad->Get("ESD_HadIntQual_DCAXY");
		histoPt = 	(TH1F*)directorySecHad->Get("ESD_HadIntQual_Pt");
	
		histoChi2PerDOF_AftCuts = 	(TH1F*)directorySecHad->Get("ESD_HadIntQualAftCuts_Chi2PerDOF");
		histoErr2D_AftCuts = 	(TH1F*)directorySecHad->Get("ESD_HadIntQualAftCuts_Err2D");
		histoErr3D_AftCuts = 	(TH1F*)directorySecHad->Get("ESD_HadIntQualAftCuts_Err3D");
		histoErrX_AftCuts = 	(TH1F*)directorySecHad->Get("ESD_HadIntQualAftCuts_ErrX");
		histoErrY_AftCuts = 	(TH1F*)directorySecHad->Get("ESD_HadIntQualAftCuts_ErrY");
		histoErrZ_AftCuts = 	(TH1F*)directorySecHad->Get("ESD_HadIntQualAftCuts_ErrZ");
		histoDCA_AftCuts = 	(TH1F*)directorySecHad->Get("ESD_HadIntQualAftCuts_DCA");
		histoDCAXY_AftCuts = 	(TH1F*)directorySecHad->Get("ESD_HadIntQualAftCuts_DCAXY");
		histoPt_AftCuts = 	(TH1F*)directorySecHad->Get("ESD_HadIntQualAftCuts_Pt");


		histoEventQual->Sumw2();
		histoNTrackSecVtx->Sumw2();
		histoNTrackSecVtx_AftCuts->Sumw2();
		histoRPhi->Sumw2();
		histoXY->Sumw2();
		histoZR->Sumw2();
		histoZPhi->Sumw2();
		histoZPhi_SPD->Sumw2();
		histoZPhi_SPDTh->Sumw2();
		histoZPhi_SSD->Sumw2();
		histoZPhi_SDD->Sumw2();
		histoZPhi_ITSTPC->Sumw2();
		histoZPhi_HotZone->Sumw2();
		histoXY_HotZone->Sumw2();
		histoXY_Beampipe->Sumw2();
		histoEta->Sumw2();
		histoChi2PerDOF->Sumw2();
		histoErr2D->Sumw2();
		histoErr3D->Sumw2();
		histoErrX->Sumw2();
		histoErrY->Sumw2();
		histoErrZ->Sumw2();
		histoDCA->Sumw2();
		histoDCAXY->Sumw2();
		histoPt->Sumw2();
		histoChi2PerDOF_AftCuts->Sumw2();
		histoErr2D_AftCuts->Sumw2();
		histoErr3D_AftCuts->Sumw2();
		histoErrX_AftCuts->Sumw2();
		histoErrY_AftCuts->Sumw2();
		histoErrZ_AftCuts->Sumw2();
		histoDCA_AftCuts->Sumw2();
		histoDCAXY_AftCuts->Sumw2();
		histoPt_AftCuts->Sumw2();
		
	}
   ULong64_t nEntriesVertex = Vertices->GetEntries();
   ULong64_t nbytesV = 0;
	ULong64_t nEntriesAccepted = 0;
	ULong64_t nEntriesAftCuts = 0;
	for (ULong64_t i=0; i<nEntriesVertex;i++) {
		nbytesV += Vertices->GetEntry(i);
		histoNTrackSecVtx->Fill(nTracksPerVtx);
		if ( nTracksPerVtx >= minSecHadTrackCut){
			histoErrX->Fill(fErrX);
			histoErrY->Fill(fErrY);
			histoErrZ->Fill(fErrZ);
			histoErr2D->Fill(fErr2D);
			histoErr3D->Fill(fErr3D);
			histoChi2PerDOF->Fill(fChiPerDOF);	
			
			Double_t dca[3];
			dca[0] = (*dcaTrack)[0];
			dca[1] = (*dcaTrack)[1];
			dca[2] = (*dcaTrack)[2];
			Int_t sort = nTracksPerVtx;
			while (sort != 1) {
				for (Int_t ii = 0  ; ii < sort; ii++){
					if ((ii+1) != nTracksPerVtx){
						if ( dca[ii+1] < dca[ii]) {
							Int_t intermed = dca[ii+1];
							dca[ii+1] = dca[ii];
							dca[ii] = intermed;
						}
					}
				}
				sort--;
			}
			if (TMath::Abs(fEta)>=etaMinCut && TMath::Abs(fEta)<etaMaxCut ){
				nEntriesAccepted++;
				histoErrX->Fill(fErrX);
				histoErrY->Fill(fErrY);
				histoErrZ->Fill(fErrZ);
				histoErr2D->Fill(fErr2D);
				histoErr3D->Fill(fErr3D);
				histoChi2PerDOF->Fill(fChiPerDOF);
				histoDCA->Fill(dca[1]);
				for (Int_t j = 0; j < nTracksPerVtx && j < 3; j++){
					histoPt->Fill((*pTTrack)[j]);
				}
				Double_t maximumDCAXY= dca[1]/(TMath::Sqrt(TMath::Power(fX-(*primVtx)[0],2) + TMath::Power(fY-(*primVtx)[1],2)))*fMultTrack;
				histoDCAXY->Fill(maximumDCAXY);
				if ( (*pTTrack)[1] > minPtTrackCut && maximumDCAXY < dcaTrackCut && dca[1] < dcaTrackAbsCut ){
					if (fChiPerDOF < chi2PerDOFCut && fErr2D < err2DCut){ //
						nEntriesAftCuts++;
						histoZR->Fill(fZ,fR);
						histoXY->Fill(fX,fY);
						histoZPhi->Fill(fZ,fPhi);
						histoRPhi->Fill(fR,fPhi);
						histoEta->Fill(fEta);
						if ( fR < rSPD){
							histoZPhi_SPD->Fill(fZ,fPhi);
						} else if (fR < rSPDTh){
							histoZPhi_SPDTh->Fill(fZ,fPhi);
						} else if (fR < rSDD){
							histoZPhi_SDD->Fill(fZ,fPhi);
						} else if (fR < rSSD){
							histoZPhi_SSD->Fill(fZ,fPhi);
						} else if (fR < rTPC){
							histoZPhi_ITSTPC->Fill(fZ,fPhi);
						}
						if (fR>rHotZoneMin && fR < rHotZoneMax){
							histoZPhi_HotZone->Fill( fZ,fPhi);
							if (fZ < zHotZoneMax && fZ > zHotZoneMin){
								histoXY_HotZone->Fill( fX,fY);
							}
						}	
						histoErrX_AftCuts->Fill(fErrX);
						histoErrY_AftCuts->Fill(fErrY);
						histoErrZ_AftCuts->Fill(fErrZ);
						histoErr2D_AftCuts->Fill(fErr2D);
						histoErr3D_AftCuts->Fill(fErr3D);
						histoChi2PerDOF_AftCuts->Fill(fChiPerDOF);			
						histoDCAXY_AftCuts->Fill(maximumDCAXY);
						histoNTrackSecVtx_AftCuts->Fill(nTracksPerVtx);
						histoDCA_AftCuts->Fill(dca[1]);
						for (Int_t j = 0; j < nTracksPerVtx && j < 3; j++){
							histoPt_AftCuts->Fill((*pTTrack)[j]);
						}
					}
				}
			}
			if (fChiPerDOF < chi2PerDOFCut && fErr2D < err2DCut && fZ < zBeamPipeInner && fZ > - zBeamPipeInner){
				histoXY_Beampipe->Fill( fX,fY);
			}		
		}
	}

	cout << "Vertices reconstructed:   " << nEntriesVertex << "   in eta range:  " << nEntriesAccepted << "  after cuts:   "  << nEntriesAftCuts << endl;			
	ULong64_t nbytesE = 0;
	for (ULong64_t i=0; i<nEntriesEvent;i++) {
     nbytesE += Event->GetEntry(i);
	  histoEventQual->Fill(5);
	}

	TFile* fileMappingWrite = new TFile(fileNameOutput.Data(),"UPDATE");
	fileMappingWrite->mkdir(Form("SecHadInt_%s",  CutSelection.Data()));
	fileMappingWrite->cd(Form("SecHadInt_%s",  CutSelection.Data()));
		histoEventQual->Write("ESD_EventQuality",TObject::kWriteDelete);
		histoNTrackSecVtx->Write("ESD_HadIntQual_nTracksSecVtx",TObject::kWriteDelete);
		histoNTrackSecVtx_AftCuts->Write("ESD_HadIntQualAftCuts_nTracksSecVtx",TObject::kWriteDelete);
		
		histoRPhi->Write("ESD_HadIntMap_RPhi",TObject::kWriteDelete);
		histoXY->Write("ESD_HadIntMap_XY",TObject::kWriteDelete);
		histoZR->Write("ESD_HadIntMap_ZR",TObject::kWriteDelete);
		histoZPhi->Write("ESD_HadIntMap_ZPhi",TObject::kWriteDelete);
		histoZPhi_SPD->Write("ESD_HadIntMap_SPD_ZPhi",TObject::kWriteDelete); 
		histoZPhi_SPDTh->Write("ESD_HadIntMap_SPDTh_ZPhi",TObject::kWriteDelete); 
		histoZPhi_SSD->Write("ESD_HadIntMap_SSD_ZPhi",TObject::kWriteDelete); 
		histoZPhi_SDD->Write("ESD_HadIntMap_SDD_ZPhi",TObject::kWriteDelete); 
		histoZPhi_ITSTPC->Write("ESD_HadIntMap_ITSTPC_ZPhi",TObject::kWriteDelete); 
		histoZPhi_HotZone->Write("ESD_HadIntMap_HotZone_ZPhi",TObject::kWriteDelete); 
		histoXY_HotZone->Write("ESD_HadIntMap_HotZone_XY",TObject::kWriteDelete); 
		histoXY_Beampipe->Write("ESD_HadIntMapInnerBeampipe_XY",TObject::kWriteDelete); 
		histoEta->Write("ESD_HadIntMap_Eta",TObject::kWriteDelete);

		histoChi2PerDOF->Write("ESD_HadIntQual_Chi2PerDOF",TObject::kWriteDelete);
		histoErr2D->Write("ESD_HadIntQual_Err2D",TObject::kWriteDelete);
		histoErr3D->Write("ESD_HadIntQual_Err3D",TObject::kWriteDelete);
		histoErrX->Write("ESD_HadIntQual_ErrX",TObject::kWriteDelete);
		histoErrY->Write("ESD_HadIntQual_ErrY",TObject::kWriteDelete);
		histoErrZ->Write("ESD_HadIntQual_ErrZ",TObject::kWriteDelete);
		histoDCA->Write("ESD_HadIntQual_DCA",TObject::kWriteDelete);
		histoDCAXY->Write("ESD_HadIntQual_DCAXY",TObject::kWriteDelete);
		histoPt->Write("ESD_HadIntQual_Pt",TObject::kWriteDelete);
		
		histoChi2PerDOF_AftCuts->Write("ESD_HadIntQualAftCuts_Chi2PerDOF",TObject::kWriteDelete);
		histoErr2D_AftCuts->Write("ESD_HadIntQualAftCuts_Err2D",TObject::kWriteDelete);
		histoErr3D_AftCuts->Write("ESD_HadIntQualAftCuts_Err3D",TObject::kWriteDelete);
		histoErrX_AftCuts->Write("ESD_HadIntQualAftCuts_ErrX",TObject::kWriteDelete);
		histoErrY_AftCuts->Write("ESD_HadIntQualAftCuts_ErrY",TObject::kWriteDelete);
		histoErrZ_AftCuts->Write("ESD_HadIntQualAftCuts_ErrZ",TObject::kWriteDelete);
		histoDCA_AftCuts->Write("ESD_HadIntQualAftCuts_DCA",TObject::kWriteDelete);
		histoDCAXY_AftCuts->Write("ESD_HadIntQualAftCuts_DCAXY",TObject::kWriteDelete);
		histoPt_AftCuts->Write("ESD_HadIntQualAftCuts_Pt",TObject::kWriteDelete);
		
	fileMappingWrite->Write();
	fileMappingWrite->Close();

}