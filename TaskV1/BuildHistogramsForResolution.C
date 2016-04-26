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
#include "TArc.h"

typedef TVectorT<double> TVectorD;

void BuildHistogramsForResolution(TString fileName = "GammaConvV1_Resolution_0000010020092663003800000.root", TString CutSelection="0199" , Bool_t kMC=0, Bool_t merge=0, TString nameOutputBase = "ResolutionGammaConv"){
	//Reset ROOT and connect tree file
	
//********************************************************************************
//* 				Definition of Cuts																*
//********************************************************************************
	
	TString V0Finder = CutSelection(0,1);
	TString etaCutNumber = CutSelection(1,1);
	TString minPtCutNumber = CutSelection(2,1);
	TString chi2CutNumber = CutSelection(3,1);
	
	if ( V0Finder.CompareTo("0") == 0){
		cout << "V0Finder: " <<  "Onfly" << endl;
	} else if (	V0Finder.CompareTo("1") == 0){
		cout << "V0Finder: " <<  "Offline" << endl;
	}
	
	Double_t etaMinCut = -0.1;
	Double_t etaMaxCut = 0.;
	if (etaCutNumber.CompareTo("0") == 0){
		etaMaxCut = 0.9;
 	} else if (etaCutNumber.CompareTo("1") == 0 ){
		etaMaxCut = 0.1;
 	} else if (etaCutNumber.CompareTo("2") == 0 ){
		etaMaxCut = 1.4;
		etaMinCut = 0.9;
 	}else if (etaCutNumber.CompareTo("3") == 0 ){
		etaMaxCut = 1.8;
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
 	} else {
		etaMaxCut = 10.;
		etaMinCut = -0.1;
 	}
 	if (etaMinCut != -0.1){
		cout << "Eta range: " <<   etaMinCut  <<  " < |eta| < " << etaMaxCut << endl;
	} else {
		cout << "Eta range: " <<  "|eta| < " << etaMaxCut << endl;
	}
	
	Double_t minPtPhotonCut = 0.0;
	if (minPtCutNumber.CompareTo("0") == 0){
		minPtPhotonCut = 0.05;
 	} else if (minPtCutNumber.CompareTo("1") == 0 ){
		minPtPhotonCut = 0.1;
 	} else if (minPtCutNumber.CompareTo("2") == 0 ){
		minPtPhotonCut = 0.15;
	} else if (minPtCutNumber.CompareTo("3") == 0 ){
		minPtPhotonCut = 0.2;
 	} else if (minPtCutNumber.CompareTo("4") == 0 ){
		minPtPhotonCut = 0.3;
	}
	cout<< "Min pt for tracks: " << minPtPhotonCut << endl;
	
	Double_t chi2PerDOFCut = 100000.0;
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
	
   gROOT->Reset();

//********************************************************************************
//* 				File definition/ loading 														*
//********************************************************************************
	
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fileName.Data());
   if (!f) {
      f = new TFile(fileName.Data());
   }


//********************************************************************************
//* 				Definition of Reconstructed Conversion Points							*
//********************************************************************************
	
	   // Declaration of leaf types
   Float_t         ESDpt;
   Float_t         ESDphi;
   Float_t         ESDeta;
   Float_t         ESDr;
   Float_t         ESDz;
   Float_t         MCpt;
   Float_t         MCphi;
   Float_t         MCeta;
   Float_t         MCr;
   Float_t         MCz;
   Float_t         chi2ndf;

   
	TTree *treeResolution = (TTree*)gDirectory->Get("Resolution");

   // Set branch addresses.
   treeResolution->SetBranchAddress("ESDpt",&ESDpt);
   treeResolution->SetBranchAddress("ESDphi",&ESDphi);
   treeResolution->SetBranchAddress("ESDeta",&ESDeta);
   treeResolution->SetBranchAddress("ESDr",&ESDr);
   treeResolution->SetBranchAddress("ESDz",&ESDz);
   treeResolution->SetBranchAddress("MCpt",&MCpt);
   treeResolution->SetBranchAddress("MCphi",&MCphi);
   treeResolution->SetBranchAddress("MCeta",&MCeta);
   treeResolution->SetBranchAddress("MCr",&MCr);
	treeResolution->SetBranchAddress("MCz",&MCz);
	treeResolution->SetBranchAddress("chi2ndf",&chi2ndf);
	
	
//********************************************************************************
//* 				Definition of Boundaries for Histograms									*
//********************************************************************************

	//ResolutionPlots
	//RESdPt
	Int_t knBinsResdPt		= 200;
	Int_t kfirstBinResdPt	= 0;
	Int_t klastBinResdPt		= 25;
	Int_t knYBinsResdPt		= 200;
	Int_t kfirstYBinResdPt	= -10;
	Int_t klastYBinResdPt		= 10;

	//RESdR
	Int_t knBinsResdR			= 720;
	Int_t kfirstBinResdR		= 0;
	Int_t klastBinResdR		= 180;
	Int_t knBinsResdZ			= 1600;
	Int_t kfirstBinResdZ		= -200;
	Int_t klastBinResdZ		= 200;
	Int_t knYBinsResdR			= 100;
	Int_t kfirstYBinResdR		= -25;
	Int_t klastYBinResdR		= 25;

	//RESdZ

	Int_t knBinsEta 			= 400;
	Double_t kfirstBinEta 	= -2.;
	Double_t klastBinEta 	= 2.;

	
//********************************************************************************
//* 		Definition of histograms for reconstructed Conversion Points 				*
//********************************************************************************	
	TH2F* histodAbsRVsR; 
	TH2F* histodAbsRVsZ; 
	TH2F* histodAbsRVsPhi;
	TH2F* histodAbsRVsPt; 
	TH2F* histodAbsRVsEta; 

	TH2F* histodAbsZVsR; 
	TH2F* histodAbsZVsZ; 
	TH2F* histodAbsZVsPhi;
	TH2F* histodAbsZVsPt; 
	TH2F* histodAbsZVsEta;
	
	TH2F* histodAbsPhiVsR; 
	TH2F* histodAbsPhiVsZ; 
	TH2F* histodAbsPhiVsPhi;
	TH2F* histodAbsPhiVsPt; 
	TH2F* histodAbsPhiVsEta; 

	TH2F* histodAbsEtaVsR; 
	TH2F* histodAbsEtaVsZ; 
	TH2F* histodAbsEtaVsPhi;
	TH2F* histodAbsEtaVsPt; 
	TH2F* histodAbsEtaVsEta; 

	TH2F* histodPtVsR; 
	TH2F* histodPtVsZ; 
	TH2F* histodPtVsPhi;
	TH2F* histodPtVsPt; 
	TH2F* histodPtVsEta; 

	TString fileNameOutput = Form("%s.root",nameOutputBase.Data()) ;
	
	TFile* fileMappingDetailedConv = new TFile(fileNameOutput.Data());
	TDirectory*  directoryConv = 		(TDirectory*)fileMappingDetailedConv->Get(Form("GammaConv_%s",  CutSelection.Data())); 			
	if (!merge || directoryConv==0) {	
		
		histodAbsRVsR = new 	TH2F("Resolution_dRAbs_VS_R","" ,knBinsResdR, kfirstBinResdR, klastBinResdR,knYBinsResdR,kfirstYBinResdR, klastYBinResdR);
		histodAbsRVsZ= new 	TH2F("Resolution_dRAbs_VS_Z","" ,knBinsResdZ, kfirstBinResdZ, klastBinResdZ, knYBinsResdR, kfirstYBinResdR, klastYBinResdR);
		histodAbsRVsPhi= new 	TH2F("Resolution_dRAbs_VS_Phi","" ,knBinsResdPt, 0, 2*TMath::Pi(), knYBinsResdR, kfirstYBinResdR, klastYBinResdR);
		histodAbsRVsPt= new 	TH2F("Resolution_dRAbs_VS_Pt","" ,knYBinsResdPt, kfirstBinResdPt, klastBinResdPt, knYBinsResdR, kfirstYBinResdR, klastYBinResdR);
		histodAbsRVsEta= new 	TH2F("Resolution_dRAbs_VS_Eta","" ,knBinsEta, kfirstBinEta, klastBinEta, knYBinsResdR, kfirstYBinResdR, klastYBinResdR);
		
		histodAbsZVsR= new 	TH2F("Resolution_dZAbs_VS_R","" ,knBinsResdR, kfirstBinResdR, klastBinResdR,knYBinsResdR,kfirstYBinResdR, klastYBinResdR);
		histodAbsZVsZ= new 	TH2F("Resolution_dZAbs_VS_Z","" ,knBinsResdZ, kfirstBinResdZ, klastBinResdZ, knYBinsResdR, kfirstYBinResdR, klastYBinResdR);
		histodAbsZVsPhi= new 	TH2F("Resolution_dZAbs_VS_Phi","" ,knBinsResdPt, 0, 2*TMath::Pi(), knYBinsResdR, kfirstYBinResdR, klastYBinResdR);
		histodAbsZVsPt= new 	TH2F("Resolution_dZAbs_VS_Pt","" ,knYBinsResdPt, kfirstBinResdPt, klastBinResdPt, knYBinsResdR, kfirstYBinResdR, klastYBinResdR);
		histodAbsZVsEta= new 	TH2F("Resolution_dZAbs_VS_Eta","" ,knBinsEta, kfirstBinEta, klastBinEta, knYBinsResdR, kfirstYBinResdR, klastYBinResdR);
		
		histodAbsPhiVsR= new 	TH2F("Resolution_dPhiAbs_VS_R","" ,knBinsResdR, kfirstBinResdR, klastBinResdR,knBinsResdPt, -TMath::Pi()/30., TMath::Pi()/30.);
		histodAbsPhiVsZ= new 	TH2F("Resolution_dPhiAbs_VS_Z","" ,knBinsResdZ, kfirstBinResdZ, klastBinResdZ, knBinsResdPt, -TMath::Pi()/30., TMath::Pi()/30.);
		histodAbsPhiVsPhi= new 	TH2F("Resolution_dPhiAbs_VS_Phi","" ,knBinsResdPt, 0, 2*TMath::Pi(), knBinsResdPt, -TMath::Pi()/30., TMath::Pi()/30.);
		histodAbsPhiVsPt= new 	TH2F("Resolution_dPhiAbs_VS_Pt","" ,knYBinsResdPt, kfirstBinResdPt, klastBinResdPt, knBinsResdPt, -TMath::Pi()/30., TMath::Pi()/30.);
		histodAbsPhiVsEta= new 	TH2F("Resolution_dPhiAbs_VS_Eta","" ,knBinsEta, kfirstBinEta, klastBinEta, knBinsResdPt, -TMath::Pi()/30., TMath::Pi()/30.);

		histodAbsEtaVsR = new 	TH2F("Resolution_dEtaAbs_VS_R","" ,knBinsResdR, kfirstBinResdR, klastBinResdR,knBinsEta, kfirstBinEta, klastBinEta);
		histodAbsEtaVsZ= new 	TH2F("Resolution_dEtaAbs_VS_Z","" ,knBinsResdZ, kfirstBinResdZ, klastBinResdZ ,knBinsEta, kfirstBinEta, klastBinEta);
		histodAbsEtaVsPhi= new 	TH2F("Resolution_dEtaAbs_VS_Phi","" ,knBinsResdPt, 0, 2*TMath::Pi(),knBinsEta, kfirstBinEta, klastBinEta);
		histodAbsEtaVsPt= new 	TH2F("Resolution_dEtaAbs_VS_Pt","" ,knYBinsResdPt, kfirstBinResdPt, klastBinResdPt,knBinsEta, kfirstBinEta, klastBinEta);
		histodAbsEtaVsEta= new 	TH2F("Resolution_dEtaAbs_VS_Eta","" ,knBinsEta, kfirstBinEta, klastBinEta,knBinsEta, kfirstBinEta, klastBinEta);

		histodPtVsR = new 	TH2F("Resolution_dPt_VS_R","" ,knBinsResdR, kfirstBinResdR, klastBinResdR,knYBinsResdPt, kfirstYBinResdPt, klastYBinResdPt);
		histodPtVsZ= new 	TH2F("Resolution_dPt_VS_Z","" ,knBinsResdZ, kfirstBinResdZ, klastBinResdZ ,knYBinsResdPt, kfirstYBinResdPt, klastYBinResdPt);
		histodPtVsPhi= new 	TH2F("Resolution_dPt_VS_Phi","" ,knBinsResdPt, 0, 2*TMath::Pi(),knYBinsResdPt, kfirstYBinResdPt, klastYBinResdPt);
		histodPtVsPt= new 	TH2F("Resolution_dPt_VS_Pt","" ,knYBinsResdPt, kfirstBinResdPt, klastBinResdPt,knYBinsResdPt, kfirstYBinResdPt, klastYBinResdPt);
		histodPtVsEta= new 	TH2F("Resolution_dPt_VS_Eta","" ,knBinsEta, kfirstBinEta, klastBinEta,knYBinsResdPt, kfirstYBinResdPt, klastYBinResdPt);

	} else { 
		histodAbsRVsR = 	(TH2F*)directoryConv->Get("Resolution_dRAbs_VS_R");
		histodAbsRVsZ = 	(TH2F*)directoryConv->Get("Resolution_dRAbs_VS_Z");
		histodAbsRVsPhi = 	(TH2F*)directoryConv->Get("Resolution_dRAbs_VS_Phi");
		histodAbsRVsPt = 	(TH2F*)directoryConv->Get("Resolution_dRAbs_VS_Pt");
		histodAbsRVsEta = 	(TH2F*)directoryConv->Get("Resolution_dRAbs_VS_Eta");
		
		histodAbsZVsR = 	(TH2F*)directoryConv->Get("Resolution_dZAbs_VS_R");
		histodAbsZVsZ = 	(TH2F*)directoryConv->Get("Resolution_dZAbs_VS_Z");
		histodAbsZVsPhi = 	(TH2F*)directoryConv->Get("Resolution_dZAbs_VS_Phi");
		histodAbsZVsPt = 	(TH2F*)directoryConv->Get("Resolution_dZAbs_VS_Pt");
		histodAbsZVsEta = 	(TH2F*)directoryConv->Get("Resolution_dZAbs_VS_Eta");
		
		histodAbsPhiVsR = 	(TH2F*)directoryConv->Get("Resolution_dPhiAbs_VS_R");
		histodAbsPhiVsZ = 	(TH2F*)directoryConv->Get("Resolution_dPhiAbs_VS_Z");
		histodAbsPhiVsPhi = 	(TH2F*)directoryConv->Get("Resolution_dPhiAbs_VS_Phi");
		histodAbsPhiVsPt = 	(TH2F*)directoryConv->Get("Resolution_dPhiAbs_VS_Pt");
		histodAbsPhiVsEta = 	(TH2F*)directoryConv->Get("Resolution_dPhiAbs_VS_Eta");

		histodAbsEtaVsR = 	(TH2F*)directoryConv->Get("Resolution_dEtaAbs_VS_R");
		histodAbsEtaVsZ = 	(TH2F*)directoryConv->Get("Resolution_dEtaAbs_VS_Z");
		histodAbsEtaVsPhi = 	(TH2F*)directoryConv->Get("Resolution_dEtaAbs_VS_Phi");
		histodAbsEtaVsPt = 	(TH2F*)directoryConv->Get("Resolution_dEtaAbs_VS_Pt");
		histodAbsEtaVsEta = 	(TH2F*)directoryConv->Get("Resolution_dEtaAbs_VS_Eta");

		histodPtVsR = 	(TH2F*)directoryConv->Get("Resolution_dPtAbs_VS_R");
		histodPtVsZ = 	(TH2F*)directoryConv->Get("Resolution_dPtAbs_VS_Z");
		histodPtVsPhi = 	(TH2F*)directoryConv->Get("Resolution_dPtAbs_VS_Phi");
		histodPtVsPt = 	(TH2F*)directoryConv->Get("Resolution_dPtAbs_VS_Pt");
		histodPtVsEta = 	(TH2F*)directoryConv->Get("Resolution_dPtAbs_VS_Eta");

		histodAbsRVsR->Sumw2();
		histodAbsRVsZ->Sumw2();
		histodAbsRVsPhi->Sumw2();
		histodAbsRVsPt->Sumw2();
		histodAbsRVsEta->Sumw2();
		histodAbsZVsR->Sumw2();
		histodAbsZVsZ->Sumw2();
		histodAbsZVsPhi->Sumw2();
		histodAbsZVsPt->Sumw2();
		histodAbsZVsEta->Sumw2();
		histodAbsPhiVsR->Sumw2();
		histodAbsPhiVsZ->Sumw2();
		histodAbsPhiVsPhi->Sumw2();
		histodAbsPhiVsPt->Sumw2();
		histodAbsPhiVsEta->Sumw2();
		histodAbsEtaVsR->Sumw2();
		histodAbsEtaVsZ->Sumw2();
		histodAbsEtaVsPhi->Sumw2();
		histodAbsEtaVsPt->Sumw2();
		histodAbsEtaVsEta->Sumw2();
		histodPtVsR->Sumw2();
		histodPtVsZ->Sumw2();
		histodPtVsPhi->Sumw2();
		histodPtVsPt->Sumw2();
		histodPtVsEta->Sumw2();
	}	
	
//********************************************************************************
//* 		Reading of Tree with reconstructed gammas/ filling of histograms 			*
//********************************************************************************
	
   Long64_t nEntriesRecGam = treeResolution->GetEntries();
	cout << "Number of Gammas to be processed: " << nEntriesRecGam << endl;
	
	
   Long64_t nbytesMCGamma = 0;
	for (Long64_t i=0; i<nEntriesRecGam;i++) {
		nbytesMCGamma += treeResolution->GetEntry(i);
		if ( ESDpt > minPtPhotonCut ){ //fChiPerDOF < chi2PerDOFCut && 
			if (TMath::Abs(ESDeta)>=etaMinCut && TMath::Abs(ESDeta)<etaMaxCut ){ //
				Float_t dEtaAbs = ESDeta - MCeta;
				Float_t dRAbs = ESDr - MCr;
				Float_t dZAbs = ESDz - MCz;
				Float_t dPhiAbs = ESDphi - MCphi;
				Float_t dPt;
				if (MCpt != 0){
					dPt= (ESDpt - MCpt)/MCpt;
				} else {
					dPt = 200;
				}
				histodAbsRVsR->Fill(MCr,dRAbs);		
				histodAbsRVsZ->Fill(MCz,dRAbs);		
				histodAbsRVsPhi->Fill(MCphi,dRAbs);		
				histodAbsRVsPt->Fill(MCpt,dRAbs);		
				histodAbsRVsEta->Fill(MCeta,dRAbs);		
				
				histodAbsZVsR->Fill(MCr,dZAbs);		
				histodAbsZVsZ->Fill(MCz,dZAbs);		
				histodAbsZVsPhi->Fill(MCphi,dZAbs);
				histodAbsZVsPt->Fill(MCpt,dZAbs);
				histodAbsZVsEta->Fill(MCeta,dZAbs);
				
				histodAbsPhiVsR->Fill(MCr,dPhiAbs);
				histodAbsPhiVsZ->Fill(MCz,dPhiAbs);
				histodAbsPhiVsPhi->Fill(MCphi,dPhiAbs);
				histodAbsPhiVsPt->Fill(MCpt,dPhiAbs);
				histodAbsPhiVsEta->Fill(MCeta,dPhiAbs);

				histodAbsEtaVsR->Fill(MCr,dEtaAbs);
				histodAbsEtaVsZ->Fill(MCz,dEtaAbs);
				histodAbsEtaVsPhi->Fill(MCphi,dEtaAbs);
				histodAbsEtaVsPt->Fill(MCpt,dEtaAbs);
				histodAbsEtaVsEta->Fill(MCeta,dEtaAbs);

				histodPtVsR->Fill(MCr,dPt);
				histodPtVsZ->Fill(MCz,dPt);
				histodPtVsPhi->Fill(MCphi,dPt);
				histodPtVsPt->Fill(MCpt,dPt);
				histodPtVsEta->Fill(MCeta,dPt);
			}
		}
	}

	
//********************************************************************************
//* 						Writing histograms to outputfile 			 						*
//********************************************************************************		
	TFile* fileMappingWrite = new TFile(fileNameOutput.Data(),"UPDATE");
	fileMappingWrite->mkdir(Form("GammaConv_%s",  CutSelection.Data()));
	fileMappingWrite->cd(Form("GammaConv_%s",  CutSelection.Data()));
		histodAbsRVsR->Write("Resolution_dRAbs_VS_R",TObject::kWriteDelete);		
		histodAbsRVsZ->Write("Resolution_dRAbs_VS_Z",TObject::kWriteDelete);		
		histodAbsRVsPhi->Write("Resolution_dRAbs_VS_Phi",TObject::kWriteDelete);		
		histodAbsRVsPt->Write("Resolution_dRAbs_VS_Pt",TObject::kWriteDelete);		
		histodAbsRVsEta->Write("Resolution_dRAbs_VS_Eta",TObject::kWriteDelete);		
		
		histodAbsZVsR->Write("Resolution_dZAbs_VS_R",TObject::kWriteDelete);		
		histodAbsZVsZ->Write("Resolution_dZAbs_VS_Z",TObject::kWriteDelete);	
		histodAbsZVsPhi->Write("Resolution_dZAbs_VS_Phi",TObject::kWriteDelete);		
		histodAbsZVsPt->Write("Resolution_dZAbs_VS_Pt",TObject::kWriteDelete);		
		histodAbsZVsEta->Write("Resolution_dZAbs_VS_Eta",TObject::kWriteDelete);		
		
		histodAbsPhiVsR->Write("Resolution_dPhiAbs_VS_R",TObject::kWriteDelete);		
		histodAbsPhiVsZ->Write("Resolution_dPhiAbs_VS_Z",TObject::kWriteDelete);		
		histodAbsPhiVsPhi->Write("Resolution_dPhiAbs_VS_Phi",TObject::kWriteDelete);
		histodAbsPhiVsPt->Write("Resolution_dPhiAbs_VS_Pt",TObject::kWriteDelete);		
		histodAbsPhiVsEta->Write("Resolution_dPhiAbs_VS_Eta",TObject::kWriteDelete);		

		histodAbsEtaVsR->Write("Resolution_dEtaAbs_VS_R",TObject::kWriteDelete);		
		histodAbsEtaVsZ->Write("Resolution_dEtaAbs_VS_Z",TObject::kWriteDelete);		
		histodAbsEtaVsPhi->Write("Resolution_dEtaAbs_VS_Phi",TObject::kWriteDelete);		
		histodAbsEtaVsPt->Write("Resolution_dEtaAbs_VS_Pt",TObject::kWriteDelete);		
		histodAbsEtaVsEta->Write("Resolution_dEtaAbs_VS_Eta",TObject::kWriteDelete);		

		histodPtVsR->Write("Resolution_dPtAbs_VS_R",TObject::kWriteDelete);		
		histodPtVsZ->Write("Resolution_dPtAbs_VS_Z",TObject::kWriteDelete);		
		histodPtVsPhi->Write("Resolution_dPtAbs_VS_Phi",TObject::kWriteDelete);		
		histodPtVsPt->Write("Resolution_dPtAbs_VS_Pt",TObject::kWriteDelete);		
		histodPtVsEta->Write("Resolution_dPtAbs_VS_Eta",TObject::kWriteDelete);		
	fileMappingWrite->Write();
	fileMappingWrite->Close();

	delete histodAbsRVsR;
	delete histodAbsRVsZ;
	delete histodAbsRVsPhi;
	delete histodAbsRVsPt;
	delete histodAbsRVsEta;
		
	delete histodAbsZVsR;
	delete histodAbsZVsZ;
	delete histodAbsZVsPhi;
	delete histodAbsZVsPt;
	delete histodAbsZVsEta;
		
	delete histodAbsPhiVsR;
	delete histodAbsPhiVsZ;
	delete histodAbsPhiVsPhi;
	delete histodAbsPhiVsPt;
	delete histodAbsPhiVsEta;

	delete histodAbsEtaVsR;
	delete histodAbsEtaVsZ;
	delete histodAbsEtaVsPhi;
	delete histodAbsEtaVsPt;
	delete histodAbsEtaVsEta;

	delete histodPtVsR;
	delete histodPtVsZ;
	delete histodPtVsPhi;
	delete histodPtVsPt;
	delete histodPtVsEta;
	
}
