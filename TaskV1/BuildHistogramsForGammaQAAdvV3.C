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
#include "TH3F.h"
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
typedef TVectorT<float> TVectorF;

using namespace std;

void SetLogBinningTH3(TH3* histoRebin){
   TAxis *axisafter 	= histoRebin->GetZaxis(); 
   Int_t bins 			= axisafter->GetNbins();
   Double_t from 		= axisafter->GetXmin();
   Double_t to 			= axisafter->GetXmax();
   Double_t *newbins 	= new Double_t[bins+1];
   newbins[0] 			= from;
   Double_t factor 		= TMath::Power(to/from, 1./bins);
   for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
   axisafter->Set(bins, newbins);
   delete [] newbins;

}

void SetLogBinningTH2(TH2* histoRebin){
   TAxis *axisafter 	= histoRebin->GetYaxis(); 
   Int_t bins 			= axisafter->GetNbins();
   Double_t from 		= axisafter->GetXmin();
   Double_t to 			= axisafter->GetXmax();
   Double_t *newbins 	= new Double_t[bins+1];
   newbins[0] 			= from;
   Double_t factor 		= TMath::Power(to/from, 1./bins);
   for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
   axisafter->Set(bins, newbins);
   delete [] newbins;
}

void SetLogBinningXTH2(TH2* histoRebin){
   TAxis *axisafter 	= histoRebin->GetXaxis(); 
   Int_t bins 			= axisafter->GetNbins();
   Double_t from 		= axisafter->GetXmin();
   Double_t to 			= axisafter->GetXmax();
   Double_t *newbins 	= new Double_t[bins+1];
   newbins[0] 			= from;
   Double_t factor 		= TMath::Power(to/from, 1./bins);
   for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
   axisafter->Set(bins, newbins);
   delete [] newbins;

}


void BuildHistogramsForGammaQAAdvV3( TString fileName 				= "GammaConvV1_QA_5460001022092970003190000000.root", 
									 TString fCutSelection 			= "5460001022092970003190000000", 
									 TString specificCutSelection 	= "0004314141",  
									 Bool_t kMC						= 0, 
									 Bool_t merge 					= 0,  
									 TString nameOutputBase 		= "PhotonQA",
									 Bool_t addExtFolderToOutput	= kFALSE
		){
   
	//Reset ROOT and connect tree file
	TString V0Finder 					= specificCutSelection(0,1);
	TString etaCutNumber 				= specificCutSelection(1,1);
	TString minPtCutNumber 				= specificCutSelection(2,1);
	TString chi2CutNumber 				= specificCutSelection(3,1);
	TString psiPairCutNumber 			= specificCutSelection(4,1);
	TString twoDimChi2PsiPairCutNumber 	= specificCutSelection(5,1);
	TString qtCutNumber 				= specificCutSelection(6,1);
	TString twoDimQtCutNumber 			= specificCutSelection(7,1);
	TString cosPointingCutNumber 		= specificCutSelection(8,1);
	TString asymmetryCutNumber 			= specificCutSelection(9,1);
	
	if ( V0Finder.CompareTo("0") == 0){
		cout << "V0Finder: " <<  "Onfly" << endl;
	} else if ( V0Finder.CompareTo("1") == 0){
		cout << "V0Finder: " <<  "Offline" << endl;
	}
	Double_t etaMinCut 				= -0.1;
	Double_t etaMaxCut 				= 0.;
	if (etaCutNumber.CompareTo("0") == 0){
		etaMaxCut 					= 0.9;
	} else if (etaCutNumber.CompareTo("1") == 0 ){
		etaMaxCut 					= 0.1;
	}else if (etaCutNumber.CompareTo("2") == 0 ){
		etaMaxCut 					= 1.4;
		etaMinCut 					= 0.9;
	}else if (etaCutNumber.CompareTo("3") == 0 ){
		etaMaxCut 					= 1.8;
		etaMinCut 					= 1.4;
	} else if (etaCutNumber.CompareTo("4") == 0 ){
		etaMaxCut 					= 1.4;
		etaMinCut 					= -0.1;
	} else if (etaCutNumber.CompareTo("5") == 0 ){
		etaMaxCut 					= 10.;
		etaMinCut 					= 2.5;
	} else if (etaCutNumber.CompareTo("6") == 0 ){
		etaMaxCut 					= 10.;
		etaMinCut 					= -0.1;
	}
	if (etaMinCut != -0.1){
		cout << "Eta range: " <<   etaMinCut  <<  " < |eta| < " << etaMaxCut << endl;
	} else {
		cout << "Eta range: " <<  "|eta| < " << etaMaxCut << endl;
	}
	
	Double_t minPtPhotonCut 		= 0.0;
	if (minPtCutNumber.CompareTo("0") == 0){
		minPtPhotonCut 				= 0.05;
	} else if (minPtCutNumber.CompareTo("1") == 0 ){
		minPtPhotonCut 				= 0.1;
	} else if (minPtCutNumber.CompareTo("2") == 0 ){
		minPtPhotonCut 				= 0.15;
	} else if (minPtCutNumber.CompareTo("3") == 0 ){
		minPtPhotonCut 				= 0.2;
	} else if (minPtCutNumber.CompareTo("4") == 0 ){
		minPtPhotonCut 				= 0.3;
	}
	cout<< "Min pt for tracks: " << minPtPhotonCut << endl;
	
	Double_t chi2PerDOFCut 			= 100000.0;
	if (chi2CutNumber.CompareTo("0") == 0){
		chi2PerDOFCut 				= 100000.;
	} else if (chi2CutNumber.CompareTo("1") == 0 ){
		chi2PerDOFCut 				= 500.;
	} else if (chi2CutNumber.CompareTo("2") == 0 ){
		chi2PerDOFCut 				= 200.;
	} else if (chi2CutNumber.CompareTo("3") == 0 ){
		chi2PerDOFCut 				= 100.;
	} else if (chi2CutNumber.CompareTo("4") == 0 ){
		chi2PerDOFCut 				= 50.;
	} else if (chi2CutNumber.CompareTo("5") == 0 ){
		chi2PerDOFCut 				= 30.;
	} else if (chi2CutNumber.CompareTo("6") == 0 ){
		chi2PerDOFCut 				= 20.;
	} else if (chi2CutNumber.CompareTo("7") == 0 ){
		chi2PerDOFCut 				= 15.;
	} else if (chi2CutNumber.CompareTo("8") == 0 ){
		chi2PerDOFCut 				= 10.;
	} else if (chi2CutNumber.CompareTo("9") == 0 ){
		chi2PerDOFCut 				= 5.;
	}
	cout<< "Chi2 per DOF Cut: " << chi2PerDOFCut << endl;
	
	Double_t psiPairCut 			= 10000.0;
	if (psiPairCutNumber.CompareTo("0") == 0){
		psiPairCut 					= 10000.;
	} else if (psiPairCutNumber.CompareTo("1") == 0 ){
		psiPairCut 					= 0.5;
	} else if (psiPairCutNumber.CompareTo("2") == 0 ){
		psiPairCut 					= 0.2;
	} else if (psiPairCutNumber.CompareTo("3") == 0 ){
		psiPairCut 					= 0.1;
	} else if (psiPairCutNumber.CompareTo("4") == 0 ){
		psiPairCut 					= 0.05;
	} else if (psiPairCutNumber.CompareTo("5") == 0 ){
		psiPairCut 					= 0.035;
	} else if (psiPairCutNumber.CompareTo("6") == 0 ){
		psiPairCut 					= 0.02;
	}
	cout<< "Psi Pair Cut: " << psiPairCut << endl;
	Bool_t f2DChi2PsiPair 			= kFALSE;
	if (twoDimChi2PsiPairCutNumber.CompareTo("0") == 0){
		f2DChi2PsiPair 				= kFALSE;
	} else if (twoDimChi2PsiPairCutNumber.CompareTo("1") == 0 ){
		f2DChi2PsiPair 				= kTRUE;
	}
	if (f2DChi2PsiPair) cout << "making 2D cut in Psi Pair and Chi2" << endl;
	
	Double_t qtCut 					= 1;
	if (qtCutNumber.CompareTo("0") == 0){
		qtCut 						= 1;
	} else if (qtCutNumber.CompareTo("1") == 0 ){
		qtCut 						= 0.15;
	} else if (qtCutNumber.CompareTo("2") == 0 ){
		qtCut 						= 0.1;
	} else if (qtCutNumber.CompareTo("3") == 0 ){
		qtCut 						= 0.07;
	} else if (qtCutNumber.CompareTo("4") == 0 ){
		qtCut 						= 0.05;
	} else if (qtCutNumber.CompareTo("5") == 0 ){
		qtCut 						= 0.03;
	} else if (qtCutNumber.CompareTo("6") == 0 ){
		qtCut 						= 0.02;
	}
	cout<< "Qt Cut: " << qtCut << endl;
	Bool_t f2DQt 					= kFALSE;
	if (twoDimQtCutNumber.CompareTo("0") == 0){
		f2DQt 						= kFALSE;
	} else if (twoDimQtCutNumber.CompareTo("1") == 0 ){
		f2DQt 						= kTRUE;
	}
	if (f2DQt) cout << "making 2D cut in qT and Alpha" << endl;
	
	Double_t cosPointingCut 		= -1;
	if (cosPointingCutNumber.CompareTo("0") == 0){
		cosPointingCut 				= -1;
	} else if (cosPointingCutNumber.CompareTo("1") == 0 ){
		cosPointingCut 				= 0;
	} else if (cosPointingCutNumber.CompareTo("2") == 0 ){
		cosPointingCut 				= 0.5;
	} else if (cosPointingCutNumber.CompareTo("3") == 0 ){
		cosPointingCut 				= 0.75;
	} else if (cosPointingCutNumber.CompareTo("4") == 0 ){
		cosPointingCut 				= 0.85;
	} else if (cosPointingCutNumber.CompareTo("5") == 0 ){
		cosPointingCut 				= 0.88;
	} else if (cosPointingCutNumber.CompareTo("6") == 0 ){
		cosPointingCut 				= 0.9;
	} else if (cosPointingCutNumber.CompareTo("7") == 0 ){
		cosPointingCut 				= 0.95;
	}
	cout<< "cos pointing Cut: " << cosPointingCut << endl;
	
	Bool_t fAsymmCut 				= kFALSE;
	if (asymmetryCutNumber.CompareTo("0") == 0){
		fAsymmCut 					= kFALSE;
	} else if (asymmetryCutNumber.CompareTo("1") == 0 ){
		fAsymmCut 					= kTRUE;
	} 
	
	//********************************************************************************
	//*            Definition of Cuts                                                *
	//********************************************************************************
	
	
	gROOT->Reset();

	//********************************************************************************
	//*            File definition/ loading                                          *
	//********************************************************************************
	
	TFile *f 						= (TFile*)gROOT->GetListOfFiles()->FindObject(fileName.Data());
	if (!f) {
		f 							= new TFile(fileName.Data());
	}
	if (!f) cout << "main List not found" << endl;
	if (f->IsZombie()) {
		cout <<fileName.Data() <<" file does not exist" << endl;
		f->Close();
		delete f;
		return;
	}   
	TString nameDirectory 			= Form("GammaConvV1_QA_%s", fCutSelection.Data());
	TDirectory* fGammaDir 			= (TDirectory*)f->Get(nameDirectory.Data());
	
	cout << Form("GammaConv_V1QA_%s",fCutSelection.Data()) << endl;
	TList* mainList 				= NULL;
	if (fGammaDir){
		mainList 					= (TList*)fGammaDir->Get(Form("GammaConv_V1QA_%s",fCutSelection.Data()));
	} else {
		mainList 					= (TList*)f->Get(Form("GammaConv_V1QA_%s",fCutSelection.Data()));
	}	
	if (!mainList) cout << "main List not found" << endl;
	TList* treeList 				= (TList*)mainList->FindObject("TreeList");
	if (!treeList) cout << "Tree List not found" << endl;
	TList* ESDQADirectory 			= (TList*)mainList->FindObject("ESD QA");					//added june 18
	if (!ESDQADirectory){															//added june 18
	cout << "ESDQADirectory List not found" << endl;}							//added june 18																		
	
	//Declaration of leaves types
	Float_t         photonPt;
	Float_t         photonTheta;
	Float_t         photonChi2Ndf;
	UChar_t          kind;
	TVectorF * daugtherProp 		= new TVectorF;
	TVectorF * gammaConvCoord 		= new TVectorF;
	TVectorF * gammaPhotonProp 		= new TVectorF;
	
	TTree *PhotonQA 				= (TTree*)treeList->FindObject("PhotonQA");
	PhotonQA->SetBranchAddress("pt",&photonPt);
	PhotonQA->SetBranchAddress("theta",&photonTheta);
	PhotonQA->SetBranchAddress("chi2ndf",&photonChi2Ndf);
	PhotonQA->SetBranchAddress("daughterProp",&daugtherProp);
	PhotonQA->SetBranchAddress("recCords",&gammaConvCoord);
	PhotonQA->SetBranchAddress("photonProp",&gammaPhotonProp);
	if (kMC) {
		PhotonQA->SetBranchAddress("kind",&kind);
	}   

	//********************************************************************************
	//*            Definition of Boundaries for Histograms                           *
	//********************************************************************************

	Int_t nBinsSigmaTPC 			= 200;
	Double_t firstBinSigmaTPC  		= -10.;
	Double_t lastBinSigmaTPC   		= 10.;

	Int_t nBinsRdetailed 			= 400;
	Double_t firstBinRdetailed  	= 0.;
	Double_t lastBinRdetailed   	= 100.;

	Int_t nBinsZdetailed 			= 960;
	Double_t firstBinZdetailed  	= -120.;
	Double_t lastBinZdetailed   	= 120.;

	Int_t nBinsX 					= 250;
	Double_t firstBinX  			= -250.;
	Double_t lastBinX				= 250.;

	Int_t nBinsY 					= 250;
	Double_t firstBinY				= -250.;
	Double_t lastBinY				= 250.;

	Int_t nBinsR 					= 180;
	Double_t firstBinR				= 0.;
	Double_t lastBinR				= 180.;
	
	Int_t nBinsSigmaTOF 			= 200;
	Double_t firstBinSigmaTOF		= -10.;
	Double_t lastBinSigmaTOF		= 10.;

	Int_t nBinsSigmaITS 			= 200;
	Double_t firstBinSigmaITS 		= -10.;
	Double_t lastBinSigmaITS		= 10.;

	//P-plots
	Int_t nBinsP					= 100;
	Double_t firstBinP				= 0.05;
	Double_t lastBinP				= 20.;

	//dEdx-plots
	Int_t nBinsdEdx					= 200;
	Double_t firstBindEdx			= 0.;
	Double_t lastBindEdx			= 200.;

	//Qt-plots
	Int_t nBinsQt					= 250;
	Double_t firstBinQt				= 0.;
	Double_t lastBinQt				= 0.25;

	//Pt-plots
	Int_t nBinsPt					= 250;
	Double_t firstBinPt				= 0.01;
	Double_t lastBinPt				= 50.;

	//TOF-plots 
	Int_t nBinsTOFsignal			= 600;
	Double_t firstBinTOFsignal		= -1000.;
	Double_t lastBinTOFsignal		= 29000.;

	//Eta-plots
	Int_t nBinsEta					= 200;
	Double_t firstBinEta			= -2.;
	Double_t lastBinEta				= 2.;

	//Eta-plots
	Int_t nBinsEtaTOF				= 80;
	Double_t firstBinEtaTOF			= -0.8;
	Double_t lastBinEtaTOF			= 0.8;

	//Phi-plots
	Int_t nBinsPhi					= 72;
	Double_t firstBinPhi			= 0;
	Double_t lastBinPhi				= 2*TMath::Pi();

	//TPCcluster to Findable-plots
	Int_t nBinsClsToF				= 150;
	Double_t firstBinClsToF			= 0.;
	Double_t lastBinClsToF			= 1.5;

	//TPCcluster 
	Int_t nBinsCls					= 200;
	Double_t firstBinCls			= 0.;
	Double_t lastBinCls				= 200;

	//ITSCluster
	Int_t nBinsClsITS				= 7;
	Double_t firstBinClsITS			= -0.5;
	Double_t lastBinClsITS			= 6.5;

	Int_t nBinsGammaR				= 240;
	Double_t firstBinGammaR			= 0;
	Double_t lastBinGammaR 			= 120.;
	
	Int_t nBinsGammaChi2 			= 500;
	Double_t firstBinGammaChi2 		= 0;
	Double_t lastBinGammaChi2 		= 200.;

	Int_t nBinsAlpha 				= 200;
	Double_t firstBinAlpha 			= -1.;
	Double_t lastBinAlpha 			= 1.;
	
	Int_t nBinsPsiPair 				= 500;
	Double_t firstBinPsiPair	 	= -0.5;
	Double_t lastBinPsiPair 		= 0.5;
	Int_t nBinsCosPoint 			= 1000;
	Double_t firstBinCosPoint 		= -1.;
	Double_t lastBinCosPoint 		= 1.;
	
	Int_t nBinsAsym 				= 200;
	Double_t firstBinAsym 			= 0.;
	Double_t lastBinAsym 			= 1.;
	
	Int_t nBinsInvMass 				= 200;
	Double_t firstBinInvMass 		= 0.;
	Double_t lastBinInvMass 		= 0.2;
	
	//********************************************************************************
	//*      Definition of histograms for reconstructed Conversion Points            *
	//********************************************************************************
	TH3F* histoElectrondEdxEtaP 				= NULL;
	TH3F* histoPositrondEdxEtaP 				= NULL;
	TH3F* histoElectronNSigmadEdxEtaP 			= NULL;
	TH3F* histoPositronNSigmadEdxEtaP 			= NULL; 
	TH3F* histoElectronTOFEtaP 					= NULL;
	TH3F* histoPositronTOFEtaP 					= NULL; 
	TH3F* histoElectronNSigmaTOFEtaP 			= NULL;
	TH3F* histoPositronNSigmaTOFEtaP 			= NULL;
	TH3F* histoElectronNSigmaITSEtaP 			= NULL; 
	TH3F* histoPositronNSigmaITSEtaP 			= NULL; 
	TH3F* histoElectronITSdEdxEtaP 				= NULL; 
	TH3F* histoPositronITSdEdxEtaP 				= NULL; 
	TH3F* histonSigdEdxElnSigdEdxPosGammaP		= NULL; // new sept 7 2015
	TH2F* histoElectronNSigmadEdxPhi 			= NULL;
	TH2F* histoPositronNSigmadEdxPhi 			= NULL;
	TH2F* histoElectronNSigmadEdxPhiEtaNeg 		= NULL;
	TH2F* histoPositronNSigmadEdxPhiEtaNeg 		= NULL;
	TH2F* histoElectronNSigmadEdxPhiEtaPos 		= NULL;
	TH2F* histoPositronNSigmadEdxPhiEtaPos 		= NULL;
	TH2F* histoElectronFClPt 					= NULL;
	TH2F* histoPositronFClPt 					= NULL;
	TH2F* histoElectronClPt 					= NULL;
	TH2F* histoPositronClPt 					= NULL;
	TH2F* histoGammaChi2NDFPt 					= NULL;
	TH2F* histoGammaChi2NDFPtEtaNeg 			= NULL;
	TH2F* histoGammaChi2NDFPtEtaPos 			= NULL;
	TH2F* histoGammaChi2NDFR 					= NULL;
	TH2F* histoElectronEtaPt 					= NULL;
	TH2F* histoPositronEtaPt 					= NULL;
	TH2F* histoElectronITSClPt 					= NULL;
	TH2F* histoPositronITSClPt 					= NULL;
	TH2F* histoElectronITSClEta 				= NULL;
	TH2F* histoPositronITSClEta 				= NULL;
	TH2F* histoPositronITSClR 					= NULL;
	TH2F* histoElectronITSClR 					= NULL;
	TH2F* histoPositronITSClPhi 				= NULL;
	TH2F* histoElectronITSClPhi 				= NULL;
	TH2F* histoPositronClR 						= NULL;
	TH2F* histoElectronClR 						= NULL;
	TH2F* histoPositronClPhi 					= NULL;
	TH2F* histoElectronClPhi 					= NULL;
	TH1F* histoGammaPhiEtaNeg 					= NULL;
	TH1F* histoGammaPhiEtaPos 					= NULL;
	TH1F* histoPositronNSigmadEdxCut 			= NULL;
	TH1F* histoElectronNSigmadEdxCut 			= NULL;
	TH1F *histoVertexZ 							= NULL;
	TH1F *histoContrVertexZ 					= NULL;
	TH1F *histoGoodESDTracks 					= NULL;
	TH1F *histoVertexZdummy 					= NULL;
	TH1F *histoContrVertexZdummy 				= NULL;
	TH1F *histoGoodESDTracksdummy 				= NULL;
	TH2F* histoGammaAsymP 						= NULL;
	TH2F* histoGammaAsymR 						= NULL;
	TH2F* histoGammaAlphaQt 					= NULL;
	TH2F* histoGammaAlphaR 						= NULL;
	TH2F* histoGammaQtPt 						= NULL;
	TH2F* histoGammaQtR 						= NULL;
	TH2F* histoGammaEtaPt 						= NULL;
	TH2F* histoGammaEtaR 						= NULL;
	TH2F* histoGammaPsiPairPt 					= NULL;
	TH2F* histoGammaPsiPairPtEtaNeg 			= NULL;
	TH2F* histoGammaPsiPairPtEtaPos 			= NULL;
	TH2F* histoGammaPsiPairR 					= NULL;
	TH2F* histoGammaCosPointPt 					= NULL;
	TH2F* histoGammaCosPointR 					= NULL;
	TH2F* histoGammaCosPointChi2 				= NULL;
	TH2F* histoGammaInvMassPt 					= NULL;
	TH2F* histoGammaInvMassR 					= NULL;
	TH2F* histoGammaCosPointPsiPair 			= NULL;
	TH2F* histoGammaChi2PsiPair 				= NULL;
	TH1F* histoGammaPhi 						= NULL;
	TH2F* histoGammaPhiR 						= NULL;
	TH2F* histoGammaZR 							= NULL;
	TH2F* histoGammaXY 							= NULL;
	
	TH2F* histoTruePrimGammaChi2NDFPt 			= NULL;
	TH2F* histoTruePrimGammaAlphaQt 			= NULL;
	TH2F* histoTruePrimGammaQtPt 				= NULL;
	TH2F* histoTruePrimGammaPsiPairPt 			= NULL;
	TH2F* histoTruePrimGammaCosPointPt 			= NULL;
	TH2F* histoTruePrimGammaAsymP 				= NULL;
	TH2F* histoTrueSecGammaChi2NDFPt 			= NULL;
	TH2F* histoTrueSecGammaAlphaQt 				= NULL;
	TH2F* histoTrueSecGammaQtPt 				= NULL;
	TH2F* histoTrueSecGammaPsiPairPt 			= NULL;
	TH2F* histoTrueSecGammaCosPointPt 			= NULL;
	TH2F* histoTrueSecGammaAsymP 				= NULL;
	TH2F* histoTrueDalitzGammaChi2NDFPt 		= NULL;
	TH2F* histoTrueDalitzGammaAlphaQt 			= NULL;
	TH2F* histoTrueDalitzGammaQtPt 				= NULL;
	TH2F* histoTrueDalitzGammaPsiPairPt 		= NULL;
	TH2F* histoTrueDalitzGammaCosPointPt 		= NULL;
	TH2F* histoTrueDalitzGammaAsymP 			= NULL;
	TH2F* histoTrueGammaCosPointChi2 			= NULL;
	TH2F* histoTrueGammaInvMassPt 				= NULL;
	TH2F* histoTrueGammaCosPointPsiPair 		= NULL;
	TH2F* histoTrueGammaChi2PsiPair 			= NULL;
	TH2F* histoTrueDalitzGammaCosPointChi2 		= NULL;
	TH2F* histoTrueDalitzGammaInvMassPt 		= NULL;
	TH2F* histoTrueDalitzGammaCosPointPsiPair 	= NULL;
	TH2F* histoTrueDalitzGammaChi2PsiPair 		= NULL;
	TH2F* histoTrueElectronNSigmadEdxTPCP 		= NULL;
	TH2F* histoTruePositronNSigmadEdxTPCP 		= NULL;
	TH2F* histoTrueElectronNSigmadEdxITSP 		= NULL;
	TH2F* histoTruePositronNSigmadEdxITSP 		= NULL;
	TH2F* histoTrueElectronNSigmaTOFP 			= NULL;
	TH2F* histoTruePositronNSigmaTOFP 			= NULL;
	TH2F* histoTruePositronITSClR 				= NULL;
	TH2F* histoTrueElectronITSClR 				= NULL;
	
	TString fileNameOutput;
	if (kMC){
		fileNameOutput = Form("%s_MC.root",nameOutputBase.Data()) ;
	} else {
		fileNameOutput = Form("%s_Data.root",nameOutputBase.Data()) ;
	}
	TFile* fileMappingDetailedConv = new TFile(fileNameOutput.Data());
	TDirectory*  directoryConv =     (TDirectory*)fileMappingDetailedConv->Get(Form("GammaConv_%s",  specificCutSelection.Data()));          
		
	if (!merge || directoryConv==0) {   
		if (ESDQADirectory)histoVertexZ 		= (TH1F*) ESDQADirectory->FindObject("Vertex_Z");			//added june 18
		if (ESDQADirectory)histoContrVertexZ 	= (TH1F*) ESDQADirectory->FindObject("ContrVertex_Z");	//added june 18
		if (ESDQADirectory)histoGoodESDTracks 	= (TH1F*) ESDQADirectory->FindObject("GoodESDTracks");

		histoElectrondEdxEtaP 					= new TH3F("histoElectrondEdxEtaP","", nBinsdEdx, firstBindEdx, lastBindEdx, nBinsEta, firstBinEta, lastBinEta, 
														   nBinsP, firstBinP, lastBinP);
		SetLogBinningTH3(histoElectrondEdxEtaP);
		histoPositrondEdxEtaP 					= new TH3F("histoPositrondEdxEtaP","", nBinsdEdx, firstBindEdx, lastBindEdx, nBinsEta, firstBinEta, lastBinEta, 
														   nBinsP, firstBinP, lastBinP);
		SetLogBinningTH3(histoPositrondEdxEtaP);
		histoElectronNSigmadEdxEtaP 			= new TH3F("histoElectronNSigmadEdxEtaP","", nBinsSigmaTPC, firstBinSigmaTPC, lastBinSigmaTPC, nBinsEta, firstBinEta, lastBinEta, 
														   nBinsP, firstBinP, lastBinP);
		SetLogBinningTH3(histoElectronNSigmadEdxEtaP);
		histoPositronNSigmadEdxEtaP 			=  new TH3F("histoPositronNSigmadEdxEtaP","", nBinsSigmaTPC, firstBinSigmaTPC, lastBinSigmaTPC, nBinsEta, firstBinEta, lastBinEta, 
															nBinsP, firstBinP, lastBinP);
		SetLogBinningTH3(histoPositronNSigmadEdxEtaP);
		histoElectronTOFEtaP 					= new TH3F("histoElectronTOFEtaP","", nBinsTOFsignal, firstBinTOFsignal, lastBinTOFsignal,nBinsEtaTOF, firstBinEtaTOF, lastBinEtaTOF, 
														   nBinsP, firstBinP, lastBinP);
		SetLogBinningTH3(histoElectronTOFEtaP);
		histoPositronTOFEtaP 					= new TH3F("histoPositronTOFEtaP","", nBinsTOFsignal, firstBinTOFsignal, lastBinTOFsignal, nBinsEtaTOF, firstBinEtaTOF, lastBinEtaTOF, 
														   nBinsP, firstBinP, lastBinP);
		SetLogBinningTH3(histoPositronTOFEtaP);
		histoElectronNSigmaTOFEtaP 				= new TH3F("histoElectronNSigmaTOFEtaP","", nBinsSigmaTOF, firstBinSigmaTOF, lastBinSigmaTOF, nBinsEtaTOF, firstBinEtaTOF, lastBinEtaTOF, 
														   nBinsP, firstBinP, lastBinP);
		SetLogBinningTH3(histoElectronNSigmaTOFEtaP);
		histoPositronNSigmaTOFEtaP 				= new TH3F("histoPositronNSigmaTOFEtaP","", nBinsSigmaTOF, firstBinSigmaTOF, lastBinSigmaTOF, nBinsEtaTOF, firstBinEtaTOF, lastBinEtaTOF, 
														   nBinsP, firstBinP, lastBinP);
		SetLogBinningTH3(histoPositronNSigmaTOFEtaP);
		histoElectronITSdEdxEtaP 				= new TH3F("histoElectronITSdEdxEtaP","", nBinsdEdx, firstBindEdx, lastBindEdx, nBinsEta, firstBinEta, lastBinEta, 
														   nBinsP, firstBinP, lastBinP);
		SetLogBinningTH3(histoElectronITSdEdxEtaP);
		histoPositronITSdEdxEtaP 				= new TH3F("histoPositronITSdEdxEtaP","", nBinsdEdx, firstBindEdx, lastBindEdx, nBinsEta, firstBinEta, lastBinEta, 
														   nBinsP, firstBinP, lastBinP);
		SetLogBinningTH3(histoPositronITSdEdxEtaP);
		histonSigdEdxElnSigdEdxPosGammaP 		= new TH3F("histonSigdEdxElnSigdEdxPosGammaP","", nBinsSigmaTPC, firstBinSigmaTPC, lastBinSigmaTPC, nBinsSigmaTPC, firstBinSigmaTPC, 												  lastBinSigmaTPC, nBinsPt, firstBinPt, lastBinPt);  // new sept 7 2015
		histoElectronNSigmaITSEtaP 				= new TH3F("histoElectronNSigmaITSEtaP","", nBinsSigmaITS, firstBinSigmaITS, lastBinSigmaITS, nBinsEta, firstBinEta, lastBinEta, 
														   nBinsP, firstBinP, lastBinP);
		SetLogBinningTH3(histoElectronNSigmaITSEtaP);
		histoPositronNSigmaITSEtaP 				= new TH3F("histoPositronNSigmaITSEtaP","", nBinsSigmaITS, firstBinSigmaITS, lastBinSigmaITS, nBinsEta, firstBinEta, lastBinEta, 
														   nBinsP, firstBinP, lastBinP);
		SetLogBinningTH3(histoPositronNSigmaITSEtaP);
		histoElectronNSigmadEdxPhi 				= new TH2F("histoElectronNSigmadEdxPhi","",nBinsSigmaTPC,firstBinSigmaTPC, lastBinSigmaTPC, nBinsPhi, firstBinPhi, lastBinPhi);
		histoPositronNSigmadEdxPhi 				= new TH2F("histoPositronNSigmadEdxPhi","",nBinsSigmaTPC,firstBinSigmaTPC, lastBinSigmaTPC, nBinsPhi, firstBinPhi, lastBinPhi);
		histoElectronNSigmadEdxPhiEtaNeg 		= new TH2F("histoElectronNSigmadEdxPhiEtaNeg","",nBinsSigmaTPC,firstBinSigmaTPC, lastBinSigmaTPC, nBinsPhi, firstBinPhi, lastBinPhi);
		histoPositronNSigmadEdxPhiEtaNeg 		= new TH2F("histoPositronNSigmadEdxPhiEtaNeg","",nBinsSigmaTPC,firstBinSigmaTPC, lastBinSigmaTPC, nBinsPhi, firstBinPhi, lastBinPhi);
		histoElectronNSigmadEdxPhiEtaPos 		= new TH2F("histoElectronNSigmadEdxPhiEtaPos","",nBinsSigmaTPC,firstBinSigmaTPC, lastBinSigmaTPC, nBinsPhi, firstBinPhi, lastBinPhi);
		histoPositronNSigmadEdxPhiEtaPos 		= new TH2F("histoPositronNSigmadEdxPhiEtaPos","",nBinsSigmaTPC,firstBinSigmaTPC, lastBinSigmaTPC, nBinsPhi, firstBinPhi, lastBinPhi);
		histoElectronITSClPt 					= new TH2F("histoElectronITSClPt","", nBinsClsITS, firstBinClsITS,    lastBinClsITS, nBinsPt, firstBinPt, lastBinPt);
		SetLogBinningTH2(histoElectronITSClPt);
		histoPositronITSClPt 					= new TH2F("histoPositronITSClPt","", nBinsClsITS, firstBinClsITS,    lastBinClsITS, nBinsPt, firstBinPt, lastBinPt);
		SetLogBinningTH2(histoPositronITSClPt);
		histoElectronITSClEta 					= new TH2F("histoElectronITSClEta","", nBinsClsITS, firstBinClsITS, lastBinClsITS, nBinsEta, firstBinEta, lastBinEta);
		histoPositronITSClEta 					= new TH2F("histoPositronITSClEta","", nBinsClsITS, firstBinClsITS, lastBinClsITS, nBinsEta, firstBinEta, lastBinEta);
		histoElectronFClPt 						= new TH2F("histoElectronFClPt","", nBinsClsToF, firstBinClsToF, lastBinClsToF,  nBinsPt, firstBinPt, lastBinPt);
		SetLogBinningTH2(histoElectronFClPt);
		histoPositronFClPt 						= new TH2F("histoPositronFClPt","", nBinsClsToF, firstBinClsToF, lastBinClsToF,  nBinsPt, firstBinPt, lastBinPt);
		SetLogBinningTH2(histoPositronFClPt);
		histoElectronClPt 						= new TH2F("histoElectronClPt","", nBinsCls, firstBinCls,    lastBinCls, nBinsPt, firstBinPt, lastBinPt);
		SetLogBinningTH2(histoElectronClPt);
		histoPositronClPt 						= new TH2F("histoPositronClPt","", nBinsCls, firstBinCls,    lastBinCls, nBinsPt, firstBinPt, lastBinPt);
		SetLogBinningTH2(histoPositronClPt);
		histoGammaChi2NDFPt 					= new TH2F("histoGammaChi2NDFPt","", nBinsGammaChi2, firstBinGammaChi2, lastBinGammaChi2, nBinsPt, firstBinPt, lastBinPt);
		SetLogBinningTH2(histoGammaChi2NDFPt);
		histoGammaChi2NDFPtEtaNeg 				= new TH2F("histoGammaChi2NDFPtEtaNeg","", nBinsGammaChi2, firstBinGammaChi2, lastBinGammaChi2, nBinsPt, firstBinPt, lastBinPt);
		SetLogBinningTH2(histoGammaChi2NDFPtEtaNeg);
		histoGammaChi2NDFPtEtaPos 				= new TH2F("histoGammaChi2NDFPtEtaPos","", nBinsGammaChi2, firstBinGammaChi2, lastBinGammaChi2, nBinsPt, firstBinPt, lastBinPt);
		SetLogBinningTH2(histoGammaChi2NDFPtEtaPos);
		histoGammaChi2NDFR 						= new TH2F("histoGammaChi2NDFR","", nBinsGammaChi2, firstBinGammaChi2, lastBinGammaChi2, nBinsR, firstBinR, lastBinR);
		histoElectronEtaPt 						= new TH2F("histoElectronEtaPt","", nBinsEta, firstBinEta, lastBinEta, nBinsPt, firstBinPt, lastBinPt);
		SetLogBinningTH2(histoElectronEtaPt);
		histoPositronEtaPt 						= new TH2F("histoPositronEtaPt","", nBinsEta, firstBinEta, lastBinEta, nBinsPt, firstBinPt, lastBinPt);
		SetLogBinningTH2(histoPositronEtaPt);
		histoGammaEtaPt 						= new TH2F("histoGammaEtaPt","", nBinsEta, firstBinEta, lastBinEta, nBinsPt, firstBinPt, lastBinPt);
		SetLogBinningTH2(histoGammaEtaPt);	
		histoGammaEtaR 							= new TH2F("histoGammaEtaR","", nBinsEta, firstBinEta, lastBinEta, nBinsRdetailed, firstBinRdetailed, lastBinRdetailed);
		histoElectronITSClR 					= new TH2F("histoElectronITSClR","",nBinsClsITS, firstBinClsITS, lastBinClsITS, nBinsGammaR, firstBinGammaR, lastBinGammaR);
		histoPositronITSClR 					= new TH2F("histoPositronITSClR","",nBinsClsITS, firstBinClsITS, lastBinClsITS, nBinsGammaR, firstBinGammaR, lastBinGammaR);
		histoElectronITSClPhi 					= new TH2F("histoElectronITSClPhi","",nBinsClsITS, firstBinClsITS, lastBinClsITS, nBinsPhi, firstBinPhi, lastBinPhi);
		histoPositronITSClPhi 					= new TH2F("histoPositronITSClPhi","",nBinsClsITS, firstBinClsITS, lastBinClsITS, nBinsPhi, firstBinPhi, lastBinPhi);
		histoElectronClR 						= new TH2F("histoElectronClR","",nBinsCls, firstBinCls, lastBinCls, nBinsGammaR, firstBinGammaR, lastBinGammaR);
		histoPositronClR 						= new TH2F("histoPositronClR","",nBinsCls, firstBinCls, lastBinCls, nBinsGammaR, firstBinGammaR, lastBinGammaR);
		histoElectronClPhi 						= new TH2F("histoElectronClPhi","",nBinsCls, firstBinCls, lastBinCls, nBinsPhi, firstBinPhi, lastBinPhi);
		histoPositronClPhi 						= new TH2F("histoPositronClPhi","",nBinsCls, firstBinCls, lastBinCls, nBinsPhi, firstBinPhi, lastBinPhi);
		histoGammaPhiEtaNeg 					= new TH1F("histoGammaPhiEtaNeg","", nBinsPhi, firstBinPhi, lastBinPhi);		
		histoGammaPhiEtaPos 					= new TH1F("histoGammaPhiEtaPos","", nBinsPhi, firstBinPhi, lastBinPhi);      	
		histoPositronNSigmadEdxCut 				= new TH1F("histoPositronNSigmadEdxCut","", nBinsSigmaTPC, firstBinSigmaTPC, lastBinSigmaTPC);	
		histoElectronNSigmadEdxCut 				= new TH1F("histoElectronNSigmadEdxCut","", nBinsSigmaTPC, firstBinSigmaTPC, lastBinSigmaTPC);	
		histoGammaAlphaQt 						= new TH2F("histoGammaAlphaQt","", nBinsAlpha, firstBinAlpha, lastBinAlpha, nBinsQt, firstBinQt, lastBinQt);
		histoGammaAlphaR 						= new TH2F("histoGammaAlphaR","", nBinsAlpha, firstBinAlpha, lastBinAlpha, nBinsR, firstBinR, lastBinR);
		histoGammaQtPt 							= new TH2F("histoGammaQtPt","", nBinsQt, firstBinQt, lastBinQt,nBinsPt, firstBinPt, lastBinPt);
		SetLogBinningTH2(histoGammaQtPt);
		histoGammaQtR 							= new TH2F("histoGammaQtR","", nBinsQt, firstBinQt, lastBinQt,nBinsR, firstBinR, lastBinR);
		histoGammaPhi 							= new TH1F("histoGammaPhi","", nBinsPhi, firstBinPhi, lastBinPhi);      
		histoGammaPhiR 							= new TH2F("histoGammaPhiR","", nBinsPhi, firstBinPhi, lastBinPhi,nBinsRdetailed, firstBinRdetailed, lastBinRdetailed);
		histoGammaPsiPairPt 					= new TH2F("histoGammaPsiPairPt","", nBinsPsiPair, firstBinPsiPair, lastBinPsiPair,nBinsPt, firstBinPt, lastBinPt);
		SetLogBinningTH2(histoGammaPsiPairPt);
		histoGammaPsiPairPtEtaNeg 				= new TH2F("histoGammaPsiPairPtEtaNeg","", nBinsPsiPair, firstBinPsiPair, lastBinPsiPair,nBinsPt, firstBinPt, lastBinPt);
		SetLogBinningTH2(histoGammaPsiPairPtEtaNeg);
		histoGammaPsiPairPtEtaPos 				= new TH2F("histoGammaPsiPairPtEtaPos","", nBinsPsiPair, firstBinPsiPair, lastBinPsiPair,nBinsPt, firstBinPt, lastBinPt);
		SetLogBinningTH2(histoGammaPsiPairPtEtaPos);
		histoGammaPsiPairR 						= new TH2F("histoGammaPsiPairR","", nBinsPsiPair, firstBinPsiPair, lastBinPsiPair,nBinsR, firstBinR, lastBinR);
		
		histoGammaCosPointPt 					= new TH2F("histoGammaCosPointPt","", nBinsCosPoint, firstBinCosPoint, lastBinCosPoint,nBinsPt, firstBinPt, lastBinPt);        
		SetLogBinningTH2(histoGammaCosPointPt);
		histoGammaChi2PsiPair 					= new TH2F("histoGammaChi2PsiPair","",nBinsGammaChi2, firstBinGammaChi2, lastBinGammaChi2, nBinsPsiPair, firstBinPsiPair, lastBinPsiPair);
		histoGammaCosPointChi2 					= new TH2F("histoGammaCosPointChi2","",nBinsCosPoint, firstBinCosPoint, lastBinCosPoint,nBinsGammaChi2, firstBinGammaChi2, lastBinGammaChi2);
		histoGammaCosPointR 					= new TH2F("histoGammaCosPointR","",nBinsCosPoint, firstBinCosPoint, lastBinCosPoint,nBinsR, firstBinR, lastBinR);
		histoGammaCosPointPsiPair 				= new TH2F("histoGammaCosPointPsiPair","",nBinsCosPoint, firstBinCosPoint, lastBinCosPoint,nBinsPsiPair, firstBinPsiPair, lastBinPsiPair);
		histoGammaInvMassPt 					= new TH2F( "histoGammaInvMassPt","",nBinsInvMass, firstBinInvMass, lastBinInvMass, nBinsPt, firstBinPt, lastBinPt);
		SetLogBinningTH2(histoGammaInvMassPt);
		histoGammaInvMassR 						= new TH2F( "histoGammaInvMassR","",nBinsInvMass, firstBinInvMass, lastBinInvMass, nBinsR, firstBinR, lastBinR);
		histoGammaAsymP 						= new TH2F("histoGammaAsymP", "",nBinsAsym,firstBinAsym,lastBinAsym,nBinsPt, firstBinPt, lastBinPt); 
		SetLogBinningTH2(histoGammaAsymP);
		histoGammaAsymR 						= new TH2F("histoGammaAsymR", "",nBinsAsym,firstBinAsym,lastBinAsym,nBinsR, firstBinR, lastBinR); 
		histoGammaZR 							= new TH2F("histoGammaZR", "",nBinsZdetailed,firstBinZdetailed,lastBinZdetailed,nBinsRdetailed, firstBinRdetailed, lastBinRdetailed); 
		histoGammaXY 							= new TH2F("histoGammaXY", "",nBinsX,firstBinX,lastBinX,nBinsY, firstBinY, lastBinY); 
		
		if (kMC){
			histoTrueElectronNSigmadEdxTPCP 	= new TH2F("histoTrueElectronNSigmadEdxTPCP","", nBinsP, firstBinP, lastBinP, nBinsSigmaTPC, firstBinSigmaTPC, lastBinSigmaTPC);
			SetLogBinningXTH2(histoTrueElectronNSigmadEdxTPCP);
			histoTruePositronNSigmadEdxTPCP 	= new TH2F("histoTruePositronNSigmadEdxTPCP","", nBinsP, firstBinP, lastBinP, nBinsSigmaTPC, firstBinSigmaTPC, lastBinSigmaTPC);
			SetLogBinningXTH2(histoTruePositronNSigmadEdxTPCP);
			histoTrueElectronNSigmadEdxITSP 	= new TH2F("histoTrueElectronNSigmadEdxITSP","", nBinsP, firstBinP, lastBinP, nBinsSigmaITS, firstBinSigmaITS, lastBinSigmaITS);
			SetLogBinningXTH2(histoTrueElectronNSigmadEdxITSP);
			histoTruePositronNSigmadEdxITSP 	= new TH2F("histoTruePositronNSigmadEdxITSP","", nBinsP, firstBinP, lastBinP, nBinsSigmaITS, firstBinSigmaITS, lastBinSigmaITS);
			SetLogBinningXTH2(histoTruePositronNSigmadEdxITSP);
			histoTrueElectronNSigmaTOFP 		= new TH2F("histoTrueElectronNSigmaTOFP","", nBinsP, firstBinP, lastBinP, nBinsSigmaTOF, firstBinSigmaTOF, lastBinSigmaTOF);
			SetLogBinningXTH2(histoTrueElectronNSigmaTOFP);
			histoTruePositronNSigmaTOFP 		= new TH2F("histoTruePositronNSigmaTOFP","", nBinsP, firstBinP, lastBinP, nBinsSigmaTOF, firstBinSigmaTOF, lastBinSigmaTOF);
			SetLogBinningXTH2(histoTruePositronNSigmaTOFP);
			histoTrueElectronITSClR 			= new TH2F("histoTrueElectronITSClR","",nBinsClsITS, firstBinClsITS, lastBinClsITS, nBinsGammaR, firstBinGammaR, lastBinGammaR);
			histoTruePositronITSClR 			= new TH2F("histoTruePositronITSClR","",nBinsClsITS, firstBinClsITS, lastBinClsITS, nBinsGammaR, firstBinGammaR, lastBinGammaR);
			histoTruePrimGammaChi2NDFPt 		= new TH2F("histoTruePrimGammaChi2NDFPt","", nBinsGammaChi2, firstBinGammaChi2, lastBinGammaChi2, nBinsPt, firstBinPt, lastBinPt);
			SetLogBinningTH2(histoTruePrimGammaChi2NDFPt);
			histoTruePrimGammaAlphaQt 			= new TH2F("histoTruePrimGammaAlphaQt","", nBinsAlpha, firstBinAlpha, lastBinAlpha, nBinsQt, firstBinQt, lastBinQt);
			histoTruePrimGammaQtPt 				= new TH2F("histoTruePrimGammaQtPt","", nBinsQt, firstBinQt, lastBinQt,nBinsPt, firstBinPt, lastBinPt);
			SetLogBinningTH2(histoTruePrimGammaQtPt);
			histoTruePrimGammaPsiPairPt 		= new TH2F( "histoTruePrimGammaPsiPairPt","", nBinsPsiPair, firstBinPsiPair, lastBinPsiPair,nBinsPt, firstBinPt, lastBinPt);        
			SetLogBinningTH2(histoTruePrimGammaPsiPairPt);
			histoTruePrimGammaCosPointPt 		= new TH2F( "histoTruePrimGammaCosPointPt","", nBinsCosPoint, firstBinCosPoint, lastBinCosPoint,nBinsPt, firstBinPt, lastBinPt);        
			SetLogBinningTH2(histoTruePrimGammaCosPointPt);
			histoTruePrimGammaAsymP 			= new TH2F("histoTruePrimGammaAsymP", "",nBinsAsym,firstBinAsym,lastBinAsym,nBinsPt, firstBinPt, lastBinPt); 
			SetLogBinningTH2(histoTruePrimGammaAsymP);
			histoTrueSecGammaChi2NDFPt 			= new TH2F("histoTrueSecGammaChi2NDFPt","", nBinsGammaChi2, firstBinGammaChi2, lastBinGammaChi2, nBinsPt, firstBinPt, lastBinPt);
			SetLogBinningTH2(histoTrueSecGammaChi2NDFPt);
			histoTrueSecGammaAlphaQt 			= new TH2F("histoTrueSecGammaAlphaQt","", nBinsAlpha, firstBinAlpha, lastBinAlpha, nBinsQt, firstBinQt, lastBinQt);
			histoTrueSecGammaQtPt 				= new TH2F("histoTrueSecGammaQtPt","", nBinsQt, firstBinQt, lastBinQt,nBinsPt, firstBinPt, lastBinPt);
			SetLogBinningTH2(histoTrueSecGammaQtPt);
			histoTrueSecGammaPsiPairPt 			= new TH2F( "histoTrueSecGammaPsiPairPt","", nBinsPsiPair, firstBinPsiPair, lastBinPsiPair,nBinsPt, firstBinPt, lastBinPt);        
			SetLogBinningTH2(histoTrueSecGammaPsiPairPt);
			histoTrueSecGammaCosPointPt 		= new TH2F( "histoTrueSecGammaCosPointPt","", nBinsCosPoint, firstBinCosPoint, lastBinCosPoint,nBinsPt, firstBinPt, lastBinPt);        
			SetLogBinningTH2(histoTrueSecGammaCosPointPt);
			histoTrueSecGammaAsymP 				= new TH2F("histoTrueSecGammaAsymP", "",nBinsAsym,firstBinAsym,lastBinAsym,nBinsPt, firstBinPt, lastBinPt); 
			SetLogBinningTH2(histoTrueSecGammaAsymP);
			histoTrueGammaChi2PsiPair 			= new TH2F( "histoTrueGammaChi2PsiPair","",nBinsGammaChi2, firstBinGammaChi2, lastBinGammaChi2, nBinsPsiPair, firstBinPsiPair, lastBinPsiPair);
			histoTrueGammaCosPointChi2 			= new TH2F( "histoTrueGammaCosPointChi2","",nBinsCosPoint, firstBinCosPoint, lastBinCosPoint,nBinsGammaChi2, firstBinGammaChi2, lastBinGammaChi2);
			histoTrueGammaCosPointPsiPair 		= new TH2F( "histoTrueGammaCosPointPsiPair","",nBinsCosPoint, firstBinCosPoint, lastBinCosPoint,nBinsPsiPair, firstBinPsiPair, lastBinPsiPair);
			histoTrueGammaInvMassPt 			= new TH2F( "histoTrueGammaInvMassPt","",nBinsInvMass, firstBinInvMass, lastBinInvMass, nBinsPt, firstBinPt, lastBinPt);
			SetLogBinningTH2(histoTrueGammaInvMassPt);
			histoTrueDalitzGammaChi2NDFPt 		= new TH2F("histoTrueDalitzGammaChi2NDFPt","", nBinsGammaChi2, firstBinGammaChi2, lastBinGammaChi2, nBinsPt, firstBinPt, lastBinPt);
			SetLogBinningTH2(histoTrueDalitzGammaChi2NDFPt);
			histoTrueDalitzGammaAlphaQt 		= new TH2F("histoTrueDalitzGammaAlphaQt","", nBinsAlpha, firstBinAlpha, lastBinAlpha, nBinsQt, firstBinQt, lastBinQt);
			histoTrueDalitzGammaQtPt 			= new TH2F("histoTrueDalitzGammaQtPt","", nBinsQt, firstBinQt, lastBinQt,nBinsPt, firstBinPt, lastBinPt);
			SetLogBinningTH2(histoTrueDalitzGammaQtPt);
			histoTrueDalitzGammaPsiPairPt 		= new TH2F( "histoTrueDalitzGammaPsiPairPt","", nBinsPsiPair, firstBinPsiPair, lastBinPsiPair,nBinsPt, firstBinPt, lastBinPt);        
			SetLogBinningTH2(histoTrueDalitzGammaPsiPairPt);
			histoTrueDalitzGammaCosPointPt 		= new TH2F( "histoTrueDalitzGammaCosPointPt","", nBinsCosPoint, firstBinCosPoint, lastBinCosPoint,nBinsPt, firstBinPt, lastBinPt);        
			SetLogBinningTH2(histoTrueDalitzGammaCosPointPt);
			histoTrueDalitzGammaAsymP 			= new TH2F("histoTrueDalitzGammaAsymP", "",nBinsAsym,firstBinAsym,lastBinAsym,nBinsPt, firstBinPt, lastBinPt); 
			SetLogBinningTH2(histoTrueDalitzGammaAsymP);
			histoTrueDalitzGammaChi2PsiPair 	= new TH2F( "histoTrueDalitzGammaChi2PsiPair","",nBinsGammaChi2, firstBinGammaChi2, lastBinGammaChi2, nBinsPsiPair, firstBinPsiPair, lastBinPsiPair);
			histoTrueDalitzGammaCosPointChi2 	= new TH2F( "histoTrueDalitzGammaCosPointChi2","",nBinsCosPoint, firstBinCosPoint, lastBinCosPoint,nBinsGammaChi2, firstBinGammaChi2, lastBinGammaChi2);
			histoTrueDalitzGammaCosPointPsiPair = new TH2F( "histoTrueDalitzGammaCosPointPsiPair","",nBinsCosPoint, firstBinCosPoint, lastBinCosPoint,nBinsPsiPair, firstBinPsiPair, lastBinPsiPair);
			histoTrueDalitzGammaInvMassPt 		= new TH2F( "histoTrueDalitzGammaInvMassPt","",nBinsInvMass, firstBinInvMass, lastBinInvMass, nBinsPt, firstBinPt, lastBinPt);
			SetLogBinningTH2(histoTrueDalitzGammaInvMassPt);
		}   
	} else {
		// reading additional histograms from histogram branch of QA
		if (ESDQADirectory){
			histoVertexZ 						= (TH1F*) directoryConv->Get("histoVertexZ");						
			histoVertexZ->Sumw2();
			histoContrVertexZ 					= (TH1F*) directoryConv->Get("histoContrVertexZ");			
			histoContrVertexZ->Sumw2();
			histoGoodESDTracks 					= (TH1F*) directoryConv->Get("histoGoodESDTracks");						
			histoGoodESDTracks->Sumw2();
			histoVertexZdummy 					= (TH1F*) ESDQADirectory->FindObject("Vertex_Z");				
			histoVertexZdummy->Sumw2();
			histoContrVertexZdummy 				= (TH1F*) ESDQADirectory->FindObject("ContrVertex_Z");	
			histoContrVertexZdummy->Sumw2();
			histoGoodESDTracksdummy 			= (TH1F*) ESDQADirectory->FindObject("GoodESDTracks");				
			histoGoodESDTracksdummy->Sumw2();
			histoVertexZ->Add(histoVertexZdummy);											
			histoContrVertexZ->Add(histoContrVertexZdummy);	
			histoGoodESDTracks->Add(histoGoodESDTracksdummy);								
		}	
		histoElectrondEdxEtaP 					= (TH3F*)directoryConv->Get("histoElectrondEdxEtaP");
		histoPositrondEdxEtaP 					= (TH3F*)directoryConv->Get("histoPositrondEdxEtaP");
		histoElectronNSigmadEdxEtaP 			= (TH3F*)directoryConv->Get("histoElectronNSigmadEdxEtaP");
		histoPositronNSigmadEdxEtaP 			= (TH3F*)directoryConv->Get("histoPositronNSigmadEdxEtaP");
		histoElectronTOFEtaP 					= (TH3F*)directoryConv->Get("histoElectronTOFEtaP");
		histoPositronTOFEtaP 					= (TH3F*)directoryConv->Get("histoPositronTOFEtaP");
		histoElectronNSigmaTOFEtaP 				= (TH3F*)directoryConv->Get("histoElectronNSigmaTOFEtaP");
		histoPositronNSigmaTOFEtaP 				= (TH3F*)directoryConv->Get("histoPositronNSigmaTOFEtaP");
		histoElectronNSigmadEdxPhi 				= (TH2F*)directoryConv->Get("histoElectronNSigmadEdxPhi");
		histoPositronNSigmadEdxPhi 				= (TH2F*)directoryConv->Get("histoPositronNSigmadEdxPhi");
		histoElectronNSigmadEdxPhiEtaNeg 		= (TH2F*)directoryConv->Get("histoElectronNSigmadEdxPhiEtaNeg");
		histoPositronNSigmadEdxPhiEtaNeg 		= (TH2F*)directoryConv->Get("histoPositronNSigmadEdxPhiEtaNeg");
		histoElectronNSigmadEdxPhiEtaPos 		= (TH2F*)directoryConv->Get("histoElectronNSigmadEdxPhiEtaPos");
		histoPositronNSigmadEdxPhiEtaPos 		= (TH2F*)directoryConv->Get("histoPositronNSigmadEdxPhiEtaPos");
		histoElectronFClPt 						= (TH2F*)directoryConv->Get("histoElectronFClPt");
		histoPositronFClPt 						= (TH2F*)directoryConv->Get("histoPositronFClPt");
		histoElectronClPt 						= (TH2F*)directoryConv->Get("histoElectronClPt");
		histoPositronClPt 						= (TH2F*)directoryConv->Get("histoPositronClPt");
		histoGammaChi2NDFPt 					= (TH2F*)directoryConv->Get("histoGammaChi2NDFPt");
		histoGammaChi2NDFPtEtaNeg 				= (TH2F*)directoryConv->Get("histoGammaChi2NDFPtEtaNeg");
		histoGammaChi2NDFPtEtaPos 				= (TH2F*)directoryConv->Get("histoGammaChi2NDFPtEtaPos");
		histoGammaChi2NDFR 						= (TH2F*)directoryConv->Get("histoGammaChi2NDFR");
		histoElectronEtaPt 						= (TH2F*)directoryConv->Get("histoElectronEtaPt");
		histoPositronEtaPt 						= (TH2F*)directoryConv->Get("histoPositronEtaPt");
		histoElectronITSClR 					= (TH2F*)directoryConv->Get("histoElectronITSClR");
		histoPositronITSClR 					= (TH2F*)directoryConv->Get("histoPositronITSClR");
		histoElectronITSClPhi 					= (TH2F*)directoryConv->Get("histoElectronITSClPhi");
		histoPositronITSClPhi 					= (TH2F*)directoryConv->Get("histoPositronITSClPhi");
		histoElectronClR 						= (TH2F*)directoryConv->Get("histoElectronClR");
		histoPositronClR 						= (TH2F*)directoryConv->Get("histoPositronClR");
		histoElectronClPhi 						= (TH2F*)directoryConv->Get("histoElectronClPhi");
		histoPositronClPhi 						= (TH2F*)directoryConv->Get("histoPositronClPhi");
		histoGammaPhiEtaNeg 					= (TH1F*)directoryConv->Get("histoGammaPhiEtaNeg");
		histoGammaPhiEtaPos 					= (TH1F*)directoryConv->Get("histoGammaPhiEtaPos");
		histoPositronNSigmadEdxCut 				= (TH1F*)directoryConv->Get("histoPositronNSigmadEdxCut");
		histoElectronNSigmadEdxCut 				= (TH1F*)directoryConv->Get("histoElectronNSigmadEdxCut");
		histoGammaEtaPt 						= (TH2F*)directoryConv->Get("histoGammaEtaPt");
		histoGammaEtaR 							= (TH2F*)directoryConv->Get("histoGammaEtaR");
		histoGammaAlphaQt		 				= (TH2F*)directoryConv->Get("histoGammaAlphaQt");
		histoGammaAlphaR 						= (TH2F*)directoryConv->Get("histoGammaAlphaR");
		histoGammaQtPt 							= (TH2F*)directoryConv->Get("histoGammaQtPt");
		histoGammaQtR 							= (TH2F*)directoryConv->Get("histoGammaQtR");
		histoGammaPhi 							= (TH1F*)directoryConv->Get( "histoGammaPhi");
		histoGammaPhiR 							= (TH2F*)directoryConv->Get( "histoGammaPhiR");
		histoGammaPsiPairPt 					= (TH2F*)directoryConv->Get( "histoGammaPsiPairPt");
		histoGammaPsiPairPtEtaNeg 				= (TH2F*)directoryConv->Get( "histoGammaPsiPairPtEtaNeg");
		histoGammaPsiPairPtEtaPos 				= (TH2F*)directoryConv->Get( "histoGammaPsiPairPtEtaPos");
		histoGammaPsiPairR 						= (TH2F*)directoryConv->Get( "histoGammaPsiPairR");
		histoGammaCosPointPt 					= (TH2F*)directoryConv->Get( "histoGammaCosPointPt");
		histoGammaCosPointR 					= (TH2F*)directoryConv->Get( "histoGammaCosPointR");
		histoGammaAsymP 						= (TH2F*)directoryConv->Get("histoGammaAsymP");
		histoGammaAsymR 						= (TH2F*)directoryConv->Get("histoGammaAsymR");
		histoGammaChi2PsiPair 					= (TH2F*)directoryConv->Get("histoGammaChi2PsiPair");
		histoGammaCosPointChi2 					= (TH2F*)directoryConv->Get("histoGammaCosPointChi2");
		histoGammaCosPointPsiPair 				= (TH2F*)directoryConv->Get("histoGammaCosPointPsiPair");
		histoGammaInvMassPt 					= (TH2F*)directoryConv->Get("histoGammaInvMassPt");
		histoGammaInvMassR 						= (TH2F*)directoryConv->Get("histoGammaInvMassR");
		histoGammaZR 							= (TH2F*)directoryConv->Get("histoGammaZR");
		histoGammaXY 							= (TH2F*)directoryConv->Get("histoGammaXY");
		
		//ITS add
		histoElectronITSdEdxEtaP 				= (TH3F*)directoryConv->Get("histoElectronITSdEdxEtaP");
		histoPositronITSdEdxEtaP 				= (TH3F*)directoryConv->Get("histoPositronITSdEdxEtaP");
		histonSigdEdxElnSigdEdxPosGammaP 		= (TH3F*)directoryConv->Get("histonSigdEdxElnSigdEdxPosGammaP"); // new sept 7 2015
		histoElectronNSigmaITSEtaP 				= (TH3F*)directoryConv->Get("histoElectronNSigmaITSEtaP");
		histoPositronNSigmaITSEtaP 				= (TH3F*)directoryConv->Get("histoPositronNSigmaITSEtaP");
		histoElectronITSClPt 					= (TH2F*)directoryConv->Get("histoElectronITSClPt");
		histoPositronITSClPt 					= (TH2F*)directoryConv->Get("histoPositronITSClPt");
		histoElectronITSClEta 					= (TH2F*)directoryConv->Get("histoElectronITSClEta");
		histoPositronITSClEta 					= (TH2F*)directoryConv->Get("histoPositronITSClEta");
		
		if (kMC){
			histoTrueElectronNSigmadEdxITSP 	= (TH2F*)directoryConv->Get("histoTrueElectronNSigmadEdxITSP");
			histoTrueElectronNSigmadEdxTPCP 	= (TH2F*)directoryConv->Get("histoTrueElectronNSigmadEdxTPCP");
			histoTrueElectronNSigmaTOFP 		= (TH2F*)directoryConv->Get("histoTrueElectronNSigmaTOFP");
			histoTrueElectronITSClR 			= (TH2F*)directoryConv->Get("histoTrueElectronITSClR");
			histoTruePositronNSigmadEdxITSP 	= (TH2F*)directoryConv->Get("histoTruePositronNSigmadEdxITSP");
			histoTruePositronNSigmadEdxTPCP 	= (TH2F*)directoryConv->Get("histoTruePositronNSigmadEdxTPCP");
			histoTruePositronNSigmaTOFP 		= (TH2F*)directoryConv->Get("histoTruePositronNSigmaTOFP");
			histoTruePositronITSClR 			= (TH2F*)directoryConv->Get("histoTruePositronITSClR");
			histoTruePrimGammaChi2NDFPt 		= (TH2F*)directoryConv->Get("histoTruePrimGammaChi2NDFPt");
			histoTruePrimGammaAlphaQt 			= (TH2F*)directoryConv->Get("histoTruePrimGammaAlphaQt");
			histoTruePrimGammaQtPt 				= (TH2F*)directoryConv->Get("histoTruePrimGammaQtPt");
			histoTruePrimGammaPsiPairPt 		= (TH2F*)directoryConv->Get("histoTruePrimGammaPsiPairPt");
			histoTruePrimGammaCosPointPt 		= (TH2F*)directoryConv->Get("histoTruePrimGammaCosPointPt");
			histoTrueSecGammaChi2NDFPt 			= (TH2F*)directoryConv->Get("histoTrueSecGammaChi2NDFPt");
			histoTrueSecGammaAlphaQt 			= (TH2F*)directoryConv->Get("histoTrueSecGammaAlphaQt");
			histoTrueSecGammaQtPt 				= (TH2F*)directoryConv->Get("histoTrueSecGammaQtPt");
			histoTrueSecGammaPsiPairPt 			= (TH2F*)directoryConv->Get("histoTrueSecGammaPsiPairPt");
			histoTrueSecGammaCosPointPt 		= (TH2F*)directoryConv->Get("histoTrueSecGammaCosPointPt");
			histoTrueDalitzGammaChi2NDFPt 		= (TH2F*)directoryConv->Get("histoTrueDalitzGammaChi2NDFPt");
			histoTrueDalitzGammaAlphaQt 		= (TH2F*)directoryConv->Get("histoTrueDalitzGammaAlphaQt");
			histoTrueDalitzGammaQtPt 			= (TH2F*)directoryConv->Get("histoTrueDalitzGammaQtPt");
			histoTrueDalitzGammaPsiPairPt 		= (TH2F*)directoryConv->Get("histoTrueDalitzGammaPsiPairPt");
			histoTrueDalitzGammaCosPointPt 		= (TH2F*)directoryConv->Get("histoTrueDalitzGammaCosPointPt");
			histoTruePrimGammaAsymP 			= (TH2F*)directoryConv->Get("histoTruePrimGammaAsymP");
			histoTrueSecGammaAsymP 				= (TH2F*)directoryConv->Get("histoTrueSecGammaAsymP");
			histoTrueDalitzGammaAsymP 			= (TH2F*)directoryConv->Get("histoTrueDalitzGammaAsymP");
			histoTrueGammaChi2PsiPair 			= (TH2F*)directoryConv->Get("histoTrueGammaChi2PsiPair");
			histoTrueGammaCosPointChi2 			= (TH2F*)directoryConv->Get("histoTrueGammaCosPointChi2");
			histoTrueGammaCosPointPsiPair 		= (TH2F*)directoryConv->Get("histoTrueGammaCosPointPsiPair");
			histoTrueGammaInvMassPt 			= (TH2F*)directoryConv->Get("histoTrueGammaInvMassPt");
			histoTrueDalitzGammaChi2PsiPair 	= (TH2F*)directoryConv->Get("histoTrueDalitzGammaChi2PsiPair");
			histoTrueDalitzGammaCosPointChi2 	= (TH2F*)directoryConv->Get("histoTrueDalitzGammaCosPointChi2");
			histoTrueDalitzGammaCosPointPsiPair = (TH2F*)directoryConv->Get("histoTrueDalitzGammaCosPointPsiPair");
			histoTrueDalitzGammaInvMassPt 		= (TH2F*)directoryConv->Get("histoTrueDalitzGammaInvMassPt");
		}   

		// store sum of square of weights for the histograms
		histoElectrondEdxEtaP->Sumw2();
		histoPositrondEdxEtaP->Sumw2();
		histoElectronNSigmadEdxEtaP->Sumw2();
		histoPositronNSigmadEdxEtaP->Sumw2();
		histoElectronTOFEtaP->Sumw2();
		histoPositronTOFEtaP->Sumw2();
		histoElectronNSigmaTOFEtaP->Sumw2();
		histoPositronNSigmaTOFEtaP->Sumw2();
		histoElectronNSigmadEdxPhi->Sumw2();
		histoPositronNSigmadEdxPhi->Sumw2();
		histoElectronNSigmadEdxPhiEtaNeg->Sumw2();
		histoPositronNSigmadEdxPhiEtaNeg->Sumw2();
		histoElectronNSigmadEdxPhiEtaPos->Sumw2();
		histoPositronNSigmadEdxPhiEtaPos->Sumw2();
		histoElectronFClPt->Sumw2();
		histoPositronFClPt->Sumw2();
		histoElectronClPt->Sumw2();
		histoPositronClPt->Sumw2();
		histoGammaChi2NDFPt->Sumw2();
		histoGammaChi2NDFPtEtaNeg->Sumw2();
		histoGammaChi2NDFPtEtaPos->Sumw2();
		histoGammaChi2NDFR->Sumw2();
		histoElectronEtaPt->Sumw2();
		histoPositronEtaPt->Sumw2();
		histoElectronITSClR->Sumw2();
		histoPositronITSClR->Sumw2();
		histoElectronITSClPhi->Sumw2();		
		histoPositronITSClPhi->Sumw2();		
		histoElectronClR->Sumw2(); 			
		histoPositronClR->Sumw2();			
		histoElectronClPhi->Sumw2(); 		
		histoPositronClPhi->Sumw2();		
		histoGammaPhiEtaNeg->Sumw2();		
		histoGammaPhiEtaPos->Sumw2();		
		histoPositronNSigmadEdxCut->Sumw2();
		histoElectronNSigmadEdxCut->Sumw2();
		histoGammaEtaPt->Sumw2();
		histoGammaEtaR->Sumw2();
		histoGammaAlphaQt->Sumw2();
		histoGammaAlphaR->Sumw2();
		histoGammaQtPt->Sumw2();
		histoGammaQtR->Sumw2();
		histoGammaPhi->Sumw2();
		histoGammaPhiR->Sumw2();
		histoGammaPsiPairPt->Sumw2();
		histoGammaPsiPairPtEtaNeg->Sumw2();
		histoGammaPsiPairPtEtaPos->Sumw2();
		histoGammaPsiPairR->Sumw2();
		histoGammaCosPointPt->Sumw2();
		histoGammaCosPointR->Sumw2();
		histoGammaAsymP->Sumw2();
		histoGammaAsymR->Sumw2();
		histoGammaChi2PsiPair->Sumw2();
		histoGammaCosPointChi2->Sumw2();
		histoGammaCosPointPsiPair->Sumw2();
		histoGammaInvMassPt->Sumw2();
		histoGammaInvMassR->Sumw2();
		histoGammaZR->Sumw2();
		histoGammaXY->Sumw2();
		//ITS add
		histoElectronNSigmaITSEtaP->Sumw2();
		histoPositronNSigmaITSEtaP->Sumw2();
		histoElectronITSdEdxEtaP->Sumw2();
		histoPositronITSdEdxEtaP->Sumw2();
		histonSigdEdxElnSigdEdxPosGammaP->Sumw2();  // new sept 7 2015
		histoElectronITSClPt->Sumw2();
		histoPositronITSClPt->Sumw2();
		histoElectronITSClEta->Sumw2();
		histoPositronITSClEta->Sumw2();
		
		if (kMC){
			histoTrueElectronNSigmadEdxITSP->Sumw2();
			histoTrueElectronNSigmadEdxTPCP->Sumw2();
			histoTrueElectronNSigmaTOFP->Sumw2();
			histoTrueElectronITSClR->Sumw2();
			histoTruePositronNSigmadEdxITSP->Sumw2();
			histoTruePositronNSigmadEdxTPCP->Sumw2();
			histoTruePositronNSigmaTOFP->Sumw2();
			histoTruePositronITSClR->Sumw2();
			histoTruePrimGammaChi2NDFPt->Sumw2();
			histoTruePrimGammaAlphaQt->Sumw2();
			histoTruePrimGammaQtPt->Sumw2();
			histoTruePrimGammaCosPointPt->Sumw2();
			histoTruePrimGammaPsiPairPt->Sumw2();
			histoTrueSecGammaChi2NDFPt->Sumw2();
			histoTrueSecGammaAlphaQt->Sumw2();
			histoTrueSecGammaQtPt->Sumw2();
			histoTrueSecGammaCosPointPt->Sumw2();
			histoTrueSecGammaPsiPairPt->Sumw2();
			histoTrueDalitzGammaChi2NDFPt->Sumw2();
			histoTrueDalitzGammaAlphaQt->Sumw2();
			histoTrueDalitzGammaQtPt->Sumw2();
			histoTrueDalitzGammaCosPointPt->Sumw2();
			histoTrueDalitzGammaPsiPairPt->Sumw2();
			histoTrueDalitzGammaAsymP->Sumw2();
			histoTruePrimGammaAsymP->Sumw2();
			histoTrueSecGammaAsymP->Sumw2();
			histoTrueGammaChi2PsiPair->Sumw2();
			histoTrueGammaCosPointChi2->Sumw2();
			histoTrueGammaCosPointPsiPair->Sumw2();
			histoTrueGammaInvMassPt->Sumw2();
			histoTrueDalitzGammaChi2PsiPair->Sumw2();
			histoTrueDalitzGammaCosPointChi2->Sumw2();
			histoTrueDalitzGammaCosPointPsiPair->Sumw2();
			histoTrueDalitzGammaInvMassPt->Sumw2();
			
		}   
	}
	

	//********************************************************************************
	//*      Reading of Tree with reconstructed gammas/ filling of histograms        *
	//********************************************************************************
	
	Long64_t nEntriesRecGam 				= PhotonQA->GetEntries();
	Int_t nGammas 							= 0;
	cout << "Number of Gammas to be processed: " << nEntriesRecGam << endl;
	//    (*ptGammaAssosiatedPi0)[k]
	Long64_t nbytesRecGam 					= 0;
	Int_t negTracksWithoutTOF 				= 0;
	Int_t posTracksWithoutTOF 				= 0;
	Int_t negTracksWithoutnSigmaTOF 		= 0;
	Int_t posTracksWithoutnSigmaTOF 		= 0;
	Int_t negTracksWithoutITS 				= 0;
	Int_t posTracksWithoutITS 				= 0;
	Int_t negTracksWithoutnSigmaITS	 		= 0;
	Int_t posTracksWithoutnSigmaITS 		= 0;

	for (Long64_t i=0; i<nEntriesRecGam;i++) {
		nbytesRecGam 						+= PhotonQA->GetEntry(i); 
		Float_t photonEta 					= -TMath::Log(TMath::Tan(photonTheta/2.)); 
		Float_t photonP 					= photonPt * TMath::CosH(photonEta);
		Float_t photonX 					= (*gammaConvCoord)[0];
		Float_t photonY 					= (*gammaConvCoord)[1];
		Float_t photonZ 					= (*gammaConvCoord)[2];
		Float_t photonPhi 					= (*gammaConvCoord)[4];
		Float_t photonR 					= (*gammaConvCoord)[3];
		Float_t photonQt 					= (*gammaPhotonProp)[0];
		Float_t photonAlpha 				= (*gammaPhotonProp)[1];
		Float_t photonPsiPair 				= (*gammaPhotonProp)[2];
		Float_t photonCosPoint 				= (*gammaPhotonProp)[3];
		Float_t photonInvMass 				= (*gammaPhotonProp)[4];
		Float_t ptPositron 					= (*daugtherProp)[0];
		Float_t thetaPositron 				= (*daugtherProp)[1];
		Float_t etaPositron 				= -TMath::Log(TMath::Tan(thetaPositron/2.));
		Float_t pPositron 					= ptPositron * TMath::CosH(etaPositron);
		Float_t dEdxPositronTPC 			= (*daugtherProp)[2];						// <---- TPC instead of IROC
		Float_t nSigmaTPCPositron 			= (*daugtherProp)[3];
		Float_t tofPositron 				= (*daugtherProp)[4];
		Float_t nSigmaTOFPositron 			= (*daugtherProp)[5];
		Float_t fracClsTPCPositron 			= (*daugtherProp)[6];
		Float_t clsITSPositron 				= (*daugtherProp)[14];						// <----- cluster ITS instead of dE/dx Middle
		Float_t dEdxPositronITS 			= (*daugtherProp)[16];						// <----- dE/dx ITS instead of dE/dx Long
		Float_t clsTPCPositron 				= (*daugtherProp)[18];
		Float_t nSigmaITSPositron 			= (*daugtherProp)[20];						// <----- added
		Float_t ptElectron 					= (*daugtherProp)[7];
		Float_t thetaElectron 				= (*daugtherProp)[8];
		Float_t etaElectron 				= -TMath::Log(TMath::Tan(thetaElectron/2.));
		Float_t pElectron 					= ptElectron * TMath::CosH(etaElectron);
		Float_t dEdxElectronTPC 			= (*daugtherProp)[9];						// <---- TPC instead of IROC
		Float_t nSigmaTPCElectron 			= (*daugtherProp)[10];
		Float_t tofElectron 				= (*daugtherProp)[11];
		Float_t nSigmaTOFElectron 			= (*daugtherProp)[12];
		Float_t fracClsTPCElectron 			= (*daugtherProp)[13];
		Float_t clsITSElectron 				= (*daugtherProp)[15];						// <----- cluster ITS instead of dE/dx Middle
		Float_t dEdxElectronITS 			= (*daugtherProp)[17];						// <----- dE/dx ITS instead of dE/dx Long
		Float_t clsTPCElectron 				= (*daugtherProp)[19];
		Float_t nSigmaITSElectron 			= (*daugtherProp)[21];						// <----- added
		Float_t photonAsymE 				= pElectron/photonP;
		Float_t photonAsymP 				= pPositron/photonP;
		Bool_t kPassedAsymE 				= 0;
		Bool_t kPassedAsymP 				= 0;
		if (fAsymmCut){
			if (photonP > 0){
				if (photonAsymE > 0.5 && photonAsymE< 0.94) {
					if (photonP > 0.175*TMath::Tan(TMath::Pi()/0.94* (photonAsymE-0.5))) kPassedAsymE = 1;
				} else if (photonAsymE< 0.5 && photonAsymE> 0.06){
					if (photonP > 0.175*1/TMath::Tan(TMath::Pi()/0.94* (photonAsymE-0.03))) kPassedAsymE = 1;
				} else if (photonAsymE== 0.5) {
					kPassedAsymE 			= 1;
				}   
			}   
			
			if (photonP > 0){
				if (photonAsymP > 0.5 && photonAsymP< 0.94) {
					if (photonP > 0.175*TMath::Tan(TMath::Pi()/0.94* (photonAsymP-0.5))) kPassedAsymP = 1;
				} else if (photonAsymP< 0.5 && photonAsymP> 0.06){
					if (photonP > 0.175*1/TMath::Tan(TMath::Pi()/0.94* (photonAsymP-0.03))) kPassedAsymP = 1;
				} else if (photonAsymP== 0.5) {
					kPassedAsymP 			= 1;
					kPassedAsymE 			= 0;
					
				}   
			}   
		} else {
			kPassedAsymE 					= kTRUE;
			kPassedAsymP 					= kTRUE;
		}	
			
		Bool_t kPassedQt 					= 0;
		if (f2DQt) {
			if (TMath::Power(photonAlpha/0.95,2)+TMath::Power(photonQt/qtCut,2) < 1 ) kPassedQt = kTRUE;
		} else {
			if (photonQt < qtCut ) kPassedQt = kTRUE;
		}	
			
		Bool_t kPassedCosPointing 			= 0;
		if (photonCosPoint > cosPointingCut) kPassedCosPointing = kTRUE;
	
		Bool_t kPassedPsiPair 				= 0;
		Bool_t kPassedChi2 					= 0;
		if (f2DChi2PsiPair){
			if (TMath::Abs(photonPsiPair) < -psiPairCut/chi2PerDOFCut*photonChi2Ndf + psiPairCut){
				kPassedPsiPair	 			= kTRUE;
				kPassedChi2 				= kTRUE;
			}	
		} else {
			if (photonChi2Ndf < chi2PerDOFCut) kPassedChi2 				= kTRUE;
			if (TMath::Abs(photonPsiPair) < psiPairCut) kPassedPsiPair 	= kTRUE;
		}	
		histonSigdEdxElnSigdEdxPosGammaP->Fill(nSigmaTPCElectron,nSigmaTPCPositron,photonPt);  // new sept 7 2015
		
		if ( photonPt > minPtPhotonCut ){ 
			if (TMath::Abs(photonEta) >= etaMinCut && TMath::Abs(photonEta) <= etaMaxCut){ 
				if ( kPassedQt && kPassedAsymE && kPassedAsymP && kPassedCosPointing
				&& kPassedPsiPair && kPassedChi2){
					nGammas++;
					histoPositrondEdxEtaP->Fill(dEdxPositronTPC,etaPositron,pPositron);
					if (dEdxPositronITS == 1000){                                          //ITS
						posTracksWithoutITS++;
					} else {
						histoPositronITSdEdxEtaP->Fill(dEdxPositronITS,etaPositron,pPositron);
						if(nSigmaITSPositron == 20){
							posTracksWithoutnSigmaITS++;                    
						} else {
							histoPositronNSigmaITSEtaP->Fill(nSigmaITSPositron,etaPositron,pPositron);
							if (kMC){
								if (kind == 0 || kind == 5 || kind == 4 || kind == 3){
								histoTruePositronNSigmadEdxITSP->Fill(pPositron,nSigmaITSPositron);
								}   
							}   
						}
					}
					
					histoPositronNSigmadEdxEtaP->Fill(nSigmaTPCPositron,etaPositron,pPositron);
					if(pPositron<0.26){histoPositronNSigmadEdxCut->Fill(nSigmaTPCPositron);}	

					histoPositronNSigmadEdxPhi->Fill(nSigmaTPCPositron,photonPhi);  //dedx photon phi
					if(etaPositron < 0.) {histoPositronNSigmadEdxPhiEtaNeg->Fill(nSigmaTPCPositron,photonPhi); }
					if(etaPositron > 0.) {histoPositronNSigmadEdxPhiEtaPos->Fill(nSigmaTPCPositron,photonPhi); }

					histoPositronEtaPt->Fill(etaPositron,ptPositron);
					histoPositronFClPt->Fill(fracClsTPCPositron,ptPositron);
					histoPositronClPt->Fill(clsTPCPositron,ptPositron);
					histoPositronITSClPt->Fill(clsITSPositron,ptPositron);     //ITS
					histoPositronITSClEta->Fill(clsITSPositron,etaPositron);   //ITS
					
					if (kMC){
						if (kind == 0 || kind == 5 || kind == 4 || kind == 3){
							histoTruePositronNSigmadEdxTPCP->Fill(pPositron,nSigmaTPCPositron);
						}   
					}   
					if (tofPositron == 20000){
						posTracksWithoutTOF++;
					} else {
						histoPositronTOFEtaP->Fill(tofPositron,etaPositron,pPositron);
						if (nSigmaTOFPositron == -20){
							posTracksWithoutnSigmaTOF++;
						} else {
							histoPositronNSigmaTOFEtaP->Fill(nSigmaTOFPositron,etaPositron,pPositron);
							if (kMC){
								if (kind == 0 || kind == 5 || kind == 4 || kind == 3){
								histoTruePositronNSigmaTOFP->Fill(pPositron,nSigmaTOFPositron);
								}   
							}
						}
					}
					
					histoElectrondEdxEtaP->Fill(dEdxElectronTPC,etaElectron,pElectron);
					if (dEdxElectronITS == 1000){                                          //ITS
						negTracksWithoutITS++;
					} else {
						histoElectronITSdEdxEtaP->Fill(dEdxElectronITS,etaElectron,pElectron);
						if(nSigmaITSElectron == 20){
							negTracksWithoutnSigmaITS++; 
						} else {
							histoElectronNSigmaITSEtaP->Fill(nSigmaITSElectron,etaElectron,pElectron);
							if (kMC){
								if (kind == 0 || kind == 5 || kind == 4 || kind == 3){
								histoTrueElectronNSigmadEdxITSP->Fill(pElectron,nSigmaITSElectron);
								}   
							}
						}
					}
					
					histoElectronNSigmadEdxEtaP->Fill(nSigmaTPCElectron,etaElectron,pElectron);
					if (kMC){
						if (kind == 0 || kind == 5 || kind == 4 || kind == 3){
							histoTrueElectronNSigmadEdxTPCP->Fill(pElectron,nSigmaTPCElectron);
						}   
					}
					histoElectronEtaPt->Fill(etaElectron,ptElectron);
					if(pElectron<0.26){histoElectronNSigmadEdxCut->Fill(nSigmaTPCElectron);}	

					histoElectronNSigmadEdxPhi->Fill(nSigmaTPCElectron,photonPhi);    //dedx photon phi
					if(etaElectron < 0.) { histoElectronNSigmadEdxPhiEtaNeg->Fill(nSigmaTPCElectron,photonPhi);}   
					if(etaElectron > 0.) { histoElectronNSigmadEdxPhiEtaPos->Fill(nSigmaTPCElectron,photonPhi);}   //dedx photon phi

					histoElectronFClPt->Fill(fracClsTPCElectron,ptElectron);
					histoElectronClPt->Fill(clsTPCElectron,ptElectron);
					histoElectronITSClPt->Fill(clsITSElectron,ptElectron);     //ITS
					histoElectronITSClEta->Fill(clsITSElectron,etaElectron);   //ITS
					if (tofElectron == 20000){
						negTracksWithoutTOF++;
					} else {
						histoElectronTOFEtaP->Fill(tofElectron,etaElectron,pElectron);
						if (nSigmaTOFElectron == -20){
							negTracksWithoutnSigmaTOF++;
						} else {
							histoElectronNSigmaTOFEtaP->Fill(nSigmaTOFElectron,etaElectron,pElectron);
							if (kMC){
								if (kind == 0 || kind == 5 || kind == 4 || kind == 3){
								histoTrueElectronNSigmaTOFP->Fill(pElectron,nSigmaTOFElectron);
								}   
							}
						}
					}

					histoElectronITSClR->Fill(clsITSElectron,photonR);
					histoPositronITSClR->Fill(clsITSPositron,photonR);
					histoElectronITSClPhi->Fill(clsITSElectron,photonPhi); 		// new ITS-Phi
					histoPositronITSClPhi->Fill(clsITSPositron,photonPhi); 		// new ITS-Phi
					histoElectronClR->Fill(clsTPCElectron,photonR); 			// new TPC-R
					histoPositronClR->Fill(clsTPCPositron,photonR); 			// new TPC-R
					histoElectronClPhi->Fill(clsTPCElectron,photonPhi); 		// new TPC-Phi
					histoPositronClPhi->Fill(clsTPCPositron,photonPhi); 		// new TPC-Phi
					if(photonEta < 0.) {
						histoGammaPhiEtaNeg->Fill(photonPhi);
						histoGammaPsiPairPtEtaNeg->Fill(photonPsiPair, photonPt);
						histoGammaChi2NDFPtEtaNeg->Fill(photonChi2Ndf,photonPt);

					} 	// new EtaNeg
					if(photonEta > 0.) {
						histoGammaPhiEtaPos->Fill(photonPhi);
						histoGammaPsiPairPtEtaPos->Fill(photonPsiPair, photonPt);
						histoGammaChi2NDFPtEtaPos->Fill(photonChi2Ndf,photonPt);
					} 	// new EtaPos

					if (kMC){
						if (kind == 0 || kind == 5 || kind == 4 || kind == 3){
							histoTrueElectronITSClR->Fill(clsITSElectron,photonR);
							histoTruePositronITSClR->Fill(clsITSPositron,photonR);
						}   
					}
					histoGammaAlphaQt->Fill(photonAlpha,photonQt);
					histoGammaAlphaR->Fill(photonAlpha,photonR);					
					histoGammaQtPt->Fill(photonQt,photonPt);
					histoGammaQtR->Fill(photonQt,photonR);
					histoGammaPhi->Fill(photonPhi);
					histoGammaPhiR->Fill(photonPhi,photonR);
					histoGammaChi2NDFPt->Fill(photonChi2Ndf,photonPt);
					histoGammaChi2NDFR->Fill(photonChi2Ndf,photonR);
					histoGammaEtaPt->Fill(photonEta, photonPt);
					histoGammaEtaR->Fill(photonEta, photonR);
					histoGammaPsiPairPt->Fill(photonPsiPair, photonPt);
					histoGammaPsiPairR->Fill(photonPsiPair, photonR);
					histoGammaCosPointPt->Fill(photonCosPoint, photonPt);
					histoGammaCosPointR->Fill(photonCosPoint, photonR);
					histoGammaChi2PsiPair->Fill(photonChi2Ndf,photonPsiPair);
					histoGammaCosPointChi2->Fill(photonCosPoint, photonChi2Ndf);
					histoGammaCosPointPsiPair->Fill(photonCosPoint,photonPsiPair);
					histoGammaInvMassPt->Fill(photonInvMass, photonPt);
					histoGammaInvMassR->Fill(photonInvMass, photonR);
					if (photonP != 0){
						histoGammaAsymP->Fill(photonAsymE,photonP);
						histoGammaAsymR->Fill(photonAsymE,photonR);
					}
					histoGammaZR->Fill(photonZ, photonR);
					histoGammaXY->Fill(photonX, photonY);
					
					if (kMC){
						if (kind ==0){
							histoTruePrimGammaChi2NDFPt->Fill(photonChi2Ndf,photonPt);
							histoTruePrimGammaAlphaQt->Fill(photonAlpha,photonQt);
							histoTruePrimGammaQtPt->Fill(photonQt,photonPt);
							histoTruePrimGammaPsiPairPt->Fill(photonPsiPair, photonPt);
							histoTruePrimGammaCosPointPt->Fill(photonCosPoint, photonPt);
							if (photonP != 0){
								histoTruePrimGammaAsymP->Fill(photonAsymE,photonP);
							}
							histoTrueGammaChi2PsiPair->Fill(photonChi2Ndf,photonPsiPair);
							histoTrueGammaCosPointChi2->Fill(photonCosPoint, photonChi2Ndf);
							histoTrueGammaCosPointPsiPair->Fill(photonCosPoint,photonPsiPair);
							histoTrueGammaInvMassPt->Fill(photonInvMass, photonPt);
						}
						if (kind ==5){
							histoTrueSecGammaChi2NDFPt->Fill(photonChi2Ndf,photonPt);
							histoTrueSecGammaAlphaQt->Fill(photonAlpha,photonQt);
							histoTrueSecGammaQtPt->Fill(photonQt,photonPt);
							histoTrueSecGammaPsiPairPt->Fill(photonPsiPair, photonPt);
							histoTrueSecGammaCosPointPt->Fill(photonCosPoint, photonPt);
							if (photonP != 0){
								histoTrueSecGammaAsymP->Fill(photonAsymE,photonP);
							}
							histoTrueGammaChi2PsiPair->Fill(photonChi2Ndf,photonPsiPair);
							histoTrueGammaCosPointChi2->Fill(photonCosPoint, photonChi2Ndf);
							histoTrueGammaCosPointPsiPair->Fill(photonCosPoint,photonPsiPair);
							histoTrueGammaInvMassPt->Fill(photonInvMass, photonPt);
						}
						if (kind ==3 || kind == 4){
							histoTrueDalitzGammaChi2NDFPt->Fill(photonChi2Ndf,photonPt);
							histoTrueDalitzGammaAlphaQt->Fill(photonAlpha,photonQt);
							histoTrueDalitzGammaQtPt->Fill(photonQt,photonPt);
							histoTrueDalitzGammaPsiPairPt->Fill(photonPsiPair, photonPt);
							histoTrueDalitzGammaCosPointPt->Fill(photonCosPoint, photonPt);
							if (photonP != 0){
								histoTrueDalitzGammaAsymP->Fill(photonAsymE,photonP);
							}
							histoTrueDalitzGammaChi2PsiPair->Fill(photonChi2Ndf,photonPsiPair);
							histoTrueDalitzGammaCosPointChi2->Fill(photonCosPoint, photonChi2Ndf);
							histoTrueDalitzGammaCosPointPsiPair->Fill(photonCosPoint,photonPsiPair);
							histoTrueDalitzGammaInvMassPt->Fill(photonInvMass, photonPt);
						}
					}
				}
			}
		}
	}
	Double_t fracTOFsucNeg = (Double_t)negTracksWithoutTOF/nGammas*100. ;
	Double_t fracTOFsucPos = (Double_t)posTracksWithoutTOF/nGammas*100. ;
	Double_t fracITSsucNeg = (Double_t)negTracksWithoutITS/nGammas*100. ;
	Double_t fracITSsucPos = (Double_t)posTracksWithoutITS/nGammas*100. ;
	cout << "**************** TOF *****************" << endl;
	cout << nGammas  << "\t" << negTracksWithoutTOF << "\t" << negTracksWithoutnSigmaTOF << "\t"<< fracTOFsucNeg << "\t" << posTracksWithoutnSigmaTOF << "\t" << posTracksWithoutTOF<< "\t" << fracTOFsucPos << endl;
	cout << "**************** ITS *****************" << endl;
	cout << nGammas  << "\t" << negTracksWithoutITS<< "\t" << negTracksWithoutnSigmaITS << "\t"<< fracITSsucNeg << "\t" << posTracksWithoutnSigmaITS << "\t" << posTracksWithoutITS<< "\t" << fracITSsucPos << endl;
		
		
	//********************************************************************************
	//*                  Writing histograms to outputfile                            *
	//********************************************************************************     
	TFile* filePhotoQAWrite = new TFile(fileNameOutput.Data(),"UPDATE");
	if(addExtFolderToOutput){
		TDirectory *mainDir = filePhotoQAWrite->mkdir(Form("GammaConvV1_QA_%s", fCutSelection.Data()));
		mainDir->cd();
	}else{
		filePhotoQAWrite->mkdir(Form("GammaConv_%s",  specificCutSelection.Data()));
		filePhotoQAWrite->cd(Form("GammaConv_%s",  specificCutSelection.Data()));
	}
		histoElectrondEdxEtaP->Write("histoElectrondEdxEtaP",TObject::kWriteDelete);
		histoPositrondEdxEtaP->Write("histoPositrondEdxEtaP",TObject::kWriteDelete);
		
		histoElectronNSigmadEdxEtaP->Write("histoElectronNSigmadEdxEtaP",TObject::kWriteDelete);
		histoPositronNSigmadEdxEtaP->Write("histoPositronNSigmadEdxEtaP",TObject::kWriteDelete);
		histoElectronTOFEtaP->Write("histoElectronTOFEtaP",TObject::kWriteDelete);
		histoPositronTOFEtaP->Write("histoPositronTOFEtaP",TObject::kWriteDelete);
		histoElectronNSigmaTOFEtaP->Write("histoElectronNSigmaTOFEtaP",TObject::kWriteDelete);
		histoPositronNSigmaTOFEtaP->Write("histoPositronNSigmaTOFEtaP",TObject::kWriteDelete);
		histoElectronNSigmadEdxPhi->Write("histoElectronNSigmadEdxPhi",TObject::kWriteDelete);
		histoPositronNSigmadEdxPhi->Write("histoPositronNSigmadEdxPhi",TObject::kWriteDelete);
		histoElectronNSigmadEdxPhiEtaNeg->Write("histoElectronNSigmadEdxPhiEtaNeg",TObject::kWriteDelete);
		histoPositronNSigmadEdxPhiEtaNeg->Write("histoPositronNSigmadEdxPhiEtaNeg",TObject::kWriteDelete);
		histoElectronNSigmadEdxPhiEtaPos->Write("histoElectronNSigmadEdxPhiEtaPos",TObject::kWriteDelete);
		histoPositronNSigmadEdxPhiEtaPos->Write("histoPositronNSigmadEdxPhiEtaPos",TObject::kWriteDelete);
		histoElectronFClPt->Write("histoElectronFClPt",TObject::kWriteDelete);
		histoPositronFClPt->Write("histoPositronFClPt",TObject::kWriteDelete);
		histoElectronClPt->Write("histoElectronClPt",TObject::kWriteDelete);
		histoPositronClPt->Write("histoPositronClPt",TObject::kWriteDelete);
		histoGammaChi2NDFPt->Write("histoGammaChi2NDFPt",TObject::kWriteDelete);
		histoGammaChi2NDFPtEtaNeg->Write("histoGammaChi2NDFPtEtaNeg",TObject::kWriteDelete);
		histoGammaChi2NDFPtEtaPos->Write("histoGammaChi2NDFPtEtaPos",TObject::kWriteDelete);
		histoGammaChi2NDFR->Write("histoGammaChi2NDFR",TObject::kWriteDelete);
		histoElectronEtaPt->Write("histoElectronEtaPt",TObject::kWriteDelete);
		histoPositronEtaPt->Write("histoPositronEtaPt",TObject::kWriteDelete);
		histoElectronITSClR->Write("histoElectronITSClR",TObject::kWriteDelete);
		histoPositronITSClR->Write("histoPositronITSClR",TObject::kWriteDelete);
		histoElectronITSClPhi->Write("histoElectronITSClPhi",TObject::kWriteDelete);		
		histoPositronITSClPhi->Write("histoPositronITSClPhi",TObject::kWriteDelete);		
		histoElectronClR->Write("histoElectronClR",TObject::kWriteDelete);					
		histoPositronClR->Write("histoPositronClR",TObject::kWriteDelete);					
		histoElectronClPhi->Write("histoElectronClPhi",TObject::kWriteDelete);				
		histoPositronClPhi->Write("histoPositronClPhi",TObject::kWriteDelete);				
		histoGammaPhiEtaNeg->Write("histoGammaPhiEtaNeg",TObject::kWriteDelete);			
		histoGammaPhiEtaPos->Write("histoGammaPhiEtaPos",TObject::kWriteDelete);			
		histoPositronNSigmadEdxCut->Write("histoPositronNSigmadEdxCut",TObject::kWriteDelete);		
		histoElectronNSigmadEdxCut->Write("histoElectronNSigmadEdxCut",TObject::kWriteDelete);		
		if (ESDQADirectory)histoVertexZ->Write("histoVertexZ",TObject::kWriteDelete);				
		if (ESDQADirectory)histoContrVertexZ->Write("histoContrVertexZ",TObject::kWriteDelete);	
		if (ESDQADirectory)histoGoodESDTracks->Write("histoGoodESDTracks",TObject::kWriteDelete);		
		histoGammaEtaPt->Write("histoGammaEtaPt",TObject::kWriteDelete);
		histoGammaEtaR->Write("histoGammaEtaR",TObject::kWriteDelete);
		histoGammaAlphaQt->Write("histoGammaAlphaQt",TObject::kWriteDelete);
		histoGammaAlphaR->Write("histoGammaAlphaR",TObject::kWriteDelete);
		histoGammaQtPt->Write("histoGammaQtPt",TObject::kWriteDelete);
		histoGammaQtR->Write("histoGammaQtR",TObject::kWriteDelete);
		histoGammaPhi->Write("histoGammaPhi",TObject::kWriteDelete);
		histoGammaPhiR->Write("histoGammaPhiR",TObject::kWriteDelete);
		histoGammaPsiPairPt->Write("histoGammaPsiPairPt",TObject::kWriteDelete);
		histoGammaPsiPairPtEtaNeg->Write("histoGammaPsiPairPtEtaNeg",TObject::kWriteDelete);
		histoGammaPsiPairPtEtaPos->Write("histoGammaPsiPairPtEtaPos",TObject::kWriteDelete);
		histoGammaPsiPairR->Write("histoGammaPsiPairR",TObject::kWriteDelete);
		histoGammaCosPointPt->Write("histoGammaCosPointPt",TObject::kWriteDelete);
		histoGammaCosPointR->Write("histoGammaCosPointR",TObject::kWriteDelete);
		histoGammaAsymP->Write("histoGammaAsymP",TObject::kWriteDelete);
		histoGammaAsymR->Write("histoGammaAsymR",TObject::kWriteDelete);
		histoGammaChi2PsiPair->Write("histoGammaChi2PsiPair",TObject::kWriteDelete);
		histoGammaCosPointChi2->Write("histoGammaCosPointChi2",TObject::kWriteDelete);
		histoGammaCosPointPsiPair->Write("histoGammaCosPointPsiPair",TObject::kWriteDelete);
		histoGammaInvMassPt->Write("histoGammaInvMassPt",TObject::kWriteDelete);
		histoGammaInvMassR->Write("histoGammaInvMassR",TObject::kWriteDelete);
		histoGammaZR->Write("histoGammaZR",TObject::kWriteDelete);
		histoGammaXY->Write("histoGammaXY",TObject::kWriteDelete);
		//ITS add
		histoElectronNSigmaITSEtaP->Write("histoElectronNSigmaITSEtaP",TObject::kWriteDelete);
		histoPositronNSigmaITSEtaP->Write("histoPositronNSigmaITSEtaP",TObject::kWriteDelete);
		histoElectronITSdEdxEtaP->Write("histoElectronITSdEdxEtaP",TObject::kWriteDelete);
		histoPositronITSdEdxEtaP->Write("histoPositronITSdEdxEtaP",TObject::kWriteDelete);
		histonSigdEdxElnSigdEdxPosGammaP->Write("histonSigdEdxElnSigdEdxPosGammaP",TObject::kWriteDelete);
		histoElectronITSClPt->Write("histoElectronITSClPt",TObject::kWriteDelete);
		histoPositronITSClPt->Write("histoPositronITSClPt",TObject::kWriteDelete);
		histoElectronITSClEta->Write("histoElectronITSClEta",TObject::kWriteDelete);
		histoPositronITSClEta->Write("histoPositronITSClEta",TObject::kWriteDelete);
		
		if (kMC){
			histoTrueElectronNSigmadEdxITSP->Write("histoTrueElectronNSigmadEdxITSP",TObject::kWriteDelete);
			histoTrueElectronNSigmadEdxTPCP->Write("histoTrueElectronNSigmadEdxTPCP",TObject::kWriteDelete);
			histoTrueElectronNSigmaTOFP->Write("histoTrueElectronNSigmaTOFP",TObject::kWriteDelete);
			histoTrueElectronITSClR->Write("histoTrueElectronITSClR",TObject::kWriteDelete);
			histoTruePositronNSigmadEdxITSP->Write("histoTruePositronNSigmadEdxITSP",TObject::kWriteDelete);
			histoTruePositronNSigmadEdxTPCP->Write("histoTruePositronNSigmadEdxTPCP",TObject::kWriteDelete);
			histoTruePositronNSigmaTOFP->Write("histoTruePositronNSigmaTOFP",TObject::kWriteDelete);
			histoTruePositronITSClR->Write("histoTruePositronITSClR",TObject::kWriteDelete);
			histoTruePrimGammaChi2NDFPt->Write("histoTruePrimGammaChi2NDFPt",TObject::kWriteDelete);
			histoTruePrimGammaAlphaQt->Write("histoTruePrimGammaAlphaQt",TObject::kWriteDelete);
			histoTruePrimGammaQtPt->Write("histoTruePrimGammaQtPt",TObject::kWriteDelete);
			histoTruePrimGammaPsiPairPt->Write("histoTruePrimGammaPsiPairPt",TObject::kWriteDelete);
			histoTruePrimGammaCosPointPt->Write("histoTruePrimGammaCosPointPt",TObject::kWriteDelete);
			histoTruePrimGammaAsymP->Write("histoTruePrimGammaAsymP",TObject::kWriteDelete);
			histoTrueSecGammaChi2NDFPt->Write("histoTrueSecGammaChi2NDFPt",TObject::kWriteDelete);
			histoTrueSecGammaAlphaQt->Write("histoTrueSecGammaAlphaQt",TObject::kWriteDelete);
			histoTrueSecGammaQtPt->Write("histoTrueSecGammaQtPt",TObject::kWriteDelete);
			histoTrueSecGammaPsiPairPt->Write("histoTrueSecGammaPsiPairPt",TObject::kWriteDelete);
			histoTrueSecGammaCosPointPt->Write("histoTrueSecGammaCosPointPt",TObject::kWriteDelete);
			histoTrueSecGammaAsymP->Write("histoTrueSecGammaAsymP",TObject::kWriteDelete);
			histoTrueGammaChi2PsiPair->Write("histoTrueGammaChi2PsiPair",TObject::kWriteDelete);
			histoTrueGammaCosPointChi2->Write("histoTrueGammaCosPointChi2",TObject::kWriteDelete);
			histoTrueGammaCosPointPsiPair->Write("histoTrueGammaCosPointPsiPair",TObject::kWriteDelete);
			histoTrueGammaInvMassPt->Write("histoTrueGammaInvMassPt",TObject::kWriteDelete);
			histoTrueDalitzGammaChi2NDFPt->Write("histoTrueDalitzGammaChi2NDFPt",TObject::kWriteDelete);
			histoTrueDalitzGammaAlphaQt->Write("histoTrueDalitzGammaAlphaQt",TObject::kWriteDelete);
			histoTrueDalitzGammaQtPt->Write("histoTrueDalitzGammaQtPt",TObject::kWriteDelete);
			histoTrueDalitzGammaPsiPairPt->Write("histoTrueDalitzGammaPsiPairPt",TObject::kWriteDelete);
			histoTrueDalitzGammaCosPointPt->Write("histoTrueDalitzGammaCosPointPt",TObject::kWriteDelete);
			histoTrueDalitzGammaAsymP->Write("histoTrueDalitzGammaAsymP",TObject::kWriteDelete);
			histoTrueDalitzGammaChi2PsiPair->Write("histoTrueDalitzGammaChi2PsiPair",TObject::kWriteDelete);
			histoTrueDalitzGammaCosPointChi2->Write("histoTrueDalitzGammaCosPointChi2",TObject::kWriteDelete);
			histoTrueDalitzGammaCosPointPsiPair->Write("histoTrueDalitzGammaCosPointPsiPair",TObject::kWriteDelete);
			histoTrueDalitzGammaInvMassPt->Write("histoTrueDalitzGammaInvMassPt",TObject::kWriteDelete);
		}   
	filePhotoQAWrite->Write();
	filePhotoQAWrite->Close();
	delete filePhotoQAWrite;
	
	delete   histoElectrondEdxEtaP;
	delete   histoPositrondEdxEtaP;
	delete   histoElectronNSigmadEdxEtaP;
	delete   histoPositronNSigmadEdxEtaP;
	delete   histoElectronTOFEtaP;
	delete   histoPositronTOFEtaP;
	delete   histoElectronNSigmaTOFEtaP;
	delete   histoPositronNSigmaTOFEtaP;
	delete   histoElectronNSigmadEdxPhi;
	delete   histoPositronNSigmadEdxPhi;
	delete   histoElectronNSigmadEdxPhiEtaNeg;
	delete   histoPositronNSigmadEdxPhiEtaNeg;
	delete   histoElectronNSigmadEdxPhiEtaPos;
	delete   histoPositronNSigmadEdxPhiEtaPos;
	delete   histoElectronFClPt;
	delete   histoPositronFClPt;
	delete   histoElectronClPt;
	delete   histoPositronClPt;
	delete   histoGammaChi2NDFPt;
	delete   histoGammaChi2NDFPtEtaNeg;
	delete   histoGammaChi2NDFPtEtaPos;
	delete   histoGammaChi2NDFR;
	delete   histoElectronEtaPt;
	delete   histoPositronEtaPt;
	delete   histoElectronITSClR;
	delete   histoPositronITSClR;
	delete   histoElectronITSClPhi;		
	delete   histoPositronITSClPhi;		
	delete   histoElectronClR;			
	delete   histoPositronClR;			
	delete   histoElectronClPhi;		
	delete   histoPositronClPhi;		
	delete   histoGammaPhiEtaNeg;		
	delete   histoGammaPhiEtaPos;		
	delete	 histoPositronNSigmadEdxCut;	
	delete	 histoElectronNSigmadEdxCut;	
	if (histoVertexZ) 				delete	 histoVertexZ;				
	if (histoContrVertexZ) 			delete	 histoContrVertexZ;			
	if (histoVertexZdummy) 			delete	 histoVertexZdummy;			
	if (histoContrVertexZdummy)		delete	 histoContrVertexZdummy;
	if (histoGoodESDTracks) 		delete	 histoGoodESDTracks;			
	if (histoGoodESDTracksdummy)	delete	 histoGoodESDTracksdummy;	
	delete   histoGammaEtaPt;
	delete   histoGammaEtaR;
	delete   histoGammaAlphaQt;
	delete   histoGammaAlphaR;
	delete   histoGammaQtPt;
	delete   histoGammaQtR;
	delete   histoGammaPhi;
	delete   histoGammaPhiR;
	delete   histoGammaPsiPairPt;
	delete   histoGammaPsiPairPtEtaNeg;
	delete   histoGammaPsiPairPtEtaPos;
	delete   histoGammaPsiPairR;
	delete   histoGammaCosPointPt;
	delete   histoGammaCosPointR;
	delete   histoGammaAsymP;
	delete   histoGammaAsymR;
	delete   histoGammaChi2PsiPair;
	delete   histoGammaCosPointChi2;
	delete   histoGammaCosPointPsiPair;
	delete   histoGammaInvMassPt;
	delete   histoGammaInvMassR;
	delete   histoGammaZR;
	delete   histoGammaXY;
	delete 	 histoElectronNSigmaITSEtaP;
	delete   histoPositronNSigmaITSEtaP;
	delete   histoElectronITSdEdxEtaP;
	delete   histoPositronITSdEdxEtaP;
	delete   histonSigdEdxElnSigdEdxPosGammaP;
	delete   histoElectronITSClPt;
	delete   histoPositronITSClPt;
	delete   histoElectronITSClEta;
	delete   histoPositronITSClEta;
	delete   histoTrueElectronNSigmadEdxITSP;
	delete   histoTrueElectronNSigmadEdxTPCP;
	delete   histoTrueElectronNSigmaTOFP;
	delete   histoTrueElectronITSClR;
	delete   histoTruePositronNSigmadEdxITSP;
	delete   histoTruePositronNSigmadEdxTPCP;
	delete   histoTruePositronNSigmaTOFP;
	delete   histoTruePositronITSClR;
	delete   histoTrueDalitzGammaChi2NDFPt;
	delete   histoTrueDalitzGammaAlphaQt;
	delete   histoTrueDalitzGammaQtPt;
	delete   histoTrueDalitzGammaCosPointPt;
	delete   histoTrueDalitzGammaPsiPairPt;
	delete   histoTrueDalitzGammaAsymP;
	delete   histoTrueSecGammaChi2NDFPt;
	delete   histoTrueSecGammaAlphaQt;
	delete   histoTrueSecGammaQtPt;
	delete   histoTrueSecGammaCosPointPt;
	delete   histoTrueSecGammaPsiPairPt;
	delete   histoTrueSecGammaAsymP;
	delete   histoTruePrimGammaChi2NDFPt;
	delete   histoTruePrimGammaAlphaQt;
	delete   histoTruePrimGammaQtPt;
	delete   histoTruePrimGammaCosPointPt;
	delete   histoTruePrimGammaPsiPairPt;
	delete   histoTruePrimGammaAsymP;
	delete   histoTrueGammaChi2PsiPair;
	delete   histoTrueGammaCosPointChi2;
	delete   histoTrueGammaCosPointPsiPair;
	delete   histoTrueGammaInvMassPt;
	delete   histoTrueDalitzGammaChi2PsiPair;
	delete   histoTrueDalitzGammaCosPointChi2;
	delete   histoTrueDalitzGammaCosPointPsiPair;
	delete   histoTrueDalitzGammaInvMassPt;

	delete daugtherProp;
	delete gammaConvCoord;
	delete gammaPhotonProp;

	if(directoryConv){
		directoryConv->Clear();
		delete directoryConv;
	}

	if(fileMappingDetailedConv){
		fileMappingDetailedConv->Close();
		delete fileMappingDetailedConv;
	}

	if(treeList){
		treeList->Clear();
		delete treeList;
	}

	if(fGammaDir){
		fGammaDir->Clear();
		delete fGammaDir;
	}

	f->Close();
	delete f;
	
}
