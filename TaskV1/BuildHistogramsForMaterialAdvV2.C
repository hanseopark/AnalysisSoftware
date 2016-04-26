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
typedef TVectorT<float> TVectorF;

TGraph* grInnerVessel ;
TGraph* grInnerVesselMC ;
TGraph* grThermalShieldSDDMC ;
TGraph* grThermalShieldSPDMC ;
TGraph* grBPMC ;

void SetLogBinningTH1(TH1* histoRebin){
	TAxis *axisafter = histoRebin->GetXaxis(); 
	Int_t bins = axisafter->GetNbins();
	Double_t from = axisafter->GetXmin();
	Double_t to = axisafter->GetXmax();
	Double_t *newbins = new Double_t[bins+1];
	newbins[0] = from;
	Double_t factor = TMath::Power(to/from, 1./bins);
	for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
	axisafter->Set(bins, newbins);
	delete [] newbins;

}

void SetLogBinningXTH2(TH2* histoRebin){
	TAxis *axisafter = histoRebin->GetXaxis(); 
	Int_t bins = axisafter->GetNbins();
	Double_t from = axisafter->GetXmin();
	Double_t to = axisafter->GetXmax();
	Double_t *newbins = new Double_t[bins+1];
	newbins[0] = from;
	Double_t factor = TMath::Power(to/from, 1./bins);
	for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
	axisafter->Set(bins, newbins);
	delete [] newbins;

}

void SetLogBinningYTH2(TH2* histoRebin){
	TAxis *axisafter = histoRebin->GetYaxis(); 
	Int_t bins = axisafter->GetNbins();
	Double_t from = axisafter->GetXmin();
	Double_t to = axisafter->GetXmax();
	Double_t *newbins = new Double_t[bins+1];
	newbins[0] = from;
	Double_t factor = TMath::Power(to/from, 1./bins);
	for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
	axisafter->Set(bins, newbins);
	delete [] newbins;

}


//____________________________________________________________________
void circleFittingInnerVesselRec(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
	//minimisation function computing the sum of squares of residuals
	Int_t np = grInnerVessel->GetN();
	f = 0;
	Double_t *x = grInnerVessel->GetX();
	Double_t *y = grInnerVessel->GetY();
	for (Int_t i=0;i<np;i++) {
		Double_t u = x[i] - par[0];
		Double_t v = y[i] - par[1];
		Double_t dr = par[2] - TMath::Sqrt(u*u+v*v);
		f += dr*dr;
	}
}

//____________________________________________________________________
void circleFittingInnerVesselMC(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
	//minimisation function computing the sum of squares of residuals
	Int_t np = grInnerVesselMC->GetN();
	f = 0;
	Double_t *x = grInnerVesselMC->GetX();
	Double_t *y = grInnerVesselMC->GetY();
	for (Int_t i=0;i<np;i++) {
		Double_t u = x[i] - par[0];
		Double_t v = y[i] - par[1];
		Double_t dr = par[2] - TMath::Sqrt(u*u+v*v);
		f += dr*dr;
	}
}

//____________________________________________________________________
void circleFittingTSSDDMC(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
	//minimisation function computing the sum of squares of residuals
	Int_t np = grThermalShieldSDDMC->GetN();
	f = 0;
	Double_t *x = grThermalShieldSDDMC->GetX();
	Double_t *y = grThermalShieldSDDMC->GetY();
	for (Int_t i=0;i<np;i++) {
		Double_t u = x[i] - par[0];
		Double_t v = y[i] - par[1];
		Double_t dr = par[2] - TMath::Sqrt(u*u+v*v);
		f += dr*dr;
	}
}

//____________________________________________________________________
void circleFittingTSSPDMC(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
	//minimisation function computing the sum of squares of residuals
	Int_t np = grThermalShieldSPDMC->GetN();
	f = 0;
	Double_t *x = grThermalShieldSPDMC->GetX();
	Double_t *y = grThermalShieldSPDMC->GetY();
	for (Int_t i=0;i<np;i++) {
		Double_t u = x[i] - par[0];
		Double_t v = y[i] - par[1];
		Double_t dr = par[2] - TMath::Sqrt(u*u+v*v);
		f += dr*dr;
	}
}

//____________________________________________________________________
void circleFittingBPMC(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
   //minimisation function computing the sum of squares of residuals
	Int_t np = grBPMC->GetN();
	f = 0;
	Double_t *x = grBPMC->GetX();
	Double_t *y = grBPMC->GetY();
	for (Int_t i=0;i<np;i++) {
		Double_t u = x[i] - par[0];
		Double_t v = y[i] - par[1];
		Double_t dr = par[2] - TMath::Sqrt(u*u+v*v);
		f += dr*dr;
	}
}


void BuildHistogramsForMaterialAdvV2(TString fileName = "kkoch_GammaConvMaterialTree.root", TString gammaCutSelection ="sdsdd", TString CutSelection="000400" , Bool_t kMC=0, Bool_t merge=0, TString nameOutputBase = "MappingDetailedGammaConv"){
   //Reset ROOT and connect tree file
   
	//********************************************************************************
	//*            Definition of Cuts                                                *
	//********************************************************************************
   
	TString V0Finder = CutSelection(0,1);
	TString etaCutNumber = CutSelection(1,1);
	TString minPtCutNumber = CutSelection(2,1);
	TString chi2CutNumber = CutSelection(3,1);
	TString minMultCutNumber = CutSelection(4,1);
	TString maxMultCutNumber = CutSelection(5,1);
	
	if ( V0Finder.CompareTo("0") == 0){
		cout << "V0Finder: " <<  "Onfly" << endl;
	} else if ( V0Finder.CompareTo("1") == 0){
		cout << "V0Finder: " <<  "Offline" << endl;
	}
	Double_t fLineCutZVtx = 5.;
	Double_t fLineCutZVtxMin = 5.;
	Double_t etaMinCut = -0.1;
	Double_t etaMaxCut = 0.;
	if (etaCutNumber.CompareTo("0") == 0){
		etaMaxCut = 0.9;
	} else if (etaCutNumber.CompareTo("1") == 0 ){
		etaMaxCut = 0.1;
	}else if (etaCutNumber.CompareTo("2") == 0 ){
		etaMaxCut = 1.4;
		etaMinCut = 0.9;
		fLineCutZVtx = 10.;
		fLineCutZVtxMin = 5.;
	}else if (etaCutNumber.CompareTo("3") == 0 ){
		etaMaxCut = 1.8;
		etaMinCut = 1.4;
		fLineCutZVtx = 20.;
		fLineCutZVtxMin = 10.;
	} else if (etaCutNumber.CompareTo("4") == 0 ){
		etaMaxCut = 1.4;
		etaMinCut = -0.1;
	} else if (etaCutNumber.CompareTo("5") == 0 ){
		etaMaxCut = 10.;
		etaMinCut = 2.5;
		fLineCutZVtx = 20.;
		fLineCutZVtxMin = 20.;
	} else if (etaCutNumber.CompareTo("6") == 0 ){
		etaMaxCut = 10.;
		etaMinCut = -0.1;
	} else if (etaCutNumber.CompareTo("7") == 0 ){
		etaMaxCut = 0.65;
	}
	if (etaMinCut != -0.1){
		cout << "Eta range: " <<   etaMinCut  <<  " < |eta| < " << etaMaxCut << endl;
	} else {
		cout << "Eta range: " <<  "|eta| < " << etaMaxCut << endl;
	}
	Double_t fLineCutZRSlope = tan (2*atan(exp(-etaMaxCut)));
	Double_t fLineCutZRSlopeMin = tan (2*atan(exp(-etaMinCut)));
	
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
	} else if (minPtCutNumber.CompareTo("5") == 0 ){
		minPtPhotonCut = 0.5;
	} else if (minPtCutNumber.CompareTo("6") == 0 ){
		minPtPhotonCut = 0.75;	
	} else if (minPtCutNumber.CompareTo("7") == 0 ){
		minPtPhotonCut = 1.;	
		
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
		
	Double_t midPtPhotonMin = 1.;
	Double_t midPtPhotonMax = 3.;
	cout << "Mid-Pt-Range : "<< midPtPhotonMin << " - " << midPtPhotonMax << endl;

	Int_t minMult = 0;
	Int_t maxMult = 2000;
	if (minMultCutNumber.CompareTo("0") == 0){
		minMult = 0.;
	} else if (minMultCutNumber.CompareTo("1") == 0 ){
		minMult = 2.;
	} else if (minMultCutNumber.CompareTo("2") == 0 ){
		minMult = 5.;
	} else if (minMultCutNumber.CompareTo("3") == 0 ){
		minMult = 10.;
	} else if (minMultCutNumber.CompareTo("4") == 0 ){
		minMult = 15.;
	} else if (minMultCutNumber.CompareTo("5") == 0 ){
		minMult = 20.;
	} else if (minMultCutNumber.CompareTo("6") == 0 ){
		minMult = 30.;	  
	} else if (minMultCutNumber.CompareTo("7") == 0 ){
		minMult = 40.;	  
	} else if (minMultCutNumber.CompareTo("8") == 0 ){
		minMult = 50.;	  
	} else if (minMultCutNumber.CompareTo("9") == 0 ){
		minMult = 60.;	  
	}
	
	if (maxMultCutNumber.CompareTo("0") == 0){
		maxMult = 2000.;
	} else if (maxMultCutNumber.CompareTo("1") == 0 ){
		maxMult = 100.;
	} else if (maxMultCutNumber.CompareTo("2") == 0 ){
		maxMult = 60.;
	} else if (maxMultCutNumber.CompareTo("3") == 0 ){
		maxMult = 50.;
	} else if (maxMultCutNumber.CompareTo("4") == 0 ){
		maxMult = 40.;
	} else if (maxMultCutNumber.CompareTo("5") == 0 ){
		maxMult = 30.;
	} else if (maxMultCutNumber.CompareTo("6") == 0 ){
		maxMult = 20.;	  
	} else if (maxMultCutNumber.CompareTo("7") == 0 ){
		maxMult = 15.;	  
	} else if (maxMultCutNumber.CompareTo("8") == 0 ){
		maxMult = 10.;	  
	} else if (maxMultCutNumber.CompareTo("9") == 0 ){
		maxMult = 5.;	  
	}
	cout<< "Mult range Cut: " << minMult << "-" << maxMult << endl;
	
	//********************************************************************************
	//*      Setting of boundaries for different regions of the material             *
	//********************************************************************************

	Double_t rMin=2.;
	Double_t rSPD=9.5;
	Double_t rSPDth=13;
	Double_t rSDD=35;
	Double_t rSSD=55;
	Double_t rTPC=72.5;
	Double_t rHotZoneMin = 5.7;
	Double_t rHotZoneMax = 50.;
	Double_t zHotZoneMin = 45.;
	Double_t zHotZoneMax = 75.;
	Double_t zBeamPipeInner = 30.;   
	
	Double_t rDetailsSDD1Min = 12.;
	Double_t rDetailsSDD1Max = 20.;
	Double_t rDetailsSDD2Min = 20.;
	Double_t rDetailsSDD2Max = 28.;
	Double_t rDetailsSDDThermalShieldMin = 28.;
	Double_t rDetailsSDDThermalShieldMax = 32.;

	gROOT->Reset();

	//********************************************************************************
	//*            File definition/ loading                                          *
	//********************************************************************************
	
	TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fileName.Data());
	if (!f) {
		f = new TFile(fileName.Data());
	}

	//********************************************************************************
	//*            Event extraction                                                  *
	//********************************************************************************
	TList *listGlobal = (TList*)f->Get(Form("GammaConvMaterial_%s",gammaCutSelection.Data()));
	TList *listEvent = (TList*)listGlobal->FindObject("EventList");
	TTree *Event = (TTree*)listEvent->FindObject("Event");
	//Declaration of leaves types
	Float_t        primVtxZ;
	Int_t           nContrVtx;
	Int_t           fMultiplicity14;
	Int_t           fMultiplicity09;
	
	// Set branch addresses.
	Event->SetBranchAddress("primVtxZ",&primVtxZ);
	Event->SetBranchAddress("nContrVtx",&nContrVtx);
	Event->SetBranchAddress("nGoodTracksEta09",&fMultiplicity09);
	Event->SetBranchAddress("nGoodTracksEta14",&fMultiplicity14);

	Long64_t nEntriesEvent = Event->GetEntries();
	cout << "number of Events: " << nEntriesEvent<< endl;
	

	//********************************************************************************
	//*            Definition of Reconstructed Conversion Points                     *
	//********************************************************************************
	
	//Declaration of leaves types
	Float_t        fTheta;
	Float_t        fChiPerDOF;
	Float_t        fPt;
	//    Float_t        fP;
	UChar_t        fKind;
	TVectorF * fDaugtherProp = new TVectorF;
	TVectorF * fCoord = new TVectorF;
	Int_t         fMultiplicity09Photon;
	
	
	TList* listRecGamma = (TList*)listGlobal->FindObject("RecGammaList");
	TTree *ConvPointsRec = (TTree*)listRecGamma->FindObject("ConvPointRec");

	// Set branch addresses.
	ConvPointsRec->SetBranchAddress("recCords",&fCoord);
	ConvPointsRec->SetBranchAddress("daughterProp",&fDaugtherProp);
	ConvPointsRec->SetBranchAddress("theta",&fTheta);
	ConvPointsRec->SetBranchAddress("chi2ndf",&fChiPerDOF);
	ConvPointsRec->SetBranchAddress("pt",&fPt);
	ConvPointsRec->SetBranchAddress("nGoodTracksEta09",&fMultiplicity09Photon);
	if (kMC) ConvPointsRec->SetBranchAddress("kind",&fKind);
	
	//********************************************************************************
	//*            Definition of MC input Conversion Points                          *
	//********************************************************************************
	/*   
	//Declaration of leaves types
	Float_t        gTheta;
	Float_t        gPt;
	TVectorF * gDaugtherProp = new TVectorF;
	TVectorF * gCoord = new TVectorF;
	
	TList* listConvGammaMC = (TList*)listGlobal->FindObject("AllMCGammaConvList");
	TTree *ConvPointsMCInput =NULL;

	if ( kMC){
		ConvPointsMCInput = (TTree*)listConvGammaMC->FindObject("ConvGammaMC");
		// Set branch addresses.
		ConvPointsMCInput->SetBranchAddress("Cords",&gCoord);
		ConvPointsMCInput->SetBranchAddress("daughterProp",&gDaugtherProp);
		ConvPointsMCInput->SetBranchAddress("Theta",&gTheta);
		ConvPointsMCInput->SetBranchAddress("Pt",&gPt);
	}  */

	//********************************************************************************
	//*            Definition of MC input Photons                                    *
	//********************************************************************************
	
	//    //Declaration of leaves types
	//    Float_t        lTheta;
	//    Float_t        lPt;
	//    
	//    TList* listAllGammaMC = (TList*)listGlobal->FindObject("AllMCGammaList");
	//    TTree *PhotonMCInput = NULL;
	// 
	//    if ( kMC){
	//       PhotonMCInput = (TTree*)listAllGammaMC->FindObject("AllGamma");
	//    // Set branch addresses.
	//       PhotonMCInput->SetBranchAddress("theta",&lTheta);
	//       PhotonMCInput->SetBranchAddress("pt",&lPt);
	//    }  
	
	
	//********************************************************************************
	//*            Definition of Boundaries for Histograms                           *
	//********************************************************************************

	//_______________EventQuality-plot____________________
		Int_t nXBinsEvtQ        = 7;
		Double_t firstXBinEvtQ  = -0.5;
		Double_t lastXBinEvtQ   = 6.5;

	//_______________Z vtx dist____________________
		Int_t nXBinsZvtx        = 400;
		Double_t firstXBinZvtx  = -20;
		Double_t lastXBinZvtx   = 20;

		
	//_______________Number of ESD track-Plot_____________
		Int_t nBinsESDtrk = 150;
		Double_t firstBinESDtrk= -0.5;
		Double_t lastBinESDtrk = 149.5;


	//______________binning fine conversions _________________
		Int_t nMappingZConv        = 9000;
		Double_t firstBinMappingZConv    = -300.;
		Double_t lastBinMappingZConv  = 300.;
		Int_t nMappingRConv        = 2700;
		Double_t firstBinMappingRConv    = 0.;
		Double_t lastBinMappingRConv  = 180.;
		Int_t nMappingXConv        = 4000;
		Double_t firstBinMappingXConv    = -200.;
		Double_t lastBinMappingXConv  = 200.;
		Int_t nMappingYConv        = 4000;
		Double_t firstBinMappingYConv    = -200.;
		Double_t lastBinMappingYConv  = 200.;
		Int_t nMappingPhiConv      = 400;
		Double_t firstBinMappingPhiConv= 0;
		Double_t lastBinMappingPhiConv = 2*TMath::Pi();
		Int_t nMappingEtaConv         = 400;
		Double_t firstBinMappingEtaConv  = -2.;
		Double_t lastBinMappingEtaConv   = 2.;
		Int_t nMappingXConvHotZone          = 1000;
		Double_t firstBinMappingXConvHotZone   = -50.;
		Double_t lastBinMappingXConvHotZone    = 50.;
		Int_t nMappingYConvHotZone          = 1000;
		Double_t firstBinMappingYConvHotZone   = -50.;
		Double_t lastBinMappingYConvHotZone    = 50.;
		Int_t nMappingZConvHotZone          = 300;
		Double_t firstBinMappingZConvHotZone   = 45.;
		Double_t lastBinMappingZConvHotZone    = 75.;
		Int_t nMappingXBPConv         = 600;
		Double_t firstBinMappingXConvBP  = -15.;
		Double_t lastBinMappingXConvBP   = 15.;
		Int_t nMappingYBPConv         = 600;
		Double_t firstBinMappingYConvBP  = -15.;
		Double_t lastBinMappingYConvBP   = 15.;

		Int_t nMappingPhiSDD = 1200;
		Double_t firstBinMappingPhiSDD= 0;
		Double_t lastBinMappingPhiSDD = 2*TMath::Pi();
		
		Int_t nMappingRSDD1 = 160;
		Double_t firstBinMappingRSDD1    = 12.;
		Double_t lastBinMappingRSDD1  = 20.;
		Int_t nMappingRSDD2 = 160;
		Double_t firstBinMappingRSDD2    = 20.;
		Double_t lastBinMappingRSDD2  = 28.;
		Int_t nMappingRSDDthermal = 120;
		Double_t firstBinMappingRSDDthermal    = 28.;
		Double_t lastBinMappingRSDDthermal  = 32.;
		
		Int_t nMappingZSDD = 160;
		Double_t firstBinMappingZSDD    = -40.;
		Double_t lastBinMappingZSDD  = 40.;		
		
	//______________Quality plots______________________________ 
		Int_t nChi2Conv         = 200;
		Double_t firstBinChi2Conv  = 0;
		Double_t lastBinChi2Conv      = 50.;

		Int_t nPtConv        = 200;
		Double_t firstBinPtConv    = 0.01;
		Double_t lastBinPtConv     = 20.;

	//********************************************************************************
	//*      Definition of histograms for reconstructed Conversion Points            *
	//********************************************************************************
	TH1I* histoEventQual;
	TH1F* histoGoodESDtrk09;
	TH1F* histoGoodESDtrk14;
	TH1F* histoGoodESDtrk0914;
	TH1F* histoGoodESDtrk09Vtx;
	TH1F* histoGoodESDtrk14Vtx;
	TH1F* histoGoodESDtrk0914Vtx;
	TH1F* histoZVtxDist;
	
	TH2F* histoZR; 
	TH2F* histoXY; 
	TH2F* histoZPhi;
	TH2F* histoRPhi; 
	TH2F* histoZRMidPt; 
	TH2F* histoZPhiMidPt;
	TH2F* histoRPhiMidPt; 
	TH2F* histoZPhi_SPD;
	TH2F* histoZPhi_SPDTh;
	TH2F* histoZPhi_SSD;
	TH2F* histoZPhi_SDD;
	TH2F* histoZPhi_ITSTPC;
	TH2F* histoZPhi_HotZone;
	TH2F* histoZPhi_SDD1;
	TH2F* histoZPhi_SDD2;
	TH2F* histoZPhi_SDDThermal;
	TH2F* histoRPhi_SDD1;
	TH2F* histoRPhi_SDD2;
	TH2F* histoRPhi_SDDThermal;
	TH2F* histoZR_SDD1;
	TH2F* histoZR_SDD2;
	TH2F* histoZR_SDDThermal;
	TH2F* histoXY_HotZone;
	TH2F* histoXY_Beampipe;
	TH1F* histoEta;

	TH1F* histoChi2PerDOF;
	TH1F* histoPt;
	//       
	TH1F* histoChi2PerDOF_AftCuts;
	TH1F* histoPt_AftCuts;

	TString fileNameOutput;
	if (kMC){
		fileNameOutput = Form("%s_MC.root",nameOutputBase.Data()) ;
	} else {
		fileNameOutput = Form("%s_Data.root",nameOutputBase.Data()) ;
	}
	TFile* fileMappingDetailedConv = new TFile(fileNameOutput.Data());
	TDirectory*  directoryConv =     (TDirectory*)fileMappingDetailedConv->Get(Form("GammaConv_%s",  CutSelection.Data()));          
	if (!merge || directoryConv==0) {   
		histoEventQual = new TH1I("ESD_EventQuality","ESD_EventQuality",nXBinsEvtQ,firstXBinEvtQ,lastXBinEvtQ);
		histoEventQual->SetXTitle("EventQuality");
		histoZVtxDist = new TH1F("Z_Vertex_distribution","Z_Vertex_distribution",nXBinsZvtx,firstXBinZvtx,lastXBinZvtx);
		histoZVtxDist->SetXTitle("Z_{prim Vtx} (cm)");
		histoGoodESDtrk09 = new TH1F("ESD_NumberOfGoodESDTracks09","Number of Good ESD tracks",nBinsESDtrk, firstBinESDtrk, lastBinESDtrk);
		histoGoodESDtrk09->SetXTitle("# good track TPC, |#eta| < 0.9");
		histoGoodESDtrk09Vtx = new TH1F("ESD_NumberOfGoodESDTracks09Vtx","Number of Good ESD tracks",nBinsESDtrk, firstBinESDtrk, lastBinESDtrk);
		histoGoodESDtrk09Vtx->SetXTitle("# good track TPC, |#eta| < 0.9");
		histoGoodESDtrk14 = new TH1F("ESD_NumberOfGoodESDTracks14","Number of Good ESD tracks",nBinsESDtrk, firstBinESDtrk, lastBinESDtrk);
		histoGoodESDtrk14->SetXTitle("# good track TPC, |#eta| < 1.4");
		histoGoodESDtrk14Vtx = new TH1F("ESD_NumberOfGoodESDTracks14Vtx","Number of Good ESD tracks",nBinsESDtrk, firstBinESDtrk, lastBinESDtrk);
		histoGoodESDtrk14Vtx->SetXTitle("# good track TPC, |#eta| < 1.4");
		histoGoodESDtrk0914 = new TH1F("ESD_NumberOfGoodESDTracks0914","Number of Good ESD tracks",nBinsESDtrk, firstBinESDtrk, lastBinESDtrk);
		histoGoodESDtrk0914->SetXTitle("# good track TPC, 0.9 < |#eta| < 1.4");
		histoGoodESDtrk0914Vtx = new TH1F("ESD_NumberOfGoodESDTracks0914Vtx","Number of Good ESD tracks",nBinsESDtrk, firstBinESDtrk, lastBinESDtrk);
		histoGoodESDtrk0914Vtx->SetXTitle("# good track TPC, 0.9 < |#eta| < 1.4");
		
		histoZR =  new TH2F("ESD_ConversionMapping_ZR","", nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingRConv, firstBinMappingRConv, lastBinMappingRConv);
		histoZR->SetXTitle("Z (cm)");
		histoZR->SetYTitle("R (cm)");
		histoXY =  new TH2F("ESD_ConversionMapping_XY","", nMappingXConv, firstBinMappingXConv, lastBinMappingXConv, nMappingYConv, firstBinMappingYConv, lastBinMappingYConv);
		histoXY->SetXTitle("X (cm)");
		histoXY->SetYTitle("Y (cm)");
		histoZPhi = new   TH2F("ESD_ConversionMapping_ZPhi","", nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
		histoZPhi->SetXTitle("Z (cm)");
		histoZPhi->SetYTitle("Phi (rad)");
		histoRPhi = new   TH2F("ESD_ConversionMapping_RPhi","",nMappingRConv, firstBinMappingRConv, lastBinMappingRConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
		histoRPhi->SetXTitle("R (cm)");
		histoRPhi->SetYTitle("Phi (rad)");
		histoZRMidPt =  new TH2F("ESD_ConversionMappingMidPt_ZR","", nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingRConv, firstBinMappingRConv, lastBinMappingRConv);
		histoZRMidPt->SetXTitle("Z (cm)");
		histoZRMidPt->SetYTitle("R (cm)");
		histoZPhiMidPt = new    TH2F("ESD_ConversionMappingMidPt_ZPhi","", nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
		histoZPhiMidPt->SetXTitle("Z (cm)");
		histoZPhiMidPt->SetYTitle("Phi (rad)");
		histoRPhiMidPt = new    TH2F("ESD_ConversionMappingMidPt_RPhi","",nMappingRConv, firstBinMappingRConv, lastBinMappingRConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
		histoRPhiMidPt->SetXTitle("R (cm)");
		histoRPhiMidPt->SetYTitle("Phi (rad)");
		
		histoZPhi_SPD = new  TH2F( "ESD_ConversionMapping_SPD_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
		histoZPhi_SPD->SetXTitle("Z (cm)");
		histoZPhi_SPD->SetYTitle("Phi (rad)");
		histoZPhi_SPDTh = new   TH2F( "ESD_ConversionMapping_SPDTh_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);  
		histoZPhi_SPDTh->SetXTitle("Z (cm)");
		histoZPhi_SPDTh->SetYTitle("Phi (rad)");
		histoZPhi_SSD = new  TH2F("ESD_ConversionMapping_SSD_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
		histoZPhi_SSD->SetXTitle("Z (cm)");
		histoZPhi_SSD->SetYTitle("Phi (rad)");
		histoZPhi_SDD = new  TH2F("ESD_ConversionMapping_SDD_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);  
		histoZPhi_SDD->SetXTitle("Z (cm)");
		histoZPhi_SDD->SetYTitle("Phi (rad)");
		histoZPhi_ITSTPC = new  TH2F("ESD_ConversionMapping_ITSTPC_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
		histoZPhi_ITSTPC->SetXTitle("Z (cm)");
		histoZPhi_ITSTPC->SetYTitle("Phi (rad)");
		histoZPhi_HotZone = new    TH2F("ESD_ConversionMapping_HotZone_ZPhi","",nMappingZConvHotZone, firstBinMappingZConvHotZone, lastBinMappingZConvHotZone, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
		histoZPhi_HotZone->SetXTitle("Z (cm)");
		histoZPhi_HotZone->SetYTitle("Phi (rad)");

		histoXY_HotZone = new   TH2F("ESD_ConversionMapping_HotZone_XY" ,"" , nMappingXConvHotZone, firstBinMappingXConvHotZone, lastBinMappingXConvHotZone, nMappingYConvHotZone, firstBinMappingYConvHotZone, lastBinMappingYConvHotZone);
		histoXY_HotZone->SetXTitle("X (cm)");
		histoXY_HotZone->SetYTitle("Y (cm)");

		histoXY_Beampipe = new  TH2F("ESD_ConversionMappingInnerBeampipe_XY" ,"" , nMappingXBPConv, firstBinMappingXConvBP, lastBinMappingXConvBP, nMappingYBPConv, firstBinMappingYConvBP, lastBinMappingYConvBP); 
		histoXY_Beampipe->SetXTitle("X (cm)");
		histoXY_Beampipe->SetYTitle("Y (cm)");

		histoEta = new TH1F( "ESD_ConversionMapping_Eta","", nMappingEtaConv, firstBinMappingEtaConv, lastBinMappingEtaConv);   
		histoEta->SetXTitle("#eta");

		histoZPhi_SDD1 = new  TH2F("ESD_ConversionMapping_SDD1_ZPhi","",nMappingZSDD, firstBinMappingZSDD, lastBinMappingZSDD, nMappingPhiSDD, firstBinMappingPhiSDD, lastBinMappingPhiSDD);  
		histoZPhi_SDD1->SetXTitle("Z (cm)");
		histoZPhi_SDD1->SetYTitle("Phi (rad)");
		histoZPhi_SDD2 = new  TH2F("ESD_ConversionMapping_SDD2_ZPhi","",nMappingZSDD, firstBinMappingZSDD, lastBinMappingZSDD, nMappingPhiSDD, firstBinMappingPhiSDD, lastBinMappingPhiSDD);  
		histoZPhi_SDD2->SetXTitle("Z (cm)");
		histoZPhi_SDD2->SetYTitle("Phi (rad)");
		histoZPhi_SDDThermal = new  TH2F("ESD_ConversionMapping_SDDThermal_ZPhi","",nMappingZSDD, firstBinMappingZSDD, lastBinMappingZSDD, nMappingPhiSDD, firstBinMappingPhiSDD, lastBinMappingPhiSDD);  
		histoZPhi_SDDThermal->SetXTitle("Z (cm)");
		histoZPhi_SDDThermal->SetYTitle("Phi (rad)");

		histoRPhi_SDD1 = new  TH2F("ESD_ConversionMapping_SDD1_RPhi","",nMappingRSDD1, firstBinMappingRSDD1, lastBinMappingRSDD1, nMappingPhiSDD, firstBinMappingPhiSDD, lastBinMappingPhiSDD);  
		histoRPhi_SDD1->SetXTitle("R (cm)");
		histoRPhi_SDD1->SetYTitle("Phi (rad)");
		histoRPhi_SDD2 = new  TH2F("ESD_ConversionMapping_SDD2_RPhi","",nMappingRSDD2, firstBinMappingRSDD2, lastBinMappingRSDD2, nMappingPhiSDD, firstBinMappingPhiSDD, lastBinMappingPhiSDD);  
		histoRPhi_SDD2->SetXTitle("R (cm)");
		histoRPhi_SDD2->SetYTitle("Phi (rad)");
		histoRPhi_SDDThermal = new  TH2F("ESD_ConversionMapping_SDDThermal_RPhi","",nMappingRSDDthermal, firstBinMappingRSDDthermal, lastBinMappingRSDDthermal, nMappingPhiSDD, firstBinMappingPhiSDD, lastBinMappingPhiSDD);  
		histoRPhi_SDDThermal->SetXTitle("R (cm)");
		histoRPhi_SDDThermal->SetYTitle("Phi (rad)");

		histoZR_SDD1 = new  TH2F("ESD_ConversionMapping_SDD1_ZR","",nMappingZSDD, firstBinMappingZSDD, lastBinMappingZSDD, nMappingRSDD1, firstBinMappingRSDD1, lastBinMappingRSDD1);  
		histoZR_SDD1->SetXTitle("Z (cm)");
		histoZR_SDD1->SetYTitle("R (cm)");
		histoZR_SDD2 = new  TH2F("ESD_ConversionMapping_SDD2_ZR","",nMappingZSDD, firstBinMappingZSDD, lastBinMappingZSDD, nMappingRSDD2, firstBinMappingRSDD2, lastBinMappingRSDD2);  
		histoZR_SDD2->SetXTitle("Z (cm)");
		histoZR_SDD2->SetYTitle("R (cm)");
		histoZR_SDDThermal = new  TH2F("ESD_ConversionMapping_SDDThermal_ZR","",nMappingZSDD, firstBinMappingZSDD, lastBinMappingZSDD, nMappingRSDDthermal, firstBinMappingRSDDthermal, lastBinMappingRSDDthermal);  
		histoZR_SDDThermal->SetXTitle("Z (cm)");
		histoZR_SDDThermal->SetYTitle("R (cm)");
		
		histoChi2PerDOF = new TH1F( "ESD_ConvQual_Chi2PerDOF","", nChi2Conv, firstBinChi2Conv, lastBinChi2Conv);
		histoChi2PerDOF->SetXTitle("#chi^{2}/ndf");
		
		histoPt = new TH1F( "ESD_ConvQual_Pt","", nPtConv, firstBinPtConv, lastBinPtConv);
		histoPt->SetXTitle("p_{T, #gamma}");
		SetLogBinningTH1(histoPt);
	//       
		histoChi2PerDOF_AftCuts = new TH1F( "ESD_ConvQualAftCuts_Chi2PerDOF","", nChi2Conv, firstBinChi2Conv, lastBinChi2Conv);
		histoChi2PerDOF_AftCuts->SetXTitle("#chi^{2}/ndf");
		histoPt_AftCuts = new TH1F( "ESD_ConvQualAftCuts_Pt","", nPtConv, firstBinPtConv, lastBinPtConv);
		histoPt_AftCuts->SetXTitle("p_{T, #gamma}");
		SetLogBinningTH1(histoPt_AftCuts);
	} else { 
		histoEventQual =  (TH1I*)directoryConv->Get("ESD_EventQuality");
		histoZVtxDist =   (TH1F*)directoryConv->Get("Z_Vertex_distribution");
		histoGoodESDtrk09 =  (TH1F*)directoryConv->Get("ESD_NumberOfGoodESDTracks09");
		histoGoodESDtrk09Vtx =  (TH1F*)directoryConv->Get("ESD_NumberOfGoodESDTracks09Vtx");
		histoGoodESDtrk14 =  (TH1F*)directoryConv->Get("ESD_NumberOfGoodESDTracks14");
		histoGoodESDtrk14Vtx =  (TH1F*)directoryConv->Get("ESD_NumberOfGoodESDTracks14Vtx");
		histoGoodESDtrk0914 =   (TH1F*)directoryConv->Get("ESD_NumberOfGoodESDTracks0914");
		histoGoodESDtrk0914Vtx =   (TH1F*)directoryConv->Get("ESD_NumberOfGoodESDTracks0914Vtx");
		
		histoRPhi =    (TH2F*)directoryConv->Get("ESD_ConversionMapping_RPhi");
		histoXY =   (TH2F*)directoryConv->Get("ESD_ConversionMapping_XY");
		histoZR =   (TH2F*)directoryConv->Get("ESD_ConversionMapping_ZR");
		histoZPhi =    (TH2F*)directoryConv->Get("ESD_ConversionMapping_ZPhi");
		histoRPhiMidPt =  (TH2F*)directoryConv->Get("ESD_ConversionMappingMidPt_RPhi");
		histoZRMidPt =    (TH2F*)directoryConv->Get("ESD_ConversionMappingMidPt_ZR");
		histoZPhiMidPt =  (TH2F*)directoryConv->Get("ESD_ConversionMappingMidPt_ZPhi");
		histoZPhi_SPD =   (TH2F*)directoryConv->Get("ESD_ConversionMapping_SPD_ZPhi"); 
		histoZPhi_SPDTh =    (TH2F*)directoryConv->Get("ESD_ConversionMapping_SPDTh_ZPhi"); 
		histoZPhi_SSD =   (TH2F*)directoryConv->Get("ESD_ConversionMapping_SSD_ZPhi"); 
		histoZPhi_SDD =   (TH2F*)directoryConv->Get("ESD_ConversionMapping_SDD_ZPhi"); 
		histoZPhi_ITSTPC =   (TH2F*)directoryConv->Get("ESD_ConversionMapping_ITSTPC_ZPhi"); 
		histoZPhi_HotZone =  (TH2F*)directoryConv->Get("ESD_ConversionMapping_HotZone_ZPhi"); 
		histoXY_HotZone =    (TH2F*)directoryConv->Get("ESD_ConversionMapping_HotZone_XY"); 
		histoXY_Beampipe =   (TH2F*)directoryConv->Get("ESD_ConversionMappingInnerBeampipe_XY"); 
		histoEta =  (TH1F*)directoryConv->Get("ESD_ConversionMapping_Eta");
		histoZPhi_SDD1 = 		(TH2F*)directoryConv->Get("ESD_ConversionMapping_SDD1_ZPhi"); 
		histoZPhi_SDD2 = 		(TH2F*)directoryConv->Get("ESD_ConversionMapping_SDD2_ZPhi"); 
		histoZPhi_SDDThermal = 	(TH2F*)directoryConv->Get("ESD_ConversionMapping_SDDThermal_ZPhi"); 
		histoRPhi_SDD1 = 		(TH2F*)directoryConv->Get("ESD_ConversionMapping_SDD1_RPhi"); 
		histoRPhi_SDD2 = 		(TH2F*)directoryConv->Get("ESD_ConversionMapping_SDD2_RPhi"); 
		histoRPhi_SDDThermal = 	(TH2F*)directoryConv->Get("ESD_ConversionMapping_SDDThermal_RPhi"); 
		histoZR_SDD1 = 			(TH2F*)directoryConv->Get("ESD_ConversionMapping_SDD1_ZR"); 
		histoZR_SDD2 = 			(TH2F*)directoryConv->Get("ESD_ConversionMapping_SDD2_ZR"); 
		histoZR_SDDThermal = 	(TH2F*)directoryConv->Get("ESD_ConversionMapping_SDDThermal_ZR"); 
		
		histoChi2PerDOF =    (TH1F*)directoryConv->Get("ESD_ConvQual_Chi2PerDOF");
		histoPt =   (TH1F*)directoryConv->Get("ESD_ConvQual_Pt");
		
		histoChi2PerDOF_AftCuts =  (TH1F*)directoryConv->Get("ESD_ConvQualAftCuts_Chi2PerDOF");
		histoPt_AftCuts =    (TH1F*)directoryConv->Get("ESD_ConvQualAftCuts_Pt");
		
		histoEventQual->Sumw2();
		histoZVtxDist->Sumw2();
		histoGoodESDtrk09->Sumw2();
		histoGoodESDtrk09Vtx->Sumw2();
		histoGoodESDtrk14->Sumw2();
		histoGoodESDtrk14Vtx->Sumw2();
		histoGoodESDtrk0914->Sumw2();
		histoGoodESDtrk0914Vtx->Sumw2();
		histoRPhi->Sumw2();
		histoXY->Sumw2();
		histoZR->Sumw2();
		histoZPhi->Sumw2();
		histoRPhiMidPt->Sumw2();
		histoZRMidPt->Sumw2();
		histoZPhiMidPt->Sumw2();
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
		histoPt->Sumw2();
		histoChi2PerDOF_AftCuts->Sumw2();
		histoPt_AftCuts->Sumw2();
		
		histoZPhi_SDD1->Sumw2();
		histoZPhi_SDD2->Sumw2();
		histoZPhi_SDDThermal->Sumw2();
		histoRPhi_SDD1->Sumw2();
		histoRPhi_SDD2->Sumw2();
		histoRPhi_SDDThermal->Sumw2();
		histoZR_SDD1->Sumw2();
		histoZR_SDD2->Sumw2();
		histoZR_SDDThermal->Sumw2();

	}  
	//********************************************************************************
	//*      Reading of Tree with event information/ filling of histograms        *
	//********************************************************************************     
	
	Long64_t nbytesE = 0;
	for (Long64_t i=0; i<nEntriesEvent;i++) {
		nbytesE += Event->GetEntry(i);
		if (fMultiplicity09 > maxMult || fMultiplicity09 < minMult) continue;
		histoGoodESDtrk09->Fill(fMultiplicity09);
		histoGoodESDtrk14->Fill(fMultiplicity14);
		histoGoodESDtrk0914->Fill(fMultiplicity14-fMultiplicity09);
		if (nContrVtx>=1){
			histoZVtxDist->Fill(primVtxZ);
			histoGoodESDtrk09Vtx->Fill(fMultiplicity09);
			histoGoodESDtrk14Vtx->Fill(fMultiplicity14);
			histoGoodESDtrk0914Vtx->Fill(fMultiplicity14-fMultiplicity09);
		} 
	}
	if (histoEventQual->GetBinContent(1)>0){
		histoEventQual->SetBinContent(1,histoEventQual->GetBinContent(1)+nEntriesEvent);
	} else {
		histoEventQual->SetBinContent(1,nEntriesEvent);
	}
	
	//********************************************************************************
	//*      Reading of Tree with reconstructed gammas/ filling of histograms        *
	//********************************************************************************
	
	Long64_t nEntriesRecGam = ConvPointsRec->GetEntries();
	cout << "Number of Gammas to be processed: " << nEntriesRecGam << endl;
	
	//    Long64_t nbytesRecGamForFit = 0;
	//    Long64_t entriesInnerVessel =0;
	//    for (Long64_t i=0; i<nEntriesRecGam;i++) {
	//       nbytesRecGamForFit += ConvPointsRec->GetEntry(i);
	//       if (fR>59.5 && fR < 64){
	//          entriesInnerVessel++;
	//       }
	//    }
	//    cout << "Number of Gammas in inner vessel: " << entriesInnerVessel << endl;
	//    grInnerVessel = new TGraph(entriesInnerVessel);
	//    nbytesRecGamForFit = 0;
	//    Long64_t entriesDummy = 0;
	//    for (Long64_t i=0; i<nEntriesRecGam;i++) {
	//       nbytesRecGamForFit += ConvPointsRec->GetEntry(i);
	//       if (fR>59.5 && fR < 64){
	//          grInnerVessel->SetPoint(entriesDummy, fX,fY);
	//          entriesDummy++;
	//       }
	//    }
	//       
	//    
	//    TVirtualFitter::SetDefaultFitter("Minuit");  //default is Minuit
	//    TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 3);
	//    fitter->SetFCN(circleFittingInnerVesselRec);
	// 
	//    fitter->SetParameter(0, "x0",   0, 0.1, 0,0);
	//    fitter->SetParameter(1, "y0",   0, 0.1, 0,0);
	//    fitter->SetParameter(2, "R",    1, 0.1, 0,0);
	// 
	//    Double_t arglist[1] = {0};
	//    fitter->ExecuteCommand("MIGRAD", arglist, 0);
	
	Long64_t nbytesRecGam = 0;
	for (Long64_t i=0; i<nEntriesRecGam;i++) {
		nbytesRecGam += ConvPointsRec->GetEntry(i);
		if (fMultiplicity09Photon > maxMult || fMultiplicity09Photon < minMult) continue;
		histoChi2PerDOF->Fill(fChiPerDOF);  
		histoPt->Fill(fPt);
		Float_t fR = (*fCoord)[3];
		Float_t fPhi = (*fCoord)[4];
		Float_t fX = (*fCoord)[0];
		Float_t fY = (*fCoord)[1];
		Float_t fZ = (*fCoord)[2];
		Float_t fEta = TMath::Log(TMath::Tan(fTheta/2.));
		Float_t fPosiEta = TMath::Log(TMath::Tan((*fDaugtherProp)[1]/2.));
		Float_t fElecEta = TMath::Log(TMath::Tan((*fDaugtherProp)[3]/2.));      
		if ( fPt > minPtPhotonCut ){ //fChiPerDOF < chi2PerDOFCut &&
			if ( TMath::Abs(fPosiEta) > etaMaxCut) continue;
			if ( TMath::Abs(fElecEta) > etaMaxCut) continue;
			if ( TMath::Abs(fEta) > etaMaxCut) continue;

	//          cout << fR << "\t" << fZ << endl;
			if (!(fR <= (TMath::Abs(fZ) * fLineCutZRSlope - fLineCutZVtx))){ 
	//             if (TMath::Abs(fEta) >=etaMinCut && TMath::Abs(fEta) <= etaMaxCut){
				if ( !(etaMinCut != -0.1 &&   fR >= (TMath::Abs(fZ) * fLineCutZRSlopeMin- fLineCutZVtxMin))){
					if (fR > rMin ){ //
	//                      cout << "Eta positron: " << fPosiEta << "\t Eta electron: " << fElecEta << "\t Photon Eta: " << fEta<< endl;
						if (fPt > midPtPhotonMin && fPt < midPtPhotonMax){
							histoZRMidPt->Fill(fZ,fR);
							histoZPhiMidPt->Fill(fZ,fPhi);
							histoRPhiMidPt->Fill(fR,fPhi);
						}
						histoZR->Fill(fZ,fR);
						histoXY->Fill(fX,fY);
						histoZPhi->Fill(fZ,fPhi);
						histoRPhi->Fill(fR,fPhi);
						histoEta->Fill(fEta);
						if ( fR < rSPD){
							histoZPhi_SPD->Fill(fZ,fPhi);
						} else if (fR < rSPDth){
							histoZPhi_SPDTh->Fill(fZ,fPhi);
						} else if (fR < rSDD){
							histoZPhi_SDD->Fill(fZ,fPhi);
						} else if (fR < rSSD){
							histoZPhi_SSD->Fill(fZ,fPhi);
						} else if (fR < rTPC){
							histoZPhi_ITSTPC->Fill(fZ,fPhi);
						}
						
						if (fR > rDetailsSDD1Min && fR < rDetailsSDD1Max ){
							histoRPhi_SDD1->Fill(fR,fPhi);
							histoZPhi_SDD1->Fill(fZ,fPhi);
							histoZR_SDD1->Fill(fZ,fR);
						}
						if (fR > rDetailsSDD2Min && fR < rDetailsSDD2Max ){
							histoRPhi_SDD2->Fill(fR,fPhi);
							histoZPhi_SDD2->Fill(fZ,fPhi);
							histoZR_SDD2->Fill(fZ,fR);
						}
						if (fR > rDetailsSDDThermalShieldMin && fR < rDetailsSDDThermalShieldMax ){
							histoRPhi_SDDThermal->Fill(fR,fPhi);
							histoZPhi_SDDThermal->Fill(fZ,fPhi);
							histoZR_SDDThermal->Fill(fZ,fR);
						}
						
						
						if (fR>rHotZoneMin && fR < rHotZoneMax){
							histoZPhi_HotZone->Fill( fZ,fPhi);
							if (fZ < zHotZoneMax && fZ > zHotZoneMin){
								histoXY_HotZone->Fill( fX,fY);
							}
						}
						if (fZ < zBeamPipeInner && fZ > - zBeamPipeInner) {
							histoXY_Beampipe->Fill( fX,fY);
						}  
						histoChi2PerDOF_AftCuts->Fill(fChiPerDOF);         
						histoPt_AftCuts->Fill(fPt);
					}
				}
			}
	//          }
		}
	}
	
	//********************************************************************************
	//*                  Writing histograms to outputfile                            *
	//********************************************************************************     
	TFile* fileMappingWrite = new TFile(fileNameOutput.Data(),"UPDATE");
	fileMappingWrite->mkdir(Form("GammaConv_%s",  CutSelection.Data()));
	fileMappingWrite->cd(Form("GammaConv_%s",  CutSelection.Data()));
		histoEventQual->Write("ESD_EventQuality",TObject::kWriteDelete);     
		histoZVtxDist->Write("Z_Vertex_distribution",TObject::kWriteDelete);    
		histoGoodESDtrk09->Write("ESD_NumberOfGoodESDTracks09",TObject::kWriteDelete);      
		histoGoodESDtrk09Vtx->Write("ESD_NumberOfGoodESDTracks09Vtx",TObject::kWriteDelete);      
		histoGoodESDtrk14->Write("ESD_NumberOfGoodESDTracks14",TObject::kWriteDelete);      
		histoGoodESDtrk14Vtx->Write("ESD_NumberOfGoodESDTracks14Vtx",TObject::kWriteDelete);      
		histoGoodESDtrk0914->Write("ESD_NumberOfGoodESDTracks0914",TObject::kWriteDelete);     
		histoGoodESDtrk0914Vtx->Write("ESD_NumberOfGoodESDTracks0914Vtx",TObject::kWriteDelete);     
		histoRPhi->Write("ESD_ConversionMapping_RPhi",TObject::kWriteDelete);
		histoXY->Write("ESD_ConversionMapping_XY",TObject::kWriteDelete);
		histoZR->Write("ESD_ConversionMapping_ZR",TObject::kWriteDelete);
		histoZPhi->Write("ESD_ConversionMapping_ZPhi",TObject::kWriteDelete);
		histoRPhiMidPt->Write("ESD_ConversionMappingMidPt_RPhi",TObject::kWriteDelete);
		histoZRMidPt->Write("ESD_ConversionMappingMidPt_ZR",TObject::kWriteDelete);
		histoZPhiMidPt->Write("ESD_ConversionMappingMidPt_ZPhi",TObject::kWriteDelete);
		histoZPhi_SPD->Write("ESD_ConversionMapping_SPD_ZPhi",TObject::kWriteDelete); 
		histoZPhi_SPDTh->Write("ESD_ConversionMapping_SPDTh_ZPhi",TObject::kWriteDelete); 
		histoZPhi_SSD->Write("ESD_ConversionMapping_SSD_ZPhi",TObject::kWriteDelete); 
		histoZPhi_SDD->Write("ESD_ConversionMapping_SDD_ZPhi",TObject::kWriteDelete); 
		histoZPhi_ITSTPC->Write("ESD_ConversionMapping_ITSTPC_ZPhi",TObject::kWriteDelete); 
		histoZPhi_HotZone->Write("ESD_ConversionMapping_HotZone_ZPhi",TObject::kWriteDelete); 
		histoXY_HotZone->Write("ESD_ConversionMapping_HotZone_XY",TObject::kWriteDelete); 
		histoXY_Beampipe->Write("ESD_ConversionMappingInnerBeampipe_XY",TObject::kWriteDelete); 
		histoEta->Write("ESD_ConversionMapping_Eta",TObject::kWriteDelete);
		histoZPhi_SDD1->Write("ESD_ConversionMapping_SDD1_ZPhi",TObject::kWriteDelete);
		histoZPhi_SDD2->Write("ESD_ConversionMapping_SDD2_ZPhi",TObject::kWriteDelete); 
		histoZPhi_SDDThermal->Write("ESD_ConversionMapping_SDDThermal_ZPhi",TObject::kWriteDelete);
		histoRPhi_SDD1->Write("ESD_ConversionMapping_SDD1_RPhi",TObject::kWriteDelete);
		histoRPhi_SDD2->Write("ESD_ConversionMapping_SDD2_RPhi",TObject::kWriteDelete);
		histoRPhi_SDDThermal->Write("ESD_ConversionMapping_SDDThermal_RPhi",TObject::kWriteDelete);
		histoZR_SDD1->Write("ESD_ConversionMapping_SDD1_ZR",TObject::kWriteDelete);
		histoZR_SDD2->Write("ESD_ConversionMapping_SDD2_ZR",TObject::kWriteDelete);
		histoZR_SDDThermal->Write("ESD_ConversionMapping_SDDThermal_ZR",TObject::kWriteDelete);
		histoChi2PerDOF->Write("ESD_ConvQual_Chi2PerDOF",TObject::kWriteDelete);
		histoPt->Write("ESD_ConvQual_Pt",TObject::kWriteDelete);
		histoChi2PerDOF_AftCuts->Write("ESD_ConvQualAftCuts_Chi2PerDOF",TObject::kWriteDelete);
		histoPt_AftCuts->Write("ESD_ConvQualAftCuts_Pt",TObject::kWriteDelete);
	fileMappingWrite->Write();
	fileMappingWrite->Close();
	
	delete histoEventQual;
	delete histoZVtxDist;
	delete histoGoodESDtrk09;
	delete histoGoodESDtrk09Vtx;
	delete histoGoodESDtrk14;
	delete histoGoodESDtrk14Vtx;
	delete histoGoodESDtrk0914;
	delete histoGoodESDtrk0914Vtx;
	delete histoRPhi;
	delete histoXY;
	delete histoZR;
	delete histoZPhi;
	delete histoRPhiMidPt;
	delete histoZRMidPt;
	delete histoZPhiMidPt;
	delete histoZPhi_SPD; 
	delete histoZPhi_SPDTh; 
	delete histoZPhi_SSD; 
	delete histoZPhi_SDD; 
	delete histoZPhi_ITSTPC; 
	delete histoZPhi_HotZone; 
	delete histoXY_HotZone; 
	delete histoXY_Beampipe; 
	delete histoEta;
	delete histoZPhi_SDD1;
	delete histoZPhi_SDD2;
	delete histoZPhi_SDDThermal;
	delete histoRPhi_SDD1;
	delete histoRPhi_SDD2;
	delete histoRPhi_SDDThermal;
	delete histoZR_SDD1;
	delete histoZR_SDD2;
	delete histoZR_SDDThermal;
	delete histoChi2PerDOF;
	delete histoPt;
	delete histoChi2PerDOF_AftCuts;
	delete histoPt_AftCuts;
	
	if (!kMC) return;

	
	//********************************************************************************
	//*      Definition of histograms for MC validated Conversion Points             *
	//********************************************************************************     
	//_________________true conversion__________________________________________
	TH2F* histoTrueZR =0x00;
	TH2F* histoTrueXY =0x00;
	TH2F* histoTrueZPhi =0x00;
	TH2F* histoTrueRPhi =0x00;
	TH2F* histoTrueZPhi_SPD =0x00;
	TH2F* histoTrueZPhi_SPDTh =0x00;
	TH2F* histoTrueZPhi_SSD =0x00;
	TH2F* histoTrueZPhi_SDD =0x00;   
	TH2F* histoTrueZPhi_ITSTPC =0x00;
	TH2F* histoTrueZPhi_HotZone =0x00;
	TH2F* histoTrueXY_HotZone =0x00;
	TH2F* histoTrueXY_Beampipe =0x00;   
	TH1F* histoTrueEta =0x00;  
	TH1F* histoTruePt =0x00;
	TH1F* histoTruePt_AftCuts =0x00;
	//________________________ true combinatorics _____________________________________
	TH2F* histoTrueCombZR = 0x00;
	TH2F* histoTrueCombElecZR = 0x00;
	TH2F* histoTrueCombPionZR = 0x00;
	TH2F* histoTrueCombElecPionZR = 0x00;
	TH2F* histoTrueCombPionProtonZR = 0x00;
	TH2F* histoTrueCombKaonXZR = 0x00;
	TH2F* histoTrueCombZPhi =0x00;
	TH2F* histoTrueCombRPhi =0x00;
	TH2F* histoTrueCombZPhi_SPD =0x00;
	TH2F* histoTrueCombZPhi_SPDTh =0x00;
	TH2F* histoTrueCombZPhi_SSD =0x00;
	TH2F* histoTrueCombZPhi_SDD =0x00;  
	TH2F* histoTrueCombZPhi_ITSTPC =0x00;
	TH2F* histoTrueCombZPhi_HotZone =0x00;
	TH1F* histoTrueCombEta =0x00; 
	TH1F* histoTrueCombElecEta =0x00;   
	TH1F* histoTrueCombPionEta =0x00;   
	TH1F* histoTrueCombElecPionEta =0x00;  
	TH1F* histoTrueCombPionProtonEta =0x00;   
	TH1F* histoTrueCombKaonXEta =0x00;  
	//________________________ true primary conversions _________________________________
	TH2F* histoTruePrimZR = 0x00;
	TH2F* histoTruePrimZPhi =0x00;
	TH2F* histoTruePrimRPhi =0x00;
	TH2F* histoTruePrimZPhi_SPD =0x00;
	TH2F* histoTruePrimZPhi_SPDTh =0x00;
	TH2F* histoTruePrimZPhi_SSD =0x00;
	TH2F* histoTruePrimZPhi_SDD =0x00;  
	TH2F* histoTruePrimZPhi_ITSTPC =0x00;
	TH2F* histoTruePrimZPhi_HotZone =0x00;
	TH1F* histoTruePrimEta =0x00; 
	//________________________ true hadronic interactions _________________________________
	TH2F* histoTrueHadZR = 0x00;
	TH1F* histoTrueHadEta =0x00;  
	//________________________ true Dalitz Pi0 _________________________________
	TH2F* histoTrueDalPi0ZR = 0x00;
	TH2F* histoTrueDalPi0ZPhi_SPD =0x00;
	TH2F* histoTrueDalPi0ZPhi_SPDTh =0x00;
	TH2F* histoTrueDalPi0ZPhi_SDD =0x00;
	TH1F* histoTrueDalPi0Eta =0x00;  
	//________________________ true Dalitz Eta _________________________________
	TH2F* histoTrueDalEtaZR = 0x00;
	TH2F* histoTrueDalEtaZPhi_SPD =0x00;
	TH2F* histoTrueDalEtaZPhi_SPDTh =0x00;
	TH2F* histoTrueDalEtaZPhi_SDD =0x00;
	TH1F* histoTrueDalEtaEta =0x00;
	
	if (kMC){
		if (!merge || directoryConv==0) {      
			histoTrueZR =  new TH2F("ESD_TrueConversionMapping_ZR","", nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingRConv, firstBinMappingRConv, lastBinMappingRConv);
			histoTrueXY =  new TH2F("ESD_TrueConversionMapping_XY","", nMappingXConv, firstBinMappingXConv, lastBinMappingXConv, nMappingYConv, firstBinMappingYConv, lastBinMappingYConv);
			histoTrueZPhi = new  TH2F("ESD_TrueConversionMapping_ZPhi","", nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTrueRPhi = new  TH2F("ESD_TrueConversionMapping_RPhi","",nMappingRConv, firstBinMappingRConv, lastBinMappingRConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTrueZPhi_SPD = new    TH2F( "ESD_TrueConversionMapping_SPD_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTrueZPhi_SPDTh = new  TH2F( "ESD_TrueConversionMapping_SPDTh_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTrueZPhi_SSD = new    TH2F("ESD_TrueConversionMapping_SSD_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTrueZPhi_SDD = new    TH2F("ESD_TrueConversionMapping_SDD_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv); 
			histoTrueZPhi_ITSTPC = new    TH2F("ESD_TrueConversionMapping_ITSTPC_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTrueZPhi_HotZone = new   TH2F("ESD_TrueConversionMapping_HotZone_ZPhi","",nMappingZConvHotZone, firstBinMappingZConvHotZone, lastBinMappingZConvHotZone, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTrueXY_HotZone = new  TH2F("ESD_TrueConversionMapping_HotZone_XY" ,"" , nMappingXConvHotZone, firstBinMappingXConvHotZone, lastBinMappingXConvHotZone, nMappingYConvHotZone, firstBinMappingYConvHotZone, lastBinMappingYConvHotZone);
			histoTrueXY_Beampipe = new    TH2F("ESD_TrueConversionMappingInnerBeampipe_XY" ,"" , nMappingXBPConv, firstBinMappingXConvBP, lastBinMappingXConvBP, nMappingYBPConv, firstBinMappingYConvBP, lastBinMappingYConvBP);   
			histoTrueEta = new TH1F( "ESD_TrueConversionMapping_Eta","", nMappingEtaConv, firstBinMappingEtaConv, lastBinMappingEtaConv); 
			histoTruePt = new TH1F( "ESD_TrueConvQual_Pt","", nPtConv, firstBinPtConv, lastBinPtConv);
			SetLogBinningTH1(histoTruePt);
			histoTruePt_AftCuts = new TH1F( "ESD_TrueConvQualAftCuts_Pt","", nPtConv, firstBinPtConv, lastBinPtConv);
			SetLogBinningTH1(histoTruePt_AftCuts);
			//________________________ true combinatorics _____________________________________
			histoTrueCombZR =  new TH2F("ESD_TrueCombMapping_ZR","", nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingRConv, firstBinMappingRConv, lastBinMappingRConv);
			histoTrueCombElecZR =  new TH2F("ESD_TrueCombElecMapping_ZR","", nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingRConv, firstBinMappingRConv, lastBinMappingRConv);
			histoTrueCombPionZR =  new TH2F("ESD_TrueCombPionMapping_ZR","", nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingRConv, firstBinMappingRConv, lastBinMappingRConv);
			histoTrueCombElecPionZR =  new TH2F("ESD_TrueCombElecPionMapping_ZR","", nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingRConv, firstBinMappingRConv, lastBinMappingRConv);
			histoTrueCombPionProtonZR =  new TH2F("ESD_TrueCombPionProtonMapping_ZR","", nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingRConv, firstBinMappingRConv, lastBinMappingRConv);
			histoTrueCombKaonXZR =  new TH2F("ESD_TrueCombKaonXMapping_ZR","", nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingRConv, firstBinMappingRConv, lastBinMappingRConv);
			histoTrueCombZPhi = new    TH2F("ESD_TrueCombMapping_ZPhi","", nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTrueCombRPhi = new    TH2F("ESD_TrueCombMapping_RPhi","",nMappingRConv, firstBinMappingRConv, lastBinMappingRConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTrueCombZPhi_SPD = new   TH2F( "ESD_TrueCombMapping_SPD_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTrueCombZPhi_SPDTh = new    TH2F( "ESD_TrueCombMapping_SPDTh_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTrueCombZPhi_SSD = new   TH2F("ESD_TrueCombMapping_SSD_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTrueCombZPhi_SDD = new   TH2F("ESD_TrueCombMapping_SDD_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv); 
			histoTrueCombZPhi_ITSTPC = new   TH2F("ESD_TrueCombMapping_ITSTPC_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTrueCombZPhi_HotZone = new  TH2F("ESD_TrueCombMapping_HotZone_ZPhi","",nMappingZConvHotZone, firstBinMappingZConvHotZone, lastBinMappingZConvHotZone, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTrueCombEta = new TH1F( "ESD_TrueCombMapping_Eta","", nMappingEtaConv, firstBinMappingEtaConv, lastBinMappingEtaConv);   
			histoTrueCombElecEta = new TH1F( "ESD_TrueCombElecMapping_Eta","", nMappingEtaConv, firstBinMappingEtaConv, lastBinMappingEtaConv); 
			histoTrueCombPionEta = new TH1F( "ESD_TrueCombPionMapping_Eta","", nMappingEtaConv, firstBinMappingEtaConv, lastBinMappingEtaConv); 
			histoTrueCombElecPionEta = new TH1F( "ESD_TrueCombElecPionMapping_Eta","", nMappingEtaConv, firstBinMappingEtaConv, lastBinMappingEtaConv);  
			histoTrueCombPionProtonEta = new TH1F( "ESD_TrueCombPionProtonMapping_Eta","", nMappingEtaConv, firstBinMappingEtaConv, lastBinMappingEtaConv); 
			histoTrueCombKaonXEta = new TH1F( "ESD_TrueCombKaonXMapping_Eta","", nMappingEtaConv, firstBinMappingEtaConv, lastBinMappingEtaConv);  
			//________________________ true primary conversions _________________________________
			histoTruePrimZR =  new TH2F("ESD_TruePrimMapping_ZR","", nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingRConv, firstBinMappingRConv, lastBinMappingRConv);
			histoTruePrimZPhi = new    TH2F("ESD_TruePrimMapping_ZPhi","", nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTruePrimRPhi = new    TH2F("ESD_TruePrimMapping_RPhi","",nMappingRConv, firstBinMappingRConv, lastBinMappingRConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTruePrimZPhi_SPD = new   TH2F( "ESD_TruePrimMapping_SPD_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTruePrimZPhi_SPDTh = new    TH2F( "ESD_TruePrimMapping_SPDTh_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTruePrimZPhi_SSD = new   TH2F("ESD_TruePrimMapping_SSD_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTruePrimZPhi_SDD = new   TH2F("ESD_TruePrimMapping_SDD_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv); 
			histoTruePrimZPhi_ITSTPC = new   TH2F("ESD_TruePrimMapping_ITSTPC_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTruePrimZPhi_HotZone = new  TH2F("ESD_TruePrimMapping_HotZone_ZPhi","",nMappingZConvHotZone, firstBinMappingZConvHotZone, lastBinMappingZConvHotZone, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTruePrimEta = new TH1F( "ESD_TruePrimMapping_Eta","", nMappingEtaConv, firstBinMappingEtaConv, lastBinMappingEtaConv);   
			//________________________ true hadronic interactions _________________________________
			histoTrueHadZR =  new TH2F("ESD_TrueHadMapping_ZR","", nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingRConv, firstBinMappingRConv, lastBinMappingRConv);
			histoTrueHadEta = new TH1F( "ESD_TrueHadMapping_Eta","", nMappingEtaConv, firstBinMappingEtaConv, lastBinMappingEtaConv);  
			//________________________ true Dalitz Pi0 _________________________________
			histoTrueDalPi0ZR =  new TH2F("ESD_TrueDalPi0Mapping_ZR","", nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingRConv, firstBinMappingRConv, lastBinMappingRConv);
			histoTrueDalPi0ZPhi_SPD = new    TH2F( "ESD_TrueDalPi0Mapping_SPD_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTrueDalPi0ZPhi_SPDTh = new  TH2F( "ESD_TrueDalPi0Mapping_SPDTh_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTrueDalPi0ZPhi_SDD = new    TH2F("ESD_TrueDalPi0Mapping_SDD_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTrueDalPi0Eta = new TH1F( "ESD_TrueDalPi0Mapping_Eta","", nMappingEtaConv, firstBinMappingEtaConv, lastBinMappingEtaConv);  
			//________________________ true Dalitz Eta _________________________________
			histoTrueDalEtaZR =  new TH2F("ESD_TrueDalEtaMapping_ZR","", nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingRConv, firstBinMappingRConv, lastBinMappingRConv);
			histoTrueDalEtaZPhi_SPD = new    TH2F( "ESD_TrueDalEtaMapping_SPD_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTrueDalEtaZPhi_SPDTh = new  TH2F( "ESD_TrueDalEtaMapping_SPDTh_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTrueDalEtaZPhi_SDD = new    TH2F("ESD_TrueDalEtaMapping_SDD_ZPhi","",nMappingZConv, firstBinMappingZConv, lastBinMappingZConv, nMappingPhiConv, firstBinMappingPhiConv, lastBinMappingPhiConv);
			histoTrueDalEtaEta = new TH1F( "ESD_TrueDalEtaMapping_Eta","", nMappingEtaConv, firstBinMappingEtaConv, lastBinMappingEtaConv);  
		} else {
			histoTrueRPhi =   (TH2F*)directoryConv->Get("ESD_TrueConversionMapping_RPhi");
			histoTrueXY =  (TH2F*)directoryConv->Get("ESD_TrueConversionMapping_XY");
			histoTrueZR =  (TH2F*)directoryConv->Get("ESD_TrueConversionMapping_ZR");
			histoTrueZPhi =   (TH2F*)directoryConv->Get("ESD_TrueConversionMapping_ZPhi");
			histoTrueZPhi_SPD =  (TH2F*)directoryConv->Get("ESD_TrueConversionMapping_SPD_ZPhi"); 
			histoTrueZPhi_SPDTh =   (TH2F*)directoryConv->Get("ESD_TrueConversionMapping_SPDTh_ZPhi"); 
			histoTrueZPhi_SSD =  (TH2F*)directoryConv->Get("ESD_TrueConversionMapping_SSD_ZPhi"); 
			histoTrueZPhi_SDD =  (TH2F*)directoryConv->Get("ESD_TrueConversionMapping_SDD_ZPhi"); 
			histoTrueZPhi_ITSTPC =  (TH2F*)directoryConv->Get("ESD_TrueConversionMapping_ITSTPC_ZPhi"); 
			histoTrueZPhi_HotZone =    (TH2F*)directoryConv->Get("ESD_TrueConversionMapping_HotZone_ZPhi"); 
			histoTrueXY_HotZone =   (TH2F*)directoryConv->Get("ESD_TrueConversionMapping_HotZone_XY"); 
			histoTrueXY_Beampipe =  (TH2F*)directoryConv->Get("ESD_TrueConversionMappingInnerBeampipe_XY"); 
			histoTrueEta =    (TH1F*)directoryConv->Get("ESD_TrueConversionMapping_Eta");
			histoTruePt =  (TH1F*)directoryConv->Get("ESD_TrueConvQual_Pt");
			histoTruePt_AftCuts =   (TH1F*)directoryConv->Get("ESD_TrueConvQualAftCuts_Pt");
		
			histoTruePrimRPhi =  (TH2F*)directoryConv->Get("ESD_TruePrimMapping_RPhi");
			histoTruePrimZR =    (TH2F*)directoryConv->Get("ESD_TruePrimMapping_ZR");
			histoTruePrimZPhi =  (TH2F*)directoryConv->Get("ESD_TruePrimMapping_ZPhi");
			histoTruePrimZPhi_SPD =    (TH2F*)directoryConv->Get("ESD_TruePrimMapping_SPD_ZPhi"); 
			histoTruePrimZPhi_SPDTh =  (TH2F*)directoryConv->Get("ESD_TruePrimMapping_SPDTh_ZPhi"); 
			histoTruePrimZPhi_SSD =    (TH2F*)directoryConv->Get("ESD_TruePrimMapping_SSD_ZPhi"); 
			histoTruePrimZPhi_SDD =    (TH2F*)directoryConv->Get("ESD_TruePrimMapping_SDD_ZPhi"); 
			histoTruePrimZPhi_ITSTPC =    (TH2F*)directoryConv->Get("ESD_TruePrimMapping_ITSTPC_ZPhi"); 
			histoTruePrimZPhi_HotZone =   (TH2F*)directoryConv->Get("ESD_TruePrimMapping_HotZone_ZPhi"); 
			histoTruePrimEta =   (TH1F*)directoryConv->Get("ESD_TruePrimMapping_Eta");

			histoTrueCombRPhi =  (TH2F*)directoryConv->Get("ESD_TrueCombMapping_RPhi");
			histoTrueCombZR =    (TH2F*)directoryConv->Get("ESD_TrueCombMapping_ZR");
			histoTrueCombElecZR =   (TH2F*)directoryConv->Get("ESD_TrueCombElecMapping_ZR");
			histoTrueCombPionZR =   (TH2F*)directoryConv->Get("ESD_TrueCombPionMapping_ZR");
			histoTrueCombElecPionZR =  (TH2F*)directoryConv->Get("ESD_TrueCombElecPionMapping_ZR");
			histoTrueCombPionProtonZR =   (TH2F*)directoryConv->Get("ESD_TrueCombPionProtonMapping_ZR");
			histoTrueCombKaonXZR=   (TH2F*)directoryConv->Get("ESD_TrueCombKaonXMapping_ZR");
			histoTrueCombZPhi =  (TH2F*)directoryConv->Get("ESD_TrueCombMapping_ZPhi");
			histoTrueCombZPhi_SPD =    (TH2F*)directoryConv->Get("ESD_TrueCombMapping_SPD_ZPhi"); 
			histoTrueCombZPhi_SPDTh =  (TH2F*)directoryConv->Get("ESD_TrueCombMapping_SPDTh_ZPhi"); 
			histoTrueCombZPhi_SSD =    (TH2F*)directoryConv->Get("ESD_TrueCombMapping_SSD_ZPhi"); 
			histoTrueCombZPhi_SDD =    (TH2F*)directoryConv->Get("ESD_TrueCombMapping_SDD_ZPhi"); 
			histoTrueCombZPhi_ITSTPC =    (TH2F*)directoryConv->Get("ESD_TrueCombMapping_ITSTPC_ZPhi"); 
			histoTrueCombZPhi_HotZone =   (TH2F*)directoryConv->Get("ESD_TrueCombMapping_HotZone_ZPhi"); 
			histoTrueCombEta =   (TH1F*)directoryConv->Get("ESD_TrueCombMapping_Eta");
			histoTrueCombElecEta =  (TH1F*)directoryConv->Get("ESD_TrueCombElecMapping_Eta");
			histoTrueCombPionEta =  (TH1F*)directoryConv->Get("ESD_TrueCombPionMapping_Eta");
			histoTrueCombElecPionEta =    (TH1F*)directoryConv->Get("ESD_TrueCombElecPionMapping_Eta");
			histoTrueCombPionProtonEta =  (TH1F*)directoryConv->Get("ESD_TrueCombPionProtonMapping_Eta");
			histoTrueCombKaonXEta =    (TH1F*)directoryConv->Get("ESD_TrueCombKaonXMapping_Eta");
			
			histoTrueHadZR =  (TH2F*)directoryConv->Get("ESD_TrueHadMapping_ZR");
			histoTrueHadEta =    (TH1F*)directoryConv->Get("ESD_TrueHadMapping_Eta");

			histoTrueDalPi0ZR =  (TH2F*)directoryConv->Get("ESD_TrueDalPi0Mapping_ZR");
			histoTrueDalPi0ZPhi_SPD =  (TH2F*)directoryConv->Get("ESD_TrueDalPi0Mapping_SPD_ZPhi"); 
			histoTrueDalPi0ZPhi_SPDTh =   (TH2F*)directoryConv->Get("ESD_TrueDalPi0Mapping_SPDTh_ZPhi"); 
			histoTrueDalPi0ZPhi_SDD =  (TH2F*)directoryConv->Get("ESD_TrueDalPi0Mapping_SDD_ZPhi"); 
			histoTrueDalPi0Eta =    (TH1F*)directoryConv->Get("ESD_TrueDalPi0Mapping_Eta");

			histoTrueDalEtaZR =  (TH2F*)directoryConv->Get("ESD_TrueDalEtaMapping_ZR");
			histoTrueDalEtaZPhi_SPD =  (TH2F*)directoryConv->Get("ESD_TrueDalEtaMapping_SPD_ZPhi"); 
			histoTrueDalEtaZPhi_SPDTh =   (TH2F*)directoryConv->Get("ESD_TrueDalEtaMapping_SPDTh_ZPhi"); 
			histoTrueDalEtaZPhi_SDD =  (TH2F*)directoryConv->Get("ESD_TrueDalEtaMapping_SDD_ZPhi"); 
			histoTrueDalEtaEta =    (TH1F*)directoryConv->Get("ESD_TrueDalEtaMapping_Eta");

			
	//       histoTrueEventQual->Sumw2();
			histoTrueRPhi->Sumw2();
			histoTrueXY->Sumw2();
			histoTrueZR->Sumw2();
			histoTrueZPhi->Sumw2();
			histoTrueZPhi_SPD->Sumw2();
			histoTrueZPhi_SPDTh->Sumw2();
			histoTrueZPhi_SSD->Sumw2();
			histoTrueZPhi_SDD->Sumw2();
			histoTrueZPhi_ITSTPC->Sumw2();
			histoTrueZPhi_HotZone->Sumw2();
			histoTrueXY_HotZone->Sumw2();
			histoTrueXY_Beampipe->Sumw2();
			histoTrueEta->Sumw2();
			histoTruePt->Sumw2();
			histoTruePt_AftCuts->Sumw2();
			histoTruePrimRPhi->Sumw2();
			histoTruePrimZR->Sumw2();
			histoTruePrimZPhi->Sumw2();
			histoTruePrimZPhi_SPD->Sumw2();
			histoTruePrimZPhi_SPDTh->Sumw2();
			histoTruePrimZPhi_SSD->Sumw2();
			histoTruePrimZPhi_SDD->Sumw2();
			histoTruePrimZPhi_ITSTPC->Sumw2();
			histoTruePrimZPhi_HotZone->Sumw2();
			histoTruePrimEta->Sumw2();
			histoTrueCombRPhi->Sumw2();
			histoTrueCombZR->Sumw2();
			histoTrueCombElecZR->Sumw2();
			histoTrueCombPionZR->Sumw2();
			histoTrueCombElecPionZR->Sumw2();
			histoTrueCombPionProtonZR->Sumw2();
			histoTrueCombKaonXZR->Sumw2();
			histoTrueCombZPhi->Sumw2();
			histoTrueCombZPhi_SPD->Sumw2();
			histoTrueCombZPhi_SPDTh->Sumw2();
			histoTrueCombZPhi_SSD->Sumw2();
			histoTrueCombZPhi_SDD->Sumw2();
			histoTrueCombZPhi_ITSTPC->Sumw2();
			histoTrueCombZPhi_HotZone->Sumw2();
			histoTrueCombEta->Sumw2();
			histoTrueCombElecEta->Sumw2();
			histoTrueCombPionEta->Sumw2();
			histoTrueCombElecPionEta->Sumw2();
			histoTrueCombPionProtonEta->Sumw2();
			histoTrueCombKaonXEta->Sumw2();
			histoTrueHadZR->Sumw2();
			histoTrueHadEta->Sumw2();
			histoTrueDalPi0ZR->Sumw2();
			histoTrueDalPi0ZPhi_SPD->Sumw2();
			histoTrueDalPi0ZPhi_SPDTh->Sumw2();
			histoTrueDalPi0ZPhi_SDD->Sumw2();
			histoTrueDalPi0Eta->Sumw2();
			histoTrueDalEtaZR->Sumw2();
			histoTrueDalEtaZPhi_SPD->Sumw2();
			histoTrueDalEtaZPhi_SPDTh->Sumw2();
			histoTrueDalEtaZPhi_SDD->Sumw2();
			histoTrueDalEtaEta->Sumw2();

		}
	//********************************************************************************
	//*      Reading of Tree with MC validated gammas/ filling of histograms         *
	//********************************************************************************     
		Long64_t nEntriesTrueGamma = ConvPointsRec->GetEntries();
		cout << "Number of Gammas to be processed: " << nEntriesTrueGamma << endl;
		Long64_t nbytesTrueGam = 0;
		Int_t nComb = 0;
		Int_t nCombElec = 0;
		Int_t nCombPion = 0;
		Int_t nCombPionProton = 0;
		Int_t nCombPionElectron = 0;
		Int_t nCombKaonX = 0;
		Int_t nCombDirectElectron = 0;
		
		for (Long64_t i=0; i<nEntriesTrueGamma;i++) {
			nbytesTrueGam += ConvPointsRec->GetEntry(i);
			if (fMultiplicity09Photon > maxMult || fMultiplicity09Photon < minMult) continue;
			Float_t fR = (*fCoord)[3];
			Float_t fPhi = (*fCoord)[4];
			Float_t fX = (*fCoord)[0];
			Float_t fY = (*fCoord)[1];
			Float_t fZ = (*fCoord)[2];
			Float_t fEta = TMath::Log(TMath::Tan(fTheta/2.));
			Float_t fPosiEta = TMath::Log(TMath::Tan((*fDaugtherProp)[1]/2.));
			Float_t fElecEta = TMath::Log(TMath::Tan((*fDaugtherProp)[3]/2.));      

			if ( fPt > minPtPhotonCut ){
				if (!(fR <= (TMath::Abs(fZ) * fLineCutZRSlope- fLineCutZVtx))){ 
				if ( TMath::Abs(fPosiEta) > etaMaxCut) continue;
				if ( TMath::Abs(fElecEta) > etaMaxCut) continue;
				if ( TMath::Abs(fEta) > etaMaxCut) continue;

					if ( !(etaMinCut != -0.1 &&   fR >= (TMath::Abs(fZ) * fLineCutZRSlopeMin- fLineCutZVtxMin))){
						if (fR > rMin ){ //
							if (fKind == 0 || fKind == 5  ){
								histoTruePt->Fill(fPt);
								histoTrueZR->Fill(fZ,fR);
								histoTrueXY->Fill(fX,fY);
								histoTrueZPhi->Fill(fZ,fPhi);
								histoTrueRPhi->Fill(fR,fPhi);
								histoTrueEta->Fill(fEta);
								if ( fR < rSPD){
									histoTrueZPhi_SPD->Fill(fZ,fPhi);
								} else if (fR < rSPDth){
									histoTrueZPhi_SPDTh->Fill(fZ,fPhi);
								} else if (fR < rSDD){
									histoTrueZPhi_SDD->Fill(fZ,fPhi);
								} else if (fR < rSSD){
									histoTrueZPhi_SSD->Fill(fZ,fPhi);
								} else if (fR < rTPC){
									histoTrueZPhi_ITSTPC->Fill(fZ,fPhi);
								}
								if (fR>rHotZoneMin && fR < rHotZoneMax){
									histoTrueZPhi_HotZone->Fill( fZ,fPhi);
									if (fZ < zHotZoneMax && fZ > zHotZoneMin){
										histoTrueXY_HotZone->Fill( fX,fY);
									}
								}
								if (fZ < zBeamPipeInner && fZ > - zBeamPipeInner) {
									histoTrueXY_Beampipe->Fill( fX,fY);
								}  
								histoTruePt_AftCuts->Fill(fPt);
							}
							if (fKind == 0){
								histoTruePrimZR->Fill(fZ,fR);
								histoTruePrimZPhi->Fill(fZ,fPhi);
								histoTruePrimRPhi->Fill(fR,fPhi);
								histoTruePrimEta->Fill(fEta);
								if ( fR < rSPD){
									histoTruePrimZPhi_SPD->Fill(fZ,fPhi);
								} else if (fR < rSPDth){
									histoTruePrimZPhi_SPDTh->Fill(fZ,fPhi);
								} else if (fR < rSDD){
									histoTruePrimZPhi_SDD->Fill(fZ,fPhi);
								} else if (fR < rSSD){
									histoTruePrimZPhi_SSD->Fill(fZ,fPhi);
								} else if (fR < rTPC){
									histoTruePrimZPhi_ITSTPC->Fill(fZ,fPhi);
								}
								if (fR>rHotZoneMin && fR < rHotZoneMax){
									histoTruePrimZPhi_HotZone->Fill( fZ,fPhi);
								}
							} else if (fKind == 1 || fKind == 10 || fKind == 11 || fKind == 12 || fKind == 13 || fKind == 14 || fKind == 15){
								if (fKind == 1) nComb++;
								if (fKind == 10){
									histoTrueCombElecZR->Fill(fZ,fR);
									histoTrueCombElecEta->Fill(fEta);
									nCombElec++;
								}
								if (fKind == 11){
									histoTrueCombPionZR->Fill(fZ,fR);
									histoTrueCombPionEta->Fill(fEta);
									nCombPion++;
								}
								if (fKind == 12){
									histoTrueCombPionProtonZR->Fill(fZ,fR);
									histoTrueCombPionProtonEta->Fill(fEta);
									nCombPionProton++;
								}
								if (fKind == 13){
									histoTrueCombElecPionZR->Fill(fZ,fR);
									histoTrueCombElecPionEta->Fill(fEta);
									nCombPionElectron++;
								}
								if (fKind == 14){
									histoTrueCombKaonXZR->Fill(fZ,fR);
									histoTrueCombKaonXEta->Fill(fEta);
									nCombKaonX++;
								}
								if (fKind == 15) nCombDirectElectron++;
								
								histoTrueCombZR->Fill(fZ,fR);
								histoTrueCombZPhi->Fill(fZ,fPhi);
								histoTrueCombRPhi->Fill(fR,fPhi);
								histoTrueCombEta->Fill(fEta);
								if ( fR < rSPD){
									histoTrueCombZPhi_SPD->Fill(fZ,fPhi);
								} else if (fR < rSPDth){
									histoTrueCombZPhi_SPDTh->Fill(fZ,fPhi);
								} else if (fR < rSDD){
									histoTrueCombZPhi_SDD->Fill(fZ,fPhi);
								} else if (fR < rSSD){
									histoTrueCombZPhi_SSD->Fill(fZ,fPhi);
								} else if (fR < rTPC){
									histoTrueCombZPhi_ITSTPC->Fill(fZ,fPhi);
								}
								if (fR>rHotZoneMin && fR < rHotZoneMax){
									histoTrueCombZPhi_HotZone->Fill( fZ,fPhi);
								}
							} else if (fKind == 2){
								histoTrueHadZR->Fill(fZ,fR);
								histoTrueHadEta->Fill(fEta);
							} else if (fKind == 3){
								histoTrueDalPi0ZR->Fill(fZ,fR);
								histoTrueDalPi0Eta->Fill(fEta);
								if ( fR < rSPD){
									histoTrueDalPi0ZPhi_SPD->Fill(fZ,fPhi);
								} else if (fR < rSPDth){
									histoTrueDalPi0ZPhi_SPDTh->Fill(fZ,fPhi);
								} else if (fR < rSDD){
									histoTrueDalPi0ZPhi_SDD->Fill(fZ,fPhi);
								} 
							} else if (fKind == 4){
								histoTrueDalEtaZR->Fill(fZ,fR);
								histoTrueDalEtaEta->Fill(fEta);
								if ( fR < rSPD){
									histoTrueDalEtaZPhi_SPD->Fill(fZ,fPhi);
								} else if (fR < rSPDth){
									histoTrueDalEtaZPhi_SPDTh->Fill(fZ,fPhi);
								} else if (fR < rSDD){
									histoTrueDalEtaZPhi_SDD->Fill(fZ,fPhi);
								} 
							}
						}
	//                   }
				}
				}
			}
		}
		TFile* fileMappingDetailedConvTrueWrite = new TFile(fileNameOutput.Data(),"UPDATE");
		
		cout << nComb << "\t" << nCombElec << "\t" << nCombPion << "\t" << nCombPionProton << "\t" << nCombPionElectron <<"\t" << nCombKaonX << "\t" << nCombDirectElectron << endl;
		
		fileMappingDetailedConvTrueWrite->mkdir(Form("GammaConv_%s",  CutSelection.Data()));
		fileMappingDetailedConvTrueWrite->cd(Form("GammaConv_%s",  CutSelection.Data()));
			histoTrueRPhi->Write("ESD_TrueConversionMapping_RPhi",TObject::kWriteDelete);
			histoTrueXY->Write("ESD_TrueConversionMapping_XY",TObject::kWriteDelete);
			histoTrueZR->Write("ESD_TrueConversionMapping_ZR",TObject::kWriteDelete);
			histoTrueZPhi->Write("ESD_TrueConversionMapping_ZPhi",TObject::kWriteDelete);
			histoTrueZPhi_SPD->Write("ESD_TrueConversionMapping_SPD_ZPhi",TObject::kWriteDelete); 
			histoTrueZPhi_SPDTh->Write("ESD_TrueConversionMapping_SPDTh_ZPhi",TObject::kWriteDelete); 
			histoTrueZPhi_SSD->Write("ESD_TrueConversionMapping_SSD_ZPhi",TObject::kWriteDelete); 
			histoTrueZPhi_SDD->Write("ESD_TrueConversionMapping_SDD_ZPhi",TObject::kWriteDelete); 
			histoTrueZPhi_ITSTPC->Write("ESD_TrueConversionMapping_ITSTPC_ZPhi",TObject::kWriteDelete); 
			histoTrueZPhi_HotZone->Write("ESD_TrueConversionMapping_HotZone_ZPhi",TObject::kWriteDelete); 
			histoTrueXY_HotZone->Write("ESD_TrueConversionMapping_HotZone_XY",TObject::kWriteDelete); 
			histoTrueXY_Beampipe->Write("ESD_TrueConversionMappingInnerBeampipe_XY",TObject::kWriteDelete); 
			histoTrueEta->Write("ESD_TrueConversionMapping_Eta",TObject::kWriteDelete);
			histoTruePt->Write("ESD_TrueConvQual_Pt",TObject::kWriteDelete);
			histoTruePt_AftCuts->Write("ESD_TrueConvQualAftCuts_Pt",TObject::kWriteDelete);
			histoTruePrimRPhi->Write("ESD_TruePrimMapping_RPhi",TObject::kWriteDelete);
			histoTruePrimZR->Write("ESD_TruePrimMapping_ZR",TObject::kWriteDelete);
			histoTruePrimZPhi->Write("ESD_TruePrimMapping_ZPhi",TObject::kWriteDelete);
			histoTruePrimZPhi_SPD->Write("ESD_TruePrimMapping_SPD_ZPhi",TObject::kWriteDelete); 
			histoTruePrimZPhi_SPDTh->Write("ESD_TruePrimMapping_SPDTh_ZPhi",TObject::kWriteDelete); 
			histoTruePrimZPhi_SSD->Write("ESD_TruePrimMapping_SSD_ZPhi",TObject::kWriteDelete); 
			histoTruePrimZPhi_SDD->Write("ESD_TruePrimMapping_SDD_ZPhi",TObject::kWriteDelete); 
			histoTruePrimZPhi_ITSTPC->Write("ESD_TruePrimMapping_ITSTPC_ZPhi",TObject::kWriteDelete); 
			histoTruePrimZPhi_HotZone->Write("ESD_TruePrimMapping_HotZone_ZPhi",TObject::kWriteDelete); 
			histoTruePrimEta->Write("ESD_TruePrimMapping_Eta",TObject::kWriteDelete);
			histoTrueCombRPhi->Write("ESD_TrueCombMapping_RPhi",TObject::kWriteDelete);
			histoTrueCombZR->Write("ESD_TrueCombMapping_ZR",TObject::kWriteDelete);
			histoTrueCombElecZR->Write("ESD_TrueCombElecMapping_ZR",TObject::kWriteDelete);
			histoTrueCombPionZR->Write("ESD_TrueCombPionMapping_ZR",TObject::kWriteDelete);
			histoTrueCombElecPionZR->Write("ESD_TrueCombElecPionMapping_ZR",TObject::kWriteDelete);
			histoTrueCombPionProtonZR->Write("ESD_TrueCombPionProtonMapping_ZR",TObject::kWriteDelete);
			histoTrueCombKaonXZR->Write("ESD_TrueCombKaonXMapping_ZR",TObject::kWriteDelete);
			histoTrueCombZPhi->Write("ESD_TrueCombMapping_ZPhi",TObject::kWriteDelete);
			histoTrueCombZPhi_SPD->Write("ESD_TrueCombMapping_SPD_ZPhi",TObject::kWriteDelete); 
			histoTrueCombZPhi_SPDTh->Write("ESD_TrueCombMapping_SPDTh_ZPhi",TObject::kWriteDelete); 
			histoTrueCombZPhi_SSD->Write("ESD_TrueCombMapping_SSD_ZPhi",TObject::kWriteDelete); 
			histoTrueCombZPhi_SDD->Write("ESD_TrueCombMapping_SDD_ZPhi",TObject::kWriteDelete); 
			histoTrueCombZPhi_ITSTPC->Write("ESD_TrueCombMapping_ITSTPC_ZPhi",TObject::kWriteDelete); 
			histoTrueCombZPhi_HotZone->Write("ESD_TrueCombMapping_HotZone_ZPhi",TObject::kWriteDelete); 
			histoTrueCombEta->Write("ESD_TrueCombMapping_Eta",TObject::kWriteDelete);
			histoTrueCombElecEta->Write("ESD_TrueCombElecMapping_Eta",TObject::kWriteDelete);
			histoTrueCombPionEta->Write("ESD_TrueCombPionMapping_Eta",TObject::kWriteDelete);
			histoTrueCombElecPionEta->Write("ESD_TrueCombElecPionMapping_Eta",TObject::kWriteDelete);
			histoTrueCombPionProtonEta->Write("ESD_TrueCombPionProtonMapping_Eta",TObject::kWriteDelete);
			histoTrueCombKaonXEta->Write("ESD_TrueCombKaonXMapping_Eta",TObject::kWriteDelete);

			histoTrueHadZR->Write("ESD_TrueHadMapping_ZR",TObject::kWriteDelete);
			histoTrueHadEta->Write("ESD_TrueHadMapping_Eta",TObject::kWriteDelete);

			histoTrueDalPi0ZR->Write("ESD_TrueDalPi0Mapping_ZR",TObject::kWriteDelete);
			histoTrueDalPi0ZPhi_SPD->Write("ESD_TrueDalPi0Mapping_SPD_ZPhi",TObject::kWriteDelete); 
			histoTrueDalPi0ZPhi_SPDTh->Write("ESD_TrueDalPi0Mapping_SPDTh_ZPhi",TObject::kWriteDelete); 
			histoTrueDalPi0ZPhi_SDD->Write("ESD_TrueDalPi0Mapping_SDD_ZPhi",TObject::kWriteDelete); 
			histoTrueDalPi0Eta->Write("ESD_TrueDalPi0Mapping_Eta",TObject::kWriteDelete);

			histoTrueDalEtaZR->Write("ESD_TrueDalEtaMapping_ZR",TObject::kWriteDelete);
			histoTrueDalEtaZPhi_SPD->Write("ESD_TrueDalEtaMapping_SPD_ZPhi",TObject::kWriteDelete); 
			histoTrueDalEtaZPhi_SPDTh->Write("ESD_TrueDalEtaMapping_SPDTh_ZPhi",TObject::kWriteDelete); 
			histoTrueDalEtaZPhi_SDD->Write("ESD_TrueDalEtaMapping_SDD_ZPhi",TObject::kWriteDelete); 
			histoTrueDalEtaEta->Write("ESD_TrueDalEtaMapping_Eta",TObject::kWriteDelete);
			
		fileMappingDetailedConvTrueWrite->Write();
		fileMappingDetailedConvTrueWrite->Close();
	
		delete histoTrueRPhi;
		delete histoTrueXY;
		delete histoTrueZR;
		delete histoTrueZPhi;
		delete histoTrueZPhi_SPD; 
		delete histoTrueZPhi_SPDTh; 
		delete histoTrueZPhi_SSD; 
		delete histoTrueZPhi_SDD; 
		delete histoTrueZPhi_ITSTPC; 
		delete histoTrueZPhi_HotZone; 
		delete histoTrueXY_HotZone; 
		delete histoTrueXY_Beampipe; 
		delete histoTrueEta;
		delete histoTruePt;
		delete histoTruePt_AftCuts;
		delete histoTruePrimRPhi;
		delete histoTruePrimZR;
		delete histoTruePrimZPhi;
		delete histoTruePrimZPhi_SPD;
		delete histoTruePrimZPhi_SPDTh;
		delete histoTruePrimZPhi_SSD;
		delete histoTruePrimZPhi_SDD;
		delete histoTruePrimZPhi_ITSTPC;
		delete histoTruePrimZPhi_HotZone;
		delete histoTruePrimEta;
		delete histoTrueCombRPhi;
		delete histoTrueCombZR;
		delete histoTrueCombElecZR;
		delete histoTrueCombPionZR;
		delete histoTrueCombElecPionZR;
		delete histoTrueCombPionProtonZR;
		delete histoTrueCombKaonXZR;
		delete histoTrueCombZPhi;
		delete histoTrueCombZPhi_SPD;
		delete histoTrueCombZPhi_SPDTh;
		delete histoTrueCombZPhi_SSD;
		delete histoTrueCombZPhi_SDD;
		delete histoTrueCombZPhi_ITSTPC;
		delete histoTrueCombZPhi_HotZone;
		delete histoTrueCombEta;
		delete histoTrueCombElecEta;
		delete histoTrueCombPionEta;
		delete histoTrueCombElecPionEta;
		delete histoTrueCombPionProtonEta;
		delete histoTrueCombKaonXEta;
		delete histoTrueHadZR;
		delete histoTrueHadEta;
		delete histoTrueDalPi0ZR;
		delete histoTrueDalPi0ZPhi_SPD;
		delete histoTrueDalPi0ZPhi_SPDTh;
		delete histoTrueDalPi0ZPhi_SDD;
		delete histoTrueDalPi0Eta;
		delete histoTrueDalEtaZR;
		delete histoTrueDalEtaZPhi_SPD;
		delete histoTrueDalEtaZPhi_SPDTh;
		delete histoTrueDalEtaZPhi_SDD;
		delete histoTrueDalEtaEta;
	}

}
