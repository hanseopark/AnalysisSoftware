// provided by Gamma Conversion Group, $ALICE_ROOT/PWGGA/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion

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
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
//#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/ExtractSignalBinning.h"
#include "AnalyseDCATest.h"
#include "THnSparse.h"
#include "TTree.h"
                            
                            
                            
// Main Function
void AnalyseDCATestV1(TString meson="", TString fileData="",  TString fileData2="" , TString cutSelection="", TString suffix = "", TString optionEnergy = "", TString optionPeriod = "", Bool_t kMC = kFALSE, Int_t numberOfBins = 10, Int_t mode = 9) {

	gROOT->Reset();

	// Set default plotting styles
	StyleSettingsThesis();
	SetPlotStyle();

	// Reading out cutSelection
   	TString fEventCutSelection				="";
	TString fGammaCutSelection				="";
	TString fClusterCutSelection			="";
	TString fElectronCutSelection			="";
	TString fMesonCutSelection				="";
	fCutSelection 							= cutSelection;
	if (mode == 9){
		ReturnSeparatedCutNumber(fCutSelection, fGammaCutSelection, fElectronCutSelection, fMesonCutSelection);
		fEventCutSelection 					= fGammaCutSelection(0, 7);
		fGammaCutSelection 					= fGammaCutSelection(7, fGammaCutSelection.Length()-7);
		cout << fEventCutSelection.Data() << "\t" << fGammaCutSelection.Data() << endl;
	} else {
		ReturnSeparatedCutNumberAdvanced(cutSelection, fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);
	}

	// Read out centrality string 
	TString intermediate 					= GetCentralityString(fEventCutSelection);
	fTextCent 								= "";
	if (intermediate.CompareTo("pp")==0){
		fTextCent 							= "MinBias";  
		intermediate 						= "";
	} else {
		fTextCent 							= Form("%s central", intermediate.Data());
	}
	// Set energy
	fEnergyText 							= ReturnFullCollisionsSystem(optionEnergy);
	if (optionPeriod.CompareTo("") != 0){
		fEnergyText 						= Form("%s, %s", fEnergyText.Data(), optionPeriod.Data());
	} 
	fEnergyFlag 							= optionEnergy;
		
	// Set MC/Data flag
	if (kMC){
		fMCFlag 							= "MC";
	} else {
		fMCFlag 							= "Data";
	}  
	
	// Set DCAz cut string
	TString fDCAZCut 						= fGammaCutSelection(GetPhotonDcaZPrimVtxCutPosition(fGammaCutSelection), 1);
	fMaxDcaZPhoton 							= AnalyseDCAZPhotonCutValue(fDCAZCut.Atoi());   
	if (fMaxDcaZPhoton> 10) fMaxDcaZPhoton 	= 10.;
	
	// Set meson type
	fMesonType 								= "#pi^{0}";
	if (meson.CompareTo("Eta") == 0) 
		fMesonType 							= "#eta";
	
	// Set date
	fdate 									= ReturnDateString();
	
	// Set Output directory
	TString outputDir 						= Form("%s/%s/%s/%s/AnalyseDCATests", cutSelection.Data(), optionEnergy.Data(), optionPeriod.Data(), suffix.Data());
	gSystem->Exec("mkdir -p "+outputDir);

	// read out files 
	TFile f(fileData.Data());
	TFile f2(fileData2.Data());
	TList *TopDir 							= (TList*)f.Get("GammaConvV1");
	TList *TopDir2 							= (TList*)f2.Get("GammaConvV1");
	if(TopDir == NULL){
		cout<<"ERROR: TopDir not Found"<<endl;
		return;
	}
	TList *HistosGammaConversion 			= (TList*)TopDir->FindObject(Form("Cut Number %s", fCutSelection.Data()));
	TList *HistosGammaConversion2 			= (TList*)TopDir2->FindObject(Form("Cut Number %s", fCutSelection.Data()));
	if(HistosGammaConversion == NULL){
		cout<<"ERROR: " << Form("Cut Number %s",fCutSelection.Data()) << " not Found in File"<<endl;
		return;
	}
	TList *ESDContainer 					= (TList*) HistosGammaConversion2->FindObject(Form("%s ESD histograms", fCutSelection.Data()));
	TList *DCAContainer 					= (TList*) HistosGammaConversion->FindObject(Form("%s Meson DCA tree", fCutSelection.Data()));
	if (!DCAContainer){
		cout<<"ERROR: " << Form("%s Meson DCA tree",fCutSelection.Data()) << " not Found in File"<<endl;
		return;
	}   
	
	TH1D* fEventQuality 					= (TH1D*)ESDContainer->FindObject("NEvents");
	if (optionEnergy.CompareTo("PbPb_2.76TeV") == 0 || optionEnergy.CompareTo("pPb_5.023TeV") == 0){
		fNEvents 							= fEventQuality->GetBinContent(1);
	} else {
		fNEvents 							=  GetNEvents(fEventQuality);
	}
	
	// Initialize histogram arrays & different ranges for fitting, extraction ...
	Initialize(meson.Data(), intermediate, optionPeriod, numberOfBins);
	
	// Set expected mass
	fMesonMassExpect 						= TDatabasePDG::Instance()->GetParticle(fMesonId)->Mass();
	
	// Read DCA tree
	TTree* dcaTree 							= (TTree*)DCAContainer->FindObject("ESD_Mesons_InvMass_Pt_DcazMin_DcazMax_Flag");
	Float_t dcaZMin, dcaZMax, pt, invMass;
	UChar_t quality, mesonMCInfo;
	dcaTree->SetBranchAddress("InvMass",&invMass);
	dcaTree->SetBranchAddress("Pt",&pt);
	dcaTree->SetBranchAddress("DcaZMin",&dcaZMin);
	dcaTree->SetBranchAddress("DcaZMax",&dcaZMax);
	dcaTree->SetBranchAddress("kind",&quality);
	if (kMC) dcaTree->SetBranchAddress("mesonMCInfo",&mesonMCInfo);

	Int_t numberMeson					[6][fNBinsPt];
	Int_t numberMesonFromK0s			[6][fNBinsPt];
	Int_t numberMesonFromEta			[6][fNBinsPt];
	Int_t numberMesonFromSomething		[6][fNBinsPt];
	Int_t numberMesonPrimary			[6][fNBinsPt];
	Int_t numberMesonDalitz				[6][fNBinsPt];
	Int_t numberMesonBackground			[6][fNBinsPt];
	Int_t numberMesonGarbage			[6][fNBinsPt];
	Double_t intTotal					[6][fNBinsPt];
	Double_t intErrTotal				[6][fNBinsPt];
	Double_t intBGHist					[3][6][fNBinsPt];
	Double_t intErrBGHist				[3][6][fNBinsPt];
	Double_t fractionIntHist			[3][6][fNBinsPt];
	Double_t fractionErrIntHist			[3][6][fNBinsPt];
	
	Double_t intTotal_AllCat			[fNBinsPt];
	Double_t intErrTotal_AllCat			[fNBinsPt];
	Double_t intBGHist_AllCat			[fNBinsPt];
	Double_t intErrBGHist_AllCat		[fNBinsPt];
	Double_t intBGFit_AllCat			[fNBinsPt];
	Double_t intErrBGFit_AllCat			[fNBinsPt];
	Double_t fractionIntFit_AllCat		[fNBinsPt];
	Double_t fractionErrIntFit_AllCat	[fNBinsPt];
	Double_t fractionIntHist_AllCat		[fNBinsPt];
	Double_t fractionErrIntHist_AllCat	[fNBinsPt];

	for (Int_t i = 0; i < 6 ; i++){
		for (Int_t j = fStartPtBin; j < fNBinsPt ; j++){
			numberMeson[i][j] 						= 0;
			numberMesonFromK0s[i][j] 				= 0;
			numberMesonFromEta[i][j] 				= 0;
			numberMesonFromSomething[i][j] 			= 0;
			numberMesonPrimary[i][j] 				= 0;
			numberMesonDalitz[i][j] 				= 0;
			numberMesonBackground[i][j] 			= 0;
			numberMesonGarbage[i][j] 				= 0;
			fHistDCAZUnderMeson_MesonPt[i][j] 		= new TH1F(Form("HistDCAZUnderMesonCat_%i_MesonPt_%3.2f-%3.2f", i+1, fBinsPt[j], fBinsPt[j+1]), 
														   Form("HistDCAZUnderMesonCat_%i_MesonPt_%3.2f-%3.2f", i+1, fBinsPt[j], fBinsPt[j+1]), 
															   201,-10,10);
			fHistDCAZUnderMesonBG1_MesonPt[i][j] 	= new TH1F(Form("HistDCAZUnderMesonBG1Cat_%i_MesonPt_%3.2f-%3.2f", i+1, fBinsPt[j], fBinsPt[j+1]),
															   Form("HistDCAZUnderMesonCat_%i_MesonPt_%3.2f-%3.2f", i+1, fBinsPt[j], fBinsPt[j+1]), 
															   201,-10,10);
			fHistDCAZUnderMesonBG2_MesonPt[i][j] 	= new TH1F(Form("HistDCAZUnderMesonBG2Cat_%i_MesonPt_%3.2f-%3.2f", i+1, fBinsPt[j], fBinsPt[j+1]),
															   Form("HistDCAZUnderMesonCat_%i_MesonPt_%3.2f-%3.2f", i+1, fBinsPt[j], fBinsPt[j+1]), 
															   201,-10,10);
			fHistInvMass_MesonPt[i][j] 				= new TH1F(Form("HistInvMassCat_%i_MesonPt_%3.2f-%3.2f", i+1, fBinsPt[j], fBinsPt[j+1]), 
															   Form("HistInvMassCat_%i_MesonPt_%3.2f-%3.2f", i+1, fBinsPt[j], fBinsPt[j+1]), 
															   numberInvMassBins ,fMesonIntRange[0], fMesonIntRange[1]);
	
			if (kMC){
				fHistDCAZTruePrimaryMesonGammaGamma_MesonPt[i][j]			= new TH1F(Form("HistDCAZTruePrimaryMesonGammaGammaCat_%i_MesonPt_%3.2f-%3.2f", i+1, fBinsPt[j], fBinsPt[j+1]), 
																					   Form("HistDCAZTruePrimaryMesonGammaGammaCat_%i_MesonPt_%3.2f-%3.2f", i+1, fBinsPt[j], fBinsPt[j+1]), 
																					   201, -10, 10);
				fHistDCAZTruePrimaryMesonDalitz_MesonPt[i][j] 				= new TH1F(Form("HistDCAZTruePrimaryMesonDalitzCat_%i_MesonPt_%3.2f-%3.2f", i+1, fBinsPt[j], fBinsPt[j+1]),
																					   Form("HistDCAZTruePrimaryMesonDalitzCat_%i_MesonPt_%3.2f-%3.2f", i+1, fBinsPt[j], fBinsPt[j+1]), 
																					   201, -10, 10);
				fHistDCAZTrueSecondaryMesonFromK0s_MesonPt[i][j] 			= new TH1F(Form("HistDCAZTrueSecondaryMesonFromK0sCat_%i_MesonPt_%3.2f-%3.2f", i+1, fBinsPt[j], fBinsPt[j+1]),
																					   Form("HistDCAZTrueSecondaryMesonFromK0sCat_%i_MesonPt_%3.2f-%3.2f", i+1, fBinsPt[j], fBinsPt[j+1]),
																					   201, -10, 10);
				fHistDCAZTrueSecondaryMesonFromEta_MesonPt[i][j] 			= new TH1F(Form("HistDCAZTrueSecondaryMesonFromEtaCat_%i_MesonPt_%3.2f-%3.2f", i+1, fBinsPt[j], fBinsPt[j+1]),
																					   Form("HistDCAZTrueSecondaryMesonFromEtaCat_%i_MesonPt_%3.2f-%3.2f", i+1, fBinsPt[j], fBinsPt[j+1]),
																					   201, -10, 10);
				fHistDCAZTrueSecondaryMesonFromSomething_MesonPt[i][j] 		= new TH1F(Form("HistDCAZTrueSecondaryMesonFromSomethingCat_%i_MesonPt_%3.2f-%3.2f", i+1, fBinsPt[j], fBinsPt[j+1]),
																					   Form("HistDCAZTrueSecondaryMesonFromSomethingCat_%i_MesonPt_%3.2f-%3.2f", i+1, fBinsPt[j], fBinsPt[j+1]),
																					   201, -10, 10);
				fHistDCAZTrueBackground_MesonPt[i][j] 						= new TH1F(Form("HistDCAZTrueBackgroundCat_%i_MesonPt_%3.2f-%3.2f", i+1, fBinsPt[j], fBinsPt[j+1]),
																					   Form("HistDCAZTrueBackgroundCat_%i_MesonPt_%3.2f-%3.2f", i+1, fBinsPt[j], fBinsPt[j+1]), 
																					   201, -10, 10);
				fHistDCAZGarbage_MesonPt[i][j] 								= new TH1F(Form("HistDCAZGarbageCat_%i_MesonPt_%3.2f-%3.2f", i+1, fBinsPt[j], fBinsPt[j+1]), 
																					   Form("HistDCAZGarbageCat_%i_MesonPt_%3.2f-%3.2f", i+1, fBinsPt[j], fBinsPt[j+1]), 
																					   201, -10, 10);
			}   
			
			if (i== 0){
				fHistDCAZUnderMeson_MesonPt_AllCat[j] 						= new TH1F(Form("HistDCAZUnderMesonAllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]),
																					   Form("HistDCAZUnderMesonAllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]), 
																					   201, -10, 10);
				fHistDCAZUnderMesonBG1_MesonPt_AllCat[j] 					= new TH1F(Form("HistDCAZUnderMesonBG1AllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]),
																					   Form("HistDCAZUnderMesonAllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]), 
																					   201, -10, 10);
				fHistDCAZUnderMesonBG2_MesonPt_AllCat[j] 					= new TH1F(Form("HistDCAZUnderMesonBG2AllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]),
																					   Form("HistDCAZUnderMesonAllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]), 
																					   201, -10, 10);
				fHistInvMass_MesonPt_AllCat[j] 								= new TH1F(Form("HistInvMassAllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]),
																					   Form("HistInvMassAllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]),
																					   numberInvMassBins, fMesonIntRange[0], fMesonIntRange[1]);
				fHistDCAZUnderMeson_MesonPt_GoodCat[j] 						= new TH1F(Form("HistDCAZUnderMesonGoodCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]),
																					   Form("HistDCAZUnderMesonGoodCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]),
																					   201, -10, 10);
				
				if (kMC){
					fHistDCAZTruePrimaryMesonGammaGamma_MesonPt_AllCat[j] 		= new TH1F(Form("HistDCAZTruePrimaryMesonGammaGammaAllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]),
																						   Form("HistDCAZTruePrimaryMesonGammaGammaAllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]),
																						   201, -10, 10);
				fHistDCAZTruePrimaryMesonDalitz_MesonPt_AllCat[j] 				= new TH1F(Form("HistDCAZTruePrimaryMesonDalitzAllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]),
																						   Form("HistDCAZTruePrimaryMesonDalitzAllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]), 
																						   201, -10, 10);
				fHistDCAZTrueSecondaryMesonFromK0s_MesonPt_AllCat[j] 			= new TH1F(Form("HistDCAZTrueSecondaryMesonFromK0sAllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]),
																						   Form("HistDCAZTrueSecondaryMesonFromK0sAllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]),
																						   201, -10, 10);
				fHistDCAZTrueSecondaryMesonFromEta_MesonPt_AllCat[j] 			= new TH1F(Form("HistDCAZTrueSecondaryMesonFromEtaAllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]),
																						   Form("HistDCAZTrueSecondaryMesonFromEtaAllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]),
																						   201, -10, 10);
				fHistDCAZTrueSecondaryMesonFromSomething_MesonPt_AllCat[j]	 	= new TH1F(Form("HistDCAZTrueSecondaryMesonFromSomethingAllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]),
																						   Form("HistDCAZTrueSecondaryMesonFromSomethingAllCat_MesonPt_%3.2f-%3.2f",  fBinsPt[j], fBinsPt[j+1]),
																						   201, -10, 10);
				fHistDCAZTrueBackground_MesonPt_AllCat[j] 						= new TH1F(Form("HistDCAZTrueBackgroundAllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]),
																						   Form("HistDCAZTrueBackgroundAllCat_MesonPt_%3.2f-%3.2f",  fBinsPt[j], fBinsPt[j+1]), 
																						   201, -10,10);
				fHistDCAZGarbage_MesonPt_AllCat[j]	 							= new TH1F(Form("HistDCAZGarbageAllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]),
																						   Form("HistDCAZGarbageAllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]), 
																						   201, -10, 10);
				}
			}   
		}
	}   
	TH1F* fHistDCAZUnderMesonAllCat_AllPt						= new TH1F("HistDCAZUnderMesonAllCat_AllPt", "HistDCAZUnderMesonAllCat_AllPt", 
																		   201, -10, 10);
	fHistDCAZUnderMesonAllCat_AllPt->Sumw2();
	TH1F* fHistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt		= new TH1F("HistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt", "HistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt", 
																		   201, -10, 10);
	fHistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt->Sumw2();
	TH1F* fHistDCAZTruePrimaryMesonDalitzAllCat_AllPt			= new TH1F("HistDCAZTruePrimaryMesonDalitzAllCat_AllPt", "HistDCAZTruePrimaryMesonDalitzAllCat_AllPt", 
																		   201, -10, 10);
	fHistDCAZTruePrimaryMesonDalitzAllCat_AllPt->Sumw2();
	TH1F* fHistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt		= new TH1F("HistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt", "HistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt", 
																		   201, -10, 10);
	fHistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt->Sumw2();
	TH1F* fHistDCAZTrueSecondaryMesonFromEtaAllCat_AllPt 		= new TH1F("HistDCAZTrueSecondaryMesonFromEtaAllCat_AllPt", "HistDCAZTrueSecondaryMesonFromEtaAllCat_AllPt", 
																		   201, -10, 10);
	fHistDCAZTrueSecondaryMesonFromEtaAllCat_AllPt->Sumw2();
	TH1F* fHistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt	= new TH1F("HistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt", "HistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt", 
																		   201, -10, 10);
	fHistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt->Sumw2();
	TH1F* fHistDCAZTrueBackgroundAllCat_AllPt					= new TH1F("HistDCAZTrueBackgroundAllCat_AllPt", "HistDCAZTrueBackgroundAllCat_AllPt", 
																		   201, -10, 10);
	fHistDCAZTrueBackgroundAllCat_AllPt->Sumw2();
	TH1F* fHistDCAZGarbageAllCat_AllPt							= new TH1F("HistDCAZGarbageAllCat_AllPt", "HistDCAZGarbageAllCat_AllPt", 
																		   201, -10, 10);
	fHistDCAZGarbageAllCat_AllPt->Sumw2();
	TH1F* fHistDCAZUnderMeson_AllPt						[6];
	TH1F* fHistDCAZTrueBackground_AllPt					[6];
	TH1F* fHistDCAZGarbage_AllPt						[6];
	TH1F* fHistDCAZTruePrimaryMesonDalitz_AllPt			[6];
	TH1F* fHistDCAZTruePrimaryMesonGammaGamma_AllPt		[6];
	TH1F* fHistDCAZTrueSecondaryMesonFromEta_AllPt		[6];
	TH1F* fHistDCAZTrueSecondaryMesonFromK0s_AllPt		[6];
	TH1F* fHistDCAZTrueSecondaryMesonFromSomething_AllPt[6];
	for (Int_t i = 0; i < 6 ; i++){
		fHistDCAZUnderMeson_AllPt[i] 							= new TH1F(Form("HistDCAZUnderMesonCat_%i_AllPt", i+1), Form("HistDCAZUnderMesonCat_%i_AllPt", i+1), 
																		   201,-10,10);
		fHistDCAZUnderMeson_AllPt[i]->Sumw2();
		fHistDCAZTrueBackground_AllPt[i] 						= new TH1F(Form("HistDCAZTrueBackgroundCat_%i_AllPt", i+1), Form("HistDCAZTrueBackgroundCat_%i_AllPt", i+1), 
																		   201,-10,10);
		fHistDCAZTrueBackground_AllPt[i]->Sumw2();
		fHistDCAZGarbage_AllPt[i] 								= new TH1F(Form("HistDCAZGarbageCat_%i_AllPt", i+1), Form("HistDCAZGarbageCat_%i_AllPt", i+1), 
																		   201,-10,10);
		fHistDCAZGarbage_AllPt[i]->Sumw2();
		fHistDCAZTruePrimaryMesonDalitz_AllPt[i] 				= new TH1F(Form("HistDCAZTruePrimaryMesonDalitzCat_%i_AllPt", i+1), Form("HistDCAZTruePrimaryMesonDalitzCat_%i_AllPt", i+1),
																		   201,-10,10);
		fHistDCAZTruePrimaryMesonDalitz_AllPt[i]->Sumw2();
		fHistDCAZTruePrimaryMesonGammaGamma_AllPt[i] 			= new TH1F(Form("HistDCAZTruePrimaryMesonGammaGammaCat_%i_AllPt", i+1), Form("HistDCAZTruePrimaryMesonGammaGammaCat_%i_AllPt", i+1),
																		   201,-10,10);
		fHistDCAZTruePrimaryMesonGammaGamma_AllPt[i]->Sumw2();
		fHistDCAZTrueSecondaryMesonFromEta_AllPt[i] 			= new TH1F(Form("HistDCAZTrueSecondaryMesonFromEtaCat_%i_AllPt", i+1), Form("HistDCAZTrueSecondaryMesonFromEtaCat_%i_AllPt", i+1),
																		   201,-10,10);
		fHistDCAZTrueSecondaryMesonFromEta_AllPt[i]->Sumw2();
		fHistDCAZTrueSecondaryMesonFromK0s_AllPt[i] 			= new TH1F(Form("HistDCAZTrueSecondaryMesonFromK0sCat_%i_AllPt", i+1), Form("HistDCAZTrueSecondaryMesonFromK0sCat_%i_AllPt", i+1),
																		   201,-10,10);
		fHistDCAZTrueSecondaryMesonFromK0s_AllPt[i]->Sumw2();
		fHistDCAZTrueSecondaryMesonFromSomething_AllPt[i] 		= new TH1F(Form("HistDCAZTrueSecondaryMesonFromSomethingCat_%i_AllPt", i+1),
																		   Form("HistDCAZTrueSecondaryMesonFromSomethingCat_%i_AllPt", i+1), 
																		   201,-10,10);
		fHistDCAZTrueSecondaryMesonFromSomething_AllPt[i]->Sumw2();
	}
	

	Long64_t nentries = dcaTree->GetEntries();
	for (Long64_t l=0;l<nentries;l++) {
		dcaTree->GetEntry(l);
		if (invMass > 0 && invMass < 0.8){ 
			for (Int_t i = 0 ; i < 6 ; i++){
				if (quality == i+1){ 
				for (Int_t j = fStartPtBin; j < fNBinsPt ; j++){
					if ( pt > fBinsPt[j] && pt <= fBinsPt[j+1] ){
						if (pt > 2){
							fHistDCAZUnderMeson_MesonPt[i][j]->Fill(dcaZMin);
							fHistDCAZUnderMeson_MesonPt[i][j]->Fill(dcaZMax);
						} else {
							if (invMass > fMesonIntRange[0] && invMass < fMesonIntRange[1]){ 
							fHistDCAZUnderMeson_MesonPt[i][j]->Fill(dcaZMin);
							fHistDCAZUnderMeson_MesonPt[i][j]->Fill(dcaZMax);
							}
						}   
					}
				}
				}
			}
		}
		if (invMass > fMesonIntRange[0] && invMass < fMesonIntRange[1]){ 
			for (Int_t i = 0 ; i < 6 ; i++){
			if (quality == i+1){
				for (Int_t j = fStartPtBin; j < fNBinsPt ; j++){
					if ( pt > fBinsPt[j] && pt <= fBinsPt[j+1] ){
	//                      fHistDCAZUnderMeson_MesonPt[i][j]->Fill(dcaZMin);
	//                      fHistDCAZUnderMeson_MesonPt[i][j]->Fill(dcaZMax);
						fHistInvMass_MesonPt[i][j]->Fill(invMass);
						fHistInvMass_MesonPt_AllCat[j]->Fill(invMass);
						fHistDCAZUnderMeson_MesonPt_AllCat[j]->Fill(dcaZMin);
						fHistDCAZUnderMeson_MesonPt_AllCat[j]->Fill(dcaZMax);
						if ( quality == 4 || quality == 5 || quality == 6){
							fHistDCAZUnderMeson_MesonPt_GoodCat[j]->Fill(dcaZMin);
							fHistDCAZUnderMeson_MesonPt_GoodCat[j]->Fill(dcaZMax);
						}   
						if (kMC){
							if (mesonMCInfo == 6){
								numberMesonPrimary[i][j]++;
								fHistDCAZTruePrimaryMesonGammaGamma_MesonPt[i][j]->Fill(dcaZMin);
								fHistDCAZTruePrimaryMesonGammaGamma_MesonPt[i][j]->Fill(dcaZMax);
								fHistDCAZTruePrimaryMesonGammaGamma_MesonPt_AllCat[j]->Fill(dcaZMin);
								fHistDCAZTruePrimaryMesonGammaGamma_MesonPt_AllCat[j]->Fill(dcaZMax);
							}
							if (mesonMCInfo == 5){
								numberMesonDalitz[i][j]++;
								fHistDCAZTruePrimaryMesonDalitz_MesonPt[i][j]->Fill(dcaZMin);
								fHistDCAZTruePrimaryMesonDalitz_MesonPt[i][j]->Fill(dcaZMax);
								fHistDCAZTruePrimaryMesonDalitz_MesonPt_AllCat[j]->Fill(dcaZMin);
								fHistDCAZTruePrimaryMesonDalitz_MesonPt_AllCat[j]->Fill(dcaZMax);
							}
							if (mesonMCInfo == 4){
								numberMesonFromK0s[i][j]++;
								fHistDCAZTrueSecondaryMesonFromK0s_MesonPt[i][j]->Fill(dcaZMin);
								fHistDCAZTrueSecondaryMesonFromK0s_MesonPt[i][j]->Fill(dcaZMax);
								fHistDCAZTrueSecondaryMesonFromK0s_MesonPt_AllCat[j]->Fill(dcaZMin);
								fHistDCAZTrueSecondaryMesonFromK0s_MesonPt_AllCat[j]->Fill(dcaZMax);
							}
							if (mesonMCInfo == 3){
								numberMesonFromEta[i][j]++;
								fHistDCAZTrueSecondaryMesonFromEta_MesonPt[i][j]->Fill(dcaZMin);
								fHistDCAZTrueSecondaryMesonFromEta_MesonPt[i][j]->Fill(dcaZMax);
								fHistDCAZTrueSecondaryMesonFromEta_MesonPt_AllCat[j]->Fill(dcaZMin);
								fHistDCAZTrueSecondaryMesonFromEta_MesonPt_AllCat[j]->Fill(dcaZMax);
							}
							if (mesonMCInfo == 2){
								numberMesonFromSomething[i][j]++;
								fHistDCAZTrueSecondaryMesonFromSomething_MesonPt[i][j]->Fill(dcaZMin); 
								fHistDCAZTrueSecondaryMesonFromSomething_MesonPt[i][j]->Fill(dcaZMax); 
								fHistDCAZTrueSecondaryMesonFromSomething_MesonPt_AllCat[j]->Fill(dcaZMin); 
								fHistDCAZTrueSecondaryMesonFromSomething_MesonPt_AllCat[j]->Fill(dcaZMax); 
							}
							if (mesonMCInfo == 1){
								numberMesonBackground[i][j]++;
								fHistDCAZTrueBackground_MesonPt[i][j]->Fill(dcaZMin); 
								fHistDCAZTrueBackground_MesonPt[i][j]->Fill(dcaZMax); 
								fHistDCAZTrueBackground_MesonPt_AllCat[j]->Fill(dcaZMin); 
								fHistDCAZTrueBackground_MesonPt_AllCat[j]->Fill(dcaZMax); 
							}
							if (mesonMCInfo == 0){
								numberMesonGarbage[i][j]++;
								fHistDCAZGarbage_MesonPt[i][j]->Fill(dcaZMin); 
								fHistDCAZGarbage_MesonPt[i][j]->Fill(dcaZMax); 
								fHistDCAZGarbage_MesonPt_AllCat[j]->Fill(dcaZMin); 
								fHistDCAZGarbage_MesonPt_AllCat[j]->Fill(dcaZMax); 
							}
						}
						numberMeson[i][j]++;
					}
				}
				}
			}
		}
		if (invMass > fMesonIntRangeBG1[0] && invMass < fMesonIntRangeBG1[1]){ 
			for (Int_t i = 0 ; i < 6 ; i++){
				if (quality == i+1){
					for (Int_t j = fStartPtBin; j < fNBinsPt ; j++){
						if ( pt > fBinsPt[j] && pt <= fBinsPt[j+1] ){
							fHistDCAZUnderMesonBG1_MesonPt[i][j]->Fill(dcaZMin);
							fHistDCAZUnderMesonBG1_MesonPt[i][j]->Fill(dcaZMax);
							fHistDCAZUnderMesonBG1_MesonPt_AllCat[j]->Fill(dcaZMin);
							fHistDCAZUnderMesonBG1_MesonPt_AllCat[j]->Fill(dcaZMax);
						}
					}
				}
			}
		}

		if (invMass > fMesonIntRangeBG2[0] && invMass < fMesonIntRangeBG2[1]){ 
			for (Int_t i = 0 ; i < 6 ; i++){
				if (quality == i+1){
					for (Int_t j = fStartPtBin; j < fNBinsPt ; j++){
						if ( pt > fBinsPt[j] && pt <= fBinsPt[j+1] ){
							fHistDCAZUnderMesonBG2_MesonPt[i][j]->Fill(dcaZMin);
							fHistDCAZUnderMesonBG2_MesonPt[i][j]->Fill(dcaZMax);
							fHistDCAZUnderMesonBG2_MesonPt_AllCat[j]->Fill(dcaZMin);
							fHistDCAZUnderMesonBG2_MesonPt_AllCat[j]->Fill(dcaZMax);
						}
					}
				}
			}
		}
	}


	Int_t mesonsCat					[6] 		= {0, 0, 0, 0, 0, 0};
	Int_t mesonsVsPt				[fNBinsPt];
	Int_t mesonsFromK0sVsPt			[fNBinsPt];
	Int_t mesonsFromEtaVsPt			[fNBinsPt];
	Int_t mesonsFromSomethingVsPt	[fNBinsPt];
	Int_t mesonsPrimaryVsPt			[fNBinsPt];
	Int_t mesonsDalitzVsPt			[fNBinsPt];
	Int_t mesonsBackgroundVsPt		[fNBinsPt];
	Int_t mesonsGarbageVsPt			[fNBinsPt];
	Int_t mesonsAllCat 							= 0;
	
	for (Int_t j = 0; j < fNBinsPt; j++){
		mesonsVsPt[j] 							= 0;
		mesonsFromK0sVsPt[j] 					= 0;
		mesonsFromEtaVsPt[j] 					= 0;
		mesonsFromSomethingVsPt[j] 				= 0;
		mesonsPrimaryVsPt[j] 					= 0;
		mesonsDalitzVsPt[j] 					= 0;
		mesonsBackgroundVsPt[j] 				= 0;
		mesonsGarbageVsPt[j] 					= 0;
	}   
	
	fHoleRadius 								= 1.5;
	TF1* fitWithHole 							= new TF1("fitWithHole", GausWithHole, -fMaxDcaZPhoton, fMaxDcaZPhoton, 3);
	fitWithHole->SetParameter(0,0.001);
	fitWithHole->SetParameter(1,0);
	fitWithHole->SetParameter(2,3);
	fitWithHole->SetParLimits(2,1,6);
		
	TVirtualFitter::SetMaxIterations(130);
	//   TVirtualFitter::SetPrecision(1.e-8);

	cout << "meson statistics" << endl;
	for (Int_t i = 0; i < 6 ; i++){
		cout <<  i+1 << "\t";
		fitWithHole->SetParameter(1,0);
		fitWithHole->SetParameter(2,3);
		fitWithHole->SetParLimits(2,1,6);
		for (Int_t j = fStartPtBin; j < fNBinsPt ; j++){
			if (i == 0){
				if (!fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 && optionPeriod.CompareTo("") == 0 && !meson.Contains("Eta")){
				fitWithHole->SetParameter(0,fHistDCAZUnderMeson_MesonPt_AllCat[j]->GetMaximum()/1e2);
				if (i == 0) fHoleRadius 		= 1.5;
				else fHoleRadius = 2;
					fHistDCAZUnderMeson_MesonPt_AllCat[j]->Fit(fitWithHole,"RME0WLQ");
					fFitDCAZUnderMesonBGEstimate_MesonPt_AllCat[j] 	= new TF1(Form("fFitDCAZUnderMesonBGEstimateAllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]), "gaus",
																			  -fMaxDcaZPhoton,fMaxDcaZPhoton);
					fFitDCAZUnderMesonBGEstimate_MesonPt_AllCat[j]->SetParameter(0,fitWithHole->GetParameter(0));
					fFitDCAZUnderMesonBGEstimate_MesonPt_AllCat[j]->SetParError(0,fitWithHole->GetParError(0));
					fFitDCAZUnderMesonBGEstimate_MesonPt_AllCat[j]->SetParameter(1,fitWithHole->GetParameter(1));
					fFitDCAZUnderMesonBGEstimate_MesonPt_AllCat[j]->SetParError(1,fitWithHole->GetParError(1));
					fFitDCAZUnderMesonBGEstimate_MesonPt_AllCat[j]->SetParameter(2,fitWithHole->GetParameter(2));
					fFitDCAZUnderMesonBGEstimate_MesonPt_AllCat[j]->SetParError(2,fitWithHole->GetParError(2));
				}
				fHistDCAZUnderMesonBGEstimate_MesonPt_AllCat[j] = (TH1F*)fHistDCAZUnderMeson_MesonPt_AllCat[j]->ShowBackground(nIterBGFit,optionBGSmoothingStandard.Data());
				fHistDCAZUnderMesonBGEstimate_MesonPt_AllCat[j]->SetName(Form("fHistDCAZUnderMesonBGEstimateAllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]));
				intBGHist_AllCat[j] 							= fHistDCAZUnderMesonBGEstimate_MesonPt_AllCat[j]->IntegralAndError(fHistDCAZUnderMesonBGEstimate_MesonPt_AllCat[j]->FindBin(
-fMaxDcaZPhoton), fHistDCAZUnderMesonBGEstimate_MesonPt_AllCat[j]->FindBin(fMaxDcaZPhoton), intErrBGHist_AllCat[j], "width");
				cout << j << "\t" << intBGHist_AllCat[j] << endl;
				if (!fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 && optionPeriod.CompareTo("") == 0 && !meson.Contains("Eta")){
					Double_t integrationWindowLow 					= fHistDCAZUnderMesonBGEstimate_MesonPt_AllCat[j]->GetBinLowEdge(fHistDCAZUnderMesonBGEstimate_MesonPt_AllCat[j]->FindBin(-fMaxDcaZPhoton));
					Double_t integrationWindowUp 					= fHistDCAZUnderMesonBGEstimate_MesonPt_AllCat[j]->GetXaxis()->GetBinUpEdge(fHistDCAZUnderMesonBGEstimate_MesonPt_AllCat[j]->FindBin(fMaxDcaZPhoton));
				intBGFit_AllCat[j] 									= fFitDCAZUnderMesonBGEstimate_MesonPt_AllCat[j]->Integral(integrationWindowLow,integrationWindowUp);
					if (kMC){
						intErrBGFit_AllCat[j] 						= TMath::Sqrt(intBGFit_AllCat[j]);
					} else {
						intErrBGFit_AllCat[j] 						= fFitDCAZUnderMesonBGEstimate_MesonPt_AllCat[j]->IntegralError(integrationWindowLow,integrationWindowUp);
					}   
				}
				intTotal_AllCat[j] 									= fHistDCAZUnderMeson_MesonPt_AllCat[j]->IntegralAndError(fHistDCAZUnderMeson_MesonPt_AllCat[j]->FindBin(-fMaxDcaZPhoton),fHistDCAZUnderMeson_MesonPt_AllCat[j]->FindBin(fMaxDcaZPhoton),intErrTotal_AllCat[j],"width");
				if (!fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 && optionPeriod.CompareTo("") == 0 && !meson.Contains("Eta")){
				if (intTotal_AllCat[j] != 0 && intBGFit_AllCat[j]!= 0 ){
					fractionIntFit_AllCat[j] 						= intBGFit_AllCat[j]/intTotal_AllCat[j];
					fractionErrIntFit_AllCat[j] 					= TMath::Sqrt(TMath::Power(intErrBGFit_AllCat[j]/intBGFit_AllCat[j],2) +
																				  TMath::Power(intErrTotal_AllCat[j]/intTotal_AllCat[j],2))*fractionIntFit_AllCat[j];
				}  else {
					fractionIntFit_AllCat[j] 						= 0.;   
					fractionErrIntFit_AllCat[j] 					= 0.;
				}
				
				}
				if (intTotal_AllCat[j] != 0 && intBGHist_AllCat[j]!= 0){
					fractionIntHist_AllCat[j] 						= intBGHist_AllCat[j]/intTotal_AllCat[j];
					fractionErrIntHist_AllCat[j] 					= TMath::Sqrt(TMath::Power(intErrBGHist_AllCat[j]/intBGHist_AllCat[j],2) +
																			      TMath::Power(intErrTotal_AllCat[j]/intTotal_AllCat[j],2))*fractionIntHist_AllCat[j];
				} else {
					fractionIntHist_AllCat[j] 						= 0.;
					fractionErrIntHist_AllCat[j] 					= 0.;
				}  
			}   
			mesonsCat[i] 					= mesonsCat[i]+numberMeson[i][j] ;
			mesonsVsPt[j] 					= mesonsVsPt[j]+numberMeson[i][j];
			mesonsAllCat 					= mesonsAllCat+numberMeson[i][j] ;
			fHistDCAZUnderMeson_AllPt[i]->Add(fHistDCAZUnderMeson_MesonPt[i][j]);
			fHistDCAZUnderMesonAllCat_AllPt->Add(fHistDCAZUnderMeson_MesonPt[i][j]);
			if (kMC){
				fHistDCAZGarbage_AllPt[i]->Add(fHistDCAZGarbage_MesonPt[i][j]);
				fHistDCAZTrueBackground_AllPt[i]->Add(fHistDCAZTrueBackground_MesonPt[i][j]);
				fHistDCAZTrueSecondaryMesonFromSomething_AllPt[i]->Add(fHistDCAZTrueSecondaryMesonFromSomething_MesonPt[i][j]);
				fHistDCAZTrueSecondaryMesonFromEta_AllPt[i]->Add(fHistDCAZTrueSecondaryMesonFromEta_MesonPt[i][j]);
				fHistDCAZTrueSecondaryMesonFromK0s_AllPt[i]->Add(fHistDCAZTrueSecondaryMesonFromK0s_MesonPt[i][j]);
				fHistDCAZTruePrimaryMesonDalitz_AllPt[i]->Add(fHistDCAZTruePrimaryMesonDalitz_MesonPt[i][j]);
				fHistDCAZTruePrimaryMesonGammaGamma_AllPt[i]->Add(fHistDCAZTruePrimaryMesonGammaGamma_MesonPt[i][j]);
				fHistDCAZGarbageAllCat_AllPt->Add(fHistDCAZGarbage_MesonPt[i][j]);
				fHistDCAZTrueBackgroundAllCat_AllPt->Add(fHistDCAZTrueBackground_MesonPt[i][j]);
				fHistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt->Add(fHistDCAZTrueSecondaryMesonFromSomething_MesonPt[i][j]);
				fHistDCAZTrueSecondaryMesonFromEtaAllCat_AllPt->Add(fHistDCAZTrueSecondaryMesonFromEta_MesonPt[i][j]);
				fHistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt->Add(fHistDCAZTrueSecondaryMesonFromK0s_MesonPt[i][j]);
				fHistDCAZTruePrimaryMesonDalitzAllCat_AllPt->Add(fHistDCAZTruePrimaryMesonDalitz_MesonPt[i][j]);
				fHistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt->Add(fHistDCAZTruePrimaryMesonGammaGamma_MesonPt[i][j]);
				
			}
			
			intTotal[i][j] 					= fHistDCAZUnderMeson_MesonPt[i][j]->IntegralAndError(fHistDCAZUnderMeson_MesonPt[i][j]->FindBin(-fMaxDcaZPhoton),fHistDCAZUnderMeson_MesonPt[i][j]->FindBin(fMaxDcaZPhoton),intErrTotal[i][j],"width");
			for (Int_t k = 0; k < 3; k ++){
				if ( k == 0){
					if (i == 1) fHistDCAZUnderMesonBGEstimate_MesonPt[k][i][j] 	= (TH1F*)fHistDCAZUnderMeson_MesonPt[i][j]->ShowBackground(nIterBGFit+2,optionBGSmoothingStandard.Data());
					else fHistDCAZUnderMesonBGEstimate_MesonPt[k][i][j] 		= (TH1F*)fHistDCAZUnderMeson_MesonPt[i][j]->ShowBackground(nIterBGFit,optionBGSmoothingStandard.Data());
				} else if (k == 1){
					if (i == 1) fHistDCAZUnderMesonBGEstimate_MesonPt[k][i][j] 	= (TH1F*)fHistDCAZUnderMeson_MesonPt[i][j]->ShowBackground(nIterBGFit+2,optionBGSmoothingVar2.Data());
					else fHistDCAZUnderMesonBGEstimate_MesonPt[k][i][j] 		= (TH1F*)fHistDCAZUnderMeson_MesonPt[i][j]->ShowBackground(nIterBGFit,optionBGSmoothingVar2.Data());
				} else if (k == 2){
					if (i == 1) fHistDCAZUnderMesonBGEstimate_MesonPt[k][i][j] 	= (TH1F*)fHistDCAZUnderMeson_MesonPt[i][j]->ShowBackground(nIterBGFit+2,optionBGSmoothingVar1.Data());
					else fHistDCAZUnderMesonBGEstimate_MesonPt[k][i][j] 		= (TH1F*)fHistDCAZUnderMeson_MesonPt[i][j]->ShowBackground(nIterBGFit,optionBGSmoothingVar1.Data());
				}   
				fHistDCAZUnderMesonBGEstimate_MesonPt[k][i][j]->SetName(Form("fHistDCAZUnderMesonBGEstimateCat_%i_MesonPt_%3.2f-%3.2f",i+1, fBinsPt[j], fBinsPt[j+1]));
				intBGHist[k][i][j] 			= fHistDCAZUnderMesonBGEstimate_MesonPt[k][i][j]->IntegralAndError(fHistDCAZUnderMesonBGEstimate_MesonPt[k][i][j]->FindBin(-fMaxDcaZPhoton),
																											   fHistDCAZUnderMesonBGEstimate_MesonPt[k][i][j]->FindBin(fMaxDcaZPhoton),
																											   intErrBGHist[k][i][j], "width");
			
				if (intTotal[i][j] != 0 && intBGHist[k][i][j]!= 0){
					fractionIntHist[k][i][j] 	= intBGHist[k][i][j]/intTotal[i][j];
					fractionErrIntHist[k][i][j] = TMath::Sqrt(TMath::Power(intErrBGHist[k][i][j]/intBGHist[k][i][j],2) + TMath::Power(intErrTotal[i][j]/intTotal[i][j],2))*fractionIntHist[k][i][j];
				} else {
					fractionIntHist[k][i][j] 	= 0.;
					fractionErrIntHist[k][i][j] = 0.;
				}
				
			}   
		}
	}   
	
	
	if (kMC){
		cout << "primary mesons statistics" << endl;
		for (Int_t i = 0; i < 6 ; i++){
		cout << i+1 << "\t";
		for (Int_t j = fStartPtBin; j < fNBinsPt ; j++){
			cout << numberMesonPrimary[i][j] << "\t" ;
			mesonsPrimaryVsPt[j] = mesonsPrimaryVsPt[j]+numberMesonPrimary[i][j];
		}
		cout<< endl;
		}   
		cout << "dalitz mesons statistics" << endl;
		for (Int_t i = 0; i < 6 ; i++){
		cout << i+1 << "\t";
		for (Int_t j = fStartPtBin; j < fNBinsPt ; j++){
			cout << numberMesonDalitz[i][j] << "\t" ;
			mesonsDalitzVsPt[j] = mesonsDalitzVsPt[j]+numberMesonDalitz[i][j];
		}
		cout<< endl;
		}
		cout << "secondary from K0s statistics" << endl;
		for (Int_t i = 0; i < 6 ; i++){
		cout << i+1 << "\t";
		for (Int_t j = fStartPtBin; j < fNBinsPt ; j++){
			cout << numberMesonFromK0s[i][j] << "\t" ;
			mesonsFromK0sVsPt[j] = mesonsFromK0sVsPt[j]+numberMesonFromK0s[i][j];
		}
		cout<< endl;
		}
		cout << "secondary from Eta statistics" << endl;
		for (Int_t i = 0; i < 6 ; i++){
		cout << i+1 << "\t";
		for (Int_t j = fStartPtBin; j < fNBinsPt ; j++){
			cout << numberMesonFromEta[i][j] << "\t" ;
			mesonsFromEtaVsPt[j] = mesonsFromEtaVsPt[j]+numberMesonFromEta[i][j];
		}
		cout<< endl;
		}
		cout << "secondary from something" << endl;
		for (Int_t i = 0; i < 6 ; i++){
		cout << i+1 << "\t";
		for (Int_t j = fStartPtBin; j < fNBinsPt ; j++){
			cout << numberMesonFromSomething[i][j] << "\t" ;
			mesonsFromSomethingVsPt[j] = mesonsFromSomethingVsPt[j]+numberMesonFromSomething[i][j];
		}
		cout<< endl;
		}
		cout << "gamma-gamma BG statistics" << endl;
		for (Int_t i = 0; i < 6 ; i++){
		cout << i+1 << "\t";
		for (Int_t j = fStartPtBin; j < fNBinsPt ; j++){
			cout << numberMesonBackground[i][j] << "\t" ;
			mesonsBackgroundVsPt[j] = mesonsBackgroundVsPt[j]+numberMesonBackground[i][j];
		}
		cout<< endl;
		}
		cout << "garbage statistics" << endl;
		for (Int_t i = 0; i < 6 ; i++){
		cout << i+1 << "\t";
		for (Int_t j = fStartPtBin; j < fNBinsPt ; j++){
			cout << numberMesonGarbage[i][j] << "\t" ;
			mesonsGarbageVsPt[j] = mesonsGarbageVsPt[j]+numberMesonGarbage[i][j];
		}
		cout<< endl;
		}
	}  
	
	cout << "fraction from primary: \t"; 
	for (Int_t j = fStartPtBin; j < fNBinsPt ; j++){
		cout << (Double_t)mesonsPrimaryVsPt[j]/(Double_t)mesonsVsPt[j] *100 << "\t";
	}   
	cout << endl;
	cout << "fraction from dalitz: \t"; 
	for (Int_t j = fStartPtBin; j < fNBinsPt ; j++){
		cout << (Double_t)mesonsDalitzVsPt[j]/(Double_t)mesonsVsPt[j] *100 << "\t";
	}   
	cout << endl;
	cout << "fraction from K0s: \t"; 
	for (Int_t j = fStartPtBin; j < fNBinsPt ; j++){
		cout << (Double_t)mesonsFromK0sVsPt[j]/(Double_t)mesonsVsPt[j] *100 << "\t";
	}   
	cout << endl; 
	cout << "fraction from Eta: \t"; 
	for (Int_t j = fStartPtBin; j < fNBinsPt ; j++){
		cout << (Double_t)mesonsFromEtaVsPt[j]/(Double_t)mesonsVsPt[j] *100 << "\t";
	}   
	cout << endl;
	cout << "fraction from Some.: \t"; 
	for (Int_t j = fStartPtBin; j < fNBinsPt ; j++){
		cout << (Double_t)mesonsFromSomethingVsPt[j]/(Double_t)mesonsVsPt[j] *100 << "\t";
	}   
	cout << endl;
	cout << "fraction from BG: \t"; 
	for (Int_t j = fStartPtBin; j < fNBinsPt ; j++){
		cout << (Double_t)mesonsBackgroundVsPt[j]/(Double_t)mesonsVsPt[j] *100 << "\t";
	}   
	cout << endl;
	cout << "fraction from Garbage: \t"; 
	for (Int_t j = fStartPtBin; j < fNBinsPt ; j++){
		cout << (Double_t)mesonsGarbageVsPt[j]/(Double_t)mesonsVsPt[j] *100 << "\t";
	}   
	cout << endl;
	cout << "fraction from dalitz/Primary Pi0: \t"; 
	for (Int_t j = fStartPtBin; j < fNBinsPt ; j++){
		cout << (Double_t)mesonsDalitzVsPt[j]/(Double_t)mesonsPrimaryVsPt[j] *100 << "\t";
	}   
	cout << endl; 
	
	
	cout << endl << endl << "all pT -- Category 1: \t" << (Double_t)mesonsCat[0]/(Double_t)mesonsAllCat*100 << "\t Category 2: \t" <<  (Double_t)mesonsCat[1]/(Double_t)mesonsAllCat*100 << "\t Category 3: \t" <<  (Double_t)mesonsCat[2]/(Double_t)mesonsAllCat*100 << "\t Category 4: \t" <<  (Double_t)mesonsCat[3]/(Double_t)mesonsAllCat*100<< "\t Category 5: \t" <<  (Double_t)mesonsCat[4]/(Double_t)mesonsAllCat*100 << "\t Category 6: \t" <<  (Double_t)mesonsCat[5]/(Double_t)mesonsAllCat*100<< endl << endl;
	

	PlotDCADistPtBinWithFitAndEstimateCat(Form("%s/%s_%s_DCAProjections_Cat_%i.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), 1, suffix.Data()),
										  "canvas_cat1", "pad_cat1", fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, kMC , fTextCent, fdate , 1);
	PlotDCADistPtBinWithFitAndEstimateCat(Form("%s/%s_%s_DCAProjections_Cat_%i.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), 2, suffix.Data()), 
										  "canvas_cat2", "pad_cat2", fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, kMC , fTextCent, fdate , 2);
	PlotDCADistPtBinWithFitAndEstimateCat(Form("%s/%s_%s_DCAProjections_Cat_%i.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), 3, suffix.Data()), 
										  "canvas_cat3", "pad_cat3",   fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, kMC , fTextCent, fdate , 3);
	PlotDCADistPtBinWithFitAndEstimateCat(Form("%s/%s_%s_DCAProjections_Cat_%i.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), 4, suffix.Data()), 
										  "canvas_cat4", "pad_cat4",   fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, kMC , fTextCent, fdate , 4);
	PlotDCADistPtBinWithFitAndEstimateCat(Form("%s/%s_%s_DCAProjections_Cat_%i.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), 5, suffix.Data()), 
										  "canvas_cat5", "pad_cat5",   fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, kMC , fTextCent, fdate , 5);
	PlotDCADistPtBinWithFitAndEstimateCat(Form("%s/%s_%s_DCAProjections_Cat_%i.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), 6, suffix.Data()),
										  "canvas_cat6", "pad_cat6",   fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, kMC , fTextCent, fdate , 6);
	PlotDCADistPtBinWithFitAndEstimate(Form("%s/%s_%s_DCAProjections_AllCat.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), suffix.Data()), 
									   "canvas_Allcat", "pad_Allcat",   fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, kMC , fTextCent, fdate);
	PlotDCADistPtBinWithFitAndEstimateGoodCat(Form("%s/%s_%s_DCAProjections_GoodCat.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), suffix.Data()), 
											  "canvas_Allcat", "pad_Allcat",   fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, kMC , fTextCent, fdate);
	
	PlotInvMassPtBinCat(Form("%s/%s_%s_InvMass_Cat_%i.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), 1, suffix.Data()),
						"canvas_cat1", "pad_cat1",   fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, kMC , fTextCent, fdate , 1);
	PlotInvMassPtBinCat(Form("%s/%s_%s_InvMass_Cat_%i.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), 2, suffix.Data()), 
						"canvas_cat1", "pad_cat1",   fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, kMC , fTextCent, fdate , 2);
	PlotInvMassPtBinCat(Form("%s/%s_%s_InvMass_Cat_%i.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), 3, suffix.Data()),
						"canvas_cat1", "pad_cat1",   fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, kMC , fTextCent, fdate , 3);
	PlotInvMassPtBinCat(Form("%s/%s_%s_InvMass_Cat_%i.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), 4, suffix.Data()), 
						"canvas_cat1", "pad_cat1",   fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, kMC , fTextCent, fdate , 4);
	PlotInvMassPtBinCat(Form("%s/%s_%s_InvMass_Cat_%i.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), 5, suffix.Data()), 
						"canvas_cat1", "pad_cat1",   fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, kMC , fTextCent, fdate , 5);
	PlotInvMassPtBinCat(Form("%s/%s_%s_InvMass_Cat_%i.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), 6, suffix.Data()), 
						"canvas_cat1", "pad_cat1",   fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, kMC , fTextCent, fdate , 6);
	PlotInvMassPtBin(Form("%s/%s_%s_InvMass.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), suffix.Data()),
					 "canvas_cat1", "pad_cat1",   fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, kMC , fTextCent, fdate);
	
	if (kMC){
		PlotDCADistPtBinWithMCSplitCat(Form("%s/%s_%s_DCAProjectionsWithMC_Cat_%i.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), 1, suffix.Data()), 
									   "canvas_cat1", "pad_cat1",   fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, kMC , fTextCent, fdate , 1);
		PlotDCADistPtBinWithMCSplitCat(Form("%s/%s_%s_DCAProjectionsWithMC_Cat_%i.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), 2, suffix.Data()), 
									   "canvas_cat2", "pad_cat2",   fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, kMC , fTextCent, fdate , 2);
		PlotDCADistPtBinWithMCSplitCat(Form("%s/%s_%s_DCAProjectionsWithMC_Cat_%i.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), 3, suffix.Data()), 
									   "canvas_cat3", "pad_cat3",   fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, kMC , fTextCent, fdate , 3);
		PlotDCADistPtBinWithMCSplitCat(Form("%s/%s_%s_DCAProjectionsWithMC_Cat_%i.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), 4, suffix.Data()), 
									   "canvas_cat4", "pad_cat4",   fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, kMC , fTextCent, fdate , 4);
		PlotDCADistPtBinWithMCSplitCat(Form("%s/%s_%s_DCAProjectionsWithMC_Cat_%i.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), 5, suffix.Data()), 
									   "canvas_cat5", "pad_cat5",   fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, kMC , fTextCent, fdate , 5);
		PlotDCADistPtBinWithMCSplitCat(Form("%s/%s_%s_DCAProjectionsWithMC_Cat_%i.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), 6, suffix.Data()), 
									   "canvas_cat6", "pad_cat6",   fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, kMC , fTextCent, fdate , 6);
		PlotDCADistPtBinWithMCSplit(Form("%s/%s_%s_DCAProjectionsWithMC_AllCat.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), suffix.Data()),
									"canvas_cat6", "pad_cat6",   fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, kMC , fTextCent, fdate );
	}
	
	TH1D* fHistFracIntHistBGvsPt				[3][6] ;
	TH1D* fHistFracIntHistBGFittedvsPt			[3][6] ;
	TH1D* fHistFracCatvsPt						[6] ;	
	TH1D* fHistCorrectionFactorsHistCat			[3][6];
	TH1D* fHistCorrectionFactorsHistCatFitted	[6];
	TH1D* fHistCorrectionFactorsFitAllCat		= new TH1D(Form("fHistCorrectionFactorsFitAllCat_vsPt"),"",fNBinsPt,fBinsPt);
	TH1D* fHistCorrectionFactorsHistAllCat		= new TH1D(Form("fHistCorrectionFactorsHistAllCat_vsPt"),"",fNBinsPt,fBinsPt);
	for (Int_t i = 0; i< 6; i++){
		fHistFracCatvsPt[i]						= new TH1D(Form("fHistFracCat_%i_vsPt",i+1),"",fNBinsPt,fBinsPt);
		for (Int_t k = 0; k < 3; k++){
			fHistFracIntHistBGvsPt[k][i]=  new TH1D(Form("fHistFracIntHistBGvsPt_Cat_%i_Variant_%i",i+1,k+1),"",fNBinsPt,fBinsPt);
			fHistFracIntHistBGFittedvsPt[k][i]	= new TH1D(Form("fHistFracIntHistBGFittedvsPt_Cat_%i_Variant_%i",i+1,k+1),"",fNBinsPt,fBinsPt);
		}
		for (Int_t j = fStartPtBin; j < fNBinsPt; j++){
			if (i == 0){
				if (!fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 && optionPeriod.CompareTo("") == 0 && !meson.Contains("Eta")){
				if (isfinite(fractionErrIntFit_AllCat[j]) && isfinite(fractionIntFit_AllCat[j])){
					fHistCorrectionFactorsFitAllCat->SetBinContent(j+1,fractionIntFit_AllCat[j]*100);
					fHistCorrectionFactorsFitAllCat->SetBinError(j+1,fractionErrIntFit_AllCat[j]*100);
				} else {
					fHistCorrectionFactorsFitAllCat->SetBinContent(j+1,0.);
					fHistCorrectionFactorsFitAllCat->SetBinError(j+1,0.);
				}
				}
				if (isfinite(fractionErrIntHist_AllCat[j]) && isfinite(fractionIntHist_AllCat[j])){
					fHistCorrectionFactorsHistAllCat->SetBinContent(j+1,fractionIntHist_AllCat[j]*100);
					fHistCorrectionFactorsHistAllCat->SetBinError(j+1,fractionErrIntHist_AllCat[j]*100);
				} else {
					fHistCorrectionFactorsHistAllCat->SetBinContent(j+1,0.);
					fHistCorrectionFactorsHistAllCat->SetBinError(j+1,0.);
				}
			}
			if (isfinite((Double_t)numberMeson[i][j]/mesonsVsPt[j]) && numberMeson[i][j] > 0){
				fHistFracCatvsPt[i]->SetBinContent(j+1,(Double_t)numberMeson[i][j]/mesonsVsPt[j]*100);
				Double_t error 					= TMath::Sqrt(TMath::Power(TMath::Sqrt((Double_t)numberMeson[i][j])/(Double_t)numberMeson[i][j],2) +
															  TMath::Power(TMath::Sqrt((Double_t)mesonsVsPt[j])/mesonsVsPt[j],2))*(Double_t)numberMeson[i][j]/mesonsVsPt[j];
				fHistFracCatvsPt[i]->SetBinError(j+1,error*100);
			} else {
				fHistFracCatvsPt[i]->SetBinError(j+1,0.);
				fHistFracCatvsPt[i]->SetBinContent(j+1,0.);
			}
			for (Int_t k = 0; k < 3; k++){   
				if (isfinite(fractionIntHist[k][i][j])){
					fHistFracIntHistBGvsPt[k][i]->SetBinContent(j+1,fractionIntHist[k][i][j]*100);
				} else {
					fHistFracIntHistBGvsPt[k][i]->SetBinContent(j+1,0.);
				}   
				if (isfinite(fractionErrIntHist[k][i][j])){
					fHistFracIntHistBGvsPt[k][i]->SetBinError(j+1,fractionErrIntHist[k][i][j]*100);
				} else {
					fHistFracIntHistBGvsPt[k][i]->SetBinError(j+1,0.);
				}      
			}   
		}   
	}
	
	Color_t  colorCat[6]    	= { kRed+1, 807, 800, kGreen+2, kCyan+2, kBlue+1};
	Style_t  styleCat[6]    	= { kFullCircle, kFullSquare, kFullStar, kFullDiamond, kFullCircle, kFullCross};
	
	TF1* fitFracIntHistBGvsPt[6];
	TCanvas* canvasCorrFrac 	= new TCanvas("canvasCorrFrac","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasCorrFrac, 0.08, 0.02, 0.02, 0.09);
	
	canvasCorrFrac->cd();
		DrawAutoGammaMesonHistos( fHistFracIntHistBGvsPt[0][0], 
								  "", Form("p_{T,%s} (GeV/c)",fMesonType.Data()), "BG/Total (%)", 
								  kFALSE, 2.5,1e-8, kFALSE,
								  kTRUE, 0, fMaxYFracBGOverIntHist, 
								  kFALSE, 0., 10.);
		fHistFracIntHistBGvsPt[0][0]->GetYaxis()->SetTitleOffset(0.8);
			
		TLegend* legendFracIntHistBG = new TLegend(0.65,0.7,0.85,0.95);
		legendFracIntHistBG->SetTextSize(0.04);         
		legendFracIntHistBG->SetFillColor(0);
		legendFracIntHistBG->SetLineColor(0);  
		
		for (Int_t j = 0; j< 3; j++){
			for (Int_t i = 0; i< 3; i++){
				DrawGammaSetMarker(fHistFracIntHistBGvsPt[j][i], styleCat[i], 1., colorCat[i], colorCat[i]);
				fitFracIntHistBGvsPt[i] 				= new TF1(Form("fitFracIntHistBGvsPt_%i",i),"[0]/pow(x,[1])");
				fitFracIntHistBGvsPt[i]->SetRange(fBinsPt[fStartPtBin], fBinsPt[fNBinsPt]);
				TFitResultPtr resultFracIntHistBGvsPt1 	= fHistFracIntHistBGvsPt[j][i]->Fit(fitFracIntHistBGvsPt[i],"SINRME+","",fBinsPt[fStartPtBin],4);
				fitFracIntHistBGvsPt[i]->SetLineColor(colorCat[i]);
				
				for (Int_t k = 2; k < fHistFracIntHistBGFittedvsPt[j][i]->GetNbinsX()+1; k++){
					Double_t ptStart 					= fHistFracIntHistBGFittedvsPt[j][i]->GetXaxis()->GetBinLowEdge(k);
					Double_t ptEnd 						= fHistFracIntHistBGFittedvsPt[j][i]->GetXaxis()->GetBinUpEdge(k);
					Double_t binWidth 					= ptEnd-ptStart;
					Double_t bgEstimate 				= (fitFracIntHistBGvsPt[i]->Integral(ptStart, ptEnd, resultFracIntHistBGvsPt1->GetParams()) / binWidth );
					Double_t errorBGEstimate		 	= (fitFracIntHistBGvsPt[i]->IntegralError(ptStart, ptEnd, resultFracIntHistBGvsPt1->GetParams(), 
																								  resultFracIntHistBGvsPt1->GetCovarianceMatrix().GetMatrixArray() ) / binWidth );
					fHistFracIntHistBGFittedvsPt[j][i]->SetBinContent(k, bgEstimate);
					fHistFracIntHistBGFittedvsPt[j][i]->SetBinError(k, errorBGEstimate);
					if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 && ptEnd> 4.) {
						fHistFracIntHistBGvsPt[j][i]->SetBinContent(k, bgEstimate);
						fHistFracIntHistBGvsPt[j][i]->SetBinError(k, errorBGEstimate);
					}
				}
				if (i == 0 && j==0)fHistFracIntHistBGvsPt[j][i]->DrawCopy("p,e1"); 
				else if (j == 0)fHistFracIntHistBGvsPt[j][i]->DrawCopy("same,p,e1"); 
				if (j == 0 && fEnergyFlag.CompareTo("pPb_5.023TeV") != 0 )fitFracIntHistBGvsPt[i]->Draw("same");
				if (j == 0)legendFracIntHistBG->AddEntry(fHistFracIntHistBGvsPt[j][i],Form("Category %i",i+1),"p");
			}   
		}   
		legendFracIntHistBG->Draw();
		
		TLatex *labelEnergy = new TLatex(0.11,0.9,Form("%s %s",intermediate.Data(),fEnergyText.Data()));
		SetStyleTLatex( labelEnergy, 0.04,4);
		labelEnergy->Draw();
	canvasCorrFrac->Update(); 
	canvasCorrFrac->SaveAs(Form("%s/%s_%s_FracBGOverIntHist.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), suffix.Data()));

	for (Int_t i = 0; i< 6; i++){
		for (Int_t k = 0; k < 3; k++){   
			fHistCorrectionFactorsHistCat[k][i] = (TH1D*)fHistFracCatvsPt[i]->Clone(Form("fHistCorrectionFactorsHistCat_%i_Variant_%i_vsPt",i+1,k+1));
			fHistCorrectionFactorsHistCat[k][i]->Sumw2();
			fHistCorrectionFactorsHistCat[k][i]->Multiply(fHistFracIntHistBGvsPt[k][i]);
			fHistCorrectionFactorsHistCat[k][i]->Scale(1./100.);
		}
		fHistCorrectionFactorsHistCatFitted[i] 	= (TH1D*)fHistFracCatvsPt[i]->Clone(Form("fHistCorrectionFactorsHistCat_%i_Variant_vsPt",i+1));
		fHistCorrectionFactorsHistCatFitted[i]->Sumw2();
		fHistCorrectionFactorsHistCatFitted[i]->Multiply(fHistFracIntHistBGFittedvsPt[0][i]);
		fHistCorrectionFactorsHistCatFitted[i]->Scale(1./100.);
	}

	
	TH1D* fHistCorrectionFactorsHist[3];
	TH1D* fHistCorrectionFactorsHistFitted 		= new TH1D(Form("fHistCorrectionFactorsHistvsPtFitted"),"",fNBinsPt,fBinsPt);
	for (Int_t i = 0; i< 3; i++){
		fHistCorrectionFactorsHistFitted->Add(fHistCorrectionFactorsHistCatFitted[i]);
		for (Int_t k = 0; k < 3; k++){
			if (i== 0){
				fHistCorrectionFactorsHist[k] 	= new TH1D(Form("fHistCorrectionFactorsHistvsPt_%i",k),"",fNBinsPt,fBinsPt);
			}
			fHistCorrectionFactorsHist[k]->Add(fHistCorrectionFactorsHistCat[k][i]);
		}   
	}
		
	Color_t colorMethod[3] 						= {kRed+2, kGreen+2, kBlue+2};
	Style_t styleMethod[3] 						= {21,29,33};
	TString nameMethod[3] 						= {"A","C","D"};
		
	canvasCorrFrac->cd();
		DrawAutoGammaMesonHistos( fHistCorrectionFactorsHistAllCat, 
								  "", Form("p_{T,%s} (GeV/c)",fMesonType.Data()), "Correction factor (%)", 
								  kFALSE, 2.,1e-8, kFALSE,
								  kTRUE, -1, 30, 
								  kFALSE, 0., 10.);
		DrawGammaSetMarker(fHistCorrectionFactorsHistAllCat, 20, 0.8, kBlack, kBlack);      
		fHistCorrectionFactorsHistAllCat->GetYaxis()->SetTitleOffset(0.8);
		fHistCorrectionFactorsHistAllCat->DrawCopy("e1,p"); 

		if (!fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 && optionPeriod.CompareTo("") == 0 && !meson.Contains("Eta")){
			DrawGammaSetMarker(fHistCorrectionFactorsFitAllCat, 25, 0.8, kCyan+2, kCyan+2);
			fHistCorrectionFactorsFitAllCat->DrawCopy("same,e1,p"); 
		}
		
		TLegend* legendCorrFraction 			= new TLegend(0.5,0.8,0.85,0.95);
		legendCorrFraction->SetTextSize(0.04);         
		legendCorrFraction->SetFillColor(0);
		legendCorrFraction->SetLineColor(0);
		
		legendCorrFraction->AddEntry(fHistCorrectionFactorsHistAllCat,"Extraction Method A","p");
		if (!fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 && optionPeriod.CompareTo("") == 0 && !meson.Contains("Eta")) 
			legendCorrFraction->AddEntry(fHistCorrectionFactorsFitAllCat,"Extraction Method B","p");
		
		for (Int_t k = 0; k < 3 ; k++){
			DrawGammaSetMarker(fHistCorrectionFactorsHist[k], styleMethod[k], 0.8, colorMethod[k], colorMethod[k]);
			fHistCorrectionFactorsHist[k]->DrawCopy("same,e1,p"); 
			legendCorrFraction->AddEntry(fHistCorrectionFactorsHist[k],Form("Extraction Method %s (sep Cat)",nameMethod[k].Data()),"p");      
		}        
		legendCorrFraction->Draw();
		
		labelEnergy->Draw();

		
	canvasCorrFrac->Update(); 
	canvasCorrFrac->SaveAs(Form("%s/%s_%s_CorrectionFactorTotal.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), suffix.Data()));
		
	canvasCorrFrac->cd();
	
		fHistCorrectionFactorsHist[0]->DrawCopy("e1,p"); 
	
		TLegend* legendCorrFractionCatHist = new TLegend(0.65,0.7,0.85,0.95);
		legendCorrFractionCatHist->SetTextSize(0.04);         
		legendCorrFractionCatHist->SetFillColor(0);
		legendCorrFractionCatHist->SetLineColor(0);
		legendCorrFractionCatHist->AddEntry(fHistCorrectionFactorsHist[0],"Total","p");
		
		
		for (Int_t i = 0; i< 3; i++){
			DrawGammaSetMarker(fHistCorrectionFactorsHistCat[0][i], styleCat[i], 0.8, colorCat[i], colorCat[i]);
			fHistCorrectionFactorsHistCat[0][i]->DrawCopy("same,e1,p"); 
			legendCorrFractionCatHist->AddEntry(fHistCorrectionFactorsHistCat[0][i],Form("Category %i",i+1),"p");
		}   
		legendCorrFractionCatHist->Draw();
		labelEnergy->Draw();
	canvasCorrFrac->Update(); 
	canvasCorrFrac->SaveAs(Form("%s/%s_%s_CorrectionFactorVsCatHist.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), suffix.Data()));
	
	
	
	canvasCorrFrac->cd();
		DrawAutoGammaMesonHistos( fHistFracCatvsPt[0], 
								"", "p_{T,meson} (GeV/c)", "N_{meson per cat}/(N_{meson}) (%)", 
								kFALSE, 2.,1e-8, kFALSE,
								kTRUE, 0, 100, 
								kFALSE, 0., 10.);
		DrawGammaSetMarker(fHistFracCatvsPt[0], styleCat[0], 1., colorCat[0], colorCat[0]);      
		fHistFracCatvsPt[0]->GetYaxis()->SetTitleOffset(0.8);
		fHistFracCatvsPt[0]->DrawCopy("p,e1"); 
	
		TLegend* legendFractionCat = new TLegend(0.65,0.7,0.85,0.95);
		legendFractionCat->SetTextSize(0.04);         
		legendFractionCat->SetFillColor(0);
		legendFractionCat->SetLineColor(0);
		legendFractionCat->AddEntry(fHistFracCatvsPt[0],"Category 1","p");
		
		
		for (Int_t i = 1; i< 6; i++){
			DrawGammaSetMarker(fHistFracCatvsPt[i], styleCat[i], 1., colorCat[i], colorCat[i]);
			fHistFracCatvsPt[i]->DrawCopy("same,p,e1"); 
			legendFractionCat->AddEntry(fHistFracCatvsPt[i],Form("Category %i",i+1),"p");
		}   
		legendFractionCat->Draw();
		labelEnergy->Draw();
	canvasCorrFrac->Update(); 
	canvasCorrFrac->SaveAs(Form("%s/%s_%s_FractionPerCatVsPt.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), suffix.Data()));
	
	for (Int_t j = fStartPtBin; j < fNBinsPt; j++){
		TCanvas* canvasDCABGCheck = new TCanvas("canvasDCABGCheck","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasDCABGCheck, 0.08, 0.02, 0.02, 0.09);
		canvasDCABGCheck->SetLogy();
		
		TH1D* fHistDCAZUnderMeson_Visual_AllCat =  (TH1D*)fHistDCAZUnderMeson_MesonPt_AllCat[j]->Clone("fHistDCAZUnderMeson_Visual_AllCat");
		fHistDCAZUnderMeson_Visual_AllCat->Sumw2();
		fHistDCAZUnderMeson_Visual_AllCat->Scale(1./fHistDCAZUnderMeson_Visual_AllCat->GetEntries());
		TH1D* fHistDCAZUnderMesonBG1_Visual_AllCat =  (TH1D*)fHistDCAZUnderMesonBG1_MesonPt_AllCat[j]->Clone("fHistDCAZUnderMesonBG1_Visual_AllCat");
		fHistDCAZUnderMesonBG1_Visual_AllCat->Sumw2();
		fHistDCAZUnderMesonBG1_Visual_AllCat->Scale(1./fHistDCAZUnderMesonBG1_Visual_AllCat->GetEntries());
		TH1D* fHistDCAZUnderMesonBG2_Visual_AllCat =  (TH1D*)fHistDCAZUnderMesonBG2_MesonPt_AllCat[j]->Clone("fHistDCAZUnderMesonBG2_Visual_AllCat");
		fHistDCAZUnderMesonBG2_Visual_AllCat->Sumw2();
		fHistDCAZUnderMesonBG2_Visual_AllCat->Scale(1./fHistDCAZUnderMesonBG2_Visual_AllCat->GetEntries());
			if (fHistDCAZUnderMeson_Visual_AllCat){
				DrawAutoGammaMesonHistos( fHistDCAZUnderMeson_Visual_AllCat, 
								"","dca_{z} #gamma (cm)", "d(dca_{z})", 
								kTRUE, 2.,5e-5, kFALSE,
								kFALSE, 0, 50, 
								kTRUE, -6., 6.);
				DrawGammaSetMarker(fHistDCAZUnderMeson_Visual_AllCat, 20, 1., kBlack, kBlack);
				fHistDCAZUnderMeson_Visual_AllCat->GetYaxis()->SetTitleOffset(0.8);
				fHistDCAZUnderMeson_Visual_AllCat->DrawCopy("p,e1"); 

				}
				if (fHistDCAZUnderMesonBG1_Visual_AllCat){
				DrawGammaSetMarker(fHistDCAZUnderMesonBG1_Visual_AllCat, 24, 1., kRed+2, kRed+2);
				fHistDCAZUnderMesonBG1_Visual_AllCat->DrawCopy("same,p,e1"); 

				} 
				if (fHistDCAZUnderMesonBG2_Visual_AllCat){
				DrawGammaSetMarker(fHistDCAZUnderMesonBG2_Visual_AllCat, 24, 1., kGreen+2, kGreen+2);
				fHistDCAZUnderMesonBG2_Visual_AllCat->DrawCopy("same,p,e1"); 
				}
				
			TLatex *labelPtBin = new TLatex(0.63,0.9,Form("%1.2f < p_{T} (GeV/c) < %1.2f", fBinsPt[j], fBinsPt[j+1]));
			SetStyleTLatex( labelPtBin, 0.04,4);   
			TLegend* legendDCADataGGBGComp = new TLegend(0.57,0.73,0.8,0.88);
			legendDCADataGGBGComp->SetTextSize(0.035);         
			legendDCADataGGBGComp->SetFillColor(0);
			legendDCADataGGBGComp->SetLineColor(0);
			legendDCADataGGBGComp->SetNColumns(1);
	//          legendDCADataGGBGComp->AddEntry((TObject*)0, ,"");
			if (fHistDCAZUnderMeson_Visual_AllCat)
				legendDCADataGGBGComp->AddEntry(fHistDCAZUnderMeson_Visual_AllCat,Form("S + BG, %1.2f < M_{#gamma#gamma} (GeV/c^{2}) < %1.2f ",fMesonIntRange[0], fMesonIntRange[1]),"p");
			if (fHistDCAZUnderMesonBG1_Visual_AllCat)
				legendDCADataGGBGComp->AddEntry(fHistDCAZUnderMesonBG1_Visual_AllCat,Form("#gamma#gamma, BG %1.2f < M_{#gamma#gamma} (GeV/c^{2}) < %1.2f ", fMesonIntRangeBG1[0], 
																						  fMesonIntRangeBG1[1]), "p");
			if (fHistDCAZUnderMesonBG2_Visual_AllCat)
				legendDCADataGGBGComp->AddEntry(fHistDCAZUnderMesonBG2_Visual_AllCat,Form("#gamma#gamma, BG %1.2f < M_{#gamma#gamma} (GeV/c^{2}) < %1.2f", fMesonIntRangeBG2[0], 
																						  fMesonIntRangeBG2[1]), "p");
			
			legendDCADataGGBGComp->Draw();         
			labelEnergy->Draw();
			labelPtBin->Draw();
		canvasDCABGCheck->Update(); 
		canvasDCABGCheck->SaveAs(Form("%s/%s_%s_DCAzGGBGComp_%i.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), j,suffix.Data()));
		delete canvasDCABGCheck;
	}
	
	
	if (kMC){
		TCanvas* canvasDCAMCComponents = new TCanvas("canvasDCAMCComponents","",200,10,1350,900);  // gives the page size
		DrawGammaCanvasSettings( canvasDCAMCComponents, 0.08, 0.02, 0.02, 0.09);
		canvasDCAMCComponents->SetLogy();
			if (fHistDCAZUnderMesonAllCat_AllPt){
				DrawAutoGammaMesonHistos( fHistDCAZUnderMesonAllCat_AllPt, 
								"","dca_{z} #gamma (cm)", "d(dca_{z})", 
								kTRUE, 2.,1, kFALSE,
								kFALSE, 0, 50, 
								kTRUE, -6., 6.);
				DrawGammaSetMarker(fHistDCAZUnderMesonAllCat_AllPt, 20, 1., kBlack, kBlack);
				fHistDCAZUnderMesonAllCat_AllPt->GetYaxis()->SetTitleOffset(0.8);
				fHistDCAZUnderMesonAllCat_AllPt->DrawCopy("p,e1"); 
			}
			if (fHistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt){
				DrawGammaSetMarker(fHistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt, 20, 1., kRed+2, kRed+2);
				fHistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt->DrawCopy("same,p,e1"); 

			} 
			if (fHistDCAZTruePrimaryMesonDalitzAllCat_AllPt){
				DrawGammaSetMarker(fHistDCAZTruePrimaryMesonDalitzAllCat_AllPt, 20, 1., kGreen+2, kGreen+2);
				fHistDCAZTruePrimaryMesonDalitzAllCat_AllPt->DrawCopy("same,p,e1"); 
			}
			if (fHistDCAZTrueBackgroundAllCat_AllPt){
				DrawGammaSetMarker(fHistDCAZTrueBackgroundAllCat_AllPt, 20, 1., kPink+2, kPink+2);
				fHistDCAZTrueBackgroundAllCat_AllPt->DrawCopy("same,p,e1"); 
			}
			if (fHistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt){
				DrawGammaSetMarker(fHistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt, 20, 1., 807, 807);
				fHistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt->DrawCopy("same,p,e1"); 
			}
			if (fHistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt){
				DrawGammaSetMarker(fHistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt, 20, 1., kViolet+2, kViolet+2);
				fHistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt->DrawCopy("same,p,e1"); 
			}
			if (fHistDCAZGarbageAllCat_AllPt){
				DrawGammaSetMarker(fHistDCAZGarbageAllCat_AllPt, 20, 1., kCyan+2, kCyan+2);
				fHistDCAZGarbageAllCat_AllPt->DrawCopy("same,p,e1"); 
			}
				
				
			TLegend* legendDCAMCComponents0 = new TLegend(0.7,0.7,0.85,0.95);
			legendDCAMCComponents0->SetTextSize(0.04);         
			legendDCAMCComponents0->SetFillColor(0);
			legendDCAMCComponents0->SetLineColor(0);
			legendDCAMCComponents0->SetNColumns(1);
			if (fHistDCAZUnderMesonAllCat_AllPt)
				legendDCAMCComponents0->AddEntry(fHistDCAZUnderMesonAllCat_AllPt,"Total","p");
			if (fHistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt)
				legendDCAMCComponents0->AddEntry(fHistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt,Form("Prim %s",fMesonType.Data()),"p");
			if (fHistDCAZTruePrimaryMesonDalitzAllCat_AllPt)
				legendDCAMCComponents0->AddEntry(fHistDCAZTruePrimaryMesonDalitzAllCat_AllPt,Form("Dalitz %s",fMesonType.Data()),"p");
			if (fHistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt && fHistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt->GetEntries() > 0)
				legendDCAMCComponents0->AddEntry(fHistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt,"Sec. #pi^{0} from K^{0}_{s}","p");
			if (fHistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt && fHistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt->GetEntries() > 0)
				legendDCAMCComponents0->AddEntry(fHistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt,"Sec. #pi^{0} from X","p");
			if (fHistDCAZTrueBackgroundAllCat_AllPt)
				legendDCAMCComponents0->AddEntry(fHistDCAZTrueBackgroundAllCat_AllPt,"#gamma#gamma BG","p");
			if (fHistDCAZGarbageAllCat_AllPt)
				legendDCAMCComponents0->AddEntry(fHistDCAZGarbageAllCat_AllPt,"other","p"); //"garbage","p");
			legendDCAMCComponents0->Draw();         
			labelEnergy->Draw();
		canvasDCAMCComponents->Update(); 
		canvasDCAMCComponents->SaveAs(Form("%s/%s_%s_DCAzMCComponentsAll.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), suffix.Data()));
		
		
		canvasDCAMCComponents->cd();
			if (fHistDCAZUnderMesonAllCat_AllPt){
				fHistDCAZUnderMesonAllCat_AllPt->Scale(1./fHistDCAZUnderMesonAllCat_AllPt->GetEntries());
				DrawAutoGammaMesonHistos( fHistDCAZUnderMesonAllCat_AllPt, 
								"","dca_{z} #gamma (cm)", "d(dca_{z})", 
								kTRUE, 2.,0.5e-4, kFALSE,
								kFALSE, 0, 50, 
								kTRUE, -6., 6.);
				DrawGammaSetMarker(fHistDCAZUnderMesonAllCat_AllPt, 20, 1., kBlack, kBlack);
				fHistDCAZUnderMesonAllCat_AllPt->GetYaxis()->SetTitleOffset(0.8);
				fHistDCAZUnderMesonAllCat_AllPt->DrawCopy("p,e1"); 

			}
			if (fHistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt){
				fHistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt->Scale(1./fHistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt->GetEntries());
				fHistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt->DrawCopy("same,p,e1"); 
			} 
			if (fHistDCAZTruePrimaryMesonDalitzAllCat_AllPt){
				fHistDCAZTruePrimaryMesonDalitzAllCat_AllPt->Scale(1./fHistDCAZTruePrimaryMesonDalitzAllCat_AllPt->GetEntries());
				fHistDCAZTruePrimaryMesonDalitzAllCat_AllPt->DrawCopy("same,p,e1"); 
			}
			if (fHistDCAZTrueBackgroundAllCat_AllPt){
				fHistDCAZTrueBackgroundAllCat_AllPt->Scale(1./fHistDCAZTrueBackgroundAllCat_AllPt->GetEntries());
				fHistDCAZTrueBackgroundAllCat_AllPt->DrawCopy("same,p,e1"); 
			}
			
				
			TLegend* legendDCAMCComponents = new TLegend(0.7,0.7,0.85,0.95);
			legendDCAMCComponents->SetTextSize(0.04);         
			legendDCAMCComponents->SetFillColor(0);
			legendDCAMCComponents->SetLineColor(0);
			legendDCAMCComponents->SetNColumns(1);
			if (fHistDCAZUnderMesonAllCat_AllPt)
				legendDCAMCComponents->AddEntry(fHistDCAZUnderMesonAllCat_AllPt,"Total","p");
			if (fHistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt)
				legendDCAMCComponents->AddEntry(fHistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt,Form("Prim. %s",fMesonType.Data()),"p");
			if (fHistDCAZTruePrimaryMesonDalitzAllCat_AllPt)
				legendDCAMCComponents->AddEntry(fHistDCAZTruePrimaryMesonDalitzAllCat_AllPt,Form("Dalitz %s",fMesonType.Data()),"p");
			if (fHistDCAZTrueBackgroundAllCat_AllPt)
				legendDCAMCComponents->AddEntry(fHistDCAZTrueBackgroundAllCat_AllPt,"#gamma#gamma BG","p");
			legendDCAMCComponents->Draw();         
			labelEnergy->Draw();
		canvasDCAMCComponents->Update(); 
		canvasDCAMCComponents->SaveAs(Form("%s/%s_%s_DCAzMCComponentsSignalAndBackground.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), suffix.Data()));
		
		canvasDCAMCComponents->cd();
			if (fHistDCAZUnderMesonAllCat_AllPt){
				fHistDCAZUnderMesonAllCat_AllPt->DrawCopy("p,e1"); 
			}
			if (fHistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt){
				fHistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt->DrawCopy("same,p,e1"); 
			} 
			if (fHistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt){
				fHistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt->Scale(1./fHistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt->GetEntries());
				fHistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt->DrawCopy("same,p,e1"); 
			}
			if (fHistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt){
				fHistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt->Scale(1./fHistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt->GetEntries());
				fHistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt->DrawCopy("same,p,e1"); 
			}
			if (fHistDCAZGarbageAllCat_AllPt){
				fHistDCAZGarbageAllCat_AllPt->Scale(1./fHistDCAZGarbageAllCat_AllPt->GetEntries());
				fHistDCAZGarbageAllCat_AllPt->DrawCopy("same,p,e1"); 
			}
			
			TLegend* legendDCAMCComponents2 = new TLegend(0.7,0.7,0.85,0.95);
			legendDCAMCComponents2->SetTextSize(0.04);         
			legendDCAMCComponents2->SetFillColor(0);
			legendDCAMCComponents2->SetLineColor(0);
			legendDCAMCComponents2->SetNColumns(1);
			if (fHistDCAZUnderMesonAllCat_AllPt)
				legendDCAMCComponents2->AddEntry(fHistDCAZUnderMesonAllCat_AllPt,"Total","p");
			if (fHistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt)
				legendDCAMCComponents2->AddEntry(fHistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt,Form("Prim %s",fMesonType.Data()),"p");
			if (fHistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt && fHistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt->GetEntries() > 0 )
				legendDCAMCComponents2->AddEntry(fHistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt,"Sec. #pi^{0} from K^{0}_{s}","p");
			if (fHistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt && fHistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt->GetEntries() > 0)
				legendDCAMCComponents2->AddEntry(fHistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt,"Sec. #pi^{0} from X","p");
			if (fHistDCAZGarbageAllCat_AllPt)
				legendDCAMCComponents2->AddEntry(fHistDCAZGarbageAllCat_AllPt,"garbage","p");
			legendDCAMCComponents2->Draw();         
			labelEnergy->Draw();
		canvasDCAMCComponents->Update(); 
		canvasDCAMCComponents->SaveAs(Form("%s/%s_%s_DCAzMCComponentsSecondaries.%s", outputDir.Data(), meson.Data(), fMCFlag.Data(), suffix.Data()));
		
		
	}
	
	
	const char* nameOutput = Form("%s/%s/%s_%s_GammaConvV1DCATestAnalysed%s.root",fCutSelection.Data(),optionEnergy.Data(),meson.Data(),fMCFlag.Data(),optionPeriod.Data());
	TFile*   fOutput1 = new TFile(nameOutput,"RECREATE");
	
	for (Int_t i = 0; i < 6 ; i++){
		for (Int_t j = fStartPtBin; j < fNBinsPt ; j++){
			if (i == 0){
				if (fHistDCAZUnderMeson_MesonPt_AllCat[j]) fHistDCAZUnderMeson_MesonPt_AllCat[j]->Write();
				if (fHistDCAZUnderMesonBG1_MesonPt_AllCat[j]) fHistDCAZUnderMesonBG1_MesonPt_AllCat[j]->Write();
				if (fHistDCAZUnderMesonBG2_MesonPt_AllCat[j]) fHistDCAZUnderMesonBG2_MesonPt_AllCat[j]->Write();
				if (fHistDCAZUnderMeson_MesonPt_GoodCat[j]) fHistDCAZUnderMeson_MesonPt_GoodCat[j]->Write();
				if (kMC){
				if (fHistDCAZTruePrimaryMesonGammaGamma_MesonPt_AllCat[j]) fHistDCAZTruePrimaryMesonGammaGamma_MesonPt_AllCat[j]->Write();
				if (fHistDCAZTruePrimaryMesonDalitz_MesonPt_AllCat[j]) fHistDCAZTruePrimaryMesonDalitz_MesonPt_AllCat[j]->Write();
				if (fHistDCAZTrueSecondaryMesonFromK0s_MesonPt_AllCat[j]) fHistDCAZTrueSecondaryMesonFromK0s_MesonPt_AllCat[j]->Write();
				if (fHistDCAZTrueSecondaryMesonFromEta_MesonPt_AllCat[j]) fHistDCAZTrueSecondaryMesonFromEta_MesonPt_AllCat[j]->Write();
				if (fHistDCAZTrueSecondaryMesonFromSomething_MesonPt_AllCat[j]) fHistDCAZTrueSecondaryMesonFromSomething_MesonPt_AllCat[j]->Write();
				if (fHistDCAZTrueBackground_MesonPt_AllCat[j]) fHistDCAZTrueBackground_MesonPt_AllCat[j]->Write();
				if (fHistDCAZGarbage_MesonPt_AllCat[j]) fHistDCAZGarbage_MesonPt_AllCat[j]->Write();
				}  
			}   
			
			if (fHistDCAZUnderMeson_MesonPt[i][j]) fHistDCAZUnderMeson_MesonPt[i][j]->Write();
			if (fHistDCAZUnderMesonBG1_MesonPt[i][j]) fHistDCAZUnderMesonBG1_MesonPt[i][j]->Write();
			if (fHistDCAZUnderMesonBG2_MesonPt[i][j]) fHistDCAZUnderMesonBG2_MesonPt[i][j]->Write();
			if (kMC){
				if (fHistDCAZTruePrimaryMesonGammaGamma_MesonPt[i][j]) fHistDCAZTruePrimaryMesonGammaGamma_MesonPt[i][j]->Write();
				if (fHistDCAZTruePrimaryMesonDalitz_MesonPt[i][j]) fHistDCAZTruePrimaryMesonDalitz_MesonPt[i][j]->Write();
				if (fHistDCAZTrueSecondaryMesonFromK0s_MesonPt[i][j]) fHistDCAZTrueSecondaryMesonFromK0s_MesonPt[i][j]->Write();
				if (fHistDCAZTrueSecondaryMesonFromEta_MesonPt[i][j]) fHistDCAZTrueSecondaryMesonFromEta_MesonPt[i][j]->Write();
				if (fHistDCAZTrueSecondaryMesonFromSomething_MesonPt[i][j]) fHistDCAZTrueSecondaryMesonFromSomething_MesonPt[i][j]->Write();
				if (fHistDCAZTrueBackground_MesonPt[i][j]) fHistDCAZTrueBackground_MesonPt[i][j]->Write();
				if (fHistDCAZGarbage_MesonPt[i][j]) fHistDCAZGarbage_MesonPt[i][j]->Write();
			}  
			if (fHistDCAZUnderMesonBGEstimate_MesonPt[0][i][j]) fHistDCAZUnderMesonBGEstimate_MesonPt[0][i][j]->Write();
		}
	}   

	for (Int_t i = 0; i < 6 ; i++){
		if (fHistDCAZUnderMeson_AllPt[i]) fHistDCAZUnderMeson_AllPt[i]->Write();
		if (fHistFracCatvsPt[i]) fHistFracCatvsPt[i]->Write();
		for (Int_t k = 0; k < 3; k++){
			if (fHistFracIntHistBGvsPt[k][i]) fHistFracIntHistBGvsPt[k][i]->Write();
			if (fHistCorrectionFactorsHistCat[k][i]) fHistCorrectionFactorsHistCat[k][i]->Write();
		}   
		if (fHistFracIntHistBGFittedvsPt[0][i]) fHistFracIntHistBGFittedvsPt[0][i]->Write();
		if (fHistCorrectionFactorsHistCatFitted[i]) fHistCorrectionFactorsHistCatFitted[i]->Write();
		if (kMC){
			if (fHistDCAZGarbage_AllPt[i]) fHistDCAZGarbage_AllPt[i]->Write();
			if (fHistDCAZTrueBackground_AllPt[i]) fHistDCAZTrueBackground_AllPt[i]->Write();
			if (fHistDCAZTrueSecondaryMesonFromSomething_AllPt[i]) fHistDCAZTrueSecondaryMesonFromSomething_AllPt[i]->Write();
			if (fHistDCAZTrueSecondaryMesonFromEta_AllPt[i]) fHistDCAZTrueSecondaryMesonFromEta_AllPt[i]->Write();
			if (fHistDCAZTrueSecondaryMesonFromK0s_AllPt[i]) fHistDCAZTrueSecondaryMesonFromK0s_AllPt[i]->Write();
			if (fHistDCAZTruePrimaryMesonDalitz_AllPt[i]) fHistDCAZTruePrimaryMesonDalitz_AllPt[i]->Write();
			if (fHistDCAZTruePrimaryMesonGammaGamma_AllPt[i]) fHistDCAZTruePrimaryMesonGammaGamma_AllPt[i]->Write();
		}
	}   
	if (fHistDCAZUnderMesonAllCat_AllPt){
		if (kMC) fHistDCAZUnderMesonAllCat_AllPt->Scale(fHistDCAZUnderMesonAllCat_AllPt->GetEntries()/fNEvents);
		else fHistDCAZUnderMesonAllCat_AllPt->Scale(1./fNEvents);
		fHistDCAZUnderMesonAllCat_AllPt->Write();
	}
	if (kMC){
		if (fHistDCAZGarbageAllCat_AllPt){
			fHistDCAZGarbageAllCat_AllPt->Scale(fHistDCAZGarbageAllCat_AllPt->GetEntries()/fNEvents);
			fHistDCAZGarbageAllCat_AllPt->Write();
		}
		if (fHistDCAZTrueBackgroundAllCat_AllPt){
			fHistDCAZTrueBackgroundAllCat_AllPt->Scale(fHistDCAZTrueBackgroundAllCat_AllPt->GetEntries()/fNEvents);
			fHistDCAZTrueBackgroundAllCat_AllPt->Write();
		}
		if (fHistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt){
			fHistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt->Scale(fHistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt->GetEntries()/fNEvents);
			fHistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt->Write();
		}   
		if (fHistDCAZTrueSecondaryMesonFromEtaAllCat_AllPt){
			fHistDCAZTrueSecondaryMesonFromEtaAllCat_AllPt->Scale(1./fNEvents);
			fHistDCAZTrueSecondaryMesonFromEtaAllCat_AllPt->Write();
		}
		if (fHistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt){
			fHistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt->Scale(fHistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt->GetEntries()/fNEvents);
			fHistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt->Write();  
		}
		if (fHistDCAZTruePrimaryMesonDalitzAllCat_AllPt){
			fHistDCAZTruePrimaryMesonDalitzAllCat_AllPt->Scale(fHistDCAZTruePrimaryMesonDalitzAllCat_AllPt->GetEntries()/fNEvents);
			fHistDCAZTruePrimaryMesonDalitzAllCat_AllPt->Write();
		}   
		if (fHistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt){
			fHistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt->Scale(fHistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt->GetEntries()/fNEvents);
			fHistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt->Write();
			
		}    
	}   
	for (Int_t k = 0; k < 3; k++){
		if (fHistCorrectionFactorsHist[k]) fHistCorrectionFactorsHist[k]->Write();
	} 
	if (fHistCorrectionFactorsHistFitted) fHistCorrectionFactorsHistFitted->Write();
	if (!fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 && optionPeriod.CompareTo("") == 0 && !meson.Contains("Eta")) 
		if (fHistCorrectionFactorsFitAllCat) fHistCorrectionFactorsFitAllCat->Write();
	if (fHistCorrectionFactorsHistAllCat) fHistCorrectionFactorsHistAllCat->Write();
	if (fEventQuality) fEventQuality->Write();
	fOutput1->Write();
	fOutput1->Close();

  
}

Double_t GausWithHole(Double_t *x, Double_t *par){
	Double_t f=0;

	if (x[0] > -fHoleRadius && x[0] < fHoleRadius ){
		TF1::RejectPoint();
		return 0;
	}
	
	f = par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2], 2.));
	
	return f; 
   
}



Double_t AsymmGaus(Double_t *x, Double_t *par){
	Double_t f=0;
	
	if(x[0]<=par[1])
		f = par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2], 2.));
	else
		f = par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[3], 2.));

	return f; 
   
}


void Initialize(TString setPi0, TString cent, TString optPeriod, Int_t numberOfBins){
	if (setPi0.CompareTo("Pi0") == 0){
		fNBinsPt 			= numberOfBins;
		fBinsPt 						= new Double_t[33];
                 if (fEnergyFlag.CompareTo("13TeV") == 0) {
			fStartPtBin					= 1;
			fColumn 					= 5;
			fRow 						= 4;
			if (fNBinsPt < 17) 
			        fColumn                                = 4;
			if (fNBinsPt > 20) {
				cout << "You have chosen to have more than 20 bins, this is not possible, it will be reduced to 20" << endl;
				fNBinsPt 				= 20;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] 				= fBinsPi013TeVPtDCA[i];
			}
			fExampleBin 			         	= 7;
			nIterBGFit 					= 12;
			fMaxYFracBGOverIntHist                          = 30;
			optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing3";
			optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing5";
			optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing7";
			Int_t nBinsPlot 			= fColumn*fRow -1;
			if (fNBinsPt-fStartPtBin > nBinsPlot) fColumn++;
			nBinsPlot 					= fColumn*fRow -1;
			if (fNBinsPt-fStartPtBin > nBinsPlot) fRow++;			
		}
                else if (fEnergyFlag.CompareTo("7TeV") == 0) {
			fStartPtBin					= 1;
			fColumn 					= 5;
			fRow 						= 4;

			if (fNBinsPt > 21) {
				cout << "You have chosen to have more than 21 bins, this is not possible, it will be reduced to 21" << endl;
				fNBinsPt 				= 21;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] 				= fBinsPi07TeVPtDCA[i];
			}
			fExampleBin 				= 10;
			nIterBGFit 					= 13;
		} else if (fEnergyFlag.CompareTo("8TeV") == 0) {
			fStartPtBin					= 1;
			fColumn 					= 5;
			fRow 						= 4;

			if (fNBinsPt > 21) {
				cout << "You have chosen to have more than 21 bins, this is not possible, it will be reduced to 21" << endl;
				fNBinsPt 				= 21;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] 				= fBinsPi08TeVPtDCA[i];
			}
			fExampleBin 				= 10;
			nIterBGFit 					= 12;
			optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing5";
			optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing3";
			optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing7";
			Int_t nBinsPlot 			= fColumn*fRow -1;
			if (fNBinsPt-fStartPtBin > nBinsPlot) fColumn++;
			nBinsPlot 					= fColumn*fRow -1;
			if (fNBinsPt-fStartPtBin > nBinsPlot) fRow++;			
		}	else if (fEnergyFlag.CompareTo("2.76TeV") == 0) {
			fStartPtBin 				= 1;
			fColumn 					= 5;
			fRow 						= 3;

			if (fNBinsPt > 14) {
				cout << "You have chosen to have more than 14 bins, this is not possible, it will be reduced to 19" << endl;
				fNBinsPt 				= 14;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] 				= fBinsPi02760GeVPtDCA[i];
			}
			fExampleBin 				= 7;
			nIterBGFit 					= 13;
		} else if (fEnergyFlag.CompareTo("900GeV") == 0) {
			fStartPtBin					= 1;
			fColumn 					= 4;
			fRow 						= 3;

			if (fNBinsPt > 11) {
				cout << "You have chosen to have more than 11 bins, this is not possible, it will be reduced to 11" << endl;
				fNBinsPt 				= 11;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] 				= fBinsPi0900GeVPt[i];
			}
			fExampleBin 				= 4;
			nIterBGFit 					= 13;
		} else if( fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0) { 
			fStartPtBin					= 1;
			fColumn 					= 5;
			fRow 						= 3;
			
			if (fNBinsPt > 15) {
				cout << "You have chosen to have more than 15 bins, this is not possible, it will be reduced to 15" << endl;
				fNBinsPt 				= 15;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				if (cent.CompareTo("60-80%")==0 || cent.CompareTo("70-80%")==0 || cent.CompareTo("75-90%")==0 ){
					fBinsPt[i] 			= fBinsPi0HIPtDCAPer[i];
				} else {
					fBinsPt[i] 			= fBinsPi0HIPtDCA[i];
				}   	
			}
			fExampleBin 	= 4;
			cout << cent.Data() << endl;
			if (cent.CompareTo("60-80%")==0){
				fColumn 					= 4;
				optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing5";
				optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing3";
				optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing7";
				nIterBGFit 					= 15;
				fMaxYFracBGOverIntHist		= 15;
			} else if (cent.CompareTo("75-90%")==0){
				fColumn 					= 4;
				optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing9";
				optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing7";
				optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing11";
				nIterBGFit 					= 14;
				fMaxYFracBGOverIntHist		= 50;
			} else if (cent.CompareTo("70-80%")==0){
				fColumn 					=   4;
				optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing9";
				optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing7";
				optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing11";
				nIterBGFit 					= 14;
				fMaxYFracBGOverIntHist		= 50;   
			} else if (cent.CompareTo("60-70%")==0){
				fColumn 					= 5;
				if (optPeriod.CompareTo("LHC11h")==0){
					optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing7";
					optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing5";
					optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing9";
					nIterBGFit 					= 13;
					fMaxYFracBGOverIntHist		= 15;
				} else {
					optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing7";
					optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing5";
					optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing9";
					nIterBGFit 					= 17;
					fMaxYFracBGOverIntHist		= 15;	      
				}  
			} else if (cent.CompareTo("50-60%")==0){
				optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing5";
				optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing3";
				optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing7";
				nIterBGFit 					= 14;
				fMaxYFracBGOverIntHist		= 12;
			} else if (cent.CompareTo("30-40%")==0) {
				optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing5";
				optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing3";
				optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing7";
				nIterBGFit 					= 18;
				fMaxYFracBGOverIntHist		= 12; 
			} else if (cent.CompareTo("20-50%")==0) {
				optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing5";
				optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing3";
				optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing7";
				nIterBGFit 					= 18;
				fMaxYFracBGOverIntHist		= 12;
			} else if (cent.CompareTo("20-30%")==0) {
				optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing5";
				optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing3";
				optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing7";
				nIterBGFit 					= 18;
				fMaxYFracBGOverIntHist		= 12;
			} else if (cent.CompareTo("30-50%")==0) {
				optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing5";
				optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing3";
				optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing7";
				nIterBGFit 					= 18;
				fMaxYFracBGOverIntHist		= 12;  
			} else if (cent.CompareTo("40-50%")==0) {
				optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing5";
				optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing3";
				optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing7";
				nIterBGFit 					= 16;
				fMaxYFracBGOverIntHist		= 12;   
			} else if (cent.CompareTo("40-60%") == 0){
				optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing5";
				optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing3";
				optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing7";
				nIterBGFit 					= 17;
				fMaxYFracBGOverIntHist		= 10;
			} else if (cent.CompareTo("20-40%") == 0){
				optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing5";
				optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing3";
				optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing7";
				nIterBGFit 					= 19;      
				fMaxYFracBGOverIntHist		= 8;
			} else {
				fMaxYFracBGOverIntHist		= 4;
				nIterBGFit 					= 21;
			}   
		} else if( fEnergyFlag.CompareTo("pPb_5.023TeV") == 0) { 
			fStartPtBin					= 1;
			fColumn 					= 3;
			fRow 						= 4;
			
			if (fNBinsPt > 12) {
				cout << "You have chosen to have more than 12 bins, this is not possible, it will be reduced to 12" << endl;
				fNBinsPt 				= 12;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] 				= fBinsPi0pPbPtDCA[i];
			}
			
			fExampleBin 			 	= 7;
			optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing5";
			optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing3";
			optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing7";
			nIterBGFit 					= 13;
			fMaxYFracBGOverIntHist		= 20;
		}
		fMesonIntRange 					=  new Double_t[2];	fMesonIntRange[0] 			= 0.08;		fMesonIntRange[1]			= 0.2;
		fMesonIntRangeBG1 				=  new Double_t[2];	fMesonIntRangeBG1[0] 		= 0.08;		fMesonIntRangeBG1[1]		= 0.1;
		fMesonIntRangeBG2 				=  new Double_t[2];	fMesonIntRangeBG2[0] 		= 0.15;		fMesonIntRangeBG2[1]		= 0.2;
		numberInvMassBins 				= 50;
		fMesonId						= 111;
		fMesonWidthExpect 				= 0.003;
		fMesonLambdaTail 				= 0.012;
		fMesonWidthRange 				= new Double_t[2];	fMesonWidthRange[0]			= 0.001;	fMesonWidthRange[1] 		= 0.009;
		fMesonLambdaTailRange 			= new Double_t[2];	fMesonLambdaTailRange[0]	= 0.001;	fMesonLambdaTailRange[1] 	= 0.02;

	} else if (setPi0.CompareTo("Eta") == 0){
		fNBinsPt 						= numberOfBins;
		fBinsPt							= new Double_t[20];
		if (fEnergyFlag.CompareTo("13TeV") == 0) {
			fStartPtBin					= 1;
			fColumn 					= 3;
			fRow                = 3;
			if ((fNBinsPt-fStartPtBin) < 6) 
			        fRow            = 2;
			if (fNBinsPt > 13) {
				cout << "You have chosen to have more than 13 bins for Eta, this is not possible, it will be reduced to 13" << endl;
				fNBinsPt 				= 13;
			}
			for (Int_t i = 0; i < fNBinsPt+2; i++) {
				fBinsPt[i] 				= fBinsEta13TeVPtDCA[i];
			}
			fExampleBin 			        	= 3;
			nIterBGFit 					= 12;
			optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing3";
			optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing5";
			optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing7";
			fMaxYFracBGOverIntHist                          = 30;
			Int_t nBinsPlot 			= fColumn*fRow -1;
			if (fNBinsPt-fStartPtBin > nBinsPlot) fColumn++;
			nBinsPlot 					= fColumn*fRow -1;
			if (fNBinsPt-fStartPtBin > nBinsPlot) fRow++;
		}
                else if (fEnergyFlag.CompareTo("7TeV") == 0) {
			fStartPtBin					= 1;
			fColumn 					= 5;
			fRow 						= 3;
			
			if (fNBinsPt > 12) {
				cout << "You have chosen to have more than 12 bins for Eta, this is not possible, it will be reduced to 12" << endl;
				fNBinsPt 				= 12;
			}
			for (Int_t i = 0; i < fNBinsPt+2; i++) {
				fBinsPt[i] 				= fBinsEta7TeVPt[i];
			}
			fExampleBin 				= 6;
			nIterBGFit 					= 13;
		} else if (fEnergyFlag.CompareTo("8TeV") == 0) {
			fStartPtBin					= 1;
			fColumn 					= 5;
			fRow 						= 3;
			
			if (fNBinsPt > 12) {
				cout << "You have chosen to have more than 12 bins for Eta, this is not possible, it will be reduced to 12" << endl;
				fNBinsPt 				= 12;
			}
			for (Int_t i = 0; i < fNBinsPt+2; i++) {
				fBinsPt[i] 				= fBinsEta8TeVPt[i];
			}
			fExampleBin 				= 6;
			nIterBGFit 					= 12;
			optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing5";
			optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing3";
			optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing7";
		} else if (fEnergyFlag.CompareTo("2.76TeV") == 0) {
			fStartPtBin					= 1;
			fColumn 					= 4;
			fRow 						= 4;

			if (fNBinsPt > 14) {
				cout << "You have chosen to have more than 15 bins for Eta, this is not possible, it will be reduced to 15" << endl;
				fNBinsPt 				= 14;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] 				= fBinsPi02760GeVPtDCA[i];
			}
			fExampleBin 				= 4;
			nIterBGFit 					= 13;
		} else if (fEnergyFlag.CompareTo("900GeV") == 0) {
			fStartPtBin					= 1;
			fColumn 					= 2;
			fRow 						= 2;

			if (fNBinsPt > 3) {
				cout << "You have chosen to have more than 3 bins for Eta, this is not possible, it will be reduced to 3" << endl;
				fNBinsPt 				= 3;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] 				= fBinsEta900GeVPt[i];
			}
			fExampleBin 				= 2;
			nIterBGFit 					= 13;
		} else if( fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0) { 
			fStartPtBin					= 4; //otherwise usually 3
			fColumn 					= 5;
			fRow 						= 3;
			
			if (fNBinsPt > 16) {
				cout << "You have chosen to have more than 16 bins, this is not possible, it will be reduced to 16" << endl;
				fNBinsPt 				= 16;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] 				= fBinsEtaHIPtDCA[i];
			}
			
			fExampleBin = 4;
			cout << cent.Data() << endl;
			if (cent.CompareTo("60-80%")==0 || cent.CompareTo("60-70%")==0 || cent.CompareTo("70-80%")==0 || cent.CompareTo("75-90%")==0 ){
				nIterBGFit 				= 15;
				fMaxYFracBGOverIntHist	= 15;
			} else if (cent.CompareTo("50-60%")==0){
				optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing5";
				optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing3";
				optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing7";
				nIterBGFit 				= 14;
				fMaxYFracBGOverIntHist	= 12;
			} else if (cent.CompareTo("40-60%") == 0 || cent.CompareTo("40-50%")==0 ){
				optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing5";
				optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing3";
				optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing7";
				nIterBGFit 				= 16;
				fMaxYFracBGOverIntHist	= 10;
			} else if (cent.CompareTo("30-40%") == 0){
				optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing5";
				optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing3";
				optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing7";
				nIterBGFit 				= 17;      
				fMaxYFracBGOverIntHist	= 8;			
			} else if (cent.CompareTo("20-50%") == 0){
				optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing5";
				optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing3";
				optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing7";	
				nIterBGFit 				= 18;      
				fMaxYFracBGOverIntHist	= 8;
			} else if (cent.CompareTo("20-40%") == 0){
				optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing5";
				optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing3";
				optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing7";	
				nIterBGFit 				= 16;      
				fMaxYFracBGOverIntHist	= 8;
			} else if (cent.CompareTo("20-30%") == 0){
				optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing5";
				optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing3";
				optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing7";
				nIterBGFit 				= 18;      
				fMaxYFracBGOverIntHist	= 8;				
			} else {
				fMaxYFracBGOverIntHist	= 4;
				nIterBGFit 				= 21;
			}   
		} else if( fEnergyFlag.CompareTo("pPb_5.023TeV") == 0) { 
			fStartPtBin					= 1;
			fColumn 					= 3;
			fRow 						= 4;

			if (fNBinsPt > 12) {
				cout << "You have chosen to have more than 12 bins, this is not possible, it will be reduced to 12" << endl;
				fNBinsPt 				= 12;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] 				= fBinsPi0pPbPtDCA[i];
			}

			fExampleBin = 7;
			optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing5";
			optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing3";
			optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing7";
			nIterBGFit 					= 13;
			fMaxYFracBGOverIntHist		= 20;
		}
		fMesonIntRange 					= new Double_t[2];	fMesonIntRange[0] 			= 0.45;		fMesonIntRange[1]			= 0.6;
		fMesonIntRangeBG1 				= new Double_t[2];	fMesonIntRangeBG1[0] 		= 0.45;		fMesonIntRangeBG1[1]		= 0.50;
		fMesonIntRangeBG2 				= new Double_t[2];	fMesonIntRangeBG2[0] 		= 0.57;		fMesonIntRangeBG2[1]		= 0.6;
		numberInvMassBins 				= 70;
		fMesonId						= 221;
		fMesonWidthExpect 				= 0.005;
		fMesonLambdaTail 				= 0.007;
		fMesonWidthRange 				= new Double_t[2];	fMesonWidthRange[0]			= 0.002;	fMesonWidthRange[1]			= 0.010;
		fMesonLambdaTailRange 			= new Double_t[2];	fMesonLambdaTailRange[0]	= 0.0005;	fMesonLambdaTailRange[1]	= 0.026;
	} else {
		fNBinsPt 						= numberOfBins;
		fBinsPt 						= new Double_t[15];
		if (fEnergyFlag.CompareTo("13TeV") == 0) {
			fStartPtBin					= 1;
			fColumn 					= 3;
			fRow 						= 2;
			
			if (fNBinsPt > 13) {
				cout << "You have chosen to have more than 13 bins for Eta, this is not possible, it will be reduced to 13" << endl;
				fNBinsPt 				= 13;
			}
			for (Int_t i = 0; i < fNBinsPt+2; i++) {
				fBinsPt[i] 				= fBinsEta13TeVPtDCA[i];
			}
			fExampleBin 				        = 3;
			Int_t nBinsPlot 			= fColumn*fRow -1;
			if (fNBinsPt-fStartPtBin > nBinsPlot) fColumn++;
			nBinsPlot 					= fColumn*fRow -1;
			if (fNBinsPt-fStartPtBin > nBinsPlot) fRow++;
			nIterBGFit 					= 12;
			optionBGSmoothingStandard 	= "BackDecreasingWindow,BackSmoothing3";
			optionBGSmoothingVar1 		= "BackDecreasingWindow,BackSmoothing5";
			optionBGSmoothingVar2 		= "BackDecreasingWindow,BackSmoothing7";
			fMaxYFracBGOverIntHist                          = 30;
		}
		else if (fEnergyFlag.CompareTo("7TeV") == 0) {
			fStartPtBin					= 1;
			fColumn 					= 5;
			fRow 						= 3;

			if (fNBinsPt > 12) {
				cout << "You have chosen to have more than 12 bins for Eta, this is not possible, it will be reduced to 12" << endl;
				fNBinsPt 				= 12;
			}
			for (Int_t i = 0; i < fNBinsPt+2; i++) {
				fBinsPt[i] 				= fBinsEta7TeVPt[i];
			}
			fExampleBin 				= 6;
			nIterBGFit 					= 13;
		}	else if (fEnergyFlag.CompareTo("8TeV") == 0) {
			fStartPtBin					= 1;
			fColumn 					= 5;
			fRow 						= 3;
			
			if (fNBinsPt > 12) {
				cout << "You have chosen to have more than 12 bins for Eta, this is not possible, it will be reduced to 12" << endl;
				fNBinsPt 				= 12;
			}
			for (Int_t i = 0; i < fNBinsPt+2; i++) {
				fBinsPt[i] 				= fBinsEta8TeVPt[i];
			}
			fExampleBin 				= 6;
			nIterBGFit 					= 13;
		} else if (fEnergyFlag.CompareTo("2.76TeV") == 0) {
			fStartPtBin					= 1;
			fColumn 					= 3;
			fRow 						= 3;

			if (fNBinsPt > 7) {
				cout << "You have chosen to have more than 7 bins for Eta, this is not possible, it will be reduced to 9" << endl;
				fNBinsPt 				= 7;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] 				= fBinsEta2760GeVPt[i];
			}
			fExampleBin 				= 4;
			nIterBGFit 					= 13;
		} else if (fEnergyFlag.CompareTo("900GeV") == 0) {
			fStartPtBin					= 1;
			fColumn 					= 2;
			fRow 						= 2;

			if (fNBinsPt > 3) {
				cout << "You have chosen to have more than 3 bins for Eta, this is not possible, it will be reduced to 3" << endl;
				fNBinsPt 				= 3;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] 				= fBinsEta900GeVPt[i];
			}
			fExampleBin 				= 2;
			nIterBGFit 					= 13;
		} else if( fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0) { 
			fStartPtBin					= 1;
			fColumn 					= 2;
			fRow 						= 2;

			if (fNBinsPt > 4) {
				cout << "You have chosen to have more than 4 bins, this is not possible, it will be reduced to 4" << endl;
				fNBinsPt 				= 4;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] 				= fBinsEtaHIPt[i];
			}

			fExampleBin = 2;
			cout << cent.Data() << endl;
			if (cent.CompareTo("60-80%")==0){
				nIterBGFit 				= 15;
			} else if (cent.CompareTo("40-60%") == 0){
				nIterBGFit 				= 17;
			} else if (cent.CompareTo("20-40%") == 0){
				nIterBGFit 				= 19;   
			} else {
				nIterBGFit 				= 21;
			}   
		} else if( fEnergyFlag.CompareTo("pPb_5.023TeV") == 0) { 
			fStartPtBin 				= 1;
			fColumn 					= 3;
			fRow 						= 3;

			if (fNBinsPt > 8) {
				cout << "You have chosen to have more than 8 bins, this is not possible, it will be reduced to 8" << endl;
				fNBinsPt 				= 8;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] 				= fBinsEtapPbPt[i];
			}
			fExampleBin 				= 4;
			nIterBGFit 					= 15;
		}

		fMesonIntRange 					=  new Double_t[2];	fMesonIntRange[0] 			= 0.08;		fMesonIntRange[1]			= 0.2;
		fMesonIntRangeBG1 				=  new Double_t[2];	fMesonIntRangeBG1[0] 		= 0.08;		fMesonIntRangeBG1[1]		= 0.1;
		fMesonIntRangeBG2 				=  new Double_t[2];	fMesonIntRangeBG2[0] 		= 0.15;		fMesonIntRangeBG2[1]		= 0.2;	
		numberInvMassBins 				= 50;
		fMesonId						= 111;
		fMesonWidthExpect 				= 0.003;
		fMesonLambdaTail 				= 0.012;
		fMesonWidthRange 				= new Double_t[2];		fMesonWidthRange[0] 		= 0.001;	fMesonWidthRange[1]			= 0.009;
		fMesonLambdaTailRange 			= new Double_t[2]; 		fMesonLambdaTailRange[0] 	= 0.001;	fMesonLambdaTailRange[1]	= 0.02;
	}
}


//______________________________ DCAz photon under Meson Peak together with Fit and Histo estimate of BG __________________________________
void PlotDCADistPtBinWithFitAndEstimateCat(TString namePlot, TString nameCanvas, TString namePad, 
										   Int_t fRowPlot, Int_t fColumnPlot, Int_t fStartBinPtRange, Int_t fNumberPtBins, 
										   Double_t* fRangeBinsPt,  Bool_t fMonteCarloInfo,  TString textCent, TString dateDummy, Int_t category){
	TCanvas *canvasDataFit 			= new TCanvas(nameCanvas.Data(),"",2800,1800);  // gives the page size
	canvasDataFit->SetTopMargin(0.0);
	canvasDataFit->SetBottomMargin(0.0);
	canvasDataFit->SetRightMargin(0.0);
	canvasDataFit->SetLeftMargin(0.0);

	TPad * padDataFit 				= new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
	padDataFit->SetFillColor(0);
	padDataFit->GetFrame()->SetFillColor(0);
	padDataFit->SetBorderMode(0);
	padDataFit->Divide(fColumnPlot,fRowPlot,0.0,0.0);
// 	padDataFit->SetLeftMargin(0.2);
	padDataFit->Draw();

	Int_t place 					= 0;
	for(Int_t iPt = fStartBinPtRange; iPt < fNumberPtBins; iPt++){
		Double_t startPt 			= fRangeBinsPt[iPt];
		Double_t endPt 				= fRangeBinsPt[iPt+1];

		place = place + 1;                  //give the right place in the page
		if (place == fColumnPlot){
			iPt--;
			padDataFit->cd(place);
			padDataFit->cd(place)->SetTopMargin(0.12);
			padDataFit->cd(place)->SetBottomMargin(0.15);
			padDataFit->cd(place)->SetRightMargin(0.0);


			string textAlice 		= "ALICE performance";
			string textEvents		= "";
			if(fMonteCarloInfo){
				textEvents 			= "MC";
			} else {
				textEvents 			= "Data";
			}
			Double_t textHeight 	= 0.055;
			Double_t startTextX 	= 0.05;
			Double_t startTextY 	= 0.75;
			Double_t differenceText = textHeight*1.25;
			
			TLatex *alice 			= new TLatex(startTextX,(startTextY+(2*differenceText)),Form("%s",textAlice.c_str()));
			TLatex *energy 			= new TLatex(startTextX,(startTextY-differenceText),Form("%s %s",textCent.Data(),fEnergyText.Data()));
			TLatex *events 			= new TLatex(startTextX,startTextY,Form("%s: %2.1e events",textEvents.c_str(), fNEvents));
			TLatex *latexDate 		= new TLatex(startTextX,startTextY+differenceText,dateDummy.Data());
			TLatex *latexCategory 	= new TLatex(startTextX,startTextY-(3*differenceText),Form("Meson Category %i",category));

			alice->SetNDC();
			alice->SetTextColor(1);
			alice->SetTextSize(textHeight);
			alice->Draw();

			energy->SetNDC();
			energy->SetTextColor(1);
			energy->SetTextSize(textHeight);
			energy->Draw();
			
			events->SetNDC();
			events->SetTextColor(1);
			events->SetTextSize(textHeight);
			events->Draw();

			latexDate->SetNDC();
			latexDate->SetTextColor(1);
			latexDate->SetTextSize(textHeight);
			latexDate->Draw();
			
			latexCategory->SetNDC();
			latexCategory->SetTextColor(1);
			latexCategory->SetTextSize(textHeight*1.5);
			latexCategory->Draw();
			
			TLegend* legend 		= new TLegend(0.,0.1,1.,0.5);
			legend->SetTextSize(textHeight);         
			legend->SetFillColor(0);
			legend->SetLineColor(0);
			legend->SetNColumns(1);
			legend->SetMargin(0.1);
			if (fHistDCAZUnderMeson_MesonPt[category-1][fStartBinPtRange])
				legend->AddEntry(fHistDCAZUnderMeson_MesonPt[category-1][fStartBinPtRange],"Data","p");
			if (fHistDCAZUnderMesonBGEstimate_MesonPt[0][category-1][fStartBinPtRange])
				legend->AddEntry(fHistDCAZUnderMesonBGEstimate_MesonPt[0][category-1][fStartBinPtRange],"Pileup Estimate Method A","l");
			if (fHistDCAZUnderMesonBGEstimate_MesonPt[1][category-1][fStartBinPtRange])
				legend->AddEntry(fHistDCAZUnderMesonBGEstimate_MesonPt[1][category-1][fStartBinPtRange],"Pileup Estimate Method C","l");
			if (fHistDCAZUnderMesonBGEstimate_MesonPt[2][category-1][fStartBinPtRange])
				legend->AddEntry(fHistDCAZUnderMesonBGEstimate_MesonPt[2][category-1][fStartBinPtRange],"Pileup Estimate Method D","l");
			legend->Draw();

		} else {
			padDataFit->cd(place);
			padDataFit->cd(place)->SetTopMargin(0.12);
			padDataFit->cd(place)->SetBottomMargin(0.15);
			padDataFit->cd(place)->SetRightMargin(0.0);

			padDataFit->cd(place)->SetLogy();
	
			int remaining = (place-1)%fColumnPlot;
			if (remaining > 0) padDataFit->cd(place)->SetLeftMargin(0.15);
			else padDataFit->cd(place)->SetLeftMargin(0.25);

			
			if (fEnergyFlag.CompareTo("PbPb_2.76TeV")==0){
				if (fHistDCAZUnderMeson_MesonPt[category-1][iPt]){
					DrawGammaDCAHisto( fHistDCAZUnderMeson_MesonPt[category-1][iPt],
									   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
									   "dca_{z} #gamma (cm)", Form("d(dca_{z}), cat %i", category),
									   -10, 10, 0, 0.5e5);
				}
			} else {
				if (fHistDCAZUnderMeson_MesonPt[category-1][iPt]){
					DrawGammaDCAHisto( fHistDCAZUnderMeson_MesonPt[category-1][iPt],
									   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
									   "dca_{z} #gamma (cm)", Form("d(dca_{z}), cat %i", category),
									   -10, 10, 0, 0.5e4);
				}
			}
			if (fHistDCAZUnderMesonBGEstimate_MesonPt[0][category-1][iPt]){
				DrawGammaDCAHisto( fHistDCAZUnderMesonBGEstimate_MesonPt[0][category-1][iPt],
								   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
								   "dca_{z} #gamma (cm)", Form("d(dca_{z}), cat %i", category),
								   -fMaxDcaZPhoton, fMaxDcaZPhoton, 1, 0.5e4);
			}
			if (fHistDCAZUnderMesonBGEstimate_MesonPt[1][category-1][iPt]){
				DrawGammaDCAHisto( fHistDCAZUnderMesonBGEstimate_MesonPt[1][category-1][iPt],
								   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
								   "dca_{z} #gamma (cm)", Form("d(dca_{z}), cat %i", category),
								   -fMaxDcaZPhoton, fMaxDcaZPhoton, 2, 0.5e4);
			} 
			if (fHistDCAZUnderMesonBGEstimate_MesonPt[2][category-1][iPt]){
				DrawGammaDCAHisto( fHistDCAZUnderMesonBGEstimate_MesonPt[2][category-1][iPt],
								   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
								   "dca_{z} #gamma (cm)", Form("d(dca_{z}), cat %i", category),
								   -fMaxDcaZPhoton, fMaxDcaZPhoton, 3, 0.5e4);
			} 
		}
	}
	canvasDataFit->Print(namePlot.Data());
	delete padDataFit;
	delete canvasDataFit;
}


void PlotInvMassPtBinCat(TString namePlot, TString nameCanvas, TString namePad, 
						 Int_t fRowPlot, Int_t fColumnPlot, Int_t fStartBinPtRange, Int_t fNumberPtBins, 
						 Double_t* fRangeBinsPt,  Bool_t fMonteCarloInfo,  TString textCent, TString dateDummy, Int_t category){
	TCanvas *canvasDataFit 			= new TCanvas(nameCanvas.Data(),"",2800,1800);  // gives the page size
	canvasDataFit->SetTopMargin(0.0);
	canvasDataFit->SetBottomMargin(0.0);
	canvasDataFit->SetRightMargin(0.0);
	canvasDataFit->SetLeftMargin(0.0);

	TPad * padDataFit 				= new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
	padDataFit->SetFillColor(0);
	padDataFit->GetFrame()->SetFillColor(0);
	padDataFit->SetBorderMode(0);
	padDataFit->Divide(fColumnPlot,fRowPlot,0.0,0.0);
// 	padDataFit->SetLeftMargin(0.2);
	padDataFit->Draw();

	Int_t place = 0;
	for(Int_t iPt = fStartBinPtRange; iPt < fNumberPtBins; iPt++){
		Double_t startPt 			= fRangeBinsPt[iPt];
		Double_t endPt 				= fRangeBinsPt[iPt+1];

		place = place + 1;                  //give the right place in the page
		if (place == fColumnPlot){
			iPt--;
			padDataFit->cd(place);
			padDataFit->cd(place)->SetTopMargin(0.12);
			padDataFit->cd(place)->SetBottomMargin(0.15);
			padDataFit->cd(place)->SetRightMargin(0.0);


			string textAlice 		= "ALICE performance";
			string textEvents		= "";
			if(fMonteCarloInfo){
				textEvents 			= "MC";
			} else {
				textEvents 			= "Data";
			}
			Double_t textHeight 	= 0.055;
			Double_t startTextX 	= 0.05;
			Double_t startTextY 	= 0.75;
			Double_t differenceText = textHeight*1.25;
	
			TLatex *alice 			= new TLatex(startTextX,(startTextY+(2*differenceText)),Form("%s",textAlice.c_str()));
			TLatex *energy 			= new TLatex(startTextX,(startTextY-differenceText),Form("%s %s",textCent.Data(),fEnergyText.Data()));
			TLatex *events 			= new TLatex(startTextX,startTextY,Form("%s: %2.1e events",textEvents.c_str(), fNEvents));
			TLatex *latexDate 		= new TLatex(startTextX,startTextY+differenceText,dateDummy.Data());
			TLatex *latexCategory 	= new TLatex(startTextX,startTextY-(3*differenceText),Form("Meson Category %i",category));

			alice->SetNDC();
			alice->SetTextColor(1);
			alice->SetTextSize(textHeight);
			alice->Draw();

			energy->SetNDC();
			energy->SetTextColor(1);
			energy->SetTextSize(textHeight);
			energy->Draw();
			
			events->SetNDC();
			events->SetTextColor(1);
			events->SetTextSize(textHeight);
			events->Draw();

			latexDate->SetNDC();
			latexDate->SetTextColor(1);
			latexDate->SetTextSize(textHeight);
			latexDate->Draw();
			
			latexCategory->SetNDC();
			latexCategory->SetTextColor(1);
			latexCategory->SetTextSize(textHeight*1.5);
			latexCategory->Draw();
			
		} else {
			padDataFit->cd(place);
			padDataFit->cd(place)->SetTopMargin(0.12);
			padDataFit->cd(place)->SetBottomMargin(0.15);
			padDataFit->cd(place)->SetRightMargin(0.0);

			int remaining = (place-1)%fColumnPlot;
			if (remaining > 0) padDataFit->cd(place)->SetLeftMargin(0.15);
			else padDataFit->cd(place)->SetLeftMargin(0.25);

			if (fHistInvMass_MesonPt[category-1][iPt]){
				DrawGammaInvMassHisto( fHistInvMass_MesonPt[category-1][iPt],
									   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
									   "M_{#gamma#gamma} (cm)", Form("d(M_{#gamma#gamma}), cat %i", category),
									   fMesonIntRange[0], fMesonIntRange[1], 0, 0);
			} 
		}
	}
	canvasDataFit->Print(namePlot.Data());
	delete padDataFit;
	delete canvasDataFit;
}

//______________________________ DCAz photon under Meson Peak together with Fit and Histo estimate of BG __________________________________
void PlotDCADistPtBinWithMCSplitCat(TString namePlot, TString nameCanvas, TString namePad,  
									Int_t fRowPlot, Int_t fColumnPlot, Int_t fStartBinPtRange, Int_t fNumberPtBins, 
									Double_t* fRangeBinsPt,  Bool_t fMonteCarloInfo,  TString textCent, TString dateDummy, Int_t category){
	TCanvas *canvasDataFit 			= new TCanvas(nameCanvas.Data(),"",2800,1800);  // gives the page size
	canvasDataFit->SetTopMargin(0.0);
	canvasDataFit->SetBottomMargin(0.0);
	canvasDataFit->SetRightMargin(0.0);
	canvasDataFit->SetLeftMargin(0.0);

	TPad * padDataFit 				= new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
	padDataFit->SetFillColor(0);
	padDataFit->GetFrame()->SetFillColor(0);
	padDataFit->SetBorderMode(0);
	padDataFit->Divide(fColumnPlot,fRowPlot,0.0,0.0);
// 	padDataFit->SetLeftMargin(0.2);
	padDataFit->Draw();

	Int_t place = 0;
	for(Int_t iPt = fStartBinPtRange; iPt < fNumberPtBins; iPt++){
		Double_t startPt 			= fRangeBinsPt[iPt];
		Double_t endPt 				= fRangeBinsPt[iPt+1];

		place = place + 1;                  //give the right place in the page
		if (place == fColumnPlot){
			iPt--;
			padDataFit->cd(place);
			padDataFit->cd(place)->SetTopMargin(0.12);
			padDataFit->cd(place)->SetBottomMargin(0.15);
			padDataFit->cd(place)->SetRightMargin(0.0);


			string textAlice 		= "ALICE performance";
			string textEvents		= "";
			if(fMonteCarloInfo){
				textEvents = "MC";	
			} else {
				textEvents = "Data";
			}
			Double_t textHeight 	= 0.055;
			Double_t startTextX 	= 0.05;
			Double_t startTextY 	= 0.75;
			Double_t differenceText = textHeight*1.25;
		
			TLatex *alice 			= new TLatex(startTextX, (startTextY+(2*differenceText)), Form("%s", textAlice.c_str()));
			TLatex *energy 			= new TLatex(startTextX, (startTextY-differenceText), Form("%s %s", textCent.Data(), fEnergyText.Data()));
			TLatex *events 			= new TLatex(startTextX, startTextY, Form("%s: %2.1e events", textEvents.c_str(), fNEvents));
			TLatex *latexDate 		= new TLatex(startTextX, startTextY+differenceText, dateDummy.Data());
			TLatex *latexCategory 	= new TLatex(startTextX, startTextY-(3*differenceText), Form("Meson Category %i", category));

			alice->SetNDC();
			alice->SetTextColor(1);
			alice->SetTextSize(textHeight);
			alice->Draw();

			energy->SetNDC();
			energy->SetTextColor(1);
			energy->SetTextSize(textHeight);
			energy->Draw();
			
			events->SetNDC();
			events->SetTextColor(1);
			events->SetTextSize(textHeight);
			events->Draw();

			latexDate->SetNDC();
			latexDate->SetTextColor(1);
			latexDate->SetTextSize(textHeight);
			latexDate->Draw();
			
			latexCategory->SetNDC();
			latexCategory->SetTextColor(1);
			latexCategory->SetTextSize(textHeight*1.5);
			latexCategory->Draw();

			TLegend* legendMCSplit 	= new TLegend(0.,0.1,1.,0.5);
			legendMCSplit->SetTextSize(textHeight);         
			legendMCSplit->SetFillColor(0);
			legendMCSplit->SetLineColor(0);
			legendMCSplit->SetMargin(0.1);
			legendMCSplit->SetNColumns(2);
			if (fHistDCAZUnderMeson_MesonPt[category-1][fStartBinPtRange])
				legendMCSplit->AddEntry(fHistDCAZUnderMeson_MesonPt[category-1][fStartBinPtRange], "Total", "p");
			if (fHistDCAZUnderMesonBGEstimate_MesonPt[0][category-1][fStartBinPtRange])
				legendMCSplit->AddEntry(fHistDCAZUnderMesonBGEstimate_MesonPt[0][category-1][fStartBinPtRange], "BG Estimate", "l");
			if (fHistDCAZTruePrimaryMesonGammaGamma_MesonPt[category-1][fStartBinPtRange])
				legendMCSplit->AddEntry(fHistDCAZTruePrimaryMesonGammaGamma_MesonPt[category-1][fStartBinPtRange], Form("Prim. %s",fMesonType.Data()), "l");
			if (fHistDCAZTruePrimaryMesonDalitz_MesonPt[category-1][fStartBinPtRange])
				legendMCSplit->AddEntry(fHistDCAZTruePrimaryMesonDalitz_MesonPt[category-1][fStartBinPtRange], Form("Dalitz %s",fMesonType.Data()), "l");
			if (fHistDCAZTrueSecondaryMesonFromK0s_MesonPt[category-1][fStartBinPtRange] && fHistDCAZTrueSecondaryMesonFromK0s_MesonPt[category-1][fStartBinPtRange]->GetEntries() > 0)
				legendMCSplit->AddEntry(fHistDCAZTrueSecondaryMesonFromK0s_MesonPt[category-1][fStartBinPtRange], "Sec. #pi^{0} from K^{0}_{s}", "l");
			if (fHistDCAZTrueSecondaryMesonFromSomething_MesonPt[category-1][fStartBinPtRange] && fHistDCAZTrueSecondaryMesonFromSomething_MesonPt[category-1][fStartBinPtRange]->GetEntries() > 0)
				legendMCSplit->AddEntry(fHistDCAZTrueSecondaryMesonFromSomething_MesonPt[category-1][fStartBinPtRange],"Sec. #pi^{0} from X", "l");
			if (fHistDCAZTrueBackground_MesonPt[category-1][fStartBinPtRange])
				legendMCSplit->AddEntry(fHistDCAZTrueBackground_MesonPt[category-1][fStartBinPtRange], "#gamma#gamma BG", "l");
			if (fHistDCAZGarbage_MesonPt[category-1][fStartBinPtRange])
				legendMCSplit->AddEntry(fHistDCAZGarbage_MesonPt[category-1][fStartBinPtRange], "garbage", "l");
			legendMCSplit->Draw();
		} else {
			padDataFit->cd(place);
			padDataFit->cd(place)->SetTopMargin(0.12);
			padDataFit->cd(place)->SetBottomMargin(0.15);
			padDataFit->cd(place)->SetRightMargin(0.0);

			padDataFit->cd(place)->SetLogy();
			int remaining = (place-1)%fColumnPlot;
			if (remaining > 0) padDataFit->cd(place)->SetLeftMargin(0.15);
			else padDataFit->cd(place)->SetLeftMargin(0.25);

			
			if (fEnergyFlag.CompareTo("PbPb_2.76TeV")==0){
				if (fHistDCAZUnderMeson_MesonPt[category-1][iPt]){
					DrawGammaDCAHisto( fHistDCAZUnderMeson_MesonPt[category-1][iPt],
									   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
									   "dca_{z} #gamma (cm)", Form("d(dca_{z}), cat %i", category),
									   -10, 10, 0, 0.5e5);
				}
			} else {
				if (fHistDCAZUnderMeson_MesonPt[category-1][iPt]){
					DrawGammaDCAHisto( fHistDCAZUnderMeson_MesonPt[category-1][iPt],
									   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
									   "dca_{z} #gamma (cm)", Form("d(dca_{z}), cat %i",category),
									   -10, 10, 0, 0.5e4);
				}
			}
			
			if (fHistDCAZUnderMesonBGEstimate_MesonPt[0][category-1][iPt]){
				DrawGammaDCAHisto( fHistDCAZUnderMesonBGEstimate_MesonPt[0][category-1][iPt],
								   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
								   "dca_{z} #gamma (cm)", Form("d(dca_{z}), cat %i", category),
								   -fMaxDcaZPhoton, fMaxDcaZPhoton, 1, 0.5e4);
			}
			if (fHistDCAZTruePrimaryMesonGammaGamma_MesonPt[category-1][iPt]){
				DrawGammaDCAHisto( fHistDCAZTruePrimaryMesonGammaGamma_MesonPt[category-1][iPt],
								   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
								   "dca_{z} #gamma (cm)", Form("d(dca_{z}), cat %i", category),
								   -fMaxDcaZPhoton, fMaxDcaZPhoton, 2, 0.5e4);
			} 
			if (fHistDCAZTruePrimaryMesonDalitz_MesonPt[category-1][iPt]){
				DrawGammaDCAHisto( fHistDCAZTruePrimaryMesonDalitz_MesonPt[category-1][iPt],
								   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
								   "dca_{z} #gamma (cm)", Form("d(dca_{z}), cat %i", category),
								   -fMaxDcaZPhoton, fMaxDcaZPhoton, 3, 0.5e4);
			}
			if (fHistDCAZTrueSecondaryMesonFromK0s_MesonPt[category-1][iPt]){
				DrawGammaDCAHisto( fHistDCAZTrueSecondaryMesonFromK0s_MesonPt[category-1][iPt],
								   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
								   "dca_{z} #gamma (cm)", Form("d(dca_{z}), cat %i", category),
								   -fMaxDcaZPhoton, fMaxDcaZPhoton, 4, 0.5e4);
			}
			if (fHistDCAZTrueSecondaryMesonFromSomething_MesonPt[category-1][iPt]){
				DrawGammaDCAHisto( fHistDCAZTrueSecondaryMesonFromSomething_MesonPt[category-1][iPt],
								   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
								   "dca_{z} #gamma (cm)", Form("d(dca_{z}), cat %i", category),
								   -fMaxDcaZPhoton, fMaxDcaZPhoton, 5, 0.5e4);
			}
			if (fHistDCAZTrueBackground_MesonPt[category-1][iPt]){
				DrawGammaDCAHisto( fHistDCAZTrueBackground_MesonPt[category-1][iPt],
								   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
								   "dca_{z} #gamma (cm)", Form("d(dca_{z}), cat %i", category),
								   -fMaxDcaZPhoton, fMaxDcaZPhoton, 6, 0.5e4);
			}
			if (fHistDCAZGarbage_MesonPt[category-1][iPt]){
				DrawGammaDCAHisto( fHistDCAZGarbage_MesonPt[category-1][iPt],
								   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
								   "dca_{z} #gamma (cm)", Form("d(dca_{z}), cat %i", category),
								   -fMaxDcaZPhoton, fMaxDcaZPhoton, 7, 0.5e4);
			}
		}
	}
	canvasDataFit->Print(namePlot.Data());
	delete padDataFit;
	delete canvasDataFit;
}

//______________________________ DCAz photon under Meson Peak together with Fit and Histo estimate of BG __________________________________
void PlotDCADistPtBinWithFitAndEstimate(TString namePlot, TString nameCanvas, TString namePad, 
										Int_t fRowPlot, Int_t fColumnPlot, Int_t fStartBinPtRange, Int_t fNumberPtBins, 
										Double_t* fRangeBinsPt,  Bool_t fMonteCarloInfo,  TString textCent, TString dateDummy){
	TCanvas *canvasDataFit 			= new TCanvas(nameCanvas.Data(),"",2800,1800);  // gives the page size
	canvasDataFit->SetTopMargin(0.0);
	canvasDataFit->SetBottomMargin(0.0);
	canvasDataFit->SetRightMargin(0.0);
	canvasDataFit->SetLeftMargin(0.0);

	TPad * padDataFit 				= new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
	padDataFit->SetFillColor(0);
	padDataFit->GetFrame()->SetFillColor(0);
	padDataFit->SetBorderMode(0);
	padDataFit->Divide(fColumnPlot,fRowPlot,0.0,0.0);
// 	padDataFit->SetLeftMargin(0.2);
	padDataFit->Draw();

	Int_t place = 0;
	for(Int_t iPt = fStartBinPtRange; iPt < fNumberPtBins; iPt++){
		Double_t startPt 			= fRangeBinsPt[iPt];
		Double_t endPt 				= fRangeBinsPt[iPt+1];

		place = place + 1;                  //give the right place in the page
		if (place == fColumnPlot){
			iPt--;
			padDataFit->cd(place);
			padDataFit->cd(place)->SetTopMargin(0.12);
			padDataFit->cd(place)->SetBottomMargin(0.15);
			padDataFit->cd(place)->SetRightMargin(0.0);


			string textAlice 			= "ALICE performance";
			string textEvents			= "";
			if(fMonteCarloInfo){
					textEvents 			= "MC";
			} else {
				textEvents 				= "Data";
			}
			Double_t textHeight 	= 0.055;
			Double_t startTextX 	= 0.05;
			Double_t startTextY 	= 0.75;
			Double_t differenceText = textHeight*1.25;
		
			TLatex *alice 			= new TLatex(startTextX, (startTextY+(2*differenceText)), Form("%s", textAlice.c_str()));
			TLatex *energy 			= new TLatex(startTextX, (startTextY-differenceText), Form("%s %s", textCent.Data(), fEnergyText.Data()));
			TLatex *events 			= new TLatex(startTextX, startTextY, Form("%s: %2.1e events", textEvents.c_str(), fNEvents));
			TLatex *latexDate 		= new TLatex(startTextX, startTextY+differenceText, dateDummy.Data());
			TLatex *latexCategory 	= new TLatex(startTextX, startTextY-(3*differenceText), Form("All Meson Categories"));

			alice->SetNDC();
			alice->SetTextColor(1);
			alice->SetTextSize(textHeight);
			alice->Draw();

			energy->SetNDC();
			energy->SetTextColor(1);
			energy->SetTextSize(textHeight);
			energy->Draw();
			
			events->SetNDC();
			events->SetTextColor(1);
			events->SetTextSize(textHeight);
			events->Draw();

			latexDate->SetNDC();
			latexDate->SetTextColor(1);
			latexDate->SetTextSize(textHeight);
			latexDate->Draw();
			
			latexCategory->SetNDC();
			latexCategory->SetTextColor(1);
			latexCategory->SetTextSize(textHeight*1.5);
			latexCategory->Draw();
			
			TLegend* legend 		= new TLegend(0., 0.1, 1., 0.5);
			legend->SetTextSize(textHeight);         
			legend->SetFillColor(0);
			legend->SetLineColor(0);
			legend->SetNColumns(1);
			legend->SetMargin(0.1);
			if (fHistDCAZUnderMeson_MesonPt_AllCat[fStartBinPtRange])
				legend->AddEntry(fHistDCAZUnderMeson_MesonPt_AllCat[fStartBinPtRange], "Data", "p");
			if (fHistDCAZUnderMesonBGEstimate_MesonPt_AllCat[fStartBinPtRange])
				legend->AddEntry(fHistDCAZUnderMesonBGEstimate_MesonPt_AllCat[fStartBinPtRange], "Pileup Estimate Method A", "l");
			if (fFitDCAZUnderMesonBGEstimate_MesonPt_AllCat[fStartBinPtRange])
				legend->AddEntry(fFitDCAZUnderMesonBGEstimate_MesonPt_AllCat[fStartBinPtRange], "Pileup Estimate Method B", "l");
			legend->Draw();
	
		} else {
			padDataFit->cd(place);
			padDataFit->cd(place)->SetTopMargin(0.12);
			padDataFit->cd(place)->SetBottomMargin(0.15);
			padDataFit->cd(place)->SetRightMargin(0.0);

			padDataFit->cd(place)->SetLogy();

			int remaining = (place-1)%fColumnPlot;
			if (remaining > 0) padDataFit->cd(place)->SetLeftMargin(0.15);
			else padDataFit->cd(place)->SetLeftMargin(0.25);
			
			if (fEnergyFlag.CompareTo("PbPb_2.76TeV")==0){
				if (fHistDCAZUnderMeson_MesonPt_AllCat[iPt]){
					DrawGammaDCAHisto( fHistDCAZUnderMeson_MesonPt_AllCat[iPt],
									   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
									   "dca_{z} #gamma (cm)", Form("d(dca_{z}), all cat"),
									   -10, 10, 0, 0.5e5);
				}
			} else {
				if (fHistDCAZUnderMeson_MesonPt_AllCat[iPt]){
					DrawGammaDCAHisto( fHistDCAZUnderMeson_MesonPt_AllCat[iPt],
									   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
									   "dca_{z} #gamma (cm)", Form("d(dca_{z}), all cat"),
									   -10, 10, 0, 0.5e4);
				}
			}
			if (fHistDCAZUnderMesonBGEstimate_MesonPt_AllCat[iPt]){
				DrawGammaDCAHisto( fHistDCAZUnderMesonBGEstimate_MesonPt_AllCat[iPt],
								   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
								   "dca_{z} #gamma (cm)", Form("d(dca_{z}), all cat"),
								   -fMaxDcaZPhoton, fMaxDcaZPhoton, 1, 0.5e4);
			}
			if (fFitDCAZUnderMesonBGEstimate_MesonPt_AllCat[iPt]){
				fFitDCAZUnderMesonBGEstimate_MesonPt_AllCat[iPt]->SetRange(-fMaxDcaZPhoton, fMaxDcaZPhoton);
				fFitDCAZUnderMesonBGEstimate_MesonPt_AllCat[iPt]->SetLineColor(kCyan+3);
				fFitDCAZUnderMesonBGEstimate_MesonPt_AllCat[iPt]->SetLineWidth(0.99);
				fFitDCAZUnderMesonBGEstimate_MesonPt_AllCat[iPt]->DrawCopy("same");
			} 
		}
	}
	canvasDataFit->Print(namePlot.Data());
	delete padDataFit;
	delete canvasDataFit;
}

//______________________________ DCAz photon under Meson Peak together with Fit and Histo estimate of BG __________________________________
void PlotDCADistPtBinWithFitAndEstimateGoodCat(TString namePlot, TString nameCanvas, TString namePad, 
											   Int_t fRowPlot, Int_t fColumnPlot, Int_t fStartBinPtRange, Int_t fNumberPtBins, 
											   Double_t* fRangeBinsPt,  Bool_t fMonteCarloInfo,  TString textCent, TString dateDummy){
	TCanvas *canvasDataFit 			= new TCanvas(nameCanvas.Data(), "", 2800, 1800);  // gives the page size
	canvasDataFit->SetTopMargin(0.0);
	canvasDataFit->SetBottomMargin(0.0);
	canvasDataFit->SetRightMargin(0.0);
	canvasDataFit->SetLeftMargin(0.0);

	TPad * padDataFit 				= new TPad(namePad.Data(), "", 0.0, 0.0, 1., 1., 0);   // gives the size of the histo areas
	padDataFit->SetFillColor(0);
	padDataFit->GetFrame()->SetFillColor(0);
	padDataFit->SetBorderMode(0);
	padDataFit->Divide(fColumnPlot,fRowPlot,0.0,0.0);
// 	padDataFit->SetLeftMargin(0.2);
	padDataFit->Draw();

	Int_t place = 0;
	for(Int_t iPt = fStartBinPtRange; iPt < fNumberPtBins; iPt++){
		Double_t startPt 			= fRangeBinsPt[iPt];
		Double_t endPt 				= fRangeBinsPt[iPt+1];

		place = place + 1;                  //give the right place in the page
		if (place == fColumnPlot){
			iPt--;
			padDataFit->cd(place);
			padDataFit->cd(place)->SetTopMargin(0.12);
			padDataFit->cd(place)->SetBottomMargin(0.15);
			padDataFit->cd(place)->SetRightMargin(0.0);


			string textAlice 		= "ALICE performance";
			string textEvents;
			if(fMonteCarloInfo){
				textEvents 			= "MC";	
			} else {
				textEvents 			= "Data";
			}
			Double_t textHeight 	= 0.055;
			Double_t startTextX 	= 0.05;
			Double_t startTextY 	= 0.75;
			Double_t differenceText = textHeight*1.25;
			
			TLatex *alice 			= new TLatex(startTextX, (startTextY+(2*differenceText)), Form("%s",textAlice.c_str()));
			TLatex *energy 			= new TLatex(startTextX, (startTextY-differenceText), Form("%s %s",textCent.Data(), fEnergyText.Data()));
			TLatex *events 			= new TLatex(startTextX, startTextY, Form("%s: %2.1e events", textEvents.c_str(), fNEvents));
			TLatex *latexDate 		= new TLatex(startTextX, startTextY+differenceText, dateDummy.Data());
			TLatex *latexCategory 	= new TLatex(startTextX, startTextY-(3*differenceText), Form("Meson Categories 4+5+6"));

			alice->SetNDC();
			alice->SetTextColor(1);
			alice->SetTextSize(textHeight);
			alice->Draw();

			energy->SetNDC();
			energy->SetTextColor(1);
			energy->SetTextSize(textHeight);
			energy->Draw();
			
			events->SetNDC();
			events->SetTextColor(1);
			events->SetTextSize(textHeight);
			events->Draw();

			latexDate->SetNDC();
			latexDate->SetTextColor(1);
			latexDate->SetTextSize(textHeight);
			latexDate->Draw();
			
			latexCategory->SetNDC();
			latexCategory->SetTextColor(1);
			latexCategory->SetTextSize(textHeight*1.5);
			latexCategory->Draw();
			
		} else {
			padDataFit->cd(place);
			padDataFit->cd(place)->SetTopMargin(0.12);
			padDataFit->cd(place)->SetBottomMargin(0.15);
			padDataFit->cd(place)->SetRightMargin(0.0);

			padDataFit->cd(place)->SetLogy();

			int remaining = (place-1)%fColumnPlot;
			if (remaining > 0) padDataFit->cd(place)->SetLeftMargin(0.15);
			else padDataFit->cd(place)->SetLeftMargin(0.25);
			
			if (fEnergyFlag.CompareTo("PbPb_2.76TeV")==0){
				if (fHistDCAZUnderMeson_MesonPt_GoodCat[iPt]){
					DrawGammaDCAHisto( fHistDCAZUnderMeson_MesonPt_GoodCat[iPt],
									   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
									   "dca_{z} #gamma (cm)", Form("d(dca_{z}), all cat"),
									   -10, 10, 0, 0.5e5);
				}
			} else {
				if (fHistDCAZUnderMeson_MesonPt_GoodCat[iPt]){
					DrawGammaDCAHisto( fHistDCAZUnderMeson_MesonPt_GoodCat[iPt],
									   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
									   "dca_{z} #gamma (cm)", Form("d(dca_{z}), all cat"),
									   -10, 10, 0, 0.5e4);
				}
			}
		}
	}
	canvasDataFit->Print(namePlot.Data());
	delete padDataFit;
	delete canvasDataFit;
}


void PlotInvMassPtBin(TString namePlot, TString nameCanvas, TString namePad,  
					  Int_t fRowPlot, Int_t fColumnPlot, Int_t fStartBinPtRange, Int_t fNumberPtBins, 
					  Double_t* fRangeBinsPt,  Bool_t fMonteCarloInfo,  TString textCent, TString dateDummy){
	TCanvas *canvasDataFit 		= new TCanvas(nameCanvas.Data(), "", 2800, 1800);  // gives the page size
	canvasDataFit->SetTopMargin(0.0);
	canvasDataFit->SetBottomMargin(0.0);
	canvasDataFit->SetRightMargin(0.0);
	canvasDataFit->SetLeftMargin(0.0);

	TPad * padDataFit 			= new TPad(namePad.Data(), "", 0.0, 0.0, 1., 1., 0);   // gives the size of the histo areas
	padDataFit->SetFillColor(0);
	padDataFit->GetFrame()->SetFillColor(0);
	padDataFit->SetBorderMode(0);
	padDataFit->Divide(fColumnPlot,fRowPlot,0.0,0.0);
	padDataFit->Draw();

	Int_t place = 0;
	for(Int_t iPt = fStartBinPtRange; iPt < fNumberPtBins; iPt++){
		Double_t startPt 		= fRangeBinsPt[iPt];
		Double_t endPt 			= fRangeBinsPt[iPt+1];

		place = place + 1;                  //give the right place in the page
		if (place == fColumnPlot){
			iPt--;
			padDataFit->cd(place);
			padDataFit->cd(place)->SetTopMargin(0.12);
			padDataFit->cd(place)->SetBottomMargin(0.15);
			padDataFit->cd(place)->SetRightMargin(0.0);


			string textAlice 				= "ALICE performance";
			string textEvents 				= "";
			if(fMonteCarloInfo){
				textEvents 					= "MC";
			} else {
				textEvents 					= "Data";
			}
			Double_t textHeight 			= 0.055;
			Double_t startTextX 			= 0.0;
			Double_t startTextY 			= 0.75;
			Double_t differenceText 		= textHeight*1.25;
			
			TLatex *alice 					= new TLatex(startTextX, (startTextY+(2*differenceText)), Form("%s",textAlice.c_str()));
			TLatex *energy 					= new TLatex(startTextX, (startTextY-differenceText), Form("%s %s",textCent.Data(), fEnergyText.Data()));
			TLatex *events 					= new TLatex(startTextX, startTextY, Form("%s: %2.1e events", textEvents.c_str(), fNEvents));
			TLatex *latexDate 				= new TLatex(startTextX, startTextY+differenceText, dateDummy.Data());
			TLatex *latexCategory 			= new TLatex(startTextX, startTextY-(3*differenceText), Form("All Meson Categories"));

			alice->SetNDC();
			alice->SetTextColor(1);
			alice->SetTextSize(textHeight);
			alice->Draw();

			energy->SetNDC();
			energy->SetTextColor(1);
			energy->SetTextSize(textHeight);
			energy->Draw();
			
			events->SetNDC();
			events->SetTextColor(1);
			events->SetTextSize(textHeight);
			events->Draw();

			latexDate->SetNDC();
			latexDate->SetTextColor(1);
			latexDate->SetTextSize(textHeight);
			latexDate->Draw();
			
			latexCategory->SetNDC();
			latexCategory->SetTextColor(1);
			latexCategory->SetTextSize(textHeight*1.5);
			latexCategory->Draw();
			
		} else {
			padDataFit->cd(place);
			padDataFit->cd(place)->SetTopMargin(0.12);
			padDataFit->cd(place)->SetBottomMargin(0.15);
			padDataFit->cd(place)->SetRightMargin(0.0);


			int remaining = (place-1)%fColumnPlot;
			if (remaining > 0) padDataFit->cd(place)->SetLeftMargin(0.15);
			else padDataFit->cd(place)->SetLeftMargin(0.25);

			if (fHistInvMass_MesonPt_AllCat[iPt]){
				DrawGammaInvMassHisto( fHistInvMass_MesonPt_AllCat[iPt],
									   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
									   "M_{#gamma#gamma} (cm)", Form("d(M_{#gamma#gamma})"),
									   fMesonIntRange[0], fMesonIntRange[1], 0, 0);
			}
		}
	}
	canvasDataFit->Print(namePlot.Data());
	delete padDataFit;
	delete canvasDataFit;
}


//______________________________ DCAz photon under Meson Peak together with Fit and Histo estimate of BG __________________________________
void PlotDCADistPtBinWithMCSplit(TString namePlot, TString nameCanvas, TString namePad,  
								 Int_t fRowPlot, Int_t fColumnPlot, Int_t fStartBinPtRange, Int_t fNumberPtBins, 
								 Double_t* fRangeBinsPt,  Bool_t fMonteCarloInfo,  TString textCent, TString dateDummy){
	
	TCanvas *canvasDataFit 			= new TCanvas(nameCanvas.Data(), "", 2800, 1800);  // gives the page size
	canvasDataFit->SetTopMargin(0.0);
	canvasDataFit->SetBottomMargin(0.0);
	canvasDataFit->SetRightMargin(0.0);
	canvasDataFit->SetLeftMargin(0.0);

	TPad * padDataFit 				= new TPad(namePad.Data(), "", 0.0, 0.0, 1., 1., 0);   // gives the size of the histo areas
	padDataFit->SetFillColor(0);
	padDataFit->GetFrame()->SetFillColor(0);
	padDataFit->SetBorderMode(0);
	padDataFit->Divide(fColumnPlot,fRowPlot,0.0,0.0);
// 	padDataFit->SetLeftMargin(0.2);
	padDataFit->Draw();

	Int_t place = 0;
	for(Int_t iPt = fStartBinPtRange; iPt < fNumberPtBins; iPt++){
		Double_t startPt 			= fRangeBinsPt[iPt];
		Double_t endPt 				= fRangeBinsPt[iPt+1];

		place = place + 1;                  //give the right place in the page
		if (place == fColumnPlot){
			iPt--;
			padDataFit->cd(place);
			padDataFit->cd(place)->SetTopMargin(0.12);
			padDataFit->cd(place)->SetBottomMargin(0.15);
			padDataFit->cd(place)->SetRightMargin(0.0);


			string textAlice 		= "ALICE performance";
			string textEvents		= "";
			if(fMonteCarloInfo){
				textEvents 			= "MC";
			} else {
				textEvents 			= "Data";
			}
			Double_t textHeight 	= 0.055;
			Double_t startTextX 	= 0.05;
			Double_t startTextY 	= 0.75;
			Double_t differenceText = textHeight*1.25;
	
			TLatex *alice 			= new TLatex(startTextX, (startTextY+(2*differenceText)), Form("%s",textAlice.c_str()));
			TLatex *energy 			= new TLatex(startTextX, (startTextY-differenceText), Form("%s %s",textCent.Data(),fEnergyText.Data()));
			TLatex *events 			= new TLatex(startTextX, startTextY, Form("%s: %2.1e events",textEvents.c_str(), fNEvents));
			TLatex *latexDate 		= new TLatex(startTextX, startTextY+differenceText, dateDummy.Data());
			
			alice->SetNDC();
			alice->SetTextColor(1);
			alice->SetTextSize(textHeight);
			alice->Draw();

			energy->SetNDC();
			energy->SetTextColor(1);
			energy->SetTextSize(textHeight);
			energy->Draw();
			
			events->SetNDC();
			events->SetTextColor(1);
			events->SetTextSize(textHeight);
			events->Draw();

			latexDate->SetNDC();
			latexDate->SetTextColor(1);
			latexDate->SetTextSize(textHeight);
			latexDate->Draw();

			TLegend* legendMCSplit 	= new TLegend(0., 0.1, 1., 0.5);
			legendMCSplit->SetTextSize(textHeight);         
			legendMCSplit->SetFillColor(0);
			legendMCSplit->SetLineColor(0);
			legendMCSplit->SetNColumns(2);
			legendMCSplit->SetMargin(0.1);
			if (fHistDCAZUnderMeson_MesonPt_AllCat[fStartBinPtRange])legendMCSplit->AddEntry(fHistDCAZUnderMeson_MesonPt_AllCat[fStartBinPtRange], "Total","p");
// 			if (fHistDCAZUnderMesonBGEstimate_MesonPt_AllCat[0][fStartBinPtRange])
// 				legendMCSplit->AddEntry(fHistDCAZUnderMesonBGEstimate_MesonPt_AllCat[0][fStartBinPtRange], "BG Estimate","l");
			if (fHistDCAZTruePrimaryMesonGammaGamma_MesonPt_AllCat[fStartBinPtRange])
				legendMCSplit->AddEntry(fHistDCAZTruePrimaryMesonGammaGamma_MesonPt_AllCat[fStartBinPtRange], Form("Prim. %s",fMesonType.Data()),"l");
			if (fHistDCAZTruePrimaryMesonDalitz_MesonPt_AllCat[fStartBinPtRange])
				legendMCSplit->AddEntry(fHistDCAZTruePrimaryMesonDalitz_MesonPt_AllCat[fStartBinPtRange],Form("Dalitz %s",fMesonType.Data()),"l");
			if (fHistDCAZTrueSecondaryMesonFromK0s_MesonPt_AllCat[fStartBinPtRange] && fHistDCAZTrueSecondaryMesonFromK0s_MesonPt_AllCat[fStartBinPtRange]->GetEntries() > 0)
				legendMCSplit->AddEntry(fHistDCAZTrueSecondaryMesonFromK0s_MesonPt_AllCat[fStartBinPtRange],"Sec. #pi^{0} from K^{0}_{s}","l");
			if (fHistDCAZTrueSecondaryMesonFromSomething_MesonPt_AllCat[fStartBinPtRange] && fHistDCAZTrueSecondaryMesonFromSomething_MesonPt_AllCat[fStartBinPtRange]->GetEntries() > 0) 
				legendMCSplit->AddEntry(fHistDCAZTrueSecondaryMesonFromSomething_MesonPt_AllCat[fStartBinPtRange],"Sec. #pi^{0} from X","l");
			if (fHistDCAZTrueBackground_MesonPt_AllCat[fStartBinPtRange])
				legendMCSplit->AddEntry(fHistDCAZTrueBackground_MesonPt_AllCat[fStartBinPtRange],"#gamma#gamma BG","l");
			if (fHistDCAZGarbage_MesonPt_AllCat[fStartBinPtRange])
				legendMCSplit->AddEntry(fHistDCAZGarbage_MesonPt_AllCat[fStartBinPtRange],"garbage","l");
			legendMCSplit->Draw();
		
			
		} else {
			padDataFit->cd(place);
			padDataFit->cd(place)->SetTopMargin(0.12);
			padDataFit->cd(place)->SetBottomMargin(0.15);
			padDataFit->cd(place)->SetRightMargin(0.0);

			padDataFit->cd(place)->SetLogy();
			
			int remaining = (place-1)%fColumnPlot;
			if (remaining > 0) padDataFit->cd(place)->SetLeftMargin(0.15);
			else padDataFit->cd(place)->SetLeftMargin(0.25);
			
			if (fEnergyFlag.CompareTo("PbPb_2.76TeV")==0){
				if (fHistDCAZUnderMeson_MesonPt_AllCat[iPt]){
					DrawGammaDCAHisto( fHistDCAZUnderMeson_MesonPt_AllCat[iPt], 
									   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
									   "dca_{z} #gamma (cm)", Form("d(dca_{z})"), 
									   -10, 10, 0, 0.5e5);
				}
			} else {
				if (fHistDCAZUnderMeson_MesonPt_AllCat[iPt]){
					DrawGammaDCAHisto( fHistDCAZUnderMeson_MesonPt_AllCat[iPt],
									   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
									   "dca_{z} #gamma (cm)", Form("d(dca_{z})"),
									   -10, 10, 0, 0.5e4);
				}
			}
			
// 			if (fHistDCAZUnderMesonBGEstimate_MesonPt_AllCat[0][iPt]){
// 				DrawGammaDCAHisto( fHistDCAZUnderMesonBGEstimate_MesonPt_AllCat[0][iPt],
// 								   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
// 								   "dca_{z} #gamma (cm)", Form("d(dca_{z})"),
// 								   -fMaxDcaZPhoton, fMaxDcaZPhoton, 1, 0.5e4);
// 			}
			if (fHistDCAZTruePrimaryMesonGammaGamma_MesonPt_AllCat[iPt]){
				DrawGammaDCAHisto( fHistDCAZTruePrimaryMesonGammaGamma_MesonPt_AllCat[iPt],
								   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c",startPt, fMesonType.Data(), endPt),
								   "dca_{z} #gamma (cm)", Form("d(dca_{z})"),
								   -fMaxDcaZPhoton, fMaxDcaZPhoton, 2, 0.5e4);
			} 
			if (fHistDCAZTruePrimaryMesonDalitz_MesonPt_AllCat[iPt]){
				DrawGammaDCAHisto( fHistDCAZTruePrimaryMesonDalitz_MesonPt_AllCat[iPt],
								   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
								   "dca_{z} #gamma (cm)", Form("d(dca_{z})"),
								   -fMaxDcaZPhoton, fMaxDcaZPhoton, 3, 0.5e4);
			}
			if (fHistDCAZTrueSecondaryMesonFromK0s_MesonPt_AllCat[iPt]){
				DrawGammaDCAHisto( fHistDCAZTrueSecondaryMesonFromK0s_MesonPt_AllCat[iPt],
								   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
								   "dca_{z} #gamma (cm)", Form("d(dca_{z})"),
								   -fMaxDcaZPhoton, fMaxDcaZPhoton, 4, 0.5e4);
			}
			if (fHistDCAZTrueSecondaryMesonFromSomething_MesonPt_AllCat[iPt]){
				DrawGammaDCAHisto( fHistDCAZTrueSecondaryMesonFromSomething_MesonPt_AllCat[iPt],
								   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
								   "dca_{z} #gamma (cm)", Form("d(dca_{z})"),
								   -fMaxDcaZPhoton, fMaxDcaZPhoton, 5, 0.5e4);
			}
			if (fHistDCAZTrueBackground_MesonPt_AllCat[iPt]){
				DrawGammaDCAHisto( fHistDCAZTrueBackground_MesonPt_AllCat[iPt],
								   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
								   "dca_{z} #gamma (cm)", Form("d(dca_{z})"),
								   -fMaxDcaZPhoton, fMaxDcaZPhoton, 6, 0.5e4);
			}
			if (fHistDCAZGarbage_MesonPt_AllCat[iPt]){
				DrawGammaDCAHisto( fHistDCAZGarbage_MesonPt_AllCat[iPt],
								   Form("%3.2f GeV/c < p_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
								   "dca_{z} #gamma (cm)", Form("d(dca_{z})"),
								   -fMaxDcaZPhoton, fMaxDcaZPhoton, 7, 0.5e4);
			}
		}
	}
	canvasDataFit->Print(namePlot.Data());
	delete padDataFit;
	delete canvasDataFit;
}




void DrawGammaDCAHisto( TH1* histo1,
						TString Title, TString XTitle, TString YTitle,
						Float_t xMin, Float_t xMax,Int_t bck, Double_t numberOfOrders ) {
	
	histo1->GetXaxis()->SetRangeUser(xMin, xMax);
	if (bck != 1 || bck != 2){
		Double_t yMin = 0;
		Double_t yMax = 0;
		for (Int_t i = histo1->GetXaxis()->FindBin(xMin); i < histo1->GetXaxis()->FindBin(xMax); i++){
			if (histo1->GetBinContent(i) > yMax){
				yMax = histo1->GetBinContent(i);
			}
		}   
		yMin = yMax/numberOfOrders;
		if (yMin < 1e-1) yMin = 1e-1;
		histo1->GetYaxis()->SetRangeUser(yMin, 2*yMax);

	}

	if(XTitle.Length() > 0){
		histo1->SetXTitle(XTitle.Data());
	}
	if(YTitle.Length() > 0){
		histo1->SetYTitle(YTitle.Data());
	}
	histo1->GetYaxis()->SetLabelSize(0.02);
	histo1->GetYaxis()->SetTitleSize(0.025);
	histo1->GetYaxis()->SetDecimals();
	histo1->GetYaxis()->SetTitleOffset(0.5);
	histo1->GetXaxis()->SetTitleSize(0.025);
	histo1->GetXaxis()->SetLabelSize(0.02);
	histo1->SetMarkerStyle(20);
	histo1->SetMarkerColor(1);
	histo1->SetLineColor(1);
	histo1->SetLineWidth(0.5);
	histo1->SetMarkerSize(0.2);
	histo1->SetTitleOffset(1.2,"xy");      
	histo1->SetTitleSize(0.05,"xy");    
	histo1->GetYaxis()->SetLabelSize(0.05);
	histo1->GetXaxis()->SetLabelSize(0.05);
	histo1->GetXaxis()->SetNdivisions(507,kTRUE);
	if(bck==1){
		histo1->SetLineStyle(1);      
		histo1->SetLineColor(4);
		histo1->SetMarkerColor(4);
		histo1->SetMarkerStyle(24);
		histo1->SetLineWidth(0.9);
		histo1->DrawCopy("hist,same");
	}
	if(bck==2){
		histo1->SetLineStyle(1);      
		histo1->SetLineColor(kRed+2);
		histo1->SetMarkerColor(kRed+2);
		histo1->SetLineWidth(0.9);
		histo1->DrawCopy("hist,same");
	}
	if(bck==3){
		histo1->SetLineStyle(1);      
		histo1->SetLineColor(kGreen+2);
		histo1->SetMarkerColor(kGreen+2);
		histo1->SetLineWidth(0.9);
		histo1->DrawCopy("hist,same");
	}
	if(bck==4){
		histo1->SetLineStyle(1);      
		histo1->SetLineColor(807);
		histo1->SetMarkerColor(807);
		histo1->SetLineWidth(0.9);
		histo1->DrawCopy("hist,same");
	}
	if(bck==5){
		histo1->SetLineStyle(1);      
		histo1->SetLineColor(kViolet+2);
		histo1->SetMarkerColor(kViolet+2);
		histo1->SetLineWidth(0.9);
		histo1->DrawCopy("hist,same");
	}

	if(bck==6){
		histo1->SetLineStyle(1);      
		histo1->SetLineColor(kPink+2);
		histo1->SetMarkerColor(kPink+2);
		histo1->SetLineWidth(0.9);
		histo1->DrawCopy("hist,same");
	}
	if(bck==7){
		histo1->SetLineStyle(1);      
		histo1->SetLineColor(kCyan+2);
		histo1->SetMarkerColor(kCyan+2);
		histo1->SetLineWidth(0.9);
		histo1->DrawCopy("hist,same");
	}

	if (bck == 0){
		histo1->SetTitle("");
		histo1->DrawCopy("e1,p");   
		if(Title.Length() > 0){
			histo1->SetTitle("");
			TLatex *alice = new TLatex(0.1,0.95,Form("%s",Title.Data())); // Bo: this was 
			alice->SetNDC();
			alice->SetTextColor(1);
			alice->SetTextSize(0.062);
			alice->Draw();    
		}
	}
}

void DrawGammaInvMassHisto( TH1* histo1,
							TString Title, TString XTitle, TString YTitle,
							Float_t xMin, Float_t xMax,Int_t bck, Double_t yMin ) {
	
	histo1->GetXaxis()->SetRangeUser(xMin, xMax);
	if (bck != 1 || bck != 2){
		Double_t yMax = 0;
		for (Int_t i = histo1->GetXaxis()->FindBin(xMin); i < histo1->GetXaxis()->FindBin(xMax); i++){
			if (histo1->GetBinContent(i) > yMax){
				yMax = histo1->GetBinContent(i);
			}
		}   
		histo1->GetYaxis()->SetRangeUser(yMin, 1.5*yMax);
	}

	if(XTitle.Length() > 0){
		histo1->SetXTitle(XTitle.Data());
	}
	if(YTitle.Length() > 0){
		histo1->SetYTitle(YTitle.Data());
	}
	histo1->GetYaxis()->SetLabelSize(0.02);
	histo1->GetYaxis()->SetTitleSize(0.025);
	histo1->GetYaxis()->SetDecimals();
	histo1->GetYaxis()->SetTitleOffset(0.5);
	histo1->GetXaxis()->SetTitleSize(0.025);
	histo1->GetXaxis()->SetLabelSize(0.02);
	histo1->SetMarkerStyle(20);
	histo1->SetMarkerColor(1);
	histo1->SetLineColor(1);
	histo1->SetLineWidth(0.5);
	histo1->SetMarkerSize(0.2);
	histo1->SetTitleOffset(1.2,"xy");      
	histo1->SetTitleSize(0.05,"xy");    
	histo1->GetYaxis()->SetLabelSize(0.05);
	histo1->GetXaxis()->SetLabelSize(0.05);
	histo1->GetXaxis()->SetNdivisions(507,kTRUE);
	if(bck==1){
		histo1->SetLineStyle(1);      
		histo1->SetLineColor(4);
		histo1->SetMarkerColor(4);
		histo1->SetMarkerStyle(24);
		histo1->SetLineWidth(0.9);
		histo1->DrawCopy("hist,same");
	}
	if(bck==2){
		histo1->SetLineStyle(1);      
		histo1->SetLineColor(kRed+2);
		histo1->SetMarkerColor(kRed+2);
		histo1->SetLineWidth(0.9);
		histo1->DrawCopy("hist,same");
	}
	if(bck==3){
		histo1->SetLineStyle(1);      
		histo1->SetLineColor(kGreen+2);
		histo1->SetMarkerColor(kGreen+2);
		histo1->SetLineWidth(0.9);
		histo1->DrawCopy("hist,same");
	}
	if(bck==4){
		histo1->SetLineStyle(1);      
		histo1->SetLineColor(807);
		histo1->SetMarkerColor(807);
		histo1->SetLineWidth(0.9);
		histo1->DrawCopy("hist,same");
	}
	if(bck==5){
		histo1->SetLineStyle(1);      
		histo1->SetLineColor(kViolet+2);
		histo1->SetMarkerColor(kViolet+2);
		histo1->SetLineWidth(0.9);
		histo1->DrawCopy("hist,same");
	}

	if(bck==6){
		histo1->SetLineStyle(1);      
		histo1->SetLineColor(kPink+2);
		histo1->SetMarkerColor(kPink+2);
		histo1->SetLineWidth(0.9);
		histo1->DrawCopy("hist,same");
	}
	if(bck==7){
		histo1->SetLineStyle(1);      
		histo1->SetLineColor(kCyan+2);
		histo1->SetMarkerColor(kCyan+2);
		histo1->SetLineWidth(0.9);
		histo1->DrawCopy("hist,same");
	}

	if (bck == 0){
		histo1->SetTitle("");
		histo1->DrawCopy("e1,p");   
		if(Title.Length() > 0){
			histo1->SetTitle("");
			TLatex *alice = new TLatex(0.1,0.95,Form("%s",Title.Data())); // Bo: this was 
			alice->SetNDC();
			alice->SetTextColor(1);
			alice->SetTextSize(0.062);
			alice->Draw();    
		}
	}
}

