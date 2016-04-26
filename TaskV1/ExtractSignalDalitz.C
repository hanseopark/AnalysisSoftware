// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion

//This file is not supposed to be run on outputfiles of the GammaConv-Software before the 30th Sept 2010.

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
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "THnSparse.h"

// Main Function
void ExtractSignalDalitz(TString meson="", TString file="", TString cutSelection="", TString fSuffix="", TString isMC="", TString option="", TString crystal="", TString optionUseMinBiasEff="",TString period="", TString thesis="",Int_t numberOfBins = 30,Bool_t addSig=kFALSE,Int_t mode = 1) {
  
	gROOT->Reset();
	gSystem->Load("libCore.so");
	gSystem->Load("libTree.so");
	gSystem->Load("libGeom.so");
	gSystem->Load("libVMC.so");
	gSystem->Load("libPhysics.so");
	gSystem->Load("libMinuit.so");
	gSystem->Load("libSTEERBase");
	gSystem->Load("libESD");
	gSystem->Load("libAOD");
	gSystem->Load("libANALYSIS");
	gSystem->Load("libANALYSISalice");
	gSystem->Load("libCORRFW.so");
        gSystem->Load("libPWGGAGammaConv.so");
	
	
	TString fEventCutSelection    = "";
	TString fGammaCutSelection    = "";
	TString fClusterCutSelection  = "";
	TString fElectronCutSelection = "";
	TString fMesonCutSelection    = "";
	
	Double_t Pi0DalitzBR = 0.0;
	Double_t Pi0GGBR     = 0.0;
	

	
	fCutSelection = cutSelection;
	TString fCutSelectionRead = cutSelection;
	fMode = mode;

	
	//Int_t mode = 9;   old output Dalitz
	//      mode = 1    new output Dalitz
	
	
	
	
	if ( fMode == 9  ){
	  
		ReturnSeparatedCutNumber(fCutSelection, fGammaCutSelection, fElectronCutSelection,fMesonCutSelection,kTRUE);
		fEventCutSelection = fGammaCutSelection(0,7);
		fGammaCutSelection = fGammaCutSelection(7,fGammaCutSelection.Length()-7);
		cout << fEventCutSelection.Data() << "\t" << fGammaCutSelection.Data() << endl;
		
		
	} else {
	
		ReturnSeparatedCutNumberAdvanced(fCutSelection, fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, fMode);
		cout<<"Testing 0: " << fEventCutSelection.Data() << "\t" << fGammaCutSelection.Data() <<"\t"<<fClusterCutSelection.Data()<<"\t"<< fElectronCutSelection.Data()<<"\t"<<fMesonCutSelection.Data()<<endl;
	}    
	
	
	TString fEventCutSelectionRead    = fEventCutSelection.Data();
	TString fGammaCutSelectionRead    = fGammaCutSelection.Data();
	TString fClusterCutSelectionRead  = fClusterCutSelection.Data();
	TString fElectronCutSelectionRead = fElectronCutSelection.Data();
	TString fMesonCutSelectionRead    = fMesonCutSelection.Data();
	
	
	if ( fMode == 6 || fMode == 7 ){
	  
	  
	    fCutSelectionRead =  Form("%s_%s_%s_%s_%s", fEventCutSelection.Data(), fGammaCutSelection.Data(), fElectronCutSelection.Data(),fClusterCutSelection.Data(), fMesonCutSelection.Data());
	    cout<<"Testing: "<<fCutSelectionRead.Data()<<endl;
	  
	  
	}
	
	
		
	if ( addSig ) {
	  
		cout << "running added Signal" << endl;
		cout << fEventCutSelection.Data() << endl;
		fEventCutSelection.Replace(6,1,"2");
		cout << fEventCutSelection.Data() << endl;
		fEventCutSelectionRead 	  = fEventCutSelection;
		fGammaCutSelectionRead	  = fGammaCutSelection;
		fElectronCutSelectionRead = fElectronCutSelection;
		fMesonCutSelectionRead 	  = fMesonCutSelection;
		if( fMode == 0 )fCutSelectionRead = Form("%s%s_%s_%s",       fEventCutSelection.Data(), fGammaCutSelection.Data(), fElectronCutSelection.Data(), fMesonCutSelection.Data());
		else if( fMode == 1 ) fCutSelectionRead = Form("%s_%s_%s_%s",fEventCutSelection.Data(), fGammaCutSelection.Data(), fElectronCutSelection.Data(), fMesonCutSelection.Data());
		cout << fCutSelectionRead.Data() << endl;
	}
	
	if( optionUseMinBiasEff.CompareTo("MinBiasEffOnly") == 0 && isMC.CompareTo("kTRUE") == 0 ){
	  
	 // if( isMC.CompareTo("kTRUE") == 0 ){
	  
		cout << "calculating MinBias Eff" << endl;
		cout << fEventCutSelection.Data() << endl;
		fEventCutSelection.Replace(1,2,"00");
		fEventCutSelectionRead 		= fEventCutSelection;
		fGammaCutSelectionRead 		= fGammaCutSelection;

		//fElectronCutSelection.Replace(18,1,"0");
		
		cout <<"Calculation "<<fElectronCutSelection.Data()<<endl;
		
		fElectronCutSelectionRead 	= fElectronCutSelection;
		
		fMesonCutSelectionRead 		= fMesonCutSelection;      
		cout << fGammaCutSelection.Data() << endl;
		if (fMode==0)fCutSelectionRead = Form("%s%s_%s_%s", 	     fEventCutSelection.Data(), fGammaCutSelection.Data(), fElectronCutSelection.Data(), fMesonCutSelection.Data());
		else if (fMode==1)fCutSelectionRead = Form("%s_%s_%s_%s",     fEventCutSelection.Data(), fGammaCutSelection.Data(), fElectronCutSelection.Data(), fMesonCutSelection.Data());
		cout << fCutSelectionRead.Data() << endl;
	}
	
	
	
	StyleSettingsThesis();
	SetPlotStyle();

	fEnergyFlag = option;
	fPrefix=meson;
	fPeriodFlag = period;
	fThesis = thesis;

	TString outputDir = Form("%s/%s/%s/ExtractSignalDalitz",fCutSelection.Data(),option.Data(),fSuffix.Data());
	gSystem->Exec("mkdir -p "+outputDir);
	cout << fSuffix.Data()<< endl;
	

	date = ReturnDateString();
	
	//TString	fCutSelectionRead = cutSelection;
	
	

	//****************************** Specification of collision system ************************************************
	
	TString textProcess = ReturnMesonString(fPrefix);
	if(textProcess.CompareTo("") == 0 ){
	    cout << "Meson unknown" << endl;
	    return ;
	}
	
	fTextMeasurement = ReturnFullTextMeson(fEnergyFlag, textProcess,kTRUE);
	fCollisionSystem = ReturnFullCollisionsSystem(fEnergyFlag);
	if (fCollisionSystem.CompareTo("") == 0){
	      cout << "No correct collision system specification, has been given" << endl;
	      return;
	}

	fDetectionProcess = ReturnFullTextReconstructionProcess(1); //Dalitz
	
	
	//****************************** Choice of Fitting procedure ******************************************************
	if(crystal.CompareTo("CrystalBall") == 0){// means we want to plot values for the pi0
		fCrysFitting=1;
		cout << "CrystalBall fit chosen ..." << endl;
	}
	else	{
		fCrysFitting=0;
		cout << "Gaussian fit chosen ..." << endl;
	}

	if(cutSelection.Length() == 0){
		cout<<"ERROR: Cut selection is not set, please do!"<<endl;
		return;
	}

	//***************************** Specification Data/MC ************************************************************
	if(isMC.CompareTo("kTRUE") == 0){
		fIsMC = 1;
		fPrefix2 = "MC";
	} else {
		fIsMC = 0;
		fPrefix2 = "data";
	}

	//***************************** Initialization of variables according to meson type ******************************
	if(meson.CompareTo("Pi0") == 0){
		Initialize("Pi0",numberOfBins);
	} else if (meson.CompareTo("Eta") == 0) {
		Initialize("Eta",numberOfBins);
	} else if(meson.CompareTo("Pi0EtaBinning") == 0) {
		Initialize("Pi0EtaBinning",numberOfBins);
	} else	{
		cout<<"ERROR: First argument in the ExtractSignalDalitz(....) has to be either Pi0 or Eta or Pi0EtaBinning"<<endl;
		return;
	}

	
	
	//************************* Start of Main routine ***************************************************************
	const char* fFileErrLogDatname = Form("%s/%s/%s_%s_FileErrLogDalitz%s_%s.dat",cutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelectionRead.Data());
	fFileErrLog.open(fFileErrLogDatname, ios::out);

	

	TFile f(file.Data());
	
	
		   
	
	
	cout <<"EventSelection: "<<fEventCutSelection.Data()<<", PhotonCut: " << fGammaCutSelection.Data() << ", ElectronCut: " << fElectronCutSelection.Data() << ", MesonCut: " << fMesonCutSelection.Data() << endl;
	
	TString nameOutputDir = "";
	
	if( fMode == 9 || fMode == 1 ){
	  
	     nameOutputDir = "GammaConvDalitzV1";
	}
	
	
	else if( fMode == 6 || fMode == 7 ){
	  
	    nameOutputDir = "GammaConvDalitzCalo";
	  
	}
	
	TList* GammaConvDalitzV1 = (TList*)f.Get(nameOutputDir.Data());

      
	  
	
	
        if( ! GammaConvDalitzV1 ) {

                    Int_t iTrain = 1;

                    while( ! GammaConvDalitzV1 && iTrain <= 100 ) {
                            
                          GammaConvDalitzV1 = (TList*)f.Get(Form("%s_%d",nameOutputDir.Data(),iTrain));
                           
                          iTrain++;
                    
                    }         
                   if( ! GammaConvDalitzV1 ) {
		     
			 GammaConvDalitzV1 = (TList*)f.Get("GammaConvDalitzCalo");
                    }
            
                    if( ! GammaConvDalitzV1 ) {
                        cout<<"ERROR: GammaConvDalitzV1 is not found in the file"<<endl;

                        return;
                    }
            
        }
  	

	TList* fHistosGammaConversionDalitz = (TList*) GammaConvDalitzV1->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()) );
	if( ! fHistosGammaConversionDalitz ){		
		cout<<Form("Cut Number %s",fCutSelectionRead.Data())<<" is not found in the file"<<endl;
		return;
	}
	
	TList *fESDContainer        = (TList*) fHistosGammaConversionDalitz->FindObject( Form("%s ESD histograms",   fCutSelectionRead.Data()));
	TList *fBackgroundContainer = (TList*) fHistosGammaConversionDalitz->FindObject( Form("%s Back histograms",  fCutSelectionRead.Data()));
        TList *fMotherContainer     = (TList*) fHistosGammaConversionDalitz->FindObject( Form("%s Mother histograms",fCutSelectionRead.Data()));
	TList *fQAContainer         = (TList*) fHistosGammaConversionDalitz->FindObject( Form("%s QA histograms",    fCutSelectionRead.Data()));
       
      
	
        fNumberOfGoodESDTracksVtx = (TH1F*)fESDContainer->FindObject("GoodESDTracks");
        fEventQuality = (TH1F*)fESDContainer->FindObject("NEvents");
	
	TString  nameInvMassPi0 = "ESD_DalitzBackground_InvMass_Pt";
	TString  nameInvBackPi0 = "ESD_DalitzMother_InvMass_Pt";
	
	if ( fMode == 6 || fMode == 7 ){
	  
	      nameInvMassPi0 = "ESD_Background_InvMass_Pt";
	      nameInvBackPi0 = "ESD_Mother_InvMass_Pt";
	  
	}
	
	

	TH2D *fBck2D         = (TH2D*)fESDContainer->FindObject(nameInvMassPi0.Data());
	TH2D *fGammaGamma2D  = (TH2D*)fESDContainer->FindObject(nameInvBackPi0.Data());
	TH1D *fBck           = (TH1D*) fBck2D->ProjectionX("ESD_DalitzBackground_InvMass");
	TH1D *fGammaGamma = (TH1D*) fGammaGamma2D->ProjectionX("ESD_DalitzMother_InvMass");

	

	
	TString intermediate = GetCentralityString(fGammaCutSelection);

	if (intermediate.CompareTo("pp")==0){
	  fTextCent = "MinBias";	
	} else {
	  fTextCent = Form("%s central", intermediate.Data());
	}
	
	TString rapidityRange;
	fYMaxMeson =  ReturnRapidityStringAndDouble(fMesonCutSelection, rapidityRange);

	
	TString fDecayChannel = "e^{+}e^{-}#gamma";
	fBackgroundMultCutNumber = fElectronCutSelection(14,1);
	
	if (fBackgroundMultCutNumber.CompareTo("0") == 0){
        fBackgroundMultNumber=5;
        cout << "using number of events for BG 5" << endl;
         }
        if (fBackgroundMultCutNumber.CompareTo("1") == 0){
            fBackgroundMultNumber=10;
            cout << "using number of events for BG 10" << endl;
        }
        if (fBackgroundMultCutNumber.CompareTo("2") == 0){
            fBackgroundMultNumber=15;
            cout << "using number of events for BG 15" << endl;
        }
        if (fBackgroundMultCutNumber.CompareTo("3") == 0){
            fBackgroundMultNumber=20;
            cout << "using number of events for BG 20" << endl;
        }
        if (fBackgroundMultCutNumber.CompareTo("4") == 0){
            fBackgroundMultNumber=2;
            cout << "using number of events for BG 2" << endl;
        }
        if (fBackgroundMultCutNumber.CompareTo("5") == 0){
            fBackgroundMultNumber=50;
            cout << "using number of events for BG 50" << endl;
        }
        if (fBackgroundMultCutNumber.CompareTo("6") == 0){
            fBackgroundMultNumber=80;
            cout << "using number of events for BG 80" << endl;
        }
        if (fBackgroundMultCutNumber.CompareTo("7") == 0){
            fBackgroundMultNumber=100;
            cout << "using number of events for BG 100" << endl;
        }	


	TString fObjectNameESD;
	TString fObjectNameBck;
	TString fObjectNameTrue;
        TString fObjectNameTruePi0Cont;
	TString fObjectNameTrueGGBck;
        TString fObjectNameTrueContBck;
	
		
	fObjectNameESD         = 	"ESD_DalitzMother_InvMass_Pt";
	fObjectNameBck         = 	"ESD_DalitzBackground_InvMass_Pt";
	fObjectNameTrue        = 	"ESD_TrueMother_InvMass_Pt";
        fObjectNameTruePi0Cont =        "ESD_TrueMotherPi0GG_InvMass_Pt";
        fObjectNameTrueGGBck   =        "ESD_TrueDalitzBckGG_InvMass_Pt";
        fObjectNameTrueContBck =        "ESD_TrueDalitzBckCont_InvMass_Pt";
	
	if( fMode == 6 || fMode == 7 ){
	  	  
	    fObjectNameESD         =        "ESD_Mother_InvMass_Pt";
	    fObjectNameBck         =        "ESD_Background_InvMass_Pt";
	    //fObjectNameTrue	   =        "ESD_TruePi0_InvMass_Pt";
	    fObjectNameTrue        =        "ESD_TruePi0NoShower_InvMass_Pt";
	    fObjectNameTruePi0Cont =        "ESD_TruePi0GG_InvMass_Pt";
	    fObjectNameTrueGGBck   =        "ESD_TrueBckGG_InvMass_Pt";
	    fObjectNameTrueContBck =        "ESD_TrueBckCont_InvMass_Pt";
	    
	
	    if ( meson.CompareTo("Eta") == 0 ) {
	
		fObjectNameTrue        =    "ESD_TrueEta_InvMass_Pt";
		fObjectNameTruePi0Cont =    "ESD_TrueEtaGG_InvMass_Pt";
        
	    }
	    
	  	  
	}
	
	
	
	TH2F *fGammaGammaInvMassVSPt= (TH2F*)fESDContainer->FindObject(fObjectNameESD.Data());
	fGammaGammaInvMassVSPt->Sumw2();
	
	TH2F *fBckInvMassVSPt = (TH2F*)fESDContainer->FindObject(fObjectNameBck.Data());
	fBckInvMassVSPt->Sumw2();
	
	fESDEposEnegInvMassPt = 0x0;
	
	if ( fQAContainer ) {
	
	     fESDEposEnegInvMassPt = (TH2F*)fQAContainer->FindObject("ESD_EposEneg_InvMassPt");
	     
	}
	
	
	
	
	const char* fFileDataLogname = Form("%s/%s/%s_%s_EffiCheck_RAWDATA_Dalitz%s_%s.dat",cutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelectionRead.Data());
	fFileDataLog.open(fFileDataLogname, ios::out);

   
	cout<<"LLego "<<endl; 	

        ProduceBckProperWeighting(fBackgroundContainer,fMotherContainer,fESDContainer);
	
	cout<<"Salio "<<endl; 	
	


	if(fIsMC){

		TString namePi0InAcc    = "MC_Pi0DalitzInAcc_Pt";
		TString namePi0Pt       = "MC_Pi0_Pt";
		TString namePi0DalitzPt = "MC_Pi0_Dalitz_Pt";
		TString namePi0GGPt     = "MC_Pi0_GG_Pt";
		
		TString nameEtaInAcc    = "MC_EtaDalitzInAcc_Pt";
		TString nameEtaPt       = "MC_Eta_Pt";
		TString nameEtaDalitzPt = "MC_Eta_Dalitz_Pt";
		TString nameEtaGGPt     = "MC_Eta_GG_Pt";
		
		
		if( fMode == 6 || fMode == 7 ){
		  
		     namePi0InAcc    = "MC_Pi0InAcc_Pt";
		     namePi0Pt       = "MC_Pi0_Pt";
		     namePi0DalitzPt = "MC_Pi0_Dalitz_Pt";
		     namePi0GGPt     = "MC_Pi0_GG_Pt";
		     
		     
		     nameEtaInAcc    = "MC_EtaInAcc_Pt";
		     nameEtaPt       = "MC_Eta_Pt";
		     nameEtaDalitzPt = "MC_Eta_Dalitz_Pt";
		     nameEtaGGPt     = "MC_Eta_GG_Pt";
		  
		}
	  
	  
		TList *fTrueHistograms = (TList*)fHistosGammaConversionDalitz->FindObject(Form("%s True histograms",fCutSelectionRead.Data()) );
		if( !fTrueHistograms ) {
			cout<<Form("%s  True histograms",fCutSelectionRead.Data())<<" is not found in the file"<<endl;
			return;
		}

		TList *fMChistograms = (TList*)fHistosGammaConversionDalitz->FindObject(Form("%s MC histograms",fCutSelectionRead.Data()) );
		if( !fMChistograms ) {
			cout<<Form("%s  MC histograms",fCutSelectionRead.Data())<<" is not found in the file"<<endl;
			return;
		}

		if( fMesonId == 111){
			fHistoMCMesonPtWithinAcceptance = (TH1D*)fMChistograms->FindObject(namePi0InAcc.Data());
			fHistoMCMesonPt			= (TH1D*)fMChistograms->FindObject(namePi0Pt.Data());   // Not the best; better having a 2D Pt_vs_Rapid in case we change limits
			fHistoMCMesonDalitzPt           = (TH1D*) fHistoMCMesonPt->Clone(namePi0DalitzPt.Data());
			fHistoMCMesonGGPt               = (TH1D*) fMChistograms->FindObject(namePi0GGPt.Data());
		}
		
		if( fMesonId == 221){
			fHistoMCMesonPtWithinAcceptance = (TH1D*)fMChistograms->FindObject(nameEtaInAcc.Data());
			fHistoMCMesonPt 		= (TH1D*)fMChistograms->FindObject(nameEtaPt.Data());   // Not the best; better having a 2D Pt_vs_Rapid in case we change limits
			fHistoMCMesonDalitzPt           = (TH1D*)fHistoMCMesonPt->Clone(nameEtaDalitzPt.Data());
			fHistoMCMesonGGPt               = (TH1D*)fMChistograms->FindObject(nameEtaGGPt.Data());
			
		}
		
		
		
		
		
		if( meson.CompareTo("Pi0") == 0 || meson.CompareTo("Pi0EtaBinning") == 0 ){
		  
		  
		  /////////////////////////Saving Branching Ratio//////////////////////////////
	
		  cout<<"Saving Branching Ration from DPG"<<endl;
	
		  fArrayBRPi0Meson->AddAt( fPi0GGBRDPG,	  0);  //  0.98823 DPG
		  fArrayBRPi0Meson->AddAt( fPi0DalitzBRDPG, 1);  //  0.01174 DPG
		  cout<<"Pi0->gg: "<<fArrayBRPi0Meson->GetAt(0)<< "  Pi0->eeg "<<fArrayBRPi0Meson->GetAt(1)<<endl;
	

		  ////////////////////////////////////////////////////////////////////////////
			   
		  cout<<"Computing BR from MC histos"<<endl;
		  Pi0DalitzBR =  fHistoMCMesonDalitzPt->GetEntries() / ( fHistoMCMesonDalitzPt->GetEntries() + fHistoMCMesonGGPt->GetEntries() );
		  Pi0GGBR     =  fHistoMCMesonGGPt->GetEntries()     / ( fHistoMCMesonDalitzPt->GetEntries() + fHistoMCMesonGGPt->GetEntries() );
			  
	          fArrayBRPi0Meson->AddAt( Pi0GGBR,	2);  // 
	          fArrayBRPi0Meson->AddAt( Pi0DalitzBR, 3);  // 
		
		  cout<<"The MC Branching ratios are set as follow:"<<endl;
		  
		  cout<<"Pi0->gg: "<<fArrayBRPi0Meson->GetAt(2)<<"  Pi0->e+e-G: "<<fArrayBRPi0Meson->GetAt(3)<<endl;
		    
			     
		} 
		
		
		
		
		fHistoTrueMesonInvMassVSPt = (TH2F*)fTrueHistograms->FindObject(fObjectNameTrue.Data());
		FillMassMCTrueMesonHistosArray(fHistoTrueMesonInvMassVSPt);
		
		fHistoTrueGGMesonInvMassVSPt = (TH2F*)fTrueHistograms->FindObject(fObjectNameTruePi0Cont.Data());
		FillMassMCTrueGGMesonHistosArray(fHistoTrueGGMesonInvMassVSPt);

                //fHistoTrueGGBckInvMassVSPt = (TH2F*)fTrueHistograms->FindObject(fObjectNameTrueGGBck.Data());
                //FillMassMCTrueMesonGGBckHistosArray( fHistoTrueGGBckInvMassVSPt );
            
                //fHistoTrueContBckInvMassVSPt = (TH2F*)fTrueHistograms->FindObject(fObjectNameTrueContBck.Data());
                //FillMassMCTrueMesonBckContHistosArray(  fHistoTrueContBckInvMassVSPt );

	}

	fMesonMassExpect = TDatabasePDG::Instance()->GetParticle(fMesonId)->Mass();
	
	
	
	 if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 || fEnergyFlag.CompareTo("pPb_5.023TeV") == 0){ 
	
	   fNEvents = fEventQuality->GetBinContent(1);
	   
	 } else {
	   
	      fNEvents =  GetNEvents(fEventQuality);
	 }
	

   
	cout<< "The mass of the meson is: "<< fMesonMassExpect<< " Events analysed:  "<< fNEvents<< endl;

	// Process the 1D invariant mass histos

	fGammaGamma->SetTitle(Form("%s %s",fGammaGamma->GetTitle(),fCutSelection.Data()));
	fGammaGamma->Rebin(fNRebinGlobal);

	//fGammaGamma->Scale(1./fNRebinGlobal);
	fBck->Rebin(fNRebinGlobal);
	//fBck->Scale(1./fNRebinGlobal);
	ProcessEM( fGammaGamma , fBck, fBGFitRange);
	fHistoMappingBackNormInvMass = fBckNorm;
	fHistoMappingSignalInvMass = fSignal;

	fGammaGamma->DrawCopy();
	fHistoMappingBackNormInvMass->DrawCopy("same");
	fHistoMappingSignalInvMass->DrawCopy("same");

	
	//Function to project the 2D histos of InvariantMass ePluseMinus to Invariant Mass spectrum
	
	if ( fESDEposEnegInvMassPt ){
	  
	    FillEposEnegHistosArray( fESDEposEnegInvMassPt );
	    
	}
	
	
	// Function to Project the 2D histos InvariantMass VS Pt into Invariant Mass spectrum
	FillMassHistosArray(fGammaGammaInvMassVSPt);


	ProcessEM( fMesonFullPtSignal, fMesonFullPtBackground, fBGFitRange);
	fMesonFullPtBackNorm = fBckNorm;

	ProcessEM( fFittingHistMidPtSignal, fFittingHistMidPtBackground, fBGFitRange);
	fFittingHistMidPtSignalSub = fSignal;
	if(fCrysFitting==0){
		fFileErrLog << "Using exp fit"<<endl;
		FitSubtractedInvMassInPtBins(fFittingHistMidPtSignalSub, fMesonIntDeltaRange,200,kTRUE);
		fFitSignalInvMassMidPt = fFitReco;
	} else {
		fFileErrLog << "Using Crystal Ball function"<<endl;
		FitCBSubtractedInvMassInPtBins(fFittingHistMidPtSignalSub, fMesonIntDeltaRange,200,kTRUE,"SinglefitfunctionMidPt");
		fFitSignalInvMassMidPt = fFitReco;
	}

	TString nameMesonFittingMidPt= Form("%s/%s_%s_MesonSubtractedFittingMidPt%s_%s_%02d.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(), fCutSelection.Data(), 200, fSuffix.Data());
	TString nameCanvasFittingMidPt= "MesonCanvasSubtractedFittingMidPt";
	TString namePadFittingMidPt= "MesonPadSubtractedFittingMidPt";
	//PlotWithFitSubtractedInvMassSinglePtBin2( fFittingHistMidPtSignalSub, fFitSignalInvMassMidPt, nameMesonFittingMidPt, nameCanvasFittingMidPt, fMesonMassRange, fIsMC,fDecayChannel );

	delete fMidPt;

	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){ // BEGIN ANALYSIS for each Pt bin
		cout << "Begin Analysis Pt Bin " << iPt <<endl;
		// Function to subtract GG minus Bck

		ProcessEM( fHistoMappingGGInvMassPtBin[iPt], fHistoMappingBackInvMassPtBin[iPt], fBGFitRange);
		fHistoMappingSignalInvMassPtBin[iPt] = fSignal;
		fHistoMappingBackNormInvMassPtBin[iPt] = fBckNorm;

		// Integrate the 2g histo

		cout<< "iPt"<< iPt<< " "<< "standard range"<<endl;
		// Fitting the subtracted spectra
		fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "normal range/right normalization" << endl;

		fFitSignalInvMassPtBin[iPt]=0x00;

		if(fCrysFitting==0){
			fFileErrLog << "Using exp fit"<<endl;
			FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE);
			fFitSignalInvMassPtBin[iPt] = fFitReco;
			fFitSignalPeakPosInvMassPtBin[iPt] = fFitGausExp;
			fFitBckInvMassPtBin[iPt] = fFitLinearBck;
			fMesonYieldsResidualBckFunc[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncError[iPt] = fIntLinearBckError;

		} else {
			fFileErrLog << "Using Crystal Ball function"<<endl;
			FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE,Form("CBFitFuncNormalBin%02d",iPt));
			fFitSignalInvMassPtBin[iPt] = fFitReco;
			fFitSignalPeakPosInvMassPtBin[iPt] = fFitGausExp;
			fFitBckInvMassPtBin[iPt] = fFitLinearBck;
			fMesonYieldsResidualBckFunc[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncError[iPt] = fIntLinearBckError;
		}
		
		if (fFitSignalInvMassPtBin[iPt] !=0x00){
		  
			fMesonMass[iPt] = fFitSignalInvMassPtBin[iPt]->GetParameter(1);
			fMesonMassError[iPt] = fFitSignalInvMassPtBin[iPt]->GetParError(1);
			
			fMesonWidth[iPt] = fFitSignalInvMassPtBin[iPt]->GetParameter(2);
			fMesonWidthError[iPt] = fFitSignalInvMassPtBin[iPt]->GetParError(2);
			
			fMesonCurIntRange[0] 			= fMesonMass[iPt] + fMesonIntDeltaRange[0];
			fMesonCurIntRangeWide[0] 		= fMesonMass[iPt] + fMesonIntDeltaRangeWide[0];
			fMesonCurIntRangeNarrow[0]		= fMesonMass[iPt] + fMesonIntDeltaRangeNarrow[0];
			fMesonCurIntRange[1] 			= fMesonMass[iPt] + fMesonIntDeltaRange[1];
			fMesonCurIntRangeWide[1] 		= fMesonMass[iPt] + fMesonIntDeltaRangeWide[1];
			fMesonCurIntRangeNarrow[1] 		= fMesonMass[iPt] + fMesonIntDeltaRangeNarrow[1];
			
			/*
			fMesonCurIntRange[0] = fMesonIntDeltaRange[0] - (fMesonMassExpect-fMesonMass[iPt]);
			fMesonCurIntRangeWide[0] = fMesonIntDeltaRangeWide[0] - (fMesonMassExpect-fMesonMass[iPt]);
			fMesonCurIntRangeNarrow[0] = fMesonIntDeltaRangeNarrow[0] - (fMesonMassExpect-fMesonMass[iPt]);
			fMesonCurIntRange[1] = fMesonIntDeltaRange[1] - (fMesonMassExpect-fMesonMass[iPt]);
			fMesonCurIntRangeWide[1] = fMesonIntDeltaRangeWide[1] - (fMesonMassExpect-fMesonMass[iPt]);
			fMesonCurIntRangeNarrow[1] = fMesonIntDeltaRangeNarrow[1] - (fMesonMassExpect-fMesonMass[iPt]);
			*/
			
			
		} else {
		        fMesonMass[iPt] 				= fMesonMassExpect;
			fMesonMassError[iPt] 			= 0.;
			fMesonWidth[iPt] = 0.;
			fMesonWidthError[iPt] = 0.;
			fMesonCurIntRange[0] 			= fMesonMassExpect + fMesonIntDeltaRange[0];
			fMesonCurIntRangeWide[0] 		= fMesonMassExpect + fMesonIntDeltaRangeWide[0];
			fMesonCurIntRangeNarrow[0] 		= fMesonMassExpect + fMesonIntDeltaRangeNarrow[0];
			fMesonCurIntRange[1] 			= fMesonMassExpect + fMesonIntDeltaRange[1];
			fMesonCurIntRangeWide[1] 		= fMesonMassExpect + fMesonIntDeltaRangeWide[1];
			fMesonCurIntRangeNarrow[1] 		= fMesonMassExpect + fMesonIntDeltaRangeNarrow[1];	   
		   
		        /*
			fMesonMass[iPt] = 0.;
			fMesonMassError[iPt] = 0.;
			fMesonWidth[iPt] = 0.;
			fMesonWidthError[iPt] = 0.;
			fMesonCurIntRange[0] = fMesonIntDeltaRange[0];
			fMesonCurIntRangeWide[0] = fMesonIntDeltaRangeWide[0];
			fMesonCurIntRangeNarrow[0] = fMesonIntDeltaRangeNarrow[0];
			fMesonCurIntRange[1] = fMesonIntDeltaRange[1];
			fMesonCurIntRangeWide[1] = fMesonIntDeltaRangeWide[1];
			fMesonCurIntRangeNarrow[1] = fMesonIntDeltaRangeNarrow[1];*/
			
		}
		
		IntegrateHistoInvMass( fHistoMappingGGInvMassPtBin[iPt], fMesonCurIntRange);
		fGGYields[iPt] = fYields;
		fGGYieldsError[iPt] = fYieldsError;
		IntegrateHistoInvMass( fHistoMappingGGInvMassPtBin[iPt], fMesonCurIntRangeWide);
		fGGYieldsWide[iPt] = fYields;
		fGGYieldsWideError[iPt] = fYieldsError;
		IntegrateHistoInvMass( fHistoMappingGGInvMassPtBin[iPt], fMesonCurIntRangeNarrow);
		fGGYieldsNarrow[iPt] = fYields;
		fGGYieldsNarrowError[iPt] = fYieldsError;

		// Integrate the bck histo
		IntegrateHistoInvMass( fHistoMappingBackNormInvMassPtBin[iPt], fMesonCurIntRange);
		fBckYields[iPt] = fYields;
		fBckYieldsError[iPt] = fYieldsError;
		IntegrateHistoInvMass( fHistoMappingBackNormInvMassPtBin[iPt], fMesonCurIntRangeWide);
		fBckYieldsWide[iPt] = fYields;
		fBckYieldsWideError[iPt] = fYieldsError;
		IntegrateHistoInvMass( fHistoMappingBackNormInvMassPtBin[iPt], fMesonCurIntRangeNarrow);
		fBckYieldsNarrow[iPt] = fYields;
		fBckYieldsNarrowError[iPt] = fYieldsError;

		// Integrate the signal histo
		fFileDataLog<< endl <<"Signal histo normal range, right norm "<< fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
		IntegrateHistoInvMassStream( fHistoMappingSignalInvMassPtBin[iPt], fMesonCurIntRange);
		fMesonYields[iPt] = fYields;
		fMesonYieldsError[iPt] = fYieldsError;
		fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
		fFileDataLog<< endl <<"Signal histo wide range, right norm " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
		IntegrateHistoInvMassStream( fHistoMappingSignalInvMassPtBin[iPt], fMesonCurIntRangeWide);
		fMesonYieldsWide[iPt] = fYields;
		fMesonYieldsWideError[iPt] = fYieldsError;
		fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
		fFileDataLog<< endl <<"Signal histo narrow range, right norm" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
		IntegrateHistoInvMassStream( fHistoMappingSignalInvMassPtBin[iPt], fMesonCurIntRangeNarrow);
		fMesonYieldsNarrow[iPt] = fYields;
		fMesonYieldsNarrowError[iPt] = fYieldsError;
		fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
		if(fIsMC){
			fFileDataLog<< endl <<"True histo normal range" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			fFitTrueSignalInvMassPtBin[iPt]=0x00;
			if(fCrysFitting==0){
				fFileErrLog << "Using exp fit"<<endl;
				FitTrueInvMassInPtBins(fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonIntDeltaRange,iPt);
			} else {
				fFileErrLog << "Using Crystal Ball function"<<endl;
				FitCBSubtractedInvMassInPtBins(fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonIntDeltaRange,iPt,kFALSE,Form("CBFitFuncMCTrueBin%02d",iPt));
			}

			//	FitSubtractedInvMassInPtBins(fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonIntDeltaRange,iPt,kFALSE);
			if (fHistoMappingTrueMesonInvMassPtBins[iPt]->GetEntries() !=0){
				fFitTrueSignalInvMassPtBin[iPt] = fFitReco;
				if (fFitTrueSignalInvMassPtBin[iPt] != 0x00){
					fMesonTrueMass[iPt] = fFitTrueSignalInvMassPtBin[iPt]->GetParameter(1);
					fMesonTrueMassError[iPt] = fFitTrueSignalInvMassPtBin[iPt]->GetParError(1);
					CalculateFWHM(fFitTrueSignalInvMassPtBin[iPt]);
					fMesonTrueFWHM[iPt] = fFWHMFunc;
					fMesonTrueFWHMError[iPt] = fFWHMFuncError;
					fFileDataLog << "TrueFWHM \t" << fMesonTrueFWHM[iPt] << "\t +-" << fMesonTrueFWHMError[iPt] << endl;
					
					fMesonTrueIntRange[0] 		= fMesonTrueMass[iPt] + fMesonIntDeltaRange[0];
					fMesonTrueIntRangeWide[0] 	= fMesonTrueMass[iPt] + fMesonIntDeltaRangeWide[0];
					fMesonTrueIntRangeNarrow[0] 	= fMesonTrueMass[iPt] + fMesonIntDeltaRangeNarrow[0];
					fMesonTrueIntRange[1] 		= fMesonTrueMass[iPt] + fMesonIntDeltaRange[1] ;
					fMesonTrueIntRangeWide[1] 	= fMesonTrueMass[iPt] + fMesonIntDeltaRangeWide[1];
					fMesonTrueIntRangeNarrow[1] 	= fMesonTrueMass[iPt] + fMesonIntDeltaRangeNarrow[1];
					
					/*
					fMesonTrueIntRange[0] = fMesonIntDeltaRange[0] - (fMesonMassExpect-fMesonTrueMass[iPt]);
					fMesonTrueIntRangeWide[0] = fMesonIntDeltaRangeWide[0] - (fMesonMassExpect-fMesonTrueMass[iPt]);
					fMesonTrueIntRangeNarrow[0] = fMesonIntDeltaRangeNarrow[0] - (fMesonMassExpect-fMesonTrueMass[iPt]);
					fMesonTrueIntRange[1] = fMesonIntDeltaRange[1] - (fMesonMassExpect-fMesonTrueMass[iPt]);
					fMesonTrueIntRangeWide[1] = fMesonIntDeltaRangeWide[1] - (fMesonMassExpect-fMesonTrueMass[iPt]);
					fMesonTrueIntRangeNarrow[1] = fMesonIntDeltaRangeNarrow[1] - (fMesonMassExpect-fMesonTrueMass[iPt]);
					*/

				} else {
				        /*
					fMesonTrueMass[iPt] = 0.;
					fMesonTrueMassError[iPt] = 1.;
					fMesonTrueFWHM[iPt] = 0.;
					fMesonTrueFWHMError[iPt] = 0.;
					fMesonTrueIntRange[0] = fMesonIntDeltaRange[0];
					fMesonTrueIntRangeWide[0] = fMesonIntDeltaRangeWide[0];
					fMesonTrueIntRangeNarrow[0] = fMesonIntDeltaRangeNarrow[0];
					fMesonTrueIntRange[1] = fMesonIntDeltaRange[1];
					fMesonTrueIntRangeWide[1] = fMesonIntDeltaRangeWide[1];
					fMesonTrueIntRangeNarrow[1] = fMesonIntDeltaRangeNarrow[1];
					*/
					fMesonTrueMass[iPt] 		= 0.;
					fMesonTrueMassError[iPt] 	= 1.;
					fMesonTrueFWHM[iPt] 		= 0.;
					fMesonTrueFWHMError[iPt] 	= 0.;
					fMesonTrueIntRange[0]	 	= fMesonMassExpect + fMesonIntDeltaRange[0];
					fMesonTrueIntRangeWide[0] 	= fMesonMassExpect + fMesonIntDeltaRangeWide[0];
					fMesonTrueIntRangeNarrow[0] = fMesonMassExpect + fMesonIntDeltaRangeNarrow[0];
					fMesonTrueIntRange[1] 		= fMesonMassExpect + fMesonIntDeltaRange[1];
					fMesonTrueIntRangeWide[1] 	= fMesonMassExpect + fMesonIntDeltaRangeWide[1];
					fMesonTrueIntRangeNarrow[1] = fMesonMassExpect + fMesonIntDeltaRangeNarrow[1];

				}

			}

			IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonTrueIntRange);
			fMesonTrueYields[iPt] = fYields;
			fMesonTrueYieldsError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

			fFileDataLog<< endl <<"True histo wide range" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonTrueIntRangeWide);
			fMesonTrueYieldsWide[iPt] = fYields;
			fMesonTrueYieldsWideError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

			fFileDataLog<< endl <<"True histo narrow range" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
			IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonTrueIntRangeNarrow);
			fMesonTrueYieldsNarrow[iPt] = fYields;
			fMesonTrueYieldsNarrowError[iPt] = fYieldsError;
			fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;


                        fFileDataLog << endl <<" TrueGG histo "<< fBinsPt[iPt] <<"-"<<fBinsPt[iPt+1] <<endl;
                        IntegrateHistoInvMassStream( fHistoMappingTrueGGMesonInvMassPtBins[iPt], fMesonTrueIntRange);
                        fMesonTrueGGYields[iPt] = fYields;
                        fMesonTrueGGYieldsError[iPt] =  fYieldsError;


                        fFileDataLog<< endl <<" TrueGG histo wide range "<< fBinsPt[iPt] <<"-"<<fBinsPt[iPt+1] <<endl;
                        IntegrateHistoInvMassStream( fHistoMappingTrueGGMesonInvMassPtBins[iPt], fMesonTrueIntRangeWide);
                        fMesonTrueGGYieldsWide[iPt] = fYields;
                        fMesonTrueGGYieldsWideError[iPt] = fYieldsError;


                        fFileDataLog<< endl <<" TrueGG histo narrow range "<< fBinsPt[iPt] <<"-"<<fBinsPt[iPt+1] <<endl;
                        IntegrateHistoInvMassStream( fHistoMappingTrueGGMesonInvMassPtBins[iPt], fMesonTrueIntRangeNarrow);
                        fMesonTrueGGYieldsNarrow[iPt] = fYields;
                        fMesonTrueGGYieldsNarrowError[iPt] = fYieldsError;

			
			if( ( fGGYields[iPt] - fMesonTrueYields[iPt]) > 0) {
				fMesonTrueSB[iPt]   = fMesonTrueYields[iPt] / ( fGGYields[iPt] - fMesonTrueYields[iPt] );
				fMesonTrueSign[iPt] = fMesonTrueYields[iPt] / pow( ( fGGYields[iPt] - fMesonTrueYields[iPt] ) , 0.5);
				fMesonTrueSBError[iPt] = 0;
				fMesonTrueSignError[iPt] = 0;
			}
			else {
				fMesonTrueSB[iPt] = 0.;
				fMesonTrueSign[iPt] = 0.;
				fMesonTrueSBError[iPt] = 0.;
				fMesonTrueSignError[iPt] = 0.;
			}
		}
		fFileDataLog<< "Residual Background leftover norm integration/right norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFunc[iPt] << "\t +- \t" << fMesonYieldsResidualBckFuncError[iPt] << endl<< endl;
		fTotalBckYields[iPt] = fBckYields[iPt] + fMesonYieldsResidualBckFunc[iPt];
		fFileDataLog<<"iPt "<<iPt<<" fBckYields "<< fBckYields[iPt] <<" fMesonYieldsResidualBckFunc "<<fMesonYieldsResidualBckFunc[iPt]<<" fTotalBckYields "<<fTotalBckYields[iPt]<<endl;
		
		
		fTotalBckYieldsError[iPt] = pow(fBckYieldsError[iPt]*fBckYieldsError[iPt] + fMesonYieldsResidualBckFuncError[iPt]*fMesonYieldsResidualBckFuncError[iPt],0.5);

		fMesonYieldsCorResidualBckFunc[iPt] = fMesonYields[iPt]- fMesonYieldsResidualBckFunc[iPt];
		
		fFileDataLog<<"iPt "<<iPt<<" fMesonYields "<< fMesonYields[iPt] <<" fMesonYieldsResidualBckFunc "<<fMesonYieldsResidualBckFunc[iPt]<<" fMesonYieldsCorResidualBckFunc "<<fMesonYieldsCorResidualBckFunc[iPt]<<" fMesonYieldsResidualBckFunc/fMesonYields "<<fMesonYieldsResidualBckFunc[iPt]/fMesonYields[iPt]<<endl;
		
		
		fMesonYieldsCorResidualBckFuncError[iPt] =
		pow((fMesonYieldsError[iPt]*fMesonYieldsError[iPt]+
		fMesonYieldsResidualBckFuncError[iPt]*fMesonYieldsResidualBckFuncError[iPt]),0.5);
		fMesonYieldsPerEvent[iPt]= fMesonYieldsCorResidualBckFunc[iPt]/fNEvents;
		fMesonYieldsPerEventError[iPt]= fMesonYieldsCorResidualBckFuncError[iPt]/fNEvents;

		//Integrate Fit Function
		IntegrateFitFunc( fFitSignalPeakPosInvMassPtBin[iPt], fHistoMappingSignalInvMassPtBin[iPt], fMesonCurIntRange);
		fMesonYieldsFunc[iPt]=fYieldsFunc;

		//GetFWHM
		CalculateFWHM( fFitSignalInvMassPtBin[iPt]);
		fMesonFWHM[iPt] = fFWHMFunc;
		fMesonFWHMError[iPt] = fFWHMFuncError;
		
		Double_t MesonCurIntRange[2];
		MesonCurIntRange[0] = fMesonMass[iPt] - fMesonFWHM[iPt];
		MesonCurIntRange[1] = fMesonMass[iPt] + fFitSignalInvMassPtBin[iPt]->GetParameter(2);

		
		if( fTotalBckYields[iPt] > 0 ){

			fMesonSB[iPt]        = fMesonYieldsCorResidualBckFunc[iPt]/fTotalBckYields[iPt];
			fMesonSBError[iPt]   = pow(pow(fMesonYieldsCorResidualBckFuncError[iPt]/fTotalBckYields[iPt],2.)+pow(fMesonYieldsCorResidualBckFunc[iPt]/(fTotalBckYields[iPt]*fTotalBckYields[iPt])*fTotalBckYieldsError[iPt],2.) ,0.5);
			fMesonSign[iPt]      = fMesonYieldsCorResidualBckFunc[iPt]/pow(fTotalBckYields[iPt],0.5);
			fMesonSignError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncError[iPt]/pow(fTotalBckYields[iPt],0.5),2.)+pow(0.5*fMesonYieldsCorResidualBckFunc[iPt]/pow(fTotalBckYields[iPt],1.5)*fTotalBckYieldsError[iPt],2.) ,0.5);
			fFileDataLog <<"Pt "<<iPt <<" signal "<<fMesonYieldsCorResidualBckFunc[iPt]<<" background "<<fTotalBckYields[iPt]<<" fMesonSB[iPt] "<<fMesonSB[iPt]  << " fMesonSign[iPt] "<<fMesonSign[iPt]<<endl;
			
			
		}else{
			
                        fMesonSB[iPt]        = 0.;
			fMesonSBError[iPt]   = 0.;
			fMesonSign[iPt]      = 0.;
			fMesonSignError[iPt] = 0.;

		}
		
		cout<< "iPt"<< iPt<< " "<< "FWHM done"<<endl;

		// Wide integration mass window
		cout<< "iPt"<< iPt<< " "<< "wide range"<<endl;
		fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "wide range/right normalization" << endl;

		if(fCrysFitting==0){
			fFileErrLog << "Using exp fit"<<endl;
			FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntDeltaRangeWide,iPt,kFALSE);
			fMesonYieldsResidualBckFuncWide[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncWideError[iPt] = fIntLinearBckError;
		} else {
			fFileErrLog << "Using Crystal Ball function"<<endl;
			FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntDeltaRangeWide,iPt,kFALSE,Form("CBFitFuncNormalWideBin%02d",iPt));
			fMesonYieldsResidualBckFuncWide[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncWideError[iPt] = fIntLinearBckError;

		}

		//FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt],fMesonIntDeltaRangeWide,iPt,kFALSE);
		fFileDataLog<< "Residual Background leftover wide integration/right norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFuncWide[iPt] <<"\t +- \t" << fMesonYieldsResidualBckFuncWideError[iPt] <<endl<< endl;
		fMesonYieldsCorResidualBckFuncWide[iPt] = fMesonYieldsWide[iPt]- fMesonYieldsResidualBckFuncWide[iPt];

		fTotalBckYieldsWide[iPt] = fBckYieldsWide[iPt] + fMesonYieldsResidualBckFuncWide[iPt];
		fTotalBckYieldsWideError[iPt] = pow(fBckYieldsWideError[iPt]*fBckYieldsWideError[iPt] + fMesonYieldsResidualBckFuncWideError[iPt]*fMesonYieldsResidualBckFuncWideError[iPt],0.5);

		fMesonYieldsCorResidualBckFuncWideError[iPt] =
		pow((fMesonYieldsWideError[iPt]*fMesonYieldsWideError[iPt]+
		fMesonYieldsResidualBckFuncWideError[iPt]*fMesonYieldsResidualBckFuncWideError[iPt]),0.5);
		fMesonYieldsPerEventWide[iPt]= fMesonYieldsCorResidualBckFuncWide[iPt]/fNEvents;
		fMesonYieldsPerEventWideError[iPt]= fMesonYieldsCorResidualBckFuncWideError[iPt]/fNEvents;

		if( fTotalBckYieldsWide[iPt]>0){
			fMesonSBWide[iPt] = fMesonYieldsCorResidualBckFuncWide[iPt]/fTotalBckYieldsWide[iPt];
			fMesonSBWideError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncWideError[iPt]/fTotalBckYieldsWide[iPt],2.)+pow(fMesonYieldsCorResidualBckFuncWide[iPt]/(fTotalBckYieldsWide[iPt]*fTotalBckYieldsWide[iPt])*fTotalBckYieldsWideError[iPt],2.) ,0.5);
			fMesonSignWide[iPt] = fMesonYieldsCorResidualBckFuncWide[iPt]/pow(fTotalBckYieldsWide[iPt],0.5);
			fMesonSignWideError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncWideError[iPt]/pow(fTotalBckYieldsWide[iPt],0.5),2.)+pow(0.5*fMesonYieldsCorResidualBckFuncWide[iPt]/pow(fTotalBckYieldsWide[iPt],1.5)*fTotalBckYieldsWideError[iPt],2.) ,0.5);

		}else{
			fMesonSBWide[iPt] = 0.;
			fMesonSBWideError[iPt] = 0.;
			fMesonSignWide[iPt] = 0.;
			fMesonSignWideError[iPt] = 0.;
		}


		// Narrow integration mass window
		fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "narrow range/right normalization" << endl; ;
		cout<< "iPt"<< iPt<< " "<< "narrow range"<<endl;

		if(fCrysFitting==0){
			fFileErrLog << "Using exp fit"<<endl;
			FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntDeltaRangeNarrow,iPt,kFALSE);
			fMesonYieldsResidualBckFuncNarrow[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncNarrowError[iPt] = fIntLinearBckError;

		} else {
			fFileErrLog << "Using Crystal Ball function"<<endl;
			FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntDeltaRangeNarrow,iPt,kFALSE, Form("CBFitFuncNormalNarrowBin%02d",iPt));
			fMesonYieldsResidualBckFuncNarrow[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncNarrowError[iPt] = fIntLinearBckError;

		}

		//FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt],fMesonIntDeltaRangeNarrow,iPt,kFALSE);
		fFileDataLog<< "Residual Background leftover narrow integration/right norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFuncNarrow[iPt] <<"\t +- \t" << fMesonYieldsResidualBckFuncNarrowError[iPt]<<endl<< endl;
		fMesonYieldsCorResidualBckFuncNarrow[iPt] = fMesonYieldsNarrow[iPt]- fMesonYieldsResidualBckFuncNarrow[iPt];
		fMesonYieldsCorResidualBckFuncNarrowError[iPt] =
		pow((fMesonYieldsNarrowError[iPt]*fMesonYieldsNarrowError[iPt]+
		fMesonYieldsResidualBckFuncNarrowError[iPt]*fMesonYieldsResidualBckFuncNarrowError[iPt]),0.5);
		fMesonYieldsPerEventNarrow[iPt]= fMesonYieldsCorResidualBckFuncNarrow[iPt]/fNEvents;
		fMesonYieldsPerEventNarrowError[iPt]= fMesonYieldsCorResidualBckFuncNarrowError[iPt]/fNEvents;

		fTotalBckYieldsNarrow[iPt] = fBckYieldsNarrow[iPt] + fMesonYieldsResidualBckFuncNarrow[iPt];
		fTotalBckYieldsNarrowError[iPt] = pow(fBckYieldsNarrowError[iPt]*fBckYieldsNarrowError[iPt] + fMesonYieldsResidualBckFuncNarrowError[iPt]*fMesonYieldsResidualBckFuncNarrowError[iPt],0.5);


		if( fTotalBckYieldsNarrow[iPt]>0){
			fMesonSBNarrow[iPt] = fMesonYieldsCorResidualBckFuncNarrow[iPt]/fTotalBckYieldsNarrow[iPt];
			fMesonSBNarrowError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncNarrowError[iPt]/fTotalBckYieldsNarrow[iPt],2.)+pow(fMesonYieldsCorResidualBckFuncNarrow[iPt]/(fTotalBckYieldsNarrow[iPt]*fTotalBckYieldsNarrow[iPt])*fTotalBckYieldsNarrowError[iPt],2.) ,0.5);
			fMesonSignNarrow[iPt] = fMesonYieldsCorResidualBckFuncNarrow[iPt]/pow(fTotalBckYieldsNarrow[iPt],0.5);
			fMesonSignNarrowError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncNarrowError[iPt]/pow(fTotalBckYieldsNarrow[iPt],0.5),2.)+pow(0.5*fMesonYieldsCorResidualBckFuncNarrow[iPt]/pow(fTotalBckYieldsNarrow[iPt],1.5)*fTotalBckYieldsNarrowError[iPt],2.) ,0.5);

		}else{
			fMesonSBNarrow[iPt] = 0.;
			fMesonSBNarrowError[iPt] = 0.;
			fMesonSignNarrow[iPt] = 0.;
			fMesonSignNarrowError[iPt] = 0.;
		}

		//////////////////////////////// Start Analysis with  Normalization at the left of the Meson Peak

		// Function to subtract GG minus Bck
		ProcessEM( fHistoMappingGGInvMassPtBin[iPt], fHistoMappingBackInvMassPtBin[iPt], fBGFitRangeLeft);
		fHistoMappingSignalInvMassLeftPtBin[iPt] = fSignal;
		fHistoMappingBackNormInvMassLeftPtBin[iPt] = fBckNorm;


		cout<< "iPt"<< iPt<< " "<< "standard range"<<endl;
		// Fitting the subtracted spectra
		fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "normal range/left normalization" << endl;

		fFitInvMassLeftPtBin[iPt] =0x00;
		if(fCrysFitting==0){
			fFileErrLog << "Using exp fit"<<endl;
			FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE);
			fFitInvMassLeftPtBin[iPt] = fFitReco;
			fFitSignalPeakPosInvMassLeftPtBin[iPt] = fFitGausExp;
			fFitBckInvMassLeftPtBin[iPt] = fFitLinearBck;
			fMesonYieldsResidualBckFuncLeft[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncLeftError[iPt] = fIntLinearBckError;

		} else {
			fFileErrLog << "Using Crystal Ball function"<<endl;
			FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE, Form("CBFitFuncLeftBin%02d",iPt));
			fFitInvMassLeftPtBin[iPt] = fFitReco;
			fFitSignalPeakPosInvMassLeftPtBin[iPt] = fFitGausExp;
			fFitBckInvMassLeftPtBin[iPt] = fFitLinearBck;
			fMesonYieldsResidualBckFuncLeft[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncLeftError[iPt] = fIntLinearBckError;

		}
		//FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE);

		if (fFitInvMassLeftPtBin[iPt] !=0x00){
			fMesonMassLeft[iPt] = fFitInvMassLeftPtBin[iPt]->GetParameter(1);
			fMesonMassLeftError[iPt] = fFitInvMassLeftPtBin[iPt]->GetParError(1);
			
			fMesonWidthLeft[iPt] = fFitInvMassLeftPtBin[iPt]->GetParameter(2);
			fMesonWidthLeftError[iPt] = fFitInvMassLeftPtBin[iPt]->GetParError(2);
			/*
			fMesonCurLeftIntRange[0] = fMesonIntDeltaRange[0] - (fMesonMassExpect-fMesonMassLeft[iPt]);
			fMesonCurLeftIntRangeWide[0] = fMesonIntDeltaRangeWide[0] - (fMesonMassExpect-fMesonMassLeft[iPt]);
			fMesonCurLeftIntRangeNarrow[0] = fMesonIntDeltaRangeNarrow[0] - (fMesonMassExpect-fMesonMassLeft[iPt]);
			fMesonCurLeftIntRange[1] = fMesonIntDeltaRange[1] - (fMesonMassExpect-fMesonMassLeft[iPt]);
			fMesonCurLeftIntRangeWide[1] = fMesonIntDeltaRangeWide[1] - (fMesonMassExpect-fMesonMassLeft[iPt]);
			fMesonCurLeftIntRangeNarrow[1] = fMesonIntDeltaRangeNarrow[1] - (fMesonMassExpect-fMesonMassLeft[iPt]);
			*/
			fMesonCurLeftIntRange[0] 		= fMesonMassLeft[iPt] + fMesonIntDeltaRange[0];
			fMesonCurLeftIntRangeWide[0] 	= fMesonMassLeft[iPt] + fMesonIntDeltaRangeWide[0];
			fMesonCurLeftIntRangeNarrow[0] 	= fMesonMassLeft[iPt] + fMesonIntDeltaRangeNarrow[0];
			fMesonCurLeftIntRange[1] 		= fMesonMassLeft[iPt] + fMesonIntDeltaRange[1];
			fMesonCurLeftIntRangeWide[1] 	= fMesonMassLeft[iPt] + fMesonIntDeltaRangeWide[1];
			fMesonCurLeftIntRangeNarrow[1]	= fMesonMassLeft[iPt] + fMesonIntDeltaRangeNarrow[1];
			
		} else {
			fMesonMassLeft[iPt] = 0.;
			fMesonMassLeftError[iPt] = 0.;
			fMesonWidthLeft[iPt] = 0.;
			fMesonWidthLeftError[iPt] = 0.;
			/*
			fMesonCurLeftIntRange[0] = fMesonIntDeltaRange[0];
			fMesonCurLeftIntRangeWide[0] = fMesonIntDeltaRangeWide[0];
			fMesonCurLeftIntRangeNarrow[0] = fMesonIntDeltaRangeNarrow[0];
			fMesonCurLeftIntRange[1] = fMesonIntDeltaRange[1];
			fMesonCurLeftIntRangeWide[1] = fMesonIntDeltaRangeWide[1];
			fMesonCurLeftIntRangeNarrow[1] = fMesonIntDeltaRangeNarrow[1];
			*/
			fMesonCurLeftIntRange[0] 		= fMesonMassExpect + fMesonIntDeltaRange[0];
			fMesonCurLeftIntRangeWide[0] 	= fMesonMassExpect + fMesonIntDeltaRangeWide[0];
			fMesonCurLeftIntRangeNarrow[0] 	= fMesonMassExpect + fMesonIntDeltaRangeNarrow[0];
			fMesonCurLeftIntRange[1] 		= fMesonMassExpect + fMesonIntDeltaRange[1];
			fMesonCurLeftIntRangeWide[1] 	= fMesonMassExpect + fMesonIntDeltaRangeWide[1];
			fMesonCurLeftIntRangeNarrow[1] 	= fMesonMassExpect + fMesonIntDeltaRangeNarrow[1];
		}

		// Integrate the bck histo
		IntegrateHistoInvMass( fHistoMappingBackNormInvMassLeftPtBin[iPt], fMesonCurLeftIntRange);
		fBckYieldsLeft[iPt] = fYields;
		fBckYieldsLeftError[iPt] = fYieldsError;
		IntegrateHistoInvMass( fHistoMappingBackNormInvMassLeftPtBin[iPt], fMesonCurLeftIntRangeWide);
		fBckYieldsLeftWide[iPt] = fYields;
		fBckYieldsLeftWideError[iPt] = fYieldsError;
		IntegrateHistoInvMass( fHistoMappingBackNormInvMassLeftPtBin[iPt], fMesonCurLeftIntRangeNarrow);
		fBckYieldsLeftNarrow[iPt] = fYields;
		fBckYieldsLeftNarrowError[iPt] = fYieldsError;

		// Integrate the signal histo
		fFileDataLog<< endl <<"Signal histo normal range, left norm " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
		IntegrateHistoInvMassStream( fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonCurLeftIntRange);
		fMesonYieldsLeft[iPt] = fYields;
		fMesonYieldsLeftError[iPt] = fYieldsError;
		fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
		fFileDataLog<< endl <<"Signal histo wide range, left norm " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
		IntegrateHistoInvMassStream( fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonCurLeftIntRangeWide);
		fMesonYieldsLeftWide[iPt] = fYields;
		fMesonYieldsLeftWideError[iPt] = fYieldsError;
		fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
		fFileDataLog<< endl <<"Signal histo narrow range, left norm " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
		IntegrateHistoInvMassStream( fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonCurLeftIntRangeNarrow);
		fMesonYieldsLeftNarrow[iPt] = fYields;
		fMesonYieldsLeftNarrowError[iPt] = fYieldsError;
		fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;

		fFileDataLog<< "Residual Background leftover norm integration/left norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFuncLeft[iPt] <<"\t +- \t" << fMesonYieldsResidualBckFuncLeftError[iPt]<<endl<< endl;
		fMesonYieldsCorResidualBckFuncLeft[iPt] = fMesonYieldsLeft[iPt]- fMesonYieldsResidualBckFuncLeft[iPt];
		fMesonYieldsCorResidualBckFuncLeftError[iPt] =
		pow((fMesonYieldsLeftError[iPt]*fMesonYieldsLeftError[iPt]+
		fMesonYieldsResidualBckFuncLeftError[iPt]*fMesonYieldsResidualBckFuncLeftError[iPt]),0.5);
		fMesonYieldsLeftPerEvent[iPt]= fMesonYieldsCorResidualBckFuncLeft[iPt]/fNEvents;
		fMesonYieldsLeftPerEventError[iPt]= fMesonYieldsCorResidualBckFuncLeftError[iPt]/fNEvents;

		fTotalBckYieldsLeft[iPt] = fBckYieldsLeft[iPt] + fMesonYieldsResidualBckFuncLeft[iPt];
		fTotalBckYieldsLeftError[iPt] = pow(fBckYieldsLeftError[iPt]*fBckYieldsLeftError[iPt] + fMesonYieldsResidualBckFuncLeftError[iPt]*fMesonYieldsResidualBckFuncLeftError[iPt],0.5);

		//Integrate Fit Function
		IntegrateFitFunc( fFitSignalPeakPosInvMassLeftPtBin[iPt], fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonCurLeftIntRange);
		fMesonYieldsFuncLeft[iPt]=fYieldsFunc;

		//GetFWHM
		CalculateFWHM(fFitInvMassLeftPtBin[iPt]);
		cout<< "iPt left "<< iPt<< " "<< "B"<<endl;
		cout << fFWHMFunc << "     " << fFWHMFuncError << endl;
		fMesonFWHMLeft[iPt] = fFWHMFunc;
		cout<< "iPt left "<< iPt<< " "<< "C"<<endl;
		fMesonFWHMLeftError[iPt] = fFWHMFuncError;
		
		if( fFitBckInvMassLeftPtBin[iPt]->Integral(fMesonMassLeft[iPt]-fMesonFWHMLeft[iPt], fMesonMassLeft[iPt]+fFitInvMassLeftPtBin[iPt]->GetParameter(2))!=0){
			Double_t background = fFitBckInvMassLeftPtBin[iPt]->Integral(fMesonMassLeft[iPt]-fMesonFWHMLeft[iPt], fMesonMassLeft[iPt]+fFitInvMassLeftPtBin[iPt]->GetParameter(2));
			Double_t backgroundErr = fFitBckInvMassLeftPtBin[iPt]->IntegralError(fMesonMassLeft[iPt]-fMesonFWHMLeft[iPt], fMesonMassLeft[iPt]+fFitInvMassLeftPtBin[iPt]->GetParameter(2));
			Double_t signal = fFitInvMassLeftPtBin[iPt]->Integral(fMesonMassLeft[iPt]-fMesonFWHMLeft[iPt], fMesonMassLeft[iPt]+fFitInvMassLeftPtBin[iPt]->GetParameter(2)) - background;
			Double_t signalErr =pow( pow(fFitInvMassLeftPtBin[iPt]->IntegralError(fMesonMassLeft[iPt]-fMesonFWHMLeft[iPt], fMesonMassLeft[iPt]+fFitInvMassLeftPtBin[iPt]->GetParameter(2)),2 )+ pow(backgroundErr,2),0.5);
			fMesonSBLeft[iPt] = signal/ background;
			fMesonSBLeftError[iPt] = pow( pow(signalErr/background,2.)+pow(signal/(background *background )*backgroundErr ,2.) ,0.5); 
			fMesonSignLeft[iPt] = signal/ pow(background + signal,0.5);
			fMesonSignLeftError[iPt] = pow(pow( (pow(background + signal ,0.5) - 0.5*signal*pow(background+signal,-0.5))/(background+signal) * signalErr ,2) + pow( 0.5*pow(signal+background,-1.5),2) ,0.5); 
		}else{
			fMesonSBLeft[iPt] = 0.;
			fMesonSBLeftError[iPt] = 0.;
			fMesonSignLeft[iPt] = 0.;
			fMesonSignLeftError[iPt] = 0.;
		}

		// Wide integration mass window
		cout<< "iPt"<< iPt<< " "<< "wide range"<<endl;
		fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "wide range/left normalization" << endl;

		if(fCrysFitting==0){
			fFileErrLog << "Using exp fit"<<endl;
			FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntDeltaRangeWide,iPt,kFALSE);
			fMesonYieldsResidualBckFuncLeftWide[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncLeftWideError[iPt] = fIntLinearBckError;

		} else {
			fFileErrLog << "Using Crystal Ball function"<<endl;
			FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntDeltaRangeWide,iPt,kFALSE, Form("CBFitFuncLeftWideBin%02d",iPt));
			fMesonYieldsResidualBckFuncLeftWide[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncLeftWideError[iPt] = fIntLinearBckError;

		}

		//FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt],fMesonIntDeltaRangeWide,iPt,kFALSE);
		fFileDataLog<< "Residual Background leftover wide integration/left norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFuncLeftWide[iPt]<<"\t +- \t" << fMesonYieldsResidualBckFuncLeftWideError[iPt] <<endl<< endl;
		fMesonYieldsCorResidualBckFuncLeftWide[iPt] = fMesonYieldsLeftWide[iPt]- fMesonYieldsResidualBckFuncLeftWide[iPt];
		fMesonYieldsCorResidualBckFuncLeftWideError[iPt] =
		pow((fMesonYieldsLeftWideError[iPt]*fMesonYieldsLeftWideError[iPt]+ fMesonYieldsResidualBckFuncLeftWideError[iPt]*fMesonYieldsResidualBckFuncLeftWideError[iPt]),0.5);
		fMesonYieldsLeftPerEventWide[iPt]= fMesonYieldsCorResidualBckFuncLeftWide[iPt]/fNEvents;
		fMesonYieldsLeftPerEventWideError[iPt]= fMesonYieldsCorResidualBckFuncLeftWideError[iPt]/fNEvents;

		fTotalBckYieldsLeftWide[iPt] = fBckYieldsLeftWide[iPt] + fMesonYieldsResidualBckFuncLeftWide[iPt];
		fTotalBckYieldsLeftWideError[iPt] = pow(fBckYieldsLeftWideError[iPt]*fBckYieldsLeftWideError[iPt] + fMesonYieldsResidualBckFuncLeftWideError[iPt]*fMesonYieldsResidualBckFuncLeftWideError[iPt],0.5);


		if( fTotalBckYieldsLeftWide[iPt]!=0){
			fMesonSBLeftWide[iPt] = fMesonYieldsCorResidualBckFuncLeftWide[iPt]/fTotalBckYieldsLeftWide[iPt];
			fMesonSBLeftWideError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncLeftWideError[iPt]/fTotalBckYieldsLeftWide[iPt],2.)+pow(fMesonYieldsCorResidualBckFuncLeftWide[iPt]/(fTotalBckYieldsLeftWide[iPt]*fTotalBckYieldsLeftWide[iPt])*fTotalBckYieldsLeftWideError[iPt],2.) ,0.5);
			fMesonSignLeftWide[iPt] = fMesonYieldsCorResidualBckFuncLeftWide[iPt]/pow(fTotalBckYieldsLeftWide[iPt],0.5);
			fMesonSignLeftWideError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncLeftWideError[iPt]/pow(fTotalBckYieldsLeftWide[iPt],0.5),2.)+pow(0.5*fMesonYieldsCorResidualBckFuncLeftWide[iPt]/pow(fTotalBckYieldsLeftWide[iPt],1.5)*fTotalBckYieldsLeftWideError[iPt],2.) ,0.5);
		}else{
			fMesonSBLeftWide[iPt] = 0.;
			fMesonSBLeftWideError[iPt] = 0.;
			fMesonSignLeftWide[iPt] = 0.;
			fMesonSignLeftWideError[iPt] = 0.;
		}


		// Narrow integration mass window
		cout<< "iPt"<< iPt<< " "<< "narrow range"<<endl;
		fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "narrow range/left normalization" << endl;

		if(fCrysFitting==0){
			fFileErrLog << "Using exp fit"<<endl;
			FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntDeltaRangeNarrow,iPt,kFALSE);
			fMesonYieldsResidualBckFuncLeftNarrow[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncLeftNarrowError[iPt] = fIntLinearBckError;

		} else {
			fFileErrLog << "Using Crystal Ball function"<<endl;
			FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntDeltaRangeNarrow,iPt,kFALSE, Form("CBFitFuncLeftNarrowBin%02d",iPt));
			fMesonYieldsResidualBckFuncLeftNarrow[iPt] = fIntLinearBck;
			fMesonYieldsResidualBckFuncLeftNarrowError[iPt] = fIntLinearBckError;

		}

		//FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt],fMesonIntDeltaRangeNarrow,iPt,kFALSE);
		fFileDataLog<< "Residual Background leftover narrow integration/left norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFuncLeftNarrow[iPt]<<"\t +- \t" << fMesonYieldsResidualBckFuncLeftNarrowError[iPt] <<endl<< endl;
		fMesonYieldsCorResidualBckFuncLeftNarrow[iPt] = fMesonYieldsLeftNarrow[iPt]- fMesonYieldsResidualBckFuncLeftNarrow[iPt];
		fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt] =
		pow((fMesonYieldsLeftNarrowError[iPt]*fMesonYieldsLeftNarrowError[iPt]+
		fMesonYieldsResidualBckFuncLeftNarrowError[iPt]*fMesonYieldsResidualBckFuncLeftNarrowError[iPt]),0.5);
		fMesonYieldsLeftPerEventNarrow[iPt]= fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]/fNEvents;
		fMesonYieldsLeftPerEventNarrowError[iPt]= fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt]/fNEvents;

		fTotalBckYieldsLeftNarrow[iPt] = fBckYieldsLeftNarrow[iPt] + fMesonYieldsResidualBckFuncLeftNarrow[iPt];
		fTotalBckYieldsLeftNarrowError[iPt] = pow(fBckYieldsLeftNarrowError[iPt]*fBckYieldsLeftNarrowError[iPt] + fMesonYieldsResidualBckFuncLeftNarrowError[iPt]*fMesonYieldsResidualBckFuncLeftNarrowError[iPt],0.5);


		if( fTotalBckYieldsLeftNarrow[iPt]!=0){
			fMesonSBLeftNarrow[iPt] = fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]/fTotalBckYieldsLeftNarrow[iPt];
			fMesonSBLeftNarrowError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt]/fTotalBckYieldsLeftNarrow[iPt],2.)+pow(fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]/(fTotalBckYieldsLeftNarrow[iPt]*fTotalBckYieldsLeftNarrow[iPt])*fTotalBckYieldsLeftNarrowError[iPt],2.) ,0.5);
			fMesonSignLeftNarrow[iPt] = fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]/pow(fTotalBckYieldsLeftNarrow[iPt],0.5);
			fMesonSignLeftNarrowError[iPt] = pow(pow(fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt]/pow(fTotalBckYieldsLeftNarrow[iPt],0.5),2.)+pow(0.5*fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]/pow(fTotalBckYieldsLeftNarrow[iPt],1.5)*fTotalBckYieldsLeftNarrowError[iPt],2.) ,0.5);

		}else{
			fMesonSBLeftNarrow[iPt] = 0.;
			fMesonSBLeftNarrowError[iPt] = 0.;
			fMesonSignLeftNarrow[iPt] = 0.;
			fMesonSignLeftNarrowError[iPt] = 0.;
		}

	}



	//******************** Data OUTPUTFILE ***************************************************
	const char* fileNameSysErrDat = Form("%s/%s/%s_%s_SystematicErrorYieldExtractionDalitz_RAWDATA_%s.dat",fCutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPrefix2.Data(), cutSelection.Data());
	fstream fileSysErrDat;
	fileSysErrDat.open(fileNameSysErrDat, ios::out);
	fileSysErrDat << "Calculation of the systematic error due to the yield extraction RAWDATA " << endl;
	fileSysErrDat <<  endl;
	fileSysErrDat << "fGGYields" << endl;
	fileSysErrDat << "Bin \t Right \t Right Wide \t Right Narr" << endl;
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << fGGYields[iPt] << "+-" << fGGYieldsError[iPt] << "\t" <<
		fGGYieldsWide[iPt] << "+-" << fGGYieldsWideError[iPt] << "\t" <<
		fGGYieldsNarrow[iPt] << "+-" << fGGYieldsNarrowError[iPt] << endl;

	}
	fileSysErrDat <<  endl;
	fileSysErrDat << "fTotalBckYields" << endl;
	fileSysErrDat << "Bin \t Right \t Right Wide \t Right Narr \t Left \t Left Wide \t Left Narr" << endl;
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" <<
		fTotalBckYields[iPt] << "+-" << fTotalBckYieldsError[iPt] << "\t" <<
		fTotalBckYieldsWide[iPt] << "+-" << fTotalBckYieldsWideError[iPt] << "\t" <<
		fTotalBckYieldsNarrow[iPt] << "+-" << fTotalBckYieldsNarrowError[iPt] << "\t" <<
		fTotalBckYieldsLeft[iPt] << "+-" << fTotalBckYieldsLeftError[iPt]<< "\t" <<
		fTotalBckYieldsLeftWide[iPt]<< "+-" << fTotalBckYieldsLeftWideError[iPt]<< "\t" <<
		fTotalBckYieldsLeftNarrow[iPt]<< "+-" << fTotalBckYieldsLeftNarrowError[iPt] << endl;
	}
	fileSysErrDat <<  endl;
	fileSysErrDat << "fMesonYields" << endl;
	fileSysErrDat << "Bin \t Right \t Right Wide \t Right Narr \t Left \t Left Wide \t Left Narr" << endl;
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" <<
		fMesonYieldsCorResidualBckFunc[iPt] << "+-" << fMesonYieldsCorResidualBckFuncError[iPt] << "\t" <<
		fMesonYieldsCorResidualBckFuncWide[iPt] << "+-" << fMesonYieldsCorResidualBckFuncWideError[iPt] << "\t" <<
		fMesonYieldsCorResidualBckFuncNarrow[iPt] << "+-" << fMesonYieldsCorResidualBckFuncNarrowError[iPt] << "\t" <<
		fMesonYieldsCorResidualBckFuncLeft[iPt]<< "+-" << fMesonYieldsCorResidualBckFuncLeftError[iPt]<< "\t" <<
		fMesonYieldsCorResidualBckFuncLeftWide[iPt]<< "+-" << fMesonYieldsCorResidualBckFuncLeftWideError[iPt]<< "\t" <<
		fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]<< "+-" << fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt] << endl;
	}
	if(fIsMC){
		fileSysErrDat <<  endl;
		fileSysErrDat << "TrueYields" << endl;
		fileSysErrDat << "Bin \t True \t True Wide \t True Narr " << endl;
		for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
			fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" <<
			fMesonTrueYields[iPt] << "\t" <<
			fMesonTrueYieldsWide[iPt] << "\t" <<
			fMesonTrueYieldsNarrow[iPt] << endl;
		}
	}
	fileSysErrDat.close();
	//******************************** OUTPUT END ******************************************************

	TString nameMeson = Form("%s/%sDalitz_%s_MesonWithBck%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),fSuffix.Data());
	TString nameCanvas = "MesonWithBckCanvas";
	TString namePad = "MesonWithBckPad";
	//PlotInvMassInPtBins( fHistoMappingGGInvMassPtBin, fHistoMappingBackNormInvMassPtBin, nameMeson, nameCanvas, namePad, fMesonMassRange, date, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel,fTextCent,fCollisionSystem,kTRUE );
	PlotInvMassInPtBins( fHistoMappingGGInvMassPtBin, fHistoMappingBackNormInvMassPtBin, nameMeson, nameCanvas, namePad, fMesonMassRange, date, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel,fTextCent,fCollisionSystem,mode);
	
	
	TString nameMesonSub= Form("%s/%sDalitz_%s_MesonSubtracted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(), fPeriodFlag.Data(), fCutSelection.Data(),fSuffix.Data());
	TString nameCanvasSub= "MesonCanvasSubtracted";
	TString namePadSub= "MesonPadSubtracted";
	//PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalInvMassPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitSignalInvMassPtBin, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, date, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fTextCent);
	PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalInvMassPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitSignalInvMassPtBin, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, date, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem);

	
	

	nameMeson= Form("%s/%sDalitz_%s_MesonWithBckLeft%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),fSuffix.Data());
	nameCanvas = "MesonWithBckCanvasLeft";
	namePad = "MesonWithBckPadLeft";
	//PlotInvMassInPtBins( fHistoMappingGGInvMassPtBin, fHistoMappingBackNormInvMassLeftPtBin, nameMeson, nameCanvas, namePad,  fMesonMassRange, date, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC , fDecayChannel,fTextCent,fCollisionSystem,kTRUE);
	PlotInvMassInPtBins( fHistoMappingGGInvMassPtBin, fHistoMappingBackNormInvMassLeftPtBin, nameMeson, nameCanvas, namePad,  fMesonMassRange, date, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC , fDecayChannel,fTextCent,fCollisionSystem,mode);


	nameMesonSub= Form("%s/%sDalitz_%s_MesonSubtractedLeft%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),fSuffix.Data());
	nameCanvasSub= "MesonCanvasSubtractedLeft";
	namePadSub= "MesonPadSubtractedLeft";
	
	//PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalInvMassLeftPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitInvMassLeftPtBin, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, date, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel,fTextCent);
	PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalInvMassLeftPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitInvMassLeftPtBin, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, date, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem);


	PlotExampleInvMassBins(fHistoMappingGGInvMassPtBin[fExampleBin], fHistoMappingSignalInvMassPtBin[fExampleBin], fHistoMappingBackNormInvMassPtBin[fExampleBin] , fFitSignalInvMassPtBin[fExampleBin], fExampleBin, outputDir.Data(),fSuffix.Data(), fMesonMassRange, pictDrawingCoordinatesFWHM, fNEvents, date, fPrefix, fPrefix2, fThesis, fCollisionSystem, fBinsPt, fDecayChannel);
        //PlotExampleInvMassBinsWithandW0Background(fHistoMappingGGInvMassPtBin[fExampleBin], fHistoMappingSignalInvMassPtBin[fExampleBin], fHistoMappingBackNormInvMassPtBin[fExampleBin] , fFitSignalInvMassPtBin[fExampleBin], fExampleBin, outputDir.Data(),fSuffix.Data(), fMesonMassRange, pictDrawingCoordinatesFWHM, fNEvents, date, fPrefix, fPrefix2, fThesis, fCollisionSystem, fBinsPt, fDecayChannel);
       
	
	if(fIsMC){
	  
		TString nameMesonTrue= Form("%s/%sDalitz_%s_TrueMesonFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),fSuffix.Data());
		TString nameCanvasTrue= "TrueMesonCanvasFitted";
		TString namePadTrue= "TrueMesonPadFitted";
		PlotWithFitSubtractedInvMassInPtBins( fHistoMappingTrueMesonInvMassPtBins, fHistoMappingTrueMesonInvMassPtBins, fFitTrueSignalInvMassPtBin, nameMesonTrue, nameCanvasTrue, namePadTrue, fMesonMassRange, date, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem);

		
		nameMesonTrue= Form("%s/%sDalitz_%s_TrueGammaGammaContamination%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),fSuffix.Data());
		nameCanvasTrue= "TrueMesonCanvasGG";
		namePadTrue= "TrueMesonPadGG";
		PlotInvMassBckGGInPtBins( fHistoMappingTrueMesonInvMassPtBins, fHistoMappingTrueGGMesonInvMassPtBins, nameMesonTrue, nameCanvasTrue, namePadTrue, fMesonMassRange, date, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem );
		
		
		FillSignalInvMassW0TruePi0HistosArray();
		FillGGInvMassW0TruePi0HistosArray();
		
	        nameMesonTrue = Form("%s/%sDalitz_%s_SignalInvMass_TrueMesonSubstracted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),fSuffix.Data());
		nameCanvasTrue = "SignalInvMassTrueMesonSubstractedCanvas";
		namePadTrue = "SignalInvMassTrueMesonSubstractedPad";
		PlotSignalInvMassW0TrueMesonInPtBins(fHistoMappingSignalInvMassW0TruePi0PtBins, nameMesonTrue, nameCanvasTrue, namePadTrue, fMesonMassRange, date, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem);
	       
		
		nameMesonTrue = Form("%s/%sDalitz_%s_InvMass_TrueMesonSubstracted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),fSuffix.Data());
		nameCanvasTrue = "InvMassTrueMesonSubstractedCanvas";
		namePadTrue = "InvMassTrueMesonSubstractedPad";
		PlotInvMassW0TrueMesonInPtBins( fHistoMappingGGInvMassW0TruePi0PtBins, fHistoMappingBackNormInvMassPtBin, nameMesonTrue, nameCanvasTrue, namePadTrue, fMesonMassRange, date, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel,fTextCent,fCollisionSystem);
	
		
	        
		
		
	}

	CreatePtHistos();
	FillPtHistos();

	if(fIsMC){
	  
	  
		//FillHistos
	  
		FillHistosArrayMC(fHistoMCMesonPtWithinAcceptance, fHistoMCMesonPt, fDeltaPt);
	
		CalculateMesonAcceptance();
		cout << "Calculated MesonAcceptance" << endl;


                fNameHistoFrac="TrueGGFrac";
                fHistoYieldTrueGGFracMeson = CalculateSecondaryFractions(fHistoYieldTrueMeson, fHistoYieldTrueGGMeson, fNameHistoFrac);
                
                fNameHistoFrac="TrueGGFracNarrow";
                fHistoYieldTrueGGFracMesonNarrow= CalculateSecondaryFractions(fHistoYieldTrueMesonNarrow, fHistoYieldTrueGGMesonNarrow, fNameHistoFrac);
                
                fNameHistoFrac="TrueGGFracWide";
                fHistoYieldTrueGGFracMesonWide= CalculateSecondaryFractions(fHistoYieldTrueMesonWide, fHistoYieldTrueGGMesonWide, fNameHistoFrac);
		
		
		//////////////////////////////////////////////Create the scaling factor/////////////////////////////////////////////
		
		TH1D*  fHistoYieldTrueMesonScaled   = (TH1D*)fHistoYieldTrueMeson->Clone("fHistoYieldTrueMesonScaled");
		TH1D*  fHistoYieldTrueGGMesonScaled = (TH1D*)fHistoYieldTrueGGMeson->Clone("fHistoYieldTrueGGMesonScaled");
		
		fHistoYieldTrueMesonScaled->Scale(fPi0DalitzBRDPG/Pi0DalitzBR);
		fHistoYieldTrueGGMesonScaled->Scale(fPi0GGBRDPG/Pi0GGBR);
				  
		fNameHistoFrac="TrueGGFracForData";
                fHistoYieldTrueGGFracMesonForData = CalculateSecondaryFractions(fHistoYieldTrueMesonScaled, fHistoYieldTrueGGMesonScaled, fNameHistoFrac);
		
		
		TH1D* fHistoYieldTrueMesonNarrowScaled   = (TH1D*)fHistoYieldTrueMesonNarrow->Clone("fHistoYieldTrueMesonNarrowScaled");
		TH1D* fHistoYieldTrueGGMesonNarrowScaled = (TH1D*)fHistoYieldTrueGGMesonNarrow->Clone("fHistoYieldTrueGGMesonNarrowScaled");
		
		fHistoYieldTrueMesonNarrowScaled->Scale(fPi0DalitzBRDPG/Pi0DalitzBR);
		fHistoYieldTrueGGMesonNarrowScaled->Scale(fPi0GGBRDPG/Pi0GGBR);
                
                fNameHistoFrac="TrueGGFracNarrowForData";
                fHistoYieldTrueGGFracMesonNarrowForData= CalculateSecondaryFractions(fHistoYieldTrueMesonNarrowScaled, fHistoYieldTrueGGMesonNarrowScaled, fNameHistoFrac);
		
		TH1D* fHistoYieldTrueMesonWideScaled   = (TH1D*) fHistoYieldTrueMesonWide->Clone("fHistoYieldTrueMesonWideScaled");
		TH1D* fHistoYieldTrueGGMesonWideScaled = (TH1D*) fHistoYieldTrueGGMesonWide->Clone("fHistoYieldTrueGGMesonWideScaled");
		
		fHistoYieldTrueMesonWideScaled->Scale(fPi0DalitzBRDPG/Pi0DalitzBR);
		fHistoYieldTrueGGMesonWideScaled->Scale(fPi0GGBRDPG/Pi0GGBR );
		
                fNameHistoFrac="TrueGGFracWideForData";
                fHistoYieldTrueGGFracMesonWideForData= CalculateSecondaryFractions(fHistoYieldTrueMesonWideScaled, fHistoYieldTrueGGMesonWideScaled, fNameHistoFrac);
				
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		
		

		/////////////////////////////////////////
		fNameHistoEffi="MesonEffiPt";
		CalculateMesonEfficiencyWithoutGGCont(fHistoYieldMeson,fHistoYieldTrueGGFracMeson,fNameHistoEffi,"RAWYieldTrueGGFrac");
		fHistoMonteMesonEffiPt = fHistoMCMesonEffiPt;
		
		fNameHistoEffi="MesonWideEffiPt";
		CalculateMesonEfficiencyWithoutGGCont(fHistoYieldMesonWide,fHistoYieldTrueGGFracMesonWide,fNameHistoEffi,"RAWYieldWideTrueGGFrac");
		fHistoMonteMesonWideEffiPt = fHistoMCMesonEffiPt;
                
		fNameHistoEffi="MesonNarrowEffiPt";
		CalculateMesonEfficiencyWithoutGGCont(fHistoYieldMesonNarrow,fHistoYieldTrueGGFracMesonNarrow,fNameHistoEffi,"RAWYieldNarrowTrueGGFrac");
		fHistoMonteMesonNarrowEffiPt = fHistoMCMesonEffiPt;
		//////////////////////////////////////////
		
                

		fNameHistoEffi="MesonLeftEffiPt";
		CalculateMesonEfficiency(fHistoYieldMesonLeft,fNameHistoEffi);
		fHistoMonteMesonLeftEffiPt = fHistoMCMesonEffiPt;

		fNameHistoEffi="MesonLeftNarrowEffiPt";
		CalculateMesonEfficiency(fHistoYieldMesonLeftNarrow,fNameHistoEffi);
		fHistoMonteMesonLeftNarrowEffiPt = fHistoMCMesonEffiPt;

		fNameHistoEffi="MesonLeftWideEffiPt";
		CalculateMesonEfficiency(fHistoYieldMesonLeftWide,fNameHistoEffi);
		fHistoMonteMesonLeftWideEffiPt = fHistoMCMesonEffiPt;

		// True Meson (only once case, because no normalization
		fNameHistoEffi="TrueMesonEffiPt";
		CalculateMesonEfficiency(fHistoYieldTrueMeson,fNameHistoEffi);
		fHistoMCTrueMesonEffiPt = fHistoMCMesonEffiPt; 

		fNameHistoEffi="TrueMesonWideEffiPt";
		CalculateMesonEfficiency(fHistoYieldTrueMesonWide,fNameHistoEffi);
		fHistoMCTrueMesonWideEffiPt = fHistoMCMesonEffiPt;

		fNameHistoEffi="TrueMesonNarrowEffiPt";
		CalculateMesonEfficiency(fHistoYieldTrueMesonNarrow, fNameHistoEffi);
		fHistoMCTrueMesonNarrowEffiPt = fHistoMCMesonEffiPt;

		
		SaveCorrectionHistos(fCutSelection, fPrefix2);
	}
	SaveHistos(fIsMC, fCutSelection, fPrefix2);

	fFileErrLog.close();
	fFileDataLog.close();
	Delete();

}


TH1D* CalculateSecondaryFractions(TH1D* histoRawYield, TH1D* histoRawYieldSec, TString nameHistoFrac/*, Bool_t scaleGGCont*/){

        TH1D* histoFracSec = (TH1D*)histoRawYieldSec->Clone(nameHistoFrac.Data());
        TH1D* histoRawYieldTotal = (TH1D*)histoRawYield->Clone();
        histoRawYieldTotal->Add(histoFracSec,1.);
        histoFracSec->Divide(histoFracSec,histoRawYieldTotal,1.,1.,"B");
	
	return histoFracSec;

}




void ProcessEM(TH1D* fGammaGamma, TH1D* fBck, Double_t * fBGFitRangeEM)
{
	
	for (Int_t binx= 0; binx < fGammaGamma->GetNbinsX()+1; binx++){
		if(fGammaGamma->GetBinContent(binx) == 0){
			fGammaGamma->SetBinError(binx,1.);
			fGammaGamma->SetBinContent(binx,0.);
		}
	}

		  
	fBckNorm = (TH1D*)fBck->Clone("fBckNorm");

		
	fGammaGamma->Sumw2();
	fBck->Sumw2();
	fBckNorm->Sumw2();

	Double_t 	r= fGammaGamma->Integral(fGammaGamma->GetXaxis()->FindBin(fBGFitRangeEM[0]),fGammaGamma->GetXaxis()->FindBin(fBGFitRangeEM[1]));
	Double_t 	b= fBck->Integral(fBck->GetXaxis()->FindBin(fBGFitRangeEM[0]),fBck->GetXaxis()->FindBin(fBGFitRangeEM[1]));
	Double_t 	norm = 1;

	if(b != 0) norm = r/b;
	fBckNorm->Scale(norm);
	cout<<"r="<<r<<" b="<<b<<" r/b="<<r/b<< " " << endl;

	Int_t numberOfZeros = 0;
	for (Int_t i = 1; i < fBckNorm->GetNbinsX()+1; i++){
		if (fBckNorm->GetBinContent(i) == 0){
			numberOfZeros++;
			if (norm > 1.){
				fBckNorm->SetBinError(i,1.);
				fBckNorm->SetBinContent(i,0.);
			}
		}
	}
	cout<<"counted " << numberOfZeros << " in the normalized BG" << endl;

	fSignal = (TH1D*)fGammaGamma->Clone("fSignal");
	fSignal->Sumw2();
	
	fSignal->Add(fBckNorm,-1.);
	
}


void ProduceBckProperWeighting(TList* fBackgroundContainer,TList* fMotherContainer, TList *fESDContainer){
   
   THnSparseF* fSparseMotherZM = 0;
   THnSparseF* fSparseBckZM = 0;


   if( !fMotherContainer && fBackgroundContainer ) {
   
      fSparseMotherZM 	= (THnSparseF*)fBackgroundContainer->FindObject("Back_Mother_InvMass_Pt_z_m");
      fSparseBckZM 	= (THnSparseF*)fBackgroundContainer->FindObject("Back_Back_InvMass_Pt_z_m");
   }

   else if( fMotherContainer && fBackgroundContainer ){

      fSparseMotherZM 	= (THnSparseF*)fMotherContainer->FindObject("Back_Mother_InvMass_Pt_z_m");
      fSparseBckZM 	= (THnSparseF*)fBackgroundContainer->FindObject("Back_Back_InvMass_Pt_z_m");

   }
   
   cout<<"Entro ruido"<<endl;
   
   
   if( fSparseMotherZM &&  fSparseBckZM ) {
 
   
   for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
     
      fHistoWeightsBGZbinVsMbin[iPt] = new  TH2F("BGWeights", "", fSparseMotherZM->GetAxis(2)->GetNbins(),  0, fSparseMotherZM->GetAxis(2)->GetNbins(),fSparseMotherZM->GetAxis(3)->GetNbins(),  0, fSparseMotherZM->GetAxis(3)->GetNbins());
      fHistoWeightsBGZbinVsMbin[iPt]->GetYaxis()->SetTitle("M-bins");
      fHistoWeightsBGZbinVsMbin[iPt]->GetXaxis()->SetTitle("Z-bins");
      fHistoWeightsBGZbinVsMbin[iPt]->Sumw2();
      fHistoFillPerEventBGZbinVsMbin[iPt] = new  TH2F("BGPoolsFillstatus", "", fSparseMotherZM->GetAxis(2)->GetNbins(),  0, fSparseMotherZM->GetAxis(2)->GetNbins(),fSparseMotherZM->GetAxis(3)->GetNbins(),  0, fSparseMotherZM->GetAxis(3)->GetNbins());
      fHistoFillPerEventBGZbinVsMbin[iPt]->GetYaxis()->SetTitle("M-bins");
      fHistoFillPerEventBGZbinVsMbin[iPt]->GetXaxis()->SetTitle("Z-bins");
      fHistoFillPerEventBGZbinVsMbin[iPt]->Sumw2();
   
      for (Int_t z=0;z < fSparseMotherZM->GetAxis(2)->GetNbins()-1;z++){
         for (Int_t m = 0; m < fSparseMotherZM->GetAxis(3)->GetNbins(); m++){ 
            // pt
            fSparseMotherZM->GetAxis(1)->SetRange((fSparseMotherZM->GetAxis(1))->FindBin(fBinsPt[iPt]+0.001),(fSparseMotherZM->GetAxis(1))->FindBin(fBinsPt[iPt+1]-0.001));
            fSparseBckZM->GetAxis(1)->SetRange((fSparseBckZM->GetAxis(1))->FindBin(fBinsPt[iPt]+0.001),(fSparseBckZM->GetAxis(1))->FindBin(fBinsPt[iPt+1]-0.001));
            // z
            fSparseMotherZM->GetAxis(2)->SetRange(z, z);
            fSparseBckZM->GetAxis(2)->SetRange(z, z);
            // m
            fSparseMotherZM->GetAxis(3)->SetRange(m,m);
            fSparseBckZM->GetAxis(3)->SetRange(m,m);
            
            fHistoMotherZMProj = (TH1D*)fSparseMotherZM->Projection(0);
            fHistoMotherZMProj->Sumw2();
            fHistoBckZMProj = (TH1D*)fSparseBckZM->Projection(0);
            fHistoBckZMProj->Sumw2();
        
            fScalingFactorBck[z][m]= 1./fBackgroundMultNumber;
            if (m==0 && z ==0){
               if(fHistoMappingBackInvMassPtBin[iPt]!= NULL){
                  delete fHistoMappingBackInvMassPtBin[iPt];
                  fHistoMappingBackInvMassPtBin[iPt]=NULL;
               }
               fNameHistoBack = Form("Mapping_Back_InvMass_in_Pt_Bin%02d", iPt);
               fHistoMappingBackInvMassPtBin[iPt]= (TH1D*)fHistoBckZMProj->Clone(fNameHistoBack);
               fHistoMappingBackInvMassPtBin[iPt]->Sumw2();
               for (Int_t ii = 0; ii < fHistoBckZMProj->GetNbinsX()+1; ii++){
                  fHistoMappingBackInvMassPtBin[iPt]->SetBinContent(ii,0.);
                  fHistoMappingBackInvMassPtBin[iPt]->SetBinError(ii,0.);
               }
            }
            Int_t startBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[0]);
            Int_t endBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[1]);
            if (fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral) != 0) {
               fScalingFactorBck[z][m] = fHistoMotherZMProj->Integral(startBinIntegral,endBinIntegral)/fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral);
//                cout << z << "\t" << m << "\t"  << fScalingFactorBck[z][m] << "\t" << 20./fBackgroundMultNumber << endl;
               if ( fScalingFactorBck[z][m]> (20./fBackgroundMultNumber) ){
//                   cout << "fail safe entered" << endl;
                  fScalingFactorBck[z][m]=1./fBackgroundMultNumber;
//                   cout << "\t" <<  z << "\t" << m << "\t"  << fScalingFactorBck[z][m] << "\t" << 20./fBackgroundMultNumber << endl;
               }
            }
            fHistoMappingBackInvMassPtBin[iPt]->Add(fHistoBckZMProj,fScalingFactorBck[z][m]);
            fHistoWeightsBGZbinVsMbin[iPt]->Fill(z+0.5,m+0.5,fScalingFactorBck[z][m]);
            fHistoFillPerEventBGZbinVsMbin[iPt]->Fill(z+0.5,m+0.5,fHistoBckZMProj->GetEntries());
            fHistoMotherZMProj->Clear();
            fHistoBckZMProj->Clear();

         }
      }
      fHistoMappingBackInvMassPtBin[iPt]->Rebin(fNRebin[iPt]);
      //fHistoMappingBackInvMassPtBin[iPt]->Scale(1./fNRebin[iPt]);
      for (Int_t ii = 0; ii < fHistoMappingBackInvMassPtBin[iPt]->GetNbinsX()+1; ii++){
         if(fHistoMappingBackInvMassPtBin[iPt]->GetBinContent(ii) == 0){
            fHistoMappingBackInvMassPtBin[iPt]->SetBinContent(ii,0.);
            fHistoMappingBackInvMassPtBin[iPt]->SetBinError(ii,1.);
         }
      }
      fFileDataLog << "Scaling Background factors for Pt bin " << iPt << " z m " << endl;
      // 		cout << "Scaling Background factors for Pt bin " << iPt << " z m " << endl;
      for (Int_t z=0; z < fSparseMotherZM->GetAxis(2)->GetNbins()-1; z++){
         fFileDataLog << fScalingFactorBck[z][0] << "\t" << fScalingFactorBck[z][1] << "\t" << fScalingFactorBck[z][2] << "\t" << fScalingFactorBck[z][3] << endl;
         // 			cout << fScalingFactorBck[z][0] << "\t" << fScalingFactorBck[z][1] << "\t" << fScalingFactorBck[z][2] << "\t" << fScalingFactorBck[z][3] << endl;
      }
   }
   
   for (Int_t z=0;z < fSparseMotherZM->GetAxis(2)->GetNbins()-1;z++){
      for (Int_t m = 0; m < fSparseMotherZM->GetAxis(3)->GetNbins(); m++) {

         // pt
         fSparseMotherZM->GetAxis(1)->SetRange((fSparseMotherZM->GetAxis(1))->FindBin(fFullPt[0]+0.001),(fSparseMotherZM->GetAxis(1))->FindBin(fFullPt[1]-0.001));
         fSparseBckZM->GetAxis(1)->SetRange((fSparseBckZM->GetAxis(1))->FindBin(fFullPt[0]+0.001),(fSparseBckZM->GetAxis(1))->FindBin(fFullPt[1]-0.001));
         // z
         fSparseMotherZM->GetAxis(2)->SetRange(z+1, z+1);
         fSparseBckZM->GetAxis(2)->SetRange(z+1, z+1);
         // m
         fSparseMotherZM->GetAxis(3)->SetRange(m+1,m+1);
         fSparseBckZM->GetAxis(3)->SetRange(m+1,m+1);
         
         fHistoMotherZMProj = (TH1D*)fSparseMotherZM->Projection(0);
         fHistoMotherZMProj->Sumw2();
         fHistoBckZMProj = (TH1D*)fSparseBckZM->Projection(0);
         fHistoBckZMProj->Sumw2();
         
         fScalingFactorBck[z][m]= 1./fBackgroundMultNumber;
         if (m==0 && z ==0){
            fNameHistoBack = "Mapping_Back_InvMass_FullPt";
            fMesonFullPtBackground = (TH1D*)fHistoBckZMProj->Clone(fNameHistoBack);
            fMesonFullPtBackground->Sumw2();
            for (Int_t ii = 0; ii < fHistoBckZMProj->GetNbinsX()+1; ii++){
               fMesonFullPtBackground->SetBinContent(ii,0.);
               fMesonFullPtBackground->SetBinError(ii,0.);
            }
         }
         Int_t startBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[0]);
         Int_t endBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[1]);
         if (fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral) != 0) {
            fScalingFactorBck[z][m] = fHistoMotherZMProj->Integral(startBinIntegral,endBinIntegral)/fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral);
            if ( fScalingFactorBck[z][m]>20./fBackgroundMultNumber ){
               fScalingFactorBck[z][m]=1./fBackgroundMultNumber;
            }
         }
         fMesonFullPtBackground->Add(fHistoBckZMProj,fScalingFactorBck[z][m]);

         fHistoMotherZMProj->Clear();
         fHistoBckZMProj->Clear();
      }
   }
   fMesonFullPtBackground->Rebin(fNRebin[4]);
   //fMesonFullPtBackground->Scale(1./fNRebin[4]);
   

   for (Int_t z=0;z < fSparseMotherZM->GetAxis(2)->GetNbins()-1;z++){
      for (Int_t m = 0; m < fSparseMotherZM->GetAxis(3)->GetNbins(); m++) {
        
         // pt
         fSparseMotherZM->GetAxis(1)->SetRange((fSparseMotherZM->GetAxis(1))->FindBin(fMidPt[0]+0.001),(fSparseMotherZM->GetAxis(1))->FindBin(fMidPt[1]-0.001));
         fSparseBckZM->GetAxis(1)->SetRange((fSparseBckZM->GetAxis(1))->FindBin(fMidPt[0]+0.001),(fSparseBckZM->GetAxis(1))->FindBin(fMidPt[1]-0.001));
         // z
         fSparseMotherZM->GetAxis(2)->SetRange(z+1, z+1);
         fSparseBckZM->GetAxis(2)->SetRange(z+1, z+1);
         // m
         fSparseMotherZM->GetAxis(3)->SetRange(m+1,m+1);
         fSparseBckZM->GetAxis(3)->SetRange(m+1,m+1);
         
         fHistoMotherZMProj = (TH1D*)fSparseMotherZM->Projection(0);
         fHistoMotherZMProj->Sumw2();
         fHistoBckZMProj = (TH1D*)fSparseBckZM->Projection(0);
         fHistoBckZMProj->Sumw2();
         
         fScalingFactorBck[z][m]= 1./fBackgroundMultNumber;
         if (m==0 && z ==0){
            fNameHistoBack = "Mapping_Back_InvMass_MidPt";
            fFittingHistMidPtBackground = (TH1D*)fHistoBckZMProj->Clone(fNameHistoBack);
            fFittingHistMidPtBackground->Sumw2();
            for (Int_t ii = 0; ii < fHistoBckZMProj->GetNbinsX()+1; ii++){
               fFittingHistMidPtBackground->SetBinContent(ii,0.);
               fFittingHistMidPtBackground->SetBinError(ii,0.);
            }
         }
         Int_t startBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[0]);
         Int_t endBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[1]);
         if (fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral) != 0) {
            fScalingFactorBck[z][m] = fHistoMotherZMProj->Integral(startBinIntegral,endBinIntegral)/fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral);
            if ( fScalingFactorBck[z][m]>20./fBackgroundMultNumber ){
               fScalingFactorBck[z][m]=1./fBackgroundMultNumber;
            }
         }
         fFittingHistMidPtBackground->Add(fHistoBckZMProj,fScalingFactorBck[z][m]);
         
         fHistoMotherZMProj->Clear();
         fHistoBckZMProj->Clear();
      }
   }
   fFittingHistMidPtBackground->Rebin(fNRebin[4]);
   //fFittingHistMidPtBackground->Scale(1./fNRebin[4]);
   
   }  else {
   
		cout << "Using TH2 for the background" << endl;
		for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		  
			cout<<"Entro bucle"<<endl;
			fHistoMotherZM = (TH2D*)fESDContainer->FindObject("ESD_Mother_InvMass_Pt");
			fHistoBckZM = (TH2D*)fESDContainer->FindObject("ESD_Background_InvMass_Pt");
			
			cout<<"Salio bucle"<<endl;

			Int_t startBin = fHistoMotherZM->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
			Int_t endBin = fHistoMotherZM->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

			fHistoMotherZMProj = fHistoMotherZM->ProjectionX("ProjectMother",startBin,endBin);
			fHistoBckZMProj = fHistoBckZM->ProjectionX("ProjectBck",startBin,endBin);

			fScalingFactorBck[0][0]= 1./fBackgroundMultNumber;
				if(fHistoMappingBackInvMassPtBin[iPt]!= NULL){
					delete fHistoMappingBackInvMassPtBin[iPt];
					fHistoMappingBackInvMassPtBin[iPt]=NULL;
				}
			fNameHistoBack = Form("Mapping_Back_InvMass_in_Pt_Bin%02d", iPt);
			fHistoMappingBackInvMassPtBin[iPt]= (TH1D*)fHistoBckZMProj->Clone(fNameHistoBack);
			fHistoMappingBackInvMassPtBin[iPt]->Sumw2();
			
			Int_t startBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[0]);
			Int_t endBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[1]);
			if (fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral) != 0) {
				fScalingFactorBck[0][0] = fHistoMotherZMProj->Integral(startBinIntegral,endBinIntegral)/fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral);
				if ( fScalingFactorBck[0][0]>20./fBackgroundMultNumber ){
					fScalingFactorBck[0][0]=1./fBackgroundMultNumber;
				}
			}
			fHistoMappingBackInvMassPtBin[iPt]->Add(fHistoBckZMProj,fScalingFactorBck[0][0]);
			fHistoMappingBackInvMassPtBin[iPt]->Rebin(fNRebin[iPt]);
			//fHistoMappingBackInvMassPtBin[iPt]->Scale(1./fNRebin[iPt]);
			for (Int_t ii = 0; ii < fHistoMappingBackInvMassPtBin[iPt]->GetNbinsX()+1; ii++){
				if(fHistoMappingBackInvMassPtBin[iPt]->GetBinContent(ii) == 0){
					fHistoMappingBackInvMassPtBin[iPt]->SetBinContent(ii,0.);
					fHistoMappingBackInvMassPtBin[iPt]->SetBinError(ii,1.);
				}
			}			

			fFileDataLog << "Scaling Background factors for Pt bin " << iPt << endl;
			fFileDataLog << fScalingFactorBck[0][0] << endl;
			
		}

		fHistoMotherZM = (TH2D*)fESDContainer->FindObject("ESD_Mother_InvMass_Pt");
		fHistoBckZM = (TH2D*)fESDContainer->FindObject("ESD_Background_InvMass_Pt");

		Int_t startBinFullPt = fHistoMotherZM->GetYaxis()->FindBin(fFullPt[0]+0.001);
		Int_t endBinFullPt = fHistoMotherZM->GetYaxis()->FindBin(fFullPt[1]-0.001);

		fHistoMotherZMProjFullPt = fHistoMotherZM->ProjectionX("ProjectMother",startBinFullPt,endBinFullPt);
		fHistoBckZMProjFullPt = fHistoBckZM->ProjectionX("ProjectBck",startBinFullPt,endBinFullPt);

		fScalingFactorBckFullPt= 1./fBackgroundMultNumber;
		fNameHistoBack = "Mapping_Back_InvMass_FullPt";
		fMesonFullPtBackground = (TH1D*)fHistoBckZMProjFullPt->Clone(fNameHistoBack);
		fMesonFullPtBackground->Sumw2();
		for (Int_t ii = 0; ii < fHistoBckZMProjFullPt->GetNbinsX()+1; ii++){
				fMesonFullPtBackground->SetBinContent(ii,0.);
				fMesonFullPtBackground->SetBinError(ii,0.);
		}
		Int_t startBinIntegralFullPt = fHistoMotherZMProjFullPt->GetXaxis()->FindBin(fBGFitRange[0]);
		Int_t endBinIntegralFullPt = fHistoMotherZMProjFullPt->GetXaxis()->FindBin(fBGFitRange[1]);
		if (fHistoBckZMProjFullPt->Integral(startBinIntegralFullPt,endBinIntegralFullPt) != 0) {
			fScalingFactorBckFullPt = fHistoMotherZMProjFullPt->Integral(startBinIntegralFullPt,endBinIntegralFullPt)/fHistoBckZMProjFullPt->Integral(startBinIntegralFullPt,endBinIntegralFullPt);
			if ( fScalingFactorBckFullPt>20./fBackgroundMultNumber ){
				fScalingFactorBckFullPt=1./fBackgroundMultNumber;
			}
		}
		fMesonFullPtBackground->Add(fHistoBckZMProjFullPt,fScalingFactorBckFullPt);

		fMesonFullPtBackground->Rebin(fNRebin[4]);
		//fMesonFullPtBackground->Scale(1./fNRebin[4]);

		Int_t startBinMidPt = fHistoMotherZM->GetYaxis()->FindBin(fMidPt[0]+0.001);
		Int_t endBinMidPt = fHistoMotherZM->GetYaxis()->FindBin(fMidPt[1]-0.001);

		fHistoMotherZMProjMidPt = fHistoMotherZM->ProjectionX("ProjectMother",startBinMidPt,endBinMidPt);
		fHistoBckZMProjMidPt = fHistoBckZM->ProjectionX("ProjectBck",startBinMidPt,endBinMidPt);

		fScalingFactorBckMidPt= 1./fBackgroundMultNumber;
				
		fNameHistoBack = "Mapping_Back_InvMass_MidPt";
		fFittingHistMidPtBackground = (TH1D*)fHistoBckZMProjMidPt->Clone(fNameHistoBack);
		fFittingHistMidPtBackground->Sumw2();
		for (Int_t ii = 0; ii < fHistoBckZMProjMidPt->GetNbinsX()+1; ii++){
			fFittingHistMidPtBackground->SetBinContent(ii,0.);
			fFittingHistMidPtBackground->SetBinError(ii,0.);
		}
				
		Int_t startBinIntegralMidPt = fHistoMotherZMProjMidPt->GetXaxis()->FindBin(fBGFitRange[0]);
		Int_t endBinIntegralMidPt = fHistoMotherZMProjMidPt->GetXaxis()->FindBin(fBGFitRange[1]);
		if (fHistoBckZMProjMidPt->Integral(startBinIntegralMidPt,endBinIntegralMidPt) != 0) {
			fScalingFactorBckMidPt = fHistoMotherZMProjMidPt->Integral(startBinIntegralMidPt,endBinIntegralMidPt)/fHistoBckZMProjMidPt->Integral(startBinIntegralMidPt,endBinIntegralMidPt);
			if ( fScalingFactorBckMidPt>20./fBackgroundMultNumber ){
				fScalingFactorBckMidPt=1./fBackgroundMultNumber;
			}
		}
		fFittingHistMidPtBackground->Add(fHistoBckZMProjMidPt,fScalingFactorBckMidPt);

		fFittingHistMidPtBackground->Rebin(fNRebin[4]);
		//fFittingHistMidPtBackground->Scale(1./fNRebin[4]);
	}
   
   
      
   
   
}

void ProduceBckWithoutWeighting(TH2F *fBckInvMassVSPt){
	//calculation background for midPt without weighting

        cout<<"Entro Axes"<<endl;
	Int_t startBinMidPt = fBckInvMassVSPt->GetYaxis()->FindBin(fMidPt[0]+0.001);
	Int_t endBinMidPt = fBckInvMassVSPt->GetYaxis()->FindBin(fMidPt[1]-0.001);
 	fFittingHistMidPtBackground = new TH1D("Mapping_Back_InvMass_MidPt","Mapping_Back_InvMass_MidPt",fBckInvMassVSPt->GetNbinsX(),0.,.8);
	fBckInvMassVSPt->ProjectionX("Mapping_Back_InvMass_MidPt",startBinMidPt,endBinMidPt);
        cout<<"Entro1 Axes"<<endl;

	fFittingHistMidPtBackground=(TH1D*)gDirectory->Get("Mapping_Back_InvMass_MidPt");
	fFittingHistMidPtBackground->Rebin(fNRebin[4]);

	cout << " created BG MidPt" << endl;
	
	//calulation background for fullPt without weighting
 	fMesonFullPtBackground = new TH1D("Mapping_Back_InvMass_FullPt","Mapping_Back_InvMass_FullPt",fBckInvMassVSPt->GetNbinsX(),0.,.8);
	Int_t startBinFullPt = fBckInvMassVSPt->GetYaxis()->FindBin(fFullPt[0]+0.001);
	Int_t endBinFullPt = fBckInvMassVSPt->GetYaxis()->FindBin(fFullPt[1]-0.001);
	fBckInvMassVSPt->ProjectionX("Mapping_Back_InvMass_FullPt",startBinFullPt,endBinFullPt);
	fMesonFullPtBackground=(TH1D*)gDirectory->Get("Mapping_Back_InvMass_FullPt");
	fMesonFullPtBackground->Rebin(fNRebin[4]);

	cout << " created BG FullPt" << endl;
	
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fNameHistoBack = Form("Mapping_Back_InvMass_in_Pt_Bin%02d", iPt);
		if(fHistoMappingBackInvMassPtBin[iPt]!= NULL){
			delete fHistoMappingBackInvMassPtBin[iPt];
			fHistoMappingBackInvMassPtBin[iPt]=NULL;
		}
 		fHistoMappingBackInvMassPtBin[iPt]=new TH1D(fNameHistoBack.Data(),fNameHistoBack.Data(),fBckInvMassVSPt->GetNbinsX(),0.,.8);
		Int_t startBin = fBckInvMassVSPt->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
		Int_t endBin = fBckInvMassVSPt->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

		fBckInvMassVSPt->ProjectionX(fNameHistoBack.Data(),startBin,endBin);
		fHistoMappingBackInvMassPtBin[iPt]=(TH1D*)gDirectory->Get(fNameHistoBack.Data());
		if(fNRebin[iPt]>1){
			fHistoMappingBackInvMassPtBin[iPt]->Rebin(fNRebin[iPt]);
		}

		cout << " created BG "<< fNameHistoBack.Data() << endl;
		
	}
}


void FillMassHistosArray(TH2F* fGammaGammaInvMassVSPt)
{

 	fFittingHistMidPtSignal = new TH1D("Mapping_GG_InvMass_MidPt","Mapping_GG_InvMass_MidPt",fGammaGammaInvMassVSPt->GetNbinsX(),0.,0.8);
	
	Int_t startBinMidPt = fGammaGammaInvMassVSPt->GetYaxis()->FindBin(fMidPt[0]+0.001);
	Int_t endBinMidPt = fGammaGammaInvMassVSPt->GetYaxis()->FindBin(fMidPt[1]-0.001);

	fGammaGammaInvMassVSPt->ProjectionX("Mapping_GG_InvMass_MidPt",startBinMidPt,endBinMidPt);
	
	fFittingHistMidPtSignal=(TH1D*)gDirectory->Get("Mapping_GG_InvMass_MidPt");
	fFittingHistMidPtSignal->Rebin(fNRebin[4]);
	cout << "Mid pt geschrieben" << endl;

 	fMesonFullPtSignal = new TH1D("Mapping_GG_InvMass_FullPt","Mapping_GG_InvMass_FullPt",fGammaGammaInvMassVSPt->GetNbinsX(),0.,0.8);

	Int_t startBinFullPt = fGammaGammaInvMassVSPt->GetYaxis()->FindBin(fFullPt[0]+0.001);
	Int_t endBinFullPt   = fGammaGammaInvMassVSPt->GetYaxis()->FindBin(fFullPt[1]-0.001);

	fGammaGammaInvMassVSPt->ProjectionX("Mapping_GG_InvMass_FullPt",startBinFullPt,endBinFullPt);

	fMesonFullPtSignal=(TH1D*)gDirectory->Get("Mapping_GG_InvMass_FullPt");
	fMesonFullPtSignal->Rebin(fNRebin[4]);
	//fMesonFullPtSignal->Scale(1./fNRebin[4]);
	cout << "Full pt geschrieben" << endl;
	
	Double_t totIntegral = 0.0;

	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){

		fNameHistoGG = Form("Mapping_GG_InvMass_in_Pt_Bin%02d", iPt);

		if(fHistoMappingGGInvMassPtBin[iPt]!= NULL){
			delete fHistoMappingGGInvMassPtBin[iPt];
			fHistoMappingGGInvMassPtBin[iPt]=NULL;
		}
 		fHistoMappingGGInvMassPtBin[iPt]=new TH1D(fNameHistoGG.Data(),fNameHistoGG.Data(),fGammaGammaInvMassVSPt->GetNbinsX(),0.,0.8);


		Int_t startBin = fGammaGammaInvMassVSPt->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
		Int_t endBin = fGammaGammaInvMassVSPt->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

		cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
		cout<< "bin values::"<< fGammaGammaInvMassVSPt->GetYaxis()->GetBinCenter(startBin)<< " "
		<< fGammaGammaInvMassVSPt->GetYaxis()->GetBinCenter(endBin)<< endl;

		fGammaGammaInvMassVSPt->ProjectionX(fNameHistoGG.Data(),startBin,endBin);
		

		fHistoMappingGGInvMassPtBin[iPt]=(TH1D*)gDirectory->Get(fNameHistoGG.Data());
		

		if(fNRebin[iPt]>1){
			fHistoMappingGGInvMassPtBin[iPt]->Rebin(fNRebin[iPt]);
			//fHistoMappingGGInvMassPtBin[iPt]->Scale(1./fNRebin[iPt]);
		}
		
		Double_t tempIntegral = fHistoMappingGGInvMassPtBin[iPt]->Integral();
	
		cout<< "Integral:: "<<tempIntegral<<endl;
		totIntegral += tempIntegral;

		
	}
	cout << "each pt written" << endl;
	cout << "Total sum:: "<<totIntegral<<endl;

}

void FillEposEnegHistosArray(TH2F* fEPosEnegInvMassVSPt)
{

          
 		
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){

		TString fNameHistoEE = Form("Mapping_EposEneg_InvMass_in_Pt_Bin%02d", iPt);

		if(fHistoMappingEEInvMassPtBin[iPt]!= NULL){
		  
			delete fHistoMappingEEInvMassPtBin[iPt];
			       fHistoMappingEEInvMassPtBin[iPt]=NULL;
			
		}
		
 		fHistoMappingEEInvMassPtBin[iPt]=new TH1D(fNameHistoEE.Data(),fNameHistoEE.Data(),fEPosEnegInvMassVSPt->GetNbinsX(),0.,5.0);


		Int_t startBin = fEPosEnegInvMassVSPt->GetYaxis()->FindBin( fBinsPt[iPt]   + 0.001);
		Int_t endBin   = fEPosEnegInvMassVSPt->GetYaxis()->FindBin( fBinsPt[iPt+1] - 0.001);

		cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
		cout<< "bin values::"<< fEPosEnegInvMassVSPt->GetYaxis()->GetBinCenter(startBin)<< " "<< fEPosEnegInvMassVSPt->GetYaxis()->GetBinCenter(endBin)<< endl;

		
		fEPosEnegInvMassVSPt->ProjectionX(fNameHistoEE.Data(),startBin,endBin);


		fHistoMappingEEInvMassPtBin[iPt]=(TH1D*)gDirectory->Get(fNameHistoEE.Data());

		//if(fNRebin[iPt]>1){
		  
		//	fHistoMappingEEInvMassPtBin[iPt]->Rebin(fNRebin[iPt]);
		//}
        }
	
}

void FillMassMCTrueMesonHistosArray(TH2F* fHistoTrueMesonInvMassVSPtFill)
{
	Double_t totIntegral = 0.0;
      
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fNameHistoTrue = Form("Mapping_TrueMeson_InvMass_in_Pt_Bin%02d", iPt);
		if(fHistoMappingTrueMesonInvMassPtBins[iPt]!= NULL){
			delete fHistoMappingTrueMesonInvMassPtBins[iPt];
			fHistoMappingTrueMesonInvMassPtBins[iPt]=NULL;
		}
		fHistoMappingTrueMesonInvMassPtBins[iPt] = new TH1D(fNameHistoTrue.Data(),fNameHistoTrue.Data(),fHistoTrueMesonInvMassVSPtFill->GetNbinsX(),0.,0.8);

		Int_t startBin = fHistoTrueMesonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
		Int_t endBin = fHistoTrueMesonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

		cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
		cout<< "bin values::"<< fHistoTrueMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(startBin)<< " "
		<< fHistoTrueMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(endBin)<< endl;

		fHistoTrueMesonInvMassVSPtFill->ProjectionX(fNameHistoTrue.Data(),startBin,endBin);
		fHistoMappingTrueMesonInvMassPtBins[iPt]=(TH1D*)gDirectory->Get(fNameHistoTrue.Data());
		if(fNRebin[iPt]>1){
			fHistoMappingTrueMesonInvMassPtBins[iPt]->Rebin(fNRebin[iPt]);
		}
		Double_t tempIntegral = fHistoMappingTrueMesonInvMassPtBins[iPt]->Integral();
		
		cout<<"Integral:: "<<tempIntegral<<endl;
		totIntegral+=tempIntegral;
		
		fHistoMappingTrueMesonInvMassPtBins[iPt]->SetLineWidth(1);
		fHistoMappingTrueMesonInvMassPtBins[iPt]->SetLineColor(2);
	}
	cout<<"Total sum true: "<<totIntegral<<endl;
}

void FillMassMCTrueGGMesonHistosArray(TH2F* fHistoTrueGGMesonInvMassVSPtFill)
{
	for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
		fNameHistoTrueGG = Form("Mapping_TrueGGMeson_InvMass_in_Pt_Bin%02d", iPt);
		if(fHistoMappingTrueGGMesonInvMassPtBins[iPt]!= NULL){
			delete fHistoMappingTrueGGMesonInvMassPtBins[iPt];
			fHistoMappingTrueGGMesonInvMassPtBins[iPt]=NULL;
		}
		fHistoMappingTrueGGMesonInvMassPtBins[iPt] = new TH1D(fNameHistoTrueGG.Data(),fNameHistoTrueGG.Data(),fHistoTrueGGMesonInvMassVSPtFill->GetNbinsX(),0.,0.8);

		Int_t startBin = fHistoTrueGGMesonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
		Int_t endBin = fHistoTrueGGMesonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

		cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
		cout<< "bin values::"<< fHistoTrueGGMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(startBin)<< " "
		<< fHistoTrueGGMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(endBin)<< endl;

		fHistoTrueGGMesonInvMassVSPtFill->ProjectionX(fNameHistoTrueGG.Data(),startBin,endBin);
		fHistoMappingTrueGGMesonInvMassPtBins[iPt]=(TH1D*)gDirectory->Get(fNameHistoTrueGG.Data());
		if(fNRebin[iPt]>1){
			fHistoMappingTrueGGMesonInvMassPtBins[iPt]->Rebin(fNRebin[iPt]);
			//fHistoMappingTrueGGMesonInvMassPtBins[iPt]->Scale(1./fNRebin[iPt]);

		}
		fHistoMappingTrueGGMesonInvMassPtBins[iPt]->SetLineWidth(1);
		fHistoMappingTrueGGMesonInvMassPtBins[iPt]->SetLineColor(2);
	}
}

void FillMassMCTrueMesonGGBckHistosArray(TH2F* fHistoTrueGGBckInvMassVSPtFill)
{
        for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
                fNameHistoTrueGGBck = Form("Mapping_TrueGGBck_InvMass_in_Pt_Bin%02d", iPt);

                if( fHistoMappingTrueGGBckInvMassPtBins[iPt]!= NULL ){
                        delete fHistoMappingTrueGGBckInvMassPtBins[iPt];
                        fHistoMappingTrueGGBckInvMassPtBins[iPt]=NULL;
                }

                fHistoMappingTrueGGBckInvMassPtBins[iPt] = new TH1D(fNameHistoTrueGGBck.Data(),fNameHistoTrueGGBck.Data(),fHistoTrueGGBckInvMassVSPtFill->GetNbinsX(),0.,0.8);

                Int_t startBin = fHistoTrueGGBckInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
                Int_t endBin   = fHistoTrueGGBckInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

                cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
                cout<< "bin values::"<< fHistoTrueGGBckInvMassVSPtFill->GetYaxis()->GetBinCenter(startBin)<< " "
                << fHistoTrueGGBckInvMassVSPtFill->GetYaxis()->GetBinCenter(endBin)<< endl;

                fHistoTrueGGBckInvMassVSPtFill->ProjectionX(fNameHistoTrueGGBck.Data(),startBin,endBin);
                fHistoMappingTrueGGBckInvMassPtBins[iPt]=(TH1D*)gDirectory->Get(fNameHistoTrueGGBck.Data());
                if(fNRebin[iPt]>1){
                        fHistoMappingTrueGGBckInvMassPtBins[iPt]->Rebin(fNRebin[iPt]);

                }
                fHistoMappingTrueGGBckInvMassPtBins[iPt]->SetLineWidth(1);
                fHistoMappingTrueGGBckInvMassPtBins[iPt]->SetLineColor(2);
                fHistoMappingTrueGGBckInvMassPtBins[iPt]->GetXaxis()->SetRangeUser(fMesonMassRange[0],fMesonMassRange[1]);
        }
}

void FillMassMCTrueMesonBckContHistosArray(TH2F* fHistoTrueContBckInvMassVSPtFill)
{
        for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
                fNameHistoTrueContBck = Form("Mapping_TrueContBck_InvMass_in_Pt_Bin%02d", iPt);
                //fHistoMappingTrueContBckInvMassPtBins
                if(fHistoMappingTrueContBckInvMassPtBins[iPt]!= NULL){
                        delete fHistoMappingTrueContBckInvMassPtBins[iPt];
                        fHistoMappingTrueContBckInvMassPtBins[iPt]=NULL;
                }
                fHistoMappingTrueContBckInvMassPtBins[iPt] = new TH1D(fNameHistoTrueContBck.Data(),fNameHistoTrueContBck.Data(),fHistoTrueContBckInvMassVSPtFill->GetNbinsX(),0.,0.8);

                Int_t startBin = fHistoTrueContBckInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
                Int_t endBin = fHistoTrueContBckInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

                cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
                cout<< "bin values::"<< fHistoTrueContBckInvMassVSPtFill->GetYaxis()->GetBinCenter(startBin)<< " "
                << fHistoTrueContBckInvMassVSPtFill->GetYaxis()->GetBinCenter(endBin)<< endl;

                fHistoTrueContBckInvMassVSPtFill->ProjectionX(fNameHistoTrueContBck.Data(),startBin,endBin);
                fHistoMappingTrueContBckInvMassPtBins[iPt]=(TH1D*)gDirectory->Get(fNameHistoTrueContBck.Data());
                if(fNRebin[iPt]>1){
                        fHistoMappingTrueContBckInvMassPtBins[iPt]->Rebin(fNRebin[iPt]);
                }
                fHistoMappingTrueContBckInvMassPtBins[iPt]->SetLineWidth(1);
                fHistoMappingTrueContBckInvMassPtBins[iPt]->SetLineColor(2);
                fHistoMappingTrueContBckInvMassPtBins[iPt]->GetXaxis()->SetRangeUser(fMesonMassRange[0],fMesonMassRange[1]);
        }
}

void FillSignalInvMassW0TruePi0HistosArray()
{
  
      for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
	
                fNameHistoSignalInvMassW0TruePi0 = Form("Mapping_Signal_InvMass_WOTruePi0_Pt_Bin%02d", iPt);
               
                if(fHistoMappingSignalInvMassW0TruePi0PtBins[iPt]!= NULL){
                        delete fHistoMappingSignalInvMassW0TruePi0PtBins[iPt];
                        fHistoMappingSignalInvMassW0TruePi0PtBins[iPt]=NULL;
                }
                 
		fHistoMappingSignalInvMassW0TruePi0PtBins[iPt] = (TH1D*)fHistoMappingSignalInvMassPtBin[iPt]->Clone(fNameHistoSignalInvMassW0TruePi0.Data());
		fHistoMappingSignalInvMassW0TruePi0PtBins[iPt]->Add(fHistoMappingTrueMesonInvMassPtBins[iPt],-1);
		fHistoMappingSignalInvMassW0TruePi0PtBins[iPt]->Add(fHistoMappingTrueGGMesonInvMassPtBins[iPt],-1); 
		

      }
    
}

void FillGGInvMassW0TruePi0HistosArray()
{
  
      for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
	
                fNameHistoGGInvMassW0TruePi0 = Form("Mapping_GG_InvMass_WOTruePi0_Pt_Bin%02d", iPt);
               
                if(fHistoMappingGGInvMassW0TruePi0PtBins[iPt]!= NULL){
                        delete fHistoMappingGGInvMassW0TruePi0PtBins[iPt];
                        fHistoMappingGGInvMassW0TruePi0PtBins[iPt]=NULL;
                }
                 
		fHistoMappingGGInvMassW0TruePi0PtBins[iPt] = (TH1D*)fHistoMappingGGInvMassPtBin[iPt]->Clone(fNameHistoGGInvMassW0TruePi0.Data());
		fHistoMappingGGInvMassW0TruePi0PtBins[iPt]->Add(fHistoMappingTrueMesonInvMassPtBins[iPt],-1);
		fHistoMappingGGInvMassW0TruePi0PtBins[iPt]->Add(fHistoMappingTrueGGMesonInvMassPtBins[iPt],-1); 
		

      }
    
}

void CreatePtHistos(){

	fDeltaPt =			 	new TH1D("deltaPt","",fNBinsPt,fBinsPt);

	fHistoYieldMeson = 			new TH1D("histoYieldMeson","",fNBinsPt,fBinsPt);
	fHistoYieldTrueMeson = 		        new TH1D("histoYieldTrueMeson","",fNBinsPt,fBinsPt);
	fHistoYieldTrueGGMeson = 		new TH1D("histoYieldTrueGGMeson","",fNBinsPt,fBinsPt);
	
	fHistoYieldMesonPerEvent = 	new TH1D("histoYieldMesonPerEvent","",fNBinsPt,fBinsPt);
	fHistoSignMeson = 			new TH1D("histoSignMeson","",fNBinsPt,fBinsPt);
	fHistoSBMeson = 			new TH1D("histoSBMeson","",fNBinsPt,fBinsPt);
	fHistoMassMeson = 			new TH1D("histoMassMeson","",fNBinsPt,fBinsPt);
	fHistoWidthMeson = 			new TH1D("histoWidthMeson","",fNBinsPt,fBinsPt);
	fHistoFWHMMeson = 			new TH1D("histoFWHMMeson","",fNBinsPt,fBinsPt);
	fHistoTrueMassMeson = 		new TH1D("histoTrueMassMeson","",fNBinsPt,fBinsPt);
	fHistoTrueFWHMMeson = 		new TH1D("histoTrueFWHMMeson","",fNBinsPt,fBinsPt);
	fHistoTrueSignMeson = 		new TH1D("histoTrueSignMeson","",fNBinsPt,fBinsPt);
	fHistoTrueSBMeson = 		new TH1D("histoTrueSBMeson","",fNBinsPt,fBinsPt);

	fHistoYieldMesonNarrow = 	new TH1D("histoYieldMesonNarrow","",fNBinsPt,fBinsPt);
	fHistoYieldTrueMesonNarrow = 	new TH1D("histoYieldTrueMesonNarrow","",fNBinsPt,fBinsPt);
	fHistoYieldTrueGGMesonNarrow = 		new TH1D("histoYieldTrueGGMesonNarrow","",fNBinsPt,fBinsPt);
	fHistoYieldMesonPerEventNarrow = new TH1D("histoYieldMesonPerEventNarrow","",fNBinsPt,fBinsPt);
	fHistoSignMesonNarrow = 		new TH1D("histoSignMesonNarrow","",fNBinsPt,fBinsPt);
	fHistoSBMesonNarrow = 		new TH1D("histoSBMesonNarrow","",fNBinsPt,fBinsPt);

	fHistoYieldMesonWide = 		new TH1D("histoYieldMesonWide","",fNBinsPt,fBinsPt);
	fHistoYieldTrueMesonWide = 	new TH1D("histoYieldTrueMesonWide","",fNBinsPt,fBinsPt);
	fHistoYieldTrueGGMesonWide = 		new TH1D("histoYieldTrueGGMesonWide","",fNBinsPt,fBinsPt);
	fHistoYieldMesonPerEventWide = new TH1D("histoYieldMesonPerEventWide","",fNBinsPt,fBinsPt);
	fHistoSignMesonWide = 		new TH1D("histoSignMesonWide","",fNBinsPt,fBinsPt);
	fHistoSBMesonWide = 		new TH1D("histoSBMesonWide","",fNBinsPt,fBinsPt);

	// Histos for normalization at the left of the peak

	fHistoYieldMesonLeft = 		new TH1D("histoYieldMesonLeft","",fNBinsPt,fBinsPt);
	fHistoYieldMesonLeftPerEvent = new TH1D("histoYieldMesonLeftPerEvent","",fNBinsPt,fBinsPt);
	fHistoSignMesonLeft = 		new TH1D("histoSignMesonLeft","",fNBinsPt,fBinsPt);
	fHistoSBMesonLeft = 		new TH1D("histoSBMesonLeft","",fNBinsPt,fBinsPt);
	fHistoMassMesonLeft = 		new TH1D("histoMassMesonLeft","",fNBinsPt,fBinsPt);
	fHistoWidthMesonLeft = 		new TH1D("histoWidthMesonLeft","",fNBinsPt,fBinsPt);
	fHistoFWHMMesonLeft = 		new TH1D("histoFWHMMesonLeft","",fNBinsPt,fBinsPt);


	fHistoYieldMesonLeftNarrow = 	new TH1D("histoYieldMesonLeftNarrow","",fNBinsPt,fBinsPt);
	fHistoYieldMesonLeftPerEventNarrow = new TH1D("histoYieldMesonLeftPerEventNarrow","",fNBinsPt,fBinsPt);
	fHistoSignMesonLeftNarrow = 	new TH1D("histoSignMesonLeftNarrow","",fNBinsPt,fBinsPt);
	fHistoSBMesonLeftNarrow = 	new TH1D("histoSBMesonLeftNarrow","",fNBinsPt,fBinsPt);


	fHistoYieldMesonLeftWide = 	new TH1D("histoYieldMesonLeftWide","",fNBinsPt,fBinsPt);
	fHistoYieldMesonLeftPerEventWide = new TH1D("histoYieldMesonLeftPerEventWide","",fNBinsPt,fBinsPt);
	fHistoSignMesonLeftWide = 	new TH1D("histoSignMesonLeftWide","",fNBinsPt,fBinsPt);
	fHistoSBMesonLeftWide = 		new TH1D("histoSBMesonLeftWide","",fNBinsPt,fBinsPt);

}

void FillPtHistos()
{
	for(Int_t iPt=fStartPtBin+1;iPt<fNBinsPt+1;iPt++){

		fDeltaPt->SetBinContent(iPt,fBinsPt[iPt]-fBinsPt[iPt-1]);
		fDeltaPt->SetBinError(iPt,0);


		fHistoMassMeson->SetBinContent(iPt,fMesonMass[iPt-1]);
		fHistoMassMeson->SetBinError(iPt,fMesonMassError[iPt-1]);
		// fHistoWidthMeson->SetBinContent(iPt);
		fHistoFWHMMeson->SetBinContent(iPt,fMesonFWHM[iPt-1]);
		fHistoFWHMMeson->SetBinError(iPt,fMesonFWHMError[iPt-1]);

		if (fIsMC) {
			fHistoTrueMassMeson->SetBinContent(iPt,fMesonTrueMass[iPt-1]);
			fHistoTrueMassMeson->SetBinError(iPt,fMesonTrueMassError[iPt-1]);
			fHistoTrueFWHMMeson->SetBinContent(iPt,fMesonTrueFWHM[iPt-1]);
			fHistoTrueFWHMMeson->SetBinError(iPt,fMesonTrueFWHMError[iPt-1]);
			fHistoTrueSignMeson->SetBinContent(iPt,fMesonTrueSign[iPt-1]);
			fHistoTrueSignMeson->SetBinError(iPt,fMesonTrueSignError[iPt-1]);
			fHistoTrueSBMeson->SetBinContent(iPt,fMesonTrueSB[iPt-1]);
			fHistoTrueSBMeson->SetBinError(iPt,fMesonTrueSBError[iPt-1]);
		}

		fHistoSignMeson->SetBinContent(iPt,fMesonSign[iPt-1]);
		fHistoSignMeson->SetBinError(iPt,fMesonSignError[iPt-1]);
		fHistoSBMeson->SetBinContent(iPt,fMesonSB[iPt-1]);
		fHistoSBMeson->SetBinError(iPt,fMesonSBError[iPt-1]);

		fHistoYieldMeson->SetBinContent(iPt,fMesonYieldsCorResidualBckFunc[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMeson->SetBinError(iPt,fMesonYieldsCorResidualBckFuncError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		if (fIsMC) {


			fHistoYieldTrueMeson->SetBinContent(iPt,fMesonTrueYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueMeson->SetBinError(iPt,fMesonTrueYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueGGMeson->SetBinContent(iPt,fMesonTrueGGYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueGGMeson->SetBinError(iPt,fMesonTrueGGYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));


		}
		fHistoYieldMesonPerEvent->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFunc[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonPerEvent->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		// Narrow integration window
		fHistoSignMesonNarrow->SetBinContent(iPt,fMesonSignNarrow[iPt-1]);
		fHistoSignMesonNarrow->SetBinError(iPt,fMesonSignNarrowError[iPt-1]);
		fHistoSBMesonNarrow->SetBinContent(iPt,fMesonSBNarrow[iPt-1]);
		fHistoSBMesonNarrow->SetBinError(iPt,fMesonSBNarrowError[iPt-1]);

		fHistoYieldMesonNarrow->SetBinContent(iPt,fMesonYieldsCorResidualBckFuncNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonNarrow->SetBinError(iPt,fMesonYieldsCorResidualBckFuncNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		if (fIsMC) {
			fHistoYieldTrueMesonNarrow->SetBinContent(iPt,fMesonTrueYieldsNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueMesonNarrow->SetBinError(iPt,fMesonTrueYieldsNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueGGMesonNarrow->SetBinContent(iPt,fMesonTrueGGYieldsNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueGGMesonNarrow->SetBinError(iPt,fMesonTrueGGYieldsNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		}
		fHistoYieldMesonPerEventNarrow->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonPerEventNarrow->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		// Wide integration window
		fHistoSignMesonWide->SetBinContent(iPt,fMesonSignWide[iPt-1]);
		fHistoSignMesonWide->SetBinError(iPt,fMesonSignWideError[iPt-1]);
		fHistoSBMesonWide->SetBinContent(iPt,fMesonSBWide[iPt-1]);
		fHistoSBMesonWide->SetBinError(iPt,fMesonSBWideError[iPt-1]);

		fHistoYieldMesonWide->SetBinContent(iPt,fMesonYieldsCorResidualBckFuncWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonWide->SetBinError(iPt,fMesonYieldsCorResidualBckFuncWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		if (fIsMC) {
			fHistoYieldTrueMesonWide->SetBinContent(iPt,fMesonTrueYieldsWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueMesonWide->SetBinError(iPt,fMesonTrueYieldsWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueGGMesonWide->SetBinContent(iPt,fMesonTrueGGYieldsWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
			fHistoYieldTrueGGMesonWide->SetBinError(iPt,fMesonTrueGGYieldsWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		}
		fHistoYieldMesonPerEventWide->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonPerEventWide->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		// Histos for integration at the left of the peak
		fHistoMassMesonLeft->SetBinContent(iPt,fMesonMassLeft[iPt-1]);
		fHistoMassMesonLeft->SetBinError(iPt,fMesonMassLeftError[iPt-1]);
		fHistoFWHMMesonLeft->SetBinContent(iPt,fMesonFWHMLeft[iPt-1]);
		fHistoFWHMMesonLeft->SetBinError(iPt,fMesonFWHMLeftError[iPt-1]);

		fHistoSignMesonLeft->SetBinContent(iPt,fMesonSignLeft[iPt-1]);
		fHistoSignMesonLeft->SetBinError(iPt,fMesonSignLeftError[iPt-1]);
		fHistoSBMesonLeft->SetBinContent(iPt,fMesonSBLeft[iPt-1]);
		fHistoSBMesonLeft->SetBinError(iPt,fMesonSBLeftError[iPt-1]);

		fHistoYieldMesonLeft->SetBinContent(iPt,fMesonYieldsCorResidualBckFuncLeft[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeft->SetBinError(iPt,fMesonYieldsCorResidualBckFuncLeftError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftPerEvent->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeft[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftPerEvent->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeftError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		// Narrow integration window
		fHistoSignMesonLeftNarrow->SetBinContent(iPt,fMesonSignLeftNarrow[iPt-1]);
		fHistoSignMesonLeftNarrow->SetBinError(iPt,fMesonSignLeftNarrowError[iPt-1]);
		fHistoSBMesonLeftNarrow->SetBinContent(iPt,fMesonSBLeftNarrow[iPt-1]);
		fHistoSBMesonLeftNarrow->SetBinError(iPt,fMesonSBLeftNarrowError[iPt-1]);

		fHistoYieldMesonLeftNarrow->SetBinContent(iPt,fMesonYieldsCorResidualBckFuncLeftNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftNarrow->SetBinError(iPt,fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftPerEventNarrow->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeftNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftPerEventNarrow->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

		// Wide integration window
		fHistoSignMesonLeftWide->SetBinContent(iPt,fMesonSignLeftWide[iPt-1]);
		fHistoSignMesonLeftWide->SetBinError(iPt,fMesonSignLeftWideError[iPt-1]);
		fHistoSBMesonLeftWide->SetBinContent(iPt,fMesonSBLeftWide[iPt-1]);
		fHistoSBMesonLeftWide->SetBinError(iPt,fMesonSBLeftWideError[iPt-1]);

		fHistoYieldMesonLeftWide->SetBinContent(iPt,fMesonYieldsCorResidualBckFuncLeftWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftWide->SetBinError(iPt,fMesonYieldsCorResidualBckFuncLeftWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftPerEventWide->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeftWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
		fHistoYieldMesonLeftPerEventWide->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeftWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
	}
}


void FitSubtractedInvMassInPtBins(TH1D* fHistoMappingSignalInvMassPtBinSingle, Double_t* fMesonIntDeltaRangeFit, Int_t ptBin, Bool_t vary)
{

	cout<<"Start Fitting spectra"<<endl;
	fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassRange[0],fMesonMassRange[1]);
	Double_t mesonAmplitude =fHistoMappingSignalInvMassPtBinSingle->GetMaximum();
	Double_t mesonAmplitudeMin;
	Double_t mesonAmplitudeMax;
	//NOTE ask ANA
	if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
		mesonAmplitudeMin = mesonAmplitude*98./100.;
		mesonAmplitudeMax = mesonAmplitude*115./100.;
		if( fEnergyFlag.CompareTo("HI") == 0) mesonAmplitudeMin = mesonAmplitude*92./100.;
	} else {
		mesonAmplitudeMin = mesonAmplitude*85./100.;
		mesonAmplitudeMax = mesonAmplitude*115./100.;
	}
	fFitReco = new TF1("GaussExpLinear","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)",fMesonFitRange[0],fMesonFitRange[1]);


	fFitGausExp = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);

	fFitLinearBck = new TF1("Linear","[0]+[1]*x",fMesonFitRange[0],fMesonFitRange[1]);


	fFitReco->SetParameter(0,mesonAmplitude);
	fFitReco->SetParameter(1,fMesonMassExpect);
	fFitReco->SetParameter(2,fMesonWidthExpect);
	
	
	if(vary){
		fFitReco->SetParameter(3,fMesonLambdaTail);
	} else {
		fFitReco->FixParameter(3,fMesonLambdaTail);
	}
	fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
	
	fFitReco->SetParLimits(1,fMesonMassRange[0],fMesonMassRange[1]);
	//fFitReco->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.15);
	
	
	fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);
	
	if(vary){
	  fFitReco->SetParLimits(3,fMesonLambdaTailRange[0],fMesonLambdaTailRange[1]);
	 }

	fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
	fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");

	fFitReco->SetLineColor(3);
	fFitReco->SetLineWidth(1);
	fFitReco->SetLineStyle(1);

	if (vary) fMesonLambdaTail = fFitReco->GetParameter(3);

	fFitGausExp->SetParameter(0,fFitReco->GetParameter(0));
	fFitGausExp->SetParameter(1,fFitReco->GetParameter(1));
	fFitGausExp->SetParameter(2,fFitReco->GetParameter(2));
	fFitGausExp->SetParameter(3,fFitReco->GetParameter(3));

	fFitGausExp->SetParError(0,fFitReco->GetParError(0));
	fFitGausExp->SetParError(1,fFitReco->GetParError(1));
	fFitGausExp->SetParError(2,fFitReco->GetParError(2));
	fFitGausExp->SetParError(3,fFitReco->GetParError(3));

	fFitLinearBck->SetParameter(0,fFitReco->GetParameter(4));
	fFitLinearBck->SetParameter(1,fFitReco->GetParameter(5));

	fFitLinearBck->SetParError(0,fFitReco->GetParError(4));
	fFitLinearBck->SetParError(1,fFitReco->GetParError(5));

	Int_t binCenterStart;
	Double_t startBinEdge;
	Int_t binCenterEnd;
	Double_t endBinEdge;

	TVirtualFitter * fitter = TVirtualFitter::GetFitter();

	if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("PROBLEMS  ") == 0 ){
	  
	  
	  
		//binCenterStart = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntDeltaRangeFit[0]-(fMesonMassExpect-fFitReco->GetParameter(1)));
		//startBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterStart)- 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
		//binCenterEnd = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntDeltaRangeFit[1]-(fMesonMassExpect-fFitReco->GetParameter(1)));
		//endBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterEnd)+ 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
		
		
		
		binCenterStart = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+fMesonIntDeltaRangeFit[0]);
		startBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterStart)- 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
		binCenterEnd = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+fMesonIntDeltaRangeFit[1]);
		endBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterEnd)+ 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);

			
		

		Int_t nFreePar = fFitReco->GetNumberFreeParameters();
		double * covMatrix = fitter->GetCovarianceMatrix();

		Float_t intLinearBack = fFitLinearBck->GetParameter(0)*(endBinEdge-startBinEdge)+
		0.5*fFitLinearBck->GetParameter(1)*(endBinEdge*endBinEdge-startBinEdge*startBinEdge);

		Float_t errorLinearBck = pow((pow( (endBinEdge-startBinEdge)*fFitReco->GetParError(4),2)+pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fFitReco->GetParError(5),2)+2*covMatrix[nFreePar*nFreePar-2]*(endBinEdge-startBinEdge)*0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5);

		fFileDataLog << "Parameter for bin " << ptBin << endl;
		fFileDataLog << "Gausexp: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<< "\t "<< fFitReco->GetParameter(3) <<"+-" << fFitReco->GetParError(3)<<endl;
		fFileDataLog << "Linear: \t"<<fFitReco->GetParameter(4)<<"+-" << fFitReco->GetParError(4) << "\t "<<fFitReco->GetParameter(5) <<"+-" << fFitReco->GetParError(5)<< endl;

		fIntLinearBck = intLinearBack/fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
		fIntLinearBckError = errorLinearBck/fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
	} else {
		fFileErrLog << "Fitting failed in test " << ptBin << " with status" << gMinuit->fCstatu.Data()<<"p" <<endl << endl;
	}
	fFitReco->DrawCopy("same");
}

void FitTrueInvMassInPtBins(TH1D* fHistoMappingSignalInvMassPtBinSingle, Double_t* fMesonIntDeltaRangeFit, Int_t ptBin)
{
	cout<<"Start Fitting spectra"<<endl;
	fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassRange[0],fMesonMassRange[1]);
	Double_t mesonAmplitude =fHistoMappingSignalInvMassPtBinSingle->GetMaximum();
	Double_t mesonAmplitudeMin;
	Double_t mesonAmplitudeMax;
	mesonAmplitudeMin = mesonAmplitude*98./100.;
	mesonAmplitudeMax = mesonAmplitude*115./100.;
	
	if (fMode == 6 || fMode == 7){
			mesonAmplitudeMin = mesonAmplitude*20./100.;
			mesonAmplitudeMax = mesonAmplitude*1000./100.;
	}
	

	fFitReco = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);

	fFitReco->SetParameter(0,mesonAmplitude);
	fFitReco->SetParameter(1,fMesonMassExpect);
	fFitReco->SetParameter(2,fMesonWidthExpect);
	fFitReco->SetParameter(3,fMesonLambdaTail);

	fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
	fFitReco->SetParLimits(1,fMesonMassRange[0],fMesonMassRange[1]);
	fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);
	fFitReco->SetParLimits(3,fMesonLambdaTailRange[0],fMesonLambdaTailRange[1]);

	fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
	fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");

	fFitReco->SetLineColor(3);
	fFitReco->SetLineWidth(1);
	fFitReco->SetLineStyle(1);

	//Int_t binCenterStart;
	//Double_t startBinEdge;
	//Int_t binCenterEnd;
	//Double_t endBinEdge;

	if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 ){
	  
		//binCenterStart = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntDeltaRangeFit[0]-(fMesonMassExpect-fFitReco->GetParameter(1)));
		//startBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterStart)- 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
		//binCenterEnd = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntDeltaRangeFit[1]-(fMesonMassExpect-fFitReco->GetParameter(1)));
		//endBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterEnd)+ 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);

		fFileDataLog << "Parameter for bin " << ptBin << endl;
		fFileDataLog << "Gausexp: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<< "\t "<< fFitReco->GetParameter(3) <<"+-" << fFitReco->GetParError(3)<<endl;
		
	} else {
		fFileErrLog << "Fitting failed in " << ptBin << " with status " << gMinuit->fCstatu.Data() <<endl << endl;
	}
	fFitReco->DrawCopy("same");
}


void FitCBSubtractedInvMassInPtBins(TH1D* fHistoMappingSignalInvMassPtBinSingle,Double_t * fMesonIntDeltaRangeFit, Int_t ptBin,Bool_t vary ,TString functionname)
{

	fFileErrLog<<"Start Fitting spectra with CB fit"<<endl;
	fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassRange[0],fMesonMassRange[1]);
	Double_t mesonAmplitude =fHistoMappingSignalInvMassPtBinSingle->GetMaximum();
	Double_t mesonAmplitudeMin = mesonAmplitude*50./100.;
	Double_t mesonAmplitudeMax = mesonAmplitude*115./100.;


	fFitReco = new TF1(functionname,CrystalBallBck,fMesonFitRange[0],fMesonFitRange[1],7);

	fFitGausExp = new TF1("crystal",CrystalBall,fMesonFitRange[0],fMesonFitRange[1],5);

	fFitLinearBck = new TF1("Linear","[0]+[1]*x",fMesonFitRange[0],fMesonFitRange[1]);

	fFitReco->SetParameter(0,mesonAmplitude);
	fFitReco->SetParameter(1,fMesonMassExpect);
	fFitReco->SetParameter(2,fMesonWidthExpect);
	fFitReco->SetParameter(3,2.);
	fFitReco->SetParameter(4,0.7);
	fFitReco->SetParameter(5,0.);
	fFitReco->SetParameter(6,1.);

	if (!vary) {
		fFitReco->FixParameter(3,fCBn);
		fFitReco->FixParameter(4,fCBAlpha);
	}

	fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
	fFitReco->SetParLimits(1,fMesonMassRange[0],fMesonMassRange[1]);
	fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);

	fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QR0");
	fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRE0");

	fFitReco->SetLineColor(3);
	fFitReco->SetLineWidth(1);
	fFitReco->SetLineStyle(1);

	if (vary) {
		fCBAlpha = fFitReco->GetParameter(4);
		fCBn = fFitReco->GetParError(3);
	}

	fFitGausExp->SetParameter(0,fFitReco->GetParameter(0));
	fFitGausExp->SetParameter(1,fFitReco->GetParameter(1));
	fFitGausExp->SetParameter(2,fFitReco->GetParameter(2));
	fFitGausExp->SetParameter(3,fFitReco->GetParameter(3));
	fFitGausExp->SetParameter(4,fFitReco->GetParameter(4));

	fFitGausExp->SetParError(0,fFitReco->GetParError(0));
	fFitGausExp->SetParError(1,fFitReco->GetParError(1));
	fFitGausExp->SetParError(2,fFitReco->GetParError(2));
	fFitGausExp->SetParError(3,fFitReco->GetParError(3));
	fFitGausExp->SetParError(4,fFitReco->GetParError(4));

	fFitLinearBck->SetParameter(0,fFitReco->GetParameter(5));
	fFitLinearBck->SetParameter(1,fFitReco->GetParameter(6));

	fFitLinearBck->SetParError(0,fFitReco->GetParError(5));
	fFitLinearBck->SetParError(1,fFitReco->GetParError(6));

	Int_t binCenterStart;
	Double_t startBinEdge;
	Int_t binCenterEnd;
	Double_t endBinEdge;

	TVirtualFitter * fitter = TVirtualFitter::GetFitter();

	if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 ){
		binCenterStart = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntDeltaRangeFit[0]-(fMesonMassExpect-fFitReco->GetParameter(1)));
		startBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterStart)- 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
		binCenterEnd = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntDeltaRangeFit[1]-(fMesonMassExpect-fFitReco->GetParameter(1)));
		endBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterEnd)+ 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);

		Int_t nFreePar = fFitReco->GetNumberFreeParameters();
		double * covMatrix = fitter->GetCovarianceMatrix();

		Float_t intLinearBack = fFitLinearBck->GetParameter(0)*(endBinEdge-startBinEdge)+
		0.5*fFitLinearBck->GetParameter(1)*(endBinEdge*endBinEdge-startBinEdge*startBinEdge);

		Float_t errorLinearBck = pow((pow( (endBinEdge-startBinEdge)*fFitReco->GetParError(5),2)+pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fFitReco->GetParError(6),2)+2*covMatrix[nFreePar*nFreePar-2]*(endBinEdge-startBinEdge)*0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5);

		fFileDataLog << "Parameter for bin " << ptBin << endl;
		fFileDataLog << "CrystalBall: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<< "\t "<< fFitReco->GetParameter(3) <<"+-" << fFitReco->GetParError(3)<< "\t "<< fFitReco->GetParameter(4) <<"+-" << fFitReco->GetParError(4)<<endl;
		fFileDataLog << "Linear: \t"<<fFitReco->GetParameter(5)<<"+-" << fFitReco->GetParError(5) << "\t "<<fFitReco->GetParameter(6) <<"+-" << fFitReco->GetParError(6)<< endl;

		fIntLinearBck = intLinearBack/fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
		fIntLinearBckError = errorLinearBck/fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
	} else {
		fFileErrLog << "Fitting failed in " << ptBin << " with status::" << gMinuit->fCstatu.Data() <<"why failed?"<<endl << endl;
	}
	fFitReco->DrawCopy("same");
}

void IntegrateHistoInvMass(TH1D * fHistoMappingSignalInvMassPtBinSingle, Double_t * fMesonIntDeltaRangeInt)
{
	Int_t binLowMassMeson = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntDeltaRangeInt[0]);
	Int_t binHighMassMeson = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntDeltaRangeInt[1]);
	fYields = fHistoMappingSignalInvMassPtBinSingle->IntegralAndError(binLowMassMeson,binHighMassMeson,fYieldsError);
}

void IntegrateHistoInvMassStream(TH1D * fHistoMappingSignalInvMassPtBinSingle, Double_t * fMesonIntDeltaRangeInt)
{
	Int_t binLowMassMeson = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntDeltaRangeInt[0]);
	Int_t binHighMassMeson = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntDeltaRangeInt[1]);
	fYields = fHistoMappingSignalInvMassPtBinSingle->IntegralAndError(binLowMassMeson,binHighMassMeson,fYieldsError);
	for ( Int_t M = binLowMassMeson; M < binHighMassMeson+1; M++){
		fFileDataLog << M << "\t" << fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(M) <<"\t" <<fHistoMappingSignalInvMassPtBinSingle->GetBinContent(M)<< "+-"<< fHistoMappingSignalInvMassPtBinSingle->GetBinError(M)<< endl;
	}
}


void IntegrateFitFunc(TF1 * fFunc, TH1D *  fHistoMappingSignalInvMassPtBinSingle,Double_t * fMesonIntDeltaRangeInt)
{

	fYieldsFunc = fFunc->Integral(fMesonIntDeltaRangeInt[0],fMesonIntDeltaRangeInt[1])/fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);

}

void FillHistosArrayMC(TH1D* fHistoMCMesonPtWithinAcceptanceFill, TH1D * fHistoMCMesonPtFill, TH1D * fDeltaPtFill)
{

   //Char_t nameHisto[100] = "fHistoMCMesonPtEtaWithinAcceptance";
   cout<<"Fallo aqui"<<endl;
   fHistoMCMesonWithinAccepPt = (TH1D*)fHistoMCMesonPtWithinAcceptanceFill->Rebin(fNBinsPt,"",fBinsPt); // Proper bins in Pt
   fHistoMCMesonWithinAccepPt->Divide(fDeltaPtFill);
   fHistoMCMesonPt1 = (TH1D*)fHistoMCMesonPtFill->Rebin(fNBinsPt,"",fBinsPt); // Proper bins in Pt
   fHistoMCMesonPt1->Divide(fDeltaPtFill);

}


void CalculateMesonAcceptance()
{

	fHistoMCMesonAcceptPt = new TH1D("fMCMesonAccepPt","",fNBinsPt,fBinsPt);
	fHistoMCMesonAcceptPt->Sumw2();
	fHistoMCMesonAcceptPt->Divide(fHistoMCMesonWithinAccepPt,fHistoMCMesonPt1,1.,1.,"B");
	fHistoMCMesonAcceptPt->DrawCopy();
	fFileDataLog << endl << "Calculation of the Acceptance" << endl;
	for ( Int_t i = 1; i < fHistoMCMesonAcceptPt->GetNbinsX()+1 ; i++){
		fFileDataLog << "Bin " << i << "\t"<< fHistoMCMesonAcceptPt->GetBinCenter(i)<< "\t" << fHistoMCMesonAcceptPt->GetBinContent(i) << "\t" << fHistoMCMesonAcceptPt->GetBinError(i) <<endl;
	}

}

void CalculateMesonEfficiency(TH1D* fMC_fMesonYieldsPt, TString nameEfi )
{
	fHistoMCMesonEffiPt = new TH1D(nameEfi.Data(),"",fNBinsPt,fBinsPt);
	fHistoMCMesonEffiPt->Sumw2();
	fHistoMCMesonEffiPt->Add(fMC_fMesonYieldsPt,1.);
	fHistoMCMesonEffiPt->Divide(fHistoMCMesonEffiPt,fHistoMCMesonWithinAccepPt,1.,1.,"B");
	fFileDataLog << endl << "Calculation of the Efficiency" << nameEfi.Data()<< endl;
	for ( Int_t i = 1; i < fHistoMCMesonEffiPt->GetNbinsX()+1 ; i++){
		fFileDataLog << "Bin " << i << "\t" << fHistoMCMesonEffiPt->GetBinCenter(i)<< "\t"<< fHistoMCMesonEffiPt->GetBinContent(i) << "\t" << fHistoMCMesonEffiPt->GetBinError(i) <<endl;
	}

}


void CalculateMesonEfficiencyWithoutGGCont(TH1D* fMC_fMesonYieldsPt, TH1D* fMC_fMesonYieldsTrueGGFrac, TString nameEfi, TString nameGGFrac)
{
        fHistoMCMesonEffiPt = new TH1D(nameEfi.Data(),"",fNBinsPt,fBinsPt);
        fHistoMCMesonEffiPt->Sumw2();
        fHistoMCMesonEffiPt->Add(fMC_fMesonYieldsPt,1.);
        fHistoMCRAWYieldTrueGGMeson = (TH1D*)fHistoMCMesonEffiPt->Clone(nameGGFrac.Data());
        fHistoMCRAWYieldTrueGGMeson->Multiply(fMC_fMesonYieldsTrueGGFrac);
        fHistoMCMesonEffiPt->Add(fHistoMCRAWYieldTrueGGMeson,-1.0); 
        fHistoMCMesonEffiPt->Divide(fHistoMCMesonEffiPt,fHistoMCMesonWithinAccepPt,1.,1.,"B");
        fFileDataLog << endl << "Calculation of the Efficiency" << nameEfi.Data()<< endl;
        for ( Int_t i = 1; i < fHistoMCMesonEffiPt->GetNbinsX()+1 ; i++){
                fFileDataLog << "Bin " << i << "\t" << fHistoMCMesonEffiPt->GetBinCenter(i)<< "\t"<< fHistoMCMesonEffiPt->GetBinContent(i) << "\t" << fHistoMCMesonEffiPt->GetBinError(i) <<endl;
        }

}

void SaveHistos(Int_t isMC, TString fCutID, TString fPrefix3)
{
	const char* nameOutput = Form("%s/%s/%s_%s_GammaConvV1DalitzWithoutCorrection%s_%s.root",fCutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPrefix3.Data(),fPeriodFlag.Data(),fCutID.Data());
	fOutput1 = new TFile(nameOutput,"RECREATE");
	
	fOutput1->WriteObject(fArrayDBinsPt,"fArrayDBinsPt");
	fOutput1->WriteObject(fArrayIParametersBins,"fArrayIParametersBins");
	fOutput1->WriteObject(fArrayDMesonMassRange,"fArrayDMesonMassRange");
	
	
	fHistoYieldMeson->Write();
	fHistoYieldMesonPerEvent->Write();
	fHistoSignMeson->Write();
	fHistoSBMeson->Write();

	fHistoYieldMesonNarrow->Write();
	fHistoYieldMesonPerEventNarrow->Write();
	fHistoSignMesonNarrow->Write();
	fHistoSBMesonNarrow->Write();

	fHistoYieldMesonWide->Write();
	fHistoYieldMesonPerEventWide->Write();
	fHistoSignMesonWide->Write();
	fHistoSBMesonWide->Write();

	fHistoMassMeson->Write();
	fHistoWidthMeson->Write();
	fHistoFWHMMeson->Write();
	fDeltaPt->Write();

	fHistoYieldMesonLeft->Write();
	fHistoYieldMesonLeftPerEvent->Write();
	fHistoSignMesonLeft->Write();
	fHistoSBMesonLeft->Write();

	fHistoYieldMesonLeftNarrow->Write();
	fHistoYieldMesonLeftPerEventNarrow->Write();
	fHistoSignMesonLeftNarrow->Write();
	fHistoSBMesonLeftNarrow->Write();

	fHistoYieldMesonLeftWide->Write();
	fHistoYieldMesonLeftPerEventWide->Write();
	fHistoSignMesonLeftWide->Write();
	fHistoSBMesonLeftWide->Write();

	fHistoMassMesonLeft->Write();
	fHistoWidthMesonLeft->Write();
	fHistoFWHMMesonLeft->Write();
	fMesonFullPtSignal->Write();
	fMesonFullPtBackground->Write();
	fMesonFullPtBackNorm->SetName("Mapping_BackNorm_InvMass_FullPt");
	fMesonFullPtBackNorm->Write();
	fNumberOfGoodESDTracksVtx->Write();
	fEventQuality->Write();
	if ( fESDEposEnegInvMassPt ) fESDEposEnegInvMassPt->Write();

	TString nameHistoSignal;
	TString nameHistoPeakPos;
	TString nameHistoBckNorm;
	TString fitnameSignal;
	TString nameHistoSignalPos;
	
	 for(Int_t ii =fStartPtBin;ii<fNBinsPt;ii++){
	   
	    if( fHistoWeightsBGZbinVsMbin[ii]    !=0x00) fHistoWeightsBGZbinVsMbin[ii]->Write(Form("BGWeights_%02d", ii));
	    if(fHistoFillPerEventBGZbinVsMbin[ii]!=0x00){
	    fHistoFillPerEventBGZbinVsMbin[ii]->Scale(1./fNEvents);
	    fHistoFillPerEventBGZbinVsMbin[ii]->Write(Form("BGPoolsFillstatus_%02d", ii));
	    }
   
	    fHistoMappingGGInvMassPtBin[ii]->Write();
	    nameHistoBckNorm = Form("Mapping_BckNorm_InvMass_in_Pt_Bin%02d", ii);
	    fHistoMappingBackNormInvMassPtBin[ii]->Write(nameHistoBckNorm.Data());
	    nameHistoSignal = Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d", ii);
	    fHistoMappingSignalInvMassPtBin[ii]->Write(nameHistoSignal.Data());
	    fitnameSignal = Form("Signal_InvMassFit_in_Pt_Bin%02d", ii);
	    if( fFitSignalInvMassPtBin[ii]      != 0x00 )fFitSignalInvMassPtBin[ii]->Write(fitnameSignal.Data());
	    if( fHistoMappingEEInvMassPtBin[ii] != 0x00 )fHistoMappingEEInvMassPtBin[ii]->Write();
	    
	  }

	
	
	if(isMC){
	  
		
		
		fHistoTrueSignMeson->Write();
		fHistoTrueSBMeson->Write();
		fHistoMCMesonPtWithinAcceptance->Write();
		fHistoMCMesonWithinAccepPt->Write(); // Proper bins in Pt
		fHistoMCMesonPt1->Write(); // Proper bins in Pt
		fHistoYieldTrueMeson->Write();
		fHistoYieldTrueMesonWide->Write();
		fHistoYieldTrueMesonNarrow->Write();
		fHistoTrueMesonInvMassVSPt->Write();
		fHistoYieldTrueGGMeson->Write();
		fHistoYieldTrueGGMesonWide->Write();
		fHistoYieldTrueGGMesonNarrow->Write();
		for(Int_t ii =fStartPtBin;ii<fNBinsPt;ii++){

			fHistoMappingTrueMesonInvMassPtBins[ii]->Write();
			if (fFitTrueSignalInvMassPtBin[ii]!=0x00) fFitTrueSignalInvMassPtBin[ii]->Write();
                        fHistoMappingTrueGGMesonInvMassPtBins[ii]->Write();
                        if ( fHistoMappingTrueGGBckInvMassPtBins[ii]!=0x00 ) fHistoMappingTrueGGBckInvMassPtBins[ii]->Write();
                        if ( fHistoMappingTrueContBckInvMassPtBins[ii]!= 0x00) fHistoMappingTrueContBckInvMassPtBins[ii]->Write();
			fHistoMappingSignalInvMassW0TruePi0PtBins[ii]->Write();
			fHistoMappingGGInvMassW0TruePi0PtBins[ii]->Write();
		}
	}
	fOutput1->Write();
	fOutput1->Close();
}

void SaveCorrectionHistos(TString fCutID, TString fPrefix3)
{
	const char* nameOutput = Form("%s/%s/%s_%s_GammaConvV1DalitzCorrectionHistos%s_%s.root",fCutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPrefix3.Data(),fPeriodFlag.Data(),fCutID.Data());
	fOutput2 = new TFile(nameOutput,"RECREATE");
	fHistoMCMesonAcceptPt->Write();
	fHistoTrueMesonInvMassVSPt->Write();
	fHistoMonteMesonEffiPt->Write();
	fHistoMonteMesonNarrowEffiPt->Write();
	fHistoMonteMesonWideEffiPt->Write();
	fHistoMonteMesonLeftEffiPt->Write();
	fHistoMonteMesonLeftNarrowEffiPt->Write();
	fHistoMonteMesonLeftWideEffiPt->Write();
	fHistoMCTrueMesonEffiPt->Write();
	fHistoMCTrueMesonNarrowEffiPt->Write();
	fHistoMCTrueMesonWideEffiPt->Write();

        //fHistoRAWYieldTrueGGFracMeson->Write();
        //fHistoRAWYieldTrueGGFracMesonWide->Write();
        //fHistoRAWYieldTrueGGFracMesonNarrow->Write();

	fHistoYieldTrueGGMeson->Write();
	fHistoYieldTrueGGMesonWide->Write();
	fHistoYieldTrueGGMesonNarrow->Write();
	fHistoYieldTrueGGFracMeson->Write();
	fHistoYieldTrueGGFracMesonWide->Write();
	fHistoYieldTrueGGFracMesonNarrow->Write();
	
	fHistoYieldTrueGGFracMesonForData->Write();
	fHistoYieldTrueGGFracMesonWideForData->Write();
	fHistoYieldTrueGGFracMesonNarrowForData->Write();

	fHistoTrueMassMeson->Write();
	fHistoTrueFWHMMeson->Write();
	fHistoMCMesonPt1->SetName("MC_Meson_genPt");
	fHistoMCMesonPt1->Write(); // Proper bins in Pt
	fHistoMCMesonPt->SetName("MC_Meson_genPt_oldBin");
	fHistoMCMesonPt->Scale(1./fHistoMCMesonPt->GetBinWidth(5));
	fHistoMCMesonPt->GetXaxis()->SetRangeUser(0.,25.);
	fHistoMCMesonPt->Write(); // Proper bins in Pt
	fHistoMCMesonDalitzPt->Write();
	fHistoMCMesonGGPt->Write();
	
	fOutput2->WriteObject(fArrayBRPi0Meson,"fArrayBRPi0Meson");
	fEventQuality->Write();
	fOutput2->Write();
	fOutput2->Close();
}

void Initialize(TString setPi0,Int_t numberOfBins)
{
	
	      if (setPi0.CompareTo("Pi0") == 0){
	  
		fNBinsPt = 		numberOfBins;
		fBinsPt= 		new Double_t[33];
		fNRebin = 		new Int_t[32];

		if ( fEnergyFlag.CompareTo("7TeV") == 0 ){
				fStartPtBin=	1;
				fColumn = 	5;
				fRow = 		5;

				if (fNBinsPt > 22) {
					cout << "You have chosen conference Plots and more than 22 bins, this is not possible, it will be reduced to 22 bins." << endl;
					fNBinsPt = 22;
				}
				for (Int_t i = 0; i < fNBinsPt+1; i++) {
					fBinsPt[i] = fBinsPi07TeVPt[i];
					if (i < fNBinsPt+1) fNRebin[i] = fBinsPi07TeVPtRebin[i];
				}
				fExampleBin = 5;
		}
                else if( fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0) {
                        
			fStartPtBin=        1;
                        fColumn =           4;
                        fRow =              3;
			
                        if (fNBinsPt > 8) {
                        cout << "You have chosen to have more than 8 bins, this is not possible, it will be reduced to 8" << endl;
                        fNBinsPt = 8;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
                        fBinsPt[i] = fBinsPi0HIPt[i];
			  if (i < fNBinsPt+1) fNRebin[i] = fBinsPi0HIPtRebin[i];
			}
            	        fExampleBin = 3;
                }
		else if( fEnergyFlag.CompareTo("pPb_5.023TeV") == 0 ){
		  
		  
		  
			        if( fMode == 6 ) {  //EMCAL Dalitz Binning

				fStartPtBin=	1;
				fColumn = 	5;
				fRow = 		5;
				
				if (fNBinsPt > 22) {
					cout << "You have chosen conference Plots and more than 22 bins, this is not possible, it will be reduced to 22 bins." << endl;
					fNBinsPt = 22;
				}
				for (Int_t i = 0; i < fNBinsPt+1; i++) {
					fBinsPt[i] = fBinsPi0EMCALDalitz5023GeVPt[i];
					if (i < fNBinsPt+1) fNRebin[i] = fBinsPi0EMCALDalitz5023GeVPtRebin[i];
				}
				fExampleBin = 7;



				} else {
		  
				fStartPtBin=	1;
				fColumn = 	6;
				fRow = 		4;
				
				if (fNBinsPt > 22) {
					cout << "You have chosen conference Plots and more than 22 bins, this is not possible, it will be reduced to 22 bins." << endl;
					fNBinsPt = 22;
				}
				for (Int_t i = 0; i < fNBinsPt+1; i++) {
					fBinsPt[i] = fBinsPi05023GeVPt[i];
					if (i < fNBinsPt+1) fNRebin[i] = fBinsPi05023GeVPtRebin[i];
				}
				fExampleBin = 5;
				
				}
				
		}
		else if (fEnergyFlag.CompareTo("2.76TeV") == 0) {
			fStartPtBin=	1;
			fColumn = 	3;
			fRow = 		3;

			if (fNBinsPt > 7) {
				cout << "You have chosen to have more than 10 bins, this is not possible, it will be reduced to 10" << endl;
				fNBinsPt = 7;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] = fBinsPi02760GeVPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsPi02760GeVPtRebin[i];
			}
			fExampleBin = 3;
		} else if (fEnergyFlag.CompareTo("900GeV") == 0) {
			fStartPtBin=	1;
			fColumn = 	4;
			fRow = 		3;

			if (fNBinsPt > 11) {
				cout << "You have chosen to have more than 11 bins, this is not possible, it will be reduced to 11" << endl;
				fNBinsPt = 11;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] = fBinsPi0900GeVPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsPi0900GeVPtRebin[i];
			}
			fExampleBin = 4;
		} else {
			fStartPtBin=	1;
			fColumn = 	5;
			fRow = 		4;

			if (fNBinsPt > 20) {
				cout << "You have chosen to have more than 20 bins, this is not possible, it will be reduced to 20" << endl;
				fNBinsPt = 20;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] = fBinsPi0HIPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsPi0HIPtRebin[i];
			}

			fExampleBin = 4;
		}
		fBGFitRange = 		new Double_t[2]; 	fBGFitRange[0]=	0.17; 	fBGFitRange[1]=	0.3; //eta 0.9
		fBGFitRangeLeft = 	new Double_t[2]; 	fBGFitRangeLeft[0]=	0.05; 	fBGFitRangeLeft[1]=	0.08;  // eta 09
		fMesonPlotRange = 	new Double_t[2]; 	fMesonPlotRange[0]=	0.13; 	fMesonPlotRange[1]=	0.138;
		
		
		
		if( fMode == 6 || fMode == 7 ) {
		  
		  cout<<"Entro al Mode "<<fMode<<" Setting fits ranges"<<endl;
		
		//fMesonIntDeltaRange       = new Double_t[2]; fMesonIntDeltaRange[0] = -0.045; 	fMesonIntDeltaRange[1]=	0.160; //NOTE *****
		//fMesonIntDeltaRangeWide   = new Double_t[2]; fMesonIntDeltaRangeWide[0]=0.075; 	fMesonIntDeltaRangeWide[1]=0.170;
		//fMesonIntDeltaRangeNarrow = new Double_t[2]; fMesonIntDeltaRangeNarrow[0]=0.1;  fMesonIntDeltaRangeNarrow[1]=0.155;
		  
		  fMesonIntDeltaRange       = new Double_t[2]; fMesonIntDeltaRange[0] = -0.045; 	fMesonIntDeltaRange[1]=	0.025; //NOTE *****
		  fMesonIntDeltaRangeWide   = new Double_t[2]; fMesonIntDeltaRangeWide[0]=-0.06; 	fMesonIntDeltaRangeWide[1]=0.035;
		  fMesonIntDeltaRangeNarrow = new Double_t[2]; fMesonIntDeltaRangeNarrow[0]=-0.035;     fMesonIntDeltaRangeNarrow[1]=0.02;
		
		
		} else {
		
		
		//fMesonIntDeltaRange 	    = 	  new Double_t[2]; 	fMesonIntDeltaRange[0] = 0.1; 	    fMesonIntDeltaRange[1]=	0.145;
		//fMesonIntDeltaRangeWide   =     new Double_t[2]; 	fMesonIntDeltaRangeWide[0]=0.08;    fMesonIntDeltaRangeWide[1]=0.160;
		//fMesonIntDeltaRangeNarrow =     new Double_t[2];        fMesonIntDeltaRangeNarrow[0]=0.12;  fMesonIntDeltaRangeNarrow[1]=0.14;
		
		/*
		fMesonIntDeltaRange 	  = 	new Double_t[2]; 	fMesonIntDeltaRange[0] = -0.035;    fMesonIntDeltaRange[1]=	0.01;
		fMesonIntDeltaRangeWide   =     new Double_t[2]; 	fMesonIntDeltaRangeWide[0]=-0.055;    fMesonIntDeltaRangeWide[1]=0.025;
		fMesonIntDeltaRangeNarrow =     new Double_t[2];        fMesonIntDeltaRangeNarrow[0]=-0.015;  fMesonIntDeltaRangeNarrow[1]=0.005;
		*/
		
		fMesonIntDeltaRange 	  = 	new Double_t[2]; 	fMesonIntDeltaRange[0] = -0.035;    fMesonIntDeltaRange[1]=	0.015; //0.015 Old
		fMesonIntDeltaRangeWide   =     new Double_t[2]; 	fMesonIntDeltaRangeWide[0]=-0.055;    fMesonIntDeltaRangeWide[1]=0.025;
		fMesonIntDeltaRangeNarrow =     new Double_t[2];        fMesonIntDeltaRangeNarrow[0]=-0.015;  fMesonIntDeltaRangeNarrow[1]=0.010; //0.010 old
		
		
		
		}
		
		
		
		fMesonMassRange = 	new Double_t[2]; 	fMesonMassRange[0]=	0.; 		fMesonMassRange[1]=	0.3;
		if( fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0) { 
			fMesonFitRange = new Double_t[2]; 	fMesonFitRange[0]=	0.05; 	fMesonFitRange[1]=	0.2; 
		} else { 
		  
			if( fMode == 6 || fMode == 7 ){
			  
			  if( fEnergyFlag.CompareTo("2.76TeV") == 0 ){
			  
			  fMesonFitRange = new Double_t[2]; 	fMesonFitRange[0]=	0.05; 	fMesonFitRange[1]=	0.25;
			  
			  } else if ( fEnergyFlag.CompareTo("pPb_5.023TeV") == 0 ){
			    
			      fMesonFitRange = new Double_t[2]; 	fMesonFitRange[0]=	0.05; 	fMesonFitRange[1]=	0.20;  
			    
			  }
			  
			} else {
		    
			  fMesonFitRange = new Double_t[2]; 	fMesonFitRange[0]=	0.07; 	fMesonFitRange[1]=	0.20; //0.25 //0.17
			
			}
		}
		//		fMesonFitRangeWithoutPeak = new Double_t[2]; fMesonFitRangeWithoutPeak[0] = 0.05; fMesonFitRangeWithoutPeak[0] = 0.3;
		fMesonId=			111;
		fMesonWidthExpect = 0.003;
		fMesonLambdaTail = 	0.012;
		fMesonWidthRange = 	new Double_t[2]; 	fMesonWidthRange[0]=0.001; 	fMesonWidthRange[1]=0.009;
		fMesonLambdaTailRange = new Double_t[2]; fMesonLambdaTailRange[0]=0.001; fMesonLambdaTailRange[1]=0.02;
		fMidPt = 							new Double_t[2]; fMidPt[0] = 0.8; fMidPt[1] = 2.5;
	} else if (setPi0.CompareTo("Eta") == 0){
		fNBinsPt = 		numberOfBins;
		fBinsPt= 		new Double_t[15];
		fNRebin = 		new Int_t[14];

		if ( fEnergyFlag.CompareTo("7TeV") == 0 ) {
			fStartPtBin=	1;
			fColumn = 	4;
			fRow = 		3;

			if (fNBinsPt > 9) {
				cout << "You have chosen to have more than 9 bins for Eta, this is not possible, it will be reduced to 9" << endl;
				fNBinsPt = 9;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] = fBinsEta7TeVPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsEta7TeVPtRebin[i];
			}
			fExampleBin = 4;
		} 
		else if( fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0) { 
                        
			fStartPtBin=        1;
                        fColumn =           4;
                        fRow =              2;
			
                        if (fNBinsPt > 4) {
                        cout << "You have chosen to have more than 4 bins, this is not possible, it will be reduced to 4" << endl;
                        fNBinsPt = 4;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
                        fBinsPt[i] = fBinsEtaHIPt[i];
			  if (i < fNBinsPt+1) fNRebin[i] = fBinsEtaHIPtRebin[i];
			}
            	        fExampleBin = 2;      
                }
		else if( fEnergyFlag.CompareTo("pPb_5.023TeV") == 0 ) {
		 	fStartPtBin=	1;
			fColumn = 	4;
			fRow = 		3;

			if (fNBinsPt > 9) {
				cout << "You have chosen to have more than 9 bins for Eta, this is not possible, it will be reduced to 9" << endl;
				fNBinsPt = 9;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] = fBinsEta5023GeVPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsEta5023GeVPtRebin[i];
			}
			fExampleBin = 4;
		}  
		else if ( fEnergyFlag.CompareTo("2.76TeV") == 0 ) {
			fStartPtBin=	1;
			fColumn = 	3;
			fRow = 		3;

			if (fNBinsPt > 7) {
				cout << "You have chosen to have more than 9 bins for Eta, this is not possible, it will be reduced to 9" << endl;
				fNBinsPt = 7;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] = fBinsEta2760GeVPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsEta2760GeVPtRebin[i];
			}
			fExampleBin = 4;
		} else if (fEnergyFlag.CompareTo("900GeV") == 0) {
			fStartPtBin=	1;
			fColumn = 	2;
			fRow = 		2;

			if (fNBinsPt > 3) {
				cout << "You have chosen to have more than 3 bins for Eta, this is not possible, it will be reduced to 3" << endl;
				fNBinsPt = 3;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] = fBinsEta900GeVPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsEta900GeVPtRebin[i];
			}
			fExampleBin = 2;
		} else {
			fStartPtBin=	1;
			fColumn = 	2;
			fRow = 		2;

			if (fNBinsPt > 4) {
				cout << "You have chosen to have more than 4 bins, this is not possible, it will be reduced to 4" << endl;
				fNBinsPt = 4;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] = fBinsEtaHIPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsEtaHIPtRebin[i];
			}

			fExampleBin = 2;
		}
		fBGFitRange = 		new Double_t[2]; 	fBGFitRange[0]=	0.58; 	fBGFitRange[1]=	0.8; //eta 0.9
		fBGFitRangeLeft = 	new Double_t[2]; 	fBGFitRangeLeft[0]=	0.35; 	fBGFitRangeLeft[1]=	0.48;  // eta 09
		fMesonPlotRange = 	new Double_t[2]; 	fMesonPlotRange[0]=	0.53; 	fMesonPlotRange[1]=	0.560;
		
		fMesonIntDeltaRange = 	new Double_t[2]; 	fMesonIntDeltaRange[0] = -0.048;     fMesonIntDeltaRange[1]=	0.022;
		fMesonIntDeltaRangeWide = new Double_t[2]; 	fMesonIntDeltaRangeWide[0]=-0.68;     fMesonIntDeltaRangeWide[1]=0.032;
		fMesonIntDeltaRangeNarrow = new Double_t[2];    fMesonIntDeltaRangeNarrow[0]=-0.033;  fMesonIntDeltaRangeNarrow[1]=0.012;
		
		
		fMesonMassRange = 	new Double_t[2]; 	fMesonMassRange[0]=	0.35; 	fMesonMassRange[1]=	0.8;;
		fMesonFitRange = 	new Double_t[2]; 	fMesonFitRange[0]=	0.4; 	fMesonFitRange[1]=	0.7;
		//		fMesonFitRangeWithoutPeak = new Double_t[2]; fMesonFitRangeWithoutPeak[0] = 0.4; fMesonFitRangeWithoutPeak[0] = 0.7;
		fMesonId=			221;
		fMesonWidthExpect = 0.005;
		fMesonLambdaTail = 	0.007;
		fMesonWidthRange = 	new Double_t[2]; 	 fMesonWidthRange[0]=0.002;     fMesonWidthRange[1]=0.010;
		fMesonLambdaTailRange = new Double_t[2]; fMesonLambdaTailRange[0]=0.0005; fMesonLambdaTailRange[1]=0.026;
		fMidPt = 							new Double_t[2]; fMidPt[0] = 1.2; fMidPt[1] = 2.5;
	} else {
		fNBinsPt = 		numberOfBins;
		fBinsPt= 		new Double_t[15];
		fNRebin = 		new Int_t[14];

		if ( fEnergyFlag.CompareTo("7TeV") == 0 ) { 
			fStartPtBin=	1;
			fColumn = 	4;
			fRow = 		3;

			if (fNBinsPt > 9) {
				cout << "You have chosen to have more than 9 bins for Eta, this is not possible, it will be reduced to 9" << endl;
				fNBinsPt = 9;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] = fBinsEta7TeVPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsPi0EtaBinning7TeVPtRebin[i];
			}
			fExampleBin = 3;
		}
		else if( fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0) { 
                        
			fStartPtBin=        1;
                        fColumn =           4;
                        fRow =              2;
			
                        if (fNBinsPt > 4) {
                        cout << "You have chosen to have more than 4 bins, this is not possible, it will be reduced to 4" << endl;
                        fNBinsPt = 4;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
                        fBinsPt[i] = fBinsEtaHIPt[i];
			  if (i < fNBinsPt+1) fNRebin[i] = fBinsPi0EtaBinningHIPtRebin[i];
			}
            	        fExampleBin = 2;      
                }
		else if ( fEnergyFlag.CompareTo("pPb_5.023TeV") == 0 ) {
		        fStartPtBin=	1;
			fColumn = 	4;
			fRow = 		3;

			if (fNBinsPt > 9) {
				cout << "You have chosen to have more than 9 bins for Eta, this is not possible, it will be reduced to 9" << endl;
				fNBinsPt = 9;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] = fBinsEta5023GeVPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsPi0EtaBinning5023GeVPtRebin[i];
			}
			fExampleBin = 4;    
		}
		else if (fEnergyFlag.CompareTo("2.76TeV") == 0) {
			fStartPtBin=	1;
			fColumn = 	3;
			fRow = 		3;

			if (fNBinsPt > 7) {
				cout << "You have chosen to have more than 9 bins for Eta, this is not possible, it will be reduced to 9" << endl;
				fNBinsPt = 7;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] = fBinsEta2760GeVPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsPi0EtaBinning2760GeVPtRebin[i];
			}
			fExampleBin = 3;
		} else if (fEnergyFlag.CompareTo("900GeV") == 0) {
			fStartPtBin=	1;
			fColumn = 	2;
			fRow = 		2;

			if (fNBinsPt > 3) {
				cout << "You have chosen to have more than 3 bins for Eta, this is not possible, it will be reduced to 3" << endl;
				fNBinsPt = 3;
			}
			for (Int_t i = 0; i < fNBinsPt+1; i++) {
				fBinsPt[i] = fBinsEta900GeVPt[i];
				if (i < fNBinsPt+1) fNRebin[i] = fBinsPi0EtaBinning900GeVPtRebin[i];
			}
			fExampleBin = 2;
		}
		fBGFitRange = 		new Double_t[2]; 	fBGFitRange[0]=	0.17; 	fBGFitRange[1]=	0.3; //eta 0.9
		fBGFitRangeLeft = 	new Double_t[2]; 	fBGFitRangeLeft[0]=	0.05; 	fBGFitRangeLeft[1]=	0.08;  // eta 09
		fMesonPlotRange = 	new Double_t[2]; 	fMesonPlotRange[0]=	0.13; 	fMesonPlotRange[1]=	0.138;
		fMesonIntDeltaRange = 	new Double_t[2]; 	fMesonIntDeltaRange[0] = 0.1; 	fMesonIntDeltaRange[1]=	0.145;
		fMesonIntDeltaRangeWide = new Double_t[2]; 	fMesonIntDeltaRangeWide[0]=0.08; 	fMesonIntDeltaRangeWide[1]=0.160;
		fMesonIntDeltaRangeNarrow = new Double_t[2]; fMesonIntDeltaRangeNarrow[0]=0.12; fMesonIntDeltaRangeNarrow[1]=0.14;
		fMesonMassRange = 	new Double_t[2]; 	fMesonMassRange[0]=	0.; 		fMesonMassRange[1]=	0.3;
		if( fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0) { 
			fMesonFitRange = 	new Double_t[2]; 	fMesonFitRange[0]=	0.05; 	fMesonFitRange[1]=	0.2; 
		}
		else { 
			fMesonFitRange = new Double_t[2]; 	fMesonFitRange[0]=	0.07; 	fMesonFitRange[1]=	0.20; //0.25 //0.17
		}
			//		fMesonFitRangeWithoutPeak = new Double_t[2]; fMesonFitRangeWithoutPeak[0] = 0.05; fMesonFitRangeWithoutPeak[0] = 0.3;
		fMesonId=			111;
		fMesonWidthExpect = 0.003;
		fMesonLambdaTail = 	0.012;
		fMesonWidthRange = 	new Double_t[2]; 	fMesonWidthRange[0]=0.001; 	fMesonWidthRange[1]=0.009;
		fMesonLambdaTailRange = new Double_t[2]; fMesonLambdaTailRange[0]=0.001; fMesonLambdaTailRange[1]=0.02;
		fMidPt = 			new Double_t[2]; fMidPt[0] = 0.8; fMidPt[1] = 2.5;
	}
	
	fArrayDBinsPt		=  new TArrayD(fNBinsPt+1,fBinsPt);
	fArrayDMesonMassRange   =  new TArrayD(2,fMesonMassRange);
	fArrayIParametersBins 	=  new TArrayI(4);
	fArrayIParametersBins->AddAt(fStartPtBin,0);
	fArrayIParametersBins->AddAt(fNBinsPt,1);
	fArrayIParametersBins->AddAt(fColumn,2);
	fArrayIParametersBins->AddAt(fRow,3);
	
	fArrayBRPi0Meson	= new TArrayD(4); //0: Pi0->GG BR DPG, 1: Pi0->eeG BR DPG, 2: Pi0->GG BR MC , 3: Pi0->eeG BR MC

	fFullPt = 						new Double_t[2]; fFullPt[0] = 0.4; fFullPt[1] = 15.;
	fMesonCurIntRange = 					new Double_t[2];
	fMesonCurIntRangeWide = 				new Double_t[2];
	fMesonCurIntRangeNarrow = 				new Double_t[2];
	fMesonCurLeftIntRange = 				new Double_t[2];
	fMesonCurLeftIntRangeWide = 				new Double_t[2];
	fMesonCurLeftIntRangeNarrow = 				new Double_t[2];
	fMesonTrueIntRange = 					new Double_t[2];
	fMesonTrueIntRangeWide = 				new Double_t[2];
	fMesonTrueIntRangeNarrow = 				new Double_t[2];

	fGGYields = 						new Double_t[fNBinsPt];
	fBckYields = 						new Double_t[fNBinsPt];
	fTotalBckYields = 					new Double_t[fNBinsPt];
	fMesonYields = 						new Double_t[fNBinsPt];
	fMesonTrueYields = 					new Double_t[fNBinsPt];
	fMesonTrueGGYields = 					new Double_t[fNBinsPt];
	fMesonYieldsFunc = 					new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFunc = 				new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFunc = 			new Double_t[fNBinsPt];
	fMesonYieldsPerEvent = 					new Double_t[fNBinsPt];
	fMesonMass = 						new Double_t[fNBinsPt];
	fMesonWidth = 						new Double_t[fNBinsPt];
	fMesonSB = 						new Double_t[fNBinsPt];
	fMesonSign = 						new Double_t[fNBinsPt];
	fMesonFWHM = 						new Double_t[fNBinsPt];
	fMesonFWHMAlpha01 =					new Double_t[fNBinsPt];
	fMesonTrueMass = 					new Double_t[fNBinsPt];
	fMesonTrueFWHM = 					new Double_t[fNBinsPt];
	fMesonTrueSB = 						new Double_t[fNBinsPt];
	fMesonTrueSign = 					new Double_t[fNBinsPt];


	// Normalization at the left of the peak
	fGGYieldsLeft = 					new Double_t[fNBinsPt];
	fBckYieldsLeft = 					new Double_t[fNBinsPt];
	fTotalBckYieldsLeft = 				new Double_t[fNBinsPt];
	fMesonYieldsLeft = 					new Double_t[fNBinsPt];
	fMesonYieldsFuncLeft = 				new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncLeft = 		new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncLeft = 	new Double_t[fNBinsPt];
	fMesonYieldsLeftPerEvent = 			new Double_t[fNBinsPt];
	fMesonMassLeft = 					new Double_t[fNBinsPt];
	fMesonWidthLeft = 					new Double_t[fNBinsPt];
	fMesonSBLeft = 					new Double_t[fNBinsPt];
	fMesonSignLeft = 					new Double_t[fNBinsPt];
	fMesonFWHMLeft = 					new Double_t[fNBinsPt];

	// Narrow Integration Window
	fGGYieldsNarrow = 					new Double_t[fNBinsPt];
	fBckYieldsNarrow = 					new Double_t[fNBinsPt];
	fTotalBckYieldsNarrow =	 			new Double_t[fNBinsPt];
	fMesonYieldsNarrow = 				new Double_t[fNBinsPt];
	fMesonTrueYieldsNarrow = 			new Double_t[fNBinsPt];
	fMesonTrueGGYieldsNarrow = 					new Double_t[fNBinsPt];
	fMesonYieldsFuncNarrow = 			new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncNarrow = 	new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncNarrow = 	new Double_t[fNBinsPt];
	fMesonYieldsPerEventNarrow = 			new Double_t[fNBinsPt];
	fMesonSBNarrow = 					new Double_t[fNBinsPt];
	fMesonSignNarrow = 					new Double_t[fNBinsPt];

	fGGYieldsLeftNarrow = 				new Double_t[fNBinsPt];
	fBckYieldsLeftNarrow = 				new Double_t[fNBinsPt];
	fTotalBckYieldsLeftNarrow = 			new Double_t[fNBinsPt];
	fMesonYieldsLeftNarrow = 			new Double_t[fNBinsPt];
	fMesonYieldsFuncLeftNarrow = 			new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncLeftNarrow = new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncLeftNarrow = new Double_t[fNBinsPt];
	fMesonYieldsLeftPerEventNarrow = 		new Double_t[fNBinsPt];
	fMesonSBLeftNarrow = 				new Double_t[fNBinsPt];
	fMesonSignLeftNarrow = 				new Double_t[fNBinsPt];

	// Wide Integration Window
	fGGYieldsWide = 					new Double_t[fNBinsPt];
	fBckYieldsWide = 					new Double_t[fNBinsPt];
	fTotalBckYieldsWide =	 			new Double_t[fNBinsPt];
	fMesonYieldsWide = 					new Double_t[fNBinsPt];
	fMesonTrueYieldsWide = 				new Double_t[fNBinsPt];
	fMesonTrueGGYieldsWide = 					new Double_t[fNBinsPt];
	fMesonYieldsFuncWide = 				new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncWide = 		new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncWide =	new Double_t[fNBinsPt];
	fMesonYieldsPerEventWide = 			new Double_t[fNBinsPt];
	fMesonSBWide = 					new Double_t[fNBinsPt];
	fMesonSignWide = 					new Double_t[fNBinsPt];

	fGGYieldsLeftWide = 				new Double_t[fNBinsPt];
	fBckYieldsLeftWide = 				new Double_t[fNBinsPt];
	fTotalBckYieldsLeftWide = 			new Double_t[fNBinsPt];
	fMesonYieldsLeftWide = 				new Double_t[fNBinsPt];
	fMesonYieldsFuncLeftWide = 			new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncLeftWide =	new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncLeftWide = new Double_t[fNBinsPt];
	fMesonYieldsLeftPerEventWide = 		new Double_t[fNBinsPt];
	fMesonSBLeftWide = 					new Double_t[fNBinsPt];
	fMesonSignLeftWide = 				new Double_t[fNBinsPt];

	fGGYieldsError = 					new Double_t[fNBinsPt];
	fBckYieldsError = 					new Double_t[fNBinsPt];
	fTotalBckYieldsError = 				new Double_t[fNBinsPt];
	fMesonYieldsError = 				new Double_t[fNBinsPt];
	fMesonTrueYieldsError = 				new Double_t[fNBinsPt];
	fMesonTrueGGYieldsError = 					new Double_t[fNBinsPt];
	fMesonYieldsFuncError = 				new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncError = 	new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncError =	new Double_t[fNBinsPt];
	fMesonYieldsPerEventError = 			new Double_t[fNBinsPt];
	fMesonMassError = 					new Double_t[fNBinsPt];
	fMesonWidthError = 					new Double_t[fNBinsPt];
	fMesonSBError = 					new Double_t[fNBinsPt];
	fMesonSignError = 					new Double_t[fNBinsPt];
	fMesonFWHMError = 					new Double_t[fNBinsPt];
	fMesonFWHMAlpha01Error = 				new Double_t[fNBinsPt];
	fMesonTrueMassError = 				new Double_t[fNBinsPt];
	fMesonTrueFWHMError = 				new Double_t[fNBinsPt];
	fMesonTrueSBError = 				new Double_t[fNBinsPt];
	fMesonTrueSignError = 				new Double_t[fNBinsPt];

	fGGYieldsLeftError = 				new Double_t[fNBinsPt];
	fBckYieldsLeftError =	 			new Double_t[fNBinsPt];
	fTotalBckYieldsLeftError = 			new Double_t[fNBinsPt];
	fMesonYieldsLeftError = 				new Double_t[fNBinsPt];
	fMesonYieldsFuncLeftError = 			new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncLeftError = 	new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncLeftError = new Double_t[fNBinsPt];
	fMesonYieldsLeftPerEventError = 		new Double_t[fNBinsPt];
	fMesonMassLeftError = 				new Double_t[fNBinsPt];
	fMesonWidthLeftError = 				new Double_t[fNBinsPt];
	fMesonSBLeftError = 				new Double_t[fNBinsPt];
	fMesonSignLeftError = 				new Double_t[fNBinsPt];
	fMesonFWHMLeftError = 				new Double_t[fNBinsPt];

	// Narrow integration Window
	fGGYieldsNarrowError = 				new Double_t[fNBinsPt];
	fBckYieldsNarrowError = 			new Double_t[fNBinsPt];
	fTotalBckYieldsNarrowError = 			new Double_t[fNBinsPt];
	fMesonYieldsNarrowError = 			new Double_t[fNBinsPt];
	fMesonTrueYieldsNarrowError = 			new Double_t[fNBinsPt];
	fMesonTrueGGYieldsNarrowError = 		new Double_t[fNBinsPt];
	fMesonYieldsFuncNarrowError = 			new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncNarrowError 	= 	new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncNarrowError = 	new Double_t[fNBinsPt];
	fMesonYieldsPerEventNarrowError = 		new Double_t[fNBinsPt];
	fMesonSBNarrowError = 				new Double_t[fNBinsPt];
	fMesonSignNarrowError = 			new Double_t[fNBinsPt];

	fGGYieldsLeftNarrowError = 			new Double_t[fNBinsPt];
	fBckYieldsLeftNarrowError = 			new Double_t[fNBinsPt];
	fTotalBckYieldsLeftNarrowError = 		new Double_t[fNBinsPt];
	fMesonYieldsLeftNarrowError = 		new Double_t[fNBinsPt];
	fMesonYieldsFuncLeftNarrowError = 		new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncLeftNarrowError = new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncLeftNarrowError = new Double_t[fNBinsPt];
	fMesonYieldsLeftPerEventNarrowError = 	new Double_t[fNBinsPt];
	fMesonSBLeftNarrowError = 			new Double_t[fNBinsPt];
	fMesonSignLeftNarrowError = 			new Double_t[fNBinsPt];

	// Wide integration Window
	fGGYieldsWideError = 				new Double_t[fNBinsPt];
	fBckYieldsWideError = 				new Double_t[fNBinsPt];
	fTotalBckYieldsWideError = 			new Double_t[fNBinsPt];
	fMesonYieldsWideError = 				new Double_t[fNBinsPt];
	fMesonTrueYieldsWideError = 			new Double_t[fNBinsPt];
	fMesonTrueGGYieldsWideError = 					new Double_t[fNBinsPt];
	fMesonYieldsFuncWideError = 			new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncWideError = 	new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncWideError = new Double_t[fNBinsPt];
	fMesonYieldsPerEventWideError = 		new Double_t[fNBinsPt];
	fMesonSBWideError = 				new Double_t[fNBinsPt];
	fMesonSignWideError = 				new Double_t[fNBinsPt];

	fGGYieldsLeftWideError = 			new Double_t[fNBinsPt];
	fBckYieldsLeftWideError = 			new Double_t[fNBinsPt];
	fTotalBckYieldsLeftWideError = 		new Double_t[fNBinsPt];
	fMesonYieldsLeftWideError = 			new Double_t[fNBinsPt];
	fMesonYieldsFuncLeftWideError = 		new Double_t[fNBinsPt];
	fMesonYieldsResidualBckFuncLeftWideError = new Double_t[fNBinsPt];
	fMesonYieldsCorResidualBckFuncLeftWideError = new Double_t[fNBinsPt];
	fMesonYieldsLeftPerEventWideError = 	new Double_t[fNBinsPt];
	fMesonSBLeftWideError = 				new Double_t[fNBinsPt];
	fMesonSignLeftWideError = 			new Double_t[fNBinsPt];

	fHistoMappingTrueMesonInvMassPtBins 	= 	new TH1D*[fNBinsPt];    
	fHistoMappingTrueGGMesonInvMassPtBins 	= 	new TH1D*[fNBinsPt];
        fHistoMappingTrueGGBckInvMassPtBins     =       new TH1D*[fNBinsPt];
        fHistoMappingTrueContBckInvMassPtBins   =       new TH1D*[fNBinsPt];
	fHistoMappingSignalInvMassW0TruePi0PtBins =     new TH1D*[fNBinsPt];
	fHistoMappingGGInvMassW0TruePi0PtBins 	=     	new TH1D*[fNBinsPt];
	
	fHistoWeightsBGZbinVsMbin 		= new TH2F*[fNBinsPt];    
	fHistoFillPerEventBGZbinVsMbin 		= new TH2F*[fNBinsPt];    
	
	fHistoMappingGGInvMassPtBin = 		new TH1D*[fNBinsPt];    
	fHistoMappingEEInvMassPtBin =           new TH1D*[fNBinsPt];
	fHistoMappingBackInvMassPtBin = 		new TH1D*[fNBinsPt];
	fHistoMappingBackNormInvMassPtBin = 	new TH1D*[fNBinsPt];
	fHistoMappingSignalInvMassPtBin =		new TH1D*[fNBinsPt];
	fHistoMappingRatioSBInvMassPtBin= 		new TH1D*[fNBinsPt];

	fFitSignalInvMassPtBin = 			new TF1*[fNBinsPt];
	fFitTrueSignalInvMassPtBin = 			new TF1*[fNBinsPt];
	fFitSignalPeakPosInvMassPtBin = 		new TF1*[fNBinsPt];
	fFitBckInvMassPtBin = 				new TF1*[fNBinsPt];
	fFitRatioInvMassPtBin = 				new TF1*[fNBinsPt];
	// Histograms for normalization on the left of the peak
	fHistoMappingBackNormInvMassLeftPtBin = new TH1D*[fNBinsPt];
	fHistoMappingSignalInvMassLeftPtBin = 	new TH1D*[fNBinsPt];
	fHistoMappingPeakPosInvMassPtBin = 	new TH1D*[fNBinsPt];

	fFitInvMassLeftPtBin = 				new TF1*[fNBinsPt];
	fFitSignalPeakPosInvMassLeftPtBin = 	new TF1*[fNBinsPt];
	fFitBckInvMassLeftPtBin = 			new TF1*[fNBinsPt];
	fFitPeakPosPtBin = 					new TF1*[fNBinsPt];
	fMesonMassPeakPos =					new Double_t[fNBinsPt];
	fMesonMassPeakPosError = 			new Double_t[fNBinsPt];

	for(Int_t i = 0;i<fNBinsPt; i++){
	   fHistoMappingTrueMesonInvMassPtBins[i]       = NULL;	// array of histos for pt slices
	   fHistoMappingTrueGGMesonInvMassPtBins[i]     = NULL;
           fHistoMappingTrueGGBckInvMassPtBins[i]       = NULL;
           fHistoMappingTrueContBckInvMassPtBins[i]     = NULL;
	   fHistoMappingSignalInvMassW0TruePi0PtBins[i] = NULL;
		
	   fHistoMappingGGInvMassPtBin[i] = NULL;
	   fHistoMappingEEInvMassPtBin[i] = NULL;
	   fHistoMappingBackInvMassPtBin[i] = NULL;
	   fHistoMappingBackNormInvMassPtBin[i] = NULL;
	   fHistoMappingSignalInvMassPtBin[i] = NULL;
	   fHistoMappingRatioSBInvMassPtBin[i] = NULL;
	   
	   fFitSignalInvMassPtBin[i] = NULL;
	   fFitTrueSignalInvMassPtBin[i] = NULL;
	   fFitSignalPeakPosInvMassPtBin[i] = NULL;
	   fFitBckInvMassPtBin[i] = NULL;
	   fFitRatioInvMassPtBin[i] = NULL;
	// Histograms for normalization on the left of the peak
	   fHistoMappingBackNormInvMassLeftPtBin[i] = NULL;
	   fHistoMappingSignalInvMassLeftPtBin[i] = NULL;
	   fHistoMappingPeakPosInvMassPtBin[i] = NULL;
	   
	   fFitInvMassLeftPtBin[i] = NULL;
	   fFitSignalPeakPosInvMassLeftPtBin[i] = NULL;
	   fFitBckInvMassLeftPtBin[i] = NULL;
	   fFitPeakPosPtBin[i] = NULL;
	}
}

void CalculateFWHM(TF1 * fFunc)
{
	// Default function
	TF1* fFunc_def;
	fFunc_def = new TF1("fFunc_def","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);
	fFunc_def->SetParameter(0,fFunc->GetParameter(0));
	fFunc_def->SetParameter(1,fFunc->GetParameter(1));
	fFunc_def->SetParameter(2,fFunc->GetParameter(2));
	fFunc_def->SetParameter(3,fFunc->GetParameter(3));


	//FWHM
	fFWHMFunc = fFunc_def->GetX(fFunc_def->GetParameter(0)*0.5,fFunc_def->GetParameter(1), fMesonFitRange[1]) - fFunc_def->GetX(fFunc_def->GetParameter(0)*0.5,fMesonFitRange[0],fFunc_def->GetParameter(1));

	//FWHM error +
	TF1* fFunc_plus;
	//	fFunc_plus = fFunc;
	fFunc_plus = new TF1("fFunc_plus","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);

	fFunc_plus->SetParameter(0,fFunc->GetParameter(0) + fFunc->GetParError(0));
	fFunc_plus->SetParameter(1,fFunc->GetParameter(1) + fFunc->GetParError(1));
	fFunc_plus->SetParameter(2,fFunc->GetParameter(2) + fFunc->GetParError(2));
	fFunc_plus->SetParameter(3,fFunc->GetParameter(3) + fFunc->GetParError(3));
	Double_t FWHM_plus = fFunc_plus->GetX(fFunc_plus->GetParameter(0)*0.5,fFunc_plus->GetParameter(1), fMesonFitRange[1]) - fFunc_plus->GetX(fFunc_plus->GetParameter(0)*0.5,fMesonFitRange[0],fFunc_plus->GetParameter(1));

	//FWHM error -
	TF1* fFunc_minus;
	//	fFunc_minus = fFunc;
	fFunc_minus = new TF1("fFunc_minus","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);
	fFunc_minus->SetParameter(0,fFunc->GetParameter(0) - fFunc->GetParError(0));
	fFunc_minus->SetParameter(1,fFunc->GetParameter(1) - fFunc->GetParError(1));
	fFunc_minus->SetParameter(2,fFunc->GetParameter(2) - fFunc->GetParError(2));
	fFunc_minus->SetParameter(3,fFunc->GetParameter(3) - fFunc->GetParError(3));
	Double_t FWHM_minus =  fFunc_minus->GetX(fFunc_minus->GetParameter(0)*0.5,fFunc_minus->GetParameter(1), fMesonFitRange[1]) -fFunc_minus->GetX(fFunc_minus->GetParameter(0)*0.5,fMesonFitRange[0],fFunc_minus->GetParameter(1));

	Double_t Error1 = TMath::Abs(fFWHMFunc-FWHM_plus);
	Double_t Error2 = TMath::Abs(fFWHMFunc-FWHM_minus);

	if(Error1>=Error2) fFWHMFuncError = Error1;
	if(Error1<Error2) fFWHMFuncError = Error2;
}

//Crystal ball function for signal +linear background, parameters are 0:normalization,1:mean,2:sigma,3:n,4:alpha;
Double_t CrystalBallBck(Double_t *x,Double_t *par) {

	Double_t t = (x[0]-par[1])/par[2];
	if (par[4] < 0) t = -t;

	Double_t absAlpha = fabs((Double_t)par[4]);

	if (t >= -absAlpha) {
		return par[0]*exp(-0.5*t*t)+par[5]+par[6]*x[0];
	}
	else {
		Double_t a =  TMath::Power(par[3]/absAlpha,par[3])*exp(-0.5*absAlpha*absAlpha);
		Double_t b= par[3]/absAlpha - absAlpha;

		return par[0]*(a/TMath::Power(b - t, par[3]))+par[5]+par[6]*x[0];
	}
}


//Crystal ball function for signal, parameters are 0:normalization,1:mean,2:sigma,3:n,4:alpha;
Double_t CrystalBall(Double_t *x,Double_t *par) {

	Double_t t = (x[0]-par[1])/par[2];
	if (par[4] < 0) t = -t;

	Double_t absAlpha = fabs((Double_t)par[4]);

	if (t >= -absAlpha) {
		return par[0]*exp(-0.5*t*t);
	}
	else {
		Double_t a =  TMath::Power(par[3]/absAlpha,par[3])*exp(-0.5*absAlpha*absAlpha);
		Double_t b= par[3]/absAlpha - absAlpha;

		return par[0]*(a/TMath::Power(b - t, par[3]));
	}
}


void Delete(){
	delete fBinsPt;
	delete fBGFitRange;
	delete fBGFitRangeLeft;
	delete fMesonPlotRange;
	delete fMesonIntDeltaRange;
	delete fMesonIntDeltaRangeWide;
	delete fMesonIntDeltaRangeNarrow;
	delete fMesonMassRange;
	delete fMesonFitRange;
	delete fMesonWidthRange;
	delete fMesonLambdaTailRange;
	delete fNRebin;
	delete fGGYields;
	delete fBckYields;
	delete fMesonYields;
	delete fMesonTrueYields;
	delete fMesonTrueYieldsError;
	delete fMesonTrueGGYields;
	delete fMesonTrueGGYieldsError;
	delete fMesonYieldsFunc;
	delete fMesonYieldsResidualBckFunc ;
	delete fMesonYieldsCorResidualBckFunc;
	delete fMesonYieldsPerEvent;
	delete fMesonMass;
	delete fMesonWidth;
	delete fMesonSB;
	delete fMesonSign;
	delete fMesonTrueSB;
	delete fMesonTrueSign;
	delete fMesonFWHM;
	delete fGGYieldsLeft;
	delete fBckYieldsLeft;
	delete fMesonYieldsLeft;
	delete fMesonYieldsFuncLeft;
	delete fMesonYieldsResidualBckFuncLeft;
	delete fMesonYieldsCorResidualBckFuncLeft;
	delete fMesonYieldsLeftPerEvent;
	delete fMesonMassLeft;
	delete fMesonWidthLeft;
	delete fMesonSBLeft;
	delete fMesonSignLeft;
	delete fMesonFWHMLeft;
	delete fGGYieldsNarrow;
	delete fBckYieldsNarrow;
	delete fMesonYieldsNarrow;
	delete fMesonTrueYieldsNarrow;
	delete fMesonTrueYieldsNarrowError;
	delete fMesonTrueGGYieldsNarrow;
	delete fMesonTrueGGYieldsNarrowError;
	delete fMesonYieldsFuncNarrow;
	delete fMesonYieldsResidualBckFuncNarrow;
	delete fMesonYieldsCorResidualBckFuncNarrow;
	delete fMesonYieldsPerEventNarrow;
	delete fMesonSBNarrow;
	delete fMesonSignNarrow;
	delete fGGYieldsLeftNarrow;
	delete fBckYieldsLeftNarrow;
	delete fMesonYieldsLeftNarrow;
	delete fMesonYieldsFuncLeftNarrow;
	delete fMesonYieldsResidualBckFuncLeftNarrow;
	delete fMesonYieldsCorResidualBckFuncLeftNarrow;
	delete fMesonYieldsLeftPerEventNarrow;
	delete fMesonSBLeftNarrow;
	delete fMesonSignLeftNarrow;
	delete fGGYieldsWide;
	delete fBckYieldsWide;
	delete fMesonYieldsWide;
	delete fMesonTrueYieldsWide;
	delete fMesonTrueYieldsWideError;
	delete fMesonTrueGGYieldsWide;
	delete fMesonTrueGGYieldsWideError;
	delete fMesonYieldsFuncWide;
	delete fMesonYieldsResidualBckFuncWide;
	delete fMesonYieldsCorResidualBckFuncWide;
	delete fMesonYieldsPerEventWide;
	delete fMesonSBWide;
	delete fMesonSignWide;
	delete fGGYieldsLeftWide;
	delete fBckYieldsLeftWide;
	delete fMesonYieldsLeftWide;
	delete fMesonYieldsFuncLeftWide;
	delete fMesonYieldsResidualBckFuncLeftWide;
	delete fMesonYieldsCorResidualBckFuncLeftWide;
	delete fMesonYieldsLeftPerEventWide;
	delete fMesonSBLeftWide;
	delete fMesonSignLeftWide;
	delete fGGYieldsError;
	delete fBckYieldsError;
	delete fMesonYieldsError;
	delete fMesonYieldsFuncError;
	delete fMesonYieldsResidualBckFuncError;
	delete fMesonYieldsCorResidualBckFuncError;
	delete fMesonYieldsPerEventError;
	delete fMesonMassError;
	delete fMesonWidthError;
	delete fMesonSBError;
	delete fMesonSignError;
	delete fMesonTrueSBError;
	delete fMesonTrueSignError;
	delete fMesonFWHMError;
	delete fGGYieldsLeftError;
	delete fBckYieldsLeftError;
	delete fMesonYieldsLeftError;
	delete fMesonYieldsFuncLeftError;
	delete fMesonYieldsResidualBckFuncLeftError;
	delete fMesonYieldsCorResidualBckFuncLeftError;
	delete fMesonYieldsLeftPerEventError;
	delete fMesonMassLeftError;
	delete fMesonWidthLeftError;
	delete fMesonSBLeftError;
	delete fMesonSignLeftError;
	delete fMesonFWHMLeftError;
	delete fGGYieldsNarrowError;
	delete fBckYieldsNarrowError;
	delete fMesonYieldsNarrowError;
	delete fMesonYieldsFuncNarrowError;
	delete fMesonYieldsResidualBckFuncNarrowError;
	delete fMesonYieldsCorResidualBckFuncNarrowError;
	delete fMesonYieldsPerEventNarrowError;
	delete fGGYieldsLeftNarrowError;
	delete fBckYieldsLeftNarrowError;
	delete fMesonYieldsLeftNarrowError;
	delete fMesonYieldsFuncLeftNarrowError;
	delete fMesonYieldsResidualBckFuncLeftNarrowError;
	delete fMesonYieldsCorResidualBckFuncLeftNarrowError;
	delete fMesonYieldsLeftPerEventNarrowError;
	delete fMesonSBLeftNarrowError;
	delete fMesonSignLeftNarrowError;
	delete fGGYieldsWideError;
	delete fBckYieldsWideError;
	delete fMesonYieldsWideError;
	delete fMesonYieldsFuncWideError;
	delete fMesonYieldsResidualBckFuncWideError;
	delete fMesonYieldsCorResidualBckFuncWideError;
	delete fMesonYieldsPerEventWideError;
	delete fMesonSBWideError;
	delete fMesonSignWideError;
	delete fGGYieldsLeftWideError;
	delete fBckYieldsLeftWideError;
	delete fMesonYieldsLeftWideError;
	delete fMesonYieldsFuncLeftWideError;
	delete fMesonYieldsResidualBckFuncLeftWideError;
	delete fMesonYieldsCorResidualBckFuncLeftWideError;
	delete fMesonYieldsLeftPerEventWideError;
	delete fMesonSBLeftWideError;
	delete fMesonSignLeftWideError;
	delete fHistoMappingTrueMesonInvMassPtBins;
	delete fHistoMappingGGInvMassPtBin;
	delete fHistoMappingEEInvMassPtBin;
	delete fHistoMappingBackInvMassPtBin;
	delete fHistoMappingBackNormInvMassPtBin;
	delete fHistoMappingSignalInvMassPtBin;
	delete fHistoMappingRatioSBInvMassPtBin;
	delete fFitSignalInvMassPtBin;
	delete fFitSignalPeakPosInvMassPtBin;
	delete fFitBckInvMassPtBin;
	delete fHistoMappingBackNormInvMassLeftPtBin;
	delete fHistoMappingSignalInvMassLeftPtBin;
	delete fFitInvMassLeftPtBin;
	delete fFitSignalPeakPosInvMassLeftPtBin;
	delete fFitBckInvMassLeftPtBin;
	delete fFitRatioInvMassPtBin;
	delete fMesonFWHMAlpha01;
	delete fMesonFWHMAlpha01Error;
	delete fHistoWeightsBGZbinVsMbin;
	delete fHistoFillPerEventBGZbinVsMbin;
	//
}


